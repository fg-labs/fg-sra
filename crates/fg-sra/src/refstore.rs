//! In-memory reference sequences, preloaded single-threaded, for reconstructing
//! aligned reads on worker threads without touching the VDB `REFERENCE` table.
//!
//! Reading the virtual `READ` column makes libncbi-vdb reconstruct through an
//! internal reference sub-cursor whose blob cache is not thread-safe. Instead we
//! load every needed reference into memory once (on the main thread, via the
//! serial [`ReferenceObj`] reader), converted to the `(ascii)READ` alphabet, and
//! hand workers a `&ReferenceStore` (plain `Vec<u8>`s — `Sync`) to slice.

use std::fmt;

use anyhow::{Context, Result};
use fg_sra_vdb::reference::{ReferenceList, ReferenceObj};
use rustc_hash::FxHashMap;

/// The `INSDC:4na:map:CHARSET`: index by 4na code (`0..=15`); code 0 renders `.`.
pub(crate) const CHARSET_4NA: &[u8; 16] = b".ACMGRSVTWYHKDBN";

/// Map `INSDC:4na:bin` codes to the CHARSET ASCII alphabet, in place into `dst`.
pub fn map_4na_to_ascii(src: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(src.len());
    dst.extend(src.iter().map(|&b| CHARSET_4NA[(b & 0x0F) as usize]));
}

/// Reference window length used by `READ` reconstruction (`get_ref_len_2`):
/// `has_ref_offset.len() + sum(ref_offset)`. May be negative for malformed input.
#[must_use]
pub fn ref_window_len(has_ref_offset_len: usize, ref_offset: &[i32]) -> i64 {
    let sum: i64 = ref_offset.iter().map(|&o| i64::from(o)).sum();
    has_ref_offset_len as i64 + sum
}

/// Failure extracting a reference window.
#[derive(Debug, PartialEq, Eq)]
pub enum WindowError {
    /// No preloaded reference for this index.
    UnknownReference(u32),
    /// Computed reference length was negative.
    NegativeLength(i64),
    /// A wrapped window longer than the reference itself (malformed offsets).
    WindowTooLong {
        /// Requested window length, in bases.
        requested: usize,
        /// Maximum accepted length (the reference's own length).
        max: usize,
    },
}

impl fmt::Display for WindowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownReference(idx) => write!(f, "no preloaded reference for index {idx}"),
            Self::NegativeLength(len) => write!(f, "negative reference window length {len}"),
            Self::WindowTooLong { requested, max } => write!(
                f,
                "wrapped reference window length {requested} exceeds reference length {max}"
            ),
        }
    }
}

impl std::error::Error for WindowError {}

/// A preloaded reference sequence in the CHARSET (ASCII) alphabet.
pub struct RefSeq {
    /// Bases mapped through `CHARSET_4NA`.
    pub bases: Vec<u8>,
    /// Whether the sub-select wraps at the end (circular or locally stored).
    pub wrap: bool,
}

/// Preloaded references keyed by `ReferenceObj` index. `window` is on the
/// per-record hot path, so the map uses `FxHashMap` (as `matecache.rs` does for
/// its analogous per-record lookup) rather than the SipHash-backed std map.
pub struct ReferenceStore {
    seqs: FxHashMap<u32, RefSeq>,
}

impl ReferenceStore {
    /// Test/constructor helper.
    #[must_use]
    pub fn from_parts(seqs: FxHashMap<u32, RefSeq>) -> Self {
        Self { seqs }
    }

    /// One store holding the references of every store in `stores`, which hold different
    /// references.
    #[must_use]
    pub fn merge(stores: impl IntoIterator<Item = ReferenceStore>) -> Self {
        let mut seqs = FxHashMap::default();
        for store in stores {
            seqs.extend(store.seqs);
        }
        Self { seqs }
    }

    /// Number of references held.
    #[must_use]
    pub fn num_references(&self) -> usize {
        self.seqs.len()
    }

    /// Bases held across all references.
    #[must_use]
    pub fn total_bases(&self) -> u64 {
        self.seqs.values().map(|seq| seq.bases.len() as u64).sum()
    }

    /// Extract the `ref_len` reference bases starting at 0-based `ref_pos`.
    ///
    /// Borrows directly from the stored sequence when the window fits before the
    /// end; wraps into `scratch` when the reference is circular/local and the
    /// window runs off the end (mirroring `ref_sub_select`). A window that runs
    /// off a non-wrapping reference is truncated (the restoration loop then
    /// reports the out-of-range base exactly as VDB does).
    ///
    /// When `ref_pos` is at or past the end of the reference the result is an
    /// empty window, not an error: `ref_sub_select` skips offset validation and
    /// returns zero bases here, leaving `restore_read` to report any base that
    /// actually needs a reference position (a fully soft-clipped read needs
    /// none, and reconstructs successfully from its mismatch bases alone).
    ///
    /// Returns [`WindowError::WindowTooLong`] when a wrapped window would exceed
    /// the reference's own length: a single read never legitimately needs more
    /// reference bases than the reference contains, so such a length can only come
    /// from malformed offsets, and materializing it would allocate an unbounded span.
    pub fn window<'a>(
        &'a self,
        ref_idx: u32,
        ref_pos: i32,
        ref_len: i64,
        scratch: &'a mut Vec<u8>,
    ) -> Result<&'a [u8], WindowError> {
        let seq = self.seqs.get(&ref_idx).ok_or(WindowError::UnknownReference(ref_idx))?;
        if ref_len < 0 {
            return Err(WindowError::NegativeLength(ref_len));
        }
        let ref_len = ref_len as usize;
        let n = seq.bases.len();
        let start = usize::try_from(ref_pos).unwrap_or(usize::MAX);
        let avail = n.saturating_sub(start);
        if !seq.wrap || ref_len <= avail {
            // `start` may sit at/past the end (avail == 0), yielding an empty
            // borrow; `start.min(n)` keeps the slice indices in range.
            let begin = start.min(n);
            return Ok(&seq.bases[begin..begin + ref_len.min(avail)]);
        }
        // Circular/local reference: fill scratch by wrapping around. `start` is
        // taken modulo the length so a position at/past the end wraps rather
        // than indexing out of range.
        //
        // A single read never legitimately needs more reference bases than the
        // reference itself contains. A window longer than the reference can only
        // come from malformed offsets (e.g. a corrupt positive REF_OFFSET); the
        // reserve/loop below would otherwise allocate and iterate over a
        // multi-gibibyte span. Reject it rather than materialize it.
        if ref_len > n {
            return Err(WindowError::WindowTooLong { requested: ref_len, max: n });
        }
        scratch.clear();
        scratch.reserve(ref_len);
        if n > 0 {
            let mut pos = start % n;
            for _ in 0..ref_len {
                scratch.push(seq.bases[pos]);
                pos += 1;
                if pos == n {
                    pos = 0;
                }
            }
        }
        Ok(&scratch[..])
    }
}

/// Which reference each row of the `REFERENCE` table belongs to, so an alignment's `REF_ID`
/// (a `REFERENCE` row) and `REF_START` (an offset within that row) can be turned into its
/// [`ReferenceList`] index and position without reading the `REFERENCE` table on a worker
/// thread.
pub struct ReferenceRows {
    /// `(first_row, last_row, reference index)`, sorted by first row.
    ranges: Vec<(i64, i64, u32)>,
    /// Bases per `REFERENCE` row (`MAX_SEQ_LEN`); 0 when the table has no such column, in
    /// which case `REF_START` is already the position on the reference.
    max_seq_len: u32,
}

impl ReferenceRows {
    /// The row ranges of every reference in `reflist` (single-threaded), whose rows hold
    /// `max_seq_len` bases each.
    pub fn new(reflist: &ReferenceList, max_seq_len: u32) -> Result<Self> {
        let count = reflist.count()?;
        let mut ranges = Vec::with_capacity(count as usize);
        for idx in 0..count {
            let (first, last) =
                reflist.get(idx)?.id_range().with_context(|| format!("id range {idx}"))?;
            ranges.push((first, last, idx));
        }
        Ok(Self::from_ranges(ranges, max_seq_len))
    }

    /// Rows from `(first_row, last_row, reference index)` ranges, in any order.
    #[must_use]
    pub fn from_ranges(mut ranges: Vec<(i64, i64, u32)>, max_seq_len: u32) -> Self {
        ranges.sort_unstable();
        Self { ranges, max_seq_len }
    }

    /// The reference index and 0-based position of offset `row_offset` within `REFERENCE`
    /// row `row`, as the schema's `REF_POS` computes it; `None` if no reference covers the row
    /// or the position doesn't fit.
    #[must_use]
    pub fn locate(&self, row: i64, row_offset: i32) -> Option<(u32, i32)> {
        let (first, idx) = self.range_of(row)?;
        let position = (row - first) * i64::from(self.max_seq_len) + i64::from(row_offset);
        Some((idx, i32::try_from(position).ok()?))
    }

    /// The first row and index of the reference covering `row`.
    fn range_of(&self, row: i64) -> Option<(i64, u32)> {
        let after = self.ranges.partition_point(|&(first, _, _)| first <= row);
        let &(first, last, idx) = self.ranges.get(after.checked_sub(1)?)?;
        (row <= last).then_some((first, idx))
    }
}

/// Split references, given as `(index, length)`, into at most `groups` groups of similar total
/// length, for loading in parallel: longest first, each to the group with the least so far.
#[must_use]
pub fn split_by_length(references: &[(u32, u64)], groups: usize) -> Vec<Vec<u32>> {
    let mut by_length = references.to_vec();
    by_length.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut split: Vec<(u64, Vec<u32>)> =
        vec![(0, Vec::new()); groups.clamp(1, references.len().max(1))];
    for (idx, len) in by_length {
        let lightest =
            split.iter_mut().min_by_key(|(total, _)| *total).expect("at least one group");
        lightest.0 += len;
        lightest.1.push(idx);
    }
    split.into_iter().map(|(_, indices)| indices).filter(|g| !g.is_empty()).collect()
}

/// Preload the distinct references in `ref_indices` into memory (single-threaded).
///
/// Duplicate indices within a single call are loaded once; a reference that
/// recurs across separate calls (e.g. overlapping `--aligned-region` batches)
/// is loaded again per call.
pub fn preload_references(reflist: &ReferenceList, ref_indices: &[u32]) -> Result<ReferenceStore> {
    let mut seqs: FxHashMap<u32, RefSeq> = FxHashMap::default();
    let mut ascii = Vec::new();
    for &idx in ref_indices {
        if seqs.contains_key(&idx) {
            continue;
        }
        let obj = reflist.get(idx).with_context(|| format!("get reference {idx}"))?;
        let label = reference_label(&obj, idx);
        let raw = obj.read_all().with_context(|| format!("failed to read reference {label}"))?;
        let expected = obj.seq_length().with_context(|| format!("seq_length {label}"))? as usize;
        anyhow::ensure!(
            raw.len() == expected,
            "reference {label}: read {} of {expected} bases",
            raw.len()
        );
        map_4na_to_ascii(&raw, &mut ascii);
        let wrap = obj.circular().with_context(|| format!("circular {label}"))?
            || !obj.external().with_context(|| format!("external {label}"))?;
        seqs.insert(idx, RefSeq { bases: std::mem::take(&mut ascii), wrap });
    }
    Ok(ReferenceStore::from_parts(seqs))
}

/// A reference's accession and, when different, its name in the archive (e.g.
/// `NC_000001.11 (chr1)`), so a user can tell which reference to fetch; `#idx` if unreadable.
fn reference_label(obj: &ReferenceObj, idx: u32) -> String {
    match (obj.seq_id(), obj.name()) {
        (Ok(seq_id), Ok(name)) if name != seq_id && !name.is_empty() => {
            format!("{seq_id} ({name})")
        }
        (Ok(seq_id), _) => seq_id,
        (Err(_), _) => format!("#{idx}"),
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn reference_rows_find_the_reference_covering_a_row() {
        let rows = ReferenceRows::from_ranges(vec![(11, 20, 1), (1, 10, 0), (21, 21, 2)], 5000);
        let reference_of = |row| rows.locate(row, 0).map(|(idx, _)| idx);
        assert_eq!(reference_of(1), Some(0));
        assert_eq!(reference_of(10), Some(0));
        assert_eq!(reference_of(11), Some(1));
        assert_eq!(reference_of(21), Some(2));
    }

    #[test]
    fn reference_rows_outside_every_range_have_no_reference() {
        let rows = ReferenceRows::from_ranges(vec![(5, 10, 0), (20, 30, 1)], 5000);
        assert_eq!(rows.locate(4, 0), None);
        assert_eq!(rows.locate(15, 0), None);
        assert_eq!(rows.locate(31, 0), None);
    }

    #[test]
    fn locate_offsets_by_whole_rows_from_the_reference_start() {
        let rows = ReferenceRows::from_ranges(vec![(1, 10, 0), (11, 20, 1)], 5000);
        assert_eq!(rows.locate(1, 17), Some((0, 17)));
        assert_eq!(rows.locate(3, 17), Some((0, 10_017)));
        assert_eq!(rows.locate(11, 0), Some((1, 0)));
        assert_eq!(rows.locate(14, 4999), Some((1, 19_999)));
    }

    #[test]
    fn locate_without_a_row_length_uses_the_offset_as_the_position() {
        let rows = ReferenceRows::from_ranges(vec![(1, 10, 0)], 0);
        assert_eq!(rows.locate(4, 123), Some((0, 123)));
    }

    #[test]
    fn split_by_length_balances_total_length() {
        let groups = split_by_length(&[(0, 100), (1, 60), (2, 50), (3, 10)], 2);
        assert_eq!(groups, vec![vec![0, 3], vec![1, 2]]);
    }

    #[test]
    fn split_by_length_makes_no_more_groups_than_references() {
        let groups = split_by_length(&[(0, 5), (1, 5)], 8);
        assert_eq!(groups.len(), 2);
    }

    #[test]
    fn split_by_length_of_nothing_is_empty() {
        assert!(split_by_length(&[], 4).is_empty());
    }

    #[test]
    fn merged_stores_hold_every_reference() {
        let store = |idx: u32, bases: &[u8]| {
            let mut seqs = FxHashMap::default();
            seqs.insert(idx, RefSeq { bases: bases.to_vec(), wrap: false });
            ReferenceStore::from_parts(seqs)
        };
        let merged = ReferenceStore::merge([store(0, b"ACGT"), store(3, b"GG")]);
        assert_eq!(merged.num_references(), 2);
        assert_eq!(merged.total_bases(), 6);
    }

    #[test]
    fn charset_maps_common_codes() {
        let mut dst = Vec::new();
        map_4na_to_ascii(&[0, 1, 2, 4, 8, 15], &mut dst);
        assert_eq!(dst, b".ACGTN");
    }

    #[test]
    fn charset_maps_every_code_including_iupac() {
        // Pin the whole 4na alphabet (INSDC:4na:map:CHARSET), so a transposition
        // of any IUPAC ambiguity code is caught, not just A/C/G/T/N.
        let all: Vec<u8> = (0u8..16).collect();
        let mut dst = Vec::new();
        map_4na_to_ascii(&all, &mut dst);
        assert_eq!(dst, b".ACMGRSVTWYHKDBN");
    }

    #[test]
    fn charset_masks_high_bits() {
        let mut dst = Vec::new();
        map_4na_to_ascii(&[0x11, 0x12], &mut dst); // & 0x0F => 1, 2
        assert_eq!(dst, b"AC");
    }

    #[rstest]
    #[case(4, &[], 4)]
    #[case(4, &[2], 6)]
    #[case(5, &[-1], 4)]
    #[case(4, &[-2], 2)]
    #[case(245, &[4, -140], 109)]
    fn window_len(#[case] hro_len: usize, #[case] ref_offset: &[i32], #[case] expected: i64) {
        assert_eq!(ref_window_len(hro_len, ref_offset), expected);
    }

    fn store_of(bases: &[u8], wrap: bool) -> ReferenceStore {
        let mut m = FxHashMap::default();
        m.insert(0u32, RefSeq { bases: bases.to_vec(), wrap });
        ReferenceStore::from_parts(m)
    }

    #[test]
    fn window_in_range_borrows() {
        let store = store_of(b"ACGTACGT", false);
        let mut scratch = Vec::new();
        assert_eq!(store.window(0, 2, 3, &mut scratch).unwrap(), b"GTA");
        assert!(scratch.is_empty(), "no-wrap path must not touch scratch");
    }

    #[test]
    fn window_truncates_at_end_when_not_wrapping() {
        let store = store_of(b"ACGTACGT", false);
        let mut scratch = Vec::new();
        // pos 6 len 5 -> only "GT" available.
        assert_eq!(store.window(0, 6, 5, &mut scratch).unwrap(), b"GT");
    }

    #[test]
    fn window_wraps_when_circular_or_local() {
        let store = store_of(b"ACGTACGT", true);
        let mut scratch = Vec::new();
        // pos 6 len 5 -> "GT" + wrap "ACG".
        assert_eq!(store.window(0, 6, 5, &mut scratch).unwrap(), b"GTACG");
    }

    #[test]
    fn window_wraps_from_position_at_or_past_end() {
        // A circular position exactly at the end wraps to the start rather than
        // panicking; pos 8 (== len) reads "ACG" from the beginning.
        let store = store_of(b"ACGTACGT", true);
        let mut scratch = Vec::new();
        assert_eq!(store.window(0, 8, 3, &mut scratch).unwrap(), b"ACG");
    }

    #[test]
    fn window_at_or_past_end_is_empty_not_error() {
        // `ref_sub_select` returns an empty window (not an error) once ref_pos is
        // at/past the end of a non-wrapping reference; restore_read then reports
        // any base that truly needs a reference position.
        let store = store_of(b"AC", false);
        let mut scratch = Vec::new();
        assert_eq!(store.window(0, 2, 1, &mut scratch).unwrap(), b"");
        assert_eq!(store.window(0, 5, 3, &mut scratch).unwrap(), b"");
    }

    #[test]
    fn window_errors() {
        let store = store_of(b"AC", false);
        let mut scratch = Vec::new();
        assert_eq!(store.window(9, 0, 1, &mut scratch), Err(WindowError::UnknownReference(9)));
        assert_eq!(store.window(0, 0, -1, &mut scratch), Err(WindowError::NegativeLength(-1)));
    }

    #[test]
    fn window_rejects_wrap_longer_than_reference() {
        // A malformed offset can make the wrapped window exceed the reference
        // length; reject it rather than reserve/loop over an unbounded span.
        let store = store_of(b"ACGT", true);
        let mut scratch = Vec::new();
        // len == n is still allowed (a full turn around the reference).
        assert_eq!(store.window(0, 1, 4, &mut scratch).unwrap(), b"CGTA");
        // len > n is rejected.
        assert_eq!(
            store.window(0, 1, 5, &mut scratch),
            Err(WindowError::WindowTooLong { requested: 5, max: 4 })
        );
    }
}
