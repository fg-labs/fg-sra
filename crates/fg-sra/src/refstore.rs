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
use fg_sra_vdb::reference::ReferenceList;
use rustc_hash::FxHashMap;

/// The `INSDC:4na:map:CHARSET`: index by 4na code (`0..=15`); code 0 renders `.`.
const CHARSET_4NA: &[u8; 16] = b".ACMGRSVTWYHKDBN";

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
        let raw = obj.read_all().with_context(|| format!("read reference {idx}"))?;
        let expected = obj.seq_length().with_context(|| format!("seq_length {idx}"))? as usize;
        anyhow::ensure!(
            raw.len() == expected,
            "reference {idx}: read {} of {expected} bases",
            raw.len()
        );
        map_4na_to_ascii(&raw, &mut ascii);
        let wrap = obj.circular().with_context(|| format!("circular {idx}"))?
            || !obj.external().with_context(|| format!("external {idx}"))?;
        seqs.insert(idx, RefSeq { bases: std::mem::take(&mut ascii), wrap });
    }
    Ok(ReferenceStore::from_parts(seqs))
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
