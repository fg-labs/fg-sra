//! Reconstruct an aligned read from the reference window and the stored
//! alignment deltas — a Rust port of libncbi-vdb's `align_restore_read_impl`.
//!
//! This lets fg-sra rebuild `PRIMARY_ALIGNMENT.READ` itself, from physically
//! stored columns, instead of reading the virtual `READ` column (which makes
//! libncbi-vdb reconstruct through an internal, unsynchronized MRU blob cache
//! that is not safe for concurrent reads).
//!
//! Ported from `libs/axf/align-restore-read.c` (`align_restore_read_impl`) in
//! NCBI ncbi-vdb (vendored 3.3.0). That source is in the **public domain** — a
//! "United States Government Work" by the National Center for Biotechnology
//! Information / National Library of Medicine, with a request to cite the
//! author. The schema function it implements is:
//!
//! ```text
//! ALIGN:align_restore_read( INSDC:4na:bin ref_read, bool has_mismatch,
//!     INSDC:4na:bin mismatch, bool has_ref_offset, I32 ref_offset [, U32 read_len] )
//! ```
//!
//! Because the CHARSET map (`.ACMGRSVTWYHKDBN`) is a bijection on 0..16 and the
//! restoration loop only *copies* elements, running it over ASCII bases (the
//! reference and `MISMATCH` mapped through CHARSET) yields exactly the bytes the
//! `(ascii)READ` column produces — no per-read conversion pass is needed.

use std::fmt;

use anyhow::{Result, bail, ensure};

use crate::refstore::CHARSET_4NA;

/// `READ_TYPE` bit (`INSDC:SRA:xread_type`) of a read aligned in sequencing orientation.
const READ_TYPE_FORWARD: u8 = 2;

/// `READ_TYPE` bit of a read aligned reverse-complemented.
const READ_TYPE_REVERSE: u8 = 4;

/// The complement of each `INSDC:dna:text` base: its 4na code bit-reversed (as the `map`
/// in `seq-restore-read.c`) and mapped back through CHARSET, so `A`↔`T`, `C`↔`G`, `M`↔`K`,
/// `R`↔`Y`, `V`↔`B`, `H`↔`D`, and `S`, `W`, `N` and `.` are their own. Other bytes are
/// left as they are.
static COMPLEMENT: [u8; 256] = complement_table();

/// Build [`COMPLEMENT`] at compile time.
const fn complement_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] = byte as u8;
        byte += 1;
    }
    let mut code = 0;
    while code < 16 {
        let reversed =
            ((code & 1) << 3) | ((code & 2) << 1) | ((code & 4) >> 1) | ((code & 8) >> 3);
        table[CHARSET_4NA[code] as usize] = CHARSET_4NA[reversed];
        code += 1;
    }
    table
}

/// Error reconstructing a read; each variant mirrors an `rcInconsistent` return
/// in the C implementation.
#[derive(Debug, PartialEq, Eq)]
pub enum RestoreReadError {
    /// `has_mismatch` and `has_ref_offset` have different lengths (checked
    /// before the loop in C).
    LengthMismatch { has_mismatch: usize, has_ref_offset: usize },
    /// A `has_ref_offset` bit is set but `ref_offset` has no more entries.
    RefOffsetExhausted { pos: usize },
    /// A `has_mismatch` bit is set but `mismatch` has no more entries.
    MismatchExhausted { pos: usize },
    /// A matching base needs a reference index outside the window (`< 0` or
    /// `>= ref_read.len()`).
    RefIndexOutOfRange { pos: usize, ref_index: i64, ref_len: usize },
}

impl fmt::Display for RestoreReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch { has_mismatch, has_ref_offset } => write!(
                f,
                "has_mismatch length {has_mismatch} != has_ref_offset length {has_ref_offset}"
            ),
            Self::RefOffsetExhausted { pos } => {
                write!(f, "ref_offset exhausted at read position {pos}")
            }
            Self::MismatchExhausted { pos } => {
                write!(f, "mismatch bases exhausted at read position {pos}")
            }
            Self::RefIndexOutOfRange { pos, ref_index, ref_len } => write!(
                f,
                "read position {pos} needs reference index {ref_index} outside window of length {ref_len}"
            ),
        }
    }
}

impl std::error::Error for RestoreReadError {}

/// Reconstruct a read into `out`, replacing its contents; returns `Ok(())` on
/// success. `out` ends up holding `has_mismatch.len()` bytes.
///
/// - `ref_read`: reference bases spanning the alignment (same alphabet as the
///   desired output; ASCII here).
/// - `has_mismatch` / `has_ref_offset`: one 0/1 byte per read base; must be the
///   same length (the read length).
/// - `mismatch`: one base per set bit in `has_mismatch`, in order.
/// - `ref_offset`: one value per set bit in `has_ref_offset` (subject to the
///   `bi >= 0` skip rule), in order.
/// - `read_len`: empty for ploidy 1 (PRIMARY/SECONDARY alignments); `N` entries
///   drive `N`-ploidy reconstruction (EVIDENCE alignments, which fg-sra does not
///   read — kept for fidelity with the C source).
///
/// `out` is cleared and filled with `has_mismatch.len()` bytes.
#[inline]
pub fn restore_read(
    ref_read: &[u8],
    has_mismatch: &[u8],
    mismatch: &[u8],
    has_ref_offset: &[u8],
    ref_offset: &[i32],
    read_len: &[u32],
    out: &mut Vec<u8>,
) -> Result<(), RestoreReadError> {
    if has_mismatch.len() != has_ref_offset.len() {
        return Err(RestoreReadError::LengthMismatch {
            has_mismatch: has_mismatch.len(),
            has_ref_offset: has_ref_offset.len(),
        });
    }
    let dst_len = has_mismatch.len();
    // Most reads have no reference offsets and lie within their window; they are the reference
    // bases with the mismatches patched in, which the general loop below finds one base at a
    // time. A read that overhangs its window is left to the loop, which reports the first
    // position that needs a reference base outside it.
    if ref_offset.is_empty()
        && read_len.is_empty()
        && dst_len <= ref_read.len()
        && all_zero(has_ref_offset)
    {
        return restore_ungapped(ref_read, has_mismatch, mismatch, out);
    }
    out.clear();
    out.reserve(dst_len);

    // Empty `read_len` => ploidy 1 (the aligned path); otherwise one read per entry.
    let mut ploidy = if read_len.is_empty() { 1 } else { read_len.len() };
    let mut read_len_idx = 0usize;

    let ref_len = ref_read.len() as i64;
    let mut mmi = 0usize; // next index into `mismatch`
    let mut roi = 0usize; // next index into `ref_offset`
    let mut rri: i64 = 0; // reference-window index (may go negative)
    let mut bi: i64 = 0; // last applied ref offset; skip `has_ref_offset` while < 0
    let mut rl: u32 = 1; // 1-based position within the current read (for ploidy)

    for (di, (&hro, &hmm)) in has_ref_offset.iter().zip(has_mismatch).enumerate() {
        if hro != 0 && bi >= 0 {
            let &off =
                ref_offset.get(roi).ok_or(RestoreReadError::RefOffsetExhausted { pos: di })?;
            bi = i64::from(off);
            rri += bi;
            roi += 1;
        }

        if hmm != 0 {
            let &base = mismatch.get(mmi).ok_or(RestoreReadError::MismatchExhausted { pos: di })?;
            out.push(base);
            mmi += 1;
        } else if rri < 0 || rri >= ref_len {
            return Err(RestoreReadError::RefIndexOutOfRange {
                pos: di,
                ref_index: rri,
                ref_len: ref_read.len(),
            });
        } else {
            out.push(ref_read[usize::try_from(rri).expect("checked 0 <= rri < ref_len")]);
        }

        if ploidy > 1 && rl == read_len[read_len_idx] {
            rri = -1;
            rl = 0;
            ploidy -= 1;
            read_len_idx += 1;
        }

        rri += 1;
        rl += 1;
        bi += 1;
    }

    Ok(())
}

/// Rebuild a read without reference offsets, whose window `ref_read` holds at least
/// `has_mismatch.len()` bases: the first bases of the window, with each mismatch written over
/// the position of its set `has_mismatch` flag. Gives what [`restore_read`]'s general loop
/// gives for the same input, including its error for too few mismatches.
fn restore_ungapped(
    ref_read: &[u8],
    has_mismatch: &[u8],
    mismatch: &[u8],
    out: &mut Vec<u8>,
) -> Result<(), RestoreReadError> {
    out.clear();
    out.extend_from_slice(&ref_read[..has_mismatch.len()]);
    let mut next_mismatch = 0;
    let mut pos = 0;
    for flags in has_mismatch.chunks(8) {
        if all_zero(flags) {
            pos += flags.len();
            continue;
        }
        for &flag in flags {
            if flag != 0 {
                let &base = mismatch
                    .get(next_mismatch)
                    .ok_or(RestoreReadError::MismatchExhausted { pos })?;
                out[pos] = base;
                next_mismatch += 1;
            }
            pos += 1;
        }
    }
    Ok(())
}

/// Whether every byte of `flags` is zero, tested eight bytes at a time.
fn all_zero(flags: &[u8]) -> bool {
    let mut words = flags.chunks_exact(8);
    let whole = words.by_ref().all(|word| word.iter().fold(0, |acc, &b| acc | b) == 0);
    whole && words.remainder().iter().all(|&b| b == 0)
}

/// Rebuild a spot's bases, as a cSRA `SEQUENCE.READ`, into `out` — a port of libncbi-vdb's
/// `seq_restore_read_impl` (`libs/axf/seq-restore-read.c`, public domain like
/// `align-restore-read.c` above).
///
/// - `cmp_read`: `CMP_READ`, the bases of the spot's unaligned reads, end to end.
/// - `align_ids`: each read's `PRIMARY_ALIGNMENT_ID`; 0 for an unaligned read.
/// - `read_lens`, `read_types`: each read's `READ_LEN` and `READ_TYPE`.
/// - `aligned_read(id, buf)`: fills `buf` with alignment `id`'s bases in reference
///   orientation, as `PRIMARY_ALIGNMENT.READ` would give them.
///
/// Reads go in slot order. An unaligned read takes the next bases of `cmp_read`. An aligned
/// read must be exactly `READ_LEN` long, and is copied if its type has the FORWARD bit or
/// reverse-complemented back to sequencing orientation if it has REVERSE; neither is an
/// error. As in the C code, when `cmp_read` holds every base of the spot it is used as is.
pub fn restore_spot(
    cmp_read: &[u8],
    align_ids: &[i64],
    read_lens: &[u32],
    read_types: &[u8],
    mut aligned_read: impl FnMut(i64, &mut Vec<u8>) -> Result<()>,
    aligned: &mut Vec<u8>,
    out: &mut Vec<u8>,
) -> Result<()> {
    out.clear();
    let total: usize = read_lens.iter().map(|&len| len as usize).sum();
    if total == cmp_read.len() {
        out.extend_from_slice(cmp_read);
        return Ok(());
    }
    ensure!(
        align_ids.len() == read_lens.len() && read_types.len() == read_lens.len(),
        "PRIMARY_ALIGNMENT_ID has {} values, READ_TYPE {} and READ_LEN {}",
        align_ids.len(),
        read_types.len(),
        read_lens.len()
    );
    out.reserve(total);
    let mut unaligned = cmp_read;
    for ((&align_id, &len), &read_type) in align_ids.iter().zip(read_lens).zip(read_types) {
        let len = len as usize;
        if align_id <= 0 {
            ensure!(
                unaligned.len() >= len,
                "CMP_READ has {} bases left but the next unaligned read has {len}",
                unaligned.len()
            );
            let (read, rest) = unaligned.split_at(len);
            out.extend_from_slice(read);
            unaligned = rest;
            continue;
        }
        aligned_read(align_id, aligned)?;
        ensure!(
            aligned.len() == len,
            "alignment {align_id} has {} bases but its read has {len}",
            aligned.len()
        );
        if read_type & READ_TYPE_FORWARD != 0 {
            out.extend_from_slice(aligned);
        } else if read_type & READ_TYPE_REVERSE != 0 {
            out.extend(aligned.iter().rev().map(|&base| COMPLEMENT[usize::from(base)]));
        } else {
            bail!(
                "alignment {align_id}'s read has neither orientation bit (READ_TYPE {read_type})"
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::all_match_4m(b"ACGT", &[0, 0, 0, 0], b"", &[0, 0, 0, 0], &[], b"ACGT")]
    #[case::single_mismatch(b"ACGT", &[0, 1, 0, 0], b"T", &[0, 0, 0, 0], &[], b"ATGT")]
    #[case::left_soft_clip_2s2m(b"AC", &[1, 1, 0, 0], b"GG", &[1, 0, 0, 0], &[-2], b"GGAC")]
    #[case::insertion_2m1i2m(b"ACGT", &[0, 0, 1, 0, 0], b"T", &[0, 0, 1, 0, 0], &[-1], b"ACTGT")]
    #[case::deletion_2m2d2m(b"ACGTAC", &[0, 0, 0, 0], b"", &[0, 0, 1, 0], &[2], b"ACAC")]
    #[case::right_soft_clip_2m2s(b"AC", &[0, 0, 1, 1], b"NN", &[0, 0, 1, 0], &[-2], b"ACNN")]
    #[case::bi_negative_skips_hro(b"AC", &[1, 1, 1, 0, 0], b"GGG", &[1, 1, 0, 0, 0], &[-3, 5], b"GGGAC")]
    #[case::dot_for_4na_zero(b".A", &[0, 0], b"", &[0, 0], &[], b".A")]
    #[case::empty_read(b"", &[], b"", &[], &[], b"")]
    fn reconstructs(
        #[case] ref_read: &[u8],
        #[case] has_mismatch: &[u8],
        #[case] mismatch: &[u8],
        #[case] has_ref_offset: &[u8],
        #[case] ref_offset: &[i32],
        #[case] expected: &[u8],
    ) {
        let mut out = Vec::new();
        restore_read(ref_read, has_mismatch, mismatch, has_ref_offset, ref_offset, &[], &mut out)
            .expect("reconstruction should succeed");
        assert_eq!(out, expected);
    }

    #[test]
    fn ploidy_two_resets_reference_index() {
        let mut out = Vec::new();
        restore_read(b"AC", &[0, 0, 0, 0], b"", &[0, 0, 0, 0], &[], &[2, 2], &mut out)
            .expect("reconstruction should succeed");
        assert_eq!(out, b"ACAC");
    }

    #[rstest]
    #[case::length_mismatch(
        &[0, 0, 0], &[0, 0, 0, 0], &[],
        RestoreReadError::LengthMismatch { has_mismatch: 3, has_ref_offset: 4 }
    )]
    #[case::ref_offset_exhausted(
        &[0, 0, 0, 0], &[1, 0, 0, 0], &[],
        RestoreReadError::RefOffsetExhausted { pos: 0 }
    )]
    #[case::ref_index_negative(
        &[0, 0, 0, 0], &[1, 0, 0, 0], &[-1],
        RestoreReadError::RefIndexOutOfRange { pos: 0, ref_index: -1, ref_len: 2 }
    )]
    #[case::ref_index_past_end(
        &[0, 0, 0, 0], &[0, 0, 0, 0], &[],
        RestoreReadError::RefIndexOutOfRange { pos: 2, ref_index: 2, ref_len: 2 }
    )]
    fn errors_without_mismatch(
        #[case] has_mismatch: &[u8],
        #[case] has_ref_offset: &[u8],
        #[case] ref_offset: &[i32],
        #[case] expected: RestoreReadError,
    ) {
        let mut out = Vec::new();
        let err = restore_read(b"AC", has_mismatch, b"", has_ref_offset, ref_offset, &[], &mut out)
            .expect_err("reconstruction should fail");
        assert_eq!(err, expected);
    }

    const B_FWD: u8 = 1 | READ_TYPE_FORWARD;
    const B_REV: u8 = 1 | READ_TYPE_REVERSE;

    /// Restore a spot whose alignment `id` has bases `"ACGGT"` followed by `id - 1` `N`s.
    fn restore(cmp_read: &[u8], align_ids: &[i64], lens: &[u32], types: &[u8]) -> Result<Vec<u8>> {
        let aligned_read = |id: i64, buf: &mut Vec<u8>| {
            buf.clear();
            buf.extend_from_slice(b"ACGGT");
            buf.extend(std::iter::repeat_n(b'N', (id - 1) as usize));
            Ok(())
        };
        let (mut aligned, mut out) = (Vec::new(), Vec::new());
        restore_spot(cmp_read, align_ids, lens, types, aligned_read, &mut aligned, &mut out)?;
        Ok(out)
    }

    #[test]
    fn complement_follows_the_4na_bit_reversal() {
        let complemented: Vec<u8> =
            CHARSET_4NA.iter().map(|&b| COMPLEMENT[usize::from(b)]).collect();
        assert_eq!(complemented, b".TGKCYSBAWRDMHVN");
    }

    #[test]
    fn spot_of_unaligned_reads_is_cmp_read() {
        let out = restore(b"AAAACC", &[0, 0], &[4, 2], &[B_FWD, B_REV]).unwrap();
        assert_eq!(out, b"AAAACC");
    }

    #[test]
    fn forward_aligned_read_is_copied() {
        let out = restore(b"", &[1], &[5], &[B_FWD]).unwrap();
        assert_eq!(out, b"ACGGT");
    }

    #[test]
    fn reverse_aligned_read_is_reverse_complemented() {
        let out = restore(b"", &[1], &[5], &[B_REV]).unwrap();
        assert_eq!(out, b"ACCGT");
    }

    #[test]
    fn aligned_and_unaligned_reads_are_spliced_in_slot_order() {
        let out = restore(b"TTT", &[0, 1], &[3, 5], &[B_FWD, B_REV]).unwrap();
        assert_eq!(out, b"TTTACCGT");
        let out = restore(b"TTT", &[1, 0], &[5, 3], &[B_FWD, B_REV]).unwrap();
        assert_eq!(out, b"ACGGTTTT");
    }

    #[test]
    fn aligned_read_of_the_wrong_length_is_an_error() {
        assert!(restore(b"", &[2], &[5], &[B_FWD]).is_err());
    }

    #[test]
    fn aligned_read_without_an_orientation_is_an_error() {
        assert!(restore(b"", &[1], &[5], &[1]).is_err());
    }

    #[test]
    fn cmp_read_too_short_for_the_unaligned_reads_is_an_error() {
        assert!(restore(b"TT", &[0, 1], &[3, 5], &[B_FWD, B_FWD]).is_err());
    }

    #[test]
    fn mismatch_exhausted() {
        let mut out = Vec::new();
        let err = restore_read(b"AC", &[1, 0, 0, 0], b"", &[0, 0, 0, 0], &[], &[], &mut out)
            .expect_err("reconstruction should fail");
        assert_eq!(err, RestoreReadError::MismatchExhausted { pos: 0 });
    }

    /// A 40-base reference window of `A`s, and `has_mismatch` flags set at `positions`.
    fn ungapped_input(positions: &[usize]) -> (Vec<u8>, Vec<u8>) {
        let mut has_mismatch = vec![0u8; 40];
        for &pos in positions {
            has_mismatch[pos] = 1;
        }
        (vec![b'A'; 40], has_mismatch)
    }

    #[test]
    fn ungapped_read_without_mismatches_is_the_reference_window() {
        let (reference, has_mismatch) = ungapped_input(&[]);
        let mut out = Vec::new();
        restore_read(&reference[..37], &has_mismatch[..37], b"", &[0; 37], &[], &[], &mut out)
            .unwrap();
        assert_eq!(out, vec![b'A'; 37]);
    }

    #[test]
    fn ungapped_mismatches_are_patched_at_their_flags() {
        let (reference, has_mismatch) = ungapped_input(&[0, 7, 8, 15, 16, 39]);
        let mut out = Vec::new();
        restore_read(&reference, &has_mismatch, b"CGTCGT", &[0; 40], &[], &[], &mut out).unwrap();
        let mut expected = vec![b'A'; 40];
        for (pos, base) in [0, 7, 8, 15, 16, 39].into_iter().zip(b"CGTCGT") {
            expected[pos] = *base;
        }
        assert_eq!(out, expected);
    }

    #[test]
    fn ungapped_read_shorter_than_its_window_takes_the_first_bases() {
        let mut out = Vec::new();
        restore_read(b"ACGTACGT", &[0, 1, 0], b"T", &[0, 0, 0], &[], &[], &mut out).unwrap();
        assert_eq!(out, b"ATG");
    }

    #[test]
    fn ungapped_read_reports_mismatch_shortage_at_its_flag() {
        let (reference, has_mismatch) = ungapped_input(&[3, 20]);
        let mut out = Vec::new();
        let err = restore_read(&reference, &has_mismatch, b"C", &[0; 40], &[], &[], &mut out)
            .expect_err("one mismatch for two flags");
        assert_eq!(err, RestoreReadError::MismatchExhausted { pos: 20 });
    }

    #[test]
    fn ungapped_read_reuses_the_output_buffer() {
        let mut out = b"stale".to_vec();
        restore_read(b"ACGT", &[0, 0], b"", &[0, 0], &[], &[], &mut out).unwrap();
        assert_eq!(out, b"AC");
    }

    #[test]
    fn all_zero_tests_every_byte() {
        assert!(all_zero(&[]));
        assert!(all_zero(&[0; 17]));
        for pos in [0, 7, 8, 15, 16] {
            let mut flags = [0u8; 17];
            flags[pos] = 1;
            assert!(!all_zero(&flags), "flag at {pos}");
        }
    }
}
