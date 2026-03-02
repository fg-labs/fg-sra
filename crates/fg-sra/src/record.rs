//! SAM/BAM record building and formatting.
//!
//! Builds complete SAM lines in a per-record `Vec<u8>` buffer, writing
//! all fields in a single pass rather than multiple formatted-print calls.

use std::io::Write;

use crate::matecache::MateInfo;

/// Column data read from VDB for an aligned record.
///
/// All fields correspond to VDB alignment table columns.
pub struct AlignedColumns {
    /// Read name from VDB (`SEQ_NAME`).
    pub seq_name: String,
    /// SAM flags (`SAM_FLAGS`).
    pub sam_flags: u32,
    /// CIGAR string (short or long form).
    pub cigar: String,
    /// Mate alignment ID (`MATE_ALIGN_ID`), 0 if unpaired.
    pub mate_align_id: i64,
    /// Mate reference name (`MATE_REF_NAME`).
    pub mate_ref_name: String,
    /// Mate reference position, 0-based (`MATE_REF_POS`).
    pub mate_ref_pos: i32,
    /// Template length (`TEMPLATE_LEN`).
    pub template_len: i32,
    /// Read sequence (`SAM_QUALITY` column provides pre-encoded quality).
    pub read: String,
    /// Quality string (already Phred+33 encoded from `SAM_QUALITY`).
    pub quality: String,
    /// Edit distance (`EDIT_DISTANCE`) for NM tag.
    pub edit_distance: u32,
    /// Read group (`SEQ_SPOT_GROUP`).
    pub spot_group: String,
    /// Alignment count (`ALIGNMENT_COUNT`) for NH tag, 0 if unavailable.
    pub alignment_count: u8,
    /// Read filter value (`READ_FILTER`), `None` if column unavailable.
    pub read_filter: Option<u8>,
}

impl AlignedColumns {
    /// Strip paired-end flags when the mate has no alignment (`mate_align_id == 0`).
    ///
    /// Matches sam-dump behavior: reads without a mate alignment are output as
    /// unpaired by clearing PAIRED, PROPER_PAIR, MATE_UNMAPPED, MATE_REVERSE,
    /// FIRST_IN_PAIR, and LAST_IN_PAIR flag bits.
    pub fn strip_paired_flags(&mut self) {
        self.sam_flags &= !sam_flags::PAIRED_MASK;
    }
}

/// Column data for an unaligned record from the SEQUENCE table.
///
/// Uses borrowed slices to avoid per-record allocations in the hot loop.
/// The caller owns the data (from VDB cursor reads) and this struct borrows it.
pub struct UnalignedColumns<'a> {
    /// Read name from VDB (`NAME`).
    pub name: &'a str,
    /// Read sequence.
    pub read: &'a str,
    /// Quality (raw Phred values, needs +33 for SAM output).
    pub quality: &'a [u8],
    /// Read group (`SPOT_GROUP`).
    pub spot_group: &'a str,
    /// Read type value (`READ_TYPE`).
    pub read_type: u8,
    /// Read filter value (`READ_FILTER`).
    pub read_filter: u8,
    /// Number of biological reads in the spot (for paired-end flag logic).
    pub num_bio_reads: u32,
    /// Which biological read this is (0 = first, 1 = second, etc.).
    pub bio_read_index: u32,
}

/// Read filter constants (from `INSDC:SRA:read_filter`).
#[allow(dead_code)]
pub const READ_FILTER_PASS: u8 = 0;
pub const READ_FILTER_REJECT: u8 = 1;
pub const READ_FILTER_CRITERIA: u8 = 2;
#[allow(dead_code)]
pub const READ_FILTER_REDACTED: u8 = 3;

/// Read type bitmask (from `INSDC:SRA:xread_type`).
pub const READ_TYPE_BIOLOGICAL: u8 = 1;
pub const READ_TYPE_REVERSE: u8 = 2;

/// SAM flag bit constants.
mod sam_flags {
    pub const PAIRED: u32 = 0x1;
    pub const PROPER_PAIR: u32 = 0x2;
    pub const UNMAPPED: u32 = 0x4;
    pub const MATE_UNMAPPED: u32 = 0x8;
    pub const REVERSE: u32 = 0x10;
    pub const MATE_REVERSE: u32 = 0x20;
    pub const FIRST_IN_PAIR: u32 = 0x40;
    pub const LAST_IN_PAIR: u32 = 0x80;
    pub const QC_FAIL: u32 = 0x200;
    pub const SUPPLEMENTARY_FILTER: u32 = 0x400;

    /// Mask of all paired-end related flags, cleared when mate has no alignment.
    pub const PAIRED_MASK: u32 =
        PAIRED | PROPER_PAIR | MATE_UNMAPPED | MATE_REVERSE | FIRST_IN_PAIR | LAST_IN_PAIR;
}

/// Options controlling SAM record formatting.
pub struct FormatOptions<'a> {
    /// Prepend this prefix to QNAME (e.g. "PRE.name").
    pub prefix: Option<&'a str>,
    /// Append spot group to QNAME (e.g. "name.RG1").
    pub spot_group_in_name: bool,
    /// Output alignment ID as XI:i tag.
    pub xi_tag: bool,
    /// Reverse unaligned reads according to read type.
    pub reverse_unaligned: bool,
}

/// Format a SAM line for an aligned record.
///
/// Writes a complete tab-delimited SAM line (with trailing newline) into `buf`.
#[allow(clippy::too_many_arguments)]
pub fn format_aligned_record(
    buf: &mut Vec<u8>,
    cols: &AlignedColumns,
    ref_name: &str,
    ref_pos: i32,
    mapq: i32,
    align_id: i64,
    mate_info: Option<&MateInfo>,
    opts: &FormatOptions<'_>,
) {
    buf.clear();

    // QNAME
    write_qname(buf, opts.prefix, &cols.seq_name, &cols.spot_group, '.', opts.spot_group_in_name);

    // FLAG — apply read filter to flags.
    let flags = apply_read_filter(cols.sam_flags, cols.read_filter);
    buf.push(b'\t');
    write_u32(buf, flags);

    // RNAME
    buf.push(b'\t');
    buf.extend_from_slice(ref_name.as_bytes());

    // POS (1-based)
    buf.push(b'\t');
    write_i32(buf, ref_pos + 1);

    // MAPQ
    buf.push(b'\t');
    write_i32(buf, mapq);

    // CIGAR
    buf.push(b'\t');
    if cols.cigar.is_empty() {
        buf.push(b'*');
    } else {
        buf.extend_from_slice(cols.cigar.as_bytes());
    }

    // RNEXT, PNEXT, TLEN — use mate cache if available, else column data.
    let (rnext, pnext, tlen) = if let Some(mate) = mate_info {
        (mate.ref_name.as_str(), mate.ref_pos, mate.tlen)
    } else {
        (cols.mate_ref_name.as_str(), cols.mate_ref_pos, cols.template_len)
    };

    buf.push(b'\t');
    if rnext.is_empty() || rnext == "*" {
        buf.push(b'*');
    } else if rnext == ref_name {
        buf.push(b'=');
    } else {
        buf.extend_from_slice(rnext.as_bytes());
    }

    // PNEXT (1-based, 0 if unavailable)
    buf.push(b'\t');
    if rnext.is_empty() || rnext == "*" {
        buf.push(b'0');
    } else {
        write_i32(buf, pnext + 1);
    }

    // TLEN
    buf.push(b'\t');
    write_i32(buf, tlen);

    // SEQ
    buf.push(b'\t');
    if cols.read.is_empty() {
        buf.push(b'*');
    } else {
        buf.extend_from_slice(cols.read.as_bytes());
    }

    // QUAL
    buf.push(b'\t');
    if cols.quality.is_empty() {
        buf.push(b'*');
    } else {
        buf.extend_from_slice(cols.quality.as_bytes());
    }

    // Optional tags.
    write_tag_u32(buf, b"NM", cols.edit_distance);

    if cols.alignment_count > 0 {
        write_tag_u32(buf, b"NH", cols.alignment_count as u32);
    }

    if !cols.spot_group.is_empty() {
        write_tag_str(buf, b"RG", &cols.spot_group);
    }

    if opts.xi_tag {
        write_tag_i64(buf, b"XI", align_id);
    }

    buf.push(b'\n');
}

/// Format a SAM line for an unaligned record.
pub fn format_unaligned_record(
    buf: &mut Vec<u8>,
    cols: &UnalignedColumns<'_>,
    opts: &FormatOptions<'_>,
) {
    buf.clear();

    // QNAME — unaligned uses '#' as spot group separator.
    write_qname(buf, opts.prefix, cols.name, cols.spot_group, '#', opts.spot_group_in_name);

    // FLAG
    let mut flags: u32 = sam_flags::UNMAPPED;
    if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
        flags |= sam_flags::REVERSE;
    }
    if cols.read_filter == READ_FILTER_REJECT {
        flags |= sam_flags::QC_FAIL;
    } else if cols.read_filter == READ_FILTER_CRITERIA {
        flags |= sam_flags::SUPPLEMENTARY_FILTER;
    }
    // Paired-end flags.
    if cols.num_bio_reads > 1 {
        flags |= sam_flags::PAIRED | sam_flags::MATE_UNMAPPED;
        if cols.bio_read_index == 0 {
            flags |= sam_flags::FIRST_IN_PAIR;
        }
        if cols.bio_read_index == cols.num_bio_reads - 1 {
            flags |= sam_flags::LAST_IN_PAIR;
        }
    }
    buf.push(b'\t');
    write_u32(buf, flags);

    // RNAME = *
    buf.extend_from_slice(b"\t*");

    // POS = 0
    buf.extend_from_slice(b"\t0");

    // MAPQ = 0
    buf.extend_from_slice(b"\t0");

    // CIGAR = *
    buf.extend_from_slice(b"\t*");

    // RNEXT = *
    buf.extend_from_slice(b"\t*");

    // PNEXT = 0
    buf.extend_from_slice(b"\t0");

    // TLEN = 0
    buf.extend_from_slice(b"\t0");

    // SEQ
    buf.push(b'\t');
    if cols.read.is_empty() {
        buf.push(b'*');
    } else {
        // Reverse complement if needed.
        if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
            write_reverse_complement(buf, cols.read.as_bytes());
        } else {
            buf.extend_from_slice(cols.read.as_bytes());
        }
    }

    // QUAL (raw phred → phred+33)
    buf.push(b'\t');
    if cols.quality.is_empty() {
        buf.push(b'*');
    } else if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
        // Reverse the quality scores.
        for &q in cols.quality.iter().rev() {
            buf.push(q + 33);
        }
    } else {
        for &q in cols.quality {
            buf.push(q + 33);
        }
    }

    // RG tag
    if !cols.spot_group.is_empty() {
        write_tag_str(buf, b"RG", cols.spot_group);
    }

    buf.push(b'\n');
}

/// Write QNAME: `[prefix.]{name}[{sep}{spot_group}]`
fn write_qname(
    buf: &mut Vec<u8>,
    prefix: Option<&str>,
    name: &str,
    spot_group: &str,
    sep: char,
    include_spot_group: bool,
) {
    if let Some(pfx) = prefix {
        buf.extend_from_slice(pfx.as_bytes());
        buf.push(b'.');
    }
    buf.extend_from_slice(name.as_bytes());
    if include_spot_group && !spot_group.is_empty() {
        buf.push(sep as u8);
        buf.extend_from_slice(spot_group.as_bytes());
    }
}

/// Apply read filter to SAM flags.
fn apply_read_filter(mut flags: u32, read_filter: Option<u8>) -> u32 {
    if let Some(filt) = read_filter {
        flags &= !sam_flags::QC_FAIL;
        if filt == READ_FILTER_REJECT {
            flags |= sam_flags::QC_FAIL;
        } else if filt == READ_FILTER_CRITERIA {
            flags |= sam_flags::SUPPLEMENTARY_FILTER;
        }
    }
    flags
}

/// Write the reverse complement of a DNA sequence.
fn write_reverse_complement(buf: &mut Vec<u8>, seq: &[u8]) {
    for &b in seq.iter().rev() {
        buf.push(complement(b));
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Write a u32 as decimal ASCII.
fn write_u32(buf: &mut Vec<u8>, val: u32) {
    let _ = write!(buf, "{val}");
}

/// Write an i32 as decimal ASCII.
fn write_i32(buf: &mut Vec<u8>, val: i32) {
    let _ = write!(buf, "{val}");
}

/// Write a SAM tag with integer value: `\tXX:i:val`
fn write_tag_u32(buf: &mut Vec<u8>, tag: &[u8; 2], val: u32) {
    buf.push(b'\t');
    buf.extend_from_slice(tag);
    buf.extend_from_slice(b":i:");
    write_u32(buf, val);
}

/// Write a SAM tag with i64 value: `\tXX:i:val`
fn write_tag_i64(buf: &mut Vec<u8>, tag: &[u8; 2], val: i64) {
    buf.push(b'\t');
    buf.extend_from_slice(tag);
    buf.extend_from_slice(b":i:");
    let _ = write!(buf, "{val}");
}

/// Write a SAM tag with string value: `\tXX:Z:val`
fn write_tag_str(buf: &mut Vec<u8>, tag: &[u8; 2], val: &str) {
    buf.push(b'\t');
    buf.extend_from_slice(tag);
    buf.extend_from_slice(b":Z:");
    buf.extend_from_slice(val.as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_opts() -> FormatOptions<'static> {
        FormatOptions {
            prefix: None,
            spot_group_in_name: false,
            xi_tag: false,
            reverse_unaligned: false,
        }
    }

    fn default_aligned_cols() -> AlignedColumns {
        AlignedColumns {
            seq_name: "read1".to_string(),
            sam_flags: 99,
            cigar: "50M".to_string(),
            mate_align_id: 2,
            mate_ref_name: "chr1".to_string(),
            mate_ref_pos: 500,
            template_len: 300,
            read: "ACGTACGT".to_string(),
            quality: "IIIIIIII".to_string(),
            edit_distance: 1,
            spot_group: "RG1".to_string(),
            alignment_count: 1,
            read_filter: None,
        }
    }

    #[test]
    fn test_aligned_record_basic() {
        let cols = default_aligned_cols();
        let opts = default_opts();
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();

        assert_eq!(fields[0], "read1"); // QNAME
        assert_eq!(fields[1], "99"); // FLAG
        assert_eq!(fields[2], "chr1"); // RNAME
        assert_eq!(fields[3], "101"); // POS (1-based)
        assert_eq!(fields[4], "60"); // MAPQ
        assert_eq!(fields[5], "50M"); // CIGAR
        assert_eq!(fields[6], "="); // RNEXT (same as RNAME)
        assert_eq!(fields[7], "501"); // PNEXT (1-based)
        assert_eq!(fields[8], "300"); // TLEN
        assert_eq!(fields[9], "ACGTACGT"); // SEQ
        assert_eq!(fields[10], "IIIIIIII"); // QUAL
    }

    #[test]
    fn test_aligned_record_with_prefix_and_spot_group() {
        let cols = default_aligned_cols();
        let opts = FormatOptions {
            prefix: Some("PRE"),
            spot_group_in_name: true,
            xi_tag: false,
            reverse_unaligned: false,
        };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[0], "PRE.read1.RG1");
    }

    #[test]
    fn test_aligned_record_with_mate_cache() {
        let cols = default_aligned_cols();
        let opts = default_opts();
        let mate = MateInfo { ref_name: "chr2".to_string(), ref_pos: 200, flags: 147, tlen: -300 };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, Some(&mate), &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[6], "chr2"); // RNEXT from mate cache
        assert_eq!(fields[7], "201"); // PNEXT from mate cache (1-based)
        assert_eq!(fields[8], "-300"); // TLEN from mate cache
    }

    #[test]
    fn test_aligned_record_tags() {
        let cols = default_aligned_cols();
        let opts = FormatOptions { xi_tag: true, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        assert!(line.contains("NM:i:1"));
        assert!(line.contains("NH:i:1"));
        assert!(line.contains("RG:Z:RG1"));
        assert!(line.contains("XI:i:42"));
    }

    #[test]
    fn test_aligned_record_read_filter_reject() {
        let mut cols = default_aligned_cols();
        cols.sam_flags = 99;
        cols.read_filter = Some(READ_FILTER_REJECT);
        let opts = default_opts();
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let flags: u32 = fields[1].parse().unwrap();
        assert!(flags & sam_flags::QC_FAIL != 0, "QC fail flag should be set");
    }

    #[test]
    fn test_unaligned_record_basic() {
        let cols = UnalignedColumns {
            name: "spot1",
            read: "ACGT",
            quality: &[30, 30, 30, 30],
            spot_group: "",
            read_type: READ_TYPE_BIOLOGICAL,
            read_filter: READ_FILTER_PASS,
            num_bio_reads: 1,
            bio_read_index: 0,
        };
        let opts = default_opts();
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[0], "spot1"); // QNAME
        assert_eq!(fields[1], "4"); // FLAG = unmapped
        assert_eq!(fields[2], "*"); // RNAME
        assert_eq!(fields[9], "ACGT"); // SEQ
        // Quality should be phred+33: 30+33 = 63 = '?'
        assert_eq!(fields[10], "????");
    }

    #[test]
    fn test_unaligned_record_paired() {
        let cols = UnalignedColumns {
            name: "spot1",
            read: "ACGT",
            quality: &[30, 30, 30, 30],
            spot_group: "",
            read_type: READ_TYPE_BIOLOGICAL,
            read_filter: READ_FILTER_PASS,
            num_bio_reads: 2,
            bio_read_index: 0,
        };
        let opts = default_opts();
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let flags: u32 = fields[1].parse().unwrap();
        assert!(flags & sam_flags::PAIRED != 0, "paired flag");
        assert!(flags & sam_flags::UNMAPPED != 0, "unmapped flag");
        assert!(flags & sam_flags::MATE_UNMAPPED != 0, "mate unmapped flag");
        assert!(flags & sam_flags::FIRST_IN_PAIR != 0, "first in pair flag");
    }

    #[test]
    fn test_unaligned_record_reverse() {
        let cols = UnalignedColumns {
            name: "spot1",
            read: "ACGT",
            quality: &[10, 20, 30, 40],
            spot_group: "",
            read_type: READ_TYPE_BIOLOGICAL | READ_TYPE_REVERSE,
            read_filter: READ_FILTER_PASS,
            num_bio_reads: 1,
            bio_read_index: 0,
        };
        let opts = FormatOptions { reverse_unaligned: true, ..default_opts() };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let flags: u32 = fields[1].parse().unwrap();
        assert!(flags & sam_flags::REVERSE != 0, "reverse flag");
        assert_eq!(fields[9], "ACGT"); // reverse complement of ACGT
        // Quality reversed: 40+33, 30+33, 20+33, 10+33 = 'I', '?', '5', '+'
        assert_eq!(fields[10], "I?5+");
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_write_qname_variants() {
        let mut buf = Vec::new();

        // Plain name.
        write_qname(&mut buf, None, "read1", "", '.', false);
        assert_eq!(&buf, b"read1");

        // With prefix.
        buf.clear();
        write_qname(&mut buf, Some("PRE"), "read1", "", '.', false);
        assert_eq!(&buf, b"PRE.read1");

        // With spot group.
        buf.clear();
        write_qname(&mut buf, None, "read1", "RG1", '.', true);
        assert_eq!(&buf, b"read1.RG1");

        // With prefix and spot group.
        buf.clear();
        write_qname(&mut buf, Some("X"), "name", "grp", '#', true);
        assert_eq!(&buf, b"X.name#grp");
    }

    #[test]
    fn test_apply_read_filter() {
        // No filter.
        assert_eq!(apply_read_filter(99, None), 99);

        // PASS filter — clears existing QC_FAIL.
        assert_eq!(apply_read_filter(99 | sam_flags::QC_FAIL, Some(READ_FILTER_PASS)), 99);

        // REJECT filter.
        assert_eq!(apply_read_filter(99, Some(READ_FILTER_REJECT)), 99 | sam_flags::QC_FAIL);

        // CRITERIA filter.
        assert_eq!(
            apply_read_filter(99, Some(READ_FILTER_CRITERIA)),
            99 | sam_flags::SUPPLEMENTARY_FILTER
        );
    }

    #[test]
    fn test_strip_paired_flags() {
        let mut cols = default_aligned_cols();

        // FLAG=73 (paired + mate_unmapped + first_in_pair) → 0
        cols.sam_flags = 0x49;
        cols.strip_paired_flags();
        assert_eq!(cols.sam_flags, 0);

        // FLAG=89 (paired + mate_unmapped + reverse + first_in_pair) → 16 (reverse only)
        cols.sam_flags = 0x59;
        cols.strip_paired_flags();
        assert_eq!(cols.sam_flags, 0x10);

        // FLAG=585 (paired + mate_unmapped + first_in_pair + QC_FAIL) → 512 (QC_FAIL only)
        cols.sam_flags = 0x249;
        cols.strip_paired_flags();
        assert_eq!(cols.sam_flags, 0x200);

        // FLAG=99 (paired + proper_pair + mate_reverse + first_in_pair) → should strip
        cols.sam_flags = 99;
        cols.strip_paired_flags();
        assert_eq!(cols.sam_flags, 99 & !(0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80));
    }

    #[test]
    fn test_unmapped_rnext_pnext() {
        let mut cols = default_aligned_cols();
        cols.mate_ref_name = "*".to_string();
        cols.mate_ref_pos = 0;
        let opts = default_opts();
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[6], "*"); // RNEXT
        assert_eq!(fields[7], "0"); // PNEXT
    }
}
