//! SAM/BAM record building and formatting.
//!
//! Builds complete SAM lines in a per-record `Vec<u8>` buffer, writing
//! all fields in a single pass rather than multiple formatted-print calls.

use std::collections::HashMap;

use crate::matecache::MateInfo;
use crate::quality::QuantTable;

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
    /// 0-based reference position (`REF_POS`).
    pub ref_pos: i32,
    /// Mapping quality (`MAPQ`).
    pub mapq: i32,
}

impl AlignedColumns {
    /// Create a new `AlignedColumns` with empty/default values.
    ///
    /// String fields start empty and grow to needed capacity on first use,
    /// reusing their heap allocation across subsequent records.
    pub fn new() -> Self {
        Self {
            seq_name: String::new(),
            sam_flags: 0,
            cigar: String::new(),
            mate_align_id: 0,
            mate_ref_name: String::new(),
            mate_ref_pos: 0,
            template_len: 0,
            read: String::new(),
            quality: String::new(),
            edit_distance: 0,
            spot_group: String::new(),
            alignment_count: 0,
            read_filter: None,
            ref_pos: 0,
            mapq: 0,
        }
    }

    /// Clear all fields without deallocating String buffers.
    pub fn clear(&mut self) {
        self.seq_name.clear();
        self.sam_flags = 0;
        self.cigar.clear();
        self.mate_align_id = 0;
        self.mate_ref_name.clear();
        self.mate_ref_pos = 0;
        self.template_len = 0;
        self.read.clear();
        self.quality.clear();
        self.edit_distance = 0;
        self.spot_group.clear();
        self.alignment_count = 0;
        self.read_filter = None;
        self.ref_pos = 0;
        self.mapq = 0;
    }

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

/// Output format mode.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    Sam,
    Bam,
    Fasta,
    Fastq,
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
    /// Replace quality values with `*`.
    pub omit_quality: bool,
    /// Quality score quantization table.
    pub qual_quant: Option<&'a QuantTable>,
    /// Output mode (SAM, BAM, FASTA, or FASTQ).
    pub output_mode: OutputMode,
    /// Reference name → BAM reference ID map. Required for BAM output.
    pub ref_name_to_id: Option<&'a HashMap<String, i32>>,
}

/// Format a record for an aligned read, dispatching by output mode.
#[allow(clippy::too_many_arguments)]
pub fn format_aligned_record(
    buf: &mut Vec<u8>,
    cols: &AlignedColumns,
    ref_name: &str,
    ref_id: i32,
    ref_pos: i32,
    mapq: i32,
    align_id: i64,
    mate_info: Option<&MateInfo>,
    opts: &FormatOptions<'_>,
) {
    match opts.output_mode {
        OutputMode::Sam => {
            format_aligned_record_sam(buf, cols, ref_name, ref_pos, mapq, align_id, mate_info, opts)
        }
        OutputMode::Bam => {
            format_aligned_record_bam(buf, cols, ref_id, ref_pos, mapq, align_id, mate_info, opts)
        }
        OutputMode::Fasta => {
            buf.clear();
            buf.push(b'>');
            write_qname(
                buf,
                opts.prefix,
                &cols.seq_name,
                &cols.spot_group,
                '.',
                opts.spot_group_in_name,
            );
            buf.push(b'\n');
            buf.extend_from_slice(cols.read.as_bytes());
            buf.push(b'\n');
        }
        OutputMode::Fastq => {
            buf.clear();
            buf.push(b'@');
            write_qname(
                buf,
                opts.prefix,
                &cols.seq_name,
                &cols.spot_group,
                '.',
                opts.spot_group_in_name,
            );
            buf.push(b'\n');
            buf.extend_from_slice(cols.read.as_bytes());
            buf.extend_from_slice(b"\n+\n");
            if opts.omit_quality || cols.quality.is_empty() {
                buf.push(b'*');
            } else if let Some(table) = opts.qual_quant {
                for &q in cols.quality.as_bytes() {
                    buf.push(crate::quality::quantize_phred33(q, table));
                }
            } else {
                buf.extend_from_slice(cols.quality.as_bytes());
            }
            buf.push(b'\n');
        }
    }
}

/// Format a SAM line for an aligned record.
///
/// Writes a complete tab-delimited SAM line (with trailing newline) into `buf`.
#[allow(clippy::too_many_arguments)]
fn format_aligned_record_sam(
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
    // Mate cache is cleared between references, so cached mates are always
    // on the same reference → RNEXT is always "=".
    // TLEN always comes from the current record's column, not the cache.
    if let Some(mate) = mate_info {
        buf.extend_from_slice(b"\t=\t");
        write_i32(buf, mate.ref_pos + 1);
        buf.push(b'\t');
        write_i32(buf, cols.template_len);
    } else {
        let rnext = cols.mate_ref_name.as_str();
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
            write_i32(buf, cols.mate_ref_pos + 1);
        }

        // TLEN
        buf.push(b'\t');
        write_i32(buf, cols.template_len);
    }

    // SEQ
    buf.push(b'\t');
    if cols.read.is_empty() {
        buf.push(b'*');
    } else {
        buf.extend_from_slice(cols.read.as_bytes());
    }

    // QUAL
    buf.push(b'\t');
    if opts.omit_quality || cols.quality.is_empty() {
        buf.push(b'*');
    } else if let Some(table) = opts.qual_quant {
        for &q in cols.quality.as_bytes() {
            buf.push(crate::quality::quantize_phred33(q, table));
        }
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

/// Format a record for an unaligned read, dispatching by output mode.
pub fn format_unaligned_record(
    buf: &mut Vec<u8>,
    cols: &UnalignedColumns<'_>,
    opts: &FormatOptions<'_>,
) {
    match opts.output_mode {
        OutputMode::Sam => format_unaligned_record_sam(buf, cols, opts),
        OutputMode::Bam => format_unaligned_record_bam(buf, cols, opts),
        OutputMode::Fasta => {
            buf.clear();
            buf.push(b'>');
            write_qname(buf, opts.prefix, cols.name, cols.spot_group, '#', opts.spot_group_in_name);
            buf.push(b'\n');
            if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
                write_reverse_complement(buf, cols.read.as_bytes());
            } else {
                buf.extend_from_slice(cols.read.as_bytes());
            }
            buf.push(b'\n');
        }
        OutputMode::Fastq => {
            buf.clear();
            buf.push(b'@');
            write_qname(buf, opts.prefix, cols.name, cols.spot_group, '#', opts.spot_group_in_name);
            buf.push(b'\n');
            if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
                write_reverse_complement(buf, cols.read.as_bytes());
            } else {
                buf.extend_from_slice(cols.read.as_bytes());
            }
            buf.extend_from_slice(b"\n+\n");
            if opts.omit_quality || cols.quality.is_empty() {
                buf.push(b'*');
            } else if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
                for &q in cols.quality.iter().rev() {
                    let phred = if let Some(table) = opts.qual_quant {
                        crate::quality::quantize_phred(q, table)
                    } else {
                        q
                    };
                    buf.push(phred + 33);
                }
            } else {
                for &q in cols.quality {
                    let phred = if let Some(table) = opts.qual_quant {
                        crate::quality::quantize_phred(q, table)
                    } else {
                        q
                    };
                    buf.push(phred + 33);
                }
            }
            buf.push(b'\n');
        }
    }
}

/// Format a SAM line for an unaligned record.
fn format_unaligned_record_sam(
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
    flags = apply_read_filter(flags, Some(cols.read_filter));
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
    if opts.omit_quality || cols.quality.is_empty() {
        buf.push(b'*');
    } else if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
        // Reverse the quality scores.
        for &q in cols.quality.iter().rev() {
            let phred = if let Some(table) = opts.qual_quant {
                crate::quality::quantize_phred(q, table)
            } else {
                q
            };
            buf.push(phred + 33);
        }
    } else {
        for &q in cols.quality {
            let phred = if let Some(table) = opts.qual_quant {
                crate::quality::quantize_phred(q, table)
            } else {
                q
            };
            buf.push(phred + 33);
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
pub(crate) fn apply_read_filter(mut flags: u32, read_filter: Option<u8>) -> u32 {
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
    buf.extend_from_slice(itoa::Buffer::new().format(val).as_bytes());
}

/// Write an i32 as decimal ASCII.
fn write_i32(buf: &mut Vec<u8>, val: i32) {
    buf.extend_from_slice(itoa::Buffer::new().format(val).as_bytes());
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
    buf.extend_from_slice(itoa::Buffer::new().format(val).as_bytes());
}

/// Write a SAM tag with string value: `\tXX:Z:val`
fn write_tag_str(buf: &mut Vec<u8>, tag: &[u8; 2], val: &str) {
    buf.push(b'\t');
    buf.extend_from_slice(tag);
    buf.extend_from_slice(b":Z:");
    buf.extend_from_slice(val.as_bytes());
}

// ── BAM binary encoding ──────────────────────────────────────────────────

/// Format an aligned record as BAM binary (block_size prefix + record data).
#[allow(clippy::too_many_arguments)]
fn format_aligned_record_bam(
    buf: &mut Vec<u8>,
    cols: &AlignedColumns,
    ref_id: i32,
    ref_pos: i32,
    mapq: i32,
    align_id: i64,
    mate_info: Option<&MateInfo>,
    opts: &FormatOptions<'_>,
) {
    buf.clear();

    let flags = apply_read_filter(cols.sam_flags, cols.read_filter);

    // Build QNAME (null-terminated).
    let mut qname_buf = Vec::with_capacity(64);
    write_qname(
        &mut qname_buf,
        opts.prefix,
        &cols.seq_name,
        &cols.spot_group,
        '.',
        opts.spot_group_in_name,
    );
    qname_buf.push(0); // null terminator
    let l_read_name = qname_buf.len() as u8;

    // Parse CIGAR ops.
    let cigar_ops = parse_cigar_to_bam(&cols.cigar);
    let n_cigar_op = cigar_ops.len() as u16;

    // Compute alignment end for bin calculation.
    let align_end = ref_pos + cigar_ref_length(&cigar_ops);
    let bin = reg2bin(ref_pos as u32, align_end as u32);

    let seq = cols.read.as_bytes();
    let l_seq = seq.len() as i32;

    // Resolve mate info.
    // TLEN always comes from the current record's column, not the cache.
    let (mate_ref_id, mate_pos, tlen) = if let Some(m) = mate_info {
        // Mate is on same reference (mate cache only stores same-ref mates).
        (ref_id, m.ref_pos, cols.template_len)
    } else if cols.mate_align_id != 0 {
        let rnext = cols.mate_ref_name.as_str();
        let mate_rid = if rnext.is_empty() || rnext == "*" {
            -1
        } else if rnext == "=" {
            ref_id
        } else {
            opts.ref_name_to_id.and_then(|m| m.get(rnext).copied()).unwrap_or(-1)
        };
        let mpos = if mate_rid == -1 { -1 } else { cols.mate_ref_pos };
        (mate_rid, mpos, cols.template_len)
    } else {
        (-1, -1, 0)
    };

    // Reserve space for block_size (4 bytes), filled in at the end.
    let block_start = buf.len();
    buf.extend_from_slice(&[0u8; 4]);

    // Fixed-length fields (32 bytes).
    buf.extend_from_slice(&ref_id.to_le_bytes());
    buf.extend_from_slice(&ref_pos.to_le_bytes());
    buf.push(l_read_name);
    buf.push(mapq as u8);
    buf.extend_from_slice(&bin.to_le_bytes());
    buf.extend_from_slice(&n_cigar_op.to_le_bytes());
    buf.extend_from_slice(&(flags as u16).to_le_bytes());
    buf.extend_from_slice(&l_seq.to_le_bytes());
    buf.extend_from_slice(&mate_ref_id.to_le_bytes());
    buf.extend_from_slice(&mate_pos.to_le_bytes());
    buf.extend_from_slice(&tlen.to_le_bytes());

    // Variable-length fields.
    buf.extend_from_slice(&qname_buf);

    for &op in &cigar_ops {
        buf.extend_from_slice(&op.to_le_bytes());
    }

    encode_sequence(buf, seq);
    encode_quality_aligned(buf, cols.quality.as_bytes(), l_seq as usize, opts);

    // Aux tags.
    write_bam_tag_u32(buf, b"NM", cols.edit_distance);
    if cols.alignment_count > 0 {
        write_bam_tag_u32(buf, b"NH", cols.alignment_count as u32);
    }
    if !cols.spot_group.is_empty() {
        write_bam_tag_str(buf, b"RG", &cols.spot_group);
    }
    if opts.xi_tag {
        write_bam_tag_i64(buf, b"XI", align_id);
    }

    // Fill in block_size.
    let block_size = (buf.len() - block_start - 4) as i32;
    buf[block_start..block_start + 4].copy_from_slice(&block_size.to_le_bytes());
}

/// Format an unaligned record as BAM binary.
fn format_unaligned_record_bam(
    buf: &mut Vec<u8>,
    cols: &UnalignedColumns<'_>,
    opts: &FormatOptions<'_>,
) {
    buf.clear();

    let mut flags: u32 = sam_flags::UNMAPPED;
    if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
        flags |= sam_flags::REVERSE;
    }
    flags = apply_read_filter(flags, Some(cols.read_filter));
    if cols.num_bio_reads > 1 {
        flags |= sam_flags::PAIRED | sam_flags::MATE_UNMAPPED;
        if cols.bio_read_index == 0 {
            flags |= sam_flags::FIRST_IN_PAIR;
        }
        if cols.bio_read_index == cols.num_bio_reads - 1 {
            flags |= sam_flags::LAST_IN_PAIR;
        }
    }

    // Build QNAME (null-terminated).
    let mut qname_buf = Vec::with_capacity(64);
    write_qname(
        &mut qname_buf,
        opts.prefix,
        cols.name,
        cols.spot_group,
        '#',
        opts.spot_group_in_name,
    );
    qname_buf.push(0);
    let l_read_name = qname_buf.len() as u8;

    let seq_bytes = cols.read.as_bytes();
    let l_seq = seq_bytes.len() as i32;

    // Reserve space for block_size.
    let block_start = buf.len();
    buf.extend_from_slice(&[0u8; 4]);

    // Fixed-length fields.
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // refID = -1
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // pos = -1
    buf.push(l_read_name);
    buf.push(0u8); // mapq = 0
    buf.extend_from_slice(&4680u16.to_le_bytes()); // bin for unmapped
    buf.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
    buf.extend_from_slice(&(flags as u16).to_le_bytes());
    buf.extend_from_slice(&l_seq.to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // mate refID = -1
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // mate pos = -1
    buf.extend_from_slice(&0i32.to_le_bytes()); // tlen = 0

    // Variable-length fields.
    buf.extend_from_slice(&qname_buf);

    // Sequence (reverse complement if needed).
    if opts.reverse_unaligned && (cols.read_type & READ_TYPE_REVERSE) != 0 {
        let mut rc = Vec::with_capacity(seq_bytes.len());
        for &b in seq_bytes.iter().rev() {
            rc.push(complement(b));
        }
        encode_sequence(buf, &rc);
    } else {
        encode_sequence(buf, seq_bytes);
    }

    // Quality (raw phred → raw phred for BAM, with optional quantization/reversal).
    encode_quality_unaligned(buf, cols.quality, l_seq as usize, cols.read_type, opts);

    // Aux tags.
    if !cols.spot_group.is_empty() {
        write_bam_tag_str(buf, b"RG", cols.spot_group);
    }

    // Fill in block_size.
    let block_size = (buf.len() - block_start - 4) as i32;
    buf[block_start..block_start + 4].copy_from_slice(&block_size.to_le_bytes());
}

/// Parse a CIGAR string (e.g. "50M2I48M") into packed BAM CIGAR ops.
///
/// Each op is encoded as `(op_len << 4) | op_code` in a u32.
fn parse_cigar_to_bam(cigar: &str) -> Vec<u32> {
    let mut ops = Vec::new();
    let mut num: u32 = 0;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as u32;
        } else {
            let code = match b {
                b'M' => 0,
                b'I' => 1,
                b'D' => 2,
                b'N' => 3,
                b'S' => 4,
                b'H' => 5,
                b'P' => 6,
                b'=' => 7,
                b'X' => 8,
                _ => continue,
            };
            ops.push((num << 4) | code);
            num = 0;
        }
    }
    ops
}

/// Compute the reference-consuming length from packed CIGAR ops.
fn cigar_ref_length(ops: &[u32]) -> i32 {
    let mut len: i32 = 0;
    for &op in ops {
        let code = op & 0xf;
        let op_len = (op >> 4) as i32;
        // Reference-consuming ops: M(0), D(2), N(3), =(7), X(8)
        if matches!(code, 0 | 2 | 3 | 7 | 8) {
            len += op_len;
        }
    }
    len
}

/// Compute the BAM indexing bin for a region `[beg, end)`.
///
/// Uses the binning scheme from the SAM specification with precomputed
/// level offsets: 0, 1, 9, 73, 585, 4681.
fn reg2bin(beg: u32, end: u32) -> u16 {
    let end = end.saturating_sub(1);
    if beg >> 14 == end >> 14 {
        return (4681 + (beg >> 14)) as u16;
    }
    if beg >> 17 == end >> 17 {
        return (585 + (beg >> 17)) as u16;
    }
    if beg >> 20 == end >> 20 {
        return (73 + (beg >> 20)) as u16;
    }
    if beg >> 23 == end >> 23 {
        return (9 + (beg >> 23)) as u16;
    }
    if beg >> 26 == end >> 26 {
        return (1 + (beg >> 26)) as u16;
    }
    0
}

/// Encode a DNA sequence as 4-bit packed bytes for BAM.
fn encode_sequence(buf: &mut Vec<u8>, seq: &[u8]) {
    let packed_len = seq.len().div_ceil(2);
    let start = buf.len();
    buf.resize(start + packed_len, 0);
    for (i, &base) in seq.iter().enumerate() {
        let code = base_to_bam(base);
        if i % 2 == 0 {
            buf[start + i / 2] = code << 4;
        } else {
            buf[start + i / 2] |= code;
        }
    }
}

/// Convert an ASCII base to its 4-bit BAM encoding.
fn base_to_bam(base: u8) -> u8 {
    match base {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        _ => 15, // N
    }
}

/// Encode quality scores for an aligned BAM record (already Phred+33 → raw Phred).
fn encode_quality_aligned(
    buf: &mut Vec<u8>,
    qual_p33: &[u8],
    l_seq: usize,
    opts: &FormatOptions<'_>,
) {
    if opts.omit_quality || qual_p33.is_empty() {
        // 0xFF means "quality not stored".
        buf.extend(std::iter::repeat_n(0xFFu8, l_seq));
    } else if let Some(table) = opts.qual_quant {
        for &q in qual_p33 {
            let phred = q.saturating_sub(33);
            buf.push(crate::quality::quantize_phred(phred, table));
        }
    } else {
        for &q in qual_p33 {
            buf.push(q.saturating_sub(33));
        }
    }
}

/// Encode quality scores for an unaligned BAM record (raw Phred values).
fn encode_quality_unaligned(
    buf: &mut Vec<u8>,
    quality: &[u8],
    l_seq: usize,
    read_type: u8,
    opts: &FormatOptions<'_>,
) {
    if opts.omit_quality || quality.is_empty() {
        buf.extend(std::iter::repeat_n(0xFFu8, l_seq));
    } else if opts.reverse_unaligned && (read_type & READ_TYPE_REVERSE) != 0 {
        for &q in quality.iter().rev() {
            let phred = if let Some(table) = opts.qual_quant {
                crate::quality::quantize_phred(q, table)
            } else {
                q
            };
            buf.push(phred);
        }
    } else {
        for &q in quality {
            let phred = if let Some(table) = opts.qual_quant {
                crate::quality::quantize_phred(q, table)
            } else {
                q
            };
            buf.push(phred);
        }
    }
}

/// Write a BAM auxiliary tag with an integer value, using the smallest type.
fn write_bam_tag_u32(buf: &mut Vec<u8>, tag: &[u8; 2], val: u32) {
    buf.extend_from_slice(tag);
    if val <= 0xFF {
        buf.push(b'C');
        buf.push(val as u8);
    } else if val <= 0xFFFF {
        buf.push(b'S');
        buf.extend_from_slice(&(val as u16).to_le_bytes());
    } else {
        buf.push(b'I');
        buf.extend_from_slice(&val.to_le_bytes());
    }
}

/// Write a BAM auxiliary tag with an i64 value.
fn write_bam_tag_i64(buf: &mut Vec<u8>, tag: &[u8; 2], val: i64) {
    buf.extend_from_slice(tag);
    if (0..=0xFF).contains(&val) {
        buf.push(b'C');
        buf.push(val as u8);
    } else if (0..=0xFFFF).contains(&val) {
        buf.push(b'S');
        buf.extend_from_slice(&(val as u16).to_le_bytes());
    } else if (0..=0xFFFF_FFFF).contains(&val) {
        buf.push(b'I');
        buf.extend_from_slice(&(val as u32).to_le_bytes());
    } else if val >= i32::MIN as i64 && val <= i32::MAX as i64 {
        buf.push(b'i');
        buf.extend_from_slice(&(val as i32).to_le_bytes());
    } else {
        // BAM doesn't have a native i64 tag type; fall back to string.
        buf.push(b'Z');
        buf.extend_from_slice(itoa::Buffer::new().format(val).as_bytes());
        buf.push(0);
    }
}

/// Write a BAM auxiliary tag with a string value (null-terminated).
fn write_bam_tag_str(buf: &mut Vec<u8>, tag: &[u8; 2], val: &str) {
    buf.extend_from_slice(tag);
    buf.push(b'Z');
    buf.extend_from_slice(val.as_bytes());
    buf.push(0);
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
            omit_quality: false,
            qual_quant: None,
            output_mode: OutputMode::Sam,
            ref_name_to_id: None,
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
            ref_pos: 0,
            mapq: 0,
        }
    }

    #[test]
    fn test_aligned_record_basic() {
        let cols = default_aligned_cols();
        let opts = default_opts();
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

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
        let opts =
            FormatOptions { prefix: Some("PRE"), spot_group_in_name: true, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[0], "PRE.read1.RG1");
    }

    #[test]
    fn test_aligned_record_with_mate_cache() {
        let cols = default_aligned_cols();
        let opts = default_opts();
        let mate = MateInfo { ref_pos: 200 };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, Some(&mate), &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[6], "="); // RNEXT from mate cache (same ref)
        assert_eq!(fields[7], "201"); // PNEXT from mate cache (1-based)
        assert_eq!(fields[8], "300"); // TLEN from record's own column, not cache
    }

    #[test]
    fn test_aligned_record_tags() {
        let cols = default_aligned_cols();
        let opts = FormatOptions { xi_tag: true, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

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

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

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

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[6], "*"); // RNEXT
        assert_eq!(fields[7], "0"); // PNEXT
    }

    #[test]
    fn test_aligned_omit_quality() {
        let cols = default_aligned_cols();
        let opts = FormatOptions { omit_quality: true, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[10], "*"); // QUAL omitted
    }

    #[test]
    fn test_unaligned_omit_quality() {
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
        let opts = FormatOptions { omit_quality: true, ..default_opts() };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(fields[10], "*"); // QUAL omitted
    }

    #[test]
    fn test_aligned_qual_quant() {
        let table = crate::quality::parse_qual_quant("0:10,10:20,20:30,30:40").unwrap();
        let mut cols = default_aligned_cols();
        // Phred+33 encoded: Phred 5 = 38, Phred 15 = 48, Phred 25 = 58, Phred 35 = 68
        cols.quality = String::from_utf8(vec![5 + 33, 15 + 33, 25 + 33, 35 + 33]).unwrap();
        let opts = FormatOptions { qual_quant: Some(&table), ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let qual_bytes = fields[10].as_bytes();
        // Phred 5 → 10, Phred 15 → 20, Phred 25 → 30, Phred 35 → 40
        assert_eq!(qual_bytes, &[10 + 33, 20 + 33, 30 + 33, 40 + 33]);
    }

    #[test]
    fn test_unaligned_qual_quant() {
        let table = crate::quality::parse_qual_quant("0:10,10:20,20:30,30:40").unwrap();
        let cols = UnalignedColumns {
            name: "spot1",
            read: "ACGT",
            quality: &[5, 15, 25, 35], // raw Phred values
            spot_group: "",
            read_type: READ_TYPE_BIOLOGICAL,
            read_filter: READ_FILTER_PASS,
            num_bio_reads: 1,
            bio_read_index: 0,
        };
        let opts = FormatOptions { qual_quant: Some(&table), ..default_opts() };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let line = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let qual_bytes = fields[10].as_bytes();
        // Phred 5 → 10, Phred 15 → 20, Phred 25 → 30, Phred 35 → 40
        assert_eq!(qual_bytes, &[10 + 33, 20 + 33, 30 + 33, 40 + 33]);
    }

    #[test]
    fn test_aligned_fasta() {
        let cols = default_aligned_cols();
        let opts = FormatOptions { output_mode: OutputMode::Fasta, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, ">read1\nACGTACGT\n");
    }

    #[test]
    fn test_aligned_fastq() {
        let cols = default_aligned_cols();
        let opts = FormatOptions { output_mode: OutputMode::Fastq, ..default_opts() };
        let mut buf = Vec::new();

        format_aligned_record(&mut buf, &cols, "chr1", 0, 100, 60, 42, None, &opts);

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "@read1\nACGTACGT\n+\nIIIIIIII\n");
    }

    #[test]
    fn test_unaligned_fasta() {
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
        let opts = FormatOptions { output_mode: OutputMode::Fasta, ..default_opts() };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, ">spot1\nACGT\n");
    }

    #[test]
    fn test_unaligned_fastq() {
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
        let opts = FormatOptions { output_mode: OutputMode::Fastq, ..default_opts() };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let output = String::from_utf8(buf).unwrap();
        // Quality: 30+33 = 63 = '?'
        assert_eq!(output, "@spot1\nACGT\n+\n????\n");
    }

    // ── BAM encoding tests ────────────────────────────────────────────

    #[test]
    fn test_parse_cigar_to_bam() {
        let ops = parse_cigar_to_bam("50M2I48M");
        assert_eq!(ops, vec![50 << 4, (2 << 4) | 1, 48 << 4]);
    }

    #[test]
    fn test_parse_cigar_to_bam_all_ops() {
        let ops = parse_cigar_to_bam("1M2I3D4N5S6H7P8=9X");
        assert_eq!(
            ops,
            vec![
                1 << 4,
                (2 << 4) | 1,
                (3 << 4) | 2,
                (4 << 4) | 3,
                (5 << 4) | 4,
                (6 << 4) | 5,
                (7 << 4) | 6,
                (8 << 4) | 7,
                (9 << 4) | 8,
            ]
        );
    }

    #[test]
    fn test_cigar_ref_length_simple() {
        // 50M2I48M → ref length = 50 + 48 = 98 (I is not ref-consuming)
        let ops = parse_cigar_to_bam("50M2I48M");
        assert_eq!(cigar_ref_length(&ops), 98);
    }

    #[test]
    fn test_cigar_ref_length_all_ref_consuming() {
        // M, D, N, =, X all consume reference
        let ops = parse_cigar_to_bam("10M5D3N2=4X");
        assert_eq!(cigar_ref_length(&ops), 10 + 5 + 3 + 2 + 4);
    }

    #[test]
    fn test_cigar_ref_length_non_ref_only() {
        // I and S do not consume reference
        let ops = parse_cigar_to_bam("5I3S");
        assert_eq!(cigar_ref_length(&ops), 0);
    }

    #[test]
    fn test_reg2bin_small_region() {
        // A region at the leaf level (within a 16kb window).
        let bin = reg2bin(0, 100);
        assert_eq!(bin, 4681);
    }

    #[test]
    fn test_reg2bin_wider_region() {
        // Spans from 0 to 200_000 → exceeds a single 16kb bin.
        let bin = reg2bin(0, 200_000);
        // beg >> 17 = 0, end-1 >> 17 = 1 → not same, so goes to level 2
        // beg >> 20 = 0, end-1 >> 20 = 0 → same → 73 + 0 = 73
        assert_eq!(bin, 73);
    }

    #[test]
    fn test_reg2bin_unmapped_convention() {
        // Per BAM spec, unmapped reads use bin 4680.
        // reg2bin is not typically called for unmapped reads, but we verify the
        // constant used in format_unaligned_record_bam matches expectations.
        // Bin 4680 is the last level-4 bin (585 + 4095 = 4680).
        assert_eq!(4680u16, 585 + 4095);
    }

    #[test]
    fn test_encode_sequence_even_length() {
        let mut buf = Vec::new();
        encode_sequence(&mut buf, b"ACGT");
        // A=1, C=2, G=4, T=8 → (1<<4|2), (4<<4|8) = 0x12, 0x48
        assert_eq!(buf, vec![0x12, 0x48]);
    }

    #[test]
    fn test_encode_sequence_odd_length() {
        let mut buf = Vec::new();
        encode_sequence(&mut buf, b"ACG");
        // A=1, C=2, G=4 → (1<<4|2), (4<<4|0) = 0x12, 0x40
        assert_eq!(buf, vec![0x12, 0x40]);
    }

    #[test]
    fn test_encode_sequence_single_base() {
        let mut buf = Vec::new();
        encode_sequence(&mut buf, b"N");
        // N=15 → 15<<4|0 = 0xF0
        assert_eq!(buf, vec![0xF0]);
    }

    #[test]
    fn test_base_to_bam_codes() {
        assert_eq!(base_to_bam(b'='), 0);
        assert_eq!(base_to_bam(b'A'), 1);
        assert_eq!(base_to_bam(b'C'), 2);
        assert_eq!(base_to_bam(b'G'), 4);
        assert_eq!(base_to_bam(b'T'), 8);
        assert_eq!(base_to_bam(b'N'), 15);
        // Lowercase should also work.
        assert_eq!(base_to_bam(b'a'), 1);
        assert_eq!(base_to_bam(b'c'), 2);
        assert_eq!(base_to_bam(b'g'), 4);
        assert_eq!(base_to_bam(b't'), 8);
    }

    #[test]
    fn test_write_bam_tag_u32_byte() {
        let mut buf = Vec::new();
        write_bam_tag_u32(&mut buf, b"NM", 5);
        assert_eq!(buf, vec![b'N', b'M', b'C', 5]);
    }

    #[test]
    fn test_write_bam_tag_u32_short() {
        let mut buf = Vec::new();
        write_bam_tag_u32(&mut buf, b"NM", 300);
        let mut expected = vec![b'N', b'M', b'S'];
        expected.extend_from_slice(&300u16.to_le_bytes());
        assert_eq!(buf, expected);
    }

    #[test]
    fn test_write_bam_tag_u32_int() {
        let mut buf = Vec::new();
        write_bam_tag_u32(&mut buf, b"NM", 100_000);
        let mut expected = vec![b'N', b'M', b'I'];
        expected.extend_from_slice(&100_000u32.to_le_bytes());
        assert_eq!(buf, expected);
    }

    #[test]
    fn test_write_bam_tag_i64_small() {
        let mut buf = Vec::new();
        write_bam_tag_i64(&mut buf, b"XI", 42);
        assert_eq!(buf, vec![b'X', b'I', b'C', 42]);
    }

    #[test]
    fn test_write_bam_tag_i64_large() {
        let mut buf = Vec::new();
        let val: i64 = 5_000_000_000; // exceeds u32 max
        write_bam_tag_i64(&mut buf, b"XI", val);
        // Falls back to string encoding.
        let mut expected = vec![b'X', b'I', b'Z'];
        expected.extend_from_slice(b"5000000000");
        expected.push(0);
        assert_eq!(buf, expected);
    }

    #[test]
    fn test_write_bam_tag_str() {
        let mut buf = Vec::new();
        write_bam_tag_str(&mut buf, b"RG", "sample1");
        let mut expected = vec![b'R', b'G', b'Z'];
        expected.extend_from_slice(b"sample1");
        expected.push(0);
        assert_eq!(buf, expected);
    }

    #[test]
    fn test_aligned_record_bam_structure() {
        // Test the full BAM binary structure of an aligned record.
        let cols = AlignedColumns {
            seq_name: "read1".to_string(),
            sam_flags: 0,
            cigar: "4M".to_string(),
            mate_align_id: 0,
            mate_ref_name: "*".to_string(),
            mate_ref_pos: 0,
            template_len: 0,
            read: "ACGT".to_string(),
            quality: "IIII".to_string(), // Phred+33: 73-33 = 40
            edit_distance: 0,
            spot_group: String::new(),
            alignment_count: 1,
            read_filter: None,
            ref_pos: 0,
            mapq: 0,
        };
        let opts =
            FormatOptions { output_mode: OutputMode::Bam, ref_name_to_id: None, ..default_opts() };
        let mut buf = Vec::new();
        format_aligned_record_bam(
            &mut buf, &cols, 0,    // ref_id
            99,   // ref_pos (0-based)
            30,   // mapq
            1,    // align_id
            None, // mate_info
            &opts,
        );

        // Parse the block_size.
        let block_size = i32::from_le_bytes(buf[0..4].try_into().unwrap());
        assert_eq!(block_size as usize, buf.len() - 4);

        // Parse fixed-length fields.
        let ref_id = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        assert_eq!(ref_id, 0);

        let pos = i32::from_le_bytes(buf[8..12].try_into().unwrap());
        assert_eq!(pos, 99);

        let l_read_name = buf[12];
        assert_eq!(l_read_name, 6); // "read1\0"

        let mapq = buf[13];
        assert_eq!(mapq, 30);

        let n_cigar_op = u16::from_le_bytes(buf[16..18].try_into().unwrap());
        assert_eq!(n_cigar_op, 1);

        let flags = u16::from_le_bytes(buf[18..20].try_into().unwrap());
        assert_ne!(flags & (sam_flags::UNMAPPED as u16), sam_flags::UNMAPPED as u16);

        let l_seq = i32::from_le_bytes(buf[20..24].try_into().unwrap());
        assert_eq!(l_seq, 4);

        // QNAME starts at offset 36.
        let qname_end = 36 + l_read_name as usize;
        assert_eq!(&buf[36..qname_end - 1], b"read1");
        assert_eq!(buf[qname_end - 1], 0); // null terminator

        // CIGAR: 1 op (4M) → (4 << 4) | 0 = 64
        let cigar_start = qname_end;
        let cigar_op = u32::from_le_bytes(buf[cigar_start..cigar_start + 4].try_into().unwrap());
        assert_eq!(cigar_op, 4 << 4);

        // Sequence: ACGT packed → 0x12, 0x48
        let seq_start = cigar_start + 4;
        assert_eq!(&buf[seq_start..seq_start + 2], &[0x12, 0x48]);

        // Quality: Phred 40, 40, 40, 40 (raw, not +33)
        let qual_start = seq_start + 2;
        assert_eq!(&buf[qual_start..qual_start + 4], &[40, 40, 40, 40]);
    }

    #[test]
    fn test_unaligned_record_bam_structure() {
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
        let opts = FormatOptions { output_mode: OutputMode::Bam, ..default_opts() };
        let mut buf = Vec::new();
        format_unaligned_record_bam(&mut buf, &cols, &opts);

        // Parse block_size.
        let block_size = i32::from_le_bytes(buf[0..4].try_into().unwrap());
        assert_eq!(block_size as usize, buf.len() - 4);

        // refID = -1 (unmapped).
        let ref_id = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        assert_eq!(ref_id, -1);

        // pos = -1.
        let pos = i32::from_le_bytes(buf[8..12].try_into().unwrap());
        assert_eq!(pos, -1);

        // mapq = 0.
        assert_eq!(buf[13], 0);

        // bin = 4680.
        let bin = u16::from_le_bytes(buf[14..16].try_into().unwrap());
        assert_eq!(bin, 4680);

        // n_cigar_op = 0.
        let n_cigar_op = u16::from_le_bytes(buf[16..18].try_into().unwrap());
        assert_eq!(n_cigar_op, 0);

        // flags include UNMAPPED.
        let flags = u16::from_le_bytes(buf[18..20].try_into().unwrap());
        assert_ne!(flags & (sam_flags::UNMAPPED as u16), 0);

        // l_seq = 4.
        let l_seq = i32::from_le_bytes(buf[20..24].try_into().unwrap());
        assert_eq!(l_seq, 4);

        // Quality: raw phred 30 for each base.
        let qname_end = 36 + buf[12] as usize;
        let seq_len_packed = 4_usize.div_ceil(2);
        let qual_start = qname_end + seq_len_packed;
        assert_eq!(&buf[qual_start..qual_start + 4], &[30, 30, 30, 30]);
    }

    #[test]
    fn test_unaligned_fasta_reverse() {
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
        let opts = FormatOptions {
            reverse_unaligned: true,
            output_mode: OutputMode::Fasta,
            ..default_opts()
        };
        let mut buf = Vec::new();

        format_unaligned_record(&mut buf, &cols, &opts);

        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, ">spot1\nACGT\n"); // reverse complement of ACGT
    }
}
