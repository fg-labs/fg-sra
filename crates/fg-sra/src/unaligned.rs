//! Unaligned read processing via SEQUENCE table scan.
//!
//! Scans the SEQUENCE table for spots that are fully or partially unaligned,
//! producing SAM records for unaligned reads.

use anyhow::{Context, Result};
use fg_sra_vdb::cursor::VCursor;
use fg_sra_vdb::database::VDatabase;

use crate::output::OutputWriter;
use crate::record::{
    FormatOptions, READ_TYPE_BIOLOGICAL, UnalignedColumns, format_unaligned_record,
};

/// VDB column names for the SEQUENCE table.
mod col {
    pub const READ: &str = "(INSDC:dna:text)READ";
    pub const QUALITY: &str = "(INSDC:quality:phred)QUALITY";
    pub const SPOT_GROUP: &str = "SPOT_GROUP";
    pub const READ_START: &str = "READ_START";
    pub const READ_LEN: &str = "READ_LEN";
    pub const READ_TYPE: &str = "READ_TYPE";
    pub const READ_FILTER: &str = "READ_FILTER";
    pub const NAME: &str = "NAME";
    pub const PRIMARY_ALIGNMENT_ID: &str = "PRIMARY_ALIGNMENT_ID";
}

/// Column indices for the SEQUENCE table cursor.
struct SeqColumnIndices {
    read: u32,
    quality: u32,
    spot_group: u32,
    read_start: u32,
    read_len: u32,
    read_type: u32,
    read_filter: u32,
    name: u32,
    primary_alignment_id: u32,
}

/// Set up a cursor on the SEQUENCE table.
fn setup_seq_cursor(db: &VDatabase) -> Result<(VCursor, SeqColumnIndices)> {
    let table = db.open_table_read("SEQUENCE").context("failed to open SEQUENCE table")?;
    let cursor = table.create_cursor_read().context("failed to create SEQUENCE cursor")?;

    let indices = SeqColumnIndices {
        read: cursor.add_column(col::READ).context("READ")?,
        quality: cursor.add_column(col::QUALITY).context("QUALITY")?,
        spot_group: cursor.add_column(col::SPOT_GROUP).context("SPOT_GROUP")?,
        read_start: cursor.add_column(col::READ_START).context("READ_START")?,
        read_len: cursor.add_column(col::READ_LEN).context("READ_LEN")?,
        read_type: cursor.add_column(col::READ_TYPE).context("READ_TYPE")?,
        read_filter: cursor.add_column(col::READ_FILTER).context("READ_FILTER")?,
        name: cursor.add_column(col::NAME).context("NAME")?,
        primary_alignment_id: cursor
            .add_column(col::PRIMARY_ALIGNMENT_ID)
            .context("PRIMARY_ALIGNMENT_ID")?,
    };

    cursor.open().context("failed to open SEQUENCE cursor")?;
    Ok((cursor, indices))
}

/// Process unaligned reads from the SEQUENCE table.
///
/// For each spot, iterates over reads and outputs those that are:
/// - Biological (READ_TYPE has biological bit set)
/// - Unaligned (PRIMARY_ALIGNMENT_ID == 0 for that read)
/// - Non-empty (READ_LEN > 0)
pub fn process_unaligned_reads(
    db: &VDatabase,
    writer: &mut OutputWriter,
    opts: &FormatOptions<'_>,
    unaligned_spots_only: bool,
) -> Result<()> {
    let (cursor, idx) = setup_seq_cursor(db)?;

    let (first_row, row_count) =
        cursor.id_range(idx.read).context("failed to get SEQUENCE row range")?;

    let mut buf = Vec::with_capacity(512);

    for row_id in first_row..first_row + row_count as i64 {
        let primary_ids = cursor.read_i64_slice(row_id, idx.primary_alignment_id)?;
        let read_types = cursor.read_u8_slice(row_id, idx.read_type)?;
        let read_starts = cursor.read_i32_slice(row_id, idx.read_start)?;
        let read_lens = cursor.read_u32_slice(row_id, idx.read_len)?;
        let read_filters = cursor.read_u8_slice(row_id, idx.read_filter)?;

        let nreads = primary_ids.len();
        if nreads == 0 {
            continue;
        }

        // If --unaligned-spots-only, skip spots that have any alignments.
        if unaligned_spots_only && primary_ids.iter().any(|&id| id != 0) {
            continue;
        }

        // Pre-compute which reads to emit (unaligned, biological, non-empty).
        // This avoids iterating twice — once to count, once to emit.
        let mut emit_indices: Vec<usize> = Vec::with_capacity(nreads);
        for (i, primary_id) in primary_ids.iter().enumerate().take(nreads) {
            if *primary_id != 0 {
                continue;
            }
            let read_type = read_types.get(i).copied().unwrap_or(0);
            if (read_type & READ_TYPE_BIOLOGICAL) == 0 {
                continue;
            }
            let read_len = read_lens.get(i).copied().unwrap_or(0);
            if read_len == 0 {
                continue;
            }
            emit_indices.push(i);
        }

        if emit_indices.is_empty() {
            continue;
        }

        let num_bio_reads = emit_indices.len() as u32;

        // Read full spot data only when we have reads to emit.
        let full_read = cursor.read_str(row_id, idx.read)?;
        let full_quality = cursor.read_u8_slice(row_id, idx.quality)?;
        let name = cursor.read_str(row_id, idx.name)?;
        let spot_group = cursor.read_str(row_id, idx.spot_group)?;

        for (bio_index, &i) in emit_indices.iter().enumerate() {
            let read_start = read_starts.get(i).copied().unwrap_or(0) as usize;
            let read_len = read_lens.get(i).copied().unwrap_or(0) as usize;
            let read_end = read_start + read_len;

            let read_seq = safe_slice_str(&full_read, read_start, read_end);
            let qual_slice = safe_slice(&full_quality, read_start, read_end);

            let cols = UnalignedColumns {
                name: &name,
                read: read_seq,
                quality: qual_slice,
                spot_group: &spot_group,
                read_type: read_types.get(i).copied().unwrap_or(0),
                read_filter: read_filters.get(i).copied().unwrap_or(0),
                num_bio_reads,
                bio_read_index: bio_index as u32,
            };

            format_unaligned_record(&mut buf, &cols, opts);
            writer.write_record(&buf)?;
        }
    }

    Ok(())
}

/// Safely slice a byte slice, clamping to available bounds.
fn safe_slice(data: &[u8], start: usize, end: usize) -> &[u8] {
    if start >= data.len() {
        &[]
    } else if end <= data.len() {
        &data[start..end]
    } else {
        &data[start..]
    }
}

/// Safely slice a string, clamping to available bounds.
fn safe_slice_str(data: &str, start: usize, end: usize) -> &str {
    if start >= data.len() {
        ""
    } else if end <= data.len() {
        &data[start..end]
    } else {
        &data[start..]
    }
}
