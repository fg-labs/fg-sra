//! Aligned read processing via PlacementSetIterator.
//!
//! Walks references → windows → positions → records, building SAM/BAM
//! output for each aligned read. Supports parallel processing across references.

use anyhow::{Context, Result};
use fg_sra_vdb::cursor::VCursor;
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::iterator::{AlignIdSrc, AlignMgr, PlacementIterator, PlacementSetIterator};
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

use crate::matecache::{MateCache, MateInfo};
use crate::output::OutputWriter;
use crate::record::{AlignedColumns, FormatOptions, format_aligned_record};

/// VDB column names for aligned records (with type casts).
mod col {
    pub const SAM_FLAGS: &str = "(U32)SAM_FLAGS";
    pub const CIGAR_SHORT: &str = "(ascii)CIGAR_SHORT";
    pub const CIGAR_LONG: &str = "(ascii)CIGAR_LONG";
    pub const MATE_ALIGN_ID: &str = "(I64)MATE_ALIGN_ID";
    pub const MATE_REF_NAME: &str = "(ascii)MATE_REF_NAME";
    pub const MATE_REF_POS: &str = "(INSDC:coord:zero)MATE_REF_POS";
    pub const TEMPLATE_LEN: &str = "(I32)TEMPLATE_LEN";
    pub const READ: &str = "(ascii)READ";
    pub const SAM_QUALITY: &str = "(INSDC:quality:text:phred_33)SAM_QUALITY";
    pub const EDIT_DISTANCE: &str = "(U32)EDIT_DISTANCE";
    pub const SEQ_SPOT_GROUP: &str = "(ascii)SEQ_SPOT_GROUP";
    pub const SEQ_NAME: &str = "(ascii)SEQ_NAME";
    pub const ALIGNMENT_COUNT: &str = "(U8)ALIGNMENT_COUNT";
    pub const READ_FILTER: &str = "(INSDC:SRA:read_filter)READ_FILTER";
}

/// Column indices for aligned read cursor.
struct AlignColumnIndices {
    sam_flags: u32,
    cigar: u32,
    mate_align_id: u32,
    mate_ref_name: u32,
    mate_ref_pos: u32,
    template_len: u32,
    read: u32,
    sam_quality: u32,
    edit_distance: u32,
    seq_spot_group: u32,
    seq_name: u32,
    alignment_count: Option<u32>,
    read_filter: Option<u32>,
}

/// Set up a VDB cursor with all alignment columns.
///
/// The cursor is NOT opened here — the `PlacementIterator` creation code
/// adds its own columns (`REF_POS`, `REF_LEN`, `MAPQ`, `SPOT_GROUP`) via
/// `TableReader_MakeCursor` and then opens the cursor.  Since columns cannot
/// be added to an already-open cursor, we must defer the open.
fn setup_align_cursor(
    db: &VDatabase,
    table_name: &str,
    use_long_cigar: bool,
) -> Result<(VCursor, AlignColumnIndices)> {
    let table = db.open_table_read(table_name).context("failed to open alignment table")?;
    let cursor = table.create_cursor_read().context("failed to create alignment cursor")?;

    let cigar_col = if use_long_cigar { col::CIGAR_LONG } else { col::CIGAR_SHORT };

    let indices = AlignColumnIndices {
        sam_flags: cursor.add_column(col::SAM_FLAGS).context("SAM_FLAGS")?,
        cigar: cursor.add_column(cigar_col).context("CIGAR")?,
        mate_align_id: cursor.add_column(col::MATE_ALIGN_ID).context("MATE_ALIGN_ID")?,
        mate_ref_name: cursor.add_column(col::MATE_REF_NAME).context("MATE_REF_NAME")?,
        mate_ref_pos: cursor.add_column(col::MATE_REF_POS).context("MATE_REF_POS")?,
        template_len: cursor.add_column(col::TEMPLATE_LEN).context("TEMPLATE_LEN")?,
        read: cursor.add_column(col::READ).context("READ")?,
        sam_quality: cursor.add_column(col::SAM_QUALITY).context("SAM_QUALITY")?,
        edit_distance: cursor.add_column(col::EDIT_DISTANCE).context("EDIT_DISTANCE")?,
        seq_spot_group: cursor.add_column(col::SEQ_SPOT_GROUP).context("SEQ_SPOT_GROUP")?,
        seq_name: cursor.add_column(col::SEQ_NAME).context("SEQ_NAME")?,
        alignment_count: cursor.add_column_optional(col::ALIGNMENT_COUNT),
        read_filter: cursor.add_column_optional(col::READ_FILTER),
    };

    // Do NOT call cursor.open() here — PlacementIterator::make will open it.
    Ok((cursor, indices))
}

/// Read all column data for one aligned record from the cursor.
fn read_aligned_columns(
    cursor: &VCursor,
    idx: &AlignColumnIndices,
    row_id: i64,
) -> Result<AlignedColumns> {
    Ok(AlignedColumns {
        seq_name: cursor.read_str(row_id, idx.seq_name)?,
        sam_flags: cursor.read_u32(row_id, idx.sam_flags)?,
        cigar: cursor.read_str(row_id, idx.cigar)?,
        mate_align_id: cursor.read_i64(row_id, idx.mate_align_id)?,
        mate_ref_name: cursor.read_str(row_id, idx.mate_ref_name)?,
        mate_ref_pos: cursor.read_coord_zero(row_id, idx.mate_ref_pos)?,
        template_len: cursor.read_i32(row_id, idx.template_len)?,
        read: cursor.read_str(row_id, idx.read)?,
        quality: cursor.read_str(row_id, idx.sam_quality)?,
        edit_distance: cursor.read_u32(row_id, idx.edit_distance)?,
        spot_group: cursor.read_str(row_id, idx.seq_spot_group)?,
        alignment_count: match idx.alignment_count {
            Some(col) => cursor.read_u8(row_id, col)?,
            None => 0,
        },
        read_filter: match idx.read_filter {
            Some(col) => Some(cursor.read_u8(row_id, col)?),
            None => None,
        },
    })
}

/// Process all aligned reads for one alignment table, writing SAM records.
///
/// Iterates references → windows → positions → records via the VDB
/// PlacementSetIterator, populating and consulting the mate cache for
/// paired-end information.
pub fn process_aligned_table(
    db: &VDatabase,
    writer: &mut OutputWriter,
    use_seqid: bool,
    use_long_cigar: bool,
    primary_only: bool,
    min_mapq: Option<u32>,
    opts: &FormatOptions<'_>,
) -> Result<()> {
    let mut reflist_opts = reflist_options::USE_PRIMARY_IDS;
    if !primary_only {
        reflist_opts |= reflist_options::USE_SECONDARY_IDS;
    }
    let reflist = ReferenceList::make_database(db, reflist_opts, 0)
        .context("failed to create ReferenceList")?;
    let ref_count = reflist.count().context("failed to get reference count")?;
    if ref_count == 0 {
        return Ok(());
    }

    let min_mapq_val = min_mapq.map(|m| m as i32).unwrap_or(0);

    // Process primary alignments.
    process_alignment_table_with_id_src(
        db,
        "PRIMARY_ALIGNMENT",
        &reflist,
        ref_count,
        writer,
        use_seqid,
        use_long_cigar,
        min_mapq_val,
        AlignIdSrc::Primary,
        opts,
    )?;

    // Process secondary alignments if requested and table exists.
    if !primary_only && db.has_table("SECONDARY_ALIGNMENT") {
        process_alignment_table_with_id_src(
            db,
            "SECONDARY_ALIGNMENT",
            &reflist,
            ref_count,
            writer,
            use_seqid,
            use_long_cigar,
            min_mapq_val,
            AlignIdSrc::Secondary,
            opts,
        )?;
    }

    Ok(())
}

/// Process a single alignment table (primary or secondary).
#[allow(clippy::too_many_arguments)]
fn process_alignment_table_with_id_src(
    db: &VDatabase,
    table_name: &str,
    reflist: &ReferenceList,
    ref_count: u32,
    writer: &mut OutputWriter,
    use_seqid: bool,
    use_long_cigar: bool,
    min_mapq_val: i32,
    id_src: AlignIdSrc,
    opts: &FormatOptions<'_>,
) -> Result<()> {
    let (cursor, col_idx) = setup_align_cursor(db, table_name, use_long_cigar)?;

    let align_mgr = AlignMgr::make_read().context("failed to create AlignMgr")?;
    let mut psi =
        align_mgr.make_placement_set_iterator().context("failed to create PlacementSetIterator")?;

    for i in 0..ref_count {
        let ref_obj = reflist.get(i).context("failed to get reference")?;
        let ref_len = ref_obj.seq_length().context("failed to get reference length")?;

        // Skip zero-length references (the VDB API requires ref_window_len >= 1).
        if ref_len == 0 {
            continue;
        }

        // PlacementIterator::make returns rcDone for references with no alignments
        // in this table (e.g. unplaced contigs). Skip those gracefully.
        let pi = match PlacementIterator::make(&ref_obj, 0, ref_len, min_mapq_val, &cursor, id_src)
        {
            Ok(pi) => pi,
            Err(e) if e.is_done() => continue,
            Err(e) => return Err(e).context("failed to create PlacementIterator"),
        };
        // add_placement_iterator returns Ok(false) when the iterator has no
        // placements (rcDone), matching sam-dump's behavior of continuing.
        match psi.add_placement_iterator(pi) {
            Ok(_) => {}
            Err(e) => return Err(e).context("failed to add PlacementIterator"),
        }
    }

    process_placement_set(&mut psi, &cursor, &col_idx, writer, use_seqid, opts)
}

/// Walk a PlacementSetIterator and emit SAM records.
fn process_placement_set(
    psi: &mut PlacementSetIterator,
    cursor: &VCursor,
    col_idx: &AlignColumnIndices,
    writer: &mut OutputWriter,
    use_seqid: bool,
    opts: &FormatOptions<'_>,
) -> Result<()> {
    let mut buf = Vec::with_capacity(1024);
    let mut mate_cache = MateCache::new();

    while let Some(next_ref) = psi.next_reference()? {
        let ref_name = if use_seqid { next_ref.seq_id()? } else { next_ref.ref_name()? };

        mate_cache.clear();

        while psi.next_window()?.is_some() {
            while let Some((_pos, _len)) = psi.next_avail_pos()? {
                while let Some(rec) = psi.next_record_at(_pos)? {
                    let align_id = rec.id();
                    let ref_pos = rec.pos();
                    let mapq = rec.mapq();

                    let mut cols = read_aligned_columns(cursor, col_idx, align_id)?;

                    // Resolve mate: look up cached mate info and store ours.
                    // When mate has no alignment, strip paired-end flags to match
                    // sam-dump's behavior (the read is output as unpaired).
                    let mate_info = if cols.mate_align_id != 0 {
                        let info = mate_cache.take(cols.mate_align_id);
                        mate_cache.insert(
                            align_id,
                            MateInfo {
                                ref_name: ref_name.clone(),
                                ref_pos,
                                flags: cols.sam_flags,
                                tlen: cols.template_len,
                            },
                        );
                        info
                    } else {
                        cols.strip_paired_flags();
                        None
                    };

                    format_aligned_record(
                        &mut buf,
                        &cols,
                        &ref_name,
                        ref_pos,
                        mapq,
                        align_id,
                        mate_info.as_ref(),
                        opts,
                    );

                    writer.write_record(&buf)?;
                }
            }
        }
    }

    Ok(())
}
