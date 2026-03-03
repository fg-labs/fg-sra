//! Aligned read processing via row-range iteration.
//!
//! Splits alignment row ID ranges across worker threads for parallel
//! processing. Each worker directly iterates rows from a VDB cursor,
//! making parallelism independent of the number of references.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use fg_sra_vdb::cursor::VCursor;
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

use crate::matecache::{MateCache, MateInfo};
use crate::output::OutputWriter;
use crate::progress::ProgressLogger;
use crate::record::{AlignedColumns, FormatOptions, format_aligned_record};

/// VDB table name for primary alignments.
const PRIMARY_ALIGNMENT_TABLE: &str = "PRIMARY_ALIGNMENT";
/// VDB table name for secondary alignments.
const SECONDARY_ALIGNMENT_TABLE: &str = "SECONDARY_ALIGNMENT";

/// Minimum number of alignment rows per work item.
const MIN_CHUNK_SIZE: i64 = 10_000;

/// Configuration for aligned read processing, bundling CLI-derived options
/// that are threaded through multiple functions.
pub struct AlignConfig<'a> {
    pub use_seqid: bool,
    pub use_long_cigar: bool,
    pub primary_only: bool,
    pub min_mapq: Option<u32>,
    pub num_threads: usize,
    pub opts: &'a FormatOptions<'a>,
    /// Genomic regions to restrict output to. Empty means all references.
    pub regions: &'a [String],
}

impl AlignConfig<'_> {
    /// Compute the `ReferenceList` option flags from `primary_only`.
    fn reflist_opts(&self) -> u32 {
        let mut opts = reflist_options::USE_PRIMARY_IDS;
        if !self.primary_only {
            opts |= reflist_options::USE_SECONDARY_IDS;
        }
        opts
    }

    /// Convert `min_mapq` to `i32` for post-read filtering.
    fn min_mapq_i32(&self) -> i32 {
        self.min_mapq.map(|m| m as i32).unwrap_or(0)
    }

    /// Compute the resource pool size: `min(num_threads, num_work_items)`, clamped to at least 1.
    fn pool_size(&self, num_work_items: usize) -> usize {
        self.num_threads.min(num_work_items).max(1)
    }
}

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
    pub const REF_POS: &str = "(INSDC:coord:zero)REF_POS";
    pub const MAPQ: &str = "(I32)MAPQ";
    pub const REF_ID: &str = "(I64)REF_ID";
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
    ref_pos: u32,
    mapq: u32,
}

/// Cache capacity for the data cursor (32 MB). Bounds VDB's MRU blob cache
/// per cursor, preventing unbounded memory growth across references.
const DATA_CURSOR_CACHE_BYTES: usize = 32 * 1024 * 1024;

/// Set up a VDB cursor with all data columns for aligned records.
///
/// Creates a cached cursor (bounded at `DATA_CURSOR_CACHE_BYTES`) and opens
/// it immediately.
fn setup_data_cursor(
    db: &VDatabase,
    table_name: &str,
    use_long_cigar: bool,
) -> Result<(VCursor, AlignColumnIndices)> {
    let table = db.open_table_read(table_name).context("failed to open alignment table")?;
    let cursor = table
        .create_cached_cursor_read(DATA_CURSOR_CACHE_BYTES)
        .context("failed to create data cursor")?;

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
        ref_pos: cursor.add_column(col::REF_POS).context("REF_POS")?,
        mapq: cursor.add_column(col::MAPQ).context("MAPQ")?,
    };

    cursor.open().context("failed to open data cursor")?;

    Ok((cursor, indices))
}

/// Read all column data for one aligned record into a reusable `AlignedColumns`.
///
/// Clears and repopulates `cols` in-place, reusing existing String allocations.
fn read_aligned_columns(
    cursor: &VCursor,
    idx: &AlignColumnIndices,
    row_id: i64,
    cols: &mut AlignedColumns,
) -> Result<()> {
    cols.clear();
    cursor.read_str_into(row_id, idx.seq_name, &mut cols.seq_name)?;
    cols.sam_flags = cursor.read_u32(row_id, idx.sam_flags)?;
    cursor.read_str_into(row_id, idx.cigar, &mut cols.cigar)?;
    cols.mate_align_id = cursor.read_i64(row_id, idx.mate_align_id)?;
    cursor.read_str_into(row_id, idx.mate_ref_name, &mut cols.mate_ref_name)?;
    cols.mate_ref_pos = cursor.read_coord_zero(row_id, idx.mate_ref_pos)?;
    cols.template_len = cursor.read_i32(row_id, idx.template_len)?;
    cursor.read_str_into(row_id, idx.read, &mut cols.read)?;
    cursor.read_str_into(row_id, idx.sam_quality, &mut cols.quality)?;
    cols.edit_distance = cursor.read_u32(row_id, idx.edit_distance)?;
    cursor.read_str_into(row_id, idx.seq_spot_group, &mut cols.spot_group)?;
    cols.alignment_count = match idx.alignment_count {
        Some(col) => cursor.read_u8(row_id, col)?,
        None => 0,
    };
    cols.read_filter = match idx.read_filter {
        Some(col) => Some(cursor.read_u8(row_id, col)?),
        None => None,
    };
    cols.ref_pos = cursor.read_coord_zero(row_id, idx.ref_pos)?;
    cols.mapq = cursor.read_i32(row_id, idx.mapq)?;
    Ok(())
}

// ── Reference boundary discovery ─────────────────────────────────────────

/// Per-reference alignment boundary: ref index, name, first/last alignment row IDs.
struct RefBoundary {
    ref_idx: u32,
    ref_name: String,
    first_row: i64,
    last_row: i64,
}

/// Discover per-reference alignment row ranges by probing the `REF_ID` column.
///
/// Since PRIMARY_ALIGNMENT rows are contiguous per reference and ordered by
/// reference index, we can find boundaries with a linear scan over reference
/// transitions using binary search to find each transition point.
///
/// Returns boundaries in reference order, only for references that have alignments.
fn find_ref_boundaries(
    db: &VDatabase,
    table_name: &str,
    reflist: &ReferenceList,
    use_seqid: bool,
) -> Result<Vec<RefBoundary>> {
    let table = db.open_table_read(table_name).context("failed to open alignment table")?;
    let cursor = table.create_cursor_read().context("failed to create boundary cursor")?;
    let ref_id_col = cursor.add_column(col::REF_ID).context("REF_ID")?;
    cursor.open().context("failed to open boundary cursor")?;

    let (first_row, total_count) = cursor.id_range(0).context("failed to get id range")?;
    if total_count == 0 {
        return Ok(Vec::new());
    }
    let last_row = first_row + total_count as i64 - 1;

    let mut boundaries = Vec::new();
    let mut current_start = first_row;

    while current_start <= last_row {
        let ref_id = cursor.read_i64(current_start, ref_id_col)?;

        // Binary search for the last row with this ref_id.
        let mut lo = current_start;
        let mut hi = last_row;
        while lo < hi {
            let mid = lo + (hi - lo + 1) / 2;
            let mid_ref_id = cursor.read_i64(mid, ref_id_col)?;
            if mid_ref_id == ref_id {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        let boundary_end = lo;

        let ref_idx = ref_id as u32;
        let ref_obj = reflist.get(ref_idx).context("failed to get reference for boundary")?;
        let ref_name = if use_seqid { ref_obj.seq_id()? } else { ref_obj.name()? };

        boundaries.push(RefBoundary {
            ref_idx,
            ref_name: ref_name.to_string(),
            first_row: current_start,
            last_row: boundary_end,
        });

        current_start = boundary_end + 1;
    }

    Ok(boundaries)
}

// ── Work item construction ───────────────────────────────────────────────

/// A unit of work: a contiguous range of alignment row IDs from one reference.
#[derive(Clone)]
struct RowRangeWorkItem {
    /// Contiguous index into the work list, used for ordered output.
    order_idx: usize,
    /// Reference index (for BAM ref_id encoding).
    ref_idx: u32,
    /// Reference name (pre-resolved, avoids per-worker ReferenceList).
    ref_name: String,
    /// First alignment row ID (inclusive).
    start_row: i64,
    /// Last alignment row ID (inclusive).
    end_row: i64,
    /// Optional coordinate filter for --aligned-region (0-based start, exclusive end).
    region_filter: Option<(i32, i32)>,
}

/// Compute target chunk size based on total rows and thread count.
fn target_chunk_size(total_rows: i64, num_threads: usize) -> i64 {
    let target = total_rows / (num_threads as i64 * 8);
    target.max(MIN_CHUNK_SIZE)
}

/// Chunk a sequence of `(boundary, region_filter)` pairs into `RowRangeWorkItem`s.
///
/// Shared core for both whole-table and region-filtered work item construction.
fn chunk_boundaries(
    entries: &[(&RefBoundary, Option<(i32, i32)>)],
    num_threads: usize,
) -> Vec<RowRangeWorkItem> {
    let total_rows: i64 = entries.iter().map(|(b, _)| b.last_row - b.first_row + 1).sum();
    if total_rows == 0 {
        return Vec::new();
    }
    let chunk_size = target_chunk_size(total_rows, num_threads);

    let mut work_items = Vec::new();
    for &(boundary, filter) in entries {
        let mut start = boundary.first_row;
        while start <= boundary.last_row {
            let end = (start + chunk_size - 1).min(boundary.last_row);
            work_items.push(RowRangeWorkItem {
                order_idx: work_items.len(),
                ref_idx: boundary.ref_idx,
                ref_name: boundary.ref_name.clone(),
                start_row: start,
                end_row: end,
                region_filter: filter,
            });
            start = end + 1;
        }
    }
    work_items
}

/// Split reference boundaries into row-range work items for parallel processing.
fn collect_row_range_work_items(
    boundaries: &[RefBoundary],
    num_threads: usize,
) -> Vec<RowRangeWorkItem> {
    let entries: Vec<_> = boundaries.iter().map(|b| (b, None)).collect();
    chunk_boundaries(&entries, num_threads)
}

/// A parsed genomic region: reference name with optional coordinate window.
struct Region {
    name: String,
    /// 0-based start (parsed from 1-based input).
    start: Option<u32>,
    /// 0-based exclusive end.
    end: Option<u32>,
}

/// Parse a region string like `"chr1:1000-2000"` or `"chr2"`.
///
/// Coordinates in the input are 1-based; the returned start is converted to
/// 0-based. End remains as-is (1-based end == 0-based exclusive end).
fn parse_region(s: &str) -> Result<Region> {
    if let Some((name, range)) = s.split_once(':') {
        let (start_str, end_str) =
            range.split_once('-').context("invalid region format: expected name:from-to")?;
        let start: u32 =
            start_str.parse::<u32>().context("invalid region start")?.saturating_sub(1);
        let end: u32 = end_str.parse().context("invalid region end")?;
        Ok(Region { name: name.to_owned(), start: Some(start), end: Some(end) })
    } else {
        Ok(Region { name: s.to_owned(), start: None, end: None })
    }
}

/// Build row-range work items for specific genomic regions.
///
/// For each region, finds the matching reference in `boundaries` and creates
/// work items covering that reference's full row range, with a coordinate
/// filter applied post-read.
fn collect_row_range_region_work_items(
    boundaries: &[RefBoundary],
    regions: &[String],
    reflist: &ReferenceList,
    num_threads: usize,
) -> Result<Vec<RowRangeWorkItem>> {
    // Build a map from ref_name → boundary index for quick lookup.
    let name_to_boundary: HashMap<&str, usize> =
        boundaries.iter().enumerate().map(|(i, b)| (b.ref_name.as_str(), i)).collect();

    let mut entries: Vec<(&RefBoundary, Option<(i32, i32)>)> = Vec::new();

    for spec in regions {
        let region = parse_region(spec)?;

        // First check the boundary map by name, then fall back to reflist.find()
        // to handle name/seqid aliasing.
        let boundary_idx = if let Some(&idx) = name_to_boundary.get(region.name.as_str()) {
            idx
        } else {
            // The region name might be the alternate form (name vs seqid). Look up
            // via reflist to get the ref_idx, then find the matching boundary.
            let ref_obj = reflist
                .find(&region.name)
                .with_context(|| format!("reference not found: {}", region.name))?;
            let ref_idx = ref_obj.idx().context("failed to get reference index")?;
            boundaries
                .iter()
                .position(|b| b.ref_idx == ref_idx)
                .with_context(|| format!("no alignments found for reference: {}", region.name))?
        };

        let filter = match (region.start, region.end) {
            (Some(s), Some(e)) => Some((s as i32, e as i32)),
            _ => None,
        };

        entries.push((&boundaries[boundary_idx], filter));
    }

    Ok(chunk_boundaries(&entries, num_threads))
}

// ── Row-range processing core ────────────────────────────────────────────

/// Process a contiguous range of alignment rows, emitting formatted records.
///
/// Iterates rows `start_row..=end_row`, applying MAPQ and region filters,
/// resolving mate information via the mate cache, and formatting each record
/// via the `emit` callback.
#[allow(clippy::too_many_arguments)]
fn process_row_range(
    cursor: &VCursor,
    col_idx: &AlignColumnIndices,
    item: &RowRangeWorkItem,
    min_mapq: i32,
    opts: &FormatOptions<'_>,
    state: &mut WorkerState,
    mut emit: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    state.mate_cache.clear();

    for row_id in item.start_row..=item.end_row {
        read_aligned_columns(cursor, col_idx, row_id, &mut state.cols)?;

        // Post-read MAPQ filter.
        if state.cols.mapq < min_mapq {
            continue;
        }

        // Post-read region coordinate filter.
        if let Some((rs, re)) = item.region_filter {
            if state.cols.ref_pos < rs || state.cols.ref_pos >= re {
                continue;
            }
        }

        // Resolve mate: look up cached mate info and store ours.
        // When mate has no alignment, strip paired-end flags to match
        // sam-dump's behavior (the read is output as unpaired).
        let mate_info = if state.cols.mate_align_id != 0 {
            let info = state.mate_cache.take(state.cols.mate_align_id);
            state.mate_cache.insert(
                row_id,
                MateInfo { ref_pos: state.cols.ref_pos, tlen: state.cols.template_len },
            );
            info
        } else {
            state.cols.strip_paired_flags();
            None
        };

        format_aligned_record(
            &mut state.record_buf,
            &state.cols,
            &item.ref_name,
            item.ref_idx as i32,
            state.cols.ref_pos,
            state.cols.mapq,
            row_id,
            mate_info.as_ref(),
            opts,
        );

        emit(&state.record_buf)?;
    }

    Ok(())
}

// ── Table-level dispatch ─────────────────────────────────────────────────

/// Process all aligned reads, writing SAM/BAM records.
///
/// Discovers per-reference alignment boundaries, splits into row-range work
/// items, and dispatches to sequential or parallel processing.
pub fn process_aligned_table(
    db: &VDatabase,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    progress_interval: u64,
) -> Result<()> {
    let mut process_table = |table_name: &str| -> Result<()> {
        // Build work items in a block scope so the ReferenceList is dropped
        // before we create per-worker cursors.
        let work_items = {
            let reflist = ReferenceList::make_database(db, config.reflist_opts(), 0)
                .context("failed to create ReferenceList")?;
            let boundaries = find_ref_boundaries(db, table_name, &reflist, config.use_seqid)?;
            if config.regions.is_empty() {
                collect_row_range_work_items(&boundaries, config.num_threads)
            } else {
                collect_row_range_region_work_items(
                    &boundaries,
                    config.regions,
                    &reflist,
                    config.num_threads,
                )?
            }
            // reflist dropped here
        };

        if work_items.is_empty() {
            return Ok(());
        }

        let progress = ProgressLogger::new(work_items.len() as u32, progress_interval);

        if config.num_threads <= 1 {
            process_table_sequential(db, table_name, &work_items, writer, config, &progress)?;
        } else {
            process_table_parallel(db, table_name, &work_items, writer, config, &progress)?;
        }

        progress.complete();
        Ok(())
    };

    process_table(PRIMARY_ALIGNMENT_TABLE)?;
    if !config.primary_only && db.has_table(SECONDARY_ALIGNMENT_TABLE) {
        process_table(SECONDARY_ALIGNMENT_TABLE)?;
    }

    Ok(())
}

/// Process a single alignment table sequentially using row-range iteration.
fn process_table_sequential(
    db: &VDatabase,
    table_name: &str,
    work_items: &[RowRangeWorkItem],
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    progress: &ProgressLogger,
) -> Result<()> {
    let (cursor, col_idx) = setup_data_cursor(db, table_name, config.use_long_cigar)?;
    let min_mapq = config.min_mapq_i32();

    let mut state = WorkerState {
        record_buf: Vec::with_capacity(1024),
        mate_cache: MateCache::new(),
        cols: AlignedColumns::new(),
    };

    for item in work_items {
        process_row_range(&cursor, &col_idx, item, min_mapq, config.opts, &mut state, |rec| {
            progress.record(1);
            writer.write_bytes(rec)
        })?;
        progress.reference_done();
    }

    Ok(())
}

// ── Parallel processing ──────────────────────────────────────────────────

/// Maximum bytes a worker buffers before sending a chunk to the collector.
const CHUNK_SIZE: usize = 8 * 1024 * 1024; // 8 MB

/// A chunk of formatted SAM output from a worker thread.
struct ResultChunk {
    /// Output ordering index (matches `RowRangeWorkItem::order_idx`).
    order_idx: usize,
    /// Chunk sequence number within this work item (0, 1, 2, ...).
    chunk_seq: usize,
    /// The formatted bytes (complete SAM lines only).
    data: Vec<u8>,
    /// True if this is the last chunk for this work item.
    is_last: bool,
}

/// VDB resource set: cursor + column indices.
/// Checked out from the bounded pool by workers, returned after processing.
struct ResourceSet {
    cursor: VCursor,
    col_idx: AlignColumnIndices,
}

/// Process one alignment table in parallel across worker threads.
///
/// Creates a bounded pool of K VDB resource sets (K = `pool_size`), distributes
/// row-range work items via a work channel, and collects chunked results in order.
fn process_table_parallel(
    db: &VDatabase,
    table_name: &str,
    work_items: &[RowRangeWorkItem],
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    progress: &ProgressLogger,
) -> Result<()> {
    let effective_threads = config.num_threads.min(work_items.len());
    let pool_size = config.pool_size(work_items.len());

    // Create bounded resource pool with pool_size sets.
    let (resource_pool_tx, resource_pool_rx) = bounded::<ResourceSet>(pool_size);
    for _ in 0..pool_size {
        let (cursor, col_idx) = setup_data_cursor(db, table_name, config.use_long_cigar)?;
        resource_pool_tx
            .send(ResourceSet { cursor, col_idx })
            .expect("resource pool channel should not be full during init");
    }

    let (work_tx, work_rx) = bounded::<RowRangeWorkItem>(effective_threads * 2);
    let (result_tx, result_rx) = bounded::<ResultChunk>(effective_threads * 2);
    // Buffer pool: collector returns emptied Vec<u8>s for workers to reuse,
    // capping total 8MB allocations to ~2× the number of workers.
    let (buf_pool_tx, buf_pool_rx) = bounded::<Vec<u8>>(effective_threads * 2);

    std::thread::scope(|s| -> Result<()> {
        // Spawn worker threads — they share the bounded resource pool.
        let mut worker_handles = Vec::with_capacity(effective_threads);
        for _ in 0..effective_threads {
            let channels = PoolWorkerChannels {
                work_rx: work_rx.clone(),
                result_tx: result_tx.clone(),
                buf_pool_rx: buf_pool_rx.clone(),
                resource_pool_rx: resource_pool_rx.clone(),
                resource_pool_tx: resource_pool_tx.clone(),
            };
            worker_handles.push(
                s.spawn(move || -> Result<()> { pool_worker_loop(channels, config, progress) }),
            );
        }
        // Drop our copies so only workers hold channel ends.
        drop(work_rx);
        drop(result_tx);
        drop(buf_pool_rx);
        drop(resource_pool_rx);
        drop(resource_pool_tx);

        // Sender thread — feeds work items to workers.
        let work_items_owned: Vec<RowRangeWorkItem> = work_items.to_vec();
        s.spawn(move || {
            for item in work_items_owned {
                if work_tx.send(item).is_err() {
                    break; // Workers died — stop sending.
                }
            }
            // work_tx dropped here, closing the work channel.
        });

        // Collector — runs on the main thread, writes chunks in reference order.
        collect_ordered_chunks(&result_rx, &buf_pool_tx, writer)?;

        // Workers have finished (result channel closed). Check for errors.
        for handle in worker_handles {
            handle.join().expect("worker thread panicked")?;
        }

        Ok(())
    })
}

/// Collect ResultChunks from workers and write them in reference order.
///
/// Chunks arrive out of order from multiple workers. This function buffers
/// them and writes in strict (order_idx, chunk_seq) order, flushing as soon
/// as the next expected chunk becomes available. Written buffers are returned
/// to workers via `pool_tx` for reuse.
fn collect_ordered_chunks(
    result_rx: &Receiver<ResultChunk>,
    pool_tx: &Sender<Vec<u8>>,
    writer: &mut impl Write,
) -> Result<()> {
    let mut next_order: usize = 0;
    // Per-work-item: next chunk_seq we expect to write.
    let mut next_chunk_seq: BTreeMap<usize, usize> = BTreeMap::new();
    // Per-work-item: chunk_seq of the is_last chunk (once seen).
    let mut last_chunk_seq: BTreeMap<usize, usize> = BTreeMap::new();
    // Buffered chunks waiting to be written, keyed by (order_idx, chunk_seq).
    let mut pending: BTreeMap<(usize, usize), Vec<u8>> = BTreeMap::new();

    // Return a written buffer to the pool for worker reuse (best-effort).
    let recycle = |mut buf: Vec<u8>| {
        buf.clear();
        let _ = pool_tx.try_send(buf);
    };

    for chunk in result_rx {
        let expected_seq = next_chunk_seq.get(&next_order).copied().unwrap_or(0);

        // Fast path: if this chunk is exactly the next one we need, write it
        // directly without inserting into the BTreeMap.
        if chunk.order_idx == next_order && chunk.chunk_seq == expected_seq {
            if !chunk.data.is_empty() {
                writer.write_all(&chunk.data)?;
            }
            recycle(chunk.data);
            if chunk.is_last {
                next_chunk_seq.remove(&next_order);
                next_order += 1;
            } else {
                next_chunk_seq.insert(next_order, expected_seq + 1);
            }
        } else {
            // Out-of-order chunk — buffer it.
            if chunk.is_last {
                last_chunk_seq.insert(chunk.order_idx, chunk.chunk_seq);
            }
            pending.insert((chunk.order_idx, chunk.chunk_seq), chunk.data);
        }

        // Flush any buffered chunks that are now in order.
        loop {
            let expected_seq = next_chunk_seq.get(&next_order).copied().unwrap_or(0);
            if let Some(data) = pending.remove(&(next_order, expected_seq)) {
                if !data.is_empty() {
                    writer.write_all(&data)?;
                }
                recycle(data);
                if last_chunk_seq.get(&next_order) == Some(&expected_seq) {
                    // This work item is complete — advance to the next one.
                    next_chunk_seq.remove(&next_order);
                    last_chunk_seq.remove(&next_order);
                    next_order += 1;
                } else {
                    next_chunk_seq.insert(next_order, expected_seq + 1);
                }
            } else {
                break;
            }
        }
    }

    Ok(())
}

/// Channels used by a pool-based worker thread.
struct PoolWorkerChannels {
    work_rx: Receiver<RowRangeWorkItem>,
    result_tx: Sender<ResultChunk>,
    buf_pool_rx: Receiver<Vec<u8>>,
    resource_pool_rx: Receiver<ResourceSet>,
    resource_pool_tx: Sender<ResourceSet>,
}

/// Mutable per-worker state reused across work items.
struct WorkerState {
    record_buf: Vec<u8>,
    mate_cache: MateCache,
    cols: AlignedColumns,
}

/// Worker loop: check out VDB resources from the pool per work item, process,
/// then return them. This bounds total VDB memory to pool_size × per-set cost.
fn pool_worker_loop(
    channels: PoolWorkerChannels,
    config: &AlignConfig<'_>,
    progress: &ProgressLogger,
) -> Result<()> {
    let mut state = WorkerState {
        record_buf: Vec::with_capacity(1024),
        mate_cache: MateCache::new(),
        cols: AlignedColumns::new(),
    };
    let min_mapq = config.min_mapq_i32();

    // Take a buffer from the pool or allocate a new one.
    let take_buf = || -> Vec<u8> {
        channels.buf_pool_rx.try_recv().unwrap_or_else(|_| Vec::with_capacity(CHUNK_SIZE))
    };

    while let Ok(item) = channels.work_rx.recv() {
        // Check out a VDB resource set (blocks if none available).
        let resources = channels
            .resource_pool_rx
            .recv()
            .map_err(|_| anyhow::anyhow!("resource pool closed unexpectedly"))?;

        // Process this row range with the checked-out resources.
        let mut output_buf = take_buf();
        let mut chunk_seq = 0usize;
        let order_idx = item.order_idx;

        let result = process_row_range(
            &resources.cursor,
            &resources.col_idx,
            &item,
            min_mapq,
            config.opts,
            &mut state,
            |rec| {
                progress.record(1);
                output_buf.extend_from_slice(rec);
                if output_buf.len() >= CHUNK_SIZE {
                    channels
                        .result_tx
                        .send(ResultChunk {
                            order_idx,
                            chunk_seq,
                            data: std::mem::replace(&mut output_buf, take_buf()),
                            is_last: false,
                        })
                        .map_err(|_| anyhow::anyhow!("result channel closed"))?;
                    chunk_seq += 1;
                }
                Ok(())
            },
        );

        // Send final chunk for this work item (may be empty).
        channels
            .result_tx
            .send(ResultChunk { order_idx, chunk_seq, data: output_buf, is_last: true })
            .map_err(|_| anyhow::anyhow!("result channel closed"))?;

        // Return resources to pool BEFORE propagating errors.
        channels
            .resource_pool_tx
            .send(resources)
            .map_err(|_| anyhow::anyhow!("resource pool return channel closed"))?;

        progress.reference_done();
        result?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use crossbeam_channel::bounded;

    use super::*;
    use crate::record::FormatOptions;

    /// Build an `AlignConfig` with the given thread count for testing.
    fn test_config(num_threads: usize) -> AlignConfig<'static> {
        static OPTS: FormatOptions<'static> = FormatOptions {
            prefix: None,
            spot_group_in_name: false,
            xi_tag: false,
            reverse_unaligned: false,
            omit_quality: false,
            qual_quant: None,
            output_mode: crate::record::OutputMode::Sam,
            ref_name_to_id: None,
        };
        AlignConfig {
            use_seqid: false,
            use_long_cigar: false,
            primary_only: true,
            min_mapq: None,
            num_threads,
            opts: &OPTS,
            regions: &[],
        }
    }

    #[test]
    fn test_data_cursor_cache_bytes() {
        assert_eq!(DATA_CURSOR_CACHE_BYTES, 32 * 1024 * 1024);
    }

    #[test]
    fn test_pool_size_equals_threads() {
        // With many work items, pool_size equals num_threads.
        let config = test_config(8);
        assert_eq!(config.pool_size(100), 8);
    }

    #[test]
    fn test_pool_size_clamped_by_threads() {
        // With 1 thread, pool_size is 1 regardless of work items.
        let config = test_config(1);
        assert_eq!(config.pool_size(100), 1);
    }

    #[test]
    fn test_pool_size_clamped_by_work_items() {
        // With 1 work item, pool_size is 1 regardless of threads.
        let config = test_config(8);
        assert_eq!(config.pool_size(1), 1);
    }

    #[test]
    fn test_pool_size_minimum_one() {
        let config = test_config(0);
        assert_eq!(config.pool_size(5), 1);
    }

    /// Helper: send chunks into a channel and collect the output via `collect_ordered_chunks`.
    fn run_collector(chunks: Vec<ResultChunk>) -> Vec<u8> {
        let (tx, rx) = bounded::<ResultChunk>(chunks.len() + 1);
        let (pool_tx, _pool_rx) = bounded::<Vec<u8>>(16);
        for chunk in chunks {
            tx.send(chunk).unwrap();
        }
        drop(tx);

        let mut output = Vec::new();
        collect_ordered_chunks(&rx, &pool_tx, &mut output).unwrap();
        output
    }

    #[test]
    fn test_single_ref_single_chunk() {
        let output = run_collector(vec![ResultChunk {
            order_idx: 0,
            chunk_seq: 0,
            data: b"line1\n".to_vec(),
            is_last: true,
        }]);
        assert_eq!(output, b"line1\n");
    }

    #[test]
    fn test_single_ref_multiple_chunks() {
        let output = run_collector(vec![
            ResultChunk { order_idx: 0, chunk_seq: 0, data: b"aaa\n".to_vec(), is_last: false },
            ResultChunk { order_idx: 0, chunk_seq: 1, data: b"bbb\n".to_vec(), is_last: false },
            ResultChunk { order_idx: 0, chunk_seq: 2, data: b"ccc\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"aaa\nbbb\nccc\n");
    }

    #[test]
    fn test_multiple_refs_one_chunk_each() {
        let output = run_collector(vec![
            ResultChunk { order_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 1, chunk_seq: 0, data: b"ref1\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 2, chunk_seq: 0, data: b"ref2\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref1\nref2\n");
    }

    #[test]
    fn test_out_of_order_refs() {
        // ref 1 arrives before ref 0 — should buffer ref 1 and write ref 0 first.
        let output = run_collector(vec![
            ResultChunk { order_idx: 1, chunk_seq: 0, data: b"ref1\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref1\n");
    }

    #[test]
    fn test_multiple_refs_multiple_chunks_out_of_order() {
        // Interleaved chunks from two references, arriving out of order.
        let output = run_collector(vec![
            ResultChunk { order_idx: 1, chunk_seq: 0, data: b"r1c0\n".to_vec(), is_last: false },
            ResultChunk { order_idx: 0, chunk_seq: 1, data: b"r0c1\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 1, chunk_seq: 1, data: b"r1c1\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 0, chunk_seq: 0, data: b"r0c0\n".to_vec(), is_last: false },
        ]);
        assert_eq!(output, b"r0c0\nr0c1\nr1c0\nr1c1\n");
    }

    #[test]
    fn test_empty_ref() {
        // A reference with only an empty final chunk should produce no output.
        let output = run_collector(vec![
            ResultChunk { order_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
            ResultChunk { order_idx: 1, chunk_seq: 0, data: Vec::new(), is_last: true },
            ResultChunk { order_idx: 2, chunk_seq: 0, data: b"ref2\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref2\n");
    }

    #[test]
    fn test_large_chunk_sequence() {
        // 10 chunks per reference × 3 refs.
        let mut chunks = Vec::new();
        for order_idx in 0..3 {
            for seq in 0..10 {
                chunks.push(ResultChunk {
                    order_idx,
                    chunk_seq: seq,
                    data: format!("r{order_idx}c{seq}\n").into_bytes(),
                    is_last: seq == 9,
                });
            }
        }
        let output = run_collector(chunks);
        let expected: String =
            (0..3).flat_map(|r| (0..10).map(move |c| format!("r{r}c{c}\n"))).collect();
        assert_eq!(output, expected.as_bytes());
    }

    #[test]
    fn test_parse_region_name_only() {
        let r = parse_region("chr1").unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, None);
        assert_eq!(r.end, None);
    }

    #[test]
    fn test_parse_region_with_coordinates() {
        let r = parse_region("chr1:1000-2000").unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, Some(999)); // 1-based → 0-based
        assert_eq!(r.end, Some(2000)); // 1-based end == 0-based exclusive
    }

    #[test]
    fn test_parse_region_start_one() {
        let r = parse_region("chr2:1-500").unwrap();
        assert_eq!(r.name, "chr2");
        assert_eq!(r.start, Some(0)); // 1 → 0
        assert_eq!(r.end, Some(500));
    }

    #[test]
    fn test_parse_region_invalid_format() {
        // Missing dash in coordinate range.
        assert!(parse_region("chr1:1000").is_err());
    }

    #[test]
    fn test_parse_region_invalid_numbers() {
        assert!(parse_region("chr1:abc-2000").is_err());
        assert!(parse_region("chr1:1000-xyz").is_err());
    }

    // ── Row-range work item tests ────────────────────────────────────────

    #[test]
    fn test_collect_row_range_single_ref_small() {
        // Single ref with fewer rows than min chunk size → 1 work item.
        let boundaries = vec![RefBoundary {
            ref_idx: 0,
            ref_name: "chr1".to_string(),
            first_row: 1,
            last_row: 5000,
        }];
        let items = collect_row_range_work_items(&boundaries, 4);
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].start_row, 1);
        assert_eq!(items[0].end_row, 5000);
        assert_eq!(items[0].ref_name, "chr1");
        assert!(items[0].region_filter.is_none());
    }

    #[test]
    fn test_collect_row_range_single_ref_large() {
        // Single ref with 1M rows at 4 threads → chunk_size = max(1M/32, 10K) = 31250.
        let boundaries = vec![RefBoundary {
            ref_idx: 0,
            ref_name: "chr1".to_string(),
            first_row: 1,
            last_row: 1_000_000,
        }];
        let items = collect_row_range_work_items(&boundaries, 4);
        // 1M rows / 31250 per chunk = 32 items.
        assert_eq!(items.len(), 32);
        // First item starts at 1.
        assert_eq!(items[0].start_row, 1);
        assert_eq!(items[0].end_row, 31250);
        // Last item ends at 1M.
        assert_eq!(items[31].end_row, 1_000_000);
        // order_idx is sequential.
        for (i, item) in items.iter().enumerate() {
            assert_eq!(item.order_idx, i);
            assert_eq!(item.ref_idx, 0);
        }
    }

    #[test]
    fn test_collect_row_range_multiple_refs() {
        let boundaries = vec![
            RefBoundary {
                ref_idx: 0,
                ref_name: "chr1".to_string(),
                first_row: 1,
                last_row: 50_000,
            },
            RefBoundary {
                ref_idx: 1,
                ref_name: "chr2".to_string(),
                first_row: 50_001,
                last_row: 100_000,
            },
        ];
        let items = collect_row_range_work_items(&boundaries, 2);
        // total = 100K, chunk_size = max(100K/16, 10K) = 10K.
        // chr1: 50K / 10K = 5 items, chr2: 50K / 10K = 5 items.
        assert_eq!(items.len(), 10);
        // First 5 are chr1.
        for item in &items[..5] {
            assert_eq!(item.ref_idx, 0);
            assert_eq!(item.ref_name, "chr1");
        }
        // Last 5 are chr2.
        for item in &items[5..] {
            assert_eq!(item.ref_idx, 1);
            assert_eq!(item.ref_name, "chr2");
        }
    }

    #[test]
    fn test_collect_row_range_empty_boundaries() {
        let items = collect_row_range_work_items(&[], 8);
        assert!(items.is_empty());
    }

    #[test]
    fn test_target_chunk_size_minimum() {
        // Very few rows should still produce at least MIN_CHUNK_SIZE.
        assert_eq!(target_chunk_size(100, 8), MIN_CHUNK_SIZE);
    }

    #[test]
    fn test_target_chunk_size_scales_with_threads() {
        // 10M rows at 8 threads → 10M / 64 = 156250.
        assert_eq!(target_chunk_size(10_000_000, 8), 156_250);
        // Same rows at 1 thread → 10M / 8 = 1250000.
        assert_eq!(target_chunk_size(10_000_000, 1), 1_250_000);
    }
}
