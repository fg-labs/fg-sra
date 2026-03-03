//! Aligned read processing via PlacementSetIterator.
//!
//! Walks references → windows → positions → records, building SAM/BAM
//! output for each aligned read. Supports parallel processing across references.

use std::collections::BTreeMap;
use std::io::Write;

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use fg_sra_vdb::cursor::VCursor;
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::iterator::{AlignIdSrc, AlignMgr, PlacementIterator, PlacementSetIterator};
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

use crate::matecache::{MateCache, MateInfo};
use crate::output::OutputWriter;
use crate::progress::ProgressLogger;
use crate::record::{AlignedColumns, FormatOptions, format_aligned_record};

/// Default genomic window size for intra-reference parallelism (5 Mbp).
/// References larger than this are split into windows in parallel mode.
const DEFAULT_WINDOW_SIZE: u32 = 5_000_000;

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

    /// Convert `min_mapq` to the `i32` expected by the VDB placement API.
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

/// Cache capacity for the data cursor (32 MB). Bounds VDB's MRU blob cache
/// per cursor, preventing unbounded memory growth across references.
const DATA_CURSOR_CACHE_BYTES: usize = 32 * 1024 * 1024;

/// Set up a VDB cursor with all data columns for aligned records.
///
/// Creates a cached cursor (bounded at `DATA_CURSOR_CACHE_BYTES`) and opens
/// it immediately. This cursor is independent of the PlacementIterator's
/// internal cursor — PI receives NULL and creates its own lightweight cursor.
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
    };

    // Open immediately — PI will NOT open this cursor (it gets NULL).
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
    Ok(())
}

/// Scan reference metadata and build work items, skipping zero-length references.
///
/// When `window_size` is `Some(ws)`, references longer than `ws` are split into
/// multiple work items with sequential `order_idx` values. When `None`, each
/// reference becomes a single whole-reference work item.
fn collect_work_items(
    reflist: &ReferenceList,
    ref_count: u32,
    window_size: Option<u32>,
) -> Result<Vec<WorkItem>> {
    let mut work_items = Vec::with_capacity(ref_count as usize);
    for i in 0..ref_count {
        let ref_obj = reflist.get(i).context("failed to get reference")?;
        let ref_len = ref_obj.seq_length().context("failed to get reference length")?;
        if ref_len == 0 {
            continue;
        }
        match window_size {
            Some(ws) if ref_len > ws => {
                let mut start: u32 = 0;
                while start < ref_len {
                    let len = ws.min(ref_len - start);
                    work_items.push(WorkItem {
                        order_idx: work_items.len(),
                        ref_idx: i,
                        window_start: start,
                        window_len: len,
                    });
                    start += ws;
                }
            }
            _ => {
                work_items.push(WorkItem {
                    order_idx: work_items.len(),
                    ref_idx: i,
                    window_start: 0,
                    window_len: ref_len,
                });
            }
        }
    }
    Ok(work_items)
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

/// Build work items for specific genomic regions.
///
/// Looks up each region name in the `ReferenceList`, restricting to the
/// requested coordinate window (or the full reference if unspecified).
fn collect_region_work_items(
    reflist: &ReferenceList,
    regions: &[String],
    window_size: Option<u32>,
) -> Result<Vec<WorkItem>> {
    let mut work_items = Vec::with_capacity(regions.len());
    for spec in regions {
        let region = parse_region(spec)?;
        let ref_obj = reflist
            .find(&region.name)
            .with_context(|| format!("reference not found: {}", region.name))?;
        let ref_idx = ref_obj.idx().context("failed to get reference index")?;
        let ref_len = ref_obj.seq_length().context("failed to get reference length")?;
        if ref_len == 0 {
            continue;
        }

        let start = region.start.unwrap_or(0).min(ref_len);
        let end = region.end.unwrap_or(ref_len).min(ref_len);
        if start >= end {
            continue;
        }
        let total_len = end - start;

        match window_size {
            Some(ws) if total_len > ws => {
                let mut pos = start;
                while pos < end {
                    let len = ws.min(end - pos);
                    work_items.push(WorkItem {
                        order_idx: work_items.len(),
                        ref_idx,
                        window_start: pos,
                        window_len: len,
                    });
                    pos += ws;
                }
            }
            _ => {
                work_items.push(WorkItem {
                    order_idx: work_items.len(),
                    ref_idx,
                    window_start: start,
                    window_len: total_len,
                });
            }
        }
    }
    Ok(work_items)
}

/// Process all aligned reads for one alignment table, writing SAM records.
///
/// When `num_threads > 1`, references are processed in parallel across worker
/// threads using a bounded resource pool.  When `num_threads <= 1`, the
/// sequential path is used (no threading overhead).
pub fn process_aligned_table(
    db: &VDatabase,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    progress_interval: u64,
) -> Result<()> {
    // Build work items in a block scope so the ReferenceList is dropped
    // before we create per-worker resource sets, avoiding N+1 ReferenceLists.
    let window_size = if config.num_threads > 1 { Some(DEFAULT_WINDOW_SIZE) } else { None };
    let work_items = {
        let reflist = ReferenceList::make_database(db, config.reflist_opts(), 0)
            .context("failed to create ReferenceList")?;
        if config.regions.is_empty() {
            let ref_count = reflist.count().context("failed to get reference count")?;
            collect_work_items(&reflist, ref_count, window_size)?
        } else {
            collect_region_work_items(&reflist, config.regions, window_size)?
        }
        // reflist dropped here
    };

    if work_items.is_empty() {
        return Ok(());
    }

    let progress = ProgressLogger::new(work_items.len() as u32, progress_interval);

    // Dispatch to sequential or parallel processing per table.
    let mut process_table = |table_name: &str, id_src: AlignIdSrc| -> Result<()> {
        if config.num_threads <= 1 {
            process_alignment_table_sequential(
                db,
                table_name,
                &work_items,
                writer,
                config,
                id_src,
                &progress,
            )
        } else {
            process_table_parallel(db, table_name, &work_items, writer, config, id_src, &progress)
        }
    };

    process_table(AlignIdSrc::Primary.table_name(), AlignIdSrc::Primary)?;
    if !config.primary_only && db.has_table(AlignIdSrc::Secondary.table_name()) {
        process_table(AlignIdSrc::Secondary.table_name(), AlignIdSrc::Secondary)?;
    }

    progress.complete();
    Ok(())
}

/// Process a single alignment table sequentially, one reference at a time.
///
/// Creates a fresh `PlacementSetIterator` per reference (matching the parallel
/// path's per-reference approach) to avoid loading all reference alignment ID
/// arrays simultaneously.
fn process_alignment_table_sequential(
    db: &VDatabase,
    table_name: &str,
    work_items: &[WorkItem],
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    id_src: AlignIdSrc,
    progress: &ProgressLogger,
) -> Result<()> {
    let (cursor, col_idx) = setup_data_cursor(db, table_name, config.use_long_cigar)?;
    let reflist = ReferenceList::make_database(db, config.reflist_opts(), 0)
        .context("failed to create ReferenceList")?;
    let min_mapq = config.min_mapq_i32();
    let align_mgr = AlignMgr::make_read().context("failed to create AlignMgr")?;

    let ctx = EmitContext {
        cursor: &cursor,
        col_idx: &col_idx,
        use_seqid: config.use_seqid,
        opts: config.opts,
    };
    let mut record_buf = Vec::with_capacity(1024);
    let mut mate_cache = MateCache::new();
    let mut cols = AlignedColumns::new();

    for item in work_items {
        let ref_obj = reflist.get(item.ref_idx).context("failed to get reference")?;

        let pi = match PlacementIterator::make(
            &ref_obj,
            item.window_start as i32,
            item.window_len,
            min_mapq,
            None,
            id_src,
        ) {
            Ok(pi) => pi,
            Err(e) if e.is_done() => {
                progress.reference_done();
                continue;
            }
            Err(e) => return Err(e).context("failed to create PlacementIterator"),
        };

        let mut psi = align_mgr
            .make_placement_set_iterator()
            .context("failed to create PlacementSetIterator")?;
        match psi.add_placement_iterator(pi) {
            Ok(_) => {}
            Err(e) => return Err(e).context("failed to add PlacementIterator"),
        }

        emit_placement_records(
            &mut psi,
            &ctx,
            &mut record_buf,
            &mut mate_cache,
            &mut cols,
            |rec| {
                progress.record(1);
                writer.write_bytes(rec)
            },
        )?;
        progress.reference_done();
        // psi + pi dropped here, freeing per-reference alignment ID arrays
    }

    Ok(())
}

/// Shared read-only context for emitting aligned records.
///
/// Bundles the cursor, column indices, and formatting options that are
/// constant across all records within a processing run.
struct EmitContext<'a> {
    cursor: &'a VCursor,
    col_idx: &'a AlignColumnIndices,
    use_seqid: bool,
    opts: &'a FormatOptions<'a>,
}

/// Walk a PlacementSetIterator and emit formatted SAM records via the `emit` callback.
///
/// Shared core between the sequential path (emits directly to writer) and the
/// parallel path (emits into a `Vec<u8>` buffer).
fn emit_placement_records(
    psi: &mut PlacementSetIterator,
    ctx: &EmitContext<'_>,
    record_buf: &mut Vec<u8>,
    mate_cache: &mut MateCache,
    cols: &mut AlignedColumns,
    mut emit: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    while let Some(next_ref) = psi.next_reference()? {
        let ref_name = if ctx.use_seqid { next_ref.seq_id()? } else { next_ref.ref_name()? };

        mate_cache.clear();

        while psi.next_window()?.is_some() {
            while let Some((pos, _len)) = psi.next_avail_pos()? {
                while let Some(rec) = psi.next_record_at(pos)? {
                    let align_id = rec.id();
                    let ref_pos = rec.pos();
                    let mapq = rec.mapq();

                    read_aligned_columns(ctx.cursor, ctx.col_idx, align_id, cols)?;

                    // Resolve mate: look up cached mate info and store ours.
                    // When mate has no alignment, strip paired-end flags to match
                    // sam-dump's behavior (the read is output as unpaired).
                    let mate_info = if cols.mate_align_id != 0 {
                        let info = mate_cache.take(cols.mate_align_id);
                        mate_cache.insert(align_id, MateInfo { ref_pos, tlen: cols.template_len });
                        info
                    } else {
                        cols.strip_paired_flags();
                        None
                    };

                    format_aligned_record(
                        record_buf,
                        cols,
                        &ref_name,
                        ref_pos,
                        mapq,
                        align_id,
                        mate_info.as_ref(),
                        ctx.opts,
                    );

                    emit(record_buf)?;
                }
            }
        }
    }

    Ok(())
}

// ── Parallel processing ──────────────────────────────────────────────────

/// Maximum bytes a worker buffers before sending a chunk to the collector.
const CHUNK_SIZE: usize = 8 * 1024 * 1024; // 8 MB

/// A unit of work sent to a worker thread: one reference window to process.
#[derive(Clone, Copy)]
struct WorkItem {
    /// Contiguous index into the filtered work list, used for ordered output.
    order_idx: usize,
    /// Index into the ReferenceList (each worker has its own list).
    ref_idx: u32,
    /// 0-based start position within the reference.
    window_start: u32,
    /// Length of this window (may be less than the full reference).
    window_len: u32,
}

/// A chunk of formatted SAM output from a worker thread.
struct ResultChunk {
    /// Output ordering index (matches `WorkItem::order_idx`).
    order_idx: usize,
    /// Chunk sequence number within this reference (0, 1, 2, ...).
    chunk_seq: usize,
    /// The formatted bytes (complete SAM lines only).
    data: Vec<u8>,
    /// True if this is the last chunk for this reference.
    is_last: bool,
}

/// VDB resource set: cursor + column indices + ReferenceList.
/// Checked out from the bounded pool by workers, returned after processing.
struct ResourceSet {
    cursor: VCursor,
    col_idx: AlignColumnIndices,
    reflist: ReferenceList,
}

/// Process one alignment table in parallel across worker threads.
///
/// Creates a bounded pool of K VDB resource sets (K = `pool_size`), distributes
/// references via a work channel, and collects chunked results in reference order.
/// Workers check out a resource set for each reference and return it when done,
/// bounding peak VDB memory to K × per-set cost.
fn process_table_parallel(
    db: &VDatabase,
    table_name: &str,
    work_items: &[WorkItem],
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    id_src: AlignIdSrc,
    progress: &ProgressLogger,
) -> Result<()> {
    let effective_threads = config.num_threads.min(work_items.len());
    let pool_size = config.pool_size(work_items.len());
    let reflist_opts = config.reflist_opts();

    // Create bounded resource pool with pool_size sets.
    let (resource_pool_tx, resource_pool_rx) = bounded::<ResourceSet>(pool_size);
    for _ in 0..pool_size {
        let (cursor, col_idx) = setup_data_cursor(db, table_name, config.use_long_cigar)?;
        let reflist = ReferenceList::make_database(db, reflist_opts, 0)
            .context("failed to create per-pool ReferenceList")?;
        resource_pool_tx
            .send(ResourceSet { cursor, col_idx, reflist })
            .expect("resource pool channel should not be full during init");
    }

    let (work_tx, work_rx) = bounded::<WorkItem>(effective_threads * 2);
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
            worker_handles.push(s.spawn(move || -> Result<()> {
                pool_worker_loop(channels, config, id_src, progress)
            }));
        }
        // Drop our copies so only workers hold channel ends.
        drop(work_rx);
        drop(result_tx);
        drop(buf_pool_rx);
        drop(resource_pool_rx);
        drop(resource_pool_tx);

        // Sender thread — feeds work items to workers.
        let work_items_owned: Vec<WorkItem> = work_items.to_vec();
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
/// them and writes in strict (ref_idx, chunk_seq) order, flushing as soon
/// as the next expected chunk becomes available. Written buffers are returned
/// to workers via `pool_tx` for reuse.
fn collect_ordered_chunks(
    result_rx: &Receiver<ResultChunk>,
    pool_tx: &Sender<Vec<u8>>,
    writer: &mut impl Write,
) -> Result<()> {
    let mut next_order: usize = 0;
    // Per-reference: next chunk_seq we expect to write.
    let mut next_chunk_seq: BTreeMap<usize, usize> = BTreeMap::new();
    // Per-reference: chunk_seq of the is_last chunk (once seen).
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
                    // This reference is complete — advance to the next one.
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
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<ResultChunk>,
    buf_pool_rx: Receiver<Vec<u8>>,
    resource_pool_rx: Receiver<ResourceSet>,
    resource_pool_tx: Sender<ResourceSet>,
}

/// Mutable per-worker state reused across references.
struct WorkerState {
    record_buf: Vec<u8>,
    mate_cache: MateCache,
    cols: AlignedColumns,
    align_mgr: AlignMgr,
}

/// Worker loop: check out VDB resources from the pool per-reference, process,
/// then return them. This bounds total VDB memory to pool_size × per-set cost.
fn pool_worker_loop(
    channels: PoolWorkerChannels,
    config: &AlignConfig<'_>,
    id_src: AlignIdSrc,
    progress: &ProgressLogger,
) -> Result<()> {
    let mut state = WorkerState {
        record_buf: Vec::with_capacity(1024),
        mate_cache: MateCache::new(),
        cols: AlignedColumns::new(),
        align_mgr: AlignMgr::make_read().context("worker: failed to create AlignMgr")?,
    };

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

        // Process this reference with the checked-out resources.
        let result = process_one_reference(
            &resources,
            &item,
            id_src,
            config,
            &mut state,
            &take_buf,
            &channels.result_tx,
            progress,
        );

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

/// Process a single reference using borrowed VDB resources.
///
/// All VDB borrows (EmitContext, PlacementSetIterator, PlacementIterator,
/// ReferenceObj) are scoped to this function call and dropped on return,
/// allowing the ResourceSet to be returned to the pool.
#[allow(clippy::too_many_arguments)]
fn process_one_reference(
    resources: &ResourceSet,
    item: &WorkItem,
    id_src: AlignIdSrc,
    config: &AlignConfig<'_>,
    state: &mut WorkerState,
    take_buf: &dyn Fn() -> Vec<u8>,
    result_tx: &Sender<ResultChunk>,
    progress: &ProgressLogger,
) -> Result<()> {
    let ctx = EmitContext {
        cursor: &resources.cursor,
        col_idx: &resources.col_idx,
        use_seqid: config.use_seqid,
        opts: config.opts,
    };
    let min_mapq = config.min_mapq_i32();

    let ref_obj = resources.reflist.get(item.ref_idx).context("worker: failed to get reference")?;

    let pi = match PlacementIterator::make(
        &ref_obj,
        item.window_start as i32,
        item.window_len,
        min_mapq,
        None,
        id_src,
    ) {
        Ok(pi) => pi,
        Err(e) if e.is_done() => {
            // No alignments on this reference — send empty final chunk.
            result_tx
                .send(ResultChunk {
                    order_idx: item.order_idx,
                    chunk_seq: 0,
                    data: Vec::new(),
                    is_last: true,
                })
                .map_err(|_| anyhow::anyhow!("result channel closed"))?;
            return Ok(());
        }
        Err(e) => {
            return Err(e).context("worker: failed to create PlacementIterator");
        }
    };

    let mut psi = state
        .align_mgr
        .make_placement_set_iterator()
        .context("worker: failed to create PlacementSetIterator")?;
    match psi.add_placement_iterator(pi) {
        Ok(_) => {}
        Err(e) => {
            return Err(e).context("worker: failed to add PlacementIterator");
        }
    }

    // Process records, sending chunks when the buffer exceeds CHUNK_SIZE.
    let mut output_buf = take_buf();
    let mut chunk_seq = 0usize;
    emit_placement_records(
        &mut psi,
        &ctx,
        &mut state.record_buf,
        &mut state.mate_cache,
        &mut state.cols,
        |rec| {
            progress.record(1);
            output_buf.extend_from_slice(rec);
            if output_buf.len() >= CHUNK_SIZE {
                result_tx
                    .send(ResultChunk {
                        order_idx: item.order_idx,
                        chunk_seq,
                        data: std::mem::replace(&mut output_buf, take_buf()),
                        is_last: false,
                    })
                    .map_err(|_| anyhow::anyhow!("result channel closed"))?;
                chunk_seq += 1;
            }
            Ok(())
        },
    )?;

    // Send final chunk for this reference (may be empty).
    result_tx
        .send(ResultChunk { order_idx: item.order_idx, chunk_seq, data: output_buf, is_last: true })
        .map_err(|_| anyhow::anyhow!("result channel closed"))?;

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
}
