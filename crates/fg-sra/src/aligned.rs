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
use fg_sra_vdb::reference::{ReferenceList, ReferenceObj, reflist_options};

use crate::matecache::{MateCache, MateInfo};
use crate::output::OutputWriter;
use crate::record::{AlignedColumns, FormatOptions, format_aligned_record};

/// Configuration for aligned read processing, bundling CLI-derived options
/// that are threaded through multiple functions.
pub struct AlignConfig<'a> {
    pub use_seqid: bool,
    pub use_long_cigar: bool,
    pub primary_only: bool,
    pub min_mapq: Option<u32>,
    pub num_threads: usize,
    pub opts: &'a FormatOptions<'a>,
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
/// When `num_threads > 1`, references are processed in parallel across worker
/// threads, each with its own VDB cursor.  When `num_threads <= 1`, the existing
/// sequential path is used (no threading overhead).
pub fn process_aligned_table(
    db: &VDatabase,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
) -> Result<()> {
    let mut reflist_opts = reflist_options::USE_PRIMARY_IDS;
    if !config.primary_only {
        reflist_opts |= reflist_options::USE_SECONDARY_IDS;
    }
    let reflist = ReferenceList::make_database(db, reflist_opts, 0)
        .context("failed to create ReferenceList")?;
    let ref_count = reflist.count().context("failed to get reference count")?;
    if ref_count == 0 {
        return Ok(());
    }

    // Dispatch to sequential or parallel processing per table.
    let mut process_table = |table_name: &str, id_src: AlignIdSrc| -> Result<()> {
        if config.num_threads <= 1 {
            process_alignment_table_with_id_src(
                db, table_name, &reflist, ref_count, writer, config, id_src,
            )
        } else {
            process_table_parallel(db, table_name, &reflist, ref_count, writer, config, id_src)
        }
    };

    process_table("PRIMARY_ALIGNMENT", AlignIdSrc::Primary)?;
    if !config.primary_only && db.has_table("SECONDARY_ALIGNMENT") {
        process_table("SECONDARY_ALIGNMENT", AlignIdSrc::Secondary)?;
    }

    Ok(())
}

/// Process a single alignment table (primary or secondary).
fn process_alignment_table_with_id_src(
    db: &VDatabase,
    table_name: &str,
    reflist: &ReferenceList,
    ref_count: u32,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    id_src: AlignIdSrc,
) -> Result<()> {
    let (cursor, col_idx) = setup_align_cursor(db, table_name, config.use_long_cigar)?;
    let min_mapq_val = config.min_mapq.map(|m| m as i32).unwrap_or(0);

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

    process_placement_set(&mut psi, &cursor, &col_idx, writer, config.use_seqid, config.opts)
}

/// Walk a PlacementSetIterator and write SAM records to the output.
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
    emit_placement_records(
        psi,
        cursor,
        col_idx,
        &mut buf,
        &mut mate_cache,
        use_seqid,
        opts,
        |rec| writer.write_bytes(rec),
    )
}

/// Walk a PlacementSetIterator and emit formatted SAM records via the `emit` callback.
///
/// Shared core between the sequential path (emits directly to writer) and the
/// parallel path (emits into a `Vec<u8>` buffer).
#[allow(clippy::too_many_arguments)]
fn emit_placement_records(
    psi: &mut PlacementSetIterator,
    cursor: &VCursor,
    col_idx: &AlignColumnIndices,
    record_buf: &mut Vec<u8>,
    mate_cache: &mut MateCache,
    use_seqid: bool,
    opts: &FormatOptions<'_>,
    mut emit: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
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
                        record_buf,
                        &cols,
                        &ref_name,
                        ref_pos,
                        mapq,
                        align_id,
                        mate_info.as_ref(),
                        opts,
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

/// A unit of work sent to a worker thread: one reference to process.
struct WorkItem {
    /// Contiguous index for ordered output collection.
    idx: usize,
    /// The reference object (moved to the worker, dropped after PI creation).
    ref_obj: ReferenceObj,
    /// Reference sequence length.
    ref_len: u32,
}

/// A chunk of formatted SAM output from a worker thread.
struct ResultChunk {
    /// Reference index for ordered output.
    ref_idx: usize,
    /// Chunk sequence number within this reference (0, 1, 2, ...).
    chunk_seq: usize,
    /// The formatted bytes (complete SAM lines only).
    data: Vec<u8>,
    /// True if this is the last chunk for this reference.
    is_last: bool,
}

/// Process one alignment table in parallel across worker threads.
///
/// Pre-creates one VDB cursor per worker, distributes references via a work
/// channel, and collects chunked results in reference order.
fn process_table_parallel(
    db: &VDatabase,
    table_name: &str,
    reflist: &ReferenceList,
    ref_count: u32,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    id_src: AlignIdSrc,
) -> Result<()> {
    // Pre-fetch reference metadata and build work items (skip zero-length refs).
    let mut work_items = Vec::with_capacity(ref_count as usize);
    for i in 0..ref_count {
        let ref_obj = reflist.get(i).context("failed to get reference")?;
        let ref_len = ref_obj.seq_length().context("failed to get reference length")?;
        if ref_len == 0 {
            continue;
        }
        let idx = work_items.len();
        work_items.push(WorkItem { idx, ref_obj, ref_len });
    }

    if work_items.is_empty() {
        return Ok(());
    }

    // Pre-create one cursor per worker (cursors cannot be shared across threads).
    let effective_threads = config.num_threads.min(work_items.len());
    let cursors: Vec<(VCursor, AlignColumnIndices)> = (0..effective_threads)
        .map(|_| setup_align_cursor(db, table_name, config.use_long_cigar))
        .collect::<Result<Vec<_>>>()?;

    let min_mapq_val = config.min_mapq.map(|m| m as i32).unwrap_or(0);
    let (work_tx, work_rx) = bounded::<WorkItem>(effective_threads * 2);
    let (result_tx, result_rx) = bounded::<ResultChunk>(effective_threads * 2);

    std::thread::scope(|s| -> Result<()> {
        // Spawn worker threads — each owns one cursor.
        let mut worker_handles = Vec::with_capacity(effective_threads);
        for (cursor, col_idx) in cursors {
            let rx = work_rx.clone();
            let tx = result_tx.clone();
            worker_handles.push(s.spawn(move || -> Result<()> {
                worker_loop(cursor, col_idx, rx, tx, config, min_mapq_val, id_src)
            }));
        }
        // Drop our copies so only workers hold channel ends.
        drop(work_rx);
        drop(result_tx);

        // Sender thread — feeds work items to workers.
        s.spawn(move || {
            for item in work_items {
                if work_tx.send(item).is_err() {
                    break; // Workers died — stop sending.
                }
            }
            // work_tx dropped here, closing the work channel.
        });

        // Collector — runs on the main thread, writes chunks in reference order.
        collect_ordered_chunks(&result_rx, writer)?;

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
/// as the next expected chunk becomes available.
fn collect_ordered_chunks(
    result_rx: &Receiver<ResultChunk>,
    writer: &mut impl Write,
) -> Result<()> {
    let mut next_ref_idx: usize = 0;
    // Per-reference: next chunk_seq we expect to write.
    let mut next_chunk_seq: BTreeMap<usize, usize> = BTreeMap::new();
    // Per-reference: chunk_seq of the is_last chunk (once seen).
    let mut last_chunk_seq: BTreeMap<usize, usize> = BTreeMap::new();
    // Buffered chunks waiting to be written, keyed by (ref_idx, chunk_seq).
    let mut pending: BTreeMap<(usize, usize), Vec<u8>> = BTreeMap::new();

    for chunk in result_rx {
        if chunk.is_last {
            last_chunk_seq.insert(chunk.ref_idx, chunk.chunk_seq);
        }
        pending.insert((chunk.ref_idx, chunk.chunk_seq), chunk.data);

        // Flush all chunks that are next in order.
        loop {
            let expected_seq = next_chunk_seq.get(&next_ref_idx).copied().unwrap_or(0);
            if let Some(data) = pending.remove(&(next_ref_idx, expected_seq)) {
                if !data.is_empty() {
                    writer.write_all(&data)?;
                }
                if last_chunk_seq.get(&next_ref_idx) == Some(&expected_seq) {
                    // This reference is complete — advance to the next one.
                    next_chunk_seq.remove(&next_ref_idx);
                    last_chunk_seq.remove(&next_ref_idx);
                    next_ref_idx += 1;
                } else {
                    next_chunk_seq.insert(next_ref_idx, expected_seq + 1);
                }
            } else {
                break;
            }
        }
    }

    Ok(())
}

/// Worker loop: process references from the work channel, send chunked results back.
fn worker_loop(
    cursor: VCursor,
    col_idx: AlignColumnIndices,
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<ResultChunk>,
    config: &AlignConfig<'_>,
    min_mapq_val: i32,
    id_src: AlignIdSrc,
) -> Result<()> {
    let mut record_buf = Vec::with_capacity(1024);
    let mut mate_cache = MateCache::new();
    let align_mgr = AlignMgr::make_read().context("worker: failed to create AlignMgr")?;

    while let Ok(item) = work_rx.recv() {
        let mut psi = align_mgr
            .make_placement_set_iterator()
            .context("worker: failed to create PlacementSetIterator")?;

        // Create a PlacementIterator for this single reference.
        let pi = match PlacementIterator::make(
            &item.ref_obj,
            0,
            item.ref_len,
            min_mapq_val,
            &cursor,
            id_src,
        ) {
            Ok(pi) => pi,
            Err(e) if e.is_done() => {
                // No alignments on this reference — send empty final chunk.
                result_tx
                    .send(ResultChunk {
                        ref_idx: item.idx,
                        chunk_seq: 0,
                        data: Vec::new(),
                        is_last: true,
                    })
                    .map_err(|_| anyhow::anyhow!("result channel closed"))?;
                continue;
            }
            Err(e) => {
                return Err(e).context("worker: failed to create PlacementIterator");
            }
        };

        match psi.add_placement_iterator(pi) {
            Ok(_) => {}
            Err(e) => {
                return Err(e).context("worker: failed to add PlacementIterator");
            }
        }

        // Process records, sending chunks when the buffer exceeds CHUNK_SIZE.
        let mut output_buf = Vec::with_capacity(CHUNK_SIZE);
        let mut chunk_seq = 0usize;
        emit_placement_records(
            &mut psi,
            &cursor,
            &col_idx,
            &mut record_buf,
            &mut mate_cache,
            config.use_seqid,
            config.opts,
            |rec| {
                output_buf.extend_from_slice(rec);
                if output_buf.len() >= CHUNK_SIZE {
                    result_tx
                        .send(ResultChunk {
                            ref_idx: item.idx,
                            chunk_seq,
                            data: std::mem::replace(
                                &mut output_buf,
                                Vec::with_capacity(CHUNK_SIZE),
                            ),
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
            .send(ResultChunk { ref_idx: item.idx, chunk_seq, data: output_buf, is_last: true })
            .map_err(|_| anyhow::anyhow!("result channel closed"))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use crossbeam_channel::bounded;

    use super::*;

    /// Helper: send chunks into a channel and collect the output via `collect_ordered_chunks`.
    fn run_collector(chunks: Vec<ResultChunk>) -> Vec<u8> {
        let (tx, rx) = bounded::<ResultChunk>(chunks.len() + 1);
        for chunk in chunks {
            tx.send(chunk).unwrap();
        }
        drop(tx);

        let mut output = Vec::new();
        collect_ordered_chunks(&rx, &mut output).unwrap();
        output
    }

    #[test]
    fn test_single_ref_single_chunk() {
        let output = run_collector(vec![ResultChunk {
            ref_idx: 0,
            chunk_seq: 0,
            data: b"line1\n".to_vec(),
            is_last: true,
        }]);
        assert_eq!(output, b"line1\n");
    }

    #[test]
    fn test_single_ref_multiple_chunks() {
        let output = run_collector(vec![
            ResultChunk { ref_idx: 0, chunk_seq: 0, data: b"aaa\n".to_vec(), is_last: false },
            ResultChunk { ref_idx: 0, chunk_seq: 1, data: b"bbb\n".to_vec(), is_last: false },
            ResultChunk { ref_idx: 0, chunk_seq: 2, data: b"ccc\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"aaa\nbbb\nccc\n");
    }

    #[test]
    fn test_multiple_refs_one_chunk_each() {
        let output = run_collector(vec![
            ResultChunk { ref_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 1, chunk_seq: 0, data: b"ref1\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 2, chunk_seq: 0, data: b"ref2\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref1\nref2\n");
    }

    #[test]
    fn test_out_of_order_refs() {
        // ref 1 arrives before ref 0 — should buffer ref 1 and write ref 0 first.
        let output = run_collector(vec![
            ResultChunk { ref_idx: 1, chunk_seq: 0, data: b"ref1\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref1\n");
    }

    #[test]
    fn test_multiple_refs_multiple_chunks_out_of_order() {
        // Interleaved chunks from two references, arriving out of order.
        let output = run_collector(vec![
            ResultChunk { ref_idx: 1, chunk_seq: 0, data: b"r1c0\n".to_vec(), is_last: false },
            ResultChunk { ref_idx: 0, chunk_seq: 1, data: b"r0c1\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 1, chunk_seq: 1, data: b"r1c1\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 0, chunk_seq: 0, data: b"r0c0\n".to_vec(), is_last: false },
        ]);
        assert_eq!(output, b"r0c0\nr0c1\nr1c0\nr1c1\n");
    }

    #[test]
    fn test_empty_ref() {
        // A reference with only an empty final chunk should produce no output.
        let output = run_collector(vec![
            ResultChunk { ref_idx: 0, chunk_seq: 0, data: b"ref0\n".to_vec(), is_last: true },
            ResultChunk { ref_idx: 1, chunk_seq: 0, data: Vec::new(), is_last: true },
            ResultChunk { ref_idx: 2, chunk_seq: 0, data: b"ref2\n".to_vec(), is_last: true },
        ]);
        assert_eq!(output, b"ref0\nref2\n");
    }

    #[test]
    fn test_large_chunk_sequence() {
        // 10 chunks per reference × 3 refs.
        let mut chunks = Vec::new();
        for ref_idx in 0..3 {
            for seq in 0..10 {
                chunks.push(ResultChunk {
                    ref_idx,
                    chunk_seq: seq,
                    data: format!("r{ref_idx}c{seq}\n").into_bytes(),
                    is_last: seq == 9,
                });
            }
        }
        let output = run_collector(chunks);
        let expected: String =
            (0..3).flat_map(|r| (0..10).map(move |c| format!("r{r}c{c}\n"))).collect();
        assert_eq!(output, expected.as_bytes());
    }
}
