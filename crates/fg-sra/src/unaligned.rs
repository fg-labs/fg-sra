//! Unaligned read processing via SEQUENCE table scan.
//!
//! Scans the SEQUENCE table for spots that are fully or partially unaligned,
//! producing SAM records for unaligned reads.
//!
//! When the table stores its unaligned reads' bases (`CMP_READ`, as aligned
//! archives do), spots are converted on several threads, each reading
//! through its own cursor a blob at a time, and written in spot order. The
//! bases come from `CMP_READ`, not the virtual `READ`, which for a partly
//! aligned spot reaches into the reference, whose cache isn't thread-safe.
//! Other tables are converted on one thread, through `READ`.

use std::ops::RangeInclusive;

use anyhow::{Context, Result, anyhow};
use crossbeam_channel::{Receiver, Sender, bounded};
use fg_sra_vdb::cursor::{BlobColumn, VCursor};
use fg_sra_vdb::database::{VDatabase, VTable};

use crate::archive::marks_rna;
use crate::output::OutputWriter;
use crate::progress::ProgressLogger;
use crate::record::{
    FormatOptions, READ_TYPE_BIOLOGICAL, UnalignedColumns, format_unaligned_record,
};
use crate::refstore::{CHARSET_4NA, CHARSET_4NA_RNA};

/// Spots per batch converted by one worker.
const BATCH_SPOTS: i64 = 16 * 1024;

/// Batches each worker may have queued or in flight, bounding memory.
const BATCHES_PER_WORKER: usize = 4;

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
    /// The unaligned reads' bases, end to end, as 4na codes.
    pub const CMP_READ_4NA: &str = "(INSDC:4na:bin)CMP_READ";
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

/// Process unaligned reads from the SEQUENCE table, on up to `num_threads`
/// threads when it stores `CMP_READ` (see the module docs).
///
/// For each spot, iterates over reads and outputs those that are:
/// - Biological (`READ_TYPE` has biological bit set)
/// - Unaligned (`PRIMARY_ALIGNMENT_ID` == 0 for that read)
/// - Non-empty (`READ_LEN` > 0)
pub fn process_unaligned_reads(
    db: &VDatabase,
    writer: &mut OutputWriter,
    opts: &FormatOptions<'_>,
    unaligned_spots_only: bool,
    num_threads: usize,
    progress: &ProgressLogger,
) -> Result<()> {
    let table = db.open_table_read("SEQUENCE").context("failed to open SEQUENCE table")?;
    if stores_cmp_read(&table)? {
        let batch = Batching { threads: num_threads, spots: BATCH_SPOTS };
        process_in_parallel(&table, writer, opts, unaligned_spots_only, batch, progress)
    } else {
        process_serially(db, writer, opts, unaligned_spots_only, progress)
    }
}

/// Whether `table` stores its unaligned reads' bases as `CMP_READ`: a physical column, not
/// just one the schema declares (a table that stores `READ` declares it too).
fn stores_cmp_read(table: &VTable) -> Result<bool> {
    let physical = table.list_physical_columns().context("failed to list SEQUENCE columns")?;
    Ok(physical.iter().any(|column| column == "CMP_READ"))
}

/// Fill `emitted` with the indices of the spot's reads to write: unaligned, biological and
/// non-empty; none if `unaligned_spots_only` and any of its reads is aligned.
fn reads_to_emit(
    primary_ids: &[i64],
    read_types: &[u8],
    read_lens: &[u32],
    unaligned_spots_only: bool,
    emitted: &mut Vec<usize>,
) {
    emitted.clear();
    if unaligned_spots_only && primary_ids.iter().any(|&id| id != 0) {
        return;
    }
    emitted.extend((0..primary_ids.len()).filter(|&i| {
        primary_ids[i] == 0
            && (read_types.get(i).copied().unwrap_or(0) & READ_TYPE_BIOLOGICAL) != 0
            && read_lens.get(i).copied().unwrap_or(0) != 0
    }));
}

/// Convert unaligned reads on one thread, reading bases through `READ`.
fn process_serially(
    db: &VDatabase,
    writer: &mut OutputWriter,
    opts: &FormatOptions<'_>,
    unaligned_spots_only: bool,
    progress: &ProgressLogger,
) -> Result<()> {
    let (cursor, idx) = setup_seq_cursor(db)?;

    let (first_row, row_count) =
        cursor.id_range(idx.read).context("failed to get SEQUENCE row range")?;

    let mut buf = Vec::with_capacity(512);
    let mut emit_indices = Vec::new();

    for row_id in first_row..first_row + row_count as i64 {
        let primary_ids = cursor.read_i64_slice(row_id, idx.primary_alignment_id)?;
        let read_types = cursor.read_u8_slice(row_id, idx.read_type)?;
        let read_starts = cursor.read_i32_slice(row_id, idx.read_start)?;
        let read_lens = cursor.read_u32_slice(row_id, idx.read_len)?;
        let read_filters = cursor.read_u8_slice(row_id, idx.read_filter)?;

        reads_to_emit(
            &primary_ids,
            &read_types,
            &read_lens,
            unaligned_spots_only,
            &mut emit_indices,
        );
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
            writer.write_bytes(&buf)?;
            progress.record(1);
        }
    }

    progress.complete();
    Ok(())
}

/// A converted batch: its formatted records and how many there are.
type BatchResult = Result<(Vec<u8>, u64)>;

/// A batch of spots, and where to send its result.
type Job = (RangeInclusive<i64>, Sender<BatchResult>);

/// How [`process_in_parallel`] divides the work: worker threads, and spots per batch.
#[derive(Clone, Copy)]
struct Batching {
    threads: usize,
    spots: i64,
}

/// Convert unaligned reads on worker threads, in batches of spots written in spot order.
fn process_in_parallel(
    table: &VTable,
    writer: &mut OutputWriter,
    opts: &FormatOptions<'_>,
    unaligned_spots_only: bool,
    batching: Batching,
    progress: &ProgressLogger,
) -> Result<()> {
    let (first_row, row_count) = {
        let cursor = table.create_cursor_read().context("failed to create SEQUENCE cursor")?;
        let read = cursor.add_column(col::READ).context("READ")?;
        cursor.open().context("failed to open SEQUENCE cursor")?;
        cursor.id_range(read).context("failed to get SEQUENCE row range")?
    };
    let last_row = first_row + row_count as i64 - 1;
    // Each worker's cursor is made here, on the table's thread, and moved to the worker.
    // Bases are rendered as `READ` renders them, with `U` for `T` in a table marked as RNA.
    let charset = if marks_rna(table) { CHARSET_4NA_RNA } else { CHARSET_4NA };
    let readers = (0..batching.threads.max(1))
        .map(|_| SpotReader::new(table, charset))
        .collect::<Result<Vec<_>>>()?;
    let num_threads = readers.len();

    std::thread::scope(|scope| -> Result<()> {
        let (job_tx, job_rx) = bounded::<Job>(num_threads);
        let (order_tx, order_rx) = bounded(num_threads * BATCHES_PER_WORKER);
        for reader in readers {
            let job_rx = job_rx.clone();
            std::thread::Builder::new()
                .stack_size(crate::archive::VDB_THREAD_STACK_BYTES)
                .spawn_scoped(scope, move || {
                    convert_batches(reader, &job_rx, opts, unaligned_spots_only);
                })
                .context("failed to start an unaligned-read worker")?;
        }
        drop(job_rx);
        scope.spawn(move || {
            let mut start = first_row;
            while start <= last_row {
                let end = (start + batching.spots - 1).min(last_row);
                let (done_tx, done_rx) = bounded(1);
                // A closed channel means the writer stopped on an error.
                if order_tx.send(done_rx).is_err() || job_tx.send((start..=end, done_tx)).is_err() {
                    return;
                }
                start = end + 1;
            }
        });
        // Returning (on an error) drops `order_rx`, which stops the dispatcher and so the
        // workers, before the scope waits for them.
        write_in_order(&order_rx, writer, progress)
    })?;
    progress.complete();
    Ok(())
}

/// Write each batch's records in dispatch order.
fn write_in_order(
    order: &Receiver<Receiver<BatchResult>>,
    writer: &mut OutputWriter,
    progress: &ProgressLogger,
) -> Result<()> {
    for batch in order {
        let (bytes, records) = batch
            .recv()
            .map_err(|_| anyhow!("an unaligned-read worker exited without converting a batch"))??;
        writer.write_bytes(&bytes)?;
        progress.record(records);
    }
    Ok(())
}

/// Worker loop: convert each batch received, sending back its records.
///
/// A panic converting a batch is sent back as that batch's error, and ends the worker, so
/// the run fails with an error rather than aborting.
fn convert_batches(
    mut reader: SpotReader,
    jobs: &Receiver<Job>,
    opts: &FormatOptions<'_>,
    unaligned_spots_only: bool,
) {
    for (spots, done) in jobs {
        let converted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut out = Vec::new();
            let mut records = 0u64;
            for row in spots {
                records += reader.convert_spot(row, opts, unaligned_spots_only, &mut out)?;
            }
            Ok((out, records))
        }));
        let panicked = converted.is_err();
        let result =
            converted.unwrap_or_else(|_| Err(anyhow!("an unaligned-read worker panicked")));
        // The writer may have stopped on another batch's error; nothing waits for this one.
        let _ = done.send(result);
        if panicked {
            return;
        }
    }
}

/// One worker's cursor on the SEQUENCE table, its columns read a blob at a time,
/// and buffers reused from spot to spot.
struct SpotReader {
    cursor: VCursor,
    cmp_read: BlobColumn,
    quality: BlobColumn,
    spot_group: BlobColumn,
    read_start: BlobColumn,
    read_len: BlobColumn,
    read_type: BlobColumn,
    read_filter: BlobColumn,
    name: BlobColumn,
    primary_alignment_id: BlobColumn,
    /// Text for each 4na code: [`CHARSET_4NA`], or [`CHARSET_4NA_RNA`] for RNA.
    charset: &'static [u8; 16],
    buffers: SpotBuffers,
}

/// A spot's cells.
#[derive(Default)]
struct SpotBuffers {
    cmp_codes: Vec<u8>,
    cmp_text: Vec<u8>,
    quality: Vec<u8>,
    spot_group: Vec<u8>,
    read_starts: Vec<i32>,
    read_lens: Vec<u32>,
    read_types: Vec<u8>,
    read_filters: Vec<u8>,
    name: Vec<u8>,
    primary_ids: Vec<i64>,
    /// Indices of the reads to write.
    emitted: Vec<usize>,
    /// Where each read's bases start in `cmp_text`.
    offsets: Vec<usize>,
    record: Vec<u8>,
}

impl SpotReader {
    /// A reader with its own (uncached, as blob reads require) cursor on `table`, rendering
    /// bases with `charset`.
    fn new(table: &VTable, charset: &'static [u8; 16]) -> Result<Self> {
        let cursor = table.create_cursor_read().context("failed to create SEQUENCE cursor")?;
        let add = |name: &str| {
            cursor.add_column(name).map(BlobColumn::new).with_context(|| name.to_string())
        };
        let reader = Self {
            cmp_read: add(col::CMP_READ_4NA)?,
            quality: add(col::QUALITY)?,
            spot_group: add(col::SPOT_GROUP)?,
            read_start: add(col::READ_START)?,
            read_len: add(col::READ_LEN)?,
            read_type: add(col::READ_TYPE)?,
            read_filter: add(col::READ_FILTER)?,
            name: add(col::NAME)?,
            primary_alignment_id: add(col::PRIMARY_ALIGNMENT_ID)?,
            charset,
            buffers: SpotBuffers::default(),
            cursor,
        };
        reader.cursor.open().context("failed to open SEQUENCE cursor")?;
        Ok(reader)
    }

    /// Append the records of spot `row`'s unaligned reads to `out`, as
    /// [`process_serially`] writes them, returning how many were written.
    fn convert_spot(
        &mut self,
        row: i64,
        opts: &FormatOptions<'_>,
        unaligned_spots_only: bool,
        out: &mut Vec<u8>,
    ) -> Result<u64> {
        let Self {
            cursor,
            cmp_read,
            quality,
            spot_group,
            read_start,
            read_len,
            read_type,
            read_filter,
            name,
            primary_alignment_id,
            charset,
            buffers: b,
        } = self;
        let cursor = &*cursor;
        primary_alignment_id.read_i64_slice_into(cursor, row, &mut b.primary_ids)?;
        read_type.read_u8_slice_into(cursor, row, &mut b.read_types)?;
        read_start.read_i32_slice_into(cursor, row, &mut b.read_starts)?;
        read_len.read_u32_slice_into(cursor, row, &mut b.read_lens)?;
        read_filter.read_u8_slice_into(cursor, row, &mut b.read_filters)?;

        reads_to_emit(
            &b.primary_ids,
            &b.read_types,
            &b.read_lens,
            unaligned_spots_only,
            &mut b.emitted,
        );
        if b.emitted.is_empty() {
            return Ok(0);
        }

        cmp_read.read_u8_slice_into(cursor, row, &mut b.cmp_codes)?;
        b.cmp_text.clear();
        b.cmp_text.extend(b.cmp_codes.iter().map(|&code| charset[usize::from(code & 0x0F)]));
        quality.read_u8_slice_into(cursor, row, &mut b.quality)?;
        name.read_u8_slice_into(cursor, row, &mut b.name)?;
        spot_group.read_u8_slice_into(cursor, row, &mut b.spot_group)?;
        let bases = std::str::from_utf8(&b.cmp_text).context("CMP_READ is not text")?;
        let spot_name = String::from_utf8_lossy(&b.name);
        let group = String::from_utf8_lossy(&b.spot_group);
        let offsets = &mut b.offsets;
        unaligned_read_offsets(&b.primary_ids, &b.read_starts, &b.read_lens, bases.len(), offsets)
            .with_context(|| format!("SEQUENCE row {row}"))?;

        let num_bio_reads = b.emitted.len() as u32;
        for (bio_index, &i) in b.emitted.iter().enumerate() {
            let len = b.read_lens[i] as usize;
            let start = b.read_starts.get(i).copied().unwrap_or(0) as usize;
            let cols = UnalignedColumns {
                name: &spot_name,
                read: safe_slice_str(bases, offsets[i], offsets[i] + len),
                quality: safe_slice(&b.quality, start, start + len),
                spot_group: &group,
                read_type: b.read_types.get(i).copied().unwrap_or(0),
                read_filter: b.read_filters.get(i).copied().unwrap_or(0),
                num_bio_reads,
                bio_read_index: bio_index as u32,
            };
            format_unaligned_record(&mut b.record, &cols, opts);
            out.extend_from_slice(&b.record);
        }
        Ok(u64::from(num_bio_reads))
    }
}

/// Fill `offsets` with where each read's bases start in `CMP_READ` (of `cmp_len`
/// bases): at its `READ_START` when `CMP_READ` holds every base of the spot, and
/// otherwise, for an unaligned read, after the earlier unaligned reads, as ncbi-vdb's
/// `seq_restore_read` takes them. Aligned reads' offsets are unused.
///
/// Fails, as `seq_restore_read` does, if the unaligned reads need more bases than
/// `CMP_READ` has.
fn unaligned_read_offsets(
    primary_ids: &[i64],
    read_starts: &[i32],
    read_lens: &[u32],
    cmp_len: usize,
    offsets: &mut Vec<usize>,
) -> Result<()> {
    offsets.clear();
    let total: usize = read_lens.iter().map(|&len| len as usize).sum();
    if total == cmp_len {
        offsets.extend(
            (0..primary_ids.len()).map(|i| read_starts.get(i).copied().unwrap_or(0) as usize),
        );
        return Ok(());
    }
    let mut next = 0;
    for (i, &id) in primary_ids.iter().enumerate() {
        offsets.push(next);
        if id <= 0 {
            next += read_lens.get(i).copied().unwrap_or(0) as usize;
        }
    }
    anyhow::ensure!(
        next <= cmp_len,
        "the unaligned reads have {next} bases but CMP_READ only {cmp_len}; the archive is \
         inconsistent"
    );
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

#[cfg(test)]
mod tests {
    use super::*;

    /// [`unaligned_read_offsets`] into a new vector.
    fn offsets_of(ids: &[i64], starts: &[i32], lens: &[u32], cmp_len: usize) -> Result<Vec<usize>> {
        let mut offsets = vec![99]; // Cleared first.
        unaligned_read_offsets(ids, starts, lens, cmp_len, &mut offsets)?;
        Ok(offsets)
    }

    #[test]
    fn test_unaligned_read_offsets_when_cmp_read_holds_every_base() {
        // CMP_READ holds the whole spot: reads start at READ_START, aligned or not.
        assert_eq!(offsets_of(&[7, 0], &[0, 4], &[4, 3], 7).unwrap(), [0, 4]);
    }

    #[test]
    fn test_unaligned_read_offsets_skip_aligned_reads() {
        // Mate 1 aligned (not in CMP_READ), a technical read and mate 2 unaligned.
        let offsets = offsets_of(&[7, 0, 0], &[0, 5, 8], &[5, 3, 4], 7).unwrap();
        assert_eq!(offsets[1..], [0, 3]);
        // Only the second read aligned.
        let offsets = offsets_of(&[0, 9, 0], &[0, 2, 6], &[2, 4, 3], 5).unwrap();
        assert_eq!((offsets[0], offsets[2]), (0, 2));
    }

    #[test]
    fn test_unaligned_read_offsets_fail_when_cmp_read_is_short() {
        let err = offsets_of(&[7, 0, 0], &[0, 5, 8], &[5, 3, 4], 6).unwrap_err();
        assert!(err.to_string().contains("CMP_READ only 6"), "{err}");
    }

    #[test]
    fn test_reads_to_emit() {
        let mut emitted = vec![99];
        let bio = READ_TYPE_BIOLOGICAL;
        // Aligned, unaligned, technical, empty, unaligned.
        let (ids, types, lens) = ([5, 0, 0, 0, 0], [bio, bio, 0, bio, bio], [3, 3, 3, 0, 3]);
        reads_to_emit(&ids, &types, &lens, false, &mut emitted);
        assert_eq!(emitted, [1, 4]);
        reads_to_emit(&ids, &types, &lens, true, &mut emitted);
        assert!(emitted.is_empty(), "the spot has an aligned read");
        reads_to_emit(&[0, 0], &[bio, bio], &[3, 3], true, &mut emitted);
        assert_eq!(emitted, [0, 1]);
    }

    /// Record-formatting options for SAM output.
    fn sam_options() -> FormatOptions<'static> {
        FormatOptions {
            prefix: None,
            spot_group_in_name: false,
            xi_tag: false,
            reverse_unaligned: false,
            omit_quality: false,
            qual_quant: None,
            output_mode: crate::record::OutputMode::Sam,
            ref_name_to_id: None,
        }
    }

    /// Conversions so far: tests run on threads of one process, so each conversion writes
    /// its own file.
    static CONVERSIONS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

    /// Convert `archive`'s unaligned reads with `convert`, returning the SAM written.
    fn convert_archive(
        archive: &str,
        name: &str,
        convert: impl FnOnce(&VDatabase, &mut OutputWriter, &FormatOptions<'_>) -> Result<()>,
    ) -> Vec<u8> {
        let manager = fg_sra_vdb::manager::VdbManager::make_read().unwrap();
        let db = manager.open_db_read(archive).unwrap();
        let n = CONVERSIONS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = std::env::temp_dir()
            .join(format!("fg-sra-unaligned-{}-{name}-{n}.sam", std::process::id()));
        let mut writer = OutputWriter::from_path_with_compression(
            &path,
            crate::output::CompressionMode::None,
            1,
        )
        .unwrap();
        convert(&db, &mut writer, &sam_options()).unwrap();
        writer.finish().unwrap();
        let sam = std::fs::read(&path).unwrap();
        std::fs::remove_file(&path).ok();
        sam
    }

    /// Unaligned SAM from the serial `READ` path and the parallel `CMP_READ` path (with
    /// `batching`), for full and spots-only output.
    fn serial_and_parallel(archive: &str, batching: Batching) -> [(Vec<u8>, Vec<u8>); 2] {
        let progress = ProgressLogger::new(0, 0);
        [false, true].map(|spots_only| {
            let serial = convert_archive(archive, "serial", |db, writer, opts| {
                process_serially(db, writer, opts, spots_only, &progress)
            });
            let parallel = convert_archive(archive, "parallel", |db, writer, opts| {
                let table = db.open_table_read("SEQUENCE")?;
                assert!(stores_cmp_read(&table)?, "the archive stores CMP_READ");
                process_in_parallel(&table, writer, opts, spots_only, batching, &progress)
            });
            (serial, parallel)
        })
    }

    #[test]
    fn parallel_conversion_from_cmp_read_matches_serial_conversion_from_read() {
        let archive = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../vendor/ncbi-vdb/test/vdb/db/VDB-3418.sra");
        // Batches of 7 spots on 4 threads: many batches, finished out of order.
        let batching = Batching { threads: 4, spots: 7 };
        for (spots_only, (serial, parallel)) in
            [false, true].into_iter().zip(serial_and_parallel(archive.to_str().unwrap(), batching))
        {
            assert!(!serial.is_empty(), "no unaligned reads (spots only: {spots_only})");
            assert!(serial == parallel, "outputs differ (spots only: {spots_only})");
        }
    }

    /// The same on a real aligned archive, whose partly aligned spots take bases from the
    /// middle of `CMP_READ`. Opt-in: set `FG_SRA_TEST_ALIGNED_SRA` to an aligned archive
    /// with some unaligned mates.
    #[test]
    fn parallel_conversion_matches_serial_conversion_on_an_aligned_archive() {
        let Ok(archive) = std::env::var("FG_SRA_TEST_ALIGNED_SRA") else {
            eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA to an aligned SRA to run this test");
            return;
        };
        let [(full, full_parallel), (spots_only, spots_only_parallel)] =
            serial_and_parallel(&archive, Batching { threads: 8, spots: BATCH_SPOTS });
        assert!(full == full_parallel, "outputs differ");
        assert!(spots_only == spots_only_parallel, "spots-only outputs differ");
        assert!(full != spots_only, "the archive has no partly aligned spots to compare");
    }
}
