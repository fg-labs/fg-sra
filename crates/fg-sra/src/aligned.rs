//! Aligned read processing via row-range iteration.
//!
//! Splits alignment row ID ranges across worker threads for parallel
//! processing. Each worker directly iterates rows from a VDB cursor,
//! making parallelism independent of the number of references.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use fg_sra_vdb::cursor::VCursor;
use fg_sra_vdb::database::{VDatabase, VTable};
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

use crate::matecache::{MateCache, MateInfo};
use crate::output::OutputWriter;
use crate::progress::ProgressLogger;
use crate::record::{AlignedColumns, FormatOptions, format_aligned_record};
use crate::refstore::{MAX_REFERENCE_LOADERS, ReferenceLoaders, ReferenceStore, ref_window_len};
use crate::restore_read::restore_read;

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
    /// Explicit cursor pool size override, or `None` for the default heuristic.
    pub pool_size_override: Option<usize>,
    pub opts: &'a FormatOptions<'a>,
    /// Genomic regions to restrict output to. Empty means all references.
    pub regions: &'a [String],
    /// Opens the database again, for each thread that preloads references in
    /// parallel.
    pub reopen: &'a (dyn Fn() -> Result<VDatabase> + Sync),
}

impl AlignConfig<'_> {
    /// Compute the `ReferenceList` option flags from `primary_only`.
    fn reflist_opts(&self) -> u32 {
        // READ_4NA: reference bases are read as 4na so we can apply the CHARSET
        // map ourselves and match the `(ascii)READ` column exactly.
        let mut opts = reflist_options::READ_4NA | reflist_options::USE_PRIMARY_IDS;
        if !self.primary_only {
            opts |= reflist_options::USE_SECONDARY_IDS;
        }
        opts
    }

    /// Convert `min_mapq` to `i32` for post-read filtering.
    fn min_mapq_i32(&self) -> i32 {
        self.min_mapq.map_or(0, |m| m as i32)
    }

    /// Compute the resource pool size (number of VDB cursors).
    ///
    /// When `--pool-size` is set, uses the explicit override.  Otherwise
    /// defaults to one cursor per thread (1:1).  The result is clamped to
    /// `[1, num_work_items]`.
    fn pool_size(&self, num_work_items: usize) -> usize {
        let n = self.pool_size_override.unwrap_or(self.num_threads);
        n.min(num_work_items).max(1)
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
    // Physical columns used to reconstruct READ in fg-sra (avoids the virtual
    // READ column, whose reference sub-select is not thread-safe).
    pub const HAS_MISMATCH: &str = "(bool)HAS_MISMATCH";
    pub const MISMATCH: &str = "(INSDC:dna:text)MISMATCH";
    pub const HAS_REF_OFFSET: &str = "(bool)HAS_REF_OFFSET";
    pub const REF_OFFSET: &str = "(I32)REF_OFFSET";
    pub const SAM_QUALITY: &str = "(INSDC:quality:text:phred_33)SAM_QUALITY";
    pub const EDIT_DISTANCE: &str = "(U32)EDIT_DISTANCE";
    pub const SEQ_SPOT_GROUP: &str = "(ascii)SEQ_SPOT_GROUP";
    pub const SEQ_NAME: &str = "(ascii)SEQ_NAME";
    pub const ALIGNMENT_COUNT: &str = "(U8)ALIGNMENT_COUNT";
    pub const READ_FILTER: &str = "(INSDC:SRA:read_filter)READ_FILTER";
    pub const REF_POS: &str = "(INSDC:coord:zero)REF_POS";
    pub const MAPQ: &str = "(I32)MAPQ";
    pub const REF_NAME: &str = "(ascii)REF_NAME";
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
    has_mismatch: u32,
    mismatch: u32,
    has_ref_offset: u32,
    ref_offset: u32,
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
        has_mismatch: cursor.add_column(col::HAS_MISMATCH).context("HAS_MISMATCH")?,
        mismatch: cursor.add_column(col::MISMATCH).context("MISMATCH")?,
        has_ref_offset: cursor.add_column(col::HAS_REF_OFFSET).context("HAS_REF_OFFSET")?,
        ref_offset: cursor.add_column(col::REF_OFFSET).context("REF_OFFSET")?,
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
    cursor.read_u8_slice_into(row_id, idx.has_mismatch, &mut cols.has_mismatch)?;
    cursor.read_u8_slice_into(row_id, idx.mismatch, &mut cols.mismatch)?;
    cursor.read_u8_slice_into(row_id, idx.has_ref_offset, &mut cols.has_ref_offset)?;
    cursor.read_i32_slice_into(row_id, idx.ref_offset, &mut cols.ref_offset)?;
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

/// A contiguous run of alignment rows that all belong to one reference.
///
/// A reference normally has exactly one run, but a table loaded as several
/// separately sorted batches stores a reference's alignments in several runs
/// (one per batch).
struct RefBoundary {
    ref_idx: u32,
    ref_name: String,
    first_row: i64,
    last_row: i64,
}

/// The per-reference alignment runs of one table, in output order.
struct RefLayout {
    /// Runs grouped by reference (in `ReferenceList` order), and in table order
    /// within a reference.
    boundaries: Vec<RefBoundary>,
    /// Number of maximal row ranges sorted by `(REF_ID, REF_POS)`, i.e. the
    /// number of separately sorted batches the table was loaded as.
    num_sorted_segments: usize,
    /// Names of references whose alignments span more than one run. Output for
    /// these references is grouped but not coordinate-sorted.
    split_references: Vec<String>,
}

/// A row's sort key in a sorted alignment table: `(REF_ID, REF_POS)`.
///
/// `REF_ID` is the `REFERENCE` table row (a fixed-size chunk of one reference)
/// holding the alignment's start, so it orders rows by reference and coarsely
/// by position; `REF_POS` orders rows within a chunk.
type SortKey = (i64, i32);

/// Minimum number of rows per slice when scanning sort keys in parallel.
const MIN_SCAN_SLICE_ROWS: i64 = 1_000_000;

/// Maximum number of sorted segments accepted before a table is rejected as not
/// sorted by reference, bounding the memory spent describing its layout.
const MAX_SORTED_SEGMENTS: usize = 100_000;

/// Split `first_row..=last_row` into at most `num_slices` contiguous slices of
/// at least `min_slice_rows` rows each (fewer when the range is small).
fn slice_row_range(
    first_row: i64,
    last_row: i64,
    num_slices: usize,
    min_slice_rows: i64,
) -> Vec<(i64, i64)> {
    let total_rows = last_row - first_row + 1;
    let max_slices = (total_rows / min_slice_rows).max(1) as usize;
    let num_slices = num_slices.clamp(1, max_slices) as i64;
    let slice_rows = (total_rows + num_slices - 1) / num_slices;
    (0..num_slices)
        .map(|i| first_row + i * slice_rows)
        .take_while(|&start| start <= last_row)
        .map(|start| (start, (start + slice_rows - 1).min(last_row)))
        .collect()
}

/// Return the rows in the slice `start..=end` whose sort key is less than the
/// previous row's, for a table whose first row is `first_row`.
///
/// The row before the slice is probed too, so a descent at `start` itself is
/// found; the table's first row has no predecessor and is never a descent.
/// Fails if the slice has more than `max_descents` descents.
fn find_descents_in_slice(
    first_row: i64,
    start: i64,
    end: i64,
    max_descents: usize,
    mut probe: impl FnMut(i64) -> Result<SortKey>,
) -> Result<Vec<i64>> {
    let mut prev = if start > first_row { Some(probe(start - 1)?) } else { None };
    let mut descents = Vec::new();
    for row in start..=end {
        let key = probe(row)?;
        if prev.is_some_and(|p| key < p) {
            anyhow::ensure!(
                descents.len() < max_descents,
                "alignment table is not sorted by reference (more than {max_descents} sorted \
                 row ranges); converting such tables is not supported"
            );
            descents.push(row);
        }
        prev = Some(key);
    }
    Ok(descents)
}

/// Find the maximal row ranges of `first_row..=last_row` that are sorted by
/// [`SortKey`], scanning one slice per probe in parallel.
///
/// `slices` pairs each slice (from [`slice_row_range`]) with a probe returning a
/// row's sort key; each probe runs on its own thread, so it must own any VDB
/// cursor it reads.
fn find_sorted_segments<P>(
    first_row: i64,
    last_row: i64,
    slices: Vec<((i64, i64), P)>,
) -> Result<Vec<(i64, i64)>>
where
    P: FnMut(i64) -> Result<SortKey> + Send,
{
    let descents = std::thread::scope(|s| -> Result<Vec<i64>> {
        let handles: Vec<_> = slices
            .into_iter()
            .map(|((start, end), probe)| {
                s.spawn(move || {
                    find_descents_in_slice(first_row, start, end, MAX_SORTED_SEGMENTS, probe)
                })
            })
            .collect();
        // Slices are ascending and joined in order, so the result is sorted.
        let mut all = Vec::new();
        for handle in handles {
            all.extend(handle.join().expect("sort-key scan thread panicked")?);
            anyhow::ensure!(
                all.len() < MAX_SORTED_SEGMENTS,
                "alignment table is not sorted by reference (more than {MAX_SORTED_SEGMENTS} \
                 sorted row ranges); converting such tables is not supported"
            );
        }
        Ok(all)
    })?;

    let mut segments = Vec::with_capacity(descents.len() + 1);
    let mut start = first_row;
    for descent in descents {
        segments.push((start, descent - 1));
        start = descent;
    }
    segments.push((start, last_row));
    Ok(segments)
}

/// Find the per-reference runs within each sorted segment.
///
/// `probe` returns a row's reference name. Within a sorted segment each
/// reference's rows are contiguous, so each run's end is found by binary search.
/// Returns `(name, first_row, last_row)` runs in row order.
fn find_runs_in_segments(
    segments: &[(i64, i64)],
    mut probe: impl FnMut(i64) -> Result<String>,
) -> Result<Vec<(String, i64, i64)>> {
    let mut runs = Vec::new();
    for &(start, end) in segments {
        let mut current_start = start;
        while current_start <= end {
            let current_name = probe(current_start)?;
            let mut lo = current_start;
            let mut hi = end;
            while lo < hi {
                let mid = lo + (hi - lo + 1) / 2;
                if probe(mid)? == current_name {
                    lo = mid;
                } else {
                    hi = mid - 1;
                }
            }
            runs.push((current_name, current_start, lo));
            current_start = lo + 1;
        }
    }
    Ok(runs)
}

/// Group runs by reference in `ReferenceList` order, keeping table order within
/// a reference. Returns the grouped runs and the names of references with more
/// than one run.
///
/// A sorted table is already in this order, and for a table of several sorted
/// batches this restores reference order, so only references whose alignments
/// span several batches are left unsorted.
fn group_runs_by_reference(mut runs: Vec<RefBoundary>) -> (Vec<RefBoundary>, Vec<String>) {
    // Stable sort: runs of the same reference keep their table order.
    runs.sort_by_key(|run| run.ref_idx);
    let split = runs
        .chunk_by(|a, b| a.ref_idx == b.ref_idx)
        .filter(|group| group.len() > 1)
        .map(|group| group[0].ref_name.clone())
        .collect();
    (runs, split)
}

/// Discover the per-reference alignment runs of `table`.
///
/// Scans every row's [`SortKey`] (split across up to `num_threads` cursors) for
/// the table's sorted segments, then binary-searches `REF_NAME` within each
/// segment for the per-reference runs. The concurrent scan is safe because each
/// thread owns a cursor that was created and opened on this thread, and
/// `REF_ID`/`REF_POS` are physical or computed per row from physical columns
/// without touching the `REFERENCE` table's sequence.
fn find_ref_boundaries(
    table: &VTable,
    reflist: &ReferenceList,
    use_seqid: bool,
    num_threads: usize,
) -> Result<RefLayout> {
    let cursor = table.create_cursor_read().context("failed to create boundary cursor")?;
    let ref_name_col = cursor.add_column(col::REF_NAME).context("REF_NAME")?;
    cursor.open().context("failed to open boundary cursor")?;

    let (first_row, total_count) = cursor.id_range(0).context("failed to get id range")?;
    if total_count == 0 {
        return Ok(RefLayout {
            boundaries: Vec::new(),
            num_sorted_segments: 0,
            split_references: Vec::new(),
        });
    }
    let last_row = first_row + total_count as i64 - 1;

    let mut slices = Vec::new();
    for slice in slice_row_range(first_row, last_row, num_threads, MIN_SCAN_SLICE_ROWS) {
        let scan = table.create_cursor_read().context("failed to create sort-key cursor")?;
        let ref_id_col = scan.add_column(col::REF_ID).context("REF_ID")?;
        let ref_pos_col = scan.add_column(col::REF_POS).context("REF_POS")?;
        scan.open().context("failed to open sort-key cursor")?;
        let probe = move |row: i64| -> Result<SortKey> {
            Ok((scan.read_i64(row, ref_id_col)?, scan.read_coord_zero(row, ref_pos_col)?))
        };
        slices.push((slice, probe));
    }
    let segments = find_sorted_segments(first_row, last_row, slices)?;
    let runs = find_runs_in_segments(&segments, |row| Ok(cursor.read_str(row, ref_name_col)?))?;

    let mut boundaries = Vec::with_capacity(runs.len());
    for (name, first, last) in runs {
        let ref_obj = reflist
            .find(&name)
            .with_context(|| format!("failed to find reference '{name}' in reference list"))?;
        let ref_idx = ref_obj.idx()?;
        let ref_name = if use_seqid { ref_obj.seq_id()? } else { name };
        boundaries.push(RefBoundary { ref_idx, ref_name, first_row: first, last_row: last });
    }
    let (boundaries, split_references) = group_runs_by_reference(boundaries);

    Ok(RefLayout { boundaries, num_sorted_segments: segments.len(), split_references })
}

/// The BAM `ref_id` of an output reference: its index among the output header's
/// `@SQ` lines, looked up by name in `ref_name_to_id` (see
/// [`crate::header::build_ref_id_map`], which also maps each reference's
/// alternate name). Without a header map (no header written) falls back to the
/// `ReferenceList` index.
fn resolve_bam_ref_id(
    ref_name_to_id: Option<&HashMap<String, i32>>,
    ref_name: &str,
    ref_idx: u32,
) -> Result<i32> {
    let Some(map) = ref_name_to_id else {
        return Ok(ref_idx as i32);
    };
    map.get(ref_name).copied().with_context(|| {
        format!("reference '{ref_name}' has alignments but no @SQ line in the output header")
    })
}

/// Set each work item's BAM `ref_id` (see [`resolve_bam_ref_id`]).
///
/// Runs on the final work items, so only references that are actually output
/// (e.g. those selected by `--aligned-region`) need an `@SQ` line.
fn assign_bam_ref_ids(
    work_items: &mut [RowRangeWorkItem],
    ref_name_to_id: Option<&HashMap<String, i32>>,
) -> Result<()> {
    for item in work_items {
        item.bam_ref_id = resolve_bam_ref_id(ref_name_to_id, &item.ref_name, item.ref_idx)?;
    }
    Ok(())
}

// ── Work item construction ───────────────────────────────────────────────

/// A unit of work: a contiguous range of alignment row IDs from one reference.
#[derive(Clone, Debug)]
struct RowRangeWorkItem {
    /// Contiguous index into the work list, used for ordered output.
    order_idx: usize,
    /// Reference index into the `ReferenceList` (for the preloaded reference).
    ref_idx: u32,
    /// Reference name (pre-resolved, avoids per-worker `ReferenceList`).
    ref_name: String,
    /// BAM `ref_id`: the reference's index among the output header's `@SQ` lines.
    bam_ref_id: i32,
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
                // The no-header fallback; `assign_bam_ref_ids` resolves it against
                // the output header.
                bam_ref_id: boundary.ref_idx as i32,
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
/// work items covering all of that reference's row runs, with a coordinate
/// filter applied post-read. A region naming a reference by its other name
/// (name vs seqid) is resolved to a `ReferenceList` index by `ref_idx_of`.
fn collect_row_range_region_work_items(
    boundaries: &[RefBoundary],
    regions: &[String],
    mut ref_idx_of: impl FnMut(&str) -> Result<u32>,
    num_threads: usize,
) -> Result<Vec<RowRangeWorkItem>> {
    // A reference may have several runs (one per sorted batch); index them all.
    let mut runs_by_name: HashMap<&str, Vec<&RefBoundary>> = HashMap::new();
    let mut runs_by_idx: HashMap<u32, Vec<&RefBoundary>> = HashMap::new();
    for b in boundaries {
        runs_by_name.entry(b.ref_name.as_str()).or_default().push(b);
        runs_by_idx.entry(b.ref_idx).or_default().push(b);
    }

    let mut entries: Vec<(&RefBoundary, Option<(i32, i32)>)> = Vec::new();

    for spec in regions {
        let region = parse_region(spec)?;

        // First match boundaries by name, then fall back to reflist.find() to
        // handle name/seqid aliasing.
        let runs = if let Some(runs) = runs_by_name.get(region.name.as_str()) {
            runs
        } else {
            // The region name might be the alternate form (name vs seqid).
            let ref_idx = ref_idx_of(&region.name)?;
            runs_by_idx
                .get(&ref_idx)
                .with_context(|| format!("no alignments found for reference: {}", region.name))?
        };

        let filter = match (region.start, region.end) {
            (Some(s), Some(e)) => Some((s as i32, e as i32)),
            _ => None,
        };

        entries.extend(runs.iter().map(|&b| (b, filter)));
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
    store: &ReferenceStore,
    state: &mut WorkerState,
    mut emit: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    for row_id in item.start_row..=item.end_row {
        read_aligned_columns(cursor, col_idx, row_id, &mut state.cols)?;

        // Post-read MAPQ filter.
        if state.cols.mapq < min_mapq {
            continue;
        }

        // Post-read region coordinate filter.
        if let Some((rs, re)) = item.region_filter
            && (state.cols.ref_pos < rs || state.cols.ref_pos >= re)
        {
            continue;
        }

        // Reconstruct READ from the preloaded reference and stored deltas (done
        // after the filters so skipped rows pay nothing).
        reconstruct_read(&mut state.cols, store, item.ref_idx, &mut state.scratch, row_id)?;

        // Resolve mate: look up cached mate info and store ours.
        // When mate has no alignment, strip paired-end flags to match
        // sam-dump's behavior (the read is output as unpaired).
        let mate_info = if state.cols.mate_align_id != 0 {
            let info = state.mate_cache.take(state.cols.mate_align_id);
            state.mate_cache.insert(row_id, MateInfo { ref_pos: state.cols.ref_pos });
            info
        } else {
            state.cols.strip_paired_flags();
            None
        };

        format_aligned_record(
            &mut state.record_buf,
            &state.cols,
            &item.ref_name,
            item.bam_ref_id,
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

/// Default reference-preload budget: bound peak memory by processing work items
/// in reference batches whose sequences sum to at most this many bytes.
///
/// The budget counts the preloaded (ASCII-mapped) reference bytes; each reference
/// is mapped in place as it is read, so it is never held twice. Loading on several
/// threads adds each loader's own reference reader and its caches on top.
const REF_PRELOAD_BUDGET_BYTES: usize = 1 << 30; // 1 GiB

/// The reference-preload budget in bytes, overridable via the
/// `FG_SRA_REF_PRELOAD_BUDGET_MB` environment variable (documented in the
/// README; useful to cap memory or to force multi-batch behavior in tests). A
/// value of 0, or one that does not parse, falls back to the default.
fn ref_preload_budget_bytes() -> usize {
    std::env::var("FG_SRA_REF_PRELOAD_BUDGET_MB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&mb| mb > 0)
        .map_or(REF_PRELOAD_BUDGET_BYTES, |mb| mb.saturating_mul(1024 * 1024))
}

/// The warning for a table whose references have alignments in several
/// separately sorted row ranges, so its output is not coordinate-sorted.
fn split_layout_warning(table_name: &str, layout: &RefLayout) -> String {
    format!(
        "{table_name} is stored as {} separately sorted row ranges; alignments for {} \
         reference(s) ({}) are grouped by reference but are not coordinate-sorted within \
         those references",
        layout.num_sorted_segments,
        layout.split_references.len(),
        abbreviate_names(&layout.split_references, 5),
    )
}

/// Join up to `max` names with commas, noting how many more were omitted.
fn abbreviate_names(names: &[String], max: usize) -> String {
    let shown = names.iter().take(max).map(String::as_str).collect::<Vec<_>>().join(", ");
    if names.len() > max { format!("{shown}, and {} more", names.len() - max) } else { shown }
}

/// Group consecutive work items into batches whose distinct references sum to at
/// most `budget_bytes`. Work items are already reference-ordered, so each batch
/// is a contiguous range covering whole references; a single reference larger
/// than the budget forms its own batch. Returns index ranges into `work_items`.
fn batch_work_items(
    work_items: &[RowRangeWorkItem],
    ref_len_of: &HashMap<u32, usize>,
    budget_bytes: usize,
) -> Vec<std::ops::Range<usize>> {
    let mut batches = Vec::new();
    let mut start = 0usize;
    let mut batch_bytes = 0usize;
    // Work items are reference-ordered and contiguous per reference, so a new
    // reference is one whose index differs from the previous item's. (In region
    // mode, overlapping regions could repeat a reference non-adjacently; that
    // only over-counts its length here, never affects correctness — preload
    // de-duplicates.)
    let mut last_ref: Option<u32> = None;
    let mut refs_in_batch = 0usize;
    for (i, w) in work_items.iter().enumerate() {
        if last_ref == Some(w.ref_idx) {
            continue;
        }
        let len = ref_len_of.get(&w.ref_idx).copied().unwrap_or(0);
        if refs_in_batch > 0 && batch_bytes + len > budget_bytes {
            batches.push(start..i);
            start = i;
            batch_bytes = 0;
            refs_in_batch = 0;
        }
        last_ref = Some(w.ref_idx);
        refs_in_batch += 1;
        batch_bytes += len;
    }
    if start < work_items.len() {
        batches.push(start..work_items.len());
    }
    batches
}

/// Whether the SECONDARY table physically stores `TMP_HAS_MISMATCH` (so
/// `HAS_MISMATCH` resolves to it rather than being generated from `READ`).
/// Returns `false` (forcing the safe serial path) if the table cannot be probed.
fn secondary_has_physical_mismatch(db: &VDatabase) -> bool {
    let Ok(table) = db.open_table_read(SECONDARY_ALIGNMENT_TABLE) else {
        return false;
    };
    let Ok(cursor) = table.create_cursor_read() else {
        return false;
    };
    // Both HAS_MISMATCH and MISMATCH are read and each falls back independently
    // to a READ-derived generator when its TMP_* physical column is absent, so
    // require both to be physically stored.
    cursor.add_column_optional("(bool)TMP_HAS_MISMATCH").is_some()
        && cursor.add_column_optional("(INSDC:dna:text)TMP_MISMATCH").is_some()
}

/// One alignment table's planned conversion: its references and work items.
struct TablePlan {
    table_name: &'static str,
    /// Lives for the whole table: its single-threaded reader is used only on the
    /// main thread, to preload each batch of references between (never during)
    /// the parallel worker phases.
    reflist: ReferenceList,
    work_items: Vec<RowRangeWorkItem>,
    /// Names of references whose alignments span several sorted row ranges.
    split_references: HashSet<String>,
}

/// Plan the conversion of one alignment table: discover its reference runs,
/// build its work items and resolve their BAM `ref_id`s. Returns `None` when
/// the table has nothing to output.
fn plan_table(
    db: &VDatabase,
    table_name: &'static str,
    config: &AlignConfig<'_>,
) -> Result<Option<TablePlan>> {
    let reflist = ReferenceList::make_database(db, config.reflist_opts(), 0)
        .context("failed to create ReferenceList")?;
    let table = db.open_table_read(table_name).context("failed to open alignment table")?;
    let layout = find_ref_boundaries(&table, &reflist, config.use_seqid, config.num_threads)?;
    if !layout.split_references.is_empty() {
        eprintln!("[layout] {}", split_layout_warning(table_name, &layout));
    }
    let mut work_items = if config.regions.is_empty() {
        collect_row_range_work_items(&layout.boundaries, config.num_threads)
    } else {
        let ref_idx_of = |name: &str| -> Result<u32> {
            let ref_obj =
                reflist.find(name).with_context(|| format!("reference not found: {name}"))?;
            ref_obj.idx().context("failed to get reference index")
        };
        collect_row_range_region_work_items(
            &layout.boundaries,
            config.regions,
            ref_idx_of,
            config.num_threads,
        )?
    };
    if work_items.is_empty() {
        return Ok(None);
    }
    assign_bam_ref_ids(&mut work_items, config.opts.ref_name_to_id)?;
    let split_references = layout.split_references.into_iter().collect();
    Ok(Some(TablePlan { table_name, reflist, work_items, split_references }))
}

/// The planned conversion of every alignment table that has output.
pub struct AlignedPlan {
    tables: Vec<TablePlan>,
}

impl AlignedPlan {
    /// Whether the aligned records will be written in coordinate order with
    /// respect to a header whose `@SQ` order is `header_ref_ids` (reference
    /// name to `@SQ` index, e.g. from [`crate::header::build_ref_id_map`]).
    ///
    /// Conservative: true only when a single table has output and its work
    /// items are coordinate-sorted (see [`work_items_are_coordinate_sorted`]).
    pub fn is_coordinate_sorted(&self, header_ref_ids: &HashMap<String, i32>) -> bool {
        tables_are_coordinate_sorted(
            self.tables.iter().map(|t| (&t.split_references, t.work_items.as_slice())),
            header_ref_ids,
        )
    }
}

/// See [`AlignedPlan::is_coordinate_sorted`]; `tables` yields each planned
/// table's `(split_references, work_items)`.
fn tables_are_coordinate_sorted<'a>(
    mut tables: impl Iterator<Item = (&'a HashSet<String>, &'a [RowRangeWorkItem])>,
    header_ref_ids: &HashMap<String, i32>,
) -> bool {
    match (tables.next(), tables.next()) {
        (None, _) => true,
        (Some((split_references, work_items)), None) => {
            work_items_are_coordinate_sorted(work_items, split_references, header_ref_ids)
        }
        _ => false, // Secondary alignments follow all primary alignments.
    }
}

/// Whether `work_items` emit records in coordinate order: no item covers a
/// reference in `split_references` (whose runs are grouped, not merged), each
/// reference's items are contiguous and visit rows in strictly increasing order
/// (within one sorted run rows are ordered by position; a repeated or
/// overlapping region revisits rows), and references appear in strictly
/// increasing `@SQ` order. A reference missing from `header_ref_ids` makes the
/// order unknown (false).
fn work_items_are_coordinate_sorted(
    work_items: &[RowRangeWorkItem],
    split_references: &HashSet<String>,
    header_ref_ids: &HashMap<String, i32>,
) -> bool {
    let mut prev: Option<(i32, i64)> = None; // (header ref id, end_row)
    for item in work_items {
        if split_references.contains(&item.ref_name) {
            return false;
        }
        let Some(&ref_id) = header_ref_ids.get(&item.ref_name) else {
            return false;
        };
        if let Some((prev_ref_id, prev_end)) = prev {
            let in_order = if ref_id == prev_ref_id {
                item.start_row > prev_end
            } else {
                ref_id > prev_ref_id
            };
            if !in_order {
                return false;
            }
        }
        prev = Some((ref_id, item.end_row));
    }
    true
}

/// Plan every alignment table that will be output (PRIMARY, then SECONDARY
/// unless `primary_only`), so a table that cannot be converted fails before any
/// output is written and the output order is known in advance.
pub fn plan_aligned_tables(db: &VDatabase, config: &AlignConfig<'_>) -> Result<AlignedPlan> {
    let mut tables = Vec::new();
    tables.extend(plan_table(db, PRIMARY_ALIGNMENT_TABLE, config)?);
    if !config.primary_only && db.has_table(SECONDARY_ALIGNMENT_TABLE) {
        tables.extend(plan_table(db, SECONDARY_ALIGNMENT_TABLE, config)?);
    }
    Ok(AlignedPlan { tables })
}

/// Process all aligned reads planned by [`plan_aligned_tables`], writing SAM/BAM
/// records, one table after another.
pub fn process_aligned_plan(
    db: &VDatabase,
    plan: &AlignedPlan,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    progress_interval: u64,
) -> Result<()> {
    // One set of loaders serves every table and batch, so the archive is opened again
    // at most once per loader.
    let mut loaders =
        ReferenceLoaders::new(config.num_threads.min(MAX_REFERENCE_LOADERS), config.reopen);
    for table in &plan.tables {
        process_table_plan(db, table, writer, config, &mut loaders, progress_interval)?;
    }
    Ok(())
}

/// Reference loaders reading through the archive reopened by [`AlignConfig::reopen`].
type TosamLoaders<'a> = ReferenceLoaders<&'a (dyn Fn() -> Result<VDatabase> + Sync)>;

/// Process one planned alignment table, preloading references batch by batch.
fn process_table_plan(
    db: &VDatabase,
    plan: &TablePlan,
    writer: &mut OutputWriter,
    config: &AlignConfig<'_>,
    loaders: &mut TosamLoaders<'_>,
    progress_interval: u64,
) -> Result<()> {
    let TablePlan { table_name, reflist, work_items, .. } = plan;
    let table_name = *table_name;

    // Reference lengths for the distinct references, to size preload batches.
    let mut ref_len_of: HashMap<u32, usize> = HashMap::new();
    for w in work_items {
        if let std::collections::hash_map::Entry::Vacant(entry) = ref_len_of.entry(w.ref_idx) {
            let len = reflist.get(w.ref_idx)?.seq_length()? as usize;
            entry.insert(len);
        }
    }
    let batches = batch_work_items(work_items, &ref_len_of, ref_preload_budget_bytes());
    // Only the SECONDARY table can force serial processing. In the align
    // schema (align.vschema), PRIMARY_ALIGNMENT stores HAS_MISMATCH/MISMATCH
    // as physical columns (`physical column <INSDC:4na:bin> ... .MISMATCH`),
    // so reading them never re-enters the reference sub-select — PRIMARY is
    // race-free on the parallel path by construction. SECONDARY, by
    // contrast, backs those columns with a TMP_* "hack" column that, when
    // absent, is generated from a RAW_READ sub-select into PRIMARY, which
    // does re-enter the unsynchronized reference reconstruction. Detect that
    // case and process such a SECONDARY table single-threaded.
    let force_serial =
        table_name == SECONDARY_ALIGNMENT_TABLE && !secondary_has_physical_mismatch(db);
    if force_serial && config.num_threads > 1 {
        eprintln!(
            "[preload] SECONDARY_ALIGNMENT lacks physical mismatch columns; \
             processing it single-threaded to avoid a libncbi-vdb data race"
        );
    }
    // One concise summary per table (orchestration layer, not the data
    // module): distinct references and their total size, plus batch count.
    let total_ref_bytes: usize = ref_len_of.values().sum();
    eprintln!(
        "[preload] {table_name}: {} reference(s), {:.1} MiB, {} batch(es)",
        ref_len_of.len(),
        total_ref_bytes as f64 / (1024.0 * 1024.0),
        batches.len(),
    );
    let progress = ProgressLogger::new(work_items.len() as u32, progress_interval);

    for range in batches {
        let batch = &work_items[range];
        let ref_indices: Vec<u32> = batch.iter().map(|w| w.ref_idx).collect();
        let store = loaders.preload(reflist, &ref_indices)?;
        // order_idx is globally contiguous, so the batch's first item's index
        // is the base the ordered collector starts from — no re-basing needed.
        let base_order = batch.first().map_or(0, |w| w.order_idx);
        if config.num_threads <= 1 || force_serial {
            process_table_sequential(db, table_name, batch, &store, writer, config, &progress)?;
        } else {
            process_table_parallel(
                db, table_name, batch, base_order, &store, writer, config, &progress,
            )?;
        }
    }

    progress.complete();
    Ok(())
}

/// Process a single alignment table sequentially using row-range iteration.
fn process_table_sequential(
    db: &VDatabase,
    table_name: &str,
    work_items: &[RowRangeWorkItem],
    store: &ReferenceStore,
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
        scratch: ReadScratch::default(),
        current_ref_idx: u32::MAX,
    };

    for item in work_items {
        if item.ref_idx != state.current_ref_idx {
            state.mate_cache.clear();
            state.current_ref_idx = item.ref_idx;
        }
        process_row_range(
            &cursor,
            &col_idx,
            item,
            min_mapq,
            config.opts,
            store,
            &mut state,
            |rec| {
                progress.record(1);
                writer.write_bytes(rec)
            },
        )?;
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
#[allow(clippy::too_many_arguments)]
fn process_table_parallel(
    db: &VDatabase,
    table_name: &str,
    work_items: &[RowRangeWorkItem],
    base_order: usize,
    store: &ReferenceStore,
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
            worker_handles.push(s.spawn(move || -> Result<()> {
                pool_worker_loop(&channels, config, store, progress)
            }));
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
        collect_ordered_chunks(&result_rx, &buf_pool_tx, writer, base_order)?;

        // Workers have finished (result channel closed). Check for errors.
        for handle in worker_handles {
            handle.join().expect("worker thread panicked")?;
        }

        Ok(())
    })
}

/// Collect `ResultChunk`s from workers and write them in reference order.
///
/// Chunks arrive out of order from multiple workers. This function buffers
/// them and writes in strict (`order_idx`, `chunk_seq`) order, flushing as soon
/// as the next expected chunk becomes available. Written buffers are returned
/// to workers via `pool_tx` for reuse.
fn collect_ordered_chunks(
    result_rx: &Receiver<ResultChunk>,
    pool_tx: &Sender<Vec<u8>>,
    writer: &mut impl Write,
    base_order: usize,
) -> Result<()> {
    let mut next_order: usize = base_order;
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
    /// Reusable buffers for reconstructing READ from the reference window.
    scratch: ReadScratch,
    /// Tracks which reference the mate cache belongs to; cleared on change.
    current_ref_idx: u32,
}

/// Reusable buffer for [`reconstruct_read`], one per worker.
#[derive(Default)]
struct ReadScratch {
    /// Reference window, when the sub-select wraps (borrowed directly otherwise).
    window: Vec<u8>,
}

/// Reconstruct `cols.read` from the reference window and the stored deltas,
/// replacing what the virtual `READ` column would have produced.
fn reconstruct_read(
    cols: &mut AlignedColumns,
    store: &ReferenceStore,
    ref_idx: u32,
    scratch: &mut ReadScratch,
    row_id: i64,
) -> Result<()> {
    let ref_len = ref_window_len(cols.has_ref_offset.len(), &cols.ref_offset);
    let window = store
        .window(ref_idx, cols.ref_pos, ref_len, &mut scratch.window)
        .with_context(|| format!("row {row_id}: reference window"))?;
    restore_read(
        window,
        &cols.has_mismatch,
        &cols.mismatch,
        &cols.has_ref_offset,
        &cols.ref_offset,
        &[],
        &mut cols.read,
    )
    .with_context(|| format!("row {row_id}: reconstruct READ"))?;
    Ok(())
}

/// Worker loop: check out VDB resources from the pool per work item, process,
/// then return them. This bounds total VDB memory to `pool_size` × per-set cost.
fn pool_worker_loop(
    channels: &PoolWorkerChannels,
    config: &AlignConfig<'_>,
    store: &ReferenceStore,
    progress: &ProgressLogger,
) -> Result<()> {
    let mut state = WorkerState {
        record_buf: Vec::with_capacity(1024),
        mate_cache: MateCache::new(),
        cols: AlignedColumns::new(),
        scratch: ReadScratch::default(),
        current_ref_idx: u32::MAX,
    };
    let min_mapq = config.min_mapq_i32();

    // Take a buffer from the pool or allocate a new one.
    let take_buf = || -> Vec<u8> {
        channels.buf_pool_rx.try_recv().unwrap_or_else(|_| Vec::with_capacity(CHUNK_SIZE))
    };

    while let Ok(item) = channels.work_rx.recv() {
        if item.ref_idx != state.current_ref_idx {
            state.mate_cache.clear();
            state.current_ref_idx = item.ref_idx;
        }
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
            store,
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
    use crate::refstore::preload_references;

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
            pool_size_override: None,
            opts: &OPTS,
            regions: &[],
            reopen: &|| anyhow::bail!("tests do not preload references"),
        }
    }

    #[test]
    fn test_data_cursor_cache_bytes() {
        assert_eq!(DATA_CURSOR_CACHE_BYTES, 32 * 1024 * 1024);
    }

    /// Build a work item for reference `ref_idx` (row range/ordering irrelevant here).
    fn wi(ref_idx: u32) -> RowRangeWorkItem {
        wi_ordered(0, ref_idx)
    }

    /// Build a work item with an explicit `order_idx` (for ordering-sensitive tests).
    fn wi_ordered(order_idx: usize, ref_idx: u32) -> RowRangeWorkItem {
        RowRangeWorkItem {
            order_idx,
            ref_idx,
            ref_name: String::new(),
            bam_ref_id: ref_idx as i32,
            start_row: 0,
            end_row: 0,
            region_filter: None,
        }
    }

    #[test]
    fn batch_empty_work_items() {
        let lens = HashMap::new();
        assert!(batch_work_items(&[], &lens, 1000).is_empty());
    }

    #[test]
    fn batch_groups_references_under_budget() {
        // refs 0,1 (400 each) then 2 (400): budget 1000 -> [0,1] then [2].
        let items = vec![wi(0), wi(0), wi(1), wi(2)];
        let lens = HashMap::from([(0, 400), (1, 400), (2, 400)]);
        assert_eq!(batch_work_items(&items, &lens, 1000), vec![0..3, 3..4]);
    }

    #[test]
    fn batch_never_splits_a_reference() {
        // Two consecutive items of ref 0 stay in the same batch even though a
        // second reference would exceed budget.
        let items = vec![wi(0), wi(0), wi(1)];
        let lens = HashMap::from([(0, 800), (1, 800)]);
        assert_eq!(batch_work_items(&items, &lens, 1000), vec![0..2, 2..3]);
    }

    #[test]
    fn batch_oversized_reference_is_own_batch() {
        let items = vec![wi(0), wi(1), wi(2)];
        let lens = HashMap::from([(0, 100), (1, 5000), (2, 100)]);
        // ref 1 alone exceeds budget -> its own batch; 0 and 2 batched separately.
        assert_eq!(batch_work_items(&items, &lens, 1000), vec![0..1, 1..2, 2..3]);
    }

    #[test]
    fn batch_region_mode_repeated_reference_only_overcounts() {
        // Region mode can repeat a reference non-adjacently (overlapping regions).
        // The repeat re-counts ref 0's length toward the budget but must not
        // affect batching correctness — items stay in reference-adjacency order.
        let items = vec![wi(0), wi(1), wi(0)];
        let lens = HashMap::from([(0, 400), (1, 400)]);
        // 0 (400) + 1 (400) = 800 <= 1000; then the non-adjacent 0 re-counts 400
        // -> 1200 > 1000 -> new batch. Only over-counts; boundaries stay valid.
        assert_eq!(batch_work_items(&items, &lens, 1000), vec![0..2, 2..3]);
    }

    /// End-to-end contract check: the in-Rust reconstruction must reproduce the
    /// bytes of libncbi-vdb's virtual `(ascii)READ` column exactly. Opt-in via
    /// `FG_SRA_TEST_ALIGNED_SRA` (a reference-compressed aligned SRA), read
    /// single-threaded so touching the virtual `READ` is safe.
    #[test]
    fn reconstructed_read_matches_virtual_read_column() {
        use fg_sra_vdb::manager::VdbManager;
        let Ok(sra) = std::env::var("FG_SRA_TEST_ALIGNED_SRA") else {
            eprintln!(
                "skipping: set FG_SRA_TEST_ALIGNED_SRA to a reference-compressed aligned SRA"
            );
            return;
        };
        let mgr = VdbManager::make_read().expect("make_read");
        mgr.disable_pagemap_thread().ok();
        let db = mgr.open_db_read(&sra).expect("open_db");
        let opts = reflist_options::READ_4NA | reflist_options::USE_PRIMARY_IDS;
        let reflist = ReferenceList::make_database(&db, opts, 0).expect("reflist");

        let table = db.open_table_read(PRIMARY_ALIGNMENT_TABLE).expect("table");
        let cursor = table.create_cursor_read().expect("cursor");
        let has_mismatch_col = cursor.add_column(col::HAS_MISMATCH).unwrap();
        let mismatch_col = cursor.add_column(col::MISMATCH).unwrap();
        let has_ref_offset_col = cursor.add_column(col::HAS_REF_OFFSET).unwrap();
        let ref_offset_col = cursor.add_column(col::REF_OFFSET).unwrap();
        let ref_pos_col = cursor.add_column(col::REF_POS).unwrap();
        let ref_name_col = cursor.add_column(col::REF_NAME).unwrap();
        let read_col = cursor.add_column("(ascii)READ").unwrap();
        cursor.open().expect("open");

        let (first, count) = cursor.id_range(read_col).expect("id_range");
        let last = first + count.min(5000) as i64;

        // Resolve each sampled row's reference to a ReferenceObj index and preload.
        let mut name_to_idx: HashMap<String, u32> = HashMap::new();
        let mut idxs = Vec::new();
        for row in first..last {
            let name = cursor.read_str(row, ref_name_col).unwrap();
            if let std::collections::hash_map::Entry::Vacant(slot) = name_to_idx.entry(name.clone())
            {
                let idx = reflist.find(&name).expect("find ref").idx().expect("idx");
                slot.insert(idx);
                idxs.push(idx);
            }
        }
        let store = preload_references(&reflist, &idxs).expect("preload");

        let mut has_mismatch = Vec::new();
        let mut mismatch = Vec::new();
        let mut has_ref_offset = Vec::new();
        let mut ref_offset = Vec::new();
        let mut window_scratch = Vec::new();
        let mut reconstructed = Vec::new();
        for row in first..last {
            let name = cursor.read_str(row, ref_name_col).unwrap();
            let ref_idx = name_to_idx[&name];
            let ref_pos = cursor.read_coord_zero(row, ref_pos_col).unwrap();
            cursor.read_u8_slice_into(row, has_mismatch_col, &mut has_mismatch).unwrap();
            cursor.read_u8_slice_into(row, mismatch_col, &mut mismatch).unwrap();
            cursor.read_u8_slice_into(row, has_ref_offset_col, &mut has_ref_offset).unwrap();
            cursor.read_i32_slice_into(row, ref_offset_col, &mut ref_offset).unwrap();
            let ref_len = ref_window_len(has_ref_offset.len(), &ref_offset);
            let window = store.window(ref_idx, ref_pos, ref_len, &mut window_scratch).unwrap();
            restore_read(
                window,
                &has_mismatch,
                &mismatch,
                &has_ref_offset,
                &ref_offset,
                &[],
                &mut reconstructed,
            )
            .unwrap();
            let expected = cursor.read_str(row, read_col).unwrap();
            assert_eq!(reconstructed, expected.as_bytes(), "row {row} ref {name} pos {ref_pos}");
        }
    }

    /// End-to-end check of boundary discovery on a real table: the runs must
    /// tile `PRIMARY_ALIGNMENT` and each run's name must match the `REF_NAME` of
    /// its rows (checked at its ends and at evenly spaced rows). Opt-in via
    /// `FG_SRA_TEST_ALIGNED_SRA`; most informative on a table stored as several
    /// separately sorted batches.
    #[test]
    fn ref_boundaries_match_ref_name_column() {
        use fg_sra_vdb::manager::VdbManager;
        let Ok(sra) = std::env::var("FG_SRA_TEST_ALIGNED_SRA") else {
            eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA to an aligned SRA");
            return;
        };
        let mgr = VdbManager::make_read().expect("make_read");
        mgr.disable_pagemap_thread().ok();
        let db = mgr.open_db_read(&sra).expect("open_db");
        let reflist = ReferenceList::make_database(&db, 0, 0).expect("reflist");
        let table = db.open_table_read(PRIMARY_ALIGNMENT_TABLE).expect("table");
        let layout = find_ref_boundaries(&table, &reflist, false, 4).expect("find_ref_boundaries");

        let cursor = table.create_cursor_read().expect("cursor");
        let ref_name_col = cursor.add_column(col::REF_NAME).unwrap();
        cursor.open().expect("open");
        let (first, count) = cursor.id_range(0).expect("id_range");

        let mut runs: Vec<_> =
            layout.boundaries.iter().map(|b| (b.first_row, b.last_row)).collect();
        runs.sort_unstable();
        let mut expected_start = first;
        for (start, end) in runs {
            assert_eq!(start, expected_start, "runs must tile the table");
            expected_start = end + 1;
        }
        assert_eq!(expected_start, first + count as i64, "runs must cover every row");

        for b in &layout.boundaries {
            let step = ((b.last_row - b.first_row) / 100).max(1);
            let rows = (b.first_row..=b.last_row).step_by(step as usize).chain([b.last_row]);
            for row in rows {
                let name = cursor.read_str(row, ref_name_col).unwrap();
                assert_eq!(name, b.ref_name, "row {row}");
            }
        }
    }

    #[test]
    fn test_pool_size_equals_threads() {
        // Default: 1:1 cursors to threads.
        let config = test_config(1);
        assert_eq!(config.pool_size(100), 1);

        let config = test_config(4);
        assert_eq!(config.pool_size(100), 4);

        let config = test_config(8);
        assert_eq!(config.pool_size(100), 8);
    }

    #[test]
    fn test_pool_size_clamped_by_work_items() {
        // 8 threads but only 3 work items → 3.
        let config = test_config(8);
        assert_eq!(config.pool_size(3), 3);
    }

    #[test]
    fn test_pool_size_minimum_one() {
        let config = test_config(0);
        assert_eq!(config.pool_size(5), 1);
    }

    #[test]
    fn test_pool_size_override_used() {
        // Explicit override of 6 with 100 work items → 6.
        let mut config = test_config(16);
        config.pool_size_override = Some(6);
        assert_eq!(config.pool_size(100), 6);
    }

    #[test]
    fn test_pool_size_override_clamped_by_work_items() {
        // Override of 10 but only 3 work items → 3.
        let mut config = test_config(16);
        config.pool_size_override = Some(10);
        assert_eq!(config.pool_size(3), 3);
    }

    #[test]
    fn test_pool_size_override_zero_clamped_to_one() {
        // Override of 0 → clamped to 1.
        let mut config = test_config(8);
        config.pool_size_override = Some(0);
        assert_eq!(config.pool_size(100), 1);
    }

    /// Helper: send chunks into a channel and collect the output via `collect_ordered_chunks`.
    fn run_collector(chunks: Vec<ResultChunk>) -> Vec<u8> {
        run_collector_base(chunks, 0)
    }

    /// As [`run_collector`], but starting collection from `base_order` (the value
    /// a non-first preload batch passes in).
    fn run_collector_base(chunks: Vec<ResultChunk>, base_order: usize) -> Vec<u8> {
        let (tx, rx) = bounded::<ResultChunk>(chunks.len() + 1);
        let (pool_tx, _pool_rx) = bounded::<Vec<u8>>(16);
        for chunk in chunks {
            tx.send(chunk).unwrap();
        }
        drop(tx);

        let mut output = Vec::new();
        collect_ordered_chunks(&rx, &pool_tx, &mut output, base_order).unwrap();
        output
    }

    #[test]
    fn test_nonzero_base_order_orders_from_base() {
        // A later preload batch numbers its work items from a non-zero global
        // order_idx; the collector must emit them in order starting at that base
        // (not from 0). Deliver them out of order to prove the ordering holds.
        let output = run_collector_base(
            vec![
                ResultChunk { order_idx: 6, chunk_seq: 0, data: b"r6\n".to_vec(), is_last: true },
                ResultChunk { order_idx: 5, chunk_seq: 1, data: b"r5c1\n".to_vec(), is_last: true },
                ResultChunk {
                    order_idx: 5,
                    chunk_seq: 0,
                    data: b"r5c0\n".to_vec(),
                    is_last: false,
                },
            ],
            5,
        );
        assert_eq!(output, b"r5c0\nr5c1\nr6\n");
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

    // ── Reference boundary discovery tests ───────────────────────────────

    /// A synthetic alignment table: each row's `(SortKey, REF_NAME)`, row IDs from 1.
    ///
    /// Each batch is a sorted run of `(reference name, row count)` entries. Rows of
    /// one reference share a `REF_ID` chunk and have increasing positions, so keys
    /// ascend within a batch (the loader's sort order) and drop at the start of a
    /// batch that begins at an earlier reference or position.
    fn synthetic_table(batches: &[&[(&str, usize)]]) -> Vec<(SortKey, String)> {
        let ordinal = |name: &str| -> i64 {
            ["1", "2", "3", "4", "X", "GL1"].iter().position(|n| *n == name).unwrap() as i64
        };
        let mut rows = Vec::new();
        for batch in batches {
            for &(name, count) in *batch {
                for pos in 0..count {
                    rows.push(((ordinal(name), pos as i32), name.to_string()));
                }
            }
        }
        rows
    }

    /// Run the production segment and run discovery (without VDB) over a
    /// synthetic table, scanning sort keys in `num_slices` threads.
    fn discover_runs(rows: &[(SortKey, String)], num_slices: usize) -> Vec<(String, i64, i64)> {
        let (first, last) = (1, rows.len() as i64);
        let slices = slice_row_range(first, last, num_slices, 1)
            .into_iter()
            .map(|slice| (slice, |row: i64| Ok(rows[(row - 1) as usize].0)))
            .collect();
        let segments = find_sorted_segments(first, last, slices).unwrap();
        find_runs_in_segments(&segments, |row| Ok(rows[(row - 1) as usize].1.clone())).unwrap()
    }

    /// Assert the runs tile rows `1..=len` in order and label every row correctly.
    fn assert_runs_label_every_row(rows: &[(SortKey, String)], runs: &[(String, i64, i64)]) {
        let mut expected_start = 1;
        for (name, first, last) in runs {
            assert_eq!(*first, expected_start, "runs must tile the table");
            for row in *first..=*last {
                assert_eq!(&rows[(row - 1) as usize].1, name, "row {row} mislabeled");
            }
            expected_start = last + 1;
        }
        assert_eq!(expected_start, rows.len() as i64 + 1, "runs must cover every row");
    }

    #[test]
    fn test_find_descents_ignores_equal_and_rising_keys() {
        let keys = [(0, 1), (0, 1), (0, 2), (1, 0), (1, 0), (2, 5)];
        let probe = |row: i64| Ok(keys[row as usize]);
        assert!(find_descents_in_slice(0, 0, 5, 10, probe).unwrap().is_empty());
    }

    #[test]
    fn test_find_descents_reports_drops_in_ref_id_and_position() {
        // Row 2 drops in REF_ID; row 4 drops in position within the same REF_ID.
        let keys = [(3, 0), (4, 9), (1, 2), (1, 7), (1, 3)];
        let probe = |row: i64| Ok(keys[row as usize]);
        assert_eq!(find_descents_in_slice(0, 0, 4, 10, probe).unwrap(), vec![2, 4]);
    }

    #[test]
    fn test_find_descents_in_slice_detects_drop_at_slice_start() {
        // Row 3 drops below row 2; a slice starting at row 3 must still see it.
        let keys = [(0, 0), (5, 0), (6, 0), (1, 0), (2, 0)];
        let probe = |row: i64| Ok(keys[row as usize]);
        assert_eq!(find_descents_in_slice(0, 3, 4, 10, probe).unwrap(), vec![3]);
        // The table's first row has no predecessor and is never a descent.
        assert!(find_descents_in_slice(0, 0, 2, 10, probe).unwrap().is_empty());
    }

    #[test]
    fn test_find_descents_rejects_unsorted_table() {
        let keys = [(3, 0), (2, 0), (1, 0), (0, 0)];
        let probe = |row: i64| Ok(keys[row as usize]);
        // Two descents are allowed; the third fails.
        assert_eq!(find_descents_in_slice(0, 0, 2, 2, probe).unwrap(), vec![1, 2]);
        let err = find_descents_in_slice(0, 0, 3, 2, probe).unwrap_err();
        assert!(err.to_string().contains("not sorted by reference"), "{err}");
    }

    #[test]
    fn test_slice_row_range_tiles_range() {
        assert_eq!(slice_row_range(1, 10, 3, 1), vec![(1, 4), (5, 8), (9, 10)]);
        assert_eq!(slice_row_range(1, 10, 1, 1), vec![(1, 10)]);
        // More slices than rows: one row per slice.
        assert_eq!(slice_row_range(5, 7, 8, 1), vec![(5, 5), (6, 6), (7, 7)]);
        // The minimum slice size caps the number of slices.
        assert_eq!(slice_row_range(1, 10, 8, 5), vec![(1, 5), (6, 10)]);
        assert_eq!(slice_row_range(1, 9, 8, 5), vec![(1, 9)]);
        assert_eq!(slice_row_range(1, 3, 8, 5), vec![(1, 3)]);
        assert_eq!(slice_row_range(1, 10, 0, 1), vec![(1, 10)]);
    }

    #[test]
    fn test_find_sorted_segments_splits_at_descents() {
        let keys = [(0, 0), (0, 1), (1, 0), (0, 5), (1, 1), (0, 0)];
        for num_slices in 1..=keys.len() {
            let slices = slice_row_range(1, 6, num_slices, 1)
                .into_iter()
                .map(|slice| (slice, |row: i64| Ok(keys[(row - 1) as usize])))
                .collect();
            let segments = find_sorted_segments(1, 6, slices).unwrap();
            assert_eq!(segments, vec![(1, 3), (4, 5), (6, 6)], "{num_slices} slices");
        }
    }

    #[test]
    fn test_discover_runs_single_batch() {
        let rows = synthetic_table(&[&[("1", 5), ("2", 3), ("X", 1)]]);
        let runs = discover_runs(&rows, 2);
        assert_eq!(runs, [("1".into(), 1, 5), ("2".into(), 6, 8), ("X".into(), 9, 9)]);
    }

    #[test]
    fn test_discover_runs_two_sorted_batches() {
        // A table loaded as two separately sorted batches stores each reference
        // twice. Every row must be labeled with its own reference. Batch 2's "1"
        // spans the table's midpoint, so a binary search over the whole table
        // would land in it and label batch 1's "2" and "3" rows as "1".
        let rows = synthetic_table(&[
            &[("1", 10), ("2", 10), ("3", 10)],
            &[("1", 60), ("2", 5), ("3", 5), ("GL1", 2)],
        ]);
        for num_slices in 1..=8 {
            let runs = discover_runs(&rows, num_slices);
            let names: Vec<_> = runs.iter().map(|r| r.0.as_str()).collect();
            assert_eq!(names, ["1", "2", "3", "1", "2", "3", "GL1"], "{num_slices} slices");
            assert_runs_label_every_row(&rows, &runs);
        }
    }

    #[test]
    fn test_discover_runs_varied_batch_layouts() {
        let layouts: &[&[&[(&str, usize)]]] = &[
            // Later batch covers only later references (no descent).
            &[&[("1", 10), ("2", 10)], &[("3", 10), ("X", 10)]],
            // Later batch covers only earlier references.
            &[&[("3", 10), ("X", 10)], &[("1", 10), ("2", 10)]],
            // Single-row references and a one-row batch.
            &[&[("1", 1), ("2", 1), ("3", 1)], &[("2", 1)], &[("1", 1), ("X", 1)]],
            // The same reference ends one batch and restarts the next at an
            // earlier position (a descent in position only).
            &[&[("1", 10), ("2", 10)], &[("2", 10), ("3", 10)]],
            // Three batches.
            &[&[("1", 7), ("4", 3)], &[("1", 2), ("2", 9)], &[("1", 4), ("X", 6)]],
        ];
        for layout in layouts {
            let rows = synthetic_table(layout);
            for num_slices in 1..=rows.len() {
                assert_runs_label_every_row(&rows, &discover_runs(&rows, num_slices));
            }
        }
    }

    #[test]
    fn test_discover_runs_detects_batch_restart_within_a_reference() {
        // Batch 2 restarts reference "2" at position 0, which REF_ID alone would
        // not reveal; "2" must come back as two runs.
        let rows = synthetic_table(&[&[("1", 4), ("2", 4)], &[("2", 4), ("3", 4)]]);
        let runs = discover_runs(&rows, 3);
        assert_eq!(
            runs,
            [("1".into(), 1, 4), ("2".into(), 5, 8), ("2".into(), 9, 12), ("3".into(), 13, 16)]
        );
    }

    /// A boundary for reference `ref_idx` named `chr{ref_idx}`.
    fn boundary(ref_idx: u32, first_row: i64, last_row: i64) -> RefBoundary {
        RefBoundary { ref_idx, ref_name: format!("chr{ref_idx}"), first_row, last_row }
    }

    #[test]
    fn test_group_runs_by_reference() {
        // Table order: 2, 1, 2, 3, 1, 2 -> grouped 1, 1, 2, 2, 2, 3; 1 and 2 split.
        let runs = vec![
            boundary(2, 1, 10),
            boundary(1, 11, 20),
            boundary(2, 21, 30),
            boundary(3, 31, 40),
            boundary(1, 41, 50),
            boundary(2, 51, 60),
        ];
        let (grouped, split) = group_runs_by_reference(runs);
        let got: Vec<_> = grouped.iter().map(|b| (b.ref_idx, b.first_row)).collect();
        assert_eq!(got, [(1, 11), (1, 41), (2, 1), (2, 21), (2, 51), (3, 31)]);
        assert_eq!(split, ["chr1", "chr2"], "each split reference is listed once");
    }

    #[test]
    fn test_group_runs_by_reference_restores_reference_order() {
        // Disjoint batches out of reference order: no reference is split, so the
        // grouped output is fully in reference order.
        let runs = vec![boundary(3, 1, 10), boundary(4, 11, 20), boundary(0, 21, 30)];
        let (grouped, split) = group_runs_by_reference(runs);
        let got: Vec<_> = grouped.iter().map(|b| b.ref_idx).collect();
        assert_eq!(got, [0, 3, 4]);
        assert!(split.is_empty());
    }

    #[test]
    fn test_resolve_bam_ref_id_uses_header_map() {
        // Header order differs from the ReferenceList index.
        let map = HashMap::from([("chr2".to_string(), 0), ("chr1".to_string(), 1)]);
        assert_eq!(resolve_bam_ref_id(Some(&map), "chr1", 0).unwrap(), 1);
        assert_eq!(resolve_bam_ref_id(Some(&map), "chr2", 1).unwrap(), 0);
    }

    #[test]
    fn test_resolve_bam_ref_id_missing_from_header_is_error() {
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let err = resolve_bam_ref_id(Some(&map), "chrZ", 7).unwrap_err();
        assert!(err.to_string().contains("'chrZ' has alignments but no @SQ line"), "{err}");
    }

    #[test]
    fn test_resolve_bam_ref_id_without_header_uses_reflist_index() {
        assert_eq!(resolve_bam_ref_id(None, "chr1", 7).unwrap(), 7);
    }

    #[test]
    fn test_assign_bam_ref_ids() {
        let mut items = vec![wi(0), wi(1), wi(0)];
        for item in &mut items {
            item.ref_name = format!("chr{}", item.ref_idx);
        }
        let map = HashMap::from([("chr1".to_string(), 0), ("chr0".to_string(), 1)]);
        assign_bam_ref_ids(&mut items, Some(&map)).unwrap();
        let ids: Vec<_> = items.iter().map(|w| w.bam_ref_id).collect();
        assert_eq!(ids, [1, 0, 1]);
    }

    /// Boundaries for references `chr0..chr{n-1}`, 10 rows each.
    fn boundaries_for(n: u32) -> Vec<RefBoundary> {
        (0..n).map(|i| boundary(i, i64::from(i) * 10 + 1, i64::from(i) * 10 + 10)).collect()
    }

    /// A `ref_idx_of` for region tests: resolves `NC_<i>` to index `i`.
    fn seqid_to_idx(name: &str) -> Result<u32> {
        name.strip_prefix("NC_")
            .and_then(|i| i.parse().ok())
            .with_context(|| format!("reference not found: {name}"))
    }

    #[test]
    fn test_assign_bam_ref_ids_needs_only_output_references() {
        // The table has chr0..chr4 but only chr4 is output (one --aligned-region);
        // a header listing only chr4 must suffice.
        let boundaries = boundaries_for(5);
        let regions = ["chr4".to_string()];
        let mut items =
            collect_row_range_region_work_items(&boundaries, &regions, seqid_to_idx, 1).unwrap();
        let map = HashMap::from([("chr4".to_string(), 0)]);
        assign_bam_ref_ids(&mut items, Some(&map)).unwrap();
        assert_eq!(items.iter().map(|w| w.bam_ref_id).collect::<Vec<_>>(), [0]);
        // Whereas resolving every reference in the table fails.
        let mut all = collect_row_range_work_items(&boundaries, 1);
        assert!(assign_bam_ref_ids(&mut all, Some(&map)).is_err());
    }

    #[test]
    fn test_region_work_items_cover_every_run_of_a_split_reference() {
        let boundaries = vec![
            boundary(0, 1, 10),
            boundary(0, 31, 40),
            boundary(1, 11, 20),
            boundary(1, 41, 50),
            boundary(2, 21, 30),
        ];
        let regions = ["chr1:5-100".to_string()];
        let items =
            collect_row_range_region_work_items(&boundaries, &regions, seqid_to_idx, 1).unwrap();
        let got: Vec<_> =
            items.iter().map(|w| (w.ref_idx, w.start_row, w.end_row, w.region_filter)).collect();
        assert_eq!(got, [(1, 11, 20, Some((4, 100))), (1, 41, 50, Some((4, 100)))]);
    }

    #[test]
    fn test_region_work_items_resolve_alternate_name() {
        let boundaries = vec![boundary(0, 1, 10), boundary(1, 11, 20), boundary(1, 21, 30)];
        let regions = ["NC_1".to_string()];
        let items =
            collect_row_range_region_work_items(&boundaries, &regions, seqid_to_idx, 1).unwrap();
        let got: Vec<_> = items.iter().map(|w| (w.start_row, w.end_row)).collect();
        assert_eq!(got, [(11, 20), (21, 30)]);
    }

    #[test]
    fn test_region_work_items_errors() {
        let boundaries = boundaries_for(2);
        // Known reference without alignments.
        let err =
            collect_row_range_region_work_items(&boundaries, &["NC_7".into()], seqid_to_idx, 1)
                .unwrap_err();
        assert!(err.to_string().contains("no alignments found for reference: NC_7"), "{err}");
        // Unknown reference.
        let err =
            collect_row_range_region_work_items(&boundaries, &["chrZ".into()], seqid_to_idx, 1)
                .unwrap_err();
        assert!(err.to_string().contains("reference not found: chrZ"), "{err}");
    }

    /// A work item for `chr{ref_idx}` covering rows `start..=end`.
    fn item(ref_idx: u32, start_row: i64, end_row: i64) -> RowRangeWorkItem {
        let mut w = wi(ref_idx);
        w.ref_name = format!("chr{ref_idx}");
        w.start_row = start_row;
        w.end_row = end_row;
        w
    }

    /// Header `@SQ` ids for `chr0..chr{n-1}` in that order.
    fn header_ids(n: u32) -> HashMap<String, i32> {
        (0..n).map(|i| (format!("chr{i}"), i as i32)).collect()
    }

    #[test]
    fn test_work_items_are_coordinate_sorted() {
        let ids = header_ids(3);
        let none = HashSet::new();
        let sorted = [item(0, 1, 10), item(0, 11, 20), item(2, 21, 30)];
        assert!(work_items_are_coordinate_sorted(&sorted, &none, &ids));
        assert!(work_items_are_coordinate_sorted(&[], &none, &ids));
        // Rows revisited within a reference (overlapping or repeated regions),
        // including by a single row at the boundary.
        for repeated in [[item(0, 1, 10), item(0, 5, 10)], [item(0, 1, 10), item(0, 10, 12)]] {
            assert!(!work_items_are_coordinate_sorted(&repeated, &none, &ids));
        }
        // References go against header order, or a reference is revisited.
        let header_swapped = HashMap::from([("chr0".to_string(), 1), ("chr1".to_string(), 0)]);
        let items = [item(0, 1, 10), item(1, 11, 20)];
        assert!(!work_items_are_coordinate_sorted(&items, &none, &header_swapped));
        let revisit = [item(0, 1, 10), item(1, 11, 20), item(0, 21, 30)];
        assert!(!work_items_are_coordinate_sorted(&revisit, &none, &ids));
        // A reference missing from the header makes the order unknown.
        assert!(!work_items_are_coordinate_sorted(&[item(5, 1, 10)], &none, &ids));
    }

    #[test]
    fn test_split_references_break_order_only_when_output() {
        let ids = header_ids(3);
        let split = HashSet::from(["chr1".to_string()]);
        // chr1's two runs are grouped: rows increase but positions restart.
        let with_split = [item(0, 1, 10), item(1, 11, 20), item(1, 31, 40), item(2, 21, 30)];
        assert!(!work_items_are_coordinate_sorted(&with_split, &split, &ids));
        // An --aligned-region excluding chr1 is sorted despite chr1 being split.
        let without_split = [item(0, 1, 10), item(2, 21, 30)];
        assert!(work_items_are_coordinate_sorted(&without_split, &split, &ids));
    }

    #[test]
    fn test_planned_work_items_of_reordered_batches_are_sorted() {
        // Two batches covering disjoint references, stored out of order. The
        // planner groups them into reference order, which is coordinate order
        // even though rows no longer increase across references.
        let runs = vec![boundary(2, 1, 10), boundary(0, 11, 20), boundary(1, 21, 30)];
        let (grouped, split) = group_runs_by_reference(runs);
        let items = collect_row_range_work_items(&grouped, 1);
        let split: HashSet<String> = split.into_iter().collect();
        assert!(work_items_are_coordinate_sorted(&items, &split, &header_ids(3)));
    }

    #[test]
    fn test_planned_region_work_items_classification() {
        let boundaries = boundaries_for(3);
        let classify = |regions: &[&str]| {
            let regions: Vec<String> = regions.iter().map(ToString::to_string).collect();
            let items = collect_row_range_region_work_items(&boundaries, &regions, seqid_to_idx, 1)
                .unwrap();
            work_items_are_coordinate_sorted(&items, &HashSet::new(), &header_ids(3))
        };
        assert!(classify(&["chr0", "chr2"]));
        assert!(!classify(&["chr2", "chr0"]), "regions out of @SQ order");
        assert!(!classify(&["chr1:1-5", "chr1:3-9"]), "overlapping regions");
    }

    #[test]
    fn test_tables_are_coordinate_sorted() {
        let ids = header_ids(2);
        let none = HashSet::new();
        let items = [item(0, 1, 10), item(1, 11, 20)];
        assert!(tables_are_coordinate_sorted(std::iter::empty(), &ids));
        assert!(tables_are_coordinate_sorted([(&none, &items[..])].into_iter(), &ids));
        // Secondary alignments follow all primary alignments.
        let two = [(&none, &items[..]), (&none, &items[..])];
        assert!(!tables_are_coordinate_sorted(two.into_iter(), &ids));
    }

    #[test]
    fn test_split_layout_warning() {
        let layout = |split: &[&str]| RefLayout {
            boundaries: Vec::new(),
            num_sorted_segments: 3,
            split_references: split.iter().map(ToString::to_string).collect(),
        };
        assert_eq!(
            split_layout_warning("T", &layout(&["1", "2"])),
            "T is stored as 3 separately sorted row ranges; alignments for 2 reference(s) \
             (1, 2) are grouped by reference but are not coordinate-sorted within those \
             references"
        );
    }

    #[test]
    fn test_abbreviate_names() {
        let names: Vec<String> = ["a", "b", "c"].iter().map(ToString::to_string).collect();
        assert_eq!(abbreviate_names(&names, 5), "a, b, c");
        assert_eq!(abbreviate_names(&names, 3), "a, b, c");
        assert_eq!(abbreviate_names(&names, 2), "a, b, and 1 more");
        assert_eq!(abbreviate_names(&[], 2), "");
    }

    #[test]
    fn test_collect_row_range_split_reference_keeps_run_labels() {
        // Two runs of chr0 (grouped by reference) then a run of chr1.
        let boundaries = vec![boundary(0, 1, 100), boundary(0, 201, 300), boundary(1, 101, 200)];
        let items = collect_row_range_work_items(&boundaries, 1);
        let got: Vec<_> =
            items.iter().map(|w| (w.ref_name.as_str(), w.start_row, w.end_row)).collect();
        assert_eq!(got, [("chr0", 1, 100), ("chr0", 201, 300), ("chr1", 101, 200)]);
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
