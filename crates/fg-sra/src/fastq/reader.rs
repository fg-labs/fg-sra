//! Reading spots from an archive's reads table.

use std::ops::RangeInclusive;

use anyhow::{Context, Result, anyhow};
use fg_sra_vdb::cursor::{BlobColumn, VCursor};
use fg_sra_vdb::database::VTable;
use fg_sra_vdb::error::VdbError;

use super::prefetch::PrefetchedReads;
use super::spot::Spot;
use crate::archive::marks_rna;
use crate::refstore::{
    CHARSET_4NA, CHARSET_4NA_RNA, ReferenceRows, ReferenceStore, ref_window_len,
};
use crate::restore_read::{restore_read, restore_spot};

/// Bytes of VDB blob cache for each alignment cursor. Each batch's alignments are read in
/// ascending order, so little is revisited, but without a cache 2-9% more bytes are read.
const ALIGNMENT_CURSOR_CACHE_BYTES: usize = 8 * 1024 * 1024;

/// Something spots can be read from by id; each worker thread has its own.
pub trait SpotSource {
    /// Prepare to read the spots `ids`, which are read next, in order.
    fn start_batch(&mut self, _ids: RangeInclusive<i64>) -> Result<()> {
        Ok(())
    }

    /// Read spot `id`. The spot borrows from the source until the next read.
    fn read(&mut self, id: i64) -> Result<Spot<'_>>;
}

/// Which optional columns to read.
#[derive(Debug, Clone, Copy)]
pub struct ReadColumns {
    /// `QUALITY`, for FASTQ output.
    pub qualities: bool,
    /// `NAME`, for `$sn`.
    pub names: bool,
    /// `SPOT_GROUP`, for `$sg`.
    pub spot_groups: bool,
}

/// What rebuilding an aligned archive's reads needs, shared by every reader: its references,
/// preloaded, and which reference each `REFERENCE` row belongs to.
#[derive(Clone, Copy)]
pub struct References<'a> {
    /// Every reference's bases.
    pub store: &'a ReferenceStore,
    /// The reference and position of each `REFERENCE` row.
    pub rows: &'a ReferenceRows,
}

/// Reads spots through a VDB cursor on the archive's reads table, copying each spot's cells
/// into buffers reused from spot to spot.
///
/// For an aligned archive the reads table holds only the unaligned reads' bases
/// (`CMP_READ`); the aligned reads are rebuilt from their alignments and the preloaded
/// references, never through the virtual `READ` column, whose reference lookups are not
/// safe on concurrent cursors.
///
/// Each reader has its own cursors, so readers on one table can be used on different
/// threads. A cursor holds its own reference to its table, but the table must stay open for
/// as long as the cursor's data is being read.
pub struct VdbSpotReader<'a> {
    cursor: VCursor,
    columns: ColumnIds,
    buffers: SpotBuffers,
    aligned: Option<AlignedReads<'a>>,
}

impl<'a> VdbSpotReader<'a> {
    /// A reader of a table that stores bases directly, with its own cursor on `table`.
    pub fn new(table: &VTable, columns: ReadColumns) -> Result<Self> {
        Self::open(table, columns, None)
    }

    /// A reader of an aligned archive's reads table `table`, rebuilding aligned reads from
    /// `alignments` (its `PRIMARY_ALIGNMENT` table) and `references`, or taking them from
    /// `prefetched`.
    pub fn new_aligned(
        table: &VTable,
        alignments: &VTable,
        columns: ReadColumns,
        references: References<'a>,
        prefetched: Option<&'a PrefetchedReads>,
    ) -> Result<Self> {
        let aligned = AlignedReads {
            alignments: AlignmentReader::new(alignments, references)?,
            prefetched,
            batch: BatchReads::default(),
            cmp_read: Vec::new(),
            aligned: Vec::new(),
        };
        Self::open(table, columns, Some(aligned))
    }

    /// A reader with its own cursor on `table`, adding the columns needed: `CMP_READ` when
    /// `aligned` is given, otherwise `READ`.
    fn open(
        table: &VTable,
        columns: ReadColumns,
        aligned: Option<AlignedReads<'a>>,
    ) -> Result<Self> {
        // Uncached, as blob reads (`BlobColumn`) require; spots are read in row order, so the
        // one blob each column holds serves its whole run of rows.
        let cursor = table.create_cursor_read()?;
        let add = |name: &str| {
            cursor
                .add_column(name)
                .map(BlobColumn::new)
                .with_context(|| format!("failed to read column {name}"))
        };
        let optional = |wanted: bool, name: &str| wanted.then(|| add(name)).transpose();
        // Bases are read as 4na codes and mapped to text here, which is cheaper than
        // libncbi-vdb's own per-base map; a table that can't give 4na is read as text. The map
        // matches the text the schema serves: an unaligned table marked as RNA gives `U` for
        // `T`, while an aligned archive's `READ` never does (though its text `CMP_READ` would).
        let (name, charset) = if aligned.is_some() {
            ("CMP_READ", CHARSET_4NA)
        } else if marks_rna(table) {
            ("READ", CHARSET_4NA_RNA)
        } else {
            ("READ", CHARSET_4NA)
        };
        let (bases, charset) = match cursor.add_column(&format!("(INSDC:4na:bin){name}")) {
            Ok(column) => (BlobColumn::new(column), Some(charset)),
            Err(_) => (add(&format!("(INSDC:dna:text){name}"))?, None),
        };
        let columns = ColumnIds {
            bases,
            charset,
            qualities: optional(columns.qualities, "(INSDC:quality:phred)QUALITY")?,
            read_starts: add("(INSDC:coord:zero)READ_START")?,
            read_lens: add("(INSDC:coord:len)READ_LEN")?,
            read_types: add("(INSDC:SRA:xread_type)READ_TYPE")?,
            read_filters: add("(INSDC:SRA:read_filter)READ_FILTER")?,
            names: optional(columns.names, "(ascii)NAME")?,
            spot_groups: optional(columns.spot_groups, "(ascii)SPOT_GROUP")?,
            align_ids: optional(aligned.is_some(), "(I64)PRIMARY_ALIGNMENT_ID")?,
        };
        cursor.open().context("failed to open a cursor on the reads table")?;
        Ok(Self { cursor, columns, buffers: SpotBuffers::default(), aligned })
    }

    /// Copy spot `id`'s cells into the buffers; optional columns not read are left empty.
    /// For an aligned archive, the bases read are `CMP_READ`, still to be restored.
    fn read_cells(&mut self, id: i64) -> Result<(), VdbError> {
        let (cursor, columns, buffers) = (&self.cursor, &mut self.columns, &mut self.buffers);
        let read_optional = |column: &mut Option<BlobColumn>, buffer: &mut Vec<u8>| {
            if let Some(column) = column {
                column.read_u8_slice_into(cursor, id, buffer)
            } else {
                buffer.clear();
                Ok(())
            }
        };
        columns.bases.read_u8_slice_into(cursor, id, &mut buffers.bases)?;
        if let Some(charset) = columns.charset {
            for base in &mut buffers.bases {
                *base = charset[usize::from(*base & 0x0F)];
            }
        }
        read_optional(&mut columns.qualities, &mut buffers.qualities)?;
        columns.read_starts.read_i32_slice_into(cursor, id, &mut buffers.read_starts)?;
        columns.read_lens.read_u32_slice_into(cursor, id, &mut buffers.read_lens)?;
        columns.read_types.read_u8_slice_into(cursor, id, &mut buffers.read_types)?;
        columns.read_filters.read_u8_slice_into(cursor, id, &mut buffers.read_filters)?;
        read_optional(&mut columns.names, &mut buffers.name)?;
        read_optional(&mut columns.spot_groups, &mut buffers.spot_group)?;
        truncate_at_nul(&mut buffers.name);
        truncate_at_nul(&mut buffers.spot_group);
        if let Some(column) = &mut columns.align_ids {
            column.read_i64_slice_into(cursor, id, &mut buffers.align_ids)?;
        }
        Ok(())
    }
}

impl SpotSource for VdbSpotReader<'_> {
    fn start_batch(&mut self, ids: RangeInclusive<i64>) -> Result<()> {
        match (&mut self.aligned, &mut self.columns.align_ids) {
            (Some(aligned), Some(column)) => aligned.restore_batch(&self.cursor, column, ids),
            _ => Ok(()),
        }
    }

    fn read(&mut self, id: i64) -> Result<Spot<'_>> {
        self.read_cells(id).with_context(|| format!("failed to read spot {id}"))?;
        if let Some(aligned) = &mut self.aligned {
            aligned
                .restore(&mut self.buffers)
                .with_context(|| format!("failed to restore spot {id}"))?;
        }
        let buffers = &self.buffers;
        Ok(Spot {
            id,
            name: &buffers.name,
            group: &buffers.spot_group,
            bases: &buffers.bases,
            qualities: &buffers.qualities,
            read_starts: &buffers.read_starts,
            read_lens: &buffers.read_lens,
            read_types: &buffers.read_types,
            read_filters: &buffers.read_filters,
        })
    }
}

/// Rebuilds an aligned archive's aligned reads for a reader: from the prefetched reads or the
/// batch's when they hold one, and otherwise through its own alignment cursor.
struct AlignedReads<'a> {
    alignments: AlignmentReader<'a>,
    prefetched: Option<&'a PrefetchedReads>,
    /// The current batch's aligned reads that weren't prefetched.
    batch: BatchReads,
    /// The spot's `CMP_READ`.
    cmp_read: Vec<u8>,
    /// The aligned read being rebuilt, in reference orientation.
    aligned: Vec<u8>,
}

impl AlignedReads<'_> {
    /// Rebuild the aligned reads of spots `ids` that weren't prefetched, in alignment order,
    /// reading the spots' `PRIMARY_ALIGNMENT_ID` through `column` of `cursor`.
    ///
    /// Some archives number their spots far out of alignment order, so a batch's spots, taken
    /// in turn, revisit hundreds of blobs of the alignment columns: more than the cursor's
    /// cache holds. In alignment order each blob is decoded once per batch.
    fn restore_batch(
        &mut self,
        cursor: &VCursor,
        column: &mut BlobColumn,
        ids: RangeInclusive<i64>,
    ) -> Result<()> {
        let Self { alignments, prefetched, batch, .. } = self;
        batch.clear();
        for spot in ids {
            column
                .read_i64_slice_into(cursor, spot, &mut batch.cell)
                .with_context(|| format!("failed to read spot {spot}'s alignment ids"))?;
            let not_prefetched = |id: &i64| prefetched.is_none_or(|p| p.get(*id).is_none());
            batch
                .ids
                .extend(batch.cell.iter().copied().filter(|&id| id > 0).filter(not_prefetched));
        }
        batch.ids.sort_unstable();
        batch.ids.dedup();
        let BatchReads { ids, starts, bases, read, .. } = batch;
        for &id in ids.iter() {
            alignments.restore(id, read)?;
            starts.push(bases.len());
            bases.extend_from_slice(read);
        }
        starts.push(bases.len());
        Ok(())
    }

    /// Replace `spot.bases`, which holds the spot's `CMP_READ`, with all of its bases.
    fn restore(&mut self, spot: &mut SpotBuffers) -> Result<()> {
        let Self { alignments, prefetched, batch, cmp_read, aligned } = self;
        std::mem::swap(&mut spot.bases, cmp_read);
        let aligned_read = |align_id: i64, out: &mut Vec<u8>| -> Result<()> {
            if let Some(bases) =
                prefetched.and_then(|p| p.get(align_id)).or_else(|| batch.get(align_id))
            {
                out.clear();
                out.extend_from_slice(bases);
                return Ok(());
            }
            alignments.restore(align_id, out)
        };
        restore_spot(
            cmp_read,
            &spot.align_ids,
            &spot.read_lens,
            &spot.read_types,
            aligned_read,
            aligned,
            &mut spot.bases,
        )
    }
}

/// Aligned reads rebuilt for one batch, looked up by alignment id; reused from batch to batch.
#[derive(Default)]
struct BatchReads {
    /// Alignment ids, ascending.
    ids: Vec<i64>,
    /// `bases[starts[i]..starts[i + 1]]` is alignment `ids[i]`'s read, in reference orientation.
    starts: Vec<usize>,
    bases: Vec<u8>,
    /// Scratch for a spot's `PRIMARY_ALIGNMENT_ID` cell.
    cell: Vec<i64>,
    /// Scratch for the read being rebuilt.
    read: Vec<u8>,
}

impl BatchReads {
    /// Forget the previous batch's reads.
    fn clear(&mut self) {
        self.ids.clear();
        self.starts.clear();
        self.bases.clear();
    }

    /// Alignment `id`'s read, in reference orientation, if it was rebuilt for the batch.
    fn get(&self, id: i64) -> Option<&[u8]> {
        let index = self.ids.binary_search(&id).ok()?;
        Some(&self.bases[self.starts[index]..self.starts[index + 1]])
    }
}

/// A cursor on an aligned archive's alignments, and what rebuilding one alignment's read
/// needs.
pub(super) struct AlignmentReader<'a> {
    cursor: VCursor,
    columns: AlignmentColumnIds,
    references: References<'a>,
    cells: AlignmentCells,
    /// Scratch for a reference window that wraps around a circular reference.
    window: Vec<u8>,
}

impl<'a> AlignmentReader<'a> {
    /// A reader with its own cursor on `alignments`, the `PRIMARY_ALIGNMENT` table.
    pub(super) fn new(alignments: &VTable, references: References<'a>) -> Result<Self> {
        let cursor = alignments.create_cached_cursor_read(ALIGNMENT_CURSOR_CACHE_BYTES)?;
        let add = |name: &str| {
            cursor
                .add_column(name)
                .with_context(|| format!("failed to read alignment column {name}"))
        };
        let columns = AlignmentColumnIds {
            has_mismatch: add("(bool)HAS_MISMATCH")?,
            mismatch: add("(INSDC:dna:text)MISMATCH")?,
            has_ref_offset: add("(bool)HAS_REF_OFFSET")?,
            ref_offset: add("(I32)REF_OFFSET")?,
            ref_start: add("(INSDC:coord:zero)REF_START")?,
            ref_id: add("(I64)REF_ID")?,
        };
        cursor.open().context("failed to open a cursor on the alignments")?;
        Ok(Self {
            cursor,
            columns,
            references,
            cells: AlignmentCells::default(),
            window: Vec::new(),
        })
    }

    /// Rebuild alignment `align_id`'s read, in reference orientation, into `out`.
    pub(super) fn restore(&mut self, align_id: i64, out: &mut Vec<u8>) -> Result<()> {
        let Self { cursor, columns, references, cells, window } = self;
        cells
            .read(cursor, columns, align_id)
            .with_context(|| format!("failed to read alignment {align_id}"))?;
        let (ref_idx, ref_pos) =
            references.rows.locate(cells.ref_row, cells.ref_start).ok_or_else(|| {
                anyhow!("alignment {align_id}: no reference holds row {}", cells.ref_row)
            })?;
        let ref_len = ref_window_len(cells.has_ref_offset.len(), &cells.ref_offset);
        let reference = references
            .store
            .window(ref_idx, ref_pos, ref_len, window)
            .with_context(|| format!("alignment {align_id}: reference window"))?;
        restore_read(
            reference,
            &cells.has_mismatch,
            &cells.mismatch,
            &cells.has_ref_offset,
            &cells.ref_offset,
            &[],
            out,
        )
        .with_context(|| format!("alignment {align_id}: rebuilding its read"))
    }
}

/// The reads table's columns, each holding its current blob; `None` for optional columns not
/// read.
struct ColumnIds {
    bases: BlobColumn,
    /// The text of each 4na code, when `bases` gives 4na codes rather than text.
    charset: Option<&'static [u8; 16]>,
    qualities: Option<BlobColumn>,
    read_starts: BlobColumn,
    read_lens: BlobColumn,
    read_types: BlobColumn,
    read_filters: BlobColumn,
    names: Option<BlobColumn>,
    spot_groups: Option<BlobColumn>,
    /// `PRIMARY_ALIGNMENT_ID`, for an aligned archive.
    align_ids: Option<BlobColumn>,
}

/// One spot's cells, reused from spot to spot.
#[derive(Default)]
struct SpotBuffers {
    bases: Vec<u8>,
    qualities: Vec<u8>,
    read_starts: Vec<i32>,
    read_lens: Vec<u32>,
    read_types: Vec<u8>,
    read_filters: Vec<u8>,
    name: Vec<u8>,
    spot_group: Vec<u8>,
    align_ids: Vec<i64>,
}

/// Alignment-cursor column indices.
struct AlignmentColumnIds {
    has_mismatch: u32,
    mismatch: u32,
    has_ref_offset: u32,
    ref_offset: u32,
    ref_start: u32,
    ref_id: u32,
}

/// One alignment's stored cells.
#[derive(Default)]
struct AlignmentCells {
    has_mismatch: Vec<u8>,
    mismatch: Vec<u8>,
    has_ref_offset: Vec<u8>,
    ref_offset: Vec<i32>,
    /// `REF_ID`: the `REFERENCE` table row the alignment starts in.
    ref_row: i64,
    /// `REF_START`: 0-based start within that row. The whole-reference position (`REF_POS`)
    /// is computed from it rather than read, as the schema derives `REF_POS` one row at a time.
    ref_start: i32,
}

impl AlignmentCells {
    /// Read alignment `align_id`'s cells through `cursor`.
    fn read(
        &mut self,
        cursor: &VCursor,
        columns: &AlignmentColumnIds,
        align_id: i64,
    ) -> Result<(), VdbError> {
        cursor.read_u8_slice_into(align_id, columns.has_mismatch, &mut self.has_mismatch)?;
        cursor.read_u8_slice_into(align_id, columns.mismatch, &mut self.mismatch)?;
        cursor.read_u8_slice_into(align_id, columns.has_ref_offset, &mut self.has_ref_offset)?;
        cursor.read_i32_slice_into(align_id, columns.ref_offset, &mut self.ref_offset)?;
        self.ref_row = cursor.read_i64(align_id, columns.ref_id)?;
        self.ref_start = cursor.read_coord_zero(align_id, columns.ref_start)?;
        Ok(())
    }
}

/// Cut a text cell at its first NUL. Some loaders store `NAME` and `SPOT_GROUP` NUL-padded to a
/// fixed width, and sra-tools reads them as C strings.
fn truncate_at_nul(text: &mut Vec<u8>) {
    if let Some(end) = text.iter().position(|&b| b == 0) {
        text.truncate(end);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rna_charset_differs_from_dna_only_at_t() {
        let differing: Vec<usize> =
            (0..16).filter(|&i| CHARSET_4NA[i] != CHARSET_4NA_RNA[i]).collect();
        assert_eq!(differing, [8]);
        assert_eq!((CHARSET_4NA[8], CHARSET_4NA_RNA[8]), (b'T', b'U'));
    }

    #[test]
    fn nul_padding_is_cut_from_text() {
        let mut text = b"tagged_866\0\0\0".to_vec();
        truncate_at_nul(&mut text);
        assert_eq!(text, b"tagged_866");
    }

    #[test]
    fn text_without_nul_is_unchanged() {
        let mut text = b"B00EBABXX:4:1:8249:1985".to_vec();
        truncate_at_nul(&mut text);
        assert_eq!(text, b"B00EBABXX:4:1:8249:1985");
    }

    #[test]
    fn all_nul_text_becomes_empty() {
        let mut text = vec![0; 8];
        truncate_at_nul(&mut text);
        assert!(text.is_empty());
    }
}
