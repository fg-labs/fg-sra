//! Safe wrapper for `VCursor` with typed column reads.
//!
//! Provides methods for adding columns, opening the cursor, and reading
//! typed cell data via `VCursorCellDataDirect`, or through each column's
//! blobs with [`BlobColumn`].

use std::ops::RangeInclusive;
use std::ptr;
use std::slice;

use crate::error::{LITERAL_FORMAT, VdbError, check_rc, to_cstring};
use crate::retry::retry_on_network_error;

/// Sentinel value for columns that were never added or are not available.
pub const INVALID_COLUMN: u32 = 0xFFFF_FFFF;

/// Reject a cell whose element width differs from what a typed read expects.
///
/// The typed readers interpret the raw cell as `T` using `row_len` as the
/// element count. If the cell's `elem_bits` is narrower than `T`, the derived
/// byte length exceeds the cell, so `slice::from_raw_parts` would read out of
/// bounds. Validate in every build (not just debug) before casting.
fn check_elem_bits(actual: u32, expected: u32) -> Result<(), VdbError> {
    if actual == expected { Ok(()) } else { Err(VdbError::ElemBitsMismatch { expected, actual }) }
}

/// Safe wrapper around the VDB `VCursor` opaque type.
///
/// A cursor is opened on a table and reads column data row-by-row.
/// Columns must be added before opening, and data is read via
/// `VCursorCellDataDirect` (zero-copy access into VDB page cache).
pub struct VCursor {
    ptr: *const fg_sra_vdb_sys::VCursor,
    /// Whether the cursor keeps a cache of recently read blobs.
    cached: bool,
}

unsafe impl Send for VCursor {}

impl VCursor {
    pub(crate) fn from_raw(ptr: *const fg_sra_vdb_sys::VCursor, cached: bool) -> Self {
        Self { ptr, cached }
    }

    /// Add a column to an unopened cursor.
    ///
    /// Returns the column index used for subsequent reads.
    /// Column names should include the type cast, e.g. `"(I64)SEQ_SPOT_ID"`.
    pub fn add_column(&self, name: &str) -> Result<u32, VdbError> {
        let c_name = to_cstring(name)?;
        let mut idx: u32 = 0;
        // Safety: VCursorAddColumn is printf-style variadic; the name is the single argument
        // to a literal "%s" format.
        let rc = unsafe {
            fg_sra_vdb_sys::VCursorAddColumn(
                self.ptr,
                &raw mut idx,
                LITERAL_FORMAT.as_ptr(),
                c_name.as_ptr(),
            )
        };
        check_rc(rc)?;
        Ok(idx)
    }

    /// Add a column, returning `None` if the column does not exist.
    ///
    /// This is equivalent to `add_opt_column` in the C implementation.
    #[must_use]
    pub fn add_column_optional(&self, name: &str) -> Option<u32> {
        self.add_column(name).ok()
    }

    /// Open the cursor after adding columns.
    pub fn open(&self) -> Result<(), VdbError> {
        let rc = unsafe { fg_sra_vdb_sys::VCursorOpen(self.ptr) };
        check_rc(rc)
    }

    /// Get the row ID range for a column (or all columns if `col_idx` is 0).
    ///
    /// Returns `(first_row_id, row_count)`.
    pub fn id_range(&self, col_idx: u32) -> Result<(i64, u64), VdbError> {
        let mut first: i64 = 0;
        let mut count: u64 = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::VCursorIdRange(self.ptr, col_idx, &raw mut first, &raw mut count)
        };
        check_rc(rc)?;
        Ok((first, count))
    }

    /// Read a single scalar value from a cell.
    ///
    /// Returns the default value (zero) if the cell is empty.
    fn read_scalar<T: Copy + Default>(
        &self,
        row_id: i64,
        col_idx: u32,
        expected_bits: u32,
    ) -> Result<T, VdbError> {
        self.cell_data_direct(row_id, col_idx)?.scalar(expected_bits)
    }

    /// Read a slice of values from a multi-element cell.
    ///
    /// Returns an empty `Vec` if the cell is empty.
    fn read_slice<T: Copy>(
        &self,
        row_id: i64,
        col_idx: u32,
        expected_bits: u32,
    ) -> Result<Vec<T>, VdbError> {
        let data = self.cell_data_direct(row_id, col_idx)?;
        if data.row_len == 0 {
            return Ok(Vec::new());
        }
        check_elem_bits(data.elem_bits, expected_bits)?;
        let values = unsafe { slice::from_raw_parts(data.base.cast::<T>(), data.row_len as usize) };
        Ok(values.to_vec())
    }

    /// Like [`read_slice`](Self::read_slice), but into an existing `Vec`,
    /// reusing its allocation. Clears `buf` and appends the cell data.
    fn read_slice_into<T: Copy>(
        &self,
        row_id: i64,
        col_idx: u32,
        expected_bits: u32,
        buf: &mut Vec<T>,
    ) -> Result<(), VdbError> {
        self.cell_data_direct(row_id, col_idx)?.copy_into(expected_bits, buf)
    }

    /// Read a single `i64` value from a cell.
    pub fn read_i64(&self, row_id: i64, col_idx: u32) -> Result<i64, VdbError> {
        self.read_scalar(row_id, col_idx, 64)
    }

    /// Read a single `i32` value from a cell.
    pub fn read_i32(&self, row_id: i64, col_idx: u32) -> Result<i32, VdbError> {
        self.read_scalar(row_id, col_idx, 32)
    }

    /// Read a single `u32` value from a cell.
    pub fn read_u32(&self, row_id: i64, col_idx: u32) -> Result<u32, VdbError> {
        self.read_scalar(row_id, col_idx, 32)
    }

    /// Read a single `u8` value from a cell.
    pub fn read_u8(&self, row_id: i64, col_idx: u32) -> Result<u8, VdbError> {
        self.read_scalar(row_id, col_idx, 8)
    }

    /// Read a single `bool` value from a cell.
    pub fn read_bool(&self, row_id: i64, col_idx: u32) -> Result<bool, VdbError> {
        Ok(self.read_scalar::<u8>(row_id, col_idx, 8)? != 0)
    }

    /// Read a string (ASCII) cell.
    ///
    /// The underlying data points into VDB's page cache; this method copies
    /// it into a new `String`.
    pub fn read_str(&self, row_id: i64, col_idx: u32) -> Result<String, VdbError> {
        let data = self.cell_data_direct(row_id, col_idx)?;
        if data.row_len == 0 {
            return Ok(String::new());
        }
        check_elem_bits(data.elem_bits, 8)?;
        let bytes = unsafe { slice::from_raw_parts(data.base.cast::<u8>(), data.row_len as usize) };
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }

    /// Read a string cell into an existing `String`, reusing its allocation.
    ///
    /// Clears `buf` and appends the cell data. The existing heap allocation
    /// is reused, avoiding a new allocation per call once `buf` has grown to
    /// the typical cell size.
    pub fn read_str_into(
        &self,
        row_id: i64,
        col_idx: u32,
        buf: &mut String,
    ) -> Result<(), VdbError> {
        buf.clear();
        let data = self.cell_data_direct(row_id, col_idx)?;
        if data.row_len == 0 {
            return Ok(());
        }
        check_elem_bits(data.elem_bits, 8)?;
        let bytes = unsafe { slice::from_raw_parts(data.base.cast::<u8>(), data.row_len as usize) };
        // VDB ASCII columns are always valid UTF-8; from_utf8_lossy borrows
        // without allocation in the common (valid) case.
        buf.push_str(&String::from_utf8_lossy(bytes));
        Ok(())
    }

    /// Read a slice of `u8` values (e.g., quality scores, read data).
    pub fn read_u8_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<u8>, VdbError> {
        self.read_slice(row_id, col_idx, 8)
    }

    /// Read a slice of `i64` values.
    pub fn read_i64_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<i64>, VdbError> {
        self.read_slice(row_id, col_idx, 64)
    }

    /// Read a slice of `u32` values (e.g., `INSDC_coord_len`).
    pub fn read_u32_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<u32>, VdbError> {
        self.read_slice(row_id, col_idx, 32)
    }

    /// Read a slice of `i32` values (e.g., `INSDC_coord_zero`).
    pub fn read_i32_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<i32>, VdbError> {
        self.read_slice(row_id, col_idx, 32)
    }

    /// Read a slice of `u8` values into an existing `Vec`, reusing its allocation.
    pub fn read_u8_slice_into(
        &self,
        row_id: i64,
        col_idx: u32,
        buf: &mut Vec<u8>,
    ) -> Result<(), VdbError> {
        self.read_slice_into(row_id, col_idx, 8, buf)
    }

    /// Read a slice of `i64` values into an existing `Vec`, reusing its allocation.
    pub fn read_i64_slice_into(
        &self,
        row_id: i64,
        col_idx: u32,
        buf: &mut Vec<i64>,
    ) -> Result<(), VdbError> {
        self.read_slice_into(row_id, col_idx, 64, buf)
    }

    /// Read a slice of `u32` values into an existing `Vec`, reusing its allocation.
    pub fn read_u32_slice_into(
        &self,
        row_id: i64,
        col_idx: u32,
        buf: &mut Vec<u32>,
    ) -> Result<(), VdbError> {
        self.read_slice_into(row_id, col_idx, 32, buf)
    }

    /// Read a slice of `i32` values into an existing `Vec`, reusing its allocation.
    pub fn read_i32_slice_into(
        &self,
        row_id: i64,
        col_idx: u32,
        buf: &mut Vec<i32>,
    ) -> Result<(), VdbError> {
        self.read_slice_into(row_id, col_idx, 32, buf)
    }

    /// Read an `INSDC_coord_zero` value (0-based coordinate, i32).
    pub fn read_coord_zero(&self, row_id: i64, col_idx: u32) -> Result<i32, VdbError> {
        self.read_i32(row_id, col_idx)
    }

    /// Read an `INSDC_coord_len` value (length, u32).
    pub fn read_coord_len(&self, row_id: i64, col_idx: u32) -> Result<u32, VdbError> {
        self.read_u32(row_id, col_idx)
    }

    /// Low-level: read raw cell data via `VCursorCellDataDirect`.
    ///
    /// Returns a `CellData` with a pointer into VDB's page cache.
    /// Retries on transient network errors with exponential backoff.
    fn cell_data_direct(&self, row_id: i64, col_idx: u32) -> Result<CellData, VdbError> {
        retry_on_network_error("cell_data_direct", || {
            let mut elem_bits: u32 = 0;
            let mut base: *const std::ffi::c_void = ptr::null();
            let mut boff: u32 = 0;
            let mut row_len: u32 = 0;
            let rc = unsafe {
                fg_sra_vdb_sys::VCursorCellDataDirect(
                    self.ptr,
                    row_id,
                    col_idx,
                    &raw mut elem_bits,
                    &raw mut base,
                    &raw mut boff,
                    &raw mut row_len,
                )
            };
            check_rc(rc)?;
            Ok(CellData { elem_bits, base, _boff: boff, row_len })
        })
    }

    /// The blob of column `col_idx` holding row `row_id`, or `None` for a row the column has no
    /// value for.
    ///
    /// Refused on a cursor with a blob cache: there `VCursorGetBlobDirect` hands out one
    /// reference too many whenever the blob was not already cached (`VCursorReadColumnDirectInt`
    /// passes on the blob it produced without releasing it, then `VTableCursorGetBlobDirect`
    /// adds its own), so the blob would never be freed.
    pub(crate) fn blob_direct(&self, row_id: i64, col_idx: u32) -> Result<Option<VBlob>, VdbError> {
        if self.cached {
            return Err(VdbError::BlobReadOnCachedCursor);
        }
        let mut blob: *const fg_sra_vdb_sys::VBlob = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::VCursorGetBlobDirect(self.ptr, &raw mut blob, row_id, col_idx)
        };
        check_rc(rc)?;
        if blob.is_null() {
            return Ok(None);
        }
        let mut first: i64 = 0;
        let mut count: u64 = 0;
        let rc = unsafe { fg_sra_vdb_sys::VBlobIdRange(blob, &raw mut first, &raw mut count) };
        match check_rc(rc).and_then(|()| blob_rows(row_id, first, count)) {
            Ok(rows) => Ok(Some(VBlob { ptr: blob, rows })),
            Err(e) => {
                unsafe { fg_sra_vdb_sys::VBlobRelease(blob) };
                Err(e)
            }
        }
    }

    /// Get the raw cursor pointer (for passing to `PlacementIterator` creation).
    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::VCursor {
        self.ptr
    }
}

impl Drop for VCursor {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::VCursorRelease(self.ptr) };
        }
    }
}

/// Reads one column's cells through its blobs, keeping the blob of the last row read so that
/// rows of the same blob are served from it directly rather than through the cursor's column
/// production chain, which `VCursorCellDataDirect` walks on every call.
///
/// Reads fail with [`VdbError::BlobReadOnCachedCursor`] on a cursor made with a blob cache.
pub struct BlobColumn {
    col_idx: u32,
    blob: Option<VBlob>,
}

impl BlobColumn {
    /// A reader of the column at `col_idx` of the cursor it is later read through.
    #[must_use]
    pub fn new(col_idx: u32) -> Self {
        Self { col_idx, blob: None }
    }

    /// Read a slice of `u8` values from row `row_id`'s cell into `buf`, replacing its contents.
    pub fn read_u8_slice_into(
        &mut self,
        cursor: &VCursor,
        row_id: i64,
        buf: &mut Vec<u8>,
    ) -> Result<(), VdbError> {
        self.cell(cursor, row_id)?.copy_into(8, buf)
    }

    /// Read a slice of `i32` values from row `row_id`'s cell into `buf`, replacing its contents.
    pub fn read_i32_slice_into(
        &mut self,
        cursor: &VCursor,
        row_id: i64,
        buf: &mut Vec<i32>,
    ) -> Result<(), VdbError> {
        self.cell(cursor, row_id)?.copy_into(32, buf)
    }

    /// Read a slice of `u32` values from row `row_id`'s cell into `buf`, replacing its contents.
    pub fn read_u32_slice_into(
        &mut self,
        cursor: &VCursor,
        row_id: i64,
        buf: &mut Vec<u32>,
    ) -> Result<(), VdbError> {
        self.cell(cursor, row_id)?.copy_into(32, buf)
    }

    /// Read a slice of `i64` values from row `row_id`'s cell into `buf`, replacing its contents.
    pub fn read_i64_slice_into(
        &mut self,
        cursor: &VCursor,
        row_id: i64,
        buf: &mut Vec<i64>,
    ) -> Result<(), VdbError> {
        self.cell(cursor, row_id)?.copy_into(64, buf)
    }

    /// Row `row_id`'s cell, from the blob held if it has the row and otherwise from the blob
    /// `cursor` gives for the row, which is then held. An empty cell stands for a row the
    /// column has no value for.
    fn cell(&mut self, cursor: &VCursor, row_id: i64) -> Result<CellData, VdbError> {
        if !self.blob.as_ref().is_some_and(|blob| blob.rows.contains(&row_id)) {
            self.blob = cursor.blob_direct(row_id, self.col_idx)?;
        }
        match &self.blob {
            Some(blob) => blob.cell_data(row_id),
            None => Ok(CellData::EMPTY),
        }
    }
}

/// The rows of a blob reported to start at row `first` and hold `count` rows, which must
/// include `row_id`, the row it was read for.
fn blob_rows(row_id: i64, first: i64, count: u64) -> Result<RangeInclusive<i64>, VdbError> {
    let last = i64::try_from(count).ok().and_then(|count| first.checked_add(count - 1));
    match last {
        Some(last) if (first..=last).contains(&row_id) => Ok(first..=last),
        _ => Err(VdbError::BlobMissesRow { row: row_id, first, count }),
    }
}

/// A reference to one of a column's blobs: the cells of a run of rows, decoded.
pub(crate) struct VBlob {
    ptr: *const fg_sra_vdb_sys::VBlob,
    /// The rows the blob holds.
    rows: RangeInclusive<i64>,
}

// A blob is reference-counted data that libncbi-vdb no longer touches once handed out, so it
// can move to the thread that reads it, like the cursor it came from.
unsafe impl Send for VBlob {}

impl VBlob {
    /// Raw cell data of row `row_id`, which must be one of the blob's rows.
    fn cell_data(&self, row_id: i64) -> Result<CellData, VdbError> {
        let mut elem_bits: u32 = 0;
        let mut base: *const std::ffi::c_void = ptr::null();
        let mut boff: u32 = 0;
        let mut row_len: u32 = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::VBlobCellData(
                self.ptr,
                row_id,
                &raw mut elem_bits,
                &raw mut base,
                &raw mut boff,
                &raw mut row_len,
            )
        };
        check_rc(rc)?;
        Ok(CellData { elem_bits, base, _boff: boff, row_len })
    }
}

impl Drop for VBlob {
    fn drop(&mut self) {
        unsafe { fg_sra_vdb_sys::VBlobRelease(self.ptr) };
    }
}

/// Raw cell data from `VCursorCellDataDirect` or `VBlobCellData`.
struct CellData {
    elem_bits: u32,
    base: *const std::ffi::c_void,
    _boff: u32,
    row_len: u32,
}

impl CellData {
    /// The cell of a row a column has no value for.
    const EMPTY: Self = Self { elem_bits: 0, base: ptr::null(), _boff: 0, row_len: 0 };

    /// The cell's first value as `T`, whose width is `expected_bits`; the default for an
    /// empty cell.
    fn scalar<T: Copy + Default>(&self, expected_bits: u32) -> Result<T, VdbError> {
        if self.row_len == 0 {
            return Ok(T::default());
        }
        check_elem_bits(self.elem_bits, expected_bits)?;
        Ok(unsafe { *self.base.cast::<T>() })
    }

    /// Replace `buf`'s contents with the cell's values as `T`, whose width is `expected_bits`.
    fn copy_into<T: Copy>(&self, expected_bits: u32, buf: &mut Vec<T>) -> Result<(), VdbError> {
        buf.clear();
        if self.row_len == 0 {
            return Ok(());
        }
        check_elem_bits(self.elem_bits, expected_bits)?;
        let values = unsafe { slice::from_raw_parts(self.base.cast::<T>(), self.row_len as usize) };
        buf.extend_from_slice(values);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_elem_bits_accepts_matching_width() {
        assert!(check_elem_bits(32, 32).is_ok());
        assert!(check_elem_bits(8, 8).is_ok());
    }

    #[test]
    fn check_elem_bits_rejects_narrower_cell() {
        // An 8-bit cell read as a 32-bit type would over-read: reject it.
        assert_eq!(
            check_elem_bits(8, 32),
            Err(VdbError::ElemBitsMismatch { expected: 32, actual: 8 })
        );
    }

    #[test]
    fn blob_rows_run_from_the_first_row_for_the_count_given() {
        assert_eq!(blob_rows(12, 10, 5), Ok(10..=14));
        assert_eq!(blob_rows(7, 7, 1), Ok(7..=7));
    }

    #[test]
    fn blob_without_the_row_read_is_refused() {
        let misses = |row, first, count| Err(VdbError::BlobMissesRow { row, first, count });
        assert_eq!(blob_rows(10, 10, 0), misses(10, 10, 0));
        assert_eq!(blob_rows(15, 10, 5), misses(15, 10, 5));
        assert_eq!(blob_rows(9, 10, 5), misses(9, 10, 5));
        assert_eq!(blob_rows(i64::MAX, 2, u64::MAX), misses(i64::MAX, 2, u64::MAX));
        assert_eq!(blob_rows(i64::MAX, i64::MAX, 2), misses(i64::MAX, i64::MAX, 2));
    }

    #[test]
    fn check_elem_bits_rejects_wider_cell() {
        assert_eq!(
            check_elem_bits(32, 8),
            Err(VdbError::ElemBitsMismatch { expected: 8, actual: 32 })
        );
    }
}
