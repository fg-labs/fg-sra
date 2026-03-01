//! Safe wrapper for `VCursor` with typed column reads.
//!
//! Provides methods for adding columns, opening the cursor, and reading
//! typed cell data via `VCursorCellDataDirect`.

use std::ptr;
use std::slice;

use crate::error::{VdbError, check_rc, to_cstring};

/// Sentinel value for columns that were never added or are not available.
pub const INVALID_COLUMN: u32 = 0xFFFFFFFF;

/// Safe wrapper around the VDB `VCursor` opaque type.
///
/// A cursor is opened on a table and reads column data row-by-row.
/// Columns must be added before opening, and data is read via
/// `VCursorCellDataDirect` (zero-copy access into VDB page cache).
pub struct VCursor {
    ptr: *const fg_sra_vdb_sys::VCursor,
}

unsafe impl Send for VCursor {}

impl VCursor {
    pub(crate) fn from_raw(ptr: *const fg_sra_vdb_sys::VCursor) -> Self {
        Self { ptr }
    }

    /// Add a column to an unopened cursor.
    ///
    /// Returns the column index used for subsequent reads.
    /// Column names should include the type cast, e.g. `"(I64)SEQ_SPOT_ID"`.
    pub fn add_column(&self, name: &str) -> Result<u32, VdbError> {
        let c_name = to_cstring(name)?;
        let mut idx: u32 = 0;
        let rc = unsafe { fg_sra_vdb_sys::VCursorAddColumn(self.ptr, &mut idx, c_name.as_ptr()) };
        check_rc(rc)?;
        Ok(idx)
    }

    /// Add a column, returning `None` if the column does not exist.
    ///
    /// This is equivalent to `add_opt_column` in the C implementation.
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
        let rc = unsafe { fg_sra_vdb_sys::VCursorIdRange(self.ptr, col_idx, &mut first, &mut count) };
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
        let data = self.cell_data_direct(row_id, col_idx)?;
        if data.row_len == 0 {
            return Ok(T::default());
        }
        debug_assert_eq!(data.elem_bits, expected_bits);
        Ok(unsafe { *(data.base as *const T) })
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
        debug_assert_eq!(data.elem_bits, expected_bits);
        let values =
            unsafe { slice::from_raw_parts(data.base as *const T, data.row_len as usize) };
        Ok(values.to_vec())
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
        debug_assert_eq!(data.elem_bits, 8);
        let bytes =
            unsafe { slice::from_raw_parts(data.base as *const u8, data.row_len as usize) };
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }

    /// Read a slice of `u8` values (e.g., quality scores, read data).
    pub fn read_u8_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<u8>, VdbError> {
        self.read_slice(row_id, col_idx, 8)
    }

    /// Read a slice of `i64` values.
    pub fn read_i64_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<i64>, VdbError> {
        self.read_slice(row_id, col_idx, 64)
    }

    /// Read a slice of `u32` values (e.g., INSDC_coord_len).
    pub fn read_u32_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<u32>, VdbError> {
        self.read_slice(row_id, col_idx, 32)
    }

    /// Read a slice of `i32` values (e.g., INSDC_coord_zero).
    pub fn read_i32_slice(&self, row_id: i64, col_idx: u32) -> Result<Vec<i32>, VdbError> {
        self.read_slice(row_id, col_idx, 32)
    }

    /// Read an INSDC_coord_zero value (0-based coordinate, i32).
    pub fn read_coord_zero(&self, row_id: i64, col_idx: u32) -> Result<i32, VdbError> {
        self.read_i32(row_id, col_idx)
    }

    /// Read an INSDC_coord_len value (length, u32).
    pub fn read_coord_len(&self, row_id: i64, col_idx: u32) -> Result<u32, VdbError> {
        self.read_u32(row_id, col_idx)
    }

    /// Low-level: read raw cell data via `VCursorCellDataDirect`.
    ///
    /// Returns a `CellData` with a pointer into VDB's page cache.
    fn cell_data_direct(&self, row_id: i64, col_idx: u32) -> Result<CellData, VdbError> {
        let mut elem_bits: u32 = 0;
        let mut base: *const std::ffi::c_void = ptr::null();
        let mut boff: u32 = 0;
        let mut row_len: u32 = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::VCursorCellDataDirect(
                self.ptr,
                row_id,
                col_idx,
                &mut elem_bits,
                &mut base,
                &mut boff,
                &mut row_len,
            )
        };
        check_rc(rc)?;
        Ok(CellData { elem_bits, base, _boff: boff, row_len })
    }

    /// Get the raw cursor pointer (for passing to PlacementIterator creation).
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

/// Raw cell data from VCursorCellDataDirect.
struct CellData {
    elem_bits: u32,
    base: *const std::ffi::c_void,
    _boff: u32,
    row_len: u32,
}
