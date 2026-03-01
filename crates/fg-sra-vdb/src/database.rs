//! Safe wrappers for `VDatabase` and `VTable`.
//!
//! RAII types that call `VDatabaseRelease` / `VTableRelease` on drop.

use std::ptr;

use crate::cursor::VCursor;
use crate::error::{VdbError, check_rc, to_cstring};

/// Safe wrapper around the VDB `VDatabase` opaque type.
///
/// Provides access to tables and metadata within a VDB database.
pub struct VDatabase {
    ptr: *const fg_sra_vdb_sys::VDatabase,
}

unsafe impl Send for VDatabase {}

impl VDatabase {
    /// Create a `VDatabase` from a raw pointer (takes ownership).
    pub(crate) fn from_raw(ptr: *const fg_sra_vdb_sys::VDatabase) -> Self {
        Self { ptr }
    }

    /// Open a read-only table within this database.
    pub fn open_table_read(&self, name: &str) -> Result<VTable, VdbError> {
        let c_name = to_cstring(name)?;
        let mut tbl: *const fg_sra_vdb_sys::VTable = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::VDatabaseOpenTableRead(self.ptr, &mut tbl, c_name.as_ptr())
        };
        check_rc(rc)?;
        Ok(VTable::from_raw(tbl))
    }

    /// Check if a table exists in this database.
    pub fn has_table(&self, name: &str) -> bool {
        self.open_table_read(name).is_ok()
    }

    /// List all table names in this database.
    pub fn list_tables(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        let rc = unsafe { fg_sra_vdb_sys::VDatabaseListTbl(self.ptr, &mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    /// Open the database metadata for reading.
    pub fn open_metadata_read(&self) -> Result<KMetadata, VdbError> {
        let mut meta: *const fg_sra_vdb_sys::KMetadata = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::VDatabaseOpenMetadataRead(self.ptr, &mut meta) };
        check_rc(rc)?;
        Ok(KMetadata { ptr: meta })
    }

    /// Get the raw pointer (for passing to C APIs).
    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::VDatabase {
        self.ptr
    }
}

impl Drop for VDatabase {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::VDatabaseRelease(self.ptr) };
        }
    }
}

/// Safe wrapper around the VDB `VTable` opaque type.
pub struct VTable {
    ptr: *const fg_sra_vdb_sys::VTable,
}

unsafe impl Send for VTable {}

impl VTable {
    pub(crate) fn from_raw(ptr: *const fg_sra_vdb_sys::VTable) -> Self {
        Self { ptr }
    }

    /// Create a read cursor on this table.
    pub fn create_cursor_read(&self) -> Result<VCursor, VdbError> {
        let mut curs: *const fg_sra_vdb_sys::VCursor = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::VTableCreateCursorRead(self.ptr, &mut curs) };
        check_rc(rc)?;
        Ok(VCursor::from_raw(curs))
    }

    /// Create a cached read cursor with the given cache capacity in bytes.
    pub fn create_cached_cursor_read(&self, capacity: usize) -> Result<VCursor, VdbError> {
        let mut curs: *const fg_sra_vdb_sys::VCursor = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::VTableCreateCachedCursorRead(self.ptr, &mut curs, capacity)
        };
        check_rc(rc)?;
        Ok(VCursor::from_raw(curs))
    }

    /// List readable column names.
    pub fn list_readable_columns(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        let rc = unsafe { fg_sra_vdb_sys::VTableListReadableColumns(self.ptr, &mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    #[allow(dead_code)]
    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::VTable {
        self.ptr
    }
}

impl Drop for VTable {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::VTableRelease(self.ptr) };
        }
    }
}

/// Safe wrapper around `KMetadata`.
pub struct KMetadata {
    ptr: *const fg_sra_vdb_sys::KMetadata,
}

impl KMetadata {
    /// Open a metadata node by path.
    pub fn open_node_read(&self, path: &str) -> Result<KMDataNode, VdbError> {
        let c_path = to_cstring(path)?;
        let mut node: *const fg_sra_vdb_sys::KMDataNode = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::KMetadataOpenNodeRead(self.ptr, &mut node, c_path.as_ptr())
        };
        check_rc(rc)?;
        Ok(KMDataNode { ptr: node })
    }
}

impl Drop for KMetadata {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::KMetadataRelease(self.ptr) };
        }
    }
}

/// Safe wrapper around `KMDataNode`.
pub struct KMDataNode {
    ptr: *const fg_sra_vdb_sys::KMDataNode,
}

impl KMDataNode {
    /// Read data from this metadata node.
    ///
    /// Returns the bytes read. `offset` is the byte position to start reading.
    pub fn read(&self, offset: usize, buffer: &mut [u8]) -> Result<(usize, usize), VdbError> {
        let mut num_read: usize = 0;
        let mut remaining: usize = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::KMDataNodeRead(
                self.ptr,
                offset,
                buffer.as_mut_ptr() as *mut std::ffi::c_void,
                buffer.len(),
                &mut num_read,
                &mut remaining,
            )
        };
        check_rc(rc)?;
        Ok((num_read, remaining))
    }

    /// Read the entire node content as a String.
    pub fn read_all(&self) -> Result<String, VdbError> {
        let mut result = Vec::new();
        let mut offset = 0usize;
        let mut buf = [0u8; 4096];
        loop {
            let (num_read, _remaining) = self.read(offset, &mut buf)?;
            if num_read == 0 {
                break;
            }
            result.extend_from_slice(&buf[..num_read]);
            offset += num_read;
        }
        Ok(String::from_utf8_lossy(&result).into_owned())
    }
}

impl Drop for KMDataNode {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::KMDataNodeRelease(self.ptr) };
        }
    }
}

/// Read all strings from a KNamelist.
fn read_namelist(names: *const fg_sra_vdb_sys::KNamelist) -> Result<Vec<String>, VdbError> {
    let mut count: u32 = 0;
    let rc = unsafe { fg_sra_vdb_sys::KNamelistCount(names, &mut count) };
    check_rc(rc)?;

    let mut result = Vec::with_capacity(count as usize);
    for i in 0..count {
        let mut name: *const std::os::raw::c_char = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::KNamelistGet(names, i, &mut name) };
        check_rc(rc)?;
        if !name.is_null() {
            let s = unsafe { std::ffi::CStr::from_ptr(name) };
            result.push(s.to_string_lossy().into_owned());
        }
    }
    Ok(result)
}
