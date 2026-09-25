//! Safe wrappers for `VDatabase` and `VTable`.
//!
//! RAII types that call `VDatabaseRelease` / `VTableRelease` on drop.

use std::ptr;

use crate::cursor::VCursor;
use crate::dependencies::VdbDependencies;
use crate::error::{LITERAL_FORMAT, VdbError, check_rc, to_cstring};

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
        // Safety: VDatabaseOpenTableRead is printf-style variadic; the name is the single
        // argument to a literal "%s" format.
        let rc = unsafe {
            fg_sra_vdb_sys::VDatabaseOpenTableRead(
                self.ptr,
                &raw mut tbl,
                LITERAL_FORMAT.as_ptr(),
                c_name.as_ptr(),
            )
        };
        check_rc(rc)?;
        Ok(VTable::from_raw(tbl))
    }

    /// Check if a table exists in this database.
    #[must_use]
    pub fn has_table(&self, name: &str) -> bool {
        self.open_table_read(name).is_ok()
    }

    /// List all table names in this database.
    pub fn list_tables(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        let rc = unsafe { fg_sra_vdb_sys::VDatabaseListTbl(self.ptr, &raw mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    /// Open the database metadata for reading.
    pub fn open_metadata_read(&self) -> Result<KMetadata, VdbError> {
        let mut meta: *const fg_sra_vdb_sys::KMetadata = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::VDatabaseOpenMetadataRead(self.ptr, &raw mut meta) };
        check_rc(rc)?;
        Ok(KMetadata { ptr: meta })
    }

    /// List dependencies of this database, triggering reference cache population.
    ///
    /// When `missing_only` is `false`, returns all dependencies.
    /// When `missing_only` is `true`, returns only those not yet cached locally.
    pub fn list_dependencies(&self, missing_only: bool) -> Result<VdbDependencies, VdbError> {
        VdbDependencies::list(self, missing_only)
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
        let rc = unsafe { fg_sra_vdb_sys::VTableCreateCursorRead(self.ptr, &raw mut curs) };
        check_rc(rc)?;
        Ok(VCursor::from_raw(curs, false))
    }

    /// Create a cached read cursor with the given cache capacity in bytes.
    pub fn create_cached_cursor_read(&self, capacity: usize) -> Result<VCursor, VdbError> {
        let mut curs: *const fg_sra_vdb_sys::VCursor = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::VTableCreateCachedCursorRead(self.ptr, &raw mut curs, capacity)
        };
        check_rc(rc)?;
        Ok(VCursor::from_raw(curs, capacity > 0))
    }

    /// List readable column names.
    pub fn list_readable_columns(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        let rc = unsafe { fg_sra_vdb_sys::VTableListReadableColumns(self.ptr, &raw mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    /// List the columns physically stored in this table.
    ///
    /// Unlike [`list_readable_columns`](Self::list_readable_columns), this excludes columns
    /// the schema computes from others, so it shows how the data was stored: e.g. whether a
    /// cSRA `SEQUENCE` table holds `CMP_READ` (bases of unaligned reads only) or `READ`.
    pub fn list_physical_columns(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        // Safety: `self.ptr` is a valid table; on success `names` is a namelist we release.
        let rc = unsafe { fg_sra_vdb_sys::VTableListPhysColumns(self.ptr, &raw mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    /// Open this table's metadata for reading.
    ///
    /// Run statistics (`STATS/TABLE/SPOT_COUNT`, `BASE_COUNT`, …) live here, in the metadata
    /// of the flat table or of a database's `SEQUENCE` table, not in the database's.
    pub fn open_metadata_read(&self) -> Result<KMetadata, VdbError> {
        let mut meta: *const fg_sra_vdb_sys::KMetadata = ptr::null();
        // Safety: `self.ptr` is a valid table; on success `meta` is owned by the `KMetadata`.
        let rc = unsafe { fg_sra_vdb_sys::VTableOpenMetadataRead(self.ptr, &raw mut meta) };
        check_rc(rc)?;
        Ok(KMetadata { ptr: meta })
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
        // Safety: KMetadataOpenNodeRead is printf-style variadic; the path is the single
        // argument to a literal "%s" format.
        let rc = unsafe {
            fg_sra_vdb_sys::KMetadataOpenNodeRead(
                self.ptr,
                &raw mut node,
                LITERAL_FORMAT.as_ptr(),
                c_path.as_ptr(),
            )
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
                buffer.as_mut_ptr().cast::<std::ffi::c_void>(),
                buffer.len(),
                &raw mut num_read,
                &raw mut remaining,
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

    /// List the names of this node's children (e.g. `PHRED_30` under `STATS/QUALITY`).
    pub fn list_children(&self) -> Result<Vec<String>, VdbError> {
        let mut names: *mut fg_sra_vdb_sys::KNamelist = ptr::null_mut();
        // Safety: `self.ptr` is a valid node; on success `names` is a namelist we release.
        let rc = unsafe { fg_sra_vdb_sys::KMDataNodeListChildren(self.ptr, &raw mut names) };
        check_rc(rc)?;
        let result = read_namelist(names);
        unsafe { fg_sra_vdb_sys::KNamelistRelease(names) };
        result
    }

    /// Read this node's value as an unsigned integer.
    ///
    /// Handles values stored as 1, 2, 4 or 8 bytes in either byte order, as the stats
    /// nodes (`STATS/TABLE/SPOT_COUNT`, …) are; fails on a node of any other size.
    pub fn read_u64(&self) -> Result<u64, VdbError> {
        let mut value: u64 = 0;
        // Safety: `self.ptr` is a valid node and `value` a writable u64.
        let rc = unsafe { fg_sra_vdb_sys::KMDataNodeReadAsU64(self.ptr, &raw mut value) };
        check_rc(rc)?;
        Ok(value)
    }

    /// Read the attribute `name` of this node (e.g. `name` on `SOFTWARE/delite`).
    ///
    /// Fails with an error for which [`VdbError::is_not_found`] is true if the node has no
    /// such attribute.
    pub fn read_attr(&self, name: &str) -> Result<String, VdbError> {
        let c_name = to_cstring(name)?;
        let mut buffer = vec![0u8; 256];
        loop {
            let mut size: usize = 0;
            // Safety: `self.ptr` is a valid node, and the buffer pointer and length describe
            // `buffer`, which the call writes at most `buffer.len()` bytes of.
            let rc = unsafe {
                fg_sra_vdb_sys::KMDataNodeReadAttr(
                    self.ptr,
                    c_name.as_ptr(),
                    buffer.as_mut_ptr().cast::<std::os::raw::c_char>(),
                    buffer.len(),
                    &raw mut size,
                )
            };
            if rc == 0 {
                buffer.truncate(size);
                return Ok(String::from_utf8_lossy(&buffer).into_owned());
            }
            // A buffer too small for the value and its NUL terminator reports the value's
            // length in `size`; every other failure (e.g. no such attribute) reports zero.
            if size < buffer.len() {
                return Err(VdbError::new(rc));
            }
            buffer.resize(size + 1, 0);
        }
    }
}

impl Drop for KMDataNode {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::KMDataNodeRelease(self.ptr) };
        }
    }
}

/// Read all strings from a `KNamelist`.
fn read_namelist(names: *const fg_sra_vdb_sys::KNamelist) -> Result<Vec<String>, VdbError> {
    let mut count: u32 = 0;
    let rc = unsafe { fg_sra_vdb_sys::KNamelistCount(names, &raw mut count) };
    check_rc(rc)?;

    let mut result = Vec::with_capacity(count as usize);
    for i in 0..count {
        let mut name: *const std::os::raw::c_char = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::KNamelistGet(names, i, &raw mut name) };
        check_rc(rc)?;
        if !name.is_null() {
            let s = unsafe { std::ffi::CStr::from_ptr(name) };
            result.push(s.to_string_lossy().into_owned());
        }
    }
    Ok(result)
}
