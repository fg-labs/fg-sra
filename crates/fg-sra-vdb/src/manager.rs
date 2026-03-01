//! Safe wrapper for `VDBManager`.
//!
//! The VDB manager is the entry point for opening databases and tables.
//! It is created once and used to open all VDB resources.

use std::ptr;

use crate::database::VDatabase;
use crate::error::{VdbError, check_rc, to_cstring};

/// Safe wrapper around the VDB `VDBManager` opaque type.
///
/// Created via [`VdbManager::make_read`]. The manager is reference-counted
/// and released on drop.
pub struct VdbManager {
    ptr: *const fg_sra_vdb_sys::VDBManager,
}

// VDBManager is internally reference-counted and thread-safe for read operations.
unsafe impl Send for VdbManager {}
unsafe impl Sync for VdbManager {}

impl VdbManager {
    /// Create a new read-only VDB manager.
    ///
    /// This is the main entry point for accessing VDB databases. Pass `None`
    /// for the working directory to use the default.
    pub fn make_read() -> Result<Self, VdbError> {
        let mut mgr: *const fg_sra_vdb_sys::VDBManager = ptr::null();
        // Safety: VDBManagerMakeRead initializes mgr; we pass NULL for default directory.
        let rc = unsafe { fg_sra_vdb_sys::VDBManagerMakeRead(&mut mgr, ptr::null()) };
        check_rc(rc)?;
        Ok(Self { ptr: mgr })
    }

    /// Disable the background pagemap pre-computation thread.
    ///
    /// This is useful when running our own thread pool to avoid contention.
    pub fn disable_pagemap_thread(&self) -> Result<(), VdbError> {
        // Safety: self.ptr is valid as long as self is alive.
        let rc = unsafe { fg_sra_vdb_sys::VDBManagerDisablePagemapThread(self.ptr) };
        check_rc(rc)
    }

    /// Open a read-only database by accession or path.
    pub fn open_db_read(&self, path: &str) -> Result<VDatabase, VdbError> {
        let c_path = to_cstring(path)?;
        let mut db: *const fg_sra_vdb_sys::VDatabase = ptr::null();
        // Safety: VDBManagerOpenDBRead is variadic; we pass the path as the only
        // variadic arg (format string with no format specifiers).
        let rc = unsafe {
            fg_sra_vdb_sys::VDBManagerOpenDBRead(
                self.ptr,
                &mut db,
                ptr::null(), // schema (NULL = use default)
                c_path.as_ptr(),
            )
        };
        check_rc(rc)?;
        Ok(VDatabase::from_raw(db))
    }

    /// Get the raw pointer (for passing to C APIs that need the manager).
    #[allow(dead_code)]
    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::VDBManager {
        self.ptr
    }
}

impl Drop for VdbManager {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // Safety: we own the reference.
            unsafe { fg_sra_vdb_sys::VDBManagerRelease(self.ptr) };
        }
    }
}
