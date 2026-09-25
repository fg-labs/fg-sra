//! Safe wrapper for `VDBManager`.
//!
//! The VDB manager is the entry point for opening databases and tables.
//! It is created once and used to open all VDB resources.

use std::ptr;

use crate::database::{VDatabase, VTable};
use crate::error::{LITERAL_FORMAT, VdbError, check_rc, path_to_cstring};

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
        let rc = unsafe { fg_sra_vdb_sys::VDBManagerMakeRead(&raw mut mgr, ptr::null()) };
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

    /// Report what kind of VDB object `path` names: a database, a flat table, or neither.
    ///
    /// Like the `open_*` methods this resolves accessions, so a bare accession may cost a
    /// network lookup. Opening a database as a table, or a table as a database, fails, so
    /// call this first when the input could be either. Fails with
    /// [`VdbError::PercentInPath`] for a path containing `%`.
    pub fn path_type(&self, path: &str) -> Result<PathType, VdbError> {
        let c_path = path_to_cstring(path)?;
        // Safety: VDBManagerPathType is printf-style variadic; the path is the single
        // argument to a literal "%s" format.
        let raw = unsafe {
            fg_sra_vdb_sys::VDBManagerPathType(self.ptr, LITERAL_FORMAT.as_ptr(), c_path.as_ptr())
        };
        Ok(PathType::from_raw(raw))
    }

    /// Open a read-only database by accession or path.
    ///
    /// Fails with [`VdbError::PercentInPath`] for a path containing `%`, which ncbi-vdb
    /// cannot handle.
    pub fn open_db_read(&self, path: &str) -> Result<VDatabase, VdbError> {
        let c_path = path_to_cstring(path)?;
        let mut db: *const fg_sra_vdb_sys::VDatabase = ptr::null();
        // Safety: VDBManagerOpenDBRead is printf-style variadic; the path is the single
        // argument to a literal "%s" format.
        let rc = unsafe {
            fg_sra_vdb_sys::VDBManagerOpenDBRead(
                self.ptr,
                &raw mut db,
                ptr::null(), // schema (NULL = use default)
                LITERAL_FORMAT.as_ptr(),
                c_path.as_ptr(),
            )
        };
        check_rc(rc)?;
        Ok(VDatabase::from_raw(db))
    }

    /// Open a read-only flat table (e.g. an older unaligned SRA run) by accession or path.
    ///
    /// This fails on a database; open its `SEQUENCE` table through [`VDatabase`] instead.
    /// Like [`open_db_read`](Self::open_db_read), it refuses a path containing `%`.
    pub fn open_table_read(&self, path: &str) -> Result<VTable, VdbError> {
        let c_path = path_to_cstring(path)?;
        let mut tbl: *const fg_sra_vdb_sys::VTable = ptr::null();
        // Safety: VDBManagerOpenTableRead is printf-style variadic; the path is the single
        // argument to a literal "%s" format.
        let rc = unsafe {
            fg_sra_vdb_sys::VDBManagerOpenTableRead(
                self.ptr,
                &raw mut tbl,
                ptr::null(), // schema (NULL = use the table's own)
                LITERAL_FORMAT.as_ptr(),
                c_path.as_ptr(),
            )
        };
        check_rc(rc)?;
        Ok(VTable::from_raw(tbl))
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

/// The kind of object at a path, from `VDBManagerPathType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathType {
    /// A VDB database: a cSRA archive, or a database holding only a `SEQUENCE` table.
    Database,
    /// A flat VDB table, as many older unaligned SRA runs are.
    Table,
    /// Nothing exists at the path, or the accession could not be resolved.
    NotFound,
    /// Anything else (a plain file or directory, a bad path, a lone column, …), with the
    /// raw `KDBPathType` value and its alias bit cleared.
    Other(u32),
}

impl PathType {
    /// Interpret a raw `VDBManagerPathType` result, ignoring the alias bit (set when the
    /// path is reached through a symbolic link).
    fn from_raw(raw: std::os::raw::c_int) -> Self {
        let Ok(raw) = u32::try_from(raw) else {
            return Self::Other(fg_sra_vdb_sys::kptBadPath);
        };
        match raw & !fg_sra_vdb_sys::kptAlias {
            fg_sra_vdb_sys::kptDatabase => Self::Database,
            fg_sra_vdb_sys::kptTable => Self::Table,
            fg_sra_vdb_sys::kptNotFound => Self::NotFound,
            other => Self::Other(other),
        }
    }
}

/// Stop VDB from resolving or reading anything over the network, for the rest of the process.
///
/// Accessions and references are then found only locally (the working directory, local
/// repositories, or beside the archive), and anything else fails. The setting is
/// process-wide (a global in ncbi-vdb's resolver), so it applies to managers made before
/// and after the call.
pub fn disable_remote_access() -> Result<(), VdbError> {
    let mut vfs: *mut fg_sra_vdb_sys::VFSManager = ptr::null_mut();
    // Safety: VFSManagerMake initializes `vfs`, which is released below.
    check_rc(unsafe { fg_sra_vdb_sys::VFSManagerMake(&raw mut vfs) })?;
    let mut resolver: *mut fg_sra_vdb_sys::VResolver = ptr::null_mut();
    // Safety: `vfs` is valid; VFSManagerGetResolver initializes `resolver`.
    let rc = unsafe { fg_sra_vdb_sys::VFSManagerGetResolver(vfs, &raw mut resolver) };
    if rc == 0 {
        // Safety: `resolver` is valid. The call returns the prior state, not an rc.
        unsafe {
            fg_sra_vdb_sys::VResolverRemoteEnable(resolver, fg_sra_vdb_sys::vrAlwaysDisable);
            fg_sra_vdb_sys::VResolverRelease(resolver);
        }
    }
    // Safety: we own the reference made above.
    unsafe { fg_sra_vdb_sys::VFSManagerRelease(vfs) };
    check_rc(rc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn database_path_type_is_recognised() {
        let raw = std::os::raw::c_int::try_from(fg_sra_vdb_sys::kptDatabase).unwrap();
        assert_eq!(PathType::from_raw(raw), PathType::Database);
    }

    #[test]
    fn table_reached_through_a_symlink_is_still_a_table() {
        let raw = fg_sra_vdb_sys::kptTable | fg_sra_vdb_sys::kptAlias;
        let raw = std::os::raw::c_int::try_from(raw).unwrap();
        assert_eq!(PathType::from_raw(raw), PathType::Table);
    }

    #[test]
    fn plain_directory_is_other() {
        let raw = std::os::raw::c_int::try_from(fg_sra_vdb_sys::kptDir).unwrap();
        assert_eq!(PathType::from_raw(raw), PathType::Other(fg_sra_vdb_sys::kptDir));
    }

    #[test]
    fn negative_path_type_is_a_bad_path() {
        assert_eq!(PathType::from_raw(-1), PathType::Other(fg_sra_vdb_sys::kptBadPath));
    }
}
