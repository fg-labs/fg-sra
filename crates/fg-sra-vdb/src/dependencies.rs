//! Safe wrappers for `VDBDependencies`.
//!
//! Provides access to VDB database dependency information and triggers
//! reference sequence cache population as a side effect.

use std::ptr;

use crate::database::VDatabase;
use crate::error::{VdbError, check_rc};
use crate::retry::retry_on_network_error;

/// Information about a single dependency (reference sequence).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DependencyInfo {
    /// The sequence ID of this dependency (e.g., "CM000667.1").
    pub seq_id: String,
    /// Whether the dependency is available locally (cached).
    pub local: bool,
}

/// Safe wrapper around the VDB `VDBDependencies` opaque type.
///
/// Created via [`VDatabase::list_dependencies`]. Listing dependencies triggers
/// VDB's internal cache population for reference sequences.
pub struct VdbDependencies {
    ptr: *const fg_sra_vdb_sys::VDBDependencies,
}

unsafe impl Send for VdbDependencies {}

impl VdbDependencies {
    /// List dependencies of a database, optionally only those not yet cached.
    ///
    /// When `missing_only` is `false`, lists all dependencies.
    /// When `missing_only` is `true`, lists only those not yet cached locally.
    ///
    /// As a side effect, this call resolves and caches reference sequences
    /// if VDB caching is enabled.
    pub(crate) fn list(db: &VDatabase, missing_only: bool) -> Result<Self, VdbError> {
        retry_on_network_error("VDatabaseListDependencies", || {
            let mut dep: *const fg_sra_vdb_sys::VDBDependencies = ptr::null();
            let rc = unsafe {
                fg_sra_vdb_sys::VDatabaseListDependencies(db.as_ptr(), &raw mut dep, missing_only)
            };
            check_rc(rc)?;
            Ok(Self { ptr: dep })
        })
    }

    /// Get the number of dependencies.
    pub fn count(&self) -> Result<u32, VdbError> {
        let mut count: u32 = 0;
        let rc = unsafe { fg_sra_vdb_sys::VDBDependenciesCount(self.ptr, &raw mut count) };
        check_rc(rc)?;
        Ok(count)
    }

    /// Get the sequence ID of the dependency at `idx`.
    pub fn seq_id(&self, idx: u32) -> Result<String, VdbError> {
        let mut seq_id: *const std::os::raw::c_char = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::VDBDependenciesSeqId(self.ptr, &raw mut seq_id, idx) };
        check_rc(rc)?;
        if seq_id.is_null() {
            return Ok(String::new());
        }
        let s = unsafe { std::ffi::CStr::from_ptr(seq_id) };
        Ok(s.to_string_lossy().into_owned())
    }

    /// Get whether the dependency at `idx` is available locally.
    pub fn local(&self, idx: u32) -> Result<bool, VdbError> {
        let mut local: bool = false;
        let rc = unsafe { fg_sra_vdb_sys::VDBDependenciesLocal(self.ptr, &raw mut local, idx) };
        check_rc(rc)?;
        Ok(local)
    }

    /// Collect information about all dependencies.
    pub fn all_info(&self) -> Result<Vec<DependencyInfo>, VdbError> {
        let count = self.count()?;
        let mut infos = Vec::with_capacity(count as usize);
        for idx in 0..count {
            infos.push(DependencyInfo { seq_id: self.seq_id(idx)?, local: self.local(idx)? });
        }
        Ok(infos)
    }
}

impl Drop for VdbDependencies {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            let _lock = crate::manager::lifecycle_lock();
            unsafe { fg_sra_vdb_sys::VDBDependenciesRelease(self.ptr) };
        }
    }
}
