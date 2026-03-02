//! Safe wrappers for `ReferenceList` and `ReferenceObj`.
//!
//! Provides access to reference sequence metadata needed for SAM header
//! generation and coordinate-sorted alignment iteration.

use std::ptr;

use crate::database::VDatabase;
use crate::error::{VdbError, check_rc, to_cstring};

/// Option flags for [`ReferenceList::make_database`].
///
/// These control which alignment ID columns are read from the REFERENCE table.
pub mod reflist_options {
    /// Include `PRIMARY_ALIGNMENT_IDS` column (required for primary alignment iteration).
    pub const USE_PRIMARY_IDS: u32 = 0x02;
    /// Include `SECONDARY_ALIGNMENT_IDS` column (required for secondary alignment iteration).
    pub const USE_SECONDARY_IDS: u32 = 0x04;
    /// Include `EVIDENCE_INTERVAL_IDS` column (required for evidence interval iteration).
    pub const USE_EVIDENCE_IDS: u32 = 0x08;
}

/// Safe wrapper around the VDB `ReferenceList` opaque type.
///
/// Created from a database, provides access to reference sequences.
pub struct ReferenceList {
    ptr: *const fg_sra_vdb_sys::ReferenceList,
}

unsafe impl Send for ReferenceList {}

impl ReferenceList {
    /// Create a `ReferenceList` from a database.
    ///
    /// - `options`: bitmask of options (0 for defaults)
    /// - `cache`: cursor cache size in bytes (0 for default)
    pub fn make_database(db: &VDatabase, options: u32, cache: usize) -> Result<Self, VdbError> {
        let mut reflist: *const fg_sra_vdb_sys::ReferenceList = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::ReferenceList_MakeDatabase(
                &mut reflist,
                db.as_ptr(),
                options,
                cache,
                ptr::null(), // name filter
                0,           // numbins
            )
        };
        check_rc(rc)?;
        Ok(Self { ptr: reflist })
    }

    /// Get the number of references.
    pub fn count(&self) -> Result<u32, VdbError> {
        let mut count: u32 = 0;
        let rc = unsafe { fg_sra_vdb_sys::ReferenceList_Count(self.ptr, &mut count) };
        check_rc(rc)?;
        Ok(count)
    }

    /// Get a reference by index.
    pub fn get(&self, idx: u32) -> Result<ReferenceObj, VdbError> {
        let mut obj: *const fg_sra_vdb_sys::ReferenceObj = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::ReferenceList_Get(self.ptr, &mut obj, idx) };
        check_rc(rc)?;
        Ok(ReferenceObj { ptr: obj })
    }

    /// Find a reference by name.
    pub fn find(&self, name: &str) -> Result<ReferenceObj, VdbError> {
        let c_name = to_cstring(name)?;
        let mut obj: *const fg_sra_vdb_sys::ReferenceObj = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::ReferenceList_Find(
                self.ptr,
                &mut obj,
                c_name.as_ptr(),
                c_name.as_bytes().len(),
            )
        };
        check_rc(rc)?;
        Ok(ReferenceObj { ptr: obj })
    }

    /// Iterate over all references.
    pub fn iter(&self) -> Result<ReferenceIter<'_>, VdbError> {
        let count = self.count()?;
        Ok(ReferenceIter { reflist: self, idx: 0, count })
    }

    #[allow(dead_code)]
    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::ReferenceList {
        self.ptr
    }
}

impl Drop for ReferenceList {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // ReferenceList_Release returns void.
            unsafe { fg_sra_vdb_sys::ReferenceList_Release(self.ptr) };
        }
    }
}

/// Iterator over references in a `ReferenceList`.
pub struct ReferenceIter<'a> {
    reflist: &'a ReferenceList,
    idx: u32,
    count: u32,
}

impl<'a> Iterator for ReferenceIter<'a> {
    type Item = Result<ReferenceObj, VdbError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.count {
            return None;
        }
        let result = self.reflist.get(self.idx);
        self.idx += 1;
        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = (self.count - self.idx) as usize;
        (remaining, Some(remaining))
    }
}

/// Safe wrapper around a single `ReferenceObj`.
///
/// Provides access to reference name, sequence ID, and length.
/// Note: `ReferenceObj_Release` returns `void`, not `rc_t`.
pub struct ReferenceObj {
    ptr: *const fg_sra_vdb_sys::ReferenceObj,
}

unsafe impl Send for ReferenceObj {}

impl ReferenceObj {
    /// Get the reference name (e.g., "chr1").
    pub fn name(&self) -> Result<String, VdbError> {
        ref_obj_name(self.ptr)
    }

    /// Get the sequence ID (e.g., "NC_000001.11").
    pub fn seq_id(&self) -> Result<String, VdbError> {
        ref_obj_seq_id(self.ptr)
    }

    /// Get the sequence length.
    pub fn seq_length(&self) -> Result<u32, VdbError> {
        let mut len: u32 = 0;
        let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_SeqLength(self.ptr, &mut len) };
        check_rc(rc)?;
        Ok(len)
    }

    /// Get the reference index.
    pub fn idx(&self) -> Result<u32, VdbError> {
        let mut idx: u32 = 0;
        let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_Idx(self.ptr, &mut idx) };
        check_rc(rc)?;
        Ok(idx)
    }

    pub(crate) fn as_ptr(&self) -> *const fg_sra_vdb_sys::ReferenceObj {
        self.ptr
    }
}

impl Drop for ReferenceObj {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // ReferenceObj_Release returns void.
            unsafe { fg_sra_vdb_sys::ReferenceObj_Release(self.ptr) };
        }
    }
}

/// Get the name from a raw `ReferenceObj` pointer.
///
/// Shared by `ReferenceObj::name()` and `NextReference::ref_name()`.
pub(crate) fn ref_obj_name(ptr: *const fg_sra_vdb_sys::ReferenceObj) -> Result<String, VdbError> {
    let mut name: *const std::os::raw::c_char = ptr::null();
    let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_Name(ptr, &mut name) };
    check_rc(rc)?;
    if name.is_null() {
        return Ok(String::new());
    }
    let s = unsafe { std::ffi::CStr::from_ptr(name) };
    Ok(s.to_string_lossy().into_owned())
}

/// Get the sequence ID from a raw `ReferenceObj` pointer.
///
/// Shared by `ReferenceObj::seq_id()` and `NextReference::seq_id()`.
pub(crate) fn ref_obj_seq_id(ptr: *const fg_sra_vdb_sys::ReferenceObj) -> Result<String, VdbError> {
    let mut seqid: *const std::os::raw::c_char = ptr::null();
    let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_SeqId(ptr, &mut seqid) };
    check_rc(rc)?;
    if seqid.is_null() {
        return Ok(String::new());
    }
    let s = unsafe { std::ffi::CStr::from_ptr(seqid) };
    Ok(s.to_string_lossy().into_owned())
}
