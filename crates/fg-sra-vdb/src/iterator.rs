//! Safe wrappers for `AlignMgr` and `PlacementSetIterator`.
//!
//! These types support coordinate-sorted iteration over aligned reads
//! via the VDB placement iterator API.

use std::ptr;

use crate::cursor::VCursor;
use crate::error::{VdbError, check_rc};
use crate::reference::ReferenceObj;

/// Check an `rc_t` for the done/exhausted state, returning `Ok(None)` for done,
/// `Err(...)` for real errors, or `Ok(Some(()))` for success.
fn check_rc_done(rc: u32) -> Result<Option<()>, VdbError> {
    if rc == 0 {
        return Ok(Some(()));
    }
    let err = VdbError::new(rc);
    if err.is_done() { Ok(None) } else { Err(err) }
}

/// Alignment ID source selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum AlignIdSrc {
    /// Primary alignments.
    Primary = 0,
    /// Secondary alignments.
    Secondary = 1,
    /// Evidence alignments.
    Evidence = 2,
}

/// Safe wrapper around `AlignMgr`.
///
/// Used to create `PlacementSetIterator` instances.
pub struct AlignMgr {
    ptr: *const fg_sra_vdb_sys::AlignMgr,
}

unsafe impl Send for AlignMgr {}

impl AlignMgr {
    /// Create a new read-only alignment manager.
    pub fn make_read() -> Result<Self, VdbError> {
        let mut mgr: *const fg_sra_vdb_sys::AlignMgr = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::AlignMgrMakeRead(&mut mgr) };
        check_rc(rc)?;
        Ok(Self { ptr: mgr })
    }

    /// Create a new `PlacementSetIterator`.
    pub fn make_placement_set_iterator(&self) -> Result<PlacementSetIterator, VdbError> {
        let mut iter: *mut fg_sra_vdb_sys::PlacementSetIterator = ptr::null_mut();
        let rc = unsafe {
            fg_sra_vdb_sys::AlignMgrMakePlacementSetIterator(self.ptr, &mut iter)
        };
        check_rc(rc)?;
        Ok(PlacementSetIterator { ptr: iter })
    }
}

impl Drop for AlignMgr {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::AlignMgrRelease(self.ptr) };
        }
    }
}

/// Safe wrapper around `PlacementSetIterator`.
///
/// Provides coordinate-sorted iteration over aligned reads across references.
/// The iteration pattern is: references → windows → positions → records.
pub struct PlacementSetIterator {
    ptr: *mut fg_sra_vdb_sys::PlacementSetIterator,
}

unsafe impl Send for PlacementSetIterator {}

impl PlacementSetIterator {
    /// Add a placement iterator (for one alignment table on one reference).
    ///
    /// Takes ownership of the `PlacementIterator`.
    pub fn add_placement_iterator(
        &mut self,
        iter: PlacementIterator,
    ) -> Result<(), VdbError> {
        let rc = unsafe {
            fg_sra_vdb_sys::PlacementSetIteratorAddPlacementIterator(self.ptr, iter.ptr)
        };
        // Don't drop the iterator - ownership transferred to the set iterator.
        std::mem::forget(iter);
        check_rc(rc)
    }

    /// Advance to the next reference.
    ///
    /// Returns `Ok(Some(ref_obj))` with the new reference, or `Ok(None)` when
    /// all references are exhausted.
    pub fn next_reference(
        &mut self,
    ) -> Result<Option<NextReference>, VdbError> {
        let mut first_pos: i32 = 0;
        let mut len: u32 = 0;
        let mut ref_obj: *const fg_sra_vdb_sys::ReferenceObj = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::PlacementSetIteratorNextReference(
                self.ptr,
                &mut first_pos,
                &mut len,
                &mut ref_obj,
            )
        };
        match check_rc_done(rc)? {
            Some(()) => Ok(Some(NextReference { first_pos, len, ref_obj })),
            None => Ok(None),
        }
    }

    /// Advance to the next window within the current reference.
    ///
    /// Returns `Ok(Some((first_pos, len)))` or `Ok(None)` when windows exhausted.
    pub fn next_window(&mut self) -> Result<Option<(i32, u32)>, VdbError> {
        let mut first_pos: i32 = 0;
        let mut len: u32 = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::PlacementSetIteratorNextWindow(self.ptr, &mut first_pos, &mut len)
        };
        match check_rc_done(rc)? {
            Some(()) => Ok(Some((first_pos, len))),
            None => Ok(None),
        }
    }

    /// Get the next available position within the current window.
    ///
    /// Returns `Ok(Some((pos, len)))` or `Ok(None)` when positions exhausted.
    pub fn next_avail_pos(&self) -> Result<Option<(i32, u32)>, VdbError> {
        let mut pos: i32 = 0;
        let mut len: u32 = 0;
        let rc = unsafe {
            fg_sra_vdb_sys::PlacementSetIteratorNextAvailPos(self.ptr, &mut pos, &mut len)
        };
        match check_rc_done(rc)? {
            Some(()) => Ok(Some((pos, len))),
            None => Ok(None),
        }
    }

    /// Get the next record at the given position.
    ///
    /// Returns `Ok(Some(record))` or `Ok(None)` when no more records at this position.
    pub fn next_record_at(
        &mut self,
        pos: i32,
    ) -> Result<Option<PlacementRecordRef>, VdbError> {
        let mut rec: *const fg_sra_vdb_sys::PlacementRecord = ptr::null();
        let rc = unsafe {
            fg_sra_vdb_sys::PlacementSetIteratorNextRecordAt(self.ptr, pos, &mut rec)
        };
        match check_rc_done(rc)? {
            Some(()) if rec.is_null() => Ok(None),
            Some(()) => Ok(Some(PlacementRecordRef { ptr: rec })),
            None => Ok(None),
        }
    }
}

impl Drop for PlacementSetIterator {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::PlacementSetIteratorRelease(self.ptr) };
        }
    }
}

/// Reference info returned by `PlacementSetIterator::next_reference`.
///
/// The contained `ReferenceObj` is borrowed from the iterator and must not
/// outlive the current reference iteration step.
pub struct NextReference {
    /// First position on this reference (0-based).
    pub first_pos: i32,
    /// Length of the reference span.
    pub len: u32,
    /// Pointer to the reference object (borrowed from the iterator).
    ref_obj: *const fg_sra_vdb_sys::ReferenceObj,
}

impl NextReference {
    /// Get the reference name.
    pub fn ref_name(&self) -> Result<String, VdbError> {
        crate::reference::ref_obj_name(self.ref_obj)
    }

    /// Get the reference sequence ID.
    pub fn seq_id(&self) -> Result<String, VdbError> {
        crate::reference::ref_obj_seq_id(self.ref_obj)
    }
}

/// A reference to a `PlacementRecord` from the iterator.
///
/// This is borrowed from the iterator and should not outlive the current
/// iteration step.
pub struct PlacementRecordRef {
    ptr: *const fg_sra_vdb_sys::PlacementRecord,
}

impl PlacementRecordRef {
    /// Get the alignment row ID.
    pub fn id(&self) -> i64 {
        unsafe { (*self.ptr).id }
    }

    /// Get the reference position (0-based).
    pub fn pos(&self) -> i32 {
        unsafe { (*self.ptr).pos }
    }

    /// Get the alignment length on the reference.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> u32 {
        unsafe { (*self.ptr).len }
    }

    /// Get the mapping quality.
    pub fn mapq(&self) -> i32 {
        unsafe { (*self.ptr).mapq }
    }

    /// Get the spot group (read group) if available.
    pub fn spot_group(&self) -> Option<String> {
        unsafe {
            let rec = &*self.ptr;
            if rec.spot_group.is_null() || rec.spot_group_len == 0 {
                return None;
            }
            let bytes = std::slice::from_raw_parts(
                rec.spot_group as *const u8,
                rec.spot_group_len as usize,
            );
            Some(String::from_utf8_lossy(bytes).into_owned())
        }
    }
}

/// Safe wrapper around a single `PlacementIterator`.
///
/// Created via `ReferenceObj_MakePlacementIterator` and added to a
/// `PlacementSetIterator`.
pub struct PlacementIterator {
    ptr: *mut fg_sra_vdb_sys::PlacementIterator,
}

unsafe impl Send for PlacementIterator {}

impl PlacementIterator {
    /// Create a placement iterator for a reference.
    ///
    /// This is a complex operation that sets up iteration over alignments
    /// on a specific reference using a given cursor.
    pub fn make(
        ref_obj: &ReferenceObj,
        ref_window_start: i32,
        ref_window_len: u32,
        min_mapq: i32,
        align_cursor: &VCursor,
        id_src: AlignIdSrc,
    ) -> Result<Self, VdbError> {
        let mut iter: *mut fg_sra_vdb_sys::PlacementIterator = ptr::null_mut();
        let rc = unsafe {
            fg_sra_vdb_sys::ReferenceObj_MakePlacementIterator(
                ref_obj.as_ptr(),
                &mut iter,
                ref_window_start,
                ref_window_len,
                min_mapq,
                ptr::null(),       // ref_cursor (NULL = use internal)
                align_cursor.as_ptr(),
                id_src as u8,
                ptr::null(),       // ext_0
                ptr::null(),       // ext_1
                ptr::null(),       // rd_group filter
                ptr::null_mut(),   // placement_ctx
            )
        };
        check_rc(rc)?;
        Ok(Self { ptr: iter })
    }

    /// Create from a raw pointer (takes ownership).
    #[allow(dead_code)]
    pub(crate) fn from_raw(ptr: *mut fg_sra_vdb_sys::PlacementIterator) -> Self {
        Self { ptr }
    }
}

impl Drop for PlacementIterator {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { fg_sra_vdb_sys::PlacementIteratorRelease(self.ptr) };
        }
    }
}
