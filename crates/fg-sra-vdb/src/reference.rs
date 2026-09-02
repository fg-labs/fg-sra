//! Safe wrappers for `ReferenceList` and `ReferenceObj`.
//!
//! Provides access to reference sequence metadata needed for SAM header
//! generation and coordinate-sorted alignment iteration.

use std::ptr;

use crate::database::VDatabase;
use crate::error::{VdbError, check_rc, to_cstring};
use crate::retry::retry_on_network_error;

/// Option flags for [`ReferenceList::make_database`].
///
/// These control which alignment ID columns are read from the REFERENCE table.
pub mod reflist_options {
    /// Make `ReferenceObj::read*` return `INSDC:4na:bin` (4-bit codes 0..=15)
    /// instead of DNA text. Used so callers can apply the CHARSET mapping
    /// themselves and match the `(ascii)READ` column exactly.
    pub const READ_4NA: u32 = 0x01;
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
        retry_on_network_error("ReferenceList_MakeDatabase", || {
            let mut reflist: *const fg_sra_vdb_sys::ReferenceList = ptr::null();
            let rc = unsafe {
                fg_sra_vdb_sys::ReferenceList_MakeDatabase(
                    &raw mut reflist,
                    db.as_ptr(),
                    options,
                    cache,
                    ptr::null(), // name filter
                    0,           // numbins
                )
            };
            check_rc(rc)?;
            Ok(Self { ptr: reflist })
        })
    }

    /// Get the number of references.
    pub fn count(&self) -> Result<u32, VdbError> {
        retry_on_network_error("ReferenceList_Count", || {
            let mut count: u32 = 0;
            let rc = unsafe { fg_sra_vdb_sys::ReferenceList_Count(self.ptr, &raw mut count) };
            check_rc(rc)?;
            Ok(count)
        })
    }

    /// Get a reference by index.
    pub fn get(&self, idx: u32) -> Result<ReferenceObj, VdbError> {
        retry_on_network_error("ReferenceList_Get", || {
            let mut obj: *const fg_sra_vdb_sys::ReferenceObj = ptr::null();
            let rc = unsafe { fg_sra_vdb_sys::ReferenceList_Get(self.ptr, &raw mut obj, idx) };
            check_rc(rc)?;
            Ok(ReferenceObj { ptr: obj })
        })
    }

    /// Find a reference by name.
    pub fn find(&self, name: &str) -> Result<ReferenceObj, VdbError> {
        let c_name = to_cstring(name)?;
        retry_on_network_error("ReferenceList_Find", || {
            let mut obj: *const fg_sra_vdb_sys::ReferenceObj = ptr::null();
            let rc = unsafe {
                fg_sra_vdb_sys::ReferenceList_Find(
                    self.ptr,
                    &raw mut obj,
                    c_name.as_ptr(),
                    c_name.as_bytes().len(),
                )
            };
            check_rc(rc)?;
            Ok(ReferenceObj { ptr: obj })
        })
    }

    /// Iterate over all references.
    #[allow(clippy::iter_not_returning_iterator)]
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

impl Iterator for ReferenceIter<'_> {
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

    /// Get the sequence ID (e.g., "`NC_000001.11`").
    pub fn seq_id(&self) -> Result<String, VdbError> {
        ref_obj_seq_id(self.ptr)
    }

    /// Get the sequence length.
    pub fn seq_length(&self) -> Result<u32, VdbError> {
        retry_on_network_error("ReferenceObj_SeqLength", || {
            let mut len: u32 = 0;
            let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_SeqLength(self.ptr, &raw mut len) };
            check_rc(rc)?;
            Ok(len)
        })
    }

    /// Get the reference index.
    pub fn idx(&self) -> Result<u32, VdbError> {
        retry_on_network_error("ReferenceObj_Idx", || {
            let mut idx: u32 = 0;
            let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_Idx(self.ptr, &raw mut idx) };
            check_rc(rc)?;
            Ok(idx)
        })
    }

    /// Read up to `buf.len()` reference bases starting at 0-based `offset`,
    /// returning the number of bases written into `buf`.
    ///
    /// With [`reflist_options::READ_4NA`] the bytes are 4na codes (`0..=15`);
    /// otherwise they are DNA text. A read past the end of a non-circular
    /// reference writes fewer bytes than requested.
    ///
    /// # Thread safety
    /// This uses the owning [`ReferenceList`]'s single shared reader. Never call
    /// it concurrently with any other read on the same list; the types are
    /// `Send` but not `Sync`, so the compiler enforces single-threaded use.
    pub fn read_into(&self, offset: i32, buf: &mut [u8]) -> Result<u32, VdbError> {
        let len = u32::try_from(buf.len()).unwrap_or(u32::MAX);
        retry_on_network_error("ReferenceObj_Read", || {
            let mut written: u32 = 0;
            let rc = unsafe {
                fg_sra_vdb_sys::ReferenceObj_Read(
                    self.ptr,
                    offset,
                    len,
                    buf.as_mut_ptr(),
                    &raw mut written,
                )
            };
            check_rc(rc)?;
            Ok(written)
        })
    }

    /// Read the entire reference sequence into a freshly allocated `Vec`.
    ///
    /// Reads sequentially in chunks until the sequence length is reached (or a
    /// chunk returns zero bases). The result may be shorter than
    /// [`seq_length`](Self::seq_length) only if the underlying reader stops
    /// early; callers that require the full sequence should check the length.
    pub fn read_all(&self) -> Result<Vec<u8>, VdbError> {
        let len = self.seq_length()? as usize;
        let mut out = vec![0u8; len];
        let mut filled = 0usize;
        let mut offset: i32 = 0;
        while filled < len {
            let written = self.read_into(offset, &mut out[filled..])? as usize;
            if written == 0 {
                break;
            }
            filled += written;
            offset = i32::try_from(filled).unwrap_or(i32::MAX);
        }
        out.truncate(filled);
        Ok(out)
    }

    /// Whether this reference is circular.
    pub fn circular(&self) -> Result<bool, VdbError> {
        retry_on_network_error("ReferenceObj_Circular", || {
            let mut circular = false;
            let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_Circular(self.ptr, &raw mut circular) };
            check_rc(rc)?;
            Ok(circular)
        })
    }

    /// Whether this reference is stored externally (a separate refseq object)
    /// rather than locally inside the run.
    pub fn external(&self) -> Result<bool, VdbError> {
        retry_on_network_error("ReferenceObj_External", || {
            let mut external = false;
            let rc = unsafe {
                fg_sra_vdb_sys::ReferenceObj_External(self.ptr, &raw mut external, ptr::null_mut())
            };
            check_rc(rc)?;
            Ok(external)
        })
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
    retry_on_network_error("ReferenceObj_Name", || {
        let mut name: *const std::os::raw::c_char = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_Name(ptr, &raw mut name) };
        check_rc(rc)?;
        if name.is_null() {
            return Ok(String::new());
        }
        let s = unsafe { std::ffi::CStr::from_ptr(name) };
        Ok(s.to_string_lossy().into_owned())
    })
}

/// Get the sequence ID from a raw `ReferenceObj` pointer.
///
/// Shared by `ReferenceObj::seq_id()` and `NextReference::seq_id()`.
pub(crate) fn ref_obj_seq_id(ptr: *const fg_sra_vdb_sys::ReferenceObj) -> Result<String, VdbError> {
    retry_on_network_error("ReferenceObj_SeqId", || {
        let mut seqid: *const std::os::raw::c_char = ptr::null();
        let rc = unsafe { fg_sra_vdb_sys::ReferenceObj_SeqId(ptr, &raw mut seqid) };
        check_rc(rc)?;
        if seqid.is_null() {
            return Ok(String::new());
        }
        let s = unsafe { std::ffi::CStr::from_ptr(seqid) };
        Ok(s.to_string_lossy().into_owned())
    })
}
