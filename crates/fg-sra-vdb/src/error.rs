//! VDB error handling: `rc_t` to `Result` conversion.
//!
//! The VDB C API returns `rc_t` (a 32-bit packed error code) from most functions.
//! This module provides conversion to idiomatic Rust `Result` types.

use std::fmt;

/// A VDB error wrapping the C library's `rc_t` return code.
///
/// The `rc_t` encodes five bit fields:
///   - Module  (bits 27-31, 5 bits)
///   - Target  (bits 21-26, 6 bits)
///   - Context (bits 14-20, 7 bits)
///   - Object  (bits  6-13, 8 bits)
///   - State   (bits  0-5,  6 bits)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VdbError {
    /// An error returned by the VDB C API via `rc_t`.
    Rc(u32),
    /// A string argument contained an interior nul byte.
    InvalidNulByte,
}

/// The `rcDone` state value, indicating iteration is complete (not an error).
const RC_STATE_DONE: u32 = 1;

impl VdbError {
    /// Create a new `VdbError` from a raw `rc_t` value.
    pub fn new(rc: u32) -> Self {
        Self::Rc(rc)
    }

    /// Returns `true` if this error represents the "done" state (end of iteration).
    pub fn is_done(&self) -> bool {
        matches!(self, Self::Rc(rc) if rc & 0x3F == RC_STATE_DONE)
    }

    /// Extract the raw `rc_t` value, or `None` for non-rc errors.
    pub fn rc(&self) -> Option<u32> {
        match self {
            Self::Rc(rc) => Some(*rc),
            Self::InvalidNulByte => None,
        }
    }

    /// Extract the State field (bits 0-5).
    pub fn state(&self) -> u32 {
        match self {
            Self::Rc(rc) => rc & 0x3F,
            _ => 0,
        }
    }

    /// Extract the Object field (bits 6-13).
    pub fn object(&self) -> u32 {
        match self {
            Self::Rc(rc) => (rc >> 6) & 0xFF,
            _ => 0,
        }
    }

    /// Extract the Context field (bits 14-20).
    pub fn context(&self) -> u32 {
        match self {
            Self::Rc(rc) => (rc >> 14) & 0x7F,
            _ => 0,
        }
    }

    /// Extract the Target field (bits 21-26).
    pub fn target(&self) -> u32 {
        match self {
            Self::Rc(rc) => (rc >> 21) & 0x3F,
            _ => 0,
        }
    }

    /// Extract the Module field (bits 27-31).
    pub fn module(&self) -> u32 {
        match self {
            Self::Rc(rc) => (rc >> 27) & 0x1F,
            _ => 0,
        }
    }
}

impl fmt::Display for VdbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Rc(rc) => write!(
                f,
                "VDB error (rc=0x{rc:08x}, mod={}, tgt={}, ctx={}, obj={}, state={})",
                self.module(),
                self.target(),
                self.context(),
                self.object(),
                self.state()
            ),
            Self::InvalidNulByte => write!(f, "VDB error: string argument contains interior nul byte"),
        }
    }
}

impl std::error::Error for VdbError {}

/// Convert a VDB `rc_t` return value to a `Result`.
///
/// Returns `Ok(())` if `rc == 0`, or `Err(VdbError)` otherwise.
pub fn check_rc(rc: u32) -> Result<(), VdbError> {
    if rc == 0 { Ok(()) } else { Err(VdbError::new(rc)) }
}

/// Convert a `&str` to a `CString`, returning `VdbError::InvalidNulByte` on failure.
pub fn to_cstring(s: &str) -> Result<std::ffi::CString, VdbError> {
    std::ffi::CString::new(s).map_err(|_| VdbError::InvalidNulByte)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_rc_ok() {
        assert!(check_rc(0).is_ok());
    }

    #[test]
    fn test_check_rc_err() {
        assert!(check_rc(1).is_err());
    }

    #[test]
    fn test_error_is_done() {
        let err = VdbError::new(1); // state=1 is rcDone
        assert!(err.is_done());
    }

    #[test]
    fn test_error_fields() {
        // Construct an rc_t: module=10, target=12, context=7, object=3, state=1
        let rc = (10 << 27) | (12 << 21) | (7 << 14) | (3 << 6) | 1;
        let err = VdbError::new(rc);
        assert_eq!(err.module(), 10);
        assert_eq!(err.target(), 12);
        assert_eq!(err.context(), 7);
        assert_eq!(err.object(), 3);
        assert_eq!(err.state(), 1);
        assert!(err.is_done());
    }

    #[test]
    fn test_invalid_nul_byte() {
        let err = VdbError::InvalidNulByte;
        assert!(!err.is_done());
        assert_eq!(err.rc(), None);
        assert!(err.to_string().contains("nul byte"));
    }

    #[test]
    fn test_to_cstring() {
        assert!(to_cstring("hello").is_ok());
        assert_eq!(to_cstring("hello\0world").unwrap_err(), VdbError::InvalidNulByte);
    }
}
