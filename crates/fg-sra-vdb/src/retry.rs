//! Retry-with-backoff for transient VDB network errors.
//!
//! VDB makes implicit network calls (accession resolution, reference metadata,
//! lazy page fetches). When connectivity is flaky, these fail with `module == rcNS`.
//! This module provides a retry helper that transparently retries on network errors.

use std::thread;
use std::time::Duration;

use crate::error::VdbError;

/// Maximum number of retry attempts after the initial call.
const MAX_RETRIES: u32 = 3;

/// Base backoff in milliseconds. Multiplied by 4^attempt (100ms, 400ms, 1600ms).
const BASE_BACKOFF_MS: u64 = 100;

/// Retry `op` on transient VDB network errors with exponential backoff.
///
/// Calls `op` up to `MAX_RETRIES + 1` times. On success or a non-network error,
/// returns immediately. On a network error, logs a warning and sleeps before retrying.
///
/// `context` is included in the warning message to identify which operation is retrying.
#[inline]
pub fn retry_on_network_error<T>(
    context: &str,
    mut op: impl FnMut() -> Result<T, VdbError>,
) -> Result<T, VdbError> {
    let mut last_err;
    match op() {
        Ok(val) => return Ok(val),
        Err(e) if !e.is_network_error() => return Err(e),
        Err(e) => last_err = e,
    }

    for attempt in 0..MAX_RETRIES {
        let backoff_ms = BASE_BACKOFF_MS * 4u64.pow(attempt);
        let rc_hex = last_err.rc().unwrap_or(0);
        eprintln!(
            "[fg-sra-vdb] {context}: network error (rc=0x{rc_hex:08x}), \
             retry {}/{MAX_RETRIES} after {backoff_ms}ms",
            attempt + 1,
        );
        thread::sleep(Duration::from_millis(backoff_ms));

        match op() {
            Ok(val) => return Ok(val),
            Err(e) if !e.is_network_error() => return Err(e),
            Err(e) => last_err = e,
        }
    }

    Err(last_err)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::*;
    use crate::error::VdbError;

    /// Build a network error with module=18.
    fn network_error() -> VdbError {
        VdbError::new(0x900995d8)
    }

    /// Build a non-network error with module=10.
    fn non_network_error() -> VdbError {
        VdbError::new((10 << 27) | 1)
    }

    #[test]
    fn test_succeeds_on_first_try() {
        let calls = Cell::new(0u32);
        let result = retry_on_network_error("test", || {
            calls.set(calls.get() + 1);
            Ok::<_, VdbError>(42)
        });
        assert_eq!(result.unwrap(), 42);
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn test_succeeds_after_two_failures() {
        let calls = Cell::new(0u32);
        let result = retry_on_network_error("test", || {
            calls.set(calls.get() + 1);
            if calls.get() <= 2 { Err(network_error()) } else { Ok(99) }
        });
        assert_eq!(result.unwrap(), 99);
        assert_eq!(calls.get(), 3);
    }

    #[test]
    fn test_exhausts_all_retries() {
        let calls = Cell::new(0u32);
        let result: Result<i32, VdbError> = retry_on_network_error("test", || {
            calls.set(calls.get() + 1);
            Err(network_error())
        });
        assert!(result.is_err());
        assert!(result.unwrap_err().is_network_error());
        // 1 initial + 3 retries = 4 total calls
        assert_eq!(calls.get(), 4);
    }

    #[test]
    fn test_non_network_error_not_retried() {
        let calls = Cell::new(0u32);
        let result: Result<i32, VdbError> = retry_on_network_error("test", || {
            calls.set(calls.get() + 1);
            Err(non_network_error())
        });
        assert!(result.is_err());
        assert!(!result.unwrap_err().is_network_error());
        assert_eq!(calls.get(), 1);
    }
}
