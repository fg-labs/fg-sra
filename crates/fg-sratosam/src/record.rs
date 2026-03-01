//! SAM/BAM record building and formatting.
//!
//! Builds complete SAM lines in a per-record `Vec<u8>` buffer, writing
//! all fields in a single pass rather than multiple formatted-print calls.
