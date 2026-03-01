//! Safe Rust wrappers over the NCBI VDB C library.
//!
//! Provides RAII types with `Drop` implementations, `Result`-based error handling,
//! and typed column reads for working with SRA/VDB databases.

pub mod error;
pub mod manager;
pub mod database;
pub mod cursor;
pub mod reference;
pub mod iterator;
