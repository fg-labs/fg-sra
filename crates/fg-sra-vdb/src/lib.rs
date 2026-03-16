//! Safe Rust wrappers over the NCBI VDB C library.
//!
//! Provides RAII types with `Drop` implementations, `Result`-based error handling,
//! and typed column reads for working with SRA/VDB databases.

pub mod cursor;
pub mod database;
pub mod error;
pub mod iterator;
pub mod manager;
pub mod reference;
pub mod retry;
