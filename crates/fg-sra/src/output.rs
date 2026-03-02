//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip/bzip2 compression
//! and multi-threaded BGZF for BAM output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::Result;

/// Default buffer size for the output writer (256 KB).
const OUTPUT_BUF_SIZE: usize = 256 * 1024;

/// Abstraction over output destinations (stdout or file, with optional compression).
///
/// Phase 2 supports SAM text output only. BAM and compression (gzip/bzip2)
/// will be added in Phase 4.
pub struct OutputWriter {
    inner: BufWriter<Box<dyn Write>>,
}

impl OutputWriter {
    /// Create a writer to stdout.
    pub fn stdout() -> Self {
        Self { inner: BufWriter::with_capacity(OUTPUT_BUF_SIZE, Box::new(io::stdout().lock())) }
    }

    /// Create a writer to a file.
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::create(path)?;
        Ok(Self { inner: BufWriter::with_capacity(OUTPUT_BUF_SIZE, Box::new(file)) })
    }

    /// Write a complete SAM header (already formatted with newlines).
    pub fn write_header(&mut self, header: &str) -> Result<()> {
        self.inner.write_all(header.as_bytes())?;
        Ok(())
    }

    /// Write a pre-formatted SAM record line (must include trailing newline).
    pub fn write_record(&mut self, record: &[u8]) -> Result<()> {
        self.inner.write_all(record)?;
        Ok(())
    }

    /// Flush and finalize the output.
    pub fn finish(mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}
