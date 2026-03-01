//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip/bzip2 compression
//! and multi-threaded BGZF for BAM output.
