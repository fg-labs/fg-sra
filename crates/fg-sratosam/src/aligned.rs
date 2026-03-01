//! Aligned read processing via PlacementSetIterator.
//!
//! Walks references → windows → positions → records, building SAM/BAM
//! output for each aligned read. Supports parallel processing across references.
