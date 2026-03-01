//! FxHashMap-based mate-pair cache.
//!
//! Stores mate information keyed by alignment ID for resolving mate fields
//! (RNEXT, PNEXT, TLEN) in paired-end reads. Cleared between references
//! in default (non-region) mode.
