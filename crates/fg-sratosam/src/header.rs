//! SAM header generation from VDB metadata.
//!
//! Reads the `BAM_HEADER` metadata node and `ReferenceList` to produce
//! `@HD`, `@SQ`, `@RG`, and `@CO` header lines.
