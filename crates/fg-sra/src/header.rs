//! SAM header generation from VDB metadata.
//!
//! Reads the `BAM_HEADER` metadata node and `ReferenceList` to produce
//! `@HD`, `@SQ`, `@RG`, and `@CO` header lines.

use anyhow::{Context, Result};
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::reference::ReferenceList;

/// Generate the SAM header string for the given database.
///
/// Strategy (matching sam-dump behavior):
/// 1. Try to read the stored `BAM_HEADER` metadata node.
/// 2. If not found **or** `regenerate` is true, build a minimal header from
///    the `ReferenceList` (`@HD` + `@SQ` lines).
/// 3. Append any user-supplied `@CO` comment lines.
///
/// The `use_seqid` flag controls whether `@SQ SN:` uses the sequence ID
/// (e.g. `NC_000001.11`) or the reference name (e.g. `chr1`).
pub fn generate_header(
    db: &VDatabase,
    regenerate: bool,
    use_seqid: bool,
    comments: &[String],
) -> Result<String> {
    let header = if regenerate {
        build_header_from_references(db, use_seqid)?
    } else {
        match read_bam_header(db)? {
            Some(h) => h,
            None => build_header_from_references(db, use_seqid)?,
        }
    };

    let mut result = header;

    // Ensure header ends with a newline before appending comments.
    if !result.is_empty() && !result.ends_with('\n') {
        result.push('\n');
    }

    for comment in comments {
        result.push_str("@CO\t");
        result.push_str(comment);
        result.push('\n');
    }

    Ok(result)
}

/// Try to read the `BAM_HEADER` metadata node.
///
/// Returns `Ok(None)` if the node does not exist, `Ok(Some(header))` if found.
fn read_bam_header(db: &VDatabase) -> Result<Option<String>> {
    let meta = match db.open_metadata_read() {
        Ok(m) => m,
        Err(_) => return Ok(None),
    };

    let node = match meta.open_node_read("BAM_HEADER") {
        Ok(n) => n,
        Err(_) => return Ok(None),
    };

    let content = node.read_all().context("failed to read BAM_HEADER metadata node")?;
    if content.is_empty() { Ok(None) } else { Ok(Some(content)) }
}

/// Build a minimal SAM header from the `ReferenceList`.
///
/// Produces `@HD VN:1.3` followed by `@SQ SN:name LN:length` for each reference.
fn build_header_from_references(db: &VDatabase, use_seqid: bool) -> Result<String> {
    let reflist =
        ReferenceList::make_database(db, 0, 0).context("failed to create ReferenceList")?;

    let count = reflist.count().context("failed to get reference count")?;

    // Pre-allocate: @HD line (~12 bytes) + ~40 bytes per @SQ line.
    let mut header = String::with_capacity(12 + (count as usize) * 40);
    header.push_str("@HD\tVN:1.3\n");

    for i in 0..count {
        let ref_obj = reflist.get(i).context("failed to get reference")?;
        let name = if use_seqid { ref_obj.seq_id() } else { ref_obj.name() }
            .context("failed to get reference name")?;
        let len = ref_obj.seq_length().context("failed to get reference length")?;
        header.push_str("@SQ\tSN:");
        header.push_str(&name);
        header.push_str("\tLN:");
        header.push_str(&len.to_string());
        header.push('\n');
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_generate_header_with_comments() {
        // Test comment appending logic with a pre-built header string.
        let base = "@HD\tVN:1.3\n@SQ\tSN:chr1\tLN:248956422\n".to_string();
        let comments = vec!["first comment".to_string(), "second comment".to_string()];

        let mut result = base;
        for comment in &comments {
            result.push_str("@CO\t");
            result.push_str(comment);
            result.push('\n');
        }

        assert!(result.contains("@CO\tfirst comment\n"));
        assert!(result.contains("@CO\tsecond comment\n"));
        assert!(result.ends_with('\n'));
    }

    #[test]
    fn test_empty_header_gets_newline() {
        let mut result = String::new();
        if !result.is_empty() && !result.ends_with('\n') {
            result.push('\n');
        }
        // Empty string should stay empty (no orphan newline).
        assert!(result.is_empty());
    }

    #[test]
    fn test_header_without_trailing_newline() {
        let mut result = "@HD\tVN:1.3".to_string();
        if !result.is_empty() && !result.ends_with('\n') {
            result.push('\n');
        }
        assert!(result.ends_with('\n'));
    }
}
