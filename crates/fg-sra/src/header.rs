//! SAM header generation from VDB metadata.
//!
//! Reads the `BAM_HEADER` metadata node and `ReferenceList` to produce
//! `@HD`, `@SQ`, `@RG`, and `@CO` header lines.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::reference::ReferenceList;

/// Generate the SAM header string for the given database.
///
/// Strategy (matching sam-dump behavior):
/// 1. If `header_file` is provided, read the file as the header.
/// 2. Otherwise, if `regenerate` is true, build from `ReferenceList`.
/// 3. Otherwise, try the stored `BAM_HEADER` metadata node, falling back to
///    `ReferenceList` if not found. A stored header lacking `@SQ` lines for some
///    references (e.g. from a multi-BAM load, where each input BAM's header
///    overwrites the last) gets `@SQ` lines for those references appended.
/// 4. Append any user-supplied `@CO` comment lines.
///
/// The `use_seqid` flag controls whether `@SQ SN:` uses the sequence ID
/// (e.g. `NC_000001.11`) or the reference name (e.g. `chr1`).
pub fn generate_header(
    db: &VDatabase,
    regenerate: bool,
    use_seqid: bool,
    comments: &[String],
    header_file: Option<&Path>,
) -> Result<String> {
    let header = if let Some(path) = header_file {
        std::fs::read_to_string(path)
            .with_context(|| format!("failed to read header file: {}", path.display()))?
    } else if regenerate {
        build_header_from_references(db, use_seqid)?
    } else {
        match read_bam_header(db)? {
            Some(h) => {
                let refs = reference_names_and_lengths(db)?;
                let (h, num_added) = add_missing_sq_lines(&h, &refs, use_seqid);
                if num_added > 0 {
                    eprintln!(
                        "[header] stored BAM_HEADER lacks {num_added} reference(s); \
                         appended @SQ lines for them"
                    );
                }
                h
            }
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
    let Ok(meta) = db.open_metadata_read() else {
        return Ok(None);
    };

    let Ok(node) = meta.open_node_read("BAM_HEADER") else {
        return Ok(None);
    };

    let content = node.read_all().context("failed to read BAM_HEADER metadata node")?;
    if content.is_empty() { Ok(None) } else { Ok(Some(content)) }
}

/// A reference's name, sequence ID and length, from the `ReferenceList`.
struct RefNames {
    name: String,
    seq_id: String,
    length: u32,
}

/// The name, sequence ID and length of every reference in the database.
fn reference_names_and_lengths(db: &VDatabase) -> Result<Vec<RefNames>> {
    let reflist =
        ReferenceList::make_database(db, 0, 0).context("failed to create ReferenceList")?;
    let mut refs = Vec::new();
    for ref_obj in reflist.iter().context("failed to iterate references")? {
        let ref_obj = ref_obj.context("failed to get reference")?;
        refs.push(RefNames {
            name: ref_obj.name().context("failed to get reference name")?,
            seq_id: ref_obj.seq_id().context("failed to get reference seq_id")?,
            length: ref_obj.seq_length().context("failed to get reference length")?,
        });
    }
    Ok(refs)
}

/// The `SN` values of a header's `@SQ` lines.
fn sq_names(header: &str) -> HashSet<&str> {
    header
        .lines()
        .filter(|line| line.starts_with("@SQ\t"))
        .filter_map(|line| line.split('\t').find_map(|field| field.strip_prefix("SN:")))
        .collect()
}

/// Append an `@SQ` line (after the last existing one) for each reference whose
/// name and sequence ID both lack one. `use_seqid` picks the `SN` value, as for
/// a regenerated header. Returns the new header and the number of lines added.
fn add_missing_sq_lines(header: &str, refs: &[RefNames], use_seqid: bool) -> (String, usize) {
    let present = sq_names(header);
    let missing: Vec<String> = refs
        .iter()
        .filter(|r| !present.contains(r.name.as_str()) && !present.contains(r.seq_id.as_str()))
        .map(|r| {
            let sn = if use_seqid { &r.seq_id } else { &r.name };
            format!("@SQ\tSN:{sn}\tLN:{}", r.length)
        })
        .collect();
    if missing.is_empty() {
        return (header.to_owned(), 0);
    }

    let mut lines: Vec<&str> = header.lines().collect();
    // Insert after the last @SQ line, or after the leading @HD line, or first.
    let insert_at = lines
        .iter()
        .rposition(|line| line.starts_with("@SQ\t"))
        .or_else(|| lines.iter().position(|line| line.starts_with("@HD")))
        .map_or(0, |i| i + 1);
    lines.splice(insert_at..insert_at, missing.iter().map(String::as_str));
    let mut result = lines.join("\n");
    result.push('\n');
    (result, missing.len())
}

/// Build the map from reference name to BAM reference ID used for both RNAME
/// and RNEXT: each `@SQ` line's index keyed by its `SN`, plus each reference's
/// other name (name vs sequence ID) keyed to the same index, so a reference
/// resolves whichever form the header uses.
pub fn build_ref_id_map(db: &VDatabase, header: &str) -> Result<HashMap<String, i32>> {
    let mut map = crate::output::build_ref_name_to_id(header);
    add_name_aliases(&mut map, &reference_names_and_lengths(db)?);
    Ok(map)
}

/// For each reference with exactly one of its name/sequence ID in `map`, map the
/// other to the same ID. Never overwrites an existing entry, so an alias cannot
/// shadow another reference's `@SQ` name.
fn add_name_aliases(map: &mut HashMap<String, i32>, refs: &[RefNames]) {
    for r in refs {
        match (map.get(&r.name).copied(), map.get(&r.seq_id).copied()) {
            (Some(id), None) => {
                map.entry(r.seq_id.clone()).or_insert(id);
            }
            (None, Some(id)) => {
                map.entry(r.name.clone()).or_insert(id);
            }
            _ => {}
        }
    }
}

/// Build a minimal SAM header from the `ReferenceList`.
///
/// Produces `@HD VN:1.3` followed by `@SQ SN:name LN:length` for each reference.
fn build_header_from_references(db: &VDatabase, use_seqid: bool) -> Result<String> {
    let reflist =
        ReferenceList::make_database(db, 0, 0).context("failed to create ReferenceList")?;

    let mut header = String::from("@HD\tVN:1.3\n");

    for ref_obj in reflist.iter().context("failed to iterate references")? {
        let ref_obj = ref_obj.context("failed to get reference")?;
        let name = if use_seqid { ref_obj.seq_id() } else { ref_obj.name() }
            .context("failed to get reference name")?;
        let len = ref_obj.seq_length().context("failed to get reference length")?;
        header.push_str("@SQ\tSN:");
        header.push_str(&name);
        header.push_str("\tLN:");
        header.push_str(itoa::Buffer::new().format(len));
        header.push('\n');
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn refs(entries: &[(&str, &str, u32)]) -> Vec<RefNames> {
        entries
            .iter()
            .map(|&(name, seq_id, length)| RefNames {
                name: name.to_string(),
                seq_id: seq_id.to_string(),
                length,
            })
            .collect()
    }

    #[test]
    fn test_add_missing_sq_lines_appends_after_last_sq() {
        let header = "@HD\tVN:1.4\n@SQ\tSN:1\tLN:100\n@RG\tID:a\n";
        let refs = refs(&[("1", "NC_1", 100), ("2", "NC_2", 200), ("3", "NC_3", 300)]);
        let (out, added) = add_missing_sq_lines(header, &refs, false);
        assert_eq!(added, 2);
        assert_eq!(
            out,
            "@HD\tVN:1.4\n@SQ\tSN:1\tLN:100\n@SQ\tSN:2\tLN:200\n@SQ\tSN:3\tLN:300\n@RG\tID:a\n"
        );
    }

    #[test]
    fn test_add_missing_sq_lines_matches_either_name_form() {
        // "NC_1" is present by seq_id, so only reference 2 is added, by seq_id.
        let header = "@HD\tVN:1.4\n@SQ\tSN:NC_1\tLN:100\n";
        let refs = refs(&[("1", "NC_1", 100), ("2", "NC_2", 200)]);
        let (out, added) = add_missing_sq_lines(header, &refs, true);
        assert_eq!(added, 1);
        assert_eq!(out, "@HD\tVN:1.4\n@SQ\tSN:NC_1\tLN:100\n@SQ\tSN:NC_2\tLN:200\n");
    }

    #[test]
    fn test_add_missing_sq_lines_complete_header_unchanged() {
        let header = "@HD\tVN:1.4\n@SQ\tSN:1\tLN:100\n";
        let (out, added) = add_missing_sq_lines(header, &refs(&[("1", "NC_1", 100)]), false);
        assert_eq!((out.as_str(), added), (header, 0));
    }

    #[test]
    fn test_add_missing_sq_lines_without_sq_goes_after_hd() {
        let header = "@HD\tVN:1.4\n@RG\tID:a\n";
        let (out, _) = add_missing_sq_lines(header, &refs(&[("1", "NC_1", 100)]), false);
        assert_eq!(out, "@HD\tVN:1.4\n@SQ\tSN:1\tLN:100\n@RG\tID:a\n");
    }

    #[test]
    fn test_add_name_aliases_maps_other_form_to_same_id() {
        let mut map = HashMap::from([("NC_1".to_string(), 0), ("chr2".to_string(), 1)]);
        add_name_aliases(&mut map, &refs(&[("chr1", "NC_1", 1), ("chr2", "NC_2", 1)]));
        assert_eq!(map.get("chr1"), Some(&0));
        assert_eq!(map.get("NC_2"), Some(&1));
        assert_eq!(map.len(), 4);
    }

    #[test]
    fn test_add_name_aliases_never_shadows_existing_name() {
        // Reference A's seq_id equals reference B's @SQ name; B keeps its ID.
        let mut map = HashMap::from([("a".to_string(), 0), ("b".to_string(), 1)]);
        add_name_aliases(&mut map, &refs(&[("a", "b", 1)]));
        assert_eq!(map.get("b"), Some(&1));
    }

    #[test]
    fn test_header_file_is_read() {
        let mut tmp = std::env::temp_dir();
        tmp.push("fg_sra_test_header.sam");
        {
            let mut f = std::fs::File::create(&tmp).unwrap();
            f.write_all(b"@HD\tVN:1.6\n@SQ\tSN:custom\tLN:100\n").unwrap();
        }

        // Exercise the header-file branch logic without a live VDB database.
        let header = std::fs::read_to_string(&tmp).unwrap();
        let mut result = header;
        if !result.is_empty() && !result.ends_with('\n') {
            result.push('\n');
        }
        result.push_str("@CO\textra\n");

        assert!(result.starts_with("@HD\tVN:1.6\n"));
        assert!(result.contains("@SQ\tSN:custom\tLN:100\n"));
        assert!(result.contains("@CO\textra\n"));

        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_header_file_missing_returns_error() {
        let result = std::fs::read_to_string("/nonexistent/header.sam");
        assert!(result.is_err());
    }

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
