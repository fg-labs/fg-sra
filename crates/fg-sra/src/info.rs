//! `fg-sra info`: describe SRA archives, to decide how to convert them.

use std::collections::BTreeMap;
use std::fmt::{self, Write as _};

use anyhow::{Context, Result, bail};
use clap::Parser;
use fg_sra_vdb::database::{VDatabase, VTable};
use fg_sra_vdb::manager::disable_remote_access;
use fg_sra_vdb::reference::ReferenceList;

use crate::archive::{Archive, QualitySource, default_accession};
use crate::fastq::spot::ReadFilter;
use crate::progress::format_count;
use crate::record::READ_TYPE_BIOLOGICAL;

/// Spot groups listed by name before the rest are only counted.
const MAX_SPOT_GROUPS_LISTED: usize = 10;

/// Describe SRA archives: what they are, what they hold, and how their spots are laid out.
///
/// For each archive, prints its kind, platform, loader, stored totals, whether qualities and
/// original names were kept, its spot groups, the alignments and references of an aligned
/// archive, and the read layout of its first spots: each read slot's type and lengths, and how
/// many biological and technical reads spots have. The layout shows which `fg-sra fastq`
/// outputs an archive needs, and which technical read is which.
#[derive(Debug, Parser)]
pub struct Info {
    /// SRA archive(s): local `.sra`/`.sralite` files, accessions or https URLs.
    #[arg(required = true)]
    pub inputs: Vec<String>,

    /// Spots to sample, from the first, for the read layout.
    #[arg(long, default_value_t = 10_000, value_name = "N")]
    pub layout_spots: u64,

    /// Never use the network: find archives only locally.
    #[arg(long)]
    pub offline: bool,
}

impl Info {
    /// Describe each input on stdout; a failure is reported and the rest still described.
    pub fn execute(&self) -> Result<()> {
        if self.offline {
            disable_remote_access().context("failed to turn off remote access")?;
        }
        let mut num_failed = 0;
        for (i, input) in self.inputs.iter().enumerate() {
            if i > 0 {
                println!();
            }
            match describe(input, self.layout_spots) {
                Ok(report) => print!("{report}"),
                Err(e) => {
                    eprintln!("[info] {input}: ERROR: {e:#}");
                    num_failed += 1;
                }
            }
        }
        if num_failed > 0 {
            bail!("{num_failed} of {} archive(s) could not be described", self.inputs.len());
        }
        Ok(())
    }
}

/// Everything reported about one archive.
struct Report<'a> {
    accession: String,
    archive: &'a Archive,
    layout: ReadLayout,
    spot_groups: Vec<String>,
    alignments: Option<AlignmentSummary>,
}

impl Report<'_> {
    /// Render a report as indented `label  value` lines.
    fn render(&self) -> String {
        let archive = self.archive;
        let mut out = String::new();
        let mut line = |label: &str, value: &str| {
            let _ = writeln!(out, "  {label:<18}{value}");
        };
        let count = format_count;

        line("location", &archive.location);
        line("kind", &kind(archive));
        if !archive.tables.is_empty() {
            line("tables", &archive.tables.join(", "));
        }
        line("platform", archive.platform.unwrap_or("unknown"));
        line("loader", archive.loader.as_deref().unwrap_or("unknown"));
        let spots = count(archive.spot_count);
        if archive.spot_count == 0 {
            line("spots", &spots);
        } else {
            let (first, last) = (archive.first_spot, archive.last_spot());
            line("spots", &format!("{spots} (spot ids {first}-{last})"));
        }
        let totals = archive.totals;
        let stored = [
            totals.spots.map(|n| format!("{} spots", count(n))),
            totals.bases.map(|n| format!("{} bases", count(n))),
            totals.biological_bases.map(|n| format!("{} biological bases", count(n))),
        ];
        let stored: Vec<String> = stored.into_iter().flatten().collect();
        line(
            "stored totals",
            &if stored.is_empty() { "none".to_owned() } else { stored.join(", ") },
        );
        line("qualities", &qualities(archive));
        line("original names", if archive.has_names { "kept" } else { "not kept" });
        line("spot groups", &spot_group_summary(&self.spot_groups));
        if let Some(alignments) = self.alignments {
            line("alignments", &alignments.to_string());
        }

        let layout = &self.layout;
        let sampled = if layout.spots == archive.spot_count {
            format!("all {} spots", count(layout.spots))
        } else {
            format!(
                "first {} of {} spots only (--layout-spots samples more)",
                count(layout.spots),
                count(archive.spot_count)
            )
        };
        line("read layout", &sampled);
        let _ = writeln!(
            out,
            "    {:<6}{:<12}{:>12}{:>10}{:>8}{:>9}{:>8}",
            "slot", "type", "non-empty", "empty", "min", "mean", "max"
        );
        for (slot, stats) in layout.slots.iter().enumerate() {
            let non_empty = stats.spots - stats.empty;
            let (min, mean, max) = if non_empty == 0 {
                ("-".to_owned(), "-".to_owned(), "-".to_owned())
            } else {
                let mean = stats.total_len as f64 / non_empty as f64;
                (stats.min_len.to_string(), format!("{mean:.1}"), stats.max_len.to_string())
            };
            let _ = writeln!(
                out,
                "    {:<6}{:<12}{:>12}{:>10}{:>8}{:>9}{:>8}",
                slot + 1,
                stats.read_type(),
                count(non_empty),
                count(stats.empty),
                min,
                mean,
                max
            );
        }
        let mut line = |label: &str, value: &str| {
            let _ = writeln!(out, "  {label:<18}{value}");
        };
        line("reads per spot", &distribution(&layout.spots_by_read_count, |_| ""));
        line(
            "biological reads",
            &distribution(&layout.spots_by_bio_reads, |reads| match reads {
                1 => " (unpaired)",
                2 => " (pairs)",
                _ => " (dropped)",
            }),
        );
        line("technical reads", &distribution(&layout.spots_by_technical_reads, |_| ""));
        let filters: Vec<String> = ReadFilter::ALL
            .iter()
            .zip(layout.reads_by_filter)
            .filter(|&(_, reads)| reads > 0)
            .map(|(filter, reads)| format!("{} {}", filter.name(), count(reads)))
            .collect();
        line("read filters", &if filters.is_empty() { "-".to_owned() } else { filters.join(", ") });

        format!("{}\n{out}", self.accession)
    }
}

/// The read layout of a sample of spots.
#[derive(Debug, Default, PartialEq, Eq)]
struct ReadLayout {
    spots: u64,
    /// Statistics of each read slot, in slot order.
    slots: Vec<SlotStats>,
    /// Spots by number of read slots.
    spots_by_read_count: BTreeMap<usize, u64>,
    /// Spots by number of non-empty biological reads.
    spots_by_bio_reads: BTreeMap<usize, u64>,
    /// Spots by number of non-empty technical reads.
    spots_by_technical_reads: BTreeMap<usize, u64>,
    /// Non-empty reads by `READ_FILTER`.
    reads_by_filter: [u64; 4],
}

impl ReadLayout {
    /// Sample spots `first..=last` of `table`.
    fn sample(table: &VTable, first: i64, last: i64) -> Result<Self> {
        let cursor = table.create_cursor_read()?;
        let read_lens = cursor.add_column("(INSDC:coord:len)READ_LEN")?;
        let read_types = cursor.add_column("(INSDC:SRA:xread_type)READ_TYPE")?;
        let read_filters = cursor.add_column("(INSDC:SRA:read_filter)READ_FILTER")?;
        cursor.open().context("failed to open a cursor on the reads table")?;
        let (mut lens, mut types, mut filters) = (Vec::new(), Vec::new(), Vec::new());
        let mut layout = Self::default();
        for id in first..=last {
            cursor.read_u32_slice_into(id, read_lens, &mut lens)?;
            cursor.read_u8_slice_into(id, read_types, &mut types)?;
            cursor.read_u8_slice_into(id, read_filters, &mut filters)?;
            layout.add(&lens, &types, &filters);
        }
        Ok(layout)
    }

    /// Add one spot, given its reads' lengths, types and filters.
    fn add(&mut self, read_lens: &[u32], read_types: &[u8], read_filters: &[u8]) {
        self.spots += 1;
        *self.spots_by_read_count.entry(read_lens.len()).or_default() += 1;
        if self.slots.len() < read_lens.len() {
            self.slots.resize(read_lens.len(), SlotStats::new());
        }
        let (mut bio_reads, mut technical_reads) = (0, 0);
        for (slot, &len) in read_lens.iter().enumerate() {
            let biological = read_types
                .get(slot)
                .is_some_and(|&read_type| read_type & READ_TYPE_BIOLOGICAL != 0);
            let stats = &mut self.slots[slot];
            stats.spots += 1;
            stats.biological += u64::from(biological);
            if len == 0 {
                stats.empty += 1;
                continue;
            }
            stats.min_len = stats.min_len.min(len);
            stats.max_len = stats.max_len.max(len);
            stats.total_len += u64::from(len);
            if biological {
                bio_reads += 1;
            } else {
                technical_reads += 1;
            }
            if let Some(filter) = read_filters.get(slot).and_then(|&raw| ReadFilter::from_raw(raw))
            {
                self.reads_by_filter[filter.index()] += 1;
            }
        }
        *self.spots_by_bio_reads.entry(bio_reads).or_default() += 1;
        *self.spots_by_technical_reads.entry(technical_reads).or_default() += 1;
    }
}

/// Statistics of one read slot over the spots that have it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SlotStats {
    spots: u64,
    biological: u64,
    empty: u64,
    /// Shortest non-empty read; `u32::MAX` while there are none.
    min_len: u32,
    max_len: u32,
    total_len: u64,
}

impl SlotStats {
    /// Statistics of a slot no spot has yet.
    fn new() -> Self {
        Self { spots: 0, biological: 0, empty: 0, min_len: u32::MAX, max_len: 0, total_len: 0 }
    }

    /// The slot's read type: `biological` or `technical`, or `mixed` if spots disagree.
    fn read_type(&self) -> &'static str {
        match self.biological {
            0 => "technical",
            n if n == self.spots => "biological",
            _ => "mixed",
        }
    }
}

/// The alignments and references of an aligned archive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AlignmentSummary {
    primary: u64,
    secondary: Option<u64>,
    references: u32,
    /// References whose bases are in other archives, rather than in this one.
    external_references: u32,
}

impl AlignmentSummary {
    /// Count the alignments and references of `database`.
    fn read(database: &VDatabase) -> Result<Self> {
        let rows = |name: &str| -> Result<u64> {
            let table = database.open_table_read(name)?;
            let cursor = table.create_cursor_read()?;
            let column = cursor.add_column("SEQ_SPOT_ID")?;
            cursor.open()?;
            Ok(cursor.id_range(column)?.1)
        };
        let secondary = if database.has_table("SECONDARY_ALIGNMENT") {
            Some(rows("SECONDARY_ALIGNMENT")?)
        } else {
            None
        };
        let references = ReferenceList::make_database(database, 0, 0)
            .context("failed to list the archive's references")?;
        let count = references.count()?;
        let mut external_references = 0;
        for idx in 0..count {
            external_references += u32::from(references.get(idx)?.external()?);
        }
        Ok(Self {
            primary: rows("PRIMARY_ALIGNMENT")?,
            secondary,
            references: count,
            external_references,
        })
    }
}

impl fmt::Display for AlignmentSummary {
    /// E.g. `1,000 primary, 20 secondary, against 25 references (3 external)`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} primary", format_count(self.primary))?;
        if let Some(secondary) = self.secondary {
            write!(f, ", {} secondary", format_count(secondary))?;
        }
        write!(
            f,
            ", against {} references ({} external)",
            self.references, self.external_references
        )
    }
}

/// Inspect `input` and render its report.
fn describe(input: &str, layout_spots: u64) -> Result<String> {
    let archive = Archive::inspect(input, None)?;
    let opened = archive.open_table()?;
    let sampled = layout_spots.min(archive.spot_count);
    let layout = if sampled == 0 {
        ReadLayout::default()
    } else {
        let last = archive.first_spot + sampled as i64 - 1;
        ReadLayout::sample(&opened.table, archive.first_spot, last)?
    };
    let alignments = match opened.database() {
        Some(database) if archive.aligned => Some(AlignmentSummary::read(database)?),
        _ => None,
    };
    let report = Report {
        accession: default_accession(input),
        archive: &archive,
        layout,
        spot_groups: spot_groups(&opened.table),
        alignments,
    };
    Ok(report.render())
}

/// Names of the spot groups in `STATS/SPOT_GROUP`, other than the unnamed `default`.
fn spot_groups(table: &VTable) -> Vec<String> {
    let Ok(metadata) = table.open_metadata_read() else { return Vec::new() };
    let Ok(node) = metadata.open_node_read("STATS/SPOT_GROUP") else { return Vec::new() };
    let children = node.list_children().unwrap_or_default();
    children
        .into_iter()
        .filter(|child| child != "default")
        .map(|child| {
            metadata
                .open_node_read(&format!("STATS/SPOT_GROUP/{child}"))
                .and_then(|group| group.read_attr("name"))
                .unwrap_or(child)
        })
        .collect()
}

/// The archive's kind: a flat table, or a database and how it holds its reads.
fn kind(archive: &Archive) -> String {
    match (&archive.table, archive.aligned) {
        (None, _) => "flat table".to_owned(),
        (Some(_), true) => "database, aligned (cSRA)".to_owned(),
        (Some(table), false) if archive.stores_bases => format!("database, reads in {table}"),
        (Some(table), false) => format!("database, reads in {table} (bases not stored directly)"),
    }
}

/// Where the archive's qualities come from, and their one value if every base has the same.
fn qualities(archive: &Archive) -> String {
    let source = match archive.quality_source {
        QualitySource::Stored => "stored",
        QualitySource::Lite => "SRA Lite: synthesised (Q30, or Q3 for reads failing the filter)",
        QualitySource::None if archive.has_qualities => "not stored: synthesised",
        QualitySource::None => "none",
    };
    match archive.constant_quality {
        Some(phred) => format!("{source}; every base Q{phred}"),
        None => source.to_owned(),
    }
}

/// The number of spot groups, and the first [`MAX_SPOT_GROUPS_LISTED`] names.
fn spot_group_summary(groups: &[String]) -> String {
    if groups.is_empty() {
        return "none".to_owned();
    }
    let listed = groups.iter().take(MAX_SPOT_GROUPS_LISTED).cloned().collect::<Vec<_>>().join(", ");
    match groups.len() {
        n if n > MAX_SPOT_GROUPS_LISTED => format!("{n}: {listed}, …"),
        n => format!("{n}: {listed}"),
    }
}

/// `key: spots` pairs of a distribution, e.g. `2 (pairs): 9,998 spots; 1: 2 spots`, in
/// descending order of spots, each key followed by `note(key)`.
fn distribution(
    spots_by_key: &BTreeMap<usize, u64>,
    note: impl Fn(usize) -> &'static str,
) -> String {
    let mut entries: Vec<(usize, u64)> = spots_by_key.iter().map(|(&k, &v)| (k, v)).collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    if entries.is_empty() {
        return "-".to_owned();
    }
    entries
        .iter()
        .map(|&(key, spots)| format!("{key}{}: {} spots", note(key), format_count(spots)))
        .collect::<Vec<_>>()
        .join("; ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archive::StoredTotals;

    const T: u8 = 0;
    const B: u8 = READ_TYPE_BIOLOGICAL;

    fn layout_of(spots: &[(&[u32], &[u8], &[u8])]) -> ReadLayout {
        let mut layout = ReadLayout::default();
        for (lens, types, filters) in spots {
            layout.add(lens, types, filters);
        }
        layout
    }

    fn archive() -> Archive {
        Archive {
            location: "/data/SRR1.sra".to_owned(),
            table: None,
            tables: Vec::new(),
            aligned: false,
            stores_bases: true,
            stores_cmp_read: false,
            platform: Some("ILLUMINA"),
            loader: Some("fastq-load.2.5.2 (2.5.2, Jun 19 2015)".to_owned()),
            first_spot: 1,
            spot_count: 2,
            quality_source: QualitySource::Stored,
            has_qualities: true,
            constant_quality: None,
            has_names: true,
            has_spot_groups: false,
            totals: StoredTotals { spots: Some(2), bases: Some(12), biological_bases: Some(12) },
        }
    }

    #[test]
    fn slot_lengths_are_summarised_over_non_empty_reads() {
        let layout = layout_of(&[(&[4, 0], &[B, B], &[0, 0]), (&[6, 3], &[B, B], &[0, 0])]);
        let first = &layout.slots[0];
        assert_eq!((first.min_len, first.max_len, first.total_len, first.empty), (4, 6, 10, 0));
        let second = &layout.slots[1];
        assert_eq!((second.min_len, second.max_len, second.empty), (3, 3, 1));
    }

    #[test]
    fn spots_are_counted_by_biological_and_technical_reads() {
        let layout = layout_of(&[
            (&[8, 50, 16, 48], &[T, B, T, B], &[0, 0, 0, 0]),
            (&[8, 50, 0, 0], &[T, B, T, B], &[0, 0, 0, 0]),
        ]);
        assert_eq!(layout.spots_by_bio_reads, BTreeMap::from([(1, 1), (2, 1)]));
        assert_eq!(layout.spots_by_technical_reads, BTreeMap::from([(1, 1), (2, 1)]));
        assert_eq!(layout.spots_by_read_count, BTreeMap::from([(4, 2)]));
    }

    #[test]
    fn slot_types_are_biological_technical_or_mixed() {
        let layout =
            layout_of(&[(&[4, 4, 4], &[B, T, B], &[0; 3]), (&[4, 4, 4], &[B, T, T], &[0; 3])]);
        let types: Vec<_> = layout.slots.iter().map(SlotStats::read_type).collect();
        assert_eq!(types, ["biological", "technical", "mixed"]);
    }

    #[test]
    fn spots_with_fewer_slots_count_only_the_slots_they_have() {
        let layout = layout_of(&[(&[4, 4], &[B, B], &[0, 0]), (&[4], &[B], &[0])]);
        assert_eq!(layout.slots[0].spots, 2);
        assert_eq!(layout.slots[1].spots, 1);
        assert_eq!(layout.spots_by_read_count, BTreeMap::from([(1, 1), (2, 1)]));
    }

    #[test]
    fn non_empty_reads_are_counted_by_filter() {
        let layout = layout_of(&[(&[4, 4, 0], &[B, B, B], &[0, 1, 3])]);
        assert_eq!(layout.reads_by_filter, [1, 1, 0, 0]);
    }

    #[test]
    fn report_shows_the_archive_and_its_layout() {
        let report = Report {
            accession: "SRR1".to_owned(),
            archive: &archive(),
            layout: layout_of(&[(&[3, 3], &[B, B], &[0, 0]), (&[3, 3], &[B, B], &[0, 0])]),
            spot_groups: Vec::new(),
            alignments: None,
        };
        let text = report.render();
        assert!(text.starts_with("SRR1\n"), "{text}");
        assert!(text.contains("  kind              flat table\n"), "{text}");
        assert!(text.contains("  platform          ILLUMINA\n"), "{text}");
        assert!(
            text.contains("  stored totals     2 spots, 12 bases, 12 biological bases\n"),
            "{text}"
        );
        assert!(text.contains("  biological reads  2 (pairs): 2 spots\n"), "{text}");
        assert!(text.contains("  technical reads   0: 2 spots\n"), "{text}");
        assert!(
            text.contains(
                "    1     biological             2         0       3      3.0       3\n"
            ),
            "{text}"
        );
    }

    #[test]
    fn report_shows_no_spot_range_for_an_empty_archive() {
        let archive = Archive { spot_count: 0, ..archive() };
        let report = Report {
            accession: "SRR1".to_owned(),
            archive: &archive,
            layout: ReadLayout::default(),
            spot_groups: Vec::new(),
            alignments: None,
        };
        let text = report.render();
        assert!(text.contains("  spots             0\n"), "{text}");
    }

    #[test]
    fn report_shows_alignments_of_an_aligned_archive() {
        let archive = Archive {
            table: Some("SEQUENCE".to_owned()),
            tables: vec![
                "PRIMARY_ALIGNMENT".to_owned(),
                "REFERENCE".to_owned(),
                "SEQUENCE".to_owned(),
            ],
            aligned: true,
            stores_bases: false,
            stores_cmp_read: true,
            ..archive()
        };
        let alignments = AlignmentSummary {
            primary: 1_000,
            secondary: Some(5),
            references: 25,
            external_references: 3,
        };
        let report = Report {
            accession: "SRR1".to_owned(),
            archive: &archive,
            layout: ReadLayout::default(),
            spot_groups: Vec::new(),
            alignments: Some(alignments),
        };
        let text = report.render();
        assert!(text.contains("  kind              database, aligned (cSRA)\n"), "{text}");
        assert!(text.contains("  alignments        1,000 primary, 5 secondary, against 25 references (3 external)\n"), "{text}");
    }

    #[test]
    fn biological_read_counts_are_annotated_with_where_fastq_sends_them() {
        let spots = BTreeMap::from([(0, 1), (1, 5), (2, 10), (3, 2)]);
        let text = distribution(&spots, |reads| match reads {
            1 => " (unpaired)",
            2 => " (pairs)",
            _ => " (dropped)",
        });
        assert_eq!(
            text,
            "2 (pairs): 10 spots; 1 (unpaired): 5 spots; 3 (dropped): 2 spots; 0 (dropped): 1 spots"
        );
    }

    #[test]
    fn many_spot_groups_are_truncated() {
        let groups: Vec<String> = (0..12).map(|i| format!("G{i}")).collect();
        let text = spot_group_summary(&groups);
        assert!(text.starts_with("12: G0, G1,"), "{text}");
        assert!(text.ends_with("G9, …"), "{text}");
    }

    #[test]
    fn qualities_mention_a_constant_value() {
        let archive = Archive { constant_quality: Some(93), ..archive() };
        assert_eq!(qualities(&archive), "stored; every base Q93");
    }
}
