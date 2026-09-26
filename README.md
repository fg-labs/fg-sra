[![Build](https://github.com/fg-labs/fg-sra/actions/workflows/ci.yml/badge.svg)](https://github.com/fg-labs/fg-sra/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/fg-sra.svg)](https://crates.io/crates/fg-sra)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/fg-sra.svg?label=bioconda)](https://bioconda.github.io/recipes/fg-sra/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/fg-sra/blob/main/LICENSE)
[![codecov](https://codecov.io/gh/fg-labs/fg-sra/branch/main/graph/badge.svg)](https://codecov.io/gh/fg-labs/fg-sra)

# fg-sra

High-performance SRA-to-SAM/BAM converter, replacing NCBI's `sam-dump` with
multi-threaded processing for significantly higher throughput.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fg-labs/fg-sra/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fg-labs/fg-sra/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fg-labs/fg-sra/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fg-sra and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Overview

fg-sra converts NCBI SRA archives to SAM or BAM format. It uses FFI
bindings to the NCBI VDB C library (`libncbi-vdb`) for reading SRA data and
processes references in parallel for high throughput.

Key features:
- **Multi-threaded** reference processing with ordered output
- **SAM and BAM** output (BAM via multi-threaded BGZF compression)
- **gzip/bzip2** compression for SAM output
- **FASTA/FASTQ** output modes
- **Region filtering** by genomic coordinates
- **Quality quantization**
- **Reference cache warming** via `cache-refs` to avoid resolver failures under load
- **Mate cache** for proper SAM flag and mate-pair information
- **Archive descriptions** via `info`: kind, totals, qualities, names and read layout
- **Paired FASTQ** via `fastq`: spot order, mates paired, multi-threaded BGZF, byte-identical at any thread count, checked against the archive's stored totals, for unaligned and aligned (cSRA) archives

For aligned runs, `fg-sra` reconstructs each read from the stored alignment
deltas rather than reading the virtual `READ` column: it preloads each aligned
reference sequence into memory (single-threaded, once per preload batch), then
rebuilds `READ` on worker threads with a port of ncbi-vdb's own reconstruction
routine. This keeps
worker threads off libncbi-vdb's `REFERENCE` sub-select, whose blob cache is not
thread-safe. The preload costs roughly one byte per reference base; references
are processed in batches to cap this (default 1 GiB, overridable with the
`FG_SRA_REF_PRELOAD_BUDGET_MB` environment variable). The budget targets the
retained reference-store bytes, not strict peak memory: a single reference
larger than the budget still forms its own batch, and loading a reference
transiently holds both its raw 4na and mapped bases, so peak memory can exceed
the configured value by roughly the largest single reference.

The following `sam-dump` options are accepted but **not yet supported**:
- `--hide-identical` — output `=` for bases matching reference
- `--with-md-flag` — compute and output the MD tag
- `--rna-splicing` / `--rna-splice-level` / `--rna-splice-log` — RNA splice detection

These require reference sequence access via VDB FFI that has not yet been implemented.

## Installation

### Building from source (vendored)

```bash
git clone --recurse-submodules https://github.com/fg-labs/fg-sra
cd fg-sra
cargo build --release
```

#### Prerequisites

- Rust (stable toolchain)
- CMake (for building the vendored ncbi-vdb C library)

The `vendored` feature is enabled by default, building ncbi-vdb from the
git submodule automatically during `cargo build`.

So is `zlib-ng`, which links zlib-ng (built from source through `libz-sys`) in place of the zlib bundled with ncbi-vdb, for faster decompression of archives; `--no-default-features --features vendored` keeps the bundled zlib.

### Building with a pre-built VDB

To link against a system-installed ncbi-vdb instead of building from source:

```bash
export VDB_INCDIR=/path/to/ncbi-vdb/interfaces
export VDB_LIBDIR=/path/to/lib/containing/libncbi-vdb.a
cargo build --release --no-default-features
```

## Usage

```bash
# Convert an SRA accession to SAM
fg-sra tosam SRR390728

# Convert to BAM
fg-sra tosam --output-format bam --output-file output.bam SRR390728

# Primary alignments only, with unaligned reads
fg-sra tosam -1 -u SRR390728

# Filter by region
fg-sra tosam --aligned-region chr1:1000000-2000000 SRR390728

# Multi-threaded with 8 threads
fg-sra tosam -t 8 SRR390728
```

For full usage, run:

```bash
fg-sra tosam --help
```

### Converting to FASTQ

`fg-sra fastq` writes an archive's reads in spot order, with mates paired. Each spot's non-empty biological reads decide where it goes: two go to `--r1`/`--r2` (or `--interleaved`), and one goes to `--unpaired`. Spots with no biological reads or more than two are dropped and counted, as are spots whose output wasn't given; it is an error if nothing at all is written. Outputs ending `.gz` or `.bgz` are BGZF-compressed, and any output may be `-` for stdout.

```bash
# Pairs only, BGZF-compressed, with a metrics TSV
fg-sra fastq SRR13232999.sra -1 r1.fq.gz -2 r2.fq.gz -m metrics.tsv

# fasterq-dump/ENA-shaped output: pairs plus orphans and single-end reads
fg-sra fastq SRR000001.sra -1 SRR000001_1.fastq.gz -2 SRR000001_2.fastq.gz -u SRR000001.fastq.gz

# Stream interleaved pairs to an aligner
fg-sra fastq SRR2584863.sra -p - | bwa mem -p ref.fa -

# Technical reads (e.g. barcodes and indexes), one file per technical read, in step with the biological reads
fg-sra fastq SRR13450125.sra -1 r1.fq.gz -2 r2.fq.gz --technical 'tech.{i}.fq.gz'

# Original read names, fasterq-dump style, for the first million spots
fg-sra fastq SRR2584863.sra --defline '$ac.$si $sn length=$rl' --max-spot-id 1000000 -1 a.fq -2 b.fq

# A random 10% of spots, the same on every run and at any thread count
fg-sra fastq SRR2584863.sra --subsample-fraction 0.1 -1 r1.fq.gz -2 r2.fq.gz
```

Reads are named `<accession>.<spot>` by default, identically in every file. `--defline` takes a template with `$ac` (accession), `$si` (spot id), `$sn` (original name), `$sg` (spot group), `$ri` (read number within its type) and `$rl` (read length); the `+` line is always bare. `--min-read-len` and `--read-filter` test biological reads and drop the whole spot, so the outputs never go out of step. Runs loaded from BAM mark duplicates `criteria` and QC failures `reject`, so `--read-filter pass` drops duplicates too. `--subsample-fraction` keeps a random fraction of spots, again whole spots with all their reads. Which spots are kept depends only on `--subsample-seed` (default 42) and each spot's id, so the output is the same at any thread count, any spot range keeps about the fraction, and with one seed a smaller fraction keeps a subset of a larger one's spots. Skipped spots aren't read, so the totals aren't checked against the archive's.

Porting from sra-tools:

| sra-tools | fg-sra fastq |
|---|---|
| `fasterq-dump` (split-3) | `-1 X_1.fq -2 X_2.fq -u X.fq` |
| `--split-spot -Z` | `-p -` for pairs only: single-read spots need their own `-u X.fq`, since only one output can be stdout |
| `--split-files --include-technical` | `-1`/`-2` or `-u`, plus `--technical 'X_tech{i}.fq'` |
| `-M N` | `--min-read-len N` |
| `-R pass` | `--read-filter pass` |
| `-N A -X B` | `--min-spot-id A --max-spot-id B` |
| `--seq-defline T --qual-defline '+'` | `--defline T` |
| `--gzip`, pigz | a `.gz` output path |

Aligned (cSRA) archives are converted too. Their references are loaded into memory first (about a byte per reference base, e.g. ~3 GB for a human genome), and each aligned read is rebuilt from its stored alignment, thread-safely, rather than through libncbi-vdb's virtual `READ` column. External references must be available locally (e.g. beside the archive, as `prefetch` puts them) or over the network; `--offline` turns network access off entirely. SRA Lite archives are converted with a warning, since their qualities are synthesised. Colour-space (SOLiD) runs are written in base space, as `fasterq-dump` writes them; `fastq-dump` needs `-B` for the same. PacBio and Oxford Nanopore native databases that have a `CONSENSUS` table are read from it by default (`--table auto`), as `fasterq-dump` reads them, rather than from `SEQUENCE`, which holds each molecule's subreads or strands as separate reads; `--table SEQUENCE` reads those instead, and `fastq-dump` needs `--table CONSENSUS` for the same.

### Describing an Archive

`fg-sra info` describes one or more archives: kind (flat table or database, aligned or not), platform, loader, stored totals, whether qualities and original read names were kept, spot groups, alignments and references, and the read layout of the first spots (10,000 by default; `--layout-spots` sets how many), so a layout that changes later in a run isn't seen. The layout lists each read slot's type and lengths, and how many biological and technical reads spots have, which shows the `fg-sra fastq` outputs an archive needs and which technical read is which.

```bash
fg-sra info SRR000001.sra
```

```
SRR000001
  kind              flat table
  platform          454
  ...
  read layout       first 10,000 of 470,985 spots only (--layout-spots samples more)
    slot  type           non-empty     empty     min     mean     max
    1     technical         10,000         0       4      4.0       4
    2     biological        10,000         0       1    168.8     685
    3     technical          5,648     4,352      44     44.0      44
    4     biological         5,634     4,366       1    110.0     392
  reads per spot    4: 10,000 spots
  biological reads  2 (pairs): 5,634 spots; 1 (unpaired): 4,366 spots
  technical reads   2: 5,648 spots; 1: 4,352 spots
```

### Pre-caching References

When running many `fg-sra tosam` conversions concurrently, the VDB reference
resolver can fail under heavy load. Use `cache-refs` to serially pre-populate
the local reference cache before launching concurrent conversions:

```bash
# Cache references for a list of accessions
fg-sra cache-refs SRR622461 SRR765989 SRR341578

# Then run conversions concurrently
parallel fg-sra tosam {} ::: SRR622461 SRR765989 SRR341578
```

The `cache-refs` command processes accessions sequentially, resolving all
reference sequence dependencies via VDB and caching them locally. Subsequent
`tosam` runs will find these references in the local cache without network access.

## Performance

SRR20022182 converted to coordinate-sorted BAM (piped through `samtools sort`)
completes in ~5s wall-clock time with ~400 MB peak memory. Use `--threads` to
enable multi-threaded aligned read processing.

`fg-sra fastq` converts SRR13232999 (9.5M paired spots, 2.9 Gbases, unaligned) to BGZF R1 and R2 at level 1 in ~2.8 s wall-clock with 8 threads on an Apple-silicon laptop; `fasterq-dump` takes ~7 s to write the same reads uncompressed. On the aligned SRR1574798 (26M spots, 5.3 Gbases, 39 embedded references) it takes ~11.4 s and 3.8 GB, against ~25 s and 4.3 GB for `fasterq-dump`, with identical output. At 16 threads on an M4 Max, ERR17774045 (69.8M paired spots, 14.1 Gbases, unaligned) takes ~6.8 s, against ~99 s for `fasterq-dump` writing uncompressed, and SRR7107873 (70.3M spots, 14.1 Gbases, aligned to 39 external RefSeq references) takes ~20 s and 5.2 GB, against ~53 s and 6.3 GB.

## Workspace Structure

```
fg-sra/
├── crates/
│   ├── fg-sra-vdb-sys/    # Raw FFI bindings to libncbi-vdb
│   ├── fg-sra-vdb/        # Safe Rust wrappers over VDB
│   └── fg-sra/            # Binary crate (the converter)
└── vendor/
    └── ncbi-vdb/           # Vendored VDB library (git submodule)
```

## Resources

- [Issues](https://github.com/fg-labs/fg-sra/issues): Report a bug or request a feature
- [Pull requests](https://github.com/fg-labs/fg-sra/pulls): Submit a patch or new feature
- [Contributors guide](https://github.com/fg-labs/fg-sra/blob/main/CONTRIBUTING.md)
- [License](https://github.com/fg-labs/fg-sra/blob/main/LICENSE): Released under the MIT license

## Authors

- [Nils Homer](https://github.com/nh13)

## Sponsors

Development of fg-sra is supported by [Fulcrum Genomics](https://www.fulcrumgenomics.com).

[Become a sponsor](https://github.com/sponsors/fulcrumgenomics)

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fg-labs/fg-sra/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fg-labs/fg-sra/issues), and better yet a [pull request](https://github.com/fg-labs/fg-sra/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.
