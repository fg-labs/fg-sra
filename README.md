[![Build](https://github.com/fg-labs/fg-sra/actions/workflows/ci.yml/badge.svg)](https://github.com/fg-labs/fg-sra/actions/workflows/ci.yml)
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
