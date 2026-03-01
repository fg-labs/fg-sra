//! fg-sratosam: High-performance SRA-to-SAM/BAM converter.

mod cli;
mod header;
mod aligned;
mod unaligned;
mod matecache;
mod record;
mod output;
mod cigar;
mod quality;
mod md_tag;

use anyhow::Result;
use clap::Parser;

use cli::Cli;

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.execute()
}
