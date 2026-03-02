//! fg-sra: High-performance SRA toolkit.

mod aligned;
mod cigar;
mod cli;
mod header;
mod matecache;
mod md_tag;
mod output;
mod quality;
mod record;
mod unaligned;

use anyhow::Result;
use clap::Parser;

use cli::Cli;

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.execute()
}
