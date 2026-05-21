mod benchmark;
mod cli;
mod frontend;
mod io;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_millis()
        .init();

    let cli = Cli::parse();
    match cli.command {
        Commands::Assemble(args) => frontend::assemble::run(args),
        Commands::Stats(args) => benchmark::run(args),
    }
}
