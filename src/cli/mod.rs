use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "lineabact",
    version,
    about = "Short-read-first assembly front-end for bacterial isolates"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Parser, Debug)]
pub enum Commands {
    /// Assemble paired-end Illumina reads with a SPAdes-backed front-end
    Assemble(AssembleArgs),

    /// Generate or materialize benchmark plans and reports
    Stats(StatsArgs),
}

#[derive(Parser, Debug)]
pub struct AssembleArgs {
    /// Forward reads (R1) FASTQ.gz
    #[arg(long)]
    pub r1: PathBuf,

    /// Reverse reads (R2) FASTQ.gz
    #[arg(long)]
    pub r2: PathBuf,

    /// Output directory
    #[arg(long, default_value = "lineabact_out")]
    pub outdir: PathBuf,

    /// Backend used to materialize SPAdes-style outputs
    #[arg(long, value_enum, default_value_t = AssembleBackend::Spades)]
    pub backend: AssembleBackend,

    /// Emit read/backend plans without running the backend
    #[arg(long, default_value_t = false)]
    pub dry_run: bool,

    /// k-mer size (must be odd, ≤ 63)
    #[arg(long, default_value_t = 55)]
    pub k: u8,

    /// Number of threads
    #[arg(long, default_value_t = 4)]
    pub threads: usize,

    /// SPAdes memory limit in GB
    #[arg(long)]
    pub spades_memory_gb: Option<usize>,

    /// SPAdes temporary directory
    #[arg(long)]
    pub spades_tmp_dir: Option<PathBuf>,

    /// Request SPAdes GFA v1.1 output
    #[arg(long, default_value_t = false)]
    pub spades_gfa11: bool,

    /// Deterministic downsampling fraction applied to read pairs before k-mer counting
    #[arg(long)]
    pub downsample_fraction: Option<f64>,

    /// Target short-read depth used to derive downsampling fraction
    #[arg(long)]
    pub target_coverage: Option<f64>,

    /// Estimated genome size in base pairs, used for depth estimation/downsampling
    #[arg(long)]
    pub genome_size_bp: Option<u64>,

    /// Seed used for deterministic downsampling
    #[arg(long, default_value_t = 0)]
    pub downsample_seed: u64,

    /// Trim adapters before backend assembly
    #[arg(long, default_value_t = false)]
    pub trim_adapters: bool,

    /// Preferred trimming tool/executable
    #[arg(long, default_value = "fastp")]
    pub trim_tool: String,

    /// Compute SHA256 digests for input/output read files
    #[arg(long, default_value_t = false)]
    pub read_sha256: bool,

    /// SPAdes executable name recorded in the backend plan
    #[arg(long, default_value = "spades.py")]
    pub spades_executable: String,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum AssembleBackend {
    Spades,
    Mock,
}

#[derive(Parser, Debug)]
pub struct StatsArgs {
    #[command(subcommand)]
    pub command: StatsCommands,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum BenchmarkScheduler {
    Local,
    Slurm,
}

#[derive(Subcommand, Debug)]
pub enum StatsCommands {
    #[command(name = "benchmark-plan")]
    /// Build a benchmark plan spanning LineaBact, SPAdes, Shovill and Unicycler
    Plan(BenchmarkPlanArgs),

    #[command(name = "benchmark-run")]
    /// Convert a benchmark plan into a report scaffold
    Run(BenchmarkRunArgs),

    #[command(name = "benchmark-merge")]
    /// Merge runtime comparisons from a secondary report into a formal report
    Merge(BenchmarkMergeArgs),
}

#[derive(Parser, Debug)]
pub struct BenchmarkPlanArgs {
    /// Tab-delimited manifest for multi-sample benchmark planning
    #[arg(long)]
    pub manifest: Option<PathBuf>,

    /// Maximum number of samples loaded from `--manifest`
    #[arg(long, default_value_t = 3)]
    pub sample_limit: usize,

    /// Unicycler fixture root used as the current-stage benchmark source
    #[arg(long, default_value = "reference_tools/Unicycler-main")]
    pub fixture_root: PathBuf,

    /// Output directory for the generated plan and report scaffold
    #[arg(long, default_value = "benchmarks/generated/reference_validation")]
    pub outdir: PathBuf,

    /// Sample name recorded in the benchmark plan
    #[arg(long, default_value = "unicycler_short_read_fixture")]
    pub sample_id: String,

    /// Number of threads recorded in planned commands
    #[arg(long, default_value_t = 4)]
    pub threads: usize,

    /// k-mer size recorded in SPAdes-backed commands
    #[arg(long, default_value_t = 55)]
    pub k: u8,

    /// LineaBact executable name recorded in the plan
    #[arg(long, default_value = "lineabact")]
    pub lineabact_executable: String,

    /// SPAdes executable name recorded in the plan
    #[arg(long, default_value = "spades.py")]
    pub spades_executable: String,

    /// Shovill executable name recorded in the plan
    #[arg(long, default_value = "shovill")]
    pub shovill_executable: String,

    /// Unicycler executable name recorded in the plan
    #[arg(long, default_value = "unicycler")]
    pub unicycler_executable: String,
}

#[derive(Parser, Debug, Clone)]
pub struct BenchmarkRunArgs {
    /// Benchmark plan JSON produced by `benchmark-plan`
    #[arg(long)]
    pub plan: PathBuf,

    /// Output directory for the benchmark report
    #[arg(long, default_value = "benchmarks/generated/reference_validation")]
    pub outdir: PathBuf,

    /// Stop execution on first failed benchmark case
    #[arg(long, default_value_t = false)]
    pub stop_on_error: bool,

    /// Number of repeated executions per benchmark case
    #[arg(long, default_value_t = 3)]
    pub repeat_count: usize,

    /// Execution backend used to materialize benchmark cases
    #[arg(long, value_enum, default_value_t = BenchmarkScheduler::Local)]
    pub scheduler: BenchmarkScheduler,

    /// Slurm partition used when `--scheduler slurm`; use `auto` or a comma-separated list to select from multiple partitions
    #[arg(long, default_value = "auto")]
    pub slurm_partition: String,

    /// Conda base directory sourced by the generated Slurm script
    #[arg(long, default_value = "/hpcfs/fpublic/app/miniforge3/conda")]
    pub slurm_conda_base: String,

    /// Conda environment activated by the generated Slurm script
    #[arg(long, default_value = "LineaBact")]
    pub slurm_conda_env: String,

    /// CPU count requested per Slurm benchmark case
    #[arg(long, default_value_t = 4)]
    pub slurm_cpus_per_task: usize,

    /// Memory in GB requested per Slurm benchmark case
    #[arg(long)]
    pub slurm_mem_gb: Option<usize>,

    /// Wall-clock time requested per Slurm benchmark case
    #[arg(long, default_value = "12:00:00")]
    pub slurm_time: String,

    /// Write Slurm scripts without submitting them
    #[arg(long, default_value_t = false)]
    pub slurm_dry_run: bool,
}

#[derive(Parser, Debug)]
pub struct BenchmarkMergeArgs {
    /// Primary benchmark report JSON to update in place
    #[arg(long)]
    pub report: PathBuf,

    /// Secondary benchmark report JSON used to compute runtime comparisons
    #[arg(long)]
    pub compare_report: PathBuf,
}
