use anyhow::Result;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::io;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadReport {
    pub input_pairs: u64,
    pub retained_pairs: u64,
    pub dropped_pairs: u64,
    pub input_reads_r1: u64,
    pub input_reads_r2: u64,
    pub retained_reads_r1: u64,
    pub retained_reads_r2: u64,
    pub input_bases_r1: u64,
    pub input_bases_r2: u64,
    pub retained_bases_r1: u64,
    pub retained_bases_r2: u64,
    pub mean_read_length_r1: f64,
    pub mean_read_length_r2: f64,
    pub input_total_bases: u64,
    pub retained_total_bases: u64,
    pub genome_size_bp: Option<u64>,
    pub estimated_depth: Option<f64>,
    pub target_coverage: Option<f64>,
    pub derived_downsample_fraction: Option<f64>,
    pub applied_downsample_fraction: f64,
    pub downsample_fraction: Option<f64>,
    pub downsample_seed: u64,
    pub trim_adapters: bool,
    pub trim_tool: String,
    pub read_sha256: bool,
    pub input_sha256_r1: Option<String>,
    pub input_sha256_r2: Option<String>,
    pub output_sha256_r1: Option<String>,
    pub output_sha256_r2: Option<String>,
    pub optimized_r1: String,
    pub optimized_r2: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendPlan {
    pub backend: String,
    pub executable: String,
    pub command: String,
    pub version: Option<String>,
    pub dry_run: bool,
    pub k: u8,
    pub threads: usize,
    pub memory_gb: Option<usize>,
    pub tmp_dir: Option<String>,
    pub gfa11: bool,
    pub input_r1: String,
    pub input_r2: String,
    pub output_dir: String,
    pub stdout_log: String,
    pub stderr_log: String,
    pub diagnostics_json: String,
    pub generated_files: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphSummary {
    pub segment_count: usize,
    pub link_count: usize,
    pub path_record_count: usize,
    pub total_segment_bases: usize,
    pub mean_segment_bases: f64,
    pub max_segment_bases: usize,
    pub min_segment_bases: usize,
    pub mean_depth: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BridgeSummary {
    pub path_record_count: usize,
    pub forward_path_records: usize,
    pub reverse_path_records: usize,
    pub candidate_count: usize,
    pub conflict_count: usize,
    pub applied_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyStats {
    pub num_contigs: usize,
    pub total_length: usize,
    pub n50: usize,
    pub longest_contig: usize,
    pub shortest_contig: usize,
    pub mean_length: f64,
    pub gc_content: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamsSnapshot {
    pub version: String,
    pub command: String,
    pub backend: String,
    pub dry_run: bool,
    pub input_r1: String,
    pub input_r2: String,
    pub outdir: String,
    pub spades_executable: String,
    pub k: u8,
    pub threads: usize,
    pub spades_memory_gb: Option<usize>,
    pub spades_tmp_dir: Option<String>,
    pub spades_gfa11: bool,
    pub downsample_fraction: Option<f64>,
    pub target_coverage: Option<f64>,
    pub genome_size_bp: Option<u64>,
    pub downsample_seed: u64,
    pub trim_adapters: bool,
    pub trim_tool: String,
    pub read_sha256: bool,
    pub read_summary: String,
    pub backend_plan: String,
    pub graph_summary: Option<String>,
    pub bridge_summary: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunSummary {
    pub version: String,
    pub command: String,
    pub backend: String,
    pub dry_run: bool,
    pub input_r1: String,
    pub input_r2: String,
    pub outdir: String,
    pub spades_executable: String,
    pub k: u8,
    pub threads: usize,
    pub spades_memory_gb: Option<usize>,
    pub spades_tmp_dir: Option<String>,
    pub spades_gfa11: bool,
    pub downsample_fraction: Option<f64>,
    pub target_coverage: Option<f64>,
    pub genome_size_bp: Option<u64>,
    pub downsample_seed: u64,
    pub trim_adapters: bool,
    pub trim_tool: String,
    pub read_sha256: bool,
    pub elapsed_seconds: f64,
    pub read_report: ReadReport,
    pub backend_plan: BackendPlan,
    pub graph_summary: Option<GraphSummary>,
    pub bridge_summary: Option<BridgeSummary>,
    pub assembly_stats: Option<AssemblyStats>,
    pub output_files: Vec<String>,
}

impl RunSummary {
    pub fn write_json(&self, path: &Path) -> Result<()> {
        io::write_json(path, self)
    }
}

impl ParamsSnapshot {
    pub fn write_toml(&self, path: &Path) -> Result<()> {
        io::write_toml(path, self)
    }
}

impl AssemblyStats {
    pub fn from_fasta(path: &Path) -> Result<Self> {
        let mut lengths = Vec::new();
        let mut gc_count = 0usize;

        let mut reader = parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            let record = record?;
            let seq = record.seq();
            let seq = seq.as_ref();
            let length = seq.len();
            if length == 0 {
                continue;
            }
            lengths.push(length);
            gc_count += seq
                .iter()
                .filter(|&&base| matches!(base, b'G' | b'g' | b'C' | b'c'))
                .count();
        }

        let num_contigs = lengths.len();
        let total_length: usize = lengths.iter().sum();
        let n50 = n50(&lengths);
        let longest_contig = lengths.iter().copied().max().unwrap_or(0);
        let shortest_contig = lengths.iter().copied().min().unwrap_or(0);
        let mean_length = if num_contigs == 0 {
            0.0
        } else {
            total_length as f64 / num_contigs as f64
        };
        let gc_content = if total_length == 0 {
            0.0
        } else {
            gc_count as f64 / total_length as f64
        };

        Ok(Self {
            num_contigs,
            total_length,
            n50,
            longest_contig,
            shortest_contig,
            mean_length,
            gc_content,
        })
    }

    pub fn write_tsv(&self, path: &Path) -> Result<()> {
        io::ensure_parent_dir(path)?;
        let file = fs::File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "metric\tvalue")?;
        writeln!(writer, "num_contigs\t{}", self.num_contigs)?;
        writeln!(writer, "total_length\t{}", self.total_length)?;
        writeln!(writer, "n50\t{}", self.n50)?;
        writeln!(writer, "longest_contig\t{}", self.longest_contig)?;
        writeln!(writer, "shortest_contig\t{}", self.shortest_contig)?;
        writeln!(writer, "mean_length\t{:.1}", self.mean_length)?;
        writeln!(writer, "gc_content\t{:.4}", self.gc_content)?;
        writer.flush()?;
        Ok(())
    }
}

fn n50(lengths: &[usize]) -> usize {
    if lengths.is_empty() {
        return 0;
    }

    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let target = sorted.iter().sum::<usize>() / 2;
    let mut cumulative = 0usize;
    for length in sorted {
        cumulative += length;
        if cumulative >= target {
            return length;
        }
    }
    0
}
