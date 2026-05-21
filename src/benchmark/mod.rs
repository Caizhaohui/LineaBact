use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::thread;
use std::time::{Duration, Instant};

use crate::cli::{
    BenchmarkMergeArgs, BenchmarkPlanArgs, BenchmarkRunArgs, BenchmarkScheduler, StatsArgs,
    StatsCommands,
};
use crate::frontend::contract;
use crate::io;

const SCHEMA_VERSION: u32 = 1;
const SHOVILL_RUNTIME_TARGET_RATIO: f64 = 1.0;
const SLURM_SUBMISSION_RETRIES: usize = 3;
const SLURM_QUERY_RETRIES: usize = 3;
const SLURM_COMMAND_TIMEOUT: Duration = Duration::from_secs(20);
const DEFAULT_SLURM_PARTITION_CANDIDATES: &[&str] =
    &["qcpu_23if", "qcpu_23i", "qcpu_23a", "qcpu_18i"];
const QUALITY_N50_MIN_RATIO: f64 = 0.90;
const QUALITY_N50_PREFERRED_RATIO: f64 = 0.95;
const QUALITY_SEED_K: usize = 19;
const QUALITY_REFERENCE_SEED_STEP: usize = 17;
const QUALITY_SEED_STEP: usize = 97;
const QUALITY_MAX_KMER_HITS: usize = 8;
const QUALITY_CLUSTER_BUCKET_BP: isize = 128;
const QUALITY_SECONDARY_CLUSTER_RATIO: f64 = 0.5;
const QUALITY_MISASSEMBLY_GAP_BP: usize = 1000;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkSample {
    pub sample_id: String,
    pub r1: String,
    pub r2: String,
    pub reference: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkCase {
    pub sample_id: String,
    pub tool: String,
    pub executable: String,
    pub command: String,
    pub outdir: String,
    pub expected_dirs: Vec<String>,
    pub expected_files: Vec<String>,
    pub notes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkPlan {
    pub schema_version: u32,
    pub fixture_root: String,
    pub samples: Vec<BenchmarkSample>,
    pub cases: Vec<BenchmarkCase>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkCaseReport {
    pub sample_id: String,
    pub repeat_index: usize,
    pub repeat_count: usize,
    pub scheduler: String,
    pub tool: String,
    pub executable: String,
    pub command: String,
    pub outdir: String,
    pub expected_dirs: Vec<String>,
    pub expected_files: Vec<String>,
    pub status: String,
    pub exit_code: Option<i32>,
    pub elapsed_seconds: f64,
    #[serde(default = "default_elapsed_seconds_source")]
    pub elapsed_seconds_source: String,
    pub missing_dirs: Vec<String>,
    pub missing_files: Vec<String>,
    pub stdout_log: String,
    pub stderr_log: String,
    pub submission_script: Option<String>,
    pub submission_command: Option<String>,
    pub job_id: Option<String>,
    #[serde(default)]
    pub quality_metrics: Option<CaseQualityMetrics>,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkReport {
    pub schema_version: u32,
    pub plan: String,
    pub case_count: usize,
    #[serde(default)]
    pub samples: Vec<BenchmarkSample>,
    pub cases: Vec<BenchmarkCaseReport>,
    #[serde(default)]
    pub runtime_comparisons: Vec<RuntimeComparison>,
    #[serde(default)]
    pub quality_comparisons: Vec<QualityComparison>,
    #[serde(default)]
    pub sample_acceptance: Vec<SampleAcceptance>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RuntimeComparison {
    pub sample_id: String,
    pub baseline_tool: String,
    pub candidate_tool: String,
    pub baseline_repeat_count: usize,
    pub candidate_repeat_count: usize,
    pub baseline_median_elapsed_seconds: f64,
    pub candidate_median_elapsed_seconds: f64,
    pub median_runtime_ratio_vs_baseline: f64,
    pub median_speedup_vs_baseline: f64,
    pub target_runtime_ratio: f64,
    pub meets_target: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaseQualityMetrics {
    pub assembly_path: String,
    pub reference_path: String,
    pub assembly_source: String,
    pub reference_bases: usize,
    pub covered_reference_bases: usize,
    pub aligned_assembly_bases: usize,
    pub genome_fraction: f64,
    pub contig_count: usize,
    pub total_length: usize,
    pub longest_contig: usize,
    pub n50: usize,
    pub mismatch_bases: usize,
    pub indel_bases: usize,
    pub mismatch_rate: f64,
    pub indel_rate: f64,
    pub misassembly_count: usize,
    pub evaluation_method: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityComparison {
    pub sample_id: String,
    pub baseline_tool: String,
    pub candidate_tool: String,
    pub baseline_repeat_count: usize,
    pub candidate_repeat_count: usize,
    pub baseline_median_n50: f64,
    pub candidate_median_n50: f64,
    pub n50_ratio_vs_baseline: f64,
    pub target_n50_ratio: f64,
    pub preferred_n50_ratio: f64,
    pub baseline_median_genome_fraction: f64,
    pub candidate_median_genome_fraction: f64,
    pub baseline_median_misassemblies: f64,
    pub candidate_median_misassemblies: f64,
    pub baseline_median_mismatch_rate: f64,
    pub candidate_median_mismatch_rate: f64,
    pub baseline_median_indel_rate: f64,
    pub candidate_median_indel_rate: f64,
    pub passes_n50: bool,
    pub passes_genome_fraction: bool,
    pub passes_misassemblies: bool,
    pub passes_mismatch_rate: bool,
    pub passes_indel_rate: bool,
    pub meets_target: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleAcceptance {
    pub sample_id: String,
    pub runtime_meets_target: bool,
    pub quality_meets_target: bool,
    pub accepted: bool,
    pub failure_reasons: Vec<String>,
}

fn default_elapsed_seconds_source() -> String {
    "unknown".to_string()
}

pub fn run(args: StatsArgs) -> Result<()> {
    match args.command {
        StatsCommands::Plan(args) => benchmark_plan(args),
        StatsCommands::Run(args) => benchmark_run(args),
        StatsCommands::Merge(args) => benchmark_merge(args),
    }
}

pub fn benchmark_plan(args: BenchmarkPlanArgs) -> Result<()> {
    fs::create_dir_all(&args.outdir)?;
    let plan = build_plan(&args)?;

    let plan_json = contract::root_path(&args.outdir, "benchmark_plan.json");
    let plan_tsv = contract::root_path(&args.outdir, "benchmark_plan.tsv");

    io::write_json(&plan_json, &plan)?;
    write_plan_tsv(&plan, &plan_tsv)?;

    Ok(())
}

pub fn benchmark_run(args: BenchmarkRunArgs) -> Result<()> {
    let args = normalize_benchmark_run_args(args)?;
    fs::create_dir_all(&args.outdir)?;
    let plan_text = fs::read_to_string(&args.plan)
        .with_context(|| format!("failed to read benchmark plan {}", args.plan.display()))?;
    let plan: BenchmarkPlan = serde_json::from_str(&plan_text)
        .with_context(|| format!("failed to parse benchmark plan {}", args.plan.display()))?;

    let mut report = BenchmarkReport {
        schema_version: SCHEMA_VERSION,
        plan: args.plan.display().to_string(),
        case_count: plan.cases.len(),
        samples: plan.samples.clone(),
        cases: run_cases(&plan.cases, &args)?,
        runtime_comparisons: Vec::new(),
        quality_comparisons: Vec::new(),
        sample_acceptance: Vec::new(),
    };
    hydrate_case_reports(&report.samples, &mut report.cases);
    report.runtime_comparisons = build_runtime_comparisons(&report.cases);
    report.quality_comparisons = build_quality_comparisons(&report.samples, &report.cases);
    report.sample_acceptance =
        build_sample_acceptance(&report.runtime_comparisons, &report.quality_comparisons);

    let report_json = contract::root_path(&args.outdir, "benchmark_report.json");
    let report_tsv = contract::root_path(&args.outdir, "benchmark_report.tsv");

    io::write_json(&report_json, &report)?;
    write_report_tsv(&report, &report_tsv)?;

    Ok(())
}

fn normalize_benchmark_run_args(mut args: BenchmarkRunArgs) -> Result<BenchmarkRunArgs> {
    if matches!(args.scheduler, BenchmarkScheduler::Slurm) {
        args.slurm_partition = resolve_slurm_partition(&args.slurm_partition)?;
    }
    Ok(args)
}

fn ensure_report_samples(report: &mut BenchmarkReport) -> Result<()> {
    if !report.samples.is_empty() {
        return Ok(());
    }

    let plan_path = Path::new(&report.plan);
    let plan_text = fs::read_to_string(plan_path)
        .with_context(|| format!("failed to read benchmark plan {}", plan_path.display()))?;
    let plan: BenchmarkPlan = serde_json::from_str(&plan_text)
        .with_context(|| format!("failed to parse benchmark plan {}", plan_path.display()))?;
    report.samples = plan.samples;
    Ok(())
}

pub fn benchmark_merge(args: BenchmarkMergeArgs) -> Result<()> {
    let report_text = fs::read_to_string(&args.report)
        .with_context(|| format!("failed to read benchmark report {}", args.report.display()))?;
    let compare_text = fs::read_to_string(&args.compare_report).with_context(|| {
        format!(
            "failed to read comparison benchmark report {}",
            args.compare_report.display()
        )
    })?;

    let mut report: BenchmarkReport = serde_json::from_str(&report_text)
        .with_context(|| format!("failed to parse benchmark report {}", args.report.display()))?;
    let mut compare_report: BenchmarkReport =
        serde_json::from_str(&compare_text).with_context(|| {
            format!(
                "failed to parse comparison benchmark report {}",
                args.compare_report.display()
            )
        })?;

    ensure_report_samples(&mut report)?;
    ensure_report_samples(&mut compare_report)?;
    hydrate_case_reports(&report.samples, &mut report.cases);
    hydrate_case_reports(&compare_report.samples, &mut compare_report.cases);

    let mut cases = report.cases.clone();
    cases.extend(compare_report.cases.iter().cloned());
    let mut samples = report.samples.clone();
    for sample in &compare_report.samples {
        if !samples
            .iter()
            .any(|existing| existing.sample_id == sample.sample_id)
        {
            samples.push(sample.clone());
        }
    }
    report.runtime_comparisons = build_runtime_comparisons(&cases);
    report.quality_comparisons = build_quality_comparisons(&samples, &cases);
    report.sample_acceptance =
        build_sample_acceptance(&report.runtime_comparisons, &report.quality_comparisons);
    report.samples = samples;

    let report_dir = args.report.parent().with_context(|| {
        format!(
            "benchmark report {} has no parent directory",
            args.report.display()
        )
    })?;
    let report_json = report_dir.join("benchmark_report.json");
    let report_tsv = report_dir.join("benchmark_report.tsv");
    io::write_json(&report_json, &report)?;
    write_report_tsv(&report, &report_tsv)?;

    Ok(())
}

fn run_cases(cases: &[BenchmarkCase], args: &BenchmarkRunArgs) -> Result<Vec<BenchmarkCaseReport>> {
    let mut reports = Vec::with_capacity(cases.len());
    let mut halted = false;

    for case in cases {
        if halted {
            reports.push(skipped_case_report(
                case,
                &args.scheduler,
                1,
                args.repeat_count,
            ));
            continue;
        }

        let repeat_count = args.repeat_count.max(1);
        for repeat_index in 1..=repeat_count {
            let repeated_case = if repeat_count == 1 {
                case.clone()
            } else {
                repeated_case(case, repeat_index)
            };
            let report = match args.scheduler {
                BenchmarkScheduler::Local => {
                    execute_case_locally(&repeated_case, repeat_index, repeat_count)?
                }
                BenchmarkScheduler::Slurm => {
                    submit_case_to_slurm(&repeated_case, args, repeat_index, repeat_count)?
                }
            };
            let failed = report.status != "ok" && report.status != "submitted";
            reports.push(report);
            if failed && args.stop_on_error {
                halted = true;
                break;
            }
        }
    }

    Ok(reports)
}

fn resolve_slurm_partition(spec: &str) -> Result<String> {
    let spec = spec.trim();
    if spec.is_empty() {
        bail!("--slurm-partition cannot be empty");
    }

    if !is_partition_selection_spec(spec) {
        return Ok(spec.to_string());
    }

    let candidates = parse_partition_candidates(spec)?;
    select_best_slurm_partition(&candidates)
}

fn is_partition_selection_spec(spec: &str) -> bool {
    spec.eq_ignore_ascii_case("auto") || spec.contains(',') || spec.chars().any(char::is_whitespace)
}

fn parse_partition_candidates(spec: &str) -> Result<Vec<String>> {
    if spec.eq_ignore_ascii_case("auto") {
        return Ok(DEFAULT_SLURM_PARTITION_CANDIDATES
            .iter()
            .map(|partition| (*partition).to_string())
            .collect());
    }

    let candidates = spec
        .split(|ch: char| ch == ',' || ch.is_whitespace())
        .map(str::trim)
        .filter(|candidate| !candidate.is_empty())
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    if candidates.is_empty() {
        bail!("no Slurm partition candidates were provided in {spec:?}");
    }
    Ok(candidates)
}

fn select_best_slurm_partition(candidates: &[String]) -> Result<String> {
    let mut best_partition = None;
    let mut best_idle_nodes = 0_u64;
    let mut best_mixed_nodes = 0_u64;

    for candidate in candidates {
        let Some((idle_nodes, mixed_nodes)) = query_slurm_partition_capacity(candidate)? else {
            continue;
        };
        if idle_nodes == 0 && mixed_nodes == 0 {
            continue;
        }
        let better_partition = best_partition.is_none()
            || idle_nodes > best_idle_nodes
            || (idle_nodes == best_idle_nodes && mixed_nodes > best_mixed_nodes);
        if better_partition {
            best_partition = Some(candidate.clone());
            best_idle_nodes = idle_nodes;
            best_mixed_nodes = mixed_nodes;
        }
    }

    best_partition.ok_or_else(|| {
        anyhow::anyhow!(
            "no usable Slurm partition found among: {}",
            candidates.join(", ")
        )
    })
}

fn query_slurm_partition_capacity(partition: &str) -> Result<Option<(u64, u64)>> {
    let sinfo_args = vec![
        "-h".to_string(),
        "-p".to_string(),
        partition.to_string(),
        "-o".to_string(),
        "%a|%D|%t".to_string(),
    ];
    let mut last_retryable_error = None;

    for attempt in 1..=SLURM_QUERY_RETRIES {
        let output = run_command_with_timeout(
            "sinfo",
            &sinfo_args,
            SLURM_COMMAND_TIMEOUT,
            &format!("while querying Slurm partition {partition}"),
        )?;
        if output.status.success() {
            let stdout_text = String::from_utf8_lossy(&output.stdout);
            let mut idle_nodes = 0_u64;
            let mut mixed_nodes = 0_u64;
            for line in stdout_text.lines() {
                let mut fields = line.split('|');
                let availability = fields.next().unwrap_or_default().trim();
                let nodes = fields.next().unwrap_or_default().trim();
                let state = fields
                    .next()
                    .unwrap_or_default()
                    .trim()
                    .to_ascii_lowercase();
                if availability != "up" {
                    continue;
                }
                let node_count = nodes.parse::<u64>().unwrap_or(0);
                if state.starts_with("idle") {
                    idle_nodes += node_count;
                } else if state.starts_with("mix") {
                    mixed_nodes += node_count;
                }
            }
            return Ok(Some((idle_nodes, mixed_nodes)));
        }

        let stderr_text = String::from_utf8_lossy(&output.stderr).trim().to_string();
        if is_retryable_slurm_error(&stderr_text) {
            last_retryable_error = Some(stderr_text.clone());
            if attempt < SLURM_QUERY_RETRIES {
                thread::sleep(Duration::from_secs(attempt as u64 * 5));
                continue;
            }
        }
        if stderr_text.contains("Invalid partition") {
            return Ok(None);
        }
        bail!(
            "failed to query Slurm partition {partition}: {}",
            if stderr_text.is_empty() {
                format!("sinfo exited with status {}", output.status)
            } else {
                stderr_text
            }
        );
    }

    bail!(
        "failed to query Slurm partitions after retries: {}",
        last_retryable_error.unwrap_or_else(|| "unknown Slurm controller error".to_string())
    )
}

fn is_retryable_slurm_error(stderr_text: &str) -> bool {
    stderr_text.contains("Unable to contact slurm controller")
        || stderr_text.contains("connect failure")
}

fn run_command_with_timeout(
    program: &str,
    args: &[String],
    timeout: Duration,
    context: &str,
) -> Result<Output> {
    let mut child = Command::new(program)
        .args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .with_context(|| format!("failed to execute {program} {context}"))?;
    let started = Instant::now();

    loop {
        if child
            .try_wait()
            .with_context(|| format!("failed while waiting for {program} {context}"))?
            .is_some()
        {
            return child
                .wait_with_output()
                .with_context(|| format!("failed to capture output from {program} {context}"));
        }

        if started.elapsed() >= timeout {
            let _ = child.kill();
            let output = child.wait_with_output().with_context(|| {
                format!("failed to collect timed out output from {program} {context}")
            })?;
            let stderr_text = String::from_utf8_lossy(&output.stderr).trim().to_string();
            if stderr_text.is_empty() {
                bail!("{program} timed out after {}s {context}", timeout.as_secs());
            }
            bail!(
                "{program} timed out after {}s {context}: {stderr_text}",
                timeout.as_secs()
            );
        }

        thread::sleep(Duration::from_millis(200));
    }
}

fn build_plan(args: &BenchmarkPlanArgs) -> Result<BenchmarkPlan> {
    let samples = if let Some(manifest) = args.manifest.as_ref() {
        load_manifest_samples(manifest, args.sample_limit)?
    } else {
        vec![build_fixture_sample(args)?]
    };

    let mut cases = Vec::new();
    for sample in &samples {
        let sample_outdir = args.outdir.join(&sample.sample_id);
        cases.push(build_lineabact_case(args, sample, &sample_outdir)?);
        cases.push(build_spades_case(args, sample, &sample_outdir));
        cases.push(build_shovill_case(args, sample, &sample_outdir));
        cases.push(build_unicycler_case(args, sample, &sample_outdir));
    }

    Ok(BenchmarkPlan {
        schema_version: SCHEMA_VERSION,
        fixture_root: args.fixture_root.display().to_string(),
        samples,
        cases,
    })
}

fn build_fixture_sample(args: &BenchmarkPlanArgs) -> Result<BenchmarkSample> {
    let sample_r1 = require_file(
        &args
            .fixture_root
            .join("sample_data")
            .join("short_reads_1.fastq.gz"),
    )?;
    let sample_r2 = require_file(
        &args
            .fixture_root
            .join("sample_data")
            .join("short_reads_2.fastq.gz"),
    )?;
    let sample_reference = require_file(
        &args
            .fixture_root
            .join("sample_data")
            .join("reference.fasta"),
    )?;

    Ok(BenchmarkSample {
        sample_id: args.sample_id.clone(),
        r1: sample_r1.display().to_string(),
        r2: sample_r2.display().to_string(),
        reference: sample_reference.display().to_string(),
    })
}

fn load_manifest_samples(manifest: &Path, sample_limit: usize) -> Result<Vec<BenchmarkSample>> {
    let text = fs::read_to_string(manifest)
        .with_context(|| format!("failed to read benchmark manifest {}", manifest.display()))?;
    let mut lines = text.lines();
    let header_line = lines
        .next()
        .with_context(|| format!("benchmark manifest {} is empty", manifest.display()))?;
    let headers = header_line.split('\t').collect::<Vec<_>>();

    let mut samples = Vec::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        if let Some(sample) = parse_manifest_sample(&headers, &fields)? {
            samples.push(sample);
            if sample_limit > 0 && samples.len() >= sample_limit {
                break;
            }
        }
    }

    if samples.is_empty() {
        bail!(
            "benchmark manifest {} did not yield any valid paired short-read samples",
            manifest.display()
        );
    }

    Ok(samples)
}

fn parse_manifest_sample(headers: &[&str], fields: &[&str]) -> Result<Option<BenchmarkSample>> {
    let sample_id = first_nonempty_field(headers, fields, &["sample_id", "slug"])
        .map(ToString::to_string)
        .unwrap_or_default();
    let r1 = first_nonempty_field(headers, fields, &["r1", "illumina_fastq_1"])
        .map(ToString::to_string)
        .unwrap_or_default();
    let r2 = first_nonempty_field(headers, fields, &["r2", "illumina_fastq_2"])
        .map(ToString::to_string)
        .unwrap_or_default();
    let reference = first_nonempty_field(headers, fields, &["reference", "reference_fasta"])
        .map(ToString::to_string)
        .unwrap_or_default();

    if sample_id.is_empty() || r1.is_empty() || r2.is_empty() || reference.is_empty() {
        return Ok(None);
    }

    let r1_path = Path::new(&r1);
    let r2_path = Path::new(&r2);
    let reference_path = Path::new(&reference);
    if !r1_path.exists() || !r2_path.exists() || !reference_path.exists() {
        return Ok(None);
    }

    Ok(Some(BenchmarkSample {
        sample_id,
        r1,
        r2,
        reference,
    }))
}

fn first_nonempty_field<'a>(
    headers: &[&str],
    fields: &'a [&str],
    candidates: &[&str],
) -> Option<&'a str> {
    for candidate in candidates {
        if let Some(index) = headers.iter().position(|header| header == candidate)
            && let Some(value) = fields.get(index)
        {
            let trimmed = value.trim();
            if !trimmed.is_empty() {
                return Some(trimmed);
            }
        }
    }
    None
}

fn build_lineabact_case(
    args: &BenchmarkPlanArgs,
    sample: &BenchmarkSample,
    sample_outdir: &Path,
) -> Result<BenchmarkCase> {
    let outdir = sample_outdir.join("lineabact");
    let genome_size_bp = reference_genome_size_bp(Path::new(&sample.reference))?;
    let command = format!(
        "{exe} assemble --backend spades --r1 {r1} --r2 {r2} --outdir {outdir} --threads {threads} --k {k} --spades-executable {spades} --trim-adapters --target-coverage 150 --genome-size-bp {genome_size_bp}",
        exe = shell_quote(&args.lineabact_executable),
        r1 = shell_quote(&sample.r1),
        r2 = shell_quote(&sample.r2),
        outdir = shell_quote(&outdir.display().to_string()),
        threads = args.threads,
        k = args.k,
        spades = shell_quote(&args.spades_executable),
        genome_size_bp = genome_size_bp,
    );

    Ok(BenchmarkCase {
        sample_id: sample.sample_id.clone(),
        tool: "lineabact".to_string(),
        executable: args.lineabact_executable.clone(),
        command,
        outdir: outdir.display().to_string(),
        expected_dirs: contract::STABLE_OUTPUT_DIRS
            .iter()
            .map(|value| value.to_string())
            .collect(),
        expected_files: lineabact_expected_files(),
        notes: "Rust front-end baseline: deterministic read optimization, adapter trimming, coverage capping, SPAdes backend, contig postprocess, and Unicycler-like graph/bridging export.".to_string(),
    })
}

fn build_spades_case(
    args: &BenchmarkPlanArgs,
    sample: &BenchmarkSample,
    sample_outdir: &Path,
) -> BenchmarkCase {
    let outdir = sample_outdir.join("spades");
    let command = format!(
        "{exe} --isolate -1 {r1} -2 {r2} -o {outdir} -t {threads} -k {k}",
        exe = shell_quote(&args.spades_executable),
        r1 = shell_quote(&sample.r1),
        r2 = shell_quote(&sample.r2),
        outdir = shell_quote(&outdir.display().to_string()),
        threads = args.threads,
        k = args.k,
    );

    BenchmarkCase {
        sample_id: sample.sample_id.clone(),
        tool: "spades".to_string(),
        executable: args.spades_executable.clone(),
        command,
        outdir: outdir.display().to_string(),
        expected_dirs: Vec::new(),
        expected_files: contract::SPADES_OUTPUT_FILES
            .iter()
            .filter(|file| {
                !matches!(
                    **file,
                    contract::SPADES_COMMAND_SH | contract::SPADES_PLAN_JSON
                )
            })
            .map(|file| file.to_string())
            .collect(),
        notes: "SPAdes backbone assembly and graph baseline.".to_string(),
    }
}

fn build_shovill_case(
    args: &BenchmarkPlanArgs,
    sample: &BenchmarkSample,
    sample_outdir: &Path,
) -> BenchmarkCase {
    let outdir = sample_outdir.join("shovill");
    let command = format!(
        "{exe} --R1 {r1} --R2 {r2} --outdir {outdir} --cpus {threads} --depth 150 --assembler spades --force",
        exe = shell_quote(&args.shovill_executable),
        r1 = shell_quote(&sample.r1),
        r2 = shell_quote(&sample.r2),
        outdir = shell_quote(&outdir.display().to_string()),
        threads = args.threads,
    );

    BenchmarkCase {
        sample_id: sample.sample_id.clone(),
        tool: "shovill".to_string(),
        executable: args.shovill_executable.clone(),
        command,
        outdir: outdir.display().to_string(),
        expected_dirs: Vec::new(),
        expected_files: vec![
            "contigs.fa".to_string(),
            "shovill.log".to_string(),
            "shovill.corrections".to_string(),
        ],
        notes: "Shovill-style read thinning, depth normalization and assembly workflow baseline."
            .to_string(),
    }
}

fn build_unicycler_case(
    args: &BenchmarkPlanArgs,
    sample: &BenchmarkSample,
    sample_outdir: &Path,
) -> BenchmarkCase {
    let outdir = sample_outdir.join("unicycler");
    let command = format!(
        "{exe} -1 {r1} -2 {r2} -o {outdir} --threads {threads} --keep 1 --spades_path {spades}",
        exe = shell_quote(&args.unicycler_executable),
        r1 = shell_quote(&sample.r1),
        r2 = shell_quote(&sample.r2),
        outdir = shell_quote(&outdir.display().to_string()),
        threads = args.threads,
        spades = shell_quote(&args.spades_executable),
    );

    BenchmarkCase {
        sample_id: sample.sample_id.clone(),
        tool: "unicycler".to_string(),
        executable: args.unicycler_executable.clone(),
        command,
        outdir: outdir.display().to_string(),
        expected_dirs: vec!["spades_assembly".to_string()],
        expected_files: vec![
            "assembly.gfa".to_string(),
            "assembly.fasta".to_string(),
            "unicycler.log".to_string(),
        ],
        notes: "Unicycler short-read-first SPAdes-optimising baseline.".to_string(),
    }
}

fn execute_case_locally(
    case: &BenchmarkCase,
    repeat_index: usize,
    repeat_count: usize,
) -> Result<BenchmarkCaseReport> {
    fs::create_dir_all(&case.outdir)?;
    let stdout_log = Path::new(&case.outdir).join("benchmark.stdout.log");
    let stderr_log = Path::new(&case.outdir).join("benchmark.stderr.log");
    let started = Instant::now();

    let output = Command::new("bash")
        .arg("-lc")
        .arg(&case.command)
        .output()
        .with_context(|| format!("failed to execute benchmark case command: {}", case.command));

    let elapsed_seconds = started.elapsed().as_secs_f64();
    let mut exit_code = None;
    let mut error = None;

    match output {
        Ok(output) => {
            exit_code = output.status.code();
            io::ensure_parent_dir(&stdout_log)?;
            fs::write(&stdout_log, &output.stdout)?;
            fs::write(&stderr_log, &output.stderr)?;
        }
        Err(err) => {
            io::ensure_parent_dir(&stdout_log)?;
            fs::write(&stdout_log, b"")?;
            fs::write(&stderr_log, format!("{err:#}\n"))?;
            error = Some(format!("{err:#}"));
        }
    }

    let missing_dirs = case
        .expected_dirs
        .iter()
        .filter(|dir| !Path::new(&case.outdir).join(dir).is_dir())
        .cloned()
        .collect::<Vec<_>>();
    let missing_files = case
        .expected_files
        .iter()
        .filter(|file| !Path::new(&case.outdir).join(file).is_file())
        .cloned()
        .collect::<Vec<_>>();

    let status = if error.is_some() {
        "error".to_string()
    } else if exit_code == Some(0) && missing_dirs.is_empty() && missing_files.is_empty() {
        "ok".to_string()
    } else if exit_code != Some(0) {
        "command_failed".to_string()
    } else {
        "missing_outputs".to_string()
    };

    Ok(BenchmarkCaseReport {
        sample_id: case.sample_id.clone(),
        repeat_index,
        repeat_count,
        scheduler: "local".to_string(),
        tool: case.tool.clone(),
        executable: case.executable.clone(),
        command: case.command.clone(),
        outdir: case.outdir.clone(),
        expected_dirs: case.expected_dirs.clone(),
        expected_files: case.expected_files.clone(),
        status,
        exit_code,
        elapsed_seconds,
        elapsed_seconds_source: "process_wallclock".to_string(),
        missing_dirs,
        missing_files,
        stdout_log: stdout_log.display().to_string(),
        stderr_log: stderr_log.display().to_string(),
        submission_script: None,
        submission_command: None,
        job_id: None,
        quality_metrics: None,
        error,
    })
}

fn submit_case_to_slurm(
    case: &BenchmarkCase,
    args: &BenchmarkRunArgs,
    repeat_index: usize,
    repeat_count: usize,
) -> Result<BenchmarkCaseReport> {
    fs::create_dir_all(&case.outdir)?;
    let slurm_dir = Path::new(&case.outdir).join("slurm");
    fs::create_dir_all(&slurm_dir)?;

    let stdout_log = slurm_dir.join("benchmark.stdout.log");
    let stderr_log = slurm_dir.join("benchmark.stderr.log");
    let submission_script = slurm_dir.join("benchmark.sbatch.sh");
    let repo_root = env::current_dir().context("failed to resolve current repository root")?;
    write_slurm_script(&submission_script, &repo_root, case, args)?;

    let submission_command =
        build_sbatch_command(case, args, &stdout_log, &stderr_log, &submission_script);
    if args.slurm_dry_run {
        return Ok(BenchmarkCaseReport {
            sample_id: case.sample_id.clone(),
            repeat_index,
            repeat_count,
            scheduler: "slurm".to_string(),
            tool: case.tool.clone(),
            executable: case.executable.clone(),
            command: case.command.clone(),
            outdir: case.outdir.clone(),
            expected_dirs: case.expected_dirs.clone(),
            expected_files: case.expected_files.clone(),
            status: "scripted".to_string(),
            exit_code: None,
            elapsed_seconds: 0.0,
            elapsed_seconds_source: "none".to_string(),
            missing_dirs: Vec::new(),
            missing_files: Vec::new(),
            stdout_log: stdout_log.display().to_string(),
            stderr_log: stderr_log.display().to_string(),
            submission_script: Some(submission_script.display().to_string()),
            submission_command: Some(submission_command),
            job_id: None,
            quality_metrics: None,
            error: None,
        });
    }

    let started = Instant::now();
    let output =
        submit_sbatch_with_retries(case, args, &stdout_log, &stderr_log, &submission_script)
            .with_context(|| format!("failed to submit Slurm benchmark case {}", case.command))?;
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let stdout_text = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let stderr_text = String::from_utf8_lossy(&output.stderr).trim().to_string();

    Ok(BenchmarkCaseReport {
        sample_id: case.sample_id.clone(),
        repeat_index,
        repeat_count,
        scheduler: "slurm".to_string(),
        tool: case.tool.clone(),
        executable: case.executable.clone(),
        command: case.command.clone(),
        outdir: case.outdir.clone(),
        expected_dirs: case.expected_dirs.clone(),
        expected_files: case.expected_files.clone(),
        status: if output.status.success() {
            "submitted".to_string()
        } else {
            "submission_failed".to_string()
        },
        exit_code: output.status.code(),
        elapsed_seconds,
        elapsed_seconds_source: "sbatch_submission".to_string(),
        missing_dirs: Vec::new(),
        missing_files: Vec::new(),
        stdout_log: stdout_log.display().to_string(),
        stderr_log: stderr_log.display().to_string(),
        submission_script: Some(submission_script.display().to_string()),
        submission_command: Some(submission_command),
        job_id: if output.status.success() {
            parse_sbatch_job_id(&stdout_text)
        } else {
            None
        },
        quality_metrics: None,
        error: if output.status.success() {
            None
        } else if !stderr_text.is_empty() {
            Some(stderr_text)
        } else {
            Some(format!("sbatch exited with status {}", output.status))
        },
    })
}

fn submit_sbatch_with_retries(
    case: &BenchmarkCase,
    args: &BenchmarkRunArgs,
    stdout_log: &Path,
    stderr_log: &Path,
    submission_script: &Path,
) -> Result<std::process::Output> {
    let sbatch_args = build_sbatch_args(case, args, stdout_log, stderr_log, submission_script);
    let mut last_output = None;

    for attempt in 1..=SLURM_SUBMISSION_RETRIES {
        let output = run_command_with_timeout(
            "sbatch",
            &sbatch_args,
            SLURM_COMMAND_TIMEOUT,
            &format!("for benchmark case {} (attempt {attempt})", case.command),
        )?;
        if output.status.success() {
            return Ok(output);
        }

        let stderr_text = String::from_utf8_lossy(&output.stderr);
        let retryable = is_retryable_slurm_error(&stderr_text);
        last_output = Some(output);
        if !retryable || attempt == SLURM_SUBMISSION_RETRIES {
            break;
        }
        thread::sleep(Duration::from_secs(attempt as u64 * 5));
    }

    match last_output {
        Some(output) => Ok(output),
        None => bail!("sbatch did not produce any submission attempts"),
    }
}

fn skipped_case_report(
    case: &BenchmarkCase,
    scheduler: &BenchmarkScheduler,
    repeat_index: usize,
    repeat_count: usize,
) -> BenchmarkCaseReport {
    let stdout_log = Path::new(&case.outdir).join("benchmark.stdout.log");
    let stderr_log = Path::new(&case.outdir).join("benchmark.stderr.log");
    BenchmarkCaseReport {
        sample_id: case.sample_id.clone(),
        repeat_index,
        repeat_count,
        scheduler: match scheduler {
            BenchmarkScheduler::Local => "local".to_string(),
            BenchmarkScheduler::Slurm => "slurm".to_string(),
        },
        tool: case.tool.clone(),
        executable: case.executable.clone(),
        command: case.command.clone(),
        outdir: case.outdir.clone(),
        expected_dirs: case.expected_dirs.clone(),
        expected_files: case.expected_files.clone(),
        status: "skipped".to_string(),
        exit_code: None,
        elapsed_seconds: 0.0,
        elapsed_seconds_source: "none".to_string(),
        missing_dirs: Vec::new(),
        missing_files: Vec::new(),
        stdout_log: stdout_log.display().to_string(),
        stderr_log: stderr_log.display().to_string(),
        submission_script: None,
        submission_command: None,
        job_id: None,
        quality_metrics: None,
        error: None,
    }
}

fn repeated_case(case: &BenchmarkCase, repeat_index: usize) -> BenchmarkCase {
    let repeated_outdir = Path::new(&case.outdir)
        .join("repeats")
        .join(format!("run_{repeat_index:03}"));
    let original_outdir = shell_quote(&case.outdir);
    let repeated_outdir_text = repeated_outdir.display().to_string();
    let repeated_outdir_quoted = shell_quote(&repeated_outdir_text);
    let command = if case.command.contains(&original_outdir) {
        case.command
            .replace(&original_outdir, &repeated_outdir_quoted)
    } else {
        case.command.replace(&case.outdir, &repeated_outdir_text)
    };

    BenchmarkCase {
        sample_id: case.sample_id.clone(),
        tool: case.tool.clone(),
        executable: case.executable.clone(),
        command,
        outdir: repeated_outdir_text,
        expected_dirs: case.expected_dirs.clone(),
        expected_files: case.expected_files.clone(),
        notes: case.notes.clone(),
    }
}

fn write_slurm_script(
    script_path: &Path,
    repo_root: &Path,
    case: &BenchmarkCase,
    args: &BenchmarkRunArgs,
) -> Result<()> {
    io::ensure_parent_dir(script_path)?;
    let file = fs::File::create(script_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "#!/usr/bin/env bash")?;
    writeln!(writer, "set -euo pipefail")?;
    write_slurm_env_bootstrap(&mut writer, args)?;
    writeln!(
        writer,
        "cd {}",
        shell_quote(&repo_root.display().to_string())
    )?;
    writeln!(writer, "exec {}", case.command)?;
    writer.flush()?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(script_path)?.permissions();
        perms.set_mode(0o755);
        fs::set_permissions(script_path, perms)?;
    }
    Ok(())
}

fn write_slurm_env_bootstrap<W: Write>(writer: &mut W, args: &BenchmarkRunArgs) -> Result<()> {
    let conda_base = Path::new(&args.slurm_conda_base);
    let conda_sh = conda_base.join("etc").join("profile.d").join("conda.sh");
    if conda_sh.exists() {
        writeln!(
            writer,
            "source {}",
            shell_quote(&conda_sh.display().to_string())
        )?;
        writeln!(
            writer,
            "conda activate {}",
            shell_quote(&args.slurm_conda_env)
        )?;
        return Ok(());
    }

    let env_prefix = resolve_slurm_env_prefix(args).with_context(|| {
        format!(
            "could not resolve Slurm conda environment from base={} env={}",
            args.slurm_conda_base, args.slurm_conda_env
        )
    })?;
    let env_bin = env_prefix.join("bin");
    writeln!(
        writer,
        "export CONDA_PREFIX={}",
        shell_quote(&env_prefix.display().to_string())
    )?;
    writeln!(
        writer,
        "export PATH={}:$PATH",
        shell_quote(&env_bin.display().to_string())
    )?;
    Ok(())
}

fn resolve_slurm_env_prefix(args: &BenchmarkRunArgs) -> Result<PathBuf> {
    let env_path = Path::new(&args.slurm_conda_env);
    if env_path.is_absolute() && env_path.join("bin").is_dir() {
        return Ok(env_path.to_path_buf());
    }

    let base_path = Path::new(&args.slurm_conda_base);
    if base_path.join("bin").is_dir() {
        return Ok(base_path.to_path_buf());
    }

    let named_env = base_path.join("envs").join(&args.slurm_conda_env);
    if named_env.join("bin").is_dir() {
        return Ok(named_env);
    }

    bail!(
        "neither {} nor {} resolves to a conda environment prefix",
        args.slurm_conda_base,
        args.slurm_conda_env
    );
}

fn build_sbatch_args(
    case: &BenchmarkCase,
    args: &BenchmarkRunArgs,
    stdout_log: &Path,
    stderr_log: &Path,
    script_path: &Path,
) -> Vec<String> {
    let mut sbatch_args = vec![
        "--parsable".to_string(),
        "--partition".to_string(),
        args.slurm_partition.clone(),
        "--cpus-per-task".to_string(),
        args.slurm_cpus_per_task.to_string(),
        "--time".to_string(),
        args.slurm_time.clone(),
        "--job-name".to_string(),
        slurm_job_name(case),
        "--output".to_string(),
        stdout_log.display().to_string(),
        "--error".to_string(),
        stderr_log.display().to_string(),
    ];
    if let Some(mem_gb) = args.slurm_mem_gb {
        sbatch_args.push("--mem".to_string());
        sbatch_args.push(format!("{mem_gb}G"));
    }
    sbatch_args.push(script_path.display().to_string());
    sbatch_args
}

fn build_sbatch_command(
    case: &BenchmarkCase,
    args: &BenchmarkRunArgs,
    stdout_log: &Path,
    stderr_log: &Path,
    script_path: &Path,
) -> String {
    let args = build_sbatch_args(case, args, stdout_log, stderr_log, script_path);
    let mut command = String::from("sbatch");
    for arg in args {
        command.push(' ');
        command.push_str(&shell_quote(&arg));
    }
    command
}

fn slurm_job_name(case: &BenchmarkCase) -> String {
    let sample = sanitize_slurm_token(&case.sample_id);
    let tool = sanitize_slurm_token(&case.tool);
    format!("lb_{sample}_{tool}")
}

fn sanitize_slurm_token(text: &str) -> String {
    let sanitized = text
        .chars()
        .map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' })
        .collect::<String>();
    sanitized.chars().take(64).collect()
}

fn parse_sbatch_job_id(stdout_text: &str) -> Option<String> {
    let first_line = stdout_text.lines().next()?.trim();
    if first_line.is_empty() {
        None
    } else {
        Some(first_line.to_string())
    }
}

fn write_plan_tsv(plan: &BenchmarkPlan, path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "sample_id\ttool\texecutable\toutdir\tcommand\texpected_dirs\texpected_files\tnotes"
    )?;
    for case in &plan.cases {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            case.sample_id,
            case.tool,
            case.executable,
            case.outdir,
            case.command,
            case.expected_dirs.join(","),
            case.expected_files.join(","),
            case.notes
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_report_tsv(report: &BenchmarkReport, path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "sample_id\trepeat_index\trepeat_count\tscheduler\ttool\texecutable\toutdir\tstatus\texit_code\tjob_id\telapsed_seconds\telapsed_seconds_source\tmissing_dirs\tmissing_files\texpected_dirs\texpected_files\tstdout_log\tstderr_log\tsubmission_script\tsubmission_command\tcommand"
    )?;
    for case in &report.cases {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            case.sample_id,
            case.repeat_index,
            case.repeat_count,
            case.scheduler,
            case.tool,
            case.executable,
            case.outdir,
            case.status,
            case.exit_code
                .map(|value| value.to_string())
                .unwrap_or_default(),
            case.job_id.clone().unwrap_or_default(),
            case.elapsed_seconds,
            case.elapsed_seconds_source,
            case.missing_dirs.join(","),
            case.missing_files.join(","),
            case.expected_dirs.join(","),
            case.expected_files.join(","),
            case.stdout_log,
            case.stderr_log,
            case.submission_script.clone().unwrap_or_default(),
            case.submission_command.clone().unwrap_or_default(),
            case.command
        )?;
    }
    writer.flush()?;

    let runtime_tsv = path.with_file_name("benchmark_runtime_comparison.tsv");
    write_runtime_comparison_tsv(&report.runtime_comparisons, &runtime_tsv)?;
    let quality_tsv = path.with_file_name("benchmark_quality_comparison.tsv");
    write_quality_comparison_tsv(&report.quality_comparisons, &quality_tsv)?;
    let acceptance_tsv = path.with_file_name("benchmark_sample_acceptance.tsv");
    write_sample_acceptance_tsv(&report.sample_acceptance, &acceptance_tsv)?;
    Ok(())
}

fn build_runtime_comparisons(cases: &[BenchmarkCaseReport]) -> Vec<RuntimeComparison> {
    let mut comparisons = Vec::new();
    let sample_ids = cases
        .iter()
        .map(|case| case.sample_id.as_str())
        .collect::<Vec<_>>();

    for sample_id in sample_ids {
        let candidate_times = collect_ok_elapsed_seconds(cases, sample_id, "lineabact");
        let baseline_times = collect_ok_elapsed_seconds(cases, sample_id, "shovill");
        if candidate_times.is_empty() || baseline_times.is_empty() {
            continue;
        }

        let Some(candidate_median) = median(&candidate_times) else {
            continue;
        };
        let Some(baseline_median) = median(&baseline_times) else {
            continue;
        };
        if baseline_median <= 0.0 || candidate_median <= 0.0 {
            continue;
        }

        let runtime_ratio = candidate_median / baseline_median;
        let speedup = baseline_median / candidate_median;
        comparisons.push(RuntimeComparison {
            sample_id: sample_id.to_string(),
            baseline_tool: "shovill".to_string(),
            candidate_tool: "lineabact".to_string(),
            baseline_repeat_count: baseline_times.len(),
            candidate_repeat_count: candidate_times.len(),
            baseline_median_elapsed_seconds: baseline_median,
            candidate_median_elapsed_seconds: candidate_median,
            median_runtime_ratio_vs_baseline: runtime_ratio,
            median_speedup_vs_baseline: speedup,
            target_runtime_ratio: SHOVILL_RUNTIME_TARGET_RATIO,
            meets_target: runtime_ratio <= SHOVILL_RUNTIME_TARGET_RATIO,
        });
    }

    comparisons.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    comparisons.dedup_by(|a, b| a.sample_id == b.sample_id);
    comparisons
}

fn build_quality_comparisons(
    samples: &[BenchmarkSample],
    cases: &[BenchmarkCaseReport],
) -> Vec<QualityComparison> {
    let mut comparisons = Vec::new();
    for sample in samples {
        let candidate_metrics = collect_ok_quality_metrics(cases, &sample.sample_id, "lineabact");
        let baseline_metrics = collect_ok_quality_metrics(cases, &sample.sample_id, "shovill");
        if candidate_metrics.is_empty() || baseline_metrics.is_empty() {
            continue;
        }

        let candidate_n50 = median_usize_field(&candidate_metrics, |metrics| metrics.n50);
        let baseline_n50 = median_usize_field(&baseline_metrics, |metrics| metrics.n50);
        let candidate_genome_fraction =
            median_f64_field(&candidate_metrics, |metrics| metrics.genome_fraction);
        let baseline_genome_fraction =
            median_f64_field(&baseline_metrics, |metrics| metrics.genome_fraction);
        let candidate_misassemblies =
            median_usize_field(&candidate_metrics, |metrics| metrics.misassembly_count);
        let baseline_misassemblies =
            median_usize_field(&baseline_metrics, |metrics| metrics.misassembly_count);
        let candidate_mismatch_rate =
            median_f64_field(&candidate_metrics, |metrics| metrics.mismatch_rate);
        let baseline_mismatch_rate =
            median_f64_field(&baseline_metrics, |metrics| metrics.mismatch_rate);
        let candidate_indel_rate =
            median_f64_field(&candidate_metrics, |metrics| metrics.indel_rate);
        let baseline_indel_rate = median_f64_field(&baseline_metrics, |metrics| metrics.indel_rate);

        if baseline_n50 <= 0.0 {
            continue;
        }

        let n50_ratio = candidate_n50 / baseline_n50;
        let passes_n50 = n50_ratio >= QUALITY_N50_MIN_RATIO;
        let passes_genome_fraction = candidate_genome_fraction + 1e-12 >= baseline_genome_fraction;
        let passes_misassemblies = candidate_misassemblies <= baseline_misassemblies + 1e-12;
        let passes_mismatch_rate = candidate_mismatch_rate <= baseline_mismatch_rate + 1e-12;
        let passes_indel_rate = candidate_indel_rate <= baseline_indel_rate + 1e-12;
        let meets_target = passes_n50
            && passes_genome_fraction
            && passes_misassemblies
            && passes_mismatch_rate
            && passes_indel_rate;

        comparisons.push(QualityComparison {
            sample_id: sample.sample_id.clone(),
            baseline_tool: "shovill".to_string(),
            candidate_tool: "lineabact".to_string(),
            baseline_repeat_count: baseline_metrics.len(),
            candidate_repeat_count: candidate_metrics.len(),
            baseline_median_n50: baseline_n50,
            candidate_median_n50: candidate_n50,
            n50_ratio_vs_baseline: n50_ratio,
            target_n50_ratio: QUALITY_N50_MIN_RATIO,
            preferred_n50_ratio: QUALITY_N50_PREFERRED_RATIO,
            baseline_median_genome_fraction: baseline_genome_fraction,
            candidate_median_genome_fraction: candidate_genome_fraction,
            baseline_median_misassemblies: baseline_misassemblies,
            candidate_median_misassemblies: candidate_misassemblies,
            baseline_median_mismatch_rate: baseline_mismatch_rate,
            candidate_median_mismatch_rate: candidate_mismatch_rate,
            baseline_median_indel_rate: baseline_indel_rate,
            candidate_median_indel_rate: candidate_indel_rate,
            passes_n50,
            passes_genome_fraction,
            passes_misassemblies,
            passes_mismatch_rate,
            passes_indel_rate,
            meets_target,
        });
    }

    comparisons.sort_by(|left, right| left.sample_id.cmp(&right.sample_id));
    comparisons
}

fn build_sample_acceptance(
    runtime_comparisons: &[RuntimeComparison],
    quality_comparisons: &[QualityComparison],
) -> Vec<SampleAcceptance> {
    let mut sample_ids = runtime_comparisons
        .iter()
        .map(|comparison| comparison.sample_id.clone())
        .collect::<Vec<_>>();
    for comparison in quality_comparisons {
        if !sample_ids
            .iter()
            .any(|sample_id| sample_id == &comparison.sample_id)
        {
            sample_ids.push(comparison.sample_id.clone());
        }
    }
    sample_ids.sort();
    sample_ids.dedup();

    let runtime_map = runtime_comparisons
        .iter()
        .map(|comparison| (comparison.sample_id.as_str(), comparison))
        .collect::<HashMap<_, _>>();
    let quality_map = quality_comparisons
        .iter()
        .map(|comparison| (comparison.sample_id.as_str(), comparison))
        .collect::<HashMap<_, _>>();

    sample_ids
        .into_iter()
        .map(|sample_id| {
            let runtime = runtime_map.get(sample_id.as_str()).copied();
            let quality = quality_map.get(sample_id.as_str()).copied();
            let runtime_meets_target = runtime.is_some_and(|value| value.meets_target);
            let quality_meets_target = quality.is_some_and(|value| value.meets_target);
            let accepted = runtime_meets_target && quality_meets_target;
            let mut failure_reasons = Vec::new();
            if runtime.is_none() {
                failure_reasons.push("missing_runtime_comparison".to_string());
            } else if !runtime_meets_target {
                failure_reasons.push("runtime_ratio_above_target".to_string());
            }
            if quality.is_none() {
                failure_reasons.push("missing_quality_comparison".to_string());
            } else if let Some(quality) = quality {
                if !quality.passes_n50 {
                    failure_reasons.push("n50_below_threshold".to_string());
                }
                if !quality.passes_genome_fraction {
                    failure_reasons.push("genome_fraction_below_shovill".to_string());
                }
                if !quality.passes_misassemblies {
                    failure_reasons.push("misassemblies_above_shovill".to_string());
                }
                if !quality.passes_mismatch_rate {
                    failure_reasons.push("mismatch_rate_worse_than_shovill".to_string());
                }
                if !quality.passes_indel_rate {
                    failure_reasons.push("indel_rate_worse_than_shovill".to_string());
                }
            }

            SampleAcceptance {
                sample_id,
                runtime_meets_target,
                quality_meets_target,
                accepted,
                failure_reasons,
            }
        })
        .collect()
}

fn write_runtime_comparison_tsv(comparisons: &[RuntimeComparison], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "sample_id\tbaseline_tool\tcandidate_tool\tbaseline_repeat_count\tcandidate_repeat_count\tbaseline_median_elapsed_seconds\tcandidate_median_elapsed_seconds\tmedian_runtime_ratio_vs_baseline\tmedian_speedup_vs_baseline\ttarget_runtime_ratio\tmeets_target"
    )?;
    for comparison in comparisons {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}",
            comparison.sample_id,
            comparison.baseline_tool,
            comparison.candidate_tool,
            comparison.baseline_repeat_count,
            comparison.candidate_repeat_count,
            comparison.baseline_median_elapsed_seconds,
            comparison.candidate_median_elapsed_seconds,
            comparison.median_runtime_ratio_vs_baseline,
            comparison.median_speedup_vs_baseline,
            comparison.target_runtime_ratio,
            comparison.meets_target
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_quality_comparison_tsv(comparisons: &[QualityComparison], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "sample_id\tbaseline_tool\tcandidate_tool\tbaseline_repeat_count\tcandidate_repeat_count\tbaseline_median_n50\tcandidate_median_n50\tn50_ratio_vs_baseline\ttarget_n50_ratio\tpreferred_n50_ratio\tbaseline_median_genome_fraction\tcandidate_median_genome_fraction\tbaseline_median_misassemblies\tcandidate_median_misassemblies\tbaseline_median_mismatch_rate\tcandidate_median_mismatch_rate\tbaseline_median_indel_rate\tcandidate_median_indel_rate\tpasses_n50\tpasses_genome_fraction\tpasses_misassemblies\tpasses_mismatch_rate\tpasses_indel_rate\tmeets_target"
    )?;
    for comparison in comparisons {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.6}\t{:.6}\t{:.4}\t{:.4}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{}\t{}\t{}\t{}\t{}\t{}",
            comparison.sample_id,
            comparison.baseline_tool,
            comparison.candidate_tool,
            comparison.baseline_repeat_count,
            comparison.candidate_repeat_count,
            comparison.baseline_median_n50,
            comparison.candidate_median_n50,
            comparison.n50_ratio_vs_baseline,
            comparison.target_n50_ratio,
            comparison.preferred_n50_ratio,
            comparison.baseline_median_genome_fraction,
            comparison.candidate_median_genome_fraction,
            comparison.baseline_median_misassemblies,
            comparison.candidate_median_misassemblies,
            comparison.baseline_median_mismatch_rate,
            comparison.candidate_median_mismatch_rate,
            comparison.baseline_median_indel_rate,
            comparison.candidate_median_indel_rate,
            comparison.passes_n50,
            comparison.passes_genome_fraction,
            comparison.passes_misassemblies,
            comparison.passes_mismatch_rate,
            comparison.passes_indel_rate,
            comparison.meets_target
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_sample_acceptance_tsv(rows: &[SampleAcceptance], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "sample_id\truntime_meets_target\tquality_meets_target\taccepted\tfailure_reasons"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            row.sample_id,
            row.runtime_meets_target,
            row.quality_meets_target,
            row.accepted,
            row.failure_reasons.join(",")
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn reference_genome_size_bp(reference: &Path) -> Result<u64> {
    let mut reader = parse_fastx_file(reference)
        .with_context(|| format!("failed to open reference FASTA {}", reference.display()))?;
    let mut total_bases = 0u64;
    while let Some(record) = reader.next() {
        let record = record?;
        total_bases += record.seq().len() as u64;
    }
    if total_bases == 0 {
        bail!(
            "reference FASTA {} contains no sequence",
            reference.display()
        );
    }
    Ok(total_bases)
}

fn hydrate_case_reports(samples: &[BenchmarkSample], cases: &mut [BenchmarkCaseReport]) {
    let sample_map = samples
        .iter()
        .map(|sample| (sample.sample_id.as_str(), sample))
        .collect::<HashMap<_, _>>();
    for case in cases {
        hydrate_case_report(&sample_map, case);
    }
}

fn hydrate_case_report(
    sample_map: &HashMap<&str, &BenchmarkSample>,
    case: &mut BenchmarkCaseReport,
) {
    normalize_case_expectations(case);
    let (missing_dirs, missing_files) = collect_missing_outputs(case);
    case.missing_dirs = missing_dirs;
    case.missing_files = missing_files;

    if let Some((elapsed_seconds, source)) = recover_elapsed_seconds(case) {
        case.elapsed_seconds = elapsed_seconds;
        case.elapsed_seconds_source = source.to_string();
    }

    if case.error.is_some() {
        return;
    }
    if !case.missing_dirs.is_empty() || !case.missing_files.is_empty() {
        if case.status == "ok" {
            case.status = "missing_outputs".to_string();
        }
        return;
    }
    if case.elapsed_seconds > 0.0 && case.elapsed_seconds_source != "sbatch_submission" {
        case.status = "ok".to_string();
    }
    case.quality_metrics = recover_case_quality(sample_map, case);
}

fn normalize_case_expectations(case: &mut BenchmarkCaseReport) {
    if case.tool == "lineabact" {
        case.expected_dirs = contract::STABLE_OUTPUT_DIRS
            .iter()
            .map(|value| value.to_string())
            .collect();
        case.expected_files = lineabact_expected_files();
    }
}

fn collect_missing_outputs(case: &BenchmarkCaseReport) -> (Vec<String>, Vec<String>) {
    let missing_dirs = case
        .expected_dirs
        .iter()
        .filter(|dir| !Path::new(&case.outdir).join(dir).is_dir())
        .cloned()
        .collect::<Vec<_>>();
    let missing_files = case
        .expected_files
        .iter()
        .filter(|file| !Path::new(&case.outdir).join(file).is_file())
        .cloned()
        .collect::<Vec<_>>();
    (missing_dirs, missing_files)
}

fn recover_elapsed_seconds(case: &BenchmarkCaseReport) -> Option<(f64, &'static str)> {
    match case.tool.as_str() {
        "lineabact" => recover_lineabact_elapsed_seconds(Path::new(&case.outdir)),
        "shovill" => recover_shovill_elapsed_seconds(Path::new(&case.stderr_log)),
        _ => None,
    }
}

fn recover_case_quality(
    sample_map: &HashMap<&str, &BenchmarkSample>,
    case: &BenchmarkCaseReport,
) -> Option<CaseQualityMetrics> {
    if case.status != "ok" || !case.missing_dirs.is_empty() || !case.missing_files.is_empty() {
        return None;
    }
    if !matches!(case.tool.as_str(), "lineabact" | "shovill") {
        return None;
    }

    let sample = sample_map.get(case.sample_id.as_str())?;
    let (assembly_path, assembly_source) = resolve_case_assembly_path(case)?;
    evaluate_reference_quality(
        Path::new(&sample.reference),
        &assembly_path,
        assembly_source,
    )
    .ok()
}

fn resolve_case_assembly_path(case: &BenchmarkCaseReport) -> Option<(PathBuf, &'static str)> {
    let outdir = Path::new(&case.outdir);
    match case.tool.as_str() {
        "lineabact" => Some((
            outdir.join(contract::ROOT_FINISHED_CONTIGS),
            "finished_contigs",
        )),
        "shovill" => Some((outdir.join("contigs.fa"), "contigs.fa")),
        "spades" => Some((outdir.join(contract::SPADES_CONTIGS_FASTA), "contigs.fasta")),
        "unicycler" => Some((outdir.join("assembly.fasta"), "assembly.fasta")),
        _ => None,
    }
}

#[derive(Debug, Clone)]
struct LoadedSequence {
    seq: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
struct SeedHit {
    ref_id: usize,
    ref_pos: usize,
}

#[derive(Debug, Clone, Copy)]
struct Anchor {
    contig_pos: usize,
    ref_pos: usize,
}

#[derive(Debug, Clone)]
struct Placement {
    ref_id: usize,
    anchors: Vec<Anchor>,
    delta: isize,
    reversed: bool,
    secondary_cluster_support: usize,
    gap_breaks: usize,
}

#[derive(Default)]
struct ContigQualityTally {
    covered_intervals: Vec<(usize, usize)>,
    mismatch_bases: usize,
    indel_bases: usize,
    aligned_assembly_bases: usize,
    misassembly_count: usize,
}

fn evaluate_reference_quality(
    reference_path: &Path,
    assembly_path: &Path,
    assembly_source: &str,
) -> Result<CaseQualityMetrics> {
    let reference_records = load_fasta_sequences(reference_path)?;
    let assembly_records = load_fasta_sequences(assembly_path)?;
    let assembly_stats = crate::frontend::report::AssemblyStats::from_fasta(assembly_path)?;
    let seed_index = build_reference_seed_index(&reference_records);

    let mut covered_by_reference = vec![Vec::<(usize, usize)>::new(); reference_records.len()];
    let mut mismatch_bases = 0usize;
    let mut indel_bases = 0usize;
    let mut aligned_assembly_bases = 0usize;
    let mut misassembly_count = 0usize;

    for contig in &assembly_records {
        let Some(placement) = select_best_placement(&contig.seq, &seed_index, &reference_records)
        else {
            indel_bases += contig.seq.len();
            aligned_assembly_bases += contig.seq.len();
            misassembly_count += 1;
            continue;
        };

        let oriented = if placement.reversed {
            reverse_complement(&contig.seq)
        } else {
            contig.seq.clone()
        };
        let tally = evaluate_contig_against_reference(
            &oriented,
            &reference_records[placement.ref_id].seq,
            &placement,
        );
        covered_by_reference[placement.ref_id].extend(tally.covered_intervals);
        mismatch_bases += tally.mismatch_bases;
        indel_bases += tally.indel_bases;
        aligned_assembly_bases += tally.aligned_assembly_bases;
        misassembly_count += tally.misassembly_count;
    }

    let reference_bases = reference_records
        .iter()
        .map(|record| record.seq.len())
        .sum::<usize>();
    let covered_reference_bases = covered_by_reference
        .iter_mut()
        .map(|intervals| merge_interval_bases(intervals))
        .sum::<usize>();
    let mismatch_rate = if aligned_assembly_bases == 0 {
        0.0
    } else {
        mismatch_bases as f64 / aligned_assembly_bases as f64
    };
    let indel_rate = if aligned_assembly_bases == 0 {
        0.0
    } else {
        indel_bases as f64 / aligned_assembly_bases as f64
    };

    Ok(CaseQualityMetrics {
        assembly_path: assembly_path.display().to_string(),
        reference_path: reference_path.display().to_string(),
        assembly_source: assembly_source.to_string(),
        reference_bases,
        covered_reference_bases,
        aligned_assembly_bases,
        genome_fraction: if reference_bases == 0 {
            0.0
        } else {
            covered_reference_bases as f64 / reference_bases as f64
        },
        contig_count: assembly_stats.num_contigs,
        total_length: assembly_stats.total_length,
        longest_contig: assembly_stats.longest_contig,
        n50: assembly_stats.n50,
        mismatch_bases,
        indel_bases,
        mismatch_rate,
        indel_rate,
        misassembly_count,
        evaluation_method: "seed_anchored_reference_comparison".to_string(),
    })
}

fn load_fasta_sequences(path: &Path) -> Result<Vec<LoadedSequence>> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open FASTA {}", path.display()))?;
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let mut seq = record.seq().as_ref().to_vec();
        seq.make_ascii_uppercase();
        records.push(LoadedSequence { seq });
    }
    Ok(records)
}

fn build_reference_seed_index(reference_records: &[LoadedSequence]) -> HashMap<u64, Vec<SeedHit>> {
    let mut index = HashMap::<u64, Vec<SeedHit>>::new();
    for (ref_id, record) in reference_records.iter().enumerate() {
        if record.seq.len() < QUALITY_SEED_K {
            continue;
        }
        let mut ref_pos = 0usize;
        while ref_pos + QUALITY_SEED_K <= record.seq.len() {
            if let Some(kmer) = encode_kmer(&record.seq[ref_pos..ref_pos + QUALITY_SEED_K]) {
                index
                    .entry(kmer)
                    .or_default()
                    .push(SeedHit { ref_id, ref_pos });
            }
            ref_pos = ref_pos.saturating_add(QUALITY_REFERENCE_SEED_STEP);
        }
    }
    index
}

fn select_best_placement(
    contig: &[u8],
    seed_index: &HashMap<u64, Vec<SeedHit>>,
    reference_records: &[LoadedSequence],
) -> Option<Placement> {
    let forward = map_contig_orientation(contig, false, seed_index, reference_records);
    let reverse = map_contig_orientation(
        &reverse_complement(contig),
        true,
        seed_index,
        reference_records,
    );

    match (forward, reverse) {
        (Some(left), Some(right)) => {
            if left.anchors.len() >= right.anchors.len() {
                Some(left)
            } else {
                Some(right)
            }
        }
        (Some(left), None) => Some(left),
        (None, Some(right)) => Some(right),
        (None, None) => None,
    }
}

fn map_contig_orientation(
    contig: &[u8],
    reversed: bool,
    seed_index: &HashMap<u64, Vec<SeedHit>>,
    _reference_records: &[LoadedSequence],
) -> Option<Placement> {
    if contig.len() < QUALITY_SEED_K {
        return None;
    }

    let mut clusters = HashMap::<(usize, isize), Vec<Anchor>>::new();
    let mut contig_pos = 0usize;
    while contig_pos + QUALITY_SEED_K <= contig.len() {
        if let Some(kmer) = encode_kmer(&contig[contig_pos..contig_pos + QUALITY_SEED_K])
            && let Some(hits) = seed_index.get(&kmer)
            && hits.len() <= QUALITY_MAX_KMER_HITS
        {
            for hit in hits {
                let delta = hit.ref_pos as isize - contig_pos as isize;
                let bucket = delta.div_euclid(QUALITY_CLUSTER_BUCKET_BP);
                clusters
                    .entry((hit.ref_id, bucket))
                    .or_default()
                    .push(Anchor {
                        contig_pos,
                        ref_pos: hit.ref_pos,
                    });
            }
        }
        contig_pos = contig_pos.saturating_add(QUALITY_SEED_STEP);
    }

    if clusters.is_empty() {
        return None;
    }

    let mut ranked_clusters = clusters.into_iter().collect::<Vec<_>>();
    ranked_clusters.sort_by(|left, right| right.1.len().cmp(&left.1.len()));
    let ((ref_id, _bucket), anchors) = ranked_clusters.remove(0);
    let secondary_cluster_support = ranked_clusters
        .first()
        .map(|cluster| cluster.1.len())
        .unwrap_or(0);

    let anchors = filter_monotonic_anchors(anchors);
    if anchors.is_empty() {
        return None;
    }

    let delta = median_isize(
        &anchors
            .iter()
            .map(|anchor| anchor.ref_pos as isize - anchor.contig_pos as isize)
            .collect::<Vec<_>>(),
    );
    let gap_breaks = count_anchor_gap_breaks(&anchors);
    Some(Placement {
        ref_id,
        anchors,
        delta,
        reversed,
        secondary_cluster_support,
        gap_breaks,
    })
}

fn filter_monotonic_anchors(mut anchors: Vec<Anchor>) -> Vec<Anchor> {
    anchors.sort_by(|left, right| {
        left.contig_pos
            .cmp(&right.contig_pos)
            .then(left.ref_pos.cmp(&right.ref_pos))
    });
    let mut filtered = Vec::new();
    let mut last_ref_end = 0usize;
    let mut last_contig_end = 0usize;
    for anchor in anchors {
        if filtered.is_empty()
            || (anchor.contig_pos >= last_contig_end && anchor.ref_pos >= last_ref_end)
        {
            last_contig_end = anchor.contig_pos + QUALITY_SEED_K;
            last_ref_end = anchor.ref_pos + QUALITY_SEED_K;
            filtered.push(anchor);
        }
    }
    filtered
}

fn count_anchor_gap_breaks(anchors: &[Anchor]) -> usize {
    anchors
        .windows(2)
        .filter(|pair| {
            let contig_gap = pair[1]
                .contig_pos
                .saturating_sub(pair[0].contig_pos + QUALITY_SEED_K);
            let ref_gap = pair[1]
                .ref_pos
                .saturating_sub(pair[0].ref_pos + QUALITY_SEED_K);
            contig_gap.abs_diff(ref_gap) > QUALITY_MISASSEMBLY_GAP_BP
        })
        .count()
}

fn evaluate_contig_against_reference(
    contig: &[u8],
    reference: &[u8],
    placement: &Placement,
) -> ContigQualityTally {
    let mut tally = ContigQualityTally {
        misassembly_count: placement.gap_breaks,
        ..Default::default()
    };
    if placement.secondary_cluster_support > 0
        && (placement.secondary_cluster_support as f64)
            >= (placement.anchors.len() as f64 * QUALITY_SECONDARY_CLUSTER_RATIO)
    {
        tally.misassembly_count += 1;
    }

    let Some(first_anchor) = placement.anchors.first().copied() else {
        tally.indel_bases += contig.len();
        tally.aligned_assembly_bases += contig.len();
        tally.misassembly_count += 1;
        return tally;
    };
    let Some(last_anchor) = placement.anchors.last().copied() else {
        tally.indel_bases += contig.len();
        tally.aligned_assembly_bases += contig.len();
        tally.misassembly_count += 1;
        return tally;
    };

    let predicted_ref_start = placement.delta.max(0) as usize;
    compare_segment(
        &mut tally,
        &contig[..first_anchor.contig_pos],
        segment(reference, predicted_ref_start, first_anchor.ref_pos),
        predicted_ref_start,
    );

    for anchor in &placement.anchors {
        let start = anchor.ref_pos.min(reference.len());
        let end = (anchor.ref_pos + QUALITY_SEED_K).min(reference.len());
        if start < end {
            tally.covered_intervals.push((start, end));
        }
        tally.aligned_assembly_bases +=
            QUALITY_SEED_K.min(contig.len().saturating_sub(anchor.contig_pos));
    }

    for pair in placement.anchors.windows(2) {
        let left = pair[0];
        let right = pair[1];
        let contig_start = left.contig_pos + QUALITY_SEED_K;
        let contig_end = right.contig_pos.min(contig.len());
        let ref_start = (left.ref_pos + QUALITY_SEED_K).min(reference.len());
        let ref_end = right.ref_pos.min(reference.len());
        compare_segment(
            &mut tally,
            segment(contig, contig_start, contig_end),
            segment(reference, ref_start, ref_end),
            ref_start,
        );
    }

    let predicted_ref_end = if placement.delta >= 0 {
        (placement.delta as usize)
            .saturating_add(contig.len())
            .min(reference.len())
    } else {
        contig
            .len()
            .saturating_sub((-placement.delta) as usize)
            .min(reference.len())
    };
    compare_segment(
        &mut tally,
        segment(
            contig,
            last_anchor.contig_pos + QUALITY_SEED_K,
            contig.len(),
        ),
        segment(
            reference,
            last_anchor.ref_pos + QUALITY_SEED_K,
            predicted_ref_end,
        ),
        last_anchor.ref_pos + QUALITY_SEED_K,
    );

    tally
}

fn compare_segment(
    tally: &mut ContigQualityTally,
    query: &[u8],
    target: &[u8],
    target_start: usize,
) {
    let compared = query.len().min(target.len());
    tally.mismatch_bases += query
        .iter()
        .zip(target.iter())
        .take(compared)
        .filter(|(left, right)| left != right)
        .count();
    tally.indel_bases += query.len().abs_diff(target.len());
    tally.aligned_assembly_bases += query.len();
    if !target.is_empty() {
        tally
            .covered_intervals
            .push((target_start, target_start + target.len()));
    }
}

fn segment(seq: &[u8], start: usize, end: usize) -> &[u8] {
    if start >= end || start >= seq.len() {
        return &[];
    }
    &seq[start..end.min(seq.len())]
}

fn merge_interval_bases(intervals: &mut [(usize, usize)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    let mut total = 0usize;
    let mut current = intervals[0];
    for interval in intervals.iter().copied().skip(1) {
        if interval.0 <= current.1 {
            current.1 = current.1.max(interval.1);
        } else {
            total += current.1.saturating_sub(current.0);
            current = interval;
        }
    }
    total + current.1.saturating_sub(current.0)
}

fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let mut value = 0_u64;
    for &base in seq {
        let bits = match base {
            b'A' => 0_u64,
            b'C' => 1_u64,
            b'G' => 2_u64,
            b'T' => 3_u64,
            _ => return None,
        };
        value = (value << 2) | bits;
    }
    Some(value)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => *other,
        })
        .collect()
}

fn median_isize(values: &[isize]) -> isize {
    if values.is_empty() {
        return 0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2]
}

fn recover_lineabact_elapsed_seconds(outdir: &Path) -> Option<(f64, &'static str)> {
    #[derive(Deserialize)]
    struct SummaryElapsed {
        elapsed_seconds: f64,
    }

    let summary_path = outdir.join(contract::ROOT_RUN_SUMMARY);
    let text = fs::read_to_string(summary_path).ok()?;
    let summary: SummaryElapsed = serde_json::from_str(&text).ok()?;
    (summary.elapsed_seconds > 0.0).then_some((summary.elapsed_seconds, "lineabact_run_summary"))
}

fn recover_shovill_elapsed_seconds(stderr_log: &Path) -> Option<(f64, &'static str)> {
    let text = fs::read_to_string(stderr_log).ok()?;
    let line = text
        .lines()
        .find(|line| line.contains("Walltime used:"))?
        .trim();
    let walltime = line.split("Walltime used:").nth(1)?.trim();
    parse_human_elapsed_seconds(walltime).map(|value| (value, "shovill_walltime_log"))
}

fn parse_human_elapsed_seconds(text: &str) -> Option<f64> {
    let mut total_seconds = 0u64;
    let mut tokens = text.split_whitespace();
    while let Some(number_text) = tokens.next() {
        let number = number_text.parse::<u64>().ok()?;
        let unit = tokens.next()?.to_ascii_lowercase();
        let delta = if unit.starts_with("hour") || unit.starts_with("hr") {
            number * 3600
        } else if unit.starts_with("min") {
            number * 60
        } else if unit.starts_with("sec") {
            number
        } else {
            return None;
        };
        total_seconds += delta;
    }
    (total_seconds > 0).then_some(total_seconds as f64)
}

fn collect_ok_elapsed_seconds(
    cases: &[BenchmarkCaseReport],
    sample_id: &str,
    tool: &str,
) -> Vec<f64> {
    let mut values = cases
        .iter()
        .filter(|case| case.sample_id == sample_id && case.tool == tool && case.error.is_none())
        .map(|case| case.elapsed_seconds)
        .filter(|elapsed| *elapsed > 0.0)
        .collect::<Vec<_>>();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    values
}

fn collect_ok_quality_metrics<'a>(
    cases: &'a [BenchmarkCaseReport],
    sample_id: &str,
    tool: &str,
) -> Vec<&'a CaseQualityMetrics> {
    cases
        .iter()
        .filter(|case| case.sample_id == sample_id && case.tool == tool && case.error.is_none())
        .filter_map(|case| case.quality_metrics.as_ref())
        .collect()
}

fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}

fn median_usize_field(
    metrics: &[&CaseQualityMetrics],
    selector: impl Fn(&CaseQualityMetrics) -> usize,
) -> f64 {
    let mut values = metrics
        .iter()
        .map(|metrics| selector(metrics) as f64)
        .collect::<Vec<_>>();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    median(&values).unwrap_or(0.0)
}

fn median_f64_field(
    metrics: &[&CaseQualityMetrics],
    selector: impl Fn(&CaseQualityMetrics) -> f64,
) -> f64 {
    let mut values = metrics
        .iter()
        .map(|metrics| selector(metrics))
        .collect::<Vec<_>>();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    median(&values).unwrap_or(0.0)
}

fn prefixed_expected_files(prefix: &str, files: &[&str]) -> Vec<String> {
    files
        .iter()
        .map(|file| Path::new(prefix).join(file).display().to_string())
        .collect()
}

fn lineabact_expected_files() -> Vec<String> {
    let mut files = prefixed_expected_files(contract::READS_DIR, contract::READS_OUTPUT_FILES);
    files.push(contract::SPADES_COMMAND_SH.to_string());
    files.push(contract::SPADES_PLAN_JSON.to_string());
    files.extend(prefixed_expected_files(
        contract::SPADES_DIR,
        &[
            contract::SPADES_LOG,
            contract::SPADES_STDOUT_LOG,
            contract::SPADES_STDERR_LOG,
            contract::SPADES_DIAGNOSTICS_JSON,
            contract::SPADES_VERSION_TXT,
            contract::SPADES_PARAMS_TXT,
            contract::SPADES_CONTIGS_FASTA,
            contract::SPADES_SCAFFOLDS_FASTA,
            contract::SPADES_GFA,
            contract::SPADES_FASTG,
            contract::SPADES_CONTIGS_PATHS,
            contract::SPADES_SCAFFOLDS_PATHS,
        ],
    ));
    files.extend(prefixed_expected_files(
        contract::GRAPH_DIR,
        contract::GRAPH_OUTPUT_FILES,
    ));
    files.extend(prefixed_expected_files(
        contract::POSTPROCESS_DIR,
        contract::POSTPROCESS_OUTPUT_FILES,
    ));
    files.extend(prefixed_expected_files(
        contract::BRIDGING_DIR,
        contract::BRIDGING_OUTPUT_FILES,
    ));
    files.extend(
        contract::ROOT_OUTPUT_FILES
            .iter()
            .map(|value| value.to_string()),
    );
    files
}

fn require_file(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        Ok(path.to_path_buf())
    } else {
        bail!(
            "expected benchmark fixture file is missing: {}",
            path.display()
        )
    }
}

fn shell_quote(text: &str) -> String {
    if text.is_empty() {
        return "''".to_string();
    }

    if !text.contains('\'')
        && !text.contains(' ')
        && !text.contains('\t')
        && !text.contains('\n')
        && !text.contains('"')
        && !text.contains('$')
        && !text.contains('`')
        && !text.contains('\\')
    {
        return text.to_string();
    }

    format!("'{}'", text.replace('\'', "'\"'\"'"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn plan_covers_the_four_short_read_tools() {
        let temp = tempdir().unwrap();
        let args = BenchmarkPlanArgs {
            manifest: None,
            sample_limit: 3,
            fixture_root: PathBuf::from("reference_tools/Unicycler-main"),
            outdir: temp.path().join("bench"),
            sample_id: "fixture".to_string(),
            threads: 4,
            k: 55,
            lineabact_executable: "lineabact".to_string(),
            spades_executable: "spades.py".to_string(),
            shovill_executable: "shovill".to_string(),
            unicycler_executable: "unicycler".to_string(),
        };

        let plan = build_plan(&args).unwrap();
        assert_eq!(plan.samples.len(), 1);
        assert_eq!(plan.cases.len(), 4);
        assert!(plan.cases.iter().any(|case| case.tool == "lineabact"));
        assert!(plan.cases.iter().any(|case| case.tool == "spades"));
        assert!(plan.cases.iter().any(|case| case.tool == "shovill"));
        assert!(plan.cases.iter().any(|case| case.tool == "unicycler"));
        assert!(plan.cases[0].command.contains("assemble"));
        assert!(plan.cases[0].command.contains("--trim-adapters"));
        assert!(plan.cases[0].command.contains("--target-coverage 150"));
        assert!(plan.cases[0].command.contains("--genome-size-bp"));
        assert!(plan.cases[1].command.contains("--isolate"));
        assert!(plan.cases[2].command.contains("--assembler spades"));
        assert!(plan.cases[3].command.contains("--spades_path"));
    }

    #[test]
    fn benchmark_run_defaults_to_three_repeats() {
        use clap::Parser;

        let cli = crate::cli::Cli::parse_from([
            "lineabact",
            "stats",
            "benchmark-run",
            "--plan",
            "benchmarks/generated/reference_validation/benchmark_plan.json",
        ]);

        let crate::cli::Commands::Stats(stats_args) = cli.command else {
            panic!("expected stats command");
        };
        let crate::cli::StatsCommands::Run(run_args) = stats_args.command else {
            panic!("expected benchmark-run subcommand");
        };
        assert_eq!(run_args.repeat_count, 3);
    }

    #[test]
    fn benchmark_run_executes_cases_and_records_results() {
        let temp = tempdir().unwrap();
        let success_outdir = temp.path().join("success");
        let failure_outdir = temp.path().join("failure");
        let cases = vec![
            BenchmarkCase {
                sample_id: "synthetic_a".to_string(),
                tool: "ok-tool".to_string(),
                executable: "bash".to_string(),
                command: format!(
                    "mkdir -p {} {}/subdir && touch {}/result.txt",
                    shell_quote(&success_outdir.display().to_string()),
                    shell_quote(&success_outdir.display().to_string()),
                    shell_quote(&success_outdir.display().to_string()),
                ),
                outdir: success_outdir.display().to_string(),
                expected_dirs: vec!["subdir".to_string()],
                expected_files: vec!["result.txt".to_string()],
                notes: "success".to_string(),
            },
            BenchmarkCase {
                sample_id: "synthetic_b".to_string(),
                tool: "bad-tool".to_string(),
                executable: "bash".to_string(),
                command: "exit 7".to_string(),
                outdir: failure_outdir.display().to_string(),
                expected_dirs: Vec::new(),
                expected_files: vec!["missing.txt".to_string()],
                notes: "failure".to_string(),
            },
        ];

        let run_args = BenchmarkRunArgs {
            plan: temp.path().join("plan.json"),
            outdir: temp.path().join("report"),
            stop_on_error: false,
            repeat_count: 1,
            scheduler: BenchmarkScheduler::Local,
            slurm_partition: "qcpu_23if".to_string(),
            slurm_conda_base: "/hpcfs/fpublic/app/miniforge3/conda".to_string(),
            slurm_conda_env: "LineaBact".to_string(),
            slurm_cpus_per_task: 4,
            slurm_mem_gb: None,
            slurm_time: "12:00:00".to_string(),
            slurm_dry_run: false,
        };

        let reports = run_cases(&cases, &run_args).unwrap();
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].status, "ok");
        assert_eq!(reports[0].sample_id, "synthetic_a");
        assert_eq!(reports[0].repeat_index, 1);
        assert_eq!(reports[0].repeat_count, 1);
        assert_eq!(reports[0].scheduler, "local");
        assert_eq!(reports[0].exit_code, Some(0));
        assert!(Path::new(&reports[0].stdout_log).exists());
        assert!(Path::new(&reports[0].stderr_log).exists());

        assert_eq!(reports[1].status, "command_failed");
        assert_eq!(reports[1].sample_id, "synthetic_b");
        assert_eq!(reports[1].exit_code, Some(7));
        assert_eq!(reports[1].missing_files, vec!["missing.txt".to_string()]);
    }

    #[test]
    fn runtime_comparison_reports_lineabact_vs_shovill_ratio() {
        let cases = vec![
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 1,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "lineabact".to_string(),
                executable: "lineabact".to_string(),
                command: "lineabact assemble".to_string(),
                outdir: "out/lineabact".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 8.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 2,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "lineabact".to_string(),
                executable: "lineabact".to_string(),
                command: "lineabact assemble".to_string(),
                outdir: "out/lineabact".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 6.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 3,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "lineabact".to_string(),
                executable: "lineabact".to_string(),
                command: "lineabact assemble".to_string(),
                outdir: "out/lineabact".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 8.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 1,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "shovill".to_string(),
                executable: "shovill".to_string(),
                command: "shovill".to_string(),
                outdir: "out/shovill".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 10.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 2,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "shovill".to_string(),
                executable: "shovill".to_string(),
                command: "shovill".to_string(),
                outdir: "out/shovill".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 12.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
            BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 3,
                repeat_count: 3,
                scheduler: "local".to_string(),
                tool: "shovill".to_string(),
                executable: "shovill".to_string(),
                command: "shovill".to_string(),
                outdir: "out/shovill".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 10.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            },
        ];

        let comparisons = build_runtime_comparisons(&cases);
        assert_eq!(comparisons.len(), 1);
        assert_eq!(comparisons[0].sample_id, "sample_a");
        assert_eq!(comparisons[0].baseline_tool, "shovill");
        assert_eq!(comparisons[0].candidate_tool, "lineabact");
        assert_eq!(comparisons[0].baseline_repeat_count, 3);
        assert_eq!(comparisons[0].candidate_repeat_count, 3);
        assert!((comparisons[0].baseline_median_elapsed_seconds - 10.0).abs() < 1e-9);
        assert!((comparisons[0].candidate_median_elapsed_seconds - 8.0).abs() < 1e-9);
        assert!((comparisons[0].median_runtime_ratio_vs_baseline - 0.8).abs() < 1e-9);
        assert!((comparisons[0].median_speedup_vs_baseline - 1.25).abs() < 1e-9);
        assert!(comparisons[0].meets_target);
    }

    #[test]
    fn benchmark_merge_writes_runtime_comparisons_from_two_reports() {
        let temp = tempdir().unwrap();
        let report_dir = temp.path().join("report_lineabact");
        fs::create_dir_all(&report_dir).unwrap();

        let primary_report = BenchmarkReport {
            schema_version: SCHEMA_VERSION,
            plan: "primary-plan.json".to_string(),
            case_count: 1,
            samples: vec![BenchmarkSample {
                sample_id: "sample_a".to_string(),
                r1: "r1".to_string(),
                r2: "r2".to_string(),
                reference: "ref.fa".to_string(),
            }],
            cases: vec![BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 1,
                repeat_count: 3,
                scheduler: "slurm".to_string(),
                tool: "lineabact".to_string(),
                executable: "lineabact".to_string(),
                command: "lineabact assemble".to_string(),
                outdir: "out/lineabact".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 8.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            }],
            runtime_comparisons: Vec::new(),
            quality_comparisons: Vec::new(),
            sample_acceptance: Vec::new(),
        };
        let compare_report = BenchmarkReport {
            schema_version: SCHEMA_VERSION,
            plan: "compare-plan.json".to_string(),
            case_count: 1,
            samples: vec![BenchmarkSample {
                sample_id: "sample_a".to_string(),
                r1: "r1".to_string(),
                r2: "r2".to_string(),
                reference: "ref.fa".to_string(),
            }],
            cases: vec![BenchmarkCaseReport {
                sample_id: "sample_a".to_string(),
                repeat_index: 1,
                repeat_count: 3,
                scheduler: "slurm".to_string(),
                tool: "shovill".to_string(),
                executable: "shovill".to_string(),
                command: "shovill".to_string(),
                outdir: "out/shovill".to_string(),
                expected_dirs: Vec::new(),
                expected_files: Vec::new(),
                status: "submitted".to_string(),
                exit_code: Some(0),
                elapsed_seconds: 10.0,
                elapsed_seconds_source: "process_wallclock".to_string(),
                missing_dirs: Vec::new(),
                missing_files: Vec::new(),
                stdout_log: "stdout".to_string(),
                stderr_log: "stderr".to_string(),
                submission_script: None,
                submission_command: None,
                job_id: None,
                quality_metrics: None,
                error: None,
            }],
            runtime_comparisons: Vec::new(),
            quality_comparisons: Vec::new(),
            sample_acceptance: Vec::new(),
        };

        let report_path = report_dir.join("benchmark_report.json");
        let compare_path = temp.path().join("benchmark_compare.json");
        io::write_json(&report_path, &primary_report).unwrap();
        io::write_json(&compare_path, &compare_report).unwrap();

        benchmark_merge(BenchmarkMergeArgs {
            report: report_path.clone(),
            compare_report: compare_path.clone(),
        })
        .unwrap();

        let merged_text = fs::read_to_string(&report_path).unwrap();
        let merged_report: BenchmarkReport = serde_json::from_str(&merged_text).unwrap();
        assert_eq!(merged_report.runtime_comparisons.len(), 1);
        assert_eq!(merged_report.runtime_comparisons[0].sample_id, "sample_a");
        assert!(
            (merged_report.runtime_comparisons[0].median_runtime_ratio_vs_baseline - 0.8).abs()
                < 1e-9
        );
        assert!(report_dir.join("benchmark_runtime_comparison.tsv").exists());
    }

    #[test]
    fn slurm_dry_run_writes_submission_script() {
        let temp = tempdir().unwrap();
        let env_prefix = temp.path().join("LineaBact");
        fs::create_dir_all(env_prefix.join("bin")).unwrap();
        let case = BenchmarkCase {
            sample_id: "synthetic_slurm".to_string(),
            tool: "lineabact".to_string(),
            executable: "lineabact".to_string(),
            command: "lineabact assemble --r1 a --r2 b --outdir x".to_string(),
            outdir: temp.path().join("slurm_case").display().to_string(),
            expected_dirs: vec!["reads".to_string()],
            expected_files: vec!["run_summary.json".to_string()],
            notes: "slurm dry run".to_string(),
        };
        let run_args = BenchmarkRunArgs {
            plan: temp.path().join("plan.json"),
            outdir: temp.path().join("report"),
            stop_on_error: false,
            repeat_count: 1,
            scheduler: BenchmarkScheduler::Slurm,
            slurm_partition: "qcpu_23if".to_string(),
            slurm_conda_base: env_prefix.display().to_string(),
            slurm_conda_env: "LineaBact".to_string(),
            slurm_cpus_per_task: 8,
            slurm_mem_gb: Some(64),
            slurm_time: "08:00:00".to_string(),
            slurm_dry_run: true,
        };

        let report = submit_case_to_slurm(&case, &run_args, 1, 1).unwrap();
        assert_eq!(report.scheduler, "slurm");
        assert_eq!(report.status, "scripted");
        assert!(report.job_id.is_none());
        let script = report.submission_script.unwrap();
        assert!(Path::new(&script).exists());
        let command = report.submission_command.unwrap();
        assert!(command.contains("--partition"));
        assert!(command.contains("qcpu_23if"));
        assert!(command.contains("--cpus-per-task"));
        assert!(command.contains("8"));
        assert!(command.contains("--mem"));
        assert!(command.contains("64G"));
        let script_text = fs::read_to_string(script).unwrap();
        assert!(script_text.contains("export CONDA_PREFIX="));
        assert!(script_text.contains("export PATH="));
    }

    #[test]
    fn repeated_case_rewrites_outdir_once() {
        let case = BenchmarkCase {
            sample_id: "sample".to_string(),
            tool: "lineabact".to_string(),
            executable: "lineabact".to_string(),
            command: "lineabact assemble --outdir '/tmp/outdir' --r1 a --r2 b".to_string(),
            outdir: "/tmp/outdir".to_string(),
            expected_dirs: Vec::new(),
            expected_files: Vec::new(),
            notes: String::new(),
        };

        let repeated = repeated_case(&case, 2);
        assert_eq!(repeated.outdir, "/tmp/outdir/repeats/run_002");
        assert!(
            repeated
                .command
                .contains("--outdir '/tmp/outdir/repeats/run_002'")
        );
        assert!(!repeated.command.contains("run_002/repeats/run_002"));
    }

    #[test]
    fn parse_partition_candidates_supports_auto_and_lists() {
        let auto = parse_partition_candidates("auto").unwrap();
        assert_eq!(
            auto,
            vec![
                "qcpu_23if".to_string(),
                "qcpu_23i".to_string(),
                "qcpu_23a".to_string(),
                "qcpu_18i".to_string()
            ]
        );

        let explicit = parse_partition_candidates("qcpu_23i, qcpu_23a qcpu_18i").unwrap();
        assert_eq!(
            explicit,
            vec![
                "qcpu_23i".to_string(),
                "qcpu_23a".to_string(),
                "qcpu_18i".to_string()
            ]
        );
    }

    #[test]
    fn manifest_plan_filters_to_valid_paired_samples() {
        let temp = tempdir().unwrap();
        let sample_root = temp.path().join("samples");
        fs::create_dir_all(&sample_root).unwrap();

        let valid_r1 = sample_root.join("valid_1.fastq.gz");
        let valid_r2 = sample_root.join("valid_2.fastq.gz");
        let valid_ref = sample_root.join("valid.fna");
        fs::write(&valid_r1, b"x").unwrap();
        fs::write(&valid_r2, b"x").unwrap();
        fs::write(&valid_ref, b">chr1\nACGTACGT\n").unwrap();

        let manifest = temp.path().join("manifest.tsv");
        fs::write(
            &manifest,
            format!(
                concat!(
                    "species\tslug\treference_fasta\tillumina_fastq_1\tillumina_fastq_2\n",
                    "Valid sample\tvalid\t{}\t{}\t{}\n",
                    "Missing r2\tmissing_r2\t{}\t{}\t\n",
                    "Missing ref\tmissing_ref\t\t{}\t{}\n"
                ),
                valid_ref.display(),
                valid_r1.display(),
                valid_r2.display(),
                valid_ref.display(),
                valid_r1.display(),
                valid_r1.display(),
                valid_r2.display(),
            ),
        )
        .unwrap();

        let args = BenchmarkPlanArgs {
            manifest: Some(manifest),
            sample_limit: 10,
            fixture_root: PathBuf::from("reference_tools/Unicycler-main"),
            outdir: temp.path().join("bench"),
            sample_id: "fixture".to_string(),
            threads: 4,
            k: 55,
            lineabact_executable: "lineabact".to_string(),
            spades_executable: "spades.py".to_string(),
            shovill_executable: "shovill".to_string(),
            unicycler_executable: "unicycler".to_string(),
        };

        let plan = build_plan(&args).unwrap();
        assert_eq!(plan.samples.len(), 1);
        assert_eq!(plan.samples[0].sample_id, "valid");
        assert_eq!(plan.cases.len(), 4);
        assert!(plan.cases.iter().all(|case| case.sample_id == "valid"));
    }
}
