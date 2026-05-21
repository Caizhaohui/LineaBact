use anyhow::{Context, Result, bail};
use flate2::{Compression, write::GzEncoder};
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::ffi::OsStr;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::io;

use super::contract;
use super::report::ReadReport;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadArtifacts {
    pub report: ReadReport,
    pub optimized_r1: PathBuf,
    pub optimized_r2: PathBuf,
    pub summary_json: PathBuf,
    pub reads_stats_tsv: PathBuf,
    pub pair_check_tsv: PathBuf,
    pub downsample_plan_tsv: PathBuf,
    pub trim_plan_tsv: PathBuf,
}

#[derive(Debug, Clone)]
pub struct ReadPreparationConfig<'a> {
    pub downsample_fraction: Option<f64>,
    pub target_coverage: Option<f64>,
    pub genome_size_bp: Option<u64>,
    pub downsample_seed: u64,
    pub trim_adapters: bool,
    pub trim_tool: &'a str,
    pub read_sha256: bool,
}

pub fn prepare_reads(
    r1: &Path,
    r2: &Path,
    outdir: &Path,
    config: &ReadPreparationConfig<'_>,
) -> Result<ReadArtifacts> {
    fs::create_dir_all(outdir)?;
    let optimized_r1 = contract::root_path(outdir, contract::READS_OPTIMIZED_R1);
    let optimized_r2 = contract::root_path(outdir, contract::READS_OPTIMIZED_R2);
    let staged_r1 = outdir.join("staged_R1.fastq.gz");
    let staged_r2 = outdir.join("staged_R2.fastq.gz");
    let summary_json = contract::root_path(outdir, contract::READS_SUMMARY_JSON);
    let reads_stats_tsv = contract::root_path(outdir, contract::READS_STATS_TSV);
    let pair_check_tsv = contract::root_path(outdir, contract::READS_PAIR_CHECK_TSV);
    let downsample_plan_tsv = contract::root_path(outdir, contract::READS_DOWNSAMPLE_PLAN_TSV);
    let trim_plan_tsv = contract::root_path(outdir, contract::READS_TRIM_PLAN_TSV);

    let inspection = inspect_read_pair(r1, r2, config.read_sha256)?;
    let input_total_bases = inspection.input_bases_r1 + inspection.input_bases_r2;
    let estimated_depth = config
        .genome_size_bp
        .filter(|value| *value > 0)
        .map(|value| input_total_bases as f64 / value as f64);

    if let Some(coverage) = config.target_coverage {
        if coverage <= 0.0 {
            bail!("target_coverage must be > 0, got {coverage}");
        }
        if config.genome_size_bp.is_none() {
            bail!("target_coverage requires --genome-size-bp");
        }
    }

    let mut retained_pairs = inspection.input_pairs;
    let mut retained_bases_r1 = inspection.input_bases_r1;
    let mut retained_bases_r2 = inspection.input_bases_r2;
    let derived_downsample_fraction = derive_target_coverage_fraction(
        input_total_bases,
        config.genome_size_bp,
        config.target_coverage,
    );
    let applied_downsample_fraction = select_applied_downsample_fraction(
        config.downsample_fraction,
        derived_downsample_fraction,
    )?;

    let mut working_r1 = r1.to_path_buf();
    let mut working_r2 = r2.to_path_buf();
    let staged_reads_written = applied_downsample_fraction < 1.0;

    if staged_reads_written {
        let (pairs, bases_r1, bases_r2, input_pairs) = write_downsampled_pair(
            r1,
            r2,
            &staged_r1,
            &staged_r2,
            applied_downsample_fraction,
            config.downsample_seed,
        )?;
        if input_pairs != inspection.input_pairs {
            bail!("internal read sampling error changed the input pair count");
        }
        retained_pairs = pairs;
        retained_bases_r1 = bases_r1;
        retained_bases_r2 = bases_r2;
        working_r1 = staged_r1.clone();
        working_r2 = staged_r2.clone();
    }

    if config.trim_adapters {
        execute_trimming(
            config.trim_tool,
            &working_r1,
            &working_r2,
            &optimized_r1,
            &optimized_r2,
        )?;
    } else {
        io::link_or_copy_file(&working_r1, &optimized_r1)?;
        io::link_or_copy_file(&working_r2, &optimized_r2)?;
    }

    if let Some(fraction) = config.downsample_fraction
        && !(0.0 < fraction && fraction <= 1.0)
    {
        bail!("downsample_fraction must be in the interval (0, 1], got {fraction}");
    }

    if retained_pairs > inspection.input_pairs {
        bail!("internal read sampling error produced more retained pairs than input pairs");
    }

    let output_sha256_r1 = if config.read_sha256 {
        Some(sha256_hex(&optimized_r1)?)
    } else {
        None
    };
    let output_sha256_r2 = if config.read_sha256 {
        Some(sha256_hex(&optimized_r2)?)
    } else {
        None
    };
    let mean_read_length_r1 = if inspection.input_reads_r1 == 0 {
        0.0
    } else {
        inspection.input_bases_r1 as f64 / inspection.input_reads_r1 as f64
    };
    let mean_read_length_r2 = if inspection.input_reads_r2 == 0 {
        0.0
    } else {
        inspection.input_bases_r2 as f64 / inspection.input_reads_r2 as f64
    };
    let retained_total_bases = retained_bases_r1 + retained_bases_r2;

    let report = ReadReport {
        input_pairs: inspection.input_pairs,
        retained_pairs,
        dropped_pairs: inspection.input_pairs.saturating_sub(retained_pairs),
        input_reads_r1: inspection.input_reads_r1,
        input_reads_r2: inspection.input_reads_r2,
        retained_reads_r1: retained_pairs,
        retained_reads_r2: retained_pairs,
        input_bases_r1: inspection.input_bases_r1,
        input_bases_r2: inspection.input_bases_r2,
        retained_bases_r1,
        retained_bases_r2,
        mean_read_length_r1,
        mean_read_length_r2,
        input_total_bases,
        retained_total_bases,
        genome_size_bp: config.genome_size_bp,
        estimated_depth,
        target_coverage: config.target_coverage,
        derived_downsample_fraction,
        applied_downsample_fraction,
        downsample_fraction: config.downsample_fraction,
        downsample_seed: config.downsample_seed,
        trim_adapters: config.trim_adapters,
        trim_tool: config.trim_tool.to_string(),
        read_sha256: config.read_sha256,
        input_sha256_r1: inspection.sha256_r1,
        input_sha256_r2: inspection.sha256_r2,
        output_sha256_r1,
        output_sha256_r2,
        optimized_r1: optimized_r1.display().to_string(),
        optimized_r2: optimized_r2.display().to_string(),
    };

    io::write_json(&summary_json, &report)?;
    write_reads_stats_tsv(&report, &reads_stats_tsv)?;
    write_pair_check_tsv(&pair_check_tsv)?;
    write_downsample_plan_tsv(
        &downsample_plan_tsv,
        &report,
        config.downsample_fraction,
        derived_downsample_fraction,
    )?;
    write_trim_plan_tsv(
        &trim_plan_tsv,
        config.trim_adapters,
        config.trim_tool,
        &working_r1,
        &working_r2,
        &optimized_r1,
        &optimized_r2,
    )?;

    if staged_reads_written {
        let _ = fs::remove_file(&staged_r1);
        let _ = fs::remove_file(&staged_r2);
    }

    Ok(ReadArtifacts {
        report,
        optimized_r1,
        optimized_r2,
        summary_json,
        reads_stats_tsv,
        pair_check_tsv,
        downsample_plan_tsv,
        trim_plan_tsv,
    })
}

struct ReadInspection {
    input_pairs: u64,
    input_reads_r1: u64,
    input_reads_r2: u64,
    input_bases_r1: u64,
    input_bases_r2: u64,
    sha256_r1: Option<String>,
    sha256_r2: Option<String>,
}

fn inspect_read_pair(r1: &Path, r2: &Path, read_sha256: bool) -> Result<ReadInspection> {
    let input_reads_r1 = count_fastq_records(r1)?;
    let input_reads_r2 = count_fastq_records(r2)?;
    if input_reads_r1 != input_reads_r2 {
        bail!(
            "paired FASTQ files have different record counts: R1={} R2={}",
            input_reads_r1,
            input_reads_r2
        );
    }

    let input_bases_r1 = count_fastq_bases(r1)?;
    let input_bases_r2 = count_fastq_bases(r2)?;
    validate_pair_names(r1, r2)?;

    Ok(ReadInspection {
        input_pairs: input_reads_r1,
        input_reads_r1,
        input_reads_r2,
        input_bases_r1,
        input_bases_r2,
        sha256_r1: if read_sha256 {
            Some(sha256_hex(r1)?)
        } else {
            None
        },
        sha256_r2: if read_sha256 {
            Some(sha256_hex(r2)?)
        } else {
            None
        },
    })
}

fn derive_target_coverage_fraction(
    input_total_bases: u64,
    genome_size_bp: Option<u64>,
    target_coverage: Option<f64>,
) -> Option<f64> {
    let genome_size_bp = genome_size_bp?;
    let target_coverage = target_coverage?;
    if genome_size_bp == 0 || input_total_bases == 0 {
        return Some(1.0);
    }
    let target_bases = target_coverage * genome_size_bp as f64;
    Some((target_bases / input_total_bases as f64).clamp(0.0, 1.0))
}

fn select_applied_downsample_fraction(
    user_fraction: Option<f64>,
    derived_fraction: Option<f64>,
) -> Result<f64> {
    if let Some(value) = user_fraction
        && !(0.0 < value && value <= 1.0)
    {
        bail!("downsample_fraction must be in the interval (0, 1], got {value}");
    }

    let selected = match (user_fraction, derived_fraction) {
        (Some(user), Some(derived)) => user.min(derived),
        (Some(user), None) => user,
        (None, Some(derived)) => derived,
        (None, None) => 1.0,
    };

    if selected <= 0.0 {
        bail!("applied downsample fraction must be > 0, got {selected}");
    }
    Ok(selected.min(1.0))
}

fn write_reads_stats_tsv(report: &ReadReport, path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "input_pairs\t{}", report.input_pairs)?;
    writeln!(writer, "retained_pairs\t{}", report.retained_pairs)?;
    writeln!(writer, "input_reads_r1\t{}", report.input_reads_r1)?;
    writeln!(writer, "input_reads_r2\t{}", report.input_reads_r2)?;
    writeln!(writer, "input_bases_r1\t{}", report.input_bases_r1)?;
    writeln!(writer, "input_bases_r2\t{}", report.input_bases_r2)?;
    writeln!(writer, "input_total_bases\t{}", report.input_total_bases)?;
    writeln!(
        writer,
        "retained_total_bases\t{}",
        report.retained_total_bases
    )?;
    writeln!(
        writer,
        "mean_read_length_r1\t{:.2}",
        report.mean_read_length_r1
    )?;
    writeln!(
        writer,
        "mean_read_length_r2\t{:.2}",
        report.mean_read_length_r2
    )?;
    writeln!(
        writer,
        "genome_size_bp\t{}",
        report
            .genome_size_bp
            .map(|value| value.to_string())
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "estimated_depth\t{}",
        report
            .estimated_depth
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "target_coverage\t{}",
        report
            .target_coverage
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "applied_downsample_fraction\t{:.6}",
        report.applied_downsample_fraction
    )?;
    writer.flush()?;
    Ok(())
}

fn write_pair_check_tsv(path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "check\tstatus\tdetail")?;
    writeln!(
        writer,
        "read_count_match\tpass\tR1 and R2 record counts are equal"
    )?;
    writeln!(
        writer,
        "pair_name_match\tpass\tR1 and R2 read names match after normalization"
    )?;
    writer.flush()?;
    Ok(())
}

fn write_downsample_plan_tsv(
    path: &Path,
    report: &ReadReport,
    user_fraction: Option<f64>,
    derived_fraction: Option<f64>,
) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "key\tvalue")?;
    writeln!(
        writer,
        "user_downsample_fraction\t{}",
        user_fraction
            .map(|value| format!("{value:.6}"))
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "derived_downsample_fraction\t{}",
        derived_fraction
            .map(|value| format!("{value:.6}"))
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "applied_downsample_fraction\t{:.6}",
        report.applied_downsample_fraction
    )?;
    writeln!(
        writer,
        "estimated_depth\t{}",
        report
            .estimated_depth
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default()
    )?;
    writeln!(
        writer,
        "estimated_retained_depth\t{}",
        match (report.estimated_depth, report.applied_downsample_fraction) {
            (Some(depth), fraction) => format!("{:.4}", depth * fraction),
            _ => String::new(),
        }
    )?;
    writeln!(writer, "retained_pairs\t{}", report.retained_pairs)?;
    writer.flush()?;
    Ok(())
}

fn write_trim_plan_tsv(
    path: &Path,
    trim_adapters: bool,
    trim_tool: &str,
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "key\tvalue")?;
    writeln!(writer, "trim_adapters\t{}", trim_adapters)?;
    writeln!(writer, "trim_tool\t{}", trim_tool)?;
    writeln!(writer, "input_r1\t{}", input_r1.display())?;
    writeln!(writer, "input_r2\t{}", input_r2.display())?;
    writeln!(writer, "planned_output_r1\t{}", output_r1.display())?;
    writeln!(writer, "planned_output_r2\t{}", output_r2.display())?;
    if trim_adapters {
        let tool_kind = trim_tool_kind(trim_tool);
        let command = match tool_kind {
            TrimToolKind::Fastp => format!(
                "{tool} -i {r1} -I {r2} -o {o1} -O {o2} --thread 1",
                tool = trim_tool,
                r1 = input_r1.display(),
                r2 = input_r2.display(),
                o1 = output_r1.display(),
                o2 = output_r2.display()
            ),
            TrimToolKind::Seqtk => format!(
                "{tool} trimfq {r1} | pigz > {o1} && {tool} trimfq {r2} | pigz > {o2}",
                tool = trim_tool,
                r1 = input_r1.display(),
                r2 = input_r2.display(),
                o1 = output_r1.display(),
                o2 = output_r2.display()
            ),
        };
        writeln!(writer, "planned_command\t{command}")?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum TrimToolKind {
    Fastp,
    Seqtk,
}

fn trim_tool_kind(trim_tool: &str) -> TrimToolKind {
    let file_name = Path::new(trim_tool)
        .file_name()
        .and_then(OsStr::to_str)
        .unwrap_or(trim_tool)
        .to_ascii_lowercase();
    if file_name.contains("seqtk") {
        TrimToolKind::Seqtk
    } else {
        TrimToolKind::Fastp
    }
}

fn execute_trimming(
    trim_tool: &str,
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
) -> Result<()> {
    match trim_tool_kind(trim_tool) {
        TrimToolKind::Fastp => {
            if Command::new(trim_tool).arg("--version").output().is_ok() {
                run_fastp(trim_tool, input_r1, input_r2, output_r1, output_r2)
            } else {
                let fallback_seqtk = "/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/seqtk";
                if Path::new(fallback_seqtk).exists() {
                    run_seqtk_trimfq(fallback_seqtk, input_r1, input_r2, output_r1, output_r2)
                } else {
                    bail!(
                        "trim_adapters requested but neither {} nor seqtk fallback is available",
                        trim_tool
                    );
                }
            }
        }
        TrimToolKind::Seqtk => {
            run_seqtk_trimfq(trim_tool, input_r1, input_r2, output_r1, output_r2)
        }
    }
}

fn run_fastp(
    executable: &str,
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
) -> Result<()> {
    let status = Command::new(executable)
        .arg("-i")
        .arg(input_r1)
        .arg("-I")
        .arg(input_r2)
        .arg("-o")
        .arg(output_r1)
        .arg("-O")
        .arg(output_r2)
        .arg("--thread")
        .arg("1")
        .arg("--disable_length_filtering")
        .arg("--dont_overwrite")
        .status()
        .with_context(|| format!("failed to execute trimming command {}", executable))?;
    if !status.success() {
        bail!(
            "trimming command {} exited with status {}",
            executable,
            status
        );
    }
    Ok(())
}

fn run_seqtk_trimfq(
    executable: &str,
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
) -> Result<()> {
    trim_one_with_seqtk(executable, input_r1, output_r1)?;
    trim_one_with_seqtk(executable, input_r2, output_r2)?;
    Ok(())
}

fn trim_one_with_seqtk(executable: &str, input: &Path, output: &Path) -> Result<()> {
    io::ensure_parent_dir(output)?;
    let pigz = if Path::new("/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/pigz").exists() {
        "/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/pigz"
    } else {
        "pigz"
    };

    let mut seqtk = Command::new(executable)
        .arg("trimfq")
        .arg(input)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| format!("failed to start seqtk trimfq on {}", input.display()))?;
    let seqtk_stdout = seqtk
        .stdout
        .take()
        .context("seqtk trimfq stdout pipe is unavailable")?;
    let output_file = fs::File::create(output)?;
    let mut pigz_child = Command::new(pigz)
        .arg("-c")
        .stdin(Stdio::from(seqtk_stdout))
        .stdout(Stdio::from(output_file))
        .spawn()
        .with_context(|| format!("failed to start pigz for {}", output.display()))?;

    let seqtk_status = seqtk.wait()?;
    let pigz_status = pigz_child.wait()?;
    if !seqtk_status.success() {
        bail!("seqtk trimfq exited with status {}", seqtk_status);
    }
    if !pigz_status.success() {
        bail!("pigz exited with status {}", pigz_status);
    }
    Ok(())
}

fn count_fastq_records(path: &Path) -> Result<u64> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open FASTQ/FASTA file {}", path.display()))?;
    let mut count = 0u64;
    while let Some(record) = reader.next() {
        record?;
        count += 1;
    }
    Ok(count)
}

fn count_fastq_bases(path: &Path) -> Result<u64> {
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("failed to open FASTQ/FASTA file {}", path.display()))?;
    let mut bases = 0u64;
    while let Some(record) = reader.next() {
        let record = record?;
        bases += record.seq().len() as u64;
    }
    Ok(bases)
}

fn validate_pair_names(r1: &Path, r2: &Path) -> Result<()> {
    let mut left = parse_fastx_file(r1)?;
    let mut right = parse_fastx_file(r2)?;
    let mut index = 0u64;

    loop {
        match (left.next(), right.next()) {
            (None, None) => break,
            (None, Some(_)) | (Some(_), None) => {
                bail!("paired FASTQ files have different lengths")
            }
            (Some(left_record), Some(right_record)) => {
                let left_record = left_record?;
                let right_record = right_record?;
                let left_name = normalize_read_name(left_record.id());
                let right_name = normalize_read_name(right_record.id());
                if left_name != right_name {
                    bail!(
                        "paired read names differ at pair {}: {} vs {}",
                        index + 1,
                        left_name,
                        right_name
                    );
                }
                index += 1;
            }
        }
    }

    Ok(())
}

fn normalize_read_name(id: &[u8]) -> String {
    let id = String::from_utf8_lossy(id);
    let first = id.split_whitespace().next().unwrap_or(&id);
    if let Some(stripped) = first.strip_suffix("/1") {
        stripped.to_string()
    } else if let Some(stripped) = first.strip_suffix("/2") {
        stripped.to_string()
    } else {
        first.to_string()
    }
}

fn write_downsampled_pair(
    r1: &Path,
    r2: &Path,
    out_r1: &Path,
    out_r2: &Path,
    fraction: f64,
    seed: u64,
) -> Result<(u64, u64, u64, u64)> {
    io::ensure_parent_dir(out_r1)?;
    io::ensure_parent_dir(out_r2)?;

    let mut left = parse_fastx_file(r1)?;
    let mut right = parse_fastx_file(r2)?;

    let writer_r1 = fs::File::create(out_r1)?;
    let writer_r2 = fs::File::create(out_r2)?;
    let mut writer_r1 = GzEncoder::new(BufWriter::new(writer_r1), Compression::default());
    let mut writer_r2 = GzEncoder::new(BufWriter::new(writer_r2), Compression::default());

    let mut retained_pairs = 0u64;
    let mut retained_bases_r1 = 0u64;
    let mut retained_bases_r2 = 0u64;

    let mut index = 0u64;
    loop {
        match (left.next(), right.next()) {
            (None, None) => break,
            (None, Some(_)) | (Some(_), None) => {
                bail!("paired FASTQ files have different lengths")
            }
            (Some(left_record), Some(right_record)) => {
                let left_record = left_record?;
                let right_record = right_record?;
                if normalize_read_name(left_record.id()) != normalize_read_name(right_record.id()) {
                    bail!("paired read names differ at pair {}", index + 1);
                }
                if keep_pair(index, fraction, seed) {
                    write_fastq_record(
                        &mut writer_r1,
                        left_record.id(),
                        left_record.seq().as_ref(),
                        left_record.qual(),
                    )?;
                    write_fastq_record(
                        &mut writer_r2,
                        right_record.id(),
                        right_record.seq().as_ref(),
                        right_record.qual(),
                    )?;
                    retained_pairs += 1;
                    retained_bases_r1 += left_record.seq().len() as u64;
                    retained_bases_r2 += right_record.seq().len() as u64;
                }
                index += 1;
            }
        }
    }

    writer_r1.finish()?;
    writer_r2.finish()?;

    Ok((retained_pairs, retained_bases_r1, retained_bases_r2, index))
}

fn keep_pair(index: u64, fraction: f64, seed: u64) -> bool {
    if fraction >= 1.0 {
        return true;
    }

    let mut hasher = Sha256::new();
    hasher.update(seed.to_le_bytes());
    hasher.update(index.to_le_bytes());
    let digest = hasher.finalize();
    let mut bytes = [0u8; 8];
    bytes.copy_from_slice(&digest[..8]);
    let value = u64::from_le_bytes(bytes);
    let normalized = (value as f64) / (u64::MAX as f64);
    normalized < fraction
}

fn write_fastq_record<W: Write>(
    writer: &mut W,
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> Result<()> {
    let qual = qual.context("FASTQ record is missing quality values")?;
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn sha256_hex(path: &Path) -> Result<String> {
    let bytes = fs::read(path)?;
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    Ok(hex::encode(hasher.finalize()))
}
