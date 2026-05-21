use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use crate::cli::{AssembleArgs, AssembleBackend};
use crate::io;

use super::contract;
use super::reads::ReadArtifacts;
use super::report::BackendPlan;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendArtifacts {
    pub plan: BackendPlan,
    pub plan_script: PathBuf,
    pub plan_json: PathBuf,
    pub spades_dir: PathBuf,
    pub log_path: PathBuf,
    pub stdout_log_path: PathBuf,
    pub stderr_log_path: PathBuf,
    pub diagnostics_json: PathBuf,
    pub version_txt: PathBuf,
    pub contigs_fasta: PathBuf,
    pub scaffolds_fasta: PathBuf,
    pub assembly_graph_gfa: PathBuf,
    pub assembly_graph_fastg: PathBuf,
    pub contigs_paths: PathBuf,
    pub scaffolds_paths: Option<PathBuf>,
    pub params_txt: PathBuf,
    pub materialized: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpadesDiagnostics {
    executable: String,
    version: Option<String>,
    command: String,
    output_dir: String,
    stdout_log: String,
    stderr_log: String,
    elapsed_seconds: f64,
    exit_code: Option<i32>,
    success: bool,
    error: Option<String>,
}

#[derive(Debug, Clone, Copy)]
struct SpadesRunConfig<'a> {
    r1: &'a Path,
    r2: &'a Path,
    outdir: &'a Path,
    output_dir: &'a str,
    stdout_log_path: &'a Path,
    stderr_log_path: &'a Path,
    memory_gb: Option<usize>,
    tmp_dir: Option<&'a Path>,
    gfa11: bool,
}

pub fn stage_backend(
    args: &AssembleArgs,
    reads: &ReadArtifacts,
    outdir: &Path,
) -> Result<BackendArtifacts> {
    fs::create_dir_all(outdir)?;

    let spades_dir = contract::dir_path(outdir, contract::SPADES_DIR);
    let precreate_spades_dir = args.dry_run || matches!(args.backend, AssembleBackend::Mock);
    if precreate_spades_dir {
        fs::create_dir_all(&spades_dir)?;
    }
    let stdout_log_path = spades_dir.join(contract::SPADES_STDOUT_LOG);
    let stderr_log_path = spades_dir.join(contract::SPADES_STDERR_LOG);
    let log_path = spades_dir.join(contract::SPADES_LOG);
    let diagnostics_json = spades_dir.join(contract::SPADES_DIAGNOSTICS_JSON);
    let version_txt = spades_dir.join(contract::SPADES_VERSION_TXT);
    let temp_stdout_log_path = outdir.join(".spades.stdout.log.tmp");
    let temp_stderr_log_path = outdir.join(".spades.stderr.log.tmp");
    let temp_log_path = outdir.join(".spades.log.tmp");
    let temp_diagnostics_json = outdir.join(".spades_diagnostics.json.tmp");
    let temp_version_txt = outdir.join(".spades.version.txt.tmp");
    let contigs_fasta = spades_dir.join(contract::SPADES_CONTIGS_FASTA);
    let scaffolds_fasta = spades_dir.join(contract::SPADES_SCAFFOLDS_FASTA);
    let assembly_graph_gfa = spades_dir.join(contract::SPADES_GFA);
    let assembly_graph_fastg = spades_dir.join(contract::SPADES_FASTG);
    let contigs_paths = spades_dir.join(contract::SPADES_CONTIGS_PATHS);
    let scaffolds_paths = Some(spades_dir.join(contract::SPADES_SCAFFOLDS_PATHS));
    let params_txt = spades_dir.join(contract::SPADES_PARAMS_TXT);

    let executable = match args.backend {
        AssembleBackend::Spades => args.spades_executable.clone(),
        AssembleBackend::Mock => "mock-spades".to_string(),
    };
    let version = match args.backend {
        AssembleBackend::Spades => detect_spades_version(&args.spades_executable).ok(),
        AssembleBackend::Mock => Some("mock-spades 0.1.0".to_string()),
    };
    if let Some(version) = version.as_ref() {
        let version_target = if precreate_spades_dir {
            &version_txt
        } else {
            &temp_version_txt
        };
        io::ensure_parent_dir(version_target)?;
        fs::write(version_target, format!("{version}\n"))?;
    }

    let plan = BackendPlan {
        backend: match args.backend {
            AssembleBackend::Spades => "spades".to_string(),
            AssembleBackend::Mock => "mock".to_string(),
        },
        executable: executable.clone(),
        command: build_spades_command_line(&executable, args, reads, "spades"),
        version: version.clone(),
        dry_run: args.dry_run,
        k: args.k,
        threads: args.threads,
        memory_gb: args.spades_memory_gb,
        tmp_dir: args
            .spades_tmp_dir
            .as_ref()
            .map(|value| value.display().to_string()),
        gfa11: args.spades_gfa11,
        input_r1: reads.optimized_r1.display().to_string(),
        input_r2: reads.optimized_r2.display().to_string(),
        output_dir: spades_dir.display().to_string(),
        stdout_log: stdout_log_path.display().to_string(),
        stderr_log: stderr_log_path.display().to_string(),
        diagnostics_json: diagnostics_json.display().to_string(),
        generated_files: vec![
            contract::SPADES_COMMAND_SH.to_string(),
            contract::SPADES_PLAN_JSON.to_string(),
            contract::SPADES_LOG.to_string(),
            contract::SPADES_STDOUT_LOG.to_string(),
            contract::SPADES_STDERR_LOG.to_string(),
            contract::SPADES_DIAGNOSTICS_JSON.to_string(),
            contract::SPADES_VERSION_TXT.to_string(),
            contract::SPADES_PARAMS_TXT.to_string(),
            contract::SPADES_CONTIGS_FASTA.to_string(),
            contract::SPADES_SCAFFOLDS_FASTA.to_string(),
            contract::SPADES_GFA.to_string(),
            contract::SPADES_FASTG.to_string(),
            contract::SPADES_CONTIGS_PATHS.to_string(),
            contract::SPADES_SCAFFOLDS_PATHS.to_string(),
        ],
    };

    let plan_script = contract::root_path(outdir, contract::SPADES_COMMAND_SH);
    write_command_script(&plan_script, &plan.command, outdir)?;

    let plan_json = contract::root_path(outdir, contract::SPADES_PLAN_JSON);
    io::write_json(&plan_json, &plan)?;

    if args.dry_run {
        fs::write(&stdout_log_path, "")?;
        fs::write(&stderr_log_path, "")?;
        write_combined_log(&stdout_log_path, &stderr_log_path, &log_path)?;
        let diagnostics = SpadesDiagnostics {
            executable: plan.executable.clone(),
            version: plan.version.clone(),
            command: plan.command.clone(),
            output_dir: plan.output_dir.clone(),
            stdout_log: stdout_log_path.display().to_string(),
            stderr_log: stderr_log_path.display().to_string(),
            elapsed_seconds: 0.0,
            exit_code: None,
            success: true,
            error: None,
        };
        io::write_json(&diagnostics_json, &diagnostics)?;
        return Ok(BackendArtifacts {
            plan,
            plan_script,
            plan_json,
            spades_dir,
            log_path,
            stdout_log_path,
            stderr_log_path,
            diagnostics_json,
            version_txt,
            contigs_fasta,
            scaffolds_fasta,
            assembly_graph_gfa,
            assembly_graph_fastg,
            contigs_paths,
            scaffolds_paths,
            params_txt,
            materialized: false,
        });
    }

    let started = Instant::now();
    let diagnostics = match args.backend {
        AssembleBackend::Spades => {
            let run_config = SpadesRunConfig {
                r1: &reads.optimized_r1,
                r2: &reads.optimized_r2,
                outdir,
                output_dir: "spades",
                stdout_log_path: &temp_stdout_log_path,
                stderr_log_path: &temp_stderr_log_path,
                memory_gb: args.spades_memory_gb,
                tmp_dir: args.spades_tmp_dir.as_deref(),
                gfa11: args.spades_gfa11,
            };
            let status_result = run_spades(&plan, &run_config)?;
            write_combined_log(&temp_stdout_log_path, &temp_stderr_log_path, &temp_log_path)?;
            let diagnostics = SpadesDiagnostics {
                executable: plan.executable.clone(),
                version: plan.version.clone(),
                command: plan.command.clone(),
                output_dir: plan.output_dir.clone(),
                stdout_log: stdout_log_path.display().to_string(),
                stderr_log: stderr_log_path.display().to_string(),
                elapsed_seconds: started.elapsed().as_secs_f64(),
                exit_code: status_result.code(),
                success: status_result.success(),
                error: None,
            };
            if !status_result.success() {
                fs::rename(&temp_stdout_log_path, &stdout_log_path).ok();
                fs::rename(&temp_stderr_log_path, &stderr_log_path).ok();
                fs::rename(&temp_log_path, &log_path).ok();
                fs::rename(&temp_version_txt, &version_txt).ok();
                io::write_json(&temp_diagnostics_json, &diagnostics)?;
                fs::rename(&temp_diagnostics_json, &diagnostics_json).ok();
                bail!(
                    "SPAdes exited with status {} after {:.2}s; see {}, {} and {}",
                    status_result,
                    diagnostics.elapsed_seconds,
                    stdout_log_path.display(),
                    stderr_log_path.display(),
                    diagnostics_json.display()
                );
            }
            diagnostics
        }
        AssembleBackend::Mock => {
            materialize_mock_outputs(&spades_dir, args.spades_gfa11)?;
            fs::write(
                &stdout_log_path,
                "mock backend materialized SPAdes-like outputs\n",
            )?;
            fs::write(&stderr_log_path, "")?;
            write_combined_log(&stdout_log_path, &stderr_log_path, &log_path)?;
            SpadesDiagnostics {
                executable: plan.executable.clone(),
                version: plan.version.clone(),
                command: plan.command.clone(),
                output_dir: plan.output_dir.clone(),
                stdout_log: stdout_log_path.display().to_string(),
                stderr_log: stderr_log_path.display().to_string(),
                elapsed_seconds: started.elapsed().as_secs_f64(),
                exit_code: Some(0),
                success: true,
                error: None,
            }
        }
    };
    io::write_json(&temp_diagnostics_json, &diagnostics)?;
    if matches!(args.backend, AssembleBackend::Spades) {
        fs::rename(&temp_stdout_log_path, &stdout_log_path).ok();
        fs::rename(&temp_stderr_log_path, &stderr_log_path).ok();
        fs::rename(&temp_log_path, &log_path).ok();
        fs::rename(&temp_diagnostics_json, &diagnostics_json).ok();
        fs::rename(&temp_version_txt, &version_txt).ok();
    } else {
        fs::rename(&temp_diagnostics_json, &diagnostics_json).ok();
    }

    let scaffolds_paths = optional_file(&spades_dir.join(contract::SPADES_SCAFFOLDS_PATHS))?;
    Ok(BackendArtifacts {
        plan,
        plan_script,
        plan_json,
        spades_dir,
        log_path,
        stdout_log_path,
        stderr_log_path,
        diagnostics_json,
        version_txt,
        contigs_fasta: require_file(&contigs_fasta)?,
        scaffolds_fasta: require_file(&scaffolds_fasta)?,
        assembly_graph_gfa: require_file(&assembly_graph_gfa)?,
        assembly_graph_fastg: require_file(&assembly_graph_fastg)?,
        contigs_paths: require_file(&contigs_paths)?,
        scaffolds_paths,
        params_txt: require_file(&params_txt)?,
        materialized: true,
    })
}

fn build_spades_command_line(
    executable: &str,
    args: &AssembleArgs,
    reads: &ReadArtifacts,
    output_dir: &str,
) -> String {
    let r1 = absolute_path(&reads.optimized_r1).unwrap_or_else(|_| reads.optimized_r1.clone());
    let r2 = absolute_path(&reads.optimized_r2).unwrap_or_else(|_| reads.optimized_r2.clone());
    let mut command = format!(
        "{exe} --isolate -1 {r1} -2 {r2} -o {out} -t {threads} -k {k}",
        exe = shell_quote(executable),
        r1 = shell_quote(&r1.display().to_string()),
        r2 = shell_quote(&r2.display().to_string()),
        out = shell_quote(output_dir),
        threads = args.threads,
        k = args.k
    );
    if let Some(memory) = args.spades_memory_gb {
        command.push_str(&format!(" -m {memory}"));
    }
    if let Some(tmp_dir) = args.spades_tmp_dir.as_ref() {
        command.push_str(&format!(
            " --tmp-dir {}",
            shell_quote(&tmp_dir.display().to_string())
        ));
    }
    if args.spades_gfa11 {
        command.push_str(" --gfa11");
    }
    command
}

fn write_command_script(script_path: &Path, command: &str, outdir: &Path) -> Result<()> {
    io::ensure_parent_dir(script_path)?;
    let file = fs::File::create(script_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "#!/usr/bin/env bash")?;
    writeln!(writer, "set -euo pipefail")?;
    writeln!(writer, "cd {}", shell_quote(&outdir.display().to_string()))?;
    writeln!(writer, "exec {}", command)?;
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

fn run_spades(
    plan: &BackendPlan,
    config: &SpadesRunConfig<'_>,
) -> Result<std::process::ExitStatus> {
    io::ensure_parent_dir(config.stdout_log_path)?;
    io::ensure_parent_dir(config.stderr_log_path)?;
    let stdout_file = fs::File::create(config.stdout_log_path)?;
    let stderr_file = fs::File::create(config.stderr_log_path)?;

    let r1 = absolute_path(config.r1)?;
    let r2 = absolute_path(config.r2)?;
    let mut command = Command::new(&plan.executable);
    command
        .current_dir(config.outdir)
        .arg("--isolate")
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(config.output_dir)
        .arg("-t")
        .arg(plan.threads.to_string())
        .arg("-k")
        .arg(plan.k.to_string())
        .stdout(Stdio::from(stdout_file))
        .stderr(Stdio::from(stderr_file));
    if let Some(memory) = config.memory_gb {
        command.arg("-m").arg(memory.to_string());
    }
    if let Some(tmp_dir) = config.tmp_dir {
        command.arg("--tmp-dir").arg(tmp_dir);
    }
    if config.gfa11 {
        command.arg("--gfa11");
    }

    command
        .status()
        .with_context(|| format!("failed to execute SPAdes command: {}", plan.command))
}

fn detect_spades_version(executable: &str) -> Result<String> {
    let output = Command::new(executable)
        .arg("--version")
        .output()
        .with_context(|| format!("failed to query SPAdes version from {executable}"))?;
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
    let text = if !stdout.is_empty() { stdout } else { stderr };
    if text.is_empty() {
        bail!("SPAdes version command produced empty output");
    }
    Ok(text)
}

fn materialize_mock_outputs(spades_dir: &Path, gfa11: bool) -> Result<()> {
    fs::create_dir_all(spades_dir)?;
    let gfa_version = if gfa11 { "1.1" } else { "1.2" };
    let seq1 = "ACGT".repeat(60);
    let seq2 = "CGTA".repeat(60);
    fs::write(
        spades_dir.join(contract::SPADES_CONTIGS_FASTA),
        format!(">NODE_1_length_480_cov_1.0\n{}{}\n", seq1, seq2),
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_SCAFFOLDS_FASTA),
        format!(">NODE_1_length_480_cov_1.0\n{}{}\n", seq1, seq2),
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_GFA),
        format!(
            "H\tVN:Z:{gfa_version}\nS\t1\t{seq1}\tLN:i:240\tKC:i:240\tdp:f:1.0\nS\t2\t{seq2}\tLN:i:240\tKC:i:240\tdp:f:1.0\nL\t1\t+\t2\t+\t3M\nP\tNODE_1_length_480_cov_1.0\t1+,2+\t*\n"
        ),
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_FASTG),
        format!(">EDGE_1_length_240_cov_1.0\n{}\n", seq1),
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_CONTIGS_PATHS),
        "NODE_1_length_480_cov_1.0\n1+,2+;2-,1-\n",
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_SCAFFOLDS_PATHS),
        "NODE_1_length_480_cov_1.0\n1+,2+\n",
    )?;
    fs::write(
        spades_dir.join(contract::SPADES_PARAMS_TXT),
        "mock params\n",
    )?;
    Ok(())
}

fn absolute_path(path: &Path) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

fn write_combined_log(stdout_path: &Path, stderr_path: &Path, combined_path: &Path) -> Result<()> {
    io::ensure_parent_dir(combined_path)?;
    let stdout_text = fs::read_to_string(stdout_path).unwrap_or_default();
    let stderr_text = fs::read_to_string(stderr_path).unwrap_or_default();
    let file = fs::File::create(combined_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "## stdout")?;
    write!(writer, "{stdout_text}")?;
    if !stdout_text.ends_with('\n') && !stdout_text.is_empty() {
        writeln!(writer)?;
    }
    writeln!(writer, "## stderr")?;
    write!(writer, "{stderr_text}")?;
    if !stderr_text.ends_with('\n') && !stderr_text.is_empty() {
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn require_file(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        Ok(path.to_path_buf())
    } else {
        bail!("expected SPAdes output file is missing: {}", path.display())
    }
}

fn optional_file(path: &Path) -> Result<Option<PathBuf>> {
    if path.exists() {
        Ok(Some(path.to_path_buf()))
    } else {
        Ok(None)
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
