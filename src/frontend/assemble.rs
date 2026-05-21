use anyhow::Result;
use log::info;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::cli::AssembleArgs;
use crate::io;

use super::backend::{BackendArtifacts, stage_backend};
use super::bridging::stage_bridging;
use super::contract;
use super::graph::{GraphArtifacts, stage_graph};
use super::postprocess::{PostprocessArtifacts, stage_postprocess};
use super::reads::{ReadPreparationConfig, prepare_reads};
use super::report::{AssemblyStats, ParamsSnapshot, RunSummary};

#[derive(Debug, Clone)]
struct OutputPaths {
    backbone_contigs: PathBuf,
    finished_contigs: PathBuf,
    assembly_graph: PathBuf,
    assembly_stats: PathBuf,
    params_toml: PathBuf,
    run_summary_json: PathBuf,
}

pub fn run(args: AssembleArgs) -> Result<()> {
    let start = Instant::now();
    info!(
        "LineaBact short-read front-end v{}",
        env!("CARGO_PKG_VERSION")
    );

    fs::create_dir_all(&args.outdir)?;
    let reads_dir = contract::dir_path(&args.outdir, contract::READS_DIR);
    let graph_dir = contract::dir_path(&args.outdir, contract::GRAPH_DIR);
    let bridging_dir = contract::dir_path(&args.outdir, contract::BRIDGING_DIR);
    let postprocess_dir = contract::dir_path(&args.outdir, contract::POSTPROCESS_DIR);

    let read_config = ReadPreparationConfig {
        downsample_fraction: args.downsample_fraction,
        target_coverage: args.target_coverage,
        genome_size_bp: args.genome_size_bp,
        downsample_seed: args.downsample_seed,
        trim_adapters: args.trim_adapters,
        trim_tool: &args.trim_tool,
        read_sha256: args.read_sha256,
    };
    let read_artifacts = prepare_reads(&args.r1, &args.r2, &reads_dir, &read_config)?;
    info!(
        "Prepared reads: {} retained pairs -> {}",
        read_artifacts.report.retained_pairs,
        read_artifacts.optimized_r1.display()
    );

    let backend_artifacts = stage_backend(&args, &read_artifacts, &args.outdir)?;
    info!(
        "Backend plan recorded at {}",
        backend_artifacts.plan_script.display()
    );
    if args.dry_run {
        return write_dry_run_outputs(&args, start, &read_artifacts, &backend_artifacts);
    }
    info!(
        "Backend outputs staged in {}",
        backend_artifacts.spades_dir.display()
    );
    info!("Backend log at {}", backend_artifacts.log_path.display());

    let graph_artifacts = stage_graph(&backend_artifacts.spades_dir, &graph_dir)?;
    info!(
        "Materialized graph artifacts at {}",
        graph_artifacts.normalized_gfa.display()
    );

    let bridge_artifacts = stage_bridging(
        &graph_artifacts.normalized_gfa,
        &graph_artifacts.contigs_paths,
        graph_artifacts.scaffolds_paths.as_deref(),
        &bridging_dir,
    )?;
    info!(
        "Generated {} bridge candidates",
        bridge_artifacts.summary.candidate_count
    );

    let postprocess_artifacts = stage_postprocess(
        &backend_artifacts.contigs_fasta,
        &backend_artifacts.scaffolds_fasta,
        &postprocess_dir,
    )?;
    info!(
        "Postprocessed contigs written to {}",
        postprocess_artifacts.finished_fasta.display()
    );

    let output_paths = write_root_outputs(&args, &postprocess_artifacts, &graph_artifacts)?;

    let assembly_stats = AssemblyStats::from_fasta(&output_paths.finished_contigs)?;
    assembly_stats.write_tsv(&output_paths.assembly_stats)?;

    let params = ParamsSnapshot {
        version: env!("CARGO_PKG_VERSION").to_string(),
        command: "lineabact assemble".to_string(),
        backend: backend_artifacts.plan.backend.clone(),
        dry_run: false,
        input_r1: args.r1.display().to_string(),
        input_r2: args.r2.display().to_string(),
        outdir: args.outdir.display().to_string(),
        spades_executable: args.spades_executable.clone(),
        k: args.k,
        threads: args.threads,
        spades_memory_gb: args.spades_memory_gb,
        spades_tmp_dir: args
            .spades_tmp_dir
            .as_ref()
            .map(|value| value.display().to_string()),
        spades_gfa11: args.spades_gfa11,
        downsample_fraction: args.downsample_fraction,
        target_coverage: args.target_coverage,
        genome_size_bp: args.genome_size_bp,
        downsample_seed: args.downsample_seed,
        trim_adapters: args.trim_adapters,
        trim_tool: args.trim_tool.clone(),
        read_sha256: args.read_sha256,
        read_summary: read_artifacts.summary_json.display().to_string(),
        backend_plan: backend_artifacts.plan_json.display().to_string(),
        graph_summary: Some(graph_artifacts.summary_json.display().to_string()),
        bridge_summary: Some(bridge_artifacts.summary_json.display().to_string()),
    };
    params.write_toml(&output_paths.params_toml)?;

    let elapsed = start.elapsed().as_secs_f64();
    let output_files = collect_full_output_files(
        &args.outdir,
        &output_paths,
        &read_artifacts,
        &backend_artifacts,
        &graph_artifacts,
        &postprocess_artifacts,
        &bridge_artifacts,
    );

    let summary = RunSummary {
        version: env!("CARGO_PKG_VERSION").to_string(),
        command: "lineabact assemble".to_string(),
        backend: backend_artifacts.plan.backend.clone(),
        dry_run: false,
        input_r1: args.r1.display().to_string(),
        input_r2: args.r2.display().to_string(),
        outdir: args.outdir.display().to_string(),
        spades_executable: args.spades_executable.clone(),
        k: args.k,
        threads: args.threads,
        spades_memory_gb: args.spades_memory_gb,
        spades_tmp_dir: args
            .spades_tmp_dir
            .as_ref()
            .map(|value| value.display().to_string()),
        spades_gfa11: args.spades_gfa11,
        downsample_fraction: args.downsample_fraction,
        target_coverage: args.target_coverage,
        genome_size_bp: args.genome_size_bp,
        downsample_seed: args.downsample_seed,
        trim_adapters: args.trim_adapters,
        trim_tool: args.trim_tool.clone(),
        read_sha256: args.read_sha256,
        elapsed_seconds: elapsed,
        read_report: read_artifacts.report.clone(),
        backend_plan: backend_artifacts.plan.clone(),
        graph_summary: Some(graph_artifacts.summary.clone()),
        bridge_summary: Some(bridge_artifacts.summary.clone()),
        assembly_stats: Some(assembly_stats.clone()),
        output_files: output_files.clone(),
    };
    summary.write_json(&output_paths.run_summary_json)?;

    info!("Wrote {} output files", output_files.len());
    info!("Assembly complete in {:.2}s", elapsed);
    info!("N50: {} bp", assembly_stats.n50);
    info!("Total length: {} bp", assembly_stats.total_length);

    Ok(())
}

fn write_dry_run_outputs(
    args: &AssembleArgs,
    start: Instant,
    read_artifacts: &super::reads::ReadArtifacts,
    backend_artifacts: &BackendArtifacts,
) -> Result<()> {
    let output_paths = output_paths(&args.outdir);
    let params = ParamsSnapshot {
        version: env!("CARGO_PKG_VERSION").to_string(),
        command: "lineabact assemble".to_string(),
        backend: backend_artifacts.plan.backend.clone(),
        dry_run: true,
        input_r1: args.r1.display().to_string(),
        input_r2: args.r2.display().to_string(),
        outdir: args.outdir.display().to_string(),
        spades_executable: args.spades_executable.clone(),
        k: args.k,
        threads: args.threads,
        spades_memory_gb: args.spades_memory_gb,
        spades_tmp_dir: args
            .spades_tmp_dir
            .as_ref()
            .map(|value| value.display().to_string()),
        spades_gfa11: args.spades_gfa11,
        downsample_fraction: args.downsample_fraction,
        target_coverage: args.target_coverage,
        genome_size_bp: args.genome_size_bp,
        downsample_seed: args.downsample_seed,
        trim_adapters: args.trim_adapters,
        trim_tool: args.trim_tool.clone(),
        read_sha256: args.read_sha256,
        read_summary: read_artifacts.summary_json.display().to_string(),
        backend_plan: backend_artifacts.plan_json.display().to_string(),
        graph_summary: None,
        bridge_summary: None,
    };
    params.write_toml(&output_paths.params_toml)?;

    let output_files = vec![
        relative(&args.outdir, &read_artifacts.summary_json),
        relative(&args.outdir, &read_artifacts.optimized_r1),
        relative(&args.outdir, &read_artifacts.optimized_r2),
        relative(&args.outdir, &read_artifacts.reads_stats_tsv),
        relative(&args.outdir, &read_artifacts.pair_check_tsv),
        relative(&args.outdir, &read_artifacts.downsample_plan_tsv),
        relative(&args.outdir, &read_artifacts.trim_plan_tsv),
        relative(&args.outdir, &backend_artifacts.plan_script),
        relative(&args.outdir, &backend_artifacts.plan_json),
        relative(&args.outdir, &backend_artifacts.log_path),
        relative(&args.outdir, &backend_artifacts.stdout_log_path),
        relative(&args.outdir, &backend_artifacts.stderr_log_path),
        relative(&args.outdir, &backend_artifacts.diagnostics_json),
        relative(&args.outdir, &backend_artifacts.version_txt),
        relative(&args.outdir, &output_paths.params_toml),
        relative(&args.outdir, &output_paths.run_summary_json),
    ];

    let summary = RunSummary {
        version: env!("CARGO_PKG_VERSION").to_string(),
        command: "lineabact assemble".to_string(),
        backend: backend_artifacts.plan.backend.clone(),
        dry_run: true,
        input_r1: args.r1.display().to_string(),
        input_r2: args.r2.display().to_string(),
        outdir: args.outdir.display().to_string(),
        spades_executable: args.spades_executable.clone(),
        k: args.k,
        threads: args.threads,
        spades_memory_gb: args.spades_memory_gb,
        spades_tmp_dir: args
            .spades_tmp_dir
            .as_ref()
            .map(|value| value.display().to_string()),
        spades_gfa11: args.spades_gfa11,
        downsample_fraction: args.downsample_fraction,
        target_coverage: args.target_coverage,
        genome_size_bp: args.genome_size_bp,
        downsample_seed: args.downsample_seed,
        trim_adapters: args.trim_adapters,
        trim_tool: args.trim_tool.clone(),
        read_sha256: args.read_sha256,
        elapsed_seconds: start.elapsed().as_secs_f64(),
        read_report: read_artifacts.report.clone(),
        backend_plan: backend_artifacts.plan.clone(),
        graph_summary: None,
        bridge_summary: None,
        assembly_stats: None,
        output_files: output_files.clone(),
    };
    summary.write_json(&output_paths.run_summary_json)?;
    info!("Dry-run complete in {:.2}s", summary.elapsed_seconds);
    Ok(())
}

fn output_paths(outdir: &Path) -> OutputPaths {
    OutputPaths {
        backbone_contigs: contract::root_path(outdir, contract::ROOT_BACKBONE_CONTIGS),
        finished_contigs: contract::root_path(outdir, contract::ROOT_FINISHED_CONTIGS),
        assembly_graph: contract::root_path(outdir, contract::ROOT_ASSEMBLY_GRAPH),
        assembly_stats: contract::root_path(outdir, contract::ROOT_ASSEMBLY_STATS),
        params_toml: contract::root_path(outdir, contract::ROOT_PARAMS),
        run_summary_json: contract::root_path(outdir, contract::ROOT_RUN_SUMMARY),
    }
}

fn collect_full_output_files(
    outdir: &Path,
    output_paths: &OutputPaths,
    read_artifacts: &super::reads::ReadArtifacts,
    backend_artifacts: &BackendArtifacts,
    graph_artifacts: &GraphArtifacts,
    postprocess_artifacts: &PostprocessArtifacts,
    bridge_artifacts: &super::bridging::BridgeArtifacts,
) -> Vec<String> {
    let mut output_files = vec![
        relative(outdir, &output_paths.backbone_contigs),
        relative(outdir, &output_paths.finished_contigs),
        relative(outdir, &output_paths.assembly_graph),
        relative(outdir, &output_paths.assembly_stats),
        relative(outdir, &output_paths.params_toml),
        relative(outdir, &output_paths.run_summary_json),
        relative(outdir, &read_artifacts.summary_json),
        relative(outdir, &read_artifacts.optimized_r1),
        relative(outdir, &read_artifacts.optimized_r2),
        relative(outdir, &read_artifacts.reads_stats_tsv),
        relative(outdir, &read_artifacts.pair_check_tsv),
        relative(outdir, &read_artifacts.downsample_plan_tsv),
        relative(outdir, &read_artifacts.trim_plan_tsv),
        relative(outdir, &backend_artifacts.plan_script),
        relative(outdir, &backend_artifacts.plan_json),
        relative(outdir, &backend_artifacts.log_path),
        relative(outdir, &backend_artifacts.stdout_log_path),
        relative(outdir, &backend_artifacts.stderr_log_path),
        relative(outdir, &backend_artifacts.diagnostics_json),
        relative(outdir, &backend_artifacts.version_txt),
        relative(outdir, &backend_artifacts.params_txt),
        relative(outdir, &backend_artifacts.contigs_fasta),
        relative(outdir, &backend_artifacts.scaffolds_fasta),
        relative(outdir, &backend_artifacts.assembly_graph_gfa),
        relative(outdir, &backend_artifacts.assembly_graph_fastg),
        relative(outdir, &backend_artifacts.contigs_paths),
        relative(outdir, &graph_artifacts.normalized_gfa),
        relative(outdir, &graph_artifacts.fastg),
        relative(outdir, &graph_artifacts.contigs_paths),
        relative(outdir, &graph_artifacts.summary_json),
        relative(outdir, &graph_artifacts.segments_tsv),
        relative(outdir, &graph_artifacts.links_tsv),
        relative(outdir, &graph_artifacts.depth_tsv),
        relative(outdir, &graph_artifacts.paths_tsv),
        relative(outdir, &graph_artifacts.anchor_segments_tsv),
        relative(outdir, &graph_artifacts.graph_qc_tsv),
        relative(outdir, &graph_artifacts.paths_summary_tsv),
        relative(outdir, &postprocess_artifacts.backbone_fasta),
        relative(outdir, &postprocess_artifacts.finished_fasta),
        relative(outdir, &postprocess_artifacts.contig_stats_tsv),
        relative(outdir, &postprocess_artifacts.rename_map_tsv),
        relative(outdir, &bridge_artifacts.manifest_json),
        relative(outdir, &bridge_artifacts.candidates_tsv),
        relative(outdir, &bridge_artifacts.candidates_jsonl),
        relative(outdir, &bridge_artifacts.conflicts_tsv),
        relative(outdir, &bridge_artifacts.decisions_jsonl),
        relative(outdir, &bridge_artifacts.bridged_graph_gfa),
        relative(outdir, &bridge_artifacts.summary_json),
    ];
    if let Some(path) = backend_artifacts.scaffolds_paths.as_ref() {
        output_files.push(relative(outdir, path));
    }
    if let Some(path) = graph_artifacts.scaffolds_paths.as_ref() {
        output_files.push(relative(outdir, path));
    }
    output_files
}

fn write_root_outputs(
    args: &AssembleArgs,
    postprocess_artifacts: &PostprocessArtifacts,
    graph_artifacts: &GraphArtifacts,
) -> Result<OutputPaths> {
    let output_paths = output_paths(&args.outdir);
    io::copy_file(
        &postprocess_artifacts.backbone_fasta,
        &output_paths.backbone_contigs,
    )?;
    io::copy_file(
        &postprocess_artifacts.finished_fasta,
        &output_paths.finished_contigs,
    )?;
    io::copy_file(
        &graph_artifacts.normalized_gfa,
        &output_paths.assembly_graph,
    )?;
    Ok(output_paths)
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{AssembleArgs, AssembleBackend};
    use crate::frontend::contract;
    use crate::frontend::fixture::extract_zip_entry;
    use serde_json::Value;
    use tempfile::tempdir;

    #[test]
    fn stable_output_contract_is_explicit() {
        assert_eq!(
            contract::ROOT_OUTPUT_FILES,
            [
                contract::ROOT_BACKBONE_CONTIGS,
                contract::ROOT_FINISHED_CONTIGS,
                contract::ROOT_ASSEMBLY_GRAPH,
                contract::ROOT_ASSEMBLY_STATS,
                contract::ROOT_PARAMS,
                contract::ROOT_RUN_SUMMARY,
            ]
        );
        assert_eq!(
            contract::STABLE_OUTPUT_DIRS,
            [
                contract::READS_DIR,
                contract::SPADES_DIR,
                contract::GRAPH_DIR,
                contract::BRIDGING_DIR
            ]
        );
    }

    #[test]
    fn assemble_smoke_uses_unicycler_fixtures() {
        let temp = tempdir().unwrap();
        let inputs = temp.path().join("inputs");
        let outdir = temp.path().join("out");
        fs::create_dir_all(&inputs).unwrap();

        let archive = PathBuf::from("reference_tools/Unicycler-main.zip");
        let r1 = inputs.join("short_reads_1.fastq.gz");
        let r2 = inputs.join("short_reads_2.fastq.gz");
        extract_zip_entry(
            &archive,
            "Unicycler-main/sample_data/short_reads_1.fastq.gz",
            &r1,
        )
        .unwrap();
        extract_zip_entry(
            &archive,
            "Unicycler-main/sample_data/short_reads_2.fastq.gz",
            &r2,
        )
        .unwrap();

        let args = AssembleArgs {
            r1,
            r2,
            outdir: outdir.clone(),
            backend: AssembleBackend::Mock,
            dry_run: false,
            k: 55,
            threads: 2,
            spades_memory_gb: Some(8),
            spades_tmp_dir: Some(outdir.join("tmp")),
            spades_gfa11: true,
            downsample_fraction: None,
            target_coverage: Some(80.0),
            genome_size_bp: Some(7_000_000),
            downsample_seed: 0,
            trim_adapters: true,
            trim_tool: "fastp".to_string(),
            read_sha256: false,
            spades_executable: "spades.py".to_string(),
        };

        run(args).unwrap();

        for relative_path in [
            "backbone_contigs.fasta",
            "finished_contigs.fasta",
            "assembly_graph.gfa",
            "assembly_stats.tsv",
            "params.toml",
            "run_summary.json",
            "reads/reads_summary.json",
            "reads/reads_stats.tsv",
            "reads/pair_check.tsv",
            "reads/downsample_plan.tsv",
            "reads/trim_plan.tsv",
            "spades.command.sh",
            "spades_plan.json",
            "spades/spades.log",
            "spades/spades.stdout.log",
            "spades/spades.stderr.log",
            "spades/spades_diagnostics.json",
            "spades/spades.version.txt",
            "spades/params.txt",
            "spades/contigs.fasta",
            "spades/scaffolds.fasta",
            "spades/assembly_graph_with_scaffolds.gfa",
            "spades/assembly_graph.fastg",
            "spades/contigs.paths",
            "spades/scaffolds.paths",
            "graph/assembly_graph.gfa",
            "graph/contigs.paths",
            "graph/scaffolds.paths",
            "graph/graph_summary.json",
            "graph/segments.tsv",
            "graph/links.tsv",
            "graph/depth.tsv",
            "graph/paths.tsv",
            "graph/anchor_segments.tsv",
            "graph/graph_qc.tsv",
            "graph/paths_summary.tsv",
            "postprocess/backbone_contigs.filtered.fasta",
            "postprocess/finished_contigs.filtered.fasta",
            "postprocess/contig_stats.tsv",
            "postprocess/rename_map.tsv",
            "bridging/bridge_manifest.json",
            "bridging/spades_path_bridge_candidates.tsv",
            "bridging/bridge_evidence.jsonl",
            "bridging/bridge_conflicts.tsv",
            "bridging/bridge_decisions.jsonl",
            "bridging/bridged_graph.gfa",
        ] {
            assert!(
                outdir.join(relative_path).exists(),
                "missing {relative_path}"
            );
        }

        let summary: Value =
            serde_json::from_str(&fs::read_to_string(outdir.join("run_summary.json")).unwrap())
                .unwrap();
        assert_eq!(summary["backend"], "mock");
        assert_eq!(summary["backend_plan"]["backend"], "mock");
        assert_eq!(summary["dry_run"], false);
        assert!(
            summary["bridge_summary"]["candidate_count"]
                .as_u64()
                .unwrap()
                > 0
        );
        assert!(summary["assembly_stats"]["num_contigs"].as_u64().unwrap() > 0);
    }

    #[test]
    fn assemble_dry_run_records_plan_without_materializing_backend_outputs() {
        let temp = tempdir().unwrap();
        let inputs = temp.path().join("inputs");
        let outdir = temp.path().join("out");
        fs::create_dir_all(&inputs).unwrap();

        let archive = PathBuf::from("reference_tools/Unicycler-main.zip");
        let r1 = inputs.join("short_reads_1.fastq.gz");
        let r2 = inputs.join("short_reads_2.fastq.gz");
        extract_zip_entry(
            &archive,
            "Unicycler-main/sample_data/short_reads_1.fastq.gz",
            &r1,
        )
        .unwrap();
        extract_zip_entry(
            &archive,
            "Unicycler-main/sample_data/short_reads_2.fastq.gz",
            &r2,
        )
        .unwrap();

        let args = AssembleArgs {
            r1,
            r2,
            outdir: outdir.clone(),
            backend: AssembleBackend::Spades,
            dry_run: true,
            k: 55,
            threads: 2,
            spades_memory_gb: Some(8),
            spades_tmp_dir: Some(outdir.join("tmp")),
            spades_gfa11: true,
            downsample_fraction: None,
            target_coverage: Some(80.0),
            genome_size_bp: Some(7_000_000),
            downsample_seed: 0,
            trim_adapters: true,
            trim_tool: "fastp".to_string(),
            read_sha256: false,
            spades_executable: "spades.py".to_string(),
        };

        run(args).unwrap();

        assert!(outdir.join("params.toml").exists());
        assert!(outdir.join("run_summary.json").exists());
        assert!(outdir.join("spades.command.sh").exists());
        assert!(outdir.join("spades_plan.json").exists());
        assert!(!outdir.join("spades/contigs.fasta").exists());
        assert!(!outdir.join("graph/assembly_graph.gfa").exists());

        let summary: Value =
            serde_json::from_str(&fs::read_to_string(outdir.join("run_summary.json")).unwrap())
                .unwrap();
        assert_eq!(summary["backend"], "spades");
        assert_eq!(summary["dry_run"], true);
        assert!(summary["graph_summary"].is_null());
        assert!(summary["bridge_summary"].is_null());
        assert!(summary["assembly_stats"].is_null());
    }
}
