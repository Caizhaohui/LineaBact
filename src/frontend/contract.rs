use std::path::{Path, PathBuf};

pub const ROOT_BACKBONE_CONTIGS: &str = "backbone_contigs.fasta";
pub const ROOT_FINISHED_CONTIGS: &str = "finished_contigs.fasta";
pub const ROOT_ASSEMBLY_GRAPH: &str = "assembly_graph.gfa";
pub const ROOT_ASSEMBLY_STATS: &str = "assembly_stats.tsv";
pub const ROOT_PARAMS: &str = "params.toml";
pub const ROOT_RUN_SUMMARY: &str = "run_summary.json";

pub const READS_DIR: &str = "reads";
pub const READS_OPTIMIZED_R1: &str = "optimized_R1.fastq.gz";
pub const READS_OPTIMIZED_R2: &str = "optimized_R2.fastq.gz";
pub const READS_SUMMARY_JSON: &str = "reads_summary.json";
pub const READS_STATS_TSV: &str = "reads_stats.tsv";
pub const READS_PAIR_CHECK_TSV: &str = "pair_check.tsv";
pub const READS_DOWNSAMPLE_PLAN_TSV: &str = "downsample_plan.tsv";
pub const READS_TRIM_PLAN_TSV: &str = "trim_plan.tsv";
pub const READS_OUTPUT_FILES: &[&str] = &[
    READS_OPTIMIZED_R1,
    READS_OPTIMIZED_R2,
    READS_SUMMARY_JSON,
    READS_STATS_TSV,
    READS_PAIR_CHECK_TSV,
    READS_DOWNSAMPLE_PLAN_TSV,
    READS_TRIM_PLAN_TSV,
];

pub const SPADES_DIR: &str = "spades";
pub const SPADES_COMMAND_SH: &str = "spades.command.sh";
pub const SPADES_PLAN_JSON: &str = "spades_plan.json";
pub const SPADES_LOG: &str = "spades.log";
pub const SPADES_STDOUT_LOG: &str = "spades.stdout.log";
pub const SPADES_STDERR_LOG: &str = "spades.stderr.log";
pub const SPADES_DIAGNOSTICS_JSON: &str = "spades_diagnostics.json";
pub const SPADES_VERSION_TXT: &str = "spades.version.txt";
pub const SPADES_PARAMS_TXT: &str = "params.txt";
pub const SPADES_CONTIGS_FASTA: &str = "contigs.fasta";
pub const SPADES_SCAFFOLDS_FASTA: &str = "scaffolds.fasta";
pub const SPADES_GFA: &str = "assembly_graph_with_scaffolds.gfa";
pub const SPADES_FASTG: &str = "assembly_graph.fastg";
pub const SPADES_CONTIGS_PATHS: &str = "contigs.paths";
pub const SPADES_SCAFFOLDS_PATHS: &str = "scaffolds.paths";
pub const SPADES_OUTPUT_FILES: &[&str] = &[
    SPADES_COMMAND_SH,
    SPADES_PLAN_JSON,
    SPADES_LOG,
    SPADES_STDOUT_LOG,
    SPADES_STDERR_LOG,
    SPADES_DIAGNOSTICS_JSON,
    SPADES_VERSION_TXT,
    SPADES_PARAMS_TXT,
    SPADES_CONTIGS_FASTA,
    SPADES_SCAFFOLDS_FASTA,
    SPADES_GFA,
    SPADES_FASTG,
    SPADES_CONTIGS_PATHS,
    SPADES_SCAFFOLDS_PATHS,
];

pub const GRAPH_DIR: &str = "graph";
pub const GRAPH_GFA: &str = "assembly_graph.gfa";
pub const GRAPH_FASTG: &str = "assembly_graph.fastg";
pub const GRAPH_CONTIGS_PATHS: &str = "contigs.paths";
pub const GRAPH_SCAFFOLDS_PATHS: &str = "scaffolds.paths";
pub const GRAPH_SEGMENTS_TSV: &str = "segments.tsv";
pub const GRAPH_LINKS_TSV: &str = "links.tsv";
pub const GRAPH_DEPTH_TSV: &str = "depth.tsv";
pub const GRAPH_PATHS_TSV: &str = "paths.tsv";
pub const GRAPH_ANCHOR_SEGMENTS_TSV: &str = "anchor_segments.tsv";
pub const GRAPH_QC_TSV: &str = "graph_qc.tsv";
pub const GRAPH_PATHS_SUMMARY_TSV: &str = "paths_summary.tsv";
pub const GRAPH_SUMMARY_JSON: &str = "graph_summary.json";
pub const GRAPH_OUTPUT_FILES: &[&str] = &[
    GRAPH_GFA,
    GRAPH_FASTG,
    GRAPH_CONTIGS_PATHS,
    GRAPH_SCAFFOLDS_PATHS,
    GRAPH_SEGMENTS_TSV,
    GRAPH_LINKS_TSV,
    GRAPH_DEPTH_TSV,
    GRAPH_PATHS_TSV,
    GRAPH_ANCHOR_SEGMENTS_TSV,
    GRAPH_QC_TSV,
    GRAPH_PATHS_SUMMARY_TSV,
    GRAPH_SUMMARY_JSON,
];

pub const BRIDGING_DIR: &str = "bridging";
pub const BRIDGING_MANIFEST_JSON: &str = "bridge_manifest.json";
pub const BRIDGING_CANDIDATES_TSV: &str = "spades_path_bridge_candidates.tsv";
pub const BRIDGING_CANDIDATES_JSONL: &str = "bridge_evidence.jsonl";
pub const BRIDGING_SUMMARY_JSON: &str = "bridge_summary.json";
pub const BRIDGING_CONFLICTS_TSV: &str = "bridge_conflicts.tsv";
pub const BRIDGING_DECISIONS_JSONL: &str = "bridge_decisions.jsonl";
pub const BRIDGING_GRAPH_GFA: &str = "bridged_graph.gfa";
pub const BRIDGING_LEGACY_CANDIDATES_TSV: &str = "bridge_candidates.tsv";
pub const BRIDGING_LEGACY_CANDIDATES_JSONL: &str = "bridge_candidates.jsonl";
pub const BRIDGING_OUTPUT_FILES: &[&str] = &[
    BRIDGING_MANIFEST_JSON,
    BRIDGING_CANDIDATES_TSV,
    BRIDGING_CANDIDATES_JSONL,
    BRIDGING_SUMMARY_JSON,
    BRIDGING_CONFLICTS_TSV,
    BRIDGING_DECISIONS_JSONL,
    BRIDGING_GRAPH_GFA,
    BRIDGING_LEGACY_CANDIDATES_TSV,
    BRIDGING_LEGACY_CANDIDATES_JSONL,
];

pub const POSTPROCESS_DIR: &str = "postprocess";
pub const POSTPROCESS_CONTIG_STATS_TSV: &str = "contig_stats.tsv";
pub const POSTPROCESS_RENAME_MAP_TSV: &str = "rename_map.tsv";
pub const POSTPROCESS_BACKBONE_FASTA: &str = "backbone_contigs.filtered.fasta";
pub const POSTPROCESS_FINISHED_FASTA: &str = "finished_contigs.filtered.fasta";
pub const POSTPROCESS_OUTPUT_FILES: &[&str] = &[
    POSTPROCESS_CONTIG_STATS_TSV,
    POSTPROCESS_RENAME_MAP_TSV,
    POSTPROCESS_BACKBONE_FASTA,
    POSTPROCESS_FINISHED_FASTA,
];

pub const ROOT_OUTPUT_FILES: &[&str] = &[
    ROOT_BACKBONE_CONTIGS,
    ROOT_FINISHED_CONTIGS,
    ROOT_ASSEMBLY_GRAPH,
    ROOT_ASSEMBLY_STATS,
    ROOT_PARAMS,
    ROOT_RUN_SUMMARY,
];

pub const STABLE_OUTPUT_DIRS: &[&str] = &[READS_DIR, SPADES_DIR, GRAPH_DIR, BRIDGING_DIR];

pub fn root_path(root: &Path, file_name: &str) -> PathBuf {
    root.join(file_name)
}

pub fn dir_path(root: &Path, dir_name: &str) -> PathBuf {
    root.join(dir_name)
}
