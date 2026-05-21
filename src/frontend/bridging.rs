use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::io;

use super::contract;
use super::graph::{PathRecord, parse_paths_file};
use super::report::BridgeSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BridgeArtifacts {
    pub summary: BridgeSummary,
    pub manifest_json: PathBuf,
    pub candidates_tsv: PathBuf,
    pub candidates_jsonl: PathBuf,
    pub conflicts_tsv: PathBuf,
    pub decisions_jsonl: PathBuf,
    pub bridged_graph_gfa: PathBuf,
    pub summary_json: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BridgeCandidate {
    pub bridge_id: String,
    pub source: String,
    pub left_segment: String,
    pub right_segment: String,
    pub path_segments: Vec<String>,
    pub bridge_length: usize,
    pub depth_agreement: f64,
    pub path_self_contained: bool,
    pub insert_size_penalty: f64,
    pub quality: f64,
    pub status: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BridgeManifest {
    contigs_paths: String,
    scaffolds_paths: Option<String>,
    source_graph_gfa: String,
    bridged_graph_gfa: String,
    candidate_tsv: String,
    evidence_jsonl: String,
    conflict_tsv: String,
    decision_jsonl: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BridgeDecision {
    bridge_id: String,
    status: String,
    reason: String,
    applied: bool,
}

pub fn stage_bridging(
    graph_gfa: &Path,
    contigs_paths: &Path,
    scaffolds_paths: Option<&Path>,
    outdir: &Path,
) -> Result<BridgeArtifacts> {
    fs::create_dir_all(outdir)?;
    let manifest_json = contract::root_path(outdir, contract::BRIDGING_MANIFEST_JSON);
    let candidates_tsv = contract::root_path(outdir, contract::BRIDGING_CANDIDATES_TSV);
    let candidates_jsonl = contract::root_path(outdir, contract::BRIDGING_CANDIDATES_JSONL);
    let conflicts_tsv = contract::root_path(outdir, contract::BRIDGING_CONFLICTS_TSV);
    let decisions_jsonl = contract::root_path(outdir, contract::BRIDGING_DECISIONS_JSONL);
    let bridged_graph_gfa = contract::root_path(outdir, contract::BRIDGING_GRAPH_GFA);
    let summary_json = contract::root_path(outdir, contract::BRIDGING_SUMMARY_JSON);
    let legacy_candidates_tsv =
        contract::root_path(outdir, contract::BRIDGING_LEGACY_CANDIDATES_TSV);
    let legacy_candidates_jsonl =
        contract::root_path(outdir, contract::BRIDGING_LEGACY_CANDIDATES_JSONL);

    let mut path_records = parse_paths_file(contigs_paths)?;
    if let Some(scaffolds_paths) = scaffolds_paths {
        path_records.extend(parse_paths_file(scaffolds_paths)?);
    }
    let candidates = generate_candidates(&path_records);

    write_candidates_tsv(&candidates, &candidates_tsv)?;
    write_candidates_jsonl(&candidates, &candidates_jsonl)?;
    write_conflicts_tsv(&conflicts_tsv)?;
    write_decisions_jsonl(&candidates, &decisions_jsonl)?;
    io::copy_file(graph_gfa, &bridged_graph_gfa)?;

    let summary = BridgeSummary {
        path_record_count: path_records.len(),
        forward_path_records: path_records
            .iter()
            .filter(|record| !record.is_reverse)
            .count(),
        reverse_path_records: path_records
            .iter()
            .filter(|record| record.is_reverse)
            .count(),
        candidate_count: candidates.len(),
        conflict_count: 0,
        applied_count: 0,
    };
    io::write_json(&summary_json, &summary)?;

    let manifest = BridgeManifest {
        contigs_paths: contigs_paths.display().to_string(),
        scaffolds_paths: scaffolds_paths.map(|path| path.display().to_string()),
        source_graph_gfa: graph_gfa.display().to_string(),
        bridged_graph_gfa: bridged_graph_gfa.display().to_string(),
        candidate_tsv: candidates_tsv.display().to_string(),
        evidence_jsonl: candidates_jsonl.display().to_string(),
        conflict_tsv: conflicts_tsv.display().to_string(),
        decision_jsonl: decisions_jsonl.display().to_string(),
    };
    io::write_json(&manifest_json, &manifest)?;

    io::copy_file(&candidates_tsv, &legacy_candidates_tsv)?;
    io::copy_file(&candidates_jsonl, &legacy_candidates_jsonl)?;

    Ok(BridgeArtifacts {
        summary,
        manifest_json,
        candidates_tsv,
        candidates_jsonl,
        conflicts_tsv,
        decisions_jsonl,
        bridged_graph_gfa,
        summary_json,
    })
}

fn generate_candidates(records: &[PathRecord]) -> Vec<BridgeCandidate> {
    let mut candidates = Vec::new();
    for record in records {
        for chunk in &record.chunks {
            for pair in chunk.steps.windows(2) {
                candidates.push(BridgeCandidate {
                    bridge_id: format!("bridge_{:05}", candidates.len() + 1),
                    source: record.name.clone(),
                    left_segment: pair[0].segment_id.clone(),
                    right_segment: pair[1].segment_id.clone(),
                    path_segments: vec![
                        format!(
                            "{}{}",
                            pair[0].segment_id,
                            if pair[0].forward { "+" } else { "-" }
                        ),
                        format!(
                            "{}{}",
                            pair[1].segment_id,
                            if pair[1].forward { "+" } else { "-" }
                        ),
                    ],
                    bridge_length: 0,
                    depth_agreement: 1.0,
                    path_self_contained: true,
                    insert_size_penalty: 0.0,
                    quality: 1.0,
                    status: "recorded".to_string(),
                    reason: "adjacent_path_steps".to_string(),
                });
            }
        }
    }
    candidates
}

fn write_candidates_tsv(candidates: &[BridgeCandidate], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "bridge_id\tsource\tleft_segment\tright_segment\tpath_segments\tbridge_length\tdepth_agreement\tpath_self_contained\tinsert_size_penalty\tquality\tstatus\treason"
    )?;
    for candidate in candidates {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            candidate.bridge_id,
            candidate.source,
            candidate.left_segment,
            candidate.right_segment,
            candidate.path_segments.join(","),
            candidate.bridge_length,
            candidate.depth_agreement,
            candidate.path_self_contained,
            candidate.insert_size_penalty,
            candidate.quality,
            candidate.status,
            candidate.reason
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_candidates_jsonl(candidates: &[BridgeCandidate], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    for candidate in candidates {
        serde_json::to_writer(&mut writer, candidate)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(())
}

fn write_conflicts_tsv(path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "bridge_id\tconflict_type\tdetail")?;
    writer.flush()?;
    Ok(())
}

fn write_decisions_jsonl(candidates: &[BridgeCandidate], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    for candidate in candidates {
        let decision = BridgeDecision {
            bridge_id: candidate.bridge_id.clone(),
            status: "not_applied".to_string(),
            reason: "current_stage_records_evidence_only".to_string(),
            applied: false,
        };
        serde_json::to_writer(&mut writer, &decision)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn stages_bridging_from_fixture_paths_golden() {
        let temp = tempdir().unwrap();
        let graph_gfa = temp.path().join("graph.gfa");
        let contigs_paths = temp.path().join("contigs.paths");
        fs::write(
            &graph_gfa,
            "H\tVN:Z:1.2\nS\t1\tACGT\tLN:i:4\nS\t2\tTGCA\tLN:i:4\n",
        )
        .unwrap();
        fs::write(
            &contigs_paths,
            "NODE_1_length_34015_cov_41.9818\n115+,123+,143+,205+,202+,304+,87+,278+,125+,88+,129-,131+,92+,190-,189+\n",
        )
        .unwrap();

        let outdir = temp.path().join("bridging");
        let artifacts = stage_bridging(&graph_gfa, &contigs_paths, None, &outdir).unwrap();
        assert_eq!(artifacts.summary.path_record_count, 1);
        assert_eq!(artifacts.summary.forward_path_records, 1);
        assert_eq!(artifacts.summary.reverse_path_records, 0);
        assert_eq!(artifacts.summary.candidate_count, 14);
        assert_eq!(artifacts.summary.conflict_count, 0);
        assert_eq!(artifacts.summary.applied_count, 0);

        let tsv = fs::read_to_string(outdir.join(contract::BRIDGING_CANDIDATES_TSV)).unwrap();
        assert!(tsv.starts_with(
            "bridge_id\tsource\tleft_segment\tright_segment\tpath_segments\tbridge_length\tdepth_agreement\tpath_self_contained\tinsert_size_penalty\tquality\tstatus\treason\n"
        ));
        assert!(tsv.contains(
            "bridge_00001\tNODE_1_length_34015_cov_41.9818\t115\t123\t115+,123+\t0\t1.0000\ttrue\t0.0000\t1.0000\trecorded\tadjacent_path_steps\n"
        ));
        assert!(outdir.join(contract::BRIDGING_CONFLICTS_TSV).exists());
        assert!(outdir.join(contract::BRIDGING_DECISIONS_JSONL).exists());
        assert!(outdir.join(contract::BRIDGING_GRAPH_GFA).exists());
        assert!(
            outdir
                .join(contract::BRIDGING_LEGACY_CANDIDATES_TSV)
                .exists()
        );
    }
}
