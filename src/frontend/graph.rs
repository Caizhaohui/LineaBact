use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::io;

use super::contract;
use super::report::GraphSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphArtifacts {
    pub summary: GraphSummary,
    pub normalized_gfa: PathBuf,
    pub fastg: PathBuf,
    pub contigs_paths: PathBuf,
    pub scaffolds_paths: Option<PathBuf>,
    pub segments_tsv: PathBuf,
    pub links_tsv: PathBuf,
    pub depth_tsv: PathBuf,
    pub paths_tsv: PathBuf,
    pub anchor_segments_tsv: PathBuf,
    pub graph_qc_tsv: PathBuf,
    pub paths_summary_tsv: PathBuf,
    pub summary_json: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphSegment {
    pub id: String,
    pub length: usize,
    pub kmer_count: Option<u64>,
    pub depth: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphLink {
    pub from_segment: String,
    pub from_forward: bool,
    pub to_segment: String,
    pub to_forward: bool,
    pub overlap: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathStep {
    pub segment_id: String,
    pub forward: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathChunk {
    pub steps: Vec<PathStep>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathRecord {
    pub name: String,
    pub base_name: String,
    pub is_reverse: bool,
    pub declared_length: Option<usize>,
    pub coverage: Option<f64>,
    pub chunks: Vec<PathChunk>,
}

pub fn stage_graph(spades_dir: &Path, outdir: &Path) -> Result<GraphArtifacts> {
    fs::create_dir_all(outdir)?;

    let normalized_gfa = contract::root_path(outdir, contract::GRAPH_GFA);
    let fastg = contract::root_path(outdir, contract::GRAPH_FASTG);
    let contigs_paths = contract::root_path(outdir, contract::GRAPH_CONTIGS_PATHS);
    let scaffolds_paths = contract::root_path(outdir, contract::GRAPH_SCAFFOLDS_PATHS);
    let segments_tsv = contract::root_path(outdir, contract::GRAPH_SEGMENTS_TSV);
    let links_tsv = contract::root_path(outdir, contract::GRAPH_LINKS_TSV);
    let depth_tsv = contract::root_path(outdir, contract::GRAPH_DEPTH_TSV);
    let paths_tsv = contract::root_path(outdir, contract::GRAPH_PATHS_TSV);
    let anchor_segments_tsv = contract::root_path(outdir, contract::GRAPH_ANCHOR_SEGMENTS_TSV);
    let graph_qc_tsv = contract::root_path(outdir, contract::GRAPH_QC_TSV);
    let paths_summary_tsv = contract::root_path(outdir, contract::GRAPH_PATHS_SUMMARY_TSV);
    let summary_json = contract::root_path(outdir, contract::GRAPH_SUMMARY_JSON);

    io::copy_file(&spades_dir.join(contract::SPADES_GFA), &normalized_gfa)?;
    io::copy_file(&spades_dir.join(contract::SPADES_FASTG), &fastg)?;
    io::copy_file(
        &spades_dir.join(contract::SPADES_CONTIGS_PATHS),
        &contigs_paths,
    )?;
    let scaffolds_paths = copy_optional_file(
        &spades_dir.join(contract::SPADES_SCAFFOLDS_PATHS),
        &scaffolds_paths,
    )?;

    let (segments, links, summary) = summarize_gfa(&normalized_gfa)?;
    write_segments_tsv(&segments, &segments_tsv)?;
    write_links_tsv(&links, &links_tsv)?;
    write_depth_tsv(&segments, &depth_tsv)?;

    let mut path_records = parse_paths_file(&contigs_paths)?;
    if let Some(scaffolds_paths) = scaffolds_paths.as_ref() {
        path_records.extend(parse_paths_file(scaffolds_paths)?);
    }
    write_paths_tsv(&path_records, &paths_tsv)?;
    write_paths_summary_tsv(&path_records, &paths_summary_tsv)?;
    write_anchor_segments_tsv(&path_records, &segments, &anchor_segments_tsv)?;
    write_graph_qc_tsv(&summary, &segments, &path_records, &graph_qc_tsv)?;

    io::write_json(&summary_json, &summary)?;

    Ok(GraphArtifacts {
        summary,
        normalized_gfa,
        fastg,
        contigs_paths,
        scaffolds_paths,
        segments_tsv,
        links_tsv,
        depth_tsv,
        paths_tsv,
        anchor_segments_tsv,
        graph_qc_tsv,
        paths_summary_tsv,
        summary_json,
    })
}

#[cfg(test)]
pub fn write_contigs_from_gfa(gfa_path: &Path, fasta_path: &Path) -> Result<()> {
    io::ensure_parent_dir(fasta_path)?;
    let text = fs::read_to_string(gfa_path)?;
    let file = fs::File::create(fasta_path)?;
    let mut writer = BufWriter::new(file);

    for line in text.lines() {
        if !line.starts_with("S\t") {
            continue;
        }
        let mut fields = line.split('\t');
        let _ = fields.next();
        let id = fields
            .next()
            .context("missing segment id in GFA segment line")?;
        let seq = fields
            .next()
            .context("missing segment sequence in GFA segment line")?;
        let sequence = if seq == "*" {
            "N".to_string()
        } else {
            seq.to_string()
        };

        writeln!(writer, ">{id}")?;
        for chunk in sequence.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }

    writer.flush()?;
    Ok(())
}

pub fn parse_paths_file(path: &Path) -> Result<Vec<PathRecord>> {
    let text = fs::read_to_string(path)?;
    let mut records = Vec::new();
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .collect();
    let mut index = 0usize;

    while index < lines.len() {
        let name_line = lines[index];
        if is_path_steps_line(name_line) {
            bail!(
                "paths file {} starts a record with a steps line: {name_line}",
                path.display()
            );
        }
        index += 1;

        let mut steps_lines = Vec::new();
        while index < lines.len() && is_path_steps_line(lines[index]) {
            steps_lines.push(lines[index]);
            index += 1;
        }
        if steps_lines.is_empty() {
            bail!(
                "paths file {} record {name_line} has no step lines",
                path.display()
            );
        }

        let steps_line = steps_lines.concat();
        records.push(parse_path_record(name_line, &steps_line)?);
    }

    Ok(records)
}

fn is_path_steps_line(line: &str) -> bool {
    let line = line.trim();
    if line.is_empty() {
        return false;
    }
    line.contains(',') || line.ends_with('+') || line.ends_with('-') || line.ends_with(';')
}

fn parse_path_record(name_line: &str, steps_line: &str) -> Result<PathRecord> {
    let is_reverse = name_line.ends_with('\'');
    let name = name_line.trim_end_matches('\'').to_string();
    let (base_name, declared_length, coverage) = parse_path_name(&name)?;
    let chunks = steps_line
        .split(';')
        .filter_map(|chunk| {
            let chunk = chunk.trim();
            if chunk.is_empty() { None } else { Some(chunk) }
        })
        .map(|chunk| {
            chunk
                .split(',')
                .filter(|token| !token.trim().is_empty())
                .map(parse_path_step)
                .collect::<Result<Vec<_>>>()
                .map(|steps| PathChunk { steps })
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(PathRecord {
        name,
        base_name,
        is_reverse,
        declared_length,
        coverage,
        chunks,
    })
}

impl PathRecord {
    pub fn all_steps(&self) -> Vec<&PathStep> {
        self.chunks
            .iter()
            .flat_map(|chunk| chunk.steps.iter())
            .collect()
    }
}

fn parse_path_name(name: &str) -> Result<(String, Option<usize>, Option<f64>)> {
    let Some((base_name, length_part)) = name.split_once("_length_") else {
        return Ok((name.to_string(), None, None));
    };
    let Some((length_text, cov_part)) = length_part.split_once("_cov_") else {
        return Ok((base_name.to_string(), None, None));
    };
    let declared_length = length_text
        .parse::<usize>()
        .with_context(|| format!("failed to parse declared path length from {name}"))?;
    let coverage = cov_part.parse::<f64>().ok();
    Ok((base_name.to_string(), Some(declared_length), coverage))
}

fn parse_path_step(token: &str) -> Result<PathStep> {
    let token = token.trim();
    if token.len() < 2 {
        bail!("invalid path step token {token:?}");
    }
    let (segment_id, orientation) = token.split_at(token.len() - 1);
    let forward = match orientation {
        "+" => true,
        "-" => false,
        _ => bail!("invalid path step orientation in token {token:?}"),
    };
    Ok(PathStep {
        segment_id: segment_id.to_string(),
        forward,
    })
}

fn summarize_gfa(path: &Path) -> Result<(Vec<GraphSegment>, Vec<GraphLink>, GraphSummary)> {
    let text = fs::read_to_string(path)?;
    let mut segments = Vec::new();
    let mut links = Vec::new();
    let mut link_count = 0usize;
    let mut path_record_count = 0usize;

    for line in text.lines() {
        if line.starts_with("S\t") {
            let segment = parse_segment_line(line)?;
            segments.push(segment);
        } else if line.starts_with("L\t") {
            link_count += 1;
            links.push(parse_link_line(line)?);
        } else if line.starts_with("P\t") {
            path_record_count += 1;
        }
    }

    let lengths: Vec<usize> = segments.iter().map(|segment| segment.length).collect();
    let total_segment_bases: usize = lengths.iter().sum();
    let mean_segment_bases = if lengths.is_empty() {
        0.0
    } else {
        total_segment_bases as f64 / lengths.len() as f64
    };
    let max_segment_bases = lengths.iter().copied().max().unwrap_or(0);
    let min_segment_bases = lengths.iter().copied().min().unwrap_or(0);
    let depth_values: Vec<f64> = segments
        .iter()
        .filter_map(|segment| segment.depth)
        .collect();
    let mean_depth = if depth_values.is_empty() {
        0.0
    } else {
        depth_values.iter().sum::<f64>() / depth_values.len() as f64
    };

    Ok((
        segments,
        links,
        GraphSummary {
            segment_count: lengths.len(),
            link_count,
            path_record_count,
            total_segment_bases,
            mean_segment_bases,
            max_segment_bases,
            min_segment_bases,
            mean_depth,
        },
    ))
}

fn parse_link_line(line: &str) -> Result<GraphLink> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        bail!("invalid GFA link line: {line}");
    }
    Ok(GraphLink {
        from_segment: fields[1].to_string(),
        from_forward: parse_orientation(fields[2])?,
        to_segment: fields[3].to_string(),
        to_forward: parse_orientation(fields[4])?,
        overlap: fields[5].to_string(),
    })
}

fn parse_orientation(token: &str) -> Result<bool> {
    match token {
        "+" => Ok(true),
        "-" => Ok(false),
        _ => bail!("invalid orientation token: {token}"),
    }
}

fn parse_segment_line(line: &str) -> Result<GraphSegment> {
    let mut fields = line.split('\t');
    let _ = fields.next();
    let id = fields
        .next()
        .context("missing segment id in GFA segment line")?
        .to_string();
    let seq = fields
        .next()
        .context("missing segment sequence in GFA segment line")?;
    let sequence_length = if seq == "*" { 0 } else { seq.len() };

    let mut length = None;
    let mut kmer_count = None;
    let mut depth = None;
    for field in fields {
        if let Some(value) = field.strip_prefix("LN:i:") {
            length = value.parse::<usize>().ok();
        } else if let Some(value) = field.strip_prefix("KC:i:") {
            kmer_count = value.parse::<u64>().ok();
        } else if let Some(value) = field.strip_prefix("dp:f:") {
            depth = value.parse::<f64>().ok();
        }
    }

    Ok(GraphSegment {
        id,
        length: length.unwrap_or(sequence_length),
        kmer_count,
        depth,
    })
}

fn write_segments_tsv(segments: &[GraphSegment], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "segment_id\tlength\tkmer_count\tdepth")?;
    for segment in segments {
        let kmer_count = segment
            .kmer_count
            .map(|value| value.to_string())
            .unwrap_or_else(|| String::from(""));
        let depth = segment
            .depth
            .map(|value| format!("{value:.4}"))
            .unwrap_or_else(|| String::from(""));
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            segment.id, segment.length, kmer_count, depth
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_links_tsv(links: &[GraphLink], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "from_segment\tfrom_orientation\tto_segment\tto_orientation\toverlap"
    )?;
    for link in links {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            link.from_segment,
            if link.from_forward { "+" } else { "-" },
            link.to_segment,
            if link.to_forward { "+" } else { "-" },
            link.overlap
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_depth_tsv(segments: &[GraphSegment], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "segment_id\tlength\tdepth")?;
    for segment in segments {
        let depth = segment
            .depth
            .map(|value| format!("{value:.4}"))
            .unwrap_or_default();
        writeln!(writer, "{}\t{}\t{}", segment.id, segment.length, depth)?;
    }
    writer.flush()?;
    Ok(())
}

fn copy_optional_file(source: &Path, destination: &Path) -> Result<Option<PathBuf>> {
    if source.exists() {
        io::copy_file(source, destination)?;
        Ok(Some(destination.to_path_buf()))
    } else {
        Ok(None)
    }
}

fn write_paths_tsv(records: &[PathRecord], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "name\tbase_name\torientation\tchunk_index\tstep_index\tsegment_id\tsegment_orientation"
    )?;
    for record in records {
        let orientation = if record.is_reverse {
            "reverse"
        } else {
            "forward"
        };
        for (chunk_index, chunk) in record.chunks.iter().enumerate() {
            for (step_index, step) in chunk.steps.iter().enumerate() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.name,
                    record.base_name,
                    orientation,
                    chunk_index + 1,
                    step_index + 1,
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                )?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_paths_summary_tsv(records: &[PathRecord], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "name\tbase_name\torientation\tdeclared_length\tcoverage\tstep_count\tgap_count\tfirst_step\tlast_step"
    )?;
    for record in records {
        let orientation = if record.is_reverse {
            "reverse"
        } else {
            "forward"
        };
        let declared_length = record
            .declared_length
            .map(|value| value.to_string())
            .unwrap_or_else(|| String::from(""));
        let coverage = record
            .coverage
            .map(|value| format!("{value:.4}"))
            .unwrap_or_else(|| String::from(""));
        let steps = record.all_steps();
        let first_step = steps
            .first()
            .map(|step| {
                format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                )
            })
            .unwrap_or_else(|| String::from(""));
        let last_step = steps
            .last()
            .map(|step| {
                format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                )
            })
            .unwrap_or_else(|| String::from(""));
        let gap_count = record.chunks.len().saturating_sub(1);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.name,
            record.base_name,
            orientation,
            declared_length,
            coverage,
            steps.len(),
            gap_count,
            first_step,
            last_step
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_anchor_segments_tsv(
    records: &[PathRecord],
    segments: &[GraphSegment],
    path: &Path,
) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "path_name\tpath_orientation\tanchor_kind\tsegment_id\tsegment_orientation\tsegment_length\tsegment_depth"
    )?;

    for record in records {
        let orientation = if record.is_reverse {
            "reverse"
        } else {
            "forward"
        };
        let steps = record.all_steps();
        for (anchor_kind, step) in [("first", steps.first()), ("last", steps.last())] {
            let Some(step) = step else {
                continue;
            };
            let segment = segments
                .iter()
                .find(|segment| segment.id == step.segment_id);
            let segment_length = segment
                .map(|value| value.length.to_string())
                .unwrap_or_default();
            let segment_depth = segment
                .and_then(|value| value.depth)
                .map(|value| format!("{value:.4}"))
                .unwrap_or_default();
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                record.name,
                orientation,
                anchor_kind,
                step.segment_id,
                if step.forward { "+" } else { "-" },
                segment_length,
                segment_depth
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_graph_qc_tsv(
    summary: &GraphSummary,
    segments: &[GraphSegment],
    records: &[PathRecord],
    path: &Path,
) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    let segments_with_depth = segments
        .iter()
        .filter(|segment| segment.depth.is_some())
        .count();
    let path_steps: usize = records.iter().map(|record| record.all_steps().len()).sum();
    let path_gaps: usize = records
        .iter()
        .map(|record| record.chunks.len().saturating_sub(1))
        .sum();

    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "segment_count\t{}", summary.segment_count)?;
    writeln!(writer, "link_count\t{}", summary.link_count)?;
    writeln!(
        writer,
        "gfa_path_record_count\t{}",
        summary.path_record_count
    )?;
    writeln!(writer, "paths_file_record_count\t{}", records.len())?;
    writeln!(
        writer,
        "total_segment_bases\t{}",
        summary.total_segment_bases
    )?;
    writeln!(
        writer,
        "mean_segment_bases\t{:.4}",
        summary.mean_segment_bases
    )?;
    writeln!(writer, "max_segment_bases\t{}", summary.max_segment_bases)?;
    writeln!(writer, "min_segment_bases\t{}", summary.min_segment_bases)?;
    writeln!(writer, "segments_with_depth\t{}", segments_with_depth)?;
    writeln!(writer, "mean_depth\t{:.4}", summary.mean_depth)?;
    writeln!(writer, "path_step_count\t{}", path_steps)?;
    writeln!(writer, "path_gap_count\t{}", path_gaps)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frontend::contract;
    use crate::frontend::fixture::extract_zip_entry;
    use tempfile::tempdir;

    #[test]
    fn parses_unicycler_paths_fixture() {
        let temp = tempdir().unwrap();
        let archive = Path::new("reference_tools/Unicycler-main.zip");
        let paths_path = temp.path().join("assembly_graph.fastg.paths");
        extract_zip_entry(
            archive,
            "Unicycler-main/test/test_assembly_graph.fastg.paths",
            &paths_path,
        )
        .unwrap();

        let records = parse_paths_file(&paths_path).unwrap();
        assert_eq!(records.len(), 280);

        let first = &records[0];
        assert_eq!(first.name, "NODE_1_length_34015_cov_41.9818");
        assert!(!first.is_reverse);
        assert_eq!(first.base_name, "NODE_1");
        assert_eq!(first.declared_length, Some(34015));
        assert!((first.coverage.unwrap() - 41.9818).abs() < 1e-6);
        assert_eq!(first.chunks.len(), 1);
        assert_eq!(first.all_steps().len(), 15);
        assert_eq!(
            first
                .all_steps()
                .iter()
                .map(|step| format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                ))
                .collect::<Vec<_>>(),
            vec![
                "115+", "123+", "143+", "205+", "202+", "304+", "87+", "278+", "125+", "88+",
                "129-", "131+", "92+", "190-", "189+",
            ]
        );

        let reverse = &records[1];
        assert_eq!(reverse.name, "NODE_1_length_34015_cov_41.9818");
        assert!(reverse.is_reverse);
        assert_eq!(reverse.base_name, "NODE_1");
        assert_eq!(reverse.chunks.len(), 1);
        assert_eq!(reverse.all_steps().len(), 15);
        assert_eq!(
            reverse
                .all_steps()
                .iter()
                .map(|step| format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                ))
                .collect::<Vec<_>>(),
            vec![
                "189-", "190+", "92-", "131-", "129+", "88-", "125-", "278-", "87-", "304-",
                "202-", "205-", "143-", "123-", "115-",
            ]
        );
    }

    #[test]
    fn parses_paths_chunks_with_gaps() {
        let temp = tempdir().unwrap();
        let paths = temp.path().join("gap.paths");
        fs::write(&paths, "NODE_gap\n1+,2+;3-,4-\n").unwrap();

        let records = parse_paths_file(&paths).unwrap();
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(record.chunks.len(), 2);
        assert_eq!(record.all_steps().len(), 4);
        assert_eq!(
            record
                .all_steps()
                .iter()
                .map(|step| format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                ))
                .collect::<Vec<_>>(),
            vec!["1+", "2+", "3-", "4-"]
        );
    }

    #[test]
    fn parses_multiline_spades_path_steps() {
        let temp = tempdir().unwrap();
        let paths = temp.path().join("multiline.paths");
        fs::write(
            &paths,
            concat!(
                "NODE_148_length_14519_cov_37.617326\n",
                "3466708+;\n",
                "3460484+\n",
                "NODE_148_length_14519_cov_37.617326'\n",
                "3460484-;\n",
                "3466708-\n"
            ),
        )
        .unwrap();

        let records = parse_paths_file(&paths).unwrap();
        assert_eq!(records.len(), 2);

        let forward = &records[0];
        assert_eq!(forward.name, "NODE_148_length_14519_cov_37.617326");
        assert_eq!(forward.chunks.len(), 2);
        assert_eq!(
            forward
                .all_steps()
                .iter()
                .map(|step| format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                ))
                .collect::<Vec<_>>(),
            vec!["3466708+", "3460484+"]
        );

        let reverse = &records[1];
        assert!(reverse.is_reverse);
        assert_eq!(reverse.chunks.len(), 2);
        assert_eq!(
            reverse
                .all_steps()
                .iter()
                .map(|step| format!(
                    "{}{}",
                    step.segment_id,
                    if step.forward { "+" } else { "-" }
                ))
                .collect::<Vec<_>>(),
            vec!["3460484-", "3466708-"]
        );
    }

    #[test]
    fn stages_graph_from_unicycler_fixture_golden() {
        let temp = tempdir().unwrap();
        let spades_dir = temp.path().join("spades");
        let outdir = temp.path().join("graph");
        fs::create_dir_all(&spades_dir).unwrap();

        let archive = Path::new("reference_tools/Unicycler-main.zip");
        extract_zip_entry(
            archive,
            "Unicycler-main/test/test_assembly_graph.gfa",
            &spades_dir.join(contract::SPADES_GFA),
        )
        .unwrap();
        extract_zip_entry(
            archive,
            "Unicycler-main/test/test_assembly_graph.fastg",
            &spades_dir.join(contract::SPADES_FASTG),
        )
        .unwrap();
        extract_zip_entry(
            archive,
            "Unicycler-main/test/test_assembly_graph.fastg.paths",
            &spades_dir.join(contract::SPADES_CONTIGS_PATHS),
        )
        .unwrap();

        let artifacts = stage_graph(&spades_dir, &outdir).unwrap();
        assert_eq!(artifacts.summary.segment_count, 336);
        assert_eq!(artifacts.summary.link_count, 452);
        assert_eq!(artifacts.summary.path_record_count, 140);
        assert_eq!(artifacts.summary.total_segment_bases, 187_896);
        assert_eq!(artifacts.summary.max_segment_bases, 32_060);
        assert_eq!(artifacts.summary.min_segment_bases, 26);

        let segments = fs::read_to_string(outdir.join(contract::GRAPH_SEGMENTS_TSV)).unwrap();
        assert!(segments.starts_with("segment_id\tlength\tkmer_count\tdepth\n"));
        assert!(segments.contains("1\t449\t37087\t82.5991\n"));

        let paths = fs::read_to_string(outdir.join(contract::GRAPH_PATHS_SUMMARY_TSV)).unwrap();
        assert!(paths.starts_with(
            "name\tbase_name\torientation\tdeclared_length\tcoverage\tstep_count\tgap_count\tfirst_step\tlast_step\n"
        ));
        assert!(paths.contains(
            "NODE_1_length_34015_cov_41.9818\tNODE_1\tforward\t34015\t41.9818\t15\t0\t115+\t189+\n"
        ));
        assert!(paths.contains(
            "NODE_1_length_34015_cov_41.9818\tNODE_1\treverse\t34015\t41.9818\t15\t0\t189-\t115-\n"
        ));
    }

    #[test]
    fn writes_contigs_from_gfa_segments() {
        let temp = tempdir().unwrap();
        let gfa = temp.path().join("graph.gfa");
        fs::write(
            &gfa,
            "H\tVN:Z:1.0\nS\tseg1\tACGT\tLN:i:4\tKC:i:2\tdp:f:1.0\nS\tseg2\tGGCC\tLN:i:4\tKC:i:2\tdp:f:2.0\n",
        )
        .unwrap();
        let fasta = temp.path().join("contigs.fasta");
        write_contigs_from_gfa(&gfa, &fasta).unwrap();

        let text = fs::read_to_string(&fasta).unwrap();
        assert!(text.contains(">seg1"));
        assert!(text.contains("ACGT"));
        assert!(text.contains(">seg2"));
    }
}
