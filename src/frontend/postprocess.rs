use anyhow::Result;
use needletail::parse_fastx_file;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::io;

use super::contract;

const MIN_CONTIG_LENGTH_BP: usize = 200;

#[derive(Debug, Clone)]
pub struct PostprocessArtifacts {
    pub backbone_fasta: PathBuf,
    pub finished_fasta: PathBuf,
    pub contig_stats_tsv: PathBuf,
    pub rename_map_tsv: PathBuf,
}

pub fn stage_postprocess(
    contigs_fasta: &Path,
    scaffolds_fasta: &Path,
    outdir: &Path,
) -> Result<PostprocessArtifacts> {
    fs::create_dir_all(outdir)?;
    let backbone_fasta = contract::root_path(outdir, contract::POSTPROCESS_BACKBONE_FASTA);
    let finished_fasta = contract::root_path(outdir, contract::POSTPROCESS_FINISHED_FASTA);
    let contig_stats_tsv = contract::root_path(outdir, contract::POSTPROCESS_CONTIG_STATS_TSV);
    let rename_map_tsv = contract::root_path(outdir, contract::POSTPROCESS_RENAME_MAP_TSV);

    let mut stats = Vec::new();
    let mut rename_map = Vec::new();
    rewrite_fasta(
        contigs_fasta,
        &backbone_fasta,
        "backbone",
        &mut stats,
        &mut rename_map,
    )?;
    rewrite_fasta(
        scaffolds_fasta,
        &finished_fasta,
        "finished",
        &mut stats,
        &mut rename_map,
    )?;
    write_contig_stats_tsv(&stats, &contig_stats_tsv)?;
    write_rename_map_tsv(&rename_map, &rename_map_tsv)?;

    Ok(PostprocessArtifacts {
        backbone_fasta,
        finished_fasta,
        contig_stats_tsv,
        rename_map_tsv,
    })
}

#[derive(Debug)]
struct ContigStatRow {
    stage: String,
    source_name: String,
    renamed_name: String,
    length_bp: usize,
    gc_content: f64,
    retained: bool,
    reason: String,
}

#[derive(Debug)]
struct RenameMapRow {
    stage: String,
    source_name: String,
    renamed_name: String,
}

fn rewrite_fasta(
    input_path: &Path,
    output_path: &Path,
    stage: &str,
    stats: &mut Vec<ContigStatRow>,
    rename_map: &mut Vec<RenameMapRow>,
) -> Result<()> {
    io::ensure_parent_dir(output_path)?;
    let mut reader = parse_fastx_file(input_path)?;
    let file = fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    let mut retained_index = 0usize;

    while let Some(record) = reader.next() {
        let record = record?;
        let source_name = String::from_utf8_lossy(record.id()).to_string();
        let seq_storage = record.seq();
        let seq = seq_storage.as_ref();
        let length_bp = seq.len();
        let gc_content = compute_gc_content(seq);
        let retained = length_bp >= MIN_CONTIG_LENGTH_BP;
        let reason = if retained {
            "retained".to_string()
        } else {
            format!("dropped_below_{}bp", MIN_CONTIG_LENGTH_BP)
        };
        let renamed_name = if retained {
            retained_index += 1;
            format!(
                "{stage}_contig_{:06} len={} gc={:.4} source={}",
                retained_index, length_bp, gc_content, source_name
            )
        } else {
            String::new()
        };

        stats.push(ContigStatRow {
            stage: stage.to_string(),
            source_name: source_name.clone(),
            renamed_name: renamed_name.clone(),
            length_bp,
            gc_content,
            retained,
            reason,
        });
        if retained {
            rename_map.push(RenameMapRow {
                stage: stage.to_string(),
                source_name,
                renamed_name: renamed_name.clone(),
            });
            write_fasta_record(&mut writer, &renamed_name, seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_fasta_record<W: Write>(writer: &mut W, header: &str, seq: &[u8]) -> Result<()> {
    writeln!(writer, ">{header}")?;
    for chunk in seq.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn compute_gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&base| matches!(base, b'G' | b'g' | b'C' | b'c'))
        .count();
    gc as f64 / seq.len() as f64
}

fn write_contig_stats_tsv(rows: &[ContigStatRow], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "stage\tsource_name\trenamed_name\tlength_bp\tgc_content\tretained\treason"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
            row.stage,
            row.source_name,
            row.renamed_name,
            row.length_bp,
            row.gc_content,
            row.retained,
            row.reason
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_rename_map_tsv(rows: &[RenameMapRow], path: &Path) -> Result<()> {
    io::ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "stage\tsource_name\trenamed_name")?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}",
            row.stage, row.source_name, row.renamed_name
        )?;
    }
    writer.flush()?;
    Ok(())
}
