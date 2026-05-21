#!/usr/bin/env bash
set -euo pipefail

manifest="${1:-rawdata/final_9_species_hybrid_runs.tsv}"
out_root="${2:-rawdata/fastq}"
tmp_root="${3:-rawdata/fastq_tmp}"

tool="/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/fasterq-dump-orig.3.4.1"
fallback_tool="/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/fastq-dump-orig.3.4.1"

mkdir -p "$out_root" "$tmp_root"

awk -F '\t' 'NR > 1 { print $3 }' "$manifest" | while IFS= read -r run; do
  [ -n "$run" ] || continue

  outdir="$out_root/$run"
  should_skip=false
  for candidate in "$outdir/$run.fastq.gz" "$outdir/${run}_1.fastq.gz" "$outdir/${run}_2.fastq.gz"; do
    if [ -f "$candidate" ] && gzip -t "$candidate" >/dev/null 2>&1; then
      should_skip=true
      break
    fi
  done
  if [ "$should_skip" = true ]; then
    printf '[%s] skip %s\n' "$(date +%F\ %T)" "$run"
    continue
  fi

  mkdir -p "$outdir" "$tmp_root/$run"
  printf '[%s] start %s\n' "$(date +%F\ %T)" "$run"

  if ! "$tool" -f --split-files -e 8 -t "$tmp_root/$run" -O "$outdir" "rawdata/sra/$run/$run.sra"; then
    printf '[%s] retry %s with fastq-dump\n' "$(date +%F\ %T)" "$run"
    "$fallback_tool" --split-files --skip-technical --gzip -O "$outdir" "rawdata/sra/$run/$run.sra"
  fi

  shopt -s nullglob
  for fastq in "$outdir"/*.fastq; do
    pigz -p 4 -1 -f "$fastq"
  done
  shopt -u nullglob

  rmdir "$tmp_root/$run" 2>/dev/null || true
  printf '[%s] done %s\n' "$(date +%F\ %T)" "$run"
done
