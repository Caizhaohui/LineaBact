#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUTDIR="benchmarks/generated/reference_validation_runtime_manual"
RESULTS_TSV="$OUTDIR/runtime_runs.tsv"
ENV_PREFIX="/hpcfs/fhome/caizhh/.conda/envs/LineaBact"

LINEABACT="$ROOT_DIR/target/release/lineabact"
SPADES="$ENV_PREFIX/bin/spades.py"
SHOVILL="$ENV_PREFIX/bin/shovill"
SEQTK="$ENV_PREFIX/bin/seqtk"
export CONDA_PREFIX="$ENV_PREFIX"
export PATH="$ENV_PREFIX/bin:$PATH"

mkdir -p "$OUTDIR"
printf "sample\ttool\trepeat\tstatus\telapsed_seconds\toutdir\n" > "$RESULTS_TSV"

run_case() {
  local sample="$1"
  local tool="$2"
  local repeat="$3"
  local case_outdir="$4"
  shift 4
  local -a cmd=( "$@" )

  mkdir -p "$case_outdir"
  local start end elapsed status
  start="$(date +%s.%N)"
  status="ok"
  if ! "${cmd[@]}" >"$case_outdir/benchmark.stdout.log" 2>"$case_outdir/benchmark.stderr.log"; then
    status="failed"
  fi
  end="$(date +%s.%N)"
  elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN{printf "%.4f", e-s}')"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "$tool" "$repeat" "$status" "$elapsed" "$case_outdir" >> "$RESULTS_TSV"
  if [[ "$status" != "ok" ]]; then
    return 1
  fi
}

for repeat in 1 2 3; do
  repeat_tag="$(printf 'run_%03d' "$repeat")"

  run_case \
    "streptomyces_clavuligerus" \
    "lineabact" \
    "$repeat" \
    "$OUTDIR/streptomyces_clavuligerus/lineabact/repeats/$repeat_tag" \
    "$LINEABACT" assemble \
    --backend spades \
    --r1 rawdata/fastq/SRR16805542/SRR16805542_1.fastq.gz \
    --r2 rawdata/fastq/SRR16805542/SRR16805542_2.fastq.gz \
    --outdir "$OUTDIR/streptomyces_clavuligerus/lineabact/repeats/$repeat_tag" \
    --threads 4 \
    --k 55 \
    --spades-executable "$SPADES" \
    --trim-adapters \
    --trim-tool "$SEQTK" \
    --target-coverage 150 \
    --genome-size-bp 9161304

  run_case \
    "streptomyces_clavuligerus" \
    "shovill" \
    "$repeat" \
    "$OUTDIR/streptomyces_clavuligerus/shovill/repeats/$repeat_tag" \
    "$SHOVILL" \
    --R1 rawdata/fastq/SRR16805542/SRR16805542_1.fastq.gz \
    --R2 rawdata/fastq/SRR16805542/SRR16805542_2.fastq.gz \
    --outdir "$OUTDIR/streptomyces_clavuligerus/shovill/repeats/$repeat_tag" \
    --cpus 4 \
    --depth 150 \
    --assembler spades \
    --force

  run_case \
    "streptomyces_rimosus" \
    "lineabact" \
    "$repeat" \
    "$OUTDIR/streptomyces_rimosus/lineabact/repeats/$repeat_tag" \
    "$LINEABACT" assemble \
    --backend spades \
    --r1 rawdata/fastq/SRR24413071/SRR24413071_1.fastq.gz \
    --r2 rawdata/fastq/SRR24413071/SRR24413071_2.fastq.gz \
    --outdir "$OUTDIR/streptomyces_rimosus/lineabact/repeats/$repeat_tag" \
    --threads 4 \
    --k 55 \
    --spades-executable "$SPADES" \
    --trim-adapters \
    --trim-tool "$SEQTK" \
    --target-coverage 150 \
    --genome-size-bp 9643891

  run_case \
    "streptomyces_rimosus" \
    "shovill" \
    "$repeat" \
    "$OUTDIR/streptomyces_rimosus/shovill/repeats/$repeat_tag" \
    "$SHOVILL" \
    --R1 rawdata/fastq/SRR24413071/SRR24413071_1.fastq.gz \
    --R2 rawdata/fastq/SRR24413071/SRR24413071_2.fastq.gz \
    --outdir "$OUTDIR/streptomyces_rimosus/shovill/repeats/$repeat_tag" \
    --cpus 4 \
    --depth 150 \
    --assembler spades \
    --force
done
