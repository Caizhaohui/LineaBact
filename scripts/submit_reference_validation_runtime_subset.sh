#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

ENV_PREFIX="/hpcfs/fhome/caizhh/.conda/envs/LineaBact"
PLAN_OUTDIR="benchmarks/generated/reference_validation_runtime_slurm"
PLAN_JSON="$PLAN_OUTDIR/benchmark_plan.json"
PLAN_FILTERED_JSON="$PLAN_OUTDIR/benchmark_plan.lineabact_shovill.json"
REPORT_OUTDIR="$PLAN_OUTDIR/report_lineabact_shovill"
MANIFEST="benchmarks/generated/reference_validation_runtime_subset.tsv"
CANDIDATE_PARTITIONS=(qcpu_23if qcpu_23i qcpu_23a qcpu_18i)

export CONDA_PREFIX="$ENV_PREFIX"
export PATH="$ENV_PREFIX/bin:$PATH"

mkdir -p "$PLAN_OUTDIR" "$REPORT_OUTDIR"

SELECTED_PARTITION="$("$ROOT_DIR/scripts/select_slurm_partition.sh" "${CANDIDATE_PARTITIONS[@]}")"
echo "Selected Slurm partition: $SELECTED_PARTITION"

cargo build --release

target/release/lineabact stats benchmark-plan \
  --manifest "$MANIFEST" \
  --sample-limit 2 \
  --fixture-root reference_tools/Unicycler-main \
  --outdir "$PLAN_OUTDIR" \
  --threads 4 \
  --k 55 \
  --lineabact-executable "$ROOT_DIR/target/release/lineabact" \
  --spades-executable "$ENV_PREFIX/bin/spades.py" \
  --shovill-executable "$ENV_PREFIX/bin/shovill" \
  --unicycler-executable "$ENV_PREFIX/bin/unicycler"

python - <<'PY'
import json
from pathlib import Path

plan_path = Path("benchmarks/generated/reference_validation_runtime_slurm/benchmark_plan.json")
filtered_path = Path("benchmarks/generated/reference_validation_runtime_slurm/benchmark_plan.lineabact_shovill.json")
plan = json.loads(plan_path.read_text())
plan["cases"] = [case for case in plan["cases"] if case["tool"] in {"lineabact", "shovill"}]
filtered_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n")
PY

target/release/lineabact stats benchmark-run \
  --plan "$PLAN_FILTERED_JSON" \
  --outdir "$REPORT_OUTDIR" \
  --repeat-count 3 \
  --scheduler slurm \
  --slurm-partition "$SELECTED_PARTITION" \
  --slurm-conda-base "$ENV_PREFIX" \
  --slurm-conda-env LineaBact \
  --slurm-cpus-per-task 4 \
  --slurm-mem-gb 64 \
  --slurm-time 12:00:00
