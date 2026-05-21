#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "usage: $0 <partition> [<partition> ...]" >&2
  exit 2
fi

best_partition=""
best_idle=-1
best_mixed=-1
query_timeout_seconds=20

for partition in "$@"; do
  if command -v timeout >/dev/null 2>&1; then
    query_cmd=(timeout "${query_timeout_seconds}s" sinfo -h -p "$partition" -o '%a|%D|%t')
  else
    query_cmd=(sinfo -h -p "$partition" -o '%a|%D|%t')
  fi
  if ! output="$("${query_cmd[@]}" 2>/dev/null)"; then
    continue
  fi
  idle_nodes=0
  mixed_nodes=0
  while IFS='|' read -r avail nodes state; do
    [ -n "$avail" ] || continue
    if [ "$avail" != "up" ]; then
      continue
    fi
    nodes="${nodes:-0}"
    case "${state,,}" in
      idle*)
        idle_nodes=$((idle_nodes + nodes))
        ;;
      mix*)
        mixed_nodes=$((mixed_nodes + nodes))
        ;;
    esac
  done <<< "$output"

  if [ "$idle_nodes" -gt "$best_idle" ] || { [ "$idle_nodes" -eq "$best_idle" ] && [ "$mixed_nodes" -gt "$best_mixed" ]; }; then
    best_partition="$partition"
    best_idle="$idle_nodes"
    best_mixed="$mixed_nodes"
  fi
done

if [ -z "$best_partition" ]; then
  echo "no usable Slurm partition found among: $*" >&2
  exit 1
fi

printf '%s\n' "$best_partition"
