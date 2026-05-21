# Testing & Benchmark Development Playbook

## Purpose
Establish rigorous testing and benchmarking standards for LineaBact, with special emphasis on terminal rescue correctness and reproducibility.

## Testing Hierarchy (Must Cover)

### 1. Unit Tests
- K-mer encoding, canonicalization, and reverse-complement
- DBG construction and compaction on synthetic data
- Terminal read recruitment logic
- Consensus building
- Attach / Reject decision logic
- Support-based trimming
- Left vs Right end independent processing

### 2. Integration Tests
- End-to-end backbone assembly on small datasets
- Terminal rescue on **synthetic terminal set**
- Full pipeline on `scl_L1` 50x / 100x downsampled data

### 3. Synthetic Terminal Set Testing (Highest Priority)
- Must create or use a dedicated synthetic dataset to evaluate terminal rescue.
- Test cases should include:
  - Successful terminal extension (true positive)
  - Low-support cases (should report `low_support` or `unresolved`)
  - Repeat-induced false extension (should be rejected)
  - Trimming behavior after extension
- Target: terminal extension base-level identity ≥ 99% on positive cases.

### 4. Negative Control Testing
- Include at least one **circular genome** dataset as negative control.
- LineaBact should **not** produce high-confidence terminal extensions on circular genomes.
- This is critical to validate the conservative nature of terminal rescue.

### 5. Regression & Golden Tests
- Maintain golden output files for key outputs (`assembly_stats.tsv`, `terminal_events.jsonl`, etc.).
- Any change that alters output schema or content must update golden files intentionally.

## Benchmark Requirements

### Reference Tools
- SPAdes (backbone assembly quality reference)
- Shovill (pipeline behavior reference)
- sparrowhawk (Rust assembler comparison)
- SKESA (conservative short-read assembler)
- Telomore (terminal finishing behavior reference)

### Key Evaluation Metrics
- Wall-clock runtime & peak memory
- `lineabact_vs_shovill` runtime ratio: `LineaBact elapsed_seconds / Shovill elapsed_seconds`, with `<= 1.0` treated as reaching Shovill speed level
- N50, contig number, total assembly length, genome fraction
- Misassemblies, mismatch rate, indel rate, and read mapping rate where reference exists
- **Terminal-specific metrics**:
  - Terminal bases recovered
  - Left/Right arm completeness
  - Candidate telomere recovery rate
  - Extension false-positive rate
  - Unresolved end count
  - Trimming ratio

### Benchmark Best Practices
- Always use **fixed random seed** for downsampling.
- Record all parameters in `params.toml` or `run_summary.json`.
- Run benchmark on multiple coverage levels: 50x / 100x / 200x / full.
- Separate **backbone assembly** evaluation from **terminal rescue** evaluation.
- Compare terminal recovery capability, not just contiguity (N50).
- For `lineabact_vs_shovill` runtime acceptance, do not use a single run. Execute at least 3 repeats per sample and evaluate the median runtime ratio.
- Treat backbone quality acceptance against `shovill` as a hard gate:
  - `lineabact` N50 must be at least 90% of `shovill`, with >= 95% as the preferred target.
  - `lineabact` genome fraction must not be lower than `shovill`.
  - `lineabact` misassembly count must not exceed `shovill`.
  - `lineabact` mismatch and indel rates must not be worse than `shovill`.
  - A faster sample is still rejected if misassemblies materially increase.
- Treat `lineabact stats benchmark-run --repeat-count 3` as the default recommended invocation; the CLI default should remain aligned with this rule.
- For the current short-read front-end phase, submit runtime benchmarks from the `LineaBact` conda environment and use the `qcpu_23if` Slurm partition. The maintained entry point is `scripts/submit_reference_validation_runtime_subset.sh`.

## Reproducibility Requirements

- All benchmark commands must be reproducible.
- Use fixed seeds for any random processes (downsampling, etc.).
- Output `run_summary.json` containing:
  - Git commit hash
  - Random seed used
  - Software versions
  - Resource usage (runtime, memory)
  - Key parameters

## Output Schema Testing

- Validate that all output files (`terminal_events.jsonl`, `terminal_qc.tsv`, `assembly_stats.tsv`, etc.) conform to expected schema.
- Changes to output format must be intentional and documented.

## Recommended Test Data Strategy

| Dataset Type           | Purpose                          | Coverage Levels     | Priority |
|------------------------|----------------------------------|---------------------|----------|
| Synthetic terminal set | Terminal rescue correctness      | 30x / 60x / 100x    | Highest  |
| `scl_L1` (downsampled) | Main development & regression    | 50x / 100x / 200x   | High     |
| `sclw_L1`              | High coverage & performance test | 50x / 100x / full   | Medium   |
| Circular negative control | Validate conservative behavior | Full                | High     |

## Common Pitfalls to Avoid

- Only optimizing for N50 while ignoring terminal recovery quality.
- Accepting a faster run even though misassemblies, mismatch rate, or indel rate got worse.
- Lack of negative control testing (leading to over-extension).
- Non-reproducible benchmark results.
- Insufficient testing of edge cases in terminal rescue (low support, repeats, trimming).
