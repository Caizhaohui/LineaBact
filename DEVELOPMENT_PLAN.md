# LineaBact 开发计划

_面向链霉菌等线性细菌基因组的 Rust 版 short-read-first assembly front-end：SPAdes-backed、Shovill-style reads 优化、Unicycler-like graph/bridging 输出，并逐步扩展 terminal-aware finishing。_

---

## 1. 项目定位

LineaBact 的目标是开发一个针对细菌线性基因组的 Rust 版 short-read-first assembly front-end。近期主线不是重写 assembler 内核，而是围绕 Illumina paired-end reads 建立可复现的前处理、SPAdes backend 编排、Shovill-style 优化、Unicycler-like graph/bridging 标准输出，并为后续线性 terminal finishing 准备稳定输入。

核心策略：

- **核心组装算法参考 SPAdes**：以 SPAdes 的 Illumina assembly、k-mer/DBG、graph simplification、paired-read evidence 和 assembly graph 输出作为 backbone assembly 的算法基线。[^spades]
- **流程优化参考 Shovill**：借鉴 Shovill 对 bacterial isolate paired-end reads 的输入检查、资源控制、adapter trimming、depth reduction、read correction、k-mer 选择、contig filtering/renaming 和可复现日志管理。[^shovill]
- **graph/bridging 输出参考 Unicycler**：借鉴 Unicycler short-read-first 思路，将 SPAdes graph、paths、depth、anchor segment 和 bridge candidate/evidence 标准化，供后续 long-read bridging、terminal bridging 和人工 Bandage 检查使用。[^unicycler]
- **Rust 工程实现参考 Sparrowhawk**：参考 Sparrowhawk 的 Rust crate 组织、bacterial short-read 约束、rayon 并行和图结构工程经验，但不把自研组装器作为 v1 产品主线。[^sparrowhawk]
- **线性基因组优化参考 Telomore**：借鉴 Telomore 对 actinomycete linear replicon 的用户指定线性 contig、terminal read recruitment、extension、QC map 和 finalized assembly 思路。[^telomore]
- **保守组装对照参考 SKESA/skesa-rs**：将 SKESA/skesa-rs 作为 microbial short-read conservative assembly 对照和 Rust library API 参考，而不是替代 SPAdes backbone。[^skesa]

当前阶段成功标准不是“重新实现 SPAdes”，也不是立刻解决线性基因组 terminal rescue，而是先做成稳定的 Rust front-end：输入 reads 优化、SPAdes 后端运行、assembly graph/paths/depth 标准化、bridging-ready 中间文件和可复现报告。

当前阶段范围：

- 开发测试数据只使用 `reference_tools/Unicycler-main.zip` 中的样例 reads、reference 和 graph fixtures。
- 只处理 short-read-first assembly front-end，不做线性基因组识别、terminal candidate detection、telomere/TIR、terminal rescue 或 Telomore-style finishing。
- 线性基因组相关设计保留为后续阶段，但不得阻塞当前 front-end 里程碑。

## 2. 参考工具分工

| 工具 | 在 LineaBact 中的角色 | 采用内容 | 明确不做 |
|---|---|---|---|
| SPAdes | 主 backbone assembler 和算法基线 | isolate Illumina assembly、assembly graph、multi-k/DBG 思路、graph simplification、hybrid evidence 参考 | v1 不重写 SPAdes 内核 |
| Shovill | 流程工程基线 | PE reads 检查、depth/coverage 控制、临时目录、日志、contig rename/filter、标准输出 | 不复制 Shovill 的 Perl 实现 |
| Unicycler | graph/bridging 输出基线 | SPAdes graph 使用、paths、depth normalization、anchor segments、bridge candidates、bridge quality/conflict 思路 | v1 不实现完整 hybrid assembler |
| Sparrowhawk | Rust 实现参考 | Rust module 边界、small bacterial genome 假设、rayon 并行、图结构数据设计 | 不把 Sparrowhawk-style 自研组装器设为 v1 主线 |
| Telomore | linear finishing 参考 | 用户指定 linear contig、左右端独立处理、terminal extension、QC、finalized assembly | 不在证据不足时强行延伸 |
| SKESA/skesa-rs | benchmark 与保守对照 | deterministic microbial assembly、k-mer extension 思想、Rust API 对照 | 不替代 SPAdes 主 backend |

工程原则：

- SPAdes 是 v1 backbone assembly 的主路径。
- Rust 负责 CLI、pipeline orchestration、二代 reads 前处理、SPAdes command planning、I/O 标准化、graph/bridging 文件生成、QC、benchmark 和可解释报告。
- 当前仓库中已有的 Rust k-mer/DBG 代码只可保留为 QC、测试夹具或后续研究资产，不作为 v1 backend 或产品路线。
- 所有参考工具只作为算法、流程和 benchmark 依据；不得直接复制不兼容许可证代码。
- `reference_tools/` 中的压缩源码归档作为只读参考材料；实现时只提炼接口、文件格式和流程设计，不直接 vendor 或复制参考工具代码。

## 3. 数据集与 truth 组织

当前阶段的开发测试数据来自 `reference_tools/Unicycler-main.zip`。7 个链霉菌样本和线性 terminal truth 暂不作为当前阶段验收输入，等 short-read-first front-end 和 graph/bridging 输出稳定后再接入。

当前阶段 fixture 来源：

```text
reference_tools/Unicycler-main.zip
├── Unicycler-main/sample_data/
│   ├── short_reads_1.fastq.gz
│   ├── short_reads_2.fastq.gz
│   ├── reference.fasta
│   ├── long_reads_low_depth.fastq.gz
│   └── long_reads_high_depth.fastq.gz
└── Unicycler-main/test/
    ├── test_assembly_graph.gfa
    ├── test_assembly_graph.fastg
    ├── test_assembly_graph.fastg.paths
    ├── test_assembly_graph_no_paths.gfa
    ├── test_bad_reads_1.fastq
    ├── test_bad_reads_2.fastq
    └── test_*_graph*.gfa
```

当前阶段使用规则：

- `sample_data/short_reads_1.fastq.gz` 和 `sample_data/short_reads_2.fastq.gz` 是 primary paired-end assembly fixture。
- `sample_data/reference.fasta` 只用于 assembly/front-end sanity check，不用于线性 terminal truth。
- `test/test_assembly_graph.gfa`、`test/test_assembly_graph.fastg` 和 `test/test_assembly_graph.fastg.paths` 用于 graph/path/parser/golden tests。
- `test/test_bad_reads_1.fastq` 和 `test/test_bad_reads_2.fastq` 用于 FASTQ/pair validation negative tests。
- 测试可以从 zip 中抽取到 `tests/fixtures/unicycler/`，但抽取脚本必须可复现，并记录 zip 内原始路径和 checksum。

后续线性阶段再启用 7 个链霉菌样本作为 Dev-core dataset，同时承担开发、验证和 pseudo-truth 角色。

推荐数据层级：

```text
rawdata/
├── sets/
│   ├── dev_test/
│   └── reference_validation/
├── fastq/
├── sra/
└── download_logs/
```

长期可整理为：

```text
datasets/
├── metadata/
│   ├── samples.tsv
│   ├── references.tsv
│   ├── reads.tsv
│   └── truth_telomeres.tsv
├── raw_reads/
├── references/
├── assemblies/
├── truth/
└── runs/
```

`samples.tsv` 最小字段：

```text
sample_id	species	strain	r1	r2	long_reads	long_platform	reference_fasta	reference_gff	provided_assembly	notes
```

`truth_telomeres.tsv` 最小字段：

```text
sample_id	replicon_id	side	expected_telomere_present	telomere_class	terminal_window	evidence_source	confidence
```

数据要求：

- 每个样本至少有 Illumina paired-end reads。
- reference replicon topology 需要能区分 linear、circular、unknown。
- terminal truth 必须按 left/right end 分开记录。
- 降采样和 benchmark 必须记录 seed、coverage、输入 checksum 和软件版本。

## 4. 目标架构

目录职责：

```text
src/
├── main.rs        # CLI parse and dispatch only
├── cli/           # clap command definitions
├── pipeline/      # assemble, inspect, finish, benchmark orchestration
├── runner/        # external tool runner: spades, shovill, telomore, skesa
├── preprocess/    # Shovill-style reads QC, trimming, depth reduction, read correction planning
├── graphio/       # SPAdes/Unicycler-like graph, path, depth and bridge-ready exports
├── bridge/        # bridge candidate/evidence schemas and later bridge application
├── qc/            # reads QC, coverage, k-mer evidence
├── io/            # FASTA/FASTQ/GFA/TSV/JSON I/O
├── postprocess/   # contig filtering, rename, stats
├── linear/        # terminal candidates, telodb, TIR, rescue
├── stats/         # run summary, benchmark plan, benchmark report
└── utils/         # command, checksum, tempdir, version capture
```

Dependency direction:

```text
cli -> pipeline -> runner/domain/io/stats
```

Rules:

- `main.rs` only parses CLI and dispatches.
- `runner/` owns external command construction, dry-run, logging and exit handling.
- `preprocess/` owns reads optimization policy; it must produce auditable files, never hidden temporary state only.
- `graphio/` owns conversion from SPAdes outputs to LineaBact standard graph/path/depth files.
- `bridge/` owns bridge schemas and evidence records; bridge application remains conservative and later-phase.
- `linear/` owns terminal detection/rescue decisions, not filesystem layout.
- `io/` and `stats/` serialize outputs; they should not drive domain decisions.
- `anyhow` is acceptable at orchestration boundaries; library-style errors use `thiserror`.

## 5. CLI 设计

### 5.1 `assemble`

Primary v1 command. Runs a Rust short-read-first front-end: Shovill-style reads optimization, SPAdes backend assembly, and Unicycler-like graph/bridging-ready export.

```bash
lineabact assemble \
  --r1 rawdata/fastq/S001/S001_1.fastq.gz \
  --r2 rawdata/fastq/S001/S001_2.fastq.gz \
  --backend spades \
  --preset streptomyces \
  --threads 16 \
  --memory-gb 64 \
  --target-depth 150 \
  --trim-adapters \
  --keep-graph-intermediates 1 \
  --downsample-seed 42 \
  --outdir runs/lineabact/S001_50x
```

Required behavior:

- Validate paired FASTQ input.
- Estimate depth and optionally reduce reads deterministically.
- Optionally trim adapters and plan read correction.
- Select SPAdes k-mers from read length and preset.
- Build and run SPAdes command with fixed resources.
- Capture command, stdout, stderr, versions and exit code.
- Normalize contigs, scaffolds, graph, path and bridging-ready files to the LineaBact output contract.

### 5.2 `inspect`

Runs linear-aware inspection on an existing assembly.

Current-stage status: deferred. Do not implement or require this command while using only Unicycler short-read fixtures.

```bash
lineabact inspect \
  --assembly runs/lineabact/S001_50x/backbone_contigs.fasta \
  --graph runs/lineabact/S001_50x/assembly_graph.gfa \
  --reference rawdata/sets/reference_validation/S001/reference.fasta \
  --outdir runs/lineabact/S001_inspect
```

Required behavior:

- Identify terminal candidates.
- Score left/right ends independently.
- Report graph dead-end evidence, k-mer support, telomere DB hits and TIR signals when available.

### 5.3 `finish`

Runs Telomore-inspired terminal rescue on user-specified linear contigs.

Current-stage status: deferred. Do not implement or test terminal rescue until linear fixtures and terminal truth are introduced.

```bash
lineabact finish \
  --assembly runs/lineabact/S001_50x/backbone_contigs.fasta \
  --graph runs/lineabact/S001_50x/assembly_graph.gfa \
  --r1 rawdata/fastq/S001/S001_1.fastq.gz \
  --r2 rawdata/fastq/S001/S001_2.fastq.gz \
  --linear-contig chromosome \
  --terminal-window 5000 \
  --terminal-min-support 3 \
  --outdir runs/lineabact/S001_finish
```

Required behavior:

- Process left and right ends independently.
- Prefer `unresolved` over unsupported extension.
- Record every attach, reject and trim decision in `terminal_events.jsonl`.

### 5.4 `stats benchmark-*`

Builds and executes benchmark plans.

```bash
lineabact stats benchmark-plan \
  --manifest rawdata/sets/reference_validation/manifest.tsv \
  --outdir benchmarks/generated/reference_validation

lineabact stats benchmark-run \
  --plan benchmarks/generated/reference_validation/benchmark_summary.json \
  --repeat-count 3
```

Full-roadmap benchmark plans should eventually include LineaBact, SPAdes, Shovill, Unicycler, SKESA/skesa-rs, Sparrowhawk and Telomore where applicable. Current-stage benchmark plans include only LineaBact, SPAdes, Shovill and Unicycler short-read-only behavior.
Current-stage recommended benchmark execution uses `benchmark-run --repeat-count 3`, and the CLI default should match that recommendation.

## 6. Output contract

Stable root outputs:

```text
lineabact_out/
├── backbone_contigs.fasta
├── finished_contigs.fasta
├── assembly_graph.gfa
├── assembly_stats.tsv
├── params.toml
├── run_summary.json
├── qc/
├── reads/
├── spades/
├── postprocess/
├── graph/
├── bridging/
├── linear/
├── logs/
└── benchmark/
```

Short-read front-end outputs:

```text
reads/reads_stats.tsv
reads/pair_check.tsv
reads/depth_estimate.tsv
reads/downsample_plan.tsv
reads/optimized_R1.fastq.gz
reads/optimized_R2.fastq.gz
```

SPAdes-backed backbone outputs:

```text
backbone_contigs.fasta
assembly_graph.gfa
postprocess/contig_stats.tsv
postprocess/rename_map.tsv
spades/contigs.fasta
spades/scaffolds.fasta
spades/assembly_graph_with_scaffolds.gfa
spades/assembly_graph.fastg
spades/contigs.paths
spades/scaffolds.paths
spades/before_rr.fasta
logs/backend.cmd.txt
logs/backend.stdout.log
logs/backend.stderr.log
```

Unicycler-like graph/bridging-ready outputs:

```text
graph/raw_spades_graph.gfa
graph/cleaned_graph.gfa
graph/overlap_trimmed_graph.gfa
graph/segments.tsv
graph/links.tsv
graph/paths.tsv
graph/depth.tsv
graph/anchor_segments.tsv
graph/graph_qc.tsv
bridging/bridge_manifest.json
bridging/spades_path_bridge_candidates.tsv
bridging/bridge_evidence.jsonl
bridging/bridge_conflicts.tsv
bridging/bridge_decisions.jsonl
bridging/bridged_graph.gfa
```

Linear-aware outputs:

Current-stage status: deferred. These files are listed as the later linear-genome output contract and are not required for Unicycler-fixture development.

```text
linear/terminal_candidates.tsv
linear/terminal_kmer_support.tsv
linear/telomere_db_hits.tsv
linear/tir_candidates.tsv
linear/topology_report.tsv
linear/terminal_extensions.tsv
linear/terminal_qc.tsv
linear/terminal_events.jsonl
linear/unresolved_ends.tsv
```

Compatibility note:

- Existing prototype outputs may currently write some terminal files at the output root.
- v1 should migrate terminal-specific files into `linear/`.
- Any migration must update `run_summary.json`, golden tests and documentation in the same change.

## 7. Roadmap

### Phase 0: Unicycler fixture 标准化与参考基线

Goal: turn `reference_tools/Unicycler-main.zip` into reproducible current-stage development fixtures.

Tasks:

- Define `tests/fixtures/unicycler/manifest.tsv` from zip internal paths.
- Extract or reference `sample_data/short_reads_1.fastq.gz`, `sample_data/short_reads_2.fastq.gz` and `sample_data/reference.fasta`.
- Extract or reference graph fixtures: `test_assembly_graph.gfa`, `test_assembly_graph.fastg`, `test_assembly_graph.fastg.paths` and selected small GFA negative cases.
- Run SPAdes, Shovill and Unicycler short-read-only behavior checks where available.
- Record zip filename, zip checksum, internal path, extracted checksum, command, tool version and resource usage.

Acceptance:

- The current-stage fixture manifest can recreate all test inputs from `reference_tools/Unicycler-main.zip`.
- The paired-end fixture can drive `assemble --dry-run`, mock SPAdes output normalization and graph/path parsing tests.
- No line in this phase requires linear topology, telomere truth or terminal rescue.

### Phase 1: Rust short-read-first front-end MVP

Goal: implement a reliable Rust front-end around Illumina paired-end reads and SPAdes.

Tasks:

- Implement `lineabact assemble --backend spades` as the main path.
- Add `--dry-run` to emit planned commands and parameter files without running SPAdes.
- Validate R1/R2 existence, gzip readability, pair count and basic read stats.
- Add deterministic depth reduction with fixed seed.
- Generate `reads/` optimization plan and optimized reads outputs.
- Run SPAdes with explicit `--isolate`, `--threads`, `--memory`, `--tmp-dir` and k-mer policy.
- Copy raw SPAdes outputs into `spades/`.
- Write `params.toml`, `run_summary.json`, `assembly_stats.tsv`, `reads/*` and backend logs.

Acceptance:

- Small fixture test passes without SPAdes using mock/dry-run mode.
- At least 2 reference-validation samples run at 50x with SPAdes backend.
- Re-running the same input with the same seed produces stable read-selection, command, params and output metadata.

### Phase 1.5: Shovill-style reads optimization and QC

Goal: close the engineering gap against Shovill before adding graph/bridging and linear-specific logic.

Tasks:

- Add adapter trimming policy and optional read correction planning.
- Estimate genome size/depth and reduce reads to target depth.
- Select SPAdes k-mers from read length, preset and user override.
- Add Shovill-style contig filtering, rename and parseable FASTA headers.
- Add read QC, k-mer spectrum QC and assembly k-mer recall.
- Capture wall time, peak memory and temp disk use.
- Compare LineaBact SPAdes-backed output to raw SPAdes and Shovill.

Acceptance:

- 50x and 100x benchmark reports exist for the dev-core dataset subset.
- N50, contig count, total length and genome fraction are within expected range versus SPAdes/Shovill.
- `reads/*`, QC TSV/JSON and contig postprocess schemas have golden tests.

### Phase 2: Unicycler-like graph and bridging-ready exports

Goal: standardize SPAdes graph outputs into downstream files usable by later graph cleanup, bridging and manual Bandage inspection. Current-stage work stops at graph/bridge-ready outputs and does not perform linear terminal analysis.

Tasks:

- Preserve SPAdes raw `assembly_graph_with_scaffolds.gfa`, `assembly_graph.fastg`, `contigs.paths`, `scaffolds.paths` and `before_rr.fasta`.
- Convert graph segments, links and paths into stable `graph/*.tsv` tables.
- Normalize segment depth and record copy/depth evidence for later anchor selection.
- Emit `graph/cleaned_graph.gfa` and `graph/overlap_trimmed_graph.gfa` when cleanup is applied.
- Identify anchor segments using length, depth consistency and graph connectivity.
- Derive SPAdes-path bridge candidates from path records and write bridge evidence without applying low-confidence bridges.

Bridge candidate fields:

```text
bridge_id
source
left_segment
right_segment
path_segments
bridge_length
depth_agreement
path_self_contained
insert_size_penalty
quality
status
reason
```

Acceptance:

- Graph exports can be loaded without parsing SPAdes-specific headers directly.
- Bridge candidates are deterministic and sorted by quality.
- Conflicting, low-quality or long-repeat bridges are recorded but not silently applied.
- GFA outputs remain viewable in Bandage/Bandage-NG.

### Phase 2.5: Deferred terminal candidate detection

Goal: identify likely linear replicon ends from assembly and graph exports.

Status: deferred until after the short-read-first front-end and Unicycler-like graph/bridging outputs are stable. This phase must not be part of current-stage acceptance.

Tasks:

- Parse `backbone_contigs.fasta`, `graph/segments.tsv`, `graph/links.tsv` and `graph/depth.tsv`.
- Compute left/right graph degree, terminal window GC, N rate, local k-mer support and depth proxy.
- Score each contig end independently.
- Output `linear/terminal_candidates.tsv` and `linear/terminal_kmer_support.tsv`.

Terminal score v1:

```text
terminal_score =
    graph_dead_end_score
  + terminal_kmer_support_score
  + coverage_consistency_score
  + telodb_match_score
  + tir_partner_score
  - repeat_or_error_penalty
```

Acceptance:

- Terminal candidate recall can be measured against reference terminal truth.
- Obvious internal contig ends are reported as low confidence or non-terminal.
- Every candidate has machine-readable reasons.

### Phase 3: Deferred telomere DB, TIR and topology

Goal: add actinomycete-specific terminal priors.

Status: deferred. Requires linear genome fixtures and terminal truth, which are explicitly out of current-stage scope.

Tasks:

- Build a minimal local telomere DB from reference terminal windows.
- Scan terminal windows for telomere-like k-mer containment.
- Detect terminal inverted repeats by comparing left/right terminal windows.
- Classify topology as `linear_with_tir`, `linear_candidate`, `linear_plasmid_candidate`, `circular_candidate`, `fragmented` or `ambiguous`.

Acceptance:

- Synthetic `linear_with_tir` and circular negative controls are covered.
- Circular controls do not produce high-confidence linear-with-TIR calls.

### Phase 4: Deferred Telomore-inspired terminal rescue

Goal: conservatively recover missing terminal bases only when evidence is strong.

Status: deferred. Do not implement or test terminal rescue in the current Unicycler-fixture stage.

Tasks:

- Accept user-specified linear contigs from CLI or file.
- Recruit terminal-supporting reads for each end independently.
- Build end-specific consensus.
- Attach, reject or trim extensions based on support, overlap quality and repeat-like evidence.
- Write `finished_contigs.fasta`, `terminal_extensions.tsv`, `terminal_qc.tsv`, `terminal_events.jsonl` and `unresolved_ends.tsv`.

Terminal states:

```text
recovered
unresolved
low_support
rejected_repeat_like
```

Acceptance:

- Synthetic positive cases reach >= 99% terminal extension identity.
- Circular negative controls do not produce high-confidence extensions.
- Every terminal decision is explainable in JSONL.

### Phase 5: Deferred long-read-assisted bridging and validation

Goal: use ONT/PacBio as auxiliary evidence for graph bridging and terminal validation, not as a required v1 input.

Status: deferred. The Unicycler archive includes low/high depth long-read fixtures, but current-stage development ignores long reads unless they are needed only for future schema planning.

Tasks:

- Support `--long` and `--long-platform ont|pacbio-hifi|pacbio-clr`.
- Produce long-read bridge evidence compatible with `bridging/bridge_evidence.jsonl`.
- Validate terminal boundaries and TIR signals with long reads.
- Report Illumina-only versus long-read-supported evidence.

Acceptance:

- Long reads improve confidence where available.
- Low-quality long-read evidence cannot force terminal extension by itself.

## 8. Benchmark design

Backbone benchmark tools:

```text
SPAdes
Shovill
LineaBact SPAdes-backed workflow
Unicycler short-read-only mode / graph output behavior
SKESA/skesa-rs
Sparrowhawk
```

Terminal benchmark tools:

```text
LineaBact finish
Telomore
manual/reference truth
```

Current-stage benchmark scope:

- Run only short-read assembly/front-end benchmarks using Unicycler sample data.
- Compare LineaBact front-end behavior against SPAdes, Shovill and Unicycler short-read-only graph behavior where tools are available.
- Do not run terminal recovery, Telomore or linear-genome truth benchmarks in the current stage.

Metrics:

- Runtime, peak memory and temp disk.
- `lineabact_vs_shovill` runtime ratio, defined as `LineaBact elapsed_seconds / Shovill elapsed_seconds`; `<= 1.0` means LineaBact has reached Shovill speed level on that sample.
- Contig count, total length, longest contig, N50 and genome fraction.
- Misassemblies and read mapping rate where reference exists.
- Assembly k-mer recall and missing solid k-mers.
- Graph dead-end count, component count and depth-filtered graph size.
- Anchor segment count and bridge candidate quality distribution.
- SPAdes path bridge recall where reference paths are known.
- Deferred linear metrics: terminal candidate recall/precision, left/right arm completeness, telomere/TIR recovery, extension identity, false-positive extension rate and unresolved end count.

Benchmark rules:

- Use fixed downsampling seed.
- Separate backbone quality from terminal recovery quality.
- Record command, software version, input checksum and parameters in machine-readable outputs.
- For runtime acceptance against Shovill, record `lineabact_vs_shovill` on at least 3 repeats per sample and evaluate the median runtime ratio rather than a single run.
- Treat backbone assembly quality acceptance as a hard gate, not descriptive context:
  - `lineabact` N50 must be at least 90% of `shovill`, with a preferred target of >= 95%.
  - `lineabact` genome fraction must not be lower than `shovill`.
  - `lineabact` misassembly count must not exceed `shovill`.
  - `lineabact` mismatch and indel rates must not be worse than `shovill`.
  - A sample is not accepted if runtime improves but misassemblies materially increase.

## 9. Testing strategy

Unit tests:

- FASTA/FASTQ/GFA parsing.
- FASTQ pair validation.
- Downsampling determinism.
- SPAdes output discovery and normalization.
- SPAdes `contigs.paths` / `scaffolds.paths` parsing.
- Segment depth normalization and anchor segment selection.
- Bridge candidate quality, conflict detection and deterministic sorting.
- Contig filtering and rename.
- Deferred linear tests: terminal window extraction, terminal score, telomere DB containment, TIR detection and attach/reject/trim decisions.

Integration tests:

- `assemble --dry-run`.
- `assemble --backend mock`.
- `assemble` mock SPAdes output normalization into `spades/`, `graph/` and `bridging/`.
- Unicycler fixture paired-end smoke test using `sample_data/short_reads_1.fastq.gz` and `sample_data/short_reads_2.fastq.gz`.
- Unicycler graph fixture parser test using `test/test_assembly_graph.gfa` and `test/test_assembly_graph.fastg.paths`.
- `stats benchmark-plan` and `stats benchmark-run` smoke tests.

Golden tests:

```text
run_summary.json
params.toml
assembly_stats.tsv
graph/segments.tsv
graph/links.tsv
graph/paths.tsv
graph/anchor_segments.tsv
bridging/bridge_manifest.json
bridging/spades_path_bridge_candidates.tsv
bridging/bridge_evidence.jsonl
```

Deferred golden tests:

- `linear/terminal_candidates.tsv`
- `linear/telomere_db_hits.tsv`
- `linear/tir_candidates.tsv`
- `linear/topology_report.tsv`
- `linear/terminal_events.jsonl`

Current-stage negative controls:

- Bad paired-end FASTQ fixtures from `Unicycler-main/test/test_bad_reads_1.fastq` and `test_bad_reads_2.fastq`.
- GFA without paths from `Unicycler-main/test/test_assembly_graph_no_paths.gfa`.
- Bridge conflict or low-quality bridge cases from selected small Unicycler GFA fixtures.

Deferred negative controls:

- Circular genome control.
- Repeat-rich terminal-like control.
- Low-support terminal evidence control.

## 10. Immediate priorities

### Priority 0: Restore implementation baseline

- Fix current compile errors before adding features.
- Synchronize CLI args, `RunSummary` schema and benchmark runner argument construction.
- Run `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings` and `cargo test` in the `LineaBact` conda environment.

### Priority 1: Implement the short-read-first front-end

- Add explicit backend selection: `spades` and `mock`.
- Make `spades` the documented default backend for v1.
- Add Shovill-style reads QC, target-depth reduction, optional trimming policy and k-mer selection.
- Add SPAdes dry-run, command logging, version capture and output discovery.
- Remove the current Rust DBG path from the documented main command path; keep reusable k-mer code only where it supports QC or terminal evidence.
- Use only `reference_tools/Unicycler-main.zip` fixtures for current-stage development tests.

### Priority 2: Standardize graph and bridging-ready files

- Copy raw SPAdes outputs into `spades/`.
- Parse GFA/path/depth outputs into `graph/segments.tsv`, `graph/links.tsv`, `graph/paths.tsv` and `graph/depth.tsv`.
- Produce `graph/anchor_segments.tsv` and `bridging/spades_path_bridge_candidates.tsv`.
- Record bridge evidence, conflicts and decisions even when bridges are not applied.
- Validate parsers and golden outputs against Unicycler graph fixtures.

### Priority 3: Stabilize output contract

- Ensure all outputs are represented in `run_summary.json`.
- Add schema/golden tests for reads, graph and bridging outputs.
- Do not require linear/terminal outputs for current-stage completion.

### Priority 4: Build reference-tool benchmarks

- Generate current-stage benchmark plans for SPAdes, Shovill, Unicycler short-read-only behavior and LineaBact.
- Defer SKESA/skesa-rs, Sparrowhawk and Telomore benchmark entries unless needed as optional references.
- Keep benchmark-only tools phase-gated in the conda environment.
- Do not let benchmark dependencies become mandatory for basic development.

## 11. Current-stage success definition

Current-stage milestone is complete when:

- The Unicycler paired-end sample fixture can drive `assemble --dry-run`, mock backend, and SPAdes-backed runs where SPAdes is installed.
- Shovill-style reads QC, logs, normalized contigs and machine-readable summaries are stable.
- Unicycler-like graph and bridging-ready files are produced from SPAdes outputs.
- SPAdes graph paths can generate deterministic bridge candidates without applying unsafe bridges.
- Unicycler GFA/FASTG/path fixtures are covered by parser and golden tests.
- `run_summary.json` records fixture source paths, checksums, commands, versions and generated files.
- `lineabact_vs_shovill` is stably recorded for the Unicycler fixture and the later reference-validation subset, using at least 3 repeats per sample and the median runtime ratio for acceptance.
- For the same accepted samples, backbone assembly quality against `shovill` passes the hard gates:
  - N50 >= 90% of `shovill`, with >= 95% as the preferred steady-state target.
  - Genome fraction is not lower than `shovill`.
  - Misassembly count is not higher than `shovill`.
  - Mismatch and indel rates are not worse than `shovill`.
  - Faster runtime does not override a quality regression.
- No linear genome, telomere, terminal candidate or terminal rescue output is required.

## 12. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Rewriting assembler internals consumes the project | Keep SPAdes as v1 backbone; restrict Rust k-mer/graph code to QC, tests or research modules |
| Output drifts between modules | Define one output contract and enforce it with golden tests |
| Graph/bridge files become incompatible with SPAdes or Bandage | Preserve raw SPAdes files and test normalized GFA/path outputs against small fixtures |
| Bridge candidates introduce misassemblies | Emit bridge evidence first; apply only later with quality/conflict thresholds |
| Terminal rescue over-extends repeats | Prefer `unresolved`; require support, overlap quality and repeat checks |
| Shovill parity is poor | Fix preprocessing, filtering, naming and resource handling before biology-specific extensions |
| Telomore assumptions do not fit fragmented assemblies | Require user-specified linear contigs and report low-confidence fragmented cases |
| Long reads introduce platform-specific errors | Treat ONT/PacBio as auxiliary evidence and report platform separately |
| Benchmark is not reproducible | Record seed, versions, checksums, resources and full commands |

## 13. Documentation governance

- Product scope, roadmap, output contract and benchmark strategy live here.
- Repository execution rules live in `AGENTS.md`.
- Task-specific implementation requirements live in `docs/playbooks/`.
- If this plan conflicts with a playbook, update the higher-level decision here first, then sync the playbook.

## 14. References

[^spades]: SPAdes Genome Assembler. https://github.com/ablab/spades

[^shovill]: Shovill: assemble bacterial isolate genomes from Illumina paired-end reads. https://github.com/tseemann/shovill

[^unicycler]: Unicycler local reference archive: `reference_tools/Unicycler-main.zip`; upstream project: https://github.com/rrwick/Unicycler

[^sparrowhawk]: Sparrowhawk: short-read assembler for bacterial genomics based on a de Bruijn graph written in Rust. https://github.com/bacpop/sparrowhawk

[^telomore]: Telomore: tool for extending linear bacterial replicons to capture telomeres. https://github.com/dalofa/telomore

[^skesa]: skesa-rs: Rust port of NCBI's SKESA genome assembler. https://docs.rs/skesa-rs/latest/skesa_rs/
