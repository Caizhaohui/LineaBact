# K-mer & De Bruijn Graph Playbook

## Purpose
Standardize the implementation of k-mer processing and de Bruijn Graph (DBG) construction in LineaBact, ensuring correctness, performance, and compatibility with terminal rescue requirements.

## Core Requirements

### 1. K-mer Encoding & Canonicalization
- Use **packed 2-bit encoding** (u64 or u128) for k-mers when k ≤ 64.
- Strictly implement **canonical k-mer** representation (forward vs reverse-complement, always store the lexicographically smaller one).
- **Critical**: Preserve orientation/direction information. Do **not** lose strand information during canonicalization, as it is required for correct left/right terminal end identification.
- Provide clear functions:
  - `canonical(kmer)`
  - `reverse_complement(kmer)`
  - `is_canonical(kmer)`

### 2. K-mer Counting
- Support efficient k-mer counting (recommend using a hashmap or specialized structure like `hashbrown` or custom implementation).
- Record both **forward** and **canonical** counts when necessary.
- Support downsampling and k-mer spectrum analysis for coverage estimation.
- Handle high-coverage data gracefully (avoid memory explosion on `scl_L1` full coverage).

### 3. De Bruijn Graph Construction
- Build DBG from canonical k-mers.
- Each edge should preferably store:
  - The k-mer itself
  - Coverage / count
  - Direction / orientation metadata (important for terminal rescue)
- Support single-k mode in v1 (multi-k can be considered later).

### 4. Graph Simplification (Phase 1.5+)
- Implement **conservative** simplification strategies:
  - Tip removal (low coverage tips)
  - Bulge removal (coverage-aware)
  - Low-frequency edge filtering
- **Important Constraint**: Graph simplification must **not** destroy terminal evidence needed for later terminal rescue. Be conservative on contig ends.

### 5. Direction & Orientation Handling (Critical for LineaBact)
- Because LineaBact needs to distinguish **left and right ends** of linear contigs, the graph must retain sufficient orientation information.
- Unitig/contig extraction must record:
  - First k-mer
  - Last k-mer
  - Overall path direction
- Avoid designs that lose strand information after canonicalization.

## Implementation Guidelines

- Prefer clarity and correctness over premature optimization in early phases.
- Use `rayon` for parallel k-mer counting when appropriate.
- Add comprehensive unit tests for:
  - K-mer encoding / decoding
  - Canonicalization correctness (including palindromic k-mers)
  - Reverse complement
  - DBG construction on synthetic data
- For high GC content (*Streptomyces* ~72%), pay attention to k-mer distribution bias.

## Common Pitfalls to Avoid

- Losing strand/orientation information after canonical k-mer conversion.
- Over-aggressive graph simplification that removes terminal supporting edges.
- Poor memory management on high-coverage datasets.
- Non-deterministic behavior in k-mer ordering or graph traversal.

## Recommended Testing

- Unit tests with synthetic k-mer sets.
- Small real-data smoke tests on `scl_L1` 50x downsampled data.
- Verify that terminal ends still have sufficient supporting edges after simplification.
