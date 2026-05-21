# Terminal Rescue Development Playbook

## Purpose
Standardized implementation rules for the terminal rescue functionality in LineaBact.

## Core Requirements

1. **Left/Right Independence**
   - Always process left and right ends of a contig separately
   - Avoid shared state between ends unless absolutely necessary

2. **Read Recruitment**
   - Recruit terminal-supporting reads based on mapping to contig ends
   - Define clear window size for recruitment (configurable via `--terminal-window`)

3. **Consensus Building**
   - Build end-specific consensus from recruited reads
   - Record coverage and quality information

4. **Attach / Reject Decision**
   - Only attach consensus if overlap quality, sequence consistency, and read support all meet thresholds
   - Clearly log the reason for rejection

5. **Support-based Trimming**
   - After attaching, perform trimming based on read support depth
   - Record how many bases were trimmed and why

6. **State Reporting**
   Every terminal end must be assigned one of the following states:
   - `recovered`: Successfully extended with high confidence
   - `unresolved`: Insufficient evidence to extend
   - `low_support`: Evidence exists but below threshold
   - `rejected_repeat_like`: Likely caused by repeats, rejected

7. **Logging & Output**
   - Every decision must be written to `terminal_events.jsonl`
   - Update `terminal_extensions.tsv` and `terminal_qc.tsv`
   - Maintain explainability for every operation

8. **Safety**
   - Be conservative: when in doubt, prefer `unresolved` over incorrect extension
   - Always consider negative control (circular genomes) behavior
