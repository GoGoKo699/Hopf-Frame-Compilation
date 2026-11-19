# Research evidence, 19 November 2026

[Current question](../../../WORKSPACE.md) · [Resource frontier](../../../docs/OPEN_PROBLEM.md) · [Research index](../../README.md)

This archive preserves the original arguments, finite checks and negative results behind the November research handoff. The dated reports use the scientific baseline at [ae613bf](https://github.com/GoGoKo699/Hopf-Frame-Compilation/tree/ae613bf8e8a38f6d0407c25910b83d69f01a2733). Consult the current question and resource frontier above for theorem scopes and the active research direction. Historical proposals for the next experiment can predate stronger results in the same archive.

| Source path prefix | Contents |
| --- | --- |
| `early/`, `encoder/`, `joint/` | Augmented field decoder, selected matching, two-batch composition, two-layer encoders, joint readout, and read-only equality. The early field checks include the actual inverse and arbitrary occupied banks. |
| `haar/`, `prefix/` | Complete positive Haar frames, rotation banks, prefix-shared rotations, sparse Haar defects, and the hypotheses on node signs and nested defects. |
| `scalar/`, `prefix/proofs/` | Returned prescribed scalar tables and their transfer to complete frames, table density, the marked-edge construction and the counting argument outside the two widest layers. |
| `mask/`, `zero/` | Mask and moving-address zero-test transfers, trace inversion, finite native examples, and the quartic basis-mixing witness. The transfer proofs specify the conditions on the hidden sectors and bad events. |
| `four/`, `supplements/` | Four-layer zero-clean scalar and rank gaps, one-clean precision and five-layer restrictions, the unabridged scalar-gap proof, Moore quartic witness, and terminal mask cuts. |
| `unreviewed/` | Index of 24 exploratory leads with surviving calculations or finite fixtures. Their presence here records which directions were tried. |
| `audits/`, `receipts/portable/` | Six closeout reviews, historical replay outputs and command records. Source pointers in the reviews use the originating handoff's file names. |

## Read and replay

Run the commands below from this directory. Python 3.11 or later can inspect, verify and extract all sources. Full numerical replay also requires NumPy. The recorded environment used Python 3.12.14 and NumPy 2.3.5.

```bash
python evidence.py list scalar/proofs
python evidence.py show scalar/proofs/single_feedback_pencil.md
python evidence.py show unreviewed/README.md
python evidence.py verify
python evidence.py replay prefix field_independent
python evidence.py replay
```

Use `python evidence.py extract DIR` to restore all 263 original paths to a new directory, including the historical aliases that are needed for neighboring imports. The [inventory](inventory.json) identifies each original package by filename and SHA-256, maps all aliases, and records the digest of each unique source. [sources.zip](sources.zip) stores the 247 distinct contents once each. Package prefixes in the table stand in for the corresponding original ZIP filenames. Seven relative links in the historical scalar dirty-sum appendix resolve at its [repository home](../../../docs/DIRTY_SUM_COMPRESSION.md).

The [runner](evidence.py) verifies the source digests and uses a new temporary directory for each suite. The eight original packet runners invoke 81 finite checks, including shared fixtures. Six more entry points reproduce the early field, matching and composition checks. Two supplemental checks cover the Moore witness and terminal cuts. A seventeenth entry reruns the archived four-layer checker against the unabridged proof's historical fixture. The proofs carry the all-size statements.
