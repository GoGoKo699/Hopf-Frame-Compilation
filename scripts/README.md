# Executable entry points

[← Repository landing page](../README.md) · [Verification map](../docs/VERIFICATION.md)

The scripts provide a readable orientation, resource ledgers, source checks,
exact verification receipts, and rendered presentation checks. They do not
replace the analytic proof.

## Research validation

Run commands from the repository root. During an edit, select the affected
modules below; the [verification catalogue](../docs/reference/VERIFICATION_CATALOGUE.md#catalogue-index)
owns their detailed evidence contracts and maps the other research components.
Use [WORKSPACE](../WORKSPACE.md) for the current question, not this command map.

| Changed material | Focused checks |
| --- | --- |
| Working documents, links and mathematical syntax | [`test_repository_navigation`](../tests/test_repository_navigation.py), [`test_reviewer_narrative`](../tests/test_reviewer_narrative.py), [`test_presentation`](../tests/test_presentation.py) |
| [Augmented field recoupler](../research/depth/AUGMENTED_FIELD_RECOUPLER.md) | [`test_augmented_field_recoupler`](../tests/test_augmented_field_recoupler.py) |
| [Projective field recoupler](../research/depth/PROJECTIVE_FIELD_RECOUPLER.md), both native allocations | [`test_projective_field_recoupler`](../tests/test_projective_field_recoupler.py), [`test_coupled_projective_recoupler`](../tests/test_coupled_projective_recoupler.py) |
| [Weak complete-probe domain](../research/depth/WEAK_PROBE_PAULI_DOMAIN.md) | [`test_weak_probe_pauli_domain`](../tests/test_weak_probe_pauli_domain.py) |
| [Clifford-sector degree barrier](../research/depth/FIELD_SECTOR_DEGREE_BARRIER.md) and signed helper tableaux | [`test_field_sector_degree`](../tests/test_field_sector_degree.py), [`test_field_sector_tableaux`](../tests/test_field_sector_tableaux.py) |

For example, the working-document checks run with:

```bash
python -m unittest tests.test_repository_navigation tests.test_reviewer_narrative tests.test_presentation
```

Use the same `python -m unittest tests.<module>` form for the affected scientific
suites. A bounded fixture checks its stated identities and implementation;
the linked proof retains responsibility for all-size claims and scope.

Before merging, run the complete local gate on the final reviewed tree:

```bash
python -m compileall -q compiler_robust_hopf scripts tests validate.py
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python scripts/complex_coarse_native_example.py
python scripts/residual_qbp_native_example.py --q 16
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
python scripts/unified_resource_ledger.py --n 12 --format json
python scripts/strict_zero_echo_ledger.py --n 12 --format json
```

Complete focused [visual checks](#visual-checks) for changed mathematics and
pages before the memory-heavy numerical suite; full presentation CI covers
the complete corpus. Focused tests support iteration and do not replace this
final gate or required CI. Record the tested commit/tree and command results
in the review record. Expected scientific receipts stay unchanged unless the
underlying evidence is intentionally regenerated and reviewed.

## Technical walkthrough

```bash
python scripts/reviewer_walkthrough.py
```

This runs ten representative checks in the same order as the narrative:

1. recursive and addressed two-qubit frames agree;
2. canonical domains and singular coordinates are interpreted consistently;
3. one prepared column is insufficient;
4. the strict-zero four-sector identity holds;
5. the borrowed suffix bit is restored;
6. the explicit coherent router matches the tail direct sum;
7. the phase diagonal is one UCG;
8. the correct workspace schedule is selected;
9. the simultaneous workspace peak respects the supplied budget;
10. the low-workspace and maximal-cut inequalities hold.

## All-workspace ledger

```bash
python scripts/unified_resource_ledger.py --n 12
```

This displays the selected schedule and transparent integer resource terms over
representative clean-workspace budgets.  Use `--format json` for a
machine-readable record.

## Strict-zero ledger

```bash
python scripts/strict_zero_echo_ledger.py --n 12
```

This displays the per-depth borrowed-suffix UCG widths, predicate terms, and the
complete strict-zero size/depth proxies.

## Native coarse-frame example

```bash
python scripts/coarse_frame_native_example.py
python scripts/coarse_frame_native_example.py --case singular --format json
```

This bounded two-qubit example emits elementary Clifford+T words, checks
complete preparation and coherent dirty-input gradient means, and reports
literal counts against the original protocol with the same observable.
Its exact finite-size target is not a generic fine-precision compiler or an
advantage experiment. See the [scope and word proof](../supplements/state_based_qbp/NATIVE_COARSE_QBP.md).

## Native complex coarse-frame example

```bash
python scripts/complex_coarse_native_example.py
python scripts/complex_coarse_native_example.py --case singular --format json
```

This bounded two-qubit complex example propagates native words through both
magnitude and phase gradient streams, checks complete preparation and dirty-input
score operators, and reports gate counts against the original protocol. It also
checks fixed-record histogram decoding, without sampling. Numerical errors are
finite floating checks; the exact target and word construction are described in
the [complex fixture proof](../supplements/state_based_qbp/NATIVE_COMPLEX_COARSE_QBP.md).
It is not a general fine-residual compiler or an advantage experiment.

## Native residual QBP ledger

```bash
python scripts/residual_qbp_native_example.py
python scripts/residual_qbp_native_example.py --q 16 --format json
```

This emits both fine-residual gradient streams for a fixed complex
one-qubit target and reports exact coefficient, coarse-distance, gate,
and gradient-bias certificates. It also demonstrates the exact histogram
decoders. It performs no statevector propagation or sampling; the bounded
native numerical checks live in the test suite. See the
[integration proof and scope](../supplements/state_based_qbp/NATIVE_RESIDUAL_QBP.md).

## Native ordered-bank components

```bash
python scripts/ordered_banks_native_example.py --n 2
python scripts/ordered_banks_native_example.py --n 3
```

This regenerates native right-readout and coarse-prefix words and reports
their gate counts, hashes, and floating errors on two vectors per component.
The script supports only n=2 and n=3; these vector checks are not all-column
verification or endpoint synthesis. The
[ordered-bank proof and tests](../research/endpoint/ORDERED_BLOCK_BANKS.md#10-remaining-native-theorem-and-reproducibility)
supply the component identities and wider evidence. The fine left source is
not emitted. Add `--write-gates /tmp/hopf-ordered-banks` to export local gate
lists; generated lists are not additional canonical evidence.

## Upstream synchronization

```bash
python scripts/check_upstream_sync.py --offline
```

The offline mode validates the recorded source schema and local lineage without
network access.  The tracked upstream commits and scientific reconciliation are
listed in [`SYNC.md`](../SYNC.md) and
[`provenance/upstream.json`](../provenance/upstream.json).

## Complete deterministic suite

```bash
python validate.py
```

The complete test map is in [`tests/README.md`](../tests/README.md).

## Fault-tolerant checks

Run `python scripts/verify_fault_tolerant.py` to reproduce four focused exact
source/kernel and rational-resource suites. See the
[scope and receipt guide](../verification/fault_tolerant/README.md). Temporary
outputs are used by default, leaving the expected evidence files unchanged.

## Visual checks

The optional `check_presentation.py` checker renders the diagrams and complete
Markdown pages, including tables, inline mathematics, and display equations.
It checks MathJax SVG and native MathML, including unsupported numbered rows
and vertically stacked equation glyphs, and saves desktop and narrow-screen
previews with measured layout checks.
The [rendering guide](../assets/README.md#rendering-checks) gives the separate
browser and MathJax dependencies and the reproducible command. These previews
model GitHub-style rendering; they do not reproduce GitHub's private client.
