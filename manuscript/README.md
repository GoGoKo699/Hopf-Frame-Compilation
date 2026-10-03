# Manuscript guide

[Publication scope](PUBLICATION_SCOPE.md) · [Read the argument](../REVIEW.md) · [Verification](../docs/VERIFICATION.md)

Working title: **Exact and Fault-Tolerant Compilation of Hopf Differential Frames**.

The [publication scope](PUBLICATION_SCOPE.md) fixes four principal
results and their corollaries. The exact theorem covers every $m\geq0$
clean-workspace budget; the sufficient-clean and one-clean T-count
constructions keep their distinct clean/dirty reservations and return
guarantees. Result D gives simultaneous T-count and T-depth bounds for
complete real frames, matching both resources in explicit
precision/workspace ranges. The [bounded internal audit](../docs/CORE_CLAIM_AUDIT.md)
is complete without an unresolved claim-level blocker, and the selected
scientific scope is frozen. Final manuscript writing remains on hold.

## Scientific contribution

The paper develops a compiler capability: preserving an entire
operationally required frame while controlling exact entangling cost,
fault-tolerant precision cost, and initialized versus borrowed workspace.
The one-clean diagonal and general one-target multiplexor results also cover
standard operator families beyond Hopf frames. The earlier two-clean
constructions remain useful at their slightly smaller dirty reservations.

The principal synthesis comparisons are Yuan–Zhang's
[exact state-preparation tradeoffs](https://quantum-journal.org/papers/q-2023-03-20-956/)
(Quantum 7, 956), Low–Kliuchnikov–Schaeffer's
[dirty-workspace tradeoffs](https://quantum-journal.org/papers/q-2024-06-17-1375/)
(Quantum 8, 1375), and Gosset–Kothari–Wu's
[optimal T-count](https://quantum-journal.org/papers/q-2026-07-22-2168/)
(Quantum 10, 2168). The manuscript must explain the new complete-operator
guarantees relative to these results, with the attribution in the source map.

The [related-work comparison](../docs/RELATED_WORK.md) includes newer synthesis
results and distinguishes their input contracts and workspace assumptions.
The opening pages should explain why the prescribed frame
requires more than state preparation, identify the new proof mechanisms,
and state Results A–D with their resource assumptions. QBP provides the
operational motivation and consequence. The exact and fault-tolerant models
remain parts of one operator-compilation story. The constant-clean count
endpoint and unrestricted large-width T-depth remain discussion questions;
neither is a prerequisite for the selected paper.

## Paper order

1. Prescribed completion versus one-column state preparation, with the
   two-qubit readout example.
2. Shared Hopf-frame and complete-input compiler contract.
3. Exact all-workspace size/depth theorem and its proof mechanisms.
4. Sufficient-clean matching T theorem, followed by the one-clean real-frame
   construction and its grouped, banked, diagonal, general multiplexor,
   and phase-dressed complex-magnitude corollaries. Keep the distinct dirty
   reservations and the separately proved two-clean T-depth assumptions.
5. Uniform same-circuit T-count and T-depth tradeoffs for complete real
   frames, including matching intervals and the low-precision improvement.
   Keep their two-clean reservation and full-input error contract explicit;
   distinguish T-depth from total depth and retain older schedules as
   applicable fallbacks.
6. Exact and approximate fixed-parameter QBP consequences, including complete
   complex gradients and the costs of classical preprocessing and output.
7. The two remaining resource gaps and the limits of the circuit model.

Detailed schedules, source preparation, reversible lookup, error estimates,
and resource sums form the technical appendices. The exact toolkit credits
Yuan and Zhang; the [source map](../docs/SOURCE_MAP.md) gives the complete
attribution. Later component refinements need not all enter the manuscript.
The changed-decoder state-based QBP result remains separate; the necessity
claim here concerns the fixed inverse-frame decoder.

The established scientific package has proof homes, explicit resource/error
contracts, attributed comparisons, and scoped executable evidence. Readiness
requires their consistency and passing verification and presentation gates;
it does not require closing the open gaps or supplying a general native
emitter. Follow the [research stopping rules](../WORKSPACE.md#research-decision-and-stopping-rules)
for any new construction pass. Otherwise retain the unresolved questions
with their actual bounds and stop adding exploratory fixtures.

When final writing resumes, assemble Results A–D and their retained
corollaries into one argument. External technical feedback and submission
follow the complete draft. Update `CITATION.cff` when a manuscript identifier
exists.
