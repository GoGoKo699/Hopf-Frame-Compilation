# Manuscript guide

[Publication scope](PUBLICATION_SCOPE.md) · [Read the argument](../REVIEW.md) · [Verification](../docs/VERIFICATION.md)

Working title: **Exact and Fault-Tolerant Compilation of Hopf Differential Frames**.

The [publication scope](PUBLICATION_SCOPE.md) is the canonical ledger of
Results A–D and retained corollaries, with literal workspace reservations
and error contracts. Result C is the optimal one-clean count theorem; the
sharper two-clean complex endpoint retains its separate allocation. The source
and verification maps identify the premises and executable evidence for each claim.

## The argument to assemble

The paper asks what it costs to preserve a prescribed Hopf completion rather
than only its first state column. Exact entangling costs and fault-tolerant
precision costs describe the same operator under two resource models. QBP
explains why the additional columns matter to its fixed inverse-frame decoder.
The one-clean diagonal and one-target multiplexor corollaries also establish
compiler capabilities beyond the motivating frame family. Result A covers
every integer clean-workspace budget $`m\ge0`$.

The principal comparisons are Yuan–Zhang's exact state-preparation tradeoffs,
Low–Kliuchnikov–Schaeffer's dirty-workspace tradeoffs, and
Gosset–Kothari–Wu's optimal T-count. Use the
[related-work comparison](../docs/RELATED_WORK.md) and
[source map](../docs/SOURCE_MAP.md) for their exact statements, newer results,
and attribution. Explain the complete-operator guarantees through their proof
mechanisms and input/workspace contracts.

## Paper order

1. Prescribed completion versus one-column preparation, using the two-qubit
   readout example and the fixed-decoder necessity theorem.
2. The Hopf operator and common complete-input contract: literal phases,
   arbitrary dirty/reference inputs, actual inverses, and workspace return.
3. Result A: exact all-workspace size, CNOT count, and depth, with matching
   lower bounds and the phase-dressed complex magnitude extension.
4. Result B: matching sufficient-clean T-count through residual dictionaries,
   one shared geometric source, and retained failure history.
5. Result C: the one-clean count frontier. Explain growing sources with
   separate precision and dimension scales, borrowed conjugate-signal recovery,
   exact dirty-bank queries, the layerwise regime splice, and the Clifford
   ledger. Compose literal leaf phases under the shared-workspace error contract.
   Retain diagonal, multiplexor, zero-clean, and sharper two-clean endpoint
   corollaries at their distinct reservations.
6. Result D: simultaneous T-count and T-depth, its two-clean reservation,
   explicit precision/workspace matching intervals, and applicable fallbacks.
7. Exact and approximate fixed-parameter QBP consequences, complete complex
   gradients, charged executions and reconstruction, and the remaining
   large-width depth question and limits of the circuit model.

Keep the formula ledger in the scope and the technical proofs in their primary
chapters. Appendices supply source preparation, echoes, decoders, routing,
dirty queries, residual composition, amplification, error/resource sums,
Euler/diagonal details, lower bounds, and finite classical preprocessing.

## Boundaries and stopping point

Distinguish exact dirty return from return included in the initialized-isometry
error. A logical zero suffix supplies clean work only on its proved active
sector. T-depth allows arbitrary Clifford interlayers and differs from total
depth. Only D claims simultaneous T-count/T-depth guarantees; its real
frame theorem does not inherit the complex extensions of A–C automatically.

Result C combines the [one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
and [banked proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier).
Its matching statements retain the exact-base precision condition and banked
width threshold. The [sharp two-clean endpoint](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
has a distinct complex dirty reservation. These count results include no
new depth or general native-emitter guarantee.
The [research catalogue](../research/README.md) preserves their dependencies,
the remaining large-width depth question, and earlier attempts with their
scoped findings. The completed
[state-based QBP supplement](../supplements/state_based_qbp/README.md) uses a
changed decoder and remains separate from the selected paper. Closing the
depth gap or extending that supplement is not a manuscript-readiness requirement.

Finite checks support fragile identities and resource ledgers; they do not
replace asymptotic proofs, a scalable native emitter, or external peer review.
Maintain consistent claims, attribution, and verification. Repair discovered
defects; otherwise follow the [research stopping rules](../WORKSPACE.md#research-decision-and-stopping-rules).
Final writing begins only on request. External technical feedback and submission
preparation follow the complete draft; update `CITATION.cff` when a manuscript
identifier exists.
