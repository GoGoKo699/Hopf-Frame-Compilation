# Manuscript guide

[Publication scope](PUBLICATION_SCOPE.md) · [Read the argument](../REVIEW.md) · [Verification](../docs/VERIFICATION.md)

Working title: **Exact and Fault-Tolerant Compilation of Hopf Differential Frames**.

The [publication scope](PUBLICATION_SCOPE.md) is the canonical ledger of
Results A–D, retained corollaries, literal workspace reservations, and error
contracts. The [bounded internal audit](../docs/CORE_CLAIM_AUDIT.md) found no
unresolved claim-level blocker. The selected scientific scope is frozen;
final manuscript writing remains on hold.

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
5. Result C: one-clean rotations and conditional-suffix grouping, followed by
   banked, diagonal, multiplexor, complex-magnitude, and zero-clean corollaries.
   Retain the earlier two-clean constructions at their smaller dirty budgets.
6. Result D: simultaneous T-count and T-depth, its two-clean reservation,
   explicit precision/workspace matching intervals, and applicable fallbacks.
7. Exact and approximate fixed-parameter QBP consequences, complete complex
   gradients, and charged quantum executions and classical reconstruction.
8. The two open resource questions and the limits of the circuit model.

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

The [research archive](../research/README.md) preserves attempts to close the
two gaps and their scoped findings. The completed
[state-based QBP supplement](../supplements/state_based_qbp/README.md) uses a
changed decoder and remains separate from the selected paper. Neither branch
adds a default requirement for manuscript readiness.

Finite checks support fragile identities and resource ledgers; they do not
replace asymptotic proofs, a scalable native emitter, or external peer review.
Maintain consistent claims, attribution, and verification. Repair discovered
defects; otherwise follow the [research stopping rules](../WORKSPACE.md#research-decision-and-stopping-rules).
Final writing begins only on request. External technical feedback and submission
preparation follow the complete draft; update `CITATION.cff` when a manuscript
identifier exists.
