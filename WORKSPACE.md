# Continuing research workspace

This is the entry point when a previous conversation or execution
workspace is unavailable. Proofs and research decisions live in the
repository; conversation summaries are retrieval aids.

The 2026-10-01 complex-state pass starts from verified main
`5a3669ec21ccbd13a53aa69a960c6d66aef90a79`. Check the current branch
and later commits before continuing. The state-only compiler prepares the
real Hopf state with two clean flags and one precision charge. A new
coarse-frame measurement and corrected classical decoder remove the
leaf-reference protocol's angle-dependent shot penalty, including at
singular angles. A complete small native example now verifies the changed
protocol with an active arbitrary dirty helper and complex observables.
The complex extension fixes a common state phase, compiles its coarse
prefix phase tables with exact dirty return, and estimates both magnitude
and leaf-phase gradients with the same two-flag allocation.
The complete-frame endpoint remains open; the task-specific result does
not require its resolution.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized for this repository.
Keep manuscript writing and release work outside this research pass.
Use small analytic examples and finite checks; no large simulations,
QRAM, resets, supplied catalysts, or hidden initialized work.

The user has emphasized the Hopf QBP task. The retained frame problem
concerns its prescribed tree completion; arbitrary unitary synthesis is
only a comparison. Supporting all allowed observables with the existing decoder
still requires the designated marker directions. The new reference-state
route explicitly changes the decoder while retaining all allowed observables
and real raw gradient coordinates, with its own sample and precision bounds.
The complex coarse-frame extension includes the separate leaf-phase stream.

1. Read the [complex coarse compiler](docs/COMPLEX_COARSE_COMPILER.md) and
   [complete complex gradient protocol](docs/COMPLEX_COARSE_QBP.md), then
   the [real coarse-frame protocol](docs/COARSE_FRAME_QBP.md) and
   the [state-only construction](docs/STATE_ONLY_COMPILER.md), including its
   common-coarse reference corollary. The earlier
   [leaf-reference protocol](docs/REFERENCE_STATE_QBP.md) remains a separate
   option under its explicit sampling tradeoff.
   The [native integration example](docs/NATIVE_COARSE_QBP.md) specifies
   what is implemented and what remains an analytic general construction.
2. Read the [task-specific assessment](docs/OPEN_PROBLEM.md#a-state-only-route-for-raw-hopf-gradients)
   and [next bounded question](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
3. For the separate complete-frame question, read the
   [grouped compiler](docs/CONDITIONAL_SUFFIX_COMPILER.md), its completed
   canonical audit, and the
   [scattering restriction](docs/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer).
4. Use [verification](docs/VERIFICATION.md), [attribution](docs/SOURCE_MAP.md),
   and [related work](docs/RELATED_WORK.md) before extending a claim.
   The [overview](README.md), [technical narrative](REVIEW.md), and
   [publication scope](manuscript/PUBLICATION_SCOPE.md) retain the frame story.

## What is established

| Result | Capability and boundary |
|---|---|
| Exact compilation | Matching size and CNOT/depth tradeoffs for the prescribed complete frame, for every clean-work budget |
| Fault-tolerant frontier | Matching T-count under the stated sufficient-clean reservation; complete-frame and fixed-parameter QBP guarantees stand |
| Constant-clean compiler | One clean qubit gives $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ at $`b\ge L+n+7`$; T-depth is a separate open problem |
| Standard operator families | Literal diagonals and complete one-target U(2) multiplexors have matched one-clean banked T-counts at their separate dirty thresholds; their use in a growing product still needs a precision ledger |
| Promised updates | Antichain and sparse ancestor-closed changes have $`O(N+L)`$ T-count under their distinct assumptions; generic rounding supplies neither promise |
| Compact residuals and transport | Linear classical Cayley/weighted data, complete coupled boundaries, and cheap source-width transitions; generic coherent conversion and interior programming remain charged |
| Small native programs | Fixed-address quaternion compression, four-mode magic-basis factors, and eight-/sixteen-mode changing-target words with full borrowed-signal return; fixed-size cost savings do not improve the generic bound |
| Hopf scattering | One packed local-rotation step has $`T=O(N+L)`$, $`G=O(NL)`$ and explicit port permutations; its feedback boundary equals the complete frame, but that feedback is not supplied by the step compiler |
| State-only Hopf QBP | At $`L\ge\max\{6,n\}`$, two compiler flags and $`b\ge L+n+7`$ give $`T=O(N+L)`$, $`G=O(NL)`$ for preparation; a common-coarse reference preserves relative branch phase |
| Coarse-frame Hopf QBP | All real angles admit depth-record norm at most five, the original logarithmic depth-block shot order, and $`O(N+L')`$ T-count per execution apart from the observable; histogram reconstruction costs $`O(S+Nn)`$ arithmetic plus preprocessing and bit costs |
| Complex state-based Hopf QBP | Gauge-fixed prefix phase tables supply an exact-return logical coarse C with $`O(N)`$ T gates; the same two flags and dirty threshold give magnitude and phase streams, $`O(N+L')`$ T per execution, and $`O(S+Nn)`$ histogram arithmetic |

For the selected complete real-frame endpoint,

```math
N=2^n,\qquad a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

```math
\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n)),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

The upper bound also uses only one clean qubit. Recent passes have not
narrowed this gap or established that its resolution is close. The
established publication results do not depend on closing it.

Ordinary fixed-accuracy QBP does not require $`L=N`$. For fixed observable
coefficient scale and raw-gradient accuracy, the
[per-execution choices](docs/QBP_APPROXIMATION.md#8-per-execution-and-complete-gradient-t-costs)
already give optimal-order $`\Theta(\sqrt N)`$ per-frame T-count with
zero compiler clean qubits and $`\Theta(\sqrt N)`$ dirty work. The
protocol still reserves its interference qubit and observable work, and
its total gradient cost is not claimed optimal.

The best compiler already replaces n ordinary
layers by $`R=O(\ell_*(n))`$ unequal groups while keeping total table
work linear. A better constant in a layerwise word does not remove that
remaining group factor. The canonical audit now transfers the paired-source
real-Y circuit to that grouped interface, including complex phases, distinct rejection
flags, and the private-suffix reflection. Its common-width precision term
is still $`(120R+4)q`$, before separately charged queries and controls.

## What the completed attempts settle

The canonical coefficient rotation uses the existing scalar flag sigma,
with the native synthesis signal borrowed outside the initialized-work
reflection. Direction-controlled Z gates select its actual inverse, so a
single coefficient program handles forward and reverse atoms. The old
one-tail source has the same symmetry; this consolidation improves both
baselines equally.

Three inner rotation appearances implement one amplified group. The new
full-operator estimate fits its existing error budget at
$`q_g=m_g+4`$. The core and borrowed signal add six dirty slots before
fixed helpers, absorbed by the established group slack; the external
clean budget stays two. Source/core return is approximate within the
complete error, while query selectors and helper subroutines return
exactly. A zero coefficient remains a nonidentity canonical rotation.

The resulting word still pays for thirty programmable-mask appearances
per group, compared with six in the consolidated old scalar route.
Common-source cancellation does not remove that precision term. Keep the
proved compatibility result and its finite fixture, and **stop optimizing
this canonical word as an endpoint mechanism**.

The width-only source transition, rejected-branch repair, shared native
word, and canonical group interface have each been audited. None gives a
variable-group precision saving. Their scoped restrictions do not make
the group factor necessary for arbitrary circuits.

The new [Frobenius-coarse specialization](docs/ENDPOINT_TREE_TRANSPORT.md#frobenius-small-residuals-also-fit-the-endpoint-budget)
also rules out lack of small residual norm as the next missing ingredient:
an $`O(N)`$-T coarse frame with exact dirty return can make
$`\|C^\dagger W-I\|_F`$ a fixed small constant. Generic near-identity
synthesis, state preparation plus ordinary Householder reflections, and
counter-only clean-work reduction do not meet the endpoint as supplied.
The [selection audit](docs/OPEN_PROBLEM.md#selection-audit-after-canonical-completion)
records their precise costs and contracts. This comparison selects no
construction and changes no retained frontier.

## Next bounded decision

The selected state-only route has an explicit program: use an actual native
coarse C, classically compute $`C^\dagger\psi`$, load its bounded off-root
coefficients into one addressed SU(2) table, and amplify the normalized
accepted state. The two flags and initialized system define the initial
reflection; every dirty wire is excluded from that reflection and returned
within the joint error. This does not implement the other frame columns.

For local branch probabilities in $`[p_0,1-p_0]`$ with fixed positive
$`p_0`$, the target state itself is a valid reference with uniformly bounded
scores. Arbitrary angles instead use the positive derivative-envelope
reference and a sufficient shot factor $`Z+1\le n+1`$. Singular examples
attain $`Z=n`$; do not treat the factor as constant without a promise.
The protocol reserves its separate interference qubit and observable work.

The new coarse-frame decoder settles this sampling question. Prepare the
actual native $`C|0\rangle`$ as reference, use the charged coarse inverse
after the controlled observable, and spread with Hadamards. Scores from
$`H^{\otimes n}C^\dagger\partial_j\psi`$ have uniformly bounded depth
norm; randomized X/Y measurements handle their complex entries. The
classical coefficients include the actual C, so its coarse discrepancy
does not become bias. A signed histogram and adjoint traversal avoid a
dense Jacobian. The new proof supplies the complete precision, reused-dirty,
shot, quantum-work, and classical-output ledger.

The bounded native integration is complete. Its two-qubit exact target uses
the state compiler's finite-size fallback, and its per-reflection echoes
actively borrow one arbitrary helper. Complete initialized isometries and
dirty-input score operators verify phases and gradient means. A reusable
integer-histogram decoder applies the actual logical tree blocks and a
reverse traversal, with explicitly floating-point final contractions.
The supplied words cost more than the original exact-frame protocol on this
fixture; they do not demonstrate a practical advantage. The general
fine-precision residual-table emitter remains outside this implementation.

The complex coarse-interface audit is complete. For supplied real leaf
phases, subtract their arithmetic mean mu. The standard prefix phase
cascade is determinant one and the borrowed reflection interpreter compiles
an actual logical $`C\approx e^{-i\mu}D_\varphi W_{\mathbb R}`$ with
exact dirty return and $`O(N)`$ T gates. Its actual phase rows can be
nondiagonal and act on every suffix; classical application costs $`O(Nn)`$.
The same residual table prepares the gauged target or its common coarse
reference. The corrected magnitude decoder removes the gauge at its leaf
weights, and the established direct phase-Y stream supplies the other N
coordinates. The physical energy gradients do not depend on the gauge.

The next bounded implementation task is to join these complex phase tables
to the existing small elementary native fixture and verify both streams
against arbitrary dirty input with literal branch phases. The present new
checks already use native one-qubit row words, but the complex prefix
selection and fine residual table are still analytic constructions.
Require complete gate counts and same-observable comparison before claiming
an implemented complex protocol. The general fine-precision emitter and
the fine complete-frame endpoint remain separate.

For the separate frame problem, globally combined coefficients with additive
O(N) tables and fully charged basis changes remain eligible. The unchanged
scattering coin and canonical shared-source word are closed as the stated
amortization candidates. Require a new algebraic template and symbolic
precision recurrence before further fixtures. Preserve literal phases,
the exact dirty allocation, arbitrary dirty/reference return, and the
complete prescribed frame. See the
[retained stopping rule](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).

## Restore and verify

Run from the repository root with `requirements.txt` installed:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

The local suite includes elementary native integration, executable integer
histogram checks, and bounded complex gauge/gradient checks alongside the
retained ideal decoder and state-amplification checks. Merge gates include Python
3.11/3.13 CI, all four exact-receipt suites, and rendered presentation.
The retained scattering checks cover complete small frames and literal port
permutations on matrices of dimension 16 and 32; their coins remain ideal.
The retained q=2 native table has dimension 256; its outer group propagates
64 logical/borrowed input columns through its 2048-dimensional space without
a dense group matrix. These are finite interface diagnostics, not
fine-precision or asymptotic certification.
Internal algebra and resource reviews are not external peer review.

This pass supplies the complex task-specific proof and its histogram
implementation. Complete native integration currently covers the real
target fixture; the frame endpoint is unchanged.
No human action is required to continue. Focused feedback from a fault-tolerant synthesis
specialist could help assess an explicit new factorization; it is not an
unstated dependency or a claim of external validation. Keep conclusions in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md) and proofs in their existing
chapters so continuation does not depend on an old chat.
