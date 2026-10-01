# Continuing research workspace

This is the entry point when a previous conversation or execution
workspace is unavailable. Proofs and research decisions live in the
repository; conversation summaries are retrieval aids.

The 2026-10-01 Hopf-scattering pass starts from the verified main checkpoint
`5bff4239a3547e555824c8cf31dfda11c97a008d`. Check the current branch
and later commits before continuing. The complete Hopf frame now has an
explicit tree-scattering identity with one cheaply compiled local-rotation
table. Turning that table into its boundary action remains necessary;
a constant number of unchanged-table queries cannot do so uniformly at
fine precision. This is a query-model restriction, not a T-count lower
bound. No endpoint construction is selected. The next selection question
is a native program with coefficients that already combine tree levels.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized for this repository.
Keep manuscript writing and release work outside this research pass.
Use small analytic examples and finite checks; no large simulations,
QRAM, resets, supplied catalysts, or hidden initialized work.

The user has emphasized the Hopf QBP task. The target is its prescribed
tree frame, not arbitrary unitary synthesis. General theorems are only
comparisons. Supporting all allowed observables with the existing decoder
still requires the designated marker directions; a restricted observable
or a different decoder must be identified as a different task.

1. Read the [current assessment](docs/OPEN_PROBLEM.md#current-assessment-what-the-results-establish)
   and [Hopf-specific scattering audit](docs/OPEN_PROBLEM.md#the-hopf-specific-scattering-test),
   then its [next task](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
2. Read [conditional-suffix Sections 3--7](docs/CONDITIONAL_SUFFIX_COMPILER.md#3-a-small-coefficient-table-and-a-streamed-coarse-circuit)
   and its one-clean extension, followed by the completed
   [canonical group audit](docs/CONDITIONAL_SUFFIX_COMPILER.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge).
   These give the actual scalar/atom interface, private work, error budget,
   and group sums.
3. Read the [shared source body](docs/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
   and [source-width analysis](docs/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation).
   Their different source families and scope must remain explicit.
4. Use [verification](docs/VERIFICATION.md), [attribution](docs/SOURCE_MAP.md),
   and [related work](docs/RELATED_WORK.md) before extending a claim.
   The [overview](README.md), [technical narrative](REVIEW.md), and
   [publication scope](manuscript/PUBLICATION_SCOPE.md) retain the main story.

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

The [scattering test](docs/OPEN_PROBLEM.md#the-hopf-specific-scattering-test)
implements the original node rotations as one addressed table. Its exact
feedback identity is a concrete Hopf representation, not a free inverse.
The path amplitude shows that a bounded-query conversion of this unchanged
table is insufficient even with arbitrary fixed interleaving operations.
Do not schedule another walk-power, fixed-routing, or feedback-query test
with that same table as the proposed endpoint mechanism.

Seek instead an explicit whole-frame or whole-residual factorization whose
programmed coefficients already combine levels, with additive table size
$`O(N)`$ and every change of basis charged. A bounded number of complete
diagonal or multiplexor primitives remains one sufficient representation;
an expanded native word with many calls but one precision charge is also
allowed. The query restriction does not price that word. No such identity
is currently known here; compact classical data or a small matrix norm is
not that identity.

The next deliverable is one actual algebraic template with its symbolic
cost, or a concise no-candidate finding. Do not schedule more fixtures
without that template. A new complete rule for carrying precision also
remains eligible. A stronger unrestricted lower-bound approach remains
open, but needs an invariant beyond the current interface restrictions.

Require an explicit operation or encoding rule and a symbolic precision
recurrence before starting another construction pass. For the retained
grouped route, a sufficient target is

```math
T_{\rm joint}\le c_RL+O\!\left(\sum_gQ_g+RP(n)\right),
```

with $`c_R`$ bounded independently of R and P a fixed polynomial
independent of L. Table and coarse-program sizes must combine additively.
Other whole-frame constructions remain allowed. An unpriced resolvent,
Hamiltonian, coordinate transform, or target-dependent decoder is not a
native construction.

Count every source, native control, reflection, and inverse. Allocate the
full error before checking the exact $`b=N+n+7`$ workspace threshold;
constant extra accuracy bits are not automatically free dirty space.
Final clean return, arbitrary dirty/reference inputs, and literal phase
must satisfy the retained complete-isometry contract. The T-only endpoint
allows a larger fully charged Clifford count than the stronger joint goal.

Use small examples only after a concrete rule is specified. For a grouped
proposal, test unequal groups and require the same rule to survive a
third. If no mechanism meets that selection standard, record that no lead
is selected rather than returning to generic mask sweeps or the completed
boundary audits. The [research decision](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule)
retains the live-work contract and stopping conditions.

## Restore and verify

Run from the repository root with `requirements.txt` installed:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

The local suite passes all 292 tests. Merge gates also include Python
3.11/3.13 CI, all four exact-receipt suites, and rendered presentation.
Four new scattering checks cover complete small frames and literal port
permutations on matrices of dimension 16 and 32; their coins remain ideal.
The retained q=2 native table has dimension 256; its outer group propagates
64 logical/borrowed input columns through its 2048-dimensional space without
a dense group matrix. These are finite interface diagnostics, not
fine-precision or asymptotic certification.
Internal algebra and resource reviews are not external peer review.

This pass adds a Hopf-specific scattering identity, a one-step native ledger,
and a scoped query restriction, not a new endpoint bound. No human action is
required to continue. Focused feedback from a fault-tolerant synthesis
specialist could help assess an explicit new factorization; it is not an
unstated dependency or a claim of external validation. Keep conclusions in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md) and proofs in their existing
chapters so continuation does not depend on an old chat.
