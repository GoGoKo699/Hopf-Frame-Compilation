# Continuing research workspace

This is the entry point when a previous conversation or execution
workspace is unavailable. Proofs and research decisions live in the
repository; conversation summaries are retrieval aids.

The 2026-10-01 canonical-group pass starts from main commit
`da0980a5ca10f8f9c47af18ba39a13418741f36b`. Check the current branch
and later commits before continuing. The canonical scalar substitution
now fits the actual grouped compiler, with a complete error/workspace
ledger. Its precision cost still grows with the group count, so this
candidate is closed for endpoint amortization. No next construction is
selected; mechanism selection is the remaining research action.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized for this repository.
Keep manuscript writing and release work outside this research pass.
Use small analytic examples and finite checks; no large simulations,
QRAM, resets, supplied catalysts, or hidden initialized work.

1. Read the [current assessment](docs/OPEN_PROBLEM.md#current-assessment-what-the-results-establish)
   and [revised next task](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
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
| Promised updates | Antichain and sparse ancestor-closed changes have $`O(N+L)`$ T-count under their distinct assumptions; generic rounding supplies neither promise |
| Compact residuals and transport | Linear classical Cayley/weighted data, complete coupled boundaries, and cheap source-width transitions; generic coherent conversion and interior programming remain charged |
| Small native programs | Fixed-address quaternion compression, four-mode magic-basis factors, and eight-/sixteen-mode changing-target words with full borrowed-signal return; fixed-size cost savings do not improve the generic bound |

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

The latest shared word reduces declared source appearances from 135 to
87 at eight modes and from 180 to 114 at sixteen. Hoisting the loader
and commuting the fixed mask through its tail leaves both comparison
words the same leading precision term $`(40g+4)q`$, for g ordinary
stages at source width q. The remaining savings concern fixed scalar
work. Its precision charge still grows with the number of stages.

This matters strategically: the best compiler already replaces n ordinary
layers by $`R=O(\ell_*(n))`$ unequal groups while keeping total table
work linear. A better constant in a layerwise word does not remove that
remaining group factor. The canonical audit now transfers the paired-source
real-Y circuit to that grouped interface, including complex phases, distinct rejection
flags, and the private-suffix reflection. Its common-width precision term
is still $`(120R+4)q`$, before separately charged queries and controls.

## Completed compatibility test and next decision

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

No variable-group fusion rule has qualified. The next task is to select a
materially different native mechanism for the coupled residual or retained
precision. A focused primary-source comparison may help, but any imported
result must be checked for its initialization, complete-unitary error,
and arbitrary dirty-work return before adaptation.

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

The baseline has 285 passing tests, Python 3.11/3.13 CI, all four
exact-receipt suites, and rendered presentation. This pass adds three
canonical-group tests, making 288. The new q=2 native table has dimension
256; the outer group propagates 64 logical/borrowed input columns through
its 2048-dimensional space without a dense group matrix. These are finite
interface diagnostics, not fine-precision or asymptotic certification.
Internal algebra and resource reviews are not external peer review.

The canonical pass adds the group compatibility proof and its focused
finite checks. Keep detailed conclusions in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md) and proofs in their existing
chapters so continuation does not depend on an old chat.
