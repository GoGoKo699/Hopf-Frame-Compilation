# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

Baseline: [PR #255](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/255),
`40cd701d47b6e2e0144bbc5c79184b0a3607a63b`, verified 26 October 2026.
Its reviewed tree passed 1,864 tests, the full local gate and all pre-merge CI.
Fetch current main before continuing; this baseline is not a pin.
Exact commit/tree receipts belong in the PR record.

The latest work sharpens complete-source symmetry, count and recovery constraints.
The smaller tree source improves component count/width, not the depth frontier.
Start with these four routes:

1. [Resource frontier](docs/OPEN_PROBLEM.md): strongest general theorems and exact domains.
2. [Clifford-interleaved lookup](research/depth/CLIFFORD_INTERLEAVED_LOOKUP.md): two-block and later-key-mixing boundaries; [central transport](research/depth/COHERENT_CENTRAL_TRANSPORT.md) isolates a partial character basis, while the [dirty Choi source](research/depth/DIRTY_CHOI_OBSERVABLE_SOURCE.md#one-bank-polarization) gives a smaller one-bank extraction from the residual probe.
3. [Clifford projector Grams](research/depth/CLIFFORD_PROJECTOR_GRAM.md#interleaved-clifford-gram-rigidity): the accepted-source boundary; [paid recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) supplies the separate return interface.
4. [Route decisions](research/ROUTE_DECISIONS.md): check the precise failed mechanism and its reopening condition before choosing another word.

<a id="current-depth-frontier"></a>

## Destination and strongest established statements

For $`N=2^n`$, integers $`n\ge1,L\ge6`$, the strongest unconditional
research upper theorem remains [R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md):

```math
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL).
```

One clean qubit suffices at real width $`b\ge17(L+n+7)`$ or
complex-magnitude width $`b\ge17(L+n+8)`$, with complete initialized-isometry
error at most $`2^{-L}`$, including arbitrary dirty/reference inputs and work return.
[Selected Result D remains R79](manuscript/PUBLICATION_SCOPE.md).
At fixed accuracy and square-root dirty width, logarithmic upper versus
constant lower native T-depth remains open. The [frontier](docs/OPEN_PROBLEM.md)
owns all other general bounds, thresholds and special-family qualifications.

## Main unfinished route and exact missing operation

For p address bits and $`m=2^p`$, construct the complete dirty indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle
```

at every fixed requested tolerance tau, with zero clean inputs,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`D_T=o(\log(p+2))`$. Every address, output, work and reference input is
arbitrary. Intermediate register motion is permitted; the complete operator
error includes every return error.

[Polynomial-overhead fusion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
also permits a smaller complete source with fixed polynomial count/width overhead.
[Bounded-clean promotion](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
permits bounded-clean discovery with paid leakage and off-code cancellation.
Neither supplies the missing shallow source. The [R98 transfer](research/depth/QUERY_DEPTH_FRAME_TRANSFER.md)
and its [constant-accuracy extension](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
retain their own log-log contribution: constant source depth alone does not
settle optimal frame depth or give a new all-precision theorem.

### The two useful source interfaces

| Interface | What remains to be constructed |
| --- | --- |
| [Complete Pauli derivative](research/depth/BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) | A native V with the specified full-space derivative, all helpers included. Two actual source calls give the indicator. The canonical diagonal phase is resource-equivalent to the indicator; its truth table is not a cheaper primitive. |
| [Accepted projector Gram](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) | A normalized native source U whose accepted block K has $`K^\dagger K\approx\Pi_{\rm eq}/2`$ on every data/dirty input. Accepted fibers and intermediate address motion may be arbitrary; selecting the correct complete lookup sector is unpaid. |

[One-clean recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
and [joint two-clean recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout)
pay actual inverses, clean leakage and dirty promotion. Their proofs own the
required Gram accuracy and full error/call/width recurrences; use those constants
when comparing a candidate with the complete indicator and frame contracts.

### Latest substantive outcomes and resumption question

The [Clifford-involution source](research/depth/CLIFFORD_INVOLUTION_SOURCE.md)
realizes every traceless Hermitian Clifford as the image of a specified Pauli,
with at most four T bands, $`3w-2`$ T gates, $`O(w^2)`$ Cliffords and no
helpers on w arbitrary wires. Literal phases and actual inverses are paid.
This replaces the first third-level block in a three-block probe source.
For $`A=C_3D_3C_2D_2`$, exact completion is equivalent to
$`A^\dagger Q A`$ being Clifford. Two eight-band cubic-sign outer phases
would give a twenty-band source; the outer factors remain unconstructed.
Arbitrary late address mixing is allowed in this search.
The [joint address/helper trial](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md)
creates the quartic sign but needs a paid carrier repair; no final Clifford alone suffices.
The [complete Clifford-symmetry normal form](research/depth/LOOKUP_CLIFFORD_SYMMETRIES.md)
forces address-algebra preservation at p at least three when a Clifford's
lookup conjugate lies strictly within $`1/\sqrt2`$ of a Clifford; the two-bit exception remains.

The [two-block theorem](research/depth/CLIFFORD_INTERLEAVED_LOOKUP.md) retains
its $`\sin(\pi/8)`$ floor with arbitrary Clifford address/work mixing.
Its [three-block extension](research/depth/CLIFFORD_INTERLEAVED_LOOKUP.md#three-block-late-key-mixing)
has a one-quarter floor above label degree six when the later Clifford walls
preserve the key algebra. Neither is a growing native-depth lower bound.

The [reduced central gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
needs only $`q=\binom p2`$ arbitrary pair wires and the m-bit carrier.
Writing $`U_F|s,c\rangle=|s,F_sc\rangle`$, its complete echo sends
$`|s,c\rangle`$ to $`|s+K,F_Kc\rangle`$. A four-band root shear completes
the probe with $`6p+(p\bmod2)`$ T gates and no helpers; the full linear
factor is unnecessary. Each pair factor has a helper-free four-band compiler, but their overlapping carrier gates have no shallow joint schedule.
The general carrier-linear gauge has an exact root-row criterion and a
strictly larger family than this representation. Its low-row degree bound
is a gauge-class constraint, not a native-depth lower bound.

The [ordered-matching source](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md)
replaces each Pfaffian sum by one matching. Its [minimal-tree specialization](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source)
needs $`6m-5`$ T gates, at most $`4p`$ bands and no private work. Its complete
forward XOR query takes eight bands, but does not supply the in-place inverse.
The retained indicator already has logarithmic depth; faster joint synthesis remains open.
The [complete-probe support bound](research/depth/NATIVE_INDICATOR_SHARING.md#complete-probe-support-count)
forces linear source T count even with arbitrary register motion. It makes the
tree count optimal in order, but cannot force growing depth at the allowed width.

The earlier [weak pair gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#weak-pair-gauge-source)
retains its distinct full-space identity, omitted affine correction, and
failed quadratic truncation and partial-Hadamard dressings. The
[HRT compiler](research/depth/COHERENT_FEEDBACK_SOURCES.md) pays individual
encoded Hadamards in sixteen bands and intersecting pair stars in twenty-four;
disjoint-pair and feedback witnesses remain in that proof home.
The [dirty cubic compiler](research/depth/DIRTY_CUBIC_PHASE_COMPILER.md)
and [half-band exchange](research/depth/COEFFICIENT_SHEAR_NATIVE_FUSION.md)
remain useful paid components. None completes the general indicator.

The [joint precision primitive](research/depth/JOINT_PRECISION_ROTATION.md)
retains the one-clean chart and projector-amplitude recovery. A Clifford image
of the actual accepted signal forces unit acceptance below error $`2\sin(\pi/16)`$.
Paid balanced branches have a direct alternating-word bypass; the paid 3:1 pilot
has a direct branch within $`3\sqrt\epsilon`$ of any recovered unitary target.
Both comparisons require separately paid branches. General precision depth remains linear.

**Next question:** choose the two outer conjugations so the complete peeled
lookup observable is Clifford, with the simultaneous polynomial discovery
budgets and arbitrary address mixing. The reduced exterior and ordered-matching
gauges are concrete alternative targets for joint synthesis. A finite-size identity must extend to a
uniform all-size native rule; restoring logarithmically many serial grade
computations does not improve the frontier. In precision, seek joint recovery
or another accepted source beyond the separately paid two-path pilots.
Any success needs a fresh complete-circuit ledger, including the
transfer's log-log term and every remaining precision-sized word. These are
sufficient routes, not necessary representations of an unrestricted compiler.

## Stopped routes, reopening conditions and superseded instructions

Use [route decisions](research/ROUTE_DECISIONS.md) for the current assessment
of each mechanism. Detailed proofs and unique failed attempts remain in their
canonical homes, indexed by the [research catalogue](research/README.md).
Old instructions to add another mixer, improve a scalar filter, or pay the
rank-one phase do not override the complete-action screens and reductions above.
A scoped failure does not retire the unrestricted destination.

<a id="research-decision-and-stopping-rules"></a>
<a id="continuing-toward-the-general-answer"></a>
<a id="read-by-claim"></a>

## Working contract and reading map

Prefer analytic identities and a few exact fixtures. Preserve literal phases,
every logical column, coherent address and arbitrary dirty/reference inputs,
actual native adjoints, and all clean leakage and work return. No QRAM, free
loading or precision states, unpriced clean history, resets or postselection.
A suffix is clean only on its active sector; inactive action needs its own proof.

Separate general theorems, conditional interfaces, useful components and
representation limits. Change mechanisms when complete synthesis restores the
old cost. Internal audits and finite checks do not replace the analytic proof
or external technical review.

| Need | Authoritative home |
| --- | --- |
| Selected A–D claims and allocations | [Publication scope](manuscript/PUBLICATION_SCOPE.md) |
| Core argument and dependencies | [Review](REVIEW.md), [documentation map](docs/README.md) |
| Strongest general resource statements | [Resource frontier](docs/OPEN_PROBLEM.md) |
| Mechanism decisions and proof locations | [Route decisions](research/ROUTE_DECISIONS.md), [research catalogue](research/README.md) |
| Attribution and executable evidence | [Source catalogue](docs/reference/SOURCE_CATALOGUE.md), [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md) |

<a id="refurnishing-policy"></a>

Each claim has one proof home. Keep navigation as short routes into those homes;
do not accumulate session summaries. Preserve unique attempts, stable identifiers,
inbound anchors, scientific tests and receipt hashes when reorganizing material.

## Restore and verify

Check worktree status and remote, fetch main, and inspect changes since the
baseline. Preserve unrelated work. [SYNC](SYNC.md) owns upstream versions.
For a proof move, compare the moved mathematics and audit old and new links.
Run the complete local gate on the final reviewed tree:

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

Use the [presentation guide](assets/README.md#rendering-checks) for changed
mathematics and reader-facing pages. Complete focused visual checks before the
memory-heavy numerical suite; use full presentation CI for the complete corpus.
Record exact-tree receipts in the review record and merge after required checks pass.
A fixture is not an asymptotic proof or a scalable native emitter.
