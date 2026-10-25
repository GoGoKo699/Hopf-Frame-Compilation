# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

Baseline: [PR #249](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/249),
`b8eafb48576b29ee17059e6ced0ec3191cdd6b62`, fetched 25 October 2026.
Its reviewed tree passed 1,778 tests, the full local gate, PR CI and post-merge CI.
Fetch current main before continuing; this baseline is not a pin.
Exact commit/tree receipts belong in the PR record.

The latest results pay native phase components and identify limits of source classes.
They supply neither a faster complete indicator nor a stronger general frame bound.
Start with these four routes:

1. [Resource frontier](docs/OPEN_PROBLEM.md): strongest general theorems and exact domains.
2. [Clifford-interleaved lookup](research/depth/CLIFFORD_INTERLEAVED_LOOKUP.md): the two-block boundary and the growing block count with Clifford fibers; [central transport](research/depth/COHERENT_CENTRAL_TRANSPORT.md) isolates a partial character basis, while the [dirty Choi source](research/depth/DIRTY_CHOI_OBSERVABLE_SOURCE.md) reduces a single residual probe to the complete indicator.
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

The [dirty cubic-phase compiler](research/depth/DIRTY_CUBIC_PHASE_COMPILER.md)
pays arbitrary shared-control weighted cubic phases in at most 13 T bands,
with returned arbitrary helpers. The [half-band exchange](research/depth/COEFFICIENT_SHEAR_NATIVE_FUSION.md)
pays its actual non-diagonal correction in four bands and $`3m`$ T gates.
These are reusable components; neither completes lookup.

The [two-block theorem](research/depth/CLIFFORD_INTERLEAVED_LOOKUP.md) retains
its $`\sin(\pi/8)`$ floor with arbitrary Clifford address/work mixing.
A separate sharp screen now requires at least p globally third-level blocks
below derivative error $`\sqrt2`$ when each block holds the address and has
Clifford bank/work fibers. Arbitrary address-algebra-preserving Clifford
interlayers are included. This is a source-block count, not native T-depth;
genuine non-Clifford fibers and address-superposition dynamics remain open.

The [complete dirty Choi calculation](research/depth/DIRTY_CHOI_OBSERVABLE_SOURCE.md)
reduces the standard corrected teleportation word to one quadratic XOR query.
Its remaining HRT probe correction is not a harmless work action. Four calls
to that residual phase, or eight calls to a source transporting one Pauli to
it, give the complete indicator with paid dirty arithmetic and constant
additional depth. Full in-place HRT synthesis is unnecessary for this interface.
The remaining diagonal phase has an explicit cubic bent dual; its Fourier
kernel does not by itself implement that diagonal phase.

[Central transport](research/depth/COHERENT_CENTRAL_TRANSPORT.md) still pays
$`I+N_x`$ in sixteen bands. An added-generator embedding now reduces its
remaining factor to a fixed partial character basis on $`\binom p2`$ commuting
affine Clifford observables, with arbitrary residual fibers and a paid cubic
middle. The basis is not supplied: if it lies in the Clifford hierarchy, its
level is at least $`\lfloor p/2\rfloor+2`$. This does not bound native depth.

**Next question:** can the single residual HRT probe, or the partial character
basis and its cubic middle jointly, be synthesized at sublogarithmic native
depth within the polynomial address overhead? These are sufficient interfaces,
not necessary representations of every direct source. Exploit the unrestricted
probe-commuting fibers; do not reconstruct the full HRT map by default.
A three-cubic-block attempt must use genuinely non-Clifford fixed-address
fibers or leave the preserved address algebra. Standard dirty teleportation,
bank Hadamards alone, and angle tuning of the failed Clifford-zeta frame do
not supply that change. Require an all-size complete operator identity before
extending native schedules. The accepted-Gram alternative still owes complete
sector selection.

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
