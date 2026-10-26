# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

Verified baseline: [PR #257](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/257),
`5a357c0c4cdb6ebc6f17b1bb3fdfe3aa3afc4ce9`, 26 October 2026.
Its reviewed tree passed 1,886 tests, the full local gate and qualifying pre-merge CI.
Fetch current main before continuing; this baseline is not a pin.
Exact commit/tree receipts belong in the PR record.

Read the [resource frontier](docs/OPEN_PROBLEM.md) for general theorems,
then the [route decisions](research/ROUTE_DECISIONS.md) for the chosen mechanism.
Use the [research catalogue](research/README.md) to find its proof dependencies.
The latest nonmonomial-source and readout results are integrated; none improves R96.

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
owns the other bounds, matching ranges and special-family qualifications.

## Main unfinished route and exact missing operation

For p address bits and $`m=2^p`$, construct the complete dirty indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle
```

at every fixed requested tolerance tau, with zero clean inputs,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`D_T=o(\log(p+2))`$. Every address, output, work and reference input is
arbitrary. Intermediate register motion is permitted; full operator error
includes every return error. This is a sufficient route, not a necessary
representation of an unrestricted frame compiler.

[Polynomial-overhead fusion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
allows discovery at fixed polynomial overhead on a smaller address;
[bounded-clean promotion](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
pays leakage and off-code cancellation. The
[R98 transfer](research/depth/QUERY_DEPTH_FRAME_TRANSFER.md) and its
[constant-accuracy extension](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
retain a log-log contribution. Constant source depth alone does not settle
optimal frame depth or give a new all-precision theorem.

<a id="the-two-useful-source-interfaces"></a>

### Choose a complete source interface

| Interface | First unpaid operation |
| --- | --- |
| [Complete Pauli derivative](research/depth/BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) | Synthesize a native source with the required full-space derivative. Two actual source/adjoint calls extract the indicator; its canonical diagonal phase is an equivalent target. |
| [Accepted projector Gram](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) | Synthesize a normalized source with $`K^\dagger K\approx\Pi_{\rm eq}/2`$ on every data/dirty input. Accepted fibers and intermediate address motion may be arbitrary. |

Use the proof home's [one-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
or [joint two-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout)
recovery ledger at the actual requested Gram accuracy. The wrappers are paid;
they do not construct the source.

### Latest substantive outcomes and resumption question

The [four-band Clifford-involution compiler](research/depth/CLIFFORD_INVOLUTION_SOURCE.md)
removes the first source block as an independent search variable. Its
[exact Galois screen](research/depth/CLIFFORD_INVOLUTION_SOURCE.md#galois-outer-symmetry-screen)
constrains the combined outer word; fourth-root phases require even Clifford
displacement rank. This arithmetic screen supplies no approximate error floor.
The [Hadamard-pair trial](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md#nonmonomial-hadamard-pair-screen)
with its specified first coupling has error at least 1/16 for every later
Clifford wall and cubic-sign batch. Weighted phases escape that proof.
The [Fourier-swap source](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md#fourier-swap-fixed-coefficient-face)
completes the fixed-top-coefficient face, but reduces to a quadratic finite
difference. The direct phase costs 28 T gates in eight bands without helpers.
Its common address-independent readout has a separate 1/2 error floor when
that top coefficient becomes arbitrary.

**Next question:** choose and pay the two outer conjugations
$`A=C_3D_3C_2D_2`$ so the complete peeled lookup observable
$`A^\dagger Q A`$ is Clifford, using genuinely later address participation
and retaining the simultaneous polynomial discovery budgets. A successful
fixed-coefficient face must also work when that coefficient is arbitrary.
The [Clifford-symmetry normal form](research/depth/LOOKUP_CLIFFORD_SYMMETRIES.md)
screens qualifying Clifford-to-Clifford images; it does not force arbitrary
three-block walls to preserve the address algebra.

Two concrete alternatives are joint synthesis of the
[reduced central gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
or the [minimal-tree probe](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source).
The latter has optimal-order count but linear address depth; its eight-band
forward XOR query does not supply the in-place inverse. The
[syndrome bound](research/depth/FIXED_READOUT_SYNDROME_BOUND.md)
also excludes an arbitrary fixed final readout after low-degree Pauli
transport. The readout must participate in the address-dependent computation
or leave that transport grammar. These components do not improve the retained
logarithmic indicator. In the separate
[precision route](research/depth/JOINT_PRECISION_ROTATION.md), joint recovery
beyond the separately paid two-path pilots remains open.

A finite-size identity must extend to a uniform all-size native rule.
Restoring serial grade computations retains the old cost. Any proposed gain
needs a fresh complete-circuit ledger, including transfer, precision words,
actual inverses, leakage and work return.

## Stopped routes, reopening conditions and superseded instructions

The [route decisions](research/ROUTE_DECISIONS.md) own each failed mechanism
and its reopening condition. Detailed proofs and unique counterexamples stay
in their linked homes. Old instructions to add another mixer, improve a scalar
filter, or pay the rank-one phase do not override those complete-action screens.
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
representation limits. Internal audits and finite checks do not replace
analytic proof or external technical review.

| Need | Authoritative home |
| --- | --- |
| Selected A–D claims and allocations | [Publication scope](manuscript/PUBLICATION_SCOPE.md) |
| Core argument and dependencies | [Review](REVIEW.md), [documentation map](docs/README.md) |
| Strongest general resource statements | [Resource frontier](docs/OPEN_PROBLEM.md) |
| Mechanism decisions and proof locations | [Route decisions](research/ROUTE_DECISIONS.md), [research catalogue](research/README.md) |
| Attribution and executable evidence | [Source catalogue](docs/reference/SOURCE_CATALOGUE.md), [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md) |

<a id="refurnishing-policy"></a>

Each claim has one proof home. Entry pages route to those homes; preserve unique
attempts, stable identifiers, inbound anchors, scientific tests and receipt hashes.
Change mechanisms when complete synthesis restores the old cost.

## Restore and verify

Check worktree status and remote, fetch main, and inspect changes since the
baseline. Preserve unrelated work. [SYNC](SYNC.md) owns upstream versions.
For a proof move, compare the mathematics and audit old and new links.
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
