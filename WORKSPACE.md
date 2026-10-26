# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

Verified baseline: [PR #260](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/260),
`9257cbe48e30b9b12fef03def2bb9217757d0133`, 26 October 2026.
Its reviewed tree passed 1,944 tests, the full local gate and qualifying pre-merge CI.
Fetch current main before continuing; this baseline is not a pin.
Exact commit/tree receipts belong in the PR record.

Read the [resource frontier](docs/OPEN_PROBLEM.md) for general theorems,
then the [route decisions](research/ROUTE_DECISIONS.md) for the chosen mechanism.
Use the [research catalogue](research/README.md) to find its proof dependencies.
The latest source and decoder screens are integrated; none improves R96.

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

The [interleaved tree probe](research/depth/INTERLEAVED_TREE_PROBE_SOURCE.md)
pays a later address query and a fixed Fourier wall. Pulling lookup through
that word gives a Pauli family with full commutation rank $`m=2^p`$.
This is an explicit complete target algebra, including literal signs and
arbitrary dirty inputs; canonical target coordinates do not pay a source decoder.

The [one-query gaps](research/depth/INTERLEAVED_TREE_PROBE_SOURCE.md#one-query-orbit-gaps)
exclude a single linear XOR query at error below $`\gamma=\sqrt{2-\sqrt2}`$,
even with an arbitrary non-Clifford starting reflection and fixed decoder.
Affine offsets retain a covariance gap $`\gamma/2`$ for $`p\ge2`$.
Two opposite interleaved queries escape these premises. For a fixed
**Clifford** starting reflection, however, the
[faithful phase-lift bound](research/depth/CLIFFORD_ORBIT_QUERY_RANK.md)
requires at least $`\lceil p/8\rceil`$ literal affine queries at error
below $`\gamma/4`$. This is a query-grammar bound, not native-depth optimality.

The [two-address pilot](research/depth/TWO_ADDRESS_INTERLEAVED_PILOT.md)
constructs the whole small family with controlled Hadamards and a paid fixed
Clifford decoder. Its complete native source and actual inverse provide a
benchmark, not an all-size depth improvement.

The [non-Clifford seed-spread bound](research/depth/NONCLIFFORD_SEED_SPREAD.md)
requires an exponentially broad Pauli-label span for bounded affine-query
count, even at constant operator error and with an arbitrary fixed decoder.
A paid two-band tensor seed supplies that spread, but its plain local
inverse decoder has error at least one half for the tree family. A localized
nonlinear patch or more elaborate fixed Clifford spreading does not meet
the requirement.
For a third-level-dressed Clifford seed, the
[decoder condition](research/depth/INTERLEAVED_TREE_PROBE_SOURCE.md#nonclifford-seed-decoder-condition)
also requires every selected target pullback to be non-Clifford below
error $`\sqrt2/4`$. A single weighted-cubic decoder cannot meet that condition.

The [field-shear construction](research/depth/FIELD_SHEAR_INTERTWINER.md)
now has an [explicit all-size recursive decoder](research/depth/RECURSIVE_FIELD_SHEAR_DECODER.md),
in addition to its exact native two-bit word. The recursion uses finite-field
coordinate changes and Clifford Fourier transforms, including every zero
sector. Dense tower queries and multiplication-XOR batches have paid
eight-band implementations. The nonlinear in-place decoder maps remain unpaid.
Preserving the physical Bell zero fiber fails even the weaker lookup probe;
exact probe freedom is a helper-valued gauge commuting with all lookup X's.

**Next question:** jointly compile or replace the recursive decoder,
including its nonlinear coordinate controls, actual inverse and returned work.
Literal level-by-level synthesis already uses at least p minus one T bands;
making each cell faster cannot meet the discovery depth target. Its
read-only affine-query factorization also has a growing degree obstruction.
Use the recursion as an exact target for joint synthesis, not as a shallow
algorithm. A different joint probe gauge or source pencil remains possible.
The alternative remains a jointly designed extensive non-Clifford seed and
nonlocal decoder for the interleaved tree family. Serial address controls
or a generic spectral decomposition retain an unpaid or excessive cost.

The [four-band Clifford front end](research/depth/CLIFFORD_INVOLUTION_SOURCE.md),
[weighted outer-phase route](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md),
[reduced central gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
and [minimal-tree probe](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source)
remain component alternatives under their separate scoped screens. In the separate
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
