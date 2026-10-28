# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

Verified baseline: [PR #265](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/265),
`97ade196941b6d8780556bbfa1585d311beee6f2`, 28 October 2026.
Its reviewed tree passed 1,996 tests in both qualifying CI versions; the local gate exited zero.
Fetch current main before continuing; this baseline is not a pin.
Exact commit/tree receipts belong in the PR record.

Read the [resource frontier](docs/OPEN_PROBLEM.md) for general theorems,
then the [route decisions](research/ROUTE_DECISIONS.md) for the chosen mechanism.
Use the [research catalogue](research/README.md) to find its proof dependencies.
The source, decoder and copy screens are integrated; none improves R96.

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
| [Encoded dirty copy](research/depth/ENCODED_DIRTY_COPY_INTERFACE.md) | Jointly synthesize the encoded bank transition and its complete copy interaction with a physical output bank. Their echo returns the carrier; no physical implementation of its abstract decoder is required. |
| [Accepted projector Gram](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) | Synthesize a normalized source with $`K^\dagger K\approx\Pi_{\rm eq}/2`$ on every data/dirty input. Accepted fibers and intermediate address motion may be arbitrary. |

Use the proof home's [one-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
or [joint two-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout)
recovery ledger at the actual requested Gram accuracy. The wrappers are paid;
they do not construct the source.

### Latest substantive outcomes and resumption question

The [tree-probe algebra](research/depth/INTERLEAVED_TREE_PROBE_SOURCE.md),
[Clifford-orbit rank bound](research/depth/CLIFFORD_ORBIT_QUERY_RANK.md) and
[non-Clifford seed-spread bound](research/depth/NONCLIFFORD_SEED_SPREAD.md)
retain their exact query/seed scopes. Target coordinates and support do not pay decoding.

The [field shear](research/depth/FIELD_SHEAR_INTERTWINER.md) has a paid
eight-band query and an [exact recursive decoder](research/depth/RECURSIVE_FIELD_SHEAR_DECODER.md).
Dense tower queries and multiplication-XOR batches are paid; its nonlinear
in-place maps remain unpaid. Bell-zero-fiber preservation fails even the
weaker lookup probe. For p at least two, a seed diagonal in any query-invariant
stabilizer basis fails, including arbitrary helpers and a non-Clifford decoder.

The [three-band source theorem](research/depth/THREE_BAND_LOOKUP_SOURCE_OBSTRUCTION.md)
excludes every zero-clean native source of depth at most three for p at least
seven below complete observable error $`\sin(\pi/8)/(2\sqrt2)`$.
The proof forces the final two absorbed Pauli bands to preserve the address,
with arbitrary Clifford interlayers, register motion and dirty width.
This is a constant source-depth bound, not a growing frame-depth lower bound.

The [copy-grading criterion](research/depth/DIRTY_COPY_GRADING.md) now
characterizes exact copying of arbitrary unitary families: every word with
nonzero occurrence parity must have zero trace. Spectral multiplicity yields
constant robust word tests. Field-orbit reflection pairs have a common
same-sign sector and cannot approximate an anticommuting pair below
$`\sqrt{2-\sqrt2}`$. Nevertheless a copyable all-size noncommuting family
exists. Its wire-translation query retains linear depth; moving it to the
field query leaves its seed and copy unpaid.

[Native copy components](research/depth/NATIVE_COPY_INTERACTION.md) retain
their one-band component, dual-support bounds and product/coded-family screens.
The [affine Pauli-query screen](research/depth/PAULI_ORBIT_COPY_SCREEN.md)
has a quarter-error floor with Clifford dual probes; it does not cover the
non-Pauli field shear or exclude all one-band field copies.

The [copy-axis bound](research/depth/COPY_AXIS_SPAN.md) requires extensive
physical-X span, compatible with linear count. The
[seed-mixing bound](research/depth/MONOMIAL_QUERY_SEED_OBSTRUCTION.md)
requires strong cross-commutators on every independent address pair;
monomial query-eigenbasis seeds have a half-error floor. The full-Hadamard
field seed passes the mixing test but still has a half-error spectral floor.

The [four-address pilot](research/depth/FIELD_COPY_CLIFFORD_PILOT.md) pays
the seed, copy and inverses: 288 T gates, at most 96 bands and seven dirty
wires. Its cheaper query gives the direct decoder 300 T, 160 bands and three dirty wires.
The [symplectic-orbit theorem](research/depth/CLIFFORD_SEED_COPY_STRUCTURE.md)
gives error at least 1/64 for p at least seven when an involutory Clifford seed's orbit
has commuting symplectic actions, even for an arbitrary copy. The pilot lies
inside this class, so retaining its Pauli-valued commutators cannot scale.

[One-band copy feasibility](research/depth/ONE_BAND_COPY_FEASIBILITY.md)
gives signed Walsh slices, their spectrum and the complete carrier fixed point.
All output probes alone can miss wrong carrier action. The triple-axis joint walls
pass the tensor screen but have zero complete seed space for p at least two,
and full copy error at least $`1/\sqrt6`$ for every unitary seed. This family is stopped.
**Next question:** solve the complete equations for changed joint walls and $`R=VPV^\dagger`$.
Two-band V and a one-band copy would give at most 42 indicator bands; no such family is known.
Every inverse and common budget is due. Retire exact ansatzes with no compatible seed;
approximations require their own error certificates. Small fixtures must check a derived rule.
If only necessary screens remain, reassess the [whole two-bank table word](research/depth/TABLE_SPECIFIC_INDICATOR_INTERFACES.md#5-what-this-changes-and-what-it-leaves-unpaid).

The [four-band Clifford front end](research/depth/CLIFFORD_INVOLUTION_SOURCE.md),
[weighted outer-phase route](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md),
[reduced central gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
and [minimal-tree probe](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source)
remain component alternatives under their separate screens. In the separate
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
