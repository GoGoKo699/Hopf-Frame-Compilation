# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames. Manuscript drafting is a separate task.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint: complete signal, scalar-kernel and Bell source trials

**Recovered main:** [PR #220](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/220),
`6f5702291016eb664d75c1ae589abd69ec31f2bf`, fetched 21 October 2026.
This checkpoint constructs signal-interleaved, scalar-kernel and coherent
Bell source families, with complete operator and native ledgers. Their
lookup-sector failures are proved below; no qualifying source is supplied.
Exact-tree receipts belong to the PR record. Fetch main before resuming; this is not a pin.

Read in this order:

1. [Resource frontier](docs/OPEN_PROBLEM.md): general theorems, literal allocations and lower-bound scopes.
2. [Signal-interleaved sources](research/depth/SIGNAL_INTERLEAVED_KERNEL_SOURCE.md): actual controlled phases, accepted polynomials, all-size failures and the changed final coupling.
3. [Scalar-kernel partial source](research/depth/SCALAR_KERNEL_FILTER_SOURCE.md) and [corank-one Bell source](research/depth/CORANK_FOURIER_BELL_SOURCE.md): correct support versus normalization, with complete arbitrary-input actions.
4. [Four-call Gram recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) and [closed six-wire scout](research/depth/FOUR_BAND_CENTER_TRANSPORT.md#coupled-six-wire-endpoint-obstruction).
5. [Route decisions](research/ROUTE_DECISIONS.md) before reopening an older mechanism; the [research catalogue](research/README.md) locates every research file.

Historical prompts and session endings are records, not current instructions.

<a id="current-depth-frontier"></a>

### Destination and strongest established statements

For $`N=2^n`$, integers $`n\ge1,L\ge6`$, the strongest unconditional
general upper theorem remains [R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md):

```math
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL).
```

One clean qubit suffices at real width $`b\ge17(L+n+7)`$ or
complex-magnitude width $`b\ge17(L+n+8)`$, with complete
initialized-isometry error at most $`2^{-L}`$. These are simultaneous
bounds with absolute constants. [R79](research/depth/UNIFORM_SUBTREE_DEPTH.md)
remains selected Result D. At fixed accuracy and square-root dirty width,
the unrestricted depth interval is still logarithmic upper versus constant
lower. The latest source and representation results do not close it.

The exhaustive at-most-one-clean three-layer Boolean theorem has universal
sign degree at most $`D_{\rm np}=D_3+576=480\,2^{43}+903`$ and nonconstant
XOR degree at most $`D_{\rm np}-1`$, with rounding radius strictly below
$`2^{-530690}`$; see the [exact case split](research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md).
The paired six/five bounds are a separate branch, not universal bounds.

### Main unfinished route and exact missing operation

For p address bits and $`m=2^p`$, construct the complete dirty indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle
```

at every fixed requested tolerance tau, with zero clean work,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`D_T=o(\log(p+2))`$. All address, output, dirty and reference inputs
are arbitrary. Intermediate address/output motion is allowed; the
complete operator and every return error must be certified.

[Polynomial-overhead fusion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
permits fixed polynomial count/width overhead in a smaller complete source.
Its eight actual source/adjoint calls and logarithmic-remainder tail are
paid. [Restored-address R98](research/depth/QUERY_DEPTH_FRAME_TRANSFER.md)
and its [constant-accuracy extension](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
transfer qualifying queries to frames, retaining their log-log contribution.
Assess success on that complete ledger, not source depth alone.

**Constructive question.** Build a normalized one-clean source U with
accepted block $`K=J^\dagger UJ`$ and $`K^\dagger K\approx\Pi_{\rm eq}/2`$,
where $`\Pi_{\rm eq}=(I+\Omega_{x,C}\Omega_{x,D})/2`$, using higher-band
native interference, or find another complete wide non-Clifford readout.
The [four-call Gram interface](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
pays recovery and serial dirty promotion, including arbitrary accepted fibers,
address motion and actual inverses. Gram error $`\delta\le\tau^2/4096`$
suffices. The six-call linear block-error interface remains available.
**The qualifying lookup source is unpaid.** The
[signal-interleaved families](research/depth/SIGNAL_INTERLEAVED_KERNEL_SOURCE.md)
now pay the controlled phases and actual accepted polynomials. The
attenuated two-clean word fails for every inner unitary; the direct
one-clean word of depth at most 314 fails for its flat Fourier mixer.
Both have Gram error at least one quarter at every nonzero address.
The latter proof also identifies the exact phase-eigenspace condition
needed to reject its offending inputs; another flat mixer is insufficient.

The [scalar-kernel source](research/depth/SCALAR_KERNEL_FILTER_SOURCE.md)
has constant fixed-accuracy depth. Its ideal Gram has the correct
half-projector normalization and rejects the entire wrong charge, but
misses most of the correct sector and has gap one half. Scalar spectral
processing cannot distinguish its missing-even and odd zero eigenvectors.
The [corank-one Fourier matrix](research/depth/CORANK_FOURIER_BELL_SOURCE.md)
has the desired abstract Gram, but shared Bell outcomes impose extra
stabilizers and dilute the accepted weight. Its actual one-clean
coherent-SWAP repair has depth at most 42 and exact gap one half at every
size. All these source budgets fit polynomial fusion; their sector
contracts fail. These are analytic, scoped failures, not general depth bounds.

**Next concrete question.** Change the final phase coupling in the
one-clean signal word by the partial-Hadamard bank rotation specified in
its [proof home](research/depth/SIGNAL_INTERLEAVED_KERNEL_SOURCE.md).
This preserves the correct zero-address block and adds only Clifford
walls. First determine the intersection of its rotated fixed-phase
space with the wrong lookup sector; then derive the complete Gram if
that screen is passed. Its full Gram identity is unresolved. A new
coupling must change this sector condition, not merely extraction
constants, flat mixing or precision. Pause the displayed failed words.

| Read before proposing a source | Exact scope |
| --- | --- |
| [Native read-only two-band bound](research/depth/READ_ONLY_TWO_BAND_SOURCE_GAP.md) | At amplitude $`1/\sqrt2`$ and relative error at most $`1/8`$, $`k\ge\lceil(p-3)/2\rceil`$ clean inputs are necessary when every actual Clifford boundary fixes the address. No flatness or exact return is assumed. |
| [One-clean two-band Gram screen](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-band-projector-gram-gap) | For p at least three, Gram error at most $`1/256`$ is impossible even with address motion. This target needs at least three source bands or a different clean allocation. |
| [Actual Clifford signal-image screen](research/depth/CROSSED_CUBIC_KERNEL_SOURCE.md#clifford-signal-gram-gap) | If $`U^\dagger Z_aU`$ is Clifford, lookup Gram error is at least $`1/(2\sqrt2)`$ for $`p\ge1`$. Arbitrary Clifford address motion is allowed; a general higher-band source need not satisfy this premise. |
| [Sharp Clifford-filter limit](research/depth/CLIFFORD_LOOKUP_FIXED_SPACE.md) | Arbitrary wide Clifford lookup correlation is at most $`(p+1)/2^p`$. Low-budget Clifford averaging is not a normalized source; this is not a T-depth bound. |
| [General two-band rank–norm screen](research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen) | Exact sources have a growing clean lower bound. The approximate floor shrinks with clean count; fixed-error moving-address growing-clean compression remains open. |

The [spectral-promotion window](research/depth/MSS_SPECTRAL_PROMOTION.md)
$`k=o(\log(p+2))`$ still requires a source, uniform precision and native
address sharing. Completed address return does not justify parallel native
calls. The paid [scalar zero reflection](research/depth/EXPANDER_ZERO_REFLECTION.md)
has zero clean work and constant fixed-accuracy depth. Conjunction alone
cannot witness dimension-growing depth at fixed accuracy.

The [paid fermionic lift](research/depth/BANK_MOTION_WALSH_INDICATOR.md#expander-fermionic-lift)
remains paused: complete spectral accuracy forces filter length
$`N=\Omega_\tau(p)`$ on unpadded lengths at constant auxiliary dimension,
regardless of selector precision; displayed coherent centering also
restores growing depth. This is a restricted sparse-word obstruction,
not a general native-depth bound; the scalar zero test is unchanged.

**Bounded lower scout completed.** The [genuinely coupled six-wire middle](research/depth/FOUR_BAND_CENTER_TRANSPORT.md#coupled-six-wire-endpoint-obstruction)
is paid (ten T gates, two bands, no helpers), but its transported probe sends
$`X_1`$ to four Paulis of magnitude one half. Every one-band output endpoint
has complete clean-return error at least one quarter, even with spectator
Clifford mixing and arbitrary logical action. Pause this candidate. Reopen
only with changed actual phase algebra or a paid non-Clifford endpoint,
then prove full Boolean return. The coefficient test can inspect more
layers; its witness is not invariant under extra non-Clifford layers.

The PR #220 signal-interleaving question and PR #219 mixer proposal are
superseded by the actual words and changed final coupling above. The
held PR #217 queue is closed: four-call recovery and the exact Gram-sector
conditions are proved; the six-wire candidate is stopped. Fourth-order
closure supplies neither lookup sectors nor a free native square.

### Stopped routes, reopening conditions and superseded instructions

[Route decisions](research/ROUTE_DECISIONS.md) records stopped mechanisms,
their reopening conditions and settled questions. Unique attempts and proofs
remain in their homes; a representation limit does not retire the destination.

<a id="research-decision-and-stopping-rules"></a>
<a id="continuing-toward-the-general-answer"></a>
<a id="read-by-claim"></a>

## Working contract and reading map

Prefer analytic identities and a few exact fixtures. Preserve every logical
column, literal phase, arbitrary dirty/reference input, actual inverse and
clean leakage/return error. A suffix is clean only on its active sector;
inactive action needs its own proof. No QRAM, free data loading, supplied
precision states, unpriced clean history, resets or postselection.
The selected manuscript and changed-decoder supplement remain separate.

Before a fixture, state the native identity and complete resource/error/return
recurrence. Change mechanisms when completion restores the old asymptotic
cost. Distinguish general theorems, conditional interfaces and
representation-specific limits. Internal audits and finite checks do not
replace an analytic proof or external technical review.

| Need | Authoritative home |
| --- | --- |
| Selected claims and exact allocations | [Publication scope](manuscript/PUBLICATION_SCOPE.md), [core audit](docs/CORE_CLAIM_AUDIT.md) |
| Argument and proof dependencies | [Review](REVIEW.md), [documentation map](docs/README.md) |
| Current general resource frontier | [Resource questions](docs/OPEN_PROBLEM.md) |
| All research proofs and outcomes | [Research catalogue](research/README.md) |
| Attribution and finite evidence | [Source catalogue](docs/reference/SOURCE_CATALOGUE.md), [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md) |

<a id="refurnishing-policy"></a>

Each result has one proof home. Navigation pages summarize their own role
and link; do not append a session narrative to every front page.
Preserve unique failed attempts, stable identifiers, inbound anchors,
scientific tests and receipt hashes. Navigation regressions check catalogue
coverage and anchor compatibility.

## Restore and verify

Check worktree status and remote, fetch main, then inspect changes since
the baseline. Preserve unrelated dirty worktrees; never reset them.
[SYNC](SYNC.md) is authoritative for tracked upstream versions; this
scientific checkpoint has no pending provenance change.

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

Use the [presentation guide](assets/README.md#rendering-checks) for the
rendered-document gate. Record exact-tree receipts in the review record;
merge only after required checks pass. A fixture is not an asymptotic proof
or a scalable native emitter.
