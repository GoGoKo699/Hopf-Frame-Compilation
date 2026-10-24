# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames. Manuscript drafting is a separate task.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint: joint native transport and the complete-bank gap

**Verified starting checkpoint:** [PR #238](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/238),
`3b04fc2f0655fbbbfff7c3a3ef56f57107b8dbb1`, fetched 24 October 2026.
A three-layer source now pays joint transport of two address coordinates.
A four-layer source pays a nested pair of bank coefficients, including a quartic term.
A bank-independent final block and one specified full-bank X-motion encoder
cannot complete the three-layer general source.
The complete-bank construction remains open. Exact-tree receipts are in the PR record.
Fetch main before resuming; this is not a pin.

Read in this order:

1. [Resource frontier](docs/OPEN_PROBLEM.md): general theorems, literal allocations and lower-bound scopes.
2. [MOD3 zeta transport](research/depth/MOD3_ZETA_TRANSPORT.md): complete observable contract and three-layer exactness; [joint transport](research/depth/WEAK_PAULI_ORBIT_COMPRESSION.md#joint-two-address-source): the paid two-coordinate source; [coupled final encoders](research/depth/COUPLED_FINAL_ENCODER.md): complete operators and two scoped three-layer boundaries.
3. [Bank response](research/depth/BANK_BIT_RESPONSE_AND_ROW_ISOLATION.md), [co-moving stages](research/depth/CO_MOVING_LOOKUP_STAGE.md), [selector readout](research/depth/CO_MOVING_SELECTOR_READOUT.md), [entangling acceptance](research/depth/ENTANGLING_SELECTOR_ACCEPTANCE.md), and [bank-coupled dynamics](research/depth/BANK_COUPLED_SELECTOR_DYNAMICS.md): retained components and source screens. The [coupled-source](research/depth/COUPLED_QUARTER_PHASE_SOURCE.md) and [matching](research/depth/QUARTER_PHASE_MATCHING_LIMITS.md) proofs retain their earlier scopes.
4. [One-clean recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) and [joint two-clean recovery](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout): the paid return and error interfaces.
5. [Route decisions](research/ROUTE_DECISIONS.md) before reopening an older mechanism; the [research catalogue](research/README.md) locates every research file.

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
logarithmic upper versus constant lower depth remains open.

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

**Retained accepted-source interface.** A normalized bounded-clean source U with
accepted block $`K=J^\dagger UJ`$ and $`K^\dagger K\approx\Pi_{\rm eq}/2`$,
where $`\Pi_{\rm eq}=(I+\Omega_{x,C}\Omega_{x,D})/2`$, remains an alternative
to direct indicator transport. Its complete source is still unpaid.
The [four-call Gram interface](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
pays recovery and serial dirty promotion, including arbitrary accepted fibers,
address motion and actual inverses. Gram error $`\delta\le\tau^2/4096`$
suffices for one clean signal. Joint two-clean recovery and full dirty promotion
give error at most $`128\sqrt\delta`$; request $`\delta\le\tau^2/2^{14}`$,
or $`\tau^2/2^{20}`$ before eight-call fusion. No extra clean recovery flag is used.
**The qualifying lookup source is unpaid.** Use the
[quantitative source screens](docs/OPEN_PROBLEM.md#paid-source-interfaces-and-their-limits)
and [route decisions](research/ROUTE_DECISIONS.md#source-specific-representation-limits)
before proposing a new word; their hypotheses are not general depth bounds.

**Retained components.** [Stage and row echoes](research/depth/BANK_BIT_RESPONSE_AND_ROW_ISOLATION.md)
pay background cancellation and fixed-row compression; coherent row selection
still requires the missing indicator. The [co-moving construction](research/depth/CO_MOVING_LOOKUP_STAGE.md)
pays complete lookup-charge preservation with zero clean work and depth at most
116. The [fixed-row central word](research/depth/ENTANGLING_SELECTOR_ACCEPTANCE.md#fixed-row-central-word)
reads only $`Z_{Y_0}`$; its paid address transport restores linear depth.
The [bank-coupled revision](research/depth/BANK_COUPLED_SELECTOR_DYNAMICS.md)
pays swap/Hadamard controls and arbitrary principal t movement. Its specified
transpose-zeta and bank-Hadamard sources retain their quarter and $`1/64`$
Gram floors. These components and scoped acceptance failures do not exclude
general nonlocal sources or provide free addressed row flips.

**Latest substantive outcome.** [Joint two-coordinate transport](research/depth/WEAK_PAULI_ORBIT_COMPRESSION.md#two-coordinate-lift)
uses twelve T gates, three layers and no extra work at every bank size.
Fixed Clifford walls cost $`O(m^2)`$. It is the complete source at p=2;
for larger p it pays two coordinates only. Its actual echo includes all
bank, work and reference inputs. The old one-coordinate product still has
error two, and its prefix-rank requirement remains representation-specific.

The [four-band nested-coefficient source](research/depth/WEAK_PAULI_ORBIT_COMPRESSION.md#four-band-nested-coefficient-source)
pays $`(-1)^{d+bt(a+cv)}`$ with twelve T gates, 34 Cliffords and no
helpers. Here b,c,t are addresses and d,a,v are arbitrary bank coefficients.
Its fixed-bank lift costs $`34+pm`$ Cliffords; it still omits independent
lookup coefficients, so it is a partial construction with full lookup error two.

The [final-block theorem](research/depth/COUPLED_FINAL_ENCODER.md#bank-independent-final-boundary)
excludes a three-layer source at p at least three when its completed final
one-band block is independent of the nonconstant coefficient bank, even
with arbitrary middle Clifford boundaries and dirty width. A fixed output
Clifford does not repair it. A separate [full-bank X-motion encoder](research/depth/COUPLED_FINAL_ENCODER.md#full-bank-x-motion-boundary)
also fails with any one-band middle, despite moving every bank coordinate.
Neither result is an unrestricted three-layer exclusion.

**Next concrete question: complete the joint coefficient action.**
Can the four-band mechanism include the five missing independent coefficients
at three addresses in one native identity, then combine growing coefficient
blocks without serial address processing? Pay the
[complete observable contract](research/depth/MOD3_ZETA_TRANSPORT.md#pauli-transport-interface)
and its actual source/adjoint ledger. The two excluded final encoders are paused;
an unrestricted exact three-layer source remains another possible mechanism. The
[three-layer rounding theorem](research/depth/MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness)
makes any error below $`\tau_*=\sin(\pi/8)/(2\sqrt2)`$ exact; two layers
are excluded for p at least three. Four or more layers may instead use the
approximate contract at every fixed tolerance with depth $`o(\log(p+2))`$.
Every address, bank, work and reference input remains arbitrary. Source work
may move; its echo and actual inverse must return it. Controlled P_0 suffices
for controlled indicator. The two-active-coordinate question is settled;
the two excluded final-encoder mechanisms are paused. Constant source depth still
retains the transfer's log-log term; fixed tolerance gives no all-precision claim.

The [syndrome-star alternative](research/depth/COUPLED_FINAL_ENCODER.md#syndrome-star-coupling)
pays a scalar reflection and a formal bank-transition identity, but normalized
complete readout remains unpaid: its leaf-to-center map is many-to-one.
Neither exponentiating that formal sum nor selecting its bright direction is free.

### Stopped routes, reopening conditions and superseded instructions

[Route decisions](research/ROUTE_DECISIONS.md) records stopped mechanisms,
their reopening conditions and settled questions. Unique attempts and proofs
remain in their homes; a representation limit does not retire the destination.
The [closed six-wire scout](research/depth/FOUR_BAND_CENTER_TRANSPORT.md#coupled-six-wire-endpoint-obstruction)
needs changed phase algebra or a paid non-Clifford endpoint before reopening.

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
| Selected claims and exact allocations | [Publication scope](manuscript/PUBLICATION_SCOPE.md) |
| Argument and proof dependencies | [Review](REVIEW.md), [documentation map](docs/README.md) |
| Current general resource frontier | [Resource questions](docs/OPEN_PROBLEM.md) |
| All research proofs and outcomes | [Research catalogue](research/README.md) |
| Attribution and finite evidence | [Core audit](docs/CORE_CLAIM_AUDIT.md), [source catalogue](docs/reference/SOURCE_CATALOGUE.md), [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md) |

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
