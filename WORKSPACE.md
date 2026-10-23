# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames. Manuscript drafting is a separate task.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint: coupled-source coherence and matching invariants

**Verified scientific checkpoint:** [PR #228](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/228),
`ea8f2ed40a8c85b465296703cd6675072b77eb03`, checked 23 October 2026.
Coupled two-path and four-path sources are paid but fail the complete Gram.
Matching invariants and a one-row spectral escape are proved; no qualifying
complete source is supplied. Exact-tree receipts belong to the PR record.
Fetch main before resuming; this is not a pin.

Read in this order:

1. [Resource frontier](docs/OPEN_PROBLEM.md): general theorems, literal allocations and lower-bound scopes.
2. [Coupled quarter-phase sources](research/depth/COUPLED_QUARTER_PHASE_SOURCE.md): complete accepted polynomials, the marked/unmarked eigenvalue collision, and the four-path coherent Gram defect.
3. [Matching limits and spectral escape](research/depth/QUARTER_PHASE_MATCHING_LIMITS.md), then [balanced two-branch readout](research/depth/COUPLED_ORDER_GRAM_SOURCE.md#balanced-two-unitary-gram-readout): exact hypotheses and the direct zero-clean square reduction.
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

**Constructive question.** Build a normalized bounded-clean source U with
accepted block $`K=J^\dagger UJ`$ and $`K^\dagger K\approx\Pi_{\rm eq}/2`$,
where $`\Pi_{\rm eq}=(I+\Omega_{x,C}\Omega_{x,D})/2`$, using higher-band
native interference, or find another complete wide non-Clifford readout.
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

**Latest substantive outcome.** The [matched-flip two-stage word](research/depth/COUPLED_QUARTER_PHASE_SOURCE.md)
has a complete arbitrary-bank relative operator, but a marked/unmarked
common eigenvalue defeats scalar filtering. Reversing branch order also
fails. A paid four-path word with bank-sensitive phases passes the empty-bank
and one-row basis-probability checks, yet has complete Gram error
exactly one half. A fixed-row selector flip does give an exact
selector-independent partial Gram on arbitrary banks.
[Matching stages](research/depth/QUARTER_PHASE_MATCHING_LIMITS.md) cover all
coordinates at constant depth, but fixed partitions and fixed bipartitions
have proved invariant-vector limits. The full dressed expander family has
a uniform one-row spectral gap; that is not an arbitrary-bank Gram theorem.

**Next concrete question: complete bank-bit response.** For a bank-preserving
candidate write $`G_Y=K(Y)^\dagger K(Y)`$. In the single-bank coordinates
of the coupled-source proof, the target is $`(I+\Omega_Y)/4`$, with
$`\Omega_Y|x\rangle=(-1)^{Y_x}|x\rangle`$. Seek an actual accepted word with

```math
G_0=I/2,\qquad
G_{Y\oplus e_v}-G_Y
=-\frac{(-1)^{Y_v}}2|v\rangle\langle v|\otimes I_{\rm work}
```

for every row v and every background bank Y. These exact identities telescope
to the complete target; they are a test for a proposed source, not a new source.
This diagnostic is scoped to bank-preserving candidates; the general source
model may move the bank.
Start from the paid row-selective partial Gram and derive its background
conjugations when distinct dressed matching stages are combined. The first
calculation is the full two-active-row operator, retaining every selector and
off-diagonal address term. Do not infer arbitrary-bank sector invariance from
the one-row gap: another active row can move a marked address.

Proceed only if an exact cancellation identity or a global norm estimate
survives at every bank. A per-flip error epsilon alone accumulates to
$`|Y|\varepsilon\le m\varepsilon`$; accuracy of order delta/m can restore the
old precision/depth cost. Pay any coherent row selection, signal mixing and
actual inverse before claiming a complete improvement. Stop a proposed repair
if it assumes a free address-controlled row flip, $`\Theta(p)`$ clean history
or a correction that reconstructs the known logarithmic-depth indicator.
A scalar filter of the failed word or another basis-only success is insufficient.
If the word reduces to two balanced unitary branches, use the
[direct-square criterion](research/depth/COUPLED_ORDER_GRAM_SOURCE.md#balanced-two-unitary-gram-readout).
Recovery is already paid. The earlier instruction to merely couple stages
or transfer the one-row spectral gap is superseded by this response test.

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
