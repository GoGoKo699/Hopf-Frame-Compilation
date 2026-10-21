# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for
complete Hopf frames. Manuscript drafting is a separate task.

## Current checkpoint: native source scope and complete readout

**Verified scientific baseline:** [PR #217](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/217),
merged 21 October 2026 as `6ba339245eda3f7473774caf80e843610c0c5dd3`,
tree `d5f0291cc6a4e6bbe4f54d6ab5d17f1bfc5c9c73`.
All ten local commands below, 1,464 deterministic tests, both Python CI
matrices and the rendered-document gate passed on that tree.
Fetch current main before resuming; this baseline is not a pin.

Read in this order:

1. [Resource frontier](docs/OPEN_PROBLEM.md): general theorems, literal allocations and lower-bound scopes.
2. [Latest substantive outcome](research/depth/READOUT_AND_NATIVE_SOURCE_REVISION_RESULT.md): revision assessment and complete paid readout ledger.
3. The proof homes for the two active questions below.
4. [Route decisions](research/ROUTE_DECISIONS.md) before reopening an older mechanism; the [research catalogue](research/README.md) locates every research file.

No conversation replay is needed. Historical prompts and session endings
are records, not current instructions.

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
The [paired projector-Gram interface](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout)
pays six-call recovery and serial dirty promotion, including arbitrary
accepted private fibers and intermediate address motion.
**The normalized shallow source is unpaid.** Paired-Pauli SELECT loses
exponential normalization; a forward XOR query or central-sign preservation
does not repair it.

| Read before proposing a source | Exact scope |
| --- | --- |
| [Native read-only two-band bound](research/depth/READ_ONLY_TWO_BAND_SOURCE_GAP.md) | At amplitude $`1/\sqrt2`$ and relative error at most $`1/8`$, $`k\ge\lceil(p-3)/2\rceil`$ clean inputs are necessary when every actual Clifford boundary fixes the address. No flatness or exact return is assumed. |
| [One-clean two-band Gram screen](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-band-projector-gram-gap) | For p at least three, Gram error at most $`1/256`$ is impossible even with address motion. This target needs at least three source bands or a different clean allocation. |
| [Sharp Clifford-filter limit](research/depth/CLIFFORD_LOOKUP_FIXED_SPACE.md) | Arbitrary wide Clifford lookup correlation is at most $`(p+1)/2^p`$. Low-budget Clifford averaging is not a normalized source; this is not a T-depth bound. |
| [General two-band rank–norm screen](research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen) | Exact sources have a growing clean lower bound. The approximate floor shrinks with clean count; fixed-error moving-address growing-clean compression remains open. |

The [spectral-promotion window](research/depth/MSS_SPECTRAL_PROMOTION.md)
$`k=o(\log(p+2))`$ still requires a source, uniform precision and native
address sharing. Completed address return does not justify parallel native
calls. The paid [scalar zero reflection](research/depth/EXPANDER_ZERO_REFLECTION.md)
does not automatically compose into a uniformly accurate whole bank.

**Focused lower scout.** Continue [center transport](research/depth/FOUR_BAND_CENTER_TRANSPORT.md)
and the [original physical-pair scout](research/depth/FOUR_BAND_ORIGINAL_PAIR_SCOUT.md).
The affected-silent reduction leaves at most $`\lfloor4k/3\rfloor`$
actual axes per band at pairing nullity k; k at most two is Clifford.
It never raises depth or adds workspace. Next test genuinely coupled
center transport outside that factorization: pay an actual output probe
and full Boolean return. Sector-compatible symplectic maps, a common
center and a smaller T-count do not supply a cheaper complete circuit.

The [held continuation notes](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/217#issuecomment-6307823212)
specify a six-wire candidate, exact Gram-sector matching conditions and
a four-call Gram-only observation. They are **not certified additions
to the proof homes**. Audit and test before promotion; fourth-order
spectral closure does not make its square a free Clifford.
The three-band paired sign/XOR bounds remain six/five (five/four at
active rank at least four), with a separate no-pair ceiling.
None is a general four-band or dimension-growing depth theorem.

### Stopped routes, reopening conditions and superseded instructions

[Route decisions](research/ROUTE_DECISIONS.md) is the compact map of why
mechanisms stopped, what must change to reopen them, and which earlier
questions are settled. Unique attempts and proofs remain in their original
homes. A representation-specific obstruction does not retire the destination.

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
