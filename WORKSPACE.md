# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge are authorized. Optimal general T-depth remains the destination;
final manuscript drafting is a separate task.

## Current checkpoint: overlapping two-band sources

**Fetched starting main:**
`0f973550fba2996b484e3dc8c67b688a5be5862f`, tree
`c5e347bd8abc9eabd891bcbc9a10a08d5870c643`, on 19 October 2026.
Read the [current session outcome](research/depth/OVERLAPPING_SOURCE_SESSION_RESULT.md)
for this continuation and its verification boundary. The new results are
source-specific overlap screens, an exact-flatness error improvement,
and a restricted original-pair degree theorem. No faster indicator or
complete-frame frontier improvement is claimed. Fetch current main
before resuming; the checkpoint's starting commit is not a pin to an
older branch.

<a id="current-depth-frontier"></a>

### Destination and strongest established statements

Write $`N=2^n`$, with integer $`n\ge1,L\ge6`$.

| Role | Statement and exact scope |
| --- | --- |
| Strongest unconditional general upper | [R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md): $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, with absolute constants on the same circuit |
| Allocation and error | One clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; complete initialized-isometry error at most $`2^{-L}`$, including arbitrary dirty inputs and references |
| Selected manuscript theorem | [R79](research/depth/UNIFORM_SUBTREE_DEPTH.md) remains Result D. The selected count frontier and literal allocations are settled; no practical crossover, efficient preprocessing or general native emitter is claimed |
| Paid constructive transfer | [Restored-address R98](research/depth/QUERY_DEPTH_FRAME_TRANSFER.md) and its [approximate extension](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md) need only constant source accuracy at each fixed final L, with uniform simultaneous count/width/Clifford/depth budgets |
| Latest positive component | [Clean-source promotion](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) removes a uniformly bounded clean bank, permitting arbitrary off-code action; error is at most $`8\cdot2^k\epsilon`$. One-signal complete blocks at amplitude $`1/2`$ or $`1/\sqrt2`$ suffice through 24 source/actual-inverse calls. [Exactly flat two-band recovery](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#61-exact-flatness-sharpens-the-complete-error-and-return-recurrence) improves promoted error to $`2\sqrt2\epsilon`$ and returns the borrowed signal/toggle exactly |
| Unrestricted indicator count | [Pauli nullity](research/depth/INDICATOR_PAULI_NULLITY.md): for $`p\ge2`$, a zero-clean full-operator implementation with error below $`1/\sqrt2`$ needs at least $`2^p-1`$ T gates, with arbitrary dirty width and intermediate mixing. This is not a depth bound |
| Current source screens | [Joint-word limits](research/depth/JOINT_FOURIER_WORD_LIMITS.md), with one signal: one physical source band fails at both paid amplitudes for $`p\ge3`$. For $`p\ge3`$, the new two-band half-amplitude floor is $`\sin(\pi/8)/(4\sqrt2)`$, with one signal and arbitrary dirty width/mixing. Two-band inverse-square-root-two remains open; [new overlap screens](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md) exclude parity-copy encoders for $`p\ge4`$ and any sufficiently accurate source with even one scalar-commuting output probe |
| Lower-bound scope | [Two-layer Boolean rigidity](research/depth/TWO_LAYER_DIAGONAL_RIGIDITY.md) gives XOR degree at most two with at most one clean qubit and robust floor $`\sin(\pi/8)/(2\sqrt2)`$ for higher degree. [Zero-clean three-layer rigidity](research/depth/THREE_LAYER_NATIVE_STRUCTURE.md) gives XOR degree at most four. The unrestricted one-clean three-layer degree and dimension-growing depth questions remain open |

The [resource frontier](docs/OPEN_PROBLEM.md) is authoritative for the
full theorem map, including the zero-clean five-layer lower bound and
R97's scoped outer-axis obstruction. At fixed accuracy and square-root
dirty width, the general upper/lower interval remains logarithmic versus
constant. An interface failure does not close that unrestricted gap.

### Main unfinished route and exact missing operation

For p address bits and $`m=2^p`$, construct the complete indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle
```

at every fixed requested tolerance tau, with zero clean work,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$, and
$`D_T=o(\log(p+2))`$. Require full-operator accuracy on arbitrary
address, output, dirty-work and reference inputs, and charge the actual
inverse and all return errors. Address/output motion within a completed
query is allowed. Success would improve the fixed-accuracy frame law
through the paid R98 transfer; the resulting frame still has its retained
$`\log_2\log_2(n+4)`$ contribution.

**Resume with genuinely signal-mixing two-band sources at inverse-square-root-two
amplitude.** Read [the complete transitions and probe identity](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md),
then the [exact-flatness recurrence](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#61-exact-flatness-sharpens-the-complete-error-and-return-recurrence).
Actual one-band encoders U_A,U_B with source error at most 0.1 have
exactly anticommuting clean-image Clifford reflections A,B. The first
unpaid operation is still their **complete normalized logical transition**
approximating I_p with uniform linear source T-count/width and quadratic
Clifford count. Alignment, recovery and clean removal are paid; their
conditional depth is at most 80 and error at most $`2\sqrt2\epsilon`$.

The next concrete question is whether actual encoders can satisfy all
indicator output-probe identities simultaneously when every original
output probe has a nonscalar commutator with A. The exact residual
$`T_j=C_{j,+}+C_{j,-}B`$ is unitary, but is generally not Clifford;
replacing its four-term Clifford sum by a free gate is the first
unsupported inference. A complete transition, its actual inverse and
all return errors are required. More than two source bands and other
bounded-clean mechanisms remain allowed.

Keep the [original physical-pair scout](research/depth/OCCUPIED_BOOLEAN_PAIR_CENTER.md#8-a-proved-bounded-span-reopening-with-shared-feature-coordinates)
active. With s original pairs and additional actual coupled silent-axis
span t, the new theorem bounds Boolean degree by $`D_K+4A_K`$,
$`K=3s+2+t`$, uniformly in width and endpoint active rank. Its shared
computational coordinates keep the physical clean bit until after the
feature calculation. **Missing general implication:** bound the actual
dressed X-mask span and an adequate feature ledger from s alone, without
assuming bounded t. Bounded coupling rank is still not bounded axis
span. The [clean-preserving Clifford-sum reduction](research/depth/OCCUPIED_BOOLEAN_PAIR_CENTER.md#9-a-clean-preserving-clifford-sum-without-bounding-the-silent-span)
keeps a bounded scalar-feature ledger even without bounded t: every
summand is a full logical Clifford generated by the same actual silent
band. Its displacement maps are linear in the bounded residual label
and have pairwise zero products, but this does not bound individual ranks.
Use that shared native structure in the missing Boolean matching step.
Fixed positive s alone, growing s, and higher layers remain open.

### Stopped routes, reopening conditions and superseded instructions

| Route or historical next step | Current decision |
| --- | --- |
| Two-band half-amplitude block | [Stopped by the new screen](research/depth/JOINT_FOURIER_WORD_LIMITS.md#51-two-source-bands-at-half-amplitude-reduce-to-two-clean-bands): near-half flatness has a literal Clifford aligner, so recovery still has only two T bands. Reopen with a different normalization, more bands or a different source interface |
| Parity-copy overlap repair | [Stopped for p at least four](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md): undoing the signal contraction gives a signed permutation whose output probe is still a two-band circuit. Reopen outside that parity-copy signal architecture |
| Growing signal-mixing rank-one-product family | [Stopped by its complete transition](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md#4-a-growing-signal-mixing-family-with-its-complete-transition): the logical unitary squares to a nontrivial quadratic Clifford, giving phase-optimized source error at least one half against every involution. Reopen with a changed transition; arbitrary independent Clifford dressing is not excluded by this square argument |
| Separated high-rank anticommuting reflections | [Stopped example](research/depth/JOINT_FOURIER_WORD_LIMITS.md#52-paused-discovery-and-the-next-complete-transition): disjoint source bands merge into one physical band. Reopen with actual overlapping dynamics; rank alone does not help |
| Independent row/pulse synthesis | [Native sharing](research/depth/NATIVE_INDICATOR_SHARING.md) forces superlinear emitted count for too-small completed disjoint bank blocks. Reopen with growing shared bank action or joint synthesis across those boundaries |
| Earlier filter, Gaussian and phase-signal routes | Read [filter limits](research/depth/APPROXIMATE_INDICATOR_FILTER_LIMITS.md), [approximate Gaussian limit](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md#4-constant-error-does-not-repair-the-gaussian-conjugator), and [signal/frame limits](research/depth/PHASE_SIGNAL_QUERY_LIMITS.md). Their failures retain full-bank error, separate-source cost, or unpaid coherent correction; reopening must explicitly leave or overcome the stated hypotheses |
| PR #195 growing-total-displacement question | Superseded by PR #196's [actual no-middle-pair closure](research/depth/GROWING_DISPLACEMENT_NATIVE_ORBITS.md). Original physical pairs are a different remaining case; do not reopen the closed no-pair branch |
| Zero-clean discovery only; one-band inverse-square-root-two; unpaid flat-sector alignment | Superseded by PR #198's clean promotion, one-band screen and paid alignment. The remaining task is the logical transition, not those operations |
| Old precision-envelope/opening/local-compiler next steps | Superseded by R94–R96. R96 is all-precision; earlier headline bounds remain historical with their useful component proofs |

The [latest session outcome](research/depth/OVERLAPPING_SOURCE_SESSION_RESULT.md)
records the overlap continuation. The previous proposed bounded-silent-span
reopening is now proved; the unrestricted pair-count question remains open.
Older session “next steps” are historical unless
reaffirmed here. PR #196 already reconciled the documentation-only
Hopf-QBP difference; [SYNC](SYNC.md) retains tracked upstream
`99e5ee5fc77004f662462e5d088e2a583d032f6a`. No further provenance task
is pending.

<a id="research-decision-and-stopping-rules"></a>
<a id="continuing-toward-the-general-answer"></a>
<a id="read-by-claim"></a>

## Working contract and reading map

Prefer analytic identities and few exact fixtures. Do not add QRAM,
free data loaders, supplied precision states, unpriced clean history,
large simulations, or an application-advantage project. Literal phases,
all logical columns, arbitrary work/reference inputs, actual inverses
and leakage/return must be charged. A suffix is clean only on its active
sector; inactive action needs a separate proof. The selected manuscript
and completed state-based QBP supplement remain separate.

| Need | Authoritative home |
| --- | --- |
| Selected claims and literal allocations | [Publication scope](manuscript/PUBLICATION_SCOPE.md), [core audit](docs/CORE_CLAIM_AUDIT.md) |
| Claim-led explanation and dependency map | [Review](REVIEW.md), [proof map](docs/README.md) |
| Full frontier and complete-circuit requirement | [Resource questions](docs/OPEN_PROBLEM.md) |
| Research proofs and historical outcomes | [Research index](research/README.md), [current session record](research/depth/OVERLAPPING_SOURCE_SESSION_RESULT.md) |
| Attribution and evidence | [Source catalogue](docs/reference/SOURCE_CATALOGUE.md), [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md) |

Separate discovery from certification. Before a new fixture, state the
native identity and resource/error/cleanup recurrence. Stop a candidate
that retains its old asymptotic cost, needs unpaid clean history, or
cannot return work; record the exact reopening condition. A session end
or a representation-specific obstruction does not retire optimal general
depth. No manuscript drafting until requested.

<a id="refurnishing-policy"></a>

Keep one proof home per claim. Front pages summarize and link. Preserve
unique failed attempts, scoped positive results, stable identifiers and
inbound anchors; do not duplicate proofs or weaken scientific tests and
receipt hashes. Internal audits do not replace external technical review.

## Restore and verify

Start with `git status --short --branch`, verify the repository remote,
fetch main, and inspect any change from the baseline above. Read this
checkpoint and the linked current source homes before older prompts.
No private scratch file or conversation transcript is needed to resume.
For substantive theorem/circuit changes, run:

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
complete rendered-document gate. Record checks on the actual reviewed
tree and merge only after required gates pass. Finite checks support
literal identities; they are not asymptotic proofs or scalable emitters.
