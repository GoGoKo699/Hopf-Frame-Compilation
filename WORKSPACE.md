# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge are authorized. Optimal general T-depth remains the destination;
final manuscript drafting is a separate task.

## Current checkpoint: two-band source gap and native profile scout

**Fetched starting main:**
`eb368d2fcccfac2c1dc22d66d0e8c8d979891384`, tree
`715b96aa427daa89ce1fffd15a024b4494959740`, on 19 October 2026.
Read the [current session outcome](research/depth/TWO_BAND_SOURCE_GAP_SESSION_RESULT.md)
for this continuation and its verification boundary. A fixed Clifford-sum
gap excludes the complete two-band indicator source at sufficiently small
fixed error and sufficiently large address size, for fixed amplitude and
uniformly bounded clean count. The original physical-pair scout gains
native profile constraints without assuming bounded silent span. No faster
indicator or unrestricted depth bound is claimed. Fetch current main
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
| Current source obstruction | [Fixed Clifford-sum gap](research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md): for fixed positive amplitude and uniformly bounded clean count, actual two-band sources fail at some fixed positive accuracy for all sufficiently large p, regardless of count or width. No flatness premise is needed. The general constants are conservative; [finite-p overlap screens](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md) remain stronger on their particular families |
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

**Resume outside the stopped two-band source class.** Read the
[fixed-sum proof](research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md), the
[complete transitions and dressed screens](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md),
and the [three-band continuation recurrence](research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#62-the-next-source-band-retains-the-general-leakage-recurrence).
The first unpaid operation is a **complete normalized logical transition**
at amplitude $`1/\sqrt2`$ from an actual one-signal source with at least
three physical T bands, or another explicitly paid interface. It must
have uniform linear source T-count/dirty width and quadratic Clifford
count at each requested accuracy. Recovery and clean removal are paid.
For three bands their conditional depth is at most 104, and promoted
error is at most $`16(3\sqrt2\epsilon+2\epsilon^2)`$; source tolerance
$`\epsilon\le\tau/72`$ suffices. Signal/toggle return is approximate
and included in this complete operator error unless exact flatness is
proved separately. For R98 use $`\epsilon_L=2^{-L}/(1152K_*)`$.

The new obstruction uses a fixed-length, fixed-coefficient sum of actual
Cliffords for a projected output-probe defect. Its normalized Hilbert–Schmidt norm
cannot tend to zero without vanishing exactly, contradicting the specified
indicator action. The theorem stops the every-fixed-accuracy interface;
its conservative constants do not exclude every modest numerical tolerance.
Allowing clean count to grow without a uniform bound, or changing block
amplitude, requires a new quantitative recurrence. Uniformly bounded
variation of clean count is already covered by padding. The theorem does
not exclude unrestricted constant depth or prove a dimension-growing bound.

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
The [central-module profile theorem](research/depth/CENTRAL_MODULE_NATIVE_PROFILES.md)
now gives three simultaneous native profile distortions bounded by
$`8c,4c,4c`$, with $`c\le3s+1`$, and a per-element escaped-axis span
at most $`8c`$ modulo the genuine center. These hold without bounded t.
The next missing operation is simultaneous native matching across the
actual square-zero displacement family, including axes that stay inside
the original middle span plus center. Per-element escape does not bound
the union of escaped axes. Fixed positive s alone, growing s, and higher
layers remain open; no bounded-t or degree conclusion is inferred.

### Stopped routes, reopening conditions and superseded instructions

| Route or historical next step | Current decision |
| --- | --- |
| General overlapping two-band source at fixed amplitude | [Stopped at the every-fixed-accuracy scope](research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md), with arbitrary dirty width and uniformly bounded clean count. Reopen with at least three source bands, a different interface, or a quantitatively paid change in unbounded clean count/amplitude |
| Two-band half-amplitude block | [Stopped by the new screen](research/depth/JOINT_FOURIER_WORD_LIMITS.md#51-two-source-bands-at-half-amplitude-reduce-to-two-clean-bands): near-half flatness has a literal Clifford aligner, so recovery still has only two T bands. Reopen with a different normalization, more bands or a different source interface |
| Parity-copy overlap repair | [Stopped for p at least four](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md): undoing the signal contraction gives a signed permutation whose output probe is still a two-band circuit. Reopen outside that parity-copy signal architecture |
| Growing signal-mixing rank-one-product family | A paid Clifford repair defeats the old wrong-square screen. The [quarter-turn coefficient screen](research/depth/OVERLAPPING_TWO_BAND_SOURCES.md#6-a-clifford-invariant-coefficient-screen-for-quarter-turn-transitions) excludes every independent Clifford dressing at explicit finite-p error. Reopen with a transition outside that quarter-turn class, also respecting the general two-band obstruction |
| Simultaneous nonscalar output probes | [Parity support constraints](research/depth/SIMULTANEOUS_OUTPUT_PARITY_LIMITS.md) force linear native support, but a two-band identity transition saturates both directions. Rank saturation alone cannot supply the indicator's specified action; the fixed-sum theorem resolves the general two-band route |
| Separated high-rank anticommuting reflections | [Stopped example](research/depth/JOINT_FOURIER_WORD_LIMITS.md#52-paused-discovery-and-the-next-complete-transition): disjoint source bands merge into one physical band. Reopen with actual overlapping dynamics; rank alone does not help |
| Independent row/pulse synthesis | [Native sharing](research/depth/NATIVE_INDICATOR_SHARING.md) forces superlinear emitted count for too-small completed disjoint bank blocks. Reopen with growing shared bank action or joint synthesis across those boundaries |
| Earlier filter, Gaussian and phase-signal routes | Read [filter limits](research/depth/APPROXIMATE_INDICATOR_FILTER_LIMITS.md), [approximate Gaussian limit](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md#4-constant-error-does-not-repair-the-gaussian-conjugator), and [signal/frame limits](research/depth/PHASE_SIGNAL_QUERY_LIMITS.md). Their failures retain full-bank error, separate-source cost, or unpaid coherent correction; reopening must explicitly leave or overcome the stated hypotheses |
| PR #195 growing-total-displacement question | Superseded by PR #196's [actual no-middle-pair closure](research/depth/GROWING_DISPLACEMENT_NATIVE_ORBITS.md). Original physical pairs are a different remaining case; do not reopen the closed no-pair branch |
| Zero-clean discovery only; one-band inverse-square-root-two; unpaid flat-sector alignment | Superseded by PR #198's clean promotion, one-band screen and paid alignment. The remaining task is the logical transition, not those operations |
| Old precision-envelope/opening/local-compiler next steps | Superseded by R94–R96. R96 is all-precision; earlier headline bounds remain historical with their useful component proofs |

The [latest session outcome](research/depth/TWO_BAND_SOURCE_GAP_SESSION_RESULT.md)
records the source gap and native profile continuation. Earlier prompts to
search general two-band transitions or merely nonscalar output probes are
superseded. The unrestricted pair-count question remains open.
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
| Research proofs and historical outcomes | [Research index](research/README.md), [current session record](research/depth/TWO_BAND_SOURCE_GAP_SESSION_RESULT.md) |
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
