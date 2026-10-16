# Preparation map for the remaining depth problem

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Research index](../README.md) · [Continuing workspace](../../WORKSPACE.md)

**Historical baseline, superseded at fixed accuracy by
[R78](SUBTREE_PROJECTOR_DECODER.md) and uniformly in precision by
[R79](UNIFORM_SUBTREE_DEPTH.md).** The
[executed session](RESEARCH_SESSION_RESULT.md) replaced overlap recovery
with native subtree-projector decoding and proved polylogarithmic depth.

At fixed accuracy and square-root width, the later
[joint-depth session](JOINT_DEPTH_SESSION_RESULT.md) proves
logarithmic depth for a structured magnitude family, a complete
sub-square-logarithmic circuit under bounded relative density, and an
improved general terminal separately. It also records the remaining early-decoder
cost and exact cache/fingerprint counterexamples. The unrestricted
square-logarithmic upper bound at that checkpoint was R79's.

The subsequent [branch-orientation session](BRANCH_ORIENTATION_SESSION_RESULT.md)
proves fixed-accuracy logarithmic depth for affine branch choices (R83)
and exact log-times-iterated-log depth for arbitrary signed corners (R84),
with one clean qubit and the same count budgets. R83 also supplies the
priced positive affine reference in R82's extension. The later
[in-place conversion checkpoint](PERMUTATION_BRIDGE_SESSION_RESULT.md)
stopped without a qualifying construction and preserved that general
frontier. Its scoped obstructions inform the continuing search; the
[current mandate](../../WORKSPACE.md#continuing-toward-the-general-answer)
keeps optimal general depth as the destination. The later
[joint orientation session](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
proves [R85](GEOMETRIC_BAND_ORIENTATION_DEPTH.md), with log-times-log-log
depth for arbitrary nonlinear orientations of shared unordered depth
magnitudes. Geometric-band unloading exposes the clean output space
needed for local orientation conversion. This bypasses the standalone
permutation interface and retained R79's general frontier at that
checkpoint. The later
[general-magnitude session](GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
proves [R86](INTERVAL_DILATION_DEPTH.md): for fixed $`L\ge6`$,
$`D_T=O_L(N/b^2+\ell^{3/2})`$, where $`\ell=\log_2(n+2)`$,
for arbitrary node-dependent magnitudes, signs, zeros, and prescribed
complex leaf phases. One clean qubit, the original literal widths,
$`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and the full initialized
contract are retained. The subsequent
[factorization session](INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
proves [R87](COIN_DILATION_DEPTH.md), improving the general depth to
$`O_L(N/b^2+\ell\log_2(\ell+2))`$ through product-coin spectral
readout. It matches R85's order without that magnitude promise. R79
remains the absolute all-precision theorem. This map preserves the
earlier preparation history below.

The later [terminal refinement R88](PREFIX_SELECTED_COIN_TERMINAL.md)
handles arbitrary real deepest-M layers selected by a preserved external
prefix, with one clean qubit, the original real width, the same count
budgets, and full initialized-isometry error. For every $`0\le M\le n`$,
its depth is $`O_L(N/b^2+\ell\log_2\log_2(M+4))`$.
A logarithmic terminal contributes
$`O_L(\ell\log_2\log_2(\ell+4))`$, and the
[contiguous coin layout](SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
now routes in $`O(\ell)`$ depth. These component gains initially left
R87's $`O_L(N/b^2+\ell\log_2(\ell+2))`$ complete bound unchanged.

The subsequent [R89 theorem](CAUSAL_COMPACT_COIN_DEPTH.md) removes the
serial coarse search through compact fixed coins, causal local matchings,
and parallel reference transforms. Its prefix-selected recursion gives
the general fixed-accuracy bound at that checkpoint

```math
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(n+2))\right),
```

with the same one-clean allocation, real/complex literal widths, count
budgets, and full initialized-isometry contract. Here $`\log_2^*x`$
counts repeated base-two logarithms until the value is at most 2.
For a real terminal of height M, R89 gives
$`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$. R87's fine decoder and
R88's sharp-query and fallback interfaces remain proof dependencies.

The [whole-terminal session](TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
then proves [R90](DIRTY_CATALOGUE_TERMINAL.md), with general fixed-L
depth $`O_L(N/b^2+\ell)`$ and the same one-clean allocation, literal
widths, count budgets, and full operator contract. Three R89 coarse
reductions leave a small terminal whose complete rounded program selects
fixed catalogue factors. Exact dirty program queries and XOR slot
translations retain the unknown program input; a commutator on two
disjoint borrowed arrays cancels their unwanted actions. The proof
charges native local synthesis, selected-factor error, catalogue volume,
and the original-width splice. R79 remains the absolute all-precision
theorem; optimal general depth remains open.

The [precision-terminal follow-up](PRECISION_TERMINAL_SESSION_RESULT.md)
proves R91's parallel word query, R92's absolute catalogue envelope, and
R93's growing-precision group theorem. At inverse-polynomial accuracy in n,
R93 improves the square-root-width bound to $`O(\ell^2/\log\ell)`$.
Its next bottleneck is the number of long-prefix group reads. Larger
groups must improve both local native execution and the initialized
matrix-unit labels or the retained bottom schedule. These are costs of
the proved construction, not unrestricted lower bounds.

This is a research decision aid, not a new theorem or a manuscript draft.
It reconstructs the pre-session upper-bound costs, isolates the lower-bound
proof gap, and records what a new mechanism must change. The audit starts
from main commit `b5d22d8`; linked chapters remain the proof homes.

## 1. A diagnostic benchmark and its fixed contract

Use fixed error $`\eta_0=1/64`$, $`L=6`$, $`N=2^n`$, one external
clean qubit, and sufficient $`b=\lceil C_b\sqrt N\rceil`$ arbitrary
dirty qubits, for a fixed constant $`C_b>0`$ and sufficiently large n.
The initial target is the **entire complete real-frame family**, including
unequal, signed, and singular parameter tuples. Seek a better worst-case
T-depth while retaining $`T=O(\sqrt N)`$ and $`G=O(N)`$ on the same
circuit. This benchmark is a lens inside the all-precision real/complex
project, not a replacement for its scope.

The [one-clean theorem R72](ONE_CLEAN_DEPTH.md) gives
$`D_T=O(\sqrt n\log(n+2))`$ here. The unrestricted lower bound is
only $`\Omega(1)`$. Count optimality and the existence of sublinear
depth are already proved. A new constant-depth exclusion alone would not
close this asymptotic gap.

The immutable task contract is

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\mathrm{op}}\le\eta_0.
```

It covers every logical input, dirty/reference correlations, literal
phase, clean leakage, and work return. Native Clifford+T gates, actual
inverses, sources, queries, reflections, guards, and precision are charged.
There is no measurement, reset, QRAM, supplied precision state, or free
coherent evaluator. T-depth permits arbitrary Clifford interlayers, whose
count is still charged. Certified classical preprocessing terminates;
efficient preprocessing and a general native emitter are not established.

The tag basis, entry-indexed factorization, amplitude-amplification word,
early/tail split, and choice of lower-bound witness are **method choices**.
They may change while the full contract remains fixed.

<a id="2-where-the-current-upper-bound-spends-depth"></a>

## 2. Where the pre-R78 upper bound spent depth

The applicable branch is [R66 §10](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision),
as adapted by [R72 §5](ONE_CLEAN_DEPTH.md#5-refined-circuits-above-a-new-width-boundary).
Write $`\ell=\log_2(n+2)`$. It protects
$`M=\lceil A[L+\log_2(n+L+2)]\rceil=\Theta(\log n)`$
terminal bits and compiles $`h=n-M=\Theta(n)`$ early levels together.
This was the audited split; R66's earlier fixed-accuracy exposition used a
different, larger protected suffix.

| Component | Fixed-accuracy allowance | Multiplicity and reason |
|---|---|---|
| Complete-column recovery | $`\nu=2\lceil\pi\sqrt K/4\rceil+1=\Theta(\sqrt n)`$ source/inverse calls, where $`K=2^{\lceil\log_2(h+1)\rceil}`$ | The undamped row/source overlap is $`W_h/\sqrt K`$; this is the multiplicity of the chosen damped recovery word |
| Large-address source queries and source reflections | $`O(\log n)`$ per call | Boolean phase, calibration, marker extraction, and source reflections use tables with $`O(n)`$ address bits |
| Arithmetic, short coefficient/label banks, and outer work reflections | $`O(\log\log n)`$ per call | Precision $`r=L+\lceil\log_2\nu\rceil+16=\Theta(\log n)`$; small predicates and arithmetic use logarithmic depth in r |
| Shared gradient preparation and actual return | $`O(L+\log r)=O(\log\log n)`$ total | Prepared once outside recovery; the recovery reflection excludes the live gradient |
| Terminal stream and deepest one-clean layers | $`O(\ell^2)`$ total | A logarithmic number of remaining layers, with their own precision/query schedule; no early recovery multiplier |
| Width penalty | $`O(N/b^2)=O(1)`$ | Already constant-order at this benchmark |
| Optional complex phase dressing | $`O(\log n)`$ additional | The one-clean literal diagonal and half-error composition retain the shared work pool |

The resulting ledger, with subordinate terms shown, is

```math
D_T=O\!\left(\sqrt n\log n+\sqrt n\log\log n+
\log^2 n+N/b^2\right).
```

These are **upper allowances**, not lower bounds on their components.
The leading term comes from large queries repeated during early recovery.
Repeated precision-state preparation has already been removed by
[shared-gradient recovery](SHARED_GRADIENT_RECOVERY.md#7-resources-and-the-uniform-precision-range).
Its preparation error is charged twice at the boundaries. Source/damping
errors instead telescope over $`\nu`$ calls, using
$`\delta=\eta/(64\nu)`$ and a full-unitary comparison within the
exactly invariant gradient sector. The exterior preparation hybrid
separately handles leakage into other gradient sectors.

The width/count proof matters when changing this schedule. The early
data tables are suppressed by $`2^{-M/2}`$, whereas standalone arithmetic
tables are not. At fixed accuracy the latter have polynomial size in n
and fit the exponential budgets above an absolute dimension cutoff.
R72 fits early/stream work below fixed fractions of
$`b_c=\sqrt{NL/[n(L+n)]}=\Theta(\sqrt N/n)`$; sequential work
reservations take a maximum, not a sum. It separately reserves the bottom
one-clean core and checks the shallow-query width threshold. Small
dimensions use the proved fallback. These analytic cutoffs are not
evaluated practical crossover estimates.

**Residual-cost check.** Improving the early encoding alone leaves the
$`O(\log^2 n)`$ terminal schedule. Improving large-query depth alone
leaves repeated arithmetic and the current recovery multiplicity.
A sublogarithmic real-frame result would also need new phase handling
before the present complex extension could inherit that order.

### The sharp barrier applies to an encoding, not to all circuits

[R67](ENTRY_INDEXED_ENCODING_LIMIT.md#4-robustness-to-operator-norm-approximation)
considers row/column preparations sharing an entry-indexed basis
$`(i,j,k)`$, with each prepared row and column having norm at most one.
For the balanced frame, represented-operator error at most
$`\epsilon\lt1`$ implies

```math
\alpha\ge(1-\epsilon)\sqrt{h+1},\qquad
\|\alpha B-W_h\|_{\mathrm{op}}\le\epsilon.
```

Reweighting and adding entry-local interference channels cannot improve
this normalization order. For exact encoding of the balanced instance,
the pre-R78 tag construction is within a factor below $`\sqrt2`$ of
the architecture's optimum. R67 supplies **neither a
recovery-query lower bound nor an unrestricted native-depth lower bound**.
A coherent matrix-basis change, non-entry-indexed preparation, different
factorization, direct native circuit, or different recovery remains eligible;
its entire native implementation must be priced.

## 3. The lower-bound route and its exact missing step

The R75–R77 route uses
$`W_n=I_{2^{n-1}}\otimes R_y(\theta)`$ with
$`\theta=\tfrac12\arccos(1/3)`$. The following stages have distinct scope.

| Stage | Status |
|---|---|
| Terminal clean return | [R75 §6](FOUR_LAYER_OPERATOR_RIGIDITY.md#6-one-clean-qubit-terminal-sectors-only): initialized error below $`\tau=\sin(\pi/8)/(2\sqrt2)`$ forces an at-most-four-layer circuit to preserve the clean sectors exactly. The unused sector's target action remains unspecified. |
| Correlated midpoint compression | Split into two native encoders with at most two layers each. They send the same input clean Pauli to S and the same data/dirty Pauli to A and B; $`Q=(I+S)/2`$ and $`\|(A-B)Q\|\le2\eta`$. These are actual shared-input encoder images, not arbitrary commuting involutions. |
| Freeze compressed conjugations | **Open for unrestricted native non-Clifford midpoints.** For every test Pauli R centralizing the first outer axis space, small occupied error must force $`Q(ARA-BRB)Q=0`$. R76 proves this for Clifford S; R77 handles bounded inner clean spread. |
| Lift the compressed commutant | [R76](ONE_CLEAN_MIDPOINT_RIGIDITY.md#lifting-the-compressed-commutant): after freezing, $`BAQ=D_0Q=QD_0`$ for an abelian-algebra unitary $`D_0`$. This already allows non-Clifford S and growing spread. |
| Discrete occupied phases | [R76](ONE_CLEAN_MIDPOINT_RIGIDITY.md#discrete-occupied-phases): native entry arithmetic forces occupied scalars to eighth roots. Small enough error forces identity; target symmetries then exclude the witness. No additional spread restriction enters here. |

For Clifford S, compressed conjugation expands into eight Clifford
boundaries. [R74's finite-boundary lemma](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
gives a positive gap. For non-Clifford S, R77 expands the inner core into
$`2^t`$ Pauli terms and obtains $`K_t=2^{t+3}`$ equal-coefficient
boundaries. Its gap $`\kappa_t=2^{-\mathcal M_t/2-1}`$ tends to zero
as t grows. The missing estimate must control cancellation in this
**correlated compressed sum** uniformly in both encoder spreads.
Renaming t as support size, boundary count, or crossing rank does not
remove the unbounded quantity. Pairwise full-space separation alone does
not provide that estimate.

The sufficient invariant is a width-independent $`c>0`$ such that, for
native $`D_T(V)\le4`$, the initialized projector $`\Pi`$, and a
data/dirty Pauli P,

```math
[V,\Pi]=0,\qquad \|[V,\Pi P]\|\lt c
\quad\Longrightarrow\quad [V,\Pi P]=0.
```

It is equivalent to a uniform occupied-image gap **on actual shared-input
native encoder pairs**. It is a sufficient statement covering more circuits
than the witness-approximation problem; no converse is proved. A counterexample
must give native exact-return four-layer circuits with
$`0\lt\|[V_j,\Pi P_j]\|\to0`$. Such a counterexample would still not
be a compiler for the witness.

[R77's common-flip example](ONE_CLEAN_INNER_SPREAD.md#5-a-common-compatible-flip-need-not-exist)
disproves a universal sign-compatible sector-flip shortcut. Its occupied
defect is the fixed constant $`\sqrt{1-1/\sqrt2}`$; padding does not
shrink operator norm. Neither that example nor a small scalar trace deficit
disproves the positive operator gap. The necessary growth conditions
$`r_{\min}=\Omega(\sqrt L)`$ and
$`t_{\min}=\Omega(\log\log L)`$ concern increasing precision L,
not increasing n at fixed accuracy. Here $`r_{\min}`$ minimizes over
native four-layer representations, and $`t_{\min}`$ minimizes both encoder
spreads over all compatible representations and splits.

### Why this witness cannot establish dimension-growing depth

This is a direct consequence of the existing compiler, not a new theorem.
At fixed error, apply the [one-clean baseline](../../docs/ONE_CLEAN_COMPILER.md#6-complete-real-frames-and-the-grouped-extension)
to the active qubit alone: $`n'=1`$, $`b_0=L+8`$,
$`T,G=O(L)`$, hence $`D_T\le T=O_\eta(1)`$.
Tensor identity on the other logical wires and all unused dirty wires.
The initialized-isometry and literal-phase guarantees survive, and
$`b_0`$ fits $`b=\Theta(\sqrt N)`$ for sufficiently large n.
At fixed L, these costs also satisfy the benchmark's $`T=O(\sqrt N)`$
and $`G=O(N)`$ budgets.

Consequently, this spectator-padded witness cannot prove depth growing
with n at fixed accuracy, even if stronger fixed-layer rigidity is proved.
Such rigidity can address precision dependence or provide a reusable
invariant. A dimension-growth objective needs a genuinely
dimension-dependent target family and an invariant that detects it.

## 4. Previous routes and what would justify reopening them

These are restrictions of recorded mechanisms, not universal impossibility
claims. The [research index](../README.md) retains all positive cases and
proof/evidence links.

| Route | Surviving limitation | A meaningful change |
|---|---|---|
| [Historical endpoint reductions](../ROUTE_HISTORY.md#the-earlier-sufficient-construction-target) | The old count target is resolved; an unpaid encoder in an alternative construction is not a hole in Result C | A new simultaneous depth ledger for the complete family |
| [Majorana source resynthesis](SOURCE_T_DEPTH.md#4-scope-and-finite-evidence) and [two-layer source gaps](SHALLOW_SOURCE_OBSTRUCTION.md#5-scope) | Growing depth bounds restrict the source architecture; the unrestricted approximate exclusion is only two layers | Leave that architecture, or find an invariant of the whole unrestricted compiler; primitive lower bounds cannot simply be added |
| [Leakage diagnostics](SOURCE_LEAKAGE_DIAGNOSTICS.md#4-what-this-closes-and-what-it-leaves-open), [echoes](HOPF_FLAG_ECHO.md#4-what-a-compiler-improvement-would-still-require), [radial filtering](HOPF_RADIAL_FILTER.md#6-resource-frontier-and-evidence-boundary) | Ideal angular stability misses physical leakage; generic echoes retain first-order error; the successful filter retains the old depth order | A full-isometry error mechanism that lowers charged native depth |
| [Common-source hoisting](COMMON_SOURCE_REUSE.md#11-a-common-source-identity-and-the-remaining-reflection) | Conjugated reflections remain; cached monitors become stale and one-use banks leave rejection | Price a shallow reflection, dynamically valid monitor, or a different complete return identity |
| [Protected sources](PROTECTED_UNARY_SOURCE.md#4-source-error-depth-and-the-remaining-complete-frame-ledger) and [windowed predicates](WINDOWED_GROUP_PREDICATES.md#4-predicate-depth-counts-and-the-complete-frame-scope) | Their grouped schedule retains separate linear logical and query allowances; current column compression already bypasses it | Change both surviving allowances or leave the sequential group interface |
| [Retained-source fusion](RETAINED_SOURCE_FUSION.md#charge-the-program-workspace-and-remaining-transport) | Stabilizer completion is batched, but the ordered transport and changed-basis selectors remain | A generic shallow transport factor with all queries and return charged |
| [Shared-prefix caching and differences](SHARED_PREFIX_QUERY_AUDIT.md#8-a-correction-needs-a-proved-structural-saving) | An arbitrary row difference can encode a fresh lookup; promised cheap refreshes do not cover generic tables | A proved generic saving or a changed, explicitly paid reader interface |
| [Dirty-query scheduling](NONUNIFORM_DIRTY_INDICATOR.md#6-the-tail-improves-the-complete-frame-frontier-is-unchanged) | The width penalty is already constant here; the logarithmic local query bound is only an upper bound | Improve queries in the current complete recovery schedule, or change global query reuse |
| [Precision streams](SHARED_BINARY_TERMINAL.md#8-bottom-and-terminal-composition) and [one-clean relocation](ONE_CLEAN_DEPTH.md#5-refined-circuits-above-a-new-width-boundary) | These constructions are completed; terminal savings and high-precision splices do not remove fixed-accuracy early recovery | Change early recovery or the full large-width schedule |
| [Entry-indexed columns](ENTRY_INDEXED_ENCODING_LIMIT.md#5-relation-to-the-current-tag-construction-and-the-remaining-freedom) | Reweighting and entry-local channels cannot improve normalization order | Change the support/factorization contract or replace recovery, with a full native ledger |
| [Arbitrary branch orientations](BRANCH_ORIENTATION_SESSION_RESULT.md#7-remaining-mechanism-and-stopping-point) | A standalone dirty-output conversion was unpaid; a long initially protected output retained square-logarithmic depth | [R85](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) bypasses that interface with geometric-band unloading; [R89](CAUSAL_COMPACT_COIN_DEPTH.md) later gives logarithmic early coarse recovery for arbitrary magnitudes through causal batching |
| [In-place conversion checkpoint](PERMUTATION_BRIDGE_SESSION_RESULT.md) | Unary relabeling needs a recursively selected permutation; four-shear and exact Clifford-encoding shortcuts have scoped obstructions, with no unrestricted native conclusion | Discover a conversion or joint recovery mechanism beyond these restrictions, then establish its full resource and return proof |
| [Midpoint rigidity](ONE_CLEAN_MIDPOINT_RIGIDITY.md#6-the-remaining-one-clean-question) and [inner spread](ONE_CLEAN_INNER_SPREAD.md#4-precision-and-the-four-layer-sizes) | Unbounded spread escapes present constants; a universal common flip is false; the witness has fixed-accuracy constant depth | Uniform native compressed freezing, a native vanishing-defect counterexample, or a dimension-sensitive witness/invariant |

## 5. What the next research plan must decide

The [joint orientation/recovery prompt](JOINT_ORIENTATION_RECOVERY_PROMPT.md)
has been executed. Its [outcome](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
is R85's complete arbitrary-orientation construction under shared unordered
depth magnitudes, with log-times-log-log depth and the retained budgets.

The [permutation checkpoint](PERMUTATION_BRIDGE_SESSION_RESULT.md#8-decision-and-verification-boundary)
ended that bounded conversion search without a general improvement.
R85 subsequently changed the complete source/recovery interface. The
[general-magnitude design](GENERAL_MAGNITUDE_RECOVERY_PROMPT.md) has also
been executed; its [outcome](GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
is R86's complete local frame compiler through interval dilation.
The balanced-Haar conversion has $`O_L(\sqrt\ell)`$ normalization,
and paid recovery yields the $`O_L(\ell^{3/2})`$ general depth term.
Its [routing support](SELECTED_DIRTY_ROUTING.md) prices selected physical
wires and variable-size query families. The subsequent
[factorization design](INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md) has
also been executed. Its [outcome](INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
is R87's product-coin readout: a paid embedding, tensor Hadamards,
and a last-nonzero-block XOR erase the fine tag with the actual source
inverse. The old interval-overlap amplification is absent.

The [continuation](COIN_CONTINUATION_SESSION_RESULT.md) first improves
the terminal and contiguous router, then proves R89's complete bound.
[Compact fixed coins](COMPACT_COIN_FLAGS.md) preserve every coarse flag;
[subset primitives](SHORT_LABEL_SUBSET_ROUTING.md) supply selected
compressions and matchings. A causal batched word realizes the embedding
without treating an arbitrary global weight-preserving permutation as
flag-preserving. Parallel reference sources include their entire
initialized work in the nonzero-block test. Prefix recursion applies
the resulting logarithmic early decoder to successively smaller tails.

The [whole-terminal session](TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
subsequently replaces the unbounded recursive terminal with R90's dirty
catalogue commutator. This reaches logarithmic general depth at fixed
accuracy; it does not prove that depth optimal. R92 and R93 subsequently
give absolute growing-precision refinements. A next plan must target
optimal general depth, the stronger growing-precision benchmark, or a
dimension-growing lower bound. Complete native resource and return
proofs remain the criterion. Proved affine, corner, and R85–R93 circuits
should not be repeated as unresolved tasks. R79 remains the absolute
general all-precision bound.

The fixed-L constants and unbounded classical preprocessing in R90 must
be retained when considering shrinking error or an effective compiler.
No practical crossover or scalable native emitter follows from the
finite catalogue. A further constructive gain needs a new complete
operator identity and resource proof.
A lower-bound attack on dimension growth first needs a different witness. The next
one-clean constant-depth obstruction is a separate valid target, with the
freezing step above precisely localized; it should be selected for that
purpose rather than described as closing the dimension-growth gap.

At the start of a bounded exploration, state the question, the assumption
or representation it changes, and the uncertainty to resolve. Before
claiming an improvement or expanding finite validation, supply the exact
native identity or encoding, its domain of valid inputs, and symbolic
depth/count/width/error/return recurrences. Price basis changes,
all table access, actual inverses, recovery multiplicity, inactive occupied
work, and the terminal tail. For a lower-bound candidate, state all
quantifiers, the norm, native-realizability conditions, and whether its
consequence grows with n, with L, or only strengthens a constant.
Stop if the proposed gain leaves the old leading order, imports unpriced
clean history, or fails work return. A finite example should test a
specific surviving premise after this analysis, not substitute for it.
A pause records what was learned and the next unresolved question;
failure of that candidate does not end the general research objective.

The preparation audit found no defect in the inspected resource composition
or the R75–R77 proof chain. It reconstructed their interfaces and the
critical imported lemmas; it did not independently reprove every synthesis
premise, construct a general emitter, or establish a new asymptotic bound.
The [research prompt](RESEARCH_SESSION_PROMPT.md) designed the subsequent
attack from this map. Its preferred hypothesis was native depth-label erasure
through nested subtree projectors, with an independent query-cost reserve
and explicit qualification gates. The [execution result](RESEARCH_SESSION_RESULT.md)
records the proved construction and the reserve routes that did not qualify.
