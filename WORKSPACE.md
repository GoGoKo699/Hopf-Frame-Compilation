# Continuing research workspace

This is the entry point when a previous conversation or execution workspace
is unavailable. Proofs and decisions live in the repository.

The 2026-10-03 uniform-precision pass starts from verified main
`9ae9685154b498e3350f8a89262dfc1f02fb07d6`, after the width-sensitive
blocked bilinear theorem. Check later commits before
continuing.
The selected state-based Hopf QBP construction and its bounded-input audit are complete; the
[consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) is their entry point.
It prepares a state and changes the gradient decoder while retaining all
raw coordinates, singular angles, and both complex-gradient streams.
The active scientific targets are the constant-clean complete-frame
endpoint and optimal T-depth. Application-level advantage is outside the
current research scope by the author's decision; retain the existing cost
comparisons and classical baselines as boundaries, not a pending task.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized. Keep manuscript writing
and release work outside this research pass. Use small analytic examples
and finite checks, without large simulations, QRAM, resets inside a
compiler execution, supplied catalysts, or hidden initialized work.

1. Begin active depth work with the
   [uniform-precision theorem](docs/UNIFORM_PRECISION_DEPTH.md),
   then the [width-sensitive bilinear theorem](docs/BLOCKED_BILINEAR_LOOKUP.md),
   then the [unary phase-source theorem](docs/UNARY_PHASE_GRADIENT.md#7-complete-frame-theorem-at-fixed-accuracy),
   then the [grouped complete-frame theorem](docs/GROUPED_PROGRAM_PREFETCH.md#8-complete-frame-theorem-at-fixed-accuracy),
   [chunked dirty indicator](docs/CHUNKED_DIRTY_INDICATOR.md),
   [conditional geometric source](docs/CONDITIONAL_GEOMETRIC_SOURCE.md),
   [two-layer obstruction](docs/SHALLOW_SOURCE_OBSTRUCTION.md),
   [radial filter](docs/HOPF_RADIAL_FILTER.md),
   [short-echo audit](docs/HOPF_FLAG_ECHO.md), and
   [error geometry](docs/HOPF_ERROR_ACCUMULATION.md),
   [masked dirty sums](docs/DIRTY_SUM_COMPRESSION.md)
   and the [hybrid tradeoff](docs/PARALLEL_DIRTY_LOOKUP.md), then the
   [amortized baseline](docs/AMORTIZED_DIRTY_LOOKUP.md),
   then its [batched allocation](docs/BATCHED_DIRTY_LOOKUP.md),
   [depth proof](docs/T_DEPTH_COMPILER.md), and
   [source-depth certificate](docs/SOURCE_T_DEPTH.md), then the
   [routed indicator](docs/PARALLEL_DIRTY_LOOKUP.md). The completed
   [state-based QBP theorem](docs/STATE_BASED_QBP_THEOREM.md) supplies its
   separate task contract, resource table, proof map, and limits.
2. Use the [cost comparison](docs/QBP_COST_COMPARISON.md) and
   [bounded-input audit](docs/BOUNDED_INPUT_QBP.md) for precision regimes,
   classical construction, and explicit Pauli baselines. The
   [residual coefficient proof](docs/RESIDUAL_TABLE_PREPROCESSING.md)
   documents the classical interval helper; the
   [native rows and tables](docs/NATIVE_RESIDUAL_ROTATION.md) connect it to
   exact masks, borrowed-signal rotations, and enabled two- and four-row tables.
   The [bounded state integration](docs/NATIVE_RESIDUAL_STATE.md) composes
   two tables and amplification with the actual inverse for one or two
   system qubits, including the enlarged initial reflection.
   The [coherent selector](docs/NATIVE_RESIDUAL_BRANCH.md) adds an arbitrary
   protocol branch to the one-system-qubit residual preparation.
   The [native residual QBP integration](docs/NATIVE_RESIDUAL_QBP.md)
   connects that selector to the actual coarse circuit, controlled
   observable, and both raw gradient decoders for a fixed complex target.
   The [coverage map](docs/VERIFICATION.md#state-based-qbp-coverage) is the
   current entry point for what is proved, implemented, and optional.
3. For the separate frame question, read the
   [research status](docs/OPEN_PROBLEM.md),
   [grouped compiler](docs/CONDITIONAL_SUFFIX_COMPILER.md), and
   [retained stopping rule](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
4. Use [verification](docs/VERIFICATION.md), [attribution](docs/SOURCE_MAP.md),
   and [related work](docs/RELATED_WORK.md) before extending a claim.
   [README](README.md) and [REVIEW](REVIEW.md) retain the established frame
   results; the existing [publication scope](manuscript/PUBLICATION_SCOPE.md)
   is unchanged.

## Current results and their proof homes

| Result | Status and proof |
|---|---|
| Prescribed complete frame, exact gates | Matching size and CNOT/depth tradeoffs for every clean-work budget; [exact theorem](docs/COMPILER_THEOREM.md) |
| Prescribed complete frame, Clifford+T | Matching T-count under the sufficient-clean reservation; [fault-tolerant theorem](docs/FAULT_TOLERANT_COMPILER.md). Smaller clean budgets retain the [borrowed](docs/BORROWED_WORKSPACE_COMPILER.md) and [one-clean grouped](docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) bounds |
| State-based real and complex QBP | Two compiler flags, fine state preparation, an exact-return coarse word, and corrected magnitude/phase decoders; [consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) |
| Additional dirty banks | Improve the state preparation T bound while charging the exact coarse circuit; [banked proof](docs/COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) |
| State-based T-depth | Two complete schedules, one retaining the sharper count at a stronger dirty reservation; [depth proof](docs/STATE_QBP_DEPTH.md) and [fair comparison](docs/QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison) |
| Complete-frame T-count and T-depth | Uniform depth O(NL/b²+nL), retaining the same-circuit count bound at b at least 17(L+n+7); matching through sqrt(N/n); [uniform theorem](docs/UNIFORM_PRECISION_DEPTH.md) |
| Slowly growing precision | For 6 ≤ L ≤ log₂(n+2)/16, depth O(NL/b²+n) and count O(sqrt(NL)+NL/b), with absolute constants; matching through sqrt(NL/n); [uniform theorem](docs/UNIFORM_PRECISION_DEPTH.md) |
| Fixed-accuracy width-dependent depth | Two external flags, T-count O(sqrt N+N/b), and depth O(N/b²+n) at b at least 17(L+n+7); both orders match through sqrt(N/n); [blocked bilinear theorem](docs/BLOCKED_BILINEAR_LOOKUP.md) |
| Complete-frame error accumulation | Sharp ideal-angle stability and finite relative spectra; coherent linear leakage in actual shared-flag source layers; [scoped error audit](docs/HOPF_ERROR_ACCUMULATION.md) |
| Filtered complete-frame source | Quadratic radial error on the same two flags, charged native selective phases, and a smaller source-precision cap; [filter proof](docs/HOPF_RADIAL_FILTER.md) |
| Conditional precision depth | Logarithmic source/reflection depth using an active zero suffix and two external flags; [conditional source](docs/CONDITIONAL_GEOMETRIC_SOURCE.md). Query/predicate depth remains charged |
| Unrestricted two-layer source obstruction | Robust constant error floor with arbitrary Clifford interlayers and dirty helpers; [transfer proof](docs/SHALLOW_SOURCE_OBSTRUCTION.md). No initialized-clean frame lower bound |
| Bounded-input construction | Polynomial construction for the listed grouped/state alternatives, explicit program output, and separate fine-search caveats; [computational audit](docs/BOUNDED_INPUT_QBP.md) |
| Implemented evidence | Certified residual rows, bounded preparations, and [complete bounded residual QBP streams](docs/NATIVE_RESIDUAL_QBP.md), alongside the earlier exact-target examples; [claim-to-proof coverage](docs/VERIFICATION.md#state-based-qbp-coverage) |

The theorem's original accuracy bits K and state precision
$`P=\max(n,K)`$ must remain distinct. The observable needs accuracy K,
not the state's dimension floor. Compiler flags, the interference branch,
initialized system, and observable work are separate reservations. The
coarse word returns dirty work exactly; fine preparation includes all work
return and leakage in its initialized-isometry error.

At fixed accuracy the original frame route retains the better available
bound. Additional banks improve the state-only T expression in the stated
high-precision regimes, but sampling still grows as $`4^K`$ up to its
confidence factor. Clifford order and classical decoding order do not
improve. Explicit Pauli inputs also admit deterministic classical gradients
and term-only classical sampling. These results establish a task-specific
compiler improvement, not a general end-to-end gradient speedup.

## Current depth frontier

The [uniform-precision theorem](docs/UNIFORM_PRECISION_DEPTH.md) gives,
for every $`L\ge6`$, two clean flags, and $`b\ge17(L+n+7)`$,
one prescribed complete real-frame circuit with absolute constants:

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

When nonempty, $`17(L+n+7)\le b\le\sqrt{N/n}`$ is a simultaneous
matching interval, with $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$. For $`L=\Theta(n)`$, or
inverse-polynomial error in N, sufficiently large $`b=\Theta(n)`$
gives $`T^\star=\Theta(N)`$ and $`D_T^\star=\Theta(N/n)`$.
Outside the matching interval, retain the nL count term; the sufficient
condition $`L\le N/n^2`$ makes it absorb into $`\sqrt{NL}`$.

There is a stronger uniform corollary for slowly growing precision:

```math
6\le L\le\frac{\log_2(n+2)}{16},\qquad
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+n\right),\qquad G=O(NL).
```

Its matching interval extends through $`b\le\sqrt{NL/n}`$.
The absolute constants and thresholds are independent of eta and L.
For fixed L this recovers the [blocked bilinear curve](docs/BLOCKED_BILINEAR_LOOKUP.md),
whose constants could depend on eta. At square-root-scale width the
fixed-accuracy depth upper bound is O(n), while the lower bound remains
Omega(1). No unrestricted depth optimality or endpoint improvement follows.

The new proof uses the charged [unary source](docs/UNARY_PHASE_GRADIENT.md)
in early groups, with cutoff sublinear in n uniformly in the displayed
precision range. Rectangular indicator allocation keeps the tail's
weighted count at $`O(\sqrt{NL})`$ and its entire depth at
$`O(NL/b^2+n)`$. For larger L, the existing
[hybrid theorem](docs/PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision)
already absorbs its $`n\log_2(n+2)`$ term into nL. No new native
query identity is needed. Source preparation and actual inverses are
charged; all dirty work returns, including arbitrary reference inputs.
The original sufficient threshold remains exactly $`17(L+n+7)`$.
The two-clean high-precision endpoint remains outside that threshold.

Older schedules remain useful when their count or workspace requirements
are sharper. In those ledgers, $`\chi(t)=\log_2(t+2)`$. The [amortized compiler](docs/AMORTIZED_DIRTY_LOOKUP.md)
removed the per-batch logarithms but retained quadratic query depth.
Its capped precision allocation has total source depth
$`O(nL+n\log(n+1))`$ and is within $`5n`$ bits of its additive
error-certificate optimum. The dirty-counter hybrid, grouped program
reuse, unary source, and blocked bilinear query successively improve the
depth upper bound. Their scoped obstructions are not frame-depth lower
bounds. The separate state-based and complex-frame schedules retain
their existing contracts.

The [state-based depth theorem](docs/STATE_QBP_DEPTH.md) closes the selected
composition audit for both real and gauge-fixed complex Hopf QBP. Put
$`B_0=P+n+7`$. With the same two compiler flags, its schedules give:

| Dirty reservation | Compiler T-depth | Compiler T-count |
|---|---|---|
| $`b\ge2B_0`$ | $`O(NP/b+P+n^3)`$ | $`O(NP)`$ |
| $`b\ge16(B_0+\sqrt{NP})`$ | $`O(P+n^3)`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ |

Both have $`G=O(NP)`$ and charge the actual coarse C, its inverse in
magnitude readout, the fine residual table, occupied flags, predicates,
and all returned work. Each row's bounds hold for the same circuit.
The [complete task ledger](docs/QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison)
adds observable depth at precision K and all S or $`2S`$ executions.
A sufficiently large common pool permits an improved real-chart depth
upper expression even where the original and state T-count orders agree.
The original real schedules and eligible complex alternatives remain
in that comparison.

The routed-indicator refinement conjugates a single X by the existing
exact bank router. Its one-hot XOR uses no extra scratch and has linear,
rather than quadratic, address T-depth. This improves the state B bound
from $`O(P+n^4)`$ to $`O(P+n^3)`$ at the same sufficient reservation.
The complete real-frame schedules also use the established two-dirty-helper
MCX construction while query storage is idle. They now give

```math
D_T=O\!\left(\frac{NL}{b}
+\min\{nL+n^2,L\ell_*(n)+n^3\}\right),\qquad T,G=O(NL),
```

at $`b\ge2(L+n+7)`$. With the stronger sufficient pool
$`b\ge C(L+n+7+\sqrt{NL})`$, the same minimum depth expression
holds without the $`NL/b`$ term and with the sharper count
$`T=O(\sqrt{NL}+L\ell_*(n))`$. At fixed L this is
$`T=O(\sqrt N)`$ and $`D_T=O(n^2)`$ in one circuit, improving the
previous cubic depth bound. Literal routed-indicator/query checks include
arbitrary dirty inputs, exact phases, and the actual inverse orientation.

These remain upper schedules. The [source-depth certificate](docs/SOURCE_T_DEPTH.md)
parallelizes the paired source's two tails and proves matching exact depths
within the specified Majorana-layer architecture. This changes constants,
not the asymptotic Hopf bounds. It is not a lower bound for arbitrary
Clifford interlayers, controlled sources, or approximate replacement sources.
The [current literature audit](docs/RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026)
records why shallow clean-workspace synthesis and supplied catalysts do not
directly replace the dirty source. Large-width depth, precision regimes
outside the proved matching window, and the selected $`b=N+n+7,L=N`$
complete-frame endpoint remain open. Obtaining sublinear
depth for these exact uncontrolled sources requires leaving the certified
architecture or changing the target. Do not repeat tail rescheduling or
denominator checks as an unrestricted lower bound.

## Completed native residual components

The bounded implementation sequence is complete. Detailed identities,
literal gate counts, certificates, and test scopes remain in their
canonical chapters; the handoff does not duplicate those ledgers.

| Component | Canonical implementation and scope |
|---|---|
| Certified coefficient programming and borrowed-signal rotations | [Native residual rows](docs/NATIVE_RESIDUAL_ROTATION.md): rational intervals to literal native words with full-operator error |
| Enabled two- and four-row tables | [Native table proof](docs/NATIVE_RESIDUAL_ROTATION.md): exact inactive identity, literal address phases, and no extra helper at these arities |
| One- and two-system-qubit preparation | [State integration](docs/NATIVE_RESIDUAL_STATE.md): two flags, actual-inverse amplification, all dirty-input return error, and the charged 28-T enlarged reflection |
| Coherent reference/target selection | [Branch integration](docs/NATIVE_RESIDUAL_BRANCH.md): arbitrary branch, preserved relative phase, and branch-independent initial reflection |
| Both complete gradient streams | [Residual QBP integration](docs/NATIVE_RESIDUAL_QBP.md): certified complex one-qubit target, actual coarse/inverse and controlled observable, all-outcome readout, and exact rational histogram decoders |

The small native preparations use q+2 arbitrary dirty wires. At q=L+10,
this exceeds the basic n=1 allocation by four wires and the basic n=2
allocation by three; these emitters fit the banked pool and do not replace
the minimum-budget fallbacks. All initialized system, compiler-flag, and
protocol-branch inputs stay separately counted. No intermediate reset,
postselection, or ideal-inverse substitution is allowed.

The [coverage audit](docs/VERIFICATION.md#state-based-qbp-coverage)
separates uniform proofs from the numerical and exact finite evidence.
In particular, exact rational decoders for a fixed fixture do not turn
the general floating-point contraction utilities into certified decoders.
The existing bounded integration supplies neither a same-accuracy native
fine-frame benchmark nor an end-to-end advantage claim.

## Remaining work and continuation criteria

The task-specific theorem, real/complex decoding, quantum resource ledger,
bounded-input construction, selected depth composition, and coverage
consolidation are complete. These are not pending research tasks. The
software extensions below are not prerequisites for the stated theorem.

| Remaining question | Concrete boundary |
|---|---|
| Certified bounded-input front end | Optional software: admitted Hopf inputs to an actual coarse word and certified residual data. The bounded-input proof already supplies construction and output bounds with its explicit small-system exception |
| Variable-size native schedule | Optional software: general tables, predicates, reflections, and banked count/depth scheduling. The bounded component-to-gradient pass is complete |
| General guarded decoder | Optional software: replace the general floating-point contractions with the proved certified arithmetic. Exact fixture decoders cover only their fixed target |
| Constant-clean complete-frame endpoint | Still open independently of the state-based task. A new candidate must supply an explicit complete native identity and symbolic precision/workspace ledger before another fixture pass |
| Optimal T-depth | Uniform matching through sqrt(N/n) above the literal threshold; for 6 ≤ L ≤ log₂(n+2)/16 this extends through sqrt(NL/n). Larger-width optimality remains open |

Do not repeat the completed coefficient-to-row, two- and four-row lookup,
one- and two-system-qubit amplification, or coherent residual-selection passes.
The selected bounded residual-to-QBP integration is also complete. No
further fixture expansion or general software API is selected. Begin a
future pass at the coverage map: specify either a concrete Hopf-QBP
input/output requirement and its missing interface, or a new scientific
claim and its proof obligation. A general compiler is not required to
close this selected Hopf-QBP task.
The guarded-batch logarithm is now amortized, and the variable-accuracy
composition is complete; do not repeat either task.
The selected fixed-accuracy depth milestone at sufficiently large
$`b=\Theta(\sqrt N)`$ now has $`D_T=O(n)`$, retaining
optimal-order $`T=\Theta(\sqrt N)`$. Its available depth bounds
are $`\Omega(1)`$ and $`O(n)`$; depth optimality remains open.
The capped precision allocation has removed the accumulated source widths
as a quadratic contribution: sources and suffix predicates now cost
$`O(n\log(n+1))`$ depth at fixed L. The retained per-layer routing
sums to $`O(n^2)`$ in that older schedule. The new hybrid prices
complete queries, including the replacement indicator work.

Two attempted shortcuts do not yet supply such a circuit. Keeping an
address-controlled router open conjugates an intervening table matrix A
to $`AP_x`$ or $`P_xA`$; this is another address-dependent operation
that must be implemented. Running independent equality tests in parallel
requires read-only shared address controls and exactly returned private
dirty work. The audited Khattar–Gidney logarithmic-depth dirty MCX circuit temporarily
borrows its controls, so its row instances cannot simply overlap them.
Dirty CNOT copies also carry unknown masks. These are gaps in the proposed
schedules, not unrestricted impossibility results. The successful route
uses a shallow exact dirty-indicator batch with a full live-width ledger.
Polynomial overhead is affordable on early layers, whose table sizes
decay geometrically, while retaining the existing queries near the leaves.
The counter construction below now makes that hybrid unconditional.

The [bilinear query reduction](docs/PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction)
now removes the word-bank router from this proposed route. Split a table
address into two parts and allocate two arbitrary dirty indicator words.
A rank-reduced bilinear map, interleaved with two calls to each indicator,
extracts the selected table entry and returns all work exactly. With the
current routed indicators, an asymmetric split retains count-efficient
queries but still has linear address depth; it does not improve the
unconditional frame bound.

The polynomial-overhead indicator is supplied by a
[dirty-counter construction](docs/PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters)
with [masked compression](docs/DIRTY_SUM_COMPRESSION.md).
A reversible sum echo adds the Hamming weight into an arbitrary dirty
counter; a cyclic-rotation echo removes its unknown offset. A final
indicator echo extracts the all-literals-true predicate. Every helper
returns exactly. Original address bits enter only through Clifford
CNOTs; all non-Clifford gates act on private row work and parallelize.
For k address bits the construction has
$`T,G,w=O(2^k(k+1)^3)`$ and
$`D_T=O(\chi(k))`$ with no clean work.

The early/late hybrid now works at every eligible width and precision.
It gives depth $`O(NL/b^2+nL+n\chi(n))`$ without changing count,
error, or the $`17B_0`$ threshold. In particular, fixed accuracy gives
$`T=O(\sqrt N)`$, $`D_T=O(n\chi(n))`$, and $`G=O(N)`$
at sufficient $`b=\Theta(\sqrt N)`$. Polynomial overhead on early
layers is absorbed by their geometrically smaller tables; the last
$`O(\log n)`$ layers retain the existing count-efficient queries.
This is an unconditional depth improvement with optimal-order count,
not a matching depth theorem or a practical crossover estimate.

The [two-bit base case](docs/PARALLEL_DIRTY_LOOKUP.md#12-the-two-bit-indicator-has-optimal-t-depth-two)
also has optimal T-depth two, with eight T gates and no extra helper.
A Pauli-transfer argument excludes one T layer for every nonaffine
all-input classical permutation even with arbitrary returned dirty
helpers. This constant lower bound does not apply to initialized-clean
subspace contracts and does not close the frame-depth gap.

The signed-translation echo and conditional complement supply a linear
TTK implementation of a read-only controlled increment. A separate
involution refinement uses two dirty bits and Vandaele's controlled
increment to reduce its depth to $`O(\log(m+2))`$ for m counter bits.
All public-address interactions remain CNOTs.

The masked sum now uses overlapping carry pipelines. Keep each column's
raw bits fixed while expressing its parity forest as linear forms.
Increment each fresh upper word in doubling blocks, testing the already
updated lower prefix for zero. A completed bit is never written or
borrowed again during the pipeline. Exact dirty-helper phase gadgets
let later columns read those bits while higher blocks continue. A bit
at distance d settles in $`O(d)`$ T layers, so launching column j at
time proportional to j finishes every carry in $`O(m)`$ T-depth.
Only then apply the deferred parity CNOTs. Two linear TTK additions
read the final words, and the actual inverse returns all work.
The helper-only offsets still cancel in the outer echo. The geometric
width bound and balanced-forest Clifford count remain polynomial.
The pipeline itself needs no fast increment or fast-adder import.

This completes the selected $`O(\log(k+2))`$ indicator milestone and
gives $`O(n\log(n+2))`$ fixed-accuracy frame depth with optimal-order
T-count. The completed [error audit](docs/HOPF_ERROR_ACCUMULATION.md)
now separates ideal angular error from physical source leakage. Exact
angle perturbations have a sharp square-sum bound with asymptotic constant
$`\sqrt2`$. For equal error magnitudes within each layer, the entire
finite relative spectrum is independent of the base angles and signs.
Nearest rounding on angular grids with a common spacing therefore still
requires half a logarithm of n in that precision at fixed accuracy, even with
coordinated midpoint tie choices. This is a restricted angle-grid result.

For the actual common-precision source stages, local errors at most
$`4\Delta`$ can yield full-frame error at least $`n\Delta/8`$ when
$`\Delta\le(256n)^{-2}`$. Their full initialized-isometry errors have
no uniform square-sum bound. This counterfamily is not a lower bound for
the capped allocation, a replacement compiler, or fixed-accuracy depth.
The count/depth frontier and high-precision endpoint remain unchanged.

The [short-echo audit](docs/HOPF_FLAG_ECHO.md) is now complete. The
proposed half-angle reflection composite is exactly the source square,
and all four diagonal Pauli flag echoes retain linear radial error at
generic angles. The two-flag-Z echo has an exact full-space identity
when the programmed cosine and sine masks coincide, but this does not
give uniform cancellation.

The [radial filter](docs/HOPF_RADIAL_FILTER.md) supplies a positive
alternative: the standard pi-over-three fixed-point sequence, with
literal phase correction and six native phase approximants, makes the
full error relative to the encoded polar rotation quadratic in the
radial defect. It uses the same two clean flags and returned dirty pool.
The modified frame obeys angular square-sum plus quadratic remainder
control. It permits

```math
m_d=L+4+\min\{n-d,\lceil\tfrac12\log_2(8n)\rceil\},
\qquad \sum_dm_d=nL+\tfrac12n\log_2n+O(n).
```

This is a reduction in source widths, not a factor-two native gate
saving: the filter triples the amplified calls and its phase words
also cost $`O(m_d)`$. The existing same-circuit T/G/depth frontier and
matching interval remain valid, and both large scientific gaps remain
open. Native phase-word existence is proved; the bounded fixtures do
not emit a variable-size filter or assert efficient fine-word search.

The [conditional geometric source](docs/CONDITIONAL_GEOMETRIC_SOURCE.md)
now supplies a sublinear precision-depth replacement on eligible layers.
A source of width m reserves at most 7m zero suffix bits on the active
sector. A reversible prefix-OR preparation has O(m) native count and
O(log(m+2)) T-depth; its tree scratch is erased before the controlled-H
batch, and its prefix outputs are then invariant. Exact inactive
cancellation permits arbitrary suffix inputs outside the active sector.
The same two external flags hold the suffix predicate and branch. A
charged enlarged reflection and literal conditional amplification phase
complete the full-isometry contract, including final predicate uncomputation.
Use the original unfiltered additive precision cap for this schedule.
At fixed accuracy, its source/reflection depth totals
O(n log log(n+2)+log²(n+2)), including the late serial fallback.
Queries and suffix predicates remain O(n log(n+2)) in that source-only
schedule. The grouped composition below now improves the full fixed-L
depth frontier; the high-precision endpoint is unchanged.

The [two-layer obstruction](docs/SHALLOW_SOURCE_OBSTRUCTION.md) also
extends beyond Majorana-preserving Clifford stages. A full-space
two-T-layer transfer entry is zero or a signed inverse power of sqrt(2).
Geometric-source witnesses therefore force operator error at least 1/16
for m≥5, and for controlled sources at m≥4, at every returned dirty
width. This is a constant bound for an isolated full-input primitive;
initialized-clean source interfaces and frame circuits are excluded.

The [grouped program construction](docs/GROUPED_PROGRAM_PREFETCH.md)
is now complete. A height-g group reserves at most 16m2^g conditional
suffix bits, including 2m(2^g−1) program bits, source, selectors, and
private phase-mask work. Every program bit is read-only through actual
source leakage, so the final query inverse erases it exactly. The inner
enable depends only on the remaining local group bits. Its amplification
minus is CZ(h,u), not Z_h. Prefetch bits use private dirty indicator pools
in parallel, charging both the polynomial T/width cost and Qw Cliffords.

The old late fallback would have retained n log n depth. The new
[chunked indicator](docs/CHUNKED_DIRTY_INDICATOR.md) avoids that trap:
a dirty tree at chunk boundaries propagates a root flip by F X F†,
implemented chronologically by F†, X, F. A leaf-readout echo returns
every dirty node. Chunk length ell gives O(2^r(ell+2)^3) count/work
and O(ceil(r/ell) log(ell+2)) T-depth. For a layer k steps from the
leaves, use ell=min(address-half length,2^floor(k/12)). Its count is
summable within O(sqrt N), while all final O(log n) queries have O(n)
total depth. No unknown dirty program is treated as an initialized cache.

The earlier geometric grouped theorem uses the unfiltered additive precision:
m=L+4+ceil log2(8n), cutoff R=min(n,256m), and early height
g=floor log2(k/(64m)). These groups fit the literal suffix reservation,
have O(n/log n) prefetches, and cost O(n log log n) depth internally.
Late source/predicate costs are O(log²n). The
[incremental selector schedule](docs/GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
now removes the internal selector/enable logarithm: consume suffix enables
before their controls change, grow prefix nodes after each target is final,
and erase the retained prefix tree in reverse. It has O(g) group depth,
fits the existing 16m2^g suffix reservation, and returns its work exactly
through source leakage and on arbitrary inactive inputs. Across early
groups this contribution is O(n). Source/reflection depth still retains
the log log n factor in that geometric construction. The unary phase
source below gives the current complete-frame frontier.

The [common-source audit](docs/GROUPED_PROGRAM_PREFETCH.md#11-a-common-source-identity-and-the-remaining-reflection)
gives an exact two-layer conjugation identity, but two conjugated success
reflections remain per layer. A legal exact zero-angle row leaves a stale
success monitor with rejected norm sqrt(3)/2. Replacing individual Q calls
by one-use fresh banks also leaves constant rejection unless a later
operation mixes that bank's success and failure sectors. These are scoped
failures of specified substitutions, not general depth lower bounds.

The [unary phase source](docs/UNARY_PHASE_GRADIENT.md) supplies a different
route. A Karatsuba rank decomposition implements a coherent one-hot cyclic
shift with private conditional work and literal full inactive identity.
The ideal Fourier source supplies signed target phases; its charged native
preparation and actual inverse cost at most twice the preparation error
for the entire group. No geometric success reflection is retained.
The natural suffix cutoff reserves the convolution pool. Early groups
have O(log n) depth and number O(n/log n). The entire remaining tail uses
chunked queries at the original capped source precision, giving O(n)
total depth while preserving O(sqrt N) T-count and O(N) Clifford count.
The six bounded tests cover the rank identity, native guard, arbitrary
source shift, preparation/inversion, unequal rows, and full source-return
error. They do not constitute a scalable native group compiler.

The [blocked bilinear query](docs/BLOCKED_BILINEAR_LOOKUP.md) now proves
the width-dependent extension: depth O(N/b²+n) and T-count
O(sqrt N+N/b), with two flags and the original 17(L+n+7) threshold.
A naive use of the earlier width-sensitive query fails because its
per-layer route has depth O(n). Instead, one returned dirty helper
implements each selected bilinear leaf in constant depth per output
bit. Two dirty traversals select its matrix block; the four-corner
indicator echo returns all masks. Chunk lengths and block sizes are
chosen together, and every tail resource sum is charged. This extends
the simultaneous fixed-accuracy matching interval to sqrt(N/n).

The [uniform precision pass](docs/UNIFORM_PRECISION_DEPTH.md) now exposes
the eta-dependent source width, conditional suffix cutoff, and weighted
query sums. Its rectangular allocation avoids the balanced query's
extra square root of precision in T-count. It proves the uniform bounds
and the stronger slowly growing precision corollary above. Bounded
integer checks cover allocation floors and workspace constraints; they
do not replace the asymptotic proof or implement a scalable compiler.

A next depth pass should target the remaining large-width gap: identify
a concrete way to share work between unary groups or a full-frame
lower-bound invariant that permits arbitrary Clifford interlayers.
The linear upper bound is not itself a lower bound. The high-precision
endpoint remains a separate question. Do not repeat the completed
selector, source-reuse, unary source, blocked-query, or precision-splice
proofs. Retain literal phases, actual inverses, full work return, and
charged preparation in any new schedule. The completed modest-width
matching theorem does not require solving the high-precision endpoint.
For any new component, keep literal phases and actual inverses, declare
all initialized inputs, include borrowed-work return in its isometry
error, and charge each reflection before composing a larger state compiler.

For the selected complete real-frame endpoint,

```math
N=2^n,\qquad a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

```math
\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n)),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

No new construction is selected for this gap, and recent task-specific
results do not show that its resolution is close. The established frame
results and the changed-decoder QBP theorem do not require its closure.
Do not replace the current task by arbitrary-unitary synthesis or invent
an uncharged observable, initialized history, or coherent evaluator.

## Completed routes that should not be repeated

The [canonical grouped audit](docs/CONDITIONAL_SUFFIX_COMPILER.md) retains
literal inverse branches, occupied flags, and full-work return, but still
pays precision per group. Its paired-source consolidation also improves
the old scalar baseline, so it is closed as the proposed amortization
mechanism. The [source-carry analysis](docs/SOURCE_REUSE_LIMITS.md) solves
width transport without pricing a cheaper joint interior program.

The [Hopf scattering step](docs/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
has one precision charge, but its unchanged-coin feedback conversion needs
a growing number of queries at fine accuracy. This is a query restriction,
not an additive native T lower bound. Globally programmed coefficients and
fully charged basis changes remain eligible for a different construction.
The [selection audit](docs/OPEN_PROBLEM.md#selection-audit-after-canonical-completion)
records the remaining complete-frame requirements.

The earlier [leaf-reference decoder](docs/REFERENCE_STATE_QBP.md) remains
an option under its explicit sampling tradeoff. The current coarse-frame
decoder removes that angle-dependent factor. Its general sampling bound,
complex gauge, and bounded-input preprocessing are settled; repeating a
special exact fixture would not strengthen those analytic claims.

## Restore and verify

Run from the repository root with `requirements.txt` installed:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python scripts/complex_coarse_native_example.py
python scripts/residual_qbp_native_example.py --q 16
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

Merge gates include Python 3.11/3.13 validation, all four exact-receipt
suites, and rendered presentation. The [verification map](docs/VERIFICATION.md)
distinguishes analytic proofs, finite native examples, exact rational
certificates, and floating-point decoder checks. Internal review and
finite tests are not external peer review or asymptotic proof.

No human action is required to resume repository work. Preserve the theorem
contracts and record any genuinely new question in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md), with its proof in the appropriate
chapter, so continuation does not depend on an old chat.
