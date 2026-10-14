# Documentation map

[Landing page](../README.md) · [Read the argument](../REVIEW.md) · [Publication scope](../manuscript/PUBLICATION_SCOPE.md)

The selected compiler package has four principal results: exact compilation,
sufficient-clean matching T-count, the optimal one-clean banked count
frontier, and simultaneous T-count/T-depth tradeoffs for complete real frames.
Start with the landing page and `REVIEW.md`, then follow the
relevant claim below.

## Proof chapters

### Shared target and error contract

State preparation fixes one column. The compiler target fixes the state and
prescribed differential-frame columns, with the declared workspace returned.
Approximate compilation includes all logical inputs, arbitrary dirty inputs
and their references, literal phases, and clean-work leakage.

| Read | Purpose |
|---|---|
| [Hopf interface](HOPF_INTERFACE.md) | Coordinates, marker columns, singular charts, addressed layers, and full-frame tangent metric and operator stability |
| [Frame-safe compilation](FRAME_SAFE_COMPILATION.md) | Complete-input error, actual adjoints, fixed-decoder necessity, and singular exceptions |
| [Compiler boundaries](COMPILER_BOUNDARIES.md) | Concrete distinctions between state preparation, checkpoints, and prescribed complete frames |

### A. Exact compilation at every clean-workspace budget

The [exact compiler theorem](COMPILER_THEOREM.md) proves matching elementary
size and depth, and a separate CNOT-count bound, for every number of clean
work qubits. It covers real and phase-dressed complex magnitude frames in
the arbitrary-one-qubit/CNOT model.

Read its echo, decoder, router, and maximal-cut constructions together with
the lower bounds. The proof returns workspace exactly and preserves every
required column; a state-preparation theorem alone does not establish that
contract. The exact gate model is distinct from the Clifford+T model below.

### B. Matching T-count with sufficient clean workspace

The [fault-tolerant compiler theorem](FAULT_TOLERANT_COMPILER.md) gives the
matching worst-case precision/workspace frontier under its sufficient clean
reservation. Real and phase-dressed complex magnitude families have their
stated extensions. Follow the complete source, residual dictionaries,
compressed failure history, amplification, and resource ledger.

Section 10 gives the real-frame reduction to diagonal synthesis, the
fixed-width counting argument, and a literal-phase root witness giving
$`T\ge2N-1`$ at $`\eta=2^{-N}`$, valid at arbitrary width.
The same probability gap gives a threshold below which approximate workspace
return is exact. These lower bounds also support C and D;
they are not obtained by adding costs of particular source calls.

### C. One-clean compilation and its corollaries

The [growing-source count theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives $`T=O(N+L)`$ and $`G=O(NL)`$ for both frame families,
every $`n\ge1`$, and integer $`L\ge6`$, with error $`2^{-L}`$
and fixed $`a=1`$.
The base dirty reservation is $`B=L+n+7`$ for real frames and
$`B=L+n+8`$ for complex magnitude frames. At $`b\ge2B`$,
the [banked theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
proves the matching frontier

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right).
```

At exactly $`b=B`$, the unbanked count is matching when
$`L\ge\max\{6,n\}`$. The proof chain is:

1. [One-clean source and layerwise baseline](ONE_CLEAN_COMPILER.md#6-complete-real-frames-and-the-grouped-extension): exact conditioning, full-output amplification, the low-precision banked bound, and the bounded-dimension fallback.
2. [Growing sources](../research/endpoint/GROWING_SOURCE_COMPILER.md): mixed local gauges, sealed modes, actual dirty replay, compensated overlap, and the high-precision regime.
3. [Conjugate recovery](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md#2-lift-the-whole-high-level-word): borrow only the phase signal, retain the initialized anchor, and charge every native interlayer.
4. [Banked block/source interfaces](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces), [literal diagonals](ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors), and the [existing lower bounds](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage): sum native costs, compose phases, and match at the stated widths.

The low-precision baseline and high-precision source meet at
$`L=N/n^2`$; their constants are uniform. Conditional logical work
uses the [exact inactive guards](CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally).
The grouped theorem is a retained corollary, not a premise of this splice.

| Supporting proof or corollary | Why it remains in the package |
|---|---|
| [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) | Full-input two-clean source, packed diagonal precision/workspace bound, amplification/error certificate, layerwise frame baseline, and separate multiplexor reservation |
| [Conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) | Exact conditional-work guards used by C; residual columns and the separately retained grouped count theorem |
| [One-clean compiler](ONE_CLEAN_COMPILER.md) | Native one-clean primitive; literal diagonals, complete one-target U(2) multiplexors, phase-dressed magnitude frames, and the separate zero-clean real-frame corollary |
| [Borrowed-workspace compiler](BORROWED_WORKSPACE_COMPILER.md) | Arbitrary clean/dirty allocations, exact dirty traversal, and the restricted matching splice |

Keep each corollary's literal dirty reservation. Reducing the clean count
does not supersede every earlier resource point. The complex magnitude
extension is a real frame followed by supplied leaf phases, not arbitrary
complex-unitary synthesis.

### D. Simultaneous real-frame T-count and T-depth

Start with the [all-precision theorem](../research/depth/ALL_PRECISION_DEPTH.md).
For every $`n\ge1`$, integer $`L\ge6`$, two clean flags, and
$`b\ge17(L+n+7)`$, it gives one circuit with full error at most
$`2^{-L}`$ and absolute constants. Define

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
```

Then

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

Count is optimal at every eligible width. Count and depth match through
$`b\le\sqrt{NL/F(n,L)}`$, subject to the same lower threshold.
For $`L>N/n^2`$, the stronger depth bound is $`O(NL/b^2+L)`$.
The [compressed-column special case](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains $`\log_2(n+2)`$ in place of $`\log_2(n+L+2)`$ for
$`6\le L\le cn`$, with some absolute $`c>0`$.
At sufficient $`b=\Theta(\sqrt{NL})`$, shrinking error retains
sublinear depth for $`L=o(n/\log(n+2))`$; fixed accuracy retains
$`O_\eta(\sqrt n\log(n+2))`$ depth. Optimal depth remains open.

| Selected branch | Canonical proof and role |
|---|---|
| All-precision composition | [All-precision depth](../research/depth/ALL_PRECISION_DEPTH.md), R70: simultaneous counts, literal widths, precision splice, and inherited matching bounds |
| Whole-frame binary schedule | [Whole-frame binary depth](../research/depth/WHOLE_FRAME_BINARY_DEPTH.md), R68: charged complete-frame schedule with precision-dependent source work |
| High-precision branch | [High-precision depth](../research/depth/HIGH_PRECISION_DEPTH.md), R69: $`O(NL/b^2+L)`$ depth for $`L>N/n^2`$ with optimal count |
| Retained baseline | [Uniform precision and depth](UNIFORM_PRECISION_DEPTH.md): the $`nL`$ depth alternative, literal-width fallback, and earlier matching intervals |

The retained uniform proof has an early unary-group branch, an exact rectangular-query tail,
and an all-precision fallback. The following chapters are mathematical
dependencies, not a sequence of research tasks to repeat.

| Part of the argument | Proof dependency and role |
|---|---|
| Uniform composition | [Uniform precision and depth](UNIFORM_PRECISION_DEPTH.md): rectangular allocation, uniform cutoff, weighted sums, literal-width fallback, and matching inequalities |
| Selected table queries | [Blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md): controlled bilinear helper echo, selected-block dirty traversal, four-call cancellation, and fixed-accuracy composition |
| Charged early source | [Unary phase-source groups](UNARY_PHASE_GRADIENT.md): exact selected cyclic shifts, conditional work, source preparation/actual inverse, and complete group error |
| Program and selector lifetime | [Grouped program prefetch](GROUPED_PROGRAM_PREFETCH.md), especially Section 10: preserved programs, incremental prefix selectors, and consume-before-change suffix enables |
| Geometric predecessor | [Conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md): the complete conditional-source interface used by the retained geometric grouped construction |
| Tail indicators | [Chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md): full-input dirty-tree echo, live-work bound, and the summable chunk allocation |
| Shallow dirty arithmetic | [Masked dirty-sum compression](DIRTY_SUM_COMPRESSION.md): compressor invariant, weighted-count recurrence, stable-prefix pipeline, and exact offset cancellation |
| Literal-width fallback | [Amortized dirty lookup](AMORTIZED_DIRTY_LOOKUP.md): exact selected-shear query, capped source precision, and the same-circuit resource/error ledger |
| General-precision hybrid | [Parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md): bilinear queries, dirty-counter primitives, and the hybrid used outside the low-precision branch |
| Native scheduling and depth bounds | [T-depth schedule](T_DEPTH_COMPILER.md): exact shared-control phase/Fredkin batches, routing, two-dirty predicates, and count-to-depth lower-bound transfer |
| Angular approximation | [Hopf error accumulation](HOPF_ERROR_ACCUMULATION.md), Section 1: full-frame square-sum bound for ideal angle perturbations |

The stronger linear-precision special case retains its separate dependency chain:

| Part of the sublinear-depth argument | Proof dependency and role |
|---|---|
| Complete early operator | [Sublinear T-depth](../research/depth/SUBLINEAR_T_DEPTH.md): orthogonal column sectors, dimension-independent source completion, literal-phase recovery, and early/tail resource closure |
| Coherent descendant source | [Descendant loader](../research/depth/COHERENT_DESCENDANT_LOADER.md): certified Boolean-state selection, preserved tag/prefix sectors, exact guarded arithmetic, and full source-isometry error |
| Early shared precision source | [Shared-gradient recovery](../research/depth/SHARED_GRADIENT_RECOVERY.md): one native source preparation/actual return, priced arithmetic and coefficient tables, invariant gradient sector, and recovery reflection excluding the live core |
| Native early data access | [Shallow dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md): private T-layer scheduling with shared read-only addresses, Boolean phase echoes, and guarded variable-suffix Hadamards |
| Shared binary terminal stream | [Shared binary terminal compilation](../research/depth/SHARED_BINARY_TERMINAL.md): a separate capped geometric source, retained failure history, one complete-stream recovery, a protected logarithmic bank, and capped bottom layers |
| Coarse-program interface | [Whole-word prefetch](../research/depth/GROUPED_TERMINAL_DEPTH.md#2-prefetch-the-exact-coarse-native-word): exactly returned native programs; the chapter's grouped terminal theorem remains a standalone alternative |
| Native terminal queries | [Nonuniform indicators](../research/depth/NONUNIFORM_DIRTY_INDICATOR.md) and [capped source precision](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision): linear indicator count/work, logarithmic address depth, and retained complete-input layer error |

The [entry-indexed encoding bound](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
separately identifies the sharp worst-case square-root normalization cost
of that source architecture, including operator-norm approximation. It is
not an unrestricted circuit-depth lower bound or an input to the upper bound.

Older constructions remain where they supply a primitive, error certificate,
or fallback with a distinct reservation. The masked compressor and its count
recurrence are needed by the later pipeline; the grouped selector schedule
is needed by the unary source. Their proofs must not be replaced by a link
to a stronger headline bound.

D is a real complete-frame theorem. The complex extensions of A–C and the
state-based protocol have separate contracts. T-depth allows arbitrary
Clifford interlayers; their gate counts are charged, and their physical
depth is not declared zero. The large-width depth optimum remains open.

### Fixed-decoder QBP consequence

The [QBP consequence](QBP_CONSEQUENCE.md) explains why preserving the
additional columns matters for the original inverse-frame decoder.
[Approximate QBP](QBP_APPROXIMATION.md) proves fixed-parameter robustness,
complete complex gradients, observable-sum and classical-weight errors,
and concentration with correlated dirty-bank reuse.

These statements separately charge quantum executions, observable access,
initialized protocol work, and classical output. They do not establish
optimality among all gradient algorithms or an end-to-end advantage.

<a id="separate-supplements-and-optional-research"></a>
<a id="separate-endpoint-theorem"></a>

<a id="selected-endpoint-count-theorem-and-complex-magnitude-corollary"></a>
<a id="selected-precision-uniform-two-clean-count-theorem"></a>
<a id="selected-growing-source-count-theorems"></a>

### Retained two-clean count corollaries

Result C above is the authoritative one-clean count statement. The
[two-clean precision-uniform theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
supplies its precision/mode separation. The [fixed complex endpoint](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
retains a distinct reservation: $`a=2,b=N+n+7,\eta=2^{-N},n\ge3`$,
with worst-case $`T^\star=\Theta(N)`$ and lower bound $`2N-1`$.
It saves one dirty qubit relative to C's one-clean complex reservation.
The [research index](../research/README.md) maps supporting constructions
and method-specific limits; the [endpoint statement](OPEN_PROBLEM.md)
records the remaining depth-optimality question.

### Completed changed-decoder supplement

The [state-based QBP supplement](../supplements/state_based_qbp/README.md)
uses a changed decoder and an initialized state-preparation interface. Its
real/complex theorem, computational input audit, task costs, and bounded
native integrations are grouped there. It does not replace the complete-frame
contract or become an extra prerequisite for A–D.

## Evidence and sources

| Read | Purpose |
|---|---|
| [Core-claim audit](CORE_CLAIM_AUDIT.md) | Recorded internal claim-to-proof review and its evidence limits |
| [Verification](VERIFICATION.md) | Analytic proofs, exact certificates, native finite checks, and missing general-emitter coverage |
| [Fault-tolerant receipts](../verification/fault_tolerant/README.md) | Four standalone exact source/kernel/resource checks |
| [Source map](SOURCE_MAP.md) | Imported premises, inherited interfaces, and local constructions |
| [Related work](RELATED_WORK.md) | Primary-source comparisons by target, gate model, precision, and workspace |
| [Provenance](../provenance/README.md) | Source versions and upstream lineage |

Finite checks test their specified contracts; they do not prove asymptotic
theorems, establish priority, or supply independent peer review. No general
elementary Clifford+T emitter for the asymptotic grouped compiler is claimed.

The [implementation](../compiler_robust_hopf/README.md),
[tests](../tests/README.md), [scripts](../scripts/README.md), and
[diagrams](../assets/README.md) provide reproduction maps.
The [manuscript guide](../manuscript/README.md) gives the selected argument's
section order.
