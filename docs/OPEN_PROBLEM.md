# Research status and the constant-clean endpoint

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [One-clean compiler](ONE_CLEAN_COMPILER.md) · [Grouped refinement](CONDITIONAL_SUFFIX_COMPILER.md)

This page gathers the retained results, the limits of explored routes, and
the next construction to test. It is a research checkpoint, not an additional
compiler theorem. The [publication scope](../manuscript/PUBLICATION_SCOPE.md)
gives the selected publication results and the [verification map](VERIFICATION.md)
separates analytic proofs from finite checks.

## Established frontier

Write $`N=2^n`$, $`q=n+a+b`$,
$`h=1+\lceil\log_2(L+n+2)\rceil`$, and
$`\ell_*(n)=1+\log_2^*(n+2)`$. The exact model has arbitrary one-qubit
gates and CNOTs; the approximate model has coherent Clifford+T gates.
The clean budget in the exact model is m; in the approximate model a and b
count clean and arbitrary dirty qubits. The precision parameter L is defined
below. Lower bounds are worst-case over the stated frame family.

| Question | Retained result | Status and proof |
|---|---|---|
| Exact size and depth versus clean workspace | Size $`\Theta(N)`$ and depth $`\Theta(n+N/(n+m))`$ for every $`m\ge0`$; CNOT count $`\Theta(N)`$ for $`n\ge2`$, zero for $`n=1`$ | Matching for the complete real and phase-dressed complex magnitude frames; [exact theorem](COMPILER_THEOREM.md) |
| T-count with sufficient clean workspace | $`T^\star=\Theta(\sqrt{NL}+L+NL/q)`$, $`G=O(NL)`$, when $`a\ge C(n+h)`$ for sufficiently large fixed C | Matching for those same families; the clean reservation is sufficient, not proved necessary; [fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md) |
| Real-frame T-count at arbitrary workspace budgets | $`T=O(NL/q+L\sqrt N)`$, $`G=O(NL)`$, for every $`a,b\ge0`$ | Splicing with the sufficient-clean theorem gives the matching frontier above for every a when $`h+b\le c\sqrt N`$, for fixed $`c>0`$; [borrowed-workspace proof](BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements) |
| Zero-clean layerwise real-frame T-count | $`T=O(N+nL)`$, $`G=O(NL)`$, at $`a=0`$, $`b\ge L+n+7`$ | Full-operator approximation using a borrowed amplification signal; [signal-symmetry corollary](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) |
| Grouped one-clean real-frame T-count | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge L+n+7`$ | Precision-uniform grouped bound without additional word banks; [one-clean extension](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean T-count with additional dirty banks | $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge2(L+n+7)`$ | Matches the lower bound if $`L\ell_*(n)^2\le N`$ or $`b\le N/\ell_*(n)`$; these are sufficient regimes; [banked one-clean proof](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean phase-dressed complex magnitude frame | The same grouped and banked T-counts, at $`b\ge L+n+8`$ and $`b\ge2(L+n+8)`$, respectively | Compose the real compiler and literal phase diagonal in the same workspace; [composition corollary](ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames) |
| T-depth with additional dirty banks | $`D_T=O(NL/b+\min\{nL+n^3,L\ell_*(n)+n^4\})`$ at $`a=2`$, $`b\ge2(L+n+7)`$, with $`T,G=O(NL)`$ | Choose between the layerwise and grouped [schedules](T_DEPTH_COMPILER.md); real frames; optimizing depth may increase T-count; no matching frontier established |
| Simultaneous T-count and T-depth | $`T=O(\sqrt{NL}+L\ell_*(n))`$, $`D_T=O(\min\{nL+n^3,L\ell_*(n)+n^4\})`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge C(L+n+7+\sqrt{NL})`$ | Same real-frame circuit, for sufficiently large fixed C; [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md); T-depth optimality remains open |

Take the best applicable construction. For fixed L, $`a=2`$ and
$`b=L+n+7=\Theta(n)`$, the arbitrary-budget matching splice gives
$`T^\star=\Theta(N/n)`$, sharper than the grouped estimate. Extra clean
qubits may be left unused; the depth schedules retain their separately
proved two-clean allocation.

All frame constructions preserve the prescribed completion and the
[fixed-parameter QBP error contract](QBP_APPROXIMATION.md). They do not
differentiate discrete synthesis. Complex leaf-phase derivatives remain a
separate stream. Literal diagonals and one-target U(2) multiplexors retain
their independent matched banked frontiers; earlier weaker frame bounds
are dependencies or fallbacks, not the current frontier.

## The count and depth gaps are different

The exact ancilla-depth theorem is matching in its elementary-gate model;
T-count is matching under its stated clean reservation. The two-clean
T-depth schedule has no matching lower bound. Complete-frame safety and
QBP substitution are already established; the gaps concern resources.

For the exactly two-clean banked regime, put $`B_0=L+n+7`$ and assume
$`b\ge2B_0`$. Then $`q=n+2+b=\Theta(b)`$. The
[inherited depth lower bound](T_DEPTH_COMPILER.md#4-lower-bounds-and-the-remaining-depth-gap)
simplifies to

```math
D_T^\star=\Omega(1+NL/b^2).
```

This follows from $`L/b=O(1)`$ and
$`\sqrt{NL}/b\le1+NL/b^2`$, not a new lower-bound argument. The table
uses $`a=2`$ throughout and respects each sufficient allocation threshold.

| Regime | Depth lower bound | Available depth upper bound | Remaining issue |
|---|---|---|---|
| Fixed L, $`b=\Theta(n)`$ | $`\Omega(N/n^2)`$ | $`O(N/n)`$ | Factor-n gap |
| Fixed L, sufficiently large $`b=\Theta(\sqrt N)`$ | $`\Omega(1)`$ | $`O(n^3)`$ with $`T=O(\sqrt N)`$ | Simultaneous count/depth capability; depth lower bound remains unmatched |
| Fixed L, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(n^3)`$ with $`T=O(\sqrt N)`$ | Extra width is not needed by this schedule; depth optimality remains open |
| $`L=N`$, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ | Serial precision cost remains |
| Selected endpoint $`L=N,b=B_0`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ from $`D_T\le T`$ | The larger-bank depth theorem does not apply |

The square-root allocation gives optimal-order count and polynomial-logarithmic
T-depth in one circuit, but does not settle its depth exponent or smaller widths.

T-depth permits Clifford circuits of nonzero depth between its T layers.
It is not total circuit depth or elapsed QBP execution time. Neither the
exact CNOT light-cone bound nor the source's linear exact T-count minimum
supplies an additional depth lower bound in this model.

### What this already gives Hopf QBP

The endpoint below is a high-precision compiler question, separate from the
precision needed for a fixed raw-gradient accuracy. For a reflection-sum
observable with coefficient one-norm Lambda, the existing
[QBP error allocation](QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights)
permits

```math
L=\max\{6,\lceil\log_2(32\Lambda/\varepsilon_\infty)\rceil\}.
```

Thus fixed observable scale and coordinate accuracy give fixed L. In this
regime, the arbitrary-budget real-frame compiler already attains
worst-case optimal-order $`T=\Theta(\sqrt N)`$ with zero compiler clean
qubits and $`b=\Theta(\sqrt N)`$. The separate parallel schedule gives
this count together with $`D_T=O(n^3)`$ using two compiler clean qubits
and a sufficiently large square-root dirty allocation. Both preserve the
prescribed frame used by QBP; their clean allocations cannot be interchanged.

The protocol reserves one additional clean interference qubit and the
controlled observable's work, beyond the initialized system and compiler
reservation. Its fixed-accuracy, fixed-confidence execution count is
$`O(1+\log n)`$; observable costs and classical gradient output remain
separately charged. The unresolved $`L=N`$ endpoint does not prevent
these existing QBP guarantees. Optimal T-depth and end-to-end gradient
optimality remain open.

## The remaining endpoint

The publication's established compiler results do not depend on resolving
this question. Let $`T^\star_{F,\mathbb R}`$ be the worst-case minimum
T-count for the prescribed complete real Hopf frame. At

```math
a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

the retained frontier remains

```math
\Omega(N)\le T^\star_{F,\mathbb R}
\le O(N\ell_*(n))=O(N\log_2^*N).
```

Here $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ for
$`0\lt\eta\le1/64`$. Literal phases, clean-work leakage, dirty-work
return, and arbitrary reference correlations are included in the
[complete-input error contract](FAULT_TOLERANT_COMPILER.md#1-target-resources-and-theorem).
The same upper bound holds with only one clean qubit and the same dirty
allocation. Zero-clean layerwise synthesis gives $`O(N\log N)`$ there;
a sufficiently large $`\Theta(n)`$ clean reservation instead attains
$`\Theta(N)`$. These are sufficient constructions, not necessary clean
reservations.

The question is whether the displayed two-clean allocation admits an
$`O(N)`$ T-count construction, or whether a stronger general lower bound
holds. The optimum minimizes T-count: a valid linear-T construction closes
this gap even with a larger fully charged Clifford count. Preserving
$`G=O(NL)`$ is the stronger joint goal of the retained constructions.
The current upper bound does not cover every prefactor in
$`b=\Theta(N)`$. Source-call minima and restrictions on particular
intermediate interfaces do not strengthen the unrestricted lower bound.
Closing this count endpoint would still leave the T-depth question open.

## Current assessment: what the results establish

The antichain and sparse-update results prove linear bounds for promised
families. The subsequent passes establish complete boundary actions, an
explicit repair, cheap source-width transitions, compact Cayley data, and
native four-, eight-, and sixteen-mode benchmarks. The latest joint word
has fewer declared source calls, but the comparison words have the same
leading precision cost after fixed-mask simplification.
**None has narrowed the generic upper/lower gap.** Wide independent changes
and deep sparse nesting each admit one precision charge; the remaining
question concerns precision reuse for unrestricted comparable updates.

| Level | Established scope |
|---|---|
| General theorems | Exact matching resources, sufficient-clean matching T-count, and one-clean grouped bound; constant-clean endpoint and optimal T-depth remain open |
| Incomparable promised families | Antichain and sparse ancestor-closed updates have $`T=O(N+L)`$, $`G=O(NL)`$; generic rounding satisfies neither promise |
| Reusable blocks | Cheap coarse frame, linear residual data, complete weighted dilations, two-flag assembly, and a coupled repair with exact child-call cancellation; the surviving target precision remains charged |
| Completed small-example pass | Fixed-address product compression, four-mode native factors, and a full-port changing-target word with borrowed-signal return; no precision recurrence independent of growing support or group count |
| Scoped diagnostics | Specified source, query, repacking, truncation, and shared-conjugator failures; no additive full-frame lower bound follows |

### A promised antichain class has a linear endpoint bound

The [antichain theorem](ANTICHAIN_COMPILER.md) requires a supplied native
baseline C with determinant-one local words of length $`O(n-d+1)`$
at depth d. The target agrees **literally** outside a prefix-free changed
set, possibly of size $`N/2`$. No closeness is needed. With zero clean qubits,

```math
b\ge L+n+7,\qquad T=O(N+L),\qquad G=O(NL),
\qquad
\|\widetilde W-W\otimes I_b\|\le\eta.
```

The exact factorization $`W=VMV^\dagger C`$ removes the native
descendant forest; charged dirty swaps pack the corrections into one SU(2)
multiplexor. The full-operator error includes core and borrowed-signal
return. Forest and packing helpers return exactly.

### Sparse nested updates also have a linear endpoint bound

The [sparse-update compiler](SPARSE_UPDATE_COMPILER.md) uses the same
literal native-baseline promise. Write $`\mathcal S`$ for the ancestor
closure of the changed nodes, including ancestors whose local words did
not change. This support set is distinct from the affine residual operator
S used below. Put

```math
m=|\mathcal S|+1,\qquad s=\lceil\log_2m\rceil.
```

Under the sufficient condition $`n\ge2s+32`$,

```math
a=1,\qquad b\ge L+n+7,\qquad
T=O(N+L),\qquad G=O(NL),
```

```math
\|\widetilde WJ_1-J_1(W\otimes I_b)\|\le\eta.
```

An exact off-support forest and charged basis packing leave m active
modes on s logical bits. Conditional logical zeros support one dense
residual dictionary, scalar source, and half-block amplification. Only the
external predicate is initialized unconditionally. Its isometry error
includes leakage, work return, and dirty references; it is not a
full-operator guarantee for arbitrary input on that external qubit.

A root-to-leaf path has $`|\mathcal S|\le n`$ and satisfies the
bound for all n: use this construction at $`n\ge44`$ and the proof's
finite fallback otherwise. The bound is asymptotic, not a practical estimate.

The two promises are incomparable. A wide antichain can have a large
ancestor closure; a sparse closure can contain arbitrarily long nested
paths and some branching. For both theorems, a concrete real Hopf family
fixes unmarked angles at $`\pi/4`$ and allows arbitrary real changes
at the promised nodes. Its complete-frame approximation preserves the
existing QBP error interface. A generic native baseline may have complex
local words: literal agreement with it does not automatically define a
real Hopf frame or the separate phase-dressed complex magnitude family.

Generic target rounding may change every internal node. Then
$`|\mathcal S|=N-1`$ and $`s=n`$, leaving no packed zero sector;
the sparse dictionary is also too large for its current accounting.
Parameter count and support dimension diagnose why this proof does not
apply. They are not hardness results. The unresolved structure is
widespread changes across comparable tree nodes, with no proved reduction
to a constant number of the two solved families.
Simply applying their theorems in sequence is not such a reduction:
each requires literal agreement with a short native baseline outside its
own marked set, a promise not automatically preserved after a nonnative
update.

## A sufficient construction to seek

A sufficient endpoint construction is an actual unitary Q using two
initialized flags and at most $`N+n+7`$ dirty qubits, satisfying

```math
\left\|2J_2^\dagger QJ_2-(W\otimes I_b)\right\|\le\eta/4,
\qquad T(Q)=O(N),\qquad G(Q)=O(N^2).
```

These displayed bounds are a sufficient joint target, not the definition
of the T-only optimum. Every source, query, control, inverse, and helper
belongs to the circuit accounting. The accepted action must hold on all logical and dirty
inputs, including reference correlations. There is no assumed return
condition on rejected branches. The
[normalization-two amplification lemma](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
then gives the complete-isometry compiler with three charged calls and the
same two flags.

The [two-flag assembly](RESIDUAL_ASSEMBLY.md) realizes this interface at
$`O(N+nL)`$ T cost. A faster whole-residual block followed by the
charged coarse frame also suffices. The uniform target $`T=O(N+L)`$,
$`G=O(NL)`$ is stronger than settling the endpoint alone.

## Reusable ingredients and the resource bottleneck

These are retained constructions; their linked proofs are the primary homes.

| Ingredient | Established capability and remaining cost |
|---|---|
| Native coarse frame C | Fixed-accuracy $`T=O(\sqrt N)`$, $`G=O(N)`$ in the endpoint pool, no initialized work, exact helper return; any new control or mask must be charged |
| Compact residual | $`O(N)`$ classical local generators and path products for $`C^\dagger W'-I`$; coherent evaluation and target transport are not free |
| Weighted forward block | Complete two-sector dilation, fixed normalization, $`T=O(N+nL)`$, $`G=O(NL)`$ at $`b\ge L+n+7`$; the native synthesis signal is borrowed, distinct from its initialized logical dilation signal |
| Exact-return alternative | Same component with $`T=O(L\sqrt N)`$, $`G=O(NL)`$ at $`b\ge n+\lceil\sqrt N\rceil+7`$; different width/return point, not the best endpoint count |
| Affine/reverse assembly | [Two-flag proof](RESIDUAL_ASSEMBLY.md): diagonal-plus-forward and inverse-based reverse blocks give normalization two, with controls/inverses charged; cost still contains $`nL`$ |
| Full-port hierarchy | [Fusion audit](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths): recursively closed complete unitary, linear classical generators; faster native synthesis unproved |
| Coupled completion | [Whole-residual boundary](RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost): one-signal normalization-two recursion; a commutator implements the complete rank-at-most-four repair without normalizing transported differences; its child calls cancel to the original target wrappers, whose precision remains charged |
| Grouped frame | Best general endpoint bound, with conditional suffix and core return in the complete error |
| Product-first residual coordinates | [Cayley recursion](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) gives constant-size local data for the complex-coarse residual; the [four-mode native benchmark](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) uses two one-target programs, while the eight-mode controlled coupling and general coherent conversion remain charged |
| Source-width transport | [Reverse-order loader](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation): total boundary T-count $`2(m_{\max}-1)`$ with linear native loaders; transformed group bodies remain charged |
| Joint source body | [Changing-target word](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) shares one fixed scalar conjugator and returns a borrowed signal through two parity boundaries; the eight-/sixteen-mode words are priced, but fair fixed-mask simplification leaves the same leading precision cost in both comparisons |

The [weighted norm proof](ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
avoids a height penalty via uniform coarse subtree accuracy. The
[forward block](WEIGHTED_TRANSPORT_BLOCK.md) charges marker gathering,
queries, and certified preprocessing at zero defects. Real-rotation X
symmetry permits a borrowed native signal; the scalar-phase word lacks
that established symmetry.

For the common algebraic target approximation W', the exact residual is

```math
C^\dagger W'=S+R,\qquad S=A+F,\qquad
\alpha=4\varepsilon_0,\qquad \sigma=2-\alpha.
```

For fixed $`0\lt\varepsilon_0\le1/64`$, the existing normalized
branches satisfy

```math
\|S/\sigma\|\le
\frac{1+2\varepsilon_0}{2-4\varepsilon_0}\lt0.54,
\qquad \|R^\dagger/\alpha\|\le\frac12.
```

Their local generators obey the additional coupled relations

```math
\begin{pmatrix}g_v&k_v\\h_v&d_v\end{pmatrix}
=(U_v^C)^\dagger
  \mathrm{diag}(g_{2v},g_{2v+1})U_v^{W'}.
```

Together they make $`S+R`$ unitary. The current affine and reverse
dilations use these data as separate contractions. Their assembly gives a
sufficient bound

```math
T=O(T_S+T_R+N+L),\qquad
G=O(G_S+G_R+NL).
```

Thus reducing their combined cost to $`O(N+L)`$ would work. It is one
sufficient route, not a requirement on a solution. A new word may use the
coupled unitary directly and choose different rejected completions. No
constant-query conversion of an opaque forward-only block to the affine
block has been proved.

### Why a single larger group is not already the answer

In the retained grouped construction, a group of s levels above a suffix
of length r, with $`e=n-r`$, has

```math
w_g=O(\log(s+2)),\qquad
Q_g=\Theta(s2^{n-r}),\qquad
m_g=L+\lfloor r/4\rfloor+8.
```

These are private initialized work, padded coefficient-table rows, and
precision-source width. Exponentially growing groups give
$`\sum_gQ_g=O(N)`$ and $`O(\ell_*(n))`$ precision charges.
A single group spanning all but a constant suffix has
$`\Theta(nN)`$ current table rows and an $`O(nN^2)`$ lookup
Clifford estimate at $`L=N`$. Its current streamed coarse-program
estimate also grows to $`O(nN)`$ under that representation's uniform
star-normalization condition. Its current term label needs
$`\Theta(\log(n+2))`$ private initialized bits, which a constant
suffix cannot supply. Removing only its depth label or private
workspace does not establish a linear-T theorem. The expanded Clifford
estimate separately misses the stronger joint goal; it is not a T lower bound.
These are costs of the current representation, not gate lower bounds.
The selected dirty allocation also falls below the proved banked
refinement's sufficient threshold.

The weighted representation can retain the cheap coarse frame instead,
but its downward path factors still use target transport. Calling the
target frame to supply those factors would be circular; replacing it by
coarse transport alone misses mixed corrections. The selected joint goal therefore requires controlling table expansion
and coarse programs as well as precision cost; a T-only improvement must
still charge those operations but may have a larger Clifford count.

## What the failure diagnostics actually rule out

None of these statements makes separate source or mask costs additive for
a jointly synthesized circuit or settles the unrestricted endpoint.

| Shortcut | Established limit and scope |
|---|---|
| Make the same exact source sublinear | Minimum exact T-count on $`m\ge2`$ dirty wires is $`2m-4`$, or $`2m-2`$ when controlled, allowing returned helpers; [primitive bound only](OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit) |
| Hoist the source basis and use cheap masks | Valid transformed masks have linear exact/fine-accuracy cost; the actual paired source has a mask requiring $`T\ge q/2-6`$ at error $`2^{-q}`$ with full return; [separate-mask bound only](SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks) |
| Carry a source code and renew it after each query | Width changes are cheap, but a legal scalar query leaves the flag-correlated code by constant norm; complete syndrome renewal costs $`\Omega(L)`$ at compilation accuracy; [specified interface only](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost) |
| Multiply accepted blocks sharing flags | Rejected components return coherently; [full word required](SOURCE_REUSE_LIMITS.md) |
| Telescope one global conjugator through a fork | Cancellation is valid, but the remaining native word has rejected returns; a coefficient-ellipse bound excludes even arbitrary mask retuning of that word at fine precision; [fork audit](SOURCE_REUSE_LIMITS.md#5-a-shared-conjugator-does-not-close-a-branching-fork) |
| Connect coupled blocks through their zero-defect word | The accepted error is exactly $`-3(A-I)(B-I)/8`$; the [complete commutator repair](RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost) cancels to the original fine-precision wrappers; its small norm does not suppress the whole merge's local synthesis error |
| Unload a query after changing its address | The inverse can leave dirty-dependent action; [query-scheduling counterexample](SOURCE_REUSE_LIMITS.md) |
| Replace initialized nilpotent source by dirty encoding | The specified full-output scalar relation needs $`\log_2L-O(1)`$ initialized width, regardless of dirty width; [that interface only](SOURCE_REUSE_LIMITS.md) |
| Compress a whole frontier to one scalar | Generic depth-d cut rank is $`2^d`$, although each edge has rank one; existing logical modes carry it, so this is [not an ancilla lower bound](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Permute a merged band to bounded-size blocks | The chosen completion has connected support growing with height; other completions or nonpermutation bases remain allowed; [fixed-completion restriction](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Truncate propagation using its norm margin | Constant Riccati messages coexist with undamped continuation; local-depth truncation misses a fixed bottom component, but the witness has a cheap global circuit; [word-specific failure](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Use fixed-order coarse transport | Mixed corrections survive at fine precision; [specified-order failure](ENDPOINT_TREE_TRANSPORT.md) |

For clarity, the rejected-return term is an exact operator identity. With
$`P=JJ^\dagger`$ and $`B_i=J^\dagger Q_iJ`$,

```math
J^\dagger Q_2Q_1J
=B_2B_1+J^\dagger Q_2(I-P)Q_1J.
```

For $`Q_1=Q_2=H\otimes H`$ on two flags, the separate accepted
blocks are $`1/2`$, yet the product's accepted block is one, not
$`1/4`$. A successful global word may exploit rejected returns, but may
not discard them or assume uncharged initialized history prevents them.
Similarly, the current chosen forward completion need not approach its
zero-defect completion when its accepted map becomes small.

The complete fusion hierarchy retains every input column: a parent and two
nonterminal children form a ten-mode unitary, and a closed subtree with m
internal nodes uses $`2m`$ modes. These are basis modes, not additional
clean wires. Its linear list of local scatterers is compact, whereas eager
entrywise expansion of an h-level affine map has $`h2^h+1`$ generic
nonzero entries. Neither linear classical storage nor dense expansion
supplies the missing native synthesis by itself.

## Revision decision and next bounded pass

The coupled-merge passes have completed their structural task. The
[commutator repair](RESIDUAL_ASSEMBLY.md#9-a-repair-word-without-an-ill-conditioned-transported-basis)
works on every signal port without normalizing transported differences.
Its actual inverse calls cancel even for a noncanonical approximate child
word. The reduced circuit retains one child call and the original local
target wrappers. The complete merge has parent error
$`\|\widehat B-B\|/2`$, with no small residual factor, so the repair
does not discount their precision. Its layerwise emission remains
$`O(N+nL)`$. The separate native shared-conjugator fork fails even
after arbitrary mask retuning. Another proof of these boundaries or
cancellations would not improve the resource frontier.

The source-carry pass tested **precision reuse between adjacent existing
groups of the best grouped compiler**. It resolves the width-only part of
the proposed transition, while exposing the missing program cost. A fixed
fork or a fixed number of fused groups still changes only constants; the
needed gain must survive a variable number of unequal groups.

### What the source-carry pass established

The [native width analysis](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
keeps literal phases and proves full-space identities, including occupied
flags and correlated dirty inputs.

| Candidate | Result | Consequence |
|---|---|---|
| Original chain loader, extracted from each group | The required one-bit eigenbasis bridge has exact T-count $`2m-1`$ and needs at least $`L-3`$ T gates at error $`2^{-L}`$ in the stated width range | The cheap opposite-order product is not this bridge |
| Reverse-order star loader | Same certified coefficient grid; monotone source boundaries cost exactly $`2(m_{\max}-1)`$ T gates and $`O(Rm_{\max})`$ Clifford gates | Width changes are solved for this choice, but transformed group programs are excluded from that count |
| Separately synthesized transformed mask | A realizable grouped coefficient $`1/(9s)`$ requires $`T\ge\max\{0,m-2\log_2s-9\}`$ at error $`2^{-m}`$ | The old table-row cost cannot simply be assigned to masks in the new basis; this is not an additive compiler lower bound |
| Flag-correlated chain-source code | Preparation is charged; width changes cost their size difference. A realizable $`c=1/8`$ query has leakage norm $`\sqrt7/4`$ | A source cannot remain a cheap flag X after that query without a changed boundary |
| Complete syndrome renewal | Needs $`T\ge(L-4)/2`$ at compilation accuracy under its stated width condition | This particular renewal interface reintroduces precision cost; code-restricted or deferred alternatives remain open |

The [flag-correlated proof](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost)
does not assume the released source tail becomes clean. It also does not
supply a complete group program on the code. The retained generic resource
frontier is unchanged.

### Revised sufficient ledger: jointly compile the interior program

Keep the groups and precision allocation of
[conditional-suffix Sections 6--7](CONDITIONAL_SUFFIX_COMPILER.md#6-exponentially-growing-groups-and-the-explicit-workspace-ledger).
For $`n\leq r_0`$, the fixed-depth fallback already costs $`O(N+L)`$.
For each remaining group, define

```math
\begin{aligned}
Q_g&=\Theta(s_g2^{n-r_g}),& m_g&=L+\lfloor r_g/4\rfloor+8,\\
w_g&=O(\log(s_g+2)),& k_g&=n-r_g+O(\log(s_g+2)).
\end{aligned}
```

Q counts table rows and the same order of coarse-program work; w is
private initialized work inside the active suffix, and k counts dirty
selectors. The existing sums are

```math
\begin{aligned}
\sum_gQ_g&=O(N),&\sum_gQ_gm_g&=O(NL),\\
R&=O(\ell_*(n)),& R&\leq n,\\
m_{\max}&=L+O(n),&\sum_g|m_{g+1}-m_g|&=O(n).
\end{aligned}
```

The new loader realizes the last line's width ledger. What remains is a
**joint interior circuit**, including transformed programming, with total
T-count bounded by

```math
O\!\left(\sum_gQ_g+L+RP(n)\right),
```

where P is a fixed polynomial independent of L. The bound permits a
charged $`O(L)`$ use *inside* that joint circuit. Dirty-only outer
conjugation leaves the logical error contract unchanged, by
[the source-hoisting argument](SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks).
Thus its boundary loaders alone cannot supply precision to independently
accurate, precision-free group bodies. The target is a global cost bound;
it is not a claim that every group separately costs only $`O(Q_g)`$.

Combining such a proved interior with the known boundary ledger would give

```math
T=O\!\left(\sum_gQ_g+L+m_{\max}+RP(n)\right)=O(N+L).
```

This remains an **unproved sufficient construction**, not a theorem or a
necessary architecture. A changed encoding correlated with logical data
may avoid separate group outputs entirely. Leave the fixed deepest-layer
tail with its existing $`O(N+L)`$ compiler. Retaining $`G=O(NL)`$
requires separately charging all interior Clifford gates; the new loader
and bridge costs themselves fit that budget.

### Small examples: completed results and remaining cost

The small-example pass supplied the following reusable results. Their
linked chapters retain the derivations and native ledgers.

| Result | Established capability | Remaining cost or limitation |
|---|---|---|
| [Fixed-address products](SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis) | Arbitrarily many noncommuting one-qubit factors compress into four quaternion coordinates and compile once | Changing the quantum address or growing logical support invalidates that fixed two-mode table |
| [Cayley residual](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) | Complete complex-coarse residual, linear classical data, two-dimensional local updates, stable inverse conversion | Generic coherent evaluation is unpriced; direct Woodbury emission returns to the fine target wrappers |
| [Four-mode native benchmark](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) | Two magic-basis factors give $`T=O(2^k+L)`$ and $`G=O(2^kL)`$ at an unchanged prefix, with one clean flag and the stated borrowed spectator | A fixed-support corollary; the separately charged complex coarse inverse remains necessary |
| [Eight-mode coupling](ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling) | Explicit controlled Bell-projector rotation and a four-Pauli comparison | Independent prefix/child factors cannot absorb the root; this does not exclude joint synthesis |
| [Shared changing-target word](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) | Exact common scalar conjugation, full SU(2) signal error, and two parity boundaries returning an arbitrary borrowed signal | Fixed-mask simplification gives the comparison words the same leading precision cost; variable-depth allocation and joint programming remain unproved |

The last construction reduces declared source counts from $`45g`$ to
$`27g+6`$: 135 to 87 at eight modes and 180 to 114 at sixteen. Its
fixed mask commutes with the loader tail, leaving a two-T seed; each
hoisted fixed-A interior then costs at most eight T gates. Both words
retain five programmed B words per stage and the same leading term
$`(40g+4)q`$. The fixed-A remainder changes from 80g to
$`8(4g+2)`$ before further optimization, with queries and controls
separately charged. These are declared-word bounds, not minima.

For local depth three or four, its stated allocation fits
$`a=2,b=L+n+7`$ and gives $`T=O(2^k+L)`$ at a k-bit unchanged
prefix. The source-width margin grows with variable depth. Neither the
fixed-depth cost nor its full borrowed-signal return narrows the generic
endpoint gap. The source-appearance reduction is not a leading precision
improvement or an additive lower-bound argument.

### Revision after the joint native audit

The small-example pass is complete as structural infrastructure. The
phrase "joint synthesis of the programmable masks" identifies a missing
result, not yet a construction. Continuing to simplify two ordinary
layers without specifying how the saving scales risks another constant
improvement to a weaker baseline.

| Comparison | What is already available | What would constitute progress |
|---|---|---|
| Fixed logical support | A fixed number of modes costs $`O(L)`$; one unchanged-address two-mode product compiles once | A representation whose charged program remains controlled as the support grows |
| Ordinary layers | The retained baseline costs $`O(N+nL)`$; the shared body changes its displayed constants | A precision recurrence that improves the best grouped construction, or a separate complete construction beating it |
| Existing unequal groups | $`T=O(N+LR)`$, $`R=O(\ell_*(n))`$, with linear total table work | A jointly emitted program with one global precision charge and a valid allocation throughout |

The paired-source real-Y word is not a replacement for a grouped scalar
SELECT. The latter keeps forward and actual-inverse scalar branches,
four complex phases, selected column maps, a private term label, and a
reflection on the initialized active suffix. Its source width and suffix
predicate change between groups. Even the one-clean extension retains
this scalar/atom separation; it relocates a flag and changes the deepest
tail compiler. A transfer from the small-example word to these grouped
interfaces has not been proved.

Keep the recent proofs and fixtures; do not add another fixed-size
optimization as an endpoint advance. Cheap width transport is available,
fixed-mask simplification is accounted for, and the remaining general
lower bound is still only the retained one. A failure of a selected
architecture would not make the grouped factor necessary.

### Next bounded task and stopping rule

The next pass is **construction selection at the actual group interface**.
Before another numerical mask search, specify one native identity or
changing encoded boundary that could remove repeated precision from the
grouped forward/reverse SELECT. State its cost hypothesis for a variable
number R of groups, not only two ordinary layers. For the sufficient
route already defined above, the target is

```math
T_{\rm joint}=O\!\left(\sum_gQ_g+L+RP(n)\right),
```

with P a fixed polynomial independent of L. Equivalently, a proposed
ledger of the form

```math
T_{\rm joint}\le c_RL+O\!\left(\sum_gQ_g+RP(n)\right)
```

must bound $`c_R`$ independently of R for this sufficient route. This
is a proposed budget, not an attained bound. A smaller constant multiplying
LR does not meet it. Table and coarse-program work must combine additively,
not through a Cartesian product of the group addresses. Other complete
constructions may use a different representation and are judged against
the same endpoint, rather than required to adopt SELECT.

The first specific compatibility test uses the retained nonnegative
coefficients $`c\in[0,1/4]`$. On the existing scalar rejection flag sigma,
consider the ideal completion

```math
R_c=\exp[-i\arccos(c)Y_\sigma]
=\begin{pmatrix}c&-\sqrt{1-c^2}\\\sqrt{1-c^2}&c\end{pmatrix}.
```

Its accepted entry is c. Replace the scalar factor in one forward atom
$`\mathcal S_c\mathcal D_\nu`$ by this rotation, and use the actual
inverse for the reverse atom. Keep the distinct atom flag and external
phase components. This is a proposed native substitution, not a free
rotation oracle or a proved grouped improvement.

Use a fully amplified native approximation to $`R_c`$, not its
unamplified scaled block. Derive a fresh full-operator perturbation ledger
against the ideal canonical group word; the previous exact rounded-scalar
identity is not inherited. An inner error of order $`2^{-m_g}`$ would
require a sufficient constant shift $`q_g=m_g+O(1)`$, which must be
included in the workspace audit. A structural zero coefficient means
$`R_0=\exp(-i\pi Y_\sigma/2)`$, not identity. Only the inactive
suffix sector is required to have exact identity action.

Use the permitted two-clean grouped layout: h and sigma are external;
the new rotation-synthesis signal would be in the arbitrary dirty pool,
outside the initialized-work reflection. Protect sigma from the coarse
interpreter's temporary buffers. Audit whether every remaining
atom, phase, SELECT, reflection, and coarse subroutine commutes with
$`Z_\sigma`$ as a complete returned-work word. If so, the shared-body
parity argument has a specific route to this interface. Specify certified
coefficient preprocessing separately, as in the retained compiler. Charge
every emitted rotation program and the extra borrowed signal/helper
reservation; any constant increase must fit the group slack. Do not assume
that dirty helpers preserve this commutation gate by gate.

Name the source family and table representation explicitly. The cheap
reverse-order one-tail width bridges do not automatically apply to the
paired-source replacement. Record private suffix work, both rejection
flags, term label, dirty core, selectors, query unloading, and reflection
in a register-lifetime table. Explain what can persist when the next
suffix predicate is computed. This focused symbolic derivation does not
require a full compiler implementation or large matrices.

Even a successful compatibility test still pays for the programmed
rotations. Continue to two unequal groups only after specifying a fusion
rule that could change $`c_R`$. Track their widths symbolically; small
finite analogues and illustrative group constants do not certify the
asymptotic allocation thresholds.

A proposed carried boundary may correlate with the logical data and need
not expose an accurate standalone output after every group. It must give
a complete transition rule and a charged final decoding. Merely defining
that boundary using the target frame or its partial products is circular
unless their native implementation is supplied. If the rule closes for
two unequal groups, apply that same rule to a third and derive its
variable-R recurrence before adding a new resource claim. Retain actual
complex coarse words, zero defects, independent branch angles, and every
signal/dirty port in the small checks.

If no such identity or boundary can be specified, report that result and
change the mechanism. Do not replace it with another generic mask sweep,
fixed-depth synthesis, or repeated audit of a known failed interface. A
focused primary-literature check is useful once a concrete alternative
has identifiable assumptions; an untargeted expansion of the literature
chapter is not the next deliverable.

The live-wire contract remains $`b=L+n+7`$ arbitrary dirty wires and
at most two external clean flags. Holding $`m_{\max}`$ throughout
cannot be assumed to fit uniformly in n: together with the deepest
group's $`n-r_0+O(1)`$ selectors it can exceed the allocation. A proposed
schedule must prove valid tail release or reassignment, or another valid
allocation. The proved transitions permit arbitrary tail correlations,
not initialized tails.
Changing suffix predicates retains only its existing active-sector
$`w_g`$ reservation. Include actual query unloading, literal inverses,
complete dirty/reference return, and all transition approximation errors.
If using an accepted-block construction, preserve normalization two and
account for rejected action in amplification; a direct complete-frame
isometry proof is also allowed.

Stop a candidate if its expanded word still pays a length-L source or
separately synthesized mask per group, requires full syndrome renewal per
query, expands the total table/coarse-program work beyond $`O(N)`$, or
uses an unproved clean-work or return assumption. These stop conditions
apply to this sufficient route; they are not general lower bounds. Record
failure in the existing proof home and change the native candidate. Do not
repeat the completed width, mode-closure, or commutator-repair audits.

Small fixtures can suggest or falsify a recurrence; they cannot prove its
asymptotic cost. No target-frame oracle, free evaluator, reset, initialized
history, or supplied catalyst is available. A bounded classical inverse
condition number does not make its quantum implementation free.

Separately improving the affine and reverse blocks remains sufficient,
and a linear-T endpoint circuit with larger fully charged Clifford cost
would still settle the T-only question. No joint interior with one global
precision charge has yet been constructed. The latest examples supply an
exact native word and a fair simplified baseline; they do not show that
the generic gap is close to resolution.

Optimal depth, practical constants, and a full elementary emitter are
separate tasks. The established publication scope is unchanged; manuscript
writing and release work are outside this pass.

## Evidence and remaining implementation work

Internal review of source conditioning, reverse-word order, workspace, and
error sums has not identified a defect; this is not independent peer
review. The [verification map](VERIFICATION.md) separates analytic proofs
from finite checks of native sources, literal phases, inverses, rejected
returns, support packing, and small complete dirty-input blocks. Tests do
not establish asymptotic theorems, optimality, or literature priority.

The strongest grouped construction has no end-to-end elementary Clifford+T
emitter. Its table fixture uses the table action directly, grouping tests
use illustrative constants, and some workspace constants and crossover
thresholds remain existential. Selected-atom pseudocode and a register-lifetime
table would improve auditability before practical resource estimates.
Proof chapters retain the arguments; this checkpoint records the next decision.
