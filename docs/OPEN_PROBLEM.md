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
families. The subsequent merge passes establish complete boundary actions,
an explicit repair, and scoped failures of proposed source sharing.
**None has narrowed the generic upper/lower gap.** Wide independent changes
and deep sparse nesting each admit one precision charge; the remaining
question concerns precision reuse for unrestricted comparable updates.

| Level | Established scope |
|---|---|
| General theorems | Exact matching resources, sufficient-clean matching T-count, and one-clean grouped bound; constant-clean endpoint and optimal T-depth remain open |
| Incomparable promised families | Antichain and sparse ancestor-closed updates have $`T=O(N+L)`$, $`G=O(NL)`$; generic rounding satisfies neither promise |
| Reusable blocks | Cheap coarse frame, linear residual data, complete weighted dilations, two-flag assembly, and a coupled repair with exact child-call cancellation; the surviving target precision remains charged |
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

The next bounded task is **precision reuse between adjacent existing
groups of the best grouped compiler**. A fixed two-depth fork already
admits $`O(L)`$ synthesis. Saving a constant number of precision charges
there, or fusing a fixed number of existing groups, would only improve
constants. The needed gain must survive a variable number of groups of
unequal, growing heights.

### A sufficient transition lemma to seek

Use the groups and precision allocation of
[conditional-suffix Sections 6--7](CONDITIONAL_SUFFIX_COMPILER.md#6-exponentially-growing-groups-and-the-explicit-workspace-ledger).
For $`n\leq r_0`$, the existing fixed-depth fallback already costs
$`O(N+L)`$. Otherwise there is at least one group.
For group g of height $`s_g`$ above $`r_g`$ suffix bits, put

```math
\begin{aligned}
Q_g&=\Theta(s_g2^{n-r_g}),& m_g&=L+\lfloor r_g/4\rfloor+8,\\
w_g&=O(\log(s_g+2)),& k_g&=n-r_g+O(\log(s_g+2)).
\end{aligned}
```

Here Q counts padded coefficient-table rows and the same order of native
coarse-program work; m is source width, w is private initialized work
inside the active suffix, and k counts dirty query selectors. The retained
construction has $`\sum_gQ_g=O(N)`$ and
$`R=O(\ell_*(n))`$ groups, with $`R\leq n`$. Its excess precision cost is
$`\sum_gm_g=O(L\ell_*(n)+N)`$.

The proposed lemma would give these groups a common, explicitly defined
boundary action and a literal native transition between adjacent groups.
It must retain their $`O(Q_g)`$ table and coarse-program costs and their
private-work budgets, while charging initial and final source work only
$`O(m_{\max})`$ in total and each transition only

```math
O\!\left(|m_{g+1}-m_g|+P(n)\right)
```

T gates. P is a fixed polynomial, independent of L and group height; it
may cover predicates and routing at a boundary. All other non-source
group work must remain $`O(Q_g+P(n))`$. This is an **unproved sufficient
mechanism**, not a new compiler theorem or a necessary form of every
solution.

The intended accounting is concrete. In the prescribed shallow-to-deep
execution order the $`r_g`$, and hence the $`m_g`$, are monotone.
Therefore

```math
m_{\max}=L+O(n),\qquad
\sum_{g=1}^{R-1}|m_{g+1}-m_g|=O(n).
```

If the transition lemma holds for the actual expanded words, their cost is

```math
\begin{aligned}
T&=O\!\left(\sum_gQ_g+m_{\max}
 +\sum_{g=1}^{R-1}|m_{g+1}-m_g|+RP(n)\right)\\
 &=O(N+L).
\end{aligned}
```

The last step uses $`R\leq n`$ and the fixed degree of P. Leave the
fixed $`r_0`$ deepest layers with their existing $`O(N+L)`$ compiler;
their different source need not be fused into this interface. The retained
table sum $`\sum_gQ_gm_g=O(NL)`$ already prices its Clifford work.
Additional transition and source Clifford gates must also be charged;
$`G=O(NL)`$ remains the stronger joint goal, not a consequence of the
T-count transition estimate alone.

### The boundary must include the actual live workspace

Simply keeping a maximum-width source live through every group is not a
valid general allocation. The deepest grouped band has fixed suffix
length $`r_0`$ and needs $`n-r_0+O(1)`$ selectors. Combining these
with $`m_{\max}=L+\lfloor r_{\max}/4\rfloor+8`$ exceeds
$`b=L+n+7`$ as $`r_{\max}`$ grows. The transition must explicitly
release or reassign source-tail wires as selector demand grows. Such
reassignment may carry arbitrary correlations; it cannot assume that a
source tail or a reused core has been reset or returned independently.
The conditional suffix is available only under its existing active-sector
promise, with its existing $`w_g`$ reservation.

The first deliverable is the full shared boundary contract, a parametric
native transition for two unequal groups, and the resulting symbolic
ledger. Establish closure under a third group and arbitrary further
merges **before adding new small matrix tests**. State every live register,
actual inverse, query address and unloading operation, and how the complete
word controls dirty return and reference correlations. At most two external
clean flags and $`L+n+7`$ arbitrary dirty wires are available. A fused
accepted-block interface must retain normalization two and include its
rejected action in the amplification argument; alternatively prove the
complete-frame isometry contract directly. The final error must include
any transition approximations, not only the old per-group digit errors.

The existing source-hoisting and changed-address restrictions still apply.
Moving a loader outside the word does not make its transformed masks cheap,
and a changed logical address does not unload a prior query. The stopped
shared-conjugator word requires a changed construction. No target-frame
oracle, free evaluator, intermediate reset, initialized history, or padded
zero sector is supplied. Subsequent finite checks should test the new
native transition and its inverse on complete dirty inputs, including a
real target with a close complex native coarse frame; checking the already
proved dense merge identities again would not test precision reuse.

Stop the proposed transition if its expanded circuit recharges a length-L
source per group or merge, expands the total table/program work beyond
$`O(N)`$, or loses its normalization, error, unloading, or live-workspace
contract. Record that failure in the existing proof home. These are
go/no-go criteria for this sufficient mechanism, not lower bounds against
other compilers. A different proved ledger may still succeed, and a valid
$`O(N)`$-T endpoint circuit with a larger fully charged Clifford count
would still close the T-only question.

Separately improving the affine and reverse blocks remains sufficient.
The current generic gap has not narrowed; resolving it still requires a
precision-amortization identity, not another structural completion.

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
