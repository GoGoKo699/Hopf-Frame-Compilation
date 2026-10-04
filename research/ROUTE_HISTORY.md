# Endpoint route assessments

These are retained assessments from completed research passes. Their proposed next steps are historical, not current work orders. Consult the [research index](README.md) for outcomes and the [current gap statement](../docs/OPEN_PROBLEM.md) for the active boundary. Repeated frontier summaries and the completed state-based protocol history have been removed from this log; their authoritative homes are the selected claim ledger and the completed supplement.

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
[normalization-two amplification lemma](../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
then gives the complete-isometry compiler with three charged calls and the
same two flags.

The [two-flag assembly](endpoint/RESIDUAL_ASSEMBLY.md) realizes this interface at
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
| Affine/reverse assembly | [Two-flag proof](endpoint/RESIDUAL_ASSEMBLY.md): diagonal-plus-forward and inverse-based reverse blocks give normalization two, with controls/inverses charged; cost still contains $`nL`$ |
| Full-port hierarchy | [Fusion audit](endpoint/RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths): recursively closed complete unitary, linear classical generators; faster native synthesis unproved |
| Coupled completion | [Whole-residual boundary](endpoint/RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost): one-signal normalization-two recursion; a commutator implements the complete rank-at-most-four repair without normalizing transported differences; its child calls cancel to the original target wrappers, whose precision remains charged |
| Grouped frame | Best general endpoint bound, with conditional suffix and core return in the complete error |
| Product-first residual coordinates | [Cayley recursion](endpoint/ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) gives constant-size local data for the complex-coarse residual; the [four-mode native benchmark](endpoint/ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) uses two one-target programs, while the eight-mode controlled coupling and general coherent conversion remain charged |
| Source-width transport | [Reverse-order loader](endpoint/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation): total boundary T-count $`2(m_{\max}-1)`$ with linear native loaders; transformed group bodies remain charged |
| Joint source body | [Changing-target word](endpoint/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) shares one fixed scalar conjugator and returns a borrowed signal through two parity boundaries; the eight-/sixteen-mode words are priced, but fair fixed-mask simplification leaves the same leading precision cost in both comparisons |

The [weighted norm proof](endpoint/ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
avoids a height penalty via uniform coarse subtree accuracy. The
[forward block](endpoint/WEIGHTED_TRANSPORT_BLOCK.md) charges marker gathering,
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
| Make the same exact source sublinear | Minimum exact T-count on $`m\ge2`$ dirty wires is $`2m-4`$, or $`2m-2`$ when controlled, allowing returned helpers; [primitive bound only](../docs/OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit) |
| Hoist the source basis and use cheap masks | Valid transformed masks have linear exact/fine-accuracy cost; the actual paired source has a mask requiring $`T\ge q/2-6`$ at error $`2^{-q}`$ with full return; [separate-mask bound only](endpoint/SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks) |
| Carry a source code and renew it after each query | Width changes are cheap, but a legal scalar query leaves the flag-correlated code by constant norm; complete syndrome renewal costs $`\Omega(L)`$ at compilation accuracy; [specified interface only](endpoint/SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost) |
| Multiply accepted blocks sharing flags | Rejected components return coherently; [full word required](endpoint/SOURCE_REUSE_LIMITS.md) |
| Telescope one global conjugator through a fork | Cancellation is valid, but the remaining native word has rejected returns; a coefficient-ellipse bound excludes even arbitrary mask retuning of that word at fine precision; [fork audit](endpoint/SOURCE_REUSE_LIMITS.md#5-a-shared-conjugator-does-not-close-a-branching-fork) |
| Connect coupled blocks through their zero-defect word | The accepted error is exactly $`-3(A-I)(B-I)/8`$; the [complete commutator repair](endpoint/RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost) cancels to the original fine-precision wrappers; its small norm does not suppress the whole merge's local synthesis error |
| Unload a query after changing its address | The inverse can leave dirty-dependent action; [query-scheduling counterexample](endpoint/SOURCE_REUSE_LIMITS.md) |
| Replace initialized nilpotent source by dirty encoding | The specified full-output scalar relation needs $`\log_2L-O(1)`$ initialized width, regardless of dirty width; [that interface only](endpoint/SOURCE_REUSE_LIMITS.md) |
| Compress a whole frontier to one scalar | Generic depth-d cut rank is $`2^d`$, although each edge has rank one; existing logical modes carry it, so this is [not an ancilla lower bound](endpoint/RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Permute a merged band to bounded-size blocks | The chosen completion has connected support growing with height; other completions or nonpermutation bases remain allowed; [fixed-completion restriction](endpoint/RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Truncate propagation using its norm margin | Constant Riccati messages coexist with undamped continuation; local-depth truncation misses a fixed bottom component, but the witness has a cheap global circuit; [word-specific failure](endpoint/RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Use fixed-order coarse transport | Mixed corrections survive at fine precision; [specified-order failure](endpoint/ENDPOINT_TREE_TRANSPORT.md) |

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

This is a retained endpoint-route decision from an earlier pass. Its
construction candidates remain conditional; the current decision follows
[the completed audit](../WORKSPACE.md#research-decision-and-stopping-rules)
and freezes the selected scientific scope.

The current selection is the
[revision checkpoint above](../docs/OPEN_PROBLEM.md#current-decision).
This section retains the earlier source-carry decision and its boundaries.

The coupled-merge passes have completed their structural task. The
[commutator repair](endpoint/RESIDUAL_ASSEMBLY.md#9-a-repair-word-without-an-ill-conditioned-transported-basis)
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

The [native width analysis](endpoint/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
keeps literal phases and proves full-space identities, including occupied
flags and correlated dirty inputs.

| Candidate | Result | Consequence |
|---|---|---|
| Original chain loader, extracted from each group | The required one-bit eigenbasis bridge has exact T-count $`2m-1`$ and needs at least $`L-3`$ T gates at error $`2^{-L}`$ in the stated width range | The cheap opposite-order product is not this bridge |
| Reverse-order star loader | Same certified coefficient grid; monotone source boundaries cost exactly $`2(m_{\max}-1)`$ T gates and $`O(Rm_{\max})`$ Clifford gates | Width changes are solved for this choice, but transformed group programs are excluded from that count |
| Separately synthesized transformed mask | A realizable grouped coefficient $`1/(9s)`$ requires $`T\ge\max\{0,m-2\log_2s-9\}`$ at error $`2^{-m}`$ | The old table-row cost cannot simply be assigned to masks in the new basis; this is not an additive compiler lower bound |
| Flag-correlated chain-source code | Preparation is charged; width changes cost their size difference. A realizable $`c=1/8`$ query has leakage norm $`\sqrt7/4`$ | A source cannot remain a cheap flag X after that query without a changed boundary |
| Complete syndrome renewal | Needs $`T\ge(L-4)/2`$ at compilation accuracy under its stated width condition | This particular renewal interface reintroduces precision cost; code-restricted or deferred alternatives remain open |

The [flag-correlated proof](endpoint/SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost)
does not assume the released source tail becomes clean. It also does not
supply a complete group program on the code. The retained generic resource
frontier is unchanged.

### Revised sufficient ledger: jointly compile the interior program

Keep the groups and precision allocation of
[conditional-suffix Sections 6--7](../docs/CONDITIONAL_SUFFIX_COMPILER.md#6-exponentially-growing-groups-and-the-explicit-workspace-ledger).
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
[the source-hoisting argument](endpoint/SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks).
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
| [Fixed-address products](endpoint/SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis) | Arbitrarily many noncommuting one-qubit factors compress into four quaternion coordinates and compile once | Changing the quantum address or growing logical support invalidates that fixed two-mode table |
| [Cayley residual](endpoint/ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) | Complete complex-coarse residual, linear classical data, two-dimensional local updates, stable inverse conversion | Generic coherent evaluation is unpriced; direct Woodbury emission returns to the fine target wrappers |
| [Four-mode native benchmark](endpoint/ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) | Two magic-basis factors give $`T=O(2^k+L)`$ and $`G=O(2^kL)`$ at an unchanged prefix, with one clean flag and the stated borrowed spectator | A fixed-support corollary; the separately charged complex coarse inverse remains necessary |
| [Eight-mode coupling](endpoint/ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling) | Explicit controlled Bell-projector rotation and a four-Pauli comparison | Independent prefix/child factors cannot absorb the root; this does not exclude joint synthesis |
| [Shared changing-target word](endpoint/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) | Exact common scalar conjugation, full SU(2) signal error, and two parity boundaries returning an arbitrary borrowed signal | Fixed-mask simplification gives the comparison words the same leading precision cost; variable-depth allocation and joint programming remain unproved |

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

The unmodified paired-source real-Y word did not supply a grouped scalar
SELECT. That interface keeps forward and actual-inverse scalar branches,
four complex phases, selected column maps, a private term label, and a
reflection on the initialized active suffix. Its source width and suffix
predicate change between groups. Even the one-clean extension retains
this scalar/atom separation; it relocates a flag and changes the deepest
tail compiler. That transfer was the unresolved compatibility question at the revision.
The canonical audit below now supplies it, without improving the precision
recurrence.

Keep the recent proofs and fixtures; do not add another fixed-size
optimization as an endpoint advance. Cheap width transport is available,
fixed-mask simplification is accounted for, and the remaining general
lower bound is still only the retained one. A failure of a selected
architecture would not make the grouped factor necessary.

### The canonical group test is complete

The [canonical scalar audit](endpoint/CANONICAL_SCALAR_COMPLETION.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge)
transfers the small native rotation to the actual grouped interface. It
uses $`R_c=\exp[-i\arccos(c)Y_\sigma]`$ on the existing scalar
flag, with a separate atom flag and an arbitrary borrowed synthesis
signal outside the initialized-work reflection. The canonical accepted
entry is c, including structural zero. Forward and reverse atoms keep
their literal phases and query order.

Direction-controlled $`Z_\sigma`$ gates select the actual native
inverse of a single positive-angle program, even after rounding. The old
one-tail scalar has the same inverse symmetry, so this one-program SELECT
is available to both baselines. The canonical substitution is not uniquely
responsible for that constant improvement.

The new full-operator error ledger gives one inner program per half-block
and three per outer amplification. With $`q_g=m_g+4`$, the group error
is below $`(45/8)2^{-m_g}`$, inside its retained budget. The paired core
and borrowed signal add six dirty slots before the fixed helper constants;
the existing group slack absorbs this by increasing its fixed thresholds.
No additional external clean qubit is needed. Final h and private-work
return remain part of the complete approximation error.

The exact common-source identity also passes through the non-scalar
group operations. But its cost fails the endpoint test: at a legal common
q, the paired comparison words have the same leading precision term
$`(120R+4)q`$. The old scalar uses six programmable mask appearances
per group, and the canonical replacement uses thirty. Unequal widths do
not acquire a free bridge. The proved totals remain
$`T=O(N+L\ell_*(n))`$ and $`G=O(NL)`$.

**Close this canonical-completion candidate for precision amortization.**
It supplies a correct group interface, exact inverse routing, and a
complete error/workspace ledger. It does not supply a variable-group
fusion rule or evidence that the endpoint is close. Further mask sweeps,
larger instances, and constant cancellation in this word are not the
next task. This is a conclusion about the displayed construction, not
an unrestricted lower bound.

### Selection audit after canonical completion

The 2026-10-01 revision checks three alternatives before selecting another
construction. The comparison is at the actual endpoint, where the accuracy
bit count is $`L=N`$; factors polynomial in L cannot be hidden in a
soft-O estimate.

| Mechanism | What it supplies | Why it is not yet an endpoint construction |
|---|---|---|
| Whole-residual Hamiltonian synthesis | A generic route from a Frobenius-small residual to a complete unitary channel | The explicit bound below still charges precision polynomially and uses initialized work; a tree-specific adaptation is missing |
| State preparation followed by Householder reflections | A prepared state determines its rank-one reflection, irrespective of the preparer's other columns | The prescribed frame needs more than its first column; neither constant reflection count nor a two-clean implementation follows |
| Reduce the shared-source proof's clean work | The sufficient-clean construction already charges precision once globally | Its source column, private SELECT labels/flags, and failure history each use initialization; dirty-safe lookup banks do not replace those roles |

**Small residual norm is already affordable.** The
[Frobenius-coarse specialization](endpoint/ENDPOINT_TREE_TRANSPORT.md#frobenius-small-residuals-also-fit-the-endpoint-budget)
uses no clean qubits and returns the dirty pool exactly, with
$`T(C)=O(n\sqrt N)=O(N)`$ and $`G(C)=O(Nn)`$, while making
$`\|C^\dagger W-I\|_F\le\varepsilon_0`$ for fixed
$`0\lt\varepsilon_0\le1/64`$. This is a specialization of the
existing borrowed compiler, not a new frame bound. It does not price a
coherent logarithm or matrix-function evaluator.

[Fang–Heunen–Wang, Corollary 3.9](https://arxiv.org/html/2607.12907v1)
give, for bounded Frobenius distance to the Clifford group and accuracy
$`\varepsilon`$,

```math
\begin{aligned}
T&=O\!\left((N+\log\log(1/\varepsilon))
             (n+\log(1/\varepsilon))^2\right),\\
a&=O\!\left(N(n+\log(1/\varepsilon))\right).
\end{aligned}
```

At $`\log(1/\varepsilon)=\Theta(N)`$, these displayed bounds are
$`O(N^3)`$ T gates and $`O(N^2)`$ initialized ancillas. Definition 1.1
traces out that initialized environment; it does not certify literal phase
and arbitrary dirty-work return. This screens out direct substitution of
that theorem, not a structure-sensitive adaptation or an optimality claim.

[Gosset–Kothari–Wu, Theorem 1.1 and Section 1.1](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/)
give optimal state preparation with initialized ancillas and the
Householder-based K-column upper bound

```math
T=O\!\left(K\sqrt{N\log(K/\varepsilon)}
             +K\log(K/\varepsilon)\right).
```

Its direct complete-column specialization $`K=N,L=N`$ is $`O(N^2)`$,
before any clean-work adaptation. There is also a simple structural check
on the proposed constant-reflection shortcut: set every upper-tree angle
to zero and all bottom pair angles to one small nonzero value. Then
$`\mathrm{rank}(W-I)=N`$, whereas a product of k rank-one
reflections differs from identity by rank at most k. For accuracy below
the smallest singular value of $`W-I`$, the same witness excludes such
a k<N approximation. This elementary rank argument concerns that literal
reflection product; packed higher-rank reflections and fully charged
changes of basis remain allowed. The witness itself is an antichain and
already has a linear compiler.

The clean-work audit is likewise specific. In
[shared-source Sections 7–9](../docs/FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel),
the prepared source supplies the full-output relation, the private labels
and separate zero flags define the selected contraction, and a known-zero
counter prevents rejected paths from returning to the accepted block.
Changing the counter alone leaves the other initialized roles intact.
The existing nilpotent-source restriction remains scoped to its stated
interface; it is not a full-kernel or compiler lower bound.

These checks select no construction. They focus the next question on
preserving the tree structure in the native program, rather than merely
making the residual small, preparing its state column, or renaming clean
work as borrowed work.

### The Hopf-specific scattering test

The [packed scattering construction](endpoint/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
uses the actual Hopf tree. One n-bit address selects its node rotation;
one additional mode qubit accommodates all internal continuation ports.
Explicit port permutations produce a unitary step

```math
\mathscr S=\begin{pmatrix}D&C_{\rm out}\\ B&A\end{pmatrix},
\qquad W=D+C_{\rm out}(I-A)^{-1}B.
```

The second identity gives every prescribed real-frame column, including
singular angles. The genuine internal block is nilpotent; two decoupled
dummy ports are padded by minus identity, so the padded A is not nilpotent.
The inverse in this algebraic identity is not a free circuit operation.

One step is nevertheless fully priced: the returned-work port permutations
cost $`O(n^4)=O(N)`$, and the retained real-rotation compiler, with four
invariant address sectors, gives $`T=O(N+L)`$, $`G=O(NL)`$, at
$`b=L+n+7`$. Its synthesis uses no initialized work; one available clean
qubit selects external ports. This is a full-operator approximation to the
step, not a compiler for its feedback boundary.

The conversion cannot use only a fixed number of queries to that unchanged
local-angle step. On the all-right path, set all n angles to theta. The
frame entry is $`\sin^n\theta`$, while any k-query word in the step and
its actual inverse, with parameter-independent interleaves, has amplitude
degree at most k. For $`k\lt n`$, uniform error on the canonical interval
$`[0,\pi/2]`$ is at least $`(8e^2)^{-n}>2^{-6n}`$.
Thus the fine endpoint rules out this constant-query conversion for
$`n\ge5`$, even with arbitrary fixed basis changes and extra work.
The small [scattering checks](../tests/test_hopf_scattering.py) test the
complete transfer, dummy handling, port routing, and the path coefficient.

This closes the unchanged-coin conversion as a bounded-call candidate.
It does not establish a native T-count lower bound: global coefficient
preprocessing, target-dependent interleaves, and shared native synthesis
of many calls remain outside the query argument. The Hopf-specific
representation is useful, but it has not narrowed the endpoint gap.
