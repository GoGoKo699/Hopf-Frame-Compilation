# Research status and the open two-clean endpoint

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [Two-clean baseline](OPERATOR_SOURCE_COMPILER.md) · [Grouped refinement](CONDITIONAL_SUFFIX_COMPILER.md)

This page gathers the retained results, the limits of explored routes, and
the next construction to test. It is a research checkpoint, not an additional
compiler theorem. The [publication scope](../manuscript/PUBLICATION_SCOPE.md)
gives the complete result list and the [verification map](VERIFICATION.md)
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
| Grouped two-clean real-frame T-count | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge L+n+7`$ | Precision-uniform grouped bound without additional word banks; [conditional-suffix theorem](CONDITIONAL_SUFFIX_COMPILER.md) |
| Two-clean T-count with additional dirty banks | $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge2(L+n+7)`$ | Matches the lower bound if $`L\ell_*(n)^2\le N`$ or $`b\le N/\ell_*(n)`$; these are sufficient regimes; [banked proof](CONDITIONAL_SUFFIX_COMPILER.md#8-additional-dirty-banks-give-a-width-tradeoff) |
| T-depth with additional dirty banks | $`D_T=O(NL/b+\min\{nL+n^3,L\ell_*(n)+n^4\})`$ at $`a=2`$, $`b\ge2(L+n+7)`$, with $`T,G=O(NL)`$ | Choose between the layerwise and grouped [schedules](T_DEPTH_COMPILER.md); real frames; optimizing depth may increase T-count; no matching frontier established |
| Simultaneous T-count and T-depth | $`T=O(\sqrt{NL}+L\ell_*(n))`$, $`D_T=O(\min\{nL+n^3,L\ell_*(n)+n^4\})`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge C(L+n+7+\sqrt{NL})`$ | Same real-frame circuit, for sufficiently large fixed C; [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md); T-depth optimality remains open |

The retained frontier is the best applicable bound among these constructions.
For example, at fixed L, $`a=2`$, and $`b=L+n+7=\Theta(n)`$,
the arbitrary-budget theorem and its matching splice give
$`T^\star=\Theta(N/n)`$, sharper than the grouped $`O(N)`$ estimate.
This does not change the high-precision endpoint below. Two-clean upper
constructions also work with additional unused clean qubits, but the banked
matching comparison above is for $`a=2`$.

Every retained frame construction preserves the prescribed completion used
by Hopf QBP. Finite-precision substitution has the fixed-parameter bias and
dirty-reference guarantees of the [QBP approximation theorem](QBP_APPROXIMATION.md).
It is not a claim about differentiating the discrete synthesis algorithm.
The strongest grouped theorem is stated for real frames. The separately
proved two-clean complex magnitude extension retains $`O(N+nL)`$ at
$`b\ge L+n+8`$, or $`O(\sqrt{NL}+nL+NL/b)`$ at
$`b\ge2(L+n+8)`$; it does not inherit the real grouped improvement.

The layerwise operator-source compiler remains a dependency and a useful
fallback. Its literal-diagonal and complete one-target U(2) multiplexor
corollaries retain matching banked frontiers as independent capabilities. Earlier, weaker
endpoint bounds are superseded as frontiers; exploratory routes are
recoverable from the [archived research snapshot](../provenance/README.md#earlier-research-snapshot).

## The count and depth gaps are different

The exact ancilla-depth theorem is matching in its elementary-gate model.
The T-count frontier is also matching under its stated clean reservation.
The two-clean depth theorem supplies an explicit upper schedule; its lower
bound is substantially smaller. Complete-frame safety and the QBP error
contract are already established in all these constructions. The remaining
questions concern resource scaling.

For the exactly two-clean banked regime, put $`B_0=L+n+7`$ and assume
$`b\ge2B_0`$. Then $`q=n+2+b=\Theta(b)`$. The
[inherited depth lower bound](T_DEPTH_COMPILER.md#4-lower-bounds-and-the-remaining-depth-gap)
simplifies to

```math
D_T^\star=\Omega(1+NL/b^2).
```

Indeed, $`L/b=O(1)`$ and $`\sqrt{NL}/b\le1+NL/b^2`$.
This is an algebraic restatement of the retained count-to-depth reduction,
not a stronger lower-bound argument. The following are consequences of the
existing bounds, with $`a=2`$ throughout. Allocations must satisfy the
theorem's sufficient threshold; the square-root row now uses the
parallel-indicator construction with a sufficiently large fixed prefactor.

| Regime | Depth lower bound | Available depth upper bound | Remaining issue |
|---|---|---|---|
| Fixed L, $`b=\Theta(n)`$ | $`\Omega(N/n^2)`$ | $`O(N/n)`$ | Factor-n gap |
| Fixed L, sufficiently large $`b=\Theta(\sqrt N)`$ | $`\Omega(1)`$ | $`O(n^3)`$ with $`T=O(\sqrt N)`$ | Simultaneous count/depth capability; depth lower bound remains unmatched |
| Fixed L, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(n^3)`$ with $`T=O(\sqrt N)`$ | Extra width is not needed by this schedule; depth optimality remains open |
| $`L=N`$, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ | Serial precision cost remains |
| Selected endpoint $`L=N,b=B_0`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ from $`D_T\le T`$ | The larger-bank depth theorem does not apply |

At fixed L, sufficiently large $`\Theta(\sqrt N)`$ dirty workspace now
gives worst-case optimal-order $`T=O(\sqrt N)`$ together with
$`D_T=O(n^3)`$ in one circuit. The parallel lookup audit therefore
resolves the previous incompatibility between the retained count and depth
schedules in this regime. It does not settle the depth exponent or the
tradeoff below the new sufficient-width threshold.

T-depth permits Clifford circuits of nonzero depth between its T layers.
It is not total circuit depth or elapsed QBP execution time. Neither the
exact CNOT light-cone bound nor the source's linear exact T-count minimum
supplies an additional depth lower bound in this model.

## The remaining endpoint

The publication establishes its compiler theorems without resolving this
endpoint. Let $`T^\star_{F,\mathbb R}`$ be the worst-case
minimum T-count for the prescribed real Hopf frame under the complete-input
error contract. At

```math
a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\geq3,
```

the retained results give

```math
\Omega(N)\le T^\star_{F,\mathbb R}
\le O(N\log_2^*N).
```

Here $`\log_2^*`$ counts repeated base-two logarithms until the value is
at most one. The quantities $`a`$ and $`b`$ count initialized clean and
arbitrary dirty qubits, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ for
$`0\lt \eta\le1/64`$. Literal phases, clean-work leakage, and dirty/reference
return error are included in the same operator norm as the
[fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md#1-target-resources-and-theorem).

The [layer-by-layer operator-source compiler](OPERATOR_SOURCE_COMPILER.md) has
$`T=O(N+nL)`$ with $`b\ge L+n+7`$. It pays the precision cost at each
tree depth. Neither its dirty-bank refinement nor its literal-diagonal
corollary removes that repeated cost for a general real frame. The
[conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) instead
groups consecutive depths. Its residual is a sum of ancestor-column maps,
their adjoints, and a diagonal; the logical register supplies the table
address. Only $`O(\log(s+2))`$ temporary initialized bits are needed for
a group of $`s`$ levels. A known-zero logical suffix supplies those bits
on the active sector, while every inactive input is preserved exactly.
Exponentially growing groups need only $`O(1+\log_2^*(n+2))`$ precision
charges. The general bound is $`O(N+L[1+\log_2^*(n+2)])`$ T gates with
$`O(NL)`$ Clifford gates and $`b\ge L+n+7`$.
With $`b\ge2(L+n+7)`$, the same grouped construction also gives
$`O(\sqrt{NL}+L[1+\log_2^*(n+2)]+NL/b)`$ T gates. This banked
refinement retains the iterated-logarithm precision term at the endpoint.

With a sufficiently large $`a=\Theta(n)`$ clean reservation and
$`b=\Theta(N)`$, the shared-source compiler instead attains
$`T^\star=\Theta(N)`$ at $`L=N`$. This sufficient clean reservation is
not proved necessary.

The open question is whether two clean qubits permit a jointly charged
$`O(N)`$ construction at the explicit allocation above, or whether a stronger
general lower bound holds. The current upper bound does not establish the
same cost for fewer clean qubits or every prefactor in $`b=\Theta(N)`$.
Restrictions proved for particular source-processing interfaces do not
settle the unrestricted frame problem. The grouped improvement is a
T-count theorem; it does not establish an optimal T-depth tradeoff.

## What the current research resolves

The geometric operator source itself admits a sharp exact synthesis statement:
on $`m\ge2`$ dirty qubits, its minimum exact T-count is $`2m-4`$;
its controlled version requires exactly $`2m-2`$. The
[source construction and Pauli-transfer argument](OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit)
allow arbitrary returned clean and dirty helpers. Thus better exact synthesis of
that same controlled source cannot make an individual call sublinear in the
precision. This is a primitive bound. Costs of separate calls cannot be added
to infer a lower bound for an unrestricted frame compiler.

Moving amplitude amplification to the end also needs a new construction.
Writing $`P=JJ^\dagger`$ and $`B_i=J^\dagger Q_iJ`$, reuse of the same flags gives

```math
J^\dagger Q_2Q_1J
=B_2B_1+J^\dagger Q_2(I-P)Q_1J.
```

The second term is a coherent return from the rejected space. It need not be
small even when the individual accepted blocks are exact. For example,
$`Q_1=Q_2=H\otimes H`$ on the two flags has $`B_1=B_2=1/2`$, but the
product's accepted block is one, not $`1/4`$. The existing layers therefore
cannot be concatenated as unamplified blocks and treated as a product of their
accepted actions. Extra initialized history would also change the two-clean
budget. This observation rules out that inference, not every global design.

## Limits of two source-reuse shortcuts

The [source-reuse limits](SOURCE_REUSE_LIMITS.md) give two more precise
restrictions on natural proposals for removing the repeated precision cost:

- A nilpotent carried contraction cannot act as an exponentially accurate
  scalar on an encoding of every arbitrary dirty input using only two
  initialized qubits. The required initialized width is at least
  $`\log_2 L-O(1)`$ for that interface, regardless of dirty width.
- Pulling the common operator-source basis change outside the stream
  transforms the programmed masks. An allowed one-bit transformed mask
  already has linear exact and fine-accuracy T cost.

Both statements concern specified intermediate interfaces. They neither
make separate source costs additive nor strengthen the unrestricted
$`\Omega(N)`$ full-frame lower bound. A jointly synthesized global
block can avoid those interfaces.

## A sufficient construction to seek

At the stated endpoint, it would suffice to construct one actual unitary $`Q`$
using two initialized flags and at most $`N+n+7`$ dirty qubits such that

```math
\left\|2J_2^\dagger QJ_2-(W\otimes I_b)\right\|\le\eta/4,
\qquad T(Q)=O(N),\qquad G(Q)=O(N^2).
```

The source, all table queries, and every helper must be included in these
budgets. The accepted action must hold on every logical and dirty input,
including reference correlations; no return condition is assumed on rejected
branches. The normalization-two amplification lemma then gives the desired
complete-isometry compiler with three charged calls, using the same two flags.
No such jointly charged block is currently established. A proof of this
sufficient interface, or a different full-frame construction, would close the
upper-bound side of the endpoint without requiring intermediate source return.

## What still costs more than linear

The current grouped proof pays for three distinct operations. For a group of
s levels above a suffix of length r, put $`e=n-r`$. Its private initialized
work, coefficient-table rows, and precision-source width satisfy

```math
w_g=O(\log(s+2)),\qquad
Q_g=\Theta(s2^{n-r}),\qquad
m_g=L+\lfloor r/4\rfloor+8.
```

The table count is for the current padded star representation, not a lower
bound on representations of the same operator. Conditional suffix work
supplies the first resource. Exponentially growing groups make
$`\sum_g Q_g=O(N)`$ while leaving
$`O(\ell_*(n))`$ precision-sized source uses. This is how the current
construction obtains its bound.

Removing the depth-label register alone would not give a linear compiler.
If a single group had $`s=n-r_0`$ with fixed suffix length $`r_0`$, its
current table would contain $`\Theta(nN)`$ rows. At $`L=N`$ this
m-bit-word table has $`\Theta(nN^2)`$ bit positions, and the existing
$`O(Q_gm_g)`$ lookup accounting becomes $`O(nN^2)`$. It therefore no
longer certifies the desired $`O(N^2)`$ Clifford bound; the table size
alone is not a gate lower bound. The explicit endpoint dirty allocation does
not meet the sufficient threshold for the proved banked refinement, so that
refinement does not apply directly. The current streamed coarse-program
bound also becomes $`O(s2^e)=O(nN)`$ for that giant group under its uniform
star-normalization accuracy condition. The weighted representation below
can instead reuse a cheap coarse frame from the borrowed-workspace theorem.
A successful replacement must
control **both** the repeated precision cost and the expanded table cost;
it must also charge coarse symbol queries and their interpreter. Constant
private width alone does not settle either resource theorem.
The residual entries are functions of the original tree angles. The
[tree-generator factorization](SOURCE_REUSE_LIMITS.md#3-tree-generators-compress-the-residual-classically)
now represents them using $`O(N)`$ classical coefficients and path products.
The downward factors still use target-frame transport. Applying that
transport by calling the target frame would be circular; replacing it by
coarse transport leaves mixed errors. Classical compression therefore does
not yet supply a cheaper coherent query implementation.

The following distinctions keep the next search from repeating shortcuts:

| Route | What has been established | What remains open |
|---|---|---|
| Better synthesis of the same source | One exact controlled source already has linear cost in its width | Joint synthesis across uses, or a different source |
| Move the source basis change outside all programming | Individually transformed masks can also have linear T cost | A jointly synthesized source/program word |
| Reuse flags and amplify only once | Rejected branches can return coherently to the accepted subspace | A designed global block with a proved rejection-space action |
| Replace an initialized nilpotent source by dirty encoding | The specified full-output relation needs growing initialized dimension | Nonnilpotent kernels, conditional sectors, or another block representation |
| Transfer the old precision-state compiler into a conditional sector | That earlier construction had a source-fit limitation | The present operator-source construction already bypasses that limitation; conditional suffix work is not ruled out |

The first four statements have their analytic homes above and in
[source-reuse limits](SOURCE_REUSE_LIMITS.md). The last is an archived
construction limitation, not a general impossibility result. None supplies
an additive lower bound for the full frame.

## Revision decision and next bounded pass

The revision retains the theorem frontier above and narrows the next task.
The parallel lookup work has produced a new same-circuit count/depth
guarantee. The tree-transport work has produced operator identities and norm
control, but no improved endpoint compiler. These are different levels of
completion.

| Ingredient | What is available | What it does not supply |
|---|---|---|
| Parallel dirty lookup | Exact returned indicators and the simultaneous theorem in [the lookup chapter](PARALLEL_DIRTY_LOOKUP.md) | Optimal T-depth, interpolation below its sufficient width, or the selected endpoint allocation |
| Residual data | $`O(N)`$ scalar tree generators and exact path products | A charged coherent evaluator for those products |
| Tree transport | Exact Gram matrix and a complete logical unitary for normalized supported columns, without a depth register | An elementary Clifford+T cost improvement |
| Forward weighted residual | A complete one-signal-flag block with constant normalization and batched native cost $`T=O(L\sqrt N)`$, $`G=O(NL)`$ under its sufficient dirty width | Linear endpoint T-count; this implementation costs $`O(N^{3/2})`$ at $`L=N`$ |
| Coarse circuit | The retained borrowed compiler gives $`T(C)=O(\sqrt N)`$, $`G(C)=O(N)`$ at fixed accuracy and the endpoint dirty allocation; both clean flags can remain untouched | The history unitary, weighted transport, or an arbitrary controlled version of C |
| High-precision residual circuit | The proved grouped compiler and its complete-input error budget | A constant total precision charge for one global weighted block |

The [weighted norm proof](ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
gives $`\|D_h\mathcal P_W\|,\|\mathcal P_C^\dagger D_k\|\le2\sqrt\gamma`$,
where $`\gamma=\max_v(1-|g_v|^2)`$. This removes an intrinsic height
penalty for those operators. The available history unitary, however,
realizes column-normalized transport. Undoing its normalization by separate
input and output filters still carries a possible $`\sqrt n`$ factor.
The correlations that make the weighted norm small must be used by an actual
circuit. Also, gamma alone does not control the literal phases in the
diagonal term D.

The scoped exclusions remain useful. Fixed-order coarse transport misses
mixed terms at exponential precision; a right-spine entry proves this
without a simulation. The two-bit Pauli routing of rejected components
removes pair returns but leaves a cubic return. Source minima, transformed
masks, nilpotent encodings, and these witnesses constrain their specified
interfaces only. None strengthens the unrestricted frame lower bound.

### Reuse: the coarse-frame prerequisite is already available

The revision identifies a useful existing ingredient. At fixed coarse
accuracy $`0\lt\varepsilon_0\le1/64`$, the borrowed-workspace theorem
already gives $`T(C)=O(\sqrt N)`$ and $`G(C)=O(N)`$ using the endpoint's
dirty pool. Its depth-dependent error schedule, native local words, exact
inactive action, and exact helper return are proved in that chapter.
No initialized-work assumption is required, so both clean flags can remain
untouched, even when occupied by an outer computation.

The [reuse argument](ENDPOINT_TREE_TRANSPORT.md#the-retained-borrowed-compiler-already-supplies-a-cheap-coarse-frame)
checks the needed tree structure and uniform subtree accuracy. It gives
$`\|D\|\le\varepsilon_0`$ and bounds each weighted term by
$`2\varepsilon_0`$. Thus a separate coarse-program bottleneck does not
remain for this representation. The older $`O(nN)`$ estimate belongs to
the giant uniform-star construction, whose coefficient normalization
requires much finer coarse words. Its accuracy condition cannot simply be
replaced inside that proof.

Do not spend the next pass rebuilding this coarse circuit. Its unconditional
application and actual inverse are priced already. The unresolved gate
cost is that of the weighted transport itself, including any controlled
versions needed by a later composition.

### Completed: a forward block with batched native synthesis

The [weighted-block construction](WEIGHTED_TRANSPORT_BLOCK.md) completes the
operator and workspace parts of the proposed forward-component audit. For
$`F=\iota D_h\mathcal P_W`$ and fixed $`\alpha=4\varepsilon_0`$, it
constructs an actual one-signal-flag unitary with accepted block
$`F/\alpha`$. The unused second clean flag supplies the same two-flag
interface as the proposed target. Rejected inputs have an explicit unitary
action; no intermediate flag is reset or assumed zero.

The first stopping proposal failed because its weighted marker columns
were not orthogonal. A recursive Schur-complement correction adds the
required rejection amplitudes. The network contains $`O(N)`$ local
four-mode and two-mode factors on $`2N`$ basis modes, so one signal flag
suffices. Certified algebraic preprocessing handles zero defects without
requiring exact zero tests on arbitrary supplied angle evaluators.

Packing disjoint factors by tree level and using the retained borrowed
reflection interpreter gives, for $`b\ge n+\lceil\sqrt N\rceil+7`$,

```math
\left\|J_2^\dagger Q_hJ_2-
\frac{F\otimes I_b}{4\varepsilon_0}\right\|\le\eta/128,
\qquad T=O(L\sqrt N),\qquad G=O(NL).
```

Each level uses at most 18 fixed addressed rotation templates. Its bank
count is proportional to the square root of its number of nodes, and the
precision allowance varies geometrically with depth. Exact packing,
gathering, and the prescribed physical marker reindex are all charged;
before absorbing polynomial terms the counts are
$`T=O(L\sqrt N+n^4)`$ and $`G=O(NL+n^4)`$.

At $`L=N`$, this improves the forward component from $`O(N^2)`$ to
$`O(N^{3/2})`$ T gates, with $`O(N^2)`$ Clifford gates. The selected
endpoint allocation $`b=N+n+7`$ satisfies the sufficient width. The signal
flag is arbitrary data throughout native synthesis; the other clean flag
is untouched, and all borrowed helpers return exactly. The component meets
the accepted-block error, workspace, and Clifford requirements, but still
does not establish linear T-count. Neither implementation upper bound is
a lower bound for this component or the frame.

### Remaining: native precision sharing with an occupied signal flag

The remaining loss is the product of word precision and square-root table
size in the borrowed interpreter. Another norm identity would not remove
that product. The two-clean stage compiler is not a direct substitute:
the dilation already occupies one signal flag, leaving only one fresh
initialized bit. The [direct flag-merging word](SOURCE_REUSE_LIMITS.md)
has an accepted-block error of at least $`1/2`$ for exact sine and cosine.
That excludes this substitution, not every one-clean construction.

The cheap coarse frame remains available and needs no reconstruction.
The next attempt should synthesize the unnormalized weighted operator
jointly. It may choose a different unitary completion: the required
interface is the accepted-block estimate, an actual unitary word and its
literal inverse, and a charged workspace ledger. Closeness to the selected
Gram–Schmidt completion is unnecessary. Indeed, a
[one-level witness](ENDPOINT_TREE_TRANSPORT.md) shows that this completion
can stay distance two from its zero-defect value even while its accepted
block tends to zero. No result permits adding costs as a lower bound,
treating a dirty register as initialized, or invoking a target-frame oracle.

Even a linear-cost forward block would be a component result. Combining
it with the reverse weighted term, the diagonal, coarse C, and amplification
must still meet the [global-block contract](#a-sufficient-construction-to-seek).
Independent accepted blocks cannot be multiplied while ignoring returns
from rejected spaces; all required controls and flag lifetimes remain charged.

Optimal T-depth and stronger unrestricted lower bounds remain separate
projects. A new lower-bound attempt needs an invariant preserved under
arbitrary Clifford interlayers and returned dirty helpers. Practical
constants and end-to-end emission also remain separate from this endpoint
decision. No improved full-frame resource theorem is claimed by this component audit.

## Evidence and remaining implementation work

The current proofs retain the complete-input guarantee; revision of their
source conditioning, reverse-word ordering, workspace ledgers, and error
sums has not identified a defect. This is an internal audit, not independent
peer review. The [finite checks](VERIFICATION.md) cover native source words,
star-support partitions, literal phases, actual inverses, rejected-space
composition, and small complete dirty-input blocks.

The strongest grouped construction is not emitted end to end as an
elementary Clifford+T circuit. Its coefficient-table fixture represents the
table action directly; its integer grouping tests use illustrative constants.
The proof's fixed workspace constants and crossover threshold remain
existential. Explicit selected-atom pseudocode and a register-lifetime table
would improve auditability before practical resource estimates or a full
emitter are attempted. Small tests do not establish the asymptotic theorem,
its optimality, or its literature priority.
