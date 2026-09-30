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
| Two-clean real-frame T-count | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$, at $`b\ge L+n+7`$ | Current strongest unbanked bound; [conditional-suffix theorem](CONDITIONAL_SUFFIX_COMPILER.md) |
| Two-clean T-count with additional dirty banks | $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$, at $`b\ge2(L+n+7)`$ | Matches the lower bound if $`L\ell_*(n)^2\le N`$ or $`b\le N/\ell_*(n)`$; these are sufficient regimes; [banked proof](CONDITIONAL_SUFFIX_COMPILER.md#8-additional-dirty-banks-give-a-width-tradeoff) |
| T-depth with additional dirty banks | $`D_T=O(NL/b+\min\{nL+n^3,L\ell_*(n)+n^4\})`$ at $`a=2`$, $`b\ge2(L+n+7)`$, with $`T,G=O(NL)`$ | Choose between the layerwise and grouped [schedules](T_DEPTH_COMPILER.md); real frames; optimizing depth may increase T-count; no matching frontier established |
| Simultaneous T-count and T-depth | $`T=O(\sqrt{NL}+L\ell_*(n))`$, $`D_T=O(\min\{nL+n^3,L\ell_*(n)+n^4\})`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge C(L+n+7+\sqrt{NL})`$ | Same real-frame circuit, for sufficiently large fixed C; [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md); T-depth optimality remains open |

Every retained frame construction preserves the prescribed completion used
by Hopf QBP. Finite-precision substitution has the fixed-parameter bias and
dirty-reference guarantees of the [QBP approximation theorem](QBP_APPROXIMATION.md).
It is not a claim about differentiating the discrete synthesis algorithm.
The strongest grouped theorem is stated for real frames; the separately
proved two-clean complex extension retains its baseline bound.

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
bound also becomes $`O(s2^e)=O(nN)`$ for that giant group. Thus a successful replacement must
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

## Next research checkpoint

### Completed: parallel dirty lookup at a fixed T-count budget

The [parallel-loader proof](PARALLEL_DIRTY_LOOKUP.md) passes the proposed
audit. Nested exact bilinear echoes construct an H-output dirty indicator
using at most $`4H`$ total dirty indicator work, $`O(H)`$ T and Clifford
gates, and $`O((1+\log H)^2)`$ depth. It uses literal exact Toffolis
and no initialized selectors. The completed whole-word query preserves
arbitrary dirty/reference inputs and is identity on inactive zero rows.

At bank count mu, $`H=Q/\mu`$, the query requires
$`\mu m+4H`$ extra dirty bits and retains
$`T=O(H+\mu m)`$. Choosing the count-efficient banks gives the new
simultaneous full-frame theorem above. The local lifetime ledger charges
all borrowed work, actual inverses, and linear table maps; those Clifford
maps still have nonzero depth. The general parallel-lookup idea is
attributed to Low, Kliuchnikov, and Schaeffer, Appendix C.

The first revision target is therefore achieved, already at sufficiently
large square-root dirty width for fixed accuracy. A further depth project
could seek a count-preserving interpolation below that width, but it would
need a separately charged partial-indicator construction. The source's
high-precision term is unchanged; this result does not close the selected
$`O(N)`$ T-count endpoint.

### Next: one bounded joint-block attempt at the linear endpoint

The closest endpoint question remains whether the residual construction
admits a **joint source and table implementation** that avoids its remaining
group multiplicity. Before attempting another full-frame theorem, specify
one candidate interface and account for:

1. The complete accepted operator and the coherent action on rejected work,
   including any intermediate reuse of flags or addresses.
2. Live initialized and dirty registers throughout the actual circuit, with
   no logical input treated as freely initialized outside its active sector.
3. Total charged source uses, table rows and word widths, coarse symbol queries
   and their interpreter, routing, reflections, and actual inverses. At the
   endpoint the targets are $`O(N)`$ T gates and
   $`O(N^2)`$ Clifford gates within the stated allocation.
4. Complete-input error after composition, including every approximate work
   return and reference correlation.

A proposal that only removes a clean label, cancels a formal basis change,
or multiplies accepted blocks does not yet meet this checkpoint. A failed
candidate should produce a restriction with its exact interface stated,
while preserving the current compiler as a fallback. A general stronger
frame lower bound remains a separate route and requires a new invariant.

The two-to-three-level diagnostic has resolved two limited questions.
Tree-node data compresses the residual classically, but substituting coarse
transport leaves nonzero mixed terms. Separately, routing scalar-source
rejections into the three nonzero two-bit Pauli labels removes pair returns
but leaves a cubic coherent return. The [scoped proofs](SOURCE_REUSE_LIMITS.md)
and finite tests do not rule out non-Pauli histories or another global block.
That routing proposal also retains separate source calls and is not a
half-unitary compiler by itself.

The next endpoint target is a charged implementation of the tree path
factors, or a full-operator invariant that avoids implementing them
separately. A proposed replacement must retain a constant total precision
charge and control the mixed terms to exponential accuracy; polynomially
small coarse errors alone do not suffice.

A bounded attempt should specify both a three-level full-unitary block and
its extension to arbitrary height. It must account for coefficient and
coarse-program queries together. Expanding all ancestor pairs, invoking
the target frame to generate its own path amplitudes, or retaining a fresh
precision charge per depth does not meet the endpoint target. Park that
candidate if its cost or mixed-error recurrence fails; keep the current
compiler as the fallback.

### Later: stronger lower bounds and practical constants

A matching depth theorem needs a new lower-bound invariant or a further
upper construction; an improved polynomial allowance alone cannot settle
the width gap above. A lower-bound attempt should start by naming the
property that arbitrary Clifford interlayers and returned dirty helpers
cannot erase. Explicit constants and complete circuit emission are useful
implementation work, but remain separate from this conceptual gap and the
linear T-count endpoint.

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
