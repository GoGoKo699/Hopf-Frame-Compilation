# Research status and the constant-clean endpoint

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [One-clean compiler](ONE_CLEAN_COMPILER.md) · [Grouped refinement](CONDITIONAL_SUFFIX_COMPILER.md)

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
| Zero-clean layerwise real-frame T-count | $`T=O(N+nL)`$, $`G=O(NL)`$, at $`a=0`$, $`b\ge L+n+7`$ | Full-operator approximation using a borrowed amplification signal; [signal-symmetry corollary](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) |
| Grouped one-clean real-frame T-count | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge L+n+7`$ | Precision-uniform grouped bound without additional word banks; [one-clean extension](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean T-count with additional dirty banks | $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge2(L+n+7)`$ | Matches the lower bound if $`L\ell_*(n)^2\le N`$ or $`b\le N/\ell_*(n)`$; these are sufficient regimes; [banked one-clean proof](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean phase-dressed complex magnitude frame | The same grouped and banked T-counts, at $`b\ge L+n+8`$ and $`b\ge2(L+n+8)`$, respectively | Compose the real compiler and literal phase diagonal in the same workspace; [composition corollary](ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames) |
| T-depth with additional dirty banks | $`D_T=O(NL/b+\min\{nL+n^3,L\ell_*(n)+n^4\})`$ at $`a=2`$, $`b\ge2(L+n+7)`$, with $`T,G=O(NL)`$ | Choose between the layerwise and grouped [schedules](T_DEPTH_COMPILER.md); real frames; optimizing depth may increase T-count; no matching frontier established |
| Simultaneous T-count and T-depth | $`T=O(\sqrt{NL}+L\ell_*(n))`$, $`D_T=O(\min\{nL+n^3,L\ell_*(n)+n^4\})`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge C(L+n+7+\sqrt{NL})`$ | Same real-frame circuit, for sufficiently large fixed C; [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md); T-depth optimality remains open |

The retained frontier is the best applicable bound among these constructions.
For example, at fixed L, $`a=2`$, and $`b=L+n+7=\Theta(n)`$,
the arbitrary-budget theorem and its matching splice give
$`T^\star=\Theta(N/n)`$, sharper than the grouped $`O(N)`$ estimate.
This does not change the high-precision endpoint below. The one-clean count
constructions also work with additional unused clean qubits; the banked
matching comparison above is for fixed $`a=1`$. The depth schedules retain
their separately proved two-clean allocations.

Every retained frame construction preserves the prescribed completion used
by Hopf QBP. Finite-precision substitution has the fixed-parameter bias and
dirty-reference guarantees of the [QBP approximation theorem](QBP_APPROXIMATION.md).
It is not a claim about differentiating the discrete synthesis algorithm.
The phase-dressed complex magnitude frame now inherits the grouped count
bound by full-isometry composition at its separate dirty threshold. Leaf-phase
derivatives still use their own measurement stream; they are not additional
columns of the magnitude frame. The older two-clean complex theorem remains
a valid baseline.

The layerwise operator-source compiler remains a dependency and a useful
fallback. Its literal-diagonal and complete one-target U(2) multiplexor
corollaries retain matching banked frontiers as independent capabilities. Earlier, weaker
endpoint bounds are superseded as frontiers; exploratory routes are
recoverable from the [archived research snapshot](../provenance/README.md#earlier-research-snapshot).

## One clean qubit now suffices for the grouped count bound

The [new primitive](ONE_CLEAN_COMPILER.md) conjugates one programmed
scalar-source block by another, with different Pauli routing on the
logical target. Its accepted block contains two signed overlaps; their
unwanted dirty terms cancel by anticommutation. Two independently
programmable geometric streams fit on paired Majorana generators in one
precision-sized dirty core. A five-call amplification word returns the
clean flag and core within the complete-input error bound.

For a tree layer, core, selectors, and one arbitrary helper occupy exactly
$`L+n+7`$ dirty wires. In the grouped compiler, the sole external clean
qubit stores the active-suffix predicate. Moving the old scalar flag into
the known-zero active suffix changes its private reservation by one bit;
the new primitive compiles the fixed deepest tail. This establishes the
one-clean grouped and banked bounds above. It is a constructive reduction
of initialized workspace, not a new lower bound or a T-depth theorem.

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

The [one-clean layer compiler](ONE_CLEAN_COMPILER.md) has
$`T=O(N+nL)`$ with $`b\ge L+n+7`$, extending the original
[two-flag source construction](OPERATOR_SOURCE_COMPILER.md). It pays the precision cost at each
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

The same grouped upper bound now holds at $`a=1`$ with the identical
dirty allocation. This reduces the sufficient clean count without closing
the selected $`a=2`$ endpoint or proving anything impossible at $`a=0`$.

With a sufficiently large $`a=\Theta(n)`$ clean reservation and
$`b=\Theta(N)`$, the shared-source compiler instead attains
$`T^\star=\Theta(N)`$ at $`L=N`$. This sufficient clean reservation is
not proved necessary.

The open question is whether two clean qubits permit a jointly charged
$`O(N)`$ construction at the explicit allocation above, or whether a stronger
general lower bound holds. The grouped endpoint upper bound does not yet extend to zero clean qubits
or to every prefactor in $`b=\Theta(N)`$. The zero-clean layerwise
corollary gives $`O(N\log N)`$ at the stated dirty allocation.
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
| Reuse flags and amplify only once | Rejected branches can return coherently; the new [affine assembly](RESIDUAL_ASSEMBLY.md) explicitly handles them | Lowering the assembled block's repeated precision cost |
| Replace an initialized nilpotent source by dirty encoding | The specified full-output relation needs growing initialized dimension | Nonnilpotent kernels, conditional sectors, or another block representation |
| Transfer the old precision-state compiler into a conditional sector | That earlier construction had a source-fit limitation | The present operator-source construction already bypasses that limitation; conditional suffix work is not ruled out |

The first four statements have their analytic homes above and in
[source-reuse limits](SOURCE_REUSE_LIMITS.md). The last is an archived
construction limitation, not a general impossibility result. None supplies
an additive lower bound for the full frame.

## Revision decision and next bounded pass

The revision distinguishes complete-frame guarantees from ingredients of a
possible endpoint construction. The exact CNOT/depth theorem and the
sufficient-clean matching T-count theorem are complete. One-clean grouped
compilation now covers the real frame and, by sequential composition, the
phase-dressed complex magnitude frame at their respective dirty thresholds.
The new real-rotation signal symmetry also gives the zero-clean layerwise
bound above. None removes the grouped endpoint's iterated logarithm.

### Reusable ingredients and their remaining obligations

| Ingredient | Established interface | Remaining obligation for a global endpoint block |
|---|---|---|
| Coarse frame C | The borrowed compiler gives $`T(C)=O(\sqrt N)`$, $`G(C)=O(N)`$ at fixed accuracy in the endpoint dirty pool, with exact helper return and no initialized work | Its actual inverse is available; any new controlled version must still be implemented and charged |
| Compact residual data | $`O(N)`$ local generators and exact path products represent $`C^\dagger W-I`$ | Classical compression alone is not a coherent evaluator |
| Weighted forward block F | A complete one-signal-flag dilation with fixed normalization; $`T=O(N+nL)`$, $`G=O(NL)`$, $`b\ge L+n+7`$; native synthesis uses no additional clean flag | This is one term, and still pays precision at every depth |
| Exact-dirty-return alternative | The same component has $`T=O(L\sqrt N)`$, $`G=O(NL)`$ at $`b\ge n+\lceil\sqrt N\rceil+7`$ | Useful when exact work return or its different width matters; it is not the best endpoint component bound |
| Affine forward and reverse blocks | The [assembly proof](RESIDUAL_ASSEMBLY.md) combines the diagonal with forward transport into one signal block, and constructs the reverse block by a swapped-pair actual inverse | Their current native T-count still contains $`nL`$ |
| Two-flag global assembly | A selector and a shared signal combine those two blocks with exact normalization two; native controls, inverses, and workspace are charged | This closes the composition step, not the linear-precision synthesis step |
| One-clean grouped frame compiler | Complete-frame error including approximate core and suffix-work return | It already gives the best retained endpoint upper bound; its precision charge is not a lower bound |

The coarse frame needs no reconstruction. The [weighted norm proof](ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
uses its uniform subtree accuracy to bound both weighted terms independently
of height. The [forward-block proof](WEIGHTED_TRANSPORT_BLOCK.md) supplies
orthogonal rejection amplitudes, a complete unitary on both signal sectors,
and certified preprocessing at zero defects. Its physical marker permutations
and lookup work are charged.

The real-rotation word's exact X symmetry is particularly useful here:
its initialized-signal error lifts to full-operator error, so native synthesis
can borrow that signal from the dirty pool. Only the forward block's own
logical dilation signal needs initialization when its accepted block is used.
The scalar-phase word does not have the same symmetry. The grouped compiler
also still needs its one external clean predicate; the zero-clean layerwise
corollary does not remove that requirement.

At $`L=N`$, the current forward component costs $`O(N\log N)`$ T gates,
whereas the complete frame already costs $`O(N\ell_*(n))`$. Thus an
improvement to that component is useful only as part of a cheaper complete
assembly. The assembly interface is now established. The remaining construction
task is a smaller joint precision charge for its two compatible terms.

### Completed: two flags assemble the whole residual

Use the exact identity

```math
C^\dagger W=A+F+R,\qquad A=I+D,\qquad
F=\iota D_h\mathcal P_W,\qquad
R=\mathcal P_C^\dagger D_k\iota^\dagger.
```

A separate selector branch for each of A, F, and R would require a third
initialized selector/signal role. The [new construction](RESIDUAL_ASSEMBLY.md)
first incorporates the diagonal into the forward tree:

```math
S=A+F,\qquad \alpha=4\varepsilon_0,\qquad s=2-\alpha,
\qquad \|S\|\le1+2\varepsilon_0\lt s.
```

The scalar Schur recursion now includes the marker's diagonal contribution
in its stop output. Root mixing and the terminal rejection phase are adjusted
so that one logical signal still supplies a complete unitary dilation.
This gives an actual block $`S/s`$. A swapped-pair construction followed by
its actual inverse gives $`R/\alpha`$.

Prepare the selector with probabilities $`s/2`$ and $`\alpha/2`$,
select between these two one-signal blocks, and undo that preparation.
Because $`s+\alpha=2`$, its two-flag accepted block is

```math
\frac{s}{2}\frac{S}{s}
+\frac{\alpha}{2}\frac{R}{\alpha}
=\frac{C^\dagger W}{2}.
```

The implementation uses one common algebraic target approximation for both
terms. Controlled templates are exactly inactive on the other selector
sector, including all dirty inputs. Final gathering is explicitly controlled;
no T gate is silently promoted to a controlled T. The proof gives the
complete flag-lifetime table and includes coefficient preparation, physical
routing, actual inverses, and normalization-two amplification. Appending the
already priced coarse C yields the prescribed frame.

At $`a=2`$ and $`b\ge L+n+7`$, the resulting complete circuit has
$`T=O(N+nL)`$ and $`G=O(NL)`$. This is an explicit assembly route;
it is weaker than the retained grouped endpoint upper bound. In particular,
the proof does not turn an opaque existing F block into an S block with a
constant number of calls: it recompiles the affine tree's local data.

### Next task: share precision in the affine and reverse blocks

The assembly now supplies a charged reduction. For compatible actual
component circuits with the stated accepted blocks and their implemented
selector controls, its T-count is bounded by

```math
O(T_S+T_R+N+L),
```

and its Clifford count by $`O(G_S+G_R+NL)`$. At the selected endpoint,
$`L=N`$ and $`b=N+n+7`$, it therefore suffices to achieve

```math
T_S+T_R=O(N),\qquad G_S+G_R=O(N^2)
```

within that same dirty pool. A precision-uniform sufficient target is
$`T_S+T_R=O(N+L)`$, $`G_S+G_R=O(NL)`$. These improved component
costs remain unproved. The current templates pay $`O(L)`$ precision work
at each of n depths.

The next bounded pass should synthesize the **affine** transport
$`S=A+F`$ together with the reverse term's compatible interface. A cheaper
F-only routine does not automatically price S. The local coefficients and
scalar storage recursion contain only $`O(N)`$ classical data; their
coherent implementation, masks, routing, and precision work must all be
charged. Replacing the repeated native source calls by a shared word is a
construction task, not a consequence of that classical data count.

The assembly does not require the future blocks to approximate the particular
Gram–Schmidt completions used in the current proof. They may choose other
unitary completions, provided their accepted-block estimates hold on all
logical and dirty inputs, their controlled actions and literal inverses are
implemented, and the two-flag allocation is respected. Amplification then
controls the complete output, including rejected-space returns.

A successful pass must give an actual native circuit and summed resource
ledger improving the repeated precision term. A new norm identity, a free
coherent evaluator, or a bound for the forward accepted block alone would
leave that task incomplete. Any obstruction must remain scoped to the
specific source or assembly interface it analyzes.

### Boundaries to carry into that pass

- Rejected components can return coherently; the exact counterexamples in
  [source-reuse limits](SOURCE_REUSE_LIMITS.md) still apply to their specified
  words. The successful conjugated source is a different construction.
- Fixed-order coarse transport misses required mixed terms at exponential
  accuracy. The cheap coarse frame does not supply free target transport.
- The chosen forward dilation can remain far from its zero-defect completion
  even when its accepted block is small. A new assembly may choose another
  completion, but must specify its unitary action and actual inverse.
- Source-call minima, transformed-mask costs, and failed merges do not make
  separate costs additive as a lower bound for unrestricted compilation.

The depth tradeoff, its lower-bound gap, and explicit practical grouping
constants remain separate tasks. They need not be reopened to test this
joint precision construction. The publication's established compiler results remain valid
independently of whether this endpoint route succeeds.

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
