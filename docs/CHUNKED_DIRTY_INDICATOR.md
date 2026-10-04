# A chunked dirty indicator and the late-layer query budget

[Dirty conjunctions](DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences) · [Bilinear queries](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction) · [Current frontier](OPEN_PROBLEM.md)

An exact dirty indicator can trade a polynomial in the address chunk
length for a shorter sequence of address-dependent stages. The tree
below interpolates between linear address T-depth with linear-size work
and logarithmic address T-depth with polynomial overhead. It uses no
initialized work, and all shared address and tree-parent wires are
read-only controls during a completed stage.

For an r-bit address, an arbitrary $`2^r`$-bit output Y, and any integer
$`1\le\ell\le r`$, the circuit implements

```math
|x\rangle|Y\rangle|h\rangle
\longmapsto |x\rangle|Y\oplus e_x\rangle|h\rangle
```

exactly on every input, with

```math
T,G,w=O\!\left(2^r(\ell+2)^3\right),\qquad
D_T=O\!\left(\left\lceil\frac r\ell\right\rceil
\log_2(\ell+2)\right).
```

Here w counts additional dirty work; including Y changes only the
constant. For r equal to zero, apply X to its single output and use no
work. All logarithms below are binary. The bound concerns T-depth;
shared-control Clifford gates may be serialized and their count is
included in G.

## 1. A reversible tree on arbitrary dirty inputs

Partition the address, in its fixed bit order, into t consecutive
nonempty chunks of lengths $`b_i\le\ell`$, using chunks of length ell
except possibly the last. Put $`s_0=0`$ and
$`s_i=b_1+\cdots+b_i`$, so $`s_t=r`$ and
$`t=\lceil r/\ell\rceil`$.

Reserve a dirty root $`w_\emptyset`$ and one distinct dirty wire
$`w_p`$ for each binary prefix p of length $`s_i`$, at every boundary
i. The final boundary has one leaf for every complete address. These
tree wires are disjoint from the original address and Y.

Define a chronological top-down word F. At stage i, for every prefix
p of length $`s_{i-1}`$ and every next-chunk value v, perform

```math
w_{pv}\longmapsto w_{pv}\oplus
w_p[\,x^{(i)}=v\,].
```

Each update is a conjunction of $`b_i+1`$ literals into its child.
Use the [exact read-only conjunction](DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
with a private dirty helper pool for that edge. For k literal controls,
its count and width are $`O((k+1)^\alpha)`$, where
$`\alpha=\log_2 3`$, its Clifford count is
$`O((k+1)^\alpha\log(k+2)\log\log(k+4))`$, and its T-depth is
$`O(\log(k+2))`$. It restores all its helpers exactly. Its original
controls occur only as Clifford CNOT controls; all non-Clifford gates
act on private work.

At one tree stage, parents belong to the preceding boundary and targets
to the next boundary. No parent is a target in that stage. Child targets
and private helper pools are distinct, so the conjunctions' T layers
can overlap even when they share a parent or address wire. Their Clifford
interactions need not have constant depth. Finish the entire stage before
starting the next, and reuse its returned helper pool afterward.

For a fixed address x, the completed F is an invertible linear map on
the tree bits. Write it as $`w\mapsto A_xw`$. Its action on the root
basis vector is particularly simple:

```math
A_x e_\emptyset
=e_\emptyset+\sum_{i=1}^{t}e_{x_{1:s_i}}.
```

Indeed, a unit root propagates at each stage to exactly the child whose
chunk agrees with x; every other child receives zero. This statement
describes the difference of two arbitrary dirty inputs, and assumes no
tree wire was initialized.

## 2. Extracting the selected leaf and restoring the tree

Let X at the root denote its Pauli bit flip. Define the algebraic
conjugate, with operators on the right acting first,

```math
K=F X_\emptyset F^\dagger.
```

Its chronological native word is therefore
$`F^\dagger,X_\emptyset,F`$. Reversing this order generally fails
to propagate the root flip to the deepest leaves. On any tree input,

```math
K:\ w\longmapsto w+A_xe_\emptyset.
```

Thus K flips the selected boundary path, including its leaf, independently
of every initial dirty bit. It is an exact involution, though its actual
reversed native word is used whenever an inverse is requested.

Let C be the Clifford word consisting of one CNOT from each leaf
$`w_p`$ to its corresponding $`Y_p`$. Execute chronologically

```math
C,\quad K,\quad C,\quad K^\dagger.
```

The first C adds the original leaf word to Y. The second adds that same
word plus $`e_x`$. Their arbitrary offsets cancel. The final inverse
of K restores the entire tree, including its root and intermediate
boundary wires. Each edge call already returned its own helpers. Hence
the complete operation is exactly the stated indicator tensored with
identity on all work, including when the work is entangled with a
reference.

For example, with two one-bit chunks, the first root conjugation flips
the root, the selected one-bit prefix, and the selected two-bit leaf.
The two leaf readouts retain only the latter flip. Omitting the final
inverse would leave the other two flips behind.

## 3. Live width and gate ledger

The boundary prefix lengths increase by at least one, so

```math
1+\sum_{i=1}^{t}2^{s_i}\le 2^{r+1}-1.
```

In particular, the number of tree edges across all stages is less than
$`2^{r+1}`$. Counting each edge once, rather than multiplying the
largest stage by the number of stages, gives

```math
T(F)=O\!\left(2^r(\ell+2)^\alpha\right),\qquad
G(F)=O\!\left(2^r(\ell+2)^\alpha
\log(\ell+3)\log\log(\ell+5)\right).
```

All tree wires remain live. Private conjunction helpers are reserved for
the largest stage and reused only after exact return. Therefore their
maximum simultaneous width, together with the tree, is
$`O(2^r(\ell+2)^\alpha)`$. No factor t occurs in this width.
There are t sequential stages, each with T-depth
$`O(\log(\ell+2))`$.

The output echo calls F or its actual inverse four times and adds
$`2^{r+1}`$ leaf CNOTs. These are fixed factors in the preceding
ledger. The safe common majorant $`P(\ell)=(\ell+2)^3`$ bounds
T-count, Clifford count, and width, proving the introductory theorem.
The endpoints ell equal to one and ell equal to r give, respectively,
$`O(2^r)`$ count/work with $`O(r)`$ T-depth, and
$`O(2^r(r+2)^3)`$ count/work with $`O(\log(r+2))`$ T-depth.

## 4. Balanced whole-word queries

Consider an arbitrary m-bit table with $`Q=2^q`$ rows. Split its address
into halves of lengths $`r_a=\lfloor q/2\rfloor`$ and
$`r_b=\lceil q/2\rceil`$, with indicator sizes
$`H=2^{r_a}`$ and $`J=2^{r_b}`$. Then
$`H+J=O(\sqrt Q)`$ and $`HJ=Q`$.

Use the [bilinear query](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction),
with one chunked indicator for each half. Its chronological word is

```math
\mathcal B,\ I_a,\ \mathcal B^\dagger,\ I_b,\
\mathcal B,\ I_a^\dagger,\ \mathcal B^\dagger,\ I_b^\dagger.
```

The rank-reduced bilinear word has T-count $`O(m\sqrt Q)`$,
Clifford count $`O(Qm)`$, and T-depth at most $`4m`$. It restores
its coordinate changes after every output bit. Choose a positive chunk
cap a and put $`\ell_i=\min(r_i,a)`$ for each nonempty half;
a zero-length half uses its single X instead. Writing
$`r=\max(r_a,r_b)`$ and $`\ell=\min(r,a)`$ when r is nonzero,
the exact full-input XOR query consequently has

```math
T=O\!\left(\sqrt Q\,[m+(\ell+2)^3]\right),\qquad
G=O\!\left(Qm+\sqrt Q(\ell+2)^3\right),
```

```math
w=O\!\left(\sqrt Q(\ell+2)^3\right),\qquad
D_T=O\!\left(m+
\left\lceil\frac r\ell\right\rceil\log(\ell+2)\right).
```

Here w includes both indicator outputs and their work, but excludes the
query's m-bit output and original address. The two helper pools can be
reused between completed indicator calls; summing them would also change
only a constant. The zero-address query is a Clifford translation.
The $`Qm`$ Clifford term is essential: arbitrary table matrices do not
have a general $`O(\sqrt Q\,\mathrm{poly}(q,m))`$ elimination bound.

## 5. A summable budget for the final Hopf layers

Fix the accuracy parameter $`L\ge6`$, let $`N=2^n`$, and index a
layer by its distance $`k=n-d`$ from the leaves. Its padded coefficient
query has

```math
Q_k=4N2^{-k},\qquad q_k=n-k+2,\qquad m_k\le L+4+k.
```

These are the existing coefficient tables; their query outputs and
source reservations are unchanged. Each layer uses a fixed number of
these completed queries, including actual inverses. On the final R
layers, where $`1\le k\le R\le n`$, choose

```math
a_k=2^{\lfloor k/12\rfloor},\qquad
\ell_{i,k}=\min\{r_{i,k},a_k\}
```

for each nonempty address half. Put $`r_k=\lceil q_k/2\rceil`$.
The polynomial overhead obeys
$`(\min(r_k,a_k)+2)^3\le27\,2^{k/4}`$. Since
$`\sqrt{Q_k}=2\sqrt N\,2^{-k/2}`$, the total T-count is at most
a constant times

```math
\sqrt N\sum_{k=1}^{R}
\left[(L+4+k)2^{-k/2}+2^{-k/4}\right]=O(\sqrt N).
```

The constant may depend on the fixed L. The maximum live query width
is also $`O(\sqrt N)`$; queries execute successively and reuse the
same pool. Clifford counts sum to $`O(N)`$, because the first term
is bounded by $`4N\sum_k2^{-k}(L+4+k)`$ and the polynomial
overhead has the same summable bound as T-count. The existing base
reservation and query outputs fit in a sufficiently large constant
multiple of $`\sqrt N`$ for fixed L.

For the depth estimate, $`a_k\ge2^{k/12}/2`$. If the cap is active,
the indicator depth is bounded by a constant times
$`n(k+1)2^{-k/12}+\log(n+4)`$. If the chunk already contains the
entire half-address, its depth is $`O(\log(n+4))`$, so the same
bound holds. Summing the two halves and the m-bit bilinear work gives

```math
\sum_{k=1}^{R}D_{T,k}
=O\!\left(n+RL+R^2+R\log(n+4)\right).
```

Indeed, $`\sum_{k\ge1}(k+1)2^{-k/12}`$ converges. Thus for
$`R\le C_0\log(n+2)`$, with fixed $`C_0`$ and L, these final
queries have total T-depth $`O(n+\log^2(n+2))=O(n)`$, T-count
$`O(\sqrt N)`$, Clifford count $`O(N)`$, and sufficient dirty
width $`C\sqrt N`$ for a fixed constant C.

This last statement prices the final coefficient queries, including
their return. It does not by itself price the source, suffix predicates,
or queries on earlier layers. Those operations need their own schedules
before asserting a new complete-frame depth bound. The
[grouped program composition](GROUPED_PROGRAM_PREFETCH.md#8-complete-frame-theorem-at-fixed-accuracy)
supplies those schedules and combines them with this late-query bound.
No approximation, new clean qubit, or altered coefficient table enters
this replacement.

## 6. Attribution and finite evidence

The read-only conjunction is the existing
[masked dirty-counter construction](DIRTY_SUM_COMPRESSION.md), and the
whole-word extraction is the established
[bilinear indicator reduction](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction),
whose attribution remains there. Conjugating a bit flip and cancelling
arbitrary offsets are standard reversible-circuit operations. The
contribution proved here is their chunk-boundary tree composition and
the summable choice of chunk lengths for the late Hopf queries; no
general synthesis-priority claim is made.

The [bounded checks](../tests/test_chunked_dirty_indicator.py) audit
literal small conjunction phases, arbitrary dirty tree inputs, actual
inverse orientation, and the selected-leaf echo. Their small conjunction
fixtures are not used to infer the asymptotic logarithmic row depth;
that bound comes from the previously proved counter construction.
