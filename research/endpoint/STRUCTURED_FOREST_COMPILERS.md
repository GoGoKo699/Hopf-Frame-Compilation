# Linear-T compilation of structured Hopf forests

[Disjoint-block primitive](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks) · [Precision gauge](PRECISION_GAUGE_COMPILER.md) · [Same-root star transfer](CANONICAL_STAR_RESOLVENT.md#1-same-root-column-transfer)

Two structural conditions give linear-T compilation of the complete prescribed
frame: subbinary growth in computational degree windows, and bounded exceptional
path depth outside a uniformly compilable angle family. The first includes
arbitrary angles on every fixed-outdegree active computational forest. Both
results retain the literal two-clean-qubit contract and arbitrary dirty inputs.
They are restricted-family theorems; the unrestricted A–D comparison is unchanged.

## 1. Model and native primitive

Let $`N=2^n`$ and $`\epsilon=2^{-N}`$. The computational-label tree has vertices
$`0,\ldots,N-1`$; the parent of a nonzero vertex clears its least significant
one. Write $`d(v)`$ for its number of trailing zero bits, with $`d(0)=n`$.
Its children have distinct degrees $`0,\ldots,d(v)-1`$.

An **active forest** is a subforest of this particular tree. Its edge rotations
are applied in the original Hopf layer order; omitted edges act as identities.
At each vertex the outgoing children occur in decreasing degree. An active
support certificate may include rotations whose angles happen to be zero.
This is a condition on the computational-label graph, not on the binary
partition tree used to define the frame.

The [disjoint-block theorem](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks)
supplies the following interface. For fixed $`\kappa,c>0`$, a direct sum of
$`\mathrm{SU}(m_j)`$ blocks with $`\sum_jm_j^2\le\kappa N`$ has a native
$`O_{\kappa,c}(N)`$-T approximation to error $`\epsilon/c`$, with two clean
qubits and at most $`N+n+7`$ arbitrary dirty qubits. Arbitrary prescribed
computational supports are included. Real determinant-one blocks are a special
case. All uses below are uses of this full native interface, not an inference
from matrix size alone.

If $`J`$ initializes the two clean qubits, each final compiler satisfies

```math
\left\|\widetilde U J-J(W\otimes I_{\rm dirty})\right\|\le\epsilon.
```

The norm covers every logical and dirty input and remains valid with an arbitrary
reference system. Factors reuse the same clean and dirty wires.

## 2. Exact grouping of suffix cuts

Choose some active edges as cuts, with every vertex's cuts forming an outgoing
suffix. Deleting them gives internal components $`C`$. Let $`U_C`$ be the full
internal component circuit and $`U_{\rm all}`$ their direct sum.

A cut from $`l\in C`$ to the root $`r`$ of a child component becomes the rotation
on the plane spanned by $`U_C^\dagger e_l`$ and $`e_r`$ when all internal
rotations are moved to the left of the matrix product. The completed $`U_C`$
is valid: the suffix condition ensures that no later internal outgoing edge
moves $`e_l`$, and the child component has no internal edge before its incoming
cut. Thus

```math
W=U_{\rm all}\,\mathcal D,
```

where $`\mathcal D`$ is the original chronological word of these dressed cuts,
written with later operations on the left.

Suppose the quotient forest of components has at most $`k`$ cut edges on each
path. Group dressed cuts by their parent component's quotient depth. This
regrouping is exact. Across different groups the only possibly noncommuting
pair is an incoming and outgoing cut of one component. Their chronological
order already agrees with increasing depth, so the grouping never inverts
such a pair. All inverted pairs commute. Within a group retain the original
order, including cuts sharing a boundary vertex.

Let $`V_r`$ be the direct sum of internal components at quotient depth $`r`$,
and $`K_r`$ the original cut rotations leaving those components. Then

```math
F_r=V_r^\dagger K_rV_r,
\qquad
W=U_{\rm all}F_{k-1}\cdots F_1F_0.
```

Ordinary matrix multiplication acts rightmost first. Increasing depth is the
chronological order, hence the reversed displayed product.

Each $`K_r`$ is a direct sum of ordinary stars on a boundary vertex and its cut
child roots. The supports are disjoint at fixed quotient depth. Universally,

```math
\sum_v(d(v)+1)^2
=(n+1)^2+\frac N2\sum_{d=0}^{n-1}(d+1)^2 2^{-d}
=6N-2n-5\le6N.
```

Consequently every cut-star factor has a linear squared-size budget, regardless
of the angles or the maximum active outdegree. There is no small-commutator
approximation in this factorization.

## 3. Subbinary degree-window growth

For a vertex $`v`$ of degree $`d`$ and $`1\le w\le d+1`$, consider the active
descendants reachable from $`v`$ while keeping degrees at least $`d-w+1`$.
Assume fixed constants $`C\ge1`$ and $`1\lt \gamma\lt 2`$ bound their number by
$`C\gamma^{w-1}`$, for every such vertex and window.

**Theorem.** Under this finite combinatorial certificate, every choice of real
edge angles has a complete-frame compiler with $`O_{C,\gamma}(N)`$ T gates,
two clean qubits, at most $`N+n+7`$ dirty qubits, and error $`\epsilon`$ in the
full contract of Section 1.

*Proof.* Choose a fixed integer $`s`$ satisfying

```math
\rho=\frac{\gamma^{1+1/s}}2\lt 1.
```

For each original active component, let $`d`$ be its root degree and set
$`w(d)=\lceil(d+1)/s\rceil`$. Partition its degrees, from $`d`$ down to zero,
into at most $`s`$ consecutive bands of width at most $`w(d)`$. Retain an edge
inside its component exactly when both endpoints lie in the same band.

The original child order already places within-band children before children
in lower-degree bands. Cuts are therefore suffixes without any star transfer,
angle comparison, or approximation. Every cut enters a strictly lower-degree
band, so a quotient path contains at most $`s-1`$ cuts.

Let $`m_j`$ be the resulting band-component sizes inside one original component.
They partition its vertices. The growth hypothesis gives

```math
\sum_jm_j^2
\le(\max_jm_j)\sum_jm_j
\le C^2\gamma^{w(d)-1}\gamma^d
\le C^2\gamma^{1/s}\gamma^{(1+1/s)d}.
```

There are $`N/2^{d+1}`$ nonzero computational labels of degree $`d\lt n`$.
Original component roots are a subset of these labels. The global root has
degree $`n`$ and contributes at most one additional term. Hence

```math
\begin{aligned}
\sum_{\rm all\ band\ components}m_j^2
&\le C^2\gamma^{1/s}
 \left(\frac N2\sum_{d=0}^{n-1}\rho^d
       +\gamma^{(1+1/s)n}\right)\\
&\le C^2\gamma^{1/s}
 \left(\frac1{2(1-\rho)}+1\right)N.
\end{aligned}
```

Apply Section 2. Each $`U_{\rm all}`$ or $`V_r`$ is a direct sum drawn from
this one collection of band components; each $`K_r`$ has squared-size sum at
most $`6N`$. There are at most $`1+3(s-1)=3s-2`$ native block factors. The
primitive in Section 1 prices all of them with the asserted resource bound. ∎

This theorem allows some unbounded-outdegree forests. For example, an arbitrary
star forest has at most $`w`$ vertices in a degree window and therefore admits
fixed $`C,\gamma\lt 2`$. The window condition is a finite graph property; checking
it does not require deciding equality of angles.

## 4. Fixed computational outdegree

**Corollary.** For each fixed integer $`b\ge1`$, arbitrary angles on an active
computational forest with at most $`b`$ children per vertex satisfy the theorem
of Section 3. In particular, no restriction on the angle exponent hierarchy
is required.

For $`b\ge2`$, let $`F_b(d)`$ be the maximum size of such a subtree rooted at
degree $`d`$. The distinct child degrees give

```math
F_b(d)\le1+\sum_{j=1}^{\min(b,d)}F_b(d-j).
```

Choose $`1\lt \gamma\lt 2`$ with

```math
\Delta=1-\sum_{j=1}^b\gamma^{-j}>0,
\qquad C\ge\Delta^{-1}.
```

Such a choice exists because the sum at $`\gamma=2`$ is $`1-2^{-b}`$.
Strong induction yields $`F_b(d)\le C\gamma^d`$: the proper-descendant bound
is at most $`1+C\gamma^d(1-\Delta)\le C\gamma^d`$. Subtracting the lower
endpoint of a degree window gives the same bound there. Section 3 applies.

For $`b=2`$, explicit choices are

```math
\gamma=\frac53,\qquad C=25,\qquad s=3,
\qquad \left(\frac53\right)^4=\frac{625}{81}\lt 8.
```

Thus seven disjoint-block factors suffice. Their constants are uniform in the
angles and in $`n`$; this is a factor count, not an assertion of seven T gates.

For $`b=1`$, each component is a chain with at most $`d+1`$ vertices at a root
of degree $`d`$. Its squared-size sum is already bounded by

```math
(n+1)^2+\frac N2\sum_{d=0}^{n-1}(d+1)^2 2^{-d}\le6N.
```

The complete target is therefore one directly priced block sum.

Exact zero, negative, and corner angles are harmless: the factorization uses
ordinary real determinant-one blocks and never divides by a branch amplitude.
If support is inferred from an oracle, only rotations certified sufficiently
close to the **identity gate** may be omitted within the error budget. A
nonidentity $`\pi/2`$ or $`\pi`$ corner must remain active unless an explicitly
tracked signed permutation is applied.

The condition is specific to the original computational-label tree. Replacing
it by an arbitrary binary topology does not preserve the distinct-child-degree
count. A full binomial component rooted at degree $`h`$ has $`2^h`$ vertices
and $`h`$ active root children; growing such components do not satisfy a fixed
outdegree or fixed subbinary window-growth certificate.

## 5. Bounded exceptional path depth

A different criterion permits unbounded computational branching. Work after
signed normalization with $`c_e\ge1/\sqrt2`$ and
$`0\le s_e\le1/\sqrt2`$. Delete known identity edges. Mark remaining edges as
base or exceptional, and suppose every computational-tree path contains at
most a fixed $`k`$ exceptional edges.

Either of the following base promises suffices:

| Base-angle condition | Native bound |
|---|---|
| $`s_e\ge2^{-E_L}`$, with $`E_L\le A N/n`$ for fixed $`A`$ | $`O_{A,k}(N)`$ T |
| Smaller-branch exponents in $`[E_{\min},E_{\max}]`$, with $`E_{\max}\le K E_{\min}`$ for fixed $`K`$ | $`O_{K,k}(N)`$ T |

Exceptional angles have no magnitude restriction. Their exponents may occupy
arbitrarily many intermediate bands. Both bounds retain the full Section 1
contract.

At each vertex stably place base children before exceptional children. The
[same-root star transfer](CANONICAL_STAR_RESOLVENT.md#1-same-root-column-transfer)
gives $`W=W'C_{\rm in}`$, where the input correction is a direct sum on
immediate-child labels. Its squared-size sum is at most
$`\sum_vd(v)^2=3N-2n-3`$, so it is one priced factor.

Here are the coefficient bounds needed to price the base components. In a
normalized degree-$`d`$ star the center root coefficient is at least
$`2^{-d/2}`$. Each original leaf coefficient is its sine times a product of
earlier cosines, hence lies between $`2^{-d/2}s_e`$ and $`s_e`$. During inverse
stick breaking every remaining center norm lies between $`2^{-d/2}`$ and one.
Therefore the reordered base coefficients obey

```math
2^{-d/2}s_e\le s'_e\le2^{d/2}s_e,
\qquad c'_e\ge2^{-d/2}.
```

The order is stable within base children. After exceptional edges are deleted,
each retained child has no more preceding siblings than before. Its inherited
binary decision depth therefore remains at most $`n`$. The underlying
computational parent-child graph is unchanged, so exceptional path depth
remains at most $`k`$.

Delete exceptional suffixes and apply the exact factorization in Section 2 to
$`W'`$. Each $`K_r`$ is a priced star sum. For the first base promise, each
retained forest needs per-depth coarse precision at most

```math
e_d\le E_L+\frac n2+O(\log n),
\qquad \sum_de_d=O_A(N).
```

The [precision-trace theorem](PRECISION_GAUGE_COMPILER.md#5-precision-trace-and-native-resource-bound)
then prices the entire forest. In particular, with $`\sum_dS_d\le N`$,

```math
\sum_d\sqrt{S_de_d}
\le\sqrt{\left(\sum_dS_d\right)\left(\sum_de_d\right)}
=O_A(N).
```

For the second promise, the retained sine exponents lie in
$`[E_{\min}-n/2,E_{\max}+n/2]`$. If $`E_{\min}`$ exceeds a sufficiently large
fixed multiple of $`n`$, their ratio remains bounded in terms of $`K`$ and the
retained cosines stay regular. The [fixed-ratio theorem](PRECISION_GAUGE_COMPILER.md#6-fixed-exponent-ratio-families)
applies. Otherwise $`E_{\max}=O_K(n)`$ and both retained branches are at least
$`2^{-O_K(n)}`$, so the regular-chart theorem applies.

There are at most $`2+3k`$ priced factors including $`C_{\rm in}`$. All-column
star transfer and full initialized-isometry composition preserve the prescribed
frame, rather than only its root state.

The marks can be supplied as a finite certificate. For the first promise a
useful robust criterion is that nonzero edges with $`s_e\lt t=2^{-AN/n}`$ have
path count at most $`k`$. Overlapping tests assign base when $`s_e>t/2`$ or
exceptional when $`s_e\lt t`$. One test always succeeds; selected exceptional
edges belong to the promised set and base precision loses at most one bit.
Certified near-identity rounding handles zero branches. If further finite-grid
rounding is needed for normalization, classify first and preserve the labels
under a sufficiently small perturbation. The combinatorial promise is on the
chosen normalized labeled tree.

## 6. Error, workspace, and verification

Allocate a fixed fraction of $`\epsilon`$ to each native factor and to any
required input rounding. For a fixed number of factors the fine-source
precision is $`q=N+O(1)`$, with constants depending on the stated structural
parameters. The primitive's two-clean and $`N+n+7`$-dirty reservation applies
for sufficiently large $`n`$; fixed smaller dimensions use its finite fallback.

If $`\|U_jJ-JW_j\|\le\eta_j`$, unitarity and a telescoping hybrid give
$`\|U_m\cdots U_1J-JW_m\cdots W_1\|\le\sum_j\eta_j`$. Actual adjoints
satisfy the corresponding inverse bound. This proof does not treat approximately
returned work as fresh clean input. All phases and dirty/reference states are
included in the operator norm.

The [combinatorial helper](../../compiler_robust_hopf/structured_forests.py)
returns degree bands, suffix-cut components, chronological cut groups, and
exact finite squared-size budgets. It does not emit native gates. The
[tests](../../tests/test_structured_forests.py) compare independently assembled
complete edge matrices with the grouped factors, include signed and exact
corner angles, verify same-root input correction on every column, and guard
against reversing the matrix-product order. Graph tests check the original
degree counts and the binary-growth certificate. These diagnostics support the
identities; the analytical estimates above establish the asymptotic counts.
