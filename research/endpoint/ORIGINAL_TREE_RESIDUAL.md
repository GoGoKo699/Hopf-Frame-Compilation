# Original-tree Cayley identities and a small collective residual

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Exact dirty interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter) · [Disjoint-block compiler](DISJOINT_BLOCK_COMPILER.md)

The prescribed, unmodified Hopf word has an explicit Cayley transform with
rank-one tree cuts. A native coarse word reduces it to a small unitary
right residual with rank-one cuts, using linear T count and exact
dirty return.
Neither statement supplies a linear-T compiler for the fine residual.
The source-completion theorem still needed is stated in Section 6.
The existing [complex left-residual recursion](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation)
already gives nested rank-two sibling blocks and constant-size local
generators. Section 4 gives a separate all-edge-cut argument for the
right residual, whose outward ordering sharpens the cut rank to one;
it does not replace that recursion.

This chapter concerns the **original complete operator**, with every
supplied angle retained. It makes neither the exterior phase/permutation
change of the [fixed-tree Cayley reduction](TREE_CAYLEY_REDUCTION.md) nor
the local-star replacement and input correction of the
[canonical-star resolvent](CANONICAL_STAR_RESOLVENT.md#1-same-root-column-transfer).
The [coupled resolvent](COUPLED_TREE_RESOLVENT.md) treats a different
coarse/fine transport representation.

## 1. The original ordered tree word

Let a tree have $`N`$ computational labels. Orient each edge
$`e=(a_e,b_e)`$ and define

```math
K_e=|b_e\rangle\langle a_e|-|a_e\rangle\langle b_e|,
\qquad R_e=\exp(\theta_eK_e).
```

Use every edge once, at distinct application time $`\tau_e`$, and set
$`W=R_{\rm last}\cdots R_{\rm first}`$. Products act rightmost first.
The supplied Hopf word is this construction on the tree
$`p(j)=j-\mathrm{lowbit}(j)`$, with rotations on
$`p0\,0^{n-d-1},p1\,0^{n-d-1}`$ at depth $`d`$. Disjoint gates within
a depth may be given any order.

**Original-tree Cayley theorem.** For every real angle tuple,

```math
\boxed{\det(I+W)=2^N\prod_e\cos^2(\theta_e/2).}
\qquad\text{(1)}
```

Thus $`I+W`$ is singular exactly when an edge angle is
$`\pi\bmod 2\pi`$. When it is invertible, put

```math
S=(W-I)(W+I)^{-1},\qquad t_e=\tan(\theta_e/2).
```

The matrix $`S`$ is real skew-symmetric. For the unique path
$`u=v_0,\ldots,v_k=v`$, write
$`\sigma_i=(K_{e_i})_{v_{i-1},v_i}\in\{-1,1\}`$. Then

```math
\boxed{
S_{uv}=
\left(\prod_{i=1}^k\sigma_i t_{e_i}\right)
\left(\prod_{i=1}^{k-1}
\mathrm{sgn}(\tau_{e_i}-\tau_{e_{i+1}})\right),
\qquad S_{uu}=0.
}
\qquad\text{(2)}
```

In particular, $`\lvert S_{uv}\rvert=\prod_{e\in\mathrm{path}(u,v)}\lvert t_e\rvert`$.
Across every tree-edge cut $`A\mid B`$, both off-diagonal blocks of $`S`$
have rank at most one. Zero angles are included.

**Proof.** A prefix product $`Q`$ is block diagonal on the components
of its edge forest. The next edge joins distinct components. Suppose
first that all half-angle tangents used are finite, and use $`S`$ for
the prefix Cayley transform. For $`F=[e_a,e_b]`$,

```math
F^TSF=0,\qquad
(I+Q)^{-1}=(I-S)/2,\qquad
Q(I+Q)^{-1}=(I+S)/2.
```

Write $`R=I+F(r-I_2)F^T`$. The determinant lemma gives

```math
\frac{\det(I+RQ)}{\det(I+Q)}
=\det\!\left(I_2+(r-I_2)/2\right)
=\cos^2(\theta/2).
```

Induction starts at $`Q=I`$, and every displayed ratio is nonzero.
This proves (1) on that chart. Continuity extends the determinant
identity to all angles and proves the asserted singular set.
Woodbury gives, on the chart,

```math
(I+RQ)^{-1}
=(I+Q)^{-1}-\tfrac12(I-S)tK(I+S),
\qquad
\boxed{S_{\rm new}=S+t(I-S)K(I+S).}
\qquad\text{(3)}
```

The correction has only cross blocks between the two components, each
of rank at most one. To obtain (2), split a path at its latest edge and
apply (3); the two remaining subpaths lie in the earlier components.
The signs at the two joins are exactly the adjacent-time comparisons
in (2). This induction also covers zero $`t_e`$. Reversing the path
changes the displayed sign by minus one, as required by skewness.

For a fixed cut edge, split every crossing path at that edge. All edge
factors and adjacent-time comparisons split into a factor depending only
on $`u\in A`$ and one depending only on $`v\in B`$. Hence
$`S_{AB}`$ is an outer product, including when edge times interleave.
This proves the cut-rank assertion. The singular-set characterization
is a mathematical identity, not an instruction to test exact equality
of computable inputs; the residual construction below needs no such test.

## 2. A regular-angle conditioning adversary

Take homogeneous Hopf angles with $`t_e=t\gt0`$, and denote their
Cayley matrix by $`S_n`$. Splitting the first-bit halves gives

```math
S_n=
\begin{pmatrix}S_{n-1}&0\\0&S_{n-1}\end{pmatrix}
+t\begin{pmatrix}0&-vv^T\\vv^T&0\end{pmatrix},
\qquad v=(I+S_{n-1})e_0.
\qquad\text{(4)}
```

Indeed, the root rotation is applied first. Apply (3) to the transposed
word, in which it is last, and transpose back. The root-to-label paths
in an $`m`$-bit binomial tree are indexed by subsets of its $`m`$ depths.
Equation (2) therefore gives

```math
v_j=t^{\mathrm{popcount}(j)},\qquad
\|v\|^2=(1+t^2)^{n-1}.
```

The norm of the root cross block and the triangle inequality in (4)
give, respectively,

```math
t(1+t^2)^{n-1}\le\|S_n\|_2
\le\frac{(1+t^2)^n-1}{t}.
\qquad\text{(5)}
```

Let $`s_1,\ldots,s_{N/2}`$ contain one singular value from each
skew-symmetric pair. Equation (1) implies

```math
\det(I-S_n)=(1+t^2)^{N-1}
=\prod_{j=1}^{N/2}(1+s_j^2),
\qquad
1+s_{\min}^2\le(1+t^2)^{2(N-1)/N}\le(1+t^2)^2.
```

Since the singular values of $`I\pm S_n`$ are
$`\sqrt{1+s_j^2}`$, this proves the actual condition-number bounds

```math
\frac{\sqrt{1+t^2(1+t^2)^{2n-2}}}{1+t^2}
\le\kappa_2(I\pm S_n)
\le\sqrt{1+\left(\frac{(1+t^2)^n-1}{t}\right)^2}.
\qquad\text{(6)}
```

For fixed $`t\gt0`$, the condition number is consequently
$`\Theta_t(N^{\log_2(1+t^2)})`$. For $`0\lt t\le1`$, there is also
the exact absolute-row-sum formula

```math
\|S_n\|_\infty=(1+t)^n-1.
\qquad\text{(7)}
```

The root row attains this value by (2). Inductively, any other row in
(4) is bounded by the previous maximum plus
$`t\lvert v_i\rvert\sum_j\lvert v_j\rvert\le t(1+t)^{n-1}`$, proving the upper bound.
At $`\theta=\pi/4`$, (7) is $`\sqrt N-1`$, while (6) grows as
$`N^{\log_2(4-2\sqrt2)}`$. This regular-angle family already has a
[linear-T compiler](PRECISION_GAUGE_COMPILER.md); these bounds obstruct
an unnormalized use of this Cayley representation, not the target
itself. Contractivity of $`(I-S)^{-1}`$ does not price a normalized
block encoding, a source, or its native inverse.

## 3. An exactly logical native coarse word

Set $`N=2^n`$, $`n\ge3`$, and

```math
\delta=\frac1{8\sqrt N},\qquad
\eta=\min\{1/128,\delta/2\},\qquad
\epsilon_d=\eta\,2^{d-n}.
\qquad\text{(8)}
```

There is a native word $`\widehat Q`$ using neither initialized flag,
such that, on its entire logical-and-dirty input space,

```math
\widehat Q=Q\otimes I_{\rm dirty},\qquad
\|W-Q\|\lt\delta.
\qquad\text{(9)}
```

Here $`Q`$ is a product of exact, possibly complex, two-level SU(2)
blocks on the original edges in the original layer order. Exact reality
of those native blocks is not assumed.

Import the [borrowed-workspace theorem](../../docs/BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements)
with zero clean work, $`b=N+n+7`$, and the tolerance $`\eta`$ in (8).
Its hypothesis $`0\lt\eta\le1/64`$ holds. With
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}=O(n)`$, it gives

```math
T(\widehat Q)
=O\!\left(\frac{NL}{n+b}+L\sqrt N\right)
=O(n\sqrt N)=O(N),\qquad G(\widehat Q)=O(Nn).
\qquad\text{(10)}
```

The theorem applies to every such $`n`$ and dirty budget, so no separate
small-dimension fallback or unquantified workspace threshold is needed.
Both supplied clean flags remain untouched. Its full-space return
identity includes arbitrary references.

The required logical support follows from that theorem's
[sector construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector),
not merely from its norm guarantee. For an addressed quarter-angle word
$`P_x`$, the actual determinant-one approximant is

```math
C_x=XP_x^\dagger XP_x,\qquad D_x=C_x^2,
\qquad \det C_x=1,\quad XC_xX=C_x^\dagger.
```

The scalar phase of $`P_x`$ cancels. The echo emits $`D_x`$ on the
selected original two-dimensional sector and exact identity on every
inactive sector, even when $`D_x`$ is complex. With each depth's block
error at most $`\epsilon_d`$, disjoint-block errors take a maximum,
and layer errors telescope to
$`\|W-Q\|\le\sum_d\epsilon_d\lt\eta\lt\delta`$.
The construction reserves a logical suffix control when one exists;
at the final depth it reserves a prefix control and processes both
sectors. These guards therefore do not enlarge the table to all
logical bits or introduce an initialized helper. The
[disjoint-block chapter](DISJOINT_BLOCK_COMPILER.md#exactly-logical-determinant-one-approximants)
uses the same exact logical coarse-word interface.

Certified strict-error search for the finite one-qubit words uses the
supplied sine/cosine enclosures. It need not select a quarter-angle
branch or decide a zero or rounding tie exactly.

## 4. Small residual and tree-cut rank

The [complex left-residual recursion](ENDPOINT_TREE_TRANSPORT.md#the-actual-complex-coarse-residual)
retains its explicit nested generators and local inverse bounds. The
following right-residual proof instead applies to every original edge
cut without computing those generators.

Let $`Q`$ be any product using one two-level unitary on
each original tree edge and suppose $`\|W-Q\|\le\delta\lt2`$.
It need not be real or use the same edge order. Define

```math
E=WQ^\dagger,\qquad S_E=(E-I)(E+I)^{-1}.
```

Then

```math
S_E^\dagger=-S_E,\qquad
\|S_E\|\le\frac{\delta}{2-\delta},\qquad
\mathrm{rank}(E_{AB}),\quad
\mathrm{rank}((S_E)_{AB})\le2
\qquad\text{(11)}
```

for every tree-edge cut $`A\mid B`$, also in the reverse direction.

**Proof.** Isolate the sole gate crossing a fixed cut. Every gate before
and after it is block diagonal across that cut. Thus both off-diagonal
blocks of $`W`$, and separately of $`Q`$, have rank at most one. In

```math
E_{AB}=W_{AA}Q_{BA}^\dagger+W_{AB}Q_{BB}^\dagger
```

each summand has rank at most one. Since
$`E+I=2I+(E-I)`$ and $`\|E-I\|\le\delta\lt2`$, all its principal
diagonal blocks are invertible. The block-inverse formula expresses
each inverse off-diagonal block as invertible left and right factors
times the original off-diagonal block, preserving the rank bound.
Now use $`S_E=I-2(E+I)^{-1}`$. Unitarity of $`E`$ gives
anti-Hermiticity, and $`\|(E+I)^{-1}\|\le(2-\delta)^{-1}`$ gives
the norm bound.

**Outward-order sharpening.** For the actual $`W,Q`$ of Section 3,
the cut ranks in (11) are at most **one**, including with complex native
blocks. More generally, this holds whenever both edge orders are outward
from the same rooted tree: every parent edge precedes all edges below it.

To prove this, let $`B`$ be the descendant component of a cut edge,
with root $`b`$, and let $`A`$ be its complement. Before the crossing
rotation, no gate has acted on $`B`$. Thus

```math
W=D_W R_W(P_W\oplus I_B),\qquad
Q=D_Q R_Q(P_Q\oplus I_B),
```

where $`D_W,D_Q`$ are block diagonal across the cut, and $`R_W,R_Q`$
act only on its two endpoint labels. Consequently

```math
E=D_W R_W(P_WP_Q^\dagger\oplus I_B)R_Q^\dagger D_Q^\dagger.
```

The middle product acts as identity on
$`\mathrm{span}\{e_j:j\in B,\ j\ne b\}`$. Its $`AB`$ block
therefore has only the $`b`$ column, and its $`BA`$ block only the
$`b`$ row. Multiplication by the outer block-diagonal factors preserves
rank at most one. The preceding block-inverse argument gives the same
bound for $`S_E`$. This is an orientation-specific property of
$`WQ^\dagger`$; it is not inferred for the left residual
$`Q^\dagger W`$ or for arbitrary reordered edge words.

For the choice (8), row and column Cauchy--Schwarz additionally give

```math
\|S_E\|\le\frac1{16\sqrt N-1},\qquad
\|E\|_1,\|E\|_\infty\le1+\sqrt N\delta=\frac98.
\qquad\text{(12)}
```

All these matrices are computable from the supplied enclosures and the
finite exact word $`Q`$. A uniformly convergent Neumann series around
$`2I`$ computes $`(E+I)^{-1}`$. Crossing-gate factors and the block-inverse
formula provide the cut factorizations without selecting nonzero minors.
For anti-Hermitian $`S,T`$, Cayley reconstruction is globally stable:

```math
\mathcal C(S)=(I+S)(I-S)^{-1},\qquad
\|\mathcal C(S)-\mathcal C(T)\|\le2\|S-T\|,
\qquad\text{(13)}
```

by the resolvent identity and $`\|(I-S)^{-1}\|,\|(I-T)^{-1}\|\le1`$.

This is the exact logical reduction $`W=EQ`$. The coarse discrepancy
is part of $`E`$, not an additional final error. A complete-isometry
approximation of $`E`$ followed on the right by $`\widehat Q`$ has the
same approximation error for $`W`$. The rank statements concern logical
matrices. Physical errors of a later approximate native circuit may
contain dirty operators and need not have these ranks; no recursive
closure of physical corrections is asserted.

## 5. Shared factors, contraction slack, and the completion gap

The small residual has an explicit low-dimensional algebraic source
factorization, but this does not yet satisfy the native resource contract.
Split the computational labels by their binary prefixes. Each sibling
off-diagonal block is a submatrix of an original tree-edge cut, so its
rank is at most one for the actual outward-ordered residual. Decompose
both directed blocks at each split into balanced singular-value factors,
and handle the diagonal separately.
Every matrix entry then belongs to exactly one such block or to the
diagonal. This constructs

```math
E=AB^\dagger,\qquad A,B\in\mathbb C^{N\times m},\qquad m\le3N-2.
\qquad\text{(14)}
```

The mode count is $`N`$ diagonal terms and two directed rank-one
terms for each of the $`N-1`$ binary splits. Zero terms may simply be
padded with zero columns; no nonzero-rank test is required.

For a diagonal entry $`z=E_{ii}`$, take factors
$`z/\sqrt{\lvert z\rvert}`$ and $`\sqrt{\lvert z\rvert}`$, defining
both as zero at $`z=0`$. Its Gram contribution is
$`\lvert z\rvert\le1`$ on each side. At one hierarchy level,
supports on either factor side are disjoint. Each
balanced off-diagonal Gram block has norm at most
$`\|E-I\|\le\delta`$. Hence

```math
AA^\dagger,BB^\dagger\preceq(1+n\delta)I,
\qquad n\delta\le\frac3{16\sqrt2}\lt\frac14.
\qquad\text{(15)}
```

Put $`F=A/\sqrt2`$ and $`G=B/\sqrt2`$. Thus
$`FG^\dagger=E/2`$ and both Gram deficits are bounded below by
$`3I/8`$. The standard orthogonal completions are the isometries

```math
V_A=\begin{pmatrix}F^\dagger\\(I-FF^\dagger)^{1/2}\\0\end{pmatrix},
\qquad
V_B=\begin{pmatrix}G^\dagger\\0\\(I-GG^\dagger)^{1/2}\end{pmatrix},
\qquad V_A^\dagger V_B=E/2.
\qquad\text{(16)}
```

They require up to $`m+2N=5N-2`$ modes. This explicit reservation
exceeds the $`4N`$ logical-and-two-clean mode space. Arbitrary dirty
qubits are not additional initialized source modes. Moreover, the two
Gram-deficit square roots and unitary extensions in (16) have not been
assigned native costs. The bound is on this construction, not a lower
bound on possible dilation dimension; abstract smaller dilations exist.

For the arbitrary-order rank-two generalization, the corresponding
shared count is $`4N-4`$: $`N`$ diagonal terms, $`N`$ bottom-level
directed terms, and at most $`2N-4`$ higher-level terms. Its separate
failure-bank completion uses $`6N-4`$ modes. The actual coarse word
retains the sharper $`3N-2`$ and $`5N-2`$ reservations above.

The factorization is an algebraic existence statement. It does not
assume an exact singular-value or zero-test oracle: finite-accuracy
factors can be found by certified search with the strict Gram slack in
(15). Such classical computation still supplies no quantum source.

## 6. Sufficient native source theorem, still unproved

A sufficient remaining theorem would produce a **joint native source
overlap** for every outward-ordered residual in (8)--(12). At the ideal
logical level, two source unitaries on $`4N`$ modes would satisfy

```math
J_2^\dagger A_E^\dagger B_EJ_2=E/2.
\qquad\text{(17)}
```

The required implementation statement is stronger than this matrix
identity. With $`\varepsilon=2^{-N}`$, its actual unitary overlap word
$`\widehat U`$ must obey

```math
T(\widehat U)=O(N),\qquad b\le N+n+7,\qquad
\left\|2J_2^\dagger\widehat UJ_2-E\otimes I_b\right\|
\le\varepsilon/4.
\qquad\text{(18)}
```

Only the stated two clean qubits may be initialized. The cost includes
source preparation, both orthogonal complements, all controls and
lookups, and the actual inverse words. The norm in (18) is on every
logical and dirty input, with arbitrary references. A constant number
of existing priced diagonal or disjoint-block slots would suffice only
with a legal flag and dirty allocation for every slot; an occupied mode
bit cannot also be assumed to be a fresh clean compiler flag.

If (18) holds, use the existing
[robust normalization-two lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
with

```math
\mathcal R=I-2J_2J_2^\dagger,\qquad
\widehat{\mathcal A}
=-\widehat U\mathcal R\widehat U^\dagger\mathcal R\widehat U.
```

It gives $`\|\widehat{\mathcal A}J_2-J_2(E\otimes I_b)\|\le\varepsilon`$,
including rejected-space leakage and dirty return. The reflection about
clean $`00`$ is an X-conjugated CZ, and the leading minus is a literal
Clifford scalar. Two forward overlap calls and one actual inverse have
linear total cost. Composing with (9) then gives the requested complete
frame at the same error and width.

No construction establishing (18) is proved here. Rank-one tree cuts,
constant row sums, classical evaluation, and the contraction slack in
(15) do not establish the native source theorem.

## 7. Verification scope

[Reference routines](../../compiler_robust_hopf/original_tree_residual.py)
and [focused tests](../../tests/test_original_tree_residual.py) exercise
the determinant and signed path identities, singular corners, tree-cut
ranks, homogeneous norm formulas, complex nearby residuals, and the
shared-factor identities, Gram bounds, and completion mode counts.
These finite matrix
checks supplement the proofs. They do not emit the coarse native word,
prove a source compiler, or change the endpoint bounds.
