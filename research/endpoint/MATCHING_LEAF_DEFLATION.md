# Matching-controlled canonical frames and leaf deflation

[Canonical-star identities](CANONICAL_STAR_RESOLVENT.md) · [Rank-weighted compiler](RANK_WEIGHTED_BOUNDARY_COMPILER.md) · [Borrowed real operations](CONJUGATE_SIGNAL_RECOVERY.md) · [Structured forests](STRUCTURED_FOREST_COMPILERS.md)

A matching certificate identifies the identity-relative active rank of the
canonical frame. A nonidentity leaf-matching baseline then exposes a smaller
canonical core with explicit local updates. Its algebra closes under
iteration; separate fine synthesis of every level still costs O(nN).
Neither result improves the unrestricted iterated-logarithmic upper bound.

Use the original computational-label forest $`p(j)=j-\mathrm{lowbit}(j)`$,
not the binary partition tree. Set $`N=2^n`$, $`\eta=2^{-N}`$. All native
claims use literal phases, arbitrary dirty/reference inputs, full return and
leakage error, and at most $`N+n+7`$ dirty qubits. Only T-count is optimized.

## 1. Paid canonicalization and exact active rank

Replace each original local star, with carrier column $`(c_v,u_v)^T`$ and
$`c_v\gt-1`$, by

```math
S'_v=\begin{pmatrix}c_v&-u_v^T\\u_v&I-u_vu_v^T/(1+c_v)\end{pmatrix}.
```

The [same-root transfer](CANONICAL_STAR_RESOLVENT.md#1-same-root-column-transfer)
proves $`W=W_{\rm can}C_{\rm in}`$. The correction acts on globally
disjoint immediate-child label sets, whose squared sizes total at most 3N.
It has a charged linear-T implementation; it cannot be discarded as a new
decoder. Rank statements below concern $`W_{\rm can}`$, not W. An original
four-label star can have $`\mathrm{rank}(W-I)=4`$ but canonical rank two.

Define $`L_{j,p(j)}=(u_{p(j)})_j/(1+c_{p(j)})`$, with other entries zero,
and $`F=I-L`$, $`R=F^{-1}`$, $`X=R^T-R`$. The canonical identities are

```math
\mathfrak C(L)=W_{\rm can}=(I-X)(I+X)^{-1},\qquad
F^TXF=L^T-L.
```

L is nilpotent and X is skew-symmetric, so both denominators are invertible.
Consequently

```math
\mathrm{rank}(W_{\rm can}-I)=\mathrm{rank}(L^T-L),\qquad
\ker(W_{\rm can}-I)=F\ker(L^T-L).
```

The last formula is an algebraic description, not a free implementation of
the generally nonunitary F. The algebra extends to finite real weights;
nonnegative-carrier conditioning bounds are not silently extended to all
negative carriers.

For a forest with nonzero edge weights, the real skew adjacency has rank
$`2\nu`$, where $`\nu`$ is its maximum matching size. This is the established
forest case of the [skew-rank theorem](https://arxiv.org/abs/1206.3409).
For completeness, eliminate a leaf u and its neighbor v. The ordered matrix is

```math
\begin{pmatrix}0&a&0\\-a&0&b^T\\0&-b&C\end{pmatrix},\qquad a\ne0.
```

Its leading two-by-two block has rank two and zero Schur correction to C.
A maximum matching can contain the leaf edge, by swapping it for another
edge at v when necessary. Both rank and twice the matching number therefore
decrease by two. Induction, deleting isolated vertices, gives

```math
\boxed{\mathrm{rank}(W_{\rm can}-I)=2\nu,\qquad
\dim\ker(W_{\rm can}-I)=m-2\nu.}
```

A declared support may include zero weights. Its matching size is then an
upper bound, requiring no exact-zero test on the input.

## 2. Matching-budget native compiler and its limits

For declared components C of size $`m_C`$ and matching number $`\nu_C`$,
assume $`\Gamma=\sum_Cm_C\nu_C\le\kappa N`$ for fixed $`\kappa`$.
The canonical rank-weight is at most $`2\Gamma`$. The rank-weighted theorem
and the paid input correction give complete real-frame synthesis with
$`O_\kappa(N)`$ T gates. Its real refinement uses one clean mode/recovery
flag and the original dirty budget.

Signed and singular tuples remain included. Approximate each active rotation
by rational unit-circle coefficients
$`\cos\alpha=(1-t^2)/(1+t^2)`$, $`\sin\alpha=2t/(1+t^2)`$, $`t\in\mathbb Q`$.
Finite t avoids cosine -1. Dense certified search gives an outward word
within $`\eta/4`$, preserves omitted edges, and has every local carrier
strictly greater than -1. Its canonicalization is computable without an
input equality oracle. Compile the canonical factor and input correction
each to $`\eta/4`$; total error is below $`3\eta/4`$. Small parameter
denominators can increase classical work, not the uniform native bound.

For example, take a top binomial tree with t vertices and attach every
available lower-degree child as a private leaf. The top vertices are a
vertex cover, each with a private leaf, so $`\nu=t`$. Choosing
$`t=\Theta(\sqrt{N/n})`$ gives $`m=\Theta(\sqrt{Nn})`$, hence
$`m\nu=\Theta(N)`$ while $`m^2=\Theta(Nn)`$. This improves the dense-block
certificate, not necessarily every older cut-and-dressing certificate.

### Suffix cuts cannot give a uniform constant-depth matching certificate

For outgoing-suffix cuts, the [exact dressing identity](STRUCTURED_FOREST_COMPILERS.md#2-exact-grouping-of-suffix-cuts)
is $`W=U_{\rm all}F_{k-1}\cdots F_0`$, with
$`F_r=V_r^TK_rV_r`$. If the internal matching budget is linear and the
quotient path depth k is fixed, internal components and disjoint cut stars
are each priced by a fixed number of the preceding compilers.

This is a restricted sufficient condition, not a missing universal
combinatorial lemma. The original tree has maximum undirected degree at
most n. Matching endpoints cover every edge, so
$`m_C-1\le2n\nu_C`$ and $`m_C^2\le(2n+1)m_C\nu_C`$ for nontrivial
components. The dressed word has at most $`J=3k+1`$ block stages and total
squared size at most $`3(2n+1)\Gamma+6N`$. The
[block/monomial support bound](ORDERED_BLOCK_BANKS.md#sharpness-in-the-literal-blockmonomial-model)
then forces, on the balanced complete frame,

```math
3k+1\ge\log_2n-\log_2\log_2n-O_\kappa(1).
```

This excludes a fixed-k certificate in that representation. It is not a
T-count lower bound against other mixing or a fused native circuit.

### Smallness does not create endpoint-negligible fixed directions

The full computational tree has the perfect matching of even vertices with
their odd leaves. Put $`s=\lceil\sqrt N\rceil+6`$, $`\delta=2^{-s}`$,
and give every shift edge weight $`\tau=\delta/(64N)`$. This is realized
by regular local stars with
$`c_v=(1-d_v\tau^2)/(1+d_v\tau^2)`$ and child coefficients
$`2\tau/(1+d_v\tau^2)`$. Disjoint shift columns give
$`\|L\|=\tau\sqrt n`$, hence $`\|\mathfrak C(L)-I\|\lt\delta/2`$.
In even/odd order,

```math
L^T-L=\begin{pmatrix}A_0&\tau I\\-\tau I&0\end{pmatrix},\qquad
(L^T-L)^{-1}=\begin{pmatrix}0&-I/\tau\\I/\tau&A_0/\tau^2\end{pmatrix}.
```

Since $`\|A_0\|\le2\tau\sqrt n`$ and
$`X=F^{-T}(L^T-L)F^{-1}`$, the inverse norm and Cayley singular-value
formula imply the conservative bound

```math
\sigma_{\min}(\mathfrak C(L)-I)
\ge\frac{\delta}{768N\sqrt n}\gt2^{-N}\quad(n\ge8).
```

Indeed $`\sigma_{\min}(L^T-L)\ge\tau/(1+2\sqrt n)`$,
$`\|F\|\le1+\tau\sqrt n`$, and
$`\|X\|\le2\tau\sqrt n/(1-\tau\sqrt n)`$ give that bound directly.
Any logical approximation fixing a nonzero vector has error at least this
minimum singular value. This restricts identity-relative retirement of this
canonical remainder, not another baseline or a cheap full-rank circuit.

## 3. Explicit nonidentity leaf-matching deflation

In even/odd order write

```math
L=\begin{pmatrix}A&0\\D&0\end{pmatrix},\quad D=\mathrm{diag}(d_i),
\quad T=(I+D^2)^{-1/2},\quad
G=\begin{pmatrix}T&DT\\-DT&T\end{pmatrix},\quad B=TAT.
```

A is strictly outward with at most one parent per row. G is a matching bank
of $`R_y(-\arctan d_i)`$. The leaf-only canonical baseline is $`P_0=G^{-2}`$.
Define $`\widetilde L=\begin{pmatrix}B&0\\DB&0\end{pmatrix}`$. Then

```math
\boxed{G\mathfrak C(L)G=\mathfrak C(\widetilde L).}
```

The G factors have the same orientation, not opposite orientations. To
verify it put $`M_0=I+D^2+A^TA-2A`$. The triangular formula gives

```math
\mathfrak C(L)=\mathrm{diag}(-I,I)
 +2\begin{pmatrix}I-A\\D\end{pmatrix}M_0^{-1}
   \begin{pmatrix}I-A^T&-D\end{pmatrix}.
```

Now $`G\mathrm{diag}(-I,I)G=\mathrm{diag}(-I,I)`$,
$`G[I-A;D]T=[I-B;DB]`$, and
$`TM_0T=I+B^TB+(DB)^T(DB)-2B`$. The corresponding right factor is
$`[I-B^T,-B^TD]`$. Substitution proves the identity without a physical
inverse call.

Let $`\mathcal P`$ be the declared parents of A and $`p=|\mathcal P|`$.
Only their upper coordinates remain source columns. With $`h_i=1+d_i^2`$,
the terminal part of parent v before its column scaling is

```math
w_v=\sum_{i\in\mathrm{ch}_A(v)}\left[
\mathbf1_{i\notin\mathcal P}\frac{A_{iv}}{\sqrt{h_i}}e_{i,u}
+\frac{d_iA_{iv}}{\sqrt{h_i}}e_{i,l}\right],
\quad
\gamma_v^2=\sum_{i\notin\mathcal P}A_{iv}^2
+\sum_{i\in\mathcal P}\frac{d_i^2A_{iv}^2}{1+d_i^2}.
```

The last two sums run over v's children. These vectors have disjoint
supports because each row has a unique parent. They exclude all upper
parent coordinates. A disjoint real state extension R maps one fixed
representative $`z_v`$ to $`w_v/\gamma_v`$ and fixes retained parents.
If $`\gamma_v=0`$, any completion suffices. Effective selection is addressed
below, not by an exact input-zero test.

The compact shift on the p parents and their p representatives is

```math
L'_{u,v}=\frac{A_{uv}}{\sqrt{(1+d_u^2)(1+d_v^2)}}\quad(u,v\in\mathcal P),
\qquad L'_{z_v,v}=\frac{\gamma_v}{\sqrt{1+d_v^2}},
```

with no other entries. For its computational embedding J,
$`R^T\widetilde L R=JL'J^T`$. Orthogonal covariance yields

```math
\boxed{\mathfrak C(L)=G^TR[\mathfrak C(L')\oplus I]R^TG^T,}
\qquad
P_0^T\mathfrak C(L)=(GR)[\mathfrak C(L')\oplus I](GR)^T.
```

The retired subspace undergoes the known nonidentity baseline P0. It is not
discarded. If all declared parent and matching edges are nonzero, the core
has the parent-representative perfect matching, and relative rank is exactly
2p; with zeros it is at most 2p.

## 4. Closure, native implementation, and the remaining cost

For the full computational tree, A is the same tree on N/2 vertices. Its
parents are its even vertices, so the core has N/2 coordinates. Interleave
parents with representatives and L' has exactly the smaller computational
tree topology. Recursion terminates at a matching on two modes. No generic
smaller-unitary or rank-selection oracle is hidden in this closure.

The updates obey

```math
\sum_j(L'_{jv})^2=\frac{\sum_iA_{iv}^2}{1+d_v^2},\qquad
\|L'\|_F^2+\sum_i d_i^2\le\|L\|_F^2.
```

Matching weights have summable squared norm, not a proved summable native
precision cost. The original uncanonicalized frame also permits elementary
last-layer stripping; dimension halving alone is not a complexity improvement.

G is a borrowed-signal real rotation bank. R is a disjoint real state-extension
bank, now available with the [zero-clean full-operator contract](CONJUGATE_SIGNAL_RECOVERY.md#zero-clean-real-reflection-and-state-extension-banks).
All padding and actual inverse words are charged. At one step,

```math
\|\widehat G^\dagger\widehat R\widehat C
\widehat R^\dagger\widehat G^\dagger J_2-J_2\mathfrak C(L)\|
\le2\epsilon_G+2\epsilon_R+\epsilon_C.
```

This propagates actual leakage by unitarity and does not refresh flags.
If components of the even forest have parent counts $`p_C`$ with
$`\sum_Cp_C^2\le\kappa N`$, the core is a direct sum of $`SO(2p_C)`$
blocks, so one step and the paid original input correction give a complete
linear-T compiler for that promise. This need not define a family outside
all older suffix-cut constructions.

For singular inputs, use dense certified rational-circle approximants of
the original angles and strict full-operator error margins. Declared absent
nonleaf edges stay absent; matching edges may be perturbed within budget.
Generic candidates avoid undefined normalizations, and all resulting
parameters are computable. The state-extension native bound is uniform
in the normalized vector, though classical effort can grow near singularities.

Separately emitting n levels gives $`O(n)`$ fine banks at
$`q=N+O(\log n)`$ precision, hence $`T=O(nN)`$. The literal dirty budget
fits for sufficiently large n; finitely many smaller dimensions use direct
synthesis. This bound is worse than the retained grouped compiler and is
not substituted for it. Equal-orientation G factors do not cancel, nor do
arbitrary intervening target-dependent state extensions commute.

[Focused tests](../../tests/test_boundary_refinement.py) check forest
matching, canonical rank, paid transfer examples, signed and zero-weight
one-step identities, full recursive reconstruction, relative retirement,
and the energy identity. Exact leaf Schur identities use rational arithmetic;
floating checks are matrix diagnostics, not endpoint-precision gate
certificates. The outstanding requirement is shared fine-precision work
in the actual complete word, not another rank-retirement identity.
