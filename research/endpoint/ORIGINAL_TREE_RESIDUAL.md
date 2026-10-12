# Original-tree Cayley identities and a small collective residual

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Exact dirty interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter) · [Ordered banks](ORDERED_BLOCK_BANKS.md)

The prescribed Hopf word has explicit Cayley identities and rank-one tree
cuts. An exactly returned native coarse word leaves an exponentially small
right residual at linear T count. A regularized joint-Gram construction
completes its canonical factors on exactly $`4N`$ modes, with uniformly
bounded conditioning and cut ranks. The fine native source and its inverse
remain unpriced; the unrestricted endpoint bounds are unchanged.

This chapter preserves every supplied angle and logical column. It does not
make the exterior phase/permutation change of the
[fixed-tree Cayley reduction](TREE_CAYLEY_REDUCTION.md) or the local-star
replacement of the [canonical-star resolvent](CANONICAL_STAR_RESOLVENT.md#1-same-root-column-transfer).
The [complex left-residual recursion](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation)
has its separate nested rank-two generators; the right-residual sharpening
below depends on outward ordering and does not replace that recursion.

## 1. The original ordered tree word

Let a tree have $`N`$ computational labels. Orient each edge
$`e=(a_e,b_e)`$ and set

```math
K_e=|b_e\rangle\langle a_e|-|a_e\rangle\langle b_e|,
\qquad R_e=\exp(\theta_eK_e).
```

Use every edge once at distinct application time $`\tau_e`$ and put
$`W=R_{\rm last}\cdots R_{\rm first}`$. Products act rightmost first.
The Hopf word uses the tree $`p(j)=j-\mathrm{lowbit}(j)`$, with depth-d
rotations on $`p0\,0^{n-d-1},p1\,0^{n-d-1}`$. Disjoint gates within a
depth may be ordered arbitrarily.

**Original-tree Cayley theorem.** For every real angle tuple,

```math
\boxed{\det(I+W)=2^N\prod_e\cos^2(\theta_e/2).} \tag{1}
```

Thus the singular set is exactly an edge angle $`\pi\bmod2\pi`$. Off
that set, let $`S=(W-I)(W+I)^{-1}`$ and $`t_e=\tan(\theta_e/2)`$.
For the unique path $`u=v_0,\ldots,v_k=v`$, put
$`\sigma_i=(K_{e_i})_{v_{i-1},v_i}\in\{-1,1\}`$. Then

```math
S_{uv}=\left(\prod_{i=1}^k\sigma_it_{e_i}\right)
       \left(\prod_{i=1}^{k-1}\mathrm{sgn}(\tau_{e_i}-\tau_{e_{i+1}})\right),
\qquad S_{uu}=0. \tag{2}
```

The matrix is real skew-symmetric, and both crossing blocks have rank at
most one across every tree-edge cut. Zero tangents are included.

**Proof.** A prefix product $`Q`$ is block diagonal on the components of
its edge forest. The next edge joins distinct components. With
$`F=[e_a,e_b]`$ and the prefix Cayley matrix $`S`$,

```math
F^TSF=0,\qquad (I+Q)^{-1}=(I-S)/2,
\qquad Q(I+Q)^{-1}=(I+S)/2.
```

Writing $`R=I+F(r-I_2)F^T`$, the determinant lemma gives

```math
\frac{\det(I+RQ)}{\det(I+Q)}
=\det(I_2+(r-I_2)/2)=\cos^2(\theta/2).
```

Induction starts at $`Q=I`$. On the finite-tangent chart each ratio is
nonzero; continuity gives (1) everywhere and the claimed singular set.
Woodbury gives

```math
(I+RQ)^{-1}=(I+Q)^{-1}-\tfrac12(I-S)tK(I+S),
\qquad S_{\rm new}=S+t(I-S)K(I+S). \tag{3}
```

The correction has only rank-one cross blocks between the joined
components. Split a path at its latest edge and apply (3). Induction on
the two subpaths yields (2), including the adjacent-time comparisons at
the joins and zero tangents. Reversing the path reverses its sign.
For a fixed cut edge, split every crossing path there. Both edge products
and the adjacent-time signs split into a factor depending only on its
endpoint on each side, giving an outer product. The singular-set identity
is not an instruction to test exact equality of computable inputs.

## 2. A regular-angle conditioning adversary

For homogeneous $`t_e=t>0`$, split the first-bit halves:

```math
S_n=\begin{pmatrix}S_{n-1}&0\\0&S_{n-1}\end{pmatrix}
+t\begin{pmatrix}0&-vv^T\\vv^T&0\end{pmatrix},
\qquad v=(I+S_{n-1})e_0. \tag{4}
```

The root rotation is first; applying (3) to the transposed word, where
it is last, proves (4). Root paths are indexed by subsets of the remaining
depths, so (2) gives $`v_j=t^{\mathrm{popcount}(j)}`$ and
$`\|v\|^2=(1+t^2)^{n-1}`$. The root cross block and a geometric
triangle bound imply

```math
t(1+t^2)^{n-1}\le\|S_n\|_2\le\frac{(1+t^2)^n-1}{t}. \tag{5}
```

Let $`s_1,\ldots,s_{N/2}`$ contain one singular value from each skew
pair. Equation (1) implies

```math
\prod_j(1+s_j^2)=\det(I-S_n)=(1+t^2)^{N-1},
\qquad 1+s_{\min}^2\le(1+t^2)^{2(N-1)/N}\le(1+t^2)^2.
```

Since the singular values of $`I\pm S_n`$ are
$`\sqrt{1+s_j^2}`$,

```math
\frac{\sqrt{1+t^2(1+t^2)^{2n-2}}}{1+t^2}
\le\kappa_2(I\pm S_n)
\le\sqrt{1+\left(\frac{(1+t^2)^n-1}{t}\right)^2}. \tag{6}
```

For fixed t this is $`\Theta_t(N^{\log_2(1+t^2)})`$. If
$`0<t\le1`$, also

```math
\|S_n\|_\infty=(1+t)^n-1. \tag{7}
```

The root row attains (7). Inductively any other row in (4) is bounded by
the previous maximum plus
$`t|v_i|\sum_j|v_j|\le t(1+t)^{n-1}`$. At $`\theta=\pi/4`$,
(7) is $`\sqrt N-1`$ and (6) grows as
$`N^{\log_2(4-2\sqrt2)}`$. This regular family already has a
[linear-T compiler](PRECISION_GAUGE_COMPILER.md). The example restricts
an unnormalized Cayley implementation, not the target. Contractivity of
the inverse does not supply a charged block encoding or source.

## 3. An exactly logical native coarse word

The [borrowed-workspace theorem](../../docs/BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements)
applies at arbitrary accuracy, using no initialized work. At an integer
precision exponent $`s\ge6`$, with $`b=N+n+7`$, it emits

```math
\widehat Q_s=Q_s\otimes I_b,\qquad
\|W-Q_s\|<2^{-s},\qquad
T(\widehat Q_s)=O\!\left(\frac{Ns}{N+2n+7}+s\sqrt N\right),
\qquad G(\widehat Q_s)=O(Ns). \tag{8}
```

The first equality is exact on the entire input space, including dirty
references. Both supplied clean flags remain unused. The actual logical
word uses possibly complex determinant-one blocks on the original edges
in their original outward layer order.

This support property follows from the theorem's
[sector echo](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector),
not merely from its error bound. Its quarter-angle word $`P_x`$ produces

```math
C_x=XP_x^\dagger XP_x,\quad D_x=C_x^2,\quad
\det C_x=1,\quad XC_xX=C_x^\dagger.
```

The scalar phase of $`P_x`$ cancels. The actual echo emits $`D_x`$ on
the selected two-dimensional sector and exact identity elsewhere. A
logical suffix supplies the guard when present; at the final depth a
prefix bit is reserved and both sectors are processed. No initialized
helper or all-logical-bit table expansion is introduced. Certified finite
word search uses strict sine/cosine error enclosures, not a quarter-angle
branch decision, exact zero test, or rounding-tie oracle.

Choose the layer tolerances $`\epsilon_d=2^{-s}2^{d-n}`$. Disjoint
block errors take a maximum, and the layer errors sum to less than
$`2^{-s}`$. The former inverse-polynomial neighborhood is obtained by
$`s=\max\{7,\lceil\log_2(16\sqrt N)\rceil\}`$, at
$`T=O(n\sqrt N)`$, $`G=O(Nn)`$.

**Stronger coarse-precision corollary.** For $`n\ge3`$, take

```math
s=\lceil\sqrt N\rceil+6,\qquad\delta=2^{-s}.
\quad T(\widehat Q_s)=O(N),\qquad G(\widehat Q_s)=O(N^{3/2}). \tag{9}
```

This is a parameter substitution in the existing native theorem, not a
new independently synthesized fine correction. In particular
$`\delta<1/(8\sqrt N)`$ for all permitted dimensions.

## 4. Small residual and tree-cut rank

Let Q use one two-level unitary on each original edge and suppose
$`\|W-Q\|\le\delta<2`$. For now its order may differ. Set

```math
E=WQ^\dagger,\qquad S_E=(E-I)(E+I)^{-1}.
```

Then, across every original edge cut $`A\mid B`$ in either direction,

```math
S_E^\dagger=-S_E,\qquad
\|S_E\|\le\frac{\delta}{2-\delta},\qquad
\mathrm{rank}E_{AB},\ \mathrm{rank}(S_E)_{AB}\le2. \tag{10}
```

**Proof.** Isolate the only gate crossing that cut. Earlier and later
gates are block diagonal, so the cross blocks of W and Q each have rank
at most one. In
$`E_{AB}=W_{AA}Q_{BA}^\dagger+W_{AB}Q_{BB}^\dagger`$ the two terms
have rank at most one. All principal blocks of $`E+I=2I+(E-I)`$ are
invertible. The block-inverse formula expresses an inverse crossing block
as invertible left and right factors times the original crossing block.
Use $`S_E=I-2(E+I)^{-1}`$. Unitarity gives anti-Hermiticity and the
norm estimate follows from $`\|(E+I)^{-1}\|\le(2-\delta)^{-1}`$.

**Outward-order sharpening.** When W and Q both have parent edges before
every descendant edge, the ranks in (10) are at most **one**, including
complex native blocks. For a cut whose descendant component B starts at
b, no gate has acted in B before the crossing gate. Thus

```math
W=D_WR_W(P_W\oplus I_B),\qquad Q=D_QR_Q(P_Q\oplus I_B),
\qquad
E=D_WR_W(P_WP_Q^\dagger\oplus I_B)R_Q^\dagger D_Q^\dagger.
```

The D matrices are cut-block-diagonal. The middle product fixes every
basis label in $`B\setminus\{b\}`$, so its AB block has only the b
column and its BA block only the b row. Outer factors preserve rank one;
the preceding inverse argument gives the same bound for $`S_E`$.
This is specific to $`WQ^\dagger`$, not the left residual or arbitrary
reorderings.

For $`\delta\le1/(8\sqrt N)`$, Cauchy--Schwarz also gives

```math
\|S_E\|\le\frac1{16\sqrt N-1},\qquad
\|E\|_1,\|E\|_\infty\le1+\sqrt N\delta\le9/8. \tag{11}
```

These matrices are computable from certified inputs and the finite Q.
A convergent Neumann series around $`2I`$ computes $`(E+I)^{-1}`$;
crossing-gate factors give the ranks without nonzero-minor selection.
For anti-Hermitian S,T, Cayley reconstruction is globally stable:

```math
\mathcal C(S)=(I+S)(I-S)^{-1},\qquad
\|\mathcal C(S)-\mathcal C(T)\|\le2\|S-T\|. \tag{12}
```

This follows from the resolvent identity and both inverse norms at most
one. The reduction $`W=EQ`$ is exact: approximating E and composing on
the right with the exactly returned native Q adds no coarse error again.
The ranks concern logical matrices, not the dirty-valued errors of a
later approximate physical circuit. No iterative physical closure follows.

## 5. Shared factors, contraction slack, and the completion gap

Split labels by binary prefixes. Decompose each directed sibling block
of E into balanced rank-one singular-value factors and treat the diagonal
separately. Every entry belongs to exactly one block or the diagonal:

```math
E=AB^\dagger,\qquad A,B\in\mathbb C^{N\times m},\qquad m\le3N-2.
\tag{13}
```

There are N diagonal terms and two directed terms per split. Zero terms
can be padded; no nonzero-rank decision is required. A diagonal entry z
uses $`z/\sqrt{|z|}`$ and $`\sqrt{|z|}`$, both defined as zero at z=0.
Its Gram contribution is $`|z|\le1`$. At each hierarchy level, factor
supports are disjoint and every balanced off-diagonal Gram block has norm
at most $`\delta`$. Therefore

```math
AA^\dagger,BB^\dagger\preceq(1+n\delta)I,
\qquad n\delta\le3/(16\sqrt2)<1/4. \tag{14}
```

Write $`F=A/\sqrt2`$, $`G=B/\sqrt2`$, $`P_F=FF^\dagger`$,
$`P_G=GG^\dagger`$. Then $`FG^\dagger=E/2`$ and
$`P_F,P_G\preceq5I/8`$. The unchanged separate completion

```math
\begin{pmatrix}F^\dagger\\(I-P_F)^{1/2}\\0\end{pmatrix},\qquad
\begin{pmatrix}G^\dagger\\0\\(I-P_G)^{1/2}\end{pmatrix}
```

reserves $`5N-2`$ modes. For the arbitrary-order rank-two variant,
the shared count is $`4N-4`$ (N diagonal, N bottom-level, at most
$`2N-4`$ higher terms), and the separate completion uses $`6N-4`$.
Those counts describe these specific constructions, not general dilation
lower bounds. The following construction changes the leading realization.

### 5.1 Regularized joint completion on four flag blocks

Set $`P_F^+=P_F+I/8`$, $`P_G^+=P_G+I/8`$ and

```math
K=\begin{pmatrix}P_F^+&E/2\\E^\dagger/2&P_G^+\end{pmatrix}
=\begin{pmatrix}F\\G\end{pmatrix}
 \begin{pmatrix}F\\G\end{pmatrix}^{\dagger}+I_{2N}/8.
\tag{15}
```

The bounds are dimension-uniform:

```math
I/8\preceq K\preceq5I/4,\qquad
I/4\preceq I-P_F^+,I-P_G^+\preceq7I/8. \tag{16}
```

For the upper bound, the diagonal block norm is at most 3/4 and the cross
norm is 1/2. Choose upper Cholesky factors
$`R^\dagger R=K`$, $`C_F^\dagger C_F=I-P_F^+`$ and similarly for G.
Let $`R=[R_F\ R_G]`$ in the original two N-column groups. Then

```math
S_F=\begin{pmatrix}R_F\\C_F\\0\end{pmatrix},\qquad
S_G=\begin{pmatrix}R_G\\0\\C_G\end{pmatrix},\qquad
S_F^\dagger S_F=S_G^\dagger S_G=I,\quad S_F^\dagger S_G=E/2.
\tag{17}
```

Both sources have exactly $`4N=2^{n+2}`$ rows: these are modes on n
logical qubits and two flags, not 4N qubits. Every Cholesky pivot has
strict spectral slack. Finite-accuracy factor data can be obtained by
certified search with the Gram slack; once these data are chosen, no exact
rank or pivot-zero test is required. This is effective matrix construction,
not a quantum source-preparation algorithm.

Let $`J_{ab}`$ embed into flag block ab. Since $`J_{11}^\dagger S_F=0`$
and $`J_{10}^\dagger S_G=0`$, the full matrices

```math
H_F=I-(S_F-J_{11})(S_F-J_{11})^\dagger,\qquad
H_G=I-(S_G-J_{10})(S_G-J_{10})^\dagger
```

are Hermitian involutions taking $`J_{11}`$ to $`S_F`$ and
$`J_{10}`$ to $`S_G`$. Indeed each difference has Gram $`2I`$.
With fixed flag flips $`P_{11},P_{10}`$, put
$`A_E=H_FP_{11}`$, $`B_E=H_GP_{10}`$. Their literal overlap is

```math
J_{00}^\dagger A_E^\dagger B_EJ_{00}=E/2. \tag{18}
```

Neither reflection is asserted to be a native elementary gate.

### 5.2 Tree-cut ranks and the ordering of Cholesky

The canonical Gram matrices retain rank-one original tree cuts. Fix a
descendant interval B. Dyadic factor supports are disjoint, contained, or
containing B. Only containing supports contribute across the cut. The
restrictions of their F columns to B lie in the same one-dimensional
column space of $`E_{BA}`$; the corresponding G columns lie in the
column space of $`E_{AB}^\dagger`$. Thus the two rows of the joint cut
$`[(P_F)_{BA},E_{BA}/2]`$ and
$`[E_{AB}^\dagger/2,(P_G)_{BA}]`$ each have rank at most one, giving
joint cut rank at most two.

Order vertices in depth-first postorder, grouping the two coordinates
of K at each vertex. Every descendant component is a contiguous interval.
For $`M=R^\dagger R`$ in an order $`A_-,B,A_+`$,

```math
M_{B,A_-}=R_{A_-,B}^\dagger R_{A_-,A_-},\qquad
R_{B,B}^\dagger R_{B,A_+}
=M_{B,A_+}-R_{A_-,B}^\dagger R_{A_-,A_+}.
```

The column spaces on the right lie in the original crossing range.
Invertible diagonal factors show that upper Cholesky preserves its
at-most-r rank in both directed crossing blocks. The other triangular
blocks vanish. Restoring the original coordinates is just an explicit
permutation. Accordingly R has cut rank at most two; $`C_F,C_G`$ at
most one; the stacked sources at most three. A source-reflection cross
block has ranks at most six from $`SS^\dagger`$, three from each of
$`SJ^\dagger,JS^\dagger`$, and zero from $`JJ^\dagger`$. The full
unitaries therefore have cut rank at most twelve. These are ideal matrix
ranks, not ranks of later physical approximation errors or native costs.

### 5.3 Why merging unchanged failure banks cannot work

Keep $`F^\dagger,G^\dagger`$ fixed and append X,Y in the same r-mode
bank. Isometry requires $`X^\dagger X=I-P_F`$,
$`Y^\dagger Y=I-P_G`$. Preserving overlap requires $`X^\dagger Y=0`$.
Both ranges have dimension N, so r must be at least 2N.
More robustly, if the two deficits exceed $`d_FI,d_GI`$ and r<2N,

```math
\|X^\dagger Y\|\ge\sqrt{d_Fd_G}. \tag{19}
```

Their ranges intersect. Polar decompositions express X,Y as isometries
times positive factors with smallest singular values at least
$`\sqrt{d_F},\sqrt{d_G}`$; the two range isometries have overlap norm
one. Multiplying by their invertible positive factors gives (19). The
present lower bound is 3/8. Thus tiny approximate error does not rescue
this fixed-leading-block proposal. Equation (17) escapes by changing those
self-Grams while retaining the cross-Gram.

### 5.4 Principal square roots need not retain rank one

A concrete 16-label outward residual is obtained by changing only the
root angle to $`2\arctan u`$. The unchanged first-eight-label subtree
has root column

```math
w=(1/2,0,1/(2\sqrt2),1/(2\sqrt2),1/(2\sqrt2),1/(2\sqrt2),1/2,0)^T.
```

Its three upper rotations are $`\pi/4`$ and its four bottom rotations
are $`0,\pi/4,\pi/4,0`$. With the other root $`e_8`$,

```math
E=I-d(ww^\dagger+e_8e_8^\dagger)+s_u(e_8w^\dagger-we_8^\dagger),
\quad s_u=\frac{2u}{1+u^2},\quad d=\frac{2u^2}{1+u^2}.
```

For its canonical Gram, $`I-P_F`$ restricted to
$`e_0,(e_2+e_3)/\sqrt2,(e_4+e_5)/\sqrt2,e_6`$ is

```math
D(u)=\tfrac12[I-2uJ/(1+u^2)-2u^2M/(1+u^2)],\quad
J=\tfrac14\mathbf1\mathbf1^T,\quad
M=\tfrac14\begin{pmatrix}1&1&0&0\\1&2&0&0\\0&0&2&1\\0&0&1&1\end{pmatrix}.
```

Its 2-by-2 crossing block has rank one. For that block C of the principal
root, the convergent series around the identity gives

```math
\det C(u)=-u^6/2048+O(u^7),\qquad
\sigma_2(C(u))=\sqrt2\,u^5/1024+O(u^6). \tag{20}
```

These coefficients follow recursively by writing
$`2D(u)=\sum H_ku^k`$, $`\sqrt{2D(u)}=\sum R_ku^k`$,
$`R_0=I`$, and
$`R_k=(H_k-\sum_{j=1}^{k-1}R_jR_{k-j})/2`$.
The leading singular value is $`u/(2\sqrt2)+O(u^2)`$; its product
with the second is the absolute determinant. Exact rational series tests
check (20). Embedding in the first four decisions and setting later
angles zero preserves the witness. At $`u=1/(64\sqrt N)`$ it satisfies
$`\|E-I\|<1/(8\sqrt N)`$, but the lost singular value is
$`\Theta(N^{-5/2})`$, not endpoint-small. A rank-one approximation to
this particular root therefore fails uniformly. This is a legal residual
family, not a claim about the output of a fixed coarse-search algorithm.
Other square-root factors, triangular completions and native compilers
are not excluded.

## 6. Sufficient native source theorem, still unproved

The ideal sources in (17)--(18) settle this matrix-completion question.
They do not give the native overlap word required at
$`\varepsilon=2^{-N}`$:

```math
T(\widehat U)=O(N),\quad b\le N+n+7,\qquad
\|2J_2^\dagger\widehat UJ_2-E\otimes I_b\|\le\varepsilon/4.
\tag{21}
```

The cost includes both source preparations, completions, controls,
lookups and actual inverses. The two flags are occupied source coordinates,
not extra freshly initialized helpers inside another compiler. The norm
covers all logical and dirty inputs and arbitrary references.

If (21) holds, the existing
[normalization-two lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
applies to the actual unitary word:

```math
\mathcal R=I-2J_2J_2^\dagger,\qquad
\widehat{\mathcal A}=-\widehat U\mathcal R\widehat U^\dagger\mathcal R\widehat U,
\qquad
\|\widehat{\mathcal A}J_2-J_2(E\otimes I_b)\|\le\varepsilon.
```

The reflection is an X-conjugated CZ and the minus a literal Clifford
scalar. Two forward calls and one actual inverse retain linear total
cost; composing with (8) then gives W with the same error. No construction
establishing (21) is proved here. Matrix evaluation, bounded rank, slack,
and computable Cholesky entries are not native source operations.

### Joint words and freedom in the rejected space

Separate full-space source approximation is sufficient, not necessary.
If $`\widehat A=C_AA`$, $`\widehat B=C_BB`$, then exactly

```math
\|\widehat A^\dagger\widehat B-A^\dagger B\|
=\|C_A^\dagger C_B-I\|=\|C_A-C_B\|.
```

A common actual left gate word cancels; a cached classical table does not.
For unitaries $`D_oJ=J=D_iJ`$, let $`U'=D_oUD_i`$. They commute with
$`R=I-2JJ^\dagger`$, so
$`-U'R(U')^\dagger RU'=D_o(-URU^\dagger RU)D_i`$.
The accepted block and full amplified error to $`JE`$ are unchanged.
Thus a cheaper rejected-space completion is eligible; it need not reproduce
the Cholesky extension. These identities price no new gates by themselves.

Principal angles alone do not determine the labelled target: for
$`S=J_0`$ and $`T=(J_0/2+\sqrt3J_1/2)E`$, both range projectors
are independent of E although $`S^\dagger T=E/2`$. The abstract
$`2N`$-mode dilation containing E and $`E^\dagger`$ explicitly is
likewise not a compiler for them.

## 7. Verification scope

The original [matrix routines](../../compiler_robust_hopf/original_tree_residual.py)
and [tests](../../tests/test_original_tree_residual.py) retain the determinant,
signed-path, conditioning and residual checks. The
[joint-completion diagnostics](../../compiler_robust_hopf/source_completions.py)
and [ordered-bank tests](../../tests/test_ordered_block_banks.py) add the
regularized isometries, literal full-unitary extensions, spectral/cut bounds,
rejected-space identities and exact rational root-series witness. Floating
checks are finite diagnostics, not exact-real equality oracles, emitted fine
sources, or a new endpoint bound.
