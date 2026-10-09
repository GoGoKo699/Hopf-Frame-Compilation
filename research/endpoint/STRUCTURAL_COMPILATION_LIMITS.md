# Structural limits of complete-frame compilation interfaces

**Proved interface restrictions.** This chapter collects the strongest
whole-operator restrictions established in the endpoint investigation. It
uses the complete real frame and literal initialized-isometry contract in
[the Hopf interface](../../docs/HOPF_INTERFACE.md), with $`N=2^n`$ and
$`R_y(\alpha)=e^{-i\alpha Y}`$. None of these results is an unrestricted
T-count lower bound or a resolution of the
[linear-count question](../../docs/OPEN_PROBLEM.md).

The distinctions matter: an entrywise loader loses useful interference;
a fixed diagonal factorization restricts the available algebras; an
original-angle query interface restricts classical reparameterization.
The allowed native compiler need not use any of these interfaces.

## 1. The sharp entrywise row-column normalization

Write $`|A|`$ for entrywise absolute value, not
$`\sqrt{A^\dagger A}`$. Consider the following interface for a finite
matrix A, allowing several local paths k for each entry:

```math
\frac{A_{ij}}{\alpha}=\sum_k a_{ijk}\overline{b_{ijk}},\qquad
\sum_{j,k}|a_{ijk}|^2\le1,\qquad
\sum_{i,k}|b_{ijk}|^2\le1,\qquad \alpha>0.
```

**Proposition 1.** Every such representation has
$`\alpha\ge\|\,|A|\,\|_{\rm op}`$. For nonzero A this bound is attained
in the abstract interface, using one path per entry.

For nonnegative unit vectors x and y, triangle inequality and
Cauchy--Schwarz give

```math
\begin{aligned}
\frac{x^{\mathsf T}|A|y}{\alpha}
&\le\sum_{i,j,k}x_i y_j|a_{ijk}b_{ijk}|\\
&\le\left(\sum_{i,j,k}x_i^2|a_{ijk}|^2\right)^{1/2}
     \left(\sum_{i,j,k}y_j^2|b_{ijk}|^2\right)^{1/2}\le1.
\end{aligned}
```

A nonnegative matrix has nonnegative maximizing singular vectors, proving
the lower bound. To see sharpness, let $`M=|A|`$ and
$`\alpha=\|M\|_{\rm op}`$. On each nonzero connected component of the
bipartite support graph, choose positive singular vectors satisfying
$`My=\alpha_c x`$ and $`M^{\mathsf T}x=\alpha_c y`$, with
$`\alpha_c\le\alpha`$. For nonzero entries in that component set

```math
a_{ij}=\frac{A_{ij}}{|A_{ij}|}
       \sqrt{\frac{M_{ij}y_j}{\alpha x_i}},\qquad
b_{ij}=\sqrt{\frac{M_{ij}x_i}{\alpha y_j}}.
```

Their product is $`A_{ij}/\alpha`$, and the corresponding row and column
squared norms are $`\alpha_c/\alpha\le1`$. Set all remaining entries to
zero. This also handles reducible support and isolated zero rows or columns.

The sharpness statement does **not** price state preparation, selectors,
or their return. The lower bound applies when paths belonging to different
pairs $`(i,j)`$ remain separate as above. It is not a restriction on every
block encoding or every dilation of A.

## 2. Complete frames attain the growing normalization

**Proposition 2.** Every complete real Hopf frame satisfies

```math
\|\,|W|\,\|_{\rm op}\le\sqrt{n+1}.
```

Equality holds for the balanced frame $`\mathsf H_n`$, whose angles are
all $`\pi/4`$. Here $`\mathsf H_n`$ denotes the adaptive Haar frame;
H below denotes the ordinary one-qubit Hadamard.

Each marker column is supported on its node's subtree. Therefore a row
has nonzero entries only in the root column and its n ancestral marker
columns, including when some angles vanish. Since W is orthogonal, for
every vector z,

```math
\|\,|W|z\|_2^2
\le(n+1)\sum_{i,j}|W_{ij}|^2|z_j|^2
=(n+1)\|z\|_2^2.
```

For $`\mathsf H_n`$, the root column is uniform, and each marker column
has constant absolute value $`|S|^{-1/2}`$ on its subtree S. In
$`|\mathsf H_n||\mathsf H_n|^{\mathsf T}`$, the root column contributes
one to every row sum. The marker columns at each depth also contribute
one, since their supports partition the rows. Thus all row sums are
$`n+1`$, attaining the bound.

Unitarity alone therefore cannot justify constant normalization after
passing to this entrywise loading interface.

### Linear total entry mass does not give constant normalization

The same subtree supports give the sharp global bound

```math
\sum_{i,j}|W_{ij}|
\le\sqrt N+\sum_{d=0}^{n-1}2^d\sqrt{2^{n-d}}
=(1+\sqrt2)N-\sqrt{2N}.
```

Each unit column has entry sum at most the square root of its support
size, and the balanced frame attains equality column by column. This
linear total mass coexists with a growing worst-case row sum. Set
$`\sin\theta_{d,0}=1/\sqrt{d+2}`$ and
$`\cos\theta_{d,0}=\sqrt{(d+1)/(d+2)}`$ on the zero path. The root
entry and all n marker entries in row zero then have magnitude
$`1/\sqrt{n+1}`$, by telescoping the later cosine factors.

Consequently, any signed-monomial decomposition
$`W=\sum_r\alpha_r P_r`$ has
$`\sum_r|\alpha_r|\ge\sqrt{n+1}`$ on this family: each monomial
unitary has row entry sum one. This restriction concerns that LCU
interface; the total entry-mass bound alone supplies neither coherent
loading nor a constant-normalization block encoding.

## 3. Arbitrarily small residuals with an inexpensive native coarse frame

Let C have angle $`\pi/4`$ at depths $`0,\ldots,n-2`$ and angle zero
at depth $`n-1`$. Let W change every final-depth angle to the same
$`\delta\in(0,\pi/4)`$. In original computational order,
$`W=(I_{N/2}\otimes R_y(\delta))C`$. In even-index/odd-index order,
a passive relabeling used only in this calculation,

```math
C=\mathrm{diag}(\mathsf H_{n-1},I),\qquad
W=\begin{pmatrix}\cos\delta\,I&-\sin\delta\,I\\
                    \sin\delta\,I&\cos\delta\,I\end{pmatrix}C.
```

Putting $`H_* =\mathsf H_{n-1}`$, $`c=\cos\delta`$ and
$`s=\sin\delta`$ gives the complete residual

```math
E=C^\dagger W-I
=\begin{pmatrix}(c-1)I&-sH_*^{\mathsf T}\\
                 sH_*&(c-1)I\end{pmatrix}.
```

**Proposition 3.** This pair obeys

```math
\begin{aligned}
\|E\|_{\rm op}&=2\sin(\delta/2),\\
\|\,|E|\,\|_{\rm op}&=1-\cos\delta+\sin\delta\sqrt n,\\
\frac{\|\,|E|\,\|_{\rm op}}{\|E\|_{\rm op}}
&=\sin(\delta/2)+\cos(\delta/2)\sqrt n
\longrightarrow\sqrt n.
\end{aligned}
```

Indeed $`C^\dagger W`$ is an orthogonal conjugate of a direct sum of
$`R_y(\delta)`$ blocks. The absolute-value matrix has diagonal blocks
$`(1-c)I`$ and off-diagonal blocks $`s|H_*|`$ and its transpose. Its norm
is $`1-c+s\|\,|H_*|\,\|=1-c+s\sqrt n`$. Restricting the construction
to any subtree gives discrepancy at most $`2\sin(\delta/2)`$ there as
well. The same residual identity permits unequal final angles in C if W
adds delta to each; the zero-final-angle choice supplies the native
resource guarantee below.

### Exact coarse circuit, including the reused logical helper

The chosen C has an exact measurement-free circuit with **one clean flag,
no supplied dirty ancillas**, and T-count, elementary Clifford count, and
sequential depth all $`O(n^3)=O(N)`$.

At depth d, the angle is independent of the prefix. The whole layer is
therefore $`R_y(\pi/4)=HZ`$ on its logical target t, conditioned only on
the $`s=n-d-1`$ suffix bits being zero. Define the native word

```math
B=SHTHS^\dagger=e^{i\pi/8}R_y(\pi/8),\qquad BZB^\dagger=H.
```

Consequently, with rightmost factors acting first,

```math
\mathrm{c}H=(I\otimes B)\mathrm{CZ}(I\otimes B^\dagger),
\qquad
\mathrm{c}(HZ)=\mathrm{c}H\mathrm{CZ}.
```

The scalar phases of B and its actual inverse cancel on both control
sectors. This is a literal two-T implementation of controlled
$`R_y(\pi/4)`$; each CZ is $`H_t\mathrm{CNOT}H_t`$, and every
remaining gate is elementary Clifford.

Compute the suffix-zero predicate into the clean flag using the exact
[borrowed-MCX recursion](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector),
borrowing t as its arbitrary helper. Apply the controlled rotation, then
reverse the predicate computation. The predicate macro equals the desired
toggle tensored with identity on t **on its full input space**. Its inverse
therefore remains valid after the middle rotation has changed t. The suffix
is unchanged, so the flag returns exactly to zero, including for arbitrary
logical inputs and references. The flag is reused only after exact return.

For s literals the recursion uses $`O(s^2)`$ exact Toffolis and paired X
gates for negative controls. Each exact Toffoli costs seven T gates and
a constant number of elementary Cliffords. Compute, controlled rotation,
and uncompute thus cost $`O(s^2+1)`$ of each gate type. Summing over
$`s=1,\ldots,n-1`$ gives the stated $`O(n^3)`$ bounds, at physical width
$`n+1`$. The second supplied clean qubit and every supplied dirty qubit
remain untouched.

### The witness targets themselves have linear T-count

These W are an **easy promised subfamily**, not hard instances for the
compiler. After exact C, only one physical final-qubit rotation remains.
The retained [one-qubit synthesis input](../../docs/FAULT_TOLERANT_COMPILER.md#primary-references)
approximates that determinant-one rotation with literal error
$`\eta=2^{-N}`$ using $`O(N)`$ T gates and Cliffords. Composition with C
gives the **same circuit** T-count, Clifford count, and sequential depth
$`O(N)`$, one exactly returned clean flag, no dirty use, and full
initialized-isometry error at most eta. Certified evaluation and
strict-slack finite-word search give terminating gate generation; no
uniform evaluator running-time bound is asserted.

Thus Proposition 3 exposes an unnecessary $`\sqrt n`$ penalty in the
entrywise interface, even for inexpensive native targets and arbitrarily
accurate inexpensive coarse frames. It does not contradict the signed
[weighted-transport bounds](ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
or [residual assembly](RESIDUAL_ASSEMBLY.md), which retain interference
before taking absolute values.

## 4. Three diagonal layers fail even with discontinuous phase choices

**Proposition 4.** For $`n\ge2`$, no fixed unitaries $`C_1,C_2`$ can
give, throughout a neighborhood of the zero-angle tuple,

```math
W(\theta)=D_3(\theta)C_2D_2(\theta)C_1D_1(\theta),
```

with all three D matrices computational-basis diagonal unitaries.
The conclusion does not assume continuous, differentiable, computable,
or uniquely chosen diagonal factors. The mixing unitaries need not be
Clifford.

At zero, W is identity. Its coordinate derivatives are the oriented
edge generators of the tree joining each $`v>0`$ to
$`v-\mathrm{lowbit}(v)`$. For a direction x with every
$`x_e>0`$,

```math
W(tx)=I+tA(x)+o(t),\qquad
A(x)=\sum_{e=\{u,v\}}x_e
       (|v\rangle\langle u|-|u\rangle\langle v|).
```

Suppose the factorization exists. Along a sequence $`t_j\downarrow0`$,
compactness allows all three chosen diagonal factors to converge to
$`D_{10},D_{20},D_{30}`$. The limit equation implies
$`Z=C_2D_{20}C_1`$ is diagonal. Define

```math
\begin{aligned}
M_j&=C_2[D_2(t_jx)D_{20}^\dagger]C_2^\dagger
    =[C_2D_2(t_jx)C_1]Z^\dagger,\\
\mathcal B&=C_2\{\text{complex diagonal matrices}\}C_2^\dagger.
\end{aligned}
```

Every $`M_j\in\mathcal B`$, and its entrywise magnitudes equal those
of $`W(t_jx)`$. The convergent exterior phases show that
$`\mathrm{offdiag}(M_j)/t_j`$ has a limit supported exactly on
the tree edges, with both oriented magnitudes equal to $`x_e`$.
The linear space $`\mathrm{offdiag}(\mathcal B)`$ is closed.
Hence some $`B_x\in\mathcal B`$ has precisely these off-diagonal
entries. No convergence of its diagonal difference quotients is needed.

Repeat for a different strictly positive direction y. Although the phase
limits may differ, $`B_x,B_y`$ lie in the **same fixed abelian algebra**
$`\mathcal B`$. For adjacent edges $`e=\{u,v\}`$ and
$`f=\{v,w\}`$, the unique common neighbor of u and w is v. Since u,w
are not adjacent, the $`(u,w)`$ commutator entry is

```math
0=(B_x)_{uv}(B_y)_{vw}-(B_y)_{uv}(B_x)_{vw}.
```

Taking magnitudes forces $`x_e y_f=y_e x_f`$. Connectivity of the tree's
line graph forces all ratios $`y_e/x_e`$ equal, contradicted by choosing
nonproportional positive directions.

This is an **exact neighborhood obstruction** for fixed arbitrary mixers.
The next result gives a constant approximation bound and allows
target-dependent mixers and arbitrary work, provided the two mixers are
flat Cliffords. Neither statement covers four or more diagonal layers
or a different whole-operator compiler.

### Three masks with two flat Cliffords: a width-independent gap

Call a physical Clifford on m qubits **flat** if all its matrix entries
have modulus $`Q^{-1/2}`$, where $`Q=2^m`$.

**Theorem.** For every logical $`n\ge3`$ there is a prescribed real Hopf
frame W such that every physical word

```math
V=D_3C_2D_2C_1D_1
```

with arbitrary computational unitary diagonals and flat Cliffords
$`C_1,C_2`$ satisfies

```math
\|VJ-J(W\otimes I_b)\|\ge\frac16.
```

Here J appends any number of clean zero qubits and b arbitrary dirty
qubits are allowed. Both Cliffords and all masks may depend on the
target. Their action on occupied-clean sectors is unrestricted.

**Three-point row lemma.** If a unit row f is supported on three distinct
binary labels, with modulus $`1/\sqrt3`$ at each, then for any phase
diagonal D and flat Clifford C,

```math
\inf_{h:\ |h_x|=Q^{-1/2}}\|fDC-h\|_2\ge\frac16.
```

A flat Clifford has the form

```math
C_{xy}=Q^{-1/2}r(x)c(y)(-1)^{x^{\mathsf T}Ay},
```

where r,c are unit phases and A is invertible over $`\mathbb F_2`$.
Indeed, conjugating the diagonal Paulis by C gives Pauli translation
relations between its columns. The translation map has trivial kernel:
otherwise a nontrivial diagonal Pauli would relate two diagonal sign
patterns through a matrix with no zero entries. The map is thus
invertible, and its relations give the displayed column characters.
Row and column phases and a binary relabeling reduce the modulus
calculation to the Walsh transform.

Absorb r and D into the three phases $`\alpha,\beta,\gamma`$ of f.
The three pairwise XOR differences are distinct and nonzero. For
$`p=|fDC|^2`$, its probability Fourier coefficients at those differences
are the three numbers $`(2/3)\cos(\alpha-\beta)`$ and their analogues.
At least one has modulus at least $`1/3`$, because

```math
\sum_{\rm pairs}\cos^2(\alpha-\beta)
=\frac{3+|e^{2i\alpha}+e^{2i\beta}+e^{2i\gamma}|^2}{4}
\ge\frac34.
```

Writing $`u_Q`$ for the uniform probability vector, this gives
$`\|p-u_Q\|_1\ge1/3`$. For every flat unit row h,
Cauchy--Schwarz gives

```math
\bigl\||fDC|^2-|h|^2\bigr\|_1\le2\|fDC-h\|_2,
```

proving the lemma in every physical dimension.

**Full-isometry witness.** Set all angles to zero except the edge
$`(0,2)`$ at angle $`\pi/4`$ and the final-depth edge $`(0,1)`$
with sine $`1/\sqrt3`$ and cosine $`\sqrt{2/3}`$. The prescribed
two-rotation product has zeroth row

```math
W_{0,:}=\frac{e_0^{\mathsf T}-e_1^{\mathsf T}-e_2^{\mathsf T}}{\sqrt3}.
```

Put $`U=W\otimes I_b`$ and suppose the initialized-isometry error is
epsilon. Actual unitarity gives

```math
V^\dagger J-JU^\dagger=-V^\dagger(VJ-JU)U^\dagger,
\qquad \|J^\dagger V-UJ^\dagger\|\le\epsilon.
```

Choose the physical output row with logical label zero, any fixed dirty
label and clean label zero. Its desired row f has three equal-modulus
entries. Its actual V row is within epsilon of f. Right-multiplication
by $`D_1^\dagger C_1^\dagger`$ makes the actual row a row of
$`D_3C_2D_2`$, hence flat. The desired row satisfies the row lemma
with $`D=D_1^\dagger`$ and $`C=C_1^\dagger`$, proving
$`\epsilon\ge1/6`$.

This proof permits arbitrary intermediate leakage and arbitrary
reference-entangled dirty inputs. It excludes the endpoint tolerance for
this three-mask architecture, including target-dependent full-support
Cliffords. Non-flat target-dependent mixers and larger factor counts
are outside its scope.

### Four alternating Walsh masks: every differentiable identity baseline

Let $`F=H^{\otimes n}`$ be the Walsh transform. For $`n\ge3`$, the
ansatz

```math
W(\theta)=D_4(\theta)F D_3(\theta)F D_2(\theta)F D_1(\theta)F
```

cannot hold near zero with all four diagonal masks differentiable there.
This includes every cancelling non-Clifford identity baseline and the
displayed exterior Clifford. It is specific to these fixed mixers.

To classify the baselines, write $`G=\mathbb F_2^n`$ and
$`D_a U D_b=V`$ at zero, where U and V are XOR-circulant unitaries,
$`U_{xy}=u(x+y)`$ and $`V_{xy}=v(x+y)`$. Their common support S obeys

```math
a(x)b(x+s)=r(s):=v(s)/u(s),\qquad s\in S.
```

Choose $`s_0\in S`$ and
$`K=\mathrm{span}\{s+s_0:s\in S\}`$. Ratios of these equations
give $`b(y+t)/b(y)=\chi(t)`$ for $`t\in K`$. Translation composition
makes chi a character, and exponent two gives $`\chi(t)\in\{1,-1\}`$.
Extend chi to $`(-1)^{z\cdot y}`$ on G. Precisely the scalar solutions
are therefore

```math
b(y)=(-1)^{z\cdot y}e(y),\quad e(y+t)=e(y)\ (t\in K),
\qquad a(x)=r(s_0)/b(x+s_0),
```

with $`r(s)=r(s_0)\chi(s+s_0)`$ on S. Thus all additional phase
freedom lives on $`G/K`$, and U,V are supported on $`s_0+K`$.

Let $`\mathcal Z`$ denote the diagonal Hermitian algebra and
$`\mathcal X=F\mathcal ZF`$. The right-logarithmic tangent is contained
in the complexification of

```math
\mathcal Z+\mathcal X+\mathrm{Ad}_{D_a}(\mathcal X)
+\mathrm{Ad}_{V}(\mathcal Z).
```

The following three cases exhaust the baseline classification.

- If $`0\lt\dim K\lt n`$, choose a coordinate vector $`e_j\notin K`$.
  At displacement $`e_j`$, the first three summands have K-invariant
  entry profiles, because $`a(x)/a(x+e_j)`$ is K-invariant. The final
  summand has support only on displacements in K. A single tree-edge
  generator at $`e_j`$ has profile support $`\{p,p+e_j\}`$; invariance
  would force $`K\subseteq\mathrm{span}\{e_j\}`$, a contradiction.
- If $`K=\{0\}`$, V is a scalar X Pauli, so the final summand is
  $`\mathcal Z`$. At any nonzero displacement the remaining profiles
  span at most $`1`$ and $`a(x)/a(x+s)`$. The final-depth tree edges
  give $`N/2\ge4`$ independent profiles at their common displacement.
- If $`K=G`$, $`D_a`$ is a scalar Z Pauli. In the Fourier basis the
  profiles at a nonzero displacement s span at most the constant function
  and $`\widehat v(k)\overline{\widehat v(k+s)}`$. At
  $`s=(1,\ldots,1)`$, the n tree edges $`\{0,e_j\}`$ have profiles
  proportional to the independent characters $`(-1)^{k_j}`$.

Each case excludes the full frame derivative. Differentiability is
essential: a commutator loop with phases proportional to $`\sqrt t`$
can produce an order-t direction outside this ordinary tangent space.
No assertion about singular factor selections, different mixers, or
approximate coverage follows. The
[bounded-diagonal reduction](BOUNDED_DIAGONAL_FACTORIZATION.md) allows
such alternatives without requiring a differentiable extraction map.

### Global obstruction for alternating balanced-Haar masks

Let $`Q=Q_n`$ be the balanced frame, with all angles $`\pi/4`$,
and $`u=N^{-1/2}(1,\ldots,1)^{\mathsf T}`$. Its first column is u;
the other columns are normalized signed Haar wavelets on dyadic
intervals. A **Haar phase multiplier** is

```math
A=Q\,\mathrm{diag}(z)Q^\dagger,\qquad |z_j|=1.
```

Unlike the preceding Walsh tangent obstruction, the following result
allows arbitrary nonsmooth masks and gives a quantitative approximation
bound. Both orientations and both parities have one common witness.

**Theorem.** For $`n\ge3`$, put $`M=N/2`$, $`H_*=Q_{n-1}`$, and
choose the full Hopf frame

```math
W_*=I_M\oplus H_*.
```

Its root and left-subtree angles are zero; its right-subtree angles are
$`\pi/4`$. Consider any word

```math
P=D_0Q^{s_1}D_1\cdots Q^{s_r}D_r,
```

where all D matrices are unitary diagonal and the signs $`s_j\in\{1,-1\}`$
alternate strictly. If $`\|P-W_*\|_{\rm op}\le\varepsilon\lt1/\sqrt2`$,
then

```math
9^{\lfloor r/2\rfloor+1}
\ge(1/\sqrt2-\varepsilon)2^{(n-1)/4}.
```

Thus every such architecture requires $`r=\Omega(n)`$ at fixed error
below $`1/\sqrt2`$, including error $`2^{-N}`$. No finite menu of
these words with an absolute bound on r covers every frame. No
continuity, differentiability, measurability, or effective exact
selection of the masks is assumed.

**Haar multiplier bound.** We first prove, for complex phases and inputs,

```math
\|A\|_{4\to4}\le9,\qquad \|A^\dagger\|_{4\to4}\le9.
```

Give the leaves uniform probability and let $`\mathcal F_k`$ be the
dyadic filtration through depth k. For complex leaf data f, set

```math
f_k=\mathbb E[f\mid\mathcal F_k],\quad
d_0=f_0,\quad d_k=f_k-f_{k-1}\ (k\ge1),\quad
S_f^2=\sum_{k=0}^n|d_k|^2,\quad f^*=\max_k|f_k|.
```

The finite maximal inequality
$`\|f^*\|_4\le(4/3)\|f\|_4`$ follows directly by stopping at the
first crossing of a level $`\lambda`$. Conditional Jensen gives

```math
\lambda\Pr(f^*\ge\lambda)
\le\mathbb E[|f|\,1_{\{f^*\ge\lambda\}}].
```

Multiply by $`4\lambda^2`$, integrate over positive levels, and apply
Hölder to obtain

```math
\mathbb E[(f^*)^4]
\le\frac43\mathbb E[|f|(f^*)^3]
\le\frac43\|f\|_4\|f^*\|_4^3.
```

The pointwise square identity is

```math
|f|^2=S_f^2+2Z,\qquad
Z=\sum_{k=1}^n\mathrm{Re}(\overline{f_{k-1}}d_k).
```

The real summands defining Z are martingale differences, hence are
orthogonal in $`L_2`$. Consequently

```math
\|Z\|_2^2
\le\mathbb E\sum_k|f_{k-1}|^2|d_k|^2
\le\mathbb E[(f^*)^2S_f^2]
\le\|f^*\|_4^2\|S_f\|_4^2.
```

Writing $`a=\|f\|_4`$ and $`s=\|S_f\|_4`$, triangle inequality
in both directions now gives

```math
s^2\le a^2+\frac83as,\qquad
a^2\le s^2+\frac83as.
```

The positive root of $`x^2-(8/3)x-1`$ is 3, proving
$`a/3\le s\le3a`$. Multiplication of each Haar coefficient by a
unit phase is a martingale transform: on each dyadic parent, it
multiplies the unique corresponding difference by that phase. Thus
$`S_{Af}=S_f`$ pointwise, including the constant coefficient, and

```math
\|Af\|_4\le3\|S_{Af}\|_4=3\|S_f\|_4\le9\|f\|_4.
```

Changing from probability-normalized $`L_4`$ to vector $`\ell_4`$
does not change the operator bound. Conjugating the phases proves the
adjoint bound. In particular, a product E of m Haar multipliers and
ordinary phase diagonals satisfies

```math
\|E\|_{4\to4},\ \|E^\dagger\|_{4\to4}\le9^m.
```

**The common witness.** Let $`u_L=(u_M,0)`$ and
$`u_R=(0,u_M)`$, where $`u_M=M^{-1/2}(1,\ldots,1)^{\mathsf T}`$.
They have Euclidean norm one and $`\ell_4`$ norm $`M^{-1/4}`$.
Also $`W_*e_0=W_*^\dagger e_0=e_0`$ and
$`W_*^\dagger u_R=e_M`$. For Euclidean operator error at most
epsilon, each coordinate error on a Euclidean unit input is at most
epsilon. This supplies the following lower bounds without a
dimension-dependent norm conversion.

- If $`r=2m`$ and the word starts with Q, pairing consecutive
  $`QD Q^\dagger`$ factors gives m multipliers. Apply its adjoint
  to $`u_R`$ to get
  $`1-\varepsilon\le9^mM^{-1/4}`$.
- If $`r=2m+1`$ and the word starts with Q, write
  $`P=E QD_r`$, with m multipliers in E. Since
  $`D_re_0=\zeta e_0`$ for a unit scalar zeta,
  $`Pe_0=\zeta Eu`$. Its zero coordinate gives
  $`1-\varepsilon\le9^mN^{-1/4}`$.
- For an odd word starting with $`Q^\dagger`$, take its adjoint.
  The preceding argument applies because $`W_*^\dagger e_0=e_0`$.
- For $`r=2m`$ starting with $`Q^\dagger`$, conjugate by Q:

```math
E=QPQ^\dagger
=(QD_0Q^\dagger)D_1(QD_2Q^\dagger)\cdots
 D_{2m-1}(QD_{2m}Q^\dagger).
```

This has $`m+1`$ multipliers. Put $`T=QW_*Q^\dagger`$,
$`a=M^{-1/2}`$, and $`c=1/\sqrt2`$. The recursion
$`Q=(H_*\oplus H_*)R`$, with R the balanced rotation on
coordinates $`0,M`$, gives

```math
Q^\dagger u_R=c(e_0+e_M),\qquad
T^\dagger u_R
=\frac{1-a}{2}u_L+c e_M+
 \left(\frac12+a\left(\frac12-c\right)\right)u_R.
```

Indeed, apply $`W_*^\dagger`$ to the first displayed vector and
then R and $`H_*\oplus H_*`$, using $`(H_*^\dagger e_0)_0=a`$.
The last coefficient is nonnegative, so the $`e_M`$ coordinate is
at least c. Since $`\|E-T\|_{\rm op}\le\varepsilon`$,

```math
1/\sqrt2-\varepsilon
\le\|E^\dagger u_R\|_4\le9^{m+1}M^{-1/4}.
```

All four cases imply the theorem's common inequality. This proof retains
literal scalar phases; taking a phase quotient is unnecessary.

**Robustness and scope.** Replace each zero angle in $`W_*`$ by
$`\kappa=\delta/N`$, retaining the right-subtree angles, where
$`0\lt\delta\lt1/\sqrt2`$. The resulting frame $`W_\delta`$ has
strictly positive angles. Telescoping at most $`N-1`$ rotations and
$`\|R_y(t)-R_y(s)\|_{\rm op}\le|t-s|`$ gives
$`\|W_\delta-W_*\|_{\rm op}\lt\delta`$. Approximation of
$`W_\delta`$ to epsilon therefore obeys the same lower bound with
$`\varepsilon+\delta`$ in place of epsilon, giving the same
$`\Omega(n)`$ conclusion when $`\varepsilon+\delta\lt1/\sqrt2`$.
The obstruction is not
confined to singular state charts.

This stronger linear call bound uses strict alternation. The next theorem
covers arbitrary orientations with a weaker asymptotic bound. Exact
finite witness certificates are in
[`test_global_haar.py`](../../tests/test_global_haar.py); the
dimension-uniform bounds follow from the proofs.

### Global obstruction for unrestricted Haar orientations

**Theorem.** For every $`n\ge3`$ there is a complete Hopf frame
$`W_n`$, with all angles in $`\{0,\pi/2\}`$, such that every word

```math
P=D_0Q^{s_1}D_1\cdots Q^{s_r}D_r,\qquad s_j\in\{1,-1\},
```

with arbitrary unitary diagonal masks and
$`\|P-W_n\|_{\rm op}\le\varepsilon\lt1`$ satisfies

```math
[2\sqrt{n+1}]^{r+1}\ge(1-\varepsilon)\sqrt N.
```

Orientations need not alternate, and the masks may be discontinuous,
noncomputable selections with arbitrary cancelling baselines. The same
target works for every word in this alphabet. Thus
$`r=\Omega(n/\log n)`$ at any fixed error below one, including the
endpoint error $`2^{-N}`$.

**Proof.** Put $`A=|Q|+|Q|^{\mathsf T}`$. Its graph is connected:
the flat root column of Q makes row and column zero of A positive.
There is a strictly positive Perron vector w with

```math
Aw=\lambda w,\qquad
\lambda=\rho(A)=\|A\|_{\rm op}\le2\sqrt{n+1}.
```

For completeness, take a maximizing eigenvector of this real symmetric
nonnegative matrix, replace it by its componentwise absolute value in
the Rayleigh quotient, and use connectedness to make all components
positive. The last bound follows from Proposition 2.

In the weighted norm $`\|x\|_{\infty,w}=\max_i|x_i|/w_i`$, every
unitary diagonal is an isometry, and both Q and its adjoint have induced
norm at most lambda: their absolute row sums against w are bounded by
$`Aw`$. Hence $`\|P\|_{\infty,w\to\infty,w}\le\lambda^r`$.
Choose j minimizing $`w_j`$. The flat root column gives

```math
\lambda w_0=(Aw)_0\ge N^{-1/2}\sum_iw_i\ge\sqrt N\,w_j.
```

Choose the frame angles on the path to leaf j to be zero for a left
step and $`\pi/2`$ for a right step, and zero off that path. Each
rotation is a literal signed permutation, and each right step moves
the root amplitude with coefficient +1. The prescribed complete frame
therefore satisfies $`W_ne_0=e_j`$; all other columns remain those
of the same Hopf product. Euclidean operator error gives
$`|(Pe_0)_j|\ge1-\varepsilon`$. Comparing weighted norms now yields

```math
\lambda^r\ge(1-\varepsilon)\frac{w_0}{w_j}
\ge\frac{(1-\varepsilon)\sqrt N}{\lambda},
```

as required. The proof retains literal scalar phase. The witness has
exact certified trigonometric values; its existence invokes no equality
oracle for the supplied angles.

**Five-mask consequence.** The equation

```math
W=D_0V D_1V^\dagger D_2V D_3V^\dagger D_4,\qquad V=Q^2,
```

has eight Haar occurrences and cannot cover all dimensions even with
all five masks free. For example, at $`n=128`$ and error at most
$`1/2`$, the necessary inequality would give
$`(2\sqrt{129})^9\ge2^{63}`$, whereas its left side is less than
$`32^9=2^{45}`$. The conclusion also excludes every bounded word in powers
of Q with exponents bounded independently of n, their adjoints, and
computational phase masks.

**Fixed-alphabet extension.** Let $`U_1,\ldots,U_s`$ be a fixed
logical mixer alphabet, each available in either orientation. Include
Q as an auxiliary comparison matrix, whether or not it is used, and set

```math
A=|Q|+|Q|^{\mathsf T}
  +\sum_{k=1}^s(|U_k|+|U_k|^{\mathsf T}),\qquad
\lambda\le2\sqrt{n+1}+2\sum_{k=1}^s\|\,|U_k|\,\|_{\rm op}.
```

The same proof supplies one Hopf witness for which
$`\lambda^{r+1}\ge(1-\varepsilon)\sqrt N`$, where r counts every
actual mixer occurrence. The unused comparison matrix contributes no
occurrence; omit its duplicate term if Q is already in the alphabet.
In particular, adding s fixed permutations to Q gives
$`\lambda\le2\sqrt{n+1}+2s`$. More generally, an alphabet with total
absolute-matrix norm polynomial in n still needs
$`\Omega(n/\log n)`$ calls at fixed error.

For a menu of architectures, this norm sum must include the union of
its mixers. The result does not treat an arbitrary target-dependent
Clifford or permutation as one member of a fixed small alphabet.
The Walsh Clifford $`H^{\otimes n}`$ has absolute-matrix norm
$`\sqrt N`$ and escapes the useful polynomial bound. These are
representation restrictions, not unrestricted T-count bounds: arbitrary
finite Clifford work and growing words with jointly amortized source
cost remain within the native compiler contract.

### Balanced-Haar boundary identities

Fixed balanced frames are another available mixer. Write
$`Q_n=B R`$, where $`B=Q_{n-1}\oplus Q_{n-1}`$ and
$`R=R_{y,\mathrm{root}}(\pi/4)`$ acts on the distinguished coordinates
0 and $`N/2`$. These fixed mixers have exact $`O(n^3)=O(N)`$ T-count
by the controlled-rotation construction in Section 3, including the
uncontrolled final layer.

For a paired child word

```math
W_{\rm child}=D_0B^\dagger D_1B D_2\cdots
 B^\dagger D_{2m-1}B D_{2m},\qquad m\ge1,
```

scalar phases on each half commute with B and can be transferred between
masks. Choose them so every even mask has equal distinguished-pair
entries, and hence commutes with R. Substitution then gives
$`W_{\rm child}=R A R^\dagger`$, where A is the same mask word with
B replaced by $`Q_n`$. The frame recursion becomes
$`W_n=R A R_{y,\mathrm{root}}(\theta-\pi/4)`$.

The variable boundary has a literal diagonal identity. Let S(t) be
$`e^{it}`$ on the left half and $`e^{-it}`$ on the right; let
$`Z_r(t)`$ agree with S(t) on the distinguished pair and equal identity
elsewhere; and put $`P=S(-\pi/4)`$. Since B commutes with S(t),

```math
R_{y,\mathrm{root}}(t)
=P Q_n^\dagger S(t)Q_nP^\dagger S(t)^\dagger Z_r(t),
```

using $`R^\dagger ZR=-X`$ and $`PXP^\dagger=Y`$ on the pair.
The last two diagonal factors cancel the unwanted spectator phases.
Only valid commuting moves give

```math
W_{\rm child}R_{y,\mathrm{root}}(t)
=P S(t)^\dagger[W_{\rm child}Z_r(t)]Q_n^\dagger S(t)Q_nP^\dagger.
```

The remaining $`Z_r(t)`$ changes each distinguished coordinate relative
to its own half, so it is not a transferable half-block scalar. Absorbing
it into a final mask still appends a new $`Q_n^\dagger/Q_n`$ pair.
These all-angle, division-free identities retain the exact boundary
phases. The preceding theorems exclude bounded-mask closure using only
Q and its adjoint, in any orientation. The
[tree Cayley reduction](TREE_CAYLEY_REDUCTION.md) supplies a separate
global normal form.

## 5. Original-angle parallel coins need at least n queries

Let the original-angle oracle be

```math
S(\theta)=\bigoplus_v R_y(\theta_v),
```

padded by identity. A query may be S, its adjoint, a controlled version,
or a fixed unitary conjugate. A circuit makes k queries with arbitrary
finite **parameter-independent** interlayers and fixed input/output
encodings. Its workspace may exceed the endpoint allowance. A tensor
product of r oracle occurrences counts as r queries.

The error-one claim uses the full signed-angle family. Under a separate
canonical-angle promise $`\theta_v\in[0,\pi/2]`$, the
[canonical-interval bound](ENDPOINT_TREE_TRANSPORT.md#a-fixed-number-of-unchanged-coin-calls-does-not-give-the-frame)
applies instead; the proof below does not assert the same threshold there.

**Proposition 5.** If $`k\lt n`$, the worst-case literal initialized-isometry
error in implementing the complete frame is at least one.

Use the all-right path $`v_d=(d,2^d-1)`$. Set its angles to
$`\theta_{v_d}=x_d\pi/2`$, $`x_d\in\{-1,1\}`$, and all other
angles to zero. These certified, possibly singular inputs are admitted.
Each one-query entry is a constant or a scalar times one Boolean sign.
After k queries, every amplitude q is a Boolean polynomial of degree at
most k, even after contraction with fixed initialized work and any fixed
dirty input/output. The prescribed root-to-rightmost-leaf amplitude is

```math
w(x)=\langle1^n|W(x)|0^n\rangle
=\prod_{d=0}^{n-1}x_d=\chi(x).
```

If $`k\lt n`$, Boolean Fourier orthogonality gives
$`\mathbb E_x[\chi(x)q(x)]=0`$. Therefore

```math
1=\left|\mathbb E_x[\chi(x)(w(x)-q(x))]\right|
\le\max_x|w(x)-q(x)|
\le\sup_\theta\|Q(\theta)J-J(W(\theta)\otimes I)\|_{\rm op}.
```

The literal scalar phase is essential: a parity decision protocol with
only a measurement guarantee is a different task. This strengthens the
[unchanged-coin frequency obstruction](ENDPOINT_TREE_TRANSPORT.md#a-fixed-number-of-unchanged-coin-calls-does-not-give-the-frame)
for the full signed-angle contract. Classical reparameterization into
derived-angle oracles, parameter-dependent interlayers, and joint native
compilation escape this model. Query costs are not automatically additive
T-count lower bounds.

## 6. A quantitative limit on one fixed-basis block compression

**Proposition 6.** Let $`n\ge2`$. Fix an exact orthogonal projector P with
$`0\lt \mathrm{rank}P\lt N`$. If every frame is within epsilon of a
unitary commuting with P, then

```math
\epsilon\ge\frac1{16n(N-1)}.
```

Every individual tree-edge rotation occurs in the frame family. Thus
$`\|P-R_ePR_e^\dagger\|\le2\epsilon`$, and a word of length ell in
these rotations changes P by at most $`2\ell\epsilon`$.
Edge rotations by pi generate every even diagonal sign pattern with at
most $`N-1`$ factors, since a tree's binary incidence columns span the
even-parity vectors. Their conjugation average kills all off-diagonal
entries for $`N\ge3`$.

The coordinate tree has diameter at most $`2n`$. A coordinate permutation
uses at most $`N-1`$ arbitrary transpositions, each realized as at most
$`4n-1`$ edge transpositions. Edge rotations by $`\pi/2`$ realize these
with signs, which do not affect diagonal entries. Averaging next over
permutations produces $`(\mathrm{rank}P/N)I`$. Every word in this
combined average has length at most $`4n(N-1)`$, so

```math
\frac12\le
\left\|P-\frac{\mathrm{rank}P}{N}I\right\|
\le8n(N-1)\epsilon.
```

For $`n\ge4`$ this excludes endpoint-accurate compression into any one
proper fixed block decomposition, including a direct sum of two-dimensional
blocks after a fixed arbitrary basis change. It does not exclude products
of noncommuting block operators or parameter-dependent encodings. The
displayed constant alone does not exclude the $`n=3`$ endpoint.

### One U(2) bank fails even with arbitrary fixed encoders and workspace

A stronger constant bound holds for blocks of dimension at most two.
Let A,B be fixed unitaries on any auxiliary width, J the initialized-work
isometry, and every $`M_X`$ belong to the **same** direct-sum-of-U(2)
block algebra. If

```math
\|A M_X B J-J(X\otimes I_{\rm dirty})\|\le\epsilon
```

holds for all complete frames X, then for $`n\ge2`$,

```math
\epsilon\ge\frac{\sqrt2}{30}.
```

For any two 2-by-2 matrices the commutator has trace zero; Cayley--Hamilton
therefore makes its square scalar. The polynomial
$`p(U,V)=[[U,V]^2,U]`$ vanishes on the common block algebra. Two legal
frame instances, each rotating one of two edges sharing a vertex by
$`\pi/2`$, restrict to

```math
U=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&1\end{pmatrix},\qquad
V=\begin{pmatrix}0&0&-1\\0&1&0\\1&0&0\end{pmatrix},\qquad
p(U,V)=\begin{pmatrix}0&-2&0\\-2&0&-2\\0&2&0\end{pmatrix}.
```

The final matrix has norm $`2\sqrt2`$. Put $`E=BJ`$ and
$`R_X=M_I^\dagger M_X`$. The identity instance and unitarity imply
$`\|R_XE-E(X\otimes I_{\rm dirty})\|\le2\epsilon`$.
The six surviving signed words in p all have length five, so telescoping
each word on E gives

```math
2\sqrt2=\|E(p(U,V)\otimes I_{\rm dirty})\|
\le60\epsilon,
```

since $`p(R_U,R_V)=0`$. This uses full initialized columns and remains
valid for arbitrary dirty inputs and references; it makes no assumption
that the encoder's image is an exact invariant subspace. It excludes the
endpoint error for every $`n\ge3`$. Multiple banks and parameter-dependent
encoders remain outside its scope; Proposition 6 separately covers larger
fixed blocks without arbitrary encoding.

## 7. Low-rank defects after restricted preconditioning

### Monomial preconditioning and reflection rank

Let $`W_*`$ have zero upper angles and final-depth angles
$`\pi/4`$, so $`W_*`$ is a direct sum of $`R_y(\pi/4)`$ blocks.
If M is any monomial unitary, R is unitary, and
$`\|W_*-MR\|\le\varepsilon\lt1`$, then

```math
\mathrm{rank}(R-I)\ge
N\frac{1-\varepsilon^2/2-1/\sqrt2}{2-\varepsilon^2/2}
```

whenever the numerator is positive. M may depend arbitrarily on the
target, including its permutation and phases.

Every entry of $`W_*`$ has modulus at most $`1/\sqrt2`$, so
$`\mathrm{Re}\,\mathrm{tr}(W_*^\dagger M)\le N/\sqrt2`$.
Write $`r=\mathrm{rank}(R-I)`$. On its kernel, of dimension $`N-r`$,
the approximation implies
$`\mathrm{Re}\langle W_*x,Mx\rangle\ge1-\varepsilon^2/2`$
for each unit x. Complete an orthonormal basis and bound each remaining
real correlation below by -1. Thus

```math
N/\sqrt2\ge(N-r)(1-\varepsilon^2/2)-r,
```

proving the claim. Products satisfy
$`\mathrm{rank}(R_k\cdots R_1-I)\le\sum_j\mathrm{rank}(R_j-I)`$.
Interleaved monomials can be moved outside by conjugating these factors,
which preserves rank. Hence bounded total reflection rank cannot cover
the family after arbitrary monomial preconditioning. This does not
exclude a bounded number of reflections of rank $`\Theta(N)`$ or
$`O(N)`$ low-rank factors with a jointly linear precision cost.

### A fixed dense baseline menu

Fix any menu $`C_1,\ldots,C_s`$ of logical unitaries, with an absolute
constant $`s\ge1`$. There is a final-depth-only Hopf frame W such that

```math
\max_{a,i,j}|(W^\dagger C_a)_{ij}|\le c_s,
\qquad c_s=\cos\frac{\pi}{16s}\lt1.
```

For each row pair, choose its real rotation independently. A column of
$`C_a`$ has pair energy $`\rho_j=|u_j|^2+|v_j|^2\le1`$, with
$`\sum_j\rho_j=2`$. Only columns with $`\rho_j>c_s^2`$ can violate
the desired bound; there are at most two since $`c_s^2>2/3`$. The
squared moduli of the rotated outputs are

```math
\frac{\rho_j}{2}\ \pm\bigl(A_j\cos(2\theta)+B_j\sin(2\theta)\bigr),
\qquad \sqrt{A_j^2+B_j^2}\le\rho_j/2\le1/2.
```

Either output exceeding $`c_s^2`$ requires
$`|\cos(2\theta-\phi_j)|>2c_s^2-1`$ for a suitable phase. The bad
angles occupy at most $`(4/\pi)\arccos c_s=1/(4s)`$ of a period.
Across at most $`2s`$ relevant columns, their union has measure at most
one half. Choosing an allowed angle for each pair proves the claim,
also for complex baselines.

For every a and every monomial M this gives
$`\mathrm{Re}\,\mathrm{tr}(W^\dagger C_aM)\le c_sN`$. The same
kernel argument therefore yields

```math
\|W-C_aMR\|\le\varepsilon\lt1
\quad\Longrightarrow\quad
\mathrm{rank}(R-I)\ge
N\frac{1-\varepsilon^2/2-c_s}{2-\varepsilon^2/2}.
```

This statement has the displayed factor order and a constant-size
baseline menu. It does not cover arbitrary interleavings of dense
mixers and target-dependent masks or a growing menu.

### Fixed diagonal displacement

Let A,B be fixed diagonal matrices. If
$`\mathrm{rank}(AW-WB)\le r`$ for every complete Hopf frame W,
and A has $`K\ge2`$ distinct diagonal entries, then $`K\le2r`$.
In particular, a simple-spectrum A forces $`r\ge N/2`$.

The identity frame gives $`\mathrm{rank}(A-B)\le r`$. Set every
angle to $`\pi/2`$ to obtain a signed permutation P whose underlying
permutation is an N-cycle. Indeed, the rotation edges form a tree, and
the product of its edge transpositions, each used once, is an N-cycle.
An induction removes a leaf edge: all remaining factors fix that leaf;
conjugating the removed transposition past one prefix attaches the leaf
to the single cycle on the other vertices.
Consequently

```math
\mathrm{rank}(AP-PA)
\le\mathrm{rank}(AP-PB)+\mathrm{rank}(P(B-A))\le2r.
```

The commutator is a monomial times a diagonal, with a nonzero entry at
each change of A's labels around the cycle. A cyclic sequence with
$`K\ge2`$ labels has at least K such transitions, proving the bound.
This is a fixed diagonal Sylvester-displacement restriction.
Nondiagonal or target-dependent displacement operators are outside its
scope. The
[hierarchical rank-one cuts](BOUNDED_DIAGONAL_FACTORIZATION.md#hierarchical-rank-one-cuts)
show why highly degenerate cut projections can still give low rank.

## 8. Phase matching and terminal Clifford gauges

The [nested-projector criterion](BOUNDED_DIAGONAL_FACTORIZATION.md#4-complete-coverage-through-nested-projectors)
requires simultaneous complete-frame identities. Flattening each local
root separately does not ensure them. The following calculations locate
two restrictions on transporting such local phase corrections.

### Unequal children obstruct one simultaneous leaf gauge

For positive four-leaf angles $`a,b,c\in(0,\pi/2)`$, the root vector is

```math
(\cos a\cos b,\ \cos a\sin b,\ \sin a\cos c,\ \sin a\sin c).
```

A computational phase gauge makes the two children Walsh-flat only if
their respective phases are in quadrature. Write the gauged child
vectors as

```math
e^{i\alpha}(\cos b,i\epsilon_L\sin b),\qquad
e^{i\beta}(\cos c,i\epsilon_R\sin c),\qquad
\epsilon_L,\epsilon_R\in\{1,-1\}.
```

Their Walsh outputs have phases $`\pm\epsilon_L b`$ and
$`\pm\epsilon_R c`$. Put
$`\Delta=\epsilon_Rc-\epsilon_Lb`$ and $`\delta=\beta-\alpha`$.
The parent's two frequency pairs are flat exactly when

```math
\cos(\delta+\Delta)=\cos(\delta-\Delta)=0.
```

This requires $`\Delta\in(\pi/2)\mathbb Z`$. For
$`b=\pi/12,c=\pi/6`$ no sign choice meets that requirement.
More quantitatively, retaining exact child flatness and minimizing the
parent's largest probability deviation from $`1/4`$ over delta gives

```math
\frac{|\cos a\sin a|}{2}
\min\{|\sin\Delta|,|\cos\Delta|\}.
```

At $`a=\pi/4`$, the minimum over both signs is
$`\sin(\pi/12)/4>0.0647`$. These identities apply inside a larger
tree without dividing by its global subtree mass.

There is an exact frequency-dependent repair. Choose
$`\alpha=\beta=0`$ and insert, between the two child Walsh transforms
and the parent Hadamard, the diagonal

```math
B(b,c)=\mathrm{diag}
 (e^{-i\epsilon_Lb},e^{i\epsilon_Lb},
  i e^{-i\epsilon_Rc},i e^{i\epsilon_Rc}).
```

It cancels the child frequency phases and places the two children in
quadrature, making all four parent output moduli $`1/2`$. Repeating
this repair adds a phase stage at each height. Its direct packed-mask
cost is $`O(nN)`$, since each fine phase bank retains its precision
charge. This root-vector calculation alone supplies no complete-frame
factorization.

### Transport into one terminal Clifford gauge is rigid

Call a unitary **Clifford-diagonalizable** if it equals
$`CDC^\dagger`$ for a Clifford C and a computational unitary diagonal D.
Write a complete frame recursively as $`W=(W_L\oplus W_R)R`$, with R
the root rotation on the two child anchors. If the child words equal
$`F_s=W_sG_s`$, then their full-operator phase mismatch is exactly

```math
(F_L\oplus F_R)R=W H_G,\qquad
H_G=R^\dagger G R,\qquad G=G_L\oplus G_R.
```

**Theorem.** Let $`G_L,G_R`$ be any four-dimensional
Clifford-diagonalizable unitaries. Let R rotate anchors 0 and 4 with
cosine $`3/5`$ and sine $`4/5`$. If $`H_G`$ is
Clifford-diagonalizable, then both child gauges fix their respective
anchors with the same eigenvalue. In particular, $`[G,R]=0`$ and
$`H_G=G`$. All three Cliffords may be target-dependent and non-flat;
arbitrary eigenvalue degeneracies are allowed.

**Proof.** Every spectral projector of an m-qubit
Clifford-diagonalizable unitary has entries in
$`2^{-m}\mathbb Z[i]`$. To see this, a Clifford column is a joint
eigenvector of commuting Paulis. Their translation constraints give an
affine binary support with constant modulus $`2^{-r/2}`$, with relative
phases in $`\{1,i,-1,-i\}`$. The other columns are Pauli translates,
with the same r. A sum of their rank-one projectors therefore has
entries in $`2^{-r}\mathbb Z[i]\subseteq2^{-m}\mathbb Z[i]`$.

Fix an eigenvalue lambda of G. Its spectral projector is
$`P=P_L\oplus P_R`$, allowing either block to be zero. The corresponding
projector of $`H_G`$ is $`Q=R^\dagger P R`$ and must have Gaussian
dyadic entries. For $`j=1,2,3`$,

```math
Q_{0j}=\frac35(P_L)_{0j}.
```

Write $`(P_L)_{0j}=(u+iv)/4`$. Projector entries have modulus at most
one, so $`|u|,|v|\le4`$. Dyadicity forces 5 to divide both integers,
hence both vanish. The corresponding right-child entries also vanish.
Each anchor is now an isolated row and column of its projector, so
its diagonal entry $`p_s`$ is zero or one. But

```math
Q_{00}=\frac{9p_L+16p_R}{25}
```

is dyadic only when $`p_L=p_R`$. For every eigenvalue, both anchors
therefore belong to the same eigenspace, proving the theorem.

For child dimension $`M=2^m`$, the same conclusion holds for a
primitive Pythagorean root rotation with cosine $`a/q`$ and sine
$`b/q`$, where q is odd and $`q>M`$. Child projector numerators now
have real and imaginary parts bounded by M. Since
$`\gcd(a,q)=1`$, the same divisibility proof isolates the anchors;
the transported diagonal then forces their eigenspaces to agree.
Such angles exist, for example from
$`(a,b,q)=(k^2-1,2k,k^2+1)`$ with even $`k\ge2`$ and $`k^2+1>M`$.

### An unequal eight-mode example with a finite gap

Take the repaired profiles above with positive quadrature signs and

```math
e^{ib}=\frac{4+3i}{5},\quad
e^{ic}=\frac{12+5i}{13},\quad
e^{id}=\frac{15+8i}{17},\quad
e^{ie}=\frac{24+7i}{25}.
```

With $`C=H\otimes H`$, set
$`G_L=CB(b,c)C^\dagger`$ and $`G_R=CB(d,e)C^\dagger`$ and use the
$`3/5,4/5`$ root rotation. The eight eigenvalues of $`H_G`$ are
pairwise distinct. Every normalized eigenvector has the same multiset
of absolute amplitudes,

```math
\left\{\frac3{10},\frac12,\frac12,\frac12,\frac25,0,0,0\right\}.
```

Thus none is a stabilizer vector: its nonzero amplitudes are unequal
and its support has size five. The first eigenprojector also has the
non-dyadic entries $`Q_{00}=9/100`$ and $`Q_{01}=3/20`$.
These profiles use four unequal bottom angles of an eight-leaf frame;
its two child root angles can also be unequal. The identity for
$`F_s=W_sG_s`$ concerns their complete unitary matrices.

The minimum separation of the eight eigenvalues is
$`g=2\sqrt{13}/65`$. A stabilizer vector has uniform modulus on a
support of size 1,2,4 or 8. Even ignoring its affine-support and phase
restrictions, its overlap with a target eigenvector is at most
$`19/20`$, attained as a modulus bound by summing the four largest
amplitudes and dividing by two. The distance between their rank-one
projectors is therefore at least $`d=\sqrt{39}/20`$.

Suppose a Clifford-diagonalizable K satisfies
$`\|H_G-K\|=\epsilon`$. Since both matrices are normal, for any target
eigenvalue lambda some eigenvalue mu of K lies within epsilon: apply
$`K-\lambda I`$ to a target eigenvector and use the normal resolvent
bound. Choose a Clifford-basis eigenvector phi of K at mu. It obeys
$`\|(H_G-\mu I)\phi\|\le\epsilon`$, while its component orthogonal
to the lambda eigenvector has norm at least d. If $`\epsilon\lt g/2`$,
all other target eigenvalues are at distance at least
$`g-\epsilon`$ from mu. Hence
$`\epsilon\ge(g-\epsilon)d`$. The same resulting bound is automatic
when $`\epsilon\ge g/2`$. Consequently

```math
\|H_G-K\|\ge\frac{gd}{1+d}>0.026>2^{-8}.
```

This excludes exact or endpoint-accurate replacement of this terminal
gauge by one Clifford-diagonalizable bank. It does not exclude a joint
rewrite using several banks and cancellations in the complete word.

### The transported correction has rank at most four

For any child dimension, define

```math
K_G=H_GG^\dagger=R^\dagger G R G^\dagger,\qquad
\mathcal S=\mathrm{span}\{e_0,e_M,Ge_0,Ge_M\}.
```

Then $`\mathrm{rank}(K_G-I)\le4`$ and $`K_G`$ fixes
$`\mathcal S^\perp`$ pointwise. Indeed,

```math
H_G-G=R^\dagger G(R-I)+(R^\dagger-I)G
```

has range in $`\mathcal S`$: $`R-I`$ has range in the anchor pair,
and $`R^\dagger`$ changes a vector only by an anchor-pair vector.
Right-multiplication by $`G^\dagger`$ preserves this range bound.
Unitarity then makes $`\mathcal S^\perp`$ a fixed subspace.

The displayed commutator is an exact constructive correction, including
degenerate child spectra. Even when both children share a Clifford
skeleton, so G itself is one packed phase bank, this word calls both G
and its actual inverse. The rank bound supplies neither their amortized
precision cost nor a simultaneous block decomposition of the overlapping
ancestor correction spaces.

## 9. Verification and scope

The coarse circuit uses the retained full-space borrowed-MCX identity,
not a zero promise on a logical helper. The
[structural checks](../../tests/test_endpoint_structural_limits.py) exercise
the complete-matrix and native-gate identities at bounded sizes. These
finite checks do not establish the asymptotic statements; the propositions
above give their analytical proofs.

The [refinement-route checks](../../tests/test_endpoint_refinement_routes.py)
also exercise the flat-Clifford witness, unequal-child phase equations,
complete eight-mode gauge identity, spectral gap, and rank-four support
calculation. Gaussian-rational identities certify the eight-mode example
without relying on numerical eigenvectors.

These results constrain specific representations. The signed residual
approach still has the separate native precision ledger documented in
[residual assembly](RESIDUAL_ASSEMBLY.md), and no argument here makes the
unrestricted ratio $`\tau_F^\star(n)/N`$ diverge.
