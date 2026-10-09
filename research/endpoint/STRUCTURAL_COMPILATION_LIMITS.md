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

This is an **exact neighborhood obstruction**. It supplies no quantified
approximation lower bound at $`2^{-N}`$. Additional exterior basis
changes, four or more diagonal layers, parameter-dependent mixing, and
workspace encodings require separate arguments. In particular, it does
not rule out a different whole-operator compiler.

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

## 7. Verification and scope

The coarse circuit uses the retained full-space borrowed-MCX identity,
not a zero promise on a logical helper. The
[structural checks](../../tests/test_endpoint_structural_limits.py) exercise
the complete-matrix and native-gate identities at bounded sizes. These
finite checks do not establish the asymptotic statements; the propositions
above give their analytical proofs.

These results constrain specific representations. The signed residual
approach still has the separate native precision ledger documented in
[residual assembly](RESIDUAL_ASSEMBLY.md), and no argument here makes the
unrestricted ratio $`\tau_F^\star(n)/N`$ diverge.
