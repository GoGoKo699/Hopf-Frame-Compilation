# A weighted transport block with one signal flag

[Tree transport](ENDPOINT_TREE_TRANSPORT.md) · [Research status](OPEN_PROBLEM.md) · [Borrowed compiler](BORROWED_WORKSPACE_COMPILER.md)

The forward weighted transport has an explicit constant-normalization
unitary dilation on one signal flag and the n-qubit logical register.
A recursive scalar correction makes its marker inputs orthogonal without
an extra depth register. Its complete rejected-space action is specified.
A conservative native implementation costs
$`T,G=O(N(L+n^3))`$, hence $`O(N^2)`$ at $`L=N`$.
It meets the endpoint's workspace and Clifford budgets, but misses the
$`O(N)`$ T-count target. The existing full-frame bound is unchanged.

This chapter completes one bounded candidate audit. It does not claim that
quadratic cost is necessary, nor that a block for one residual term already
compiles the frame. The finite checks use matrices of dimension at most 32;
the native resource bound is analytic, not an end-to-end circuit fixture.

## 1. Why defect-weighted stopping needs a rejection correction

Use the heap indexing and overlap generators from the
[tree-transport chapter](ENDPOINT_TREE_TRANSPORT.md). Put
$`M=D_h\mathcal P_W`$, with internal-node rows and logical root/marker
columns. Write $`u_0,u_1`$ for the columns of $`U_v^W`$ and
$`\epsilon_v=1-|g_v|^2`$, with zero defects at leaves. Define

```math
z_v=\sum_{b=0}^1|u_1[b]|^2\epsilon_{2v+b},
\qquad
\tau_v=\sum_{b=0}^1\overline{u_0[b]}u_1[b]\epsilon_{2v+b}.
```

The exact weighted Gram entries are

```math
\begin{aligned}
(M^\dagger M)_{\ast,\ast}&=\epsilon_1,&
(M^\dagger M)_{v,v}&=z_v,\\
(M^\dagger M)_{\ast,v}&=\overline{a_W(v)}\tau_v,&
(M^\dagger M)_{x,v}&=\overline{\beta_W(x,v)}\tau_v
\quad(x\text{ a strict ancestor of }v).
\end{aligned}
```

Incomparable marker entries vanish; the other entries follow by adjoint.
Apply subtree defect telescoping below each child of v to prove these
identities. The local overlap matrix also gives
$`\tau_v=-(\overline g_vk_v+\overline h_vd_v)`$.
Unequal child defects therefore destroy the orthogonality used by the
unweighted history unitary.

For a two-level example, take the native coarse frame $`C=I`$ and target
words $`U_1^W=U_2^W=R_y(t)`$, $`U_3^W=I`$. With
$`c=\cos t`$, $`s=\sin t`$, and columns ordered as
$`(\ast,1,2,3)`$, the weighted matrix is

```math
M=\begin{pmatrix}s&0&0&0\\cs&-s^2&0&0\\0&0&0&0\end{pmatrix}.
```

Its two nonzero Gram eigenvalues are $`s^2(1+c)`$ and $`s^2(1-c)`$.
For $`0\lt t\lt\pi/2`$, this gives
$`\|M\|>\sqrt\gamma`$, where $`\gamma=s^2(1+c^2)`$.
The independently normalized columns have inner product
$`-c/\sqrt{1+c^2}`$, tending to $`-1/\sqrt2`$ as t tends to zero.
Every subtree discrepancy is at most $`2|t|`$, so the obstruction remains
for arbitrarily accurate coarse frames. It rules out this normalization
shortcut, not a dilation at the safe larger normalization.

## 2. A scalar recursion supplies the missing orthogonality

Use $`N=2^n`$, $`n\ge1`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
Fix a dyadic constant $`0\lt\varepsilon_0\le1/64`$ and set
$`\alpha=4\varepsilon_0`$. Choose the actual coarse C using the retained
borrowed compiler at accuracy $`\varepsilon_0/2`$. For the moment assume
$`\gamma\le\varepsilon_0^2`$; Section 4 preserves this allowance during
certified preprocessing.

For the subtree at v, split its weighted transport as
$`H_v=[t_v\ \ M_v]`$, where $`t_v`$ is the root-state column and
$`M_v`$ contains that subtree's marker columns. The retained norm lemma gives
$`\|H_v\|\le2\sqrt\gamma`$ and $`\|t_v\|^2=\epsilon_v`$.
Define the nonnegative scalar

```math
\begin{aligned}
\rho_v
&=\sup_x\left(\|t_v+M_vx\|^2-\alpha^2\|x\|^2\right)\\
&=t_v^\dagger(I-M_vM_v^\dagger/\alpha^2)^{-1}t_v
\le\frac{\epsilon_v}{1-4\gamma/\alpha^2}
\le\frac43\epsilon_v.
\end{aligned}
```

This is a strictly concave quadratic maximization in x. Completing the
square gives the inverse formula; the norm bound gives its positivity and
estimate. These are standard Schur-complement operations; see
[Boyd–Vandenberghe, Appendix A.5.5](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).
Their application here supplies a tree-local dilation, not a free matrix
inverse subroutine.

The scalars can be evaluated bottom up without forming $`H_v`$.
Set leaf values to zero and define

```math
\begin{gathered}
R_v=\mathrm{diag}(\rho_{2v},\rho_{2v+1}),\qquad
A_v=u_0^\dagger R_vu_0,\quad B_v=u_1^\dagger R_vu_1,
\quad C_v=u_0^\dagger R_vu_1,\\
\beta_v=\alpha^2-B_v,\qquad
\rho_v=|h_v|^2+A_v+\frac{|C_v|^2}{\beta_v}.
\end{gathered}
```

Indeed, eliminate the independent child-marker inputs in the supremum.
The remaining scalar marker input z contributes
$`|h_v|^2+\sum_b\rho_{2v+b}|u_0[b]+u_1[b]z|^2-\alpha^2|z|^2`$.
Maximizing over z gives the recursion. Its denominators obey

```math
\beta_v\ge\alpha^2-\frac43\gamma
\ge\frac{11}{12}\alpha^2>0,
\qquad \rho_1/\alpha^2\le1/12.
```

The variational definition proves finiteness before the recursion is used;
there is no assumption that the recursion happens to remain bounded.

### Local complete unitaries

For a nonterminal internal node, order four output modes as stop, left
continuation, right continuation, rejection. If $`\rho_v>0`$, prescribe
the first two columns

```math
p_v=\frac1{\sqrt{\rho_v}}
\begin{pmatrix}
h_v\\ \sqrt{\rho_{2v}}u_0[0]\\
\sqrt{\rho_{2v+1}}u_0[1]\\-\overline{C_v}/\sqrt{\beta_v}
\end{pmatrix},
\qquad
q_v=\frac1\alpha
\begin{pmatrix}
0\\ \sqrt{\rho_{2v}}u_1[0]\\
\sqrt{\rho_{2v+1}}u_1[1]\\\sqrt{\beta_v}
\end{pmatrix}.
```

The recursion gives unit norm to both columns. Their inner product is
$`C_v/(\alpha\sqrt{\rho_v})-C_v/(\alpha\sqrt{\rho_v})=0`$.
Complete them by Gram–Schmidt on the standard basis in order, skipping
zero projections and using positive normalizing norms. Phase the last completion
column to make the determinant one. Thus all rejected input modes have a
defined unitary action as well.

If $`\rho_v=0`$, take $`p_v`$ to be the stop basis vector; the displayed
$`q_v`$ is still defined and orthogonal to it. The input continuation
amplitude on the accepted sector will be zero, so this choice does not
alter the required block. At a terminal internal node use the two-mode
unitary $`\mathrm{diag}(\omega_v,\overline{\omega_v})`$, where
$`\omega_v=h_v/|h_v|`$ for nonzero h and $`\omega_v=1`$ otherwise.
Its second output is rejection. These conventions include singular charts.

## 3. One flag supplies all unused modes

The register has $`2N`$ basis modes, implemented by n logical qubits and
one signal flag. Initially the flag-zero half contains all N logical input
amplitudes, including arbitrary superpositions of root and marker columns.
The flag-one half has zero amplitude. These are unused basis modes, not
N individually initialized work qubits.

Use the following explicit allocation algorithm.

1. Mix the root input mode and one unused mode by the real determinant-one
   rotation with first column $`(s_0,r_0)^{\mathsf T}`$,
   where $`s_0=\sqrt{\rho_1}/\alpha`$ and
   $`r_0=\sqrt{1-s_0^2}`$. Label its outputs root continuation and root rejection.
2. Process internal nodes from shallow to deep. At each nonterminal node,
   take its continuation mode, its untouched logical marker input, and two
   previously unused modes. Apply the four-mode unitary above; label its
   outputs stop, two child continuations, and rejection.
3. At each terminal node, apply the two-mode unitary to its continuation
   and untouched marker input. Retain the stop and rejection outputs.
4. Permute the N−1 stop modes into the corresponding flag-zero marker
   outputs. Send the one remaining unused mode to the flag-zero root output,
   and the N rejection modes to the flag-one half.

There are $`N/2-1`$ nonterminal internal nodes. The total number of unused
modes consumed is $`1+2(N/2-1)=N-1`$; one remains for padding.
All factors are complete unitaries, and the final map is a permutation.
For $`n=1`$ the same count holds with no four-mode factors.

Let $`x_v`$ be the marker input amplitude and $`a_v`$ the incoming
path amplitude, so $`a_1=x_\ast`$ and
$`a_{2v+b}=u_0[b]a_v+u_1[b]x_v`$. Inductively, the actual continuation
amplitude at v is $`\sqrt{\rho_v}a_v/\alpha`$. Its stop output is
$`h_va_v/\alpha`$, while the two continuation outputs preserve the
induction. If $`\rho_v=0`$, nonnegativity of the recursion gives
$`h_v=0`$ and $`\sqrt{R_v}u_0=0`$, so the same conclusion holds.

Consequently the resulting complete unitary U obeys

```math
(\langle0|\otimes I_N)U(|0\rangle\otimes I_N)
=\iota D_h\mathcal P_W/\alpha.
```

The flag-one input action is supplied by the actual unitary factors, not
reset or discarded. There are $`O(N)`$ local factors with classically
specified scalar data. This is not yet an $`O(N)`$ Clifford+T statement.

## 4. Finite preprocessing at zero defects

Exact zero tests on arbitrary computable angles are unnecessary. First
replace the target local rotations by algebraic real rotations $`U_v'`$
with operator error at most

```math
\delta=\frac{\alpha\eta}{512n^{3/2}}.
```

Rational stereographic coordinates on the unit circle, using either chart,
give certified finite choices. The fixed native coarse words are algebraic.
Every generator, recursive scalar, square root, and selected completion can
therefore be computed with algebraic zero tests and certified enclosures.
Classical running time and bit complexity are separate; no efficiency claim
for these algebraic calculations is used in the quantum counts.

For fixed C, subtree telescoping gives
$`\max_v|h_v(W)-h_v(W')|\le n\delta`$.
The depth-j row block of $`\mathcal P_W`$ is a truncated depth-j frame,
so stacking the block differences gives
$`\|\mathcal P_W-\mathcal P_{W'}\|\le n^{3/2}\delta`$.
Together with $`\|\mathcal P_W\|=\sqrt n`$ and $`|h_v(W')|\le1`$,

```math
\|D_h(W)\mathcal P_W-D_h(W')\mathcal P_{W'}\|
\le2n^{3/2}\delta.
```

Thus the normalized accepted block changes by at most $`\eta/256`$.
Since $`n\delta\lt\varepsilon_0/2`$, the chosen coarse accuracy ensures
all modified subtree discrepancies remain below $`\varepsilon_0`$.
The recursion and denominator bounds therefore still apply. A completion
may change discontinuously at a zero; its accepted block has the uniform
error estimate just proved. The native circuit approximates the selected
complete algebraic unitary.

## 5. A charged native implementation, and its limitation

The following is a conservative implementation, not a T-count lower bound.
Decompose each determinant-one local factor into a constant number of
two-level SU(2) matrices and then Euler rotations. The pair decomposition
and basis routing are standard; see
[Barenco et al., Lemma 4.1 and Section 8](https://arxiv.org/pdf/quant-ph/9503016v1).
There are $`K=O(N)`$ rotation factors overall. Approximating each to
$`\eta/(256K)`$ takes $`O(L+n)`$ native word length with the retained
phase-calibrated one-qubit synthesis primitive.

An arbitrary pair of $`(n+1)`$-bit mode labels can be routed to adjacent
labels by a Gray path and routed back. It uses $`O(n)`$ multiply-controlled
X gates, each charged at $`O(n^2)`$ exact Toffolis using an arbitrary
borrowed helper. For the adjacent pair, use the
[borrowed-sector echo](BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
to implement the conditional rotation. Its actual word
$`C_x=XQ_x^\dagger XQ_x`$ has the exact cancellation symmetry;
inactive inputs and borrowed work return exactly. The same construction
handles an Rz factor by Clifford conjugation. The total per pair is
$`O(L+n^3)`$ T and Clifford gates. The final output permutation needs at
most $`2N`$ transpositions, charged by the same routing method.

Hence this native realization has

```math
T,G=O\!\left(N(L+n^3)\right).
```

The signal flag is part of the data being transformed; it is never assumed
zero inside a subroutine. The second clean flag remains untouched. Routing
uses a fixed number of arbitrary borrowed helpers, well within the selected
dirty allocation, and returns them exactly even with references. The small
cases with fewer controls use fixed-size controlled circuits directly.

Unitary telescoping bounds the error on the entire signal-flag input space
by $`\eta/256`$. Adding preprocessing error bounds the accepted-block
error by $`\eta/128`$, below the proposed $`\eta/64`$ allowance.
Appending the unused second flag gives the requested two-flag interface.
At $`L=N`$, $`n^3=O(N)`$, so both counts are $`O(N^2)`$.
The Clifford target is met; the linear T target is not established.

The existing two-clean stage compiler cannot simply replace this pricing:
the signal flag already contains accepted and rejected amplitudes, leaving
only one fresh initialized flag. That is a limitation of this proposed
substitution, not a general one-clean impossibility theorem. No claim of
additive source costs or optimality of the quadratic implementation follows.

The next endpoint improvement must synthesize this structured unitary more
cheaply, or use a different joint construction. Normalization and a complete
flag action are now explicit for this component. Combining it with the reverse
weighted term, diagonal, coarse C, and amplification still requires a full
composition proof within the same two-clean allocation.
