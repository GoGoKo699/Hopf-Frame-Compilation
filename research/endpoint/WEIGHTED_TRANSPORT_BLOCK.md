# A weighted transport block with one signal flag

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Tree transport](ENDPOINT_TREE_TRANSPORT.md) · [Research status](../../docs/OPEN_PROBLEM.md) · [Borrowed compiler](../../docs/BORROWED_WORKSPACE_COMPILER.md)

The forward weighted transport has an explicit constant-normalization
unitary dilation on one signal flag and the n-qubit logical register.
A recursive scalar correction makes its marker inputs orthogonal without
an extra depth register. Its complete rejected-space action is specified.
The [borrowed-signal extension](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
of the real-rotation compiler gives this block
$`T=O(N+nL)`$ and $`G=O(NL)`$, using $`b\ge L+n+7`$ dirty wires.
Its error includes approximate return of the dirty core and borrowed
amplification signal. No additional initialized flag is needed, and the
dilation signal may be occupied throughout the compilation. An alternative
with exact dirty-work return needs no additional initialized flag and has
$`T=O(L\sqrt N)`$, $`G=O(NL)`$ when
$`b\ge n+\lceil\sqrt N\rceil+7`$.
At $`L=N`$ the first implementation costs $`O(N\log N)`$ T gates.
The linear-T target and the existing full-frame endpoint bound remain
unchanged.

This is a block for one weighted residual term, not a complete frame
compiler. Both native bounds include the physical marker permutations
and the action on both signal-flag sectors, with their respective
dirty-work return guarantees stated explicitly.
The finite matrix checks have dimension at most 32; the native resource
bound is analytic, not an end-to-end circuit fixture.

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

## 5. A native implementation by depth batching

**Proposition.** With the preprocessing and constant normalization above,
suppose $`b\ge n+\lceil\sqrt N\rceil+7`$. There is a coherent
Clifford+T implementation of this weighted block with

```math
T=O(L\sqrt N),\qquad G=O(NL).
```

Its accepted-block error is at most $`\eta/128`$. The signal flag is part
of the data, with a specified complete unitary action; every borrowed wire
returns exactly, jointly with arbitrary references. The second clean flag
is untouched. The proof below first compiles a particular complete unitary
for the algebraic target chosen in Section 4, then uses that section's
accepted-block perturbation bound. It does not assume that completions
vary continuously with the target.

### An allocation that can be packed uniformly

For now use heap logical labels $`0,1,\ldots,N-1`$, with 0 the root input,
and write a mode as $`(F,v)`$, where F is the signal bit. Physical marker
reindexing is charged separately below. Choose the continuation slots

```math
c_1=(0,0),\qquad c_v=(1,v)\quad(v>1).
```

The root rotation acts on $`(0,0),(1,0)`$. At each nonterminal v, order
its four physical slots as

```math
c_v,\quad(0,v),\quad(1,2v),\quad(1,2v+1).
```

The last two slots have not been used before this node. Keep the stop in
$`c_v`$, put the rejection in $`(0,v)`$, and send the two child
continuations to $`(1,2v),(1,2v+1)`$. In this physical slot order, use
the Section 2 completion with output rows reordered as $`(0,3,1,2)`$.
This is an even three-cycle of the last three rows, so the local factor
still has determinant one. Its first two columns give exactly the same
stop, child, and rejection amplitudes as before. Terminal factors use
$`c_v,(0,v)`$, with stop and rejection in that order. Thus the induction
of Section 3 applies unchanged. The unused padding mode is $`(1,1)`$.

At a nonterminal depth $`1\le d\le n-2`$, write each node
$`v=2^d+x`$, with x a d-bit string. Separate the register into signal F,
$`n-d-2`$ leading logical bits, two bits A and B, and a d-bit suffix.
The active parent slots have logical form $`0^{n-d-2}01x`$, whereas
the child slots have $`0^{n-d-2}1xb`$, with signal one. Apply an
A-controlled right cyclic shift to the last $`d+1`$ logical bits. It
uses d Fredkins and sends the child suffix $`xb`$ to $`bx`$, while
leaving parent suffixes unchanged. The four slots now have common suffix
x and three-bit labels

```math
(F,A,B)=101,\quad001,\quad110,\quad111.
```

Apply the fixed three-bit permutation
$`\pi=(0\ 4\ 7\ 3\ 6\ 2\ 5)`$, fixing label 1. Its images of
these four labels are $`000,001,010,011`$, respectively. This constant
permutation has an exact constant-size Toffoli realization: decompose its
seven-cycle into transpositions and route each pair along a three-bit
Gray path. Consequently each node is now one two-mode-bit block, addressed
by x, with all outer bits zero. Pack, apply the block table, and undo the
packing. The complete packing circuit is a permutation on every input;
its action outside the selected blocks need not be interpreted as fresh
workspace. At depth zero, the active labels already are
$`000,001,110,111`$, so CNOT from A to F suffices. The total packing and
unpacking cost over all depths is $`O(n^2)`$ exact gates.

### Fixed determinant-one rotation templates

There is a common constant-length two-level decomposition for every local
SU(4) factor, including singular choices of the prescribed columns. Perform
QR elimination in the fixed pair order

```math
(2,3),\ (1,2),\ (0,1),\ (2,3),\ (1,2),\ (2,3).
```

For a column pair $`(a,b)^{\mathsf T}`$ with
$`r=\sqrt{|a|^2+|b|^2}>0`$, the determinant-one matrix

```math
\frac1r\begin{pmatrix}\overline a&\overline b\\-b&a\end{pmatrix}
```

sends it to $`(r,0)^{\mathsf T}`$. Use identity for a zero pair.
This eliminates the lower triangle with positive first three diagonal
pivots. A triangular unitary is diagonal; the final diagonal entry is also
one because the determinant is one. Reversing the six eliminations
therefore reconstructs the entire factor, including its phase. Express
each two-level SU(2) matrix as three Euler rotations
$`R_zR_yR_z`$, with $`R_z(\theta)=e^{-i\theta Z}`$. This gives at most
18 fixed rotation templates, whose angles depend on the address x.
Zero angles pad missing operations. This is the standard two-level/Euler
construction; see
[Barenco et al., Lemma 4.1 and Section 8](https://arxiv.org/pdf/quant-ph/9503016v1).
Here the determinant-one choice removes any separate address-dependent
scalar-phase gate.

For each template, a fixed permutation of the two mode bits maps its
selected pair to adjacent binary labels. Such a two-bit permutation is
Clifford: the affine permutations on two bits are all 24 permutations.
One mode bit becomes the rotation target, the other the sector control
$`\beta`$; the d address bits remain unchanged. The outer fixed bits
form the predicate for the
[borrowed-sector echo](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector).
For a requested $`R_y(\theta_x)`$ or $`R_z(\theta_x)`$, synthesize
$`Q_x`$ for the corresponding quarter-angle rotation up to scalar phase
and use the actual native word

```math
C_x=XQ_x^\dagger XQ_x,\qquad
\det C_x=1,\qquad XC_xX=C_x^\dagger.
```

Conjugation by X reverses either requested rotation. Thus the same echo
implements $`C_x^2`$ on the selected sector and exact identity on the
others. Scalar phases of $`Q_x`$ cancel. The reflection-table interpreter
applies these determinant-one words with exact dirty-work return; it does
not use the signal as an initialized helper. The fixed mode permutations
are undone after each template.

The terminal depth has $`2^{n-1}`$ factors. For $`n\ge2`$, they are
$`R_z`$ rotations on the signal, addressed by the lower $`n-1`$ logical
bits and selected by the top logical bit being one. The same echo handles
this layer; its outer predicate is empty. Root mixing is one selected
$`R_y`$ rotation. For $`n=1`$, the constant-size two-qubit network can be
decomposed directly by the same two-level and echo construction, with
$`O(L)`$ T and Clifford gates.

### Gathering and physical marker reindexing

After the tree, the stop slots are $`(0,0)`$ for node 1 and $`(1,v)`$
for the other nodes; the padding slot is $`(1,1)`$. Flip F whenever
$`v\ge2`$, and on the two remaining logical values apply the cycle

```math
(0,0)\longmapsto(0,1)\longmapsto(1,1)\longmapsto(0,0).
```

This puts every stop at $`(0,v)`$, padding at $`(0,0)`$, and every
rejection in the signal-one half. The first operation is an unconditional
X on F followed by its inverse conditioned on the upper $`n-1`$ logical
bits being zero. The three-cycle, conditioned on that same zero prefix,
is CNOT from the low logical bit to F, then X on the low bit controlled
on F being zero. Each conditional toggle uses the retained
$`O(n^2)`$ exact-Toffoli construction with an arbitrary borrowed helper.
Thus gathering costs $`O(n^2)`$ and specifies the entire output
permutation, including rejected inputs.

The heap labels are not the prescribed physical marker labels. A node
with binary form $`v=0^r1x`$, where $`|x|=d`$ and $`r=n-d-1`$, has
physical marker

```math
\kappa(v)=x10^r,\qquad\kappa(0)=0.
```

This permutation also has an explicit circuit without initialized work.
For every d, reverse its d suffix bits conditioned on the prefix being
$`0^r1`$. These conditions select disjoint invariant subspaces: none of
the swaps changes the first nonzero bit. Then reverse all n bits
unconditionally. The result is
$`0^r1x\mapsto0^r1\mathrm{rev}(x)\mapsto x10^r`$; zero stays zero.
A prefix-controlled swap of suffix bits a and b is CNOT from a to b,
then X on a controlled by the prefix and b, then the same CNOT. The
middle gate costs $`O(n^2)`$ exact Toffolis with one borrowed helper.
There are $`O(n^2)`$ such swaps, so the reindexing costs $`O(n^4)`$
T and Clifford gates. Apply its inverse before the heap circuit and its
forward version after it, independently of the signal. Both permutations
are therefore charged on the full signal space, not silently absorbed
into the definition of the block.

### Error, workspace, and resource sums

Let J be an absolute upper bound on the number of rotation templates in
one depth layer, also allowing for root mixing. For each rotation template
at depth d choose full selected-sector error at most

```math
\xi_d=\frac{\eta\,2^{d-n}}{512J},\qquad
w_d=O(L+n-d).
```

The retained one-qubit synthesis and exact echo give this error with native
word length $`w_d`$. Within a template, addresses have orthogonal invariant
supports, so the table error is the maximum of its row errors. All packing,
mode routing, and gathering operations are exact. Unitary telescoping over
the constant number of templates per depth and all depths consequently
gives a full-unitary error below $`\eta/256`$, on both signal input
sectors. Together with Section 4's $`\eta/256`$ accepted-block change,
this is below $`\eta/128`$. No correctness claim depends on the rejected
sector being unoccupied inside a subroutine.

For a table with $`S=2^d`$ rows, choose a power-of-two number of dirty
banks within a factor two of $`\sqrt S`$. The exact reflection-table
interpreter costs

```math
T_d=O\!\left((L+n-d)2^{d/2}+n^2\right),\qquad
G_d=O\!\left((L+n-d)2^d+n^2\right).
```

The additive terms charge the sector predicate toggles. The interpreter
uses at most d dirty selectors and $`\lceil\sqrt S\rceil`$ dirty
banks; routing and predicate helpers can be reused after each operation
returns them. The stated $`n+\lceil\sqrt N\rceil+7`$ dirty allocation
therefore suffices throughout. The signal and all logical bits are data;
the second clean flag is never used. Each complete reflection-table interpreter and sector echo returns its
borrowed wires as an exact identity factor; predicate helpers are restored
before the next operation. The guarantee therefore extends to inputs
entangled with external references.

Including the two physical reindexings and all exact routing gives

```math
\begin{aligned}
T&=O\!\left(\sum_{d=0}^{n-1}(L+n-d)2^{d/2}+n^4\right)
  =O(L\sqrt N+n^4)=O(L\sqrt N),\\
G&=O\!\left(\sum_{d=0}^{n-1}(L+n-d)2^d+n^4\right)
  =O(NL+n^4)=O(NL).
\end{aligned}
```

The sums follow after substituting $`j=n-d`$ and summing the convergent
geometric series with coefficients $`L+j`$. The last equalities use
$`N=2^n`$, $`L\ge6`$, and the finite supremum of
$`n^4/2^{n/2}`$; they do not discard an uncharged permutation. All
classical table construction and angle evaluation remain separate
preprocessing, as in Section 4. ∎

At the selected endpoint $`L=N`$, $`b=N+n+7`$, the sufficient dirty
allocation holds and this block has $`T=O(N^{3/2})`$, $`G=O(N^2)`$.
This improves the earlier separate-factor implementation
$`T,G=O(N(L+n^3))`$ by sharing the native word tables across each depth.
It supplies no lower bound on the cost of this block and no improvement
to the existing full-frame endpoint bound.

The implementation above leaves the second clean flag unused and returns
all dirty work exactly. The following implementation reduces the T-count
using only borrowed synthesis work, with approximate core and amplification-
signal return included in its error.


## 6. A faster implementation using a borrowed signal

The [borrowed-signal lemma](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
strengthens the real addressed primitive from an initialized-isometry
estimate to a full-unitary estimate. Its exact X symmetry lets the
amplification signal be arbitrary dirty work. With precision q, k free
address bits, and p unchanged predicate bits, the primitive uses

```math
b_{\rm primitive}=q+k+3
```

arbitrary dirty qubits and no initialized work. Its full-operator error
is less than $`43\,2^{-q}`$, its T-count is $`O(2^k+q+p^2)`$, and its
Clifford count is $`O(2^kq+q+p^2)`$. The reservation consists of the
$`q+1`$-qubit core, k selectors, one predicate helper, and one borrowed
amplification signal. The selectors and helper return exactly; approximate
return of the core and borrowed signal is included in the operator norm.
Inactive predicates give exact identity on the complete input space.

Apply this primitive to the fixed Euler templates of Section 5. A logical
rotation target may be a mode bit or the occupied dilation signal; it is
distinct from the borrowed amplification signal. The unchanged predicate
contains the outer zero bits and the other mode bit's selected value.
Target-only Clifford conjugation gives Rz with the same X symmetry. All
packing, gathering, and physical marker permutations remain exact and
fully charged. No rejected signal sector is assumed empty.

### A fixed sector split preserves the exact dirty threshold

The constant number of templates requires a constant precision allowance.
Increasing the core alone would exceed the specified
$`b=L+n+7`$ threshold. A fixed sector split supplies the needed room
without losing accuracy.

Choose an absolute constant J bounding the number of rotation templates
per depth, including root mixing at depth zero. Set

```math
C=\left\lceil\log_2(43\cdot512J)\right\rceil,\qquad
K=C-4,\qquad q_d=L+n-d+C.
```

These are fixed construction constants except for $`q_d`$; take
$`J\ge19`$, so $`K\ge0`$. For a depth $`d\ge K`$, select K of its
address bits as sector literals and use the other
$`k=d-K`$ bits as the free table address. Process the $`2^K`$ sectors
in turn, adding the selected literals to the template's logical predicate.
The dirty requirement in every sector is exactly

```math
(q_d+1)+(d-K)+1+1
=L+n+C-K+3=L+n+7.
```

The core, selector, helper, and borrowed amplification-signal roles fit
the stated pool without borrowing the logical dilation signal.
The number of sector literals is constant, and the total predicate
length remains $`O(n)`$.

These sector errors take a maximum, rather than a sum. To see why, the
primitive preserves each sector's address bits, and its actual action on
every other sector is exactly identity on the complete flag/core space.
On any one invariant address sector, only its own circuit acts, even if
that circuit leaves small borrowed-signal or core disturbance. Consequently the
product over the sectors has full-operator error bounded by the largest
sector error. This argument uses the complete inactive action, not merely
an identity accepted block. It also explains why the sector split needs no
additional factor $`2^K`$ in precision.

For $`d\lt K`$, use the direct native borrowed-sector echo from Section 5
on each row, to error at most
$`\eta 2^{d-n}/(512J)`$ per template. There are fewer than $`2^K`$
rows in all such shallow layers, an absolute constant. Their native word
length is $`O(L+n)`$, their predicate cost is $`O(n^2)`$, and only a
fixed number of exactly returned borrowed helpers are needed. Root mixing
is included in this direct implementation. This handles all small n as
well, without weakening the dirty threshold. The constants need not make
this split advantageous for small registers.

### Full-input error and total cost

For every compiled deep template,

```math
43\,2^{-q_d}
\le\frac{\eta\,2^{d-n}}{512J}.
```

Let $`U'`$ be the chosen exact dilation of the algebraic target from
Section 4. On the entire logical-signal and dirty-work space, the native
implementation V satisfies

```math
\bigl\|V-(U'\otimes I_b)\bigr\|\lt\frac{\eta}{256}.
```

Sum the template errors over depth, using their maximum over address
sectors. The ordinary unitary hybrid includes arbitrary dilation-signal
inputs, dirty/reference correlations, and approximate return of both the
core and borrowed amplification signal. It needs no initialized synthesis
work or intermediate exact work return. Inverses are actual circuit
inverses. Adding Section 4's preprocessing error gives accepted-block error
below $`\eta/128`$ when the single dilation signal is initialized to zero.

The fixed number of sectors changes only absolute constants. At each
deep depth the costs are

```math
T_d=O(2^d+L+n-d+n^2),\qquad
G_d=O\!\left(2^d(L+n-d)+L+n-d+n^2\right).
```

Including shallow direct synthesis, root mixing, exact packing, and the
$`O(n^4)`$ physical reindexings yields

```math
\begin{aligned}
T&=O(N+nL+n^4)=O(N+nL),\\
G&=O(NL+n^4)=O(NL).
\end{aligned}
```

Here $`n^4=O(2^n)=O(N)`$. The complete result is therefore

```math
b\ge L+n+7,\qquad a=1,\qquad
T=O(N+nL),\qquad G=O(NL).
```

The one initialized flag is the block's dilation signal. Native synthesis
uses only arbitrary dirty work and approximates the complete unitary action
on both signal sectors. Unlike Section 5, return of its source core and
borrowed amplification signal is approximate within the full-operator norm.

At $`L=N`$, the bound is $`T=O(N\log N)`$, $`G=O(N^2)`$, within
$`b=N+n+7`$. It improves the native cost of this component while retaining
a separate precision charge at each depth. It does not improve the
existing $`O(N\log^*N)`$ full-frame bound or close the $`O(N)`$
endpoint. The [two-flag assembly](RESIDUAL_ASSEMBLY.md) now supplies a
complete composition: incorporate the diagonal into a new affine forward
dilation, select it against the reverse block, and amplify before applying
the actual coarse C. Its cost remains $`O(N+nL)`$. A cheaper isolated
forward block would not automatically price that affine component.
