# Precision gauges and native forest compilation

The complete Hopf frame admits a linear-T compiler whenever its required
coarse precision has a linear sum over tree depths. This gives native
$`O(N)`$ synthesis on polynomially regular charts, on every fixed-ratio
interval of nonzero branch exponents, and on the separated families below.
The construction retains the actual rejected operators of native scalar
words. No prepared dirty state or scalar model of a rejection is assumed.

The remaining unrestricted cost is an additive precision sum, identified
in Section 5. These results do not establish an all-instance linear bound.

## 1. Contract and effective normalization

Set $`N=2^n`$, $`n\ge3`$, and $`\epsilon=2^{-N}`$. All native
circuits use H,S,CNOT,T and their adjoints, two clean qubits, and at most
$`N+n+7`$ arbitrary dirty qubits. The approximation is the complete
[initialized-isometry norm](../../docs/FAULT_TOLERANT_COMPILER.md), with
literal phase, every logical column, both clean outputs, arbitrary dirty
inputs, and reference systems included. Classical computation need only
terminate; no efficient classical search is assumed.

A pruned forest means a direct sum of the original binary decision trees
with selected rotations replaced by known corners and the corresponding
zero edges cut after signed-permutation normalization. Every logical
label and remaining subtree is retained in a direct-sum component,
including components absent from the prepared root column. Its maximum
binary depth is at most n. At depth d let
$`S_d`$ count its active internal nodes; always $`\sum_d S_d\le N-1`$.
This inherited decision depth differs from graph distance in the
computational-label tree, whose parent clears the least significant one.

Signed row and input permutations normalize a real frame to
$`0\le s_v\le c_v`$ without losing any columns. Here is the local
identity. Write a binary node as $`\mathrm{diag}(A,B)G_\theta`$,
let m label the right root, let S exchange the two output subtrees, and
let $`S_{\rm desc}`$ exchange corresponding descendant input columns
while fixing the two root columns. Then

```math
S\,\mathrm{diag}(A,B)G_\theta S_{\rm desc}
=\mathrm{diag}(B,A)G_{\pi/2-\theta} Z_m.
```

For quadrant signs $`a=\mathrm{sgn}(c)`$, $`b=\mathrm{sgn}(s)`$,
choose the sign of a zero coefficient to be +1. Then
multiply the output subtrees by a,b and the descendant input columns by
the same respective signs. Give the root input sign 1 and the detail
input sign ab. This changes the node coefficients to $`|c|,|s|`$.
Every child input correction fixes its root, so these identities compose
recursively to $`PWQ=W_+`$. The
[native permutation and sign compiler](DISJOINT_BLOCK_COMPILER.md)
charges both permutations in $`O(N)`$ T gates.

Effective implementation does not decide equality of real oracle values.
Use overlapping certified tests to replace an angle near a Clifford
corner when its smaller branch is below $`2h`$, or retain it when that
branch exceeds h, with $`h=\epsilon/(c_0n)`$ for a sufficiently large
fixed $`c_0`$. One test eventually succeeds. At a single binary depth
the rotations have disjoint supports, so its perturbation is the maximum
local error; the total corner-rounding error is $`O(nh)`$. Retained
branches have exponent at most $`N+O(\log n)`$. Signs and child order
can alternatively be chosen on a certified odd angular grid: a mesh
$`2\pi/M`$ with $`M=2^{N+\lceil\log_2 n\rceil+c_1}`$ perturbs
each depth by at most $`4/M`$ and avoids axes and diagonal ties. Preserve
known exact corners when using a promised nonzero exponent interval.

## 2. The full-column regular gauge and its Parseval source

At node v let the real target first column be $`(c,s)`$, and let the
actual emitted coarse SU(2) word have first column $`(a,b)`$. Its second
column is $`(-\bar b,\bar a)`$; the complex phases are retained. Define

```math
d_{\rm root}=1,\qquad d_L=d_vc/a,\qquad d_R=d_vs/b,
\qquad
\beta_v=\frac{cs}{ab},\quad
\alpha_v=-\frac{c^2\bar b}{a}+\frac{s^2\bar a}{b}.
```

Let $`D`$ contain the leaf gauges, $`W_0`$ the complete coarse frame,
and $`T_0`$ the matrix of its internal subtree-root states. A descendant
gauge maps each coarse subtree root to $`d_v`$ times the target root.
At the local detail column, direct two-coordinate multiplication gives

```math
\mathrm{diag}(d_L,d_R)(-\bar b,\bar a)^T
=d_v\bigl[\alpha_v(c,s)^T+\beta_v(-s,c)^T\bigr].
```

Consequently, with the original detail input label $`j_v`$ retained,

```math
W=D(W_0C+T_0A),\qquad
C_0=1,\quad C_{j_v}=\frac1{d_v\beta_v},\quad
A_{v,j_v}=-\frac{\alpha_v}{d_v\beta_v}.
\qquad\text{(1)}
```

All other entries of A vanish. This proves (1) on the root and every
detail column. If $`a=c(1+e)`$, $`b=s(1+f)`$ with
$`|e|,|f|\le\eta=1/(1000n^2)`$, then $`|D|,|C|\lt 1.01`$ and
$`|A|\le3\eta`$. The cancellation that removes a spurious inverse-small-
branch factor is

```math
\frac{\alpha_v}{\beta_v}
=cs\bigl(|1+e|^2-|1+f|^2\bigr).
```

Known exact corners use their exact coarse gates, equal child gauges,
$`\alpha=0`$, $`\beta=1`$. No division by a zero branch occurs.

An explicit unitary network supplies the internal states in (1). Choose
exact native splitters with first columns $`(r_d,t_d)`$ and
$`1/(16n)\le|t_d|^2\le1/(8n)`$. Start with N leaf modes and N−1
empty feature modes, indexed by the original detail labels. Process
depths bottom-up: merge child carriers with the actual coarse gate
adjoint, retain its detail output, then split the parent carrier against
its fresh feature mode. Put

```math
\lambda_{d+1}=\prod_{k=d+1}^{n-1}r_k,\qquad
\gamma_{j_v}=\lambda_{d+1},\qquad
\gamma_0=\prod_{k=0}^{n-1}r_k,\qquad
\omega_v=t_d\lambda_{d+1}.
```

Both child carriers have the same complex factor before each merge.
The exact initialized source is therefore

```math
P_0x=\begin{bmatrix}\Gamma W_0^\dagger x\\
\Omega T_0^\dagger x\end{bmatrix},\qquad P_0^\dagger P_0=I,
\qquad |\gamma_j|^2\ge7/8,\quad
|\omega_v|^2\ge7/(128n).
\qquad\text{(2)}
```

Parseval follows from the physical merge/split network, rather than an
uncharged square-root completion. The N−1 empty features are modes of
one initialized binary flag, not N−1 initialized qubits.

For an uneven pruned forest of depth D, set
$`\kappa_d=\prod_{k=d}^{D-1}r_k`$ and initially multiply a terminal
leaf at depth d by $`\kappa_d`$. Its sibling carrier then has the same
factor. Fold this multiplier into the row attenuation below. The network
itself remains exactly unitary. Its product with the initial diagonal
contraction has the displayed gamma/omega formulas with D in place of n;
this contracted shared map alone is not Parseval. Its missing norm stays
in Q's attenuation failure bank. There is one root coefficient for each
component and at most N−1 feature coefficients in total.

Take $`\Lambda=5/4`$, $`\kappa=8/5`$. Pre-attenuating
$`D^\dagger/\Lambda`$ and applying (2) gives a row source with shared
part $`[\Gamma W_0^\dagger D^\dagger;\Omega T_0^\dagger D^\dagger]/\Lambda`$.
The column source has shared coefficients
$`C_j/(\kappa\bar\gamma_j)`$ and
$`A_{v,j}/(\kappa\bar\omega_v)`$. They have uniform contraction
slack. Prepare the feature amplitude by an SU(2) bank while the compiler
flag is clean, then scalar-attenuate only the other branch. The branch
intended to retain amplitude 1 uses the exact identity mask of the
[one-clean scalar word](../../docs/ONE_CLEAN_COMPILER.md#literal-diagonal-phases-use-the-same-flag),
including its occupied-flag action. The two sources place rejection in
opposite N-mode banks. Their overlap is exactly $`W/2`$ at the ideal
table level. Section 4 supplies the common recovery argument.

On a fixed regular chart $`\min(|c_v|,|s_v|)\ge N^{-K}`$,
polynomial(n)-length exact native coarse words give the required relative
accuracy. The
[exact borrowed-word interpreter and routing](DISJOINT_BLOCK_COMPILER.md#2-exact-routing-and-coarse-words)
cost
$`\mathrm{poly}(n)\sqrt N=O_K(N)`$ in total. It returns dirty work
exactly; an approximately returned native word is not used as an exact
coarse logical gate. A constant number of final precision-N table calls
therefore proves the regular-chart $`O_K(N)`$ theorem.

## 3. Lossy native scalars supply their own correction features

The following construction replaces exact coarse SU(2) gates when a
branch is very small. At a noncorner node put

```math
\begin{gathered}
t=e^{2i\theta},\qquad a=(1+z)/2,\qquad b=(z-1)/(2i),\\
B=\begin{bmatrix}a&-b\\b&a\end{bmatrix},\qquad
\ell=(1,-i)^T/\sqrt2,\qquad l=1-|z|^2.
\end{gathered}
```

Then $`B\ell=z\ell`$, its orthogonal eigenvalue is 1, and
$`I-BB^\dagger=l\ell\ell^\dagger`$. Use the same node gauges
$`d_L=d_vc/a`$, $`d_R=d_vs/b`$. Since
$`\det[(-b,a)^T,\ell]=-1/\sqrt2`$, direct inversion gives

```math
\begin{aligned}
\mathrm{diag}(d_L,d_R)^{-1}(-s,c)^T
 &=C_v(-b,a)^T+A_v\ell,\\
C_v&=\frac{z-\cos(2\theta)}{2id_vcs},\\
A_v&=-\frac{(z-t)(z-\bar t)}{2\sqrt2\,d_vcs}\\
 &=-\frac{\sqrt2 i(z-t)}{d_v}
   -\frac{(z-t)^2}{2\sqrt2\,d_vcs}.
\end{aligned}
\qquad\text{(3)}
```

Choose a dyadic $`0\lt \tau\le\eta\min(|c|,|s|)/100`$ and an
actual native accepted scalar z satisfying
$`|z-(1-\tau)t|\le\tau/100`$. The elementary bounds in (3) give

```math
\begin{gathered}
l\ge1.9\tau,\qquad |z-t|\le1.01\tau,\qquad
|d|,|d|^{-1}\lt1.01,\\
.98\lt|C_v|\lt1.02,\qquad |A_v|/l\lt.8,\qquad |z|\gt.99.
\end{gathered}
\qquad\text{(4)}
```

For example the leading ratio in the A bound is at most
$`\sqrt2(1.01)/(1.9\cdot.99)`$ and the quadratic correction is
$`O(\eta)`$. Gauges retain their possibly large complex phases.

Let the actual one-clean scalar word satisfy

```math
U_v|0\rangle=|0\rangle\bar z_v I+|1\rangle F_v,
\qquad F_v^\dagger F_v=l_vI.
```

The rejected operator $`F_v`$ may act on all dirty wires. Conjugating
the child pair by the fixed Clifford eigenbasis of Y and applying this
word only to the ell component gives the exact three-port merge
$`[B_v^\dagger;F_v\ell^\dagger]`$. Bottom-up composition yields

```math
Q_0x=\begin{bmatrix}W_B^\dagger x\\
(F_v(L_B^\dagger x)_v)_v\end{bmatrix},\qquad Q_0^\dagger Q_0=I.
\qquad\text{(5)}
```

Here $`W_B`$ is the complete, generally nonunitary coarse frame, and
$`L_B`$ consists of ell injected at each node and propagated through
its coarse descendants. Formula (3) proves the complete-column identity

```math
W=D(W_BC+L_BA).
\qquad\text{(6)}
```

The source has N coarse root/detail modes and N−1 leakage modes. It does
not introduce a second subtree-feature bank. The physical label of a
feature is the original midpoint $`j_v`$. At depth d, CNOT from the
feature bit to logical bit d, the suffix-controlled scalar word, and
the inverse CNOT pair the carrier with that feature. Previously emitted
deeper features have a later nonzero suffix bit and are exactly inactive.
For a pruned forest, the priced mode permutations pack each disjoint
matching into this layout and undo the packing. All predicates act only
on unchanged address bits.

The exact native scalar is obtained from the literal phase source before
amplification. Its accepted polynomial is

```math
p_5(w)=(5-20|w|^2+16|w|^4)w.
```

For $`\rho=1-\tau`$ choose
$`b_\rho=\sin(\arcsin(\rho)/5)`$ and program
$`w=b_\rho e^{-2i\theta}`$. Then $`p_5(w)=\rho e^{-2i\theta}`$.
The polynomial has an absolute Lipschitz bound on the unit disk, so
$`e=\lceil\log_2(1/\tau)\rceil+O(1)`$ programmed bits suffice.
This uses the actual scalar accepted by the emitted word; its rejection
is retained in (5). Inactive predicates condition only central source
Paulis, so the entire word is exactly identity there, including occupied
signal inputs. Loaders cancel with their actual inverses.

At a known exact corner use the exact Clifford
$`B=e^{i\theta}R_y(\theta)`$, put both child gauges equal to
$`e^{-i\theta}d_v`$, and set $`C_v=1/d_v`$, $`A_v=0`$.
Such a node needs no feature query. Its root/detail identity holds
directly, so neither (3) nor a later division by $`l_v=0`$ is evaluated.

## 4. Matching arbitrary dirty rejection and recovering the frame

The column source must reproduce the same $`F_v`$. For any unitary
$`U|0\rangle=[z_0I;F]`$, unitarity gives the full operator identity

```math
U\begin{bmatrix}e^{i\phi}I&0\\0&I\end{bmatrix}U^\dagger|0\rangle
=\begin{bmatrix}
1+(e^{i\phi}-1)|z_0|^2\\
(e^{i\phi}-1)\bar z_0F
\end{bmatrix}.
\qquad\text{(7)}
```

Thus two actual queries, followed by separate literal sector phases,
prepare every normalized pair $`(xI,yF)`$ with
$`|x|^2+l|y|^2=1`$ and $`|y|\le2|z_0|`$. Choose
$`2|z_0||\sin(\phi/2)|=|y|`$ and match the two output phases.
The algebraic order in (7) is essential: its chronological order is
$`U^\dagger`$, phase, U. Reversing it need not reproduce F.

Use outer flag a for shared/failure and inner flag b for root/detail
versus feature modes. For each detail label put

```math
x_v=C_v/\kappa,\qquad y_v=A_v/(\kappa l_v),\qquad
r_v=(|x_v|^2+l_v|y_v|^2)^{1/2},\qquad \kappa=8/5.
```

At a root use $`x=1/\kappa,y=0`$. Bounds (4) imply
$`.6\lt r_v\lt .7`$ and $`|y_v|/r_v\lt 1\lt 2|z_v|`$. First attenuate by
$`r_v`$ on a, placing failure in $`a1b0`$. On $`a=0`$, use (7)
on b to prepare $`(x_v/r_v,(y_v/r_v)F_v)`$, with the same physical
core and the same scalar word as (5). Grouping feature labels by node
depth costs only a constant multiple of the coarse ledger.

The fine signal phases use a as their compiler flag and preserve the
entire address $`(j,b)`$ exactly. A needed occupied-flag fact is as
follows. If a unitary V preserves address and
$`\|VJ_0-J_0D\|\le\delta`$ for a unitary diagonal D, then V
is within $`3\delta`$ of $`D\oplus E`$, with E an address-preserving
unitary on flag 1, for $`\delta\le1/2`$. Indeed the two off-diagonal
blocks have norm at most delta; the lower diagonal block has singular
values at least $`\sqrt{1-\delta^2}`$. Replace it by its polar
unitary and bound the resulting block error. Therefore old P failure
stays in $`b=0`$ up to the charged full operator error. No ideal action
on an occupied flag is assumed.

For Q, first attenuate $`D^\dagger/\Lambda`$ on a, route its failure
to $`a1b1`$, and apply (5) only on $`a=0`$, with $`\Lambda=5/4`$.
At ideal fine-table values the shared columns and their overlap are

```math
\begin{gathered}
P_{\rm shared}=
\begin{bmatrix}C/\kappa\\(A_vF_v/(\kappa l_v))_v\end{bmatrix},
\qquad
Q_{\rm shared}=\frac1\Lambda
\begin{bmatrix}W_B^\dagger D^\dagger\\
(F_v(L_B^\dagger D^\dagger)_v)_v\end{bmatrix},\\
Q^\dagger P=\frac12D(W_BC+L_BA)=W/2.
\end{gathered}
\qquad\text{(8)}
```

The middle equality uses $`F_v^\dagger F_v=l_vI`$ on all dirty
inputs. The two failure banks are orthogonal. Their mode count is
$`(2N-1)+N+N=4N-1`$, within precisely two clean flags.

Let $`J`$ initialize those flags, $`U=V_Q^\dagger V_P`$ be the
actual source-overlap unitary, and $`R=I-2JJ^\dagger`$. The literal
word $`-URU^\dagger RU`$ has, when $`UJ=[A;B]`$, accepted and
rejected blocks

```math
3A-4AA^\dagger A,\qquad B(I-4A^\dagger A).
\qquad\text{(9)}
```

If $`\|A-W/2\|\le\delta`$, these are within $`12\delta`$ of
$`[W;0]`$ for sufficiently small delta. The error is linear, including
rejection. A constant number of fine table calls, each at error
$`O(2^{-q})`$, therefore suffices with $`q=N+O(1)`$. Physical source
inverses are the literal reverse native circuits. They are never replaced
by accepted-block inverses or by an unknown dirty polar factor.

## 5. Precision trace and native resource bound

Let $`e_d`$ be the largest scalar precision needed at forest depth d,
including the relative tolerance $`\eta`$. The
[banked one-clean query](../../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
implements all scalars at that depth. Routing matchings and computing
unchanged suffix predicates cost $`\mathrm{poly}(n)\sqrt N`$ in
total, absorbed in $`O(N)`$. The resulting ledger is

```math
T=O(N)+O\!\left(\sum_{d=0}^{D-1}
\left[\sqrt{S_de_d}+e_d+\mathrm{poly}(n)\right]\right),
\qquad D\le n.
\qquad\text{(10)}
```

The dirty-width reservation in (10) is explicit. For
$`e_d\le N/2-O(n)`$ the banked threshold fits the available pool;
its additional $`S_de_d/N`$ term is at most $`\sqrt{S_de_d}`$.
For larger $`e_d\le N+O(\log n)`$, use the unbanked source instead:
$`S_d+e_d=O(e_d)`$. Self-borrowed masks and central predicate helpers
require only the precision core plus $`O(n)`$ overhead, and the
large-precision implementation uses the core itself for those helpers.
The final fine calls need $`q+1=N+O(1)`$ dirty core wires. Thus all
stages fit $`N+n+7`$ for sufficiently large n. For the finitely many
remaining n, direct native synthesis with the same two flags and no
additional work has a finite worst-case cost; enlarge the absolute
constant. All other gate counts and classical searches are finite.

**Precision-trace theorem.** If $`\sum_d e_d\le C N`$ for a fixed
C, the complete forest has $`T=O_C(N)`$ in the stated endpoint model.
Indeed

```math
\sum_d\sqrt{S_de_d}
\le\sqrt{\left(\sum_dS_d\right)\left(\sum_de_d\right)}=O_C(N).
\qquad\text{(11)}
```

For example, exponents at most $`A N/n+O(n)`$ at every retained
branch give $`e_d\le A N/n+O(n+\log n)`$ and a linear trace. This
includes the cheap forests produced by stable local star reordering.
Classical evaluation of $`A_v/l_v`$ to final precision may need
$`N+\log_2(1/l_v)+O(\log n)`$ oracle bits. This is offline
enclosure precision; it does not enlarge the physical fine source q.

After effective corner rounding every $`e_d\le N+O(\log n)`$.
For the complete balanced tree the square-root sum is already linear,
since $`S_d=2^d`$. The additive $`\sum_de_d`$ can still be
$`\Theta(nN)`$. Reusing a loader across successive scalar words has
not been proved: intermediate source bridges are target-dependent
operators, not automatically cancelling masks. Formula (10) records
this remaining cost rather than claiming an unrestricted compiler.

### Mixed heights and the limit of rescheduling

A pathwise precision budget cannot replace the trace in (10) merely by
rescheduling the same elementary merge banks. The following family makes
the distinction precise. Take integers $`t\ge2`$ and put

```math
n=2^{2t},\qquad N=2^n,\qquad k=t,\qquad
h_r=\frac{n}{2^{r+1}},\quad E_r=\frac{N}{h_r},\quad 0\le r\lt k.
```

Let $`M=N/2^{\lceil\log_2 k\rceil}`$. Give each class one aligned
M-label interval, partitioned into $`M/2^{h_r}`$ aligned blocks of
size $`2^{h_r}`$. Their total size $`kM`$ lies between $`N/2`$ and N.
Set all rotations crossing those blocks to zero and retain a complete
height-$`h_r`$ Hopf component inside each block. This is an actual
pruned frame on the original computational labels. Every smaller branch
in class r may vary independently within
$`2^{-E_r-1}\le s_v\le2^{-E_r}`$.

Choose the standard coarse precision $`q_r=E_r+c_n`$, where one common
$`c_n=O(\log n)`$ supplies the relative tolerance in Section 3.
Every binary root-to-leaf path then has total assigned precision at most

```math
P=\max_r h_rq_r=N+c_nh_0=O(N).
```

The same bound holds on computational parent-child paths inside a block.
These are distinct path notions: the binary leftmost carrier chain below
corresponds to a computational root star, not a computational path.

The preceding sufficient criteria do not directly settle this family with
uniform constants. Compiling whole components as arbitrary SU blocks
would have squared-dimension charge

```math
\sum_j m_j^2=M\sum_{r=0}^{k-1}2^{h_r}\ge M2^{h_{k-1}},
\qquad h_{k-1}=\sqrt n.
```

Even the smallest class exceeds every fixed linear block budget. A fixed
exponent-ratio interval covers only a bounded number of classes, as does
any fixed cheap threshold $`A N/n`$. No branch is ultraweak at the
threshold of Section 7. The complete components have unbounded
computational degree and binary descendant growth, and outside any fixed
cheap or fixed-ratio base their exceptional path depth tends to infinity.
These observations concern the stated sufficient criteria; they do not
exclude a different factorization of the same family.

**Unfused-bank scheduling proposition.** Retain the elementary three-port
merges of (5), their carrier dependencies, and their assigned precisions.
Allow arbitrary ordering and grouping into disjoint matchings. A bank b
has one common precision $`q_b`$ at least that of each merge it contains.
Then

```math
\sum_b q_b\ge q_0h_0+
\sum_{r=1}^{k-1}(q_r-q_{r-1})h_r
=N\frac{k+1}{2}+c_nh_0.
\qquad\text{(11a)}
```

To prove this, select one component from each class. Its binary leftmost
carrier participates in $`h_r`$ successive merges; no disjoint matching
can contain two of those merges. Therefore at least $`h_r`$ banks have
precision at least $`q_r`$. Write each bank's charge as its contributions
above the successive thresholds $`q_0,\ldots,q_{k-1}`$, and sum these
inequalities. Since $`h_rE_r=N`$ and $`E_r=2E_{r-1}`$, the right side
is the displayed expression. Packing or relabeling the modes preserves
the carrier dependencies.

Equality is attained at the scheduling level: align the bottom-up stages
of all components, with one class-r stage in each of the first $`h_r`$
banks. Each such bank is a disjoint matching. Thus the minimum precision
charge is $`\Theta(N\log n)`$, although P is linear. The unbanked query
ledger also gives total $`O(N\log n)`$, since each active node contributes
only a fixed number of query entries per stream and routing costs $`O(N)`$.

This is a bound for the specified scheduling model. The actual paired
loader has $`2q_b`$ one-T splitting rotations per invocation, so separately
emitting the unchanged bank templates retains this order of loader cost.
It is not a lower bound after joint native simplification, replacement of
the sources, or fusion of several merges. In particular, the existing
unrestricted $`O(N[1+\log_2^*(n+2)])`$ compiler is already better than
this direct schedule.

The row-source merge stream alone has the charge (11a). Replaying the
same dirty rejection in the column source and using the actual inverses
in (9) adds only fixed multiplicities; fine precision, flags and dirty
reservation remain those of Sections 4–5. Hoisting a common loader leaves
the [transformed masks and width-transition bodies](SOURCE_REUSE_LIMITS.md#2-factoring-out-the-source-can-move-its-cost-into-the-masks)
to be compiled, so exterior cancellation does not supply a new ledger.
This rescheduling route stops at (11a). The
[heterogeneous-source compiler](HETEROGENEOUS_PRECISION_SOURCE.md) supplies
a different native identity: a prepared geometric source with a tiny
class-dependent ratio is carried through exact dyadic Cayley kernels.
It proves linear cost for this mixed-height family, including preparation,
actual inverses, and full work return. It changes the sources and merge
implementation, so it does not contradict the unchanged-bank proposition.
The unrestricted endpoint question remains open.

## 6. Fixed exponent-ratio families

Suppose every nonzero normalized smaller branch satisfies

```math
2^{-E_{\max}}\le s_v\le2^{-E_{\min}},\qquad
E_{\max}\le K E_{\min},
\qquad\text{(12)}
```

for one fixed K. The interval endpoints and the number of distinct
exponents can depend on n and on the input. This family has
$`T=O_K(N)`$; the same statement holds for a pruned forest.

After corner rounding, the regime $`E_{\min}=O(\log n)`$ is
regular. Otherwise put $`\delta=2^{-E_{\min}}`$ and
$`p=\lceil(N+2n+c_2)/E_{\min}\rceil`$, with a fixed sufficiently
large $`c_2`$. Partition graph depths of the computational-label tree
into slabs of p levels. Moving internal edges left gives
$`W=U_0`$ times dressed cut rotations. A cut from boundary vertex l
of component C to child root r rotates
$`\mathrm{span}\{U_C^\dagger e_l,e_r\}`$. Completed-component
dressing is valid because later internal edges fix that boundary vertex.

For an incoming and outgoing cut through C, their only plane overlap is
$`\langle e_l,U_Ce_{\rm root}\rangle`$, bounded by
$`\delta^{p-1}`$. Writing a planar rotation as
$`I+sA+(c-1)P`$ shows its commutator with the adjacent cut is at most
$`4\delta^{p+1}`$: each generator/projection commutator is bounded
by the plane overlap, and $`1-c\le|s|`$. Stable partition by parent
component-depth parity uses at most $`N^2`$ swaps, hence

```math
\|W-U_0F_{\rm even}F_{\rm odd}\|
\le4N^2\delta^{p+1}=O(2^{-c_2}\epsilon).
\qquad\text{(13)}
```

Each colored factor is a disjoint sum of components with their outgoing
cut roots. A block rooted at a label of trailing degree d has size at most
$`m(d,p)=\sum_{r=0}^p\binom dr`$. Split roots into
$`d\le32p`$ and $`d>32p`$. For high degree, the binary entropy
bound $`H_2(1/32)\lt 1/4`$ gives $`m(d,p)^2\le2^{d/2}`$.
There are at most $`N/2^{d+1}`$ nonroot labels of degree d; the root
contributes at most $`2^{n/2}\le N`$. Thus each high-degree direct
sum has squared block-size sum $`O(N)`$ and is compiled by the
[disjoint SU-block theorem](DISJOINT_BLOCK_COMPILER.md).

Every low-degree internal component has binary depth
$`D\le\min(n,32p)`$. A dressed cut block is implemented as
$`U_C^\dagger K_CU_C`$, where both $`U_C`$ and the boundary-star
forest $`K_C`$ have that depth bound. They are three actual forest
calls; the dressed block is not identified with a single Hopf tree.
Their precisions obey

```math
D(E_{\max}+O(\log n))=O_K(N),
```

so (10)–(11) compile each in linear T. There are at most ten factors:
the low/high parts of $`U_0`$, and three low forest calls plus one
high block sum for each cut color. Allocate constant fractions of
epsilon to their native errors and (13). A unitary hybrid on the same
two flags bounds their complete initialized-isometry error; no reset
between factors is used.

For pruned forests, descendant counts only decrease, and the same degree
distribution and binary-depth bounds apply. Overlapping exponent bins
and the certified corner test of Section 1 avoid exact dyadic or zero
comparisons. Deleting the smallest branches and fixed additive changes
to exponents preserve (12) up to a constant depending on K.

## 7. An arbitrary ultraweak tail

There is a complementary family whose exponent ratio can grow with n.
Fix A and suppose every normalized branch is either cheap, with exponent
$`E\le A N/n`$, or ultraweak, with exponent
$`E\ge N/2+3n/2+c_3`$, for a sufficiently large fixed $`c_3`$.
Their positions and the number encountered on a computational path are
unrestricted.

Use the exact
[stable base-first star reorder](STRUCTURED_FOREST_COMPILERS.md#5-bounded-exceptional-path-depth).
Its input correction is one disjoint sum with squared sizes $`O(N)`$.
Cheap branches lose at most $`n/2+O(1)`$ exponent bits; their
inherited binary depth remains at most n. All ultraweak edges become
local suffixes with sine at most
$`\delta_H=2^{n/2+1-E_H}`$. Delete these suffixes to form cheap
components. Completed-component dressing is exact. Adjacent cuts now
have arbitrary plane overlap at most 1, so stable parity reordering has
error at most

```math
4N^2\delta_H^2=O(2^{-2c_3}\epsilon).
```

The cheap component forests satisfy (11). Each cut-color factor is
$`U_{\rm color}^\dagger K_{\rm color}U_{\rm color}`$, with cheap
forests U and disjoint boundary stars K. The original trailing-degree
distribution gives $`\sum_v(d_v+1)^2=O(N)`$, so K is compiled by
the disjoint-block theorem. Seven forest/star factors, the one input
correction, and the signed permutations therefore give $`O_A(N)`$ T.

The same proof allows the cheap portion to be one fixed-ratio interval
from (12). Reordering changes exponents by $`O(n)`$. If its lower
endpoint is much larger than n, it retains a fixed ratio; otherwise all
retained branches are polynomially regular in N. Section 6 or Section 2
then applies. This proves $`O_K(N)`$ for one such interval mixed with
an arbitrary ultraweak tail. It does not cover unrestricted intermediate
scales, and it does not require bounded exceptional path depth.

## 8. Verification boundary

The accompanying `tests/test_precision_gauge_compiler.py` checks the
complete regular and lossy column identities, complex Parseval source
recurrence, arbitrary-dirty two-query replay, and literal recovery
polynomial. Algebraic fixtures use exact $`\mathbb Q(\sqrt2,i)`$
arithmetic where practical; numerical dirty completions check the full
operator. Exact packing checks and exhaustive small chain schedules check
the mixed-height construction and (11a), without asserting a native
compiler lower bound. These are bounded diagnostics of the displayed identities.
The native count follows from the charged interfaces and estimates above,
not from finite matrix agreement.
