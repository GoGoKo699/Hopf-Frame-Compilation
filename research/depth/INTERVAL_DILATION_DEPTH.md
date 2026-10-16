# General frame depth by interval dilation

[General resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Geometric recovery](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) ·
[Native family sources](COHERENT_DESCENDANT_LOADER.md) ·
[Recursive terminal](RECURSIVE_TERMINAL_DEPTH.md) ·
[Selected routing](SELECTED_DIRTY_ROUTING.md) ·
[Session outcome](GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)

**Local result R86.** At fixed accuracy, arbitrary node-dependent
magnitudes admit a complete-frame compiler of T-depth
$`O_L(N/b^2+\log^{3/2}(n+2))`$, with one external clean qubit and the
retained optimal T-count and linear Clifford-count bounds. The new
mechanism dilates a rounded weighted tree into contiguous intervals of
uniform atoms. Its change of basis to a balanced Haar frame has
entrywise-absolute norm $`O(\sqrt{\log(n+2)})`$. The exposed suffix in
the geometric recovery supplies the dilation registers. Only the coarse
geometric band is recovered by spectral search.

This chapter is the unique proof home for the interval-dilation mechanism.
It does not assert an efficient classical preprocessing algorithm or a
scalable native emitter. Finite matrix checks support its identities;
the asymptotic statement is analytic.

## 1. Theorem and the local interface

Fix an integer $`L\ge6`$, put $`N=2^n`$ and
$`\ell=\log_2(n+2)`$. For every complete prescribed real frame and
integer dirty width $`b\ge17(L+n+7)`$, there is a native Clifford+T
circuit V using one external initialized qubit and at most b arbitrary
dirty qubits such that

```math
\|VJ_1-J_1(W\otimes I_b)\|\le2^{-L},
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell^{3/2}).
```

The corresponding prescribed complex magnitude completion, with supplied
leaf phases, has the same bounds at $`b\ge17(L+n+8)`$. All bounds refer
to the same circuit. The norm includes every logical column, literal
phase, work return and leakage, arbitrary dirty input, and an external
reference. No measurement, reset, supplied precision state, QRAM, or
uncharged coherent evaluator is used. Constants may depend on L; no
precision-uniform strengthening is claimed here.

The inherited geometric recovery preserves a computational prefix u,
a coarse bin label R, and the saved activity bit H. A bin has endpoints
$`a\lt d`$ and width $`q=d-a`$. After its descendant unloader, the ideal
active image is a q-bit local frame band, with an exposed zero suffix.
Write its complete local frame as $`W_{u,q}`$ and its scaling projectors
as $`A_j^{u,q}`$. Its bands are

```math
\Pi_0^{u,q}=A_0^{u,q},\qquad
\Pi_j^{u,q}=A_j^{u,q}-A_{j-1}^{u,q}\quad(1\le j\le q).
```

The old fine tag is $`t=a+j`$. A local compiler for the full
$`W_{u,q}`$, preserving u and R, is sufficient: its inverse converts
$`\mathrm{Ran}\Pi_j^{u,q}`$ into the computational last-one
sector j. XOR that label into the old tag and restore the local frame.
The first bin includes j zero; all other bins have j positive. This
chapter constructs that sufficient compiler. It does not assume a
full-input label-XOR primitive as a premise.

## 2. Finite rounding and the atom embedding

Put $`\eta=2^{-L}`$. Round all node pairs once, before the geometric
decomposition, with

```math
\tau=\left\lceil2\log_2(n/\eta)\right\rceil+12,
\qquad \Delta=2^{-\tau},\qquad K_0=2^\tau,
\qquad \nu=2^{-\lceil\tau/2\rceil-3}\le\sqrt\Delta/8.
```

The integer tau may equivalently be found by finite rational comparisons
ensuring $`\Delta\le(\eta/(64n))^2`$. At each node with original
normalized pair $`(c,s)`$, perform the following finite preprocessing.

1. Obtain a certified rational v with $`|v-c^2|\le\Delta/4`$.
2. Round $`K_0v`$ to a nearest integer, using a fixed rational tie rule,
   then clip that integer to $`1\le k\le K_0-1`$.
3. Obtain rational approximants to c and s with errors at most nu. Select
   each sign to be positive when its approximant is nonnegative, negative
   otherwise. Call the selected signs $`e_0,e_1`$.
4. Set $`c'=e_0\sqrt{k/K_0}`$ and
   $`s'=e_1\sqrt{1-k/K_0}`$.

This rule requires no exact-zero or exact-sign test on the original
coefficients. All selected coefficients are finite algebraic data with
strictly positive magnitudes and exact normalization. The retained
certified coefficient interface makes every step terminate; no
preprocessing running-time claim is made.

Nearest rounding and clipping give $`|k/K_0-c^2|\le\Delta`$. The
square-root inequality bounds each magnitude error by $`\sqrt\Delta`$.
If a selected sign is wrong, the corresponding original magnitude is at
most nu. Therefore

```math
|c'-c|,|s'-s|\le\sqrt\Delta+2\nu
\le\frac54\sqrt\Delta,
\qquad
\|R(c,s)-R(c',s')\|
=\sqrt{(c-c')^2+(s-s')^2}\lt 2\sqrt\Delta.
```

The [addressed-layer product](../../docs/HOPF_INTERFACE.md) has n layers.
Each layer is a direct sum of its node rotations and identity on the
complementary zero-suffix predicate. Thus its norm perturbation is the
maximum node perturbation at that depth. Unitary telescoping gives

```math
\|W-W_{\rm rnd}\|\le2n\sqrt\Delta\le\eta/32.
\qquad(1a)
```

This is a full operator-norm estimate. Every node below an original zero
incoming amplitude is still rounded from its own prescribed pair; none
of its marker columns is deleted. Tensoring arbitrary dirty systems or
references preserves the estimate.

Define the rounded frame once. Every descendant source, local atom tree,
scaling projector, and terminal factor below uses that same tuple of
rounded coefficients. Its finite signs have the exact retained gauge

```math
W_{\rm rnd}=D_{\rm sign}W_{\rm mag}S_{\rm sign}.
\qquad(1b)
```

Here the output leaf sign is the product of selected signs on its path;
the state-input sign is one, and the marker sign at v is the product of
its ancestor edge signs times $`e_0(v)e_1(v)`$. The all-column gauge
is proved in the [common-flag chapter](HOMOGENEOUS_MAGNITUDE_DEPTH.md#2-separate-signs-without-a-sign-oracle).
Section 10 prices both sign diagonals and the prescribed complex leaf
phase. In particular phases at original zero-state leaves are retained.

Fix one positive rounded q-level local frame W. Put
$`r=q\tau`$ and $`D=2^r`$. Every leaf probability is
$`\mu_x=m_x/D`$ for an integer $`m_x\ge1`$. In lexicographic leaf
order, partition the integer atom interval $`[0,D)`$ into consecutive
half-open intervals

```math
I_x=[L_x,L_x+m_x),\qquad L_x=\sum_{z\lt x}m_z.
```

For each tree prefix v, its descendant union $`I_v`$ is also an
interval. Its two children have the rounded prescribed probability
ratio. All lengths are integers: a depth-e prefix numerator is
multiplied by $`2^{(q-e)\tau}`$ in the common denominator.

Define the isometry from q leaf bits to r atom bits by

```math
E|x\rangle=\frac1{\sqrt{m_x}}\sum_{z\in I_x}|z\rangle.
\qquad(1)
```

Disjoint interval supports give $`E^\dagger E=I`$. If $`v_v`$ is the
normalized positive scaling vector of prefix v, then

```math
Ev_v=|I_v|^{-1/2}\mathbf1_{I_v}.
\qquad(2)
```

Consequently EW consists of the constant atom vector and the normalized
interval-tree Haar wavelets. For a parent interval I with child lengths
$`l_0,l_1`$ and $`L_I=l_0+l_1`$, its wavelet has, up to the fixed
marker sign convention, the two values

```math
-\sqrt{\frac{l_1}{L_Il_0}}\quad\hbox{on }I_0,
\qquad
+\sqrt{\frac{l_0}{L_Il_1}}\quad\hbox{on }I_1.
\qquad(3)
```

### Native realization and actual inverse

On leaf X, atom Z, and a zero source bank, use the generic real family
source to prepare $`|I_x|^{-1/2}\mathbf1_{I_x}`$ into Z while
preserving X. Then apply the exact arbitrary-output XOR lookup

```math
X\longmapsto X\mathbin\oplus F_u(Z),
\qquad F_u(z)=\text{the unique leaf x with }z\in I_x.
\qquad(4)
```

On its initialized image this clears X and realizes (1), including its
positive phase. Unlike quantum row-label removal, (4) has a deterministic
classical inverse relation on every atom string. All finite lookup data
are computed in preprocessing. No coherent interval arithmetic is free
or needed.

The physical inverse is the actual reversed word: membership lookup,
then the actual source inverse. Its compression from X zero to atom and
source-work zero is $`E^\dagger`$, regardless of its completion on
other active inputs. Approximate source return is included in the
source-isometry norm. Every lookup returns arbitrary dirty work exactly.
Including H in each source and lookup enable gives exact identity at
$`H=0`$ on all occupied local registers.

## 3. A square-root normalization bound

Let B be the complete balanced Haar frame on the r atom bits, and put

```math
K=B^\dagger EW.
\qquad(5)
```

This rectangular matrix is an isometry. Its root row and root column form
an isolated unit entry: EW's root is the constant atom vector, and all
its other columns have zero mean. The remaining entries are
$`\langle h_J,w_I\rangle`$, where $`h_J`$ is a balanced Haar
wavelet on a dyadic atom interval J and $`w_I`$ is (3). Write

```math
s_J=|J|,\qquad m_I=\min\{|I_0|,|I_1|\},\qquad L_I=|I|.
```

Rounded positivity gives $`L_I/m_I\le2^\tau`$. The following proof
is uniform in the tree height, atom dimension, and all selected pairs.

**Sparse row-energy estimate.** If a matrix M has at most d nonzero
entries in each column and every row has squared Euclidean norm at most
e, then $`\|M\|\le\sqrt{de}`$. For each row, Cauchy--Schwarz
bounds $`|(Mx)_j|^2`$ by e times the squared norm of x restricted to
that row's support. Summing counts every coordinate at most d times.
This applies equally to the entrywise absolute matrix.

A piecewise constant wavelet $`w_I`$ has at most three jumps: its
endpoints and child split. At a fixed dyadic scale, only the cells
containing those jumps can have nonzero Haar overlap. There are at most
three such cells. Half-open integer intervals avoid boundary ambiguity;
a jump exactly at a cell endpoint contributes no additional cell.

### Middle scales

Restrict entries to $`m_I\le s_J\le L_I`$. Each column has at most
$`3(\tau+2)`$ entries. Every row has squared norm at most one by
Bessel's inequality for the orthonormal family $`\{w_I\}`$. Hence

```math
\|\,|K_{\rm mid}|\,\|\le\sqrt{3(\tau+2)}.
\qquad(6)
```

### Fine scales

For each atom x and each positive threshold A,

```math
\sum_{I:m_I\ge A}|w_I(x)|^2\le\frac1A.
\qquad(7)
```

Only ancestors of x contribute. At such an ancestor,
$`|w_I(x)|^2=1/|I_{\rm child}(x)|-1/|I|`$. Add any omitted
nonnegative ancestor terms through the deepest qualifying node, then
telescope. The child at that node has length at least A. The qualifying
nodes need not form an initial segment; this subset argument is why (7)
remains valid on a heavy-child path.

Minkowski's inequality and $`\|h_J\|_1=\sqrt{s_J}`$ now give

```math
\sum_{I:m_I\ge A}|\langle h_J,w_I\rangle|^2\le s_J/A.
\qquad(8)
```

Slice the fine range $`s_J\lt m_I`$ by
$`2^k\le m_I/s_J\lt 2^{k+1}`$, $`k\ge0`$. Every column meets at
most two dyadic scales, hence at most six rows. By (8), each row's squared
norm is at most $`2^{-k}`$. The slice norm is at most
$`\sqrt6\,2^{-k/2}`$, and their sum is bounded by
$`\sqrt6/(1-2^{-1/2})`$.

### Coarse scales

Fix $`h_J`$ and set $`A=s_J2^{-k}`$. If
$`L_I\le A`$ and $`\langle h_J,w_I\rangle\ne0`$, I must
cross a jump of $`h_J`$; otherwise $`h_J`$ is constant on I and the
zero mean of $`w_I`$ kills the overlap. All these I lie in the union S
of at most three half-open integer neighborhoods of radius
$`\lfloor A\rfloor`$ around those jumps. Each contributing interval
has integer length at most $`\lfloor A\rfloor`$ and crosses an
integer boundary. Thus $`|S|\le6\lfloor A\rfloor\le6A`$; when
$`A\lt1`$ there are no contributing intervals. They remain orthonormal
and supported in S, so

```math
\sum_{I:L_I\le A}|\langle h_J,w_I\rangle|^2
\le\|h_J\mathbf1_S\|_2^2\le6A/s_J.
\qquad(9)
```

In a slice $`2^k\le s_J/L_I\lt 2^{k+1}`$ of the coarse range,
each column again has at most six entries, and (9) bounds row energy by
$`6\,2^{-k}`$. The slice norm is at most $`6\,2^{-k/2}`$.

Combining the isolated root entry and these three ranges proves

```math
\|\,|K|\,\|
\le 1+\sqrt{3(\tau+2)}
 +\frac{\sqrt6+6}{1-2^{-1/2}}
 =O(\sqrt{\tau+2}).
\qquad(10)
```

The constants are deliberately loose. This is a bound on a specified
change of basis, not a general unitary synthesis assumption.

## 4. A native rectangular overlap

All compression and amplification identities in Sections 4–6 are on a
fixed preserved prefix u, bin R, and active $`H=1`$ sector. The physical
inactive action is separately proved to be identity on every occupied
input. Thus scaled blocks below are not assertions about $`H=0`$.

Choose one rational $`\alpha`$ strictly larger than twice the right
side of (10), uniformly across all prefixes and bins. Let
$`M=|K|`$ and define the positive Schur weights

```math
\mathcal B_M=\begin{pmatrix}0&M\\M^{\mathsf T}&0\end{pmatrix},
\qquad
\binom fg=(I-\mathcal B_M/\alpha)^{-1}\mathbf1.
\qquad(11)
```

The Neumann series converges with a known strict gap. All entries of f
and g are at least one, and

```math
Mg=\alpha(f-1),\qquad M^{\mathsf T}f=\alpha(g-1).
```

This is the rectangular version of the certified Schur construction in
[R82](BOUNDED_DENSITY_DEPTH.md#5-certified-schur-weights-and-the-sparse-source-interface).
All matrix entries are effectively specified finite algebraic numbers
after dyadic rounding, so signs and exact zeros are decidable in
preprocessing. No coherent eigensolver or Schur solver is invoked.
Large classical values do not require quantum registers storing those
values; the native source prepares normalized states to absolute norm
error.

For atom y and leaf x, put

```math
c_{yx}=\mathrm{sign}(K_{yx})
 \sqrt{\frac{|K_{yx}|f_y}{\alpha g_x}},
\qquad
l_{yx}=\sqrt{\frac{|K_{yx}|g_x}{\alpha f_y}}.
\qquad(12)
```

Both are zero at a zero entry. Their overlap is $`K_{yx}/\alpha`$,
while the squared column norm of c is $`1-1/g_x`$ and the squared
row norm of l is $`1-1/f_y`$.

Use registers X of q bits, Z of r bits, two failure bits F, and a common
padded source bank A. The column source preserves X and prepares

```math
|x\rangle\longmapsto
 \sum_y c_{yx}|x,y,00\rangle
 +g_x^{-1/2}|x,0,01\rangle.
\qquad(13)
```

The row source preserves Z and prepares

```math
|y\rangle\longmapsto
 \sum_x l_{yx}|x,y,00\rangle
 +f_y^{-1/2}|0,y,10\rangle.
\qquad(14)
```

These are normalized real-vector family sources, supplied by the charged
native source theorem. Their failure supports are orthogonal. With

```math
J_{\rm in}|x\rangle=|x,0,00,0^A\rangle,
\qquad
J_{\rm out}|y\rangle=|0,y,00,0^A\rangle,
```

and physical source unitaries $`U_R,U_L`$, the exact-source comparison
has

```math
J_{\rm out}^\dagger U_L^\dagger U_RJ_{\rm in}=K/\alpha.
\qquad(15)
```

If each actual initialized-source error is at most delta, the actual
block differs by at most $`2\delta`$. Indeed the block is the overlap
$`(U_LJ_{\rm out})^\dagger(U_RJ_{\rm in})`$. This estimate needs no
specified full-space action on occupied source work. The sources preserve
their read-only labels and are exact identity at $`H=0`$ on every input.

## 5. Square compression without rejected-history return

Add a separate q-bit register Y, a separate interval-source bank
$`A_E`$, and a balanced-frame tag $`T_B`$, all initially zero.
Chronologically apply

```math
U_K=U_L^\dagger U_R,\qquad
B\text{ on Z},\qquad
E^\dagger\text{ on }(Y,Z),\qquad
\mathrm{SWAP}(X,Y).
\qquad(16)
```

Every operation is guarded by H. B is the exact guarded balanced native
word from
[R82 Section 6](BOUNDED_DENSITY_DEPTH.md#6-the-balanced-factor-is-an-explicitly-guarded-native-word),
with its own returned short tag. It is not a supplied frame oracle.
The swap bank is an exact returned-dirty controlled-involution echo.
The completed word V is identity at $`H=0`$ on arbitrary occupied
registers.

Let J initialize all work other than X. Then

```math
J^\dagger VJ=W/\alpha
\qquad(17)
```

in the exact-source comparison. To see why rejection cannot contaminate
(17), a nonzero old X after $`U_K`$ is untouched by B and
$`E^\dagger`$, then moves into final Y. Failure F and source work A
are never touched by later factors. Any of those rejections is therefore
excluded by the final anchor. On the remaining block, B converts K into
EW, and the inverse embedding on the fresh Y converts this into W with
Z and $`A_E`$ zero. The final swap places W on X. The separate
$`A_E`$ bank is essential: its rejected histories also remain outside
J. Every inverse is the actual reversed native word.

The same compression argument with actual source columns gives

```math
\|J^\dagger V_{\rm act}J-W/\alpha\|\le C\delta
\qquad(18)
```

for an absolute constant C. The balanced frame is exact; all other
continuous preparations are included in delta. Equivalently one may use
the dimension-free source-isometry extension lemma and unitary
comparison. Neither proof assumes exact actual source-work return or
an initialized-tail statement on leaked inputs.

## 6. Uniform calibrated amplification

Choose an integer $`k>\alpha`$, for example
$`k=\lceil\alpha\rceil+1`$, and set

```math
\theta=\frac\pi{4k+2},\qquad
\zeta=\sin\theta,\qquad \rho=\alpha\zeta\lt 1.
```

Append a damping contraction rho using a fresh flag and source bank.
Including their errors, the complete actual unitary U has

```math
\|J^\dagger UJ-\zeta W\|\le\epsilon_{\rm blk},
\qquad \epsilon_{\rm blk}\le C\delta.
\qquad(19)
```

All new work belongs to J. U remains exact identity at $`H=0`$.
The full-isometry recovery argument of
[R82 Section 7](BOUNDED_DENSITY_DEPTH.md#7-a-complete-guarded-constant-normalization-early-circuit)
is uniform in this growing k. In fact $`\zeta\le1/2`$; for
$`\epsilon_{\rm blk}\le1/16`$, the rejection Gram operator
$`I-A^\dagger A`$, $`A=J^\dagger UJ`$, has an absolute positive
lower bound. The square-root Lipschitz estimate and the isometry-extension
constant are therefore absolute, not proportional to $`1/\zeta`$.
There is a comparison unitary $`U_*`$ with exact block $`\zeta W`$
and

```math
\|U-U_*\|\le C'\epsilon_{\rm blk}.
```

Take that comparison sectorwise in the preserved prefix, bin, and H;
on the inactive sector choose $`U_*=U=I`$. Define the exact anchored
reflection and scalar phase by

```math
R_H=I-2HJJ^\dagger,\qquad
Z_H=I-2H.
```

The physical word is

```math
\mathcal F(U)=[Z_HUR_HU^\dagger R_H]^kU.
\qquad(20)
```

It contains $`2k+1`$ actual forward/inverse occurrences. The scalar
$`Z_H`$ removes the literal recovery sign. Since
$`(2k+1)\theta=\pi/2`$, the ideal accepted amplitude is one and
the rejected amplitude is zero. Thus

```math
\|\mathcal F(U)J-JW\|
\le C'(2k+1)\epsilon_{\rm blk}.
\qquad(21)
```

At $`H=0`$, (20) is identity on the whole occupied-work space. Set
$`\delta=c\epsilon_{\rm loc}/(k+1)`$ for sufficiently small
absolute c. The source precision is

```math
p_{\rm src}=O\!\left(\log\frac{k+1}{\epsilon_{\rm loc}}\right),
\qquad k=O(\sqrt{\tau+2}).
\qquad(22)
```

At fixed local accuracy this is $`O_L(\log(\ell+2))`$. The same
physical source banks are reused: no reset or exact intermediate return
is assumed. The all-work reflection is a charged returned-helper
conjunction, with logarithmic depth and polynomial count in the number
of tested wires. It does not require an arbitrary truth table on every
initialized wire.

## 7. Reached-band erasure and selected native families

Let $`\widehat W_{u,q}`$ denote the actual compiler (20). On the ideal
unloaded band, apply its actual inverse, XOR
$`a+\mathrm{lastone}(X)`$ into the old fine tag T, then apply
$`\widehat W_{u,q}`$. In the exact comparison this gives

```math
|\psi_j\rangle_X|a+j\rangle_T|0^{\rm work}\rangle
\longmapsto
|\psi_j\rangle_X|0\rangle_T|0^{\rm work}\rangle,
\qquad \psi_j\in\mathrm{Ran}\Pi_j^{u,q}.
\qquad(23)
```

The compiler uses the retained bin and prefix, not the erased T, so the
last inverse dependence is explicit. Equation (23) holds on the whole
orthogonal sum of reached bands and prefixes. Its error is at most twice
the local compiler error, by ideal-source-image telescoping through actual
unitaries. The source leakage is never treated as a clean suffix.
The effective enable on the label XOR includes H. On inactive input the
complete local eraser is therefore literal identity on arbitrary work.

The local extra clean allocation is

| Register | Initialized bits beyond the local input |
| --- | ---: |
| Atom register Z | $`r=q\tau`$ |
| Temporary output Y | q |
| Overlap failure bits F | 2 |
| Common row/column source bank A | $`O(p_{\rm src})`$ |
| Separate interval-source bank $`A_E`$ | $`O(p_{\rm src})`$ |
| Balanced-frame tag $`T_B`$ | $`O(\log(r+2))`$ |
| Separate damping flag and bank | $`O(p_{\rm src})`$ |

A fixed external-bank layout may reserve one additional q-bit donor
block. We reserve $`r+2q`$ exposed donor wires and a separate
$`O_L(\ell)`$ protected metadata bank for the fixed ancillary fields
in the table. The old fine tag, coarse bin, saved predicate, and outer
source records have their own disjoint reservation; none is counted
again as an available donor.

For a bin with preserved prefix length a, the significant large source
addresses are

| Operation | Relevant data-address bits |
| --- | ---: |
| Schur column source | $`a+q+r`$ |
| Schur row source | $`a+r+q`$ |
| Interval source | $`a+q+r`$ |
| Atom-to-leaf membership | $`a+r`$, with q outputs |

Short bin, expansion-index, calibration, and enable fields multiply row
counts only by a fixed polynomial in n. There are a constant number of
source words in one U. All labels, output lengths, masks, and invalid
padding targets are finite preprocessing data. The family-source theorem
allows arbitrary normalized real vectors including failure coordinates.
Its error is a supremum over preserved sectors, not a sum over the number
of prefixes.

### Coherent variable sizes

Use one common source grammar across all bins, not a serial call to a
whole source for each bin. Choose the same expansion length and the same
source-index/calibration registers. At bin i, each finite coarse source
vector is a vector on only its active target prefix of length $`q_i`$
or $`r_i`$, together with the constant failure fields. Each of the
constant number of grammar positions from the
[native source proof](COHERENT_DESCENDANT_LOADER.md#3-constant-query-flagging-and-literal-inner-amplification)
is compiled jointly: coefficient rotations act on the common source
bank; selected Hadamards act on the active target prefix and failures;
Boolean phases use the relevant prefix, row, active target, and index;
calibration reads the preserved label; and the inner zero reflection
tests only the active target, failures, index, and attenuation field.
The other inner reflection tests only its common source core.

Thus no maximum-size zero-padded source is assumed. Inactive bank padding
is a tensor identity after each completed selected source macro. Raw
gate prefixes inside a returned query may transiently correlate public
inputs and private dirty helpers; their exact unquery is included. The
native source's variable-suffix masks are replaced by variable-prefix
masks selected by the same short bin. Source-isometry errors are uniform
over the selected sectors and over arbitrary states of the padding.

Do not pad the large lookup address with the union of all logical prefix,
row, and atom banks. Compile bin i's table on only its
$`a_i+q_i+r_i`$ relevant inputs. All bin queries can run as one parallel
family with private dirty outputs and private returned helpers; public
address wires are used only as read-only CNOT controls. Disjoint bin
predicates make their selected phases multiply correctly. To combine
private outputs into a shared XOR target, use the chronological echo

```math
Q_{\rm all},\quad \mathrm{CNOT}_{\rm private\to target},\quad
Q_{\rm all}^{\dagger},\quad \mathrm{CNOT}_{\rm private\to target}.
\qquad(24)
```

Unknown dirty offsets cancel. Shared targets must be disjoint from every
public address wire. Non-Clifford supports remain private, so the table
family does not introduce one serial source layer per bin. Its total
volume is bounded by the sum of the relevant table volumes.

Selected reversible layouts use a fixed Beneš network with classically
specified switch tables depending on the short bin R. Each stage is an
exact dirty-control swap echo with private parallel helper pools.
There are $`O(\log n)`$ stages; for $`p=O_L(\ell^2)`$ bins their
T-depth is $`O(\ell\log\log(\ell+3))`$, with polynomial count and
dirty width. The layout is established around the complete local
amplification and reversed afterwards. Padding wires are preserved;
only selected live fields enter the local anchor. A variable-size
balanced frame uses the exact guarded balanced word on its selected atom
prefix, with returned short-tag masks. Fixed tensor-product factors
cancel on the unused padding. These layout and mask operations have no
extra external clean input. These complete routing and selected-query
identities, including every private simultaneous pool, are proved in
[the supporting routing lemma](SELECTED_DIRTY_ROUTING.md).

The amplification anchor tests only the selected live Y and Z fields
and the local failure, source, damping, and balanced-tag fields. It does
not test the original logical slots holding parked arbitrary dirty
values, inactive bank padding, preserved prefix, old fine tag, or coarse
bin. For each bin, form the conjunction of H, its bin equality, and its
own initialized-field zero tests, using private returned helpers. The
parallel phase echo combines these disjoint cases in $`O(\ell)`$
T-depth and polynomial count. No arbitrary truth table on the union of
all work fields is used. The helper pools exclude the entire union of
tested and preserved public wires.

## 8. Local resource and error envelope

For a coherent family of bins, let

```math
\mathcal Q=\sum_i2^{a_i+q_i+r_i}.
```

The native source and sharp dirty-query interfaces, including their
calibration output widths, private helper pools, exact membership,
returned predicates, masks, and actual inverses, give a fixed-degree polynomial
$`P_L(n)`$ such that each complete local U obeys

```math
T,\ w\le P_L(n)\sqrt{\mathcal Q},\qquad
G\le P_L(n)\mathcal Q,
\qquad D_T=O_L(\ell+p_{\rm src})
```

at its standalone sufficiently wide dirty allocation. The square-root
bound may absorb the polynomial number of separate-bin pools using
$`\sum_i\sqrt{Q_i}\le\sqrt{p\sum_iQ_i}`$. The routed width-dependent
query terms are charged in the outer ledger. In particular, neither the
selected family nor its output width is free.

Amplification uses $`O(\sqrt{\tau+2})`$ copies, so at fixed local
error, after the selected layout is established,

```math
D_{\rm local}=O_L(\ell\sqrt{\tau+2})=O_L(\ell^{3/2}).
\qquad(25)
```

The exact balanced word, zero-work reflections, and selected last-one
XOR fit the same depth allowance. Layout depth is
$`O_L(\ell\log\log(\ell+3))=O_L(\ell^{3/2})`$. The repetition
count is only polynomial in n and is included in the fixed polynomial
count and workspace envelope. The global count margin must absorb that
polynomial; it cannot be omitted merely because each individual bin
fits its donor bank.

The local error in (23) is a single initialized-isometry error uniform in
all selected prefixes, bins, and arbitrary dirty references. It is
allocated as a fixed fraction of the outer error, while individual
source errors use (22). Global rounding contributes once. Outer source
and coarse-reflection errors are handled by the inherited geometric
word with actual inverses; Section 9 supplies their allocation and
preserved-sector proof.

## 9. Outer geometric allocation and global ledger

We now price the complete positive rounded early circuit. All choices
below depend on fixed L but not on the node coefficients. First fix the
short metadata allocation and a free reserve
$`m_0=\lceil c_L\ell\rceil`$ with $`m_0\ge4(\tau+2)`$.
Metadata includes the old fine tag, coarse selector, search and threshold
records, outer source work, and the local failure, source, damping, and
balanced-frame tags. Their total is $`O_L(\ell)`$ and is disjoint
from the free reserve. Only afterwards choose

```math
M=\lceil A_L\ell\rceil,\qquad h=n-M.
```

The large-dimension branch has $`h\ge1`$. The remaining finite range
is handled in Section 10. Extra protected wires may remain unused;
increasing $`A_L`$ does not enlarge $`c_L`$.

### Thin bins and actual donor space

Put $`\kappa=\tau+2`$. Construct bins backwards from d=h. At each
step set

```math
v=h-d+m_0,\qquad
q=\min\left\{d,\left\lfloor\frac{v}{2\kappa}\right\rfloor\right\},
\qquad a=d-q,
```

record $`(a,d]`$, and continue from d=a. Every nonempty step has
$`q\ge1`$, and

```math
v\ge2q(\tau+2)=2(r+2q),\qquad
 a+q+r=d+r\le h+m_0.
\qquad(26)
```

After the selected endpoint unloader, the early tail of length h-d is
zero on the ideal reached band. Together with the free reserve it gives
v clean donor wires. Equation (26) fits the atom and temporary banks
with a constant spare factor; persistent metadata is not borrowed.
Every nonterminal step grows v by at least $`1+1/(4\kappa)`$, so

```math
p=O_L\!\left(\kappa\log(1+h/m_0)\right)=O_L(\ell^2),
\qquad k_{\rm bin}=\lceil\log_2 p\rceil=O_L(\log(\ell+2)).
\qquad(27)
```

The largest q is $`O_L(n/\ell)`$ and the largest r is O(n). The
[supporting geometry and layout proof](SELECTED_DIRTY_ROUTING.md#3-thin-bin-geometry-for-an-atom-register-of-size-q-tau)
gives fixed external arbitrary-dirty banks with the following selected
placement. Move the active local data into the active part of the fixed
leaf bank, and swap donor zeros into the active atom and temporary banks.
Park their former arbitrary contents in the original local and donor
slots. The active body preserves those parked values; the completed
actual inverse router returns them exactly. Any approximation residual
on the initialized fields returns to the corresponding logical exposed
slots and belongs to the full operator error.

Keep this layout fixed throughout the local compiler, its actual inverse,
all amplification reflections, and the fine-tag XOR. Routing is performed
only at the two boundaries of that entire local episode. The completed
router preserves the selected prefix and all metadata. Coherent sums of
bins and arbitrary dirty reference correlations are included.

### Complete early recovery word

Use the rounded positive early-frame tagged source of R85. With increasing
endpoints $`0=d_0\lt d_1\lt \cdots\lt d_p=h`$, define the zero-based bin by

```math
g(t)=\min\{i:0\le i\lt p,\ t\le d_{i+1}\}.
```

Thus the first bin contains the root t=0. Compute $`R\mathrel{\oplus}=Hg(T)`$ from the
old fine tag by an exact finite lookup. The selected descendant inverse
at $`d_{R+1}`$ exposes the local frame and donor suffix. In the selected
layout, run the local compiler inverse, XOR
$`H(a+\mathrm{lastone}(X_{\rm loc}))`$ into T, and run the compiler.
Reverse the layout and reload the selected descendants. In the ideal
comparison T is now zero, while R still records the coarse band and the
original early-frame vector is restored with its literal phase.

Erase R using the inherited coherent binary search over the p coarse
bands, a fresh search register, XOR into R, and the actual reverse search.
Its threshold comparison $`g(t)\le i`$ is the scaling-projector test
$`t\le d_{i+1}`$. Padded search thresholds use endpoint h and hence
an always-true comparison. The transported reflection still includes
its initialized source-work projector: its ideal full-space form is
$`I-2(A_h\otimes|0\rangle\langle0|_A)`$, which is minus identity
on initialized source work, not identity. This is exactly the full-input coarse eraser
of [R85](GEOMETRIC_BAND_ORIENTATION_DEPTH.md); it uses
$`k_{\rm bin}`$ comparison rounds, and no initial coarse spectral
search is needed. All transient search and threshold fields return on
the ideal output.

The geometric source and unloader identity uses no shared-magnitude
promise. That promise entered R85 only in its local flag transporter,
which Sections 2–8 replace. Set the completed call tolerances, in the
same norm interfaces as R85, to

```math
\delta_0=\delta_U=\delta_{\rm loc}=\eta/128,
\qquad
\delta_R=\frac{\eta}{512(k_{\rm bin}+1)}.
```

The inherited ideal-source-image hybrid gives

```math
\epsilon_{\rm early}
\le\delta_0+2\delta_U+2\delta_{\rm loc}
       +8k_{\rm bin}\delta_R
\lt \eta/8.
\qquad(28)
```

Here $`\delta_{\rm loc}`$ is the complete local compiler error from
(21), not the error of one of its sources. Its individual source errors
are chosen by (22). Exact classical labels, routers, balanced transforms,
and structural predicates add no approximation error. All inverses are
actual inverses. The hybrid compares with ideal reached states when it
invokes exposed zeros and propagates every actual leakage component by
unitaries; no leaked input is presumed initialized.

Save H from the protected zero predicate before using any protected
work. Every completed guarded source and local recovery word is identity
on the full H-zero sector, including occupied data and work. Threshold
preparations meet their actual inverses; the coarse eraser is its actual
search/inverse pair with the middle XOR disabled. Thus the entire early
word has literal inactive identity. Finally reverse the actual
saved-predicate computation. On the ideal active image all protected
fields have returned, so H clears. On actual leakage that last inverse
merely propagates the existing norm error. The resulting early word has
the required one-external-clean full-isometry contract.

### Simultaneous count and dirty-width margin

By (26), all relevant large table volumes, including all selected bins,
obey

```math
\mathcal Q\le p\,2^{h+m_0},\qquad
2^{h+m_0}=N2^{-(M-m_0)}.
\qquad(29)
```

Short source indices, calibration precisions, masks, outputs, reflected
predicates, repeated amplification calls, and simultaneous private pools
contribute fixed polynomial factors. The outer descendant sources and
coarse reflections also fit this envelope: their relevant early-data
address has at most h bits, with polynomially many tags and source-index
values. Thus there is a fixed degree c, chosen before $`A_L`$, such
that the entire early circuit has simultaneous bounds

```math
T_E,w_E\le C_L(n+2)^c\left(2^{(h+m_0)/2}+1\right),
```

```math
G_E\le C_L(n+2)^c\left(2^{h+m_0}+1\right),
\qquad D_E=O_L(\ell^{3/2}).
\qquad(30)
```

The width $`w_E`$ includes the fixed external data/work banks, routing
padding, and all parallel selected-query helpers. Independent simultaneous
pools are summed; fully returned pools are reused between calls. The
early depth consists of $`O_L(\ell\sqrt\tau)`$ for the local word,
$`O_L(\ell\log(\ell+2))`$ for coarse recovery, and
$`O_L(\ell\log\log(\ell+4))`$ for boundary routing, all within
(30).

Choose $`A_L>c_L+2c+20`$, increasing it if necessary for the separately
fixed metadata reservation. Equations (29)–(30) then imply

```math
T_E=O_L(\sqrt N),\qquad G_E=O_L(N),\qquad
w_E=o(\sqrt N/n).
\qquad(31)
```

The independent count margin $`M-m_0`$ is essential; the same wires
have not both enlarged the table and been counted as untouched suppression.
No efficient preprocessing or practical crossover follows from this
asymptotic choice.

## 10. Terminal factor, widths, signs, and prescribed phases

Use the wide construction only beyond a fixed L-dependent dimension
cutoff and when $`b\ge\sqrt N/n`$. Equation (31) makes the entire
early dirty peak fit in a fixed fraction of b. All polynomial metadata,
routing, predicate, and external-bank reservations fit there as well.
The early circuit does not need a width tradeoff inside its large tables
in this branch.

The exact early/terminal factorization is applied to the same rounded
positive coefficient tuple throughout. Its terminal height is
$`M=\lceil A_L\ell\rceil`$. The
[uniform terminal theorem R81](RECURSIVE_TERMINAL_DEPTH.md) gives, at
error $`\eta/4`$, one clean qubit and the simultaneous bounds

```math
T_T=O_L(\sqrt N+N/b),\qquad G_T=O_L(N),
```

```math
D_T^{\rm terminal}
=O_L\!\left(N/b^2+\ell\log(M+2)\right)
=O_L\!\left(N/b^2+\ell\log(\ell+2)\right).
\qquad(32)
```

The terminal's stronger internal constant-accuracy reservation fits
b automatically beyond the selected cutoff. It reuses the returned
arbitrary-dirty pool. Sequential composition uses initialized-isometry
hybrids and actual unitaries, so no reset or perfect actual intermediate
clean return is assumed.

The two finite sign gauges in (1b) are exact Boolean phase tables. Their
retained sharp/routed dirty-query implementations have
$`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and
$`D_T=O_L(N/b^2+\ell)`$, with exact arbitrary-dirty return. In the
prescribed complex magnitude family, retain the given diagonal
$`D_\phi`$ at every leaf, including leaves of original zero state
amplitude. It multiplies the entire frame, so dropping those phases
would change other columns. The existing one-clean diagonal compiler
implements it at error $`\eta/4`$ with the same count envelope and
$`O_L(N/b^2+\ell)`$ depth. This is the prescribed phase-dressed
completion, not an arbitrary complex-unitary synthesis theorem.

Rounding, early compilation, terminal compilation, and the optional
complex phase contribute at most

```math
\eta/32+\eta/8+\eta/4+\eta/4\lt \eta.
\qquad(33)
```

The exact sign gauges preserve these errors. Combining (30)–(33) proves
the displayed simultaneous bounds in the wide branch, with all dirty
systems and external references included.

If $`b\lt \sqrt N/n`$, use the complete original
[R72 one-clean compiler](ONE_CLEAN_DEPTH.md) at the original total
accuracy eta, rather than at a stricter internal allocation. Its depth
is $`O_L(N/b^2+nL)=O_L(N/b^2)`$, since $`N/b^2>n^2`$. It has
the required T-count and Clifford-count bounds and preserves the literal
reservations $`17(L+n+7)`$ and $`17(L+n+8)`$. For the finitely many
dimensions below the chosen cutoff, use that same original-accuracy
compiler. Its angle-uniform bounds are absorbed by L-dependent constants.
No compactness argument, exact sign decision, or new small-width
assumption is needed. This completes the theorem at every stated width.

The result is a general upper bound at fixed L. It does not establish
absolute all-precision constants, logarithmic optimality, or a growing
lower bound. Those questions remain separate from this strict improvement
below square-logarithmic depth.

## 11. Evidence and remaining scope

[The bounded tests](../../tests/test_interval_dilation_depth.py) check
heterogeneous integer interval embedding, every prescribed local scaling
flag and marker sign, rectangular overlap and an isolated root block,
classical membership erasure and its actual inverse, square compression
with all rejection fields retained, calibrated amplification with literal
phase and source-work leakage, and fine-tag erasure with an external
reference. They use small matrix completions to check the interfaces;
those completions are not native circuit emitters. A finite norm
diagnostic does not prove (10), and numerical work is not used to infer
the asymptotic resource bounds.

The analytic native ingredients are the existing family-source compiler,
the charged selected-query/routing and active-anchor lemmas, the balanced
frame word, and the actual-inverse amplification above. The
[session record](GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md) distinguishes
this result from the stopped row-label and sequential relative-transport
attempts. Scientific receipts and the existing selected count theorems
are unchanged. The continuing destination includes sharper general depth;
this construction closes the present sub-square-logarithmic target.
