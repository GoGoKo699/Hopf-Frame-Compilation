<a id="a-growing-source-compiler-for-complete-real-hopf-frames"></a>

# A growing-source compiler for complete Hopf frames

[Mixed local gauges](HETEROGENEOUS_PRECISION_SOURCE.md#2-mixed-local-merges-with-one-feature-per-node) · [Exact dirty scalars](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger) · [Disjoint state extensions](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)

The carried source can grow as the bottom-up row network retires its
detail and feature outputs. Fresh mode arenas supply the larger source
and all private fields. Every rejected mode is retained permanently,
and every later completed word fixes it. An exact dyadic refinement
changes the source ratio without another endpoint-precision synthesis.
The resulting precision charges are summable. A single final column
source compensates the known attenuation of each retired output.
The real construction leaves a half-error margin; the
[complex-magnitude corollary](#9-complex-magnitude-endpoint) uses it to
append a literal diagonal within the same workspace and linear T-count.
The [precision-uniform extension](#10-precision-uniform-two-clean-count)
separates the logical dimension from the requested accuracy and gives
$`O(N+L)`$ count with two clean qubits.
The [one-clean extension](#11-precision-uniform-one-clean-count)
borrows only the fine-phase signal, retaining the initialized recovery
anchor. The [banked frontier](#12-optimal-banked-count-frontier)
then matches the worst-case count lower bound at its stated dirty widths.

## 1. Statement, clipping, and the fixed anchor

Let $`n\ge3`$, $`N=2^n`$, and $`\epsilon=2^{-N}`$. Let W be a
prescribed complete real Hopf frame whose rotation data admit terminating
certified evaluation. No regularity, sparsity, or exponent promise is
imposed.

**Theorem.** There is an effective coherent Clifford+T circuit V with

```math
T(V)=O(N),\qquad
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-N},\qquad b=N+n+7.
\qquad\text{(1)}
```

The constant is universal. The norm includes every logical column,
both initialized outputs, arbitrary dirty inputs, and external
references. All inverses are actual circuit inverses. Classical
coefficient evaluation and finite table construction are preprocessing;
there are no supplied precision states, uncharged coherent evaluators,
measurements, or resets.

Use the [effective normalization and clipping convention](HETEROGENEOUS_PRECISION_SOURCE.md#1-family-contract-and-the-upper-sector).
After an optional separately charged signed normalization, the retained
normalized branches have selected integer certificates

```math
2^{-E_v-2}\le s_v\le2^{-E_v},\qquad
0\le E_v\le N+\lceil\log_2n\rceil+5.
```

The clipping rule uses a certified dyadic approximation and overlapping
tests, not a zero oracle or a floor of a real logarithm. It changes the
complete frame by at most $`\epsilon/16`$. At a clipped or known zero
node use the exact identity chart. The original signed permutations are
priced in $`O(N)`$. Reserve another $`\epsilon/16`$ for any rounding
needed by effective normalization.

Choose one sufficiently large universal integer $`d_0`$. Its size is
fixed by the mode-volume bound in Section 4, not by n or the target.
The working construction is for n above a universal threshold with
$`n\gt d_0`$; Section 8 treats every smaller height by a finite
fallback before any $`d_0`$-bit suffix is used.
Separate the bottom $`d_0`$ physical binary levels of the clipped frame:

```math
\widetilde W=W_{\rm bottom}W_{\rm upper}.
```

The bottom part is a direct sum of $`N/2^{d_0}`$ blocks in
$`\mathrm{SU}(2^{d_0})`$, with sum of squared block dimensions
$`2^{d_0}N`$. The [disjoint-block compiler](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks),
at its fixed constants $`\kappa=2^{d_0}`$ and $`c=4`$, therefore gives
error at most $`\epsilon/4`$ and $`O(N)`$ T-count with the same
two clean qubits and $`N+n+7`$ dirty pool. The forward
frame circuit executes the upper compiler before this bottom compiler.
Only the row-source construction inside the upper compiler runs
bottom-up.

The upper frame acts on the original anchor set
$`\mathcal A=\{x:\text{the final }d_0\text{ logical bits are zero}\}`$,
of size $`N_0=N/2^{d_0}`$, and is identity on its complement.
Compute $`h=[x\in\mathcal A]`$ in external clean qubit 1 with the exact
borrowed predicate word. Preserve h throughout. External clean qubit 2,
denoted c, is the only fine-diagonal compiler flag. All other source,
feature, attenuation, scalar, and arithmetic fields below are coordinates
of the original n logical mode bits. They are not extra qubits.

Every completed upper subroutine is exactly identity on $`h=0`$,
including on occupied private fields and arbitrary dirty inputs. The
[conditional-suffix interface](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally)
then permits the original anchor suffix to supply the initial work.
After recovery the ideal output again lies in $`\mathcal A`$ and has $`c=0`$.
Reversing the actual predicate computation returns h; the full-isometry
error also bounds the actual suffix and flag leakage.

## 2. A polynomial schedule and the mixed local gauge

The physical height d is measured upward from the original leaves.
Initial upper carriers are the roots at height $`d_0`$. For
$`d_0\le d\le n`$ set

```math
\begin{gathered}
P_d=2^{2\lceil\log_2(d+1)\rceil},\qquad
M_d=4P_d^2,\qquad q_d=N/P_d^2,\qquad\beta_d=2^{-q_d},\\
p_d=\frac1{8P_d},\qquad B_d=6+2\log_2P_d,\qquad
\Delta_d=q_d+B_d.
\end{gathered}
\qquad\text{(2)}
```

Thus $`(d+1)^2\le P_d\lt4(d+1)^2`$. The sequence P changes only
by a factor four. For all n beyond a universal threshold, every
$`q_d`$ is a positive integer: $`P_d^2`$ is a power of two bounded
by a fixed polynomial in n, while $`N=2^n`$. Enlarge that threshold
so that $`q_d\ge B_d`$ and every source ratio used below is at most
the corresponding $`p_d`$. All smaller n have a finite fallback in
Section 8. No logical dimension is padded.

At height d call a retained node cheap if $`E_v\lt4N/P_d`$ and
expensive otherwise. Use local relative tolerance and splitter loss

```math
\eta_d=\frac1{1000P_d^2},\qquad
l_d\in[1/(1600P_d),1/(800P_d)].
\qquad\text{(3)}
```

The [one-feature mixed construction](HETEROGENEOUS_PRECISION_SOURCE.md#2-mixed-local-merges-with-one-feature-per-node)
supplies at each node the exact physical local column

```math
\begin{bmatrix}B_v^\dagger\\F_vL_v^\dagger\end{bmatrix},
\qquad B_vB_v^\dagger+l_vL_vL_v^\dagger=I,
\qquad F_v^\dagger F_v=l_vI.
\qquad\text{(4)}
```

At a cheap node this is the actual native lossy scalar word with its
actual rejected operator on the dirty core. Fix one precision
$`q_{{\rm ch},d}=\lceil4N/P_d\rceil+O(\log P_d)`$ at that
height, and fix the emitted masks and physical core before computing
the gauge. Source-index bits and future arena coordinates are never
borrowed as that core: on an enabled node the word is exactly
$`I_{\rm source}\otimes U_v`$. The dyadic loss parameter can be
$`2^{-E_v-\lceil\log_2(400/\eta_d)\rceil}`$.
At an expensive node the mathematical coarse merge is the Cayley
rotation $`R(t_v)`$, with

```math
t_v=2^{-E_v}u_v,\qquad 1/16\le u_v\le1,
\qquad m_d=O(\log P_d)
```

dyadic mantissa bits and relative entry error at most $`\eta_d`$.
For example $`2^{-m_d}\le1/(16000P_d^2)`$ suffices. Its row
implementation is $`R(t_v)^\dagger`$ followed by one actual native
splitter with literal first column $`(r_d,t_d^{\rm sp})`$ and loss
(3). Consequently
$`B_v=R(t_v)\mathrm{diag}(\bar r_d,1)`$,
$`L_v=R(t_v)e_0`$, and $`F_v=t_d^{\rm sp}I`$.
No complex phase is discarded. A zero node has $`B_v=I,F_v=0`$ and
no splitter or scalar query.

Apply the same local two-coordinate gauge identity from that chapter
to obtain the coarse matrices $`W_B,L_B`$, leaf diagonal D, and
column coefficients C,A. It gives on every upper-frame column

```math
\begin{gathered}
W_{\rm upper}=D(W_BC+L_BA),\\
Q_0x=\begin{bmatrix}W_B^\dagger x\\
 (F_v(L_B^\dagger x)_v)_v\end{bmatrix},\qquad
Q_0^\dagger Q_0=I.
\end{gathered}
\qquad\text{(5)}
```

Here $`Q_0`$ denotes the mathematical mixed row network before source
attenuation. Its accepted carrier blocks are scalar on the dirty space;
each feature leaves the carrier network once and is fixed afterward.
Thus its feature contains exactly its own $`F_v`$, with no product
of rejected dirty operators from different nodes.

The variable tolerances retain the previous uniform gauge slack. Along
any path,

```math
\sum_{d>d_0}\eta_d\le\frac1{3000d_0^3},\qquad
\sum_{d>d_0}l_d\le\frac1{800d_0}.
```

Increasing fixed $`d_0`$ therefore ensures
$`.999\lt|D|\lt1.002`$ and $`.98\lt|C_v|\lt1.02`$.
The cheap bound is $`|A_v|/l_v\lt.8`$; the expensive bound is
$`|A_v|/l_d\le4.8/P_d\le.3`$. These statements retain the
complex local phases in the gauge. Root and zero-node columns have
zero feature coefficient, without evaluating a zero denominator.

## 3. Exact source refinement and carried kernels

For $`0\lt\beta\lt1`$ define the normalized finite source

```math
g_{\beta,M}=a_{\beta,M}\sum_{j=0}^{M-1}\beta^j|j\rangle,
\qquad a_{\beta,M}=\sqrt{\frac{1-\beta^2}{1-\beta^{2M}}}.
```

Prepare $`g_{\beta_{d_0},M_{d_0}}`$ once on every initial carrier,
using the one-clean disjoint state-extension interface with compiler
flag c. This has $`O(N)`$ cost at any fixed error fraction of
$`\epsilon`$. Its supports and padding fit the mode allocation in
Section 4. Its emitted native inverse is the only inverse ever used.
No intermediate source is prepared again or unprepared.

### A source square root with an exact dyadic scalar

Appending a zero coordinate bit and interleaving
$`|j\rangle|b\rangle\mapsto|2j+b\rangle`$ gives the exact identity

```math
(I+\beta X)(g_{\beta^2,M}\otimes|0\rangle)
=\sqrt{1+\beta^2}\,g_{\beta,2M}.
\qquad\text{(6)}
```

For $`p=2^{-r}`$ with positive integer r and dyadic
$`\beta\le p`$, set
$`f=(1-p)\beta/p`$. Prepare an exact rare event of probability p,
apply the exact raw dyadic scalar f on that event with a separate
retained scalar flag, apply X on the appended source bit on that same
rare event, and reverse
the rare-event preparation. Its zero-private compression is

```math
T_{\rm ref}=(1-p)I+pfX=(1-p)(I+\beta X).
```

It is a contraction because $`(1-p)(1+\beta)\le1`$. On the
intended normalized source it has the scalar attenuation

```math
\lambda_{\rm ref}=(1-p)\sqrt{1+\beta^2},\qquad
1-p\le\lambda_{\rm ref}\le1.
\qquad\text{(7)}
```

The square root in (7) is not omitted. With $`\beta=2^{-q}`$ and
$`p=p_d`$, f is the exact dyadic $`(8P_d-1)2^{-q}`$, requiring
only $`q+O(1)`$ scalar precision. The source scalar has accepted
action fI on every dirty input; its actual rejected operator is kept.

To process nodes at height $`d>d_0`$, first refine all still-live
child carriers from $`g_{d-1}`$ to $`g_d`$. There is no refinement
when P is unchanged. At a jump $`P_d=4P_{d-1}`$, use four binary
substeps, each halving q and doubling M, with rare probability $`p_d`$.
Their product changes $`q_{d-1}`$ to $`q_d=q_{d-1}/16`$ and
$`M_{d-1}`$ to $`M_d=16M_{d-1}`$. Every substep retains
$`qM=4N`$. A new zero coordinate comes from an actual fresh mode
embedding, not an assumed reset.

### The expensive merge kernel

For an expensive node at height d compute offline

```math
E_v=k_v\Delta_d+e_v,\qquad 0\le e_v\lt\Delta_d.
```

For sufficiently large n, $`q_d\ge B_d`$ and the expensive cutoff
imply $`k_v\ge2P_d`$. Clipping gives $`k_v\le2P_d^2`$.
Only this $`O(\log d)`$-bit multiplier is loaded. The remainder,
which can require order n classical bits, stays entirely in the
offline scalar masks. For $`1\le\nu\le P_d`$ define

```math
f_{v,\nu}=2P_d(8P_d-1)u_v^\nu
                 2^{-\nu(k_vB_d+e_v)}.
\qquad\text{(8)}
```

These are exact dyadics in $`[0,1/4]`$. Set them to zero on cheap,
zero, and frozen correction rows. Their scalar precision is bounded by

```math
q_{{\rm coef},d}
\le P_d(m_d+2P_d^2B_d+\Delta_d)+O(\log P_d)
=\frac N{P_d}+O(P_d^3\log P_d).
\qquad\text{(9)}
```

Let $`S_k=\sum_{j=k}^{M_d-1}|j-k\rangle\langle j|`$ when
$`k\lt M_d`$, and $`S_k=0`$ otherwise. Its native dilation first
compares the full product $`k_v\nu`$ with the source index, toggles
a distinct underflow flag, and only then subtracts modulo $`M_d`$.
Thus a shift at least $`M_d`$ has accepted block zero. Metadata and
arithmetic are uncomputed before the logical quarter-turn changes the
carrier endpoint. The quarter-turn is $`J_v^\nu`$, not
$`J_v^{k_v\nu}`$, with $`J_v=-\begin{bmatrix}0&-1\\1&0\end{bmatrix}`$
on the carrier pair. Source-underflow and scalar-rejection flags are
distinct and retained.

Use an exact uniform $`P_d`$-valued power label and a rare event of
probability $`p_d`$. The zero-private block of the full unitary kernel is

```math
T_d=(1-p_d)I+\frac{p_d}{P_d}
\sum_{\nu=1}^{P_d}\sum_{v\text{ at height }d}
S_{k_v\nu}\otimes f_{v,\nu}J_v^\nu.
\qquad\text{(10)}
```

On the common source, use
$`t_v=\beta_d^{k_v}u_v2^{-(k_vB_d+e_v)}`$. The finite-shift
residual is at most $`2\beta_d^{M_d}`$ for every shift, including
shifts past the cap. The same direct-sum argument as in the
[fixed-source kernel](HETEROGENEOUS_PRECISION_SOURCE.md#4-exact-dyadic-local-kernels)
therefore gives

```math
\begin{gathered}
\|T_dE_{g_d}-E_{g_d}(1-p_d)\mathcal R_d\|
\le\delta_d,\\
\delta_d=2p_d\beta_d^{M_d}
 +\frac{2(1-p_d)t_{\max,d}^{P_d+1}}{1-t_{\max,d}}
\le5\,2^{-4N}.
\end{gathered}
\qquad\text{(11)}
```

Here $`\mathcal R_d`$ is the direct sum of expensive Cayley merge
adjoints and identities on the other pairs. Indeed
$`\beta_d^{M_d}=2^{-4N}`$ and
$`t_{\max,d}^{P_d+1}\le2^{-4N}`$. Disjoint pairs give a maximum,
not a sum over nodes. The intertwiner $`E_{g_d}`$ tensors the same
source with identity on the complete system and arbitrary dirty
register. It thus commutes with the actual dirty-dependent cheap
interlayers, including on occupied feature modes.
Every surviving carrier, including one at a cheap or zero node,
receives the same baseline $`1-p_d`$ and the same required source
refinements. Zero corrections do not skip that baseline.

## 4. Full-unitary mode arenas and exact guards

### Fixed external flags and packed compact banks

On $`h=1`$, first apply an exact permutation $`\pi_0`$ sending the
strided original anchor set $`\mathcal A`$ to a compact
$`N_0`$-mode input bank.
After this packing take a,b as two logical mode coordinates. Reserve
four compact $`N_0`$-mode subbanks: ordinary representatives in
$`a0b0`$, feature representatives in $`a0b1`$, column attenuation
failure in $`a1b0`$, and row attenuation failure in $`a1b1`$.
Unused root-feature positions are harmless padding. Every growing
arena and every final shared source tag lies in $`a=0`$. The coordinate
b need not keep its feature interpretation inside those arenas.

Both source words include the same initial packing:
$`Q=\widetilde Q\pi_0`$ and $`P=\widetilde P\pi_0`$.
Their actual inverses undo it. Consequently the overlap is conjugated
back to the original anchors, and the recovery reflection below tests
the original suffix, not the packed bank or later arena fields.

### A summable full-space reservation

At height d there are at most $`2N/2^d`$ entering carriers. A complete
kernel arena includes the node/port label, all $`M_d=O(d^4)`$ source
coordinates, the power and rare labels, multiplier and full product,
scalar and underflow flags, node type, and reversible arithmetic work.
The multiplier, product, and scratch have $`O(\log d)`$ bits. No
mantissa or remainder register is initialized. Standard reversible
addition, comparison, and multiplication use $`O(\log d)`$ work
bits and polynomially many gates, returning the scratch before the
stage boundary. A refinement has no larger type of reservation.

It follows directly that the complete arenas at height d, including
every private-field value and a constant number of refinements, have
total size at most

```math
C_1N2^{-d}(d+1)^D
\qquad\text{(12)}
```

for universal constants $`C_1,D`$. Round every reservation to a
dyadic size; absorb that factor and the padding for final tag
extensions into $`C_1,D`$. The bound counts whole unitary dilations,
not only initialized source columns or accepted modes. Choose fixed
$`d_0`$ so large that

```math
4\,2^{-d_0}+C_1\sum_{d\ge d_0}2^{-d}(d+1)^D\lt1/4.
\qquad\text{(13)}
```

The initial source arena is included at $`d=d_0`$. All compact
representatives, complete stage arenas, sealed rejections, retired
outputs, and conservative tag padding now fit inside N logical modes.
No dirty qubit is used as one of these mode coordinates.

An explicit aligned placement avoids an implicit fragmentation
assumption. Put the two compact shared subbanks at the starts of the
two $`a0b`$ quarters. The upper half of each such quarter has size
$`N/8`$ and is free. Every
dyadic arena is at most $`N/8`$, and their total is below $`N/4`$;
descending-size dyadic packing places them in those two free halves.
The two compact attenuation-failure banks lie in their respective
$`a1b`$ quarters. All arena predicates are unchanged high-prefix
conditions on their internal source and private operations.

### Fresh work is supplied by a permutation

Before a stage, distinguish live source cells L, frozen shared cells,
and sealed row failures. Choose a fresh arena X disjoint from every
previous reservation and the compact failure banks. Inject every basis
label in L into X with its source coordinate unchanged and every new
coordinate/private bit zero. For a refinement this includes the new
source bit. If the injection is f, the disjoint-swap permutation

```math
\pi_X=\prod_{j\in L}(j,f(j))
\qquad\text{(14)}
```

is a full-space unitary. It fixes all old frozen and sealed cells and
supplies the zero fields on the entire live subspace, including dirty
and reference entanglement. It does not reset a register or assume
anything about arbitrary inputs elsewhere in the mode space.

Apply the full refinement/kernel word within X and identity elsewhere.
Seal every output with nonzero acceptance-private fields. Guarded on
those fields being zero, apply the exact cheap merges and native splitters,
freeze each new detail and feature source block, and carry only the
parent cells onward. Future permutations and completed words fix all
old frozen and sealed cells, on their entire source and dirty spaces.
Old evacuated live cells remain reserved; none is silently recycled as
fresh work. Each operation is a unitary on every physical input.
Any fourth port padding of a three-port merge, or other non-designated
arithmetic and private output coordinate, also remains permanently
reserved. The assertion that only live, detail, and feature columns
occur concerns the ideal initialized input, not arbitrary occupied
padding inputs.

The ideal initialized trajectory has no amplitude in a fresh arena
before (14). Actual fine-preparation leakage may have such amplitude;
the complete-isometry hybrid propagates it through the same full
unitaries. It is never declared zero or discarded.

### A zero coefficient does not provide an arena guard

The full-word guard must freeze an old arena even when all its private
fields are occupied. Let e be the preserved conjunction of h, the
shared sector, and the current aligned arena prefix. Write an exact
dirty source as $`M=U_mX_0U_m^\dagger`$ and its unitary mask as
$`P_f`$. The raw scalar word has the form

```math
\mathcal S_f=H_z C_{z=0}(M)P_fMP_f^\dagger C_{z=1}(M)H_z.
```

Replace its three source occurrences by
$`C_{e\wedge z=0}(M)`$, $`C_e(M)`$, and
$`C_{e\wedge z=1}(M)`$, respectively. On $`e=1`$ this is the
original scalar word. On $`e=0`$ it is exactly
$`H_zP_fP_f^\dagger H_z=I`$ on every signal and dirty input.
The second mask is its actual inverse. Real-kernel masks are
involutions, but paired masks used by cheap complex scalar constituents
need not square to I; their literal phases are covered by this form.
Each guarded
source is implemented as

```math
C_G(M)=U_mC_G(X_0)U_m^\dagger.
\qquad\text{(15)}
```

Only its central X receives the long prefix control; the long source
loaders are not controlled gate by gate. Borrowed predicate helpers
return before the inverse loader. This is the full-space
[central-source guard](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger),
and it applies to every raw constituent of the cheap scalar word as
well. Assigning f=0 alone would not prove this identity.

Short metadata queries can also retain their compressed addresses.
Let Q be an involutory XOR query, let z be a dirty control disjoint
from its targets, and let $`T_e`$ toggle z by e. The chronological echo
$`T_e,C_zQ,T_e,C_zQ`$ implements $`Q^e`$ and returns z, since
$`(z\mathbin\oplus e)\mathbin\oplus z=e`$. The controlled query
has only one extra short address bit. Each completed query returns
borrowed selectors before dirty-changing work; its inverse works on
the possibly changed arbitrary dirty state.

Coefficient masks read only the node and power labels, not the source
index or all ambient N mode labels. Thus the guard does not expand the
table. Rare/power Hadamards, arithmetic, and quarter-turn SELECTs have
only polynomially many small-register gates and use the same exact
borrowed-prefix conditioning. All inactive private inputs are covered.

The constant fine calls obey the same whole-word condition. Use the
completed conditional one-clean primitive for each literal diagonal
and attenuation bank, including its conditional leading sign when
it contains amplification. Condition the state-extension Walsh,
coarse, center, and routing interlayers by the exact borrowed-control
interface. All these completed words are identity on $`h=0`$ even
when c is occupied. An identity phase table alone is not used as a
substitute for that full-word guard. Preserved block predicates also
fix the appropriate compact failure banks during later source steps.

## 5. Retired source tags and the compensated overlap

Include frozen successful outputs and current live cells in the
accepted direct sum. A new stage's accepted compression is a
contraction; it is identity on all earlier frozen coordinates.
Equation (6) is exact on a live source column, and (11) bounds each
kernel's error on that column, on the entire dirty space. Therefore
contractive telescoping gives

```math
e_d\le e_{d-1}+\delta_d,\qquad
e_n\le5n\,2^{-4N}.
\qquad\text{(16)}
```

Source-cap residuals remain in their accepted source coordinates.
Neither a live source nor a failure register is asserted to have
returned at an intermediate boundary.

For each final output r let $`\Gamma_r`$ be the product of every
kernel factor $`1-p_d`$ and every normalized refinement factor (7)
on its live path until it freezes. A node's detail and feature have
the same factor $`\Gamma_d`$. The final root has the product along
the full path. At most one kernel and four refinement events occur at
each height, so

```math
\sum_{\rm events}p\le\frac5{8d_0},\qquad
7/8\le\Gamma_r\le1
\qquad\text{(17)}
```

after increasing fixed $`d_0`$ if necessary. The upper bound uses
$`\lambda_{\rm ref}\le1`$, not a discarded source normalization.
All these scalar factors are known by terminating classical evaluation.

Let $`J_{\rm tag}`$ send each compact ordinary/feature representative
to the normalized source $`g_d`$ in its final frozen cell, and the
root to $`g_n`$. Those target cells are pairwise disjoint. Let
$`\Gamma`$ be their positive diagonal attenuation. The ideal shared
part of the row source, after attenuation by $`D^\dagger/\Lambda`$
with $`\Lambda=5/4`$, is

```math
Q_{\rm shared}=\frac1\Lambda J_{\rm tag}\Gamma
\begin{bmatrix}W_B^\dagger D^\dagger\\
(F_v(L_B^\dagger D^\dagger)_v)_v\end{bmatrix}.
\qquad\text{(18)}
```

Equation (18) is a reference shared column. The actual source is
obtained from the native stream, not from an uncharged completion of
that reference. Its accepted-column error is (16), plus
its charged fine-preparation errors. Its other amplitudes remain in
the sealed row-failure cells or the row attenuation-failure bank.
There is no intermediate amplification and no final row-source
unpreparation.

### The column source repeats the actual dirty rejection

Set $`\kappa=8/5`$. At the compact detail representative for node v,
prepare shared coefficients

```math
x_v=\frac{C_v}{\kappa\Gamma_d},\qquad
y_v=\frac{A_v}{\kappa\Gamma_d l_v},\qquad
\rho_v=\sqrt{|x_v|^2+l_v|y_v|^2}.
\qquad\text{(19)}
```

The unscaled mixed construction has $`\rho_v\lt.7`$; (17) gives
$`\rho_v\lt.8`$ here. The normalized replay ratio
$`|y_v|/\rho_v`$ is unchanged by compensation. At root and zero
nodes use zero feature coefficient and no division by l.

For the actual local word $`U_v|0\rangle=[z_{0,v}I;F_v]`$, use
the exact identity

```math
U_v\mathrm{diag}(e^{i\phi},1)U_v^\dagger|0\rangle
=\begin{bmatrix}
1+(e^{i\phi}-1)|z_{0,v}|^2\\
(e^{i\phi}-1)\bar z_{0,v}F_v
\end{bmatrix}.
\qquad\text{(20)}
```

It prepares the required normalized pair after two sector phases.
For cheap nodes replay the exact emitted row word, with the same
physical dirty core and finite masks at that height. Relabeling its
node address to its unique compact detail label does not replace it
by another synthesis with the same accepted scalar. Expensive nodes
repeat their actual native splitter. All these words have disjoint
compact detail labels. Execute all blocks of $`U_*^\dagger`$, one
global fine phase bank, all blocks of $`U_*`$, and the two global
sector-phase banks. There is no fine phase call inside the depth loop.
Completed routing and guards may differ, but on the active node they
give that exact $`U_v`$ and identity on returned helpers and selectors.
Only logical modes are relabeled. Thus the rejected $`F_v`$ is
identical even when earlier row words have changed the arbitrary dirty
state. Different heights may use different core sizes; their actual
words and cores are retained in the replay.

First attenuate by $`\rho_v`$ and put failure in $`a1b0`$.
Every subsequent replay, phase, and tag operation is enabled only in
the shared sector and fixes that failure bank. The fine compiler flag
is the independent external c, never an occupied attenuation mode.
Choose replay phases by certified amplitude-norm approximants, as in
the mixed-source construction, including when a requested component
vanishes. No exact test whether a nonzero node's correction $`A_v`$
is zero is needed.

After this compact replay, apply one disjoint state extension realizing
$`J_{\rm tag}`$. Its representative $`e_r`$ lies outside its target
tag cell, so an exact extension is the product of the two reflections
with normals $`(e_r+g_d)/\sqrt2`$ and $`e_r`$. Supports for distinct
r are disjoint. Their padded sizes, at most $`2(M_d+1)`$ each, were
included in (13). Invoke the full disjoint-state interface, including
its priced coarse preparation; the first normal is not presumed close
to $`e_r`$. This costs $`O(N)`$ with a constant number of fine calls
and ideally fixes the column-failure bank on its whole space.

The resulting ideal shared column is

```math
P_{\rm shared}=\frac1\kappa J_{\rm tag}\Gamma^{-1}
\begin{bmatrix}C\\(A_vF_v/l_v)_v\end{bmatrix}.
\qquad\text{(21)}
```

The row's sealed failures and attenuation failure are disjoint from
the column attenuation failure. The ideal column has no amplitude in
any sealed row-failure cell. Equations (18)--(21), with
$`F_v^\dagger F_v=l_vI`$, give the literal complete-column overlap

```math
Q_{\rm shared}^\dagger P_{\rm shared}
=\frac1{\Lambda\kappa}D(W_BC+L_BA)
=W_{\rm upper}/2.
\qquad\text{(22)}
```

The tags and their positive attenuations cancel. Distinct node features
exclude cross-node dirty products. No information about the size or
polar part of a sealed rejected operator is used.

## 6. Actual inverses and the original-anchor recovery

Let J initialize only the original suffix $`d_0`$ and external c on
the active $`h=1`$ sector. It is identity on the original reduced
logical input and on all dirty inputs. It does not test fields that
are only later interpreted as source, metadata, or private coordinates.
After the initial packing has been undone by the actual source
inverses, (22) is expressed on precisely this original anchor.

Let $`\widehat P,\widehat Q`$ be the full emitted source words and
$`U=\widehat Q^\dagger\widehat P`$. Fine-call replacements are
bounded by a complete-isometry hybrid; earlier leakage continues
unitarily. Compare the overlap directly as
$`J^\dagger UJ=(\widehat QJ)^\dagger(\widehat PJ)`$.
Both source columns are independently initialized. The ideal column
source is orthogonal to every sealed row failure, so only the row's
shared projection requires (16). No approximation for
$`\widehat Q^\dagger`$ on the arbitrary output of
$`\widehat P`$ is assumed. From these column bounds,

```math
\|J^\dagger UJ-W_{\rm upper}/2\|
\le e_n+\epsilon_P+\epsilon_Q=:\delta.
\qquad\text{(23)}
```

The source budgets $`\epsilon_P,\epsilon_Q`$ include certified
approximation of attenuation coefficients, tag states, local-chart
parameters and replay phases, as well as native errors of every fine
call and its actual inverse. The
[disjoint-state chart](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)
uses amplitude-norm approximation, including at zero components.
Allocate half of each source budget to certified data and half to native
realization. The exact mixed gauge already corrects the actual cheap
scalars, mantissas and splitter losses; these are not additional final
frame errors.

Only a constant number of fine calls occurs in each source. Choose
$`q=N+c_0`$ for one fixed sufficiently large $`c_0`$ so that the
total fine hybrid errors satisfy
$`\epsilon_P,\epsilon_Q\le\epsilon/384`$. Beyond the universal
threshold, (16) is at most $`\epsilon/192`$, giving
$`\delta\le\epsilon/96`$.

For $`R=I-2JJ^\dagger`$ the actual word $`-URU^\dagger RU`$
has accepted and rejected columns

```math
3A_U-4A_UA_U^\dagger A_U,\qquad
B_f(I-4A_U^\dagger A_U),
```

where $`A_U=J^\dagger UJ`$ and $`B_f=(I-JJ^\dagger)UJ`$.
The [normalization-two lemma](../../docs/FAULT_TOLERANT_COMPILER.md#4-from-a-half-unitary-block-to-the-complete-output)
bounds the full output error by $`12\delta\le\epsilon/8`$.
All displayed identities hold on arbitrary dirty and reference inputs;
the reflection excludes them.

Use the literal conditional word

```math
Z_hU_hR_hU_h^\dagger R_hU_h,\qquad
R_h=I-2\Pi_{h=1,\,x\in\mathcal A,\,c=0}.
\qquad\text{(24)}
```

Its leading sign is correct on $`h=1`$. Every factor is exactly
identity on $`h=0`$. The zero-work test is the original suffix and c;
it never tests all of the later growing arena coordinates. Reversing
the actual h predicate then returns both external qubits. Compose the
bottom compiler at error at most $`\epsilon/4`$. Together with
clipping and optional normalization, the complete error is at most

```math
\frac{\epsilon}{16}+\frac{\epsilon}{16}
+\frac{\epsilon}{8}+\frac{\epsilon}{4}
=\frac{\epsilon}{2}.
```

Exact predicates, routing and their actual reversals add no error.
This proves (1) with a half-error margin and the same literal workspace.

## 7. Summable native charges

Write $`S_d=N/2^d`$ for the number of nodes at physical height d.
At one kernel the coefficient table has $`O(S_dP_d)`$ rows; its
values are identical across source indices and unused private fields.
The raw scalar and exact guards cost

```math
O(S_dP_d+q_{{\rm coef},d}+\mathrm{poly}(n)).
```

Short multiplier/type tables have no larger row charge. Their small
outputs do not supply their own selector bank: completed lookups borrow
and return $`O(n)`$ selectors from the much larger current scalar
core. This happens before dirty-changing work, and inverse lookups are
valid on the resulting arbitrary dirty state.

From $`P_d\asymp(d+1)^2`$ and (9),

```math
\begin{gathered}
\sum_{d>d_0}S_dP_d=O(N),\\
\sum_{d>d_0}q_{{\rm coef},d}
\le N\sum_{d>d_0}P_d^{-1}+\mathrm{poly}(n)=O(N).
\end{gathered}
\qquad\text{(25)}
```

The additive polynomial is absorbed beyond a universal threshold;
the precision term is retained explicitly. Cheap row words and their
exact column replays similarly cost

```math
O\left(\sum_{d>d_0}[S_d+4N/P_d+O(\log P_d)]
       +\mathrm{poly}(n)\right)=O(N).
```

The common native splitter at a height has polynomial-in-log(d) word
length. Load the preserved node-type bit once, apply that word under
the bit, and unload it. Its table is not read separately for every
native symbol. Even a conservative $`S_d\mathrm{poly}(d)`$ row
allowance is summable. The mantissas remain offline mask data.

Across every refinement substep the precision q halves strictly.
Hence

```math
\sum_{\rm refinements}q_{\rm new}\lt q_{d_0}=N/P_{d_0}^2.
\qquad\text{(26)}
```

Their controls and arithmetic contribute only $`\mathrm{poly}(n)`$.
The initial source preparation, final tag extension, attenuations,
and global replay phase banks each occur a constant number of times
and cost $`O(N)`$. The single outer amplification multiplies these
charges by a constant, not by the tree height.

Exact packings and fresh-cell permutations use the
[dirty-only permutation interface](DISJOINT_BLOCK_COMPILER.md#2-exact-routing-and-coarse-words)
on at most $`2N`$ modes including h. Each costs
$`O(n\sqrt N+\mathrm{poly}(n))`$ T gates and needs no clean
helper. There are only polynomially many such whole-mode operations,
not one per table row or precision digit. Their total is
$`O(\mathrm{poly}(n)\sqrt N)=O(N)`$. The central source guards,
metadata echoes, and small-register exact guards have additive
polynomial cost; they never multiply a long precision loader by n.
The fixed bottom part and signed normalization add $`O(N)`$.
This proves the count in (1).

## 8. Literal workspace and scope

The only external initialized wires are h and c. All growing private
bits and failure records are logical mode coordinates counted in the
complete volume (13). No separate clean history register is used:
sealing full rejected mode spaces prevents their later re-entry.

The dirty reservations are sequential. Increasing fixed $`d_0`$ and
the universal n threshold gives

```math
\begin{gathered}
\max_d q_{{\rm coef},d}\le N/P_{d_0}+O(n^6\log n)\lt N/2,\\
\max_d q_{{\rm ch},d}\le4N/P_{d_0}+O(\log n)\lt N/2,\qquad
\max q_{\rm ref}\le N/P_{d_0}^2\lt N/2.
\end{gathered}
```

Each such core supplies returned query selectors and a fixed number of
central-guard helpers. Permutations and coarse reflection routing use
$`O(\sqrt N+\mathrm{poly}(n))`$ dirty space instead. Fine calls
use one shared core of $`q+1=N+c_0+1`$ arbitrary dirty wires.
Their largest address includes h and at most n logical bits, so
self-borrowing applies once $`q+1\ge2(n+1)`$. Returned core helpers
handle the remaining central predicates; any fixed external helpers
fit inside the remaining $`n+6-c_0`$ wires. These are the same literal
reservations as the one-clean fine and disjoint-state interfaces.
At no time is a logical private bit also counted as dirty storage.
The peak is therefore $`N+n+7`$.

For the finitely many n below the universal threshold, the existing
literal-phase determinant-one finite-net fallback compiles every frame
at accuracy $`\epsilon/2`$ with the same two flags and dirty pool.
The finite maximum of T/N is absorbed into the universal constant.
Neither the fixed anchor size nor the asymptotic threshold is asserted
to be a practical crossover.

The real endpoint construction uses exact refinement, full-space sealed
arena allocation, compensated overlap, and the summable native ledger.
The next corollary extends its count to the phase-dressed complex
magnitude frame. Section 10 separately extends the count to arbitrary
precision. Sections 11--12 separately supply the one-clean extension
and a charged $`O(NL)`$ Clifford bound. No result here supplies an
optimal depth bound, efficient classical preprocessing algorithm, or
scalable native emitter.

## 9. Complex-magnitude endpoint

Let $`n\ge3`$, $`N=2^n`$, and $`\epsilon=2^{-N}`$. In addition
to the certified real rotation data, independently supply every unit
diagonal entry of $`D_\phi`$ with terminating certified real and
imaginary evaluation. This includes phases at zero-amplitude leaves;
they are part of the prescribed operator and are not inferred from its
state column.

**Corollary.** For the complete phase-dressed magnitude frame
$`W_{\mathbb C,\mathrm{mag}}=D_\phi W_{\mathbb R}`$, there is an
effective coherent Clifford+T circuit V satisfying

```math
\begin{gathered}
a=2,\qquad b=N+n+7,\qquad T(V)=O(N),\\
\|VJ_2-J_2(W_{\mathbb C,\mathrm{mag}}\otimes I_b)\|_{\rm op}
\le2^{-N}.
\end{gathered}
```

The constant is universal. Every logical column, literal common and
relative phase, dirty/reference input, and returned clean flag is
covered. The worst-case minimum count for each of the real and complex
magnitude families therefore obeys

```math
2N-1\le T^\star\le O(N),\qquad T^\star=\Theta(N).
```

**Proof.** Sections 6 and 8 give a real-frame circuit $`V_R`$ with
error at most $`\epsilon/2`$ at this same allocation, including the
finite-size fallback. For the
[packed literal-diagonal compiler](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames),
choose $`\ell=N+1`$. Since N is even,

```math
q_{\rm diag}=\left\lceil\frac{N+5}{2}\right\rceil=N/2+3,\qquad
b_{\rm diag}=q_{\rm diag}+n+1=N/2+n+4\le N+n+7.
```

Its full-isometry error is at most $`2^{-\ell}=\epsilon/2`$.
It preserves the literal diagonal phase and has the explicit count

```math
T(V_D)\le1344N+36q_{\rm diag}-18=1362N+90.
```

The remaining dirty wires are spectators; both circuits use the same
two clean flags and the same dirty pool sequentially.

Put $`J=J_2`$, $`R=W_{\mathbb R}\otimes I_b`$ and
$`D=D_\phi\otimes I_b`$. For the actual emitted product
$`V=V_DV_R`$,

```math
V_DV_RJ-JDR
=V_D(V_RJ-JR)+(V_DJ-JD)R.
```

Unitarity bounds the norm by $`\epsilon/2+\epsilon/2=\epsilon`$.
In particular, $`V_D`$ propagates the earlier clean leakage and dirty
disturbance unitarily; its approximation guarantee is used only on
$`JR`$. No intermediate reset or exact work return is assumed.
Tensoring any reference preserves this operator-norm estimate. The
actual inverse has the same initialized-isometry error, since

```math
V^\dagger J-J(DR)^\dagger
=V^\dagger(JDR-VJ)(DR)^\dagger.
```

The two T-counts add to $`O(N)`$. Setting $`D_\phi=I`$ includes
the real family and its
[root-angle witness](../../docs/FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness),
so the lower bound $`2N-1`$ transfers without a new premise. ∎

This corollary concerns the prescribed complex magnitude frame, not
arbitrary complex unitaries. Leaf-phase derivatives keep their separate
measurement stream. It changes neither the one-clean reservations nor
the real-frame count/depth guarantees of Result D.

## 10. Precision-uniform two-clean count

Let $`n\ge3`$, $`N=2^n`$, and let L be any integer with $`L\ge6`$.
The target data have the same certified-evaluation contract as above.
Put $`\epsilon=2^{-L}`$.

**Theorem.** Every prescribed complete real frame and every prescribed
phase-dressed complex magnitude frame has an effective coherent
Clifford+T circuit with two clean qubits and

```math
\begin{gathered}
T=O(N+L),\qquad
\|VJ_2-J_2(W\otimes I_b)\|_{\rm op}\le2^{-L},\\
b\ge L+n+7\quad\text{(real)},\qquad
b\ge L+n+8\quad\text{(complex magnitude)}.
\end{gathered}
```

The constants are universal and independent of both n and L. The
contract includes every logical column, literal phase, arbitrary dirty
and reference input, and returned clean work. All inverses are actual
emitted inverses. The fixed-precision corollary in Section 9 still gives
the tighter complex reservation $`N+n+7`$ at $`L=N`$.

### Separate dimension from precision

Choose one universal threshold $`n_0`$, large enough for the estimates
below and the fixed $`d_0`$ of Section 4. First suppose

```math
n>n_0,\qquad L>N/n,\qquad K=2^{\lceil\log_2L\rceil},
\qquad L\le K\lt 2L.
```

Retain the logical dimension N, the same anchor, and the same
$`P_d,M_d,p_d,B_d`$ from (2). Change only the precision scale:

```math
q_d=K/P_d^2,\qquad \beta_d=2^{-q_d},\qquad
\Delta_d=q_d+B_d,\qquad q_dM_d=4K.
```

Since $`L>N/n=2^n/n`$, a single $`n_0`$ makes K dominate every
fixed polynomial needed here, uniformly over this whole regime.
Thus each $`q_d`$ is a positive integer, $`q_d\ge B_d`$, and the
source ratios satisfy the rare-event bounds. At a jump of P, the same
four exact binary refinements halve q and double M each time.
Neither their identity nor their normalized attenuation depends on
$`L=N`$.

Apply the clipping convention with
$`H=L+\lceil\log_2n\rceil+6`$. Its complete-frame error is at most
$`\epsilon/16`$ and its retained certificates obey
$`E_v\le L+\lceil\log_2n\rceil+5`$. The optional signed
normalization has its separate $`\epsilon/16`$ budget. Declare a
node cheap when $`E_v\lt 4K/P_d`$; keep the same relative tolerances
and splitter losses. Use the same native-mask construction at the new
precision, fixing each emitted word before computing its mixed local gauge.
The cheap precision is now
$`q_{{\rm ch},d}=\lceil4K/P_d\rceil+O(\log P_d)`$.

For an expensive node write $`E_v=k_v\Delta_d+e_v`$.
The same inequalities give $`2P_d\le k_v\le2P_d^2`$ beyond the
universal threshold. Only this bounded multiplier is loaded. The
remainder can require $`O(\log(L+2))`$ classical bits; it stays in
the offline scalar masks, so it adds no logical mode register.
Formula (8) and its dyadic scalar realization are unchanged, with

```math
q_{{\rm coef},d}\le K/P_d+O(P_d^3\log P_d).
```

The bound $`f_{v,\nu}\le1/4`$, exact full-space guards, dirty
replay, and attenuation compensation depend on the unchanged
$`P_d,p_d,B_d`$ and these multiplier bounds. Hence their original
proofs apply on all occupied private inputs as well as initialized
ones. In particular, large L adds dirty precision width, not private
history or logical modes. The arena volume (12) and packing (13)
remain functions of N and d alone.

### Uniform error, native cost, and workspace

The finite-source and Cayley-tail bounds become

```math
\delta_d\le5\,2^{-4K},\qquad
e_n\le5n\,2^{-4K}\le\epsilon/192.
```

The last inequality holds uniformly for $`L>N/n`$ after increasing
$`n_0`$. The
[high-precision disjoint-block and state-source interfaces](DISJOINT_BLOCK_COMPILER.md#precision-parametrized-high-precision-form)
compile the fixed bottom blocks, initial source, and final tag
extensions at any fixed fraction of $`2^{-L}`$, with
$`O(N+L)`$ T gates and $`L+n+7`$ dirty qubits. Their coarse words
and mode dimensions still depend on N; their constant number of fine
calls uses precision $`q=L+c_0`$, with fixed $`c_0`$.
The same replacement applies to the attenuation and replay phase banks.
Their certified data and actual fine words, including inverses, split
the two source budgets exactly as in Section 6:
$`\epsilon_P,\epsilon_Q\le\epsilon/384`$.
Consequently the actual overlap has $`\delta\le\epsilon/96`$,
and the same normalization-two recovery has full-output error at most
$`\epsilon/8`$. Adding the bottom, clipping, and normalization errors
again gives $`\epsilon/2`$. No new recovery or reset is introduced.

All variable native charges remain explicit:

```math
\begin{aligned}
\sum_{d>d_0}S_dP_d&=O(N),\\
\sum_{d>d_0}q_{{\rm coef},d}
&\le K\sum_{d>d_0}P_d^{-1}+\mathrm{poly}(n)=O(L+N),\\
\sum_{d>d_0}q_{{\rm ch},d}&=O(L+\mathrm{poly}(n)),\\
\sum_{\rm refinements}q_{\rm new}&\lt K/P_{d_0}^2.
\end{aligned}
```

The exact routing and coarse state-source words still cost
$`O(\mathrm{poly}(n)\sqrt N)=O(N)`$. There are a constant number
of fine banks, each costing $`O(N+L)`$, and a single outer recovery.
Thus the total is $`O(N+L)`$.

For the literal workspace, enlarge the fixed $`d_0`$ and $`n_0`$
once so that all cheap, coefficient, and refinement dirty cores are
smaller than $`L/2`$. This follows from $`K\lt 2L`$, the displayed
precision bounds, and $`\mathrm{poly}(n)=o(N/n)`$.
Their returned selectors and guard helpers have only polynomial-in-n
size. Exact routing and coarse construction use
$`O(\sqrt N+\mathrm{poly}(n))=o(N/n)\lt L`$ dirty work in separate
stages. The fine core uses $`L+c_0+1`$ dirty wires; self-borrowing
holds because $`L+c_0+1\ge2(n+1)`$. The remaining fixed helpers
fit in $`n+6-c_0`$ wires, as in Section 8. Thus the peak is
$`L+n+7`$, independent of how far L exceeds N. Unused dirty wires
are spectators. No logical arena coordinate is also counted as dirty
storage.

### Splice the remaining regimes and compose phases

For $`L\le N/n`$, the existing
[two-clean real baseline](../../docs/OPERATOR_SOURCE_COMPILER.md)
uses $`b=L+n+7`$ and $`T=O(N+nL)=O(N)`$ at error $`2^{-L}`$.
For every $`3\le n\le n_0`$ and every L, use that same baseline:
$`N+nL\le\max\{1,n_0\}(N+L)`$. This is one universal
constant. In particular, the endpoint's fixed-accuracy finite-net
fallback is not extended to unbounded precision. These two regimes
and the growing-source regime prove the real theorem.

For the complex frame, invoke this uniform real theorem at precision
$`L+1`$ on $`b=L+n+8`$ dirty qubits. Compose the packed literal
diagonal at $`\ell=L+1`$, whose reservation is

```math
b_{\rm diag}=n+1+\left\lceil\frac{L+5}{2}\right\rceil
\le L+n+8.
```

Both counts are $`O(N+L)`$. Their two half-error full-isometry
bounds compose on the same two flags by the exact identity in
Section 9, including prior leakage, dirty/reference correlations,
literal phases, and the actual inverse. This proves the complex
claim at the displayed allocation. ∎

**Matching base-width corollary.** At the exact displayed reservations,
the total physical width is $`q_{\rm phys}=L+2n+9`$ for real frames
and $`q_{\rm phys}=L+2n+10`$ for complex magnitude frames.
When $`L\ge\max\{6,n\}`$, each is $`\Theta(L)`$. The existing
[worst-case lower bound](../../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
therefore gives

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{q_{\rm phys}}\right)
=\Omega(N+L),\qquad T^\star=\Theta(N+L).
```

The complex family contains the real one by taking $`D_\phi=I`$.
The quoted lower bound already combines the fixed-width counting
regime with its complementary square-root bound; no extra restriction
on $`q_{\rm phys}^2/(NL)`$ is needed here. This matching statement
concerns the exact base widths, not every larger dirty allocation.
Separately the width-independent root-angle witness gives
$`T^\star\ge2L-1`$, so the same order is optimal even at larger
widths when $`L\ge N`$.

For $`L\lt n`$ the uniform result is an upper bound. The existing
banked and sufficient-clean frontiers retain their stronger applicable
guarantees. No new Clifford-count, one-clean, depth,
classical-efficiency, or native-emitter bound follows here.

## 11. Precision-uniform one-clean count

Under the same certified-input contract, for every $`n\ge3`$ and
integer $`L\ge6`$, both frame families have an effective coherent
Clifford+T compiler satisfying

```math
\begin{gathered}
a=1,\qquad T=O(N+L),\qquad G=O(NL),\\
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L},\\
b\ge L+n+7\quad\text{(real)},\qquad
b\ge L+n+8\quad\text{(complex magnitude)}.
\end{gathered}
```

Constants are universal. The guarantee covers the complete operator on
initialized clean work, arbitrary dirty and reference inputs, literal
phases, and the actual inverse. It is the fine-phase flag c that is
removed; the anchor predicate h remains initialized. In particular,
this does not assert that the old two-clean complex endpoint in
Section 9 keeps its sharper $`N+n+7`$ width with one clean qubit.

### Lift before compiling the fine slots

Work first in the high regime of Section 10: $`n>n_0`$ and
$`L>N/n`$, after increasing a universal $`n_0`$ as needed. Keep
its precision scale K, mode allocation, clipping, and actual cheap
words. Before implementing the ideal fine diagonal slots, remove their
compiler flag c from the high-level words. It is not a logical mode
or a recovery coordinate of those ideal words. Retain h and all the
original logical suffix coordinates.

Reserve one arbitrary dirty signal r disjoint from every dirty core,
selector set, and bank throughout the source and its replay. Use the
[initialized-real-target lifting identity](CONJUGATE_SIGNAL_RECOVERY.md#a-real-initialized-target-is-sufficient)
and the [overlap recovery lemma](CONJUGATE_SIGNAL_RECOVERY.md#lift-an-overlap-before-its-recovery):

```math
\Pi_r(U)=|0\rangle\langle0|_r\otimes U
       +|1\rangle\langle1|_r\otimes\overline U.
```

Lift every exact native interlayer, including the coarse words, cheap
scalar words, actual inverse replays, predicate signs, and routing.
Real H and CNOT gates lift unchanged; S remains Clifford and a T
lift uses two T gates. The operation preserves products and adjoints.
Fix each actual cheap word and core before evaluating its mixed gauge;
the column replays that same word, now lifted, on that same core.
Thus the full rejected operator identity, including arbitrary dirty
input, holds in each r sector.

For every remaining ideal diagonal D, implement $`\Pi_r(D)`$ by
the actual finite phase word of
[the full-operator phase lemma](CONJUGATE_SIGNAL_RECOVERY.md#1-the-actual-finite-phase-word).
Its signal is r itself; there is no second borrowed rotation signal.
At precision $`q=L+c_0`$ the error is $`O(2^{-q})`$ on the
entire input space, including both signal sectors and their mixing.
Exact whole-word conditioning controls only the central source gates,
so the completed word is exactly identity on an inactive predicate
even for occupied r and arbitrary dirty work. All h guards survive.
Choose the fixed $`c_0`$ for the same constant number of fine slots
and source error budgets as in Section 10.

### Recover on the original anchor, excluding the signal

Let $`J_A`$ initialize only the original logical suffix on the active
h sector, with identity on the reduced input and all other dirty
work. Put $`J_r=I_r\otimes J_A`$. Before the fine replacements,
the two lifted source columns have compression

```math
J_r^\dagger\Pi_r(Q^\dagger P)J_r
=A\oplus\overline A,\qquad
\|A-W_{\rm upper}/2\|\le e_n.
```

The target $`W_{\rm upper}`$ and embedding $`J_A`$ are real.
Consequently both signal sectors approximate the same half-unitary,
regardless of whether the source words or their complementary blocks
are real. Certified-data and native fine-slot errors give

```math
\left\|J_r^\dagger\widehat Q^\dagger\widehat P J_r
       -I_r\otimes W_{\rm upper}/2\right\|
\le e_n+\epsilon_P+\epsilon_Q
\le\epsilon/96.
```

Here $`\epsilon_P,\epsilon_Q\le\epsilon/384`$ are complete
initialized-source-column bounds. Certified-data errors retain their
original amplitude-norm budgets: lifting a source column gives the
direct sum of that column and its conjugate, preserving its norm even
when the source target is complex. Only the native fine-slot
replacements require full-operator bounds. Their unitary telescope
retains all earlier leakage and any r mixing. Comparing the two actual
source columns gives the displayed overlap estimate without an
operator-norm approximation to the entire data-dependent source.
The normalization-two recovery
uses the actual $`U=\widehat Q^\dagger\widehat P`$ and its
actual inverse, with the guarded reflection

```math
R_h=I-2\Pi_{h=1,\,x\in\mathcal A}\otimes I_r\otimes I_{\rm dirty},
\qquad Z_hU_hR_hU_h^\dagger R_hU_h.
```

There is no test on r or on later private fields. This gives complete
output error at most $`\epsilon/8`$, including return of the signal
to its arbitrary input state. Every completed factor is exactly
identity on $`h=0`$. Reverse the actual anchor-predicate computation;
the ideal recovered output belongs to the original anchor, and the
same error bounds actual suffix and h leakage.

The real bottom blocks lie in $`\mathrm{SO}(2^{d_0})`$. Invoke the
[one-clean real disjoint-block interface](CONJUGATE_SIGNAL_RECOVERY.md#precision-parametrized-real-blocks-and-sources)
with $`\kappa=2^{d_0}`$, error $`\epsilon/4`$, and the same
clean wire after the upper circuit. Its complete isometry bound
composes with prior leakage unitarily. Initial real state extensions
and final real tag extensions use that interface's zero-clean source
form, or equivalently the same lifted ideal fine-slot words. The
clipping, normalization, upper recovery and bottom budgets still sum
to at most $`\epsilon/2`$.

### Literal width, counts, and remaining regimes

The fine reservation is $`q+1`$ core wires and one independent r.
Its $`L+O(1)`$ total, the fixed returned guard helpers, and at most
$`n+2`$ self-borrowed address selectors fit in $`L+n+7`$ after
increasing the fixed threshold on n. In particular
$`q+1\ge2(n+2)`$. Selectors are portions of the dirty core;
they are not extra initialized work. The smaller cheap, coefficient,
and refinement cores and coarse routing reservations from Section 10
fit with the same separate r in this pool in their respective stages.
No logical private coordinate is counted as dirty storage.

Lifting at most doubles the native T-count outside the constant number
of fine calls. Each such call is $`O(N+L)`$, so the total remains
$`O(N+L)`$. The Clifford ledger in Section 12 below also applies
with unbanked exact queries, giving $`G=O(NL)`$ at this base width.

For $`L\le N/n`$, use the existing
[one-clean real baseline](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger):
$`T=O(N+nL)=O(N)`$, $`G=O(NL)`$, and $`b=L+n+7`$.
For each $`3\le n\le n_0`$ and arbitrary L, that same baseline
has $`N+nL=O(N+L)`$ with one universal enlarged constant.
This covers all dimensions without a finite-net fallback at variable
precision.

For the complex magnitude frame, compose the real theorem at precision
$`L+1`$ with the existing
[one-clean literal diagonal](../../docs/ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames)
at precision $`L+1`$. Both use the same initialized wire and
$`b=L+n+8`$ dirty pool sequentially. Their half-error full-isometry
bounds compose by the identity in Section 9; counts add, preserving
$`T=O(N+L)`$ and $`G=O(NL)`$. If $`J=J_1`$ and the dirty
identity is included in W, the actual inverse satisfies

```math
V^\dagger J-JW^\dagger=V^\dagger(JW-VJ)W^\dagger.
```

Its error is the same. All these operator-norm identities remain valid
after tensoring an arbitrary reference. ∎

At the exact base dirty widths, the physical widths are now
$`q_{\rm phys}=L+2n+8`$ for real and $`L+2n+9`$ for complex
magnitude frames. For $`L\ge\max\{6,n\}`$ the existing lower
bound gives $`T^\star=\Theta(N+L)`$, as in Section 10. This
base-width matching statement is not extended to arbitrary larger b.

## 12. Optimal banked count frontier

For every $`n\ge3`$ and integer $`L\ge6`$, retain one clean
qubit and allow

```math
b\ge2(L+n+7)\quad\text{(real)},\qquad
b\ge2(L+n+8)\quad\text{(complex magnitude)}.
```

Under the same complete-input and certified-data contract there are
effective compilers with

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad G=O(NL),
\qquad \|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L}.
```

These match the worst-case lower bound throughout the stated widths:

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right).
```

### Bank exact native tables and sum their costs

In the high regime $`n>n_0`$, $`L>N/n`$, retain the construction
of Section 11. Use the
[banked real block and source interfaces](CONJUGATE_SIGNAL_RECOVERY.md#banked-real-interfaces)
for the fixed bottom blocks and constant source extensions. Their
literal base reservation includes the separate borrowed signal r.
For each fine, cheap or coefficient query, leave r outside both the
core and the additional word banks.

The [banked self-borrowing lemma](DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces)
implements the complete m-bit XOR query on an S-row table with

```math
T_{\rm query}=O\!\left(\sqrt{Sm}+m+\frac{Sm}{b}\right),
\qquad G_{\rm query}=O(Sm).
```

The base pool is $`B_0=L+n+7`$; the additional pool has at least
$`b/2`$ wires. Two halves of a long core provide one another's
returned address selectors, while at most $`b-B_0`$ bank wires
are used. For shorter metadata queries, selectors come from a
disjoint returned portion of the base core. Every completed query is
the original mask tensored with identity on the banks, on all inputs.
Hence the scalar source, its rejected operator, full-space guards,
and actual inverses are unchanged. Fix the banked row word before
forming its mixed gauge, and replay that exact word and core in the
column. Banks need not be initially zero or independent of the input.

At height d the coefficient tables have
$`Q_d=O(S_dP_d)`$ rows, $`S_d=N/2^d`$, and word widths

```math
m_d=O\!\left(L/P_d+P_d^3\log P_d\right),\qquad P_d\asymp(d+1)^2.
```

There are only constantly many whole coefficient masks per kernel;
small multiplier/type lookups have no larger row charge. The source
gates and central predicates add $`O(m_d+\mathrm{poly}(d))`$
per height, rather than a predicate charge per precision digit. The
convergent geometric row sums give

```math
\begin{aligned}
\sum_d Q_dm_d&=O(NL),\\
\sum_d\sqrt{Q_dm_d}&=O(\sqrt{NL}),\\
\sum_d m_d&=O(L+\mathrm{poly}(n)).
\end{aligned}
```

For example, the leading product is $`O(NL/2^d)`$, while each
remaining product is $`N2^{-d}`$ times a fixed polynomial in d.
Its square root is likewise summable. Cheap row words and their
actual column replays have $`O(S_d)`$ rows and
$`O(L/P_d+\log P_d)`$ word width, satisfying the same bounds.
The strictly halving refinement widths sum to $`O(L)`$ and use no
long table. The source, attenuation, tag, and replay fine banks occur
only constantly many times and each has $`O(N)`$ rows and
$`O(L)`$ word width. Their totals have the claimed three terms.

All exact routing, coarse blocks, short splitter words, and predicate
echoes have T-count at most $`N^{3/4}\mathrm{poly}(n)`$
and Clifford count at most $`N\mathrm{poly}(n)`$, by the
independent coarse ledger in the banked interface. For the growing
node-type words, load the preserved type bit once with an exact banked
one-bit query, execute its $`\mathrm{poly}(d)`$-length controlled
splitter word once, and unload. The T-count per height is
$`O(\sqrt{S_d}+1+S_d/b+\mathrm{poly}(d))`$; summing gives
$`O(\sqrt N+N/b+\mathrm{poly}(n))`$, within the claimed frontier.
Their Clifford cost is at most $`O(S_d\mathrm{poly}(d))`$ and
sums to $`O(N)`$. Since $`L>N/n`$,

```math
N^{3/4}\mathrm{poly}(n)=O(\sqrt{NL}),\qquad
N\mathrm{poly}(n)=O(NL),\qquad
\mathrm{poly}(n)=O(L).
```

Thus exact logical transport and the coarse compiler are charged, as
well as every precision table. The lift sends $`T\mapsto O(T)`$
and $`G\mapsto O(G+T)`$; $`T=O(NL)`$ throughout this regime.
One outer recovery multiplies the sums by a constant. The original
error and cleanup proof is unchanged because full mask operators are
unchanged. This proves the high-regime banked upper bound.

For the unbanked assertion of Section 11, use the same row and width
products with ordinary exact queries. Their T-count is $`O(S+m)`$
and Clifford count $`O(Sm)`$. The same sums give $`T=O(N+L)`$
and $`G=O(NL)`$ without extra banks. The coarse ledger is identical.

### Splice, compose, and match the lower bound

For $`L\le N/n`$, the existing one-clean grouped banked theorem
in [the one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md#6-complete-real-frames-and-the-grouped-extension)
gives

```math
T=O\!\left(\sqrt{NL}+L\ell_*(n)+\frac{NL}{b}\right),
\qquad G=O(NL).
```

Here $`\ell_*(n)`$ is its iterated grouping overhead. Its ratio
to $`\sqrt n`$ is bounded by a universal constant. Therefore
$`L\ell_*(n)\le\sqrt{NL}\,\ell_*(n)/\sqrt n=O(\sqrt{NL})`$
in this regime. For each remaining bounded $`3\le n\le n_0`$
and arbitrary L, use the existing banked baseline with count
$`O(\sqrt{NL}+nL+NL/b)`$. The factor n is uniformly bounded,
so it gives the desired formula. Both routes use one clean qubit and
the stated literal banked width, covering all regimes.

Compose the real result at $`L+1`$ with the banked one-clean literal
diagonal at $`L+1`$. Their required dirty threshold is
$`2(L+n+8)`$ and their counts add to the displayed complex bound.
The half-error hybrid and actual-inverse identity are those of
Section 11. Setting $`D_\phi=I`$ includes the real family.

Finally the existing
[complete-frame lower bound](../../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
at total width $`q_{\rm phys}=n+1+b`$ gives

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{q_{\rm phys}}\right).
```

At either stated banked threshold, $`q_{\rm phys}=\Theta(b)`$.
The lower and upper bounds match, with universal constants, for both
families. This is a count frontier at sufficient dirty width. It gives
no new T-depth bound and makes no claim about the full transition
between the base and doubled-base reservations. ∎

<a id="9-bounded-implementation-evidence"></a>
<a id="10-bounded-implementation-evidence"></a>
<a id="11-bounded-implementation-evidence"></a>

## 13. Bounded implementation evidence

The [exact helper](../../compiler_robust_hopf/growing_precision_source.py)
records the dyadic height schedule and normalized refinement factors,
with an optional requested precision L and a separate rounded source
scale K. Omitting L preserves the symbolic endpoint schedule.
The [focused tests](../../tests/test_growing_precision_source.py) check
the schedule identity $`M_dq_d=4N`$ and the four binary refinements at
a factor-four P jump. Exact checks at non-dyadic precisions, including
$`L>N`$, retain the same logical mode lengths and verify
$`M_dq_d=4K`$ and refinement compatibility. Small-L symbolic schedules
can have nonintegral q and do not certify a native construction; the
uniform proof uses the baseline in that regime. A separate rare-event
fixture uses an actual
scalar flag, nontrivial dirty rejection, controlled X, and inverse
rare preparation to check its full accepted operator.

The staged finite-mode fixture uses explicit fresh-cell swap
permutations, retained rejection spaces, frozen features, and dirty
reference entanglement. It checks the normalized $`\Gamma`$ factors,
their inverse compensation in the actual dirty replay, the shared
overlap, and recovery with actual circuit inverses and the original
anchor projector. Negative checks expose a missing refinement
normalization, reusing a rejection cell, and testing dirty or later
source fields in the recovery reflection. A separate perturbed-overlap
fixture checks the rejected column and the linear recovery error.

The five-wire phase-dressing check composes a diagnostic real recovery
with finite native packed-diagonal words, reusing both flags and the
dirty core. It retains intermediate leakage and tests literal phases,
the actual inverse and reference-entangled inputs. The addressed
direct sum is supplied; this is a composition check, not a native
growing-source or general table emitter.

The [conjugate recovery diagnostic](../../tests/test_conjugate_recovery.py)
adds one 16-dimensional initialized-real-target check. It lifts an
actual native interlayer and a supplied source, recovers on the
original mode anchor while excluding the arbitrary signal, and tests
dirty/reference entanglement, the actual inverse, retained leakage,
and the full-operator fine-slot hybrid. Negative controls detect a
signal-dependent recovery projector and an unlifted interlayer. The
uninitialized complement remains complex, so the test does not assume
whole-space reality. Its source matrix is supplied; the separate
finite-phase tests check the actual native phase word.

These fixtures use bounded dimensions and diagnostic dilations. They
do not emit the asymptotic compiler or simulate its exponentially fine
endpoint parameters. The theorem rests on the native identities,
complete-space allocation, error recurrence, and literal resource
ledger above.
