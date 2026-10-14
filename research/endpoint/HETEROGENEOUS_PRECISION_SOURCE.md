# A common carried source for nodewise precisions

[Precision gauges](PRECISION_GAUGE_COMPILER.md) · [Disjoint reflection sources](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources) · [Carried-source history](../../docs/FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization)

Complete connected Hopf trees admit a linear native compiler when every
nonzero branch exponent is at most a fixed multiple of $`N/\log_2n`$.
The exponent can change independently at every node. Cheap native
scalar merges and a common-source polynomial compiler share one correction
feature per node. One small-ratio geometric source is retained through the entire coarse row
network; a short integer shift multiplier and an exactly dyadic scalar
encode each node's precision. The source returns in one final
amplification, with the original two clean qubits and dirty reservation.

The earlier mixed-height family is a corollary. The construction does
not reschedule its unchanged high-precision merge banks, and does not
establish an unrestricted linear-count compiler.

## 1. Family, contract, and the upper sector

Fix a constant $`A_\star\gt0`$, an integer $`t\ge2`$, and

```math
n=2^{2t},\qquad N=2^n.
```

At every node of the complete height-n binary tree, allow either a known
zero rotation or a normalized rotation with smaller branch

```math
2^{-E_v-1}\le s_v\le2^{-E_v},\qquad
0\le E_v\le\frac{A_\star N}{\log_2n},\qquad E_v\in\mathbb Z.
```

The integer bins are supplied effective classical data, independently at
every nonzero node. No floor of an unknown logarithm or equality test
at a dyadic boundary is required. Exact-zero status may likewise be
supplied, or certified using the promised gap below the smallest allowed
nonzero branch. The supplied normalized chart has positive cosine and
nonnegative sine; the already priced signed normalizations may also be
attached if needed. The theorem in particular permits all rotations to
be nonzero, so there is no component-label invariance hypothesis.

**Theorem.** Every such prescribed complete real frame W has an effective
coherent Clifford+T compiler satisfying

```math
T=O_{A_\star}(N),\qquad
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-N},
\qquad b=N+n+7.
\qquad\text{(1)}
```

The constant depends only on $`A_\star`$, not on t, the exponent
assignment, or the branch values. The norm includes
every logical column, both initialized outputs, arbitrary dirty inputs,
and external references. All inverses are actual circuit inverses.
Certified classical evaluation and table construction are preprocessing;
no prepared precision state, quantum evaluator, reset, or measurement is
supplied.

Choose one sufficiently large universal positive integer C and put

```math
K=C\lceil\log_2(n+2)\rceil,\qquad N'=N/2^K.
\qquad\text{(2)}
```

All asymptotic arguments below are for a fixed threshold beyond which
$`K\lt n`$ and the displayed workspace inequalities hold.
The finitely many smaller n are handled at the end. Increase C once to
cover the fixed clean-work constants and the polynomial table factors
identified below; no particular numerical value of C is asserted.

Separate the bottom K binary depths from the upper frame, so
$`W=W_{\rm bottom}W_{\rm up}`$ and the upper part acts first. The bottom
part is a forest whose precision trace is at most

```math
O\bigl(K(E_{max}+\log n)\bigr)
=O_{A_\star}(N+\log^2n)=O_{A_\star}(N).
\qquad\text{(3)}
```

The precision-trace theorem therefore compiles it in linear T count.
The upper frame acts only when the last K logical bits are zero. On that
sector it is a complete height-$`n-K`$ frame on $`N'`$ labels, with the
same known zero rotations where present. Standard aligned labels are
retained throughout.

Compute $`h=[\text{last K bits}=0]`$ into the first external clean qubit
using the exact borrowed predicate construction. Keep h unchanged until
the entire upper compiler, including its fine-gauge recovery, is complete.
On $`h=1`$ the suffix supplies K initialized private wires. The second
external clean qubit is available as a feature-mode flag. Every completed
subroutine is conditioned to be exactly identity on $`h=0`$, on the
whole space of the remaining registers, using the
[conditional-suffix rules](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally).
Finally reverse the actual predicate word. A complete-isometry hybrid
charges any suffix leakage; no intermediate reset is used.

## 2. Mixed local merges with one feature per node

Write a target normalized node rotation as

```math
R_y(\theta)=\begin{bmatrix}c&-s\\s&c\end{bmatrix},
\qquad \tau=\tan(\theta/2)=\frac{s}{1+c}.
```

Call a nonzero node cheap when $`E_v\lt2N/n`$ and expensive otherwise.
This comparison is on the supplied integer bin. Both constructions below
use $`\eta=1/(1000n^2)`$ and retain literal complex phases.

### Expensive nodes: a Cayley rotation and a native splitter

For an expensive node v choose a positive dyadic $`u_v`$, with m fractional
bits, approximating $`2^{E_v}\tau_v`$, where

```math
m=O(\log n),\qquad 2^{-m}\le\frac1{16000n^2},
\qquad \frac18\le u_v\le1.
```

These bounds hold uniformly for sufficiently large n, by certified
rounding of $`2^{E_v}\tau_v\in[1/4,1]`$. Put $`t_v=2^{-E_v}u_v`$.
For the polynomial compiler only, put $`t_v=u_v=0`$ at every cheap or
known zero node; its actual cheap merge will be a separate native interlayer.
The exact mathematical coarse gate is the dyadic Cayley rotation

```math
J=\begin{bmatrix}0&-1\\1&0\end{bmatrix},\qquad
R(t)=(I+tJ)(I-tJ)^{-1}
=\frac1{1+t^2}\begin{bmatrix}1-t^2&-2t\\2t&1-t^2\end{bmatrix}.
\qquad\text{(4)}
```

The derivatives of its first column, together with
$`s_v\ge2^{-E_v-1}`$, give relative error at most
$`\eta`$ in both nonzero entries, after the fixed
rounding margin above. All coefficients in (4) are exact computable
rationals; (4) is a mathematical target, not an assumed exact native gate.
Write its first column as $`(a_v,b_v)`$, with real positive entries.
First apply $`R(t_v)^\dagger`$ to the child carriers, then apply a fixed
exact native splitter $`U^{\rm sp}`$ between the parent carrier and its
fresh feature. Choose its literal first column $`(r,t)`$ so that

```math
l=|t|^2=1-|r|^2\in[1/(1600n),1/(800n)].
```

A strict interior target and a polynomial-in-log(n) native word suffice.
Its entries may be complex; neither r nor t is replaced by its magnitude.
The resulting local column on arbitrary dirty inputs is

```math
\begin{bmatrix}B_v^\dagger\\F_vL_v^\dagger\end{bmatrix},\qquad
B_v=R(t_v)\,\mathrm{diag}(\bar r,1),\quad
L_v=R(t_v)e_0=(a_v,b_v)^T,\quad F_v=tI.
```

Thus $`B_vB_v^\dagger+lL_vL_v^\dagger=I`$ and
$`F_v^\dagger F_v=lI`$. Only one feature mode is used.

### Cheap nodes: the actual lossy native scalar

Use the [lossy scalar construction](PRECISION_GAUGE_COMPILER.md#3-lossy-native-scalars-supply-their-own-correction-features).
Choose the dyadic loss parameter

```math
\tau_v=2^{-E_v-\lceil\log_2(200/\eta)\rceil},\qquad
\eta s_v/400\le\tau_v\le\eta s_v/100.
```

Fix one common precision $`q_{\rm ch}=\lceil2N/n\rceil+O(\log n)`$
and one physical dirty source core for all cheap nodes, with a fixed
sufficient margin. Let $`z_v`$ be the actual accepted scalar of the
word emitted at that precision,
with $`|z_v-(1-\tau_v)e^{2i\theta_v}|\le\tau_v/100`$. The
exact emitted word $`U_v`$ and its actual rejection satisfy

```math
U_v|0\rangle=|0\rangle\bar z_vI+|1\rangle F_v,\qquad
F_v^\dagger F_v=l_vI,\qquad l_v=1-|z_v|^2.
```

Take

```math
a_v=(1+z_v)/2,\qquad b_v=(z_v-1)/(2i),\qquad
B_v=\begin{bmatrix}a_v&-b_v\\b_v&a_v\end{bmatrix},\qquad
L_v=\ell=(1,-i)^T/\sqrt2.
```

The fixed Clifford eigenbasis construction implements exactly
$`[B_v^\dagger;F_vL_v^\dagger]`$, with
$`B_vB_v^\dagger+l_vL_vL_v^\dagger=I`$. Its rejected operator
is retained, not replaced by a scalar completion. Choose and fix each
actual row word before computing the gauge. These words are exact native
interlayers even though their accepted carrier blocks are contractions;
the precision-trace theorem is not being used as an exact coarse unitary.

At a known normalized zero use $`B_v=I`$, $`F_v=0`$, and no feature
query or splitter. One may fix $`L_v=e_0`$ there; its coefficient below
is zero. Every node retains its allocated feature label, even when unused.

### Full-column gauge for the mixed network

For either nonzero chart, let $`B_v=[b_{0,v},b_{1,v}]`$ and set
$`d_{root}=1`$ and

```math
d_L=d_vc_v/(b_{0,v})_0,\qquad d_R=d_vs_v/(b_{0,v})_1.
```

Define the two scalar correction coefficients by the local identity

```math
\mathrm{diag}(d_L,d_R)^{-1}(-s_v,c_v)^T
=b_{1,v}C_v+L_vA_v.
```

For a cheap node these are exactly the coefficients in the lossy gauge:

```math
C_v=\frac{z_v-\cos(2\theta_v)}{2id_vc_vs_v},\qquad
A_v=-\frac{(z_v-e^{2i\theta_v})(z_v-e^{-2i\theta_v})}
                 {2\sqrt2\,d_vc_vs_v}.
```

For an expensive node, direct expansion in the orthonormal vectors
$`(-b_v,a_v)^T`$ and $`(a_v,b_v)^T`$ gives

```math
d_L=\frac{d_vc_v}{\bar r a_v},\qquad
d_R=\frac{d_vs_v}{\bar r b_v},\qquad
C_v=\frac{\bar r a_vb_v}{d_vc_vs_v},\qquad
A_v=\frac{\bar r(c_v^2b_v^2-s_v^2a_v^2)}{d_vc_vs_v}.
```

At a zero node use equal child gauges $`d_L=d_R=d_v`$,
$`C_v=1/d_v`$, and $`A_v=0`$, without evaluating a zero denominator.

Let $`W_B`$ be the complete coarse frame obtained by inserting the
matrices $`B_v`$ in the tree recursion. Let $`L_B`$ contain the vector
$`L_v`$ injected at node v and propagated through its coarse descendants.
Let D contain the leaf gauges, let $`C_{root}=1`$, and place $`C_v`$
and $`A_v`$ at v's original detail column. Descendant induction gives

```math
\begin{gathered}
W_{\rm up}=D(W_BC+L_BA),\\
Q_0|x,0\rangle=
\begin{bmatrix}W_B^\dagger x\\
 (F_v(L_B^\dagger x)_v)_v\end{bmatrix},\qquad Q_0^\dagger Q_0=I.
\end{gathered}
\qquad\text{(5)}
```

Indeed each coarse root is mapped by its descendant diagonal gauge to
$`d_v`$ times the target root. The displayed two-coordinate detail
identity then supplies every detail column. The row identity follows by
composing the actual local columns: accepted carrier blocks are scalar
on the dirty register, while a feature leaves the carrier network at its
unique node. Later stages fix that feature. Thus its one surviving
rejection is exactly $`F_v`$. No common attenuation of sibling carriers
is needed; their differing losses are already contained in $`B_v`$.

The physical network is a full unitary on all $`2N'`$ modes and all
dirty inputs, including occupied features. This full dirty-dependent
unitary $`Q_0`$ is the target of the carried-source compiler.

The constants have uniform slack. Every path has at most n nodes, and
the relative coarse errors are at most $`\eta`$ after removing the
cheap node's local unit-modulus phase $`e^{i\theta_v}`$. That phase
remains in the complex gauge; the following bound concerns its modulus.
For $`n\ge16`$,

```math
(1+\eta)^{-n}\le|d_v|
\le(1-\eta)^{-n}(1-1/(800n))^{-n/2},\qquad
.999\lt|d_v|\lt1.002.
```

Consequently $`.98\lt|C_v|\lt1.02`$ at both node types. The cheap
bound $`|A_v|/l_v\lt.8`$ is unchanged. At an expensive node, writing
$`a_v=c_v(1+e)`$ and $`b_v=s_v(1+f)`$ gives

```math
A_v=\frac{\bar r}{d_v}c_vs_v\bigl[(1+f)^2-(1+e)^2\bigr],\qquad
|A_v|\le3\eta,\qquad |A_v|/l\le4.8/n\le.3.
```

This cancellation is why a very small branch does not enlarge its
feature coefficient. The root coefficient and all zero-node coefficients
obey the same uniform bounds.

For clarity, if p and m are the left and right carrier labels of a node,
the merge acts on $`(0,p),(0,m)`$ and its splitter on
$`(0,p),(1,m)`$, where the first bit is the feature flag. In standard
aligned labels, p and m differ
only in that node's target bit. Conjugating by CNOT from the feature flag
to that target bit maps the splitter pair to a feature-bit rotation
guarded by target zero and the unchanged logical suffix. This supplies
an explicit native matching layout on all modes, with the conjugating
CNOTs cancelling on other pairs. No whole logical endpoint or exponent
is copied to a fresh clean register.

Let $`\mathcal R_d`$ apply only the expensive Cayley merge adjoints at
depth d, and let $`C_d`$ apply the cheap native three-port merges and the
expensive fixed splitters on their disjoint node modes. It is an exact
native unitary on the logical, feature, and dirty registers. On expensive
carrier pairs the polynomial generator is $`J_d=-J`$; every other row
has identity polynomial merge. There is one list of
$`R=n-K\le n`$ depths. Thus

```math
Q_0=C_R\mathcal R_R\cdots C_1\mathcal R_1.
\qquad\text{(6)}
```

The same construction can include a preserved row-attenuation flag: on
its failure value set all polynomial coefficients to zero and disable
every native merge and splitter. Its target is exactly I, while the normalization
below stays uniform.

## 3. One common finite source with a small ratio

Because n is dyadic, all of the following choices are exact:

```math
\Delta=\frac{N}{n^2},\qquad P=n,\qquad M=4n^2,\qquad
p=\frac1{8n},\qquad B=6+2\log_2n,\qquad \beta=2^{B-\Delta}.
\qquad\text{(7)}
```

The spacing $`\Delta`$ is an integer. For sufficiently large n,
$`0\lt\beta\le1/2`$. On an M-mode source register define

```math
\begin{gathered}
g_\beta=a_\beta\sum_{j=0}^{M-1}\beta^j|j\rangle,
\qquad
a_\beta=\sqrt{\frac{1-\beta^2}{1-\beta^{2M}}},\\
S_d=\sum_{j=d}^{M-1}|j-d\rangle\langle j|\quad(0\le d\lt M),
\qquad S_d=0\quad(d\ge M).
\end{gathered}
\qquad\text{(8)}
```

For $`0\le d\lt M`$, only the final d coordinates contribute to the
residual. Direct summation gives

```math
\|S_dg_\beta-\beta^dg_\beta\|
=\beta^M\sqrt{\frac{1-\beta^{2d}}{1-\beta^{2M}}}
\le2\beta^M.
\qquad\text{(9)}
```

For $`d\ge M`$ the residual is $`\beta^d\le\beta^M`$, so the final
bound in (9) holds for every nonnegative integer d. This is a full-output
relation on the source column, not just a projected overlap. Define
$`E_g=|g_\beta\rangle\otimes I`$ with exactly the same source on every
carrier, feature, and row-attenuation-control value, and tensor it with
identity on the entire arbitrary dirty register. In this chapter the
system for the inner compiler includes those dirty inputs. Thus every
actual cheap merge, including its dirty action, intertwines with
$`E_g`$. No unchanged class label or intermediate dirty return is required.

### Charged preparation and its full-space extension

Choose the real source extension

```math
v=\frac{e_0+g_\beta}{\|e_0+g_\beta\|},
\qquad A=R_vR_{e_0},\qquad Ae_0=g_\beta.
\qquad\text{(10)}
```

It acts only on the M source modes, with identical copies on all
reduced system labels. Including
feature and row-attenuation bits, the complete mode dimension is at most
$`4N'M\le N`$ once C and the asymptotic threshold are fixed. Including
the external conditioning bit still gives at most $`2N`$ addressed modes. The
[one-clean disjoint-reflection interface](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)
therefore implements A with four fine diagonal calls and one exact center
sign. Denote its emitted word by $`\widehat A`$. It has cost $`O(N)`$
and full initialized-isometry error $`\delta_A`$
of any fixed chosen multiple of $`2^{-N}`$.

There is no need for a coarse preparation. Indeed
$`\|v-e_0\|=O(\beta)\le1/(16M)`$ for sufficiently large n,
so the reflection interface's local state chart takes its coarse word
to be identity. Its Walsh interlayers act only on source bits, and its
diagonals and center signs preserve every reduced logical, feature, and
row-attenuation address label. The emitted word preserves those labels
exactly, even when its private compiler flag leaks. Disturbance of its
arbitrary dirty compiler core is included in $`\delta_A`$; no exact
dirty-return claim is imposed on this approximate preparation.
Its inverse below is its actual native inverse. The flag is distinct
from all private flags used by the carried stream.

## 4. Exact dyadic local kernels

The Cayley series in a carrier pair is

```math
R(-t_v)=I+2\sum_{\nu=1}^{\infty}t_v^\nu J_d^\nu.
\qquad\text{(11)}
```

For every expensive node compute classically

```math
E_v=k_v\Delta+e_v,\qquad
k_v=\lfloor E_v/\Delta\rfloor,\qquad 0\le e_v\lt\Delta,
\qquad 2n\le k_v\le k_{max}\le\frac{A_\star n^2}{\log_2n}.
```

This division is on a supplied integer, not an unknown real exponent.
At a cheap or known zero node put $`k_v=e_v=u_v=0`$. These multiplier
symbols are distinct from the local coarse matrix entries. For
$`1\le\nu\le P`$, define

```math
f_{v,\nu}
=\frac{2P(1-p)}p\,u_v^\nu2^{-\nu(k_vB+e_v)}
=2n(8n-1)u_v^\nu2^{-\nu(k_vB+e_v)}.
\qquad\text{(12)}
```

Set f to zero on every frozen or feature mode and on the disabled
row-attenuation sector. On a carrier pair the same f is assigned to both
labels. These numbers are exact dyadics, obey $`0\le f_{v,\nu}\le1/4`$,
and need at most

```math
q_{\rm coef}\le P(m+k_{max}B+\Delta)+O(1)
=\frac Nn+O_{A_\star}(n^3).
\qquad\text{(13)}
```

binary fractional digits. The
[raw operator-source scalar](../../docs/OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients)
therefore realizes an exact accepted block $`f_{v,\nu}I`$ with one
retained scalar flag and $`q_{\rm coef}+O(1)`$ arbitrary dirty core
wires. Enlarge the finite sign grid by a fixed number of places when
necessary. This is the unamplified scalar word; no inverse polynomial
or unrounded real coefficient is being encoded.

The long remainder $`e_v`$ appears only in the offline scalar-mask
table. It is never loaded into a quantum register: unlike $`k_v`$, its
ordinary binary description can have order n bits. Only the multiplier
$`k_v`$, of $`O(\log n)`$ bits for fixed $`A_\star`$, is loaded
coherently.

Prepare a uniform P-valued label $`\ell`$ and interpret
$`\nu=\ell+1`$. Prepare $`\log_2(8n)`$ Hadamard bits and designate
their all-zero string as a rare event, of probability p. In the common
event SELECT is identity. In the rare event apply the scalar word for
f, the shift dilation for $`S_{k_v\nu}`$, and the literal quarter-turn
$`J_d^\nu`$ on the selected carrier pair. Its phase is retained; modulo
four powers suffice. The quarter-turn exponent is still $`\nu`$, not
$`k_v\nu`$.

Load $`k_v`$ by an exact whole-word table, multiply by $`\nu`$ into
$`O(\log n)`$ clean arithmetic work, and let $`d=k_v\nu`$ be the
full nonnegative product. Toggle a distinct underflow flag by
$`[j\lt d]`$ **before** reducing d modulo M, then subtract d modulo M
from the source index. If $`d\ge M`$, the predicate is true for every
source index, so the flag-zero block is exactly zero. Otherwise it is
exactly $`S_d`$. This gives the required contraction on every source
input, including a selected shift beyond the source cap. Reverse the
product arithmetic and multiplier lookup before the quarter-turn changes
the logical endpoint. Their controls are unchanged, so their scratch
returns exactly on every scalar and underflow failure branch.

The source-underflow and scalar-rejection flags are distinct. Both are
retained. The rare bits and power label are unprepared by the actual
Hadamard inverses. All coefficient queries return their borrowed
selectors exactly. Their address values are unchanged during the query;
the subsequent quarter-turn stays in the same carrier pair, on which f
is constant. Every SELECT case specifies a full unitary, including
rejected and inactive inputs.

Let F consist of these rare bits, the power label, both retained flags,
and any retained private work. Exact arithmetic scratch is restored at
the boundary. The actual kernel $`K_d`$ has F-zero block

```math
T_d=(1-p)I+\frac pP\sum_{\nu=1}^P\sum_{v\text{ at depth }d}
S_{k_v\nu}\otimes f_{v,\nu}J_{d,v}^\nu.
\qquad\text{(14)}
```

Here $`J_{d,v}`$ is $`-J`$ on v's carrier pair, extended by zero
outside that pair. The pairs at a depth are disjoint. The multiplier
and f tables take the same values on both labels of each pair; feature,
frozen, and disabled rows have zero correction. Formula (14) follows by
projecting the
two separate dilation flags, not by dropping a rejected term of a
later operation on the same flag. Since it is a compression of a
unitary, $`\|T_d\|\le1`$.

Use the exact identity
$`t_v=\beta^{k_v}u_v2^{-(k_vB+e_v)}`$ in (11). With
$`D_d=(1-p)\mathcal R_d`$, equations (9), (12), and (14) give

```math
\|T_dE_g-E_gD_d\|
\le 2p\beta^M
 +\frac{2(1-p)t_{max}^{P+1}}{1-t_{max}}
=:\delta_{\rm ker}.
\qquad\text{(15)}
```

The bound holds on every system and dirty input, with references.
Disjoint carrier pairs make the nodewise contribution a maximum,
so no factor equal to the number of nodes occurs in (15).
Frozen and disabled sectors have exactly $`(1-p)I`$ as their kernel;
they do not acquire a different attenuation. Also
$`\|D_d\|=1-p`$. Every expensive node has $`E_v\ge2N/n`$;
all other polynomial rows have $`t_v=0`$. Increasing the fixed
threshold so that $`4n^2B\le N`$ gives

```math
M(\Delta-B)\ge3N,\qquad
t_{max}^{P+1}\le2^{-2N},\qquad
R\delta_{\rm ker}=o(2^{-N}).
\qquad\text{(16)}
```

## 5. Whole-network history, normalization, and return

Share F among all R kernels. Initialize a history counter H of
$`\lceil\log_2R\rceil`$ bits, omitting it if $`R=1`$. After each
kernel except the last, increment H if and only if $`F\ne0`$. There
are at most $`R-1`$ increments. Hence a path beginning at H=0 cannot
wrap back to zero. Later kernels never act on H.

Expand in the F and H computational bases. Final H=F=0 retains exactly
those paths whose F value was zero at every intervening boundary. It
therefore inserts the intermediate projectors, giving the exact accepted
operator on carried source and system

```math
V_{\rm acc}=(I\otimes C_R)T_R\cdots(I\otimes C_1)T_1.
\qquad\text{(17)}
```

Each $`C_d`$ is defined on every occupied feature and arbitrary dirty
input. It acts on those system registers and never on the carried source
or private history. The source is the same on every system input. Thus
$`(I\otimes C_d)E_g=E_gC_d`$, and contractive telescoping of (15)
proves

```math
\|V_{\rm acc}E_g-E_g\gamma Q_0\|
\le R\delta_{\rm ker},\qquad
\gamma=(1-p)^R\ge1-Rp\ge7/8.
\qquad\text{(18)}
```

This estimate includes the carried source's residual. It makes no claim
that an intermediate source or F register has returned on the actual
state. Failed branches continue through defined full unitaries.

Apply $`\widehat A`$ once, run the complete stream, and apply the actual
native inverse $`\widehat A^\dagger`$.
Comparison with the ideal source extension (10), using a unitary hybrid,
adds at most $`2\delta_A`$ to the final zero-source block error. The
source-preparation compiler flag is included among final tested work.

Use one additional retained scalar flag, outside F and its history test,
to encode

```math
w=\frac1{2\gamma}\in[1/2,4/7]
```

by a raw operator-source scalar word. With $`q+1`$ dirty generators its
signed grid has spacing $`2^{1-q}`$. Since w is an exact rational,
rounding to that grid gives $`\widetilde w`$ with
$`\delta_w\le2^{-q}=O(2^{-N})`$ for $`q=N+c_0`$. This one fine
call costs $`O(N)`$. Its accepted operator is
$`\widetilde w I_{\rm dirty}`$, so it multiplies the completed stream
block even when its dirty core is shared. Its failure flag never
re-enters any other scalar subroutine. If U is the whole actual word,
and J initializes all its private work, its accepted block satisfies

```math
\|J^\dagger UJ-Q_0/2\|
\le R\delta_{\rm ker}+2\delta_A+\delta_w.
\qquad\text{(19)}
```

The [normalization-two amplification lemma](../../docs/FAULT_TOLERANT_COMPILER.md#4-from-a-half-unitary-block-to-the-complete-output)
applied to $`-URU^\dagger RU`$, with $`R=I-2JJ^\dagger`$, now gives
the full initialized-isometry approximation to $`Q_0`$. It includes
every private output, all dirty disturbance, and arbitrary occupied
feature inputs. There are two forward words, one actual inverse, and
two work reflections, so only a constant number of fine calls.

For the external upper-sector marker use the literal conditional form

```math
Z_h U_h R_h U_h^\dagger R_h U_h,
\qquad R_h=I-2\Pi_{h=1,\,\mathrm{private}=0}.
\qquad\text{(20)}
```

On $`h=1`$ this is the required leading-minus word. On $`h=0`$ every
factor, including its complete output, is exactly identity. The
reflection does not test arbitrary logical, feature, or dirty inputs.
When compiling the row-controlled version, its preserved attenuation
control is likewise system data, not initialized work.

## 6. Replaying the mixed features and recovering the frame

Take $`\Lambda=5/4`$ and $`\kappa=8/5`$, so
$`\Lambda\kappa=2`$. The ideal row source first attenuates by
$`D^\dagger/\Lambda`$, routes its failure to one bank, and applies
the full mixed $`Q_0`$ on the shared sector. This is a contraction
because $`|D|\lt1.002`$. For the column source put, at every nonzero
detail label,

```math
x_v=C_v/\kappa,\qquad y_v=A_v/(\kappa l_v),\qquad
\rho_v=(|x_v|^2+l_v|y_v|^2)^{1/2}.
```

At the root and zero nodes use $`x=C/\kappa`$, $`y=0`$, and omit
all divisions by l and feature queries. Section 2 gives
$`.6\lt\rho_v\lt.7`$ and $`|y_v|/\rho_v\lt1`$. For a cheap
node the accepted scalar of $`U_v`$ is $`z_{0,v}=\bar z_v`$; for
an expensive node use the very same fixed splitter, with
$`z_{0,v}=r`$ and $`F_v=tI`$. In both cases
$`2|z_{0,v}|\gt1`$.

The [two-query replay identity](PRECISION_GAUGE_COMPILER.md#4-matching-arbitrary-dirty-rejection-and-recovering-the-frame)
is a full-operator consequence of unitarity:

```math
U_v\begin{bmatrix}e^{i\phi_v}I&0\\0&I\end{bmatrix}U_v^\dagger|0\rangle
=\begin{bmatrix}
1+(e^{i\phi_v}-1)|z_{0,v}|^2\\
(e^{i\phi_v}-1)\bar z_{0,v}F_v
\end{bmatrix}.
```

Choose $`2|z_{0,v}|\,|\sin(\phi_v/2)|=|y_v|/\rho_v`$ and
apply separate literal output-sector phases. This prepares the normalized
pair $`(x_v/\rho_v,(y_v/\rho_v)F_v)`$ while retaining its actual
dirty rejection. After an initial attenuation by $`\rho_v`$, its
shared column is $`(x_v,y_vF_v)`$. Failure stays in the opposite
bank from the row source.

The chronology and fine-call count matter. Construct
$`U_* = \bigoplus_v U_v`$ on the preserved original detail label and
feature flag. Each cheap block replays the exact row scalar word on the
exact same physical dirty wires. The midpoint-to-parent relabeling of
the row bank, or its identical finite mask table indexed by the detail
label, supplies that word. It is not replaced by a fresh synthesis
having merely the same accepted scalar. The common cheap source width
and physical source circuit were fixed before defining $`z_v,F_v`$.
Expensive blocks likewise use the exact row splitter word. All root and
zero blocks are identity.

Execute all blocks of $`U_*^\dagger`$, then **one global** fine
diagonal $`\Phi(j,b)`$, then all blocks of $`U_*`$, followed by the
two output-sector phase banks. The algebraic middle word is
$`U_*\Phi U_*^\dagger`$. Its chronological order is not reversed.
After the initial attenuation, this entire replay and every phase bank
are enabled only on its shared flag value. They are exactly identity on
the attenuation-failure sector, on every feature and dirty input.
Consequently column failure remains in its original $`a=1,b=0`$ bank.
The detail labels are disjoint, so grouping the actual queries by their
original depth changes no block. In particular there is no fine
precision-N phase call inside each depth loop. The phase and attenuation
banks use a separate initialized suffix compiler flag; no occupied
attenuation flag is presumed clean. Certified approximate phase data
can be chosen in amplitude norm also when a requested component vanishes.

The ideal shared source columns are exactly

```math
\begin{gathered}
P_{\rm shared}=\frac1\kappa
\begin{bmatrix}C\\(A_vF_v/l_v)_v\end{bmatrix},\qquad
Q_{\rm shared}=\frac1\Lambda
\begin{bmatrix}W_B^\dagger D^\dagger\\
(F_v(L_B^\dagger D^\dagger)_v)_v\end{bmatrix},\\
Q^\dagger P=\frac12D(W_BC+L_BA)=W_{\rm up}/2.
\end{gathered}
\qquad\text{(21)}
```

The zero-feature entries in this display are defined as zero without
dividing by l. Its overlap uses $`F_v^\dagger F_v=l_vI`$ on all
dirty inputs. Distinct features prevent cross-node dirty products. The
shared bank has $`2N'-1`$ modes and the two failure banks have $`N'`$
modes each, totaling $`4N'-1`$ within the feature and attenuation flags.

The feature flag is the second external clean qubit. The row/column
attenuation flag and separate fine-compiler flag are active-suffix wires.
The first external flag remains h throughout. Source, history,
normalization, and preparation flags are distinct private suffix wires.
Replace the ideal mixed row-network call by Section 5's approximation
on all feature and dirty inputs, controlled on the shared attenuation
value. Its other sector has target I. The inner work reflection excludes
the arbitrary dirty input; it tests only its initialized private work.

Compile the constant number of remaining fine banks at
$`q=N+c_0`$ for a fixed sufficiently large $`c_0`$. A complete-isometry
hybrid bounds these replacements, propagating earlier leakage unitarily
and using actual inverse words. Apply the outer normalization-two
recovery with the conditional sign convention (20). Its reflection tests
every initialized private and suffix wire and the feature flag; h is
the preserved enable, excluded from the zero-work test. It does not test
reduced logical or dirty input. Fixed error fractions give upper-frame
error below $`2^{-N}/2`$.

Finally reverse the upper-sector predicate and compose with the bottom
compiler from (3), also at error below $`2^{-N}/2`$. This proves (1).
The cheap native words themselves, including their dirty action, are
executed as exact interlayers and exact replays. Only the expensive
Cayley gates require the carried polynomial compiler. No approximately
returned low-precision frame is used as an exact coarse layer.

## 7. Native count and literal workspace

At one stage the complete mode-and-power coefficient table has at most
$`O(N'P)`$ rows, including fixed feature and attenuation controls.
Whole-word self-borrowed queries implement its dyadic sign masks with
T cost $`O(N'P+q_{\rm coef})`$. The output mask width affects the
Clifford work, not the row contribution to T. The source-word constant
is fixed. The short multiplier table costs $`O(N')`$ per stage, or
$`O(N'P)`$ if conservatively addressed by the power as well. Its
$`O(\log n)`$-bit output does not supply enough self-borrowed selectors;
instead this completed query borrows and exactly returns $`O(n)`$
selectors inside the much larger scalar-source dirty core. Exact shift
arithmetic, rare predicates, quarter-turn powers, and history operations
have a fixed polynomial cost in n. Crucially, (13) gives

```math
nq_{\rm coef}\le n^2\Delta+O_{A_\star}(n^4)
=N+O_{A_\star}(n^4).
```

The long residual-exponent charge is therefore retained. The
conservative full-table ledger over $`R\le n`$ stages is

```math
T_{\rm kernels}
=O\left(n^2\frac{N}{2^K}+nq_{\rm coef}+\mathrm{poly}(n)\right)
=O_{A_\star}(N).
\qquad\text{(22)}
```

The exactly native splitters have polynomial-length coarse descriptions.
Their controlled matching tables and suffix predicates cost
$`N'\mathrm{poly}(n)+\mathrm{poly}(n)`$, absorbed by a
single sufficiently large choice of C. The cheap/expensive/zero selector
is an arbitrary node-type table. At a fixed depth it depends on the
unchanged node prefix, excluding the pair target. Compute its constant
number of bits, execute the selected carrier/feature word, and reverse
that lookup. Each completed lookup returns its borrowed dirty selectors
before the selected word runs. Its inverse lookup then works on the
possibly changed arbitrary dirty state and returns its selectors exactly.
A conservative total $`O(nN')`$ covers these type tables.

The cheap row banks and all their column replays use the fixed
$`q_{\rm ch}=2N/n+O(\log n)`$ source core and their actual finite
masks. If $`S_d`$ counts addressed nodes at depth d, their combined
native cost is

```math
T_{\rm cheap+replay}
=O\left(\sum_d S_d+nq_{\rm ch}
              +N'\mathrm{poly}(n)+\mathrm{poly}(n)\right)=O(N).
```

There are only a constant number of global fine phase banks in the
replay, not one per depth. Source preparation and its
inverse, global normalization, and all fine gauge banks each have
$`O(N)`$ cost and occur a constant number of times, even after the two
amplifications. Equation (3) charges the bottom part. This proves the
T count in (1).

The simultaneous active clean reservation consists of a source register
of $`\log_2M=O(\log n)`$ bits, a power register and its increment
scratch, a multiplier and full-product register of $`O(\log n)`$
bits, $`\log_2(8n)`$ rare bits, $`\lceil\log_2R\rceil`$ history
bits, arithmetic/predicate scratch of $`O(\log n)`$ bits, and a fixed
number of distinct scalar, underflow, preparation, normalization, and
fine-gauge flags, including the separate fine compiler and the returned
node-type flags. This totals $`O(\log n)`$ with a fixed constant,
covered by K. Long logical suffix predicates use the existing returned
borrowed construction, not an n-bit clean register.

At a fixed depth, carrier and splitter predicates use the unchanged
logical suffix, fixed flags, and the returned node-type table just
charged. Quarter-turn and splitter controls
therefore use polynomially many exact gates and returned borrowed
predicate work, rather than initialized logical-endpoint metadata.

During fine calls reserve one shared core of $`q+1=N+c_0+1`$ arbitrary
dirty wires. The largest fine address has at most n+1 bits because its
mode dimension is at most $`2N`$. Self-borrowing applies once
$`q+1\ge2(n+1)`$. All bounded extra controls and exact center predicates
reuse returned core helpers in completed subroutines. A fixed number of
external dirty helpers, if used, fits in the remaining $`n+6-c_0`$
wires for sufficiently large n. Coarse reflection routing uses
$`O(\sqrt N+\mathrm{poly}(n))`$ dirty space instead, and the
dyadic kernels use at most $`N/n+O_{A_\star}(n^3)`$ dirty source
wires, which is below $`N/2`$ after enlarging the threshold.
The cheap source uses $`2N/n+O(\log n)`$ wires, also below $`N/2`$.
Their short-output metadata queries reuse returned selectors inside
that core. These reservations are sequential and reuse the same pool.
No active-suffix bit is also
counted as dirty storage. Thus the literal maximum is $`N+n+7`$.
The coefficient of K can be chosen universally: for each fixed
$`A_\star`$, the $`O(\log A_\star)`$ metadata margin is at most
$`O(\log n)`$ beyond an $`A_\star`$-dependent threshold. Only that
threshold and the count constant depend on $`A_\star`$.

All finitely many n below the chosen threshold have fixed dimension and
fixed required accuracy. The existing direct determinant-one,
literal-phase finite-net fallback compiles every frame in those compact
families with the same two flags and dirty pool. The finite maximum of
T/N is absorbed into the constant depending on $`A_\star`$. This is an asymptotic
theorem, not a proposed practical crossover threshold.

## 8. Scope and verification

**Nodewise square-root-window corollary.** The theorem covers arbitrary
independent bins $`0\le E_v\le N/\sqrt n`$ on a fully connected
tree, because $`\sqrt n\ge\log_2n`$ for the permitted n. This is
the case $`A_\star=1`$.

**Mixed-height corollary.** The earlier
[mixed-height family](PRECISION_GAUGE_COMPILER.md#mixed-heights-and-the-limit-of-rescheduling)
has $`h_r=n/2^{r+1}`$ and $`E_r=N/h_r`$, with
$`E_{min}=2N/n`$ and $`E_{max}=N/\sqrt n`$. All of its zero
crossings are allowed by the theorem. The native T count is therefore
$`O(N)`$ uniformly in its growing number of classes. Its unchanged-bank
scheduling proposition remains a statement about those original banks;
the present circuit uses a common source and different local kernels.

The exponentially small shift eigenvalue in (9) permits endpoint
accuracy with only $`M=4n^2`$ initialized source modes. Nodewise scale
changes on expensive nodes are encoded by the short multiplier $`k_v`$
and the charged dyadic remainder. The cheap nodes retain their actual
native rejected operators in the same single-feature gauge. The upper
proof needs no class preservation. The threshold $`2N/n`$ controls the
expensive Taylor order and cheap source sum; it is no longer a lower
endpoint of the theorem. The upper envelope
$`A_\star N/\log_2n`$ pays for the bottom K depths. Frames whose
remaining bottom precision trace is uncontrolled are outside this
theorem.

The [exact helper](../../compiler_robust_hopf/heterogeneous_precision_source.py)
records finite source residuals, dyadic Cayley coefficients, and symbolic
parameter margins, including the common-source spacing and multiplier
ledger. The
[focused tests](../../tests/test_heterogeneous_precision_source.py)
check finite shift residuals including shifts beyond the cap, exact
dyadic quotient/remainder coefficients, the retained linear precision
charge, truncated kernel columns, and exact large-integer metadata.
Full-input checks cover separate scalar/underflow flags, multiplier
cleanup on occupied inputs, interleaved node scales, and no-wrap history
with uniform attenuation. Mixed-chart fixtures compose the actual
three-port row words with noncommuting dirty rejection operators,
complex splitter phases, and an exact-zero variant. They check the
mixed gauge, the common-source intertwiner on occupied features and
dirty inputs, actual column replay, the P/Q overlap, and full recovery.
Negative checks detect dropped splitter phase, scalarized rejection,
reversed replay order, and a reflection that tests the arbitrary dirty
register. Earlier class-dependent source fixtures and the all-node
zero-splitter fixture remain scoped to their earlier charts. The mixed
chart here uses $`B_v=I,F_v=0`$ and no splitter at exact-zero nodes.
Those bounded checks verify the displayed interfaces; they do not
establish the asymptotic count, replace the source-preparation charge,
or constitute an unrestricted native frame emitter. The native theorem
follows from the exact kernels, contractive recurrence, fixed fine-call
count, and workspace ledger above.
