# Coarse flags with a final-accuracy sparse correction

[R94 baseline](PORTABLE_WORKSPACE_CATALOGUE.md) ·
[Causal compact coins](CAUSAL_COMPACT_COIN_DEPTH.md) ·
[Shared binary gradient](SHARED_GRADIENT_RECOVERY.md) ·
[Sparse row and column sources](BOUNDED_DENSITY_DEPTH.md#5-certified-schur-weights-and-the-sparse-source-interface)

**Opening component of [R95](LOGARITHMIC_PRECISION_DEPTH.md).** This chapter changes
an early compiler's precision allocation. Its coarse mathematical tree is
rounded only to inverse-polynomial accuracy in the local height, while its
physical circuit and a sparse additive correction are implemented at the
requested final accuracy. It does not assume that a low-accuracy native
word returns its work exactly.

## 1. Prefix-selected early interface

Let the global logical register have n bits. An unchanged prefix P has
$`n-m`$ bits, the local register has m bits, and its bottom r original
logical bits form Z. Put $`a=m-r`$. The exact early factor is

```math
E=\sum_P|P\rangle\langle P|\otimes
 \left(W_{P,a}\otimes\Pi_Z+I\otimes(I-\Pi_Z)\right),
\qquad \Pi_Z=|0^r\rangle\langle0^r|.
```

The prescribed real frames include every column, signs and singular pairs.
Let $`l\ge6`$ be the requested stage precision and define

```math
K=m+l+2,\qquad k=\lceil\log_2 K\rceil,\qquad v=l+k^3.
\qquad(1)
```

The construction below uses $`O(v)`$ conditionally initialized bits, a
constant number of sequential P-dependent table positions, and one external
clean qubit. All constants and polynomial degrees are independent of m, l,
P and the frame. Once those constants are fixed, choose
$`r=\lceil Av\rceil`$, for a sufficiently large absolute A, and use this
stage only when $`r\lt m`$.

Save $`H=[Z=0^r]`$ in the external clean bit before using Z as work. Every
completed body word is exactly identity on H zero, also on occupied work.
On H one its complete initialized-isometry error is at most
$`2^{-l}/16`$. The actual saved-predicate inverse is applied last. Dirty
query pools are arbitrary, may be entangled with references, and return
exactly after each completed query. Approximate source and gradient return
is included in the stage error.

The target application is after R94's first uniform R89 stage, where
$`m=O(\ell^3)`$, $`l=L+O(1)`$, $`6\le L\le3\ell`$ and
$`\ell=\log_2(n+2)`$. There (1) is $`O(\ell)`$ after one absolute
cutoff. This is sufficient for the complete-depth assembly; no
linear-in-$`l+\log m`$ reserve is asserted for all m,l.

## 2. A coarse exact tree and a small positive residual norm

For each fixed P, use the certified finite rounding of
[R86 Section 2](INTERVAL_DILATION_DEPTH.md#2-finite-rounding-and-the-atom-embedding)
on every early node. Choose an absolute

```math
\epsilon_0=\frac1{64(m+1)^2}
```

as a bound on each two-dimensional rotation error. A dyadic probability
precision $`\tau_0=O(\log_2(m+2))`$ suffices: take the R86 bound
$`2\,2^{-\tau_0/2}\le\epsilon_0`$. Its clipping makes the rounded
magnitudes positive, and its finite sign choices do not require deciding
whether an original coefficient is zero. All nodes below originally zero
incoming amplitudes are retained. The rounded pair and its finite signs
define an exact mathematical real Hopf frame $`W^0_P`$.

Put $`D_P=W_{P,a}-W^0_P`$. Unitary layer telescoping gives

```math
\max_j\|(D_P)_{\cdot j}\|_2\le a\epsilon_0\le m\epsilon_0.
```

Both frame matrices have the same allowed support: a computational leaf
meets only the state column and its ancestor marker columns. Thus every row
of $`D_P`$ has at most $`a+1\le m+1`$ nonzero entries. For every vector z,
Cauchy--Schwarz in each row gives

```math
\begin{aligned}
\|\,|D_P|z\|_2^2
&\le(m+1)\sum_j|z_j|^2\sum_x|(D_P)_{xj}|^2,\\
\|\,|D_P|\|
&\le\sqrt{m+1}\max_j\|(D_P)_{\cdot j}\|_2
 \le m\sqrt{m+1}\epsilon_0\lt \frac1{16}.
\end{aligned}\qquad(2)
```

Here $`|D_P|`$ is entrywise absolute value. Equation (2) is uniform over
singular angles and requires no relative branch-accuracy promise.

The physical implementation of $`W^0_P`$ in Sections 3–5 is at FINAL
accuracy. Its native leakage is therefore an ordinary final-error term;
it is not silently included in the sparse matrix $`D_P`$.

## 3. Causal compact coins at coarse geometric precision

Use the exact positive magnitude tree underlying $`W^0_P`$ in R89's
local construction. Every probability and compact atom has geometric
precision $`\tau_0`$, while every physical source has final precision.
Keep the exact sign diagonals outside this positive construction; each is
a charged Boolean-phase query on P and the appropriate logical label.

The retained geometric bins satisfy

```math
p=O(\tau_0\log_2(m+2))=O(\log_2^2(m+2)).
\qquad(3)
```

The free donor reserve in the bin construction excludes every live source,
label and gradient register. Fixing a sufficiently large constant multiple
of v accommodates that reserve as well. The original tagged source,
endpoint inverse, fine product-coin erasure and endpoint reload have a
constant number of P-dependent source positions. Their sources are
implemented at final tolerance, and their actual inverse words are used.
The selected contiguous routers retain the logarithmic-depth schedule of
R89/R91; their parameters, donors and output banks are charged separately.

For block i of height $`q_i`$, the compact fixed coin has
$`2^{q_i}(q_i\tau_0+1)`$ positive atoms. Its coin register uses
$`q_i+O(\log K)`$ bits. The R89 staged shear compressor, one joint rank
source and three-query causal matching are unchanged. Hash and shift
catalogues are chosen separately for each immutable P and uniformly over
that P's internal original prefixes. Their selectors have $`O(\log K)`$
bits. Matching remains causal at the original prefixes, which is the
hypothesis preserving every coarse band; no arbitrary global
weight-preserving matching is substituted.

All fixed code padding, hash records, matching indices, rank-window
metadata and live router parameters require $`O(p\log K)=O(k^3)`$
additional conditional coordinates. The single joint rank source and the
constant number of other source banks use $`O(l+\log K)`$ work. The only
remaining obstacle in the old allocation is the separate final-precision
bank for each of p simultaneous reference sources. The next section
replaces those banks.

## 4. One gradient for all fixed reference sources

### 4.1 Explicit simultaneous grammar

Let $`C_i`$ be the i-th compact coin and $`\gamma_i`$ its fixed positive
reference state. These references depend on the block sizes and
$`\tau_0`$, not on P or the node probabilities. For every i use the
[real Boolean-phase source grammar](COHERENT_DESCENDANT_LOADER.md#3-constant-query-flagging-and-literal-inner-amplification)
with one common padded expansion length

```math
J=O(l+\log_2(p+2)),\qquad q=\log_2J,
```

where J is a power of two. Each local initialized bank $`A_i`$ has its
q-bit coefficient index and one attenuation bit. An ordinary source has
five appearances of its flagging word or actual inverse, separated by
exact local anchor and success reflections. The simultaneous sources use
this same constant chronological grammar. Their Boolean phases and
reflections act on disjoint $`C_i,A_i`$ and private returned dirty pools.

Reserve ONE g-bit binary gradient F, one g-bit program B, and one g-bit
returned addition register C, where

```math
g=l+O(\log_2(p+2))+c_g.
\qquad(4)
```

The constant is fixed by the error allocation below. The paid arithmetic,
its actual inverse, the unsigned guarded translation echo and the
coefficient eigenbasis are exactly
[shared-gradient Sections 2–4](SHARED_GRADIENT_RECOVERY.md#2-a-paid-logarithmic-depth-modular-adder).
No gradient or initialized addition bank is supplied externally.

At a coefficient-bank or attenuation position, first conjugate the
participating targets into their Y eigenbasis. For a local coefficient
bank, round its WHOLE eigenphase $`\Theta_i(x_i)`$ once to a multiple
of $`4\pi/2^g`$, as in the shared-gradient coefficient-bank construction.
For an attenuation bank, round its one angle. Let $`m_i(x_i)`$ be the
resulting signed integer eigenphase. There is one rounding error per
local bank, with no coefficient-index-length factor. Classical
preprocessing computes the finite combined table

```math
M(x_1,\ldots,x_p)=\sum_i m_i(x_i)\pmod{2^g}.
\qquad(5)
```

Load M into zero B by an exact H-enabled XOR query; apply one unsigned
guarded translation echo; reverse that query while all its addresses
are unchanged; undo the target basis changes. This is one completed
rotation bank. Its actual inverse is the reversed word.

On the exact gradient

```math
|\Gamma_g\rangle=2^{-g/2}\sum_z e^{2\pi iz/2^g}|z\rangle,
```

the bank applies the product of the individually rounded eigenphases.
The integers in (5) are summed EXACTLY. Rounding a whole joint phase
after summation is not used: preserving this exact tensor factorization
is needed for the reference-flag argument.

The combined table reads at most $`O(pq+p)=O(k^3)`$ local bits and has
g outputs. It does not read P. Each reference Boolean phase reads its own
coin and local index; anchors test its own coin and every bit of $`A_i`$.
They exclude the shared F,B,C registers.

### 4.2 The invariant mode space and complete band readout

Let

```math
\mathcal K_0=
|\Gamma_g\rangle_F\otimes|0^g\rangle_B\otimes|0^g\rangle_C
\otimes\mathcal H_{\mathrm{coins,local\ work,other\ systems}}.
\qquad(6)
```

Every completed rotation bank and its actual inverse preserves
$`\mathcal K_0`$ for ALL states of the coins and local source work.
The addition word preserves C zero for every gradient input; each
completed program query/body/unquery preserves each computational B
value; translations preserve $`\Gamma_g`$. Other source operations do
not touch F,B,C. Thus the restriction of the physical reference word to
(6) is a genuine tensor product $`\bigotimes_i\widetilde U_i`$ on
$`C_i,A_i`$, with arbitrary dirty systems as identity factors.

No desired-angle accuracy is asserted on occupied B. It is unnecessary:
(6) is an invariant finite-dimensional mode space, and both the nearby
completion lemma and all comparison algebra can be performed there.

Assign each local reference isometry error at most
$`\varepsilon_i=\delta_{\rm ref}/(p+1)`$. Increasing the absolute
constant in (4) and J supplies this error. The nearby-completion lemma
within (6) gives ideal comparison unitaries $`\widehat U_i`$ with
exact zero columns $`\gamma_i|0^{A_i}\rangle`$ and

```math
\|\widetilde U_i-\widehat U_i\|\le6\varepsilon_i.
```

The per-block computational nonzero test includes $`C_i`$ AND every bit
of $`A_i`$. Fixed encoding padding is zero on the ideal matched image
and untouched by the references. R89's complete-band identity therefore
maps band i to the last nonzero augmented block i, independent of the
other columns chosen for $`\widehat U_i`$.

After the ideal compact embedding, prepare the common gradient by a
native $`P_F`$, apply the simultaneous reference inverses, XOR the
classical coarse label into its tag, apply the actual forward reference
word, and apply $`P_F^\dagger`$. The XOR is guarded by H and leaves
F,B,C unchanged. The original compact embedding is then reversed.
The coarse spectral XOR has literal phase +1. The reference sandwich's
body error is at most $`12\sum_i\varepsilon_i`$; if
$`P_F|0\rangle`$ has error $`e_F`$ from its ideal gradient including
its common phase, its boundary pair contributes at most $`2e_F`$.
The common preparation phase cancels through the actual inverse.

On H zero every completed reference macro and the middle XOR is exactly
identity on its full occupied input space. Thus
$`P_F^\dagger I P_F=I`$ there, without an initialization assumption.
The same argument handles prior source leakage: compare the earlier
embedding to its exact isometry and propagate its error through the
actual unitary sandwich. No intermediate reset is used.

### 4.3 Simultaneous initialized work

The concurrent allocation is

| Bank | Size |
|---|---:|
| Permanent code padding, hashes, matching indices, router fields | $`O(p\log K)`$ |
| Reference-local index and attenuation banks | $`O(p\log(l+\log(p+2)))`$ |
| One reference gradient, program and addition bank | $`O(l+\log(p+2))`$ |
| Original/fine/joint source work and short tags | $`O(l+\log K)`$ |
| Free fine-decoder donor reserve | $`O(v)`$ |

The total is $`O(v)`$. These are distinct simultaneous reservations;
none is counted as arbitrary dirty query work. Serial episodes may reuse
returned pools. This replaces the old $`p\,O(l)`$ reference allocation,
while retaining the complete coarse flag rather than merely its root.

## 5. Final-accuracy coarse implementation

Set every nonreference source tolerance, total reference tolerance and
gradient preparation tolerance to a sufficiently small constant multiple
of a requested $`\delta_0`$. R89's initial tagged source, fine erasure,
endpoint operations, compact embedding and final coarse readout use a
constant number of such completed components. Equations (3)–(6) and the
retained actual-inverse inequalities imply

```math
\|U^0_PJ-J(W^0_P\otimes I_{\rm dirty})\|\le\delta_0.
\qquad(7)
```

The exact sign diagonals are composed before and after the positive
coarse word. They preserve (7). All comparison and error statements are
uniform over P, so their coherent direct sum takes a maximum, not a sum.
The effective enable may additionally include a constant-size LCU branch
predicate. The native enabled-query and echoed-mask construction of
[R82 Section 7](BOUNDED_DENSITY_DEPTH.md#7-a-complete-guarded-constant-normalization-early-circuit)
prices this control without generically controlling a large circuit.

## 6. Native sparse residual sources and normalization-two recovery

Write $`K_P=|D_P|`$. By (2), $`\|K_P\|\lt 1/16`$. The positive finite
Schur weights

```math
B_P=\begin{pmatrix}0&K_P\\K_P^{\mathsf T}&0\end{pmatrix},
\qquad \binom uv=(I-B_P)^{-1}\mathbf1
```

are effectively computable with a known strict norm gap. Their Neumann
series is positive and gives $`u_x,v_j\ge1`$,
$`K_Pv=u-1`$ and $`K_P^{\mathsf T}u=v-1`$. Define

```math
L_{xj}=\sqrt{(K_P)_{xj}v_j/u_x},\qquad
R_{xj}=\mathrm{sgn}((D_P)_{xj})
       \sqrt{(K_P)_{xj}u_x/v_j}.
```

The function $`t\mapsto\mathrm{sgn}(t)\sqrt{|t|}`$, with value
zero at zero, is continuous and effectively computable. This formula
requires no exact sign or zero decision. Then

```math
\sum_jL_{xj}^2=1-1/u_x,\qquad
\sum_xR_{xj}^2=1-1/v_j,\qquad L_{xj}R_{xj}=(D_P)_{xj}.
\qquad(8)
```

Use one logical leaf register X, a short depth tag T, two failure bits F
and a complete source bank. The ideal column isometry is

```math
\mathcal R|j\rangle
 =\sum_xR_{xj}|x,t(j),F=00\rangle
   +v_j^{-1/2}|\mathrm{anchor}(j),t(j),F=01\rangle.
```

For a marker j, exact R78 packing exposes its unchanged node prefix and
zero suffix. The generic signed descendant source prepares this whole
normalized vector, including its failure coordinate. The root column
has tag zero, anchor zero, and its residual sum ranges over ALL leaves;
it is not silently set to pure failure. Distinct columns have distinct
preserved tag/prefix sectors, hence are orthogonal.

The ideal row isometry preserves X:

```math
\mathcal L|x\rangle
 =\sum_{j:\ x\text{ in the support of }j}
 L_{xj}|x,t(j),F=00\rangle
 +u_x^{-1/2}|x,T=0,F=10\rangle.
```

Its short source target has at most a+1 ancestor/root tags and the failure
coordinate. The generic real-vector source supplies both isometries at
final tolerance. Orthogonal failure supports and (8) give
$`\mathcal L^\dagger\mathcal R=D_P`$. Their native word is the
actual column source followed by the actual inverse row source on the
same padded source bank. No second initialized leaf index is used.

Reserve disjoint conditional factor banks for this residual word and the
coarse word (7), plus one LCU branch bit. Exact Hadamard PREP, branch
SELECT, and actual PREP inverse have ideal all-work-zero compression

```math
\frac{W^0_P+D_P}{2}=\frac{W_{P,a}}2.
\qquad(9)
```

Every branch enable includes H. On H zero SELECT is identity and the PREP
pair cancels on arbitrary occupied branch and private work. Constant many
source/coarse errors give block error $`\varepsilon\le C\delta_0`$.

Let J initialize ALL private work in this combined word U and put
$`R_H=I-2HJJ^\dagger`$. This reflection excludes P, X and dirty pools.
The normalization-two special case of
[R82's full-isometry recovery](BOUNDED_DENSITY_DEPTH.md#full-isometry-recovery-with-a-linear-error-constant)
uses the physical word

```math
Z_H U R_H U^\dagger R_H U,
\qquad Z_H=(-1)^H.
\qquad(10)
```

At exact compression (9), its accepted amplitude is
$`\sin(3\pi/6)=1`$, its rejection is zero, and $`Z_H`$ corrects the
literal recovery sign. The block-error-to-isometry argument has an
absolute linear constant because the ideal success amplitude is 1/2,
strictly separated from one. Thus (10) has full initialized-isometry
error at most $`C'\delta_0`$. Choose
$`\delta_0=2^{-l}/(16C')`$, increasing source precision by an absolute
constant. The actual saved-H predicate inverse propagates all remaining
leakage unitarily and returns H on the ideal output. This proves the
stage error stated in Section 1.

## 7. Native query envelope and the intended application

All large P-dependent tables read at most

```math
n-r+C v
```

bits. This includes tagged/endpoint/fine sources, compact metadata,
high-part erasure, the joint rank source, causal matching, all large
source anchors, the two sign phases, and both residual sources. The
per-P hash/shear and matching-shift catalogues also read P, even when
their selectors are short; they belong to this same global envelope.
Only a constant number of these positions are sequential. Outputs,
simultaneous instances and standalone masks have a fixed polynomial
bound in K. In particular, source factors exponential in l are included
in $`Cv`$ address bits, not called polynomial in K.

For each simultaneous one-bit-output query use the retained
[R88 read-only sharp interface](PREFIX_SELECTED_COIN_TERMINAL.md#4-parallel-sharp-queries-with-a-read-only-public-address).
A query with s address bits requires its own pool of at least
$`16c_I(s+2)`$ dirty bits. P appears only as read-only Clifford controls
at the non-Clifford boundaries. Splitting a fixed fraction of b among at
most $`K^{d_0}`$ simultaneous instances therefore suffices when

```math
b\ge C K^{d_w}(n+v+2),
\qquad(11)
```

for one fixed absolute degree $`d_w`$. Outputs and helper pools are
disjoint from every live logical, conditional-zero and address register.
This is a sufficient wide-branch condition, not a new literal minimum.

Put $`\overline Q=2^{n-r+Cv}`$. Fixing an envelope degree d after all
native source and output degrees gives the safe large-query bounds

```math
\begin{aligned}
T&\le C K^d(\sqrt{\overline Q}+\overline Q/b+1),\\
G&\le C K^d(\overline Q+1),\\
D_T&\le C\left[\log_2(n+Cv+2)+l+\log_2K
                    +K^d\overline Q/b^2\right].
\end{aligned}\qquad(12)
```

Standalone fixed-reference Boolean tables have at most
$`2^{m+O(v)}`$ rows; combined angle tables have $`2^{O(k^3)}`$ rows;
shared-gradient arithmetic has $`2^{O(l+\log(p+2))}`$ rows. A common
safe count/Clifford/width majorant is $`2^{C(m+v)}`$, with depth
$`O(l+\log_2K)`$. These costs are charged separately from (12), not
suppressed by a nonexistent global-prefix factor. The native gradient
preparation has polynomial count and depth $`O(l+\log K)`$.

Now specialize to $`m\le C_0\ell^3`$, $`l=L+c_0`$ and
$`6\le L\le3\ell`$, with fixed absolute constants. Then
$`v=O(\ell)`$, $`m+v=O(\ell^3)=o(n)`$ and the standalone majorant
is $`N^{o(1)}`$. For every wide-branch allocation
$`b\ge\sqrt N/n`$, these terms and (11) fit after one absolute
cutoff. Fix A after C,d and all clean-reservation constants, so

```math
r-Cv\ge(2d+80)v.
```

Since $`v\ge\log_2K`$, the large-table terms in (12) have summable
inverse-power K slack. They fit $`O(\sqrt{NL}+NL/b)`$ T gates,
$`O(NL)`$ elementary Clifford gates, and

```math
D_T=O(NL/b^2+L+\ell).
```

The remaining bottom height is $`r=O(\ell)`$. This chapter proves only
this early component; a complete improved depth bound additionally
requires a valid final-precision compiler for that bottom. Original-width
fallbacks, the complex leaf phase and the simultaneous complete-circuit
bounds belong to that separate assembly. No new general depth lower
bound, efficient preprocessing, practical crossover or scalable emitter
is asserted here.
