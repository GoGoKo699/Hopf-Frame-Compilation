# One charged binary gradient across complete-column recovery

[Complete-frame assembly](SUBLINEAR_T_DEPTH.md) · [Coherent descendant source](COHERENT_DESCENDANT_LOADER.md) · [Exact shallow dirty queries](SHALLOW_DIRTY_QUERIES.md)

The binary gradient can remain prepared throughout the early frame's
complete-column recovery. Every source coefficient bank, calibration,
and outer damping operation is compiled into exact modular translations
on that same register. A paid two-query construction implements addition
with logarithmic T-depth and one returned zero register. Its exponential
precision cost fits the frame budgets in a sufficiently small absolute
linear-precision range.

The recovery reflection omits the retained gradient. The full source
comparison and recovery proof take place inside its invariant eigenspace;
the native preparation and its actual inverse are charged once each at
the two exterior boundaries. This is a source-preparation saving, not an
improvement to the complete-column normalization. The final layers have their separately charged
[grouped terminal schedule](GROUPED_TERMINAL_DEPTH.md).

## 1. Registers and the early-segment contract

Use the complete-frame split with h active logical bits X and M protected
terminal bits Z, and put $`Q=2^h`$. The initialized external flag H records
the original predicate $`[Z=0^M]`$ before Z is used as workspace. On
$`H=1`$, the designated Z registers are zero; on $`H=0`$, their state
is arbitrary, including references. The early target is

```math
W_E=W_h\otimes\Pi_Z+I_X\otimes(I_Z-\Pi_Z).
```

Write $`K=2^{\lceil\log_2(h+1)\rceil}`$ for the padded tag dimension,
and choose the complete-frame recovery parameters

```math
k=\left\lceil\frac{\pi\sqrt K}{4}\right\rceil,\qquad
\nu=2k+1,\qquad
\vartheta=\frac{\pi}{4k+2},\qquad
\gamma=\sqrt K\sin\vartheta\le1.
```

For precision $`\eta=2^{-L}`$, $`L\ge6`$, set

```math
\delta=\frac{\eta}{64\nu},\qquad
r=L+\lceil\log_2\nu\rceil+16,\qquad D=2^r.
```

Then $`\nu\le5\sqrt{n+2}`$ for $`h\le n`$, and
$`r\le L+\tfrac12\log_2(n+2)+20`$. Reserve three separate r-bit
registers inside Z: the angle program B, the gradient F, and the returned
addition register C. The source's $`q=O(\log r)`$ index bits, its inner
attenuation bit, the $`s=\lceil\log_2(h+1)\rceil`$ tag, and one outer
attenuation bit are separate. The simultaneous conditional-zero footprint
is therefore

```math
3r+q+s+2=O(L+\log(n+2)).
```

All additional query helper pools and temporary mask outputs are
supplied arbitrary dirty work. They are disjoint from every live address,
target and conditional-zero register, and return exactly after each
completed query or mask echo. The designated B, C, tag and logical data
outputs retain the allocations stated here and in the complete-frame
source. No operation retests Z while this protected work is occupied.

## 2. A paid logarithmic-depth modular adder

The shallow dirty-query theorem implements every fixed Boolean output
table exactly, not just a promised state map. Apply it to the two
r-bit-output functions on the 2r-bit address $`(m,x)`$,

```math
f_+(m,x)=x+m\pmod D,\qquad
f_-(m,x)=x-m\pmod D.
```

Their native XOR queries are $`Q_+`$ and $`Q_-`$. Define the addition
word chronologically by

```math
\mathsf A:\quad Q_+(B,F;C),\ \mathrm{SWAP}(F,C),\ Q_-(B,F;C).
```

On $`|m,x,0\rangle_{BFC}`$ the three stages give

```math
(m,x,0)\longmapsto(m,x,x+m)
\longmapsto(m,x+m,x)\longmapsto(m,x+m,0),
```

with all arithmetic modulo D. This holds for every m and x, hence for
arbitrary superpositions and external references. B is preserved on
the entire input space. The $`C=0`$ subspace is exactly invariant, and
the actual reversed native word implements subtraction there:

```math
\mathsf A^\dagger|m,x,0\rangle=|m,x-m,0\rangle.
```

On occupied C the word is still a full-space unitary permutation; no
addition promise is made there. In particular, its internal queries and
swap are unconditional. Exact inactive behavior is supplied by the
complete echo below, not by incorrectly zeroing query rows while leaving
an unconditional swap.

The query theorem gives, including both queries and all elementary swaps,

```math
T,w=O((r+1)^4 2^r),\qquad
G=O(r4^r+(r+1)^4 2^r),\qquad
D_T=O(\log(r+2)).
```

Here w is additional arbitrary dirty width. There is no extra initialized
helper beyond C. Each query and every actual inverse returns its dirty
pool on the full input space; the same pool can be reused after every
completed query. The tables have $`4^r`$ rows each, and their costs are
retained explicitly. This is not an uncharged arithmetic oracle.

## 3. Exact inactive guards and a preserved gradient

For a preserved bit c, let $`\mathsf J_c`$ complement every F bit
controlled by c. It is a Clifford CNOT bank; at $`c=1`$ its arithmetic
action is $`x\mapsto-x-1\pmod D`$. For a target t put
$`\mathsf V_t=\mathsf J_t\mathsf A\mathsf J_t`$. On $`C=0`$ this
shifts F by $`(-1)^t m`$ and returns C. The chronological word

```math
\mathsf V_t,\ \mathsf J_H,\ \mathsf V_t^\dagger,\ \mathsf J_H
```

has operator $`\mathsf J_H\mathsf V_t^\dagger\mathsf J_H\mathsf V_t`$.
On $`H=0`$ it is exactly identity for arbitrary B, F, C and t, since
the entire actual word cancels against its adjoint. On $`H=1,C=0`$,
each completed adder preserves C zero, and complementation conjugates
translations to their inverses. The word therefore shifts F by
$`2(-1)^t m`$. The same echo with $`\mathsf V=\mathsf A`$ gives the
unsigned shift $`2m`$. Both versions preserve B and C exactly on the
active zero-C sector, for every gradient input.

Define the positive-exponent gradient

```math
|\Gamma_r\rangle=D^{-1/2}\sum_{x=0}^{D-1}e^{2\pi i x/D}|x\rangle.
```

Translation by d has eigenvalue $`e^{-2\pi i d/D}`$ on this state.
Thus the signed echo gives
$`\mathrm{diag}(e^{-i\theta_m},e^{i\theta_m})`$ on t, with
$`\theta_m=4\pi m/D`$. Let $`Q_Y=SH`$, so
$`Q_YZQ_Y^\dagger=Y`$. Chronological conjugation by
$`Q_Y^\dagger`$, the echo, and $`Q_Y`$ implements the literal
$`R_y(\theta_m)=e^{-i\theta_mY}`$. There is no row-dependent scalar.

Load m into zero B by the exact H-enabled query, apply the complete
rotation macro, and reverse the actual query while its address is
unchanged. B returns exactly. On $`H=0`$ the query row is zero and the
middle echo is full-space identity, so both target basis changes cancel
even on occupied B, F and C. No preparation gate occurs inside a macro.

All completed active macros preserve the subspace

```math
\mathcal K=
\mathrm{span}\{|\Gamma_r\rangle_F|0^r\rangle_C\}
\otimes\mathcal H_{\rm remaining}.
```

This invariance holds for every state of the remaining registers,
including occupied programs and source work, not just source columns.
Since the restrictions are finite-dimensional unitaries, their actual
inverses preserve the same subspace. The property of C also holds with
an arbitrary F input, covering all subsequent preparation error.

## 4. One short table for every coefficient bank

The descendant source has $`J=2^q=O(r)`$ expansion terms and coefficient
preparation

```math
U_\chi=\bigotimes_{\ell=0}^{q-1}R_y(\theta_\ell),\qquad
\theta_\ell=\arctan\!\left(\beta^{2^\ell/2}\right),\qquad
\beta=4/5.
```

In the Y eigenbasis this whole bank is one diagonal:

```math
U_\chi=Q_Y^{\otimes q}\,
 \mathrm{diag}_x(e^{-i\Theta_x})\,
 (Q_Y^\dagger)^{\otimes q},\qquad
\Theta_x=\sum_{\ell=0}^{q-1}(-1)^{x_\ell}\theta_\ell.
```

For every x, choose an integer $`m_x\in\{0,\ldots,D/2-1\}`$ with
certified circular error

```math
\mathrm{dist}_{2\pi}(\Theta_x,4\pi m_x/D)\lt4\pi/D.
```

This selection terminates: the grid has spacing $`4\pi/D`$, so some
point is within $`2\pi/D`$, and finite interval certificates for the
looser strict test can be dovetailed. No exact midpoint decision is
required. Compute the whole angle sum classically when constructing this
finite table; its quantum evaluation remains a charged query.

Apply $`(Q_Y^\dagger)^{\otimes q}`$, load $`m_x`$ into B by its
H-enabled query, use the unsigned guarded echo, unload B with the actual
query inverse, and apply $`Q_Y^{\otimes q}`$. The index bits x remain
unchanged between query and unquery. On $`\mathcal K`$ with B initialized
to zero, this has operator error less than $`4\pi/D`$ from $`U_\chi`$
on every index input, tensored with identity on the remaining registers.
Every completed word and its actual inverse preserves each computational
B value exactly. On occupied B it instead implements an XOR-offset
program, so no desired-angle accuracy is claimed on that larger space.
There is no additional factor q, because the complete eigenphase was
rounded once. Every phase, including a common scalar, is retained.

The table has $`O(J)=O(r)`$ rows, r output bits and
$`O(\log(r+2))`$ address bits. Its T-count, dirty width and elementary
Clifford count are polynomial in r; its T-depth is
$`O(\log\log(r+4))`$. The arithmetic adds $`O(\log(r+2))`$ depth.
The inverse coefficient bank is the actual reversed native word, not a
freshly rounded negative-angle table.

The source-label-dependent attenuation uses the signed echo from
Section 3 and the existing preserved-tag/prefix calibration table. The
outer damping rotation needs a stronger occupied-program contract. For
its fixed rounded angle word $`m_\gamma`$, hardwire that constant into
the arithmetic tables $`f_{\pm,\gamma}(x)=x\pm m_\gamma\pmod D`$.
The same compute/swap/uncompute construction and signed echo then ignore
B entirely, and give $`R_y(4\pi m_\gamma/D)\otimes I_B`$ on
$`\mathcal K`$ even when B is occupied. These smaller r-address tables
fit the generic arithmetic bounds. Loading a fixed word by XOR into an
occupied B would not give this full-operator damping contract.

Thus no remaining coefficient or damping operation is synthesized to
precision eta inside each recovery call. The original source angles and
literal inner amplification formula are those of the descendant proof;
the outer damping angle satisfies $`\cos\theta_\gamma=\gamma`$.

## 5. Source comparison and recovery inside the invariant sector

The source's Boolean phases, marker packing, Hadamards and zero tests
do not act on F or C. All inner reflections act as identity on these
registers; they retain their original index, attenuation and local-target
tests. Combined with Sections 3–4, every completed native source and its
actual inverse restrict to a unitary on $`\mathcal K`$. Supplied dirty
systems return exactly, so this restricted unitary tensors with identity
on those systems and on external references.

Let $`\widehat R_\Gamma`$ be the restricted complete-column source,
including marker packing. Let $`J_0`$ initialize its tag and other source
work, and let R be the ideal complete-column isometry. The latter has
the row-overlap identity $`P^\dagger R=W_h/\sqrt K`$, where P is the
tag-Hadamard row isometry. These are the same disjoint tag/prefix columns
as in the complete-frame proof.

Choose the common power-of-two expansion length J so that
$`2\beta^J\le\delta/2`$. Its size is $`O(r)`$, and the strict
residual-selection rule from the descendant source remains terminating.
One flagging word contains two coefficient-bank occurrences and one
calibration. Two inner amplification rounds contain five forward/inverse
flagging words, hence at most fifteen such angle banks. Our choice of r
gives

```math
\rho=4\pi2^{-r}\le\frac{\pi\delta}{256}\lt\frac\delta{64}.
```

For these angle-error comparisons, restrict further to
$`\mathcal K_0=\mathcal K\cap\{B=0\}`$. Every completed programmed
macro and its actual inverse preserves each B value: its query address
is unchanged, the adder preserves B, and the actual query inverse removes
the loaded word. This holds even for arbitrary gradient inputs and
occupied C. The other source macros do not change B, and the desired
continuous-angle comparison operations also preserve B zero. Thus
$`\mathcal K_0`$ is an exactly invariant subspace for the source-angle
hybrid, on which both coefficient and calibration error bounds apply.
Unitary telescoping there bounds their total angular error by
$`15\rho\lt\delta/2`$. The data error is at most
$`\delta/2`$. Preserved tag/prefix sectors make the full source norm
the maximum branch error, not a sum over columns. Consequently

```math
\|\widehat R_\Gamma J_0-R\|\le\delta.
```

Apply the complete-frame dimension-independent isometry-extension lemma
on the larger reduced mode space $`\mathcal K`$, which permits occupied
B. Only the displayed source-isometry bound on initialized anchors is
needed; intended-angle accuracy on occupied B is not assumed. The lemma
supplies a mathematical unitary
$`C_\Gamma`$ with $`\|C_\Gamma-I\|\le6\delta`$ and

```math
R_{*,\Gamma}=C_\Gamma\widehat R_\Gamma,\qquad
R_{*,\Gamma}J_0=R.
```

The correction is never implemented. If a physical extension is desired,
extend it by identity on $`\mathcal K^\perp`$ and on $`H=0`$. It
preserves $`\mathcal K`$, has the same norm bound, and makes no desired
action claim about the other gradient sectors. Actual inverse source
words obey the same restricted full-operator comparison.

For physical outer recovery, let $`\Pi_{\rm work}`$ test all initialized
live work except F: tag, index, attenuation, B, C and the outer damping
bit. It excludes X, H and the supplied dirty pool. Use

```math
S_H=I-2[H=1]\Pi_{\rm work},\qquad
\widehat A_H=Z_H\widehat U_HS_H\widehat U_H^\dagger S_H.
```

F is not tested for zero or for any other state. The same is true of
unused Z spectators. On $`\mathcal K`$, the C-zero factor is identity,
and this reflection restricts to exactly the initialized-work reflection
needed by the complete-frame recovery. It is an exact returned-helper
predicate on $`O(r+\log(n+2))`$ bits, with logarithmic T-depth and
charged polynomial count. It neither projects onto a supplied phase
state nor re-prepares the gradient.

The outer damping word has restricted full-operator error at most
$`\rho\lt\delta`$. The ideal comparison
$`U_{*,\Gamma}=G\,H_D^{\otimes s}R_{*,\Gamma}`$ therefore satisfies

```math
J_0^\dagger U_{*,\Gamma}J_0
=\sin\vartheta\,W_h,\qquad
\|\widehat U_\Gamma-U_{*,\Gamma}\|\le7\delta.
```

Here $`J_0`$ is enlarged to initialize the outer attenuation bit.
The literal recovery formula $`(-U_*SU_*^\dagger S)^kU_*J_0=J_0W_h`$
holds on this invariant sector with positive phase. The physical $`Z_H`$
supplies the active minus; every factor is identity on arbitrary work
when $`H=0`$. Telescoping across the $`\nu=2k+1`$ actual forward/inverse
source occurrences gives body error at most

```math
7\nu\delta=7\eta/64.
```

No source call is assumed to return its tag, index or attenuation exactly.
Those errors are covered by the full restricted-unitary comparison.
Only the exactly invariant gradient-and-C sector is used as a promise.

## 6. One native preparation boundary pair

For $`\alpha_j=2\pi2^j/D`$, use the product comparison preparation

```math
P_0=\bigotimes_{j=0}^{r-1}(e^{-i\alpha_jZ/2}H),\qquad
P_0|0^r\rangle=\lambda_r|\Gamma_r\rangle,\qquad
\lambda_r=e^{-i\pi(1-1/D)}.
```

Approximate its r determinant-one rotations in parallel, using the
retained one-qubit synthesis premise, with total literal operator error
$`e_P\le\eta/64`$. This gives a native P with

```math
T,G=O(r[L+\log(r+2)]),\qquad
D_T=O(L+\log(r+2)).
```

Apply P unconditionally to F after computing H, then the entire native
recovery body V, then the actual $`P^\dagger`$. Let $`E_0`$ embed the
active input with F and all other work zero, and let $`E_\Gamma`$ use
$`\lambda_r\Gamma_r`$ on F instead. For arbitrary logical, dirty and
reference inputs,

```math
\|PE_0-E_\Gamma\|\le e_P,\qquad
\|VE_\Gamma-E_\Gamma W_h\|\le7\eta/64.
```

Unitarity and the actual exit inverse give

```math
\|P^\dagger VPE_0-E_0W_h\|
\le2e_P+7\eta/64\le9\eta/64\lt\eta/4.
```

The common preparation phase cancels. This inequality includes all
gradient, source-work and eventual flag leakage. It has no factor equal
to the number of recovery calls. On $`H=0`$, V is exactly identity on
the full arbitrary-work space, including occupied C, so
$`P^\dagger VP=I`$ there with no preparation assumption. Finally reverse
the original exact H predicate. Applying that same inverse to actual and
ideal outputs preserves the combined norm bound and returns H on the
ideal output. No physical reset is used.

## 7. Resources and the uniform-precision range

The source has a constant number of large phase/calibration queries and
short coefficient-bank queries per call. Its long tables retain at most
$`O((h+1)Qr)`$ padded rows; neither F nor C is appended to their address.
The original tag masks, Boolean phases and reflections have depth
$`O(\log(n+2))`$. Every new angle bank adds
$`O(\log(r+2))`$ depth. Including the single exterior preparation pair,
the early segment therefore has

```math
D_{T,E}=O\!\left(L+\log(r+2)
 +\sqrt n[\log(n+2)+\log(r+2)]\right).
```

All constants are absolute. For $`u=n+L+2`$, the native query bounds and
$`O(\sqrt n)`$ recovery calls give a conservative common majorant

```math
T_E,w_E\le C u^6(2^{h/2}+2^r+1),\qquad
G_E\le C u^6(2^h+4^r+1).
```

The width bound deliberately overestimates repeated calls: the same dirty
pool is reused after exact return. The terms $`2^r`$ and $`4^r`$ are the
paid arithmetic tables; they are not suppressed by the early-frame factor
$`2^{-M}`$ and cannot be discarded at unrestricted precision.

To fit the precision-uniform frame budgets, choose an absolute A large
enough for all conditional reservations and polynomial suppression, and
take $`M=\lceil A[L+\log_2(n+L+2)]\rceil`$. Choose an absolute
$`c\gt0`$ small enough that, throughout $`6\le L\le cn`$ and all
sufficiently large n, $`M\le n/2`$ and $`r\le n/4`$. Such a c and
an absolute threshold exist by the explicit r bound in Section 1.

The $`2^{h/2}`$ and $`2^h`$ terms are suppressed by M as in the
complete-frame precision ledger. The standalone arithmetic terms fit
independently: for every fixed degree d,

```math
\frac{n^d2^r}{\sqrt{N/n}}
\le n^{d+1/2}2^{-n/4}\longrightarrow0,\qquad
\frac{n^d4^r}{N}
\le n^d2^{-n/2}\longrightarrow0.
```

Thus all early queries and arithmetic fit every
$`b\ge\sqrt{N/n}`$ in that large-n branch, while preserving the
desired $`O(\sqrt{NL})`$ T and $`O(NL)`$ Clifford allowances. The
one-time preparation count is polynomial and fits the same budgets.
The finite-dimension and smaller-width cases require the retained
compiler splice in the complete-frame proof, with absolute constants.

For $`L\le cn`$, the displayed early depth simplifies to
$`O(L+\sqrt n\log(n+2))`$. The separate
[grouped terminal construction](GROUPED_TERMINAL_DEPTH.md) supplies the
remaining layers and their precision-dependent depth bound.
The full theorem, width splice, and simultaneous count statements belong
to the complete-frame assembly. No unrestricted precision range,
efficient classical preprocessing, general native emitter, or optimal
large-width depth is asserted here.

The [bounded shared-gradient checks](../../tests/test_shared_gradient_recovery.py)
test literal arithmetic, occupied-work guards, coefficient phases,
program-dependent versus fixed damping, source-boundary error, and the
reflection mask. Their small recovery-shaped matrices are local
diagnostics; the invariant-sector argument and uniform resource absorption
above supply the analytic proof.
