# Coupled final encoders and scoped three-layer boundaries

[Complete Pauli transport](MOD3_ZETA_TRANSPORT.md#pauli-transport-interface) ·
[Joint two-coordinate source](WEAK_PAULI_ORBIT_COMPRESSION.md#joint-two-address-source) ·
[Final-block boundary](WEAK_PAULI_ORBIT_COMPRESSION.md#final-band-commutation-boundary)

A single native band can mix every address coordinate with the probe and
arbitrary dirty tags while keeping the conjugated target's columns exactly
flat. This supplies an explicit family outside the earlier final-block
commutation condition. It is a final encoder, not a complete source.
At three or more address bits, an analytic obstruction now rules out every
one-band middle completion of this family. The argument covers any final
one-band block whose completed action is independent of the residual bank,
even if its native gates temporarily move those wires. A second construction
below moves every bank direction in one band; a different Pauli-centralizer
argument excludes that particular transversal family as well. Other
bank-dependent final blocks and sources with more layers remain open.
These are representation boundaries, not an unrestricted native-depth theorem.

<a id="native-coupled-encoder"></a>

## 1. An actual one-band family

Let a have p address bits, with $`p\ge1`$. Let d be a principal bank
coordinate, Y the remaining principal bank coordinates, and r a register
of p arbitrary dirty tags. Consider the Boolean diagonal target

```math
Q=Z_dR,\qquad R|a,Y\rangle=(-1)^{f(a,Y)}|a,Y\rangle,
```

where f is independent of d. This includes lookup in its fixed binary
bank coordinates: d is the constant coefficient and the other bank
coordinates are the coefficients of the nonconstant address monomials.
All inputs, including r and external references, are arbitrary.

Define commuting Hermitian Pauli axes and literal native rotations by

```math
A=X_d,\qquad B_i=X_{a_i},\qquad C_i=X_dX_{a_i}Z_{r_i},\qquad
T(P)=e^{i\pi/8}e^{-i\pi P/8},
```

and put

```math
F=T(A)^p\prod_{i=1}^p T(B_i)T(C_i),\qquad G=H_dF.       \qquad(1)
```

Reduce $`T(A)^p`$ to a Clifford power and at most one odd T gate.
The remaining $`2p+(p\bmod2)`$ axes are independent: the unique
$`Z_{r_i}`$ label first forces every C_i coefficient to vanish,
then the distinct $`X_{a_i}`$ labels force the B_i coefficients to
vanish, and the optional A coefficient vanishes last. The dirty tags
are essential to this independence. Dropping them would leave dependent
axes and would not certify one physical band.

Here is a literal chronological word for G. Write
$`k=\lfloor p/2\rfloor\bmod4`$.

1. Apply H to d and every address bit.
2. For every i, apply CX(d,r_i), then CX(a_i,r_i).
3. Apply $`S_d^k`$, one T_d if p is odd, and one T to every a_i and r_i.
   All odd gates act on distinct physical wires in one band.
4. Reverse the CNOT word and apply H to every address bit.

The missing final H_d in this word is precisely the extra H_d in (1),
cancelled against the last basis-change H_d. The common forward and
actual-reverse ledger is

```math
T_G=2p+(p\bmod2),\qquad D_T(G)=1,\qquad
G_G\le6p+4,\qquad w_G=p,\qquad a_{\rm clean}=0.        \qquad(2)
```

The respective exact Clifford counts are $`6p+1+k`$ and
$`6p+1+((-k)\bmod4)`$, after reducing adjacent S powers. All
inverse occurrences reverse the emitted word and exchange T with
T-dagger. The family has zero implementation error. G can correlate the
dirty tags with the other inputs; their final return remains an obligation
of a completed source echo.

<a id="complete-hybrid-address-kernel"></a>

## 2. The complete conjugated operator

Let $`\omega=e^{i\pi/4}`$, $`\eta_i=(-1)^{r_i}`$ and
$`\alpha\in\{+1,-1\}`$. In the A=alpha eigenspace, (1) has the
exact block, apart from its one common native scalar,

```math
F_\alpha=e^{-i\pi p\alpha/8}
 \prod_{i:\alpha\eta_i=1}e^{-i\pi X_{a_i}/4}.          \qquad(3)
```

The common scalar is $`e^{3i\pi p/8}`$ and cancels with the actual
adjoint below. Q exchanges the two A eigenspaces. The products for
alpha and minus alpha act on complementary sets of address coordinates.
Consequently each matrix element has a unique intermediate address,
rather than a sum over possible intermediate addresses.

For input and output addresses a,a', define

```math
z_i=\begin{cases}a'_i,&r_i=d,\\a_i,&r_i\ne d,\end{cases}
\qquad
e(a,a',r,d)=\sum_i(-1)^{r_i\oplus d}(a_i\oplus a'_i),
```

where e is an ordinary integer. The only nonzero entries of
$`M=GQG^\dagger`$ are

```math
\begin{aligned}
&\langle1-d,a',Y,r|M|d,a,Y,r\rangle\\
&\qquad=2^{-p/2}\omega^{p(-1)^d}
             i^{e(a,a',r,d)}(-1)^{f(z,Y)}.             \qquad(4)
\end{aligned}
```

To verify the orientation, the right factor $`F_{(-1)^d}^\dagger`$
mixes the coordinates with r_i=d, while the left factor mixes their
complement. Their preserved coordinates force the displayed z. A changed
coordinate contributes +i on the right or -i on the left, giving e.
The H_d walls turn the exchanged A eigenspaces into the literal d flip.

Every column therefore has exactly $`2^p`$ entries of magnitude
$`2^{-p/2}`$, with the same affine support: d flips, a' is free, and
Y,r are unchanged. Formula (4) is a full matrix identity, including
coherent dirty tags and reference extensions. It is not an identity only
on one selected tag sector. The encoder generally does not commute with
R, and the residual sign in (4) depends on both input and output addresses.

<a id="bank-independent-final-boundary"></a>

## 3. A bank-independent final block cannot finish three-layer lookup

For the displayed encoder, an exact three-layer source would require
one-band native blocks B,C satisfying

```math
E=G^\dagger BC,\qquad K=CZ_dC^\dagger,\qquad
BKB^\dagger=GQG^\dagger.                              \qquad(5)
```

K is a literal Hermitian Clifford. An arbitrary Hermitian Clifford is
not assumed to have the required first-band realization C. The theorem
below excludes (5) even when K is allowed to be any Hermitian Clifford.

**Theorem.** Let $`p\ge3`$ and let the independent coefficient bank
contain one bit $`y_S`$ for each nonempty subset S of the address positions.
Write

```math
g_y(a)=\sum_{\varnothing\ne S\subseteq[p]}y_S\prod_{i\in S}a_i
 \pmod2,\qquad R_y=(-1)^{g_y(a)},\qquad Q_y=Z_dR_y.
```

Let G be any native circuit of at most one T layer whose **completed
operator** factors as $`G_{\rm work}\otimes I_y`$, where work includes
address, d and arbitrary additional inputs. Then
$`M=\sum_y|y\rangle\langle y|\otimes G_{\rm work}Q_yG_{\rm work}^\dagger`$
cannot equal $`BKB^\dagger`$ for any one-band native B and any
Hermitian Clifford K on the complete space. No restriction is imposed on
B's address, bank or work movement. G's native boundaries need not preserve
the factorization; the hypothesis concerns its full operator.

### 3.1. Bank fibers force a Pauli centralizer

Call a binary Pauli label good for M when its actual conjugate by M is
Clifford. Good labels form a linear subspace, because conjugation preserves
products and the Clifford operators form a group. Every bank Z is good,
since M commutes with it. Thus a good Pauli's bank-Z factor can be removed,
leaving, up to its irrelevant constant phase,
$`X_y^\delta\otimes P`$, where P is a work Pauli.

Set $`M_y=G_{\rm work}Q_yG_{\rm work}^\dagger`$. If
$`V=M(X_y^\delta\otimes P)M`$ is Clifford, its fixed bank translation
is delta and its work fibers are

```math
V_y=M_{y+\delta}PM_y.                                  \qquad(6)
```

For every coefficient shift s, $`V_sV_0^\dagger`$ is a scalar times a
Pauli. Strip V's fixed bank translation by multiplying by $`X_y^\delta`$.
Conjugating $`X_y^s`$ by this bank-preserving Clifford produces a Pauli
with bank translation s; its work fiber at bank zero is precisely the
displayed ratio, up to the ordering's constant phase.

Choose s to encode the point function $`g_s(a)=[a=x]`$ for a nonzero
address x. Such a coefficient vector exists because this function vanishes
at address zero. Put

```math
S_s=G_{\rm work}(I-2|x\rangle\langle x|_a\otimes I)G_{\rm work}^\dagger.
```

If D is the work Hilbert-space dimension and $`m=2^p`$, then
$`\mathrm{rank}(I-S_s)=D/m`$. All the M_y and S_s commute, and
$`M_{y+s}=M_yS_s`$. Therefore

```math
\begin{aligned}
V_s-V_0&=M_\delta(S_sPS_s-P)M_0,\\
\mathrm{rank}(V_s-V_0)&\le2D/m\lt D/2.                 \qquad(7)
\end{aligned}
```

A nonidentity scalar-Pauli operator differs from I in rank at least D/2.
For a nonscalar Pauli, an eigenvalue equal to one can occupy at most half
the space; a nontrivial scalar multiple of I has full difference rank.
Applying this to $`V_sV_0^\dagger`$ forces $`V_s=V_0`$ literally,
so P commutes with S_s. This covers every nonzero address projector; their
complement supplies address zero. Hence P commutes with every

```math
J_i=G_{\rm work}Z_{a_i}G_{\rm work}^\dagger.             \qquad(8)
```

The argument did not assume that the good Pauli has zero bank-X shift.
That is why arbitrary bank motion in a putative middle block is covered.

### 3.2. A one-band middle requires isotropic Pauli support

Suppose $`M=BKB^\dagger`$. Write the one-band B as a Clifford-conjugated
physical band followed by a Clifford, absorbing the latter into K. Let L
be the isotropic binary span of its independent odd axes. For every Pauli
P in $`L^\perp`$, the band commutes with P, and its conjugation of the
Pauli $`KPK^\dagger`$ is Clifford. Consequently

```math
L^\perp\subseteq\mathrm{Good}(M),\qquad
\mathrm{Good}(M)^\perp\subseteq L.                     \qquad(9)
```

A Pauli commuting with J_i must commute with every nonzero Pauli term in
its expansion: conjugation only changes each term's sign and cannot mix
its distinct labels. By (8), the span of all these support labels therefore
lies in $`\mathrm{Good}(M)^\perp`$, and hence must be isotropic.

But J_i is a one-band Pauli image. The retained
[physical-axis expansion](TWO_LAYER_DIAGONAL_RIGIDITY.md#1-model-and-full-space-diagonal-theorem)
has support $`P_i+L_i`$, with all coefficients nonzero and each active
axis anticommutes with the Pauli core P_i. If J_i is non-Pauli, its support
contains both P_i and P_i plus an active axis; these two terms anticommute.
This contradicts (9). Therefore every J_i must be a Pauli.

The resulting commuting independent Pauli family can be mapped back to
$`Z_{a_1},\ldots,Z_{a_p}`$, with the literal signs, by a Clifford C_0
on work alone. Put $`A=C_0G`$ and $`B'=C_0B`$. Then A commutes
with every address Z and with R. Its image
$`K_A=AZ_dA^\dagger`$ is Clifford, and (5) would imply

```math
R=B'K B'^\dagger K_A.                                 \qquad(10)
```

This is an actual two-band diagonal circuit. The
[two-layer diagonal rigidity theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md)
forces R to be diagonal C3, whose Boolean exponent has degree at most
three. Here $`g_y(a)`$ contains
$`y_{[p]}a_1\cdots a_p`$ and has degree p+1, a contradiction.
This completes the proof without assuming a first-band realization of K.

### 3.3. Scope and the displayed family

The encoder (1) is independent of the coefficient bank, so the theorem
closes its three-layer completion for p at least three. Its address image
makes the nonisotropic support particularly explicit:

```math
GZ_{a_i}G^\dagger
 =\frac{Z_{a_i}-Y_{a_i}-Z_dY_{a_i}Z_{r_i}
                    -Z_dZ_{a_i}Z_{r_i}}2.            \qquad(11)
```

The terms $`Z_{a_i}`$ and $`Y_{a_i}`$ already anticommute. Adding
an arbitrary fixed Clifford on the left of G, including a bank Hadamard,
does not repair the failure: it can be absorbed into B in (5). Additional
identity helpers also cannot repair it; all dimension and Pauli arguments
above include them.

The [three-layer exactness theorem](MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness)
extends the exclusion to complete error strictly below tau_*. This is a
zero-clean, full-operator statement. It does not cover a promised input
sector, a bank-dependent final block, or a source with more T layers.
In particular, it is not an unrestricted three-layer lower bound. The
source in the linked two-address construction is not retired by it.

<a id="diagonal-middle-boundary"></a>

### 3.4. A separate edge-phase screen also applies at two addresses

For the particular G in (1), a computationally diagonal one-band B in
(5) forces $`\deg(f)\le2`$, even for a general f rather than the full
coefficient bank. This retains a useful smaller-case distinction.

On the always nonzero edge $`a'=a`$, (4) has phase
$`\omega^{p(-1)^d}(-1)^{f(a,Y)}`$. A Clifford matrix has quadratic
phase on its affine support: its Choi state's diagonal Pauli stabilizers
fix the support, and its translation stabilizers give linear phase
differences. Restricting to this fixed affine edge keeps phase degree at
most two modulo eight. A diagonal one-band B is C3, with weighted phase
$`\sum l_ju_j+2\sum q_{jk}u_ju_k+4\sum t_{jkl}u_ju_ku_l\pmod8`$.
Its difference on flipping only d is quadratic modulo eight; the possible
remaining cubic coefficient is divisible by eight. Comparing edge phases
and subtracting $`p(-1)^d=p-2pd`$ would therefore make $`4f`$
quadratic modulo eight, forcing Boolean degree at most two. Thus this
particular diagonal-middle repair already fails for two-address lookup.

<a id="full-bank-x-encoder"></a>
<a id="native-full-bank-encoder"></a>

## 4. A full-bank X-motion encoder and its separate boundary

This candidate moves every bank direction, so Section 3's completed
bank-independence hypothesis does not apply. Use m original bank bits
$`Y_j`$, indexed by $`j\in\mathbb F_2^p`$, and one additional
arbitrary probe d. The target is $`Q=Z_d\Omega_Y`$, where
$`\Omega_Y|a,Y\rangle=(-1)^{Y_a}|a,Y\rangle`$.

This probe interface is paid and sufficient. Let L be the m-CNOT word
$`\prod_j\mathrm{CX}(d,Y_j)`$. It obeys $`LQL=\Omega_Y`$.
If an emitted E transported Z_d to Q, then
$`LE\,\mathrm{SWAP}(d,Y_0)`$ would transport Z_Y0 to Omega_Y
with only m+3 additional Cliffords. Its echo would return d on every
input. Bank Hadamards then give the indicator. No clean probe is used.

### 4.1. Native family and complete rewritten operator

Define the commuting CNOT word

```math
W=\prod_j\mathrm{CX}(Y_j,d)
       \prod_{i:j_i=1}\mathrm{CX}(Y_j,a_i),\qquad
T_X=\prod_j H_{Y_j}T_{Y_j}H_{Y_j},\qquad
G_{\rm bank}=H_dWT_XW^\dagger.                       \qquad(12)
```

W has m+pm/2 CNOTs and is an involution. Its column action is
$`(d,a,Y)\mapsto(d+\pi(Y),a+HY,Y)`$, where
$`\pi(Y)=\sum_jY_j\pmod2`$ and
$`HY=\bigoplus_{j:Y_j=1}j`$. The axes of G_bank's band are
$`X_dX_a^jX_{Y_j}`$, independent by their distinct bank X labels.
Emit W, bank H gates, one T on every bank wire, bank H gates,
the actual reverse of W, and H_d. Forward and actual inverse have

```math
T_{G_{\rm bank}}=m,\qquad D_T(G_{\rm bank})=1,\qquad
G_{G_{\rm bank}}=pm+4m+1,\qquad w_{G_{\rm bank}}=1,
\qquad a_{\rm clean}=0.                              \qquad(13)
```

The one work bit counted here is d; there are no tags. The actual
adjoint retains the emitted T phases and reverses every gate.
The fixed Clifford rewrite is

```math
WQW=Q_h=(-1)^h,\qquad
h=d+\pi(Y)+Y_{a+HY},\qquad
M'=WH_d(G_{\rm bank}QG_{\rm bank}^\dagger)H_dW
   =T_XQ_hT_X^\dagger.                              \qquad(14)
```

There is an explicit complete operator expansion. Let b_u be the
bank mask with $`(b_u)_j=1+u\cdot j`$, all binary. Fourier expansion
of the selected bank sign gives

```math
M'=\frac{Z_d}{m}\sum_{u,j\in\mathbb F_2^p}
 (-1)^{u\cdot j}Z_a^u
 \prod_{k:(b_u+e_j)_k=1}\frac{Z_{Y_k}-Y_{Y_k}}{\sqrt2}.
                                                               \qquad(15)
```

This is a full-space identity, including every bank, probe and reference
input. It follows by expanding Q_h as
$`m^{-1}Z_d\sum_{u,j}(-1)^{u\cdot j}Z_a^uZ_Y^{b_u+e_j}`$
and conjugating each bank Z by T_X. The common native scalar cancels
with its actual inverse.

<a id="full-bank-encoder-boundary"></a>
<a id="full-bank-x-motion-boundary"></a>

### 4.2. A Clifford centralizer lemma

Put $`R_j=(Z_{Y_j}-Y_{Y_j})/\sqrt2`$, a Hermitian single-qubit
Clifford. If a full Clifford C commutes with every R_j, then

```math
C=\left(\prod_j R_j^{e_j}\right)\otimes V,             \qquad(16)
```

where V is a Clifford on all nonbank wires, including its common phase.
No register-preservation premise is needed. Indeed the Pauli expansion of
$`CR_jC^\dagger=R_j`$ forces exactly one signed case on each wire:
$`(CZ_jC^\dagger,CY_jC^\dagger)=(Z_j,Y_j)`$ or
$`(-Y_j,-Z_j)`$. The latter is conjugation by R_j. Multiplying by
the corresponding product of R_j makes C fix the entire bank Pauli
algebra pointwise, hence factor as bank identity times V. V is Clifford
because the remaining Pauli conjugation action is Clifford. This also
excludes a bank permutation or a controlled R_j; scalar phases remain
in V.

Consequently, if D is Boolean diagonal and $`T_XDT_X^\dagger`$
is Clifford, then D is a bank-Z Pauli times a Boolean diagonal
Clifford on the remaining wires. Its bank coefficients are linear and
independent of address/work, and its nonbank exponent has degree at
most two. This follows from (16), since D commutes with every bank Z.

### 4.3. Good Paulis cannot supply a one-band middle

Use the good-label definition and necessary isotropy condition (9).
M' commutes with every R_j. If P is good and has Z or Y component
on bank wire j, then $`PR_jP=\pm R_j^+`$, where
$`R_j^+=(Z_j+Y_j)/\sqrt2=iX_jR_j`$. Thus the Clifford
$`C=M'PM'`$ satisfies

```math
CR_jC^\dagger=\pm i(M'X_jM')R_j.                      \qquad(17)
```

It follows that X_j would also be good. But X_j commutes with T_X,
so goodness requires $`T_X(-1)^{\Delta_jh}T_X^\dagger`$ to be
Clifford. At Y=0, for every j including j=0,
$`\Delta_jh=1+[a=0]`$, of Boolean degree p. Section 4.2 excludes
this for p at least three. Therefore every good Pauli has zero
bank-Z label: its bank factor is X_Y^v.

Write the address translation of such a Pauli as u, its probe
translation as epsilon, and set $`t=u+Hv`$, $`k=a+HY`$.
The Pauli commutes with T_X. Its goodness therefore requires the
following derivative to satisfy Section 4.2's bank-linearity condition:

```math
\Delta h=\epsilon+\pi(v)+Y_k+Y_{k+t}+v_{k+t}.          \qquad(18)
```

Compare Y=0 with Y=e_0. Their difference is $`[a=0]+[a=t]`$,
which is nonconstant when t is nonzero and p is at least two. Hence
t=0. Equation (18) becomes $`\epsilon+\pi(v)+v_{a+HY}`$.
Comparing Y=e_j with Y=0 now shows that every derivative
$`v_{a+j}+v_a`$ is independent of a. Thus v, viewed as a bank
truth table on row labels, must be affine. Its possible masks lie
in the p+1-dimensional span of the constant mask and the coordinate
masks $`(j_i)_j`$.

These necessary containments already contradict (9). Good(M') has
no bank-Z projection and only affine bank-X projections, so its
orthogonal complement contains both

```math
X_{Y_0},\qquad
Z_{Y_0}Z_{Y_{e_1}}Z_{Y_{e_2}}Z_{Y_{e_1+e_2}}.         \qquad(19)
```

Every affine truth table has even parity on this four-point plane.
The two displayed Paulis anticommute, so the complement is not
isotropic. No arbitrary one-band B and Hermitian Clifford K can
therefore satisfy $`M'=BKB^\dagger`$ for p at least three.
This argument does not assume that K has a first-band Pauli realization.

### 4.4. Exact scope and reopening condition

The conclusion closes this particular full-bank X-motion final encoder
at three source layers. It permits arbitrary bank/address motion in a
putative middle, arbitrary identity helpers, and a fixed output Clifford
on G_bank: the last is absorbed into B. The retained three-layer
rounding theorem extends it to complete error strictly below tau_*.
It does not exclude a different bank-dependent final encoder, a changed
nonmonomial input Clifford wall, or a source with four or more layers.

Decorating the axes with private Z tags does not change this result:
conjugation by $`\prod_j\mathrm{CZ}(Y_j,r_j)`$ adds those tags
and commutes with Q. It is only a Clifford change of coordinates on
the same candidate. Reopening needs a different complete operator,
rather than an assumed middle or an unpriced tag initialization.

## 5. Research consequence and evidence

The family pays a noncommuting final encoder at every address length and
identifies its complete interference pattern. Its proposed three-layer
completion is now paused at p at least three. The same obstruction covers
all completed bank-independent one-band final blocks, not just its
triangle axes. Merely choosing another non-diagonal middle does not reopen
that class.

The full-bank X-motion family also supplies an actual native encoder and
complete operator rewrite, but its separate centralizer obstruction stops
three-layer completion. Bank motion by itself does not certify the missing
source. The next constructive operation must change the final encoder's
interaction with the bank beyond the two stopped families, or move to
more source layers and pay a complete exact or approximate mechanism.
The normalization, lookup sectors, arbitrary work/reference inputs,
precision and actual inverses remain part of that obligation. A complete
source would still need the simultaneous resource ledger and the retained
fusion and frame transfer.

No improved complete-frame bound follows here. R96 remains the general
upper frontier. The constant-depth approximate scalar zero reflection is
unaffected, and no conjunction depth witness is being reinstated.

The [exact native fixtures](../../tests/test_coupled_final_encoder.py)
check (4), every dirty-tag input and actual inverses for a complete
two-address coefficient bank and a cubic three-address sign. They also
check physical-axis independence and the native ledgers. Bounded full-bank
fixtures also check the self-addressed charge through p=3, literal inverses
and expansion columns at p=1,2, and ledgers through p=6. The all-size kernel
and both scoped final-block obstructions are analytic proofs;
these small fixtures do not establish the obstruction or a general source.
The proofs have internal audits, not external mathematical review.

<a id="syndrome-star-coupling"></a>

## 6. A syndrome-star coupling still needs complete readout

There is a concrete higher-layer direction using the paid approximate scalar
zero reflection. Use the raw bank Y indexed by all binary addresses j, and
put

```math
s=a+HY,\qquad HY=\bigoplus_{j:Y_j=1}j,\qquad
R_S|a,Y\rangle=(-1)^{[s=0]}|a,Y\rangle.                \qquad(20)
```

The syndrome part of the Clifford W in Section 4 computes s into the
existing address with pm/2 CNOTs. Apply the
[complete scalar zero reflection](EXPANDER_ZERO_REFLECTION.md) there and
reverse those CNOTs. At requested full-operator error epsilon, this gives
an actual zero-clean approximation to R_S with
$`T,w=O_\epsilon(p)`$, $`G=pm+O_\epsilon(p)`$ and
$`D_T=O_\epsilon(1)`$. Every address, bank, implementation helper and
reference input is covered. The actual reversed word supplies its inverse;
no approximate reflection is assumed to be exactly involutory.

Consider the **formal Hermitian sums**

```math
H_X=\sum_{j\ne0}X_{Y_j},\qquad
D=\frac{H_X-R_SH_XR_S}{2}.                             \qquad(21)
```

Flipping Y_j changes s to s+j. Consequently R_S changes the sign of this
edge precisely when s is zero or j, and the complete basis action is

```math
D|a,Y\rangle=
\begin{cases}
\displaystyle\sum_{j\ne0}|a,Y+e_j\rangle,&s=0,\\
|a,Y+e_s\rangle,&s\ne0.
\end{cases}                                           \qquad(22)
```

This has a simple full-space decomposition. For each fixed address a,
every syndrome-zero bank column Y_c is the center of one star with leaves
$`Y_c+e_j`$, j nonzero. Each leaf has syndrome j and has exactly that
center as its neighbor in (22), so the stars are disjoint and cover the
whole bank. There are $`2^{m-p}`$ stars, each with m vertices. Define

```math
|b_c\rangle=\frac1{\sqrt{m-1}}\sum_{j\ne0}|Y_c+e_j\rangle.
```

On this star, (21) is exactly

```math
D_c=\sqrt{m-1}\bigl(|b_c\rangle\langle Y_c|
                         +|Y_c\rangle\langle b_c|\bigr).\qquad(23)
```

Thus $`\|D\|=\sqrt{m-1}`$. Each star also has m-2 orthogonal leaf
directions annihilated by D. The normalized coupling's matrix element
from any one leaf to its center is $`1/\sqrt{m-1}`$. The apparently
addressed X action in the second line of (22) is not a complete unitary
lookup: all m-1 leaves map to the same center. Its first line and its dark
directions cannot be discarded on arbitrary bank inputs.

Equations (21)--(23) are operator identities, not native implementations of
a linear combination or of its exponential. Even the direct precision
estimate for this formal combination has a dimension factor. If an emitted
$`\widetilde R_S`$ has error at most epsilon, use its actual adjoint in
$`\widetilde D=(H_X-\widetilde R_SH_X\widetilde R_S^\dagger)/2`$;
then $`\|\widetilde D-D\|\le\epsilon(m-1)`$. This is a sufficient
telescoping estimate, not a necessary precision bound on every changed
mechanism.

The first unpaid operation is a native unitary or accepted-block readout
with the required constant normalization on **every** star direction and
the correct lookup sectors. The scalar reflection component and the formal
star action do not supply $`K^\dagger K\approx\Pi_{\rm eq}/2`$ or a
complete indicator. A proposed pulse word must prove that operator identity
and pay its rotations, precision, actual inverses and work return together.
No lower bound on arbitrary finite pulse words or higher-layer sources is
claimed by this decomposition.
