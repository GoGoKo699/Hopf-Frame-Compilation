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
The later sections pay two syndrome readouts and a shifted-syndrome
four-band commutator, deriving their complete actions and separate failures.

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

```math
\sum l_ju_j+2\sum q_{jk}u_ju_k+4\sum t_{jkl}u_ju_ku_l\pmod8.
```

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
\begin{aligned}
WQW&=Q_h=(-1)^h,\qquad h=d+\pi(Y)+Y_{a+HY},\\
M'&=WH_d(G_{\rm bank}QG_{\rm bank}^\dagger)H_dW\\
  &=T_XQ_hT_X^\dagger.                              \qquad(14)
\end{aligned}
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

The scalar reflection and the formal star action do not yet supply a
native readout with constant normalization on **every** star direction.
The following paid pulse candidate tests that missing operation directly.

<a id="syndrome-star-pulse-boundaries"></a>

### 6.1 The original pulse algebra has a blind translation

For a bank mask v with Hv=0, the translation $`X^v`$ commutes with R_S
and with every product of bank-X rotations, at arbitrary angles. This
remains true after arbitrary address-only, selector-only or helper-only
interlayers, including selector-controlled R_S calls. Thus every ideal
word in this grammar, its accepted block and its accepted Gram commute
with $`X^v`$. The claim concerns completed operators, not the native
boundaries of an approximate scalar-zero compiler.

The complete lookup sign $`\Omega_{a,Y}=(-1)^{Y_a}`$ does not have
this symmetry: conjugation by $`X^v`$ multiplies it by $`(-1)^{v_a}`$.
Any nonzero kernel mask therefore gives direct-phase error at least one.
For the canonical two-bank target, a translation on one bank gives

```math
\left\|G-\frac{I+\Omega_{a,C}\Omega_{a,D}}4\right\|
 \ge\frac14                                                     \qquad(24)
```

whenever G has the stated symmetry. Indeed the target and its translated
copy differ by norm one half, whereas G is unchanged. The address may be
coherent; restricting to one address with v_a=1 establishes the norm gap.
Adding a fixed $`Z_{Y_0}`$ baseline still leaves, at p at least two,
the kernel mask supported on three distinct rows i,j,i+j. An actual
implementation within full-word error eta of this ideal grammar retains
direct gap at least $`1-\eta`$ and Gram gap at least $`1/4-2\eta`$.
No exact commutation is asserted for its approximate native subwords.

### 6.2 A wide-T changed word pays control but remains thin

Use the fixed Clifford wall $`L:(C,D)\mapsto(C,B=C+D)`$, with
m CNOTs. It turns the canonical target into
$`(I+\Omega_{a,B})/4`$; C remains an arbitrary spectator. On B put

```math
\begin{aligned}
P&=[a+HB=0],\qquad R=I-2P,\\
T_X&=H_B T_B H_B,\qquad F=T_X T_B T_X,\\
U&=R F^\dagger R F,\qquad K=\frac{I-U}{2}.
\end{aligned}                                                   \qquad(25)
```

Here T_B applies one literal T on each bank wire. F is an emitted
three-band word with 3m T gates and 4m Hadamards. Its middle T_B prevents
the cancellation that would occur if a bank-Z dressing merely conjugated
R. This is a changed interference mechanism: the exact p=3 fixture below
shows that its accepted Gram no longer commutes with the kernel
translation supported on rows 1,2,3.

There is no unpaid controlled F in (25). With one clean acceptance bit q,
put

```math
C_R=|0\rangle\langle0|_q\otimes I+
|1\rangle\langle1|_q\otimes R.
```

Emit, chronologically,

```math
L,\ H_q,\ F,\ C_R,\ F^\dagger,\ C_R^\dagger,
\ Z_q,\ H_q,\ L.                                                \qquad(26)
```

The middle word is controlled U: its q=0 branch is the actual
$`F^\dagger F=I`$, and its q=1 branch is (25). The zero-to-zero
accepted block of (26) is $`L K L`$. To implement C_R, compute
s=a+HB into the existing address, flip q, apply the paid scalar zero
reflection on (s,q), undo that flip, and undo the syndrome CNOTs. This
tests q=1 and s=0 without a second clean bit. All scalar-zero helpers are
arbitrary dirty inputs. This construction also specifies the actual
adjoint in (26); an approximate emitted reflection is never replaced by
an assumed involution.

Write T_0,G_0,w_0,d_0 for the scalar-zero compiler's costs on p+1 bits
at full error epsilon. One common safe source/inverse ledger is

```math
\begin{aligned}
T&=6m+2T_0,& w_{\rm dirty}&=w_0,\\
G&\le2pm+10m+20+2G_0,& D_T&\le6+2d_0.
\end{aligned}                                                   \qquad(27)
```

The 20 term permits X and Z to be expanded using H and positive S.
There is one clean acceptance bit, no initialized bank, and no hidden
clean history. Thus at fixed epsilon this candidate has linear T-count,
dirty width $`O_\epsilon(p)`$, $`O(pm)+O_\epsilon(p)`$ Cliffords,
and constant depth. These are costs of the failed accepted source, not a
complete returned indicator. The full-operator estimates include arbitrary
bank, address, dirty work and reference inputs.

Nevertheless, the exact complete Gram is

```math
\begin{aligned}
K^\dagger K&=\frac{I-\mathrm{Re}\,U}{2}
            =P+P'-PP'-P'P,\qquad P'=F^\dagger PF,\\
\mathrm{ran}(K)&\subseteq\mathrm{ran}(P)+\mathrm{ran}(P'),\qquad
\mathrm{rank}(K^\dagger K)\le\frac{2D}{m}.
\end{aligned}                                                   \qquad(28)
```

D denotes the entire unaccepted input dimension, including the arbitrary
C bank and any identity helpers. The syndrome map has rank p, so P has
rank D/m. The target $`(I+\Omega_{a,B})/4`$ has eigenvalue one half
on a subspace of dimension D/2. For m greater than four, that subspace
has a nonzero intersection with the kernel of (28). Consequently

```math
\left\|K^\dagger K-\frac{I+\Omega_{a,B}}4\right\|
 \ge\frac12\qquad(p\ge3).                                    \qquad(29)
```

This failure is independent of the old translation symmetry. Replacing
the two controlled reflections by an emitted word of error at most
epsilon and its actual adjoint changes the complete source and K by at
most 2 epsilon, and its Gram by at most 4 epsilon. The actual native
candidate therefore retains gap at least $`1/2-4\epsilon`$.

The same rank argument pauses a fixed-call extension with zero accepted
baseline. If a complete word V contains k such thin reflections and
deleting them leaves I, telescoping gives
$`\mathrm{rank}(V-I)\le kD/m`$, even with arbitrary intervening
unitaries. Acceptance $`(I-V)/2`$ cannot have error below one half at
k less than m/2. This is a limitation of that baseline and call grammar;
it is neither a native T-depth bound nor an exclusion of a broad-baseline
source. Combining many calls into fewer physical bands would need an
actual native schedule of its own.

The [exact fixtures](../../tests/test_syndrome_star_readout.py) check all
columns and the actual inverse of the p=1 native two-bank acceptance,
including its normalization, and two p=3 Gram entries witnessing the
changed symmetry. The all-size support and error statements are analytic;
the finite checks do not prove them. This internal assessment supplies no
external mathematical review. The first unpaid operation is now a broad
accepted-sector mechanism with the correct lookup action, not control of
(25), its native synthesis, or a repair of the stopped fixed-call word.

<a id="broad-baseline-syndrome-source"></a>

### 6.3 A normalized broad baseline does not finish the wide-T trial

The next actual source starts with a broad fixed-row accepted sector,
so its failure is not inferred from zero acceptance. Retain L, F, R and
$`C_R`$ from (25)--(26), and define the literal signal rotation

```math
T_Y=S_qH_qT_qH_qS_q^\dagger
   =e^{i\pi/8}e^{-i\pi Y_q/8}.
```

Emit, chronologically,

```math
L,\ T_Y,\ F,\ C_R,\ F^\dagger,\ C_R^\dagger,
\ \mathrm{CZ}(q,B_0),\ T_Y^\dagger,\ H_q,\ L.           \qquad(30)
```

All inverses are actual native reversals. The signal is the only clean
input. C, B and the scalar-zero compiler's helpers are arbitrary, and
intermediate motion of the signal is allowed in the approximate compiler.
In the ideal word, the middle controlled commutator has signal branches
I and $`U=RF^\dagger RF`$. Put $`c=\cos(\pi/8)`$ and
$`s=\sin(\pi/8)`$. The two signal boundaries give

```math
T_Y|0\rangle=e^{i\pi/8}(c|0\rangle+s|1\rangle),\qquad
\langle0|H_qT_Y^\dagger=e^{-i\pi/8}(s\langle0|+c\langle1|).
```

Their literal scalars cancel and $`cs=1/(2\sqrt2)`$. Hence the
complete zero-to-zero block and Gram, before the two L walls, are

```math
\begin{aligned}
K_{\rm b}&=\frac{I+Z_{B_0}U}{2\sqrt2},\\
G_{\rm b}=K_{\rm b}^\dagger K_{\rm b}
 &=\frac{I+\mathrm{Re}(Z_{B_0}U)}4,
 &0\le G_{\rm b}\le I/2.                              \qquad(31)
\end{aligned}
```

The canonical block is $`LK_{\rm b}L`$. Deleting both controlled
syndrome calls makes the actual F and its adjoint cancel and leaves

```math
K_0=\frac{I+Z_{B_0}}{2\sqrt2},\qquad
K_0^\dagger K_0=\frac{I+Z_{B_0}}4.                    \qquad(32)
```

Thus the trial has correct half normalization and a baseline of rank
half the complete input dimension. It does not initialize a row bank
or assume a clean syndrome.

Each $`T_Y`$ or its actual adjoint uses one T-type gate and six
Cliffords when S-dagger is expanded as three positive S gates. Its
signal envelope in (30), including CZ and H, costs two T gates, two
T layers and sixteen Cliffords. The controlled syndrome word costs
$`T_0`$ T gates, at most $`pm+G_0+8`$ Cliffords, depth $`d_0`$
and $`w_0`$ dirty helpers. Consequently one common source/inverse
ledger, at scalar-zero error epsilon, is

```math
\begin{aligned}
T_{\rm b}&\le6m+2T_0+2,&w_{\rm b}&\le w_0,\\
G_{\rm b,native}&\le2pm+10m+32+2G_0,
&D_T&\le8+2d_0.                                      \qquad(33)
\end{aligned}
```

There is exactly one clean acceptance bit. At fixed epsilon these
are linear T-count, $`O_\epsilon(p)`$ dirty width,
$`O(pm)+O_\epsilon(p)`$ Cliffords and constant T-depth. Two
approximate controlled-reflection calls change the full source by at
most $`2\epsilon`$ and its accepted Gram by at most
$`4\epsilon`$. These estimates use the emitted reflection and its
actual adjoint; they include dirty return and arbitrary references,
without asserting that the emitted signal branches remain exact.

The paid source nevertheless fails on the correct lookup sectors.
The syndrome ignores $`B_0`$, and F is a tensor product of identical
single-wire words. Its $`B_0`$ factor therefore cancels inside U:

```math
U=I_{B_0}\otimes U_{\rm rest}.                       \qquad(34)
```

Fix any nonzero address a, and let D denote the dimension of the
remaining arbitrary inputs, including the whole bank B. On the
$`B_0=1`$ sector, each of P and $`F^\dagger PF`$ has rank
$`D/(2m)`$. The two-reflection difference therefore has

```math
\mathrm{rank}\bigl((I-U)|_{B_0=1}\bigr)\le D/m. \qquad(35)
```

The target-positive subspace $`B_0=1,B_a=0`$ has dimension D/4.
When $`m\gt4`$, it intersects the kernel in (35) nontrivially.
On a unit vector in this intersection, (31) vanishes while
$`(I+\Omega_{a,B})/4`$ acts as one half. Since both positive
operators are bounded above by I/2, the complete ideal distance is
exactly

```math
\left\|G_{\rm b}-\frac{I+\Omega_{a,B}}4\right\|
 =\frac12\qquad(p\ge3).                              \qquad(36)
```

The L walls preserve this distance to the canonical two-bank target.
The emitted source retains gap at least $`1/2-4\epsilon`$. The
rank argument grants arbitrary identity helpers and a coherent address;
restriction to the chosen address fiber already witnesses operator norm.

There is a useful **baseline-relative** statement, without the special
$`B_0`$ factorization. Let the unaccepted input have dimension D,
and insert k controlled syndrome reflections into any specified native
baseline. Deleting those calls defines the baseline unitary and its
accepted block K_0. Each physical controlled reflection has difference
rank D/m, on the full signal-plus-input space of dimension 2D.
Telescoping and compression give

```math
\mathrm{rank}(K-K_0)\le kD/m.                   \qquad(37)
```

If a subspace S lies in $`\ker K_0`$ and the desired half-projector
acts there as I/2, then $`\dim S\gt kD/m`$ forces Gram error at
least one half. For the fixed-row baseline (32), across all addresses,
S can be chosen with $`a\ne0,B_0=1,B_a=0`$ and has dimension
$`D(m-1)/(4m)`$. Thus fewer than $`(m-1)/4`$ such insertions
cannot repair that baseline. The sharper p=3 result (36) uses (34);
it is not obtained by silently halving the general physical rank.
A baseline with no such kernel subspace is outside this stated test.

This pauses the particular broad-baseline trial (30). It does not turn
thin-call count into native depth: a genuinely joint schedule for a
growing number of calls, or a different broad baseline already carrying
the required lookup information, remains a constructive possibility.
The two-band scalar rotation envelope is paid; the missing operation is
that complete lookup-dependent joint action. Neither this source failure
nor the earlier thin-acceptance result changes R96 or excludes unrestricted
higher-band sources.

The [exact broad-source fixtures](../../tests/test_broad_syndrome_readout.py)
check all physical columns and the actual inverse of the p=1 canonical
source, its full accepted block, the half normalization and the two-T
baseline. This p=1 controlled reflection is an exact Clifford
specialization, not an execution of the general approximate compiler.
A p=2 matrix checks the $`B_0`$ spectator, and an exact
p=3 dimension count checks the two subspaces used in (35). These are
finite internal checks. The all-size operator, resource and rank claims
above are analytic and internally audited, not externally reviewed.

<a id="shifted-syndrome-wide-commutator"></a>

## 7. A four-band broad word retains the wrong Pauli coefficient

A different native trial uses no thin reflection or accepted compression.
Let $`p\ge1`$, $`m=2^p`$, and let a and B be the arbitrary address
and m-bit bank. Fix the nonzero row label $`v=e_0`$. Define

```math
\begin{aligned}
P&=Z_{B_0},&
P_j&=X_a^{j+v}X_{B_j},& Q_j&=Z_a^jZ_{B_j},\\
A_v&=\prod_{j\in\mathbb F_2^p}e^{-i\pi P_j/8},&
B_0&=\prod_{j\in\mathbb F_2^p}e^{-i\pi Q_j/8},&
E&=A_v^\dagger B_0^\dagger A_vB_0.                    \qquad(38)
\end{aligned}
```

Here the subscript on the operator $`B_0`$ denotes the unshifted
syndrome word; $`B_j`$ inside a Pauli denotes a bank wire. The
shift matters: the unshifted syndrome never couples the physical
wire $`B_0`$ to the address. Alternating shifted and unshifted
words removes that particular spectator argument. Both operators
in (38) act on the whole Hilbert space, without a clean sector.

### Actual native family and complete action

The Clifford word $`W_v`$ consists of CNOTs from bank wire j to
address bit i whenever $`(j+v)_i=1`$. Its basis action is

```math
W_v|a,B\rangle
 =\left|a+\sum_j(j+v)B_j,B\right\rangle.
```

It is its own actual inverse and costs pm/2 CNOTs. Emit W_v, H on
every bank wire, one transversal T band, the same bank Hadamards,
and the reversed W_v. The literal operator is
$`\widetilde A_v=e^{im\pi/8}A_v`$ and costs m T gates and
$`pm+2m`$ Cliffords. For $`\widetilde B_0`$, instead compute each
address parity $`a\cdot j`$ into bank wire j, apply one transversal
T band, and undo the parity CNOTs. This has literal operator
$`\widetilde B_0=e^{im\pi/8}B_0`$, m T gates and pm Cliffords.
Every band has m distinct physical T targets.

The chronological source is

```math
\widetilde B_0,\quad\widetilde A_v,\quad\widetilde B_0^\dagger,\quad
\widetilde A_v^\dagger.
```

Its four literal scalars cancel, giving
exactly E in (38). The source and its actual reversed native word have

```math
T_E=4m,\qquad D_T(E)\le4,\qquad G_E=4pm+4m,
\qquad w_E=a_{\rm clean}=0.                           \qquad(39)
```

All inputs and external references remain arbitrary. The Pauli echo
$`EPE^\dagger`$ costs $`8m`$ T gates, at most eight T layers and
$`8pm+8m+2`$ Cliffords. These are ledgers for an actual candidate,
not a claim that its echo implements lookup or returns the desired
address and bank state.

Put $`R=-iP_0P`$, a Hermitian Pauli. Only P_0 anticommutes with P,
so the first conjugation is exactly

```math
A_vPA_v^\dagger=\frac{P+R}{\sqrt2}.
```

All Q_j commute with P. The axes anticommute with R precisely for

```math
\mathcal J=\{j:v\cdot j+\mathbf1_{j=0}=1\},
\qquad r=|\mathcal J|=m/2+1.
```

The m/2 labels satisfying $`v\cdot j=1`$ are nonzero, and j=0
adds one more. Consequently the complete echo has the operator formula

```math
\begin{aligned}
R'&=B_0^\dagger R B_0
   =R\prod_{j\in\mathcal J}\frac{I-iQ_j}{\sqrt2},\\
EPE^\dagger&=A_v^\dagger\frac{P+R'}{\sqrt2}A_v.
                                                               \qquad(40)
\end{aligned}
```

This is a whole-space identity, not an expectation on a prepared bank.
The Q_j are independent Pauli labels, because each has its own bank
Z label. Thus only the empty product contributes the original R in
R', with normalized Pauli coefficient $`2^{-r/2}`$.

### A complete norm witness and its scope

For operators on the address and bank, write
$`\langle F,G\rangle=2^{-(p+m)}\mathrm{Tr}(F^\dagger G)`$.
The P component of R' is zero, since conjugation by B_0 preserves P
and $`\langle P,R\rangle=0`$. Unitary invariance and (40) therefore
give

```math
\begin{aligned}
\langle P,EPE^\dagger\rangle
 &=\frac12\langle P+R,P+R'\rangle\\
 &=\frac12+2^{-r/2-1}.                                  \qquad(41)
\end{aligned}
```

The required lookup sign is
$`\Omega_{a,B}=\sum_a|a\rangle\langle a|\otimes Z_{B_a}`$.
Only its address-zero block contributes to its P coefficient, so
$`\langle P,\Omega_{a,B}\rangle=1/m`$. Since P is unitary,
normalized trace is bounded by operator norm. It follows that

```math
\left\|EPE^\dagger-\Omega_{a,B}\right\|
 \ge\frac12+2^{-(m/2+1)/2-1}-\frac1m.                    \qquad(42)
```

The bound is one quarter at p=1 and strictly greater than one quarter
for every p at least two. It also holds after tensoring arbitrary
reference or identity-work systems. This exact native candidate fails
at a dimension-independent requested accuracy even though all four
bands, Clifford walls and actual inverses are paid. No approximation
or thin-reflection rank argument enters the proof.

The same coefficient and gap hold with independently chosen T or
T-dagger orientations on every axis of A_v and B_0, using the actual
adjoints in E. If the orientation of P_0 is sigma, the first image is
$`(P+\sigma R)/\sqrt2`$. For B_0 each retained R factor is still
$`\cos(\pm\pi/4)=1/\sqrt2`$, and both cross coefficients with P
are zero. The two occurrences of sigma multiply to one in (41).
Literal native scalars change with the orientations but still cancel
in the commutator. Thus sign tuning within these four bands does not
repair the mismatch; (39) is unchanged.

Pause this commutator as a lookup source. The retained half-weight
P component is the concrete mismatch: the middle B_0 preserves P,
while its independent quarter-phase axes leave a positive R component.
A changed middle action on P, changed angles, a different source word
or a different accepted Gram is outside (42). In particular this is
not an obstruction to the whole shifted-syndrome grammar, an
unrestricted four-layer theorem, or a growing-depth lower bound.
It does not reverse the constant-depth approximate scalar zero test
or change the complete-frame frontier.

The [exact commutator fixtures](../../tests/test_shifted_syndrome_commutator.py)
check every native source column and actual inverse at p=1 and p=2,
the Pauli coefficients $`3/4`$ and $`1/2+\sqrt2/8`$, and the
distinct-target bands and ledgers. A p=1 full echo checks the operator
orientation in (40), and a signed p=2 word checks the unchanged
coefficient under mixed orientations. Those finite checks challenge the
literal phase, orientation and coefficient calculation; (38)--(42)
are the analytic all-size argument. They are internal evidence, not
external mathematical review.
