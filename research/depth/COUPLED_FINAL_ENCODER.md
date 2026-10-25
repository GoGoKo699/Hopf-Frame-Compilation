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

## 7. A four-band broad word misses complete bank rows

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

For p at least two, the same complete operator has a stronger witness
that also survives arbitrary axis-angle changes. Choose a nonzero row
label k with $`v\cdot k=0`$. Such a label exists because this
hyperplane has dimension p minus one. Then k is absent from
$`\mathcal J`$. Both P and R commute with $`X_{B_k}`$, and so
does every Q_j appearing in (40), since j differs from k. All the
axes of A_v are X type, hence A_v also commutes with $`X_{B_k}`$.
Equation (40) gives

```math
[EPE^\dagger,X_{B_k}]=0.                               \qquad(43)
```

Take a normalized input with address k, bank wire k in its plus
X eigenstate and all other bank wires zero. The actual unitary echo
keeps the output in the plus eigenspace of $`X_{B_k}`$, regardless
of its address motion. The required lookup applies $`Z_{B_k}`$ on
that address and sends the input to the minus eigenspace. The actual
and desired outputs are normalized and orthogonal. Therefore

```math
\left\|EPE^\dagger-\Omega_{a,B}\right\|\ge\sqrt2
\qquad(p\ge2).                                        \qquad(44)
```

For the angle extension, replace the exponent on each P_j or Q_j
by any real angle, retaining the same axes and actual commutator.
The first image $`A_vPA_v^\dagger`$ still lies in the span of P
and R. The Q_k factor commutes with both and drops out of their
conjugation; all remaining factors commute with $`X_{B_k}`$.
Thus (43)--(44) still hold. This is an operator statement about that
fixed-axis word, not a claim that arbitrary rotations have the native
ledger (39). The native T/T-dagger family already has the gap (44).

There is also a broader midpoint test with the native A_v fixed.
Replace B_0 by any unitary B commuting with P, including arbitrary
address, bank and work motion consistent with that commutation, and
put $`E_B=A_v^\dagger B^\dagger A_vB`$. The Hermitian unitary
$`V=A_v^\dagger PA_v=(P-R)/\sqrt2`$ witnesses the failure.
Cyclic invariance of normalized trace and $`BP=PB`$ give

```math
\begin{aligned}
\langle V,E_BPE_B^\dagger\rangle
 &=\langle P,B^\dagger(A_vPA_v^\dagger)B\rangle
   =\frac1{\sqrt2},\\
\langle V,\Omega_{a,B}\rangle&=\frac1{m\sqrt2},\\
\left\|E_BPE_B^\dagger-\Omega_{a,B}\right\|
 &\ge\frac{1-1/m}{\sqrt2}.                              \qquad(45)
\end{aligned}
```

In the second line the R term vanishes because R is off-diagonal and
lookup is diagonal. The notation B in E_B denotes the replacement
unitary; the bank in $`\Omega_{a,B}`$ is unchanged. Arbitrary extra
work is included by identity extension and normalized trace on the
full space. In particular, a diagonal middle product with many tags
$`Z_a^{M_\ell j}Z_{B_j}`$ still commutes with P. Spanning the
address labels with those tags cannot repair this midpoint mismatch,
even if it removes every particular missing-row witness. No compiler
for that proposed product is needed to apply (45).
The numerical constant in (45) uses the native outer angles of A_v;
changing that outer operation is outside this stated midpoint bound.
For the original single-tag axes, (44) separately survives changes
to all angles, including the outer ones.

Pause this commutator as a lookup source. Its coefficient mismatch
at p=1 and unseen bank row at p at least two are both paid-source
failures. A changed middle action on P, changed axes, a different
source word or a different accepted Gram can leave these hypotheses;
changing only the axis angles does not. In particular this is not
an obstruction to the whole shifted-syndrome grammar, an
unrestricted four-layer theorem, or a growing-depth lower bound.
It does not reverse the constant-depth approximate scalar zero test
or change the complete-frame frontier.

The [exact commutator fixtures](../../tests/test_shifted_syndrome_commutator.py)
check every native source column and actual inverse at p=1 and p=2,
the Pauli coefficients $`3/4`$ and $`1/2+\sqrt2/8`$, and the
distinct-target bands and ledgers. A p=1 full echo checks the operator
orientation in (40), and a signed p=2 word checks the unchanged
coefficient under mixed orientations. A separate p=2 full native echo
checks (43) and the orthogonal witness with squared distance two.
Those finite checks challenge the literal phase, orientation, coefficient
and missing-row calculation; (38)--(45)
are the analytic all-size argument. They are internal evidence, not
external mathematical review.

<a id="common-probe-wide-source"></a>

## 8. Common-probe bands and a paired invariant state

The next source changes the first image as well as the middle action.
Let d be one arbitrary probe, let a be the p-bit address, and let B
be the complete bank of $`m=2^p`$ arbitrary bits. Put

```math
\begin{aligned}
P&=Z_d,& P_j&=X_dX_a^jX_{B_j},& Q_j&=Y_dZ_a^jZ_{B_j},\\
A&=\prod_j e^{-i\pi P_j/8},&
B_{\rm mid}&=\prod_j e^{-i\pi Q_j/8},\\
C|B\rangle&=(-1)^{\binom{\mathrm{wt}(B)}2}|B\rangle,&
E_C&=C^\dagger A^\dagger B_{\rm mid}^\dagger C A B_{\rm mid}.
                                                               \qquad(46)
\end{aligned}
```

Every axis in each band anticommutes with P. Each list is independent
by its distinct bank labels. The desired complete image is
$`E_CP E_C^\dagger=Z_d\Omega_{a,B}`$. This interface is sufficient
without preparing d: the fixed word
$`L=\prod_j\mathrm{CX}(d,B_j)`$ obeys
$`L(Z_d\Omega_{a,B})L=\Omega_{a,B}`$. Thus
$`LE_C\,\mathrm{SWAP}(d,B_0)`$ would be a source for the ordinary
bank probe, adding m+3 Cliffords in each direction. Its complete echo
would return d. The identities below test the full operator before
making that inference.

### 8.1. Literal native words and their simultaneous ledger

For A, emit the involutory CNOT word

```math
W=\prod_j\mathrm{CX}(B_j,d)
       \prod_{i:j_i=1}\mathrm{CX}(B_j,a_i),
```

then bank H gates, one transversal T band, the bank H gates and the
actual reverse of W. This emits $`\widetilde A=e^{im\pi/8}A`$
with m T gates and $`pm+4m`$ Cliffords. For the middle, let
$`R_d=S_dH_d`$, whose chronological word is H_d,S_d. Emit
$`R_d^\dagger`$, CNOTs from d to every bank wire and from a_i to
B_j when $`j_i=1`$, one transversal bank T band, the reversed
CNOT list and R_d. It emits
$`\widetilde B_{\rm mid}=e^{im\pi/8}B_{\rm mid}`$, with m T
gates and $`pm+2m+6`$ Cliffords. The six probe Cliffords expand
S-dagger as three positive S gates. Actual reversals have the same
common bounds. All physical T targets within each band are distinct.

The complete-graph bank phase C is also paid without a quadratic
Clifford cost. If $`\pi(B)=\mathrm{wt}(B)\bmod2`$, then

```math
C=i^{\mathrm{wt}(B)-\pi(B)}
  =\left(\prod_jS_{B_j}\right)S_{\pi(B)}^\dagger.       \qquad(47)
```

Compute the parity into one existing bank wire, apply S-dagger and
undo the parity CNOTs. This needs no clean or borrowed parity bit.
Including the initial bank S gates, C costs $`3m+1`$ elementary
Cliffords. Its literal reversed native word costs $`5m-1`$ when
inverse S gates are expanded; that actual inverse is used below,
even though the completed operator C is an involution.

The chronological source is

```math
\widetilde B_{\rm mid},\widetilde A,C,
\widetilde B_{\rm mid}^\dagger,\widetilde A^\dagger,C^\dagger.
```

The four T-band scalars cancel literally. The simultaneous common
forward/actual-inverse ledger is

```math
\begin{aligned}
T_{E_C}&=4m,&D_T(E_C)&\le4,&G_{E_C}&=4pm+20m+12,\\
w_{E_C}&=1,&a_{\rm clean}&=0,&\varepsilon_{E_C}&=0.
                                                               \qquad(48)
\end{aligned}
```

The one arbitrary work bit is d. Its motion inside the source is
allowed. The emitted full Pauli echo has $`8m`$ T gates, at most
eight T layers and $`8pm+40m+26`$ Cliffords, including
$`P=S_d^2`$. These are exact ledgers for the candidate; complete
lookup and work return remain the operator claim to be checked.

<a id="two-row-local-frame-covariance"></a>

### 8.2. Two row-local bank frames retain an address translation

There is an obstruction to using only two such frames, before adding
C. It permits arbitrary row-dependent address and helper Paulis,
not merely the linear tags in (46). Suppose the available axes have
the form

```math
\mathcal P_j=A_j\otimes H_j\otimes X_{B_j},\qquad
\mathcal Q_j=B_j^{\rm addr}\otimes K_j\otimes Z_{B_j}.
                                                               \qquad(49)
```

Here A_j and $`B_j^{\rm addr}`$ are arbitrary address Paulis;
H_j and K_j are arbitrary Paulis on probes and other work. Each
axis has no action on any other bank wire. Allow arbitrary angles,
orders and repetitions of these rotations, selector-controlled
versions, and arbitrary unitaries on the work and selectors alone.
There is no layer-count hypothesis and no assumed preparation of the
arbitrary bank.

Fix any nonzero address mask u. Let alpha_j and beta_j be the binary
commutation pairings of $`X_a^u`$ with A_j and
$`B_j^{\rm addr}`$, respectively. The Hermitian Pauli

```math
F_u=X_a^u\otimes I_{\rm work}\otimes
       \prod_j i^{\alpha_j\beta_j}
                      X_{B_j}^{\beta_j}Z_{B_j}^{\alpha_j}
                                                               \qquad(50)
```

commutes with both axes for every j. Its local bank-Z pairing cancels
alpha_j, and its local bank-X pairing cancels beta_j. It also
commutes with every permitted work-only operation. Hence every
completed word in this algebra commutes with F_u, regardless of
its native implementation cost.

Let the initial probe P be a Hermitian Pauli on work only, as in
(46). Then every transported probe also commutes with F_u. This
prevents the complete target $`P\Omega_{a,B}`$ with norm error
less than square root two. To see the exact witness, choose an
address a and a computational bank B such that

```math
B_a+B_{a+u}+\beta_{a+u}=1\pmod2.                    \qquad(51)
```

The two selected bank bits are distinct, so this is always possible.
Choose the probe in a plus eigenstate of P and let b denote the full
product input. The normalized vector
$`(|b\rangle+F_u|b\rangle)/\sqrt2`$ is in the plus eigenspace
of F_u. The two terms have opposite target signs by (51), so the
target sends it to the minus eigenspace. The actual unitary probe
image keeps it in the plus eigenspace. Their normalized outputs
are orthogonal, giving direct error at least $`\sqrt2`$.

The same covariance gives a Gram obstruction with arbitrary common
clean insertion and readout on work only. Since F_u acts trivially
there, the accepted block and its Gram commute with F_u. On the
same principal plus-eigenstate, the target
$`(I+\Omega_{a,B})/4`$ has a minus-eigenspace component of norm
one quarter; the actual Gram has none. Thus its complete Gram
error is at least one quarter. The fixed two-bank XOR wall converts
this statement to $`\Pi_{\rm eq}/2`$ in the canonical interface.
This is an alphabet restriction: new bank frames, coupled bank axes
or interlayers outside the stated algebra are not excluded.

For the two frames in (46), (50) reduces to
$`F_u=X_a^uX_B^{h_u}`$, where $`(h_u)_j=u\cdot j`$.
The third frame supplied by C genuinely removes these particular
symmetries at the native angles. For $`p\ge2`$, h_u has even
weight. Conjugating F_u by the complete-graph C adds bank-Z mask
h_u, up to its constant Pauli phase. It then anticommutes with
P_j whenever $`u\cdot j=1`$. Those axes are independent, so
$`C^\dagger A C`$ does not commute with F_u. Since
$`E_C=(C^\dagger A^\dagger C)B_{\rm mid}^\dagger A B_{\rm mid}`$
and the remaining factors commute with F_u, neither does E_C.
The two-frame argument alone therefore does not reject (46).

<a id="common-probe-paired-state"></a>

### 8.3. The complete third-frame word still has a paired invariant state

A separate witness decides this actual repair. For $`p\ge2`$,
choose $`u=e_0+e_1`$ and pair every row j with j+u. Let the bank
be a product of singlets on those disjoint pairs, and set

```math
|\chi\rangle_a=\frac{|0\rangle+|u\rangle}{\sqrt2},
\qquad
|\Psi\rangle=|\chi\rangle_a|0\rangle_d
       \bigotimes_{\{j,j+u\}}|\psi^-\rangle_{B_jB_{j+u}}.
                                                               \qquad(52)
```

The address has $`X_a^u=Z_a^u=+1`$, since $`u\cdot u=0`$.
Every singlet has both XX and ZZ eigenvalue minus one. Therefore

```math
(P_j+P_{j+u})|\Psi\rangle=0,\qquad
(Q_j+Q_{j+u})|\Psi\rangle=0.                          \qquad(53)
```

The same pair cancellation holds with any probe state: it does not
initialize d inside either band. Within each frame all axes commute,
so the ideal A and $`B_{\rm mid}`$ fix this vector exactly.
Their literal native words multiply it only by $`e^{im\pi/8}`$.
The bank has definite Hamming weight m/2, so C acts by the scalar
$`(-1)^{\binom{m/2}2}`$. These phases cancel between the actual
source calls and their adjoints. Hence E_C and its actual adjoint
fix $`|\Psi\rangle`$ exactly, as does the full echo
$`E_CP E_C^\dagger`$.

Every bank wire has zero Z expectation in the singlet product. The
address is independent of the bank and the probe in (52) has Z
eigenvalue plus one. Thus

```math
\langle\Psi|Z_d\Omega_{a,B}|\Psi\rangle=0,
\qquad
\left\|E_CP E_C^\dagger-Z_d\Omega_{a,B}\right\|
 \ge\sqrt2\quad(p\ge2).                              \qquad(54)
```

This is a complete-state operator witness, not a basis-only test or
an assumed clean-bank sector. The source is native and all four
bands are paid; its complete lookup image is wrong.

The cancellation (53) survives arbitrary uniform angles on each
frame, or more generally angles equal within every paired row
j,j+u in each frame. C may be replaced by any unitary phase
depending only on total bank Hamming weight: it is still scalar on
(52). Actual adjoints cancel that scalar. Arbitrary row-dependent
angles or phases are not covered; equality within the pairs is an
essential hypothesis. In particular, this is not an arbitrary-angle
extension to a general third-frame source or a native-depth bound.

This pauses the uniform common-probe word with total-weight phase
dressing. It shows why removing the two-frame Pauli covariance was
insufficient: a different invariant state survives the complete
circuit. Row-dependent phases can leave this particular witness;
the next subsection tests that larger diagonal-dressing class.

There is a separate Gram extension within this paired-state grammar.
Let $`\Pi`$ project the address onto its common plus eigenspace of
$`X_a^u,Z_a^u`$ and the bank onto the fixed product of singlets,
with every probe and helper unrestricted. Each paired pulse is identity
on this subspace; every total-weight bank phase is scalar there.
Consequently arbitrary repetitions, with probe/helper/selector-only
interlayers and clean insertion or acceptance on those registers,
give an accepted Gram G commuting with the principal projector Pi.
Every $`Z_{B_a}`$ takes the singlet product orthogonally out of it.
Thus $`\Pi\Omega\Pi=0`$ and
$`\|(I-\Pi)(I+\Omega)\Pi/4\|=1/4`$, proving Gram error at least
one quarter. The same statement holds with an arbitrary probe factor
$`Z_d`$ in the target. This is a paired-code restriction, not an
extension of the next result to arbitrary accepted sources.

<a id="static-diagonal-common-probe-spectrum"></a>

### 8.4. A static diagonal dressing cannot repair the relative spectrum

Keep the native angles of A and $`B_{\rm mid}`$, but now allow C
in (46) to be **any** bank-diagonal unitary, without assuming uniform
row phases, a Clifford implementation or a resource bound. Its two
probe fibers share the same change of bank basis. Consequently even
a diagonal C that breaks (52) cannot give the required complete
image at p at least three.

Here is the complete fiber calculation, with all banks arbitrary.
The Clifford W in Section 8.1 changes coordinates to
$`r=d+\pi(B)`$, $`s=a+HB`$, where the p-by-m binary matrix H
has column j at row label j. Let $`\boldsymbol1_m`$ be the
all-ones bank mask and put

```math
M=I+\boldsymbol1_m\boldsymbol1_m^{\mathsf T}+H^{\mathsf T}H,
\qquad T=\prod_j e^{-i\pi X_{B_j}/8}.
```

Then $`WAW=T`$. In the Y_r eigenfiber eta, at fixed co-moving
address s, the middle word is

```math
D_\eta(s)=\prod_j
 \exp\!\left[-\frac{i\pi\eta(-1)^{s\cdot j}}8
                         Z_B^{Me_j}\right],
\qquad D_-(s)=D_+(s)^\dagger.                         \qquad(55)
```

The native T scalars in A and the middle cancel in each completed
source; T and D in (55) use their stated determinant-one rotations.
C is unchanged by W because it is bank-diagonal. Write
$`D=D_+(s)`$, $`Z_\pi=\prod_jZ_{B_j}`$. The complete source
fibers and the off-diagonal probe block of its Pauli image are

```math
\begin{aligned}
E_{C,\eta}&=C^\dagger T^\dagger D_\eta^\dagger C T D_\eta,\\
F_C&=E_{C,+}Z_\pi E_{C,-}^\dagger,\\
F_C&=L_C F_I L_C^\dagger,
\qquad L_C=C^\dagger T^\dagger C T.                  \qquad(56)
\end{aligned}
```

Indeed $`E_{C,\eta}=L_C E_{I,\eta}`$, since C commutes with
D_eta. Thus C changes the common bank basis but leaves the spectrum
of F_C unchanged. The target, in the same Y_r fibers, has both
off-diagonal blocks equal to the Hermitian unitary
$`Q_s=(-1)^{\pi(B)+B_{s+HB}}`$. The norm of the complete image
error on this address fiber is exactly $`\|F_C-Q_s\|`$.

For $`p\ge3`$, binary orthogonality gives $`M^{\mathsf T}M=M^2=I`$
and $`Mh_t=h_t`$, where $`(h_t)_j=t\cdot j`$. The integer weight
difference modulo four is

```math
\mathrm{wt}(My)-\mathrm{wt}(y)
 =2(\boldsymbol1_m+h_{\boldsymbol1_p})\cdot y\pmod4.
```

To verify it, orthogonality makes the difference
additive modulo four; on each basis vector e_j, its value is
$`2(1+\boldsymbol1_p\cdot j)`$. It follows that, literally,

```math
D^2=\left(\prod_j S_{B_j}\right)
             Z_B^{\boldsymbol1_m+h_{s+\boldsymbol1_p}},
\qquad
R_s=Z_\pi D^2=\bigotimes_j r_j,\quad r_j\in\{S,S^\dagger\}.
                                                               \qquad(57)
```

The potential constant phase here is
$`\exp[-i\pi m[s=0]/4]=1`$ at p at least three. No scalar is
dropped. Direct multiplication of (56), using $`D_-=D^\dagger`$,
now gives

```math
(DT)F_I(T^\dagger D^\dagger)
 =T R_s T^\dagger D^{-2}
 =\bigotimes_j G_j,
\qquad G_j=T_Xr_jT_X^\dagger r_j,                    \qquad(58)
```

where $`T_X=e^{-i\pi X/8}`$. In the last equality,
$`D^{-2}=R_s^\dagger Z_\pi`$ and
$`r_j^\dagger Z=r_j`$. Either choice of r_j has
$`G_j^2`$ eigenvalues $`e^{\pm i\beta}`$, with

```math
\cos\beta=\frac{1+2\sqrt2}{4}.
```

This follows by multiplying the two single-qubit matrices: their
squared product has determinant one and trace
$`1/2+\sqrt2`$. Since $`m\ge8`$ is even, choose four unmatched
positive spectral signs and cancel the remaining signs in pairs.
Then $`F_C^2`$ has an eigenvalue $`e^{4i\beta}`$. For any
Hermitian unitary Q, unitarity and $`Q^2=I`$ imply

```math
\begin{aligned}
\|F_C-Q\|&\ge\frac12\|F_C^2-I\|\\
 &\ge|\sin(2\beta)|\gt\frac12.
                                                               \qquad(59)
\end{aligned}
```

The first inequality uses
$`F_C^2-I=F_C(F_C-Q)+(F_C-Q)Q`$. Applying it to Q_s proves
complete Pauli-image error at least one half for every static
bank-diagonal C at p at least three. This includes nonuniform and
nonlinear row phase choices regardless of their synthesis cost.
The stronger square-root-two result (54) retains its separate
uniform-angle and total-weight-phase scope. Equation (59) assumes
the native A and middle angles and makes no new accepted-Gram claim.

Pause the entire static bank-diagonal repair in (46) with the fixed
native A and middle angles. Its first
unpaid operation is a changed relative probe-fiber unitary: a common
bank-basis change cannot turn this spectrum into the required
Hermitian lookup sign. A different interleaving or accepted interface
must supply its own complete operator identity, native schedule and
return ledger. R96 and the complete-frame frontier remain unchanged;
the approximate scalar zero reflection is unaffected.

The [common-probe native fixtures](../../tests/test_common_probe_source.py)
check every p=1 native source column and its actual inverse, and every
p=1 co-moving address/bank fiber column with the specified probe
orientation. At p=2 and p=3 they check the normalized singlet states,
literal pulse phases, source and inverse return, complete Pauli echo
and orthogonal target witness. They also check all 2048 p=3 entries
of (57), both single-qubit spectral orientations in (58), and native
bands and ledgers through p=6. These finite checks challenge the
operator formulas; the all-size covariance, paired-state and
relative-spectrum arguments above are analytic. The proofs and
fixtures have internal audits, not external mathematical review.

<a id="relative-probe-involution-repair"></a>

## 9. A paid relative-block repair still needs the complete lookup action

The first unpaid relative-spectrum operation in Section 8.4 can be
performed by changing the middle pulse. This is a partial construction:
it makes the **relative probe block** an exact Hermitian involution,
but does not give the lookup sectors. The full transported Pauli
$`EPE^\dagger`$ was already a Hermitian involution for every unitary E;
that tautology is different from the additional relative-block property
proved here.

### 9.1. Correct the two middle fibers and resynthesize the result

Retain A, W, T, D and $`Z_\pi=\prod_jZ_{B_j}`$ from Section 8.
Replace the physical middle by

```math
\widehat B=\prod_j
 \exp\!\left[-\frac{i\pi}{8}Z_\pi Z_a^j Z_{B_j}\right].
                                                               \qquad(60)
```

In the co-moving coordinates $`r=d+\pi(B),\ s=a+HB`$, this is
D_s from (55) on **both** Y_r fibers. Indeed the extra bank parity
and the address shift give bank mask
$`\boldsymbol1_m+e_j+H^{\mathsf T}j=Me_j`$. Equivalently, the
old fibers $`D_+=D_s,D_-=D_s^\dagger`$ have been multiplied by

```math
|+Y\rangle\langle+Y|_r\otimes I+
|-Y\rangle\langle-Y|_r\otimes D_s^2.                \qquad(61)
```

For p at least three, (57) identifies the correction in the second
fiber as an address-dependent diagonal Clifford. The implementation
below resynthesizes the corrected middle jointly; it does not assume
that a controlled Clifford has zero native cost.

Here is an exact one-band implementation of (60), valid for every
$`p\ge1`$. On the bank implement
$`V:y\mapsto y+\pi(y)\boldsymbol1_m`$. Since m is even, V is an
involution. Gather all other bank bits into B_0, fan B_0 into those
other bits, and repeat the gather. This uses $`3(m-1)`$ CNOTs and
realizes V on arbitrary inputs. Next add $`a\cdot j`$ into B_j
for every row j, apply one transversal bank T band, and reverse the
entire exposure word. The literal operator is
$`\widetilde{\widehat B}=e^{im\pi/8}\widehat B`$. It costs m T
gates, one T layer and $`pm+6m-6`$ Cliffords, with no helper or
clean input.

More generally, let C be a bank unitary in co-moving coordinates and
let $`C_{\rm phys}=W^\dagger C W`$. Define

```math
\begin{aligned}
\widehat E_C&=C_{\rm phys}^\dagger A^\dagger\widehat B^\dagger
                     C_{\rm phys}A\widehat B,\\
E_s(C)&=C^\dagger T^\dagger D_s^\dagger C T D_s,\\
F_s(C)&=E_s(C)Z_\pi E_s(C)^\dagger,\qquad F_s(C)^2=I.
                                                               \qquad(62)
\end{aligned}
```

The source is now independent of the co-moving probe r. Its transported
physical probe has off-diagonal Y_r blocks both equal to F_s(C).
The required block remains
$`Q_s(y)=(-1)^{\pi(y)+y_{s+Hy}}`$. Thus (62) pays the relative
involution, while $`F_s(C)\approx Q_s`$ remains the complete lookup
obligation. All addresses, banks, probe states and references are
arbitrary in this statement.

For C=I, emit chronologically

```math
\widetilde{\widehat B},\widetilde A,
\widetilde{\widehat B}^{\dagger},\widetilde A^\dagger.
```
The four literal scalars cancel. A common forward/actual-inverse
ledger is

```math
T=4m,\qquad D_T\le4,\qquad G=4pm+20m-12,
\qquad w_{\rm dirty}=1,\quad a_{\rm clean}=0,\quad\varepsilon=0.
                                                               \qquad(63)
```

The one arbitrary work bit is d. The safe complete-echo ledger is
$`8m`$ T gates, depth at most eight, and $`2G+2`$ Cliffords,
using the actual source inverse and $`Z_d=S_d^2`$. For an emitted
Clifford $`C_{\rm phys}`$, add its forward and actual-inverse
Clifford costs to G. No free implementation is granted to an arbitrary
C in the analytic statements below. These are candidate ledgers,
not completed lookup or work-return claims.

Without C, the repair still commutes with
$`X_a^uX_B^{h_u}`$, $`(h_u)_j=u\cdot j`$, for every nonzero u
when $`p\ge2`$. The parity of h_u is even, and its rowwise signs
cancel those of the address in (60). The witness from Section 8.2
therefore gives complete image error at least square root two. This
does not undo the relative-block repair; it identifies a separate
missing sector action.

<a id="relative-repair-diagonal-modulus"></a>

### 9.2. Diagonal dressing cannot repair the transition magnitudes

Keep the fixed native angles in (60) and A. Allow C to be any static
bank-diagonal unitary, with no cost restriction. Put
$`H_0=T Z_\pi T^\dagger`$. Conjugating the desired identity in
(62) by TC reduces its two sides to

```math
(D_s^\dagger C)H_0(C^\dagger D_s),\qquad
R_s=T Q_s T^\dagger.                                \qquad(64)
```

Every matrix entry of H_0 has magnitude $`2^{-m/2}`$: each local
factor is $`(Z-Y)/\sqrt2`$. The left side of (64) has exactly
these magnitudes, regardless of C. The target has a different
amount of squared matrix magnitude at large bank displacement.
This gives a constant error floor, not merely a small individual-entry
defect.

For $`p\ge3`$, the exact Walsh expansion is

```math
Q_s=\frac1m\sum_{a,t\in\mathbb F_2^p}
 (-1)^{t\cdot(s+a)}Z_B^{\boldsymbol1_m+h_t+e_a}.
                                                               \qquad(65)
```

All $`m^2`$ masks are distinct: a nonzero h_t has weight m/2,
which exceeds the weight of $`e_a+e_b`$. For t=0 the m masks
have weight m-1. For each nonzero t, half have weight m/2-1 and
half have weight m/2+1.

Normalize Hilbert--Schmidt norm by dividing its square by $`2^m`$.
Restrict matrix entries to $`|y+y'|\ge k`$, where the addition is
binary and the bars denote Hamming weight. A conjugated Pauli mask
of weight w distributes its squared norm over displacements as
$`\mathrm{Bin}(w,1/2)`$. Distinct masks in (65) remain orthogonal
after this restriction: their expanded Z/Y Pauli terms have different
supports, even when their Y displacement masks coincide. Consequently
the exact restricted squared norms of the left and right sides of
(64) are

```math
\begin{aligned}
P_H&=\Pr[\mathrm{Bin}(m,1/2)\ge k],\\
P_R&=\frac1m\Pr[\mathrm{Bin}(m-1,1/2)\ge k]\\
&\quad+\frac{m-1}{2m}\left(
 \Pr[\mathrm{Bin}(m/2-1,1/2)\ge k]
 +\Pr[\mathrm{Bin}(m/2+1,1/2)\ge k]\right),\\
\|F_s(C)-Q_s\|&\ge\sqrt{P_H}-\sqrt{P_R}.
                                                               \qquad(66)
\end{aligned}
```

The last inequality uses the reverse triangle inequality after the
entry restriction, and the fact that operator norm bounds normalized
Hilbert--Schmidt norm. Entry restriction is used only for that norm;
it is not asserted to contract operator norm.

For $`m\ge128`$, choose k=3m/8. Chebyshev's inequality gives
$`P_H\ge1-16/m`$. Each of the two smaller-binomial tails is at
most
$`(8m+16)/(m-4)^2\le16/m`$, while the first term of P_R is at
most 1/m. Thus $`P_R\le17/m`$ and (66) is greater than one half.
The four remaining powers of two are finite rational binomial sums:

| m | k | Exact lower bound on P_H | Exact upper bound on P_R |
|---|---|---|---|
| 8 | 4 | $`3/5`$ | $`3/20`$ |
| 16 | 7 | $`3/5`$ | $`3/20`$ |
| 32 | 12 | $`3/5`$ | $`3/20`$ |
| 64 | 24 | $`3/5`$ | $`3/20`$ |

All eight comparisons are strict and follow by integer cross
multiplication of the displayed binomial sums. For example at m=8,
$`P_H=163/256`$ and $`P_R=37/256`$. The same rational comparisons
are checked by the fixtures. These bounds give
$`\sqrt{P_H}-\sqrt{P_R}>\sqrt{3/20}>1/3`$. Therefore

```math
\|F_s(C)-Q_s\|>\frac13\qquad(p\ge3)                 \qquad(67)
```

for every static bank-diagonal C, at every co-moving address s.
The relative block in (62) is already an involution; this new
obstruction concerns transition magnitudes, not its spectrum.
It assumes the stated native pulse angles and source word, and
makes no claim about arbitrary higher-band circuits or accepted
Gram interfaces.

<a id="relative-repair-hadamard-mask"></a>

### 9.3. A non-diagonal Hadamard mask also has a complete witness

A concrete change outside Section 9.2 is
$`C=\prod_{j\in S}H_{B_j}`$ for any nonempty set S of bank rows.
The physical interlayer $`W^\dagger C W`$ makes genuine joint
address, probe and bank motion. Its native cost is
$`pm+2m+|S|`$ Cliffords, with no T gates or additional work.
Both its forward word and actual inverse have this cost. The source
in (62) therefore has

```math
T=4m,\qquad D_T\le4,\qquad
G=6pm+24m+2|S|-12,\qquad w_{\rm dirty}=1,
\quad a_{\rm clean}=0,\quad\varepsilon=0.             \qquad(68)
```

This includes masks selected by a nonperiodic quadratic function of
the row label. It leaves the relative image Hermitian, but it does
not supply the complete lookup action.

To see this for every nonempty S, choose any q in S. Conjugate the
actual F_s(C) and target Q_s by TC. The two operators become

```math
A_s=D_s^\dagger C H_0 C^\dagger D_s,\qquad
R_s=T C Q_s C^\dagger T^\dagger.
```

The local q factor in $`C H_0 C^\dagger`$ is
$`(X_q+Y_q)/\sqrt2`$. Hence A_s anticommutes with Z_q;
the diagonal D_s does not change that anticommutation.
At s=0, both bank strings 0 and e_q have Q_0 sign plus one,
including q=0. In fact the singleton bank e_q has odd parity and
its selected bit at address $`He_q=q`$ equals one. Set

```math
|\psi\rangle=\frac{|0\rangle+e^{-i\pi/4}|e_q\rangle}{\sqrt2},
\qquad |\phi\rangle=T C|\psi\rangle.                \qquad(69)
```

The q factor of psi is a plus eigenstate of
$`C^\dagger T^\dagger Z_q T C=(X_q-Y_q)/\sqrt2`$;
all other bank bits in psi are zero. Therefore phi has Z_q
eigenvalue plus one and $`R_0|\phi\rangle=|\phi\rangle`$.
The unitary A_0 sends phi into the minus eigenspace of Z_q.
Their normalized outputs are orthogonal, proving
$`\|F_0(C)-Q_0\|\ge\sqrt2`$ for every nonempty S and every
$`p\ge1`$. Equivalently, psi itself is the complete relative-block
witness before the TC change of basis. Taking the co-moving probe
in minus Y turns this into a full native Pauli-image witness with
plus-Y output; no bank preparation is assumed by the candidate.

<a id="relative-repair-cubic-hadamard"></a>

### 9.4. An entangling cubic/Hadamard word also collapses

A paid nonlinear, non-diagonal interlayer can leave both preceding
premises and still have the wrong complete action. For p at least two,
partition the bank into four-row faces along the first two row-label
coordinates, and set

```math
P|B\rangle=(-1)^{f(B)}|B\rangle,\qquad
f(B)=\sum_j B_jB_{j+e_0}B_{j+e_1}\pmod2,\qquad
C=P H_B.
```

Each face contributes its four distinct cubic triples exactly once.
The chronological order for C is the bank Hadamards followed by P.
Use this C in the corrected family (62). To pay it, the identity

```math
4abc=a+b+c-(a\mathbin\oplus b)-(a\mathbin\oplus c)
 -(b\mathbin\oplus c)+(a\mathbin\oplus b\mathbin\oplus c)
 \pmod8
```

emits each CCZ with seven T-type gates and ten CNOTs by separately
computing and uncomputing its three pair parities and its triple parity.
The four triples of each face run serially; distinct faces use disjoint
wires and run in parallel. Thus P and its actual native inverse have
$`T=7m`$, depth at most 28 and at most $`10m`$ Cliffords, with no
helpers. Including the bank Hadamards and both W walls gives
$`T(C_{\rm phys})=7m`$, depth at most 28 and
$`G(C_{\rm phys})\le pm+13m`$. The complete source and its actual
reverse therefore have the simultaneous safe ledger

```math
T=18m,\qquad D_T\le60,\qquad G\le6pm+46m-12,
\qquad b_{\rm extra}=1,\qquad a_{\rm clean}=0.
```

The extra bit is the arbitrary probe. The complete Pauli echo has
$`T=36m`$, depth at most 120 and at most $`12pm+92m-22`$ Cliffords.
These are exact native candidate costs, not a lookup or probe-return
certificate. The inverse of P is emitted as the actual reversed word;
its equality to P as a completed operator does not delete either call.

There is an all-size collapse for p at least three. Write
$`X_\pi=X^{\otimes m}`$ and $`S_\pi=\prod_jS_{B_j}`$. Complementing
all four bits of a face changes its cubic sum by that face's parity,
so $`f(B+\boldsymbol1_m)+f(B)=\pi(B)`$. Consequently

```math
P X_\pi P=X_\pi Z_\pi,\qquad
C H_0 C^\dagger=e^{im\pi/4}X_\pi S_\pi.
```

The latter identity uses the literal one-wire factor
$`(X+Y)/\sqrt2=e^{i\pi/4}XS^\dagger`$. Also
$`M\boldsymbol1_m=\boldsymbol1_m`$, so
$`D_s(B+\boldsymbol1_m)=\overline{D_s(B)}`$ and
$`D_s^\dagger X_\pi D_s=X_\pi D_s^2`$. Applying (57), including
its scalar, reduces the transformed operator from Section 9.3 to a Pauli:

```math
A_s=D_s^\dagger C H_0 C^\dagger D_s
 =X_\pi Z_B^{h_{s+\boldsymbol1_p}}\qquad(p\ge3).
```

Here $`e^{im\pi/4}=1`$ and $`S_\pi^2=Z_\pi`$ were retained in
the cancellation. Choose $`s=\boldsymbol1_p`$. Then A_s is $`X_\pi`$,
which commutes with T, and the complete relative image is exactly

```math
F_s(C)=C^\dagger T^\dagger X_\pi TC
      =H_B P X_\pi P H_B=Y^{\otimes m}.
```

It sends the zero bank to the all-one bank, with literal phase
$`i^m=1`$, while $`Q_s|0^m\rangle=|0^m\rangle`$. The outputs
are orthogonal, proving complete image error at least square root two
for every p at least three. This is a whole relative-block identity;
the probe-fiber embedding gives the same full native operator witness.
The exact p=2 fixture separately gives
$`\langle0|F_{11}(C)|0\rangle=1/8`$, hence column error
$`\sqrt7/2`$. That finite calculation does not invoke (57) outside
its p-at-least-three scope.

Thus this specific entangling cubic/Hadamard interlayer is also paused.
Its common-fiber correction and native controls are paid; the complete
lookup action fails. This does not exclude other non-diagonal,
entangling interlayers or general native sources.

This pauses the two dressing classes and the specified cubic word. A changed
joint bank operation or accepted interface must supply the complete
lookup-dependent transition action. Merely repairing the relative
involution, choosing diagonal row phases, or changing the Hadamard
mask does not do so. R96 and the simultaneous complete-frame
resource frontier remain unchanged. These source obstructions do
not contradict the approximate scalar zero-reflection construction
and do not give an unrestricted native-depth lower bound.

The [relative-probe fixtures](../../tests/test_relative_probe_repair.py)
check all p=1 source columns and their actual inverse for the basic
repair. Selected complete p=3 co-moving columns check the relative
block, its square and the literal native echo. At p=2 and p=3 they
check the non-diagonal source and actual inverse on selected inputs,
the normalized witness (69), native/fiber agreement, zero target
overlap and squared output distance two; a p=1 fixture also checks
the selected-row q=0 corner. The displacement test uses exact rational
binomial sums at m=8,16,32,64, without a large dense circuit fixture.
The cubic/Hadamard fixture checks every native four-row phase column
and its actual inverse, the p=2 scalar overlap, three p=3 relative
columns, complete native/fiber echo agreement at both sizes and the
emitted count ledgers. It uses the actual reversed phase word.
Native band and Clifford ledgers are checked through p=6. These
finite checks challenge the identities; the all-size fiber,
displacement, anticommutation and cubic-collapse arguments are analytic results
with internal audits, not external mathematical review.

<a id="probe-asymmetric-paid-trials"></a>

## 10. Paid probe-asymmetric repairs: three finite failures and a coordinate boundary

The next three words leave the static bank-diagonal premise of Section 8.4.
Their native controls are paid, but their tested relative blocks are not
Hermitian. The witnesses below reject the specified families at one finite
size; they are not all-size exclusions of asymmetric interlayers. In
particular, patching that size would require a new complete identity, rather
than being ruled out by these fixtures.

### 10.1. The changed fibers and their actual relative word

Retain the co-moving coordinates r,s, the rotations T and D=D_+(s), and
$`Z_\pi=\prod_j Z_{B_j}`$ from (55). Replace the dressing in (46) by a
probe-controlled bank unitary

```math
\mathcal C_V
 =\frac{I+Y_r}{2}\otimes I
  +\frac{I-Y_r}{2}\otimes V.
```

Here V is a specified full-space bank unitary. In the original physical
coordinates the dressing is $`C_V=W\mathcal C_VW`$, with the actual Clifford
W from Section 8.1. The completed source is still
$`E_{C_V}=C_V^\dagger A^\dagger B_{\rm mid}^\dagger C_VAB_{\rm mid}`$.
Thus its two literal fibers are

```math
\begin{aligned}
E_+&=T^\dagger D^\dagger TD,\\
E_-&=V^\dagger T^\dagger DVT D^\dagger,\\
F_V&=E_+Z_\pi E_-^\dagger
 =T^\dagger D^\dagger T R_s T^\dagger V^\dagger D^\dagger TV,
\qquad R_s=Z_\pi D^2.                                    \qquad(70)
\end{aligned}
```

No commutation between V and D was used. The target is the same Hermitian
lookup sign $`Q_s=(-1)^{\pi(B)+B_{s+HB}}`$. In this probe basis the two
Pauli-image blocks are $`F_V,F_V^\dagger`$; their complete distance from the
target blocks equals $`\|F_V-Q_s\|`$. The initial target probe remains
arbitrary throughout.

The rotations in (70) are determinant-one operators. A and the middle
word are emitted with their literal native scalars, as in Section 8.1;
the two scalar pairs cancel in the completed source. V and its actual
adjoint occur explicitly and retain every relative phase. None of the
following calculations replaces the inverse of a relative-phase shear
by its identical binary permutation.

### 10.2. Three native choices of V

**Controlled bank parity.** Take $`V=Z_\pi`$. The controlled operation
$`\mathcal C_V`$ is Clifford: change the probe's Y basis to Z, apply CZ
from the probe to each bank wire, and undo the basis change. The W walls
are also Clifford. This repair therefore retains four source T bands,
$`T=4m`$, $`G=O(pm)`$, one arbitrary probe and no extra helper. It changes
the two probe fibers differently and is not a static bank-diagonal C.

**Controlled syndrome involution.** For $`p\ge3`$, take
$`V=U_M`$, where $`U_M|B\rangle=|MB\rangle`$ and M is the binary matrix
in Section 8.4. Put q=p+1 and let L have first row $`\boldsymbol1_m^\mathsf T`$
and remaining rows H. Then

```math
N=M-I=L^\mathsf TL,\qquad LL^\mathsf T=0,\qquad
\mathrm{rank}L=\mathrm{rank}N=q,
\qquad N^2=0.                                           \qquad(71)
```

The rows of L are independent affine truth tables. Their pairwise inner
products are even for p at least three, including each row's square;
this proves the zero Gram matrix. Since L is onto and $`L^\mathsf T`$
is injective, N has rank q. The condition $`2q\le m`$ holds, with equality
at p=3. No claim of this involution or rank formula is made for p=1,2.

There is a concrete binary basis exposing q disjoint CNOTs. Index the rows
of L by 0,...,p. Write $`u_i=L^\mathsf Te_i`$; take $`v_0=e_0`$ and
$`v_i=e_0+e_{e_{i-1}}`$ for $`1\le i\le p`$, where the outer subscripts
name physical bank rows. Then $`Lv_i=e_i`$. Extend
$`\{u_i\}`$ to a basis of $`\ker N`$ using vectors k_j. The ordered
vectors $`u_i,v_i,k_j`$ are independent: applying N to a relation first
kills every coefficient of v_i, after which independence in the kernel
finishes the proof. In this basis
$`M:u_i\mapsto u_i,\ v_i\mapsto v_i+u_i,\ k_j\mapsto k_j`$.
Gaussian elimination emits the fixed binary basis change and its actual
reverse with $`O(m^2)`$ CNOTs.

Conditioning these q disjoint CNOTs on the probe gives q Toffolis with a
shared control. Use the retained [returned private-phase Toffoli](../../docs/DIRTY_SUM_COMPRESSION.md#read-only-affine-controls-with-private-phase-work): eight T
gates, at most twelve Cliffords, at most eight T layers, and one arbitrary
phase helper per pair. Its physical T targets are private, so all q
Toffolis have an actual parallel schedule. Including the probe-basis and
W walls gives $`T(C_V)\le8q`$, depth at most eight,
$`G(C_V)=O(m^2)`$, and q exactly returned dirty helpers. The full source
therefore has

```math
T\le4m+16(p+1),\qquad D_T\le20,\qquad
G=O(m^2),\qquad b_{\rm extra}\le p+2,
\qquad a_{\rm clean}=0.                                \qquad(72)
```

Here the extra width includes the arbitrary probe. The source may move
that probe and the principal address/bank; the private phase helpers
return exactly. This does not claim the desired final indicator return.

**Controlled paired coefficient shear.** Pair every row j with j+e_0,
choosing the lower row by $`j_0=0`$. Let A and B denote the lower and upper
banks of q=m/2 bits. Set

```math
V|A,B\rangle=i^{\mathrm{wt}(B)}|A+B,B\rangle.             \qquad(73)
```

The complete controlled operation is exactly the
[two-band coefficient shear](WEAK_PAULI_ORBIT_COMPRESSION.md#complete-coefficient-shear)
$`G_c^{(m/2)}`$ with $`c=(I-Y_r)/2`$, after the paid probe basis change.
It needs no helpers. The W walls contain $`pm+2m`$ CNOTs in total;
the shear has $`2m`$ T gates, two T bands and $`4m`$ Cliffords; the two
probe basis changes cost six Cliffords in the positive-S alphabet.
Consequently common forward/actual-inverse bounds are

```math
\begin{aligned}
T(C_V)&=2m,&D_T(C_V)&\le2,&G(C_V)&\le pm+6m+6,\\
T(E_{C_V})&=8m,&D_T(E_{C_V})&\le8,
&G(E_{C_V})&\le6pm+24m+24.                              \qquad(74)
\end{aligned}
```

There is one arbitrary probe and no other work. The complete Pauli echo
has $`16m`$ T gates, depth at most sixteen and at most
$`12pm+48m+50`$ Cliffords, including $`Z_d=S_d^2`$. These are exact
native ledgers for the candidate, with arbitrary coherent bank, address
and reference inputs; they do not certify its desired lookup action.

### 10.3. Exact finite witnesses reject the three displayed repairs

At p=3 and s=0, evaluate the complete matrix coefficient
$`\langle0^8|F_V|0^8\rangle`$ using (70). Every coefficient belongs to
$`\mathbf Q(\sqrt2,i)`$. For reproducibility, on a bank value y the
D phase is

```math
D(y)=\omega^{\mathrm{wt}(My+h_s)-m/2},\qquad
\omega=(1+i)/\sqrt2,\qquad (h_s)_j=s\cdot j.
```

The literal one-wire native X rotation is $`HTH`$. At m=8 its transversal
word is minus the determinant-one T; its four occurrences in (70) cancel
that scalar. Apply, chronologically, V, T, D-adjoint, V-adjoint,
T-adjoint, R_s, T, D-adjoint, T-adjoint. This computes the actual relative
column, including the inverse phase in (73), and gives

| V | Exact imaginary part |
| --- | --- |
| $`Z_\pi`$ | $`-21(2+\sqrt2)/4096`$ |
| $`U_M`$ | $`(846-37\sqrt2)/16384`$ |
| Paired shear (73) | $`(686-1213\sqrt2)/65536`$ |

Each imaginary part is nonzero, so each tested F_V is non-Hermitian.
In particular its distance from Q_s is at least the absolute value of
that imaginary part. These exact small fixtures refute the stated
uniform identities. They neither exclude all larger instances nor every
possible repair of the finite case. Additional changes would need a new
whole-space lookup certificate.

### 10.4. Completed classical coordinate rounds retain a logarithmic requirement

A different proposal would realize Q_s by conjugating a fixed bank Pauli
through completed low-degree monomial coordinate changes. For every
fixed s the Boolean exponent

```math
q_s(B)=\pi(B)+B_{s+HB}
```

has degree exactly p+1. The upper bound follows by expressing each
address equality as a product of p affine bank functions. For the lower
bound restrict the bank to row zero and the p basis rows e_i. Write
$`z_i=B_{e_i}`$. Then $`HB=z`$, and the selected-row-zero contribution is
$`B_0[z=s]`$, whose monomial $`B_0\prod_i z_i`$ has coefficient one.
All other selected-row contributions are independent of B_0, and the
parity term is linear, so that monomial cannot cancel.

If a common-coordinate realization $`Q_s=V_s Z^\ell V_s^\dagger`$
uses r monomial rounds whose forward and inverse binary coordinate maps
have degree at most k, affine boundaries included, then degree composition
gives $`k^r\ge p+1`$. Arbitrary monomial phases cancel in this conjugation.
For this exact diagonal-sign grammar, an operator error strictly below
two already forces the Boolean signs to agree. Thus quadratic completed
rounds still require logarithmically many rounds. This specializes the
[existing classical-round boundary](JOINT_ROUTING_INDICATOR_INTERFACES.md#5-classical-round-boundaries-retain-a-logarithmic-degree-requirement)
to the current co-moving target. This counts completed coordinate rounds;
it does not preclude native synthesis across those rounds or give an
unrestricted T-depth bound.

The argument does not apply to a nonlinear shear inserted inside a
nonmonomial interference word, nor to unrelated probe fibers whose
relative phases supply the target. The trials above test precisely such
changed fiber words; their native costs are paid and their finite
identities fail. The first unpaid operation remains a complete relative
lookup action at sublogarithmic native depth. No complete-frame frontier
moves, and the scalar zero-reflection theorem is unchanged.

The [asymmetric-fiber fixtures](../../tests/test_probe_asymmetric_trials.py)
check the three exact p=3 relative coefficients, the syndrome involution's
disjoint-CNOT coordinates through p=6, and the top-degree restriction on
every small fixed syndrome through p=6. The finite coefficient failures
are distinguished from the analytic native schedules and degree argument.
They have internal audits, not external mathematical review.

<a id="affine-fourier-complete-transition-gap"></a>

## 11. Affine bank mixing around a full Fourier layer still misses lookup

A full bank Hadamard layer can be combined with an arbitrary invertible
linear change of bank coordinates, instead of changing only a Hadamard
mask. This leaves the diagonal-dressing premise of Section 9.2 and makes
nonlocal bank transitions. Nevertheless the corrected source (62) has a
complete error floor of square root two for both orders of those layers.
The obstruction persists under arbitrary X-basis phase functions in T
and arbitrary computational-diagonal middle phases. It concerns the
displayed source word; it is not a lower bound for arbitrary native
circuits or for two genuinely nonmonomial mixing stages.

### 11.1. Complete transition support for both affine/Fourier orders

Let $`p\ge2`$, $`m=2^p`$, and let L be any phase-free affine bank
permutation

```math
L|y\rangle=|Ay+b\rangle,\qquad A\in\mathrm{GL}(m,2).
```

Take either $`C=H_B L`$ or $`C=L H_B`$, where
$`H_B=H^{\otimes m}`$. The affine offset b is arbitrary; no initialization
of bank bits is allowed. For this subsection allow T to be any unitary
diagonal in the joint X basis, and $`D_s`$ to be any computational-diagonal
unitary. Retain the same source form
$`E_s=C^\dagger T^\dagger D_s^\dagger C T D_s`$ and
$`F_s=E_sZ_\pi E_s^\dagger`$. This analytic generality supplies no free
native implementation of T or $`D_s`$. Write $`\mathbf1=\boldsymbol1_m`$
and $`H_0=T Z_\pi T^\dagger`$. Conjugating actual and target by TC gives

```math
\mathcal A_s=D_s^\dagger C H_0 C^\dagger D_s,
\qquad
\mathcal R_s=T C Q_s C^\dagger T^\dagger.                 \qquad(71)
```

Because T commutes with every X-type Pauli while $`Z_\pi`$
anticommutes with the odd-weight ones,

```math
H_0 X^u H_0=(-1)^{\mathbf1\cdot u}X^u.
```

The affine L conjugates $`X^u`$ to $`X^{Au}`$ without a phase.
It follows that $`C H_0 C^\dagger`$ is computational-monomial,
with a single fixed displacement v in every column:

| C | Displacement v |
|---|---|
| $`H_B L`$ | $`A^{-\mathsf T}\mathbf1`$ |
| $`L H_B`$ | $`A\mathbf1`$ |

For the first row, conjugation of each $`Z^u`$ by $`C H_0 C^\dagger`$
returns $`(-1)^{(A^{-\mathsf T}\mathbf1)\cdot u}Z^u`$; the
joint computational eigenspaces are one dimensional. The second row
follows by first using $`H_B`$ and then applying L. The affine offset can
change the phase of a column but not its displacement. The diagonal
$`D_s`$ in (71) likewise changes only the column phase. Thus for every y,

```math
\mathcal A_s|y\rangle=\gamma_{s,y}|y+v\rangle,
\qquad |\gamma_{s,y}|=1.                                  \qquad(72)
```

The Walsh identity (65), which remains valid even when its masks are not
distinct at p=2, has support among

```math
\mathbf1+h_t+e_a,\qquad a,t\in\mathbb F_2^p.
```

None of these masks equals $`\mathbf1`$: $`h_0`$ is zero, and each
nonzero $`h_t`$ has weight $`m/2\ge2`$, whereas $`e_a`$ has weight one.
For $`C=H_B L`$, the corresponding X-displacement masks of
$`C Q_s C^\dagger`$ are $`A^{-\mathsf T}(\mathbf1+h_t+e_a)`$;
for $`C=L H_B`$ they are $`A(\mathbf1+h_t+e_a)`$. Affine
phases affect coefficients only. Neither support contains v.

Both conjugated targets are sums of X-type Paulis, so they commute with
the X-diagonal T. Therefore $`\mathcal R_s=C Q_s C^\dagger`$
and

```math
\langle y+v|\mathcal R_s|y\rangle=0
\quad\text{for every }s,y.
```

Equations (71)--(72) give normalized orthogonal output columns. Hence

```math
\|F_s(C)-Q_s\|\ge\sqrt2
\qquad\text{for every }s,\ p\ge2.                         \qquad(73)
```

Before the TC change of basis, an explicit witness is
$`C^\dagger T^\dagger|y\rangle`$. Embedding it in the minus-Y
co-moving probe fiber gives the same complete native Pauli-image error.
These witnesses test permitted arbitrary inputs; no state preparation
resource is supplied to the source.

### 11.2. The native costs are paid, but do not repair the transition

Suppose L is emitted using $`g_L`$ elementary Cliffords in a common
forward/actual-inverse bound, with CNOTs, swaps expanded into CNOTs,
and X expanded as HSSH. Both displayed choices of C have a physical
implementation $`C_{\rm phys}=W^\dagger C W`$ with

```math
G(C_{\rm phys}),\ G(C_{\rm phys}^\dagger)
\le pm+3m+g_L.
```

It uses no T gates or extra helpers. For the native product angles and
corrected middle of Section 9.1, the original X-axis source bands
therefore give the complete candidate ledger

```math
T=4m,\qquad D_T\le4,\qquad
G\le6pm+26m-12+2g_L,
\qquad b_{\rm extra}=1,\quad a_{\rm clean}=0.              \qquad(74)
```

The extra bit is the arbitrary physical probe. The actual Pauli echo has
T-count 8m, depth at most eight, and Clifford count at most twice the
bound in (74) plus two. Gaussian elimination supplies an emitted
$`g_L=O(m^2)`$ for any specified affine map, so these candidates fit the
permitted Clifford budget. The failure in (73) is their complete
operator action, not an unpaid control or uncounted native gate.

### 11.3. No monomial dressing gives exact lookup at the native angles

Return here to the native product $`T=\prod_j e^{-i\pi X_j/8}`$
and the corrected middle $`D_s`$ of Section 9.1. A separate exact statement
covers every computational-monomial C on the bank, including nonlinear basis permutations and arbitrary diagonal
phases, without granting it a free implementation. Such a C preserves
the fact that every entry of $`H_0`$ has modulus $`2^{-m/2}`$. Thus the
actual $`\mathcal A_s`$ in (71) has that same modulus everywhere.
The transformed target is
$`\mathcal R_s=T\mathrm{diag}(f)T^\dagger`$, where
$`\mathrm{diag}(f)=C Q_s C^\dagger`$ and f takes values
in $`\{+1,-1\}`$.

Let $`\widehat f(S)`$ be its normalized Boolean Walsh coefficients.
A direct product of the one-qubit transition probabilities of T gives

```math
\langle y|\mathcal R_s|y\rangle
 =\sum_{S\subseteq[m]}\widehat f(S)\,2^{-|S|/2}(-1)^{S\cdot y}.
```

Parseval consequently yields

```math
2^{-m}\sum_y|\langle y|\mathcal R_s|y\rangle|^2
 =\sum_S\widehat f(S)^2 2^{-|S|}\ge2^{-m}.                \qquad(75)
```

Here $`\sum_S\widehat f(S)^2=1`$, and equality in (75) is
possible only when all Fourier mass is at $`S=[m]`$. If
$`F_s(C)=Q_s`$ exactly, the flat actual matrix forces equality, hence

```math
C Q_s C^\dagger=\pm Z_\pi.                              \qquad(76)
```

Thus a monomial dressing that completed this exact four-band word would
already transport the full lookup sign to parity. The remaining flat
matrix identity is stronger: it forces the dressing's basis permutation
to be affine, so no monomial dressing completes the lookup for p at least
two.

To prove this, write $`S_B=\prod_j S_{B_j}`$. The native product has
the literal identity

```math
H_0=S_B^\dagger H_B S_B.
```

Let P be the basis permutation underlying C, with
$`P|y\rangle=|\sigma(y)\rangle`$. Substituting (76) into (71)
gives $`D_s^\dagger C H_0 C^\dagger D_s=\pm H_0`$. Absorbing
all diagonal factors on both sides, including the arbitrary phases of C,
therefore gives a diagonal unitary J such that

```math
J P H_B P^\dagger J^\dagger=\pm H_B.                    \qquad(77)
```

Take the multiplicative rectangle ratio of entries at
$`(x,y),(x,0),(0,y),(0,0)`$. Every entry is nonzero; the four
row/column phases of J and the common sign cancel. The Walsh kernel
$`(H_B)_{xy}=2^{-m/2}(-1)^{x\cdot y}`$ then gives, over the
binary field,

```math
g(x)\cdot g(y)=x\cdot y,
\qquad g(x)=\sigma^{-1}(x)+\sigma^{-1}(0).              \qquad(78)
```

The map g is bijective and fixes zero. For all x,x',y, (78) implies
$`[g(x+x')+g(x)+g(x')]\cdot g(y)=0`$. Surjectivity in y makes
the bracket zero. Hence g is linear and orthogonal, and sigma is affine.
But (76) would then make the exponent of $`Q_s`$ affine as well. This
contradicts its degree p+1 from Section 10.4 whenever p is at least two.
Consequently **no computational-monomial C**, even with a nonlinear
permutation and arbitrary phases, supplies exact complete lookup in this
source word at the fixed native T angles. No synthesis-cost assumption
on C was used.

This exact screen does **not** give a dimension-independent approximation
gap for arbitrary monomial dressings. The diagonal restriction only gives

```math
\|F_s(C)-Q_s\|\ge
\left(\sum_S\widehat f(S)^2 2^{-|S|}\right)^{1/2}-2^{-m/2}.
```

This bound can be small; it must not be promoted to a fixed-accuracy exclusion.
The robust result (73) has the separate affine/Fourier hypotheses and
already excludes angle tuning or a changed diagonal middle in that word.
Neither result excludes a genuinely different nonmonomial interlayer,
a nondiagonal changed middle or an accepted-Gram source. R96 remains
unchanged.

The [affine/Fourier fixtures](../../tests/test_affine_fourier_repair.py)
check both layer orders with a nontrivial CNOT basis and affine offset,
fixed transition displacement, exact target support exclusion, actual
p=3 source inverses and native echoes on the coherent witnesses, and
the native ledger through p=6. A separate changed-word fixture uses
correlated X-basis cubic phases and a new diagonal middle. The monomial
fixture includes cubic diagonal phases and checks the exact noise
identity and flat actual columns. A nonlinear permutation that does
map the lookup sign to parity is rejected by the exact Walsh-rectangle
identity, testing the second step of the monomial proof. These finite checks challenge the
analytic statements; they are not a general synthesis or depth proof.

<a id="two-clifford-mixer-screen"></a>

## 12. Two Clifford mixing stages: a transition screen and its exact exception

This screen concerns a specified complete bank-source grammar. It does not
bound arbitrary native T-depth, the earlier four-band X-rotation family,
or sources using additional active work registers or intermediate address
motion. Fix $`p\ge3`$, $`m=2^p`$, and a co-moving address $`s`$. All $`m`$ bank
bits and any passive reference are arbitrary. Put

```math
P=Z_\pi,\qquad Q_s|y\rangle=(-1)^{q_s(y)}|y\rangle,
\qquad q_s(y)=\pi(y)+y_{s+Hy}.
```

Allow the actual bank word

```math
E=D_2 C_2 D_1 C_1 D_0,\qquad F=EPE^\dagger,             \qquad(T1)
```

where $`C_1,C_2`$ are arbitrary bank Cliffords and the $`D_i`$ are arbitrary
bank-diagonal unitaries. No native implementation of a general $`D_i`$ is
granted. The screen is therefore valid even after giving this diagonal
freedom for free. The rightmost operation acts first; the inverse is the
actual adjoint. Diagonal boundaries $`D_0,D_2`$ may have literal phases.

### 12.1. A Fourier-coset identity for the complete diagonal

Write $`C_1PC_1^\dagger`$ as a Hermitian Pauli with binary X label $`a`$. Since
$`D_1`$ is diagonal, $`D_1(C_1PC_1^\dagger)D_1^\dagger`$ has matrix entries only
between $`y`$ and $`y+a`$. Its Pauli expansion is contained in one affine Pauli
coset with fixed X label $`a`$. Let

```math
\mathcal A=C_2\langle Z_{B_0},\ldots,Z_{B_{m-1}}\rangle C_2^\dagger,
\qquad V=\{u:Z^u\text{ belongs to }\mathcal A
                         \text{ up to scalar}\}.
```

Define $`h=m-\dim V`$; equivalently $`h`$ is the rank of the binary X-label
projection of this maximal commuting Pauli group. After $`C_2`$, the Pauli
support is in one affine coset of $`\mathcal A`$. Its intersection with the diagonal
Paulis is either empty or an affine binary coset $`u_0+V`$. $`D_0`$ commutes with
$`P`$, while $`D_2`$ leaves every diagonal entry unchanged. Consequently the
complete diagonal $`f(y)=\langle y|F|y\rangle`$, which is real, has Fourier support in
$`u_0+V`$, or vanishes identically. In the nonempty case,

```math
f(y+v)=(-1)^{u_0\cdot v}f(y)
\qquad\text{for every }v\in V^\perp.                    \qquad(T2)
```

This is an identity of the actual source, not a proposed classical
approximation of its transitions. If $`\|F-Q_s\|\lt1`$, then
$`|f(y)-(-1)^{q_s(y)}|\lt1`$ on every basis column. Thus $`f`$ is nonzero and has
the target's sign everywhere. Equation (T2) forces

```math
q_s(y+v)+q_s(y)=u_0\cdot v
\qquad(v\in V^\perp),                                  \qquad(T3)
```

uniformly in $`y`$. The complete lookup sign has only one nonzero linear
structure: $`v=\boldsymbol1_m`$, with derivative one. This follows from its distinct
Walsh masks $`\boldsymbol1_m+h_t+e_a`$: constant derivative requires $`v_a`$ to be independent
of $`a`$, and then $`Hv=0`$. Since $`m`$ is even, $`v=\boldsymbol1_m`$ has odd pairing with every
one of these masks. Hence (T3) requires $`h\le1`$.

The $`h=0`$ case also fails. Here $`C_2`$ is a monomial Clifford. If $`a\ne0`$,
the conjugated operator still has zero diagonal, giving norm distance at
least one from $`Q_s`$. If $`a=0`$, $`D_1`$ commutes with the first Pauli image and $`F`$
is a diagonal Pauli; $`Q_s`$ is not a Pauli and their norm distance is two.
Therefore a necessary condition for error strictly below one in (T1) is

```math
h=1,\qquad V^\perp=\mathrm{span}\{\boldsymbol1_m\},
\qquad u_0\cdot\boldsymbol1_m=1.                                 \qquad(T4)
```

Every other case has complete operator error at least one. The statement
includes arbitrary diagonal phases, arbitrary bank Clifford choices and
all bank columns, but it retains the source grammar (T1). In particular,
it is not a general restriction on the number of T layers.

<a id="common-affine-two-mask-screen"></a>

### 12.2. Two Fourier masks in one affine bank frame cannot meet the exception

Let $`L`$ be any invertible binary $`m\times m`$ matrix, let $`t`$ be any fixed bank
translation, and let $`U_{L,t}|z\rangle=|Lz+t\rangle`$. In that one common affine frame,
consider

```math
E=U_{L,t}D_2 H_A D_1 H_B D_0U_{L,t}^\dagger,            \qquad(T5)
```

where $`A`$ and $`B`$ are arbitrary sets of bank coordinates. The two masks need
not agree, and the diagonal middle and boundaries are unrestricted.
Every such word has $`\|EPE^\dagger-Q_s\|\ge1`$.

To prove this, the input parity in $`z`$ coordinates is a scalar times
$`Z^\ell`$ with $`\ell=L^{\mathsf T}\boldsymbol1_m`$. The unique nonzero target linear structure is
$`v_*=L^{-1}\boldsymbol1_m`$, still with derivative one. For $`C_2=H_A`$ the space $`V^\perp`$
is the span of the coordinate vectors indexed by $`A`$. The empty $`A`$ case is
the excluded $`h=0`$ case. Otherwise (T3) forces $`A=\{j\}`$ and $`v_*=e_j`$. But then
$`\ell_j=\boldsymbol1_m^{\mathsf T}Le_j=m=0\pmod2`$. The first mask $`H_B`$ therefore creates no
X component on $`j`$. Any X component it creates outside $`j`$ survives the
last $`H_j`$ and gives zero diagonal. If there is no such component, $`D_1`$
commutes with the first Pauli image and the result is a Pauli. Both cases
fail as above. A change of the two masks, arbitrary diagonal phases, and
a shared invertible affine exposure cannot produce the missing lookup.

The fixed affine walls are Clifford and can be emitted with $`O(m^2)`$
CNOTs plus $`O(m)`$ X gates and their actual inverses. This pays those walls,
not an arbitrary middle $`D_1`$. The exclusion remains valid regardless of
that middle's cost. Different affine frames between the mixers, active
extra registers, non-Clifford mixing stages and general source words are
outside (T5).

<a id="rank-one-mixer-lookup-derivative"></a>

### 12.3. The aligned rank-one case is real and exposes the missing phase

Condition (T4) must not be promoted to an impossibility claim. Define the
fixed involutory linear map $`L`$ by

```math
y_0=z_0,\qquad y_j=z_j+z_0\quad(j\ne0).
```

Its first column is $`\boldsymbol1_m`$. Since the target flips under complete bank
complement, there is a Boolean function $`g_s`$ of $`z_1,\ldots,z_{m-1}`$ with

```math
q_s(Lz)=z_0+g_s(z_1,\ldots,z_{m-1}).                   \qquad(T6)
```

Let $`R`$ be the involutory CNOT gather $`z_0\mapsto z_0+\sum_{j\ne0}z_j`$ and set

```math
\begin{aligned}
C_1&=H_{B_0}R,& C_2&=U_LH_{B_0},\\
D_s|z\rangle&=(-1)^{z_0g_s(z_{\ne0})}|z\rangle,\\
E_s&=C_2D_sC_1.
\end{aligned}                                                \qquad(T7)
```

The gather maps $`P`$ to $`Z_0`$, so $`C_1PC_1^\dagger=X_0`$. The exact complete
diagonal derivative is

```math
D_sX_0D_s^\dagger=(-1)^{g_s(z_{\ne0})}X_0.
```

The last $`H_{B_0}`$ and $`U_L`$ give $`E_sPE_s^\dagger=Q_s`$ on every bank column.
There is no initialized bank or omitted inverse. This realizes precisely
the rank-one exception. The Clifford envelope uses two Hadamards and
$`2(m-1)`$ CNOTs in each direction, zero T gates and no additional work.

Its missing operation is $`D_s`$. Formula (T7) merely writes the complete
lookup information into a controlled diagonal phase. Since $`q_s`$ has
Boolean degree $`p+1`$ and an invertible linear substitution preserves degree,
(T6) has $`\deg g_s=p+1`$. The displayed sign phase $`z_0g_s`$ has degree $`p+2`$.
No native sublogarithmic compiler or simultaneous resource ledger for that
middle is supplied. A general complex diagonal with the same derivative
can multiply (T7) by an arbitrary phase independent of $`z_0`$; the required
lookup-dependent derivative is unchanged.

Thus two genuinely different Clifford mixing stages are not excluded.
The complete transition calculation identifies exactly where this route
stores the unpaid operation: in a high-degree diagonal derivative.
Replacing that phase by a charged shallow construction would be new
progress; declaring its truth table, normalized spectrum, or formal
diagonal unitary does not improve the complete-indicator or frame frontier.

The [exact two-mixer tests](../../tests/test_two_mixer_lookup_screen.py) check the $`p=3`$ Walsh
support and every bank translation, a full Clifford-coset source fixture,
the same-frame slope mismatch, and all columns plus actual inverses of
the aligned rank-one construction. Its middle is evaluated as a truth-table
diagonal, explicitly not as a scalable native implementation. These finite
checks support the identities and boundaries; they do not replace the
all-size proof or supply the missing phase compiler.
