# Uniform native transport through controlled Paulis

[Six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Hybrid diagonal/core transport](HYBRID_DIAGONAL_CORE_TRANSPORT.md) ·
[Diagonal-source transport](NATIVE_CLIFFORD_TRANSPORT.md)

Every third-level source whose computational control blocks are scalar
Paulis admits width-independent native transport bounds. A signed fixed-
Clifford normalization puts its exact Pauli-preserving domain in canonical
form; no separate domain or fiber-compatibility hypothesis is required.
Polynomial phase arithmetic then avoids interpolating an arbitrary matrix
on the whole target register. The theorem covers a known non-semi-Clifford
source and its tensor powers.

The controlled-Pauli form is a source hypothesis. It does not follow merely
from complete workspace return, arbitrary C3 membership, or a three-native-
band description. No physical synthesis saving or returned-query
construction is asserted.

## 1. Statement and native contract

In fixed analytical input/output Clifford coordinates suppose the source
has scalar-Pauli computational blocks:

```math
F=\sum_{z\in\mathbb F_2^r}|z\rangle\langle z|_R\otimes F_z,
\qquad F_z=e^{i\theta(z)}X^{a(z)}Z^{b(z)}.
```

The Boolean label functions and scalar phase may initially be arbitrary.
The control and target widths are unrestricted. Assume only that F is
third-level. Write $`\omega=e^{i\pi/4}`$ and use the tensor-factorizing
Hermitian section $`\tau(x,y)=i^{x\cdot y}X^xZ^y`$, with integer dot
product in the exponent. All equalities below are literal, including the
phase in $`X^aZ^b`$. Polynomial degrees mean Boolean algebraic degrees
over $`\mathbb F_2`$.

The signed normalization proved below gives fixed Cliffords A,B with

```math
H=AFB=I_S\otimes\sum_w|w\rangle\langle w|_{R'}\otimes H_w,
```

where every $`H_w`$ is still a scalar Pauli and the exact domain is

```math
E_H=\{p:H\tau(p)H^\dagger\text{ is a signed Hermitian Pauli}\}
    =\mathcal P_S\oplus Z_{R'}.
```

The zero label is included. This equality is a conclusion, not an extra
source premise. H remains C3. Below, relabel the normalized controls as
R, denote their string by z, and write their blocks afresh as
$`H_z=e^{i\theta(z)}X^{a(z)}Z^{b(z)}`$. Existing identity spectators
can be included in S.

After the corresponding fixed endpoint changes, let the actual native
six-layer inner family L and outer Clifford comparison K obey

```math
\|K_x-HL_xH^\dagger\|_{\rm op}\le2\eta.
```

[Retain the native ledger](NATIVE_CLIFFORD_TRANSPORT.md#4-native-degrees-including-modulo-eight-carries)
for L and its actual inverse: fixed probe labels have degree two and signs
degree six; a probe label of degree d has output label degree at most
$`d+2`$ and literal conjugation sign degree at most $`3d+6`$. Fixed
Clifford coordinates preserve this ledger and the complete operator error.

**Controlled-Pauli transport theorem.** At full query error
$`\eta\lt\sqrt2/8`$,

```math
\deg S_K\le58,\qquad\deg\kappa_K(P)\le174
```

for every fixed Hermitian probe P, including after restoring the original
fixed output coordinates. The six-layer target degree is at most 58 for
visible outer target shear and at most 174 in general. Neither the bounds
nor the radius depends on control/target width, source coefficients,
Pauli deficiency, or dirty width. No sharpness is claimed.

### Signed normalization of the exact domain

This normalization needs no C3 assumption. For the original source F,
write a binary Pauli label as $`(u,p,v)`$: u,v are its control X,Z labels
and p its target label. Let E be its exact Pauli-preserving domain.
It is a linear subspace with an injective symplectic image map; the
literal signs come from F. Since F fixes every control Z,
$`Z_R\subseteq E`$. Commutation with these fixed operators shows that
every domain image has the same control-X component as its input.

The following construction produces H and its exact domain without
initializing a sector, discarding a scalar, or adding a wire.

#### Ordered controlled-Pauli Cliffords

For any linear map $`A_0:\mathbb F_2^r\to\mathcal P_Q`$, with columns
$`A_j=(\alpha_j,\beta_j)`$, choose the actual ordered product

```math
W_{A_0}=\prod_{j=1}^r
\left(|0\rangle\langle0|_j\otimes I
      +|1\rangle\langle1|_j\otimes\tau(A_j)\right).
```

Each factor is Clifford: a target Clifford takes the Hermitian axis to
a single signed Z, reducing its controlled gate to CZ with a possible
control Z. The product fixes control Z and has scalar-Pauli blocks. Its
binary action, and that of its inverse, both have form

```math
(u,p,v)\longmapsto(u,p+A_0u,v')
```

with a linear control-Z tail. Conjugating a pure control translation
contributes $`A_0u`$ to the target label; conjugating an existing target
Pauli changes only its control-Z tail. These assertions also cover a
zero column, whose controlled gate is the identity.

The columns need not commute. For the displayed algebraic order the
literal block is

```math
(W_{A_0})_z=i^{\psi(z)}X^{\sum_j\alpha_jz_j}
                              Z^{\sum_j\beta_jz_j},
```

```math
\psi(z)=\sum_j(\alpha_j\cdot\beta_j)z_j
       +2\sum_{j\lt k}(\beta_j\cdot\alpha_k)z_jz_k\pmod4.
```

The X,Z exponents are binary sums; the displayed dots are integer dots.
This retains the order-dependent scalar. Replacing the block by
$`\tau(A_0z)`$ can change the operator. The source class is closed under
left/right multiplication by these actual Cliffords, left multiplication
by diagonal control Cliffords, and conjugation by control-only invertible
linear coordinates. These are literal block identities.

#### Straighten mixed input translations

Let V be the control-X projection of E, with $`h=\dim V`$. Choose a
linear section $`u\mapsto e(u)\in E`$. Its target component $`p_D(u)`$
is linear. Extend it to all control translations and construct the
ordered Clifford $`W_D`$ as above.

For $`u\in V`$, $`W_DX^uW_D^\dagger`$ has control-X component u
and target component $`p_D(u)`$. Its label differs from $`e(u)`$ only
by a control-Z label, which is in E. Hence

```math
F_1=FW_D,\qquad E_{F_1}=S_{W_D}^{-1}E
```

has every pure $`X^u`$, $`u\in V`$, in its exact domain. Its control-X
projection is still exactly V because $`W_D`$ preserves that component.
It still fixes $`Z_R`$ literally.

The selected section need not commute: the commuting lifts supplied
by $`W_DX^uW_D^\dagger`$ differ by available control-Z tails. Signs
are irrelevant only to domain membership; the actual $`W_D`$ signs
remain in $`F_1`$.

#### Match output translations with all signs

For $`u\in V`$, the Pauli $`F_1X^uF_1^\dagger`$ has control-X
component u and a linear target component $`p_I(u)`$. Extend the latter
and construct $`W_I`$ in the same way. Put $`F_2=W_I^\dagger F_1`$.
The displayed binary action cancels the target component even though
it was already present in the image. Thus

```math
F_2X^uF_2^\dagger=(-1)^{s(u)}\tau_R(u,c(u)),\qquad u\in V,
```

where $`c:V\to\mathbb F_2^r`$ is linear. Left multiplication by a
Clifford leaves the exact input domain unchanged. Commuting input
translations have commuting images, giving

```math
u\cdot c(v)=v\cdot c(u),\qquad u,v\in V.
```

A symmetric linear map $`C:\mathbb F_2^r\to\mathbb F_2^r`$ extends
the entire prescribed map c. To see this, choose a complement of V and
extend the bilinear form $`(w,u)\mapsto w\cdot c(u)`$. Its V-by-V
block is symmetric. Its cross block is prescribed by pairing c with
the complement; the opposite block is its transpose. Choose an arbitrary
symmetric complement-by-complement block. Thus the extension preserves
the values of c paired with all R, not only with V.

A diagonal control Clifford realizes the shear
$`(u,v)\mapsto(u,v+Cu)`$: diagonal matrix entries use S gates and
off-diagonal entries use CZ gates in the original computational basis.
For $`u\in V`$ its X image has the same label as the image under
$`F_2`$. The ratio of their signs is a character of V, because both
conjugations respect the multiplication of $`\{X^u:u\in V\}`$ with
the same Pauli cocycle. A control Z Pauli realizes any such character.
Include that correction in an actual diagonal Clifford G. It satisfies

```math
GX^uG^\dagger=F_2X^uF_2^\dagger,\qquad u\in V,
```

literally. Consequently $`F_3=G^\dagger W_I^\dagger F W_D`$
commutes all $`X^u`$, $`u\in V`$, and all control Z.

Choose a control CNOT/SWAP Clifford $`C_R`$ taking V to the translation
span of the first h control wires $`S_0`$. Then

```math
F_4=C_RF_3C_R^\dagger=I_{S_0}\otimes F_5.
```

Commutation with the full Pauli algebra on $`S_0`$ proves this exact
factorization on arbitrary entangled inputs and references. No wire is
initialized or physically discarded. On remaining controls $`R_0`$
the domain has zero control-X projection; subtract its available
control-Z tails to obtain the exact equalities

```math
E_{F_5}=E_Q\oplus Z_{R_0},\qquad
E_{F_4}=\mathcal P_{S_0}\oplus E_{F_5}.
```

No new control translation can appear after stripping the exact identity
factor, because the projection before stripping was exactly V.

#### Match the target domain and promote its radical

For $`p\in E_Q`$, blockwise Pauli conjugation gives

```math
F_5\tau(p)F_5^\dagger=(-1)^{\kappa_p(z)}\tau(p),\qquad
\kappa_p(z)=\langle P_5(z),p\rangle.
```

The target label is unchanged. Exact-domain membership makes this sign
affine in the remaining control string:

```math
\kappa_p(z)=s(p)+\sum_jf_j(p)z_j.
```

Evaluation at zero and each control basis vector shows that s and all
$`f_j`$ are linear functionals of p on $`E_Q`$, including when that
subspace is degenerate. Each extends to the full target space, whose
nondegenerate pairing represents it by a Pauli label. Choose representing
labels $`q_0,q_j`$ and the actual ordered Clifford

```math
W_0=\tau(q_0)\prod_j
\left(|0\rangle\langle0|_j\otimes I
      +|1\rangle\langle1|_j\otimes\tau(q_j)\right).
```

Its block label is $`q_0+\sum_jq_jz_j`$, so its literal conjugation
matches $`F_5`$ on every Hermitian Pauli in $`E_Q`$. The phases of the
ordered product are retained. Therefore $`F_6=W_0^\dagger F_5`$
commutes $`E_Q`$ and $`Z_{R_0}`$ literally, with the same exact domain.

Let $`2s_1`$ be the symplectic rank of $`E_Q`$ and c its radical
dimension. A target-only Clifford $`C_Q`$ takes it to
$`\mathcal P_{S_1}\oplus Z_C`$, with $`|S_1|=s_1`$ and $`|C|=c`$.
Conjugation preserves scalar-Pauli blocks including its actual signs.
Commutation with $`\mathcal P_{S_1}`$ forces identity on $`S_1`$;
commutation with $`Z_C`$ forbids X labels on C. Each block thus has form

```math
I_{S_1}\otimes Z_C^{b_C(z)}
       \otimes\widetilde\lambda_z\tau(\widetilde P(z)).
```

Regard C as additional controls, absorbing its full action into the
scalar $`(-1)^{b_C(z)\cdot c_0}\widetilde\lambda_z`$ at the joint
sector $`(z,c_0)`$. This regroups every sector; it selects no eigenvalue.
A fixed wire permutation $`\Pi`$ puts the registers in the desired
order. Then, exactly,

```math
E_H=\mathcal P_{S_0\cup S_1}\oplus Z_{R_0\cup C}.
```

Wire counts are $`|S|=h+s_1`$, $`|R'|=r-h+c`$, and
$`|Q'|=q-s_1-c`$, whose sum is unchanged. The domain dimension checks:
$`\dim E=r+h+2s_1+c=2|S|+|R'|`$.

#### Both endpoint coordinates and literal error

Extend $`W_0,C_Q`$ by identity on $`S_0`$. The complete construction is

```math
\begin{aligned}
H&=AFB,\\
A&=\Pi(I_{S_0}\otimes C_QW_0^\dagger)
       C_RG^\dagger W_I^\dagger,\\
B&=W_DC_R^\dagger(I_{S_0}\otimes C_Q^\dagger)\Pi^\dagger.
\end{aligned}
```

Each factor is its actual selected operator, including any scalar, so
$`F=A^\dagger HB^\dagger`$ literally. Given an endpoint
$`K_x\approx FL_xF^\dagger`$, set

```math
K'_x=AK_xA^\dagger,\qquad L'_x=B^\dagger L_xB.
```

Then

```math
\|K'_x-HL'_xH^\dagger\|_{\rm op}
=\|K_x-FL_xF^\dagger\|_{\rm op}.
```

The actual inverse is $`(L'_x)^\dagger=B^\dagger L_x^\dagger B`$.
The fixed input and output coordinates generally differ; using the same
coordinate conjugation on both sides would not preserve this identity.
Initially supplied source endpoint Cliffords can be absorbed into A,B.
Existing identity spectator wires may be included in S.

The normalization is analytical and uses only fixed Cliffords on either
side of F, so it preserves C3. The following phase bounds are derived
afresh for H; no polynomial degree is propagated through these possibly
wide coordinate changes.

## 2. C3 derives the polynomial phase bounds

Third-level membership forces bounded-degree control data. The raw
scalar phase may have a weighted quartic term; it need not itself be
diagonal C3.

First, conjugation of target $`X_i`$ and $`Z_i`$ gives respectively

```math
H X_i H^\dagger=(-1)^{b_i(z)}X_i,
\qquad H Z_i H^\dagger=(-1)^{a_i(z)}Z_i.
```

The left sides are Clifford. Hence the Boolean scalar phases on controls
are diagonal Cliffords, so every $`a_i,b_i`$ has degree at most two.

For a control coordinate j put $`\delta_j a(z)=a(z+e_j)+a(z)`$, and
similarly for b. The Clifford $`J_j=HX_jH^\dagger`$ is a monomial
Hermitian involution. Its permutation on a basis input $`(z,t)`$ is

```math
(z,t)\longmapsto(z+e_j,t+\delta_j a(z)),
```

and its literal phase at that input is

```math
e^{i(\theta(z+e_j)-\theta(z))}
(-1)^{\delta_j b(z)\cdot(t+a(z))}.
```

All additions inside the binary Pauli labels are binary; the dot in this
sign is parity. A monomial Clifford has affine permutation and, up to a
scalar $`\lambda_j`$, a phase $`i^{\ell(z,t)}(-1)^{q(z,t)}`$, where
$`\ell`$ is integer linear and q is Boolean quadratic. Because $`J_j`$
is a literal Hermitian involution, the product of its two nonzero entries
along a permutation cycle is one. Hence $`\lambda_j^2\in\mu_4`$, so
$`\lambda_j\in\mu_8`$, where $`\mu_k`$ denotes the k-th roots of unity.
Every nonzero entry of $`J_j`$ is therefore an eighth root. The adjacent
ratio $`e^{i(\theta(z+e_j)-\theta(z))}`$ also lies in $`\mu_8`$.
Connectedness of the Boolean cube gives

```math
e^{i\theta(z)}=e^{i\theta(0)}\omega^{p(z)}
```

for a uniquely normalized modulo-eight function p. The fixed initial
scalar cancels from all transport formulas and remains literal in H.

At target input $`t=0`$, the phase of $`J_j`$ is a weighted-quadratic
modulo-eight function $`r_j(z)`$, including its constant term. Consequently

```math
p(z+e_j)-p(z)
=r_j(z)-4\,\delta_j b(z)\cdot a(z)\pmod8.
```

The final term has Boolean degree at most three, with coefficient four:
$`\delta_j b`$ is affine and a is quadratic. Write the unique multilinear
integer modulo-eight expansion $`p(z)=\sum_S c_S\prod_{k\in S}z_k`$.
Restrict the last displayed identity to $`z_j=0`$. For $`j\in S`$,
its coefficient on $`S\setminus\{j\}`$ is exactly $`c_S`$.
Coefficients with $`|S|\ge2`$ are even, those with $`|S|\ge3`$ are
divisible by four, and those with $`|S|\ge5`$ vanish. Therefore

```math
p(z)=c+\sum_i \alpha_i z_i+2\sum_{i\lt j}\beta_{ij}z_i z_j
 +4\sum_{i\lt j\lt k}\gamma_{ijk}z_i z_j z_k
 +4\sum_{i\lt j\lt k\lt l}\delta_{ijkl}z_i z_j z_k z_l\pmod8.
```

The extra quartic terms contribute only to the high bit, at degree $`4h_0`$
when z coordinates have degree $`h_0`$. The phase p has bit-degree
bounds $`(h_0,2h_0,4h_0)`$, which are the phase bounds used below. No
source-specific phase gap or restriction on control width is needed.

## 3. Exact compatibility and fiber formula

The [C3 rounding lemma](NATIVE_CLIFFORD_TRANSPORT.md#1-a-robust-exactness-lemma-for-every-third-level-source)
makes the endpoint projectively exact, because its error is at most
$`2\eta\lt\sqrt2/4`$. It requires $`H\in\mathcal C_3`$, without an
inverse-C3 premise. The canonical domain proved in section 1 makes exactness force L
to normalize $`E_H`$ and its radical $`Z_R`$, and H commutes $`E_H`$
literally, as in the [hybrid proof](HYBRID_DIAGONAL_CORE_TRANSPORT.md#5-exact-endpoint-and-domain-compatibility).
Therefore

```math
J_x=L_x^\dagger H L_x H^\dagger
```

is an actual Clifford up to its retained scalar and commutes $`E_H`$.
The actual K and LJ differ only by a scalar, irrelevant to their Pauli
conjugation action after this exact reduction.

L permutes control sectors affinely by $`z\mapsto\phi_x(z)`$; its
fibers factor as Q and S Cliffords. Thus J is identity on S and commutes
$`Z_R`$. These properties now follow from the normalized exact domain,
without a separate compatibility premise.
At every fixed control input z,

```math
J_{x,z}=M_{x,z}H_{\phi_x(z)}M_{x,z}^\dagger H_z^\dagger,
\qquad M_{x,z}=L_{Q,x,z}^\dagger.
```

The S factors and fiber scalars cancel literally. The local inverse
Clifford M is read from the actual inverse L at fixed input sector z.
The control string $`\phi_x(z)`$ has coordinate degree at most six,
from $`L_x^\dagger Z_jL_x`$. No inverse-affine substitution is required.

## 4. Polynomial phase arithmetic independent of width

For Hermitian labels $`u=(x,y)`$ and $`v=(s,t)`$, literal multiplication gives

```math
\tau(u)\tau(v)=i^{e(u,v)}\tau(u+v),
```

```math
\begin{aligned}
e(u,v)&=\sum_i(-x_it_i+s_iy_i)\\
&\quad+2\sum_i(x_iy_it_i+s_iy_it_i+x_is_iy_i+x_is_it_i)\pmod4.
\end{aligned}
```

This follows by expanding XOR as $`x+s-2xs`$, retaining integer dot products.
For two moving labels of degree at most d, the low bit of e has degree
at most $`2d`$ and its high bit at most $`4d`$: the odd-weight terms
have degree $`2d`$, the even-weight terms degree $`3d`$, and the high-bit carry is a second
symmetric polynomial of the odd-weight terms. If one label is fixed,
the respective bounds improve to d and $`2d`$. Arbitrarily many coordinates
are summed without increasing these degree bounds.

Convert $`H_z`$ to the Hermitian section, retaining $`\lambda=e^{i\theta(0)}`$:

```math
H_z=\lambda\omega^{h(z)}\tau(P(z)),\qquad
P(z)=(a(z),b(z)),\qquad h(z)=p(z)-2a(z)\cdot b(z)\pmod8.
```

If the coordinates of z have degree at most $`h_0`$, the three bits of $`h(z)`$
have degrees at most

```math
(h_0,\ 4h_0,\ 8h_0).
```

Section 2 gives p bit degrees $`h_0,2h_0,4h_0`$, including its quartic
terms. Each $`a_i b_i`$ has degree $`4h_0`$, and multiplication of their sum by two puts its parity in the weight-two
bit and its pairwise carry in the weight-four bit. The subtraction carry
preserves the displayed bounds. This is why one must not silently treat
the scalar of the Hermitian Pauli section as weighted cubic.

For two modulo-eight expressions whose bit degrees are bounded by
$`(c_0,c_1,c_2)`$, addition and subtraction preserve those bounds whenever
$`c_1\ge2c_0`$ and $`c_2\ge2c_1`$. This is the binary carry rule; fixed constants are
allowed. Adding twice a modulo-four expression with bit degrees $`(d_0,d_1)`$
gives bit bounds

```math
(c_0,\ \max(c_1,d_0),\ \max(c_2,d_1,c_1+d_0)).
```

Any lower-bit carry involving $`c_0`$ is already bounded by the stated terms.

## 5. Every fixed J block is a low-degree literal Pauli

At a fixed sector z, $`P(\phi_x(z))`$ has label degree at most 12,
while the phase $`h(\phi_x(z))`$ has three bit degrees at most $`(6,24,48)`$. Conjugating
that moving Hermitian probe by M uses the actual native inverse:

```math
M_{x,z}\tau(P(\phi_x(z)))M_{x,z}^\dagger
=(-1)^{\sigma_z(x)}\tau(T_z(x)),
\quad \deg T_z\le14,\quad\deg\sigma_z\le42.
```

This local action is obtained by restricting
$`L_x^\dagger\tau(P(\phi_x(z)))L_x`$ to the fixed input control z.
Its control-Z tail adds at most degree 14 to the sign, already below 42. No target-width tableau interpolation is
used.

Since $`H_z`$ is fixed, multiplication by its adjoint has only one moving
Pauli label. The phase rule of section 4 therefore gives

```math
J_{x,z}=\omega^{\nu_z(x)}\tau(t_z(x)),
\quad \deg t_z\le14,
\quad\deg(\nu_{z,0},\nu_{z,1},\nu_{z,2})\le(6,24,48).
```

The fixed scalar $`\lambda`$ cancels. The remaining phase consists of
$`h(\phi_x(z))-h(z)`$, four times $`\sigma_z`$, and twice the
modulo-four multiplication phase. The latter has bit degrees 14,28;
their carry into h has degree at most $`24+14=38`$. Thus 48 bounds the
high bit, including every literal factor.

## 6. Fixed-probe recovery and the sign carry

Omit the identity S factor. Fix $`\tau(p_Q)\otimes\tau(u,v)_R`$,
with u,v fixed. Since J is Clifford and commutes $`Z_R`$, its conjugation image is

```math
(-1)^{\kappa(x)}\tau(t(x))\otimes\tau(u,v'(x)).
```

At any fixed control string z let

```math
A_z=J_{x,z+u}\tau(p_Q)J_{x,z}^\dagger.
```

The two J-block labels have degree 14, so $`A_z`$ has target label
degree 14. Its modulo-eight scalar phase has bit degrees at most
$`(6,28,56)`$: the J phases have bounds $`(6,24,48)`$, while the product
of the two moving Paulis contributes twice a modulo-four phase with
bit degrees 28,56. The fixed
middle probe does not increase these bounds.

The exact full Pauli identity is

```math
A_z=\alpha(x)(-1)^{(v'(x)-v)\cdot z}\tau(t(x)),
\quad\alpha(x)=(-1)^{\kappa(x)}i^{u\cdot(v'(x)-v)}.
```

Thus $`A_z`$ has scalar in $`\mu_4`$: its weight-one bit is identically
zero. Write its scalar as $`i^{\mu_z}`$. The two bits of $`\mu_z`$
have degrees at most 28 and 56. At $`z=0`$ the target label t has degree
14, and $`\alpha=i^{\mu_0}`$.
For each control coordinate j,

```math
A_{e_j}A_0^\dagger=(-1)^{v'_j-v_j}I.
```

The exact equal target labels cancel. Therefore $`v'_j-v_j`$ is the
high bit of $`\mu_{e_j}-\mu_0`$, of degree at most 56. At most the four
fixed blocks $`0,u,e_j,e_j+u`$ are used for this coordinate; no fixed four blocks
are claimed to reconstruct the entire wide output simultaneously.

Finally the literal sign obeys

```math
2\kappa=\mu_0-u\cdot(v'-v)\pmod4.
```

The integer sum $`u\cdot v'`$ has low bit degree 56 and high bit degree
112. The extra subtraction borrow has degree at most $`28+56=84`$. Hence

```math
\deg S_J\le56,\qquad\deg\kappa_J(P)\le112.
```

This retains the Hermitian-section carry rather than summing signs of
an unbounded collection of moving generators.

## 7. Native composition and scope

$`K=LJ`$ projectively. Apply the actual native moving-probe ledger to a
label of degree 56 and incoming sign degree 112. Output labels have
degree 58 and signs at most $`3\cdot56+6=174`$. Restoring fixed output
Clifford coordinates adds only a cubic sign in a degree-58 label, also 174.
The existing visible-shear and general fixed-probe target extractors
add degree-two and degree-six baselines respectively, so the bounds 58/174 remain.

The scalar-Pauli source form is an additional hypothesis. It does not
follow from complete workspace return, from arbitrary C3 membership, or
from a three-native-band description alone. This theorem does not make
its source semi-Clifford, nor construct a returned query. Its precision
starts at the complete operator endpoint comparison and retains all
logical columns, arbitrary dirty/reference inputs, literal phases,
leakage and full return. Every physical source, gate, control, actual
inverse, helper and simultaneous workspace remains charged. All block
and coordinate operations are analytical.

## 8. Known seven-qubit example with literal source accounting

The following source is exactly Clifford-conjugate to the
Gottesman–Mochon example reported by Salman Beigi and Peter W. Shor,
*C3, Semi-Clifford and Generalized Semi-Clifford Operations*,
Quantum Information and Computation **10** (2010), 41–59,
[Section 1.1](https://arxiv.org/pdf/0810.5108). It is an explicit scope
fixture; no novelty is claimed for the source.

Use four controls $`z_1,\ldots,z_4`$ and three targets
$`t_1,t_2,t_3`$. Define

```math
f(z)=(z_1z_2,z_1z_3,z_1z_4),\qquad
g(z)=(z_3z_4,z_2z_4,z_2z_3),\qquad F_z=X^{f(z)}Z^{g(z)}.
```

These blocks have the required source form. When both exponents on one
target equal one, the literal factor is $`XZ=-iY`$; its phase is retained.
The full basis action is

```math
F|z,t\rangle=(-1)^{g(z)\cdot t}|z,t+f(z)\rangle.
```

### C3 membership

F fixes all control Z operators. Its target generator images are

```math
F X_{t_i}F^\dagger=(-1)^{g_i(z)}X_{t_i},\qquad
F Z_{t_i}F^\dagger=(-1)^{f_i(z)}Z_{t_i}.
```

These are a target Pauli times a control CZ, hence Clifford. Write
$`\delta_j`$ for the Boolean derivative in $`z_j`$. A control-X image
has block

```math
F_{z+e_j}F_z^\dagger
 =(-1)^{\delta_j g(z)\cdot f(z)}
       X^{\delta_j f(z)}Z^{\delta_j g(z)}.
```

For every j, $`\delta_j g\cdot f=0`$ identically. The derivative vectors
are

| j | $`\delta_j f`$ | $`\delta_j g`$ |
|---|---|---|
| 1 | $`(z_2,z_3,z_4)`$ | $`(0,0,0)`$ |
| 2 | $`(z_1,0,0)`$ | $`(0,z_4,z_3)`$ |
| 3 | $`(0,z_1,0)`$ | $`(z_4,0,z_2)`$ |
| 4 | $`(0,0,z_1)`$ | $`(z_3,z_2,0)`$ |

Within each row, the X and Z target supports are disjoint. Thus the
control-X image is $`X_{z_j}`$ times controlled X and controlled Z gates,
with literal scalar +1. Every generator image is Clifford, and their
products are Clifford, proving $`F\in\mathcal C_3`$.

### Its exact Pauli-preserving domain

Let an input Pauli have control-X translation u and fixed target label p.
Its output target label contains $`P(z+u)+P(z)`$, where $`P=(f,g)`$.
The six coordinates of P are all six distinct quadratic monomials
$`z_i z_j`$, $`i\lt j`$. For that difference to be independent of z,
every $`u_i z_j+u_j z_i`$ must vanish, which forces $`u=0`$.

With $`u=0`$, the variable phase for a target label $`(x,v)`$ is

```math
g(z)\cdot x+f(z)\cdot v.
```

The six quadratic monomials are independent, so this is affine only when
$`x=v=0`$. The exact domain is consequently

```math
E_F=\mathrm{span}\{Z_{z_1},Z_{z_2},Z_{z_3},Z_{z_4}\}.
```

Its largest isotropic subspace has dimension four on seven qubits.
Therefore F has no Pauli-preserving Lagrangian and is not semi-Clifford.
For a tensor product of sources, a tensor-product Pauli has a Pauli
image exactly when each factor does: its Pauli support is the Cartesian
product of the factor supports. Disjoint tensor powers thus have exact
domain equal to all their control Z axes, with no additional directions.
They satisfy the theorem at unbounded target width while remaining
non-semi-Clifford. Every copy remains charged as a physical source.

### Literal square and failure of inverse-C3

Since $`f\cdot g`$ equals three copies of $`z_1z_2z_3z_4`$,

```math
F^2=(-1)^{z_1z_2z_3z_4}I_Q.
```

This is a quartic control phase, not a global scalar. In the inverse
conjugation of $`X_{z_1}`$, the derivative phase is
$`\delta_1 f\cdot g=z_2z_3z_4`$. Hence

```math
F^\dagger X_{z_1}F
 =X_{z_1}\left(\prod_{j=2}^4\mathrm{CX}(z_j;t_{j-1})\right)
       \mathrm{CCZ}(z_2,z_3,z_4).
```

All displayed factors commute. The first factors are Clifford, whereas
CCZ is not; their product is therefore not Clifford. In particular
$`F^\dagger\notin\mathcal C_3`$. The transport theorem needs no
inverse-C3 premise.

Target Hadamards also illustrate the quartic phase permitted by section 2:
they change $`X^fZ^g`$ into
$`(-1)^{f\cdot g}X^gZ^f`$. This Clifford-conjugate source remains C3,
although its raw scalar phase is $`p(z)=4z_1z_2z_3z_4`$ modulo eight.

### Identification with the recorded Gottesman–Mochon example

[Beigi–Shor, Section 1.1](https://arxiv.org/pdf/0810.5108) describes U as
three controlled swaps on $`(R,A_i,B_i)`$, and V as four CCZ gates with
phase

```math
A_1A_2A_3+A_1B_2B_3+B_1A_2B_3+B_1B_2A_3.
```

Use the invertible binary Clifford coordinates

```math
z_1=R,\quad z_2=A_1+B_1,\quad z_3=A_2+B_2,\quad z_4=A_3+B_3,
\qquad t_1=A_1,\quad t_2=B_2,\quad t_3=B_3.
```

The controlled swaps become $`t\mapsto t+f(z)`$. Substitution into the
four-term cubic phase cancels all terms except

```math
t_1z_3z_4+t_2z_2z_4+t_3z_2z_3=g(z)\cdot t.
```

Thus UV is literally Clifford-conjugate to F, with no residual cubic or
global phase. The attribution and exact coordinate equivalence are
separate from the new transport bound.

### Fully paid native realization

Apply, chronologically, the three CCZ gates on

```math
(z_3,z_4,t_1),\qquad(z_2,z_4,t_2),\qquad(z_2,z_3,t_3),
```

then the three CCX gates with control/target triples

```math
(z_1,z_2;t_1),\qquad(z_1,z_3;t_2),\qquad(z_1,z_4;t_3).
```

This is exactly $`X^fZ^g`$. The
[retained literal CCZ parity word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger) uses
seven T/T-dagger gates, ten CNOTs and at most five physical T layers,
based on the identity

```math
4abc=a+b+c-(a\mathbin\oplus b)-(a\mathbin\oplus c)
     -(b\mathbin\oplus c)+(a\mathbin\oplus b\mathbin\oplus c)\pmod8.
```

CCX is that word conjugated by two target Hadamards. Serial execution
therefore gives the complete source ledger

```math
\begin{aligned}
T&=42,\qquad G=66=60\ \mathrm{CNOT}+6\ H,\\
D_T&\le30,\qquad \text{data width}=7,\\
\text{extra helpers}&=\text{clean helpers}=0.
\end{aligned}
```

The error and leakage are zero. Every physical column has the displayed
sign. The actual reversed native inverse has the same ledger, and
$`F^\dagger F=I`$ literally. F acts nontrivially on its seven data wires;
this source fixture is not itself a returned-work query. The parity word
uses control coordinates as native T targets, so it supplies no strict
native-read-only-address certificate. No three-native-band realization
or resource improvement is claimed.

## 9. Finite verification

The [focused tests](../../tests/test_controlled_pauli_transport.py) check
literal phase arithmetic, the source's exact domain and its inverse,
the Clifford coordinate equivalence to the recorded example, and the
fully charged native word. Two exact endpoints built from the actual
source and Pauli inner operators exercise four-block recovery, including
imaginary scalar phases and a carry that changes the recovered sign.
Tensor-product checks retain the arbitrary-width scope of the stated
family. A five-wire normalization fixture checks all physical columns
and the full Pauli domain with mixed, noncoordinate control translations,
noncommuting section target labels, essential literal sign corrections,
and a degenerate target-domain direction promoted to a control. Its
different left and right endpoint coordinates preserve an exact Clifford
endpoint and the actual inverse, including a nonreal scalar. These
bounded fixtures support the explicit mechanisms; the degree bounds and
source-uniform radius follow from the proof above.

