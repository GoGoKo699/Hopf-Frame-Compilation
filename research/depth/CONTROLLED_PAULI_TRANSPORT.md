# Uniform native transport through controlled Paulis

[Six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Hybrid diagonal/core transport](HYBRID_DIAGONAL_CORE_TRANSPORT.md) ·
[Diagonal-source transport](NATIVE_CLIFFORD_TRANSPORT.md)

A third-level source whose computational control blocks are scalar
Paulis admits width-independent native transport bounds under an explicit
compatibility condition on its exact Pauli-preserving domain. The proof
uses polynomial phase arithmetic instead of interpolating an arbitrary
matrix on the whole target register. It applies to a known non-semi-
Clifford source and its tensor powers.

The controlled-Pauli form and domain condition are source hypotheses.
Neither follows merely from complete workspace return or an arbitrary
three-native-band description. No physical source-synthesis saving or
returned-query construction is asserted.

## 1. Statement and native contract

In fixed analytical input/output Clifford coordinates suppose

```math
H=I_S\otimes\sum_{z\in\mathbb F_2^r}|z\rangle\langle z|_R
       \otimes H_z,
\qquad H_z=e^{i\theta(z)}X^{a(z)}Z^{b(z)}.
```

The Boolean vectors a,b and scalar function $`\theta`$ may initially be
arbitrary. Write $`\omega=e^{i\pi/4}`$. The control width r and target width q are unrestricted.
Assume H is third-level and its exact signed Pauli-preserving domain is

```math
E_H=\{p:H\tau(p)H^\dagger\text{ is a signed Hermitian Pauli}\}
    =\mathcal P_S\oplus Z_R.
```

The zero label belongs to this domain. All equalities are literal,
including the phase in $`X^aZ^b`$. Use the tensor-factorizing Hermitian
section $`\tau(x,y)=i^{x\cdot y}X^xZ^y`$, where the exponent's dot
product is over the integers. All polynomial degrees below are Boolean
algebraic degrees over $`\mathbb F_2`$.

The exact-domain premise can instead be replaced by the explicit
compatibility condition that the particular inner family L preserves
the control algebra and has Clifford fibers factoring on Q and S.
The proof uses only that compatibility once the displayed source form
is supplied. Neither compatibility nor the domain equality is inferred
from a generic controlled-source representation.

Let the actual native six-layer inner family L and outer Clifford
comparison K obey the complete operator endpoint bound

```math
\|K_x-HL_xH^\dagger\|_{\rm op}\le2\eta.
```

[Retain the native ledger](NATIVE_CLIFFORD_TRANSPORT.md#4-native-degrees-including-modulo-eight-carries)
for L and its actual inverse: fixed probe labels
have degree two and signs degree six; a probe label of degree d has output
label degree at most $`d+2`$ and literal conjugation sign degree at most $`3d+6`$.
Fixed Clifford coordinates preserve the ledger.

**Controlled-Pauli transport theorem.** At full query error $`\eta\lt\sqrt2/8`$,

```math
\deg S_K\le58,\qquad\deg\kappa_K(P)\le174
```

for every fixed Hermitian probe P. Thus the six-layer target degree is at
most 58 for visible outer target shear and at most 174 in general.
Neither the bounds nor the radius depends on q,r, source coefficients,
or dirty width. The constants are conservative; no sharpness is claimed.

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
inverse-C3 premise. Under the domain hypothesis, exactness makes L
normalize $`E_H`$ and its radical $`Z_R`$, and H commutes $`E_H`$
literally, as in the [hybrid proof](HYBRID_DIAGONAL_CORE_TRANSPORT.md#5-exact-endpoint-and-domain-compatibility).
Therefore

```math
J_x=L_x^\dagger H L_x H^\dagger
```

is an actual Clifford up to its retained scalar and, under the
exact-domain premise, commutes $`E_H`$.
The actual K and LJ differ only by a scalar, irrelevant to their Pauli
conjugation action after this exact reduction.

L permutes control sectors affinely by $`z\mapsto\phi_x(z)`$; its
fibers factor as Q and S Cliffords. If the alternative explicit fiber
compatibility premise is used, this step is supplied directly. In either
case J is identity on S and commutes $`Z_R`$; the alternative premise
does not identify the radical of a possibly larger $`E_H`$.
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

The source form and exact domain are additional hypotheses. They do not
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
family. These bounded fixtures support the explicit mechanisms; the
degree bounds and source-uniform radius follow from the proof above.

