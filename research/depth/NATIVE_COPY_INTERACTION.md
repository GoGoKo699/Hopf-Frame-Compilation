# Native copy interactions and their field-orbit screens

[Complete dirty copy interface](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Noncommuting copy grading](DIRTY_COPY_GRADING.md) ·
[Pauli-query copy screen](PAULI_ORBIT_COPY_SCREEN.md) ·
[Field-shear source](FIELD_SHEAR_INTERTWINER.md)

A non-Clifford copy can be a real native component: the reflection
$`THT^\dagger`$ has a two-T, one-band copy interaction with no helpers.
That observation does not supply a compatible all-size field-orbit seed.
This chapter gives an extensive T-count floor for arbitrary copies,
a dual-readout rank test, and exact failures of two concrete product
constructions. The unrestricted one-band field-copy problem remains open.

The [complete one-band feasibility reduction](ONE_BAND_COPY_FEASIBILITY.md)
gives exact affine Walsh support spaces, their spectrum and the remaining
carrier fixed-point equations for arbitrary joint Clifford walls. Output
dual signatures alone are insufficient; the reduced equations still need
a uniform native seed and copy.

Throughout, $`m=2^p`$ with $`p\ge1`$, Y is the physical m-qubit output bank, and
$`R_x=Q_xRQ_x^\dagger`$ is a family of Hermitian involutions on the
full arbitrary carrier A. The fixed interaction C must satisfy

```math
C(R_x\otimes I_Y)C^\dagger=R_x\otimes X_{Y_x}
\quad\text{for every }x.                              \qquad(1)
```

The complete interface owns the echo, its actual inverse and the resource
transfer. No change to the general frame-compiler bound is claimed here.

## 1. Every copy pays an extensive T count

Suppose the binary symplectic entries of the Clifford family Q have
Boolean degree at most d in x. Write
$`M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj`$.
For any fixed native copy C containing t T or T-dagger gates, define its
complete operator error by the maximum discrepancy in (1). Then

```math
\eta\ge
 \sqrt{2-2\sqrt{\min\{1,(M_d(p)+t)/m\}}}.              \qquad(2)
```

Arbitrary Clifford interlayers, carrier/output mixing and arbitrary
helpers are permitted. This is a count bound; it is not a depth bound.

The [projected-axis proof](COPY_AXIS_SPAN.md) is the authoritative home
for this bound and its refinement. It preserves Hilbert--Schmidt mass in
output-label cosets of the absorbed rotation-axis span, including arbitrary
coefficient interference. If that span has $`r_X`$ independent physical-X
directions, the stronger bound replaces t by $`r_X\le t`$. It also
allows arbitrary unitary seeds and arbitrary rotation angles.

For the one-query field shear, d is one. Exact copying needs
$`t\ge m-p-1`$; for $`\eta\lt\sqrt2`$ it needs

```math
t\ge m(1-\eta^2/2)^2-p-1.                            \qquad(4)
```

This linear count is compatible with the discovery budget and gives no
growing native-depth lower bound. The native components and distinct dual
support arguments below retain their own proofs.

## 2. A genuine one-band non-Clifford copy component

Put $`T_P=e^{i\pi/8}e^{-i\pi P/8}`$ for a Hermitian Pauli P.
On one carrier qubit a, consider

```math
A=T_aH_aT_a^\dagger=(X_a+Y_a)/2+Z_a/\sqrt2.
```

Let $`C_H=H_y\,\mathrm{CY}(y;a)\,H_y`$. This Clifford maps
$`X_a\mapsto X_aX_y`$ and $`Z_a\mapsto Z_aX_y`$, and fixes
$`Y_a`$. Hence it copies H. The interaction

```math
C_A=T_aC_HT_a^\dagger
    =T_{Z_a}T_{Z_aX_y}^\dagger C_H                  \qquad(5)
```

copies A. Its two axes commute and are independent. A Clifford D with
$`DZ_aD^\dagger=Z_a`$ and $`DZ_yD^\dagger=Z_aX_y`$ is given
chronologically by $`\mathrm{CX}(a;y),H_y`$. Therefore a literal native
word is $`C_H,D^\dagger,T_a,T_y^\dagger,D`$: two T gates on distinct
wires in one physical band, no initialized inputs and no private helpers.
The retained test checks the full matrix and the actual inverse.

The carrier directions copied by $`C_A`$ form the span of
$`T_aX_aT_a^\dagger`$ and $`Z_a`$. It fixes the complementary span
of I and $`T_aY_aT_a^\dagger`$. Tensoring the component over m
pairs gives 2m T gates in one physical band.

A Clifford cannot copy A itself. Its three nonzero Pauli components
have labels X, Y and Z, with $`Y=X+Z`$. A Clifford copying A would
have to give all three the same nonzero output label X, contradicting
linearity of the label map. Thus (5) genuinely extends the Clifford-copy
class. It does not arrange the required field covariance of its seed.

## 3. Dual support must carry the full translation rank

For the rest of the chapter Q is the additive field-shear representation,
including identity on any helpers:
$`Q_xQ_y=Q_{x+y}`$. For an m-bit table f set

```math
J_f=C^\dagger Z_Y^fC,\qquad K_x=Q_x^\dagger J_fQ_x.
```

Equation (1) implies

```math
R_yK_xR_y^\dagger=(-1)^{f(x+y)}K_x.                 \qquad(6)
```

Let $`r_f`$ be the binary rank of the translates of f. Choose
$`r_f`$ linearly independent translated sign tables. The ordered
products of their corresponding K operators have all
$`2^{r_f}`$ distinct conjugation characters in (6). They are
nonzero unitaries. Any two such products are Hilbert--Schmidt
orthogonal because conjugation by some $`R_y`$ gives them opposite
signs. Commutativity among the $`R_y`$ is not used.

If the Pauli support of $`J_f`$ spans r binary labels, the supports
of all its field conjugates lie in a space of dimension at most
$`(p+1)r`$. Indeed the field symplectic action is
$`S_x=I+\sum_i x_iN_i`$, so the orbit of a support basis v lies in
the span of v and its p images $`N_iv`$. The Pauli algebra of a
label space of dimension h has vector-space dimension $`2^h`$,
whether or not that label space is isotropic. It contains all the
ordered products above. Consequently

```math
r_f\le(p+1)r.                                       \qquad(7)
```

Every odd-weight f has $`r_f=m`$. In the group algebra
$`\mathbb F_2[\mathbb F_2^p]`$, its convolution square is the identity:
cross terms cancel and each group element squares to the identity.
Thus convolution by f, whose columns are its translates, is invertible.

In particular every individual physical output readout
$`f=e_j`$ requires Pauli-support span at least $`m/(p+1)`$.
More generally (6)--(7) hold for any output Pauli P: use its Z-label
on the physical Y bank as f. Output X factors and extra output helpers
do not change those commutation signs.

## 4. Output-only Clifford mixing does not repair the product component

Let $`C_{\rm prod}`$ be a tensor product of the components (5),
with identity on any remaining output helpers. Consider

```math
C=(K_A\otimes K_Y)\,C_{\rm prod}\,(L_A\otimes L_Y),   \qquad(8)
```

where the four walls are arbitrary Cliffords on their indicated
registers. They may mix all carrier helpers among themselves or all
output helpers among themselves. They do not mix carrier and output
registers at these walls.

Use the complete internal output Pauli basis, both X and Z on every
output/helper wire. Pullback through $`C_{\rm prod}`$ gives a Pauli
for an X probe, or a two-term operator with support span at most two
for a Z probe. For a paired Z this is

```math
C_A^\dagger Z_yC_A
 =(T_aY_aT_a^\dagger)Z_y.
```

The pre-walls in (8) preserve support-span rank. Now conjugate the
complete internal basis by $`K_Y`$ to obtain physical output probes.
Their projected Z-labels on the m specified Y wires span
$`\mathbb F_2^m`$. At least one has odd weight; otherwise all would
lie in the even-parity hyperplane. Applying (7) to that probe gives

```math
m\le2(p+1).                                         \qquad(9)
```

Thus this class cannot satisfy (1) for $`p\ge4`$, regardless of
the seed and arbitrary helpers. Both X and Z probes are essential:
the output-Z activation submatrix of a Clifford need not be invertible.

Equation (8) is a restriction on the actual native interaction.
A general Clifford wall that jointly mixes carrier and output is not
covered by this choice of low-support physical probes.

## 5. Coded two-band copies and the fixed field coordinates

There is a broader genuine native construction. Take a full-column-rank
binary matrix A of size $`w\times m`$. Choose a Clifford $`V_0`$
on w carrier wires with

```math
V_0Z_iV_0^\dagger=X^{Ae_i}\quad(0\le i\lt m),
\qquad V=T^{\otimes w}V_0.
```

Such a Clifford exists by extending the independent X-type labels to
a symplectic basis. With the canonical copy from its first m virtual
carrier wires to Y, set

```math
C=V C_{\rm copy}V^\dagger.
```

This is an exact copy for $`R_i=VX_iV^\dagger`$, with 2w T gates
in two physical bands, Clifford walls and no initialized helpers.
Its dual probes are

```math
J_f=(T^{\otimes w}X^{Af}(T^\dagger)^{\otimes w})\otimes Z_Y^f.
                                                               \qquad(10)
```

For nonzero $`Af`$, the carrier support spans
$`\mathrm{wt}(Af)+1`$ labels. A code of large distance therefore
evades a naive low-weight test on these full carrier supports.

This construction does not justify collapsing its two bands for free.
When all rows of A are nonzero, the displayed absorbed axes are the 2w
distinct labels $`Z_j`$ and $`Z_jX_Y^{\mathrm{row}_j(A)}`$.
Their span has dimension w+m. If $`w>m`$ they are dependent, so
they cannot all be the images of distinct native single-wire axes in
one simultaneous Clifford diagonalization. Zero rows give canceling
pairs and should first be removed. This observation does not exclude
a different one-band synthesis.

Now fix the physical field coordinates: the first m carrier wires
are its two field registers; the remaining wires are untouched helpers.
Split $`A=(A_{\rm core};A_{\rm help})`$, where
$`A_{\rm core}`$ is square. The form (10) is assumed in this fixed
basis. An arbitrary outer Clifford mixing core and helpers after the
product T layer is not part of the following class.

Suppose some Hermitian seed R makes this C satisfy (1) with the field
orbit. If $`A_{\rm core}f=0`$, (10) acts only on untouched helpers
and output, so all $`K_x`$ equal $`J_f`$. Equation (6), with
$`y=0`$, forces f to be constant. Hence

```math
\ker A_{\rm core}\subseteq\langle\mathbf1\rangle,
\qquad \mathrm{rank}A_{\rm core}\ge m-1.        \qquad(11)
```

If the core matrix is invertible, its inverse unit columns form a
basis of masks. At least one of these masks is odd, and its core
image has weight one.

If its rank is m-1, the kernel in (11) is exactly the all-one mask.
Since m is even, parity descends to a nonzero linear functional on
the image hyperplane. Every binary hyperplane has a basis of vectors
of weight at most two: for a defining nonzero check h, choose a pivot
j with $`h_j=1`$ and use $`e_i`$ when $`h_i=0`$ or
$`e_i+e_j`$ otherwise, for $`i\ne j`$. Some basis vector has odd
preimage under the descended parity. Thus there is an odd f with
$`\mathrm{wt}(A_{\rm core}f)\le2`$.

Write k for this core weight. The core factor in (10) has support
span at most k+1. Its helper/output factor is one fixed Hermitian
involution F, which is unchanged by Q. Every product of translated
probes belongs to the core orbit Pauli algebra tensored with
$`\mathrm{span}\{I,F\}`$. The dimension argument of Section 3
therefore gives the slightly enlarged bound

```math
m\le(p+1)(k+1)+1.                                   \qquad(12)
```

For an invertible core, k is one and (12) fails for $`p\ge4`$.
In the general case k is at most two and it fails for $`p\ge5`$.
Arbitrarily many untouched helpers cannot repair this particular
coded product-magic class. This exact obstruction does not give a
uniform approximation floor for the class.

## 6. A complete signed-Pauli pilot

A finite instance below the asymptotic exclusion is also inconsistent.
Take $`p=2`$, eight carrier wires, field core
$`(u_0,u_1,v_0,v_1)`$ on the first four, and code columns

```math
A=(255,\ 178,\ 212,\ 232)
```

written as eight-bit integers. The first four rows give an invertible
core, so this is not the trivial helper-only kernel failure.

For each code generator j and each address x, put
$`B_j=T^{\otimes8}X^{Ae_j}(T^\dagger)^{\otimes8}`$ and
$`K_{xj}=Q_x^\dagger B_jQ_x`$. Each is a Clifford Hermitian
involution. A compatible seed must obey the linear coefficient
constraints

```math
K_{xj}RK_{xj}=(-1)^{[x=j]}R.                         \qquad(13)
```

Use Hermitian labels $`\tau(a,b)=i^{\mathrm{wt}(a\mathbin{\&}b)}X^aZ^b`$.
The adjoint action of $`B_j`$ sends

```math
(a,b)\longmapsto(a,b+(a\mathbin{\&}Ae_j)),
\quad\text{sign }(-1)^{\mathrm{wt}(b\mathbin{\&}\neg a
                                      \mathbin{\&}Ae_j)}.
```

Compose this signed permutation with the field query on both sides.
A union-find calculation with parity relations applies all 16
constraints to all $`4^8=65\,536`$ coefficients. It obtains
1,048,576 relations and 2,760 connected components; every component
has a contradictory sign cycle. Thus every coefficient of R must
vanish. No nonzero operator, and hence no Hermitian seed, satisfies
(13) for this pilot.

A second finite pilot uses the same code and applies the same field
query simultaneously to both four-wire carrier blocks. The identical
coefficient calculation again has 65,536 variables and 1,048,576
relations, now with 1,097 connected components and no consistent
component. It also forces R to vanish. This repeated-query example is
outside the one-active-core hypothesis of Section 5; its finite failure
does not extend that all-size theorem to repeated query cores.

The [native fixtures](../../tests/test_native_copy_interaction.py)
retain this exact finite calculation, native checks of its signed label
conventions, the full one-band component and inverse, and a small
complete coded-copy/dual fixture. The analytic all-size statements
are Sections 1--5.

The remaining constructive target is a native interaction whose
carrier/output mixing or carrier phase axes escape these specific
support restrictions, jointly with a compatible seed. The extensive
T-count floor permits such a construction within the discovery count
budget; its complete action and depth still have to be paid.
