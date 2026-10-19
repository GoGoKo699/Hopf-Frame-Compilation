# Representation parameters of easy clean transitions

[Bounded-relation matching](BOUNDED_RELATION_NATIVE_MATCHING.md#8-one-clean-qubit-no-middle-pairs-and-bounded-active-ranks) ·
[Determinant defect](DETERMINANT_DEFECT_DEGREE.md) ·
[Current frontier](../../docs/OPEN_PROBLEM.md)

The single Pauli derivative identity for the clean transition does not
bound the relation dimension of its residual commuting parity-T list.
An actual family with three independent native lists, no middle pairs,
and unbounded endpoint clean ranks already disproves that inference.
Its residual is a product of disjoint CCZ gates. The Boolean degree is
only three, and a zero-clean three-band implementation is explicit.

The determinant-defect codimension can likewise grow for an occupied
identity target, as Section 5 records. These examples distinguish
parameters of a supplied representation from target hardness. They
do not obstruct the unrestricted circuit model or change the
complete-frame resource frontier.

## 1. Three independent clean-active lists

Let n be a positive even integer. Use one clean wire c and n disjoint
logical triples $`x_j=(x_{j,1},x_{j,2},x_{j,3})`$. Choose
$`t_j\in\{1,-1\}`$ with $`\sum_jt_j=0`$. For
$`a\in\mathbb F_2^3`$, define the positive Hermitian Pauli

```math
P_{j,a}=Z_c\prod_{i=1}^3 Z_{j,i}^{a_i}.
```

All these axes commute and anticommute with $`P_A=P_C=X_c`$.
With the literal convention
$`T_P=(I+P)/2+\omega(I-P)/2`$, set

```math
\begin{aligned}
E_A&=\prod_j T_{P_{j,100}}^{t_j}
                  T_{P_{j,010}}^{t_j}T_{P_{j,001}}^{t_j},\\
E_B&=\prod_j T_{P_{j,110}}^{-t_j}T_{P_{j,101}}^{-t_j},\\
E_C&=\prod_j T_{P_{j,011}}^{t_j}T_{P_{j,111}}^{-t_j}.
\end{aligned}                                                    \qquad(1)
```

Each displayed list is physically independent. Projecting its axes
onto the three coordinates of block j gives respectively the
independent sets

```math
\{100,010,001\},\qquad\{110,101\},\qquad\{011,111\}.
```

Different blocks have disjoint projected coordinates, so no cross-block
relation is possible within any list. Their sizes are 3n, 2n and 2n.
They therefore are actual native bands, not dependent lists being
charged as native bands.

For binary c and x, write $`\ell(c,a,x)=c\mathbin\oplus(a\cdot x)`$.
The exact integer identity

```math
\sum_{a\in\mathbb F_2^3}(-1)^{|a|}\ell(c,a,x)
=-4(-1)^c x_1x_2x_3
```

shows that one eight-axis block with coefficients
$`-t_j(-1)^{|a|}`$ has phase $`4x_{j,1}x_{j,2}x_{j,3}`$ modulo
eight. The repeated $`a=000`$ axis is Z_c; its total coefficient is
$`-\sum_jt_j=0`$, so those factors cancel literally. The remaining
seven axes are exactly (1), with C inverted. Consequently

```math
F_0=E_C^\dagger E_BE_A
=I_c\otimes D_{g_n},\qquad
g_n(x)=\sum_{j=1}^n x_{j,1}x_{j,2}x_{j,3}.        \qquad(2)
```

Every scalar in (2) is literal. In particular
$`F_0X_cF_0^\dagger=X_c`$, so the aligning Clifford may be identity
and the residual is F_0 itself.

Use $`K_A=K_C=H_c`$ to make this an actual one-clean circuit. Its two
interior clean images have active ranks 3n and 2n. Every B axis commutes
all A and C axes, and its conjugation by the first clean image is its
literal negative. Thus all B orbits are negative singletons and m=0.
The common active span is the full Z space on the 3n+1 physical wires.

## 2. The complete word has an equally shallow zero-clean baseline

Equation (2) holds on every physical column, including the unused clean
sector. Removing Z_c from every axis in (1) gives a second literal
implementation of the same logical sign: each block still has phase
$`4t_jx_{j,1}x_{j,2}x_{j,3}=4x_{j,1}x_{j,2}x_{j,3}`$ modulo eight.
The projected three lists remain independent. Thus the clean wire can
be omitted entirely for this target.

For explicit H/CNOT coordinate words, the full ledgers are:

| Complete word | T count | Native T-depth | H/CNOT count | Clean work |
| --- | ---: | ---: | ---: | ---: |
| Displayed one-clean word | $`7n`$ | $`\le3`$ | $`20n+2`$ | 1 |
| Actual inverse of that word | $`7n`$ | $`\le3`$ | $`20n+2`$ | 1 |
| Zero-clean baseline | $`7n`$ | $`\le3`$ | $`8n`$ | 0 |
| Actual inverse of baseline | $`7n`$ | $`\le3`$ | $`8n`$ | 0 |

The coordinate words are given in the exact fixture linked below.
Only H and CNOT are used for these Clifford changes, so their reversed
words have the stated counts; every inverse T orientation is reversed.
No arbitrary data wire is initialized. The target acts on the triples,
which are not being claimed as identity-return dirty helpers. Any
additional dirty wires can be left untouched, with arbitrary references.

## 3. Every single parity-T expression has growing relation dimension

The following lower bound allows arbitrary replacement parity axes,
their signs and inverse orientations, a diagonal Clifford correction,
and a retained common scalar. It also permits extra arbitrary variables
on which the target is identity. It is not restricted to the seven
axes used in (1).

Suppose a commuting parity expression for D_gn has the literal form

```math
D_{g_n}=\gamma C\prod_{i=1}^{M}T_{Z^{a_i}}^{\epsilon_i},
\qquad\epsilon_i\text{ odd},                    \qquad(3)
```

where C is diagonal Clifford. Signed axes can be converted to positive
axes with their fixed scalar retained. Identity axes contribute a
scalar; repetitions and even powers can be absorbed into C. Thus one
may use distinct nonzero odd axes in (3). Let A have those a_i as
columns, let d be its binary rank, and let $`r=M-d`$ be its relation
dimension.

The diagonal Clifford exponent has the form

```math
q_C(x)=q_0+2\sum_j u_jx_j+4\sum_{j\lt k}v_{jk}x_jx_k\pmod8.
```

Reduce the phase identity (3) modulo four. A parity expands there as

```math
a\cdot x\pmod2
\equiv\sum_j a_jx_j-2\sum_{j\lt k}a_ja_kx_jx_k\pmod4.
```

Uniqueness of the multilinear coefficients on binary inputs now makes
every first and second odd-axis moment vanish modulo two:

```math
\sum_i(a_i)_j=0,\qquad
\sum_i(a_i)_j(a_i)_k=0\quad(j\ne k).
```

The first identities supply the diagonal entries, so

```math
AA^{\mathsf T}=0.                               \qquad(4)
```

The row space of A is self-orthogonal in the nondegenerate binary dot
product on $`\mathbb F_2^M`$. Hence $`d\le M-d=r`$.

It remains to show that arbitrary new parity choices cannot make d
small. For a translation h in $`\ker A^{\mathsf T}`$, every parity
factor in (3) is unchanged. The derivative of a diagonal Clifford is
a Pauli phase, so $`g_n(x+h)+g_n(x)`$ must be affine. But its
homogeneous quadratic part is

```math
\sum_j\bigl(
h_{j,1}x_{j,2}x_{j,3}
+h_{j,2}x_{j,1}x_{j,3}
+h_{j,3}x_{j,1}x_{j,2}\bigr).                    \qquad(5)
```

These are distinct monomials. Thus every active component of h is
zero. With additional untouched variables, the kernel may use only
those variables, and its codimension is still at least 3n. Therefore

```math
d\ge3n,\qquad r\ge d\ge3n.                     \qquad(6)
```

The common scalar and Clifford signs cannot change either (4) or the
quadratic part (5). An alternative aligning Clifford inside the common
isotropic algebra is exactly a diagonal Clifford correction of this
kind after projection. An arbitrary non-diagonal aligning operation
can leave the commuting-phase class and is not covered by this claim.

## 4. The first unsupported compression step

The derivative identity has already made the residual a commuting
weighted-cubic phase on the occupied logical space. It has not made
that phase a single list of bounded relation dimension. Equations
(1)--(6) disprove such a bound even with three independent original
lists, all their axes clean-active, and complete full-space return.

The obstruction belongs to the proposed representation: the same
family already has the zero-clean three-band baseline in Section 2
and degree three. A successful rank-independent m=0 argument must
retain more of the original multi-band structure, control its
interaction with the silent bands directly, or construct a different
paid representation. It cannot invoke the bounded-relation theorem
after silently assigning a uniform relation budget to this residual.

The [exact native fixture](../../tests/test_clean_transition_relation_obstruction.py)
checks the smallest two-block member, its full physical columns and
actual inverse, both native ledgers, the moment identity, and all
nonzero logical translation directions. The arbitrary-width relation
lower bound is the argument above, not an inference from enumeration.

## 5. Unbounded determinant codimension for literal occupied identity

For every positive integer n, take one clean wire c and n arbitrary
non-clean wires. Use the literal convention
$`T_P=(I+P)/2+\omega(I-P)/2`$, where $`\omega=e^{i\pi/4}`$.
Set

```math
E_A=\prod_{j=1}^nT_{Z_cX_j}^{\dagger},\qquad
E_B=T_{X_c}\prod_{j=1}^nT_{X_j},\qquad
E_C=T_{X_c},\qquad K_A=K_C=I.
```

Each displayed list consists of independent commuting physical Pauli
axes. In particular, the n labels $`Z_cX_j`$ are independent because
their data X coordinates are independent. Thus
$`V_n=E_C^\dagger E_BE_A`$ is an actual three-native-band circuit,
with no synthesized dependent list counted as one physical band.
The two clean X rotations cancel literally, leaving

```math
V_n=\prod_{j=1}^nT_{X_j}T_{Z_cX_j}^{\dagger}
=|0\rangle\langle0|_c\otimes I
 +|1\rangle\langle1|_c\otimes
       \omega^{-n}\prod_{j=1}^nS_{X_j},
\qquad S_P=T_P^2.                              \qquad(7)
```

Indeed $`T_{-P}=\omega T_P^\dagger`$. On the clean-zero sector
the two factors for each j are mutual inverses. On the other sector
their product is $`\omega^{-1}S_{X_j}`$, including its scalar.
Equation (7) is literal occupied identity on every data input and
external reference, with full clean return. The unused block is
$`\exp(-i\pi\sum_jX_j/4)`$.

### 5.1. All determinant hypotheses and the complete defect family

In the actual split above,

```math
J_c=E_AZ_cE_A^\dagger=Z_c,\qquad
S=E_CZ_cE_C^\dagger=T_{X_c}Z_cT_{X_c}^\dagger.
```

The first interior clean image is Pauli; the second is non-Pauli of
intrinsic active rank one. Every middle axis is a singleton for this
clean comparison: J_c negates X_c and fixes every X_j. Thus the clean
pair count is zero. The common Pauli space and its probe space are

```math
K=\langle X_c,Z_c,X_1,\ldots,X_n\rangle,\qquad
N=L_B^\perp\cap\langle X_c,Z_c\rangle^\perp
 =\langle X_1,\ldots,X_n\rangle.
```

The space K has symplectic rank two and center generated by the X_j.
Use their eigenbits u as its intrinsic central coordinates. For every
complete computational-Z label $`z\in\mathbb F_2^n`$, write
$`Z(z)=\prod_jZ_j^{z_j}`$ and $`X(z)=\prod_jX_j^{z_j}`$.
With $`U_C=E_C`$ and $`U_{AB}=U_CV_n`$, direct conjugation in (7)
gives

```math
\begin{gathered}
C_z=U_CZ(z)U_C^\dagger=Z(z),\\
A_z=U_{AB}Z(z)U_{AB}^\dagger,\\
W_z=C_zA_z
 =\Pi_++i^{|z|}X(z)\Pi_-,\\
\Pi_\pm=(I\pm S)/2.
\end{gathered}                                  \qquad(8)
```

For example, the unused one-wire block sends Z_j to $`-Y_j`$,
so multiplying on the left by Z_j gives $`iX_j`$. The independent
wire factors prove (8) for every z, including the literal phase.
Both A_z and C_z conjugate each X_j to $`(-1)^{z_j}X_j`$.
They therefore have identical action on every probe in N. All complete
logical probe characters are zero, not just the characters of selected
logical generators.

The common determinant defect is consequently

```math
\alpha_z(u)=i^{|z|}(-1)^{z\cdot u}.
```

In fact the actual middle determinant and central permutation are
$`\delta(u)=\omega i^{\sum_j u_j}`$ and
$`\sigma_z(u)=u+z`$. Thus

```math
\alpha_z(u)=\frac{\delta(u+z)}{\delta(u)},\qquad
\kappa_z(u)=1.
```

The growth already occurs in a pure common determinant coboundary;
its affine-sign correction family is trivial.

For every nonzero z this is nonconstant on the full intrinsic center.
Thus the constant-defect subgroup and its codimension are

```math
H=\{0\},\qquad r=n.                              \qquad(9)
```

The family satisfies all hypotheses of the
[determinant-codimension degree bound](DETERMINANT_DEFECT_DEGREE.md)
while r grows arbitrarily and its occupied target stays identity.
Its full alpha-family also needs at least n common linear central
coordinates: if it factors through a linear map, every character
$`(-1)^{u_j}`$ factors through that map, forcing its kernel to be zero.
This rules out a width-independent codimension or common-quotient bound
from those hypotheses alone. It does not refute the fixed-r degree
theorem, or establish unbounded r when both interior clean images are
required to be non-Pauli.

### 5.2. Actual native words and inverse costs

Let $`F=\prod_{j=1}^n\mathrm{CX}_{c,j}`$, and let H_data denote
Hadamards on all n non-clean wires. An explicit chronological word is

```math
H_{\rm data};\ F;\ T_{\rm data}^{\dagger};\ F;\ H_c;
\ T_{\rm all};\ H_{\rm data};\ T_c^\dagger;\ H_c.
                                                               \qquad(10)
```

Here a semicolon separates successive operations, and each displayed
T operation is one native band on distinct physical wires. To verify
the first band, $`D=H_{\rm data}F`$ sends Z_j to Z_cX_j, so
$`E_A=DT_{\rm data}^{\dagger}D^\dagger`$. The adjacent Hadamards
in the direct implementations of E_A, E_B and E_C dagger cancel to
give exactly (10).

| Actual operation | T count | Native T-depth | H count | CNOT count | Added work |
| --- | ---: | ---: | ---: | ---: | ---: |
| Word (10) | $`2n+2`$ | $`\le3`$ | $`2n+2`$ | $`2n`$ | 0 |
| Its actual reversed inverse | $`2n+2`$ | $`\le3`$ | $`2n+2`$ | $`2n`$ | 0 |

The displayed schedule has three bands; cancelling the two clean X
rotations gives a two-band implementation of the same full operator.
The determinant parameters above refer to the displayed actual split.
The inverse is the reversed gate word with every T orientation
reversed; it is not a separately chosen identity extension. All n
non-clean wires may be dirty and reference-entangled. Since (7) is
occupied identity, they return exactly under both words.
The [exact fixtures](../../tests/test_clean_transition_parameters.py)
check the actual native words, their inverses, every occupied column,
and the complete defect family for small n. The displayed algebra
proves the arbitrary-n result.
