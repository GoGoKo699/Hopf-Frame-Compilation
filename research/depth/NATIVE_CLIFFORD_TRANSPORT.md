# Clifford transport through a third-level source

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[One-band scope](ONE_BAND_AXIS_MATCHING.md)

This note separates robust exactness from polynomial degree. Any fixed
third-level source has a width-independent rounding radius when both
comparison endpoints are Clifford. For a source diagonal in fixed Clifford
coordinates, a signed native transport argument also bounds the output
tableau degree. The source Clifford being transported may be nonmonomial;
none of the intermediate native band passages is assumed Clifford.

The resulting six-layer theorem requires the specified three-band prefix
to have the additional form below. It does not resolve the general
six-layer question or improve a complete-frame compiler.

The later [native source certificates](NATIVE_BAND_SOURCE_STRUCTURE.md)
identify additional physical factorizations satisfying this source form.
The [bounded-deficiency theorem](BOUNDED_DEFICIENCY_TRANSPORT.md) instead
handles a small Pauli-domain codimension without third-level membership;
its hypothesis is independent of the diagonal-source hypothesis here.

## 1. A robust exactness lemma for every third-level source

Write $`\tau(p)`$ for the Hermitian Pauli section, and call F third-level
when $`F\tau(p)F^\dagger`$ is Clifford for every p. Global phases are
retained in operator estimates; $`\doteq`$ denotes equality up to scalar
only in subsequent algebraic comparisons.

**Rounding lemma.** Let F be any third-level unitary and C,K be Cliffords,
possibly with arbitrary scalar phases. If

```math
d_{\rm pr}(FCF^\dagger,K)\lt\frac{\sqrt2}{4},
\qquad
d_{\rm pr}(U,V)=\inf_{|\alpha|=1}\|U-\alpha V\|_{\rm op},
```

then $`FCF^\dagger\doteq K`$ exactly.

Put $`U=FCF^\dagger`$ and $`A_p=F\tau(p)F^\dagger`$. For every p,
both

```math
UA_pU^\dagger=F[C\tau(p)C^\dagger]F^\dagger,
\qquad KA_pK^\dagger
```

are Hermitian Cliffords. Their distance is at most
$`2d_{\rm pr}(U,K)\lt1/\sqrt2`$. The retained
[Clifford separation lemma](READ_ONLY_THREE_LAYER_RIGIDITY.md#81-peeling-the-outer-band-and-quantizing-a-clifford-cocycle)
makes them projectively equal. Hermiticity restricts the scalar to plus
or minus one, and the distance is less than two, excluding the minus sign.
They are therefore literally equal. Thus $`K^\dagger U`$ commutes with
every $`A_p`$. Those operators span the full matrix algebra, so their
commutant consists only of scalars. This proves the lemma.

The proof does not assume $`F^\dagger`$ third-level, or closure of the
third level under multiplication. It uses the actual source adjoint in U.
It yields exactness, not by itself a polynomial formula for the transported
Clifford tableau.

## 2. The diagonal source and its Pauli directions

Assume first that F is diagonal with

```math
F|z\rangle=\omega^{q(z)}|z\rangle,\qquad \omega=e^{i\pi/4},
```

up to one fixed scalar, where

```math
q(z)=c+\sum_i a_i z_i
       +2\sum_{i\lt j}b_{ij}z_i z_j
       +4\sum_{i\lt j\lt k}c_{ijk}z_i z_j z_k\pmod8.       \qquad(1)
```

Here $`a_i\in\mathbb Z_8`$, $`b_{ij}\in\mathbb Z_4`$ and
$`c_{ijk}\in\mathbb F_2`$. Translation of (1) shows directly that F
is third-level: its conjugation of a Pauli is a Pauli times a diagonal
Clifford. This includes any diagonal phase polynomial formed from physical
T and T-dagger phases on affine work parities.

Define the fixed binary space

```math
\mathcal R=\{a:F X^a F^\dagger\text{ is Pauli up to scalar}\}.
```

It is a linear space, because conjugation preserves multiplication and F
commutes with every Z Pauli. There exists a **fixed diagonal Clifford** G
whose action agrees literally with F on every $`X^a`$, $`a\in\mathcal R`$.
Here is the phase-sensitive extension argument. The output label has form
$`(a,c(a))`$, with c linear on $`\mathcal R`$. Commutation gives

```math
a\cdot c(b)=b\cdot c(a)\quad(a,b\in\mathcal R).
```

Consequently c extends to a symmetric binary matrix on the full work
space: extend a basis of $`\mathcal R`$ to a basis of the space, prescribe
the corresponding rows and columns by symmetry, and choose the remaining
symmetric block arbitrarily. A diagonal Clifford built from S and CZ
has precisely this shear. The residual signs on a basis of
$`\mathcal R`$ can be matched by a fixed Z Pauli. Multiplicativity then
matches the literal action on every $`X^a`$ in the space, including its
Hermitian signs. No unsigned replacement is made.

Put $`H=G^\dagger F`$. It is diagonal, commutes literally with all
$`X^a`$, $`a\in\mathcal R`$, and its phase still has form (1).
It is therefore constant on cosets of $`\mathcal R`$. Choose fixed
independent vectors $`r_1,\ldots,r_k`$ spanning $`\mathcal R^\perp`$
and a fixed linear section of the quotient. Restricting the phase of H
to that section gives a polynomial $`\bar q`$ of form (1), and

```math
H|z\rangle
=\omega^{\bar q(r_1\cdot z,\ldots,r_k\cdot z)}|z\rangle. \qquad(2)
```

Substitution of binary linear parities preserves form (1) modulo eight;
this follows by retaining terms of degree at most three in their integer
parity expansion. All choices in this paragraph are fixed with the source
and independent of the address.

## 3. Arbitrary Clifford inputs reduce to a diagonal correction

Let $`L_x`$ be any Clifford family such that
$`F L_xF^\dagger`$ is Clifford up to scalar. For every input Z generator,
the Pauli $`L_x Z_jL_x^\dagger`$ must have its X label in
$`\mathcal R`$, because its conjugation by F is Pauli. It follows that
for every $`r\in\mathcal R^\perp`$,

```math
L_x^\dagger Z^r L_x=(-1)^{\beta_r(x)}Z^{a_r(x)}.           \qquad(3)
```

For clarity, if the symplectic matrix of L has blocks

```math
S_L=\begin{pmatrix}A&B\\C&D\end{pmatrix},
```

the first statement is $`\mathrm{im}B\subseteq\mathcal R`$.
The X label of $`L_x^\dagger Z^rL_x`$ is $`B^T r=0`$, proving (3).
This permits nonmonomial L; only the indicated quotient observables are
forced diagonal.

Functional calculus applied to (2) and (3) gives the exact identity

```math
\begin{aligned}
\Delta_x&=(L_x^\dagger H L_x)H^\dagger,\\
\Delta_x|z\rangle&=\omega^{p_x(z)}|z\rangle,\\
p_x(z)&=\bar q\bigl(\beta_{r_j}(x)+a_{r_j}(x)\cdot z\bigr)_{j=1}^k
          -\bar q\bigl(r_j\cdot z\bigr)_{j=1}^k\pmod8,\\
F L_xF^\dagger&=G L_x\Delta_x G^\dagger.                 \qquad(4)
\end{aligned}
```

The correction $`\Delta_x`$ is diagonal Clifford, because the last
operator and G,L are Clifford. Its inverse order and all phases in (4)
are literal. Under approximate Clifford endpoints, section 1 first makes
the endpoint projectively exact, so the same argument applies; no
intermediate native band is rounded or assumed Clifford.

## 4. Native degrees, including modulo-eight carries

Specialize now to the exact inner family L of the
[six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md#2-a-clifford-valued-native-pullback).
It is a product of two one-band Pauli-image families and a linear-label
Pauli. The retained
[native moving-probe ledger](READ_ONLY_SIX_LAYER_REDUCTION.md#4-native-tableau-degree-of-the-eight-factor-product)
says that a probe label of degree d has output label degree at most
$`d+2`$ and literal sign degree at most $`3d+6`$. The same holds for
the actual inverse, by reversing the two one-band factors and retaining
their orientations. Fixed Clifford conjugations preserve these bounds:
their signs are cubic in the probe labels and the final label degree is
at most $`d+2`$.

In particular, the labels $`a_{r_j}(x)`$ in (3) have degree at most two,
and the signs $`\beta_{r_j}(x)`$ have degree at most six. For each fixed
work string z, every argument of $`\bar q`$ in (4) has degree at most
six in x.

The bit of weight four in a modulo-eight phase of form (1), evaluated
on Boolean functions of degree at most h, has degree at most $`4h`$,
not necessarily $`3h`$. To see this, split the integer sum into terms
of weights one, two and four. Their Boolean monomials have degrees at
most h, $`2h`$ and $`3h`$, respectively. The weight-four bit of a sum
of weight-one bits is its fourth elementary symmetric polynomial; the
weight-two bit is its second. Carries between the first two weight
classes consequently also have degree at most $`4h`$. The weight-four
class has degree at most $`3h`$. This proves the bound, including arbitrary
width and all modulo-eight carries. The lower bits have degree at most h
and $`2h`$. Differences of finitely many such phase evaluations have
the same bounds by the same weighted-sum argument.

Apply this with h=6 to $`p_x(z)`$ at fixed z. Its weight-four bit has
degree at most 24. Since $`\Delta_x`$ is diagonal Clifford, write

```math
p_x(z)-p_x(0)
=2\sum_i\ell_i(x)z_i+4\sum_{i\lt j}m_{ij}(x)z_i z_j\pmod8.
```

Its diagonal shear entries and generator signs are obtained from

```math
\ell_i(x)=\frac{p_x(e_i)-p_x(0)}2\pmod4,\qquad
m_{ij}(x)=\frac{p_x(e_i+e_j)-p_x(e_i)-p_x(e_j)+p_x(0)}4\pmod2.
                                                               \qquad(5)
```

The promised divisibility follows from diagonal Clifford membership.
All shear entries have degree at most 24, as do the generator signs.
For a fixed Hermitian probe $`\tau(a,b)`$, put
$`c_x=B_{\Delta_x}a`$. Its image has label $`(a,b+c_x)`$, of
degree at most 24, and literal sign determined by

```math
4\kappa_{\Delta_x}(a,b)
=p_x(a)-p_x(0)+2[a\cdot b-a\cdot(b+c_x)]\pmod8.          \qquad(6)
```

The dot products here are integer-valued and $`b+c_x`$ is a binary
sum. The last term's weight-four bit is the weight-two bit of a sum of
degree-at-most-24 bits, so its degree is at most 48. Thus the fixed-probe
sign of $`\Delta_x`$ has degree at most 48. This conservative bound
does not discard the Hermitian-section carries.

Finally conjugate a fixed probe through the last line of (4). The right
fixed $`G^\dagger`$ leaves it fixed. The correction has output label
degree 24 and sign degree at most 48. L then has output label degree
at most 26 and sign degree at most $`3\cdot24+6=78`$. The left fixed
G contributes a cubic sign on a degree-26 label, also at most 78. Hence

```math
\deg_x S_{F L_xF^\dagger}\le26,\qquad
\deg_x\kappa_{F L_xF^\dagger}(q)\le78
\quad\text{for every fixed Hermitian probe }q.             \qquad(7)
```

Both bounds are independent of work width. The same bounds hold when the
fixed source has form $`\mathcal F=C_{\rm out}F C_{\rm in}`$ with
fixed Cliffords on either side: conjugate L into the input coordinates
and the output back. The native input ledger is unchanged, and the last
fixed Clifford adds at most another cubic sign on a degree-26 label.

## 5. Consequence for the actual six-layer query

Use the full native contract and actual definitions of $`\mathcal F,L_x,K_x`$
in the [six-layer proof](READ_ONLY_SIX_LAYER_REDUCTION.md). Suppose its
specified three-band prefix admits the additional factorization

```math
\mathcal F=C_{\rm out}F C_{\rm in},
```

where the Cliffords are fixed and F has form (1). This is an additional
source hypothesis, not a consequence of complete workspace return.
Such a source is third-level. The complete-circuit comparison estimate is

```math
\|\mathcal F L_x\mathcal F^\dagger-K_x\|\le2\eta.
```

At full operator error $`\eta\lt\sqrt2/8`$, section 1 makes each comparison
projectively exact. Equation (7) therefore bounds the actual K tableaus.
No Pauli-valued comparison or monomial L hypothesis is needed.

If the outer target shear is visible, the retained
[fixed linear extractor](READ_ONLY_FIVE_LAYER_REDUCTION.md#3-visible-outer-shear-a-translation-invariant-space-of-functions)
expresses $`g(x)=f(x)+f(0)`$ as a fixed linear functional of $`S_K(x)`$
plus a quadratic baseline. It follows that $`\deg f\le26`$.
If the outer target shear vanishes, a fixed anticommuting Pauli probe
has sign $`p_q(x)+g(x)`$, where the native baseline satisfies
$`\deg p_q\le6`$. Equation (7) then gives $`\deg f\le78`$.

Thus this specified-source six-layer subcase has degree at most 78,
improving to 26 for visible outer shear, below $`\sqrt2/8`$.
Sharpness or attainability is not asserted. A three-band native prefix
need not admit the displayed source factorization; its general degree
transport remains unresolved. Even the general third-level rounding
lemma supplies no unrestricted tableau-degree bound by itself.

## 6. Operator and resource meaning

Every approximate estimate starts from the complete query norm, including
every logical column, literal target phase, arbitrary dirty inputs and
external references, leakage and complete workspace return. Projective
comparisons are used only after those estimates and leave every Pauli
conjugation sign literal. No initialized state or selected column is used.

The fixed coordinate choices, diagonal Clifford extension, quotient phase
and correction in this proof are analytical identities. They add no
physical operation or helper and are not offered as free native sources.
The original circuit retains all six physical T slots, all T and T-dagger
gates, every elementary H, S and CNOT, every control and actual inverse,
and its simultaneous arbitrary-work allocation. T-count, elementary
Clifford count and finite dirty width remain unrestricted in this necessary
condition; no complete-circuit upper resource bound changes.

For example, monomial Clifford intervals between the prefix's physical
bands allow the T phases to be collected into (1), up to a fixed monomial
Clifford. This is a source identification, with all original gates still
charged. It does not assume that conjugating a nonmonomial L through each
of those bands separately stays Clifford. More general nonmonomial
prefixes require a separate proof of the stated factorization.

The concrete remaining question is a width-independent degree map for
the actual three-band prefix beyond this fixed diagonal source class,
or an actual native counterexample to such transport. The exactness lemma
isolates that degree question for any third-level prefix; outside the
third level, even that rounding conclusion needs additional justification.

The six
[native transport checks](../../tests/test_native_clifford_transport.py)
cover the signed diagonal Clifford extension, a native nonmonomial input
with a nontrivial correction, the actual inverse quotient observables,
modulo-eight carries, Hermitian Pauli signs, and finite rounding examples.
They also retain a non-Clifford intermediate with an exact Clifford
endpoint and a two-band source already outside the third level. These bounded
checks support the identities; the proofs above supply the arbitrary-width
bounds.

## 7. Third-level membership alone does not bound tableau degree

The following permutation family is due to
[He, Robitaille and Tan, Theorem 42](https://doi.org/10.4230/LIPIcs.TQC.2026.2).
The retained use here is a self-contained signed-transport stress test with
an explicit native ledger. It is not claimed as a new third-level
permutation construction. The example separates the exactness lemma from
a degree theorem: exact third-level sources can send a linear Pauli family
to Clifford families of unbounded tableau and sign degree.

No three-T-band realization of these sources is supplied. Thus the example
does not contradict the specified-source theorem above or settle the
general six-layer query problem.

### A quadratic coordinate map with a higher-degree inverse

Fix $`r\ge2`$. Let $`\mathcal A_r`$ be the binary algebra with one basis
element $`e_S`$ for every nonempty subset of $`[r]`$, with multiplication

```math
e_Ie_J=
\begin{cases}
e_{I\cup J},&I\cap J=\varnothing,\\
0,&I\cap J\ne\varnothing.
\end{cases}
```

This associative commutative algebra has dimension $`w=2^r-1`$, and every
element a satisfies $`a^2=0`$. Consequently $`a\circ b=a+b+ab`$ is an
elementary abelian group operation: it is multiplication of $`1+a`$ and
$`1+b`$ in the algebra with an identity adjoined, and every element is its
own inverse.

Choose any fixed order on the basis labels and put

```math
\nu(a)=\sum_{I\lt J}a_Ia_J e_Ie_J,\qquad k(a)=a+\nu(a).
```

The quadratic map satisfies

```math
\nu(a+b)=\nu(a)+\nu(b)+ab,\qquad \nu(ab)=0.               \qquad(8)
```

For the second identity, expand $`ab=\sum_{I,J}a_Ib_Je_Ie_J`$.
A nonzero term of its quadratic image has an unordered pair
$`\{(I,J),(K,L)\}`$ of ordered pairs, with all four sets pairwise
disjoint. Pair it with $`\{(I,L),(K,J)\}`$. The terms have equal
coefficient and basis product, and are distinct because the sets are
nonempty and disjoint. They cancel in characteristic two. Terms with an
overlap were already zero. Equation (8) and $`a^2=b^2=0`$ now give

```math
\begin{aligned}
k(a\circ b)
 &=a+b+ab+\nu(a)+\nu(b)+\nu(ab)+ab+(a+b)ab\\
 &=k(a)+k(b).
\end{aligned}                                             \qquad(9)
```

In coordinates,

```math
k(a)_S=a_S+
 \sum_{\substack{\{I,J\}:I,J\ne\varnothing\\I\sqcup J=S}}a_Ia_J.
                                                               \qquad(10)
```

Only strict subsets of S occur in the added terms, so k is a triangular
quadratic bijection. Its inverse f is

```math
f(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1.     \qquad(11)
```

Indeed k sends each $`e_S`$ to itself, and (9) sends the product in (11)
to the ordinary sum of those basis elements. Every nonzero product in (11)
uses disjoint nonempty subsets, so f has Boolean degree at most r. On
singleton-supported inputs,

```math
f\!\left(\sum_{i=1}^r x_i e_{\{i\}}\right)
 =\sum_{\varnothing\ne S\subseteq[r]}
       \left(\prod_{i\in S}x_i\right)e_S,                \qquad(12)
```

so its degree is exactly r. The actual inverse of k must therefore not be
replaced by k itself when $`r\ge3`$.

### Literal third-level transport and growing degrees

Represent the w coefficients by physical work qubits and define
$`F_r|c\rangle=|f(c)\rangle`$. For every fixed b, equation (9) gives

```math
F_rX^bF_r^\dagger|a\rangle
 =|a\circ f(b)\rangle=|a+f(b)+af(b)\rangle.              \qquad(13)
```

This is an invertible affine work permutation and hence a Clifford.
For every fixed t,

```math
F_rZ^tF_r^\dagger|a\rangle
 =(-1)^{t\cdot k(a)}|a\rangle.                           \qquad(14)
```

The exponent is Boolean quadratic, so (14) is a product of Z and CZ
gates. Equations (13)--(14) prove that $`F_r`$ is third-level. They hold
literally on every work column. No initialized input or inverse-closure
property of the Clifford hierarchy is used.

Now take r read-only address bits and

```math
b(x)=\sum_{i=1}^r x_i e_{\{i\}},\qquad
L_x=X^{b(x)},\qquad K_x=F_rL_xF_r^\dagger.
```

Write $`a(x)=f(b(x))`$. The affine map (13) has linear part
$`B_x=I+M_{a(x)}`$, where $`M_a(y)=ya`$. Since $`a^2=0`$,
$`B_x^{-1}=B_x`$ and $`B_xa=a`$, and its symplectic matrix is

```math
S_{K_x}=\begin{pmatrix}B_x&0\\0&B_x^{-T}\end{pmatrix}.
```

For nonempty T and S, an off-diagonal entry of B is
$`a_{S\setminus T}(x)`$ when $`T\subsetneq S`$, and otherwise zero.
All entries have degree at most $`r-1`$, while

```math
(B_x)_{[r],\{1\}}=\prod_{i=2}^r x_i                      \qquad(15)
```

has that degree exactly. The fixed probe $`Z_{[r]}`$ has literal image

```math
K_xZ_{[r]}K_x^\dagger
 =(-1)^{a_{[r]}(x)}Z^{B_x^Te_{[r]}},\qquad
a_{[r]}(x)=\prod_{i=1}^r x_i.                            \qquad(16)
```

Its fixed-probe sign degree is exactly r. Thus exact Clifford endpoints
and third-level source membership alone admit unbounded symplectic and
signed degree growth. The diagonal source hypothesis in the degree theorem
provides an additional restriction beyond the rounding lemma.

### Charging the native source and its actual inverse

Equation (10) gives an in-place construction of k. Process target sets S
in decreasing cardinality. For each unordered disjoint nonempty pair I,J
with union S, apply one Toffoli with controls I,J and target S. Controls
are strict subsets of S and still contain their original inputs. Gates
with the same target cardinality cannot modify one another's controls.
The word therefore computes k on every work input; its actual reversed
inverse computes f and hence $`F_r`$.

The number of Toffolis is

```math
M_r=\frac{3^r-2^{r+1}+1}{2}.                              \qquad(17)
```

There are $`3^r`$ ordered disjoint pairs including empty sets; subtract
both empty-set cases, restore the doubly subtracted empty pair, and divide
by two. Use the retained exact
[seven-T CCZ word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger)
with two target Hadamards for each Toffoli. This costs seven T/T-dagger
gates, twelve elementary Cliffords and at most five T layers, with no
helper. Its Cliffords are H and CNOT, so the actual reverse has the same
elementary count. Sequential emission supplies the explicit upper ledger
for each of $`F_r`$ and its actual inverse:

```math
T=7M_r,\qquad G=12M_r,\qquad D_T\le5M_r,\qquad
W=2^r-1\text{ arbitrary work qubits},\qquad a_{\rm clean}=0.
```

There are no additional helpers, approximation errors or leakage. A
displayed implementation of $`K_x`$ pays both sources and the r
address-controlled CNOTs, giving

```math
T=14M_r,\qquad G=24M_r+r,\qquad D_T\le10M_r
```

on the same simultaneous arbitrary-work allocation. These are direct
upper bounds, not optimality claims. Every source, control and actual
inverse is charged. Neither $`F_r`$ nor $`K_x`$ is advertised as a returned
XOR query: $`K_x`$ generally acts nontrivially on work. A source followed
by its actual inverse returns all work, including external references,
but supplies no query implementation.

No construction here places $`F_r`$ in the specified three physical T
bands of the general six-layer reduction. That missing native condition
is why this example diagnoses a hierarchy-only degree argument without
refuting native three-band transport. The
[bounded checks](../../tests/test_c3_transport_scope.py) verify the algebra,
actual inverse, growing degrees and a literal native source on every
small work column.

The later [balanced native construction](BALANCED_ALGEBRA_PERMUTATION.md)
improves this source's depth to $`O(\log r)`$ with $`O(3^r)`$
T-count, elementary Clifford count and additional arbitrary returned
dirty work. It pays the in-place multiplier and compares with a
cardinality-batched $`O(r)`$ baseline. Its separate returned degree-r
query treats every coefficient as address data. It still supplies
neither a three-band source nor a bounded-depth unbounded-degree query;
the zero-extra-helper serial ledger above keeps its distinct allocation.
