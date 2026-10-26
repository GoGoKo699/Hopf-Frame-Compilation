# Complete Clifford symmetries of lookup

[Two-outer-conjugation search](CLIFFORD_INVOLUTION_SOURCE.md#two-outer-conjugations) ·
[Clifford-interleaved lookup](CLIFFORD_INTERLEAVED_LOOKUP.md) ·
[Complete Clifford correlation](CLIFFORD_LOOKUP_FIXED_SPACE.md)

There is an exact normal form for every Clifford R whose lookup conjugate
QRQ is also Clifford. The characterization imposes no prior restriction
on address, coefficient-bank or helper mixing. Every such R factors into a Clifford fixing a specified
computational algebra and a monomial Clifford satisfying an explicit
Boolean quadratic condition. The first factor has no remaining restriction;
all nontrivial compatibility is in that affine permutation condition.
For p at least three, every such Clifford normalizes the address-Z algebra.
It may affinely relabel addresses, but cannot mix their linear forms with
coefficient or helper coordinates. The two-address exception admits genuine
address/high-coefficient swaps.

This characterizes an intersection of two Clifford groups. A separate
uniform rounding lemma gives distance at least $`1/\sqrt2`$ from the
Clifford group whenever QRQ is not Clifford. It makes the symmetry screen
robust at complete probe error below $`1/(2\sqrt2)`$. Neither result
classifies all three-block sources or gives a growing native-depth lower
bound: a general intermediate observable need not be Clifford on either
side of lookup.

<a id="full-clifford-symmetry-normal-form"></a>

## 1. The complete symmetry group and its normal form

Use the coefficient-bank Clifford coordinates, with p address bits,
$`m=2^p`$ arbitrary coefficient bits and any number of arbitrary helpers h:

```math
Q|x,c,h\rangle=(-1)^{f(x,c)}|x,c,h\rangle,
\qquad f(x,c)=\sum_{S\subseteq[p]}c_S\prod_{i\in S}x_i.
```

Assume $`p\ge2`$. Call x together with the coefficients $`c_S`$ for
$`|S|\ge2`$ the hard coordinates z. The other coordinates w consist of
$`c_\varnothing`$, all p singleton coefficients and every helper.
The names describe this computational partition only; they do not assert
individual qubit hardness or initialized status. Let $`\mathcal S`$ be
the group of Cliffords fixing every $`Z_{z_j}`$ literally. Such Cliffords
may mix all w wires, with z-dependent Clifford fiber phases, while keeping
every z value.

Define

```math
\mathcal G_Q=\{R\in\mathcal C_2:QRQ\in\mathcal C_2\}.
```

**Normal-form theorem.** Exactly the elements of $`\mathcal G_Q`$ have
factorizations

```math
R=S M,\qquad S\in\mathcal S,\qquad
M|v\rangle=\lambda i^{q(v)}|\pi(v)\rangle,
\qquad \deg\bigl(f(\pi(v))+f(v)\bigr)\le2.             \qquad(1)
```

Here $`v=(x,c,h)`$, M is a monomial Clifford, pi is its invertible affine
binary permutation, and q is its Clifford phase polynomial modulo four.
The degree is Boolean ANF degree on the entire physical input. There is
no separate restriction on M's literal Clifford phase. Every S in the
stated group is permitted. The factorization need not be unique.

For p zero or one, Q itself is Clifford, so every Clifford R belongs to
$`\mathcal G_Q`$ and the hard/soft restriction is unnecessary.

## 2. A quadratic criterion on the affine Choi support

The following observation applies to any diagonal Boolean sign Q, not
only lookup. Let R be Clifford and let L be the affine support of its
normalized Choi state, with coordinates $`(v_{\rm out},v_{\rm in})`$.
Then

```math
QRQ\text{ is Clifford}
\quad\Longleftrightarrow\quad
\bigl[f(v_{\rm out})+f(v_{\rm in})\bigr]\big|_L
\text{ is Boolean quadratic}.                         \qquad(2)
```

Both directions retain literal phases. The Choi state of a Clifford has
constant-modulus amplitudes on an affine support and a modulo-four
quadratic phase with even cross coefficients. Multiplying R on both sides
by Q changes neither support nor modulus; its amplitude ratio is exactly
$`(-1)^{f(v_{\rm out})+f(v_{\rm in})}`$. If both operators are Clifford,
the ratio of their stabilizer phases is modulo-four quadratic. Because
that ratio takes only the values plus and minus one, all coefficients
are even after fixing its value at one affine origin. Boolean Möbius
inversion therefore gives a Boolean quadratic exponent. Conversely,
multiplying a stabilizer state's affine-support phases by a Boolean
quadratic sign produces another stabilizer state. QRQ remains a unitary,
so its maximally entangled stabilizer Choi state is that of a Clifford.

The statement concerns exact equality, not small entrywise error. In
particular, no uniform rounding radius follows merely by comparing the
individual Choi amplitudes as the register width grows. Section 6 instead
proves a uniform radius using complete Pauli-image gaps.

## 3. Computational support can move only in affine-derivative directions

Let W be the common direction space of the computational column supports
of R. A vector a belongs to W precisely when $`(a,0)`$ is a direction
of L. Taking a difference of the quadratic in (2) along that direction
gives

```math
f(v_{\rm out}+a)+f(v_{\rm out})
\quad\text{affine on }L.
```

The output projection of L is surjective, since every row of a unitary
is nonzero. A surjective affine projection has an affine section. Thus
this same function is affine on the entire output space, proving

```math
W\subseteq\mathcal L_f,
\qquad
\mathcal L_f=\{a:\deg(f(v+a)+f(v))\le1\}.              \qquad(3)
```

This derivative space is the retained
[complete Pauli-subspace calculation](INDICATOR_PAULI_NULLITY.md#2-the-complete-indicators-exact-pauli-subspace),
now used to constrain an entire Clifford rather than an individual Pauli.
For coefficient lookup and $`p\ge2`$, it is exactly

```math
\mathcal L_f=
\{(a_x,a_c,a_h):a_x=0,\ (a_c)_S=0\text{ for }|S|\ge2\}.
                                                               \qquad(4)
```

Indeed, if some address displacement $`(a_x)_j`$ is one, the coefficient
of the independent $`c_{[p]}`$ in the difference is
$`\prod_i(x_i+(a_x)_i)+\prod_i x_i`$. Its degree-(p-1) monomial
omitting j has coefficient one. Including $`c_{[p]}`$ gives degree p,
which is at least two. Coefficient and helper translations cannot cancel
this independent-bank term. Once $`a_x=0`$, the difference is simply
$`\sum_S(a_c)_S x^S`$, affine precisely under (4).

Thus the column support of every R in $`\mathcal G_Q`$ has no hard
coordinate displacement. In particular, R has no address Hadamard
direction; among bank coordinates it has at most p+1 independent such
directions. Arbitrarily many helper directions remain allowed. This is
a statement about the full Clifford R's support, not a prohibition on
individual Hadamards inside a canceling native decomposition.

In the ordinary addressed bank Y, the same condition reads

```math
\mathcal L_{Y_x}=
\{(a_x,b_Y,a_h):a_x=0,\ b_Y(x)=b_0+\sum_i b_i x_i\}.
                                                               \qquad(5)
```

The fixed bank Möbius Clifford maps these affine truth-table masks to the
empty and singleton coefficient directions in (4). All helper values stay
arbitrary in either coordinate system.

## 4. Factoring the unrestricted Clifford and proving sufficiency

For every hard output Z, (4) says its computational value is constant on
each column support of R. Since its pullback is a Pauli, it follows that

```math
R^\dagger Z_{z_j}R
\quad\text{is a signed computational Z Pauli}.
```

These pullbacks have independent binary labels. Extend their affine
linear forms to an invertible affine map pi and choose the corresponding
monomial Clifford M, with signs fixed so that

```math
M^\dagger Z_{z_j}M=R^\dagger Z_{z_j}R
\quad\text{for every hard coordinate }j.
```

Set $`S=RM^\dagger`$. Then $`S^\dagger Z_{z_j}S=Z_{z_j}`$
literally, so $`S\in\mathcal S`$. Binary elimination constructs M using
$`O(n^2)`$ CNOT and X gates on n physical wires; X is paid as $`HS^2H`$
in the native alphabet. This factorization is a description of the given
Clifford, not an uncharged nonlinear operation.

Every member of $`\mathcal S`$ is automatically in $`\mathcal G_Q`$.
To see this, split

```math
Q=Q_{\rm low}Q_{\rm high},\qquad
Q_{\rm low}=(-1)^{c_\varnothing+\sum_i c_i x_i},\qquad
Q_{\rm high}=(-1)^{\sum_{|S|\ge2}c_S x^S}.
```

The low factor is Clifford. The high factor depends only on z and hence
commutes every S fixing z. Since both factors are diagonal involutions,

```math
Q S Q=Q_{\rm low}S Q_{\rm low},                       \qquad(6)
```

which is Clifford, including its actual scalar. Furthermore
$`\mathcal G_Q=\mathcal C_2\cap Q\mathcal C_2Q`$ is a group.
Consequently R is in that group if and only if M is. For a monomial M,
QRQ has the same affine permutation and its phase changes by the sign
$`(-1)^{f(\pi(v))+f(v)}`$. A diagonal sign is Clifford exactly when
its Boolean exponent has degree at most two. This proves (1) in both
directions, without constraining address motion to be absent in M.

<a id="characteristic-address-algebra"></a>

## 5. The hard algebra is characteristic, and so is the address algebra

First, every R in $`\mathcal G_Q`$ normalizes the entire hard-Z algebra.
Write the affine permutation in (1) as $`\pi(v)=Bv+b`$. Adding a Boolean
quadratic to f does not change $`\mathcal L_f`$: every derivative of that
quadratic is affine. Affine substitution gives

```math
\mathcal L_{f\circ\pi}=B^{-1}\mathcal L_f.
```

Thus (1) forces $`B\mathcal L_f=\mathcal L_f`$. By (4), B preserves
the soft translation space. Every hard output of pi consequently depends
only on hard inputs. M normalizes the hard-Z algebra, and S fixes it
pointwise. In particular, the remaining affine condition can equivalently
be written

```math
\deg\bigl(h(\alpha(z))+h(z)\bigr)\le2,
\qquad h(z)=\sum_{|S|\ge2}c_Sx^S,                      \qquad(7)
```

where alpha is pi's induced invertible affine map on the hard coordinates.
The omitted low factor is quadratic before and after affine substitution.
This is algebra normalization, not pointwise fixation of every hard value.

**Address-algebra theorem.** If $`p\ge3`$, every R in $`\mathcal G_Q`$
normalizes the address-Z algebra itself. Thus its conjugates of address
Z generators are signed products of address Z generators, with no bank or
helper factors. The induced address permutation may be any affine map
that meets the remaining condition (1); it need not be the identity.

Here is an invariant proof that keeps Boolean degree reductions explicit.
Let V be the vector space of linear coordinate forms on the full physical
register. The associated graded algebra of the Boolean-degree filtration
is the commutative square-zero, or characteristic-two exterior, algebra
on V. Indeed $`v_i^2=v_i`$ lowers degree, so a square vanishes in the
leading symbol. An invertible affine substitution acts on degree k symbols
by the exterior power of its linear pullback; affine offsets affect only
lower degree.

Put $`d=c_{[p]}`$ and $`\mathcal A=\langle x_1,\ldots,x_p\rangle`$.
The unique degree-(p+1) symbol of f is

```math
\omega=d\wedge x_1\wedge\cdots\wedge x_p,
\qquad \mathcal U=\langle d,x_1,\ldots,x_p\rangle.      \qquad(8)
```

Because the defect in (1) is quadratic, pi preserves omega and therefore
its essential linear span $`\mathcal U`$. For an exterior form, the
essential span is the span of all contractions leaving one linear form;
for a collection, take the sum of those spans. This construction is
equivariant under invertible linear maps. In particular it recovers
$`\mathcal U`$ from (8).

Let $`\mathcal E=\mathrm{ann}(\mathcal U)`$ be the translation subspace
on which all address forms and d vanish. It is B-invariant. Taking the
finite difference **before** affine substitution gives the exact identity

```math
\Delta_a f=\sum_{S\subsetneq[p]}(a_c)_Sx^S,
\qquad a\in\mathcal E.                                \qquad(9)
```

The top term has already disappeared. These differences have degree at
most p-1, and their degree-(p-1) symbols span exactly

```math
\mathcal W=\Lambda^{p-1}\mathcal A.
```

Writing $`q=f\circ\pi+f`$, finite differences obey

```math
(\Delta_{Ba}f)\circ\pi=\Delta_a f+\Delta_a q.
```

The last term is affine. For $`p\ge3`$, the degree-(p-1) layer is of
degree at least two, so it is unaffected by that error. Since B permutes
$`\mathcal E`$, its linear pullback preserves $`\mathcal W`$. The
essential span of $`\Lambda^{p-1}\mathcal A`$ is $`\mathcal A`$:
the wedges omitting one address form already have essential spans whose
sum is all of $`\mathcal A`$. Therefore pi preserves the address linear
forms. Its affine offsets supply only the signed Z factors. S fixes all
address values, completing the proof for R.

At p=3 the preserved derivative layer is quadratic and $`\Delta_a q`$
is only affine, so the threshold includes p=3. At p=2 that layer is
linear and can change under the defect; the conclusion is false, as the
two-address example below shows. Approximate consequences use the separate
rounding lemma next, rather than an approximate degree-symbol argument.

<a id="uniform-boolean-sign-rounding"></a>

## 6. A uniform Clifford rounding gap for every Boolean sign

This lemma is not specific to lookup. Let $`Q_f|v\rangle=(-1)^{f(v)}|v\rangle`$
for any Boolean function f on all physical wires. For Cliffords R and C,
including their actual phases,

```math
\|Q_f R Q_f-C\|\lt\frac1{\sqrt2}
\quad\Longrightarrow\quad Q_f R Q_f\in\mathcal C_2.     \qquad(10)
```

Thus, if $`Q_fRQ_f`$ is not Clifford, its distance from every Clifford is
at least $`1/\sqrt2`$. The statement remains true after optimizing a
common phase on C. Its radius is independent of width and helper count.
All norms are complete operator norms, including arbitrary reference systems.
The conclusion is Clifford membership, not literal equality to the particular C.

**The exact Pauli gap.** For any Hermitian Pauli P, $`Q_f P Q_f`$ has
the same basis permutation as P, and its phases are P's phases times a
Boolean derivative sign. If this image is not Pauli, its distance from
every signed Hermitian Pauli is at least square root two. A different X
mask gives orthogonal column destinations. With the same X mask, multiply
by $`P^\dagger`$: the comparison is between a Boolean sign and a diagonal
Pauli with scalar in $`\{1,-1,i,-i\}`$. An imaginary scalar gives distance
square root two in every entry; a real scalar either gives exact equality
or disagrees by two on some column. This is the general form of the
[complete Pauli-image gap](INDICATOR_PAULI_NULLITY.md#3-a-uniform-pauli-image-gap-makes-the-count-robust).

Suppose the norm in (10) is epsilon below the stated radius. Put
$`P_i=RZ_iR^\dagger`$. Since $`Q_f`$ commutes every computational Z,
unitary telescoping gives

```math
\|Q_fP_iQ_f-CZ_iC^\dagger\|\le2\varepsilon\lt\sqrt2.
```

The Pauli gap forces every $`Q_fP_iQ_f`$ to be exactly Pauli. Define
$`\mathcal L_f`$ by (3), for this general f. It is a linear translation
space: the sum of two affine-derivative directions again has an affine
derivative. The full label space of Paulis with Pauli images under $`Q_f`$
is precisely

```math
\mathcal P_f=Z_{\rm all}\oplus X_{\mathcal L_f},
\qquad \mathcal P_f^\perp=Z_{\mathcal L_f^\perp}.
```

The labels $`R\langle Z\rangle R^\dagger`$ form a Lagrangian contained
in $`\mathcal P_f`$. Its symplectic centralizer is itself, so it contains
$`Z_{\mathcal L_f^\perp}`$. Consequently R pulls those latter Z labels
back to computational Z Paulis exactly.

Choose invertible linear coordinates with soft translations $`\mathcal L_f`$
and complementary hard coordinates z. Every soft coordinate derivative of
f is affine. In the Boolean ANF, this excludes every monomial of degree at
least three that contains a soft coordinate. Hence

```math
f=q+h(z),\qquad \deg q\le2.
```

The exact pullbacks above now give the same factorization $`R=SM`$ as
Section 4, with M monomial Clifford and S fixing every hard Z literally.
The factor $`Q_h`$ commutes S and $`Q_q`$ is Clifford, so
$`C_S=Q_fSQ_f=Q_qSQ_q`$ is Clifford. It follows that the diagonal sign

```math
D=M^\dagger Q_f M Q_f
```

is within epsilon of the Clifford $`M^\dagger C_S^\dagger C`$.
The retained [diagonal-sign threshold](CLIFFORD_INTERLEAVED_LOOKUP.md#2-two-elementary-clifford-thresholds)
says that a diagonal sign within strict distance
$`\delta=2\sin(\pi/8)=\sqrt{2-\sqrt2}`$ of a Clifford, even allowing
its common scalar, is itself Clifford. Since
$`1/\sqrt2\lt\delta`$, D is Clifford. Multiplying the actual factors
back shows that $`Q_fRQ_f`$ is Clifford, proving (10).

The radius is not claimed optimal. For p=2 lookup, take $`R=H_{c_{12}}`$.
The low Clifford factor commutes R, and the exact conjugate is

```math
\begin{aligned}
A=QRQ&=\frac{Z_{c_{12}}+X_{c_{12}}\mathrm{CZ}_{x_1,x_2}}{\sqrt2},\\
(A-Z_{c_{12}})^\dagger(A-Z_{c_{12}})&=(2-\sqrt2)I.
\end{aligned}                                                   \qquad(11)
```

A is not Clifford: it maps $`Z_{c_{12}}`$ to the non-Pauli Clifford
$`X_{c_{12}}\mathrm{CZ}_{x_1,x_2}`$. Nevertheless its distance from
$`Z_{c_{12}}`$ is exactly delta. This retains the distinction between
the proved lower radius and a concrete nearby non-Clifford example.

## 7. The two-address exception and genuine nonmonomial symmetries

At p=2, the high factor is the cubic sign
$`(-1)^{c_{12}x_1x_2}`$. Swapping $`x_1`$ and $`c_{12}`$ preserves
that factor exactly. The lower factor changes only by
$`(-1)^{c_1(x_1+c_{12})}`$, a Clifford sign. This swap therefore belongs
to $`\mathcal G_Q`$: the normal form does not require the affine
permutation to preserve the address algebra.

It can also be combined with a genuinely nonmonomial S, for example a
Hadamard on a singleton coefficient and a Hadamard on an arbitrary helper,
with hard-controlled CNOTs into that helper. Equation (6) still applies.
The complete word mixes address and bank coordinates through M and has
nonmonomial helper/bank fibers through S. Its lookup conjugate remains
Clifford on every physical input. Conversely, a lone Hadamard on an
address or a coefficient of degree at least two violates (3)--(4).

## 8. The necessary screen at exact and finite probe error

Suppose the [two-outer-conjugation search](CLIFFORD_INVOLUTION_SOURCE.md#two-outer-conjugations)
succeeds exactly, so

```math
Q=A K A^\dagger,\qquad A=C_3D_3C_2D_2,\qquad K\text{ Hermitian Clifford}.
```

Let P be a computational Z Pauli for which $`KPK^\dagger`$ is also
Z type, with its literal sign retained. Both

```math
R=A P A^\dagger,\qquad
Q R Q=A(KPK^\dagger)A^\dagger                         \qquad(12)
```

are Cliffords: $`D_2`$ commutes each input probe, $`C_2`$ maps it to
a Pauli, and the third-level $`D_3`$ maps that Pauli to a Clifford.
Thus every such R must obey (1)--(4), normalize the hard-Z algebra, and,
for p at least three, normalize the address-Z algebra. This supplies a
concrete necessary condition on proposed outer factors.

The same exact conclusions hold when the complete probe has finite error

```math
\|A K A^\dagger-Q\|\le\eta\lt\frac1{2\sqrt2}.
```

The two operators $`R=APA^\dagger`$ and
$`R'=A(KPK^\dagger)A^\dagger`$ are still actual Cliffords by the same
argument. Telescoping the approximate probe gives
$`\|QRQ-R'\|\le2\eta\lt1/\sqrt2`$. Lemma (10) therefore puts
$`R\in\mathcal G_Q`$ exactly. Its hard-algebra normalization and, for
p at least three, its address-algebra normalization follow from the exact
classification. This is a statement about complete probe error, not about
an individual matrix entry or an unreturned source block.

If K is monomial, the screen applies to the whole computational Z algebra.
For a genuinely nonmonomial K, the
intersection of that algebra with its K image may be smaller or even
trivial. The argument does not replace non-Clifford intermediate images
by Cliffords or exclude arbitrary three-block sources. All outer phase counts,
Clifford gates, arbitrary work, actual inverses and complete return must
still be paid. The [exact fixtures](../../tests/test_lookup_clifford_symmetries.py)
check the affine-derivative space, a complete genuinely mixed two-address
Clifford example, the literal low/high conjugation, the characteristic
degree symbols and their threshold, excluded support directions, and the
literal non-Clifford distance witness (11).
No general frame frontier changes here.
