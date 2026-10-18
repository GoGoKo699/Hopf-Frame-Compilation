# Five-layer reduction and conditional degree bounds

[Current depth frontier](../../docs/OPEN_PROBLEM.md) ·
[Three- and four-layer proof](READ_ONLY_THREE_LAYER_RIGIDITY.md) ·
[Read-only Clifford normal form](READ_ONLY_TWO_LAYER_RIGIDITY.md)

A five-layer read-only query admits a useful reformulation: its target
comparison can remain Clifford-valued, while a two-layer Pauli image
acts by conjugation on that comparison family. The latter conjugation
need not preserve the full Pauli group. Shared native prefixes nevertheless
give bounded-degree intrinsic comparison features. The remaining general
obstacle is extracting the Boolean target relative to its outer native
baseline; no closure of the Clifford hierarchy is assumed.

This note proves a **conditional** five-layer degree bound. If the
specified middle two-layer Pauli representation is Clifford-valued,
then the Boolean XOR output has degree at most seven, independently of
dirty width. A visible outer target shear improves this to degree at
most four. The condition is additional; full workspace return does not
supply it. Neither bound is claimed sharp or sufficient for synthesis.
A separate subcase needs no Clifford-valued middle: if the comparison K
is Pauli-valued, degree is at most two at error below $`\sqrt2/12`$.
Native peeling freezes that subgroup action exactly projectively;
independent inner axes then exclude cubic degree. No exact normalizer
is assumed for a general approximate query. A further native subcase, where
K permutes the outer active axes, has a conservative degree bound of 244
(or 41 for an invisible target shear). It allows both non-Pauli K and
non-Clifford B; its extra hypothesis remains explicit.
The weaker condition that K preserves the outer axis **span** now also has
a finite degree bound: **112**, at the smaller error radius
$`\eta\lt\sqrt2/24`$. Section 5.4 constructs intrinsic group-algebra
features from the native comparison subgroup. The earlier stronger bounds
and the larger-radius axis-normalizer theorem remain available under their
own hypotheses. Section 5.5 proves faithful degree-four intrinsic comparison
features without a span hypothesis below $`\sqrt2/40`$. Combined with span
preservation, these improve the output bound to **21** at that smaller radius.
Target extraction without span preservation remains open.

## 1. Native contract and an exact one-band reduction

Let x be an a-bit address, and let W contain the arbitrary output t and
all arbitrary dirty helpers. There are no initialized work qubits. Take

```math
V=C_5D_5C_4D_4C_3D_3C_2D_2C_1D_1C_0,
\qquad F_x=X_t^{f(x)}\otimes I_{\rm dirty}.
```

Every completed Clifford interlayer preserves the computational address
pointwise. Every D is one physical layer of T, T-dagger and identity on
distinct work wires. Empty layers are allowed. Assume

```math
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}\le\eta.
```

The norm includes all output and dirty inputs, arbitrary external
references, leakage and complete work return. T-count, elementary
Clifford count and finite dirty width are unrestricted. The proof adds
no physical gate or auxiliary wire. Its fixed changes of Pauli
coordinates are analytical operations, not free compiler subroutines.

Use the retained address-preserving normal form at every boundary:

```math
C_j(x)=\lambda_j(x)P_j(x)C_j(0),\qquad
P_j(0)=I,\quad\lambda_j(0)=1.
```

Each P has binary-linear labels; its projective multiplication factors
are bilinear signs. Every address S and CZ phase remains in lambda.
Write

```math
\lambda(x)=\prod_{j=0}^{5}\lambda_j(x),\qquad
W_5=C_5(0)D_5,\qquad
R_x=C_0(0)^\dagger P_0(x)C_0(0).
```

Define three exact projective Pauli representations by

```math
A_x=C_4(0)D_4P_3(x)D_4^\dagger C_4(0)^\dagger,
```

```math
B_x=C_4(0)D_4C_3(0)D_3P_2(x)D_3^\dagger C_3(0)^\dagger
D_4^\dagger C_4(0)^\dagger,
```

```math
J_x=C_4(0)D_4C_3(0)D_3C_2(0)D_2P_1(x)D_2^\dagger C_2(0)^\dagger
D_3^\dagger C_3(0)^\dagger D_4^\dagger C_4(0)^\dagger.
```

Only A is automatically Clifford-valued. Direct multiplication gives,
for the actual work fibre $`U_x=V(x)`$,

```math
U_x=\lambda(x)P_5(x)W_5P_4(x)A_xB_xJ_xW_5^\dagger U_0R_x.
                                                               \qquad(1)
```

Define the phase-retaining target comparisons

```math
M_x=P_4(x)^\dagger W_5^\dagger
 [\lambda(x)^{-1}P_5(x)^\dagger F_xR_x^\dagger F_0^\dagger]W_5,
\qquad K_x=A_x^\dagger M_x.                                  \qquad(2)
```

Each M and K is a specified scalar times a Clifford: the bracket in
(2) is a phase times a Pauli, and only one physical T layer conjugates
it. Solving (1) with the **actual inverse** $`U_0^\dagger`$,
then replacing the two target approximations, gives

```math
\|B_xJ_x-K_x\|\le2\eta.                                    \qquad(3)
```

If $`B_xB_y=b(x,y)B_{x+y}`$ and
$`J_xJ_y=j(x,y)J_{x+y}`$, put
$`\zeta(x,y)=b(y,x)j(x,y)`$. Their exact representation laws,
followed by three uses of (3), give

```math
\|K_{x+y}-\zeta(x,y)B_yK_xB_y^\dagger K_y\|\le6\eta.
                                                               \qquad(4)
```

For an exact query, (4) is an exact equality including its specified
scalar. Thus, even when B is non-Clifford, it conjugates the particular
Clifford comparison family back into Cliffords. Let G be the projective
group generated by all K. Equation (4) at zero error implies

```math
B_yGB_y^\dagger=G.                                          \qquad(5)
```

Inclusion follows on the generators from (4); equality follows because
conjugation by $`B_y`$ is involutive. Equation (5) is a partial
normalizer statement. It does not assert that B normalizes every Pauli
or every Clifford, and no approximate version of (5) is asserted here.

## 2. Additional hypothesis and a robust conditional theorem

**Conditional theorem.** In the native setting above, suppose the
specified family $`B_y`$ is Clifford-valued for every address y.
If

```math
\eta\lt\frac{\sqrt2}{12},
```

then $`\deg_{\mathbb F_2}f\le7`$. If the target is visible in
the outer shear, in the precise sense $`E\ne0`$ of section 3,
then $`\deg_{\mathbb F_2}f\le4`$.

The additional condition can be stated on the address generators:
B is a projective representation, so Clifford-valued generator images
make every B Clifford. This is a condition on the displayed native
representation, not a claimed intrinsic characterization of five-layer
queries. No optimal general depth statement follows from it.

Under this condition both sides of (4) are projective Cliffords. The
[retained Clifford separation argument](READ_ONLY_THREE_LAYER_RIGIDITY.md#81-peeling-the-outer-band-and-quantizing-a-clifford-cocycle)
first forces their symplectic actions to agree, since
$`6\eta\lt1/\sqrt2`$. Their quotient is then a scalar times a
Pauli; its label must vanish since $`6\eta\lt\sqrt2`$.
Consequently the following cocycle is exact **projectively**:

```math
K_{x+y}\doteq B_yK_xB_y^\dagger K_y.                       \qquad(6)
```

Equality of the residual scalars is neither asserted nor needed. The
literal target phases and full operator error remain those of (1)--(4).

## 3. Visible outer shear: a translation invariant space of functions

Put $`g(x)=f(x)+f(0)`$. The outer one-T Pauli image in (2) has
symplectic matrix $`I+L(x)+E g(x)`$, where L is linear and E
is fixed. Both belong to the outer physical shear space. The
one-T family A has symplectic action $`I+H(x)`$, with H linear
and $`H(x)^2=0`$. Therefore

```math
S_K(x)=(I+H(x))(I+L(x)+E g(x)).                            \qquad(7)
```

Assume $`E\ne0`$. Let $`\mathcal P_r`$ denote the Boolean
polynomials in x of degree at most r, and let $`\mathcal W`$
be the binary span of all matrix-entry functions of $`S_K(x)`$.
Equation (7) gives

```math
\mathcal W\subseteq\mathcal P_2+g\mathcal P_1.             \qquad(8)
```

The symplectic action of (6) is

```math
S_K(x+y)=S_B(y)S_K(x)S_B(y)^{-1}S_K(y).                   \qquad(9)
```

For fixed y the right-hand side is a fixed linear transformation of
the entries of $`S_K(x)`$. Thus $`\mathcal W`$ is invariant
under translation of x.

A fixed linear functional extracts g from (7) up to a quadratic
polynomial. To see this, choose a nonzero column v of E. Every
$`H(x)v`$ lies in the fixed isotropic native-axis space
$`\mathcal L_A`$, and H vanishes on that space. If
$`v\in\mathcal L_A`$, all these vectors vanish. Otherwise v
is outside their span. In either case there is a binary linear
functional l with

```math
l(v)=1,\qquad l(H(x)v)=0\quad\text{for every }x.
```

Apply l to the corresponding column of (7). Its coefficient of g is
one and its other term is quadratic, proving

```math
g\in\mathcal W+\mathcal P_2.
```

Translation invariance and (8) now imply

```math
g(x+y)\in\mathcal P_2+g\mathcal P_1\quad\text{for every }y.
                                                               \qquad(10)
```

Let $`d=\deg g`$. Restrict x to the d variables of a nonzero
degree-d monomial, setting the other variables to zero. The degree
remains d, and (10) restricts to this address subgroup. The translates
of a d-variable full-degree Boolean polynomial span a space of
dimension $`2^d`$: the derivatives indexed by subsets S have
leading monomials $`\prod_{i\notin S}x_i`$, with only lower-degree
terms beneath them, so they form a triangular basis of all Boolean
functions. Each derivative is a sum of translates.

But the containing space in (10) has dimension at most

```math
\dim(\mathcal P_2+g\mathcal P_1)
\le1+d+\binom d2+(d+1).
```

For $`d\ge5`$ this is strictly smaller than $`2^d`$.
Constant f already satisfies the conclusion. Otherwise
$`\deg f=\deg g\le4`$. This argument uses native
independent axes for A; it does not assign such axes to B.

## 4. The sign of a fixed Clifford on Hermitian Paulis

For binary labels $`v=(p,q)`$, use the Hermitian section

```math
\tau(v)=i^{p\cdot q}X^pZ^q,
```

where the dot product in the exponent is integer-valued. For a fixed
Clifford V define

```math
V\tau(v)V^\dagger=(-1)^{\kappa_V(v)}\tau(S_Vv).
```

**Sign lemma.** The Boolean function $`\kappa_V`$ has degree
at most three in the binary entries of v.

In the ordered XZ section, the conjugation phase modulo four is an
integer linear polynomial plus twice a quadratic polynomial. Conversion
to the Hermitian section adds $`p\cdot q-p'\cdot q'`$,
where $`(p',q')=S_Vv`$. Modulo four, a binary parity has degree
at most two, with even quadratic coefficients. The product of two
such parities has degree at most three modulo four: products of two
quadratic terms have coefficients divisible by four. Thus the full
phase has a multilinear polynomial of degree at most three modulo four.
It is pointwise even by Hermiticity. Uniqueness of the multilinear
polynomial representation modulo four makes every coefficient even.
Dividing them by two gives a binary polynomial of degree at most three,
proving the lemma. Width does not enter the degree bound.

## 5. Invisible outer shear: degree at most seven

Suppose $`E=0`$. Independence of the outer physical axes implies
that the target commutes with every active outer T axis after the fixed
outer Clifford. Hence

```math
Q=W_5^\dagger X_tW_5
```

is a fixed nonidentity Pauli, with binary label v. Moving this fixed
factor to the right in (2), while retaining any induced phase in the
actual comparison, gives projectively

```math
K_x\doteq K_x^0Q^{g(x)},\qquad
K_x^0=A_x^\dagger P_xN_x.                                 \qquad(11)
```

Here P has linear Pauli labels and N is a fixed Clifford conjugate of
a one-physical-T-conjugated linear Pauli representation. The symplectic
matrix $`S_K(x)`$ is independent of g and has degree at most two
in x, by (7).

For a fixed Hermitian probe $`\tau(q)`$, let $`p_q(x)`$ be
its conjugation sign under $`K_x^0`$. A conservative bound sufficient
here is

```math
\deg p_q\le6.                                             \qquad(12)
```

To check this directly, N sends a fixed probe to a Pauli with affine
label in x. In N's physical-axis coordinates its sign is linear in x;
its fixed output Clifford adds the cubic sign from section 4 evaluated
on an affine label. P adds a pairing of two affine labels, of degree at
most two, without changing the current label. The fixed input Clifford
of A adds another cubic sign on an affine label and changes that label
linearly. On one active physical wire of A, let $`(u(x),w(x))`$
be its Pauli family label and $`(a(x),b(x))`$ the incoming probe
label. All four functions are affine. The one-T Pauli image contributes

```math
u b(1+a)+w a+\epsilon u a,
\qquad\epsilon\in\{0,1\},                               \qquad(13)
```

where epsilon distinguishes the two T orientations. This sign has
degree at most three; the resulting probe label has degree at most two.
The final fixed Clifford of A consequently adds a sign of degree at
most six by section 4. Inverse A has the same conjugation action, since
A is a Pauli image up to scalar. Adding these signs proves (12).

Let $`r_q(x)`$ be the sign of the full K on a fixed probe.
Equation (11) gives

```math
r_q(x)=p_q(x)+\langle v,q\rangle g(x).                    \qquad(14)
```

Apply the exact projective cocycle (6) by conjugation to that probe.
For fixed y, put

```math
q'_y=S_B(y)^{-1}S_K(y)q.
```

The signs contributed by $`K_y`$ and $`B_y^\dagger`$ are
constant in x. The remaining signs give

```math
r_q(x+y)=r_{q'_y}(x)
 +\kappa_{B_y}\bigl(S_K(x)q'_y\bigr)+c(y,q).              \qquad(15)
```

The middle term has degree at most six: the fixed Clifford sign is
cubic, and its binary label argument has degree at most two. Choose
q with $`\langle v,q\rangle=1`$, which exists because v is
nonzero. Substitution of (12) and (14) into (15) yields

```math
g(x+y)=\epsilon_y g(x)+h_y(x),\qquad
\epsilon_y\in\mathbb F_2,\qquad\deg h_y\le6.              \qquad(16)
```

If $`\deg g\gt6`$, translation preserves this degree, so
$`\epsilon_y=1`$ for every y. All Boolean derivatives of g
therefore have degree at most six. Every monomial of degree at least
eight would have a coordinate derivative of degree at least seven,
which is impossible. Thus $`\deg g\le7`$. If its degree was
already at most six, the conclusion holds directly. Together with
section 3, this proves the conditional theorem.

## 5.1. A Pauli-valued comparison: robust quadratic bound

**Pauli-comparison theorem.** In the native setting of section 1, suppose
that every specified comparison $`K_x`$ in (2) is a scalar times a Pauli.
If $`\eta\lt\sqrt2/12`$, then $`\deg f\le2`$, independently of dirty
width. The family B may be non-Clifford. This is a different additional
hypothesis from section 2; neither hypothesis is inferred from full return.
Attainability within this comparison subcase is not asserted.

Let $`g(x)=f(x)+f(0)`$. If g is affine there is nothing to prove. Otherwise
$`S_K(x)=I`$ and (7), multiplied on the left by $`I+H(x)`$, give

```math
L(x)+E g(x)=H(x).
```

The two matrix-valued functions L and H are linear. A nonzero entry of E
would make g linear, so $`E=0`$ and $`H=L`$. The factorization (11) becomes

```math
K_x\doteq A_x^\dagger P_xN_xQ^{g(x)},
\qquad S_A(x)=S_N(x)=I+H(x).
```

A and N are exact projective linear-label Pauli representations conjugated
by one physical T layer and fixed Cliffords. Consequently
$`C_x=A_x^\dagger N_x`$ is Pauli-valued. Its binary label c obeys

```math
c(x+y)=(I+H(y))c(x)+c(y),\qquad c(0)=0.
```

This follows from the two exact representation laws, as in (8) of the
[three-layer proof](READ_ONLY_THREE_LAYER_RIGIDITY.md#2-a-robust-pauli-label-cocycle).
The independent physical axes of A force $`H(y)c(x)=0`$ for all x,y:
in native coordinates, setting y=x and reversing x,y give, on each active
axis, $`u_j(x)p_j(x)=0`$ and
$`u_j(y)p_j(x)=u_j(x)p_j(y)`$. Choosing y with $`u_j(y)=1`$ makes
$`p_j`$ vanish; a zero $`u_j`$ contributes nothing. Hence c is linear.

Write p for the linear Pauli label of P and v for the fixed nonzero label
of Q. Conjugating P by A gives the full comparison label

```math
m(x)=r(x)+v g(x),\qquad
r(x)=(I+H(x))p(x)+c(x),\qquad \deg r\le2.                 \qquad(18)
```

The quadratic coefficient vectors of r lie in the fixed span of A's
physical T axes, since its only quadratic term is $`H(x)p(x)`$.

### Native peeling freezes the Pauli cocycle

Make the fixed analytical change of Pauli coordinates by
$`C_4(0)^\dagger`$, and use the same letters for the conjugated operators
and their labels. Put

```math
D=D_4,\qquad
A'_y=C_3(0)D_3P_2(y)D_3^\dagger C_3(0)^\dagger.
```

Then $`B_y=D A'_yD^\dagger`$, and A' is an exact native one-T
Clifford Pauli family. The axes whose span contains the quadratic
coefficients in (18) are precisely the active axes of D. Denote their
label span by $`\mathcal L_D`$. This coordinate change adds no physical
gate or resource.

For every z,y, right-multiply (4) by $`K_y^\dagger`$. This gives

```math
\|K_{z+y}K_y^\dagger-\zeta(z,y)B_yK_zB_y^\dagger\|
\le6\eta.                                               \qquad(19)
```

Conjugate this comparison by $`D^\dagger`$. Its two sides are

```math
D^\dagger K_{z+y}K_y^\dagger D,
\qquad
\zeta(z,y)A'_y(D^\dagger K_zD)A_y'^\dagger.
```

Both are scalar times Cliffords, because the input to each D-conjugation
is a Pauli. Their distance is at most $`6\eta\lt1/\sqrt2`$. The retained
Clifford gap first forces equal symplectic actions. Their quotient is
then a scalar times a Pauli; a nonidentity quotient would have
projective distance at least $`\sqrt2`$, so its Pauli label also vanishes.
Undoing the conjugation proves (19) exactly **projectively**. The
remaining literal scalar need not be equal.

Thus each B_y exactly normalizes the projective Pauli group generated
by all K_x, even though it need not normalize the full physical Pauli
group. Let W be that subgroup's binary label span, and let T_y be its
induced invertible linear action on W. The exact label cocycle is

```math
m(x+y)=T_y m(x)+m(y).                                    \qquad(20)
```

No intermediate Clifford property was assumed for a general K. Here
it follows specifically from the Pauli-valued comparison hypothesis.

### A cubic is the only remaining possibility

The quadratic baseline in (18) sums to zero on any three-dimensional
address cube. Multiplying the eight exact projective versions of (19)
in a fixed order cancels adjacent $`B_y^\dagger B_y`$ on one side and
cancels the eight Pauli factors $`K_y^\dagger`$ projectively on the other.
For $`d(x)=\Delta_{h_1}\Delta_{h_2}\Delta_{h_3}g(x)`$, this gives

```math
Q^{d(x+y)}\doteq B_yQ^{d(x)}B_y^\dagger.
```

Identity cannot be conjugate to a nonidentity Pauli up to scalar.
Hence every third Boolean derivative of g is translation invariant,
so $`\deg g\le3`$. Equivalently, this follows by comparing the highest
coefficients in (20) and then taking one Boolean derivative. Only
binary labels and necessary projective equalities are used; the
original literal error target has not changed.

### Inner independent axes exclude cubic degree

Suppose for contradiction that g has degree three. A cubic coefficient
of m is v, so $`v\in W`$. Comparing that coefficient in (20) gives
$`T_yv=v`$ for every y. Every quadratic coefficient vector w of r also
lies in W: subtract $`v g_{ij}`$ from the corresponding coefficient
of m. As observed after (18), it also lies in $`\mathcal L_D`$.

For any $`w\in W\cap\mathcal L_D`$, its Pauli commutes with D.
The image under $`B_y=D A'_yD^\dagger`$ is again a Pauli by the
exact subgroup normalization. Therefore

```math
D(A'_y\tau(w)A_y'^\dagger)D^\dagger
```

is Pauli. Independence of D's active axes implies that a Pauli has a
Pauli image under D only if it commutes with every active axis: otherwise
its expansion has at least two nonzero coefficients, as proved in the
[one-band matching theorem](ONE_BAND_AXIS_MATCHING.md#2-recovering-the-axis-list-from-pauli-coefficient-norms).
Consequently D leaves this Pauli label unchanged, and

```math
T_yw=S_{A'}(y)w\quad\text{for every }w\in W\cap\mathcal L_D.
                                                               \qquad(21)
```

Taking quadratic coefficients in (20), using (18), $`T_yv=v`$
and (21), yields

```math
v(\Delta_y g)^{[2]}
=(I+S_{A'}(y))r^{[2]}
\quad\text{coefficientwise}.                            \qquad(22)
```

The independent native axes of A' have labels a_j and linear activation
functions u_j(y), so

```math
(I+S_{A'}(y))w=\sum_j u_j(y)\langle a_j,w\rangle a_j.
```

A nonzero cubic has a nonzero derivative quadratic part. Equation (22)
therefore places v in the span of these independent axes. Choose an
axis coordinate k with coefficient one in v. For every quadratic
coefficient w_ij of r, its corresponding coefficient in the derivative
must satisfy

```math
(\Delta_y g)_{ij}=u_k(y)\langle a_k,w_{ij}\rangle.          \qquad(23)
```

Thus the linear map sending y to the entire quadratic part of
$`\Delta_y g`$ has rank at most one. But a nonzero cubic has contraction
rank at least three: choose a nonzero coefficient of $`x_i x_j x_k`$;
restrict the directions to i,j,k and the output coordinates to the
pairs jk,ik,ij. The resulting matrix is the three-by-three identity.
This contradiction excludes cubic degree, proving $`\deg f\le2`$.

The proof uses full operators and the same fixed approximation radius
as the retained four-layer bound. It adds no gate or workspace, assumes
no intermediate work return, and retains arbitrary output, dirty and
reference inputs, leakage and the actual inverse in the original norm.
It rounds a specified conjugation action on the generated Pauli subgroup;
it does not round B itself to a Clifford. Non-Pauli K still leaves the
general five-layer problem open.

## 5.2. An outer-axis-normalizer comparison: a finite degree bound

**Axis-normalizer theorem.** In the native setting of section 1, make the
fixed analytical change of coordinates by $`C_4(0)^\dagger`$ and put
$`D=D_4`$. Suppose every specified comparison K maps D's individual active
binary Pauli axes bijectively to themselves, allowing signs. At
$`\eta\lt\sqrt2/12`$, the Boolean output has degree at most **244**,
independently of dirty width. If the original outer target shear is
$`E=0`$, the bound improves to **41**. If also $`D^\dagger QD`$ is Pauli,
it improves to **31**. These conservative constants are not sharpness
claims. The stronger degree-two theorem of section 5.1 still applies
when K itself is Pauli.

This is an additional native hypothesis, not a consequence of full return.
It permits both non-Pauli K and non-Clifford B. The returned affine fixture
in section 6 has $`D=T`$ and comparison $`K_1=Z(T^\dagger XT)`$.
This comparison is non-Pauli but maps Z to minus Z; its displayed B is
non-Clifford. Thus the new hypothesis is not subsumed by either previously
proved subcase. The fixture remains a deliberately expensive realization
of a single CNOT, with no resource improvement.

### Peeling and a finite Clifford feature space

The [one-band matching theorem](ONE_BAND_AXIS_MATCHING.md) makes
$`C_x=D^\dagger K_xD`$ Clifford. Axis permuters form a group, so
$`D^\dagger K_{x+y}K_y^\dagger D`$ is Clifford too. Right-multiply (4)
by $`K_y^\dagger`$ before peeling D. Both resulting sides are Clifford;
the same separation argument as section 5.1 freezes the equality
projectively:

```math
C_{x+y}\doteq A'_y C_x A_y'^\dagger C_y,
\qquad B_y=D A'_yD^\dagger.                              \qquad(24)
```

A' is the native inner one-T Clifford Pauli family. Thus translating x
acts by fixed left and right Clifford multiplication on C. No claim is
made that B itself is Clifford.

For a Clifford C let S be its binary symplectic matrix and let h list its
Hermitian generator signs. Define its feature vector $`\mathcal V(C)`$
to contain 1, all monomials of degree at most four in entries of S, and
the entries of h. Fixed left or right Clifford multiplication acts linearly
on the span of these features. The degree four here needs a different
argument from the fixed-Clifford sign lemma in section 4.

For a fixed probe q and a varying Clifford C, its sign is linear in h
plus a polynomial of degree at most four in S. To prove this, order X
generators before Z generators and write $`S e_j=(p_j,r_j)`$.
Let $`p'=\sum_jq_jp_j`$ and $`r'=\sum_jq_jr_j`$ be binary sums.
Beyond the generator signs, the conjugation phase is half the following
even quantity modulo four; the dot products displayed here are integer:

```math
\begin{aligned}
F={}&q_X\cdot q_Z+\sum_jq_j(p_j\cdot r_j)\\
 &+2\sum_{j\lt k}q_jq_k(r_j\cdot p_k)-p'\cdot r'.
\end{aligned}                                             \qquad(25)
```

A binary parity modulo four is its integer linear sum minus twice the
sum of pairwise products. Expanding (25), odd-coefficient monomials have
degree at most two in S and even-coefficient monomials degree at most
three. The high bit of a sum modulo four is the sum of the summands'
high bits plus pairwise products of their low bits. It therefore has
degree at most four. This argument works even though F is guaranteed
even only on the symplectic domain; it does not divide odd polynomial
coefficients by two.

For fixed Cliffords L and R, the symplectic matrix of LCR is linear in S.
Its generator signs use the varying-Clifford sign just proved, together
with the fixed-Clifford cubic sign evaluated on a label linear in S.
They are consequently linear in h plus degree-four polynomials in S.
This proves the feature-space assertion. Tensor powers of the features
also transform linearly.

### A polynomial formula on the axis-normalizer domain

Use D's physical Pauli coordinates, so its active axes are unsigned Z
operators and its T orientations retain every native axis sign. Choose
dual labels to these axes. The permutation entry
$`\pi_{ji}`$ of C on those axes is a pairing with $`S_Ca_i`$, hence
linear in S. Its sign $`\epsilon_i=\kappa_C(a_i)`$ belongs to the
feature span just defined. Let $`\sigma_j`$ record the orientation of
D's T or T-dagger on axis j. Up to scalar, the correcting Clifford in
$`DCD^\dagger C^\dagger`$ is a product of S or S-dagger rotations
about these independent axes, with activations

```math
\delta_j=\sigma_j+\sum_i\pi_{ji}(\sigma_i+\epsilon_i).
                                                               \qquad(26)
```

Its orientation on each axis is fixed by D. All terms in (26) are binary.
Thus delta is degree at most two in $`\mathcal V(C)`$, and

```math
S_{DCD^\dagger}
 =\left(I+\sum_j\delta_j N_j\right)S_C                 \qquad(27)
```

is degree at most three in those features. Here N_j is the fixed
transvection increment for axis j. The same statements hold with D and
its actual inverse exchanged. They are polynomial formulas on the
axis-normalizer domain; outside it they are only polynomial expressions,
not asserted physical tableaus.

For later degree accounting, a conservative generator-sign bound is nine
in the features: the incoming Pauli label has degree one, the independent
S rotations give a label of degree at most three, and a fixed output
coordinate Clifford contributes its cubic sign. The intermediate local
sign has degree at most four. This bound allows arbitrary fixed Pauli
coordinates and does not charge an analytical coordinate change as a gate.

### The visible-target bound

Temporarily replace $`g(x)=f(x)+f(0)`$ by a fresh Boolean variable t in
the native comparison. Projectively it has the form

```math
K(x,t)=A_x^\dagger P_x N_{(x,t)}.
```

P has linear Pauli labels; N is a fixed-Clifford conjugate of a one-T
image of a linear-label Pauli family. This defines a Clifford for every
pair x,t. Its symplectic entries have total degree at most two, and its
sign on every fixed probe has degree at most six. The latter is exactly
the direct calculation in section 5, allowing t among N's affine input
variables. It does not assume axis normalization at the other value of t.

Apply the polynomial expressions (26)--(27) to this formal family.
Permutation coordinates have degree at most two, and axis signs have
degree at most six, so delta has degree at most eight. The resulting
symplectic entries have degree at most ten. Generator signs have degree
at most thirty: a fixed coordinate Clifford's cubic sign is evaluated
on a label of degree at most ten. Therefore the features of
$`C_x=D^\dagger K_xD`$ have representatives of total degree at most
forty in x,t, evaluated at $`t=g(x)`$. Only at this actual value are they
claimed to be physical tableaus.

Suppose $`E\ne0`$. Write $`K_x\doteq K_x^0Q^{g(x)}`$ with the fixed
Clifford $`Q=W_5^\dagger X_tW_5`$. The baseline symplectic matrix S_0
and its inverse have degree at most two in x; symplectic inversion is a
fixed transpose-linear operation. A nonzero entry of
$`S_0^{-1}S_K-I`$ extracts g. C also permutes these same active axes,
because D commutes with them. Recover S_K from C by (27), of degree at
most three in its features. Consequently g belongs to the space
$`\mathcal W`$ spanned by degree-at-most-two monomials in x times
feature monomials of degree at most three. Equation (24) makes this
space translation invariant. Its functions have formal total degree at
most 122 in x,t, and hence

```math
g\in\mathcal W\subseteq\mathcal P_{122}+g\mathcal P_{121}.
                                                               \qquad(28)
```

For completeness, if a translation-invariant space contains g and is
contained in $`\mathcal P_d+g\mathcal P_{d-1}`$, then
$`\deg g\le2d`$. Restrict to the n variables of a highest-degree
monomial, where $`n=\deg g`$. Its translates span all $`2^n`$ Boolean
functions, as in section 3. But the containing space has dimension at most

```math
\sum_{j=0}^{d}\binom nj+\sum_{j=0}^{d-1}\binom nj,
```

strictly smaller than $`2^n`$ when $`n\ge2d+1`$. Applying this to (28)
proves the bound 244.

### The invisible-target improvement

If $`E=0`$, Q is a fixed nonidentity Pauli. Pauli multiplication changes
axis signs without changing their permutation, so the baseline $`K^0`$ also
normalizes D's axes. Thus

```math
C_x\doteq C_x^0(Q')^{g(x)},\qquad
C_x^0=D^\dagger K_x^0D,\qquad Q'=D^\dagger QD
```

are Cliffords. The preceding native calculation, now with no formal t,
gives degree at most ten for the symplectic entries of $`C^0`$, and at most
thirty for its sign on every fixed probe.

If Q' is non-Pauli, a nonzero entry of $`S_{C^0}^{-1}S_C-I`$ extracts g.
Take the span of degree-at-most-ten x monomials times the symplectic
entries of C and the constant. It is translation invariant by (24) and
lies in $`\mathcal P_{20}+g\mathcal P_{20}`$. Its dimension is too small
to contain all translates of a full-degree polynomial on 42 or more
variables, giving degree at most 41.

If Q' is Pauli, the symplectic matrix of C is independent of g and has
degree at most ten. Its sign on a fixed probe q is
$`p_q(x)+\langle v',q\rangle g(x)`$, where v' is Q's transformed
nonzero label and $`\deg p_q\le30`$. Apply (24) to a probe detecting v'.
The fixed conjugator's cubic sign on a degree-ten label has degree at most
thirty. As in (15)--(16), this yields
$`g(x+y)=\epsilon_y g(x)+h_y(x)`$ with $`\deg h_y\le30`$.
If g has degree greater than thirty, every $`\epsilon_y=1`$ and every
first derivative has degree at most thirty. Hence $`\deg g\le31`$.
This completes the stated bounds.

The proof adds no native operation or workspace. The original full operator
error includes every logical column, literal target phase, arbitrary dirty
and reference input, leakage, the actual inverse and complete work return.
T-count, elementary Clifford count and finite dirty width remain unrestricted
in this necessary invariant. The general comparison may fail the additional
axis-normalizer hypothesis; the general five-layer and complete-frame depth
questions remain open. The [bounded checks](../../tests/test_outer_axis_normalizer.py)
test the varying-tableau sign polynomial, signed band formula, exact returned
nonredundancy witness and dimension thresholds.

## 5.3. Preserving the axis span freezes a third-level cocycle

There is a further reduction when K preserves the binary span
$`\mathcal L_D`$ of D's active axes, while possibly mixing individual
axes. Under this weaker additional hypothesis, at
$`\eta\lt\sqrt2/24`$ the same peeled cocycle (24) becomes exact
projectively. Its comparisons C need only lie in the third level of the
Clifford hierarchy. This reduction supplies the exact cocycle used by the
degree theorem in section 5.4 below.

Let $`G_{\mathcal L_D}`$ be the Clifford group that preserves this
axis span. For any member K,

```math
D^\dagger KD=\left[D^\dagger(KDK^\dagger)\right]K.
                                                               \qquad(29)
```

All Pauli axes of the two bands in brackets belong to the same isotropic
span, so they commute. Their product is a collection of commuting T or
T-dagger Pauli rotations, with signs retained. Conjugating any Pauli by
this product gives that Pauli times commuting S or S-dagger Pauli
rotations, hence a Clifford. The bracket therefore lies in the third
hierarchy level, denoted $`\mathcal C_3`$. Multiplication by K preserves
this level. Consequently

```math
D^\dagger G_{\mathcal L_D}D\subseteq\mathcal C_3.          \qquad(30)
```

The left side is an actual group, because it is a conjugate of a group.
This asserts closure only for this specified subgroup; the full third
hierarchy level is not being treated as a group.

Two third-level unitaries at projective operator distance less than
$`\sqrt2/4`$ are projectively equal. Indeed, conjugating any Hermitian
Pauli by the two operators gives Cliffords at distance less than
$`1/\sqrt2`$. The retained Clifford separation makes those images equal
projectively. Their Hermiticity allows only a residual sign, and the same
strict distance excludes the minus sign. The images are therefore equal
literally on every Pauli, so the original operators differ only by scalar.
This gap is independent of work width.

Now right-multiply (4) by $`K_y^\dagger`$ and peel D. The left comparison
$`C_{x+y}C_y^\dagger`$ lies in the subgroup (30). The right comparison
$`A'_y C_x A_y'^\dagger`$ lies in $`\mathcal C_3`$ by fixed Clifford
conjugation. Their distance is at most $`6\eta\lt\sqrt2/4`$, so the
gap proves (24) projectively. It neither rounds B to a Clifford nor makes
the intermediate C Clifford. Already a CNOT on two active T wires
preserves their Z span but mixes individual axes: its D-dagger conjugate
is third-level and non-Clifford.

Section 5.2's Clifford tableau features do not close when C_y is
third-level. Section 5.4 resolves this branch by using a truncated group
algebra of the actual comparison subgroup instead. Its bounded native
degree depends on the two specified one-T baseline families. Abstract
partial normalization alone remains insufficient, by the
[counterexample](ONE_BAND_AXIS_MATCHING.md#5-an-abstract-normalizer-can-have-arbitrarily-high-tableau-degree).

The argument adds no native gate or workspace and retains the original
full operator contract. The commuting rotations in (29) are a structural
factorization, not a free compiler call. The
[axis-span fixture](../../tests/test_outer_axis_normalizer.py) checks every
Pauli image of the mixed-axis example, including its non-Clifford nature.

### 5.4 Span-preserving comparisons: a finite degree bound

**Span-preserving theorem.** In the native setting of section 1, suppose
that every specified K preserves the binary span $`\mathcal L_D`$ of
the active D axes, after the same fixed analytical coordinate change.
If $`\eta\lt\sqrt2/24`$, then

```math
\deg_{\mathbb F_2} f\le112.
```

The constant is conservative and independent of dirty width. Individual
axes may mix, the comparisons may be non-Pauli, and B may be non-Clifford.
This is an additional native hypothesis, not a consequence of full return.
The degree-two Pauli-comparison theorem and the earlier conditional bounds
remain stronger in their own subcases. The axis-normalizer bound 244 also
retains its larger error radius $`\sqrt2/12`$.

The feature space below is an analytical invariant. Its dimension may grow
with work width; its **native polynomial degree** does not. No group-algebra
coordinate is supplied to a physical compiler as free workspace or a source.

#### Separating the native baseline from the target

Set $`g(x)=f(x)+f(0)`$. If g is affine the conclusion already holds.
Otherwise work projectively and write the formal comparison as

```math
K(x,t)=A_x^{-1}P_xN_xQ^t,\qquad
K_x=K(x,g(x)),\qquad K_x^0=K(x,0).                    \qquad(31)
```

Here A and N are the two specified native one-T Pauli-image families,
P has linear Pauli labels, and $`Q=W_5^\dagger X_tW_5`$ is the fixed
nonidentity target image. The use of t as a formal Boolean variable in
$`K(x,t)`$ is separate from the physical output-wire subscript in Q.
N and Q share the outer one-T conjugator: their projective labels form
a Pauli representation, and their symplectic increments have a common
isotropic image space $`\mathcal M`$ and multiply pairwise to zero.
All operators are understood in the fixed coordinates of this section.
The zero values are $`A_0=P_0=N_0=I`$ and $`Q^0=I`$;
Q itself remains nonidentity.

Write $`\mathcal L=\mathcal L_D`$. The symplectic A family fixes
$`\mathcal L`$ pointwise, and P has identity symplectic action. The
quotient obstruction to N followed by $`Q^{g(x)}`$ preserving
$`\mathcal L`$ is therefore linear in x plus g times a fixed matrix.
If an entry of that fixed matrix were nonzero, span preservation for all x
would force g linear. Since g is nonaffine, both terms vanish separately.
Thus **every N_x and Q preserves $`\mathcal L`$ individually**.
No span property is inferred at a counterfactual t without this argument.

Let $`\mathsf H`$ be the finite projective Clifford subgroup generated
by all address-generator images of A, P, N, together with Q. It contains
every formal $`K(x,t)`$ and both native baselines. Let
$`\mathsf G=\langle K_x:x\in\mathbb F_2^a\rangle`$. These groups
need not be equal. Only $`\mathsf G`$ is asserted to be normalized
by the actual B family.

#### A native matrix algebra of nilpotency index at most four

Denote A-family symplectic increments by h and N/Q-family increments by n.
For any choices of these increments,

```math
\begin{gathered}
hh'=0,\qquad nn'=0,\qquad
\mathrm{im}\,h\subseteq\mathcal L,\quad h(\mathcal L)=0,\\
\mathrm{im}\,n\subseteq\mathcal M,\quad n(\mathcal M)=0,
\qquad n(\mathcal L)\subseteq\mathcal L\cap\mathcal M.
\end{gathered}                                            \qquad(32)
```

In particular $`hnh'=0`$. A possible nonzero word $`nhn'`$ has image
in $`\mathcal L\cap\mathcal M`$, annihilated on the left by every
h and n. All words of length four vanish. Hence the nonunital binary
matrix algebra $`\mathcal J`$ generated by these increments satisfies

```math
\mathcal J^4=0,\qquad S(\mathsf H)\subseteq1+\mathcal J.
                                                               \qquad(33)
```

For the lower central series $`\gamma_i(\mathsf H)`$, commutators
of algebra units give $`S(\gamma_i(\mathsf H))\subseteq1+\mathcal J^i`$.
Thus $`\gamma_4(\mathsf H)`$ is projective Pauli. A commutator of
a Pauli label v with h has label $`(S_h-I)v`$ up to the immaterial choice
of commutator convention. Subsequent lower-central label spaces lie in
$`\mathcal J V_{\rm label},\mathcal J^2 V_{\rm label},\mathcal J^3 V_{\rm label}`$,
and then zero. Consequently $`\gamma_8(\mathsf H)=1`$.

Also every element of $`\mathsf H`$ has fourth power in the projective
Pauli kernel, while every element of $`\gamma_2(\mathsf H)`$ or
$`\gamma_3(\mathsf H)`$ has square in that kernel. Therefore

```math
\begin{aligned}
\mathrm{exp}\,\mathsf H&\mid8,\\
\mathrm{exp}\,\gamma_2(\mathsf H),\quad
\mathrm{exp}\,\gamma_3(\mathsf H)&\mid4,\\
\mathrm{exp}\,\gamma_i(\mathsf H)&\mid2\quad(4\le i\le7).
\end{aligned}                                             \qquad(34)
```

The finite symplectic quotient lies in the finite 2-group
$`1+\mathcal J`$ and the projective Pauli kernel is elementary abelian.
Thus $`\mathsf H`$ is a finite 2-group.

#### Augmentation features and a subgroup extension lemma

For a finite 2-group F, let $`I_F`$ be the augmentation ideal of its
binary group algebra $`\mathbb F_2[F]`$. We use two classical
Jennings–Lazard facts:

```math
\mathfrak D_r(F)=F\cap(1+I_F^r)
 =\prod_{i2^j\ge r}\gamma_i(F)^{2^j};                    \qquad(35)
```

and the weighted basis obtained by choosing representatives of bases of
$`\mathfrak D_i(F)/\mathfrak D_{i+1}(F)`$. Ordered products of the
corresponding $`(u-1)`$, each used zero or one times, form a basis of the
group algebra; the products of total weight at least r form a basis of
$`I_F^r`$. The weight of u is i. These are imported algebraic facts,
attributed to Jennings and Lazard; an accessible primary reference is
[Usefi, *Identifications in modular group algebras*](https://www.math.mun.ca/~usefi/papers/JPAA.pdf),
the dimension-subgroup formula on page 2, Theorem 3.1 and Corollary 3.4.
Its Lemma 5.1 also treats the subgroup-induced filtration used below.
The weighted basis is explicitly stated in
[Sakurai, Theorem 3.5](https://arxiv.org/pdf/1701.03799).

Equations (34)–(35) imply

```math
\mathfrak D_8(\mathsf H)=1.                              \qquad(36)
```

For example, terms of (35) with i equal to one require an eighth power;
i equal to two or three require a fourth power; i from four through seven
require a square; and all larger lower-central terms are trivial. Every
Jennings generator of $`\mathsf H`$ therefore has weight at most seven.
Equation (36) means that distinct group elements remain distinct in
$`\mathbb F_2[\mathsf H]/I_{\mathsf H}^8`$. It does **not** assert
that the whole augmentation ideal has nilpotency index eight.

We need the following subgroup extension, because the actual normalizer
acts on $`\mathsf G`$ while the native factors lie in $`\mathsf H`$.
If $`G\le H`$, $`\mathfrak D_{w+1}(H)=1`$, and $`m\ge1`$, then

```math
\mathbb F_2[G]\cap I_H^{w(m-1)+1}\subseteq I_G^m.
                                                               \qquad(37)
```

To prove it, at every weight i choose a basis for the image of
$`G\cap\mathfrak D_i(H)`$ in
$`\mathfrak D_i(H)/\mathfrak D_{i+1}(H)`$ and extend it to a full
basis. Choose the first representatives inside G. Descending through the
induced filtration gives unique ordered binary product normal forms for G.
The augmentation monomials in those selected representatives form a basis
of $`\mathbb F_2[G]`$ that is a subset of the full Jennings basis of
$`\mathbb F_2[H]`$. A selected monomial of weight at least
$`w(m-1)+1`$ has at least m factors, since each weight is at most w.
It belongs to $`I_G^m`$. Expansion in this common basis proves (37).

It follows that every linear functional on $`\mathbb F_2[G]/I_G^m`$
extends to a linear functional on
$`\mathbb F_2[H]/I_H^{w(m-1)+1}`$. This is a **linear extension**,
not an asserted algebra embedding between these two truncated quotients.

#### Exact translation action and bounded native degree

By section 5.3, the cocycle is exact projectively. In the original
comparison coordinates it reads

```math
K_{x+y}=\phi_y(K_x)K_y\quad\hbox{in }\mathsf G,
\qquad \phi_y(K)=B_yKB_y^\dagger.                        \qquad(38)
```

Conjugation is an automorphism of $`\mathsf G`$: the cocycle proves
inclusion on its generators, and the projective involution B_y proves
equality. Define the intrinsic feature vector

```math
v(x)=[K_x]\in\mathcal V
 :=\mathbb F_2[\mathsf G]/I_{\mathsf G}^8.               \qquad(39)
```

An automorphism of $`\mathsf G`$ preserves its augmentation powers.
Right multiplication by K_y also preserves the quotient. Thus (38) acts
linearly on $`\mathcal V`$, and its coordinate-function span is
translation invariant. No normalization of the larger $`\mathsf H`$
is needed.

Apply (37) with $`w=7,m=8`$. Every feature functional extends to
$`\mathbb F_2[\mathsf H]/I_{\mathsf H}^{50}`$. Each family in
(31) is a projective Pauli representation, hence an ordered product of
its fixed address-generator images to binary powers. In the group algebra,

```math
u^{x_i}=1+x_i(u-1).
```

Expand the native word $`A_x^{-1}P_xN_xQ^t`$ in these factors.
Every term containing at least 50 augmentation factors vanishes in this
quotient. Each coordinate of v therefore has a Boolean polynomial
representative in x and the fresh t of total degree at most 49, evaluated
at $`t=g(x)`$. Repeated address variables can only reduce multilinear
degree. The counterfactual word need not belong to $`\mathsf G`$;
the linear extension just proved is precisely what permits this calculation.

#### Detecting the target and bounding its degree

The natural group-algebra map induces a linear map

```math
\mathcal V\longrightarrow
\mathbb F_2[\mathsf H]/I_{\mathsf H}^8,
```

because $`I_{\mathsf G}^8`$ maps into $`I_{\mathsf H}^8`$.
In that ambient quotient,

```math
[(K_x^0)^{-1}][K_x]
 =[Q]^{g(x)}=1+g(x)([Q]-1).                             \qquad(40)
```

The native inverse baseline has polynomial degree at most seven in x in
this quotient, by the same augmentation expansion, with every inverse
in its actual reversed order. The target Q is nonidentity projectively,
so (36) gives $`[Q]-1\ne0`$. Choose a binary linear functional taking
value one on this difference. Applying it to (40) plus 1 extracts g.

Let $`\mathcal W`$ be the span of degree-at-most-seven x monomials
times feature coordinates from (39). It includes constants, is translation
invariant, and contains g. The native degree calculation gives

```math
g\in\mathcal W\subseteq\mathcal P_{56}+g\mathcal P_{55}.
                                                               \qquad(41)
```

If $`d=\deg g`$, restrict to the d variables of a highest-degree
monomial. Its translates span all $`2^d`$ Boolean functions, whereas
the right side of (41) has dimension at most

```math
\sum_{j=0}^{56}\binom dj+\sum_{j=0}^{55}\binom dj,
```

strictly less than $`2^d`$ for $`d\ge113`$. This is the same
translation-space argument used in section 5.2. Hence $`\deg g\le112`$,
which proves the theorem.

This necessary invariant adds no T gate, elementary Clifford gate, source
call, inverse call or workspace. T-count, elementary Clifford count and
finite dirty width remain unrestricted within the lower-bound statement.
Its fixed error radius comes entirely from the original full-operator
cocycle reduction, including every logical column, literal signs and phases,
arbitrary dirty and reference inputs, leakage and complete work return.
The proof uses projective groups only after that reduction; it does not
replace the original literal target by a phase-relaxed specification.

The [bounded evidence](../../tests/test_span_preserving_augmentation.py)
checks a native mixed-axis shear with nonzero cubic products, the fourth-power
nilpotence mechanism, the distinction between subgroup and ambient
augmentation filtrations, a nonabelian finite group and the degree threshold.
These checks support the analytic proof and do not extrapolate from finite
width. General comparisons need not preserve $`\mathcal L_D`$;
the general five-layer and complete-frame depth questions remain open.

### 5.5 Shared-prefix comparisons: a robust quadratic group map

There is an unconditional structural restriction on the **actual**
comparison family, even when it does not preserve the outer axis span.
It uses the common native prefix of B and J, which the independently
chosen baselines in the
[baseline-word counterexample](BASELINE_WORD_OBSTRUCTION.md) do not retain.

**Shared-prefix theorem.** In the native setting of section 1, suppose

```math
\eta\lt\frac{\sqrt2}{40}.
```

Then the projective group $`\mathsf G=\langle K_x\rangle`$ has
nilpotency class at most two, exponent dividing four, and commutator
subgroup of exponent dividing two. In particular,
$`\mathfrak D_3(\mathsf G)=1`$. The comparison family is a
quadratic group-valued function of the address, in the explicit sense
proved below. Its intrinsic features in
$`\mathbb F_2[\mathsf G]/I_{\mathsf G}^3`$ distinguish group
elements and have Boolean polynomial degree at most four in x.

This does not yet bound the Boolean output f. The formal outer baseline
and target Q need not lie in this group, or in a common ambient group
with a bounded dimension-subgroup filtration. Section 5.4 establishes
that additional bridge under span preservation; it is missing in general.

#### The common prefix places the exact native products in a small-class group

Retain every fixed Clifford and actual inverse, and define

```math
\mathcal F=C_4(0)D_4C_3(0)D_3,\qquad E=C_2(0)D_2,
```

```math
Y_x=EP_1(x)E^\dagger,\qquad
H_x=P_2(x)Y_x,\qquad
\widehat K_x=B_xJ_x=\mathcal F H_x\mathcal F^\dagger.
                                                               \qquad(42)
```

The equality is an exact identity of the specified native operators;
no independent replacement source is introduced. Equation (3) gives
$`\|K_x-\widehat K_x\|\le2\eta`$ on the full work space.

Let $`\mathsf H_{\rm in}`$ be the projective group generated by
the P_2 and Y address-generator images. Each of these two families
commutes internally projectively, since it is a fixed-unitary image
of a linear-label Pauli representation. Each generator is an involution
projectively. A cross commutator of a Pauli P with Y is again a Pauli,
whose binary label is $`(S_Y-I)p`$. This label lies in the fixed
isotropic native-axis span

```math
\mathcal M=C_2(0)\mathrm{span}_{\mathbb F_2}
 \{Z_j:D_2\text{ is active on }j\}.
```

Here C_2 acts on binary Pauli labels. Every Y shear fixes M pointwise,
and every Pauli has trivial symplectic action. Thus all cross
commutators are central projectively. The two families generate a
class-at-most-two group whose commutators are Pauli involutions.
Products have fourth power one projectively. Therefore
$`\mathsf H_{\rm in}`$ has exponent dividing four.
Conjugation by the possibly non-Clifford F preserves these intrinsic
group statements exactly.

For exact return, (42) already proves the group assertion for G.
In particular, no comparison value, nor any product of comparison
values, can have nontrivial odd projective order. The earlier independent
baseline source with arbitrary odd symplectic order cannot be an actual
five-layer returned comparison family. This is an obstruction to that
proposed counterexample, not a new restriction on freely chosen baselines.

#### A bounded word freezes the quadratic address law under full error

Work projectively for the group identities in this paragraph. For the
native pullbacks, define the ordered defect

```math
Z_H(x,y)=H_{x+y}H_y^{-1}H_x^{-1}.
```

Using the two exact projective representation laws and the centrality
just proved gives

```math
Z_H(x,y)=[Y_x,P_2(y)].                                  \qquad(43)
```

The choice of inverse commutator convention does not affect (43): the
commutator is a projective involution. It is central, has square one,
and is multiplicative under address addition in each variable.
Conjugating by F gives the same statements for the native products
$`\widehat K`$.

Now define the actual Clifford defect

```math
Z(x,y)=K_{x+y}K_y^{-1}K_x^{-1}.
```

The following relations are words in actual K values and their actual
inverses. Their lengths count every occurrence before cancellation:

| Relation | Word length |
|---|---:|
| $`[Z(x,y),K_z]=1`$ | 8 |
| $`Z(x,y)^2=1`$ | 6 |
| $`Z(x+y,z)=Z(x,z)Z(y,z)`$ | 9 |
| $`Z(x,y+z)=Z(x,y)Z(x,z)`$ | 9 |

Replacing each K by its native product changes a length-l word by at
most $`2l\eta`$ in full operator norm. The native word is scalar
by (43). The actual word is Clifford up to scalar. The projective
Clifford separation gap $`1/\sqrt2`$ therefore freezes every
displayed relation exactly projectively at the stated error radius.
No arbitrary approximate normalizer is rounded, and no error grows
with dirty width or address dimension.

Since $`K_0=I`$, the diagonal relation
$`Z(x,x)=K_x^{-2}`$ gives $`K_x^4=1`$ projectively.
Interchanging x and y in the defining defect relation shows that
$`[K_x,K_y]`$ is a product of two central involutions. Hence
the generators have fourth power one, their commutators are central
involutions, and G has class at most two and exponent dividing four.
It is finite because it is a subgroup of the finite projective
Clifford group on the given finite work space. It is a 2-group.
The dimension-subgroup formula (35) now gives

```math
\mathfrak D_3(\mathsf G)
 =\mathsf G^4\,\gamma_2(\mathsf G)^2\,\gamma_3(\mathsf G)
 =1.                                                    \qquad(44)
```

Thus intrinsic filtration depth is bounded even without an outer-span
hypothesis. This statement does not mean that the whole augmentation
ideal has nilpotency index three.

#### A degree-four intrinsic address feature, and its missing target detector

Let e_i denote the address coordinate vectors, with an arbitrary fixed
ordering. The central, bilinear defect law gives the ordered normal form

```math
K_x=
\left(\prod_i K_{e_i}^{x_i}\right)
\left(\prod_{i\lt j}Z(e_i,e_j)^{x_i x_j}\right)
\quad\text{in }\mathsf G.                              \qquad(45)
```

In the group algebra each binary power has the expansion
$`u^b=1+b(u-1)`$. Modulo $`I_{\mathsf G}^3`$, terms with
three augmentation factors vanish. Each retained factor contributes
an address monomial of degree at most two in (45). Therefore every
coordinate of the faithful intrinsic feature

```math
v(x)=[K_x]\in\mathbb F_2[\mathsf G]/I_{\mathsf G}^3
```

has Boolean degree at most four. This conclusion uses the actual
comparison group; there is no unjustified extension to an independent
outer baseline. It also does not make the symplectic tableau entries
linear functions of these intrinsic coordinates.

The same expansion shows, for every integer m at least one, that every
coordinate of $`[K_x]`$ modulo $`I_{\mathsf G}^m`$ has degree at
most $`2(m-1)`$. Only terms with fewer than m augmentation factors
survive. No faithfulness claim is needed for this degree calculation.

The remaining general question is now sharper: can one detect the fixed
target factor in $`K_x=(K_x^0)Q^{g(x)}`$ from these bounded-degree
intrinsic features together with bounded-degree native baseline data?
One must construct that detector, or a bounded ambient extension for
the actual baseline and Q, while retaining their shared sources. The
independent-baseline odd-order example rules out assuming such an ambient
extension for arbitrary pairs. Bounded intrinsic filtration of G itself
is no longer the missing invariant.

**Smaller-error span-preserving corollary.** Combining this degree
calculation with section 5.4 gives

```math
\deg f\le21
\quad\text{if K preserves the outer axis span and }
\eta\lt\sqrt2/40.
```

Indeed, take m equal to eight. The intrinsic features have degree at
most 14. Their natural linear map to the span-preserving ambient quotient
$`\mathbb F_2[\mathsf H]/I_{\mathsf H}^8`$ is well defined,
because $`I_{\mathsf G}^8`$ maps into $`I_{\mathsf H}^8`$.
The actual inverse baseline has degree at most seven there, and the
nonzero target difference $`[Q]-1`$ is detected as in (40).
The detector therefore expresses g as a polynomial of degree at most
14 plus seven. The degree-112 bound remains valid at its larger radius
$`\sqrt2/24`$; neither constant is claimed sharp.

This proof adds no native gate, source, control, inverse or workspace.
T-count, elementary Clifford count and finite dirty width remain
unrestricted in this structural lower-bound statement. The error refers
to the original literal full-operator query, including arbitrary output,
dirty and external-reference inputs, actual inverses, leakage and complete
workspace return. Projective group analysis does not relax that contract.
The [bounded evidence](../../tests/test_shared_prefix_comparison.py)
checks the native pullback, noncommutative group law, ordered defects,
robust word lengths and faithful intrinsic features. No general Boolean
degree bound or complete-frame improvement is claimed here.

## 6. Full return does not force a Clifford middle

A literal native affine query already falsifies that implication. On
one arbitrary work qubit put

```math
U=THT,\qquad K=\frac{X-Y}{\sqrt2},\qquad
B=UXU^\dagger=\frac{Z+K}{\sqrt2},\qquad
J=UYU^\dagger=\frac{K-Z}{\sqrt2}.
```

Both B and J are two-T Pauli images and are non-Clifford. Direct
multiplication gives

```math
BZB=K,\qquad BKB=Z,\qquad BJ=ZK=:D.
```

The unitary K and product D are Clifford. The family $`D_x=D^x`$
for one address bit obeys the exact cocycle

```math
D_{x+y}=B^yD_xB^yD_y.
```

Indeed $`BDB=-D`$ and $`D^2=-I`$. The same B also normalizes
the irreducible Clifford subgroup generated by $`iI,Z,K`$;
normalizing an irreducible Clifford subgroup therefore does not force
a native two-T Pauli image to be Clifford.

This example occurs in an actual returned query, not only an abstract
normalizer equation. The complete fibre word is

```math
V_x=T Z^x T H T X^x Y^x T^\dagger H S^\dagger=X^x.
                                                               \qquad(17)
```

At x=0 it is $`SS^\dagger=I`$. At x=1, the identity
$`XY=iZ`$ reduces the word to
$`iZSXS^\dagger=iZY=X`$, including its scalar. In the six-band
notation its nonconstant Clifford bands are
$`C_4=Z^x,C_2=X^x,C_1=Y^x`$, with
$`C_3=H,C_0=HS^\dagger,C_5=I`$. The physical layers are
$`D_5=D_4=D_3=T,D_2=I,D_1=T^\dagger`$.
The middle family in section 1 is exactly $`B^x`$.
If $`Y^x`$ is rewritten in the ordered XZ section, its literal
factor $`i^x`$ must stay in lambda.

A direct native expansion uses four T gates, four T layers and 14
elementary Cliffords, with zero clean qubits and zero additional dirty
helpers. The output is arbitrary. Controlled-Z uses two H and one
CNOT; controlled-Y uses S, CNOT and $`S^\dagger=S^3`$; the
remaining fixed H and S inverses are included in that count. Equation
(17) proves zero error, zero leakage and full return with external
references. Its fair simplified baseline is one CNOT, with zero T
count and depth and one elementary Clifford. The fixture supplies no
resource saving; it proves that full return alone cannot justify the
extra hypothesis in section 2. Its selected middle representation also
changes under resynthesis.

## 7. The missing general invariant and scope

The native reduction retains the actual middle action on a particular
Clifford comparison group. Sections 5.1 and 5.2 resolve Pauli comparisons
and individual-axis permutations. Sections 5.3–5.4 resolve the broader
span-preserving branch at a fixed full-error radius, by an intrinsic
augmentation feature space and a bounded ambient extension.

Section 5.5 now supplies a structural invariant **without** span preservation:
the actual comparison group is class at most two of exponent dividing four,
and its address map has faithful degree-four intrinsic features. This holds
below error $`\sqrt2/40`$. Unitary conjugation still need not act with
bounded Boolean degree on ordinary Clifford tableaus, and faithful intrinsic
features do not make every operator observable a low-degree feature functional.

The remaining step is to detect g in $`K_x=K_x^0Q^{g(x)}`$ with bounded-degree
native baseline data. The baseline and Q need not lie in the actual group.
An arbitrary ambient group generated by independent one-T baselines can have
odd order, so the span-preserving proof's ambient filtration cannot simply be
assumed. A general detector or compatibility identity for the actual shared
sources is required before a five-layer Boolean degree theorem, and much more
before a growing-depth recurrence.

The conditional theorem retains every logical column and literal target
phase, arbitrary output and dirty inputs with external references, the
actual inverse, leakage and complete workspace return. It does not
provide a new native query compiler, a smaller complete-frame upper
bound, unrestricted one-clean rigidity, or a dimension-growing general
frame lower bound. R96, selected Result D R79, and the conditional status
of R98 are unchanged.

The four bounded [five-layer checks](../../tests/test_read_only_five_layer_reduction.py)
verify the literal six-band normal form and actual inverse, the native
returned-query counterexample and its ledger, the full-degree translation
rank against the dimension bound, and the Hermitian Clifford sign formulas.
These finite identity checks support the analytic proof; they do not
establish an unrestricted degree law by enumeration.

The five [Pauli-comparison checks](../../tests/test_pauli_comparison_cubic.py)
verify native-axis label cancellation, the quadratic-baseline cube product,
the peeled Clifford comparison, a non-Clifford middle's exact partial
Pauli action, and cubic contraction rank. The cube-error fixture also
retains the weaker direct product estimate from the initial derivation;
the theorem uses the stronger native separation step. These finite checks
do not enumerate or prove the unrestricted five-layer case.
