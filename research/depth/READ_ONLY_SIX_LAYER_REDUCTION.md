# Six-layer native reduction and a Pauli-comparison degree bound

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[General five-layer theorem](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound) ·
[Six-layer boundary](READ_ONLY_FIVE_LAYER_REDUCTION.md#71-the-six-layer-native-boundary)

The general zero-clean, read-only six-layer degree question remains open.
This note proves a width-independent bound for a specified **Pauli-valued
target comparison**. The proof does not assume that the three-band prefix
or its intermediate conjugations are Clifford. The extra Pauli-comparison
hypothesis is stated on the displayed native circuit and is not inferred
from full workspace return. General comparisons need not form a nilpotent
group, as the returned fixture below demonstrates.

**Six-layer Pauli-comparison theorem.** In the native setting below,
suppose every specified comparison $`K_x`$ is a scalar times a Pauli.
If

```math
\eta\lt\frac{\sqrt2}{64},
```

then $`\deg_{\mathbb F_2}f\le51`$, independently of dirty width.
The constants are conservative; sharpness and attainability in this subcase
are not asserted. The proof freezes an eight-factor Clifford product. It
does not promote an arbitrary intermediate of a three-band conjugation to
a Clifford.

## 1. Native contract and the phase-retaining comparison

Let x be an a-bit address and W contain the arbitrary output t and all
arbitrary dirty helpers. There are zero initialized work qubits. Write

```math
V=C_6D_6C_5D_5\cdots C_1D_1C_0,\qquad
F_x=X_t^{f(x)}\otimes I_{\rm dirty}.
```

Each completed Clifford interlayer preserves the original computational
address pointwise. Each D is one physical band of T, T-dagger and identity
on distinct work wires; empty bands are allowed. Assume the complete
operator approximation

```math
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}\le\eta.
```

This norm retains every address, output and dirty input, arbitrary external
references, literal target phases, leakage and complete workspace return.
T-count, elementary Clifford count and finite dirty allocation are
unrestricted in this necessary condition. There are no measurements,
postselected inputs or initialized-only source promises.

Use the retained address-preserving Clifford normal form

```math
C_j(x)=\lambda_j(x)P_j(x)C_j(0),\qquad
P_j(0)=I,\quad\lambda_j(0)=1.
```

The binary Pauli labels of P are linear in x. Their actual projective
phases, including all address S and CZ phases, stay in the normal form.
Put

```math
\lambda(x)=\prod_{j=0}^6\lambda_j(x),\qquad
W_6=C_6(0)D_6,\qquad
R_x=C_0(0)^\dagger P_0(x)C_0(0).
```

Let A, B, J and O be the successive exact native Pauli-image families:
A conjugates $`P_4(x)`$ by $`C_5(0)D_5`$; B conjugates
$`P_3(x)`$ by $`C_5(0)D_5C_4(0)D_4`$; J adds
$`C_3(0)D_3`$ to that prefix and conjugates $`P_2(x)`$; O adds
$`C_2(0)D_2`$ and conjugates $`P_1(x)`$. All inverses in these
definitions are the actual adjoints. In particular,

```math
A_x=C_5(0)D_5P_4(x)D_5^\dagger C_5(0)^\dagger.
```

For the actual work fibre $`U_x=V(x)`$, direct multiplication gives

```math
U_x=\lambda(x)P_6(x)W_6P_5(x)A_xB_xJ_xO_x
                  W_6^\dagger U_0R_x.                    \qquad(1)
```

Define the specified phase-retaining comparisons

```math
\begin{aligned}
M_x&=P_5(x)^\dagger W_6^\dagger
 [\lambda(x)^{-1}P_6(x)^\dagger F_xR_x^\dagger F_0^\dagger]W_6,\\
K_x&=A_x^\dagger M_x.
\end{aligned}                                             \qquad(2)
```

Each M and K is a scalar times a Clifford: the bracket is a scalar
times a Pauli, and its passage through $`W_6`$ uses only one physical
T band. Solving (1) with the actual inverse $`U_0^\dagger`$, then
replacing $`U_x`$ and $`U_0^\dagger`$ by their target approximations,
gives

```math
\|B_xJ_xO_x-K_x\|\le2\eta.                              \qquad(3)
```

This is the six-layer version of the
[phase-retaining five-layer reduction](READ_ONLY_FIVE_LAYER_REDUCTION.md#1-native-contract-and-an-exact-one-band-reduction).
No source call or inverse has been omitted from V; (1)--(3) are analytical
identities and estimates for that complete circuit.

## 2. A Clifford-valued native pullback

Set

```math
\begin{aligned}
S&=C_3(0)D_3,\qquad E_{\rm in}=C_2(0)D_2,\\
\mathcal F&=C_5(0)D_5C_4(0)D_4S,\\
L_x&=(S^\dagger P_3(x)S)\,P_2(x)\,
                   (E_{\rm in}P_1(x)E_{\rm in}^\dagger).
\end{aligned}
```

Then the shared native prefixes cancel exactly:

```math
B_xJ_xO_x=\mathcal F L_x\mathcal F^\dagger,\qquad
\|L_x-\mathcal F^\dagger K_x\mathcal F\|\le2\eta.        \qquad(4)
```

L is Clifford-valued up to its actual scalar. Its first factor is
$`D_3^\dagger[C_3(0)^\dagger P_3(x)C_3(0)]D_3`$, a one-band
Pauli image with the inverse T orientation retained. Its last factor is
also a one-band Pauli image, with potentially different native axes.
Neither $`\mathcal F`$ nor its action on arbitrary Cliffords is asserted
Clifford-valued.

The [existing native checks](../../tests/test_six_layer_defect_boundary.py)
verify (1)--(4) and actual inverse order. They also exhibit a non-Pauli
ordered defect of L. Thus a Pauli defect cannot be assumed in the general
six-layer problem.

## 3. The Pauli-comparison label and its third difference

Put $`g(x)=f(x)+f(0)`$. If g has degree at most two, the theorem already
holds. Otherwise the Pauli-comparison hypothesis implies

```math
m(x)=r(x)+v g(x),\qquad \deg r\le2,\qquad v\ne0,          \qquad(5)
```

where m is the binary Pauli label of K and v is fixed.
This is precisely the outer-factor argument in the retained
[Pauli-comparison proof](READ_ONLY_FIVE_LAYER_REDUCTION.md#51-a-pauli-valued-comparison-robust-quadratic-bound),
through its equation (18). It does not use the number or nature of the
inner factors B, J or O: the outer comparison has the same two one-band
factors. Briefly, its trivial symplectic action forces the outer target
shear to vanish for nonaffine g and makes the two native shears equal.
Their Pauli-valued quotient has a linear label by the independent-axis
cocycle. The remaining conjugated linear Pauli label is quadratic, giving
(5). The target image is a nonidentity Pauli, so v is nonzero.

Fix any three address directions u, v_1 and v_2, allowing dependence and
zero directions. In any fixed ordering of the eight binary triples s,
write

```math
\begin{aligned}
x_s&=x+s_1u+s_2v_1+s_3v_2,\\
e(x)&=\delta_u\delta_{v_1}\delta_{v_2}g(x),\\
\mathcal C(x)&=\prod_{s\in\{0,1\}^3}L_{x_s}.
\end{aligned}                                             \qquad(6)
```

Here $`\delta_a g(x)=g(x+a)+g(x)`$ is the Boolean difference.
The same ordered product of the K operators has Pauli label
$`\sum_s m(x_s)=v e(x)`$, since the third difference of r is zero.
Keep its literal scalar as $`\omega(x)`$ and put

```math
\mathcal R=\mathcal F^\dagger\tau(v)\mathcal F.
```

Eight uses of (4), telescoping the actual unitary products, give

```math
\|\mathcal C(x)-\omega(x)\mathcal R^{e(x)}\|
\le16\eta.                                               \qquad(7)
```

This estimate is established before passing to projective operators.
The scalar includes every comparison phase and Pauli multiplication sign;
it is not set to one. No commutation of approximate factors is used.

## 4. Native tableau degree of the eight-factor product

Use the Hermitian Pauli section and its literal conjugation signs from the
[fixed-Clifford sign lemma](READ_ONLY_FIVE_LAYER_REDUCTION.md#4-the-sign-of-a-fixed-clifford-on-hermitian-paulis).
The symplectic entries of L have degree at most two in x, since its two
one-band Pauli images each have affine symplectic action. Therefore

```math
\deg_x S_{\mathcal C}(x)\le16.                            \qquad(8)
```

A moving-probe bound also follows from the same native factors. If the
incoming Pauli label has degree at most d, passage through L has outgoing
label degree at most $`d+2`$ and sign degree at most $`3d+6`$.
The following conservative ledger allows fixed Clifford coordinates on
both one-band factors; an identity coordinate change is permitted.

| Passage through L, right to left | Sign degree | Outgoing label degree |
|---|---:|---:|
| Fixed coordinates into the last one-band image | $`3d`$ | d |
| Its physical one-band Pauli image | $`2d+1`$ | $`d+1`$ |
| Fixed coordinates back | $`3d+3`$ | $`d+1`$ |
| The middle linear-label Pauli | $`d+2`$ | $`d+1`$ |
| Fixed coordinates into the first one-band image | $`3d+3`$ | $`d+1`$ |
| Its physical one-band Pauli image | $`2d+3`$ | $`d+2`$ |
| Fixed coordinates back | $`3d+6`$ | $`d+2`$ |

For a physical one-band image, the local sign is the retained expression
$`u b(1+a)+w a+\epsilon u a`$, with affine family labels u,w and
incoming probe labels a,b. Fixed Clifford signs are cubic in the probe
label. These facts justify every row, including inactive wires and either
T orientation. Phases of the one-band source itself cancel only in this
conjugation action; its actual adjoint remains in the operator definition.

Shifting x to any $`x_s`$ preserves these degree bounds. Starting with a
fixed probe, its label degree before the j-th factor is at most
$`2(j-1)`$, so that factor's sign degree is at most $`6j`$.
Adding the eight Boolean signs consequently gives

```math
\deg_x\kappa_{\mathcal C(x)}(q)\le48
\quad\text{for every fixed Hermitian Pauli probe }q.       \qquad(9)
```

These are polynomial degrees of the exact native Clifford family L and
its product. They do not assume a polynomial action of an arbitrary
unitary on Clifford tableaus.

## 5. One common frozen Clifford detects every third difference

Write $`d_{\rm pr}(U,V)=\inf_{|\alpha|=1}\|U-\alpha V\|`$.
The retained dimension-independent Clifford separation says that
projective Cliffords at distance strictly below $`1/\sqrt2`$ are
projectively equal.

If every third difference of g is zero, g has degree at most two and
there is nothing further to prove. Otherwise choose **one** triple of
directions and one base address with e equal to one, and take its exact
Clifford product $`\mathcal C`$ as $`\mathcal R'`$. This single
reference is used for all subsequent triples and all base addresses.
Equation (7) gives

```math
d_{\rm pr}(\mathcal R',\mathcal R)\le16\eta.
```

The operator $`\mathcal R`$ is a conjugate of a nonidentity Hermitian
Pauli and has eigenvalues plus and minus one. Hence
$`d_{\rm pr}(\mathcal R,I)=\sqrt2`$. At the stated error radius,
$`\mathcal R'`$ cannot be projectively identity.

For any cube (6), a zero value of e gives projective distance at most
$`16\eta`$ from I. A value one gives distance at most $`32\eta`$
from $`\mathcal R'`$, by the triangle inequality. Since

```math
32\eta\lt\frac1{\sqrt2},
```

Clifford separation freezes all cubes:

```math
\mathcal C(x)\doteq
\begin{cases}
I,&e(x)=0,\\
\mathcal R',&e(x)=1.
\end{cases}                                               \qquad(10)
```

No claim is made that $`\mathcal R`$ itself is Clifford. Only the exact
native Clifford reference $`\mathcal R'`$ and the other exact Clifford
products are compared using their separation gap.

If $`\mathcal R'`$ is non-Pauli, a nonzero entry of
$`S_{\mathcal R'}-I`$ detects e in (10). Equation (8) then gives
$`\deg e\le16`$ for every triple, and hence $`\deg g\le19`$.
If $`\mathcal R'`$ is Pauli up to scalar, its nonzero label admits a
fixed anticommuting Hermitian Pauli probe. Conjugating this probe in (10)
gives the literal sign $`(-1)^{e(x)}`$. Equation (9) therefore gives
$`\deg e\le48`$ for every triple.

Finally, a monomial of g of degree at least 52 would have a derivative
in three of its coordinate directions of degree at least 49. This
contradicts the last bound. Thus $`\deg g\le51`$, and adding the
constant $`f(0)`$ proves the theorem for f. The non-Pauli reference
branch gives the smaller bound 19, without asserting that either branch
is attained by a returned circuit.

## 6. A conditional defect detector for the general comparison

The remaining target-extraction step can be isolated without a nilpotence
assumption. For general Clifford-valued K define its actual ordered defect

```math
Z_K(x,y)=K_{x+y}K_y^{-1}K_x^{-1},\qquad
K_{x+y}=Z_K(x,y)K_xK_y.                                  \qquad(11)
```

Suppose that, for every fixed y, its symplectic entries have degree at most
d in x. Suppose also that its literal conjugation sign on any Pauli probe
whose label has degree at most two in x has degree at most s. These must
be bounds for the actual defect family, established without assuming a
bound on the unknown g. Then the outer comparison alone implies

```math
\begin{aligned}
\deg g&\le2d+4 &&\text{if the outer target shear is visible},\\
\deg g&\le\max(s,6)+1 &&\text{if that shear vanishes}.
\end{aligned}                                             \qquad(12)
```

Indeed, in the visible case the native symplectic expression for K has
formal total degree two in x and a new Boolean bit replacing g(x).
Equation (11) and the retained fixed linear target extractor therefore give

```math
g(x+y)\in\mathcal P_{d+2}+g\mathcal P_{d+1}.
```

The translation-dimension argument in the
[five-layer proof](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
gives the first bound in (12). In the invisible case, K has symplectic
degree at most two and its sign on a fixed probe q is
$`p_q(x)+\langle v,q\rangle g(x)`$, with $`\deg p_q\le6`$.
Apply (11) to a probe detecting v. The rightmost $`K_y`$ supplies a
constant probe label; $`K_x`$ supplies a label of degree at most two;
the assumed sign bound for $`Z_K`$ then gives
$`g(x+y)=\epsilon_y g(x)+h_y(x)`$ with
$`\deg h_y\le\max(s,6)`$. If g has higher degree, every epsilon is
one, bounding all first derivatives and proving the second assertion.

This is a conditional algebraic detector, not a new native bound on d or s.
It adds no approximation premise beyond whatever establishes those actual
defect bounds, and introduces no gate or workspace.

The six-layer input to that missing step is explicit. Set

```math
Z_L(x,y)=L_{x+y}L_y^{-1}L_x^{-1}.
```

For fixed y its symplectic degree is at most four. Its sign on a degree-D
moving probe has degree at most $`3D+12`$: the right variable factor
increases label degree by two, the middle constant Clifford contributes
a cubic sign on that label, and the left variable factor adds the same
native ledger again. In particular the sign degree is at most eighteen
for a degree-two probe. These bounds also hold for the actual inverse
factors, whose native one-band factors occur in reversed order.

Three uses of (4) give the full-operator estimate

```math
\|Z_K(x,y)-\mathcal F Z_L(x,y)\mathcal F^\dagger\|
\le6\eta.                                                \qquad(13)
```

Thus a general six-layer theorem needs controlled amplification of these
native degree-four/eighteen data through the specified three-band prefix.
The source defect is Clifford but can be non-Pauli. Its intermediate
conjugations need not be Clifford, so one-band matching cannot simply be
iterated. Neither (12) nor (13) supplies the missing degree transfer or its
width-independent error radius. A depth recurrence must track that radius
as well as both polynomial degrees.

## 7. Resource meaning and the remaining general question

All eight-factor products and reference selections above are analytical
comparisons. The proof adds no physical gate, source call, query, control,
inverse, clean qubit or dirty helper. The actual circuit retains the six
physical T bands and all its T and elementary Clifford gates; no source
or inverse is supplied as an uncharged construction. The bound is
independent of its finite simultaneous work allocation and T-count.

The argument retains the original full operator error. Projective scalars
are discarded only after the literal telescoping estimate (7), when using
Clifford separation and conjugation of Hermitian Pauli probes. Literal
query signs, phases, arbitrary dirty/reference inputs and workspace return
have not been weakened to an initialized-column promise.

The general six-layer comparison K can be a non-Pauli Clifford. For that
case the quadratic binary label (5) is unavailable, and an eight-factor
operator cube need not cancel the formal baseline. The
[baseline-word obstruction](BASELINE_WORD_OBSTRUCTION.md) prevents
silently replacing that missing argument by a universal ordering of
independent baseline factors. The next question remains a bounded native
action or target detector for the actual shared three-band source and
non-Pauli comparisons, with a controlled error radius.

This result gives no unrestricted one-clean or address-mixing rigidity,
no dimension-growing general-frame lower bound, and no improvement to
R96's complete-circuit upper law. R79 remains selected Result D and R97
keeps its stated hypotheses. The
[focused evidence](../../tests/test_six_layer_pauli_comparison.py) checks
the native sign degrees, the quadratic-label cube, reference separation
and threshold accounting; the analytic argument proves the arbitrary-width
statement.

## 8. An actual returned comparison with an odd-order defect

Full return does not force the specified six-layer comparison group to be
a 2-group. The following is a complete native counterexample to that
candidate invariant. It also prevents an unconditional recurrence that
repeatedly differentiates multiplicative defects until they become Pauli.
It is deliberately redundant and supplies no circuit saving.

On one arbitrary work/output qubit with two untouched address bits u,v,
take the literal fibre

```math
V_{u,v}=(-1)^uX^vT^\dagger HTY^uS^\dagger X^uTHTX^v=I.   \qquad(14)
```

Indeed, $`TS^\dagger T=I`$, whereas

```math
TYS^\dagger XT=(TYT^\dagger)(T^\dagger XT)=-I.
```

The explicit address phase $`(-1)^u`$ cancels this minus sign. Thus (14)
is identity on every address and work column, including their arbitrary
superpositions and external references. Its error and leakage are zero,
and it returns the complete workspace without an initialized-state premise.

Place this word into the six native slots as follows; the order remains
$`C_6D_6\cdots C_1D_1C_0`$.

| j | $`C_j(u,v)`$ | $`D_j`$ |
|---|---|---|
| 6 | $`(-1)^uX^v`$ | $`T^\dagger`$ |
| 5 | H | T |
| 4 | $`Y^uS^\dagger`$ | I |
| 3 | $`X^u`$ | T |
| 2 | H | T |
| 1 | $`X^v`$ | I |
| 0 | I | — |

Every completed Clifford interlayer preserves the address pointwise.
Choose the Pauli representatives $`P_4=Y^u`$, $`P_3=X^u`$ and
$`P_1=P_6=X^v`$. Then the retained scalar is
$`\lambda(u,v)=(-1)^u`$. Define

```math
\mathsf A=T^\dagger XT,\qquad
\mathsf B=HTXT^\dagger H.
```

The exact prefix and inner family from section 2 become

```math
\mathcal F=HTS^\dagger T=H,\qquad
L_{u,v}=\mathsf A^u\mathsf B^v.
```

The prefix is Clifford after fair simplification; its nominal bands are
not claimed irreducible. Merging empty T slots or resynthesizing the word
changes the defined comparison K; this witness with four nonempty bands
therefore does not contradict the five-layer comparison theorem. The outer
family and target comparison in (2),
for the identity target f=0, are literally

```math
\begin{aligned}
A_{u,v}&=HTY^uT^\dagger H=(-1)^uH\mathsf A^uH,\\
M_{u,v}&=(-1)^uTX^vT^\dagger,\\
K_{u,v}&=A_{u,v}^\dagger M_{u,v}
        =H\mathsf A^u\mathsf B^vH
        =\mathcal F L_{u,v}\mathcal F^\dagger.
\end{aligned}                                             \qquad(15)
```

Thus K is the actual comparison defined by the full returned circuit,
including its literal scalar cancellation. It is not an independently
assigned Clifford family. In binary X,Z coordinate order its tableaus are

```math
S_{K_{1,0}}=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad
S_{K_{0,1}}=\begin{pmatrix}1&0\\1&1\end{pmatrix},\qquad
S_{K_{1,1}}=\begin{pmatrix}0&1\\1&1\end{pmatrix}.
```

The last matrix has order three. For the literal phases, direct Pauli
multiplication gives

```math
\mathsf A=\frac{X-Y}{\sqrt2},\qquad
\mathsf B=\frac{Z-Y}{\sqrt2},\qquad
\mathsf A\mathsf B=\frac{I-i(X+Y+Z)}2.
```

Since $`(X+Y+Z)^2=3I`$, this implies
$`(\mathsf A\mathsf B)^3=-I`$ and hence $`K_{1,1}^3=-I`$.
For $`x=(0,1)`$ and $`y=(1,0)`$ the ordered defect of section 6 is

```math
Z_K(x,y)=K_{1,1}K_{1,0}^{-1}K_{0,1}^{-1}
        =H\mathsf C H,\qquad
\mathsf C=\mathsf A\mathsf B\mathsf A\mathsf B,
\qquad \mathsf C^3=I.                                    \qquad(16)
```

Its symplectic matrix is $`\left(\begin{smallmatrix}1&1\\1&0\end{smallmatrix}\right)`$,
also of order three. The actual comparison group is consequently not a
2-group, even projectively. The general five-layer class-two/exponent-four
conclusion cannot be carried over to these six-layer comparisons.

### 8.1 Repeated multiplicative differences retain the odd-order mode

Fix $`y=(1,0)`$ and restrict $`x=(0,u)`$. The inner defect is
$`W(u)=\mathsf C^u`$. For the ordinary multiplicative difference

```math
\Delta W(u)=W(u+1)W(u)^\dagger
```

the values are $`(\mathsf C,\mathsf C^{-1})`$. Differentiating that pair
again leaves it unchanged, since $`\mathsf C^{-2}=\mathsf C`$ and
$`\mathsf C^2=\mathsf C^{-1}`$. Every positive iterated difference in
this direction is therefore non-Pauli. The fixed outer conjugation by H
preserves the same conclusion for the actual comparison defects.

This refutes that specific universal Pauli-return recurrence. It does not
rule out a richer invariant, a different native representation, or a
general depth lower bound. Any successful detector for the displayed
general six-layer comparison must accommodate these odd-order modes;
full workspace return alone does not remove them.

### 8.2 Charged construction, actual inverse and fair baseline

Use elementary H, S and CNOT, with $`S^\dagger=S^3`$.
Controlled Y is $`S_t\mathrm{CNOT}_{u,t}S_t^\dagger`$, costing
five elementary Cliffords; the address phase is $`Z_u=S_u^2`$.
The displayed complete word has four T gates, fifteen elementary
Cliffords and four nonempty T layers, with zero clean qubits, zero extra
helpers, one arbitrary work/output qubit and two explicit address qubits.
Its Clifford charge is address Z:2; two controlled $`X^v`$:2;
controlled $`X^u`$:1; controlled $`Y^u`$:5; two H:2; fixed
$`S^\dagger`$:3. Every source and control is included.

The actual adjoint, retaining the reverse order, is

```math
V_{u,v}^\dagger
=X^vT^\dagger HT^\dagger X^uSY^uT^\dagger HTX^v(-1)^u.
```

This displayed inverse has four T gates, thirteen elementary Cliffords
and four nonempty T layers, on the same allocation: fixed
$`S^\dagger`$ becomes S, while controlled Y and address Z are
self-adjoint with their same charged implementations. The complete word
and its inverse both simplify to identity, whose fair baseline has zero
gates, zero depth and no additional workspace on the same logical
registers. No resource improvement is inferred from the unsimplified
witness.

The four
[returned-word checks](../../tests/test_returned_six_layer_odd_group.py)
verify the full eight-dimensional operator and actual inverse, the
elementary ledgers, the necessity of the literal address phase, the actual
normal-form comparison, both order-three calculations and the persistent
non-Pauli differences. The identities above establish the claim without
restricting any input or discarding any logical column.
