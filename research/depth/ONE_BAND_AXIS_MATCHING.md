# Matching native one-T bands through Cliffords

[Five-layer reduction](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Current depth frontier](../../docs/OPEN_PROBLEM.md)

A physical T layer has more structure than an arbitrary third-level
Clifford-hierarchy unitary. Its independent active Pauli axes can be
recovered from the Pauli coefficient norms of its conjugation action.
This gives an exact, and fixed-error, matching theorem for a single
native band. It does not assert that either intermediate passage in a
two-band cancellation is Clifford.

## 1. Statement with signs, inactive wires and error

Work on q arbitrary qubits. Let D and E be physical T/T-dagger layers,
possibly expressed in fixed Clifford coordinates. Their active Pauli
axes are respectively the independent commuting signed Hermitian
Paulis A_i and B_j. Inactive wires contribute no axis. Up to their
specified native scalar phases, these layers have the form

```math
D\doteq\prod_i e^{-i\epsilon_i\pi A_i/8},\qquad
E\doteq\prod_j e^{-i\delta_j\pi B_j/8},\qquad
\epsilon_i,\delta_j\in\{1,-1\}.
```

Axis matching below concerns binary Pauli labels. It allows the signs
of A_i and B_j, as well as the T orientations, to differ. All native
scalars remain in the actual operators D and E; projective notation
is used only for the necessary structural conclusion.

**One-band matching theorem.** Let C and C' be Clifford unitaries. If
there is a scalar omega of modulus one such that

```math
\|DC-\omega C'E\|_{\rm op}\lt\frac{\sqrt2}{4},
```

then C maps the set of active binary axes of E bijectively to the set
of active binary axes of D. In particular their active counts agree.
The threshold and conclusion are independent of q.

Conversely, if C has this axis-matching property, then DC E-dagger is
exactly a scalar times a Clifford, so some C' and omega satisfy the
displayed equality with zero error. The converse allows a new C'; it
does not say that an already specified nearby C' is exactly equal to
the native product.

The special case E=D gives an exact characterization:

```math
DCD^\dagger\text{ is projectively Clifford}
\quad\Longleftrightarrow\quad
C\text{ permutes the individual active axes of }D
\text{ up to signs}.
```

Preserving only their linear span is insufficient. For example, a
CNOT on two active T wires preserves their Z span but sends one
individual Z axis to their product; its D-conjugate is non-Clifford.

## 2. Recovering the axis list from Pauli coefficient norms

Use Hermitian Paulis tau(p), indexed by binary labels p, and the
symplectic pairing. Write a_i for the label of A_i and define

```math
k_D(p)=\#\{i:\langle a_i,p\rangle=1\}.
```

Conjugating tau(p) through an axis that commutes with it leaves it
unchanged. An anticommuting active axis replaces it by a sum of two
Paulis with coefficient magnitudes 1/sqrt(2). Since all active axes
commute, their commutation status with the current Pauli does not
change during the expansion. Since their labels are independent,
different subsets of those axes produce distinct resulting labels.
There are no collisions or cancellations. Therefore the Pauli
expansion of $`D\tau(p)D^\dagger`$ contains exactly $`2^{k_D(p)}`$
terms of equal magnitude $`2^{-k_D(p)/2}`$. Its coefficient l1 norm is

```math
\left\|D\tau(p)D^\dagger\right\|_{P,1}
=2^{k_D(p)/2}.                                             \qquad(1)
```

This norm is invariant under Clifford conjugation and under left or
right multiplication by a Pauli. These operations only permute Pauli
coefficients and multiply them by unit scalars. The expansion itself
is an operator identity, not a state-preparation argument.

The complete function k_D recovers the active binary axis list:

```math
k_D(p)=\frac r2-\frac12\sum_{i=1}^r
                 (-1)^{\langle a_i,p\rangle}.             \qquad(2)
```

The characters on the binary Pauli label space are linearly
independent. The constant coefficient gives r; each nonzero character
coefficient gives the multiplicity of its axis. Thus equality of two
such functions forces equality of their axis multisets. Independence
means these multisets are sets in the present native setting.

## 3. Exact and robust proofs

First suppose DC=omega C'E exactly. Conjugating a Pauli through both
sides and using (1) gives

```math
k_D(S_Cp)=k_E(p)\quad\text{for every }p.                  \qquad(3)
```

Equation (2) identifies the labels S_C inverse times a_i with the
labels b_j, proving the required matching.

For the error statement, put eta=norm(DC-omega C'E). For every p the
two conjugation images differ in operator norm by at most 2 eta.
Both images are Cliffords: a one-physical-T layer sends any Pauli to
a Clifford, and the other factors are Clifford conjugations.
Projective Cliffords with different symplectic actions are at least
1/sqrt(2) apart in operator norm, even after optimizing a scalar.
Indeed a Pauli witnessing the different symplectic actions has two
images at projective distance at least sqrt(2), while conjugating by
two nearby unitaries changes that image by at most twice their
operator distance.

Since 2 eta is strictly less than 1/sqrt(2), the two images have the
same symplectic action. They therefore differ by a scalar times a
Pauli. Their Pauli coefficient l1 norms are equal, so (3) follows
exactly, and the Fourier argument applies without a width-dependent
continuity estimate for the l1 norm.

For the converse, matching the axes means C E C-dagger and D rotate
the same individual axes. Their orientations differ, if at all, by
an even multiple of a T rotation. Consequently

```math
DC E^\dagger C^\dagger
```

is a scalar times a product of Pauli rotations through integer
multiples of pi/4, hence a Clifford. Multiplication by C proves the
converse. Native scalar phases may be assigned to omega; no literal
phase is discarded from a claimed circuit equality.

## 4. What this does and does not supply at five layers

In fixed outer Clifford coordinates, the two-T dressed Pauli B_y in
the [five-layer reduction](READ_ONLY_FIVE_LAYER_REDUCTION.md) is

```math
B_y=D A_y D^\dagger,
```

where D is one physical T layer and A_y is Clifford: the inner
physical T layer has already conjugated a Pauli to a Clifford.
For a Clifford K, its partial-normalizer passage is

```math
B_yKB_y^\dagger
=D\left[A_y(D^\dagger K D)A_y^\dagger\right]D^\dagger.
                                                               \qquad(4)
```

The theorem characterizes either one-band passage **when its input
and output are known to be Clifford**. If K permutes D's axes, the
first passage D-dagger K D is Clifford. If (4) is also Clifford,
the theorem then characterizes the second passage. This can be a
useful sufficient route through a specified instance.

Complete return does not force the first passage to be Clifford.
Already on one arbitrary qubit, take D=T and

```math
A=T^\dagger X T=\frac{X-Y}{\sqrt2},\qquad
B=DAD^\dagger=X,\qquad K=H.
```

A is a one-T Pauli image and hence Clifford. The full output B K
B-dagger=X H X is Clifford, but D-dagger H D is not: H sends the
active Z axis to X. Thus even exact cancellation can pass through a
non-Clifford intermediate. This native witness uses the actual
inverse layers. Its displayed B word has four T/T-dagger gates in
four sequential slots and four elementary Cliffords, using X=HSSH
in the fixed H,S,CNOT alphabet. It uses no clean qubit or additional
dirty helper; the one work qubit is arbitrary. Its fair simplified
baseline is X, with zero T count and depth and four elementary
Cliffords. The full comparison BHB-dagger has eight T calls in eight
unsimplified slots and nine elementary Cliffords; its simplified
baseline XHX has zero T count and depth and nine elementary
Cliffords before further Clifford simplification. All these
displayed equalities are exact full-operator equalities, with zero
leakage. They supply a failure of an inference, not a resource
improvement.

The [general five-layer theorem](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
subsequently resolves the output-degree question by applying this matching
lemma to a native Pauli defect. It does not prove a general degree law for
arbitrary two-band action on Clifford tableaus. The current
[six-layer question](READ_ONLY_SIX_LAYER_REDUCTION.md#6-a-conditional-defect-detector-for-the-general-comparison)
has a Clifford defect and three-band prefix, with intermediate Clifford
membership still unavailable.

No auxiliary qubit or physical operation is introduced by this proof.
The inequalities concern full operators on arbitrary inputs,
including external references. T-count, elementary Clifford count
and simultaneous finite workspace remain unrestricted in this
necessary invariant; no new upper bound or work-return compiler is
claimed. The [bounded checks](../../tests/test_one_band_axis_matching.py)
test the expansion, matching, inactive axes, literal signs and failed
intermediate inference. The analytic proof supplies the
width-independent theorem.

## 5. An abstract normalizer can have arbitrarily high tableau degree

The native two-T restriction cannot be discarded in favor of an arbitrary
unitary normalizer, even if that unitary is an involution and an abstract
Pauli image. The following construction makes this distinction on the
**full binary Clifford tableau**, including the Pauli signs.

Let z be a register of a arbitrary work qubits. Let the data bank contain
one arbitrary qubit for every u in the binary vector space with a
coordinates. For each x, let C_x permute these data wires by
$`u\mapsto u+x`$. These are Clifford operators with exact laws
$`C_xC_y=C_{x+y}`$ and $`C_x^2=I`$. Define

```math
U=\sum_z |z\rangle\!\langle z|\otimes C_z,
\qquad U=U^\dagger,\quad U^2=I.
```

Direct action on every basis vector, with all amplitudes retained, gives

```math
U(X_z^x\otimes I)U^\dagger=X_z^x\otimes C_x,
\qquad U(I\otimes C_x)U^\dagger=I\otimes C_x.
```

Thus U normalizes the abelian Clifford subgroup generated by the two
families $`X_z^x\otimes I`$ and $`I\otimes C_x`$. The input family
$`K_x=X_z^x\otimes I`$ has identity symplectic matrix; its only varying
tableau entries are Pauli signs linear in x. The conjugated family's
symplectic matrix has a data-wire permutation block with entries

```math
\delta_{v,u+x}
=\prod_{j=1}^{a}(1+u_j+v_j+x_j)
\quad\text{over }\mathbb F_2.
```

Every such entry has Boolean degree a. Consequently no absolute bound on
the degree of the full-tableau partial action follows from subgroup
normalization alone: a polynomial of bounded degree in the input tableau
entries would restrict to that same bounded degree in x.

To impose balanced multiplicities as well, add one arbitrary work qubit
and put $`B=U\otimes X`$. The extra X makes the positive and negative
eigenspaces of the Hermitian involution B equal in dimension. Choosing
orthonormal eigenbases therefore writes B as an arbitrary-unitary
conjugate of a physical Pauli. Its action on the displayed subgroup is
unchanged. This is an abstract Pauli-image statement, not a two-T native
construction.

The source does have a fully charged, growing-depth native realization.
For each coordinate j, conditionally swap the data pairs
$`u,u+e_j`$ with $`u_j=0`$, controlled by $`z_j`$. The resulting
$`a2^{a-1}`$ Fredkin gates implement U on every arbitrary work input.
Each Fredkin uses two CNOTs and one Toffoli. Using the literal Toffoli
expansion with seven T gates, ten CNOTs and two H gates gives the
conservative serial ledger

```math
T(U)=7a2^{a-1},\qquad
G(U)=14a2^{a-1},\qquad
D_T(U)\le5a2^{a-1}.
```

For B add four elementary Cliffords for $`X=HSSH`$. The simultaneous
work allocation is $`a+2^a+1`$ arbitrary qubits, with zero clean qubits
and no additional helper. Its actual inverse is the reversed inverse
native circuit; Hermiticity does not make either source call free. The
identities are exact, with zero leakage and no approximation error, and
hold with arbitrary external references. U and B themselves intentionally
permute work. The displayed conjugation implements
$`X_z^x\otimes C_x\otimes I`$ exactly; it does not return the data bank
to identity and is not a qualifying returned XOR query. Both calls of B
and its actual inverse are charged when that conjugation is used: its
T count and depth bounds are twice the source bounds, and its Clifford
count is $`2G(B)+G(P_x)`$ for the specified Pauli word $`P_x`$.
For a fixed literal X word, $`G(P_x)=4|x|`$ using HSSH; for a
coherent a-bit address controlling that word, it is a CNOTs. The latter
conjugation still permutes the bank and is not a returned XOR query.

This ledger supplies no constant-depth construction. Conversely, the
example does not prove that another constant-depth native representation
is impossible. It falsifies only the proposed shortcut from an abstract
partial-normalizer statement to a bounded-degree tableau action. The
native two-T conjugator and the special five-layer address cocycle remain
additional structure that a successful general argument must use.

The [scope checks](../../tests/test_partial_normalizer_scope.py) verify the
full-tableau degree, exact permutation identities and native Fredkin word
in bounded instances. The arbitrary-width conclusion follows from the
displayed product formula, not from enumeration.

## 6. A projected band with one relation

Restricting an independent physical band to a one-stabilizer code can
create at most one binary relation among its projected axes. The
following extension recovers the individual axes except for one explicit
three-dimensional ambiguity. It is an invariant of the projected band,
not a general one-clean three-layer query theorem.

Let E contain commuting nonzero projected axis labels
$`a_1,\ldots,a_m`$ with at most one relation. Keep all axis signs,
T orientations and native scalar phases. A projected identity contributes
only its actual scalar. A relation of length two combines its two equal
axes into an exact commuting Clifford or scalar; absorb that factor
into the adjacent Clifford before applying the result below.

For a Pauli label p define

```math
k(p)=\#\{i:\langle a_i,p\rangle=1\},\qquad
w_E(p)=2\log_2\|E\tau(p)E^\dagger\|_{P,1}.
```

An independent active subset has the usual flat expansion, so w=k.
Let R be the unique relation's support, of length r. The only dependent
active subset contains all of R. This can happen only when r is even:
pairing the relation with p makes its number of active members even.
In the expansion, subsets S and their symmetric difference with R
give the same Pauli label. Their coefficient ratio is a fixed sign
times $`(-1)^{|S\cap R|}`$. Since r is even, this parity is well
defined on each paired class; exactly half the classes cancel and
half double. Thus there are $`2^{k-2}`$ surviving terms of magnitude
$`2^{-(k-2)/2}`$. Consequently

```math
\begin{gathered}
w_E(p)=k(p)-2I_R(p),\\
I_R(p)=\mathbf1[\langle a_i,p\rangle=1\text{ for all }i\in R]
                                                               \qquad(6.1)
\end{gathered}
```

for an even relation. For an odd relation w=k everywhere. This proof
does not assume that a dependent expansion has the independent support
or normalization.

Put $`\chi_v(p)=(-1)^{\langle v,p\rangle}`$. For even r,

```math
I_R=2^{-r}\prod_{i\in R}(1-\chi_{a_i}).
```

Every label in $`L_R=\mathrm{span}\{a_i:i\in R\}`$ has
exactly two subset representations, differing by R. Their common
cardinality parity is denoted pi(v). The Fourier coefficient of I_R
at v is $`2^{1-r}(-1)^{\pi(v)}`$. Since
$`k=m/2-\tfrac12\sum_i\chi_{a_i}`$, setting
$`\delta_r=2^{2-r}`$ gives the nonconstant coefficients of w:

| Character label | Coefficient |
| --- | --- |
| An axis outside $`L_R`$ | $`-1/2`$ |
| A relation axis | $`-1/2+\delta_r`$ |
| Another nonzero label in $`L_R`$ | $`-\delta_r(-1)^{\pi(v)}`$ |
| Any other label | Zero |

No outside axis lies in $`L_R`$, since that would create a second
relation. For $`r\ge6`$, coefficients of magnitude strictly greater
than 1/4 therefore recover exactly the individual axes. Independent
lists and odd relations retain the original coefficient -1/2 recovery.

The sole non-Clifford exception is r=4. Its relation span has dimension
three. Its odd coset consists exactly of the four relation axes; all
seven nonzero labels in the span have coefficient -1/4. The profile
recovers this span, but not the original four-axis subset. All outside
axes still have coefficient -1/2 and are recovered individually.

**Matching consequence.** Suppose D,E are two such projected bands and
C,C' are Clifford on the quotient. If

```math
DC=\gamma C'E,
```

then $`w_D(F_Cp)=w_E(p)`$ for every p. Therefore C matches the
individual axes, except that a length-four relation supplies only a
matched three-dimensional span. The profile distinguishes that exception
from the independent, odd-relation and longer-even-relation cases.
Span matching in the exceptional case is necessary, not claimed
sufficient for the displayed Clifford equality.

The same necessary conclusion holds if the displayed products are at
operator distance strictly below $`\sqrt2/4`$, optimizing gamma. Each
projected phase-polynomial band still sends Paulis to Cliffords. The
symplectic gap from Section 3 makes the compared images have equal
symplectic actions; they differ by a scalar Pauli and hence have equal
coefficient l1 norms exactly. No width-dependent l1 continuity estimate
is used.

The [exact fixtures](../../tests/test_projected_band_matching.py)
check dependent Pauli expansions and the exceptional Fourier spectra.
The component controls a band that really restricts to a Pauli
stabilizer code. It supplies no such restriction for a non-Pauli
intermediate code and does not turn a general dependent phase polynomial
into a native band. That remaining distinction is essential in the
[three-layer continuation](THREE_LAYER_NATIVE_STRUCTURE.md).
