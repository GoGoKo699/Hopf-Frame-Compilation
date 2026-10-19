# Computational derivatives of the completed clean-transition phase

[Common determinant](NONPAULI_CODE_DECODING.md#5-the-zero-pair-zero-character-logical-defect-is-abelian) ·
[Rank-one interface](RANK_ONE_CLEAN_PHASE_INTERFACE.md) ·
[Two-layer diagonal rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md)

The intrinsic determinant phases can be connected to computational
derivatives without bounding their number of central variables. An
actual clean-odd Pauli gives a Clifford relation between the occupied
and unused blocks. Under the zero-character determinant hypotheses,
each associated target derivative has an actual two-native-band
logical word and therefore Boolean degree at most three. The resulting
quartic-plus-quotient factorization is a statement about the target;
removing its quartic factor at the same native depth is not proved.

For an actual no-pair clean transition with both interior clean images
non-Pauli, Boolean reality also identifies a shared computational
center. Section 10 reduces the remaining active transition to a central
cubic sign and one Pauli controlled by a quadratic function of those
actual input parities, retaining three independent logical native
lists. Section 9 closes the actual-source subcase with a noncentral
common endpoint parity, giving sign degree at most five and XOR degree
at most four under the stated sign-representation hypothesis. The
remaining question is global matching with the quadratic
control; fixing its central variables need not preserve native
independence.

## 1. Actual split and occupied Boolean action

Let Q count the clean wire c and every other physical wire. Assume
complete clean return and write the actual native three-band circuit as

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,
\qquad U_C=E_CK_C,\quad U_{AB}=E_BE_AK_A.
```

Each E has independent commuting signed physical axes, with its literal
T or T-dagger orientations. Let

```math
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger,
\qquad J=E_AP_AE_A^\dagger,\quad S=E_CP_CE_C^\dagger.
```

Suppose S is non-Pauli and the native clean comparison has no paired
middle orbit. A middle singleton commutes every C-active axis: its
binary displacement under S vanishes, and independence of the native
C list forces every active-axis pairing to vanish. This conclusion
does not require a bound on either endpoint active rank.

On the positive clean sector let

```math
G|x,w\rangle=e^{i\theta}(-1)^{g(x)}|x,w\rangle,
```

where x has n designated sign variables and w contains arbitrary dirty
variables. The common scalar is retained. Unitarity and exact return
give a full-space block decomposition $`V=\mathrm{diag}(G,U)`$.
No claim about the unused block U is part of the target specification.

## 2. An actual Clifford relation between the two clean blocks

Choose an actual signed C-active axis Q_0. It commutes every E_C axis
and, by the no-pair premise, every E_B axis. It anticommutes S. Hence

```math
q_C=U_C^\dagger Q_0U_C=K_C^\dagger Q_0K_C,
\qquad
q_A=U_{AB}^\dagger Q_0U_{AB}
    =K_A^\dagger E_A^\dagger Q_0E_AK_A.             \qquad(1)
```

The first is a clean-odd Hermitian Pauli. The second is a clean-odd
Hermitian Clifford, since a one-band inverse is in C3. Both are
actual full-space operators and satisfy $`V^\dagger q_CV=q_A`$.

Choose an actual occupied-fixing Clifford D with
$`DX_cD^\dagger=q_C`$. A clean-odd Pauli has off-diagonal blocks
equal to a logical Pauli times a fourth root, so D can be chosen as
a Pauli controlled by c, with a clean S or Z correction if required.
It costs $`O(Q)`$ elementary Cliffords, uses no T gates or helper,
and fixes every clean-zero vector literally. Set

```math
V'=D^\dagger V,\qquad K_C'=K_CD,
\qquad U_C'=E_CK_C'.
```

The actual native depth and occupied G are unchanged. Since
$`(V')^\dagger X_cV'=q_A`$, its block identity is

```math
q_A=
\begin{pmatrix}0&L\\L^\dagger&0\end{pmatrix},
\qquad
V'=\begin{pmatrix}G&0\\0&GL\end{pmatrix}.
                                                               \qquad(2)
```

Here L is an actual logical Clifford. Indeed $`X_cq_A`$ is a
clean-preserving actual Clifford, whose two block restrictions are
$`L^\dagger,L`$. Restriction of that full-space Clifford, rather
than extension of an occupied equality, supplies the assertion.

The exact operation $`D_L=\mathrm{diag}(I,L)`$ belongs to C3.
It fixes Z_c; it sends a logical Pauli P to
$`\mathrm{diag}(P,LPL^\dagger)`$, which is a Clifford because
the two blocks are Paulis. Its image of X_c is $`X_cq_AX_c`$,
also a Clifford. Thus $`V'=(I_c\otimes G)D_L`$ is an actual
block factorization. This gives no free implementation of D_L, no
bound on its support, and no right to erase it from the three-band word.
Section 2 needs no zero-character assumption.

## 3. Central defects give actual two-band computational derivatives

Now add the determinant hypotheses: S has intrinsic active rank other
than two, and every complete logical computational-Z probe character
is zero. The
[common determinant lemma](NONPAULI_CODE_DECODING.md#51-a-common-determinant-coboundary-constrains-the-phases)
then gives

```math
K=L_B+L_S+\langle p_C\rangle,
\qquad W_z=C_zA_z=\Pi_++\alpha_z\Pi_-,
\qquad \Pi_\pm=(I\pm S)/2.                         \qquad(3)
```

Every alpha_z is a full Clifford in the center of Alg(K), including
its literal scalar. These are intrinsic central Clifford phases;
their coordinates have not been identified with x or w.

The normalization D only changes
$`C_z`$ to $`C_zS^{\chi(z)}`$, because D conjugates each
logical Z(z) to $`Z(z)Z_c^{\chi(z)}`$. Thus it multiplies
alpha_z by the actual sign $`(-1)^{\chi(z)}`$. The hypotheses
and central Clifford conclusion remain valid. Use the normalized
symbols C_z, W_z and alpha_z below.

Define the actual Pauli unitary

```math
P_z=Z(z)LZ(z)L^\dagger.                            \qquad(4)
```

Its scalar may be a fourth root; do not replace it by an independently
chosen Hermitian representative. Pulling (3) through U_C' and using
(2) gives

```math
(U_C')^\dagger W_zU_C'
 =\mathrm{diag}(I,GP_zG^\dagger).                  \qquad(5)
```

There is an exact logical native word for the nontrivial block in (5).
Centrality makes alpha_z commute P_C, Q_0 and every C-active axis.
Also $`K_C'Z_c(K_C')^\dagger=P_C`$ and
$`K_C'X_c(K_C')^\dagger=Q_0`$. The entire C-silent band commutes
both P_C and Q_0. Consequently there are literal logical operators

```math
\begin{aligned}
(K_C')^\dagger E_{C,\mathrm{silent}}K_C'&=I_c\otimes E_s,\\
(K_C')^\dagger\alpha_zK_C'&=I_c\otimes A_z^{\mathrm{cen}},\\
(U_C')^\dagger\alpha_zU_C'&=
           I_c\otimes(E_s^\dagger A_z^{\mathrm{cen}}E_s).
\end{aligned}                                                     \qquad(6)
```

The logical E_s is one actual independent-axis native band. Projection
is injective here: a Pauli commuting with both clean X and clean Z
has identity on the clean wire. Thus no relation is introduced from
the physical independent C-silent list. Each $`A_z^{\mathrm{cen}}`$
is an actual logical Clifford. Equation (6) uses the original endpoint
Clifford and band; it is not a conditional choice of coordinate circuit.
Comparing (5)--(6) proves

```math
GP_zG^\dagger=E_s^\dagger A_z^{\mathrm{cen}}E_s.    \qquad(7)
```

Let u_z be the X mask of P_z, including all non-clean wires. Removing
that literal Pauli yields

```math
\mathrm{diag}\bigl((-1)^{g(x)+g(x+u_{z,x})}\bigr)
 =E_s^\dagger A_z^{\mathrm{cen}}E_sP_z^\dagger.     \qquad(8)
```

Here u_{z,x} is its projection onto the designated sign variables;
the dirty translation disappears because g is independent of w.
The common scalar of G cancels. The right side is a full logical
two-native-band circuit with only actual Clifford interlayers and
endpoints. All logical and dirty inputs are arbitrary. The
[zero-clean two-layer diagonal theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#1-model-and-full-space-diagonal-theorem)
therefore makes (8) diagonal C3, proving

```math
\deg_{\mathbb F_2}\Delta_{u_{z,x}}g\le3.            \qquad(9)
```

No assertion about the degree of intrinsic-coordinate alpha_z was
transferred directly to computational variables. Its actual Clifford
operator, and the retained native endpoint band in (8), supply that
transfer.

## 4. A quartic part and a residual quotient

Let F_L be the binary symplectic action of the actual Clifford L.
Since the initial Z(z) in (4) has zero X mask, the map

```math
z\longmapsto u_{z,x}
 =\mathrm{pr}_x\,X(F_LZ(z))
```

is linear. Let U_x be its image in the n-dimensional sign-variable
space and put $`d=\dim U_x`$. Derivatives in every direction of
U_x have degree at most three, either directly from (9), or from
the identity $`\Delta_{u+v}g=\Delta_ug+\Delta_vg(\,\cdot+u)`$.

Choose linear sign coordinates whose first d directions span U_x.
In the unique ANF, any monomial involving one of those coordinates
has degree at most four: differentiating in that coordinate removes
exactly that variable and cannot cancel with another monomial.
Hence there is an exact target factorization

```math
g(x)=q_4(x)+h(\pi x),\qquad
\deg q_4\le4,\qquad
\pi:\mathbb F_2^n\longrightarrow\mathbb F_2^n/U_x.
                                                               \qquad(10)
```

In particular,

```math
\deg g\le\max\{4,n-d\}.                           \qquad(11)
```

This is uniform in endpoint active ranks, central dimension and dirty
width whenever the computational residual codimension $`n-d`$ is
fixed. It does not require bounded determinant-defect codimension r.
The earlier identity family has arbitrarily large r and, if all its
non-clean variables are designated, full U_x.

For an XOR query, form the sign representation by Hadamards on its
output before applying these hypotheses. Its sign is
$`(-1)^{y f(x)}`$; for nonzero f the corresponding bound on its
Boolean degree is one less than (11). No invariance of U_x under that
output basis change is presumed. The retained R74 rounding consequence
first forces exact occupied Boolean action below $`2^{-530690}`$;
the exact circuit can then be analyzed by the theorem above.

The factorization (10) is not a same-depth circuit factorization.
A quartic sign need not be a Clifford, and no native word removing
q_4 within the original three bands has been supplied. Applying a
three-band theorem to h after deleting q_4 would be unsupported.

## 5. More clean flips do not enlarge this derivative space

The mask space in Section 4 is invariant under all admissible choices
of the clean flip, including choices arising from another split of
the same full operator V. Indeed let two occupied-fixing controlled
Pauli normalizations give

```math
V_q=\mathrm{diag}(G,GL_q),\qquad
V_v=\mathrm{diag}(G,GL_v).
```

Their relation has the form

```math
L_qL_v^\dagger=G^\dagger R_{qv}G,                 \qquad(12)
```

where R_{qv} is a literal logical Pauli unitary. The left side is
Clifford. Since G is diagonal, the right side has the same constant
XOR permutation of computational basis vectors as that Pauli. It is
therefore a Pauli times a diagonal Clifford. Such a Clifford preserves
the X mask of every Pauli label, giving

```math
X(F_{L_q}Z(z))=X(F_{L_v}Z(z))\quad\text{for every }z.              \qquad(13)
```

Taking a span over more clean flips cannot enlarge U_x. The assertion
concerns one fixed full V; changing its unused-sector extension while
preserving the occupied target is a different operation.

There is also an actual arbitrary-width zero-space witness. If
$`V=I_c\otimes G`$, then every admissible L_q has the form
$`G^\dagger R_qG`$, so the same argument gives U_x=0. The
[disjoint-CCZ family](CLEAN_TRANSITION_PARAMETERS.md#1-three-independent-clean-active-lists)
has this full operator, both clean interiors non-Pauli, no middle
pairs, and C active rank $`2n\ne2`$ for its positive even block
count n. Because full V is diagonal, every logical comparison has
$`A_z=C_z`$: all characters and defects vanish. Thus all determinant
hypotheses hold, but no choice of active clean flip or resplitting of
that same operator yields nonzero derivative coverage. Its target
degree is only three, so this is a limitation of a uniform-coverage
inference, not a high-degree native-query counterexample.

## 6. A broader derivative theorem for a commuting-C3 sandwich

There is a second computational derivative argument with different
hypotheses. Use the actual
[four-factor reduction](RANK_ONE_CLEAN_PHASE_INTERFACE.md) when both
interior clean images are non-Pauli and there are no middle pairs.
It does not require zero probe characters or exclude active rank two.
On the full logical register, including every arbitrary dirty wire,
write its literal expression as

```math
G=K_{\mathrm{out}}^\dagger C_s^\dagger FBAK_{\mathrm{in}}
 =K_{\mathrm{out}}^\dagger C_s^\dagger DAK_{\mathrm{in}},
\qquad D=FB.                                      \qquad(14)
```

The endpoint lists A and C_s are independent logical native bands.
The axes of F and B commute with one another, so D is an actual
commuting Pauli C3 operation. Their projected relation dimensions
remain $`(0,\le1,\le1,0)`$ in input-to-output order. D is used
as an exact C3 operation, not credited as a free physical band.
In fact the following argument applies to every full-space Boolean
sandwich of the second form in (14), with arbitrary commuting C3 D.

Let L_C be the axis span of C_s. For any signed Hermitian Pauli R
with label in $`L_C^\perp`$, put
$`Q_R=K_{\mathrm{out}}^\dagger R K_{\mathrm{out}}`$.
It has a definite X mask a_R. R commutes C_s literally, so

```math
G^\dagger Q_RG
 =K_{\mathrm{in}}^\dagger A^\dagger J_RA K_{\mathrm{in}},
\qquad J_R=D^\dagger RD.                          \qquad(15)
```

J_R is an actual Clifford. Removing the literal Hermitian Pauli Q_R
gives

```math
Q_R G^\dagger Q_RG
 =[Q_RK_{\mathrm{in}}^\dagger]A^\dagger J_RA K_{\mathrm{in}}
 =\mathrm{diag}\bigl((-1)^{\Delta_{a_{R,x}}g(x)}\bigr).
                                                               \qquad(16)
```

The common scalar of G cancels, and $`Q_R^2=I`$ fixes the Pauli
phase convention. The right side is an actual full-space two-native-band
word, so its Boolean derivative has degree at most three. Time reversal
gives the same conclusion for the X masks of
$`K_{\mathrm{in}}^\dagger L_A^\perp K_{\mathrm{in}}`$,
where L_A is the span of the actual input band. This also includes
masks of the endpoint silent axes themselves, because each native
span is isotropic and lies in its own symplectic orthogonal.

Let U_sand be the span of both mask images, projected onto the
designated sign variables x. The ANF argument of Section 4 gives

```math
g=q_4+h_{\mathrm{sand}},\qquad \deg q_4\le4,
```

where h_sand depends only on the quotient by U_sand. Equivalently
its allowed computational parity labels form

```math
\begin{aligned}
T_x=U_{\mathrm{sand}}^\perp
 =\{z_x:\;& Z(z_x)\otimes I_w\text{ has label in }\\
 &K_{\mathrm{out}}^\dagger L_CK_{\mathrm{out}}
 \cap K_{\mathrm{in}}^\dagger L_AK_{\mathrm{in}}\}.
\end{aligned}                                                     \qquad(17)
```

On the whole logical register, this is simply the intersection of
those two transformed endpoint spans with the computational Z space M.
The identity on w makes the active-variable projection explicit.
Thus (17) locates the higher-degree part in actual computational
parity coordinates. The dimension of this common intersection need
not be bounded by the one-core interface.

When both sets of hypotheses hold, U_sand and the determinant mask
space U_x may be combined by taking their linear span. Every derivative
in that span still has degree at most three, so the quartic-plus-quotient
statement holds for the smaller quotient. This does not assert that
either space, or their sum, has uniformly bounded codimension, and it
does not supply free native synthesis of q_4.

## 7. Boolean reality gives a shared computational center

This lemma concerns the actual full-space four-factor form

```math
V=K_C^\dagger C_s^\dagger FBAK_A=e^{i\theta}D_g,
```

where g is Boolean, A is an independent commuting native T/T-dagger
band, and every factor retains its literal scalar. Suppose an isotropic
Pauli-label space R is centralized pointwise by B, F and C_s. In the
rank-one interface take R to be the span of the differences between the
F axes. No initialized restriction is used below.

For a signed Hermitian Pauli R_r in this space, put
$`Q_C(r)=K_C^\dagger R_r K_C`$. Then

```math
V^\dagger Q_C(r)V=K_A^\dagger A^\dagger R_r A K_A.
```

The right-hand side is a Hermitian Clifford. The left-hand side is a
Pauli times a diagonal Boolean derivative. Consequently its binary
Clifford displacement has image in the computational Z space M.
The independent-native intrinsic-axis lemma now shows that every A
axis active for R_r lies in $`K_A M K_A^\dagger`$. Boolean reality gives
more. Let A_i range over those active axes and let
$`Z_i=K_A^\dagger A_i K_A`$ be their signed computational Z representatives.
The exact one-band conjugation formula writes the image as
$`K_A^\dagger R_r K_A`$ times a common scalar and a product of $`S_{Z_i}^{+1}`$
or $`S_{Z_i}^{-1}`$. The two endpoint Pauli cores have the same X mask,
since their quotient is diagonal. Therefore the Boolean derivative
itself is a scalar times a computational Z Pauli times this product
of odd S powers. Squaring it gives the identity. Its binary Pauli label
is the sum of the active Z_i labels, so that sum must vanish. Physical
independence of the original A list forces the active set to be empty.
This holds for every r in R. Hence **A already commutes with the entire
center R**, and

```math
V^\dagger Q_C(r)V=Q_A(r):=K_A^\dagger R_rK_A
```

is an exact Pauli identity. In particular every Boolean derivative of
g along the X-mask space U of Q_C(R) is affine, of degree at most one.
This sharpens the mere Clifford-image derivative bound and requires no
C3 correction or removed native factor.

There is an actual real diagonal Clifford C matching the entire family
simultaneously. To see this directly in computational variables, choose
an invertible linear coordinate change putting U in the first u
coordinate directions. The ANF derivative condition says that every
monomial of g involving one of those directions has degree at most two.
Collect precisely those monomials into a quadratic Boolean q, transform
back to the original coordinates, and set C=D_q. This is a concrete
Z/CZ Clifford, including its chosen literal scalar, with

```math
C^\dagger Q_C(r)C=Q_A(r)\quad\hbox{for every }r\in R.
```

This extraction is legitimate because a quadratic Boolean phase is an
actual Clifford; no higher-level correction is being charged for free.
The literal residual

```math
V_1:=VC^\dagger
 =K_C^\dagger C_s^\dagger FBA(K_AC^\dagger)
```

is still a Boolean diagonal sign up to the same common scalar and
commutes pointwise with every Q_C(r). All four native/algebraic middle
factors are unchanged.

The common center can now be made computational without losing
diagonality of the residual target. Let U be the X-mask space of Q_C(R),
with dimension u. A CNOT change of computational variables makes U the
first u coordinate directions. An S/CZ/Z Clifford removes the Z
parts and signs of a chosen commuting family of lifts of those
directions. This is possible by the same symmetric-matrix extension;
it fixes every computational Z operator. The remaining zero-mask
kernel is put in coordinate Z form by CNOTs on the other variables.
The diagonal V_1 is unchanged by the diagonal Clifford and, because
it commutes with the resulting $`X_1,\ldots,X_u`$, acts as identity on those
u factors. Applying H to precisely those factors therefore leaves
V_1 diagonal and unchanged. This converts all of Q_C(R) to a
computational Z subgroup. At the level of the residual diagonal phase,
only the preceding invertible linear input change affects its arguments.

Use the corresponding common fixed Clifford coordinates on the middle
space as well. Every middle factor and both resulting endpoint
Cliffords preserve this computational Z subgroup pointwise. This is a
shared computational center on the entire register, not a collection
of separately chosen fiber extensions. The final input phase recovers
the original Boolean target by multiplying the residual phase by C,
a real quadratic diagonal Clifford. Their literal scalars and actual
inverse words are retained.

The operation extracts an actual quadratic Clifford and uses fixed
Clifford coordinate changes; all displayed phase lists are retained.
It adds no free physical T band. Fixing the shared
center's eigenvalues can nevertheless create an unbounded relation
dimension among the remaining projected native axes, so a native
matching theorem cannot simply be applied independently to each fiber.

### 7.1. The residual F factor has two actual-input control bits

In this common-center form suppose all F axes lie in one affine Pauli
coset of R. Choose one fixed signed Hermitian representative P, commuting
with the center. Each signed axis has the literal form
$`(-1)^{a_i(z)}P`$ in the center fiber z, where a_i is an affine Boolean
parity of genuine computational center coordinates. Let
$`\sigma_i=\pm1`$ be its T orientation, and define integer functions

```math
A(z)=\sum_i\sigma_i a_i(z),\qquad s=\sum_i\sigma_i.
```

The exact identity $`T_{-P}=\omega T_P^\dagger`$ gives

```math
F(z)=\omega^{A(z)}T_P^{s-2A(z)}.
```

Let epsilon be s modulo 2, chosen in {0,1}, and let b_0,b_1 be the two
binary digits modulo four of (s-epsilon)/2-A(z). Then

```math
F(z)=\omega^{A(z)}T_P^\epsilon S_P^{b_0(z)}P^{b_1(z)}.
```

Here b_0 is affine and b_1 has Boolean degree at most two. This follows
by expanding the integer sum of signed affine Boolean parities modulo
four: the first digit is affine, and the carry digit is quadratic.
The scalar-center factor $`\omega^{A(z)}`$ is a diagonal commuting C3 phase
and commutes with every remaining factor, including the input endpoint.
It therefore peels literally as a known diagonal input phase. What
remains is a fixed $`T_P^\epsilon`$ followed by one of four fixed Cliffords
on the single Pauli direction P, selected by two bounded-degree bits
of actual computational input coordinates. No central coordinate has
been silently identified with an input bit: that identification was
proved above by Boolean rigidity and actual Clifford alignment.

This finite control alphabet does not yet prove a uniform Boolean
degree bound. The unresolved step is a matching theorem retaining the
global independent native lists when the insertion has these actual
central controls. Restriction to one control value is not a free
full-space diagonal extension and can destroy projected independence.


## 8. Scope and verification

The actual quadratic phase in Section 7 can be compiled using
$`O(Q^2)`$ Z/CZ gates. Its fixed symplectic coordinate changes
also admit $`O(Q^2)`$ elementary Clifford compilations. Each
actual inverse is obtained by reversing and inverting its displayed
word and carries its own corresponding ledger. No helper, initialized
projection, or conditional coordinate circuit is used. The central C3
phase in Section 7.1 is extracted as a literal commuting input phase;
this identity is not a claim of free synthesis within the original
three native bands.

Section 9 gives a degree bound under its stated common-parity
hypothesis, sharpened to five/four for the actual clean source in
Section 9.4. Section 10 strengthens the actual clean-derived normal
form to a central Boolean cubic and a quadratic-controlled Pauli.
The cubic is retained as an exact input phase, and the controlled
Pauli keeps its literal source word and inverse; neither is assigned
a free synthesis or a free native-band deletion.

The [exact fixtures](../../tests/test_computational_phase_derivatives.py)
check all occupied columns and the actual inverse of a native
three-band CZ example, its three admissible clean flips, every logical
Z defect, retained fourth-root phases, and every derived computational
derivative. This example is not sharp for the quartic bound. A separate
four-axis commuting-list guard demonstrates why native independence
is essential to the Boolean reality argument: a dependent list can
have a non-Pauli real quadratic derivative.

A four-wire independent-axis fixture uses different orthogonal F and
active-B bases with identical Gram matrices. It verifies the actual
four-band word and inverse, the exact aligning Clifford, and the
retained central sign $`e_2(z_1,z_2,z_3)+z_1z_2z_3`$ on every column.

The remaining general step is to control the completed target phase
while retaining the original globally independent lists. The proved
quartic quotient is not a free circuit deletion, and the actual
shared-center normalization does not make the axes independent in
every restricted fiber. Neither limitation is a degree counterexample
or a change to the established resource frontier.

## 9. A noncentral common parity closes an independent active-middle subcase

There is a concrete bounded-degree consequence of the common endpoint
parity space in Section 6. Use the noncommuting rank-one interface,
where the joint output-axis result makes F independent, and first apply
the actual shared-center normalization of Section 7. Keep the original
F list rather than replacing it by its control formula. Write

```math
V=K_C^\dagger C_s^\dagger FBAK_A=e^{i\theta}D_g,
\qquad
R_F=\mathrm{span}\{f_i+f_0\}.
```

R_F is now a shared computational center. The A and C_s lists are
independent; B has at most one relation. The input quadratic Clifford
and affine coordinate changes used for this normalization are retained.
They change Boolean degree by at most two and preserve every global
list relation bound.

**Conditional degree theorem.** Suppose there is a computational parity
label t in the common endpoint space

```math
T=M\cap K_A^\dagger L_AK_A\cap K_C^\dagger L_CK_C
```

such that $`P=K_AZ(t)K_A^\dagger`$ anticommutes with at least one
F axis. Suppose further that the sublist B_a of B axes anticommuting
with P is independent. Then the original Boolean sign has degree at
most 2048. A literal XOR predicate whose Hadamard-conjugated circuit
meets these hypotheses has degree at most 2047.

These are additional hypotheses on a specific full-space representation.
They do not assert that the common parity space contains such a probe,
or that every active middle sublist is independent. The probe is
evaluated in the normalized coordinates. Their residual Boolean phase
differs from the original only by the retained invertible affine input
change and quadratic sign; the Hadamards act only on its exact identity
modes. This suffices to transfer the resulting degree bound back to the
original input variables.

### 9.1. The active middle span is the F span

Put $`P'=K_CZ(t)K_C^\dagger`$. The endpoint phases commute
with their respective Pauli cores. The exact diagonal target therefore
gives

```math
(FB)P(FB)^\dagger=P'.
```

P commutes with every central difference in R_F, so every F axis is
active for P. Let B_0 be the complementary middle sublist. Its axes
commute P, F and B_a. Thus the commuting C3 phase

```math
D_a=F B_a
```

already sends P to P' literally. Its Clifford displacement vanishes,
giving the binary Gram identity

```math
\sum_{f\in F}ff^{\mathsf T}
 =\sum_{b\in B_a}bb^{\mathsf T}.
```

For an independent list, this Gram matrix has image equal to the
axis span. Independence of both F and B_a consequently gives

```math
W:=\mathrm{span}(F)=\mathrm{span}(B_a),\qquad
W\cap P^\perp=R_F.                                \qquad(18)
```

No individual-axis matching is claimed. Different orthogonal bases
of the same span can leave a non-Clifford Boolean cubic phase.

### 9.2. Extract the aligning Clifford while keeping a Boolean residual

Let s be the sum of the literal T/T-dagger orientations in D_a.
Every signed phase axis in D_a anticommutes with P. The native identity
$`T_{-Q}=\omega T_Q^\dagger`$ yields

```math
P D_a P=\omega^sD_a^\dagger,\qquad
D_a^2=\omega^s P'P.                               \qquad(19)
```

The second operator is a scalar Pauli. Independently, D_a squared is
a product of commuting S rotations, so this Pauli has a fourth-root
scalar. Choose a literal diagonal Clifford
$`C_1\in\mathrm{Alg}(W)`$ with $`C_1^2=D_a^2`$: take an
eighth-root scalar square root and an S rotation on the Hermitian
Pauli label of D_a squared. The identity-label case is included.
All chosen scalars are actual Clifford scalars.

C_1 commutes with D_a. Hence
$`H_1=C_1^\dagger D_a`$ satisfies $`H_1^2=I`$ exactly and is
a Boolean commuting C3 phase on W. It sends P to a Pauli, so its
Boolean derivative along the nonzero P translation in W coordinates
is affine. Choose those coordinates with that translation first.
Writing the derivative as an affine function l of the other coordinates,
the real quadratic phase $`C_2=(-1)^{x_1l}`$ has the same derivative.
Thus

```math
C_0=C_1C_2,\qquad H_2=C_0^\dagger D_a
```

obey

```math
C_0PC_0^\dagger=P',\qquad H_2^2=I,\qquad
H_2\in\mathrm{Alg}(W\cap P^\perp)=\mathrm{Alg}(R_F).             \qquad(20)
```

H_2 is a Boolean cubic phase. Since R_F is a proved shared
computational center, this is a cubic sign of affine actual input
parities. It commutes with every remaining factor and endpoint and
can be moved literally to the input.

C_0 also commutes with every B_0 axis. Absorb it into the output
encoder by

```math
C_s'=C_0^\dagger C_sC_0,\qquad K_C'=C_0^\dagger K_C.
```

After retaining the central input sign H_2, the remaining full-space
Boolean operator has the exact form

```math
V H_{2,\mathrm{in}}^\dagger
 =(K_C')^\dagger(C_s')^\dagger B_0 A K_A.          \qquad(21)
```

Its three list relation dimensions are $`(0,\le1,0)`$. The
retained [one-relation three-band argument](PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md)
gives Boolean sign degree at most 2048. One may also realize its
physical hypothesis directly: lift the sole possible middle relation
with Z on one clean qubit, choosing the lift coefficients to have odd
pairing with that relation. The lifted middle list is independent,
all bands preserve the clean Pauli, and their occupied action is
exactly (21). This is an actual Clifford synthesis with its own cost,
not an assertion that a dependent logical list was already a physical
band. No additional clean resource beyond that one is needed.

Adding back H_2's degree-three sign and the Section 7 quadratic sign
does not raise the bound. Every intermediate identity is full-space;
no unused-sector replacement or separately chosen control-fiber
extension is involved. The aligning Clifford and coordinate changes
have paid fixed Clifford words, and their actual inverses are reversed
and inverted words. The cubic phase removal is exact degree accounting,
not a free same-depth compilation of the original circuit. This closes
the stated subcase without changing the complete-frame frontier.

### 9.3. The one-relation active sublist leaves a precise extra direction

For an abstract four-factor word satisfying the displayed list bounds,
if B_a has one relation, that relation necessarily has even weight:
every B_a axis anticommutes with P. The retained rank formula then gives

```math
\mathrm{rank}\!\left(\sum_{b\in B_a}bb^{\mathsf T}\right)
 =|B_a|-2=|F|,\qquad
\dim\mathrm{span}(B_a)=|F|+1.
```

Thus the F span has codimension one in the active-middle span. A
linear functional vanishing on the F span selects precisely the
support of the unique B_a relation: the Gram identity forces the
sum of its selected axes to vanish, and the relation kernel is
one-dimensional. Its support length need not be bounded.

The same literal square-root and quadratic alignment leaves a Boolean
cubic supported on R_F plus one additional commuting Pauli direction H.
In actual center coordinates it has the form

```math
(-1)^{a(z)}H^{b(z)},\qquad \deg a\le3,\quad\deg b\le2.           \qquad(22)
```

H commutes every remaining middle axis. The central cubic sign is
known, but the quadratically controlled Pauli still requires a global
matching argument. Equation (22) is not included in the conditional
closure above, and setting its controls to constants would not supply
a free full-space native extension.

### 9.4. The actual clean origin makes the active sublist independent

For the actual no-pair clean-derived representation, the independence
hypothesis in the conditional theorem is automatic whenever its
noncentral common parity exists. This holds at every non-Pauli clean
active rank, including rank two. Consequently that actual-source
subcase satisfies the abstract theorem above. The additional physical
split below sharpens its sign bound to five, or its XOR bound to four
when the output-Hadamard sign representation has the stated witness.

Transport the normalized common parity back through the actual fixed
Clifford changes and quadratic extraction. This gives input and output
logical Paulis with exact occupied intertwining. They need not both be
computational Paulis in the original coordinates. Their endpoint-span
memberships give signed physical lifts t_A in the original A-axis span,
commuting with P_A, and t_C in the original C-silent span. The input
lift comes from the merged A list and therefore commutes the chosen
A-active q. A lift may differ from the preferred endpoint representative
by the signed clean stabilizer; it has the same occupied action.

Set

```math
J_t=E_Bt_AE_B^\dagger,\qquad W=t_CJ_t.
```

Both J_t and t_C are actual Cliffords commuting with S. Their occupied
agreement makes W fix the positive S code pointwise. The
[pointwise clean-code stabilizer theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#1-the-exact-pointwise-clifford-stabilizer)
therefore gives W equal to I or S, apart from its two possible rank-two
exceptions. Those exceptions cannot occur here. In their canonical
coordinates the S displacement image is
$`L_{C,\mathrm{active}}=\langle Z_v,Z_t\rangle`$,
while the two exceptional displacement planes are
$`\langle Z_v,Y_t\rangle`$ and $`\langle Z_v,X_t\rangle`$.
Neither is contained in $`L_{C,\mathrm{active}}^\perp`$.
Every native B axis does lie in this orthogonal space under the
no-pair premise, so the displacement image of J_t, and hence of W,
lies there. This excludes both exceptions.

If W were I, the one-band Pauli image J_t would be Pauli. Native
independence forces its active B list to be empty, giving
$`t_A=t_C`$ literally. This label commutes both q and every original
C-active axis, contradicting its noncentral pairing with F. Hence

```math
E_Bt_AE_B^\dagger=t_CS.
```

Equality of intrinsic displacement images now gives
$`\mathrm{span}(B_{\mathrm{active}})=L_{C,\mathrm{active}}`$.
Equality of their Pauli-support cosets also gives
$`t_A+t_C+p_C\in L_{C,\mathrm{active}}`$ at the label level.
Pairing this identity with every original C-active axis shows that
t_A pairs to one with each such axis. Its pairing functional on this
span is therefore the clean functional chi, so every active B axis is
clean-negative.

The physical active B list is independent. The map tau is injective
on its span: in the coupled interface q cannot belong to the original
C-active span, because some C-silent axis pairs to one with q.
Furthermore tau preserves the active subset, since t_A commutes q.
The subsequent clean projection is injective as in the rank-one
interface proof. Thus the logical B_a in Section 9 is independent,
and Sections 9.1--9.2 apply with no additional active-list assumption.
Section 9.3 remains a possible obstruction for the abstract four-factor
class, but not for this actual-source noncentral-parity branch.

There is more information before projection. Write B_a for this
physical active sublist and B_r for the complementary clean-negative
sublist. The original no-pair clean Gram identity and the comparison
just proved give

```math
N_C=N_A+N_{B,-},\qquad N_{B_a}=N_C,
\qquad N_{B_r}=N_A.
```

Both B_r and the original A-active list are physically independent,
so their spans agree. Thus the actual clean phase splits as

```math
F_0=D_CD_A,\qquad
D_C=E_{C,\mathrm{active}}^\dagger E_{B_a},\qquad
D_A=E_{B_r}E_{A,\mathrm{active}}.
```

The two commuting factors are supported respectively on the original
C-active and A-active spans. Each sends P_A to a Pauli: its active
Gram sum vanishes. Apply the literal square-root and real quadratic
alignment of Section 9.2 separately in those spans. It gives actual
commuting Cliffords C_C,C_A and Boolean cubic phases H_C,H_A with

```math
D_C=C_CH_C,\qquad D_A=C_AH_A,
\qquad
\begin{aligned}
H_C&\in\mathrm{Alg}(L_{C,\mathrm{active}}\cap P_A^\perp),\\
H_A&\in\mathrm{Alg}(L_{A,\mathrm{active}}\cap P_A^\perp).
\end{aligned}
```

The symmetric Boolean reality argument of Section 10 makes both even
active spans a shared computational center. Their product is therefore
a literal cubic input sign, after the paid quadratic and fixed Clifford
normalization. Retain that sign and remove it only for degree accounting.
C_C commutes the original C-silent band, C_A commutes the original
A-silent band, and both commute every positive B axis. The residual
physical word before those quadratic and coordinate dressings is exactly

```math
(C_C^\dagger K_C)^\dagger
 E_{C,\mathrm{silent}}^\dagger E_{B,+}E_{A,\mathrm{silent}}
 (C_AK_A).                                             \qquad(23)
```

Its two endpoints map the clean Z to the same signed Pauli

```math
P_*=C_AP_AC_A^\dagger=C_C^\dagger P_CC_C.
```

Every displayed physical band preserves P_*. Their restrictions are
independent logical native lists. Indeed the A-silent and B-positive
spans commute an actual A-active axis q_A, which remains clean-odd
for P_* because C_A fixes it. The C-silent span commutes an actual
C-active axis q_C, which is clean-odd for P_* because C_C fixes it.
None of these spans can contain P_*, so projection introduces no
relation. Fixed Clifford clean coordinates make the endpoint
restrictions actual logical Cliffords. Equation (23) is therefore an
actual zero-clean three-native-band logical circuit, including every
arbitrary dirty input and its literal scalar.

Apply the same paid Section 7 quadratic adjustment and fixed coordinate
changes to the endpoints of (23). In these dressed coordinates the
removed central factor is the computational cubic sign just identified,
and the residual is a Boolean diagonal circuit. The bare original-word
expression (23) need not itself be computationally diagonal before
these dressings. They preserve all three native lists and add only
actual Clifford operations.

The retained zero-clean theorem bounds the dressed word's Boolean sign
degree by five.
Restoring the cubic input sign and the quadratic normalization leaves
that bound unchanged. This proves the stated five/four consequence for
the actual-source noncentral-parity subcase. The additional aligners
have paid $`O(Q^2)`$ elementary Clifford words, no new helper is used,
and every actual inverse is obtained by reversing and inverting its
word. The cubic sign was retained in an exact target factorization;
its removal is not free same-depth synthesis of the original target.

### 9.5. The complementary common parities have literal shared lifts

Retain the actual no-pair clean-source setting and normalized common
endpoint space T, without the noncentral-witness premise of Section 9.4.
Suppose every t in T has its input-frame
Pauli commuting with F. Then every such parity has a common signed
physical lift commuting with all three original native bands.

Use the transported lifts t_A,t_C from Section 9.4. The same pointwise
stabilizer argument, including the rank-two exclusion, gives

```math
W=t_C E_Bt_AE_B^\dagger\in\{I,S\}.
```

If W were S, its support-coset identity would make t_A pair to one
with every original C-active axis c. Since t_A commutes the chosen
A-active q, it would also pair to one with every projected F axis
c+q. This contradicts the premise. Therefore W=I, so the native
one-band image $`E_Bt_AE_B^\dagger`$ is Pauli. Physical B-axis
independence forces its active sublist to be empty. Literal conjugation
then gives

```math
E_Bt_AE_B^\dagger=t_A=t_C.
```

This common lift commutes E_A because it belongs to the original
A-axis span, commutes E_C because it belongs to the original C-silent
span, and commutes E_B by the empty active-sublist conclusion. The
whole lifted family is commuting and retains its signs. All these
identities hold on the full physical space. The input and output
parities were transported through the actual normalization; no claim
that they were originally computational parities is needed.

Thus the complementary branch places the Section 6 high-degree
residual in a genuine shared physical center. This alone supplies no
degree bound for that residual. The quartic quotient is still a target
factorization, and fixing shared-center eigenvalues can create an
unbounded number of relations among projected native axes. The next
matching problem must retain the global lists together with this
literal common-lift property.

## 10. The actual clean transition leaves a quadratic-controlled Pauli

The actual clean-derived transition admits a sharper form than the
general F factor in Section 7.1. Retain the no-middle-pair hypotheses
with both interior clean images non-Pauli. In the original common
isotropic active span, write

```math
L_0=L_A+L_C,\qquad \chi(v)=\langle v,p_A\rangle,\qquad
L_A^0=L_A\cap\ker\chi,\quad L_C^0=L_C\cap\ker\chi,
\qquad R=L_A^0+L_C^0.
```

Here L_A,L_C denote the original clean-active spans, not the endpoint
silent spans of Section 9. Boolean reality makes R a shared center.
For an element of L_C^0, the logical four-factor reduction commutes
with that Pauli in B, F and C_s; Section 7 forces A to commute with it
as well. Conversely, an element of L_A^0 commutes with A, B and F.
The time-reversed one-band reality argument forces C_s to commute
with it. The logical quotient is injective on these active directions,
so these are also commutation statements about the original physical
bands. Thus every actual band, and every commuting active-span aligning
Clifford, centralizes R.

If a clean-odd vector lies in L_A intersect L_C, it already supplies
the shared clean flip of the
[retained criterion](RANK_ONE_CLEAN_PHASE_INTERFACE.md#23-a-full-odd-coset-clean-flip-criterion-and-its-dual).
Suppose that criterion fails. Then L_A intersect L_C is clean-even,
and dimension counting gives

```math
\dim R=\dim L_0-2,\qquad
\ker\chi=R\mathbin\oplus\langle h\rangle.
```

One may take h to be the sum of one actual A-active axis and one
actual C-active axis. Their literal commuting product is a signed
Hermitian Pauli H. It commutes every positive middle axis.

### 10.1. Choose an aligner with a literal Boolean residual

Keep the actual clean-transition word

```math
F_0=E_{C,\mathrm{active}}^\dagger E_{B,-}E_{A,\mathrm{active}},
\qquad F_0P_AF_0^\dagger=P_C.
```

Every axis in this word anticommutes with P_A. The same literal
calculation as (19) gives

```math
F_0^2=\omega^sP_CP_A,
```

where s sums its signed T orientations, including the inverted C
orientations. This is a scalar Pauli. Since the square is also a
product of commuting S rotations, its scalar relative to a Hermitian
Pauli is a fourth root. Choose an actual Clifford C_1 in Alg(L_0)
with C_1 squared equal to F_0 squared, using an eighth-root scalar
and an S rotation on that Pauli label. The scalar is an actual
Clifford scalar; its native word and inverse are retained.

The commuting phase H_1=C_1-dagger F_0 is then literally Boolean C3.
It sends P_A to a Pauli. In intrinsic coordinates where P_A translates
the first coordinate, its Boolean derivative is affine and independent
of that coordinate. A real quadratic Clifford C_2 has exactly this
derivative. Consequently the actual commuting Clifford and residual

```math
C_0=C_1C_2,\qquad G_0=C_0^\dagger F_0
```

satisfy

```math
C_0P_AC_0^\dagger=P_C,\qquad
G_0^2=I,\qquad G_0\in\mathrm{Alg}(\ker\chi).
```

G_0 is a Boolean cubic in that commuting algebra, with its literal
constant sign. This choice of aligner avoids leaving a fixed T or S
rotation in the residual; it does not change the occupied action.

### 10.2. Make the controls computational and retain the native lists

Apply the simultaneous real quadratic extraction and fixed Clifford
normalization of Section 7 to the shared R. Its eigenbits z are now
genuine computational center coordinates; Hadamards were used only
on residual identity modes. Since the residual Boolean cubic is
supported on R plus H, it has the exact form

```math
G_0=(-1)^{a(z)}H^{b(z)},\qquad
\deg a\le3,\qquad \deg b\le2.
```

This is the unique Boolean expansion in the H eigenbit: its
coefficient has degree at most two. The central cubic sign commutes
with every remaining factor and both aligned endpoints and can be
moved literally to the input. The remaining insertion is one fixed
Pauli controlled by a quadratic function of actual input parities.
It still commutes each positive middle axis.

Both silent endpoint bands and the positive middle band remain in
the exact word. All three logical lists are independent here. The
A-silent and B-positive physical spans commute an actual A-active
clean-odd axis, and the aligned C-silent span commutes an actual
C-active clean-odd axis. None can therefore contain P_A, so clean
restriction introduces no relation. Their actual Clifford-coordinate
native realizations are retained.

Only actual Clifford aligners and coordinates, with their corresponding
polynomial Clifford ledgers, were substituted. The source of G_0 is
still the literal product $`C_0^\dagger F_0`$, and its inverse
is the actual reversed inverse. No initialized helper is added and
no quadratic-controlled Pauli is declared a free native gate. Its
conjugation through a silent native band can be a Hermitian Clifford
with unbounded displacement rank. Thus this sharper actual-input
normal form still requires a global matching theorem; it does not
close the unrestricted growing-rank degree question.
