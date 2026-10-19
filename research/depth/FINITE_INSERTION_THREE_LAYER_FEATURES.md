# Three commuting phase bands with a finite Pauli-supported insertion

[Bounded-relation native matching](BOUNDED_RELATION_NATIVE_MATCHING.md) ·
[Partial-diagonal features](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md) ·
[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md)

A bounded insertion can be handled without deleting it or commuting it
through a silent endpoint band. The additional hypothesis used here is
that it commutes with **each middle-axis Pauli**. Its effect on native
Pauli coefficient norms is bounded, allowing joint axis recovery with a
finite exceptional set. The remaining exceptional operation may be
non-Clifford, non-C3, and supported on a nonisotropic Pauli space.

All comparisons in this note are full-space identities. Application to
an initialized clean code requires a separate exact reduction producing
the stated full logical operator and the axiswise commutation premise.

## 1. The full-space expression and explicit bound

Suppose
```math
V=K_C^\dagger E_C^\dagger F E_BE_AK_A               \qquad(1)
```
is projectively a Boolean diagonal sign on its entire input space.
The endpoint K operators are actual Cliffords. Each E is a product of
literal T or T-dagger rotations about commuting signed Hermitian
Paulis. Their binary relation dimensions are bounded by r_A,r_B,r_C.

Let F be any fixed unitary satisfying
```math
F\in\mathrm{Alg}(S_F),\qquad \dim S_F\le d,
\qquad [F,B]=0\quad\text{for every middle-axis Pauli }B .
                                                               \qquad(2)
```
The support space S_F need not be isotropic. F need not belong to C3.
Commutation with the product E_B alone is not the premise in (2).

Identity axes contribute their actual scalars. Repeated axes can be
combined, with even powers absorbed as commuting Clifford corrections.
For the middle list, F commutes those corrections by (2). As in the
bounded-relation theorem, this reduces to distinct nonzero axes without
increasing the three relation bounds or changing (2).

Define
```math
\begin{gathered}
K=32(r_B+2d)^2,\qquad h=4r_B+7K,\\
\rho=r_B+d+h,\qquad w=d+h+4\rho^2,\\
k_*=r_A+w+r_C,\\
D_{\rm ins}(r_A,r_B,r_C,d)
   =36\,2^{k_*}+16k_*+w+84 .
\end{gathered}                                      \qquad(3)
```

**Finite-insertion degree theorem.** For
$`V=e^{i\theta}\mathrm{diag}((-1)^{g(x)})`$,
```math
\deg_{\mathbb F_2}g
       \le D_{\rm ins}(r_A,r_B,r_C,d).              \qquad(4)
```

Only the fixed parameters in (3) enter the bound. No width, number of
middle pairs, or number of central controls appears. Dependent phase
lists are algebraic factors, not automatically one physically native
T band. The theorem makes no synthesis cost claim for F.

## 2. Bounded norm distortion replaces exact rank-profile invariance

For a computational Z label z, put
```math
R_z=K_AZ(z)K_A^\dagger,\quad
J_z=E_AR_zE_A^\dagger,\quad
C_z=E_CK_CZ(z)K_C^\dagger E_C^\dagger .
```
The diagonal target and (1) give
```math
F E_BJ_zE_B^\dagger F^\dagger=C_z.                 \qquad(5)
```
Both J_z and C_z are actual Hermitian Cliffords. The J family is an
exact commuting representation. Its binary actions satisfy
```math
\Phi_z=I+N_z,\qquad
N_zv=\sum_i\lambda_i(z)\langle a_i,v\rangle a_i,
\qquad N_zN_u=0.                                    \qquad(6)
```
Peel globally silent A axes to the input as a diagonal C3 phase;
this leaves the J family unchanged.

Use normalized Hermitian Pauli coefficients and their $`\ell_1`$
norm. A unitary supported in a d-dimensional Pauli space has at most
$`2^d`$ coefficients and coefficient squared norm one. Hence
```math
\|F\|_{P,1}\le2^{d/2}.                              \qquad(7)
```
Pauli coefficient norm is submultiplicative, and Clifford conjugation
preserves it exactly.

For any Pauli P, set $`X=E_BPE_B^\dagger`$.
Equation (5) implies
```math
E_BJ_zPJ_z^\dagger E_B^\dagger
 =(F^\dagger C_zF)X(F^\dagger C_z^\dagger F).
                                                               \qquad(8)
```
The two expressions compared in (8) are Hermitian Cliffords because
E_B is a commuting C3 phase. The right side has two conjugations by
F or its inverse, with a Clifford conjugation between them. Its
Pauli coefficient norm is at most
$`\|F\|_{P,1}^4\|X\|_{P,1}\le2^{2d}\|X\|_{P,1}`$.
The reverse comparison has the same bound.

Let w(v) be the displacement rank of $`E_BPE_B^\dagger`$.
Clifford flatness gives its coefficient norm as $`2^{w(v)/2}`$.
Therefore
```math
|w(\Phi_zv)-w(v)|\le4d.                             \qquad(9)
```

For a commuting list of relation dimension r_B, the exact rank formula
from [bounded-relation matching](BOUNDED_RELATION_NATIVE_MATCHING.md#2-a-rank-profile-recovers-all-but-boundedly-many-axes-jointly)
is
```math
w(v)=\mathrm{wt}_{A_B}(v)-\delta(v),\qquad
0\le\delta(v)\le2r_B.
```
It includes all signed orientations. Thus
```math
|\mathrm{wt}_{A_B}(\Phi_zv)-\mathrm{wt}_{A_B}(v)|
       \le2r_B+4d.                                 \qquad(10)
```
Fourier orthogonality converts (10) into
```math
|A_B\mathbin{\triangle}\Phi_z(A_B)|
       \le16(r_B+2d)^2.                             \qquad(11)
```

The bound holds uniformly for every member of the entire finite J
label group, not only its generators. Strict-majority orbit membership
therefore gives an invariant isotropic proxy B with
$`|A_B\triangle B|\le K`$ and relation dimension at most
$`r_B+K`$. Exact rank-profile invariance is not asserted.

## 3. Localize the combined finite exceptional operation

The square-zero identities (6) make every proxy orbit an affine set
$`b+W_b`$. An orbit of size $`2^s\ge4`$ contributes
$`2^s-s-1`$ independent internal relations. The argument from the
bounded-relation theorem therefore bounds all large proxy orbits by
$`4(r_B+K)`$ axes.

Remove those orbits and all singleton or pair orbits meeting
$`B\setminus A_B`$. The remaining ordinary list consists of actual
middle axes, is invariant under the whole J family, and has only
singleton or pair orbits. All remaining actual bad axes number at most
```math
|A_B\setminus B|+4(r_B+K)+2|B\setminus A_B|
 \le h.                                             \qquad(12)
```
Write $`E_B=E_{\rm bad}E_{\rm ord}`$ and put
```math
E_*=F E_{\rm bad},\qquad
S_*=S_F+\mathrm{span}(A_{\rm bad}),\qquad
\dim S_*\le d+h.                                    \qquad(13)
```
E_* is not assumed commuting C3. It does commute every ordinary axis,
because F does so and all actual B axes commute.

The signed permutation of ordinary axes makes
```math
D_{{\rm ord},z}
  =E_{\rm ord}J_zE_{\rm ord}^\dagger J_z^\dagger
```
a commuting Clifford on those axes. It commutes E_*. Equation (5)
then shows that $`E_*J_zE_*^\dagger`$ is Clifford independently.

Define the full-space intrinsic Pauli kernel
```math
\mathcal K_*=
 \{p:E_*P E_*^\dagger\text{ is a Pauli with its literal sign}\}.
                                                               \qquad(14)
```
It is a linear label space and is invariant under every Phi_z.
The latter follows by conjugating a Pauli image with the actual
Clifford $`E_*J_zE_*^\dagger`$. Moreover
$`S_*^\perp\subseteq\mathcal K_*`$, and every ordinary B axis
belongs to K_* and is fixed literally by E_*.

Actual conjugation gives a signed symplectic isometry on K_*.
Extend its binary isometry to a full symplectic automorphism, choose
a Clifford representative, and correct its residual sign character
on K_* by a Pauli. This supplies an actual Clifford C_* matching
E_*'s exact Pauli action on every member of K_*.
In particular C_* fixes $`S_*^\perp`$ and every ordinary B axis
pointwise. The fixed commutant implies
```math
C_*\in\mathrm{Alg}(S_*),\qquad
[C_*,B]=0\quad(B\text{ ordinary}).                  \qquad(15)
```
Thus membership in the same finite algebra is a consequence of the
literal fixed-Pauli action, not an assumption about a chosen extension.

Set
```math
E_{\rm eff}=C_*^\dagger E_*,\qquad
W=\mathcal K_*^\perp.
```
Then
```math
E_{\rm eff}\in\mathrm{Alg}(W),\qquad
W\subseteq S_*,\qquad
\dim W\le d+h.                                     \qquad(16)
```
W is invariant under the complete J family and commutes every ordinary
axis. It need not be isotropic. The unitary E_eff need not be C3.
Its scalar remains literal.

Absorb C_* into the C encoder. Since it commutes the ordinary band,
the middle expression is now exactly $`E_{\rm eff}E_{\rm ord}`$,
and the C list remains a Clifford conjugate of its previous list with
the same relation bound.

## 4. Pair peeling does not require an isotropic exceptional space

The ordinary singleton and pair characters obey the same exact signed
group laws as in bounded-relation matching. For a pair, denote its
swap character by lambda, its common sign character by chi, and its
axis difference by r. The Hermitian product of its two axes is fixed
by every J_z. As before, r is an A-span direction and pulls back to
a computational input Z parity.

The exceptional defect
```math
D_{{\rm eff},z}
 =E_{\rm eff}J_zE_{\rm eff}^\dagger J_z^\dagger
```
is an actual Clifford in Alg(W). It fixes every Pauli in
$`W^\perp`$ pointwise. Since J_z preserves W, its binary action
preserves $`W^\perp`$. Ordinary axes commute W, so their defect
also preserves that domain.

Restrict the displacement identity to inputs in $`W^\perp`$ and
project its outputs to the ambient label quotient $`\mathcal V/W`$.
The pairing
```math
(\mathcal V/W)\times W^\perp\longrightarrow\mathbb F_2
```
is perfect, whether or not W is isotropic. The exceptional defect
vanishes on this domain. The only potentially quadratic contribution
is consequently
```math
\sum_j\chi_j(z)\lambda_j(z)\,
                   \bar r_j\bar r_j^{\mathsf T}.    \qquad(17)
```
No quotient $`W^\perp/W`$ is used.

The projected difference list has relation dimension at most
$`r_B+\dim W\le\rho`$.
The rank-one tensor map has kernel dimension at most rho.
Every tensor relation has support at most $`2\rho`$ by
Sylvester's inequality, so a kernel basis uses at most
$`2\rho^2`$ pair coordinates. Outside that union a dual
functional isolates each tensor in (17), forcing chi lambda to have
no degree-two component. Thus chi is zero or lambda.
Orient and peel those pairs with the retained literal scalar identities.

Adjoin the axes of the at-most-$`2\rho^2`$ remaining pairs to W
and include their phase factors in the exceptional operation. Call
the enlarged space W' and the resulting unitary E'. Then
```math
\dim W'\le d+h+4\rho^2=w.                            \qquad(18)
```
W' is still J-invariant and commutes all remaining middle axes.
E' is an arbitrary fixed unitary in Alg(W'), with no C3 claim.

Every peeled pair phase is diagonal on the full logical input space.
Its exponent has ordinary degree at most five modulo eight, by the
retained parity-sector Hermitian-Clifford argument. The argument uses
the A band being C3; it makes no assumption about the removed
exceptional factor.

## 5. A common bounded quantum algebra for all residual factors

The source activation subspace has full coordinate support after
silent A peeling. For each remaining singleton b, its pairing vector
with the A axes, multiplied coordinatewise by that activation subspace,
lies in the source relation kernel. The
[product-space rank lemma](BOUNDED_RELATION_NATIVE_MATCHING.md#5-relations-also-bound-singleton-coupling-and-hidden-c-directions)
therefore bounds the A/singleton cross-pairing rank by r_A.

Their joint span has symplectic rank at most $`2r_A`$.
Adjoining W' gives a space S of rank at most $`2(r_A+w)`$.
The actual residual source
$`U=E'E_{B,\rm single}E_{A,\rm active}`$
belongs to Alg(S), although its factors need not commute with one
another. The full-space comparison still says $`UR_zU^\dagger`$
is a Clifford with the C encoder's native displacement.

The Clifford $`(UR_zU^\dagger)R_z^\dagger`$ belongs to Alg(S),
so its displacement vanishes on $`S^\perp`$.
Applying the same product-space rank lemma to the active C list
shows that its span escapes S in at most r_C directions.
After peeling silent C axes, all residual factors belong to one
space $`\widehat S`$ with
```math
\mathrm{rank}_{\rm symp}\widehat S
          \le2(r_A+w+r_C)=2k_* .                   \qquad(19)
```

There are three commuting C3 factors in this common algebra:
active A, singleton B, and active C. The fourth factor is the
arbitrary unitary E' with support dimension bounded by w.
The central radical of the common algebra may have unbounded dimension.

## 6. Exact matrix features and Boolean phase extraction

The residual diagonal operator has the form
$`V_0=K_C^\dagger W_0K_A`$, where W_0 is the ordered product
of the four residual factors. It need not be a real sign.
On $`K_A^\dagger\widehat S^\perp K_A`$, V_0 acts on Paulis
literally. Extract a diagonal Clifford C with that signed action.
Then
```math
J=K_CCK_A^\dagger\in\mathrm{Alg}(\widehat S),\qquad
V_1=C^\dagger V_0=K_A^\dagger J^\dagger W_0K_A.
                                                               \qquad(20)
```
V_1 is diagonal. Its Pauli support lies in the conjugate of
$`\widehat S\cap K_A\langle Z\rangle K_A^\dagger`$.

Use canonical central Z coordinates and k quantum coordinates for
Alg($`\widehat S`$), with $`k\le k_*`$.
The final isotropic support has commuting core projections, diagonalized
by one fixed core Clifford. Its target eigenbits are affine parities
of the original computational input.
Choose the exact affine section: prescribe the central target eigenbits,
set complementary central coordinates to zero, set each selected
core bit to its target eigenbit plus its central parity offset,
and set the other core bits to zero.
Pauli-support uniqueness proves that this section evaluates the same
diagonal phase. The canonical Clifford need not preserve all other
input computational Paulis.

Each of the three C3 factors is diagonalized in the core by a fixed
Clifford. At every fixed core string its central exponent is weighted
cubic modulo eight, whose three binary digits have degree at most four.
These factors require at most $`9\,2^k`$ such features.

For E', choose a Pauli basis of W', of size at most w.
In the common canonical coordinates, each basis Pauli has a fixed
core label and one central Z parity. Its scalar normalization is fixed.
The exact matrix E'(z) is therefore determined by at most w affine
central parity bits. This uses its Pauli expansion in Alg(W'), not
a C3 tableau or a commuting-core diagonalization of E'.

Finally the Clifford J has fixed core Clifford, affine ordered-X/Z
labels, and a central diagonal-Clifford scalar. Its labels and scalar
digits, together with the central offsets and core entry selector,
use at most $`4k+3`$ features of degree at most four.
Every fixed scalar is retained.

The aggregate peeled phase has ordinary exponent degree at most five
modulo eight. Its three digits have degree at most twenty.
The extracted diagonal Clifford C adds three digits of degree at most
four. Once those phases are restored, the original Boolean target
makes its sign a Boolean function of all the listed features.
Extend that function arbitrarily off the attained feature vectors.
The total degree is bounded by
```math
\begin{aligned}
\deg g
 &\le4(9\,2^k+4k+3)+w+3\cdot20+3\cdot4\\
 &=36\,2^k+16k+w+84\\
 &\le D_{\rm ins}(r_A,r_B,r_C,d).
\end{aligned}                                      \qquad(21)
```

Arbitrary intermediate matrix entries were not assumed Boolean or
roots of unity. Only the final target supplies the sign condition.
All inverses and common phases are literal throughout. No step deletes
the finite insertion or attributes bounded codimension to the input
constraints produced by a silent-band conjugation.

## 7. Scope and resource accounting

The fixed-space statement follows from the entire logical comparison
family and axiswise commutation in (2). It does not follow from the
insertion having a small support alone, or from an occupied comparison
being extended arbitrarily to the unused sector. Every Clifford extracted
here matches an actual signed Pauli action on a full-space Pauli subgroup.

The feature coordinates describe existing operators. They are not free
physical coordinate circuits. A dependent commuting phase list is not
credited as one native layer, and every inserted or extracted Clifford
and its actual inverse must be charged in a synthesis ledger. The
physical clean-code application below supplies the required logical
expression without changing the occupied action or allocating workspace.

The [exact finite-insertion fixtures](../../tests/test_finite_insertion_three_layer_features.py)
check a case where the insertion permits a genuine failure of raw
native rank-profile invariance, and a non-C3 insertion whose effective
support is nonisotropic. They retain the physical native words,
their complete column actions, and their actual inverses.

## 8. One clean qubit with bounded active ranks, including middle pairs

Consider an actual one-clean native three-band circuit
```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A
```
with complete clean and dirty return and occupied Boolean sign action.
Use the actual clean Paulis and interior images
```math
\begin{gathered}
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger,\\
J=E_AP_AE_A^\dagger,\qquad
S=E_CP_CE_C^\dagger=E_BJE_B^\dagger .
\end{gathered}                                      \qquad(22)
```
Suppose both interior images are non-Pauli. Let L_A,L_C be the spans
of their clean-active endpoint axes, of dimensions k_A,k_C. Then
```math
\deg g\le D_{\rm ins}(1,1,1,k_A+k_C-1).              \qquad(23)
```
No hypothesis on the number or type of middle pairs is added.
The cases with a Pauli interior image retain their stronger
[previous bounds](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md).

### 8.1. Localize the entire clean transition

Apply actual independent-axis native matching to the clean probe (22).
The middle orbits are signed singletons and pairs. Write the pairs as
$`\{b_j,d_j\}`$, and their differences as r_j=b_j+d_j. The r_j are
independent and lie in both L_A and L_C: the native displacement
formulas for J and S give this for the same pair, since conjugation
by the commuting middle band does not change its permutation.
If there are m pairs, put
```math
S_0=L_A+L_C+\mathrm{span}\{b_j,d_j:1\le j\le m\}.
```
Since each d_j=b_j+r_j and the intersection contains the m independent
differences,
```math
\dim S_0\le(k_A+k_C-m)+m=k_A+k_C.                  \qquad(24)
```
S_0 need not be isotropic.

Every negative singleton axis belongs to L_A+L_C. To see this, the
literal defect
```math
D=E_BJE_B^\dagger J^\dagger=SJ^\dagger
```
is a commuting Clifford on the independent B axes: signed permutation
changes each odd phase exponent by an even number. If J and S have
binary actions I+N_A and I+N_C, then
```math
N_D=N_C+N_A+N_CN_A,
\qquad\mathrm{im}\,N_D\subseteq L_A+L_C.
```
A negative singleton contributes an odd S power to D. For independent
commuting axes, the image of that Clifford displacement is precisely
the span of its odd-S axes. This proves the inclusion.

Let E_{B,act} contain every pair and every negative singleton; let
E_{B,+} contain the positive singletons. A singleton fixed as a binary
label by J commutes every clean-active A axis, by independence of
that endpoint list. The same argument with S gives commutation with
every clean-active C axis. All B axes commute each other. Consequently
every B-positive axis commutes S_0, and its positive signed action
under J also makes it commute P_A.

The actual clean-transition word
```math
F_0=E_{C,act}^\dagger E_{B,act}E_{A,act}
       \in\mathrm{Alg}(S_0)
```
satisfies $`F_0P_AF_0^\dagger=P_C`$, with the literal signs.
It may be non-C3. No single commuting phase expansion is assumed.

### 8.2. Extract a Clifford without flattening the residual

The actual signed action of F_0 fixes every Pauli in $`S_0^\perp`$
and sends P_A to P_C. On
```math
\mathcal K_0=S_0^\perp+\langle p_A\rangle
```
it is therefore a signed symplectic isometry. Extend it to an actual
Clifford C_0 by the signed-extension argument of Section 3. Matching
the fixed Paulis in $`S_0^\perp`$ forces
```math
C_0\in\mathrm{Alg}(S_0),\qquad
G_0=C_0^\dagger F_0\in
       \mathrm{Alg}(S_0\cap p_A^\perp).             \qquad(25)
```
The second inclusion uses the exact commutation of G_0 with P_A.
At least one active A axis in S_0 anticommutes with P_A, so
```math
\dim(S_0\cap p_A^\perp)\le k_A+k_C-1.              \qquad(26)
```
Every B-positive axis commutes C_0 and G_0.

Set
```math
K_C'=C_0^\dagger K_C,\qquad
E_{C,sil}'=C_0^\dagger E_{C,sil}C_0.
```
Reordering only the stated commuting factors gives the literal
full-space identity
```math
V=(K_C')^\dagger(E_{C,sil}')^\dagger
          G_0 E_{B,+}E_{A,sil}K_A.                  \qquad(27)
```
All internal factors now preserve P_A, and both endpoints send Z_c
to P_A. Both silent endpoint factors have been retained. The C list
has merely been conjugated by an actual Clifford, preserving its
physical independence and all signs.

Use fixed Clifford coordinates for the positive P_A code and restrict
(27) to that code. The two endpoint restrictions are actual logical
Cliffords. Each independent commuting silent or positive sublist
projects to a logical list with at most one relation: the kernel of
the Pauli quotient $`p_A^\perp\to p_A^\perp/\langle p_A\rangle`$
has dimension one. Identity axes, repeated images and their scalars
are retained. The restricted G_0 is an arbitrary fixed unitary with
Pauli support dimension at most k_A+k_C-1. It commutes **each**
projected middle axis. Thus the full logical expression satisfies
Sections 1--6 with parameters $`(1,1,1,k_A+k_C-1)`$, proving (23).

For an XOR query, Hadamards on its output turn the occupied action
into $`(-1)^{zf(x)}`$. They commute with Z_c, so the two interior
clean images, their active ranks and the clean matching are unchanged.
Therefore
```math
\deg f\le D_{\rm ins}(1,1,1,k_A+k_C-1)-1.           \qquad(28)
```
The retained
[R74 rounding consequence](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
transfers this scoped obstruction below occupied operator error
$`2^{-530690}`$, allowing an optimized common phase. It is applied
to the original native three-band circuit, not to a circuit with a
freely supplied G_0 gate.

### 8.3. Resource and model scope

The construction above is an exact algebraic reduction for a degree
proof. It neither deletes G_0 nor counts it as a native T band.
Its literal source is C_0-dagger followed by the actual F_0 word;
its actual inverse reverses that word and conjugates every scalar.
Physical Clifford coordinate changes and C_0 can be synthesized with
$`O(Q^2)`$ elementary Cliffords on the same Q wires if needed, and
their inverse costs must also be charged. No helper is initialized
or added. The original clean leakage/return, every logical column,
arbitrary dirty inputs and external references are preserved.

This closes all fixed endpoint clean-active ranks, including the
general paired case and rank two. It does not give a degree bound
uniform in growing active ranks, an unrestricted one-clean three-band
theorem, or a new upper resource frontier. The
[easy-transition counterfamilies](CLEAN_TRANSITION_PARAMETERS.md)
explain why replacing G_0 by a single bounded-relation phase list
cannot supply that missing uniformity.
