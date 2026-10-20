# A complete original-pair bound when the input silent list commutes with the middle

[Original physical pairs](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[Larger-support sector revision](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#8-larger-physical-support-a-two-band-sector-encoder-and-a-one-band-obstruction) ·
[Occupied Boolean center](OCCUPIED_BOOLEAN_PAIR_CENTER.md)

The larger-support sector code can require two native endpoint bands.
That encoding obstruction does not prevent a complete degree theorem
when the original lists supply a different matching mechanism. Here
the Boolean condition either supplies a full physical clean flip or
makes the entire retained clean transition commute with an enlarged
logical Pauli subgroup. No parity sector is compiled separately.

## 1. Statement in the original physical lists

Use an actual one-clean three-native-band word

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,\qquad
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger.
```

Each phase list consists of independent commuting signed Hermitian
Pauli axes with their literal T or T-dagger orientations. Assume exact
return of the sole clean wire and the occupied action

```math
G=e^{i\theta}\mathrm{diag}((-1)^{g(x)})
```

on the entire logical and dirty register. The scalar is fixed; the
unused clean block is retained and need not be Boolean. Suppose:

1. The original middle list has exactly one genuine paired orbit, with
   axes B,D and difference label $`r=b+d`$.
2. This r is the sum of **all actual clean-active axes** at A, and also
   the sum of all actual clean-active axes at C.
3. Every actual input-silent A axis commutes with every original middle
   axis.

Then

```math
\deg g\le2428,\qquad
\deg f\le2427\quad\text{for a complete Boolean XOR query }Q_f.
                                                               \qquad(1)
```

The time-reversed version of condition 3 also suffices: every actual
output-silent C axis may instead commute with the entire original
middle list. In particular, an empty silent list at either endpoint
suffices. The result permits both active ranks, both silent spans,
the middle singleton list and the dirty width to grow without bound.
The active-support conditions are statements about physical lists,
not a choice of basis after quotienting. The common active rank is
even; the new larger-support cases have rank at least four. The
existing rank-two theorem remains much sharper in its overlap.

This is a source-class Boolean-degree bound. It does not settle an
arbitrary original pair, fixed pair count in general, or unrestricted
depth. The original-circuit rounding consequence in Section 7 gives
the retained conservative higher-degree error floor.

The class also falls qualitatively under the earlier bounded-coupled-span
theorem after reversal: the matching proved below makes the original
input-silent list centralize S_0, so the reversed coupled-silent sublist
is empty and t=0. The present result supplies an explicit criterion in
the original lists, its complete matching mechanism and the sharper
bound (1). It does not enlarge the earlier theorem's fixed-(s,t) scope.

## 2. Orientation repair removes every negative singleton

Write $`J=E_AP_AE_A^\dagger`$ and $`JBJ=sD`$, where s is a sign.
For the original pair orientations epsilon,delta define the actual
operators

```math
F=T_B^\epsilon T_D^{s\epsilon},\qquad
C=T_D^{\delta-s\epsilon},\qquad E_B=CFE_0.
```

Here C is Clifford, F commutes J and every singleton, and E_0 retains
all original singleton axes and orientations. Absorb C into the output
encoder and output phase list exactly as in the
[physical insertion identity](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md#1-keep-the-actual-middle-axes-and-their-signs).
Use C-side notation for this corrected data below. C fixes r, so the
assumed sum of the actual output-active labels is still r. No native
list loses independence and no scalar is discarded. The clean
comparison becomes

```math
S=E_0JE_0^\dagger=E_CP_CE_C^\dagger.               \qquad(2)
```

Let $`B_-`$ be the actual negative singleton list under J. Every
singleton u has $`(F_J-I)u=0`$. Independence of the input-active
list therefore forces u to commute with every input-active axis.
For a negative singleton the literal identity is

```math
T_U^\epsilon J T_U^{-\epsilon}
 =\omega^{-\epsilon}S_U^\epsilon J,
\qquad \omega=e^{i\pi/4}.
```

The displacement cross terms vanish, since these singleton axes
commute the active span. For any independent endpoint active list A,
the diagonal quadratic expression of its displacement is

```math
\langle v,(F_J-I)v\rangle
 =\left\langle\sum_{a\in A_{\rm active}}a,v\right\rangle.
```

Applying this at both ends of (2) gives

```math
r=r+\sum_{u\in B_-}u.
```

The original middle axes are independent, so $`B_-`$ is empty. Thus
E_0 consists only of positive singletons and (2) gives $`S=J`$
literally. The intrinsic displacement images of the two independent
endpoint lists give

```math
L_A=L_C=L.
```

Their complete nonzero Pauli supports are the affine spaces
$`p_A+L`$ and $`p_C+L`$, so $`p_A+p_C\in L`$. Moreover
$`(F_J-I)b=r`$ and the all-active support hypothesis imply

```math
\langle b,a\rangle=1\quad(a\in A_{\rm active}),\qquad
\langle b,l\rangle=\langle p_A,l\rangle\quad(l\in L).
                                                               \qquad(3)
```

Write B_+ for the retained singleton list. Its axes commute L and
P_A. The corrected input/output silent lists commute L and P_A too;
for the output use $`p_A+p_C\in L`$.

## 3. Boolean transport and an exhaustive linear dichotomy

The [endpoint commutation lemma](OCCUPIED_BOOLEAN_PAIR_CENTER.md#1-actual-source-and-the-endpoint-commutation-lemma)
uses the complete occupied Boolean action. An A-silent axis commutes
P_A, the full A list and, by condition 3, the full middle list.
The lemma therefore makes it commute the full corrected C list.
Every positive singleton commutes P_A and the A-active list; condition
3 also makes it commute the A-silent list. It commutes the full middle
list, so the same lemma makes it commute every C-silent axis.
Consequently the three retained lists A-silent, B_+ and C-silent are
pairwise commuting, including their actual signed representatives.

Let M be the binary span of all three complete corrected physical
phase lists. Equivalently,

```math
M=L+\langle b\rangle+L_{A,s}+L_{B,+}+L_{C,s}.       \qquad(4)
```

**Case I: $`p_A\notin M`$.** Nondegeneracy of the ambient symplectic
form supplies a label h orthogonal to M with
$`\langle h,p_A\rangle=1`$. Since $`p_A+p_C\in L`$, its actual
Hermitian Pauli H anticommutes both clean Paulis and commutes every
original retained phase axis, including both members of the pair.
The full physical word transports the two Pauli clean flips
$`K_A^\dagger HK_A`$ and $`K_C^\dagger HK_C`$ exactly.

The [clean-controlled Pauli dressings](RANK_ONE_CLEAN_PHASE_INTERFACE.md#23-a-full-odd-coset-clean-flip-criterion-and-its-dual)
fix every clean-zero vector and align both flips to X_c. Together
with exact clean return, this gives the full-space operator
$`I_c\otimes G`$. All three original physical phase lists stay native
and unchanged. The zero-clean three-layer theorem yields the stronger
sign/XOR bounds five/four in this case.

**Case II: $`p_A\in M`$.** By (3), pairing any actual A-active axis
with an expression for p_A in (4) forces the coefficient of b to be
one. Thus

```math
p_A=b+l+u+v+w,
\quad l\in L,\ u\in L_{A,s},\ v\in L_{B,+},\ w\in L_{C,s}.
                                                               \qquad(5)
```

Pair (5) with any actual C-silent axis. Its pairings with p_A,L and
its own silent span vanish. Its pairings with u and v vanish by the
Boolean transport above. Hence it commutes b, and also d=b+r.
Input-silent axes commute the pair by hypothesis; positive singletons
commute it by original middle commutation. Thus every axis of all
three retained lists centralizes

```math
S_0=L+\langle b\rangle.                          \qquad(6)
```

In particular this proves an empty actual coupled-silent sublist,
not merely a bounded pairing rank or a bounded span assumption.

## 4. Enlarge the shared center before restricting the clean wire

Continue in Case II. The genuine center of the occupied-center
reduction is

```math
R=L\cap p_A^\perp.
```

Equation (3) makes it central in S_0. The quotient S_0/R is a
two-dimensional symplectic space. Put

```math
H_0=S_0\cap p_A^\perp.
```

This H_0 is isotropic and $`\dim(H_0/R)=1`$. Explicitly, if
$`\chi(l)=\langle p_A,l\rangle`$ and
$`\beta=\langle p_A,b\rangle`$, its elements are
$`l+e b`$ with $`\chi(l)=e\beta`$; the pairing of any two such
elements vanishes by (3).

Use the retained exact clean-transition localization, with its actual
Clifford $`C_0\in\mathrm{Alg}(S_0)`$ and

```math
G_0=C_0^\dagger E_{C,\rm active}^\dagger F E_{A,\rm active}
   \in\mathrm{Alg}(H_0),\qquad
K_C'=C_0^\dagger K_C.
```

Every retained silent/positive axis commutes S_0, so C_0 fixes those
axes literally. The full physical word is therefore

```math
V=(K_C')^\dagger E_{C,s}^\dagger
                G_0E_{B,+}E_{A,s}K_A.             \qquad(7)
```

All its internal factors preserve P_A. All commute every Hermitian
Pauli in H_0: the retained lists centralize S_0, and G_0 lies in the
commutative algebra of H_0. Equation (7) therefore transports the
entire H_0 Pauli subgroup between the two endpoint encoders on the
**full physical space**.

Restricting this identity to the complete positive clean code gives
actual logical Pauli transports. Quotient only the possible kernel
$`H_0\cap\langle p_A\rangle`$. Their Boolean derivatives are affine,
so the [shared computational-coordinate normalization](COMPUTATIONAL_PHASE_DERIVATIVES.md#7-boolean-reality-gives-a-shared-computational-center)
applies to this enlarged commuting subgroup. Extract its one quadratic
Boolean Clifford, use fixed affine coordinates, and apply Hadamards
only on factors where the residual target is already identity. The
logical image of H_0 is now computational Z at both endpoints.

This is the same global transport argument used for R, applied to the
larger subgroup whose transport was just proved. No independently
chosen coordinate word in different parity sectors is used. Its
alignment preserves the meaning of a degree conclusion in the original
input bits: the residual target is independent of the Hadamard-rotated
factors, and the remaining coordinate changes are affine.

## 5. A diagonal residual and the complete feature bound

In these aligned coordinates let F_out be the complete logical
restriction of $`(K_C')^\dagger G_0K_C'`$. It is diagonal, since
the entire logical H_0 subgroup is computational Z. Equation (7)
gives the literal full-logical-space factorization

```math
G=F_{\rm out}W_3,
\qquad W_3=(K_C')^\dagger E_{C,s}^\dagger E_{B,+}E_{A,s}K_A,
                                                               \qquad(8)
```

where the endpoints in W_3 are their actual clean-code restrictions.
One common fixed Clifford coordinate for P_A makes both restrictions
literal logical Cliffords.
The three separate projected lists remain independent by the existing
active clean-flip witnesses. No dependent union is assigned a native
band. Their union does commute, by Section 3. A product of commuting
Pauli T phases lies in $`\mathcal C_3`$, even if its full union has
relations: conjugating any Pauli leaves a product of actual commuting
S phases and that Pauli, which is Clifford. Clifford endpoints preserve
this membership. Thus W_3 is full-space $`\mathcal C_3`$. Since both
G and F_out are diagonal, W_3 is diagonal as well. Neither intermediate
factor was assumed Boolean.

The [retained coefficient ledger](OCCUPIED_BOOLEAN_PAIR_CENTER.md#9-a-clean-preserving-clifford-sum-without-bounding-the-silent-span)
for one original pair gives

```math
G_0=f_0(R)I+f_1(R)P_h,\qquad
N(1)=9\cdot2^6+12+15=603.                        \qquad(9)
```

Here h represents H_0/R and every coefficient together is determined
by at most 603 Boolean features of degree at most four. Literal signs
and the fixed scalar are retained. The enlarged normalization makes
the R eigenvalues affine computational parities, preserving these
degree bounds. The logical P_h is either identity or one additional
signed computational Z parity. Charge its eigenvalue bit as one more
feature. In particular, $`p_A\in M`$ alone does **not** make this bit
constant; that stronger simplification would require
$`p_A\in S_0`$.

A diagonal $`\mathcal C_3`$ unitary has, up to its fixed scalar, a
weighted-cubic exponent modulo eight. Its three binary exponent digits
have degrees at most one, two and four, respectively: carries among
linear terms have degree at most four, coefficients of quadratic terms
are even, and cubic coefficients are multiples of four. Charge three
more degree-at-most-four features for W_3.

The complete Boolean sign in (8) is therefore a Boolean function of
at most

```math
603+1+3=607
```

features of degree at most four. Extend this sign function arbitrarily
off its attained feature vectors and take its algebraic normal form.
This extends only a finite classical sign function, not an intermediate
quantum action. Composition proves $`\deg g\le4\cdot607=2428`$.
Restoring the extracted quadratic phase does not raise that bound.
Together with Case I this proves (1) for the entire occupied action,
including both pair-parity sectors simultaneously.

## 6. Reverse orientation and the XOR target

For the actual reversed inverse, the original middle axes commute
the original middle band. Their orbit action under S is consequently
the same signed singleton/pair action as under J. The endpoint active
lists exchange places, and the all-active sum condition is unchanged.
If the original C-silent list commutes the original middle list, the
inverse satisfies condition 3. Its Boolean sign is the same g, with
the conjugated common scalar, so the proof applies unchanged.

For a Boolean XOR query, apply output Hadamards to its actual endpoint
Cliffords. The clean Paulis, physical phase lists and all three source
hypotheses stay unchanged. Its sign polynomial is $`y f(x)`$, of
degree $`\deg f+1`$ for nonzero f. This gives the bound 2427 while
preserving every arbitrary output, dirty and reference input.

## 7. Resources, approximation and remaining scope

The source and its actual inverse retain their original T counts,
at most three native T bands, one clean qubit and the original dirty
width. Orientation correction, the actual active Clifford C_0 and
fixed global coordinate words cost $`O(Q^2)`$ elementary Cliffords
on Q existing wires if emitted; each inverse is charged as its actual
reversed word. The Case I clean-controlled Pauli dressings cost
$`O(Q)`$ Cliffords and no new wire. No helper is allocated. The
feature description is a proof device and is not a supplied synthesis
with free native operations.

All reductions are exact. They retain the occupied scalar, the full
unused clean block, all logical columns and reference correlations.
Any genuinely dirty work is arbitrary throughout and obeys the target's
complete return contract. No clean return or error statement is inferred
from individual scalar features.

The [focused exact fixture](../../tests/test_all_active_original_pair_closure.py)
uses the actual rank-four physical pair and an extension by one arbitrary
logical wire. It checks original-list independence and all-active support,
both alternatives of the total-axis-span test, the isotropic H_0, every
physical column and the actual inverse. In the outside-span case it
emits the clean-controlled Pauli dressing and verifies the full
$`I_c\otimes G`$ extension. It also checks the native gate ledger and
retains the original non-Boolean entries. These are mechanism witnesses,
not Boolean-target enumeration or a new hard family.

For approximate Boolean action, first apply the
[original-circuit approximation-to-exactness theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
to the original three-native-band source, before any logical restriction
or coordinate proof. Its conservative threshold $`2^{-530690}`$
therefore excludes sign degree above 2428 or XOR degree above 2427
within the stated original-list class, allowing one common scalar.
This is not an arbitrary fixed-tolerance obstruction.

The unresolved next step is the case where both actual silent lists
couple to the original middle pair. Then the Boolean endpoint lemma
does not make all retained lists commute, and Case II need not have
an empty coupled-silent sublist. Neither a bounded pair count nor a
rank-one coupling functional bounds that actual span. General original
pairs and unrestricted depth remain open; the general frame resource
frontier is unchanged.
