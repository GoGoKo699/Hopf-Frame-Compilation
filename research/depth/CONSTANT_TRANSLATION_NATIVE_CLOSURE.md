# Constant endpoint translations force a cubic Boolean phase

[Global native matching](QUADRATIC_CONTROL_NATIVE_MATCHING.md#53-a-larger-actual-center-gives-monomial-clifford-endpoint-images) ·
[Computational derivatives](COMPUTATIONAL_PHASE_DERIVATIVES.md) ·
[Two-layer diagonal rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md)

The remaining center phase can be bounded when the actual monomial
endpoint images translate computational inputs by constants. Their
native origin forces every activated endpoint axis to be computationally
diagonal. Those actual sublists can be retained as diagonal C3 factors,
leaving a full-space two-band diagonal problem after fixed Clifford
alignment. The Boolean sign then has degree at most three, without a
bound on the number of center coordinates or endpoint axes.

## 1. Full-space hypotheses and conclusion

On Q arbitrary logical wires, including every dirty wire, let

```math
V=K_C^\dagger E_C^\dagger D E_AK_A
  =e^{i\theta}\mathrm{diag}\bigl((-1)^{g(x)}\bigr).       \qquad(1)
```

The scalar is literal. The K operators are actual Cliffords. Each E is
one actual native band with independent commuting signed Pauli axes and
its original T/T-dagger orientations. Let S be an isotropic Pauli-label
space and assume

```math
D\in\mathrm{Alg}(S),\qquad D\in C_3.                \qquad(2)
```

Thus D is a commuting Pauli C3 phase. It need not be an independent
native band; no synthesis cost is assigned to (2). In the actual
no-middle-pair clean-transition application, D is the retained product
of the quadratically controlled Pauli and the independent middle band,
with its original literal source and inverse.

For every signed Hermitian representative P_s of S, define the actual
Hermitian Clifford endpoint images

```math
Q_A(s)=K_A^\dagger E_A^\dagger P_sE_AK_A,
\qquad
Q_C(s)=K_C^\dagger E_C^\dagger P_sE_CK_C.            \qquad(3)
```

**Constant-translation hypothesis.** Every image in both families is
monomial in the full computational basis, and its basis permutation is
a fixed translation, independent of every input coordinate. Equivalently,

```math
Q_A(s)|v\rangle=\eta_{A,s}(v)|v+a_{A,s}\rangle,
\qquad
Q_C(s)|v\rangle=\eta_{C,s}(v)|v+a_{C,s}\rangle,       \qquad(4)
```

where the a vectors are constant and the phases are those of the actual
Cliffords. Quadratic work-coordinate phases are allowed. This is the
specialization $`L_s=0`$ of the enlarged-center monomial form; a constant
a_s need not vanish. Since D commutes every P_s, (1) implies
$`V^\dagger Q_C(s)V=Q_A(s)`$. Diagonality of V consequently makes
$`a_{A,s}=a_{C,s}`$.

**Constant-translation degree theorem.** Under (1)--(4), V is a diagonal
C3 operator up to its retained common scalar. In particular,

```math
\deg_{\mathbb F_2}g\le3.                           \qquad(5)
```

There is no restriction on the size of S, center dimension, endpoint
active ranks or dirty width. The hypotheses are full-space operator
identities, not properties of selected center fibers. For an actual
one-clean circuit, the cited clean-transition reduction must first
supply (1) on the entire logical register with its independent lists.

## 2. Native independence identifies actual diagonal endpoint sublists

Write M for the computational Z-label space. A constant-translation
monomial Clifford sends every computational Z to itself with a possible
sign. Its symplectic displacement therefore vanishes on M.

Fix s. The independent native-axis formula for Q_A(s) gives, at the
binary label level,

```math
(\Phi_{Q_A(s)}-I)v
 =\sum_{i:\langle a_i,s\rangle=1}
      \langle K_A^\dagger a_iK_A,v\rangle
      K_A^\dagger a_iK_A.                          \qquad(6)
```

Orientations and Pauli signs do not change this displacement formula;
the corresponding operators still retain them. The labels on the right
are independent. Because (6) vanishes for every v in M, each activated
label is orthogonal to M, hence itself lies in M. Thus every A axis
activated by at least one s in S is computationally diagonal after its
actual input endpoint Clifford. The same reasoning applies at C.

Let A_d and C_d be precisely these unions of activated axes, and let
A_0,C_0 be the complementary sublists. These are sublists of the original
native bands, not new linear combinations of axes. Define

```math
D_A=K_A^\dagger A_dK_A,\qquad
D_C=K_C^\dagger C_dK_C.
```

Both are actual diagonal C3 native sublists. A_0 and C_0 centralize S
pointwise and remain independent. Since each endpoint band is commuting,
the exact factorization is

```math
V=D_C^\dagger W D_A,
\qquad
W=K_C^\dagger C_0^\dagger D A_0K_A.                \qquad(7)
```

W is full-space diagonal because V,D_A,D_C are diagonal. It need not be
Boolean: the removal of the actual native sublists can leave complex
diagonal phases. That causes no gap in the later two-layer step, whose
full-space diagonal theorem permits such phases. No arbitrary quartic
factor has been removed, and no dependent list has been declared native.

## 3. Match the remaining Pauli families by an actual diagonal Clifford

The residual bands centralize S, so (7) gives exact signed Pauli
intertwinings

```math
W^\dagger P_C(s)W=P_A(s),\qquad
P_C(s)=K_C^\dagger P_sK_C,\qquad
P_A(s)=K_A^\dagger P_sK_A.                          \qquad(8)
```

There is an actual diagonal Clifford C matching the whole family:

```math
C^\dagger P_C(s)C=P_A(s)\quad(s\in S).             \qquad(9)
```

Here is a full-space construction. Diagonality of W makes the two Pauli
families have the same X masks. If a member has zero X mask, the two
signed Paulis agree literally. The difference between their Z labels
therefore factors through their common X-mask space U. The commutation
relations make this linear map on U the restriction of a symmetric
binary matrix. Extend it to the whole computational X space and realize
the resulting shear by S and CZ gates. The remaining sign discrepancy
is a linear character on U, because the two families have the same
literal multiplication law. A computational Z Pauli corrects that
character. The zero-mask kernel has no discrepancy, so the correction
is consistent.

Equivalently, (8) supplies a consistent signed symplectic map on the
Pauli subgroup generated by the family and all computational Z Paulis;
its extension can be chosen to fix the latter subgroup pointwise. Both
descriptions construct a global diagonal Clifford. They do not extend
an equality known only on initialized vectors or on one center fiber.
The construction uses $`O(Q^2)`$ elementary Clifford gates and no helper.

Set

```math
W_0=WC^\dagger
 =K_C^\dagger C_0^\dagger D A_0(K_AC^\dagger).
```

W_0 is still full-space diagonal and now centralizes the entire common
Pauli family P_C(S) pointwise. The input Pauli family in this displayed
word is the same signed family as the output family. Only an actual
Clifford endpoint has changed; the two remaining native lists are
unchanged.

## 4. Normalize the common algebra and apply the complex diagonal theorem

The common Pauli family can be made computational without destroying
this residual diagonality. First use CNOTs to put its X-mask space in
coordinate form. A diagonal S/CZ/Z Clifford removes the Z decorations
and signs on commuting lifts of those X masks. The zero-mask kernel
is then placed in coordinate Z form using CNOTs on the remaining
coordinates. These are actual fixed coordinate operations.

The resulting W_0 is diagonal and commutes each of the selected X
operators, so it acts as identity on those particular tensor factors.
Hadamards on exactly those factors leave it unchanged and convert the
common family to a computational Z subgroup. This argument uses no
initialized work and holds for the complex diagonal W_0 as well as for
a Boolean sign. It is the identity action, not a chosen fiber extension,
that permits these Hadamards.

Choose matching fixed Clifford coordinates on the middle space S.
Both resulting endpoint Cliffords preserve the now-computational
subgroup pointwise, and both remaining native bands centralize it.
The transformed D is supported on that subgroup. Hence it moves
literally to the input as a diagonal C3 factor D_in. In these coordinates
write

```math
W_0'=W_2D_{\mathrm{in}},
\qquad
W_2=(K_C')^\dagger(C_0')^\dagger A_0'K_A'.         \qquad(10)
```

W_0' and D_in are diagonal, so W_2 is diagonal on the full register.
Its two phase lists are retained independent native sublists with actual
Clifford interlayers and endpoints. The
[full-space two-layer diagonal theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#1-model-and-full-space-diagonal-theorem)
makes W_2 diagonal C3, including its complex phases and common scalar.
The product with D_in is therefore diagonal C3.

Clifford conjugation preserves C3 membership, so returning from the fixed
coordinates makes W_0 a diagonal C3 operator. Restoring the actual
diagonal Clifford C gives the same conclusion for W. Finally all three
operators in $`V=D_C^\dagger W D_A`$ are diagonal C3. Diagonal C3
phases form a group, proving (5).

This is a degree factorization with actual sublists and coordinate
words. Removing D_in is not a claim of a new implementation of the
original circuit at unchanged native depth. Every literal scalar,
orientation and actual inverse remains attached to its source word.
The additional aligners and coordinate changes have $`O(Q^2)`$
elementary Clifford cost, with the same bound charged separately for
their reversed inverse words. No clean or dirty helper is added, and
no arbitrary logical input is initialized or discarded.

## 5. Application, remaining gap, and the continuing frontier

For a literal XOR query, apply the theorem to its actual output-Hadamard
sign representation. If its predicate f is nonzero, the sign is
$`(-1)^{y f(x)}`$, so (5) implies $`\deg f\le2`$. The
constant-translation hypotheses must be checked in that representation;
no invariance under the output basis change is presumed. The retained
exact-rounding theorem may first be applied at its stated accuracy floor;
it does not remove the hypotheses of this exact result.

The one-clean reduction retains a known central cubic sign and an
actual real quadratic normalization. Restoring those signs does not
raise the cubic bound. The associated affine input changes preserve
Boolean degree, and the Hadamards were used only on exact identity modes.
Thus this is a genuine growing-rank subcase of the actual no-pair
problem when the resulting endpoint family has constant translations.
It does not require those translations to be zero or demand a Clifford
decoder for the original monomial algebra.
The original Q families can remain non-Pauli until their actual diagonal
native sublists are peeled; no fixed Clifford decoder is assumed.

When the whole middle support is already a genuine shared computational
center, the theorem includes the simpler argument: D moves to the input
directly, and the remaining complex diagonal two-band word is C3. The
more general proof above first obtains that common center by retaining
the identifiable diagonal endpoint sublists and their phases.

### 5.1. An empty native middle list has a uniform bound

There is a further closure in the actual
[quadratic-control normal form](QUADRATIC_CONTROL_NATIVE_MATCHING.md#1-actual-input-controls-and-the-full-space-word)

```math
V=K_C^\dagger E_C^\dagger H^{b(z)}E_BE_AK_A
  =e^{i\theta}D_g.
```

Retain all of that form's hypotheses: the three lists are globally
independent native bands; R is a genuine shared computational center;
the signed Hermitian H commutes with R and each middle axis; and b is
Boolean quadratic in the actual center inputs. Suppose the B list is
empty. Then, uniformly in width and control rank,

```math
\deg g\le D_{\mathrm{ins}}(0,0,0,3).              \qquad(11)
```

Here D_ins is the explicit
[finite-insertion bound](FINITE_INSERTION_THREE_LAYER_FEATURES.md#1-the-full-space-expression-and-explicit-bound).
The sharper cubic bound holds when the polar rank of b is zero or at
least four; only rank two uses (11).

For the high-rank assertion, the simultaneous native matching identities
on $`S=R+\langle h\rangle`$ already place $`P_u`$ in the original
computational center R when B is empty. They give

```math
P_u=H(J_uHJ_u^\dagger)=(-1)^{p_u(z)},\qquad
C_uJ_u^\dagger=P_u^{b(z)}=(-1)^{b(z)p_u(z)},       \qquad(12)
```

where p_u is an affine center parity with its literal sign, and the
right side is an actual diagonal Clifford. The native middle defect
is exactly I because the list is empty. Thus the Boolean polynomial
$`b p_u`$ has no cubic part. If the linear part ell of p_u is nonzero,
choose affine center coordinates with ell as the first coordinate.
Vanishing of the cubic part forces the quadratic part of b to have
that coordinate in every quadratic monomial. Equivalently it is
$`\ell\wedge m`$ for a linear m, and its polar rank is at most two.
At polar rank at least four, every p_u is therefore constant.

Consequently all computational-probe actions on S have zero binary
displacement: they already fix R, and (12) now fixes the label h as
well. Each endpoint image Q_A(s),Q_C(s) sends every computational Z
to itself up to sign. Its permutation is a constant translation, so
the theorem above proves $`\deg g\le3`$.

At rank zero, $`H^{b(z)}`$ is an actual Clifford; absorbing it leaves
only the two independent endpoint bands, whose diagonal theorem gives
the same cubic bound. At rank two, the paid affine-control absorption
leaves an insertion supported on H and two actual center parities.
The finite-insertion theorem applies with relation ranks zero and
support dimension at most three, proving (11). These are actual fixed
Clifford changes with their retained phases and inverse words, not
separate choices on center fibers. Restoring the known central cubic
and real quadratic signs preserves the displayed bounds.

There is a useful sufficient condition for the empty-list hypothesis.
Let $`\mathsf P_A,\mathsf P_B,\mathsf P_C`$ be the literal products
of all signed axes in the three native lists. For a genuine quadratic
control, the
[whole-word Boolean Galois identity](NATIVE_CENTER_SATURATION_LIMITS.md#4-whole-word-galois-products-supply-at-most-two-further-directions)
is $`\mathsf P_A=\pm\mathsf P_B\mathsf P_C`$. If
$`\mathsf P_A=\pm\mathsf P_C`$, it forces the binary label of
$`\mathsf P_B`$ to vanish. Independence makes a nonempty B list
impossible, so the preceding result applies. For affine b, do not invoke
that genuine-control identity: direct Clifford absorption instead
leaves the existing three-native-band bound five. Thus matching endpoint
total products give the uniform sign bound
$`\max\{5,3,D_{\mathrm{ins}}(0,0,0,3)\}=D_{\mathrm{ins}}(0,0,0,3)`$
for every quadratic b. In the actual XOR sign representation the
corresponding predicate bound is one less. This does not identify the
individual endpoint axes or remove the distinct-total, nonempty-middle
case.

### 5.2. The remaining center-dependent case

The unresolved case has nonzero center-dependent translation matrices
in $`y\mapsto y+a_s+L_sz`$. Their nontrivial native endpoint activation
cannot be peeled as computationally diagonal sublists by (6). Although
the enlarged center is shared and the middle span is invariant, the
native middle matching defect can still carry a central cubic cocycle.
Multiplying the middle factor by a central phase does not remove that
cocycle because the whole probe family fixes the center pointwise.
A new argument must match those actual center-dependent translations
and phases while keeping all original global native lists. Separate
fiber decoders, arbitrary quartic removal, and replacement of independent
Pauli axes by non-Pauli Clifford involutions do not supply that step.

The general destination remains a uniform bound for growing ranks, or
a complete actual native family of growing Boolean degree. This scoped
closure changes neither the complete-frame resource frontier nor the
unrestricted one-clean question. It constrains a previously unpaid
center-phase step using the actual endpoint bands, rather than merely
bounding another representation parameter.
