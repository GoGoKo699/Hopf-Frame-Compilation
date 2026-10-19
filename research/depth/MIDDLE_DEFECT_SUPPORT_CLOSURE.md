# The native middle defect confines the remaining quadratic control

[Quadratic-control matching](QUADRATIC_CONTROL_NATIVE_MATCHING.md) ·
[Constant translations](CONSTANT_TRANSLATION_NATIVE_CLOSURE.md) ·
[Native center constraints](NATIVE_CENTER_SATURATION_LIMITS.md) ·
[Finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md)

The quadratic-controlled insertion need not have constant translations
on the entire middle algebra. Constant translations on its own Pauli
support already give sign degree at most five. The simultaneous native
middle defects force this weaker condition unless the quadratic control
has a specific concentration in the shared center contained in the
exceptional middle span. In particular, a bounded dimension of that
intersection gives a uniform degree bound even when endpoint ranks,
middle width and the number of other center coordinates grow.

These are full-space degree results. They retain the original independent
native lists and literal phases. The intermediate factorizations do not
give a new same-depth synthesis of a controlled gate.

## 1. Actual normal form and simultaneous phase identities

Use the retained full-space logical normal form

```math
V=K_C^\dagger E_C^\dagger F E_BE_AK_A
  =e^{i\theta}D_g,\qquad F=H^{b(z)}.                \qquad(1)
```

Every input, including each dirty input, is arbitrary. Each E is an
actual independent commuting signed Pauli T/T-dagger list; the three
lists need not have an independent union. The K's are actual Cliffords.
The signed Hermitian Pauli H commutes with every middle axis and with
the genuine shared computational center R. Every band and both
normalized endpoints fix R pointwise. The Boolean b has degree at most
two in the actual center bits z. The source of F, the known central cubic
sign and the real quadratic normalization keep their original accounting.

For a genuine quadratic b, use the already-paid simultaneous center
enlargement of
[global native matching](QUADRATIC_CONTROL_NATIVE_MATCHING.md#53-a-larger-actual-center-gives-monomial-clifford-endpoint-images).
It retains the old control coordinates while adjoining actual shared
parities, so it does not alter the polar rank of b. Write L_B for the
span of the original independent middle labels and set

```math
S=R+L_B+\langle h\rangle,\qquad
J_u=E_AK_AZ(u)K_A^\dagger E_A^\dagger,\qquad
C_u=E_CK_CZ(u)K_C^\dagger E_C^\dagger.
```

The two families are exact commuting representations of the full
computational Z group. They preserve S, fix R pointwise, agree literally
on S, and preserve L_B. Their binary displacement on S lies in R.
Consequently there are affine center bits, including their actual signs,
such that

```math
J_uB_iJ_u^\dagger=(-1)^{\ell_{u,i}(z)}B_i,
\qquad
J_uHJ_u^\dagger=(-1)^{p_u(z)}H.                    \qquad(2)
```

The independent signed B axes have original orientations
$`\epsilon_i\in\{1,-1\}`$. Their literal defect and its actual
Clifford correction satisfy

```math
\begin{aligned}
D_u&=E_BJ_uE_B^\dagger J_u^\dagger
 =\omega^{-\sum_i\epsilon_i\ell_{u,i}(z)}
   \prod_i S_{B_i}^{\epsilon_i\ell_{u,i}(z)},\\
K_u&=C_uJ_u^\dagger=D_u(-1)^{b(z)p_u(z)}.
\end{aligned}                                                   \qquad(3)
```

Here $`\omega=e^{i\pi/4}`$ and $`S_P=T_P^2`$. Equation (3)
keeps the scalar from $`T_{-P}=\omega T_P^\dagger`$. All its
factors commute in the isotropic algebra of S. In particular, K_u is
an actual diagonal Clifford in any fixed coordinates diagonalizing S.

For clarity, these are simultaneous identities, not separate choices
for individual probes. In one fixed global coordinate system on S,
let beta_i be the affine eigenbit of the signed B_i and put
$`\phi_B=\sum_i\epsilon_i\beta_i`$. If tau_u is the affine action
of J_u on those coordinates, then

```math
\phi_B-\phi_B\circ\tau_u+4bp_u=\kappa_u\pmod8.    \qquad(4)
```

The tau actions fix z; the affine flags in (2) add under probe addition.
The literal diagonal-Clifford phase kappa obeys
$`\kappa_{u+v}=\kappa_u+\kappa_v\circ\tau_u\pmod8`$.
Only z are already actual computational input coordinates. Equation
(4) does not identify the other middle eigenbits with target inputs.

## 2. Constant translations on the insertion support suffice

The following extension uses the complex diagonal version of the
[three-native-layer theorem](THREE_LAYER_NATIVE_STRUCTURE.md#1-exact-statement-and-conventions).

**Partial-support closure.** In a full-space Boolean word (1), more
generally let F be a commuting C3 phase supported in an isotropic Pauli
space T, and suppose every actual middle axis centralizes T. If both
endpoint image families

```math
K_A^\dagger E_A^\dagger P_tE_AK_A,\qquad
K_C^\dagger E_C^\dagger P_tE_CK_C\quad(t\in T)
```

have constant computational translations, then

```math
\deg g\le5.                                      \qquad(5)
```

Their translations may be nonzero, and their diagonal Clifford phases
may be non-Pauli. No corresponding hypothesis is imposed on the other
middle-axis probes.

Here is the paid reduction. Apply the native displacement argument of
[constant-translation closure](CONSTANT_TRANSLATION_NATIVE_CLOSURE.md#2-native-independence-identifies-actual-diagonal-endpoint-sublists)
only to T. Each actual endpoint axis activated by T is computationally
diagonal after its endpoint Clifford: the symplectic displacement
vanishes on computational Z, and independence prevents cancellation
between its activated axis vectors. Retain exactly these original
sublists as diagonal C3 factors D_A,D_C. The residual full word is

```math
V=D_C^\dagger W D_A,\qquad
W=K_C^\dagger C_0^\dagger F E_BA_0K_A,             \qquad(6)
```

where A_0,C_0 remain independent and centralize T pointwise. W is
diagonal, possibly with complex phases. The middle list also centralizes
T. Thus W intertwines the actual input and output Pauli families on T.
Their signed zero-X-mask kernel agrees literally. As in the cited proof,
extend their symmetric Z-label difference by a diagonal S/CZ Clifford,
and correct its residual sign character by a Z Pauli. This constructs
one actual diagonal Clifford C aligning the whole family. Then
$`W_0=WC^\dagger`$ is diagonal and centralizes one common Pauli family.

Use actual CNOTs and diagonal Clifford decorations to put that family
in X/Z coordinate form. Let W_pre denote W_0 in these coordinates.
It is diagonal and commutes every selected X, so it is exactly identity
on those particular tensor modes I. Hadamards on those modes make the
whole family computational without changing W_pre:

```math
W_{\mathrm{pre}}=H_IW_{\mathrm{pre}}H_I=W_0'.       \qquad(7)
```

Matching fixed middle coordinates make F a diagonal C3 factor F_in
which moves literally to the input. What remains is a full-space
complex diagonal word with the three retained native sublists:

```math
W_0'=W_3'F_{\mathrm{in}}',\qquad
W_3'=(K_C')^\dagger(C_0')^\dagger E_B'A_0'K_A'.    \qquad(8)
```

The coordinate undo matters. Equation (7) allows the same factorization
(8) for W_pre without conjugating either factor by H_I. Neither factor
is asserted to be identity on I separately. Undo only the preceding
affine coordinate permutation and diagonal decorations on W_3' and
F_in'. The first remains an actual full-space diagonal three-native-band
word in the original computational coordinates; the second remains
diagonal C3. The three-layer theorem gives an ordinary mod-eight phase
polynomial of degree at most five for the first word in those coordinates.
The second and the restored D_A,D_C,C have phase degree at most three.
Their exact product therefore has degree at most five. For the original
Boolean sign this proves (5).

This does not assert affine invariance for an arbitrary degree-five
phase polynomial or conjugate the separate factors through an unpaid
Hadamard. The theorem is applied to the actual diagonal native word
after its affine coordinate undo. All fixed coordinate and alignment
words use $`O(Q^2)`$ elementary Cliffords on Q logical wires, with
their separately reversed actual inverse words and literal scalars.
No helper, conditional decoder or extra native-band synthesis is used.

For (1), it is enough that every p_u in (2) be constant. Indeed J_u and
C_u then have zero binary displacement on $`T=R+\langle h\rangle`$.
The corresponding endpoint images fix every computational Z label,
possibly changing its sign, and hence have constant translations.
This proves (5) without requiring constant translations on all of S.
The earlier stronger hypothesis on all of S still gives its sharper
sign bound three; it is not replaced by (5).

## 3. Remove only true native permutation orbits from the cubic defect

The original B basis supplies an intrinsic useful pruning. Let B_perm
be the set of actual B axes whose whole orbit under the binary J action
stays in the actual B-axis set. It is a union of complete group orbits.
Let B_exc be its complementary actual sublist and put

```math
L_{\mathrm{exc}}=\mathrm{span}(B_{\mathrm{exc}}).
```

The native Gram identity from global matching makes the binary J action
orthogonal in the independent B basis. The coordinate span of B_perm
is invariant, so its coordinate complement L_exc is invariant too.
This is a partition of the original native axes, not an orthogonal
change to a more convenient independent basis.

Each B_perm orbit has size one or two. If an elementary abelian orbit
had at least four members, two independent orbit translations u,v would
make $`(\Phi_{J_u}-I)(\Phi_{J_v}-I)b`$ the sum of four distinct B
axes. Independence makes that sum nonzero. The common isotropic A
encoder, however, gives $`N_uN_v=0`$ for these displacements, a
contradiction.

Split the literal phase list as $`E_B=E_{\mathrm{perm}}E_{\mathrm{exc}}`$.
The signed-permutation defect

```math
D_{\mathrm{perm},u}
 =E_{\mathrm{perm}}J_uE_{\mathrm{perm}}^\dagger J_u^\dagger
```

is an actual Clifford: after the signed axis permutation, each odd
T exponent changes by an even amount. The scalar caused by a negative
signed axis remains attached. The residual defect is supported on the
invariant span L_exc and satisfies

```math
\begin{aligned}
D_{\mathrm{exc},u}
 &=E_{\mathrm{exc}}J_uE_{\mathrm{exc}}^\dagger J_u^\dagger,\\
\widetilde K_u
 &=D_{\mathrm{perm},u}^\dagger K_u
   =D_{\mathrm{exc},u}(-1)^{bp_u}\quad\text{is Clifford}.
\end{aligned}                                                   \qquad(9)
```

All displayed factors commute. In particular, the cubic defect can be
localized to the actual exceptional coordinate span. Keeping all B axes
instead also gives the valid weaker conclusions below with L_B in
place of L_exc.

## 4. Normalize the defect to an exact Boolean cubic

Define the shared overlap

```math
W=R\cap L_{\mathrm{exc}}.
```

Do not take the cubic coefficients of an unnormalized weighted phase as
an intrinsic tensor. A literal Boolean normalization supplies the needed
coordinate-independent statement.

Because the Boolean factor in (9) squares to I,
$`D_{\mathrm{exc},u}^2=\widetilde K_u^2`$. A Clifford inside an
isotropic Pauli algebra is a diagonal Clifford in fixed coordinates,
so its square is a Pauli unitary, including its scalar. Since the left
side is supported on L_exc, write it as $`\lambda_u P_u'`$ with
P_u' a signed Hermitian Pauli there. The literal commuting signed
T/T-dagger product D_exc has eighth-root eigenvalues, so lambda_u is a
fourth root and its square root is an allowed eighth-root Clifford
scalar. Choose the actual diagonal Clifford

```math
L_u=\sqrt{\lambda_u}\,S_{P_u'},\qquad
L_u^2=D_{\mathrm{exc},u}^2,\qquad
U_u=L_u^\dagger D_{\mathrm{exc},u}.                \qquad(10)
```

The specified square root is retained as a native Clifford scalar;
identity-axis cases use the same literal equation. Then U_u is an exact
Boolean C3 phase supported on L_exc, and

```math
L_u^\dagger\widetilde K_u=U_u(-1)^{bp_u}
```

is an exact Boolean Clifford. This is an operator identity for analyzing
the existing defect, not a free gauge inserted into V.

Let beta be the quadratic homogeneous symbol of b in the exterior
algebra of the center parity space R, and let ell_u be the linear part
of p_u. The polar rank of beta is the polar rank of b. Boolean cubic
symbols are invariant under affine changes of coordinates; this is the
degree-three associated graded part of the Boolean polynomial algebra.
Equation (10) and the Boolean Clifford identity give, simultaneously,

```math
\beta\wedge\ell_u\in\Lambda^3W
       \quad\text{for every computational probe }u.              \qquad(11)
```

To see the support intersection explicitly, use one fixed basis of
$`R+L_{\mathrm{exc}}`$ adapted to their intersection. The Boolean
cubic of U_u uses only L_exc, while that of bp_u uses only R; their
difference is Boolean quadratic. The common cubic therefore belongs to
$`\Lambda^3R\cap\Lambda^3L_{\mathrm{exc}}=\Lambda^3W`$.
No eigenvalue fiber has been selected, and none of these middle
coordinates has been silently identified with a target input.

## 5. The exterior constraint gives growing-rank closures

The following elementary consequence makes (11) useful. For a nonzero
linear ell with $`\beta\wedge\ell\in\Lambda^3W`$:

* If ell does not lie in W, then $`\beta=\ell\wedge a`$ for some a,
  and the polar rank of beta is at most two.
* If ell lies in W, then
  $`\beta\in\Lambda^2W+\ell\wedge R`$.

For the first claim, complete W and ell to a basis. Every term in
$`\beta\wedge\ell`$ has the ell coordinate, so it can belong to
$`\Lambda^3W`$ only by vanishing. The kernel of wedge multiplication
by nonzero ell on quadratic symbols is exactly $`\ell\wedge R`$.
For the second claim choose a complement of W. Terms with two complement
coordinates must vanish; for each complement coordinate, its W-valued
coefficient must be proportional to ell. These are precisely the two
asserted forms.

Assume now that the polar rank of b is at least four. Every nonzero
ell_u lies in W. Choose a linear complement $`R=W\oplus R'`$ and
decompose beta into its W-only, mixed and R'-only parts. This gives
the following exact alternatives.

* If the projection of beta to $`\Lambda^2(R/W)`$ is nonzero, every
  ell_u vanishes, so the partial-support theorem gives $`\deg g\le5`$.
* If that projection vanishes but the mixed coefficient map has rank
  at least two, the same conclusion holds. Its W-valued image would
  otherwise have to lie in the single line spanned by each nonzero ell_u.
* If the mixed rank is exactly one, every nonzero ell_u lies on that
  one fixed line in W. This is a reduction of the simultaneous defect,
  not yet a degree bound for that remaining case.
* If the mixed rank is zero as well, the quadratic symbol is entirely
  in $`\Lambda^2W`$. Equation (11) leaves the remaining native cubic
  cancellation inside W open.

The quotient projection is intrinsic. When it vanishes, the mixed
coefficient map and its rank do not depend on the chosen complement.
Affine terms of b and the literal constant parts of p_u are retained;
they do not change these homogeneous cubic implications.

One coordinate-free sufficient condition for the first two closures is

```math
\min_{\gamma\in\Lambda^2W}
    \mathrm{rank}_{\mathrm{polar}}(\beta+\gamma)\ge4.
                                                               \qquad(12)
```

Indeed any nonzero ell_u would express beta as a W-supported quadratic
plus one decomposable quadratic, contradicting (12). Condition (12) is
sufficient, not necessary: the nonzero quotient projection above can
also close a case whose minimum residual polar rank is only two.

## 6. A bounded shared overlap closes every control rank

Let $`d=\dim W`$. If a nonzero ell_u exists, Section 5 gives

```math
\mathrm{rank}_{\mathrm{polar}}b\le d+2.             \qquad(13)
```

Thus polar rank greater than d+2 forces all p_u constant and gives the
sign/XOR bounds five/four. At nonzero polar rank at most d+2,
the retained
[affine-control absorption](QUADRATIC_CONTROL_NATIVE_MATCHING.md#2-fixed-polar-rank-gives-a-uniform-degree-bound)
leaves an insertion supported on at most d+3 actual Pauli directions.
It uses an actual Clifford, conjugates the output native list without
changing its independence, and retains every scalar and inverse word.
At polar rank zero the entire insertion is Clifford and the ordinary
three-native-band bound five applies.

Consequently a fixed d gives the uniform bound

```math
\deg g\le
B(d):=\max\{5,D_{\mathrm{ins}}(0,0,0,d+3)\},       \qquad(14)
```

with the explicit D_ins from
[finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md#1-the-full-space-expression-and-explicit-bound).
No endpoint clean-active rank, number of original middle axes, total
center dimension or dirty width appears in B(d). Using the larger
$`R\cap L_B`$ gives a valid, possibly weaker value of d. If every
middle axis belongs to a true native permutation orbit, L_exc is zero;
high control rank then gives five/four, and low rank uses the already
bounded insertion. The empty-middle case retains its previously proved
sharper high-rank three/two bound.

For literal XOR queries, first form the actual output-Hadamard sign
representation and check the hypotheses there. A nonzero predicate has
degree one less than its sign bound. The known central cubic and real
quadratic signs in the clean reduction do not raise (5) or (14).
Every statement is on the full logical register with arbitrary dirty
states and references; the actual inverse and literal source of F
remain charged. These reductions are degree accounting, not new
same-native-depth implementations.

## 7. What remains after the reduction

The whole-word Boolean Galois relation and its genuine parity centers
remain in force. The proof does not replace them by probe-by-probe
phase choices. Nor does it replace a non-permutation orthogonal action
on the exceptional B coordinates by a native axis permutation.

The remaining case within this no-pair normal form has growing overlap
W and can retain a quadratic control supported wholly in W, or the
one-line mixed case of Section 5. There the original independent exceptional middle list
can have a nonzero Boolean cubic defect which cancels bp_u. Equations
(3)--(4) still constrain the whole family, but setting that defect to
zero, rebasing the native list, or deleting the remaining central phase
is unsupported. A decoder chosen separately for each center eigenvalue
does not supply the missing global native matching.

This closes new growing-rank subcases while leaving the complete-frame
resource frontier and the unrestricted one-clean degree question
unchanged. The next question is the simultaneous native cancellation
within the surviving shared overlap, with the full Boolean target and
all literal Galois and endpoint constraints retained.

The [bounded exact fixtures](../../tests/test_middle_defect_support_closure.py)
check the exterior implication for all dimension-five alternating forms,
the signed native permutation defect with its actual inverse, and a
non-permutation native defect whose Clifford normalization leaves a
Boolean cubic. The uniform conclusions are the arguments above, not
an extrapolation from those finite checks.

### 7.1. Native independence and the total-product relation do not remove the overlap

The following exact skeleton separates the remaining phase constraint
from a tempting purely linear inference. It satisfies native list
independence, the clean Gram identity and the total-product Galois
relation, but it is not a Boolean query source.

Use ten physical wires with independent computational Z labels
a,c,p,z_1,...,z_6,d. The same letters denote their eigenbits in phase
formulas. Put

```math
w=\sum_{j=1}^6z_j,\qquad h=a+c,\qquad
v=p+z_1+z_2.
```

All displayed phase orientations are positive T. The clean-active and
middle lists are

```math
\begin{aligned}
A_a&=\{a\}\cup\{a+z_j:1\le j\le6\},\\
C_a&=\{c\}\cup\{c+p+z_j:1\le j\le6\},\\
B_-&=\{a+w,c+p+w,c,c+p\},\\
B_+&=\{h+z_1,h+z_3,h+z_4,h+p+z_2+z_3+z_4\}.
\end{aligned}                                                   \qquad(15)
```

Take the input clean Pauli to be $`P_A=X_aX_c`$. The full A list
adds the silent axes $`X_c,X_cZ_v,Z_d`$. In the aligned output frame
the silent C axes are $`X_a,X_aZ_v,Z_dZ_v`$. The actual original
C-silent axes are their conjugates by

```math
C_0=Z_cZ_pS_{Z_w}\,CZ_{a,w},
```

namely $`X_aZ_w,X_aZ_wZ_v,Z_dZ_v`$. Here CZ_{a,w} is the actual
Clifford controlled by the indicated two commuting parity coordinates.
The full endpoint lists each have ten independent commuting labels;
the combined B list has eight. Choose the actual endpoint words
$`K_A=CX_{a,c}H_a`$ and $`K_C=C_0K_A`$, with a the input clean wire
and c the input logical t wire. K_A sends the input clean Z to P_A
and the input logical Z_t to Z_h while fixing p,z_1,...,z_6,d. K_C
sends the input clean Z to $`P_C=C_0P_AC_0^\dagger`$. The full
source is $`K_C^\dagger E_C^\dagger E_BE_AK_A`$. Its literal active
transition is

```math
E_{C_a}^\dagger E_{B_-}E_{A_a}
 =C_0(-1)^{(a+c+p)e_2(z)},\qquad
e_2(z)=\sum_{i\lt j}z_iz_j,                         \qquad(16)
```

and it sends P_A to $`-X_aX_cZ_w`$. Its binary clean Gram identity
is $`N_{A_a}+N_{C_a}=N_{B_-}`$. The even-active center
$`\langle p,z_1,\ldots,z_6\rangle`$ is enlarged by the shared
spectator parity to
$`R=\langle p,z_1,\ldots,z_6,d\rangle`$. The vector v lies in R
but in none of $`L_{A_a},L_{C_a},L_{B_-}`$. Thus independence does
not force this candidate silent displacement into the negative-middle
span. The fresh d padding is relevant: the retained logical endpoint
total products are respectively $`Z_vZ_d`$ and $`Z_d`$, while the
positive-middle product is $`Z_v`$. They satisfy the literal Boolean
total-product relation, not merely its scalar-free label version.

Nevertheless, the phase equation fails. The original computational
probe with middle label h shifts every positive B axis by v at the
binary level. None of those images is another axis in B_+, so this is
a genuine exceptional span. The corresponding linear part of p_u is
v. The actual residual quadratic control is $`b=e_2(z)`$, of polar
rank six; the extra term $`p e_2(z)`$ in (16) is a retained central
cubic sign. Its cubic defect $`\beta\wedge v`$ has coefficient one
on $`p z_5z_6`$, whereas the positive B span and all its images use
neither z_5 nor z_6. Equation (11) therefore cannot hold.

The failure is also visible with literal phases. On the occupied
fiber p=1,z=0,d=0, the aligned endpoint core bands are both
$`\omega I`$, the positive-middle band is $`\omega S_{Z_t}`$,
and the output padding contributes $`\omega^{-1}`$ after inversion.
Their product is exactly $`S_{Z_t}`$, not a common scalar times a
Boolean sign. This fiber is a check of the fixed full source, not a
replacement circuit. The example refutes only the proposed inference
from Gram identities, independence and total products; it does not
counterexample the full Boolean theorem or supply growing Boolean
degree. The focused fixture checks its labels and literal parity
phases without large matrix simulation.
