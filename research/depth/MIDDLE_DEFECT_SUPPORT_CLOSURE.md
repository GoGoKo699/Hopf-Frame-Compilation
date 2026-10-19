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
For actual clean sources, the original transition further excludes the
mixed-support alternative and gives a uniform bound at fixed dimension
of the original clean-active intersection; see Section 8.

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

At the level of this abstract no-pair normal form, growing overlap W
can retain a quadratic control supported wholly in W, or the one-line
mixed case of Section 5. Section 8 excludes that mixed case for actual
clean sources using their original transition. The original independent
exceptional middle list can have a nonzero Boolean cubic defect which
cancels bp_u. Equations
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

## 8. The original clean transition removes the mixed case

The [external-center continuation](EXTERNAL_CENTER_MIDDLE_CORRECTION.md)
uses the remaining source constraints below to correct all flags with
one middle Pauli. Its Boolean cubic correction remains in the middle;
it is not automatically another independent native band.

The original physical source imposes more than the abstract normal
form. Retain the actual no-middle-pair branch with both interior clean
images non-Pauli. If the original active spans have a common clean-odd
vector, the existing shared-flip theorem already gives five/four.
Here suppose that criterion fails.

In this section L_A,L_C mean the original physical **clean-active**
spans, not the logical silent endpoint spans used above. Put

```math
L_0=L_A+L_C,\qquad
R_{\mathrm{act}}=L_A^0+L_C^0,
\qquad L_A^0=L_A\cap\ker\chi,\quad
L_C^0=L_C\cap\ker\chi,
```

where $`\chi(v)=\langle v,p_A\rangle`$ is the original clean
functional. The control b initially lives on R_act. Later shared-center
enlargements may add directions outside L_0. We keep R_act distinct
from that larger computational center R throughout the argument.

### 8.1. A literal commutator determines the intrinsic quadratic symbol

There is a linear functional psi on L_0 which is zero on every A-active
axis and one on every C-active axis. Indeed, these prescriptions agree
on their intersection precisely because it contains no clean-odd
vector. Psi kills R_act and is one on the chosen h, the sum of one
A-active and one C-active axis. Choose an actual signed Hermitian Pauli
Q_psi realizing this functional; no initialized-sector assertion is used.

Retain the actual transition and its commuting Clifford aligner:

```math
F_0=E_{C,a}^\dagger E_{B,-}E_{A,a}=C_0G_0,
\qquad G_0=(-1)^{a(z)}H^{b(z)}.
```

Because Q_psi flips the H bit and fixes the old center bits, its exact
commutator with G_0 is $`(-1)^{b(z)}`$. Let s_psi be the sum of the
original signed T orientations in F_0 on the psi-active axes, including
the inverted C orientations. The literal identity is

```math
(-1)^b
 =(C_0^\dagger Q_\psi C_0Q_\psi)\,
   \omega^{-s_\psi}
   \prod_{v:\psi(v)=1} S_v^{\epsilon_v}.          \qquad(17)
```

The first factor is a Pauli unitary, with its actual scalar. All
displayed factors commute in the original active algebra. Thus their
binary symplectic displacement gives the intrinsic formula

```math
M:=N_\beta=N_C+N_{U_1}=N_A+N_{U_0},               \qquad(18)
```

where U_1 and U_0 are the spans of the original negative-middle axes
selected by psi=1 and psi=0, respectively. For an independent list X,
$`N_X(q)=\sum_{x\in X}\langle x,q\rangle x`$. The matrix M is the
polar displacement of the quadratic Boolean phase b, and
$`S_\beta:=\mathrm{im}M`$ is its intrinsic support in R_act. Signs
and T orientations are not dropped from (17); they simply do not
change the binary S-rotation displacement in (18).
The notation N_{U_j} always retains the selected original native list;
it does not permit an arbitrary basis of the same span.

The second equality follows either by using chi+psi or from the
original clean Gram identity $`N_A+N_C=N_{B,-}`$. All negative-middle
axes remain original physical labels here. Their independence gives
$`L_{B,-}=U_0\oplus U_1`$. They are not projected through a logical
clean quotient. Changing the commuting Clifford aligner changes b
only by an affine term, so it does not change beta or (18).

### 8.2. Lift the cubic-support identity faithfully to the physical space

Choose actual clean-odd axes q_A in L_A and q_C in L_C. The physical
section of the positive-clean Pauli quotient consisting of labels
commuting with q_A is faithful: the other representative differs by
p_A, which anticommutes with q_A. The old center R_act, H and every
positive-middle axis commute with q_A. Their joint logical labels
therefore lift injectively to the original physical space. Pull the
paid normalizations back through this same section before using (18).

At polar rank at least four, any nonzero linear flag ell_u from Section
5 lies in W. Its positive-middle lift p is unique and lies in the
physical positive-middle span. At either endpoint the native displacement
of h uses only silent axes: h commutes with every clean-active axis,
so those axes contribute zero to the native symplectic displacement.
It also equals the A-silent displacement
lift, since both representatives commute with q_A. For C, use the
aligned silent band $`C_0^\dagger E_{C,s}C_0`$, which preserves P_A.
The aligner fixes each original C-active axis, so both that silent lift
and the positive-middle lift commute with q_C. Their possible quotient
difference p_A does not. Hence these physical lifts agree as well.
Because the positive-middle span commutes with C_0, undoing the aligner
places p in the original C-silent span too.

Full native independence at each physical endpoint now gives

```math
p\notin L_A,\qquad p\notin L_C\quad\text{for }p\ne0.           \qquad(19)
```

In this section let W also denote its faithful physical lift from
Section 4. The original full middle list is independent, so

```math
W\cap(U_0\oplus U_1)=0.                           \qquad(20)
```

This compares only faithfully lifted positive labels with original
physical negative labels. It does not try to project the latter into
the positive-clean section.

### 8.3. A projection argument removes the one-line mixed term

Suppose the polar rank of b is at least four and a nonzero flag p
survives. The exterior result gives a physical decomposition

```math
M=M_W+p r^{\flat}+r p^{\flat},\qquad
\mathrm{im}M_W\subseteq W,\qquad
p^{\flat}(q)=\langle p,q\rangle.                  \qquad(21)
```

Here M_W is the polar map of a W-supported quadratic symbol. Project
$`N_C=M+N_{U_1}`$ to the ordinary vector-space quotient by W+U_1.
Every term vanishes except $`[r]p^{\flat}`$. If [r] were nonzero,
every vector in the kernel of N_C would pair to zero with p. Native
C-axis independence gives $`\ker N_C=L_C^\perp`$, so this would
put p in L_C, contradicting (19). Therefore r lies in W+U_1. The
same argument at A puts r in W+U_0. Equation (20) makes the
intersection of those two spaces exactly W, and hence r lies in W.

Thus M is supported wholly in W. The one-line mixed alternative of
the abstract normal form is excluded for this actual source. Moreover,
the maps M and N_{U_1} have disjoint image spaces, and each kernel is
the symplectic orthogonal complement of its image. Their sum
has kernel $`S_\beta^\perp\cap U_1^\perp`$, so (18) gives

```math
L_A=S_\beta\oplus U_0,\qquad
L_C=S_\beta\oplus U_1,
\qquad
S_\beta=L_A\cap L_C\subseteq L_{B,+}.             \qquad(22)
```

In fact S_beta lies in the exceptional positive-middle overlap W.
The equality of the active intersection follows from (20) and the
independence of U_0,U_1. In particular, a surviving nonconstant flag
forces the polar rank of b to equal the dimension of the entire
original active intersection.

This is not a contradiction. Equations (20)--(22) imply
$`W\cap L_0=S_\beta`$. By (19) the nonzero p does not lie in
S_beta and therefore lies **outside the original L_0**. Such a
direction is allowed in the enlarged shared center. Replacing that
larger R by the old R_act at this step would incorrectly claim an
unrestricted closure.

### 8.4. A fixed original active intersection gives a uniform bound

Let $`t=\dim(L_A\cap L_C)`$ for the original physical active lists.
At polar rank at least four, if the rank differs from t, (22) rules
out every nonconstant flag. Section 2 then gives sign degree at most
five. In particular this applies when the polar rank is greater than t.
Any remaining high-rank case has rank exactly t; rank two is already
bounded by finite insertion, and rank zero is absorbed as an actual
Clifford. Consequently the actual source satisfies

```math
\deg g\le
\max\left\{5,
D_{\mathrm{ins}}\bigl(0,0,0,\max\{t,2\}+1\bigr)\right\}.
                                                               \qquad(23)
```

The low-rank step uses precisely the paid affine-control absorption:
the remaining insertion has support dimension at most its polar rank
plus one. Original widths, actual native lists, literal scalars and
actual inverse words remain retained. The central cubic and real
quadratic signs do not raise (23). The corresponding literal XOR bound
is one less after checking these hypotheses in its actual output-Hadamard
sign representation. This is uniform with t fixed even if both active
ranks, the positive-middle overlap and dirty width grow.

Since polar rank is even, an odd t cannot equal a surviving high rank.
For odd t the high-rank bound is therefore five, and all control ranks
have the uniform bound $`\max\{5,D_{\mathrm{ins}}(0,0,0,3)\}`$,
even when that odd intersection dimension grows.

The abstract normal form still has the one-line mixed possibility of
Section 5. For actual sources the remaining requirement is stronger:
the full quadratic support must equal the original active intersection,
be carried by the exceptional positive-middle span, and its nonconstant
flags must point outside the original active algebra into the enlarged
center. Matching that external center transport with the full Boolean
phase remains open. Neither (22) nor the fixed-intersection bound changes
the complete-frame resource frontier or the unrestricted arbitrary-pair
one-clean question.
