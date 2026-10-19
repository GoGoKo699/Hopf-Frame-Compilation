# One middle Pauli absorbs the external center displacement

[Actual-source support](MIDDLE_DEFECT_SUPPORT_CLOSURE.md#8-the-original-clean-transition-removes-the-mixed-case) ·
[Global native matching](QUADRATIC_CONTROL_NATIVE_MATCHING.md#5-a-genuine-quadratic-control-forces-a-common-invariant-middle-span) ·
[Constant-translation closure](CONSTANT_TRANSLATION_NATIVE_CLOSURE.md)

The remaining actual no-pair branch has a uniform degree bound when
the total middle displacement space has fixed dimension. The bound
allows the original clean-active intersection, quadratic-control rank
and dirty width to grow. Its proof uses one fixed product of the
original positive-middle axes to describe every nonconstant displacement,
then normalizes the entire fixed middle Pauli family and retains an
actual two-native-band residual on a bounded quantum core.

This distinction is substantive. An explicit transvection model below
shows that the corrected middle need not admit even a one-native-band
implementation with arbitrary Clifford endpoints on the same full
register. The model is a barrier to that repair step, not a Boolean
query counterexample. The degree proof in Section 7 does not use that
repair. The actual clean-returning example in Section 6 likewise retains
a complex phase and does not satisfy the Boolean premise. No
complete-frame resource frontier or unrestricted one-clean degree bound
changes here.

The [native-orbit continuation](GROWING_DISPLACEMENT_NATIVE_ORBITS.md)
subsequently closes growing total displacement in this actual no-middle-pair
branch. It retains the growing permutation part as a third native band
rather than applying Section 7 to a smaller dimension without paying that
part. The original physical middle-pair and higher-layer questions remain
open.

## 1. Retain the actual source and its whole probe family

Use the full-space logical word and literal conventions of
[middle-defect support](MIDDLE_DEFECT_SUPPORT_CLOSURE.md):

```math
V=K_C^\dagger E_C^\dagger H^{b(z)}E_BE_AK_A
  =e^{i\theta}D_g.                                      \qquad(1)
```

Each E is an independent commuting signed Pauli T/T-dagger list on
the whole register. Dirty inputs and references are arbitrary. The
quadratic b uses genuine shared computational center coordinates, and
H commutes with that center and every middle axis. All previously
extracted central cubic signs and actual Clifford changes keep their
original accounting.

This note concerns the remaining actual no-middle-pair branch after
Section 8 of the cited proof: both interior clean images are non-Pauli,
there is no already-closing shared clean flip, and the polar rank
$`t=\mathrm{rank}\beta`$ is at least four. The intrinsic support
S_beta equals the original physical clean-active intersection and is
carried by the positive-middle span. A nonconstant displacement may
point outside the original active algebra into the enlarged genuine
center; that possibility is not identified with the old active center.

Index the independent signed positive-middle axes by B_1,...,B_m.
Their labels give a fixed coordinate space $`\mathbb F_2^m`$ with
its ordinary dot product. This coordinate metric is different from
the physical Pauli symplectic pairing. No new native basis is chosen.
The actual computational probe families are

```math
J_u=E_AK_AZ(u)K_A^\dagger E_A^\dagger,
\qquad
C_u=E_CK_CZ(u)K_C^\dagger E_C^\dagger.
```

Their common action on the middle labels is an orthogonal involution
$`O_u=I+N_u`$. Thus N_u is symmetric in these fixed middle
coordinates and $`N_u^2=0`$. All the probe actions commute. Put

```math
U=\sum_u\mathrm{im}N_u,\qquad K=U^\perp.
```

Every displacement lifts to the original A-silent span. The support
S_beta is in the original A-active span, and the whole physical
endpoint list is independent. Consequently

```math
U\cap S_\beta=0.                                  \qquad(2)
```

Here the physical-to-logical section is the faithful one used in the
actual-source theorem. In particular, (2) is not obtained by projecting
original negative-middle labels through a clean quotient.

Write p_u for the label of the linear part of the signed affine flag
in $`J_uHJ_u^\dagger=(-1)^{p_u(z)}H`$; when the distinction matters,
call the full affine bit f_u(z). Literal Boolean normalization of the
native defect gives the tensor identity

```math
c(O_u)=\beta\wedge p_u,                            \qquad(3)
```

where, in the fixed native coordinates,

```math
c(O)(x,y,z)
 =\sum_i(Ox)_i(Oy)_i(Oz)_i+\sum_i x_i y_i z_i.
```

This is the polarization of the normalized Boolean cubic, not a raw
weighted-phase coefficient. The normalization retains its actual
Clifford square root and scalar as in the cited proof. Constants and
affine parts of b do not enter beta; they reappear in the exact operator
identities below.

## 2. One fixed native-coordinate product corrects every probe

Restriction of beta to K retains its full rank t. Indeed, its radical
is $`S_\beta^\perp`$, and (2) implies that K surjects onto the
quotient by that radical. The same argument applies to
$`\ker N_u=(\mathrm{im}N_u)^\perp`$ for each u.

First restrict (3) to three vectors in ker N_u. Its left side vanishes.
If p_u had a nonzero restriction there, the vanishing of
$`\beta\wedge p_u`$ would force the restricted alternating form
to have rank at most two. Its rank is t, so

```math
p_u\in\mathrm{im}N_u.                             \qquad(4)
```

Choose one pair $`x,y\in K`$ with $`\beta(x,y)=1`$ and use
this same pair for every probe. Let x*y denote coordinatewise product
in the original B list and set $`v=x*y`$. For every z, both x and y
are fixed by O_u. Equations (3)--(4) give

```math
\begin{aligned}
c(O_u)(x,y,z)
 &=\langle x*y,N_uz\rangle
  =\langle N_uv,z\rangle,\\
(\beta\wedge p_u)(x,y,z)&=\langle p_u,z\rangle.
\end{aligned}
```

Therefore the one fixed vector v satisfies

```math
N_uv=p_u\qquad\text{for every }u.                  \qquad(5)
```

Let P_v be the literal signed product of the original B_i selected by
v. It is a Hermitian Pauli because those actual operators commute.
Define

```math
H'=HP_v.
```

The linear flags of H and P_v cancel by (5). Hence there is a constant
sign character epsilon_u, with no dependence on any input coordinate,
such that both actual endpoint families obey

```math
J_uH'J_u^\dagger=(-1)^{\epsilon_u}H',\qquad
C_uH'C_u^\dagger=(-1)^{\epsilon_u}H'.              \qquad(6)
```

Their literal actions agree on the entire retained middle algebra, so
the two signs in (6) agree as well. Multiplication in the commuting
probe family makes epsilon a character. The endpoint images of H'
therefore have constant computational translations. The same is already
true of the genuine center, which is fixed pointwise.

Equation (5) is a whole-family statement. Choosing a separate v for
each probe would not supply the fixed operator H' in (6).

## 3. Exact factorization and the remaining middle cost

All controlled factors commute, so the replacement is literal:

```math
H^{b(z)}=(H')^{b(z)}P_v^{b(z)},\qquad
\widetilde E_B=P_v^{b(z)}E_B.                       \qquad(7)
```

No scalar is omitted from (7). The new middle is a commuting C3 phase:
in fixed coordinates diagonalizing the joint isotropic algebra
$`R+L_B`$, the extra sign has Boolean phase b times the affine
eigenbit of the signed P_v. This retains any affine part of b on
center directions outside L_B; only its quadratic support was proved
to lie in L_B.

Its whole probe family has actual Clifford defects. To see the signs,
write

```math
J_uHJ_u^\dagger=(-1)^{f_u(z)}H,\qquad
J_uP_vJ_u^\dagger=(-1)^{r_u(z)}P_v.
```

The affine bits satisfy $`f_u+r_u=\epsilon_u`$. With the original
literal identities $`K_u=C_uJ_u^\dagger=D_u(-1)^{bf_u}`$, the
corrected defect is exactly

```math
\widetilde E_BJ_u\widetilde E_B^\dagger J_u^\dagger
 =D_u(-1)^{br_u}
 =K_u(-1)^{\epsilon_u b}.                          \qquad(8)
```

The last factor is a diagonal Clifford because b is quadratic.
Thus (8) retains the original scalars and holds simultaneously, not
only as a statement about cubic symbols.

The actual diagonal-sublist extraction and signed Clifford alignment
from the partial-support theorem can now be applied to the support of
(H')^b. They use original endpoint sublists, at most quadratic Clifford
coordinate cost, no extra helper, and their actual reversed inverse
words. The commuting middle need not be native for these preliminary
steps. What remains is a full-space complex diagonal sandwich with
two retained independent endpoint bands and the commuting C3 middle
\widetilde E_B, together with retained diagonal C3 factors.

The first unsupported inference is to apply the three-native-band
degree theorem to that sandwich. Equation (8) does not make
\widetilde E_B an independent native list. Nor has a uniform degree
theorem been proved here for an arbitrary commuting C3 middle between
two native endpoint bands. The following obstruction shows that simply
repairing the corrected middle as one native band cannot be assumed.

## 4. A corrected middle can require more than one native band

This is an exact middle-algebra model, not a completed Boolean source.
Take m independent computational Z axes, where $`m\equiv2\pmod4`$
and $`m\ge6`$, and put

```math
p=(1,\ldots,1),\quad O=I+pp^T,\quad
P=p\cdot x,\quad A=x_m,\quad
\kappa=\frac{m-2}{4}\pmod2.
```

The native middle is $`\tau_B(x)=\omega^{|x|}`$. Define the exact
Boolean quadratic

```math
b(x)=\kappa+e_2(x)+PA+A.                           \qquad(9)
```

It is invariant under $`x\mapsto x+p`$. Its polar symbol is
$`\beta=e_2+p\wedge e_m`$, with kernel spanned by p and e_m
and rank m-2. Its support is disjoint from the displacement line
spanned by p. Directly, O complements x only when P=1, and the native
phase ratio is

```math
\frac{\tau_B(x)}{\tau_B(Ox)}=(-1)^{P b(x)}.        \qquad(10)
```

The constant kappa matters: it is one for m congruent to six modulo
eight and zero for m congruent to two modulo eight.

Every fixed preimage v of p under $`N=pp^T`$ has $`p\cdot v=1`$.
The corrected middle associated with such a v is

```math
D(x)=\tau_B(x)(-1)^{b(x)(v\cdot x)}.               \qquad(11)
```

**One-layer obstruction.** D cannot be written as C_1 E C_2 on the
same full register, where C_1,C_2 are arbitrary actual Cliffords and
E is an independent commuting signed Pauli T/T-dagger band.

For proof, rewrite a proposed expression as $`D=E'C`$ by conjugating
the independent band through C_1. Every Pauli in the symplectic
orthogonal complement of the B span is fixed by D. A one-band Pauli
image can be Pauli only with empty native activation. Explicitly,
$`(E')^\dagger Q E'=CQC^\dagger`$ is Pauli for each such Q.
Applying inverse-band rigidity shows that every axis of E' lies in
the B span and
that C is a Clifford in its commuting algebra. This also covers any
untouched spectator wires on the same full register.

Since the Boolean factor in (11) squares to one,
$`D^2=\tau_B^2`$. The square of the diagonal Clifford C is a
Pauli up to its literal scalar. Binary displacement of this equality
therefore gives $`AA^T=I`$ for the matrix of axes of E' in the
original B coordinates. Independence makes A a full orthogonal basis.
Every such basis axis has dot product one with p. Thus complementing
x by p flips every original B eigenbit and every E' eigenbit.

The two native complement products
$`\tau_B(x+p)\tau_B(x)`$ and $`E'(x+p)E'(x)`$ are constants,
including arbitrary signed axes and T orientations. For a diagonal
Clifford, $`C(x+p)C(x)`$ is a constant times an affine Boolean sign.
Consequently the Boolean correction f in (11) would have affine
derivative $`\Delta_p f`$. But b is p-invariant and $`p\cdot v=1`$,
so

```math
\Delta_p\bigl(b(x)(v\cdot x)\bigr)=b(x),
```

which has polar rank m-2. This contradiction proves the claim.
All global scalar phases cancel only through the displayed constant
products; none is silently discarded. The result rules out one native
layer with Clifford endpoints, also after tensoring identity on any
returned dirty helpers. No initialized-helper extension is claimed.
It does not rule out a different global three-layer realization of a
completed query.

## 5. A rank-one displacement has only two quotient labels

There is an additional useful classification, still on the original
middle list. Suppose some N_u has rank one and its flag p_u is nonzero.
Symmetry and nilpotence give $`N_u=pp^T`$ for an even-weight vector
p. Equation (4) identifies its nonzero flag with p. Let I be the support
of p, of even size w.

The native cubic tensor of this transvection is
$`e_3(I)=p\wedge e_2(I)`$. Hence (3) implies

```math
\beta=e_2(I)+p\wedge a                             \qquad(12)
```

for some a. Because S_beta is disjoint from the displacement line,
a has no component outside I: otherwise an external covector would
map under beta to the nonzero vector p. Also $`p\cdot a=1`$;
if it were zero, applying beta to p would again place p in its image.
On the I coordinates its matrix is

```math
I+pp^T+pa^T+ap^T.
```

Its kernel is exactly the span of p,a. Therefore

```math
\mathrm{rank}\beta=w-2,\qquad
S_\beta=\langle p,a\rangle^\perp\cap\mathbb F_2^I.
                                                               \qquad(13)
```

Modulo S_beta, each original axis in I has one of the two nonzero
labels with coordinates $`(1,a_i)`$. Both occur an odd number of
times, since w is even and a has odd weight. If the entire simultaneous
displacement space U is one-dimensional, every nonzero N_u is this
same transvection, so (13) describes the entire affected middle
sublist.

This is a global label statement. On a fixed S_beta eigenvalue, the
axes collapse to two labels with their actual central signs. Their
integer T exponents must still be retained modulo eight. Replacing
these fibers by two free rotations would discard those phases and
would not produce one globally independent native list. Section 4
already shows that this distinction can survive even at rank one.

## 6. A literal clean-returning source still retains a complex phase

The following source realizes the external-center mechanism with
independent physical bands. It supplies a useful limit on what the
retained matching conditions alone establish. It is not Boolean.

Use nine physical wires a,c,z_1,...,z_4,r,q,d, with a the sole initialized
wire. Write a bit label also for its positive computational Z Pauli,
and put

```math
h=a+c,\qquad w=\sum_{i=1}^4z_i,\qquad
p=h+r,\qquad y=h+q.
```

The original clean-active lists and negative-middle list are

```math
\begin{aligned}
A_a&=\{a\}\cup\{a+z_i:1\le i\le4\},\\
C_a&=\{c\}\cup\{c+z_i:1\le i\le4\},\\
B_-&=\{a+w,c+w\}.
\end{aligned}
```

The six positive-middle axes are

```math
B_+=\{y+z_i:1\le i\le4\}\cup\{y,y+p+w\}.
```

Every original rotation is T except the displayed y rotation, which
is T-dagger. The fresh q is essential to the middle independence:
the positive span contains y rather than h, and h is already the sum
of the two negative-middle labels.

Let

```math
C_0=S_aS_w\,\mathrm{CZ}_{c,w},\qquad
b'=1+w+e_2(z),\qquad H=Z_h.
```

The active transition satisfies the exact full-space identity

```math
E_{C,a}^\dagger E_{B,-}E_{A,a}=C_0H^{b'}.          \qquad(14)
```

Choose original A-silent axes
$`\{X_cX_r,(X_cX_r)Z_p,Z_d\}`$. The original C-silent band is
the C_0 conjugate of the list
$`\{X_aX_r,(X_aX_r)Z_p,Z_{d+p}\}`$, including its actual Pauli
signs. Set

```math
K_A=\mathrm{CX}_{a,c}H_a,\qquad K_C=C_0K_A,
\qquad V=K_C^\dagger E_C^\dagger E_BE_AK_A.
```

Each physical band has eight independent commuting axes. Both
interior clean ranks are five, the middle has no clean pairs, and
the original active intersection has rank four. All logical inputs
and any attached references are arbitrary.

In the occupied logical coordinates write t=h and
$`Q=X_tX_r`$. Then $`p=t+r`$, $`y=t+q`$, and the two core endpoint
bands agree as the exact operator

```math
E=T_QT_{QZ_p}.
```

It is S_Q on the p=0 sector and exactly $`\omega I`$ on p=1.
For $`D=H^{b'}E_{B,+}`$, a direct parity calculation gives
$`[E,D]=0`$. The occupied operator is consequently diagonal:

```math
G=E^\dagger DE\,T_dT_{d+p}^\dagger
  =\omega^{\phi},
```

with the literal exponent

```math
\phi=\sum_i z_i+w+4q b'
       +2p(d-w-y+2wy)\pmod8.                      \qquad(15)
```

Here w,p,y in (15) are their Boolean parity values, while
$`\sum_i z_i`$ is the integer sum. The full clean decomposition is
$`V=\mathrm{diag}(G,GL)`$, where L is $`-iQ`$ on p=0 and I on
p=1. Thus clean return holds on every occupied column, and the unused
sector is specified rather than replaced by a convenient extension.

At the zero input phi is zero. At
$`t=r=q=d=0`$ and $`z=(1,1,0,0)`$, it is two. Hence the two
occupied phases include 1 and i, so no common scalar makes this a
Boolean sign. The residual failure is already visible modulo four;
the common active quadratic has not disappeared merely because the
core conjugation became diagonal.

The emitted coordinate words use 24 T/T-dagger gates in three native
bands, one clean wire and eight arbitrary wires, with no helper. The
forward and actual inverse words each contain 223 Clifford gates when
S-dagger is one gate. Expanding S-dagger as three S gates gives 231
and 227 elementary Cliffords respectively. These are counts of the
explicit unoptimized words, not a claim of minimal synthesis cost.

The [exact fixtures](../../tests/test_external_center_middle_correction.py)
check all 512 physical columns and the actual inverse using those
compiled words. They also check the six-axis transvection, its native
phase defect, the corrected square, and the quadratic complement
derivative underlying Section 4. A further fixture checks the complex
diagonal Clifford alignment and the distinction between whole-residual
identity-mode Hadamards and their action on separate factors. The
analytic arbitrary-width and
whole-family claims are the proofs above, not extrapolations from
those checks.

## 7. Bounded total displacement gives a uniform degree bound

Return to the actual Boolean source and let $`d=\dim U`$, where U
is the total displacement space in the original independent middle
coordinates of Section 1. Put

```math
s=\frac{d(d+1)}2,\qquad k=s+d+1,
\qquad C(d)=D_k+72\,2^k,                           \qquad(16)
```

with D_k from the
[partial-diagonal finite-feature theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants).
At polar rank at least four, the hypotheses of Sections 1--2 imply

```math
\deg g\le C(d).                                   \qquad(17)
```

This is independent of the number of original middle axes, original
clean-active ranks, central controls and dirty wires. In particular,
the rank-one case of Section 5 closes even when its affected native
list and its control polar rank grow. The proof uses no center-fiber
collapse and does not replace the middle by one newly native band.

### 7.1. Peel actual diagonal endpoint sublists

Peel from each endpoint every original axis whose pullback through
that endpoint Clifford has zero computational X mask. They are actual
diagonal native C3 sublists D_A,D_C. Write the resulting exact word as

```math
V=D_C^\dagger V_1D_A,
\qquad V_1=K_C^\dagger C^\dagger F E_B A K_A,
\qquad F=H^b.                                     \qquad(18)
```

V_1 is complex diagonal; Booleanity is not assumed after these peels.
The probe families J_u,C_u and their middle actions are unchanged,
since the removed input/output factors commute all computational Z
probes. Every remaining endpoint axis has a nonzero X-mask functional
lambda_i(u).

The physical Pauli space corresponding to $`K=U^\perp\subset L_B`$,
using the original native-coordinate dot product, is fixed by every
middle action. For k in that space, the native formula
is

```math
0=N_uk=\sum_i\lambda_i(u)\langle a_i,k\rangle a_i.
```

Native independence and nonzero lambda_i imply that every remaining
A axis commutes K. The same argument applies at C. Both F and E_B
also commute K. Thus V_1 intertwines the actual input and output
Pauli families on K.

Use the literal diagonal-Pauli-intertwining construction from
[constant-translation closure](CONSTANT_TRANSLATION_NATIVE_CLOSURE.md#3-match-the-remaining-pauli-families-by-an-actual-diagonal-clifford)
to obtain one diagonal S/CZ/Z Clifford L with the same transport.
Then $`V_0=V_1L^\dagger`$ is complex diagonal and centralizes the
common Pauli family. This uses actual full-space Pauli identities,
including their signs, rather than a Boolean reality condition or an
initialized-sector extension.

### 7.2. Normalize the fixed family without conjugating restored factors through Hadamards

Normalize the common K family jointly with the genuine computational
center R, using a fixed affine computational permutation and diagonal
Clifford decorations T. After this change, V_0 commutes selected X
operators and is diagonal, so it is exactly identity on their tensor
modes I. Hadamards on those modes give

```math
H_I(TV_0T^\dagger)H_I=TV_0T^\dagger.              \qquad(19)
```

The transformed native word therefore represents the same diagonal
residual in the affine T frame, and K is now a computational Pauli
subgroup. The original genuine center stays computational: it commutes
the fixed K family and has no component on the X modes being changed.

Keep the restored factors D_A,D_C,L conjugated only by T. Equation
(19) is an equality of the whole residual; it does not assert that
its separate factors are identity on I or permit their separate
Hadamard undo. In particular the three restored factors remain actual
diagonal C3 phases in this frame.

For the new computational probes, those outside I are old
computational-Z combinations after the affine/diagonal change. The
new probes on I belong to the fixed K family, possibly multiplied by
genuine-center Paulis. They activate neither remaining endpoint list
and commute H. Thus the following common matrix-entry bound and the
simultaneous identity (5) remain valid in these coordinates.

### 7.3. Bound the quantum masks while keeping every original middle axis

The symmetric linear maps N_u have image in U. Their matrix-entry
linear functions of u therefore span a space Gamma of dimension at
most $`s=d(d+1)/2`$. No orthogonal change of the native basis is
used in this count. These functions are linear because the common
isotropic native encoder gives $`N_uN_{u'}=0`$, and hence
$`N_{u+u'}=N_u+N_{u'}`$.

If an actual remaining A axis a_i is active against some middle axis
B_j, its coefficient in $`N_uB_j`$ is exactly lambda_i(u).
Independence of the A list therefore recovers lambda_i from the common
N_u matrix entries. Its input X mask lies in Gamma. At C the identical
middle action gives the same bound in the same Gamma.

An axis active against H is already active against the middle span.
Indeed, compare its coefficient in the full native identities
$`(\Phi_{J_u}-I)h=N_uv`$ and the corresponding C-family identity,
which are (5) expressed in the full physical label space. This does
not apply the middle-coordinate map N_u to h outside its domain.
If the axis commuted
every B axis, the right coefficient would vanish for all u. Since its
lambda_i is nonzero, it would commute H as well. The same argument
holds at C. Only binary displacements are used here; the literal
constant sign difference of Section 3 is retained in the operators.

Split the remaining actual lists into active and passive sublists:
passive axes commute every B axis and H. In the normalized endpoint
frames define

```math
\begin{aligned}
A_{\rm in}&=K_A^\dagger A_{\rm act}K_A,&
C_{\rm out}&=K_C^\dagger C_{\rm act}K_C,\\
D_{\rm in}&=K_A^\dagger(FE_B)K_A,&
W_2&=K_C^\dagger C_{\rm pass}^\dagger A_{\rm pass}K_A.
\end{aligned}
```

The residual factors in the exact order

```math
TV_0T^\dagger
 =C_{\rm out}^\dagger W_2 D_{\rm in}A_{\rm in}.    \qquad(20)
```

Here W_2 is an actual two-native-band word with its actual Clifford
endpoints. Moving the passive sublists uses their axiswise commutation
with F and E_B. The full F E_B is kept as one commuting C3 factor;
neither E_B nor its corrected version is split into conditional native
lists.

The pullback of K is computational. Thus the X masks of the whole
middle span are generated by at most d representatives of L_B/K.
The support of D_in adds at most the one X mask of H; its b controls
are genuine computational parities. Choose a computational-Z subgroup
M_0 annihilating these at-most-d+1 masks and Gamma. It has codimension
at most k. The factors A_in,C_out,D_in fix M_0 pointwise. Since the
left side of (20) is diagonal, W_2 also fixes M_0 pointwise on the
entire register.

### 7.4. Full matrix features give the Boolean degree bound

The existing partial-diagonal theorem applies to W_2 with an empty
third band. It determines every exact core matrix from the features
counted by D_k, including its complex phases and matrix entries that
are not roots of unity. Besides W_2, the full product contains six commuting
C3 phase factors: A_in, C_out, D_in, and the three restored factors
D_A,D_C,L, with their actual inverses and order. Each factor commutes
M_0; the factors are not asserted mutually commuting.

For any one commuting C3 factor fixing M_0, its commuting core Pauli
projections can be diagonalized by one fixed core Clifford. Its phase
is a weighted cubic in core bits and affine central parities. The
three binary digits of that exponent modulo eight, evaluated at each
of the at-most-$`2^k`$ core strings, give at most $`3\,2^k`$
Boolean features of degree at most four. This is the retained
[commuting-C3 feature calculation](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#5-the-residual-and-silent-factors-have-bounded-features-too),
with all scalar phases kept.

Multiplying the finite matrices in their actual order and selecting a
diagonal Boolean entry adds at most
$`6\cdot3\cdot4\cdot2^k=72\,2^k`$ to the D_k degree count.
This proves (17). Returning from T only composes the feature functions
with an invertible affine input change; diagonal decorations commute
the Boolean target. Their Boolean degrees do not increase. No
individual factor is conjugated back through H_I.

The actual peels, coordinate changes and diagonal Clifford alignment
have $`O(Q^2)`$ additional elementary Clifford cost, with their actual
inverse words charged separately and no added work. The proof is a
degree factorization, not a new shallow implementation. Original
central cubic and real quadratic signs do not raise this bound.

At control polar rank two, the retained finite-insertion theorem uses
support at most three after paid affine-control absorption. At rank
zero the insertion is Clifford and the ordinary three-native-band
bound is five. Consequently the actual no-pair branch has, for fixed
total displacement dimension d, the all-rank bound

```math
\deg g\le
\max\{C(d),D_{\rm ins}(0,0,0,3)\}.                \qquad(21)
```

For a literal XOR query, check the hypotheses in its actual
output-Hadamard sign representation first; its predicate degree is
then at most one less than (21). All statements retain arbitrary
dirty inputs, their references, full clean return and literal phases.

## 8. The remaining question has growing total displacement

The external displacement is absorbed by one simultaneous middle
Pauli product, and bounded total displacement now closes without
repairing the corrected middle into an independent native band.
The remaining actual no-pair case allows U to grow together with the
original common active support. Its full Boolean phase must still be
matched with the original global lists.

Constant translations of H' do not imply constant translations on
the entire middle support. Nor does the rank-one quotient give free
center-dependent rotations. Section 4 excludes a generic one-layer
repair, while Section 7 avoids that step through a fixed-dimensional
full-matrix argument. Section 6 shows why a diagonal core matching
still requires a literal Boolean phase check.

This revision does not establish the unrestricted three-layer degree
target, a new complete-frame resource frontier, a Boolean counterfamily,
or a reduction of the arbitrary-middle-pair branch. Every original
native list, actual inverse, dirty/reference contract and central phase
remains retained.
