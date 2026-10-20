# All-active original-pair closure with arbitrary singleton coupling

[Original singleton-coupling reduction](ORIGINAL_PAIR_SINGLETON_COUPLING.md) ·
[Occupied computational center](OCCUPIED_BOOLEAN_PAIR_CENTER.md#4-computational-center-and-a-finite-feature-description) ·
[Partial-diagonal native features](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md)

The current [structural diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)
sharpens this same all-active source class to sign/XOR degree six/five
uniformly, and five/four when the common active rank is at least four.
It uses Sections 2--5 below and replaces the large feature bound in
Sections 1 and 6. Those historical bounds and their proof are retained;
the generic partial-diagonal theorem is unchanged. See the
[revision checkpoint](STRUCTURAL_PAIR_REVISION_RESULT.md) for the current
frontier and remaining general-pair question.

The actual scalar branch of one original physical pair supplies native
matching on one common stabilizer code. The other branch then forces a
single dressed Clifford to have only one computational X mask. This
removes the bounded singleton-to-silent coupling hypothesis in the
all-active genuine one-pair class. It does not bound that coupling rank,
replace an original native list, or settle a general original pair.

## 1. Original physical hypotheses and uniform bound

Use the literal one-clean three-native-band word

```math
\begin{aligned}
V&=K_C^\dagger E_C^\dagger E_BE_AK_A,\\
P_A&=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger.
\end{aligned}
```

Each phase list consists of independent commuting signed Hermitian
Paulis with its actual T/T-dagger orientations. Assume exact clean
return and complete occupied action

```math
G=e^{i\theta}\,\mathrm{diag}((-1)^{g(x)})
```

on the entire logical and dirty register, with arbitrary references.
The scalar is retained. The unused clean block is unrestricted.
Suppose the original middle list has exactly one genuine physical pair
B,D and its difference label $`r=b+d`$ is the sum of **all actual
clean-active axes at each endpoint**. Original singleton lists and
both silent lists are unrestricted.

The historical feature-theorem conclusion, now superseded by the
[structural six/five bound](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md), is

```math
\begin{aligned}
\deg g&\le 37748865,\\
\deg f&\le 37748864
   &&\text{for a complete Boolean XOR query }Q_f.
\end{aligned}
```

This constant is $`D_1+36`$, using the retained
[partial-diagonal constant](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants)
$`D_1=37748829`$. At this historical feature step, the six/five bound
for common active rank two and the 2428/2427 one-sided bound were sharper
in their respective subclasses. No bound on the actual
singleton-to-silent coupling rank h is asserted or needed.

The proof keeps all original phase gates and their literal inverses.
Its restriction to the physical clean code occurs only after every
retained internal factor preserves that code. A further one-Pauli
restriction is used to prove matching, not to supply an implementation
or replace the complete occupied target by a branch promise.

## 2. A literal central-controlled quarter-turn

Use the exact orientation repair and common active-span conclusions of
[the original-list reduction](ORIGINAL_PAIR_SINGLETON_COUPLING.md#2-literal-repair-and-the-common-active-span).
Write $`J=E_AP_AE_A^\dagger`$, $`JBJ=sD`$ and let epsilon,delta
be the original pair orientations. The actual factorization is

```math
\begin{aligned}
F&=T_B^\epsilon T_D^{s\epsilon},\qquad
C=T_D^{\delta-s\epsilon},\\
E_B&=CFE,
\end{aligned}
```

where E retains every original singleton and C is Clifford. Absorb C
into the output list and encoder literally. Every singleton is then
positive, the two clean images agree, and

```math
\begin{gathered}
L_A=L_C=L,\qquad p_A+p_C\in L,\\
\langle b,l\rangle=\langle p_A,l\rangle\quad(l\in L),\qquad
R=L\cap p_A^\perp.
\end{gathered}
```

The genuine physical center R commutes both clean Paulis and every
actual retained phase axis. Let $`k=\dim L`$. It is even and at least
two: b pairs with every actual active axis, while b commutes the pair
difference r, their sum.

There is an actual Clifford $`C_0\in\mathrm{Alg}(L)`$ taking
$`P_A`$ to $`P_C`$ with its signed action. Indeed, products of
S rotations on the actual A-active axes realize any required label
shift in L of $`P_A`$; an active Pauli corrects the remaining sign.
This Clifford fixes every silent axis literally. Since all three
factors below lie in the same commuting algebra, the operator

```math
D_R=C_0^\dagger E_{C,a}^\dagger E_{A,a}
```

is diagonal C3 in fixed L coordinates. Both
$`E_{C,a}^\dagger E_{A,a}`$ and $`C_0`$ take $`P_A`$ to
$`P_C`$, so $`D_R`$ commutes $`P_A`$ and belongs to
$`\mathrm{Alg}(R)`$.

Choose one actual active A label a, put

```math
\beta=\langle p_A,b\rangle,\qquad h=b+\beta a,
```

and choose a signed Hermitian representative $`P_h`$. It commutes
$`P_A`$ and R. Define the actual central Pauli and its projectors
by

```math
\Gamma=sBD,\qquad
\Pi_\pm=\frac{I\pm\Gamma}{2},\qquad
\eta=\frac{1-s}{2}.
```

The native convention $`T_{-P}^\epsilon=\omega^\epsilon T_P^{-\epsilon}`$ gives the full physical identity

```math
F=\omega^{-\epsilon\eta}
       \left(\Pi_+S_B^\epsilon+
                         \Pi_-\omega^\epsilon I\right),
\qquad \omega=e^{i\pi/4}.
```

Thus the minus branch is literally scalar; its scalar is not omitted.
Let $`H=E_{A,a}^\dagger BE_{A,a}`$. On the plus branch H commutes
$`P_A`$. Its Pauli support there is contained in $`h+R`$, and
there is one global Hermitian diagonal Clifford
$`Q_R\in\mathrm{Alg}(R)`$ such that

```math
\Pi_+H=\Pi_+Q_RP_h.
```

Section 3 proves this assertion, including its central quadratic rank.
Only the displayed plus-branch equality is used. H is not replaced by
$`Q_RP_h`$ on the other branch.

The complete clean-preserving insertion therefore has the exact
full-physical-space form

```math
\begin{aligned}
G_0&=C_0^\dagger E_{C,a}^\dagger F E_{A,a}\\
&=D_R\omega^{-\epsilon\eta}
  \left(\Pi_+S_{Q_RP_h}^\epsilon+
                         \Pi_-\omega^\epsilon I\right).
\end{aligned}
```

It commutes every original positive singleton. A sector-dependent
Clifford frame, a free controlled phase and an extension of the target
into its unused clean block have not been introduced.

## 3. The central quadratic sign has rank k minus two

Every actual active A axis anticommutes with B. Keeping signed axes
and orientations, direct native conjugation gives

```math
H=\omega^{\sum_i\epsilon_i}
       \left(\prod_i S_{A_i}^{-\epsilon_i}\right)B.
```

Consequently $`HP_h`$ is a global diagonal Clifford in
$`\mathrm{Alg}(L)`$. It need not be Hermitian outside the plus
branch. Choose one fixed Clifford coordinate system for the actual
independent A-active axes, with their eigenbits $`x_1,\ldots,x_k`$.
The condition $`\Gamma=+1`$ is one affine parity equation in x,
and $`P_A`$ translates x by the all-ones vector. Since k is even,
that translation preserves the parity equation.

On this affine parity hyperplane, $`HP_h`$ is Hermitian and invariant
under the $`P_A`$ translation. Its sign is a quadratic Boolean
function q. Up to affine terms and a constant, it is half the
integer sum $`\sum_i\epsilon_i x_i`$ on the fixed parity class.
On the tangent space

```math
\mathcal E=\{v\in\mathbb F_2^k:\ \sum_i v_i=0\},
```

its polar form is the ordinary dot product. The signs epsilon change
only affine terms. Since k is even, the radical of this restricted
dot product is exactly the all-ones line. Descending the actual
$`P_A`$ translation therefore gives

```math
\mathrm{rank}(q_2)=k-2
```

on the $`k-2`$ remaining central coordinates of the plus branch.
Here $`q_2`$ denotes the quadratic polar form, independently of the
chosen affine section.

Take one fixed linear section of this translation quotient. The
quadratic sign on that section is a polynomial in the characters of
$`R/\langle\Gamma\rangle`$. Extending this same polynomial
independently of the Gamma bit produces the global $`Q_R`$ used
above. It is a Hermitian diagonal Clifford in $`\mathrm{Alg}(R)`$.
Changing a signed Pauli representative changes only an affine term;
the quadratic rank is independent of the fixed quotient coordinates.
In particular, when $`k\ge4`$, its
sign on the plus branch is genuinely nonaffine. The $`k=2`$ case
retains the already proved rank-two endpoint bound.

## 4. Scalar-branch matching on one common native code

After the preceding full physical operations, the word is

```math
V=(K_C')^\dagger E_{C,s}^\dagger G_0 E E_{A,s}K_A,
\qquad K_C'=C_0^\dagger K_C.
```

Every retained internal factor now preserves $`P_A`$, and both
endpoint encoders send the physical clean Z to $`P_A`$. Restrict
this exact word to its positive clean code. The retained silent and
positive singleton lists are still independent native lists: an
actual active axis commutes each such list and anticommutes
$`P_A`$, ruling out a clean-quotient relation.

Use the proved
[shared computational-center alignment](OCCUPIED_BOOLEAN_PAIR_CENTER.md#4-computational-center-and-a-finite-feature-description).
It keeps the physical clean coordinate until restriction, retains the
extracted quadratic phase, and makes the logical R a computational-Z
subgroup at both endpoints. Rename the resulting logical endpoints
and bands $`K_A,K_C,A,C,E`$. The residual target remains Boolean
diagonal on its entire logical register. For every computational
Z label z, define

```math
\begin{aligned}
J_z&=AK_AZ(z)K_A^\dagger A^\dagger,\\
C_z&=CK_CZ(z)K_C^\dagger C^\dagger.
\end{aligned}
```

These are actual Hermitian Cliffords fixing the same signed R
pointwise. Full diagonality gives

```math
EJ_zE^\dagger=G_0^\dagger C_zG_0.
```

The logical Gamma is nonzero and is outside the span of the actual
positive singleton axes. To check this without discarding the physical
clean coordinate, a putative physical relation

```math
r+t=e p_A,\qquad t\in\mathrm{span}(E),
```

pairs with any actual active a to force $`e=0`$. It would then
contradict independence of the original middle list B,D and all its
singletons. The same argument, using original endpoint active/silent
independence and $`r\in L`$, shows Gamma is outside either silent
span after clean restriction.

Now restrict to **one fixed** logical stabilizer code
$`\Gamma=-1`$. The positive singleton list remains genuinely
independent, so its restricted phase band is native. Both comparison
Cliffords descend to actual Cliffords. Since $`D_R`$ commutes
$`C_z`$, the scalar branch gives exactly

```math
E_-J_{z,-}E_-^\dagger=C_{z,-}.
```

The actual [one-band matching theorem](ONE_BAND_AXIS_MATCHING.md)
therefore supplies a permutation of the projected singleton axes.
Lifting to the full clean-reduced logical label space gives

```math
\Phi_z b_i=b_{\pi_z(i)}+e_{z,i}r,
\qquad e_{z,i}\in\mathbb F_2,
```

where $`\Phi_z`$ is the binary action of $`J_z`$. This is a
simultaneous matching modulo the single actual pair difference.
Orientations and signs are retained; native matching permits them.

The same binary projected matching holds on the other code
$`\Gamma=+1`$. Its projected E list is also independent. The
converse direction of one-band matching thus proves

```math
S_{z,+}:=E_+J_{z,+}E_+^\dagger
\quad\hbox{is an actual Clifford}.
```

Its literal scalar is retained, as in the matching theorem.

No claim that $`EJ_zE^\dagger`$ is globally Clifford follows at
this point: the possible central lifts $`e_{z,i}r`$ have not been
removed. Only the two code restrictions are Clifford. Restricting all
of R at the matching step would introduce unbounded native relations;
that different step is not used.

## 5. Quadratic signs force the actual dressed axis to be monomial

Assume $`k\ge4`$. Work in one fixed Clifford frame on the plus
code which makes the remaining genuine R coordinates computational Z.
No frame is chosen separately for individual central characters.
If $`P_h`$ becomes central after physical clean restriction, every
$`C_z`$ already fixes its label and the conclusion below is immediate.
Otherwise its nonzero quotient direction is retained throughout.
Write chi for these characters. A fixed Pauli representative $`P_h`$
has character dependence given by one affine sign. Absorb it into q;
the nonzero quadratic part from Section 3 is unchanged.

Modulo central scalar, the insertion on this fiber is

```math
F_\chi=S_{(-1)^{q(\chi)}\tau_h}^{\epsilon}.
```

An actual Clifford fixing R pointwise has a central-character-independent
quotient symplectic map, and its signs on fixed quotient Pauli
representatives are affine in chi. More explicitly, in this one fixed
frame its fiber is a central scalar times $`P_{a(\chi)}C_*`$, where
$`C_*`$ is fixed and the ordered Pauli displacement $`a(\chi)`$ is
affine. Both $`C_{z,+}`$ and
$`S_{z,+}`$ are such actual Cliffords. The elementary identity
$`S_{-P}=iPS_P`$, with its inverse version, shows that replacing the
fixed quarter-turn $`S_{\tau_h}^\epsilon`$ by $`F_\chi`$ in

```math
S_{z,+}=F_\chi^\dagger C_{z,+}F_\chi
```

adds, after a fixed Clifford change, a Pauli displacement

```math
q(\chi)\bigl(h+\overline C_z h\bigr).
```

Here $`\overline C_z`$ is the quotient symplectic map. If the
parenthesized label were nonzero, some fixed quotient Pauli would
acquire a nonzero quadratic sign q. This cannot equal the affine
central sign of an actual Clifford. Hence

```math
\overline C_z h=h\qquad\hbox{for every }z.
```

Lifting first through the plus code and then through the remaining
central quotient says that the full logical displacement of h by
$`C_z`$ belongs to R. Its actual native formula also places it in
the projected C-silent span $`U_C`$:

```math
(C_z-I)h
  =\sum_j\lambda_j(z)\langle c_j,h\rangle c_j
  \ \in\ U_C.
```

This last display denotes binary actions. Crucially,

```math
U_C\cap R=0
```

in the clean-reduced space. Physically, original output active/silent
independence gives $`U_C\cap L=0`$. A relation
$`u+r'=e p_A`$ introduced by clean restriction again pairs with an
actual active a to force $`e=0`$, so it introduces no intersection.
The Clifford $`C_0\in\mathrm{Alg}(L)`$ fixed the silent list
literally; no replacement list is involved. Therefore

```math
C_zP_hC_z^\dagger=\pm P_h
\qquad\hbox{on the entire clean-reduced logical space}.
```

Define the actual dressed Hermitian Clifford

```math
H_{\rm out}=K_C^\dagger C^\dagger P_hCK_C.
```

Conjugating the preceding equality gives

```math
Z(z)H_{\rm out}Z(z)=\pm H_{\rm out}
\qquad\hbox{for all }z.
```

Every Pauli in its nonzero Pauli expansion consequently has the same
computational X mask u. Equivalently,
$`H_{\rm out}=X(u)D`$ with D an actual diagonal Clifford and its
literal scalar retained. The case u equal to zero is allowed. This
bounds the dressed action itself, without bounding the number, rank
or span of coupled singleton or silent axes.

## 6. One shared matrix ledger completes the degree bound

This retained feature step is superseded for the present source class by
the [full-space diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md).
Its one-mask input from Section 5 remains part of the current proof;
the generic matrix-feature theorem and the calculation below are
preserved as component results.

The complete logical factorization remains

```math
\begin{aligned}
G&=F_{\rm out}W_3,\\
W_3&=K_C^\dagger C^\dagger EAK_A,\\
F_{\rm out}&=K_C^\dagger C^\dagger G_0CK_C.
\end{aligned}
```

The three lists of $`W_3`$ are the actual independent projected
native lists, with every original singleton and silent factor retained.
All central coefficient functions commute C. Put
$`R_{\rm out}=K_C^\dagger RK_C`$ and use the corresponding
output notation for $`D_R,Q_R,\Gamma`$. The exact insertion is

```math
F_{\rm out}=D_R\omega^{-\epsilon\eta}
  \left(\Pi_+S_{Q_RH_{\rm out}}^\epsilon+
                         \Pi_-\omega^\epsilon I\right).
```

Its Pauli support has only the X masks zero and u. Choose one fixed
invertible affine computational change which puts u on at most one
core bit. Let M be the computational-Z subgroup on all remaining
logical bits. Both $`F_{\rm out}`$ and the diagonal G fix M
pointwise, hence so does the actual full-space $`W_3`$. The
partial-diagonal theorem applies with core size at most one.

The insertion needs only nine degree-at-most-four Boolean features in
these same central computational coordinates:

- Three exponent digits of the diagonal-C3 phase $`D_R`$, of
  degrees at most one, two and four.
- The affine Gamma eigenbit and the quadratic sign of $`Q_R`$.
- Two mod-four exponent digits of the diagonal Clifford D at each
  of the two possible core input values, four features of degrees
  at most two each.

These central factors commute $`H_{\rm out}`$, so their functions
are invariant under the X-mask translation and depend only on the
central bits after the fixed affine change. The listed features
determine every exact insertion matrix entry, including zero entries
and the scalar minus branch. With no core bit, the same count is a
valid overestimate. No intermediate coefficient is required to be a
Boolean phase.

Combining this ledger with the retained partial-diagonal features for
$`W_3`$, and only then using the complete Boolean target, gives

```math
\deg g\le D_1+4\cdot9=37748865.
```

As in the existing feature theorem, extend the final sign function
off its attained classical feature vectors and compose its algebraic
normal form. This is a classical extension, not an extension of an
intermediate quantum operator to an unused sector. The shared-center
alignment transports the assertion back; its extracted quadratic
phase does not increase this bound. When $`k=2`$, the retained
[rank-two endpoint theorem](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#5-all-paired-sources-with-two-clean-active-axes-at-both-endpoints)
gives the stronger six/five conclusion and completes the remaining
case. Output Hadamards turn a Boolean XOR query into the sign
$`(-1)^{zf(x)}`$, giving the stated one-degree decrement.

## 7. What this removes, and what remains open

This closes the all-active genuine one-pair class with arbitrary
positive singleton-to-silent coupling. The symmetric coupling block
of the earlier h-dependent proof need not be bounded or diagonalized.
The mechanism instead uses the pair's literal scalar branch, native
matching on one fixed stabilizer code, and the nonaffine central sign
on its other branch. Its final computational X-mask bound applies to
the actual dressed insertion.

All physical coordinate Cliffords and their actual inverses, if
emitted, carry $`O((Q+1)^2)`$ elementary-Clifford costs on the same
register. No initialized helper is added, no native factor is freely
resynthesized, and all complete logical columns, arbitrary dirty and
reference inputs, exact clean return, and the unrestricted unused
clean block are retained. The algebraic proof adds no approximation
error and gives no new upper compiler.

The [small exact fixture](../../tests/test_all_active_singleton_closure.py)
uses four actual active axes and one arbitrary logical wire. It checks
the nonaffine central quadratic, the plus-branch identity without
extending it to the minus branch, both signed-swap conventions, the
literal scalar branch, the complete native word and its actual inverse
on all 32 physical columns. The fixture is not asserted to be a Boolean
query; it guards the phase-sensitive structural identity used above.

There is also the retained original-circuit approximation consequence.
If a native circuit of at most three T layers and at most one clean
qubit approximates a complete Boolean sign or XOR query with initialized
operator error strictly below $`2^{-530690}`$, the
[initialized rounding theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
forces that same circuit to have exact clean return and exact occupied
Boolean action, up to the allowed common eighth-root phase. No exact
return premise is needed before applying that rounding theorem. If the
resulting actual clean comparison and original physical lists satisfy
the all-active one-pair hypotheses of Section 1, the present degree
bound then applies. Thus a higher-degree target has the stated error
floor within this source class. The condition on the original lists
is still required; this is one conservative fixed tolerance, not a claim
excluding all approximate three-band circuits or every fixed accuracy.

A general original pair need not have its difference supported on all
actual active axes at both ends. The nonaffine quadratic rank and
single scalar-branch reduction above are proved only for the stated
all-active class. Several original pairs, unrestricted one-clean
three-layer degree, and dimension-growing native-depth bounds remain
open. The general resource frontier is unchanged.
