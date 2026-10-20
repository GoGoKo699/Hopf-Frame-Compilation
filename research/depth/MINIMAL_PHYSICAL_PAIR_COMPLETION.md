# A minimal physical pair and its required output partner

[Original physical pairs](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[Native profiles and the exported guard](CENTRAL_MODULE_NATIVE_PROFILES.md#62-the-transformed-union-can-grow-in-the-exported-comparison) ·
[Two-native-band classification](NATIVE_BAND_SOURCE_STRUCTURE.md#2-two-native-bands-exact-c3-classification)

An actual one-pair source can cancel a positive middle singleton and
implement occupied identity. Its pair residual also contains a second,
centrally controlled T phase. That partner prevents the minimal source
from realizing the exported growing-orbit cancellation guard: adding a
silent input band that moves the singleton makes every one-band output
completion impossible, even with arbitrary output Cliffords and arbitrary
action on the unused clean sector.

The same physical reduction gives an actual zero-clean three-band word
with just one middle axis for every complete clean-return target in this
fixed family. Its diagonal targets are diagonal $`\mathcal C_3`$:
Boolean sign degree is at most three and XOR degree at most two. These
are results for the specified input and middle word, not a theorem for
all fixed-pair sources or unrestricted depth. Section 5 extends the
scout to every paired source whose two endpoint clean-active ranks are
two, proving sign/XOR degree at most six/five. General inside-span
matching and growing endpoint ranks remain open.

## 1. A literal original one-pair source

Use physical clean wire c, arbitrary wire a, and a nonidentity Hermitian
Pauli B on the remaining arbitrary wires. Define

```math
\begin{aligned}
K&=\mathrm{CX}_{c\to a}H_c, & A&=T_{Z_c}T_{Z_a},\\
P&=X_cX_a, & Q&=Y_cY_a,\\
R_p&=PQ=-Z_cZ_a.
\end{aligned}
```

All phases have the literal convention

```math
T_H=\frac{1+\omega}{2}I+\frac{1-\omega}{2}H,
\qquad \omega=e^{i\pi/4}.
```

Take the actual bands and complete operator

```math
\begin{aligned}
E_A&=A,\\
E_B&=T_{PB}T_{QB}T_B^\dagger,\\
E_C&=A T_{R_pB},\\
V&=K^\dagger E_C^\dagger E_BE_AK.
\end{aligned}
```

The A list has two independent axes. The C list has three independent
commuting axes because B has nonzero support off c,a. The middle list
$`\{PB,QB,B\}`$ is likewise independent and commuting. Thus each
displayed band is physically native.

Both endpoint clean Paulis are P. The actual clean image
$`J=APA^\dagger`$ obeys

```math
JPJ=Q,\qquad JQJ=P,\qquad JBJ=B.
```

Consequently PB and QB form one original physical pair with positive
swap signs, and B is a positive singleton. The middle band commutes J.
The extra C axis is silent, so
$`E_CP E_C^\dagger=J=E_BJE_B^\dagger`$. This proves full physical
clean-observable return before restricting the input.

The pair product commutes A literally: writing it as
$`\omega\exp[-i\pi(P+Q)B/8]`$, the commutator vanishes because
$`P+Q=X_cX_a+Y_cY_a`$ commutes $`Z_c+Z_a`$. Also

```math
K^\dagger PK=Z_c,\qquad
K^\dagger QK=-Z_cZ_a,\qquad
K^\dagger R_pK=-Z_a.
```

In the positive clean code, with $`R_l=-Z_a`$, the occupied pair
residual is therefore

```math
\begin{aligned}
F&=\left.K^\dagger A^\dagger T_{PB}T_{QB}AK\right|_{c=0}\\
 &=T_B T_{R_lB}\\
 &=\Pi_{a=0}\,\omega I+\Pi_{a=1}\,S_B.
\end{aligned}
```

The complete occupied word is exactly

```math
T_{R_lB}^\dagger F T_B^\dagger=I.
```

No common scalar was discarded. The full physical operator, including
its unused sector, is

```math
V=\Pi_{c=0}I+
  \Pi_{c=1}\bigl(\Pi_{a=0}I+\Pi_{a=1}\,iB\bigr).
```

Thus every occupied input and arbitrary reference is returned literally;
the full operator is not asserted to be Boolean identity.

## 2. Genuine center and the partner that remains

In the [occupied-center notation](OCCUPIED_BOOLEAN_PAIR_CENTER.md), both
clean-active endpoint spans are $`L_A=L_C=\langle z_c,z_a\rangle`$.
Their genuine central space is the one-dimensional span of
$`z_c+z_a`$, represented by the signed $`R_p`$. The active localization
and its clean-preserving subspace are

```math
\begin{aligned}
S_0&=\langle z_c,z_a,p+b\rangle,\\
H_0&=S_0\cap p^\perp
    =\langle z_c+z_a,p+b\rangle,\\
\dim(H_0/R)&=1.
\end{aligned}
```

The actual active transition preserves P, so its matching Clifford may
be $`C_0=I`$. After clean restriction, F has support in two R-cosets.
This is an actual original-pair instance of insertion dimension one.

The positive singleton cancels $`T_B`$, leaving $`T_{R_lB}`$.
Cancelling both factors inside the original middle band would require
adding $`T_{R_pB}^\dagger`$, whose proposed four axes satisfy

```math
(p+b)+(q+b)+b+(p+q+b)=0.
```

That four-axis list is dependent and is not one native band. In the
actual identity source, the output silent band pays the remaining phase.

## 3. No one-band output completion for the moving-axis extension

Keep the middle band above and change the actual input to $`AWK`$,
where W is any actual clean-silent subband completing A to one
independent native input band. This includes axes touching c,a.
Commutation with $`Z_c,Z_a,P`$ gives the exact form

```math
W=\prod_i T_{R_p^{e_i}D_i}^{\epsilon_i},
\qquad e_i\in\{0,1\},\quad\epsilon_i\in\{1,-1\},
```

with signed Hermitian D_i acting off c,a. Their binary labels are
independent: a relation among them would place the corresponding W
product in the active span $`\langle z_c,z_a\rangle`$, contradicting
independence of the complete input band. The input clean image is
still J. In the positive clean code the projected band is

```math
W_0=\prod_iT_{R_l^{e_i}D_i}^{\epsilon_i},\qquad R_l=-Z_a,
```

an actual independent native band on the full logical space.
Allow an arbitrary actual independent native output band
$`\widehat E_C`$ and an arbitrary Clifford $`\widehat K_C`$.

**Completion obstruction.** If at least one D_i anticommutes with B,
the complete unitary

```math
\widehat V=
\widehat K_C^\dagger\widehat E_C^\dagger E_BAWK
```

cannot have exact clean return and an occupied Clifford target on the
whole logical register. In particular it cannot have occupied identity.
There is no restriction on its unused-clean-sector action.

**Proof.** Suppose such a completion exists. Full initialized-unitary
return means $`\widehat V`$ maps the entire positive clean subspace
onto itself. Unitarity then preserves its orthogonal complement too,
giving $`[\widehat V,Z_c]=0`$ on the full physical space. Hence the
literal clean-observable comparison is

```math
\begin{aligned}
\widehat P_C&=\widehat K_C Z_c\widehat K_C^\dagger,\\
\widehat E_C\widehat P_C\widehat E_C^\dagger&=J.
\end{aligned}
```

The [intrinsic active-span identity](ONE_CLEAN_MIDPOINT_RIGIDITY.md#the-essential-axes-are-intrinsic)
makes the displacement image of J exactly
$`L=\langle z_c,z_a\rangle`$. Independence of the output native
list therefore forces exactly two active axes spanning L. The full
Pauli support of J is $`p+L`$; the output expansion has support
$`\widehat p_C+L`$, so $`\widehat p_C\in p+L`$. Since L is
isotropic, only $`z_c,z_a`$ among its three nonzero labels anticommute
with $`\widehat p_C`$. The actual active output axes are thus precisely
$`Z_c,Z_a`$, up to their retained Hermitian signs.

Normalize those actual signs and orientations only. The active output
band has the literal form

```math
\widehat E_{C,\mathrm{act}}=A C_{\mathrm{act}},
\qquad
C_{\mathrm{act}}\widehat P_C C_{\mathrm{act}}^\dagger=P,
```

where $`C_{\mathrm{act}}`$ is an actual diagonal Clifford on c,a,
including its scalar: each orientation differs from A by an even T
power. Every silent output axis commutes $`Z_c,Z_a,\widehat P_C`$,
hence also P and $`C_{\mathrm{act}}`$. Put

```math
\begin{aligned}
W_C&=K^\dagger\widehat E_{C,\mathrm{sil}}K,\\
D_C&=K^\dagger C_{\mathrm{act}}\widehat K_C,\\
U_C&=K^\dagger A^\dagger\widehat E_C\widehat K_C
    =W_CD_C.
\end{aligned}
```

Both factors preserve physical $`Z_c`$ individually. The silent axes
remain independent on the positive clean code: a relation created by
restriction would put P in their physical span, but that span commutes
$`Z_c`$ while P anticommutes with it. Thus the occupied $`W_C`$ is
one actual native band, and the occupied $`D_C`$ is a logical Clifford.
The entire occupied relative encoder $`U_C`$ belongs to
$`\mathcal C_3`$, without selecting a convenient unused-sector extension.

On the other hand, the fixed input and middle word have exact occupied
relative action

```math
\begin{aligned}
U_{\mathrm{in}}
 &=\left.K^\dagger A^\dagger E_BAWK\right|_{c=0}\\
 &=T_B T_{R_lB}T_B^\dagger W_0\\
 &=T_{R_lB}W_0.
\end{aligned}
```

The existing [two-native-band theorem](NATIVE_BAND_SOURCE_STRUCTURE.md#2-two-native-bands-exact-c3-classification)
says this is in $`\mathcal C_3`$ exactly when every projected W_0 axis commutes
with $`R_lB`$. Since each projected axis is $`R_l^{e_i}D_i`$,
this is exactly commutation of D_i with B. An anticommuting axis therefore makes
$`U_{\mathrm{in}}\notin\mathcal C_3`$. An occupied Clifford target
G would require $`U_{\mathrm{in}}=(U_C|_{c=0})G`$, contradicting
invariance of $`\mathcal C_3`$ under right Clifford multiplication.
This proves the obstruction.

If all D_i commute B, the projected native probe family fixes B's binary
label, so this mechanism supplies no growing orbit. The explicit formal
completion $`A T_{R_pB}W`$ pays the partner but requires two output
bands whenever an axis moves B. The theorem above also excludes every
different one-band completion of the same fixed input and middle word.

### 3.1. Complete clean removal and diagonal closure for this family

The clean comparison and the factorization $`U_C=W_CD_C`$ above use
only full initialized-unitary return, without assuming that its occupied
target G is Clifford. Write $`W_{C,0},D_{C,0}`$ for the positive-clean
restrictions of those individually clean-preserving factors. On the
entire logical register the exact identity is

```math
\begin{aligned}
G&=D_{C,0}^\dagger W_{C,0}^\dagger T_{R_lB}W_0,\\
R_l&=-Z_a.
\end{aligned}
```

This is an actual zero-clean three-band word in the original logical
coordinates. The projected input band $`W_0`$ and output band $`W_{C,0}`$
are independent native lists, $`R_lB`$ is one nonidentity Hermitian
Pauli, and $`D_{C,0}`$ is an actual logical Clifford. Every literal
scalar is retained. No unused-sector extension has been substituted.

If G is diagonal, the
[native-remainder closure](THREE_LAYER_NATIVE_STRUCTURE.md#5-peeling-the-pairs-and-closing-the-native-remainder)
applies with an empty pair list: the one-axis logical middle band has
only singleton orbits under the simultaneous diagonal comparisons.
That closing argument proves that G is diagonal $`\mathcal C_3`$,
including its possibly complex phases. In particular,

```math
G=e^{i\theta}\,\mathrm{diag}((-1)^{g(x)})
\quad\Longrightarrow\quad
\deg_{\mathbb F_2}g\le3.
```

For a literal occupied XOR query $`Q_f`$, Hadamards on its output turn
this same complete logical word into the sign $`(-1)^{zf(x)}`$.
They change only its outer endpoint Cliffords; the logical middle band
still has one axis. Thus

```math
G=Q_f\quad\Longrightarrow\quad
\deg_{\mathbb F_2}f\le2.
```

For w input-W axes and v silent output axes, the rewrite uses
$`w+v+1`$ T/T-dagger gates in at most three bands on the original
logical wires, with no clean allocation. Fixed native coordinate and
logical Clifford synthesis costs $`O(q^2)`$ elementary Cliffords on q
logical wires. Its actual reversed inverse has the same bounds.
The full logical operator identity preserves arbitrary work and reference
inputs and the target's complete return contract, with zero new error.

This closes the diagonal-target degree question for this specified
physical family, including non-Clifford Boolean signs. Changing the
fixed active input or middle word, or adding original pairs leaves this
corollary's scope. All admissible core-coupled silent input axes are
already covered. The general
inside-span matching step and unrestricted depth problem remain open.

## 4. Native ledger and exact fixtures

For three wires $`c=0,a=1,B=X_2`$, the following coordinate words are
chronological:

```text
K   = [H0, CX(0,1)]
U_B = [CX(2,0), CX(2,1), H1, S1, S1, H1,
       CX(0,1), H0, CX(0,1), H2]
U_C = [CX(0,2), CX(1,2), H2, S2, S2]
```

Conjugation by $`U_B`$ maps $`Z_0,Z_1,Z_2`$ to $`PB,QB,B`$;
conjugation by $`U_C`$ maps them to $`Z_0,Z_1,R_pB`$. The
chronological band words are therefore

```text
E_B = U_B†, [T0,T1,T2†], U_B
E_C = U_C†, [T0,T1,T2], U_C
V   = K, [T0,T1], E_B, E_C†, K†
```

The complete source has eight T/T-dagger gates in three bands of sizes
two, three, three. Its displayed coordinate word has fourteen CNOTs,
twelve H gates, four S gates and four S-dagger gates. Expanding every
S-dagger into three S gates gives at most forty-two elementary H/S/CNOT
gates. The actual reversed inverse has the same counts. Only c is
initialized; all remaining wires and arbitrary references are returned
on the occupied code, with zero error and leakage. For B on a larger
register, a fixed Pauli-coordinate Clifford and its actual inverse cost
$`O(q^2)`$ elementary Cliffords on the existing q wires; the same
eight T gates and three bands suffice.

The [three bounded exact fixtures](../../tests/test_minimal_physical_pair_completion.py)
check the three-wire native coordinates, signed pair, genuine singleton,
every physical column, full unused sector, actual inverse and occupied
factorization. A two-logical-wire extension with $`W=T_{Z_t}`$ and
$`B=X_t`$ checks a specific non-Clifford Pauli image of
$`T_{-Z_aB}W`$, protecting the unpaid partner. These certify the stated
mechanism; they do not enumerate candidate completions or replace the
analytic arbitrary-output proof. The added four-wire fixture uses a
core-coupled input axis $`R_pZ_t`$ and an extra positive middle
singleton $`R_pZ_s`$. It retains all sixteen columns and the actual
inverse, and checks the two pair-parity blocks, including the literal
omega scalar and Clifford $`S_{X_t}`$. Its occupied target is not
asserted Boolean; it checks the sector-reduction mechanism used below.


## 5. All paired sources with two clean-active axes at both endpoints

The fixed cancellation in Sections 1--3 is not needed for the following
broader degree theorem. Consider any actual one-clean three-band word

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A
```

with complete occupied Boolean sign action, arbitrary logical/dirty
inputs and references, and exact return of its clean qubit. Assume that
both intrinsic endpoint clean-active ranks are two and that the original
physical middle list has a genuine paired orbit. Then

```math
\deg_{\mathbb F_2}g\le6,\qquad
\deg_{\mathbb F_2}f\le5
\quad\text{for an occupied Boolean XOR query }Q_f.
```

This permits arbitrary actual middle singletons, arbitrary silent
endpoint bands including core-coupled axes, arbitrary endpoint
Cliffords, and arbitrary unused-sector action. It is uniform in width
and the sizes of those silent lists. It is not a bound for growing
endpoint active ranks or a new unrestricted depth obstruction. The
no-original-pair branch retains its separate existing theorem; it is
not assigned the bound six by this argument.

### 5.1. Exactly one pair, and no negative singleton after orientation repair

Write a_1,a_2 for the independent active A-axis labels, p_A for the
input clean Pauli label, and J for its actual one-band image. Then

```math
N_A b=\langle a_1,b\rangle a_1+
       \langle a_2,b\rangle a_2.
```

For a nonsingleton orbit b,d, native matching gives
$`d=b+N_A b`$. Since the actual middle axes b,d commute,
$`\langle b,N_A b\rangle=0`$. Independence of a_1,a_2 makes
the only nonzero possibility

```math
\langle a_1,b\rangle=\langle a_2,b\rangle=1,
\qquad r=b+d=a_1+a_2.
```

Two different paired orbits would have this same difference, giving a
relation among four distinct physical middle axes. The independent
middle list therefore has exactly one pair.

Normalize its actual orientations by the literal Clifford correction of
[the original-pair factorization](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md#1-keep-the-actual-middle-axes-and-their-signs).
Thus

```math
E_B=C F E_-E_+,
```

where F is the reoriented paired phase commuting J, E_- consists of
negative singleton axes, and E_+ consists of positive singleton axes.
Absorb the actual Clifford C into the output encoder and conjugate its
native list. This preserves its intrinsic active rank and every literal
phase. The corrected clean comparison is

```math
S'=E_-JE_-^\dagger.
```

Every singleton commutes both active A axes: its fixed binary label
lies in the kernel of N_A. Its sign under J is therefore its sign
under P_A. In particular each negative singleton anticommutes P_A.
The span of the negative list intersects $`L_A=\langle a_1,a_2\rangle`$
trivially. Indeed every such combination commutes the paired axis b,
whereas a_1 and a_2 separately anticommute with b. The only possible
nonzero intersection is r. But r is the sum of the two original pair
labels, so its presence would violate independence of the complete
middle list.

Consequently the union of the two active A axes and the negative
singleton list is independent and commuting, and every member is active
for P_A. It is an actual one-band representation of S', giving

```math
\mathrm{rank}(S')=2+|B_-|.
```

Here rank denotes intrinsic symplectic displacement rank. The corrected
output image has rank two, so B_- is empty. Thus S'=J and both active
endpoint spans are the same L_A. This excludes a hidden cancellation
between negative singleton directions and the active core.

### 5.2. A common computational parity before sector restriction

Let R be the signed Hermitian product of the actual paired axes.
Its label r belongs to both endpoint active spans and to the original
middle span. It commutes every axis in all three bands and both clean
images. The original-pair transport is a full physical identity,

```math
V(K_A^\dagger R K_A)V^\dagger=K_C^\dagger R K_C.
```

Its occupied restrictions are nontrivial logical Paulis. The
[paid pair-center normalization](OCCUPIED_LOGICAL_TRANSPORT.md#3-a-paid-common-normalization-of-all-pair-differences)
uses this identity and the Boolean target to extract a quadratic sign,
make an invertible affine change of active variables, and rotate only
returned flat directions. It places the input and output restrictions
of R at the same computational Z parity. The original Boolean degree
is at most the maximum of two and the normalized sign degree.

These are actual endpoint Clifford changes. They preserve every
interior native axis list, its active rank and the original pair.
A possible physical clean-Z discrepancy between the endpoint parities
is corrected by the retained clean-controlled Pauli Clifford before
restriction. There is then one common physical computational parity
for R at both ends. Its occupied restriction is chosen as coordinate z
using only an affine change on the remaining sign variables and returned
flat wires. The physical parity may still include the clean Z factor,
which is fixed to +1 on the occupied input. No independently selected
coordinates in its two sectors have
been substituted for this global normalization.

### 5.3. Each parity sector is a complete zero-clean three-band word

Fix either eigenvalue of R. The two active endpoint axes have product
R up to their retained sign, so their images in this sector coincide
up to sign. Their two odd T powers therefore combine to an even power
and a literal scalar: an actual Clifford. The same is true of the two
axes of the original middle pair. All these scalars are retained.
Their Clifford descriptions may differ between the two parity sectors.

The surviving silent endpoint lists remain independent after quotienting
by r: a relation would put r in the span of a silent list, whereas r
belongs to its complementary active span in one independent physical
endpoint list. The positive singleton list likewise stays independent
modulo r. Otherwise a combination of its labels would equal r, the
product of the removed original pair, contradicting independence of the
complete physical middle list.

The actual clean image J in this parity sector is now a signed Pauli
J_z, because the active endpoint band has become Clifford. Every
surviving axis commutes J_z. No second relation is introduced by the
clean-code quotient. An active A axis a_1 commutes every surviving
list and with R, but anticommutes J and hence J_z. Thus J_z cannot
belong to the span of any surviving list after the R quotient.
This is the independence witness; exact clean return alone would not
have supplied it.

The pair-sector Clifford preserves J_z and commutes the positive
singletons. Move it to an endpoint by conjugating the neighboring
native list by that actual Clifford. Conjugation preserves the proved
independence and transports its clean-code witness. Absorb the active
endpoint Cliffords into the encoders in the same way. Fixed Clifford
coordinates for J_z then remove the sole clean qubit. The restricted
endpoint isometries are actual logical Cliffords, and the three
surviving lists are independent commuting logical Pauli lists.

Therefore, for each fixed parity bit z, the entire remaining logical
operator is a zero-clean circuit with at most three actual native bands.
This includes every original dirty input, not just an initialized column.
Its target is the normalized diagonal Boolean sign at that z. The
[zero-clean three-layer theorem](THREE_LAYER_NATIVE_STRUCTURE.md)
gives a Boolean polynomial g_z of degree at most five in those common
remaining computational inputs.

### 5.4. Restore the parity bit, the quadratic phase and the contract

In the common computational coordinates from Section 5.2, the normalized
sign polynomial is exactly

```math
h(z,y)=(1+z)g_0(y)+z g_1(y)
\quad\text{over }\mathbb F_2.
```

Thus h has degree at most six. Returning through affine coordinates
preserves that bound; restoring the extracted quadratic sign does not
increase it. This proves the sign theorem. For an XOR query, first
conjugate its output by Hadamards. This changes only the endpoint
Cliffords, preserves both clean-active ranks and the original middle
pair, and gives the sign $`(-1)^{t f(x)}`$. Its degree is one plus
that of a nonzero f, proving the bound five.

The two sector reductions are a degree-proof device, not a free coherent
choice between two compiled circuits. The implementation remains the
original physical word and its actual native inverse: one clean qubit,
its original dirty allocation, at most three T bands and its original
T/Clifford counts. The global normalization, if emitted, uses
$`O(Q^2)`$ elementary Cliffords on Q existing wires in each direction.
For a fixed parity sector, the comparison word removes the four active
endpoint T factors and the two paired factors, retaining no more than
$`|A_s|+|B_+|+|C_s|`$ T gates. Its coordinate Cliffords and their
actual inverses also cost $`O(Q^2)`$. No such sector compiler is
credited as a uniform upper-bound improvement.

All identities retain literal scalar phases, arbitrary references,
complete logical columns, and the exact clean/dirty return contract.
The [retained approximation-to-exactness theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
also gives its conservative floor $`2^{-530690}`$ for higher-degree
Boolean targets in this class, allowing one common phase. It does not
exclude arbitrary fixed tolerances or establish dimension-growing depth.
