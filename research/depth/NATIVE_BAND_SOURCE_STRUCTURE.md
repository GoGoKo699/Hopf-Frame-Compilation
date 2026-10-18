# Native band structure of third-level sources

[Native Clifford transport](NATIVE_CLIFFORD_TRANSPORT.md) ·
[One-band matching](ONE_BAND_AXIS_MATCHING.md) ·
[Six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md)

The physical independence of a T band's axes supplies source restrictions
that third-level membership alone does not. Two native bands are third-level
exactly when their axis lists commute. A native band added on either side
of a diagonal third-level source cannot leave that source class while
remaining third-level. For three bands, the overlap of the first two axis
spans limits how many middle axes can fail to commute with the outer band.
In particular, zero overlap certifies the source hypothesis of the retained
six-layer degree-78/26 theorem. A further commutation-cover certificate
allows arbitrarily many bad zero-character middle axes: if each middle
axis commutes with the entire first list or the entire last list, an exact
C3 source is diagonal C3 in fixed Clifford coordinates.

These are properties of the **specified native factorization**, not
invariants of the implemented unitary. They do not infer third-level
membership from complete workspace return, prove a general three-band
classification, or improve a complete-frame compiler. No intermediate
Clifford membership is assumed.

## 1. Common native coordinates and exact hypotheses

Absorb fixed Clifford intervals into fixed axis coordinates. Up to fixed
Cliffords on the left and right, a specified three-band source is

```math
F=E_CE_BE_A.
```

Each E is a Clifford conjugate of one actual physical T/T-dagger layer.
Its active binary Pauli axes form an independent isotropic set, respectively
A, B and C. Write $`U_A,U_B`$ for the first two binary spans. Inactive wires
contribute no axis. Axis signs, T orientations and native scalar phases
remain in the actual operators, although binary-axis arguments do not
depend on them. Third-level membership is unchanged by fixed Clifford
factors on either side.

A source U is third-level, or C3, when $`U\tau(p)U^\dagger`$ is Clifford
for every Hermitian Pauli $`\tau(p)`$. We call U **diagonal C3 in fixed
Clifford coordinates** when

```math
U=C_{\rm out}D C_{\rm in},\qquad
D|z\rangle=\omega^{q(z)}|z\rangle,\qquad \omega=e^{i\pi/4},
```

up to a fixed retained scalar, where

```math
q(z)=c+\sum_i a_i z_i
       +2\sum_{i\lt j}b_{ij}z_i z_j
       +4\sum_{i\lt j\lt k}c_{ijk}z_i z_j z_k\pmod8.
```

This is precisely the source class of the
[signed transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#2-the-diagonal-source-and-its-pauli-directions).
The two fixed end Cliffords need not be inverses.

We use the established [one-band matching theorem](ONE_BAND_AXIS_MATCHING.md):
for a native band E and a Clifford J,
$`EJE^\dagger`$ is projectively Clifford if and only if J permutes the
individual active binary axes of E. Its robust direction has projective
operator radius $`\sqrt2/4`$. Preserving only the axis span is insufficient.

## 2. Two native bands: exact C3 classification

**Two-band theorem.** The source $`E_BE_A`$ is C3 if and only if every
axis in A commutes with every axis in B. Consequently every two-physical-band
C3 source is diagonal C3 in fixed Clifford coordinates.

For every Pauli P, $`J_P=E_APE_A^\dagger`$ is Clifford. If the source is
C3, one-band matching says that its symplectic action permutes the B-axis
set. Independence of A permits a Pauli P whose commutation vector against
A is any specified unit vector. Isolating $`a\in A`$ gives the single
transvection

```math
T_a(v)=v+\langle a,v\rangle a.
```

A transvection cannot nontrivially permute an isotropic set. If it sends b
to the distinct member $`b+a`$, then
$`\langle b,b+a\rangle=\langle b,a\rangle=1`$, contradicting isotropy.
It therefore fixes every B axis, proving all cross-axis commutations.

Conversely the union of A and B is isotropic. A fixed Clifford makes all
its axes diagonal Paulis. Dependencies and signs between the two lists are
allowed. Their native phase becomes an integer sum of affine binary
parities modulo eight, which has the weighted cubic form above by the
integer parity expansion. Thus the source is diagonal C3. This is an
analytical identification: a dependent union is not declared to be one
physical T band.

## 3. Native right extension

**Half-angle lemma.** Let H be a traceless Hermitian involution. If

```math
R=\frac{I+iH}{\sqrt2}
```

is projectively Clifford, then H is a signed Hermitian Pauli.

A Clifford has a flat Pauli expansion. For completeness, if K has
symplectic matrix S, normalized Clifford trace orthogonality gives

```math
|\mathrm{tr}(P K)|\in\{0,2^{-\mathrm{rank}(S-I)/2}\}
```

for every Pauli P. Indeed the trace of conjugation by PK on the Pauli
basis is the sum of a sign character over $`\ker(S-I)`$; it is either
zero or the size of that kernel. Parseval therefore gives exactly
$`2^{\mathrm{rank}(S-I)}`$ nonzero Pauli coefficients, all of the displayed
magnitude. Since H is traceless, R has identity coefficient exactly
$`1/\sqrt2`$. Its flat expansion consequently has exactly two terms.
The other is the sole nonzero Pauli coefficient of H, which must be a
signed Hermitian Pauli because $`H=H^\dagger`$ and $`H^2=I`$.
A global scalar does not change this argument.

**Right-extension lemma.** Suppose U is C3, $`E_A`$ is a native band,
and $`UE_A`$ is C3. Then $`UaU^\dagger`$ is a signed Pauli for every
actual Hermitian active axis a of $`E_A`$.

Choose P isolating a among the A axes. With the actual scalar and rotation
orientation retained, its one-band image is of the form

```math
E_APE_A^\dagger=P\exp(\pm i\pi a/4).
```

Both its U-conjugate and $`UPU^\dagger`$ are Clifford. Multiplication by
the latter's actual inverse shows that
$`\exp(\pm i\pi UaU^\dagger/4)`$ is projectively Clifford. The axis
image is a traceless Hermitian involution, so the half-angle lemma applies.

**Diagonal right-extension corollary.** If U is diagonal C3 in fixed
Clifford coordinates, $`E_A`$ is native and $`UE_A`$ is C3, then
$`UE_A`$ is diagonal C3 in fixed Clifford coordinates.

Remove U's fixed end Cliffords, conjugating the native band into its input
coordinates, so that U is diagonal. Define its Pauli translation domain

```math
\mathcal R=\{r:UX^rU^\dagger\text{ is Pauli up to scalar}\}.
```

The [phase-sensitive diagonal Clifford extension](NATIVE_CLIFFORD_TRANSPORT.md#2-the-diagonal-source-and-its-pauli-directions)
gives a fixed diagonal Clifford G such that $`H=G^\dagger U`$ commutes
literally with every $`X^r`$, $`r\in\mathcal R`$, and every Z Pauli.
Its weighted cubic phase depends only on quotient Z observables with labels
in $`\mathcal R^\perp`$.

The right-extension lemma places the X label of every A axis in
$`\mathcal R`$. Hence each A axis commutes with H and with every quotient
Z observable. These observables and the A axes have a joint isotropic span.
A fixed Clifford simultaneously diagonalizes them. In those coordinates H
remains weighted cubic, since affine parity substitution preserves that
form modulo eight, and $`E_A`$ contributes ordinary T parity phases.
Thus $`HE_A`$, and consequently $`UE_A=GHE_A`$, has the stated form.
All signed corrections are retained.

### Certified native stripping

**Stripping lemma.** Suppose W is C3, $`E_L`$ is a commuting native
phase band, and $`WPW^\dagger`$ is a signed Hermitian Pauli for every
actual active axis P of $`E_L`$. Then $`WE_L^\dagger`$ is C3.

Use the actual signed-axis convention

```math
N_e(P)=e^{ie\pi/8}e^{-ie\pi P/8},\qquad e\in\{1,-1\}.
```

For any Hermitian Pauli p the exact inverse-band derivative is

```math
E_L^\dagger pE_Lp
=\prod_{P\in L:\langle P,p\rangle=1}e^{ie_P\pi P/4}.
```

The native scalars cancel in this derivative, while remaining in the
actual band and its inverse. Conjugation by W turns the displayed product
into Pauli quarter-turns, hence a Clifford. Therefore

```math
(WE_L^\dagger)p(WE_L^\dagger)^\dagger
=W(E_L^\dagger pE_Lp)W^\dagger\,WpW^\dagger
```

is Clifford for every p. This proves the lemma without assuming
$`W^\dagger`$ is C3. Independence is not needed for this implication;
the use below strips an actual independent subset of A. A native band
is removed only after its complete signed Pauli images have been proved.

## 4. Native left extension without an inverse-C3 premise

**Left-extension lemma.** Suppose U is C3, $`E_C`$ is native and
$`E_CU`$ is C3. Then $`U^\dagger cU`$ is a signed Pauli for every
active Hermitian axis c of $`E_C`$. No assumption that $`U^\dagger`$
is C3 is made.

For each Hermitian Pauli P, the operator $`J_P=UPU^\dagger`$ is Clifford.
One-band matching says that its adjoint action permutes the individual
Hermitian C axes up to signs. It therefore preserves their real linear span
V. These restricted adjoint actions are real orthogonal involutions and
commute exactly: $`J_PJ_Q`$ and $`J_QJ_P`$ differ only by the original
Pauli scalar. They have a simultaneous real orthonormal eigenbasis.

For a nonzero common eigenvector v, the pullback $`U^\dagger vU`$ is a
simultaneous eigenoperator of conjugation by every Pauli. Each such joint
eigenspace in the full matrix algebra is one-dimensional and spanned by a
Pauli. Since v is Hermitian, its pullback is a real scalar times a Hermitian
Pauli. Therefore $`v^2`$ is scalar.

Write $`v=\sum_i r_i c_i`$ with real coefficients. The C axes commute and
are independent, so all labels $`c_i+c_j`$ for distinct unordered pairs
are distinct and nonzero. In

```math
v^2=\left(\sum_i r_i^2\right)I+2\sum_{i\lt j}r_i r_j c_ic_j
```

no nonidentity terms can cancel. Thus $`r_ir_j=0`$ for every pair: each
common eigenvector is proportional to one individual axis. A simultaneous
eigenbasis consists of these coordinate vectors, so every individual C axis
is a common eigenvector. Pulling it back proves the lemma.

The pulled-back axes remain independent and isotropic, because conjugation
preserves products and commutation. With their actual signs, the original
T orientations and the retained native scalar, they define a native-axis
band $`E_{C'}`$ satisfying the literal identity

```math
E_CU=UE_{C'}.
```

**Diagonal left-extension corollary.** If U is diagonal C3 in fixed
Clifford coordinates and $`E_CU`$ is C3, then $`E_CU`$ has the same
source form. Apply the diagonal right-extension corollary to $`UE_{C'}`$.
This identity does not supply a free physical coordinate change or source.

## 5. Three-band overlap restricts noncommuting middle axes

**Overlap theorem.** Suppose $`F=E_CE_BE_A`$ is C3, and put

```math
t=\dim(U_A\cap U_B).
```

At most 2t individual active B axes fail to commute with some active C axis.

Take Pauli labels $`p\in U_A^\perp`$, so that the first native band
commutes with them. Their commutation vectors against B range over

```math
\mathcal U=\{(\langle b_i,p\rangle)_i:p\in U_A^\perp\}
\subseteq\mathbb F_2^{|B|}.
```

A coefficient vector r lies in $`\mathcal U^\perp`$ exactly when
$`\sum_i r_i b_i\in U_A`$. Independence of B gives
$`\dim\mathcal U^\perp=t`$. For every $`u\in\mathcal U`$, the middle
Pauli image has symplectic action

```math
S_u=I+N_u,\qquad
N_u(v)=\sum_i u_i\langle b_i,v\rangle b_i.
```

C3 membership and outer one-band matching force $`S_u`$ to permute the
individual C axes. Isotropy of B gives $`N_uN_v=0`$, so
$`S_uS_v=S_{u+v}`$.

Every orbit of this elementary abelian action on the independent C-axis set
has size at most two. Otherwise choose u,v and c for which
$`c,S_uc,S_vc,S_{u+v}c`$ are four distinct orbit members. Their binary
sum is zero, contradicting independence of C.

Fix $`c\in C`$. For a singleton orbit, independence of B implies that
every coordinate functional $`u\mapsto u_i`$ with
$`\langle b_i,c\rangle=1`$ vanishes on $`\mathcal U`$. For a two-element
orbit $`\{c,c'\}`$, let $`\chi:\mathcal U\to\mathbb F_2`$ be its
nonzero swap character. Then

```math
\sum_{i:\langle b_i,c\rangle=1}u_i b_i
=\chi(u)(c+c').
```

Independence of B says that each coordinate functional occurring on the
left is either zero or exactly $`\chi`$. At least two nonzero ones must
occur: with precisely one, $`c+c'=b_i`$ and
$`\langle c,c'\rangle=\langle c,b_i\rangle=1`$, contrary to isotropy
of C.

Partition all B indices by equality of their coordinate functionals on
$`\mathcal U`$. Let z be the size of the zero class, and
$`m_1,\ldots,m_s`$ the sizes of nonzero classes of size at least two.
Every B axis that anticommutes with a C axis belongs to one of these
classes. The z relations $`e_i`$ for zero coordinates and the
$`m_j-1`$ relations $`e_i+e_{\rm base}`$ within each repeated class
belong to $`\mathcal U^\perp`$. Their disjoint class supports make all
these relations jointly independent. Hence

```math
z+\sum_j(m_j-1)\le t,\qquad
z+\sum_jm_j\le2\left(z+\sum_j(m_j-1)\right)\le2t.
```

This proves the claim. Both native axis-list independence hypotheses are
essential to the argument. The theorem permits arbitrarily many commuting
axes and does not impose a nilpotence condition on comparison groups.

## 6. Source certificates and the six-layer consequence

Each of the following conditions certifies that the specified three-band
C3 source is diagonal C3 in fixed Clifford coordinates:

1. $`U_A\cap U_B=\{0\}`$. The overlap theorem makes B and C commute,
   so the outer two-band source is diagonal C3; use right extension.
2. The specified outer two-band source $`E_CE_B`$ is C3, equivalently
   all B axes commute with all C axes; again use right extension.
3. The specified inner two-band source $`E_BE_A`$ is C3, equivalently
   all A axes commute with all B axes; use left extension.
4. Every B axis commutes with the entire A list or the entire C list;
   use the [commutation-cover theorem](#61-commutation-cover-theorem) below.

No condition makes both adjacent partial sources Clifford. For example,
on one arbitrary work qubit the literal source

```math
F=T_ZT_Z^\dagger T_X=T_X,\qquad T_X=HTH,
```

has zero A/B span overlap although its inner two-band source
$`T_Z^\dagger T_X`$ is outside C3. It uses three T/T-dagger calls,
two H gates and three declared physical slots; the actual inverse has the
same charge. The fair simplified source has one T call, two H gates and
one slot. This cancellation fixture is a source identity, not a returned
query or a resource improvement.

More generally, if every specified partial source $`E_j\cdots E_1`$
remains C3, the complete source is diagonal C3. The first band has this
form, and left extension supplies induction. The same holds for a sequence
built by successive right extensions. No inverse-C3 premise is used.

For the actual [six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md),
any of these certificates therefore verifies the source hypothesis
of the [native transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query).
Under the original full operator error

```math
\eta\lt\sqrt2/8,
```

the Boolean output has degree at most **78**, or **26** when the outer
target shear is visible. This is a new sufficient test for the existing
source class, with its existing error radius; it is not an unrestricted
six-layer degree theorem.

### 6.1 Commutation-cover theorem

**Theorem.** Suppose $`F=E_CE_BE_A`$ is exactly C3 and every middle
axis commutes with the entire A list or with the entire C list. Then F
is diagonal C3 in fixed Clifford coordinates, with no bound on width,
overlap or relation count. In particular this covers every source whose
bad B axes all have zero activation character: such axes lie in $`U_A`$
and therefore commute with A.

Let $`B_0`$ be exactly the B axes that fail to commute with some C axis,
and let $`B_g`$ be the remaining B axes. The hypothesis makes $`B_0`$
commute with all of A; $`B_g`$ commutes with all of C by definition.
Set

```math
U=E_CE_{B_g},\qquad E=E_{B_0}E_A,\qquad F=UE.
```

Every phase in U commutes with every other phase in U, and the same is
true of E. Both are diagonal C3 sources by the parity expansion in
section 2. Relations between their component lists are allowed; neither
merged list is declared to be one physical band. Every A axis commutes
with E, while every $`g\in B_g`$ commutes with U and with $`B_0`$.

Choose an actual A axis $`a_i`$ and a Hermitian Pauli p isolating it
among the entire independent A list. With the signed phase convention
above, the exact left derivative is

```math
K_i=EpE^\dagger p
=e^{-ie_i\pi a_i/4}
 \prod_{b_j\in B_0:\langle p,b_j\rangle=1}
 e^{-i\delta_j\pi b_j/4}.
```

Here $`e_i,\delta_j\in\{1,-1\}`$ are the actual orientations.
Since F and U are C3, their complete Pauli images are Cliffords, and

```math
J_i=(FpF^\dagger)(UpU^\dagger)^\dagger=UK_iU^\dagger
```

is Clifford with its literal scalar retained. Suppose $`a_i`$
anticommutes with some $`g\in B_g`$. All selected $`B_0`$
quarter-turns fix g. Consequently

```math
K_igK_i^\dagger=-ie_i a_i g,\qquad
J_igJ_i^\dagger=-ie_i(Ua_iU^\dagger)g.
```

The second identity uses $`UgU^\dagger=g`$ literally. Its left side
is a signed Pauli, so multiplication by g and the displayed scalar proves
that $`Ua_iU^\dagger`$ is a signed Hermitian Pauli. Because E fixes
$`a_i`$, this is also $`Fa_iF^\dagger`$.

Let $`A_1`$ consist of the actual A axes that anticommute with at least
one $`B_g`$ axis, and put $`A_0=A\setminus A_1`$. The preceding
argument certifies every $`A_1`$ image as Pauli. The stripping lemma
therefore proves that

```math
V=FE_{A_1}^\dagger=E_CE_BE_{A_0}
```

is C3. This cancellation uses the actual inverse, including the native
scalar. Every $`A_0`$ axis commutes with $`B_g`$ by definition and
with $`B_0`$ by hypothesis. Thus $`D=E_BE_{A_0}`$ is a commuting
diagonal-C3 source. The outer $`E_C`$ remains an independent native
band. The diagonal left-extension corollary applied to $`V=E_CD`$
makes V diagonal C3. Finally $`F=VE_{A_1}`$ is C3 and $`E_{A_1}`$
is an actual native subset, so the diagonal right-extension corollary
makes F diagonal C3 too.

No original adjacent prefix is assumed C3, and no inverse-C3 hypothesis
enters. This is an exact source certificate. It does not extend the
all-Pauli robust source promises in section 7, nor infer source membership
from return. Its six-layer consequence is the existing 78/26 bound above;
it changes no physical gate, coordinate, inverse or workspace charge.

### 6.2 Removing good middle axes in a general C3 source

The certified stripping argument also gives a reduction without the
commutation-cover hypothesis. Let $`B_g`$ contain every middle axis
commuting with the entire C list, and put $`B_0=B\setminus B_g`$.
Let $`A_1`$ contain every A axis that anticommutes with some $`B_g`$
axis, and put $`A_0=A\setminus A_1`$. For any exact native three-band
C3 source F, define the actual reduced source

```math
W=E_CE_{B_0}E_{A_0}.
```

**Good-middle elimination.** W is C3, and

```math
F=E_{B_g}WE_{A_1}
```

literally. Moreover F is diagonal C3 in fixed Clifford coordinates if
and only if W is. Every remaining middle axis is bad relative to the
unchanged C list. This is an analytical reduction; the removed rotations
remain in the displayed factorization and in every original resource ledger.

Again set $`U=E_CE_{B_g}`$ and $`E=E_{B_0}E_A`$. U is diagonal
C3, but E is now allowed to lie outside C3. For p isolating $`a_i`$,
write

```math
R_i=E_ApE_A^\dagger p=e^{-ie_i\pi a_i/4},\qquad
K_B=E_{B_0}pE_{B_0}^\dagger p.
```

The exact derivative and the Clifford comparison are

```math
K_i=EpE^\dagger p
   =(E_{B_0}R_iE_{B_0}^\dagger)K_B,\qquad
J_i=UK_iU^\dagger=(FpF^\dagger)(UpU^\dagger)^\dagger.
```

Every $`g\in B_g`$ commutes with $`E_{B_0}`$ and is fixed by
$`K_B`$ and U. If $`a_i`$ anticommutes with g, these exact identities
give

```math
J_igJ_i^\dagger=-ie_i(Fa_iF^\dagger)g.
```

Since $`J_i`$ is Clifford, $`Fa_iF^\dagger`$ is a signed Pauli.
Thus every $`A_1`$ axis has a certified image even though E need not
be C3. Stripping $`E_{A_1}`$ gives the C3 source

```math
V=FE_{A_1}^\dagger=E_CE_BE_{A_0}.
```

The source V fixes every $`B_g`$ axis literally: these axes commute
with C, B and $`A_0`$. A second application of the stripping lemma
therefore makes

```math
VE_{B_g}^\dagger=E_CE_{B_0}E_{A_0}=W
```

C3. Commuting $`E_{B_g}`$ through the unchanged outer and reduced
inner factors proves the stated literal factorization of F.

If W is diagonal C3, the native left-extension corollary gives
$`V=E_{B_g}W`$ diagonal C3, and native right extension restores
$`F=VE_{A_1}`$. Conversely, if F is diagonal C3, native right
extension by the actual inverse band $`E_{A_1}^\dagger`$ gives V
diagonal C3 because V is already proved C3; extension by
$`E_{B_g}^\dagger`$ then gives W diagonal C3. This proves the
equivalence.

The companion [domain-axis phase closure](CONTROLLED_PAULI_TRANSPORT.md#closure-under-commuting-phase-gadgets-on-domain-axes)
also makes scalar-Pauli block form in fixed input/output Clifford
coordinates equivalent for F and W. Apply it in the two actual stages:
$`W\mapsto V=E_{B_g}W`$ uses the fixed $`B_g`$ axes in W's exact
output domain; $`V\mapsto F=VE_{A_1}`$ uses the certified $`A_1`$
axes in V's exact input domain. Their Pauli images persist under the
respective phase factors, giving the converse too. No assertion that
$`A_1`$ lies in W's input domain is needed. The companion proof retains
the ordered phases and both endpoint coordinates.

The reduction can remove axes from the first span, so activation
characters must be recomputed for $`A_0,B_0`$. It does not remove a
mixed core in which every remaining B axis interacts with C and some
also interact with $`A_0`$. In particular it supplies no argument that
a bad singleton C factor in such a core has Clifford sector restrictions.

## 7. Robust source certification requires all Pauli probes

The exact C3 hypotheses in the two-band theorem, overlap theorem and
zero-overlap certificate can be weakened to the following promise for the
specified source F:

```math
\text{for every }p,\qquad
\inf_{K\ {\rm Clifford}}d_{\rm pr}(F\tau(p)F^\dagger,K)
\lt\sqrt2/4.
```

For two bands, robust one-band matching on isolated-axis Pauli probes gives
the same exact cross-axis commutations. For the overlap theorem, the probes
centralizing A give exact outer-axis permutations by robust matching, so
the entire orbit and relation-counting proof is unchanged.

For three bands with zero overlap, this first makes B and C commute
exactly, giving an exact diagonal-C3 source $`U=E_CE_B`$. An isolated-A
probe writes its complete image as the product of the exact Clifford
$`UPU^\dagger`$ and $`UR_aU^\dagger`$, where
$`R_a=\exp(\pm i\pi a/4)`$ is Clifford. The promise puts the second
factor projectively within $`\sqrt2/4`$ of a Clifford. The
[third-level rounding lemma](NATIVE_CLIFFORD_TRANSPORT.md#1-a-robust-exactness-lemma-for-every-third-level-source)
for U makes this factor exactly projectively Clifford. Half-angle rigidity
then places $`UaU^\dagger`$ in the Paulis. The diagonal right-extension
proof applies directly, so F is exactly diagonal C3 in fixed coordinates.

The left-extension statement also has a robust version with U exactly C3
and the displayed promise for $`E_CU`$: robust matching gives exact axis
permutations for every $`UPU^\dagger`$, and its proof then applies
unchanged. The matching converse also makes every complete Pauli image
exactly Clifford. Thus the partial-source induction remains valid if every
partial source has the displayed robust promise. At each step its predecessor
is already exactly C3, so robust matching rounds the next step before the
left-extension corollary is used.

These are uniform **all-source-Pauli** promises. The actual returned-query
comparison only controls a particular native Clifford family $`L_x`$ and
does not supply them. No general source rounding is inferred from full return.

## 8. Dependence of axes defeats the left-extension inference

Independence cannot be replaced by an arbitrary commuting phase-polynomial
axis set. On three arbitrary work qubits a,b,t, let U be the literal
Toffoli $`t\mapsto t+ab`$, and take the four commuting axes

```math
c_{u,v}=Z_tZ_a^uZ_b^v,\qquad u,v\in\mathbb F_2.
```

They are dependent: their product is I. Let E apply a positive T phase on
each of these parities. This is a four-gadget diagonal C3 source, **not one
physical native band**.

For a Pauli P with X label $`(\alpha,\beta,\gamma)`$, direct conjugation
by the Clifford $`J_P=UPU^\dagger`$ gives

```math
J_Pc_{u,v}J_P^\dagger
=(-1)^{\gamma+u\alpha+v\beta+\alpha\beta}
 c_{u+\beta,v+\alpha}.
```

Thus the full Pauli-image frame permutes these four axes with signs.
Collecting their T rotations after this permutation changes orientations
only by Clifford rotations, so $`EJ_PE^\dagger`$ is Clifford for every P.
Both U and EU are therefore C3. Nevertheless

```math
U^\dagger Z_tU
=Z_t\,\mathrm{CZ}_{a,b}
=\frac12(c_{0,0}+c_{1,0}+c_{0,1}-c_{1,1})
```

is not Pauli. It is a Hermitian involution with four nonzero coefficients.
Dependence allows the pair-product cancellations excluded in the native
left-extension proof. This is a counterexample to removing that hypothesis,
not to any physical-band theorem.

The source has a fully charged native implementation on the same three
arbitrary wires, with no clean qubit or helper. The retained seven-T CCZ
word and two target Hadamards give U seven T calls, twelve elementary
Cliffords and at most five T layers. E uses one T on t and three parity
words requiring respectively two, two and four CNOTs. Thus E costs four
T calls, eight CNOTs and at most four layers. EU and its actual reversed
inverse each cost eleven T calls, twenty elementary Cliffords and at most
nine layers in this displayed schedule. These are conservative source
ledgers, not optimality claims. Neither source is offered as a returned
XOR query. Every basis column, literal phase and actual inverse is retained.

## 9. Remaining gap and operator meaning

The subsequent [overlap-one theorem](ONE_OVERLAP_SOURCE_STRUCTURE.md)
closes t equal to one, certifying the diagonal C3 class and its existing
78/26 endpoint-degree bounds. It also proves a separate diagonal-source
extension lemma for a commuting T list with exactly one binary dependency;
that conclusion does not restore individual-axis preservation in Section 8.
The later [overlap-two theorem](TWO_OVERLAP_SOURCE_STRUCTURE.md) closes
all t at most two. The [isolated-core extension](ISOLATED_NATIVE_CORE_STRUCTURE.md)
then removes the total-overlap bound for arbitrarily many isolated paired
cores. Its later coupled-incidence theorem covers every exact source in
which zero-character B axes commute with C, forcing all interacting
nonzero classes to be isolated pairs. A separate
[hybrid theorem](HYBRID_DIAGONAL_CORE_TRANSPORT.md) supplies transport for
a wide diagonal C3 part coupled to a bounded quantum core, and localizes
commuting lists with bounded relation count. Thus the former t-at-least-two
next step is historical. The commutation-cover theorem above now closes
the case where every bad B axis commutes with A, including arbitrarily
many bad zero-character axes. Genuinely mixed components, containing bad
zero-character axes and middle axes that interact with both A and C,
remain outside this certificate. The general
native three-band transport, its error radius outside C3, and the general
returned six-layer theorem remain open.

All source axes and overlaps refer to the specified physical factorization.
A different resynthesis can change them and change which partial sources
remain C3. No test here is an optimal-depth criterion or a unitary invariant.

Every equality is an operator identity on the full work space, with arbitrary
dirty inputs and external references. Native scalar phases are retained;
projective reasoning is used only for Clifford membership and binary axes.
The six-layer corollary starts from the original full-return comparison,
including all logical columns, literal target phases, leakage and actual
adjoints. Fixed coordinates and quotient observables are analytical and
introduce no physical gate or helper. All original T gates, elementary
Cliffords, controls, inverse calls and simultaneous workspace remain charged.

The [focused checks](../../tests/test_native_band_source_structure.py) cover
native commuting and noncommuting two-band sources, zero-overlap cancellation,
signed left/right axis transport, half-angle rigidity, the overlap action and
relation budget, strict source-certification fixtures, and the fully charged
dependent-axis Toffoli counterexample. The arbitrary-width results follow
from the proofs above, not from finite enumeration.
The [commutation-cover checks](../../tests/test_native_commutation_cover.py)
use a five-wire, overlap-four fixture with four composite bad zero-character
axes, a good B axis outside $`U_A`$, and a certified A axis to strip.
They retain all columns, literal phases and the actual inverse, and check
C3 membership before and after stripping. This bounded axis fixture makes
no synthesis-cost or optimal-depth claim. A separate three-wire fixture
checks good-middle elimination when the intermediate $`E_{B_0}E_A`$
is outside C3, retaining the literal two-sided factorization.
