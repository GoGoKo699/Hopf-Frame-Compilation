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
six-layer degree-78/26 theorem.

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
any of the three certificates above therefore verifies the source hypothesis
of the [native transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query).
Under the original full operator error

```math
\eta\lt\sqrt2/8,
```

the Boolean output has degree at most **78**, or **26** when the outer
target shear is visible. This is a new sufficient test for the existing
source class, with its existing error radius; it is not an unrestricted
six-layer degree theorem.

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
next step is historical. Bad zero-character interactions remain the
structural source gap after the nonzero-class closure. The general
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
