# Arbitrarily many isolated native paired cores

[Two-overlap source theorem](TWO_OVERLAP_SOURCE_STRUCTURE.md) ·
[One-overlap source theorem](ONE_OVERLAP_SOURCE_STRUCTURE.md) ·
[Native band structure](NATIVE_BAND_SOURCE_STRUCTURE.md) ·
[Signed Clifford transport](NATIVE_CLIFFORD_TRANSPORT.md)

An exact C3 source with three specified native T bands is diagonal C3 in
fixed Clifford coordinates when its noncommuting middle/outer axes split
into isolated paired cores, with no bound on their number. The native phase
formula controls signs across arbitrarily many controls. Commuting
control-X images then give a common Pauli-preserving Lagrangian. This is
an arbitrary-overlap source certificate, not a classification of all
three-band sources or an inference from complete workspace return. Full
C3 membership also excludes a single isolated wider activation orbit.

## 1. Exact native theorem

**Isolated-core theorem.** Let $`F=E_CE_BE_A`$ be C3 on the full Hilbert space. Each band has an
independent isotropic list of actual signed native T/T-dagger axes. Suppose
the complete B/C interaction consists of k isolated paired cores, and write
$`U_A`$ for the binary span of A:

```math
B_i=\{b_i,b_i+r_i\},\qquad
C_i=\{c_i,c_i+r_i\},\qquad
\langle b_i,c_i\rangle=1,\qquad 1\le i\le k.
```

Assume the $`r_i`$ are independent and belong to $`U_A`$. The spaces
$`W_i=\mathrm{span}\{r_i,b_i,c_i\}`$ are mutually symplectically orthogonal,
with $`r_i`$ radical. All remaining B/C axes commute with every core and
with one another. They form a commuting diagonal-C3 phase source G.

Then F is diagonal C3 in fixed Clifford coordinates, for arbitrary k.
Additional relations between A and good B axes are permitted. There is no
bound on total overlap, dirty width or nondepth gate counts in this theorem.

The conclusion is the source form
$`F=C_{\mathrm{out}}D C_{\mathrm{in}}`$ with fixed end Cliffords and D a
weighted cubic phase modulo eight, as defined in the
[native source theorem](NATIVE_BAND_SOURCE_STRUCTURE.md#1-common-native-coordinates-and-exact-hypotheses).
The end Cliffords need not be inverses.

The assertion concerns the specified native factorization. All signs,
orientations, literal scalar phases, actual inverses and physical resource
charges remain in the source. Coordinates and core extraction below are
analytical and supply no free gates or prepared workspace.

## 2. The rank obstruction localizes the A relations

Put $`R=\mathrm{span}\{r_1,\ldots,r_k\}`$. Write each control uniquely
in A coordinates. For an A axis participating in those relations, choose
a Pauli p isolating it, and let

```math
v=(\langle p,r_1\rangle,\ldots,\langle p,r_k\rangle).
```

Choose symplectic control coordinates so p is precisely the product of the
control X axes selected by v, and flips those entries.
Target representatives can be adjusted by R to commute with the chosen
control X axes. Such adjustments give affine control-dependent signs to
the target axes, but do not change the binary activation of a paired native
band: its target transvection is active in exactly one value of its own
control. Every outer paired-core branch is therefore an ordinary target
Clifford with its exact scalar retained.

The isolating A phase supplies one Pauli quarter-turn in the ratio between
control sectors u and $`u+v`$. The G ratio is a diagonal spectator Clifford D
with $`S_D^2=I`$. The full Clifford image $`FpF^\dagger`$ has residual
blocks with a common binary symplectic action and inverse blocks in paired
sectors. The
[single half-angle lemma](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension),
exactly as in the
[one-core proof](ONE_OVERLAP_SOURCE_STRUCTURE.md#43-full-c3-membership-removes-every-spectator-component),
yields

```math
S=M_vT_hS_D,\qquad S^2=I.
```

Here $`M_v`$ has an order-three binary block on every target with $`v_i=1`$
and identity on the other targets. The scalar half-angle replaces $`T_h`$
by identity. Let $`L=M_vS_D`$. Then $`L^2=M_v^2`$ and
$`\mathrm{rank}(S-L)\le1`$, so

```math
2\,\mathrm{wt}(v)
=\mathrm{rank}(M_v^2-I)
\le\mathrm{rank}((S-L)S)+\mathrm{rank}(L(S-L))\le2.
```

Thus every nonzero v is a unit vector. The A-coordinate supports of all
control relations are pairwise disjoint.

For $`v=e_i`$, apply the
[one-core symplectic block lemma](ONE_OVERLAP_SOURCE_STRUCTURE.md#43-full-c3-membership-removes-every-spectator-component)
with all inactive targets included among spectators. It forces h to have a nonzero
component only on target i. G acts only on original spectators and cannot
erase a Pauli factor on an inactive target or spectator. Thus each A axis
in the $`r_i`$ relation belongs to $`W_i+R`$.

If $`V_i`$ is their span, disjoint A supports give
$`V_i\cap R=\mathrm{span}\{r_i\}`$. Quotient by the other controls:
independence is preserved in a space with one symplectic plane and the
one-dimensional radical $`r_i`$. Its maximal isotropic dimension is two.
The relation therefore has at most two axes. A single axis would equal
$`r_i`$ and have zero target component, excluded by the block lemma. Each
relation is exactly a native A pair

```math
A_i=\{a_i,a_i+r_i\},\qquad a_i\in W_i+R.
```

This generalizes the
[two-core localization](TWO_OVERLAP_SOURCE_STRUCTURE.md#33-each-relation-is-one-native-a-pair)
without changing its full-space or native-independence hypotheses.

## 3. Extracting the core preserves the remaining A band

Choose a control $`X_i`$ isolating one member of $`A_i`$ among the entire
A list. Their pairings with the controls are the identity matrix. Adding
linear combinations of control Z axes makes these X axes mutually commuting
without changing any A-axis commutation. Complete the symplectic coordinates
and adjust target representatives by R.

Every remaining A axis commutes with all chosen control X and Z axes. The
selected A pairs act only on their own target and own control; the B/C pairs
may retain affine signs involving arbitrary controls. Let K contain the
bad C and B pairs and the selected A pairs. It acts on the k controls and
k target qubits, with exact sector form

```math
K=\sum_{u\in\mathbb F_2^k}|u\rangle\langle u|\otimes K(u),
\qquad K(u)=\bigotimes_{i=1}^k K_i(u),
```

including all sector scalars in the factors. Every $`K_i(u)`$ is an ordinary
one-qubit Clifford with retained eighth-root scalar. Its binary action
depends only on $`u_i`$, and its relative action between the two values is
a nonidentity transvection: it is conjugate to $`M_iT_{h_i}`$ with $`M_i`$
order three and $`h_i\ne0`$.

Write

```math
F=U_0E_{A,\mathrm{rest}},\qquad U_0=GK.
```

The remaining A list is still an actual independent native band and commutes
with all $`X_i`$. Consequently

```math
J_i=U_0X_iU_0^\dagger=FX_iF^\dagger
```

is Clifford. These are commuting Hermitian involutions. This identity
retains all wide A coupling: no assumption that the remaining A axes commute
with every target is made.

## 4. Exact native pair phases and the target-sign module

Let

```math
N_e(P)=\frac{I+P}{2}+\omega^e\frac{I-P}{2},\qquad
\omega=e^{i\pi/4},\qquad e\in\mathbb Z,\qquad S_P=N_2(P).
```

In fixed control sector u, a paired factor on target i can be written

```math
C(u)=N_e((-1)^{\ell(u)}P)
     N_f((-1)^{\ell(u)+u_i}P),
```

where $`e,f\in\{1,-1\}`$ and $`\ell`$ is an affine Boolean parity. A fixed discrepancy between
the two actual signed axes is absorbed using
$`N_f(-Q)=\omega^fN_{-f}(Q)`$, retaining that fixed scalar. Define the
Boolean activation

```math
a(u_i)=
\begin{cases}
1+u_i,&e=f,\\
u_i,&e=-f.
\end{cases}
```

The exact identity $`S_P^{-e}=S_P^eP`$ gives

```math
C(u)=\gamma(u)S_P^{e a(u_i)}P^{a(u_i)\ell(u)},
```

with literal scalar

```math
\gamma(u)=\omega^{e\ell(u)+f(\ell(u)\oplus u_i)}
=\omega^{f u_i+(e+f-2f u_i)\ell(u)}.
```

The second exponent is an integer expression modulo eight with the Boolean
values represented by 0 and 1. It is not a silent projective equality.
All these scalars remain in K and its control-X images; they cancel only
in the particular target-Pauli adjoint being analyzed.

Suppose the incoming target Pauli label depends only on $`u_i`$. Its
commutation bit with P is a function $`c(u_i)`$. The extra Pauli factor
contributes sign $`a(u_i)c(u_i)\ell(u)`$. The baseline Clifford
$`S_P^{e a(u_i)}`$ contributes a sign depending only on $`u_i`$ and leaves
the outgoing binary label dependent only on $`u_i`$. Since a Boolean
function of one bit is affine, the sign belongs to the module

```math
\mathcal M_i=\mathrm{Aff}(u)+u_i\mathrm{Aff}(u).
```

Composing the three paired native factors preserves this module. Thus every
target-Pauli image through K has a binary label depending only on its own
control and a sign of degree at most two. More precisely, every quadratic
monomial in that sign contains $`u_i`$; a product $`u_j u_l`$ with both
indices different from i never arises.

## 5. The extracted source U0 is C3 for arbitrary core count

Control Z axes are fixed. Their X images are the Cliffords $`J_i`$ supplied
by the full source. A target Pauli image has two Hermitian Pauli branches
selected by its own control, times the quadratic sign from section 4.
The two-branch Pauli operator is Clifford: factor out one branch and use a
controlled Pauli with its allowed fourth-root scalar. The quadratic sign
is a diagonal Clifford on all controls. G commutes with both factors.
Therefore every target Pauli has a Clifford image through $`U_0`$.

K acts trivially on spectators, so their $`U_0`$ images are their Clifford
images through the commuting diagonal-C3 source G. All Pauli generators
now have Clifford images, proving $`U_0\in C_3`$ with no bound on k.

This replaces the previous two-control observation that every Boolean sign
on two bits is quadratic. The native phase identity supplies the required
quadratic bound at arbitrary control count.

## 6. All control-X commutators fix the target signs

For each target i choose a nonzero input Pauli $`p_i`$ fixed by its relative
branch transvection. Its two branch images have the same output target
label $`q_i`$. Retaining all signs gives

```math
Q_i=U_0\tau(p_i)U_0^\dagger
=(-1)^{f_i(u)}\tau(q_i).
```

Let $`S_j`$ be the symplectic action of $`J_j`$. It fixes the control-Z
space R pointwise. Its induced action on target i is $`T_{q_i}`$ when
$`j=i`$, and identity when $`j\ne i`$. Target labels acquire only
control-Z tails, not spectator components: the native block factorization
and G's commutations give this directly, and $`J_j`$ is Clifford.

Choose a target-i label w anticommuting with $`q_i`$. Then

```math
S_iw=w+q_i+r_{iw},\qquad
S_jw=w+r_{jw}\quad(j\ne i),\qquad
r_{iw},r_{jw}\in R.
```

The involution $`S_i^2=I`$, with R fixed, gives $`S_iq_i=q_i`$. For
$`j\ne i`$, the exact commuting operators $`J_i,J_j`$ have commuting
symplectic maps. Evaluating $`S_iS_jw=S_jS_iw`$ and cancelling the fixed
R tails gives

```math
S_jq_i=q_i\qquad\text{for every }j.
```

Hence each $`J_j`$ sends $`\tau(q_i)`$ to itself up to a constant sign,
with no control-Z factor. The operators $`J_j,Q_i`$ commute exactly because
their input Paulis do. Since $`J_j`$ flips only control j, it follows that

```math
f_i(u+e_j)+f_i(u)=\text{constant}\qquad\text{for all }j.
```

Every Boolean derivative is constant, so $`f_i`$ is affine. Therefore all
$`Q_i`$ are signed Paulis. This conclusion actually uses no degree bound
on $`f_i`$: the cross-control commutators eliminate even a hypothetical
quadratic on two other controls. As a separate check, section 4's stronger
native sign module already lets the own-control derivative eliminate every
possible quadratic monomial.

The controls, all input target axes $`p_i`$, and a full spectator Z basis
form a Pauli-preserving Lagrangian for $`U_0`$. A fixed Clifford matches
their literal signed images. The resulting source is diagonal, and its
proved C3 membership gives the weighted cubic phase form with retained
scalar. Finally
[native diagonal right extension](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
applies to
$`F=U_0E_{A,\mathrm{rest}}`$. Thus F is diagonal C3.

## 7. A separate controlled-Clifford certificate

The commutator argument gives a useful exact statement independent of native
sign synthesis. Suppose $`U=GK`$ is C3, where G is diagonal on controls and
spectators and identity on targets, and K has one target qubit per control
with ordinary Clifford fibers $`K(u)=\bigotimes_iK_i(u)`$. Include all
sector scalars. Suppose the binary action of $`K_i(u)`$ depends only on
$`u_i`$ and its relative two-branch action is a nonidentity transvection.
Then U is diagonal C3 in fixed Clifford coordinates.

The fixed target axes and commuting control-X images in section 6 give
the Pauli-preserving Lagrangian, while the assumed C3 property supplies
the weighted cubic diagonal form. Arbitrary scalar phases are not ignored:
the full C3 premise constrains their control-X images. In the native source
theorem, section 5 proves this C3 premise from the exact native formula and
the full source's isolating probes.

The nonidentity relative action is essential. A constant binary target
action does not supply the transvection used to remove its control-Z tails.
The [known Gottesman–Mochon source](CONTROLLED_PAULI_TRANSPORT.md#8-known-seven-qubit-example-with-literal-source-accounting)
has controlled Pauli fibers, hence constant binary target action, and is
not semi-Clifford. The certificate therefore does not extend to arbitrary
controlled Clifford fibers with identity relative action.

## 8. An activation certificate at arbitrary overlap

The isolated-core hypotheses follow from a directly checkable condition.
For clarity, the activation code and coordinate characters are

```math
\mathcal U=\{(\langle b_j,p\rangle)_j:p\in U_A^\perp\},\qquad
\lambda_j:\mathcal U\longrightarrow\mathbb F_2,\quad \lambda_j(u)=u_j.
```

A zero character means $`b_j\in U_A`$. Equal nonzero characters form a
repeated class. Call a B axis good when it commutes with every C axis.
Assume the exact native source is C3, every zero-activation B axis is good
(commutes with every C axis), and every C axis anticommutes with at most two
B axes. Then the source is diagonal C3, with no restriction on total overlap.

Indeed a bad C neighborhood contains no zero-character axis. The
[native activation-orbit proof](NATIVE_BAND_SOURCE_STRUCTURE.md#5-three-band-overlap-restricts-noncommuting-middle-axes)
puts that neighborhood in one nonzero coordinate-character class. It is
nonempty and even by isotropy of the orbit pair, hence has
exactly two B axes. Their sum is the C-pair difference and belongs to
$`U_A`$ because the coordinate characters agree.

Two different two-element neighborhoods cannot overlap once: the difference
of the first C pair would anticommute with a C axis in the second orbit,
contradicting C isotropy. Equal neighborhoods cannot support two C orbits,
because their four distinct C axes would sum to zero, contradicting native
independence. Thus all neighborhoods are disjoint paired cores. Disjoint
B supports make their control differences independent; native B/C
commutation makes their core spaces orthogonal, and every leftover phase
centralizes the cores and other leftover phases. Section 1 applies.

In particular, if there are no zero coordinate characters and every nonzero
character class has size at most three, the certificate holds. Every bad
C neighborhood is a nonempty even subset of one such class, so has size
two. This gives an arbitrary-rank activation-family consequence, not just
an isolated finite-rank example.

## 9. Orbit incidence without an overlap bound

Within one nonzero coordinate-character class, every nontrivial C orbit
has two members $`\{c_j,c_j+r_j\}`$. If $`J_j`$ is its B neighborhood,
then

```math
r_j=\sum_{i\in J_j}b_i,\qquad |J_j|=0\pmod2.
```

The neighborhood is nonempty. For distinct orbits, C isotropy gives
$`|J_j\cap J_l|=0\pmod2`$. The differences $`r_j`$ are independent:
a relation among them would be a relation among disjoint pairs of the
independent C list. B independence then makes the binary incidence rows
independent. Thus the orbit-incidence matrix H of each repeated class obeys

```math
\mathrm{rank}(H)=\text{number of rows},\qquad HH^T=0.
```

These necessary conditions do not force weight-two rows. The row 1111
is permitted; so are the two independent rows 111100 and 001111, whose
supports overlap in two places. Such data have independent isotropic
realizations: take $`b_i=Z_i`$ and the C pairs
$`\{X^{J_j},X^{J_j}Z^{J_j}\}`$, retaining their Hermitian signs.
The X and Z components separately prove independence of the C list.
An independent basis of the even-parity Z subspace in A gives all B
coordinates the same nonzero activation character. These are axis data,
not certificates that the resulting full source is C3.

## 10. A single wider orbit is excluded by full C3 membership

**Single-orbit theorem.** Suppose all actual B/C noncommutation is one
C orbit $`\{c,c+r\}`$, whose neighborhood J lies in one nonzero
coordinate-character class. Every other B/C phase commutes with that core
and with every other good phase. For independent isotropic native A, B, C
lists, exact full-space C3 membership forces $`|J|=2`$. Consequently F is
diagonal C3, with no bound on its total A/B overlap.

In particular there are no bad zero-character B axes in this statement.
The proof excludes one isolated wider orbit; it does not separate several
wider orbits with overlapping neighborhoods.

### 10.1 Exact control-sector ratio and the half-angle axis

Write $`k=|J|`$, a positive even integer, and
$`r=\sum_{j\in J}b_j\in U_A`$. Choose an A axis $`a_i`$ in r's unique
A-coordinate expansion, and a Pauli p isolating that axis in the full A
list. Then $`\langle p,r\rangle=1`$. Use r,p as a control Z,X pair,
and choose the representative c of its C orbit to commute with p. Write

```math
b_j=\bar b_j+u_jr,\qquad
u_j=\langle p,b_j\rangle,\qquad
a_i=r+a_*,\qquad
\sum_{j\in J}u_j=1.
```

All barred labels commute with both control axes. They have exactly one
relation, $`\sum_{j\in J}\bar b_j=0`$. For $`k\ge4`$ they are nonzero
and pairwise distinct, and every one anticommutes with c.

Let $`B_+,B_-`$ be the two restrictions of the actual native J-band.
They need not be Clifford. Their literal ratio
$`Q=B_+B_-^\dagger`$ is Clifford: each sign-changing native phase
contributes a quarter-turn about its barred axis. Its binary action is

```math
S_Q=\prod_{j:u_j=1}T_{\bar b_j}.
```

Let G be the product of all good B/C phases. It centralizes every barred
B axis, c and r. Its sector ratio $`D=G_+G_-^\dagger`$ is a diagonal
Clifford with $`S_D^2=I`$, and it fixes every barred B axis and c. The
native C pair restricts to ordinary Cliffords $`C_+,C_-`$, with relative
binary action $`T_c`$. The isolating p leaves just one A phase in the
sector ratio, so the residual block of $`FpF^\dagger`$ is literally

```math
V=C_+Q\,
 \bigl[B_-G_+e^{\pm i\pi\tau(a_*)/4}
       G_+^\dagger B_-^\dagger\bigr]
 C_-^\dagger D.
```

All scalar phases of the native restrictions remain in this equality.
The full source is C3, so V is Clifford. Its paired blocks have a common
binary action and are inverses, giving $`S_V^2=I`$ on the entire remaining
space.

Every factor outside the brackets is Clifford. For $`a_*\ne0`$ define

```math
H=B_-G_+\tau(a_*)G_+^\dagger B_-^\dagger.
```

This H is a traceless Hermitian involution even though $`B_-`$ need not
be Clifford. The [half-angle lemma](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
therefore makes H a signed Pauli: the projective Clifford
$`e^{\pm i\pi H/4}=(I\pm iH)/\sqrt2`$ has identity Pauli coefficient
$`1/\sqrt2`$, and flat Clifford Pauli coefficients force exactly two
terms. No prior Clifford property of H is assumed. If $`a_*=0`$, the
bracket is scalar; the order-three obstruction below excludes this case.

### 10.2 Reduction to one target plane keeps the spectator action

Put $`v=\sum_{j:u_j=1}\bar b_j`$. The selected cardinality is odd, so
$`\langle v,c\rangle=1`$. On binary labels factor

```math
S_Q=T_vQ_{\mathrm{rad}},\qquad
Q_{\mathrm{rad}}=T_vS_Q.
```

Both $`S_Q`$ and $`T_v`$ fix v and send c to $`c+v`$. Because all
barred B labels commute, $`Q_{\mathrm{rad}}`$ is involutory and fixes
v,c pointwise. It commutes with $`S_D`$, since D fixes every barred B
axis. The fixed target plane is therefore symplectically orthogonal to
the remaining action of these two maps.

To make the conjugation direction explicit, D commutes with $`C_-`$ and

```math
C_-^\dagger V C_-
 =(C_-^\dagger C_+)Qe^{\pm i\pi H/4}D.
```

If h is the signed Pauli H's binary label, move
$`Q_{\mathrm{rad}}`$ through its transvection. The resulting involutory
binary map is

```math
S=(T_cT_v)T_{h'}D',\qquad
h'=Q_{\mathrm{rad}}h,\qquad
D'=Q_{\mathrm{rad}}S_D,\qquad
(D')^2=I.
```

The map $`D'`$ is identity on $`\mathrm{span}\{v,c\}`$ and acts only
on its symplectic complement. The target map $`T_cT_v`$ has order three.
The [one-core involution calculation](ONE_OVERLAP_SOURCE_STRUCTURE.md#43-full-c3-membership-removes-every-spectator-component)
forces $`h'`$ to be nonzero and entirely in that target plane. A scalar
half-angle cannot remove the order-three block. Since
$`Q_{\mathrm{rad}}`$ fixes the plane pointwise, $`h=h'`$ lies there too.

### 10.3 The native pullback excludes width at least four

Every sector of G fixes the target plane. Since G also commutes with
$`B_-`$, pulling back H shows that $`\tau(a_*)=B_-^\dagger H B_-`$
is a Pauli. If $`h=c`$ or $`h=c+v`$, H anticommutes with all k barred
B axes. This pullback is H times a product of quarter-turns about all
those axes, with the literal native signs retained. Their binary action is

```math
I+N,\qquad
N=\sum_{j\in J}\bar b_j\ell_{\bar b_j},\qquad
\ell_b(x)=\langle b,x\rangle.
```

For $`k\ge4`$, $`\mathrm{rank}(N)=k-2\gt0`$. To see this, let M have
the k barred labels as columns. Its kernel is the span of the all-ones
vector, and the image of $`M^T\Omega`$ is the even-parity subspace of
$`\mathbb F_2^k`$, where $`\Omega`$ is the invertible ambient symplectic
form. That image has dimension $`k-1`$ and intersects $`\ker M`$ in
one dimension, because k is even. Hence
$`\mathrm{rank}(MM^T\Omega)=k-2`$. This is exactly the asserted rank.

The quarter-turn product is consequently not Pauli, so multiplying H by
it cannot give a Pauli. The only remaining target label is $`h=v`$.
Both $`B_-`$ and $`G_+`$ fix that Pauli, giving $`a_*=v`$ at the label
level and thus

```math
a_i=r+v=\sum_{j:u_j=1}b_j.
```

This is an odd subset of one nonzero activation-character class. Its
functional on $`\mathcal U`$ is that nonzero character, so its sum cannot
belong to $`U_A`$. This contradicts $`a_i\in A`$ and excludes every
$`k\ge4`$. Only $`k=2`$ remains, where section 1 applies directly.

## 11. Scope and six-layer consequence

This establishes an arbitrary-rank source certificate for isolated paired
cores, not an arbitrary-rank classification of all activation patterns.
Multiple wider C orbits can have overlapping even neighborhoods. Although
their differences are independent conserved Paulis, one physical B phase
may then belong to several orbit cores. The isolating-probe weighted Gram
matrix has entries $`\sum_i u_iH_{ji}H_{li}`$; $`HH^T=0`$ does not make
its off-diagonal entries vanish. The proof does not classify these coupled
incidence patterns or arbitrary bad zero-activation interactions.
No inverse-C3 assumption or approximate source-membership statement enters.

Under the [actual six-layer query contract](READ_ONLY_SIX_LAYER_REDUCTION.md),
an exactly C3 specified prefix with this structure satisfies the retained
[diagonal-source transport hypothesis](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query).
At full operator error $`\eta\lt\sqrt2/8`$, the necessary
Boolean-degree bound is 78, or 26 for visible outer shear. Literal target
phase, arbitrary dirty inputs and external references, leakage, complete
return and the actual inverse remain part of that full operator contract.
All original native gates, all six T slots, elementary Cliffords and dirty
width remain charged; the theorem changes no construction's resource ledger.

## 12. Bounded checks

The [focused tests](../../tests/test_isolated_native_core_structure.py)
compare 32 native eigenvalue/orientation cases using exact
eighth-root arithmetic, and 6,144 four-control matrix blocks covering all
32 affine parities, both orientations of each native phase and all three
target Pauli axes. Every target label depends only on its own control, and
every quadratic adjoint-sign monomial contains that control. These checks
exercise signs involving more than two controls; the arbitrary-count
conclusion rests on the proof above.

The tests also check the single-orbit derivative rank at widths
2, 4, 6, 8 and 10, and eight width-four native candidates that satisfy
the activation geometry but fail full C3 membership.

An exact three-core fixture uses six arbitrary wires, zero helpers and
overlap three. Literal CZ words attach the parity of both other controls
to the first target. The tests compare every one of its 64 computational
columns, including scalar phase, with a fair one-band implementation,
and compare the actual reverse-and-adjoint word with the literal inverse.
The retained elementary ledgers are:

| Implementation | T count | Elementary Clifford count | T depth |
| --- | ---: | ---: | ---: |
| Original source | 18 | 42 | 3 |
| Actual source inverse | 18 | 42 | 3 |
| Fair baseline | 6 | 54 | 1 |
| Fair baseline inverse | 6 | 60 | 1 |

These are bounded source fixtures. They assert neither optimal depth nor
an asymptotic resource advantage or a new returned-query construction.
