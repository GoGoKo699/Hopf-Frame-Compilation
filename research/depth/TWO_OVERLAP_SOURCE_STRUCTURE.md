# Three native bands with two-dimensional overlap

[One-overlap source theorem](ONE_OVERLAP_SOURCE_STRUCTURE.md) ·
[Native band structure](NATIVE_BAND_SOURCE_STRUCTURE.md) ·
[Signed Clifford transport](NATIVE_CLIFFORD_TRANSPORT.md) ·
[Six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md)

An exact C3 source with three specified native T bands is diagonal C3 in
fixed Clifford coordinates when the first two axis spans overlap in at
most two dimensions. Two separate controlled cores do not escape that
conclusion: full C3 membership restricts their A-axis coupling and fixes
the signs needed for a common Pauli-preserving Lagrangian. The zero-axis
and mixed activation patterns reduce to absorbable duplicates; their
remaining nonduplicate cores contradict full-space Clifford involution.

## 1. Exact statement and scope

In common fixed Clifford coordinates write the specified source as

```math
F=E_CE_BE_A,\qquad t=\dim(U_A\cap U_B)=2.
```

Each factor is a Clifford conjugate of one actual physical T/T-dagger band.
Its active binary Pauli axes form an independent isotropic list, respectively
A, B and C. Axis signs, orientations and native scalar phases remain in the
actual operators. All C3 hypotheses concern the full Hilbert space: every
Hermitian Pauli has a Clifford image. No inverse-C3 or intermediate-C3
premise is made.

For the activation code of B coordinates on Paulis centralizing A, write

```math
\mathcal U=\{(\langle b_i,p\rangle)_i:p\in U_A^\perp\},\qquad
\mathcal R=\mathcal U^\perp,
\qquad \dim\mathcal R=2.
```

Let $`\lambda_i`$ be coordinate i restricted to $`\mathcal U`$, and let
$`z=|\{i:\lambda_i=0\}|=|B\cap U_A|`$. A B axis is **bad** when it
anticommutes with at least one C axis; otherwise it is good.

**Two-overlap theorem.** Suppose F is C3. Then F is diagonal C3 in fixed
Clifford coordinates. All three possibilities $`z=0,1,2`$ are covered.

The conclusion means $`F=C_{\mathrm{out}}D C_{\mathrm{in}}`$, with fixed
Cliffords and D a weighted cubic phase modulo eight, including its retained
scalar. The end Cliffords need not be inverses. The
[one-overlap theorem](ONE_OVERLAP_SOURCE_STRUCTURE.md#1-exact-statement-and-native-scope)
already gives this conclusion for $`t\le1`$.

Together these results cover $`t\le2`$. This note does not classify
$`t\ge3`$ or sources known only through complete workspace return. The
overlap concerns the specified native factorization and can change under
resynthesis. This is an exact source certificate, not a robust
source-membership test.

## 2. Activation classes and the local one-core certificates

### 2.1 Complete exceptional-class list

The [native overlap proof](NATIVE_BAND_SOURCE_STRUCTURE.md#5-three-band-overlap-restricts-noncommuting-middle-axes)
shows that each C-axis orbit has one or two elements. A two-element orbit
uses one nonzero coordinate-character class, together with any
zero-character coordinates; the latter make no contribution to its
displacement. A singleton C axis can anticommute only with zero-character
B axes.

A zero character gives $`e_i\in\mathcal R`$, and equal characters give
$`e_i+e_j\in\mathcal R`$. If the repeated nonzero classes have sizes
$`m_j`$, their independent coordinate relations give

```math
z+\sum_j(m_j-1)\le2.
```

The possible exceptional classes are therefore:

| Zero-character axes | Repeated nonzero classes |
| --- | --- |
| None | None, one pair, one triple, or two disjoint pairs |
| One | None or one pair |
| Two | None |

The table permits unused relation budget: a further relation of larger
support can create no zero or repeated coordinate character. Singleton
nonzero classes contain only good B axes.

When $`z=0`$, a nontrivial C orbit uses an even nonempty subset of its
repeated class. Its displacement is the sum of those B axes, and isotropy
of the two C members makes that subset even. A repeated pair therefore
gives an orbit $`\{c,c+b_i+b_j\}`$. C-axis independence permits at most
one such orbit for that pair.

In a repeated triple, the only nonempty even subsets are its three pairs.
Two different pairs intersect once. The displacement of one resulting C
orbit would then anticommute with a member of the other, violating C
isotropy. Thus a triple supports only one actual bad B pair; its third B
axis is good.

For two repeated pairs the nonzero characters are distinct. Each C orbit
uses one pair alone, and there is at most one orbit per pair. When both
occur, their core spaces are symplectically orthogonal. All remaining B/C
axes commute with the cores and with one another. Section 3 closes this
case.

### 2.2 Local hypotheses suffice for the earlier proofs

The one-overlap proofs use the overlap dimension only to establish their
local interaction patterns. Their subsequent arguments remain valid with
additional relations involving good axes.

If the only bad B axis is b and $`b\in U_A`$, set
$`U=E_CE_{B,\mathrm{good}}`$. This is diagonal C3, and the source is
$`UE_AE_b`$. The [weight-one proof](ONE_OVERLAP_SOURCE_STRUCTURE.md#32-closing-the-source-case)
uses only the unique A-coordinate expansion of b: a duplicate axis merges
to a retained Clifford; a composite expansion contradicts native
support-coset rigidity. Thus this local pattern is diagonal C3 at any
overlap.

Likewise suppose the full bad B/C interaction consists of

```math
\{b,b+r\},\qquad\{c,c+r\},\qquad
r\in U_A,\qquad\langle b,c\rangle=1,
```

with all remaining B/C axes commuting with this core and each other. The
[weight-two proof](ONE_OVERLAP_SOURCE_STRUCTURE.md#4-a-weight-two-relation)
then applies to the unique A expansion of r. It localizes that expansion
to the core, forces exactly one A pair, extracts a phase-sensitive
diagonal-C3 source and applies native right extension. No premise on the
dimension of other A/B relations is used in those steps.

These local certificates close the table's no-zero cases except two actual
bad pairs. They also close the one-zero case whenever either its zero axis
or its repeated pair is good. A genuinely mixed interaction is not two
orthogonal controlled cores: its paired C orbit can anticommute with the
zero B axis, and singleton C axes can also anticommute with that axis.
Section 5 treats this remaining pattern separately.

## 3. Two orthogonal controlled cores

### 3.1 Native core geometry

Suppose both repeated pairs interact. Label their B and C axes by

```math
B_i=\{b_i,b_i+r_i\},\qquad
C_i=\{c_i,c_i+r_i\},\qquad
\langle b_i,c_i\rangle=1,\qquad i=1,2.
```

The controls $`r_1,r_2`$ are independent and belong to $`U_A`$. The core
spaces $`W_i=\mathrm{span}\{r_i,b_i,c_i\}`$ are symplectically orthogonal,
with $`r_i`$ radical in $`W_i`$. Write
$`R=\mathrm{span}\{r_1,r_2\}`$. The remaining B/C phases form a commuting
diagonal-C3 source G centralizing both cores.

The [exact paired-phase formula](ONE_OVERLAP_SOURCE_STRUCTURE.md#42-isolate-one-a-coordinate-of-the-relation)
will be used throughout. Each paired band restricts in its control sectors
to ordinary Clifford quarter-turns or scalars, retaining powers of
$`\omega=e^{i\pi/4}`$. Its binary transvection occurs in exactly one
value of its own control. Common axis-sign changes alter the quarter-turn
orientation and scalar, but not this binary activation.

### 3.2 One A coordinate cannot flip both controls

Write both controls uniquely in the independent A list. For an A axis
$`a_j`$ participating in either relation, choose a Pauli p isolating it:
$`\langle p,a_k\rangle=\delta_{jk}`$. Its control-flip vector is

```math
v=(\langle p,r_1\rangle,\langle p,r_2\rangle)\in\mathbb F_2^2.
```

Choose control X coordinates so p is their product specified by v.
Representatives of the target axes can be adjusted by R to commute with
those control X axes. This preserves both target planes and introduces
only control-dependent signs in the native rotations. In a fixed control
sector $`u=(u_1,u_2)`$, the bad outer bands give a tensor product of target
Cliffords $`R_u`$. The binary map

```math
M_v=S_{R_{u+v}}^{-1}S_{R_u}
```

is order three on target i when $`v_i=1`$, and identity there when
$`v_i=0`$. Each active block is either
$`T_{b_i}T_{c_i}`$ or $`T_{c_i}T_{b_i}`$.

All A axes commute with R. Since p isolates $`a_j`$, the sector ratio
$`E_{A,u}E_{A,u+v}^{\dagger}`$ is one Pauli quarter-turn, possibly scalar.
The ratio of the G sectors is a diagonal spectator Clifford D, with
$`S_D^2=I`$. Full C3 membership makes $`FpF^\dagger`$ Clifford. Its
residual blocks $`V_u=F_uF_{u+v}^{\dagger}`$ are therefore Clifford with
a common binary action; the block at $`u+v`$ is the actual inverse, so
$`S_{V_u}^2=I`$.

Cancel the outer target Cliffords and D. The
[single half-angle lemma](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
makes the G-conjugated quarter-turn axis a Pauli. After conjugating by a
target Clifford, the necessary condition is

```math
S=M_vT_hS_D,\qquad S^2=I.
```

The scalar case uses the identity in place of $`T_h`$. Put
$`L=M_vS_D`$. The maps on targets and spectators commute, so
$`L^2=M_v^2`$, while $`\mathrm{rank}(S-L)\le1`$. Hence

```math
\mathrm{rank}(S^2-L^2)
\le\mathrm{rank}((S-L)S)+\mathrm{rank}(L(S-L))\le2.
```

For $`v=(1,1)`$, both target blocks have order three and
$`\mathrm{rank}(M_v^2-I)=4`$, a contradiction. The A-coordinate supports
expressing $`r_1`$ and $`r_2`$ are therefore disjoint.

### 3.3 Each relation is one native A pair

For $`v=e_i`$, apply the
[one-core symplectic block argument](ONE_OVERLAP_SOURCE_STRUCTURE.md#43-full-c3-membership-removes-every-spectator-component),
including the inactive target among the spectators. Its spectator map is
identity on that target and $`S_D`$ on the original spectators, so remains
involutory. The argument forces h to be nonzero on target i and zero on
every spectator and inactive target. G cannot remove a nonidentity factor
on either. Every selected A axis thus lies in

```math
W_i+R=\mathrm{span}\{r_1,r_2,b_i,c_i\}.
```

Let $`V_i`$ span the A axes in the $`r_i`$ relation. Disjoint relation
supports give $`V_i\cap R=\mathrm{span}\{r_i\}`$. Quotienting by the
other control preserves their independence. The quotient of $`W_i+R`$ has
one nondegenerate target plane and one-dimensional radical, so its maximal
isotropic dimension is two. The relation has at most two A axes. A single
axis would equal $`r_i`$ and have zero target component, already excluded.
Thus each relation is exactly

```math
A_i=\{a_i,a_i+r_i\},\qquad a_i\in W_i+R.
```

Cross-control sign factors are permitted here. In particular, $`a_i`$ may
contain the other control. They will be retained in the extracted core.

### 3.4 The extracted source is C3

Choose a control $`X_i`$ isolating one member of $`A_i`$ among the entire
A list. These choices are dual to the controls because the supports are
disjoint. They can be made commuting: adding $`r_1`$ to $`X_2`$ changes
their mutual pairing without changing any A-axis commutation. Adjust
target representatives by R and complete the symplectic coordinates.

Every remaining A axis commutes with both controls and both chosen X axes.
The selected A pairs act on their respective target and own control; the
outer pairs may retain cross-control signs. Let K be the product of the
bad C and B pairs and these two A pairs. It acts on the four-qubit core.
The literal source factors as

```math
F=U_0E_{A,\mathrm{rest}},\qquad U_0=GK.
```

The remaining A list is still an actual independent native band. Since it
commutes with both $`X_i`$,

```math
J_i=U_0X_iU_0^\dagger=FX_iF^\dagger
```

is Clifford. In each control sector K is a tensor product of two exact
one-qubit Cliffords, including its scalar. The binary action on target i
depends only on $`u_i`$. Its relative action between the two values is
conjugate to $`M_iT_{h_i}`$, with $`M_i`$ order three and $`h_i\ne0`$.
This is a nonidentity involution, hence a transvection. Identity would
require an order-three map to equal a transvection.

Control Z axes are fixed, and their X images are Clifford by the displayed
identity. A target Pauli image under K has a binary target label depending
only on its own control and an extra sign depending on both controls. An
operator with two Hermitian Pauli branches on one control is Clifford:
factoring out one branch leaves a controlled Pauli with a fourth-root
scalar. The remaining sign is $`(-1)^{f(u_1,u_2)}`$. Every Boolean function
of two bits is quadratic, so this is a diagonal Clifford. G commutes with
these controlled target-Pauli operators. All target Pauli images under
$`U_0`$ are therefore Clifford.

K acts trivially on spectators, so their images under $`U_0`$ are their
Clifford images through G. All Pauli generators now have Clifford images:
$`U_0`$ is C3. The use of two controls is essential at this step; arbitrary
signs on more controls need not be Clifford phases.

### 3.5 The fixed target-axis signs are affine

For target i choose a nonzero input Pauli $`p_i`$ fixed by its relative
branch transvection. Its branch images have the same target label $`q_i`$.
Keeping their literal signs, write

```math
Q_i=U_0\tau(p_i)U_0^\dagger
=(-1)^{f_i(u_1,u_2)}\tau(q_i).
```

The Hermitian Clifford involution $`J_i`$ fixes all control-Z binary labels.
Modulo additions of control-Z labels, it acts on target i as
$`T_{q_i}`$. No target image gains a spectator component: K has the
displayed block form and G commutes with controls and targets.

Put $`N=S_{J_i}-I`$. Involution gives $`N^2=0`$. For a target label w
anticommuting with $`q_i`$,

```math
Nw=q_i+r_w,\qquad r_w\in R.
```

Since N vanishes on R, another application gives $`Nq_i=0`$. Thus
$`J_i\tau(q_i)J_i^\dagger=\pm\tau(q_i)`$ with one constant sign and
no added control Z. The operators $`J_i,Q_i`$ commute exactly because
their input Paulis do. As $`J_i`$ flips only control i, this yields

```math
f_i(u+e_i)+f_i(u)=\text{constant}.
```

The coefficient of $`u_1u_2`$ vanishes. Hence $`f_i`$ is affine, and
$`Q_i`$ is a signed Pauli.

The two controls, the two input target axes $`p_i`$, and a full commuting
spectator Z basis form a Pauli-preserving Lagrangian for $`U_0`$. A fixed
Clifford matches its literal signed images. In fixed Clifford coordinates
$`U_0`$ is diagonal; its proved C3 membership gives the weighted cubic
phase form with retained scalar. Apply the
[native diagonal right-extension corollary](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
to $`F=U_0E_{A,\mathrm{rest}}`$. This proves the two-core case and finishes
all no-zero-character cases.

## 4. Two zero-character B axes

### 4.1 Membership classes and duplicate absorption

Let the zero-character axes be $`b_0,b_1`$. Every other B axis is good.
If either zero axis is also good, the local one-axis certificate applies.
Assume both are bad and set

```math
U=E_CE_{B,\mathrm{good}},\qquad
F=UE_AE_{b_0}E_{b_1}.
```

U is diagonal C3. Both $`b_j`$ lie in $`U_A`$, so their native phases
commute with the entire A band. Their Hermitian images
$`H_{b_j}=Ub_jU^\dagger=E_Cb_jE_C^\dagger`$ are non-Pauli, because
each anticommutes with a nonempty independent set of C axes.

The unique A expansions partition A into membership classes
$`I_{00},I_{10},I_{01},I_{11}`$, where

```math
b_0=\sum_{i\in I_{10}}a_i+\sum_{i\in I_{11}}a_i,
\qquad
b_1=\sum_{i\in I_{01}}a_i+\sum_{i\in I_{11}}a_i.
```

For $`i\in I_{10}`$, an isolating A probe activates exactly $`a_i,b_0`$
in the inner commuting list. Cancelling its Clifford U image gives the
[paired half-angle premise](ONE_OVERLAP_SOURCE_STRUCTURE.md#2-paired-half-angle-rigidity).
When $`a_i\ne b_0`$, the
[native support-coset consequence](ONE_OVERLAP_SOURCE_STRUCTURE.md#31-a-support-coset-consequence-for-a-native-image)
then forces exactly two C axes anticommuting with $`b_0`$ and

```math
a_i+b_0=c_0+c_0'.
```

There is at most one such $`a_i`$. The same argument applies to
$`I_{01}`$ and $`b_1`$.

If $`I_{11}`$ is empty, $`I_{10}`$ must consist of the single duplicate
axis $`b_0`$: two or more members would all be distinct from $`b_0`$ and
the displayed equality would make them identical. Likewise
$`I_{01}=\{b_1\}`$. The actual A/B rotations on each duplicate axis
combine into a scalar times a Clifford, commuting with the remaining A
band. Thus

```math
F=UE_{A,\mathrm{rest}}C_{\mathrm{dup}}.
```

Cancel the actual right Clifford, apply native diagonal right extension,
and restore the Clifford and scalar. F is diagonal C3. This branch permits
unbounded raw C-axis spread; no bound on the original C-axis count follows.

### 4.2 The only nonduplicate pattern has a three-qubit core

Suppose $`I_{11}`$ is nonempty. Independence rules out a duplicate in
$`I_{10}`$ or $`I_{01}`$, so each has at most one member. At least one
is nonempty because $`b_0\ne b_1`$. If only $`I_{10}=\{a_0\}`$ were
nonempty, then $`a_0+b_0=b_1`$. The paired half-angle lemma would make
$`H_{a_0}H_{b_0}=\pm H_{b_1}`$ Pauli, contradicting badness. The opposite
case is identical. Hence both classes contain one axis.

Put $`d=\sum_{i\in I_{11}}a_i\ne0`$. Then

```math
a_0=b_0+d,\qquad a_1=b_1+d.
```

The support-coset consequence identifies d as the sum of the two C axes
anticommuting with each $`b_j`$. C independence makes those unordered pairs
equal. Thus one common C pair $`\{c,c+d\}`$ anticommutes with both
$`b_0,b_1`$, and every other C axis commutes with both.

The labels $`b_0,b_1,d`$ are independent by their disjoint A membership
supports; c is independent of them because it anticommutes with $`b_0`$.
Their span has symplectic rank two and radical
$`\mathrm{span}\{d,b_0+b_1\}`$. All remaining B/C phases centralize
this space. This is a three-qubit core in analytical coordinates, without
any initialized-work premise.

### 4.3 Full C3 membership excludes the nonduplicate pattern

Choose $`a_i\in I_{11}`$ and a Pauli p isolating it among A. Then p
anticommutes with $`d,b_0,b_1`$ and commutes with $`a_0,a_1`$ and
$`r=b_0+b_1`$. Use d,p as control Z,X. Set $`b=a_0=b_0+d`$, which
commutes with p, and choose the member c of the C pair commuting with p.
The target axes b,c anticommute and commute with both control axes. In
these coordinates,

```math
B_{\mathrm{bad}}=\{b+d,b+d+r\},\qquad
A_{\mathrm{selected}}=\{b,b+r\},\qquad
C_{\mathrm{bad}}=\{c,c+d\}.
```

The label r lies on the remaining space and commutes with b,c. Write
$`a_i=d+a_*`$. Every A axis, hence $`a_*`$, commutes with b and r.
Only the isolated A phase survives between the two d sectors, giving the
retained scalar times $`\exp(\pm i\pi\tau(a_*)/4)`$. It commutes with
both restrictions of the bad B band.

Let $`C_+,C_-`$ be the bad C restrictions. They are ordinary Clifford
quarter-turns or scalars, with relative binary action $`T_c`$. Let
$`B_+,B_-`$ be the bad B restrictions. Reversing both B-axis signs between
sectors makes $`B_-`$ a retained scalar times $`B_+^\dagger`$. Thus

```math
Q=B_+B_-^\dagger\quad\text{is Clifford},\qquad
S_Q=T_bT_{b+r}.
```

This follows by squaring the two native phases, for either orientation;
the individual $`B_\pm`$ need not be Clifford. The remaining phase source
$`G=E_{C,\mathrm{good}}E_{B,\mathrm{good}}`$ centralizes b,c,r,d. Its
sector ratio $`D=G_+G_-^\dagger`$ is a diagonal Clifford fixing b,c,r
literally. The exact relative source block, with its scalar retained, is

```math
V=F_+F_-^\dagger
=C_+Q\,[G_+\exp(\pm i\pi\tau(a_*)/4)G_+^\dagger]C_-^\dagger D.
```

F commutes with d. Its C3 image of p has off-diagonal blocks V and
$`V^\dagger`$, so V is Clifford and $`S_V^2=I`$. Every factor outside
the brackets is Clifford. The single half-angle lemma therefore makes
$`H=G_+\tau(a_*)G_+^\dagger`$ a signed Pauli when $`a_*\ne0`$.
It commutes with b,r, and its label has the form

```math
h=\delta b+s,\qquad \delta\in\mathbb F_2,
```

where s has no target component and commutes with r. In the scalar case
take $`\delta=0,s=0`$ and replace its transvection by identity.

Conjugating $`S_V`$ by $`S_{C_-}^{-1}`$ gives

```math
S=T_cT_bT_{b+r}T_hS_D.
```

The map $`S_D`$ fixes b,c,r. The transvection $`T_h`$ fixes b,r, and
$`T_bT_{b+r}`$ fixes s. Direct evaluation gives

```math
Sb=b+c,\qquad
Sc=c+r+\delta(b+c+s),\qquad Sr=r,
```

and hence

```math
S^2b=b+r+\delta(b+c+s).
```

For $`\delta=0`$ this differs from b by nonzero r. For $`\delta=1`$,
equality would require $`s=b+c+r`$, impossible because s has no target
component. Thus S is not involutory, contradicting full C3 membership.
The nonduplicate pattern is excluded. Section 4.1 proves the complete
two-zero-axis case.

## 5. One zero axis and one repeated pair

Let e be the zero-character B axis and $`r=b+b'`$ the repeated-pair
relation. Assume both e and the pair are bad; section 2.2 covers the other
cases. The C list has exactly one nontrivial orbit $`\{c,c+r\}`$.
Every other C axis commutes with b,b', but may anticommute with e. Every
B axis outside $`\{e,b,b'\}`$ is good, and these three B axes are
independent.

Use the unique A expansions to define membership classes by

```math
e=\sum_{i\in I_{10}}a_i+\sum_{i\in I_{11}}a_i,
\qquad
r=\sum_{i\in I_{01}}a_i+\sum_{i\in I_{11}}a_i.
```

### 5.1 Localization and paired probes

For $`i\in I_{01}`$, an A-isolating Pauli p anticommutes with r and
commutes with e. Move the actual $`E_e`$ phase through the commuting A/B
factors to the far right. It cancels literally from $`FpF^\dagger`$.
The remaining comparison has the single controlled pair
$`\{b,b+r\},\{c,c+r\}`$, the actual A band and a commuting good phase
centralizing $`W=\mathrm{span}\{r,b,c\}`$. The
[one-core block argument](ONE_OVERLAP_SOURCE_STRUCTURE.md#43-full-c3-membership-removes-every-spectator-component)
uses only this isolating probe's Clifford image. It does not require C3
membership of the source after cancellation. Consequently

```math
a_i\in W\qquad(i\in I_{01}).
```

The zero-axis gate remains in the original source; only this adjoint
comparison cancelled it. In particular $`I_{10}`$ is nonempty. Otherwise
$`e+r=\sum_{i\in I_{01}}a_i\in W`$, making $`e\in W`$. But
$`W\cap U_B=\mathrm{span}\{r,b\}`$: c anticommutes with b, whereas
$`U_B`$ is isotropic. This contradicts independence of e,b,b+r.

For $`i\in I_{10}`$, an A-isolating p has B activation
$`u_e=1,u_b=u_{b'}`$. The centralizer activation code permits the common
pair activation and all good B activations freely. Adjust p by a Pauli
centralizing A to make it commute with every B axis except e.

Set $`V=E_CE_{B\setminus\{e\}}`$, so $`F=VE_AE_e`$. Then
$`VpV^\dagger=E_CpE_C^\dagger`$ is Clifford. Moreover
$`Va_iV^\dagger`$ and $`VeV^\dagger`$ are Clifford, because the inner
commuting list $`E_AE_e`$ fixes both axes and F is C3. Cancelling the
Clifford image of p gives the paired half-angle premise for these two
Hermitian Clifford involutions.

If $`a_i=e`$, it is an exact duplicate, handled in section 5.4. Otherwise
the axes are distinct, commuting and orthogonal in normalized trace. Since
$`VeV^\dagger=E_CeE_C^\dagger`$ is non-Pauli, paired rigidity forces
exactly two C axes to anticommute with e. The image of $`a_i`$ has the
same four-term support, and their product is a Pauli whose label d is the
sum of those two C axes, by the native support-coset argument.

These two C axes are either both orbit members or both singleton axes.
They cannot be one of each: e commutes with r and hence anticommutes with
both orbit members or neither. Their sum d commutes with every B axis.
For the orbit pair it is r; singleton axes commute with every B axis except
possibly e, whose two anticommutations cancel in d. It also commutes with
all C axes. Thus V fixes d literally, and pulling the product back gives

```math
a_i+e=d.
```

All nonduplicate members of $`I_{10}`$ would be identical. Since A is
independent and this class is nonempty, $`I_{10}=\{a_0\}`$. In the
nonduplicate case, $`I_{11}`$ is nonempty and

```math
d=a_0+e=\sum_{i\in I_{11}}a_i.
```

### 5.2 Two singleton C axes are impossible

If the two bad C axes were singletons, their sum d would centralize W.
On the other hand,

```math
d=r+\sum_{i\in I_{01}}a_i\in W.
```

The radical of W is the line spanned by r. Nonzero d must therefore equal
r. The two singleton axes would sum to the same label as the orbit pair,
contradicting independence of the complete C list.

The only nonduplicate possibility has e anticommuting with exactly the
two orbit members. Then

```math
d=r,\qquad a_0=e+r,\qquad I_{01}=\varnothing,
\qquad r=\sum_{i\in I_{11}}a_i.
```

Every other C axis now commutes with all three bad B axes. Together with
the good B phases they form a commuting diagonal-C3 source G centralizing
$`\mathrm{span}\{e,b,c,r\}`$.

### 5.3 The final nonduplicate core is not involutory

Choose $`a_i\in I_{11}`$ and an A-isolating p. It anticommutes with e,r
and commutes with $`a_0=e+r`$. Use r,p as control Z,X. Choose representatives
b,c in their actual pairs commuting with p and put $`a=a_0=e+r`$.
The labels a,b commute, each anticommutes with c, and all three commute
with both control axes. Their nonzero sum $`h=a+b`$ is radical in their
three-dimensional span. It is nonzero by independence of e,b,b+r.

Move the actual A factor $`E_a`$ next to $`E_e`$. This commutes actual
operators: a belongs to the isotropic B span and commutes with every A
axis. The literal source becomes

```math
F=G\,E_{\{c,c+r\}}E_{\{b,b+r\}}E_{\{a,a+r\}}
E_{A\setminus\{a\}}.
```

The last displayed pair contains one original A factor and the zero B
factor, with their signed axes and orientations. This is an analytical
reordering, not a free physical band merger. The exact paired-phase formula
makes each of the three pair restrictions a Clifford quarter-turn or a
scalar. Their core restrictions $`R_+,R_-`$ are therefore ordinary
Cliffords with literal eighth-root scalars retained.

Only the isolated A phase survives in the remaining A-sector ratio. Write
$`a_i=r+a_*`$ in these control coordinates. The exact relative block,
including its retained native scalar, is

```math
V=F_+F_-^\dagger
=R_+\,[G_+\exp(\pm i\pi\tau(a_*)/4)G_+^\dagger]R_-^\dagger D,
\qquad D=G_+G_-^\dagger.
```

D is diagonal Clifford and fixes a,b,c literally. F commutes with r, so its
Clifford image of p has blocks V and $`V^\dagger`$. Thus V is Clifford
and $`S_V^2=I`$ on the full remaining Pauli space; no radical eigenvalue
is fixed. The half-angle lemma makes
$`H=G_+\tau(a_*)G_+^\dagger`$ a signed Pauli when $`a_*\ne0`$.
It commutes with a because the original A axes do and G centralizes a.
The scalar case uses the identity transvection.

Write the binary branch actions as complementary powers of
$`T_c,T_b,T_a`$. If L is the product of the $`T_b,T_a`$ factors in the
minus branch, then

```math
S_{R_-}^{-1}S_{R_+}=L(T_cT_bT_a)L^{-1}.
```

Indeed $`T_a,T_b`$ commute and their two branch products multiply to
$`T_bT_a`$. L fixes a and h and commutes with $`S_D`$. Conjugating
$`S_V`$ first by $`S_{R_-}^{-1}`$ and then by $`L^{-1}`$, and renaming
the conjugated Pauli label H, gives

```math
S=M_0T_HS_D,\qquad M_0=T_cT_bT_a,
\qquad \langle H,a\rangle=0.
```

On the core labels,

```math
M_0a=a+c,\qquad M_0c=c+h,\qquad M_0h=h.
```

Set $`\delta=\langle H,c\rangle`$. Since $`S_D`$ fixes a,c,h,

```math
Sa=a+c,\qquad Sc=c+h+\delta M_0H,
\qquad S^2a=a+h+\delta M_0H.
```

Involutivity would force $`h=\delta M_0H`$. For $`\delta=0`$ this
contradicts $`h\ne0`$. For $`\delta=1`$ it gives
$`H=M_0^{-1}h=h`$, but then $`\delta=\langle h,c\rangle=0`$, another
contradiction. The scalar half-angle is the first case. This excludes the
last nonduplicate branch on the full remaining space.

### 5.4 Duplicate absorption finishes the mixed case

The surviving mixed case has e itself as an A axis. Its actual A/B rotations
combine into a scalar times a right Clifford $`C_e`$. It commutes with the
remaining A and all B factors. Remove both duplicate rotations and retain
the Clifford on the right. The remaining specified three-band source is
C3 after multiplication by $`C_e^\dagger`$.

Its first/middle span overlap is at most one. Those spans are subspaces of
the original ones, whose intersection was $`\mathrm{span}\{e,r\}`$;
the remaining B span does not contain e by B-axis independence. The
one-overlap theorem gives diagonal-C3 form. Restoring the actual right
Clifford and scalar gives that form for F, completing every case in the
activation table and proving the two-overlap theorem.

## 6. Six-layer consequence and physical contract

Use the actual prefix and comparison operators in the
[six-layer native reduction](READ_ONLY_SIX_LAYER_REDUCTION.md). If the
specified three-band prefix is exactly C3 and its axis data satisfy one of
the proved cases here, its fixed end Cliffords can be absorbed into the
diagonal-C3 source form. The retained
[signed transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query)
then applies to the actual comparison estimate

```math
\|\mathcal F L_x\mathcal F^\dagger-K_x\|\le2\eta.
```

At full operator error $`\eta\lt\sqrt2/8`$, the resulting necessary
Boolean-degree bounds are

```math
\deg f\le78,
\qquad
\deg f\le26\quad\text{when the outer target shear is visible}.
```

No Pauli-valued comparison, monomial inner Clifford, initialized state or
selected logical column is assumed. Exact source membership and its
structural certificate remain additional premises; returned action alone
does not imply them. The approximate query conclusion does not make the
source classification robust.

Every equality and norm estimate concerns the literal full operator,
arbitrary dirty inputs and external references, target phase, leakage and
complete workspace return. Actual adjoints are used in their displayed
order. Pauli signs remain literal even where a comparison is projective.
The coordinate choices, core identifications and fixed Clifford extensions
are analytical. They introduce no free source, inverse, helper or physical
gate. The original six physical T slots, all T/T-dagger gates, every
elementary Clifford, controls, actual inverse and simultaneous dirty-work
allocation remain charged. Finite dirty width and the nondepth counts are
unrestricted in this necessary condition. No complete-frame upper resource
bound or general six-band classification is changed.

## 7. Bounded evidence and fair literal fixtures

The [focused source tests](../../tests/test_two_overlap_source_structure.py)
check all 155 two-dimensional relation spaces at B width five, the repeated
triple obstruction, 272 full-space involution configurations including
radical couplings, four two-core sign fixtures, and 256 native orientations
of the nonduplicate zero and mixed cores. These are bounded mechanism
checks; the arbitrary-width conclusion is the analytical theorem above.

Two positive fixtures compare all full-space computational-basis columns
and the actual inverse exactly with a fair native baseline. Each table
tuple lists T-count, elementary Clifford count and T-depth in that order.
Every wire is arbitrary; no prepared-state premise is used.

| Fixture | Width | Specified source and actual inverse | Fair source | Fair actual inverse |
| --- | --- | --- | --- | --- |
| Two zero-axis duplicates | 3 | $`(6,8,3)`$ each | $`(2,8,1)`$ | $`(2,8,1)`$ |
| Mixed duplicate with a noncommuting paired core | 4 | $`(9,14,3)`$ each | $`(3,18,1)`$ | $`(3,20,1)`$ |

The mixed fixture retains a non-C3 adjacent two-band partial while its full
source is C3. Its inverse baseline includes the separate literal S-gate
corrections responsible for the different Clifford count. These are finite
source equalities, not returned-query constructions, depth-optimality claims
or asymptotic resource improvements. The displayed source and actual inverse
remain charged in any circuit that uses them.
