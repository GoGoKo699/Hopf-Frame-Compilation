# One-clean outer-overlap rigidity

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Midpoint rigidity](ONE_CLEAN_MIDPOINT_RIGIDITY.md) · [Inner spread](ONE_CLEAN_INNER_SPREAD.md) · [Research index](../README.md)

This result, R97, gives a width-independent one-clean obstruction when
the common outer phase-axis space has bounded misalignment with either
physical axis list in a native two-plus-two split. The common space
itself, both inner clean spreads and the outer crossing rank may be
arbitrarily large. The condition refers to the actual native
representation; it does not follow from terminal clean return alone.

For the complete-frame witness

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad n\ge1,
```

split a native circuit of T-depth at most four as
$`V=U_B^\dagger U_A`$, with each encoder of T-depth at most two.
After moving Clifford blocks into the axes, write each encoder as an
outer physical phase layer times a native circuit of T-depth at most one.
Let $`L_A,L_B`$ be the binary spans of those outer physical axes. Put
$`I=L_A\cap L_B`$ and let $`\mathcal U_A,\mathcal U_B`$ be the spans of the individual
physical axes in each list that themselves lie in I. Define

```math
\mu_A=\dim I-\dim\mathcal U_A,\qquad
\mu_B=\dim I-\dim\mathcal U_B,\qquad r=\min\{\mu_A,\mu_B\},
```

```math
k(r)=2r+(r+1)(2^r-1),\qquad t(r)=8r+8k(r).
```

Thus zero misalignment means that I is spanned by a subset of one
physical list. An arbitrary basis change inside I does not establish
that condition: the physical T axes must be retained.

Define the constants already used by R77:

```math
K_t=2^{t+3},\qquad a_t=2(K_t^2-1),\qquad c_t=7K_t+5t+17,
```

```math
\mathcal M_t=c_t\frac{a_t^{K_t}-1}{a_t-1},\qquad
\eta_r=2^{-\mathcal M_{t(r)}/2-2}.
```

**Outer-overlap theorem.** If $`J_1`$ initializes the one clean qubit,
then every such realization obeys

```math
\inf_{\phi\in\mathbb R}
\|VJ_1-e^{i\phi}J_1(W_n\otimes I_b)\|_{\rm op}
\ge\eta_r.
```

The norm includes all logical columns, arbitrary dirty inputs, external
references and clean leakage. Clifford interlayers, finite dirty width
and gate counts are otherwise unrestricted. The constants are
conservative; their useful property is independence of physical width
at each fixed r. When r=0, the threshold is exactly R76's
Clifford-midpoint threshold, with no bound on $`\dim I`$. This includes
$`L_A=L_B`$ even when their physical basis lists differ. Bounded
intersection dimension is an immediate sufficient condition because
$`r\le\dim I`$.

Consequently, an arbitrarily accurate four-layer one-clean sequence
would require both outer misalignments to diverge, even after minimizing
over native representations, in addition to the growing crossing rank
and both growing inner spreads required by R74 and R77. The radical of
a shared clean bilinear form must also grow, as proved in Section 9.
These necessary conditions do not exclude a sequence in which all those
quantities grow. The unrestricted one-clean five-layer obstruction and optimal general depth
remain open. The compiler bounds and selected manuscript claims are
unchanged.

## 1. Correlated native images

Suppose first that V returns the terminal clean sector exactly. With
$`\Pi=(I+Z_C)/2`$ and a data/dirty Hermitian Pauli P, the two encoders
supply

```math
S=U_AZ_CU_A^\dagger=U_BZ_CU_B^\dagger,
\qquad Q=(I+S)/2,
```

```math
A=U_APU_A^\dagger=E_AK_AE_A^\dagger,\qquad
B=U_BPU_B^\dagger=E_BK_BE_B^\dagger,
```

```math
S=E_AJ_AE_A^\dagger=E_BJE_B^\dagger.
```

The four cores are Hermitian Cliffords, A,B commute with S, and EA,EB
are Clifford conjugates of physical T/T-dagger layers. Each layer's
active Pauli axes are independent and isotropic. All conjugation
identities below retain Hermitian signs. Binary labels omit those signs
only when discussing symplectic maps.

We first prove an asymmetric compressed-image theorem: EA may be any
phase-polynomial layer on LA, while EB remains a physical phase layer
with independent axes. The cores retain their displayed Hermitian
Clifford properties. For $`d=\dim(L_A\cap L_B)`$, we prove

```math
\|(A-B)Q\|\lt\kappa_d
\quad\Longrightarrow\quad (A-B)Q=0,
\qquad \kappa_d=2^{-\mathcal M_{t(d)}/2-1}.             \qquad(1)
```

The proof first finds a common change of midpoint coordinates in which
S has a bounded finite Clifford expansion. The transformed outer
operations remain phase-polynomial layers. A phase-polynomial layer
means a Clifford-conjugated diagonal weighted-cubic phase operation,
as defined in [R76, Section 1](ONE_CLEAN_MIDPOINT_RIGIDITY.md#1-phase-polynomial-layers-and-operator-rigidity).
Its displayed axes may be dependent. This analytical change of
coordinates makes no new physical T-depth-one claim.

## 2. A physical-axis normalizer

**Normalizer lemma.** Let E be a product of T/T-dagger rotations about
independent commuting Hermitian Pauli axes. For a Clifford K,
$`EKE^\dagger`$ is Clifford if and only if K signed-permutes those
individual active axes.

For the forward direction, write $`K'=EKE^\dagger`$. For each active
axis a, choose a Pauli R anticommuting with a and commuting with all
other active axes. Independence and nondegeneracy of the ambient
symplectic form guarantee this choice. The identity

```math
K(E^\dagger RE)K^\dagger
=E^\dagger(K'RK'^\dagger)E
```

has exactly two Pauli terms on the left, with label difference
$`F_Ka`$. A physical phase layer conjugates a Pauli anticommuting
with s active axes into exactly $`2^s`$ distinct nonzero Pauli terms.
Thus the right side has s=1, and its label difference is an individual
active axis. K permutes the finite independent list, with possible
Hermitian signs.

Conversely, signed permutation makes
$`EKE^\dagger K^\dagger`$ a product of even T-angle rotations
on that same isotropic axis space. It is a diagonal Clifford.
The lemma uses physical independence; it is not asserted for an
arbitrary phase-polynomial layer with dependent displayed axes.

## 3. Shared-image compatibility and exceptional axes

Write B's outer axes as $`a_1,\ldots,a_m`$ and put

```math
I=L_A\cap L_B,\qquad c_i=Ja_iJ,\qquad
C_{ij}=\langle c_i,a_j\rangle.
```

The binary matrix C is symmetric because J is a symplectic involution.
For $`R\in L_A^\perp`$, set

```math
r_i=\langle a_i,R\rangle,\qquad
K_R=J(E_B^\dagger RE_B)J,\qquad G=SRS.
```

Both KR and G are Clifford. On the A side, R commutes with EA, so
$`G=E_A(J_ARJ_A)E_A^\dagger`$ is a one-layer Pauli image and

```math
\mathrm{im}(F_G-I)\subseteq L_A.
```

On the B side $`G=E_BK_RE_B^\dagger`$. The normalizer lemma makes
KR a signed permutation of B's physical axes. Its binary displacement is

```math
N_r:=F_{K_R}-I=\sum_i r_i c_ic_i^{\mathsf T}\Omega,\qquad
N_rN_s=0,
```

where Omega is the ambient symplectic matrix. The allowed activation
vectors form

```math
H=I^\perp\subseteq\mathbb F_2^m,
```

where I is expressed in the independent B-axis coordinates and this
perpendicular uses their ordinary coordinate pairing.

### A common matching

The induced permutations of the B axes form an elementary abelian
two-group, and their matrices satisfy
$`(P_r-I)(P_s-I)=0`$. Every orbit has size at most two. Otherwise two
independent moves on an axis would produce four distinct axes whose
sum is $`N_rN_sa=0`$, contradicting physical independence. Thus there
is one matching of B's axis list: every allowed r fixes or swaps each
matched pair.

For a fixed column j,

```math
N_r(a_j)=\sum_i r_iC_{ij}c_i
```

is either zero or the fixed difference of its matched pair. Independence
of the ci implies that all nonzero coordinate characters
$`r\mapsto r_i`$ appearing in that column coincide. Coordinates whose
character vanishes on H are not constrained by this statement.

Because $`N_r(L_B)\subseteq L_B`$, pairing its images with B axes gives

```math
C\mathrm{diag}(r)C=0\qquad(r\in H).
```

The diagonal entries, using symmetry of C, show that each column of C
annihilates H and hence lies in I. Therefore the single nonzero
activation character in a column must occur an even number of times.
A singleton nonzero-character coordinate cannot occur in any column.

Let X contain every zero-character coordinate and every coordinate in
a repeated nonzero-character class. If there are z zero coordinates and
the repeated classes have sizes $`q_1,q_2,\ldots`$, their independent
relations in I give

```math
z+\sum_j(q_j-1)\le d,\qquad
|X|=z+\sum_jq_j\le2d.
```

The rows of C outside X are zero, and symmetry gives zero columns too.
We call indices outside X regular. For every regular i, ci commutes
with all of LB.

### Containment in the two outer spaces

Every regular ci also belongs to $`L_A+L_B`$. To see this, the normalizer
lemma gives a diagonal Clifford

```math
D_R=E_BK_RE_B^\dagger K_R^\dagger
```

on LB. Since $`G=D_RK_R`$, the difference between the symplectic maps
of G and KR has image in LB. Projecting modulo $`L_A+L_B`$ therefore
annihilates

```math
F_{K_R}-I=\sum_i r_i c_ic_i^{\mathsf T}\Omega.
```

The ci are independent, so a dual Pauli test vector isolates every
activated ci. Each regular coordinate has a nonzero character on H
and is activated by some R. This proves the containment.

## 4. Only boundedly many regular images are mixed

Choose a d-by-m binary matrix whose row space is I in the B-axis
coordinates. Its columns have at most $`2^d`$ patterns. These patterns
are distinct from the activation-character classes used in Section 3.
Equal column patterns at i,j give $`e_i+e_j\in H`$.

Take regular i,j with equal patterns, and choose R activating precisely
those two indices. Since ci,cj commute with LB, the symmetric
displacement tensor of G is

```math
(F_G-I)\Omega^{-1}
=c_i\otimes c_i+c_j\otimes c_j+
\sum_l u_l a_l\otimes a_l,\qquad u_l=\langle c_l,R\rangle.
```

Its image is contained in LA. Because this tensor is symmetric, both
factors lie in LA. In the direct decomposition

```math
(L_A+L_B)/I=(L_A/I)\oplus(L_B/I),\qquad
\overline c_i=x_i+y_i,
```

the cross block therefore gives

```math
x_i\otimes y_i=x_j\otimes y_j.
```

If the common rank-one tensor is nonzero, its factors over the binary
field are unique. All ci in this pattern class then have the same image
modulo I. Their independence limits the class to at most d+1 members.
If the tensor is zero, every member lies wholly in LA or wholly in LB.
The zero column pattern admits the singleton activation ei; the same
argument directly gives $`x_i\otimes y_i=0`$ there. Consequently at
most

```math
(d+1)(2^d-1)
```

regular indices have mixed images, meaning ci lies in neither LA nor LB.

Partition B's actual physical axis list into three disjoint sets:

- Z contains X and all mixed regular indices.
- T contains the remaining indices with $`c_i\in L_A\setminus L_B`$.
- P contains the remaining indices with $`c_i\in L_B`$, including
  images lying in I.

Write their spans as $`L_Z,L_T,L_P`$, and put $`k=|Z|`$. Then

```math
k\le k(d)=2d+(d+1)(2^d-1),\qquad
L_B=L_P\oplus L_T\oplus L_Z.
```

For i in T, regularity gives $`Sa_iS=c_i`$, hence $`Sc_iS=a_i`$.
Since ci lies in LA, the latter identity is

```math
E_A(J_Ac_iJ_A)E_A^\dagger=a_i.
```

Here EA is allowed to be phase-polynomial, so one cannot infer that
its input equals its Pauli output. The following simultaneous extraction
supplies the exact replacement needed for all transfer axes.

## 5. Simultaneous extraction and the transfer gauge

### A diagonal Clifford matching all Pauli outputs

Let D be a phase-polynomial unitary on an isotropic space L. Suppose it
maps a commuting independent Pauli family $`p_i`$ to a Pauli family
$`a_i`$, including their literal signs. There is a diagonal Clifford C
on L with the same action on their whole generated group. Moreover
$`D_0=DC^\dagger`$ commutes with both families, and its phase support
can be chosen inside $`L\cap\mathrm{span}\{a_i\}^\perp`$.

To prove this, diagonalize L to k Z coordinates. D preserves the
X-translation v on those coordinates. Its Z-label shift Bv depends
only on v: two inputs with the same translation differ by a Pauli
commuting with D. The translations form a subspace T0. Preservation of
commutation gives

```math
v\cdot Bw=w\cdot Bv\qquad(v,w\in T_0).
```

Extend this partial symmetric map to a symmetric k-by-k matrix: choose
a basis with T0 first, use the prescribed symmetric block and cross
block, use its transpose for the opposite block, and choose the
remaining block arbitrarily. A diagonal Clifford with this shear has
the correct unsigned Pauli action. The residual literal signs form a
linear character of T0 and are corrected by a Z Pauli. Thus C exists.

Since D and C commute, D0 fixes both the inputs and outputs. Its diagonal
phase is invariant under the translation subspace T0. In binary
coordinates adapted to T0 it is independent of those coordinates.
Removing them preserves the weighted-cubic phase-polynomial form and
leaves support in the annihilator of T0, as claimed.

### Applying the extraction

Apply the lemma to EA and the transfer family
$`p_i=J_Ac_iJ_A`$, whose exact Pauli outputs are ai by Section 4.
Write $`E_A=D_0C`$ and absorb C into both A Clifford cores. Let LA0
be the support of D0, chosen so that

```math
L_{A,0}\subseteq L_A\cap L_T^\perp.
```

Factor B's physical layer as $`E_B=E_TE_PE_Z`$ and left-multiply both
encoders by $`E_T^\dagger`$. V and the occupied discrepancy norm
are unchanged. A's new outer operation is $`E_T^\dagger D_0`$, one
phase-polynomial layer on $`L_T+L_{A,0}`$; the summands commute but
may be dependent. B's new outer layer is $`E_PE_Z`$. Both clean/data
cores on each side are still Hermitian Cliffords. This is an analytical
change of midpoint coordinates and does not assert physical independence
for the new A outer word.

Use primes for these new midpoint images. Set

```math
M=L_P+JL_P.
```

Both summands lie in LB, so M is isotropic; J preserves it because
$`J^2=I`$. If $`h_0=\dim M-|P|`$, then

```math
h_0\le d+k.                                         \qquad(2)
```

Indeed, the independent vectors in JT lie in LA, so their rank modulo
LB is at least $`|T|-d`$. Hence

```math
\dim(L_B\cap JL_B)\le m-|T|+d.
```

M lies in that intersection, while $`|P|=m-|T|-k`$, proving (2).

Strip EZ for analysis and define

```math
S_0=E_PJE_P^\dagger=D_PJ,\qquad
D_P=E_P(JE_P^\dagger J),\qquad
S'=E_ZS_0E_Z^\dagger.
```

DP is a diagonal phase-polynomial operation on M: its two axis lists
lie in that same isotropic space. Thus S0 is a third-level
Clifford-hierarchy operation, and $`G_0=S_0RS_0`$ is Clifford for
every Pauli R. This assertion concerns an exact analytical factorization,
not a free implementation of DP.

## 6. Bounded central sectors give Clifford restrictions

Let

```math
L_*=L_A+L_T+L_Z,\qquad K=M\cap L_*.
```

Lstar need not be isotropic. The dimension of K obeys

```math
\dim K\le d+h_0\le2d+k.                             \qquad(3)
```

For a direct verification, M lies in LB and
$`L_*\cap L_B=I+L_T+L_Z`$. Project LB onto $`L_B/L_P`$.
The image of M has dimension h0. The kernel of this projection on K
has dimension at most d, because LP is complementary to
$`L_T+L_Z`$ and only I can contribute to that kernel. This proves (3).

Take any $`R\in L_*^\perp`$. It commutes with A's new outer
phase-polynomial layer and with EZ. Hence

```math
G'=S'RS'
```

is Clifford with displacement image in $`L_A+L_T`$. Also
$`G'=E_ZG_0E_Z^\dagger`$. Apply the physical-axis normalizer lemma
to the size-k sublayer EZ and the Clifford G0. The discrepancy
$`E_ZG_0E_Z^\dagger G_0^\dagger`$ is a diagonal Clifford on LZ.
It follows that $`\mathrm{im}(F_{G_0}-I)\subseteq L_*`$.
The expression

```math
G_0=D_P(JRJ)D_P^\dagger
```

also puts this image in M. Therefore

```math
\mathrm{im}(F_{G_0}-I)\subseteq K.              \qquad(4)
```

As R ranges over $`L_*^\perp`$, the translations of JRJ on the
commuting M coordinates fill the annihilator of JK in the dual of M.
Indeed, an element of M annihilates all those translations exactly when
its J-image lies in Lstar; because JM=M, this annihilator is JK.

### The odd-axis tensor

For a diagonal phase-polynomial word on M, define its binary symmetric
tensor by retaining the odd T coefficients:

```math
\Theta=\sum_v e_v\,v\otimes v\otimes v,
\qquad e_v\in\mathbb F_2.
```

Contracting the first factor against a Pauli translation gives the
symmetric displacement matrix of its Clifford conjugate. Equation (4)
therefore says that, for every translation s annihilating JK,

```math
\Theta(s,-,-)\in K\otimes K.                         \qquad(5)
```

Image containment gives one factor in K; symmetry gives the other.
Put

```math
A_0=K+JK,\qquad h=\dim A_0\le2\dim K\le4d+2k.
```

The projection of Theta to $`(M/A_0)^{\otimes3}`$ vanishes. All three
dual arguments there annihilate A0; the first annihilates JK, and the
other two annihilate K in (5).

For completeness, zero tensor is exactly the diagonal Clifford
condition. Write a parity phase in binary coordinates as

```math
p(x)=\sum_v\alpha_v(v\cdot x\bmod2)\pmod8.
```

Its linear, quadratic and cubic monomial coefficients are respectively

```math
\ell_i=\sum_v\alpha_vv_i,\qquad
q_{ij}=-2\sum_v\alpha_vv_iv_j,\qquad
r_{ijk}=4\sum_v\alpha_vv_iv_jv_k.
```

It is Clifford exactly when the linear coefficients are even, the
quadratic coefficients are divisible by four, and the cubic coefficients
vanish modulo eight. These are precisely the zero entries
$`\Theta_{iii},\Theta_{iij},\Theta_{ijk}`$.

### Fixing the central characters

Choose Clifford coordinates diagonalizing M, with its first h generators
spanning A0, and fix any character z of those h commuting Paulis.
A fixed character changes the remaining parity coefficients only by
signs and adds a constant. Its odd-axis tensor is the projection just
shown to vanish. The restriction of DP is therefore a diagonal
Clifford $`C_z`$ on the remaining coordinates.

Its actual scalar is an eighth root. J normalizes M and hence permutes
its joint characters, including Pauli signs. DP is the ratio of the
corresponding EP eigenvalues; the common scalar of EP cancels. Thus DP
has eighth-root entries, including after fixing z. Extend each Cz to a
Clifford on all coordinates with that same scalar. Then

```math
D_P=\sum_z\Pi_zC_z,\qquad
\Pi_z=2^{-h}\sum_{a\in A_0}\chi_z(a)a.
```

There are $`2^{2h}`$ Clifford summands, each with coefficient magnitude
$`2^{-h}`$. Multiplication by J gives the same expansion for S0.
The data images need not preserve these sectors: this decomposition is
an operator identity for the midpoint, used only as a finite boundary
expansion.

## 7. A finite Clifford expansion and the occupied gap

For a Hermitian Pauli a, use the phase convention

```math
T_a^\epsilon=\frac{1+\zeta_8^\epsilon}{2}I
             +\frac{1-\zeta_8^\epsilon}{2}a,
\qquad \epsilon\in\{1,-1\},\qquad \zeta_8=e^{i\pi/4}.
```

Expanding each numerator into two terms expresses one such gate as
four Pauli terms of coefficient magnitude one half, with eighth-root
phases. Conjugation by the size-k layer EZ consequently expresses each
Clifford as $`16^k`$ Clifford terms of magnitude $`4^{-k}`$.
Rotation-convention common phases cancel in conjugation.

Combining this with Section 6 gives a flat, possibly repeated expansion

```math
S'=2^{-t_0/2}\sum_{j=1}^{2^{t_0}}\sigma_jC_j,
\qquad t_0=2h+4k\le8d+8k\le t(d),                  \qquad(6)
```

where eighth roots are absorbed into fixed ring-valued Clifford
representatives and $`\sigma_j\in\{1,-1\}`$. Both t0 and t(d)
are even. To use exactly t(d), duplicate each term
$`2^{(t(d)-t_0)/2}`$ times and lower its coefficient to
$`2^{-t(d)/2}`$, then add cancelling Clifford pairs until there are
$`2^{t(d)}`$ terms. The deficit is nonnegative and even; the zero case
requires no padding.

Apply the compressed-discrepancy argument of
[R77, Section 2](ONE_CLEAN_INNER_SPREAD.md#2-compressed-image-theorem)
to the transformed pair, using (6). For every Pauli R commuting with
A's transformed outer axis space, both
$`G=A'RA'`$ and $`S'GS'`$ are Clifford, because each is a
phase-polynomial conjugate of a Pauli. The other side has four
Clifford boundaries sharing B's single outer phase layer. Thus the
normalized Choi discrepancy has the
[R74 finite-boundary form](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
with

```math
K=K_t=2^{t+3},\qquad b=t+4,\qquad m=0,
\qquad t=t(d).
```

The proof uses only the finite Clifford expansion in (6), not a bound
on either original inner spread. It gives

```math
Q'(A'RA'-B'RB')Q'=0
```

or normalized Hilbert--Schmidt squared norm at least
$`2^{-\mathcal M_t}`$. If $`\varepsilon=\|(A'-B')Q'\|`$,
the discrepancy's operator norm is at most $`2\varepsilon`$.
Below $`\kappa_d`$ every such compressed conjugation is therefore frozen.
No union bound over Paulis is required.

The [R76 compressed commutant lift](ONE_CLEAN_MIDPOINT_RIGIDITY.md#lifting-the-compressed-commutant)
then makes $`B'A'Q'=D_0Q'=Q'D_0`$, where D0 is a unitary in the
abelian Pauli algebra of A's transformed outer layer. In its diagonal
coordinates, that phase-polynomial layer has eighth-root entries up to
one common scalar. Every nonzero entry of $`S'`$ is consequently a
flat Clifford entry times an eighth root, hence a unit of
$`\mathbb Z[i,1/\sqrt2]`$. The
[R77 occupied-phase argument](ONE_CLEAN_INNER_SPREAD.md#3-freezing-and-occupied-phases)
places each occupied scalar of D0 in that ring and makes it an eighth
root. A nonidentity scalar would have distance at least
$`\sqrt{2-\sqrt2}`$ from one, contradicting the smaller threshold.
Thus $`(A'-B')Q'=0`$. Undoing the common midpoint change proves (1).

## 8. Native misalignment and the frame witness

Return to two native outer layers and their original intersection I.
Suppose first that the A list has $`\dim\mathcal U_A`$ individual axes in I.
Left-multiply both encoders by the inverse of that A sublayer. A's
remaining outer layer is still physical, with its coordinate complement
$`L_A'`$ satisfying

```math
\dim(L_A'\cap L_B)=\dim I-\dim\mathcal U_A=\mu_A.
```

B's new outer operation is a phase-polynomial layer on LB, because every
removed A axis lies in LB. Its actual phase support may shrink; using
LB as a containing isotropic space is sufficient. Apply the asymmetric
compressed theorem with the two sides exchanged. It gives (1) with
$`d=\mu_A`$. The same reduction using B's individually contained
axes gives $`d=\mu_B`$. We may therefore take
$`r=\min\{\mu_A,\mu_B\}`$.

In particular, when I is spanned by a subset of either physical list,
r=0 and the original R76 threshold applies. No common list or matching
basis is required on the other side. Equal outer spaces are included.

Suppose now that V approximates the initialized witness with error eta.
[R75](FOUR_LAYER_OPERATOR_RIGIDITY.md) first forces exact terminal clean
return when

```math
\eta\lt\tau=\frac{\sin(\pi/8)}{2\sqrt2}.
```

For every target-preserved data/dirty Pauli P, the encoders then have
the common clean image used above and obey

```math
\|(A-B)Q\|\le2\eta.
```

Since $`\eta_r=\kappa_r/2\lt\tau`$, error below $`\eta_r`$ freezes
all spectator and dirty X,Z observables and the active Y observable.
The [R76 witness consequence](ONE_CLEAN_MIDPOINT_RIGIDITY.md#4-initialized-error-and-the-frame-witness)
forces the occupied circuit to be a one-qubit rotation by an integer
multiple of $`\pi/8`$, modulo a common phase, tensored with the
spectator/dirty identity. Its distance from the target is at least

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\gt\tau,
```

the contradiction proving the theorem.

## 9. An intrinsic clean bilinear form

Before the analytical gauges, define a binary form on the native common
outer space I by

```math
\beta(i,j)=\langle i,J_Aj\rangle.
```

It is symmetric because JA is a symplectic involution. It has the same
value using J in place of JA. Indeed, for $`j\in I`$, the Clifford
$`G_j=SjS`$ is represented on either side by the corresponding
one-layer Pauli image. Its conjugation sign on $`i\in I`$ is
$`(-1)^{\langle i,J_Aj\rangle}`$ on A and
$`(-1)^{\langle i,Jj\rangle}`$ on B. These signs are equal. The form is therefore intrinsic to S on this
specified I; the space I itself still depends on the native representation.

Let FA span A's physical axes outside I. The essential-axis span of the
Clifford $`G_j`$ is intrinsic: a physical phase layer conjugating a Pauli
has symplectic displacement image exactly the span of its activated
axes, by [R76's essential-axis lemma](ONE_CLEAN_MIDPOINT_RIGIDITY.md#the-essential-axes-are-intrinsic).
The A and B spans must agree and hence lie in I. No A axis outside I
can be activated, so

```math
F_A\perp J_AI.
```

Physical independence gives

```math
I=\mathcal U_A\oplus(I\cap F_A),\qquad
I\cap F_A\subseteq\mathrm{rad}\beta.
```

Consequently

```math
\mu_A\le\dim\mathrm{rad}\beta,\qquad
\mu_B\le\dim\mathrm{rad}\beta.
```

A bounded radical therefore also supplies a positive width-independent
floor. A nondegenerate beta forces alignment with both physical lists
and gives the uniform R76 threshold, even at unbounded intersection
dimension. Nondegeneracy is sufficient, not necessary: the aligned
case of Section 8 permits degenerate beta as well.

## 10. Scope and the reopening condition

The misalignments concern the given native representation. A common
symplectic coordinate change preserves them, but a different native
word or split may change them. Since the theorem applies to every such
representation, a sequence whose errors tend to zero must have both
misalignments tending to infinity even after minimizing over native
representations. Its common-intersection dimension and the radical
dimension from Section 9 must then diverge as well. These invariants
are necessary conditions, not a characterization of approximation power.

The unresolved step is precise: the finite expansion and its gap constant
grow with the unaligned part of the common outer space. This proof does
not supply a constant independent of that misalignment. Reopening the
unrestricted route requires an operator-norm argument controlling those
correlated sectors uniformly, a different normal form removing that
dependence, or an actual native family with vanishing occupied defect.
Increasing finite test widths or merely renaming the unbounded space
does not resolve the step. The spectator-padded witness has constant
depth at fixed accuracy; a dimension-growing lower bound still needs
another witness. Optimal general depth remains the destination.
