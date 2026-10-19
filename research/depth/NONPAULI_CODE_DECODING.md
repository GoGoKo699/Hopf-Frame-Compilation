# Exact decoding of a one-band non-Pauli clean code

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md) ·
[Occupied logical transport](OCCUPIED_LOGICAL_TRANSPORT.md)

This component gives an exact controlled-Clifford description of a
one-band clean image. Its controls have bounded Boolean degree, even
when the intrinsic active rank is unbounded. The decoder is generally
not a full Clifford, and inserting it does not preserve a three-band
budget. The growing-rank example identifies this interface limitation explicitly.

## 1. Signed native coordinates

Let P be a signed Hermitian Pauli and let

```math
E=\prod_i T_{A_i}^{\epsilon_i},\qquad
T_A=\frac{I+A}{2}+\omega\frac{I-A}{2},\qquad
\omega=e^{i\pi/4},\qquad \epsilon_i\in\{1,-1\},
```

be a native band on independent commuting physical Pauli axes. Suppose
exactly k axes anticommute with P, with k positive. The silent factors
commute with P and with the active factors, so
$`S=EPE^\dagger=E_aPE_a^\dagger`$.

Literal inverse orientations can be removed from the active list using

```math
T_A^\dagger=\omega^{-1}T_{-A}.
```

Let n_minus be the number of negative active orientations and put
$`\widetilde A_i=\epsilon_i A_i`$. Then

```math
E_a=\omega^{-n_{\mathrm{minus}}}\prod_{i=1}^kT_{\widetilde A_i}.
                                                               \qquad(1)
```

The scalar in (1) is retained; it cancels only in the conjugation
defining S. Choose a fixed Clifford B sending the actual signed Paulis to

```math
BPB^\dagger=X_t,\qquad
B\widetilde A_1B^\dagger=Z_t,\qquad
B\widetilde A_{i+1}B^\dagger=Z_tZ_i\quad(1\le i\lt k).
                                                               \qquad(2)
```

Such coordinates exist: the commuting independent differences
$`\widetilde A_1\widetilde A_{i+1}`$ commute with P, and together
with the anticommuting pair $`(\widetilde A_1,P)`$ extend to a signed
Pauli frame. All other wires are spectators. The cost of the actual
Clifford word B and its actual inverse is paid separately below.

Write m=k-1, let z be the m central computational bits, and put

```math
w(z)=\sum_{i=1}^m z_i,\qquad e=k\bmod2,\qquad
r(z)=\frac{k-e}{2}-w(z)\pmod4.
```

In these coordinates the positive-orientation active band has the exact
fiber identity

```math
E_+(z)=\omega^{w(z)}S_t^{r(z)}T_t^e,
\qquad BE_aB^\dagger=\omega^{-n_{\mathrm{minus}}}E_+.
                                                               \qquad(3)
```

Indeed its phase on a core bit t is
$`t+\sum_i(t\mathbin\oplus z_i)=w(z)+(k-2w(z))t`$ modulo eight.

## 2. Four possible core directions and low-degree controls

Define the controlled core Clifford

```math
D=\sum_z|z\rangle\langle z|\otimes S_t^{r(z)},\qquad
Q_e=T_t^eX_tT_t^{-e}.
```

Then the full physical clean observable satisfies

```math
BSB^\dagger=DQ_eD^\dagger.                         \qquad(4)
```

For even k, Q_e is X. For odd k, it is
$`(X+Y)/\sqrt2`$. Thus only four equatorial core directions occur,
with their parity class fixed by k. They are not eight independent
choices. Both eigenspaces, including the unused clean sector, are
covered by (4).

The low binary digit of r is affine in z and its high digit has Boolean
degree at most two, by the binary-carry formula for an integer affine
sum modulo four. D is therefore a Clifford on the core for every
central character, controlled by two bounded-degree bits. Globally its
phase is

```math
D|t,z\rangle=\omega^{2t r(z)}|t,z\rangle.            \qquad(5)
```

This is a diagonal C3 phase polynomial. It is generally not a full-space
Clifford: already k=2 gives $`D=S_t\,\mathrm{CS}_{t,1}^\dagger`$.
The central scalar in (3) has not been discarded when describing the
encoder; it merely does not affect the observable in (4).

For an alternative coordinate convention retaining the original signed
axes as positive $`Z_t,Z_tZ_i`$ while retaining orientations
$`\epsilon_i`$, replace w by
$`b(z)=\sum_{i=1}^m\epsilon_{i+1}z_i`$ and k by
$`c=\sum_{i=1}^k\epsilon_i`$ in (3):

```math
E_a(z)=\omega^{b(z)}S_t^{(c-e)/2-b(z)}T_t^e.
                                                               \qquad(6)
```

Here $`c\bmod2=e`$. Equation (6) is literal and gives the same
degree bounds without a projective identification.

## 3. Native forward and inverse costs

Let

```math
C_X=\prod_{i=1}^m\mathrm{CX}_{t,i}.
```

The shared-control CNOTs commute; they are ordinary Clifford gates and
do not require extra workspace. A literal native realization of D is

```math
D=C_X\left(T_t^{1-e}\prod_{i=1}^mT_i\right)
       C_X^\dagger\left(\prod_{i=1}^mT_i^\dagger\right).
                                                               \qquad(7)
```

Within each displayed T band every T gate acts on a distinct physical
wire. The phase of (7) is
$`(1-e)t+\sum_i(t\mathbin\oplus z_i)-\sum_i z_i=2tr(z)`$
modulo eight. This proves the scalar as well as the conjugation action.

| Canonical operation | T count | Native T-depth | H/S/CNOT count | Added work |
| --- | ---: | ---: | ---: | ---: |
| D | $`2k-1-e`$ | at most 2 | $`2(k-1)`$ | 0 |
| Actual inverse of D | $`2k-1-e`$ | at most 2 | $`2(k-1)`$ | 0 |
| $`L=H_tE_+^\dagger`$ | k | at most 1 | $`2(k-1)+1`$ | 0 |
| Actual inverse of L | k | at most 1 | $`2(k-1)+1`$ | 0 |

For k=1, D is identity and its depth is zero. The inverse rows use the
actual reversed words, with every T orientation reversed. Only H and
CNOT appear in these canonical Clifford words, so their inverse counts
are unchanged. All identities hold on arbitrary states and external
references; no initialized or borrowed helper is used.

The cheaper complete code decoder L obeys

```math
L(BSB^\dagger)L^\dagger=Z_t.
```

It differs from $`H_tT_t^{-e}D^\dagger`$ by the known central phase
$`\omega^{-w(z)}`$. Thus it decodes the same code but is not the same
full isometry. This distinction matters when comparing encoders.

For the physical conjugated operation $`B^\dagger DB`$, add
$`G(B)+G(B^\dagger)`$ to the Clifford column. Its actual inverse has
the same sum. For a decoder from physical coordinates to canonical clean
coordinates, $`LB`$ has Clifford cost
$`G(B)+2(k-1)+1`$; its inverse has cost
$`G(B^\dagger)+2(k-1)+1`$. No equality of the costs of a general
Clifford word and its inverse is assumed. Fixed Clifford coordinate
changes do not increase native T-depth.

## 4. The decoder does not preserve a bounded quantum core at an interface

Having two low-degree control bits does not make conjugation by the
decoder a Clifford change of Pauli coordinates. For an explicit family,
take

```math
D_m=\prod_{i=1}^m\mathrm{CS}_{t,i}
    =\sum_z|z\rangle\langle z|\otimes S_t^{\sum_i z_i}.
```

Again the two control digits have degrees at most one and two. Yet

```math
D_mX_tD_m^\dagger
=X_t\left(\prod_iS_i\right)
       \left(\prod_i\mathrm{CZ}_{t,i}\right)        \qquad(8)
```

is a Hermitian Clifford whose symplectic displacement rank grows with m.
The diagonal Clifford on the right of (8) has displacement matrix

```math
\mathcal B_m=
\begin{pmatrix}0&\mathbf1^T\\\mathbf1&I_m\end{pmatrix},
\qquad
\mathrm{rank}\mathcal B_m=
\begin{cases}m,&m\text{ even},\\m+1,&m\text{ odd}.
\end{cases}                                                    \qquad(9)
```

To compute its kernel, the lower rows force every central component
to equal the core component; the top row then imposes m times that
component equal to zero. The Pauli X factor has trivial symplectic
action, so it does not change this rank.

Consequently a native Pauli-axis band conjugated through D need not
remain a Pauli-axis band. The first unsupported interface inference
would be to apply native axis matching, or the three-layer degree theorem,
after this decoding while retaining the original three native bands.
Equation (7) supplies an actual implementation, but it has its own band
cost. Equations (8)--(9) separately show why low-degree decoder controls
alone do not supply a bounded-rank conjugation argument. These are
limitations of that proposed reduction, not a high-degree Boolean-query
counterexample.

The [native fixtures](../../tests/test_nonpauli_code_decoding.py) check
small full-space decoding and inverse identities and recover displacement
rank from physical Pauli images. The algebra above establishes the
width-independent statements.

## 5. The zero-pair, zero-character logical defect is abelian

This section has additional joint hypotheses. Use the native logical
probe notation from [clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md#3-native-probes-of-a-compressed-logical-image):

```math
A_z=E_BJ_zE_B^\dagger,\qquad A_z\Pi=C_z\Pi,\qquad
\Pi=(I+S)/2,\qquad W_z=C_zA_z.                    \qquad(10)
```

Here z ranges over the complete logical computational-Z label group;
the A_z, J_z and C_z are exact commuting Hermitian representations,
J_z and C_z are actual full-space Cliffords, and A_z and C_z commute
with the shared non-Pauli clean image S.
Let L_B be the actual independent middle-band span, let L_S be the
intrinsic active span of S, and let its full Pauli support be p+L_S.
Assume the clean matching has no middle pairs, so m=0. Assume also
that the active rank of S is not two and that **every** logical probe
character nu_z from the cited theorem is zero. The last condition is
not implied by m=0.

The support and retained-factor lemmas give the same Pauli space for
the entire family:

```math
K=L_B+L_S+\langle p\rangle,\qquad
\mathrm{supp}_P(W_z)\subseteq K,\qquad
C_z(K)=K.                                        \qquad(11)
```

To obtain the normalization in (11), the retained-factor lemma gives
$`C_z(L_B)\subseteq K`$, and $`[C_z,S]=0`$ preserves the full support
coset p+L_S and hence its linear span. Invertibility then gives equality.
Since m=0, L_B+L_S is isotropic. The vector p pairs nontrivially with
L_S because S is a non-Pauli one-band clean image. Thus

```math
\mathrm{rank}_{\mathrm{symp}}K=2.                   \qquad(12)
```

Choose one fixed Pauli coordinate system for K. Its algebra has central
computational bits u, one quantum core qubit t, and identity spectators.
The u variables are intrinsic central coordinates and are not presumed
to be the target's computational input bits. Every Pauli in p+L_S
has nonidentity core support: choose a clean-active a in L_S with
$`\langle a,p\rangle=1`$; it anticommutes every such Pauli. Therefore
every fiber S(u) is a traceless Hermitian involution. It has one-dimensional
positive and negative eigenspaces, with projectors Pi_+(u), Pi_-(u).

By (10), W_z fixes the positive eigenspace pointwise and commutes S.
There is consequently a unit-modulus scalar alpha_z(u) with the exact
fiber identity

```math
W_z(u)=\Pi_+(u)+\alpha_z(u)\Pi_-(u).               \qquad(13)
```

In particular, all W_z commute on the whole space, including the unused
clean sector. This is an abelian defect statement, not a claim that
each W_z is a Clifford.

Let $`\beta_v(X)=C_vXC_v`$. Because $`A_z=C_zW_z`$ and both the
A and C families are exact commuting representations, multiplication
gives

```math
W_{z+v}=\beta_v(W_z)W_v=\beta_z(W_v)W_z.            \qquad(14)
```

Each C_v induces an affine permutation sigma_v of the central characters
u and preserves S. Equation (14) is thus the scalar cocycle identity

```math
\alpha_{z+v}(u)=\alpha_z(\sigma_vu)\alpha_v(u).     \qquad(15)
```

### 5.1. A common determinant coboundary constrains the phases

The actual native factors give more than commutation. On the same K,
the retained factorization is

```math
W_z=(C_zE_BC_z^\dagger)(C_zJ_z)E_B^\dagger.        \qquad(16)
```

Both native bands and the middle Clifford belong to the algebra of K.
Define the one-core determinants

```math
\delta(u)=\det E_B(u),\qquad
\kappa_z(u)=\det(C_zJ_z)(u).
```

For a native axis with nonidentity core Pauli, its two eigenvalues are
positive and negative one in every central fiber. Its literal T or
T-dagger factor therefore has determinant $`\omega^{\epsilon_i}`$,
independent of u. A central-only signed axis with eigenvalue
$`(-1)^{\ell_i(u)}`$ contributes determinant
$`i^{\epsilon_i\ell_i(u)}`$, where ell_i is affine. Hence

```math
\delta(u)=
\omega^{\sum_{i\;\mathrm{noncentral}}\epsilon_i}
i^{\sum_{i\;\mathrm{central}}\epsilon_i\ell_i(u)}. \qquad(17)
```

This is a central diagonal Clifford phase, up to its displayed fixed
eighth-root scalar. In particular, expanding affine Boolean parities
modulo four gives linear terms plus twice quadratic terms.

A full Clifford inside this one-core algebra has blocks

```math
(C_zJ_z)(u)=\gamma_z i^{a_z(u)}
             X^{v_z(u)}Z^{w_z(u)}H_{z,*},
```

where H_{z,*} is a fixed one-qubit Clifford, v_z,w_z are affine, and
a_z is a central Clifford exponent modulo four, in ordered X-then-Z
convention. Taking a two-dimensional determinant doubles a_z, removing
its quadratic terms. It follows that

```math
\kappa_z(u)=\xi_z(-1)^{b_z(u)},\qquad
\xi_z\in\{1,-1,i,-i\},\qquad b_z\text{ affine}.  \qquad(18)
```

The root phases in (18) are those of the actual native Clifford, not
independently chosen projective representatives. Since C_z normalizes
K, its conjugation permutes the central fibers and conjugates their core
matrices. Determinants are unchanged by the core conjugation. Finally,
(13) gives $`\det W_z(u)=\alpha_z(u)`$. Taking determinants of (16)
therefore yields the literal common-coboundary formula

```math
\alpha_z(u)=\kappa_z(u)\,
             \frac{\delta(\sigma_zu)}{\delta(u)}. \qquad(19)
```

Thus every unused-sector eigenphase is itself a central diagonal
Clifford phase, and the same delta from the actual middle band relates
the entire family to the affine-sign phases (18). This strengthens the
pointwise bounded-core statement to a common abelian phase constraint.
All assertions include arbitrary spectator states and references.

### 5.2. The remaining inference is not a Boolean degree theorem

Equations (15) and (19) describe phases in intrinsic central coordinates
u. They do not identify u with the target's computational inputs, or
construct a native simultaneous change of coordinates carrying the
complete occupied logical algebra to computational Z. Interpreting delta
as a cocycle gauge is an algebraic statement; a controlled operation on
the non-Pauli negative eigenspace need not be a full Clifford and has
not been supplied at zero band cost.

Accordingly this lemma gives no uniform degree bound for the occupied
Boolean sign or XOR target. The active-rank-two exceptional stabilizer,
nonzero probe characters, and failure of simultaneous native decoding
are not resolved by it. The next question in this subcase is whether
the common cocycle (19), together with both actual endpoint encoders,
admits a paid reduction that controls the occupied computational phase.
The two-low-degree-control decoder in Sections 1--4 alone does not
provide that reduction.

The later [determinant-defect degree theorem](DETERMINANT_DEFECT_DEGREE.md)
supplies a partial bridge: the labels whose entire defect phase is
constant form a subgroup H. One occupied-fixing CNOT correction makes
H a globally preserved computational-Z subgroup. Bounded codimension
of H, or one bounded common central quotient for the whole phase family,
therefore gives a uniform Boolean degree bound. Establishing that
codimension bound for every actual native representation remains open.
