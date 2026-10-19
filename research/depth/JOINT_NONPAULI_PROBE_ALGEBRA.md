# A joint algebra for nonzero clean-code probe characters

[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md) ·
[Occupied logical transport](OCCUPIED_LOGICAL_TRANSPORT.md) ·
[Common determinant defect](DETERMINANT_DEFECT_DEGREE.md)

The separate small algebras supplied by nonzero logical probe characters
can be placed in one common algebra when both endpoint clean images are
non-Pauli. Its symplectic rank is at most four, independently of the
number of logical generators or characters. The rank-two exceptional
stabilizer also admits a joint reduction: only one of its nonidentity
elements can occur. Its genuinely exceptional case separates the clean
transition from every silent band and gives an occupied-equivalent
zero-clean three-band circuit, hence Boolean sign degree at most five
and XOR-query degree at most four.

The generic common-algebra statement alone does not identify intrinsic
central controls with computational inputs or bound their joint phase
complexity. The degree closure uses the additional exceptional
separation proved below. Neither result changes the complete-frame
resource frontier.

## 1. Common setup and the supported scope

Use the actual three-band factorization

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,
\qquad U_C=E_CK_C,\qquad U_{AB}=E_BE_AK_A.
```

Each E is a native band on independent commuting physical Pauli axes.
Let L be the binary span of the middle B axes. For every computational
logical-Z label z, including arbitrary dirty wires, write

```math
\begin{aligned}
J_z&=E_AK_AZ(z)K_A^\dagger E_A^\dagger,&
A_z&=E_BJ_zE_B^\dagger,\\
C_z&=U_CZ(z)U_C^\dagger,&W_z&=C_zA_z.
\end{aligned}
```

The J_z and C_z are actual Hermitian Cliffords, and both logical
families are exact commuting representations. Assume they preserve
the common clean code and agree on its occupied sector. Denote the
clean images before and after B by J_c and S, and their intrinsic
active spans by L_J and L_S. For the output clean Pauli P_C, with
binary label p, put

```math
W_S=L_S+\langle p\rangle,\qquad
N=L^\perp\cap W_S^\perp.
```

First suppose the active rank of S differs from two. The retained
[probe theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#3-native-probes-of-a-compressed-logical-image)
gives, for each z, a linear character nu_z on N such that

```math
A_z R A_z=(C_z R C_z)S^{\nu_z(R)}
\quad(R\in N).
                                                               \qquad(1)
```

Assume at least one nu_z is nonzero, and assume **both J_c and S are
non-Pauli**. Then all W_z, together with their retained two-band and
Clifford factorizations, are supported in one Pauli space K satisfying

```math
\mathrm{rank}_{\mathrm{symp}}K\le4.               \qquad(2)
```

The first endpoint assumption is essential to this algebra proof. A
Pauli J_c is a real native possibility, but its degree question is
already covered by the retained
[one-Pauli-interior theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#7-one-pauli-interior-clean-image-suffices).

## 2. The clean-active lists become disjoint middle coordinates

A nonzero value in (1) makes the displacement image of the native
one-band Pauli image on its left equal to L_S. Native independence
therefore gives more than containment: L_S is spanned by a subset
I_C of the actual B-axis list. In particular L_S is contained in L,
so the middle clean permutation has no paired orbit: m=0.

In these independent B coordinates, the displacement of S is

```math
N_S=\sum_{b\in I_C}bb^*,\qquad b^*(x)=\langle b,x\rangle.
```

S fixes every B label. Let I_minus be the axes it negates literally.
Every axis in I_C belongs to I_minus. One way to see the signs is to
compare (1) with the actual one-band expansion: the active axes of
that expansion anticommute with its input Pauli; its Pauli multiplier
C_z R C_z commutes S and hence commutes all L_S axes. Equivalently,
Hermiticity fixes the characteristic vector of the independent
transvection form. Since the support of S is p+L_S, its sign on a B
axis is precisely its pairing with p.

Conjugating S by the inverse middle band consequently gives

```math
N_{J_c}=N_S+\sum_{b\in I_{\mathrm{minus}}}bb^*
       =\sum_{b\in I_A}bb^*,\qquad
I_A=I_{\mathrm{minus}}\setminus I_C.              \qquad(3)
```

The independent columns imply that L_J is the coordinate span of I_A.
Thus the two endpoint clean-active spans are disjoint coordinate
subspaces of L. The support of J_c is contained in p+L, so its clean
Pauli label p_A satisfies p_A+p in L. In particular

```math
\langle p,a\rangle=\langle p_A,a\rangle=1
\quad(a\in I_A).                                \qquad(4)
```

I_A is nonempty because J_c is non-Pauli. No assumption is made that
the endpoint active axes use the same basis as these middle
coordinates; only their common spans and displacement maps are used.

## 3. One common choice of all character representatives

Let L_A be the span of **all** first-endpoint axes, including its silent
ones, and write F_{J_z}=I+N_z. A one-band Pauli image has

```math
N_z=\sum_{a\text{ active on }K_AZ(z)K_A^\dagger}aa^*.
```

Therefore N_z is self-adjoint for the binary symplectic pairing, and
its image lies in the one fixed isotropic space L_A.

Define the linear functional lambda on L by giving value one on the
B axes in I_C and zero on the other B axes. Choose b_0 in I_C and set

```math
h_z=N_zb_0,\qquad H=\mathrm{span}\{h_z:z\}.      \qquad(5)
```

In (1), the active B-axis set is empty when nu_z(R)=0 and is exactly
I_C when nu_z(R)=1. Self-adjointness therefore gives

```math
\langle h_z,R\rangle=\nu_z(R),\qquad
N_zb-\lambda(b)h_z\in N^\perp=L+\langle p\rangle
\quad(b\in L).                                  \qquad(6)
```

This applies to every logical product, including zero-character
members. Choose any a in I_A. Both N_zb and h_z lie in L_A and hence
commute a, while every element of L also commutes a. Equation (4)
then removes the possible p coefficient in (6):

```math
N_zb-\lambda(b)h_z\in L.                         \qquad(7)
```

Pair (7) with b_0 and use self-adjointness again. Since lambda(b_0)=1,

```math
\langle h_z,b\rangle
=\langle N_zb_0,b\rangle
=\langle b_0,N_zb\rangle
=\lambda(b)\langle b_0,h_z\rangle.               \qquad(8)
```

Thus the entire cross-pairing between L and H has rank at most one.
H itself is isotropic because it is contained in L_A. Consequently
L+H has symplectic rank at most two, and

```math
K=L+H+\langle p\rangle
```

has rank at most four. The retained support theorem, using the
representatives (5), places every W_z in this same K.

This also retains the actual native factors. For a Pauli R in
K-perp, W_z and E_B commute R, so

```math
E_B(J_z R J_z)E_B^\dagger=C_z R C_z.
```

Independence says that a one-band image is Pauli only when none of its
axes is active. Hence J_z R J_z=C_z R C_z and the Clifford C_zJ_z
centralizes K-perp. The same equality says C_z R C_z commutes every
B axis. Varying R over K-perp proves C_z(L) is contained in K;
normalization of all K is proved in the following paragraph. The literal identity

```math
W_z=(C_zE_BC_z^\dagger)(C_zJ_z)E_B^\dagger       \qquad(9)
```

therefore has both native bands and its middle Clifford in K. A
fixed Clifford choice of coordinates describes at most two quantum
coordinates, possibly many commuting controls, and identity spectators.
Equation (9) neither compiles that coordinate change for free nor
inserts an uncharged decoder.

The common space is also normalized by the complete logical Clifford
families. Since all N_z have isotropic image in L_A, N_zN_v=0 and
F_{J_z} fixes every h_v. Equation (7) maps L into L+H. Moreover J_z
commutes J_c, so it preserves J_c's support coset p_A+L_J and maps p_A
into that coset. As p_A+p lies in L, these statements prove

```math
J_z(K)=K.
```

The Clifford C_zJ_z in (9) centralizes K-perp and therefore normalizes
K. It follows that C_z(K)=K as well. This is one common invariant
algebra, rather than separate coordinates for individual generators.

## 4. The rank-two exception has only one joint nontrivial type

For intrinsic clean rank two, use the retained exact coordinates

```math
S=\mathrm{diag}_v(Y_t,X_t),\qquad
W_0=\mathrm{diag}_v(Y_t,I_t),\qquad
W_1=\mathrm{diag}_v(I_t,X_t).
```

The pointwise Clifford stabilizer is {I,W_0,W_1,S}. Its three nonidentity
displacement planes are respectively

```math
D_{W_0}=\langle Z_v,Y_t\rangle,\quad
D_{W_1}=\langle Z_v,X_t\rangle,\quad
D_S=\langle Z_v,Z_t\rangle.                     \qquad(10)
```

Any two have a nonisotropic sum. If a probe comparison has nonidentity
value Q, the native displacement argument puts D_Q inside the one
isotropic middle span L. Therefore **only one nonidentity Q can occur
across the entire logical family and every probe**. The homomorphisms
in the retained exceptional statement consequently have image at most
one common copy of F_2, not unrelated Klein-four images.

If Q=S, Sections 2--3 apply unchanged when J_c is also non-Pauli.
Conjugating a discrepancy by C_z when forming W_z can interchange
W_0 and W_1; the assertion about one type concerns the displayed
comparison A_z R A_z versus C_z R C_z, and its kernels are unchanged.

Suppose Q is W_0 or W_1, and put P_t=Y_t or X_t respectively. The
active B subset has exactly two members and spans
D_Q=span{Z_v,P_t}. S preserves the B list and swaps P_t with Z_vP_t.
Independence excludes a list containing Z_v and only one of those
two labels: its image would require a third dependent axis. Thus the
actual middle pair is exactly

```math
\{P_t,Z_vP_t\},                                  \qquad(11)
```

up to literal signs. In particular it has no spectator dressing.
Every other middle axis is a singleton with no t-core component; it
commutes the entire space W_S=span{Z_v,X_t,Z_t}. There is exactly one
pair, with difference Z_v.

The pair (11) rotates just one v-block of S by a one-qubit Clifford
and fixes the other. Depending on orientations it either leaves S
unchanged or replaces one of its distinct Pauli blocks by a signed
Z_t. Thus J_c still has four flat Pauli terms, intrinsic rank two,
and W_J=W_S. Both endpoint active spans and both clean Pauli labels
are contained in this same three-dimensional space. Every endpoint
silent axis commutes the clean Pauli and its active span, hence
commutes W_S.

Consequently all three silent bands commute the entire clean-relevant
two-wire factor. This is a strict local separation for the genuinely
exceptional probe type, stronger than a bounded endpoint-rank claim.
It supplies no shared clean flip for the original representation.

### 4.1. Paid occupied replacement and the degree-five closure

Write a and s for each endpoint's active and silent factors, and write
E_B,pair and E_B,s for the middle pair and singleton factors. The
literal clean-relevant operator is

```math
F_0=E_{C,a}^\dagger E_{B,\mathrm{pair}}E_{A,a},
\qquad F_0P_AF_0^\dagger=P_C.
```

Both signed clean Paulis are nonradical elements of W_S. There is an
actual Clifford C_0 in its algebra with $`C_0P_AC_0^\dagger=P_C`$: a core
Clifford rotates the nonidentity core Pauli, a controlled core Pauli
with control v adjusts its Z_v factor if necessary, and a core Pauli
corrects its literal sign. Every one of these operations is supported
in W_S. Thus

```math
G_0=C_0^\dagger F_0\in
\mathrm{Alg}(W_S\cap p_A^\perp)
=\mathrm{Alg}(\langle p_A,Z_v\rangle).            \qquad(12)
```

All three silent bands commute C_0 and G_0 and preserve P_A. On its
occupied +1 sector, G_0 is a scalar unitary q(Z_v). An actual Clifford
simultaneously diagonalizing P_A and Z_v puts its two scalar values
in the native coefficient ring. The retained
[unit-circle lemma](../../docs/FAULT_TOLERANT_COMPILER.md) makes both
values eighth roots. Therefore, with its scalar retained,

```math
q(Z_v)=\gamma T_{Z_v}^{r},\qquad
\gamma\in\{\omega^j:j\in\mathbb Z\},\quad r\in\mathbb Z/8.
                                                               \qquad(13)
```

Write $`r=2s+\epsilon`$ with $`\epsilon\in\{0,1\}`$, and absorb the
commuting Clifford $`S_{Z_v}^s`$ into C_0, obtaining C_1. The new actual circuit

```math
\widetilde V=
\gamma K_C^\dagger C_1 E_{C,s}^\dagger
        (E_{B,s}T_{Z_v}^{\epsilon})E_{A,s}K_A     \qquad(14)
```

has exactly the original occupied action, with complete clean return.
The replacement of G_0 by (13) is made only on that initialized sector.
It is not a claim that (14) has the original unused-sector action.
All preceding factorizations, through (12), were full-space identities.

Equation (14) still has at most three actual native bands. The only
new candidate axis, Z_v, is independent of the entire middle singleton
list: otherwise the difference of the removed independent pair would
be a combination of its other members. Each surviving axis lies in
W_S-perp. Choose fixed Clifford coordinates sending P_A and a
symplectic partner in W_S to the clean Z and X. Every surviving axis,
including Z_v, then acts as identity on that clean coordinate. Deleting
it preserves the independence of every list. The two endpoint
Clifford isometries restrict to actual logical Cliffords.

This constructs a **zero-clean logical circuit of at most three native
bands** on all the remaining wires, including every arbitrary dirty
wire. If the original occupied action is a Boolean sign G, the
[zero-clean three-layer theorem](THREE_LAYER_NATIVE_STRUCTURE.md)
therefore gives

```math
\deg G\le5,\qquad \deg f\le4
\quad\text{for a Boolean XOR query }f.           \qquad(15)
```

For the XOR statement, first conjugate the query output by Hadamards
to form its Boolean sign representation, and require the non-S
exceptional probe witness for that representation. Its sign has
degree one plus the degree of a nonzero f. The logical-Z probe family
need not be unchanged by this output basis change.

The construction allocates no new work. It removes the six
clean-relevant T gates and inserts at most one, so its T count is at
most the displayed original count minus five. Its Clifford operations
are actual syntheses: a fixed number of width-Q Clifford words costs
$`O(Q^2)`$ elementary gates, with their actual inverse words charged
separately. The fixed scalar gamma is also realizable without T gates;
the native Clifford word $`(HS)^3`$ realizes omega times identity. The
logical circuit can be embedded with an idle returned clean wire if
the original interface retains that wire. All initialized operator
identities include arbitrary reference entanglement and dirty return.

This closure is specific to a non-S exceptional probe value. It does
not close the generic S-valued branch, all rank-two clean codes, or
arbitrary paired middle bands.

The [exact two-wire fixture](../../tests/test_joint_nonpauli_probe_algebra.py)
checks the separated native pair, the literal occupied scalar, every
original physical column, actual inverses, and the changed unused
extension. Its occupied target is a T-dagger phase, not a Boolean
sign. It tests the replacement mechanism after separation, not the
existence of a non-S logical probe witness or a Boolean lower bound.

## 5. What remains unpaid

For the generic nonzero type, the first endpoint may be Pauli. Then
I_A is empty and the elimination of p in (6) is unavailable. The
common-rank proof above does not cover that branch; its degree bound
is supplied independently by the retained one-Pauli-interior theorem.

Even in the supported two-core case, the number of central controls
can grow. A bounded intrinsic quantum core alone does not yield a
bounded computational-Z defect codimension or a Boolean degree bound.
That requires an additional common transport or phase argument, such
as the separately scoped determinant result linked above. The general
one-clean three-layer question and the complete-frame frontier remain
unchanged.
