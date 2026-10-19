# A rank-one interface for an unbounded clean transition

[Clean-transition parameters](CLEAN_TRANSITION_PARAMETERS.md) ·
[Bounded-relation matching](BOUNDED_RELATION_NATIVE_MATCHING.md#8-one-clean-qubit-no-middle-pairs-and-bounded-active-ranks) ·
[Finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md)

With no middle pairs, the clean transition can retain its original
three active lists instead of flattening them into one dependent list.
This gives four phase factors, each with commuting axes, on the occupied
logical space, with relation bounds $`(0,1,1,0)`$. The two adjacent output factors
have cross-pairing rank at most one. Their central support can still
grow without bound. Even an occupied identity circuit shows that this
one-core interface does not bound the rank-profile distortion used in
the finite-insertion proof.

## 1. Literal clean alignment and evaluation

Use the actual one-clean notation

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,\qquad
P_A=K_AZ_cK_A^\dagger,\quad P_C=K_CZ_cK_C^\dagger.
```

Assume both interior clean images are non-Pauli and their native clean
comparison has no paired middle orbit. Split A and C into clean-active
and clean-silent lists. Split B into its negative and positive
singletons. The retained no-pair reduction gives the isotropic span

```math
L_0=L_{A,\mathrm{active}}+L_{B,-}
   =L_{A,\mathrm{active}}+L_{C,\mathrm{active}},
```

and the literal commuting phase

```math
F_0=E_{C,\mathrm{active}}^\dagger E_{B,-}E_{A,\mathrm{active}},
\qquad F_0P_AF_0^\dagger=P_C.                    \qquad(1)
```

Every B axis commutes L_0. Every active axis in (1) anticommutes
with P_A; for C this also follows from $`p_A+p_C\in L_0`$.
Extract an actual commuting-axis Clifford $`C_0\in\mathrm{Alg}(L_0)`$
with the same signed action on P_A. Then

```math
G_0=C_0^\dagger F_0\in\mathrm{Alg}(L_0\cap P_A^\perp).
```

Choose q to be an actual signed A-active axis. On the commuting Pauli
algebra define $`\tau`$ by literal evaluation at the positive q
eigenvalue:

```math
\tau(Q)=I,\qquad
\tau(R)=R Q^{\langle r,p_A\rangle}\quad(r\in L_0).
                                                               \qquad(2)
```

Products in (2) are literal products of commuting signed Hermitian
Paulis. Equivalently, use a signed commuting basis consisting of Q
and a basis of $`L_0\cap P_A^\perp`$. This makes tau a unital
star-algebra homomorphism. It fixes G_0 exactly, on both clean sectors.
A diagonal Clifford restricted to one coordinate eigenvalue is again
a diagonal Clifford, including its fixed scalar. Thus

```math
D_0=\tau(C_0),\qquad
G_0=D_0^\dagger
 (E_{C,\mathrm{active}}^\tau)^\dagger
 E_{B,-}^\tau E_{A,\mathrm{active}}^\tau           \qquad(3)
```

is a literal full-space identity. In particular this step does not
replace a chosen unused-sector extension. Each projected T factor
retains its orientation and its signed axis. Any scalar from an
identity axis or Clifford restriction is retained.

Put

```math
L=D_0C_0^\dagger,\qquad K_C''=LK_C,\qquad
E_{C,s}''=L E_{C,\mathrm{silent}}L^\dagger,
```

and define

```math
\widehat E_A=E_{A,\mathrm{active}}^\tau E_{A,\mathrm{silent}},
\quad
\widehat E_B=E_{B,-}^\tau E_{B,+},
\quad
\widehat E_{C,a}=E_{C,\mathrm{active}}^\tau.
```

Commuting only factors justified by the native no-pair structure gives

```math
V=(K_C'')^\dagger (E_{C,s}'')^\dagger
       \widehat E_{C,a}^\dagger\widehat E_B\widehat E_AK_A.
                                                               \qquad(4)
```

All four factors preserve P_A, and both endpoints map Z_c to P_A.
The active C factor commutes every individual axis of the merged B
factor. Equation (4) is still full-space equality.

The clean derivative also constrains the retained scalar phase.
In signed q/central coordinates, every original active axis has the
form $`\pm QZ^u`$. If f_0(x) is the literal exponent of F_0 at
$`q=+1`$, its exponent on both q sectors is

```math
f(t,x)=f_0(x)+t[N-2f_0(x)]\pmod8,
```

where N is the signed sum of the original T orientations, with the
C-active orientations reversed. The fixed constants from negative
axes are included in f_0. Since the derivative in (1) is Pauli,
$`N-2f_0(x)=c+4\ell(x)\pmod8`$ for a constant c and a linear
central parity ell. Therefore
$`f_0(x)=c_0+2\ell_0(x)\pmod4`$. Its weighted-cubic coefficients
leave a Boolean polynomial of degree at most three after subtracting
an actual diagonal-Clifford exponent and the fixed scalar. The same
description holds for G_0. This is a statement in the common commuting
algebra's coordinates, not a degree bound in computational input
coordinates. The disjoint-CCZ family in
[the parameter obstruction](CLEAN_TRANSITION_PARAMETERS.md) shows that
the residual cubic support can still be arbitrarily large.

## 2. Projected lists and the single coupled core

The kernel of tau on L_0 is exactly the binary line spanned by q.
The merged A list stays independent after its Q factor becomes I and
is dropped: q was itself one member of the original independent A
list. The merged B list and active C list each gain at most one
relation. These three lists commute Q, so projection onto the positive
P_A code is injective on their spans: P_A anticommutes Q.

The conjugated C-silent list remains independent and also projects
injectively. Indeed L fixes each original C-active axis R, its silent
axes commute R, and P_A anticommutes R. Hence P_A cannot belong to
the C-silent span. After fixed Clifford coordinates and clean
restriction, the four logical relation bounds, in input-to-output
order, are therefore

```math
(r_A,r_B,r_{C,a},r_{C,s})\le(0,1,1,0).             \qquad(5)
```

For an original active C label r and silent C label s, L fixes both
r and q. The surviving cross-pairing is

```math
\langle r+q,Ls\rangle=\langle q,s\rangle.         \qquad(6)
```

Its rows are identical, so its rank is at most one. The union of the
two isotropic C-side spans has symplectic rank at most two: at most
one quantum core, with an unrestricted central span. Clean projection
preserves this pairing.

### 2.1. The entire output union has at most one relation

The two separate bounds in (5) conceal a stronger fact. Let L_C be the
span of the original full independent C list, and define the linear map

```math
M(v)=Lv+\langle Lv,p_A\rangle q,\qquad v\in L_C.
```

Its images are exactly the retained active and silent C labels. On the
ambient Pauli space, the map after L is a projection with kernel
$`\langle q\rangle`$. Thus its restriction to L_C has kernel dimension
at most one.

Clean restriction introduces no further relation in this entire mixed
union. If $`p_A=M(v)`$, then Lv is either p_A or $`p_A+q`$. Both
pair to one with any original active C axis r. But L fixes r, and
Lv and r belong to the isotropic space L(L_C), a contradiction. The
union therefore has relation dimension at most one before and after
clean restriction.

If the cross-pairing (6) is nonzero, the union is fully independent.
Indeed q cannot belong to L(L_C): some silent axis in that isotropic
space pairs to one with q. In this case M has no kernel on L_C.

### 2.2. Exact intrinsic coordinates for the coupled output factors

Suppose the cross-pairing is nonzero. Select an active label r_0 and
a silent label s_0 that pair to one. Every active label pairs to one
with s_0. Split the silent list into labels that pair to one with r_0
and labels that commute with it. The following are independent central
labels in the combined output algebra:

```math
r+r_0\quad(r\ne r_0\text{ active}),\qquad
s+s_0\quad(s\ne s_0\text{ coupled silent}),\qquad
s\quad(s\text{ uncoupled silent}).
```

They commute each other and both r_0,s_0. The preceding union
independence proves their independence together with this core pair.
A fixed paid Clifford can therefore put the two lists in the form

```math
\begin{array}{ll}
\text{active:}&Z_t,\ Z_tZ_1,\ldots,Z_tZ_{a-1},\\
\text{coupled silent:}&X_t,\ X_tZ_a,\ldots,X_tZ_{a+b-2},\\
\text{uncoupled silent:}&Z_{a+b-1},\ldots,Z_{a+b+h-2}.
\end{array}
```

Here a and b are positive, h may be zero, and empty index ranges are
omitted. Independence permits literal Pauli-sign correction by a
fixed Pauli; every original T orientation and common scalar is kept.
The two controlled core phases use disjoint independent intrinsic
central coordinates. These coordinates have not been identified with
computational input bits. Their existence alone supplies no Boolean
degree bound in the required input coordinates.

### 2.3. A full odd-coset clean-flip criterion and its dual

Write $`\chi(v)=\langle v,p_A\rangle`$ on L_0. Suppose there is a
label q in L_0 such that

```math
\chi(q)=1,\qquad
\langle q,L_{A,\mathrm{silent}}\rangle=0,
\qquad\langle q,L_{C,\mathrm{silent}}\rangle=0.       \qquad(11)
```

The associated Hermitian Pauli Q commutes every actual axis of A, B
and C, and anticommutes with both P_A and P_C. Thus both endpoint
preimages $`K_A^\dagger QK_A`$ and $`K_C^\dagger QK_C`$ are
Paulis flipping the physical clean Z_c. The complete circuit
transports one to the other literally.

An actual clean-controlled Pauli dressing at each endpoint makes
both these preimages X_c while fixing every clean-zero vector. Exact
clean return already makes the full circuit block diagonal in Z_c.
After dressing it also commutes X_c, so its full operator is
$`I_c\otimes G`$, with G the original occupied logical action and
its literal phase. These are endpoint Clifford changes only: all
three original physical native lists remain independent and unchanged.
The old clean wire can now be an arbitrary returned dirty wire. No
individual native list has been projected or assigned a free relation.

Consequently the retained zero-clean three-layer theorem gives sign
degree at most five, and Boolean XOR predicate degree at most four,
whenever (11) holds. The forward and actual inverse keep their T counts
and native T-depth, with $`O(Q)`$ additional elementary Clifford gates
on Q physical wires and no new wire. The inverse dressings are paid
as their actual inverse words. For XOR, apply the argument to the
output-Hadamard sign representation; its clean data and criterion
(11) are unchanged.

In particular, it is enough to find any clean-odd label in the whole
A-active span that annihilates the C-silent span. The label need not
be an actual member of the native A list. If the rank-one interface
vanishes for the q chosen in Section 1, this condition holds.

The failure of (11) has an exact linear dual description. Let
$`L_{C,s}'=C_0^\dagger L_{C,\mathrm{silent}}C_0`$, which preserves
P_A and has the same pairing rows on L_0 as the original silent span.
Then there exist $`s_A\in L_{A,\mathrm{silent}}`$ and
$`s_C\in L_{C,s}'`$ such that

```math
\chi(v)=\langle v,s_A+s_C\rangle\quad(v\in L_0).       \qquad(12)
```

Indeed a functional vanishes on the common kernel of a family of
linear functionals exactly when it belongs to their row span. Put
$`\alpha(v)=\langle v,s_A\rangle`$ and
$`\beta(v)=\langle v,s_C\rangle`$. Every active A axis has
$`(\alpha,\beta)=(0,1)`$, every active C axis has (1,0), and every
negative B axis lies in exactly one of these two affine classes.
This partitions the original independent negative B list. It does
not authorize deleting either resulting phase or supplying a
sector-dependent native operation for free. Only the no-common-flip
case inside L_0 is diagnosed; no absence claim about every physical
Pauli outside L_0 is needed.

The general four-factor reduction (4) is an exact structural reduction,
not a three-native-band compilation. The two lists with possible relations are algebraic phase
factors; dependent axes cannot be assigned a free physical band.
The coordinate Cliffords, C_0, D_0, their literal scalars and their
actual inverse words must be paid. Fixed Clifford synthesis needs
polynomial Clifford cost in the physical width; no helper is required.
Deleting the clean wire produces the full occupied logical operator
after these paid coordinate changes, but (5) alone does not bound its
Boolean degree or improve the complete-frame frontier.

## 3. One core does not bound native rank-profile distortion

Here is a clean-derived identity witness with arbitrarily large
distortion. Use physical wires c,t and $`z_1,\ldots,z_m`$, with
$`m\ge1`$, and set

```math
\begin{gathered}
Q=Z_c,\qquad R_i=Z_cZ_tZ_i,\qquad K_A=K_C=H_c,\\
E_A=T_QT_{X_t},\qquad
E_B=T_Q^\dagger\prod_iT_{R_i},\\
E_C=\left(\prod_iT_{R_i}\right)T_{X_cX_t}.
\end{gathered}                                      \qquad(13)
```

Every list is physically independent. The last C axis commutes every
R_i by two anticommutations. Here $`P_A=P_C=X_c`$, the endpoint
active ranks are $`k_A=1`$ and $`k_C=m`$, all B axes are negative
singletons, the middle pair count is $`m_{\mathrm{pair}}=0`$, and
$`F_0=I`$ literally. Both clean interiors are non-Pauli. The full
actual operator is

```math
V=H_cT_{X_cX_t}^\dagger T_{X_t}H_c
 =|0\rangle\langle0|_c\otimes I
  +|1\rangle\langle1|_c\otimes\omega^{-1}S_{X_t}.
                                                               \qquad(14)
```

All z wires and any external reference are untouched. Thus its
occupied Boolean target is identity, including literal phase.
There is no target-hardness obstruction.

Applying (2)--(4), then restricting to the clean code, gives

```math
A=C=T_{X_t},\qquad B=\prod_iT_{Z_tZ_i},\qquad
F=B^\dagger,\qquad C^\dagger FBA=I.               \qquad(15)
```

Here the logical C is the projected C-silent factor; F is the inverse
projected C-active factor. The original C clean rank remains m, and
must not be confused with the single axis of the logical silent C.

The insertion F commutes every B axis. Its axes pair with C's single
axis with rank one, and their total Pauli span has symplectic rank
two. Nevertheless take the logical Z_t probe. Its actual Hermitian
Clifford image under A is

```math
J=T_{X_t}Z_tT_{X_t}^\dagger=(Z_t-Y_t)/\sqrt2,
\qquad \Phi(Z_t)=Y_t
```

at the binary-label level. For B's native displacement-rank profile,

```math
w(Z_t)=0,\qquad w(\Phi(Z_t))=w(Y_t)=m.            \qquad(16)
```

The latter equality follows because Y_t anticommutes with all m
independent B axes. The corresponding Hermitian Clifford has Pauli
coefficient norm $`2^{m/2}`$. Thus the finite-insertion proof's
rank-profile estimate cannot replace support dimension by this
adjacent pairing rank or quantum-core count, even on an actual
clean-derived occupied identity. The cancellation $`FB=I`$ shows
precisely why this is a representation limitation.

For a concrete native ledger, let
$`N=\prod_i\mathrm{CX}_{c,i}\mathrm{CX}_{t,i}`$ and
$`U=N\mathrm{CX}_{c,t}H_c`$. These map Z_i to R_i, and U
maps Z_c to $`X_cX_t`$. Use H_t coordinates for A, N coordinates
for B, and U coordinates for C. The one clean wire is the c already
specified; the t and z wires remain arbitrary.

| Actual schedule | T count | Native T-depth | CNOT | H | Clean work | Added helper |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Displayed forward word | $`2m+4`$ | $`\le3`$ | $`8m+2`$ | 6 | 1 | 0 |
| Its actual reversed inverse | $`2m+4`$ | $`\le3`$ | $`8m+2`$ | 6 | 1 | 0 |

The [exact native fixture](../../tests/test_rank_one_clean_phase_interface.py)
checks all physical columns, the nontrivial unused block, actual
inverses, independence and the profile jump for small m. The formulas
above establish the arbitrary-width claims.

The next question is whether matching the correlated product FB can
use the native union and dual constraints without separately bounding
F's distortion. The present finite-insertion norm estimate does not
supply that step. The
[computational phase analysis](COMPUTATIONAL_PHASE_DERIVATIVES.md)
supplies an actual input-derivative bridge in its zero-character
subcase and identifies the remaining quotient and native cancellation
gap. Neither intrinsic low-degree controls nor a formal removal of
the resulting phase supplies that missing implementation.
