# Occupied logical transport and a diagonal pair-center normal form

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md) ·
[Three-layer native structure](THREE_LAYER_NATIVE_STRUCTURE.md)

A clean-swapped middle-band pair supplies a common **clean-preserving**
Pauli transport even when no common clean flip exists. This uses the
actual three native bands and gives an identity on the entire physical
space. For an occupied Boolean sign, these transports identify directions
along which its Boolean derivative is affine. Removing the resulting
quadratic part gives a paid normal form in which every pair difference
is a globally preserved computational-Z parity.

The result is an unrestricted structural reduction for exact one-clean
three-layer signs and XOR queries. It does not assert a uniform Boolean
degree bound. The surviving normal form still permits non-Pauli clean
codes and unused-sector freedom, even when there are no pairs; projecting
a nonempty pair center can also create arbitrarily many axis relations.
No complete-frame resource bound changes.

## 1. Literal model and common Pauli transports

Let Q be the total physical width, including one clean qubit c and every
arbitrary dirty helper. Suppose an actual native circuit has T-depth at
most three and returns c exactly. Absorb its endpoint and intermediate
Cliffords into the literal factorization

```math
V=K_C^\dagger E_C^\dagger E_B E_AK_A.
```

Each E is a native T/T-dagger band on independent commuting signed
Hermitian Pauli axes. Empty bands are permitted. Every scalar and inverse
in this factorization belongs to the actual operator.

Exact clean return on every other input implies $`[V,Z_c]=0`$ by
unitarity. Define

```math
\begin{aligned}
J&=E_AK_AZ_cK_A^\dagger E_A^\dagger,\\
S&=E_CK_CZ_cK_C^\dagger E_C^\dagger
  =E_BJE_B^\dagger.
\end{aligned}
```

Both J and S are Hermitian Cliffords. The
[native one-band matching theorem](ONE_BAND_AXIS_MATCHING.md) says that
J permutes the individual middle-band axes, up to signs. Since E_B
commutes those axes, S induces exactly the same signed permutation.
Their binary symplectic actions have the one-band form

```math
F_J=I+N_A,\qquad F_S=I+N_C,
```

where the displacement images are precisely the spans of the respective
clean-active endpoint axes. These displacement maps square to zero.
The independent middle-axis list therefore has only singleton and paired
orbits, as in the retained three-layer argument.

Let $`\{B_j,D_j\}`$, $`j=1,\ldots,m`$, be its paired orbits. Use the
literal signed Hermitian products and their binary labels

```math
R_j=B_jD_j,\qquad r_j=b_j+d_j.
```

The common signed permutation gives

```math
r_j=N_A b_j=N_C b_j.                              \qquad(1)
```

Thus r_j lies in both clean-active endpoint spans and in the middle-band
span. Consequently **R_j commutes every axis of all three entire bands**.
The r_j are independent because each uses its own two members of the
independent middle list. A signed swap squaring to identity has equal
signs on its two members; therefore every R_j also commutes J and S.

Define signed physical Paulis

```math
P_{A,j}=K_A^\dagger R_jK_A,\qquad
P_{C,j}=K_C^\dagger R_jK_C.
```

Commuting R_j through the three actual bands proves the full operator
identity

```math
V P_{A,j}V^\dagger=P_{C,j}.                       \qquad(2)
```

Both endpoint Paulis commute Z_c. Equation (2) holds on the unused clean
sector as well as the occupied sector; no extension of an occupied
Clifford equality has been inferred. These are clean-preserving
transports, not the shared clean flips excluded by the retained examples.

Their restrictions to c=0 remain an independent commuting logical Pauli
family. Indeed, a product restricting to identity would be a scalar or a
signed Z_c before restriction. Independence excludes a nonempty scalar
product. A product equal to signed Z_c would make the corresponding
product of R_j equal to J, or S, up to sign. That product is Pauli,
whereas any paired orbit makes the clean displacement nonzero and hence
J and S non-Pauli. If m=0, the family and subsequent normalization are
empty and the assertions remain valid.

## 2. Affine Boolean derivatives and a quotient of the active variables

Suppose the occupied action is the literal Boolean sign

```math
G|x,w\rangle=(-1)^{g(x)}|x,w\rangle,
```

where x consists of n active sign variables and w comprises every
arbitrary dirty variable. In particular, g is independent of w.
A freely retained common scalar would cancel in all conjugation identities
and does not alter the following conclusion.

Write $`p_{A,j},p_{C,j}`$ for the logical restrictions of the endpoint
Paulis. A diagonal sign preserves a Pauli's X mask. Comparing the literal
basis columns in

```math
G p_{A,j}G^\dagger=p_{C,j}
```

therefore shows that their X masks coincide and that

```math
\Delta_{a_j}g(x):=g(x)+g(x+a_j)
             =\ell_j(x)+\epsilon_j,              \qquad(3)
```

with $`a_j\in\mathbb F_2^n`$ the projection of that X mask onto the
active variables, $`\ell_j`$ linear, and $`\epsilon_j`$ constant.
The right side cannot contain a dirty variable because the left side
contains none. Literal Pauli phases determine the constant; they cannot
supply an unrecorded complex phase to a Boolean derivative.

Let

```math
U=\mathrm{span}_{\mathbb F_2}\{a_1,\ldots,a_m\},
\qquad r=\dim U.
```

Every derivative in U is affine: use
$`\Delta_{a+b}g(x)=\Delta_a g(x)+\Delta_b g(x+a)`$.
Choose one invertible linear change x=My sending the first r coordinate
directions onto U. In these common coordinates, the algebraic normal
form of $`g(My)`$ has no degree-three-or-higher monomial involving one
of the first r variables. Such a monomial would leave a distinct
nonaffine monomial after differentiation in that coordinate, contradicting
(3). Hence

```math
\begin{aligned}
g(My)&=q(y)+h(y_{r+1},\ldots,y_n),\\
\deg q&\le2,\qquad
\deg g\le\max\{2,n-r\}.                          \qquad(4)
\end{aligned}
```

Choose q to contain precisely the monomials touching the first r
coordinates; the constant and the remaining monomials belong to h.
This fixes the literal sign without introducing a scalar correction.

Equation (4) is a circuit-dependent degree certificate. It does not say
that the projected rank r must grow with the number m of paired orbits:
pair transports may be Z-type on the active variables or have X support
only on already dirty variables.

For a literal XOR query, Hadamards on its output turn its occupied action
into $`g(x,z)=zf(x)`$, with n counting the original address and output
variables, excluding all dirty helpers. Since $`\deg g=\deg f+1`$ for
nonzero f, the same circuit certificate gives

```math
\deg f\le\max\{1,n-r-1\}.                       \qquad(5)
```

The coordinate change M may mix the original output with addresses;
there is no claim that the quotient h is itself a query with that same
output designation. The sign degree and its implication (5) do not
require such a claim.

For a full-degree Boolean function on a address bits, with a at least
three, the sign $`g(x,z)=zf(x)`$ has its unique degree-$`a+1`$
monomial. Differentiating it in any nonzero direction leaves a nonzero
degree-a sum; derivatives of lower-degree terms cannot cancel that sum.
Thus U is zero. Every pair transport then has zero X mask on all address
and output wires, although it may have dirty X support. In particular,
for AND queries the normalization below needs no active-variable change
or quadratic subtraction and preserves the original query sign exactly.

## 3. A paid common normalization of all pair differences

The preceding quotient can be implemented with endpoint Cliffords.
The linear change is an actual CNOT circuit. The sign $`(-1)^{q(y)}`$
is an actual product of Z and CZ gates, with every CZ expanded into
H-CNOT-H when elementary counts are required. Applying these operations
to V gives a circuit whose occupied action is h alone. The first r
coordinates are now additional returned arbitrary dirty variables.
Total physical width and the single clean allocation are unchanged.

Keep the transformed pair transports and denote them again by
$`P_{A,j},P_{C,j}`$. Their X masks now touch only dirty coordinates;
the remaining active components are Z-type. Thus h commutes their
occupied input restrictions. Equation (2) implies

```math
P_{C,j}=P_{A,j}Z_c^{e_j},\qquad e_j\in\mathbb F_2, \qquad(6)
```

on the full physical space. Both are Paulis commuting Z_c, and their
occupied restrictions coincide literally, so no minus sign is permitted
in (6).

The logical input family is independent and commuting. The values e_j
therefore extend linearly to a character on its binary span. By
nondegeneracy of the logical symplectic pairing, choose a Hermitian
logical Pauli Q_0 such that

```math
\langle Q_0,p_{A,j}\rangle=e_j\quad\text{for every }j.
```

The actual clean-controlled Pauli

```math
D=|0\rangle\langle0|_c\otimes I
  +|1\rangle\langle1|_c\otimes Q_0
```

is Clifford and fixes every occupied vector literally. Its conjugation
sends every P_C,j in (6) to P_A,j simultaneously. Left multiplication by
D consequently makes the complete circuit commute with all P_A,j on
**both** clean sectors, while leaving the occupied action h unchanged.

The active parts of these commuting Paulis are Z-type, so their dirty
parts commute with one another. A Clifford acting only on the dirty
register simultaneously diagonalizes that abelian dirty Pauli family.
Conjugating the circuit by it leaves h unchanged. Every paired difference
now has one identical computational-Z-parity preimage at both endpoints,
and the complete physical circuit preserves each of those parities.

Only endpoint Cliffords were changed. The actual three native bands,
their clean images J,S, their active ranks, and their middle paired orbits
are retained. In particular, no conditional coordinate circuit depending
on a selected parity sector has been supplied for free.

### 3.1. Native ledger and operator scope

CNOT Gaussian elimination, a quadratic diagonal Clifford, and dirty
Clifford diagonalization each require $`O(Q^2)`$ elementary H/S/CNOT
gates with absolute constants. A clean-controlled Pauli requires
$`O(Q)`$ such gates and no helper: controlled X is CNOT, controlled Z
is H-CNOT-H, and controlled Y is an S-conjugated CNOT. A sign on Q_0 is
an additional Z on the control. Every S inverse can be emitted as three
S gates.

Writing $`\widehat V`$ for the normalized actual circuit, its simultaneous
ledgers are

```math
\begin{aligned}
T(\widehat V)&=T(V),&
D_T(\widehat V)&\le D_T(V),&
G(\widehat V)&\le G(V)+O(Q^2),\\
T(\widehat V^\dagger)&=T(V^\dagger),&
D_T(\widehat V^\dagger)&\le D_T(V^\dagger),&
G(\widehat V^\dagger)&\le G(V^\dagger)+O(Q^2).
\end{aligned}
```

The inverse is the reversed actual native word, with all inverse
Cliffords paid separately. These bounds do not identify a forward word
with its inverse or compare G(V) with an uncharged inverse count.
There are no additional wires, measurements, supplied states, errors,
or discarded branches. The normalized occupied action is h tensor the
identity on all old and newly inactive dirty wires. Exact clean return
is preserved. Full operator identities (1), (2), and (6), and the occupied
identity for h, extend to arbitrary external references.

## 4. The exact remaining diagonal pair-center question

In this normal form, the common preimage of every R_j is a preserved
computational-Z parity. Consequently the circuit can be studied on each
of their joint sectors without imposing a new initial-state promise on
the actual circuit. The middle paired phases become known Cliffords or
scalars in a fixed sector. With the native convention
$`T_P=(I+P)/2+\omega(I-P)/2`$, $`\omega=e^{i\pi/4}`$, their literal
actions are

| Pair factor | $`R_j=+1`$ | $`R_j=-1`$ |
| --- | --- | --- |
| $`T_{B_j}T_{D_j}`$ | $`S_{B_j}`$ | $`\omega I`$ |
| $`T_{B_j}T_{D_j}^\dagger`$ | $`I`$ | $`\omega^{-1}S_{B_j}`$ |

Other orientations are obtained by retaining the appropriate fixed
Clifford factors and scalars. The sector dependence of this entire
Clifford word remains part of the physical operator. In particular,
choosing a different simplifying Clifford in every sector is not an
available free native operation.

Quotienting independent endpoint axes by m pair parities can create up
to m relations. A subsequent genuine Pauli clean-code quotient can
increase this number by one. The intermediate clean code may also remain
non-Pauli. Thus the retained
[one-relation axis matching theorem](ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation)
does not automatically control the resulting endpoint factors.

The first unsupported inference would be to remove the now sectorwise
Clifford middle pairs, apply a one-relation or full-space two-layer
result to the remaining projected endpoints, and infer one common
Boolean-degree bound. Neither its relation hypothesis nor a uniform
control of the chosen sectorwise coordinates has been established.

This is not the only surviving issue. When m=0 the pair center is empty
and the normalization may do nothing. The clean-silent endpoint bands
and a potentially non-Pauli clean code still remain. An occupied logical
image agreeing with a Clifford does not become Clifford on the unused
sector; the complete-family examples retain that failure even with no
middle paired orbit. Thus pair-free does not mean proved, and the
many-relation issue must not be presented as the sole remaining gap.

The new question is therefore precise: **do the actual independent
endpoint bands, the complete occupied logical family, and their common,
globally preserved diagonal pair center bound the completed Boolean
phase degree uniformly, including its joint dependence on all pair-parity
values?** A proof must retain the conditional Clifford word, all projected
relations, clean-silent bands, and the unused-sector freedom of a
non-Pauli clean code. Conversely,
a native family in this normal form with unbounded completed degree
would settle the original exact three-layer test constructively.

The quotient removes every higher-degree dependence detected by the
non-diagonal pair transports and shows that the remaining difficulty can
be placed entirely in this diagonal pair-center form. It does not prove
that m or the number of surviving active variables is bounded, and it
does not change R96, selected R79, or the scoped R97 obstruction.

### 4.1. A pair-free closure with an explicit silent-band hypothesis

A useful subcase can be closed without extending an occupied Clifford
identity to the unused sector. Suppose m=0, and suppose every
clean-silent A axis, pulled back through K_A, is a computational-Z Pauli;
make the analogous assumption for every clean-silent C axis through K_C.
These are additional hypotheses on the actual endpoint factors.
Under them, **every occupied diagonal action of V is diagonal C3**.
In particular, an occupied Boolean sign has degree at most three and a
literal Boolean XOR query has predicate degree at most two.

To prove this, split each endpoint band into its clean-active and
clean-silent factors, $`E_A=E_{Aa}E_{As}`$ and
$`E_C=E_{Ca}E_{Cs}`$. The actual factors

```math
D_A=K_A^\dagger E_{As}K_A,\qquad
D_C=K_C^\dagger E_{Cs}^\dagger K_C
```

are full-space diagonal C3 under the stated hypothesis. Peel them to
write $`V=D_C V_*D_A`$, where

```math
V_*=K_C^\dagger E_{Ca}^\dagger E_B E_{Aa}K_A.
```

They preserve the clean sectors, so V_* still returns the clean qubit,
and its occupied action is diagonal. Removing silent endpoint axes does
not change J or S.

Every B orbit is a singleton. Its label b therefore satisfies
$`N_A b=0`$. Independence of the clean-active A axes forces b to commute
each such axis individually. Consequently
$`U=E_B E_{Aa}`$ is a C3 phase polynomial on one isotropic Pauli space L.
Its clean image is

```math
S=U(K_AZ_cK_A^\dagger)U^\dagger
 =E_{Ca}(K_CZ_cK_C^\dagger)E_{Ca}^\dagger.
```

For every Pauli with label in $`L^\perp`$, conjugation by S returns that
label, since U commutes that Pauli. Symplecticity gives
$`\mathrm{im}(F_S-I)\subseteq L`$. The independent clean-active C
axes span precisely this displacement image; hence every such axis is
in L. It follows that $`E_{Ca}^\dagger U`$ is a commuting-axis C3 phase
polynomial. Endpoint Clifford multiplication shows that **V_* is C3 on
the entire physical space**. This conclusion follows from its actual
common-axis factorization, not from its occupied logical restriction.

A C3 unitary commuting Z_c induces a logical C3 unitary on c=0. Indeed,
a logical Pauli lifts to a physical Pauli commuting Z_c; its conjugate is
a physical Clifford commuting Z_c, whose restriction is a logical
Clifford. The occupied V_* is therefore diagonal C3. Restoring the
occupied diagonal C3 factors D_A,D_C preserves diagonal C3 membership
and proves the claim, with every literal phase retained.

No new circuit, workspace, or native operation was introduced by this
proof. The hypothesis may be weakened to requiring the two clean-silent
endpoint factors to have diagonal occupied actions, since they are
already clean-preserving C3 operators. Individual Z-type pullbacks are
a transparent sufficient condition.

The silent-band condition is substantive. The retained
[paid cubic query](THREE_LAYER_NATIVE_STRUCTURE.md#8-paying-the-third-band-an-exact-cubic-xor-query)
has a three-layer quartic occupied sign. Its Clifford-conjugated Pauli
bands E,E' are Z-diagonal; the remaining endpoint and interband Cliffords
act only on non-clean wires. In that valid factorization, J=S=Z_c and
the middle pair count is zero. Thus m=0 alone does not even
force an occupied real diagonal action to be C3. Its silent endpoint
factors do not satisfy the sufficient diagonal hypothesis above. The
non-Pauli pair-free case likewise cannot discard its silent factors
without a separate argument.

## 5. Complete-family extension witnesses

Neither large clean-active rank nor a complete commuting logical family
repairs the failed full-space extension. The following literal family
also has no middle pairs and no shared Pauli clean flip.

For every integer $`k\ge2`$, use wires $`c,1,\ldots,k-1`$, with only c
initialized. Put, in operator order with the rightmost gate acting first,

```math
\begin{aligned}
P&=\prod_{i=1}^{k-1}X_i,&
K&=\left(\prod_{i=1}^{k-1}\mathrm{CX}_{c,i}\right)H_c,\\
E&=T_c\prod_{i=1}^{k-1}T_i,&
E_0&=T_P T_{Z_cP}^\dagger,\\
U_A&=EK,& U_B&=EK E_0.
\end{aligned}
```

The two axes of E_0 are physically independent and commute. Its literal
clean-zero block is I and its clean-one block is $`\exp(-i\pi P/4)`$.
Both encoders therefore agree on every occupied input, and their common
clean image is

```math
S=E(X_cP)E^\dagger
 =\prod_{j=c,1,\ldots,k-1}\frac{X_j+Y_j}{\sqrt2}.
```

Its intrinsic active rank is k. It defines one half-dimensional code,
not k separate clean conditions. The identity
$`K E_0K^\dagger=T_P T_{X_c}^\dagger`$ exhibits U_B as two actual
native bands with endpoint Clifford K.

For every logical generator Z_i, its first image is the Pauli
$`C_i=U_AZ_iU_A^\dagger=Z_cZ_i`$. Write

```math
\begin{aligned}
J_i&=(T_P T_{X_c}^\dagger)(Z_cZ_i)
                           (T_P T_{X_c}^\dagger)^\dagger,\\
A_i&=U_BZ_iU_B^\dagger=EJ_iE^\dagger.
\end{aligned}
```

Each J_i is an actual Hermitian Clifford, since it is a one-band Pauli
image. Direct signed conjugation gives

```math
J_iZ_cJ_i=Y_c,\qquad
A_iZ_cA_i=\frac{Y_c-X_c}{\sqrt2}.
```

The latter is non-Pauli, so **every logical generator A_i is non-Clifford
on the full space**. Nevertheless, for $`\Pi=(I+S)/2`$,
$`A_i\Pi=C_i\Pi`$. All A_i commute with one another and S, and the
same occupied equality holds for every logical product. Even products
agree on the full space; odd products retain the defect. Thus neither
high rank nor agreement of the whole family supplies the missing
Clifford membership.

There is also no shared Pauli clean flip. A physical Pauli with a Pauli
preimage under U_A must be Z-type, since independent outer T axes give
distinct noncancelling terms otherwise. Write it as
$`H=Z_c^s\prod_iZ_i^{z_i}`$ and put $`b=\sum_i z_i\pmod2`$.
Its preimage under K has X_c exponent $`s+b`$. Its commutation bits
with E_0's two axes P and Z_cP are b and s. Its preimage under U_B is
therefore Pauli precisely when $`b=s=0`$; every such H commutes S.
Clean flips require $`s+b=1`$, which is excluded. Here the outer E
axes span precisely the clean-active space. S fixes their binary labels
individually while negating their signed Z representatives, so the middle
clean-pair count is zero.

### 5.1. A paid easy target and actual inverses

The complete transition is $`U_B^\dagger U_A=E_0^\dagger`$. To compile
its two axes, use

```math
D=\left(\prod_{j=2}^{k-1}\mathrm{CX}_{1,j}\right)
          H_1\mathrm{CX}_{1,c},\qquad
E_0=D(T_c^\dagger T_1)D^\dagger.
```

D and K each cost k elementary H/CNOT gates. Both actual inverses use
the same counts, and no S expansion changes the ledger.

| Full native word | T count | T-depth | Elementary H/S/CNOT |
| --- | ---: | ---: | ---: |
| Displayed encoder comparison | $`2k+2`$ | at most 3 | $`4k`$ |
| Actual inverse of that comparison | $`2k+2`$ | at most 3 | $`4k`$ |
| Simplified E_0-dagger baseline | 2 | 1 | $`2k`$ |
| Actual inverse of the baseline | 2 | 1 | $`2k`$ |

There are k physical wires, exactly one clean, and no additional work.
The occupied target is identity on every arbitrary data input and
reference. It needs no gates if the unused extension is unspecified.
The one-band baseline retains that extension literally; its clean-one
block is $`\exp(i\pi P/4)`$, so its forward word must not be mistaken
for its actual inverse. Error and leakage are zero. This is a
counterexample to two proposed proof inferences, not a hard-query family
or an obstruction in the unrestricted target model.

The [exact fixtures](../../tests/test_occupied_logical_transport.py)
retain all logical products, the full unused blocks, native phases,
shared-flip exclusion, active ranks, and separately emitted inverses.
The arbitrary-width statements follow from the identities above.
