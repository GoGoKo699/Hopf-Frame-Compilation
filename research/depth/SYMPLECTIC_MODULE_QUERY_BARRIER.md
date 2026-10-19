# A source-count barrier for symplectic commutator readout

[Full-input algebra extraction](ALGEBRA_QUERY_EXTRACTION.md) ·
[Earlier translation-carrier bound](SCALAR_REPRESENTATION_PRODUCT.md#5-a-source-cost-barrier-survives-pauli-carrier-compression) ·
[Conditional query interface R98](QUERY_DEPTH_FRAME_TRANSFER.md)

Compressing a full-table extraction into an arbitrary linear Pauli
carrier does not remove its separately paid source cost. A degree-d
matrix coefficient of an elementary abelian symplectic representation
requires at least $`2^{d-1}`$ carrier qubits and that many T gates in a
source producing the representation by Pauli conjugation. The carrier
need not be an algebra quotient, a regular representation or a specified
Pauli orbit.

The same numerical barrier was already known for the narrower Clifford
translation-carrier architecture. The result here removes that orbit
hypothesis. It remains a source-architecture theorem. Two independent
half-address representations escape its single-representation premise;
section 5 gives a completely charged witness and identifies its remaining
selector cost. No unrestricted query lower bound or faster R98 query is
claimed.

## 1. The readout contract and theorem

Let $`V=\mathbb F_2^{2w}`$ be a w-qubit Pauli-label space with its usual
nondegenerate alternating pairing. Let

```math
S:\mathbb F_2^a\longrightarrow\mathrm{Sp}(V),\qquad
S_{x+y}=S_xS_y,\qquad S_0=I
```

be any binary symplectic representation. No choice of computational
coordinates, carrier subspace, native axes or quotient is imposed.
For fixed probe labels u,v, consider the Boolean function

```math
h(x)=\langle u,(S_x-I)v\rangle.                         \qquad(1)
```

This is exactly the scalar exposed by a nested Pauli commutator. If C_x
is any Clifford lift of S_x, and the Pauli probes use their actual
unitary adjoints, then

```math
K_x=C_xP_vC_x^\dagger P_v^\dagger,\qquad
K_xP_uK_x^\dagger P_u^\dagger=(-1)^{h(x)}I.           \qquad(2)
```

K_x is a Pauli up to its literal scalar. That scalar cancels with its
actual adjoint in the second commutator. Thus (1) does not borrow an
address-dependent phase from the Clifford lift. Conversely, a different
readout involving an additional signed phase or nonlinear function of the
tableau is not covered merely because it also uses Cliffords.

**Theorem.** If h has Boolean degree $`d\ge1`$, then

```math
w\ge2^{d-1}.                                          \qquad(3)
```

For each direction in a chosen degree-d monomial, after restricting the
other address coordinates to zero, one also has

```math
\mathrm{rank}(S_{e_i}-I)\ge2^{d-1}.             \qquad(4)
```

Suppose, in addition, a fixed physical source F and a linear-label Pauli
family P_x satisfy $`C_x=FP_xF^\dagger`$ with every C_x Clifford on the
entire physical carrier. Then every exact Clifford+T implementation of
F, including any extra arbitrary returned helpers, obeys

```math
T(F)\ge2^{d-1}.                                       \qquad(5)
```

It suffices that the relevant input Pauli subgroup has Clifford images;
F need not be C3 on all probes and its inverse need not be C3. The P_x
may have their fixed Pauli multiplication cocycle: conjugation preserves
it, and the symplectic matrices still form the displayed representation.
A source specified only on an initialized sector is outside this theorem.

## 2. A free cyclic module is forced by the top coefficient

Restrict to the coordinates of a degree-d monomial and relabel them
$`1,\ldots,d`$. Put $`N_i=S_{e_i}-I`$. The representation law gives

```math
N_i^2=0,\qquad N_iN_j=N_jN_i,\qquad
S_x=\prod_{i=1}^d(I+x_iN_i).
```

The top coefficient of (1) is one, so
$`\langle u,N_1\cdots N_dv\rangle=1`$. In particular its vector
$`N_{[d]}v`$ is nonzero. All $`2^d`$ vectors

```math
N_Av:=\prod_{i\in A}N_i v,\qquad A\subseteq[d]
```

are linearly independent. To prove this, take a nonzero relation and
choose a supported A of minimum cardinality. Multiply the relation by
$`N_{[d]\setminus A}`$. A term B not contained in A vanishes because a
square-zero factor is repeated. A proper subset B of A has zero
coefficient by minimality. Only the A term remains, giving the
contradiction $`N_{[d]}v=0`$.

Thus their span is a free cyclic module for
$`\mathbb F_2[e_1,\ldots,e_d]/(e_1^2,\ldots,e_d^2)`$, embedded in V.
It has dimension $`2^d`$, proving (3). On this independent basis, N_i
maps each A not containing i to $`A\cup\{i\}`$ and kills all remaining
vectors. Its rank on this subspace is $`2^{d-1}`$, proving (4).

The proof requires neither regularity of the full carrier nor a prior
algebra-quotient description. It also applies whenever any linear
functional of S_x has degree d: the corresponding top matrix coefficient
is nonzero, so choose v on which that matrix is nonzero. In particular,
XORing any fixed collection of probe readouts for the same representation
does not bypass (3)–(4).

## 3. The actual source cost, including approximation scope

A Clifford C with symplectic matrix S has exactly
$`2^{\mathrm{rank}(S-I)}`$ nonzero Hermitian-Pauli coefficients,
all of the same magnitude in normalized Hilbert–Schmidt coordinates.
The retained [Pauli-support proof](SCALAR_REPRESENTATION_PRODUCT.md#52-a-physical-t-gate-can-at-most-double-pauli-support)
applies without a translation-orbit premise: fixed labels confine support
to one coset of $`\mathrm{im}(S-I)`$, and the Clifford intertwining
relation makes the whole coset occur.

Conjugating one Pauli through F, with t physical T/T-dagger gates,
produces at most $`2^t`$ terms. Every T conjugation splits a term into
at most two, and every Clifford only permutes terms. Applying this to
$`FP_{e_i}F^\dagger=C_{e_i}`$ and using (4) proves
$`t\ge\mathrm{rank}N_i\ge2^{d-1}`$. The count here is T(F),
not twice T(F): one conjugation through a physical T gate doubles the
support once. Returned arbitrary helpers tensor the target image with
identity and change neither its support nor its displacement rank.

The same retained support argument covers full-operator source
approximation. If a native W is within delta in operator norm of F tensor
helper identity, with $`\delta\lt 1/2`$, then

```math
T(W)\ge2^{d-1}+\log_2(1-4\delta^2).
```

For $`\delta\lt 1/(2\sqrt2)`$, integer T-count is still at least
$`2^{d-1}`$. This is a guarantee about the whole source on every input.
Small error of a final commutator alone does not imply it.

For a full-degree a-address table, the eight separately implemented
source calls of the [algebra extraction](ALGEBRA_QUERY_EXTRACTION.md)
therefore require at least $`8\cdot2^{a-1}=4Q`$ T gates, where
$`Q=2^a`$. The inverse has the same minimum count by native reversal.
Even an optimal compressed carrier and an optimal implementation of each
complete source cannot meet R98's $`O(\sqrt Q)`$ count at square-root
or larger dirty width by this route. This strengthens the earlier
algebra-quotient diagnosis to arbitrary linear Pauli carriers under (1).

This lower bound must not be charged after a compiler fuses and
resynthesizes the entire scalar word. Such a circuit need not implement
any of the separately specified sources. The theorem gives no growing
T-depth lower bound and no lower bound for a different query mechanism.

## 4. The binary carrier bound is sharp

The dimension and rank statements are sharp as binary representation
claims. On the $`2^d`$-dimensional square-zero algebra R_d, define

```math
\langle b,c\rangle=[e_{[d]}](bc),\qquad
S_xb=\left(\prod_{i=1}^d(1+x_ie_i)\right)b.
```

Complementary basis subsets pair to one, so the form is nondegenerate.
It is alternating: every square b squared is its constant coefficient,
which has zero top coefficient for positive d. Multiplication by each
unit preserves the pairing because its square is one. Thus this is a
symplectic representation on exactly $`w=2^{d-1}`$ qubits' label space.
With $`u=v=1`$, (1) is precisely the degree-d conjunction, and every
N_i has rank $`2^{d-1}`$.

This witnesses the binary bound, not an optimal native source. No
simultaneous Clifford lift, Pauli-conjugating source or synthesis with
exactly that T-count is inferred from the matrices alone.

## 5. Two half-address families escape the premise

Split the address into u,v with lengths p,q. Put $`m=2^p`$ and
$`n=2^q`$, so $`mn=Q`$. For an arbitrary m-by-n truth table D, use
m arbitrary carrier qubits A_i and n arbitrary carrier qubits B_j.
Define

```math
\gamma_u=X_{A_u},\qquad
\delta_v=Z_{B_v}\prod_i Z_{A_i}^{D_{iv}}.
```

Their literal commutator is $`(-1)^{D_{uv}}I`$ on every carrier input.
Both families have explicit Clifford translations. Let P_u^A and P_v^B
be the wire permutations adding the half-address to the corresponding
bank index. Let

```math
C_D=\prod_{i,j:D_{ij}=1}\mathrm{CX}(A_i;B_j).
```

These cross-bank CNOTs commute; its actual inverse is their reversed
word. It uses exactly wt(D) elementary Cliffords, at most Q. Then

```math
\gamma_u=P_u^A X_{A_0}(P_u^A)^\dagger,\qquad
\delta_v=C_D P_v^B Z_{B_0}(P_v^B)^\dagger C_D^\dagger.
```

The two symplectic translation representations generally do not commute
with each other. Their joint product is not the representation required
in section 1. This is a real escape from that premise: for a balanced
split the carrier has $`m+n=O(\sqrt Q)`$ qubits, even for a full-degree
table. These half-translation families need not have the extra source
form $`FP_xF^\dagger`$; for example, a nontrivial wire permutation can
have nonzero trace, whereas a conjugated nonidentity Pauli has zero trace.
Thus (5) is not being applied to these sources. No claim about an
unavailable free wire translation is made.

Here is the complete query interface that would use these families.
Let U_A and U_B be actual coherent controlled wire-translation sources
on their respective half-address and bank, each with its own supplied
native word, actual inverse and arbitrary-helper-return certificate.
They may change those addresses internally but return them at each
completed source call. They do not include the original output z.
Build the controlled moving probe

```math
A_z=U_A\mathrm{CX}(z;A_0)U_A^\dagger
```

and the moving phase probe

```math
B=C_DU_B Z_{B_0}U_B^\dagger C_D^\dagger.
```

Chronologically apply $`A_z,B,A_z^\dagger,B^\dagger`$, with Hadamards
on z at both ends. Its action is the literal XOR of $`D_{uv}`$ into z,
with both carrier banks, both half-addresses and all helpers returned.
The controlled probe uses the arbitrary z itself, not a copied control.
Let $`T_A,G_A,d_A,b_A`$ and $`T_B,G_B,d_B,b_B`$ be common bounds
for each supplied source word **and its actual native inverse**. This
matters because a positive-S alphabet need not give equal raw Clifford
counts under reversal. The full native ledger is

```math
\begin{aligned}
T&\le4T_A+4T_B,\\
G&\le4G_A+4G_B+4\mathrm{wt}(D)+8,\\
D_T&\le4d_A+4d_B,\\
B_{\mathrm{extra}}&\le m+n+\max(b_A,b_B),\qquad a_{\mathrm{clean}}=0.
\end{aligned}                                                   \qquad(6)
```

The four C_D occurrences, four top-level source inverses, two controlled
probe CNOTs, two $`Z=S^2`$ words and two output Hadamards are all included.
The source-helper pool is separate from both carrier banks and reused
only after complete return. The actual inverse has the same displayed
ledger after adjacent S-power reduction. Identity actions for zero
half-addresses are included in the full-source promises, not omitted
from a coherent controlled call.

There is an existing paid but insufficient source baseline: execute one
batch of controlled disjoint swaps per half-address bit. Literal Fredkin
words give $`T_A,G_A=O(pm)`$, $`d_A=O(p)`$, and no extra work; similarly
for B. With a balanced split this has $`O(a\sqrt Q)`$ T count and O(a)
depth, although its total Clifford count is O(Q) and carrier width is
$`O(\sqrt Q)`$. The native batches share bank wires and cannot be declared
simultaneous merely because the ideal permutations commute.

The retained [R91 dirty-indicator mechanism](PARALLEL_BILINEAR_WORD_QUERY.md)
already implements the same scalar table query with
$`O(\sqrt Q)`$ T count, O(Q) Cliffords and logarithmic address depth.
That construction already supplies the appropriate fair baseline. The
first missing operation here is a complete selector or probe-transport
implementation with the same count and width and sublogarithmic depth;
rewriting the table as a commutator does not pay it. Even such a
one-output improvement would need a joint multi-output schedule satisfying
all R98 word-query budgets, with shared addresses charged correctly.

## 6. Evidence and research consequence

The [exact fixtures](../../tests/test_symplectic_module_query_barrier.py)
check the free-module ranks and top coefficient in
mixed symplectic coordinates, a sharp binary carrier, the general
Clifford-support law on small matrices, and every physical column of a
native two-half query and its actual inverse. A wider binary fixture
exhibits a full-degree table on a carrier smaller than (3), together with
noncommuting half-address representations, preventing misuse of the bound.

The source compression question has a definite answer under (1): changing
the carrier basis or module cannot supply the missing R98 counts. A
constructive continuation must change the readout/source arrangement,
fuse the native word without separately realizing that source, or solve
the smaller selector task in section 5. The unrestricted depth question
remains open.
