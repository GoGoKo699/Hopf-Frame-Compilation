# A four-address field orbit with a copyable Clifford seed

[Field shear](FIELD_SHEAR_INTERTWINER.md) ·
[Complete copy grading](DIRTY_COPY_GRADING.md) ·
[Encoded copy contract](ENCODED_DIRTY_COPY_INTERFACE.md)

The two-bit field query admits a noncommuting, exactly copyable orbit
of Hermitian Clifford operators. Its seed has a literal 26-Clifford-gate
word, with no T gates or private helpers. A paid physical copy completes
the four-row indicator with 288 T gates, at most 96 T bands and seven
returned dirty wires. The complete grading follows
from two nondegenerate Pauli algebras selected by a central physical
Pauli. In particular, the monomial and full-Hadamard seed failures do
not exclude every Clifford seed at this address size.

This is a finite complete-operator pilot. Its four addresses do not give
an all-size source or improve the general frame-depth frontier.
The later [flat-seed pilot](FLAT_FIELD_BANK_PILOT.md) gives a different
commuting bank and a cheaper 72-T, 28-band direct indicator with three
dirty wires. The noncommuting central-sector construction here remains
a distinct algebraic and native calibration.

## 1. The physical seed and query

Use four arbitrary carrier wires, numbered 0 through 3. The field-query
banks are $`u=(u_0,u_1)`$ on wires 0,1 and $`v=(v_0,v_1)`$ on wires
2,3. With address $`x\in\mathbb F_2^2`$, use

```math
Q_x|u,v\rangle=|u,v+M_xu\rangle,\qquad
M_0=0,\quad M_1=I,\quad
M_2=\begin{pmatrix}0&1\\1&1\end{pmatrix},\quad M_3=M_1+M_2.
```

This is exactly the existing query pencil over the degree-two field
extension, with modulus $`t^2+t+1`$. Address labels 0,1,2,3 are binary,
and address addition below is XOR. Define

```math
A=\frac{-X_1Z_0+Z_1Z_3}{\sqrt2},\qquad
B=\frac{X_0Z_1Z_2+Z_0}{\sqrt2},\qquad R=AB.
                                                               \qquad(1)
```

The Pauli summands in each factor anticommute. The two summands of A
commute with both summands of B, so A and B commute. Each factor is a
Hermitian Clifford involution; hence so is R. This description retains
all literal signs, not just projective Pauli labels.

For a direct native word, specify the following Clifford walls in
chronological order:

```text
K_A: Z(1), CZ(0,1), CX(3,1).
K_B: CZ(0,1), CZ(0,2).
```

Their literal Pauli actions give

```math
A=K_AH_1K_A^\dagger,\qquad B=K_BH_0K_B^\dagger.
```

Thus one chronological seed word is
$`K_B^\dagger,H_0,K_B,K_A^\dagger,H_1,K_A`$.
Expand every CZ as H,CX,H and every Z as S,S. Each conjugated-H
factor costs 13 elementary Clifford gates, for 26 in total. The
actual adjoint has the same count after reducing Z adjoints by
$`Z^\dagger=Z`$. Neither seed word initializes, measures or resets
any carrier qubit.

Put $`R_x=Q_xRQ_x^\dagger`$. Every address-conditioned query is a
Clifford, so all four orbit members are Hermitian Cliffords, even
though the coherent addressed query has its separately paid native
T cost.

## 2. A central Pauli and the complete word algebra

Let

```math
P=Z_0X_2.
```

Direct Pauli multiplication gives $`[R_x,P]=0`$ for all x and the
following literal commutators. All orbit members square to I.

| Pair $`(x,y)`$ | $`R_xR_yR_xR_y`$ |
| --- | --- |
| 0,1 | $`-P`$ |
| 0,2 | $`I`$ |
| 0,3 | $`P`$ |
| 1,2 | $`P`$ |
| 1,3 | $`I`$ |
| 2,3 | $`-P`$ |

These finite identities can be checked directly from (1) using the
shear's generator rules: X on the u bank adds the corresponding
column of $`M_x`$ to the v-bank X label, while Z on the v bank
adds the corresponding row to the u-bank Z label. X on v and Z on
u are fixed. Their products retain the displayed signs.

Since P is central, every word in the R's and their actual inverses
reduces to

```math
\pm P^e R_0^{s_0}R_1^{s_1}R_2^{s_2}R_3^{s_3},
\qquad e\in\{0,1\},\quad s\in\mathbb F_2^4.       \qquad(2)
```

This is an all-word reduction, not a cutoff in word length.
The two P eigenspaces have dimension eight. On the positive sector,
the anticommuting pairs are $`(R_0,R_1)`$ and $`(R_2,R_3)`$.
On the negative sector, they are $`(R_0,R_3)`$ and $`(R_1,R_2)`$.
All other distinct pairs commute within the respective sector.
Both resulting binary commutation matrices are nonsingular.

For every nonzero s, some generator therefore anticommutes with
the ordered word in (2), separately on each P sector. Its trace on
each sector is zero. Multiplication by P or a scalar sign does not
change that conclusion. Every word of nonzero occurrence parity
consequently has trace zero on the complete physical carrier.

The [exact grading criterion](DIRTY_COPY_GRADING.md#1-the-exact-full-space-criterion)
therefore supplies a fixed unitary C with

```math
C(R_x\otimes I_Y)C^\dagger=R_x\otimes X_{Y_x}
\qquad(x=0,1,2,3),                                  \qquad(3)
```

for an arbitrary four-qubit physical output bank Y. The criterion
includes every carrier and output input and all external references.
Tensoring identity on additional arbitrary helpers preserves it.

More explicitly, each P sector carries the full two-qubit Pauli
algebra $`M_4`$ with multiplicity two. The generated algebra is
$`M_4\oplus M_4`$, represented with that multiplicity in each block.
Its dimension is 32; each of the 16 parity components has dimension
two. This is a noncommuting copyable family, rather than an independent
commuting virtual bank.

## 3. A native copy in the physical field coordinates

Write $`\Pi_\pm=(I\pm P)/2`$. The commutator table supplies the four
Hermitian involutions

```math
\begin{aligned}
J_0&=\Pi_+R_1+\Pi_-R_3,&J_1&=\Pi_+R_0+\Pi_-R_2,\\
J_2&=\Pi_+R_3+\Pi_-R_1,&J_3&=\Pi_+R_2+\Pi_-R_0.
\end{aligned}                                                   \qquad(4)
```

Each satisfies $`J_jR_xJ_j=(-1)^{[x=j]}R_x`$. Apply controlled
$`J_0,J_1,J_2,J_3`$ in that chronological order, with controls the
four output wires, and place Hadamards on all output wires before
and after this word. In output-X sector f the carrier receives the
actual ordered product $`J_3^{f_3}J_2^{f_2}J_1^{f_1}J_0^{f_0}`$.
Its conjugation gives all required signs simultaneously, proving (3).
The J's need not commute and the resulting copy C is not an involution;
its actual reversed native word is required.

### Paying the sector controls

The following fixed Clifford K is a useful control coordinate change;
its gates are again listed chronologically:

```text
SWAP(0,1), SWAP(2,3), CX(3,1), CZ(0,1), CZ(0,2), Z(1), Z(2).
```

It has 17 elementary Clifford gates after expanding SWAP into three
CNOTs and using the CZ/Z macros above. Direct conjugation gives

```math
K(H_0H_1)K^\dagger=R,\qquad K^\dagger P K=-X_3.
```

Let $`K_x=Q_xK`$. The central P commutes with the field queries, so
$`K_x^\dagger P K_x=-X_3`$ and
$`K_x(H_0H_1)K_x^\dagger=R_x`$. After applying $`K_x^\dagger`$
and $`H_3`$, the desired positive-P sector is the ordinary bit-3
value one. A negative-P control adds X on that bit before the core
and its actual inverse after the core. All coordinate changes are
then reversed.

For the core, write $`T(A)=\exp(i\pi(I-A)/8)`$ and
$`U=T(Y_0)T(Y_1)`$. Its literal native implementation obeys
$`UZ_0Z_1U^\dagger=H_0H_1`$. Compute parity of bits 0,1 onto bit 1,
apply CCZ with controls the physical output and bit 3 and target bit 1,
and undo the parity. Conjugating this controlled-$`Z_0Z_1`$ operation
by U gives the required doubly controlled $`H_0H_1`$.

The eight-T, two-band dirty CCZ uses arbitrary carrier bit 2 as a
mask and returns it exactly before the surrounding coordinate changes
are undone. Together with U and its actual inverse the core uses
12 T gates in four bands. No new helper is introduced. The exact
[emitter and core fixtures](../../tests/test_field_copy_clifford_native.py)
retain the two CCZ parity bases, all scalar phases and actual inverses.

Each controlled J in (4) consists of two such sector factors. The
full physical copy has eight factors, with the independently counted
forward and actual-inverse bounds

```math
T(C)=96,\qquad D_T(C)\le32,\qquad G_{\rm all}(C)=792.
                                                               \qquad(5)
```

Here $`G_{\rm all}`$ includes both Clifford and T-type gates. The
copy uses precisely the four carrier and four output wires, and
no initialized workspace. It need not return the carrier at its own
boundary; complete return follows from the final echo.

## 4. Three-cube coherent query and the complete indicator

In the v-bank Fourier coordinates, the two-bit shear is the cubic
phase with exponent

```math
f=x_0u_0v_0+x_0u_1v_1+x_1u_1v_0+x_1u_0v_1+x_1u_1v_1.
```

It has the exact factorization over the binary field

```math
f=x_0u_0(v_0+v_1)+x_1u_1v_0+v_1(x_0+x_1)(u_0+u_1). \qquad(6)
```

Apply the [returned cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md)
to these three products of linear forms, with a distinct arbitrary
mask per product. In every Gray-code phase round the three T-type
gates act on different masks. Computing each linear form directly
onto its mask and closing each walk returns all three masks. There
are 24 T-type gates, 30 CNOTs and four external v-bank Hadamards:

```math
T(Q)=24,\qquad D_T(Q)\le8,\qquad G_{\rm all}(Q)=58.
                                                               \qquad(7)
```

The actual inverse has the same resources and complete action.
Every address, carrier and mask input is arbitrary. No field-element
state or clean parity wire is supplied.

Let $`\mathcal A=QRQ^\dagger`$ with coherent address Q. The
chronological word
$`\mathcal A^\dagger,C^\dagger,\mathcal A,C`$
implements the complete indicator and identity on its private work.
Combining (5), (7) and the 26-Clifford seed gives one finite native word
and its actual inverse with

```math
T=288,\qquad D_T\le96,\qquad G_{\rm all}=1868,
\qquad b_{\rm dirty}=7.                              \qquad(8)
```

The allocation is two logical address wires, four physical output
wires, four arbitrary returned carrier wires and three arbitrary
returned query masks: 13 physical wires in total. The seven dirty
wires include the carrier, not just the three query masks. The copy
reuses a carrier wire internally only with complete return of that
wire at the CCZ boundary. Every carrier and mask returns after the
full echo, with literal phases preserved on all input columns and
external references.

### The improved query also benefits the direct decoder

The existing [four-wire decoder](FIELD_SHEAR_INTERTWINER.md#6-a-paid-four-wire-decoder-for-two-address-bits)
has 63 T gates, at most 36 bands and 193 total gates per direction.
Its direct indicator uses four decoder calls, two coherent queries
and one physical X, paid as H,S,S,H. Substituting (7) therefore gives
another complete native word, with 300 T gates, at most 160 bands,
892 total gates and only the three returned query masks as dirty work.
No separate carrier bank is used in that direct construction.

| Four-address construction | T gates | T bands at most | Total gates | Dirty wires |
| --- | ---: | ---: | ---: | ---: |
| Previously recorded direct decoder and five-cube query | 332 | 160 | 944 | 5 |
| Direct decoder with the three-cube query (7) | 300 | 160 | 892 | 3 |
| Noncommuting seed/copy with the same query | 288 | 96 | 1868 | 7 |

These are finite upper bounds, not optimality claims. Against the improved
direct decoder, the copy construction saves 12 T gates and 64 bands,
using four additional dirty wires and more total gates. Both new words
and their actual inverses are included in the retained emitter and ledger
fixtures. A scalable result still needs a uniform seed/copy pair on
arbitrary address lengths.
Moreover, the [symplectic-orbit theorem](CLIFFORD_SEED_COPY_STRUCTURE.md)
prevents extending the pilot's Pauli-valued commutator structure to
p at least seven at copy error below 1/64. A Clifford extension must
change that symplectic structure, not only the size of its central sectors.
The general upper theorem R96 is unchanged.

## 5. Exact evidence

[The algebra fixtures](../../tests/test_field_copy_clifford_pilot.py)
verify the literal 26-gate seed and its inverse on every carrier
column, all six central commutators, both sector forms, and every
nonzero-parity ordered product with and without P. These finite
checks agree with the all-word reduction above.

[The native fixtures](../../tests/test_field_copy_clifford_native.py)
check every column of the doubly controlled core and all eight sector
factors, their actual inverses, and all sixteen output-X blocks of
the full copy. Equality on these complete blocks includes off-sector
coherences. The query and its actual inverse are tested on all 512
address/carrier/mask columns. Exact gate and band checks cover the
common 13-wire copy allocation and the nine-wire direct-decoder allocation,
with both directions of each indicator word.
An additional assembled-word check exercises every address with nonzero
dirty inputs in both actual native directions, guarding the register wiring.
The complete echo follows from these full-space component identities
and (3), not from a basis-state-only copy promise or a finite search
for an asymptotic circuit.
