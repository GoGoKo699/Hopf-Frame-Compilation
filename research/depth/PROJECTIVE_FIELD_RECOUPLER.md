# A projective-line recoupler and a paid eight-address decoder

[Recursive field decoder](RECURSIVE_FIELD_SHEAR_DECODER.md) ·
[Four-address flat decoder](FLAT_FIELD_BANK_PILOT.md)

The top recursive field recoupler can be organized by projective lines instead
of treating every nonzero field input as a separate Clifford fiber. A common
unconditional Clifford removes every selected Fourier transform. The remaining
selected maps are binary-linear permutations and displacements. The zero fiber
is specified by the emitted product; it needs no separate nonzero flag.

For three address bits this gives a complete native eight-wire intertwiner,
with every input arbitrary and no additional helper. It is a finite paid
calibration of the recursive construction. Enumerating the projective lines is
exponential in bank width, so the construction does not improve R96 or supply
the uniform depth bound required by the frame transfer.

[Section 6](#coupled-six-helper-allocation) gives a second complete allocation:
148 T gates in at most 39 bands, with six additional returned arbitrary helpers.
It preserves the zero-extra-helper allocation below as a distinct width/depth
tradeoff; neither finite construction is a uniform source-depth improvement.

<a id="projective-fiber-identity"></a>

## 1. Homogeneity of the complete conditional cell

Use the [recursive tower](RECURSIVE_FIELD_SHEAR_DECODER.md#1-a-compatible-tower-and-address-pencil).
At level p let $`F=K_{p-1}`$, $`K=F[\theta]`$,
$`\theta^2+\theta+\delta=0`$, and write
$`U=d_0+d_1\theta`$. Put

```math
\nu(U)=d_0^2+d_0d_1+\delta d_1^2,\qquad q_0(U)=d_0d_1.
```

For nonzero U, the linear parts of the two controlled coordinate maps are

```math
A_U=\begin{pmatrix}
d_0+d_1&d_0+(\delta+1)d_1\\
d_1/\nu&(d_0+d_1)/\nu
\end{pmatrix},\qquad
B_U=\begin{pmatrix}d_0&\delta d_1\\d_1&d_0+d_1\end{pmatrix}.
```

Let $`L_M`$ denote the literal permutation $`|v\rangle\mapsto|Mv\rangle`$,
and let H denote the trace-pairing Fourier transform on the first F coordinate.
The complete middle cell of the recoupler is

```math
K_U=X_{f_U}C_U,\qquad
C_U=L_{B_U}H L_{A_U},\qquad
f_U=B_U\binom{0}{q_0(U)/\nu(U)}.                    \qquad(1)
```

This includes the affine offset, which cannot be omitted when compiling the
literal operator. For every $`s\in F^\times`$,

```math
A_{sU}=\mathrm{diag}(s,s^{-1})A_U,\qquad
B_{sU}=sB_U,\qquad f_{sU}=s f_U.                   \qquad(2)
```

The Fourier identity $`H L_s=L_{s^{-1}}H`$ holds literally for the trace
pairing, with no input-dependent or scalar phase. The two scalings in (2)
therefore cancel across H and give

```math
C_{sU}=C_U,\qquad K_{sU}=X_{s f_U}C_U.             \qquad(3)
```

Thus the nonzero inputs have only $`|F|+1`$ Clifford linear fibers, one for each
projective line in $`F^2`$. The displacement is F-linear along each line. This
statement is about the entire conditional operator, not only its symplectic
matrix or its action on a promised state.

<a id="common-baseline-and-zero-fiber"></a>

## 2. A common baseline removes selected Fourier transforms

Choose line representatives $`\mathcal D=\{1\}\cup\{\theta+a:a\in F\}`$.
For D=1, $`B_1=I`$ and $`f_1=0`$, so $`C_1=H L_{A_1}`$. Apply $`C_1`$
unconditionally. On the line represented by $`D\ne1`$, the required linear
correction is

```math
C_D C_1^\dagger
 =L_{B_D}H L_{A_D A_1^{-1}}H.                      \qquad(4)
```

Use the fixed linear change of the existing U register
$`U/D=s+t\theta`$. Its line predicate is exactly $`t=0`$. For each nonbaseline
D, apply the following complete factor and then undo that change of U:

1. Apply H unconditionally on the first F coordinate of V.
2. Apply $`L_{A_D A_1^{-1}}`$ only when $`t=0`$.
3. Apply H unconditionally again.
4. Apply $`L_{B_D}`$ only when $`t=0`$.
5. Add $`s f_D`$ to V only when $`t=0`$.

On a disabled line the two Fourier walls cancel. On an enabled line the factor
is $`X_{s f_D}C_D C_1^\dagger`$. Every nonzero U enables exactly one factor,
or none when U lies on the baseline line. Equation (3) then proves the required
cell on every nonzero input, including every V input and literal phase.

The zero input belongs to all lines. There s=0, so all displacements vanish,
and the zero fiber is the explicitly ordered product

```math
K_0^{\rm emitted}
 =\left(\prod_{D\in\mathcal D\setminus\{1\}}^{\rm chronological}
          C_D C_1^\dagger\right) C_1.             \qquad(5)
```

Here chronological order means later operators multiply on the left. Equation
(5) is a fixed Clifford, not an identity promise. The source queries are the
identity on the entire U=0 sector. After the recoupler's fixed coordinate sum,
the target pairs satisfy $`a=c`$; the lower queries and the half swap are also
the identity on that entire intermediate sector. Consequently any unitary
zero fiber, including (5), preserves the complete intertwining proof. There
is no zero detector, clean flag, discarded sector or unreturned workspace.

Let $`\widehat R_p`$ be this cell with the same external Fourier and coordinate
walls as the recursive recoupler. Its complete identity remains

```math
\widehat R_p Q^{(p)}_{(x',b)}\widehat R_p^\dagger
 =(Q^{(p-1)}_{x'}\otimes Q^{(p-1)}_{x'})S^b.        \qquad(6)
```

The inverse is the reversed complete word, including the inverse of (5).
It need not agree with the original recoupler on U=0.

<a id="borrowed-quartic-native-word"></a>

## 3. A sixteen-T quartic phase with one returned borrow

The retained four-T relative Toffoli is the literal unitary

```math
R(a,b;h)=(-i)^{ab}X_h^{ab}.
```

Let $`U=H_hR(a,b;h)H_h=(-i)^{ab}Z_h^{ab}`$ and
$`V=R(c,d;h)=(-i)^{cd}X_h^{cd}`$. Their chronological commutator
$`U,V,U^\dagger,V^\dagger`$ is

```math
V^\dagger U^\dagger VU=(-1)^{abcd}I_h.             \qquad(7)
```

The two relative phases cancel, and the remaining sign is the Pauli X/Z
commutator. This holds for both arbitrary values of h and all coherent
control inputs. Hadamards on d turn (7) into a three-control X on d.
The emitted word has **16 T gates, at most eight physical T bands and 46 total
elementary gates**. The actual inverse has the same counts. It uses only
H, S, CNOT, T and T-dagger, and h returns literally. The p=3 construction
borrows an unused logical coordinate for h; it allocates no additional wire.

## 4. The finite p=3 emitter

For p=3, $`F=\mathbb F_4`$, $`K=\mathbb F_{16}`$ and the nonbaseline
representatives are $`D=4,5,6,7`$ in the retained little-endian tower basis.
A direct selected-CNOT compilation of (4) is possible, but the four linear
corrections admit a smaller joint reflection factorization.

Use local V-wire labels
$`P_{a+16b}=i^{\mathrm{wt}(a\mathbin\&b)}X^a Z^b`$ and
$`T(P)=\exp(i\pi(I-P)/8)`$. For commuting independent A and B that both
anticommute with P, define the Hermitian Clifford reflection

```math
J(A,B;P)=T(A)T(B)P T(B)^\dagger T(A)^\dagger.        \qquad(8)
```

Anticommutation rewrites J as P times the two quarter-turn Pauli rotations
$`\exp(i\pi A/4)`$ and $`\exp(i\pi B/4)`$, which are Clifford.
Hermiticity and the involution property follow from its displayed conjugation.

The following literal identities specify each correction. Each row lists the
earlier and later reflection; the final column means
$`C_D C_1^\dagger=\sigma J_1J_0`$.

| D | Earlier $`(A,B;P)`$ | Later $`(A,B;P)`$ | $`\sigma`$ |
| --- | --- | --- | ---: |
| 4 | (19, 50; 32) | (156, 228; 112) | +1 |
| 5 | (156, 228; 64) | (19, 50; 32) | +1 |
| 6 | (89, 247; 1) | (19, 50; 13) | +1 |
| 7 | (19, 50; 32) | (89, 247; 87) | -1 |

There is no residual Pauli correction. The finite identities can be verified
directly by multiplying the displayed signed Paulis or comparing their
Clifford action and one literal column; the exact fixture compares all 16
columns with the field matrices in (4).

For each pair the emitter constructs one joint Clifford K taking two distinct
physical Z axes to A and B, with their positive signs. One explicit procedure
first diagonalizes A, pulls B through that inverse, then diagonalizes the
remaining independent Pauli while fixing the first Z axis; an X corrects a
negative pulled sign. It emits
$`K^\dagger,[T,T],K`$, so the two T rotations occupy one physical band.
The words for K and its actual inverse contribute their elementary Clifford
gates to the ledger.

The line predicate has two binary zero controls. A Clifford-conjugated
relative Toffoli gives the selected Pauli P with the common phase
$`(-i)^{ab}`$. Conjugating it by the joint two-T word gives a selected J,
using **eight T gates in four bands**. The two selected J factors contribute
$`(-1)^{ab}`$ together. A CZ on the two controls corrects this for D=4, 5 and
6; for D=7 it is exactly the required negative sign. These phases remain
literal throughout the emitted word.

Thus the four linear corrections cost 64 T gates and 32 bands. The three
nonzero displacement maps use two three-control X gates each, paid by (7),
adding 96 T gates and 48 bands. Every such gate borrows an existing logical
coordinate and returns it. The baseline, line-coordinate maps and exterior
Fourier walls are Clifford. No oracle, initialized flag or controlled gate
remains in the elementary emitted alphabet.

| Complete native word | T count | T bands | Forward total gates | Actual inverse gates | Physical wires | Extra helpers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| $`\widehat R_3`$ | 160 | 80 | 1,104 | 1,104 | 8 | 0 |
| $`\widehat E_3=(E_2\otimes E_2)\widehat R_3`$ | 172 | 83 | 1,234 | 1,286 | 8 | 0 |

Here $`E_2`$ is the full six-T intertwiner from the
[four-address pilot](FLAT_FIELD_BANK_PILOT.md), including its final Clifford
completion. The two identical half-bank words run in parallel, contributing
12 T gates and three bands. Using the complete intertwiner rather than only
its bank-decoding prefix gives, with literal scalar phase one,

```math
\widehat E_3 Q^{(3)}_x\widehat E_3^\dagger=W_x
\quad\text{for every three-bit address }x.         \qquad(9)
```

All eight bank inputs are arbitrary. The internal borrowed logical
coordinates return before their next use; the actual complete inverse is
paid separately. Equation (9) therefore also holds on arbitrary coherent
addresses and external references.

## 5. Verification, scope and the next operation

The [native fixture](../../tests/test_projective_field_recoupler.py) checks
the quartic primitive on every five-wire column and its actual inverse,
the projective scaling identity, all paired reflection factorizations and
their simultaneous signed native axis maps, all 240 nonzero-U columns against the original recursive
recoupler, and every column of both actual native inverses. It checks (9)
for all eight addresses and all 256 bank columns, including the modified
zero fiber. Integer arithmetic in $`\mathbb Z[\omega]/2^d`$,
$`\omega=e^{i\pi/4}`$, retains literal phases; normalization and overflow
bounds are explicit. These checks do not use a numerical tolerance.

This removes the unpaid primitive at this finite size. It does not claim
that this is a better general eight-address indicator than the repository's
other constructions. The projective-fiber count at level p is

```math
|F|+1=2^{2^{p-2}}+1=2^{m/4}+1.
```

Enumerating these factors is exponential in m and outside the polynomial
discovery budgets. Separately composing levels also retains a stage cost
at each level. The useful new target is the joint projective-line action in
(4), including its allowed zero-fiber freedom: a uniform coupled synthesis
must compress these factors or change their representation while retaining
the full operator identity and paid work return. Neither a constant-depth
source nor a new frame-depth theorem is established here.

<a id="coupled-six-helper-allocation"></a>

## 6. Coupled reflections and a six-helper allocation

The four projective corrections can share one reflection selector. Their six
displacement terms can also occupy one private-helper phase batch. Both changes
retain the complete literal operator, with a newly specified zero fiber.

### 6.1 One common reflection, with its selector phase paid

Use the same Hermitian Pauli labels and reflection (8), and abbreviate

```math
\begin{aligned}
J_A&=J(19,50;32),\\
J_{B,64}&=J(156,228;64),& J_{B,112}&=J(156,228;112),\\
J_{C,44}&=J(89,247;44),& J_{C,87}&=J(89,247;87).
\end{aligned}
```

The following signed identities are exact:

```math
\begin{aligned}
C_4C_1^\dagger&=J_{B,112}J_A,&
C_5C_1^\dagger&=J_AJ_{B,64},\\
C_6C_1^\dagger&=J_AJ_{C,44},&
C_7C_1^\dagger&=-J_{C,87}J_A.                    \qquad(10)
\end{aligned}
```

In the third identity, changing the earlier C-reflection probe to 44 permits
the same A-reflection with probe 32 in every row. These identities follow by
multiplying the signed Pauli matrices; the fixture independently compares them
with the field trace-Fourier matrices. Every J is a Hermitian Clifford
involution, with the same independent commuting-axis premises as in Section 4.

Let $`f_D(U)`$ indicate the line represented by D, including its zero point.
Write the four binary U coordinates as $`u_0,u_1,u_2,u_3`$. Their common selector is

```math
g(U)=\bigoplus_{D=4}^7 f_D(U)=u_2\lor u_3.        \qquad(11)
```

For nonzero U exactly one nonbaseline line is selected, or none on the baseline
line. At zero all four line predicates are one, and g is zero. Emit the line-5
B reflection and line-6 C reflection first, the common A reflection next, then
the line-4 B and line-7 C reflections. On every nonzero line this chronological
order selects exactly the required pair in (10).

The common A selector needs its literal relative phase. For $`P=P_{32}`$,
which is the single physical $`Z_{V_1}`$, one has

```math
(iP)^{u_2}(iP)^{u_3}(-iP)^{u_2u_3}
   =i^{u_2\lor u_3}P^{u_2\lor u_3}.              \qquad(12)
```

The first two factors are Clifford controlled-iP words. The last is a
Clifford-conjugated four-T relative Toffoli. Conjugation by the common two-T
preparation from Section 4 produces $`i^gJ_A^g`$, at eight T gates in four
bands. Each of the four remaining selected reflections costs the same eight
T gates and four bands and contributes $`(-i)^{f_D}`$. A CZ in the line-7
coordinates supplies its required $`(-1)^{f_7}`$ sign. Thus, on nonzero U,

```math
i^g(-i)^{\sum_{D=4}^7f_D}=1.                     \qquad(13)
```

The linear corrections therefore cost **40 T gates in 20 bands**. All line
coordinate changes and their actual inverses are included. On the entire
zero fiber the chosen linear cell is the fixed Clifford

```math
-J_{C,87}J_{B,112}J_{C,44}J_{B,64}C_1.           \qquad(14)
```

The four relative phases multiply to one there; the minus sign in (14) comes
from the line-7 CZ. Section 2 proves that this full zero-sector action is
permitted. No zero column or scalar sign is discarded.

<a id="private-mask-quartic-batch"></a>

### 6.2 Sixteen-band quartic batching with arbitrary private helpers

The [R91 private-mask word](DIRTY_CUBIC_PHASE_COMPILER.md) gives
$`F(a,b;h)=(-i)^{ab}Z_h^{ab}`$ using four T gates on an arbitrary private
helper h. Its affine version is useful here. Write
$`a=\ell_a\oplus\kappa_a`$, $`b=\ell_b\oplus\kappa_b`$, with binary
linear forms $`\ell_a,\ell_b`$ that do not involve h. Visit helper parity
vertices $`0,\ell_a,\ell_a\oplus\ell_b,\ell_b`$ and return to zero.
At the vertex indexed by a subset J of the two forms, use T with orientation
$`(-1)^{|J|+\sum_{j\in J}\kappa_j}`$. The resulting literal operator is

```math
\widetilde F(a,b;h)
 =i^{\kappa_a\kappa_b}(-i)^{ab}Z_h^{ab}.          \qquad(15)
```

The scalar factor in (15) is essential: it is the constant term of the affine
finite difference. Equivalently, the T-phase exponent is
$`2\kappa_a\kappa_b-2ab+4hab`$ modulo eight. This identity permits repeated,
dependent, overlapping and constant affine forms, without a clean assumption.

Set $`U=\widetilde F(a,b;h)`$ and
$`V=H_h\widetilde F(c,d;h)H_h`$. The same chronological commutator as in (7)
gives

```math
V^\dagger U^\dagger VU=(-1)^{abcd}I_h.             \qquad(16)
```

Both affine scalar factors cancel against their actual inverses. Assign a
distinct private helper to each of Q quartic terms and run corresponding
Gray-vertex T gates in parallel. All T targets in a band are distinct;
principal-to-helper parity CNOTs leave the principal inputs unchanged.
Shared principal controls therefore cause no interference between helper
actions. For L equal to the sum of the four linear support sizes over all
terms, the forward and actual-inverse bounds are

```math
T=16Q,\qquad D_T\le16,\qquad G=20Q+4L,
\qquad b_{\rm dirty}=Q,\qquad b_{\rm clean}=0.    \qquad(17)
```

The gate count includes sixteen T gates and four Hadamards per helper, plus
four visits to each linear support. Both values of every helper have scalar
identity action in (16), so the return identity also holds for arbitrary
helper entanglement and external references. This parallel allocation uses
private masks instead of the compact eight-band logical borrow in Section 3.

For the p=3 displacement, encode an affine form by its binary mask and constant.
U occupies bits 0 through 3, and Fourier-V occupies bits 4 through 7. The six
quartic products are specified by the following four-form rows:

| Term | First form | Second form | Third form | Fourth form |
| --- | --- | --- | --- | --- |
| 1 | (15, 1) | (5, 1) | (4, 0) | (16, 0) |
| 2 | (15, 1) | (5, 1) | (8, 0) | (32, 0) |
| 3 | (14, 1) | (7, 1) | (13, 0) | (144, 0) |
| 4 | (14, 1) | (7, 1) | (6, 0) | (224, 0) |
| 5 | (6, 1) | (11, 1) | (15, 0) | (160, 0) |
| 6 | (6, 1) | (11, 1) | (5, 0) | (240, 0) |

These are the two binary zero predicates for each of the three nonzero-shift
lines, its scalar coordinate, and its Fourier-V displacement pairing. Their
linear supports sum to L=60. The phase batch has 96 T gates, 16 bands and
360 total gates. Hadamard walls on all four V coordinates convert it into
the complete XOR displacement, giving **368 total gates** in either direction.

All displacements can be moved after the coupled linear word: at nonzero U
only one line is active, while at U=0 every displacement vanishes. Thus this
reordering also retains the specified zero fiber. Six arbitrary helpers at
physical wires 8 through 13 suffice and return literally.

### 6.3 Complete allocation and its verification

Including the original Fourier/coordinate walls and the two parallel full
four-address intertwiners gives the following second allocation:

| Complete native word | T count | T bands | Forward total gates | Actual inverse gates | Physical wires | Extra helpers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Coupled $`\widehat R_3`$ | 136 | 36 | 980 | 984 | 14 | 6 |
| Coupled $`\widehat E_3`$ | 148 | 39 | 1,110 | 1,166 | 14 | 6 |

Both words have zero approximation error and zero clean inputs. The actual
inverse reverses the entire word, including the zero-cell Clifford and private
helper batches. Its elementary-gate count is kept distinct because the native
alphabet emits S-dagger with positive S powers. The complete identity is

```math
\widehat E_3(Q_x^{(3)}\otimes I_6)\widehat E_3^\dagger
   =W_x\otimes I_6.                              \qquad(18)
```

The [coupled native fixture](../../tests/test_coupled_projective_recoupler.py)
checks the literal affine constant in (15), all five-wire quartic columns and
actual inverses with affine shifts, and all 3,072 private-helper sector matrices
for the six-term batch. The latter checks both columns of each helper for every
principal input and both actual orientations; tensor independence proves the
complete fourteen-wire identity without constructing a dense matrix of that
size. Independent field and signed-Pauli formulas check (10), all 240 nonzero
linear-cell columns, and every zero-cell column in (14). Integer arithmetic in
$`\mathbb Z[\omega]/2^d`$ then checks every complete logical column, both actual
inverses, all eight addresses in (18), and every physical T-band target list.

This does not replace the 172-T, 83-band zero-extra-helper allocation. At a
general field size the common selector is the nonzero predicate on
$`\dim_{\mathbb F_2}F=m/4`$ binary coordinates, with that growing Boolean degree.
The finite two-bit OR in (11) does not pay it at bounded depth, and the line
list remains exponential. Coupling the recursive levels still requires a
separate uniform construction within the stated budgets.

The [augmented recoupler](AUGMENTED_FIELD_RECOUPLER.md) supplies such an
all-size native cell after explicitly changing the helper query. It retains
linear depth across levels; the [sector-degree theorem](FIELD_SECTOR_DEGREE_BARRIER.md)
limits compression by fixed-frame helper-affine stages. Neither result
changes the two identity-helper finite allocations established here.
