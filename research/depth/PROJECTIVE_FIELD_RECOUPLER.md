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
