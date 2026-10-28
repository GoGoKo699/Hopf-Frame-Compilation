# A jointly encoded dirty carrier and its copy interaction

[Field-shear source](FIELD_SHEAR_INTERTWINER.md) ·
[Universal source gauges](WEAKER_INDICATOR_INTERFACES.md) ·
[Table-specific probes](TABLE_SPECIFIC_INDICATOR_INTERFACES.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

A fixed encoding on an entirely arbitrary dirty carrier need not be
physically synthesized when it cancels from the completed operation.
This gives a different construction target for the field-shear route:
jointly implement a source seed and a copy interaction, rather than the
whole recursive decoder. The physical output bank remains unencoded.
The interaction is a genuine missing native operation. In particular,
a fixed Clifford copy cannot serve the one-query field orbit at growing
size, even with an arbitrary non-Clifford seed and arbitrary helpers.

This is a sufficient complete-operator interface and a restriction on
one companion class. It supplies no fast seed/copy pair and no improved
general depth theorem. The existing gauge and table-syndrome chapters
own the corresponding classifications with a fixed physical probe;
here the carrier and its interaction are designed together.

## 1. Cancel an encoding only on the returned dirty carrier

Let x have p bits, let $`m=2^p`$, and let A be a w-qubit arbitrary
carrier. A may include encoded spectator qubits. Let Y be a separate
arbitrary m-qubit physical output bank. The desired operation is

```math
I_p|x,Y\rangle=|x,Y+e_x\rangle.
```

Suppose a fixed unitary E identifies m virtual bank qubits inside A,
with identity spectators, and define

```math
R_a=E^\dagger(X_a\otimes I)E,
\qquad R_x=Q_xR_0Q_x^\dagger.                         \qquad(1)
```

The second equality is a required source identity. For the field route,
$`Q_x`$ is the paid affine shear; the all-size decoder proves that an
abstract E with this property exists. Existence does not pay a native
word for its seed $`R_0`$ or any interaction with A.

The canonical encoded copy would be

```math
C_0=(E^\dagger\otimes I_Y)
       \left(\prod_a\mathrm{CX}(A_a;Y_a)\right)
       (E\otimes I_Y).
```

Its useful property is only the following action. A native candidate C,
independent of x, may instead satisfy

```math
C(R_x\otimes I_Y)C^\dagger=R_x\otimes X_{Y_x}
\quad\text{for every }x.                              \qquad(2)
```

C need not equal $`C_0`$, preserve the carrier at its own boundary,
or be a computational copy. It may mix its carrier and arbitrary work.
Every such register must be included in (2), with identity in the ideal
comparison where appropriate. The output Y is a physical logical bank,
so an encoding on Y cannot be silently omitted.

Put $`\mathcal A=\sum_x|x\rangle\langle x|\otimes R_x`$.
Equation (2), on the complete space, gives

```math
C\mathcal A C^\dagger\mathcal A^\dagger
   =I_A\otimes I_p.                                  \qquad(3)
```

A chronological native word is
$`\mathcal A^\dagger,C^\dagger,\mathcal A,C`$.
The actual inverse is essential; C need not be an involution. Equation
(3) returns all of A and every helper on arbitrary inputs and references.
There is no initial application of E. Equivalently, conjugating the
whole ideal word by E changes nothing because its final action on A is
identity. The abstract encoding has canceled, not been granted a free
implementation. Every remaining physical interaction must still be paid.

In fact (2)--(3) need no virtual encoding or commutativity assumption:
any directly specified unitary family $`R_x`$ and C satisfying (2)
give the same complete indicator. Equation (1) is one constructive way
to search for that family.
The [exact parity-grading criterion](DIRTY_COPY_GRADING.md) characterizes
which unitary families admit such a C and gives robust word tests. It also
exhibits a noncommuting all-size family; an abstract copy implementer still
needs its own native synthesis.

## 2. Native resources, error and entry into the frame transfer

Suppose the coherent Q, fixed seed R and copy C have common forward
and actual-inverse bounds $`t_Q,t_R,t_C`$ for T count,
$`d_Q,d_R,d_C`$ for T depth, and $`g_Q,g_R,g_C`$ for total gates.
The source word $`Q R Q^\dagger`$ and echo (3) give

```math
T_I\le4t_Q+2t_R+2t_C,\qquad
D_T(I)\le4d_Q+2d_R+2d_C,
```
```math
G_I\le4g_Q+2g_R+2g_C,\qquad
b_I\le w+h.                                         \qquad(4)
```

Here $`b_I`$ excludes x and the specified output Y. The quantity h
counts all additional live private wires in one common allocation.
Private pools may be reused only at boundaries where their full return
has been proved; otherwise their widths must be added. All parts of A
are charged, even though A is returned. No clean input is used.

For approximation, let $`\widetilde{\mathcal A}`$ and
$`\widetilde C`$ be actual unitaries with complete operator errors
$`\epsilon_A,\epsilon_C`$ from a pair satisfying (2). Four-factor
telescoping, using the actual adjoints, gives

```math
\|\widetilde C\widetilde{\mathcal A}
  \widetilde C^\dagger\widetilde{\mathcal A}^\dagger
       -I_A\otimes I_p\|
 \le2\epsilon_A+2\epsilon_C.                         \qquad(5)
```

For exact Q and a seed error $`\epsilon_R`$,
$`\epsilon_A\le\epsilon_R`$. If Q is approximate with complete
error $`\epsilon_Q`$, the bound is
$`\epsilon_A\le2\epsilon_Q+\epsilon_R`$.
Alternatively one can certify the error of (2) directly for the actual
unitaries; multiplying by their actual source adjoint gives exactly
that error in (3). No number of carrier bits multiplies these bounds.

Thus at each fixed requested accuracy, a jointly paid pair with
$`T,b=O_\tau(mp^{O(1)})`$, $`G=O_\tau(m^2p^{O(1)})`$ and
source-plus-copy depth $`o(\log(p+2))`$ supplies the discovery
interface for polynomial-overhead fusion. The resulting physical
indicator then enters the existing R91 query assembly and R98 transfer.
Those proofs already charge output batching, traversal, actual inverses
and full error; no unencoded table wall is inserted into the carrier.
The same log-log term and fixed-precision qualifications remain.

## 3. The missing companion is coherent dual readout

For the canonical copy, put $`P_f=E^\dagger Z^fE`$ for an m-bit
mask f, including identity spectators. Hadamards on the physical Y bank
turn $`C_0`$ into

```math
H_Y C_0 H_Y
 =\sum_f|f\rangle\langle f|_Y\otimes P_f.             \qquad(6)
```

The universal coherent interaction, including every superposition of
f and every carrier input, is the missing component. Cheap separate
words for selected classical masks do not establish (6) at the required
count, depth or accuracy.

For one fixed table f, the weaker relation
$`P_fR_xP_f=(-1)^{f_x}R_x`$ suffices for its scalar phase echo.
Using it for an XOR output also requires paid coherent control, for
example a controlled $`P_f`$ word.
The full R91 bilinear construction uses coherent table-dependent
combinations, so this scalar observation does not replace its complete
interface. The linked table-specific chapter owns the resulting
syndrome and row-query reductions.

For a fixed full intertwiner E, multiplying it on the left by a gauge
commuting with every $`X_a`$ leaves every $`E^\dagger X_aE`$
unchanged. Such a gauge may simplify a word for E, but does not redesign
the seed. Right multiplication by a unitary commuting with every
$`Q_x`$ can change that source family. Neither freedom supplies the
native copy in (2).

## 4. A Clifford copy has a uniform error floor

Let R be any Hermitian involution on the full carrier, including any
arbitrary helpers, and let $`Q_x`$ be Cliffords whose binary symplectic
entries have Boolean degree at most d. Put $`R_x=Q_xRQ_x^\dagger`$
and $`M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj`$.
For every fixed, address-independent Clifford C on carrier and Y,

```math
\max_x\|C(R_x\otimes I)C^\dagger-R_x\otimes X_{Y_x}\|
 \ge\sqrt{2-2\sqrt{M_d(p)/2^p}}.                     \qquad(7)
```

Expand $`R=\sum_v r_v\tau(v)`$ in Hermitian Paulis, with
$`\sum_v r_v^2=1`$. For a fixed v, the output-bank Pauli label
after conjugation by C is a vector-valued degree-d function
$`F_v(x)`$. Its values lie in the binary span of its at most
$`M_d(p)`$ ANF coefficient vectors. It can therefore equal the
distinct labels $`(e_x,0)`$ for at most $`M_d(p)`$ addresses:
these output X labels are linearly independent. Clifford conjugation
is a signed bijection of Pauli labels, so no coefficients merge.

Let $`w_x=\sum_v r_v^2[\,F_v(x)=(e_x,0)\,]`$. Then
$`2^{-p}\sum_xw_x\le M_d(p)/2^p`$. The target operator in (7)
lies wholly in that output-label sector and has normalized
Hilbert--Schmidt norm one. Cauchy--Schwarz bounds its overlap with the
actual operator by $`\sqrt{w_x}`$. If the maximum operator error
is eta, normalized Hilbert--Schmidt distance consequently gives

```math
1-\eta^2/2
 \le2^{-p}\sum_x\sqrt{w_x}
 \le\sqrt{M_d(p)/2^p},
```

proving (7). This permits C to mix every carrier, output and helper
wire; it does not permit dependence on the external address x.
For one field shear $`d=1`$, the ratio is $`(p+1)/2^p`$.
Exact copying is impossible for $`p\ge2`$, and the displayed error
floor tends to $`\sqrt2`$. For q literal affine queries with fixed
Clifford carrier interlayers and the external address held as a
parameter, $`d\le q`$. This is a companion-class screen,
not a native-depth lower bound for a jointly non-Clifford construction.

Even a fixed-table Pauli companion is limited. The Pauli labels v for
which $`R\tau(v)R=(-1)^{\sigma(v)}\tau(v)`$ form a binary linear
space, and sigma is a linear functional there, by multiplying the
conjugation identities. For a one-shear orbit, the label of
$`Q_x^\dagger P_fQ_x`$ is affine in x. The required sign is therefore
affine as well. A nonaffine table cannot have a physical Pauli companion,
regardless of the seed's Pauli spread. Non-Pauli companions remain open.

## 5. Two stronger screens for candidate pairs

The [output-axis span bound](COPY_AXIS_SPAN.md) permits a general
non-Clifford copy. For a degree-d Clifford query orbit, every absorbed
Pauli-rotation description of the copy must have
$`r_X\ge m(1-\eta^2/2)^2-M_d(p)`$ at error $`\eta\lt\sqrt2`$,
where $`r_X`$ counts independent physical X directions in its projected
axis span. For one field shear, exact copying needs at least
$`m-p-1`$ such directions and T gates. This excludes localized
corrections but is compatible with the linear count budget; it gives
no growing depth bound.

Independently, the [monomial-seed screen](MONOMIAL_QUERY_SEED_OBSTRUCTION.md)
gives copy error at least one half for p at least two whenever the seed
is monomial in a simultaneous query eigenbasis. More generally it
suffices that the original query algebra commute with its seed-conjugated
image. The copy may be completely arbitrary and non-Clifford. The proof
uses a four-word identity and includes arbitrary unitary seeds and
helpers. A seed must therefore mix the query eigenspaces beyond this
class, as well as admit the extensive physical copy interaction.

## 6. Exact finite evidence and next construction

`tests/test_encoded_dirty_copy_interface.py` checks a native one-bit
source with a non-Pauli seed, a noncanonical copy interaction, its
actual inverse and the full returned-work indicator. A separate
support fixture tests the Clifford-copy obstruction with a genuinely
non-Clifford seed. These are finite contract checks, not an all-size
source/copy compiler.

The [native-copy chapter](NATIVE_COPY_INTERACTION.md) extends the Clifford
screen to an extensive T-count floor and a dual-support requirement. It
supplies a one-band non-Clifford component, but excludes the specified
product-magic/output-only mixing and coded two-band field-copy families
at growing size. These classes do not exhaust native one-band copies.
The separate [affine Pauli-orbit screen](PAULI_ORBIT_COPY_SCREEN.md) gives
a one-quarter error floor with Clifford dual probes, including native
one-band copies, for p at least two. The field shear is not a Pauli query.

The [four-address Clifford pilot](FIELD_COPY_CLIFFORD_PILOT.md) now pays
a noncommuting seed and physical copy in the same field coordinates, giving
a complete zero-clean finite indicator. Its central Pauli commutators also
illustrate the [symplectic-orbit restriction](CLIFFORD_SEED_COPY_STRUCTURE.md):
at p at least seven, a projectively involutive Clifford seed with commuting
symplectic orbit actions has copy error at least 1/64, even with arbitrary C.

The all-size task remains a jointly paid pair satisfying (2) within (4).
The later [flat four-address bank](FLAT_FIELD_BANK_PILOT.md) supplies a
genuinely non-Clifford seed, a 12-T six-band physical copy, and a cheaper
direct indicator through its six-T three-band decoder. Extending that
complete probe decoder uniformly is another concrete route; its finite
word does not settle either all-size interface.
Full intertwiners must also meet the
[Pauli-domain requirement](FIELD_DECODER_PAULI_DOMAIN.md), including with
moving Clifford helpers, and the
[disjoint-span two-band screen](TWO_BAND_FIELD_DECODER_OBSTRUCTION.md).
These conditions do not restrict a decoder required only to transport the
complete physical probe. For the Clifford-dual alternative, an associative
algebra preserved by scalar-field shears has a
[four-address relation](SCALAR_FIELD_RADICAL_COPY_SCREEN.md) for every
address embedding; changing the field basis does not repair that model.
For a Clifford seed, noncommuting symplectic dynamics are therefore necessary
at arbitrarily small fixed error; non-Clifford seeds are another possibility.
The complete dual readout, native inverse and every arbitrary-input sector
must still be proved. Global dirty-carrier cancellation removes the abstract
decoder only after the physical companion is paid.
