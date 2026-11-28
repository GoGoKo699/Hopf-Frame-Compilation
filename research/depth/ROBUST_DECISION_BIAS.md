# Quantitative decision bias at three native T bands

[Decision source](PRESCRIBED_SCALAR_TABLES.md#one-clean-decision-source) ·
[Two-band spectral screen](FOUR_LAYER_SCALAR_GAP.md#two-band-signed-gap-sources) ·
[Raw decision count](BOOLEAN_QUERY_FRAME_REDUCTION.md#raw-decision-count)

A uniformly random Boolean table cannot have a constant average signed bias
from an unrestricted circuit with at most three physical T bands.
The statement allows arbitrary Clifford walls, including physical address
motion, and places no bound on either width or gate count. It already holds
when the dirty bank is averaged in the computational basis. The complete
one-clean decision condition implies this weaker bias requirement.
The [four-band decision screen](#four-band-complementary-spans)
covers complementary axis spans even at large middle crossing rank.
The [complete-frame lower bound](#complete-return-four-band-floor)
uses the initialized isometry and its actual adjoint to exclude every
four-band word for a typical complete deepest-corner frame.
A [two-table complete frame](#complete-return-two-table-frame)
extends this floor to any polynomial clean count.
A [rare three-cycle frame](#zero-clean-six-band-floor) requires
six bands at small fixed error in the zero-clean model.

## 1. Statement and normalization

Fix $`m\ge1`$ computational input bits to $`x\in\mathbb F_2^m`$, and average the
other $`w`$ inputs uniformly in their computational basis. With a fixed
Hermitian output Pauli $`S`$ and a native word $`V`$, write

```math
b_V(x)=\mathsf{Tr}\left[
 (|x\rangle\langle x|\otimes I/2^w)V^\dagger S V
 \right].
```

This is a mathematical witness rather than a circuit assumption about dirty
inputs. A promise on every dirty column implies the corresponding average
promise. In particular, with one clean input, the promised bias is tested by
setting its coordinate of $`x`$ to zero. The output Pauli can be chosen with
the circuit.

**Three-band entropy bound.** For $`0\lt\epsilon\le1/2`$, the set of these
bias functions for all native words with at most three T/T-dagger bands
has an $`\epsilon`$-net in normalized $`L^2`$ of size at most

```math
2^{\mathrm{poly}(m+\log(1/\epsilon))}.
```

The polynomial is absolute and has no width dependence. With at most two
bands, the exponent can be taken to be
$`O(m^2(m+\log(1/\epsilon))^2)`$.

Consequently there is an absolute polynomial $`P`$ such that, for a
uniformly random $`a`$-bit sign table $`F`$, with probability at least
$`1-e^{-a}`$, every one-clean source with at most three native T bands
satisfies

```math
\left|\mathbb E_x F(x)b_V(x,0)\right|\le P(a)\,2^{-a/2}.
```

In the two-band case $`P(a)=O(a^2)`$ suffices. Here $`(x,0)`$ is the
address and initialized signal. If the source has a uniform signed margin
$`\gamma`$ on every computational dirty input, then the left side is at
least $`\gamma`$. The one-clean coherent sign condition makes the same
implication, even when the output leaves quantum garbage or moves the
address. Neither estimate grows with address size as a depth lower bound.

## 2. Three bands reduce to a single diagonal Choi decoration

Write the actual word as $`V=C_3 E_3 C_2 E_2 C_1 E_1 C_0`$, including
empty bands if needed. With $`W=m+w`$ and
$`\tau=2^{-W}\mathsf{Tr}`$, set

```math
\begin{aligned}
 A&=C_2^\dagger E_3^\dagger C_3^\dagger S C_3 E_3 C_2,\\
 P_z&=C_0(Z(z)\otimes I)C_0^\dagger,\\
 B_z&=C_1 E_1 P_z E_1^\dagger C_1^\dagger,\\
 \beta_z&=\mathbb E_x(-1)^{x\cdot z}b_V(x)
          =\tau(AE_2B_zE_2^\dagger).
\end{aligned}
```

The last equality is the Walsh expansion of the computational input
projector. The two end Paulis become literal Hermitian Cliffords after
transport through their physical T bands. On the normalized
operator Choi space define $`|K\rangle=(K\otimes I)|\Phi_W\rangle`$.
Then

```math
|\psi\rangle=E_{\rm ch}|A\rangle,\quad
E_{\rm ch}=E_2^\dagger\otimes E_2^{\mathsf T},\quad
\beta_z=\langle B_z|\psi\rangle.
```

Here $`|A\rangle`$ and $`|B_z\rangle`$ are stabilizer vectors.
$`E_{\rm ch}`$ is a diagonal T/T-dagger band on their doubled analytical
register. The Choi space supplies neither physical work nor initialization.

### A large overlap gives a common quadratic restriction

A stabilizer vector is flat on its affine support. Its phase there is
an eighth-root scalar times $`i`$ to an integer-linear exponent times a
quadratic sign in affine coordinates. On the affine support $`H`$ of
$`|A\rangle`$, the phase of $`|\psi\rangle`$ has exponent modulo eight
of degree at most three, with degree-two coefficients even and
degree-three coefficients divisible by four. This follows from the
integer parity expansion of each affine coordinate in the signed
T exponent of $`E_{\rm ch}`$.

Put $`N=2^m`$ and $`\alpha=\epsilon/(4\sqrt N)`$. If every
$`|\beta_z|\lt\alpha`$, the zero function already approximates the bias
by Parseval. Otherwise choose $`z_0`$ with $`|\beta_{z_0}|\ge\alpha`$.
Let $`H_{z_0}`$ be the support of $`|B_{z_0}\rangle`$,
$`I_0=H\cap H_{z_0}`$, and let $`d,d',h`$ be the respective
dimensions of $`H,H_{z_0},I_0`$. The overlap is the product of
$`2^{h-(d+d')/2}`$ and the average of the relative phase on $`I_0`$.
It follows that $`d-h\le2\log_2(1/\alpha)`$ and that the phase
average has modulus at least $`\alpha`$.

For the relative exponent $`p`$, three additive differences modulo eight
equal four times the binary trilinear form

```math
T(h_1,h_2,h_3)=\sum_j(e_j\bmod2)\ell_j(h_1)\ell_j(h_2)\ell_j(h_3).
```

Here $`e_j`$ is the signed T exponent on physical coordinate $`j`$ of the
doubled register, and $`\ell_j`$ is the linear part of that coordinate on
the support directions. This formula follows by differing each affine
parity bit three times. The Clifford phases contribute zero third
differences. Three applications of Cauchy–Schwarz, or the
$`U^3`$ norm of the relative phase, give

```math
\alpha^8\le
\mathbb E_{h_1,h_2,h_3\in\mathsf{dir}(I_0)}
                         (-1)^{T(h_1,h_2,h_3)}.
```

The average is the bias of a trilinear form and is nonnegative by
averaging its third argument. Its negative binary logarithm is its
analytic rank. [F56](../../docs/reference/SOURCE_CATALOGUE.md) bounds
partition rank polynomially in analytic rank for tensors of fixed order
over $`\mathbb F_2`$. Here each partition-rank-one term has a linear
factor on one of the three copies of $`\mathsf{dir}(I_0)`$.
Intersect the kernels of these factors in the common coordinate space.
This gives a subspace $`K\subseteq\mathsf{dir}(I_0)`$ on which
$`T|_{K^3}=0`$ and

```math
R=\mathsf{codim}_{\mathsf{dir}(H)}K
 \le 2\log_2(1/\alpha)+\mathrm{poly}\bigl(1+\log_2(1/\alpha)\bigr).
```

On each affine $`K`$-coset in $`H`$, the phase of $`|\psi\rangle`$
is a stabilizer phase, with its original scalar retained. In coset
coordinates, distinct third differences kill the cubic coefficients.
Using $`\Delta_e^2=-2\Delta_e`$ for an order-two direction, two repeated
differences force each quadratic coefficient to be divisible by four.
Three repeated differences force each linear coefficient to be even.
The remaining exponent is a constant plus twice an integer-linear
polynomial plus four times a binary quadratic polynomial.

Choose linear shifts identifying the cosets with
$`t\in\mathbb F_2^R`$. Normalize their restrictions with literal phases
as $`|S_t\rangle`$. The X masks of the stabilizers of $`|A\rangle`$
span $`\mathsf{dir}(H)`$. Indeed the Z-only stabilizers have
dimension $`2W-d`$, and the full stabilizer has dimension $`2W`$.
For each basis shift choose an actual stabilizer Pauli with that mask.
The ordered product $`P'_t`$ of those generators fixes $`|A\rangle`$
and carries the base coset to coset $`t`$. Therefore

```math
|\psi\rangle=2^{-R/2}\sum_t |S_t\rangle,\qquad
 |S_t\rangle=Q_t|S_0\rangle,\qquad Q_t=E_{\rm ch}P'_tE_{\rm ch}^\dagger.
```

Both $`P'_t`$ and the input Paulis $`P_z`$ have linear binary labels in
their respective parameters and fixed multiplication cocycles.

## 3. The large coset overlaps have a common active support

Let fixed Cliffords $`F_0,H_0`$ prepare $`|\Phi_W\rangle`$ and
$`|S_0\rangle`$ from computational zero on the doubled space, including
their literal phases. Set

```math
g(z,t)=\langle B_z|S_t\rangle=\langle0|C(z,t)|0\rangle,\qquad
 C(z,t)=F_0^\dagger(B_z\otimes I)^\dagger Q_t H_0.
```

The binary symplectic entries of $`C(z,t)`$ and its inverse have degree
at most two in $`z,t`$. A single T band conjugating a Pauli gives a
product of independent same-basis Clifford transvections. The X mask
of the Pauli says which occur, so its tableau is affine in that mask
and hence affine in its input parameters. The two such families occur
between fixed Clifford walls in $`C`$.

Put $`\delta=\alpha 2^{-R/2}`$ and let $`\mathcal D`$ be the set of
parameter pairs with $`|g(z,t)|\ge\delta`$. A significant pair's
computational stabilizer support contains zero and has direction
dimension at most $`2\log_2(1/\delta)`$. Its direction is the column
space of the X-output/Z-input symplectic block $`M(z,t)`$.
Lift $`(z,t)`$ to its monomials of degree at most two and take the
affine hull $`L`$ of the lifted $`\mathcal D`$. The empty case
contributes only discarded overlaps. Otherwise choose an affine basis
of $`L`$ from those lifted points. Since $`M`$ is affine in the lifts,
every column of $`M`$ at a lift-in-$`L`$ point belongs to the span
$`R_*`$ of the column spaces at the affine-basis points. Thus

```math
r=\dim R_*
 \le O\bigl((m+R)^2(1+\log(1/\delta))\bigr).
```

A fixed computational linear Clifford preserving zero takes $`R_*`$
into $`r`$ active wires. Include it at the output of $`C`$. The
stabilizer output now factors for every lift-in-$`L`$ pair as a
stabilizer state on the active wires and one deterministic computational
string on spectators. A nonzero spectator string makes the overlap
zero. Outside $`L`$ every overlap has modulus below $`\delta`$.

### Compact signed stabilizer data

A bounded number of parameterized diagonal Clifford phases with
bounded-degree parameter coefficients, between fixed Clifford walls,
has polynomial Pauli labels and signs of bounded degree in the
parameters and in the input Pauli label. The bound is independent of
the physical width. For a fixed Clifford, expand the ordered product of
its signed generator images. Its scalar uses one- and two-generator terms
modulo four, with the phase of the chosen Pauli section included. For a
diagonal phase, translating its exponent by the X mask gives its new Z
mask and its phase at computational zero. Modulo-four addition costs
pairwise carries, and the third bit of an integer sum modulo eight
uses at most fourth-order carries. Fixed walls and a bounded number
of these parameter rounds preserve a constant degree bound.
[Literal Pauli phase arithmetic](CONTROLLED_PAULI_TRANSPORT.md#4-polynomial-phase-arithmetic-independent-of-width)

For each spectator $`Z_j`$, pull back through $`C`$. Its X mask is zero
on the selected lift and its sign is a polynomial in $`z,t`$. The
spectator string is zero exactly when all those sign exponents vanish.
Store a basis of their polynomial span. For an arbitrary active Pauli,
introduce $`2r`$ formal label variables. Pull back that Pauli through $`C`$.
Store a basis of the span of its polynomial X coordinates and its single
sign polynomial. The active expectation vanishes when an X coordinate
survives on input zero. Otherwise it equals the signed value. Those
expectations determine the active state projectively. The polynomials
use only $`m+R+2r`$ variables of bounded degree, so both tests and the
projective state have a polynomial-size description. The next step
retains the scalar beyond this Pauli description.

## 4. Recover the literal overlap phase

### An active correction with a fixed number of variable rounds

Choose an active Clifford $`K_\xi`$ which takes the reconstructed active
stabilizer state to zero up to a scalar. Its actual inverse is known by
the decoder. This correction can use one fixed architecture with
$`O(r^2)`$ formal parameter bits $`\xi`$. In that architecture a bounded
number of diagonal S/CZ layers and Pauli shifts have variable
coefficients. All Hadamards are fixed.

To see this, straighten the stabilizer's affine support by a
computational linear map and a translation. A diagonal Clifford
undoes its relative phase, and partial Hadamards on the straightened
support take it to zero. Every invertible binary matrix $`A_0`$ is a
product of two invertible symmetric matrices. In cyclic module
coordinates $`\mathbb F_2[t]/(p)`$, the coefficient of
$`t^{\deg p-1}`$ in the product of two elements is a nondegenerate
symmetric pairing for which multiplication by $`t`$ is self-adjoint.
Nondegeneracy follows by multiplying a nonzero representative of
degree $`j`$ by $`t^{\deg p-1-j}`$. Taking the direct sum over cyclic
blocks produces symmetric invertible $`B`$ with
$`BA_0=A_0^{\mathsf T}B`$. Then $`A_0=B^{-1}(BA_0)`$ has the required
two factors.

For a symmetric invertible $`D`$, the three alternating Z and X
symplectic shears with matrices $`D,D^{-1},D`$ have product

```math
\begin{pmatrix}0&D^{-1}\\D&0\end{pmatrix}.
```

A full Hadamard on the appropriate side gives the desired in-place
linear symplectic map. Partial Hadamards themselves have a three-shear
factorization using the diagonal mask of the selected subset for each
of the three shears. Each Z shear is a parallel diagonal phase with
the symmetric matrix entries as parameters. Each X shear sandwiches
such a phase by fixed full Hadamards. One additional variable Pauli
can correct the literal Pauli action, since two Cliffords with the
same symplectic map differ by a Pauli up to scalar. This proves the
fixed architecture claim. Its actual small-register unitary is
fully determined by the decoder, including the phase it obtains
from the selected matrices and Pauli signs.

On a lift-in-$`L`$ pair passing the spectator test, let $`\xi`$ be
the selected correction parameters. Then

```math
 K_\xi C(z,t)|0\rangle=\lambda|0\rangle,\qquad
 g(z,t)=\lambda\langle0|K_\xi^\dagger|0\rangle .
```

The second factor is a matrix element on $`r`$ active wires.
The remaining task is to describe the eighth root $`\lambda`$ without
paying a full physical tableau.

### Literal phases from a bounded-round real path sum

Consider a real Clifford word with fixed Clifford walls and a bounded
number of diagonal quadratic sign rounds whose coefficients are
bounded-degree polynomials in formal Boolean parameters $`\theta`$.
On the parameter subset where the word takes $`|0\rangle`$ to
$`\pm|0\rangle`$, its literal sign agrees with one
bounded-degree Boolean polynomial in $`\theta`$. This is a pointwise
version of the [Gaussian path-marker argument](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md#4-pauli-markers-bound-the-inverse-polar-degree).
The subset may be arbitrary.

For the proof, expand the fixed walls into real H, CNOT and diagonal
Clifford gates. With $`h`$ fixed Hadamards, the amplitude from zero to
the all-plus output bra is a quadratic Gauss sum on $`h`$ path bits,
times its fixed normalization. The path quadratic representative
$`Q(\theta)`$ and the constant $`c(\theta)`$ have bounded-degree
coefficients. At a scalar-return point the normalization forces the
polar matrix $`Q+Q^{\mathsf T}`$ to be nonsingular. Put formal Z
markers after the original Hadamards. At those same points the
two-marker mixed phase coefficients are the off-diagonal entries of
the inverse polar matrix. Propagating marker Paulis through their
actual suffix shows that each such coefficient is a dot product of
one propagated X mask with one Z mask. The two masks are
bounded-degree polynomials from the product of a bounded number of
variable symplectic walls. Individual Pauli signs cancel from the
mixed coefficient. This defines one polynomial matrix
$`P(\theta)`$ with zero diagonal whose value is the inverse polar matrix
at every scalar-return point. The
[Arf sign identity](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md#2-a-basis-independent-formula-for-the-arf-sign)
gives the literal sign exponent on this subset as

```math
 c(\theta)+c_2(P(\theta)Q(\theta))+\binom{h/2}{2}\pmod2.
```

Here $`c_2`$ is the second characteristic coefficient. This is a
bounded-degree polynomial in the formal parameters, including at
points where the inverse did not exist. If $`h`$ is odd there are no
real scalar-return points. If $`h=0`$, the constant path phase gives
the sign directly.

To apply the argument to $`K_\xi C(z,t)`$, first retain the physical
phase of each conjugated-Pauli family. Write a Pauli with its ordered
product phase as $`i^{s(v)}X^{u(v)}Z^{v'(v)}`$, where the two labels
are binary linear in $`v`$, and $`s`$ uses the modulo-four cocycle.
For a physical T band $`E`$ with signed orientations
$`e_j\in\{-1,0,1\}`$, the exact expression is

```math
 E P_v E^\dagger=
 \omega^{2s(v)+\sum_j e_ju_j(v)}
 X^{u(v)}Z^{v'(v)}\prod_j S_j^{-e_j u_j(v)},
 \qquad \omega=e^{i\pi/4}.
```

Store the initial eighth-root factor by the polynomial digits of its
exponent. The modulo-eight parity expansion and carries give constant
degree. After removing those known factors, the families and their
inverses use a bounded number of variable Pauli shifts and diagonal
Clifford phase rounds. The variable active correction has this form
as well.

Realify the word using one analytical rebit, as in the
[realification of Clifford intervals](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md#5-complex-clifford-intervals-and-the-approximation-threshold).
Fixed complex Cliffords become fixed real Cliffords. In a variable
S/S-dagger layer let $`q_i`$ denote the S occurrence bit, $`w_i`$ the
data bit, $`a_i=q_iw_i`$, and $`b`$ the rebit. The layer applies
$`J^{\sum_i a_i}`$ on the rebit, where $`J=XZ`$, with one extra
minus sign for each S-dagger on an active data bit. Its phase before
the parity toggle is

```math
 (-1)^{\sum_{i\lt j}a_i a_j+b\sum_i a_i+\sum_{i\ {\rm negative}}a_i}.
```

This uses $`J^2=-I`$ and the parity of the second binary digit of
$`\sum_i a_i`$. The exponent is quadratic in data, with coefficients
of degree at most two in occurrence parameters. The rebit toggle by
$`\sum_i a_i\bmod2`$ is a diagonal quadratic round between two fixed
rebit Hadamards. Variable X shifts can likewise be implemented by
fixed Hadamards around a sign round. Thus the realified word meets
the real path-marker premise with a constant round count.

For each $`j=0,1,2,3`$, append the fixed realification of the scalar
rotation $`\omega^{-j}`$ and store its sign interpolation polynomial.
For a selected corrected complex endpoint, exactly one of these four
rotations leaves the realified rebit in computational zero. That
rotation produces a real endpoint $`\pm|0\rangle`$, because the
complex scalar is an eighth root. Identify the rotation by pulling
back the rebit Z to the all-zero input. Store a basis of the X-mask
polynomials and the sign polynomial for each rotation, as in the
active-stabilizer description. A zero X mask and positive sign gives
the required rotation. The corresponding real path polynomial
determines its literal sign. Undo the rotation and multiply back the
stored eighth-root factors to recover $`\lambda`$. The rebit and
parameter correction are proof and decoding devices, not components
of the original native source.

## 5. Count descriptions and compare to random tables

To approximate the Fourier coefficients, evaluate the exact
$`g(z,t)`$ on lift-in-$`L`$ pairs using the stored data. Set it to zero
outside $`L`$. There every removed term had modulus below $`\delta`$.
Using the normalized coset sum, each $`\beta_z`$ changes by at most

```math
 2^{-R/2}2^R\delta=\alpha .
```

The Fourier coefficients of the original real bias are real.
Take real parts of the reconstructed coefficients before inverting
the Walsh transform. Parseval then bounds the normalized ambient
$`L^2`$ error by $`\sqrt N\alpha=\epsilon/4`$.
The recorded affine hull, dimensions, bounded-degree polynomials and
constraint spans use a polynomial number of bits in
$`m+\log(1/\epsilon)`$. Indeed $`R,r,\log(1/\delta)`$ are such
polynomials, and every formal polynomial variable count is
$`O(m+R+r^2)`$. The degree bounds are absolute even for unbounded
physical width. This proves the claimed net. The decoder is allowed
arbitrary classical computation and needs no circuit synthesis on
the unbounded-width native registers.

Restricting the address-plus-signal function to clean signal zero
multiplies its normalized $`L^2`$ error by at most $`\sqrt2`$. The
approximants have norm at most $`1+\epsilon`$ by the triangle
inequality. For any fixed real approximant $`h`$, independence of
the uniformly random signs gives

```math
 \Pr\left\{\big|\mathbb E_x F(x)h(x)\big|\ge u\right\}
 \le 2\exp\!\left(-\frac{2^a u^2}{2\|h\|_2^2}\right).
```

Use the union bound over the net with
$`\epsilon=2^{-a/2}`$ (or $`1/2`$ for $`a=1`$), and use
Cauchy–Schwarz for the approximation error. The resulting failure
probability is at most $`e^{-a}`$ at a threshold
$`P(a)2^{-a/2}`$ for a suitable absolute polynomial $`P`$.

### A direct complete-frame error floor

The theorem also bounds the original three-band frame depth without
substituting the two-call XOR reduction. Use the [deepest-corner frames](BOOLEAN_QUERY_FRAME_REDUCTION.md#1-a-legal-complete-frame-subfamily),

```math
U_f=\sum_x|x\rangle\langle x|\otimes (XZ)^{f(x)}.
```

Test the target input at zero along with the initialized compiler
helpers, and measure target $`Z`$. This selects columns of the
complete-isometry contract without adding physical initialized work.
Suppose the compiler's full initialized error is $`\eta\le\sqrt2`$.
For every tested input, the real overlap with its ideal output
is at least $`1-\eta^2/2`$, after any allowed common phase is
included. The probability of the correct target outcome is at least
the squared modulus of this overlap. Consequently its averaged
signed bias is at least

```math
2(1-\eta^2/2)^2-1 .
```

Fixing $`h`$ computational input coordinates instead of one multiplies
the normalized $`L^2`$ net error by at most $`2^{h/2}`$. Invoke the
ambient net at accuracy $`2^{-(a+h)/2}/2`$ to keep its logarithmic
cardinality polynomial in $`a+h`$. For any fixed number of clean
compiler helpers, the random-table bound therefore applies to the
tested slice with a polynomial in $`a`$ times $`2^{-a/2}`$.
It also applies when the clean count is any fixed polynomial in $`a`$,
with a polynomial prefactor for that count bound.

Thus for every fixed $`\eta\lt\sqrt{2-\sqrt2}`$ and sufficiently
large $`a`$, almost every deepest-corner frame requires at least four
physical T layers at error $`\eta`$. This allows arbitrary dirty
width, counts and Clifford interlayers, including approximate work
return. It also holds for the complex magnitude family by taking zero
leaf phases. The [complete-return test](#complete-return-four-band-floor)
below strengthens the depth floor to five layers at the same error.
The limiting error floor on this corner family is sharp:
the single target Clifford $`HZ=R_y(\pi/4)`$ is a zero-T
approximation with full operator error $`\sqrt{2-\sqrt2}`$ from
both possible target rotations. The fixed-error lower bound improves
the constant threshold for this family. It leaves the asymptotic
depth order open.

## 6. A shorter two-band proof with a smaller polynomial

With two raw bands, the same Fourier normal form is
$`\beta_z=\tau(AE P_z E^\dagger)`$, where $`A`$ is a fixed
Hermitian Clifford and $`E`$ is the first physical band.
For a Clifford $`D`$ with symplectic matrix $`M_D`$,
its normalized trace is zero or has magnitude
$`2^{-\mathsf{rank}(M_D-I)/2}`$. If its trace is nonzero, it
acts up to fixed Clifford coordinates on at most
$`\mathsf{rank}(M_D-I)`$ effective qubits. The reason is that
all the fixed-label Paulis then commute with the Clifford literally,
and that fixed space contains at least $`W-\mathsf{rank}(M_D-I)`$
independent symplectic pairs.
[Clifford trace and support bound](FOUR_LAYER_SCALAR_GAP.md#native-normal-form-and-the-clifford-trace-bound)

Let $`k=\lceil m+2\log_2(1/\epsilon)\rceil`$ and retain the
affine hull of the Fourier labels with
$`|\beta_z|\ge2^{-k/2}`$. Outside that hull the $`L^2`$ truncation
cost is at most $`2^{(m-k)/2}\le\epsilon`$. If there are no
such labels the zero approximant suffices. Otherwise fix a
significant base point $`z_0`$. The symplectic matrix of
$`EP_z E^\dagger`$ is identity plus the independent physical
Z-axis transvections toggled by $`P_z`$'s X mask. The rank of
the difference of two such matrices is exactly the Hamming
distance between their toggles. For any two significant labels
this rank is at most $`2k`$, by subtracting the two product
matrices for $`AE P_z E^\dagger`$ and using the trace bound.
Consequently all differences on an affine basis about $`z_0`$
toggle a union of at most $`2km`$ T targets. Let $`E'`$ keep only
those targets and $`U`$ be the direction of the retained hull. Then
for $`u\in U`$ the exact coefficients satisfy

```math
 \beta_{z_0+u}=\tau\left(D_0 E' P_u E'^\dagger\right),
 \qquad D_0=A E P_{z_0}E^\dagger .
```

The ordered input Paulis obey
$`P_{z_0+u}=P_{z_0}P_u`$ literally. The discarded physical axes commute
with $`P_u`$, by linearity of the toggle labels. The nonzero trace of
$`D_0`$ localizes it to at most $`k`$ effective qubits. Enlarge that
symplectic support by the axes of $`E'`$ and the at most $`m`$
generators of $`P_u`$, and then add one symplectic partner per
remaining radical direction. This fits in $`O(km)`$ effective
qubits. On that many qubits one can record the actual Clifford
up to its eighth-root scalar, the signed axes, the input Pauli
generators and the label hull in $`O((km)^2)`$ bits.
The decoder computes the normalized trace with their literal signs.

This proves the sharper two-band exponent. In the random
one-clean consequence take ambient tolerance
$`2^{-a/2}/2`$ before restricting the clean input. The width-independent
description has $`O(a^4)`$ bits, giving the stated $`O(a^2)`$
factor in the margin bound.

## 7. A middle-axis rank bound for longer native words

Suppose there are $`L`$ middle T bands after taking the first and
last bands into the input and output Paulis. Move the middle
Clifford walls to one side, absorbing their product into the fixed
outer Clifford. In this common frame let $`2s`$ be the rank of the
binary symplectic form on the span of all middle axes.
A fixed Clifford sends its nondegenerate pairs to $`s`$ active wires
and its radical to spectator Z labels. For two middle bands, $`s`$
is the rank of their cross-band commutation matrix. The pairing
pulled back to the two isotropic spans has the cross-band matrix
and its transpose as its only nonzero blocks. Any kernel of the
pullback lies in the radical.

Let $`s_j`$ be the rank of the active Pauli-label span in middle
band $`j`$ in that splitting and set $`K_*=2\sum_j s_j\le 2Ls`$.
There is an $`\epsilon`$-net for the physical averaged bias with
logarithmic cardinality $`\mathrm{poly}(m+K_*+\log(1/\epsilon))`$.
It follows that, with probability at least $`1-e^{-a}`$ over the
same random sign table, simultaneously for all words and values
of $`K_*`$,

```math
 \left|\mathbb E_x F(x)b_V(x,0)\right|
 \le P_*(a+K_*)\,2^{-a/2}
```

for an absolute polynomial $`P_*`$. In particular the margin
vanishes for subexponential middle record count
$`K_*=2^{o(a)}`$ and for sufficiently small exponential rates,
depending on the degree of $`P_*`$. This is a rank condition on the
middle native axes. It gives neither an unrestricted four-band
bound nor a growing lower bound on unrestricted depth.

**A single diagonal parity layer in the proof.** Sections 2–5
also describe any Fourier coefficients
$`\langle B_z|D_{\rm par}|A_0\rangle`$, where
$`|A_0\rangle`$ is a normalized stabilizer vector, $`D_{\rm par}`$
is a diagonal product of T phases on computational parity functions,
and $`|B_z\rangle`$ is the same type of conjugated-Pauli family on
a fixed stabilizer reference. The parity functions need not be
independent. Each weighted parity has the same third-derivative
formula as a physical T coordinate. After choosing cosets, conjugating
a Pauli by the parity phase uses the phase and shear formulas
of Section 4 with $`u_i,w_i`$ replaced by parity evaluations
$`\ell_i(u),\ell_i(w)`$. Its tableau is still affine in the Pauli
parameter. In the realification, the pairwise carry products
are products of two data parities, so they are quadratic in data
with bounded-degree parameter coefficients. The large-overlap
localization and phase reconstruction therefore go through exactly.
Adding a fixed stabilizer register to the $`B_z`$ reference
also preserves the description.

**Coherent middle history.** In each band the active components
of the axes commute, since the spectator components are Z labels.
Choose $`s_j`$ independent commuting active Pauli generators
and their joint eigenprojectors $`P_{j,i}`$. An axis's
eigenvalue sign on a branch is an affine parity of the generator
eigenbits, plus its spectator Z parity. Thus band $`j`$, in
the common frame, acts as

```math
 \sum_i P_{j,i}\otimes D_{j,i},
```

where each $`D_{j,i}`$ is a diagonal spectator parity phase.
The relative global phases are exactly those given by the signed
physical rotations. On each Choi leg coherently write each
band's $`s_j`$ eigenbits into fresh mathematical record bits.
The corresponding isometry is a fixed Clifford change to the
generators' computational coordinates, CNOTs to the record, and
the inverse change. The transposed leg has transposed stabilizer
projectors with the same property. Sequential recording transforms
the original fixed Choi stabilizer into one normalized stabilizer
state $`|\Psi_0\rangle`$ on the data and $`K_*`$ record bits.

All branch-dependent phases can be moved after recording:
they are diagonal on the previous records and spectators,
which later recording acts on trivially. Each resulting exponent
uses one parity of the branch eigenbits and spectator data,
with the literal orientation from the corresponding native gate.
In total this is one product $`D_{\rm par}`$ of diagonal
parity T gates on the mathematical history. Projecting its record
onto all-plus implements the coherent sum of the original middle
branches, with uniform amplitude $`2^{-K_*/2}`$. Therefore the
exact Fourier overlap is

```math
 \beta_z=2^{K_*/2}
 \langle B_z,+^{K_*}|D_{\rm par}|\Psi_0\rangle .
```

The history and projection are proof devices. They are not extra
physical resources or assumptions about native implementations.

The record count can be reduced further. A stabilizer's
computational support is affine, so the record marginal is
supported on an affine subspace. If that subspace has dimension
$`K_{\rm eff}`$, a fixed affine computational change of record
coordinates identifies its support with $`(t,0)`$, where
$`t\in\mathbb F_2^{K_{\rm eff}}`$. The all-plus record bra is
invariant under the permutation. Each inactive record then
contributes a factor $`2^{-1/2}`$, and the remaining phase
depends on affine parities of the active records and spectators.
Thus one can replace $`K_*`$ below by $`K_{\rm eff}`$ in the
net and probability bounds.

Apply the parity-phase version of the net at tolerance
$`\epsilon\,2^{-K_*/2}`$ and scale the approximating coefficients
back. The logarithm in its accuracy argument increases only by
$`K_*/2`$. The approximation remains within $`\epsilon`$ of the
physical bias. In particular it has normalized $`L^2`$ norm
at most $`1+\epsilon`$, since the physical bias is bounded by
one. This avoids an exponentially large factor in the
random-table margin. For uniformity, use the union bound of
Section 5 at each $`K_*`$ with failure probability at most
$`e^{-a-K_*-1}`$ and sum over the nonnegative integer record counts.
Choosing the single polynomial $`P_*`$ large enough gives
the displayed probability and threshold.

### Exact cancellation at large crossing rank

A pair of actual Choi endpoints shows why high overlap by itself
does not bound the absolute cubic ranks of the two decorations.
On two qubits put

```math
 A=(X_1+Z_1X_2)/\sqrt2,\qquad E=T_1T_2,\qquad O=E^\dagger A E.
```

The anticommuting Paulis $`P=X_1`$ and $`Q=Z_1X_2`$ give
$`A=e^{-i\pi(iQP)/8}P e^{i\pi(iQP)/8}`$, a physical one-band
Pauli image after fixed Clifford transport of the rotation axis.
Write $`O=(g_1+g_2+g_3+g_4)/2`$, with

```math
(g_1,g_2,g_3,g_4)=(X_1,-Y_1,Z_1X_2,-Z_1Y_2).
```

These signed Paulis pairwise anticommute. The fixed Clifford
$`C=(i g_1g_4)(g_2+g_3)/\sqrt2`$ swaps $`g_2,g_3`$ and fixes
$`g_1,g_4`$ under conjugation. The second factor swaps the former
pair and negates the latter, while the first factor negates only
the latter. In particular $`COC^\dagger=O`$ literally.
It sends $`Z_1`$ to $`Y_1X_2`$ and $`Z_2`$ to $`-X_1Y_2`$,
so the two copies of the physical E band across the Clifford
have an all-ones cross-commutation matrix of rank one.

The affine Choi support of $`A`$ is
$`y=x+e_1+b(e_1+e_2)`$, where $`x\in\mathbb F_2^2`$
is the input and $`b\in\mathbb F_2`$ chooses its two output
shifts. The third derivative tensor of the
$`E^\dagger\otimes E^{\mathsf T}`$ decoration is nonzero on the directions
$`(e_{x_1},e_{x_1},e_b)`$, from the output coordinate $`y_1`$.
Taking $`k`$ independent copies gives crossing rank $`k`$ and
analytic rank linear in $`k`$ on each side, since the trilinear
biases multiply and a nonzero tensor has bias strictly below one.
The normalized overlap is still one.

These endpoints occur in the actual four-band normal form.
If $`A=F_1P F_1^\dagger`$ with one band $`F_1`$, the native
word $`F_1^\dagger E C E^\dagger F_1`$ preserves $`P`$ by the
displayed commutation. Fixed input and output Cliffords can
identify $`P`$ with any nontrivial input and output Pauli, giving
Fourier coefficient one for the corresponding input character.
On tensor copies each of the four bands acts on disjoint copies
in parallel. This example has a simple resulting bias. Its role
is to rule out a rank-localization argument based only on the
absolute cubic ranks in the actual native representation.

The remaining question is whether relative cancellation between
the two decorations gives a short description at large joint
middle-axis rank. The next section resolves a class in which
each side's Pauli paths at the midpoint have distinct labels.

<a id="four-band-complementary-spans"></a>

## 8. A four-band bound for complementary axis spans

Absorb the Clifford walls as in the
[axis trajectory formulation](AXIS_TRAJECTORY_BIAS.md).
Let $`L_1,L_2,L_3,L_4`$ be the four native band axis spans in
chronological order, with their actual independent axis bases.
Suppose

```math
 L_2\cap(L_1+L_3)=0,\qquad L_3\cap(L_2+L_4)=0,
```

and the two symplectic pairing maps
$`L_1\longrightarrow L_2^*`$ and
$`L_4\longrightarrow L_3^*`$ are injective. For example, these
hypotheses hold for alternating full bases of two complementary
isotropic spans with a perfect pairing, including different bases
on each revisit. Middle crossing, temporal span overlap, and
total axis nullity can all scale with the full axis count in
that example.

For $`0\lt\delta\lt1/2`$ the Fourier coefficient arrays in
Section 1 have a uniform coefficient net at tolerance delta
whose logarithmic cardinality is $`O(n^2\log(n+2))`$, where
$`n=(m+1)(1+\log_2(1/\delta))`$. The bound is independent of
physical width and gate count. In particular, with one clean
bit and $`m=a+1`$, a uniform random sign table has, with
probability at least $`1-e^{-a}`$, absolute signed margin at
most $`O((a+1)^2\sqrt{\log(a+2)}\,2^{-a/2})`$ simultaneously
over this class. An entire initialized-code signed gap implies
at least as large an averaged signed margin. The
[complete-frame measurement](#a-direct-complete-frame-error-floor)
also applies to this class.

### Midpoint marginals penalize active exterior axes

For a fixed input Fourier label z, transport its Pauli through
bands one and two, and transport the output Pauli backward
through bands four and three. Their normalized Pauli coefficient
vectors have inner product $`\beta_z`$. Each half has distinct
Pauli labels for distinct turn words, because $`L_1\cap L_2=0`$
and $`L_3\cap L_4=0`$. A turn at an anticommuting axis uses
one of two terms of equal squared modulus. The squared
coefficient arrays are thus probability distributions on the
Pauli labels, sampled by independent uniform turns at the
active axes of each successive band.

Let $`h_1(z)`$ be the number of active axes in band one.
The second band's activation vector is an affine function
of those $`h_1(z)`$ independent turn bits. Its linear part has
rank $`h_1(z)`$ by the injective pairing. Choose that many
axis coordinates in band two with independent activation bits.
The matching linear coordinate functionals on $`L_2`$ extend
by zero on $`L_1+L_3`$, by the disjointness assumption.
Project both probability distributions onto these functionals.
In the forward half each resulting bit, up to translation by
the input Pauli label, is one exactly when its uniform
activation bit and its independent uniform turn coin are both
one. The projected law is a product of Bernoulli $`1/4`$ bits,
up to those translations. In the backward half the projection
annihilates every turn in band three, so its law is uniform
on an affine subspace from the band four turns.

For probabilities P and Q their affinity is
$`\sum_v\sqrt{P(v)Q(v)}`$. Taking absolute values in the
Pauli inner product and then grouping by the projected
coordinates bounds $`|\beta_z|`$ by the projected affinity.
For a uniform affine law of dimension d on h coordinates,
choose d coordinate pivots giving a bijection from its support
to $`\mathbb F_2^d`$. In the affinity with the biased product
law, sum over these pivots and bound each remaining coordinate
factor by $`\sqrt{3}/2`$. This gives

```math
 \left(\frac{\sqrt{3}}2\right)^{h-d}
 \left(\frac{\sqrt{3}+1}{2\sqrt{2}}\right)^d
 \le \cos(\pi/12)^h.
```

Reversing the argument uses $`L_3\cap(L_2+L_4)=0`$ and the
other injective pairing. Thus a coefficient at least delta
in modulus requires both exterior activation counts at most
$`K=\lceil\log_{1/\cos(\pi/12)}(1/\delta)\rceil`$.

### A sparse description of the retained paths

Only consider addresses z with both exterior counts at most K.
For fixed z and exterior turn words, the two midpoint label
constraints are affine cosets of $`L_2`$ and $`L_3`$.
Their intersection has at most one label because
$`L_2\cap L_3=0`$. Thus there are at most $`2^{2K}`$ matched
path pairs per z. On the affine domain of feasible exterior
turn words the two interior turn words are uniquely determined
affine functions. Exclude any path whose sum of interior
activation counts exceeds $`k=\lceil4K+2\log_2(1/\delta)\rceil`$.
Each excluded contribution has magnitude at most $`2^{-k/2}`$,
so omitting them costs at most delta per coefficient.

Take the affine hull of all valid retained tuples consisting
of z and the two exterior turn words. Its dimension is
$`O((m+1)K)`$. Indeed the hull of the tested z values has an
affine basis of at most m+1 such values. Each exterior mask
is affine in z, so only the union of their at most
$`2(m+1)K`$ active indices on that basis can be used by
exterior turns. If the retained tuple hull is empty, use the
zero array. Otherwise all four turn and activation vectors
are affine on that hull. At each basis point, selected from
the retained tuples, their supports have total size
$`O(K+k)`$. Relabel their union and store the sparse index
lists as in the
[trajectory code](AXIS_TRAJECTORY_BIAS.md#total-axis-nullity-net).
The ordered Pauli phases on both halves combine into one
quadratic exponent modulo four on the hull. The sparse lists,
address values and phase polynomial use $`O(n^2\log(n+2))`$ bits.

The decoder interpolates these data and sums matched paths
with valid turns and both exterior counts at most K. Its
additional valid hull paths, if any, belong to the discarded
large interior-count class, so the omitted-term bound still
applies. Addresses with either exterior count above K have
exact coefficient magnitude below delta by the affinity
bound. Taking real parts, applying Parseval and restricting
the clean coordinate to zero gives normalized $`L^2`$ error
at most $`2^{-a/2}/2`$ at $`\delta=2^{-a-2}`$.
Clip each profile to $`[-1,1]`$ and apply the random-sign
tail bound and union count as in Section 5.

The argument uses a projection that keeps one interior band's
axis coordinates and erases the opposite interior band and
the earlier exterior band. It also uses injectivity of the
two exterior pairings. The large-crossing cancellation example
in Section 7 falls outside these conditions. The unrestricted
asymptotic depth gap remains open.

<a id="four-band-first-pair-rank"></a>

## 9. A four-band bound for small first-pair rank

In a common Clifford frame let r be the rank of the cross-commutation
matrix between the first two native T-band axis spans. The averaged
bias functions of these four-band words have an epsilon net with
logarithmic cardinality $`\mathrm{poly}(m+r+\log(1/\epsilon))`$,
independent of physical width. Here m counts the fixed computational
input coordinates, as in Section 1. At one clean input, almost every
a-bit sign table therefore has signed margin at most
$`P_{\rm ext}(a+r)2^{-a/2}`$ for an absolute polynomial, with failure
probability at most $`e^{-a}`$ simultaneously across all r. A fixed
positive margin requires exponential first-pair rank for typical tables.
This complements the middle-pair bound of Section 7.

### The early Choi vector has a uniform stabilizer expansion

Modulo the radical of their sum, the two isotropic axis spans are
complementary r-dimensional Lagrangians. Lift bases from the original
spans. A fixed Clifford makes their respective bases Z and X on r
active wires, and makes the radical diagonal on spectators. Track each
physical axis sign as an affine offset in its T exponent. At spectator
word y the first two bands become
$`U(y)=H^{\otimes r}D_B(y)H^{\otimes r}D_A(y)`$, with diagonal
modulo-eight exponents $`p_B(e,y),p_A(d,y)`$. Each exponent is a
sum of signed affine parity bits. With $`\omega=e^{i\pi/4}`$,
the active matrix entry is

```math
 U(y)_{a,d}=2^{-r}\sum_e (-1)^{e\cdot(a+d)}
                    \omega^{p_B(e,y)+p_A(d,y)} .
```

In these coordinates let the Fourier input Pauli have X and Z masks
$`u=(u_A,u_Y)`$ and $`v=(v_A,v_Y)`$, and ordered-product scalar
$`i^{s(z)}`$. Its masks are linear in z and its scalar is the
bounded-degree modulo-four Pauli cocycle. Expanding both active
matrix entries of $`UP_zU^\dagger`$ gives its normalized Choi vector
on W physical wires as

```math
 \begin{aligned}
 |UP_zU^\dagger\rangle
   &=2^{-5r/2}\sum_{h=(a,b,d,e,f)} |B(z,h)\rangle,\\
 |B(z,h)\rangle
   &=2^{-(W-r)/2}\sum_y\omega^{q(y;z,h)}
         |a,y+u_Y\rangle\,|b,y\rangle ,
 \end{aligned}
```

where all five history indices are r-bit words and

```math
 \begin{aligned}
 q(y;z,h)
   ={}&p_B(e,y+u_Y)-p_B(f,y)\\
      &+p_A(d+u_A,y+u_Y)-p_A(d,y)\\
      &+4\big[e\cdot(a+d+u_A)+f\cdot(b+d)
                   +v_A\cdot d+v_Y\cdot y\big]+2s(z).
 \end{aligned}
```

The coefficient uses the Choi normalization $`2^{-W/2}`$,
two matrix-entry factors $`2^{-r}`$, and the norm
$`2^{(W-r)/2}`$ of each unnormalized spectator sum.

In either band's paired exponent, each spectator parity
$`\ell_Y(y)`$ occurs with the same signed weight. Its difference is

```math
 (\ell_Y(y)\oplus\alpha)-(\ell_Y(y)\oplus\beta)
                 =(\alpha-\beta)(1-2\ell_Y(y)),
```

where alpha and beta are affine bits of $`(z,h)`$. This literal
integer identity removes every degree-three spectator term. After
extracting $`\omega^{q(0;z,h)}`$, the relative phase is a single
diagonal Clifford wall on fixed spectator parities, with affine
S or S-dagger occurrence bits and a linear Z mask. The scalar has
bounded-degree modulo-eight digits by the phase arithmetic in
Section 4, including the active quadratic Hadamard signs. A fixed
spectator Bell preparation, this diagonal wall, and affine Pauli
shifts on the active labels and first spectator leg prepare every
$`|B(z,h)\rangle`$ uniformly. Realifying the parity S gates uses
the same pairwise carries and fixed rebit Hadamards as in Sections
4 and 7. Those carries are quadratic in spectator data and have
bounded-degree parameter coefficients. The histories and Bell
registers are mathematical descriptions of overlaps.

### One description for all histories

Absorb the fourth band into the measured Pauli. Its Choi vector is
stabilizer, and the third band contributes one fixed T decoration
after a fixed change of Clifford frame. Denote the resulting
fixed decorated vector by $`|\psi\rangle`$. The Fourier coefficients
are $`2^{-5r/2}\sum_h\langle B(z,h)|\psi\rangle`$.
Sections 2–4 apply to the bounded-round Clifford preparation
family above. A significant overlap gives one common cubic
restriction on the support of the fixed decorated vector. After
the coset split, lift by all monomials up to the absolute degree
bound of the relative symplectic entries. The active support,
signed Pauli descriptions, and literal scalar reconstruction then
use a polynomial number of bits across all $`(z,h)`$.

The history coefficients have one-norm $`2^{5r/2}`$. Set the
term tolerance proportional to $`\epsilon 2^{-(m+5r)/2}`$ before
the coset split. The summed coefficient error and Parseval give
normalized $`L^2`$ error at most epsilon for the physical bias.
The increased log inverse tolerance is $`O(m+r+\log(1/\epsilon))`$.
The decoded bias has $`L^2`$ norm at most $`1+\epsilon`$.
Thus the random-sign comparison in Section 5 pays no history
factor in the final margin. Sum the failure probabilities
$`e^{-a-r-1}`$ over r to give the simultaneous bound.

### The actual adjoint also screens the last pair for complete corners

For the [complete deepest-corner frames](#a-direct-complete-frame-error-floor),
the actual adjoint of an initialized-isometry approximation has the
same error against the ideal adjoint, by the
[initialized-isometry inverse identity](BOOLEAN_QUERY_FRAME_REDUCTION.md#approximation-leakage-and-the-actual-inverse).
It reverses native T bands and preserves their cross-commutation
ranks under the corresponding change of Clifford frame. On target
input zero, the ideal adjoint has the same computational decision
outcome f(x) as the ideal forward frame. Apply the same error floor
and first-pair bound to the reverse to screen the original last pair.
At fixed error below $`\sqrt{2-\sqrt2}`$ and polynomial clean count,
these decision tests impose exponential cross-commutation rank at all
three adjacent pairs of a hypothetical four-band complete corner.
The [complete-isometry test below](#complete-return-four-band-floor)
excludes even that case when the clean count stays a fixed fraction below
the address width.
This inference uses the complete initialized action of the frame and
its adjoint. An arbitrary moving raw decision can rotate its sign
subspaces into different addresses under reversal.

<a id="complete-return-four-band-floor"></a>

## 10. Complete return excludes four bands

Let $`a`$ be the number of deepest address bits and let $`f`$ be a
uniform Boolean table. Write $`W_f`$ for the complete deepest-corner
frame from Section 5, tensored with the dirty identity, and let
$`J`$ append $`k`$ initialized zeros. Fix
$`\eta\lt\sqrt{2-\sqrt2}`$ and $`\delta\gt0`$. For all sufficiently
large $`a`$, with probability at least $`1-e^{-a}`$, every native word
satisfying $`\|VJ-e^{i\phi}JW_f\|\le\eta`$ with
$`k\le(1-\delta)a`$ has at least five physical T bands. This
covers arbitrary Clifford interlayers, moving addresses, unrestricted
dirty width, arbitrary dirty inputs with references, and the actual
inverse. Zero leaf phases give the same bound for the prescribed
complex magnitude-completion family. The error threshold is sharp
by the single target Clifford in Section 5.

### Three-band traces with Pauli probes

Use the target $`Y`$ eigenbasis, with label $`b=0`$ for its positive
eigenvector. The eigenvalue of $`XZ`$ is $`-i(-1)^b`$. A fixed
target Clifford implements this choice of coordinates. For signed
Hermitian Paulis $`P,Q`$ and a word $`B`$ with at most three
physical T bands, test the complex profile

```math
h(x)=\tau_{\mathrm{dirty}}\langle x,b,0^k|P^\dagger B Q B^\dagger|x,b,0^k\rangle .
```

Here $`\tau_{\mathrm{dirty}}`$ is normalized trace, and $`b`$ is
fixed. Both quadratures of these profiles have an
$`\epsilon`$-net in normalized $`L^2`$ with logarithmic cardinality
$`\mathrm{poly}(a+k+1+\log(1/\epsilon))`$. To see this, first vary
all $`a+k+1`$ computational coordinates and let $`L_z`$ be their
Pauli $`Z`$ character at Walsh label $`z`$. The ambient Fourier
coefficient is $`\tau(QB^\dagger(L_zP^\dagger)B)`$. Apply Sections 2–5
to the actual word $`B^\dagger`$, with fixed output Pauli $`Q`$ and
affine scalar Pauli probe $`L_zP^\dagger`$. In the notation of
Section 2 the Choi overlap is $`\langle B_z^\dagger|\psi\rangle`$.
The adjoint retains affine Pauli labels and bounded-degree
modulo-four phase digits, so the tableau and scalar reconstruction
apply to this ket, including when the trace is complex.
The bounds on individual overlaps use modulus and thus apply to both
quadratures. Restrict to the desired slice after approximating at
ambient accuracy $`\epsilon 2^{-(k+1)/2}`$.

Set $`r_x=(-1)^{f(x)}`$ and define the two challenges

```math
 v_{u,+}(x)=(r_x+r_{x+u})/2,\qquad
 v_{u,-}(x)=(r_{x+u}-r_x)/2 .
```

For both quadratures $`g`$ of every profile above, both $`b`$, and
every affine hyperplane $`H`$ of addresses or the full address
space, a common event gives

```math
 \big|\mathbb E_x\,1_H(x)v_{u,\pm}(x)g(x)\big|
 \le p(a+k+1)2^{-a/2}
```

for an absolute polynomial $`p`$. In fact use
$`\epsilon=2^{-a/2}`$ in the net and clip each quadrature to
$`[-1,1]`$. Multiplying an approximant by $`1_H`$ and
symmetrizing or antisymmetrizing under $`x\mapsto x+u`$ gives a
Rademacher test with normalized $`L^2`$ norm at most one. Apply
Section 5 and union over the net, all $`u,H,b`$, and both signs
and quadratures. The polynomial can absorb a union over
$`0\le k\le a`$, with failure probability at most $`e^{-a}/2`$.
Also on this event, after allocating the remaining failure budget,
each challenge used below is nonzero on at least one third of
every such $`H`$. For $`u\ne0`$ these counts use independent
address pairs, each with multiplicity at most two in $`H`$.
For $`u=0`$ the plus challenge has full support and the minus
challenge is unused.

### Complete action gives a correlation for every translating probe

A Pauli $`P`$ that commutes with every clean $`Z`$ preserves the
initialized code. Let $`u`$ be its address translation and
$`v`$ its target eigenbit translation. On target eigenbit $`b=0`$,
the ideal conjugated Pauli $`W_f(J^\dagger P J)W_f^\dagger`$
agrees on each column with $`J^\dagger P J`$ times phase

```math
 \xi_{u,v,0}(x)=
 \begin{cases}
 i^{f(x)-f(x+u)},&v=0,\\
 i^{f(x)+f(x+u)},&v=1.
 \end{cases}
```

At $`b=1`$ the phase is conjugated. Thus for $`v=1`$ its
real part is $`v_{u,+}`$, and for $`v=0,u\ne0`$ its imaginary
part on $`b=0`$ is $`v_{u,-}`$. On the respective challenge
supports the ideal phase is a signed real or imaginary unit.

The [actual initialized inverse identity](BOOLEAN_QUERY_FRAME_REDUCTION.md#approximation-leakage-and-the-actual-inverse)
gives $`\|V^\dagger J-e^{-i\phi} JW_f^\dagger\|\le\eta`$.
Apply it to a lifted column with any unit dirty input, then apply
$`P`$ and the forward approximation on the initialized code. The
two comparisons are between unit vectors at distance at most
$`\eta`$ each. Their real spherical angles are at most
$`2\arcsin(\eta/2)`$ each, so the real overlap of the actual and
ideal conjugated vectors is at least

```math
 \gamma(\eta)=2(1-\eta^2/2)^2-1\gt0 .
```

The trace with comparison Pauli $`P^\dagger`$ and conjugation
$`VPV^\dagger`$ therefore has the correct signed quadrature at
least $`\gamma(\eta)`$ on each challenge column. This remains
true after dirty averaging. Any nonzero logical translation
gives correlation at least $`\gamma(\eta)/3`$ in the full
address space.

### Sparse physical activation and the four-band exclusion

Work in binary symplectic Pauli labels. Let $`C_Z`$ be the span
of initialized $`Z`$ labels and $`L`$ the span of logical
address and target eigenbit $`Z`$ labels. Pull Clifford walls
to the output and write $`V=B E_1`$, where $`E_1`$ is the
first physical band in the common input frame. Its independent
Pauli axes span an isotropic space $`A_1`$. If $`B`$ has at
most three bands, every Pauli in
$`A_1^\perp\cap C_Z^\perp`$ commutes with $`E_1`$ and has
a three-band trace. The correlation lower bound forces it to
have zero logical translation. Taking symplectic orthogonals
gives $`L\subseteq A_1+C_Z`$.

For the exclusion, start with a target-flipping Pauli and vary
every address and dirty Pauli label, the target eigenbit
$`Z`$, and all clean $`Z`$ labels. These variations preserve
clean and still flip the target eigenbit. In the independent
axis coordinates their commutation syndromes have an annihilator
identified with $`A_1\cap\mathrm{span}(Z_b,C_Z)`$, of dimension
at most $`k+1`$. Any coset of a binary space of codimension
$`j`$ has a representative of weight at most $`j`$, by
choosing independent pivot columns for its constraints. Thus
some target-flipping Pauli $`P`$ anticommutes with at most
$`k+1`$ individual physical axes. Conjugating through an
anticommuting T rotation splits a Pauli into two signed terms
with coefficient magnitudes $`2^{-1/2}`$. The Pauli expansion
of $`E_1P E_1^\dagger`$ therefore has coefficient
$`\ell^1`$ norm at most $`2^{(k+1)/2}`$. Combining its
terms with the three-band trace bound gives correlation at most
$`2^{(k+1)/2}p(a+k+1)2^{-a/2}`$, which tends to zero under
the stated clean bound. This contradicts complete return.

<a id="complete-return-pauli-spreading"></a>

### Any longer word must spread the translating probes at both ends

For $`d\ge3`$, let $`F`$ be the first $`d-3`$ bands in
the common input frame of a $`d`$-band word. The same trace
comparison forces every clean-preserving Pauli $`P`$ with
nonzero logical translation to have

```math
 \|FPF^\dagger\|_{\mathrm{Pauli},1}
 \ge\frac{\gamma(\eta)2^{a/2}}{3p(a+k+1)} .
```

Here the norm sums absolute coefficients in the normalized
orthogonal Pauli basis. In particular the sum of the prefix
axis spans with $`C_Z`$ contains $`L`$, by applying the
bound to Paulis commuting with that sum. Also start with the
label of $`P`$ and collect axes that anticommute with some
term entering their own band. All final labels lie in the
translate of the span of these axes. If this span has binary
dimension $`j`$, then the Pauli coefficient norm above
is at most $`2^{j/2}`$ by normalized Hilbert–Schmidt norm one
and Cauchy–Schwarz. Thus $`j\ge a-O(\log(a+k+2))`$ at fixed
error. Applying the same calculation to the actual adjoint
gives the corresponding suffix constraints. The ideal adjoint
conjugates the literal phases and preserves the challenge
supports. These spreading bounds allow a single wide band to
activate many axes.

<a id="complete-return-five-band-intersection"></a>

## 11. A five-band axis-intersection requirement

Use the same typical deepest corners and fixed error as in Section 10.
With zero or one clean input, for every actual five-band decomposition
the first two axis spans in the common Clifford frame must intersect
in dimension at least $`a-O(\log(a+2))`$. The last two obey the
same bound. These are binary dimensions of intersections, distinct
from the ranks of cross-commutation matrices. The result allows
arbitrary dirty width and still leaves five-band words with large
intersections open.

Write $`V=B E_2 E_1`$ with $`B`$ using the remaining three
bands. Let $`A_1,A_2`$ be the two early axis spans,
$`r=\dim(A_1\cap A_2)`$, and suppose
$`r\le a-C\log(a+2)`$ for a sufficiently large constant
$`C`$ depending on the error. First strengthen
$`L\subseteq A_1+A_2+C_Z`$ from Section 10 to
$`L\subseteq A_1+C_Z`$. To do so use probes commuting with
$`A_1+C_Z`$. Their $`E_2`$ syndrome space has codimension
$`\dim(A_2\cap(A_1+C_Z))\le r+k`$ in the independent physical
axis coordinates. Since $`L\subseteq A_1+A_2+C_Z`$, any
logical translation of these probes is a linear function of their
$`E_2`$ syndrome. If that function were nonzero, it would
be nonzero on an echelon basis word of weight at most $`r+k+1`$.
Such a probe is unchanged by $`E_1`$, and its coefficient norm
after $`E_2`$ is at most $`2^{(r+k+1)/2}`$, contradicting
the prefix bound. The stronger inclusion follows. In particular
$`A_1\cap C_Z^\perp`$ is logical-diagonal by isotropy.

If $`k=1`$, let $`K_*`$ be the binary sum of all actual
$`E_1`$ axes anticommuting with the clean $`Z_c`$.
For $`k=0`$ set $`K_*=0`$. Choose one logical coordinate
label $`T_0`$ among the address and target eigenbit $`Z`$ labels
outside $`C_Z+\mathrm{span}(K_*)`$. Choose a Pauli $`P`$
that anticommutes with $`T_0`$ and commutes with $`C_Z`$
and $`K_*`$. The combined activity syndrome on the
concatenated physical axis lists can be chosen of weight at
most $`r+k+2`$. Indeed the linear map from the concatenated
axis coefficients to physical Pauli labels has kernel dimension
$`r`$. The annihilator of the allowed syndrome variations is
the preimage of $`C_Z+\mathrm{span}(T_0,K_*)`$, of dimension at
most $`r+k+2`$. The same pivot argument as in Section 10
then supplies the claimed weight. Write its two activity weights
as $`\ell_1,\ell_2`$.

### Removing the first band while controlling clean overlap

Put $`Q=E_1^\dagger P E_1`$, $`Q_c=J^\dagger Q J`$,
$`\Pi=JJ^\dagger`$, and $`V_\mathrm{tail}=V E_1^\dagger`$.
The exact trace identity uses the original physical inverse,

```math
 J^\dagger V_\mathrm{tail}P V_\mathrm{tail}^\dagger J
 =(V^\dagger J)^\dagger Q(V^\dagger J).
```

Only clean-preserving Pauli terms of $`Q`$ survive in $`Q_c`$.
They differ in label from $`P`$ by vectors in
$`A_1\cap C_Z^\perp`$ and hence share the same nontrivial
logical translation $`(u,v)`$. Thus on each $`(x,b)`$ column
$`W_f Q_c W_f^\dagger`$ agrees with $`Q_c`$ times the phase
$`\xi_{u,v,b}(x)`$ of Section 10.

When $`k=0`$, $`Q_c^\dagger Q_c=I`$. For $`k=1`$
write the conjugation as $`Q=P G\exp(iH_1)`$, where $`G`$
is the product of the activated clean-preserving Pauli quarter
turns and $`H_1`$ is the signed $`\pi/4`$ sum of the
activated clean-anticommuting axes. All these axes commute.
The compression of $`\exp(iH_1)`$ is the restriction of
$`\cos(H_1)`$ because $`Z_c H_1 Z_c=-H_1`$. Both
$`P`$ and $`G`$ preserve clean, so
$`Q_c^\dagger Q_c=J^\dagger\cos^2(H_1)J`$. The number of axes
in $`H_1`$ is even because $`P`$ commutes with $`K_*`$.
If there are none the square is $`I`$. Otherwise
$`\cos(2H_1)`$ is a signed Hermitian Pauli $`K`$, the
product of these axes up to sign, with label neither zero
(by axis independence) nor $`Z_c`$ (which could not belong
to $`A_1`$ in the presence of any such axis). This product
preserves clean and is logical-diagonal. Thus the square is
the projection $`(I+J^\dagger KJ)/2`$.

For each $`(x,b)`$ choose a dirty basis diagonalizing this
projection on that block. If $`K`$ has nonzero dirty support,
half of the dirty columns have projected square one at every
$`(x,b)`$. If $`K`$ has no dirty support, the active
columns are cut out by one nontrivial affine equation in
$`(x,b)`$. Choose $`b`$ so its slice is all addresses or an
affine hyperplane $`H`$. In all cases at least half the dirty
columns have projected square one at each $`x\in H`$.

The tested profile is

```math
 h(x)=\tau_{\mathrm{dirty}}
 \langle x,b,0^k|Q^\dagger\Pi V_\mathrm{tail}
     P V_\mathrm{tail}^\dagger|x,b,0^k\rangle .
```

A dirty column with square zero contributes zero. For an
initialized logical-dirty unit column $`w`$ of square one,
$`Q_cw`$ is a unit vector. Pair the actual vectors
$`V^\dagger J Q_cw`$ and $`Q V^\dagger Jw`$.
Each is within $`\eta`$ of the version with actual
$`V^\dagger J`$ replaced by $`e^{-i\phi}JW_f^\dagger`$.
The ideal inner product is exactly $`\xi_{u,v,b}(x)`$, so the
two-angle comparison gives the correct signed quadrature at
least $`\gamma(\eta)`$ on the challenge support.
After dirty averaging, the signed quadrature there is at least
$`\gamma(\eta)/2`$ on $`H`$. The masked correlation is
therefore at least $`\gamma(\eta)/12`$.

Finally expand $`Q^\dagger`$ and $`E_2 P E_2^\dagger`$
in Paulis, at coefficient norm at most
$`2^{\ell_1/2}`$ and $`2^{\ell_2/2}`$ respectively.
The clean projection has Pauli coefficient norm one. Every
remaining term has a three-band trace. Scalar Pauli products
can rotate its real and imaginary quadratures by a power of
$`i`$, already covered by the same event. The correlation
upper bound is consequently
$`2^{(r+k+2)/2}p(a+k+1)2^{-a/2}`$, which contradicts the
positive lower bound for the assumed $`r`$. Apply the
same argument to the actual adjoint to obtain the last-pair
statement. All projections in this argument are trace probes
on the initialized isometry, and incur no physical preparation
or postselection cost.

<a id="complete-return-two-table-frame"></a>

## 12. Two random Hopf levels force a nonisotropic logical span

For every fixed $`\eta\lt\sqrt{2-\sqrt2}`$ and every fixed
$`M\gt0`$, a different typical family of complete frames
requires at least five physical T bands even with up to
$`a^M`$ clean inputs. Here $`a`$ is again the number of
deepest address bits. For an actual five-band word with
$`k`$ initialized inputs, the forward and reverse first-pair
axis intersections each have dimension at least
$`a-k-O(\log(a+k+2))`$. These statements hold with probability
at least $`1-e^{-a}`$ for sufficiently large $`a`$, over
the two independent tables below. The five-band test needs only
one expansion through an exposed physical band.

### The literal two-table frame

Use computational logical coordinates $`(x,y,b)`$, where
$`x`$ has $`a-1`$ bits and $`y,b`$ are single bits.
Take independent uniform Boolean tables $`g(x)`$ and
$`f(x,y)`$. The penultimate angle is $`\pi g(x)/2`$,
with the actual complete suffix predicate $`b=0`$.
It produces a signed permutation $`F_g`$ that on that
sector applies $`K_y^{g(x)}`$ for $`K=XZ`$. The deepest
angle is $`\pi f(x,y)`$. As $`R_y(\pi)=-I`$, this
produces the same scalar $`(-1)^{f(x,y)}`$ on both
$`b`$ columns. Call this diagonal $`G_f`$. All earlier
angles are zero, so the complete product is $`W=G_fF_g`$.
For the complex frame use magnitude operator $`F_g`$ and
leaf phases $`\pi f(x,y)`$ on both $`b`$ values.

Restrict to clean-preserving Pauli probes $`P`$ commuting with
$`Z_b`$, tested on $`b=0`$. Write their address
translation as $`(u,t)`$ and their binary $`Z_y`$
label as $`z_y(P)`$. Suppose $`(u,t)\ne0`$ and
choose the largest of the four masks

```math
 H_{u,p,q}=\{(x,y):g(x)=p,\ g(x+u)=q\}.
```

Its density is at least one quarter. On this mask, conjugation
through $`F_g`$ agrees with the fixed scalar Pauli
$`Q_+=K_y^q P K_y^{-p}`$. It translates the address by
$`d=(u,t+p+q)`$, which is nonzero because $`p=q`$
on a nonempty mask if $`u=0`$. Putting $`z=(x,y)`$
and $`S_{d,s}(z)=(-1)^{f(z+s)+f(z+s+d)}`$, the ideal forward
conjugation on the mask equals $`S_{d,0}(z)Q_+`$.
For the reverse, use $`Q_-=K_y^{-q}P K_y^p`$. The
two deepest signs are evaluated at $`(x,y+p)`$ and
$`(x+u,y+p+t)`$, so ideal reverse conjugation equals
$`S_{(u,t),(0,p)}(z)Q_-`$. In either direction the phase
relative to the comparison Pauli is a real sign. Its fixed
scalar Pauli factor can be included by choosing a quadrature.

If instead $`(u,t)=0`$ and $`z_y(P)=1`$, both deepest
signs cancel. On $`b=0`$, conjugation of $`Z_y`$
by $`K_y`$ or its inverse changes sign. Thus the ideal test
compares with $`P`$ itself and has the sign
$`(-1)^{g(x)}`$ in either direction, over the full
address space.

### Common sign tests and the axis bounds

For fixed $`g`$, use Section 10's three-band trace nets
on this computational basis. Test their net points against
$`1_{H_{u,p,q}}S_{d,s}`$ for all masks and all nonzero
$`d`$, with $`s`$ a fixed translation of the bit
$`y`$. For each $`d`$, the products of $`f`$ signs
on the unordered pairs $`\{z+s,z+s+d\}`$ are independent
Rademachers. Each coefficient in the masked correlation sums
at most two bounded net values and divides by $`2^a`$.
The sum of squared coefficients is at most $`2^{1-a}`$.
The separate sign test $`(-1)^{g(x)}`$, repeated on the
two $`y`$ values, has the same bound without conditioning.
Hoeffding, the union bound and the net approximation give a
common correlation upper bound $`p(a+k+1)2^{-a/2}`$ for
both tests, uniformly over the masks, all Pauli probes and
three-band words, and every $`k\le a^M`$. The polynomial
absorbs union over the finite mask labels and clean counts.
Allocate half the failure budget to each table. In the first
test conditioning on $`g`$ fixes all masks before sampling
$`f`$.

The complete forward and actual inverse errors give the same
two-angle comparison as Section 10, now using the appropriate
comparison Pauli above. It holds pointwise on arbitrary dirty
columns. Thus every probe with $`(u,t,z_y(P))\ne0`$
gives correlation at least $`\gamma(\eta)/4`$ in one of
these tests. This holds for the approximation to $`W`$ and
the actual adjoint approximation to $`W^\dagger`$.

Put $`C_*=C_Z+\mathrm{span}(Z_b)`$ and

```math
 L_* =\mathrm{span}(Z_{x_1},\ldots,Z_{x_{a-1}},Z_y,X_y).
```

The detectable probes are precisely those in
$`C_*^\perp`$ with nonzero pairing against $`L_*`$.
As in the spreading argument, for a $`d`$-band word,
the required correlation forces

```math
 L_*\subseteq A_1+\cdots+A_{d-3}+C_*.
```

Indeed a Pauli in the orthogonal of the right side commutes
through the whole prefix, leaving a three-band trace. With at
most four bands the inclusion is impossible. The logical
$`X_y,Z_y`$ pair commutes with $`C_*`$ and still
anticommutes when either element is changed by an element of
$`C_*`$. Two lifts into the single isotropic span
$`A_1`$ would commute.

If the word has five bands, let
$`j=\dim(A_2\cap(A_1+C_*))`$. For probes commuting with
$`A_1+C_*`$ their detected labels
$`(u,t,z_y(P))`$ are linear functions of their $`E_2`$
syndrome by the prefix inclusion. At least one such function
is nonzero, since $`L_*\subseteq A_1+C_*`$ is impossible.
The $`E_2`$ syndrome space has codimension $`j`$ and
thus has an echelon basis with weight at most $`j+1`$
per vector. Some basis vector has a detectable label.
Its Pauli commutes through $`E_1`$ and the
$`E_2`$ expansion has coefficient norm at most
$`2^{(j+1)/2}`$. The three-band correlation upper bound
therefore forces $`j\ge a-O(\log(a+k+2))`$. Since
$`j\le\dim(A_1\cap A_2)+k+1`$, this gives the
claimed first-pair bound, and applying it to the actual adjoint
gives the last-pair bound. Large intersections at five bands
remain possible under these conditions.

<a id="zero-clean-six-band-floor"></a>

## 13. A rare three-cycle requires six bands with zero clean inputs

For any fixed $`\eta\lt\sin(\pi/24)/2`$ and all sufficiently
large address sizes, an explicit complete frame needs at least six
physical T bands at full operator error $`\eta`$ in the
zero-clean model. Arbitrary dirty width, gate counts, Clifford
interlayers, physical address motion and common approximation phases
are allowed. The bound also holds for the complex magnitude completion
with zero leaf phases.

Fix one marked prefix of $`p`$ bits and call the last two
logical bits $`y,b`$. Outside that prefix, set every angle to
zero. At the marked prefix set the penultimate angle on $`y`$
to $`\pi/2`$, with its actual suffix predicate $`b=0`$.
Set the deepest angle on $`b`$ to $`\pi/2`$ at
the marked prefix extended by $`y=0`$. All other angles
vanish. The ordered frame $`W`$ fixes the unmarked blocks.
On the marked $`(y,b)`$ block its literal signed action is

```math
 |00\rangle\mapsto |10\rangle,qquad
 |10\rangle\mapsto-|01\rangle,qquad
 |01\rangle\mapsto-|00\rangle,qquad
 |11\rangle\mapsto |11\rangle.
```

The sign product on the cycle is positive, so the spectrum includes
$`1,e^{2\pi i/3},e^{-2\pi i/3}`$. Write
$`\|\cdot\|_{2,\tau}`$ for the normalized Hilbert–Schmidt norm.
The marked block has trace one on four dimensions, and tensoring
the dirty identity leaves the normalized norm unchanged. Hence

```math
 \|W-I\|_{2,\tau}=\sqrt{3/2}\,2^{-p/2}.
```

Suppose $`V`$ has at most five physical bands and
$`\|V-cW\|\le\eta`$ for a common phase $`|c|=1`$.
Then $`\|V-cI\|_{2,\tau}\le\eta+\sqrt{3/2}\,2^{-p/2}`$.
Let $`A`$ be the isotropic span of its first physical
band's axes in the input frame. For any Hermitian Pauli $`P`$
in $`A^\perp`$, the conjugation $`VPV^\dagger`$ uses at
most four bands because $`P`$ commutes with the first one.
By the [four-band mean Pauli gap](FOUR_LAYER_SCALAR_GAP.md#mean-error-exceptional-rank-and-complete-operators),
either $`VPV^\dagger=P`$ or

```math
 \|VPV^\dagger-P\|_{2,\tau}
 \ge 2\sqrt{(1-\cos(\pi/12))/8}=\sin(\pi/24).
```

On the other hand, inserting and subtracting $`cI`$ in
the conjugation gives distance at most $`2\|V-cI\|_{2,\tau}`$.
For sufficiently large $`p`$ this is strictly below the gap.
Thus $`V`$ commutes with every Pauli in $`A^\perp`$.
Expanding in the Pauli basis shows that its support lies in
$`(A^\perp)^\perp=A`$. A Clifford diagonalizes this
isotropic algebra. In that basis the diagonal entries remain
in the native ring $`\mathbb Z[i,1/\sqrt2]`$ and have
unit modulus. By the
[unit-circle classification](FOUR_LAYER_OPERATOR_RIGIDITY.md#4-discrete-phases-in-the-common-commutant)
they are eighth roots of unity.

For an eigenvector of $`cW`$ with eigenvalue $`\lambda`$,
the norm of $`V-\lambda I`$ on that vector is at most
$`\eta`$. Normality of $`V`$ puts an eighth root
within $`\eta`$ of $`\lambda`$. In particular some
eighth roots $`\zeta_0,\zeta_1`$ approximate
$`c,ce^{2\pi i/3}`$, respectively. Their ratio satisfies

```math
 |\overline{\zeta_0}\zeta_1-e^{2\pi i/3}|\le2\eta.
```

The nearest eighth root to $`e^{2\pi i/3}`$ is at
distance $`2\sin(\pi/24)`$, contradicting the assumed
error. The comparison to the identity uses the full physical
operator. Initialization of even one extra clean qubit would
leave its other input sector outside that norm guarantee.

<a id="five-band-third-level-clean-boundary"></a>

## 14. The rare cycle constrains one-clean exterior images

At sufficiently small fixed complete error and large p, a five-band
one-clean compiler for Section 13 must have **both** exterior
two-band images of clean Z outside the third Clifford hierarchy
level $`\mathcal C_3`$. These are the forward image through the
first two bands and the backward image through the last two,
using the actual Clifford walls. A member of $`\mathcal C_3`$
conjugates every Pauli to a Clifford. Arbitrary dirty width and
physical address motion are allowed.

### A third-level exterior image forces exact clean return

Let J append one clean zero. An approximation $`\|VJ-cJW\|\le\eta`$
with $`|c|=1`$ gives leakage from zero at most eta. Unitarity
gives the same norm for the other off-diagonal block: their
squared singular spectra are the deficiencies of the occupied
square block in its two multiplication orders. Thus
$`\|VZV^\dagger-Z\|\le2\eta`$ on the entire physical space.

Use the actual reverse if necessary to place the known third-level
image L at the forward two-band prefix. At the cut after the
third band E compare $`A=E L E^\dagger`$ and the backward
two-band clean image $`B=E_B K_B E_B^\dagger`$, where E and
$`E_B`$ are conjugated physical bands and $`K_B`$ is a
Hermitian Clifford. Their operator distance is at most $`2\eta`$.
For every Pauli R commuting with the axes of E, the conjugation
$`ARA=E(LRL)E^\dagger`$ has one band decoration of a Hermitian
Clifford. So does BRB, since $`E_B^\dagger R E_B`$ is Clifford.
The distance of these probe images in normalized Hilbert–Schmidt
norm is at most $`4\eta`$. The
[four-band mean gap](FOUR_LAYER_SCALAR_GAP.md#mean-error-exceptional-rank-and-complete-operators)
$`g_4=\sin(\pi/24)`$ forces equality below that threshold.

Thus BA commutes with the entire Pauli centralizer of the axis
span of E. Its Pauli support lies in that isotropic span, so a
Clifford diagonalizes BA over the native ring. The
[unit-circle classification](FOUR_LAYER_OPERATOR_RIGIDITY.md#4-discrete-phases-in-the-common-commutant)
gives eighth-root eigenvalues. Since $`\|BA-I\|\le2\eta`$,
at sufficiently small eta it is identity. Consequently
$`[V,Z]=0`$ exactly. Reorient by actual reversal as needed
to put the known third-level image S on the backward suffix.

### A common phase frame for S

Write $`S=E_B J_B E_B^\dagger`$ with $`J_B`$ Hermitian
Clifford. For each actual axis $`P_i`$ of $`E_B`$, there
is a Pauli R activating only that axis, by independence.
In $`S R S`$, factor out the Clifford
$`E_B(J_B R J_B^\dagger)E_B^\dagger`$.
The remaining factor has form $`(I\pm iK)/\sqrt2`$,
where $`K=E_B(J_B P_i J_B^\dagger)E_B^\dagger`$ is traceless.
Since S is third-level this factor is Clifford. A Clifford
has flat Pauli coefficient magnitudes on its support.
The identity coefficient $`1/\sqrt2`$ forces support
of size two, hence K is Pauli. Conjugating a Pauli by
an independent-axis band gives $`2^s`$ distinct Pauli
terms if s axes activate, so $`J_B P_i J_B^\dagger`$
commutes with every axis of $`E_B`$.

The union of these two axis lists is isotropic. Thus
$`S=D_S J_B`$, with $`D_S=E_B(J_BE_B^\dagger J_B^\dagger)`$
a diagonal weighted-cubic phase in a fixed Clifford frame
shared with $`E_B`$. Conversely, that commutation condition
implies S is third-level. If it fails, some K in the test
above has all nonidentity coefficients of modulus at most
$`1/\sqrt2`$. In $`(I\pm iK)/\sqrt2`$, the identity
coefficient is $`1/\sqrt2`$ and all the others at most 1/2.
A flat Clifford Pauli vector either differs at the identity
by at least $`1/\sqrt2-1/2`$, or has support of size two
containing the identity and differs by at least that amount
at its other entry. The Pauli-image telescope shows that
a non-third-level two-band image is at normalized
Hilbert–Schmidt distance at least
$`(1/\sqrt2-1/2)/2`$ from the third level, even up to phase.

### Occupied Pauli comparisons with this phase frame

At a two-plus-two cut, write the forward and backward
Hermitian Pauli images $`A=E_A K_A E_A^\dagger`$ and
$`B=E_B K_B E_B^\dagger`$. Let $`Q=(I+S)/2`$ and
assume A and B commute with Q. This commutation will follow
from exact terminal return in the application. Put
$`\varepsilon=\|(A-B)Q\|_{2,\tau}`$ with the full physical
trace normalization.

For a Hermitian Pauli R activating at most five axes of
$`E_A`$, its $`s`$ activated axes give $`2^s`$ terms
in $`E_A^\dagger R E_A`$. Hence ARA is a sum of
$`2^s`$ Hermitian Cliffords with coefficient modulus
$`2^{-s/2}`$. For each such Clifford G, the terms
G, SG, GS and SGS of its compression are flat affine
Choi states with weighted-cubic phases in the common
diagonal frame. Indeed $`SG=D_S J_B G`$, $`GS=(SG)^\dagger`$,
and $`SGS=D_S J_B G J_B^\dagger D_S^\dagger`$.
The four terms from compressing BRB remain $`E_B`$
decorations of Cliffords because
$`Q=E_B(I+J_B)E_B^\dagger/2`$.
After absorbing literal eighth-root phases and duplicating
terms using $`\sqrt2=\omega+\omega^{-1}`$ with
$`\omega=e^{i\pi/4}`$, all coefficients have modulus
$`2^{-(s+4)/2}`$, with at most 160 terms and $`s+4\le9`$.

The [finite-boundary gap](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
at zero active doubled coordinates gives an absolute
$`g_{\rm pal}\gt0`$ for a nonzero such palette. The telescope
$`\|Q(ARA-BRB)Q\|_{2,\tau}\le2\varepsilon`$ freezes the
compressed probes below half that gap. Zero-activation
Paulis span the full operator blocks on the axis fibers z.
Their compression includes all rank-one maps on each
occupied fiber image in Q, so BA is scalar there.
A fixed activation mask e similarly spans every block
$`z\to z+e`$. Nonzero compression between two occupied
fibers identifies their scalars when $`|e|\le5`$.

To see connectivity, in the $`E_A`$ axis coordinates
an unoccupied z satisfies $`S|z,v\rangle=-|z,v\rangle`$
for every spectator v. Let F implement the fixed diagonal
frame of $`D_S`$. The states $`F|z,v\rangle`$ and
$`FJ_B|z,v\rangle`$ are Pauli translates of fixed stabilizer
states. Their support directions are fixed and their
origins affine in $`(z,v)`$. If the directions differ,
every fiber is occupied. Otherwise requiring equal supports
for every v gives an affine locus $`H_0`$ in z. In that
locus use the first affine origin as the shared origin
and a fixed linear inverse on the common support direction
to lift the other coordinates. The ratio of the two
phases, including $`D_S`$ on the second state, has a
modulo-eight weighted-cubic exponent in $`(z,v,u)`$,
where u runs over the shared support coordinates:
stabilizer phases are weighted quadratic and variable
Pauli translations add only affine and quadratic terms.

For each fixed $`(v,u)`$, a phase differing from minus
one is witnessed by a nonzero digit of the exponent minus
four. The three digit degrees are at most 1, 2 and 4 by
the [carry bound](FOUR_LAYER_SCALAR_GAP.md#a-gap-for-the-choi-vectors).
Extend these witnesses by a fixed affine retraction onto
$`H_0`$, and include the affine violations of $`H_0`$.
Their supports have union exactly the occupied set.
In a union of degree-at-most-four supports, a closest
pair in separate radius-five components would have
distance at least six. The connecting subcube has no
occupied interior by minimality. A witness at one
endpoint would restrict to one vertex or two opposite
vertices, of degree at least five. Thus the occupied
set is connected.

The scalars of BA therefore coincide on Q. The Hermitian
involutions imply ratio plus or minus one, and the
negative case has $`\varepsilon=\sqrt2`$. For
$`\varepsilon\lt\varepsilon_{\rm occ}=\min\{g_{\rm pal}/2,1\}`$
this proves exact occupied agreement.

### Apply the comparison to the rare cycle

Let $`A_1`$ be the first native band's isotropic axis span.
For each Hermitian P in $`A_1^\perp\cap Z^\perp`$, the first
band drops from the Pauli conjugation. At the cut after the
third band the compared images are two-band images as above.
Both commute with S: this is immediate for the backward
image, and exact clean return identifies S with the forward
clean image through the original three bands. The actual
reverse isometry and a commutator telescope yield

```math
 \|(VPV^\dagger-P)(I+Z)/2\|_{2,\tau}
 \le\sqrt2\bigl(\eta+\sqrt{3/2}\,2^{-p/2}\bigr).
```

The occupied normalization gives two isometry errors eta
and twice the rare-cycle distance from identity. Conversion
to physical trace normalization contributes $`1/\sqrt2`$.
For $`\eta\lt\min\{g_4/4,\varepsilon_{\rm occ}/(2\sqrt2)\}`$
and sufficiently large p, every such probe freezes.
Empty bands can be inserted for a shorter word.

The occupied unitary $`V_0=J^\dagger VJ`$ now has Pauli
support in the projection of $`A_1\cap Z^\perp`$ in the
quotient $`Z^\perp/\langle Z\rangle`$. Indeed the physical
orthogonal complement of the probe span before restriction
is $`A_1+\langle Z\rangle`$. Its projected span is isotropic.
So a logical-dirty Clifford diagonalizes $`V_0`$ over the
native ring, giving eighth-root eigenvalues. Their ratios
cannot approximate the order-three eigenvalue in the rare
cycle at the same error, by the spectral comparison of
Section 13. Actual reversal conjugates that eigenvalue and
preserves the conclusion, excluding either third-level
exterior image. The case with two non-third-level images
remains open.

<a id="native-interaction-spectral-volume"></a>

## 15. A spectral-volume bound for native interaction components

Collect the Cliffords in a d-band native word as
$`V=C E_d\cdots E_1`$. Each $`E_j`$ is the actual
transported band of independent commuting Pauli rotations.
Write $`\omega=e^{i\pi/4}`$ for the native phase root.
Connect two labels among these axes and the input clean
Pauli $`S=Z_c`$ when they anticommute. For each connected
component let $`r_i`$ be half the symplectic rank of its
span, and put $`r=\max_i r_i`$. This graph may have many
components and arbitrarily large total rank. Its labels
come from the actual native word.

At complete initialized-isometry error $`\eta\le1/64`$,
the rare cycle of Section 13 on p prefix bits satisfies

```math
p\le 2^{O((d+1)(r+1))}.
```

Thus bounded component rank requires logarithmically growing
native depth. With total symplectic half-rank R, a separate
bound is $`p+1\le 7+3R+d(3R+\binom R2)`$. Both bounds
allow arbitrary dirty width and Clifford interlayers.
Large components, even at small d, leave the general
depth gap open.

### Reduce to component programs

Let M span the labels, including S, and let
$`U=E_d\cdots E_1`$. For a Pauli P with label in
$`M^\perp`$, $`VPV^\dagger=CPC^\dagger`$ is Pauli.
On the code, $`\|VJ-cJ\|_{2,\mathrm{occ}}\le\eta+\sqrt{3/2}\,2^{-p/2}`$, where the norm uses
normalized trace on the initialized input. Comparing the
two Pauli actions first on VJ and then on cJ gives
$`\|(CPC^\dagger-P)J\|_{2,\mathrm{occ}}\le4(\eta+\sqrt{3/2}\,2^{-p/2})`$.
Two Pauli code columns are orthogonal in this inner
product or coincide up to sign, by tracing against
$`(I+S)/2`$. At large p this bound forces equality, so
$`CPC^\dagger=P S^{f(P)}`$. The bit f is linear on
$`M^\perp`$ by multiplication and vanishes at S if S
belongs to that space. Extend f with $`f(S)=0`$, and
represent it by the symplectic pairing with a Pauli
label Q that commutes with S. The Clifford

```math
D=(I+S)/2+(I-S)Q/2
```

uses a signed Hermitian representative Q. It fixes J
pointwise and conjugates $`P S^{f(P)}`$ back to P.
Thus $`C'=DC`$ belongs to the Pauli algebra of M,
and $`J^\dagger DVJ=J^\dagger VJ`$. The finitely small
p omitted by this step satisfy both stated bounds.

In a fixed Clifford frame, the radical of M consists of
computational Z coordinates z. The nondegenerate quotient
splits as a symplectic direct sum of the component spaces,
because their spans are pairwise orthogonal. The kernel
of each component's map to the quotient is exactly its
own radical. Choose $`r_i`$ active qubits in component i,
and ignore the remaining spectator factor on which the
algebra is identity. A radical S restricts z to an affine
hyperplane. Otherwise choose S itself as an active Z of
its component. The clean condition then fixes that
one active bit. The occupied z sectors have equal weight.

At fixed z each component band is diagonal in a fixed
active Clifford chart. Its phase on at most $`r_i`$ active
eigenbits t is an exponent modulo eight, a sum of signed
affine parity bits of $`(z,t)`$. Expanding a parity
$`a(z)\oplus b\cdot t`$ in t gives coefficient
$`(1-2a(z))(-2)^{|I|-1}`$ on each nonempty squarefree
monomial supported on $`I\subseteq\mathrm{supp}(b)`$.
Thus cubic coefficients are constant, each quadratic
coefficient has one varying affine bit, and each linear
coefficient has varying bits of degrees at most one
and two. The three digits of the scalar coefficient have
degrees at most one, two and four by the parity-sum
carry bound. These bounds do not count the number of
transported axes.

The Clifford $`C'`$ fixes every radical Z. On an active Pauli
its conditional symplectic action is independent of z,
and the signs on its generators are affine in z. Hence
$`C'_z=\omega^{q_0(z)}C'_0 P(z)`$, where P is an ordered
product of active Pauli generators with affine controls.
Indeed those controls are Clifford gates. Removing them
and the fixed active Clifford leaves a diagonal Clifford
on the radical, whose exponent $`q_0`$ is weighted quadratic.
Choose the order of P by component. The purely radical band
phases and this scalar may be added, with the same degree
bounds one, two and four for the total scalar digits.

A tuple determining one component's exact active matrix
therefore contains bits with sum of degree bounds at most

```math
A_i=2r_i+d\bigl(7+3r_i+\binom{r_i}{2}\bigr).
```

It includes P's affine activations and the varying phase
digits of every band. Any fixed attained tuple recurs on
at least $`2^{-A_i}`$ of the occupied radical sectors.
Multiply the appropriate bits or complements to obtain
its nonzero Boolean indicator, and apply the minimum
weight $`2^{-a}`$ of a nonzero degree-a Boolean polynomial.
The same argument for several attained tuples adds
their degree bounds. It still holds on the clean affine
hyperplane, when present.

### Factor the compression at each radical sector

The normalized Choi vector of $`U_z J_z`$ is a product
across components. Here $`J_z`$ is identity except in the
possible active clean component. The normalized Choi
vector of $`(C'_z)^\dagger J_z`$ is a stabilizer vector,
and all its z versions differ by products of local
Paulis and a scalar. Since the full initialized comparison
to cI in normalized Hilbert–Schmidt norm is small,
there is some occupied z for which these two Choi
vectors have overlap greater than $`1/\sqrt2`$.
A stabilizer Choi vector entangled across a component
cut has largest Schmidt coefficient at most $`1/\sqrt2`$:
the reduction of its stabilizer projector is a scalar
times a projector. It would have at most that overlap
with a product vector. Thus the Choi vector factors
across every cut, for all z.

The factors are stabilizer isometries, as their reference
marginals are maximally mixed. On a full component this
is a Clifford up to phase, since the maximally entangled
stabilizer pairs each input Pauli with an output Pauli.
Choose literal stabilizer phases in the factors. Their
product differs from the original by an eighth-root scalar,
by comparing a nonzero computational amplitude. Absorb
that scalar in the radical phase. The compression now
has the form

```math
B(z)=\omega^{q(z)}B_*(z)\otimes\bigotimes_{i\ne *} G_i(z).
```

Here $`B_*`$ is the contraction of dimension at most $`2^r`$
on the active clean component, or scalar one if absent.
Each full factor $`G_i`$ is unitary, a literal fixed
Clifford and affine Pauli followed by its local bands.
Their determinants are eighth roots. The factor matrices
are determined by the respective program tuples.

### A spectral outlier must recur

For $`B=J^\dagger VJ`$, the bound $`\|B-cW\|\le\eta`$
implies that at most a $`2^{-p-1}`$ fraction of singular
values of $`B-cI`$ exceed eta, since that is the rank
fraction of $`W-I`$. At least one is at least
$`\sqrt3-\eta`$ by the rare cycle's spectrum.

Each full factor that is nonscalar has spectral chord
diameter at least $`\Delta=2^{-(2d+1)r}/4`$.
In fact its matrix entries use denominator at most
$`(\sqrt2)^{(2d+1)r}`$ in $`\mathbb Z[\omega]`$:
each band is a diagonal eighth-root phase between
two fixed active Clifford charts. This ledger holds
even for dependent projected axes. The remaining
Clifford costs one such denominator. The nontrivial
field embedding fixing i sends the matrix to another
unitary. For a nonzero off-diagonal entry or a nonzero
difference of diagonal entries, the algebraic norm
bounds its absolute value below by
$`2^{-(2d+1)r}/2`$, using absolute value at most two
under the embedding. Distance to a scalar is at least
half that amount, and distance to one of its eigenvalues
is at most the spectral diameter.

If some full factor has diameter at least 1/8 at a
sector z, select it. Otherwise lift each spectrum to
a short angular interval of width $`w_i\lt1/3`$.
Each positive width is at least Delta. If their sum is
at least 1/8, select nonscalar factors up to the first
partial sum that reaches 1/8. Their product has two
eigenvalues separated by at least $`2\sin(1/16)\gt2\eta`$,
from the angular endpoints. In either case at most
$`m=1+1/(8\Delta)`$ factors are needed. Force their
program tuples to recur. For every choice of the
remaining eigenvalues, the two selected products
give two multiples of the same clean block. If both
were within eta of cI, their scalar ratio would be
within $`2\eta`$ of one by subtraction. Thus at least
one gives a singular outlier. In each recurring
z sector the fraction above eta is at least
$`2^{-mr-r}`$. Multiplying by the recurrence
fraction proves $`p+1\le m(\max_i A_i+r)+r`$.

In the remaining case $`\sum_i w_i\lt1/8`$ throughout
the occupied sectors. Replace each full factor by
the scalar with its mean angle in the chosen interval.
The scalar raised to $`2^{r_i}`$ is the determinant,
so the product center, including the radical phase,
is in $`\mu_{8\cdot2^r}`$. The product is within 1/8
of its center in every z sector. Each center digit of
one factor is a Boolean function of its program bits
and has degree at most $`A_i`$. When adding all center
exponents, bit j has degree at most
$`2^j\max\{4,\max_i A_i\}`$ by the same carry bound,
regardless of component count. Force the total center
and clean program from the sector of the largest outlier
to recur. This costs degree
$`O(2^r(1+\max_i A_i))`$. Each recurring sector has at
least one singular outlier per full-factor eigenvalue
choice, since replacing the product center changes
each such clean block by less than 1/8. The resulting
volume gives $`p\le2^{O((d+1)(r+1))}`$ in both cases.

For the total-rank estimate use a single R-qubit active
block and leave its Clifford unfactored. Force its
whole program at an outlier to recur. Only the
varying linear and quadratic coefficients, the active
Pauli controls and a single total scalar need be fixed.
Their degree cost is $`7+2R+d(3R+\binom R2)`$.
The recurring block has outlier fraction at least
$`2^{-R}`$, giving the asserted estimate.

### The same components constrain literal phase precision

The component argument has a simpler zero-clean consequence
for a full-space approximation $`\|V-\omega_3 I\|\le2^{-L}`$,
where $`\omega_3=e^{2\pi i/3}`$ is literal. Use the same
interaction graph without the clean Pauli. At sufficiently
small error each Pauli in $`M^\perp`$ is exactly fixed by
C, because the distinct Pauli operator-norm gap exceeds
the conjugation error. Hence C already lies in the
algebra of M. In every radical sector its Choi overlap
with the product of local native bands exceeds
$`1/\sqrt2`$, so the same factorization applies on the
full space. Fix all other eigenvalues to see that each
local factor has diameter at most $`2^{1-L}`$.
If $`(2d+1)r\lt L-3`$, the arithmetic diameter bound
forces each to be scalar. Then every radical block is a
scalar native unitary, hence an eighth root, a contradiction.
Thus $`d(r+1)=\Omega(L)`$ for d at least one.
This applies to a literal active control sector if its
read-only specialization gives an actual zero-clean
native word. It constrains the
[finite phase target](COHERENT_COUNTER_SOURCE.md#private-gauge-mod3)
when all interaction components are small. A large
coupled component lies within the allowed discovery budget.
