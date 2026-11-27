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
The [four-band result below](#four-band-complementary-spans)
covers complementary axis spans even at large middle crossing rank.

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
leaf phases. The limiting error floor on this corner family is sharp:
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
typical complete corners thus require exponential rank at all three
adjacent pairs if achieved in four bands. The forward and reverse
probability estimates can use half of the same failure budget.
This inference uses the complete initialized action of the frame and
its adjoint. An arbitrary moving raw decision can rotate its sign
subspaces into different addresses under reversal.
