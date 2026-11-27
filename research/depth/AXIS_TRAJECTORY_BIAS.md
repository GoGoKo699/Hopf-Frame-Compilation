# Bias entropy from Pauli axis relations and components

[Quantitative decision bias](ROBUST_DECISION_BIAS.md) ·
[Temporal axis cuts](TEMPORAL_AXIS_BIAS.md) ·
[Open depth question](../../docs/OPEN_PROBLEM.md#large-width-t-depth)

A large crossing rank can coexist with simple Pauli trajectories.
This chapter bounds the bias from the linear relations in the whole
rotation list and from its decomposition into independent commuting
components. The bounds allow arbitrary Clifford walls and address
motion. They use only the dirty-averaged measurement, so a complete
signed-gap source must satisfy the same obstructions.

## Fourier coefficients in a common Pauli frame

Move the actual Clifford walls to the output side. Each T or
T-dagger becomes a rotation $`e^{-i\varepsilon_j\pi A_j/8}`$
about a signed Hermitian Pauli, with $`\varepsilon_j\in\{1,-1\}`$.
The constant gate scalars cancel in conjugation. Rotations within
a band commute and may be ordered arbitrarily. Write U for the
ordered rotation list and S for the signed output Pauli in this
frame. For m computational input bits, Walsh expansion gives the
coefficients of the averaged bias

```math
 \beta_z=\tau(SU P_z U^\dagger),\qquad
 \tau(A)=2^{-W}\mathsf{Tr}(A).
```

Here W is the physical width and $`P_z`$ is the signed product
of the input Z generators selected by $`z\in\mathbb F_2^m`$.
In the section $`\sigma(u,v)=X^u Z^v`$, its binary labels
are linear in z. If the ordered generators have sections
$`i^{s_i}X^{u_i}Z^{v_i}`$, their product phase is

```math
 q(z)=\sum_i s_i z_i+2\sum_{i\lt j}(v_i\cdot u_j)z_i z_j
 \pmod 4 .
```

This convention also applies separately to tensor components,
whose local Pauli sections can be non-Hermitian. The coefficient
descriptions below allow an affine input label and a separately
stored quadratic scalar in z. They use uniform error on
coefficients, including when a local array is complex.

<a id="total-axis-nullity-net"></a>

## 1. Nullity of the complete rotation list

Let t be the number of rotations and $`a_j`$ their binary
Pauli labels, counting every actual occurrence. Define the
**axis nullity** of the representation by

```math
 D=t-\dim\mathrm{span}\{a_1,\ldots,a_t\}.
```

At nullity at most D, the coefficient arrays admit a uniform
$`\delta`$-net with logarithmic cardinality
$`O(n_0^2\log(n_0+2))`$, where
$`n_0=m+D+\log(1/\delta)+1`$ and
$`0\lt\delta\le1/2`$. There is no width, gate-count or
native-depth limit. Under an additional count cap T, an
**exact description** costs $`O(m^2+(m+D+1)T)`$ bits.

**Pauli paths.** A rotation acts trivially on a commuting
Pauli P. On an anticommuting one it gives the two terms
$`P/\sqrt2`$ and $`-i\varepsilon_j A_jP/\sqrt2`$.
Let $`w\in\mathbb F_2^t`$ specify which second terms
are taken, and write $`p(z),s`$ for the input and
output labels. The activation and endpoint tests are

```math
 \alpha_j(z,w)=\left\langle a_j,p(z)+\sum_{i\lt j} w_i a_i\right\rangle,
 \qquad
 p(z)+\sum_i w_i a_i=s.
```

A path is valid when every one in w occurs at an active
position. A valid endpoint-matched path has coefficient magnitude
$`2^{-A/2}`$, where $`A=\sum_j\alpha_j`$
is the integer activation count. For each z at most $`2^D`$
turn words solve the endpoint equation. Jointly in z and w
that equation defines an affine space of dimension at most
m+D.

**Literal scalar.** On the matched paths fix the
non-Hermitian Pauli section above. For
$`A_j=i^{s_j}\sigma(u_j,v_j)`$ and
$`p(z)=(u(z),v(z))`$, the internal phase exponent
modulo four is

```math
 \sum_j (s_j-\varepsilon_j)w_j
 +2\sum_{i\gt j}(v_i\cdot u_j)w_iw_j
 +2\sum_j w_j(v_j\cdot u(z)).
```

The input phase and constant trace of $`S\sigma(s)`$
supply the remaining scalar. This formula follows from
multiplying the inserted Paulis in reverse chronological order.
In any affine coordinates for the path labels, an affine
Boolean bit lifted modulo four is a degree-two polynomial
with even quadratic coefficients. Quadratic terms already
multiplied by two are computed over $`\mathbb F_2`$.
Thus the internal phase needs only a quadratic-size record
in the number of affine coordinates.

**Sparse description.** Take
$`k=\lceil 2D+2\log_2(1/\delta)\rceil`$ and let H
be the affine hull of the valid matched tuples of activation
count at most k. Use the zero array if H is empty. Otherwise
take an affine basis of H from these tuples. All activation
and turn bits are affine on H, and bits vanishing on every
basis point vanish throughout H. Thus their support can be
relabeled inside a set of size at most $`(m+D+1)k`$.
Each basis tuple uses at most k activation and k turn indices.

Record those index lists and the address values on the affine
basis, together with the phase on H. The decoder interpolates
all the bits affinely, enumerates each tuple of H once and sums
valid terms for each z. The decoded address and full turn word
identify each tuple uniquely, since omitted turn bits vanish
throughout H. Invalid paths are discarded by the
activation test. All omitted valid paths have count exceeding k,
and their total contribution per coefficient is at most
$`2^D2^{-k/2}\le\delta`$. The lists cost
$`O((m+D+1)k\log((m+D+1)k+2))`$ bits, giving
the stated sparse bound.

**Exact count-bounded description.** Use the full affine
endpoint solution space, and store the full activation and
turn vectors at each affine basis point. The same decoder
now sums every valid path exactly. The vectors cost
$`O((m+D+1)t)`$ bits. The address and phase
data cost $`O((m+D+1)^2)`$, which fits
$`O(m^2+(m+D+1)T)`$ because $`D\le t\le T`$.

### Consequences for random tables

For one clean bit put $`m=a+1`$ and
$`\delta=2^{-a-2}`$. By Parseval, taking real parts,
restricting the clean coordinate to zero and clipping to
$`[-1,1]`$, the coefficient net gives
normalized $`L^2`$ error at most $`2^{-a/2}/2`$
on the tested slice. For any fixed clipped profile g and
independent uniform table signs $`s_x`$, Hoeffding's bound gives
$`Pr[mathbb{E}_x s_xg(x)ge u]leexp(-2^a u^2/2)`$.
A union bound gives the signed-margin estimate

```math
 O\!\left((a+D+1)\sqrt{\log(a+D+2)}\,2^{-a/2}\right)
```

with probability at least $`1-e^{-a}`$, simultaneously
for all D. For example assign failure budget $`e^{-a-D-1}`$
to the class of nullity at most D and sum. Thus constant margin
needs $`D=\Omega(2^{a/2}/\sqrt a)`$, independent of
depth and width. The exact count-bounded description strengthens
this to $`D=\Omega_{C,\gamma}(2^{a/2})`$ for margin
$`\gamma\gt0`$ at $`T\le C2^{a/2}`$.
The latter union over $`D\le T`$ adds only
$`O(\log(T+1))`$ to the description exponent.

This tests a different parameter from middle crossing rank.
For instance disjoint X and Z banks on the same wires can
have full crossing rank and an independent concatenated list.
Such a list cannot achieve the typical constant margin.

<a id="orthogonal-axis-component-net"></a>

## 2. Orthogonal components of the rotation list

Split the axis span into an independent direct sum of pairwise
symplectically orthogonal spaces $`L_j`$, with every axis
assigned to one component. Write $`2s_j`$ for its symplectic
rank and $`d_j`$ for its radical dimension. The components
can be put in disjoint Clifford coordinates on
$`r_j=s_j+d_j`$ qubits each. Split off the nondegenerate
pairs, then complete the joint isotropic radical basis by
symplectic partners assigned to its components. If a component
has $`t_j`$ actual axes then $`r_j\le\dim L_j\le t_j`$.

There is a finest partition at the axis-label level. Choose
a linear basis of the labels and join the axes in each of its
fundamental dependence circuits. Add edges between any two
anticommuting axes. Distinct connected components now use disjoint
subsets of the chosen basis and have orthogonal spans. Every
admissible partition must keep each of the specified dependence
circuits and anticommuting pairs together.

In these coordinates U is a tensor product of rotation words
$`U_j`$, with identity on the spectators. Therefore

```math
 \beta_z=i^{q(z)}c_{\rm sp}(z)\prod_j c_j(z),\qquad
 c_j(z)=\tau_{r_j}(S_j U_j P_{z,j} U_j^\dagger).
```

Here the local Pauli sections have affine label maps in z,
and $`c_{\rm sp}`$ is zero unless the spectator input
label matches the output label. It is constant on that affine
match domain. Encode the domain, the spectator scalar and the
phase polynomial in $`O(m^2)`$ bits. Each component
may be encoded directly through its literal local rotations and
endpoint labels in $`O((m+r_j+1)t_j)`$ bits. With total count
cap T and $`r=\max_jr_j`$, this gives exact description
length $`O(m^2+(m+r+1)T)`$ including list delimiters.
Typical constant-margin tables at $`T\le C2^{a/2}`$ must
therefore have a component of width $`\Omega_{C,\gamma}(2^{a/2})`$.

**A width-limited variant without a total count cap.** By Pauli
orthonormality the squared moduli of $`c_j`$ over all
distinct local input labels sum to one. At most one label
has coefficient modulus greater than $`1/\sqrt2`$. At a
point with $`|\beta_z|\ge\delta`$, at most
$`2\log_2(1/\delta)`$ components can fail to have
that dominant label. Take an affine basis of the hull of
these significant points. Keep the components that fail at
any of its at most m+1 basis points. Every omitted component
has a constant label and hence contributes a constant on
this hull. Outside the hull the exact coefficient has modulus
below $`\delta`$.

Thus without a count cap only
$`O((m+1)\log(1/\delta))`$ components need explicit
descriptions for a uniform coefficient net. At native depth d
and local width cap r, encode each retained component in
$`O(dr^2+mr)`$ bits, using at most r independent
local axes per band. Without a depth cap, tabulate its
$`4^r`$ Pauli coefficients instead. Quantize them
and the omitted scalar inside the unit disk to total error
$`O(\delta)`$ using product telescoping. Either
variant is independent of the number of components.

<a id="same-component-square-root-constraints"></a>

## 3. Width, dependence and temporal overlap in the same component

Let $`D_j=t_j-\dim L_j`$ be the internal nullity.
Let $`b_j`$ be the maximum prefix-suffix axis
intersection dimension in the local word, at its actual
native band cuts. Omit bands empty on that component. In
each component choose the cheapest of the explicit rotation
code, the exact nullity code, and the
[exact temporal cut code](TEMPORAL_AXIS_BIAS.md). All three
apply to the local affine Pauli maps appearing above. This
gives total logarithmic cardinality bounded by

```math
 O\!\left(m^2(T+1)+T\max_j\min\{r_j,D_j,b_j\}\right).
```

The $`m^2(T+1)`$ term covers the global
endpoint and spectator conditions, the local address data
and model choice in each of at most T components. Supposing
$`T\le C2^{a/2}`$, this overhead is $`o(2^a)`$.
With probability at least $`1-e^{-a}`$, a random table
with a constant signed margin $`\gamma\gt0`$ therefore
requires **one same axis component** with simultaneous
$`r_j,D_j,b_j=\Omega_{C,\gamma}(2^{a/2})`$.
This conclusion holds for every physical representation
with count inside the cap. It follows by fixing a sufficiently
small common threshold on the minimum of the three quantities,
then applying the random-sign tail bound over the exact descriptions.

<a id="same-component-symplectic-strengthening"></a>

**Full operator gap.** For a source whose margin holds on
each entire initialized address sector, the same component
must have $`s_j,D_j,b_j=\Omega_{C,\gamma}(2^{a/2})`$, where
$`s_j`$ counts symplectic pairs. The statement also holds when
the complete operator guarantee implies such a sector margin.

For each component selected for the direct code, condition on
the part of its radical orthogonal to all m input stabilizers
and the measured output Pauli. This joint radical space commutes
with the rotations, the input constraints and the measurement.
Each selected radical contributes a subspace of codimension
at most m+1, so the local word, restricted to its eigensector,
uses at most $`s_j+m+1`$ qubits. Choose the conditioned
eigenvalue bits affinely in the address x. The only compatibility
constraints come from the binary label intersection of the joint
radical with the input stabilizers, where the required character
is affine in x.
An affine extension to a basis of the joint radical always exists.

On this compatible eigensector, take the normalized projector
onto the remaining input code. Its trace tests the required
margin by the full operator promise. Its Pauli expansion uses
independent projected input generators with affine eigenvalue
bits in x. The output sign and each selected rotation sign
acquire affine bits in x from the conditioning. The unselected
components are unchanged and take the exact nullity or temporal
code, now on the projected input generators. In the description
of the factorized coefficients, the remaining global phase,
input eigenvalues and output sign cost $`O(m^2)`$ bits. Each
selected component costs $`O((s_j+m+1)t_j+m^2)`$, including
the address-dependent signs. Choose the cheapest mode separately
for each component. The resulting exact profile descriptions
therefore have exponent

```math
 O\!\left(m^2(T+1)+T\max_j\min\{s_j,D_j,b_j\}\right).
```

The conditioned scalar profile still takes values in
$`[-1,1]`$, so the same random-sign count applies. The
conditioning uses the freedom to test arbitrary dirty states
in each address sector.

The result concerns representation structure. It applies
also to the complete deepest-corner frame through the same
output measurement as [the three-band proof](ROBUST_DECISION_BIAS.md#a-direct-complete-frame-error-floor).
Tensor repeats of the [exact cancellation witness](ROBUST_DECISION_BIAS.md#exact-cancellation-at-large-crossing-rank)
can have arbitrarily large aggregate crossing and nullity with
constant component widths. A single large component with
large dependence and large temporal overlap can still fit
within the square-root T budget, and remains outside these
structural screens.


<a id="linear-trace-cut-dependence"></a>

## 4. Linear traces require dependence within a half

A complete returned frame also gives a linear trace test. In this
test the small-nullity condition can be applied separately to the
two halves of any chronological cut. Define
$`e=e_L+e_R`$, where each $`e_H`$ is the number of rotation
occurrences in half H minus the rank of their binary axes. The halves
may share arbitrarily many axes and the total Clifford is arbitrary.

For a probe family $`P_z`$ indexed by $`z\in\mathbb F_2^m`$,
assume affine Pauli labels and the literal ordered-product phase
convention of Section 1. The arrays $`\tau(P_z^\dagger V)`$
admit a uniform $`\delta`$-net of logarithmic cardinality

```math
O\!\left((m+s^2+1)^2\log(m+s+2)\right),\qquad
s=e+1+\log(1/\delta),\qquad 0\lt\delta\le1/2.
```

There is no width, count or depth cap. For a typical complete
deepest-corner frame at fixed error below $`\sqrt{2-\sqrt2}`$ and
polynomially many initialized bits, every cut in a successful
representation has $`e=\Omega(2^{a/4}/a^{1/4})`$ for a-bit
addresses.

### Independent rotations give tensor products of Hadamard eigenstates

Write the word as $`V=CGF`$ in the common rotation frame, with
F the earlier half and G the later half. Retain the native scalars
from $`T^\varepsilon=e^{i\varepsilon\pi/8}e^{-i\varepsilon\pi Z/8}`$
as one overall phase. Use operator Choi space in the computational
coordinates given by the literal Pauli sections. A fixed Bell Clifford
implements this change of analytical basis.

For either F or the actual reverse $`G^\dagger`$, expand the product
of signed rotations by choosing either its identity or its Pauli
term. If t rotations have turn word u, its scalar magnitude is

```math
 \cos(\pi/8)^{t-|u|}\sin(\pi/8)^{|u|}.
```

The literal phase has exponent modulo four with even quadratic
coefficients from the ordered Pauli products, as in Section 1.
Keep an independent sublist and enumerate the omitted turns v.
The full label is $`Au+Bv`$ for an injective A on the kept
turns u, and the exponent has form
$`q_0(u)+q_1(v)+2u^{\mathsf T}Qv`$. Thus the vector for
each v is one fixed Clifford image of a tensor product of

```math
 |h\rangle=\cos(\pi/8)|0\rangle+\sin(\pi/8)|1\rangle
```

on the kept coordinates and zero on spectators, followed by a
Pauli translation affine in v. The fixed Clifford applies the
quadratic phase and the injective computational label map. The
mixed phase is a Pauli Z shift in the kept coordinates.
The omitted turns contribute known scalar weights with sum of
moduli at most $`2^{e_H/2}`$.

The normalized trace equals $`\langle G^\dagger C^\dagger P_z,F\rangle`$
in the Choi inner product. Right multiplication by the fixed
Clifford is a Clifford on this space, and right multiplication
by $`P_z`$ is a Pauli with affine labels. After pulling the fixed
coordinate Cliffords to one side, each term, or its complex
conjugate, is an overlap

```math
 a(\theta)=\langle L|P(\theta)K|R\rangle,
```

where $`\theta`$ contains z and the omitted turns, K is one fixed
Clifford, and $`P(\theta)`$ is an affine Pauli family with stored
quadratic literal phase. L and R are products of $`|h\rangle`$
and zero. K is common to the whole parameter array. It suffices
to approximate these overlaps uniformly at
$`\epsilon_0=\delta 2^{-e/2}/2`$ before summing omitted turns.

### Large overlap aligns the spectator constraints

Suppose some overlap has modulus at least
$`\epsilon_0`$. Let $`S_L`$ be the unsigned spectator Z span of L
and $`S_R`$ the image under K of the corresponding span for R.
Put $`c=\mathsf{rank}\,\omega(S_L,S_R)`$ and
$`j=\dim(S_L\cap S_R)`$. The probability in L of the
right spectator eigensector is at least $`\epsilon_0^2`$.
In the signed Pauli expansion of that projector, terms outside
$`S_L^\perp`$ vanish. The remaining subgroup has codimension c
in $`S_R`$; its projection to the active left Pauli labels
has dimension $`d=\dim S_R-c-j`$. A uniform label in a
d-dimensional binary span has d independent coordinate pivots,
and its nonidentity Pauli count is at least half their weight.
On $`|h\rangle`$ the nonidentity expectation has modulus at most
$`1/\sqrt2`$. The projector probability is therefore at most

```math
 2^{-c}\left(\frac{1+2^{-1/4}}2\right)^d .
```

Thus $`\dim S_R-j=O(s)`$, and exchanging the states gives the
same bound on the left. If all overlaps were below
$`\epsilon_0`$, use the zero array.

Use a CNOT basis change on each zero-spectator register to align
the common intersection with the same j Z wires. In these bases
K fixes each selected Z up to a fixed sign. Returning the common
zero sector under $`P(\theta)K`$ is an affine condition
on $`\theta`$. Store that condition and use zero overlap
off its domain. On its domain factor one fixed X sign correction
from K and restrict to the common zero wires. The remaining
Clifford is fixed. The residual Pauli translation is affine, and
its stored quadratic phase includes multiplication by the fixed
X correction. Its remaining common Z labels act as positive
identity on the chosen sector. Only $`O(s)`$ unmatched zero
spectators remain in the two states. Expand each using

```math
 |0\rangle=\cos(\pi/8)|h\rangle+\sin(\pi/8)iY|h\rangle .
```

This produces a sum with coefficient modulus sum
$`2^{O(s)}`$ and a common fixed Clifford between translated
full products of $`|h\rangle`$. The additional Pauli
labels are affine in the new branch bits. A significant overlap
thus supplies a full-product branch of magnitude $`2^{-O(s)}`$.

### Relative Pauli statistics compress the remaining Clifford

The squared overlap of two pure states on q qubits is
$`2^{-q}`$ times the sum of products of their Hermitian Pauli
expectations. Taking absolute values bounds it by the affinity
of their normalized squared Pauli expectations. For the full
product of $`|h\rangle`$, this is exactly the ternary product
law of the [relative product inverse theorem](RELATIVE_FOUR_BAND_BIAS.md#2-an-inverse-statement-for-the-relative-product-law).
Its one-cell masses are $`1/2,1/4,1/4,0`$ for I, X, Z and Y
respectively. Pauli translations change only the expectation signs.

Apply the inverse theorem at the significant branch. The fixed
Clifford's symplectic matrix is a whole-cell permutation with
possible X/Z swaps, plus a correction of rank $`O(s^2)`$.
The permutation and plain Hadamards fix the full product
literally. Remove that Clifford on the product-input side.
The remaining symplectic matrix $`M'`$ satisfies
$`\rho=\mathsf{rank}(M'-I)=O(s^2)`$.

The Pauli coefficients of the remaining Clifford form a
[stabilizer vector](ROBUST_DECISION_BIAS.md#2-three-bands-reduce-to-a-single-diagonal-choi-decoration)
in Bell coordinates. Its affine support
dimension is rho because the stabilizer indexed by p moves
the Bell label by $`p+M'p`$. All supported coefficients
have magnitude $`2^{-\rho/2}`$ and a weighted quadratic
Clifford phase on affine coordinates. A fixed common scalar
phase is retained, or stored to the requested tolerance.

Expand over these support coordinates as further analytical
parameters. Each product Pauli expectation has factors
$`1,1/\sqrt2,0`$ at local affine labels, and a quadratic
literal phase. At tolerance $`\zeta=2^{-O(s+\rho)}`$,
the significant terms have $`O(s+\rho)`$ nonidentity cells.
Take the affine hull of their parameter values and record
the local labels on an affine basis. Only the union of
nonidentity cells there can occur throughout the hull.
Store each sparse list, the phase polynomial and the hull.
The decoder evaluates those terms exactly and uses zero
outside the hull. The combined parameter count is
$`O(m+s^2+1)`$. The sparse lists and phases give the stated
description length. Choose the constant in the precision
exponent to cover the coefficient modulus sums of all
analytical expansions. No new initialized physical registers
are used.

### Test the complete returned frame with a linear trace

For the [deepest-corner witness](ROBUST_DECISION_BIAS.md#a-direct-complete-frame-error-floor)
put $`A=XZ`$ and $`F(x)=(-1)^{f(x)}`$. Then
$`U_f(x)=B_0+F(x)B_1`$, where $`B_0=(I+A)/2`$ and
$`B_1=(I-A)/2`$ both have norm $`1/\sqrt2`$.
Let $`v_x`$ denote the compression of V to address x and
the initialized values at both ends. Its target and dirty inputs
are unrestricted. The real overlap of the two unit-norm full
outputs gives $`\mathsf{Re}\,\tau(U_f(x)^\dagger v_x)\ge1-\eta^2/2`$
at complete error eta. Since $`v_x`$ is a contraction,
the linear profile obeys

```math
 F(x)\mathsf{Re}\,\tau(B_1^\dagger v_x)
 \ \ge\ 1-\eta^2/2-1/\sqrt2 .
```

It has positive margin below the same limiting error threshold.
The ambient profile allowing all values of the tested initialized
bits is also bounded by one. Its Walsh coefficients are real
linear combinations of the two Pauli trace arrays from
$`I+A`$. Using coefficient tolerance $`2^{-m-2}`$,
Parseval and clean-slice restriction give normalized error
at most $`2^{-a/2}/2`$. A random-sign tail bound and a union
bound over the net exclude the constant margin with probability at
least $`1-e^{-a}`$ whenever
$`e\le c_\eta\,2^{a/4}/a^{1/4}`$ for sufficiently small
$`c_\eta\gt0`$ and $`m=\mathrm{poly}(a)`$.
The cut may be chosen by the compiler, since the net already
counts all such pairs of halves. This excludes using only the
relations across a cut to explain a typical complete return.
Large dependence within the two halves can still occur at
constant depth.
