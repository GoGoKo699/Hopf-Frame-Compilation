# Basis changes for explicit quartic sign frames at fixed accuracy

Consider a unitary Clifford+T circuit with ordered native T bands. Transport all physical
T-axis labels to one Pauli coordinate frame through the fixed Clifford intervals alone. Let
$`L_j`$ be the binary span of the axes in band j. Define

```math
\Gamma_{\rm all}=\min_{\Lambda\ {\rm isotropic}}\sum_j
       [\dim L_j-\dim(L_j\cap\Lambda)] .
```

The sum includes both outer bands. The circuit may have arbitrary T-count, dirty workspace and
initialized workspace. The angles below use the prescribed
[complete real and complex frame convention](../../docs/OPEN_PROBLEM.md).

Write $`k=n-1`$ and let $`f`$ be the binary sum of $`\lfloor k/4\rfloor`$ disjoint products of
four prefix bits. Choose the upper real Hopf angles zero and at the final level choose
$`\theta_x=\pi f(x)`$. The deepest real canonical domain is $`[0,2\pi)`$. In this convention
the block is $`R_y(\theta)=\cos(\theta)I-i\sin(\theta)Y`$,
so each final block is $`(-1)^{f(x)}I_2`$. Thus the complete frame is exactly the sign on each
prefix, tensored with identity on the final physical bit. Zero incoming node weights at these
canonical boundary angles leave the addressed complete-layer product well defined. In the
complex magnitude family the same target follows independently by putting zero in every
magnitude angle and selecting leaf phases $`\pi f(x)`$.

Suppose the emitted $`V`$ has initialized-isometry error $`\epsilon\lt\sqrt2`$ to this target
even up to a single global phase, on all logical inputs with complete dirty return. In
particular fix the final physical bit and **every** other private bit to zero. For each k-bit
x the resulting input $`|s_x\rangle`$ is its own target output basis vector, with sign
$`(-1)^{f(x)}`$. Taking the real part after the allowed common phase gives

```math
\left|{\mathbb E}_x(-1)^{f(x)}\langle s_x|V|s_x\rangle\right|
   \ge 1-\epsilon^2/2 .
```

For the disjoint-block quartic, Sections 3 and 5 give

```math
\Gamma_{\rm all}\ge
 \frac{\lfloor k/4\rfloor}{8}\log_2\!\frac{1}{\,1-2(20160/65536)\,}
   +\log_2(1-\epsilon^2/2).
```

In particular the necessary rank grows linearly with n at fixed accuracy. The exact
construction in Section 4 attains this order with constant T-depth.

## 1. Diagonal frames for the native bands

Take an optimal $`\Lambda`$, extend it to a maximal isotropic space, and choose fixed Clifford
coordinates putting it in Z. On those coordinates the jth band's X projection has dimension
$`h_j=\dim L_j-\dim(L_j\cap\Lambda)`$. It can be diagonalized with at most $`h_j`$ Hadamards.
First a classical invertible linear change places its X-row space on $`h_j`$ pivot bits, and
choose basis rows with one X pivot each. Remaining Z-only rows have no pivot Z by commutation.
The pivot rows' Z entries within the pivots form a symmetric matrix, again by commutation. A
symmetric phase shear using S and CZ cancels their entire Z content, extending the pivot rows
to a symmetric square matrix on all coordinates. Hadamards on these $`h_j`$ pivots now make
every generator Z-only. The axes originally used in the T rotations keep their actual signs
and all terms. The auxiliary span basis only selects the Clifford frame.

Absorb the two end basis changes, which may be arbitrarily large, into Clifford factors
$`C_R,C_L`$. For the remaining changes of basis use the stated circuits, placing both frames
at each internal transition. This gives $`H\le2\sum_j h_j=2\Gamma_{\rm all}`$ middle Hadamards
among monomial Cliffords and weighted-cubic diagonal phases. A T band in a Z frame produces
$`\zeta_8`$ raised to a sum of signed binary parities. Written as integer multilinear
polynomials modulo eight, parity functions have degree-j coefficients divisible by
$`2^{j-1}`$, so the resulting phase is weighted cubic. Quadratic Clifford signs and
fourth-root linear phases have the same property. Affine substitutions preserve these weighted
powers of two and discard terms of degree above three modulo eight.

## 2. Exact path bound with both boundary Cliffords

Expand each of the $`H`$ middle single-qubit Hadamards by its binary displacement b, so an
input bit t goes to $`t+b`$ with coefficient $`2^{-1/2}(-1)^{t(t+b)}`$. For a fixed path of
displacements the middle word is $`2^{-H/2}`$ times an affine computational permutation with
entries $`\zeta_8^{q(v)}`$ times a constant phase, for a weighted-cubic q of its initial
physical coordinate v. The equation for the successive computational labels is affine,
including all monomial Cliffords.

For each x, the two endpoint stabilizer states are $`C_R|s_x\rangle`$ and
$`C_L^\dagger|s_x\rangle`$. Each is flat on an affine support with fixed direction and offset
affine in x, with a jointly weighted-cubic phase on x and its free coordinates. An explicit
way to see the dependence and its global scalar is to transport X(x) through each fixed
Clifford from the all-zero input, obtaining a Pauli with affine labels and fourth-root phase
of quadratic type. Its action on the fixed zero-input stabilizer wavefunction produces the
asserted joint phase, including the relative scalar between different x. The fixed zero-input
stabilizer has constant-magnitude amplitudes with quadratic Clifford phases on its affine
support. Equivalently one can use the standard affine quadratic form of Clifford matrix
elements.

For any fixed middle path, its affine permutation carries the right support to a support
having a fixed direction. An intersection with the left support is possible on an affine set S
of x. On S the intersection is affine of a fixed dimension $`t\le\min(h_R,h_L)`$, where
$`h_R,h_L`$ are the dimensions of the two support directions. An affine section in x and t
free bits u parametrizes it. After summing over the common points the path contributes to the
diagonal entry, up to a constant phase,

```math
2^{-H/2}2^{\,t-(h_L+h_R)/2}\,{\bf1}_{S}(x)\,
             {\mathbb E}_{u}\zeta_8^{P(x,u)} .
```

The second factor is at most one. The phase exponent for fixed $`u`$ is weighted cubic in x
and has such an extension to all x by an affine retraction to S; an empty S contributes zero.
This step prevents an arbitrary external Clifford or an arbitrary number of auxiliary bits
from enlarging a path's norm.

## 3. Correlation with a cubic phase

The sign $`\chi=(-1)^f`$ has correlation with any weighted-cubic phase at most its $`U^4`$
norm. Here $`\|\chi\|_{U^4}^{16}`$ is the mean product of the sixteen vertex signs over a base
point and four independent uniform directions. Four applications of Cauchy–Schwarz give the
correlation bound because the fourth multiplicative derivative of a weighted-cubic phase is
one. To verify the latter directly, a monomial of degree j multiplied by $`2^{j-1}`$ in a
phase modulo eight has additive-difference degree at most three: a difference decreasing
ordinary degree loses at least one variable, and each surviving change at the same ordinary
degree acquires a factor of two. Equivalently, polarizing an integer parity raises 2-adic
divisibility at every extra derivative. The bound also holds multiplied by any affine-support
indicator because its Fourier expansion has absolute coefficient sum one (or zero on the empty
set). Its linear characters are themselves weighted-cubic phases.

For one quartic monomial on four variables the fourth derivative of its exponent in four
independent uniform directions equals their four-by-four determinant over $`\mathbb F_2`$.
This follows by distributing four separate differences over the product; all terms using a
variable more than once vanish at order four. The probability the determinant is one equals
$`|\mathrm{GL}_4(2)|/2^{16}=20160/65536`$. Hence its fourth uniformity norm to the sixteenth
power is $`1-2\cdot20160/65536`$. Uniformity norms multiply across independent blocks by
factorization of their defining expectations. Summing the $`2^H`$ path bounds gives

```math
\left|{\mathbb E}_x\chi(x)\langle s_x|V|s_x\rangle\right|
 \le 2^{H/2}\|\chi\|_{U^4}\le2^{\Gamma_{\rm all}}\|\chi\|_{U^4}.
```

Combining with the initialized overlap proves the claim.

Some sign tables give a stronger bound. There are at most $`8^{1+k+\binom k2+\binom k3}`$
weighted-cubic phases. For a uniformly random truth-table sign χ, Hoeffding on the real and
imaginary parts plus a union bound gives a simultaneous correlation $`O(k^{3/2}2^{-k/2})`$
with all of them. The path proof used only this correlation bound and arbitrary affine-support
indicators with Fourier absolute-coefficient sum one. Choosing one such target yields
$`\Gamma_{\rm all}\ge k/2-O(\log(k+2))`$ at fixed $`\epsilon`$. There is no restriction on how
the Clifford boundaries are selected for this truth table.

## 4. Consecutive diagonal frames

An intrinsically sharper measure for the path proof is

```math
H_{\min}=\min_{\Lambda_j\supseteq L_j}
       \sum_{j=1}^{D-1}\bigl(q-\dim(\Lambda_j\cap\Lambda_{j+1})\bigr),
```

where each $`\Lambda_j`$ is a full Lagrangian in the physical $`2q`$-dimensional Pauli space.
For any chosen completions select Clifford diagonal frames for them. The active Fourier
dimension of the transition between two neighboring frames is
$`q-\dim(\Lambda_j\cap\Lambda_{j+1})`$. Indeed, the rows of either Lagrangian modulo the other
have exactly that X rank in its coordinates. The same pivot diagonalization used above sends
one full Lagrangian to Z with that many Hadamards. Any residual Clifford preserving all Z
labels as a subspace is monomial. This proves the required decomposition independently at each
transition even if the frames on both ends had already been chosen.

The overlap bound therefore holds with $`2^{H_{\min}/2}`$ in place of
$`2^{\Gamma_{\rm all}}`$. In particular

```math
H_{\min}\ \ge\ 2\log_2\big((1-\epsilon^2/2)/\|\chi\|_{U^4}\big)
```

for the explicit witness, and $`H_{\min}\ge k-O(\log(k+2))`$ for a worst-case truth table at
fixed accuracy. To see $`H_{\min}\le2\Gamma_{\rm all}`$, fix a common optimal full Lagrangian
$`\Lambda`$. For each j take the Lagrangian completion
$`\Lambda'_j=L_j+(\Lambda\cap L_j^\perp)`$. It has intersection codimension $`h_j`$ with
$`\Lambda`$, by the perfect symplectic dual pairing and the dimension formula. The triangle
inequality for intersection codimension gives
$`\mathrm{dist}(\Lambda'_j,\Lambda'_{j+1})\le h_j+h_{j+1}`$.

For the explicit quartic family the linear lower bounds on $`H_{\min}`$ and
$`\Gamma_{\rm all}`$ are both achieved in order of growth, even at **constant** native
T-depth. For each four-bit block use one arbitrary dirty bit t. Apply
$`\mathrm{CCZ}(t,x_3,x_4)`$, toggle t by $`x_1x_2`$ with a Toffoli, repeat the same CCZ and
apply the actual inverse Toffoli. The unknown t cancels in the product of these two phases,
leaving exactly the desired quartic sign with work restored. Use the
[seven-T, five-layer CCZ word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger)
for each phase and conjugate it by Hadamards on the target for each Toffoli. The disjoint
quartets run in parallel, using 28 T gates and one dirty bit per quartet and at most 20 T
layers overall. The whole family has $`O(k)`$ T gates, and $`\Gamma_{\rm all}\le T`$ and
$`H_{\min}\le2\Gamma_{\rm all}`$.

## 5. A trace-quartic witness with a larger explicit coefficient

For $`k=n-1\ge3`$ identify the prefix with $`\mathbb F_{2^k}`$ and choose the last real angles
$`\pi f(x)`$ for $`f(x)=\mathrm{Tr}(x^{15})`$. The complex leaf phases give the same sign. For
every unit-modulus nonclassical cubic phase $`\phi`$ on the binary address space, meaning its
fourth multiplicative differences vanish (equal one),

```math
|\mathbb E_x(-1)^{f(x)}\phi(x)|\le \min(1,(112\,2^{-k})^{1/8}).
```

This applies to every weighted mod-eight phase from an internal monomial path of the native
frame construction. It also applies when multiplied by the indicator of any affine subspace,
since that indicator has a linear-character expansion of total coefficient modulus at most
one. Affine retractions extend the conditional phase to all addresses with its third-degree
property.

Applying Cauchy–Schwarz three times bounds the eighth power of the correlation by the average
over $`h_1,h_2,h_3`$ of the modulus of the x-mean of the third derivative of $`(-1)^f\phi`$.
The third derivative of $`\phi`$ is constant. The multiplicative third derivative of
$`(-1)^f`$ is an affine binary phase as a function of x. Its x-mean has modulus one or zero
according as its linear part vanishes or not.

The linear part is the fourth derivative of f evaluated in $`(h_1,h_2,h_3,x)`$. Since
$`x^{15}=x x^2 x^4 x^8`$ and the four factors are binary-linear, the polarized coefficient is
the trace of the 4-by-4 Moore determinant (in characteristic two all permutation signs are
one). Expand in the last row as $`\sum_{j=0}^3 C_j x^{2^j}`$. Trace duality and the fact a
nonzero trace pairing is nondegenerate show that the linear form vanishes exactly if

```math
C_0^8+C_1^4+C_2^2+C_3=0.
```

Trace invariance under squaring writes the linear part as the trace pairing against
$`\sum_j C_j^{2^{k-j}}`$. Raising this coefficient to the eighth power and using $`z^{2^k}=z`$
gives the equation above. The displayed polynomial in the three uniform field directions is
nonzero as an ordinary polynomial: each minor $`C_j`$ has disjoint determinant monomials and
is homogeneous of degree $`15-2^j`$. After the indicated powers the four degrees are 112, 52,
22 and 7, so the leading part cannot cancel. The elementary polynomial zero bound over the
finite field gives the exceptional probability at most $`112/2^k`$. This proves the
correlation statement, trivial for bounded small k. For the disjoint-block witness in Section
3, the third derivative is also affine. Its linear part vanishes on one four-bit block exactly
when the three directions there are dependent. Its probability is
$`1-(15\cdot14\cdot12)/16^3=1-2(20160/65536)`$. Independence across blocks gives the eighth
root of the same block product in place of the sixteenth root. This strengthens the
coefficient while preserving the linear order.

For an $`\epsilon`$-accurate native frame let $`H_{\min}`$ be the minimum total
stabilizer-basis changes for band-by-band diagonalization (sum of the Lagrangian intersection
codimensions) and let $`\Gamma_{\rm all}`$ be the single common isotropic-space distance
parameter. Every path bound with $`\ge1-\epsilon^2/2`$ target correlation gives

```math
\Gamma_{\rm all}\ge (k-\log_2 112)/8+\log_2(1-\epsilon^2/2),
 \qquad
 H_{\min}\ge (k-\log_2 112)/4+2\log_2(1-\epsilon^2/2).
```

The constant-depth realization in Section 4 applies to the disjoint-block quartic.
