# Relative Pauli bias for four perfectly paired bands

[Quantitative decision bias](ROBUST_DECISION_BIAS.md) ·
[Complementary spans](ROBUST_DECISION_BIAS.md#four-band-complementary-spans) ·
[Literal Pauli trajectories](TEMPORAL_AXIS_BIAS.md)

Consider four actual native T bands, with all the Clifford walls
absorbed into one output Clifford. Let their binary axis spans be
$`L_1,L_2,L_3,L_4`$ in chronological order, with perfect
symplectic pairings between $`L_1,L_2`$ and between $`L_4,L_3`$.
Their two half-spans are $`W_L=L_1\oplus L_2`$ and
$`W_R=L_3\oplus L_4`$, with half-dimensions $`r_L,r_R`$.
These conditions allow large intersections of one
interior span with the sum of the other interior span and an
exterior span. Assume $`r_R\ge r_L`$ for the raw bias statement
and let $`k_4`$ count the axes of the fourth band commuting with
the fixed measured Pauli in the absorbed frame. Complete returned
frames can be oriented using their actual reverse, as in Section 5.

For $`m`$ fixed computational input bits and all other inputs
averaged, the bias profiles of these words have a normalized
$`L^2`$ epsilon-net of logarithmic cardinality
$`\mathrm{poly}(m+k_4+\log(1/\epsilon))`$ for
$`0\lt\epsilon\lt1/2`$. The same polynomial bound holds when
$`k_4`$ is an upper bound. No physical width or gate count
appears. In particular, $`k_4=\mathrm{poly}(a)`$ excludes
constant average signed margin on a typical $`a`$-bit table
with one clean signal. For a complete returned frame the
measured Pauli can be chosen so that $`k_4`$ is automatically
small, as shown in Section 5. Section 6 replaces perfect pairings
by a quantitative allowance for their rank deficiencies.

## 1. A perfectly paired half gives a product law with affine constraints

Use the orthonormal Pauli sections $`\sigma_q=X^{q_X}Z^{q_Z}`$
under normalized Hilbert--Schmidt inner product. Let $`P_z`$
be the input Pauli for Fourier label $`z\in\mathbb F_2^m`$
and let S be the measured Pauli in the absorbed frame. The
Fourier coefficient is the inner product of the output Pauli
transported backward through the last two bands with the input
Pauli transported forward through the first two bands. Each
half has distinct labels for distinct turn words, by the direct
sum in that half and independence of each physical band.

In the forward half, the input label specifies the active masks
of the first band. Given its turn word t, the activation vector
$`\alpha`$ of the second band is an affine bijective function
of t, by perfect pairing. A direction t which turns at one of
the $`k_1(z)`$ inactive exterior axes is excluded. Where allowed,
t is uniform in its coordinate subspace in the squared Pauli
coefficient distribution. Each interior turn bit w is zero
when its activation bit is zero, and otherwise uniform. If we
temporarily take every t uniform, the pair of coordinates
$`(w,w+\alpha)`$ at each interior axis has the probability law

```math
 \nu(00)=1/2,\qquad \nu(01)=\nu(10)=1/4,\qquad \nu(11)=0.
```

The true squared amplitudes are consequently a product
$`\nu^{r_L}`$, conditioned on $`k_1(z)`$
affine linear tests with total probability $`2^{-k_1(z)}`$.
The backward half uses $`\nu^{r_R}`$ with $`k_4`$ tests.
This padding to all t is a probability calculation.

An affine subspace of codimension c has $`\nu^r`$ probability
at most $`(3/4)^c`$. Eliminate the constraints in block columns,
in the order of the two-bit cells. If a cell has $`j=1`$ or 2
pivots, the echelon constraints on this cell, after fixing later
cells, have rank j. Its conditional maximum acceptance is $`3/4`$
for $`j=1`$ and $`1/2\le(3/4)^2`$ for $`j=2`$. Integrating over
cells gives the bound.

The midpoint labels can match only when $`p(z)+s\in W_L+W_R`$,
where p and s are the input and output labels. Put
$`\ell=\dim W_R-\dim(W_L\cap W_R)`$ and $`d=r_R-r_L`$.
The affine space containing the forward conditional support has
dimension $`2r_L-k_1(z)`$. Its intersection with the backward
label space therefore has codimension at least both ell and
$`2d+k_1(z)`$ in that space. The backward law on this space
is a conditioned product law with probability normalization
$`2^{k_4}`$. The absolute coefficient is at most the square
root of its mass on the forward support. A coefficient of
magnitude at least delta consequently requires

```math
 \max\{\ell,2d+k_1(z)\}=O(k_4+1+\log(1/\delta)).
```

Pad the forward product coordinates with d cells constrained to
00. This raises their squared-amplitude normalization by $`2^d`$.
The forward match domain is an affine subspace in the padded
coordinates of codimension ell, including the zero-padding
constraints. Its direction is identified with $`W_L\cap W_R`$.
The matching map to the backward domain has a fixed homogeneous
part that is an isomorphism between these directions. Extend it
to an invertible binary map M on the full padded coordinate
space. An affine choice of matching origins for z in the endpoint
match domain then gives $`y=Mx+b(z)`$, with b affine in z. The
overlap still enforces the ell match tests on x.

The two halves use the same $`\sigma_q`$ at every matched label.
Their relative phase is $`i^{q(x,z)}`$, where q is a degree-two
polynomial modulo four with even quadratic coefficients. Indeed
each turned gate on axis $`A=i^{s_A}\sigma_a`$ contributes its
signed scalar from $`-i\varepsilon A\sigma_p`$. The ordered
Pauli products supply only pairwise signs. Substituting affine
binary labels preserves the phase degree, and reversing the
rotation changes its literal orientation. The phase of the
initial and final Pauli is included.

## 2. An inverse statement for the relative product law

The affinity of probabilities P and Q is
$`\mathrm{Aff}(P,Q)=\sum_x\sqrt{P(x)Q(x)}`$.
Let $`X\sim\nu^r`$ and let M be invertible on its $`2r`$ bits.
If the affinity of X and $`MX+b`$ is at least $`2^{-s}`$ for
$`s\ge1`$, then

```math
 M=\Pi+R,\qquad \mathrm{rank}_{\mathbb F_2} R=O(s^2),
```

where the permutation $`\Pi`$ carries whole two-bit cells to
whole cells and may swap the two coordinates in each cell.
The bound allows arbitrary affine displacement b.

### Sparse rows up to a low-rank correction

A nonzero linear character on one cell has mean $`1/2`$ for
either coordinate and zero for their sum. Hence a character
touching h cells has mean of modulus at most $`2^{-h}`$.
Take $`K=\lceil C s\rceil`$ for a sufficiently large absolute C.
Greedily choose rows of M from distinct output cells whose
span has minimum nonzero source-cell weight greater than 2K.
If K rows are chosen, the law of those output bits has all
nonconstant Fourier biases bounded by $`2^{-2K}`$. Its total
variation distance from uniform is at most $`2^{-3K/2}/2`$
by Parseval. Meanwhile the same marginal of X has independent
Bernoulli $`1/4`$ coordinates. Its affinity with uniform is
$`\cos(\pi/12)^K`$. The square-root vectors of two laws at total
variation distance t have Euclidean distance at most $`\sqrt{2t}`$.
The projected affinity is therefore less than $`2^{-s}`$ for
large enough C. Projection can only increase affinity, so the
greedy choice must stop sooner.

Each row outside the used output cells is now within 2K source
cells of the selected span. Put the span corrections and the rows
of the used output cells into a matrix of rank at most 3K.
Write M as this correction plus S, where every row of S touches
at most 2K source cells.

### Isolate multi-cell rows and repair the pairing

Greedily choose a row of S outside the used output cells whose
restriction to unused source cells touches at least two cells,
or touches exactly one by the sum of its two coordinates. Mark
the row's output cell and mark two witness cells in the first
case and one in the second. Let L be the number of choices. Their
witness sets are disjoint. Make a conflict graph by joining two
chosen rows when either row touches the other's witness set.
There are at most $`2KL`$ edges, so an independent collection of
$`l\ge L/(4K+1)`$ rows exists by the average-degree bound.

Fix all source cells outside those witnesses. On the selected
rows the bits of $`SX`$ are independent and have means in
$`[3/8,5/8]`$, including after any fixed translation. In the target
law the same bits are independent with mean $`1/4`$. The event
that their mean is at most $`5/16`$ has complement at most
$`\exp(-c_0 l)`$ in the target law and probability at most
$`\exp(-c_0 l)`$ for every fixed translation of $`SX`$, by
Hoeffding. The actual correction plus b has at most $`2^{3K}`$
possible translation values. Use a union bound over those values,
without conditioning on the actual correction. For any event E,
affinity is at most $`\sqrt{P(E^c)}+\sqrt{Q(E)}`$.
This forces $`L=O(K(s+K))=O(s^2)`$.

Absorb the marked source columns and output rows into the
low-rank correction. At termination each other row is either
zero or a single coordinate. Invertibility of M implies that the
rank deficit of this one-hot matrix is at most the correction
rank. Change at most that many rows to complete a coordinate
permutation, by keeping one row per used source coordinate.

Let B count its output cells whose two coordinates come from
different source cells. Each source cell is used twice in total,
so at least $`B/3`$ broken pairs have disjoint sets of source
cells. On these pairs the forbidden outcome 11 has probability
at least $`1/16`$ under the permutation law with any fixed
translation, independently between selected pairs. The true law
can differ by at most $`2^{O(s^2)}`$ translations. Its probability
of the target support is therefore at most
$`2^{O(s^2)}(15/16)^{B/3}`$. Since affinity is at most the square
root of this probability, $`B=O(s^2)`$. Repair the broken cells
using the source cells unused by the unbroken ones. This alters
only $`2B`$ rows and proves the relative product statement.

## 3. A decoder for a quadratic phase averaged against the product law

Let Q be a fixed quadratic phase exponent modulo four with even
cross coefficients and zero constant. Let the parameter vector
$`\theta`$ have p bits and let the tilt coefficients
$`\ell_j(\theta)`$ be Boolean polynomials of fixed bounded degree.
The arrays

```math
 g(\theta)=\mathbb E_{x\sim\nu^r}i^{Q(x)+2\ell(\theta)\cdot x}
```

have a uniform eta-net with bit length
$`\mathrm{poly}(p+1+\log(1/\eta))`$ independent of r. A scalar
of bounded-degree phase in the parameters can be recorded
separately. The proof keeps the fixed polar form and allows
dense linear tilts on every physical cell.

### Find a low-rank completion of the cross-cell polar form

Sample a $`\nu`$ cell by choosing one coordinate uniformly and
then making it a uniform bit, with the other coordinate zero.
Conditional on these coordinate choices, a quadratic Gauss
average of phase Q with any linear tilt has modulus at most
$`2^{-\mathrm{rank}(M_0)/2}`$, where $`M_0`$ is the polar matrix
on the chosen coordinates. Squaring the average and differencing
the exponent proves this bound by averaging a linear character.

Let h be the maximum cross rank of the full polar matrix M over
rows and columns in disjoint sets of cells. Choose an invertible
$`h\times h`$ minor there. The random coordinate choice retains
a constant fraction of its independent rows except with
probability $`\exp(-\Omega(h))`$, by independence of cells.
On the retained rows choose independent columns from the minor.
The choices on those column cells are independent of the row
choices and keep a constant fraction of the column pivots with
the same exponential tail. Therefore the above conditional
average has expected modulus $`\exp(-\Omega(h))`$. Unless the
whole array is below eta, $`h=O(1+\log(1/\eta))`$.

Write U and V for the row and column coordinates of a maximum
invertible cross minor and remove their cells. The Schur
complement of $`M[U,V]`$ shows that between any distinct
remaining cells M agrees with the rank-h matrix
$`M[:,V]M[U,V]^{-1}M[U,:]`$ restricted to the remaining
coordinates. Otherwise an enlarged minor on disjoint cell sets
would contradict maximality. This completion is symmetric off
its diagonal cells. Its sum with its transpose has rank at most
2h and is supported on the diagonal cells, so remove at most h
additional asymmetric cells. The remaining symmetric completion
has rank at most h, and the residual polar form is cell-local.

Enumerate the removed cells, retaining their original $`\nu`$
weights and moving their cross terms to the linear tilt. The
quadratic phase with the symmetric rank-h completion as polar
form is a function of at most h linear checks, up to a fixed
even linear term. Indeed on the radical it is an additive
modulo-four character, necessarily even. Extend that character
linearly and subtract it. The resulting phase is constant on
radical cosets. Record its phase on a linear section of the
check space as a quadratic polynomial. No coordinates in the
full physical space are part of the record.

### Retain only cells with nonidentity phase factors

Expand the low-rank phase in Fourier characters of the check
space. Its coefficients are determined by the recorded small
polynomial and have sum of moduli at most $`2^h`$.
Each product after taking those shifts and fixing the removed
cells has local factors of the form

```math
 \frac12+\frac14 i^{u_j(\theta')}+\frac14 i^{v_j(\theta')},
```

where $`\theta'`$ contains the original parameters, the removed
cell values, and the check Fourier label. The two phase digits
are bounded-degree polynomials, since each local exponent is
constant plus an even linear tilt in the physical bits. The
factor is one when both phases equal one and has modulus at most
$`\sqrt{10}/4`$ otherwise. A product of magnitude at least
$`\eta'=\eta 2^{-h}`$ thus has only
$`O(1+\log(1/\eta'))`$ nonidentity cells.

Lift $`\theta'`$ to its monomials of the required fixed degree
and take the affine hull of the lifts where the product is at
least eta-prime in modulus. At an affine basis from those
points, keep the union of nonidentity cells. All phase digits on
each omitted cell vanish identically on the hull, so that cell
contributes one. Store the hull and only the polynomials for the
kept cells. Decode the exact product on the hull and zero off
it. The bit length is polynomial in $`p+1+\log(1/\eta)`$.
The normalized sum over the removed cells and the Fourier
coefficient norm show error at most eta in the reconstructed
average. If the significant lift set is empty, the zero product
has the same bound.

## 4. Describe the full coefficient array

Fix coefficient tolerance delta, and use a smaller constant
multiple of it for each truncation. Discard labels outside the
endpoint match domain or with $`k_1(z)\gt K_1`$, where
$`K_1=O(k_4+1+\log(1/\delta))`$ is large enough by Section 1.
The affine hull of the other labels has a basis of at most m+1
points. Its inactive first-axis flags are affine in z, and only
the union of at most $`(m+1)K_1`$ flagged directions at these
points can occur on this hull. Store those flags and still test
the cutoff when decoding. An empty hull gives zero coefficients.

If all coefficients have modulus below delta, the zero net point
suffices. Otherwise Section 1 also bounds ell and d by the same
order as $`K_1`$. At one significant label, dropping the
conditioning and match tests gives product affinity at least
$`\delta 2^{-(K_1+k_4+d)/2}`$. By Section 2 the same fixed M
used for every matched label is a cell permutation plus a
correction of rank $`b=O((k_4+1+\log(1/\delta))^2)`$.
Relabel the backward cells by the permutation, which preserves
their product law. Factor the remaining correction as DA with
b row coordinates. For each $`\lambda=Ax`$ the backward
coordinates then equal $`x+\Delta(z,\lambda)`$ for a fixed
affine function Delta.

A cell shifted by a nonzero two-bit vector has affinity at most
$`1/\sqrt2`$ with its unshifted law. This is obtained directly
from the three nonzero probabilities of $`\nu`$. A branch with
H shifted cells therefore has absolute contribution at most
$`2^{(K_1+k_4+d)/2}2^{-H/2}`$, even after the path constraints
and phase. There are at most $`2^b`$ branches. Discard the
branches with more than
$`H_0=O(b+K_1+k_4+d+1+\log(1/\delta))`$ shifted cells,
with total coefficient error at most a chosen constant fraction
of delta. The affine hull of all retained $`(z,\lambda)`$
has at most m+b+1 basis points. The union of displaced cells
on that basis has size at most $`(m+b+1)H_0`$. Record the
affine shift on these cells and zero elsewhere on the hull.
Decode only branches in the hull passing the original cutoffs.

For a retained branch enumerate the values of x on those
displaced cells. Their weights are the explicit factors
$`\sqrt{\nu(x_j)\nu(x_j+\Delta_j)}`$, with sum at most one
over each cell. Every other cell has its ordinary $`\nu`$
weight. The constraints are $`Ax=\lambda`$, the ell tests of
the forward match domain, the fixed inactive backward turns,
and the possible inactive forward turns with their recorded
activation flags. The number of constraints is
polynomial in $`m+k_4+1+\log(1/\delta)`$. Each is affine in
x at fixed parameters. The path cutoff has already removed
Fourier labels with too many exterior inactive flags. Enforce
an indicated constraint $`\ell(x,z,\lambda)=0`$ by the exact
two-character average over $`(-1)^{u e(z)\ell(x,z,\lambda)}`$,
where $`e(z)`$ is its occurrence bit and $`u\in\mathbb F_2`$.

The remaining expectation has the fixed quadratic polar form
from the relative phase on the remaining cells. All added phase
terms are scalar or even linear in these cells, with coefficients
of bounded degree in z, lambda, the enumerated cell values and
constraint Fourier bits. Thus Section 3 encodes their joint
expectation with a polynomial-size record. The constraint
averaging is normalized. The product weights of the removed
cells sum to at most one, so an internal uniform error eta
costs at most $`2^{b+(K_1+k_4+d)/2}\eta`$ per coefficient.
Choose eta to leave the remaining fraction of delta. The
description includes the exterior and padding conditioning
factors, the literal relative phase, and the affine cutoffs. The large
physical matrices A and D enter only through these compressed
phase and constraint data.

This proves a uniform coefficient net with logarithmic cardinality
polynomial in $`m+k_4+\log(1/\delta)`$. Take real parts and
apply Parseval at $`\delta\le\epsilon 2^{-m/2}`$ for the stated
normalized $`L^2`$ net. Restricting to fixed clean input bits and
counting against random signs works as in the
[three-band bias proof](ROBUST_DECISION_BIAS.md#5-count-descriptions-and-compare-to-random-tables).
The Fourier cutoff can be exponentially small in m. Its logarithm
still gives polynomial description length.

## 5. A complete returned frame supplies the active output probe

Consider the [deepest-corner frames](ROBUST_DECISION_BIAS.md#a-direct-complete-frame-error-floor)
and a complete approximation with operator error
$`\eta\lt\sqrt{2-\sqrt2}`$ in the geometry above. If needed,
use the actual native reverse to put the longer half on the
measurement side. For the initialized embedding J, multiply
$`VJ-JU_f`$ by $`V^\dagger`$ on the left and $`U_f^\dagger`$
on the right. This gives the same error for the actual reverse
on the code, with the conjugated common phase if one is used.
The ideal reverse has the same target signs at input zero.

Use target input zero and the standard target Z measurement, but multiply
the measured Pauli by a dirty-output Pauli R. Test inputs on
the $`+1`$ eigenspace of R. The ideal frame returns the dirty
bank, so its correct signal sign is also the sign of the modified
measurement on this input subspace. The correct-outcome probability
still has the lower bound $`(1-\eta^2/2)^2`$ from overlap with
the ideal state. This gives constant positive signed bias at
the stated fixed error.

To choose R, transport the independent last-band axes into the
physical output coordinates. Their projections onto the dirty
output coordinates have kernel dimension at most twice the
number of nondirty output wires. Choose a maximal independent
subset of those projections. The dirty symplectic pairing is
nondegenerate, so one can prescribe its pairing with R on all
members of this subset. Choose the values that make every selected
axis anticommute with the modified measurement. Thus $`k_4`$ is
at most twice the number of nondirty output wires. An input
Clifford on the dirty bank makes the selected R sector a single
computational constraint, adding one tested bit when R is
nonidentity. This does not add an initialized circuit wire and
conjugates all absorbed spans by the same symplectic map.

With at most polynomially many initialized work bits the
resulting m and $`k_4`$ are polynomial in the address length.
The coefficient net, restriction to the tested slice and random
sign bound therefore exclude these four-band words for almost
every deepest-corner frame at any fixed error below the threshold.
The argument uses ideal work return and allows actual work error
within the complete-isometry guarantee. A raw decision source
without dirty return has no corresponding freedom to change its
measured Pauli.

## 6. Allow rank deficiencies in the two pairings

Let $`r_L,r_R`$ now be the ranks of the two within-half pairing
matrices of the actual band bases. Write $`t_j`$ for the band
sizes and set

```math
 e=\sum_{j=1}^4t_j-2(r_L+r_R).
```

For $`r_R\ge r_L`$ the normalized $`L^2`$ net has logarithmic
cardinality $`\mathrm{poly}(m+k_4+e+\log(1/\epsilon))`$,
where $`k_4`$ may count inactive axes in the entire fourth band.
This also permits dependent concatenated axis lists within a half.

In each pairing matrix retain a square nonsingular minor of its
full rank. The resulting band sublists give perfectly paired
cores. The other e axes are actual omitted literals. Reorder
commuting axes within each native band to place omitted exterior
rotations nearest their endpoints and omitted interior rotations
nearest the midpoint. Use the same order from each endpoint
toward the midpoint. This leaves the word and its actual inverse
unchanged as operators.

Expand the omitted rotations by their activation and turn bits,
with at most $`2^{2e}`$ joint parameter values. A turn at an
inactive axis is discarded. All omitted activations are affine
tests on the incoming labels. The selected omitted exterior
turns shift the two starting labels affinely in the new parameters.
They leave the kept exterior activation masks unchanged, since
axes within the same band commute. The omitted interior turns
shift the midpoint matching equation. Their activation tests
involve only the incoming label and exterior turns, because the
kept and omitted interior axes commute within their band.

At each parameter, the two kept halves have precisely the product
laws in Section 1. The omitted activations impose at most e
additional affine tests in the product coordinates. The omitted
turns affect the affine origins and leave the two core spans
and the linear matching map fixed. Literal phases of both kept
and omitted turns are quadratic modulo four in the combined
variables, by ordered Pauli multiplication. The omitted
normalization is $`2^{-h/2}`$ for h their total activation count,
so it is at most one and is known from the parameter.

Apply Sections 1–4 to this entire parameterized array at
tolerance $`\delta 2^{-2e}`$. The match and exterior cutoff
bounds follow from the same product mass estimate, since the
extra tests only restrict the support. At a significant branch,
dropping the extra tests gives the product affinity needed to
find one relative linear map description for all branches. The
constraint Fourier average in Section 4 enforces the extra
tests. There are m+2e parameter bits and $`O(e)`$ extra bits
in the logarithmic accuracy budget. Summing the decoded branches
with their literal weights gives coefficient error at most delta.
This proves the stated net without recording path collisions
of the full halves as independent physical data.

For complete returned frames, use the actual reverse if needed
to put the larger core rank at the measurement side. The output
probe adjustment in Section 5 applies to the entire independent
last band's list. Thus for polynomial clean count, a four-band
word for a typical deepest-corner frame at the fixed error
threshold requires pairing deficiency outside the polynomial
window. The net also rules out $`e=2^{o(a)}`$ at that accuracy,
since its logarithmic cardinality is $`2^{o(a)}`$ on the tested
slice. Arbitrary four-band words can have much larger e. The
unrestricted asymptotic depth gap is open.
