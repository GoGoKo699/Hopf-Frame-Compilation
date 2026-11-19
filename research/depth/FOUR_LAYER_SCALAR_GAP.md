# Scalar and rank gaps at four native T layers

For a physical circuit $`V`$ of at most four native T layers, two Hermitian Paulis $`P,Q`$, and
physical dimension $`D`$, let $`e=(1-D^{-1}\mathrm{Tr}(QVPV^\dagger))/2`$. For data-bit
observables this is the mean disagreement probability on uniformly sampled computational
columns, with no clean initialization. Then

```math
e=0\quad\text{or}\quad e\ge\frac{1-\cos(\pi/12)}8.
```

There is a stronger conclusion when the exceptional columns are known exactly. If the two
observables agree on a subspace of relative dimension greater than $`1-\sin^2(\pi/24)/2`$, they
agree everywhere. In particular, agreement off one $`k`$-bit address forces agreement there as
well whenever $`k\ge7`$. The circuit may move the address and use arbitrary dirty work.

## Sparse Pauli probes

In the
[R75 midpoint frame](FOUR_LAYER_OPERATOR_RIGIDITY.md#1-two-layer-pauli-images-and-commuting-axis-spaces),
write $`A=E_AK_AE_A^\dagger`$ and $`B=E_BK_BE_B^\dagger`$ for the two transported observables.
Here each $`E_j`$ is a Clifford conjugate of one physical T band, and each $`K_j`$ is Hermitian
Clifford. Their full-space normalized Hilbert–Schmidt distance is $`2\sqrt e`$. The binary Pauli
axes of $`E_A`$ are independent because its original T gates act on distinct wires. The map
pairing a Pauli with these axes is onto. Choose a basis of its kernel and a lift of each
coordinate vector of the image. This gives a generating set of Hermitian Pauli probes $`H`$,
each anticommuting with at most one axis.

If $`H`$ meets one axis, $`E_A^\dagger H E_A`$ expands into two distinct signed Paulis with
coefficients $`1/\sqrt2`$. Conjugation by $`K_A`$ preserves Paulis, and conjugation by $`E_A`$
takes each Pauli to a Hermitian Clifford. The normalized Choi vector of $`AHA`$ is therefore
$`(U+V')/\sqrt2`$ for orthogonal unit stabilizer vectors $`U,V'`$ (absorb the signs). With no
anticommuting axis it is one stabilizer vector. On the other side $`BHB=E_B C_HE_B^\dagger`$
with $`C_H`$ Hermitian Clifford. In an $`E_B`$-diagonalizing Clifford basis, its Choi vector
$`C`$ is a stabilizer vector multiplied by the diagonal phase of $`E_B\otimes\overline E_B`$.

Each of these states is flat on its affine support. Ordinary stabilizer phases there are a fixed
eighth root times $`i`$ to a sum of affine bits times a quadratic sign. The eighth-root global
phase follows from Hermiticity by comparing a nonzero matrix entry with its transposed
conjugate, as in
[R75 §2](FOUR_LAYER_OPERATOR_RIGIDITY.md#2-the-one-phase-gap-with-arbitrary-hermitian-clifford-boundaries).
The extra diagonal phase is an eighth root raised to a signed sum of affine parity bits on the
support.

## A gap for the Choi vectors

Put $`b=2\sin(\pi/24)`$. We show that the two Choi vectors coincide or have distance at least
$`b`$. The single-stabilizer case follows directly from R75 §2. In the two-summand case, if
either overlap with $`C`$ has modulus at most $`1/2`$, orthogonality and Bessel's inequality
bound the real overlap of $`C`$ with the sum by $`(1+\sqrt3)/(2\sqrt2)=\cos(\pi/12)`$. This
gives the claimed distance.

Suppose both overlaps exceed $`1/2`$. For intersecting affine supports $`S,T`$ define
$`d(S,T)=\dim S+\dim T-2\dim(S\cap T)`$; set $`d=\infty`$ for disjoint supports. A flat-state
overlap has modulus at most $`|S\cap T|/\sqrt{|S||T|}=2^{-d/2}`$. Relative to
$`S=\mathrm{supp}C`$, each other support is therefore equal to $`S`$, one index-two subspace of
it, or one index-two extension of it.

The following phase estimate applies on every affine support or half-support below. A relative
phase $`\omega^p`$, with $`\omega=e^{i\pi/4}`$, has exponent digits of binary degrees at most
$`1,2,4`$. Indeed, for a sum of affine bits, digit $`j`$ is the elementary symmetric polynomial
of order $`2^j`$ in those bits, by the parity of the corresponding binomial coefficient. Repeat
affine bits to represent signed integer coefficients modulo eight. A quadratic sign contributes
$`4q`$ with $`\deg q\le2`$. The additional exponent $`-2q`$ occurring below contributes
$`2q+4q`$ modulo eight. Its middle digit has degree at most two and the carry from that digit to
the high digit has degree at most four. If the phase is nonconstant or a nonidentity constant,
its lowest nonzero digit takes value one on at least $`1/2,1/4,`$ or $`1/16`$ of the points,
respectively, by the minimum weight of a nonzero binary polynomial. The corresponding real
cosine losses are at least $`1-1/\sqrt2,1,2`$. Thus the mean loss is zero or at least $`1/8`$.

* **A support extends $`S`$.** If the other support lies in $`S`$, the outside contribution to
  the norm is $`1/2`$. Distinct extensions have disjoint outside cosets. If the two extensions
  coincide, both summands on $`S`$ have half the target magnitude. Since both states being
  compared have norm one, their squared distance is the sum of the two mean phase losses on
  $`S`$. This is zero or at least $`1/8`$.
* **A support is half of $`S`$.** Two noncomplementary halves leave at least a quarter of $`S`$
  uncovered. On complementary halves the squared distance is the sum of the two mean phase
  losses. If the other support equals $`S`$, orthogonality on their shared half gives masses
  $`3/4,1/4`$ for the sum on the two blocks. The target masses are $`1/2,1/2`$. The triangle
  inequality on each block gives squared distance at least
  $`(\sqrt3/2-1/\sqrt2)^2+(1/2-1/\sqrt2)^2=b^2`$.
* **All supports equal $`S`$.** The ratio $`w=V'/U`$ is a fixed eighth root times a Clifford
  phase. For odd eighth-root exponents, $`|(1+w)/\sqrt2|`$ differs from one by at least
  $`\sqrt{1+1/\sqrt2}-1>b`$ at every point. Otherwise $`w`$ takes fourth-root values. Its
  real-valued positions are either empty or at least half the support, because the low bit of
  its mod-four exponent is affine. On real positions the magnitude differs from one by at least
  $`\sqrt2-1`$. If there are no real positions, write $`w=i(-1)^q`$ with $`q`$ binary quadratic.
  Then $`(U+V')/\sqrt2=U\omega i^{-q}`$ throughout, and the phase estimate applies.

## Mean error, exceptional rank and complete operators

Every probe has $`\|AHA-BHB\|_{2,\mathrm{tr}}\le2\|A-B\|_{2,\mathrm{tr}}=4\sqrt e`$. Below $`b`$
all probe conjugations agree, so $`BA`$ commutes with a full Pauli basis and is a scalar. Since
$`A,B`$ are Hermitian involutions this gives $`A=\pm B`$, and the minus sign has $`e=1`$. This
proves the mean bound.

For an exact exceptional subspace put $`W=BA`$. Telescoping $`WHW^\dagger-H`$ bounds its rank by
$`2\mathrm{rank}(W-I)`$ and its operator norm by two, so

```math
\|AHA-BHB\|_{2,\mathrm{tr}}^2
\le 8\mathrm{rank}(W-I)/D.
```

When the right side is below $`b^2`$, the same probes force $`W`$ scalar. The fixed subspace
then gives $`W=I`$. For $`k\ge7`$ it is at most $`1/16\lt b^2`$. Positive flip probability may
be arbitrarily small on the exceptional address. Exact silence on every other computational
column means the corresponding Pauli eigenspaces agree on those columns.

For $`p\ge7`$ the complete dirty indicator $`|x,Y,w\rangle\mapsto|x,Y\oplus e_x,w\rangle`$ and
the fully returned $`p`$-controlled NOT have four-layer operator error at least
$`\tau=\sin(\pi/8)/(2\sqrt2)`$, even up to a common phase. Hadamards on the target bank
transform the exact operator into a diagonal sign. Below $`\tau`$, R75 operator rigidity forces
commutation with each computational $`Z`$, hence the implementing circuit is diagonal. Its
unit-modulus diagonal entries are eighth roots in the native phase convention by
[R75 §4](FOUR_LAYER_OPERATOR_RIGIDITY.md#4-discrete-phases-in-the-common-commutant). After
multiplying by the ideal signs, all entries lie in one disk of diameter
$`2\tau\lt 2\sin(\pi/8)`$, hence coincide. The resulting exact oracle contradicts the
single-address rank gap.

An explicit real frame sets its deepest rotation angle to $`\pi`$ at one prefix and all other
angles to zero. Its sign on that prefix acts as identity on the leaf bit. Conjugating a prefix
bit by Hadamards gives a controlled flip, up to a sign on the selected branch, when the other
$`n-2`$ prefix bits match. The diagonal freezing argument followed by the exact single-address
test still applies to that sign. Thus $`n\ge9`$ gives the same zero-clean error floor for a
complete real frame and for its complex realization. The
[canonical frame domains](../../docs/HOPF_INTERFACE.md#canonical-domains) admit these angles and
phases.
