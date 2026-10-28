# Clifford seed copying needs noncommuting symplectic dynamics

[Encoded copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Clifford query rank](CLIFFORD_ORBIT_QUERY_RANK.md) ·
[Exact finite Clifford pilot](FIELD_COPY_CLIFFORD_PILOT.md)

A noncommuting Clifford seed orbit can admit a complete physical copy;
the linked four-address pilot supplies one. But noncommutation only by
Pauli factors cannot persist as a successful all-size field construction.
For seven or more address bits, commuting symplectic actions of the seed
orbit force complete copy error at least one sixty-fourth. The copy may
be an arbitrary unitary, with unrestricted depth and arbitrary helpers.

This is a restriction on a seed class. It does not exclude all Clifford
seeds, all native copies, or the complete indicator by another source.
It gives no improved R96 bound or growing native-depth lower bound.

## 1. Hypotheses and uniform error floor

Let $`x\in\mathbb F_2^p`$. On the full arbitrary carrier and every
helper, let Q be a literal additive Clifford representation:

```math
Q_xQ_y=Q_{x+y},\qquad Q_0=I.
```

Assume its binary symplectic matrix $`F(x)`$ has entrywise Boolean
degree at most d, where $`d\ge1`$. The one-query field shear has d equal
to one, including identity helper blocks. Let R be a fixed Clifford
whose projective order divides two; a Hermitian Clifford involution
is the principal case. Set

```math
R_x=Q_xRQ_x^\dagger,\qquad A(x)=S(R_x).
```

Suppose all matrices $`A(x)`$ commute. Equivalently, each actual
operator commutator of two orbit elements is a Pauli up to scalar.
It need not be a scalar, and the $`R_x`$ need not commute even
projectively. This is strictly broader than the commuting projective
orbit in the linked query-rank theorem.

For an arbitrary physical output bank Y of $`2^p`$ qubits and an
arbitrary fixed unitary C on carrier, outputs and helpers, define

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                    -R_x\otimes X_{Y_x}\|.
```

**Theorem.** Under the stated hypotheses,

```math
p\ge5d+2\quad\Longrightarrow\quad
\eta\ge2^{-5d-1}.
```

In particular one field shear has

```math
p\ge7\quad\Longrightarrow\quad\eta\ge1/64.        \qquad(1)
```

The complete indicator echo has exactly this copy defect, as proved
in the encoded-copy interface. Completing the actual-adjoint echo
therefore does not remove (1). No sharpness of these constants is claimed.

## 2. A large address cube has a Pauli product

Because Q is an additive representation over a binary group,
$`F(x)^{-1}=F(x)`$. Hence

```math
A(x)=F(x)S(R)F(x),\qquad \deg A\le2d.
```

Projective order two makes every $`A(x)`$ involutive. Apply the
[commuting square-zero factorization](CLIFFORD_ORBIT_QUERY_RANK.md#multiplicative-degree-rank)
to these commuting involutions. It gives fixed commuting involutions
$`H_T`$ such that

```math
A(x)=\prod_{|T|\le4d}H_T^{x^T}.                     \qquad(2)
```

The factorization takes place in the ambient binary matrix group;
the individual factors need not be symplectic. Only the identity
for the original symplectic values is used here. In particular the
product of $`A(x)`$ over any affine cube of dimension $`4d+1`$ is I:
each degree-at-most-$`4d`$ Boolean coefficient sums to zero on that cube.
Substitution of affine coordinates cannot increase Boolean degree.

Choose independent directions $`a_1,\ldots,a_{4d+1}`$ and put
$`U=\langle a_1,\ldots,a_{4d+1}\rangle`$. Fix an order on U and define

```math
W_y=\prod_{u\in U}^{\rm fixed\ order}R_{y+u}.
```

Equation (2) makes $`S(W_y)=I`$. The kernel of the Clifford
symplectic action is the Pauli group including scalar phases, so
each $`W_y`$ is a Pauli up to scalar. No Hermiticity of this ordered
word is assumed. Exact covariance also gives

```math
W_y=Q_yW_0Q_y^\dagger.                              \qquad(3)
```

If v is the binary Pauli label of $`W_0`$, the label of $`W_y`$ is
$`F(y)v`$, of degree at most d. This step does not discard a Pauli
correction or assume that phases have low address degree.

## 3. A second cube gives a scalar relation with nonzero copy parity

Choose another $`d+1`$ independent directions, independent also of U,
and let V be their span. Use any fixed order on V. The Pauli label of

```math
\Omega=\prod_{y\in V}^{\rm fixed\ order}W_y
```

is the binary sum of $`F(y)v`$ over V, which vanishes by its degree
bound. Thus $`\Omega=\lambda I`$ for a scalar of modulus one.
Its literal phase need not be computed and need not be one.

Expanded as a word in the $`R_x`$, this product has length

```math
\ell=2^{4d+1}2^{d+1}=2^{5d+2}.
```

The direct-sum condition $`U\cap V=\{0\}`$ makes all its address
labels distinct. It contains every point of $`U+V`$ exactly once.
Its occurrence parity is therefore nonzero: the desired copied word
has the physical output factor

```math
X_{U+V}=\prod_{x\in U+V}X_{Y_x}\ne I.
```

Multiplying the approximate copy identities in exactly the same order
costs at most $`\ell\eta`$ by unitary telescoping. The actual product
is $`\lambda I`$, while the requested product is
$`\lambda I\otimes X_{U+V}`$. Their operator distance is two, since
the nonidentity output X word has a negative-one eigenspace. Hence
$`2\le\ell\eta`$, proving the theorem.

For a field shear this is a word of 128 distinct orbit elements:
four translated five-dimensional cubes. Arbitrary identity helpers
do not alter the symplectic degree, the scalar relation or the output
parity. The proof includes every carrier, output and reference input.

## 4. Perturbations and what a scalable construction must change

Let an actual seed be within epsilon of a seed covered by the theorem,
and let every actual query be within delta of the exact Q above, in
complete operator norm. Each actual orbit element is within
$`\epsilon+2\delta`$ of its reference. For any same copy C the two
copy defects consequently differ by at most
$`2\epsilon+4\delta`$. In the field case,

```math
\widetilde\eta+2\epsilon+4\delta\ge1/64
\qquad(p\ge7).                                      \qquad(4)
```

This is a distance-to-class statement. It does not infer a nearby
commuting symplectic orbit from small commutators, nor exclude
arbitrary non-Clifford seeds.

The four-address Clifford pilot has nonscalar central Pauli
commutators, and therefore lies inside the theorem's seed class.
It demonstrates a genuine finite copy mechanism, while (1) prevents
an all-size extension retaining that structure. A Clifford-seed route
at growing address size and arbitrarily small fixed error must allow
noncommuting symplectic actions, or change the query family beyond the
bounded-degree hypothesis. A non-Clifford seed is another possibility.
Each alternative still needs a paid complete copy; crossing this screen
alone is not a construction or an improved frame-depth theorem.

The [exact fixtures](../../tests/test_clifford_seed_copy_structure.py)
check the necessary degree inflation in (2), literal scalar phases of
translated Pauli words, and a sparse full-operator 128-address cube
relation on a field-shear orbit. Finite checks support the identities;
Sections 1--3 supply the all-size result.
