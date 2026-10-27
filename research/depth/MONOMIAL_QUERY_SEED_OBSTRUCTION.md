# Commuting original and transported query algebras obstruct copying

[Encoded copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Field shear](FIELD_SHEAR_INTERTWINER.md)

The encoded-copy interface allows arbitrary unitary seeds and arbitrary
non-Clifford copy interactions. Even with that freedom, a seed monomial
in a simultaneous query eigenbasis has a complete copy-error floor of
one half when the address has at least two bits. The proof uses an
alternating four-word identity; it does not require a Hermitian seed,
a Clifford eigenbasis, or a restriction on the copy circuit.

This is a seed-class restriction, not a general native-depth bound.
A nonmonomial seed can escape the identity. R96 is unchanged.

## 1. The exact four-word identity

Let x range over the binary p-dimensional vector space. On the whole
arbitrary carrier and its helpers let the query obey the literal laws

```math
Q_xQ_y=Q_{x+y},\qquad Q_0=I.
```

Let R be any unitary and put

```math
R_x=Q_xRQ_x^\dagger,\qquad
B_x=RQ_xR^\dagger,\qquad D_x=Q_xB_x^\dagger.
```

Assume the original and transported query algebras commute:

```math
[Q_x,B_y]=0\quad\hbox{for all }x,y.                 \qquad(1)
```

Both Q and B are additive commuting representations. Equation (1)
therefore gives an additive commuting representation D, and
$`R_x=D_xR`$. For every origin c and directions a,b,

```math
R_{c+a}R_c^\dagger R_{c+b}R_{c+a+b}^\dagger=I.       \qquad(2)
```

Indeed the two adjacent R/actual-adjoint pairs cancel, leaving
$`D_{c+a}D_c^\dagger D_{c+b}D_{c+a+b}^\dagger=I`$.
This ordering matters: an arbitrary seed need not square to identity.

In particular, (1) holds whenever R is monomial in a common eigenbasis
of all Q's. Conjugating a diagonal matrix by a monomial unitary gives
a diagonal matrix. The seed's permutation and column phases can be
arbitrary, and the eigenbasis may include every helper. No native
implementation of that basis is being supplied. Condition (1) is
slightly broader: it can also hold for a nonmonomial seed acting
inside degenerate joint eigenspaces.

For the field query, Hadamards on its v bank give the explicit diagonal
family with signs $`(-1)^{v\cdot M_xu}`$. Thus any seed monomial in
those coordinates falls under this theorem, even if it has an
arbitrarily complicated reversible permutation or phase function.
This differs from the existing exclusion for seeds diagonal in a
query-invariant stabilizer basis.

## 2. A uniform floor for every physical copy interaction

Let Y be an arbitrary physical output bank with one qubit per address,
and let C be any fixed unitary on carrier, Y and arbitrary helpers,
independent of the external address. Define

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                    -R_x\otimes X_{Y_x}\|.          \qquad(3)
```

For p at least two, choose independent a,b. The four points in (2)
are distinct. Multiplying their four approximate copy identities,
with the actual adjoints in the second and fourth positions, costs
at most $`4\eta`$ by telescoping. The actual product is I by (2);
the ideal product is identity on the carrier and X on exactly those
four output wires. Their operator distance is two. Hence

```math
\eta\ge\tfrac12.                                   \qquad(4)
```

The proof allows C to mix every carrier, output and helper wire. It
uses the complete operator comparison on arbitrary input sectors;
identity helpers cannot alter the four-word identity. No claim of
sharpness of (4) is made.

The corresponding indicator echo
$`C\mathcal A C^\dagger\mathcal A^\dagger`$, where
$`\mathcal A=\sum_x|x\rangle\langle x|\otimes R_x`$,
has exactly the same defect as (3), by right multiplication with its
actual source adjoint and the direct-sum norm. Thus this interface
cannot hide the floor by completing the echo.

## 3. Robustness and the escape condition

Suppose an actual seed $`\widetilde R`$ has operator distance epsilon
from a seed R satisfying (1), with the same exact queries. Its orbit
has the same distance bound at every address. For every C, the copy
defects of the two orbits differ by at most $`2\epsilon`$, so

```math
\widetilde\eta+2\epsilon\ge\tfrac12.                \qquad(5)
```

This is a distance-to-class statement, not a claim that a small
commutator implies closeness to a monomial seed. If the actual query
at each address additionally differs from Q by at most delta, orbit
distance is at most $`\epsilon+2\delta`$ and the same argument gives
$`\widetilde\eta+2\epsilon+4\delta\ge1/2`$.

### A quantitative cross-commutator requirement

There is a direct screen for an arbitrary unitary seed, without a
monomial approximation or the cross-commutation assumption (1). Fix
independent a,b and let W be the alternating word (2) with c=0.
Because Q and B are representations by Hermitian involutions,

```math
W=Q_a(B_aQ_bB_aQ_b)Q_a,\qquad
\kappa_{a,b}=\|W-I\|=\|[B_a,Q_b]\|.               \qquad(6)
```

For the norm equality, right multiplication of
$`[B_a,Q_b]=B_aQ_b-Q_bB_a`$ by $`B_aQ_b`$ gives
$`B_aQ_bB_aQ_b-I`$. All multipliers are unitary; R itself need
not be Hermitian.

Write $`\operatorname{Re}W=(W+W^\dagger)/2`$. Unitarity gives
$`(W-I)^\dagger(W-I)=2I-2\operatorname{Re}W`$, hence

```math
\operatorname{Re}W\ge(1-\kappa_{a,b}^2/2)I.
```

The actual four-word copy product is $`C(W\otimes I_Y)C^\dagger`$.
It has the same lower bound on its Hermitian part. The target is
$`W\otimes X^{\rm plane}`$; on the nonempty negative eigenspace
of the four-output X word its Hermitian part is at most
$`-(1-\kappa_{a,b}^2/2)I`$. Testing the difference on any unit
vector in that sector proves distance at least
$`2-\kappa_{a,b}^2`$ whenever that number is positive. Four-factor
telescoping therefore yields

```math
\eta\ge\max\{0,(2-\kappa_{a,b}^2)/4\}.
                                                               \qquad(7)
```

Consequently every independent pair must satisfy

```math
\eta<1/2\quad\Longrightarrow\quad
\|[RQ_aR^\dagger,Q_b]\|\ge\sqrt{2-4\eta}.        \qquad(8)
```

In particular exact copying requires each of these norms to be at
least square root two. The condition is necessary, not sufficient:
passing this two-query screen supplies neither the entire copy
interaction nor a native implementation. The argument uses the
complete operator norm, not an average commutator measure.

A successful source must thus mix the original and transported query
algebras quantitatively, while retaining a complete, cheaply synthesized
copy. Merely adding nonlinear phases to a reversible seed in the
query eigenbasis cannot do so. A finite nonmonomial word escaping (2)
is included in the fixtures; escape is not a complete construction.

## 4. Exact finite evidence

`tests/test_monomial_query_seed_obstruction.py` checks arbitrary
non-Hermitian monomial phases, shifted address planes, the physical
field-query Fourier basis, arbitrary helper mixing in the seed, and
a nonmonomial escape example. Rational and native quarter-angle seeds
check the exact cross-commutator identity, its positive-error screen
and its square-root-two boundary. The output-word witness has a literal
minus-one eigenvector, giving distance two without a numerical
singular-value estimate. These tests support the identities; the
all-size result follows from the argument above.
