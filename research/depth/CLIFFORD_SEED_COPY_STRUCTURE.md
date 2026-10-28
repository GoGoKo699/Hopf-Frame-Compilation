# Clifford-seed copy obstructions from orbit structure and native depth

<a id="clifford-seed-copying-needs-noncommuting-symplectic-dynamics"></a>

[Encoded copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Clifford query rank](CLIFFORD_ORBIT_QUERY_RANK.md) ·
[Exact finite Clifford pilot](FIELD_COPY_CLIFFORD_PILOT.md)

A noncommuting Clifford seed orbit can admit a complete physical copy;
the linked four-address pilot supplies one. But noncommutation only by
Pauli factors cannot persist as a successful all-size field construction.
For seven or more address bits, a projectively involutive seed with commuting
symplectic orbit actions has copy error at least one sixty-fourth. The copy may
be an arbitrary unitary, with unrestricted depth and arbitrary helpers.

There is a separate native-depth restriction: for nine or more address
bits, **every Clifford seed**, of any projective order, has copy error at
least $`\sin(\pi/8)`$ when the copy has one physical T/T-dagger band.
This second theorem permits noncommuting symplectic dynamics. Its proof
uses robust one-band matching and a literal phase-degree bound.

These are restrictions on specified seed/copy pairs. They do not exclude
all Clifford seeds with deeper copies, non-Clifford seeds, or the complete
indicator by another source. They give no improved R96 bound or growing
native-depth lower bound.

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
Section 5 further excludes a one-band copy for every Clifford seed at
large address size, even when its symplectic dynamics do not commute.

The [exact fixtures](../../tests/test_clifford_seed_copy_structure.py)
check the necessary degree inflation in (2), literal scalar phases of
translated Pauli words, and a sparse full-operator 128-address cube
relation on a field-shear orbit. Finite checks support the identities;
Sections 1--3 supply the all-size result.

<a id="arbitrary-clifford-seed-one-band-copy"></a>
<a id="one-band-clifford-seed"></a>

## 5. Every Clifford seed fails with a one-band copy at large address size

Now specialize Q to the literal field CNOT shear, identity on every
helper. Let R be **any Clifford**, with no involution or projective-order
hypothesis, and let the proposed copy have the native form

```math
C=F D G,\qquad D=\prod_{i\in B}T_i^{\epsilon_i},
\qquad \epsilon_i\in\{+1,-1\}.
```

F and G may jointly mix the carrier, physical output and arbitrary
helpers. The active T coordinates are distinct physical wires. The
complete copy error eta is defined exactly as in Section 1.

**One-band theorem.** In this class,

```math
p\ge9\quad\Longrightarrow\quad\eta\ge\sin(\pi/8). \qquad(5)
```

The statement allows unrestricted count and width and arbitrary
reference entanglement. In particular, it does not assume that the seed
orbit or an output-dual orbit has commuting symplectic actions.
The threshold nine is conservative; no sharpness is asserted.

### Robust exactness for this native conjugation

For arbitrary Cliffords K,L and a physical band D,

```math
\|D K D^\dagger-L\|\lt\sin(\pi/8)
\quad\Longrightarrow\quad D K D^\dagger\doteq L.   \qquad(6)
```

Here $`\doteq`$ means equality up to a scalar. This strengthens the
retained [general third-level rounding radius](NATIVE_CLIFFORD_TRANSPORT.md#1-a-robust-exactness-lemma-for-every-third-level-source)
for this particular physical-band conjugation. It does not assert the
stronger radius for every third-level source.

Let B also denote the binary span of the active Z axes. A Hermitian
Pauli P with label in $`B^\perp`$ commutes with D, so its image under
$`DKD^\dagger`$ is $`D(KPK^\dagger)D^\dagger`$. If
$`KPK^\dagger`$ anticommutes with any active axis, that image has
at least two orthogonal Pauli coefficients, all of magnitude at most
$`1/\sqrt2`$. Its normalized Hilbert--Schmidt distance from every
signed Hermitian Pauli is at least

```math
\gamma=\sqrt{2-\sqrt2}=2\sin(\pi/8).
```

The assumed operator distance in (6) would put this image within
$`2\|DKD^\dagger-L\|\lt\gamma`$ of the Pauli $`LPL^\dagger`$.
Thus K preserves $`B^\perp`$, and consequently B.

It follows that $`E=DKD^\dagger K^\dagger`$ is a diagonal function
of the active Z axes. It is the identity on their commuting full matrix
factor, and its eigenvalue ratios are eighth roots of unity. Set
$`L'=LK^\dagger`$. For every Pauli in $`B^\perp`$, conjugation
by L' differs from that by E by less than $`2\sin(\pi/8)\lt\sqrt2`$.
Distinct signed Hermitian Paulis have operator distance at least
$`\sqrt2`$. Hence L' fixes every such Pauli exactly. The Pauli
commutant makes L' a diagonal function of the active Z axes as well.
Being Clifford, its eigenvalue ratios are fourth roots of unity.

The ratios of corresponding eigenvalues of E and L' therefore lie in
one scalar multiple of the eighth roots. Any two distinct values are
separated by at least gamma. The operator-norm estimate and the triangle
inequality put every such pair at distance strictly less than gamma.
All ratios coincide, proving (6), with no chosen input state or
initialized subspace.

The retained [exact axis-matching theorem](ONE_BAND_AXIS_MATCHING.md#1-statement-with-signs-inactive-wires-and-error)
now implies that K signed-permutes the individual active Z axes.
Preservation of their span alone would not imply this conclusion.
This use of matching includes empty bands, inactive wires and either T
orientation; the matching proof remains in its original home.

### Literal phase degree of a field-conjugated Clifford

Use the ordered Pauli section $`P(a,b)=X^aZ^b`$. For any fixed
Clifford U,

```math
U P(v)U^\dagger=i^{\phi_U(v)}P(S_Uv),
```

where $`\phi_U`$ is a polynomial modulo four with linear terms
and even-coefficient quadratic terms. Multiply the images of the
ordered X and Z generators to obtain this formula; their pairwise
commutation phases supply exactly the quadratic terms.

If the coordinates v depend on x with Boolean degree at most d, this
phase has a multilinear polynomial modulo four of degree at most 2d.
Indeed, a Boolean sum of monomials lifts modulo four to their ordinary
sum minus twice the sum of pairwise products. Terms with at least three
monomial factors vanish modulo four. A quadratic term already multiplied
by two needs only its Boolean product. These are exact finite-ring
identities, independent of the number of coordinates or monomials.

The field CNOT query is phase-free in this ordered section:

```math
Q_x P(v)Q_x^\dagger=P(F(x)v),\qquad\deg F\le1.
```

Consequently, for any fixed Pauli input, the Clifford

```math
K_x=G(R_x\otimes I_Y)G^\dagger
```

has output-label degree at most two and literal phase degree at most
four modulo four. Read the native conjugation from the right: G-dagger
contributes a fixed phase; the first Q makes the label affine; R
contributes phase degree at most two; the second Q makes the label
quadratic; G contributes phase degree at most four.

### The physical output sign has degree at most eight

Suppose, for contradiction, that $`\eta\lt\sin(\pi/8)`$. Apply
(6) to K_x and the Clifford target
$`F^\dagger(R_x\otimes X_{Y_x})F`$ for every address. Exact
matching gives permutations and literal signs

```math
K_x Z_i K_x^\dagger=(-1)^{s_i(x)}Z_{\pi_x(i)}.
```

The permutation entries $`P_{ji}(x)=[\pi_x(i)=j]`$ have degree
at most two, since they are entries of the binary symplectic action of
K_x on the independent active axes. Each sign s_i has Boolean degree
at most four. To justify division by two, the corresponding phase
polynomial modulo four has only even values. Boolean Möbius inversion
then gives only even multilinear coefficients. Dividing those
coefficients by two gives the asserted Boolean sign polynomial.

Write $`\epsilon_i=(-1)^{e_i}`$. The exact diagonal correction
is, up to its scalar,

```math
\begin{aligned}
E_x&=D K_x D^\dagger K_x^\dagger\doteq\prod_j S_j^{k_j(x)},\\
k_{\pi_x(i)}
 &=\frac{\epsilon_{\pi_x(i)}-\epsilon_i(-1)^{s_i(x)}}2\\
 &=e_i-e_{\pi_x(i)}+\epsilon_i s_i(x).              \qquad(7)
\end{aligned}
```

The last line is an integer identity, with k in $`\{-1,0,1\}`$;
it avoids an unnecessary product of a permutation entry and a lifted
sign. The omitted scalar is irrelevant only to this adjoint action.

Fix a physical output probe $`Z_{Y_t}`$, and let
$`P_0=F^\dagger Z_{Y_t}F`$, retaining its fixed Hermitian scalar.
Write a_i for its active X coordinates. If the active X coordinates
after conjugation by K_x are a'_j, symplectic pairing gives
$`a'_{\pi_x(i)}=a_i`$. The S correction in (7) therefore contributes
phase

```math
\sum_i a_i\bigl(e_i-e_{\pi_x(i)}+\epsilon_i s_i(x)\bigr)
\pmod4.                                                    \qquad(8)
```

The binary function $`e_{\pi_x(i)}=\sum_j e_jP_{ji}(x)`$ has
degree at most two, so its modulo-four lift has degree at most four.
The modulo-four lift of s_i has degree at most eight. Thus (8) has
degree at most eight. K_x's own phase on P_0 has degree at most four.

Projective equality to the requested copy implies that the label after
E_x K_x is exactly the fixed label of P_0: the final F returns it to
$`Z_{Y_t}`$. The final F phase is therefore constant. Hence

```math
(CR_xC^\dagger)Z_{Y_t}(CR_x^\dagger C^\dagger)
       =(-1)^{[x=t]}Z_{Y_t}
```

has a literal phase polynomial modulo four of degree at most eight.
Its values are even, so the same coefficientwise division proves that
$`[x=t]`$ has Boolean degree at most eight. But a singleton on p
address bits has degree exactly p:

```math
[x=t]=\prod_{i=1}^p(1+x_i+t_i).
```

This contradicts $`p\ge9`$ and proves (5).

### Scope, stability and exact fixtures

A seed $`R=VPV^\dagger`$ with V containing at most one physical
T band is Clifford. Thus, for a one-band copy at arbitrarily small
fixed error, the two-band V grammar is the first preparation depth
not excluded by this theorem. Two-band V can produce a non-Clifford
seed and remains outside the proof. A deeper copy is another possible
change; its native cost and complete return must still be paid.

If a seed is within epsilon of this Clifford class and each actual
query is within delta of the exact field shear, the same one-band
copy obeys the distance-to-class bound

```math
\widetilde\eta+2\epsilon+4\delta\ge\sin(\pi/8),
\qquad p\ge9.                                         \qquad(9)
```

This follows by the full-operator triangle inequality. It does not
infer nearness to the class from small commutators or a finite fit.

The [exact fixtures](../../tests/test_one_band_clifford_seed.py)
check signed affine phase intersections, including an orthogonal
four-axis map that fails the cubic condition; literal native matching
phases and actual inverses with an arbitrary inactive helper; and a
field-shear example whose phase genuinely has degree four modulo four.
They also verify the required modulo-four lifts and the singleton
contradiction. These bounded checks support the displayed identities;
the analytic proof supplies the all-size error floor.
