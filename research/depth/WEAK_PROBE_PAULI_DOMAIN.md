# Pauli-domain bounds for a weak complete-probe decoder

[Full-decoder restriction](FIELD_DECODER_PAULI_DOMAIN.md) ·
[Four-address native cell](FLAT_FIELD_BANK_PILOT.md) ·
[Recursive field query](RECURSIVE_FIELD_SHEAR_DECODER.md)

A decoder that transports just one physical probe already has a restricted
Pauli domain. For affine Clifford queries, its image-domain bank X masks
must be affine functions of the address below complete probe error
one-quarter. A full-rank field query gives a stronger restriction on pure
bank Z masks. Together these facts give independent copies of the retained
four-wire cell an error floor of one-twelfth from three address bits onward,
even after arbitrary joint input and output Clifford walls.

These are Pauli-domain, axis-rank and count restrictions. They do not require
the complete query endpoint to be a wire permutation, and do not give a
growing native T-depth lower bound or change R96. Coupling the non-Clifford
operations remains outside the independent-cell screen.

<a id="complete-probe-contract"></a>

## 1. The weaker complete-space contract

Let $`x\in\mathbb F_2^p`$, $`m=2^p`$, and $`n=m+q`$. The target consists
of an m-qubit bank indexed by x and q other arbitrary qubits. Write
$`Z_x=Z_{Y_x}\otimes I_q`$. Let $`Q_x`$ be a Clifford family on all n
physical qubits, and let its binary symplectic matrix $`F_x`$ have entrywise
Boolean degree at most d. This first theorem needs neither additivity nor
an identity action on any source helpers.

For an arbitrary unitary E put

```math
R=E^\dagger Z_0E,\qquad R_x=Q_xRQ_x^\dagger,\qquad
\epsilon=\max_x\|ER_xE^\dagger-Z_x\|_{\rm op}.
```

This is a complete operator comparison, including every dirty input and
reference entanglement. It constrains the transported probe only; it does
not require $`EQ_xE^\dagger`$ to be Clifford or to equal a specified
wire-translation query on the whole space.

Use the [binary Pauli domain](BOUNDED_DEFICIENCY_TRANSPORT.md#2-signed-localization-of-a-possibly-degenerate-domain)
and its conjugated image:

```math
H_E=\{v:EP(v)E^\dagger\text{ is a signed Hermitian Pauli}\},
\qquad \mathcal U_E=\text{the image-label space of }H_E.
```

Both are linear spaces, and conjugation gives an injective symplectic map
between them. Let $`\mathrm{RM}(d,p)`$ be the Boolean functions of degree
at most d, evaluated on all addresses, and set

```math
M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj.
```

**Degree-sensitive theorem.** If $`\epsilon\lt2^{-d-1}`$, every bank X
mask in $`\mathcal U_E`$ belongs to $`\mathrm{RM}(d,p)`$. Consequently

```math
\dim H_E\le m+M_d(p)+2q,\qquad
\mathrm{rank}_{\rm symp}H_E\le2M_d(p)+2q.
```

If A is the binary span of the absorbed T axes in any native word for E,
then

```math
\dim A\ge m-M_d(p),\qquad T(E)\ge m-M_d(p).
```

For affine queries, $`d=1`$, $`M_1(p)=p+1`$, and the threshold is
one-quarter. No degree bound on the Clifford conjugation signs is assumed.

## 2. Scalar Pauli cubes retain every phase

Fix a Hermitian Pauli $`P\in H_E`$, put $`S=EPE^\dagger`$, and denote
its bank X mask by $`a(x)`$. Since
$`Z_xSZ_x=(-1)^{a(x)}S`$, two-factor telescoping gives

```math
\|R P_xR-(-1)^{a(x)}P_x\|\le2\epsilon,
\qquad P_x=Q_x^\dagger P Q_x.
```

The inverse-query orientation matters. Its label is $`F_x^{-1}v`$,
whose degree is still at most d because
$`F_x^{-1}=JF_x^TJ`$ for the fixed binary symplectic form J.

For $`d\lt p`$, fix any affine cube with $`d+1`$ independent directions
and an order of its $`2^{d+1}`$ vertices. The sum of the $`P_x`$ labels
on the cube is zero. Their ordered product is therefore

```math
\Omega=\prod_{x\text{ in the cube}}P_x=\lambda I,
\qquad \lambda\in\{1,-1,i,-i\}.
```

The factors need not commute, and the scalar need not be positive.
Conjugation by R fixes this scalar product exactly. Applying the preceding
approximate conjugations in the same order costs at most
$`2^{d+2}\epsilon`$. An odd cube parity of a would change the desired
product to $`-\Omega`$, at distance two. Hence

```math
\text{odd cube parity of }a
\quad\Longrightarrow\quad
2\le2^{d+2}\epsilon.
```

Below the theorem's threshold, every such derivative vanishes, so a is in
$`\mathrm{RM}(d,p)`$. If $`d\ge p`$, that code is the whole bank mask
space and the conclusion is automatic. This proof uses the degree of the
labels only; arbitrary literal Pauli signs remain inside $`\Omega`$.

## 3. Domain dimension, symplectic rank and absorbed axes

The image domain lies in

```math
V_{d,p,q}=\{(a,b,c):a\in\mathrm{RM}(d,p),\quad
 b\in\mathbb F_2^m,\quad c\in\mathbb F_2^{2q}\}.
```

Its dimension is $`m+M_d(p)+2q`$. Its radical consists of pure bank Z
masks orthogonal to every Reed--Muller X mask; it has dimension
$`m-M_d(p)`$. Thus the alternating form has rank $`2M_d(p)+2q`$.
Restriction to a subspace cannot increase that rank, proving both bounds.

Every Pauli commuting with all absorbed rotation axes stays Pauli under E,
so $`A^\perp\subseteq H_E`$. Nondegeneracy of the n-qubit pairing gives
$`2n-\dim A\le\dim H_E`$. The T count is at least the span rank of
its individual axes, giving the asserted count bound.

<a id="exact-additive-kernel"></a>

## 4. Exact additive queries also restrict the pure-Z kernel

Assume now $`\epsilon=0`$ and that the queries are additive up to scalar
phase: $`Q_yQ_x`$ is proportional to $`Q_{x+y}`$. For a bank mask v define

```math
Z^v=\prod_x Z_x^{v_x},\qquad
K_Z=\{v:E^\dagger Z^vE\text{ is a Pauli}\}.
```

Let $`\Pi_y`$ translate bank masks by y. Exact orbit covariance gives

```math
Q_yE^\dagger Z^vEQ_y^\dagger=E^\dagger Z^{\Pi_yv}E.
```

Thus $`K_Z`$ is translation-invariant. The map $`\phi`$ taking a mask
to its pulled-back Pauli label is linear and injective: a nonzero positive
Z word cannot be conjugated to a scalar. At the label level,

```math
\phi(\Pi_yv)=F_y\phi(v).
```

For any $`d+1`$ independent directions, additivity and the degree bound give

```math
\prod_{j=1}^{d+1}(F_{y_j}-I)
=\sum_{x\in\mathrm{span}\{y_1,\ldots,y_{d+1}\}}F_x=0.
```

Injectivity transfers the same translation-difference identity to
$`K_Z`$, so $`K_Z\subseteq\mathrm{RM}(d,p)`$ and
$`\dim K_Z\le M_d(p)`$. Project $`\mathcal U_E`$ onto its bank X
mask and all helper labels. Its kernel is precisely the pure-bank Z
space $`K_Z`$, and its image has dimension at most $`M_d(p)+2q`$.
Therefore exact additive queries satisfy the stronger bounds

```math
\dim H_E\le2M_d(p)+2q,\qquad
\dim A\ge2m-2M_d(p),\qquad T(E)\ge2m-2M_d(p).
```

No approximate invariance of $`H_E`$ or $`K_Z`$ is asserted. The robust
field argument below explicitly pays the weight of the Z word being moved.

<a id="field-query-and-pauli-helpers"></a>

## 5. Full-rank field queries, including Pauli helper motion

Let $`p\ge2`$. On $`2k`$ core qubits take the additive field shear

```math
Q_x^{\rm fld}|u,v\rangle=|u,v+M_xu\rangle,\qquad
M_x=\sum_{j=0}^{p-1}x_jM_j,
```

where every nonzero $`M_x`$ is invertible. This permits any injective
address subspace of the field. On h further source qubits allow a Pauli
query $`B_x`$ whose binary label is linear in x, and put

```math
Q_x=Q_x^{\rm fld}\otimes B_x,\qquad n=2k+h=m+q.
```

Identity source helpers are the case $`B_x=I`$. Arbitrary scalar phases
of $`B_x`$ are irrelevant to conjugation; its linear labels make the
joint query a projective additive representation. The Pauli factors need
not commute literally. This includes the moving Pauli helper query in an
augmented recoupling construction.

The symplectic action is still $`F_x=I+N_x`$, with
$`N_xN_y=0`$. Every nonzero $`N_x`$ has the same kernel

```math
F_{\rm fix}=\mathrm{span}\{Z_u,X_v\}\oplus V_h.
```

This is a common fixed-**label** space. With Pauli helpers, a fixed label
can acquire an address-dependent sign and need not be fixed literally.
The proof below retains that distinction.

For any Hermitian Pauli P with label v, the source orbit
$`P_y=Q_yPQ_y^\dagger`$ is pairwise commuting. Indeed,
$`\langle v,N_yv\rangle=0`$, and the derivative labels lie in the
isotropic core space $`\mathrm{span}\{Z_u,X_v\}`$, with zero helper
component. Every four-address plane product
$`P_0P_yP_zP_{y+z}`$ consequently has zero label and is a Hermitian
scalar $`\pm I`$. The helper signs do not change these facts.

At exact error, Section 4 puts $`K_Z`$ in $`\mathrm{RM}(1,p)`$.
Every affine mask has a nonzero translating stabilizer because $`p\ge2`$.
For such a stabilizer y, $`\phi(v)`$ is fixed by $`F_y`$ and hence by
every $`F_x`$. Injectivity of $`\phi`$ forces v to be translation-fixed,
so

```math
K_Z\subseteq\mathrm{span}\{\boldsymbol1\},\qquad
\dim H_E\le p+2+2q,\qquad
\dim A\ge2m-p-2,\qquad T(E)\ge2m-p-2.
```

Here q counts all physical coordinates outside the target bank. It need not
equal h or describe query-fixed wires. The theorem does not prescribe any
decoded helper action and does not require E itself to return helpers.

<a id="robust-low-weight-kernel"></a>

## 6. A robust low-weight kernel lemma

Keep the field and Pauli-helper hypotheses of Section 5, but allow
nonzero $`\epsilon`$. Define $`K_Z`$ by exact Pauli pullback as before;
it is linear even when it is not translation-invariant.

**Low-weight lemma.** If $`0\ne v\in K_Z`$ has weight w, then

```math
\epsilon\lt\frac1{3w}\quad\Longrightarrow\quad
v=\boldsymbol1.
```

Put $`P=E^\dagger Z^vE`$ and form $`C_v`$ as a fixed ordered product
of the w operators $`R_x`$ at occupied addresses. The probe comparison
gives $`\|C_v-P\|\le w\epsilon`$. Query covariance translates that
same ordered product. Comparing its factors with the commuting target Z
words gives, for every y,

```math
\|Q_yPQ_y^\dagger-E^\dagger Z^{\Pi_yv}E\|
\le2w\epsilon.
```

No approximate factors are reordered. At $`y=0`$ the comparison is exact.
Compare the four translated Paulis on the plane through $`0,y,z,y+z`$
with the corresponding target Z words. The source product is $`\pm I`$;
the zero vertex has no error and the other three cost at most
$`2w\epsilon`$ each. If $`(I+\Pi_y)(I+\Pi_z)v`$ is nonzero, the
target product is a conjugate of a nontrivial Z word, at distance two from
either scalar. Thus $`6w\epsilon\lt2`$ forces every second difference
to vanish, making v affine.

Choose a nonzero translating stabilizer y of that affine mask. The same
comparison puts $`Q_yPQ_y^\dagger`$ within $`2w\epsilon\lt2`$ of P.
These commuting Hermitian Paulis have distance either zero or two, so
their labels agree. The common fixed-label kernel then gives
$`Q_xPQ_x^\dagger=\sigma_xP`$ for every x, with
$`\sigma_x\in\{1,-1\}`$. After conjugating the comparison by E,

```math
\|\sigma_x Z^v-Z^{\Pi_xv}\|\le2w\epsilon\lt2.
```

Two canonical positive Z words, with the indicated possible sign on the
first, have distance zero only when their masks agree and $`\sigma_x=1`$;
otherwise the distance is two. Hence every translation fixes v and every
sign is positive. This proves the lemma without incorrectly promoting a
fixed helper label to a literally fixed helper Pauli.

Since $`w\le m`$, the stronger field domain and axis bounds of Section 5
also hold under the width-dependent robust threshold

```math
\epsilon\lt\frac1{3m}.
```

This threshold is below one-quarter at $`p\ge2`$, so the bank X-mask
bound used in the domain projection applies simultaneously.

<a id="independent-cell-error-floor"></a>

## 7. Independent four-wire cells have a uniform error floor

The retained [six-T cell D](FLAT_FIELD_BANK_PILOT.md) has Pauli domain

```math
H_D=\mathrm{span}\{X_0,X_2,X_3,X_1Z_0\},\qquad
\dim H_D=4,\qquad\mathrm{rank}_{\rm symp}H_D=2.
```

Let $`m=4r=2^p`$ and consider

```math
E=L(D^{\otimes r}\otimes I_q)K,
```

where K and L are arbitrary fixed joint Cliffords, allowed to mix every
cell and helper. Its domain has dimension at least $`m+2q`$ and
symplectic rank at least $`2r+2q`$. These invariants survive both walls.
The factor $`I_q`$ describes the candidate E; it does not require the
source query to be identity on those physical coordinates.

For any degree-d Clifford source, the first theorem already excludes this
class below $`2^{-d-1}`$ whenever

```math
m\gt4M_d(p).
```

In particular, affine queries have an error floor of one-quarter for
$`p\ge5`$. The field structure strengthens the statement at smaller sizes.

Suppose $`p\ge3`$ and $`\epsilon\lt1/12`$. Projection onto bank X
masks and all helper labels has dimension at most $`p+1+2q`$, while
the domain has dimension at least $`m+2q`$. Its pure-Z kernel therefore
satisfies

```math
\dim K_Z\ge m-p-1,\qquad
[\mathbb F_2^m:K_Z]\le2^{p+1}=2m.
```

There are $`1+m+\binom m2`$ masks of weight at most two, more than
$`2m`$ because their difference is $`(m-1)(m-2)/2\gt0`$. Two distinct
such masks lie in the same coset of $`K_Z`$; their difference is a
nonzero kernel word of weight at most four. Since $`m\ge8`$, that word
is not all-ones. But $`\epsilon\lt1/12\le1/(3w)`$ contradicts the
low-weight lemma. Thus

```math
\epsilon\ge\frac1{12}\quad(p\ge3),\qquad
\epsilon\ge\frac14\quad(p\ge5).
```

These floors include the Pauli-helper sources of Section 5 and arbitrary
joint outer Clifford walls. They do not exclude coupled non-Clifford bands,
non-Clifford spectator factors, or a different query family.

At $`p=2`$ the actual four-address native decoder saturates the exact
field bound: its image domain has dimension four, its X-mask projection is
all of $`\mathrm{RM}(1,2)`$ with dimension three, and its pure-Z kernel
is exactly $`\{0,\boldsymbol1\}`$. Identity spectators add $`2q`$
to the domain dimension. The permitted weight-four kernel word is then
all-ones, so the uniform-error contradiction correctly starts at $`p=3`$.

## 8. Exact fixtures and scope

The [bounded exact fixture](../../tests/test_weak_probe_pauli_domain.py)
checks literal scalar cube phases, including a noncommuting negative cube;
Reed--Muller derivative and radical dimensions; joint Clifford invariance;
the actual native four-address domain; the weight-four code boundary; and
the signed helper-label distinction in Section 6. The analytic arguments
supply the all-size conclusions. The fixture is neither a uniform native
emitter nor a proof of a growing depth lower bound.

The complete-probe freedom remains useful, but it cannot rescue a tensor
product of the finite cells by adding only outer Clifford mixing. A new
construction must alter the non-Clifford coupling or its Pauli domain while
still paying the complete probe, actual inverse and helper-return contract.
