# Clifford-orbit rank of literal affine queries

[Interleaved tree probe](INTERLEAVED_TREE_PROBE_SOURCE.md#full-rank-probe-algebra) ·
[Fixed-readout syndrome bound](FIXED_READOUT_SYNDROME_BOUND.md) ·
[Clifford involution source](CLIFFORD_INVOLUTION_SOURCE.md)

A bounded word of literal affine XOR queries cannot generate the complete
interleaved tree-probe family from a fixed Hermitian Clifford, even with an
arbitrary fixed final decoder and arbitrarily many dirty helper wires.
The obstruction counts independent elements of a projective group. It does
not identify group rank with the linear rank of matrix entries.

The theorem first uses an independently certified degree bound on a faithful
binary lift. For the actual affine XOR-query grammar below, each query has
lift degree at most two. The resulting necessary condition is
$`q\ge\lceil p/8\rceil`$, at any complete operator error
$`\delta\lt\sqrt{2-\sqrt2}/4`$. This is a bound on that query word,
not an unrestricted native-depth lower bound.

<a id="faithful-phase-lift"></a>
## 1. A binary lift retaining normal-Pauli phases

On any finite number $`w`$ of work wires, write a normal Pauli as

```math
i^{t_0+2t_1}X^aZ^b,\qquad
a,b\in\mathbb F_2^w,\quad t_0,t_1\in\mathbb F_2.
```

Put $`z=(a,b,t_0)`$. The feature column consists of the constant one,
all entries of $`z`$, their square-free pair products, and $`t_1`$:

```math
\Phi(a,b,t_0,t_1)
 =\bigl(1,(z_i)_i,(z_i z_j)_{i\lt j},t_1\bigr)^{\mathsf T}.
```

Every Clifford conjugation acts linearly on these features. To see the
phase explicitly, order the normal generators as all X generators followed
by all Z generators. Write their images as
$`i^{k_j}X^{r_j}Z^{s_j}`$, and write $`v=(a,b)`$. The image label is
linear in $`v`$, while the normal phase increment modulo four is

```math
q(v)=\sum_j k_jv_j+
 2\sum_{j\lt k}(s_j\cdot r_k)v_jv_k\pmod4.
```

Its low bit $`q_0`$ is linear and its high bit $`q_1`$ is quadratic:
the carry from a sum of odd coefficients is the sum of their pair
products. Adding the incoming phase therefore gives

```math
t'_0=t_0+q_0(v),\qquad
t'_1=t_1+q_1(v)+t_0q_0(v).
```

Thus $`z'`$ is linear in $`z`$, and $`t'_1`$ has the required form.
The features $`t_0a_j,t_0b_j`$ are included; dropping them would lose
the carry. The feature functions are distinct Boolean monomials, so are
linearly independent as functions. Their evaluation columns span the
feature space. This gives a unique invertible binary matrix $`L(C)`$ and

```math
\Phi(CPC^\dagger)=L(C)\Phi(P),\qquad
L(CD)=L(C)L(D).
```

The kernel consists exactly of scalar Cliffords: an identity lift fixes
every normal Pauli, and an operator commuting with every X and Z is scalar.
Hence this representation is faithful on the projective Clifford group.
It retains Pauli signs while discarding only the global phase of the
conjugating Clifford.

<a id="literal-affine-query-degree"></a>
## 2. What the literal query grammar certifies

One query is an actual complete affine XOR shear between two disjoint
work banks, with fixed address $`x\in\mathbb F_2^p`$:

```math
Q_x:\ |u,v\rangle\longmapsto
 |u,v+M(x)u+h(x)\rangle.
```

Every entry of $`M`$ and $`h`$ is affine in the address. Both work banks
are arbitrary inputs, and the query is its own inverse. Its exact action
on normal-Pauli coordinates is

```math
\begin{aligned}
a'_u&=a_u,&a'_v&=a_v+M(x)a_u,\\
b'_u&=b_u+M(x)^{\mathsf T}b_v,&b'_v&=b_v,\\
t'_0&=t_0,&t'_1&=t_1+b_v\cdot h(x).
\end{aligned}
```

The linear feature coordinates have affine address coefficients; their
pair products have degree at most two. Therefore
$`\deg L(Q_x)\le2`$, including all normal-section phases. Either bank
orientation and fixed relabelings obey the same bound.

Let $`C_x`$ contain $`q`$ such queries, interleaved with arbitrary fixed
address-independent Clifford walls on all work and helper wires. Its
literal inverse contains the same number of queries. For a fixed
Hermitian Clifford $`R`$, put $`R_x=C_xRC_x^\dagger`$. Matrix
multiplication over Boolean functions gives

```math
\deg L(C_x),\ \deg L(C_x^{-1})\le2q,
\qquad \deg L(R_x)\le4q.                         \qquad(1)
```

This does **not** follow from affine symplectic entries alone. For example,
$`C_x=Z^{f(x)}`$ has identity symplectic matrix for any Boolean function
$`f`$, but its phase-lift row contains $`f(x)`$. Arbitrary address-dependent
sign choices are outside the certified grammar. A more general stage is
eligible if its full lift degree is proved independently; overlapping
query products are not silently charged as one affine query.

<a id="multiplicative-degree-rank"></a>
## 3. Low degree bounds multiplicative group rank

Let $`A(x)`$ be a family of commuting involutive binary matrices, with
entrywise Boolean degree at most $`d`$. Write its ANF as

```math
A(x)=I+\sum_{|S|\le d}x^S N_S,
\qquad x^S=\prod_{i\in S}x_i.
```

Each $`N_S`$ is a binary linear combination of the values $`A(x)-I`$.
Those values commute and square to zero. In characteristic two their
linear combinations consequently commute and square to zero too.

For any commuting square-zero elements $`n_j`$ and binary scalars
$`\alpha_j`$, the following identity is exact:

```math
1+\sum_j\alpha_jn_j
 =\prod_j(1+\alpha_jn_j)
  \prod_{j\lt k}(1+\alpha_j\alpha_kn_jn_k).       \qquad(2)
```

Prove it first in the free commutative square-zero algebra and then
substitute the matrices. A coefficient on a set of $`k`$ distinct
generators counts partitions into singletons and pairs, equivalently
involutions on that set. Their counts satisfy
$`I_0=I_1=1`$, $`I_2=2`$ and
$`I_k=I_{k-1}+(k-1)I_{k-2}`$. They are even for every $`k\ge2`$.
Only the constant and singleton terms survive.

Apply (2) with $`\alpha_S=x^S`$ and $`n_S=N_S`$. Pair factors have
exponent $`x^{S\cup T}`$. Collect all factors with the same address
monomial into fixed commuting involutions $`H_T`$:

```math
A(x)=\prod_{|T|\le2d}H_T^{x^T},\qquad
\mathrm{rank}\langle A(x):x\in\mathbb F_2^p\rangle
 \le M_{2d}(p):=\sum_{j=0}^{\min(2d,p)}\binom pj. \qquad(3)
```

The rank in (3) is the rank of an elementary abelian group. No linear
independence of matrix entries is asserted. In particular, multiplicatively
independent involutions can be linearly dependent as matrices; the pair
factors in (2) are essential.

<a id="robust-projective-independence"></a>
## 4. A fixed decoder and a nondegenerate target algebra

Let $`P_x`$ be Hermitian Paulis, including identity on every arbitrary
helper, indexed by all $`m=2^p`$ addresses. Assume their commutator matrix
$`K`$ is nonsingular over $`\mathbb F_2`$:

```math
P_xP_yP_xP_y=(-1)^{K(x,y)}I.
```

Let $`R_x`$ be Hermitian Cliffords. Suppose a fixed unitary $`U`$ on the
complete work space, independent of the address and otherwise arbitrary,
satisfies

```math
\max_x\|UR_xU^\dagger-P_x\|\le\delta,
\qquad \delta\lt\frac{\gamma}{4},
\qquad \gamma=\sqrt{2-\sqrt2}.                    \qquad(4)
```

A nonscalar Clifford has normalized trace of magnitude at most
$`1/\sqrt2`$. This follows from
[flat Pauli support](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension):
a nonscalar singleton support is a nonidentity Pauli and has zero trace;
otherwise every nonzero coefficient has magnitude at most $`1/\sqrt2`$.
The normalized Hilbert--Schmidt norm therefore gives, for either sign,

```math
\|C-sI\|^2\ge
 2-2\mathrm{Re}\bigl(s\,\mathrm{tr}_{\mathrm{norm}}C\bigr)
 \ge2-\sqrt2,\qquad s\in\{1,-1\}.              \qquad(5)
```

Four-factor telescoping in (4) puts the Clifford commutator
$`R_xR_yR_xR_y`$ within $`4\delta\lt\gamma`$ of the required scalar
sign. Equation (5) forces it to be scalar. A scalar commutator of two
involutions is necessarily $`\pm I`$, since conjugating twice squares
that scalar. The wrong sign has distance two, so

```math
R_xR_yR_xR_y=(-1)^{K(x,y)}I.                       \qquad(6)
```

Consequently the faithful lifts $`L(R_x)`$ commute and square to the
identity, although the actual operators need not commute. A projective
relation $`\prod_xR_x^{c_x}`$ being scalar would imply
$`K^{\mathsf T}c=0`$ by commuting it past each generator. Nonsingularity
forces $`c=0`$. There are therefore exactly $`m`$ independent lifted
generators. If their independently certified entry degree is at most
$`d`$, (3) gives the necessary condition

```math
M_{2d}(p)\ge m=2^p,\qquad 2d\ge p.               \qquad(7)
```

<a id="tree-orbit-query-bound"></a>
## 5. Consequences for the complete tree probe

The [pulled tree-probe algebra](INTERLEAVED_TREE_PROBE_SOURCE.md#full-rank-probe-algebra)
has the commutator matrix

```math
K_p(x,y)=\sum_{i=1}^p(x_i+y_i)
                    \prod_{j>i}(1+x_jy_j).
```

It is nonsingular for every $`p\ge1`$. Its last-bit block recursion and
base case are

```math
K_p=\begin{pmatrix}K&K+J\\K+J&0\end{pmatrix},
\qquad K_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},
```

where $`K=K_{p-1}`$ and $`J`$ is the all-ones matrix.
For nonsingular alternating $`K`$,
$`\mathbf1^{\mathsf T}K^{-1}\mathbf1=0`$, making $`K+J`$ and then
the displayed block matrix invertible.

Apply (7) to the orbit in (1). Any fixed decoder achieving (4) must obey

```math
M_{8q}(p)\ge2^p,
\qquad q\ge\lceil p/8\rceil.                    \qquad(8)
```

The source Clifford $`R`$ may be non-Pauli. The width, fixed Clifford
walls, and fixed final decoder have no restriction beyond the stated
complete-space action. This permits the specified later address-dependent
tree processing only because its pulled target algebra is explicitly
known. Arbitrary later address processing, unrestricted non-Clifford
first images, and a general native synthesis remain outside the theorem.

For the original commuting complete lookup target, exact fixed decoding
also requires (8), by independence of its projective Pauli labels. Under
(4), the same commutator argument instead forces the actual source family
to commute. A commuting Clifford group of projective rank $`r`$ realizes
at most $`2^r`$ joint sign patterns: if $`-I`$ is present, its character
is fixed and supplies no extra pattern. The existing
[syndrome-covering proof](FIXED_READOUT_SYNDROME_BOUND.md) then gives the
compact approximate corollary

```math
M_{8q}(p)\ge2^p\bigl[1-H_2(\delta^2/4)\bigr].
```

The exact fixtures in
[`test_clifford_orbit_query_rank.py`](../../tests/test_clifford_orbit_query_rank.py)
check the phase lift, literal query degree, square-zero factorization,
Clifford trace gap and projective independence. The proofs above supply
the all-width and all-address statements.
