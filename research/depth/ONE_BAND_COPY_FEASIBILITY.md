# Exact one-band copy feasibility from signed Walsh blocks

[Complete copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Native copy components](NATIVE_COPY_INTERACTION.md) ·
[Copy-axis span](COPY_AXIS_SPAN.md)

A literal one-band copy admits an exact coefficient-space reduction.
Its potential input and output supports are affine binary slices. On each
physical T block, the transfer between those slices is a signed Walsh
matrix with an explicitly known spectrum. This identifies the full
coefficient space that can avoid output leakage, with a uniform gap for the
remaining coefficients. A separate carrier fixed-point condition is
essential: even all physical output Pauli tests do not enforce it.

This is a complete linear feasibility criterion and a tool for
constructing or rejecting a specified copy. It supplies neither a successful all-size
field seed nor a new frame-depth theorem. The native seed and its inverse
remain separate synthesis obligations.

Two further necessary conditions guide the remaining search: the
[Clifford-seed one-band floor](CLIFFORD_SEED_COPY_STRUCTURE.md#one-band-clifford-seed)
forces a genuinely non-Clifford seed at growing address size, while the
[dual-group grading](CLIFFORD_DUAL_COPY_GRADING.md) applies to arbitrary
unitary seeds. Neither replaces the complete fixed-point equations below.

An explicit triple-axis joint-wall family below pays 3m T gates in one
band and passes the binary tensor check, but has zero complete seed space
for $`p\ge2`$. Its full copy error is at least $`1/\sqrt6`$ for any
unitary seed, including arbitrary dirty helpers. This family-specific
result illustrates why the carrier condition cannot be omitted.

## 1. Physical one-band coordinates

Let A denote the complete arbitrary carrier, including every private
helper, and let Y be the physical m-qubit output. Write a proposed copy
as

```math
C=F D G,\qquad D=\prod_{i\in B}T_i^{\epsilon_i},
\qquad \epsilon_i\in\{+1,-1\},
```

where F and G are arbitrary fixed Cliffords on A and Y. The gates in D
act on distinct physical wires. Scalar phases remain in the native word;
they cancel only in conjugation. Each actual inverse is obtained by
reversing the emitted word and replacing each gate by its adjoint.

Use Hermitian Paulis
$`\tau(a,z)=i^{\mathrm{wt}(a\mathbin{\&}z)}X^aZ^z`$, orthonormal
in the normalized Hilbert--Schmidt inner product. A Clifford permutes
their binary labels symplectically and supplies a literal sign. For D,
fix the X label a and the Z coordinates outside

```math
J(a)=\{i\in B:a_i=1\},\qquad r=|J(a)|.
```

The free Z coordinates in J form one block. Conjugation on a single
active wire, in the ordered X,Y basis, is HZ for T and ZH for its
adjoint. Therefore the exact block is

```math
K=Z^v H^{\otimes r} Z^u,\qquad
K_{y,z}=2^{-r/2}(-1)^{y\cdot z+u\cdot z+v\cdot y},
                                                               \qquad(1)
```

where u marks positive T orientations and v marks negative ones.
They are disjoint and sum to the all-one vector. Inactive I/Z factors
are unchanged. This also covers r=0: the block is the scalar identity.

The binary carrier-only Pauli space, transformed by G, intersects this
block in either the empty set or an affine slice
$`S=s_0+L\subseteq\mathbb F_2^r`$. For a chosen physical output
Pauli $`X_{Y_j}`$, the labels of all carrier Paulis times that output,
pulled back through F, intersect it in an empty set or an affine slice
$`T=t_0+M`$. These slices follow by binary linear elimination; signs
of F and G are retained when transforming coefficients. No matrix on
the exponentially large physical Hilbert space is required to find them.

## 2. Spectrum and complete no-leakage space of one block

Assume S and T are nonempty, and put

```math
l=\dim L,\quad t=\dim M,\quad
h=\dim(L\cap M^\perp).
```

The restriction $`P_T K P_S`$ has rank and nonzero squared singular
values

```math
\mathrm{rank}(P_T K P_S)=2^{l-h},\qquad
\lambda=2^{t+h-r}.                                  \qquad(2)
```

In particular there exists a nonzero vector supported in S whose full
image is supported in T if and only if

```math
M^\perp\subseteq L.                                \qquad(3)
```

When (3) holds that space has dimension $`2^{l+t-r}`$. When it fails,
the restriction has norm at most $`1/\sqrt2`$. If either slice is
empty the restriction is zero; there is no supported nonzero vector.
All assertions include r=0 with the usual zero-dimensional conventions.

### Proof of the spectrum

The row and column signs in (1) do not change singular values. For the
unsigned Walsh matrix the restricted Gram kernel, for z,z' in S, is

```math
(P_S H P_T H P_S)_{z,z'}
 =\begin{cases}
 2^{t-r}(-1)^{t_0\cdot(z+z')},&z+z'\in M^\perp,\\
 0,&z+z'\notin M^\perp.
 \end{cases}
```

Partition S into cosets of $`L\cap M^\perp`$. Each coset gives one
rank-one signed matrix with its sole nonzero eigenvalue
$`2^{t-r}2^h`$. There are $`2^{l-h}`$ such cosets, proving (2).
Since $`h\le r-t`$, a norm-one singular direction occurs exactly
when equality holds, which is (3). Otherwise the integer exponent in
(2) is at most minus one. Unitarity of the full Walsh transform then
identifies norm-one directions precisely with absence of leakage to
the complement of T.

### Constructive real basis

When (3) holds, set $`H_0=M^\perp`$ and choose one representative s
from each coset of H_0 in S. The following real unit vectors give an
orthonormal basis of the no-leakage input space:

```math
f_s(z)=2^{-\dim H_0/2}
 (-1)^{u\cdot z+t_0\cdot(z+s)}[\,z\in s+H_0\,].
                                                               \qquad(4)
```

Their exact images are

```math
(Kf_s)(y)=2^{-t/2}(-1)^{(s+v)\cdot y}[\,y\in T\,]. \qquad(5)
```

Indeed, write $`g(z)=(-1)^{u\cdot z}f(z)`$ to remove the input
sign. The Fourier-support condition is then
$`g(z+h')=(-1)^{t_0\cdot h'}g(z)`$ for every $`h'\in H_0`$.
If a nonzero coefficient lies in S, its entire H_0 coset must lie
there, forcing (3). Conversely (4) has precisely that character and
Fourier transforms to (5). Clifford-wall signs multiply these explicit
coefficients when returning to physical carrier labels; they are not
discarded as a phase convention.

## 3. Full copy requires the carrier fixed point

Write $`\iota(A)=A\otimes I_Y`$ and use normalized Hilbert--Schmidt
norms on both spaces, so iota is isometric. Define the physical
carrier coefficient map

```math
\Theta_j(A)=2^{-m}\mathrm{Tr}_Y
 [(I_A\otimes X_{Y_j})\,C(A\otimes I_Y)C^\dagger].  \qquad(6)
```

It is a Hilbert--Schmidt contraction. Orthogonal projection of the
conjugated operator onto the output-X_j sector is exactly
$`\Theta_j(A)\otimes X_{Y_j}`$. Consequently

```math
\Theta_j(A)=A
\quad\Longleftrightarrow\quad
C(A\otimes I_Y)C^\dagger=A\otimes X_{Y_j}.          \qquad(7)
```

For the forward implication, the projected norm already equals the
norm of the full operator. Every other output-Pauli component must
therefore vanish. The converse is immediate. Equation (7) holds for
arbitrary operators, not only unitaries.

The affine-block calculation constructs the full no-leakage subspace
$`\mathcal E_j`$ and the action of Theta on it: pull (4) back through
G and push (5) through F, retaining their signs. A compatible operator
must in addition satisfy $`\Theta_j(A)=A`$. A nonzero no-leakage
space by itself is insufficient.

For a fixed Clifford query family define $`\mathcal Q_x(A)=Q_xAQ_x^\dagger`$.
The complete compatible seed space is exactly

```math
\mathcal K=\bigcap_x\ker
 (I-\mathcal Q_x^{-1}\Theta_x\mathcal Q_x).
                                                               \qquad(8)
```

The block result reduces (8) to explicit real coefficient equations
on the no-leakage bases, rather than unconstrained dense operator
matrices. In general those equations may still be exponentially large.
An efficient symbolic description or an all-size solution must be
proved for the chosen walls; the reduction does not assert one.

Equivalently, let $`T_x=\mathcal Q_x^{-1}\Theta_x\mathcal Q_x`$.
Every $`I-(T_x+T_x^*)/2`$ is positive semidefinite, and

```math
\mathcal K=\ker\sum_x[I-(T_x+T_x^*)/2].              \qquad(9)
```

Equality in the contraction and Cauchy--Schwarz bounds proves that a
zero quadratic form for a summand forces $`T_x(A)=A`$. Positivity
then identifies the kernel of the sum with the common fixed space.

The same matrix gives a quantitative certificate for the complete
error, including any residual carrier action. There are m addresses;
define

```math
\mathcal L=\frac1m\sum_x(2I-T_x-T_x^*),\qquad
\Delta_x(A)=C\iota(\mathcal Q_x(A))C^\dagger
                 -\mathcal Q_x(A)\otimes X_{Y_x}.
```

Unitary invariance and orthogonal projection give the exact identity

```math
\frac1m\sum_x\|\Delta_x(A)\|_{\rm HS,norm}^2
       =\langle A,\mathcal L A\rangle_{\rm HS,norm}. \qquad(9a)
```

Thus $`\ker\mathcal L=\mathcal K`$. If the smallest eigenvalue
of this finite matrix is positive, every unitary seed has worst-case
complete operator error at least $`\sqrt{\lambda_{\min}(\mathcal L)}`$.
Indeed its normalized Hilbert--Schmidt norm is one, and the operator
norm bounds the normalized Hilbert--Schmidt norm at each address.
This is a computable certificate for specified walls; no uniform
positive eigenvalue is asserted for a family of increasing sizes.

### Full output duals are a prefilter, not a substitute

For one-band C, every $`J_v=C^\dagger P_Y(v)C`$ is a Hermitian
Clifford. Requiring

```math
J_v(R_x\otimes I_Y)J_v^\dagger
 =(-1)^{\langle X_{Y_x},v\rangle}(R_x\otimes I_Y)
```

for all output Pauli generators gives precisely the no-leakage
condition. To work in seed coordinates, conjugate each J_v by
$`Q_x^\dagger\otimes I_Y`$. Its adjoint action is a signed
permutation of Hermitian Pauli coefficients. For each carrier-only
label u, impose the resulting relation between its coefficient and
that of the image label, including the prescribed character sign.
An image outside the carrier-only space forces that coefficient to
zero. A signed component with such a zero constraint or a contradictory
sign cycle vanishes; every other component supplies one real basis
vector with its propagated signs. Signed union-find therefore computes
this full output-dual prefilter without dense physical matrices.
It can still have exponentially many coefficient labels.

Even all output X and Z generators do not enforce the carrier fixed
point. If $`C'=Z_A C`$ for a carrier Pauli, every output dual is
unchanged, whereas $`\Theta'_j(A)=Z_A\Theta_j(A)Z_A`$. A copied
operator not fixed by that carrier conjugation no longer satisfies
(7). The retained native one-bit-address copy for the seed H gives such
an example with $`ZHZ\ne H`$. Testing only output signatures would
accept it incorrectly.

## 4. Uniform gaps and balanced native seeds

If every nonempty input slice for a fixed j has an empty target slice
or fails (3), the norm of Theta is at most $`1/\sqrt2`$.
For a unitary seed element A, both complete
operators compared in (7) have normalized Hilbert--Schmidt norm one,
and their real overlap is at most $`1/\sqrt2`$. Hence the complete
operator error is at least

```math
\gamma=\sqrt{2-\sqrt2}.                            \qquad(10)
```

There is also a partial-leakage statement. For any fixed j, let E_j be
the orthogonal projector onto the norm-one subspace of the full
output-sector restriction. Every other singular value is at most
$`1/\sqrt2`$. For arbitrary A,

```math
\|C\iota(A)C^\dagger
 -\Theta_j(A)\otimes X_{Y_j}\|_{\rm HS,norm}^2
 \ge\tfrac12\|(I-E_j)A\|_{\rm HS,norm}^2.           \qquad(11)
```

In particular a complete copy defect at most eta for a unitary A
implies distance at most $`\sqrt2\eta`$ from this allowed coefficient
space. This is a necessary full-operator test; it does not make an
approximately allowed operator into a nearby native seed.

For a Hermitian-reflection construction, (8) must contain R with
$`R^2=I`$. It is automatically traceless: every nonzero output-X
sector has trace zero, while unitary conjugation preserves the trace.
In a real Hermitian basis $`E_\alpha`$ of K, the equation
$`(\sum_\alpha r_\alpha E_\alpha)^2=I`$ is an explicit quadratic
coefficient test. A nonzero vector in K does not suffice. Nor does an
abstract balanced involution pay a native factorization.

The bounded construction target chooses that factorization from the
start, $`R=V P V^\dagger`$, with P a nonidentity Pauli and V two
literal T bands with explicit Clifford interlayers. Together with a
one-band copy and the paid eight-band query, a completed identity
would give the conservative indicator depth bound 42. Its width,
T count and Clifford count still require one common ledger. This
Hermitian grammar is an intentional test, not an exhaustive class;
higher-order unitary seeds remain allowed by the complete interface.

## 5. Exact finite checks and constructive use

[Signed Walsh fixtures](../../tests/test_one_band_walsh_blocks.py)
exhaust affine slice pairs through dimension three, include selected
four-dimensional pairs, and check (2)--(5) with exact arithmetic.
They also recover the block signs from actual T/T-dagger conjugations,
including inactive coordinates. These are checks of the all-size
linear-algebra proof, not finite evidence for a successful field copy.

The [complete-kernel fixtures](../../tests/test_one_band_copy_feasibility.py)
use exact coefficients in $`\mathbb Q(\sqrt2)`$ and audit literal
native conjugation on all 84 Paulis at widths one, two and three.
For the retained one-bit-address copy and Pauli query $`Q=(I,Z)`$, the
complete compatible space is the line spanned by X+Z. Appending
carrier Z leaves every physical output X/Z dual unchanged but makes
the complete kernel zero; the normalized H seed has squared defect
matrix $`\Delta^\dagger\Delta=2I`$. An arbitrary carrier helper
enlarges the original kernel to the expected four directions
$`H\otimes\{I,X,Y,Z\}`$. The p=2 tensor-product copy has eight
T gates scheduled in one physical band and a zero complete field-query
kernel. That last rejection concerns this specified finite copy.

For a chosen uniform wall family, first compute its binary slices and
remove blocks failing (3). Next impose the exact reduced carrier
fixed-point equations (8), and finally enforce the chosen native seed
factorization. A vanishing complete kernel rejects that exact instance;
a proved uniform positive error floor rejects the family below that
floor. Surviving support, a balanced dense matrix,
or an isolated small-size solution does not establish the desired
uniform native construction.

## 6. An overlapping joint-wall family passes the tensor check but fails full copy

The following uniform family is a concrete application of the complete
criterion. Its walls jointly mix carrier and output, and every core
Pauli survives a necessary binary tensor check. Nevertheless its complete
kernel is zero for $`p\ge2`$, even with arbitrary operators on its helpers.
The same family has complete copy error at least $`1/\sqrt6`$ for every
unitary seed. These conclusions concern the specified walls, not all
one-band copies.

### A necessary tensor check

For each active native T coordinate write the binary axis labels
$`a_i=FZ_iF^\dagger`$ and $`b_i=G^\dagger Z_iG`$; superscript Y
denotes output projection. For an output Pauli label v, the Clifford
$`J_v=C^\dagger P_Y(v)C`$ has symplectic action

```math
u\longmapsto u+
 \sum_i\langle a_i^Y,v\rangle\langle b_i,u\rangle b_i.
```

This follows from the single-coordinate T conjugation: each activated
coordinate contributes the transvection along b_i, and those axes are
independent and commute. If a carrier-only coefficient u acquires a
nonzero output label, its coefficient must vanish in the signed-dual
constraint. Taking every output v therefore gives the necessary test

```math
\sum_i\langle b_i,u\rangle\,b_i^Y\otimes a_i^Y=0. \qquad(12)
```

For a seed coefficient this must hold with u replaced by Q_x u for
every address. Equation (12) tests output labels, not their signs or
the carrier fixed point.

### Explicit native walls

Let $`m=2^p`$ and use the standard field query on its m core wires,
with u and v halves of size $`m/2`$. Index core wires j modulo m.
Let g_j be X on the j-th u wire or Z on the corresponding v wire.
These axes are independent and commute. Add two arbitrary helpers
$`h_{j1},h_{j2}`$ per core wire and define

```math
\ell_{j1}=g_jZ_{h_{j1}},\qquad
\ell_{j2}=g_jZ_{h_{j2}},\qquad
\ell_{j3}=g_jZ_{h_{j1}}Z_{h_{j2}}.
```

A literal Clifford K sends the three native Z axes of each triple
$`(j,h_{j1},h_{j2})`$ to these axes. Number its positions 0,1,2 and
apply, chronologically,

```math
\mathrm{CX}(2;1),\quad\mathrm{CX}(1;2),\quad
\mathrm{CX}(1;0),\quad\mathrm{CX}(0;2),
```

then H on every u core wire. Here CX lists control before target.
This fixes all Pauli signs; before those H gates its binary basis map is

```math
L=\begin{pmatrix}1&1&1\\0&1&1\\1&0&1\end{pmatrix},\qquad
L^{-T}=\begin{pmatrix}1&1&1\\1&0&1\\0&1&1\end{pmatrix}.
```

Define binary matrices alpha and beta with three rows per j:

```math
\alpha_j=(e_j,e_j,e_j),\qquad
\beta_j=(e_j,e_{j+1},e_j+e_{j+1}).
```

Let E_M be the commuting CNOT word from output Y_t to carrier
coordinate i whenever M_it=1, and let H_Y be H on every output. Set

```math
F=(K\otimes H_Y)E_\alpha,\qquad
G=E_\beta(K^\dagger\otimes H_Y),\qquad
D=T^{\otimes3m},\qquad C=FDG.                    \qquad(13)
```

The word has exactly 3m T gates in one physical band, 18m elementary
Clifford gates and 4m total wires, of which 2m are arbitrary helpers.
It uses zero clean inputs. The actual inverse reverses the word and
adjoints each gate. The count includes both K walls, both H_Y walls,
3m CNOTs in E_alpha and 4m in E_beta. Cyclic overlap of the output
controls makes the displayed walls jointly mix carrier and output.
No claim is needed about alternative representations of the same unitary.

Its axis projections are $`a_i=(\ell_i,X_Y^{\alpha_i})`$ and
$`b_i=(\ell_i,X_Y^{\beta_i})`$. Every tensor in (12) is nonzero.
Within each triple the three tensors have their single relation;
different triples have independent second factors e_j. Thus (12) says
the three activations $`\langle\ell_{ji},u\rangle`$ agree.
At x=0 this forces all helper X labels to vanish. Conversely every
core Pauli times every helper-Z word passes for every x, since Q_x
fixes the helpers. The test removes no core Pauli label.

### Complete coefficient action and its contradiction

In the output X-eigenbasis with eigenvalue bits y, C has carrier blocks

```math
C_y=KX^{\alpha y}T^{\otimes3m}X^{\beta y}K^\dagger.
```

Full copy requires $`C_yR_xC_y^\dagger=(-1)^{y_x}R_x`$.
In particular $`[C_0,R_x]=0`$, so the same equation holds with C_y
replaced by the Clifford $`J_y=C_yC_0^\dagger`$. This use of the
y=0 branch enforces the carrier return missing from output signatures.

Express the core in the g_j basis, by conjugating only with H on its
u half. A surviving Hermitian basis element has core labels (w,z) and
helper labels (r_1,r_2). In one triple its native labels after K's
inverse are

```math
s=(w,w,w),\qquad q=(z+r_2,z+r_1,z+r_1+r_2),
```

with Hermitian-Pauli gauge $`(-1)^{w(z+r_1+r_2+r_1r_2)}`$.
For one native coordinate, conjugation by
$`X^A(TX^BT^\dagger)`$ sends q to q+Bs with sign
$`(-1)^{Aq+ABs+Bq(1-s)}`$. Combining the three coordinates and
their gauges shows that J_y preserves the physical core label, sends

```math
(r_1,r_2)\longmapsto(r_1+w y_{j+1},r_2+w y_j),
```

and contributes sign

```math
(-1)^{y_jz+y_jr_1+y_{j+1}r_2
             +w(y_j+y_{j+1}+y_jy_{j+1})}.          \qquad(14)
```

All exponents here are binary. In particular J_{e_t} acts on the
vector of coefficients indexed by helper-Z masks as the Pauli

```math
D_t(w,z)=(-1)^{z_t+w_t+w_{t-1}}
 Z_{h_{t1}}X_{h_{t2}}^{w_t}
 Z_{h_{(t-1)2}}X_{h_{(t-1)1}}^{w_{t-1}}.          \qquad(15)
```

The X and Z symbols in (15) act on the coefficient-index space, not
on physical helper states. All displayed factors act on different
index wires. The gauge in (14) is needed for the literal sign in (15).

Take a nonzero core-Pauli coefficient vector r_u of a putative seed,
including all its helper-Z coefficients. Q_x sends its core label to
u_x and supplies one overall Pauli sign to that vector. This sign
cancels in the copy equation. Because the field symplectic map is
affine in x, the coordinates w(x),z(x) of u_x are affine. Complete
copy therefore forces

```math
D_t(u_x)r_u=(-1)^{[x=t]}r_u\qquad\hbox{for every }x,t. \qquad(16)
```

For fixed t, (15) uses fixed commuting coefficient Paulis with affine
exponents and an affine scalar sign. Its product over every affine
address two-plane is consequently literal +I. When $`p\ge2`$, choose
such a plane containing t. It contains t exactly once, so multiplying
the four equations (16) gives $`r_u=-r_u`$. Every coefficient vanishes.
Thus $`\mathcal K=\{0\}`$ for (13), for arbitrary seed depth and
arbitrary operators on the full carrier and its helpers. Extra identity
spectators only tensor the same argument with their coefficient spaces.

The [native family fixtures](../../tests/test_joint_wall_copy_family.py)
check K, the tensor survivor space and the literal branch action (14)
directly from emitted native gates. The obstruction follows from the
all-size affine-plane identity above.

### A uniform error floor for the same walls

Let R be any unitary on the complete carrier, with arbitrary helpers,
and let the worst complete operator error be

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                  -R_x\otimes X_{Y_x}\|,\qquad R_x=Q_xRQ_x^\dagger.
```

Let P_good project in coefficient space onto operators having no helper
X labels, and write $`R_g=P_{\rm good}R`$, $`R_b=(I-P_{\rm good})R`$.
All Hilbert--Schmidt norms below are normalized; unitarity gives
$`\|R_g\|^2+\|R_b\|^2=1`$.

In a native D block with active coordinates J, the input free labels
obey $`\beta_J^Tq=\kappa_{\rm in}`$ and the target free labels obey
$`\alpha_J^Tq'=e_x+\kappa_{\rm out}`$, with constants fixed by the
block. Their direction spaces are therefore
$`\ker\beta_J^T`$ and $`\ker\alpha_J^T`$. By (3), a norm-one
singular direction requires

```math
\beta_J^T\alpha_J=0.
```

For these triple walls that condition says exactly that the three
activations agree, or equivalently that the physical helper X labels
vanish. Thus the full no-leakage subspace is contained in the good
space; equality is not required. Every bad block has squared singular
norm at most one half. The leakage estimate (11), and the fact that
Q_x preserves helper labels, give

```math
\|R_b\|_{\rm HS,norm}^2\le2\eta^2.                \qquad(17)
```

For every x the branch errors at y=0 and y=e_t are at most eta in
operator norm, since C is block diagonal in the output X-eigenbasis.
Their triangle inequality gives

```math
\|J_{e_t}R_xJ_{e_t}^\dagger-(-1)^{[x=t]}R_x\|\le2\eta.
```

The good space is preserved by Q_x and by conjugation with J_y, as
(14) shows, so its orthogonal projector commutes with both actions.
Projecting and pulling back by Q_x gives residual at most $`2\eta`$
in normalized Hilbert--Schmidt norm for R_g. Write these pulled-back
maps as $`\mathcal B_x=\mathcal Q_x^{-1}\mathcal J_t\mathcal Q_x`$,
where $`\mathcal J_t(A)=J_{e_t}AJ_{e_t}^\dagger`$. For fixed t,
the four maps on the good space over an address two-plane containing t
are commuting Hermitian involutions, by (15), and their product is +I.
The four desired signs have product minus one. In every simultaneous
eigenspace at least one sign disagrees, contributing squared residual
four. Hence

```math
4\|R_g\|_{\rm HS,norm}^2
 \le\sum_{x\text{ in the plane}}
 \|\mathcal B_x(R_g)-(-1)^{[x=t]}R_g\|_{\rm HS,norm}^2
 \le16\eta^2.                                    \qquad(18)
```

Combining (17) and (18) yields the uniform full-operator floor

```math
1\le6\eta^2,\qquad \eta\ge\frac1{\sqrt6},\qquad p\ge2. \qquad(19)
```

No Hermiticity or seed-depth assumption was used. Additional identity
spectators are allowed. The native fixtures also check a local exact
half-gap block and the four-sign spectral certificate. This retires
the displayed wall family below (19); other one-band joint walls remain
open.
