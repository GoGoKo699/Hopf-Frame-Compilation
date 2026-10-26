# The Pauli algebra after an interleaved tree query

[Minimal-tree source and paid XOR query](ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source) ·
[Fixed endpoint syndrome bound](FIXED_READOUT_SYNDROME_BOUND.md) ·
[Individual Clifford involution sources](CLIFFORD_INVOLUTION_SOURCE.md)

A fixed Fourier wall after the complete tree XOR query changes the original
commuting lookup family into a full-rank noncommuting Pauli family when the
observable is pulled back through that wall. This gives an exact target
algebra for an interleaved source search. The wall and the query are both
paid on the complete input space. Normalizing this target algebra by a
fixed Clifford does not supply a cheap simultaneous decoder for an
arbitrary family of source Clifford involutions.

## 1. The paid wall and its literal pulled-back probe

Let $`p\ge1`$, $`m=2^p`$, and index each of two arbitrary m-qubit banks c,t
by subsets $`S\subseteq[p]`$. For nonempty S put
$`i(S)=\min S`$ and $`\pi(S)=S\setminus\{i(S)\}`$. Use column vectors
throughout this note:

```math
T_xe_S=x_{i(S)}e_{\pi(S)},\qquad T_xe_\varnothing=0,
\qquad A_x=I+T_x,\qquad
(v_x)_S=\prod_{i\in S}x_i.                              \qquad(1)
```

The empty product is one. Write $`X_c(a),Z_c(b)`$ for Pauli words on the
indicated bank, and similarly for t. The full query from the minimal-tree
source note is

```math
Q_{{\rm tree},x}|c,t,h\rangle
   =|c,t+A_xc,h\rangle.                                \qquad(2)
```

Here h consists of $`m-1`$ arbitrary dirty helpers returned exactly. The
physical query also preserves the coherent address x. Its actual reversed
native word implements its adjoint; the fact that (2) is an involution
does not replace that word by an unpaid inverse.

Define a fixed Clifford, in rightmost-first operator order, by

```math
C=\left(\prod_S\mathrm{CX}(c_S;t_S)\right)
  \left(\prod_S H(c_S)\right),\qquad
B_x=CQ_{{\rm tree},x}.                                 \qquad(3)
```

C uses exactly m Hadamards and m CNOTs. In particular
$`C^\dagger Z_t(v)C=X_c(v)Z_t(v)`$ for every v. The complete native B word
is the paid query followed by these two Clifford walls. For B and its
actual adjoint, common bounds are

```math
\begin{aligned}
T_B&\le8(m-1),& D_T(B)&\le8,\\
G_{{\rm Cliff},B}&\le12m-8,&G_B&\le20m-16,\\
w_{\rm private}&=m-1,&a_{\rm clean}&=0,\qquad\epsilon=0.
\end{aligned}                                                   \qquad(4)
```

$`G_B`$ counts every elementary gate, including T gates. The two banks
are allocated separately: the total width including the address is
$`p+3m-1`$. Both orientations have literal full-space action and return
all private helpers, including after tensoring an arbitrary reference.

The addressed observable to be read is $`Q_x=Z_t(v_x)`$. Its pullback is

```math
\begin{aligned}
P_x&=B_x^\dagger Q_xB_x\\
   &=X_c(v_x)X_t(A_xv_x)Z_c(A_x^{\mathsf T}v_x)Z_t(v_x).
\end{aligned}                                                   \qquad(5)
```

Equation (5) is in the literal normal order X then Z, with no suppressed
scalar. It follows directly from the full XOR conjugation identities
$`X_c(a)\mapsto X_c(a)X_t(A_xa)`$ and
$`Z_t(b)\mapsto Z_c(A_x^{\mathsf T}b)Z_t(b)`$. Thus the binary labels are

```math
a_x=(v_x,A_xv_x),\qquad b_x=(A_x^{\mathsf T}v_x,v_x).
                                                               \qquad(6)
```

In fact $`A_x^{\mathsf T}v_x=e_\varnothing`$: at each nonempty S,
the two entries $`(v_x)_S`$ and $`x_{i(S)}(v_x)_{\pi(S)}`$ cancel.
The root entry remains one. Thus the c-bank Z factor in (5) is always
the single root Z, although its X factor contains the full lookup vector.

Their binary overlap is zero because the two terms in
$`a_x\cdot b_x=v_x\cdot A_x^{\mathsf T}v_x+(A_xv_x)\cdot v_x`$
are equal. The normal word is therefore a Hermitian involution. If the
Hermitian Pauli section is instead
$`\tau(a,b)=i^{|a\wedge b|}X^aZ^b`$, (5) has the explicit sign
$`(-1)^{|a_x\wedge b_x|/2}`$ relative to that section. For example, at
$`x=0`$ it is $`-Y_{c_\varnothing}Y_{t_\varnothing}`$ on the root pair.
Using the unsigned Hermitian section would lose this literal sign.

<a id="full-rank-probe-algebra"></a>

## 2. Full commutation rank and canonical target coordinates

The m vectors $`v_x`$ form the subset-zeta matrix: its entries are one
exactly when $`S\subseteq\mathrm{supp}(x)`$. Ordering subsets by
inclusion makes it triangular with diagonal one. Consequently the m
labels in (6) are independent, already in their first X block.

Let the binary commutation matrix be

```math
K_p(x,y)=a_x\cdot b_y+b_x\cdot a_y,
\qquad P_xP_y=(-1)^{K_p(x,y)}P_yP_x.
```

The exact all-size formula is

```math
K_p(x,y)=\sum_{i=1}^p(x_i+y_i)
                    \prod_{j>i}(1+x_jy_j).             \qquad(7)
```

Indeed, for any address vector z,
$`v_y^{\mathsf T}T_zv_x=\sum_i z_ix_i\prod_{j>i}(1+x_jy_j)`$:
the children with minimum i are $`\{i\}\cup U`$ with
$`U\subseteq\{i+1,\ldots,p\}`$. Expanding (6)'s pairing cancels the
identity parts of $`A_x+A_y`$ and leaves
$`v_x^{\mathsf T}(T_{x+y}+T_{x+y}^{\mathsf T})v_y`$, giving (7).

Split rows and columns by their final address bit. With J the all-ones
$`2^{p-1}`$-square matrix, (7) gives

```math
K_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
K_p=\begin{pmatrix}
K_{p-1}&K_{p-1}+J\\
K_{p-1}+J&0
\end{pmatrix}.                                        \qquad(8)
```

Suppose K is invertible alternating, and write $`J=uu^{\mathsf T}`$.
Its inverse is alternating as well: for $`z=Kw`$,
$`z^{\mathsf T}K^{-1}z=w^{\mathsf T}Kw=0`$. Hence
$`u^{\mathsf T}K^{-1}u=0`$, and the rank-one determinant identity gives
$`\det(K+J)=\det K`$. The off-diagonal blocks in (8) are therefore
invertible. A vector in (8)'s kernel has zero first component by the
second block row and zero second component by the first. Induction proves

```math
\mathrm{rank}_{\mathbb F_2}K_p=m.                      \qquad(9)
```

In particular the addresses with $`x_p=1`$ give a commuting independent
half of the family. They span a Lagrangian in the m-dimensional target
label space; the other half pairs with them through the invertible
$`K_{p-1}+J`$. Binary symplectic elimination first inverts this cross
pairing and then adds a triangular combination of the commuting half to
remove the first half's internal pairing. This explicitly gives
$`m/2`$ logical Pauli pairs. Products used as new generators retain their
actual Pauli phases; when such a product squares to minus identity, its
required factor of $`\pm i`$ makes the chosen logical generator Hermitian.

Extend those pairs to a symplectic basis of all physical Pauli labels.
Clifford elimination realizes this fixed target-coordinate change with
$`O(m^2)`$ elementary Cliffords on the two banks, no T gates and no
additional work. Literal signs of the canonical generators can be fixed
by a Pauli correction. For any total carrier width $`w\ge2m`$, the
generated complex operator algebra is therefore exactly

```math
\mathrm{Alg}_{\mathbb C}\{P_x:x\in\mathbb F_2^p\}
 \cong M_{2^{m/2}}(\mathbb C)\otimes I_{2^{w-m/2}}.      \qquad(10)
```

Its center consists only of scalars. This is a canonical coordinate
description of the target observables. If a source family $`R_x`$ is
also Pauli with the same independent commutation pattern, symplectic
extension likewise supplies an exact fixed Clifford matching every
literal generator. For general Hermitian **Clifford** $`R_x`$, the same
relations supply an abstract unitary equivalence of operator algebras;
they do not supply a native decoder. A non-Pauli Clifford cannot be
turned into a Pauli by Clifford conjugation.

There is no universal cheap replacement hidden in that distinction.
For $`k=m/2`$, take any homogeneous cubic Boolean f on k logical qubits
and $`D_f=\mathrm{diag}((-1)^f)`$. The family
$`D_fX_iD_f^\dagger,Z_i`$ consists of Hermitian Clifford involutions with
the full canonical commutation pattern. Its first generators are X
times quadratic CZ phases. The $`2^{\binom{k}{3}}`$ choices give distinct
families separated by two in maximum generator operator norm: some
Boolean derivative of the difference of two cubics is nonzero. A fixed
binary change of generator basis gives the same example for (8)'s Gram
matrix; the separation remains two since corresponding monomial operators
have identical translations and differ only by signs.

On W physical wires, a circuit with at most t T gates has at most
$`2^{O(W^2+Wt)}`$ possibilities, after absorbing its Clifford intervals
into t signed Pauli rotation axes and one final Clifford. A circuit with
d native T bands has at most $`2^{O((d+1)W^2)}`$ possibilities. These counts
already permit arbitrary intermediate action on dirty helpers. For every
fixed polynomial allowance $`W=O(m\,\mathrm{poly}(p))`$, they cannot
cover the cubic families with both $`t=O(m\,\mathrm{poly}(p))`$ and
sublogarithmic d, even at uniform generator error below one. This is an
existential obstruction to an algebra-only decoder theorem, not a lower
bound for this specific sparse tree source. A successful source decoder
must exploit additional structure and pay its actual word and inverse.

The [fixtures](../../tests/test_interleaved_tree_probe.py) compare (5)
against the complete native word with dirty inputs, retain the literal
root sign, and independently check the label kernel and canonical-pair
reduction. The analytic rank argument is (7)--(9), not a finite-size fit.

<a id="one-query-orbit-gaps"></a>

## 3. One-query gaps without a Clifford-probe premise

The following screens allow an arbitrary Hermitian involution R and
arbitrary fixed unitary walls and decoder, including operations on every
dirty helper. Put

```math
\gamma=\sqrt{2-\sqrt2}.
```

All errors below are operator norms on the complete work space. Tensoring
an arbitrary reference changes neither the conclusions nor their constants.
These are restrictions on a specified source orbit, not general native
T-depth lower bounds.

### A majority eigenspace forces a common probe eigenvalue

Let L be any unitary on a finite-dimensional space of dimension n, with
an eigenspace of dimension strictly greater than n/2. Let R be any
Hermitian involution. At least one of R's two eigenspaces has dimension
at least n/2. It intersects that eigenspace of L nontrivially. For a
unit vector v in the intersection, write

```math
Lv=\lambda v,\qquad Rv=rv,\qquad r\in\{1,-1\}.
```

Then $`LRL^\dagger v=rv`$ too. After any fixed unitary D, the two
images $`DRD^\dagger`$ and $`DLRL^\dagger D^\dagger`$ therefore
share a unit eigenvector $`\psi=Dv`$ with the same eigenvalue r.

Suppose anticommuting Hermitian involutions P and Q approximate those
images, each within delta. On that vector,

```math
r\langle P\rangle_\psi,\ r\langle Q\rangle_\psi
 \ge1-\delta^2/2.
```

Anticommutation implies
$`\langle P\rangle_\psi^2+\langle Q\rangle_\psi^2\le1`$:
every real combination $`aP+bQ`$ has norm $`\sqrt{a^2+b^2}`$.
If delta were smaller than gamma, the two displayed expectations would
have squared sum greater than one. Thus

```math
\max\{\|DRD^\dagger-P\|,
       \|DLRL^\dagger D^\dagger-Q\|\}\ge\gamma.       \qquad(11)
```

This constant is sharp for the stated lemma. Take $`L=I`$, $`D=I`$,
$`R=H=(X+Z)/\sqrt2`$, and targets X,Z. Both distances are exactly
gamma. Neither Clifford structure nor a trace condition on R was used.

### Application to one linear XOR query

Consider a single queried orbit
$`S_x=DU_xRU_x^\dagger D^\dagger`$, where
$`U_x=FQ_xG`$ and F,G,D are fixed unitaries. Let

```math
Q_x|u,v,h\rangle=|u,v+M_xu,h\rangle
```

be a phase-free linear XOR shear; the binary matrix $`M_x`$ may have
arbitrary address dependence. For any two addresses, $`Q_xQ_y`$ is
the same kind of involution with matrix $`\Delta=M_x+M_y`$.
Its trace counts fixed basis columns. If n is the full dimension,
including all spectator helpers, then

```math
\mathrm{Tr}(Q_xQ_y)=n\,2^{-\mathrm{rank}\Delta}>0.
```

Its +1 eigenspace consequently has dimension
$`(n+\mathrm{Tr}(Q_xQ_y))/2>n/2`$. The relative unitary
$`U_x^\dagger U_y=G^\dagger Q_xQ_yG`$ has the same spectrum.
Equation (11) applies to every anticommuting target pair. For (5),
$`P_0`$ and $`P_{e_1}`$ anticommute at every $`p\ge1`$, so

```math
\max_x\|S_x-P_x\|\ge\gamma.                         \qquad(12)
```

In particular (12) excludes a single minimal-tree XOR query before a
fixed readout, even when the first probe R is genuinely non-Clifford.
The trace includes every arbitrary bank input; it is not a statement
restricted to the fixed columns used to count it.

### Affine offsets retain a separate address-covariance gap

Address-dependent offsets can remove that majority eigenspace. For a
single affine shear, however, assume the actual parameters are affine:

```math
Q_x|u,v\rangle=|u,v+M(x)u+h(x)\rangle,\qquad
M(x)=M_0+\sum_i x_iM_i,\quad
h(x)=h_0+\sum_i x_ih_i.
```

The shears commute, and $`V_x=Q_0Q_x`$ obeys
$`V_xV_z=V_{x+z}`$. More generally the next argument only requires
a projective representation of the additive address group: scalar
multiplication phases cancel in its adjoint action.
Absorb $`Q_0`$ and all fixed walls into the probe and decoder. The
source orbit then has fixed unitaries $`W_z`$ satisfying

```math
S_{x+z}=W_zS_xW_z^\dagger
\quad\hbox{for every }x,z.                          \qquad(13)
```

The target algebra does not have this covariance. For $`p\ge2`$,
equation (7) gives

```math
K_p(0,e_1)=1,\qquad K_p(e_2,e_1+e_2)=0.             \qquad(14)
```

If $`\max_x\|S_x-P_x\|\le\delta`$, (13) implies that
$`W_{e_2}P_0W_{e_2}^\dagger`$ and
$`W_{e_2}P_{e_1}W_{e_2}^\dagger`$ are each within $`2\delta`$
of $`P_{e_2}`$ and $`P_{e_1+e_2}`$, respectively. Call the first
pair E,F; it anticommutes exactly. The second pair commutes, so some
nonzero simultaneous eigenspace exists. Choose a unit vector psi in
it, with the respective eigenvalues r,s. The two norm comparisons give

```math
r\langle E\rangle_\psi,\ s\langle F\rangle_\psi
 \ge1-(2\delta)^2/2,\qquad
\langle E\rangle_\psi^2+\langle F\rangle_\psi^2\le1.
```

For $`2\delta\lt\gamma`$ these inequalities contradict each other.
No particular joint sign sector is assumed to be present. Therefore

```math
\max_x\|S_x-P_x\|\ge\gamma/2,\qquad p\ge2.           \qquad(15)
```

This covers arbitrary non-Clifford probes and fixed decoders for every
single affine-query orbit, including offsets. It does not cover
non-affine offsets or arbitrary address-dependent native words.

The offset boundary is real for a single address bit. Take
$`V_0=I,V_1=X`$ and $`R=H`$. Then

```math
R_0=(X+Z)/\sqrt2,\qquad R_1=(X-Z)/\sqrt2,
\qquad R_0R_1=-R_1R_0.
```

The relative X has no majority eigenspace. The literal native word
$`D=SH T^\dagger H S^\dagger`$ satisfies
$`DR_0D^\dagger=Z`$ and $`DR_1D^\dagger=X`$.
Here $`S^\dagger=S^3`$ if only positive S is emitted, and the actual
adjoint word is used. This example only marks the boundary of (11);
it does not evade the full-family obstruction (15) for $`p\ge2`$.

### Two genuinely interleaved queries can leave these premises

On one c wire and one t wire, set

```math
L=\mathrm{CX}(t;c)\mathrm{CX}(c;t).
```

This is an order-three permutation, with spectrum
$`\{1,1,\zeta,\zeta^2\}`$ for $`\zeta=e^{2\pi i/3}`$.
Its largest eigenspace has dimension exactly half the total, and it is
not a Hermitian involution. Literal conjugation gives

```math
Y_c\longmapsto Z_cY_t\longmapsto X_cY_t
     \longmapsto Y_c,                               \qquad(16)
```

with positive signs. These three observables anticommute pairwise.
Thus $`C_0=I,C_1=L`$ supplies an anticommuting pair from a fixed probe
using two opposite address-controlled XOR stages. The controlled L
word is exactly a Toffoli from x,c to t followed by a Toffoli from
x,t to c. The retained seven-T, four-band Toffoli word pays it with
14 T gates, at most eight T bands, no extra wires and no initialized
inputs. Its actual inverse reverses the stages; on the x=1 fiber,
$`L^\dagger=L^2\ne L`$.

This is a boundary witness for the one-query premises, not a uniform
source for (5). Changing R alone or adding affine offsets cannot repair
one queried orbit. Genuinely later address participation must instead
change its orbit structure, and a proposed all-size construction still
needs the complete native source, decoder, inverse and resource ledger.
For fixed Clifford R, the separate
[Clifford-orbit rank theorem](CLIFFORD_ORBIT_QUERY_RANK.md)
also bounds words with multiple interleaved queries.

The [one-query fixtures](../../tests/test_query_orbit_spectral_gaps.py)
check full-space fixed-point multiplicities, a non-Clifford probe with an
arbitrary fixed decoder, sharpness of (11), affine/projective covariance,
the target mismatch (14), and both boundary examples. The controlled
order-three word and its actual inverse are checked on every native
basis column.

<a id="nonclifford-seed-decoder-condition"></a>

## 4. A non-Clifford seed also constrains the actual decoder

Suppose $`R=D_0KD_0^\dagger`$ is genuinely non-Clifford up to scalar,
where $`D_0`$ is third-level and K is a Hermitian Clifford. Let
$`C_x`$ be any queried Clifford family, and let E be the actual fixed
decoder on the complete work space. If even one selected target has
$`E^\dagger P_xE`$ Clifford, put

```math
T_x=C_x^\dagger E^\dagger P_xE C_x,\qquad
\|EC_xRC_x^\dagger E^\dagger-P_x\|=\|R-T_x\|.
```

Here $`T_x`$ is Clifford. At error strictly below $`\sqrt2/4`$, the
[third-level rounding lemma](NATIVE_CLIFFORD_TRANSPORT.md#1-a-robust-exactness-lemma-for-every-third-level-source)
forces $`D_0KD_0^\dagger`$ to equal $`T_x`$ up to scalar, contradicting
the premise on R. This proves an error floor $`\sqrt2/4`$ for that
selected target, independently of the number of queried Clifford stages.

In particular, a decoder with $`E^\dagger`$ third-level is excluded.
A single weighted-cubic decoder in fixed Clifford coordinates satisfies
that adjoint premise: negating its phase polynomial preserves the same
weighted-cubic form. No inverse-closure assertion is made for a general
third-level unitary. Leaving this hierarchy class globally is also
insufficient: every selected pullback $`E^\dagger P_xE`$ must avoid the
Clifford class at the stated accuracy. A viable source must therefore
change both its intermediate reflection and the selected action of its
paid decoder. All statements concern the actual complete operators,
including arbitrary dirty/reference inputs.
