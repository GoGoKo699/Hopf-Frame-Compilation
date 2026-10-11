# Native realification, dyadic completion, and tree-pattern restrictions

[Complete-frame contract](../../docs/HOPF_INTERFACE.md) · [One-clean grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) · [Endpoint](../../docs/OPEN_PROBLEM.md)

Every native compiler for a real target has an exactly real native
implementation with one additional clean or dirty qubit, at most twice
the T-count, and the same full-isometry error. In contrast, exact
ring-valued orthogonal matrices that preserve the frame's mandatory
marker supports have a constant approximation gap. The distinction is
between the enlarged physical completion and an exactly preserved
logical sparsity pattern. Every complete frame also admits a determinant-one
dyadic unitary approximation on exactly two clean qubits with linear
denominator height, as proved in Section 3. That arithmetic result does
not provide a linear T-count factorization.

## 1. Realification preserves the full initialized-isometry error

For a complex, possibly rectangular matrix B, define

```math
\mathcal R(B)=
\begin{pmatrix}\mathrm{Re}B&-\mathrm{Im}B\\
\mathrm{Im}B&\mathrm{Re}B\end{pmatrix}.
```

This real-algebra homomorphism preserves products and adjoints. With
$`F=HS=2^{-1/2}\begin{pmatrix}1&i\\1&-i\end{pmatrix}`$,

```math
(F\otimes I)\mathcal R(B)(F^\dagger\otimes I)
=B\oplus\overline B,
\qquad \|\mathcal R(B)\|=\|B\|.
```

The identities on the two sides have the appropriate sizes for
rectangular B. Its singular values are therefore repeated twice by
realification, with no quotient by scalar phase.

**Theorem.** Suppose a finite native word V uses a clean qubits, b
arbitrary dirty qubits and t T or inverse-T gates, and satisfies

```math
\|VJ_a-J_a(W\otimes I_b)\|\le\varepsilon
```

for a real target W. There is an exactly real native word implementing
$`\mathcal R(V)`$, with the following alternatives:

| Additional register | Clean qubits | Dirty qubits | T-count | Full error |
|---|---:|---:|---:|---:|
| Arbitrary qubit | a | b+1 | at most 2t | at most epsilon |
| Initialized qubit | a+1 | b | at most 2t | at most epsilon |

Both variants include clean leakage, dirty return, arbitrary reference
correlations and the literal target phase.

*Proof of the error claim.* Put $`S_0=W\otimes I_b`$ and
$`E=VJ_a-J_aS_0`$. Since $`J_a`$ and $`S_0`$ are real,

```math
\mathcal R(V)(I_2\otimes J_a)
-(I_2\otimes J_a)(I_2\otimes S_0)=\mathcal R(E).
```

The norm is $`\|E\|`$ on the entire added-register space, including
complex superpositions. Consequently the additional qubit may be dirty
and entangled with any other input or inaccessible reference. Tensoring
with a reference preserves the operator norm. Restricting its input to
$`|0\rangle`$ instead gives the clean variant, with its return error
already included. ∎

### Exact native gate words

The extra realification register is r. The real gates H and CNOT act
unchanged on their original wires. Use the literal real-rotation convention

```math
R_y(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix}.
```

Then

```math
\begin{aligned}
\mathcal R(S_x)
&=C_x[R_y(\pi/2)_r]
=\mathrm{CNOT}_{x,r}\mathrm{CZ}_{x,r},\\
\mathcal R(T_x)
&=C_x[R_y(\pi/4)_r]
=C_x[H_r]\mathrm{CZ}_{x,r}.
\end{aligned}
```

These follow from $`R_y(\pi/2)=XZ`$ and
$`R_y(\pi/4)=HZ`$. Put $`A=SHTHS^\dagger`$. Direct multiplication gives
$`AZA^\dagger=H`$, hence

```math
C_x[H_r]=(I_x\otimes A_r)\mathrm{CZ}_{x,r}
                 (I_x\otimes A_r^\dagger).
```

This word contains exactly one T and one inverse T; CZ is
H–CNOT–H. Replace inverse gates by the actual inverse words. Products
then implement $`\mathcal R(V)`$ exactly, using only H, S, CNOT, T
and their native inverses. SWAPs can place r anywhere. No arbitrary
rotation primitive or supplied state enters this construction.

**Endpoint specialization.** The
[one-clean grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit)
leaves the second supplied clean flag available as r. At
$`L=N`$ it therefore gives an exactly real native completion with

```math
a=2,\qquad b=N+n+7,\qquad
T=O\bigl(N[1+\log_2^*(n+2)]\bigr),\qquad
\|\mathcal R(V)J_2-J_2(W\otimes I_b)\|\le2^{-N},
```

after the fixed register reordering. Thus requiring an exactly real
enlarged native completion does not change the established asymptotic
upper bound. Realification does not by itself reduce that bound to
$`O(N)`$.

## 2. Exact ring-valued marker supports have a constant gap

Let $`R=\mathbb Z[1/\sqrt2]`$. Every exact native Clifford+T matrix has
entries in $`R[i]`$, by inspection of the generators and closure under
matrix multiplication. An exactly returned real logical unitary
therefore has entries in R, even with clean or arbitrary dirty work.

**Circle lemma.** The only $`x,y\in R`$ satisfying $`x^2+y^2=1`$ are

```math
(\pm1,0),\quad(0,\pm1),\quad
(\pm2^{-1/2},\pm2^{-1/2}),
```

with all independent signs in the last pair.

*Proof.* Write $`x=(a+b\sqrt2)/2^k`$ and
$`y=(c+d\sqrt2)/2^k`$, with integer coefficients and $`k\ge0`$.
Rational and irrational parts give

```math
a^2+c^2+2b^2+2d^2=4^k,\qquad ab+cd=0.
```

For $`k\ge1`$, a and c have the same parity. If both were odd,
the first equation modulo four would force b and d to have opposite
parity, contradicting the second equation modulo two. Thus a and c
are even. If $`k\ge2`$, set $`a=2a'`$, $`c=2c'`$. Dividing the
first equation by two shows that b and d have the same parity. If both
were odd, reduction modulo four would force $`a',c'`$ to have
opposite parity, again contradicting $`a'b+c'd=0`$ modulo two.
All four coefficients are even, so divide by two and descend in k.
At $`k=0`$ only the axis points occur; $`k=1`$ adds exactly the four
diagonal points. ∎

For each tree node $`(d,p)`$, let

```math
I_{d,p}=\{p2^{n-d},\ldots,(p+1)2^{n-d}-1\},\qquad
z_{d,p}=(2p+1)2^{n-d-1}.
```

The prescribed frame's marker column $`z_{d,p}`$ is supported in
$`I_{d,p}`$. No support restriction is imposed on column zero.

**Exact-pattern theorem.** Every real orthogonal Q with these marker
supports is a Hopf frame times a possible sign on its root column.
If Q has entries in R, every frame angle is a multiple of
$`\pi/4`$ modulo $`2\pi`$. Consequently this exact-pattern class
has at most $`2\cdot8^{N-1}`$ ring-valued matrices.

*Proof.* A bottom marker is a unit vector on coordinates
$`2p,2p+1`$. Choose the real plane rotation sending the second basis
vector to that marker. Multiplying Q by its transpose sends the marker
to $`|2p+1\rangle`$; orthogonality makes every other column zero in
that odd row. Remove all odd marker columns and rows. The remaining
matrix is orthogonal with the smaller tree's marker supports. Repeating
this operation leaves a root coordinate equal to +1 or −1 and proves
the first assertion.

If Q is ring-valued, the circle lemma makes each bottom stripping
rotation ring-valued with angle a multiple of $`\pi/4`$. The residual
matrix remains ring-valued, so induction proves the second assertion.
For an original Hopf frame the bottom marker is exactly

```math
-\sin\theta_{n-1,p}|2p\rangle+
\cos\theta_{n-1,p}|2p+1\rangle.
```

This also restricts each of its prescribed angles. No ancestor amplitude is divided out, and
singular angles and all signs are included. ∎

**Uniform approximation gap.** Set only
$`\theta_{n-1,0}=\pi/8`$ nonzero. Its marker column 1 has distance
at least

```math
2\sin(\pi/16)>0.39
```

from every ring-valued unit vector supported on coordinates 0 and 1.
Hence the operator distance from every ring-valued Q in the stated
exact-pattern class has the same lower bound, independently of n.
Additional denominator bits or exactly returned work cannot remove it.

The support hypothesis is essential: small off-pattern entries and
approximate work return lie outside this theorem. In particular, the
realified physical completion of Section 1 is not required to preserve
the logical marker pattern.

## 3. Every frame has a two-clean linear-height dyadic completion

Use **dyadic height** H to mean that every matrix entry belongs to
$`2^{-H}\mathbb Z[i]`$; this is a power-of-two denominator exponent,
not the least square-root-of-two exponent of a Pauli-transfer matrix.

**Theorem.** Let $`n\ge3`$, $`N=2^n`$, and let W be any prescribed
complete real Hopf frame with computable angles, including zero and
singular branches. There is an effective construction of an exact
determinant-one unitary V on $`n+2`$ qubits such that

```math
V\in 2^{-H}M_{4N}(\mathbb Z[i]),\qquad
H\le1085N,\qquad
\|VJ_2-J_2W\|\lt2^{-N}.
```

Every dirty wire may remain untouched. The same error therefore includes
arbitrary dirty inputs and references, with literal phase. There is a
finite exact native synthesis of V on these same wires, but **no
$`O(N)`$ T-count bound is claimed**. The construction below proves the
height directly; it does not infer height from generic column elimination.

### Regular laminar charts: integer columns and an empty quadrant

First consider a direct sum of binary laminar frames on at most N modes,
each of depth at most n. For an internal node v, write its branch
coefficients as $`c_v,s_v`$ and its local root-to-node scaling amplitude
as $`g_v`$. Suppose a certified integer $`B\ge0`$ satisfies

```math
\mu=\min_v |g_v|^3|c_vs_v|\ge2^{-B}>0.
```

For a requested error $`\delta=2^{-L}`$, $`L\ge1`$, set

```math
k=L+n+B+10,\qquad K=3k+2L+12=5L+3n+3B+42.
```

In each component let q be its normalized root column. Certified
enclosures produce an integer vector z with
$`\|2^{-k}z-q\|\le\sqrt N\,2^{-k}`$. Use z as the integer root
column. For a split with child restrictions $`z_L,z_R`$, put

```math
a_v=\|z_L\|^2,\quad b_v=\|z_R\|^2,\qquad
u_v=\mathop{\rm sign}(g_vc_vs_v)(-b_vz_L,a_vz_R).
```

These integer columns are exactly orthogonal: disjoint nodes have
disjoint support, and an ancestor column restricts to a scalar multiple
of z on every descendant, whereas each displayed detail is orthogonal
to z there. For the unrounded root, the complete-column convention gives

```math
(-\|q_R\|^2q_L,\|q_L\|^2q_R)
=g_v^3c_vs_v\,W e_v.
```

Thus the signs recover the prescribed marker columns, not an arbitrary
completion of the first column. On vectors of norm at most two this
cubic map is 12-Lipschitz. Normalizing its nonzero values changes each
column by at most $`24\sqrt N\,2^{-k}/\mu`$. Consequently the
orthogonal matrix Q obtained by normalizing all integer columns satisfies

```math
\|Q-W\|\le\|Q-W\|_F
\le24N2^{-k}/\mu\le(24/1024)\delta\lt\delta/32.
```

This proof uses only laminar child restrictions and disjoint components;
unary nodes may be deleted, and isolated coordinates are identity blocks.
There is no equal-size-child requirement.

Each integer column has norm at most $`2^{3k+2}`$. With
$`\eta=\delta^2/64`$, choose

```math
\ell_j=\left\lfloor\frac{(1-\eta/2)2^K}{\|u_j\|}\right\rfloor,
\qquad A_j=2^{-K}\ell_j u_j,\qquad
\alpha_j=2^{-K}\ell_j\|u_j\|.
```

The floor is computable by integer/rational square comparisons. The
rounding loss is at most $`2^{3k+2-K}=\eta/16`$, so
$`1-\eta\le\alpha_j\le1-\eta/2`$. The Gram matrix of A is
exactly diagonal. Write each positive integer
$`2^{2K}-\ell_j^2\|u_j\|^2`$ as four integer squares, found by finite
exhaustive search. Pairing these squares into two Gaussian integers
gives diagonal dyadic matrices C,D such that

```math
X=\begin{bmatrix}A\\C\\D\\0\end{bmatrix},\qquad
X^\dagger X=I,\qquad
(X-J_2Q)^\dagger(X-J_2Q)=\mathrm{diag}(2-2\alpha_j).
```

Let E include the logical modes in the unused fourth quadrant and put

```math
R=I-(E-X)(E-X)^\dagger,\qquad
P=X_{\rm clean,1}X_{\rm clean,2},\qquad V=RP.
```

Here $`(E-X)^\dagger(E-X)=2I`$, hence R is an exact Hermitian
unitary, $`RE=X`$, and $`VJ_2=X`$. Its height is at most $`2K`$.
R has exactly N negative eigenvalues and P is an even permutation, so
$`\det V=1`$. The full-output error is

```math
\|VJ_2-J_2W\|\le\sqrt{2\eta}+\delta/32\lt\delta,
\qquad H\le10L+6n+6B+84.
```

No matrix square root, odd-prime denominator, extra clean quadrant, or
projective phase is hidden in this completion.

### Certified normalization and weighted components

Put $`\varepsilon=2^{-N}`$ and $`M=2^{N+n+8}`$. Replace every angle
by an odd grid point $`(2j+1)\pi/M`$ whose sine/cosine pair is at
Euclidean distance less than $`4/M`$. Such a point exists, and finite
dovetailing over strict certified enclosures finds one without deciding
an equality. The complete operator changes by at most
$`4N/M=\varepsilon/64`$. Grid branches are nonzero, have magnitude
at least $`2/M`$, and have no equal-magnitude ties.

Signed input/output permutations normalize each split to

```math
2/M\le s_v\le1/\sqrt2,\qquad c_v\ge1/\sqrt2.
```

For clarity, these are full-frame identities. At a node with child
frames A,B and root-pair rotation G, quadrant normalization multiplies
output halves by $`a=\mathop{\rm sign}(c)`$ and
$`b=\mathop{\rm sign}(s)`$. On the input, use signs a on left detail
columns, b on right detail columns, 1 on the root column and ab on the
marker. To exchange the larger branch, swap the output halves and the
corresponding descendant input columns, fixing both root inputs, then
put a minus sign on the marker. The new angle is $`\pi/2-\theta`$.
Child input corrections fix child roots and commute through the parent
rotation, so these identities compose over the whole frame.

View the computational labels as the tree with parent
$`p(j)=j-\mathrm{lowbit}(j)`$. Put

```math
\Lambda=N+2n+10,\qquad a=2^{-\Lambda}.
```

Grow components down this tree. At a vertex let x be the product of
internal path sines since its component root. Expand **all** outgoing
edges if $`x>a`$, or cut **all** outgoing edges if $`x\lt2a`$; the tests
overlap and at least one terminates. Every child past a cut begins a new
component with x equal to one. Each included nonroot vertex has path
product greater than $`a(2/M)`$, and each boundary vertex has path
product less than $`2a`$.

The all-outgoing rule is essential. It ensures that no later internal
edge moves a boundary axis after its cut. An edge-by-edge stopping rule
would not justify the following completed-component conjugation.

### Seven regular factors, with the reordering error charged

Let $`U_0`$ be the direct sum of all internal component words. Moving
these internal rotations left rewrites the normalized rounded frame as
$`U_0`$ times the original cut word, with a cut from boundary l of
component C to child-component root r acting in the plane
$`\mathrm{span}\{U_C^\dagger e_l,e_r\}`$. The original prefix
conjugation equals the completed $`U_C`$ conjugation because all later
internal edges fix l.

Outgoing cuts with distinct boundary vertices have orthogonal parent
axes; those sharing a boundary retain their order. The only reordered
cuts that may fail to commute are an incoming and an outgoing cut of
one component. Their sole plane overlap is
$`\gamma=\langle e_l,U_C e_{\rm root(C)}\rangle`$, bounded by
the internal path-sine product, hence by $`2a`$.

For planar rotations with nonnegative cosines,

```math
\|[R_A,R_B]\|
\le|\gamma|(|s_A|+1-c_A)(|s_B|+1-c_B)
\le4|s_As_B\gamma|\lt8a.
```

This follows by writing each rotation as
$`I+sA+(c-1)P_A`$; the four resulting commutators of plane generators
and projections have norm at most $`|\gamma|`$. Stable partition by
parent-component parity needs at most $`N^2`$ swaps. Consequently,

```math
\|W_+-U_0F_{\rm even}F_{\rm odd}\|
\lt8N^2a=\varepsilon/128.
```

For either parity, let $`U_{\rm parity}`$ be the direct sum of internal
words of parent components of that parity. Let $`K_{\rm parity}`$ be
their original cut word. Its stars have disjoint supports, centered at
boundary vertices with their cut-child roots. Thus

```math
F_{\rm parity}=U_{\rm parity}^\dagger K_{\rm parity}U_{\rm parity},
```

giving the seven factors

```math
U_0U_{\rm even}^\dagger K_{\rm even}U_{\rm even}
   U_{\rm odd}^\dagger K_{\rm odd}U_{\rm odd}.
```

An internal component is the original binary decision tree with branches
pruned and unary nodes deleted. Its laminar depth is at most n. Every
ancestor scaling coefficient has at most n cosine factors in addition
to its internal path sines, so

```math
|g_v|\ge2^{-\Lambda-\log_2 M-n/2},\qquad
|g_v|^3|c_vs_v|\ge2^{-B},\qquad
B=3\Lambda+4\log_2M+2n+1=7N+12n+63.
```

For a cut star, the continuing-center path has at most n cosines; its
chart satisfies the smaller bound $`B_{\rm star}=\log_2M+2n+1`$.
The common B therefore applies to all seven factors.

### Height, determinant, and the remaining native gap

Apply the regular-chart construction to each factor with
$`\delta=\varepsilon/32`$, so $`L=N+5`$. Each exact physical factor
has determinant one and height at most

```math
H_{\rm factor}\le10(N+5)+6n+6B+84=52N+78n+512.
```

Use actual adjoints for inverse factors. For any full unitary A with
$`\|AJ-JU\|\le\delta`$,

```math
A^\dagger J-JU^\dagger
=A^\dagger(JU-AJ)U^\dagger.
```

Thus inverse errors are unchanged, and the usual initialized-isometry
hybrid composes all seven factors without reinitializing their clean
qubits. Their heights add. Since $`n\le N/2`$ and $`N\ge8`$,

```math
H\le364N+546n+3584\le1085N,
\qquad
\frac{\varepsilon}{64}+\frac{\varepsilon}{128}
 +\frac{7\varepsilon}{32}\lt\frac{\varepsilon}{4}.
```

Undo the exact signed input/output permutations. Extended identically
to all four clean sectors, each has determinant equal to the fourth
power of its logical determinant, hence one; its entries are integers.
This preserves the height, determinant, and full error, and proves the
theorem for the original signed, possibly singular tuple.

The matrix belongs to the native ring. The determinant-one case of
[Giles--Selinger, Corollary 2](https://arxiv.org/pdf/1212.0506)
provides finite exact Clifford+T synthesis without an additional ancilla
on $`n+2\ge5`$ wires. If a scalar $`e^{i\pi/4}I`$ appears in that
gate convention, the literal Clifford word $`(HS)^3`$ supplies it.
This establishes effective native existence only. A height bound does
not control T-count, and a generic exact synthesis may enlarge intermediate
denominators. Factoring these particular integer-wavelet completions
with $`O(N)`$ T gates remains a separate unsupported step.

## 4. Verification

The [native-representation tests](../../tests/test_native_representation_limits.py) verify the literal gate words
in exact $`\mathbb Q(\sqrt2,i)`$ arithmetic, check the circle equations
over bounded integer ranges, and exercise the rectangular norm and
full-isometry identities. Small frame checks confirm the marker-column
gap. The all-size conclusions follow from the proofs above, independently
of these finite checks.

The same tests check exact integer-wavelet orthogonality and the
empty-quadrant reflection completion, including its determinant and
dyadic denominator. They independently exercise the completed-component
cut identity and its parity factorization on full small matrices. They
do not run a general exact synthesis or certify a linear T-count.
