# Canonical stars and their sparse resolvent

A complete tree frame can change each local star to a determinant-one plane
rotation, paying one disjoint input correction. The resulting frame has a
Cayley representation with a lower-triangular denominator. Its complete
local scattering oracle admits an $`O(N)`$-T native implementation with two
clean qubits. Implementing the feedback inverse remains a separate problem.

Use $`N=2^n`$ and the computational representative tree
$`p(j)=j-\mathrm{lowbit}(j)`$, for $`j>0`$. Products act rightmost first.
The local-star results assume a real nonnegative carrier component
$`c_v`$ in the root column; normalized real Hopf frames satisfy this
condition after the
[effective signed normalization](PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization).
Zero child components and $`c_v=0`$ are included. The algebra
also allows complex child components, retaining their literal phases.

This is a different representative from the
[fixed-tree phase-gauge Cayley reduction](TREE_CAYLEY_REDUCTION.md): it
uses whole-star root columns and an explicitly priced input correction.
The formulas below identify its native components and the unresolved
inverse without importing the other Cayley family's compiler condition.
The [coupled resolvent](COUPLED_TREE_RESOLVENT.md) instead tracks a graded
coarse/fine defect; no separate-power loader or feedback compiler is
imported from that construction.

## 1. Same-root-column transfer

At vertex $`v`$, let $`S_v`$ be any determinant-one unitary on the carrier
and its immediate child labels. The complete subtree frame satisfies

```math
W_v=\mathrm{diag}(1,W_{w_1},\ldots,W_{w_d})S_v.
```

Choose another determinant-one $`S'_v`$ with exactly the same first column.
Then $`S_v'^\dagger S_v`$ fixes the carrier pointwise and restricts to a
determinant-one child-space block $`C_v`$. Induction gives the exact
complete-operator identity

```math
\boxed{W=W'C_{\rm in},\qquad C_{\rm in}=\bigoplus_v C_v.}
```

Indeed, each correction below a child fixes that child's carrier, and
therefore commutes with the parent star. Moving these corrections to the
input leaves exactly the displayed direct sum. Immediate-child label
sets are globally disjoint, even though the original subtree supports
overlap. The degree count is

```math
\sum_v d(v)^2
=n^2+\sum_{d=0}^{n-1}\frac{N d^2}{2^{d+1}}
=3N-2n-3\le3N.
```

Thus the [disjoint-block compiler](DISJOINT_BLOCK_COMPILER.md) implements
$`C_{\rm in}`$ once with $`O(N)`$ T gates, under its complete-return
contract. This statement prices the correction; it does not assert that
$`W'`$ itself has a linear-T compiler.

The transfer also permits a different ordering of the same outgoing
root-column components. For a real root column $`(c,u_1,\ldots,u_d)`$
with $`c>0`$, and any chosen child order, define

```math
r_j=\sqrt{c^2+\sum_{k>j}u_k^2},\qquad r_0=1,
\qquad \sin\alpha_j=u_j/r_{j-1},\quad \cos\alpha_j=r_j/r_{j-1}.
```

The ordered carrier-child rotations with these coefficients have that
same root column. This changes the local frame, without globally
permuting the child coordinates. The additional coefficient and depth
bounds for stable cheap-first ordering belong to the
[structured forest construction](STRUCTURED_FOREST_COMPILERS.md).

## 2. A canonical local frame and tree shift

Write the local root column as $`(c_v,u_v)^T`$, where
$`c_v\ge0`$ and $`\|u_v\|^2=1-c_v^2`$. Choose

```math
S'_v=\begin{pmatrix}
c_v&-u_v^\dagger\\
u_v&I-u_vu_v^\dagger/(1+c_v)
\end{pmatrix}.
```

On the carrier and the direction $`u_v`$ this is a determinant-one plane
rotation; it fixes the orthogonal child complement. If $`u_v=0`$ it is
identity. Only $`1+c_v`$ is inverted, never a small child component.

Set $`L_{j,p(j)}=(u_{p(j)})_j/(1+c_{p(j)})`$ and all other entries to
zero. Its columns have disjoint supports, so

```math
L^\dagger L=\mathrm{diag}(\lambda_v),\qquad
\lambda_v=\frac{1-c_v}{1+c_v},\qquad L^{n+1}=0.
```

Leaf vertices have $`c_v=1`$ and $`\lambda_v=0`$. The root is included
in this shift: $`L^n`$ need not vanish. This differs from a shift whose
indices contain only nonroot internal ports.

For external input $`x`$, output $`o=W'x`$, and incoming carrier $`y_v`$
at vertex $`v`$, with $`y_0=x_0`$, the local star equations imply the
**carrier identities**

```math
\boxed{y-x=L(y+o),\qquad o-y=-L^\dagger(y+x).}
```

For example,
$`y_j-x_j=L_{j,p(j)}(y_{p(j)}+o_{p(j)})`$ follows by substituting the
two blocks of $`S'_{p(j)}`$. The second identity follows from
$`(1+c_v)L_{\mathrm{children}(v),v}=u_v`$. Absent root and leaf terms
are zero, so these identities include every external input and output
coordinate.

## 3. Cayley identity and a triangular denominator

Define

```math
F=I-L,\qquad R=F^{-1}=\sum_{r=0}^n L^r,\qquad X=R^\dagger-R.
```

Then $`X^\dagger=-X`$ and the **Cayley identity** is

```math
\boxed{W'=(I-X)(I+X)^{-1}=2(I+X)^{-1}-I.}
```

To prove this, put $`w=y+x`$ in the carrier identities. Elimination gives

```math
(I-2L+LL^\dagger)w=2(I-L)x,\qquad o=(I-L^\dagger)w-x.
```

Both congruences hold exactly:

```math
FXF^\dagger=F^\dagger XF=L^\dagger-L.
```

The first makes $`I-2L+LL^\dagger=F(I+X)F^\dagger`$ and proves the
Cayley identity.
The second supplies the sparse orientation

```math
D=I+L^\dagger L,\qquad
M=F^\dagger(I+X)F=D-2L,
\qquad D_v=\frac{2}{1+c_v}\in[1,2].
```

In particular, $`D`$ is diagonal and $`M`$ is lower triangular. The
**triangular resolvent identity** is therefore

```math
\boxed{W'=2(I-L)(D-2L)^{-1}(I-L^\dagger)-I.}
```

There is no singular zero-angle chart or unspecified scalar phase:
$`L=0`$ gives $`W'=I`$. On a single edge, the Cayley identity is the usual
half-angle formula.

Let $`A_{j,p(j)}=(u_{p(j)})_j`$. The **carrier factorization** of the same
denominator is

```math
2LD^{-1}=A,\qquad D-2L=(I-A)D,
\qquad (D-2L)^{-1}=D^{-1}(I-A)^{-1}.
```

Thus the triangular resolvent still contains the descendant scattering
inverse. A native implementation of that inverse cannot be inferred from
triangularity.
For comparison, the [complete scattering account](ENDPOINT_TREE_TRANSPORT.md)
separates local oracle construction from external feedback recovery.

## 4. Uniform conditioning

Disjoint column supports and nilpotency imply

```math
\|L\|=\max_v\sqrt{\frac{1-c_v}{1+c_v}}\le1,\qquad
\|A\|=\max_v\sqrt{1-c_v^2}\le1,
\qquad \|R\|\le n+1,\quad \|X\|\le2n.
```

The smallest singular value of $`I+X`$ is at least one, because
$`X`$ is skew-Hermitian. Its condition number is at most
$`\sqrt{1+4n^2}`$. The carrier factorization gives
$`\|M^{-1}\|\le n+1`$ and $`\|M\|\le4`$, hence a condition number
at most $`4(n+1)`$. These bounds include zero coefficients.

If every $`c_v\ge c_*>0`$, then
$`\|L\|\le\sqrt{(1-c_*)/(1+c_*)}\lt1`$ and the corresponding geometric
series bound improves $`\|R\|`$. Choosing a different detail completion
with the same root column cannot change $`c_v`$ or this contraction ratio.
These conditioning bounds price no inverse or matrix function by themselves.

## 5. Native compilation of the complete local oracle

Consider the isometry on two banks of $`N`$ modes,

```math
V=\begin{bmatrix}I\\L\end{bmatrix}D^{-1/2},\qquad V^\dagger V=I.
```

Column $`v`$ is supported on the upper mode $`v`$ and the lower modes of
its children. These supports are disjoint, have total size $`2N-1`$,
and leave the lower root unused. Its upper-anchor component is the real
number $`1/\sqrt{D_v}\ge1/\sqrt2`$.

For this column $`z`$ and its upper anchor $`e`$, set
$`w=(e+z)/\|e+z\|`$ and $`H_w=I-2ww^\dagger`$. Then
$`H_wH_e e=z`$ literally. Each completed reflection fixes every coordinate
outside the original column support, including any temporary zero padding.
Their direct sum defines a particular complete unitary $`U_V`$, rather
than an unspecified state-preparation extension.

### One-bin support bound

Leaf columns have $`z=e`$ and $`H_wH_e=I`$ exactly. They can be omitted
using topology alone, with no test for an unknown zero angle. For a
nonleaf of degree $`d\ge1`$, its next-power-of-two support size obeys

```math
B_v=2^{\lceil\log_2(d+1)\rceil}\le2d,
\qquad \sum_{v\text{ nonleaf}}B_v\le2(N-1)\lt2N.
```

Descending-size packing therefore fits every padded block as an aligned
interval in the single $`2N`$-mode space. Leaf upper coordinates may be
padding, since the completed ideal reflections fix them. There is no
additional source bank or initialized address register.

The [disjoint-reflection source theorem](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)
then implements $`U_V`$ in $`O(N)`$ T gates. The first supplied clean
bit is part of the $`n+1`$-qubit mode address and may be arbitrary during
the operation. The second is the sole initialized compiler flag. Exact
coarse words, packing permutations, and anchor signs use no further
initialized bit and return their borrowed work exactly.

### Full mode-space and dirty-return contract

Let $`J`$ initialize only the compiler flag. Its domain includes all
$`2N`$ modes and every dirty input. The imported reflection theorem gives

```math
\|\widetilde U_VJ-J(U_V\otimes I_b)\|\le\delta.
```

It is not merely an upper-anchor-column promise. Multiplying this
comparison by actual and ideal inverses gives the same initialized bound
for $`\widetilde U_V^\dagger`$. Ordinary unitary telescoping therefore
controls the literal word

```math
\widetilde{\mathcal S}
=\widetilde U_V Z\widetilde U_V^\dagger Z,
\qquad Z=\mathrm{diag}(I_N,-I_N),
```

on the entire initialized domain, with error at most $`2\delta`$.
Actual intermediate leakage is propagated by actual unitaries; neither
clean nor dirty work is reinitialized. Tensoring an arbitrary reference
does not change any of these operator-norm statements.

At the ideal level,

```math
\mathcal S=(2VV^\dagger-I_{2N})Z.
```

On each upper-carrier/lower-children support, this is exactly $`S'_v`$.
On the unused lower root, both factors contribute $`-1`$, so their
product is $`+1`$. Thus the entire direct sum of canonical local stars
is implemented, with no discarded global or support-dependent phase.

### Precision and width

For real root columns, the source has four fine diagonal occurrences;
the oracle and its actual source inverse have eight. The one-clean
literal diagonal bound gives total error below $`2^{-N}`$ with

```math
q=N+\lceil\log_2(30\cdot8)\rceil+1=N+9,
\qquad b_{\rm core}=q+1=N+10\le N+n+7\quad(n\ge3).
```

The eight fine-word errors total less than $`2^{-N}/2`$; reserve, for
example, another $`2^{-N}/4`$ for certified classical coefficient
rounding. The bounded midpoint denominator makes this approximation
effective without dividing by a small branch.

All fine calls use the same dirty core. Their addresses have $`n+1`$
bits, their T counts are $`O(2N+q)=O(N)`$, and self-borrowing requires
$`q+1\ge2(n+1)`$, satisfied here. The exact coarse and permutation
costs are also $`O(N)`$ and reuse the same dirty pool. Any additional
finite-small-dimension threshold of those imported interfaces retains
their direct-synthesis fallback; the displayed inequality alone is not
a separate small-dimension proof.

For complex child components, conjugate the real-magnitude reflection
by its diagonal phase table. This adds two fine occurrences per source,
so the oracle has at most twelve and can take $`q=N+10`$. The core fits
for $`n\ge4`$; smaller dimensions use the same finite fallback.
Certified rational state rounding avoids an equality test at zero phases.
The real-target width claim does not depend on this extension.

## 6. What remains to be recovered

The priced oracle does not directly act as $`W'`$ on the external input
space. The Cayley and triangular formulas still require a native
realization of the resolvent or a separately proved feedback recovery with aggregate
$`O(N)`$ T count. A precision-dependent sequence of $`O(N)`$-T oracle
calls does not establish that bound. The unchanged-angle query obstruction
in the [scattering account](ENDPOINT_TREE_TRANSPORT.md#a-fixed-number-of-unchanged-coin-calls-does-not-give-the-frame)
concerns its specified oracle and parameter-independent interlayers; it
does not automatically apply to this classically reparameterized oracle.

The [operator tests](../../tests/test_canonical_star_resolvent.py) compare
the resolvent with independently composed complete stars, verify the
same-root input correction, and check the full local oracle, its zero
padding, and the support/width counts. These are finite checks of signs,
orientation, and complete-operator identities, not native synthesis of
the feedback inverse.
