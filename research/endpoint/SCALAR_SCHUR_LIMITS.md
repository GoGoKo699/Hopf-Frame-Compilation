# Tree-correlation gauges and source-dimension limits

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Research index](../README.md)

A complete frame has an exact coarse-frame Schur representation with
bounded outer gauges, including zero target branches. Its correlation
factor has a support-forced rank obstruction. The scalar proof extends to
bounded operator-valued sources acting on arbitrarily large dirty work,
with an error floor independent of that dirty dimension. These results
concern logical-index-preserving sources around one prescribed coarse call;
they are not unrestricted circuit or T-count lower bounds.

## 1. Exact gauge with a complex native coarse frame

Take target splits with real $`c_v\ge1/\sqrt2`$, $`s_v\ge0`$, and
$`c_v^2+s_v^2=1`$. The signed-frame reduction and input correction are
in the [precision-gauge compiler](PRECISION_GAUGE_COMPILER.md).
Let each actual coarse SU(2) gate be

```math
G_v=\begin{pmatrix}a_v&-\overline b_v\\b_v&\overline a_v\end{pmatrix},
\qquad |a_v|^2+|b_v|^2=1,\qquad |a_v|,|b_v|>0.
```

Assume $`|b_v|\ge s_v`$ and define

```math
\lambda_v=\frac{c_v}{|a_v|},\qquad
q_v=\frac{s_v|a_v|}{c_v|b_v|}\in[0,1].
```

On binary-tree nodes set $`d_{\rm root}=1`$ and

```math
d_{v0}=d_v\lambda_v\frac{|a_v|}{a_v},\qquad
d_{v1}=d_v\lambda_v\frac{|b_v|}{b_v}.
```

Let $`D`$ be the leaf diagonal, and $`j_v`$ the right-child
representative/detail label of $`v`$. Define

```math
C_0=1,\qquad C_{j_v}=\frac{a_vb_v}{|a_v||b_v|d_v}.
```

Connect every computational label $`j>0`$ to
$`p(j)=j-\mathrm{lowbit}(j)`$. Put the corresponding $`q_v`$ on
that edge; let $`K_{ij}`$ be the product along the unique path from
$`i`$ to $`j`$, with $`K_{ii}=1`$. For the complete coarse word $`W_0`$,

**Schur identity (1).**

```math
W_\theta=D(W_0\circ K)C,
```

where $`\circ`$ is entrywise multiplication. For column zero, each right
branch contributes $`q_v`$, giving exactly the target branch products.
For a detail column, the path to the left half crosses its own edge, while
the path to the right half does not. Its two initial gauged coefficients
are respectively $`-s_v,c_v`$; descendants obey the column-zero calculation.
Outside its support both sides vanish. This proves every column of (1)
and retains all complex coarse phases.

An inverse-polynomial upward angular margin before native rounding can
ensure $`\lambda_v\le\exp(O(\delta))`$ for
$`\delta=1/\mathrm{poly}(n)`$, with both coarse branches inverse-polynomially
bounded. Then $`1\le|D_i|\le\exp(O(n\delta))`$ and $`|C_j|\le1`$.
No tiny target coefficient is inverted. Native synthesis uses the
[coarse encoder's interfaces](COARSE_PREFIX_ENCODER.md), not an assumed
exact native implementation of an arbitrary real rotation.

For orthonormal formal vectors $`e_j`$, the recursion

```math
v_0=e_0,\qquad v_j=q_jv_{p(j)}+\sqrt{1-q_j^2}\,e_j
```

gives $`\langle v_i,v_j\rangle=K_{ij}`$. This abstract Gram source does
not price a native circuit. Its aggregate support is at most
$`\sum_j(1+\mathrm{popcount}(j))=N(1+n/2)`$, with equality when all
edge correlations lie strictly between zero and one.

## 2. Every complex completion has a forced determinant

Consider the representative chain

```math
x_j=N-2^{n-j},\qquad 0\le j\le n.
```

With all coarse branches nonzero, detail column $`x_j`$ has exact support
$`[x_{j-1},x_{j-1}+2^{n-j+1})`$, and column zero has full support.
The chain submatrix therefore supports the whole lower triangle and first
superdiagonal. Farther upper entries need not be supported.

For $`m\le n+1`$, let $`G`$ be any complex $`m\times m`$ matrix
agreeing with $`K`$ on that pattern:

```math
G_{ij}=\prod_{\ell=j+1}^{i}q_\ell\quad(i\ge j),\qquad
G_{j-1,j}=q_j.
```

The unspecified upper entries are arbitrary. Without any positivity or
Hermitian assumption,

**Completion determinant (2).**

```math
\det G=\prod_{j=1}^{m-1}(1-q_j^2).
```

Indeed, right multiplication by the lower bidiagonal matrix with unit
diagonal and $`-q_{j+1}`$ in position $`(j+1,j)`$ simultaneously
replaces column $`j`$ by its original value minus $`q_{j+1}`$ times
the original column $`j+1`$. The result is upper triangular with diagonal
$`1-q_1^2,\ldots,1-q_{m-1}^2,1`$; the multiplier has determinant one.

If every selected $`q_j\lt1`$, each completion has rank at least $`m`$.
Thus a scalar two-source factorization $`G_{ij}=\langle u_i,v_j\rangle`$
needs dimension at least $`n+1`$. Different source families and
non-Hermitian completions do not evade this exact bound.

## 3. Robust completion bound

Suppose $`\mathrm{rank}(G)\lt m`$, $`|G_{ij}|\le B`$, and the maximum
entry error on the required pattern is $`\delta`$. Set
$`M=\max(B,1)`$, $`\Delta=\prod_j(1-q_j^2)`$, and replace just the
required entries by their exact values to obtain $`\widetilde G`$.
Its determinant is $`\Delta`$, whereas $`\det G=0`$. Telescoping
column by column and applying Hadamard's inequality gives

**Robust bound (3).**

```math
\Delta\le m^{1+m/2}\delta M^{m-1},\qquad
\delta\ge\frac{\Delta}{m^{1+m/2}M^{m-1}}.
```

Each difference column has norm at most $`\sqrt m\delta`$; every other
old or new column has norm at most $`\sqrt m M`$. The bounded-amplitude
hypothesis is essential for the approximate statement.

## 4. Polynomial operator-error floor at fixed scalar dimension

Fix $`r\ge1`$ and $`n\ge r`$. First take uniform real coarse angle
$`\beta=1/(10n)`$ and target $`\theta=\arctan(\tan\beta/2)`$.
All correlations are $`1/2`$. Write $`\lambda=\cos\theta/\cos\beta`$.
Then $`D=\lambda^n I`$ and $`C_j=\lambda^{-d}`$ for a detail label
at binary-node depth $`d`$. Put $`H=DW_\beta C`$. On the first
$`r+1`$ chain labels the required entries are

```math
H_{x_i,0}=\cos^n\theta\,\tan^i\beta,
```

```math
H_{x_i,x_j}=\cos^{n-j+1}\theta\,\tan^{i-j}\beta
\quad(i\ge j\ge1),\qquad
H_{x_j,x_{j+1}}=-\cos^{n-j}\theta\,\tan\beta.
```

Since $`\cos^n\theta\ge1/2`$, each magnitude is at least
$`\alpha=1/[2(10n)^r]`$. If $`\widehat W=D(W_\beta\circ G)C`$
with $`\mathrm{rank}(G)\le r`$ and $`|G_{ij}|\le B`$, an operator
error $`\eta`$ bounds the required entry errors of $`G`$ by
$`\eta/\alpha`$. Applying (3) with $`m=r+1`$ yields

**Operator gap (4).**

```math
\eta\ge
\frac{(3/4)^r}
{2(10n)^r(r+1)^{1+(r+1)/2}\max(B,1)^r}.
```

For a native-compatible instance choose an actual exact native SU(2)
gate approximating that coarse rotation closely enough that
$`1/(20n)\le|b|\le1/(5n)`$. Set
$`\tan\theta=|b|/(2|a|)`$ and repeat the actual gate at every node.
Equation (1) retains its phases; the same magnitude proof uses
$`\beta_{\rm eff}=\arcsin|b|`$. Thus (4) holds with $`20n`$ replacing
$`10n`$. Only inverse-polynomial native coarse accuracy is needed;
no exactly real native approximation is presumed.

For fixed $`r,B`$ the gap is $`\Omega_{r,B}(n^{-r})`$, exceeding
$`2^{-N}`$ for sufficiently large $`n`$. In particular a four-dimensional
scalar environment has a polynomial floor in this interface. Fixed block
normalization only changes the bounded-amplitude constant.

## 5. Operator-valued sources with arbitrary dirty dimension

Let $`\mathcal H_D`$ be the dirty input space of finite dimension $`d`$.
For every logical output $`i`$ and input $`j`$, take source columns

```math
A_i,B_j:\mathcal H_D\longrightarrow\mathbb C^r\otimes\mathcal H_D,
\qquad G_{ij}=A_i^\dagger B_j,\qquad \|G_{ij}\|\le B.
```

Isometric source columns give $`B=1`$. The two source operations preserve
logical labels; the sole intervening coarse call acts as $`W_0`$ on those
labels and as identity on the source modes and dirty space. With the fixed
outer gauges from (1), the accepted block has entries
$`H_{ij}G_{ij}`$, where $`H=DW_0C`$. Allowing these scalar gauges for free
only enlarges this restricted architecture. Any different gauge choice
must re-establish the bounded-coefficient hypothesis.

On $`m`$ selected labels, the collective source maps have domain
$`\mathcal H_D^{\oplus m}`$ and codomain
$`\mathbb C^r\otimes\mathcal H_D`$. Consequently

```math
G=A^\dagger B,\qquad \mathrm{rank}G\le rd.
```

In particular $`G`$ is singular when $`m>r`$, even though its blocks need
not commute. Suppose its required chain entries approximate
$`K_{ij}I_D`$ to operator error at most $`\delta`$, with
$`0\le q_j\lt 1`$. Replace just those entries by their exact values to obtain
$`\widetilde G`$. Set

```math
M=\max(B,1),\qquad Q=1+\max_j q_j,\qquad
\kappa=\min_j(1-q_j^2)>0.
```

Every block of $`\widetilde G`$ has norm at most $`M`$, and the block norm
estimate gives $`\|G-\widetilde G\|\le m\delta`$. Right multiply by the
same scalar lower bidiagonal matrix $`T`$ used in Section 2, tensored with
$`I_D`$. Its norm is at most $`Q`$. The block matrix
$`F=\widetilde G T`$ is upper triangular: its diagonal blocks are
$`(1-q_1^2)I_D,\ldots,(1-q_{m-1}^2)I_D,I_D`$, and every strict upper
block has norm at most $`QM`$. Only scalar multiplication is used in the
cancellations, so the argument does not assume commuting dirty operators.

For $`Fx=y`$, put $`Y=\max_i\|y_i\|`$. Back substitution gives

```math
\|x_i\|\le\frac{Y}{\kappa}
  +\frac{QM}{\kappa}\sum_{j>i}\|x_j\|
\le\frac{Y}{\kappa}
  \left(1+\frac{QM}{\kappa}\right)^{m-1-i}.
```

Thus

```math
\|\widetilde G^{-1}\|=\|TF^{-1}\|
\le\frac{Q\sqrt m}{\kappa}
  \left(1+\frac{QM}{\kappa}\right)^{m-1}.
```

Because $`G`$ is singular, its distance from the invertible matrix
$`\widetilde G`$ is at least $`1/\|\widetilde G^{-1}\|`$. Combining the
two estimates proves the following bound, independent of $`d`$.

**Operator-valued completion bound (5).** For $`m>r`$,

```math
\delta\ge
\frac{\kappa}
 {Qm^{3/2}\left(1+QM/\kappa\right)^{m-1}}.
```

At zero error, the block triangular calculation itself forces $`r\ge m`$.
For $`q_j=1/2`$ and isometric source columns, (5) becomes
$`\delta\ge1/[2m^{3/2}3^{m-1}]`$.

For $`n\ge r\ge1`$, apply this to the native-compatible witness of Section 4
with $`m=r+1`$. Its fixed gauged coarse entries on the required pattern have
magnitude at least $`\alpha=1/[2(20n)^r]`$. Here the effective angle is
$`\beta_{\rm eff}=\arcsin|b|`$, and the lower bound uses
$`\tan\beta_{\rm eff}\ge|b|\ge1/(20n)`$. If the accepted block
approximates $`W_\theta\otimes I_D`$ with operator error $`\eta`$, each
required block obeys

```math
|H_{ij}|\,\|G_{ij}-K_{ij}I_D\|\le\eta.
```

Hence $`\delta\le\eta/\alpha`$, and (5) yields

**Dirty-dimension-independent operator gap (6).** For isometric sources,

```math
\eta\ge
\frac{1}{4(20n)^r(r+1)^{3/2}3^r}.
```

For fixed $`r`$, arbitrarily many dirty qubits therefore cannot remove the
polynomial error floor of this interface. A fixed normalization changes
$`B`$ in (5), retaining a polynomial floor for fixed $`r,B`$. All norms
are unchanged upon tensoring an arbitrary reference identity. Requiring
the target blocks to be $`K_{ij}I_D`$ enforces the entire dirty input
space, not merely one dirty state or an averaged matrix element. Since
projection onto the accepted flags cannot increase an error norm, the
same lower bound applies to a full initialized-input isometry satisfying
this architecture.

## 6. Scope and evidence

This result covers one prescribed coarse call between two
logical-index-preserving sources with bounded scalar or dirty-operator
coefficients. It does not cover source operations mixing logical labels,
additional retained labels beyond the stated source space, growing source
dimension, multiple coherent coarse queries, coarse transport acting
nontrivially on the source or dirty space, or every adaptive choice of a
different coarse frame. The [canonical-star resolvent](CANONICAL_STAR_RESOLVENT.md)
uses a distinct local oracle; its descendant inverse remains unpaid.

[Finite tests](../../tests/test_scalar_schur_limits.py) check the full
complex-phase gauge, support, exact rational-complex determinants with
free upper entries, and minor-weight bounds. Those tests concern the scalar
identities; the operator-valued extension follows from the block rank and
inverse estimates above. The rank and asymptotic operator-error conclusions
follow from the analytic proofs.
