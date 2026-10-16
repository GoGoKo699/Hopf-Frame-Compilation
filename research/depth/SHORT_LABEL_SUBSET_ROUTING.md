# Short-label subset routing

This note proves two finite classical routing lemmas and their reversible
query words. They are supporting primitives, not a new numbered
one-clean depth theorem. Both words use a short initialized metadata
register, return the native query work, and have literal inverses.
Neither word supplies a frame-label eraser by itself.

The first lemma chooses a bijection between two dense sets. It does not
implement a prescribed arbitrary permutation. The second lemma compresses
a specified sparse set into a short code chosen by the construction. It
retains its stage label as part of that code; it does not produce a
prescribed rank ordering.

## 1. A short shift catalogue between dense sets

Let $`V=2^m`$, and let $`S,T\subseteq\mathbb F_2^m`$ have the same
positive cardinality $`a=\alpha V`$. There is a multiset of

```math
d=O((m+1)/\alpha^2)
```

shifts $`v_1,\ldots,v_d`$ such that the bipartite graph with left and right
vertex sets $`\mathbb F_2^m`$ and edges
$`x\longmapsto x\oplus v_i`$ has a perfect matching from $`S`$ to $`T`$.
Repeated shifts and parallel edges are allowed. A matching still gives a
bijection of vertices.

### 1.1 Existence of a suitable shift multiset

Choose the shifts independently and uniformly. Let $`A`$ be the
adjacency matrix of the resulting bipartite multigraph. Its constant
singular direction has value $`d`$. Every other Fourier eigenvalue is

```math
\lambda_z=\sum_{i=1}^d(-1)^{z\cdot v_i},\qquad z\ne0.
```

A scalar Hoeffding bound and a union bound over fewer than $`V`$
nonzero characters show that, for the stated choice of $`d`$ with a
sufficiently large absolute constant, there is positive probability that

```math
\lambda:=\max_{z\ne0}|\lambda_z|\le\alpha d/8.
```

For each vertex $`x`$ and either set $`U\in\{S,T\}`$, its number of
neighbors in $`U`$, counted with multiplicity, is a sum of $`d`$
independent Bernoulli variables with mean $`\alpha d`$. A Chernoff bound
and a union bound over these $`2V`$ choices give, simultaneously,

```math
\deg_U(x)\ge7\alpha d/8.
```

Both events hold together for $`d=C(m+1)/\alpha^2`$ with an absolute
constant $`C`$. This argument is an existence proof with finite
preprocessing; efficient construction of the catalogue is not assumed.

The Fourier bound gives the usual mixing inequality, directly by
splitting indicator vectors into their constant and orthogonal parts:

```math
e(A_0,B_0)\le d|A_0||B_0|/V+\lambda\sqrt{|A_0||B_0|}.
```

Here the repeated symbol $`A_0`$ denotes a vertex subset, not the
adjacency matrix, and edges are counted with multiplicity.

### 1.2 Hall's condition, including the large-set case

Suppose $`A_0\subseteq S`$ has $`k\le\lceil a/2\rceil`$ elements and
its neighborhood $`B_0\subseteq T`$ has fewer than $`k`$ elements.
Every edge from $`A_0`$ into $`T`$ ends in $`B_0`$, so the degree lower
bound gives

```math
e(A_0,B_0)\ge7\alpha d k/8.
```

On the other hand, $`|B_0|\le k-1\le a/2`$, and mixing gives

```math
e(A_0,B_0)\le d k(k-1)/V+\lambda k\le5\alpha d k/8,
```

a contradiction. The same argument applies with $`S,T`$ exchanged.

If a Hall violation has $`k\gt\lceil a/2\rceil`$, choose
$`a-k+1`$ vertices of $`T\setminus N(A_0)`$. This number is at most
$`\lceil a/2\rceil`$, and their neighborhood in $`S`$ lies in
$`S\setminus A_0`$, which has only $`a-k`$ elements. This is a small
Hall violation in the reverse direction. Thus no violation exists.

## 2. Three queries implement the chosen matching

Fix a matching $`f:S\to T`$, and choose for every $`x\in S`$ an index
$`c(x)\in\{1,\ldots,d\}`$ satisfying

```math
f(x)=x\oplus v_{c(x)}.
```

Set $`v_0=0`$. Extend $`c`$ by zero outside $`S`$, and define
$`b(y)=c(f^{-1}(y))`$ on $`T`$, extended by zero elsewhere.
An initialized $`s=\lceil\log_2(d+1)\rceil`$-bit register $`C`$
suffices for the word

```math
\begin{aligned}
C&\mathrel{\oplus}=c(X),\\
X&\mathrel{\oplus}=v_C,\\
C&\mathrel{\oplus}=b(X).
\end{aligned}
```

On the reached subspace it acts exactly as

```math
|x\rangle_X|0\rangle_C\longmapsto
|f(x)\rangle_X|0\rangle_C,\qquad x\in S.
```

The first and third queries read $`X`$ and target $`C`$; the middle
query reads $`C`$ and targets $`X`$. Each read/target pair is disjoint.
Every query is a unitary XOR oracle on all inputs, so their composition
has a full-input unitary extension and an actual inverse. No clean-return
claim is made outside the reached $`S\otimes|0\rangle_C`$ subspace.
In particular, the arbitrary extension need not be the identity there.

For $`\alpha^{-1}=\mathrm{poly}(m)`$, the metadata has $`O(\log(m+2))`$
bits. The table address widths are $`m,s,m`$, and their total volume is
$`2^m\mathrm{poly}(m)`$. The native read-only XOR-query compiler may
therefore be applied with only a constant number of query appearances.
Its private dirty work is returned by each completed query.

The same proof applies to finitely many disjoint source/target classes,
with a separate catalogue per class, provided the concatenated catalogue
is polynomial in the outer parameters. The source and output tables
select the class as well as its shift. This qualification matters: an
exponentially sparse class does not satisfy the polynomial-catalogue
bound merely because the number of classes is small.

## 3. A fixed-width staged shear compressor

Let $`S\subseteq\mathbb F_2^m`$ have $`r\ge1`$ elements. Define

```math
k=\min\{m,\lceil\log_2(4r)\rceil\}.
```

If $`k=m`$, the identity already uses at most
$`\log_2 r+3`$ data bits, and no compression is needed. Suppose henceforth
that $`k\lt m`$, and write the data as $`X=(H,L)`$, with $`L`$ containing
exactly $`k`$ bits. The split is fixed throughout this lemma.

For a binary matrix $`M`$ with $`k`$ rows and $`m-k`$ columns, use the
linear shear hash

```math
h_M(H,L)=L\oplus MH.
```

For two distinct inputs with different high parts, a uniformly random
$`M`$ makes their hashes collide with probability exactly $`2^{-k}`$.
If their high parts agree, their low parts differ and collision is
impossible.

### 3.1 A logarithmic-size stage label

Start with the remaining set $`S_1=S`$. At stage $`j`$, let
$`R_j=|S_j|`$ and choose a shear matrix $`M_j`$ for which at least half
the points of $`S_j`$ have a hash value unique within $`S_j`$.
Such a matrix exists: the expected number of points having at least one
collision is at most

```math
R_j(R_j-1)2^{-k}\le R_j/4.
```

Assign all isolated points to this stage and remove them. After at most
$`t\le\lceil\log_2 r\rceil+1\le m+1`$ stages every point is assigned.
Let $`j(x)\in\{1,\ldots,t\}`$ be its stage, and set $`j(x)=0`$ outside
$`S`$. Within each assigned stage set, $`h_{M_j}`$ is injective.

For stage $`j\ge1`$ and a hash value $`z`$ attained by an assigned point,
let $`g(j,z)`$ be that unique point's high part. Extend $`g`$ by zero on
unused pairs and on stage zero, and set $`M_0=0`$.

### 3.2 The reversible word and its exact domain

Initialize a stage register $`J`$ of
$`\lceil\log_2(t+1)\rceil=O(\log(m+2))`$ bits. Apply

```math
\begin{aligned}
J&\mathrel{\oplus}=j(H,L),\\
L&\mathrel{\oplus}=M_JH,\\
H&\mathrel{\oplus}=g(J,L).
\end{aligned}
```

For $`x=(H_x,L_x)\in S`$, the output is

```math
|0^{m-k}\rangle_H
|h_{M_{j(x)}}(x)\rangle_L
|j(x)\rangle_J.
```

The pair $`(j(x),h_{M_{j(x)}}(x))`$ is an injective code. Its total
length is at most $`\log_2 r+O(\log(m+2))`$. The stage register is part
of the output code, not a returned helper. Outside $`S`$, an initialized
stage register remains zero and the word is exactly the identity.

Again every line is a full-input XOR unitary. The actual reverse word
recovers $`x`$ and resets $`J`$ on the code image. On arbitrary inputs
the same reverse word is still the inverse unitary; no interpretation
of arbitrary nonzero stage inputs is required.

The read/target separation is literal:

| Query | Read address | XOR target |
|---|---|---|
| Stage | $`H,L`$ | $`J`$ |
| Shear | $`J,H`$ | $`L`$ |
| High-part erasure | $`J,L`$ | $`H`$ |

The three table volumes are bounded by
$`2^m\mathrm{poly}(m)`$. The shear can be treated as an ordinary finite
XOR table; its efficient classical evaluation is not a free coherent
oracle premise. Finite preprocessing chooses the matrices, assigned
sets, and inverse tables.

## 4. Native scope, selectors, and errors

The words above are compatible with the native read-only query
interfaces in [uniform subtree depth](UNIFORM_SUBTREE_DEPTH.md).
They use a constant number of query appearances. At fixed query
accuracy their T-depth is $`O(\log(m+2))`$ when the corresponding query
table volume and private dirty pools are paid. Their counts are the
query compiler's counts for $`2^m\mathrm{poly}(m)`$ table volume;
this note does not replace that quantity by $`O(2^m)`$.

A preserved external prefix may select the sets, matching, or shear
matrices. Its bits then occur in each relevant table address, and the
entire coherent family must be charged once at its global table volume.
A preserved activity bit may also guard every table, including the
middle shift or shear table: set every inactive row to zero. Each
completed query is then literally the identity in the inactive sector,
even when the short metadata is arbitrarily occupied.
Likewise, a variable split $`k`$ is not automatically covered by the
fixed-split word in Section 3: a layout lemma must first keep every
query's read wires disjoint from its target wires across the selected
family. Serially applying a separate full word for every selector value
is not the claimed constant-query implementation.

If the three physical queries approximate their ideal full XOR
unitaries in operator norm by $`\epsilon_1,\epsilon_2,\epsilon_3`$,
their product differs from the displayed ideal word by at most their
sum. This comparison remains valid for superpositions of the reached
inputs and with arbitrary external references. Using the actual
physical inverse costs the same norm error. Native query work is
returned according to the full query contract; the short initialized
metadata requirements stated above remain explicit.

These are conditional-workspace primitives. They do not establish a
standalone one-clean circuit, an arbitrary prescribed rank permutation,
or a complete frame recovery. A prospective compact dilation must still
prove its source normalization, reference transform, selected layout,
and all output/metadata return conditions.
