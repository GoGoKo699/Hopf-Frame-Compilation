# Compact fixed-coin representations of coarse frame flags

[R87 product coins](COIN_DILATION_DEPTH.md) · [Native family sources](COHERENT_DESCENDANT_LOADER.md) · [Deterministic-frame scope](DETERMINISTIC_FRAME_DEPTH.md)

This note gives an exact finite representation with smaller, nonuniform
coins. It makes no independent native compiler claim. The separate
[causal compact-coin theorem R89](CAUSAL_COMPACT_COIN_DEPTH.md) supplies
an in-place embedding, fixed reference transforms, and their native
resource ledger. A small ambient coin register alone does not supply
the initialized target of a label-preserving source; Section 6 explains
how the causal implementation avoids that particular requirement.

## 1. A fixed atom distribution for all finite dyadic probability vectors

Let M and K be positive integers. There is a fixed list of positive
integer atom weights $`a_1,\ldots,a_D`$, independent of any target
probability vector, with the following properties:

```math
\sum_{z=1}^D a_z=K,\qquad
D=O\bigl(M(1+\log(K+1))\bigr).
```

Every vector of nonnegative integers $`(k_1,\ldots,k_M)`$ summing
to K is the sum of a partition of this same list. All atom weights may
be chosen to be powers of two.

Here is a terminating construction. Start with remaining total T=K.
When $`T\ge M`$, append the largest power of two at most T/M;
when $`0\lt T\lt M`$, append one. Subtract that weight and continue.
At every step the appended weight is positive and at most T.
For $`T\ge M`$ it exceeds $`T/(2M)`$, so T decreases by a
factor less than $`1-1/(2M)`$. There are
$`O(M(1+\log(K+1)))`$ such steps and fewer than M final unit steps.

To assign these fixed atoms to a given target vector, maintain its
remaining integer capacities. Their sum is T. A largest remaining
capacity is at least T/M and therefore at least the next atom weight.
Assign the atom to such a capacity, breaking ties by the lowest label,
and subtract it. In the final unit regime some capacity is positive.
All capacities remain nonnegative and finish at zero. This proves the
partition assertion, including zero target probabilities.

Let $`w_z=a_z/K`$ and pad the atom list with zero weights to a
power-of-two dimension. The fixed real coin state is

```math
|\gamma\rangle=\sum_z\sqrt{w_z}|z\rangle.
```

For a positive target vector $`p_x=k_x/K`$, write $`\chi(z)=x`$
for its atom assignment and put

```math
E_p|x\rangle=
\frac1{\sqrt{p_x}}\sum_{z:\chi(z)=x}\sqrt{w_z}|z\rangle.
\qquad(1)
```

The columns have disjoint support and unit norm. Moreover,

```math
E_p^\dagger E_p=I,\qquad
E_p\sum_x\sqrt{p_x}|x\rangle=|\gamma\rangle.
\qquad(2)
```

Zero-probability outcomes are allowed in the partition lemma, but (1)
is asserted only for positive outcomes. R87's globally rounded tree has
strictly positive branches, so every outcome needed below is positive.
The membership function on padded zero-weight atoms can be assigned
arbitrarily; those atoms do not occur in these isometry images.

For M=2 and $`K=2^s`$, the construction gives the familiar finite
weight list

```math
\frac12,\frac14,\ldots,2^{-s},2^{-s}.
```

Thus one binary branch of s-bit dyadic precision requires only
$`\lceil\log_2(s+1)\rceil`$ coin qubits in this representation.
Its reference state is nonuniform; it is not a smaller uniform coin.

## 2. A second finite construction with two outcomes per column

The following construction makes membership more structured. Assume
M outcomes with capacities $`k_x`$ summing to K. Set integer residuals
$`c_x=Mk_x`$. Maintain small labels with $`c_x\lt K`$ and large
labels with $`c_x\ge K`$. If both kinds occur, choose a small
label j and a large label l, record

```math
\alpha_j=c_j/K,\qquad \mathrm{parent}(j)=l,
```

remove j, and replace $`c_l`$ by $`c_l-(K-c_j)`$. Reclassify l.
When no small label remains, every remaining residual equals K:
the residual sum equals K times the number of surviving labels.
Set $`\alpha_j=1`$ at each survivor.

Each removed label points to one that is still present. The resulting
directed graph has no cycle. For a uniformly chosen column j, output j
with probability $`\alpha_j`$ and its recorded parent otherwise.
This outputs x with probability exactly $`k_x/K`$. One proof is to
view the invariant residual $`c_x`$ as its original scaled capacity
minus all mass already assigned to x by removed columns. Processing j
assigns exactly its remaining $`c_j`$ units to j and exactly
$`K-c_j`$ units to l. The last survivor columns complete the capacities.

For dyadic K, each acceptance probability has denominator K. Use the
fixed binary coin of Section 1 and its two-outcome partition to realize
that acceptance test. The product of a uniform M-valued column and this
fixed coin has $`M(1+\log_2 K)`$ positive atoms. Its membership is
only

```math
\chi(j,z)=
\begin{cases}
j,&z\text{ is accepted at }j,\\
\mathrm{parent}(j),&z\text{ is rejected at }j.
\end{cases}
\qquad(3)
```

This is an alternative fixed atom distribution with a finite acyclic
assignment rule. It is not an in-place permutation: many removed
columns can point to the same parent. Reversing (3) coherently requires
preparing the normalized state on all incoming column/coin pairs.
An acceptance bit alone does not recover the previous column.

## 3. Products of compact coins for consecutive tree blocks

Use R87's one globally rounded positive tree, with each branch
probability an integer multiple of $`2^{-\tau}`$. Divide its q
levels into consecutive blocks of heights
$`q_1,\ldots,q_p`$, and write
$`d_0=0,d_i=q_1+\cdots+q_i`$.

At a prefix u of length $`d_{i-1}`$, the conditional distribution
on the $`M_i=2^{q_i}`$ coarse children has denominator
$`K_i=2^{q_i\tau}`$. Apply either fixed-atom construction to
this vector. Its atom weights $`w_{i,z}`$ and coin
$`\gamma_i`$ depend only on $`M_i,K_i`$, whereas its membership
$`\chi_{i,u}`$ depends on the preserved prefix and the finite tree.
Its number of qubits is

```math
r_i=\lceil\log_2 D_i\rceil
\le q_i+O\bigl(\log(q_i\tau+1)\bigr).
\qquad(4)
```

A tuple of coin values is decoded block by block: apply
$`\chi_{i,u}`$ using the already decoded preceding prefix u.
For a logical leaf x let $`\mathcal R_x`$ be the set of coin
tuples decoding to x. These sets partition the positive-weight coin
tuples. Put

```math
\omega(z)=\prod_i w_{i,z_i},\qquad
\mu_x=\sum_{z\in\mathcal R_x}\omega(z),\qquad
E|x\rangle=\frac1{\sqrt{\mu_x}}
\sum_{z\in\mathcal R_x}\sqrt{\omega(z)}|z\rangle.
\qquad(5)
```

The conditional partitions imply that $`\mu_x`$ is exactly the
original leaf probability. Thus E is a positive isometry with disjoint
columns. No typical-set truncation or root-distribution-only estimate
is used.

At a coarse boundary prefix u of length $`d_i`$, its normalized
scaling vector obeys

```math
E\rho_u=|\Gamma_u\rangle\otimes
|\gamma_{i+1}\rangle\otimes\cdots\otimes|\gamma_p\rangle,
\qquad(6)
```

where $`\Gamma_u`$ is the normalized product-coin state on the
prefix rectangle. Indeed each descendant coefficient
$`\sqrt{\mu_x/\mu_u}`$ cancels the $`\sqrt{\mu_x}`$
denominator in (5). The remaining amplitudes are
$`\sqrt{\omega(z)/\mu_u}`$ on that rectangle, with unrestricted
future coins.

Now take an original marker at depth
$`d_{i-1}\lt t\le d_i`$. It has a computational prefix u of length
$`d_{i-1}`$. Within that coarse subtree, expand it in the normalized
coarse-child scaling vectors at depth $`d_i`$; call the resulting
unit vector $`v_{u,t}`$. The marker is orthogonal to that subtree's
positive root, so

```math
\langle\rho_{u,\mathrm{block}},v_{u,t}\rangle=0.
```

Equations (1)–(2) give a unit vector
$`\zeta_{u,t}=E_{i,u}v_{u,t}`$ orthogonal to
$`\gamma_i`$. Therefore

```math
Ew=|\Gamma_u\rangle\otimes|\zeta_{u,t}\rangle\otimes
|\gamma_{i+1}\rangle\otimes\cdots\otimes|\gamma_p\rangle,
\qquad \langle\gamma_i,\zeta_{u,t}\rangle=0.
\qquad(7)
```

The expansion retains the marker's literal sign convention. By
linearity the identity holds on the entire coarse band
$`\mathrm{Ran}(A_{d_i}-A_{d_{i-1}})`$, including arbitrary
superpositions of its prefixes and original depths.

## 4. An exact coarse flag for any fixed coin completions

For each i choose any unitary $`U_i`$ whose zero column is
$`\gamma_i`$, and let $`F_i=U_i^\dagger`$.
Such finite algebraic completions exist; a Householder construction
suffices. This existence is not a native implementation claim.

Write $`F=\bigotimes_i F_i`$ and let
$`\lambda(z)=\max(\{i:z_i\ne0\}\cup\{0\})`$.
Since $`F_i\gamma_i=|0\rangle`$ and
$`\langle0|F_i\zeta_{u,t}\rangle=0`$, (6)–(7) prove

```math
FE\rho_\varnothing=|0\rangle,
\qquad
FE\mathrm{Ran}(A_{d_i}-A_{d_{i-1}})
\subseteq\mathrm{span}\{|z\rangle:\lambda(z)=i\}.
\qquad(8)
```

This is sector containment, not surjectivity onto the larger coin
sector. It is independent of the chosen action of $`U_i`$ on all
other columns, including padded zero-weight coordinates.

Consequently, with an initialized ideal embedding available, the word
E, F, XOR $`\lambda(C)`$ into a coarse tag, $`F^\dagger`$,
$`E^\dagger`$ erases a consistent old coarse tag and restores the
logical vector with phase +1. The inverse embedding is only applied on
its reached image. This statement is the same exact image identity as
R87's coin lemma, for nonuniform block coins. It does not supply the
physical embedding, its inverse, or the fixed coin transforms for free.

## 5. Why a compact coin does not give a short conditional-rank source

There is a finite obstruction to one tempting implementation of (5).
Suppose a fixed atom distribution realizes every M-outcome probability
vector with denominator K, where K is divisible by M. The uniform target
forces every atom weight to be at most $`1/M`$.

For $`K\ge2M`$, consider the positive target

```math
\left(1-\frac{M-1}{K},\frac1K,\ldots,\frac1K\right).
```

Its heavy outcome has mass greater than one half. Its fiber therefore
contains more than $`M/2`$ positive atoms. If a circuit preserves the
old $`\log_2 M`$-bit outcome label, prepares a d-bit conditional
rank, and then applies a basis permutation, that column can have at
most $`2^d`$ nonzero entries. Thus

```math
2^d>M/2.
\qquad(9)
```

For power-of-two M this already requires $`d\ge\log_2 M`$.
An $`O(\log(q\tau+1))`$-bit conditional rank does not suffice
uniformly for $`M=2^q`$. This is a support-count restriction on
that label-preserving construction. It does not exclude an in-place
unitary that changes the old label, different embeddings, or additional
coherent interference.

## 6. Native scope and bounded evidence

The direct family-source realization of (5) preserves the old logical
leaf while preparing the whole coin register, and then uses membership
to erase the leaf. Its actual inverse reconstructs the leaf before
reversing that source. The
[charged source interface](COHERENT_DESCENDANT_LOADER.md) requires
the separate initialized target bank; it does not infer that bank from
(4). Over p blocks, the ambient output length is

```math
\sum_i r_i=q+O\left(\sum_i\log(q_i\tau+1)\right),
```

but preparing this output while retaining a q-bit arbitrary input still
needs a long additional clean register in that realization. R87's
protected suffix does not supply it before a selected band is unloaded.

The acyclic construction in Section 2 makes the remaining map a
weighted incidence isometry, but an outcome may have many incoming
children. A constant-norm weighted comb completion does not by itself
price its graph-label transport or its source table addresses.
Likewise, the deterministic-frame theorem R84 explicitly does not
compile arbitrary tree automorphisms or arbitrary in-place membership
maps. None is invoked here.

Preserving atom weights alone is also insufficient for the full coarse
flag. Consider two one-level blocks, all branch probabilities one half,
and the Section 2 coin with $`\tau=1`$. Each block has a column bit
$`j_i`$ and a fair coin bit $`z_i`$. One product embedding is

```math
E|x_1x_2\rangle
=|x_1\rangle|+\rangle_{z_1}|x_2\rangle|+\rangle_{z_2},
\qquad F=H^{\otimes4}.
```

Swapping $`j_1`$ and $`j_2`$ preserves every atom weight, even
the full tuple of weights of the two blocks. It preserves the global
positive root. But it sends the top marker's transformed nonzero bit
from block one to block two, violating (8). Thus an arbitrary matching
between equal-weight atom classes cannot replace a product embedding.
A proposed compiler must preserve the prefix-rectangle identities, or
provide a separate full-band intertwining argument for its new
reference transform.

[R89](CAUSAL_COMPACT_COIN_DEPTH.md) resolves these implementation
obligations by compressing child labels according to fiber size,
preparing joint conditional ranks, and choosing the weight-class
matchings locally at each original prefix. Short matching selectors
batch those causal local assignments without changing their product
structure. Its separate proof charges fixed coin transforms, source
approximation, occupied inactive work, and all metadata and reference
return. The arbitrary global matching counterexample above remains
valid; local causal assignment is the additional structure used there.

Using compact coins only inside the already selected local episodes
would leave the retained coarse binary search and terminal costs. The
native improvement requires R89's full causal construction, rather
than just the dimension estimate (4).

The [bounded tests](../../tests/test_compact_coin_flags.py) exhaust
small integer capacity vectors, verify the acyclic two-outcome
construction by exact integer mass accounting, and check (6)–(8) on
heterogeneous trees, all marker columns, and noncanonical unitary
completions. They also exhibit the heavy-fiber rank obstruction.
The two-block swap checks the distinction between root preservation and
full coarse-flag preservation.
They are finite representation checks, not native circuit emitters or
asymptotic proofs.
