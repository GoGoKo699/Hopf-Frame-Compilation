# Weighted-chain factorization of finite Haar trees

[R86 interval dilation](INTERVAL_DILATION_DEPTH.md) · [Executed design](INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md) · [Selected routing](SELECTED_DIRTY_ROUTING.md)

This note proves an exact factorization of finite weighted-tree Haar
bases. It is an independent representation result, not a new numbered
native-depth frontier. The simpler coin-dilation construction is the
primary compiler route from the same session. No intermediate Haar basis
below is assumed to be a free physical gate.

## 1. Statement and conventions

A finite full binary tree has D leaves, each of mass one. A node's mass
is its number of descendant leaves. For a node with child masses a and b,
its positive scaling vector is the normalized indicator of its leaves;
its marker is negative on its first child and positive on its second.
The root scaling vector and all internal-node markers form a real
orthonormal basis F of the D-dimensional atom space. Leaves carry fixed
physical atom labels. Intermediate subtrees may be arbitrary nested atom
sets; they need not be intervals in the original atom order.

Suppose every child fraction is at least $`2^{-b}`$, for an integer
$`b\ge2`$. There is a finite sequence of Haar bases
$`F_0=F,F_1,\ldots,F_s`$ with

```math
s=O(\log(b+2)),\qquad
Q_j=F_{j+1}^{\dagger}F_j,\qquad
\|\,|Q_j|\,\|\lt125,
```

such that every child fraction of $`F_s`$ is at least $`2^{-8}`$.
Every Q fixes the root coordinate exactly. Here and below absolute-value
bars around a matrix mean entrywise absolute value; the outer norm is
the Euclidean operator norm.

If $`D=2^r`$, let B be the standard balanced atom-Haar basis. A final
permutation P of physical atom coordinates puts $`F_s`$ in depth-first
leaf order. Then

```math
Q_{\rm end}=B^{\dagger}PF_s,\qquad
\|\,|Q_{\rm end}|\,\|=O(1),
```

and

```math
F=P^{\dagger}B Q_{\rm end}Q_{s-1}\cdots Q_0.
```

All matrices and permutations have terminating certified finite
preprocessing from the integer tree. Their entries are real algebraic
numbers. The construction has no dependence on tree height in its
number of factorization stages.

## 2. A complete weighted comb has constant absolute norm

Let $`\mu_0,\ldots,\mu_m`$ be positive weights summing to one, including
the residual terminal weight $`\mu_m`$. Put
$`M_i=\sum_{j=i}^m\mu_j`$. The comb H has root column
$`v(j)=\sqrt{\mu_j}`$ and marker columns, for $`i\lt m`$,

```math
w_i(j)=
\begin{cases}
0,&j\lt i,\\
-\sqrt{M_{i+1}/M_i},&j=i,\\
\sqrt{\mu_i\mu_j/(M_iM_{i+1})},&j>i.
\end{cases}
```

For $`m=0`$, H is the one-dimensional identity. Otherwise,

```math
\|w_i\|^2=\frac{M_{i+1}}{M_i}+\frac{\mu_i}{M_i}=1,
\qquad
\langle v,w_i\rangle=0.
```

For $`i\lt k`$, the restriction of $`w_i`$ to the tail beginning at k
is proportional to its positive scaling vector. The marker $`w_k`$ is
orthogonal to that vector, so the displayed columns are an orthonormal
basis. The minus sign on the peeled child is the literal marker sign.

Here is a finite proof of the required positive-norm bound. Give
$`\{0,\ldots,m\}`$ mass $`\mu_j`$ at j and define

```math
T_i=\{i+1,\ldots,m\},\qquad
A_i(g)=\frac{1}{M_{i+1}}\sum_{j>i}\mu_jg(j),
\qquad
\mathcal M g(j)=\max\bigl(\{A_i(|g|):i\lt j\}\cup\{0\}\bigr).
```

For $`t>0`$, if k is the smallest index with $`A_k(|g|)>t`$, the
set $`\{\mathcal M g>t\}`$ is exactly $`T_k`$. Thus

```math
t\,\mu\{\mathcal M g>t\}
\le\sum_{j:\mathcal M g(j)>t}\mu_j|g(j)|.
```

If no such k exists, both sides vanish. Integrating over t on this finite
space and using Cauchy--Schwarz gives

```math
\|\mathcal M g\|_{L^2(\mu)}^2
\le2\sum_j\mu_j|g(j)|\mathcal M g(j)
\le2\|g\|_{L^2(\mu)}\|\mathcal M g\|_{L^2(\mu)}.
```

Consequently $`\|\mathcal M g\|_2\le2\|g\|_2`$. If nonnegative
weights satisfy $`\sum_{i=k}^{m-1}\lambda_i\le C M_{k+1}`$, the same
smallest-index argument gives

```math
\sum_{i:A_i(|g|)>t}\lambda_i
\le C\mu\{\mathcal M g>t\},
\qquad
\sum_i\lambda_i|A_i(g)|^2\le4C\|g\|_{L^2(\mu)}^2.
```

Decompose $`|H|=R+D_0+S`$ into its root column, marker diagonal, and
strict positive tails. The first two norms are at most one. Set

```math
\lambda_i=\frac{\mu_iM_{i+1}}{M_i},\qquad
g(j)=\frac{f(j)}{\sqrt{\mu_j}}.
```

The marker coordinates of $`S^{\dagger}f`$ are
$`\sqrt{\lambda_i}A_i(g)`$. Moreover,

```math
\lambda_k\le M_{k+1},\qquad
\sum_{i=k+1}^{m-1}\lambda_i
\le\sum_{i=k+1}^{m-1}\mu_i\le M_{k+1}.
```

The finite estimate applies with $`C=2`$. Since
$`\|g\|_{L^2(\mu)}=\|f\|_2`$, it proves

```math
\|\,|H|\,\|\le2+2\sqrt2\lt5.
```

Every averaging tail includes the terminal leaf. The estimate is uniform
in the number of leaves, their weights, and the terminal weight.
Transposition, row or column permutations, and literal sign changes
preserve this absolute-norm bound.

## 3. Packet reassociation and its exact affected space

Partition the abstract leaves of a comb into nonempty packets, in any
chosen order. On each packet apply its normalized weighted comb basis;
let U be the direct sum of those bases in the chosen leaf order. Let A
map that order back to the original abstract-leaf coordinates. Let V be
the outer comb on the packet-root coordinates, with identity on their
internal marker coordinates. A coordinate permutation C arranges the
final column labels. The packet-tree basis is $`H_{\rm pkt}=AUVC`$, and

```math
Q=H_{\rm pkt}^{\dagger}H,
\qquad
\|\,|Q|\,\|
\le\|\,|H_{\rm pkt}^{\dagger}|\,\|
   \|\,|H|\,\|\lt5^3.
```

The comb and packet tree have the same positive root vector. After
aligning its label, $`Q=1\oplus\widetilde Q`$ exactly. The factors U
and V separately need not fix that coordinate; only the complete
coefficient transition is used as the batched primitive.

An abstract leaf may stand for the scaling vector of an attached
subtree. Its internal markers are orthogonal to that scaling vector and
are unchanged by reassociation. Thus the complete transition acts only
on the marker space of the replaced path. This distinction permits
nested subtree supports while keeping the coefficient blocks disjoint.
The final residual abstract leaf is always present.

## 4. Full heavy paths and a sorting lemma

At each internal node choose a heavier child, breaking a tie by a fixed
rule. Maximal paths following these chosen edges partition the internal
nodes. Continue every path to an actual unit-mass atom; do not stop at
a balanced node or at an unresolved large subtree.

For one path, its abstract leaves are its successive side subtrees and
its terminal atom. Write their unnormalized integer masses as
$`\mu_0,\ldots,\mu_m`$ and $`M_i=\sum_{j=i}^m\mu_j`$. If
$`\delta=2^{-b}`$, the old tree gives

```math
\delta M_i\le\mu_i\le M_{i+1}\quad(i\lt m),
\qquad \mu_m=1.
```

If a heavy continuation is the first child, its old marker is the
negative of the displayed comb marker. Retain this literal column sign
in the coefficient transition; Section 2 bounds both conventions.

In particular the last side mass is one. Sort the complete list in
nonincreasing order, obtaining $`\nu_0,\ldots,\nu_m`$ with tails
$`P_i=\sum_{j=i}^m\nu_j`$. The sorted comb still satisfies

```math
\delta P_i\le\nu_i\le P_{i+1}\quad(i\lt m).
```

For the upper bound, take an occurring weight $`w>1`$ and the last
original index j with $`\mu_j\ge w`$. It is nonterminal. Every later
weight is below w, and their sum is at least $`\mu_j\ge w`$. Hence

```math
\sum_{k:\mu_k\lt w}\mu_k\ge w.
```

All those smaller weights remain after any sorted occurrence of w.
For $`w=1`$, every nonterminal sorted occurrence has at least one unit
weight after it. This proves the upper bound, including ties.

For the lower bound, take the earliest original index j whose mass is
at most w. It is nonterminal, because the last side mass and terminal
mass are both one. Its original suffix contains every weight at most w,
so

```math
\sum_{k:\mu_k\le w}\mu_k
\le M_j\le\mu_j/\delta\le w/\delta.
```

The sorted suffix beginning with any occurrence of w is no larger than
this sum. This proves the lower bound. Sorting introduces no unpriced
physical operation: it specifies the new abstract tree and the entries
of its complete coefficient transition.

## 5. A bias-halving packet pass without borrowed boundary paths

For $`b>8`$ set $`\epsilon=2^{-\lceil b/2\rceil}`$. Process the
sorted list from its beginning, leaving its final unit atom as the
terminal outer-comb leaf. At a current suffix of mass P:

- If its first mass is at least $`\epsilon P`$, make that leaf a
  singleton packet.
- Otherwise accumulate successive leaves until their total first reaches
  $`\epsilon P`$, and make them one packet.

The second case always finishes before the terminal atom. If there are
at least two leaves, their positive integer masses give
$`P-1\ge P/2\ge\epsilon P`$. Since the list is descending and the
first mass was below $`\epsilon P`$, every leaf accumulated in this
case is also below that threshold. Its packet mass A therefore satisfies

```math
\epsilon P\le A\lt2\epsilon P,
\qquad P-A>(1-2\epsilon)P.
```

A singleton packet has mass between $`\epsilon P`$ and $`P/2`$ by
Section 4. Thus all outer-comb child fractions are at least
$`\epsilon`$; all peeled fractions are at most one half.

For a nonsingleton packet, every constituent leaf has mass at least
$`\delta(1-2\epsilon)P`$. Indeed its sorted global suffix still
contains the entire positive tail after the packet. Any suffix inside
the packet has mass at most $`A\lt2\epsilon P`$, and each of its two
nonempty children contains at least one constituent leaf. Every inner
comb child fraction is consequently at least

```math
\frac{\delta(1-2\epsilon)}{2\epsilon}
\ge\frac{3\epsilon}{8}
\ge2^{-(\lceil b/2\rceil+2)}.
```

Here $`\epsilon\le1/8`$ and $`\delta\ge\epsilon^2`$ suffice.
Singleton packets have no inner markers. The terminal atom contributes
no split. Hence this pass replaces every split along the path by splits
with child fractions at least $`2^{-b'}`$, where

```math
b'=\lceil b/2\rceil+2.
```

Unlike a pass stopped at the first nonheavy subtree, this construction
has no subthreshold unresolved tail. Reaching unit atoms and sorting the
whole abstract list are both essential to that conclusion.

## 6. Global batching, recurrence, and the final atom order

For each maximal heavy path, take the complete transition from its old
comb to its sorted packet tree. The paths partition the old internal
markers. Each transition fixes its path-root scaling vector and leaves
the internal markers of every attached side subtree unchanged. Replacing
all paths therefore gives a well-defined new full binary tree on the
same physical atoms.

Assign each new path's markers bijectively to its old marker labels;
this is part of the basis definition, not a separate physical gate.
In those old and new coordinate labels, the global transition is
one root identity plus a direct sum of the path-marker blocks from
Section 3. Nested supports cause no overlap between these blocks:
the parent path only uses a side subtree's scaling vector, while the
nested path acts in its orthogonal marker space. It follows that

```math
\|\,|F_{j+1}^{\dagger}F_j|\,\|\lt125.
```

No local rotation schedule or boundary-path coloring is invoked.
Every node of the new tree belongs to one of the packet constructions,
so its global child-fraction bound is exactly the bound from Section 5.

Iterate while $`b_j>8`$, using
$`b_{j+1}=\lceil b_j/2\rceil+2`$. For those integer values,

```math
b_{j+1}-5\le\frac{b_j-5}{2}.
```

The process terminates after $`O(\log(b+2))`$ passes at $`b_s\le8`$.
The number of leaves and the tree height do not enter this recurrence.

The intermediate trees use arbitrary nested subsets of the original
atoms. Choose a single final depth-first permutation P. Under P every
subtree becomes a contiguous interval, and child masses are unchanged.
The finite-interval norm argument in
[R86 Sections 2–3](INTERVAL_DILATION_DEPTH.md#3-a-square-root-normalization-bound)
then gives a constant bound for $`B^{\dagger}PF_s`$: it uses
orthonormality, interval boundary sparsity, ancestor telescoping, and a
constant child-fraction bound, without requiring equal leaf depth.
All atom permutations have absolute norm one. Telescoping now proves the
full identity in Section 1, with the explicit final $`P^{\dagger}`$.

For the original interval embedding $`FR=EW`$, this specializes to

```math
W=E^{\dagger}P^{\dagger}B
Q_{\rm end}Q_{s-1}\cdots Q_0R.
```

The atom permutation is a physical coordinate permutation. Column-label
permutations used to describe individual coefficient blocks are already
included in the corresponding Q and must not be counted a second time.

## 7. An admissible sharpness witness for the unchanged overlap

The logarithmic factor count is compatible with an exact obstruction to
using just the original single positive-overlap representation. Let
$`s\ge4`$, $`q\ge2^s-1`$, $`\tau\ge s`$, and
$`p=2^{-s}`$. In a complete q-level tree put branch probabilities
$`(p,1-p)`$ on the first $`m=2^s-1`$ all-right-spine nodes and use
balanced pairs elsewhere. All probabilities lie exactly on the
$`2^{-\tau}`$ grid. With $`r=q\tau`$ and $`D=2^r`$, all prefix and
leaf atom masses are positive integers.

In unit-interval normalization the peeled cells have masses
$`m_t=p(1-p)^t`$, for $`0\le t\lt m`$. They obey
$`p/4\le m_t\le p`$ and cover $`[0,1/2]`$. For the first assertion,
$`\log(1-x)\ge-x/(1-x)`$ gives the lower bound. For the second,
Bernoulli's inequality gives
$`(1+1/(2^s-1))^{2^s-1}\ge2`$.

Let J be a dyadic interval inside $`[0,1/2]`$, of length
$`u\ge8p`$. At most two peeled cells intersect the endpoints of its
first half without being wholly contained there. Thus that half contains
at least $`u/(4p)`$ complete peels. If peel t is one of them, its
weighted-tree wavelet has absolute overlap with the dyadic wavelet
$`h_J`$ at least

```math
|\langle h_J,w_t\rangle|
\ge\sqrt{(1-p)m_t/u}\ge\sqrt{p/(8u)}.
```

There is no cancellation: after choosing the custom wavelet's positive
sign on its peel, its constant negative tail begins in J's first half,
where the integral of $`h_J`$ over the remaining tail is nonpositive.

Test $`|B^{\dagger}EW|`$ on the unit vector constant on these m
spine-marker columns. Each row J contributes at least $`u/128`$ to
the squared output norm. Each scale has $`1/(2u)`$ such rows, and the
$`s-3`$ scales from $`u=1/2`$ through $`u=8p`$ are distinct. Therefore

```math
\|\,|B^{\dagger}EW|\,\|\ge\frac{\sqrt{s-3}}{16}.
```

This is an actual rounded-tree family. For fixed L, choose
$`n=2^{2s}`$ and the R86 precision $`\tau=4s+2L+12`$. Its largest
thin bin has width at least $`h/[8(\tau+2)]`$, where
$`h=n-O_L(\log n)`$, and eventually contains the required spine.
For completeness, under the backward rule
$`q=\min\{b,\lfloor(h-b+m_0)/(c\kappa)\rfloor\}`$, with
$`h\ge4c\kappa`$ and $`m_0\ge2c\kappa`$, some bin has width at
least $`h/(4c\kappa)`$: at the first endpoint at most $`h/2`$,
either the preceding bin already crossed below $`h/4`$, or the next
step from $`[h/4,h/2]`$ has that width. Take
$`c=2,\kappa=\tau+2`$ for R86. Its certified rounding preserves the
specified exact grid probabilities.

The witness gives an $`\Omega(\sqrt\tau)`$ norm along admissible
instances. It is neither a native-depth lower bound nor an
approximation-robust obstruction. For an exact factorization of the
unchanged $`K=B^{\dagger}EW`$ whose factors all have absolute norm at
most a fixed $`C>1`$, entrywise
submultiplicativity implies $`C^{s+O(1)}\ge c\sqrt\tau`$, requiring
$`s=\Omega(\log\tau)`$ on this family. This gives the logarithmic lower
scale for that fixed-reference interface. The construction in this note
also changes the physical atom order; the witness does not establish
optimality for that larger interface, since the additional conjugate
$`B^{\dagger}P^{\dagger}B`$ need not have bounded absolute norm.

## 8. Native scope and evidence boundary

Completing the original interval tree by integer midpoint refinement
inside each leaf gives initial $`b\le\max\{\tau,2\}`$. The exact
factorization thus supplies $`O(\log(\tau+2))`$ square atom factors
with uniformly bounded absolute norms, plus one physical atom
permutation. The native square row/column-source wrapper must still be
applied to every selected factor, with its guards, actual inverse,
anchor, approximation error, and returned workspace.

In particular, full atom factors use two r-bit indices and table address
$`a+2r`$, unlike the rectangular R86 overlap's $`a+q+r`$. The exposed
suffix must pay for that address margin and constant-many atom banks.
The [selected routing proof](SELECTED_DIRTY_ROUTING.md) supplies the
ragged-table and dirty-return interfaces; it does not make the factors
free. A factorization alone is not an unqualified compiler theorem.

The present claims are finite analytic identities and inequalities.
The comb maximal estimate, global sorting and packet pass, exact
permutation placement, and admissible norm witness have separate
internal mathematical checks. No scalable native emitter, efficient
classical preprocessing bound, or new scientific receipt is claimed.
