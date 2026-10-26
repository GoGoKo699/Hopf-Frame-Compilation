# A complete probe from ordered matchings

[Weak pair gauges](COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge) ·
[Complete probe readout](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) ·
[Returned cubic batches](DIRTY_CUBIC_PHASE_COMPILER.md)

The complete lookup probe does not force the exterior-algebra pair
representation. A different representation retains just one matching for
each even subset and one distinguished singleton for each odd subset.
It has the correct probe on every carrier and address input, with no
Pfaffian sums in its matrix entries. Its exact native baseline below has
linear address depth; a sublogarithmic joint compiler is still missing.
This is a new sufficient source, not an improvement to R96.

## 1. An associative ordered-matching algebra and module

Fix ordered labels $`1,\ldots,p`$, with $`p\ge2`$,
$`m=2^p`$ and $`q=\binom p2`$. Introduce commuting variables
$`z_{ij}`$ for $`i\lt j`$. Form the binary monomial quotient algebra

```math
\mathcal A_p=\mathbb F_2[z_{ij}:i\lt j]/\mathcal I,
\qquad
\mathcal I=\langle z_{ij}^2,\ z_{ij}z_{kl}:
 [i,j]\cap[k,l]\ne\varnothing\rangle.                 \qquad(1)
```

Intervals here are closed real intervals. Thus a product is zero when
two edges share an endpoint, cross, or nest. A surviving monomial is a
collection of separated intervals. It is uniquely specified by its
even endpoint subset S: list S increasingly and pair consecutive labels.
Write $`M(S)`$ for that matching. Every even S occurs once, including
the empty subset. Associativity is inherited from the polynomial ring;
this is not a rule for choosing a matching afresh during multiplication.

For each i, take a cyclic module $`u_i\mathcal A_p`$ and additionally
annihilate every $`z_{jk}`$ with $`j\le i`$. Its surviving basis is
$`u_i\prod_{e\in M(S\setminus\{i\})}z_e`$ for odd subsets S whose
minimum is i. Use the direct sum

```math
\mathcal M_p=\mathcal A_p\ \oplus\
 \bigoplus_{i=1}^p u_i\bigl(\mathcal A_p/
                    \langle z_{jk}:j\le i\rangle\bigr).
                                                               \qquad(2)
```

Its basis has exactly one vector for every subset S, hence dimension m.
The empty-subset vector belongs to the even component. A singleton i is
the cyclic vector of its odd component. These are labels for physical
carrier coordinates, not initialized vectors supplied to a circuit.

Let $`N_{ij}`$ be the row-oriented multiplication matrix of $`z_{ij}`$:
its row U has a one in column $`U\cup\{i,j\}`$ exactly when adding that
edge to the basis monomial for U survives. Equivalently, it removes the
edge ij from the canonical matching of a column subset. Consequently

```math
N_eN_f=N_fN_e,\qquad N_e^2=0,\qquad
H_s=\prod_e(I+s_eN_e),\qquad H_sH_t=H_{s+t},\quad H_s^2=I.
                                                               \qquad(3)
```

Every factor is invertible. The carrier map $`c\mapsto H_sc`$ is thus
a phase-free permutation of every m-bit carrier, for every arbitrary
q-bit pair key s. The group law is exact on all these inputs.

## 2. One canonical monomial gives each lookup coefficient

Let $`e_S^{\mathsf T}`$ denote a carrier coordinate row, and put

```math
\ell_x=e_\varnothing^{\mathsf T}+\sum_i x_i e_{\{i\}}^{\mathsf T},
\qquad
v_x=\sum_{S\subseteq[p]}\left(\prod_{i\in S}x_i\right)e_S^{\mathsf T},
\qquad K(x)_{ij}=x_ix_j.
```

The empty row of $`H_s`$ has coefficient
$`\prod_{e\in M(S)}s_e`$ in column S when S is even, and zero when S
is odd. Singleton row i has coefficient
$`\prod_{e\in M(S\setminus\{i\})}s_e`$ when S is odd with minimum i,
and zero in all other components. There is exactly one way to obtain
each displayed basis monomial in (3). Therefore

```math
\ell_x H_{K(x)}=v_x.                                    \qquad(4)
```

For even S, substituting $`s_{ij}=x_ix_j`$ gives its full address
monomial. For odd S, its unique singleton contributes the remaining
factor. No sum over possible perfect matchings or possible singletons
is needed.

Use the [paid root-only completion](COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
$`L_x`$, which changes only the empty carrier bit by
$`c_\varnothing\mapsto c_\varnothing+\sum_i x_ic_{\{i\}}`$.
It is an involution with empty row $`\ell_x`$. Define, in rightmost-first
operator order,

```math
V_x=H_{K(x)}L_x,\qquad
V_x^{-1}=L_xH_{K(x)},\qquad
VZ_{c_\varnothing}V^\dagger
 =(-1)^{\sum_S c_S\prod_{i\in S}x_i}=Q.                \qquad(5)
```

The inverse empty row is (4). This proves the complete probe, including
coherent x, every carrier column, and any external reference. The order
in (5) matters: the new representation need not commute with $`L_x`$
or with the full linear factor used by the exterior representation.

The representation is genuinely different already at p=4. Its empty
row in the top column is $`s_{12}s_{34}`$, whereas the exterior
representation has
$`s_{12}s_{34}+s_{13}s_{24}+s_{14}s_{23}`$. For key
$`s_{13}=s_{24}=1`$ and all other pair bits zero, those entries differ.
Both give the required top monomial on clique keys. This is not just a
fixed carrier-coordinate change with pair labels held fixed: here
$`N_{13}N_{24}=0`$, whereas the corresponding exterior product is
$`N_1N_3N_2N_4\ne0`$. Simultaneous similarity would preserve that product.
No claim about arbitrary pair-key relabeling is needed. Equation (5)
continues to hold on every original address and arbitrary carrier.

## 3. The arbitrary pair-key echo and its exact native baseline

Introduce only q arbitrary pair bits s, and let
$`C|s,c\rangle=|s,H_sc\rangle`$. The controlled pair toggle
$`Q_K|x,s,c\rangle=|x,s+K(x),c\rangle`$ gives

```math
C Q_K C:\quad (s,c)\longmapsto(s+K(x),H_{K(x)}c).
                                                               \qquad(6)
```

Thus the physical source $`C Q_K C L_x`$ has (5)'s probe and an
arbitrary auxiliary translation. Its actual inverse cancels that
translation inside the completed probe. No pair bit is initialized.

### Each interval family with a common point is one quadratic shear

For a fixed edge e, all nonzero entries of $`N_e`$ pair distinct carrier
coordinates; no coordinate is both source and target. The controlled
factor $`I+s_eN_e`$ is their shared-control CCX matching.

If a family of edge intervals has a common point, then
$`N_eN_f=0`$ for every two edges in the family. Moreover there is no
two-step carrier path: a nonzero product of these partial matching
matrices would contain that unique path with coefficient one. Hence
the union of all sources is disjoint from the union of all targets.
The entire family is one quadratic XOR shear
$`I+\sum_e s_eN_e`$, rather than an inference from completed commutation.
Hadamards on its distinct targets turn it into a cubic phase batch.
For R raw CCX occurrences and t distinct targets, the returned-dirty
compiler gives

```math
T\le8R,\qquad D_T\le8,\qquad
G_{\rm all}\le16R+2t,\qquad w_{\rm private}\le R.       \qquad(7)
```

Each of the R helpers is arbitrary and returns, with literal phase one.
The actual reversed word has the same bounds. Shared controls and
targets are included in this batch statement.

### A complete paid schedule

Let $`h=\lfloor p/2\rfloor`$. Assign each interval to the first even
integer it contains. Every interval between distinct integer labels
contains such an integer, and this partitions all edges into h
common-point families. Execute their shears sequentially. The total
number of raw occurrences has the exact count

```math
E=\sum_{S\subseteq[p]}\lfloor|S|/2\rfloor
  =(p-1)2^{p-2}.                                        \qquad(8)
```

Indeed every canonical matching contributes one removable edge per
pair, with the singleton pivot preserved on odd subsets. The second
identity uses the total subset size $`p2^{p-1}`$ and the
$`2^{p-1}`$ odd subsets. Summing (7), with at most m targets per family,
gives a complete exact implementation of C and its actual inverse:

```math
T_C\le8E,\qquad D_T(C)\le8h,\qquad
G_C\le16E+2mh,\qquad w_C\le E.                         \qquad(9)
```

Reuse the returned private pool between families. Combining the
eight-band, 8q-T pair toggle and the helper-free four-band root
completion gives the literal source in (6), followed in operator
notation by $`L_x`$ on the right, with

```math
\begin{aligned}
T_V&\le16E+8q+6p+(p\bmod2),\\
D_T(V)&\le16h+12,\\
G_V&\le32E+4mh+18q+16p+(p\bmod2)+5,\\
w_{\rm extra}&\le q+\max(E,q),\qquad
a_{\rm clean}=0,\qquad\epsilon=0.
\end{aligned}                                                   \qquad(10)
```

Extra width excludes the address and principal carrier. It includes
the q pair bits and one returned pool reused by C and $`Q_K`$;
$`L_x`$ needs no helper. The same common ledger applies to the actual
inverse. The complete two-call probe readout has its usual additional
Clifford walls and doubles source T-count and T-depth. This schedule
has $`T,w=O(p2^p)`$ and $`G=O(p2^p)`$, within the discovery budgets,
but $`D_T=O(p)`$. It is not a faster indicator than the retained
logarithmic-depth construction.

The first unpaid operation is a faster native synthesis of the joint
word in (6), or any source with its complete probe, using this different
representation. Distinct interval families may have nonzero products;
their completed commutation does not parallelize their native words.
Restricting to subsets formed from disjoint consecutive pairs leaves
the usual subset-product action on $`\lfloor p/2\rfloor`$ effective
generators. That identifies a remaining serial construction problem;
it is not a native-depth lower bound. Arbitrary nonmonomial dynamics
and later address mixing remain eligible.

The [ordered-matching fixtures](../../tests/test_ordered_matching_probe.py)
check the quotient/module matrices, exact source order, arbitrary pair
keys, the first difference from the Pfaffian representation, and an
actual returned-dirty native batch. Their finite checks do not replace
the all-size algebra or supply the missing shallow compiler.

<a id="minimal-tree-source"></a>

## 4. A smaller tree source has linear count and no private work

A further specialization removes the pair representation itself. It is a
sparse version of the [affine-pencil source](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md),
with one coefficient carrier, no extra field coordinates and no Vandermonde
change. Its native baseline still has linear address depth. The gain is in
paid count and width, not in the depth frontier.

For every nonempty subset S, let $`i(S)=\min S`$ and
$`\pi(S)=S\setminus\{i(S)\}`$. On coefficient columns define

```math
T_xe_S=x_{i(S)}e_{\pi(S)},\qquad T_xe_\varnothing=0,
\qquad A_x=I+T_x.
```

There are exactly $`m-1`$ directed edges. Every edge lowers subset size, so
$`T_x^{p+1}=0`$ and

```math
A_x^{-1}=\sum_{k=0}^{p}T_x^k,
\qquad e_\varnothing^{\mathsf T}A_x^{-1}
 =\sum_{S\subseteq[p]}\left(\prod_{i\in S}x_i\right)e_S^{\mathsf T}.
                                                               \qquad(11)
```

Indeed each S has exactly one path to the root: delete its labels in increasing
order. Its path weight is the full address monomial on S. The direct source

```math
V|x,c\rangle=|x,A_xc\rangle,
\qquad VZ_{c_\varnothing}V^\dagger=Q                  \qquad(12)
```

therefore has the complete probe. V, rather than its inverse, is the quadratic
forward map: each carrier bit changes by a sum of $`x_i c_S`$ terms. All x and
c inputs are arbitrary, and (12) tensors with every reference. There is no
pair key, root-completion word or private initialized coefficient.

### Literal native schedule and actual inverse

Let $`N_i`$ delete i only when i is the smallest label of its column subset.
Set $`B_i=I+x_iN_i`$. Each B_i is a disjoint carrier matching of
$`r_i=2^{p-i}`$ shared-control Toffolis. A product $`N_iN_j`$ vanishes whenever
$`i\le j`$: after deleting j, the remaining minimum is larger than j.
Consequently, in rightmost-first operator notation,

```math
A_x=B_1B_2\cdots B_p.                                  \qquad(13)
```

The emitted chronological order is p down to 1. Reversing that actual native
word gives $`A_x^{-1}`$ and its long path products; it is not replaced by a
fresh quadratic query. The [shared-control compiler](COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge)
pays each matching in four T bands with $`6r_i+(r_i\bmod2)`$ T gates,
$`12r_i`$ non-S Cliffords, and a reduced root power of S.

Since the $`r_i`$ are powers of two and sum to $`m-1`$, only $`r_i=1`$
contributes an odd T count. Only $`r_i=2,4`$ contribute nonzero root S powers.
After reducing adjacent S powers modulo four, their total count is at most
three forward and five in the actual inverse. Safe common bounds for
$`p\ge1`$ are therefore

```math
\begin{aligned}
T_V&=6m-5,\qquad D_T(V)\le4p,\\
G_{\rm Cliff}(V)&\le12(m-1)+5,\qquad G_V\le18m-12,\\
w_{\rm private}&=0,\qquad a_{\rm clean}=0,\qquad\epsilon=0.
                                                               \qquad(14)
\end{aligned}
```

Every basis column has literal phase one. The source and its actual inverse
use no private work. Two source appearances give the coefficient-lookup
phase by (12), with the usual complete-indicator readout. For an ordinary
addressed bank, also pay the fixed Möbius bank change and its inverse,
$`O(pm)`$ Cliffords; they add no T gates or dirty work.
The [complete-probe count bound](NATIVE_INDICATOR_SHARING.md#complete-probe-support-count)
makes (14) order-optimal in T-count for this interface at $`p\ge2`$, with
zero private work; no optimal constant is asserted. Its depth remains linear,
whereas the retained general indicator already has logarithmic depth. Thus
(14) is a count/width alternative, not an improved complete-frame upper bound.

### The full forward XOR query is also paid

On a separate arbitrary m-bit output t, the complete query

```math
Q_A:\ |x,c,t,w\rangle\longmapsto|x,c,t+A_xc,w\rangle
```

adds c by m CNOTs and applies all $`m-1`$ bilinear edge terms to t.
The edge parents are precisely the $`m/2`$ subsets not containing label 1.
Hadamards on those output coordinates turn the edge terms into one cubic
phase batch with $`m-1`$ returned dirty helpers. The eight-band compiler and
both Hadamard walls give common forward/actual-inverse bounds

```math
\begin{aligned}
T_{Q_A}&\le8(m-1),\qquad D_T(Q_A)\le8,\\
G_{Q_A}&\le18m-16,\qquad w_{Q_A}\le m-1,\\
a_{\rm clean}&=0,\qquad\epsilon=0.                     \qquad(15)
\end{aligned}
```

The separate output t is excluded from $`w_{Q_A}`$ and must be allocated
when this query is used inside another source. This is exact on every
output and helper input, but does not implement A in place or supply its
inverse for free. A constant number of these affine-in-x linear shears, interleaved with
fixed binary carrier maps, has bounded address degree; broader Gaussian
restrictions retain their separately stated hypotheses. A faster in-place word or
nonmonomial complete readout still needs its own native construction.
The [tree fixtures](../../tests/test_minimal_tree_probe.py) check the unique-path
identity, source orientation, literal native schedule and actual inverse,
and the complete forward query with arbitrary dirty helpers.
