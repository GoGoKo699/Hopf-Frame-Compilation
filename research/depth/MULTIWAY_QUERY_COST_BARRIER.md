# The count cost of replacing two indicators by many

[Native bilinear query R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Conditional frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md) ·
[Paired permutation wrappers](PAIRED_PERMUTATION_ROUTING.md)

Splitting an address into many short fields appears to reduce the depth
of each dirty indicator. For a separated multilinear middle word, the
price is the number of product terms required to encode a general table.
This chapter gives that price and checks two other proposed depth
shortcuts. It proves no unrestricted native T-count or T-depth lower
bound.

## 1. The separated-reader architecture

Partition an a-bit address into k nonempty fields with lengths
$`a_1,\ldots,a_k`$, where $`a=\sum_i a_i`$. Put

```math
M_i=2^{a_i},\qquad Q=\prod_i M_i=2^a.
```

The proposed reader uses a dirty one-hot bank Y_i of length M_i for
each field. Its central Boolean polynomial is a sum of r separated
products,

```math
F(Y_1,\ldots,Y_k)
=\bigoplus_{\nu=1}^{r}
       \prod_{i=1}^{k}\langle c_{\nu,i},Y_i\rangle,
\qquad c_{\nu,i}\in\mathbb F_2^{M_i}.                 \qquad(1)
```

Arbitrary invertible linear coordinate changes within individual banks
are allowed and absorbed in the coefficients. On one-hot arguments the
represented table is

```math
g(x_1,\ldots,x_k)
=\bigoplus_{\nu=1}^{r}\prod_{i=1}^{k}c_{\nu,i}(x_i).
                                                               \qquad(2)
```

Thus r is the binary separated tensor-rank budget. Allowing an affine
form on a one-hot bank does not enlarge (2): its constant is represented
by the all-ones coefficient vector. Zero factors pad a shorter expression
to r terms.

The usual complete dirty extraction takes a finite difference in each
bank. Multilinearity gives

```math
\bigoplus_{S\subseteq\{1,\ldots,k\}}
F\bigl(Y_i+\mathbf1_{i\in S}e_{x_i}:1\le i\le k\bigr)
=F(e_{x_1},\ldots,e_{x_k}).                           \qquad(3)
```

All unknown dirty offsets cancel in (3). For a scalar phase the output
is $`(-1)^g`$; including an arbitrary output bit h as a phase control
and conjugating h by Hadamards gives the XOR query. A literal expansion
of (3) has $`2^k`$ central appearances and all required indicator
inverses. The argument below optimistically ignores that multiplicity,
indicator costs and the cost of synthesizing a k-factor controlled
phase. Those omitted costs can only make this particular route more
expensive.

The scope of (1)–(2) matters. A general Clifford+T word need not admit
such a decomposition with r bounded by its T-count. A fused native
implementation could also be cheaper than emitting the terms of (1)
separately. No inequality $`r\le T`$ is asserted for such circuits.

## 2. A counting bound on the central representation

For the fixed address partition of section 1, there are at most

```math
2^{r\sum_i M_i}                                      \qquad(4)
```

coefficient lists in (1), hence at most that many distinct Boolean
tables representable by r terms. Different lists can give the same
table; this only strengthens the bound. There are $`2^Q`$ Boolean
tables on the chosen address set. Covering all of them therefore requires

```math
r\sum_i M_i\ge Q.                                   \qquad(5)
```

More quantitatively, the fraction of tables having such a representation
is at most $`2^{r\sum_iM_i-Q}`$. This is a worst-case counting
statement, not a rank lower bound for every supplied table.

If the strategy chooses its partition from a family $`\mathcal P`$
depending on the table, there are at most $`k^a`$ labelled partitions
into k fields. Union counting therefore gives the necessary condition
$`r\max_{\pi\in\mathcal P}\sum_i M_i(\pi)+a\log_2k\ge Q`$.
For a fixed size profile this changes (5) only by the additive
$`a\log_2k`$. Since $`k\le a`$ and $`a\log_2a=o(Q)`$, this
allowance does not change the asymptotic consequences below, including
the conclusion when every permitted field has sublinear length.

If all fields have equal size $`a_i=a/k`$ (when k divides a), (5)
gives

```math
r\ge\frac1k Q^{1-1/k}.                              \qquad(6)
```

For any fixed $`k\ge3`$ this exceeds $`\sqrt Q`$ by an unbounded
factor. In contrast, the two-field matrix representation needs at most
$`\min\{M_1,M_2\}`$ rank-one terms, and R91 implements its
balanced bilinear middle at square-root T-count. The comparison already
uses that simplified constant-depth bilinear baseline.

For unequal fields let $`a_{\max}=\max_i a_i`$. Then
$`\sum_iM_i\le k2^{a_{\max}}`$. If a strategy restricts its
term budget to $`r\le C\sqrt Q`$, (5) forces

```math
a_{\max}\ge a/2-\log_2 k-\log_2 C
          \ge a/2-\log_2 a-\log_2 C.                 \qquad(7)
```

Thus a separated middle with square-root term budget must retain at
least one field containing a linear fraction of all address bits.
Making every indicator address sublinear cannot preserve that term
budget. Using the retained indicator compiler on the necessary long
field still gives its logarithmic address-depth allowance. This last
sentence is a statement about that chosen compiler, not a lower bound
on every possible implementation of an indicator.

### Translating the term budget into an actual resource claim

For an implementation that physically emits each nonzero separated
product as its own native controlled-phase word, with a constant
positive number of T gates charged per term and no cross-term
cancellation, an $`O(\sqrt Q)`$ T-count permits only
$`O(\sqrt Q)`$ such terms. Then (7) applies. For k growing, a
term may require $`\Omega(k)`$ gates or a larger cost; granting it
unit cost in this scout is deliberately favorable.

This gate-count implication is restricted to that term-by-term
implementation rule. It does not apply to a compiler that synthesizes
the whole tensor polynomial jointly without this charging rule.
Likewise, the bound does not cover nonseparated quantum branches,
intermediate nonlinear coordinate changes, general interference, or a
different representation of the table. Such a mechanism would have to
provide its own complete native circuit and returned-work proof.

## 3. Coherent Boolean-branch averaging does not remove this cost

The same representation bound survives one natural approximation
attempt. Allow an arbitrary address-independent pure seed
$`\sum_s\alpha_s|s\rangle`$, even with unlimited initialized seed
space. Each branch implements exactly the scalar Boolean phase of
some table $`g_s`$ represented by at most r terms in (2), and the
actual seed preparation is then reversed. All other work returns.

If the full initialized-isometry error from the literal phase
$`(-1)^{f(x)}`$ is at most $`\varepsilon\lt\sqrt2`$, then for
every x, norm preservation by the final actual inverse gives

```math
4\sum_s|\alpha_s|^2\mathbf1_{g_s(x)\ne f(x)}
\le\varepsilon^2.                                   \qquad(8)
```

This also applies to XOR-query branches by choosing the arbitrary
output in $`|{-}\rangle`$. Put
$`\delta=\varepsilon^2/4\lt1/2`$. Averaging (8) over the Q
addresses shows that at least one branch table is within relative
Hamming distance delta of f. A Hamming ball of that radius contains
at most $`2^{H_2(\delta)Q}`$ Boolean tables, where
$`H_2(\delta)=-\delta\log_2\delta-(1-\delta)\log_2(1-\delta)`$.
Together with (4), universal coverage implies

```math
r\sum_i M_i\ge[1-H_2(\delta)]Q.                     \qquad(9)
```

At fixed error the factor in brackets is positive and independent of
a. Thus (6)–(7) retain their asymptotic consequences. The seed and its
distribution may depend on the prescribed table, but not on the live
address. No seed-preparation cost is used in the lower estimate.
Allowing table-dependent partitions adds the same $`a\log_2k`$
capacity allowance described after (5), with the maximum taken over
permitted partitions; the asymptotic conclusions again remain.

This argument requires Boolean phase branches followed by actual
unpreparation. It does not cover address-dependent seed preparation,
branches with unreturned garbage, non-Boolean intermediate phases, or
interference that does not have this branch form. In particular it is
not a lower bound for arbitrary approximate Clifford+T queries.

## 4. Other depth shortcuts and their unpaid operations

### Routing a fixed bit flip to the selected site

An indicator has the identity
$`I_x=P_xX_0P_x^{-1}`$, where $`P_x`$ is an appropriate XOR
permutation of the $`2^a`$ output sites. Substituting the new paired
permutation wrapper does not by itself accelerate it. That router has
$`p=2^a`$ selectable permutations and its sparse-grid shear calls
an a-bit dirty indicator. The supposed construction invokes the same
indicator it is meant to improve. The direct grid also has
$`O(pw)`$ count and dirty volume, with $`w=2^a`$ here. Neither
the recursive identity nor extra supplied width removes its T-count.

### Replacing a selected site by a small linear-hash bucket

A scalar zero-test can have few hash collisions at each fixed nonzero
input. That observation does not justify replacing the complete
indicator XOR by a hash bucket. Let a branch s use a linear map
$`H_s:\mathbb F_2^a\to\mathbb F_2^{d_s}`$ of rank strictly less
than a and toggle every output site in

```math
B_s(x)=\{y:H_sy=H_sx\}=x+\ker H_s.
```

This bucket contains at least two sites. On the legitimate arbitrary
output input $`|0^Q\rangle`$, the target produces the computational
basis state $`|e_x\rangle`$, whereas this branch produces
$`|\mathbf1_{B_s(x)}\rangle`$. They are orthogonal for every x
and every such s.

To include coherent seeds, compare before the actual inverse of seed
preparation. The ideal and proposed normalized states are respectively

```math
\sum_s\alpha_s|s\rangle|e_x\rangle,
\qquad
\sum_s\alpha_s|s\rangle|\mathbf1_{B_s(x)}\rangle.
```

Their inner product is zero, so their norm distance is exactly
$`\sqrt2`$, even after optimizing a common phase. Actual seed
unpreparation preserves this distance. The argument also permits
address-dependent seed amplitudes, because x is fixed in the
comparison. If only a seed weight p uses rank-deficient hashes and the
other branches produce the correct singleton, the squared distance is
$`2p`$. Arbitrary unknown output data need not be analyzed further:
the valid input $`|0^Q\rangle`$ already witnesses failure of the
full operator contract.

This excludes only the literal bucket-toggle shortcut. A different
interference scheme, or a procedure that resolves the collisions before
returning the seed, must provide and charge that additional operation.

### Applying all sequential gates on borrowed slots at once

For fixed equal-width unitaries $`U_0,\ldots,U_{t-1}`$, let
$`A=\bigotimes_i U_i`$ act in parallel on arbitrary slots. Define
the fixed cyclic slot permutation P so that
$`P(\bigotimes_i B_i)P^{-1}=\bigotimes_i B_{i-1}`$, with indices
modulo t. The exact commutator is

```math
APA^{-1}P^{-1}
=\bigotimes_i U_iU_{i-1}^{-1}.                       \qquad(10)
```

It gives neighboring ratios on **every** slot, not the sequential
product on one designated slot with identity on the others. Requiring
all ratios except one to be identity forces all corresponding U_i to
agree; a lone nontrivial ratio cannot be obtained around the cycle.
One can choose prefix-product factors to prescribe a chain of ratios,
but synthesizing those factors then includes the original sequential
products. Alternatively iterating a shift-and-gate word t times uses t
charged appearances of the gate array. These exact identities do not
flatten an arbitrary sequential query compiler.

The paired routing lemma cancels its own inverse-permuted auxiliary
array around a body. It does not cancel the unwanted gate products in
(10). A different cancellation identity remains possible, but must
return every occupied slot and price the actual array factors.

## 5. Resumption boundary

These candidates do not provide the smaller query-depth function
required by R98. The multiway route would need a native central reader
that escapes the separated term cost without merely rebuilding a long
indicator. The borrowed-slot route would need a complete return identity
whose array factors are actually cheaper to compile than the original
prefix products. Neither missing operation is supplied by these notes.

The counting argument and exact operator identities are analytic. No
finite-rank experiment is used as evidence for a native lower bound.
R96's unconditional general upper law, R98's conditional status and
the continuing optimal-depth objective are unchanged.
