# A uniform precision-depth bound and a stronger low-precision range

[Blocked bilinear query](BLOCKED_BILINEAR_LOOKUP.md) · [Unary group compiler](UNARY_PHASE_GRADIENT.md) · [Existing all-precision hybrid](PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision) · [Capped source precision](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)

The blocked query can retain its depth bound while using unequal
indicator lengths to retain the frame's square-root precision
dependence in T-count. This changes the workspace allocation, not its native circuit
identities. The resulting tail schedule combines with the unary early
groups when precision is at most a small constant times log n. At larger
precision, the existing hybrid's logarithmic term is already absorbed
by its nL term. This gives one uniform theorem, with constants independent
of the requested accuracy.

## 1. The uniform theorem and its stronger low-precision corollary

Let $`n\ge1`$, $`N=2^n`$, and $`0\lt\eta\le1/64`$. Put

```math
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad B_0=L+n+7.
```

For every $`b\ge17B_0`$, every prescribed complete real Hopf frame W
has one coherent Clifford+T circuit V with two external clean flags and
at most b arbitrary dirty qubits such that

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

All implied constants are absolute: none depends on n, L, eta, or b.
All bounds hold on the same circuit. The error is the full
initialized-isometry norm, including arbitrary dirty inputs, work return,
and reference correlations. Every quantum query and actual inverse is
charged. There are no measurements, resets, QRAM, supplied phase states,
or uncharged compiler primitives. T-depth permits arbitrary intervening
Clifford circuits; their elementary gate count remains included in G.

The construction has a stronger bound throughout the explicit range

```math
6\le L\le\frac1{16}\log_2(n+2):
```

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+n\right).
```

These constants are also absolute. The range is an asymptotic sufficient
condition, not a practical crossover estimate. In particular it contains
every fixed accuracy for sufficiently large n. Unlike a fixed-accuracy
statement with an unspecified eta-dependent constant, this corollary
controls the constants while L varies in the displayed range.

The inherited lower bounds give simultaneous matching count and depth
on each nonempty interval

```math
17B_0\le b\le\sqrt{N/n}
```

for the uniform theorem, and on the larger interval

```math
17B_0\le b\le\sqrt{NL/n}
```

in the displayed low-precision range. On either relevant interval,

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

These are worst-case bounds over prescribed complete real frames. No
unrestricted lower bound for the additive n or nL is asserted. The
large-width depth optimum and the high-precision complete-frame count
endpoint remain separate questions; the uniform theorem retains nL
in its general T-count. The sufficient condition $`L\le N/n^2`$
absorbs that term into $`\sqrt{NL}`$, making the count optimal in
order at every eligible width in that precision range.

## 2. An unequal indicator allocation preserves the block-depth bound

The [blocked bilinear construction](BLOCKED_BILINEAR_LOOKUP.md#1-a-controlled-bilinear-output-with-one-arbitrary-dirty-helper)
implements a controlled rank-reduced bilinear output with one arbitrary
dirty helper. A two-pass dirty traversal selects its high-address block,
and the usual four-call, two-indicator echo extracts the selected table
entry. All those exact native words, phases, and arbitrary-input return
identities remain unchanged here.

For a Q-row table, $`Q=2^r`$, with $`m\ge1`$ output bits, allocate
H and J indicator outputs and $`K=Q/(HJ)`$ high-address blocks. All
three sizes are powers of two. The generic proved ledger is

```math
T=O(Km\min(H,J)+(H+J)P),\qquad
G=O(Qm+(H+J)P),\qquad
D_T=O(Km+D_H+D_J),
```

where $`P=(a+2)^3`$ for positive integer chunk cap a, and
$`D_H,D_J`$ are the two chunked indicator depths. Let the absolute
constant $`c_I\ge1`$ bound each indicator's T-count, Clifford count,
and dirty width, including its output word, by its output length times
$`c_IP`$. A sufficient simultaneous extra width is

```math
c_I(H+J)P+p+1,\qquad p=\log_2K.
```

The last p bits are arbitrary dirty traversal selectors and the final
bit is the separate controlled-bilinear helper. Keep all these registers
disjoint from the address, output, and existing compiler reservation.
Assume the same sufficient extra-width condition

```math
B\ge16c_I(P+r+1).
```

For $`Q\ge2`$, choose J as the largest power of two at most

```math
\min\!\left\{Q,\sqrt{Qm},\frac{B}{8c_IP}\right\},
```

and then set

```math
H=\min\{J,Q/J\},\qquad K=Q/(HJ).
```

Both H and K are valid powers of two, $`H\le J`$, and $`HJ\le Q`$.
Thus these define disjoint low-address fields of lengths
$`\log_2H,\log_2J`$ and a remaining high field of length p. The
indicator pools use at most B divided by four; $`p+1\le r+1`$
fits the remaining reservation. For $`Q=1`$, emit its sole row as a
Clifford X word with no work.

Since $`H\le J`$, the selected bilinear contribution is
$`KmH=Qm/J`$. Power-of-two rounding yields

```math
\frac{Qm}{J}
=O\!\left(m+\sqrt{Qm}+\frac{QmP}{B}\right).
```

The indicators have $`(H+J)P\le2JP\le2\sqrt{Qm}P`$ cost.
Consequently the same exact query has

```math
T=O\!\left(m+\sqrt{Qm}P+\frac{QmP}{B}\right),\qquad
G=O(Qm+\sqrt{Qm}P).
```

The standalone +m count term is needed when $`m>Q`$: then the
cap $`J\le Q`$ prevents $`Qm/J`$ from falling below m. It will
be charged by the unweighted tail sum.

The count cap preserves the claimed block-depth order. If
$`J\ge\sqrt Q`$, then $`H=Q/J`$ and $`K=1`$. Otherwise
$`H=J`$. If a count cap determines J, then
$`J\gt\sqrt Q/2`$, so $`K\lt4`$. If the width cap determines
J, then $`J\gt B/(16c_IP)`$ and
$`K\lt256c_I^2QP^2/B^2`$. This proves in all cases

```math
K=O\!\left(1+\frac{QP^2}{B^2}\right),\qquad
D_T=O\!\left(m\left[1+\frac{QP^2}{B^2}\right]+D_H+D_J\right).
```

No initialized indicator, copied dirty-bit promise, or new native
controlled gate is introduced. Unequal indicator address lengths are
allowed by the generic echo and by the chunked-indicator theorem. The
matrix-coordinate changes are still charged separately for every
block and output bit, so G retains its $`Qm`$ term.

## 3. Uniform tail bounds at the original capped precisions

Write $`L'=L+2`$, the precision parameter for accuracy eta divided
by four. On a tail layer k steps from the leaves, use the original
full-input operator source at

```math
m_k=L'+4+\min\{k,\lceil\log_2(8n)\rceil\},\qquad
Q_k=4N2^{-k},\qquad r_k=n-k+2.
```

Let the tail have at most M layers. Set

```math
a_k=\min\{n+2,2^{\lfloor k/12\rfloor}\},\qquad
P_k=(a_k+2)^3.
```

Then

```math
P_k\le(n+4)^3,\qquad P_k\le27\,2^{k/4},\qquad
m_k\le L'+4+k.
```

Use Section 2 for the entire m-bit output word. The two indicators are
shared by all its output-bit operations. Each indicator address length
is at most $`r_k\le n+1`$, even with unequal H and J. Its depth
therefore obeys the same bound as in the balanced case,

```math
D_H+D_J=O\!\left(n(k+1)2^{-k/12}+\log(n+2)\right).
```

When the exponential chunk cap is smaller than the address length,
its inverse is within a fixed factor of $`2^{-k/12}`$ and its
logarithm is $`O(k+1)`$. Otherwise one chunk suffices and its depth
is logarithmic in the address length. A zero-bit indicator is Clifford.
Thus each completed tail query has

```math
T_k=O\!\left(m_k+\sqrt{Q_km_k}P_k+
                         \frac{Q_km_kP_k}{B}\right),
```

```math
G_k=O(Q_km_k+\sqrt{Q_km_k}P_k),
```

```math
D_{T,k}=O\!\left(m_k+\frac{Q_km_kP_k^2}{B^2}
                  +n(k+1)2^{-k/12}+\log(n+2)\right).
```

Every query fits B if
$`B\ge16c_I((n+4)^3+n+2)`$. This sufficient condition will be used
only in the large-width branch below. Actual inverses have the same
cost and return all their arbitrary dirty work.

The weighted sums have universal constants for $`L'\ge6`$:

```math
\sum_{k=1}^{M}\sqrt{Q_km_k}P_k=O(\sqrt{NL'}),
```

```math
\sum_{k=1}^{M}Q_km_kP_k=O(NL'),\qquad
\sum_{k=1}^{M}Q_km_kP_k^2=O(NL'),\qquad
\sum_{k=1}^{M}Q_km_k=O(NL').
```

For the first sum, factor out $`\sqrt{NL'}`$ and sum
$`\sqrt{1+(k+4)/L'}\,2^{-k/4}`$. For the P-squared sum, it
suffices to sum $`(L'+4+k)2^{-k/2}`$; the other two decay faster.
All bounds hold uniformly for every $`M\le n`$.

For the unweighted terms retain
$`m_k\le L'+4+\lceil\log_2(8n)\rceil`$. The finite sum
$`\sum_{k\ge1}(k+1)2^{-k/12}`$ gives the tail-query ledger

```math
T=O\!\left(\sqrt{NL'}+\frac{NL'}B+
                         M[L'+\log(n+2)]\right),\qquad G=O(NL'),
```

```math
D_T=O\!\left(\frac{NL'}{B^2}+n+
                         M[L'+\log(n+2)]\right).
```

The constant number of queries in each amplified source layer changes
only constants. These replace the original table words exactly on the
full input space; the source's complete rejected action and error
certificate are unchanged. Serial source words and suffix predicates
add $`O(M[L'+\log(n+2)])`$ depth and polynomial local counts in
the regime where this lemma will be composed.

## 4. Prove the stronger low-precision corollary

Assume $`6\le L\le\log_2(n+2)/16`$. First consider

```math
b\le\sqrt N/n.
```

The [original amortized theorem](AMORTIZED_DIRTY_LOOKUP.md#workspace-and-the-simultaneous-resource-ledger)
already supplies the required counts and depth
$`O(NL/b^2+nL+n^2)`$. Here $`NL/b^2\ge Ln^2`$ absorbs both
additive depth terms. Also $`nL=O(\sqrt{NL})`$ uniformly in the
present low-precision range: divide by $`\sqrt{NL}`$ and use
$`L\le\log_2(n+2)/16`$. The function
$`n\sqrt{\log_2(n+2)}/2^{n/2}`$ is uniformly bounded. Thus this
branch already obeys the stronger corollary.

For the remaining case let $`b\gt\sqrt N/n`$. For sufficiently
large n, with an absolute threshold specified by the fixed circuit
constants, use the unary early groups followed by Section 3's tail.
The following estimates make the uniformity in eta explicit.

### A cutoff uniformly smaller than n

Choose q to be the least power of two satisfying

```math
q\ge\frac{4\pi\sqrt{2n}}\eta,
\qquad \delta=\frac\eta{8n},\qquad \rho=\log_2 3.
```

Since $`1/\eta\le2^L`$, and because power-of-two rounding costs
at most a factor two,

```math
q=O(2^L\sqrt n)=O((n+2)^{9/16}).
```

Use the unary local lemma's absolute reservation constant A, increased
to at least one, and its natural cutoff

```math
M=\left\lceil16A(q^\rho+q\log_2q+\log_2q+1)\right\rceil,
\qquad C=16A.
```

This satisfies, uniformly over the entire low-precision range,

```math
M=O\!\left((n+2)^{9\rho/16}
                  +(n+2)^{9/16}\log(n+2)\right),
\qquad \frac{9\rho}{16}\lt0.892\lt1,
```

```math
M[L+\log(n+2)]=o(n).
```

Thus $`M\lt n`$ for all sufficiently large n with one absolute
threshold, not an eta-dependent threshold. At remaining height
$`k>M`$, use

```math
g=\left\lfloor\log_2\frac{k}{Cq}\right\rfloor,
\qquad w=q(2^g-1).
```

The existing unary reservation proof applies unchanged: its conditional
suffix work fits within k divided by eight, leaving the required outer
suffix. Its group height has

```math
g\le\log_2q,\qquad
g\ge\lfloor(\rho-1)\log_2q\rfloor=\Omega(\log(n+2)).
```

All constants here are absolute. The upper bound uses $`q^2\ge n`$;
the lower bound uses the q-to-the-rho term in M. There are consequently
$`O(n/\log(n+2))`$ early groups.

### Charge all early work uniformly

The native phase-state preparation uses $`\ell=\log_2q`$ disjoint
one-qubit words with precision $`\delta/\ell`$. Its word-length
parameter is

```math
1+\log_2((\ell+1)/\delta)=O(L+\log(n+2))=O(\log(n+2)).
```

This follows from $`1/\delta=8n/\eta\le8n2^L`$. Preparation,
actual unpreparation, unary conversion, internal translations, and
selectors therefore have $`O(\log(n+2))`$ T-depth per early group.
Prefetch, unload, and the outer-predicate pair have that same order.
Their total early depth is $`O(n)`$.

The largest local source, convolution, selector, and preparation costs
are a fixed polynomial in n, uniformly over this precision range.
Summing them over at most n groups is still a fixed polynomial. Those
terms fit $`O(\sqrt{NL})`$ T gates and $`O(NL)`$ Clifford gates
with universal constants. This statement charges the actual precision
of every preparation word; it does not hide an eta-dependent native
synthesis cost.

The early query program has $`w\le k/C`$. Its independently returned
parallel dirty pools have the same bound as in the unary proof,

```math
O\!\left(\sqrt N(n+2)^3
                  \sum_{k>M}k2^{-k/2}\right)=o(\sqrt N/n).
```

This estimate is uniform because $`q=\Omega(\sqrt n)`$ and hence
$`M=\Omega(n^{\rho/2})`$, with absolute constants. Their total
T-count is $`O(\sqrt N)`$ and their total Clifford count is
$`O(N)`$. They fit the large-width branch's available dirty pool and
remain disjoint from the logical suffix supplying conditional-zero work.

### Allocate the tail and preserve the literal width threshold

Let $`L'=L+2`$, $`B'_0=L'+n+7`$, and reserve two separate dirty
predicate helpers. The extra query pool has size

```math
B=b-B'_0-2=b-B_0-4\ge b/2.
```

The last inequality follows from the original
$`b\ge17B_0`$. On this branch it gives
$`B\gt\sqrt N/(2n)`$. For sufficiently large n this exceeds
$`16c_I((n+4)^3+n+2)`$, with an absolute threshold. Thus every
unequal-indicator tail query satisfies Section 2's width condition.
The source base, predicate helpers, clean flags, and live query registers
are separately reserved; a returned query pool is reused only after
its completed call or actual inverse.

After early grouping, at most M individual layers remain. Section 3
therefore gives tail query depth

```math
O\!\left(\frac{NL}{b^2}+n+M[L+\log(n+2)]\right)
=O\!\left(\frac{NL}{b^2}+n\right).
```

The sources and predicates add only the displayed sublinear term.
The query counts are
$`O(\sqrt{NL}+NL/b+M[L+\log(n+2)])`$ T gates and
$`O(NL)`$ Cliffords. The unweighted term, and the polynomial local
source/predicate counts, are absorbed into the same budgets. Combining
with the early groups proves all three low-precision bounds on one
circuit.

### Full error and the universal finite-size fallback

Early angle rounding costs at most eta divided by four by the
[complete-frame angular square-sum bound](HOPF_ERROR_ACCUMULATION.md#1-a-sharp-angle-error-bound-on-the-full-frame).
The native source preparation and actual inverse cost at most
$`2\delta`$ for an entire unary group, hence at most eta divided
by four over all early groups. This includes the source's full return,
not only its accepted action. The remaining capped source layers at
accuracy eta divided by four cost at most eta divided by four.
Their queries are exact replacements and add no error. Unitary
telescoping includes all leakage and reference correlations without
resetting any intermediate work. The final error is at most eta.

Every sufficiently-large-n threshold above is absolute under the stated
low-precision restriction. For the finitely many smaller n, use the
existing all-precision hybrid at accuracy eta and width
$`b\ge17B_0`$. Its additive terms
$`nL+n\log(n+2)`$ are $`O(n)`$ on that finite low-precision set,
with one absolute constant. Its nL count term is uniformly absorbed by
$`\sqrt{NL}`$ as above. This supplies the stronger corollary for all
its stated n without an eta-dependent hidden constant or a larger
workspace threshold.

## 5. Complete the uniform theorem and its matching statements

When $`L\gt\log_2(n+2)/16`$, use the already proved
[all-precision hybrid](PARALLEL_DIRTY_LOOKUP.md#the-improved-matching-range).
Its same-circuit resources are

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
```

```math
D_T=O\!\left(\frac{NL}{b^2}+nL+n\log(n+2)\right).
```

Now $`n\log_2(n+2)\lt16nL`$, so its depth is
$`O(NL/b^2+nL)`$ with an absolute constant. On the complementary
range, Section 4's stronger result implies the uniform theorem since
$`L\ge6`$. Selecting these circuits by their declared precision
proves Section 1 for every eta and n.

Under $`b\ge17B_0`$, physical width is $`\Theta(b)`$, and the
inherited complete-frame lower bounds give

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
D_T^\star=\Omega(NL/b^2).
```

For $`b\le\sqrt{N/n}`$, the general depth term nL is at most
$`NL/b^2`$. Also $`b\le\sqrt{NL}`$, so the square-root count
term is at most $`NL/b`$, and $`nL\le NL/b^2\le NL/b`$.
This proves the general matching interval.

In the low-precision range, use its stronger depth and count bounds.
For $`b\le\sqrt{NL/n}`$, n is at most $`NL/b^2`$ and
$`\sqrt{NL}\le NL/b`$. This proves its larger simultaneous
matching interval. Both statements require their intervals to be
nonempty and concern worst-case complete frames. Neither asserts that
the additive depth term is necessary beyond that interval.

## 6. Evidence and contribution boundary

The only new local choice is the count-efficient cap on J and the
complementary definition of H. The controlled bilinear helper echo,
dirty block traversal, four-call indicator echo, and native phase
identities are exactly those of the blocked-query proof. No native
algorithm or correctness contract is altered. The tail source remains
the original full-input source at its proved capped precisions, and
the early source remains the proved unary group with its charged
preparation and actual inverse.

The [bounded allocation checks](../tests/test_rectangular_query_allocation.py)
use integer arithmetic to check power-of-two rounding, asymmetric
address fields, sufficient dirty width, and the two branches of the
block-depth bound. They include odd address lengths, output words longer
than the table, a Clifford single-row query, and a balanced-allocation
family with an exact square-root-of-precision count penalty. The
[existing native query checks](../tests/test_blocked_bilinear_lookup.py)
cover literal phases and arbitrary helper return. Neither finite checks
nor integer resource ledgers replace the analytic uniform sums and
precision split proved above; they do not emit a scalable complete-frame
compiler or establish a new unrestricted depth lower bound.

The [source map](SOURCE_MAP.md) records this allocation and composition
as R45; the [related-work scope](../research/RELATED_WORK.md#21-uniform-precision-and-rectangular-query-allocation-3-october-2026)
identifies its inherited premises. No new external construction or
priority claim is needed.
