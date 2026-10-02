# Grouped program reuse and improved complete-frame T-depth

[Conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md) · [Conditional suffix workspace](CONDITIONAL_SUFFIX_COMPILER.md) · [Parallel dirty queries](PARALLEL_DIRTY_LOOKUP.md) · [Current frontier](OPEN_PROBLEM.md)

A group can load all of its coefficient masks once, use them without
changing the stored program, and erase them with the actual inverse
lookup. The storage is a logical suffix that is zero on the group's
active sector. Two external clean flags suffice. Inactive inputs may
have arbitrary program, source, and computation work: the completed
group acts exactly as identity on that sector.

The local interface below combines with a
[chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md) to give, at fixed
accuracy and sufficiently large square-root-scale dirty width, one
complete real-frame circuit with optimal-order $`T=O(\sqrt N)`$ and
$`D_T=O(n\log\log(n+2))`$. Section 8 prices the entire schedule.
Both ingredients are needed: retaining logarithmically many old late
queries would retain the previous $`O(n\log(n+2))`$ depth allowance.
Depth optimality and the high-precision endpoint remain open.

## 1. Group, program, and local statement

Fix a group of heights $`d,\ldots,d+g-1`$, with $`g\ge1`$.
Its logical registers are a preserved d-bit prefix x, g local bits
$`t_0,\ldots,t_{g-1}`$, and an outer suffix of length
$`r=n-d-g`$. The ideal group acts only when that outer suffix is zero.
The first external zero flag stores $`h=[\mathrm{suffix}=0]`$;
the second is the source branch flag b.

Use a common source width m for the group, with
$`m\ge2`$ and $`g\le m`$. At local height ell, each prefix
$`p\in\{0,1\}^{\ell}`$ has cosine and sine mask words
$`f_{\ell,p,\beta}(x)\in\{0,1\}^m`$, where
$`\beta=0,1`$. Concatenate all these words into

```math
F(x)=\bigl(f_{\ell,p,\beta,j}(x)\bigr)_{\ell,p,\beta,j},
\qquad w=2m(2^g-1).
```

Let $`\mathcal L_F`$ be any exact arbitrary-input XOR lookup

```math
(h,x,W,z)\longmapsto(h,x,W\oplus hF(x),z),
```

returning all its dirty helper bits z. The table has at most
$`2^{d+1}`$ rows and w output bits. Its output W is part of the
conditional suffix, not an extra initialized program register. Lookup
helpers are disjoint from the logical registers and from the live program.

**Local interface.** If $`r\ge16m2^g`$, there is a two-flag
implementation of the group that uses one $`\mathcal L_F`$ and its
actual inverse. Beyond those queries and the outer predicate, its
resources obey

```math
T,G=O(m2^g),\qquad D_T=O\!\left(g\log(m+2)\right).
```

The outer predicate and its inverse use two disjoint returned dirty
helpers, $`O(r)`$ T and Clifford gates, and
$`O(\log(r+1))`$ T-depth using the existing borrowed-MCX primitive.
Thus the explicit full local ledger is

```math
T,G\le2C(\mathcal L_F)+O(r+m2^g),\qquad
D_T\le2D_T(\mathcal L_F)+O\!\left(\log(r+1)+g\log(m+2)\right),
```

where C denotes the respective count. This is not an assumption that
the prefetch itself is shallow or count optimal at every width.

If local layer ell has certified coefficient-vector error at most
$`\zeta_\ell\le1/4`$, the complete group error is at most
$`4\sum_{\ell=0}^{g-1}\zeta_\ell`$. This is a full
initialized-isometry estimate over every logical and dirty input and
every reference. The program and selector work return exactly; source
leakage and the final return of the external flags are included in the
norm estimate.

## 2. A program that can be erased despite source leakage

Compute h first. On h equal to one, all r suffix bits are zero, so
$`\mathcal L_F`$ writes F(x) into W. On h equal to zero its row is
zero and the lookup is exactly identity, regardless of the initial W.

Every subsequent local operation preserves x, h, and each original
program bit. Indeed, the implementation below uses program bits only
as controls of CNOT copies; source preparation, target rotations,
reflections, and selector computation do not act on them. This is a
property of the literal circuit, including its action on leaked source
and flag inputs.

After the internal layers, apply the actual inverse
$`\mathcal L_F^\dagger`$. Since x and h never changed, the program
is erased exactly on the active sector and is unchanged on the inactive
sector. No fresh-source or small-leakage assumption is needed for this
erasure. Finally reverse the actual outer predicate circuit.

## 3. Small-address selectors and the internal enable

At local height ell, define

```math
u_\ell=[t_{\ell+1}\cdots t_{g-1}=0],\qquad
\pi_p=[t_0\cdots t_{\ell-1}=p].
```

The empty suffix and empty prefix predicates equal one. These predicates
depend only on local group bits, not on the d-bit external prefix.
Compute them into separate conditional-zero suffix work before the
amplified stage, and reverse their actual computation afterward.

For each p, copy the ell prefix literals to private zero leaves, negate
the copies for zero literals, and compute their AND in a balanced tree
of zero targets. For ell at least one this uses
$`(2\ell-1)2^\ell`$ wires; ell equal to zero uses one X-prepared
constant-one wire. Across rows and tree levels the Toffolis act on
disjoint private wires. Shared logical bits appear only as CNOT controls.
Compute $`u_\ell`$ similarly with at most $`2g+1`$ wires.
The compute and inverse therefore cost

```math
T,G=O(\ell2^\ell+g),\qquad D_T=O(\log(g+2)).
```

The current stage changes only its target $`t_\ell`$ among the local
logical bits. It preserves all prefix bits and all later local bits,
even on rejected source inputs. Therefore these actual inverses return
their work exactly after the stage. This statement does not assume
that its source core or branch flag returned to zero.

On an outer inactive input these selector work bits need not be clean.
Their computation still has an actual inverse, and the intervening
completed stage is exactly identity when h is zero. The whole
compute-stage-uncompute word consequently remains identity there.

## 4. A literal parallel phase mask

Reserve the geometric source core and its preparation V as in the
[conditional-source proof](CONDITIONAL_GEOMETRIC_SOURCE.md). On active
zero preparation work, V preserves that work for every core input.
It prepares the same m one-hot geometric weights from the zero core.

For each term $`(p,\beta,j)`$ at the current height, privately copy
the five bits

```math
\pi_p,\quad u_\ell,\quad [b=\beta],\quad
W_{\ell,p,\beta,j},\quad {\rm core}_j.
```

Use five zero leaves and four zero AND-tree targets, whose root is q.
For the b-zero literal, negate its private copy. Each term has nine
private work bits. All such trees have depth at most three and are
computed in parallel. Write E for this copy-and-AND word. Apply
$`\mathrm{CZ}(h,q)`$ for every root, then apply the actual inverse
$`E^\dagger`$. The central gates are Clifford. In operator order the
mask is

```math
P_{h,u,W}=E^\dagger
 \left(\prod_{p,\beta,j}\mathrm{CZ}(h,q_{p,\beta,j})\right)E.
```

On zero mask work, this is the exact diagonal phase

```math
(-1)^{\displaystyle h u_\ell
 \sum_{p,\beta,j}\pi_p[b=\beta]
 W_{\ell,p,\beta,j}{\rm core}_j}.
```

Thus on the selected local row it is precisely the programmed Pauli
mask on the source core. The construction holds for every core input,
not only the prepared one-hot state. Every source, program, selector,
and branch bit is preserved by the completed mask; all private mask
work returns exactly.

There are $`q_\ell=2m2^\ell`$ terms. Their native implementation has

```math
T(P_{h,u,W})\le56q_\ell,\qquad
G(P_{h,u,W})=O(q_\ell),\qquad D_T(P_{h,u,W})\le24.
```

These bounds use four exact seven-T, four-layer Toffolis per AND tree
and their actual inverses. All simultaneous trees are private, so no
shared-control T serialization is omitted. Only the original h controls
the central CZs. Consequently h equal to zero makes the whole word
exactly $`E^\dagger E=I`$ on arbitrary mask work, program, and core.
This exact inactive identity does not rely on interpreting dirty copies
as Boolean values.

## 5. Internal normalization-two amplification

For the current height put $`K=XZ=-iY`$ on $`t_\ell`$, and define

```math
Q_{h,u}=H_b\,C_{h=u_\ell=b=1}(K)\,
 V^\dagger P_{h,u,W}V\,H_b,
```

```math
R_{h,u}=I-2[\,h=1,u_\ell=1,b=0,{\rm core}=0^m\,],
\qquad
A_{h,u}=\mathrm{CZ}(h,u_\ell)\,
 Q_{h,u}R_{h,u}Q_{h,u}^\dagger R_{h,u}Q_{h,u}.
```

The conditional K has a fixed number of controls and costs constant
native resources, with literal phase. Its three-controlled Z and X
decompositions borrow and return a fixed number of idle private mask-work
bits, disjoint from the program, selectors, and source core. No separate
helper-free three-controlled X is assumed. The prefactor CZ supplies minus
one only on the enabled sector. Replacing it by Z on h would apply an
incorrect sign when the outer group is active but the current inner
suffix predicate is false.

On h equal to zero, the phase mask is identity on the entire helper
space. The actual V and inverse cancel, as do the outer Hadamards;
the controlled K, reflection, and prefactor are all identity. Hence
$`A_{h,u}=I`$ on the full inactive sector.

On h equal to one, restrict to the invariant subspace with zero source
preparation auxiliaries and zero mask work. The computed selector values
are retained, and the program is arbitrary but preserved. If u is zero,
the mask and the entire stage are identity for every source-core and b
input. This inner inactive statement does require the stated zero-work
subspace; it is not an arbitrary-helper claim on h equal to one.

When h and u are both one, initialization of b and the source core gives
the accepted block $`(cI+sK)/2`$ for the selected program row.
The usual full-isometry amplification proof applies. Source auxiliary
and mask work remain zero for arbitrary core inputs, so the reflection
tests neither that work nor the generally nonzero program and selector
registers.

Each reflection uses the existing logarithmic-depth multi-control
construction with two dirty helpers. Between completed masks, any two
bits of the disjoint private mask-work pool can serve as these helpers;
they return exactly. They are not source, program, selector, or reflection
control bits. The cost is $`O(m)`$ T/Clifford gates and
$`O(\log(m+2))`$ T-depth. The outer predicate separately requires
two helpers outside the entire logical suffix, since all suffix bits
are its controls.

## 6. Sufficient simultaneous suffix reservation

Reuse local selector, enable, and mask computation work between completed
heights, but keep the whole program and source core live throughout the
group. A conservative simultaneous reservation is

| Live register | Suffix bits |
|---|---:|
| Complete program | $`2m(2^g-1)`$ |
| Source core and preparation auxiliaries | $`7m`$ |
| Prefix selectors | $`g2^g`$ |
| Inner suffix enable and computation work | $`2g+1`$ |
| Private mask copies and AND trees | $`9m2^g`$ |

The last row uses the maximum term count at local height $`g-1`$.
The selector bound includes the constant-one root when ell is zero.
For $`1\le g\le m`$, $`m\ge2`$, the sum is at most

```math
11m2^g+g2^g+5m+2g+1
\le12m2^g+7m+1\le16m2^g.
```

These are logical suffix bits, conditionally zero on h equal to one.
They are not also counted as arbitrary dirty lookup banks. The external
query helper pool and the two external predicate helpers are separately
reserved and returned. The only external initialized wires are h and b.

Across all heights, selector and mask count sums satisfy

```math
\sum_{\ell=0}^{g-1}O(\ell2^\ell+g+m2^\ell+ m)
=O\!\left((m+g)2^g\right)=O(m2^g).
```

Each height has a fixed number of V/V-adjoint calls, masks, and source
reflections. Its depth is $`O(\log(m+2)+\log(g+2))`$.
This proves the local resource statement, with load/unload and the
outer predicate still present in their separate terms.

## 7. Full error and the boundary of the result

For a fixed loaded program, replace each actual amplified stage in turn
by its ideal local rotation with source core and b returned. The maximum
row error is at most $`4\zeta_\ell`$. Unitary telescoping bounds
the product by their sum; it does not reinitialize the actual source
between layers or assume that later calls receive unleaked flags.
Selector and program preservation hold on those leaked inputs as well.

The ideal local rotations constitute the prescribed complete g-level
Hopf frame, including every marker column. The inactive outer sector is
exactly identity on arbitrary inputs. Finally, actual program unloading
and the outer predicate inverse preserve the norm bound. The ideal
comparison returns the outer suffix and both flags to their original
states. Source-core and flag leakage of the actual circuit is included
in the same complete-input error estimate, including references.

The local identity alone does not improve a global depth order. A global
use must choose group sizes, common precisions, and actual prefetch
schedules meeting the width and count budgets. The following composition
does this at fixed accuracy, using a different late query. It does not
establish unrestricted source-depth optimality or the high-precision
endpoint.

## 8. Complete-frame theorem at fixed accuracy

**Theorem.** Fix $`0\lt\eta\le1/64`$ and put
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. There is a constant
$`C_\eta`$ such that, for every $`n\ge1`$, $`N=2^n`$, and
$`b\ge C_\eta\sqrt N`$, every prescribed complete real Hopf frame W
has a coherent Clifford+T circuit V using two external clean flags and
at most b arbitrary dirty qubits, with

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T=O_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(n\log\log(n+2)\right).
```

All bounds hold on the same circuit. The error includes returned work
and arbitrary references. No measurements, resets, supplied catalysts,
QRAM, or uncharged quantum oracles are used. T-depth permits arbitrary
Clifford circuits between its layers; their gate count is included in G,
but their elementary depth is not claimed small.

The existing [worst-case T-count lower bound](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
is $`\Omega_\eta(\sqrt N)`$, so count has the optimal order.
The available general depth lower bound at this width remains only
$`\Omega(1)`$. Constants here may depend on the fixed accuracy. This
does not strengthen the all-precision, all-width theorem or its literal
$`17(L+n+7)`$ reservation.

### Group sizes and unchanged accuracy

Use the original unfiltered additive cap. Set

```math
h_n=\lceil\log_2(8n)\rceil,\qquad
m=L+4+h_n,\qquad R=\min\{n,256m\}.
```

At a group start d let $`k=n-d`$. While $`k\gt R`$, take

```math
g=\left\lfloor\log_2\frac{k}{64m}\right\rfloor,
\qquad d\longleftarrow d+g.
```

Here $`2\le g\le m`$. The remaining outer suffix has length
$`r=k-g`$, and

```math
16m2^g\le k/4\le k-g,\qquad
w=2m(2^g-1)\le k/32.
```

Thus every group meets the literal reservation in Section 6. Even the
last such group ends with remaining length greater than $`255m`$,
so every layer it contains has the original capped precision m.
After grouping stops, the remaining $`R'\le R`$ individual layers
use their original widths

```math
m_k=L+4+\min\{k,h_n\},\qquad 1\le k\le R'.
```

Consequently the [original additive certificate](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
applies without changing any layer's precision or error allocation.
Group unloading and the actual predicate inverses preserve the full
error estimate from Section 7. Unitary telescoping over groups and late
layers gives the displayed complete-frame bound.

### Price the early prefetches

For each of the w program bits, use its own balanced two-indicator
bilinear query on the $`d+1`$ address bits $`(h,x)`$.
Use the existing counter indicator, with a conservative cubic
count/width majorant, and allocate separate arbitrary dirty indicator
and helper pools to different output bits. Each single-bit bilinear
middle word has constant T-depth. The shared public address enters the
counter circuits only through CNOT controls, so all w queries can run
in parallel. The exact [query construction](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction)
and its actual inverse give, with $`Q=2^{d+1}`$,

```math
T_{\rm load},w_{\rm dirty}
=O\!\left(w\sqrt Q\,(n+2)^3\right),
\qquad D_{T,\rm load}=O(\log(n+2)),
```

```math
G_{\rm load}
=O\!\left(wQ+w\sqrt Q\,(n+2)^3\right).
```

In particular, the arbitrary binary coordinate changes cost
$`O(wQ)`$ Clifford gates; this term is not omitted. The w output
bits are the separately counted logical program, not additional clean
external work. All private dirty query work returns before group use.

Group starts have distinct k, and $`w\le k/32`$. When the early set
is nonempty, $`R=256m`$, so geometric summation gives

```math
\sum_{\rm early}T_{\rm load}
\le O\!\left(\sqrt N(n+2)^3
       \sum_{k>R}k2^{-k/2}\right)=O_\eta(\sqrt N),
```

```math
\sum_{\rm early}G_{\rm load}
\le O\!\left(N\sum_{k>R}k2^{-k}\right)
       +O_\eta(\sqrt N)=O_\eta(N).
```

The same first estimate bounds each peak private dirty pool. Load and
unload multiply these costs only by two. Local masks, sources, and
outer predicates cost $`O(n^2)`$ gates altogether: each group has
$`m2^g=O(k)`$, predicate length at most k, and there are at most n
groups. This polynomial fits both stated exponential count budgets.

For fixed L the number of early groups is $`O_\eta(n/\log(n+2))`$.
For sufficiently large n, when $`k\ge\sqrt n`$ the chosen g is at
least a fixed positive multiple of $`\log n`$. Below that threshold
there are at most $`\sqrt n`$ groups. The finitely many smaller n
are covered by accuracy-dependent constants. Therefore all early
prefetches and outer predicates together have depth $`O_\eta(n)`$.
Since $`\sum g\le n`$, their internal layers have total depth

```math
O\!\left(n\log(m+2)\right)
=O_\eta\!\left(n\log\log(n+2)\right).
```

### Use a different query in the late segment

The late segment must not revert to an $`O(n)`$-depth router at
every layer. At remaining height k its existing mask query has
$`Q_k=4N2^{-k}`$ rows, output width $`m_k\le L+4+k`$, and
$`n-k+2`$ address bits. Split that address as evenly as possible.
On each nonempty half of length s use the
[chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md) with

```math
\ell_k=\min\{s,2^{\lfloor k/12\rfloor}\}.
```

The zero-bit half is a Clifford X. The indicator's cubic workspace
overhead is at most a constant times $`2^{k/4}`$. Its balanced
bilinear query therefore satisfies

```math
T_k,w_k
\le O\!\left(\sqrt N\,2^{-k/2}
       [L+k+4+2^{k/4}]\right),
```

```math
G_k\le O\!\left(N2^{-k}(L+k+4)\right)+O(T_k),
```

```math
D_{T,k}\le O\!\left(
L+k+(k+1)n2^{-k/12}+\log(n+4)\right).
```

The weighted count sums converge, giving late T-count
$`O_\eta(\sqrt N)`$, Clifford count $`O_\eta(N)`$, and peak dirty
width $`O_\eta(\sqrt N)`$. Since $`R'=O_\eta(\log(n+2))`$,
the depth sum is

```math
\sum_{k=1}^{R'}D_{T,k}
=O_\eta\!\left(n+\log^2(n+2)\right).
```

The original late source calls, flag reflections, and separately priced
suffix predicates fit an additional $`O_\eta(\log^2(n+2))`$ depth.
Their counts are polynomial. The constant number of mask calls in
each amplified stage does not change these orders.

Finally keep the original $`B_0=L+n+7`$ dirty base, reserve two
disjoint predicate helpers when needed, and reuse the returned query
pools between completed calls. All early and late peaks fit
$`C_\eta\sqrt N`$ for a sufficiently large fixed constant. No
conditional suffix storage is counted again as dirty query workspace.
Combining the early and late ledgers proves the theorem.

## 9. Bounded executable evidence

The [group checks](../tests/test_grouped_program_prefetch.py) cover
literal mask circuits, a native AND fixture, native geometric preparation,
inactive arbitrary dirty scratch, and exact program erasure through
retained source leakage. Negative controls protect the actual prefetch
inverse, mask uncomputation, original-h guard, and internal-enable phase.
Exact integer checks exercise the reservations and cap agreement of
Section 8. Group reflections and mask actions also use reduced diagonal
operators; these tests do not emit the complete native frame compiler.
The [chunked-indicator checks](../tests/test_chunked_dirty_indicator.py)
separately test the late-query interface and resource sums. The
asymptotic theorem follows from the analytic ledgers above.
