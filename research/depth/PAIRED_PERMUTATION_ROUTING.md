# Short-label permutation wrappers with paired dirty arrays

[Native phase grid](PARALLEL_BILINEAR_WORD_QUERY.md#2-parallel-bilinear-outputs-and-their-controlled-version) ·
[Retained selected routing](SELECTED_DIRTY_ROUTING.md) ·
[Conditional frame transfer](QUERY_DEPTH_FRAME_TRANSFER.md)

An arbitrary selected permutation can be applied to one physical array,
with its inverse applied to a separate dirty array, using a sparse native
phase grid. The second array returns in a complete layout/body/unlayout
wrapper. This removes the wire-count contribution to that wrapper's
T-depth by paying a larger polynomial dirty allocation.

## 1. Contract and distinction from a returned router

Let $`p,w\ge1`$, let $`k=\lceil\log_2\max\{p,2\}\rceil`$, and
let a preserved k-bit register R select one of p prescribed permutations
of w physical qubit sites. Padded selector values specify the identity.
Let H be a separate preserved activity bit. Write $`P_s`$ for the
binary permutation matrix at $`s=(R,H)`$, with $`P_s=I`$ when H is
zero. The action convention is $`A\mapsto P_s A`$ on the column of
computational bits held by the physical array A.

There is an exact native Clifford+T word on two disjoint arbitrary
w-site arrays A,E implementing

```math
\mathcal R:\quad
|s,A,E,d\rangle\longmapsto
|s,P_s A,P_s^{-1}E,d\rangle                           \qquad(1)
```

with

```math
T,G,w_{\rm extra}=O(pw),\qquad
D_T=O(\log_2(k+2)).                                  \qquad(2)
```

The additional dirty width in (2) excludes A and the selector but
includes E, the private grid, the dirty indicator, and all its helpers.
There is no initialized routing work. Every additional helper except
E returns at the boundary (1); **E is permuted there**. The literal
phase is one. The equality is on the full input space, hence includes
arbitrary dirty states and external references.

Let a body B preserve s exactly and act as identity on E. It may act
on A and on any other disjoint physical registers. Then, with rightmost
factors acting first,

```math
\mathcal R^\dagger B\mathcal R
=\sum_s|s\rangle\langle s|\otimes
    \bigl(P_s^{-1}B_sP_s\bigr)\otimes I_E.             \qquad(3)
```

Thus the complete wrapper returns E exactly. The actual body must
preserve s on its entire physical input space, including occupied work
and leakage sectors. Preserving it only on an ideal initialized image
is insufficient for the exact cancellation in (3). Internal native gates
of B may change s: only the complete full-space operator must preserve
it, and no other operation is interleaved with that body. All source
anchors and tested work in B exclude E.

This is a wrapper theorem, not a standalone returned-work implementation
of $`P_s(A)`$. R91's paired XOR translations are the special case
$`P_s=P_s^{-1}`$ with an additional XOR group law. No group law among
the prescribed permutations is needed here. Conversely, (1) does not
provide R91's removal of an unknown dirty program: s is an available,
preserved selector.

## 2. A native sparse selected shear

There are $`p'=2^{k+1}=O(p)`$ possible s values, including activity.
Reserve a p'-bit arbitrary dirty indicator Y. For every pair
$`(u,j)`$, with $`0\le u\lt p'`$ and $`0\le j\lt w`$,
reserve one private arbitrary dirty grid bit. Apply the R91 phase-grid
word to the bilinear map

```math
\mathcal B:\quad
E_j\longmapsto E_j\mathbin\oplus
       \bigoplus_{u=0}^{p'-1}Y_u(P_u A)_j.             \qquad(4)
```

For each u,j its row form contains exactly one A bit. The native
linear word therefore has $`p'w`$ CNOTs, rather than the generic
dense-matrix allowance $`p'w^2`$. R91 gives at most $`8p'w`$
T gates, eight physical T layers, $`O(p'w)`$ elementary Clifford
gates, and $`p'w`$ returned private grid bits. Shared controls receive
only Clifford interactions; every T gate acts on a private grid wire.
Neither A nor Y changes in the completed map (4).

Let $`I_s`$ be the retained exact dirty indicator
$`Y\mapsto Y\oplus e_s`$. Execute chronologically

```math
\mathcal B,\quad I_s,\quad
\mathcal B^\dagger,\quad I_s^\dagger.                  \qquad(5)
```

The two output contributions differ by precisely $`P_s A`$. Hence
(5) is the full-input shear

```math
U_P:(A,E)\longmapsto(A,E\oplus P_sA),                 \qquad(6)
```

and returns Y and every private helper. The indicator's counts and
width are $`O(p')`$, and its depth is
$`O(\log_2(k+2))`$. Its controls exclude A,E and all live grid
wires. Complete sequential calls can reuse the returned grid/indicator
pool. Replacing every $`P_u`$ by $`P_u^{-1}`$ and exchanging the
array roles gives

```math
L_{P^{-1}}:(A,E)\longmapsto
                  (A\oplus P_s^{-1}E,E)              \qquad(7)
```

with exactly the same resource orders. Every dagger in (5) is the
actual reversed native word. In particular, the relative phases of
the underlying four-T words are canceled by R91's exact dirty-wire
echo; no relative-phase Toffoli is substituted for (4).

## 3. The paired word and complete work return

Execute chronologically $`U_P,L_{P^{-1}},U_P`$, then the fixed
swap of A and E. The array values are

```math
\begin{aligned}
(A,E)&\longmapsto(A,E+P_sA)\\
&\longmapsto(P_s^{-1}E,E+P_sA)\\
&\longmapsto(P_s^{-1}E,P_sA)\\
&\longmapsto(P_sA,P_s^{-1}E).
\end{aligned}                                        \qquad(8)
```

Here addition is over $`\mathbb F_2`$. The final fixed swap is
an explicit $`3w`$-CNOT word. The three shear pools are reused only
after their complete return. At the peak the additional systems are
E, Y, its helper pool and one $`p'w`$ grid, so (2) is simultaneous.
The inverse uses the actual reverse of this complete word. If H is
zero, (8) is exactly identity on both arrays even though its individual
shears need not be identity.

Since B is block diagonal in s and does not touch E, the E factors
in (3) commute through B and cancel. This proves exact E return even
when E is entangled with s,A, other work or an external reference.
Every helper already returned at each shear boundary. No clean-state
promise has been used anywhere in this identity.

An initialized-isometry error bound for the original routed body is
unchanged: conjugation by the exact layout words preserves operator
norm, and the paired auxiliary factors cancel identically. The entire
wrapper is identity on H zero whenever the original body has that
property, including when it obtains it by its own actual-inverse
sandwich. Pairing does not impose an inactive identity on an unrelated
body.

## 4. Fair baseline and applications

The retained arbitrary selected Beneš router has depth
$`O(\log(w+2)\log(k+3))`$. Its contiguous-interval specialization
already simplifies this to $`O(\log(w+2))`$; that sharper baseline
is the relevant comparison for the current frame construction.
Equation (2) removes its wire-count factor inside a complete wrapper,
at the cost of $`O(pw)`$ additional dirty volume. For
$`p=\mathrm{poly}(w)`$, the new depth is $`O(\log\log(w+4))`$.
The statement does not claim that this larger allocation fits every
literal minimum frame width.

The global early stage uses a constant number of selected layout
positions. Its bin selectors and the local interval metadata have
$`O(\log(m+2))`$ bits, there are polynomially many selectable
layouts, and their physical arrays have polynomial size in m. Each
selected body preserves its bin/metadata selector exactly: the relevant
queries read those fields, target disjoint banks, and retain them until
the actual inverse layout. Parallel block layouts have disjoint arrays
and private pools. Thus their paired E arrays and grids add only
polynomial count and dirty width, already allowed by that stage's
standalone budget, and have depth $`O(\log\log(m+4))`$.

The two uses keep their actual selector interfaces:

| Layout position | Selector excluded from the routed array | Body between the two layout boundaries |
|---|---|---|
| Selected fine-erasure episode | H, the coarse bin R, and its loaded interval starts/lengths | Fine source, tag XOR and actual inverse source; the old fine tag may change, but R and interval parameters do not |
| Each compact-compressor block | H and that block's fixed split/interval metadata | Local shear and global High erasure; both read the metadata and target separate Low/High banks |

Keep the original parameter-load words and their actual unqueries.
For the new permutation table, enumerate the old completed router's
action on **every** selector/parameter bit pattern, including malformed
active metadata and occupied values. Each such action is a physical site
permutation. Use that exact complete table; do not replace an old invalid
row by the identity unless it already was identity. The number of rows
is still polynomial in m. This preserves the original full-space
extension, not just its action on ideal loaded parameters. H-zero rows
are identity by the old exact guard. Every listed selector field is
disjoint from both A and E and remains unchanged by its actual body.

Every selected permutation fixes its prescribed padding at the completed
boundary. Any padding and the new E arrays are arbitrary, and remain
excluded from the body's initialization predicates. This replacement
does not modify a source, remove an inverse, or change the decoded
operator. Its relevance to the remaining long-query contribution is
proved separately in the conditional transfer chapter.

## 5. Selected address packing for completed XOR queries

Suppose a native exact query compiler implements, for every finite
q-bit table f, the full-space operator

```math
Q_f|x,Z,w\rangle=|x,Z\mathbin\oplus f(x),w\rangle,     \qquad(9)
```

with literal phase one and arbitrary output, helper and reference
states. The address x may change inside its native word. This section
uses only its exact return in (9), together with the paid query bounds;
it does not run multiple such words concurrently on shared address wires.

Let s be a k-bit selector and P a common physical address register. For
each s, specify an ordered subset $`A_s`$ of a local physical array A,
with $`a_s=|A_s|\le a_{\max}`$. The array contains the union of the
local address subsets and any required padding. It is disjoint from
P, s, the entire final output word Z, and every private query/router
register. The desired selected query is

```math
|s,P,A,Z\rangle\longmapsto
|s,P,A,Z\mathbin\oplus f_s(P,A_s)\rangle.             \qquad(10)
```

Every row $`f_s`$ is specified, including malformed and inactive
selector patterns. When activity H is part of s, an H-zero identity is
implemented by a zero output row. Selector-dependent output positions
can be embedded in one fixed output word by putting zeros in the other
slices. Its full physical output union must be disjoint from the entire
address union, not merely the address selected on an ideal input.

Choose a physical permutation $`\Pi_s`$ of A that places the ordered
$`A_s`$ in the first $`a_s`$ positions of a fixed
$`a_{\max}`$-site block. Complete the permutation arbitrarily on
the other sites. On a prescribed zero row no packing is needed, so the
identity permutation is allowed. The subsets depend only on the short
selector; P affects table values, not the physical permutation.

Apply the paired router to A and a separate arbitrary array E. The
sparse-grid proof of Sections 2–3 also applies to an arbitrary total
k-bit permutation table, with $`2^k`$ rows and no separate activity
bit. Run one complete query B on P, s and the fixed packed block, with
row s reading only its first $`a_s`$ positions and ignoring all
filler values. Finally apply the actual inverse router. At the complete
query boundary B preserves P and s and is identity on E and the unused
A positions. Equation (3) therefore gives

```math
\mathcal R^{\dagger}B\mathcal R
=\sum_s|s\rangle\langle s|\otimes
       \bigl(\Pi_s^{-1}B_s\Pi_s\bigr)\otimes I_E.    \qquad(11)
```

This is exactly (10), with every address, E and private helper returned.
It remains true when B temporarily changes P or s, applies non-Clifford
gates to them, or entangles them with its work. The proof is a literal
operator equality and includes arbitrary dirty/reference entanglement.
It requires exact boundary preservation on every physical sector.

There is one complete long query on at most

```math
a_{\mathrm{query}}=|P|+a_{\max}+k                         \qquad(12)
```

address bits. The union size $`w=|A|`$ enters only the router ledger:
both router appearances together use $`O(2^k w)`$ T gates,
elementary Clifford gates and additional dirty bits, and have T-depth
$`O(\log_2(k+2))`$. Their additional width includes E. The query
receives its separately charged pool, disjoint from P,s,A,Z,E and all
other live registers. Returned router-only helpers may be reused at
complete boundaries, but E remains reserved across B. Keeping the
pools disjoint gives the same additive simultaneous bound. Actual
inverse queries and routers have the same resource orders.

A Boolean phase uses the exact query/Z/actual-unquery/Z dirty echo on
a separate arbitrary bit, so only a constant number of these complete
queries is required. No controlled version of an arbitrary compiler is
assumed: activity is an address bit with specified zero rows. Existing
malformed-selector rows must be copied from the intended total function;
they cannot be changed merely because the ideal source never reaches
them.

For an unselected query position whose genuine input union X is small,
concatenate all its output bits into one table on X. This realizes the
same product of completed XOR queries, provided the full output union
is disjoint from X. Bounds on individual addresses alone do not bound
X. Selected packing supplies the alternative when ragged physical
subsets have a large union. The conditional R98 proof audits these two
cases separately and charges their common output word and live arrays.

The [six bounded fixtures](../../tests/test_paired_permutation_routing.py)
check the native sparse-grid phase echo on all 32,768 occupied inputs,
all three-site permutations including noninvolutive cycles, actual
inverse work return with every selector and activity value, and a
negative case in which a permanent selector change prevents auxiliary
return. The two packing fixtures check a native body that temporarily
changes its address and selector but restores them, and a ragged
maximum-width-two/union-width-three example with common prefix, explicit
inactive and invalid rows, arbitrary output offsets and paired arrays.
The native evaluator is shared with R91's retained tests. These exact
finite checks support the identities; the resource bounds above are
analytic and do not assert a scalable emitter.
