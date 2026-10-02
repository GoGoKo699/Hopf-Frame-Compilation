# A dirty-workspace upper bound for T-depth

[Two-clean grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md) · [Operator source](OPERATOR_SOURCE_COMPILER.md) · [QBP error contract](QBP_APPROXIMATION.md)

T-count measures how many non-Clifford gates are used. T-depth measures
how many sequential layers of those gates remain when arbitrary Clifford
circuits may be placed between them. The two resources have different
workspace tradeoffs. The construction below chooses larger dirty lookup
banks to reduce T-depth, and can use more T gates than the count-optimized
compiler.

The [parallel-indicator refinement](PARALLEL_DIRTY_LOOKUP.md) retains
count-efficient banks and uses extra dirty selectors instead. Under its
stronger sufficient-width condition, it obtains the best retained T-count
and low T-depth simultaneously. This chapter supplies the routing primitive
and the wider-range bank schedule used in that refinement.

The [amortized extension](AMORTIZED_DIRTY_LOOKUP.md) reuses a smaller indicator
pool and shares dirty traversal across the chunk addresses. For every
$`L\ge6`$, two clean qubits and $`b\ge17(L+n+7)`$ give one complete
real-frame circuit with $`T=O(\sqrt{NL}+NL/b+nL)`$, $`G=O(NL)`$,
and $`D_T=O(NL/b^2+nL+n^2)`$. Count and depth are both optimal in
order when $`b\le\sqrt{NL/(nL+n^2)}`$ above that threshold.
At fixed L the earlier matching depth range $`b\le\sqrt N/n`$
is retained. The general-precision schedules below remain available
where their source or workspace costs are sharper.
Its [capped source allocation](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
reduces source depth to $`O(nL+n\log(n+1))`$; the remaining per-layer
routing keeps the fixed-accuracy total at $`O(N/b^2+n^2)`$.

**Theorem.** Let $`n\geq1`$, $`N=2^n`$,
$`0\lt\eta\leq1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. Put

```math
B_0=L+n+7,\qquad
\ell_*(n)=1+\log_2^*(n+2).
```

For every $`a\geq2`$ and $`b\geq2B_0`$, the prescribed complete real
Hopf frame has a coherent Clifford+T compiler with

```math
\|VJ_a-J_a(W\otimes I_b)\|\leq\eta,
```

```math
D_T=O\!\left(\frac{NL}{b}+L\ell_*(n)+n^3\right),
\qquad T=O(NL),\qquad G=O(NL).
```

Only two external clean qubits are used; unused clean and dirty qubits
remain untouched. The same full-input norm includes initialized leakage,
arbitrary dirty inputs, and their reference correlations. All inverses
are actual circuit inverses. No measurements, resets, supplied magic
states, catalysts, or uncharged quantum oracles are introduced.

This is an upper bound, not an optimal T-depth theorem. It also is not a
claim about total elementary-gate depth: Clifford circuits between the T
layers still have nonzero depth. The high-precision allocation
$`a=2,b=N+n+7,L=N`$ does not satisfy the larger-bank hypothesis here.

## 1. Model and the routing primitive

A T layer applies T or T-dagger on distinct physical qubits. An arbitrary
Clifford circuit may occur before, after, or between these layers. We
assume the same all-to-all logical connectivity as the other compiler
chapters. Clifford gates are counted in G, even though they do not
contribute to $`D_T`$.

The lookup-depth mechanism is established by Low, Kliuchnikov, and
Schaeffer, *Trading T gates for dirty qubits in state preparation and
unitary synthesis*, [arXiv:1812.00954v2](https://arxiv.org/html/1812.00954v2),
Table 2 and Appendix B. Their dirty SelectSwap construction has T-depth
$`O(Q/\mu+\log\mu)`$. We give an explicit routing schedule and then
apply it to the complete-frame compiler's dirty-selector queries.
The primitive is not a new lookup result.

### A shared-control Fredkin batch has at most four T layers

Consider R controlled swaps, all controlled by c, with disjoint target
pairs $`(u_i,v_i)`$. Each Fredkin is two CNOTs surrounding a Toffoli.
Hadamards on the Toffoli targets reduce the middle batch to
$`\prod_i\mathrm{CCZ}(c,u_i,v_i)`$. These changes are Clifford and
require no ancillary qubits.

For bits, the exact phase-polynomial identity is

```math
4cuv=c+u+v-(c\oplus u)-(c\oplus v)
-(u\oplus v)+(c\oplus u\oplus v)\pmod8.
```

Summing over the disjoint pairs produces $`Rc`$ and the remaining local
terms. Each parity can be placed on one of its existing target wires by
CNOTs and then uncomputed. The phases can be scheduled in four layers:

| T layer | Parities receiving a T or T-dagger phase |
|---|---|
| 1 | $`u_i,v_i`$ for every i, together with $`T^R`$ on c |
| 2 | $`c\oplus u_i,c\oplus v_i`$, with negative signs |
| 3 | $`u_i\oplus v_i`$, with negative signs |
| 4 | $`c\oplus u_i\oplus v_i`$ |

Within each row, the phase targets are distinct. The common c wire is
only a control of the parity-computing CNOTs in rows 2 and 4; its value
is never copied into a clean ancilla. The $`T^R`$ in row 1 contributes
at most one T gate after removing its Clifford power. The literal batch
therefore has

```math
D_T\leq4,\qquad T=6R+(R\bmod2),\qquad G=O(R).
```

The phase identity is exact, including the common phase. Applying the
surrounding Clifford changes gives the same bounds for the Fredkins.
The parity CNOT circuits can have substantial elementary depth; this
argument does not count them as constant-depth circuits.

### A bank router

Let $`\mu`$ be a power of two. A binary router moves one addressed
m-bit bank into position zero. In its first level, it conditionally
swaps neighboring pairs of banks; the next level conditionally swaps
the selected representatives of neighboring groups; subsequent levels
continue to the root. Every level uses one low-address bit as its
common control and disjoint target pairs.

There are $`\log_2\mu`$ levels and $`\mu-1`$ word swaps in total.
The preceding schedule gives

```math
D_T(\mathcal R)\leq4\log_2\mu,
\qquad T(\mathcal R),G(\mathcal R)=O(\mu m).
```

This is an exact permutation of arbitrary bank states. Its actual inverse
has the same resources. No clean copies of the address are required.

## 2. Exact lookup with depth-optimized banks

Consider a Q-row, m-bit XOR table. Reserve $`\mu m`$ dirty bank bits,
separate from its m-bit output and dirty selectors, where
$`1\leq\mu\leq Q`$ is a power of two. The high-address loader
$`\mathcal L`$ writes the $`\mu`$ words for each high address into
these banks. The existing two-pass dirty traversal uses
$`O(Q/\mu)`$ Toffolis, and hence $`O(Q/\mu)`$ T-depth by any fixed
exact Toffoli decomposition. Its many leaf CNOTs are Clifford.

Let C copy the routed bank into the arbitrary output word by CNOTs.
Use the chronological sequence

```math
\mathcal L,\ \mathcal R,\ C,\ \mathcal R^\dagger,\
\mathcal L^\dagger,\ \mathcal R,\ C,\ \mathcal R^\dagger.
```

If the initially selected bank contains z, the two output contributions
are $`z\oplus f(x)`$ and z. Every bank and selector returns exactly.
The basis identity extends to arbitrary entangled inputs and references.
The loader inverse is applied after reversing routing, and its targets
exclude the output word. The resources are

```math
D_{T,\mathrm{query}}=O(Q/\mu+\log\mu),
\qquad T_{\mathrm{query}}=O(Q/\mu+\mu m),
\qquad G_{\mathrm{query}}=O(Qm).
```

If a table is inactive for some h or mode address, assign that row the
zero output word. The completed query is exactly identity on that sector.
Its internal routing need not be identity there. Thus the construction
does not add a control to every routing or leaf gate.

Reserve $`B_0`$ dirty wires for the existing source, selectors, and
returned helpers. The remaining pool has

```math
K=b-B_0\geq b/2\geq B_0.
```

Every source width m in the grouped or reserved-tail construction is at
most $`B_0`$. Choose the largest power of two satisfying

```math
\mu\leq\min\{Q,K/m\}.
```

This fits the disjoint banks and differs from the T-count-optimized
choice, which may stop near $`\sqrt{Q/m}`$. Since rounding loses at
most a factor of two,

```math
D_{T,\mathrm{query}}
=O\!\left(1+\frac{Qm}{b}+\log Q\right),
\qquad T_{\mathrm{query}},G_{\mathrm{query}}=O(Qm).
```

The same reasoning applies to the constant-size streamed coarse symbols.
Their bank pool can be reused between queries.

## 3. Composition with the grouped full-frame compiler

### Returned helpers for logarithmic-depth predicates

Khattar and Gidney, *Rise of conditionally clean ancillae for optimizing
quantum circuits*, [arXiv:2407.17966v1, Section 5.4](https://arxiv.org/html/2407.17966v1#S5.SS4),
give an exact k-controlled X with two arbitrary dirty helpers,
$`O(k)`$ Toffolis, and $`O(\log(k+1))`$ Toffoli depth.
Use literal exact Clifford+T decompositions of its Toffolis, including
their actual inverses. This gives the same asymptotic T-count, Clifford
count, and T-depth without measurements or relative-phase substitutions.
Negative controls add only X gates; an exact multiple-controlled Z is
obtained by conjugating one target with Hadamards.

The two helpers come from the extra dirty pool $`b-B_0\geq B_0`$,
not from either initialized compiler flag. Predicates, reflections, and
selected local-word operations occur between completed queries. At
these times the word banks and any indicator registers have returned
exactly and are idle. Two distinct pool wires are therefore available,
disjoint from all logical, label, flag, and suffix controls and from the
predicate target. Every such predicate returns both helpers before the
next query, source, or local-word operation. No pool wire is borrowed
while it holds a live query value. The equality is on arbitrary helper
inputs, so previous approximation leakage and reference correlations
do not invalidate this reuse.

### Group schedule

Use exactly the partition, coefficient tables, coarse circuits, and
precision choices of the [conditional-suffix construction](CONDITIONAL_SUFFIX_COMPILER.md).
For a group of height s and ending depth e, set $`r=n-e`$. Its
coefficient table has

```math
Q_g=O(s2^e),\qquad
m_g=L+\lfloor r/4\rfloor+8,
\qquad s\leq2^{r/C_1},\qquad C_1>2.
```

There are $`R=O(\ell_*(n))`$ groups. The r values are distinct, and
the group heights partition at most n logical depths. The construction
uses only a constant number of scalar-source calls per group, including
both orientations and the actual inverses in amplification. Schedule
each existing source word serially; it has T-depth $`O(m_g)`$.
No parallel source construction is assumed.

The coefficient queries contribute

```math
O\!\left(\frac{Q_gm_g}{b}+n+1\right)
```

per group. Their logarithmic routing term is $`O(n+1)`$, because
$`Q_g=O(s2^e)`$ and $`s,e\leq n`$.

The coarse interpreter streams $`O(s)`$ symbols at each of its s local
depths. At absolute depth d the table has $`O(2^d)`$ rows and constant
output width. Summing its query-depth terms within one group gives

```math
O\!\left(\frac{s2^e}{b}+s^2(n+1)\right).
```

For the selected atoms, unroll the $`O(s)`$ possible column types.
A fixed type is selected by equality tests on the unchanged term label,
together with its mode and direction. Its marker condition is either an
all-zero local word or one specified local bit followed by a zero suffix.
The corresponding atom-flag toggle needs a constant number of conjunction
tests on $`O(n)`$ wires. For a negated marker, toggle on the selected
type and then on the selected type together with the positive marker.
The new predicate circuit gives $`O(\log(n+1))`$ T-depth for these
tests. Decode a type into an existing private temporary bit, apply its
at most s controlled Hadamards and one controlled X serially, and erase
the type bit. The decoding address is unchanged. Retain the existing
h/mode controls on every consumed operation, so the completed word is
exactly identity on inactive sectors. Each native Hadamard has only a
fixed number of these controls and therefore constant native cost by
the conditional-suffix construction's bounded-control lemma.
This gives $`O(s^2+s\log(n+1))`$ atom T-depth per group, including
reverse atoms through their actual inverses. No additional clean bit or
parallel controlled-Clifford theorem is required.

The streamed coarse alphabet has constant size. Its interpretation has
bounded controls, and its selected suffix tests use the same two returned
helpers. Even charging one fresh conjunction test and its inverse for
each of the $`O(s^2)`$ symbols costs only
$`O(s^2\log(n+1))`$ T-depth. Suffix tests and amplification reflections
outside these streams cost $`O(\log(n+1))`$ each. Since

```math
\sum_gs_g\leq n,\qquad
\sum_gs_g^2\leq n^2,\qquad R\leq n,
```

the atom work is $`O(n^2)`$ in total, and the streamed predicate work
is $`O(n^2\log(n+1))`$. Both fit $`O(n^3)`$. The
$`O(s^2(n+1))`$ coarse-query routing allowance also sums to
$`O(n^3)`$, and is retained explicitly rather than claimed optimal.

The weighted table sum proved in the grouped compiler is

```math
\sum_g Q_gm_g
\leq O(N)\sum_{r\geq r_0}
2^{-(1-1/C_1)r}(L+r+8)=O(NL).
```

Also $`\sum_gs_g2^{e_g}=O(N)`$ and
$`\sum_gm_g=O(L\ell_*(n)+nR)`$. Thus the grouped part has

```math
D_T=O\!\left(\frac{NL}{b}+L\ell_*(n)+n^3\right).
```

The fixed number of deepest reserved layers uses the old operator-source
blocks with the same depth-optimized banks. Their contribution is
$`O(NL/b+L+n)`$. Each such layer has a constant number of queries
and predicates, source width $`O(L+1)`$, and logarithmic routing at
most $`O(n)`$. When n is below the fixed grouping threshold,
using that original compiler for every layer has the asserted form,
with constants depending only on the fixed threshold.

Changing an exact lookup implementation does not change its unitary,
accepted coefficient, source precision, or inactive-sector behavior.
Consequently the original amplification proof, error budget, and
complete-input telescoping remain valid without alteration. In
particular, the approximate returned-work guarantee is not being
replaced by an assumption that intermediate work is reset.

The depth-selected query T-count is at most $`O(Q_gm_g)`$, and its
Clifford count has the same bound. Their sums are $`O(NL)`$.
Coarse tables and the fixed-degree scheduling overheads also fit this
bound, since every fixed polynomial in n is $`O(2^n)`$ and $`L\geq6`$.
This proves the stated simultaneous $`T,G=O(NL)`$ guarantees.
It does not retain the sharper count-optimized
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ bound on this same schedule.
That limitation of maximal banking is removed, under a stronger workspace
condition, by the separate [parallel-loader construction](PARALLEL_DIRTY_LOOKUP.md).

For comparison, applying only the depth-optimized queries to the older
layer-by-layer compiler gives the explicit upper schedule

```math
D_T=O\!\left(\frac{NL}{b}+nL+n^2\right),
\qquad T,G=O(NL),\qquad b\geq2B_0.
```

Here the source is paid at all n depths: the widths
$`m_d=L+n-d+4`$ sum to $`O(nL+n^2)`$. The constant number of
queries per depth contributes $`O(NL/b+n^2)`$, and the constant
number of suffix tests contributes $`O(n\log(n+1))`$ using the
same idle-pool helper allocation. The two initialized-flag reflections
have constant depth. Either schedule may be used; these are upper
schedules, not optimal leading constants.
Choosing the better one gives the immediate scheduling corollary

```math
D_T=O\!\left(\frac{NL}{b}
+\min\{nL+n^2,\;L\ell_*(n)+n^3\}\right),
\qquad T,G=O(NL),\qquad b\geq2B_0.
```

## 4. Lower bounds and the remaining depth gap

Let $`q=n+a+b`$ be the allowed physical width. A T-depth-d circuit has
at most q T or T-dagger gates in each layer, hence at most $`qd`$ in
total. Applying the existing worst-case T-count lower bound therefore
gives

```math
D_T^\star
=\Omega\!\left(
\frac{\sqrt{NL}}q+\frac Lq+\frac{NL}{q^2}\right).
```

This transfers a count lower bound; it is not a new additive lower bound
for individual source calls. There is also a constant worst-case lower
bound of one T layer. Set all angles except one to zero and choose that
rotation to produce
$`(\cos(\pi/8)|0\rangle+\sin(\pi/8)|1\rangle)\otimes|0\cdots0\rangle`$
on the all-zero input, including the permitted all-zero dirty input.
A Clifford output is a stabilizer state. Projecting all but the displayed
qubit to zero leaves either zero or a subnormalized stabilizer state,
whose overlap with the displayed qubit is at most $`\cos(\pi/8)`$.
The full-output distance is therefore at least
$`\sqrt{2-2\cos(\pi/8)}`$, larger than $`1/64`$, independently of
helper width. One may write the combined lower bound with a maximum of
one and the displayed expression.

These bounds leave a substantial gap. At $`L=N`$ and width
$`q=\Theta(N)`$, they supply only $`\Omega(1)`$ T-depth. The exact
linear T-count of the operator-source primitive does not prove a linear
T-depth lower bound. Its serial schedule is an upper bound for one
implementation.

The [exact source-depth analysis](SOURCE_T_DEPTH.md) sharpens this statement
within a restricted architecture: signed-permutation Majorana Clifford
stages interleaved with disjoint-plane rotations. It certifies the geometric
source and parallelizes the paired source's two tails with a matching
lower bound in that class. These are constant-factor schedules. Arbitrary
Clifford interlayers can leave the Majorana representation, so the restricted
certificate supplies no additional asymptotic lower bound for this theorem.
It also does not constrain approximate replacement sources or establish
optimal controlled-source depth. The [workspace comparison](RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026)
records the clean-initialization and catalyst assumptions of recent shallow
rotation constructions.

At fixed accuracy, with $`L=O(1)`$ and $`b=\Theta(N)`$, the older
layer schedule with depth-optimized banks gives $`D_T=O(n^2)`$ with two clean qubits,
while the worst-case T-count remains at least $`\Omega(\sqrt N)`$.
This illustrates the difference between count and depth; it does not
make the polynomial T-depth optimal.
The [parallel-loader refinement](PARALLEL_DIRTY_LOOKUP.md) achieves
$`T=O(\sqrt N)`$ and $`D_T=O(n^2)`$ together already with sufficiently
large $`\Theta(\sqrt N)`$ dirty workspace at fixed accuracy. The
[dirty-counter hybrid](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters)
improves that depth to $`O(n\log n(\log\log n)^2)`$ with the same
count and workspace orders. Its variable-width extension also enlarges
the matching range at every precision; see the
[hybrid theorem](PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision). The depth lower-bound gap remains open.

At smaller eligible widths, the [amortized loader](AMORTIZED_DIRTY_LOOKUP.md)
instead closes the gap: for fixed L and
$`17B_0\le b\le\sqrt N/n`$, it gives
$`D_T^\star=\Theta(N/b^2)`$ while retaining optimal-order count.
For variable L, matching now also holds when
$`17(L+n+7)\le b\le\sqrt{NL/(nL+n^2)}`$. This includes
$`L=\Theta(n)`$, sufficient $`b=\Theta(n)`$, and
$`D_T^\star=\Theta(N/n)`$ with $`T^\star=\Theta(N)`$.
Large-width depth and the selected high-precision endpoint remain unresolved.

The same complete-frame norm transfers to the existing Hopf-QBP bias
and inverse guarantees. It does not by itself determine the observable
oracle's depth or the total elapsed depth of a complete QBP execution.

## 5. Checks and attribution boundaries

The [focused tests](../tests/test_t_depth.py) check the literal phase
polynomial and inverse of shared-control Fredkin batches, the four
layers' disjoint T targets, and exact dirty-query return on small full
input spaces. They do not establish optimal T-depth or replace the
asymptotic scheduling argument.

The shared-control routing and SelectSwap depth mechanism have primary
precedent in Low, Kliuchnikov, and Schaeffer, cited above. This chapter
supplies an explicit schedule for that mechanism and its composition
with the repository's complete real-frame compiler, including the
two-clean workspace and dirty/reference return contracts. The sharper
predicate schedule uses Khattar and Gidney's existing two-dirty-helper
construction, with the query-pool lifetime checked above. No claim of
a new general lookup primitive, optimal T-depth, or optimal total
elementary depth is made.
