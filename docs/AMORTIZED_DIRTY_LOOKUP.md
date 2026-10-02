# Amortized dirty lookup and matching T-depth ranges

[Batched construction](BATCHED_DIRTY_LOOKUP.md) · [Depth model and routing](T_DEPTH_COMPILER.md) · [Complete-frame source](OPERATOR_SOURCE_COMPILER.md)

The low-address indicator router can be shared across every chunk of a
dirty table query. The required change is to select a linear map on the
dirty indicator, then cancel its unknown input with one outer echo.
The existing two-pass dirty traversal selects that map with constant
T-depth per chunk. This removes both logarithms that were previously
paid for every guarded batch.

**Variable-accuracy theorem.** Let $`n\ge1`$, $`N=2^n`$, and
$`0\lt\eta\le1/64`$. Set
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and $`B_0=L+n+7`$.
For every $`b\ge17B_0`$, the prescribed complete real Hopf frame has
a coherent Clifford+T implementation using two initialized compiler
flags and at most b arbitrary dirty qubits, satisfying

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL+n^2\right).
```

These bounds hold on the same circuit. The norm includes initialized
flag leakage, arbitrary dirty inputs, and their reference correlations.
Every coherent query and actual inverse is implemented and charged.
There are no measurements, resets, initialized indicators, supplied
catalysts, or quantum memory oracles.

When the following interval is nonempty, the count and depth both have
the optimal worst-case order:

```math
17B_0\le b\le\sqrt{\frac{NL}{nL+n^2}},\qquad
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

At fixed L, the separate [fixed-accuracy corollary](#fixed-accuracy-corollary)
retains $`T=O(\sqrt N+N/b)`$ throughout $`b\ge17B_0`$ and matching
depth on the original literal interval $`b\le\sqrt N/n`$. The
[capped precision allocation](#capping-the-source-precision) gives
accumulated source depth $`O(nL+n\log(n+1))`$; the retained query
routing contributes $`O(n^2)`$. The separate
[dirty-counter hybrid](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters)
improves fixed-accuracy depth to $`O(n\log^2(n+2))`$ at sufficient
$`b=\Theta(\sqrt N)`$, retaining optimal-order T-count.
Large-workspace optimal depth and the
high-precision constant-clean endpoint remain open. T-depth permits arbitrary
Clifford circuits between T layers; their elementary depth is not
bounded by this theorem, and their gate count remains included in G.
The earlier batched construction remains a valid separate proof.

## 1. A controlled linear shear has constant T-depth

Let Y and Z be disjoint arbitrary registers of s and w bits. For a
classically specified binary matrix $`A\in\mathbb F_2^{w\times s}`$,
define the Clifford shear

```math
C_A:(Y,Z)\longmapsto(Y,Z\oplus AY).
```

It is an involution. Every pair of these shears commutes, even for
different A, since Y is unchanged and their contributions to Z add
over $`\mathbb F_2`$.

Let h be a separate, arbitrary control bit. The controlled shear is

```math
C_h(C_A):(h,Y,Z)\longmapsto(h,Y,Z\oplus hAY).
```

Put $`\rho=\mathrm{rank}(A)`$. Row and column elimination gives
invertible binary matrices P and Q with

```math
PAQ=J_\rho,
```

where $`J_\rho`$ is the rectangular matrix with its first rho diagonal
entries equal to one and all other entries zero. A pivot uses at most
$`O(s+w)`$ row and column additions and swaps. Thus P, Q, and their
inverses have CNOT circuits with total size $`O(\rho(s+w))`$.
This bound includes the three-CNOT realization of a wire swap.
No auxiliary register is needed for this elimination circuit.

Apply the separate basis changes
$`(Y,Z)\mapsto(Q^{-1}Y,PZ)`$, the controlled shear for
$`J_\rho`$, and their actual inverses. The middle operation consists
of rho Toffolis with common control h and disjoint pairs
$`(Y_i,Z_i)`$. It adds $`hJ_\rho Q^{-1}Y=hPAY`$ in the changed
Z basis, so the complete word is exactly $`C_h(C_A)`$.

The [shared-control phase-polynomial schedule](T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
implements this middle Toffoli batch with

```math
D_T\le4,\qquad T=6\rho+(\rho\bmod2),\qquad G=O(\rho).
```

It uses only h and the existing paired Y and Z wires. In particular,
it requires no initialized copies of h and no dirty fanout helpers.
The surrounding basis changes are Clifford. The complete controlled
shear therefore has

```math
T=O(\rho),\qquad D_T\le4,\qquad
G=O(\rho(s+w))=O(sw),
```

with zero work beyond h, Y, and Z. For $`\rho=0`$ it is identity.
Both the CNOT basis changes and the exact Toffoli phase polynomial
retain literal scalar phases. Although Y can change inside a leaf's
basis changes, the completed controlled shear restores it exactly.

## 2. Selecting a shear with dirty unary traversal

Let a p-bit address c select one of $`K=2^p`$ matrices $`A_t`$.
Reserve p arbitrary dirty selector bits $`W_1,\ldots,W_p`$,
disjoint from c, Y, and Z. We implement

```math
F:(c,Y,Z,W)\longmapsto(c,Y,Z\oplus A_cY,W).
```

For $`p\ge1`$, use the existing
[two-pass dirty traversal](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter),
replacing a leaf's CNOT mask by a controlled shear. At a leaf t, write
$`\ell_i=[c_i=t_i]`$. On entry to its depth-first path, update

```math
v_1=W_1\oplus\ell_1,\qquad
v_i=W_i\oplus v_{i-1}\ell_i\quad(2\le i\le p).
```

The first update is a CNOT with an optional negative literal; later
updates are exact Toffolis with an optional negative address literal.
At leaf t, apply $`C_{v_p}(C_{A_t})`$. Undo each selector update
when leaving its subtree. The complete traversal returns every selector.
Its action on Y and Z is the product, over leaves in traversal order,
of $`C_{A_t}^{v_p(t;c,W)}`$.

The recurrence expands as

```math
v_p(t;c,W)=\prod_{i=1}^p\ell_i\ \oplus\ h_t(c,W),
```

where $`h_t`$ is the same recurrence with the first-selector update
omitted. Perform that second traversal with identical leaf insertions,
but omit those first-selector updates, and take its actual circuit
inverse. Its selectors also return, and its completed leaf exponents
are $`h_t(c,W)`$. All completed shears commute and square to identity.
For p=1 this second traversal still visits both leaf labels, controlled
by the unchanged $`W_1`$; omitting the root updates does not omit any
leaf insertion.
Consequently the first traversal followed by the inverse of the second
has action

```math
\prod_t C_{A_t}^{v_p(t;c,W)\oplus h_t(c,W)}
=\prod_t C_{A_t}^{[c=t]}=C_{A_c}.
```

This is the claimed F on every basis input. The proof does not assume
any selector or any Y or Z bit is zero. Each leaf leaves c and the
selectors untouched, so its internal basis changes do not alter the
selector recurrence. Literal phase-free identities extend the equality
to arbitrary superpositions and reference-entangled inputs. F is a
semantic involution; whenever an inverse word is used below, it is the
actual inverse of its complete native circuit.

A complete binary depth-first traversal has $`O(K)`$ nodes, so it
uses $`O(K)`$ selector updates, rather than a separate p-step
predicate computation at every leaf. Each leaf occurs twice in F.
Writing $`\rho_t=\mathrm{rank}(A_t)`$, Section 1 gives

```math
T(F)=O\!\left(K+\sum_t\rho_t\right),\qquad
D_T(F)=O(K),\qquad
G(F)=O\!\left(K+\sum_t\rho_t(s+w)\right).
```

All p selectors return before the next operation. For $`p=0`$,
use the direct Clifford shear $`C_{A_0}`$, with no selectors and zero
T-count or T-depth. For $`p\ge1`$, an explicit sequential ledger
using four-layer, seven-T Toffolis has

```math
\begin{aligned}
T(F)&=7(8K-16)+2\sum_t[6\rho_t+(\rho_t\bmod2)],\\
D_T(F)&\le4(8K-16)+8\,\#\{t:\rho_t>0\}.
\end{aligned}
```

Each traversal has $`4K-8`$ nonroot selector Toffolis, including
compute and uncompute. The first-selector operations are Clifford.
These expressions describe this unoptimized literal schedule; the
asymptotic proof does not require additional parallelism between leaves.

## 3. One outer indicator echo and a completed query

Consider a Q-row, m-bit XOR table $`f`$, with $`Q=2^r`$ and
$`m\ge1`$. A non-power-of-two table can be padded by zero rows.
Choose powers of two $`1\le\mu\le Q`$ and
$`1\le s\le H=Q/\mu`$. Put

```math
w=\mu m,\qquad K=H/s,\qquad p=\log_2K.
```

Split the address as $`(c,u,v)`$, selecting the chunk, row within
its chunk, and word bank, respectively. Let Z consist of w arbitrary
dirty word-bank bits, and let Y consist of s arbitrary dirty indicator
bits. For each t, column u of $`A_t`$ concatenates
$`f(t,u,v)`$ over all banks v. These are classical binary matrices;
their implemented circuits are included in the ledger.

The unguarded routed indicator

```math
I_u:\quad\mathcal R_u,\ X_{Y_0},\ \mathcal R_u^\dagger
```

maps $`Y\mapsto Y\oplus e_u`$, where $`\mathcal R_u`$ routes
the selected one-bit bank to position zero. It needs no extra work and
has $`T,G=O(s)`$ and $`D_T=O(\log_2s)`$; for s=1 it is simply
an X gate. Use the selected shear F from Section 2 in the chronological
word

```math
\mathcal L:\quad F,\ I_u,\ F^\dagger,\ I_u^\dagger.
```

For an initial indicator word y, the two contributions to Z are
$`A_cy`$ and $`A_c(y\oplus e_u)`$. Their sum is $`A_ce_u`$.
Thus this loader writes every selected high-address table word into Z,
returns Y and every selector, and is exact on the complete input space.
The indicator router is paid only twice for the whole loader. There is
no chunk-controlled conjunction gate in this word.

Complete the query with the established SelectSwap echo. Let
$`\mathcal R_v`$ route the chosen m-bit bank to position zero, and
let C copy that bank by CNOTs into the separate arbitrary m-bit output.
Apply

```math
\mathcal L,\ \mathcal R_v,\ C,\ \mathcal R_v^\dagger,\
\mathcal L^\dagger,\ \mathcal R_v,\ C,\ \mathcal R_v^\dagger.
```

The output receives its original value XOR $`f(c,u,v)`$. All banks,
indicators, and selectors return exactly. The loader inverse reverses
the entire native word and is used only after the word-bank route has
been undone. Its targets exclude the query output. Inactive zero rows,
arbitrary dirty contents, and all reference correlations obey the same
literal equality.

Since $`\sum_t\rho_t\le Ks=H`$ and $`Ksw=Qm`$, the completed
query satisfies

```math
\begin{aligned}
T_{\rm query}&=O(H+w),\\
G_{\rm query}&=O(Qm),\\
D_{T,\rm query}&=O(K+\log_2s+\log_2\mu).
\end{aligned}
```

The depth includes both dirty traversals, both outer indicator calls,
the inverse loader, and every word-bank route, each with constant
multiplicity. The Clifford estimate charges every matrix basis change:
$`\rho_t(s+w)\le2sw`$. This is not an uncharged linear-algebra oracle.

| Additional live register | Dirty width | Lifetime and return |
|---|---|---|
| Word banks Z | $`w=\mu m`$ | Live throughout the completed query; returned exactly |
| Indicator Y | s | Live across the outer echo; returned exactly |
| Traversal selectors W | $`p=\log_2K`$ | Returned after each completed F or its inverse |

The literal additional dirty width is $`w+s+p`$. The address, separate
query output, existing source core, flags, and previous reservations
are excluded from this table. No helper wire is taken from a live bank
or a live indicator. The p=0 case needs no traversal selector.

## 4. Allocation at a prescribed dirty width

Let B be the extra dirty width, excluding the output and all existing
reservations. Assume

```math
B\ge16(m+r+1),\qquad r=\log_2Q.
```

Use the same count-efficient allocation as the batched proof. For
$`Q\ge m`$, choose mu as the largest power of two at most
$`\min\{\sqrt{Q/m},B/(4m)\}`$; for $`Q\lt m`$, set mu=1.
Then choose s as the largest power of two at most
$`\min\{H,B/4\}`$. The parameters divide Q and obey

```math
w\le B/4,\qquad s\le B/4,\qquad p\le r\le B/16,\qquad
w+s+p\le B.
```

Unused wires remain untouched. For $`Q\ge m`$, rounding gives
$`H=O(\sqrt{Qm}+Qm/B)`$ and $`w\le\sqrt{Qm}`$.
For $`Q\lt m`$, we have $`H=Q\lt m\le B/16`$, hence s=H
and K=1. Therefore

```math
T_{\rm query}=O\!\left(\sqrt{Qm}+m+\frac{Qm}{B}\right),\qquad
G_{\rm query}=O(Qm).
```

If K=1, the selected shear is Clifford and the remaining T-depth is
$`O(\log_2(Q+1))`$. If $`K\gt1`$, rounding forces
$`s\gt B/8`$ and the case $`Q\lt m`$ is impossible. If the
width cap determines mu, then $`\mu\gt B/(8m)`$ and
$`K\lt64Qm/B^2`$. If the count cap determines mu, then
$`H\lt2\sqrt{Qm}`$. Combining $`H\gt s\gt B/8`$ with this
last inequality gives $`\sqrt{Qm}/B\gt1/16`$ and again
$`K=O(Qm/B^2)`$. Thus

```math
D_{T,\rm query}=O\!\left(\log_2(Q+1)+\frac{Qm}{B^2}\right).
```

When Q=1, the sole row may instead be emitted directly as an X word,
with zero T-count and T-depth, $`O(m)`$ Clifford gates, and no work.
The query lemma allows arbitrary m. The complete-frame application below
uses that freedom for a precision-dependent output word.

## 5. Complete-frame composition and the matching range

Use the layerwise two-flag compiler's certified coefficient encoding,
native source, suffix echo, and actual-inverse amplification, with the
capped precisions proved below. At depth $`d=0,\ldots,n-1`$, its
fixed number of whole-word queries have $`Q_d=2^{d+2}`$ rows and
address width $`r_d=d+2`$. All bounds in this section hold uniformly
for $`L\ge6`$.

### Capping the source precision

The quadratic sum of source widths is not required by the complete-frame
error contract. It can be reduced by changing only the precision assigned
to each existing layer. Put $`k=n-d`$ and define

```math
h=\lceil\log_2(8n)\rceil,\qquad
\widetilde m_d=L+4+\min\{k,h\}.
```

Use the same certified coefficient encoding, native source, suffix echo,
and actual-inverse amplification at this smaller width. Equation (24) of
the [operator-source proof](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
gives the full initialized-isometry error of one layer as

```math
\epsilon_d\le10\sqrt2\,2^{-\widetilde m_d}.
```

This is the error including rejected flags, dirty-core disturbance, and
reference correlations. No square-root conversion or projection onto a
successful branch is used. Separating the uncapped geometric terms from
the capped tail gives

```math
\begin{aligned}
\sum_{d=0}^{n-1}\epsilon_d
&\le\frac{5\sqrt2}{8}\,2^{-L}
  \sum_{k=1}^{n}2^{-\min\{k,h\}}\\
&\lt\frac{5\sqrt2}{8}\,2^{-L}(1+n2^{-h})\\
&\le\frac{45\sqrt2}{64}\,2^{-L}\lt2^{-L}\le\eta.
\end{aligned}
```

The last strict inequality follows from $`4050\lt4096`$. Complete-input
telescoping therefore proves the same frame approximation and QBP
substitution guarantees. Intermediate work is never assumed to have reset.

Every new source width is at most the old $`L+n-d+4`$. In particular,
$`\widetilde m_d+r_d+1\le B_0`$. Keep the entire reserved base of
$`B_0=L+n+7`$ dirty qubits and the two separate clean flags; unused base
wires remain untouched. The extra-pool allocation, its sufficient threshold
$`b\ge17B_0`$, and each query's exact return consequently remain valid.
The old weighted bounds also hold pointwise:

```math
\sum_dQ_d\widetilde m_d=O(NL),\qquad
\sum_d\sqrt{Q_d\widetilde m_d}=O(\sqrt{NL}).
```

Writing $`s=\min\{n,h\}`$, the unweighted sum is now exactly

```math
\sum_d\widetilde m_d
=n(L+4)+ns-\frac{s(s-1)}2
=O\!\left(nL+n\log(n+1)\right).
```

Thus the retained serial source words use
$`O(nL+n\log(n+1))`$ T gates and T-depth. At fixed L, their depth
is $`O(n\log(n+1))`$; the retained routers still pay
$`O(\log Q_d)`$ at each layer, totaling $`O(n^2)`$. This identifies
the remaining quadratic contribution in this construction; it is not a
lower bound for another lookup or frame circuit.

The allocation is also near-optimal for its particular additive error
certificate. If arbitrary assigned precisions $`m_1,\ldots,m_n`$
satisfy $`\sum_k10\sqrt2\,2^{-m_k}\le2^{-L}`$, convexity gives

```math
\frac1n\sum_km_k
\ge L+\log_2 n+\log_2(10\sqrt2).
```

The cap above has average precision at most $`L+\log_2 n+8`$ and
obeys every original layerwise width cap. Its total assigned precision
is therefore within $`5n`$ bits of the best allocation certified by this
inequality, including under those width caps. This is a budgeting result,
not a lower bound on the actual accumulated error, joint source synthesis,
or unrestricted T-depth. Changing the error analysis or circuit remains
eligible to improve the full-frame bound.

### Workspace and the simultaneous resource ledger

Reserve the entire old base $`B_0=L+n+7`$, including the source core,
old selectors, and separate dirty suffix control, even when some of its
wires are unused. The two initialized flags are separate. The extra pool
has size

```math
B=b-B_0\ge16B_0\ge16(\widetilde m_d+r_d+1).
```

Thus Section 4 supplies every query's banks, indicators, and traversal
selectors outside the occupied base. Reuse the pool only after the
completed query or its actual inverse. Between queries, two returned
extra-pool wires supply the
[logarithmic-depth suffix predicates](T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
No source or predicate borrows a live register from an incomplete query.
Each query implements the exact table unitary at its assigned precision
on arbitrary inputs, so replacing its implementation changes neither
the amplification identity nor the error estimate above.

The weighted sums proved above and $`B=\Theta(b)`$ give total query
T-count

```math
O\!\left(\sum_d\sqrt{Q_d\widetilde m_d}
+\sum_d\widetilde m_d
+\frac1B\sum_dQ_d\widetilde m_d\right)
=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL+n\log(n+1)\right).
```

The sources have the same unweighted precision charge. Suffix predicates
have $`O(n^2)`$ T-count and $`O(n\log(n+1))`$ T-depth; the two-flag
reflections have constant cost per use. Since $`n^2=O(\sqrt N)`$ and
$`L\ge6`$, all of these count terms except $`nL`$ are absorbed into
$`O(\sqrt{NL})`$. The query Clifford count is
$`O(\sum_dQ_d\widetilde m_d)=O(NL)`$; source, predicate, and reflection
Cliffords also fit that bound. Consequently the same circuit has

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL).
```

For depth, Section 4 and the source and predicate schedules give

```math
\begin{aligned}
D_T
&=O\!\left(\sum_d\log_2(Q_d+1)
+\frac1{B^2}\sum_dQ_d\widetilde m_d
+nL+n\log(n+1)\right)\\
&=O\!\left(\frac{NL}{b^2}+nL+n^2\right).
\end{aligned}
```

The bank cap prevents excess available workspace from worsening the
T-count. There is no upper restriction on b for these upper bounds.
The additive $`nL`$ is retained: count optimality is not asserted for
every variable-precision allocation covered by the upper theorem.

### Matching variable-precision window

The [complete real-frame lower bound](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
is $`\Omega(\sqrt{NL}+L+NL/q)`$ for all $`L\ge6`$ and physical
widths $`q=n+2+b`$. Its fixed-width packing proof first obtains the
$`NL/q`$ term when $`q^2\le c_0NL`$; outside that range the inherited
$`\sqrt{NL}`$ lower bound already dominates $`NL/q`$. Thus the
combined bound has no additional precision or width restriction.
Here $`q=\Theta(b)`$, and every T layer uses at most q T gates, so

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
D_T^\star=\Omega\!\left(\frac{\sqrt{NL}}b+\frac Lb
+\frac{NL}{b^2}\right).
```

These are worst-case bounds over the prescribed complete real frames;
the depth lower bound permits arbitrary Clifford interlayers. Suppose

```math
17B_0\le b\le\sqrt{\frac{NL}{nL+n^2}}.
```

Then $`nL+n^2\le NL/b^2`$. Also $`b\le\sqrt{NL}`$, so
$`\sqrt{NL}\le NL/b`$, and $`nL\le NL/b^2\le NL/b`$.
The upper and lower bounds therefore give

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

The displayed upper count and depth are achieved together. There is no
separate restriction $`L\le N/n^2`$ in this matching window: its depth
condition already absorbs the extra source count. The window is a
sufficient range and may be empty for a given n and L.

### Fixed-accuracy corollary

For fixed L, $`nL=O(\sqrt N)`$, and the same circuit has

```math
T=O\!\left(\sqrt N+\frac Nb\right),\qquad G=O(N),\qquad
D_T=O\!\left(\frac N{b^2}+n^2\right)
\quad(b\ge17B_0).
```

Its count has optimal worst-case order throughout that width range.
Retain the original literal matching-depth interval separately:
$`17B_0\le b\le\sqrt N/n`$ implies $`n^2\le N/b^2`$ and
$`\sqrt N\le N/b`$. Hence the simultaneous worst-case orders there are
$`T^\star=\Theta(N/b)`$ and $`D_T^\star=\Theta(N/b^2)`$.
This corollary does not require inclusion in the displayed
variable-precision window with its particular constants.

### Inverse-polynomial accuracy at linear dirty width

Fix $`\alpha\gt0`$ and take $`\eta=N^{-\alpha}`$ for sufficiently
large n, so $`L=\Theta(n)`$. Choose $`b=\Theta(n)`$ with a sufficiently
large constant to satisfy $`b\ge17B_0`$. Exponential growth of N then
places this allocation in the matching window. The same complete-frame
circuit has the simultaneous worst-case orders

```math
T^\star=\Theta(N),\qquad D_T^\star=\Theta(N/n),\qquad G=O(Nn).
```

More generally, for fixed c and $`6\le L\le cn`$, sufficiently large
n and eligible $`b=\Theta(n)`$ give
$`T^\star=\Theta(NL/n)`$ and $`D_T^\star=\Theta(NL/n^2)`$.
These are prescribed complete real-frame results, with the same QBP
substitution contract; no new state-only or complex-frame theorem is
asserted. At large workspace the depth gap remains. In particular,
the allocation $`L=N,b=L+n+7`$ does not meet this theorem's bank
threshold, and the variable-accuracy extension does not resolve the
high-precision constant-clean endpoint. The restricted Majorana source
certificate is not promoted to an unrestricted depth lower bound.

## 6. Attribution and proof scope

Low, Kliuchnikov, and Schaeffer,
[*Trading T gates for dirty qubits in state preparation and unitary synthesis*](https://arxiv.org/html/1812.00954v2),
supply the SelectSwap and dirty-indicator framework. Khattar and Gidney,
[*Rise of conditionally clean ancillae for optimizing quantum circuits*](https://arxiv.org/html/2407.17966v1),
Section 7.3 of v1, supply dirty unary iteration; Section 5.4 supplies
the returned two-dirty-helper conjunction used by the unchanged suffix
schedule. The repository's explicit two-pass traversal gives the selector
cancellation used above, and its shared-control phase-polynomial schedule
gives the controlled-shear leaf.

Kim and Laakkonen,
[*Any Clifford+T circuit can be controlled with constant T-depth overhead*](https://arxiv.org/pdf/2512.24982v1),
Theorem 3, already give constant-Toffoli-depth control of a CNOT circuit
without ancillary qubits. Section 1 is an explicit rank-sensitive
specialization for this shear, not a new constant-depth control
principle. Boyd,
[*Low-Overhead Parallelisation of LCU via Commuting Operators*](https://arxiv.org/html/2312.00696v2),
Section III and Appendix A, uses commuting-operator grouping and basis
changes in SELECT and QROM constructions. Its copied-address registers
have initialized inputs, so that construction does not directly provide
the arbitrary-dirty workspace bound proved here.

The composition proved here selects a commuting linear shear before
the indicator echo, charges its Clifford basis changes, and applies the
resulting exact query to the prescribed complete real frame. It makes
no claim that the generic lookup identity is new in the literature.
The [focused checks](../tests/test_amortized_dirty_lookup.py) exhaust
154 small rectangular binary matrices, verify symbolic dirty traversal
and complete-query identities, and test literal native phases, actual
inverses, and emitted resource ledgers. A missing-second-traversal
counterexample protects the arbitrary-selector contract. These bounded
fixtures do not establish the asymptotic or matching-range statements;
those follow from the full-input identities and ledgers above.
