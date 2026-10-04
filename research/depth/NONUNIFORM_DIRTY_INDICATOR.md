# A linear-size dirty indicator with logarithmic address depth

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Dirty tree and echoes](../../docs/CHUNKED_DIRTY_INDICATOR.md) · [Read-only conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences) · [Rectangular blocked allocation](../../docs/UNIFORM_PRECISION_DEPTH.md#2-an-unequal-indicator-allocation-preserves-the-block-depth-bound) · [Current research checkpoint](../../docs/OPEN_PROBLEM.md#current-decision)

Unequal address chunks remove the polynomial size overhead of the earlier
uniform-chunk indicator while retaining logarithmic address T-depth. The
large conjunctions occur near the root, where there are fewer tree edges.
The final high-multiplicity stage has bounded arity. The existing dirty
root and leaf echoes still return every arbitrary helper exactly.

For every integer $`s\ge0`$, put $`S=2^s`$. There is an exact
Clifford+T circuit

```math
|x,Y,W\rangle\longmapsto|x,Y\oplus e_x,W\rangle
```

on an s-bit address x, an arbitrary S-bit output Y, and arbitrary dirty
work W, with absolute constants in

```math
T,G,w=O(S),\qquad D_T=O(\log_2(s+2)).
```

Here w is additional dirty width beyond x and Y; G counts elementary
Clifford gates. No input wire is assumed initialized. The equality keeps
literal phase and holds with arbitrary reference correlations. Every
inverse used below is the actual reversed native word. T-depth allows
intervening Clifford circuits of unrestricted elementary depth, whose
count is included in G. For s equal to zero, the indicator is a single
Clifford X and needs no work.

The associated blocked query retains its square-root count term in Qm
and its stated width penalty without a polynomial indicator overhead.
This removes the tail's
linear-depth indicator allowance in the current large-width frame
schedule. Four separately charged early contributions remain; no new
complete-frame asymptotic depth bound follows from this lemma alone.
The subsequent [protected-source theorem](PROTECTED_UNARY_SOURCE.md)
also removes repeated early source boundaries; logical stages, outer
predicates, and program queries retain their separate linear allowances.

## 1. The exact tree interface permits unequal chunks

Fix any positive chunk lengths $`b_1,\ldots,b_t`$ summing to s.
Let $`d_i=\sum_{j=1}^i b_j`$ and $`d_0=0`$. Reserve an arbitrary
dirty root and one dirty node for each binary prefix of length $`d_i`$,
including all S leaf nodes. These are distinct from Y and the address.
At stage i, for each parent prefix p and each next-chunk label v, perform

```math
w_{pv}\longmapsto w_{pv}\oplus w_p[\,x^{(i)}=v\,].
```

The conjunction has exactly $`b_i+1`$ original controls: its parent
bit and the current chunk's literals. It does not recompute a conjunction
of the entire preceding address. Use the proved read-only conjunction
with a private arbitrary dirty helper pool for each edge.

There is an absolute constant c such that an edge's T-count, elementary
Clifford count, and additional dirty width are each at most
$`c(b_i+3)^3`$, and its T-depth is at most
$`c\log_2(b_i+3)`$. These safe majorants follow from the existing
read-only conjunction theorem; the extra constant in the argument counts
the parent control. Original controls enter through its read-only
Clifford interface, and every private helper returns exactly.

At one stage all children and private pools are disjoint. Parents belong
to the preceding boundary and are not targeted at that stage. Therefore
all edge T layers can run together despite shared parent and address
controls. Their Clifford interactions may be serialized and remain
charged. Complete a stage before starting the next, and reuse its
private pool only after return. All tree nodes remain live.

Call the chronological top-down word F. For fixed x its completed action
on the entire arbitrary tree word is linear, $`w\mapsto A_xw`$,
and

```math
A_xe_\emptyset=e_\emptyset+\sum_{i=1}^t e_{x_{1:d_i}}.
```

A unit difference at the root propagates to precisely one node per
boundary. This is a difference of two arbitrary tree inputs and assumes
no clean node. Unequal boundary spacings do not change the argument.

Define the algebraic conjugate $`K=F X_\emptyset F^\dagger`$,
implemented chronologically by $`F^\dagger,X_\emptyset,F`$.
It adds $`A_xe_\emptyset`$ to any tree word. Let C copy each leaf
into its corresponding Y bit by CNOT. The chronological echo

```math
C,\quad K,\quad C,\quad K^\dagger
```

therefore adds exactly $`e_x`$ to Y and returns the entire tree.
Each edge already returned its private work. Literal phase-free basis
identities give the same equality on arbitrary superpositions and
references. The complete indicator uses F or its actual inverse four
times, not once for each chunk. Its extra Clifford operations are two
leaf-readout words and two root X gates.

## 2. Choose the remaining address sizes

Let $`r_0=s`$. Whenever $`r_i\gt64`$, define

```math
r_{i+1}=\left\lceil4\log_2(r_i+2)\right\rceil,
\qquad b_{i+1}=r_i-r_{i+1}.
```

Once $`1\le r_i\le64`$, take the final chunk of size $`r_i`$
and set its next remainder to zero. An initial zero remainder uses the
single-output Clifford case instead.

For every real $`r\ge64`$,

```math
4\log_2(r+2)+1\le r/2.
```

It holds at 64 and the right side grows faster thereafter. Thus every
nonfinal rounded remainder satisfies $`r_{i+1}\le r_i/2`$, and
every chosen chunk is positive. The process terminates; all boundaries
are valid integer address lengths. Its final chunk has at most 64 bits.

There are exactly $`2^{s-r_{i+1}}`$ edges at a stage ending with
remainder $`r_{i+1}`$. This count includes every child label and
every preceding parent. It already includes the factor $`2^{b_{i+1}}`$
for the chunk's equality tests; that factor is not multiplied in again.

## 3. Linear count and a complete simultaneous-work bound

For a nonfinal stage, abbreviate $`r=r_i`$, $`r'=r_{i+1}`$,
and $`b=r-r'`$. Its normalized cubic cost satisfies

```math
\frac{2^{s-r'}(b+3)^3}{S}
=\frac{(b+3)^3}{2^{r'}}
\le\frac1{r+2}.
```

Indeed $`2^{r'}\ge(r+2)^4`$, and $`r'\ge1`$ implies
$`b+3\le r+2`$. The nonfinal remainders are integers greater
than 64 and decrease by at least a factor two. Reading their sequence
backward from its last value gives

```math
\sum_{\rm nonfinal}\frac1{r_i+2}
\le\sum_{\rm nonfinal}\frac1{r_i}\le\frac2{65}.
```

The final stage has S edges and $`b\le64`$, so its normalized
cubic cost is at most $`67^3`$. Consequently

```math
\sum_i2^{s-r_{i+1}}(b_{i+1}+3)^3
\le\left(67^3+\frac2{65}\right)S.
```

This is a constant geometric budget. There is no uncharged factor from
the number of chunks or an iterated logarithm. The T and elementary
Clifford counts of F are at most c times this budget. Four F calls and
the $`2S+2`$ extra Clifford gates prove the stated linear indicator
counts. The constants are sufficient bounds, not an optimized crossover.

The live tree has

```math
1+\sum_{i=1}^t2^{d_i}\le2S-1
```

nodes, since the $`d_i`$ are distinct increasing positive integers
at most s. The largest private conjunction pool is at most
$`c67^3S`$: the final stage fits this bound, and each nonfinal
stage's normalized cost is at most $`1/(r_i+2)`$. A single pool
of this maximum size suffices because every edge returns its work before
the stage finishes. It is disjoint from every live tree node, output,
and address. Hence additional dirty work is at most

```math
(2+c67^3)S.
```

This counts all simultaneous work, including nodes retained from earlier
stages and the root used by the conjugation. The inverse traversal has
the same peak. No temporary clean node or erased-but-still-live pool is
assumed.

## 4. Logarithmic depth, including rounding

The depth of F is at most
$`c\sum_i\log_2(b_{i+1}+3)`$. The stronger decrease

```math
\left\lceil4\log_2(r+2)\right\rceil+2\le\sqrt{r+2}
\qquad(r\ge4096)
```

bounds this sum. To include the ceiling, compare
$`4\log_2(r+2)+3`$ with $`\sqrt{r+2}`$: the inequality holds
at 4096 and their difference grows thereafter. Thus, while the remainder
is at least 4096, $`\log_2(r_i+2)`$ decreases by at least a factor
two. Since $`b+3\le r_i+2`$ at every nonfinal stage, those stages
contribute at most $`2\log_2(s+2)`$ to the depth sum.

Once the remainder is below 4096, at most one nonfinal stage remains:
for $`64\lt r\lt4096`$, the next remainder is at most 49, hence
at most 64. It is followed by the bounded final stage. These contribute
only an absolute constant. If the initial remainder is already at most
64, its one stage has depth $`O(\log_2(s+2))`$ directly.

Therefore

```math
\sum_i\log_2(b_{i+1}+3)=O(\log_2(s+2)).
```

The four traversals in the exact echo only multiply this by four.
Root flips and leaf CNOTs have zero T-depth. This proves the introductory
indicator theorem together with its complete-input return contract.

## 5. A count-efficient rectangular blocked-query consequence

Consider an arbitrary Q-row, m-bit table with $`Q=2^r`$ and
$`m\ge1`$. B denotes additional dirty width, excluding its original
address, output, and all other compiler reservations. Choose an absolute
$`c_I\ge1`$ bounding the new indicator's T-count, Clifford count,
and width including its output by $`c_IS`$; its depth is
$`O(\log_2(s+2))`$ for s address bits.

Under the sufficient condition

```math
B\ge16c_I(r+2),
```

use the existing rectangular blocked allocation:

```math
J=\operatorname{pow2floor}
  \min\left\{Q,\sqrt{Qm},\frac{B}{8c_I}\right\},
\qquad H=\min\{J,Q/J\},\qquad K=Q/(HJ).
```

Here pow2floor is the largest power of two no greater than its argument.
For $`Q=1`$, emit the sole row as a Clifford X word instead. The
indicator pools occupy at most $`c_I(H+J)\le B/4`$ wires, and
the high-address selectors and one dirty controlled-bilinear helper use
at most $`r+1`$ more. Thus all simultaneous work fits B.

The [proved rectangular ledger](../../docs/UNIFORM_PRECISION_DEPTH.md#2-an-unequal-indicator-allocation-preserves-the-block-depth-bound)
now has constant indicator overhead. The selected middle word costs
$`O(Qm/J)`$ T gates, $`O(Qm)`$ Cliffords, and $`O(Km)`$
T-depth. Rounding gives

```math
Qm/J=O\!\left(m+\sqrt{Qm}+Qm/B\right),\qquad
K=O(1+Q/B^2).
```

Each indicator has at most r address bits, so the two indicator depths
sum to $`O(\log_2(r+2))`$. Both actual inverses are charged. Thus
one exact arbitrary-input XOR query has

```math
T=O\!\left(m+\sqrt{Qm}+\frac{Qm}{B}\right),\qquad
G=O(Qm),
```

```math
D_T=O\!\left(m\left[1+\frac Q{B^2}\right]+\log_2(r+2)\right).
```

The indicator T-count fits because $`H+J\le2J\le2\sqrt{Qm}`$.
Their Clifford count is absorbed because $`H+J\le2Q`$
and $`m\ge1`$. All returned helpers, literal phases, coherent
addresses, arbitrary outputs, and zero inactive rows retain the generic
blocked-query contract. No separate word-bank router is present.

For any fixed m and sufficiently large fixed $`B=C\sqrt Q`$, this
query has $`T=O(\sqrt Q)`$, $`G=O(Q)`$, and
$`D_T=O(\log_2(r+2))`$. Its dirty width is $`O(\sqrt Q)`$.
This is a local query result; it does not assert a matching query-depth
lower bound.

## 6. The tail improves; the complete-frame frontier is unchanged

Use the uniform low-precision unary schedule with
$`6\le L\le\log_2(n+2)/16`$, source parameter $`L'=L+2`$,
and its natural cutoff M. Its proof gives

```math
M[L+\log(n+2)]=o(n)
```

uniformly over that range. On every remaining layer,

```math
Q_k=4N2^{-k},\qquad
m_k=L'+4+\min\{k,\lceil\log_2(8n)\rceil\},\qquad
r_k=n-k+2.
```

Where the extra pool satisfies $`B\ge16c_I(n+3)`$, replace the
tail queries by Section 5. The original full-input source and capped
error certificate are unchanged. Using the existing weighted sums
$`\sum_kQ_km_k=O(NL)`$ and
$`\sum_k\sqrt{Q_km_k}=O(\sqrt{NL})`$, the tail has

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{B}+M[L+\log(n+2)]\right),
\qquad G=O(NL),
```

```math
D_{T,\rm tail}
=O\!\left(\frac{NL}{B^2}+M[L+\log(n+2)]\right).
```

The latter also includes the old serial sources, reflections, and suffix
predicates. Their other local gate counts are polynomial in n uniformly
in this low-precision range and fit the stated square-root T-count and
$`NL`$ Clifford budgets. Whenever $`B\gg\sqrt{NL/n}`$ in this eligible regime,
the tail is $`o(n)`$; sufficiently large $`B=\Theta(\sqrt{NL})`$
is one such case. In particular, the fixed-accuracy square-root-width
tail no longer contributes its previous linear indicator allowance.
Exact query replacement adds no approximation error or reset assumption.

The early sequential logical/selector stages, phase-source boundaries,
outer predicates, and program load/unload remain separately charged at
$`O(n)`$ in the current upper ledger. Applying the new indicator to
an early prefetch still gives $`O(\log(n+2))`$ depth per group;
it does not reduce the number of groups. Hence this result does not
improve the complete-frame $`O(n)`$ fixed-accuracy bound by itself,
nor either matching interval of the uniform-precision theorem.

Section 5's sufficient width constant is explicit and may exceed the
old query reservation constant. There is no claim that this replacement
fits every allocation at the literal threshold $`17(L+n+7)`$.
The established full-frame theorems retain that threshold unchanged:
where the new query does not fit, keep their proved query schedule.
At the large widths used for the tail improvement it does fit for
sufficiently large n. Existing finite-size fallbacks remain available.
The high-precision count endpoint is not changed.

## 7. Evidence boundary

The proof uses the existing all-input conjunction and exact dirty-tree
identities with a new nonuniform partition and resource sum. Its linear
Clifford count, live dirty width, and logarithmic T-depth are analytic
claims supported by the explicit bounds above, not by extrapolating
small circuits. The [five bounded checks](../../tests/test_nonuniform_dirty_indicator.py)
test unequal chunk order, the actual root conjugation, arbitrary
dirty-helper cancellation, small native edge phases, and exact integer
rounding/resource inequalities at 8,502 address lengths through one million.
The large cases use normalized arithmetic, not exponential-register simulation.
These checks do not supply an emitted large indicator,
a complete native frame compiler, or a new unrestricted depth lower
bound.
