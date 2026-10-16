# R91: parallel bilinear word queries with returned dirty work

[Previous rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence) ·
[Blocked traversal and four-corner echo](../../docs/BLOCKED_BILINEAR_LOOKUP.md) ·
[Precision-terminal design](PRECISION_TERMINAL_PROMPT.md)

The rectangular query can process its output bits together without
repeating a change of coordinates on shared indicator registers. A
relative-phase word puts every T gate on a private dirty wire. An exact
echo removes both its relative phase and that wire's unknown initial
value. The resulting query removes the previous additive output-length
term in T-depth while preserving the same count envelope.

**Exact query theorem.** Let an a-bit address select one of
$`Q=2^a`$ prescribed q-bit words, with $`q\ge1`$. There is an
absolute constant C such that, for every extra dirty allocation
$`B\ge C(a+2)`$, a Clifford+T circuit implements

```math
|x,Z,d\rangle\longmapsto
|x,Z\mathbin\oplus f(x),d\rangle
```

with simultaneous bounds

```math
T=O\!\left(\sqrt{Qq}+\frac{Qq}{B}\right),\qquad
G=O(Qq),
```

```math
D_T=O\!\left(\log_2(a+2)+\frac{Qq}{B^2}\right).
```

The constants are absolute. B excludes the original address, q-bit
output, and other caller reservations. No initialized work is used;
every borrowed register returns exactly. The literal operator equality
includes arbitrary outputs, dirty states, superpositions, and external
references. Every inverse is the actual native inverse word. Original
address wires can be restricted to Clifford controls throughout.
Elementary Clifford gates are counted in G, although their depth between
T layers is unrestricted. For Q equal to one, emit its fixed Clifford
X word directly and use no work.

This is a query component, not a complete-frame theorem or a query
lower bound. In particular, for q larger than Q, a full indicator may
beat the displayed square-root count envelope. The proof claims no
external priority for relative-phase synthesis or the retained echoes.

## 1. A phase derivative on an arbitrary dirty wire

For three bits define the integer-valued expression

```math
p(d;y,z)=d-(y\mathbin\oplus d)-(z\mathbin\oplus d)
                  +(y\mathbin\oplus z\mathbin\oplus d).
```

Checking the four values of the pair y,z gives the exact identity

```math
p(d;y,z)=(4d-2)yz.                                      \qquad(1)
```

The diagonal word F with phase $`\exp(i\pi p/4)`$ has four
T layers. Place each of the four displayed parities on d by CNOTs,
apply its T or T-dagger phase, and undo those CNOTs. No other wire
receives a T gate. In operator notation,

```math
F(y,z;d)=\mathrm{CS}^{\dagger}(y,z)\,
                         \mathrm{CCZ}(y,z,d).           \qquad(2)
```

Equation (2) describes the exact relative phase; a separate CS gate is
not emitted. For many private wires $`d_{ij}`$ and shared controls
$`y_i,z_j`$, all copies of F run in four T layers. At each layer
the phase targets are distinct. Their parity CNOTs can be serialized
at that layer boundary and contribute only a constant number of
Clifford gates per private wire.

Let X be a separate unchanged register, and let L toggle each
$`d_{ij}`$ by a fixed binary linear form $`s_{ij}(X)`$.
L is an explicit Clifford CNOT word. Execute chronologically

```math
F^{\dagger},\quad L,\quad F,\quad L^{\dagger}.           \qquad(3)
```

Here F denotes the product over all pairs. L is its own inverse as an
operator, but the last word is still its actual inverse. Each d returns.
Its phase difference, modulo eight, is

```math
p(d\mathbin\oplus s;y,z)-p(d;y,z)
       =4s(1-2d)yz=4syz\pmod8.
```

Consequently the complete word (3) is exactly

```math
(-1)^{\sum_{i,j}y_i z_j s_{ij}(X)}
       \otimes I_d.                                    \qquad(4)
```

No initialized d, copied control, or discarded scalar phase is present.
The shared-control relative phase in (2) cancels because L does not
change y or z. This cancellation is needed; F alone is not CCZ.

## 2. Parallel bilinear outputs and their controlled version

Let Y and X contain H and J arbitrary bits, respectively. For each
output bit $`Z_j`$, prescribe a binary matrix
$`D_j\in\mathbb F_2^{H\times J}`$. Use

```math
s_{ij}(X)=(D_jX)_i.
```

The linear word L has at most qHJ CNOTs. Conjugate (3) by Hadamards on
every output Z. Equation (4) then gives the exact shear

```math
(Y,X,Z,d)\longmapsto
\left(Y,X,
  \left(Z_j\mathbin\oplus Y^{\mathsf T}D_jX\right)_{j=1}^q,d
\right).                                               \qquad(5)
```

Its T-count is at most 8qH, its T-depth is at most eight, its Clifford
count is $`O(qHJ)`$, and its extra dirty width is qH. The
matrices need not have a common rank decomposition: their row forms
are computed into separate arbitrary dirty wires. Y and X are never
changed by a coordinate transformation.

If only gH private wires are available, split the outputs into batches
of at most g. A completed batch preserves X and Y and returns its
private wires, so the same pool can serve the next batch.

For a separate arbitrary control h, reserve H further dirty bits
$`e_i`$. Let $`\mathcal B_e`$ be the diagonal version (4),
with e replacing y. It returns its d pool. Set

```math
E=\prod_{i=1}^{H}\mathrm{CCX}(h,y_i;e_i).
```

The [shared-control phase schedule](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
implements E in four T layers and $`O(H)`$ gates: the pairs
$`(y_i,e_i)`$ are disjoint. Execute chronologically

```math
\mathcal B_e,\quad E,\quad
\mathcal B_e^{\dagger},\quad E^{\dagger}.
```

The resulting phase is

```math
(-1)^{h\sum_{i,j}y_i z_j(D_jX)_i}.                       \qquad(6)
```

All e and d bits return. Hadamards on Z turn (6) into the controlled
version of (5). For one g-output batch its simultaneous resources are

```math
T=O(gH),\quad G=O(gHJ),\quad D_T=O(1),\quad
w=gH+H.                                                \qquad(7)
```

This is exact identity on the entire h-zero sector, including arbitrary
helpers. It is therefore a valid leaf for the retained dirty traversal.

## 3. A selected block and the complete four-corner query

Split the address into three fields, with $`Q=KHJ`$ and K,H,J
powers of two. The high field c selects one of K matrix families.
Use the retained
[two-pass dirty traversal](../../docs/BLOCKED_BILINEAR_LOOKUP.md#2-select-a-bilinear-block-with-dirty-traversal),
with (7) at each leaf. Process all q outputs in batches before leaving
that leaf. Each completed leaf is a controlled XOR shear preserving
X, Y, and the selectors. The completed shears commute. The two passes
therefore cancel the unknown dirty-selector contribution exactly.

The selected bilinear word has

```math
T=O(KqH),\qquad G=O(Qq),\qquad
D_T=O\!\left(K\lceil q/g\rceil\right).                  \qquad(8)
```

It uses gH+H private bilinear wires, p dirty selectors with
$`p=\log_2K`$, and a fixed number of reusable helpers. The same
d/e pool is reused after every completed batch. For K equal to one,
use the unguarded bilinear word; (8) remains a valid upper bound.

If the original address must enter only as a Clifford control,
implement the traversal's one- and two-literal selector updates with
the retained [read-only conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences).
Their arity is bounded, so this changes only absolute constants in
count and depth and adds a fixed reusable helper pool. No address wire
then receives a non-Clifford gate, even during the traversal.

Surround the entire selected word by the retained four-corner echo
with the dirty indicators
$`I_a:Y\mapsto Y\mathbin\oplus e_a`$ and
$`I_b:X\mapsto X\mathbin\oplus e_b`$:

```math
F_c,\ I_a,\ F_c^{\dagger},\ I_b,\quad
F_c,\ I_a^{\dagger},\ F_c^{\dagger},\ I_b^{\dagger}.
```

The four contributions to each output sum to
$`e_a^{\mathsf T}D_{c,j}e_b`$. Thus every indicator, selector,
and private helper returns, and the output receives exactly its table
word. Each indicator is called only twice, outside all output batches;
there is no extra indicator-depth factor $`\lceil q/g\rceil`$.

Enlarge the absolute
[linear indicator constant](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence)
$`c_I\ge1`$ once to dominate the fixed helper constants. A
sufficient simultaneous width is

```math
c_I(H+J)+gH+H+p+c_0.                                   \qquad(9)
```

The whole query has

```math
\begin{aligned}
T&=O(KqH+H+J),\qquad G=O(Qq),\\
D_T&=O\!\left(K\lceil q/g\rceil+\log_2(a+2)\right).
\end{aligned}                                         \qquad(10)
```

The Clifford count includes every row-form CNOT, output Hadamard,
traversal update, and actual inverse. The Hadamards contribute at most
$`O(Kq)\le O(Qq)`$ gates. Literal basis identities prove the full
operator equality and hence its external-reference version.

## 4. An allocation at every sufficient dirty width

Let pow2floor denote the largest power of two not exceeding its
argument, and set

```math
t=\mathrm{pow2floor}\,\frac{B}{64c_I}.
```

Choose the absolute minimum-width constant C so that t is at least one
and $`a+c_0\le B/4`$. Then $`t=\Theta(B)`$ with absolute
constants.

First, if $`q\ge Q`$ and $`Q\le t`$, use a full Q-bit dirty
indicator Y. Let C be the Clifford word XORing its table-weighted
linear forms into Z. The chronological word $`C,I,C^{\dagger},I^{\dagger}`$
returns Y and adds exactly the selected row. It costs $`O(Q)`$
T gates, $`O(Qq)`$ Cliffords, and logarithmic address T-depth,
using at most $`c_IQ`$ dirty wires. Here $`Q\le\sqrt{Qq}`$,
so it meets the stated count envelope.

In every remaining case choose

```math
\begin{aligned}
J&=\mathrm{pow2floor}\,
               \min\{Q,\sqrt{Qq},t\},\qquad g=\min\{q,t\},\\
H&=\mathrm{pow2floor}\,\min\{Q/J,t/g\},\qquad K=Q/(HJ).
\end{aligned}                                         \qquad(11)
```

H and J are powers of two dividing Q, $`HJ\le Q`$, and
$`gH\le t`$. Both H and J are at most t. Therefore (9) fits B,
including both indicator pools at once. The three cases below exhaust
(11); q need not be a power of two.

| Case | Allocation consequences | Bilinear depth contribution |
|---|---|---|
| $`t\ge\sqrt{Qq}`$ | The remaining branch has $`q\lt Q`$; $`g=q`$, $`J=\Theta(\sqrt{Qq})`$, $`H=\Theta(\sqrt{Q/q})`$, $`K=O(1)`$ | $`O(1)`$ |
| $`t\lt\sqrt{Qq}`$ and $`q\le t`$ | $`J=t`$, $`H=\Theta(t/q)`$, $`K=O(Qq/t^2)`$ | $`O(Qq/t^2)`$ |
| $`q\gt t`$ and $`Q\gt t`$ | $`g=t`$, $`H=1`$, $`J=t`$, $`K=Q/t`$ | $`O(Qq/t^2)`$ |

For the first case, power-of-two rounding gives
$`J\gt\sqrt{Qq}/2`$ and $`H\ge\sqrt{Q/q}/2`$.
The bilinear count is $`KqH=Qq/J=O(\sqrt{Qq})`$;
the indicators have the same count bound.

In the second case the direct-indicator branch has already removed
q at least Q. Thus $`q\lt Q`$ and $`t\lt Q`$. The inequality
$`t^2\lt Qq`$ forces the width term t/q to determine H in (11).
Rounding gives $`H\gt t/(2q)`$. The count is $`Qq/t`$, and
the indicator cost $`O(t)`$ is absorbed by it.

In the third case $`\lceil q/t\rceil\le2q/t`$. The count is
again $`Qq/t`$, with the same absorption of indicators. If
$`q\gt t`$ but $`Q\le t`$, the direct-indicator branch would
already have applied, so no fourth case is missing.

Substitution into (10), with $`t=\Theta(B)`$, proves the theorem.
All matrix circuits are charged separately, while only their returned
dirty pools are reused. A larger physical B can simply leave unused
wires untouched. No finite cutoff depends on the table or its precision.

## 5. What this removes from the terminal problem

For $`B=\Theta(\sqrt{Qq})`$ with a sufficiently large absolute
constant, a complete q-bit word query now has logarithmic address
T-depth. The previous serial q term came from processing different
matrix coordinate changes in sequence; (3) removes that requirement.

The output is still an arbitrary word XORed with f(x). This theorem
does not make it a clean program. The query alone does not remove the $`2^q`$ catalogue or the local
synthesis cost in [R90](DIRTY_CATALOGUE_TERMINAL.md). Section 6 supplies
a separate faster slot router at a larger dirty allocation. Composing
the reader through
noncommuting target gates with changing local selectors requires its
own complete returned-work identity.

The bounded [native fixtures](../../tests/test_parallel_bilinear_word_query.py)
check literal phase on every basis input of the relative-phase and
multigrid words, controlled return, and a complete small four-corner
query. Exact integer checks also probe power-of-two allocation
boundaries. They validate these interfaces; the proof above supplies
the arbitrary-size and simultaneous-resource statements.

## 6. A faster selected-slot router with a larger dirty allocation

The same grid permits a separate tradeoff for the catalogue router.
Let $`V=2^q`$ and let A and E be two disjoint arbitrary arrays,
each consisting of V slots of m bits, with $`m\ge1`$. A separate
q-bit register r specifies the XOR permutation

```math
(P_rA)_{t,\beta}=A_{t\mathbin\oplus r,\beta},
\qquad P_r^2=I,\qquad P_rP_s=P_{r\mathbin\oplus s}.
\qquad(12)
```

There is a native exact shear

```math
(r,A,E)\longmapsto(r,A,E\mathbin\oplus P_rA)
\qquad(13)
```

with $`O(mV^2)`$ T-count, Clifford count and additional dirty
width, and T-depth $`O(\log_2(q+2))`$. The additional width
excludes r and both live arrays. Every additional helper returns
exactly; neither A nor E is initialized.

To prove (13), reserve a V-bit arbitrary dirty indicator Y and its
private helpers, and one private grid bit $`d_{u,t,\beta}`$ for
every triple. Equation (5) gives the complete XOR shear

```math
E_{t,\beta}\longmapsto E_{t,\beta}
 \mathbin\oplus\bigoplus_{u\in\mathbb F_2^q}
                   Y_uA_{t\mathbin\oplus u,\beta}.
\qquad(14)
```

Call it $`\mathcal B`$. Each row form has exactly one input A
bit, so the linear map into the private grid costs mV squared CNOTs.
The dense matrix estimate from Section 2 is unnecessary here. The grid
has mV squared private wires, eight T layers and $`O(mV^2)`$
gates. A remains unchanged throughout each completed word.

Let $`I_r`$ toggle Y by $`e_r`$. Execute chronologically

```math
\mathcal B,\quad I_r,\quad
\mathcal B^{\dagger},\quad I_r^{\dagger}.
```

Linearity in Y cancels its arbitrary initial value and leaves exactly
(13). This needs only one indicator echo: A is the quantum data input,
not another indicator whose unknown offset should be removed. The
indicator has logarithmic address depth and linear size in V. Literal
phase, arbitrary inputs and references are preserved. For q zero,
(13) is simply a fixed CNOT word.

### Pair the translations, then return the auxiliary array in the full word

Write $`U_P:E\mapsto E\mathbin\oplus P_rA`$ and
$`L_P:A\mapsto A\mathbin\oplus P_rE`$. Their completed words
return the same reusable grid and indicator pool. Apply chronologically

```math
U_P,\quad L_P,\quad U_P,\quad\mathrm{SWAP}(A,E).
```

Using $`P_r^2=I`$, the four successive array values are

```math
\begin{aligned}
(A,E)&\longmapsto(A,E+P_rA)\\
     &\longmapsto(P_rE,E+P_rA)\\
     &\longmapsto(P_rE,P_rA)\\
     &\longmapsto(P_rA,P_rE).
\end{aligned}                                         \qquad(15)
```

Here additions are over $`\mathbb F_2`$. The final fixed array
swap is an explicit Clifford word of 3mV CNOTs. Thus the completed
paired translation is

```math
\widehat T_r=P_r(A)\otimes P_r(E),
```

with the same $`O(mV^2)`$ count and additional-width bounds,
and T-depth $`O(\log_2(q+2))`$. E is permuted at this boundary;
it is not yet a returned helper. All grid and indicator work has returned.

Now let P be an immutable external prefix and let the exact query
$`Q_f`$ toggle a separate arbitrary dirty program R by f(P).
It excludes the logical m-bit target X and both arrays. The chronological
word

```math
Q_f,\quad\widehat T_R,\quad
Q_f^{\dagger},\quad\widehat T_R^{\dagger}
\qquad(16)
```

is exactly $`\widehat T_f`$, with R and all query/router helpers
returned. This follows from the XOR group law on both arrays. No clean
program is produced or assumed.

Let $`S_0`$ be the fixed swap between X and slot zero of A. Then

```math
\widehat T_f^{\dagger}S_0\widehat T_f
       =S_f(X,A)\otimes I_E.                           \qquad(17)
```

The E permutation commutes with $`S_0`$ and cancels in this
complete conjugation. Hence (17) is the exact selected swap, including
return of E, R, all grid and indicator work, and all query work. This
holds even when those registers are entangled. Every dagger denotes
the actual inverse circuit word. The paired translation alone must not
be substituted for a returned-work single-array translation elsewhere.

### Consequence for the local catalogue ledger

Use (17) for every selected swap in R90's two-array commutator. Keep
those two original arrays disjoint. A third auxiliary array E and one
grid/indicator pool can be reused between completed selected swaps,
since each such word returns them exactly.

Suppose the immutable address has a bits and $`Q=2^a`$ rows,
the program has $`q\ge1`$ bits, and each fixed catalogue factor
has T-count, Clifford count and T-depth at most S. Reserve
$`B\ge C(a+2)`$ dirty wires exclusively for the Section 4
query, outside all live programs, arrays and router work. A sufficient
whole-catalogue ledger is then

```math
\begin{aligned}
T&=O\!\left(\sqrt{Qq}+\frac{Qq}{B}+VS+mV^2\right),\\
G&=O(Qq+VS+mV^2),\\
D_T&=O\!\left(\log_2(a+2)+\frac{Qq}{B^2}
                         +S+\log_2(q+2)\right).
\end{aligned}                                         \qquad(18)
```

A sufficient simultaneous dirty reservation, excluding the logical
address and X, is

```math
q+3mV+C_{\rm router}mV^2+B,                            \qquad(19)
```

for one absolute router constant. Every query, factor, and actual
inverse is charged. The fixed factors still run in parallel on distinct
slots. Their exact dirty-array cancellation and the direct-sum maximum
error estimate are unchanged. A zero-bit program has one fixed entry;
its query and routing can be omitted.

This corollary reduces the router's q term to logarithmic q depth by
paying quadratic catalogue volume. It is useful only when that larger
simultaneous allocation and gate budget fit. It does not alter R93's
retained group synthesis or initialized-label ledgers, and supplies no
additional complete-frame depth bound.
