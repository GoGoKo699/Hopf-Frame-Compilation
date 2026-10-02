# T-depth schedules for state-based Hopf QBP

[State-based task theorem](STATE_BASED_QBP_THEOREM.md) · [Depth routing](T_DEPTH_COMPILER.md) · [Parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md)

The state compiler uses a constant number of precision-source programs.
This chapter schedules their exact lookups and the actual coarse circuit
in T layers. One schedule prioritizes depth and may increase T-count;
the second retains a simultaneous count/depth bound with more dirty work.
The sources themselves are still scheduled serially.

## 1. Statement and physical model

Let $`N=2^n`$, $`n\ge1`$, and let the integer precision satisfy
$`P\ge\max\{6,n\}`$. Put $`B_0=P+n+7`$. Supply a real Hopf
state, or the consistently gauge-fixed complex state psi' of
[the complex compiler](COMPLEX_COARSE_COMPILER.md#1-target-gauge-and-workspace).
Use its actual recorded logical coarse unitary C, with

```math
\|C-W'\|\le\min\{1/64,1/(4\sqrt N)\}.
```

The two schedules below prepare either psi' alone or the coherent pair
with reference $`C|0^n\rangle`$, with the existing complete
initialized-isometry error at most $`2^{-P}`$. Two initialized compiler
flags suffice. The n initialized system qubits, the unchanged protocol
branch for a pair, and any observable work are separate resources.

**Theorem.** The following bounds hold for the same circuit within each
row, with arbitrary dirty input and external references:

| Schedule | Sufficient dirty width | T-count | Clifford count | T-depth |
|---|---|---|---|---|
| A: depth-optimized banks | $`b\ge2B_0`$ | $`O(NP)`$ | $`O(NP)`$ | $`O(NP/b+P+n^3)`$ |
| B: count-efficient banks and parallel indicators | $`b\ge16(B_0+\sqrt{NP})`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ | $`O(NP)`$ | $`O(P+n^3)`$ |

Appending the actual coarse inverse and $`H^{\otimes n}`$, as in
magnitude readout, preserves these orders. A fixed number of preparation
or readout uses also preserves them. The stronger reservation is a
literal sufficient condition, not a necessary workspace threshold.

Use the all-to-all logical model of
[the depth theorem](T_DEPTH_COMPILER.md#1-model-and-the-routing-primitive):
a T layer acts on distinct physical wires, and arbitrary Clifford
circuits may occur between layers. Their gates are counted. This is not
a total elementary-depth or elapsed-time bound. No clean address copies,
additional initialized signals, measurements, or intermediate resets are
used. Classical preprocessing remains separately charged as in the
[task theorem](STATE_BASED_QBP_THEOREM.md#5-construction-proof-dependencies-and-scope).

The coarse **logical** C and its recorded row coefficients are unchanged
by the schedules. Its physical word changes. Schedule A can cost
$`O(Nn)`$ T gates for C and does not retain the earlier count-optimized
state bound. Schedule B gives the simultaneous counts in its row.
Both physical implementations are exactly $`C\otimes I_b`$;
the fine state circuit instead returns its work within its norm error.

## 2. Exact queries and their live registers

For a Q-row, m-bit mask query, the output is an arbitrary occupied word.
Reserve it and the existing selectors/helpers before allocating banks.
Let mu be a power of two with $`1\le\mu\le Q`$. The completed
load/route/copy/unload echo from
[the depth lookup lemma](T_DEPTH_COMPILER.md#2-exact-lookup-with-depth-optimized-banks)
implements the original XOR query tensor identity on its added dirty
registers. It has

```math
D_{T,\rm query}=O(Q/\mu+\log\mu),\qquad
T_{\rm query}=O(Q/\mu+\mu m),\qquad G_{\rm query}=O(Qm).
```

Its word banks occupy $`\mu m`$ dirty wires. The high-address selectors
are part of the existing base reservation. The inverse is the actual
reversed gate word. These are literal all-input equalities, including
phases, arbitrary output contents, and entangled references.

Alternatively, [the parallel indicator](PARALLEL_DIRTY_LOOKUP.md#3-count-efficient-whole-word-queries)
replaces the high-address loader by an exactly equal loader. It uses
$`Q/\mu`$ dirty indicator outputs, disjoint from the banks and
output, and no additional scratch. It gives

```math
D_{T,\rm query}=O\!\left(1+\log(Q/\mu)+\log\mu\right)=O(1+\log Q)
```

with the same T- and Clifford-count bounds. Thus, for a base of size
$`B_0`$, its literal live-width condition is

```math
B_0+\mu m+Q/\mu\le b.
```

Choose mu as the largest power of two at most $`\sqrt{Q/m}`$ when
$`Q\ge m`$, and choose $`\mu=1`$ otherwise. The added width is
at most $`3\sqrt{Qm}`$ in the former case and $`2m`$ in the latter.
The count is $`O(\sqrt{Qm}+m)`$. These extra pools are reused
between completed queries; their widths are not summed over calls.

No initialized indicator is assumed. Conjugating a single X by the exact
dirty-bank router produces its one-hot XOR without initializing any work,
as proved in the linked chapter. The
indicator's classical-matrix CNOT layers can have substantial elementary
depth, which is not removed by this T-depth accounting.

## 3. Fine state preparation

First suppose $`n\ge6`$. The retained residual SU(2) table has three
real-angle factors, possibly conjugated by fixed target Cliffords. Set

```math
q=P+10,\qquad m=q+1=P+11,\qquad k=n-6,\qquad Q=2^k=N/64.
```

For a coherent pair, eight fixed address literals give 256 sectors; a
single-state table needs at most that many. The address, including the
protocol branch when present, is unchanged throughout each rotation.
The target and occupied compiler flags are not borrowed as synthesis work.

| Base register | Dirty width | Return |
|---|---|---|
| Precision core, also the mask-query output | $`P+11`$ | Included in the full-operator approximation |
| Query selectors | $`n-6`$ | Exact after each query |
| Predicate helper | 1 | Exact after each predicate |
| Borrowed synthesis signal | 1 | Included in the full-operator approximation |

These widths sum to exactly $`B_0`$, and $`m\le B_0`$. The X- and
Z-mask queries run sequentially and reuse all returned extra work. The
precision core, signal, and helper are never counted as available banks
while a query is live.

Replace each mask query on its full space by the exact query from Section
2 tensored with identity on the extras. This leaves the actual rotation
unitary unchanged on its old wires. In particular, it preserves the
exact commutation with X on the borrowed signal, actual inverse words,
and the inactive-predicate identity. Therefore the
[borrowed-signal error](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
remains $`43\,2^{-q}`$, including arbitrary signal and core inputs.
Accepted-block equality alone would not justify this replacement.

The repetitions are fixed: one scalar word has three source calls;
one pre-amplified rotation has three scalar words; its five-call
amplification uses 45 source calls and at most 20 addressed whole-word
queries. The fixed scalar masks are Pauli words. Three table factors,
three appearances of the state half-amplitude word, and at most 256
sectors multiply these bounds only by constants. Each source loader has
$`2q`$ fixed Pauli splitting rotations, so a source and its controlled
versions have a serial schedule of $`O(P)`$ T-depth and native count.
Controls act only on the central X in $`M=UX_0U^\dagger`$;
the loader U and its actual inverse remain unconditional. No individual
splitting rotation or T gate is promoted to a generically controlled gate.
The sector predicates have at most eight literals and constant depth.
No parallel source implementation is assumed.

The outer initial-state reflection has the existing $`O(n^2)`$
schedule using one arbitrary returned helper. The good reflection is
Clifford; the controlled Hadamard sweep costs $`O(n)`$ if scheduled
serially. All helpers are reused only after the preceding subroutine is
complete, and their exact identities hold even on leaked work.

For schedule A, let $`B=b-B_0\ge b/2\ge B_0`$ and choose the
largest power of two at most $`\min\{Q,B/m\}`$. Then
$`B_0+\mu m\le b`$ and

```math
D_{T,\rm fine}=O(P+n^2+NP/b),\qquad
T_{\rm fine},G_{\rm fine}=O(NP).
```

For schedule B use the count-efficient choice in Section 2. Since
$`Qm\le NP`$ and $`m\le B_0`$, its live width is bounded by
$`B_0+3\sqrt{NP}`$ or $`3B_0`$, each below
$`16(B_0+\sqrt{NP})`$. Every query has depth $`O(n+1)`$, giving

```math
D_{T,\rm fine}=O(P+n^2),\qquad
T_{\rm fine}=O(\sqrt{NP}+P),\qquad G_{\rm fine}=O(NP).
```

The reflection count is absorbed because $`n^2=O(\sqrt{NP})`$ for
$`P\ge n`$; its serial depth remains explicitly included.

The finite table-rounding error, exact query replacements, sector
direct sums, and actual-inverse state amplification retain
$`390\,2^{-q}\lt2^{-P}`$. Inactive sectors remain exactly identity;
no reset or projection is inserted between the three state-word uses.

For $`1\le n\le5`$, use the
[fixed-all-address source construction](RESIDUAL_TABLE_PREPROCESSING.md#6-small-systems-at-the-banked-reservation).
There are at most 128 sectors, no free address, and only
$`P+13\le2B_0`$ dirty wires for core, helper, and signal. Emit each
known one-row mask as a Pauli word. No selector, word bank, or indicator
is allocated. The fine depth and both gate counts are $`O(P)`$;
the same error proof holds. This separate case does not assume that its
live core fits a $`B_0`$ reservation or leave room for a needless bank.

## 4. Scheduling the exact coarse circuit

The real coarse frame has one native row word per addressed tree node;
the complex extension additionally has n full prefix-addressed phase
layers. Their coarse tolerance and geometric error allocation give
row lengths $`w_d=O(n)`$ at depth $`d=0,\ldots,n-1`$. Fix these
actual words before computing the fine residual coefficients. Each
determinant-one word has the literal linear-length rewrite in the fixed
reflection alphabet from
[the borrowed interpreter](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter).
For the real rows retain their recorded half-words
$`A_x=XQ_x^\dagger XQ_x`$, so that
$`XA_xX=A_x^\dagger`$ exactly and the actual row is $`A_x^2`$.
The real-sector echo below interprets these same half-words, not an
arbitrary newly synthesized approximation or a chosen square root.

For a reflection G, a one-bit dirty loader and router select the bank
value $`z\oplus f(x)`$. Apply G controlled by that value and, for a
real sector, the existing logical control beta. Undo routing and loading;
route the original bank again and apply the matched controlled G.
The target exponent is

```math
\beta(z\oplus f(x))\oplus\beta z=\beta f(x).
```

Omit beta for the full phase-prefix multiplexor. Since $`G^2=I`$,
the completed word returns every bank and selector exactly and applies
the prescribed reflection. Complete this echo before advancing to a
different reflection; different target reflections need not commute.
The serial or parallel high-address loader can be substituted here by
its full-space equality, even though the routed bank is used directly
rather than copied into an output word. Controlled reflections have
constant-size exact words with no scratch.

For the real tree and $`n\ge2`$, choose $`r=n-2`$ in the
[logical-sector echo](BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector).
At nonfinal depth d use all d prefix bits freely and choose the first
suffix bit as beta. At the final depth use one prefix bit as beta and
the other $`n-2`$ as free address, treating its two values separately.
Thus there is one sector per nonfinal depth and two at the final depth.
The remaining suffix predicate is toggled into beta using the logical
target as the arbitrary borrowed helper.
This directly specified allocation is legal with the external selectors
reserved below; it need not equal the count-optimized address cap in the
borrowed compiler's separate arbitrary-width schedule.

Crucially, the four predicate toggles in that echo surround entire
controlled table words. They are not repeated for every reflection
symbol. Their total cost is $`O(n^2)`$ per sector, hence $`O(n^3)`$
over the real coarse frame. Every loader is undone before such a toggle.
At depth d the free table has at most $`2^d`$ rows. A complex prefix
layer has exactly $`2^d`$ rows, acts on every suffix, and needs no
logical-sector predicate.

Reserve the same $`B_0`$ base while the fine program is idle. At most
$`n-1`$ external selectors are needed; they fit in this base. The real
predicate borrows the logical target, and controlled reflections require
no extra dirty helper. Additional one-bit banks and indicators are
disjoint from these selectors and logical controls. All returned work is
reused between symbols, sectors, and layers. For $`n=1`$, emit the
constant-length coarse row directly instead.

For schedule A, choose the largest power of two at most
$`\min\{Q,b-B_0\}`$ for each one-bit reflection table. A symbol
then has depth $`O(1+Q/b+\log Q)`$. There are $`O(w_d)`$ symbols
per layer and at most two sectors. Hence

```math
\begin{aligned}
D_T(C)&=O\!\left(\sum_{d=0}^{n-1}
 w_d[1+2^d/b+d]+n^3\right)
 =O(Nn/b+n^3),\\
T(C),G(C)&=O\!\left(\sum_{d=0}^{n-1}2^dw_d+n^3\right)
 =O(Nn).
\end{aligned}
```

The final equalities use $`w_d=O(n)`$, the geometric sum, and the
elementary absorption $`n^3=O(Nn)`$. Real and complex coarse words
obey these bounds; composing their two factors changes only constants.

For schedule B, each one-bit table uses mu as the largest power of two
at most $`\sqrt Q`$. Its extra width is at most $`3\sqrt Q`$;
$`Q\le N/2`$ and the common reservation covers it. Each symbol has
T-count $`O(\sqrt Q)`$, Clifford count $`O(Q)`$, and depth
$`O(d+1)`$. Consequently

```math
T(C)=O\!\left(\sum_d w_d2^{d/2}+n^3\right)
 =O(n\sqrt N),\qquad G(C)=O(Nn),
```

```math
D_T(C)=O\!\left(\sum_d w_d(d+1)+n^3\right)=O(n^3).
```

Here $`n^3=O(n\sqrt N)`$ absorbs the predicate count with an
absolute constant. No count or depth of an approximately returned
source is substituted for C. Every complete coarse reflection and
sector word implements the same exact logical action as before,
tensored with dirty identity. The actual inverse reverses this chosen
physical schedule and has identical counts and depth.

## 5. Composition and limits

Apply the scheduled exact C after the fine preparation. It maps the
computed residual to psi' exactly, so coarse approximation error does
not add to the final fine-isometry error. The coherent construction uses
the same common C on both branches and preserves their literal phase;
it does not control an arbitrary compiled C. Adding Sections 3 and 4
proves both rows of the theorem, using $`P\ge n`$ to absorb
$`Nn/b`$ into $`NP/b`$ and $`Nn`$ into $`NP`$.

The same addition proves the readout statement with actual
$`C^\dagger`$. Every new lookup register returns exactly; core and
borrowed-signal return remain part of the fine approximation norm.
The unchanged complete-input contracts therefore include arbitrary
reference correlations and conditional dirty states in repeated QBP
executions. The observable's depth and gate counts, and the number of
executions, remain separate charges in the
[complete serial task ledger](QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison).

Schedule A does not claim the count bound of schedule B at its smaller
width. Schedule B proves all three resource bounds on one physical
schedule, rather than choosing its count and depth from different
circuits. Neither schedule proves an optimal T-depth tradeoff, reduces
the serial precision-source depth below $`O(P)`$, or compiles the fine
prescribed complete frame. This is the state-based task's composition
of the retained source, exact dirty lookup, reflection-interpreter, and
state-amplification constructions.
