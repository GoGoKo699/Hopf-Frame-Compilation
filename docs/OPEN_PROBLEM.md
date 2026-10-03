# Research status and the constant-clean endpoint

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [One-clean compiler](ONE_CLEAN_COMPILER.md) · [Grouped refinement](CONDITIONAL_SUFFIX_COMPILER.md)

This page records the complete-frame frontier and the limits of explored
routes. The separate [state-based QBP theorem](STATE_BASED_QBP_THEOREM.md)
consolidates the completed changed-decoder result, its real/complex
contracts, resource bounds, computational assumptions, and proof map.
Its two-flag state preparation does not resolve the prescribed
complete-frame endpoint.

The [publication scope](../manuscript/PUBLICATION_SCOPE.md) retains the
established frame results. The [verification map](VERIFICATION.md)
distinguishes analytic proofs from implemented and finite evidence.

## Revision checkpoint and selected next test

This checkpoint separates completed component results from the remaining
complete-frame questions. Write $`N=2^n`$. The uniform-precision and
low-precision bounds below remain established and unchanged; these gaps
are separate from the completed Hopf-QBP contract.

| Remaining gap | Current boundary | Selected treatment |
|---|---|---|
| Large-width T-depth | At fixed accuracy, two clean flags and sufficient $`b=\Theta(\sqrt N)`$, $`\Omega(1)\le D_T^\star\le O(n)`$, with $`T=O(\sqrt N)`$ and $`G=O(N)`$ | Three overhead refinements are complete; audit the remaining program-query and logical-stage costs |
| High-precision complete-frame count endpoint | At $`a=2,L=N,b=N+n+7,n\ge3`$, $`\Omega(N)\le T^\star\le O(N\ell_*(n))`$, where $`\ell_*(n)=1+\log_2^*(n+2)`$ | Park until a new complete native identity supplies its symbolic precision, workspace, and work-return ledger |

### Completed components and the current depth ledger

The [nonuniform indicator](NONUNIFORM_DIRTY_INDICATOR.md) implements the
exact dirty/reference-safe transformation

```math
|x,Y,W\rangle\longmapsto|x,Y\oplus e_x,W\rangle,
\qquad T,G,w=O(2^s),\qquad D_T=O(\log_2(s+2))
```

for an s-bit address. Its rectangular-query consequence removes the old
linear indicator allowance from the low-precision late tail at sufficient
width. The original literal full-frame reservation is retained through
its established fallbacks; no matching query-depth lower bound is claimed.

The [protected-source theorem](PROTECTED_UNARY_SOURCE.md) retains one
unary source across the entire early segment. H records the original
zero predicate of a fixed terminal logical bank. The preparation U is
unconditional; the actual $`U^\dagger`$ is used at the global exit.
Inactive inputs cancel exactly with h initially zero. On the active
sector, one global $`2\delta`$ error bound includes all final source
and H leakage. Setting $`\delta=\eta/8`$ gives a single
$`O(L+\log(n+2))`$ source-boundary depth allowance.

The [windowed-predicate theorem](WINDOWED_GROUP_PREDICATES.md) adds
$`2J+1`$ disjoint conditional-zero cache bits to that bank, where
$`J=\lceil\log_2(n+2)\rceil`$. One long outer-suffix predicate,
short block-zero bits, and a suffix-product chain supply the activity
flags for up to J consecutive groups. Block caches are erased before
their logical targets change; chain bits are consumed while their later
controls remain valid. The full cache returns exactly, including arbitrary
inactive-cache and dirty/reference inputs at $`Hh=00`$. Sequential short
predicate calls reuse the original two dirty helpers. Five bounded checks
cover unequal/partial windows, retained sources, arbitrary inactive caches,
native flag phases, and incorrect cleanup or guard choices.

Its aggregate activity depth in the uniform low-precision regime is

```math
D_{T,\rm activity}
=O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}+\log(n+2)\right)=o(n).
```

At fixed accuracy and sufficient square-root dirty width, the resulting
same-circuit ledger is:

| Contribution | Current total depth allowance |
|---|---|
| Logical stages and incremental selectors | $`O(n)`$ |
| Program prefetch and unload | $`O(n)`$ |
| Windowed activity predicates | $`O(n\log\log(n+2)/\log(n+2)+\log(n+2))`$ |
| Protected-source and initial-bank-predicate boundaries | $`O(\log(n+2))`$ |
| Late queries, sources, and predicates | $`o(n)`$ |

The two linear rows are upper bounds for the current construction, not
lower bounds. No improved complete-frame asymptotic order or high-precision
endpoint follows from the three completed overhead refinements alone.

### Query audit: what can and cannot be retained

The [shared-prefix query audit](SHARED_PREFIX_QUERY_AUDIT.md) now tests
this interface directly. With a retained dirty indicator
$`Y=d\oplus e_x`$, completing only the short-address echo reads the
selected table bit plus $`d^{\mathsf T}D e_t`$. Reusing a correction
computed at an earlier local address leaves

```math
d^{\mathsf T}D(e_{t_{\rm old}}\oplus e_{t_{\rm new}}).
```

It cancels for every dirty d exactly when the two table columns coincide.
The chapter gives a native two-group counterexample and an exact positive
example whose extra baseline bit must start clean. The latter is a
reservation demonstration, not a query-depth improvement.

A broader boundary applies to one specific reader architecture. Its
encoder uses the unchanged logical prefix x only as a read-only control
on a cache of c initialized and b arbitrary dirty bits. Subsequent exact
row readers cannot access x or an uncounted prefix-bearing register;
any final x-dependent decoder acts only on the cache. If the table has R
distinct row signatures, orthogonal encoded supports give

```math
R2^b\le2^{b+c},\qquad c\ge\lceil\log_2 R\rceil.
```

This is an elementary exact-programming dimension argument, with the
Nielsen–Chuang precedent identified in the proof and
[related-work comparison](RELATED_WORK.md#26-exact-shared-prefix-cache-capacity-3-october-2026).
It is not a T-depth lower bound. Conditional logical zeros and every
other initialized wire touched by the prefix encoder count among c.

For G consecutive target bits there are $`2^G-1`$ independent local
angle nodes per prefix. Two allowed angles per node already give legal
tables with $`R=2^{\min(d,2^G-1)}`$. Thus the necessary cache capacity is

```math
c\ge\min\{d,2^G-1\}.
```

Generously granting all k remaining logical bits and both flags to this
cache gives $`c\le k+2`$. In the latter portion where $`d=n-k\gt k+2`$,

```math
G\le\log_2(k+3).
```

Each current early group has height of order log n from below in the
uniform low-precision regime. This exact reusable-reader architecture
therefore cannot join a growing number of current groups into a window
there. The conclusion permits arbitrary cache encoding within its stated
interface; it does not assume q-bit one-hot rows. When $`c\ge d`$, a
clean copy of x is possible and this window-height consequence does not
follow. Restricted tables must use their actual R.

The earlier direct one-hot prefetch restriction remains a simpler special
case of storage accounting: its program uses $`q(2^G-1)`$ clean logical
bits for unary phase modulus q. The new dimension argument neither forces
that representation nor claims that every fused frame circuit must expose
an exact reader for every angle node.

### Charged correction: an exact identity with a remaining query

The [charged two-group word](SHARED_PREFIX_QUERY_AUDIT.md#6-a-charged-two-group-program-refresh)
completes the selected prefix-access audit. Both groups use a common padded
program bank. One prefix indicator remains live while the first body runs;
an exact table-difference query changes the program to the second body's
word at the current, changed address. The last actual inverse returns the
program and indicator. The identity holds for arbitrary dirty program and
source inputs, with the body's required one-hot promise stated separately
at its use points. A transition of the existing activity flag is included
without storing a third clean flag.

For E groups, this algebra leaves one outer indicator pair and E minus one
correction queries. Every correction and its temporary work remain charged.
The generic correction family contains ordinary fresh table queries: fix
the old legal one-hot word to a constant and vary the new row freely. One
fixed Clifford XOR recovers the new query. A constant first angle of
$`\pi/2`$ can change the next address; its real Hopf rotation
is $`XZ`$, including the sign on its one input. It imposes no relation on
the later angle rows.

There is a useful restricted positive example. For an explicitly affine
parity u of the prefix and a current local bit t, the correction
$`\Delta(x,t)=u(x)t`$ uses two literal Toffolis and Clifford parity
computations, returning one arbitrary dirty helper. But rank one alone
is insufficient: $`f(x)t`$ with arbitrary f contains a fresh query at
$`t=1`$. Independent Hopf parameters do not supply cheap prefix factors.
These are exact identities and reductions, not unrestricted depth lower
bounds. The complete-frame frontier and both linear allowances are unchanged.

**Stop rule.** Do not continue the general table-difference route without
either a better generic query schedule or a factorization proved cheap for
every admissible row family. Recounting fewer explicit indicator symbols,
assuming related rows, or moving the readout into a conjugated basis does
not pay for the missing correction. The cache-capacity and whole-window
offset-echo failures remain in the same proof chapter.

### Retained-source fusion: the exact small-block result

The [retained-source continuation](UNARY_PHASE_GRADIENT.md#9-retained-source-fusion-without-a-larger-modulus)
now lifts the existing ideal magic-basis benchmark to the complete source
operator. Its four-mode factors have opposite determinant monomials.
Their product is exact over the original Laurent ring with $`S^q=I`$;
no half-angle source enlargement is required. This is a shared-source
product, not two independent phase states. The two joint exponent tables
use only the original labels and require two selected-shift rounds.

At eight modes, jointly diagonalizing the Bell-controlled root leaves
only integer exponents zero or plus/minus the root label. One original
source shift implements the whole root; its four commuting Pauli terms
do not require four quarter-angle sources. Together with the two child
rounds, this recovers three rounds. The retained compiler still supplies
the height-linear schedule; diagonalization alone does not supply a new
transformed-selector construction.

There is also a complete Bell-stabilizer extraction. For each child
$`C_s=A_s\otimes B_s`$ over the shared source ring, put

```math
V_s=A_sB_s^{\mathsf T},\qquad
K_s=B_s^{-\mathsf T}\otimes B_s,\qquad
C_s=(V_s\otimes I)K_s.
```

The transpose affects only logical matrix indices. Each $`K_s`$ fixes
the Bell state for every source input, so the controlled K commutes with
the parent root U even for unequal children. Consequently
$`CU=VUK`$ and $`CUC^\dagger=VUV^\dagger`$. The full word
still needs K. Its direct charged schedule uses six rounds versus three
for the baseline, and the conjugated-root schedule uses seven versus
five. Derived one-hot program rows are explicitly charged. These upper
schedules establish an exact interface, not a native depth improvement.
The retained transport has noncommuting matrix coefficients; it is not
another scalar two-dimensional factor to which the same formula applies.

### A shallow completion and the remaining transport

The [growing completion schedule](UNARY_PHASE_GRADIENT.md#10-two-source-shifts-for-the-stabilizer-completion)
uses the stabilizer extraction across an even-height group. Its completion
factors act on mutually orthogonal sectors: the last nonzero target pair
and the earlier computational prefix identify the sector. Their two
paired-axis words preserve these sectors on arbitrary source inputs.
Cache all sector labels before applying the logical basis changes, use
the cached labels in both selected shifts, and erase them after restoring
the original basis. No transformed prefix is reread as an old address.

For even g in the existing range, this gives

```math
D_{T,\rm completion}=O(\log(g+2)),\qquad
T,w=O(q2^g+R),\qquad G=O(q2^g+qR),\qquad R=3^\ell.
```

The ledger excludes the charged program query pair, outer activity
predicates, and protected-source preparation/return. The program includes
derived difference rows. All conditional-zero cache and selector work is
reserved; on the inactive sector the full word cancels on arbitrary work.
The result is an exact component compiler on every logical/source column.

The complete factorization still has an ordered transport product with
three source-shift factors per target pair. Its available depth remains
O(g). A faster state preparation sharing its first column cannot replace
that product without pricing the change to every complementary column.

**Decision.** The next unresolved interface is this noncommuting transport;
no faster native rule is selected. Before another fixture pass, require
an explicit rule with a depth/count/width recurrence and complete cleanup,
or a better generic query schedule. Fixed magic-basis recursion or source
phase rearrangement alone does not supply that rule. A small next matrix
is not by itself a reason to open another construction pass.

The complete-frame frontier and both linear allowances are unchanged.
An exact rewrite must hold on all source characters; a different unitary
body that agrees only on the ideal source remains eligible with a proved
global initialized-isometry and inactive/work-return contract.

## Established frontier

Write $`N=2^n`$, $`q=n+a+b`$,
$`h=1+\lceil\log_2(L+n+2)\rceil`$, and
$`\ell_*(n)=1+\log_2^*(n+2)`$. The exact model has arbitrary one-qubit
gates and CNOTs; the approximate model has coherent Clifford+T gates.
The clean budget in the exact model is m; in the approximate model a and b
count clean and arbitrary dirty qubits. The precision parameter L is defined
below. Lower bounds are worst-case over the stated frame family. For the
hybrid depth bounds, write

```math
\chi(t)=\log_2(t+2).
```

| Question | Retained result | Status and proof |
|---|---|---|
| Exact size and depth versus clean workspace | Size $`\Theta(N)`$ and depth $`\Theta(n+N/(n+m))`$ for every $`m\ge0`$; CNOT count $`\Theta(N)`$ for $`n\ge2`$, zero for $`n=1`$ | Matching for the complete real and phase-dressed complex magnitude frames; [exact theorem](COMPILER_THEOREM.md) |
| T-count with sufficient clean workspace | $`T^\star=\Theta(\sqrt{NL}+L+NL/q)`$, $`G=O(NL)`$, when $`a\ge C(n+h)`$ for sufficiently large fixed C | Matching for those same families; the clean reservation is sufficient, not proved necessary; [fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md) |
| Real-frame T-count at arbitrary workspace budgets | $`T=O(NL/q+L\sqrt N)`$, $`G=O(NL)`$, for every $`a,b\ge0`$ | Splicing with the sufficient-clean theorem gives the matching frontier above for every a when $`h+b\le c\sqrt N`$, for fixed $`c>0`$; [borrowed-workspace proof](BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements) |
| Zero-clean layerwise real-frame T-count | $`T=O(N+nL)`$, $`G=O(NL)`$, at $`a=0`$, $`b\ge L+n+7`$ | Full-operator approximation using a borrowed amplification signal; [signal-symmetry corollary](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) |
| Grouped one-clean real-frame T-count | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge L+n+7`$ | Precision-uniform grouped bound without additional word banks; [one-clean extension](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean T-count with additional dirty banks | $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$, at $`a=1`$, $`b\ge2(L+n+7)`$ | Matches the lower bound if $`L\ell_*(n)^2\le N`$ or $`b\le N/\ell_*(n)`$; these are sufficient regimes; [banked one-clean proof](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) |
| One-clean phase-dressed complex magnitude frame | The same grouped and banked T-counts, at $`b\ge L+n+8`$ and $`b\ge2(L+n+8)`$, respectively | Compose the real compiler and literal phase diagonal in the same workspace; [composition corollary](ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames) |
| T-depth with additional dirty banks | $`D_T=O(NL/b+\min\{nL+n^2,L\ell_*(n)+n^3\})`$ at $`a=2`$, $`b\ge2(L+n+7)`$, with $`T,G=O(NL)`$ | Choose between the layerwise and grouped [schedules](T_DEPTH_COMPILER.md); real frames; optimizing depth may increase T-count; no matching frontier established |
| Simultaneous T-count and T-depth | $`T=O(\sqrt{NL}+L\ell_*(n))`$, $`D_T=O(\min\{nL+n^2,L\ell_*(n)+n^3\})`$, $`G=O(NL)`$, at $`a=2`$, $`b\ge C(L+n+7+\sqrt{NL})`$ | Same real-frame circuit, for sufficiently large fixed C; [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md); T-depth optimality remains open |
| Fixed-accuracy count and depth versus width | $`T=O(\sqrt N+N/b)`$, $`D_T=O(N/b^2+n)`$, $`G=O(N)`$, at $`a=2`$, fixed L, $`b\ge17(L+n+7)`$ | Same complete real-frame circuit; count is optimal in order, and depth is matching through $`b\le\sqrt{N/n}`$; [blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md) |
| Uniform-precision count and depth | $`T=O(\sqrt{NL}+NL/b+nL)`$, $`D_T=O(NL/b^2+nL)`$, $`G=O(NL)`$, at $`a=2`$, $`L\ge6`$, $`b\ge17(L+n+7)`$ | Same complete real-frame circuit, with absolute constants; both are matching through $`b\le\sqrt{N/n}`$ when eligible; [uniform-precision theorem](UNIFORM_PRECISION_DEPTH.md) |
| Low-precision count and depth | $`T=O(\sqrt{NL}+NL/b)`$, $`D_T=O(NL/b^2+n)`$, $`G=O(NL)`$, at $`a=2`$, $`6\le L\le\log_2(n+2)/16`$, $`b\ge17(L+n+7)`$ | Absolute constants; count is optimal in order, and both resources match through $`b\le\sqrt{NL/n}`$ when the interval is nonempty; [low-precision theorem](UNIFORM_PRECISION_DEPTH.md) |
| Fixed-accuracy large-width depth | $`T=O_\eta(\sqrt N)`$, $`G=O_\eta(N)`$, $`D_T=O_\eta(n)`$, at $`a=2`$, $`b\ge C_\eta\sqrt N`$ | Same complete real-frame circuit with charged unary source preparation/return; [unary theorem](UNARY_PHASE_GRADIENT.md#7-complete-frame-theorem-at-fixed-accuracy); T-count is optimal in order, depth lower bound remains $`\Omega(1)`$ |

Take the best applicable construction. For fixed L, $`a=2`$ and
$`b=L+n+7=\Theta(n)`$, the arbitrary-budget matching splice gives
$`T^\star=\Theta(N/n)`$, sharper than the grouped estimate. Extra clean
qubits may be left unused; the depth schedules retain their separately
proved two-clean allocation.

The uniform-precision count is optimal at every eligible width under the
sufficient condition $`L\le N/n^2`$, which absorbs its $`nL`$ term
into $`\sqrt{NL}`$. No count-optimality claim is made for all L.
The uniform and low-precision bounds include arbitrary dirty/reference
inputs and every returned work register on the same circuit.

All frame constructions preserve the prescribed completion and the
[fixed-parameter QBP error contract](QBP_APPROXIMATION.md). They do not
differentiate discrete synthesis. Complex leaf-phase derivatives remain a
separate stream. Literal diagonals and one-target U(2) multiplexors retain
their independent matched banked frontiers; earlier weaker frame bounds
are dependencies or fallbacks, not the current frontier.

## The count and depth gaps are different

The exact ancilla-depth theorem is matching in its elementary-gate model;
T-count is matching under its stated clean reservation. The two-clean
T-depth tradeoff is now matching in explicit precision/workspace ranges,
including fixed and inverse-polynomial accuracy in N. Large-width depth
and precision regimes outside those ranges remain open. Complete-frame safety and
QBP substitution are already established; the gaps concern resources.

For the exactly two-clean banked regime, put $`B_0=L+n+7`$ and assume
$`b\ge2B_0`$. Then $`q=n+2+b=\Theta(b)`$. The
[inherited depth lower bound](T_DEPTH_COMPILER.md#4-lower-bounds-and-the-remaining-depth-gap)
simplifies to

```math
D_T^\star=\Omega(1+NL/b^2).
```

This follows from $`L/b=O(1)`$ and
$`\sqrt{NL}/b\le1+NL/b^2`$, not a new lower-bound argument. The table
uses $`a=2`$ throughout and respects each sufficient allocation threshold.

| Regime | Depth lower bound | Available depth upper bound | Remaining issue |
|---|---|---|---|
| Fixed L, $`b=\Theta(n)`$ and $`b\ge17B_0`$ | $`\Omega(N/n^2)`$ | $`O(N/n^2)`$ with $`T=\Theta(N/n)`$ | Matching count and depth in one circuit; [amortized schedule](AMORTIZED_DIRTY_LOOKUP.md) |
| Fixed L, $`17B_0\le b\le\sqrt{N/n}`$ | $`\Omega(N/b^2)`$ | $`O(N/b^2)`$ with optimal-order count | Matching throughout this interval when nonempty |
| Variable L, $`17B_0\le b\le\sqrt{N/n}`$ | $`\Omega(NL/b^2)`$ | $`O(NL/b^2)`$ with $`T=\Theta(NL/b)`$ | Uniform-precision matching interval, when nonempty |
| $`6\le L\le\log_2(n+2)/16`$, $`17B_0\le b\le\sqrt{NL/n}`$ | $`\Omega(NL/b^2)`$ | $`O(NL/b^2)`$ with $`T=\Theta(NL/b)`$ | Larger low-precision matching interval, when nonempty |
| $`L=\Theta(n)`$, sufficient $`b=\Theta(n)`$ | $`\Omega(N/n)`$ | $`O(N/n)`$ with $`T=\Theta(N)`$ | Matching at inverse-polynomial error in N for sufficiently large n |
| Fixed L, $`2B_0\le b\lt17B_0`$ | $`\Omega(N/n^2)`$ | $`O(N/n)`$ | Earlier schedule remains the proved fallback at this literal reservation |
| Fixed L, sufficiently large $`b=\Theta(\sqrt N)`$ | $`\Omega(1)`$ | $`O(n)`$ with $`T=O(\sqrt N)`$ | [Unary phase-source groups](UNARY_PHASE_GRADIENT.md#7-complete-frame-theorem-at-fixed-accuracy); depth lower bound remains unmatched |
| Fixed L, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(n)`$ with $`T=O(\sqrt N)`$ | Extra width is not needed by this schedule; depth optimality remains open |
| $`L=N`$, $`b=\Theta(N)`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ | Serial precision cost remains |
| Selected endpoint $`L=N,b=B_0`$ | $`\Omega(1)`$ | $`O(N\ell_*(n))`$ from $`D_T\le T`$ | The larger-bank depth theorem does not apply |

The amortized schedule removes the earlier per-batch logarithm from
both indicator routing and chunk selection. The earlier hybrid replaced
its additive $`n^2`$ term by $`n\chi(n)`$ at every eligible width
and precision. It remains a dependency of the uniform theorem.
The [capped source precision](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
reduces source depth to $`O(n\log(n+1))`$ at fixed L without changing
the count or workspace orders. The quadratic contribution remaining in
this schedule comes from query routing; it is not an unavoidable source cost.
The exact dirty-counter indicator, [masked-sum refinement](DIRTY_SUM_COMPRESSION.md),
and bilinear-query hybrid improved
the square-root-width upper bound to $`O(n\chi(n))`$. The
[grouped-program refinement](GROUPED_PROGRAM_PREFETCH.md#8-complete-frame-theorem-at-fixed-accuracy)
gives $`O(n\log\log(n+2))`$ at fixed accuracy, with returned
work and the same optimal-order T-count. It combines cached conditional
programs with a summable chunked-indicator budget for the late layers.
The [unary phase-source refinement](UNARY_PHASE_GRADIENT.md#7-complete-frame-theorem-at-fixed-accuracy)
now gives $`O(n)`$ depth under the same external-flag and sufficient
dirty-width orders. It charges source preparation/inversion and the
coherent one-hot program interface. The [blocked bilinear refinement](BLOCKED_BILINEAR_LOOKUP.md)
extends this to $`D_T=O(N/b^2+n)`$ throughout
$`b\ge17(L+n+7)`$, retaining the optimal-order count
$`T=O(\sqrt N+N/b)`$. Both orders match through
$`b\le\sqrt{N/n}`$ when the interval is nonempty.
At square-root-scale dirty width, the depth lower bound is still constant,
so depth optimality there remains open.
The lower bound does not assume count optimality; the upper circuit also
retains optimal-order T-count. No high-precision endpoint improvement follows.
The [uniform-precision theorem](UNIFORM_PRECISION_DEPTH.md) now removes
the separate $`n\chi(n)`$ term: its depth is
$`O(NL/b^2+nL)`$ with absolute constants. For
$`L=\Omega(\log n)`$, the earlier hybrid already absorbed that term
into $`nL`$. The new regime includes slowly growing
$`L=o(\log n)`$: rectangular blocked queries preserve the
$`\sqrt{NL}`$ count scale, and a uniform unary cutoff gives the
stronger $`O(NL/b^2+n)`$ depth for
$`6\le L\le\log_2(n+2)/16`$. The general count bound retains its
$`nL`$ term; outside the sufficient condition $`L\le N/n^2`$ or a
matching interval, count optimality need not follow. Earlier grouped
and parallel constructions remain eligible where they are sharper.

T-depth permits Clifford circuits of nonzero depth between its T layers.
It is not total circuit depth or elapsed QBP execution time. Neither the
exact CNOT light-cone bound nor the source's linear exact T-count minimum
supplies an additional depth lower bound in this model.

The [source-depth certificate](SOURCE_T_DEPTH.md) now settles the exact
geometric and paired sources within the specified Majorana-layer architecture.
The paired-tail schedule improves constants only. This rules out obtaining
sublinear precision depth merely by rescheduling those exact sources inside
that class; arbitrary Clifford interlayers and different approximate sources
remain eligible. Recent shallow rotation synthesis requires either growing
clean workspace or a prepared catalyst in the constructions audited
[here](RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026).
We therefore have useful ingredients, but no complete argument closing
either the general T-depth gap or the constant-clean frame endpoint.

The [error-accumulation audit](HOPF_ERROR_ACCUMULATION.md) now rules out
one proposed shortcut. Ideal angle perturbations obey a sharp square-sum
bound, and their finite relative spectrum is independent of the base
angles when each layer's error magnitudes agree. Independent nearest-grid
rounding still requires a logarithmic precision budget. The literal
shared-flag source stages instead admit a coherent linear-leakage family,
so their full isometry errors cannot receive the same square-sum guarantee.
These are restrictions of the specified error models, not a depth lower
bound. The [short-echo audit](HOPF_FLAG_ECHO.md) closes the proposed
two-half-angle square and three diagonal Pauli alternatives: generic
radial error remains linear, despite a special exact equal-mask identity.
The [radial filter](HOPF_RADIAL_FILTER.md) instead uses the standard
pi-over-three fixed-point sequence with six charged native phase words.
It gives quadratic full polar error on the same two flags and a physical
square-sum-plus-quadratic composition bound. Source widths may then use
$`m_d=L+4+\min\{n-d,\lceil\tfrac12\log_2(8n)\rceil\}`$.
This halves the logarithmic coefficient in assigned source precision,
but additional calls and phase words preserve the same asymptotic costs.
The earlier [conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md)
reduces the precision part to logarithmic depth when the active
suffix can reserve 7m temporary clean bits for source width m.
It uses the same two external flags and an enlarged, charged reflection.
With the unfiltered additive cap, source/reflection depth at fixed
accuracy becomes $`O(n\log\log(n+2)+\log^2(n+2))`$ after the late
fallback. Its [grouped-program theorem](GROUPED_PROGRAM_PREFETCH.md)
supplies the query/predicate composition: exact cached-program
erasure survives source leakage and changing internal addresses; chunked
dirty indicators avoid the late routing bottleneck. This earlier depth is
$`O(n\log\log(n+2))`$ with optimal-order count at sufficient
square-root-scale dirty width and fixed accuracy.

The [incremental group selectors](GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
now have $`O(g)`$ total depth in a height-g group, within the original
$`16m2^g`$ suffix reservation. They erase each suffix enable before its
controls change and retain prefix nodes until the final reverse traversal.
Their work returns exactly through source leakage and on arbitrary
inactive inputs. The resulting selector contribution is $`O(n)`$;
the source/reflection contribution retains the old geometric bound.
The [common-source identity](GROUPED_PROGRAM_PREFETCH.md#11-a-common-source-identity-and-the-remaining-reflection)
moves a shared preparation to the group boundaries but retains two
conjugated success reflections per stage. A legal exact row rules out a
stale success monitor and a one-use source-bank substitution, with constant
rejected norm. These are interface restrictions, not depth lower bounds.
The [unary construction](UNARY_PHASE_GRADIENT.md) instead uses a Fourier
eigenstate and constant-T-depth programmed cyclic shifts. Conditional
bilinear work returns exactly on all source inputs, the inactive sector
is literal identity, and source preparation plus actual inversion costs
at most twice its preparation error for the entire group. Its global
allocation proves the displayed $`O(n)`$ upper bound. The
[blocked bilinear theorem](BLOCKED_BILINEAR_LOOKUP.md) now closes the
width-dependent extension. It replaces the residual bank route by
selected bilinear blocks with one returned dirty phase helper, then
charges the coupled block-size/chunk-length allocation. The
[uniform-precision chapter](UNIFORM_PRECISION_DEPTH.md) now proves the
precision extension: rectangular indicator dimensions retain the correct
precision-dependent count, and the explicit unary cutoff is sublinear
uniformly throughout its stated low-precision range. The remaining
questions concern unrestricted large-width depth and the high-precision
constant-clean endpoint; neither is resolved by these two-clean schedules.

The [two-layer obstruction](SHALLOW_SOURCE_OBSTRUCTION.md) separately
allows unrestricted Clifford interlayers: the original source at width
at least five, or its controlled version at width at least four, remains
at operator distance at least $`1/16`$ from every two-T-layer full-input
circuit, regardless of dirty width. This does not cover a conditionally
initialized source or the two-clean complete-frame isometry.

### What this already gives Hopf QBP

The synthesis target throughout is the prescribed Hopf tree frame; an
arbitrary-unitary compiler is not required. The existing decoder's
universal observable guarantee still uses its marker directions, as the
[frame-safety necessity argument](FRAME_SAFE_COMPILATION.md#necessity-for-all-observable-dependent-gradient-means)
shows. A narrower observable class or a changed decoder would define a
separate task and require its own resource comparison.

The endpoint below is a high-precision compiler question, separate from the
precision needed for a fixed raw-gradient accuracy. For a reflection-sum
observable with coefficient one-norm Lambda, the existing
[QBP error allocation](QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights)
permits

```math
L=\max\{6,\lceil\log_2(32\Lambda/\varepsilon_\infty)\rceil\}.
```

Thus fixed observable scale and coordinate accuracy give fixed L. In this
regime, the arbitrary-budget real-frame compiler already attains
worst-case optimal-order $`T=\Theta(\sqrt N)`$ with zero compiler clean
qubits and $`b=\Theta(\sqrt N)`$. The separate parallel schedule gives
this count together with $`D_T=O(n^2)`$ using two compiler clean qubits
and a sufficiently large square-root dirty allocation. Both preserve the
prescribed frame used by QBP; their clean allocations cannot be interchanged.

The protocol reserves one additional clean interference qubit and the
controlled observable's work, beyond the initialized system and compiler
reservation. Its fixed-accuracy, fixed-confidence execution count is
$`O(1+\log n)`$; observable costs and classical gradient output remain
separately charged. The unresolved $`L=N`$ endpoint does not prevent
these existing QBP guarantees. Optimal T-depth and end-to-end gradient
optimality remain open.

## The remaining endpoint

The publication's established compiler results do not depend on resolving
this question. Let $`T^\star_{F,\mathbb R}`$ be the worst-case minimum
T-count for the prescribed complete real Hopf frame. At

```math
a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

the retained frontier remains

```math
\Omega(N)\le T^\star_{F,\mathbb R}
\le O(N\ell_*(n))=O(N\log_2^*N).
```

Here $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ for
$`0\lt\eta\le1/64`$. Literal phases, clean-work leakage, dirty-work
return, and arbitrary reference correlations are included in the
[complete-input error contract](FAULT_TOLERANT_COMPILER.md#1-target-resources-and-theorem).
The same upper bound holds with only one clean qubit and the same dirty
allocation. Zero-clean layerwise synthesis gives $`O(N\log N)`$ there;
a sufficiently large $`\Theta(n)`$ clean reservation instead attains
$`\Theta(N)`$. These are sufficient constructions, not necessary clean
reservations.

The question is whether the displayed two-clean allocation admits an
$`O(N)`$ T-count construction, or whether a stronger general lower bound
holds. The optimum minimizes T-count: a valid linear-T construction closes
this gap even with a larger fully charged Clifford count. Preserving
$`G=O(NL)`$ is the stronger joint goal of the retained constructions.
The current upper bound does not cover every prefactor in
$`b=\Theta(N)`$. Source-call minima and restrictions on particular
intermediate interfaces do not strengthen the unrestricted lower bound.
Closing this count endpoint would still leave the T-depth question open.

## Current assessment: what the results establish

The antichain and sparse-update results prove linear bounds for promised
families. The subsequent passes establish complete boundary actions, an
explicit repair, cheap source-width transitions, compact Cayley data, and
native four-, eight-, and sixteen-mode benchmarks. The latest joint word
has fewer declared source calls, but the comparison words have the same
leading precision cost after fixed-mask simplification.
**None has narrowed the generic upper/lower gap.** Wide independent changes
and deep sparse nesting each admit one precision charge; the remaining
question concerns precision reuse for unrestricted comparable updates.

| Level | Established scope |
|---|---|
| General theorems | Exact matching resources, sufficient-clean matching T-count, one-clean grouped bound, and matching T-count/depth in explicit accuracy/workspace ranges; constant-clean endpoint and the unrestricted depth frontier remain open |
| Incomparable promised families | Antichain and sparse ancestor-closed updates have $`T=O(N+L)`$, $`G=O(NL)`$; generic rounding satisfies neither promise |
| Reusable blocks | Cheap coarse frame, linear residual data, complete weighted dilations, two-flag assembly, and a coupled repair with exact child-call cancellation; the surviving target precision remains charged |
| Completed small-example pass | Fixed-address product compression, four-mode native factors, and a full-port changing-target word with borrowed-signal return; no precision recurrence independent of growing support or group count |
| Scoped diagnostics | Specified source, query, repacking, truncation, and shared-conjugator failures; no additive full-frame lower bound follows |

### A promised antichain class has a linear endpoint bound

The [antichain theorem](ANTICHAIN_COMPILER.md) requires a supplied native
baseline C with determinant-one local words of length $`O(n-d+1)`$
at depth d. The target agrees **literally** outside a prefix-free changed
set, possibly of size $`N/2`$. No closeness is needed. With zero clean qubits,

```math
b\ge L+n+7,\qquad T=O(N+L),\qquad G=O(NL),
\qquad
\|\widetilde W-W\otimes I_b\|\le\eta.
```

The exact factorization $`W=VMV^\dagger C`$ removes the native
descendant forest; charged dirty swaps pack the corrections into one SU(2)
multiplexor. The full-operator error includes core and borrowed-signal
return. Forest and packing helpers return exactly.

### Sparse nested updates also have a linear endpoint bound

The [sparse-update compiler](SPARSE_UPDATE_COMPILER.md) uses the same
literal native-baseline promise. Write $`\mathcal S`$ for the ancestor
closure of the changed nodes, including ancestors whose local words did
not change. This support set is distinct from the affine residual operator
S used below. Put

```math
m=|\mathcal S|+1,\qquad s=\lceil\log_2m\rceil.
```

Under the sufficient condition $`n\ge2s+32`$,

```math
a=1,\qquad b\ge L+n+7,\qquad
T=O(N+L),\qquad G=O(NL),
```

```math
\|\widetilde WJ_1-J_1(W\otimes I_b)\|\le\eta.
```

An exact off-support forest and charged basis packing leave m active
modes on s logical bits. Conditional logical zeros support one dense
residual dictionary, scalar source, and half-block amplification. Only the
external predicate is initialized unconditionally. Its isometry error
includes leakage, work return, and dirty references; it is not a
full-operator guarantee for arbitrary input on that external qubit.

A root-to-leaf path has $`|\mathcal S|\le n`$ and satisfies the
bound for all n: use this construction at $`n\ge44`$ and the proof's
finite fallback otherwise. The bound is asymptotic, not a practical estimate.

The two promises are incomparable. A wide antichain can have a large
ancestor closure; a sparse closure can contain arbitrarily long nested
paths and some branching. For both theorems, a concrete real Hopf family
fixes unmarked angles at $`\pi/4`$ and allows arbitrary real changes
at the promised nodes. Its complete-frame approximation preserves the
existing QBP error interface. A generic native baseline may have complex
local words: literal agreement with it does not automatically define a
real Hopf frame or the separate phase-dressed complex magnitude family.

Generic target rounding may change every internal node. Then
$`|\mathcal S|=N-1`$ and $`s=n`$, leaving no packed zero sector;
the sparse dictionary is also too large for its current accounting.
Parameter count and support dimension diagnose why this proof does not
apply. They are not hardness results. The unresolved structure is
widespread changes across comparable tree nodes, with no proved reduction
to a constant number of the two solved families.
Simply applying their theorems in sequence is not such a reduction:
each requires literal agreement with a short native baseline outside its
own marked set, a promise not automatically preserved after a nonnative
update.

## A sufficient construction to seek

A sufficient endpoint construction is an actual unitary Q using two
initialized flags and at most $`N+n+7`$ dirty qubits, satisfying

```math
\left\|2J_2^\dagger QJ_2-(W\otimes I_b)\right\|\le\eta/4,
\qquad T(Q)=O(N),\qquad G(Q)=O(N^2).
```

These displayed bounds are a sufficient joint target, not the definition
of the T-only optimum. Every source, query, control, inverse, and helper
belongs to the circuit accounting. The accepted action must hold on all logical and dirty
inputs, including reference correlations. There is no assumed return
condition on rejected branches. The
[normalization-two amplification lemma](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
then gives the complete-isometry compiler with three charged calls and the
same two flags.

The [two-flag assembly](RESIDUAL_ASSEMBLY.md) realizes this interface at
$`O(N+nL)`$ T cost. A faster whole-residual block followed by the
charged coarse frame also suffices. The uniform target $`T=O(N+L)`$,
$`G=O(NL)`$ is stronger than settling the endpoint alone.

## Reusable ingredients and the resource bottleneck

These are retained constructions; their linked proofs are the primary homes.

| Ingredient | Established capability and remaining cost |
|---|---|
| Native coarse frame C | Fixed-accuracy $`T=O(\sqrt N)`$, $`G=O(N)`$ in the endpoint pool, no initialized work, exact helper return; any new control or mask must be charged |
| Compact residual | $`O(N)`$ classical local generators and path products for $`C^\dagger W'-I`$; coherent evaluation and target transport are not free |
| Weighted forward block | Complete two-sector dilation, fixed normalization, $`T=O(N+nL)`$, $`G=O(NL)`$ at $`b\ge L+n+7`$; the native synthesis signal is borrowed, distinct from its initialized logical dilation signal |
| Exact-return alternative | Same component with $`T=O(L\sqrt N)`$, $`G=O(NL)`$ at $`b\ge n+\lceil\sqrt N\rceil+7`$; different width/return point, not the best endpoint count |
| Affine/reverse assembly | [Two-flag proof](RESIDUAL_ASSEMBLY.md): diagonal-plus-forward and inverse-based reverse blocks give normalization two, with controls/inverses charged; cost still contains $`nL`$ |
| Full-port hierarchy | [Fusion audit](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths): recursively closed complete unitary, linear classical generators; faster native synthesis unproved |
| Coupled completion | [Whole-residual boundary](RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost): one-signal normalization-two recursion; a commutator implements the complete rank-at-most-four repair without normalizing transported differences; its child calls cancel to the original target wrappers, whose precision remains charged |
| Grouped frame | Best general endpoint bound, with conditional suffix and core return in the complete error |
| Product-first residual coordinates | [Cayley recursion](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) gives constant-size local data for the complex-coarse residual; the [four-mode native benchmark](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) uses two one-target programs, while the eight-mode controlled coupling and general coherent conversion remain charged |
| Source-width transport | [Reverse-order loader](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation): total boundary T-count $`2(m_{\max}-1)`$ with linear native loaders; transformed group bodies remain charged |
| Joint source body | [Changing-target word](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) shares one fixed scalar conjugator and returns a borrowed signal through two parity boundaries; the eight-/sixteen-mode words are priced, but fair fixed-mask simplification leaves the same leading precision cost in both comparisons |

The [weighted norm proof](ENDPOINT_TREE_TRANSPORT.md#the-actual-weighted-pieces-have-no-height-penalty)
avoids a height penalty via uniform coarse subtree accuracy. The
[forward block](WEIGHTED_TRANSPORT_BLOCK.md) charges marker gathering,
queries, and certified preprocessing at zero defects. Real-rotation X
symmetry permits a borrowed native signal; the scalar-phase word lacks
that established symmetry.

For the common algebraic target approximation W', the exact residual is

```math
C^\dagger W'=S+R,\qquad S=A+F,\qquad
\alpha=4\varepsilon_0,\qquad \sigma=2-\alpha.
```

For fixed $`0\lt\varepsilon_0\le1/64`$, the existing normalized
branches satisfy

```math
\|S/\sigma\|\le
\frac{1+2\varepsilon_0}{2-4\varepsilon_0}\lt0.54,
\qquad \|R^\dagger/\alpha\|\le\frac12.
```

Their local generators obey the additional coupled relations

```math
\begin{pmatrix}g_v&k_v\\h_v&d_v\end{pmatrix}
=(U_v^C)^\dagger
  \mathrm{diag}(g_{2v},g_{2v+1})U_v^{W'}.
```

Together they make $`S+R`$ unitary. The current affine and reverse
dilations use these data as separate contractions. Their assembly gives a
sufficient bound

```math
T=O(T_S+T_R+N+L),\qquad
G=O(G_S+G_R+NL).
```

Thus reducing their combined cost to $`O(N+L)`$ would work. It is one
sufficient route, not a requirement on a solution. A new word may use the
coupled unitary directly and choose different rejected completions. No
constant-query conversion of an opaque forward-only block to the affine
block has been proved.

### Why a single larger group is not already the answer

In the retained grouped construction, a group of s levels above a suffix
of length r, with $`e=n-r`$, has

```math
w_g=O(\log(s+2)),\qquad
Q_g=\Theta(s2^{n-r}),\qquad
m_g=L+\lfloor r/4\rfloor+8.
```

These are private initialized work, padded coefficient-table rows, and
precision-source width. Exponentially growing groups give
$`\sum_gQ_g=O(N)`$ and $`O(\ell_*(n))`$ precision charges.
A single group spanning all but a constant suffix has
$`\Theta(nN)`$ current table rows and an $`O(nN^2)`$ lookup
Clifford estimate at $`L=N`$. Its current streamed coarse-program
estimate also grows to $`O(nN)`$ under that representation's uniform
star-normalization condition. Its current term label needs
$`\Theta(\log(n+2))`$ private initialized bits, which a constant
suffix cannot supply. Removing only its depth label or private
workspace does not establish a linear-T theorem. The expanded Clifford
estimate separately misses the stronger joint goal; it is not a T lower bound.
These are costs of the current representation, not gate lower bounds.
The selected dirty allocation also falls below the proved banked
refinement's sufficient threshold.

The weighted representation can retain the cheap coarse frame instead,
but its downward path factors still use target transport. Calling the
target frame to supply those factors would be circular; replacing it by
coarse transport alone misses mixed corrections. The selected joint goal therefore requires controlling table expansion
and coarse programs as well as precision cost; a T-only improvement must
still charge those operations but may have a larger Clifford count.

## What the failure diagnostics actually rule out

None of these statements makes separate source or mask costs additive for
a jointly synthesized circuit or settles the unrestricted endpoint.

| Shortcut | Established limit and scope |
|---|---|
| Make the same exact source sublinear | Minimum exact T-count on $`m\ge2`$ dirty wires is $`2m-4`$, or $`2m-2`$ when controlled, allowing returned helpers; [primitive bound only](OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit) |
| Hoist the source basis and use cheap masks | Valid transformed masks have linear exact/fine-accuracy cost; the actual paired source has a mask requiring $`T\ge q/2-6`$ at error $`2^{-q}`$ with full return; [separate-mask bound only](SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks) |
| Carry a source code and renew it after each query | Width changes are cheap, but a legal scalar query leaves the flag-correlated code by constant norm; complete syndrome renewal costs $`\Omega(L)`$ at compilation accuracy; [specified interface only](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost) |
| Multiply accepted blocks sharing flags | Rejected components return coherently; [full word required](SOURCE_REUSE_LIMITS.md) |
| Telescope one global conjugator through a fork | Cancellation is valid, but the remaining native word has rejected returns; a coefficient-ellipse bound excludes even arbitrary mask retuning of that word at fine precision; [fork audit](SOURCE_REUSE_LIMITS.md#5-a-shared-conjugator-does-not-close-a-branching-fork) |
| Connect coupled blocks through their zero-defect word | The accepted error is exactly $`-3(A-I)(B-I)/8`$; the [complete commutator repair](RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost) cancels to the original fine-precision wrappers; its small norm does not suppress the whole merge's local synthesis error |
| Unload a query after changing its address | The inverse can leave dirty-dependent action; [query-scheduling counterexample](SOURCE_REUSE_LIMITS.md) |
| Replace initialized nilpotent source by dirty encoding | The specified full-output scalar relation needs $`\log_2L-O(1)`$ initialized width, regardless of dirty width; [that interface only](SOURCE_REUSE_LIMITS.md) |
| Compress a whole frontier to one scalar | Generic depth-d cut rank is $`2^d`$, although each edge has rank one; existing logical modes carry it, so this is [not an ancilla lower bound](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Permute a merged band to bounded-size blocks | The chosen completion has connected support growing with height; other completions or nonpermutation bases remain allowed; [fixed-completion restriction](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Truncate propagation using its norm margin | Constant Riccati messages coexist with undamped continuation; local-depth truncation misses a fixed bottom component, but the witness has a cheap global circuit; [word-specific failure](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths) |
| Use fixed-order coarse transport | Mixed corrections survive at fine precision; [specified-order failure](ENDPOINT_TREE_TRANSPORT.md) |

For clarity, the rejected-return term is an exact operator identity. With
$`P=JJ^\dagger`$ and $`B_i=J^\dagger Q_iJ`$,

```math
J^\dagger Q_2Q_1J
=B_2B_1+J^\dagger Q_2(I-P)Q_1J.
```

For $`Q_1=Q_2=H\otimes H`$ on two flags, the separate accepted
blocks are $`1/2`$, yet the product's accepted block is one, not
$`1/4`$. A successful global word may exploit rejected returns, but may
not discard them or assume uncharged initialized history prevents them.
Similarly, the current chosen forward completion need not approach its
zero-defect completion when its accepted map becomes small.

The complete fusion hierarchy retains every input column: a parent and two
nonterminal children form a ten-mode unitary, and a closed subtree with m
internal nodes uses $`2m`$ modes. These are basis modes, not additional
clean wires. Its linear list of local scatterers is compact, whereas eager
entrywise expansion of an h-level affine map has $`h2^h+1`$ generic
nonzero entries. Neither linear classical storage nor dense expansion
supplies the missing native synthesis by itself.

## Revision decision and next bounded pass

The current selection is the
[revision checkpoint above](#revision-checkpoint-and-selected-next-test).
This section retains the earlier source-carry decision and its boundaries.

The coupled-merge passes have completed their structural task. The
[commutator repair](RESIDUAL_ASSEMBLY.md#9-a-repair-word-without-an-ill-conditioned-transported-basis)
works on every signal port without normalizing transported differences.
Its actual inverse calls cancel even for a noncanonical approximate child
word. The reduced circuit retains one child call and the original local
target wrappers. The complete merge has parent error
$`\|\widehat B-B\|/2`$, with no small residual factor, so the repair
does not discount their precision. Its layerwise emission remains
$`O(N+nL)`$. The separate native shared-conjugator fork fails even
after arbitrary mask retuning. Another proof of these boundaries or
cancellations would not improve the resource frontier.

The source-carry pass tested **precision reuse between adjacent existing
groups of the best grouped compiler**. It resolves the width-only part of
the proposed transition, while exposing the missing program cost. A fixed
fork or a fixed number of fused groups still changes only constants; the
needed gain must survive a variable number of unequal groups.

### What the source-carry pass established

The [native width analysis](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
keeps literal phases and proves full-space identities, including occupied
flags and correlated dirty inputs.

| Candidate | Result | Consequence |
|---|---|---|
| Original chain loader, extracted from each group | The required one-bit eigenbasis bridge has exact T-count $`2m-1`$ and needs at least $`L-3`$ T gates at error $`2^{-L}`$ in the stated width range | The cheap opposite-order product is not this bridge |
| Reverse-order star loader | Same certified coefficient grid; monotone source boundaries cost exactly $`2(m_{\max}-1)`$ T gates and $`O(Rm_{\max})`$ Clifford gates | Width changes are solved for this choice, but transformed group programs are excluded from that count |
| Separately synthesized transformed mask | A realizable grouped coefficient $`1/(9s)`$ requires $`T\ge\max\{0,m-2\log_2s-9\}`$ at error $`2^{-m}`$ | The old table-row cost cannot simply be assigned to masks in the new basis; this is not an additive compiler lower bound |
| Flag-correlated chain-source code | Preparation is charged; width changes cost their size difference. A realizable $`c=1/8`$ query has leakage norm $`\sqrt7/4`$ | A source cannot remain a cheap flag X after that query without a changed boundary |
| Complete syndrome renewal | Needs $`T\ge(L-4)/2`$ at compilation accuracy under its stated width condition | This particular renewal interface reintroduces precision cost; code-restricted or deferred alternatives remain open |

The [flag-correlated proof](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost)
does not assume the released source tail becomes clean. It also does not
supply a complete group program on the code. The retained generic resource
frontier is unchanged.

### Revised sufficient ledger: jointly compile the interior program

Keep the groups and precision allocation of
[conditional-suffix Sections 6--7](CONDITIONAL_SUFFIX_COMPILER.md#6-exponentially-growing-groups-and-the-explicit-workspace-ledger).
For $`n\leq r_0`$, the fixed-depth fallback already costs $`O(N+L)`$.
For each remaining group, define

```math
\begin{aligned}
Q_g&=\Theta(s_g2^{n-r_g}),& m_g&=L+\lfloor r_g/4\rfloor+8,\\
w_g&=O(\log(s_g+2)),& k_g&=n-r_g+O(\log(s_g+2)).
\end{aligned}
```

Q counts table rows and the same order of coarse-program work; w is
private initialized work inside the active suffix, and k counts dirty
selectors. The existing sums are

```math
\begin{aligned}
\sum_gQ_g&=O(N),&\sum_gQ_gm_g&=O(NL),\\
R&=O(\ell_*(n)),& R&\leq n,\\
m_{\max}&=L+O(n),&\sum_g|m_{g+1}-m_g|&=O(n).
\end{aligned}
```

The new loader realizes the last line's width ledger. What remains is a
**joint interior circuit**, including transformed programming, with total
T-count bounded by

```math
O\!\left(\sum_gQ_g+L+RP(n)\right),
```

where P is a fixed polynomial independent of L. The bound permits a
charged $`O(L)`$ use *inside* that joint circuit. Dirty-only outer
conjugation leaves the logical error contract unchanged, by
[the source-hoisting argument](SOURCE_REUSE_LIMITS.md#the-current-paired-source-also-has-expensive-transformed-masks).
Thus its boundary loaders alone cannot supply precision to independently
accurate, precision-free group bodies. The target is a global cost bound;
it is not a claim that every group separately costs only $`O(Q_g)`$.

Combining such a proved interior with the known boundary ledger would give

```math
T=O\!\left(\sum_gQ_g+L+m_{\max}+RP(n)\right)=O(N+L).
```

This remains an **unproved sufficient construction**, not a theorem or a
necessary architecture. A changed encoding correlated with logical data
may avoid separate group outputs entirely. Leave the fixed deepest-layer
tail with its existing $`O(N+L)`$ compiler. Retaining $`G=O(NL)`$
requires separately charging all interior Clifford gates; the new loader
and bridge costs themselves fit that budget.

### Small examples: completed results and remaining cost

The small-example pass supplied the following reusable results. Their
linked chapters retain the derivations and native ledgers.

| Result | Established capability | Remaining cost or limitation |
|---|---|---|
| [Fixed-address products](SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis) | Arbitrarily many noncommuting one-qubit factors compress into four quaternion coordinates and compile once | Changing the quantum address or growing logical support invalidates that fixed two-mode table |
| [Cayley residual](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) | Complete complex-coarse residual, linear classical data, two-dimensional local updates, stable inverse conversion | Generic coherent evaluation is unpriced; direct Woodbury emission returns to the fine target wrappers |
| [Four-mode native benchmark](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) | Two magic-basis factors give $`T=O(2^k+L)`$ and $`G=O(2^kL)`$ at an unchanged prefix, with one clean flag and the stated borrowed spectator | A fixed-support corollary; the separately charged complex coarse inverse remains necessary |
| [Eight-mode coupling](ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling) | Explicit controlled Bell-projector rotation and a four-Pauli comparison | Independent prefix/child factors cannot absorb the root; this does not exclude joint synthesis |
| [Shared changing-target word](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) | Exact common scalar conjugation, full SU(2) signal error, and two parity boundaries returning an arbitrary borrowed signal | Fixed-mask simplification gives the comparison words the same leading precision cost; variable-depth allocation and joint programming remain unproved |

The last construction reduces declared source counts from $`45g`$ to
$`27g+6`$: 135 to 87 at eight modes and 180 to 114 at sixteen. Its
fixed mask commutes with the loader tail, leaving a two-T seed; each
hoisted fixed-A interior then costs at most eight T gates. Both words
retain five programmed B words per stage and the same leading term
$`(40g+4)q`$. The fixed-A remainder changes from 80g to
$`8(4g+2)`$ before further optimization, with queries and controls
separately charged. These are declared-word bounds, not minima.

For local depth three or four, its stated allocation fits
$`a=2,b=L+n+7`$ and gives $`T=O(2^k+L)`$ at a k-bit unchanged
prefix. The source-width margin grows with variable depth. Neither the
fixed-depth cost nor its full borrowed-signal return narrows the generic
endpoint gap. The source-appearance reduction is not a leading precision
improvement or an additive lower-bound argument.

### Revision after the joint native audit

The small-example pass is complete as structural infrastructure. The
phrase "joint synthesis of the programmable masks" identifies a missing
result, not yet a construction. Continuing to simplify two ordinary
layers without specifying how the saving scales risks another constant
improvement to a weaker baseline.

| Comparison | What is already available | What would constitute progress |
|---|---|---|
| Fixed logical support | A fixed number of modes costs $`O(L)`$; one unchanged-address two-mode product compiles once | A representation whose charged program remains controlled as the support grows |
| Ordinary layers | The retained baseline costs $`O(N+nL)`$; the shared body changes its displayed constants | A precision recurrence that improves the best grouped construction, or a separate complete construction beating it |
| Existing unequal groups | $`T=O(N+LR)`$, $`R=O(\ell_*(n))`$, with linear total table work | A jointly emitted program with one global precision charge and a valid allocation throughout |

The unmodified paired-source real-Y word did not supply a grouped scalar
SELECT. That interface keeps forward and actual-inverse scalar branches,
four complex phases, selected column maps, a private term label, and a
reflection on the initialized active suffix. Its source width and suffix
predicate change between groups. Even the one-clean extension retains
this scalar/atom separation; it relocates a flag and changes the deepest
tail compiler. That transfer was the unresolved compatibility question at the revision.
The canonical audit below now supplies it, without improving the precision
recurrence.

Keep the recent proofs and fixtures; do not add another fixed-size
optimization as an endpoint advance. Cheap width transport is available,
fixed-mask simplification is accounted for, and the remaining general
lower bound is still only the retained one. A failure of a selected
architecture would not make the grouped factor necessary.

### The canonical group test is complete

The [canonical scalar audit](CONDITIONAL_SUFFIX_COMPILER.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge)
transfers the small native rotation to the actual grouped interface. It
uses $`R_c=\exp[-i\arccos(c)Y_\sigma]`$ on the existing scalar
flag, with a separate atom flag and an arbitrary borrowed synthesis
signal outside the initialized-work reflection. The canonical accepted
entry is c, including structural zero. Forward and reverse atoms keep
their literal phases and query order.

Direction-controlled $`Z_\sigma`$ gates select the actual native
inverse of a single positive-angle program, even after rounding. The old
one-tail scalar has the same inverse symmetry, so this one-program SELECT
is available to both baselines. The canonical substitution is not uniquely
responsible for that constant improvement.

The new full-operator error ledger gives one inner program per half-block
and three per outer amplification. With $`q_g=m_g+4`$, the group error
is below $`(45/8)2^{-m_g}`$, inside its retained budget. The paired core
and borrowed signal add six dirty slots before the fixed helper constants;
the existing group slack absorbs this by increasing its fixed thresholds.
No additional external clean qubit is needed. Final h and private-work
return remain part of the complete approximation error.

The exact common-source identity also passes through the non-scalar
group operations. But its cost fails the endpoint test: at a legal common
q, the paired comparison words have the same leading precision term
$`(120R+4)q`$. The old scalar uses six programmable mask appearances
per group, and the canonical replacement uses thirty. Unequal widths do
not acquire a free bridge. The proved totals remain
$`T=O(N+L\ell_*(n))`$ and $`G=O(NL)`$.

**Close this canonical-completion candidate for precision amortization.**
It supplies a correct group interface, exact inverse routing, and a
complete error/workspace ledger. It does not supply a variable-group
fusion rule or evidence that the endpoint is close. Further mask sweeps,
larger instances, and constant cancellation in this word are not the
next task. This is a conclusion about the displayed construction, not
an unrestricted lower bound.

### Selection audit after canonical completion

The 2026-10-01 revision checks three alternatives before selecting another
construction. The comparison is at the actual endpoint, where the accuracy
bit count is $`L=N`$; factors polynomial in L cannot be hidden in a
soft-O estimate.

| Mechanism | What it supplies | Why it is not yet an endpoint construction |
|---|---|---|
| Whole-residual Hamiltonian synthesis | A generic route from a Frobenius-small residual to a complete unitary channel | The explicit bound below still charges precision polynomially and uses initialized work; a tree-specific adaptation is missing |
| State preparation followed by Householder reflections | A prepared state determines its rank-one reflection, irrespective of the preparer's other columns | The prescribed frame needs more than its first column; neither constant reflection count nor a two-clean implementation follows |
| Reduce the shared-source proof's clean work | The sufficient-clean construction already charges precision once globally | Its source column, private SELECT labels/flags, and failure history each use initialization; dirty-safe lookup banks do not replace those roles |

**Small residual norm is already affordable.** The
[Frobenius-coarse specialization](ENDPOINT_TREE_TRANSPORT.md#frobenius-small-residuals-also-fit-the-endpoint-budget)
uses no clean qubits and returns the dirty pool exactly, with
$`T(C)=O(n\sqrt N)=O(N)`$ and $`G(C)=O(Nn)`$, while making
$`\|C^\dagger W-I\|_F\le\varepsilon_0`$ for fixed
$`0\lt\varepsilon_0\le1/64`$. This is a specialization of the
existing borrowed compiler, not a new frame bound. It does not price a
coherent logarithm or matrix-function evaluator.

[Fang–Heunen–Wang, Corollary 3.9](https://arxiv.org/html/2607.12907v1)
give, for bounded Frobenius distance to the Clifford group and accuracy
$`\varepsilon`$,

```math
\begin{aligned}
T&=O\!\left((N+\log\log(1/\varepsilon))
             (n+\log(1/\varepsilon))^2\right),\\
a&=O\!\left(N(n+\log(1/\varepsilon))\right).
\end{aligned}
```

At $`\log(1/\varepsilon)=\Theta(N)`$, these displayed bounds are
$`O(N^3)`$ T gates and $`O(N^2)`$ initialized ancillas. Definition 1.1
traces out that initialized environment; it does not certify literal phase
and arbitrary dirty-work return. This screens out direct substitution of
that theorem, not a structure-sensitive adaptation or an optimality claim.

[Gosset–Kothari–Wu, Theorem 1.1 and Section 1.1](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/)
give optimal state preparation with initialized ancillas and the
Householder-based K-column upper bound

```math
T=O\!\left(K\sqrt{N\log(K/\varepsilon)}
             +K\log(K/\varepsilon)\right).
```

Its direct complete-column specialization $`K=N,L=N`$ is $`O(N^2)`$,
before any clean-work adaptation. There is also a simple structural check
on the proposed constant-reflection shortcut: set every upper-tree angle
to zero and all bottom pair angles to one small nonzero value. Then
$`\mathrm{rank}(W-I)=N`$, whereas a product of k rank-one
reflections differs from identity by rank at most k. For accuracy below
the smallest singular value of $`W-I`$, the same witness excludes such
a k<N approximation. This elementary rank argument concerns that literal
reflection product; packed higher-rank reflections and fully charged
changes of basis remain allowed. The witness itself is an antichain and
already has a linear compiler.

The clean-work audit is likewise specific. In
[shared-source Sections 7–9](FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel),
the prepared source supplies the full-output relation, the private labels
and separate zero flags define the selected contraction, and a known-zero
counter prevents rejected paths from returning to the accepted block.
Changing the counter alone leaves the other initialized roles intact.
The existing nilpotent-source restriction remains scoped to its stated
interface; it is not a full-kernel or compiler lower bound.

These checks select no construction. They focus the next question on
preserving the tree structure in the native program, rather than merely
making the residual small, preparing its state column, or renaming clean
work as borrowed work.

### The Hopf-specific scattering test

The [packed scattering construction](ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
uses the actual Hopf tree. One n-bit address selects its node rotation;
one additional mode qubit accommodates all internal continuation ports.
Explicit port permutations produce a unitary step

```math
\mathscr S=\begin{pmatrix}D&C_{\rm out}\\ B&A\end{pmatrix},
\qquad W=D+C_{\rm out}(I-A)^{-1}B.
```

The second identity gives every prescribed real-frame column, including
singular angles. The genuine internal block is nilpotent; two decoupled
dummy ports are padded by minus identity, so the padded A is not nilpotent.
The inverse in this algebraic identity is not a free circuit operation.

One step is nevertheless fully priced: the returned-work port permutations
cost $`O(n^4)=O(N)`$, and the retained real-rotation compiler, with four
invariant address sectors, gives $`T=O(N+L)`$, $`G=O(NL)`$, at
$`b=L+n+7`$. Its synthesis uses no initialized work; one available clean
qubit selects external ports. This is a full-operator approximation to the
step, not a compiler for its feedback boundary.

The conversion cannot use only a fixed number of queries to that unchanged
local-angle step. On the all-right path, set all n angles to theta. The
frame entry is $`\sin^n\theta`$, while any k-query word in the step and
its actual inverse, with parameter-independent interleaves, has amplitude
degree at most k. For $`k\lt n`$, uniform error on the canonical interval
$`[0,\pi/2]`$ is at least $`(8e^2)^{-n}>2^{-6n}`$.
Thus the fine endpoint rules out this constant-query conversion for
$`n\ge5`$, even with arbitrary fixed basis changes and extra work.
The small [scattering checks](../tests/test_hopf_scattering.py) test the
complete transfer, dummy handling, port routing, and the path coefficient.

This closes the unchanged-coin conversion as a bounded-call candidate.
It does not establish a native T-count lower bound: global coefficient
preprocessing, target-dependent interleaves, and shared native synthesis
of many calls remain outside the query argument. The Hopf-specific
representation is useful, but it has not narrowed the endpoint gap.

### A state-only route for raw Hopf gradients

The completed [state-based theorem](STATE_BASED_QBP_THEOREM.md) is the
current task-specific statement. It uses two compiler flags to prepare a
real or gauge-fixed complex Hopf state, or coherently select it with the
actual coarse reference. Its initialized-isometry error includes flags,
arbitrary dirty work, and external references. The prescribed other frame
columns are not compiled.

The earlier [leaf-reference protocol](REFERENCE_STATE_QBP.md) retains its
own tradeoff: its derivative-envelope factor can grow to $`n+1`$, and a
fixed local branch-probability promise makes that factor constant. This is
an optional predecessor, not a restriction on the current theorem.

### A coarse inverse removes the angle-dependent shot penalty

The [coarse-frame proof](COARSE_FRAME_QBP.md) uses the actual native
reference $`C|0\rangle`$, a charged $`C^\dagger`$ in magnitude readout,
and randomized X/Y measurements. Classical scores use that same recorded
C, so its coarse discrepancy does not become gradient bias. Each magnitude
depth block has record norm at most five, including singular tuples.
The theorem retains the logarithmic depth-block shot order and
$`O(S+Nn)`$ histogram arithmetic, with coefficient and bit costs charged.

### The complex coarse interface closes after fixing a common phase

The [complex coarse construction](COMPLEX_COARSE_COMPILER.md) fixes the
arithmetic-mean phase gauge and implements its actual logical C with exact
dirty return. The [complex decoder](COMPLEX_COARSE_QBP.md) estimates both
magnitude and leaf-phase coordinates with the same two compiler flags.
All coherent branches and classical coefficients use the same gauge;
physical energy gradients are unchanged. The consolidated theorem keeps
this phase convention separate from literal complete-frame synthesis.

### What the complete task-cost comparison establishes

The [comparison](QBP_COST_COMPARISON.md) uses original accuracy bits K,
state precision $`P=\max(n,K)`$, and observable precision K. Its fine
gauged complex borrowed baseline and literal bank reservations prevent an
unfairly expensive original comparator. Extra dirty banks improve the
available state T bound only in its stated precision/workspace regimes.
Fixed accuracy favors the original bound; no gradient optimality or
end-to-end speedup follows.

The [bounded-input audit](BOUNDED_INPUT_QBP.md) supplies polynomial
construction for the listed grouped/state alternatives and preserves the
separate fine native-search caveats. Its explicit Pauli model also has
deterministic all-gradient and classical term-sampling baselines, with
shared duplicate aggregation and coefficient scale. The
[algebraic residual helper](RESIDUAL_TABLE_PREPROCESSING.md) outputs
certified classical coefficients. The
[native residual bridge](NATIVE_RESIDUAL_ROTATION.md) now converts them
to an elementary unaddressed U(z) row or an enabled two-row table, with
arbitrary borrowed inputs and the retained full-operator error. The
table preserves address and enable exactly and is identity when disabled.
These components implement existing identities.

### Next bounded task and stopping rule

The task-specific theorem, quantum/classical cost comparison, and
bounded-input construction are complete and consolidated. The bounded
[real](NATIVE_COARSE_QBP.md) and [complex](NATIVE_COMPLEX_COARSE_QBP.md)
native examples verify full-input phases, dirty return, and gradient means
on their special targets. Their original-frame comparators are cheaper;
the fixtures establish integration, not a general advantage.

The selected state-based depth pass is also complete, as recorded below.
The unaddressed row, enabled two- and four-row tables, and bounded
[one- and two-qubit residual preparation](NATIVE_RESIDUAL_STATE.md) implementation
passes are complete. Both preparations use two declared clean flags, literal
amplification phases, the actual inverse, and a full dirty-input isometry
bound. Their dirty allocations exceed the respective basic small-system reservations;
they do not replace the minimum-budget fallbacks. The
[four-row lookup](NATIVE_RESIDUAL_ROTATION.md#6-four-rows-with-two-address-bits-and-no-additional-helper)
adds the second address without any extra helper, using exact Clifford
conjugation of quadratic-mask Toffolis. Both addresses and enable return
exactly and the inactive sector is identity. The two-system-qubit word
charges its enlarged initial reflection at 28 T gates, reusing and
returning an arbitrary core helper on every leaked input. The bounded
[one-system-qubit coherent selector](NATIVE_RESIDUAL_BRANCH.md) also
implements reference/target selection under an arbitrary protocol branch,
preserving its relative phase and excluding it from the initial reflection.
Its dirty allocation exceeds the basic n=1 reservation by four wires.
The [native residual QBP integration](NATIVE_RESIDUAL_QBP.md) now connects
this selector to an actual coarse circuit, controlled observable, actual
inverse magnitude readout, direct phase readout, and exact histogram
decoders for a fixed complex one-qubit target. Its coefficient and full
coarse-distance promises have rational certificates. The
[coverage audit](VERIFICATION.md#state-based-qbp-coverage) now pairs the
stated claims with their analytic proofs and bounded implementations.
It separates three optional software deliverables: a certified front end
for the admitted bounded Hopf inputs, a variable-size or banked native
schedule, and general guarded decoder arithmetic. Their absence limits
executable coverage without leaving a prerequisite of the stated theorem
unresolved. The small-system search and workspace exceptions remain.

The selected Hopf-QBP construction and bounded native integration are
complete. No further fixture expansion or general software API is selected.
A new implementation pass should specify its Hopf-QBP input/output
contract and the missing interface it will supply. A new scientific pass
should identify a claim beyond the established bounds and its proof
obligation before producing another fixture. The complete-frame endpoint
and depth optimality remain separate research questions; neither is
needed for this state-based task. Application-level advantage is outside
the current research scope.
The existing cost comparisons and classical baselines remain valid
boundaries; they are not pending work.

#### Selected state-based depth audit

The [state-based depth theorem](STATE_QBP_DEPTH.md) closes this audit for
both real and consistently gauged complex charts. With
$`P=\max(n,K)`$, $`B_0=P+n+7`$, and two compiler flags, it gives:

| Sufficient dirty width | Compiler T-depth | Compiler T-count |
|---|---|---|
| $`b\ge2B_0`$ | $`O(NP/b+P+n^3)`$ | $`O(NP)`$ |
| $`b\ge16(B_0+\sqrt{NP})`$ | $`O(P+n^3)`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ |

Both schedules have $`G=O(NP)`$, with count and depth holding on the
same circuit. Exact query replacement preserves the borrowed-signal and
amplification contracts. The actual coarse reflection interpreter is
charged separately; its predicates surround whole row words, yielding
an $`n^3`$ overhead in both retained schedules. Coarse work returns exactly,
while fine state work retains its full initialized-isometry error.
The proof includes the coherent common reference, complex prefix phases,
literal live-register inequalities, and the separate small-system source.

The [complete task comparison](QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison)
adds the actual coarse inverse, observable depth at K, and all S or
$`2S`$ executions. Existing original real-depth schedules and eligible
complex count-based alternatives remain charged at common physical work.
The count-preserving real upper expressions can differ even when their
T-count orders agree. This is not a lower bound for the old protocol or
a total-runtime advantage; the retained Clifford and classical work,
including explicit Pauli baselines, still matters.

The pass uses the already proved exact dirty routing and indicator
identities, whose general lookup mechanism is inherited from
[Low–Kliuchnikov–Schaeffer](https://arxiv.org/html/1812.00954v2).
No new source-parallelization or lookup primitive was required. The
selected stopping condition is met. Matching T-depth, elementary-depth
improvements, the fine native emitter, and an end-to-end advantage model
are separate questions; none is automatically selected by this result.

#### Separate complete-frame question

No new complete-frame endpoint construction is selected. For that separate
question, the remaining selection target is a tree-specific whole-residual
factorization whose programmed coefficients already combine tree levels.
The unchanged-node scattering table has now been tested; another walk-power
or fixed-routing conversion of that same table is not the next candidate.
The question is whether globally preprocessed coefficients give a bounded
number of complete native primitives already priced here, such as diagonals
or one-target multiplexors, with every change of basis charged.
This is a sufficient representation to investigate, not an available
factorization or a claim that it exists. A packed higher-rank reflection
would likewise need its own complete native implementation.
An expanded native program with many calls but one precision charge also
remains eligible; the scattering query restriction does not price it.

The next deliverable is one explicit candidate identity and its symbolic
cost, or a concise explanation that none was found. Do not begin another
fixture pass merely because a matrix representation is compact. A different
way to retain precision between groups remains eligible if it supplies
a new complete transition rule. A stronger full-frame lower bound is also
a valid resolution, but would need an invariant that applies to arbitrary
circuits; the existing interface restrictions do not provide one.

Before committing to another construction pass, require an explicit
native operation or encoding rule and a symbolic cost hypothesis. For
the retained grouped route, a sufficient target remains

```math
T_{\rm joint}\le c_RL+O\!\left(\sum_gQ_g+RP(n)\right),
```

where $`c_R`$ is bounded independently of R and P is a fixed polynomial
independent of L. Table and coarse-program work must combine additively.
A different complete construction may bypass the groups and is judged
against the actual endpoint; the displayed form is not compulsory.

For a proposed bounded factorization, count the total table entries, basis
changes, and native controls, and allocate each factor's error before
claiming a resource gain. Constant shifts in precision can consume the
fixed dirty margin; $`b=N+n+7`$ must hold literally, not only as
$`b=\Theta(N)`$. A precision-dependent number of separately synthesized
factors does not pass the selection check merely because each factor uses
the existing optimal one-target compiler.

The selection check must identify where the precision-dependent native
work is paid, how its count scales, what is initialized, and what returns
on arbitrary dirty/reference inputs. Imported channel or state-preparation
results must first be checked against this complete-unitary contract.
Classical compactness and a bounded inverse condition number do not
supply a free coherent evaluator.

Only after a concrete rule passes this check should small examples test
it. For a grouped proposal use unequal groups, retained complex coarse
words, and actual inverses; require the same rule to close on a third
group before treating its recurrence as plausible. Intermediate group
action may be deferred, but final decoding and all rejected action must
be charged. If no rule qualifies, record that there is no selected lead
rather than assigning another generic mask or boundary audit.

The live-wire contract remains $`b=L+n+7`$ arbitrary dirty wires and
at most two external clean flags. Holding $`m_{\max}`$ throughout
cannot be assumed to fit uniformly in n: together with the deepest
group's $`n-r_0+O(1)`$ selectors it can exceed the allocation. A proposed
schedule must prove valid tail release or reassignment, or another valid
allocation. The proved transitions permit arbitrary tail correlations,
not initialized tails.
Changing suffix predicates retains only its existing active-sector
$`w_g`$ reservation. Include actual query unloading, literal inverses,
complete dirty/reference return, and all transition approximation errors.
If using an accepted-block construction, preserve normalization two and
account for rejected action in amplification; a direct complete-frame
isometry proof is also allowed.

Stop a candidate if its expanded word still pays a length-L source or
separately synthesized mask per group, requires full syndrome renewal per
query, expands the total table/coarse-program work beyond $`O(N)`$, or
uses an unproved clean-work or return assumption. These stop conditions
apply to this sufficient route; they are not general lower bounds. Record
failure in the existing proof home and change the native candidate. Do not
repeat the completed width, mode-closure, or commutator-repair audits.

Small fixtures can suggest or falsify a recurrence; they cannot prove its
asymptotic cost. No target-frame oracle, free evaluator, reset, initialized
history, or supplied catalyst is available. A bounded classical inverse
condition number does not make its quantum implementation free.

Separately improving the affine and reverse blocks remains sufficient,
and a linear-T endpoint circuit with larger fully charged Clifford cost
would still settle the T-only question. No joint interior with one global
precision charge has yet been constructed. The completed examples and
group audit supply native interfaces and fair comparison costs; they do
not show that the generic gap is close to resolution.

Optimal depth, practical constants, and a full elementary emitter are
separate tasks. The established publication scope is unchanged; manuscript
writing and release work are outside this pass.

## Evidence and remaining implementation work

Internal review of source conditioning, reverse-word order, workspace, and
error sums has not identified a defect; this is not independent peer
review. The [verification map](VERIFICATION.md) separates analytic proofs
from finite checks of native sources, literal phases, inverses, rejected
returns, support packing, and small complete dirty-input blocks. Tests do
not establish asymptotic theorems, optimality, or literature priority.

The strongest grouped construction has no end-to-end elementary Clifford+T
emitter. Its table fixture uses the table action directly, grouping tests
use illustrative constants, and some workspace constants and crossover
thresholds remain existential. The canonical audit adds a local register
table and explicit scalar/atom ordering; it does not supply the complete
emitter or practical crossover estimates.
Proof chapters retain the arguments; this checkpoint records the next decision.
