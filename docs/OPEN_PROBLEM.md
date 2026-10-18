<a id="two-resource-questions-for-complete-hopf-frames"></a>

<a id="endpoint-theorem-and-remaining-depth-question"></a>

<a id="two-clean-count-theorems-and-remaining-depth-question"></a>

<a id="one-clean-count-frontier-and-remaining-depth-question"></a>

<a id="one-clean-count-frontier-and-sublinear-large-width-depth"></a>

# Resource frontiers and the remaining depth question

[Compiler claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) · [Research results](../research/README.md)

Write $`N=2^n`$ and $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
The coherent Clifford+T model
has $`a`$ initialized flags and $`b`$ arbitrary dirty qubits. Its full
initialized-isometry contract includes all logical inputs, dirty/reference
return, literal phases, and clean leakage. The supplied angles and leaf phases
admit terminating certified evaluation, including unequal and singular tuples.

## Current status

| Question | Status and reading route |
|---|---|
| Count at the selected clean/dirty allocations | Resolved by [Results B–C](../manuscript/PUBLICATION_SCOPE.md#four-principal-results); the exact-base and banked ranges remain distinct |
| Count and depth on one circuit | [Result D](#large-width-t-depth), R79, gives absolute all-precision bounds and polylogarithmic depth through $`L=\Theta(\log n)`$ with one clean qubit for both frame families; optimal large-width depth remains open |
| Fixed-accuracy general depth | [R90 dirty catalogue terminal](../research/depth/DIRTY_CATALOGUE_TERMINAL.md) gives $`O_L(N/b^2+\ell)`$, with $`\ell=\log_2(n+2)`$, for arbitrary magnitudes, retaining one clean qubit, the original literal widths, and the count budgets; optimal depth remains open |
| Stronger all-precision research bound | [R96](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md) gives $`O(NL/b^2+L+\ell)`$ for every $`n\ge1,L\ge6`$, with optimal counts, one clean qubit, original literal widths and full error; [R95](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md), [R94](../research/depth/PORTABLE_WORKSPACE_CATALOGUE.md) and [R92](../research/depth/UNIFORM_CATALOGUE_PRECISION.md) retain their independent narrower constructions |
| Supporting depth results | [R80–R85](#large-width-t-depth) retain promised-family and terminal theorems; R87's fine decoder, [R88's sharp-query/fallback interfaces](../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md), and [R89's coarse reductions](../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md) support R90 |
| Constant-depth obstruction | The [zero-clean five-layer bound](#constant-depth-frontier) is proved; [R97](#outer-axis-overlap) excludes bounded outer-axis misalignment; the unrestricted one-clean five-layer bound remains open |
| Query reformulation | [Conditional R98](../research/depth/QUERY_DEPTH_FRAME_TRANSFER.md) transfers a faster exact read-only dirty word-query interface to general fixed-accuracy frame depth; [two-call corner reductions](../research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md) separately transfer unrestricted Boolean-query lower bounds |
| Scoped returned-query obstructions | [Affine Gaussian scalar loops](../research/depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md) have Boolean degree at most $`4r+2`$ after r rounds; [two read-only T layers](../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md) force degree at most two and quadratic polar rank at most two with one clean qubit. Neither is an unrestricted dimension-growing T-depth lower bound |
| Three-layer query frontier | [Three read-only T layers](../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md) force affine XOR outputs with zero clean qubits and degree at most four with one clean qubit, at error below $`\sqrt2/6`$ independently of dirty width. [Exact one-clean two-layer queries are affine](../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-one-clean-queries-are-affine). Native AND depth is four with zero clean and three with one clean; the one-clean optimality is exact-only. This does not settle general depth |
| Four-layer zero-clean query closure | [Four read-only T layers](../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md#8-four-layers-with-zero-clean-qubits-degree-at-most-two) force XOR degree at most two at error below $`\sqrt2/12`$, independently of dirty width; literal Boolean scalar phases have degree at most three. The proof handles both visible and vanishing outer symplectic contributions; it gives no one-clean four-layer or unrestricted-address conclusion |

The preceding [native-defect round](../research/depth/NATIVE_DEFECT_SESSION_RESULT.md)
resolves the general **zero-clean read-only five-layer XOR degree** question:
[Section 5.6](../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
proves degree at most **14** at full operator error below $`\sqrt2/24`$,
independently of dirty width. An invisible outer target shear gives degree
at most 11. No comparison-span, axis-permutation, Pauli-comparison or
Clifford-middle hypothesis is added. A native Pauli pullback makes one-band
matching applicable to the ordered defect; its bounded tableau and sign
degrees then detect the Boolean output. The sharp degree remains open.

The earlier subcases retain their larger error radius $`\sqrt2/12`$:
Clifford-valued middle gives degree seven (four for a visible outer shear),
Pauli-valued comparison gives degree two, and individual-axis permutation
gives degree 244 (41 for an invisible target shear). The span-preserving
112/21 bounds are superseded by the new general theorem at their radii.
The [shared-prefix invariant](../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#55-shared-prefix-comparisons-a-robust-quadratic-group-map)
still supplies class-two, exponent-four comparison structure and faithful
degree-four intrinsic features below $`\sqrt2/40`$. Target extraction no
longer requires putting the outer baseline into that group. A
[returned affine counterexample](../research/depth/BASELINE_WORD_OBSTRUCTION.md#5-a-returned-affine-query-can-have-an-odd-order-outer-baseline)
shows why such an ambient-group inference would fail.

The [scalar representation construction](../research/depth/SCALAR_REPRESENTATION_PRODUCT.md)
still realizes arbitrary Boolean phases in four abstract factors. Its
Clifford-translation source mechanism needs at least $`2^{d-1}`$ carrier
qubits and T gates per full source for a degree-d table, even after carrier
compression and at the stated small full-source error. This is not a lower
bound for fused whole-word synthesis. Fusion removes irrelevant carrier
motion but leaves the table-bit loading task; the retained implementation
uses R91's logarithmic-depth dirty indicators. No improved reader is proved.

The [preceding checkpoint](../research/depth/SIX_LAYER_TRANSPORT_SESSION_RESULT.md)
adds a [six-layer Pauli-comparison theorem](../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md):
degree at most 51 below $`\sqrt2/64`$, independently of dirty width.
It does not cover general Clifford-valued comparisons. An actual fully
returned six-slot identity circuit has odd-order comparisons and persistent
non-Pauli multiplicative defects, ruling out that proposed recurrence even
under full return. Merging its empty slots changes the specified comparison;
its fair circuit baseline is identity.

A separate [diagonal-source theorem](../research/depth/RETURNED_DIAGONAL_ECHO_DEGREE.md)
bounds Boolean scalar degree by $`\max(2,q-1)`$ for q arbitrary diagonal
work-source occurrences joined by read-only monomial Clifford intervals,
at error below $`2\sin(\pi/2^{\max(q,3)+1})`$. Every inverse occurrence
counts. The four-call echo radius improves to $`2\sin(\pi/32)`$.
This is a source-architecture bound, with a radius decreasing in q; it is
not a general native T-depth lower bound or a faster reader.

The [diagonal-source transport theorem](../research/depth/NATIVE_CLIFFORD_TRANSPORT.md)
gives six-layer degree at most 78, or 26 for visible outer target shear,
below $`\sqrt2/8`$, when the specified prefix is diagonal C3 in fixed
Clifford coordinates. The inner Clifford need not be monomial. C3
membership alone only rounds nearby endpoints and cannot bound tableau
degree; the charged known permutation-family stress test is not a
three-band returned query.

The [preceding checkpoint](../research/depth/TWO_OVERLAP_SESSION_RESULT.md)
extended the source certificate to overlap two. The
[overlap-two theorem](../research/depth/TWO_OVERLAP_SOURCE_STRUCTURE.md)
proves that an exact three-band C3 prefix with first/middle span overlap
at most two is diagonal C3 in fixed Clifford coordinates, even when
neither adjacent partial source is C3. It therefore certifies the existing
78/26 bounds. The [earlier source theorem](../research/depth/NATIVE_BAND_SOURCE_STRUCTURE.md)
still confines middle/outer noncommutation to at most 2t middle axes for
general overlap t. Full-space involution arguments resolve the paired
cores and zero-character cases; raw outer-axis spread can remain unbounded
through duplicate cancellation. C3 membership is an extra hypothesis,
and overlap is a property of the specified physical representation.

A separate [hybrid transport theorem](../research/depth/HYBRID_DIAGONAL_CORE_TRANSPORT.md)
handles either ordering of a wide diagonal C3 source and a fixed a-qubit
core, when the full source is C3. With
$`N(a)=3\,2^{5a}+15a+2`$ and $`B(a)=24N(a)+8a^2+12a`$, it gives
label degree $`4B(a)+2`$ and literal sign degree $`12B(a)+6`$ below
$`\sqrt2/8`$, independently of dirty width and source coefficients.
Zero- and one-qubit cores retain the sharper 26/78 bounds. A commuting
T-gadget list with r binary relations following a diagonal-C3 source has
a hybrid core of at most $`4(2^r-1)`$ qubits if their product is C3.
These explicit source hypotheses are not consequences of return.

The [preceding round](../research/depth/COUPLED_INCIDENCE_SESSION_RESULT.md)
extends the [native source theorem](../research/depth/ISOLATED_NATIVE_CORE_STRUCTURE.md)
to arbitrary overlapping nonzero-character incidence. If every
zero-character B axis commutes with C, exact full-space C3 forces each
interacting class to consist of one isolated two-axis B/C pair. The
full involution defect localizes each isolated A probe; native A
isotropy forces its B-neighborhood to be uniform on each activation
class. This closes coupled wider rows without bounding total overlap
and certifies the existing 78/26 transport bounds. The
[current round](../research/depth/ZERO_CHARACTER_SESSION_RESULT.md)
adds two complementary certificates. The
[commutation-cover theorem](../research/depth/NATIVE_BAND_SOURCE_STRUCTURE.md)
closes every exact native three-band C3 source whose B axes each commute
with all of A or all of C, including arbitrarily many bad zero-character
axes. The [paired-orbit extension](../research/depth/ISOLATED_NATIVE_CORE_STRUCTURE.md)
proves that all zero-character axes are good if every singleton C orbit
is good. It uses character nonvanishing where the earlier diagonal Gram
bit can vanish. Both conclusions retain 78/26 transport at the same radius.

The separate [controlled-Pauli theorem](../research/depth/CONTROLLED_PAULI_TRANSPORT.md)
now assumes only an exact C3 source with scalar-Pauli computational
blocks in fixed input/output Clifford coordinates. A signed
normalization derives its canonical exact domain, removing the former
domain or fiber-compatibility premise. Label degree **58** and literal
sign/general target degree **174** hold below $`\sqrt2/8`$ independently
of control and target width. C3 derives quadratic labels and bounded
scalar phase degrees, retaining the weighted quartic term and all
Hermitian-section carries. The known non-semi-Clifford Gottesman–Mochon
source and its tensor powers satisfy the block hypothesis. Their paid
realization does not supply a native three-band source or a returned
query. The same proof home now gives a necessary-and-sufficient intrinsic
criterion: if E is the exact Pauli input domain and
$`N_p=S_{F\tau(p)F^\dagger}-I`$, scalar-Pauli coordinates exist exactly
when $`N_pN_q=0`$ for every $`p\in E^\perp`$ and every q.
This criterion includes literal signed normalization. A known C3
permutation fails it, so arbitrary C3 membership alone does not imply
the transport theorem. That source is not a native-three-band counterexample.

A separate [uniform transport theorem](../research/depth/BOUNDED_DEFICIENCY_TRANSPORT.md)
uses the intrinsic codimension k of the source's exact Pauli-preserving
domain, including degenerate domains and non-C3 sources. For k positive,
put $`D(k)=(2k+1)(4k+6)`$. It gives output degree at most $`3D(k)`$,
or $`D(k)`$ for visible shear, at error below
$`1/[4\sqrt2(2^k+1)]`$. Its radius is uniform in the source. Native
depth does not bound k; a wide single T band already has growing k.

The [commutative character bound](../research/depth/NONLINEAR_CHARACTER_TRANSPORT.md)
now covers every finite-dimensional commutative nilpotent binary algebra:
degree-a affine-input return needs at least $`2^{a-1}`$ coefficient
wires. A fully paid [noncommutative cycle construction](../research/depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md)
escapes that bound, realizing degree-d AND on $`d^2`$ coefficient
wires. Its exact native schedule uses $`O(d^3)`$ T and Clifford gates,
depth at most $`64(d-1)`$, $`d^2+d`$ dirty qubits and no clean work.
The retained conjunction circuits are faster. Thus the width obstruction
is not representation-independent, but no improved general reader follows.

The cycle source has an exact reformulation as the complete word query
for every entry of a binary unitriangular inverse, with only Hadamard
basis changes on the lower-triangular coefficient register. Its
controlled-source wrapper prices both source calls and all helpers.
The new commuting-pivot identity implements the full in-place inverse
with $`8\binom d3`$ T gates, $`20\binom d3`$ elementary Cliffords,
T-depth at most $`8\max(d-2,0)`$, no clean work, and at most
$`\lfloor(d-1)^2/4\rfloor`$ additional arbitrary helpers for
$`d\ge2`$. All helpers return exactly. Overlapping gates prevent
inferring a parallel schedule from whole-pivot commutativity. The circuit
changes its matrix inputs internally and gives no improved native
read-only query or complete frame.

The [balanced construction](../research/depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md)
now permits helper garbage during inverse computation and returns it
with the actual inverse after readout. For a d-by-d unitriangular input,
put $`q=\binom d2`$, $`h=\max\{\lceil\log_2d\rceil-1,0\}`$,
$`a=\lfloor d/2\rfloor`$, and $`c=\lceil d/2\rceil`$.
The exact complete inverse-word XOR query has

```math
\begin{aligned}
T&=64\binom d3,&G&\le160\binom d3+23q+4d,\\
D_T&\le128h,&W&=ac(c+3)/2.
\end{aligned}
```

Here W is additional arbitrary work beyond the q coefficient inputs
and $`d(d+1)/2`$ arbitrary output bits; d equal to one uses W zero.
All work and coefficient inputs return, no clean bit is used, and the
error is zero. Its actual reversed word and direct controlled variant
have separate complete ledgers in the proof home. The phase source has
the same T/depth/work orders. This improves the completed structured
source/query from linear to logarithmic depth at cubic rather than
quadratic work. It avoids the former repeated-child-inversion cost;
the determinant obstruction remains valid for its narrower returned
in-place, frozen-coefficient interface. The matrix inputs change
internally, so native read-only access and an arbitrary-table encoding
within R98's budgets remain missing. The [representation audit](../research/depth/GAUSS_SUM_QUERY_ENCODING.md)
clarifies that inverse-entry encodings already exist, but generic
substitution has excessive count/work and affine dirty loading leaves
offsets. Its completed-oracle counting boundary is scoped to affine
wiring and the retained generic circuit costs. An exact finite-field
interpolation repair returns all work with native read-only access,
but repeated evaluations lose to the direct path-propagation baseline.
Sparse specialization, coupled-mask extraction and quantum interference
remain outside those boundaries. The complete cycle conjunction's
separate transport and fair baseline are unchanged.

The next lower-bound question remains a **native degree-versus-depth
recurrence**. For exact three-native-band C3 sources, the immediate
remaining structural issue contains both a bad singleton C axis and
a middle B axis interacting with both outer bands. Does the native
factorization force $`N_pN_q=0`$ for $`p\in E^\perp`$ and arbitrary q?
The intrinsic criterion makes this sufficient for 174/58 transport
without requiring diagonalization. A representation-dependent parameter is
not an unrestricted depth invariant. General non-C3 prefixes, robust
radii for growing depth, one-clean work and address mixing remain
additional gaps. None of these results changes the complete-frame upper law.
A robust growing **unrestricted native** query lower bound
would transfer to a legal deepest-layer real-frame subfamily with two
calls, the same clean/dirty allocation, at most doubled error, T-count and
T-depth, and elementary Clifford budget $`G(V)+G(V^\dagger)+2\le4G(V)+2`$.
The actual inverse is charged in the fixed H, S, CNOT alphabet.
A faster exact multi-output query with R91's native read-only
interface would, by conditional R98, give

```math
D_T=O_L\!\left(N/b^2+d(2n+10)+\log_2\log_2(n+4)\right)
```

for both complete frame families at their original one-clean thresholds.
Here d is the hypothetical query's additive address-depth bound.
R91 supplies logarithmic d; no smaller d has been proved. These upper
and lower transfers have distinct query hypotheses. The
[Gaussian inverse repair](../research/depth/GAUSSIAN_INVERSE_EXTRACTION.md)
gives a paid identity and an affine-round obstruction, not an unrestricted
query or frame depth lower bound. The new
[scalar-loop proof](../research/depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
also detects the central sign, closing bounded-round extraction by actual
affine quadratic phase primitives. General native work interference is
outside that architecture. The
[multiway cost audit](../research/depth/MULTIWAY_QUERY_COST_BARRIER.md)
rules out a termwise separated-reader shortcut at the square-root count
envelope; it does not lower-bound joint native synthesis.
R96's unconditional frontier is unchanged.

The separate constant-layer question is whether four T layers with one clean
qubit can approximate the explicit witness arbitrarily well when dirty
width, inner clean spread, middle crossing rank and both outer-axis
misalignment dimensions all grow. The
[bounded-spread result](#one-clean-inner-spread) rules out fixed spread;
[R97](#outer-axis-overlap) rules out bounded misalignment on either side,
even at unbounded outer-intersection dimension. All statements below
concern the complete initialized-isometry contract, not one state column.

The [preparation map](../research/depth/DEPTH_OBSTACLE_MAP.md) separates
that precision-dependent obstruction from the dimension-growing depth
question. It reconstructs the pre-R78 upper costs and lower proof gap;
the present spectator-padded witness itself has constant depth at fixed
accuracy by the existing one-qubit compiler.

<a id="precision-uniform-two-clean-count"></a>

<a id="precision-uniform-one-clean-count"></a>
<a id="optimal-banked-count"></a>

## Result C: the one-clean count frontier

For every $`n\ge1`$ and integer $`L\ge6`$, the
[one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives error at most $`2^{-L}`$, with one clean qubit and universal
constants. Its base dirty reservation is $`B=L+n+7`$ for real frames
or $`B=L+n+8`$ for phase-dressed complex magnitude frames:

| Dirty allocation | T-count guarantee |
|---|---|
| Every $`b\ge B`$ | $`T=O(N+L)`$ |
| Exactly $`b=B`$, with $`L\ge\max\{6,n\}`$ | Worst-case $`T^\star=\Theta(N+L)`$ |
| Every $`b\ge2B`$ | Worst-case $`T^\star=\Theta(\sqrt{NL}+L+NL/b)`$ |

All upper constructions have $`G=O(NL)`$ and charge sources, queries,
actual inverses, and work return. The
[publication scope](../manuscript/PUBLICATION_SCOPE.md#result-c-the-one-clean-count-theorem)
is the complete statement; the [proof map](README.md#c-one-clean-compilation-and-its-corollaries)
gives its dependencies. At base width with $`L\lt n`$, or in the
transition $`B\lt b\lt2B`$, retain the upper bound without a general
matching claim. These count results add no depth guarantee.

## Constant-clean high-precision count

At the literal allocation

```math
a=2,\qquad \eta=2^{-N},\qquad L=N,\qquad b=N+n+7,\qquad n\ge3,
```

the selected worst-case count for each of the complete real and phase-dressed
complex magnitude families is

```math
2N-1\le T^\star\le O(N),\qquad T^\star=\Theta(N).
```

The [growing-source construction and literal phase composition](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
give the upper bounds under the full return contract. The
[literal-phase root witness](FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
gives the lower bound even at unrestricted width. This allocation saves
one dirty qubit for the complex family relative to C; the one-clean
complex theorem uses $`N+n+8`$ dirty qubits at $`L=N`$.

## Large-width T-depth

The [uniform subtree-depth theorem](../research/depth/UNIFORM_SUBTREE_DEPTH.md), R79,
applies to both families for every $`n\ge1`$, integer $`L\ge6`$,
and one clean flag.
Put $`B_D=L+n+7`$ for complete real frames and $`B_D=L+n+8`$
for phase-dressed complex magnitude frames $`D_\phi W_{\mathbb R}`$.
At every $`b\ge17B_D`$, use the following bound.

With $`\ell=\log_2(n+2)`$ and $`\lambda=\log_2(n+L+2)`$, put

```math
F(n,L)=\min\{nL,(L+\ell)\lambda\}.
```

One circuit has error at most $`2^{-L}`$ in the full initialized-isometry
norm and simultaneous bounds with absolute constants:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

Count is worst-case optimal at every eligible precision and width.
The inherited lower bounds also match depth on every nonempty interval

```math
17B_D\le b\le\sqrt{NL/F(n,L)}:
\qquad T^\star=\Theta(NL/b),\quad D_T^\star=\Theta(NL/b^2).
```

The one-clean high-precision branch gives the stronger $`D_T=O(NL/b^2+L)`$ whenever $`L\gt N/n^2`$,
with the same count, Clifford, and error bounds. R79 charges the native
subtree decoder's precision-dependent source calls and all-work predicates.
Its protected suffix allocation absorbs their table volumes with absolute
constants. The retained [one-clean baseline, R72](../research/depth/ONE_CLEAN_DEPTH.md),
supplies the terminal flag relocation, width and precision fallbacks, and
literal phase composition at the original dirty thresholds.

The complex extension uses the one-clean literal-diagonal schedule with
$`D_T=O(NL/b^2+L+\log_2(n+2))`$, then composes half-error factors
on the same flag and dirty pool. Its stated threshold includes that
precision increment; the high-precision bound applies to both families.
Every source preparation, query, actual inverse, and work-return error
is charged. Classical preprocessing terminates, with no efficiency claim.

At sufficient $`b=\Theta(\sqrt{NL})`$ and $`L=\Theta(\log n)`$,
R79 gives $`O(\log^2 n)`$ depth with error inverse polynomial in n.
More generally, $`L=o(n/\log(n+2))`$ gives sublinear depth.
The earlier two-clean
[compressed-column special case](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains its original source and recovery proofs; R79 supersedes its
resource bound.

The [subtree-projector decoder](../research/depth/SUBTREE_PROJECTOR_DECODER.md),
R78, is the proof home for native subtree reflections and coherent spectral
label erasure. Its fixed-accuracy theorem is now a corollary of R79.
For every fixed $`0\lt\eta\le1/64`$, put
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. At one clean
flag, every $`b\ge17(L+n+7)`$ for real frames, and every
$`b\ge17(L+n+8)`$ for phase-dressed complex magnitude frames, it gives

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(N/b^2+\log_2^2(n+2)\right).
```

The same circuit has full initialized-isometry error at most eta.
At sufficient $`b=\Theta_\eta(\sqrt N)`$, this is optimal-order
T-count with $`O_\eta(\log^2(n+2))`$ depth. Both count
and depth match through
$`b\le\sqrt N/\log_2(n+2)`$, whenever the
interval above the family's literal lower threshold is nonempty: there
$`T^\star=\Theta_\eta(N/b)`$ and $`D_T^\star=\Theta_\eta(N/b^2)`$.
The [interval-dilation theorem R86](../research/depth/INTERVAL_DILATION_DEPTH.md)
first improved the fixed-accuracy general bound to an additive
$`O_L(\ell^{3/2})`$. The later
[product-coin theorem R87](../research/depth/COIN_DILATION_DEPTH.md)
removes its overlap-amplification factor. For every fixed integer
$`L\ge6`$, every n, and every original eligible real or complex dirty
allocation, one clean qubit gives

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell\log_2(\ell+2)),\qquad \ell=\log_2(n+2).
```

It covers arbitrary node-dependent magnitudes, signs, zeros, and the
prescribed leaf phases, with initialized-isometry error at most
$`2^{-L}`$. Retained global finite rounding makes a positive weighted
tree. A product-coin embedding followed by tensor Hadamards exposes
the fine marker depth through the last nonzero coin block. The paid
source and its actual inverse erase the tag and restore the work.
The [selected-routing lemma](../research/depth/SELECTED_DIRTY_ROUTING.md)
supplies its physical routing and variable-size query interfaces.
The subsequent [causal compact-coin theorem R89](../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md)
retains those counts, widths, and full real/complex contract, and improves
the general depth to

```math
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(n+2))\right).
```

Here $`\log_2^*x`$ counts repeated base-two logarithms until the value
is at most 2. Its compact fixed coins, causally batched local matchings,
and parallel reference transforms erase the coarse label in
$`O_L(\ell)`$ depth. Preserved-prefix recursion supplies the remaining
iterated-logarithmic factor. Every source, actual inverse, reference-work
zero test, and query is charged.

The later [dirty catalogue theorem R90](../research/depth/DIRTY_CATALOGUE_TERMINAL.md)
retains those counts, literal widths, and complete real/complex guarantees,
and removes the iterated-logarithmic factor:

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell).
```

Three R89 coarse reductions leave a small terminal. Its complete rounded
program selects local catalogue factors through an exact arbitrary-dirty
word query and XOR slot translations. A commutator using two disjoint
borrowed arrays cancels their unwanted actions exactly; the local factors
need no initialized ancilla. The selected-factor approximation, full
initialized-isometry error, and original-width splice are proved together.
At square-root dirty width, the general upper bound is now $`O_L(\ell)`$,
while the unrestricted lower bound remains $`\Omega(1)`$. Optimal depth
remains open. R79 retains its absolute all-precision scope; R90's constants
depend on fixed L. The finite catalogue implies no efficient classical
preprocessing, practical crossover, or scalable native emitter guarantee.

The [precision-terminal session](../research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
adds two absolute-constant frame refinements. Put

```math
F_{79}(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

[R92](../research/depth/UNIFORM_CATALOGUE_PRECISION.md) proves, for all
$`n\ge1,L\ge6`$ and an absolute $`C_*`$,

```math
D_T=O\!\left(NL/b^2+
\min\{F_{79}(n,L),\ell+2^{C_*(L+1)^3}\}\right).
```

In particular, logarithmic depth holds for the unbounded precision window
$`6\le L\le c[\log_2(\ell+2)]^{1/3}`$, with one absolute positive c.
The independent [grouped-catalogue theorem R93](../research/depth/GROUPED_CATALOGUE_PRECISION.md)
proves, on $`6\le L\le2\ell`$,

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)
                     +\frac{L\ell}{\log_2(\ell+2)}\right).
```

Both retain optimal $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$,
one clean qubit, the original family-specific literal widths, and the full
initialized-isometry contract. Select a complete circuit to take the best
available depth bound; R79's stronger high-precision clause is retained.
At $`L=\Theta(\ell)`$ and $`b=\Theta(\sqrt{NL})`$, R93 gives
$`O(\ell^2/\log_2(\ell+2))`$. The later
[portable-workspace theorem R94](../research/depth/PORTABLE_WORKSPACE_CATALOGUE.md)
proves the stronger uniform bound

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right),
\qquad 6\le L\le2\ell.
```

R94 has absolute constants and retains the same optimal T-count,
$`G=O(NL)`$, one external clean qubit, real width
$`b\ge17(L+n+7)`$, complex width $`b\ge17(L+n+8)`$, and full
initialized-isometry error $`2^{-L}`$. R79 compiles larger groups
directly to final accuracy. Carrying its entire local workspace through
the selected transfer makes the borrowed slot action cancel exactly,
without initialized residual labels. A constant number of groups uses
the [R91 parallel word query and paired router](../research/depth/PARALLEL_BILINEAR_WORD_QUERY.md).
R92 and R93 remain independent of R91. R94 reaches the
$`O(\ell\log_2(\ell+2))`$ square-root-width benchmark.

The [logarithmic-precision theorem R95](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
removes the remaining logarithmic factor:

```math
D_T=O(NL/b^2+L+\ell),\qquad 6\le L\le2\ell.
```

Its constants are absolute, and it preserves the optimal T-count,
Clifford budget, original real/complex literal widths, one initialized
qubit and full error contract. The
[coarse-flag construction](../research/depth/COARSE_FLAG_SPARSE_CORRECTION.md)
implements a coarse mathematical frame to final accuracy and adds a
constant-normalization sparse correction. The
[small-workspace source](../research/depth/SMALL_WORKSPACE_DESCENDANT_SOURCE.md)
supports a [local compiler](../research/depth/RELAXED_COUNT_LOCAL_DEPTH.md)
with depth linear in height plus precision and much larger dirty counts.
R94's complete-payload transfer absorbs those local costs globally.
Thus inverse-polynomial accuracy in n has $`O(\ell)`$ depth at
$`b=\Theta(\sqrt{NL})`$.

The subsequent [all-precision theorem R96](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
removes that range restriction:

```math
D_T=O(NL/b^2+L+\ell),\qquad n\ge1,\quad L\ge6.
```

It preserves $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, one external
initialized qubit, real width $`17(L+n+7)`$, complex width
$`17(L+n+8)`$, and full initialized-isometry error $`2^{-L}`$.
All constants are absolute. A globally priced borrowed-precision stage
and fast capped references give a constant number of height reductions.
A directly programmed high-precision terminal uses the original global
prefix in every fine mask and coarse instruction query. It does not
allocate a catalogue of all precision programs. At sufficient width
$`b=\Theta(\sqrt{NL}+L+n)`$, its depth is $`O(L+\ell)`$ for every
precision. No efficient preprocessing, practical crossover or scalable
native emitter is claimed. Optimal general depth remains open; the
selected manuscript Result D remains R79.

The [structured theorem R80](../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md)
attains logarithmic depth when angle magnitudes depend only on tree depth,
including arbitrary node signs and supplied leaf phases. At the same
fixed-accuracy reservations it gives

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(N/b^2+\log_2(n+2)).
```

The signed balanced real subfamily has an exact version. This is a
restricted-family result, not a replacement for R79. The separate
[recursive terminal lemma R81](../research/depth/RECURSIVE_TERMINAL_DEPTH.md)
now covers every terminal height $`0\le M\le n`$, with depth
$`O_L(N/b^2+\ell\log_2(M+2))`$, where $`\ell=\log_2(n+2)`$.
A logarithmic protected suffix contributes only log-times-log-log depth;
the unrestricted early decoder still determined the square-logarithmic
bound at that checkpoint. The
[joint session record](../research/depth/JOINT_DEPTH_SESSION_RESULT.md)
identifies the remaining native mechanism and scoped failed routes.

The later [prefix-selected terminal theorem R88](../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md)
improves that component for every $`0\le M\le n`$ to

```math
D_T=O_L\!\left(N/b^2+\ell\log_2\log_2(M+4)\right).
```

It compiles the deepest M layers of an arbitrary real frame, including
signed and singular tuples, while the unchanged external prefix may
select unrelated local frames. It retains one clean qubit,
$`b\ge17(L+n+7)`$, $`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and
full initialized-isometry error strictly below $`2^{-L}/4`$.
For $`M\le C_0\ell`$ with fixed $`C_0`$, the additive depth is
$`O_{L,C_0}(\ell\log_2\log_2(\ell+4))`$.
The [contiguous coin-layout router](../research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
separately has $`O(\ell)`$ depth, including arbitrary dirty padding and
occupied inactive work. These component gains initially left R87's
complete bound unchanged. R89 subsequently replaces its coarse search
and gives the stronger arbitrary-height real-terminal bound
$`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$. R88 remains an independent
weaker terminal theorem and the proof home for sharp-query and fallback
interfaces used by R89.

The [bounded-density theorem R82](../research/depth/BOUNDED_DENSITY_DEPTH.md)
also improves the complete circuit when root probabilities have fixed
positive upper and lower ratios to a classically specified positive
product or affine-oriented reference. The latter uses the priced R83
reference core, including its leaf and marker corrections. With fixed
accuracy and ratio bounds, the same count and literal width guarantees
hold, with

```math
D_T=O\!\left(N/b^2+\log_2(n+2)(1+\log^*(n+2))\right).
```

The constants may depend on those fixed bounds. A native weighted
correction replaces repeated early threshold comparisons; every
recursive stage still pays for its global-prefix queries. This covers
nonhomogeneous magnitudes under the stated promise, including the
uniform-reference case. It does not require or establish an unrestricted
general improvement. R89 later attains this depth order for arbitrary
magnitudes, without the density promise.

The [affine-orientation theorem R83](../research/depth/AFFINE_ORIENTATION_DEPTH.md)
extends R80's logarithmic depth and count orders to depth-dependent
unordered magnitude pairs whose branch choices are affine functions of
the node prefix. Node signs, zeros, and supplied leaf phases are allowed.
Both R83 and the R82 extension retain, for fixed $`L\ge6`$, the literal
widths $`b\ge17(L+n+7)`$ for real frames and $`b\ge17(L+n+8)`$ for
the prescribed phase-dressed complex completion, with error at most
$`2^{-L}`$. Their constants may depend on L, and R82's also on the
fixed density bounds; these are not new all-precision claims.

The [deterministic-frame theorem R84](../research/depth/DETERMINISTIC_FRAME_DEPTH.md)
instead permits arbitrary node-dependent signed quarter turns and is
exact. For every $`b\ge17(n+13)`$, it has one clean qubit,
$`T=O(\sqrt N+N/b)`$, $`G=O(N)`$, and the same displayed
log-times-iterated-log depth order as R82, with absolute constants. Its
classical output-depth query erases the tag because every column is a
signed basis leaf. The
[branch-orientation session](../research/depth/BRANCH_ORIENTATION_SESSION_RESULT.md)
kept the general frontier separate: the factorization alone did not pay
for arbitrary nonlinear orientation at interior bias.

The [geometric-band theorem R85](../research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md)
now covers arbitrary nonlinear branch choices when the unordered
magnitude pair depends only on depth. It includes arbitrary node signs,
zeros, balanced pairs, and supplied leaf phases, with the same fixed-L
literal real and complex widths as R83. Its simultaneous bounds are

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\bigl(N/b^2+\ell\log_2(\ell+2)\bigr).
```

The source's fine tag first selects a coarse geometric band. A
selected subtree unloader then exposes enough clean suffix to implement
the local orientation permutation through two returned-output XOR
queries. A product reference removes the fine label; only the coarse
band needs spectral search. The initialized-isometry proof includes
every column, literal phases, actual inverses, leakage, dirty inputs,
and references. This reaches sub-square-logarithmic depth for the
arbitrary-orientation promise, without a general whole-register
permutation compiler. R87 subsequently matched its asymptotic depth for
arbitrary magnitudes, and R89 improves the general bound further.
R90 reaches logarithmic depth for arbitrary magnitudes, including this
promised family; R85's independent proof remains valid.

The [bounded in-place conversion session](../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
stopped at its mechanism checkpoint. Its exact word and encoding
restrictions apply to specified architectures, not unrestricted native
circuits. They supply no new general upper or lower bound. The sufficient
whole-register reduction and the remaining recovery/return obligation are
recorded there. R86 subsequently obtains a general improvement through
interval dilation rather than that standalone permutation interface.

The [joint orientation session](../research/depth/JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
records R85's discovery and review. The later
[general-magnitude session](../research/depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
records R86 and its reserve routes. The
[factorization session](../research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
then proves R87 through product-coin spectral readout. The
[continuation](../research/depth/COIN_CONTINUATION_SESSION_RESULT.md)
proves R89 through compact coins, causal batching, and prefix recursion.
The [whole-terminal session](../research/depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
proves R90 through fixed-number coarse reduction and dirty catalogue
commutators. The [precision follow-up](../research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
then proves R91–R93. The [catalogue bootstrap](../research/depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md)
proves R94 and reaches its growing-precision benchmark. The
[R95 construction](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
then reaches logarithmic depth at inverse-polynomial accuracy through
coarse flags, sparse correction and local compilation with relaxed counts.
The [precision bridge](../research/depth/BORROWED_PRECISION_BRIDGE_SESSION_RESULT.md)
then proves R96's all-precision upper law through globally priced sources,
fast references and a directly programmed terminal. The
[continuing research mandate](../WORKSPACE.md#continuing-toward-the-general-answer)
retains optimal general depth and dimension-growing lower bounds as open
questions. All further
improvements must retain the full native return and resource contract.

## Constant-depth frontier

### Zero-clean error floors and fixed-clean thresholds

The [constant-depth obstruction](../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md),
R73, strengthens that constant. The complete frame
$`W_n=I_{2^{n-1}}\otimes R_y(\theta)`$ with
$`\theta=\tfrac12\arccos(1/3)`$ has exact zero-clean error floor

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\approx0.169714
```

at T-depth at most three, uniformly over n and dirty width b. A Clifford
attains this floor, so integer $`L\ge3`$ requires at least four layers.
For every fixed positive clean budget a, a positive threshold
$`\varepsilon_a`$ exists; R74's elementary refinement below makes it
explicit. These statements
allow arbitrary Clifford interlayers and dirty/reference inputs.
Separately, an explicit one-clean witness requires at least three layers
for integer $`L\ge8`$. Zero leaf phases put both witnesses in the complex
magnitude family. These depth lower bounds are constant.

The [four-layer operator-rigidity theorem](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md),
R75, strengthens the zero-clean obstruction without restricting crossing
rank. Distinct Pauli images under two-layer circuits have operator-norm
distance at least

```math
\rho=\frac{\sin(\pi/8)}{\sqrt2}\approx0.270598.
```

The same complete frame therefore has error at least
$`\rho/2\approx0.135299`$ at T-depth at most four, with zero clean qubits
and arbitrary dirty width. Since $`1/8\lt\rho/2`$, every integer
$`L\ge3`$ requires at least five layers. This lower bound does not claim
that the four-layer error floor is sharp. With one initialized clean
qubit, error below $`\rho/2`$ forces exact terminal clean-sector
preservation.

### One-clean midpoint restrictions

The [one-clean midpoint theorem](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md),
R76, gives two further error floors with arbitrary dirty width and
crossing rank. For a split $`V=V_2V_1`$ with each half of T-depth at
most two, let $`S=V_1Z_CV_1^\dagger`$ be the encoded clean observable.

| A four-layer realization admits | One-clean consequence for the same witness |
|---|---|
| A split with S Pauli | Error at least $`\tau=\rho/2`$; integer $`L\ge3`$ requires five layers in this subclass |
| A split with S Hermitian Clifford | Error at least an explicit positive $`\eta_C`$, independent of width and crossing rank; the proof gives a very conservative constant |

Terminal return does not imply either midpoint hypothesis. R76's
[remaining question](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md#6-the-remaining-one-clean-question)
is freezing the compressed conjugations for a non-Clifford midpoint;
commutant lifting and phase discreteness already extend to that case.
A separate spectral argument gives error at least $`1/196608`$ when
the exactly returned occupied block has finite projective order; in
particular, four-layer full circuits of finite projective order are
excluded at $`L\ge18`$. Generic Clifford+T circuits need not satisfy
that condition.

### One-clean inner spread

The [inner clean-spread theorem](../research/depth/ONE_CLEAN_INNER_SPREAD.md),
R77, extends the compressed rigidity argument to non-Clifford midpoints.
In either two-layer encoder, let t count the first layer's axes that
anticommute with the incoming clean Pauli after its initial Clifford.
For every fixed t there is an explicit positive error floor, independent
of dirty width and middle crossing rank. The proof expands the midpoint
into $`2^t`$ Clifford terms and applies R74's finite-boundary argument
without crossing phases. Define $`t_{\min}(V)`$ by minimizing over
both encoders and every native four-layer representation. At error
$`2^{-L}`$, its explicit inequality implies

```math
t_{\min}(V)=\Omega(\log\log L).
```

Thus the first and last T layers must each have diverging clean-active
size. Together with R74 below, all four layers must grow in any
hypothetical arbitrarily accurate one-clean sequence. This does not
exclude unbounded spread. An exact two-qubit example also shows that
the separate encoder flips cannot always be replaced by one common
sign-compatible flip, even allowing an arbitrary unitary flip. Its
occupied defect is a positive constant, so it does not refute the
remaining gap conjecture.

### Outer-axis overlap

The [outer-overlap theorem](../research/depth/OUTER_OVERLAP_RIGIDITY.md),
R97, gives a positive initialized-isometry error floor depending only on
$`r=\min\{\mu_A,\mu_B\}`$. Here $`I=L_A\cap L_B`$ is the
intersection of the two native outer axis spaces, and $`\mu_A`$ is
$`\dim I`$ minus the number of individual A axes lying in I; define
$`\mu_B`$ similarly. The physical lists are independent, so these are
nonnegative dimensions. Both inner clean spreads, crossing rank, dirty
width and the full intersection dimension may grow while this floor
remains fixed.

At $`r=0`$, the intersection is spanned by a subset of either physical
list and the floor is R76's uniform Clifford-midpoint threshold. Literal
matching between the lists is unnecessary. An arbitrarily accurate
four-layer sequence must have both misalignment dimensions diverge in
every admissible native representation. That unbounded case remains
open. These are constant-depth obstructions, not matching
dimension-growing lower bounds for R96.

The [representation limits](../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md)
retain two distinct conclusions. Exact clean return does not guarantee an
exact same-space phase-polynomial descent, even with one additional dirty
wire. Growing center mixing in the displayed R77 tensor family also does
not obstruct an alternative native split: that family has a
Pauli-preserving three-layer realization. Neither example refutes the
remaining occupied-sector gap.

### Middle crossing rank and the remaining gap

The [four-layer crossing-rank theorem](../research/depth/FOUR_LAYER_CROSSING_RANK.md),
R74, gives a positive error gap for the same $`1/3`$ witness for each
fixed clean budget a and rank bound R, for circuits admitting a four-layer
decomposition whose middle active symplectic crossing rank is at most R.
Its elementary initialized-observable argument strengthens this to an
explicit rigidity gap $`2^{-O_a(r^2)}`$. Consequently, at error
$`2^{-L}\lt c_0`$, every such approximation satisfies

```math
r_{\min}=\Omega_a(\sqrt L),
```

where the minimum is over all four-layer decompositions of the circuit.
The proof supplies explicit constants for each fixed clean budget and
also an explicit fixed-clean three-layer threshold. Its original
polynomial-limit closure proof remains valid. With positive clean
budget, this rank restriction permits four-layer circuits whose crossing
rank grows with precision; it is not an unrestricted five-layer lower
bound. The width-independent four-layer Hilbert--Schmidt gap remains
unresolved, separately from R75's operator-norm result.
Optimal depth remains an open problem.

The [entry-indexed encoding bound](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
shows that the entry-indexed overlap architecture needs worst-case normalization
$`\Omega(\sqrt n)`$, even at fixed operator error. That is a bound
on this encoding method, with no unrestricted native-depth consequence.
R78 uses native spectral-label erasure instead of that overlap recovery.
T-depth permits arbitrary Clifford interlayers and is distinct from total
physical depth. The [depth results](../research/README.md#large-workspace-t-depth)
record the new complete construction and preserve the earlier component
and interface results under their original hypotheses.

## Earlier endpoint constructions

The [research index](../research/README.md#constant-clean-high-precision-count)
is the catalogue of earlier constructions and their exact hypotheses.
Their method-specific open steps are alternatives to the resolved endpoint,
not gaps in the growing-source theorem. Unique proofs, negative results,
and executable evidence remain at their original paths.

<a id="proved-linear-count-families"></a>
<a id="bounded-diagonal-factorization"></a>
<a id="a-regular-fixed-tree-cayley-family"></a>
<a id="source-program-resynthesis"></a>

| Retained route | Catalogue and boundary |
|---|---|
| Linear-count promised families | [Native primitives and promised families](../research/README.md#native-primitives-and-promised-families): distinct promises and workspace reservations; several support C |
| Fixed-tree, Cayley, scattering, residual, and ordered-bank reductions | [Full-frame reductions](../research/README.md#full-frame-reductions-and-remaining-native-steps): compact representation and conditioning alone do not price native recovery |
| Source-program resynthesis | [Source programs and encodings](../research/README.md#source-programs-and-encodings): paid fusion rules and scoped history/query restrictions |
| Bounded diagonal and reflection factorizations | [Diagonal and reflection architectures](../research/README.md#diagonal-and-reflection-architectures): certified selection and restrictions retain their stated mixers and regularity assumptions |

## Full-circuit requirement

Every compiler must charge source construction, tables, basis changes,
queries, actual inverses, amplification, and workspace return on the same
circuit. Result C also charges
$`G=O(NL)`$ elementary Clifford gates. Certified classical preprocessing
has no asserted uniform efficiency bound; the depth question retains its
displayed Clifford budget. Neither permits QRAM, a supplied precision state,
measurement, reset, postselection, or uncharged coherent evaluator.
A lower bound must apply to this full model; a failure of one source,
query interface, or factorization is a restriction of that method.
