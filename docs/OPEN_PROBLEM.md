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
| Supporting depth results | [R80–R84](#large-width-t-depth) give restricted-family circuits and a general terminal improvement with one clean qubit and the count budgets; the [permutation checkpoint](../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md) leaves the arbitrary-orientation gap open |
| Constant-depth obstruction | The [zero-clean five-layer bound](#constant-depth-frontier) is proved; the unrestricted one-clean five-layer bound remains open |

The immediate lower-bound question is whether four T layers with one clean
qubit can approximate the explicit witness arbitrarily well when dirty
width, inner clean spread, and middle crossing rank all grow. The
[bounded-spread result](#one-clean-inner-spread) rules out fixed spread,
and its common-flip counterexample rules out that universal shortcut.
Neither settles the unrestricted question. All statements below concern
the complete initialized-isometry contract, not one state column.

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

The one-clean high-precision branch gives the stronger $`D_T=O(NL/b^2+L)`$ whenever $`L>N/n^2`$,
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
Optimal large-width depth remains open,
with unrestricted lower bound $`\Omega(1)`$ at square-root width.

The later [structured theorem R80](../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md)
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
improves only the last logarithmically many general layers; the early
decoder still determines the general square-logarithmic bound. The
[joint session record](../research/depth/JOINT_DEPTH_SESSION_RESULT.md)
identifies the remaining native mechanism and scoped failed routes.

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
uniform-reference case, and leaves the unrestricted frontier unchanged.

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
keeps the general frontier separate: arbitrary nonlinear orientations at
interior bias still require a native in-place permutation or a complete recovery word;
the factorization alone does not pay for either.

The [bounded in-place conversion session](../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
stopped at its mechanism checkpoint. Its exact word and encoding
restrictions apply to specified architectures, not unrestricted native
circuits. They supply no new general upper or lower bound. The sufficient
whole-register reduction and the remaining recovery/return obligation are
recorded there; arbitrary node-dependent magnitudes remain a broader
question.

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
