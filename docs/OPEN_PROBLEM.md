<a id="two-resource-questions-for-complete-hopf-frames"></a>

<a id="endpoint-theorem-and-remaining-depth-question"></a>

<a id="two-clean-count-theorems-and-remaining-depth-question"></a>

<a id="one-clean-count-frontier-and-remaining-depth-question"></a>

# One-clean count frontier and sublinear large-width depth

[Compiler claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) · [Research results](../research/README.md)

Write $`N=2^n`$ and $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
The coherent Clifford+T model
has $`a`$ initialized flags and $`b`$ arbitrary dirty qubits. Its full
initialized-isometry contract includes all logical inputs, dirty/reference
return, literal phases, and clean leakage. The supplied angles and leaf phases
admit terminating certified evaluation, including unequal and singular tuples.

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

All upper constructions have $`G=O(NL)`$. The
[banked proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
charges exact word-bank queries, growing sources, coarse interlayers, actual
inverses, and all workspace return. A borrowed signal carries conjugate
source sectors while the genuine initialized recovery anchor remains.
The layerwise baseline covers low precision and bounded n; literal leaf
phases reuse the sole clean flag afterward. Banks and selectors return
exactly on arbitrary inputs, while source and signal return are included
in the full-isometry norm.

The count frontier is proved on the displayed banked range. At base width
and $`L\lt n`$, or throughout the transition $`B\lt b\lt2B`$,
retain the stated upper bound without a general matching claim. The
[publication scope](../manuscript/PUBLICATION_SCOPE.md#result-c-the-one-clean-count-theorem)
records the complete hypotheses. Result B retains its sufficient-clean
allocations, and Result D retains its distinct same-circuit depth guarantees.

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

The [growing-source compiler](../research/endpoint/GROWING_SOURCE_COMPILER.md)
gives the upper bound for every prescribed real frame, including unequal
and singular tuples. It charges source preparation, retained rejection,
actual inverses, and complete workspace return. The
[complex-magnitude corollary](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
composes a packed literal leaf-phase diagonal on the same flags and dirty
pool. Both factors' errors, including their work return, enter the complete
initialized-isometry norm. Result C also gives linear
count for real frames at this dirty allocation. The [explicit single-angle witness](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
at root angle $`\pi/16`$ gives the lower bound even at unrestricted
workspace width, retaining literal phases. Zero leaf phases include this
real-family witness in the complex magnitude family. Both two-clean count
endpoints are resolved at the displayed allocation. The complex endpoint
saves one dirty qubit relative to the uniform corollary. These are prescribed
frames; the one-clean complex theorem uses one additional dirty qubit.
Arbitrary complex-unitary synthesis and T-depth retain their separate scopes.

## Large-width T-depth

The [one-clean theorem](../research/depth/ONE_CLEAN_DEPTH.md), R72,
applies to both families for every $`n\ge1`$, integer $`L\ge6`$,
and one clean flag.
Put $`B_D=L+n+7`$ for complete real frames and $`B_D=L+n+8`$
for phase-dressed complex magnitude frames $`D_\phi W_{\mathbb R}`$.
At every $`b\ge17B_D`$, define

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
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
with the same count, Clifford, and error bounds. R72 supplies native
one-clean layer schedules, protected terminal flag relocation, and its
own width/precision splice without changing the literal dirty threshold.
Its complex extension proves a one-clean literal-diagonal schedule with
$`D_T=O(NL/b^2+L+\log_2(n+2))`$, then composes half-error factors
on the same flag and dirty pool. Its stated threshold includes that
precision increment; the high-precision bound applies to both families.
Every source preparation, query, actual inverse, and work-return error
is charged. Classical preprocessing terminates, with no efficiency claim.

The retained two-clean
[compressed-column special case](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains, for some absolute $`c>0`$, the stronger real-frame
$`\log_2(n+2)`$ dependence on
$`6\le L\le cn`$:

```math
D_T=O\!\left(NL/b^2+\min\{nL,(\sqrt n+L)\log_2(n+2)\}\right).
```

At sufficient $`b=\Theta(\sqrt{NL})`$, precision $`L=O(\sqrt n)`$
retains $`O(\sqrt n\log(n+2))`$ depth. In particular,
$`L=\lfloor\sqrt n\rfloor`$ has that depth and error at most
$`2^{-\lfloor\sqrt n\rfloor}`$, for sufficiently large n.
More generally $`L=o(n/\log(n+2))`$ gives sublinear depth.
Its early phase-gradient source and terminal binary source are separately
charged. The small absolute c is not a practical crossover estimate.

The one-clean fixed-accuracy consequence of R72 covers every $`n\ge1`$.
For every fixed $`0\lt\eta\le1/64`$, put
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. At one clean
flag and every $`b\ge17(L+n+7)`$, it gives

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(N/b^2+\min\{n,\sqrt n\log_2(n+2)\}\right).
```

At sufficient $`b=\Theta_\eta(\sqrt N)`$, this is optimal-order
T-count with $`O_\eta(\sqrt n\log(n+2))=o(n)`$ depth. Both count
and depth match through
$`b\le\sqrt N/[n^{1/4}\sqrt{\log_2(n+2)}]`$, whenever the
interval above the literal lower threshold is nonempty. The sublinear-depth
existence question is resolved; optimal large-width depth remains open,
with unrestricted lower bound $`\Omega(1)`$ at square-root width.

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
shows that the current source architecture needs worst-case normalization
$`\Omega(\sqrt n)`$, even at fixed operator error. That is a bound
on this encoding method, with no unrestricted native-depth consequence.
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

### Proved linear-count families

The [native primitives and promised families](../research/README.md#native-primitives-and-promised-families)
include disjoint and rank-weighted blocks, precision gauges, heterogeneous
sources, structured forests, and sparse updates. Several are dependencies
of the growing-source proof; their separate clean/dirty reservations and
stronger promises remain explicit in their own statements.

### Bounded diagonal factorization

[Bounded diagonal coverage](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
is a sufficient route with terminating certified selection. The
[structural limits](../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md)
apply to specified mixers and regularity assumptions, not arbitrary compilers.

### A regular fixed-tree Cayley family

The [full-frame reductions](../research/README.md#full-frame-reductions-and-remaining-native-steps)
preserve Cayley, prefix, scattering, residual, and ordered-bank constructions.
Each distinguishes its priced native operations from its remaining recovery
or fine-source step. Small rank, bounded conditioning, and compact storage
alone do not certify native cost.

### Source-program resynthesis

The [source programs and encodings](../research/README.md#source-programs-and-encodings)
and [diagonal and reflection architectures](../research/README.md#diagonal-and-reflection-architectures)
record paid fusion rules and restrictions on particular source, history,
and query interfaces. They supply neither an additive unrestricted lower
bound nor an uncharged native implementation of a different representation.

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
