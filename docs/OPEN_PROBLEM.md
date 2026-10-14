<a id="two-resource-questions-for-complete-hopf-frames"></a>

<a id="endpoint-theorem-and-remaining-depth-question"></a>

<a id="two-clean-count-theorems-and-remaining-depth-question"></a>

# One-clean count frontier and remaining depth question

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

For complete real frames at fixed $`0\lt\eta\le1/64`$, two clean flags,
and sufficient $`b=\Theta_\eta(\sqrt N)`$, the proved construction has

```math
T=\Theta_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(n).
```

The unrestricted depth lower bound is $`\Omega(1)`$. The question is
whether one circuit can attain $`D_T=o(n)`$ with the same count, Clifford,
and workspace budgets. The explicit matching intervals at smaller widths
are those of [Result D](UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary).
T-depth permits arbitrary Clifford interlayers and is distinct from total
physical depth. The [depth results](../research/README.md#large-workspace-t-depth)
separate source, predicate, query, and logical-transport costs.

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
