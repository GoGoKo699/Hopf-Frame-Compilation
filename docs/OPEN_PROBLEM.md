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

## Precision-uniform one-clean count

For $`n\ge3`$ and every integer $`L\ge6`$, the
[one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives full-isometry error at most $`2^{-L}`$ and

```math
a=1,\qquad T=O(N+L),\qquad G=O(NL).
```

The base dirty reservation is $`B=L+n+7`$ for real frames and
$`B=L+n+8`$ for phase-dressed complex magnitude frames; every
$`b\ge B`$ is allowed. At exactly $`b=B`$ and
$`L\ge\max\{6,n\}`$, the existing worst-case lower bound matches:
$`T^\star=\Theta(N+L)`$. For $`L\lt n`$ the base-width statement
is an upper bound.

A borrowed signal implements the two conjugate source sectors and returns
within the full-operator error because their recovered real target agrees.
The existing one-clean baseline covers modest precision and bounded n.
Literal leaf phases then use the sole clean compiler flag. All source
construction, queries, actual inverses, and work return are charged.

## Optimal banked count

For the same two frame families, $`n\ge3`$, every integer $`L\ge6`$,
and $`b\ge2B`$, the
[banked frontier](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
is

```math
a=1,\qquad
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O(NL).
```

Exact word-bank queries replace the growing sources' coefficient masks;
their weighted row and precision costs are summable. At modest precision,
the existing grouped compiler already attains this count. Banks and
selectors return exactly on arbitrary inputs, while source and signal
return are included in the common error norm. This is a count theorem
on its displayed dirty-width range, with no new depth guarantee.
The [earlier two-clean proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
retains its independently established allocation and proof.

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
initialized-isometry norm. The new one-clean theorem also gives linear
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
circuit. The selected one-clean count frontier also charges
$`G=O(NL)`$ elementary Clifford gates. Certified classical preprocessing
has no asserted uniform efficiency bound; the depth question retains its
displayed Clifford budget. Neither permits QRAM, a supplied precision state,
measurement, reset, postselection, or uncharged coherent evaluator.
A lower bound must apply to this full model; a failure of one source,
query interface, or factorization is a restriction of that method.
