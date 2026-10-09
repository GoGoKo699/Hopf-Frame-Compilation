# Two resource questions for complete Hopf frames

[Compiler claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) · [Research results](../research/README.md)

Write $`N=2^n`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
and $`\ell_*(n)=1+\log_2^*(n+2)`$. The coherent Clifford+T model
has $`a`$ initialized flags and $`b`$ arbitrary dirty qubits. Its full
initialized-isometry contract includes all logical inputs, dirty/reference
return, literal phases, and clean leakage. The supplied angles admit
terminating certified evaluation, including unequal and singular tuples.

## Constant-clean high-precision count

At the literal allocation

```math
a=2,\qquad \eta=2^{-N},\qquad L=N,\qquad b=N+n+7,\qquad n\ge3,
```

the worst-case real-frame bounds are

```math
2N-1\le T^\star\le O(N\ell_*(n)).
```

The [one-clean grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit)
gives the upper bound while leaving the second available flag unused.
The [explicit single-angle witness](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
at root angle $`\pi/16`$ gives the lower bound even at unrestricted
workspace width, retaining literal phases. The question
is whether the upper bound can be reduced to O(N) at this exact allocation.

### Bounded diagonal factorization

The [packed diagonal compiler](OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames)
gives error $`2^{-\ell}`$ using two flags and

```math
b_{\rm diag}=n+1+\left\lceil\frac{\ell+4}{2}\right\rceil,
\qquad T=O(N+\ell),\qquad G=O(N\ell).
```

Consequently, **any absolute constant number of diagonal factors** fits
the endpoint workspace with error left for certified factor selection.
The [factorization reduction](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
proves a terminating rational-grid search from exact coverage alone,
without requiring a continuous or computable exact phase-selection map.
It includes a finite-small-dimension fallback and propagation of existing
work leakage. The remaining hypothesis is coverage by a uniformly bounded
number of factors with charged interlayers.

For arbitrary diagonal masks and r Haar or adjoint occurrences, a
[global norm bound](../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md)
gives a common Hopf witness requiring

```math
[2\sqrt{n+1}]^{r+1}\ge(1-\epsilon)\sqrt N.
```

Thus even nonalternating Haar words need $`\Omega(n/\log n)`$ calls;
the five-mask $`Q_n^2`$ family does not give uniform coverage. Fixed
alphabets with polynomially bounded summed entrywise absolute norms obey
the same restriction. Dense Walsh mixers lie outside that bound. The
factorization chapter gives a stable, all-column nested-projector
criterion for a different mixer representation.

### A regular fixed-tree Cayley family

Every complete frame has the exact representation

```math
W=D_{\rm out}\Pi\sigma P(\phi)\sigma^\dagger D_{\rm in}^\dagger,
\qquad |\cot\phi_j|\lt\sqrt3.
```

The [tree reduction](../research/endpoint/TREE_CAYLEY_REDUCTION.md)
constructs both literal diagonals and both permutations effectively.
The permutations cost O(N) T gates exactly and return their dirty helpers.
The remaining operator acts on the same fixed incidence tree and is a
Cayley transform with a bounded real parameter table. Its variable
generator is a compressed difference of two fixed diagonal algebras.
The [boundary propagation construction](../research/endpoint/BOUNDARY_PROPAGATION.md)
gives an O(N)-T block encoding of $`(I-\mathcal B)/2`$ with two signal
ports. For every real table,
$`\sigma_{\min}(I-\mathcal B)=2\sin(\pi/(4n+2))`$.
Self-borrowed table programming and reuse of the occupied dilation port
allow R calls at precision $`q=N+h+9`$, where
$`h=\lceil\log_2R\rceil`$. Their dirty core uses $`q+1`$ wires and
their summed full-operator error is below $`43\eta/256`$ when
$`h\le n-3`$. Additional unchanged controls consume the corresponding
width slack. The direct T-count remains O(RN).

The complete feedback identity cancels an adjoint chain encoder on its
initialized entrance, leaving exactly one terminal prefix-frame encoder.
The [simultaneous coarse construction](../research/endpoint/COARSE_PREFIX_ENCODER.md)
implements all prefix frames in O(N) T gates. For the regular-core coins,
$`s=\lceil N/n\rceil`$, its full-operator error is below
$`(43/64)2^{-s}`$ and its dirty reservation is $`s+n+9\le N+n+7`$.
Conjugation around a middle gate within $`\kappa`$ of identity on the full space suppresses encoder
error to $`2\kappa\delta`$; the middle gate's native error is added
separately. A complete fine-precision correction with summed O(N) T-count
is the remaining resource question. The non-small chain reversal does not
itself receive this error suppression.

The [collective refinement theorem](../research/endpoint/COLLECTIVE_PRECISION_REFINEMENT.md)
constructs a replacement for the terminal regular-core frame with

```math
\epsilon\le40\,2^{-3s}+2^{-L},\qquad
T=O(N+3ns+n^4)+O(N+L),
```

for $`n\ge16`$, $`s=\lceil N/n\rceil`$, and $`3s\le L\le N`$,
within $`b=N+n+7`$. Exact geometric histories and midpoint conjugation
collect the local defects into one bank; opposite-sign purification
cancels its quadratic phase. A native baseline replacement accounts for
the original dirty disturbance. The complete initialized-output norm
includes literal phase, leakage, and dirty/reference return.

This fixed-order construction costs O(N), as does directly increasing
coarse precision by a fixed factor. It does not change the endpoint bounds.
Higher-order defects need not retain the rounded-coin form, and an
independent N-row program per quadratic stage would cost O(N log n).
The unresolved implication is a closed physical correction or another
global representation whose total source and table charges are O(N).

The [structural limits](../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md)
and [source interfaces](../research/endpoint/SOURCE_REUSE_LIMITS.md)
give precise restrictions on particular representations. They do not
strengthen the unrestricted lower bound beyond linear order.

## Large-width T-depth

At fixed $`0\lt\eta\le1/64`$, two clean flags, and sufficient
$`b=\Theta_\eta(\sqrt N)`$, the proved construction has

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

## Full-circuit requirement

A resolution must charge source construction, tables, basis changes,
queries, actual inverses, amplification, and workspace return on the same
circuit. The count question allows arbitrary finite Clifford work and
classical preprocessing. It allows no QRAM, supplied precision state,
measurement, reset, postselection, or uncharged coherent evaluator.
A lower bound must apply to this full model; a failure of one source,
query interface, or factorization is a restriction of that method.
