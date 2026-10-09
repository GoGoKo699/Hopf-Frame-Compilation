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
gives an O(N)-T block encoding of $`(I-\mathcal B)/2`$ using the two clean signals
and exactly $`N+n+7`$ dirty qubits. Its full-unitary error is below
$`(43/64)\eta`$. For every real parameter table,
$`\sigma_{\min}(I-\mathcal B)=2\sin(\pi/(4n+2))`$.
Completing the coupled Cayley feedback and returning the signals would
resolve this route; the block encoding alone is not the complete frame.

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
