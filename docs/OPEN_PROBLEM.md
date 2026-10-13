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

### Proved linear-count families

The following families have $`O(N)`$ T-count at the literal allocation
above. Constants may depend on fixed promises, never on $`n`$ or the
target tuple.

| Family | Exact promise and proof |
|---|---|
| Regular or common-scale branches | The [precision-gauge compiler](../research/endpoint/PRECISION_GAUGE_COMPILER.md) covers polynomially regular normalized branches, a fixed-ratio interval of branch exponents, and the stated cheap/ultraweak separation |
| Sparse computational forests | The [structured-forest compiler](../research/endpoint/STRUCTURED_FOREST_COMPILERS.md) covers fixed outgoing degree in the original tree $`p(j)=j-\mathrm{lowbit}(j)`$, more generally certified sub-binary growth within degree windows; active angles may be arbitrary, signed, or singular |
| Bounded exceptional path depth | The same forest theorem allows arbitrary-scale exceptional edges when each computational path meets only a fixed number, over a cheap or fixed-ratio base family |
| Native-baseline updates | [Antichain](../research/endpoint/ANTICHAIN_COMPILER.md) and [ancestor-closed sparse updates](../research/endpoint/SPARSE_UPDATE_COMPILER.md) retain their separate promises |

A common primitive is the [disjoint SU-block compiler](../research/endpoint/DISJOINT_BLOCK_COMPILER.md):
$`\sum_j m_j^2\le\kappa N`$ for fixed $`\kappa`$ gives
$`O_\kappa(N)`$ T-count, including the actual source construction,
inverses, and return error. Whole components can be too large for this
criterion; the forest theorem obtains a constant number of eligible
factors by exact component dressing.

The [rank-weighted compiler](../research/endpoint/RANK_WEIGHTED_BOUNDARY_COMPILER.md)
extends this primitive to supplied bounds
$`\mathrm{rank}(U_j-I)\le r_j`$ with $`\sum_jm_jr_j\le\kappa N`$.
It also compiles complete sequences of boundary corrections when their
collective basis and core satisfy the stated weighted budgets. Complex
blocks use two clean flags; real blocks use one through the
[borrowed conjugate signal](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md).
All retain the literal dirty reservation. The disjoint or collective-budget
promises do not follow from adding ranks of arbitrary overlapping factors.

These promises do not cover every supplied tuple. In particular, full
binomial components have binary descendant growth, while arbitrary
heterogeneous precision scales need not have a bounded precision sum.
Reordering a star into an arbitrary binary graph does not preserve the
original tree's degree-window condition.

### Bounded diagonal factorization

Any fixed number of literal diagonal factors with charged interlayers fits
the endpoint width. The [factorization reduction](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
proves terminating certified selection from exact coverage; uniformly
bounded coverage remains a hypothesis.

The [structural limits](../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md)
separate restrictions on Haar/flat-Clifford mixers from the general
problem. With fixed Clifford interleavings and differentiable diagonal
stages that are Clifford-valued at the origin, exact full-family synthesis
needs $`n`$ stages. An approximate version additionally requires a
physical Lipschitz bound. Non-Clifford baselines, nonsmooth selection,
and target-dependent interleavings are outside that theorem.

### A regular fixed-tree Cayley family

The [all-angle tree reduction](../research/endpoint/TREE_CAYLEY_REDUCTION.md)
expresses every frame through bounded cotangent parameters, two literal
diagonals, and exact $`O(N)`$-T permutations. The
[boundary feedback identity](../research/endpoint/BOUNDARY_PROPAGATION.md)
leaves one complete prefix encoder.
The [coarse encoder](../research/endpoint/COARSE_PREFIX_ENCODER.md) and
[fixed-order refinement](../research/endpoint/COLLECTIVE_PRECISION_REFINEMENT.md)
price coarse and cubic replacements in $`O(N)`$ T gates.
The [coupled inverse](../research/endpoint/COUPLED_TREE_RESOLVENT.md)
has an all-order norm bound and ideal scattering completion, but its
direct native precision cost remains $`O(N+nL)`$.

A complementary [canonical-star identity](../research/endpoint/CANONICAL_STAR_RESOLVENT.md)
gives a full native local scattering oracle at linear count and writes
the target using a triangular descendant inverse. The inverse is still
unpaid. Conditioning, small rank, and accepted-block identities do not
supply a complete native recovery circuit.

The [original-tree residual](../research/endpoint/ORIGINAL_TREE_RESIDUAL.md)
provides a complementary reduction without a change of tree. Its signed
path and determinant formulas locate the direct Cayley singularities.
An exactly returned native coarse frame has $`O(N)`$ T-count and leaves
a residual smaller than $`2^{-\lceil\sqrt N\rceil-6}`$, with rank at
most one across every original tree-edge cut. A regularized joint-Gram
construction completes its canonical factors on exactly $`4N`$ modes,
with bounded conditioning and cut ranks. This removes the earlier specific
completion overrun, not the unproved native cost of its fine source and
actual inverse.

The [ordered-bank and finite-history results](../research/endpoint/ORDERED_BLOCK_BANKS.md)
give an exact full-frame factorization with $`(4+o(1))N`$ total squared
block size and $`(1+o(1))\log_2n`$ ordered factors. The stored banks,
right readout, coarse prefix bank and local fine SU(2) bank have charged
native implementations under their stated interfaces. Their ordered fine
evaluation or the alternative fine left source remains unconstructed.
The sharp block-factorization and independently compiled-module bounds
restrict those representations; they do not change the unrestricted
T-count comparison.

The [scalar Schur theorem](../research/endpoint/SCALAR_SCHUR_LIMITS.md)
gives an exact complex coarse-frame gauge and a support-forced completion
rank of at least $`n+1`$ when the selected chain correlations satisfy
$`q_j\lt1`$. Its robust polynomial error floor concerns one
prescribed coarse call with fixed source-mode dimension and bounded,
logical-label-preserving sources. The operator-valued extension makes
that floor independent of arbitrary dirty width. Logical-label mixing,
multiple coarse calls, and other recovery mechanisms remain eligible.
The [native representation results](../research/endpoint/NATIVE_REALIFICATION.md)
separate linear ring height from native factorization and give a phase
cost for normalized scalar acceptance, even at arbitrarily small success
probability. The [single-echo bound](../research/endpoint/SOURCE_REUSE_LIMITS.md#whole-word-xor-cancellation-and-a-cyclic-alternative)
likewise applies only to its commuting involutory interface. These statements leave the
unrestricted endpoint bounds unchanged.

### Source-program resynthesis

The [source-program results](../research/endpoint/SOURCE_PROGRAM_RESYNTHESIS.md)
distinguish raw scalar-source programs from fully recovered groups.
After complete group recovery, a prefix-controlled dirty-core gauge that
commutes with the group target preserves its initialized-isometry error.
Different groups may therefore use independent gauges without physical
connectors. Raw prefix changes still require their connectors, and every
transformed group interior remains to be compiled.

At this endpoint, the same chapter gives a normal form using one fixed
physical core of $`\mu=N+6+\lceil\log_2n\rceil`$ dirty qubits for all
source calls, within the literal allocation and with two clean flags.
Width transitions disappear, but repeated precision charges retain the
$`O(N\ell_*(n))`$ bound. Resynthesis at fixed source algebra or fixed local
dimension does not price the growing interiors of the unrestricted family.

The [packed rejection-bank construction](../research/endpoint/PACKED_REJECTION_BANKS.md)
fits supported failure paths into the two mode flags. An input preweight
and one final recovery remove separate group recoveries, but its direct
fine-source cost remains O(nN). Its depthwise balanced schedule gives exact
four-port history compression and a full rounded-source comparison within
the literal dirty width. The preceding fine prefix-unwinding bank remains
unconstructed; a bank-preserving scalarization must approximate that prefix
on every logical and dirty input. The
[joint-word formulation](../research/endpoint/JOINT_NATIVE_WORDS.md) also
compresses returned dirty work according to the Pauli axes' symplectic rank,
without additional T gates or error. Neither result supplies a uniform
O(N) native selection rule or improves the endpoint comparison.

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
