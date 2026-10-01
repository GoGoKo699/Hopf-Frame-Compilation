# Publication scope

[Paper architecture](README.md) · [Read the argument](../REVIEW.md) · [Verification](../docs/VERIFICATION.md)

The selected publication is **one full-length theoretical compiler paper**.
Its scientific ingredients are developed and checked
in this repository before final manuscript writing.
Its subject is the prescribed Hopf differential frame, studied in two resource
models. The constant-clean endpoint is an open problem in the discussion;
resolving it is not a prerequisite for this paper.

Working title: **Exact and Fault-Tolerant Compilation of Hopf Differential Frames**.

## The central claim

Preparing a state fixes one column of a unitary. Hopf quantum backpropagation
uses a prescribed completion containing the state and its coordinate-frame
directions. This paper shows how to compile that stronger operator interface
with optimal exact size/depth and explicit fault-tolerant precision/workspace
tradeoffs. The QBP consequence explains why preserving the extra columns
matters operationally.

This is the connection between CNOT-based logical synthesis and T-count:
they price the same prescribed operator under different gate models. The
theorems use different constructions and do not assert that one circuit
jointly minimizes CNOT count, T-count, depth, and workspace.

## Three principal results

Write $`N=2^n`$, $`n\ge1`$. The exact model uses arbitrary one-qubit gates
and CNOTs with all-to-all connectivity. Its workspace budget is $`m`$ clean
qubits. The coherent Clifford+T model uses $`a`$ clean and $`b`$ dirty
qubits, $`q=n+a+b`$, and

```math
0\lt \eta\le1/64,\qquad
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad
h=1+\lceil\log_2(L+n+2)\rceil.
```

Upper bounds hold for every prescribed parameter tuple. Matching lower bounds
are worst-case statements over the specified family, not costs imposed on
every individual frame.

Write $`T^\star_{F,\mathbb R}(n,a,b,\eta)`$ for the worst-case minimum
T-count over real frames, and $`T^\star_{F,\mathbb C,\mathrm{mag}}`$ for
the corresponding complex magnitude family. A displayed $`T^\star`$ bound
applies to each family explicitly named in its row. Supplied angles and phases admit certified evaluation;
classical coefficient generation and table preprocessing are excluded from
the quantum gate counts.

| Result | Statement selected for the paper | Required scope |
|---|---|---|
| **A. Exact complete-frame compilation** | Total elementary size $`\Theta(N)`$, CNOT count $`\Theta(N)`$ for $`n\ge2`$, and depth $`\Theta(n+N/(n+m))`$ | Every integer $`m\ge0`$; real and phase-dressed complex magnitude frames; workspace returned exactly; CNOT count is zero for $`n=1`$ |
| **B. Matching fault-tolerant frontier** | $`T^\star=\Theta(\sqrt{NL}+L+NL/q)`$, with $`O(NL)`$ Clifford cost for the upper construction | Real frame and, by Corollary 7, the phase-dressed complex magnitude frame; sufficient clean reservation $`a\ge C(n+h)`$ for a sufficiently large fixed $`C`$; literal-phase, complete-input error |
| **C. One-clean compilation** | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ | Real Hopf frame; $`n\ge1`$, $`L\ge6`$, $`a=1`$, $`b\ge L+n+7`$; complete-input error including workspace return |

Result A includes a separate CNOT lower bound with arbitrary one-qubit gates
free: a circuit with $`K`$ CNOTs has at most $`n+4K`$ relevant one-qubit
slots after removing inactive clean ancillas. No optimal leading constant is
claimed. Result B is precision-uniform under its clean reservation;
Result C uses a constant clean allocation. Here
$`\ell_*(n)=1+\log_2^*(n+2)`$, where $`\log_2^*`$ counts repeated
base-two logarithms until the value is at most one. At $`L=N`$ and
$`b=N+n+7`$, it gives $`T=O(N\ell_*(n))`$ and $`G=O(N^2)`$.
The endpoint's $`\Omega(N)`$ lower bound is unchanged.

Proofs: [A](../docs/COMPILER_THEOREM.md),
[B](../docs/FAULT_TOLERANT_COMPILER.md),
[C](../docs/ONE_CLEAN_COMPILER.md), with its
[grouped extension](../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit). The
[operator-source baseline](../docs/OPERATOR_SOURCE_COMPILER.md) gives
$`O(N+nL)`$ with two clean qubits and supplies the dirty-lookup
mechanism, exact source analysis, and earlier diagonal and multiplexor
constructions. The one-clean proof establishes the corollaries below
with their stated reservations.

## A compiler capability beyond Hopf frames

The one-clean mechanism also compiles **every literal diagonal and every
complete one-target U(2) multiplexor**. These corollaries establish the
mechanism's reach beyond the motivating frame family.

For a multiplexor with $`k\ge1`$ address qubits, one target, and $`M=2^k`$
arbitrary U(2) blocks, the [proof](../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
gives

```math
a=1,\quad b\ge2(L+k+9),\qquad
T^\star_{\mathrm{mux}}
=\Theta\!\left(\sqrt{ML}+L+\frac{ML}{k+2+b}\right),
\qquad G=O(ML).
```

The unbanked construction gives $`O(M+L)`$ T gates already at
$`b\ge L+k+9`$. Four addressed Euler factors reuse the same initialized
flag and arbitrary dirty pool. Certified approximate Euler coordinates
avoid singular inversion assumptions; their classical computation is
separate and has no universal efficiency guarantee for arbitrary input
evaluators. Literal block phases and complete-input error are preserved.
The upper constructions also cover $`k=0`$; the displayed matching claim
uses $`k\ge1`$.

The earlier two-clean constructions retain smaller base dirty reservations:
$`L+k+5`$ for a literal diagonal and $`L+k+7`$ for a multiplexor,
with matching banked bounds at twice those reservations. Reducing the clean
allocation does not supersede those distinct resource points.

This is one full addressed operation, not a claim for arbitrary many-qubit
unitaries. Direct composition of the noncommuting Hopf layers incurs the
baseline $`nL`$ term; conditional-suffix grouping reduces the number of
precision charges uniformly over real-frame accuracy. The
[current literature comparison](../docs/RELATED_WORK.md#12-contemporary-comparisons-and-the-broader-compiler-contribution)
distinguishes this clean-workspace guarantee from prior complete multiplexor
and arbitrary-unitary synthesis.

## Corollaries that stay in the paper

These complete the resource picture without becoming separate storylines.

| Corollary | Statement and placement |
|---|---|
| Dirty-bank refinement | With $`a=1`$, $`b\ge2(L+n+7)`$, the grouped compiler gives $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$. It matches the existing lower bound if $`\ell_*(n)^2L\le N`$ or $`b\le N/\ell_*(n)`$, subject to that allocation. State beside C; give bank scheduling in the appendix. |
| Zero-clean real-frame baseline | The [signal-symmetry corollary](../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) gives $`T=O(N+nL)`$, $`G=O(NL)`$, with $`a=0`$ and $`b\ge L+n+7`$. Its layerwise circuit borrows the signal qubit and includes its return error. This does not extend the grouped or literal-phase constructions to zero clean qubits. |
| Constructive T-depth bound | Choose between the layerwise and grouped [depth-oriented schedules](../docs/T_DEPTH_COMPILER.md) to obtain $`D_T=O(NL/b+\min\{nL+n^3,L\ell_*(n)+n^4\})`$ with two clean qubits, $`b\ge2(L+n+7)`$, and $`T,G=O(NL)`$. Retain as an appendix scheduling corollary for real frames; this trades gate count for depth and does not establish optimal T-depth or total circuit depth. |
| Simultaneous count and depth | With $`b\ge C(L+n+7+\sqrt{NL})`$ for a sufficiently large fixed C, [parallel dirty lookup](../docs/PARALLEL_DIRTY_LOOKUP.md) gives $`T=O(\sqrt{NL}+L\ell_*(n))`$ and $`D_T=O(\min\{nL+n^3,L\ell_*(n)+n^4\})`$ in the same two-clean real-frame circuit, with $`G=O(NL)`$. At fixed accuracy this is optimal-order T-count and polynomial-logarithmic T-depth with sufficiently large square-root dirty workspace; no depth optimality or selected-endpoint closure. |
| Literal diagonal synthesis | For $`\ell\ge6`$, one clean qubit gives $`O(N+\ell)`$ T gates at error $`2^{-\ell}`$ with $`b\ge\ell+n+7`$. For $`b\ge2(\ell+n+7)`$, the bound $`O(\sqrt{N\ell}+\ell+N\ell/b)`$ matches the diagonal lower bound. Use as a supporting compiler corollary. |
| Complex magnitude frame | Compose the one-clean grouped real frame and literal diagonal, with separately supplied certified phases and error $`\eta/2`$ per factor. This gives $`O(N+L\ell_*(n))`$ T gates at $`b\ge L+n+8`$, or $`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ at $`b\ge2(L+n+8)`$, with $`G=O(NL)`$. This compiles $`D_\phi W_{\mathbb R}`$; it is not arbitrary complex-unitary synthesis. |
| Smaller clean/dirty allocations | For every $`a,b\ge0`$, the real-frame borrowed-workspace compiler gives $`T=O(NL/q+L\sqrt N)`$ and $`G=O(NL)`$. Its matching splice requires $`h+b\le c\sqrt N`$ for fixed $`c>0`$. Retain as an appendix comparison; it can beat the grouped bound at low precision. |
| Fixed-parameter QBP robustness | Exact substitution preserves the global record. Magnitude bias is at most $`4\lvert a_j\rvert(\eta+\eta_O)`$; the separate phase-vector bias is at most $`4(\eta+\eta_O)`$. The full complex gradient, reflection-sum observables, rounded classical weights, and correlated dirty-bank reuse have explicit error and cost budgets. |

The exact and sufficient-clean complex extensions have their own proofs.
The [one-clean corollary](../docs/ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames)
uses full-isometry composition, including intermediate leakage and dirty-core
disturbance. The earlier two-clean complex baseline retains its proved
$`O(N+nL)`$ count at $`b\ge L+n+8`$. Leaf-phase derivatives use a
separate measurement stream; they are not extra columns of the magnitude frame.

## The common error contract

Let $`J_a`$ append the initialized clean work. The finite-precision target is

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

This includes all system inputs, arbitrary dirty inputs and their external
references, and leakage outside the clean-work subspace. Literal phases are
retained. No measurements, resets, postselection, or free supplied precision
sources are used.

The sufficient-clean compiler returns its dirty work exactly. In the
baseline two-clean construction, lookup banks, selectors, and suffix-control work
return exactly, while the operator core returns approximately within the
displayed norm. The one-clean construction also includes core disturbance in that norm.
In the conditional-suffix refinement, final return of the
logical suffix scratch and its clean predicate is also approximate within
that norm; intermediate table and selector subroutines return their work
exactly. The paper must preserve these distinctions wherever it states
workspace return.

QBP concerns simultaneous absolute accuracy of the raw coordinate gradient
at a fixed ideal parameter tuple. It assumes phase-calibrated controlled
access to the specified Hermitian-unitary observable. The matched-program
comparison separately charges quantum executions, per-execution circuit
cost, and classical output. It neither differentiates a discretely synthesized
gate word nor proves optimality among all gradient algorithms.

## The unresolved endpoint

The selected high-precision point is

```math
a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
\qquad
\Omega(N)\le T^\star\le O(N\ell_*(n)).
```

The one-clean construction attains the same upper bound with $`a=1`$
and the same dirty allocation; the selected two-clean endpoint above
remains open as well.

This gap is stated once in the main results and revisited in the discussion.
The paper does not claim that the grouped precision charge is necessary, or
that every constant clean allocation and every linear dirty allocation has
the same upper bound. Restrictions on particular source-processing interfaces
do not supply an additive full-frame lower bound.

The tree-generator and weighted-transport studies refine this open question.
Their [two-flag residual assembly](../docs/RESIDUAL_ASSEMBLY.md) is an explicit
complete-frame route at $`O(N+nL)`$ T cost, including native controls and
actual inverses. It supports the endpoint research without changing Results
A–C or adding a fourth principal result. The remaining joint precision task
is specified in the [research revision](../docs/OPEN_PROBLEM.md#revision-decision-and-next-bounded-pass);
its proposed linear gate budget is not an established publication claim.

## Main text and appendices

| Location | Contents |
|---|---|
| Main text | Operational necessity of the full frame; shared contract; Results A–C and resource regimes; one-clean diagonal and multiplexor capability; proof mechanisms; separately scoped complex-gradient consequence; open endpoint |
| Technical appendices | Echo, decoder, and router schedules; source preparation and exact primitive costs; residual composition and inherited failure compression; amplification and error sums; dirty queries/banks; Euler and diagonal proofs; lower-bound reductions; QBP concentration and classical preprocessing |

Exploratory attenuation, modular, batching, clock, and catalysis investigations
are outside the selected manuscript and the active repository. Their earlier
versions remain recoverable through the snapshot identified in
[provenance](../provenance/README.md). The current proof chapters contain all
dependencies of the three principal results and retained corollaries.

## Contribution and evidence boundaries

The local claims concern the complete-operator factorizations, workspace
schedules, precision-uniform residual composition, one-clean geometric
operator construction, and conditional-suffix grouping that yield these
resource guarantees. The earlier two-clean source and scheduling results
retain their separately established scope. Hopf geometry
and QBP records are inherited from the earlier Hopf work. UCG synthesis,
dirty lookup, Clifford-algebra loaders, geometric weighting, oblivious
amplification, and block-product compression retain their established attribution in the
[source map](../docs/SOURCE_MAP.md).

The submission is a theoretical construction paper with executable finite
checks. Those checks support fragile identities, conventions, and resource
ledgers; they do not prove universal statements, certify priority, or supply
a general elementary Clifford+T emitter. Hardware connectivity, physical
noise thresholds, optimal T-depth, optimizer convergence, and practical
end-to-end speedups are outside the selected claims.

## Scientific ingredients before final writing

Each ingredient has one primary home; the final manuscript draws on this
package rather than introducing unsupported research claims during writing.

| Ingredient | Scientific content and primary home |
|---|---|
| Necessity of the target | [Frame-safe contract](../docs/FRAME_SAFE_COMPILATION.md): universal fixed-decoder means force the full frame at regular points, up to common phase; sharp sensitivity and singular exceptions |
| Exact resources | [Exact theorem](../docs/COMPILER_THEOREM.md): complete constructions, all clean budgets, total size, CNOT-only count, and depth lower bounds |
| Precision sharing | [Fault-tolerant proof](../docs/FAULT_TOLERANT_COMPILER.md): sufficient-clean frontier, residual-dictionary sufficient condition, all charged source/history work, and complex extension |
| Constant-clean capability | [One-clean proof](../docs/ONE_CLEAN_COMPILER.md): real-frame primitive, one-clean diagonal and multiplexor frontiers, and phase-dressed complex-magnitude composition; [conditional-suffix proof, Section 10](../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit): grouped and banked one-clean real-frame bounds; [operator-source proof](../docs/OPERATOR_SOURCE_COMPILER.md): two-clean baseline, smaller dirty reservations, exact source T minima, dirty/reference return, and finite certified preprocessing |
| Operational consequence | [Approximation proof](../docs/QBP_APPROXIMATION.md): complete complex gradients, general reflection sums, weight errors, conditional-mean concentration under dirty-bank reuse, and quantum/classical costs |
| Comparisons and attribution | [Related work](../docs/RELATED_WORK.md) and [source map](../docs/SOURCE_MAP.md): input families, precision, initialized/borrowed work, current general-unitary baselines, and inherited compression |
| Reproducible evidence | [Verification](../docs/VERIFICATION.md): native finite one-clean primitives and two-clean circuits, actual inverse and rejected-work composition, source witnesses, gradient fixtures, and exact receipts; no full grouped elementary emitter |
| Honest unresolved question | [Open endpoint](../docs/OPEN_PROBLEM.md): the remaining full-frame gap between $`\Omega(N)`$ and $`O(N\ell_*(n))`$, the limit of the exact-source lower bound, and the unproved linear cost of a global block whose operator interface is already realized at $`O(N+nL)`$ |

Final writing assembles these established statements and proofs into one
argument with consistent notation and bibliography. No general elementary
Clifford+T emitter, hardware experiment, or solution of the open endpoint is
asserted by this package. External technical feedback on the complete draft
and final submission preparation follow.
