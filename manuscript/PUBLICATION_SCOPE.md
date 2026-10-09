# Publication scope

[Paper architecture](README.md) · [Read the argument](../REVIEW.md) · [Verification](../docs/VERIFICATION.md)

The selected publication is **one full-length theoretical compiler paper**
on prescribed Hopf differential frames. This is its canonical claim ledger:
Results A–D, retained corollaries, and their resource and error contracts.
Neither open resource question below is a prerequisite for these results.

Working title: **Exact and Fault-Tolerant Compilation of Hopf Differential Frames**.

## The central claim

Preparing a state fixes one column of a unitary. The fixed inverse-frame QBP
decoder uses a prescribed completion containing the state and its coordinate
directions. The paper compiles that stronger operator interface with optimal
exact size/depth and explicit fault-tolerant precision/workspace tradeoffs.
The exact and fault-tolerant theorems price the same prescribed operator under
different gate models; they do not assert that one circuit jointly minimizes
CNOT count, T-count, depth, and workspace.

## Four principal results

Write $`N=2^n`$, $`n\ge1`$. The exact model uses arbitrary one-qubit gates
and CNOTs with all-to-all connectivity. Its workspace budget is $`m`$ clean
qubits. The coherent Clifford+T model uses $`a`$ clean and $`b`$ dirty
qubits, $`q=n+a+b`$, and

```math
0\lt \eta\le1/64,\qquad
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad
h=1+\lceil\log_2(L+n+2)\rceil.
```

Here G counts elementary Clifford gates. Upper bounds hold for every prescribed
parameter tuple. Matching lower bounds are worst-case statements over the
named family, not costs imposed on every individual frame.

Write $`T^\star_{F,\mathbb R}(n,a,b,\eta)`$ and
$`T^\star_{F,\mathbb C,\mathrm{mag}}`$ for worst-case minimum T-count in
the real and complex magnitude families. A displayed $`T^\star`$ applies to
each explicitly named family. Supplied angles/phases admit certified evaluation;
classical coefficient generation and preprocessing are excluded from gate counts.

| Result | Statement selected for the paper | Required scope |
|---|---|---|
| **A. Exact complete-frame compilation** | Total elementary size $`\Theta(N)`$, CNOT count $`\Theta(N)`$ for $`n\ge2`$, and depth $`\Theta(n+N/(n+m))`$ | Every integer $`m\ge0`$; real and phase-dressed complex magnitude frames; workspace returned exactly; CNOT count is zero for $`n=1`$ |
| **B. Matching fault-tolerant frontier** | $`T^\star=\Theta(\sqrt{NL}+L+NL/q)`$, with $`O(NL)`$ Clifford cost for the upper construction | Real frame and, by Corollary 7, the phase-dressed complex magnitude frame; sufficient clean reservation $`a\ge C(n+h)`$ for a sufficiently large fixed $`C`$; literal-phase, complete-input error |
| **C. One-clean compilation** | $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ | Real Hopf frame; $`n\ge1`$, $`L\ge6`$, $`a=1`$, $`b\ge L+n+7`$; complete-input error including workspace return |
| **D. Simultaneous T-count and T-depth** | Uniform same-circuit bounds and matching precision/workspace regimes stated below | Complete real frames; $`a=2`$, $`b\ge17(L+n+7)`$; absolute constants and the full initialized-isometry error contract |

Result A includes a separate CNOT lower bound with arbitrary one-qubit gates
free: a circuit with $`K`$ CNOTs has at most $`n+4K`$ relevant one-qubit
slots after removing inactive clean ancillas. No optimal leading constant is
claimed. Result B is precision-uniform under its clean reservation;
Result C uses a constant clean allocation. Here
$`\ell_*(n)=1+\log_2^*(n+2)`$, where $`\log_2^*`$ counts repeated
base-two logarithms until the value is at most one.

### Result D: the current count–depth tradeoff

For two clean flags and $`b\ge17(L+n+7)`$, one complete real-frame
circuit satisfies, with absolute constants,

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+nL\right),\qquad G=O(NL).
```

In the explicit slowly growing precision range
$`6\le L\le\log_2(n+2)/16`$, the same theorem gives

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+n\right),\qquad G=O(NL).
```

Both resources match their worst-case lower bounds on these intervals when nonempty:

| Precision regime | Simultaneous matching interval |
|---|---|
| Every $`L\ge6`$ | $`17(L+n+7)\le b\le\sqrt{N/n}`$ |
| $`6\le L\le\log_2(n+2)/16`$ | $`17(L+n+7)\le b\le\sqrt{NL/n}`$ |

On either interval, $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$, attained simultaneously. Outside these
intervals, retain the stated upper bounds and their actual precision terms.
The sufficient condition $`L\le N/n^2`$ absorbs the general count's
$`nL`$ term into $`\sqrt{NL}`$, giving optimal-order count at every
eligible width; it does not establish depth optimality there.

D obeys the common full-input error contract below, with source preparation,
queries, inverses, and conditional work charged. The complex extensions of A–C
do not automatically extend D. T-depth permits arbitrary Clifford interlayers
and differs from total circuit depth.

## A compiler capability beyond Hopf frames

The one-clean mechanism also compiles **every literal diagonal and every
complete one-target U(2) multiplexor**.

For a multiplexor with $`k\ge1`$ address qubits, one target, and $`M=2^k`$
arbitrary U(2) blocks, the [proof](../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
gives

```math
a=1,\quad b\ge2(L+k+9),\qquad
T^\star_{\mathrm{mux}}
=\Theta\!\left(\sqrt{ML}+L+\frac{ML}{k+2+b}\right),
\qquad G=O(ML).
```

The unbanked construction gives $`O(M+L)`$ T gates at $`b\ge L+k+9`$.
Four addressed Euler factors reuse the same flag and dirty pool. Certified
approximate coordinates avoid singular inversion assumptions and preserve
literal phases and complete-input error. Their separate classical computation
has no universal efficiency guarantee for arbitrary evaluators. Upper bounds
also cover $`k=0`$; the matching claim uses $`k\ge1`$.

With two clean flags, the
[packed literal-diagonal compiler](../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames)
has base dirty reservation $`B_{\rm diag}=k+1+\lceil(L+4)/2\rceil`$,
T-count $`O(M+L)`$, and matching banked count at $`b\ge2B_{\rm diag}`$.
The two-clean multiplexor retains its separate base reservation $`L+k+7`$
and matching banked bound at twice that reservation.

## Corollaries that stay in the paper

| Corollary | Statement and placement |
|---|---|
| Dirty-bank refinement | With $`a=1`$, $`b\ge2(L+n+7)`$, the grouped compiler gives $`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, $`G=O(NL)`$. It matches the existing lower bound if $`\ell_*(n)^2L\le N`$ or $`b\le N/\ell_*(n)`$, subject to that allocation. State beside C; give bank scheduling in the appendix. |
| Zero-clean real-frame baseline | The [signal-symmetry corollary](../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) gives $`T=O(N+nL)`$, $`G=O(NL)`$, with $`a=0`$ and $`b\ge L+n+7`$. Its layerwise circuit borrows the signal qubit and includes its return error. This does not extend the grouped or literal-phase constructions to zero clean qubits. |
| Literal diagonal synthesis | For $`\ell\ge6`$, error $`2^{-\ell}`$ and $`T=O(N+\ell)`$ use base dirty reservation $`B=\ell+n+7`$ with one clean flag, or $`B=n+1+\lceil(\ell+4)/2\rceil`$ with two. At $`b\ge2B`$, $`O(\sqrt{N\ell}+\ell+N\ell/b)`$ matches the diagonal lower bound; $`G=O(N\ell)`$. |
| Complex magnitude frame | Compose the one-clean grouped real frame and literal diagonal, with separately supplied certified phases and error $`\eta/2`$ per factor. At $`a=1`$, this gives $`O(N+L\ell_*(n))`$ T gates at $`b\ge L+n+8`$, or $`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ at $`b\ge2(L+n+8)`$, with $`G=O(NL)`$. This compiles $`D_\phi W_{\mathbb R}`$; it is not arbitrary complex-unitary synthesis. |
| Smaller clean/dirty allocations | For every $`a,b\ge0`$, the real-frame borrowed-workspace compiler gives $`T=O(NL/q+L\sqrt N)`$ and $`G=O(NL)`$. Its matching splice requires $`h+b\le c\sqrt N`$ for fixed $`c>0`$. Retain as an appendix comparison; it can beat the grouped bound at low precision. |
| Fixed-parameter QBP robustness | Exact substitution preserves the global record. Magnitude bias is at most $`4\lvert a_j\rvert(\eta+\eta_O)`$; the separate phase-vector bias is at most $`4(\eta+\eta_O)`$. The full complex gradient, reflection-sum observables, rounded classical weights, and correlated dirty-bank reuse have explicit error and cost budgets. |

Earlier [depth-oriented schedules](../docs/T_DEPTH_COMPILER.md) and
[parallel dirty lookup](../docs/PARALLEL_DIRTY_LOOKUP.md) remain proof
dependencies and alternative resource points, including smaller literal
dirty reservations or high-precision regimes where they are sharper.
Present their applicable fallbacks in the appendices; Result D states the
current uniform tradeoff and does not dominate every older construction.

The [operator-source baseline](../docs/OPERATOR_SOURCE_COMPILER.md) retains
its source/digit primitives, dirty queries, amplification, and certified table
generation. Its two-clean real-frame bound is $`O(N+nL)`$ at
$`b\ge L+n+7`$; its complex baseline has the same count at
$`b\ge L+n+8`$. Leaf-phase derivatives use a separate measurement stream;
they are not additional columns of the magnitude frame.

## The common error contract

Let $`J_a`$ append the initialized clean work. The finite-precision target is

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

This includes all system inputs, arbitrary dirty inputs and their external
references, and leakage outside the clean-work subspace. Literal phases are
retained. No measurements, resets, postselection, QRAM, free target-frame
oracle, or supplied precision source is used. Coherent queries and actual
inverse calls are charged.

The sufficient-clean compiler returns its dirty work exactly. In the
baseline two-clean construction, lookup banks, selectors, and suffix-control work
return exactly, while the operator core returns approximately within the
displayed norm. The one-clean construction also includes core disturbance in that norm.
In the conditional-suffix refinement, final return of the
logical suffix scratch and its clean predicate is also approximate within
that norm; intermediate table and selector subroutines return their work
exactly. The paper must preserve these distinctions wherever it states
workspace return.

The [frame-safe proof](../docs/FRAME_SAFE_COMPILATION.md) establishes the
need for the prescribed completion for the fixed inverse-frame decoder at
regular points, up to common phase, with stated singular exceptions.
QBP concerns simultaneous absolute accuracy of the raw coordinate gradient
at a fixed ideal parameter tuple. It assumes phase-calibrated controlled
access to the specified Hermitian-unitary observable. The matched-program
comparison separately charges quantum executions, per-execution circuit
cost, and classical output. It neither differentiates a discretely synthesized
gate word nor proves optimality among all gradient algorithms or an
end-to-end application advantage.

## Two unresolved resource questions

The selected constant-clean high-precision open problem is

```math
a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
\qquad
\Omega(N)\le T^\star\le O(N\ell_*(n)).
```

One clean qubit attains the same upper bound and $`G=O(N^2)`$ at that dirty
allocation. The selected two-clean endpoint remains open as well.

The grouped precision charge is not claimed necessary, nor is the same upper
bound claimed for every constant clean count or linear dirty prefactor.
Source-interface restrictions do not supply additive full-frame lower bounds.

At fixed accuracy and sufficient $`b=\Theta_\eta(\sqrt N)`$, the
complete real-frame schedule attains $`T=\Theta_\eta(\sqrt N)`$ with
$`D_T=O_\eta(n)`$. The unrestricted depth lower bound remains
$`\Omega(1)`$. This large-width gap is separate from the high-precision
count endpoint. Result D's literal width threshold excludes the displayed
endpoint allocation, and its unary-source mechanism does not scale to that
precision.

The [research archive](../research/README.md) records tried routes, proved
restricted cases, component improvements, and their scoped limits. They are
not additional principal results or unrestricted lower bounds. The completed
[state-based QBP supplement](../supplements/state_based_qbp/README.md) uses a
changed decoder and remains a separate result; the fixed-decoder necessity
claim here does not apply to every gradient algorithm.

## Contribution and evidence boundaries

The local claims concern complete-operator factorizations and charged workspace,
residual, one-clean, conditional-suffix, and query/source schedules. Earlier
two-clean results retain their stated scope. Hopf geometry and QBP records are
inherited from the earlier Hopf work. UCG synthesis, dirty lookup, Clifford-algebra
loaders, geometric weighting, oblivious amplification, and block-product
compression retain their attribution in the [source map](../docs/SOURCE_MAP.md).

The submission is a theoretical construction paper with executable finite
checks. Those checks support fragile identities, conventions, and resource
ledgers; they do not prove universal statements, certify priority, or supply
a general elementary Clifford+T emitter. Hardware connectivity, physical
noise thresholds, the unrestricted T-depth frontier, optimizer convergence, and practical
end-to-end speedups are outside the selected claims.

## Proof chain and manuscript placement

| Ingredient | Primary proof home and role |
|---|---|
| Target and necessity | [Hopf interface](../docs/HOPF_INTERFACE.md) and [frame-safe contract](../docs/FRAME_SAFE_COMPILATION.md): prescribed columns, chart continuation, fixed-decoder necessity |
| A | [Exact theorem](../docs/COMPILER_THEOREM.md): all clean budgets, size, CNOT-only count, and depth |
| B | [Fault-tolerant theorem](../docs/FAULT_TOLERANT_COMPILER.md): shared precision source, residual dictionaries, failure history, full error, complex extension, lower bounds |
| C and compiler corollaries | [One-clean proof](../docs/ONE_CLEAN_COMPILER.md), [conditional-suffix Sections 1–10](../docs/CONDITIONAL_SUFFIX_COMPILER.md), and [operator-source primitives](../docs/OPERATOR_SOURCE_COMPILER.md); [borrowed-workspace proof](../docs/BORROWED_WORKSPACE_COMPILER.md) supplies queries, sector echoes, and the arbitrary-budget comparison |
| D | [Uniform-precision theorem](../docs/UNIFORM_PRECISION_DEPTH.md) and its linked source/query schedules: same-circuit upper bounds and matching intervals |
| QBP consequence | [Exact record and costs](../docs/QBP_CONSEQUENCE.md), [approximation budgets](../docs/QBP_APPROXIMATION.md) |
| Attribution and evidence | [Related work](../docs/RELATED_WORK.md), [source map](../docs/SOURCE_MAP.md), [verification](../docs/VERIFICATION.md), and [internal audit](../docs/CORE_CLAIM_AUDIT.md) |

The main text follows the target, shared contract, Results A–D, compiler
corollaries, operational consequence, and two open questions. Technical
appendices contain exact echo/decoder/router schedules, source preparation,
dirty queries, residual composition, amplification, error and resource sums,
Euler/diagonal details, lower-bound reductions, and classical preprocessing.

The stopping criteria are complete proofs, charged resource/error contracts,
consistent attribution, scoped finite checks, and passing verification gates.
The completed internal audit is not external peer review. Practical crossover
constants and a scalable native emitter remain outside the selected claims.

Closing either gap or adding further fixtures is not a readiness requirement.
