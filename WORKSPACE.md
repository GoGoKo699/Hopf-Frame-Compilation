# Continuing research workspace

This is the current research checkpoint, revised **3 October 2026** by the
[bounded core-claim audit](docs/CORE_CLAIM_AUDIT.md). Its reviewed baseline is
main `67a4cf25be30ec0440dd6e63a85181ca9c38121d`, after merged PR #88.
That baseline passes 505 tests, four exact fault-tolerant receipt suites,
and all five branch/PR checks, including rendered mathematics.

**Decision:** the bounded audit found no unresolved claim-level blocker;
freeze the selected scientific package. Automatic construction passes stop,
and no new frontier construction is selected. The two open
resource gaps are research opportunities, not prerequisites of the selected
paper. The [scope](manuscript/PUBLICATION_SCOPE.md) now includes the strongest
proved count/depth theorem; final manuscript writing remains on hold.

## Mandate and model

Continue theorem-led research only in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized. Manuscript writing and
release work are outside this pass. Application-level advantage is outside
the selected scope; existing cost comparisons and classical baselines are
retained boundaries, not pending tasks.

Use small analytic examples and finite checks. Preserve literal phases,
actual inverses, and the prescribed complete-frame columns. The compiler
has no measurements, intermediate resets, QRAM, supplied catalysts,
uncharged initialized history, free target-frame oracle, or uncharged coherent evaluator.
All dirty work may be entangled with an external reference. Its return and
all clean leakage belong to the full initialized-isometry error. A logical
zero suffix supplies temporary clean work only on its active sector; the
inactive action must be proved separately on arbitrary inputs.
T-depth permits arbitrary Clifford interlayers; their elementary count
and nonzero physical depth remain separate costs.

## What is established

| Track | Analytic result and proof home | Implemented evidence and remaining boundary |
|---|---|---|
| Exact prescribed frames | Matching size and elementary depth across clean-work budgets; [exact theorem](docs/COMPILER_THEOREM.md) | Small complete matrices and reversible decoder/router ledgers. This exact-gate theorem is distinct from Clifford+T depth |
| Sufficient-clean frame T-count | Matching worst-case precision/workspace frontier for real and phase-dressed complex magnitude frames; [fault-tolerant proof](docs/FAULT_TOLERANT_COMPILER.md) | Analytic full-frame theorem under its sufficient clean reservation; it is distinct from the constant-clean constructions |
| Constant-clean frame T-count | One-clean grouped bound and separately banked refinements; [grouped proof](docs/CONDITIONAL_SUFFIX_COMPILER.md), [one-clean theorem](docs/ONE_CLEAN_COMPILER.md) | Native source/amplification and bounded group checks. The selected high-precision endpoint is still open |
| Complete real-frame T-depth | Uniform and low-precision same-circuit bounds below; [uniform theorem](docs/UNIFORM_PRECISION_DEPTH.md) | Unary source components, native bilinear leaves, symbolic complete queries, and exact allocation checks. No scalable native frame emitter or unrestricted depth-optimality theorem |
| State-based real/complex QBP | Fine state preparation, actual coarse reference/inverse, both original raw-gradient streams, and charged reconstruction; [consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) | Certified bounded rows/tables, preparations, and a fixed complex residual-to-gradient integration. The selected theorem and bounded integration are complete |

The [verification map](docs/VERIFICATION.md) distinguishes analytic
constructions, native components, exact certificates, and floating-point
checks. Tests support their stated finite contracts; they are not proofs
of asymptotic bounds or external peer review.

For QBP, retain accuracy K separately from the state precision floor
$`P=\max(n,K)`$. The initialized system, interference branch, compiler
flags, and observable work are separate reservations. The two-clean count
is not the total initialized width of a QBP execution. State preparation
does not substitute for a prescribed frame on arbitrary logical inputs.
The [cost comparison](docs/QBP_COST_COMPARISON.md) charges all executions,
observable access, preprocessing, and classical reconstruction. No general
end-to-end gradient advantage follows.

## Current complete-frame frontier

Put $`N=2^n`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
and $`B_0=L+n+7`$. For $`0\lt\eta\le1/64`$, two external clean
flags, and the literal dirty reservation $`b\ge17B_0`$, one real-frame
circuit has absolute-constant bounds

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

For the explicit slowly growing precision range, the same theorem gives

```math
6\le L\le\frac{\log_2(n+2)}{16},\qquad
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+n\right),\qquad G=O(NL).
```

| Regime, always above the stated width threshold | Simultaneous matching interval, when nonempty |
|---|---|
| Every precision | $`b\le\sqrt{N/n}`$: $`T^\star=\Theta(NL/b)`$, $`D_T^\star=\Theta(NL/b^2)`$ |
| Displayed slowly growing precision range | $`b\le\sqrt{NL/n}`$: the same matching orders |
| Fixed accuracy | [Blocked bilinear theorem](docs/BLOCKED_BILINEAR_LOOKUP.md): $`T=O_\eta(\sqrt N+N/b)`$, $`D_T=O_\eta(N/b^2+n)`$ |

The general count retains nL outside its matching interval. The sufficient
condition $`L\le N/n^2`$ absorbs this into $`\sqrt{NL}`$ and gives
optimal-order count at every eligible width. Constants in the fixed-eta
predecessor may depend on eta; those in the uniform theorem do not.
The small coefficient in the low-precision range is a sufficient
asymptotic condition, not a practical crossover estimate.

Two scientific gaps remain:

| Question | Known boundary |
|---|---|
| Fixed-accuracy, large-width T-depth | At sufficient $`b=\Theta_\eta(\sqrt N)`$, optimal-order $`T=\Theta_\eta(\sqrt N)`$ accompanies $`D_T=O_\eta(n)`$; the unrestricted depth lower bound is only $`\Omega(1)`$ |
| Constant-clean high-precision count | At $`a=2`$, $`L=N`$, $`b=N+n+7`$, $`n\ge3`$: $`\Omega(N)\le T^\star\le O(N\ell_*(n))`$, where $`\ell_*(n)=1+\log_2^*(n+2)`$ |

The uniform-depth width hypothesis excludes that literal endpoint.
Its unary-source mechanism also does not scale to endpoint precision.
The separate state-based and phase-dressed complex results retain their
own contracts; the new real-frame depth theorem is not automatically a
complex-frame theorem.

## Research decision and stopping rules

### Vision: one operator-compilation question, two resource models

The central question is whether preserving a prescribed completion costs
more than preparing its first column. We have an optimal exact
workspace–depth theorem, a sufficient-clean optimal T-count frontier,
constant-clean constructions, and simultaneous matching T-count/T-depth
in explicit real-frame regimes. QBP explains why the additional columns
matter. The completed state-based decoder is a separately scoped option;
it does not replace the prescribed-frame theorem.

The goal is a defensible account of these resource tradeoffs, with each
claim attached to its assumptions and proof. It is not a requirement to
solve every width and precision regime. None of these statements promises
acceptance, external validation, or an application-level speedup.

### What the recent refinements bought

| Result | Established gain | Present effect on the full-frame frontier |
|---|---|---|
| [Uniform precision and depth](docs/UNIFORM_PRECISION_DEPTH.md) | One circuit attains the displayed count/depth bounds and two explicit matching intervals | This is the latest complete-frame frontier advance and belongs in the principal theorem map |
| [Nonuniform indicators](docs/NONUNIFORM_DIRTY_INDICATOR.md), [protected source](docs/PROTECTED_UNARY_SOURCE.md), [windowed predicates](docs/WINDOWED_GROUP_PREDICATES.md) | Reduce late-query, source-boundary, and activity overheads below linear depth at fixed accuracy and sufficient width | Useful components; the two remaining linear allowances still dominate |
| [Query-cache and correction audit](docs/SHARED_PREFIX_QUERY_AUDIT.md) | Exact reader-capacity restriction, charged reuse identity, and the remaining generic correction query | Closes specific proposed shortcuts; supplies no unrestricted depth lower bound |
| [Retained-source completion](docs/UNARY_PHASE_GRADIENT.md#10-two-source-shifts-for-the-stabilizer-completion) | Two source shifts and logarithmic completion depth for even-height groups, with full-frame and all-source identity | The ordered transport remains linear in group height; the complete-frame order does not improve |

At fixed accuracy and sufficient square-root dirty width, the current
same-circuit accounting remains:

| Contribution | Aggregate T-depth allowance |
|---|---|
| Logical transport and incremental selectors | $`O(n)`$ |
| Program queries and actual unloads | $`O(n)`$ |
| Windowed activity predicates | $`O(n\log\log(n+2)/\log(n+2)+\log(n+2))`$ |
| Source and initial-bank-predicate boundaries | $`O(\log(n+2))`$ |
| Late queries, sources, and predicates | $`o(n)`$ |

These are construction costs, not lower bounds. The extracted transport
uses three shifts per target pair; including the shallow completion gives
$`3g/2+2`$ slots for even height g, versus g in the retained original
compiler. Its factorization is useful structure, not a faster full group.
Transport can move the pair's 00 state into its complement, so the
completion's invariant-sector cache argument does not transfer to it.
Noncommuting transported blocks also prevent direct scalar recursion.

### Completed gate: bounded claim-to-proof audit

The [audit record](docs/CORE_CLAIM_AUDIT.md) traces Results A–D to their
upper constructions, matching lower bounds, literal work reservations,
precision regimes, and global error contracts. It also checks the borrowed
source premises and the fixed-decoder QBP interface. All four results pass
this internal review. One conservative tail-work expression is now described
as an upper bound; no theorem or resource frontier changes.

The selected scientific scope is frozen. This means the claimed package has
no identified blocker from this bounded review, not that its proofs are
formally verified or independently peer reviewed. The analytic constructions,
finite certificates, native components, and missing scalable emitter remain
distinct in the evidence map. Practical crossover constants are not supplied.

Final writing begins only when requested. A general native emitter, practical
crossover study, experiment, or application advantage is not a default
prerequisite. Repair a discovered defect; otherwise apply the reopening gate
below instead of adding fresh requirements to the completed package.

### Reopening research requires a qualifying mechanism

The smallest substantial new large-width target is a complete real-frame
circuit, at fixed accuracy and sufficient $`b=\Theta_\eta(\sqrt N)`$,
with two external clean flags and

```math
D_T=o(n),\qquad T=O_\eta(\sqrt N),\qquad G=O_\eta(N).
```

This would improve the frontier without claiming optimal depth. Both
linear rows must become sublinear, or a fused circuit must bypass them.
Improving only transport or only queries is a component milestone; a pass
pursuing one must state the other remaining allowance at the outset.

Before another construction pass, require an explicit native identity or
encoding rule and symbolic depth/count/width/error recurrences. Price
program generation queries, actual unloads, conditional logical zeros,
source boundaries, and dirty/reference return. For a transport component,
a sublinear group-depth recurrence must fit the retained asymptotic group
workspace and gate budgets. A different body agreeing only on the ideal
source remains eligible with a global initialized-isometry proof; matching
only the first logical column is insufficient.

**Candidate stop:** a transport-improvement proposal fails its target if
its expanded schedule still pays a constant number of sequential shifts
per target pair. A query-improvement proposal fails if its aggregate query
cost remains unchanged after corrections and unloads are charged. Fresh
queries are allowed if a new schedule makes their aggregate cost sublinear;
a component milestone may retain the other unresolved linear allowance.
Unpriced initialized history or an unproved return invalidates either
proposal. More modes, determinant rearrangements, and constant-factor
savings do not meet this pass's asymptotic target. One explicit failure
record in the relevant proof home is sufficient. This rejects the candidate
and its stated ledger, not all possible compilers.

**Selection stop:** if no qualifying rule emerges from a bounded analytical
selection pass, record “no selected construction” and park the frontier.
That is the current status. Existing ingredients do not establish that
closure is close, nor that the missing idea lies in already available theory.

The high-precision endpoint remains separately parked until a global native
identity meets its [literal endpoint conditions](docs/OPEN_PROBLEM.md#separate-complete-frame-question).
An unrestricted lower-bound project would need a new invariant valid with
arbitrary Clifford interlayers, initialized flags, and returned dirty work.
The existing interface obstructions do not supply it. Neither route is
selected simply to keep exploration active.

**Project stop reached:** the bounded consolidation audit has no unresolved
claim-level blocker, so the current scientific package is complete for its
selected scope. Preserve both open gaps explicitly, keep all supporting
proofs accessible, and move to manuscript preparation only on request.
Reopen for a concrete qualifying idea, a discovered defect, or an explicit
change of scope—not merely another request to “continue.”

## Proof map and completed work

- Current depth: [uniform composition](docs/UNIFORM_PRECISION_DEPTH.md)
  uses [unary groups](docs/UNARY_PHASE_GRADIENT.md),
  [blocked bilinear queries](docs/BLOCKED_BILINEAR_LOOKUP.md), and the
  [older all-precision hybrid](docs/PARALLEL_DIRTY_LOOKUP.md).
- Unary groups use [ideal-angle stability](docs/HOPF_ERROR_ACCUMULATION.md),
  [incremental selectors](docs/GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables),
  and charged preparation/actual return. The late tail retains the
  [capped source certificate](docs/AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision).
- Scoped source, reflection, transport, and endpoint failures remain in
  the [route audit](docs/OPEN_PROBLEM.md#what-the-failure-diagnostics-actually-rule-out)
  and their linked proof homes. They are not unrestricted lower bounds.
- State-based QBP starts at its [consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md)
  and [coverage map](docs/VERIFICATION.md#state-based-qbp-coverage).
  Certified general input preprocessing, a variable-size native emitter,
  and general guarded decoding remain optional software; none is selected
  without a concrete input/output requirement. Do not repeat the completed
  residual-to-gradient fixture sequence.

## Restore and verify

Start with repository status and the relevant proof chapter. For substantive
circuit or theorem changes, the retained verification workflow is:

```bash
git status --short --branch
git rev-parse HEAD
python -m compileall -q compiler_robust_hopf scripts tests validate.py
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python scripts/complex_coarse_native_example.py
python scripts/residual_qbp_native_example.py --q 16
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
python scripts/unified_resource_ledger.py --n 12 --format json
python scripts/strict_zero_echo_ledger.py --n 12 --format json
```

Documentation-only revisions need focused link, math, and presentation
checks locally; retain all required remote merge gates. Those include
Python 3.11/3.13 validation, four exact-receipt suites, and rendered
presentation. Record new scientific claims in their proof home and the
research-status map so continuation does not depend on an old chat.
