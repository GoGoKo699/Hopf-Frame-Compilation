# Continuing research workspace

This is the current research checkpoint, revised **3 October 2026** from
verified main `efbb75fdf3573a9d7403e2d3afd207b5180db053` (PR #82).
That baseline passes 481 tests, four exact fault-tolerant receipt suites,
and all five CI checks. This continuation proves protected unary-source
reuse across all early groups and adds four bounded checks. The complete-frame frontier is unchanged;
the next windowed-predicate test is unproved.

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
| Constant-clean high-precision count | At $`a=2`$, $`L=N`$, $`b=N+n+7`$: $`\Omega(N)\le T^\star\le O(N\ell_*(n))`$, where $`\ell_*(n)=1+\log_2^*(n+2)`$ |

The uniform-depth width hypothesis excludes that literal endpoint.
Its unary-source mechanism also does not scale to endpoint precision.
The separate state-based and phase-dressed complex results retain their
own contracts; the new real-frame depth theorem is not automatically a
complex-frame theorem.

## Revision decision and next bounded task

The [current decision record](docs/OPEN_PROBLEM.md#revision-checkpoint-and-selected-next-test)
separates the completed indicator and source milestones from the next test.
The [nonuniform indicator](docs/NONUNIFORM_DIRTY_INDICATOR.md) now proves
an exact dirty indicator for $`S=2^s`$ outputs:

```math
|x,Y,W\rangle\longmapsto|x,Y\oplus e_x,W\rangle,
\qquad T,G,w=O(S),\qquad D_T=O(\log_2(s+2)).
```

It uses no initialized helpers and returns arbitrary dirty/reference
inputs with literal phase. The arbitrary-partition echo, rounded schedule,
linear total count, live width, and logarithmic depth are proved. Five
bounded tests cover its fragile interfaces and exact resource inequalities.
The rectangular query consequently has square-root count and logarithmic
address T-depth at sufficient square-root dirty width for a fixed word.

This removes the previous linear indicator allowance from the late tail.
In the uniform low-precision range, its depth is now
$`O(NL/B^2+M[L+\log(n+2)])`$, with the second term o(n), whenever
the new query's sufficient reservation fits. The full-frame literal width
threshold retains its established fallback. After the two refinements,
the early accounting is:

| Contribution | Current accounting |
|---|---|
| Early target-dependent shifts and selector updates | O(g) per group; total group height O(n) |
| Early source preparation, conversion, and actual inverse | One global pair of depth O(L+log n) |
| Early outer predicates and their inverses | O(log n) per group |
| Early program loading and unloading | O(log n) per group |

There are O(n/log n) early groups. These are separate upper allowances,
not lower bounds; the tail improvement alone does not make the full
frame sublinear.

The [protected-source theorem](docs/PROTECTED_UNARY_SOURCE.md) now removes
the repeated source-boundary cost. A fixed terminal logical bank supplies
the source, preparation auxiliaries, and convolution work. Its original
zero predicate is stored once in H; the other initialized flag h records
each group's outer-zero predicate excluding that bank, conjoined with H.
Actual preparation U is unconditional. The entire middle is identity on
H equal to zero with h initially zero, so U and its actual inverse cancel
on arbitrary inactive inputs. This does not assert identity for an
arbitrary initial h; both external flags belong to the initialized contract.

On the active sector every group's temporary work and h return exactly
for any source input. A single $`2\delta`$ bound covers the entire early
segment, including final bank leakage and H erasure. With
$`\delta=\eta/8`$, the two global source/predicate boundaries have
$`O(L+\log(n+2))`$ T-depth. The actual source need not return exactly
before H is erased. The proof compares the complete isometry, with no
reset or discarded rejected component.

Logical stages, outer predicates, and program queries remain separately
$`O(n)`$. The full-frame count, Clifford count, depth, matching intervals,
and literal width threshold retain the existing theorem and fallbacks.

The next selected **unproved component test** is a window of cached outer
predicates. Reserve $`2J+1`$ additional conditional-zero bits in the protected bank. For
$`J=\lceil\log_2(n+2)\rceil`$ consecutive early groups, cache each
group block's zero bit and a suffix-product chain, plus one predicate for
the suffix beyond the whole window, excluding the protected bank. Use
these bits and H to set the external group flag h.

Proceed in this order:

1. Prove a consume-before-change schedule: erase a block's cached zero
   bit before its logical targets change; erase each consumed chain bit
   while its controls are still valid. Return all cache bits at the end
   of the window, including arbitrary inactive inputs and references.
2. Keep cache storage disjoint from source, convolution, program, and
   selector work. Only the two original external flags may be initialized;
   all additional zero promises must follow from H.
3. Charge one long predicate pair per window and the short block/chain
   cleanup. Test the candidate depth $`O(n\log\log(n+2)/\log(n+2))`$
   in the uniform low-precision regime before claiming any improvement.

Stop if cleanup needs an unchanged logical input after that input has
become a target, assumes clean inactive cache bits, resets source leakage,
or overlaps live work. A successful predicate test would still leave
logical stages and program queries with linear allowances. Add only
bounded checks that resolve a new concrete interface risk after the
symbolic contract and resource ledger close.

Keep the high-precision endpoint parked until an explicit new global
native identity or encoding rule survives the
[endpoint selection conditions](docs/OPEN_PROBLEM.md#separate-complete-frame-question).
Completed source transport, canonical completion, and unchanged-angle
scattering attempts do not remove the joint interior's precision charge.
A lower-bound project remains eligible, but must specify an invariant
valid with arbitrary Clifford interlayers and returned dirty work.
The [focused literature check](docs/RELATED_WORK.md#22-scope-of-recent-depth-lower-bounds-3-october-2026)
does not supply such a growing bound in the selected large-width regime.

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
