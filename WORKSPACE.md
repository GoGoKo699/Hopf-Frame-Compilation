# Continuing research workspace

This is the current research checkpoint, revised **3 October 2026** from
verified main `88350ab15ef22bdc49c9bfec3fcd2764283edbf7` (PR #85).
That baseline passes 495 tests, four exact fault-tolerant receipt suites,
and all five CI checks. The charged query-correction audit below gives an
exact fused word, but its generic correction retains a fresh lookup.
Two early linear-depth allowances remain, and the complete-frame frontier
is unchanged.

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
consolidates three completed refinements:

- The [nonuniform dirty indicator](docs/NONUNIFORM_DIRTY_INDICATOR.md)
  has linear count/work and logarithmic address T-depth, with literal
  dirty/reference return. At sufficient width this makes the low-precision
  late tail sublinear.
- The [protected source](docs/PROTECTED_UNARY_SOURCE.md) serves every early
  group. One global $`2\delta`$ initialized-isometry bound includes final
  source/flag leakage. Source preparation and its actual inverse contribute
  $`O(L+\log(n+2))`$ total depth with $`\delta=\eta/8`$.
- The [windowed predicates](docs/WINDOWED_GROUP_PREDICATES.md) use
  $`J=\lceil\log_2(n+2)\rceil`$ groups per full window and
  $`2J+1`$ additional conditional-zero bank bits. The two original clean
  flags and two returned dirty predicate helpers still suffice.

The windowed predicates' aggregate activity depth in the uniform
low-precision regime is

```math
O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}+\log(n+2)\right)=o(n).
```

At fixed accuracy and sufficient square-root dirty width, the current
same-circuit depth accounting is:

| Contribution | Current total allowance |
|---|---|
| Logical shifts and incremental selectors | O(n) |
| Program prefetch and unload | O(n) |
| Activity predicates across group windows | O(n log log n/log n+log n) |
| Protected-source and initial-bank-predicate boundaries | O(log n) |
| Late queries, sources, and predicates | o(n) |

These are construction costs, not lower bounds. The two remaining linear
rows keep the complete-frame upper bound at O(n); the count, Clifford
count, matching intervals, and literal width threshold retain the
established theorem and fallbacks.

The window caches return exactly on active inputs with conditional-zero
cache and an arbitrary source core. On $`Hh=00`$, the full window is
identity for arbitrary cache, source, logical, dirty, and reference inputs.
Each block-zero cache is erased before its targets change. Each chain bit
is consumed while its later controls are still valid. Future logical bits
may be temporary group work, but they must return before the cache is read
again. No clean-cache assumption is made on the inactive sector, and no
identity is asserted for arbitrary initial h equal to one.

The [shared-prefix query audit](docs/SHARED_PREFIX_QUERY_AUDIT.md)
finds that a retained dirty indicator leaves an unknown output mask.
A correction computed before the local address changes is generally
stale, even for commuting XOR queries. A one-clean-bit native example
shows precisely which baseline promise repairs the small construction.

More generally, consider an encoder that preserves the original prefix x,
followed by exact reusable row readers that cannot access x at all.
With c initialized cache bits and any number of arbitrary dirty bits,
at most $`2^c`$ distinct row signatures fit. This counts every initialized
wire used by the encoder, including conditional logical zeros. For arbitrary
legal Hopf tables over G target bits, the requirement is

```math
c\ge\min\{d,2^G-1\},\qquad d=n-k.
```

Even granting all k remaining logical bits and both flags to this cache,
when $`d\gt k+2`$ its window satisfies $`G\le\log_2(k+3)`$.
It therefore cannot combine a growing number of the current early groups
in that latter portion, without assuming a one-hot encoding.
It does not constrain a reader that accesses x, even by CNOT,
or a circuit that implements a fused window without exposing row queries.

The [charged correction](docs/SHARED_PREFIX_QUERY_AUDIT.md#6-a-charged-two-group-program-refresh)
now gives an exact two-group query/body/inverse-query identity. It retains
one prefix indicator and changes the common program by the XOR of the
old and new table words, reading the second address after the first body.
A reused activity flag has an explicitly charged transition; no extra
old-flag copy is assumed. All actual inverses and arbitrary dirty/source
inputs remain in the contract.

This is not a general depth saving. Fixing the first program to a legal
constant leaves the next program arbitrary, so a fixed Clifford XOR turns
the correction into a fresh query. A rank-one difference is still hard in
this reduction sense if its prefix factor is arbitrary. An affine parity
factor does give an exact two-Toffoli correction with one returned dirty
helper, but independently prescribed Hopf rows need not have that form.
This reduction is not a depth lower bound. Stop the generic difference-only
route unless a new generic lookup schedule or uniformly cheap factorization
is supplied; do not pursue more stale-mask or dirty-cache variants.

The next selected task is the logical-stage row: write two consecutive
unary stages as an exact matrix whose entries are Laurent polynomials in
the cyclic source shift S, with $`S^q=I`$. Keep the full source action,
changed logical controls, literal phase, and actual cleanup. Start at
$`q=8`$ and four logical modes, then test closure at eight modes before
claiming any recurrence. For this exact rewrite, equality only on the
prepared Fourier source is insufficient: verify every source character
or an exact polynomial identity.
The existing ideal four-mode factorization is a starting comparison, not
a new result. A fixed pairwise saving cannot reduce the total asymptotic
depth; a growing-block recurrence must include selectors, work return,
program access, and its T/Clifford/width costs. A different unitary body
that agrees only on the ideal source remains eligible with a proved global
initialized-isometry and work-return contract; failure of an exact rewrite
does not exclude it. The query row remains separately linear unless a
later construction improves it as well.

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
