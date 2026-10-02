# Continuing research workspace

This is the entry point when a previous conversation or execution workspace
is unavailable. Proofs and decisions live in the repository.

The 2026-10-02 source-precision pass starts from verified main
`018b22f8e95162a273f4da0510550a1373a00241`, after the matching
modest-width depth theorem. Check later commits before
continuing.
The selected state-based Hopf QBP construction and its bounded-input audit are complete; the
[consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) is their entry point.
It prepares a state and changes the gradient decoder while retaining all
raw coordinates, singular angles, and both complex-gradient streams.
The active scientific targets are the constant-clean complete-frame
endpoint and optimal T-depth. Application-level advantage is outside the
current research scope by the author's decision; retain the existing cost
comparisons and classical baselines as boundaries, not a pending task.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized. Keep manuscript writing
and release work outside this research pass. Use small analytic examples
and finite checks, without large simulations, QRAM, resets inside a
compiler execution, supplied catalysts, or hidden initialized work.

1. Begin active depth work with the [amortized tradeoff](docs/AMORTIZED_DIRTY_LOOKUP.md),
   then its [batched allocation](docs/BATCHED_DIRTY_LOOKUP.md),
   [depth proof](docs/T_DEPTH_COMPILER.md), and
   [source-depth certificate](docs/SOURCE_T_DEPTH.md), then the
   [routed indicator](docs/PARALLEL_DIRTY_LOOKUP.md). The completed
   [state-based QBP theorem](docs/STATE_BASED_QBP_THEOREM.md) supplies its
   separate task contract, resource table, proof map, and limits.
2. Use the [cost comparison](docs/QBP_COST_COMPARISON.md) and
   [bounded-input audit](docs/BOUNDED_INPUT_QBP.md) for precision regimes,
   classical construction, and explicit Pauli baselines. The
   [residual coefficient proof](docs/RESIDUAL_TABLE_PREPROCESSING.md)
   documents the classical interval helper; the
   [native rows and tables](docs/NATIVE_RESIDUAL_ROTATION.md) connect it to
   exact masks, borrowed-signal rotations, and enabled two- and four-row tables.
   The [bounded state integration](docs/NATIVE_RESIDUAL_STATE.md) composes
   two tables and amplification with the actual inverse for one or two
   system qubits, including the enlarged initial reflection.
   The [coherent selector](docs/NATIVE_RESIDUAL_BRANCH.md) adds an arbitrary
   protocol branch to the one-system-qubit residual preparation.
   The [native residual QBP integration](docs/NATIVE_RESIDUAL_QBP.md)
   connects that selector to the actual coarse circuit, controlled
   observable, and both raw gradient decoders for a fixed complex target.
   The [coverage map](docs/VERIFICATION.md#state-based-qbp-coverage) is the
   current entry point for what is proved, implemented, and optional.
3. For the separate frame question, read the
   [research status](docs/OPEN_PROBLEM.md),
   [grouped compiler](docs/CONDITIONAL_SUFFIX_COMPILER.md), and
   [retained stopping rule](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
4. Use [verification](docs/VERIFICATION.md), [attribution](docs/SOURCE_MAP.md),
   and [related work](docs/RELATED_WORK.md) before extending a claim.
   [README](README.md) and [REVIEW](REVIEW.md) retain the established frame
   results; the existing [publication scope](manuscript/PUBLICATION_SCOPE.md)
   is unchanged.

## Current results and their proof homes

| Result | Status and proof |
|---|---|
| Prescribed complete frame, exact gates | Matching size and CNOT/depth tradeoffs for every clean-work budget; [exact theorem](docs/COMPILER_THEOREM.md) |
| Prescribed complete frame, Clifford+T | Matching T-count under the sufficient-clean reservation; [fault-tolerant theorem](docs/FAULT_TOLERANT_COMPILER.md). Smaller clean budgets retain the [borrowed](docs/BORROWED_WORKSPACE_COMPILER.md) and [one-clean grouped](docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) bounds |
| State-based real and complex QBP | Two compiler flags, fine state preparation, an exact-return coarse word, and corrected magnitude/phase decoders; [consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) |
| Additional dirty banks | Improve the state preparation T bound while charging the exact coarse circuit; [banked proof](docs/COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) |
| State-based T-depth | Two complete schedules, one retaining the sharper count at a stronger dirty reservation; [depth proof](docs/STATE_QBP_DEPTH.md) and [fair comparison](docs/QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison) |
| Fixed-accuracy complete-frame T-depth | Matching count and depth in the same circuit at modest dirty width; [amortized tradeoff](docs/AMORTIZED_DIRTY_LOOKUP.md) |
| Bounded-input construction | Polynomial construction for the listed grouped/state alternatives, explicit program output, and separate fine-search caveats; [computational audit](docs/BOUNDED_INPUT_QBP.md) |
| Implemented evidence | Certified residual rows, bounded preparations, and [complete bounded residual QBP streams](docs/NATIVE_RESIDUAL_QBP.md), alongside the earlier exact-target examples; [claim-to-proof coverage](docs/VERIFICATION.md#state-based-qbp-coverage) |

The theorem's original accuracy bits K and state precision
$`P=\max(n,K)`$ must remain distinct. The observable needs accuracy K,
not the state's dimension floor. Compiler flags, the interference branch,
initialized system, and observable work are separate reservations. The
coarse word returns dirty work exactly; fine preparation includes all work
return and leakage in its initialized-isometry error.

At fixed accuracy the original frame route retains the better available
bound. Additional banks improve the state-only T expression in the stated
high-precision regimes, but sampling still grows as $`4^K`$ up to its
confidence factor. Clifford order and classical decoding order do not
improve. Explicit Pauli inputs also admit deterministic classical gradients
and term-only classical sampling. These results establish a task-specific
compiler improvement, not a general end-to-end gradient speedup.

## Current depth frontier

The selected next-step recommendations were accepted: retain optimal-order
T-count in the same circuit, and treat a logarithmic depth gap as a useful
bounded milestone. The [amortized indicator proof](docs/AMORTIZED_DIRTY_LOOKUP.md)
now closes that gap for prescribed complete real frames at fixed L. With two
clean flags and $`b\ge17(L+n+7)`$, one circuit has

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O\!\left(N/b^2+n^2\right).
```

At $`b=\Theta(n)`$ above that threshold, one circuit has optimal-order
$`T=\Theta(N/n)`$ and $`D_T=\Theta(N/n^2)`$.
The matching depth tradeoff $`D_T^\star=\Theta(N/b^2)`$ holds throughout
$`17(L+n+7)\le b\le\sqrt N/n`$ when this interval is nonempty.
The loader places one low-address indicator echo around a multiplexed
family of linear shears. A two-pass dirty traversal and rank-reduced
controlled shears amortize both former per-batch logarithms.
Dirty selector stacks remain separate from live banks and indicators;
all lookup work returns exactly. Full-frame error and QBP substitution
are inherited unchanged. A [capped precision allocation](docs/AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
now reduces accumulated source depth to $`O(nL+n\log(n+1))`$ without
increasing any layer's width or the weighted lookup costs. It retains
the full error guarantee using the sharper local isometry constant.
The total assigned precision is within $`5n`$ bits of the optimum under
that additive error certificate; this is not a frame-depth lower bound.
Large-workspace depth still has the additive $`n^2`$ routing contribution.
The following general-precision and state schedules retain
their separate contracts; this pass does not extend the new bound to them.

The [state-based depth theorem](docs/STATE_QBP_DEPTH.md) closes the selected
composition audit for both real and gauge-fixed complex Hopf QBP. Put
$`B_0=P+n+7`$. With the same two compiler flags, its schedules give:

| Dirty reservation | Compiler T-depth | Compiler T-count |
|---|---|---|
| $`b\ge2B_0`$ | $`O(NP/b+P+n^3)`$ | $`O(NP)`$ |
| $`b\ge16(B_0+\sqrt{NP})`$ | $`O(P+n^3)`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ |

Both have $`G=O(NP)`$ and charge the actual coarse C, its inverse in
magnitude readout, the fine residual table, occupied flags, predicates,
and all returned work. Each row's bounds hold for the same circuit.
The [complete task ledger](docs/QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison)
adds observable depth at precision K and all S or $`2S`$ executions.
A sufficiently large common pool permits an improved real-chart depth
upper expression even where the original and state T-count orders agree.
The original real schedules and eligible complex alternatives remain
in that comparison.

The routed-indicator refinement conjugates a single X by the existing
exact bank router. Its one-hot XOR uses no extra scratch and has linear,
rather than quadratic, address T-depth. This improves the state B bound
from $`O(P+n^4)`$ to $`O(P+n^3)`$ at the same sufficient reservation.
The complete real-frame schedules also use the established two-dirty-helper
MCX construction while query storage is idle. They now give

```math
D_T=O\!\left(\frac{NL}{b}
+\min\{nL+n^2,L\ell_*(n)+n^3\}\right),\qquad T,G=O(NL),
```

at $`b\ge2(L+n+7)`$. With the stronger sufficient pool
$`b\ge C(L+n+7+\sqrt{NL})`$, the same minimum depth expression
holds without the $`NL/b`$ term and with the sharper count
$`T=O(\sqrt{NL}+L\ell_*(n))`$. At fixed L this is
$`T=O(\sqrt N)`$ and $`D_T=O(n^2)`$ in one circuit, improving the
previous cubic depth bound. Literal routed-indicator/query checks include
arbitrary dirty inputs, exact phases, and the actual inverse orientation.

These remain upper schedules. The [source-depth certificate](docs/SOURCE_T_DEPTH.md)
parallelizes the paired source's two tails and proves matching exact depths
within the specified Majorana-layer architecture. This changes constants,
not the asymptotic Hopf bounds. It is not a lower bound for arbitrary
Clifford interlayers, controlled sources, or approximate replacement sources.
The [current literature audit](docs/RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026)
records why shallow clean-workspace synthesis and supplied catalysts do not
directly replace the dirty source. Large-width and general-precision T-depth
optimality and the selected $`b=N+n+7,L=N`$ complete-frame endpoint remain open. Obtaining sublinear
depth for these exact uncontrolled sources requires leaving the certified
architecture or changing the target. Do not repeat tail rescheduling or
denominator checks as an unrestricted lower bound.

## Completed native residual components

The bounded implementation sequence is complete. Detailed identities,
literal gate counts, certificates, and test scopes remain in their
canonical chapters; the handoff does not duplicate those ledgers.

| Component | Canonical implementation and scope |
|---|---|
| Certified coefficient programming and borrowed-signal rotations | [Native residual rows](docs/NATIVE_RESIDUAL_ROTATION.md): rational intervals to literal native words with full-operator error |
| Enabled two- and four-row tables | [Native table proof](docs/NATIVE_RESIDUAL_ROTATION.md): exact inactive identity, literal address phases, and no extra helper at these arities |
| One- and two-system-qubit preparation | [State integration](docs/NATIVE_RESIDUAL_STATE.md): two flags, actual-inverse amplification, all dirty-input return error, and the charged 28-T enlarged reflection |
| Coherent reference/target selection | [Branch integration](docs/NATIVE_RESIDUAL_BRANCH.md): arbitrary branch, preserved relative phase, and branch-independent initial reflection |
| Both complete gradient streams | [Residual QBP integration](docs/NATIVE_RESIDUAL_QBP.md): certified complex one-qubit target, actual coarse/inverse and controlled observable, all-outcome readout, and exact rational histogram decoders |

The small native preparations use q+2 arbitrary dirty wires. At q=L+10,
this exceeds the basic n=1 allocation by four wires and the basic n=2
allocation by three; these emitters fit the banked pool and do not replace
the minimum-budget fallbacks. All initialized system, compiler-flag, and
protocol-branch inputs stay separately counted. No intermediate reset,
postselection, or ideal-inverse substitution is allowed.

The [coverage audit](docs/VERIFICATION.md#state-based-qbp-coverage)
separates uniform proofs from the numerical and exact finite evidence.
In particular, exact rational decoders for a fixed fixture do not turn
the general floating-point contraction utilities into certified decoders.
The existing bounded integration supplies neither a same-accuracy native
fine-frame benchmark nor an end-to-end advantage claim.

## Remaining work and continuation criteria

The task-specific theorem, real/complex decoding, quantum resource ledger,
bounded-input construction, selected depth composition, and coverage
consolidation are complete. These are not pending research tasks. The
software extensions below are not prerequisites for the stated theorem.

| Remaining question | Concrete boundary |
|---|---|
| Certified bounded-input front end | Optional software: admitted Hopf inputs to an actual coarse word and certified residual data. The bounded-input proof already supplies construction and output bounds with its explicit small-system exception |
| Variable-size native schedule | Optional software: general tables, predicates, reflections, and banked count/depth scheduling. The bounded component-to-gradient pass is complete |
| General guarded decoder | Optional software: replace the general floating-point contractions with the proved certified arithmetic. Exact fixture decoders cover only their fixed target |
| Constant-clean complete-frame endpoint | Still open independently of the state-based task. A new candidate must supply an explicit complete native identity and symbolic precision/workspace ledger before another fixture pass |
| Optimal T-depth | Fixed-accuracy depth is matching for $`17(L+n+7)\le b\le\sqrt N/n`$; the large-width $`n^2`$ upper floor and general precision remain unresolved |

Do not repeat the completed coefficient-to-row, two- and four-row lookup,
one- and two-system-qubit amplification, or coherent residual-selection passes.
The selected bounded residual-to-QBP integration is also complete. No
further fixture expansion or general software API is selected. Begin a
future pass at the coverage map: specify either a concrete Hopf-QBP
input/output requirement and its missing interface, or a new scientific
claim and its proof obligation. A general compiler is not required to
close this selected Hopf-QBP task.
The guarded-batch logarithm is now amortized; do not repeat that task.
The next bounded depth target is fixed accuracy at sufficiently large
$`b=\Theta(\sqrt N)`$, retaining optimal-order $`T=\Theta(\sqrt N)`$.
Its available depth bounds are still $`\Omega(1)`$ and $`O(n^2)`$.
The capped precision allocation has removed the accumulated source widths
as a quadratic contribution: sources and suffix predicates now cost
$`O(n\log(n+1))`$ depth at fixed L. The retained per-layer routing
still sums to $`O(n^2)`$. An improvement must price the complete queries,
not only their rank-reduced shears or a single routed batch.

Two attempted shortcuts do not yet supply such a circuit. Keeping an
address-controlled router open conjugates an intervening table matrix A
to $`AP_x`$ or $`P_xA`$; this is another address-dependent operation
that must be implemented. Running independent equality tests in parallel
requires read-only shared address controls and exactly returned private
dirty work. The audited Khattar–Gidney logarithmic-depth dirty MCX circuit temporarily
borrows its controls, so its row instances cannot simply overlap them.
Dirty CNOT copies also carry unknown masks. These are gaps in the proposed
schedules, not unrestricted impossibility results. A concrete next route
is a shallow exact dirty-indicator batch with a full live-width ledger.
Polynomial overhead could be affordable on early layers, whose table sizes
decay geometrically, while retaining the existing queries near the leaves;
that hybrid remains conditional on the missing batch construction.
Before another native fixture pass, require either a complete all-input
construction with its count/depth/live-width ledger, or a lower bound
in the unrestricted T-depth model. Rescheduling the certified exact source
does not answer this question. The completed modest-width matching theorem
does not require solving the high-precision endpoint.
For any new component, keep literal phases and actual inverses, declare
all initialized inputs, include borrowed-work return in its isometry
error, and charge each reflection before composing a larger state compiler.

For the selected complete real-frame endpoint,

```math
N=2^n,\qquad a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

```math
\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n)),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

No new construction is selected for this gap, and recent task-specific
results do not show that its resolution is close. The established frame
results and the changed-decoder QBP theorem do not require its closure.
Do not replace the current task by arbitrary-unitary synthesis or invent
an uncharged observable, initialized history, or coherent evaluator.

## Completed routes that should not be repeated

The [canonical grouped audit](docs/CONDITIONAL_SUFFIX_COMPILER.md) retains
literal inverse branches, occupied flags, and full-work return, but still
pays precision per group. Its paired-source consolidation also improves
the old scalar baseline, so it is closed as the proposed amortization
mechanism. The [source-carry analysis](docs/SOURCE_REUSE_LIMITS.md) solves
width transport without pricing a cheaper joint interior program.

The [Hopf scattering step](docs/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
has one precision charge, but its unchanged-coin feedback conversion needs
a growing number of queries at fine accuracy. This is a query restriction,
not an additive native T lower bound. Globally programmed coefficients and
fully charged basis changes remain eligible for a different construction.
The [selection audit](docs/OPEN_PROBLEM.md#selection-audit-after-canonical-completion)
records the remaining complete-frame requirements.

The earlier [leaf-reference decoder](docs/REFERENCE_STATE_QBP.md) remains
an option under its explicit sampling tradeoff. The current coarse-frame
decoder removes that angle-dependent factor. Its general sampling bound,
complex gauge, and bounded-input preprocessing are settled; repeating a
special exact fixture would not strengthen those analytic claims.

## Restore and verify

Run from the repository root with `requirements.txt` installed:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python scripts/complex_coarse_native_example.py
python scripts/residual_qbp_native_example.py --q 16
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

Merge gates include Python 3.11/3.13 validation, all four exact-receipt
suites, and rendered presentation. The [verification map](docs/VERIFICATION.md)
distinguishes analytic proofs, finite native examples, exact rational
certificates, and floating-point decoder checks. Internal review and
finite tests are not external peer review or asymptotic proof.

No human action is required to resume repository work. Preserve the theorem
contracts and record any genuinely new question in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md), with its proof in the appropriate
chapter, so continuation does not depend on an old chat.
