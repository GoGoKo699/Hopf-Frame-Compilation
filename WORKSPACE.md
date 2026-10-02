# Continuing research workspace

This is the entry point when a previous conversation or execution workspace
is unavailable. Proofs and decisions live in the repository.

The 2026-10-02 coverage-consolidation pass starts from verified main
`736c25500dbb01fc0b449aec1cf1fa0105837817`, after native residual
QBP integration. Check later commits before continuing.
The selected state-based Hopf QBP construction and its bounded-input audit are complete; the
[consolidated theorem](docs/STATE_BASED_QBP_THEOREM.md) is their entry point.
It prepares a state and changes the gradient decoder while retaining all
raw coordinates, singular angles, and both complex-gradient streams.
The prescribed complete-frame endpoint remains a separate open problem.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized. Keep manuscript writing
and release work outside this research pass. Use small analytic examples
and finite checks, without large simulations, QRAM, resets inside a
compiler execution, supplied catalysts, or hidden initialized work.

1. Read the [state-based QBP theorem](docs/STATE_BASED_QBP_THEOREM.md) for
   the current task, exact assumptions, resource table, proof map, and limits.
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

## Completed state-based depth pass

The [state-based depth theorem](docs/STATE_QBP_DEPTH.md) closes the selected
composition audit for both real and gauge-fixed complex Hopf QBP. Put
$`B_0=P+n+7`$. With the same two compiler flags, its schedules give:

| Dirty reservation | Compiler T-depth | Compiler T-count |
|---|---|---|
| $`b\ge2B_0`$ | $`O(NP/b+P+n^3)`$ | $`O(NP)`$ |
| $`b\ge16(B_0+\sqrt{NP})`$ | $`O(P+n^4)`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ |

Both have $`G=O(NP)`$ and charge the actual coarse C, its inverse in
magnitude readout, the fine residual table, occupied flags, predicates,
and all returned work. Each row's bounds hold for the same circuit.
The [complete task ledger](docs/QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison)
adds observable depth at precision K and all S or $`2S`$ executions.
A sufficiently large common pool permits an improved real-chart depth
upper expression even where the original and state T-count orders agree.
The original real schedules and eligible complex alternatives remain
in that comparison.

This completes the bounded pass using existing exact lookup and source
identities; no new circuit primitive or large simulation was needed.
These are upper schedules, not matching T-depth tradeoffs or total-runtime
claims. The source programs still run serially. Source parallelization
and depth lower bounds remain separate questions.

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
| End-to-end algorithmic advantage | No example is selected. A new claim needs a concrete observable-access model and a classical comparator; explicit Pauli inputs do not supply the high-precision advantage claimed by T-count alone |
| Constant-clean complete-frame endpoint | Still open independently of the state-based task. A new candidate must supply an explicit complete native identity and symbolic precision/workspace ledger before another fixture pass |
| Depth optimality and practical constants | Remain open after the completed upper-bound audit; optimal T-count does not imply optimal depth or a practical crossover |

Do not repeat the completed coefficient-to-row, two- and four-row lookup,
one- and two-system-qubit amplification, or coherent residual-selection passes.
The selected bounded residual-to-QBP integration is also complete. No
further fixture expansion or general software API is selected. Begin a
future pass at the coverage map: specify either a concrete Hopf-QBP
input/output requirement and its missing interface, or a new scientific
claim and its proof obligation. A general compiler is not required to
close this selected Hopf-QBP task.
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
