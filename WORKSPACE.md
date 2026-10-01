# Continuing research workspace

This is the entry point when a previous conversation or execution workspace
is unavailable. Proofs and decisions live in the repository.

The 2026-10-01 direction review starts from verified main
`204a78ee40adc8eedc07939f0a9b7082d750a6a2`, which merged the consolidated
theorem. Check the current branch and later commits before continuing.
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
   documents the implemented classical interval helper.
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
| Bounded-input construction | Polynomial construction for the listed grouped/state alternatives, explicit program output, and separate fine-search caveats; [computational audit](docs/BOUNDED_INPUT_QBP.md) |
| Implemented evidence | Complete bounded [real](docs/NATIVE_COARSE_QBP.md) and [complex](docs/NATIVE_COMPLEX_COARSE_QBP.md) native examples, histogram decoders, and certified classical residual coefficients; [verification map](docs/VERIFICATION.md) |

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

## Selected next pass: depth of the state-based gradient circuit

The mathematical Hopf QBP task is established. The next bounded question
is whether its two-flag state construction admits an explicit dirty-workspace
versus T-depth schedule for the complete gradient execution. This follows
the original ancilla–depth motivation and has a concrete starting point in
the existing [routing](docs/T_DEPTH_COMPILER.md) and
[parallel dirty lookup](docs/PARALLEL_DIRTY_LOOKUP.md) lemmas.
Those chapters currently compose schedules for complete real frames;
they do not already prove the corresponding state-based real/complex result.

The [selected audit and stopping rule](docs/OPEN_PROBLEM.md#selected-state-based-depth-audit)
requires the actual coarse C and its inverse, the fine residual table,
occupied flags, predicates, all returned work, both gradient streams, and
the observable's depth to be charged. A T-count improvement alone does
not establish a depth improvement. Depth-optimized and count-preserving
schedules must have separate resource ledgers, with all simultaneous
claims applying to the same circuit.

This is a bounded composition audit, not a new general lookup or
state-preparation program. Its success condition is a proved task-level
upper bound and a fair comparison with eligible original schedules. If
a required operation cannot be scheduled with the declared work, record
the precise obstruction and stop that route. Matching depth lower bounds,
source parallelization, and a general native emitter are outside this pass.

## Remaining work and continuation criteria

The task-specific theorem, real/complex decoding, quantum resource ledger,
and bounded-input construction are established. Consolidation is complete;
these are not pending research tasks.

| Remaining question | Concrete boundary |
|---|---|
| State-based T-depth | Selected next: audit the existing schedules on the complete real/complex state-based gradient circuit; no new depth bound is claimed by this revision |
| General fine-precision native emitter | The full gate-emission package is not implemented. The certified residual helper outputs coefficient intervals, and bounded native examples use special finite-size targets. Deferred while the depth audit runs; a later implementation pass must name the missing emitted primitive and its complete error/workspace contract |
| End-to-end algorithmic advantage | No example is selected. A new claim needs a concrete observable-access model and a classical comparator; explicit Pauli inputs do not supply the high-precision advantage claimed by T-count alone |
| Constant-clean complete-frame endpoint | Still open independently of the state-based task. A new candidate must supply an explicit complete native identity and symbolic precision/workspace ledger before another fixture pass |
| Depth optimality and practical constants | Remain separate from the selected upper-bound audit; optimal T-count does not imply optimal depth or a practical crossover |

If native emission is resumed, the first missing bridge is certified
cosine/sine intervals to exact paired-source sign masks and then a literal
borrowed-signal Ry/Rz word. Begin with one unaddressed residual U(z) row.
The present test helpers use floating-point rounding and some dense
matrices; they are not a certified, scalable implementation of that bridge.

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
