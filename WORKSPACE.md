# Continuing research workspace

This is the entry point when a previous conversation or execution workspace
is unavailable. Proofs and decisions live in the repository.

The 2026-10-02 native residual QBP pass starts from verified main
`ad312be7b4247f6534a810074263eb7d7e8fd27d`, after coherent residual
branch selection. Check later commits before continuing.
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
| Implemented evidence | Complete bounded [real](docs/NATIVE_COARSE_QBP.md) and [complex](docs/NATIVE_COMPLEX_COARSE_QBP.md) native examples, histogram decoders, [certified residual rows and bounded tables](docs/NATIVE_RESIDUAL_ROTATION.md), and [bounded residual state amplification](docs/NATIVE_RESIDUAL_STATE.md); [verification map](docs/VERIFICATION.md) |

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

The [native bridge](docs/NATIVE_RESIDUAL_ROTATION.md) now programs exact
paired-source signs from rational cosine/sine intervals, emits literal
borrowed-signal Ry/Rz words, and composes one unaddressed residual U(z)
row. Its full-operator error is below $`130\,2^{-q}`$ on arbitrary target,
core, signal, and reference inputs. It uses q+2 borrowed wires excluding
the target, no clean work, and O(q) gate storage. The unsimplified row has
540q T/TDG gates. Exact rational programming and small all-input native
checks replace the previous floating-point-only programming evidence.
The two-row extension adds one address bit and one enable literal. It
preserves both logical wires exactly and is exactly identity when
disabled. Source-center CNOT/Toffoli gates implement the predicate with
no extra helper at this control arity. Its full-operator row bound remains
below $`130\,2^{-q}`$, taking the maximum over addresses; it has q+5
total wires and the same q+2 borrowed work wires. The unsimplified table
has 540q+630 T/TDG gates. Complete small-input checks use four fixed
control sectors and retain their relative phases. These components
implement existing identities; they do not change the asymptotic
theorem or emit its complete state-preparation schedule.

The [one-system-qubit residual state emitter](docs/NATIVE_RESIDUAL_STATE.md)
now composes two enabled tables, an exact controlled Hadamard, and one
state-amplification step with the actual reversed word. Its full
initialized-isometry bound is below 390 times 2 to the power minus q,
including both returned flags and arbitrary borrowed/reference input.
The unsimplified count is 3240q+3793 T/TDG gates, with q+2 dirty wires.
At q=L+10 this bounded emitter exceeds the minimum n=1 dirty reservation
by four wires; it fits the banked pool and does not replace the
minimum-budget fallback. Exact reflections and the literal leading
minus sign are included. No coarse C or QBP branch is emitted here.

The [four-row lookup](docs/NATIVE_RESIDUAL_ROTATION.md#6-four-rows-with-two-address-bits-and-no-additional-helper)
adds a second address bit without any extra helper. Constant and linear
mask terms are Clifford; each nonempty same-axis quadratic support is
a Clifford conjugate of one exact seven-T Toffoli. The literal residual
count is 540q+630+70(2k_z+k_y), at most 540q+1050, with q+2 dirty
wires and q+6 total wires. Its full-operator error remains below
130 times 2 to the power minus q. Complete native words return both
addresses and enable exactly, including their phases; temporary
address changes inside native Toffolis are explicitly tracked in tests.
The disabled action remains exactly identity. This is a fixed four-row
component, not a general lookup or larger state-preparation schedule.

The [two-system-qubit state emitter](docs/NATIVE_RESIDUAL_STATE.md#5-two-system-qubit-preparation-and-a-returned-core-helper)
composes four-row tables and two controlled Hadamards. Its enlarged
initial reflection reuses core wire zero as an arbitrary helper, returns
it exactly on every input, and charges four seven-T Toffolis. The
complete word uses two clean flags and q+2 dirty wires, with error below
390 times 2 to the power minus q. Its literal count is
3240q+3820+210(2k_z+k_y), at most 3240q+5080. At q=L+10 the dirty
pool exceeds the minimum n=2 allocation by three wires and fits the
banked pool. The actual inverse retains all leakage; no reset or supplied
coarse C is used. General tables remain a separate integration task.

The [coherent residual selector](docs/NATIVE_RESIDUAL_BRANCH.md) now uses
the second address bit as an arbitrary protocol branch for one system
qubit. Its seven-T initial reflection excludes that branch; literal
relative phase and the actual inverse are retained. Both coefficient
normalizations are checked separately. The complete branch-and-dirty
isometry has error below 390 times 2 to the power minus q and count
3240q+3793+210(2k_z+k_y), at most 3240q+5053. It uses q+2 dirty
wires, one initialized system, two clean flags, and the separate arbitrary
branch. At q=L+10 the pool exceeds the basic n=1 allocation by four wires
and fits the banked pool. The common-coarse reference specialization uses
root one and tail zero on branch zero; no coarse C or observable is
emitted by this component.

The [native residual QBP fixture](docs/NATIVE_RESIDUAL_QBP.md) completes
that bounded integration. An exact rational coarse-distance certificate
places its complex one-qubit target inside 1/64. Certified square-root
intervals program its residual coefficients without floating angles.
Magnitude readout charges the coherent residual word, common C, the
controlled Hadamard observable, and actual C inverse; phase readout uses
the unbranched residual word, common C, and the same observable. Both
retain all flag/dirty leakage and have exact rational histogram decoders.
The extra T charges are six and four respectively, so S executions of
each stream cost S times (T_pair+T_one+10). The dirty allocation retains
the banked n=1 caveat. This is finite native evidence for the existing
theorem, without a general emitter or a same-accuracy fine-frame benchmark.

## Remaining work and continuation criteria

The task-specific theorem, real/complex decoding, quantum resource ledger,
bounded-input construction, and the selected depth composition are
established. These are not pending research tasks.

| Remaining question | Concrete boundary |
|---|---|
| General fine-precision native emitter | Certified rows, bounded tables and preparations, coherent residual selection, and a one-system-qubit native residual QBP integration with both decoders are implemented. General table sizing, larger predicates, and full variable-size fine emission remain separate. The next pass consolidates the proof-to-code map and checks which remaining implementation gaps affect any stated claim before selecting another component |
| End-to-end algorithmic advantage | No example is selected. A new claim needs a concrete observable-access model and a classical comparator; explicit Pauli inputs do not supply the high-precision advantage claimed by T-count alone |
| Constant-clean complete-frame endpoint | Still open independently of the state-based task. A new candidate must supply an explicit complete native identity and symbolic precision/workspace ledger before another fixture pass |
| Depth optimality and practical constants | Remain open after the completed upper-bound audit; optimal T-count does not imply optimal depth or a practical crossover |

Do not repeat the completed coefficient-to-row, two- and four-row lookup,
one- and two-system-qubit amplification, or coherent residual-selection passes.
The selected bounded residual-to-QBP integration is also complete. Do not
automatically enlarge its dimension or add more fixtures: first identify
a distinct interface or scientific claim that the next check would resolve.
For the next component, keep literal phases and actual inverses, declare
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
