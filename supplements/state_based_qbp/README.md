# Completed state-based Hopf QBP results

[Main compiler package](../../README.md) · [Publication scope](../../manuscript/PUBLICATION_SCOPE.md) · [Original QBP contract](../../docs/QBP_APPROXIMATION.md)

This completed supplementary package estimates the original raw Hopf gradients
using fine state preparation, an actual coarse reference and inverse, and a
changed classical decoder. It includes real and phase-dressed complex charts,
both complex-gradient streams, charged quantum and classical costs, and bounded
native integrations. These are separate task-level results outside the selected
Results A–D manuscript scope.

The package does not prescribe the fine circuit's other logical columns.
Accordingly, it neither closes the complete-frame endpoint nor changes the
[necessity statement for the fixed inverse-frame decoder](../../docs/FRAME_SAFE_COMPILATION.md).
The mathematical construction and selected bounded integrations are complete.
A general certified front end, variable-size native emitter, and general guarded
decoder implementation remain optional software extensions.

## Read the result and its proof

Start with the [consolidated theorem](STATE_BASED_QBP_THEOREM.md). Follow state
preparation to the actual coarse interface, then the real and complex decoder
proofs. The cost and bounded-input chapters state what is charged and what the
classical input assumptions permit.

| Chapter | Role |
|---|---|
| [STATE_BASED_QBP_THEOREM.md](STATE_BASED_QBP_THEOREM.md) | Complete task, accuracy, workspace, confidence, and evidence contract |
| [STATE_ONLY_COMPILER.md](STATE_ONLY_COMPILER.md) | Fine real-state preparation and coherent reference selection |
| [COMPLEX_COARSE_COMPILER.md](COMPLEX_COARSE_COMPILER.md) | Gauge-fixed complex coarse circuit, exact dirty return, complex preparation, and banked refinement |
| [COARSE_FRAME_QBP.md](COARSE_FRAME_QBP.md) | Actual-coarse X/Y scores, bounded depth records, and real histogram reconstruction |
| [COMPLEX_COARSE_QBP.md](COMPLEX_COARSE_QBP.md) | Original magnitude and leaf-phase gradients with certified reconstruction budgets |
| [REFERENCE_STATE_QBP.md](REFERENCE_STATE_QBP.md) | Separate leaf-only alternative with its own reference and sampling tradeoffs |
| [STATE_QBP_DEPTH.md](STATE_QBP_DEPTH.md) | Completed same-circuit T-count and T-depth upper schedules |
| [QBP_COST_COMPARISON.md](QBP_COST_COMPARISON.md) | Common-accuracy and workspace comparison with the original frame programs |
| [BOUNDED_INPUT_QBP.md](BOUNDED_INPUT_QBP.md) | Constructive preprocessing, explicit output costs, and classical Pauli baselines |
| [RESIDUAL_TABLE_PREPROCESSING.md](RESIDUAL_TABLE_PREPROCESSING.md) | Certified algebraic rotation coefficients and the small-system banked construction |

The proofs reuse the main [borrowed-workspace interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md),
[one-clean rotation construction](../../docs/ONE_CLEAN_COMPILER.md),
[dirty-bank query](../../docs/OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost),
and [parallel lookup](../../docs/PARALLEL_DIRTY_LOOKUP.md).
The [original QBP approximation proof](../../docs/QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights)
supplies the shared observable-access and reused-dirty-work concentration
premises. Results A–D do not rely on this supplementary decoder.

## Bounded native evidence

These notes specify implemented components and complete small integrations.
Their contracts, exact certificates, and literal gate ledgers complement the
[verification map](../../docs/reference/VERIFICATION_CATALOGUE.md#state-based-qbp-coverage);
finite checks do not establish the asymptotic construction by themselves.

| Evidence note | Implemented scope |
|---|---|
| [NATIVE_RESIDUAL_ROTATION.md](NATIVE_RESIDUAL_ROTATION.md) | Certified unaddressed residual rows and enabled two- and four-row tables |
| [NATIVE_RESIDUAL_STATE.md](NATIVE_RESIDUAL_STATE.md) | One- and two-qubit preparation, two flags, and actual-inverse amplification |
| [NATIVE_RESIDUAL_BRANCH.md](NATIVE_RESIDUAL_BRANCH.md) | One-qubit reference/target selection preserving an arbitrary branch and its relative phase |
| [NATIVE_COARSE_QBP.md](NATIVE_COARSE_QBP.md) | Complete real-target two-qubit integration using the exact finite-size fallback |
| [NATIVE_COMPLEX_COARSE_QBP.md](NATIVE_COMPLEX_COARSE_QBP.md) | Both complex-gradient streams, actual prefix rows, and a certified coarse word |
| [NATIVE_RESIDUAL_QBP.md](NATIVE_RESIDUAL_QBP.md) | Fine residual preparation inside both streams for a fixed complex one-qubit target |

The [implementation index](../../compiler_robust_hopf/README.md),
[tests](../../tests/README.md), and [example scripts](../../scripts/README.md)
provide executable entry points. The notes retain their precise small-system
workspace exceptions and distinguish floating-point checks from exact rational
certificates and analytic error bounds.

## Contracts retained throughout

Accuracy bits K remain separate from the state precision floor
$`P=\max\{n,K\}`$. The two compiler flags are additional to the initialized
system, interference branch, and observable work. The complex preparation uses
one consistent arithmetic-mean phase gauge; it does not compile the prescribed
frame's literal common phase. Actual coarse words and their actual inverses are
used in both the quantum circuit and decoder.

Error includes clean leakage and arbitrary dirty/reference inputs. No QRAM,
supplied precision source, intermediate reset, or postselection is introduced.
The tuple, decoder, and observable are fixed during an execution batch. All
executions, observable access, preprocessing, arithmetic precision, and output
costs remain charged. The results assert neither gradient-cost optimality nor
an end-to-end quantum advantage. Further work requires a concrete missing
interface or a new claim; enlarging a finite fixture is not a readiness gate.
