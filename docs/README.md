# Documentation map

[Landing page](../README.md) · [Read the argument](../REVIEW.md) · [Publication scope](../manuscript/PUBLICATION_SCOPE.md)

This is the **proof-dependency map**, not another theorem ledger.
The [publication scope](../manuscript/PUBLICATION_SCOPE.md) owns the selected
A–D statements and literal reservations; the [resource frontier](OPEN_PROBLEM.md)
owns the strongest research bounds and remaining general question.
For active work, start with [WORKSPACE](../WORKSPACE.md). The
[research catalogue](../research/README.md) preserves individual constructions,
restricted results, failed mechanisms and their proof homes.

## Proof chapters

### Shared target and error contract

The prescribed operator includes state and differential-frame columns,
literal phases, every logical input, arbitrary dirty/reference inputs,
and clean leakage. Every actual inverse and return error is charged.

| Proof home | Role |
| --- | --- |
| [Hopf interface](HOPF_INTERFACE.md) | Coordinates, marker columns, singular charts, addressed layers, tangent metric and operator stability |
| [Frame-safe compilation](FRAME_SAFE_COMPILATION.md) | Complete-input error, actual adjoints, fixed-decoder necessity and singular exceptions |
| [Compiler boundaries](COMPILER_BOUNDARIES.md) | State preparation, weaker checkpoints and prescribed complete frames are different targets |

### A. Exact compilation at every clean-workspace budget

The [exact compiler theorem](COMPILER_THEOREM.md) owns the matching size,
CNOT-count and total-depth statements for every clean allocation, for real
and phase-dressed complex magnitude frames. Its gate model is arbitrary
one-qubit gates and CNOTs, distinct from Clifford+T.

| Step | Canonical argument |
| --- | --- |
| No added work | [Borrowed-suffix echo](COMPILER_THEOREM.md#4-strict-zero-workspace), retaining all borrowed-bit sectors |
| Small positive work | [Direct flagged schedule](COMPILER_THEOREM.md#5-small-positive-workspace) |
| Larger work | [Tree cut](COMPILER_THEOREM.md#6-exact-tree-cut), [binary–one-hot decoder](COMPILER_THEOREM.md#7-conditioned-prefix-construction) and [routed parallel tails](COMPILER_THEOREM.md#8-routed-parallel-tail) |
| Optimality and phases | [Matching lower bounds](COMPILER_THEOREM.md#9-matching-lower-bounds) and [literal phase dressing](COMPILER_THEOREM.md#10-phase-dressed-complex-magnitude-frame) |

### B. Matching T-count with sufficient clean workspace

The [fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md) owns the sufficient-clean
precision/workspace frontier. Follow the full-input source, residual dictionaries,
retained failure history and full-output normalization as one charged circuit.

| Step | Canonical argument |
| --- | --- |
| Coarse/fine separation | [Grouped dictionaries](FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame) |
| Shared fine precision | [Reusable source and correction kernel](FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel) |
| Recovery and return | [Full-output normalization](FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization) |
| Lower bounds | [Family reduction, coherent counting and literal-phase witness](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage) |

The lower bounds also support C and D. Neither the sufficient clean allocation
nor sums of source-call costs establish a clean-space or unrestricted depth
lower bound.

### C. One-clean compilation and its corollaries

The [growing-source theorem, Section 11](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
owns C; [Section 12](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
owns its banked extension. The low-precision baseline and high-precision
source meet uniformly, including bounded dimension. The grouped theorem
is a retained corollary, not a premise of that splice.

| Dependency | Required role |
| --- | --- |
| [One-clean source and layerwise baseline](ONE_CLEAN_COMPILER.md#6-complete-real-frames-and-the-grouped-extension) | Exact conditioning, full-output amplification, low-precision banks and bounded-dimension fallback |
| [Growing sources](../research/endpoint/GROWING_SOURCE_COMPILER.md) | Precision/mode separation, sealed modes, actual dirty replay and compensated overlap |
| [Conjugate-signal recovery](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md#2-lift-the-whole-high-level-word) | Borrow the phase signal while retaining the genuine initialized anchor and paying native interlayers |
| [Banked block/source interfaces](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces) | Sum native costs at the stated reservations |
| [Conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally) | Exact inactive-work guards; its grouped count construction remains separately scoped |
| [Literal diagonals and U(2) multiplexors](ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors) | Literal phase composition and separate dirty reservations |
| [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) | Full-input two-clean source, packed diagonal and amplification/error certificate |
| [Borrowed-workspace compiler](BORROWED_WORKSPACE_COMPILER.md) | Mixed clean/dirty allocations and their restricted matching splice |

Reducing clean work does not supersede every earlier allocation.
Keep the base and doubled-base thresholds and the exact-base precision
qualification from the claim ledger. C adds no T-depth theorem.

<a id="d-simultaneous-real-frame-t-count-and-t-depth"></a>

### D. Simultaneous frame T-count and T-depth

**Selected Result D remains R79.** It combines the native R78 decoder
with the R72 terminal and fallback schedules, at the respective real/complex
one-clean dirty reservations. Counts and depth belong to the same circuit;
the full-isometry composition includes earlier leakage and literal phases.

| Selected dependency | Role |
| --- | --- |
| [Uniform subtree depth, R79](../research/depth/UNIFORM_SUBTREE_DEPTH.md) | Absolute all-precision source/predicate costs, protected suffix allocation and all-width composition |
| [Subtree-projector decoder, R78](../research/depth/SUBTREE_PROJECTOR_DECODER.md) | Native transported reflections and coherent spectral-label erasure |
| [One-clean depth, R72](../research/depth/ONE_CLEAN_DEPTH.md) | Protected terminal flag, width/precision fallback and one-clean literal phase composition |
| [High-precision branch, R69](../research/depth/HIGH_PRECISION_DEPTH.md) | Stronger additive-precision bound in its stated regime |
| [Uniform baseline](UNIFORM_PRECISION_DEPTH.md) | Rectangular allocation, cutoff, weighted sums and literal-width fallback |
| [Retained real R70](../research/depth/ALL_PRECISION_DEPTH.md) and [complex R71](../research/depth/COMPLEX_FRAME_DEPTH.md) | Separate two-clean proofs and their original allocations |

**R96 is the strongest unconditional research upper bound.** It is not
silently substituted for selected Result D. Its current statement, precision
range and open optimality gap are in the [resource frontier](OPEN_PROBLEM.md).

| R96 proof chain | Role |
| --- | --- |
| [All-precision linear depth, R96](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md) | Final general theorem and original-width resource/error closure |
| [Global borrowed-precision stage](../research/depth/GLOBAL_BORROWED_PRECISION_STAGE.md) and [fast capped reference](../research/depth/FAST_CAPPED_GEOMETRIC_REFERENCE.md) | Separate initialized precision scratch from coarse geometry; compare on the full space |
| [Polynomial-precision bridge](../research/depth/POLYNOMIAL_PRECISION_DEPTH.md) | Three charged coarse contractions |
| [Direct-prefix terminal](../research/depth/DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md) | Fixed fine-mask sequence and explicitly programmed shorter native coarse word |

The following primitives remain dependencies or distinct fallback interfaces,
not new tasks to repeat. Their individual hypotheses are part of the proofs.

| Interface | Proof home |
| --- | --- |
| Dirty table selection | [Blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md), [parallel word-query refinement R91](../research/depth/PARALLEL_BILINEAR_WORD_QUERY.md), [selected routing](../research/depth/SELECTED_DIRTY_ROUTING.md) |
| Early precision and program lifetime | [Unary phase-source groups](UNARY_PHASE_GRADIENT.md), [grouped prefetch](GROUPED_PROGRAM_PREFETCH.md), [conditional geometric predecessor](CONDITIONAL_GEOMETRIC_SOURCE.md) |
| Tail and literal-width fallback | [Chunked indicators](CHUNKED_DIRTY_INDICATOR.md), [dirty-sum compression](DIRTY_SUM_COMPRESSION.md), [amortized lookup](AMORTIZED_DIRTY_LOOKUP.md), [parallel lookup](PARALLEL_DIRTY_LOOKUP.md) |
| Native scheduling and errors | [T-depth schedule](T_DEPTH_COMPILER.md), [Hopf angular error accumulation](HOPF_ERROR_ACCUMULATION.md) |
| Earlier linear-precision chain | [Complete early operator](../research/depth/SUBLINEAR_T_DEPTH.md), [descendant loader](../research/depth/COHERENT_DESCENDANT_LOADER.md), [shared gradient](../research/depth/SHARED_GRADIENT_RECOVERY.md), [dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md), [binary terminal](../research/depth/SHARED_BINARY_TERMINAL.md), [whole-word prefetch](../research/depth/GROUPED_TERMINAL_DEPTH.md#2-prefetch-the-exact-coarse-native-word), [nonuniform indicators](../research/depth/NONUNIFORM_DIRTY_INDICATOR.md) |

The [research catalogue](../research/README.md#large-workspace-t-depth) indexes
intermediate R80–R95 constructions and their distinct promised families,
catalogues and workspace-transfer components. Stronger headline bounds do
not replace those component proofs. The [constant-depth frontier](OPEN_PROBLEM.md#constant-depth-frontier)
separates general lower bounds from midpoint, spread and outer-overlap
restrictions. Source normalization, logical calls and representation
obstructions are not unrestricted T-depth bounds.

### Fixed-decoder QBP consequence

[Exact QBP substitution](QBP_CONSEQUENCE.md) and [approximate QBP](QBP_APPROXIMATION.md)
retain the original inverse-frame decoder. They charge executions,
observable access, initialized protocol work and classical output separately.
They establish neither optimality among all gradient algorithms nor
end-to-end application advantage.

<a id="separate-supplements-and-optional-research"></a>
<a id="separate-endpoint-theorem"></a>
<a id="selected-endpoint-count-theorem-and-complex-magnitude-corollary"></a>
<a id="selected-precision-uniform-two-clean-count-theorem"></a>
<a id="selected-growing-source-count-theorems"></a>

### Retained two-clean count corollaries

[Growing-source Section 10](../research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
retains precision/mode separation. [Section 9's complex endpoint](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
keeps its distinct allocation: two clean qubits, one fewer dirty qubit
than the one-clean complex theorem at the stated high-precision endpoint.
The [claim ledger](../manuscript/PUBLICATION_SCOPE.md) retains its exact conditions.

### Completed changed-decoder supplement

The [state-based QBP supplement](../supplements/state_based_qbp/README.md)
has its own initialized state-preparation interface, real/complex theorem,
task costs and native checks. It does not replace the complete-frame
contract or add a prerequisite to A–D.

## Evidence and sources

| Read | Purpose |
| --- | --- |
| [Core-claim audit](CORE_CLAIM_AUDIT.md) | Internal claim-to-proof review and its evidence limits |
| [Verification](VERIFICATION.md) and [fault-tolerant receipts](../verification/fault_tolerant/README.md) | Proof-to-code coverage, exact certificates, native finite checks and reproduction |
| [Source map](SOURCE_MAP.md), [related work](RELATED_WORK.md) and [provenance](../provenance/README.md) | Imported premises, local contributions and source versions |
| [Implementation](../compiler_robust_hopf/README.md), [tests](../tests/README.md), [scripts](../scripts/README.md) and [diagrams](../assets/README.md) | Executable and rendering maps |
| [Manuscript guide](../manuscript/README.md) | Selected argument's section order |

Finite checks do not prove asymptotic theorems, establish priority or supply
external peer review. No general elementary emitter for the asymptotic
shared-source compiler is claimed.
