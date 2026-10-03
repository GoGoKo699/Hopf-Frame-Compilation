# Documentation map

[Landing page](../README.md) · [Read the argument](../REVIEW.md) · [Publication scope](../manuscript/PUBLICATION_SCOPE.md)

For the gradient task, start with the
[state-based QBP theorem](STATE_BASED_QBP_THEOREM.md). For the complete-frame
compiler argument, start with the landing page and `REVIEW.md`. Each formal
topic has one primary chapter below.

## Proof chapters

### Frame interfaces and compiler proofs

| Chapter | Role |
|---|---|
| [Hopf interface](HOPF_INTERFACE.md) | Coordinates, marker columns, singular charts, and addressed layers |
| [Frame-safe compilation](FRAME_SAFE_COMPILATION.md) | Complete-input contract, operational necessity, sharp marker sensitivity, and actual adjoints |
| [Compiler boundaries](COMPILER_BOUNDARIES.md) | Concrete distinctions between state, checkpoint, and complete-frame promises |
| [Exact compiler theorem](COMPILER_THEOREM.md) | All clean-workspace budgets; echo, decoder, router, and matching size/CNOT/depth bounds |
| [Fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md) | Real and complex sufficient-clean matching T-count; shared-source residual proposition |
| [One-clean compiler](ONE_CLEAN_COMPILER.md) | Two programmable overlaps, native packed source, five-call amplification, and real/phase-dressed frame, diagonal, and multiplexor bounds |
| [Two-clean compiler](OPERATOR_SOURCE_COMPILER.md) | Exact source costs, frame/bank bounds, matched diagonal and general multiplexor frontiers, and certified preprocessing |
| [Conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) | Precision-uniform grouped bounds and their one-clean extension, ancestor-column residuals, and complete-input return |
| [T-depth schedule](T_DEPTH_COMPILER.md) | Explicit dirty-bank depth tradeoff, literal shared-control swaps, and the remaining lower-bound gap |
| [Exact source depth](SOURCE_T_DEPTH.md) | Parallel paired tails and matching exact source depths within a specified Majorana-linear circuit class; no unrestricted optimality claim |
| [Two-layer source obstruction](SHALLOW_SOURCE_OBSTRUCTION.md) | Width-independent approximation gaps for full-input sources with arbitrary Clifford interlayers and returned dirty helpers |
| [Conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md) | Logarithmic precision depth using the active logical suffix as temporary clean work; lookup and suffix-predicate costs remain separate |
| [Grouped program reuse](GROUPED_PROGRAM_PREFETCH.md) | Complete real-frame depth O(n log log n) at fixed accuracy and sufficient square-root dirty width; O(n) selector maintenance and scoped source-reuse audit |
| [Unary phase-source groups](UNARY_PHASE_GRADIENT.md) | Complete real-frame T-depth O(n) at fixed accuracy and sufficient square-root dirty width, retaining optimal-order T-count; charged preparation and exact guarded cyclic shifts |
| [Blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md) | Full-input selected bilinear blocks with returned dirty work; fixed-accuracy depth O(N/b²+n), optimal-order T-count, and matching range through sqrt(N/n) |
| [Uniform precision and depth](UNIFORM_PRECISION_DEPTH.md) | Absolute-constant depth O(NL/b²+nL); rectangular allocation gives O(NL/b²+n) for slowly growing precision with the same-circuit count guarantee |
| [Chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md) | A tunable exact dirty-tree indicator and a summable late-query budget that removes the late routing bottleneck |
| [Hopf error accumulation](HOPF_ERROR_ACCUMULATION.md) | Sharp ideal-angle stability, finite relative spectra, and coherent leakage in the actual shared-flag sources; scoped precision boundaries |
| [Flag-echo audit](HOPF_FLAG_ECHO.md) | Exact errors of four diagonal Pauli echoes, their generic linear leakage, and an exact equal-mask exception |
| [Radial source filter](HOPF_RADIAL_FILTER.md) | Phase-calibrated fixed-point filtering, quadratic radial error, charged native phases, and a smaller complete-frame precision cap |
| [Parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md) | Exact returned dirty indicators and simultaneous count-efficient, low-T-depth full-frame compilation |
| [Masked dirty-sum compression](DIRTY_SUM_COMPRESSION.md) | Deferred parity forests and overlapping carry blocks give logarithmic-depth dirty indicators and a sharper complete-frame bound |
| [Batched dirty lookup](BATCHED_DIRTY_LOOKUP.md) | Partial indicator batches give optimal-order fixed-accuracy frame count and a logarithmic depth gap at modest dirty width |
| [Amortized dirty lookup](AMORTIZED_DIRTY_LOOKUP.md) | Same-circuit count and depth bounds at every accuracy, matching in an explicit workspace range that includes inverse-polynomial error |
| [Borrowed-workspace appendix](BORROWED_WORKSPACE_COMPILER.md) | Exact dirty lookup and predicates; arbitrary-budget bound and restricted matching splice |

### Gradient protocols and task costs

| Chapter | Role |
|---|---|
| [State-based QBP theorem](STATE_BASED_QBP_THEOREM.md) | Task-level statement, precision and workspace contract, implementation scope, and complete-cost boundary |
| [State-based QBP T-depth](STATE_QBP_DEPTH.md) | Whole real/complex gradient schedules, live dirty workspace, and simultaneous T-count/T-depth upper bounds |
| [QBP consequence](QBP_CONSEQUENCE.md) | Exact substitution, raw-coordinate accuracy, and matched-program accounting |
| [Approximate QBP](QBP_APPROXIMATION.md) | Complete complex gradient, observable sums, rounded weights, correlated dirty reuse, and quantum/classical budgets |
| [State-only compiler](STATE_ONLY_COMPILER.md) | Two clean flags and one precision charge for a known real Hopf state, including coherent reference selection |
| [Reference-state QBP](REFERENCE_STATE_QBP.md) | Leaf-only interference with explicit reference and sampling tradeoffs |
| [Coarse-frame QBP](COARSE_FRAME_QBP.md) | All real angles, bounded depth records, a charged coarse inverse, and exact classical correction without a fine inverse frame |
| [Native coarse-frame example](NATIVE_COARSE_QBP.md) | Complete elementary two-qubit integration, an active dirty helper, complex observables, and a comparison with the original protocol |
| [Complex coarse compiler](COMPLEX_COARSE_COMPILER.md) | Gauge-fixed phase tables with exact dirty return, two-flag state preparation, and its additional dirty-bank refinement |
| [Complex coarse-frame QBP](COMPLEX_COARSE_QBP.md) | Complete magnitude and leaf-phase gradients, actual native phase-row reconstruction, and the two-stream cost ledger |
| [Native complex QBP example](NATIVE_COMPLEX_COARSE_QBP.md) | Both elementary native gradient streams, full-input dirty echoes, an exact coarse certificate, and a same-observable cost comparison |
| [Hopf QBP cost comparison](QBP_COST_COMPARISON.md) | Common accuracy and workspace, the fine gauged borrowed baseline, banked state bounds, and separate quantum, sampling, and classical costs |
| [Bounded-input QBP costs](BOUNDED_INPUT_QBP.md) | Constructive preprocessing, explicit program output, deterministic and sampled classical Pauli baselines, and the end-to-end comparison boundary |
| [Residual-table preprocessing](RESIDUAL_TABLE_PREPROCESSING.md) | Certified algebraic rotation coefficients without Euler search, unchanged error constants, and the small-system banked construction |
| [Native residual rows and tables](NATIVE_RESIDUAL_ROTATION.md) | Exact coefficient-to-mask programming, elementary residual rows, and enabled two- and four-row tables with arbitrary borrowed inputs |
| [Native residual state preparation](NATIVE_RESIDUAL_STATE.md) | One- and two-system-qubit integration, two clean flags, charged reflections, actual-inverse amplification, and complete borrowed-input isometry error |
| [Coherent residual selection](NATIVE_RESIDUAL_BRANCH.md) | One-system-qubit reference/target selection under an arbitrary branch, literal relative phase, two clean flags, and complete branch-and-dirty isometry error |
| [Native residual QBP integration](NATIVE_RESIDUAL_QBP.md) | Certified complex fixture, fine residual preparation inside both native streams, charged coarse/observable words, and exact raw-gradient histogram decoders |

### Research boundaries and related constructions

| Chapter | Role |
|---|---|
| [Research status and open endpoint](OPEN_PROBLEM.md) | Reconciled hierarchy of bounds, promised update families, and supporting components; current task-specific protocol and separate complete-frame questions |
| [Source-reuse limits](SOURCE_REUSE_LIMITS.md) | Scoped source restrictions, classical tree-generator compression, and the coherent transport and leakage obstacles |
| [Endpoint tree transport](ENDPOINT_TREE_TRANSPORT.md) | Sparse path representation, explicit normalized unitary columns, weighted norm bound, and the remaining joint precision cost |
| [Weighted transport block](WEIGHTED_TRANSPORT_BLOCK.md) | Complete one-signal-flag dilation, rejection correction, singular-case preprocessing, and native synthesis with exact or approximate dirty return |
| [Complete residual assembly](RESIDUAL_ASSEMBLY.md) | Affine forward dilation, actual reverse branch, two-clean selection and amplification, and a complete-frame bound that retains the per-depth precision charge |
| [Antichain correction compiler](ANTICHAIN_COMPILER.md) | Zero-clean $`O(N+L)`$ T-count for targets differing from the specified native baseline only on a prefix-free node set; the general frame endpoint remains open |
| [Sparse-update compiler](SPARSE_UPDATE_COMPILER.md) | One-clean $`O(N+L)`$ T-count for a sufficiently small ancestor closure, including a full changed path; the unrestricted frame endpoint remains open |

## Evidence and sources

| Page | Role |
|---|---|
| [Verification](VERIFICATION.md) | Proof-to-code map, reproducible checks, and evidence limits |
| [Fault-tolerant receipts](../verification/fault_tolerant/README.md) | Four standalone exact source/kernel/resource checks |
| [Source map](SOURCE_MAP.md) | Imported results, inherited interfaces, and local constructions |
| [Related work](RELATED_WORK.md) | Primary-source comparisons and contribution boundaries |
| [Provenance](../provenance/README.md) | Source versions and upstream lineage |

The [implementation](../compiler_robust_hopf/README.md),
[tests](../tests/README.md), [scripts](../scripts/README.md), and
[diagrams](../assets/README.md) have short
maps for reproducing the paper's constructions. The
[manuscript guide](../manuscript/README.md) gives the section order.
