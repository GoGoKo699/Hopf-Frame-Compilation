# Validation suite

[← Repository landing page](../README.md) · [State-based QBP theorem](../docs/STATE_BASED_QBP_THEOREM.md) · [Verification map](../docs/VERIFICATION.md) · [Implementation map](../compiler_robust_hopf/README.md)

The tests are organized around the proof interfaces rather than around one
particular circuit library.  They are designed to expose convention, operator,
workspace, and asymptotic-accounting errors.

Finite tests support the analytic proof; they do not establish the asymptotic
theorem by numerical extrapolation.

## Test groups

| Test file | Principal questions |
|---|---|
| [`test_frames.py`](test_frames.py) | Do the recursive and addressed constructions give the same frame? Are marker columns, chart domains, metric weights, and singular continuations consistent? |
| [`test_complex_analysis.py`](test_complex_analysis.py) | Do magnitude and leaf-phase differentials, gauge relations, and zero-amplitude behavior match the frame convention? |
| [`test_compiler_boundaries.py`](test_compiler_boundaries.py) | State-only failure, sharp all-observable sensitivity, common phases, leakage, singular markers, and checkpoint interfaces |
| [`test_strict_zero_echo.py`](test_strict_zero_echo.py) | Does the four-sector echo implement every addressed layer, restore the borrowed data bit, and compose to the complete frame? |
| [`test_strict_zero_audit.py`](test_strict_zero_audit.py) | Do the strict-zero size, depth, and lower-bound inequalities hold in exact arithmetic? |
| [`test_tree_decoder.py`](test_tree_decoder.py) | Is the binary–one-hot decoder a reversible permutation with disjoint layers, exact counts, and the intended code-space frame action? |
| [`test_router.py`](test_router.py) | Does the explicit CNOT/Fredkin router work on entangled complex inputs and return data, tokens, copies, and flags clean? |
| [`test_unified_compiler.py`](test_unified_compiler.py) | Does schedule selection cover every workspace budget and does the phase diagonal equal one exact UCG? |
| [`test_resource_bounds.py`](test_resource_bounds.py) | Do the low-workspace absorption, maximal-cut, workspace-envelope, and lower-bound diagnostics hold? |
| [`test_decoders.py`](test_decoders.py) | Do direct parity, histogram, Walsh–Hadamard, and direct phase-record decoders agree? |
| [`test_approximation_contract.py`](test_approximation_contract.py) | Actual-adjoint bias, complex phase records, reflection sums, finite weights, and conditional means with correlated dirty reuse |
| [`test_reference_state_qbp.py`](test_reference_state_qbp.py) | Leaf-interference gradients, derivative envelopes, singular angles, and the two-flag state-amplification identity; ideal small matrices |
| [`test_coarse_frame_qbp.py`](test_coarse_frame_qbp.py) | Corrected X/Y means, uniform depth-record bounds, histogram/adjoint reconstruction, and bias witnesses for omitted quadratures or coarse corrections; ideal four/eight-mode matrices |
| [`test_coarse_frame_decoder.py`](test_coarse_frame_decoder.py) | Executable histogram reconstruction, large exact counters, canceled records, and input contracts without a dense derivative table |
| [`test_native_coarse_qbp.py`](test_native_coarse_qbp.py) | Expanded elementary circuits, complete preparation isometries, coherent dirty-input score operators, literal phases, and comparison with the original protocol |
| [`test_complex_coarse_qbp.py`](test_complex_coarse_qbp.py) | Phase-gauge factorization, actual local native rows on every suffix, executable histogram reconstruction, and both gradient means; residual amplification uses small ideal matrices |
| [`test_complex_coarse_certificate.py`](test_complex_coarse_certificate.py) | Exact small-field commutator trace and rational certificate for the complex fixture's coarse error |
| [`test_native_complex_coarse_qbp.py`](test_native_complex_coarse_qbp.py) | Complete elementary prefix selections, coherent preparation, both gradient streams on arbitrary dirty input, literal gate counts, and wrong-gauge/branch-phase witnesses |
| [`test_residual_table_preprocessing.py`](test_residual_table_preprocessing.py) | Rational square-root and shared half-phase enclosures, tiny-radius/unit-boundary cases, input contracts, and finite error/workspace certificates; no quantum circuit simulation |
| [`test_rotation_programming.py`](test_rotation_programming.py) | Exact geometric-tail moments, interval certificates beyond floating-point precision, rounding/head boundaries, and invalid-input rejection |
| [`test_native_residual_rotation.py`](test_native_residual_rotation.py) | Literal emitted Ry/Rz and residual words on all small-register input columns, actual adjoints, signal symmetry, fine algebraic checks, and linear gate counts |
| [`test_native_residual_table.py`](test_native_residual_table.py) | Literal two-row masks and enable controls, exact inactive identity, active row phases on every control sector, and high-precision table certificates/counts |
| [`test_native_residual_lookup.py`](test_native_residual_lookup.py) | Four-row quadratic masks on arbitrary core inputs, temporary address tracking, literal phases, complete enabled/disabled sectors, and exact gate ledgers |
| [`test_native_residual_state.py`](test_native_residual_state.py) | Two-flag native preparation, actual inverse and reflections, complete dirty-input isometry, literal phase, and rational precision/resource certificates |
| [`test_native_two_qubit_residual_state.py`](test_native_two_qubit_residual_state.py) | Two-system-qubit preparation, exact borrowed-core reflection on arbitrary inputs, actual inverse through leakage, complete dirty-input isometry, and charged gate counts |
| [`test_native_branched_residual_state.py`](test_native_branched_residual_state.py) | Coherent residual selection on an arbitrary branch, literal relative phase, branch-independent reflection, complete branch-and-dirty isometry, and separate normalization checks |
| [`test_native_residual_qbp.py`](test_native_residual_qbp.py) | Certified complex residual fixture, complete native magnitude/phase streams on coherent dirty inputs, no-postselection score means, exact histogram decoders, and full coarse/observable gate ledger |
| [`test_operator_source_compiler.py`](test_operator_source_compiler.py) | Native two-clean frame composition, optimal source words and witnesses, dirty echoes/banks, and literal U(2) multiplexor phases |
| [`test_source_reuse_limits.py`](test_source_reuse_limits.py) | Nilpotent encoded-source dimension limits, assumption counterexamples, and transformed-mask operator identities |
| [`test_conditional_suffix_compiler.py`](test_conditional_suffix_compiler.py) | Ancestor-column residuals, separate dilation flags, conditional suffix use, complete-output amplification, and resource ledgers |
| [`test_t_depth.py`](test_t_depth.py) | Literal shared-control Fredkin batches, four disjoint T layers, and native dirty-bank queries with exact return |
| [`test_source_t_depth.py`](test_source_t_depth.py) | Exact geometric and paired-source depth certificates within the Majorana-layer architecture, paired native T layers, literal phases, and denominator witnesses |
| [`test_shallow_source_obstruction.py`](test_shallow_source_obstruction.py) | Full-input two-layer transfer alphabet, exact dyadic subset grids, optimized robust gaps, and native source witnesses with dirty extensions |
| [`test_conditional_geometric_source.py`](test_conditional_geometric_source.py) | Conditional geometric preparation, prefix cleanup, native controlled-H phases, scalar masks, actual inverse, inactive sectors, and amplification with work return |
| [`test_grouped_program_prefetch.py`](test_grouped_program_prefetch.py) | Literal masks, native AND/PREP, inactive dirty scratch, coherent program unloading through leakage, enable phases, and adaptive group reservations; reduced group reflections and queries are not a full native emitter |
| [`test_unary_phase_gradient.py`](test_unary_phase_gradient.py) | Karatsuba convolution, guarded native phase gates, exact cyclic-shift return, native small phase sources, and three unequal reduced stages with actual inverse and source error charged once per group |
| [`test_retained_source_fusion.py`](test_retained_source_fusion.py) | Exact Gaussian-dyadic shared-source polynomials for four/eight/sixteen modes, every q=8 four-mode label triple, physical q=4 bank, Bell extraction, two-round orthogonal-sector completion batching, and phase/transpose/guard counterchecks; native scope is the fixed Clifford basis |
| [`test_protected_unary_source.py`](test_protected_unary_source.py) | Native full-bank q=4 preparation, two unequal group partitions, every initial bank/logical column, actual-inverse return, retained source error including final-flag leakage, and incorrect bank-retest/inverse counterchecks; group actions are reduced operators |
| [`test_windowed_group_predicates.py`](test_windowed_group_predicates.py) | Exact rational partial and unequal windows, retained source and transient future work, arbitrary inactive cache/dirty return, native constant-arity flag phases, and late-erasure/missing-guard counterchecks; completed bodies are reduced operators |
| [`test_shared_prefix_query_audit.py`](test_shared_prefix_query_audit.py) | Native dirty residues and stale-mask failures, charged two-group fusion with signed rotations and one reused flag, actual inverses, and a 14-T affine-parity refresh; bounded feature/query interfaces, not a scalable emitter |
| [`test_grouped_selector_reuse.py`](test_grouped_selector_reuse.py) | Incremental prefix growth, early suffix-enable cleanup, exact rational three-stage leakage, inactive arbitrary-work return, literal private-copy phases, and invalid cleanup orders; reduced stages test the selector interface |
| [`test_grouped_source_reuse.py`](test_grouped_source_reuse.py) | Common-source two-stage conjugation, a wrong-reflection negative control, and legal exact-row stale-monitor/one-use leakage; native PREP with reduced masks and reflections |
| [`test_chunked_dirty_indicator.py`](test_chunked_dirty_indicator.py) | Native shared-control phases and actual inverses, all-input dirty-tree return, conjugation-order regression, and exact late-query resource sums; asymptotic counter depth is analytic |
| [`test_nonuniform_dirty_indicator.py`](test_nonuniform_dirty_indicator.py) | Unequal-chunk tree echoes on all symbolic dirty inputs, small native edge interfaces, inverse/order counterchecks, exact recursion thresholds, and normalized resource/depth majorants; no exponential-register simulation |
| [`test_hopf_error_accumulation.py`](test_hopf_error_accumulation.py) | Small ideal-frame relative spectra and coherent shared-flag leakage, using one common reduced source algebra and actual inverses; not an emitted native source word |
| [`test_flag_echo.py`](test_flag_echo.py) | Literal-source square and Pauli-echo errors, generic leakage, exact equal-mask cancellation, and actual inverses in bounded reduced source algebra |
| [`test_radial_filter.py`](test_radial_filter.py) | Fixed-point source filtering with literal phase, full polar error, phase-approximation budgets, actual inverse, and smaller precision allocation; native phase-word costs are analytic |
| [`test_parallel_dirty_lookup.py`](test_parallel_dirty_lookup.py) | Scratch-free routed indicators, actual-inverse orientation, literal phases, symbolic all-input return, complete native queries, and disjoint T layers |
| [`test_bilinear_dirty_lookup.py`](test_bilinear_dirty_lookup.py) | Rectangular bilinear basis changes, shared-target native phases, the two-indicator query echo, actual inverses, and arbitrary dirty-input return; uses existing routed indicators |
| [`test_blocked_bilinear_lookup.py`](test_blocked_bilinear_lookup.py) | Native dirty-controlled bilinear phases and selected blocks, symbolic blocked queries on every input variable, actual inverses, emitted T-layer ledgers, and missing-helper/selection counterchecks; no scalable chunked-query emitter |
| [`test_rectangular_query_allocation.py`](test_rectangular_query_allocation.py) | Exact power-of-two allocation, rectangular address partition, simultaneous pool reservation, integer resource inequalities, width/precision cap boundaries, and the square-allocation penalty; normalized analytic ledgers, no new native word |
| [`test_dirty_indicator_depth.py`](test_dirty_indicator_depth.py) | Exact emitted parity phases, two-layer dirty-helper Toffoli, complete six-wire indicator, actual inverses, and no extra helper |
| [`test_counter_dirty_indicator.py`](test_counter_dirty_indicator.py) | Two-adder signed increments, TTK and shortened RV modular-adder actions, cyclic echoes, arbitrary dirty offsets, full indicator return, and shared-address native scheduling; optimized RV depth is analytic |
| [`test_readonly_dirty_increment.py`](test_readonly_dirty_increment.py) | Two-dirty-bit involution increment, both literal polarities, full native phases and actual inverses, shared-address scheduling, and wrong-order decrement witness; optimized increment depth is analytic |
| [`test_dirty_sum_interfaces.py`](test_dirty_sum_interfaces.py) | Exact masked compressors, native phases and actual inverses, helper-offset cancellation in the full sum echo, column-height/resource ledgers, and the signed-increment cycle witness; optimized increment depth is analytic |
| [`test_pipelined_dirty_sum.py`](test_pipelined_dirty_sum.py) | Deferred parity forests, emitted doubling-block carries with settled-prefix controls, exact dirty phase localization, full masked-sum echoes, and dependency deadlines |
| [`test_batched_dirty_lookup.py`](test_batched_dirty_lookup.py) | Guarded partial indicators, symbolic multi-bank queries, literal native phase and work return, disjoint T layers, and finite workspace/depth ledgers |
| [`test_amortized_dirty_lookup.py`](test_amortized_dirty_lookup.py) | Rank-reduced controlled shears, two-pass dirty selectors, native query phases/inverses and exact return, emitted resource ledgers, capped-precision error budgets, and variable-accuracy allocation/composition checks |
| [`test_tree_residual_structure.py`](test_tree_residual_structure.py) | Complete residual reconstruction from classical tree generators, complex coarse words, singular angles, and mixed transport errors |
| [`test_tree_transport.py`](test_tree_transport.py) | Sparse transport and exact Gram identities, complete history-unitary columns, weighted subtree norms, and finite-order correction witnesses |
| [`test_weighted_transport_block.py`](test_weighted_transport_block.py) | Weighted Gram obstruction, scalar recursion versus Schur solves, complete one-flag dilations and actual inverses, level packing and physical reindexing, and zero-defect behavior |
| [`test_residual_assembly.py`](test_residual_assembly.py) | Affine subtree Schur solves, complete forward/reverse branches, two-flag selection and amplification, mode packing, all-input perturbation bounds, and wrong-phase/inverse controls; operator fixtures of dimension at most 64, not a native emitter |
| [`test_affine_tree_fusion.py`](test_affine_tree_fusion.py) | Full-input ten-mode and recursive fusion, ranks and entry/support counts, undamped continuation, native paired-source correlations, changed-address masks, and shared-signal returns; matrices of dimension at most 64 and representation-specific diagnostics |
| [`test_coupled_residual_merge.py`](test_coupled_residual_merge.py) | Complete one-signal coupled recursion, close complex native baseline with real target, height three, actual inverses, arbitrary spectators, stability, anchored mixed term, and rank-four full repair; no native repair cost claim |
| [`test_shared_conjugator_merge.py`](test_shared_conjugator_merge.py) | Literal native outer cancellation on all dirty/signal ports, two returning fork paths, nonuniform normalization, coefficient-ellipse retuning bound, and the emitted source ledger; scoped word failure |
| [`test_antichain_compiler.py`](test_antichain_compiler.py) | Exact descendant-forest factorization, mixed-depth packing, native dirty-Fredkin and reflection predicate echoes, and scoped counterchecks; matrices of dimension at most 64 |
| [`test_sparse_update_compiler.py`](test_sparse_update_compiler.py) | Ancestor-closed support and forest factorization, affine mode packing, dense column dictionaries, separate rejection flags, complete amplification, inactive sectors, and dirty-core error; small component matrices and batched columns |
| [`test_one_clean_compiler.py`](test_one_clean_compiler.py) | Native paired-Majorana source and masks, conjugated scalar blocks, five-call return bounds, borrowed-signal X symmetry, the scalar-phase counterexample, addressed phases, and exact inactive sectors |
| [`test_source_merge.py`](test_source_merge.py) | Native scalar-source parity, direct single-flag merge failure, and two-flag Pauli routing with a surviving cubic return |
| [`test_provenance.py`](test_provenance.py) | Are upstream commits, source roles, and local lineage recorded consistently? |
| [`test_literature_policy.py`](test_literature_policy.py) | Does the active proof use one compiler framework and maintain the declared contribution boundary? |
| [`test_reviewer_narrative.py`](test_reviewer_narrative.py) | Do the primary reading route, diagrams, links, terminology, and source-version statements remain coherent? |
| [`test_math_typography.py`](test_math_typography.py) | Do notation, quantifiers, adjoints, and mathematical prose retain their intended meaning? |
| [`test_hyphen_inline_math.py`](test_hyphen_inline_math.py) | Does inline mathematics next to a word hyphen use the protected GitHub syntax? |
| [`test_presentation.py`](test_presentation.py) | Does the presentation checker detect damaged math handoffs, diagram overflow, overlapping labels, and obscured connectors? |

## Highest-leverage operator checks

The following checks are especially useful when modifying the scientific code:

1. recursive frame equals addressed-layer frame;
2. state-equivalent completion changes the two-qubit decoded gradient;
3. each strict-zero echo layer equals the corresponding addressed layer;
4. route followed by inverse route is identity on arbitrary complex inputs;
5. routed tail equals the exact subtree direct sum;
6. complete routed cut equals the direct frame;
7. exact-compiler work registers return to their promised inputs;
8. one-UCG phase blocks reproduce the complete leaf-phase diagonal;
9. baseline two-clean operator-source blocks satisfy full-output error bounds, including
   core return, while lookup banks and suffix controls return exactly;
10. actual inverse circuits and sequential composition retain their bounds
    after intermediate leakage;
11. grouped conditional-suffix blocks preserve every inactive sector and
    include active suffix and predicate leakage in the complete error;
12. corrected magnitude and direct phase streams reproduce analytic gradients
    as complete dirty-input score operators, including singular tuples;
13. rational residual coefficients retain their certified error and consistent
    half-phase branch at zero, tiny radius, and the unit-circle boundary;
14. exact paired-source programs and elementary unaddressed residual rows
    preserve literal phases, borrowed-signal symmetry, and full-input bounds;
15. enabled two-row tables preserve the address/enable data and every
    inactive input exactly, with coherent relative phases retained;
16. the bounded residual state word uses both clean flags, includes all
    dirty input columns and intermediate leakage, and amplifies with its
    actual inverse and literal sign;
17. four-row masks use exact core-pivot conjugations, retain temporary
    address phases, and need no additional helper;
18. two-system-qubit preparation returns the reused reflection helper on
    arbitrary leaked inputs and includes both system bits in the initial
    reflection, with all elementary gates charged;
19. coherent residual selection excludes the arbitrary protocol branch
    from the initial reflection, retains relative phase on coherent
    inputs, and checks both branch normalizations separately;
20. the residual QBP integration uses the actual coarse inverse only in
    magnitude readout, retains all leakage in its output probabilities,
    preserves both raw phase coordinates, and charges every wrapper.

The real and complex elementary integration fixtures use two logical
qubits and exact finite-size preparation. Their bounded propagated columns
do not implement the general fine residual-table emitter. Classical
coefficient certificates and floating histogram reconstruction have their
separate scopes above. The unaddressed native row is a further implemented
component, extended to addressed two- and four-row tables and bounded
one- and two-system-qubit state amplification. The one-system-qubit
coherent selector includes an arbitrary branch and its relative phase.
Its bounded QBP integration supplies certified fixture coefficients and
both charged native readout streams with exact histogram decoders.
General lookup and the full fine state schedule remain separate.

All exact-frame, resource, QBP, one-clean primitive, and two-clean compiler
checks are retained.
Earlier exploratory endpoint constructions are preserved in repository
history; their finite checks are outside this suite. The focused
source-reuse checks above accompany the current open-problem note and
preserve its explicit interface restrictions.

## Run the suite

```bash
python validate.py
```

The validation entry point prints the Python and NumPy versions and executes the
complete unittest collection. The standalone exact receipt suites are run
separately with `python scripts/verify_fault_tolerant.py`.

The shorter orientation is:

```bash
python scripts/reviewer_walkthrough.py
```

The [verification map](../docs/VERIFICATION.md) explains which statements are
proved analytically, represented as explicit schedules, imported from the
compiler literature, or checked only in finite dimensions.
