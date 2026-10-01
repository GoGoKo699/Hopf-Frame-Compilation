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
| [`test_operator_source_compiler.py`](test_operator_source_compiler.py) | Native two-clean frame composition, optimal source words and witnesses, dirty echoes/banks, and literal U(2) multiplexor phases |
| [`test_source_reuse_limits.py`](test_source_reuse_limits.py) | Nilpotent encoded-source dimension limits, assumption counterexamples, and transformed-mask operator identities |
| [`test_conditional_suffix_compiler.py`](test_conditional_suffix_compiler.py) | Ancestor-column residuals, separate dilation flags, conditional suffix use, complete-output amplification, and resource ledgers |
| [`test_t_depth.py`](test_t_depth.py) | Literal shared-control Fredkin batches, four disjoint T layers, and native dirty-bank queries with exact return |
| [`test_parallel_dirty_lookup.py`](test_parallel_dirty_lookup.py) | Exact bilinear dirty echoes, native phases, symbolic indicator return, reusable scratch, and disjoint parallel T layers |
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
    preserve literal phases, borrowed-signal symmetry, and full-input bounds.

The real and complex elementary integration fixtures use two logical
qubits and exact finite-size preparation. Their bounded propagated columns
do not implement the general fine residual-table emitter. Classical
coefficient certificates and floating histogram reconstruction have their
separate scopes above. The unaddressed native row is a further implemented
component; general address selection and state amplification remain separate.

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
