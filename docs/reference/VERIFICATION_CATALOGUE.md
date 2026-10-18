# Verification catalogue

[Core verification map](../VERIFICATION.md) · [Source catalogue](SOURCE_CATALOGUE.md) · [Tests](../../tests/README.md)

This catalogue preserves detailed fixture contracts and software boundaries.
The core map links Results A–D and the complete-frame count theorems to their proofs and
principal evidence; this page records the smaller input, precision,
workspace, and gate-emission
scopes of the individual checks. A research fixture's presence here does
not make its construction a current principal result. Finite checks do
not prove universal statements or asymptotic optimality.

## Catalogue index

| Material | Entry point | Scope |
|---|---|---|
| C: real and complex-magnitude one-clean count frontier | [Count evidence](../VERIFICATION.md#endpoint-count-evidence) | Growing-source diagnostics, exact non-dyadic precision scales, bounded conjugate-signal recovery, and a separate phase-dressing composition check |
| D: one-clean count and depth | [Detailed finite checks](#detailed-finite-checks) | R79's analytic precision ledger and splices use unchanged R78 decoder and R72 one-clean interfaces; existing decoder, source, recovery and dirty-query checks retain their finite scopes |
| Fixed-accuracy logarithmic general depth R90 | [Detailed finite checks](#detailed-finite-checks) | Bounded dirty program, catalogue commutator, literal native factor and reference-return fixtures; fixed-number reduction, native local synthesis, catalogue absorption, selected-factor error and original-width composition remain analytic |
| Precision-terminal refinements R91–R93 | [Detailed finite checks](#detailed-finite-checks) | Seven exact reader fixtures and four grouped residual fixtures; R92's uniform envelope and R93's native resource/error and width splices remain analytic |
| Portable workspace catalogue R94 | [Detailed finite checks](#detailed-finite-checks) | Six bounded full-operator checks for transfer, initialized leakage, current selectors and guard return; the simultaneous native bounds and uniform width/precision splices are analytic |
| All-precision general depth R96 | [Detailed finite checks](#detailed-finite-checks), [all-precision proof](../../research/depth/ALL_PRECISION_LINEAR_DEPTH.md) | Four bounded capped-reference fixtures check exact overlap/cutoff, full-space completion, all original bands and occupied actual-inverse cancellation; retained native controls, source/recovery and query suites support their existing interfaces. Global addresses, common fine grammar, coarse support, precision and width splices remain analytic |
| Outer-overlap rigidity R97 | [Detailed finite checks](#detailed-finite-checks), [outer-overlap proof](../../research/depth/OUTER_OVERLAP_RIGIDITY.md) | Nine exact fixtures check native axis compatibility, transverse gauge, displayed overlap and center mixing, a distinct three-layer resplitting, and the same-space descent obstruction. The bounded-misalignment gap and its uniform aligned case remain analytic |
| Query-depth reformulation R98 | [Query/frame identities](../../research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md), [conditional transfer](../../research/depth/QUERY_DEPTH_FRAME_TRANSFER.md), [paired router](../../research/depth/PAIRED_PERMUTATION_ROUTING.md) | Exact small native query/frame words, actual-inverse phase falsifier and initialized leakage hybrid; paired routing and selected-address packing are tested on complete arbitrary inputs, including a body that temporarily changes its selector. The conditional general-frame resource/error assembly is analytic and assumes an improved query interface |
| Gaussian inverse repair | [Proof and scope](../../research/depth/GAUSSIAN_INVERSE_EXTRACTION.md), [bounded fixtures](../../tests/test_gaussian_inverse_extraction.py) | Complete small binary edge maps, Weyl phases, cheap-shear counterexample and inverse polynomial degree; retained R91 fixtures check native edge substitution. The robust affine-round obstruction is analytic, not a general T-depth bound |
| Inverse-entry representation and read-only conversion | [Proof home](../../research/depth/GAUSS_SUM_QUERY_ENCODING.md), [bounded fixtures](../../tests/test_gaussian_inverse_extraction.py) | Existing arbitrary-table encoding, scoped affine-embedding capacity and completed-generic-query pricing bounds, exact dirty-offset failure and finite-field interpolation. The fully returned read-only repair loses to direct path propagation in count, depth and work; the common-mask evaluation boundary does not exclude coupled masks, quantum interference or sparse native specialization |
| Gaussian scalar-loop obstruction | [Proof home](../../research/depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md), [exact fixtures](../../tests/test_gaussian_scalar_loop_obstruction.py) | Arf characteristic-coefficient formula, literal Fourier cube and marked path-sum inverse polar. Width-independent degree and robust return are analytic; the primitive is an actual affine quadratic phase, not an arbitrary scalar lift of a symplectic matrix |
| Two-layer read-only query rigidity | [Proof home](../../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md), [zero-clean fixtures](../../tests/test_read_only_two_layer_rigidity.py), [initialized-work fixtures](../../tests/test_read_only_clean_rigidity.py) | Literal Clifford phases, actual-inverse hybrid with leakage, Pauli separation, one-clean quadratic structure, and a native two-layer Toffoli showing why unrestricted address mixing lies outside the theorem |
| Unrestricted one-clean two-layer transfer constraints | [Proof home](../../research/depth/UNRESTRICTED_TWO_LAYER_QUERY.md), [joint-transfer fixtures](../../tests/test_joint_transfer_two_layer.py), [Galois fixtures](../../tests/test_two_layer_galois_return.py) | Exact row normalization, infinite-alphabet tails, shared native clean-axis spans and trace cancellation, full-input query coefficients, native inverses and literal rational return. Two-clean and rational non-C3 counterexamples delimit the invariant; a six-bit cubic bent query passes the complete scalar alphabet but is excluded by the subsequent joint theorem. Universal bounds and clean-return rounding are analytic |
| Unrestricted two-layer Boolean degree bound | [Proof home](../../research/depth/TWO_LAYER_DIAGONAL_RIGIDITY.md), [exact diagonal fixtures](../../tests/test_two_layer_diagonal_rigidity.py), [paid bent-query baseline](../../tests/test_bent_query_baseline.py) | Full preserved-algebra proof, paid clean removal, signed one-relation quotient analysis and native row cancellation prove degree at most two with at most one clean qubit; the nonquadratic error floor is uniform in arity and dirty width. Exact complex diagonal counterexample checks retain all occupied columns, unused sector, literal inverses and ledger; its formal real correction cannot fit the same native layers. The 40-T/90-Clifford/ten-layer bent baseline restores every input and borrowed logical helper. No dimension-growing depth or complete-frame improvement is claimed |
| Unrestricted zero-clean three-layer degree and paid one-clean query | [Proof home](../../research/depth/THREE_LAYER_NATIVE_STRUCTURE.md), [native mechanism fixtures](../../tests/test_three_layer_native_mechanisms.py) | Analytic whole-algebra proof gives XOR degree at most four and a uniform error floor with zero clean qubits. Exact one-clean AND3 uses 18 T/99 Cliffords/three layers, with all occupied columns, leakage, dirty return and full actual inverse checked; the inverse has 109 Cliffords after reduction. A smaller-count six-layer baseline and a split-specific counterexample with a one-layer simplification delimit the claims. General one-clean degree and full frame depth remain open |
| Clean-code Clifford rigidity and native probes | [Proof home](../../research/depth/CLEAN_CODE_CLIFFORD_RIGIDITY.md), [exact fixtures](../../tests/test_clean_code_clifford_rigidity.py) | Literal pointwise Clifford stabilizer, exceptional rank-two group, restriction gap and native logical localization are proved analytically. Finite native comparisons include a correlated non-Clifford extension counterexample and a growing-rank no-shared-flip family with a paid one-layer simplification. Clean pair count is not an unrestricted lower-bound invariant; no one-clean three-layer degree theorem is claimed |
| Complete quadratic XOR queries | [Proof home](../../research/depth/QUADRATIC_NATIVE_QUERY.md), [exact fixtures](../../tests/test_quadratic_native_query.py) | Two physical parity bases, dirty phase cancellation, constructive polar congruence and paid affine correction give exact optimal depth two with zero clean and one dirty for every nonaffine quadratic. Tests retain complete small matrices, general finite coordinate checks, physical bands, literal native adjoints and asymmetric inverse S counts. The improved cubic baseline is checked on all occupied columns; simultaneous shared-control batching is not inferred |
| Balanced HRT in-place recovery | [Proof home](../../research/depth/BALANCED_ALGEBRA_PERMUTATION.md), [exact native and algebra fixtures](../../tests/test_balanced_algebra_permutation.py) | Six returned bilinear calls in five native stages pay each in-place multiplier; native full-column phases and adjoints, concurrent shared controls with private T targets, shared dirty coefficient banks, balanced forward/inverse maps and simultaneous ledgers support the analytic $`O(\log r)`$ source depth with $`O(3^r)`$ counts and dirty work. All coefficient inputs are arbitrary data. The fair cardinality-batched baseline is charged; no bounded-depth counterfamily or complete-frame improvement is claimed |
| One-relation projected native matching | [Proof home](../../research/depth/ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation), [exact coefficient fixtures](../../tests/test_projected_band_matching.py) | Dependent Pauli expansions, signs, exact Fourier spectra, independent outside axes and a full native stabilizer lift verify the component. Individual axes are recovered except for a four-axis relation's three-dimensional span; its matching is necessary only. The robust condition acts on an already identified Pauli stabilizer quotient, not an arbitrary intermediate clean code |
| Exact one-clean two-layer phase strengthening | [Proof home](../../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-one-clean-queries-are-affine), [four bounded checks](../../tests/test_read_only_two_layer_exact_phase.py) | Physical one-T flat-support commutation, an abstract literal-AND falsifier outside that restriction, retained imaginary phase, and cancellation of the actual complementary-sector unitary. Exact affineness is analytic; the stronger statement has no asserted positive approximation radius |
| Three-layer read-only query rigidity | [Proof home](../../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md), [four-band checks](../../tests/test_read_only_three_layer_rigidity.py), [one-clean checks](../../tests/test_read_only_three_layer_clean_rigidity.py), [native fixtures](../../tests/test_read_only_query_layer_fixtures.py) | Thirteen bounded tests: literal actual-inverse cocycle, independent T axes and generic symplectic falsifier, nonzero/opposite code intersections, orbit extension and codimension-two derivative lemma, scalar conversion, complete zero-clean four-layer and one-clean three-layer AND words, and a two-layer relative-phase failure with every gate charged. Width-independent degree bounds are analytic and retain the preserved-address hypothesis |
| Four-layer zero-clean read-only rigidity | [Proof home](../../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md#8-four-layers-with-zero-clean-qubits-degree-at-most-two), [bounded checks](../../tests/test_read_only_four_layer_rigidity.py) | Native outer-layer normalization and auxiliary Pauli-defect identities support the analytic degree-two bound at error below $`\sqrt2/12`$. Both outer-symplectic cases are covered; arbitrary address mixing, initialized work and general depth growth remain outside this theorem |
| General five-layer zero-clean read-only rigidity | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound), [native-defect checks](../../tests/test_five_layer_native_defect.py) | Native bilinear Pauli pullback, one-band matching, moving-probe signs and translation dimensions support degree at most 14 (11 for invisible target shear) below $`\sqrt2/24`$, without extra comparison hypotheses. The analytic theorem is independent of dirty width; it is not a general-frame or arbitrary-depth lower bound |
| Six-layer native reduction and Pauli comparisons | [Proof home](../../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md), [native boundary checks](../../tests/test_six_layer_defect_boundary.py), [Pauli-comparison checks](../../tests/test_six_layer_pauli_comparison.py) | Shared three-band factorization, quadratic-label cube, literal moving-probe degree 48 and common frozen reference give degree at most 51 below $`\sqrt2/64`$ under the specified Pauli-comparison hypothesis. The conditional general detector leaves native transport unresolved |
| Actual returned odd-order comparison | [Native construction](../../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md#8-an-actual-returned-comparison-with-an-odd-order-defect), [four checks](../../tests/test_returned_six_layer_odd_group.py) | Full identity word, literal address phase, actual inverse, 4T/15G forward and 4T/13G inverse, order-three comparison/defect and persistent non-Pauli differences. Fair baseline is identity with zero gates; no saving or general depth lower bound is inferred |
| Return does not imply source coverage | [Six-slot proof](../../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md#9-full-return-supplies-no-source-wide-inverse-coverage), [three exact checks](../../tests/test_returned_six_layer_no_source_coverage.py) | A literal returned identity with six nonempty T slots has specified source THTHT, zero exact Pauli domains and no nonidentity inverse-Clifford probes. Actual comparisons are identity; this refutes a source-wide inference, not a general depth bound |
| Signed transport through a diagonal third-level source | [Proof home](../../research/depth/NATIVE_CLIFFORD_TRANSPORT.md), [six mechanism checks](../../tests/test_native_clifford_transport.py) | Literal radical matching, nonmonomial inner action, quotient observables, modulo-eight carries and actual endpoint cancellation support the analytic six-layer degree bound 78 (26 for visible shear) below $`\sqrt2/8`$. The specified prefix must be diagonal third-level in fixed Clifford coordinates; general native transport remains open |
| Native band source structure | [Proof home](../../research/depth/NATIVE_BAND_SOURCE_STRUCTURE.md), [original checks](../../tests/test_native_band_source_structure.py), [stripping checks](../../tests/test_native_commutation_cover.py) | Two-band classification, signed extension, certified stripping and the 2t exceptional-axis bound. The commutation-cover theorem closes all bad-zero-only interactions; general good-middle elimination preserves C3. Bounded nonduplicate overlap-four and non-C3 intermediate fixtures retain literal inverses. Physical axis independence is essential; return does not imply C3 membership |
| Overlap-one native source closure | [Proof home](../../research/depth/ONE_OVERLAP_SOURCE_STRUCTURE.md), [focused checks](../../tests/test_one_overlap_source_structure.py) | Paired half-angle rigidity, literal two-qubit core and its adjacent-partial failures, spectator block identity, and the dependent four-axis exception support the analytic t at most one C3 source classification and separate one-relation extension. Existing six-layer bounds 78/26 apply; no general C3 or returned-source classification is claimed |
| Overlap-two native source closure | [Proof home](../../research/depth/TWO_OVERLAP_SOURCE_STRUCTURE.md), [focused checks](../../tests/test_two_overlap_source_structure.py) | Complete activation-class analysis, coupled-core rank obstruction and full-space zero/mixed-core involutions certify all exact three-native-band C3 sources with overlap at most two. Literal phases and remaining axes are retained; the existing six-layer 78/26 bounds apply. C3 is not inferred from return |
| Native cores and coupled incidence | [Proof home](../../research/depth/ISOLATED_NATIVE_CORE_STRUCTURE.md), [focused checks](../../tests/test_isolated_native_core_structure.py) | Full involution defects and native A commutation close arbitrary coupled incidence with good zero-character B axes. The new paired-orbit certificate derives that premise when every singleton C axis is good; character nonvanishing replaces a potentially zero Gram bit. Bounded overlapping-row and zero-diagonal checks retain positive binary corrections. Mixed bad-singleton cases remain unresolved |
| Controlled-Pauli transport | [Proof home](../../research/depth/CONTROLLED_PAULI_TRANSPORT.md), [focused checks](../../tests/test_controlled_pauli_transport.py) | Signed normalization and phase carries give uniform 174/58 bounds for exact C3 scalar-Pauli sources. An intrinsic support-coset criterion, equivalently a tableau-deviation product test for C3, characterizes this source class under all fixed endpoint Cliffords. Symbolic HRT checks cover its exact domain, all permutation columns, failure at four generators and success at three. Certified domain-axis phase multiplication preserves the class. None of the attributed source fixtures is a returned three-band query |
| Hybrid diagonal/core transport | [Proof home](../../research/depth/HYBRID_DIAGONAL_CORE_TRANSPORT.md), [focused checks](../../tests/test_hybrid_diagonal_core_transport.py) | Signed Pauli-domain localization, ordered controlled-Pauli phases and coordinatewise four-block recovery with the mod-four sign carry give a width-independent degree bound depending on core size. Relation-count localization preserves the wide diagonal phase. These are explicit C3 source hypotheses, not general native depth invariants |
| Uniform bounded-deficiency transport | [Proof home](../../research/depth/BOUNDED_DEFICIENCY_TRANSPORT.md), [focused checks](../../tests/test_bounded_deficiency_transport.py) | Degenerate-domain localization, literal Pauli phase reconstruction and approximate near-Clifford endpoints support a source-uniform radius. The intrinsic codimension parameter need not be bounded by T-depth |
| Permutation-source transport | [Scoped proof](../../research/depth/PERMUTATION_SOURCE_TRANSPORT.md), [focused checks](../../tests/test_permutation_source_transport.py) | Uniform monomial endpoint gap, forward-degree and phase bounds, and a fully native four-band non-C3 candidate with actual inverse. A qualifying three-band returned instance remains unconstructed; failure of one direct band merge is not a synthesis lower bound |
| Third-level-only transport counterexample | [Scope and native ledger](../../research/depth/NATIVE_CLIFFORD_TRANSPORT.md#7-third-level-membership-alone-does-not-bound-tableau-degree), [scope checks](../../tests/test_c3_transport_scope.py) | The known He–Robitaille–Tan permutation family has growing symplectic and literal sign degree under exact Clifford endpoint transport. Every source and inverse is paid; no three-band synthesis or returned-query counterexample is supplied |
| Five-layer partial-normalizer reduction | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md), [four bounded checks](../../tests/test_read_only_five_layer_reduction.py) | Six-band normalization with the actual inverse and literal cocycle phases; a full-return native affine query with non-Clifford middle representation; translation-span and Hermitian-Pauli sign checks. The degree-seven bound (four for a visible shear) requires the additional Clifford-middle hypothesis |
| Five-layer Pauli-comparison subcase | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#51-a-pauli-valued-comparison-robust-quadratic-bound), [five bounded checks](../../tests/test_pauli_comparison_cubic.py) | Native-axis cocycle kernel, quadratic baseline, actual-inverse cube product, peeled projective separation and cubic-contraction rank support the analytic degree-two bound below $`\sqrt2/12`$. The comparison must be Pauli-valued; the middle may be non-Clifford |
| Five-layer outer-axis-normalizer subcase | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#52-an-outer-axis-normalizer-comparison-a-finite-degree-bound), [bounded checks](../../tests/test_outer_axis_normalizer.py) | Varying-tableau phase carries, signed native-axis correction, translation-space dimensions and a returned query satisfying neither earlier extra hypothesis support the analytic degree-244 bound (41 for an invisible target). Individual-axis permutation remains an additional condition; a mixed-axis CNOT fixture supports the separate third-level cocycle reduction under span preservation |
| Five-layer span-preserving subcase | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#54-span-preserving-comparisons-a-finite-degree-bound), [bounded checks](../../tests/test_span_preserving_augmentation.py) | Native shear nilpotence, subgroup augmentation filtrations and the strict translation-dimension threshold support degree at most 112 below $`\sqrt2/24`$; the shared-prefix combination gives 21 below $`\sqrt2/40`$. The analytic proof imports the stated Jennings basis theorem; finite tests are not an arbitrary-width proof. K must preserve the outer axis span; the degree bounds are superseded by the general theorem |
| Shared-prefix five-layer comparison invariant | [Proof home](../../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#55-shared-prefix-comparisons-a-robust-quadratic-group-map), [bounded checks](../../tests/test_shared_prefix_comparison.py) | Central native-axis commutators, bilinear address defect, short-word robustness and faithful degree-four augmentation features. No span hypothesis is needed below $`\sqrt2/40`$; target extraction is subsequently resolved by the native-defect theorem |
| Baseline-word obstruction | [Proof home](../../research/depth/BASELINE_WORD_OBSTRUCTION.md), [five checks](../../tests/test_baseline_word_obstruction.py) | Native full-column source and inverse, unbounded odd symplectic orders, incidence Gram matrix and literal target cancellation. The exact and robust word obstruction concerns independent baselines. The added [returned-affine checks](../../tests/test_returned_affine_baseline.py) verify an actual order-two comparison group with odd-order outer baseline, literal return and complete native ledger; this affine example is not a high-degree returned query |
| Scalar representation products | [Proof home](../../research/depth/SCALAR_REPRESENTATION_PRODUCT.md), [bounded checks](../../tests/test_scalar_representation_product.py) | Literal twisted-Pauli translations, exact four-factor scalar return, source/inverse costs, commutation/displacement ranks and flat Clifford Pauli support. The added source-count barrier survives Pauli-carrier compression and small full-source error; it does not lower-bound fused whole-word synthesis or general query depth |
| Native one-band matching and partial-action scope | [Proof home](../../research/depth/ONE_BAND_AXIS_MATCHING.md), [band checks](../../tests/test_one_band_axis_matching.py), [scope checks](../../tests/test_partial_normalizer_scope.py) | Independent-axis expansion, fixed-error separation, signed/inactive matching and non-Clifford intermediate cancellation; the controlled-bank example tests full-tableau degree and a paid native source. No general two-band partial-action bound is asserted |
| Exact phase relocation | [Proof home](../../research/depth/READ_ONLY_PHASE_RELOCATION.md), [six exact checks](../../tests/test_read_only_phase_relocation.py) | Null phase identity, actual four-layer CNOT bases and inverses, affine constants, complete dirty return and simultaneous ledgers; the hybrid scalar CCZ attains four layers with 64 T gates, 218 Cliffords and 17 dirty helpers, including its shared-work inverse order. No asymptotic query improvement is claimed |
| Returned diagonal source words and translation echo | [Proof home](../../research/depth/RETURNED_DIAGONAL_ECHO_DEGREE.md), [four source-word checks](../../tests/test_diagonal_source_word.py), [eight echo checks](../../tests/test_returned_diagonal_echo_degree.py) | Joint address/work differences annihilate q arbitrary diagonal sources under affine monomial-Clifford transport, giving degree at most $`\max(2,q-1)`$ below $`2\sin(\pi/2^{\max(q,3)+1})`$. The echo is cubic below $`2\sin(\pi/32)`$; all source errors and actual inverses are charged. Source occurrences are not native T-depth |
| Nonlinear returned unit-group characters | [Proof home](../../research/depth/NONLINEAR_CHARACTER_TRANSPORT.md), [native checks](../../tests/test_nonlinear_character_transport.py) | All-character degree bound, nonlinear multiplication, literal dirty return, in-place factorization and charged controlled phases. The extended proof covers every commutative nilpotent binary coefficient algebra. A degree-a output requires at least $`2^{a-1}`$ coefficient wires in this architecture; the complete slow native construction does not improve R91 or lower-bound arbitrary queries |
| Noncommutative cycle character | [Proof home](../../research/depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md), [native checks](../../tests/test_noncommutative_cycle_character.py) | Degree-d return on d-squared arbitrary coefficient wires, full dirty/output columns, literal factor source and inverse, shared-control phase schedule and exact sequential ledgers. The added unitriangular-inverse identity, complete word-query basis change and dirty controlled-source wrapper retain every output column and charge both source calls. The added commuting-pivot inverse is exact on all matrix inputs with paid private helpers, cubic counts and linear depth; whole-pivot commutation does not parallelize individual gates or establish native read-only access. The block-pivot/XOR extension retains all helpers and exposes the inverse counter-action in a proposed balanced recurrence. Its determinant obstruction is restricted to frozen-coefficient uniform array shears. A later two-shear computation retains dirty garbage until actual readout/uncompute, giving a complete structured inverse-word query and phase source with logarithmic T-depth, cubic counts/work, direct charged controls and exact return. Bounded native columns and recursive arbitrary-input checks support the analytic recurrence. Coefficient inputs change internally; no native read-only arbitrary-table query follows. The full conjunction transport remains slower than the retained baseline; the earlier inverse-degree obstruction is only for affine/disjoint-Toffoli circuits |
| Multiway query cost audit | [Proof home](../../research/depth/MULTIWAY_QUERY_COST_BARRIER.md) | Analytic tensor-description and coherent-branch entropy counts; exact hash-bucket and cyclic-slot return failures. The T-count consequence requires termwise emission and is not a lower bound for joint native synthesis |
| Logarithmic-precision general depth R95 | [Detailed finite checks](#detailed-finite-checks) | Four bounded sparse-correction checks for live root and signed Schur sources, complete-work recovery, actual inverses and leakage; retained native/source, R91 and R94 fixtures support the inherited interfaces, while the asymptotic proof is analytic |
| Retained causal compact-coin theorem R89 | [Detailed finite checks](#detailed-finite-checks) | Four causal-batching fixtures, seven compact-coin fixtures, and two subset-routing fixtures; catalogue existence, native resources, prefix recursion, full-isometry error, and all-width composition remain analytic; its coarse reductions support R90 |
| Retained product-coin theorem R87 | [Detailed finite checks](#detailed-finite-checks) | Seven bounded product-coin fixtures; its fine decoder supports R89, and native synthesis, full dirty/reference contract, rounding, and all-width resources remain analytic |
| Retained terminal R88 and contiguous routing | [Detailed finite checks](#detailed-finite-checks) | Three selected-terminal matrix fixtures and three logical interval-routing fixtures; R88's sharp-query and fallback interfaces support R89's stronger terminal theorem |
| Interval-dilation predecessor R86 | [Detailed finite checks](#detailed-finite-checks) | Seven bounded interval-dilation fixtures; the normalization, native implementation, full dirty/reference contract, and resource bounds are analytic |
| Separate state-based QBP result | [Coverage and software boundary](#state-based-qbp-coverage) | State/reference preparation with a changed decoder |
| Endpoint source and representation fixtures | [Proof/module/test map](#endpoint-source-and-representation-checks) | Joint completions, ordered banks, recovery, boundary structure, encodings, reflections, and source resynthesis |
| Detailed construction fixtures | [Detailed finite checks](#detailed-finite-checks) | Core ingredients, supporting schedules, and research diagnostics, with individual limitations |
| Complete exact-frame checks | [Core verification §§2–8](../VERIFICATION.md#2-what-is-represented-locally) | Frame identities, explicit routing, resource ledgers, and records |
| Four standalone exact receipts | [Receipt index](../../verification/fault_tolerant/README.md) | Gray source, shift kernel, resources, and reflection |
| Reproduction and rendering | [Commands](../VERIFICATION.md#9-reproduce-the-checks), [rendering details](#rendering-checks) | Deterministic checks and presentation coverage |

The detailed paragraphs are retained as an evidence record, including
negative controls and unsuccessful routes. Their actual return guarantees
and finite ranges take precedence over any shorthand description.

### State-based QBP coverage

The [task theorem](../../supplements/state_based_qbp/STATE_BASED_QBP_THEOREM.md) and the bounded native
residual-to-gradient integration are complete within their stated scopes.
The following map separates each analytic claim from its executable
coverage. The detailed finite checks below retain their individual
input, precision, and workspace boundaries.

| Claim | Proof home | Executable coverage | Software boundary |
|---|---|---|---|
| Common coarse reference with literal phase and exact dirty return | [Complex coarse compiler §§1–6](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md) | [Real](../../supplements/state_based_qbp/NATIVE_COARSE_QBP.md), [complex](../../supplements/state_based_qbp/NATIVE_COMPLEX_COARSE_QBP.md), and [residual QBP](../../supplements/state_based_qbp/NATIVE_RESIDUAL_QBP.md) fixtures record actual coarse words | No integrated coarse synthesizer for all admitted inputs |
| Fine single-state and coherent-pair isometries, including returned work | [State compiler §§1–7](../../supplements/state_based_qbp/STATE_ONLY_COMPILER.md) and [banked extension](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) | Certified rows, two- and four-row tables, one- and two-qubit preparations, and a coherent one-qubit selector | No variable-size or banked native state emitter |
| Original raw gradient means, confidence, and reconstruction accuracy | [Real decoder](../../supplements/state_based_qbp/COARSE_FRAME_QBP.md) and [complex decoder §§2–6](../../supplements/state_based_qbp/COMPLEX_COARSE_QBP.md) | General histogram utilities test the identities; [the residual fixture](../../supplements/state_based_qbp/NATIVE_RESIDUAL_QBP.md) has exact rational decoders | General coefficient contractions still use floating point; the guarded decoder is proved, not implemented |
| Compiler T-count and simultaneous T-depth upper schedules | [Task ledger](../../supplements/state_based_qbp/STATE_BASED_QBP_THEOREM.md#4-quantum-and-classical-resource-ledger) and [depth proof](../../supplements/state_based_qbp/STATE_QBP_DEPTH.md) | Exact resource ledgers, native dirty-lookup checks, and bounded literal gate inventories | No full native program realizing the variable-size banked depth schedule |
| Polynomial construction and explicit program-output bound under the input restriction | [Bounded-input audit §§3–4](../../supplements/state_based_qbp/BOUNDED_INPUT_QBP.md#3-which-native-searches-are-polynomial-in-the-input-parameters) | [Residual intervals](../../supplements/state_based_qbp/RESIDUAL_TABLE_PREPROCESSING.md), exact sign programming, and bounded emitted words | No integrated certified front end from the admitted bounded Hopf inputs to the full program |
| Fine residual preparation inside both complete gradient streams | [Bounded residual QBP integration](../../supplements/state_based_qbp/NATIVE_RESIDUAL_QBP.md) | [Native words and rational decoders](../../compiler_robust_hopf/native_residual_qbp.py), checked by [independent oracles](../../tests/test_native_residual_qbp.py) | Fixed complex one-qubit target; no same-accuracy native fine-frame benchmark or sampled-accuracy experiment |

The polynomial construction statement applies at the basic reservation
for $`n\ge6`$, or at the banked reservation for every $`n\ge1`$.
The basic-budget $`n\le5`$ fallback retains its sufficient
$`2^{O(P)}`$ search bound. No polynomial evaluation bound is inferred
for arbitrary computable inputs. The bounded one- and two-qubit native
preparations also retain their extra-four and extra-three dirty-wire
allocations; they do not implement those minimum-budget fallbacks.

The remaining software can be grouped into three optional deliverables:

1. A certified front end for the admitted bounded Hopf inputs, including
   the recorded coarse word and certified residual data.
2. A variable-size native emitter with general tables, predicates,
   reflections, bank lifetimes, and the claimed count/depth schedule.
3. General guarded decoder arithmetic with the stated output-error
   certificate; exact integer histograms alone do not certify the later
   coefficient contractions.

The corresponding construction and error arguments are supplied by the
linked proofs. These software gaps limit executable coverage; none is an
unfinished prerequisite for the theorem as stated. The fixed rational
fixture discharges its own coefficient promises and does not certify an
arbitrary caller's data.

The bounded implementation pass stops here. No larger fixed fixture or
general software API is selected. Further implementation should begin
with a concrete Hopf-QBP use requirement and an input/output contract
identifying the missing interface. The
[Result C proof](../../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
and its [banked frontier](../../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
establish the one-clean count guarantees at the [stated reservations](../../manuscript/PUBLICATION_SCOPE.md);
the general T-depth frontier remains an
[active research question](../OPEN_PROBLEM.md#full-circuit-requirement).
Count and depth match in explicit accuracy/workspace ranges by the
[precision-uniform subtree theorem](../../research/depth/UNIFORM_SUBTREE_DEPTH.md),
R79, for both frame families. R72 retains its earlier one-clean statement;
the earlier real and complex theorems retain their two-clean proof scopes.
The existing
state-based result does not require the remaining questions' resolution.
Application-level advantage is outside the current research scope; the
existing resource comparisons and classical baselines remain documented.

### Endpoint source and representation checks

These fixtures support the endpoint research chapters. Exact arithmetic,
emitted native component words, and bounded floating matrix checks are
distinguished below; none emits the unrestricted fine-precision complete
frame. Floating rank, support, spectral, and chart decisions are diagnostics,
not certified decisions about arbitrary computable inputs. The linked proofs
supply the all-size statements within their declared interfaces.

| Construction | Proof home | Module and tests | Evidence and boundary |
|---|---|---|---|
| Joint source completion | [Original-tree residual §5](../../research/endpoint/ORIGINAL_TREE_RESIDUAL.md#5-shared-factors-contraction-slack-and-the-completion-gap) | [source_completions](../../compiler_robust_hopf/source_completions.py); [ordered-bank tests](../../tests/test_ordered_block_banks.py) | Regularized 4N-mode isometries, full-unitary extensions, cut/Gram bounds and rejected-space identities are finite matrix checks; the root-series obstruction is rational-exact. No fine source emitter. |
| Ordered banks and native readout | [Ordered block banks](../../research/endpoint/ORDERED_BLOCK_BANKS.md) | [ordered_block_banks](../../compiler_robust_hopf/ordered_block_banks.py), [native_ordered_banks](../../compiler_robust_hopf/native_ordered_banks.py); [tests](../../tests/test_ordered_block_banks.py) | Exact integer schedules and rational witnesses; bounded source matrices; native right-readout and coarse-prefix words with small full dirty interpreters and occupied/entangled-vector checks. No fine left-source or endpoint emitter. |
| Conjugate-sector recovery | [Conjugate signal recovery](../../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md) | [conjugate_recovery](../../compiler_robust_hopf/conjugate_recovery.py); [tests](../../tests/test_conjugate_recovery.py) | Exact native recovery counterexample and gate lifts; bounded full-sector recovery and phase-error matrices; 16-dimensional initialized-real-target recovery excludes the borrowed signal from its projector and retains full-operator leakage, actual inverse, and dirty/reference inputs. Its source is supplied; no fine source is emitted. |
| Boundary rank and forest deflation | [Rank-weighted compiler](../../research/endpoint/RANK_WEIGHTED_BOUNDARY_COMPILER.md), [matching/leaf deflation](../../research/endpoint/MATCHING_LEAF_DEFLATION.md) | [boundary_refinement](../../compiler_robust_hopf/boundary_refinement.py); [tests](../../tests/test_boundary_refinement.py) | Finite residual-rank, reflection-LCU/recovery, boundary cancellation, transport and closed-forest checks; exact combinatorial matching/bin checks. Floating ranks do not supply a native compiler or endpoint error certificate. |
| Paid logical encodings | [Paid encodings](../../research/endpoint/PAID_LOGICAL_ENCODINGS.md) | [paid_logical_encodings](../../compiler_robust_hopf/paid_logical_encodings.py); [tests](../../tests/test_paid_logical_encodings.py) | Complete Boolean permutations, exact six-bit native matrix, encoded controls and bounded frame/code diagnostics. Encoders act on all inputs; zero remainders require the stated code. No arbitrary-angle core is emitted. |
| Center-phase reflection family | [Centered phase reflections](../../research/endpoint/CENTERED_PHASE_REFLECTIONS.md) | [centered_phase_reflections](../../compiler_robust_hopf/centered_phase_reflections.py); [tests](../../tests/test_centered_phase_reflections.py) | One five-qubit native member: exact initialized columns and full actual inverse, with 144 T appearances. Floating checks test optimized floors and real-family identities; they do not certify general endpoint synthesis. |
| Programmed reflection orbits | [Programmed reflections §§1–4](../../research/endpoint/PROGRAMMED_REFLECTIONS.md) | [programmed_reflections](../../compiler_robust_hopf/programmed_reflections.py); [tests](../../tests/test_programmed_reflections.py) | Exact rational dual certificate and fixed native corner/inverse checks; bounded orbit, kernel and moment diagnostics. A kernel solution alone is not full programming coverage; no variable-angle program or eigenbasis is emitted. |
| Regular spectral chart | [Programmed reflections §§5–7](../../research/endpoint/PROGRAMMED_REFLECTIONS.md#5-a-regular-signed-cycle-chart) | [spectral_reflection_chart](../../compiler_robust_hopf/spectral_reflection_chart.py); [tests](../../tests/test_spectral_reflection_chart.py) | Fixed native signed-cycle and alternate-baseline words have exact dirty-return/inverse checks. Spectral, polar and quadrature calculations are bounded floating diagnostics, not emitted variable spectral reflections. |
| Raw-source resynthesis | [Source-program resynthesis §§1–4](../../research/endpoint/SOURCE_PROGRAM_RESYNTHESIS.md) | [source_resynthesis](../../compiler_robust_hopf/source_resynthesis.py); [tests](../../tests/test_source_resynthesis.py) | Exact dyadic masks, rational source ranks/denominator growth and binary Lie closures; bounded common-basis, alignment and product matrices. Arbitrary rotations in those matrices are not supplied native gates. |
| Native selected tails | [Source-program resynthesis §4](../../research/endpoint/SOURCE_PROGRAM_RESYNTHESIS.md#4-conserved-prefix-raw-programs-have-a-coherently-priced-basis) | [native_source_resynthesis](../../compiler_robust_hopf/native_source_resynthesis.py); [tests](../../tests/test_native_source_resynthesis.py) | Exact four-qubit selected-tail word and literal inverse; complete dirty Boolean queries and bounded full-matrix/reference checks. Emits selected-tail components, not arbitrary-angle fine resynthesis. |
| Borrowed selected axes | [Borrowed-axis compiler](../../research/endpoint/BORROWED_AXIS_COMPILER.md) | [exact tests](../../tests/test_borrowed_pauli_rotations.py) | Full CQ2 matrices with unequal selected angles, arbitrary borrowed/core columns, signed Pauli phases and actual reverse inverses. Checks the conjugation identity; no arbitrary-angle native emitter or endpoint bound is inferred from the fixture. |
| Paired-source rejection orientations | [Orientation bank](../../research/endpoint/BORROWED_AXIS_COMPILER.md) | [tests](../../tests/test_paired_source_orientation.py) | Exact small raw-source words and Clifford-vector identities, bounded normalized full-bank gauges, and exact sign-class dimensions. The joint O(Q+q+ell) compiler is analytic; these fixtures do not emit its fine rotations. |
| Completed-group resynthesis | [Source-program resynthesis §§5–6](../../research/endpoint/SOURCE_PROGRAM_RESYNTHESIS.md#5-completed-groups-admit-independent-gauges-without-connectors) | [source_group_resynthesis](../../compiler_robust_hopf/source_group_resynthesis.py); [tests](../../tests/test_source_group_resynthesis.py) | Exact rational coefficient/mask tables and raw-connector sign witness; bounded complete recovery/gauge words in faithful compact Clifford coordinates, retaining both central sectors and leakage. No large physical-core simulation or gate emitter. |

### Detailed finite checks

The [homogeneous-magnitude theorem R80](../../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md)
has four [operator fixtures](../../tests/test_homogeneous_depth.py):
independent layer and marker-column frames for complex local unitaries;
literal native conjugator/involution echoes on all dirty and tag inputs;
complete tag return and inactive occupied work; and rational sign gauges
on heterogeneous signed/singular trees. The bounded query action is a
supplied exact XOR fixture. Native query synthesis, recursive workspace,
all-width bounds and fixed-accuracy asymptotics remain analytic.
The [terminal theorem R81](../../research/depth/RECURSIVE_TERMINAL_DEPTH.md)
uses retained native identities with a prefix-uniform resource proof.
The uniform terminal-height extension to all $`0\le M\le n`$ is
analytic; it adds no primitive or cutoff-certifying fixture. Existing
immutable receipts are unchanged.

The [affine orientation R83](../../research/depth/AFFINE_ORIENTATION_DEPTH.md)
and [deterministic frame R84](../../research/depth/DETERMINISTIC_FRAME_DEPTH.md)
proofs have five [orientation fixtures](../../tests/test_orientation_depth.py):
affine marker correction with outside-prefix dependence; the complete
shared-tag circuit with actual complex native factors; a literal
Clifford+T guarded affine map; signed deterministic source and
last-deviation cleanup; and the failed simultaneous dirty update.
Expected frames use independent path products. Occupied inactive tags,
all small dirty columns, negative corners, and marker phases are retained.
Queries use their complete reversible truth-table action; the native
query schedule, all-width recurrences, and asymptotic bounds remain
analytic. R82's affine-reference extension reuses the proved reference
core and adds no supplied frame oracle or new immutable receipt.

The [geometric-band theorem R85](../../research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md)
has six [orientation-recovery fixtures](../../tests/test_geometric_band_orientation.py):
integer geometric capacities and bin counts; all eligible truncated-frame
columns and local last-one bands on a four-qubit nonlinear tree; the
physical save/query/swap/inverse-query/restore word, its actual inverse,
external dirty return even on occupied tails, and inactive identity;
coherent fine-tag clear/reload across bins with zero pairs, node signs,
leaf phases, and a dirty system entangled with a reference; complete
coarse search/XOR/reverse on every column of a five-qubit frame with five
zero-based bins and a
padded-threshold countercheck; and occupied inactive protected records.
Expected frames use independent path products. Rotations are bounded
ideal matrices, and queries are exact basis permutations. The coarse
fixture reuses R78's literal Hadamard/reflection/complement search, with
projectors reconstructed from descendant scaling columns. Native
synthesis, full source-error and dirty-reference guarantees, finite
cutoffs, and resource bounds follow from the analytic proof; no scalable
emitter or scientific receipt is added.

The [interval-dilation theorem R86](../../research/depth/INTERVAL_DILATION_DEPTH.md)
has seven [bounded fixtures](../../tests/test_interval_dilation_depth.py):

1. Integer atom intervals and the isometric embedding preserve every
   prescribed subtree scaling and marker column of a small weighted tree.
2. The rectangular overlap reconstructs the full frame, keeps the root
   scaling column and literal signs, and has the expected per-scale
   sparsity; its norm check is a finite diagnostic.
3. Membership XOR erases the original leaf label, and the actual inverse
   has the required accepted block on every atom input, independently
   of the chosen source completion.
4. The square anchor has the required accepted frame block while
   retaining old-label, failure, and source-work rejection coordinates.
5. Calibrated amplification preserves literal phase; a perturbed source
   and its actual inverse retain final leakage within the finite hybrid
   bound.
6. Fine-tag erasure covers every local band, with an occupied-work
   perturbation and coherent reference inputs.
7. The guarded balanced-frame word permits variable active widths and
   is exactly identity on occupied inactive tag and padding registers.

These use bounded ideal matrices, explicit basis permutations, and
independent tree-column formulas. They do not emit the variable-size
native compiler or prove the universal Haar norm, arbitrary dirty/reference
contract, global rounding, resource bounds, or all-width splice. Those
are analytic in R86 and its
[selected-routing support](../../research/depth/SELECTED_DIRTY_ROUTING.md).
No new immutable scientific receipt is added. The
[session record](../../research/depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md#5-verification-and-integration)
records the separate full-validation and presentation checks.

The [product-coin theorem R87](../../research/depth/COIN_DILATION_DEPTH.md)
has seven [bounded fixtures](../../tests/test_coin_dilation_depth.py):

1. Heterogeneous coin rectangles partition the full coin space, form an
   isometric embedding, and preserve every prescribed subtree scaling flag.
2. Every marker has its literal mean-zero split coin and uniform future
   coin blocks, with unit norm.
3. Tensor Hadamards place every complete-frame column in exactly its
   last-nonzero-coin-block sector, including the zero root label.
4. The physical source and membership XOR realize the embedding, and
   the actual inverse returns the initial column for distinct source
   completions.
5. The full word clears the fine tag on every column with literal phase at zero
   and nonzero offsets, including coherent external reference inputs.
6. A perturbed source retains nonzero final leakage within the complete
   initialized-output norm bound, including reference correlations.
7. On the inactive sector the whole word cancels exactly on occupied
   source, coin, and tag work with an external reference.

These checks use bounded ideal matrices and explicit basis permutations.
They test the product-coin identities and cleanup interfaces; they do
not emit the native family source or prove arbitrary-size resource,
rounding, dirty/reference, or all-width guarantees. Those follow from
the analytic R87 proof and its retained premises. No new immutable
scientific receipt is added. The
[session record](../../research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
records review and integration separately.

The [prefix-selected terminal theorem R88](../../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md)
has three [bounded fixtures](../../tests/test_prefix_selected_coin_terminal.py):

1. Independently constructed complete selected-frame columns agree with
   ordered layer factors for signed and singular tuples; the early
   factor is identity on occupied protected sectors, and nonzero
   external prefixes select genuinely different frames.
2. Successive stages may use different positive dyadic denominators on
   their own layers. Their product preserves the complete rounded-frame
   target, and its error obeys the sum of the individual stage bounds,
   including prescribed descendants below zero incoming amplitudes.
3. The full operator error is the maximum over the preserved prefix
   sectors, also after coherent prefix mixing and tensoring a reference.

The [contiguous-layout router](../../research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
has three [logical permutation fixtures](../../tests/test_interval_barrel_routing.py):

1. Variable intervals, including empty intervals, exchange exactly the
   selected bank prefix and return arbitrary padding and dirty controls.
2. Inactive routing is identity for arbitrary parameter words, including
   lengths and starts that would be invalid on the active sector.
3. Two disjoint interval exchanges preserve parked external-bank values;
   a nontrivial active body followed by the actual inverse layout returns
   every bank, padding wire, and control.

The terminal fixtures use small real matrices and the routing fixtures
use explicit logical permutations. They do not emit the native family
sources, sharp queries, or Fredkin T-layer schedule. The recursive
error and resource sums, arbitrary dirty/reference guarantees, and
all-width allocation are analytic. R88's terminal theorem and the
contiguous router initially left R87's complete-frame order unchanged;
R89 subsequently improves the complete and terminal bounds. R79's
absolute all-precision scope remains unchanged. No new immutable
scientific receipt is added.

The [causal compact-coin theorem R89](../../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md)
has four [bounded fixtures](../../tests/test_causal_compact_coin_depth.py):

1. Batched compression clears all High fields with prefix-dependent
   metadata and genuine staged-hash collisions; the actual inverse
   returns every logical child and metadata field.
2. The global short-shift word implements the causal local matchings
   injectively and clears its selector using the inverse selector.
   The fixture's chosen complete local permutations also permit bounded
   off-image checks; no general off-image matching contract is inferred.
3. Every complete-frame column lies in its prescribed last-nonzero-block
   band. The complete coarse-tag XOR and actual inverses retain literal
   phases, arbitrary old tag values, and an external reference.
4. An explicit reference-source completion places an orthogonal marker
   in source work while its coin remains zero, showing why the block
   predicate must test all initialized reference-source work as well.

These fixtures use explicit finite matching catalogues, logical reversible
words, and small ideal matrices. They do not establish polynomial
catalogue existence or synthesize the variable-size native sources and
queries. The simultaneous resource envelope, prefix recursion, summable
errors, arbitrary dirty/reference guarantee, and literal-width splice
are proved in R89.

The [dirty catalogue theorem R90](../../research/depth/DIRTY_CATALOGUE_TERMINAL.md)
has [bounded fixtures](../../tests/test_dirty_catalogue_terminal.py):

1. A complete selected commutator uses different noncommuting H/S/T
   words at two catalogue slots, retains their actual inverse words, and
   checks every logical, prefix, program, and array input with literal phase.
2. A coherent prefix and an explicit dirty/reference Bell pair retain
   their full joint return contract.
3. Reusing one array leaves the predicted nontrivial dirty commutator;
   two disjoint arrays cancel that residue.
4. An independently constructed four-column real Hopf frame is recovered
   by explicitly computed unitary commutator factors. This spectral
   calculation is a bounded numerical fixture, not a native compiler.
5. A nonlinear two-bit program tests the dirty XOR translation echo,
   actual inverse, and selected swap on every basis input.
6. A shared-control CCZ phase-polynomial identity checks the literal
   common-control contribution used by the native Fredkin batching lemma.

The fixtures use small matrices, explicit basis permutations, and literal
native factors. They do not emit arbitrary catalogue factors or establish
the asymptotic allocation. The fixed-number R89 reduction, zero-clean
local synthesis, selected-factor error bound, program/catalogue/query
resource absorption, full-isometry composition, and literal-width splice
are proved in R90. No efficient classical preprocessing, practical
crossover, depth optimality, scalable native emitter, or new immutable
scientific receipt is asserted. The
[session outcome](../../research/depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
records review and integration separately.

The [R91 parallel word query](../../research/depth/PARALLEL_BILINEAR_WORD_QUERY.md)
has [seven bounded fixtures](../../tests/test_parallel_bilinear_word_query.py).
Integer phases modulo eight and explicit basis permutations test the
literal private-target bilinear echo, actual inverses, controlled
arbitrary-dirty action, and four-corner cancellation. A separate native
H/CNOT/T matrix checks the shared-control Toffoli substitution. The
allocation checks cover balanced, blocked, wide-output and full-indicator
regimes. A native sparse shear and a composed masked-program selected
swap check the paired router, distinguishing its temporarily permuted
auxiliary array from exact return in the full conjugation. These fixtures
do not prove the asymptotic allocation or query optimality.

The [R93 grouped terminal](../../research/depth/GROUPED_CATALOGUE_PRECISION.md)
has [four bounded fixtures](../../tests/test_grouped_catalogue_precision.py).
The actual nine-wire Clifford+T matrix-unit dilation is checked on all
512 columns, with arbitrary labels, rejection flag and helper, exact
inactive identity, literal phase and actual inverse. An explicit
64-type, four-digit SELECT contracts to the stated accepted kernel;
independent source and atom rejection flags, common normalization, and
inactive G/H behavior are checked separately. A two-group noncommuting
example changes the later prefix. Negative controls freeze that address,
merge failure flags, or project the shared source too early. These checks
stop at the bounded operator boundaries; retained shared-source tests
cover history, padding and recovery. The finite fixture parameters do
not assert the theorem's precision schedule or emit a general circuit.

[R92](../../research/depth/UNIFORM_CATALOGUE_PRECISION.md) changes the
uniform analytic ledger for existing circuit words and needs no new
operator fixture. Its absolute cutoffs and precision window, and R93's
local decoder, full error, simultaneous widths, geometric resources and
real/complex splices are proved analytically and independently reviewed.
The [precision session record](../../research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
records validation and integration. No immutable receipt is regenerated.

The [R94 portable-workspace catalogue](../../research/depth/PORTABLE_WORKSPACE_CATALOGUE.md)
has [six bounded fixtures](../../tests/test_portable_workspace_catalogue.py).
Complete native H/S/T/CNOT words check full-payload extraction, coherent
selection between noncommuting entries, literal phase, and actual inverse.
A complete controlled Hopf matrix checks every marker column and the
forward branch; an explicit leaking extension checks the initialized
inverse estimate and twice-error commutator bound. A partial-routing
counterexample leaves a dirty-slot residue. A changing-prefix example
rejects a frozen selector. The guard fixture distinguishes exact
restoration of its original value from arbitrary-input tensor independence,
and detects failure when its predicate reads changing payload bits.
All operators have at most 256 rows and columns.

The last interfaces use small abstract operators where stated; they do
not emit the full R79 local compiler or a general frame circuit. R91's
retained fixtures cover dirty-program echoes and the complete paired
router. Exact array return and initialized logical approximation are
distinct contracts in the proof. Uniform catalogue absorption, geometric
query sums, the full error hybrid, constant bottom, and original-width
real/complex splices are analytic and independently reviewed. The
[bootstrap session record](../../research/depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md)
records validation and integration; no immutable receipt is regenerated.

The [R96 all-precision theorem](../../research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
has [four bounded capped-reference fixtures](../../tests/test_fast_capped_geometric_reference.py)
for its [fast geometric reference component](../../research/depth/FAST_CAPPED_GEOMETRIC_REFERENCE.md):

1. Original cap sizes from 2 through 128 and every smaller power-of-two
   cap check the literal overlap, squared-distance formula, exact equal-cap
   case, and precision cutoff. Renormalizing a truncated prefix instead of
   using the exact capped column gives a different vector.
2. One augmented reference block checks the analysis-only plane rotation,
   its full-operator distance, the transported initialized anchor, and its
   identity on occupied private scratch. Uniform-column factors and equal
   cap dimensions are included.
3. Two original eight-level coins test all 64 root/marker columns and all
   four incoming label values. The ideal comparison uses every private
   scratch bit and every original high coin bit in its augmented nonzero
   test. The actual shortened references obey the full-output sandwich
   bound. Omitting either set of coordinates fails. The original positive
   tail weights also distinguish the unchanged mathematical matching from
   the shortened native preparation.
4. Coherent active/inactive inputs include occupied private work, a dirty
   bit and an external reference bit. The actual inverse gives exact
   inactive cancellation and the active full-output error bound. Replacing
   it with the comparison inverse breaks that cancellation.

These are dense finite verification matrices with deliberately varied
unspecified completion columns. They are not emitted reference-source
oracles. The plane correction is analysis-only; its existence is not
priced as an available native gate. The exact recursive Gray-path
source and the full-space error argument have their analytic proof homes.

The complete R96 proof also uses the
[global stage](../../research/depth/GLOBAL_BORROWED_PRECISION_STAGE.md),
[polynomial-precision bridge](../../research/depth/POLYNOMIAL_PRECISION_DEPTH.md),
and [direct-prefix terminal](../../research/depth/DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md).
The retained [native grouped-terminal fixtures](../../tests/test_grouped_terminal_depth.py)
check literal H/T/T-dagger control words, phase corrections, actual
unloading, and occupied inactive programs and helpers. They support
bounded native-control interfaces, not the new unbounded program schedule.
Existing [growing-source diagnostics](../../tests/test_growing_precision_source.py),
[conjugate recovery](../../tests/test_conjugate_recovery.py),
[one-clean native words](../../tests/test_one_clean_compiler.py),
[R91 queries](../../tests/test_parallel_bilinear_word_query.py), and
[R94 complete-payload transfers](../../tests/test_portable_workspace_catalogue.py)
retain their separate finite contracts.

The common fine skeleton across immutable prefixes, exact completed-mask
substitution and actual gauge replay, the precision-independent coarse
support, every global table factor, simultaneous dirty pools, absolute
cutoffs and all-precision initialized-isometry composition are analytic.
The full law applies to every $`n\ge1,L\ge6`$ with one external clean
qubit and the original real/complex widths; finite tests do not certify
its asymptotic bounds or unrestricted depth optimality. No scientific
receipt, new external premise or general native emitter is added.
Selected manuscript Result D remains R79.

The [R95 logarithmic-precision theorem](../../research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
has [four bounded fixtures](../../tests/test_coarse_flag_sparse_correction.py)
for its [coarse-frame sparse correction](../../research/depth/COARSE_FLAG_SPARSE_CORRECTION.md):

1. Frames of heights one through four include singular branches and negative
   coefficients. Their complete signed columns, permitted residual support,
   live root/state-column residual, and absolute residual norm bound are
   checked against an independent frame construction.
2. Single-leaf signed Schur row and column sources at heights one through
   three have the claimed overlap, including zero residual columns and the
   root tag. Replacing the root by pure failure or dropping the coefficient
   signs gives explicit failures.
3. A height-two normalization-two recovery retains all source and branch
   work, the literal leading minus, actual inverse, and full initialized-work
   reflection anchor. Omitting the sign or a source tag from the anchor fails.
   Occupied-work inactive identity and coherent dirty/reference inputs are
   checked separately.
4. A perturbed height-two source completion retains nonzero final-work
   leakage within the full-output error bound and remains an isometry.
   Substituting the stale ideal inverse for the actual inverse changes the
   output.

These are small dense matrix audits. Their source completions are verification
objects, not emitted native source oracles. The retained
[one-clean native words](../../tests/test_one_clean_compiler.py),
[operator-source checks](../../tests/test_operator_source_compiler.py),
[conjugate recovery](../../tests/test_conjugate_recovery.py),
[shared binary sources](../../tests/test_shared_binary_terminal.py),
[R91 exact queries](../../tests/test_parallel_bilinear_word_query.py), and
[R94 complete-payload transfers](../../tests/test_portable_workspace_catalogue.py)
retain their individual scopes. The
[small-workspace source](../../research/depth/SMALL_WORKSPACE_DESCENDANT_SOURCE.md),
[two-stage local compiler](../../research/depth/RELAXED_COUNT_LOCAL_DEPTH.md),
native schedules, full-isometry error, catalogue absorption, and uniform
width/precision resource bounds are proved analytically. No new scientific
receipt, external primitive, general native emitter, or experimentally
certified asymptotic cutoff is added. Selected manuscript Result D remains
R79.

The [compact fixed-coin support](../../research/depth/COMPACT_COIN_FLAGS.md)
has [seven bounded checks](../../tests/test_compact_coin_flags.py):
universal finite atom partitions; exact mass and acyclicity of two-outcome
columns; a heavy-fiber obstruction to a short label-preserving rank source;
an equal-weight swap that fixes the root but changes a coarse flag;
every coarse scaling vector; whole coarse flags for different fixed-coin
completions; and full tag return with reference inputs. The two negative
examples delimit interfaces, not unrestricted native circuits.

The [short-label subset support](../../research/depth/SHORT_LABEL_SUBSET_ROUTING.md)
has [two logical fixtures](../../tests/test_short_label_subset_routing.py):
the dense-set shift word with its full guarded extension and actual
inverse, and staged-shear compression with collisions and return on
reached codes. These chosen matchings and codes are not prescribed
arbitrary permutations or standalone frame erasers. Catalogue existence
and native query costs remain analytic. No new immutable scientific
receipt or general native emitter is added.

The [permutation checkpoint](../../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
has two [boundary fixtures](../../tests/test_permutation_bridge_limits.py):
the exact four-shear one-clean order-eight word on all allowed inputs,
with a triangular order-sixteen contrast, and explicit linear Clifford
cycles with approximate spectral embeddings. The latter checks both the
isometry and its operator-norm intertwining error. These bounded examples
guard the scope of the analytic restrictions; the dense encoder is not a
supplied native primitive. No general depth claim or receipt is added.

The [bounded-density theorem R82](../../research/depth/BOUNDED_DENSITY_DEPTH.md)
has six [finite fixtures](../../tests/test_weighted_haar_depth.py):

1. Independently normalized weighted-frame columns agree with the retained
   Haar gauge, including uniform weights and zero correction columns.
2. Exact rational subtree variance equals the sum of martingale energies,
   and the corresponding finite Carleson bound holds.
3. Row and column sources share one leaf/depth label, have orthonormal
   columns and the prescribed Schur overlap, and retain separate failure
   padding, including the full state-input failure column.
4. A biased positive product reference checks the relative-density gauge,
   exact rational weighted variance and Carleson identities, and both
   signed correction-source overlaps. It tests reference geometry, not
   native synthesis of that reference frame.
5. A complete small three-branch LCU retains source, branch, product-factor
   and damping flags through literal amplitude amplification. Its active
   initialized columns return all work; arbitrary occupied work and
   complex amplitudes on the inactive branch see literal identity.
6. A perturbed source is paired with its actual inverse. The complete
   initialized-isometry error obeys the finite hybrid bound, while
   nonzero final leakage is retained and column normalization is checked.

The variance tests use exact rational arithmetic; frame, source-overlap,
and recovery comparisons are bounded floating-point checks. Source
completions are supplied verification matrices, not native source
oracles. These fixtures establish neither the general arbitrary-dirty/
reference contract nor asymptotic count, T-depth, workspace, precision,
all-width, or bias-independent constants. Those follow from the chapter's
native synthesis and analytic ledger. No large native emitter, efficient
reference search, or new immutable scientific receipt is asserted. The
supporting-family result does not itself improve the unrestricted
frontier; R86, R87, R89, and R90 supply the later fixed-accuracy general
improvements.

The [precision-uniform subtree theorem R79](../../research/depth/UNIFORM_SUBTREE_DEPTH.md)
reuses R78's decoder and R72's one-clean terminal and fallback words.
The upgrade retains the source's explicit precision factors, chooses an
enlarged protected bank, and proves an absolute cutoff and all-precision
composition analytically. It adds no circuit identity or finite fixture.
The five decoder checks below and the retained source, terminal and query
suites keep their existing contracts; they do not certify an unbounded
precision range or the asymptotic absorption inequalities. No receipt or
native-emitter guarantee changes.

The [subtree-projector decoder R78](../../research/depth/SUBTREE_PROJECTOR_DECODER.md)
retains its fixed-accuracy polylogarithmic theorem with one external clean qubit.
Its five [focused tests](../../tests/test_subtree_projector_decoder.py)
check exact rational signed/singular marker geometry, nested projector
bands, reversible label extraction, padded thresholds, literal phases,
occupied inactive work, guarded tag XOR, and leaked-source full-projector
bounds with arbitrary unitary completions. Wrong signs, padding, guards,
and inverse substitutions are rejected by counterchecks. Native source
realization, query costs, polynomial workspace absorption, all-width
composition and asymptotic depth remain analytic; the bounded completions
are verification fixtures, not source emitters. Existing exact receipts
are unchanged.

The [one-clean count/depth theorem](../../research/depth/ONE_CLEAN_DEPTH.md),
R72, adds an analytic native schedule, protected terminal flag relocation,
a width/precision splice, and one-clean literal phase composition. The
[one-clean suite](../../tests/test_one_clean_compiler.py) checks the native
phase and real words, complete input columns, and actual inverses. Retained
source, recovery and dirty-query suites support their component contracts.
They do not certify the relocated logical bank, unbounded schedule, or
asymptotic cutoffs by testing. R72 changes no native fixture or immutable
receipt and supplies no general native emitter.

The retained two-clean [all-precision depth theorem](../../research/depth/ALL_PRECISION_DEPTH.md)
adds two analytic schedule proofs. The
[whole-frame stream](../../research/depth/WHOLE_FRAME_BINARY_DEPTH.md)
reuses the existing shared-binary-terminal native identities with an explicit
growing source cardinality. The
[high-precision schedule](../../research/depth/HIGH_PRECISION_DEPTH.md)
uses exact equality of completed real masks after conjugate lifting,
preserving their full rejected operators and actual replays. Existing
shared-source, conjugate-recovery and blocked-query checks support those
component contracts. They do not test the new asymptotic cutoffs or establish
the all-precision theorem by extrapolation. No receipt is changed and no
general native emitter is claimed.

The [phase-dressed extension](../../research/depth/COMPLEX_FRAME_DEPTH.md)
adds an analytic depth schedule for the packed literal diagonal, including
the routed-indicator fallback at the stated workspace constant, and a
half-error full-isometry composition. The retained
[operator-source suite](../../tests/test_operator_source_compiler.py)
checks literal diagonal phases, actual inverses, arbitrary core inputs,
and shared-flag composition without discarding intermediate leakage.
The [nonuniform indicator](../../tests/test_nonuniform_dirty_indicator.py)
and [blocked-query](../../tests/test_blocked_bilinear_lookup.py) checks
support the exact query substitutions. These component checks do not emit
the complete complex circuit or establish its asymptotic bounds by testing.
The existing native fixtures and immutable receipts remain unchanged.

The [packed rejection-bank proof](../../research/endpoint/PACKED_REJECTION_BANKS.md)
has [integer-routing and full-unitary diagnostics](../../compiler_robust_hopf/packed_rejection_banks.py)
and [tests](../../tests/test_packed_rejection_banks.py). They distinguish exact
mode permutations from bounded floating checks of dirty rejected blocks,
preweight, recovery leakage, and incorrect-inverse/bank-collision controls.
No fine-precision native word is emitted.

The [balanced-history tests](../../tests/test_balanced_rejection_history.py)
check exact rational attenuation and preweight identities, complete integer
bank routes, and a bounded complex canonical cascade with prefix unwinding
and exact Hadamard compression. Omitting the prefix bank breaks scalar
compression while preserving the accepted block and ordinary final recovery.
These tests supply the prefix operator explicitly; they do not price it or
identify rounded raw-source defects with ideal balanced defects.

The [depthwise balance tests](../../tests/test_depth_balanced_rejection.py)
check the revised equal-amplitude schedule within each three-axis depth:
exact rational weights, the geometric comparison bound and literal width,
complete routed Walsh gates, and a bounded complex cascade with separate
success-only coarse operations. The proof supplies the native and asymptotic
cost statements. The finite matrices still supply the fine prefix explicitly;
they are not a prefix compiler or an endpoint gate-count certificate.

The [joint-history tests](../../tests/test_joint_history_compression.py)
check the [three-call extraction reduction](../../research/endpoint/JOINT_HISTORY_COMPRESSION.md):
complete active and inactive columns, arbitrary occupied-input completions,
dirty entanglement, the 24-times error bound, and a conserved-label
constant-gap example. The compressor is supplied as a finite matrix;
these checks do not construct a native compressor or improve its count.

The [history-reflection tests](../../tests/test_history_reflections.py)
check the anchored Householder completion and its bounded normalization,
the full-operator reflection error from approximate compressor columns,
and a promised root family whose range reflection suppresses orientation
changes. Exact scalar checks support the restricted query-bound constants;
the native T-count remains a separate question.

The [joint-word proof](../../research/endpoint/JOINT_NATIVE_WORDS.md) has
[exact phase, projected-Clifford and error-Gram witnesses](../../tests/test_joint_native_words.py).
The source-leakage fixture explicitly avoids assuming that the approximate
input lies in the projected Clifford's source code. The all-width compression
and counting statements follow from the proof, not these finite examples.

The [modular-clock source audit](SOURCE_CATALOGUE.md#modular-clock-counterexample)
has [exact rational full-circuit checks](../../tests/test_modular_history_audit.py)
of two-pulse return, repeated accumulation and a nonwrapping comparison.
The supplied-unitary counterexample is distinct from a native emitter.

The [operator-source proof](../OPERATOR_SOURCE_COMPILER.md) is accompanied by
[finite circuit checks](../../tests/test_operator_source_compiler.py).
These test the native geometric operator, its programmable scalar block,
the two-flag rotation, coherent amplification, dirty word-bank routing, the
restored suffix echo, and literal phase-diagonal composition on arbitrary dirty
inputs.
The analytic proof establishes the improved upper bound; finite matrices test
its fragile identities and conventions.

The suite now includes a complete two-layer **native Clifford+T** fixture
with nine wires and exactly two initialized flags. All 128 logical/dirty
input columns are propagated through the actual gate words, including an
exact seven-T Toffoli decomposition, the literal amplification phase,
actual reversed-word inverses, and retained intermediate leakage. Its
deliberately coarse precision tests composition rather than claiming the
theorem's asymptotic accuracy. Separate fixtures verify general U(2)
multiplexor Euler order, address-dependent phases, and four-stage composition,
as well as the optimal exact source words and Pauli-transfer witnesses.

The [packed diagonal tests](../../tests/test_packed_diagonal.py) reuse the
retained exact arithmetic and native gate evaluator. They check the one-tail
paired-Majorana source on one to three dirty core qubits, every small sign
mask, literal scalar and complex blocks, and the full-output amplification
identity. These checks accompany the reduced dirty reservation in
[operator-source Section 8](../OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames).
The bounded-factor consequence and terminating classical search have their
analytic proof in [the factorization reduction](../../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md).

The [tree phase/chart tests](../../tests/test_tree_phase_chart.py) use exact
rational-complex matrices for the all-angle phase gauge, polynomial
resolvent, singular inputs, coarse-permutation conjugacy, stable cotangent
chart, subtree flow sums, and Haar diagonal-difference identity. They
support the [fixed-tree reduction](../../research/endpoint/TREE_CAYLEY_REDUCTION.md),
not a native implementation of its remaining Cayley operator.

The [boundary propagation tests](../../tests/test_tree_boundary.py) compare
the complete sparse-preconditioned core with an independently ordered
edge-phase product, including unequal signed and zero parameters. They
check triangular coefficients, nearest-parent recurrence, root-chain
compression, the all-table singular gap, dummy states and all dilation
ports. The complete feedback entrance, chain reversal and recurrence are
checked on every logical column, including the cancellation that leaves
one prefix encoder. The local rotation factorization retains literal
phases, and the repeated-call reservation and error fraction are exact
integer/rational checks. The [resource proof](../../research/endpoint/BOUNDARY_PROPAGATION.md)
imports the charged borrowed-signal primitive; the dense reference
matrices do not emit that precision-dependent native circuit.

The [self-borrowed query tests](../../tests/test_self_borrowed_query.py)
propagate all 64 input columns in 24 native XOR/Pauli fixtures, including
the actual inverse. Exact Boolean checks cover overlapping output and
selector bits, the full-core Hadamard sign identity, and a source-core bit
used as the central MCX helper. These check the all-input contracts in the
[source compiler](../OPERATOR_SOURCE_COMPILER.md) and
[one-clean compiler](../ONE_CLEAN_COMPILER.md); output borrowing supplies no
initialized register.

The [coarse prefix tests](../../tests/test_coarse_prefix_encoder.py) check
the complete depth-frame/Wold encoder, delimiter permutation, and
full-operator small-conjugation bound, together with exact resource
ledgers. The [encoder theorem](../../research/endpoint/COARSE_PREFIX_ENCODER.md)
establishes the arbitrary-size native construction at coarse precision.

The [collective refinement tests](../../tests/test_collective_precision.py)
check the exact routed geometric history on occupied as well as initialized
inputs, all-column tangent and midpoint comparisons, generic two-flag
purification, opposite-sign cancellation with dirty work, and the same
circuit's precision/width ledger. Their
[analytic proof](../../research/endpoint/COLLECTIVE_PRECISION_REFINEMENT.md)
gives the quadratic and cubic native replacements. Reduced operators test
the correction identities; elementary source and borrowed-control synthesis
are imported from their separately checked contracts. No full variable-size
refinement emitter is asserted.

The [endpoint interface tests](../../tests/test_endpoint_refinement_routes.py)
check dense three-point, four-mask XOR and unequal gauge witnesses,
native dirty Klein dressing and unequal-address higher-grade terms,
and exact correlated history identities. The precision-six inverse-pair
fixture uses the literal full five-call source words, compares the dirty
polar operator, and verifies its scalar-defect two-flag filter. Lookup
fixtures test uniform-program invariance and changing-address transport.
These support the structural and source-interface restrictions, without
establishing source-witness reachability by the fixed global encoder.

The [coupled-tree tests](../../tests/test_coupled_tree_resolvent.py)
check graded powers, tails and computational input support, independent
descendant parameters in a rank-one column, and a complete complex
30-mode scattering unitary. Its zero-defect rejected column agrees with
the prefix-history formula. The
[resolvent proof](../../research/endpoint/COUPLED_TREE_RESOLVENT.md)
establishes the all-size bounds; these ideal matrix checks do not emit
the native precision-dependent completion.

The [disjoint-block tests](../../tests/test_disjoint_block_compiler.py)
check the literal dirty clock phase on every input, adjacent-row SU meshes
including padding, complex pair-source orientation, and full accepted/rejected
amplification with inactive labels. The [block proof](../../research/endpoint/DISJOINT_BLOCK_COMPILER.md)
supplies the source coverage and native resource theorem.

The [precision-gauge tests](../../tests/test_precision_gauge_compiler.py)
compare complete complex regular and lossy gauges, exact Parseval networks,
actual arbitrary-dirty rejection replay, source overlap, and full recovery.
The [precision-trace theorem](../../research/endpoint/PRECISION_GAUGE_COMPILER.md)
establishes its dimension-independent native schedule. The
[heterogeneous-source checks](../../tests/test_heterogeneous_precision_source.py)
test the finite geometric residual, exact exponent remainders and dyadic
coefficient grids, saturated shifts, and exponent-only precision ledgers.
Connected noncommuting planes share one source with retained failure
history; occupied feature columns, zero-merge splitters, and pair-symmetric
metadata return have independent checks and negative controls. Mixed-chart
fixtures check the full row word with noncommuting dirty rejection,
complex splitter phases and known-zero nodes, then actual source replay,
$`W/2`$ overlap and complete initialized-isometry recovery. Negative
controls replace dirty rejection by its norm, reverse replay, drop a
splitter phase or incorrectly include dirty inputs in the clean reflection.
Exact integer metadata also has an overflow regression. All-height checks
cover rational certified clipping on every frame column, dyadic scalar-order
padding without changing the logical dimension, the selected retained
bottom-depth trace, and source shifts beyond the clipped cap. The
[all-height common-source proof](../../research/endpoint/HETEROGENEOUS_PRECISION_SOURCE.md)
prices one prepared source, the complete row network with actual cheap
dirty-source action, global fine replay, and final recovery at
$`O(N+\mathcal B_K)`$ cost. It charges the actual selected bottom trace
and imposes no exponent promise on the upper tree. The former exponent
envelopes are corollaries. These finite identities are not a variable-size
native emitter.

The separate [growing-source helper](../../compiler_robust_hopf/growing_precision_source.py)
records exact rational refinement factors and symbolic source scales.
Its [tests](../../tests/test_growing_precision_source.py) check the normalized
rare-event refinement and its dirty-dependent unitary dilation, fresh-cell
permutations that fix occupied and entangled sealed failures, and a bounded
128-coordinate two-stage source word. The latter checks frozen source tags,
height-dependent attenuation compensation, the actual noncommuting dirty
rejection replay, and the complete $`W/2`$ overlap. Negative controls omit
the geometric normalization, reuse sealed failures, reflect about final
source tags or selected dirty values, and replace the actual inverse.
A separate perturbation check includes the rejected output and verifies
linear full-isometry recovery error from an accepted-block error. Exact
integer checks preserve $`q_dM_d=4N`$ across four binary refinements at
each scale change. An optional precision argument sets K to the least power
of two at least L. Exact checks at non-dyadic precision boundaries, including
$`L\gt N`$, preserve $`q_dM_d=4K`$ and leave the logical mode
schedule unchanged. Refinement accepts equal rounded source precision and
rejects different rounded precision. A nonintegral q remains a symbolic
schedule, not a claim that the growing-source construction applies there.
These are finite diagnostics, not emitted
Clifford+T source words or a certificate of asymptotic arena capacity.
The [growing-source proof](../../research/endpoint/GROWING_SOURCE_COMPILER.md)
supplies the full-space arena allocation, exact guards, actual dirty replay,
summable native charges and literal dirty reservation. Its
[one-clean recovery](../../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
and [banked composition](../../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
establish Result C; the [claim ledger](../../manuscript/PUBLICATION_SCOPE.md)
retains its exact domains and the distinct two-clean endpoint.
All-size summability, fallback dispatch and workspace feasibility come
from these proofs, not the schedule tests.

The same suite's separate composition fixture uses one logical bit, two
arbitrary dirty-core bits, and two shared clean flags: a 32-dimensional
full unitary with all eight logical/dirty input columns checked. An explicit
diagnostic real recovery retains nonzero dirty-dependent leakage; the
following packed-diagonal stage uses actual q=2 amplified native words in
a supplied addressed direct sum. Unequal rational unit-circle phases and
a common phase $`i`$ check literal composition, the full-output error sum,
the actual inverse, and a full-Schmidt reference input. Counterchecks discard
intermediate leakage, reverse the factor order, or omit the common phase.
The real factor is not emitted growing-source synthesis, and the addressed
direct sum is not native table emission. This bounded check does not realize
the requested general precision or certify its all-size workspace bound. The fixtures
and internal proof review do not constitute independent external review or
a scalable native emitter.

The [conjugate-recovery suite](../../tests/test_conjugate_recovery.py)
separately checks a 16-dimensional word on an arbitrary borrowed signal,
an initialized mode suffix, one logical bit, and one dirty bit. Its
supplied source has a real accepted target and a complex complementary
action. A native coarse interlayer and its actual lifted inverse test
the whole-word conjugate construction. The recovery projector includes
the original mode suffix and excludes the borrowed signal; both signal
sectors and all eight signal/logical/dirty input columns are checked.
A perturbation mixes signal and suffix sectors while acting on the dirty
bit; the test retains nonzero leakage and checks the full-operator
three-call hybrid, the actual inverse, and an eight-label reference with
full Schmidt support. Negative controls add the signal to the projector
or leave the coarse interlayer unlifted. This fixture supplies neither a
fine source emitter nor evidence of asymptotic optimality by simulation.

The [structured-forest tests](../../tests/test_structured_forests.py) independently
check computational degree counts and growth, stable star transfer, and exact
component dressing for signed and corner angles. A reversed-product negative
control detects an incorrect quotient-depth order. The
[exceptional-exit tests](../../tests/test_exceptional_exits.py) check stable
star transfer followed by dressed-cut degree sorting and a second band
factorization on every column. They also verify structural row envelopes,
commuting distinct-boundary cuts, repeated quotient child degrees, exact
dyadic packing, long exceptional chains, and finite block-support budgets.
The helper records combinatorial certificates without inferring support from
numerical zeros or emitting gates. The native conclusion
uses the [forest proof](../../research/endpoint/STRUCTURED_FOREST_COMPILERS.md),
not a gate count inferred from these matrices.

The [canonical-star tests](../../tests/test_canonical_star_resolvent.py)
check both Cayley congruences, the triangular descendant denominator, full
complex stars, midpoint completions, support packing and literal-width bounds.
The [oracle theorem](../../research/endpoint/CANONICAL_STAR_RESOLVENT.md)
prices local scattering, while inverse recovery remains separate.

The [Schur tests](../../tests/test_scalar_schur_limits.py) verify the full
complex-phase gauge, exact rational-complex completion determinants with
unconstrained upper entries, support and minor-weight bounds, and a bounded-rank
approximation example. The [completion proof](../../research/endpoint/SCALAR_SCHUR_LIMITS.md)
establishes rank and robust error floors for the stated scalar interface
and its bounded, logical-label-preserving dirty-operator extension. The
existing finite tests cover the scalar identities; the dimension-independent
extension uses the analytic block-rank and inverse estimates.

The [operator-stability tests](../../tests/test_frame_operator_stability.py)
use an independent product-rule derivative to check the tangent PSD sum,
sharp factor-two path, uniform-depth spectral equality at arbitrary base
angles and signs, finite perturbations, and isolated-angle precision.
The [interface proof](../HOPF_INTERFACE.md#42-complete-frame-operator-stability)
complements the existing [depth-cap theorem](../HOPF_ERROR_ACCUMULATION.md).
Neither ideal-angle estimate replaces a native source-leakage bound.

The [native representation tests](../../tests/test_native_representation_limits.py)
verify realification gate words in exact cyclotomic arithmetic, the
rectangular norm and full-isometry identity, bounded integer circle
solutions, and literal marker-column gaps. Exact character and phase-grid
enumeration checks Spin multiplicities at n=3 and n=4, retaining the
small-dimension exception. The
[realification](../../research/endpoint/NATIVE_REALIFICATION.md) and
[Spin-decoding](../../research/endpoint/SPIN_DECODING_LIMITS.md) proofs
give their arbitrary-size conclusions under the stated support and
wrapper hypotheses. Additional exact integer-wavelet and four-square
Householder fixtures check the dyadic construction, including determinant
and unused-quadrant identities. The linear-height and seven-factor
bounds are analytic; ring height is not a native T-count certificate.
Exact native two-qubit words additionally check the normalized accepted
scalar phase bound, both Galois estimates, and tiny acceptance amplitudes.
A rephased-boundary counterexample checks why the
[phase theorem](../../research/endpoint/NATIVE_REALIFICATION.md#4-normalized-accepted-amplitudes-retain-a-phase-cost)
requires literal Clifford calibration. It concerns the specified coherent
scalar amplitude, not a success probability or arbitrary retained features.
These representation checks do not establish a stronger unrestricted
endpoint bound.

The [exact Haar tests](../../tests/test_global_haar.py) use arithmetic over
$`\mathbb Q(\sqrt2)`$ to verify common signed-path witnesses and
absolute-matrix certificates. The arbitrary-order Haar obstruction, its
fixed-alphabet extension and the stronger alternating bound are
[analytic](../../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md).
The spectral encoding restriction, branch-rank bound and return rigidity
also have analytic proofs; none is inferred from a finite fitting campaign.

The [source-reuse limits](../../research/endpoint/SOURCE_REUSE_LIMITS.md) have
[separate finite checks](../../tests/test_source_reuse_limits.py) for nilpotent
contractions, arbitrary encoding bases, dirty-dimension independence, and
transformed-mask operator identities, exact source hoisting,
matched-Majorana commutation and rank-limited reflection extraction.
Single-involution echo checks retain clean leakage and test failures when
lookup unitarity or target commutation is omitted; the analytic bound
also requires the lookup to preserve the initialized subspace.
Counterexamples test why nilpotence and complete-output accuracy cannot
be omitted. These support the scoped
analytic restrictions; they do not establish a full-frame impossibility.

The [structural endpoint checks](../../tests/test_endpoint_structural_limits.py)
test the [full-frame metric](../HOPF_INTERFACE.md#41-complete-frame-parameter-geometry),
[normalization and residual witnesses](../../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md),
singular-safe nested projectors, restricted reflection and displacement
bounds, native control phases, query parity, and source-return counterexamples.
The centered-stage fixtures additionally check pairwise anticommuting target
Paulis, clean/dirty projection, the exact n-slot realization, and a scalar
finite-error witness under the separate Lipschitz hypothesis.
They use small complete matrices and independent formulas, including
unequal and singular angles. The proofs establish the all-dimension
claims; these numerical checks are not a general compiler, a proof of the
fixed-mixer theorem, or an unrestricted endpoint impossibility result.

The [representation-interface tests](../../tests/test_endpoint_representation_limits.py)
add exact integer block-polynomial witnesses and binary tree-label rank,
plus finite floating checks of sharp entry mass, full-dirty Spin covariance,
overlap and reflection distance. The Walsh baseline classification is
analytic; fitting a finite collection of masks is not used as its evidence.

The [conditional-suffix proof](../CONDITIONAL_SUFFIX_COMPILER.md) has
[focused grouped-block checks](../../tests/test_conditional_suffix_compiler.py).
They combine a native small operator source with separate scalar and
atom flags, retain every dirty input column through amplification,
and check literal phases, the disjoint ancestor-column support partition,
forward and reverse column blocks, inactive suffix sectors, and predicate
uncomputation after leakage. A separate integer ledger tests exponentially
growing groups and the dirty-width/error inequalities. Its sample
grouping constants do not certify the unspecified fixed constants of the
coarse synthesis primitive. These are finite interface checks, not a native
elementary circuit emitter for the whole asymptotic grouped compiler.

The [T-depth schedule](../T_DEPTH_COMPILER.md) has
[native routing checks](../../tests/test_t_depth.py) for literal four-layer
shared-control Fredkins, disjoint T-layer supports, actual inverses, and
whole dirty-bank queries. Clifford layers may have substantial depth;
the fixture does not treat T-depth as total execution depth.
The [exact source-depth proof](../../research/depth/SOURCE_T_DEPTH.md) has
[separate checks](../../tests/test_source_t_depth.py) for its source schedules,
literal native phases, and exact denominator witnesses. The matching lower
bounds apply only to the defined Majorana-layer architecture and exact
targets, including its arbitrary returned dirty extensions. They are not
general Hopf or controlled-source depth lower bounds.
The [unrestricted two-layer obstruction](../../research/depth/SHALLOW_SOURCE_OBSTRUCTION.md)
has [bounded checks](../../tests/test_shallow_source_obstruction.py) of the
full transfer alphabet for native words with unrestricted Clifford
interlayers, exact rational subset grids, optimized witness gaps, and
native source/controlled-source witnesses with arbitrary dirty extensions.
Its robust constant lower bounds are analytic, not inferred from sampling
the finite circuits. Initialized-clean isometries remain outside its scope.
The separate [clean-aware proof](../../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md)
has [eleven bounded checks](../../tests/test_clean_aware_depth_obstruction.py)
for the initialized-projector sum and robustness, exact complete-frame
witnesses, the explicit one-clean gap with an infinite-tail certificate,
mod-eight residue events, a native three-layer counterexample to the
two-layer alphabet, the sharp preserved-Pauli gap, the literal zero-clean
Clifford attainer, and the one-clean terminal-sector identity and its
occupied-sector limitation. The original fixed-clean three-layer proof
uses the cited polynomial-limit theorem; R74's initialized-sector
argument now supplies an independent elementary explicit threshold.
These tests do not prove the imported theorem. The separate
[four-layer proof](../../research/depth/FOUR_LAYER_CROSSING_RANK.md) has
[bounded checks](../../tests/test_four_layer_crossing_rank.py) of
the symplectic normal form with intersecting spans, native boundary words
attaining the normalization bound, and the abstract cubic-bias family.
The elementary initialized-observable proof gives an explicit
$`2^{-O_a(r^2)}`$ gap and minimum crossing rank
$`\Omega_a(\sqrt L)`$. Its phase-digit and two-sided projector identities
have finite checks in the same suite; the width-independent gap and
rank bound are analytic. No native realization of the abstract scalar
family is inferred.
The [four-layer operator proof](../../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md)
has [five finite audits](../../tests/test_four_layer_operator_rigidity.py):
the exact one-qubit image orbit, a literal four-layer word and dirty
extension, arbitrary Hermitian-Clifford boundaries, common active-axis
commutants, and the distinction between rare-sector trace deficits and
operator norm. The unrestricted zero-clean five-layer lower bound follows
from the analytic commutant and phase-quantization proof, not from finite
enumeration. Neither the tests nor that theorem establish an unrestricted positive-clean
five-layer bound or a uniform scalar gap.
The [one-clean midpoint proof](../../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md),
R76, has [bounded checks](../../tests/test_one_clean_midpoint_rigidity.py)
of dependent-axis phase layers, diagonal Clifford extraction, native
two-layer return after intermediate mixing, a non-Pauli Clifford midpoint,
compressed-commutant lifting, the eight-boundary recurrence with its
independent coefficient scale, and bounded cyclotomic orbits, polynomial
constants and phase cancellation. Its Pauli and Clifford midpoint bounds
allow arbitrary dirty width and crossing rank. The width-independent
gaps and complementary finite-projective-order obstruction are analytic;
the finite fixtures do not prove the unrestricted one-clean conjecture or
infer a midpoint condition from terminal clean return.
The [inner clean-spread proof](../../research/depth/ONE_CLEAN_INNER_SPREAD.md),
R77, has [finite audits](../../tests/test_one_clean_inner_spread.py)
of flat inner Pauli expansions, their outer Clifford images, the
compressed finite-boundary expansion, and eighth-root coefficient
homogenization. Exact arithmetic in $`\mathbb Q(\sqrt2)`$ verifies
the unequal sector product traces of a native two-qubit counterexample
to a common sign-compatible flip. Literal physical layers independently
check its non-Clifford midpoint, one clean-active inner axis in each
encoder, exact clean return, arbitrary dirty extension, and occupied
versus full-space norms. The width-independent gap and asymptotic
spread condition are analytic; the counterexample's positive constant
defect is not an asymptotic failure of rigidity.

The [outer-overlap proof](../../research/depth/OUTER_OVERLAP_RIGIDITY.md),
R97, and the [representation limits](../../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md)
have [nine exact fixtures](../../tests/test_outer_overlap_rigidity.py).
They reuse the retained $`\mathbb Q(\sqrt2)`$ Pauli arithmetic and
$`\mathbb Q(\sqrt2,i)`$ native-amplitude ring. Signed permutations are
checked with T inverses, inactive wires and changed Clifford coordinates;
an orthogonal nonpermutation checks that pair-overlap conditions alone
do not suffice. A three-qubit non-Clifford midpoint verifies literal
physical routing, shared clean images and a common gauge with exact gate
cancellation. Small R77 tensor families verify one-clean parity encoding,
growing displayed overlap and crossing rank, center sign changes and
complete projector-mixing identities. A literal alternative three-layer
word retains its exact scalar and shows that this displayed complexity
does not obstruct resplitting. A separate four-wire native circuit checks
every occupied entry, zero leakage and determinant $`-i`$, including the
obstruction after one added dirty wire and its disappearance after two.
These are checks of specific signs and scope, not numerical evidence for
an infinite-width gap. The physical normalizer, simultaneous Clifford
extraction, bounded-misalignment estimate and aligned uniform floor are
proved analytically.

The [conditional geometric source](../CONDITIONAL_GEOMETRIC_SOURCE.md)
has [separate checks](../../tests/test_conditional_geometric_source.py) for
prefix preparation, native controlled-H phases, actual inverse, scalar
programming, inactive sectors, and amplification with work return.
Its shallow prefix schedule and uniform resource bounds are proved
analytically; the fixtures do not emit the whole variable-size frame.
The [grouped-prefetch proof](../GROUPED_PROGRAM_PREFETCH.md) has
[nine bounded checks](../../tests/test_grouped_program_prefetch.py) for
literal phase-mask words, native AND and source preparation, inactive
dirty scratch, coherent program erasure through source leakage, and the
internal-enable amplification sign. Exact integer ledgers check the
adaptive group reservations and unchanged precision caps. The reduced
group fixtures use direct diagonal masks and reflections; they are not
a complete native query, reflection, or full-frame emitter.
The [incremental-selector checks](../../tests/test_grouped_selector_reuse.py)
add exact rational one-, two-, and three-stage comparisons with retained
source leakage, exhaustive two-stage inactive arbitrary-work return,
native private-copy phases, and incorrect-cleanup negative controls.
These stages are reduced rational unitaries testing the dependency
interface, not emitted amplified source circuits. The linear group-depth
and unchanged reservation are proved in
[Section 10](../GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables).
The [source-reuse checks](../../tests/test_grouped_source_reuse.py) compare the
two-stage conjugation on all small source/preparation-work inputs,
detect replacement of the conjugated reflection by the original one,
and check constant stale-monitor and one-use-bank leakage on a legal
exact zero-angle row. PREP is native; masks, monitors, and reflections
are reduced operator fixtures. They supply no new reflection emitter
or unrestricted depth lower bound.
The [conditional unary phase-source proof](../UNARY_PHASE_GRADIENT.md) has
[six bounded checks](../../tests/test_unary_phase_gradient.py). Exact integer
checks establish the Karatsuba cyclic bilinear identities on every basis
pair at $`q=1,2,4,8`$, check their rank counts and dense inputs, and verify
the in-place shift and temporary-word return on every four-bit source
string for both directions of every one-hot program. A literal 28-T
scalar guarded trilinear phase gadget checks its active phase and full
inactive arbitrary-work identity. Native source preparations at q equal
to four and eight include the product phases, reversible binary-to-unary
decode, and actual inverse; these small decode words do not implement
the asymptotic parallel tree schedule.

The group checks use reduced operators on the invariant unary source
space. Three stages with unequal row programs test the signed-shift
convention, local suffix enables, all logical columns, and exact inactive
action on every reduced source input. A perturbed source is retained
through all three stages and unprepared with its actual inverse; the
full initialized-isometry error obeys the $`2\delta`$ group bound and
includes the residual source component. Negative controls expose a
missing original-h guard, an incorrect inverse shift, a reversed phase
convention, and substitution of the ideal preparation inverse. These
tests do not emit the complete native bilinear XOR, selected-program
network, scalable group, or full-frame lookup schedule. Their uniform
work-return, count, width, depth, and error statements are analytic
proofs, not extrapolations from the finite native and reduced fixtures.
The [retained-source fusion identities](../../research/depth/RETAINED_SOURCE_FUSION.md#9-retained-source-fusion-without-a-larger-modulus)
have [seven bounded checks](../../tests/test_retained_source_fusion.py). Exact
Gaussian-dyadic coefficients modulo $`z^q-1`$ compare all 512 four-mode
label triples at q equal to eight with independent plane words, including
unitarity. Three additional label triples check the literal determinant
monomials and detect an incorrectly dropped source phase. Integer arithmetic
has checked pre-operation bounds; its tensor products multiply coefficients on one
shared source. Separate numerical checks cover every character of those
words, the existing native magic Clifford, and all 16 states of a physical
four-wire cyclic source bank. Exact eight-mode checks include unequal
children, joint integer root exponents, and the Bell-stabilizer extraction
with logical-only transpose. A noncommuting pair of transported factors
checks the boundary word on Bell inputs and detects its invalid extension
to the complete space. Negative cases retain a spurious source
phase, invert source powers during transpose, or omit the required final
factor. The inactive test checks cancellation of the fixed Clifford
sandwich, without modeling private selector work. A noncommuting pair
rules out only a common fixed basis for one diagonal round. A four-target
fixture partitions sixteen logical modes into five orthogonal completion
sectors and the vacuum. Each completed axis preserves every sector;
two rounds using the original cached labels equal the serial completion
word and its inverse. Independent upper/lower plane words verify the
complete frame factorization on all columns. Recomputing a sector in the
temporary basis, or dropping the upper transport's suffix guard, fails.
This uses exact sector blocks, without simulating arbitrary cache inputs.
General native selectors, work return, resource ledgers, and growing-block
recurrences are analytic claims, not emitted by these fixtures.
The [protected unary source bank](../../research/depth/PROTECTED_UNARY_SOURCE.md) has
[four bounded checks](../../tests/test_protected_unary_source.py). Its q=4
preparation is a literal native word on the full six-wire bank. Two
unequal group partitions are compared on every initial logical/bank
column with both flags zero. The reduced group operators include signed
shifts on every bank bitstring, the initial zero-bank predicate, and
group enables that exclude the protected bank. Inactive interior identity
is checked for $`H=0,h=0`$; arbitrary $`h=1`$ is outside that contract.
Perturbed preparation creates components outside the unary source space;
the actual inverse and final zero-bank predicate act on that live
leakage. The complete output error obeys the single $`2\delta`$ bound
across both groups, including a nonzero final-H residual and a coherent
reference extension. Negative controls retest the physical bank or
substitute the ideal inverse. Predicates and group actions are reduced
operator fixtures, not a native group or QROM emitter; scalable work
return and resource bounds remain analytic.
The [windowed group predicates](../../research/depth/WINDOWED_GROUP_PREDICATES.md) have
[five bounded checks](../../tests/test_windowed_group_predicates.py). Exact
rational comparisons cover one-group partial windows and unequal group
sizes on every small logical/source/tail/dirty-helper input column, with
zero active cache. Direct full-suffix predicates supply the independent
expected enables. Noncommuting reduced bodies retain source leakage and
temporarily use future logical bits as work, restoring them before the
next cache operation. Exhaustive $`H=h=0`$ checks include arbitrary cache
and both dirty helpers. The four-Toffoli constant-arity flag word is
expanded to native gates and checked on all six-wire inputs for literal
phase, helper return, and its actual inverse. Negative controls delay
block-zero erasure, delay chain erasure, or omit the original H guard.
Completed group bodies remain reduced operators; the scalable predicate
emitter, reservation, and asymptotic depth bounds are analytic.

The [shared-prefix query audit](../../tests/test_shared_prefix_query_audit.py)
adds eight native checks on at most eight wires, with independent signed Boolean
expected actions for every input column and literal phase. Two small
bilinear forms expose the retained-prefix dirty residue and check the
actual inverse. A stale-mask example changes the local address between
baseline and loaded evaluations: every cache returns, but the output is
wrong. Completed per-group XOR queries supply a positive reference with
arbitrary dirty programs and retained source/data inputs; a delayed
whole-word echo instead leaves an extra source X. A separate baseline
works under its explicit extra clean-bit promise and fails for arbitrary
baseline input. Full-input comparisons imply coherent/reference safety
for the positive words, with an additional explicit coherent check.
Three additional checks compare the complete original and charged-fusion
words, preserving the controlled rotation's literal sign and reading the
second address after the first target changes. A single nonlinear feature
$`Y\mapsto Y\oplus x_0x_1`$ keeps this fixture small; it is not a
full one-hot indicator. The guarded fixture reuses the actual activity
flag, tests both directions of its transition without an old-flag copy,
and checks the resulting correction on every input. Missing, premature,
and incorrectly guarded corrections are detected. A separate affine-parity
refresh checks all input phases, dirty-helper return, coherent reference
columns, exactly 14 T gates, and eight disjoint T layers in its emitted
schedule. The general fusion and query-resource bounds remain analytic.
These are bounded query-interface circuits and counterexamples, not
complete Hopf groups, scalable query emitters, or lookup lower bounds.

The [chunked dirty indicator](../CHUNKED_DIRTY_INDICATOR.md) has
[three bounded checks](../../tests/test_chunked_dirty_indicator.py) for
native shared-control phases, actual inverses, arbitrary dirty-tree
return, and the selected-path conjugation order. Exact arithmetic checks
the late-layer resource sums. The small row macros do not establish
the imported asymptotic counter depth. Together these fixtures support
the fragile interfaces of the fixed-accuracy complete-frame theorem;
its count, width, and depth bounds are proved analytically.
The [nonuniform-indicator checks](../../tests/test_nonuniform_dirty_indicator.py)
add five bounded audits of unequal chunk partitions. Exact Boolean
polynomials verify the complete leaf-readout echo, actual inverse,
and return of every dirty input wire for small unequal trees. Separate
native matrices check all input phases of the one- and two-bit edge
interfaces composed by those trees; they are not dense native matrices
of the whole tree. Negative cases reverse the root-conjugation order
or omit final cleanup. Exact ceiling and normalized rational ledgers
check the remaining-length recursion, private-pool bound, cubic count
majorant, and logarithmic depth surrogate, including zero address and
large arithmetic-only cases. No exponential register is allocated in
those resource checks. The toy reversible ladders test the edge action;
the shallow conjunction depth and uniform asymptotic bounds come from
the analytic construction.
The [parallel dirty-lookup checks](../../tests/test_parallel_dirty_lookup.py)
audit routed-indicator cancellation, literal native phases, disjoint T
layers, scratch-free width, and arbitrary-input return using symbolic
Boolean polynomials. A reversed-router negative case protects the actual
inverse orientation; a complete native query checks phases and work return. Their [analytic composition](../PARALLEL_DIRTY_LOOKUP.md) retains
the count bound while reducing T-depth under its sufficient dirty-width
condition. The linear table maps still have a charged Clifford-depth cost.
The [bilinear lookup checks](../../tests/test_bilinear_dirty_lookup.py) audit
the separate [two-indicator reduction](../PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction):
rectangular binary basis changes, literal shared-target Toffoli phases,
the four-corner echo, actual inverses, and arbitrary dirty-input return.
These bounded fixtures use the existing routed indicators. They do not
supply the shallow indicator themselves.

The [blocked bilinear checks](../../tests/test_blocked_bilinear_lookup.py)
add five bounded audits of the [selected-block query](../BLOCKED_BILINEAR_LOOKUP.md).
Native matrices cover rank-zero, rank-one, and rank-two controlled
bilinear words, rectangular coordinate changes, arbitrary helper return,
and an eight-wire selected block with a dirty traversal selector. Their
emitted schedules check literal phases, actual inverses, T-count,
T-depth, and disjoint targets within every T layer. Complete queries
with three or four address bits and one or two output bits are checked
as exact Boolean polynomials in every input wire, including both dirty
indicator words, the selector stack, the helper, and arbitrary outputs.
These larger queries use a reduced controlled-bilinear action whose
native interface is checked separately. Negative cases remove the
helper inverse or second traversal, or select the wrong block. The
fixtures do not emit the scalable chunked indicators or establish the
width-sensitive frame theorem by extrapolation; that resource and
composition argument remains analytic.

The [rectangular allocation checks](../../tests/test_rectangular_query_allocation.py)
audit the [uniform-precision allocation](../UNIFORM_PRECISION_DEPTH.md)
with five exact arithmetic checks and the indicator constant normalized
to one. An independent search over feasible powers of two checks the
floor-based choice, address partition, simultaneous indicator pools,
and traversal reservation. Integer comparisons, including squared
positive remainders, check the block-count, selected-count, and indicator
bounds without approximating square roots. Boundary cases include odd
address length, output precision exceeding the row count, width-cap
transitions, and the helper-free single-row Clifford word. An exact
family shows the square allocation's square-root precision penalty in
its selected-middle cost. These finite certificates neither assign
unit constants to physical indicators nor establish asymptotic bounds
by fitted data; they test the allocation proof's arithmetic interfaces.

The [dirty-counter checks](../../tests/test_counter_dirty_indicator.py)
audit the two-adder signed increment, both modular-adder actions, cyclic
routing, nested full-input echoes, actual inverses, and parallel native
layers. Complete counter fixtures emit the linear TTK adder; separate
fixtures check the shortened RV macro. The optimized RV ladder depth is
imported analytically and is not inferred from those serial macro checks.
The [involution-increment checks](../../tests/test_readonly_dirty_increment.py)
separately audit the two-dirty-bit replacement, both literal polarities,
native phases, actual inverses, and shared-address scheduling. They use
serial controlled increments for their bounded fixtures; the logarithmic
depth is the imported analytic contract. Reversing the complement's
position gives decrement, a negative case included in the checks.
The [masked-sum checks](../../tests/test_dirty_sum_interfaces.py) separately
audit the two-controlled-increment compressor, literal native phases,
weighted helper offsets, the complete outer translation echo, actual
inverses, and the static column schedule. Their native increment is a
slower exact MCX expansion; the optimized one-dirty-helper increment
depth is imported from Vandaele's theorem. The
[pipelined-sum checks](../../tests/test_pipelined_dirty_sum.py) separately
check deferred parity forests, emitted doubling-block increments,
updated-prefix carry conditions, private phase helpers, full sum echoes,
and the static bit-release schedule. The refined pipeline uses linear
TTK arithmetic; its overlapping schedule is proved analytically.
The
[analytic counter proof](../PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters),
[compression proof](../DIRTY_SUM_COMPRESSION.md), and hybrid composition
give the improved complete-frame depth bound;
the finite checks do not emit a variable-size complete frame.
The [two-bit depth checks](../../tests/test_dirty_indicator_depth.py) separately
verify the eight-parity CCZ phase identity, both invertible CNOT bases,
all dirty-helper inputs, actual inverses, and the complete six-wire
indicator in exactly two T layers. The matching depth lower bound is
analytic and concerns full-input implementations without initialized work.
The [batched lookup checks](../../tests/test_batched_dirty_lookup.py) cover
guarded partial indicators, symbolic arbitrary-input query return, a
complete native query with literal phases, disjoint T layers, and the
allocation ledger. The [fixed-accuracy frame theorem](../../research/depth/BATCHED_DIRTY_LOOKUP.md)
is an analytic composition; these finite checks do not implement the
general two-dirty-helper MCX or a variable-size full-frame emitter.
The [amortized lookup checks](../../tests/test_amortized_dirty_lookup.py)
audit controlled rectangular linear maps, two-pass dirty traversal,
literal complete-query phases and inverses, all dirty-register return,
and the emitted resource schedule. Their [analytic composition](../AMORTIZED_DIRTY_LOOKUP.md)
gives same-circuit bounds at every accuracy and matching count and depth
in its stated workspace range, including inverse-polynomial error in N.
The full-frame theorem remains an analytic result;
the fixture emits only its small lookup components.
An exact-rational budget check covers the capped source precisions,
including ceiling boundaries, the complete-error margin, unchanged dirty
reservation, source-count savings, and weighted query sums. It adds no
new native fixture and does not certify a smaller total depth.
The variable-accuracy resource check applies the literal power-of-two
bank and indicator allocations, worst-rank native query ledgers, and
full-frame sums, including high precision and the empty/nonempty matching
window boundary. These are finite arithmetic checks of the analytic
composition, not an emitted variable-size frame circuit.
The [Hopf error checks](../../tests/test_hopf_error_accumulation.py) compare
small complete relative-frame spectra with the analytic finite recursion,
including different bases and error signs. A single common representation
of every active and inactive source mask then checks the exact Q/OAA
algebra, its rejected-space sign, the coherent leakage witness, and actual
inverse return. These source fixtures use a joint Clifford-algebra
reduction evaluated in floating point; they are not emitted native words.
Their dense matrices have dimension at most 64, and the larger small-n
witnesses propagate one state through 32-dimensional blocks. The uniform
stability and leakage bounds are proved in the [error chapter](../HOPF_ERROR_ACCUMULATION.md).
The [flag-echo checks](../../tests/test_flag_echo.py) verify the exact square
and Pauli-echo errors, the near-zero improvement, and the full-space
equal-mask exception in one common source representation. The
[radial-filter checks](../../tests/test_radial_filter.py) verify the
fixed-point word with the literal global phase, complex accepted scalar,
full polar error, actual inverse, and six-occurrence phase-error budget.
Their source/filter matrices have dimension 16. Continuous phase
perturbations check the norm budget; they are not emitted native
approximants. The phase-word resource bound and full-frame composition
are analytic. Exact arithmetic checks cover the smaller precision cap.
The [state-based QBP depth proof](../../supplements/state_based_qbp/STATE_QBP_DEPTH.md) composes the same exact
queries and [borrowed-signal rotations](../ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations).
Existing tests above and [rotation checks](../../tests/test_one_clean_compiler.py)
support those primitives. The composition preserves literal phases and
work-return contracts; its simultaneous count/depth ledger is analytic.
No new circuit fixture or general emitter is claimed, and T-depth remains
distinct from total elementary depth.
The [tree-residual checks](../../tests/test_tree_residual_structure.py)
reconstruct complete small residuals from classical tree generators,
including complex coarse words and singular angles.
The [tree-transport checks](../../tests/test_tree_transport.py) also reconstruct
the sparse resolvent, verify the exact transport Gram matrix and complete
three-mode unitary columns, and test subtree overlap-defect telescoping,
weighted norms, and omitted perturbation terms. These support the
[endpoint candidate analysis](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md); they do not
construct a cheaper native joint block.
The [weighted-block checks](../../tests/test_weighted_transport_block.py)
verify its nonorthogonal-column witness, scalar recursion against independent
subtree Schur solves, complete one-flag unitaries and actual inverses, and
weighted-map perturbation bounds at zero and near-zero defects. They also
check the batched allocation, literal level packing and gathering, the
physical marker permutation, and the selected completion's discontinuity.
The [batched native bound](../../research/endpoint/WEIGHTED_TRANSPORT_BLOCK.md#5-a-native-implementation-by-depth-batching)
is proved analytically. These operator matrices have dimension at most 32;
the separate permutation checks enumerate basis labels without dense matrices.
They do not emit the elementary circuit or establish the linear T-count target.
The [residual-assembly checks](../../tests/test_residual_assembly.py) compare the
affine forward recursion with independent subtree Schur solves, reconstruct
the complete affine and reverse dilations, and assemble and amplify their
two-flag selection on every logical input column. They check literal mode
packing and gathering, singular and zero-defect cases, actual inverses,
and perturbations that retain rejected-flag leakage and arbitrary dirty
inputs. Negative controls detect the wrong reverse branch, a changed
relative root phase, and replacement of an actual inverse by a forward call.
The [native assembly and resource bound](../../research/endpoint/RESIDUAL_ASSEMBLY.md) are proved
analytically. These operator fixtures use matrices of dimension at most 64;
they do not emit the native compiler or demonstrate an improved endpoint
T-count.

Exact support and Walsh-column checks in the same residual-assembly suite
verify the independent-forest parametrization and its cube dimension.
Rational capacity checks support the packing ledger. The resulting
$`\Omega(Nn)`$ theorem concerns the enlarged independent family, including
accepted-block-only approximation; it is not a bound on the correlated
Hopf subfamily.
The [affine-tree fusion checks](../../tests/test_affine_tree_fusion.py) compare
the ten-mode and recursive full-input merges with literal local products,
actual inverses, and independent path maps. They test continuation ranks,
growing column support, eager entry counts, and a bottom-only perturbation
whose internal transport is undamped despite a strict external norm margin.
Native paired-source fixtures separately check transformed-mask correlations,
changed-address uncomputation, and rejected-space return when a scalar signal
is reused. Matrices have dimension at most 64. These support the
[fusion audit](../../research/endpoint/RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths)
and [scoped source-reuse arguments](../../research/endpoint/SOURCE_REUSE_LIMITS.md); the mask
T-count bound is analytic, and neither the fixtures nor the representation
witnesses establish a general gate lower bound or a cheaper native compiler.
The [coupled-merge checks](../../tests/test_coupled_residual_merge.py) compare
the complete normalization-two recursion and its unfolded target wrappers
at a fork and height three, with close complex native coarse words and real
targets. They retain every signal port, actual inverse, and correlated
spectator input; they test the stability estimate, anchored mixed term,
and the complete rank-four repair. These small operator fixtures do not
price the transported repair modes or emit a faster native compiler.
The [transported-repair checks](../../tests/test_transported_repair.py) test
the subsequent commutator circuit using only the supplied child and local
parent words. They cover intersecting and vanishing transported modes,
the weighted repair-error bound, exact cancellation with noncanonical
actual words on shared dirty work, and the failure of an unmatched ideal
inverse. The repaired full merge retains the local half-error bound;
these tests do not establish shared precision or a better T-count.
The [shared-conjugator checks](../../tests/test_shared_conjugator_merge.py)
use literal Clifford+T source and routing words on all 128 core/logical/signal
columns. They test valid outer cancellation, the surviving two fork returns,
a direct half-unitary target failure, and the precision-independent
coefficient-ellipse obstruction to retuning the same word. A second flag
is an untouched spectator. The literal source ledger concerns that emitted
word; it is not a general gate lower bound.
The [precision-carry checks](../../tests/test_precision_carry.py) compare literal
native chain and reversed-star loaders, retain their scalar phases, and
check the orientation and gate counts of unequal-width bridges. They verify
three-width full-operator composition after logical address changes,
reversed coefficient grids, and a transformed-mask correlation arising from
a legal grouped coefficient. The
[width-boundary proof](../../research/endpoint/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
supplies the analytic bridge and mask bounds. Cheap outer basis changes do
not price the transformed group bodies or improve the endpoint frontier.
The [correlated-carry checks](../../tests/test_correlated_precision_carry.py)
test a flag/source code on every dirty input, unequal-width physical code
transport, arbitrary reference correlations, and the reuse of released tail
wires as dirty work. A realizable grouped mask exposes code leakage; full
two-syndrome extraction and its actual inverse test source recovery and
its error bound when approximate extraction leaves flag leakage. The
[correlated-boundary proof](../../research/endpoint/SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost)
charges initialization and decoding and restricts that complete
renewal interface. It does not apply the same lower bound to an accepted-only
block or establish an unrestricted carried-source impossibility.
The [small-product checks](../../tests/test_small_product_compilation.py)
compare quaternion products with complete addressed matrices, including
literal phases and variable word lengths. A changed-address control rejects
invalid rowwise preprocessing; overlapping Hopf pairs retain their mixed
path amplitude. The positive compilation statement reuses the existing
one-target theorem with its stated allocation.
The [Cayley checks](../../tests/test_tree_cayley.py) start with four modes and
extend to eight. They compare direct complete residuals with the recursive
skew-Hermitian generator, including actual complex native coarse words,
zero defects, chart boundaries, small-inverse conditioning, and inverse
Cayley stability. This verifies the recursive operator representation, not a
Clifford+T emitter or a reduced precision-source count. The inverse audit
also checks Woodbury cancellation, complete scattering ports, and an
actual native child with non-scalar dirty compression.
The [original-tree checks](../../tests/test_original_tree_residual.py)
extend this evidence to arbitrary ordered trees, oriented path/time signs,
and the determinant identity. Independent complete products test both
off-diagonal directions of every tree cut for the logical residual
$`WQ^\dagger`$ and its Cayley transform, including complex coarse blocks.
The outward Hopf order has rank-one cuts; reordered words include a sharp
rank-two witness. The shared-factor checks reconstruct the Hopf residual
with $`3N-2`$ modes, test the Gram bounds, and retain the $`5N-2`$ count
for separate failure banks as a diagnostic of that construction. The
[regularized joint completion](../../research/endpoint/ORIGINAL_TREE_RESIDUAL.md#51-regularized-joint-completion-on-four-flag-blocks)
fits the legal $`4N`$-mode space; its evidence is mapped
[above](#endpoint-source-and-representation-checks). Neither the separate-bank
count nor the joint matrix construction supplies a native-cost lower or upper
bound. Numerical rank/chart tolerances are diagnostic, not certified equality
decisions.
The [four-mode native checks](../../tests/test_native_cayley.py) retain the
magic Clifford's literal phase and orientation, recover regular Cayley
factors at zero and nonzero defects, and test an addressed native word
that borrows each idle logical target as an active predicate helper.
Every dirty column and both occupied flag sectors enter the helper-return
check; the initialized-flag composition includes nonzero leakage without
a reset. Its q=2 source is a small interface diagnostic, not a certified
fine-precision instance or a full SO(4) circuit emitter. Its explicit
source and T-gate counts include actual inverse appearances.
The [eight-mode checks](../../tests/test_eight_mode_coupling.py) verify the
controlled root's four commuting Pauli factors, the product-distance
witness and attaining product, and the root/child commutator. A special
angle has a literal four-T word with its scalar correction. These tests
check scoped constructions and obstructions, not a generic precision
lower bound. Matrices in the native spectator check have dimension at
most 128; the eight-mode matrices have dimension eight.
The [joint-source checks](../../tests/test_joint_source_body.py) emit complete
eight-mode words and verify their sixteen-mode recurrence on matrices of
dimension 128 and 256. They compare the retained primitive with its new
controlled-Y route on every signal and dirty input, including inactive
predicates and exact return of a used logical helper. Independent encoded
blocks fix the intended angle sign. Tests check literal amplification
cancellation, full signal error, the final parity boundaries, and a
separately charged complex coarse inverse. Source counts are 135 versus
87 and 180 versus 114 before further native simplification. Emitted
loader/seed-simplified eight-mode words have the same leading precision
term in both comparisons; fixed-mask tail cancellation is checked at
q=2, 3, and 5. The q=2 complete words are interface diagnostics, not
fine-precision certification. The [analytic allocation](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
is restricted to local depth three or four; no generic endpoint gain or
optimality of the displayed counts is asserted.
The [canonical group checks](../../tests/test_grouped_scalar_completion.py)
test the real scalar completion, zero coefficients, and direction-controlled
Z selection of the actual inverse. The retained one-tail scalar is checked
against the same symmetry, so the comparison receives that consolidation
too. Native q=2 table words have dimension 256; complete SELECT columns
have dimension 512. The outer group is applied to 64 input columns in a
2048-dimensional space, keeping every logical, borrowed-signal, and core
input with the other work initialized. This avoids forming a large dense
group matrix. Its four-label coefficient table checks half-block algebra,
amplification perturbations, inactive sectors, and common-source/parity
identities; it is not a supplied Hopf residual with an exactly unitary
unrounded target. The [group proof](../../research/endpoint/CANONICAL_SCALAR_COMPLETION.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge)
supplies that specialization and the L-bit error/workspace ledger.
The executed program count is three inner rotations per amplified group.
Literal source counts 135 versus 87 are structural word comparisons,
not a leading precision improvement or a claim of optimality.
The [Hopf scattering checks](../../tests/test_hopf_scattering.py) reconstruct
all columns of four- and eight-mode real frames from explicit tree ports,
including singular and negative-angle cases. They distinguish the
nilpotent physical internal block from its two minus-identity dummy ports.
The input permutation is emitted as X/CNOT/Toffoli gates with each Toffoli
expanded to its phase-correct native word; complete matrices of dimension
16 and 32 verify arbitrary helper return. The output permutation is a
literal X. A deterministic Fourier grid checks the analytically known
right-path coefficient. Coin rotations remain ideal in these fixtures;
the [proof](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
prices one native step and proves the query restriction. Neither the
small checks nor that step ledger supplies a feedback compiler.
The [reference-state QBP checks](../../tests/test_reference_state_qbp.py)
compare every decoded real raw derivative against the analytic Jacobian for
small regular, signed, and singular trees with complex Hermitian-unitary
observables. They test the linear envelope recurrence and its sharp Z=1
and Z=n examples. Ideal small matrices verify the complex residual's
half-amplitude block, the literal sign of one amplification step, the
actual-inverse perturbation estimate, and coherent two-branch preparation
with arbitrary spectator input. These are state-isometry and decoder
checks; the SU(2) coins and coarse matrices in these fixtures are ideal.
The [native cost proof](../../supplements/state_based_qbp/STATE_ONLY_COMPILER.md) reuses separately verified
primitives; no new full elementary emitter or complete-frame endpoint is
claimed by these tests.
The [coarse-frame QBP checks](../../tests/test_coarse_frame_qbp.py) verify
X/Y interference against analytic gradients for complex observables and
regular, signed, and singular real charts. They check the uniform depth
bound and the histogram/adjoint reconstruction. A concrete phase fixture
detects the bias caused by dropping Y measurements or replacing actual-C
scores with ideal Walsh scores. Their coarse words are small ideal
unitaries; the native resource and complete-work claims are analytic in
the [decoder proof](../../supplements/state_based_qbp/COARSE_FRAME_QBP.md) and its preparation dependency.
The [native coarse-frame integration checks](../../tests/test_native_coarse_qbp.py)
now emit the finite-size two-qubit fallback completely in elementary gates.
They compare full preparation isometries and two-by-two decoded score
operators on arbitrary dirty input, preserving literal phase and retaining
all output rows. The helper is used and returned; the two reserved compiler
flags are untouched. The same controlled observable is priced in the
original protocol. The [executable histogram checks](../../tests/test_coarse_frame_decoder.py)
compare reconstruction with independent small Jacobians and preserve
integer cancellations beyond fixed-width and floating-point precision.
Final contractions remain NumPy arithmetic. The [native example note](../../supplements/state_based_qbp/NATIVE_COARSE_QBP.md)
separates these finite checks from the general residual-table construction
and its analytic cost proof.
The [complex coarse-frame checks](../../tests/test_complex_coarse_qbp.py)
verify arithmetic-mean gauge factorization with winding phases, actual
native one-qubit rows repeated across every suffix, and the bounded
complex residual amplification. Independent analytic derivatives check
X/Y magnitude means, executable integer-histogram reconstruction, and the
direct phase-Y stream, including signed and singular trees. Wrong-gauge
and omitted-phase corrections give explicit bias witnesses. The
[complex utility](../../compiler_robust_hopf/complex_coarse_decoder.py) uses
floating-point final coefficients; neither these checks nor the utility
emit the complete native prefix selection or fine residual-table circuit.
Their general resource contracts are proved in the
[compiler](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md) and [decoder](../../supplements/state_based_qbp/COMPLEX_COARSE_QBP.md).
The subsequent [native complex checks](../../tests/test_native_complex_coarse_qbp.py)
emit the complete two-qubit phase-prefix tables, coherent preparation,
actual inverse, and both gradient streams. Independent analytic four-mode
oracles check every active system/branch/helper input in 64-by-four batches;
two-by-two dirty score operators check magnitude and phase means, including
nonzero signals at a singular tuple. The
[exact certificate](../../tests/test_complex_coarse_certificate.py) uses
rational arithmetic in the basis $`(1,\sqrt2,i,i\sqrt2)`$ to verify the
commutator trace and coarse bound $`3/200\lt1/64`$. Complete elementary
gate counts include the same observable in the original protocol. The
[report](../../scripts/complex_coarse_native_example.py) reproduces both streams
without Monte Carlo. This finite-size fallback leaves both compiler flags
unused and does not emit the general fine residual table; see its
[scope and word ledger](../../supplements/state_based_qbp/NATIVE_COMPLEX_COARSE_QBP.md).
The [task-cost arithmetic certificates](../../tests/test_qbp_cost_comparison.py)
check the exact weighted phase-word sum, disjoint core/selector/helper/signal
and bank reservations, a conservative integer query-cost inequality, and
the one-wire complex-bank eligibility gap when the state and original
precisions coincide. They use bounded integer grids without statevectors.
The [comparison chapter](../../supplements/state_based_qbp/QBP_COST_COMPARISON.md) and
[banked state proof](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation)
supply the analytic resource bounds and precision regimes; these finite
checks do not prove asymptotics or an end-to-end gradient advantage.
The [residual preprocessing certificates](../../tests/test_residual_table_preprocessing.py)
use exact rational arithmetic to check square-root enclosures, the paired
half-phase identities, tiny and unit-circle boundary cases, both sides of
the phase cut, and the small-system banked allocation. The
[helper](../../compiler_robust_hopf/residual_table_preprocessing.py) emits
coefficient intervals, not quantum gates. The
[preprocessing proof](../../supplements/state_based_qbp/RESIDUAL_TABLE_PREPROCESSING.md) supplies the uniform
error and bit bounds, and the [bounded-input audit](../../supplements/state_based_qbp/BOUNDED_INPUT_QBP.md)
separately prices native-word searches and the classical gradient baseline.
These checks do not implement a general native emitter or certify a
polynomial-time optimal one-qubit synthesizer.

The [certified residual-row bridge](../../supplements/state_based_qbp/NATIVE_RESIDUAL_ROTATION.md) adds a
production path from these intervals to elementary gates. Its
[exact programming tests](../../tests/test_rotation_programming.py) recover
source moments independently, certify beta and block errors through
2048-bit precision, and cover head, grid, phase-cut, and input boundaries.
The [native tests](../../tests/test_native_residual_rotation.py) compare
literal emitted words against independent Majorana algebra on every input
column of eight wires at q=5; they check phases, amplified moments,
borrowed-signal symmetry, Rz orientation, and residual composition.
Fine-precision error checks use both small Clifford-algebra representations,
while output ledgers verify linear gate storage and 540q T/TDG gates per
row. Certificates and signs altered after programming are rejected.
This implements one unaddressed row, with no general table or complete
state compiler claim. The full-operator bound rests on the existing
source and amplification proofs, not numerical extrapolation.

The [two-row table checks](../../tests/test_native_residual_table.py) extend
this path to one address bit and one enable literal. They check the exact
Toffoli and Clifford mask tables, then propagate all target/core/signal
inputs in each fixed address/enable sector at q=5. The sector evaluator
retains all phases on fixed controls, verifies their exact return, and
compares active blocks to the independent Majorana algebra and inactive
blocks to identity. This covers the full ten-wire operator through four
256-dimensional blocks, including coherent address phases and the
active-zero convention. High-precision checks use rational certificates
and literal counts of 540q+630 T/TDG gates per residual table. This
control arity needs no extra helper; general dirty lookup and larger
predicates are not emitted by it.

The [four-row lookup tests](../../tests/test_native_residual_lookup.py) add
a second address bit without another helper. Small complete native mask
matrices verify constant, linear, and grouped quadratic terms, including
odd X/Z overlap and literal axis order. The q=5 residual table is checked
on all eight address/enable sectors, retaining all 256 target/signal/core
columns per sector. Fixed-to-fixed CNOTs inside the native Toffoli are
tracked as temporary classical address changes with literal phases;
quantum mixing and unrestored controls are rejected. Active words are
compared with independent Majorana algebra and disabled words with
identity. Fine-q certificates and native counts verify the conditional
ledger at most 540q+1050 T/TDG gates. These are bounded component checks,
not a general table or a new asymptotic result.

The [bounded state tests](../../tests/test_native_residual_state.py) compose
two certified tables into a one-system-qubit preparation with two clean
flags. Exact normalized coefficient fixtures cover complex amplitudes
and boundary cases. The small native blocks feed a complete 128-column
dirty-input isometry at q=5, retaining leakage through both reflections
and the actual inverse. Independent elementary propagation checks the
flattened word without phase alignment. Exact high-precision certificates
and the unsimplified 3240q+3793 T/TDG ledger are checked without growing
statevectors. The [scope and proof](../../supplements/state_based_qbp/NATIVE_RESIDUAL_STATE.md) retain the
external coefficient promises and the larger-than-minimum small-system
dirty allocation. This is not a general emitter or an advantage claim.

The [two-system-qubit tests](../../tests/test_native_two_qubit_residual_state.py)
extend preparation to four-row tables and a larger initial reflection.
Every logical/helper input of the 28-T reflection is checked, including
helper-one and reference-entangled inputs. The complete 2048-by-128
initialized isometry at q=5 retains all dirty-input columns and leakage
through the actual inverse and reflections. Native sector blocks and an
independent source-algebra oracle are compared, with direct flattened
propagation on coherent inputs. An exact complex coefficient fixture has
three nonzero tails and satisfies the actual 1/64 residual neighborhood; no
coarse circuit is inferred from that fixture. Fine-q certificates and
counts are checked without a large statevector. The enlarged reflection
borrows an existing arbitrary core wire and returns it exactly; it does
not introduce a third clean compiler flag or a reset.

The [coherent residual-selection tests](../../tests/test_native_branched_residual_state.py)
use one system qubit and a separate arbitrary protocol branch. Their
complete 2048-by-256 initialized isometry at q=5 retains both branch
values and every dirty-input column. Native table sectors and independent
source algebra agree with literal phases; flattened propagation and the
actual inverse also cover coherent reference inputs. Ideal witnesses
detect a changed relative phase and an initial reflection that wrongly
tests the branch. The two main complex fixtures satisfy the actual 1/64
residual neighborhood. A root-reference specialization and separate
normalization checks retain the common-coarse coefficient contract.
Fine-q rational certificates and counts establish the implemented
conditional ledger without large statevectors. The
[branch proof and scope](../../supplements/state_based_qbp/NATIVE_RESIDUAL_BRANCH.md) exclude a supplied
coarse circuit, observable, or full native QBP integration from this component.

The [native residual QBP tests](../../tests/test_native_residual_qbp.py)
connect the fine preparation to an actual coarse word, controlled
Hadamard observable, and both gradient readouts. Independent rational
target/residual formulas and square-root enclosures check the coefficient
promises and full logical coarse-distance certificate. Small coherent
dirty/reference column batches propagate through complete native words
and an independent source-algebra oracle, retaining all flag leakage.
Probabilities sum over every work outcome without postselection; their
score operators check magnitude and raw phase means. Ideal witnesses
detect omitted Y records, the wrong coarse inverse, and a changed
relative branch phase. Exact histogram tests retain finite raw phase
sums, and fine-q ledgers charge the observable and every coarse call.
The [proof and scope](../../supplements/state_based_qbp/NATIVE_RESIDUAL_QBP.md) distinguish these finite
checks from analytic error bounds and sampling guarantees. This pass
does not supply a native fine-frame cost benchmark or a general emitter.

The [antichain checks](../../tests/test_antichain_compiler.py) compare the exact
strict-descendant forest factorization with complete complex tree words,
pack mixed-depth disjoint updates into one last-bit multiplexor, and check
native dirty-Fredkin echoes on every input, including inactive sectors.
They also test reflection-by-reflection predicate echoes and detect the
failure of the stated factorization when updates are comparable or the
forest is omitted. Matrices have dimension at most 64. The
[antichain compiler](../../research/endpoint/ANTICHAIN_COMPILER.md) proves the native resource and
full-operator error bounds analytically; these fixtures neither emit its
fine-precision synthesis circuit nor establish the unrestricted endpoint.
The [sparse-update checks](../../tests/test_sparse_update_compiler.py) verify
ancestor-closed support, exact off-support forest factorization, affine
basis-state transpositions, and the packed operator's identity complement.
They reconstruct the dense column dictionary and detect reversed atom/filter
order and an incorrectly shared rejection flag. Small component matrices
and batched input columns test actual inverses, complete amplification,
inactive-sector identity, and dirty-core perturbations including rejection.
The [sparse-update proof](../../research/endpoint/SPARSE_UPDATE_COMPILER.md) supplies the native
resource and conditional-workspace bounds. These checks are finite interface
tests, not an elementary emitter for the precision-dependent compiler.
The [one-clean checks](../../tests/test_one_clean_compiler.py) reconstruct the
paired-Majorana source and general Pauli masks from native gates, audit the
conjugated scalar word and five-call amplification on all dirty input
columns, and check addressed relative phases and exact inactive action.
They also test X symmetry, the full-operator borrowed-signal estimate,
reference stability, and a scalar-phase counterexample to that extension.
The fine-precision checks use an independent small Clifford-algebra
representation, rather than a large precision-core simulation. The
[one-clean theorem](../ONE_CLEAN_COMPILER.md) supplies the analytic error and
workspace proof; these fixtures do not emit its asymptotic grouped circuit.
The [source-merge checks](../../tests/test_source_merge.py) use native scalar
sources and noncommuting three-level logical operations to exhibit the
Pauli-routed cubic return on every dirty input. These support the scoped
claims in [source-reuse limits](../../research/endpoint/SOURCE_REUSE_LIMITS.md); neither classical
compression nor a failed merge fixture settles the linear frame endpoint.
The direct single-flag merge is also checked against its actual native
branch word, including its constant error and identity-angle amplification
failure.

## Rendering checks

For the optional browser audit, follow the
[rendering guide](../../assets/README.md#rendering-checks). It checks all diagrams
and complete Markdown pages, including table mathematics, display equations,
missing expressions, and horizontal overflow. Equations are checked both with
MathJax SVG and the browser's native MathML layout. The native checks reject
unsupported numbered rows and detect the narrow vertical stacking that can
otherwise pass expression-count and overflow checks. Equation numbers are
printed as ordinary mathematical text so native rendering preserves them.
Comparison signs use HTML-safe TeX commands, guarded against literal
less-than characters that GitHub can misinterpret before math rendering.
Saved desktop and narrow-screen previews support human inspection. These
local checks cover both rendering paths but do not reproduce GitHub's private
renderer or guarantee identical layout across all browsers.
