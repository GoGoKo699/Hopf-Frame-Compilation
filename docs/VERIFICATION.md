# Proof and executable correspondence

[← QBP consequence](QBP_CONSEQUENCE.md) · [Complete narrative](../REVIEW.md) · [Source map →](SOURCE_MAP.md)

The repository uses executable checks to test the conventions and constructions
that are easiest to get wrong: column labels, chronological order, relative
phases, coherent routing, work-register cleanup, and simultaneous resource
peaks.  The asymptotic theorem itself is established analytically.

## 1. Evidence levels

| Level | Meaning in this repository |
|---|---|
| analytic proof | a dimension-independent operator or resource argument |
| explicit construction | a complete reversible or logical-gate schedule is supplied |
| imported exact synthesis | an elementary compiler theorem is used under its stated model |
| finite regression check | independently constructed matrices, states, schedules, or ledgers are compared |

These levels are kept distinct.  A matrix test does not prove an asymptotic
bound, and an asymptotic bound does not certify an implementation's bit order.

### State-based QBP coverage

The [task theorem](STATE_BASED_QBP_THEOREM.md) and the bounded native
residual-to-gradient integration are complete within their stated scopes.
The following map separates each analytic claim from its executable
coverage. The detailed finite checks below retain their individual
input, precision, and workspace boundaries.

| Claim | Proof home | Executable coverage | Software boundary |
|---|---|---|---|
| Common coarse reference with literal phase and exact dirty return | [Complex coarse compiler §§1–6](COMPLEX_COARSE_COMPILER.md) | [Real](NATIVE_COARSE_QBP.md), [complex](NATIVE_COMPLEX_COARSE_QBP.md), and [residual QBP](NATIVE_RESIDUAL_QBP.md) fixtures record actual coarse words | No integrated coarse synthesizer for all admitted inputs |
| Fine single-state and coherent-pair isometries, including returned work | [State compiler §§1–7](STATE_ONLY_COMPILER.md) and [banked extension](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) | Certified rows, two- and four-row tables, one- and two-qubit preparations, and a coherent one-qubit selector | No variable-size or banked native state emitter |
| Original raw gradient means, confidence, and reconstruction accuracy | [Real decoder](COARSE_FRAME_QBP.md) and [complex decoder §§2–6](COMPLEX_COARSE_QBP.md) | General histogram utilities test the identities; [the residual fixture](NATIVE_RESIDUAL_QBP.md) has exact rational decoders | General coefficient contractions still use floating point; the guarded decoder is proved, not implemented |
| Compiler T-count and simultaneous T-depth upper schedules | [Task ledger](STATE_BASED_QBP_THEOREM.md#4-quantum-and-classical-resource-ledger) and [depth proof](STATE_QBP_DEPTH.md) | Exact resource ledgers, native dirty-lookup checks, and bounded literal gate inventories | No full native program realizing the variable-size banked depth schedule |
| Polynomial construction and explicit program-output bound under the input restriction | [Bounded-input audit §§3–4](BOUNDED_INPUT_QBP.md#3-which-native-searches-are-polynomial-in-the-input-parameters) | [Residual intervals](RESIDUAL_TABLE_PREPROCESSING.md), exact sign programming, and bounded emitted words | No integrated certified front end from the admitted bounded Hopf inputs to the full program |
| Fine residual preparation inside both complete gradient streams | [Bounded residual QBP integration](NATIVE_RESIDUAL_QBP.md) | [Native words and rational decoders](../compiler_robust_hopf/native_residual_qbp.py), checked by [independent oracles](../tests/test_native_residual_qbp.py) | Fixed complex one-qubit target; no same-accuracy native fine-frame benchmark or sampled-accuracy experiment |

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
identifying the missing interface. The prescribed complete-frame
endpoint and general T-depth frontier remain
[active research questions](OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule);
count and depth now match in explicit accuracy/workspace ranges by the
[amortized construction](AMORTIZED_DIRTY_LOOKUP.md). The existing
state-based result does not require the remaining questions' resolution.
Application-level advantage is outside the current research scope; the
existing resource comparisons and classical baselines remain documented.

### Detailed finite checks

The [operator-source proof](OPERATOR_SOURCE_COMPILER.md) is accompanied by
[finite circuit checks](../tests/test_operator_source_compiler.py).
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

The [source-reuse limits](SOURCE_REUSE_LIMITS.md) have
[separate finite checks](../tests/test_source_reuse_limits.py) for nilpotent
contractions, arbitrary encoding bases, dirty-dimension independence, and
transformed-mask operator identities. Counterexamples test why nilpotence
and complete-output accuracy cannot be omitted. These support the scoped
analytic restrictions; they do not establish a full-frame impossibility.

The [conditional-suffix proof](CONDITIONAL_SUFFIX_COMPILER.md) has
[focused grouped-block checks](../tests/test_conditional_suffix_compiler.py).
They combine a native small operator source with separate scalar and
atom flags, retain every dirty input column through amplification,
and check literal phases, the disjoint ancestor-column support partition,
forward and reverse column blocks, inactive suffix sectors, and predicate
uncomputation after leakage. A separate integer ledger tests exponentially
growing groups and the dirty-width/error inequalities. Its sample
grouping constants do not certify the unspecified fixed constants of the
coarse synthesis primitive. These are finite interface checks, not a native
elementary circuit emitter for the whole asymptotic grouped compiler.

The [T-depth schedule](T_DEPTH_COMPILER.md) has
[native routing checks](../tests/test_t_depth.py) for literal four-layer
shared-control Fredkins, disjoint T-layer supports, actual inverses, and
whole dirty-bank queries. Clifford layers may have substantial depth;
the fixture does not treat T-depth as total execution depth.
The [exact source-depth proof](SOURCE_T_DEPTH.md) has
[separate checks](../tests/test_source_t_depth.py) for its source schedules,
literal native phases, and exact denominator witnesses. The matching lower
bounds apply only to the defined Majorana-layer architecture and exact
targets, including its arbitrary returned dirty extensions. They are not
general Hopf or controlled-source depth lower bounds.
The [parallel dirty-lookup checks](../tests/test_parallel_dirty_lookup.py)
audit routed-indicator cancellation, literal native phases, disjoint T
layers, scratch-free width, and arbitrary-input return using symbolic
Boolean polynomials. A reversed-router negative case protects the actual
inverse orientation; a complete native query checks phases and work return. Their [analytic composition](PARALLEL_DIRTY_LOOKUP.md) retains
the count bound while reducing T-depth under its sufficient dirty-width
condition. The linear table maps still have a charged Clifford-depth cost.
The [bilinear lookup checks](../tests/test_bilinear_dirty_lookup.py) audit
the separate [two-indicator reduction](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction):
rectangular binary basis changes, literal shared-target Toffoli phases,
the four-corner echo, actual inverses, and arbitrary dirty-input return.
These bounded fixtures use the existing routed indicators. They do not
supply the shallow indicator themselves.
The [dirty-counter checks](../tests/test_counter_dirty_indicator.py)
audit the two-adder signed increment, both modular-adder actions, cyclic
routing, nested full-input echoes, actual inverses, and parallel native
layers. Complete counter fixtures emit the linear TTK adder; separate
fixtures check the shortened RV macro. The optimized RV ladder depth is
imported analytically and is not inferred from those serial macro checks.
The [involution-increment checks](../tests/test_readonly_dirty_increment.py)
separately audit the two-dirty-bit replacement, both literal polarities,
native phases, actual inverses, and shared-address scheduling. They use
serial controlled increments for their bounded fixtures; the logarithmic
depth is the imported analytic contract. Reversing the complement's
position gives decrement, a negative case included in the checks.
The [masked-sum checks](../tests/test_dirty_sum_interfaces.py) separately
audit the two-controlled-increment compressor, literal native phases,
weighted helper offsets, the complete outer translation echo, actual
inverses, and the static column schedule. Their native increment is a
slower exact MCX expansion; the optimized one-dirty-helper increment
depth is imported from Vandaele's theorem. The
[pipelined-sum checks](../tests/test_pipelined_dirty_sum.py) separately
check deferred parity forests, emitted doubling-block increments,
updated-prefix carry conditions, private phase helpers, full sum echoes,
and the static bit-release schedule. The refined pipeline uses linear
TTK arithmetic; its overlapping schedule is proved analytically.
The
[analytic counter proof](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters),
[compression proof](DIRTY_SUM_COMPRESSION.md), and hybrid composition
give the improved complete-frame depth bound;
the finite checks do not emit a variable-size complete frame.
The [two-bit depth checks](../tests/test_dirty_indicator_depth.py) separately
verify the eight-parity CCZ phase identity, both invertible CNOT bases,
all dirty-helper inputs, actual inverses, and the complete six-wire
indicator in exactly two T layers. The matching depth lower bound is
analytic and concerns full-input implementations without initialized work.
The [batched lookup checks](../tests/test_batched_dirty_lookup.py) cover
guarded partial indicators, symbolic arbitrary-input query return, a
complete native query with literal phases, disjoint T layers, and the
allocation ledger. The [fixed-accuracy frame theorem](BATCHED_DIRTY_LOOKUP.md)
is an analytic composition; these finite checks do not implement the
general two-dirty-helper MCX or a variable-size full-frame emitter.
The [amortized lookup checks](../tests/test_amortized_dirty_lookup.py)
audit controlled rectangular linear maps, two-pass dirty traversal,
literal complete-query phases and inverses, all dirty-register return,
and the emitted resource schedule. Their [analytic composition](AMORTIZED_DIRTY_LOOKUP.md)
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
The [Hopf error checks](../tests/test_hopf_error_accumulation.py) compare
small complete relative-frame spectra with the analytic finite recursion,
including different bases and error signs. A single common representation
of every active and inactive source mask then checks the exact Q/OAA
algebra, its rejected-space sign, the coherent leakage witness, and actual
inverse return. These source fixtures use a joint Clifford-algebra
reduction evaluated in floating point; they are not emitted native words.
Their dense matrices have dimension at most 64, and the larger small-n
witnesses propagate one state through 32-dimensional blocks. The uniform
stability and leakage bounds are proved in the [error chapter](HOPF_ERROR_ACCUMULATION.md).
The [flag-echo checks](../tests/test_flag_echo.py) verify the exact square
and Pauli-echo errors, the near-zero improvement, and the full-space
equal-mask exception in one common source representation. The
[radial-filter checks](../tests/test_radial_filter.py) verify the
fixed-point word with the literal global phase, complex accepted scalar,
full polar error, actual inverse, and six-occurrence phase-error budget.
Their source/filter matrices have dimension 16. Continuous phase
perturbations check the norm budget; they are not emitted native
approximants. The phase-word resource bound and full-frame composition
are analytic. Exact arithmetic checks cover the smaller precision cap.
The [state-based QBP depth proof](STATE_QBP_DEPTH.md) composes the same exact
queries and [borrowed-signal rotations](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations).
Existing tests above and [rotation checks](../tests/test_one_clean_compiler.py)
support those primitives. The composition preserves literal phases and
work-return contracts; its simultaneous count/depth ledger is analytic.
No new circuit fixture or general emitter is claimed, and T-depth remains
distinct from total elementary depth.
The [tree-residual checks](../tests/test_tree_residual_structure.py)
reconstruct complete small residuals from classical tree generators,
including complex coarse words and singular angles.
The [tree-transport checks](../tests/test_tree_transport.py) also reconstruct
the sparse resolvent, verify the exact transport Gram matrix and complete
three-mode unitary columns, and test subtree overlap-defect telescoping,
weighted norms, and omitted perturbation terms. These support the
[endpoint candidate analysis](ENDPOINT_TREE_TRANSPORT.md); they do not
construct a cheaper native joint block.
The [weighted-block checks](../tests/test_weighted_transport_block.py)
verify its nonorthogonal-column witness, scalar recursion against independent
subtree Schur solves, complete one-flag unitaries and actual inverses, and
weighted-map perturbation bounds at zero and near-zero defects. They also
check the batched allocation, literal level packing and gathering, the
physical marker permutation, and the selected completion's discontinuity.
The [batched native bound](WEIGHTED_TRANSPORT_BLOCK.md#5-a-native-implementation-by-depth-batching)
is proved analytically. These operator matrices have dimension at most 32;
the separate permutation checks enumerate basis labels without dense matrices.
They do not emit the elementary circuit or establish the linear T-count target.
The [residual-assembly checks](../tests/test_residual_assembly.py) compare the
affine forward recursion with independent subtree Schur solves, reconstruct
the complete affine and reverse dilations, and assemble and amplify their
two-flag selection on every logical input column. They check literal mode
packing and gathering, singular and zero-defect cases, actual inverses,
and perturbations that retain rejected-flag leakage and arbitrary dirty
inputs. Negative controls detect the wrong reverse branch, a changed
relative root phase, and replacement of an actual inverse by a forward call.
The [native assembly and resource bound](RESIDUAL_ASSEMBLY.md) are proved
analytically. These operator fixtures use matrices of dimension at most 64;
they do not emit the native compiler or demonstrate an improved endpoint
T-count.
The [affine-tree fusion checks](../tests/test_affine_tree_fusion.py) compare
the ten-mode and recursive full-input merges with literal local products,
actual inverses, and independent path maps. They test continuation ranks,
growing column support, eager entry counts, and a bottom-only perturbation
whose internal transport is undamped despite a strict external norm margin.
Native paired-source fixtures separately check transformed-mask correlations,
changed-address uncomputation, and rejected-space return when a scalar signal
is reused. Matrices have dimension at most 64. These support the
[fusion audit](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths)
and [scoped source-reuse arguments](SOURCE_REUSE_LIMITS.md); the mask
T-count bound is analytic, and neither the fixtures nor the representation
witnesses establish a general gate lower bound or a cheaper native compiler.
The [coupled-merge checks](../tests/test_coupled_residual_merge.py) compare
the complete normalization-two recursion and its unfolded target wrappers
at a fork and height three, with close complex native coarse words and real
targets. They retain every signal port, actual inverse, and correlated
spectator input; they test the stability estimate, anchored mixed term,
and the complete rank-four repair. These small operator fixtures do not
price the transported repair modes or emit a faster native compiler.
The [transported-repair checks](../tests/test_transported_repair.py) test
the subsequent commutator circuit using only the supplied child and local
parent words. They cover intersecting and vanishing transported modes,
the weighted repair-error bound, exact cancellation with noncanonical
actual words on shared dirty work, and the failure of an unmatched ideal
inverse. The repaired full merge retains the local half-error bound;
these tests do not establish shared precision or a better T-count.
The [shared-conjugator checks](../tests/test_shared_conjugator_merge.py)
use literal Clifford+T source and routing words on all 128 core/logical/signal
columns. They test valid outer cancellation, the surviving two fork returns,
a direct half-unitary target failure, and the precision-independent
coefficient-ellipse obstruction to retuning the same word. A second flag
is an untouched spectator. The literal source ledger concerns that emitted
word; it is not a general gate lower bound.
The [precision-carry checks](../tests/test_precision_carry.py) compare literal
native chain and reversed-star loaders, retain their scalar phases, and
check the orientation and gate counts of unequal-width bridges. They verify
three-width full-operator composition after logical address changes,
reversed coefficient grids, and a transformed-mask correlation arising from
a legal grouped coefficient. The
[width-boundary proof](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
supplies the analytic bridge and mask bounds. Cheap outer basis changes do
not price the transformed group bodies or improve the endpoint frontier.
The [correlated-carry checks](../tests/test_correlated_precision_carry.py)
test a flag/source code on every dirty input, unequal-width physical code
transport, arbitrary reference correlations, and the reuse of released tail
wires as dirty work. A realizable grouped mask exposes code leakage; full
two-syndrome extraction and its actual inverse test source recovery and
its error bound when approximate extraction leaves flag leakage. The
[correlated-boundary proof](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost)
charges initialization and decoding and restricts that complete
renewal interface. It does not apply the same lower bound to an accepted-only
block or establish an unrestricted carried-source impossibility.
The [small-product checks](../tests/test_small_product_compilation.py)
compare quaternion products with complete addressed matrices, including
literal phases and variable word lengths. A changed-address control rejects
invalid rowwise preprocessing; overlapping Hopf pairs retain their mixed
path amplitude. The positive compilation statement reuses the existing
one-target theorem with its stated allocation.
The [Cayley checks](../tests/test_tree_cayley.py) start with four modes and
extend to eight. They compare direct complete residuals with the recursive
skew-Hermitian generator, including actual complex native coarse words,
zero defects, chart boundaries, small-inverse conditioning, and inverse
Cayley stability. This verifies the new operator representation, not a
Clifford+T emitter or a reduced precision-source count. The inverse audit
also checks Woodbury cancellation, complete scattering ports, and an
actual native child with non-scalar dirty compression.
The [four-mode native checks](../tests/test_native_cayley.py) retain the
magic Clifford's literal phase and orientation, recover regular Cayley
factors at zero and nonzero defects, and test an addressed native word
that borrows each idle logical target as an active predicate helper.
Every dirty column and both occupied flag sectors enter the helper-return
check; the initialized-flag composition includes nonzero leakage without
a reset. Its q=2 source is a small interface diagnostic, not a certified
fine-precision instance or a full SO(4) circuit emitter. Its explicit
source and T-gate counts include actual inverse appearances.
The [eight-mode checks](../tests/test_eight_mode_coupling.py) verify the
controlled root's four commuting Pauli factors, the product-distance
witness and attaining product, and the root/child commutator. A special
angle has a literal four-T word with its scalar correction. These tests
check scoped constructions and obstructions, not a generic precision
lower bound. Matrices in the native spectator check have dimension at
most 128; the eight-mode matrices have dimension eight.
The [joint-source checks](../tests/test_joint_source_body.py) emit complete
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
fine-precision certification. The [analytic allocation](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
is restricted to local depth three or four; no generic endpoint gain or
optimality of the displayed counts is asserted.
The [canonical group checks](../tests/test_grouped_scalar_completion.py)
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
unrounded target. The [group proof](CONDITIONAL_SUFFIX_COMPILER.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge)
supplies that specialization and the L-bit error/workspace ledger.
The executed program count is three inner rotations per amplified group.
Literal source counts 135 versus 87 are structural word comparisons,
not a leading precision improvement or a claim of optimality.
The [Hopf scattering checks](../tests/test_hopf_scattering.py) reconstruct
all columns of four- and eight-mode real frames from explicit tree ports,
including singular and negative-angle cases. They distinguish the
nilpotent physical internal block from its two minus-identity dummy ports.
The input permutation is emitted as X/CNOT/Toffoli gates with each Toffoli
expanded to its phase-correct native word; complete matrices of dimension
16 and 32 verify arbitrary helper return. The output permutation is a
literal X. A deterministic Fourier grid checks the analytically known
right-path coefficient. Coin rotations remain ideal in these fixtures;
the [proof](ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer)
prices one native step and proves the query restriction. Neither the
small checks nor that step ledger supplies a feedback compiler.
The [reference-state QBP checks](../tests/test_reference_state_qbp.py)
compare every decoded real raw derivative against the analytic Jacobian for
small regular, signed, and singular trees with complex Hermitian-unitary
observables. They test the linear envelope recurrence and its sharp Z=1
and Z=n examples. Ideal small matrices verify the complex residual's
half-amplitude block, the literal sign of one amplification step, the
actual-inverse perturbation estimate, and coherent two-branch preparation
with arbitrary spectator input. These are state-isometry and decoder
checks; the SU(2) coins and coarse matrices in these fixtures are ideal.
The [native cost proof](STATE_ONLY_COMPILER.md) reuses separately verified
primitives; no new full elementary emitter or complete-frame endpoint is
claimed by these tests.
The [coarse-frame QBP checks](../tests/test_coarse_frame_qbp.py) verify
X/Y interference against analytic gradients for complex observables and
regular, signed, and singular real charts. They check the uniform depth
bound and the histogram/adjoint reconstruction. A concrete phase fixture
detects the bias caused by dropping Y measurements or replacing actual-C
scores with ideal Walsh scores. Their coarse words are small ideal
unitaries; the native resource and complete-work claims are analytic in
the [decoder proof](COARSE_FRAME_QBP.md) and its preparation dependency.
The [native coarse-frame integration checks](../tests/test_native_coarse_qbp.py)
now emit the finite-size two-qubit fallback completely in elementary gates.
They compare full preparation isometries and two-by-two decoded score
operators on arbitrary dirty input, preserving literal phase and retaining
all output rows. The helper is used and returned; the two reserved compiler
flags are untouched. The same controlled observable is priced in the
original protocol. The [executable histogram checks](../tests/test_coarse_frame_decoder.py)
compare reconstruction with independent small Jacobians and preserve
integer cancellations beyond fixed-width and floating-point precision.
Final contractions remain NumPy arithmetic. The [native example note](NATIVE_COARSE_QBP.md)
separates these finite checks from the general residual-table construction
and its analytic cost proof.
The [complex coarse-frame checks](../tests/test_complex_coarse_qbp.py)
verify arithmetic-mean gauge factorization with winding phases, actual
native one-qubit rows repeated across every suffix, and the bounded
complex residual amplification. Independent analytic derivatives check
X/Y magnitude means, executable integer-histogram reconstruction, and the
direct phase-Y stream, including signed and singular trees. Wrong-gauge
and omitted-phase corrections give explicit bias witnesses. The
[complex utility](../compiler_robust_hopf/complex_coarse_decoder.py) uses
floating-point final coefficients; neither these checks nor the utility
emit the complete native prefix selection or fine residual-table circuit.
Their general resource contracts are proved in the
[compiler](COMPLEX_COARSE_COMPILER.md) and [decoder](COMPLEX_COARSE_QBP.md).
The subsequent [native complex checks](../tests/test_native_complex_coarse_qbp.py)
emit the complete two-qubit phase-prefix tables, coherent preparation,
actual inverse, and both gradient streams. Independent analytic four-mode
oracles check every active system/branch/helper input in 64-by-four batches;
two-by-two dirty score operators check magnitude and phase means, including
nonzero signals at a singular tuple. The
[exact certificate](../tests/test_complex_coarse_certificate.py) uses
rational arithmetic in the basis $`(1,\sqrt2,i,i\sqrt2)`$ to verify the
commutator trace and coarse bound $`3/200\lt1/64`$. Complete elementary
gate counts include the same observable in the original protocol. The
[report](../scripts/complex_coarse_native_example.py) reproduces both streams
without Monte Carlo. This finite-size fallback leaves both compiler flags
unused and does not emit the general fine residual table; see its
[scope and word ledger](NATIVE_COMPLEX_COARSE_QBP.md).
The [task-cost arithmetic certificates](../tests/test_qbp_cost_comparison.py)
check the exact weighted phase-word sum, disjoint core/selector/helper/signal
and bank reservations, a conservative integer query-cost inequality, and
the one-wire complex-bank eligibility gap when the state and original
precisions coincide. They use bounded integer grids without statevectors.
The [comparison chapter](QBP_COST_COMPARISON.md) and
[banked state proof](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation)
supply the analytic resource bounds and precision regimes; these finite
checks do not prove asymptotics or an end-to-end gradient advantage.
The [residual preprocessing certificates](../tests/test_residual_table_preprocessing.py)
use exact rational arithmetic to check square-root enclosures, the paired
half-phase identities, tiny and unit-circle boundary cases, both sides of
the phase cut, and the small-system banked allocation. The
[helper](../compiler_robust_hopf/residual_table_preprocessing.py) emits
coefficient intervals, not quantum gates. The
[preprocessing proof](RESIDUAL_TABLE_PREPROCESSING.md) supplies the uniform
error and bit bounds, and the [bounded-input audit](BOUNDED_INPUT_QBP.md)
separately prices native-word searches and the classical gradient baseline.
These checks do not implement a general native emitter or certify a
polynomial-time optimal one-qubit synthesizer.

The [certified residual-row bridge](NATIVE_RESIDUAL_ROTATION.md) adds a
production path from these intervals to elementary gates. Its
[exact programming tests](../tests/test_rotation_programming.py) recover
source moments independently, certify beta and block errors through
2048-bit precision, and cover head, grid, phase-cut, and input boundaries.
The [native tests](../tests/test_native_residual_rotation.py) compare
literal emitted words against independent Majorana algebra on every input
column of eight wires at q=5; they check phases, amplified moments,
borrowed-signal symmetry, Rz orientation, and residual composition.
Fine-precision error checks use both small Clifford-algebra representations,
while output ledgers verify linear gate storage and 540q T/TDG gates per
row. Certificates and signs altered after programming are rejected.
This implements one unaddressed row, with no general table or complete
state compiler claim. The full-operator bound rests on the existing
source and amplification proofs, not numerical extrapolation.

The [two-row table checks](../tests/test_native_residual_table.py) extend
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

The [four-row lookup tests](../tests/test_native_residual_lookup.py) add
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

The [bounded state tests](../tests/test_native_residual_state.py) compose
two certified tables into a one-system-qubit preparation with two clean
flags. Exact normalized coefficient fixtures cover complex amplitudes
and boundary cases. The small native blocks feed a complete 128-column
dirty-input isometry at q=5, retaining leakage through both reflections
and the actual inverse. Independent elementary propagation checks the
flattened word without phase alignment. Exact high-precision certificates
and the unsimplified 3240q+3793 T/TDG ledger are checked without growing
statevectors. The [scope and proof](NATIVE_RESIDUAL_STATE.md) retain the
external coefficient promises and the larger-than-minimum small-system
dirty allocation. This is not a general emitter or an advantage claim.

The [two-system-qubit tests](../tests/test_native_two_qubit_residual_state.py)
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

The [coherent residual-selection tests](../tests/test_native_branched_residual_state.py)
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
[branch proof and scope](NATIVE_RESIDUAL_BRANCH.md) exclude a supplied
coarse circuit, observable, or full native QBP integration from this component.

The [native residual QBP tests](../tests/test_native_residual_qbp.py)
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
The [proof and scope](NATIVE_RESIDUAL_QBP.md) distinguish these finite
checks from analytic error bounds and sampling guarantees. This pass
does not supply a native fine-frame cost benchmark or a general emitter.

The [antichain checks](../tests/test_antichain_compiler.py) compare the exact
strict-descendant forest factorization with complete complex tree words,
pack mixed-depth disjoint updates into one last-bit multiplexor, and check
native dirty-Fredkin echoes on every input, including inactive sectors.
They also test reflection-by-reflection predicate echoes and detect the
failure of the stated factorization when updates are comparable or the
forest is omitted. Matrices have dimension at most 64. The
[antichain compiler](ANTICHAIN_COMPILER.md) proves the native resource and
full-operator error bounds analytically; these fixtures neither emit its
fine-precision synthesis circuit nor establish the unrestricted endpoint.
The [sparse-update checks](../tests/test_sparse_update_compiler.py) verify
ancestor-closed support, exact off-support forest factorization, affine
basis-state transpositions, and the packed operator's identity complement.
They reconstruct the dense column dictionary and detect reversed atom/filter
order and an incorrectly shared rejection flag. Small component matrices
and batched input columns test actual inverses, complete amplification,
inactive-sector identity, and dirty-core perturbations including rejection.
The [sparse-update proof](SPARSE_UPDATE_COMPILER.md) supplies the native
resource and conditional-workspace bounds. These checks are finite interface
tests, not an elementary emitter for the precision-dependent compiler.
The [one-clean checks](../tests/test_one_clean_compiler.py) reconstruct the
paired-Majorana source and general Pauli masks from native gates, audit the
conjugated scalar word and five-call amplification on all dirty input
columns, and check addressed relative phases and exact inactive action.
They also test X symmetry, the full-operator borrowed-signal estimate,
reference stability, and a scalar-phase counterexample to that extension.
The fine-precision checks use an independent small Clifford-algebra
representation, rather than a large precision-core simulation. The
[one-clean theorem](ONE_CLEAN_COMPILER.md) supplies the analytic error and
workspace proof; these fixtures do not emit its asymptotic grouped circuit.
The [source-merge checks](../tests/test_source_merge.py) use native scalar
sources and noncommuting three-level logical operations to exhibit the
Pauli-routed cubic return on every dirty input. These support the scoped
claims in [source-reuse limits](SOURCE_REUSE_LIMITS.md); neither classical
compression nor a failed merge fixture settles the linear frame endpoint.
The direct single-flag merge is also checked against its actual native
branch word, including its constant error and identity-angle amplification
failure.

## 2. What is represented locally

| Component | Local representation | Principal executable check | Imported ingredient |
|---|---|---|---|
| Hopf states and frames | dense matrices from independent recursive and addressed-layer constructions | complete equality, orthogonality, state and marker columns, chart domains, singular coordinates | Hopf geometry inherited from the earlier work |
| strict-zero echo | dense logical gates and exact reversible permutations | all four sectors, every nonfinal depth, full frames, inverse, and complex composition | UCG and ancilla-free MCT resource theorems |
| binary–one-hot decoder | explicit X/CNOT/Toffoli layers | basis action, arbitrary-basis reversibility, disjoint layers, counts, and clean return | constant-cost elementary decompositions |
| coherent router | explicit CNOT-fanout and Fredkin layers with sparse complex-state simulation | basis routing, entangled inputs, tail direct sum, complete cut, and zero leakage | coherent copy and constant-cost controlled gates |
| controlled subtree frames | exact logical token/flag-controlled rotations | equality to the ideal block-diagonal tail and flag cleanup | UCG and MCT synthesis |
| leaf-phase diagonal | exact block-diagonal UCG matrix | complete diagonal, inverse, common phase, and complex magnitude composition | UCG synthesis |
| resource theorem | integer and exact-rational ledgers | workspace peaks, schedule dispatch, endpoints, geometric sums, and cut inequalities | asymptotic primitive bounds |
| gradient records | parity, signed-histogram, and fast Walsh–Hadamard decoders | agreement of routes, empirical means, and deterministic record norms | QBP record identities |

The repository does not reproduce the elementary UCG or multi-controlled-X
compiler.  Those exact results are imported from the all-workspace
state-preparation framework.  Toffoli, Fredkin, controlled one-qubit gates, and
fixed-width controlled Givens rotations have exact constant-size,
constant-depth decompositions in the declared arbitrary-one-qubit+CNOT model.

## 3. Geometry and compiler contract

The frame suite checks:

- equality of the recursive and addressed-layer frames;
- real-frame orthogonality;
- phase-dressed complex magnitude-frame unitarity;
- the breadth-first marker convention;
- $g_{j,j}=a_j^2$ for unrestricted angles;
- $a_j=\sqrt{g_{j,j}}$ on the canonical domains;
- tolerance-aware regularity at floating-point chart boundaries;
- zero raw derivative and a unit chart-selected marker continuation at a
  singular coordinate.

The compiler-boundary fixtures check:

- two unitaries with the same prepared state but different marker columns;
- the resulting gradient change
  ```math
  (2,0,0)\longmapsto(0,\sqrt2,0);
  ```
- a checkpoint suffix that preserves one state but changes a derivative;
- an active-interface-safe suffix that preserves designated means without
  preserving the complete distribution.

Additional boundary tests exercise the sharp worst-observable sensitivity,
common-phase cancellation, projected marker error versus physical leakage,
and singular-marker ambiguity across all two-qubit Pauli observables.

Files:

- [frame implementation](../compiler_robust_hopf/frames.py)
- [frame tests](../tests/test_frames.py)
- [complex geometry tests](../tests/test_complex_analysis.py)
- [compiler-boundary implementation](../compiler_robust_hopf/compiler_boundaries.py)
- [compiler-boundary tests](../tests/test_compiler_boundaries.py)

## 4. Strict-zero construction

The strict-zero suite checks:

1. $C^2=R_y(\theta)$ and $XCX=C^{-1}$;
2. all four $(h,b)$ sectors;
3. exact restoration of the borrowed logical suffix bit;
4. absence of hidden workspace;
5. every nonfinal addressed depth through $n=8$;
6. complete real frames through $n=8$;
7. inverse frames and phase-dressed complex magnitude frames;
8. the endpoints $n=1$, $d=0$, $d=n-2$, and the final depth.

The resource audit checks, using exact arithmetic,

```math
\sum_{q=2}^{n}\frac{2^q}{q}
\leq6\frac{2^n}{n}
```

and the absorption of the polynomial predicate terms into $O(2^n/n)$.

Files:

- [strict-zero construction](../compiler_robust_hopf/strict_zero_echo.py)
- [operator tests](../tests/test_strict_zero_echo.py)
- [exact-rational resource audit](../compiler_robust_hopf/strict_zero_audit.py)
- [resource-audit tests](../tests/test_strict_zero_audit.py)

## 5. Binary–one-hot prefix decoder

The decoder is tested as a complete reversible permutation, not only on its
intended clean input.  The suite verifies:

- $\lvert x\rangle\lvert0\rangle$ maps to $\lvert0\rangle\lvert e_x\rangle\lvert0\rangle$ and returns under the inverse;
- arbitrary computational-basis contents return after forward and inverse;
- every declared layer has disjoint wire support;
- the exact workspace formula
  ```math
  3\,2^t-2-t;
  ```
- the exact gate and depth formulas;
- the one-hot Givens network equals the complete prefix Hopf frame on the code.

Files:

- [decoder schedule](../compiler_robust_hopf/tree_decoder.py)
- [decoder tests](../tests/test_tree_decoder.py)

## 6. Coherent route–operate–unroute

The router is supplied as an explicit schedule.  The tests cover:

- branch-data, token, copy, and reusable-flag register allocation;
- balanced prefix fanout and exact uncopy;
- disjoint Fredkin layers;
- every clean basis input routed to its prefix-selected branch;
- forward route followed by inverse route on arbitrary complex inputs;
- prefix–suffix-entangled inputs;
- token-controlled action of all subtree frames;
- cleanup of every branch flag before inverse routing;
- equality to
  ```math
  \bigoplus_r W_s^{(r)};
  ```
- equality of the complete routed cut to the direct Hopf frame;
- zero probability outside the clean-workspace subspace.

The explicit schedule is cross-checked against

```math
\text{copy wires}=(2^t-1)(s+1)-t,
```

```math
\text{forward Fredkins}=(2^t-1)(s+1).
```

The simulator iterates over branches for convenience.  The circuit-depth claim
uses the declared parallel schedule: branch data, tokens, flags, and subtree
frames have disjoint physical support.

Files:

- [router and sparse-state simulator](../compiler_robust_hopf/router.py)
- [router operator tests](../tests/test_router.py)
- [tree identities](../compiler_robust_hopf/tree_structure.py)

## 7. Resource theorem

The ledgers use integer or exact-rational checks rather than fitted slopes.  They
cover:

- strict zero workspace;
- the one-clean-flag endpoint;
- the direct-to-routed transition;
- maximal feasible cuts over broad $n,m$ grids;
- copy-pool reuse as branch flags;
- the $s=1$ routed endpoint;
- arbitrarily large workspace;
- matching real-state parameter and light-cone lower bounds.

Representative checked implications are

```math
n^2
=O\left(\frac{2^n}{n+m}\right)
\qquad(1\leq m\lt 4n),
```

and

```math
\frac{2^s}{s}
=O\left(1+\frac{2^n}{n+m}\right)
```

for the maximal routed cut, including the saturated endpoint $`s=1`$.

Files:

- [all-workspace resource selection](../compiler_robust_hopf/unified_compiler.py)
- [resource diagnostics](../compiler_robust_hopf/resource_bounds.py)
- [general resource tests](../tests/test_resource_bounds.py)
- [unified compiler tests](../tests/test_unified_compiler.py)

## 8. Gradient records

The decoder suite compares:

- direct record-wise parity averages;
- dense Walsh transforms;
- the fast Walsh–Hadamard implementation;
- direct phase-stream signed one-hot records.

It verifies the deterministic norm-two record property used by the
coordinatewise concentration statement.

Files:

- [decoders](../compiler_robust_hopf/decoders.py)
- [decoder tests](../tests/test_decoders.py)

## 9. Reproduce the checks

Use Python 3.11 or 3.13.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the short orientation:

```bash
python scripts/reviewer_walkthrough.py
```

Run the deterministic unittest collection and the separate exact
fault-tolerant receipt suites:

```bash
python validate.py
python scripts/verify_fault_tolerant.py
```

The first command discovers the files in `tests/`; the second reproduces the
four retained receipts in `verification/fault_tolerant/`.

Print the two resource ledgers:

```bash
python scripts/unified_resource_ledger.py --n 12
python scripts/strict_zero_echo_ledger.py --n 12
```

Check the recorded upstream versions without network access:

```bash
python scripts/check_upstream_sync.py --offline
```

For the optional browser audit, follow the
[rendering guide](../assets/README.md#rendering-checks). It checks all diagrams
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

## 10. Evidence boundary

The repository does not use finite experiments to establish asymptotic
optimality.  It also does not test:

- device connectivity or routing overhead;
- a hardware-native gate set;
- a general elementary-gate emitter for the full asymptotic Clifford+T compiler;
- noisy execution or readout mitigation;
- application-specific controlled-observable implementations;
- arbitrary non-Hopf differential frames.

The proof should be assessed in four separate steps:

1. the Hopf operator identities;
2. the logical circuits and workspace cleanup;
3. the imported synthesis bounds and resource sums;
4. the matched QBP output and access conventions.

## 11. Fault-tolerant evidence and approximate QBP

The [fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md) separates the analytic
resource proof from the [focused finite checks](../verification/fault_tolerant/README.md).
Run `python scripts/verify_fault_tolerant.py` from the repository root. Four
stdlib-only suites reproduce the retained exact Gray-source, shift-kernel,
resource, and reflection receipts without changing their saved reference files.
The runner compares scientific receipt fields and separately verifies the
deterministic source-file hashes, whose values changed when the sources were
relocated into this repository.

| Standalone receipt suite | What it checks | Evidence boundary |
|---|---|---|
| Gray geometric source | Emitted logical source words, actual inverse, scratch return, reference witnesses, and gate-count recurrences | Finite clean-scratch fixtures; the full source is not emitted as an elementary Clifford+T circuit |
| Shift kernel | Complete finite kernel columns, actual inverses, source defects, failure tracking, and negative controls | Small fixed parameters, without the production residual construction or emitted outer amplification |
| Shift resource | Exact rational scales, coefficient rounding, scalar error budgets, and resource ledgers | Analytic proxies and inequalities; the rows are not measured gate counts |
| Geometric reflection | Small source/reflection matrices, exact Pauli-transfer arithmetic, and algebraic norm separation | Supporting source-cost evidence and a kernel-suite dependency; no additive lower bound follows from these fixtures |

The [receipt map](../verification/fault_tolerant/README.md) records the exact
fixture ranges, source dependencies, and provenance for all four suites.

The [approximation bridge](QBP_APPROXIMATION.md) has tests for full-input error,
actual-adjoint transfer, leakage, reference-entangled borrowed inputs, and
bounded estimator bias. These validate the finite examples behind the general
proof, not differentiation of a synthesized family.

The approximation suite also checks the complex phase stream directly,
reflection-term sampling, pathwise finite-weight error, and a two-shot
measurement instrument whose retained dirty state changes the next-shot law.
The last fixture tests the conditional-mean premise of the martingale proof;
it is not a fitted statistical scaling experiment.

---

[← QBP consequence](QBP_CONSEQUENCE.md) · [Complete narrative](../REVIEW.md) · [Source map →](SOURCE_MAP.md)
