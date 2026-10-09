# Verification catalogue

[Core verification map](../VERIFICATION.md) · [Source catalogue](SOURCE_CATALOGUE.md) · [Tests](../../tests/README.md)

This catalogue preserves detailed fixture contracts and software boundaries.
The core map links Results A–D to their proofs and principal evidence; this
page records the smaller input, precision, workspace, and gate-emission
scopes of the individual checks. A research fixture's presence here does
not make its construction a current principal result. Finite checks do
not prove universal statements or asymptotic optimality.

## Catalogue index

| Material | Entry point | Scope |
|---|---|---|
| Separate state-based QBP result | [Coverage and software boundary](#state-based-qbp-coverage) | State/reference preparation with a changed decoder |
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
identifying the missing interface. The prescribed complete-frame
endpoint and general T-depth frontier remain
[active research questions](../OPEN_PROBLEM.md#full-circuit-requirement);
count and depth now match in explicit accuracy/workspace ranges by the
[uniform-precision construction](../UNIFORM_PRECISION_DEPTH.md). The existing
state-based result does not require the remaining questions' resolution.
Application-level advantage is outside the current research scope; the
existing resource comparisons and classical baselines remain documented.

### Detailed finite checks

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

The [native representation tests](../../tests/test_native_representation_limits.py)
verify realification gate words in exact cyclotomic arithmetic, the
rectangular norm and full-isometry identity, bounded integer circle
solutions, and literal marker-column gaps. Exact character and phase-grid
enumeration checks Spin multiplicities at n=3 and n=4, retaining the
small-dimension exception. The
[realification](../../research/endpoint/NATIVE_REALIFICATION.md) and
[Spin-decoding](../../research/endpoint/SPIN_DECODING_LIMITS.md) proofs
give their arbitrary-size conclusions under the stated support and
wrapper hypotheses. No unrestricted circuit lower bound is inferred.

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
Counterexamples test why nilpotence
and complete-output accuracy cannot be omitted. These support the scoped
analytic restrictions; they do not establish a full-frame impossibility.

The [structural endpoint checks](../../tests/test_endpoint_structural_limits.py)
test the [full-frame metric](../HOPF_INTERFACE.md#41-complete-frame-parameter-geometry),
[normalization and residual witnesses](../../research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md),
singular-safe nested projectors, restricted reflection and displacement
bounds, native control phases, query parity, and source-return counterexamples.
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
Cayley stability. This verifies the new operator representation, not a
Clifford+T emitter or a reduced precision-source count. The inverse audit
also checks Woodbury cancellation, complete scattering ports, and an
actual native child with non-scalar dirty compression.
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
