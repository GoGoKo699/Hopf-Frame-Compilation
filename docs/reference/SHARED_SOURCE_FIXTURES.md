# Shared-source and recovery fixtures

[Verification catalogue](VERIFICATION_CATALOGUE.md#detailed-finite-checks) · [Core verification map](../VERIFICATION.md) · [Tests](../../tests/README.md)

<a id="shared-source-and-recovery-checks"></a>

The [compressed-column theorem](../../research/depth/SUBLINEAR_T_DEPTH.md)
has separate source, exact dirty-query, complete-recovery and tail proofs.
Its five [bounded recovery checks](../../tests/test_sublinear_depth_recovery.py)
use independently written two-target descendant vectors and small numerical
unitary completions. They check every frame column for unequal signed and
singular tuples, marker packing, padded-tag overlap, actual inverse use,
active-only amplification phase, and residual work/flag leakage. Negative
controls preserve the root while corrupting a marker, substitute a different
completion's inverse, or apply the active minus globally. The last error
check includes final guard erasure and exact identity on inactive columns.

Five [binary-gradient calibration checks](../../tests/test_phase_gradient_calibration.py)
reuse the emitted helper-free TTK adder. Exhaustive small register inputs
cover occupied programs, gradient offsets, both target signs, inactivity,
and actual inverses. A two-bit Clifford+T matrix checks literal phases;
larger small kickback blocks check the real rotation. A coarse native
product source produces nonzero leakage and verifies its complete-output
error bound. Negative controls drop a guard or target sign, or substitute
an ideal preparation inverse. The source is deliberately coarse, not an
implementation of arbitrary-accuracy one-qubit synthesis.

The [shared-gradient checks](../../tests/test_shared_gradient_recovery.py)
add a literal two-bit compute/swap/uncompute adder, occupied-work inactive
echoes, and a complete coefficient-bank phase table. Recovery-shaped
noncommuting words test one native preparation and its actual inverse at
the two boundaries, including nonzero final leakage and a preparation-error
bound independent of the number of middle calls. Counterchecks reject a
reflection that tests the live gradient, a substituted ideal inverse, and
a discarded coefficient-bank scalar phase. Programmable angle accuracy
requires a zero program; the fixed outer damping instead ignores program
contents. These checks use small exact permutations and matrices of
dimension at most 64. They test the local identities and source boundary
comparison, not a complete emitted amplitude-amplification compiler.

The [grouped-terminal checks](../../tests/test_grouped_terminal_depth.py)
combine literal native H/T/T-dagger interpretation with complete-word
prefetch and actual unloading. They test occupied inactive programs,
predicates and dirty helpers, with negative controls for an omitted phase
correction, a mutable query address and a missing activity guard. Integer
checks cover rounded group sizes, execution order and the source/coarse
balance. Their small column batches and illustrative allocation constants
do not establish a scalable emitter or a numerical crossover.

The [shared binary terminal checks](../../tests/test_shared_binary_terminal.py)
add an exact sparse two-stage, sixteen-type accepted product on every
source/logical column. Noncommuting coarse words produce coherent prefix
addresses for the next stage; failure amplitudes and literal complex phases
are retained. The source is unprepared only at the final boundary. Missing
history and intermediate source projection give different accepted blocks.
This algebra fixture supplies coefficient-table actions and folds the local
suffix predicate into atom rejection; it does not emit the full query or
production-precision circuit. Separate native shift/history words check all
occupied register inputs, dirty-helper return, literal phases and actual
inverses. A native source/body/inverse word checks inactive cancellation on
occupied coherent inputs. Wrong inverse ordering, a raw disabled source,
and preparation substituted for its inverse have negative controls. The
pinned generic source/shift/history receipts remain unchanged.

These fixtures do not emit the derived Boolean tables or a complete native
descendant loader. The existing native dirty-query, shared-control phase,
returned-helper predicate and nonuniform-indicator tests cover their stated
finite interfaces. Terminating source selection, polynomial and arithmetic-table overhead
absorption, and the same-circuit asymptotic theorem are proved analytically;
no numerical extrapolation or efficient preprocessing claim is made.

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

