# Research results for complete-frame compilation

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Endpoint status](../docs/OPEN_PROBLEM.md) · [Core proof map](../docs/README.md)

This catalogue separates selected proof dependencies, retained special
cases, earlier methods, and scoped obstructions. The
[core proof map](../docs/README.md#c-one-clean-compilation-and-its-corollaries)
gives the minimal dependency chain. The [claim ledger](../manuscript/PUBLICATION_SCOPE.md)
states the resource formulas and literal reservations; the
[resource frontier](../docs/OPEN_PROBLEM.md) states the remaining questions.
A method-specific obstruction is not a lower bound for every compiler.

The selected growing-source count theorem and
[precision-uniform subtree depth theorem](depth/UNIFORM_SUBTREE_DEPTH.md)
cover complete real and phase-dressed complex magnitude frames with one
clean qubit. The depth theorem combines the
[subtree-projector decoder](depth/SUBTREE_PROJECTOR_DECODER.md) with the
terminal and fallback interfaces of the retained
[one-clean predecessor](depth/ONE_CLEAN_DEPTH.md). The
[sublinear-depth construction](depth/SUBLINEAR_T_DEPTH.md) keeps its
two-clean linear-precision statement and supporting source proofs. The
[state-based QBP supplement](../supplements/state_based_qbp/README.md)
uses a different decoder and has its own task-level guarantee.

The [current research checkpoint](depth/COMPUTATIONAL_PHASE_SESSION_RESULT.md)
connects native probes to computational target derivatives and closes the
active-center coordinate step using Boolean reality and a paid Clifford.
Global-center matching remains open: fixing a center fiber can destroy
native independence. A fully native Majorana Walsh probe retains excessive
count and depth; its Gaussian source obstruction does not cover joint
unrestricted synthesis. R96 remains the strongest general bound, R79
remains selected Result D, and the complete-frame frontier is unchanged.

## Endpoint results

<a id="precision-uniform-two-clean-count"></a>

### Precision-uniform one-clean count

The [one-clean growing-source theorem](endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
is the proof home for Result C's unbanked count. Its
[two-clean predecessor](endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
supplies the precision/mode separation. The one-clean lift borrows the
phase signal while retaining the initialized anchor marker; the layerwise
baseline covers low precision and bounded dimension.

The [claim ledger](../manuscript/PUBLICATION_SCOPE.md) records each family's
base dirty reservation and the matching qualification at that exact width.
The guarantee includes every logical input, dirty/reference correlation,
literal phase, clean-work leakage, and actual inverse.

### Optimal banked count frontier

The [banked growing-source theorem](endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
proves Result C's optimal count frontier throughout the stated doubled-base
width regimes. Banked source masks return every bank and selector before
the next operation. Summed native charges and the existing width-sensitive
lower bounds establish matching; the theorem's Clifford budget and
real/complex reservations remain distinct obligations.

### Constant-clean high-precision count

The [endpoint statement](../docs/OPEN_PROBLEM.md#constant-clean-high-precision-count)
and [growing-source compiler](endpoint/GROWING_SOURCE_COMPILER.md) retain
the two-clean high-precision theorem for complete prescribed real frames.
The [complex-magnitude corollary](endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
composes the half-error real construction with the packed literal diagonal.
Its sharper complex reservation saves one dirty qubit relative to the
one-clean theorem at this endpoint. The real zero-phase witness supplies
the same lower bound. Growing sources, sealed rejection modes, frozen
output tags, actual dirty replay, and final recovery are charged. These
count results do not claim arbitrary-complex-unitary synthesis, a new
T-depth bound, or a general native emitter.

Its proof uses the [mixed local gauges and dyadic kernels](endpoint/HETEROGENEOUS_PRECISION_SOURCE.md),
[disjoint-block compilation of the fixed bottom levels](endpoint/DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks),
and [disjoint source extensions](endpoint/DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources),
together with the linked conditional-work guards, scalar, and recovery
primitives. The conditional-suffix grouped count remains a separate
corollary; C does not require its iterated grouping in the precision splice.

The remaining tables catalogue earlier methods and their individual limits.
An unconstructed source or fine step in one of these methods is not an
unresolved part of the growing-source theorem or a new work order.

<a id="native-primitives-and-promised-families"></a>

#### Retained primitives and promised families

| Study | Result | Scope |
|---|---|---|
| [Disjoint native blocks](endpoint/DISJOINT_BLOCK_COMPILER.md) | Linear endpoint count for disjoint SU blocks; unbanked and banked high-precision block/state interfaces with full Clifford ledgers | High-precision forms require $`L\gt N/n^2`$ beyond their fixed dimension threshold; each interface retains its own clean allocation, literal phase, and complete return |
| [Rank-weighted boundary compiler](endpoint/RANK_WEIGHTED_BOUNDARY_COMPILER.md) | Linear native count under a support-times-rank budget; boundary-column and collective common-basis fusion | Two-clean complex and one-clean real domains; growing sequential overlap does not inherit the disjoint budget |
| [Borrowed phase and recovery sectors](endpoint/CONJUGATE_SIGNAL_RECOVERY.md) | Full conjugate-sector phase word, real initialized-target and source-overlap lifts, zero-clean real reflection banks, and one-clean real blocks | The phase signal is borrowed outside every core; the initialized anchor marker and genuine logical recovery coordinates are retained; complementary-sector restrictions keep their original scope |
| [Borrowed selected axes](endpoint/BORROWED_AXIS_COMPILER.md) | Full-input selected Pauli/source-axis rotations with two arbitrary helpers; jointly compiled raw rejection orientations at O(Q+q+ell) T-count | Literal query phases and shared-core width are charged; raw canonicalization does not compile transported logical prefixes |
| [Matching and leaf deflation](endpoint/MATCHING_LEAF_DEFLATION.md) | Canonical rank equals twice forest matching size; matching-budget compilers and an explicitly closed nonidentity deflation | Paid input correction retained; separate-level precision still costs O(nN); no unrestricted endpoint improvement |
| [Precision-gauge compiler](endpoint/PRECISION_GAUGE_COMPILER.md) | Linear endpoint count on regular charts, fixed-ratio exponent intervals, and the stated cheap/ultraweak family; exact mixed-height scheduling bound | Rescheduling unchanged merge banks retains repeated precision cost; the next construction changes those sources |
| [Heterogeneous precision source](endpoint/HETEROGENEOUS_PRECISION_SOURCE.md) | At every integer height, $`O(N+\mathcal B_K)`$ count after certified clipping; linear retained bottom trace gives linear count | Fixed common-source construction with mixed local gauges and actual dirty replay; includes the earlier exponent-envelope and mixed-height families |
| [Structured forests](endpoint/STRUCTURED_FOREST_COMPILERS.md) | Linear endpoint count under sub-binary degree-window growth, including fixed computational outdegree; over a cheap or fixed-ratio base, either bounded exceptional path depth or bounded exits per base component suffices | Bounded exits allow arbitrarily deep exceptional paths; original computational-label tree and fixed promise constants, with exact signed and corner tuples included |
| [Antichain changes](endpoint/ANTICHAIN_COMPILER.md) | **Proved special case:** linear endpoint count for a literal native-baseline antichain promise | Generic rounded frames need not satisfy that promise |
| [Sparse nested changes](endpoint/SPARSE_UPDATE_COMPILER.md) | **Proved special case:** linear endpoint count under the stated ancestor-closure condition, including a changed path | This is a restricted update family, not the unrestricted frame family |

<a id="full-frame-reductions-and-remaining-native-steps"></a>

#### Earlier full-frame reductions and their native limits

| Study | Result | Scope |
|---|---|---|
| [Fixed-tree Cayley reduction](endpoint/TREE_CAYLEY_REDUCTION.md) | All-angle phase gauge, regular parameter chart, exact O(N)-T permutations, and a two-diagonal-algebra generator | This route's direct fine-core realization costs O(nN); it supplies no separate linear native ledger |
| [Boundary propagation](endpoint/BOUNDARY_PROPAGATION.md) | Sparse preconditioning, repeated-call precision within the exact width, and complete feedback cancellation to one prefix encoder | Repeated calls cost O(RN); the surviving fine encoder has a direct O(nN)-T certificate |
| [Coarse encoder](endpoint/COARSE_PREFIX_ENCODER.md) and [collective refinement](endpoint/COLLECTIVE_PRECISION_REFINEMENT.md) | Joint prefix compilation, exact geometric history, quadratic physical purification, and a cubic midpoint replacement at O(N) T-count | Fixed-order refinement of the terminal regular-core frame; this refinement route has no all-order linear-cost ledger |
| [Coupled tree resolvent](endpoint/COUPLED_TREE_RESOLVENT.md) | All-order graded inverse, tail and support bounds, ideal normalization-three scattering, and hidden-prefix baseline identity | Direct native synthesis retains O(N+nL) cost; small rank and static feedback do not establish a linear native ledger |
| [Canonical-star resolvent](endpoint/CANONICAL_STAR_RESOLVENT.md) | Exact root-column transfer, triangular Cayley denominator, and a linear-count full native local oracle | The descendant feedback inverse remains unpaid |
| [Original-tree residual](endpoint/ORIGINAL_TREE_RESIDUAL.md) | Signed path and determinant identities; exponentially accurate exactly returned coarse word; regularized 4N-mode source completion with bounded tree-cut ranks | Ideal source geometry is proved; fine native source preparation and its inverse remain unpaid |
| [Ordered block banks and finite histories](endpoint/ORDERED_BLOCK_BANKS.md) | Linear-volume complete-frame bands, sharp block/monomial factor-count scale, one nested-support source construction, and native right-readout, coarse-prefix and local fine banks | Stored banks are not ordered execution; fine left source remains open; the modular lower bound does not cover fused native synthesis |
| [Tree transport](endpoint/ENDPOINT_TREE_TRANSPORT.md) | Sparse generators, weighted norms, complete transport columns, and native small-mode constructions | Representation size and classical conditioning do not price the joint native operation; unchanged local scattering retains its query limitation |
| [Weighted transport block](endpoint/WEIGHTED_TRANSPORT_BLOCK.md) | Valid complete dilation with exact/approximate dirty-return variants | The constructed native word retains the repeated precision cost |
| [Packed rejection banks](endpoint/PACKED_REJECTION_BANKS.md) | Supported rejected paths fit the two mode flags; preweight and one final recovery give the complete frame; depthwise balance gives exact four-port compression and a full rounded comparison within the literal dirty width | Fine prefix transport remains unpaid; bankwise scalarization still requires a fine chronological prefix, and direct source cost remains O(nN) |
| [Joint history compression](endpoint/JOINT_HISTORY_COMPRESSION.md) | Three compressor calls recover a complete fine prefix; an anchored Householder specifies a stable full completion, and two compressor calls give its history range reflection | Native compression remains unconstructed; logical-label-preserving mixing has a constant error floor, and a reflection-only oracle needs growing queries to recover the orientation |
| [Residual assembly](endpoint/RESIDUAL_ASSEMBLY.md) | Forward/reverse and coupled repairs, independent-forest packing, orthogonal-history capacity, and damping restrictions | The selected affine branch, native repair and refinement-table costs remain charged; all-order inverse bounds are in the coupled-resolvent chapter |

#### Source programs and encodings

| Study | Result | Scope |
|---|---|---|
| [Source reuse](endpoint/SOURCE_REUSE_LIMITS.md) | Transformed-mask fusion, dirty Klein dressing, Spin closure, full-source dirty-polar obstruction, and a single involutory echo bound | The source witness is not proved reachable by the fixed global encoder; the echo bound requires invariant clean inputs and a commuting target |
| [Source-program resynthesis](endpoint/SOURCE_PROGRAM_RESYNTHESIS.md) | Joint raw-source resynthesis, prefix-selected tails, completed-group gauges, and a fixed-core endpoint normal form | Raw connectors and growing local costs remain distinct; independent gauges apply only after full group recovery; endpoint bound unchanged |
| [Joint native words](endpoint/JOINT_NATIVE_WORDS.md) | Literal Pauli-word contract and returned-dirty compression to half the dirty symplectic rank, without added T gates or error | Canonical counting gives only a linear lower bound; no short-word selection rule or unrestricted improvement |
| [Canonical grouped scalar](endpoint/CANONICAL_SCALAR_COMPLETION.md) | A canonical scalar and direction-controlled inverse fit the existing group interface | Consolidating the scalar calls also helps the old construction; precision is still paid per group |
| [Paid logical encodings](endpoint/PAID_LOGICAL_ENCODINGS.md) | Bit-optimal vacuum/one-hot permutation, a 15-T six-bit native encoder, branching support and exact sector controls | Complete input action retained; physical common-core and regular phase-stage limits are scoped, not unrestricted T bounds |
| [Native realification and ring height](endpoint/NATIVE_REALIFICATION.md) | Realification at at most twice the T-count, exact marker-support gap, linear-height dyadic approximation, and a normalized accepted-phase bound | Ring height does not price a native word; the phase bound requires literally calibrated stabilizer boundaries |
| [Spin decoding](endpoint/SPIN_DECODING_LIMITS.md) | Spectral multiplicity bounds for one doubled lift with a target-dependent inverse encoder or unrelated Clifford wrappers | The endpoint witnesses require n at least four; general noninverse native wrappers and other representations remain eligible |
| [Schur source limits](endpoint/SCALAR_SCHUR_LIMITS.md) | Exact coarse-frame correlation gauge and a polynomial approximation floor at fixed source-mode dimension, even with arbitrary dirty work | One prescribed coarse call with bounded logical-label-preserving sources; no unrestricted source or T-count lower bound |

#### Diagonal and reflection architectures

| Study | Result | Scope |
|---|---|---|
| [Bounded diagonal factorization](endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md) | Certified factor search at the literal dirty budget, stable nested-projector coverage, and an explicit nonflat n-slot factorization | The n-slot realization costs O(nN); uniformly bounded coverage by eligible mixers remains a separate premise |
| [Structural compilation limits](endpoint/STRUCTURAL_COMPILATION_LIMITS.md) | Haar/flat-Clifford restrictions, gauge-transport rigidity, and the tight n-stage bound for differentiable Clifford-centered diagonal architectures | Approximate stage bounds additionally require a physical Lipschitz constant; other baseline and selection models remain eligible |
| [Center-phase reflections](endpoint/CENTERED_PHASE_REFLECTIONS.md) | Exact all-table optimum, real-image classification, quadratic floors, and an exactly verified five-qubit member | The two slots are a table and its inverse; multiple centers and different source words remain eligible |
| [Programmed reflections](endpoint/PROGRAMMED_REFLECTIONS.md) | Exact phase-orbit and reversing-equation tests, a regular signed-cycle chart, and native fixed offsets | Classical factor selection is distinct from native synthesis; the scoped query restriction is paired with a cheap alternate baseline |

The [archived route assessments](ROUTE_HISTORY.md) retain additional
derivations and their assumptions. The independent
[single-angle precision lower bound](../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
and [full-frame geometry](../docs/HOPF_INTERFACE.md) live in their canonical
proof chapters; they are not duplicated here.

### Large-workspace T-depth

The [precision-uniform subtree theorem](depth/UNIFORM_SUBTREE_DEPTH.md),
R79, is the selected all-precision simultaneous count/depth proof for both
frame families. With absolute constants, one clean qubit and the retained
family-specific dirty thresholds, it gives optimal
$`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, and

```math
D_T=O\!\left(NL/b^2+
\min\{nL,(L+\log_2(n+2))\log_2(n+L+2)\}\right).
```

At sufficient $`b=\Theta(\sqrt{NL})`$, depth is
$`O(\log^2(n+2))`$ for $`L=O(\log(n+2))`$, including
inverse-polynomial error in n. R78's fixed-accuracy theorem remains a
corollary; R72's stronger high-precision clause is retained.
The precision splice, literal widths, matching intervals, and shrinking-error
special cases are stated in the [claim ledger](../manuscript/PUBLICATION_SCOPE.md).
The full initialized-isometry guarantee includes common phase, earlier
leakage, arbitrary dirty inputs, and their references. Optimal large-width
depth remains open. At fixed integer $`L\ge6`$, the later
[dirty catalogue theorem R90](depth/DIRTY_CATALOGUE_TERMINAL.md)
improves the general bound to

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell),\qquad \ell=\log_2(n+2).
```

It retains one clean qubit and the original literal real/complex widths,
including arbitrary node-dependent signed magnitudes, zeros, prescribed
leaf phases, and the full initialized-isometry contract. Three R89 coarse
reductions leave a small terminal; a complete rounded program selects
fixed local factors on two borrowed arrays whose unwanted actions cancel
in a commutator. R79 remains the selected manuscript result. R90
does not assert depth optimality or constants uniform in growing L.

The [R92 absolute precision envelope](depth/UNIFORM_CATALOGUE_PRECISION.md)
and [R93 grouped terminal](depth/GROUPED_CATALOGUE_PRECISION.md)
give independent growing-precision refinements. The later
[portable-workspace theorem R94](depth/PORTABLE_WORKSPACE_CATALOGUE.md)
proves, uniformly on $`6\le L\le2\ell`$,

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right),
```

with the original one-clean optimal T-count, Clifford budget, literal
real/complex widths and full operator contract. At $`L=\Theta(\ell)`$
and square-root dirty width, it reaches the
$`O(\ell\log_2(\ell+2))`$ benchmark. The
[executed design](depth/CATALOGUE_BOOTSTRAP_PROMPT.md) and
[session outcome](depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md) record the
bootstrap: R79 compiles a constant number of larger groups directly to final
accuracy, with its complete local workspace carried through selected
transfer. R91 supplies parallel queries and paired routing; no residual
endpoint labels are needed. The [frontier statement](../docs/OPEN_PROBLEM.md#large-width-t-depth)
retains the other precision ranges and the open optimal-depth question.
The [executed full-flag design](depth/FULL_FLAG_DEPTH_PROMPT.md) and
[session outcome](depth/FULL_FLAG_DEPTH_SESSION_RESULT.md) establish
[R95](depth/LOGARITHMIC_PRECISION_DEPTH.md):
$`D_T=O(NL/b^2+L+\ell)`$ on the same range, with the same simultaneous
count, width and full-isometry guarantees. A coarse mathematical frame
implemented to final accuracy admits a constant-normalization sparse
correction. Borrowed real rotations reduce a local source's initialized
scratch to logarithmic size; exact compact references and two explicit
local reductions then give depth $`O(h+\kappa)`$ with much larger
allowed local dirty resources. R94's full-payload catalogue absorbs them.
The [executed precision bridge](depth/BORROWED_PRECISION_BRIDGE_SESSION_RESULT.md)
establishes [R96](depth/ALL_PRECISION_LINEAR_DEPTH.md): the same bound
for every $`n\ge1,L\ge6`$, with absolute constants and the original
count, width and complete-error guarantees. The globally priced stage
separates logarithmic precision scratch from geometric metadata; a
shorter exact cap speeds the original reference. A directly programmed
terminal then avoids the full final-accuracy catalogue. Optimal general
depth and matching dimension-growing lower bounds remain open.

The [executed four-layer session](depth/ONE_CLEAN_FOUR_LAYER_SESSION_RESULT.md)
establishes [outer-overlap rigidity, R97](depth/OUTER_OVERLAP_RIGIDITY.md).
A positive error floor depends only on the smaller outer-axis misalignment
dimension: intersection dimension minus the number of individual axes
from that list lying in it. An intersection spanned by a subset of either
physical list gives R76's uniform floor at any dimension. The unrestricted
gap remains open when both misalignment dimensions grow.
The [executed design](depth/ONE_CLEAN_FOUR_LAYER_PROMPT.md) and outcome
retain the exact descent and center-mixing limits. This is a constant
lower-bound branch; the upper theorem and selected Result D are unchanged.

The [dirty-query round](depth/DIRTY_QUERY_SESSION_RESULT.md) establishes
[conditional R98](depth/QUERY_DEPTH_FRAME_TRANSFER.md): an improved exact
restored-address word-query interface transfers to general fixed-accuracy depth,
with a fully paid $`O_L(\log\log(n+4))`$ residual. The
[paired permutation wrapper](depth/PAIRED_PERMUTATION_ROUTING.md) is an
unconditional native component. The
[corner-frame/query reduction](depth/BOOLEAN_QUERY_FRAME_REDUCTION.md)
separately transfers lower bounds for arbitrary native Boolean queries
with two calls and the same clean/dirty allocation. These two query
interfaces have different hypotheses. The
[Gaussian inverse extraction](depth/GAUSSIAN_INVERSE_EXTRACTION.md)
prices a candidate repair and proves a tight affine-symplectic round
obstruction for its in-place map. It does not obstruct direct scalar
extraction or unrestricted T-depth. R96 remains the strongest
unconditional upper law, and R79 remains selected Result D.

The [restored-address revision](depth/RESTORED_ADDRESS_SESSION_RESULT.md)
removes gatewise read-only access from R98's sufficient hypothesis.
Joint common-bank queries and paid selected-address packing replace the
shared-address parallel schedule. Every address and helper must return
at the completed query boundary; internal address mixing is allowed.
The conditional bound is unchanged, and no improved arbitrary-table
primitive or unrestricted lower bound is inferred. The
[six-layer source counterexample](depth/READ_ONLY_SIX_LAYER_REDUCTION.md#9-full-return-supplies-no-source-wide-inverse-coverage)
separately rules out inferring source-wide C3 or inverse-Clifford coverage
from full return. The main next question is a faster complete query with
all simultaneous budgets, not another representation identity alone.

The subsequent [joint-transfer round](depth/JOINT_TRANSFER_SESSION_RESULT.md)
first finds the old indicator cost in sparse inverse extraction, then
proves an [unrestricted two-layer query theorem](depth/UNRESTRICTED_TWO_LAYER_QUERY.md).
Shared native one-band encodings force a single occupied transfer alphabet
with one returned clean qubit, without an address-preservation hypothesis.
R75 rounding yields three-control AND error floor $`(3-2\sqrt2)/8`$
and AND arity $`a\ge4`$ floor $`2^{-a}`$. Exact nonaffine queries
must have relative distance at least 1/4 from every affine function.
The same proof home gives an exact rational-target Galois constraint and
native counterexamples to its two-clean and C3 extensions. A six-bit
cubic bent XOR passes the whole scalar alphabet, demonstrating why the
subsequent degree theorem needs a joint invariant. No complete-frame
bound or selected publication claim changes.

The [joint diagonal round](depth/DIAGONAL_ENCODING_SESSION_RESULT.md)
then proves the [unrestricted two-layer Boolean degree bound](depth/TWO_LAYER_DIAGONAL_RIGIDITY.md):
with at most one clean qubit, every exact XOR query has degree at most two,
and every nonquadratic query has uniform initialized-isometry error floor
$`\sin(\pi/8)/(2\sqrt2)`$. Arbitrary dirty width and intermediate
address mixing are allowed. The proof combines full-space diagonal C3
rigidity, explicit clean removal for a non-Pauli image, and native
quotient-row cancellation for a Pauli image. It excludes the six-bit
bent query. A 12-T, two-layer, 44-Clifford one-clean complex diagonal
counterexample still disproves the blanket occupied-C3 extension; its
input-dependent phase cannot receive a free correction in the same
native layers. A paid ten-layer bent-query baseline gives no improved
upper frontier. The next lower question is an arity- and width-independent
Boolean degree bound at three unrestricted layers; R98's faster complete
query remains the parallel constructive target.

The [three-layer native round](depth/THREE_LAYER_NATIVE_SESSION_RESULT.md)
then proves an [unrestricted zero-clean degree bound](depth/THREE_LAYER_NATIVE_STRUCTURE.md):
three layers force XOR degree at most four, with the same uniform tau
error floor for higher-degree predicates and arbitrary intermediate
Clifford mixing. A paid third-band correction realizes literal AND3
using 18 T, 99 Cliffords, one clean and one dirty qubit; depth three is
optimal at this allocation. The actual inverse costs 109 Cliffords.
A fair 16-T/34-Clifford/six-layer baseline prevents mistaking the depth
improvement for count improvement. A growing-arity quartic-sign family
still has an unpaid address phase when used as a larger query.
The [projected-band extension](depth/ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation)
isolates one possible three-dimensional defect. The full one-clean
three-layer degree bound remains open because an intermediate code may
be non-Pauli; the old shared-clean-flip construction is not automatic.
No general complete-frame frontier or publication claim changes.

The [clean-code and algebra round](depth/CLEAN_CODE_ALGEBRA_SESSION_RESULT.md)
then gives three supported components. [Clifford code rigidity](depth/CLEAN_CODE_CLIFFORD_RIGIDITY.md)
classifies literal pointwise stabilizers and localizes the complete native
logical comparison family. Exact examples disprove blanket full-space
Clifford extension and show that arbitrarily many clean-paired axes can
have a one-layer full target. [Quadratic queries](depth/QUADRATIC_NATIVE_QUERY.md)
have optimal depth two with zero clean and one returned dirty helper,
with paid coordinates and actual inverses. A fair six-band cubic baseline
stops the quartic-phase echo as a saving. The
[balanced HRT source](depth/BALANCED_ALGEBRA_PERMUTATION.md) pays in-place
recovery by splitting an algebra generator and returning a dirty
correction bank; six calls occupy five native stages per merge. Its
source and inverse have $`O(\log r)`$ depth and
$`O(3^r)`$ counts and dirty allocation on $`2^r-1`$ arbitrary inputs,
improving the fair $`O(r)`$ depth baseline. Its degree-r returned query
still has growing depth. The next unrestricted lower question is the
one-clean three-layer degree bound after handling unused-sector freedom.
The retained R74 theorem already reduces sufficiently accurate Boolean
approximation to exact occupied action at a uniform small radius;
a faster general R98 query remains open, and R96 is unchanged.

The [occupied-transport revision](depth/OCCUPIED_TRANSPORT_SESSION_RESULT.md)
gives a new [full-space pair transport](depth/OCCUPIED_LOGICAL_TRANSPORT.md).
Its affine Boolean derivatives remove only a quadratic part; paid
endpoint Cliffords put all remaining pair differences into globally
preserved computational parities at the same T-depth and width.
A complete-family counterexample retains non-Clifford unused action at
arbitrarily high clean rank even with zero pairs. The pair-free C3
subcase needs its explicit diagonal silent-band premise.
The [returned algebra query](depth/ALGEBRA_QUERY_EXTRACTION.md) removes
the coefficient-zero promise through eight source calls. Its fast native
subset-product compiler improves the source to
$`O(r^2 2^r\log r)`$ T/G and $`O(r^2 2^r)`$ extra dirty width,
still at $`O(\log r)`$ depth. Neither that table query nor the scoped
arithmetic compression limits change the full-frame frontier.
The [Pauli-boundary theorem](depth/PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md)
closes one substantial subcase: actual Pauli clean images at every band
boundary force sign degree at most 2048, and query degree at most 2047,
even with arbitrary Clifford/address mixing. It retains the exceptional
four-axis relation and bounds the phase dependence of an unbounded
central register through a five-qubit quantum core. The next lower
question concerns non-Pauli intermediate codes, with joint parity-control
phases, silent endpoint bands and unused-sector freedom still present.

The [preceding non-Pauli-boundary checkpoint](depth/NONPAULI_BOUNDARY_SESSION_RESULT.md)
extends that lower closure. The [partial-diagonal theorem](depth/PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md)
proves a bounded-feature description for a full-space native three-band
unitary preserving a computational-Z subgroup of fixed codimension.
Its one-clean reduction gives sign degree at most $`2^{26}`$ and XOR
degree at most $`2^{26}-1`$ whenever either **interior** clean image is
Pauli; the input and final clean stabilizers do not suffice. The retained
both-Pauli bound remains 2048/2047. [Exact repacking](depth/PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#8-both-non-pauli-images-of-active-rank-one-can-be-repacked)
also gives those stronger bounds when both non-Pauli interior images
have active rank one. Thus an unbounded-degree family must have both
interior images non-Pauli and at least one active rank greater than one.
The [one-band non-Pauli decoder](depth/NONPAULI_CODE_DECODING.md)
has an exact controlled-Clifford form, bounded-degree control bits and
paid native forward/inverse words. Its growing displacement-rank witness
prevents using that decoder as a free Pauli-coordinate change within the
same three native bands. Its separate zero-pair lemma assumes active
rank other than two and zero character for every logical probe. The
unused-sector defects then commute, and their phases share one central
Clifford determinant coboundary from the actual middle band. These are
intrinsic central coordinates; the lemma alone gives no uniform target
degree bound or free native implementation of that gauge.

The [bounded-relation theorem](depth/BOUNDED_RELATION_NATIVE_MATCHING.md)
now controls full-space three commuting phase lists with fixed binary
relation dimensions. Its one-clean consequence closes the no-middle-pair
case when both endpoint clean-active ranks are bounded. The proof
retains the resulting dependent middle list and its literal scalar;
it does not compile that list as a free physical native band. The later
[finite-insertion theorem](depth/FINITE_INSERTION_THREE_LAYER_FEATURES.md)
handles an arbitrary bounded Pauli-supported operation that commutes with
each middle axis. Its one-clean reduction closes every fixed pair of
endpoint active ranks, including general middle pairs. For two non-Pauli
interior images it gives sign bound
$`D_{\rm ins}(1,1,1,k_A+k_C-1)`$ and XOR bound one less. Neither
the insertion nor its coordinate changes are supplied for free.
Uniformity in growing active ranks remains open.

The [determinant-defect theorem](depth/DETERMINANT_DEFECT_DEGREE.md)
supplies a different bridge in the zero-pair, all-zero-character branch.
Labels with constant defect phase form a subgroup H. An actual CNOT
correction makes H globally preserved, at unchanged width and T-depth.
Codimension r gives sign degree at most $`D_r`$, with degree five when
r=0. If the entire phase family factors through one common
t-dimensional intrinsic linear quotient, then $`r\le t(t+3)/2`$.
The proof does not identify each intrinsic coordinate with a target
input. The [identity-target family](depth/CLEAN_TRANSITION_PARAMETERS.md#5-unbounded-determinant-codimension-for-literal-occupied-identity)
proves that r and the common central quotient dimension can grow without
bound under these hypotheses, already for a pure determinant coboundary.
Its first interior clean image is Pauli. A separate family in the same
note has both images non-Pauli and degree-three CCZ targets, while its
residual parity-list relation dimension grows. These are representation
limits, not high-degree native query families.

The [rank-one clean-phase interface](depth/RANK_ONE_CLEAN_PHASE_INTERFACE.md)
retains four phase factors in the no-middle-pair case with relation
profile $`(0,\le1,\le1,0)`$. The clean-active factor commutes with
the middle factor and has rank-at-most-one coupling to the silent
endpoint. Raw displacement change can nevertheless grow on an identity
target. The coupled C union is now proved independent, and a shared
clean flip gives sign/XOR bounds five/four at unchanged native depth.
The [computational phase theorem](depth/COMPUTATIONAL_PHASE_DERIVATIVES.md)
uses actual two-band derivatives to localize higher-degree dependence
in the common computational parity part of the two endpoint spans.
Boolean reality forces the active center to transport as Paulis;
a paid Clifford makes it computational. Applied at both endpoints,
this reduces the actual clean transition to one quadratic-controlled
Pauli between three independent native bands, plus a central cubic
sign. A separate noncentral-parity subcase of the actual clean-derived
no-pair circuit has sign/XOR bounds five/four, including endpoint rank two.
Global matching remains open:
conditioning the center can destroy native independence even for an
identity target. An arbitrary quartic correction is not a free gate.

The [joint nonzero-probe theorem](depth/JOINT_NONPAULI_PROBE_ALGEBRA.md)
places the complete generic nonzero-character family and its actual
native factors in one algebra of symplectic rank at most four, assuming
both interior clean images are non-Pauli. This is a common bounded-algebra
component, not a computational phase bound. The genuinely exceptional
non-S probe type at clean active rank two separates from every silent
band. Its occupied action then has an actual zero-clean three-band
logical representation, giving sign degree at most five and XOR degree
at most four. The exceptional hypothesis does not cover every rank-two
code; the finite-insertion theorem covers the cases where both endpoint
ranks are fixed, with a larger bound.

The [symplectic-module readout theorem](depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md)
removes the earlier translation-orbit hypothesis from the source-count
barrier: a degree-d linear matrix coefficient forces a free cyclic module
of dimension $`2^d`$, hence at least $`2^{d-1}`$ carrier qubits and that
many T gates per complete Pauli-conjugating source. It covers arbitrary
modules and full-source approximation at the stated radius, while
excluding signed phase readouts and fused whole-word synthesis. Two
half-address representations escape the single-representation premise
with square-root carrier width and at most Q table Cliffords. Every
source, controlled probe and actual inverse is charged. Its
[universal-probe conversion](depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md#6-universal-pauli-probes-contain-a-complete-dirty-indicator)
uses two complete probe calls to implement an arbitrary-input dirty
indicator, including noncommuting dual labels and address-dependent
source phases. A uniform faster independent-probe source with linear
T/width and quadratic Clifford cost would then supply
[all R98 word-query budgets](depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md#7-the-complete-r98-word-query-interface-follows-conditionally)
for every output length and dirty allocation. That conditional schedule
is proved; the faster source remains unpaid. R91 supplies the unchanged
fair query baseline. The new [Walsh phase attempt](depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md#71-a-dense-walsh-phase-program-does-not-supply-the-faster-source)
pays its formal parity network, but requires generally nonnative
rotations and misses the source count/depth budgets. Its parity-only
repair is excluded; general Clifford interference is not.

The [paid Majorana Walsh source](depth/MAJORANA_WALSH_QUERY.md) replaces
those formal phases by an exact full-input native probe. Its
$`T=2mp`$, $`G=O(m^2)`$, depth at most $`4p`$ and $`m/2`$ carrier
width for $`m=2^p`$ still miss the improved-query budgets. Sharp separate
Walsh-source bounds hold only for Majorana-preserving synthesis.
A flat-source alternative has a paid linear-count word, but retains
depth $`2p`$ and a costly complete dirty dual-frame readout. Its exact
commutator gives a concrete alternate indicator interface; neither
missing resource improvement is supplied. The separate unrestricted
flat-source lower bound is constant, so R96 remains unchanged.

The [returned-query continuation](depth/RETURNED_QUERY_SESSION_RESULT.md)
closes that specific scalar escape: the
[Gaussian scalar-loop theorem](depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
bounds the extracted sign's Boolean degree by $`4r+2`$ for r actual
affine quadratic phase rounds, with fixed Clifford interlayers and even
initialized helpers. Its marker/path-sum proof detects central phases;
the bound is on completed rounds, not unrestricted native T-depth.
The separate [two-layer read-only theorem](depth/READ_ONLY_TWO_LAYER_RIGIDITY.md)
gives robust initialized-work degree and polar-rank restrictions, with
an explicit unrestricted two-layer Toffoli falsifying a broader claim.
The [multiway cost audit](depth/MULTIWAY_QUERY_COST_BARRIER.md) stops
termwise tensor splitting, rank-deficient bucket toggles and cyclic
borrowed-slot substitutes under their respective hypotheses. No improved
unconditional frame bound or growing native depth lower bound is claimed.

The [three-layer continuation](depth/THREE_LAYER_QUERY_SESSION_RESULT.md)
answers the next native query question. Its
[proof home](depth/READ_ONLY_THREE_LAYER_RIGIDITY.md) shows that, below
error $`\sqrt2/6`$, three work-only T layers with preserved address
boundaries force affine XOR outputs with zero clean qubits and degree at
most four with one clean qubit, independently of dirty width. A signed
Pauli orbit controls intermediate clean-code excursions. Exact native
AND fixtures use four zero-clean layers or three one-clean layers; a
two-layer candidate retains a non-Clifford branch phase. The
[exact one-clean two-layer phase argument](depth/READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-one-clean-queries-are-affine)
also forces affineness, making the one-clean three-layer AND depth exactly
optimal in this interface; no robust affine radius is asserted. The next issue
is degree growth with the number of layers, while fixed-error robustness
and arbitrary address mixing remain separate requirements for a general
frame lower bound. The complete-frame upper frontier is unchanged.

The [four-layer follow-up](depth/FOUR_LAYER_QUERY_SESSION_RESULT.md)
closes the first case of that handoff. In
[section 8 of the same proof home](depth/READ_ONLY_THREE_LAYER_RIGIDITY.md#8-four-layers-with-zero-clean-qubits-degree-at-most-two),
zero-clean read-only four-layer XOR queries have degree at most two
below error $`\sqrt2/12`$, independently of dirty width. Clifford
comparison targets and an auxiliary Pauli cocycle handle the vanished
symplectic contribution as well. The exact AND fixture makes this
degree bound sharp. There is no new unrestricted or one-clean four-layer
conclusion.

The [five-layer continuation](depth/FIVE_LAYER_QUERY_SESSION_RESULT.md)
keeps the comparison Clifford and identifies the partial-normalizer gap.
Its [proof home](depth/READ_ONLY_FIVE_LAYER_REDUCTION.md) gives degree at
most seven, or four for a visible shear, under the Clifford-middle condition
at error below $`\sqrt2/12`$. The
[previous follow-up](depth/NATIVE_PARTIAL_ACTION_SESSION_RESULT.md) adds
the Pauli-comparison bound of two without a Clifford middle, independent-axis
matching, and the four-layer scalar CCZ word with 64 T gates, 218 Cliffords
and 17 dirty helpers. These results keep their own hypotheses and ledgers.

The [returned-relations round](depth/RETURNED_RELATIONS_SESSION_RESULT.md)
proved the individual-axis comparison bound of 244 (41 for an invisible
target), the [baseline-word obstruction](depth/BASELINE_WORD_OBSTRUCTION.md),
the [four-factor scalar construction](depth/SCALAR_REPRESENTATION_PRODUCT.md),
and the [robust diagonal-echo bound](depth/RETURNED_DIAGONAL_ECHO_DEGREE.md).
Its proposed span-preserving continuation is now resolved by the
[span-preserving round](depth/SPAN_PRESERVING_SESSION_RESULT.md):
[Section 5.4](depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#54-span-preserving-comparisons-a-finite-degree-bound)
gives degree at most 112 below error $`\sqrt2/24`$, allowing axis mixing,
non-Pauli comparisons and a non-Clifford middle. Intrinsic augmentation
features handle the right action without assuming an ambient normalizer.

The same round proves a source-cost barrier for the scalar construction,
including compressed Pauli carriers and small full-source approximation.
A degree-d table needs at least $`2^{d-1}`$ T gates per Clifford translation
source. Whole-word synthesis remains outside this source-call bound.
A further [shared-prefix invariant](depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#55-shared-prefix-comparisons-a-robust-quadratic-group-map)
gives faithful degree-four actual-group features without span preservation,
below error $`\sqrt2/40`$. With span preservation it sharpens the output
bound to 21 at that smaller radius. The [native-defect round](depth/NATIVE_DEFECT_SESSION_RESULT.md)
resolves the missing target extraction: [Section 5.6](depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
proves the general zero-clean read-only five-layer degree bound **14**
(**11** for an invisible target shear) below $`\sqrt2/24`$, with no extra
comparison hypothesis. The defect's native Pauli pullback permits one-band
matching and bounded-degree signed tableau transport. This supersedes the
span-preserving degree bounds but preserves their structural proofs.
The [six-layer comparison round](depth/SIX_LAYER_TRANSPORT_SESSION_RESULT.md) proves a
[six-layer Pauli-comparison bound](depth/READ_ONLY_SIX_LAYER_REDUCTION.md)
of degree 51 below $`\sqrt2/64`$. The same home contains a fully returned
identity fixture with actual odd-order comparisons and persistent non-Pauli
defects; repeated Pauli-return cannot be the general invariant. The
[diagonal-source proof](depth/RETURNED_DIAGONAL_ECHO_DEGREE.md) now covers q
arbitrary diagonal sources with affine monomial-Clifford transport, giving
degree at most $`\max(2,q-1)`$ below its explicit q-dependent radius.
Its four-call echo radius improves to $`2\sin(\pi/32)`$.
Neither result improves R96 or supplies the missing R98 query.

The [diagonal-source transport proof](depth/NATIVE_CLIFFORD_TRANSPORT.md)
gives six-layer degree at most 78, or 26 for visible outer shear, below
$`\sqrt2/8`$. The [preceding round](depth/TWO_OVERLAP_SESSION_RESULT.md)
proved the [overlap-two source theorem](depth/TWO_OVERLAP_SOURCE_STRUCTURE.md):
an exact three-band C3 source with first/middle overlap at most two
belongs to the diagonal class, even when both adjacent partial sources
fail to be C3. It extends the [overlap-one proof](depth/ONE_OVERLAP_SOURCE_STRUCTURE.md)
using full-space coupled-control and zero-character involution arguments.
The [earlier source theorem](depth/NATIVE_BAND_SOURCE_STRUCTURE.md)
confines middle/outer noncommutation to at most 2t middle axes for general
overlap t. Full return alone does not supply C3 membership.

The new [hybrid transport theorem](depth/HYBRID_DIAGONAL_CORE_TRANSPORT.md)
gives explicit degree bounds through a wide diagonal source coupled to an
a-qubit core, assuming the full source is C3. The bounds depend only on a;
the error radius remains $`\sqrt2/8`$ even with unbounded Pauli deficiency.
At most four fixed blocks recover each output label coordinate; bounded
modulo-four carries recover the literal signs.
A commuting gadget list with r relations following a diagonal-C3 source
has a hybrid core of at most $`4(2^r-1)`$ qubits if their product is C3.
Monomial localization keeps the wide phase diagonal. This does not make
the coordinate changes or dependent gadget list a free native band.

The [preceding round](depth/COUPLED_INCIDENCE_SESSION_RESULT.md) extends
[source certification](depth/ISOLATED_NATIVE_CORE_STRUCTURE.md) to arbitrary
coupled nonzero-character incidence. When every zero-character B axis
commutes with C, the full involution defect and native A isotropy force
all interacting classes to be isolated two-axis pairs. Exact C3/native
hypotheses remain essential; the certificate gives 78/26 transport at
arbitrary overlap. The [zero-character round](depth/ZERO_CHARACTER_SESSION_RESULT.md)
adds [commutation-cover and certified-stripping results](depth/NATIVE_BAND_SOURCE_STRUCTURE.md):
if each B axis commutes with an entire outer band, the source is diagonal
C3, including arbitrary bad-zero-only interactions. Good middle axes can
be removed after certified A-axis stripping, preserving C3. The
[paired-orbit extension](depth/ISOLATED_NATIVE_CORE_STRUCTURE.md) also
proves that zero-character interactions cannot occur only on paired C
orbits when every singleton C axis is good. Both certificates retain
78/26 transport; mixed cases remain.
The [controlled-Pauli transport proof](depth/CONTROLLED_PAULI_TRANSPORT.md)
now derives its canonical exact Pauli domain by signed fixed-Clifford
normalization, removing the former compatibility premise. Every exact
C3 scalar-Pauli block source has uniform 174/58 bounds below
$`\sqrt2/8`$. Its necessary weighted quartic term and both endpoint
coordinate changes are retained. The known Gottesman–Mochon source
shows this class need not be semi-Clifford; its charged source schedule
is not a returned three-band query. The same proof home now characterizes
scalar-Pauli form intrinsically: for exact input domain E and
$`N_p=S_{F\tau(p)F^\dagger}-I`$, it is equivalent to
$`N_pN_q=0`$ for all $`p\in E^\perp`$ and all q. The criterion
keeps signed normalization and is invariant under fixed endpoint
Cliffords. The known four-generator HRT source fails it, separating this
class from generic C3 without asserting a native-three-band counterexample.

The [bounded-deficiency theorem](depth/BOUNDED_DEFICIENCY_TRANSPORT.md)
controls approximate transport through arbitrary sources with bounded
Pauli-domain codimension, uniformly even near Cliffords and for degenerate
domains. Its parameter can grow within a single T band, so it is not a
general depth invariant. The [permutation-source note](depth/PERMUTATION_SOURCE_TRANSPORT.md)
retains another conditional transport theorem; a qualifying three-band
non-C3 returned instance remains unconstructed.

The [commutative character theorem](depth/NONLINEAR_CHARACTER_TRANSPORT.md)
now proves that, for every commutative nilpotent binary coefficient
algebra, a degree-a affine-input character needs at least $`2^{a-1}`$
coefficient wires. This exceeds R98's square-root-width target as a grows. The
[noncommutative cycle mechanism](depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md)
produces degree-d AND on $`d^2`$ coefficient wires with literal full
return. Its native $`O(d^3)`$ counts and $`O(d)`$ depth lose to the
retained conjunction baseline; source and inverse costs remain decisive.

The cycle phase source is identified exactly with the complete
unitriangular-inverse word query after Hadamards on its lower coefficient
register. A fully charged two-call wrapper supplies its control on
arbitrary auxiliary inputs. The new exact commuting-pivot inversion
identity has cubic native counts and linear depth, with returned arbitrary
helpers and no clean work. Whole-pivot commutativity does not eliminate
cross-pivot gate dependencies. This in-place circuit does not provide the
completed XOR query within the budgets needed for a complete-frame improvement.

The balanced block construction now gives a complete structured
inverse-word query and phase source with logarithmic T-depth, cubic
T/Clifford counts, cubic arbitrary work and zero clean work. One parallel
child-inversion position suffices because helper garbage may remain live
until the actual compute/readout/uncompute finishes. Direct controlled
readouts borrow that same work and return it exactly. This leaves the
scoped frozen-coefficient determinant obstruction intact; it does not
supply returned-work in-place inversion. Coefficient inputs change
internally; strengthened R98 now allows this at complete query boundaries.
Paid arbitrary-table coefficient access and extraction within its
count/work/depth budgets remain missing. The full cycle conjunction's
transport and faster fair baseline are unchanged.

The [inverse-query cost audit](depth/GAUSS_SUM_QUERY_ENCODING.md) identifies
the remaining constructive cost: arbitrary-table inverse-entry encodings
already exist, but generic query substitution exceeds the count/work
envelope and dirty loading retains unknown offsets. Affine embeddings
and fixed completed generic-query compositions have scoped capacity
bounds. A full read-only interpolation repair cancels the offset, but
its repeated source evaluations lose to direct read-only path propagation.
Sparse native specialization, coupled masks and nonlinear work remain open.

The remaining component source question is whether native three-band structure forces
the intrinsic tableau product to vanish in the remaining mixed cases.
This would give 174/58 transport without demanding diagonalization.
Alternative native syntheses must be tested before treating a
source parameter as a depth invariant. General non-C3, one-clean and
address-mixing lower bounds and optimal complete-frame depth remain open.
The [query/frame reduction](depth/BOOLEAN_QUERY_FRAME_REDUCTION.md) keeps
the actual inverse Clifford word separately charged.

The [depth preparation map](depth/DEPTH_OBSTACLE_MAP.md) reconstructs the
pre-R78 fixed-accuracy cost ledger, the exact lower-proof bottleneck, and
the assumptions a new route must change. It separates dimension-growing
depth from the next constant-depth obstruction and preserves the scope of
earlier method-specific limits.

The [real R70 theorem](depth/ALL_PRECISION_DEPTH.md) and
[complex R71 theorem](depth/COMPLEX_FRAME_DEPTH.md) retain their two-clean
statements. The two-clean [compressed-column theorem](depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains its original linear-precision statement and source interfaces. The
[older uniform construction](../docs/UNIFORM_PRECISION_DEPTH.md) still
supplies the two-clean precision and width splices.

| Current proof | Role |
|---|---|
| [Precision-uniform subtree depth](depth/UNIFORM_SUBTREE_DEPTH.md) | Selected R79 theorem: retain explicit per-call source precision costs, enlarge the protected bank, and combine R78's decoder with R72's one-clean terminal and all-precision splices; all resource constants are absolute |
| [Dirty catalogue terminal](depth/DIRTY_CATALOGUE_TERMINAL.md) | Current R90 fixed-accuracy general bound: three R89 reductions and a zero-clean catalogue commutator give $`O_L(N/b^2+\ell)`$ depth, with original count budgets, literal widths, and the complete operator contract |
| [Parallel bilinear word query](depth/PARALLEL_BILINEAR_WORD_QUERY.md) | R91 exact component used by R94: private phase targets parallelize arbitrary output words, removing additive q depth with no clean work; independent of R92 and R93 |
| [Uniform catalogue precision](depth/UNIFORM_CATALOGUE_PRECISION.md) | R92 absolute catalogue envelope and a slow unbounded precision window with logarithmic depth; retains R79's other regimes |
| [Grouped catalogue precision](depth/GROUPED_CATALOGUE_PRECISION.md) | R93 growing-precision improvement: local early decoder, whole noncommuting coarse groups, shared matrix-unit residual source, and the original one-clean width/phase splices |
| [Portable workspace catalogue](depth/PORTABLE_WORKSPACE_CATALOGUE.md) | R94 growing-precision bound: carry the entire local R79 workspace through dirty selection, compile larger groups directly to final accuracy, and return borrowed arrays exactly with parallel R91 queries and routing |
| [All-precision linear depth](depth/ALL_PRECISION_LINEAR_DEPTH.md) | R96: $`O(NL/b^2+L+\ell)`$ for all $`n\ge1,L\ge6`$, with optimal counts, one clean bit, original literal widths and full initialized-isometry error |
| [Global borrowed-precision stage](depth/GLOBAL_BORROWED_PRECISION_STAGE.md) | Separated $`O(\log l+\log^3m)`$ reserve, actual global-prefix addresses, explicit standalone costs and full recovery |
| [Fast capped geometric reference](depth/FAST_CAPPED_GEOMETRIC_REFERENCE.md) | Shorter exact geometric cap, full-space comparison to the original reference, complete augmented flags and actual-inverse cancellation |
| [Polynomial-precision bridge](depth/POLYNOMIAL_PRECISION_DEPTH.md) | Independently closed intermediate theorem through $`L\le n^{3/4}`$, using three fast stages and R94's complete-payload catalogue |
| [Direct-prefix high-precision terminal](depth/DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md) | Fixed fine-mask grammar with global-prefix queries and a paid native coarse instruction schedule; original conditional flags, literal controls and complete work return |
| [Logarithmic precision depth](depth/LOGARITHMIC_PRECISION_DEPTH.md) | R95: uniform $`O(NL/b^2+L+\ell)`$ on $`6\le L\le2\ell`$, with optimal counts, one clean bit, original literal widths and full-isometry error |
| [Coarse flags and sparse correction](depth/COARSE_FLAG_SPARSE_CORRECTION.md) | R95 opening: coarse geometric precision, final-accuracy shared-gradient references, a small absolute residual norm, signed Schur sources and complete recovery |
| [Small-workspace descendant source](depth/SMALL_WORKSPACE_DESCENDANT_SOURCE.md) | Logarithmic initialized source scratch via borrowed real calibration; approximate dirty-core return is retained in a whole-physical-space comparison |
| [Local frame depth with relaxed counts](depth/RELAXED_COUNT_LOCAL_DEPTH.md) | One-clean depth $`O(h+\kappa)`$ with fixed-polynomial exponential local count/dirty allowance; exact reference coins, decreasing two-stage construction and complete-payload prefix slots |
| [Causal compact-coin depth](depth/CAUSAL_COMPACT_COIN_DEPTH.md) | Retained R89 theorem and coarse-reduction interface: causal local matchings and parallel fixed-coin reference transforms erase the coarse label in logarithmic depth; its original preserved-prefix recursion has an iterated-logarithmic overhead |
| [Product-coin dilation](depth/COIN_DILATION_DEPTH.md) | Retained R87 theorem and fine-decoder interface: a paid product-coin embedding and tensor Hadamards expose the fine spectral band, with actual-inverse cleanup and the retained all-width splice |
| [Prefix-selected coin terminal](depth/PREFIX_SELECTED_COIN_TERMINAL.md) | Retained R88 real terminal theorem with depth $`O_L(N/b^2+\ell\log_2\log_2(M+4))`$; its sharp-query and literal-width fallback interfaces support R89's stronger prefix recursion |
| [Compact fixed coins](depth/COMPACT_COIN_FLAGS.md) | Exact finite atom partitions, acyclic two-outcome columns, and whole coarse flags; a compact representation alone does not supply an in-place native source |
| [Short-label subset routing](depth/SHORT_LABEL_SUBSET_ROUTING.md) | Finite dense-set matching and staged-shear compression with reversible query words; neither primitive is a prescribed arbitrary-permutation compiler or a frame eraser on its own |
| [Interval-dilation depth](depth/INTERVAL_DILATION_DEPTH.md) | R86 intermediate general bound: positive interval embedding, a bounded absolute-norm Haar conversion, and paid local recovery; its rounding and geometric allocation are retained by R87 |
| [Selected dirty routing](depth/SELECTED_DIRTY_ROUTING.md) | Supporting exact selected wire permutations with arbitrary dirty controls, variable-size query families, and exposed clean-window geometry; [contiguous coin layouts](depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals) have $`O(\ell)`$ routing depth |
| [Subtree-projector decoder](depth/SUBTREE_PROJECTOR_DECODER.md) | R78 decoder proof and fixed-accuracy corollary: replace early amplitude amplification by native transported reflections and coherent spectral-label erasure; retain the one-clean terminal, width splice, literal phases, and optimal count |
| [One-clean count and depth](depth/ONE_CLEAN_DEPTH.md) | Retained R72 predecessor and selected proof interfaces: native one-clean layer words, protected terminal work, explicit width/precision splices, and literal phase composition |
| [All-precision composition](depth/ALL_PRECISION_DEPTH.md) | Retained two-clean theorem: combine complete circuits across precision ranges, retaining literal widths, optimal count, and explicit inherited matching intervals |
| [Complex magnitude-frame depth](depth/COMPLEX_FRAME_DEPTH.md) | Retained two-clean theorem: charge literal phase dressing and shared-workspace composition, with the complex dirty reservation and the same count/depth envelope |
| [Whole-frame binary stream](depth/WHOLE_FRAME_BINARY_DEPTH.md) | Reuse the binary source above a logarithmic protected bank through exponentially growing precision, with explicit precision cardinality in every query charge |
| [High-precision depth schedule](depth/HIGH_PRECISION_DEPTH.md) | Replace completed real masks after the signal lift, sum their rectangular-query depths, and serialize the separately bounded native remainder |
| [Complete sublinear-depth construction](depth/SUBLINEAR_T_DEPTH.md) | Compress all early columns into orthogonal depth/prefix sectors, complete the source error to a full-unitary comparison, recover the prescribed operator, and charge the early/tail composition |
| [Shared binary terminal source](depth/SHARED_BINARY_TERMINAL.md) | Retain one compact geometric source across single-layer residual kernels, record failures coherently, and recover the complete stream within logarithmic conditional workspace |
| [Coherent descendant loader](depth/COHERENT_DESCENDANT_LOADER.md) | Certified signed-state construction, coherent source labels, native calibration, and complete work/error contract |
| [Shared-gradient recovery](depth/SHARED_GRADIENT_RECOVERY.md) | One native phase-source preparation and actual return across early recovery, priced arithmetic and coefficient tables, and a reflection excluding the invariant live gradient |
| [Shallow dirty queries](depth/SHALLOW_DIRTY_QUERIES.md) | Parallel private T layers, arbitrary-dirty Boolean phases, and variable-suffix Hadamard echoes |
| [Nonuniform indicators](depth/NONUNIFORM_DIRTY_INDICATOR.md) | Linear count/work and logarithmic address T-depth, needed to retain the final layers' square-root T-count |

The [joint decoder/terminal session](depth/JOINT_DEPTH_SESSION_RESULT.md)
added three supporting results while keeping R79 as the general
frontier at that checkpoint. At fixed accuracy and square-root width,
[R80](depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md) gives logarithmic
depth for depth-homogeneous magnitudes with arbitrary signs, using one
common product conjugator and exact tag return. Its signed balanced real
case is exact. [R81](depth/RECURSIVE_TERMINAL_DEPTH.md) gives a uniform
terminal-height bound, including a log-times-log-log contribution for
the general logarithmic terminal;
the general early decoder still had square-logarithmic depth then. Their
full contracts, all-width splices, and scope restrictions have separate
proof homes. At the same accuracy and width scale,
[R82](depth/BOUNDED_DENSITY_DEPTH.md) supplies a complete
sub-square-logarithmic circuit for density bounded relative to a positive
product reference. Its weighted correction has constant normalization;
sparse native sources and iterated-logarithmic recursion retain the
one-clean allocation and count budgets. The constants depend on fixed
accuracy and density bounds, so this is also a promised-family result.

The later [branch-orientation session](depth/BRANCH_ORIENTATION_SESSION_RESULT.md)
adds [R83](depth/AFFINE_ORIENTATION_DEPTH.md): logarithmic depth at fixed
accuracy when each depth has a common unordered magnitude pair and its
branch choices are affine functions of the prefix. The leaf map and
signed marker correction are both implemented natively. It also proves
[R84](depth/DETERMINISTIC_FRAME_DEPTH.md), an exact compiler for arbitrary
signed corner rotations with
$`O(\log(n+2)(1+\log^*n))`$ depth at square-root width. Both preserve
one clean qubit, square-root T-count, and linear Clifford count. The
[R82 reference extension](depth/BOUNDED_DENSITY_DEPTH.md#affine-oriented-reference-extension)
uses R83's priced conditional core to admit positive affine-oriented
references with unchanged fixed-density bounds. The subsequent
[in-place permutation checkpoint](depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
stopped without a qualifying general construction. It records a unary
relabeling identity with an unpaid recursive selector, scoped four-shear
and exact encoding restrictions, and a limited corner-word obstruction.
These do not exclude general native circuits. The subsequent
[joint orientation session](depth/JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
executes its [five-hour design](depth/JOINT_ORIENTATION_RECOVERY_PROMPT.md)
and proves [R85](depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md): arbitrary
nonlinear orientations with shared unordered depth magnitudes have
$`O_L(\log(n+2)\log\log(n+4))`$ depth at square-root dirty width,
with one clean qubit and the existing count budgets. A geometric coarse
band selects an unloader that exposes clean suffix space for a local
orientation permutation. Only the coarse tag needs spectral recovery.
Signs, zeros, balanced pairs, and prescribed leaf phases are included.

The [executed general-magnitude design](depth/GENERAL_MAGNITUDE_RECOVERY_PROMPT.md)
and [session outcome](depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
advanced the fixed-accuracy general frontier through R86. A rounded
weighted tree is embedded in consecutive uniform intervals; the
balanced-Haar overlap has absolute norm $`O_L(\sqrt\ell)`$. Its
paid recovery uses the exposed suffix and
[selected dirty routing](depth/SELECTED_DIRTY_ROUTING.md), then erases
the fine label with the local frame compiler and its actual inverse.
The subsequent
[executed factorization design](depth/INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md)
and [session outcome](depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
prove R87 through a different local representation. Independent coin
words encode the branch decisions; tensor Hadamards make the last
nonzero coin block equal the marker depth. The paid source and its
actual inverse erase the fine tag without R86's interval-overlap
amplification. R87 matched R85's depth order for arbitrary magnitudes;
R85 retains its independent promised-family proof.

The subsequent [terminal refinement R88](depth/PREFIX_SELECTED_COIN_TERMINAL.md)
compiles arbitrary deepest-M real-frame layers selected by a preserved
external prefix, including signed and singular tuples. It retains the
one-clean allocation, original real width, full initialized-isometry
contract, and count budgets, with depth
$`O_L(N/b^2+\ell\log_2\log_2(M+4))`$ for every $`0\le M\le n`$.
A logarithmic terminal therefore contributes only
$`O_L(\ell\log_2\log_2(\ell+4))`$. The contiguous-layout routing
specialization has $`O(\ell)`$ depth. At that checkpoint R87's complete
bound was unchanged because coarse recovery still contributed
$`O_L(\ell\log_2(\ell+2))`$.

The [continuation record](depth/COIN_CONTINUATION_SESSION_RESULT.md)
then gives [R89](depth/CAUSAL_COMPACT_COIN_DEPTH.md). One global
compression word and causal local matchings realize the compact-coin
embedding; parallel reference transforms read the coarse band directly.
The last-nonzero test includes each reference source's initialized work.
Preserved-prefix recursion gives the stronger arbitrary-height real
terminal bound $`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$ and the
complete real/complex frontier at that checkpoint. Here $`\log_2^*x`$
counts repeated base-two logarithms until the value is at most 2. Every
source and actual inverse is paid.

The [executed whole-terminal design](depth/TERMINAL_SOURCE_GAUGE_PROMPT.md)
and [session outcome](depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
give [R90](depth/DIRTY_CATALOGUE_TERMINAL.md). Three coarse reductions
replace the unbounded recursive tail by a sufficiently small complete
rounded program. An exact dirty word query and XOR slot translations
select local catalogue factors. Two disjoint borrowed arrays transfer
their commutator to the logical target and return exactly. The selected
factor error, rather than the number of catalogue slots, controls the
complete operator error. The proof includes a zero-clean local synthesis
ledger and the original-width splice.

R79 remains the selected all-precision theorem. The
[executed precision-terminal design](depth/PRECISION_TERMINAL_PROMPT.md)
and [session record](depth/PRECISION_TERMINAL_SESSION_RESULT.md)
give R91–R93. Small complete groups share long-prefix reads while a paid
shared source corrects each actual coarse operator. R93 strictly improves
the inverse-polynomial-accuracy benchmark; R92 independently makes a
slowly growing precision range uniform. R94 uses R79's proved compiler
inside each larger catalogue entry and R91 for parallel selection. Its
complete workspace transfer removes the residual correction and permits
a constant number of groups at final accuracy. R95 improves both its
opening and local depth, reaching logarithmic depth at inverse-polynomial
accuracy while preserving the global count and width contracts.
R96 then extends $`O(NL/b^2+L+\ell)`$ to all precisions through a
global borrowed-work stage and directly programmed terminal. The
[continuing mandate](../WORKSPACE.md#continuing-toward-the-general-answer)
retains optimal general depth and dimension-growing lower bounds as
open questions. Full work return and
simultaneous count/width bounds remain part of any further improvement. The finite
catalogue gives no efficient classical preprocessing, practical crossover,
or scalable native emitter guarantee.

The retained binary-source terminal proof reuses the compact source and shared-source
composition R12/R14, with the inherited geometric weights, coherent
failure history and LCU recovery F1/F10/F11. R20 supplies exact bounded
controls, and the [grouped terminal alternative](depth/GROUPED_TERMINAL_DEPTH.md)
supplies the retained whole-word coarse-prefetch interface. That alternative
keeps its proved $`L^{3/2}`$ contribution and column-forest residuals;
the selected terminal kernel uses sixteen one-layer matrix-unit types.
In the retained compressed-column proof, the early shared phase gradient
and the terminal geometric source have separate paid preparations and return
proofs. No new external premise
or priority claim is introduced.

The [entry-indexed encoding bound](depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
proves sharp worst-case normalization $`\sqrt{n+1}`$ for that
architecture using the balanced Haar frame. At operator error epsilon,
the lower bound is $`(1-\epsilon)\sqrt{n+1}`$. Extra entry-local
channels and reweighting obey the same bound. It does not constrain
unrestricted block encodings or establish a native T-depth lower bound.
R78 replaces that overlap recovery with explicitly priced spectral-label
erasure; R79 makes the resulting resource bound uniform in precision.
The normalization theorem remains unchanged.

The following obstruction branch is independent of the constructive proof
chain. The [resource frontier](../docs/OPEN_PROBLEM.md) gives the exact
witnesses, thresholds, and remaining one-clean question. Zero leaf phases
place the witnesses in the phase-dressed family as well.

| Proof | Established obstruction | Scope that remains open |
|---|---|---|
| [Clean-aware depth obstruction](depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md), R73 | Fixed-clean error floors at three layers, uniform in dimension and dirty width; explicit zero-clean optimum and separate one-clean two-layer gap | Arbitrary Clifford interlayers are allowed; these are constant-depth bounds |
| [Four-layer crossing rank](depth/FOUR_LAYER_CROSSING_RANK.md), R74 | Elementary explicit rigidity gap and necessary middle-rank growth, minimized over all four-layer decompositions | Crossing rank may grow with precision; the qualitative polynomial-limit argument remains separately valid |
| [Four-layer operator rigidity](depth/FOUR_LAYER_OPERATOR_RIGIDITY.md), R75 | Width-independent gap between distinct two-layer Pauli images; zero-clean five-layer bound and exact one-clean terminal-sector preservation below the gap | Terminal-sector preservation alone does not settle one-clean four-layer approximation |
| [One-clean midpoint rigidity](depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md), R76 | Pauli and Hermitian-Clifford midpoint floors using phase discreteness, finite Clifford boundaries, and compressed-commutant lifting | Unrestricted non-Clifford midpoints need an additional rigidity argument |
| [Inner clean spread](depth/ONE_CLEAN_INNER_SPREAD.md), R77 | Bounded-spread non-Clifford extension; both outer layers must grow in every accurate representation; exact example excludes a common sign-compatible sector flip even among arbitrary unitaries | Unbounded inner spread remains eligible; the counterexample does not refute the remaining gap conjecture |
| [Outer-overlap rigidity](depth/OUTER_OVERLAP_RIGIDITY.md), R97 | Positive floor at bounded misalignment relative to either physical outer-axis list; uniform R76 floor at zero misalignment, without bounding inner spread or crossing rank | Both misalignment dimensions must grow; exact same-space descent fails, and displayed center mixing does not exclude another native split |

Together R74 and R77 force all four layers to grow in any hypothetical
arbitrarily accurate one-clean sequence for this witness. R97 additionally
forces both outer-axis misalignment dimensions to diverge in every
admissible native representation. Optimal depth remains open.

The following earlier studies isolate costs or restrict particular
interfaces. They remain as mathematical results and records of the
previous grouped construction; their individual limits do not obstruct
the complete-column source and recovery above.

| Study | What survives | Limit of the result |
|---|---|---|
| [Restricted source depth](depth/SOURCE_T_DEPTH.md) | Exact source schedules and matching depth in a specified Majorana-linear architecture | No unrestricted full-frame depth lower bound |
| [Two-layer source obstruction](depth/SHALLOW_SOURCE_OBSTRUCTION.md) | Full-input source approximation gaps in the stated two-layer model | Does not rule out other source architectures or frame compilers |
| [Source leakage](depth/SOURCE_LEAKAGE_DIAGNOSTICS.md) | Explicit common-precision source families with coherent leakage accumulation | Ideal-angle stability cannot simply be substituted for the full output error; other circuits remain eligible |
| [Flag echoes](depth/HOPF_FLAG_ECHO.md) | Exact errors and a same-mask cancellation exception | Generic radial leakage survives these proposed echo words |
| [Radial filter](depth/HOPF_RADIAL_FILTER.md) | A phase-correct filter with quadratic radial error and a smaller precision cap | Charged extra calls retain the earlier source schedule's depth order |
| [Batched lookup](depth/BATCHED_DIRTY_LOOKUP.md) | A valid intermediate schedule and its distinct ledger | The selected proof uses the later self-contained amortized query; this schedule is retained as a predecessor |
| [Common-source rearrangement](depth/COMMON_SOURCE_REUSE.md) | Exact conjugation and counterexamples to a stale success monitor or one-use fresh bank | The conjugated reflection or another work-return mechanism remains necessary |
| [Protected unary source](depth/PROTECTED_UNARY_SOURCE.md) | One source preparation/actual return across early groups, with a global error charge | Removes source-boundary overhead only |
| [Windowed predicates](depth/WINDOWED_GROUP_PREDICATES.md) | Consumed activity caches with exact inactive cancellation | Predicate overhead becomes sublinear; logical transport and queries remain |
| [Shared-prefix queries](depth/SHARED_PREFIX_QUERY_AUDIT.md) | Exact cache-capacity restrictions, charged refresh, and a fresh-query reduction | Special affine differences help only with an extra promise; generic correction is still charged |
| [Retained-source fusion](depth/RETAINED_SOURCE_FUSION.md) | Exact factors and a two-shift, logarithmic-depth stabilizer completion | The complementary ordered transport remains linear in group height; the original full-group schedule is still retained |

The earlier grouped schedules separately charged O(n) logical transport
and program query/unload. Their local refinements alone did not remove
both allowances. The current construction bypasses those stages by
recovering the complete early operator from its charged column source.

## Attribution and reproducibility

The [research literature notes](RELATED_WORK.md) preserve comparisons for
these studies. Stable source/result identifiers remain in the
[source catalogue](../docs/reference/SOURCE_CATALOGUE.md); the
[finite-evidence catalogue](../docs/reference/VERIFICATION_CATALOGUE.md)
identifies the corresponding checks and their limitations.

Each proof links its executable evidence in the shared package and test
directories. The complete research corpus is included in mathematics,
link and presentation validation.
