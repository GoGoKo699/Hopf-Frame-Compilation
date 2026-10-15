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
| [Disjoint native blocks](endpoint/DISJOINT_BLOCK_COMPILER.md) | Linear endpoint count for disjoint SU blocks; unbanked and banked high-precision block/state interfaces with full Clifford ledgers | High-precision forms require $`L>N/n^2`$ beyond their fixed dimension threshold; each interface retains its own clean allocation, literal phase, and complete return |
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
depth remains open.

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
adds three supporting results while keeping R79 as the selected general
theorem. At fixed accuracy and square-root width,
[R80](depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md) gives logarithmic
depth for depth-homogeneous magnitudes with arbitrary signs, using one
common product conjugator and exact tag return. Its signed balanced real
case is exact. [R81](depth/RECURSIVE_TERMINAL_DEPTH.md) recursively improves
the general logarithmic terminal to a log-times-log-log contribution;
the general early decoder still has square-logarithmic depth. Their
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
references with unchanged fixed-density bounds. The branch report is
the current constructive stopping point: arbitrary nonlinear leaf and
marker permutations remain unpaid, and R79 remains the general theorem.

The terminal proof reuses the compact binary source and shared-source
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

Together R74 and R77 force all four layers to grow in any hypothetical
arbitrarily accurate one-clean sequence for this witness. Optimal depth
remains open.

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
