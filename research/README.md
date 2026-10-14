# Research results for complete-frame compilation

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Endpoint status](../docs/OPEN_PROBLEM.md) · [Core proof map](../docs/README.md)

The growing-source proofs below establish Result C's one-clean count
frontier for complete real and phase-dressed complex magnitude frames.
The [core proof map](../docs/README.md#c-one-clean-compilation-and-its-corollaries)
gives the minimal dependency chain; the remaining studies retain their
own hypotheses and distinct resource points. The
[sublinear-depth construction](depth/SUBLINEAR_T_DEPTH.md) resolves
$`o(n)`$ T-depth existence at fixed accuracy and sufficiently large
constant square-root dirty width. Its precision-uniform extension retains
sublinear depth with shrinking error $`2^{-L}`$ for
$`L=o(n/\log(n+2))`$ at sufficient $`b=\Theta(\sqrt{NL})`$.
The [one-clean depth theorem](depth/ONE_CLEAN_DEPTH.md) combines
that branch with a whole-frame binary stream and a high-precision query
schedule, retaining optimal count for every precision for both frame
families. Optimal T-depth remains open.
Each result retains its mathematical hypotheses; a method-specific
obstruction is not a lower bound for every compiler.

The [state-based QBP supplement](../supplements/state_based_qbp/README.md)
uses a different decoder and has its own task-level guarantee.

## Endpoint results

<a id="precision-uniform-two-clean-count"></a>

### Precision-uniform one-clean count

The [one-clean growing-source theorem](endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives $`T=O(N+L)`$ and $`G=O(NL)`$ with universal constants, for
$`n\ge1`$, integer $`L\ge6`$, and error $`2^{-L}`$. It uses one clean flag and
$`b\ge L+n+7`$ dirty qubits for real frames, or $`b\ge L+n+8`$
for phase-dressed complex magnitude frames. The complete-source conjugate
lift borrows the phase signal and retains the initialized anchor marker.
The [two-clean theorem](endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
supplies the precision/mode separation. The one-clean layerwise baseline
and growing-source interfaces meet at $`L=N/n^2`$; the same baseline
covers bounded n, including $`n=1,2`$.
Every logical input, dirty/reference correlation, literal phase, clean-work
leakage, and actual inverse is included in the initialized-isometry guarantee.
At exactly $`b=L+n+7`$ for real frames or $`b=L+n+8`$ for complex
magnitude frames, and $`L\ge\max\{6,n\}`$, the existing finite-width
lower bound gives worst-case $`T^\star=\Theta(N+L)`$. The upper bound
allows additional dirty work; this matching corollary concerns the displayed
base widths.

### Optimal banked count frontier

The [banked growing-source theorem](endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
gives, with one clean flag and error $`2^{-L}`$,

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O(NL),
```

for all $`n\ge1`$ and integer $`L\ge6`$, at real
$`b\ge2(L+n+7)`$ or complex magnitude $`b\ge2(L+n+8)`$.
The banked source masks return every bank and selector before the next
operation. Summed native charges and the existing width-sensitive lower
bounds establish matching throughout these banked regimes.

### Constant-clean high-precision count

The target is the complete prescribed frame at the literal allocation in
the [endpoint statement](../docs/OPEN_PROBLEM.md#constant-clean-high-precision-count).
The [growing-source compiler](endpoint/GROWING_SOURCE_COMPILER.md) attains
$`T=\Theta(N)`$ in the worst case at $`\eta=2^{-N}`$, $`n\ge3`$,
with two clean flags and $`N+n+7`$ dirty qubits. Its upper bound covers
every prescribed real frame. Growing sources, sealed rejection modes,
frozen output tags, actual dirty replay, and final recovery are charged.
The [complex-magnitude corollary](endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
composes the half-error real construction with the packed literal diagonal,
covering $`D_\phi W_{\mathbb R}`$ at the same allocation and accuracy.
The real witness at zero leaf phases gives the same $`2N-1`$ lower bound.
This sharper complex reservation remains a two-clean result; the one-clean
complex theorem uses $`N+n+8`$ dirty qubits at $`L=N`$. None of the
selected count theorems claims arbitrary-complex-unitary synthesis, a new
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

The [one-clean theorem](depth/ONE_CLEAN_DEPTH.md), R72, gives, for every
$`n\ge1`$, integer $`L\ge6`$, one external clean flag and every
$`b\ge17(L+n+7)`$,

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\},
```

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

All constants are absolute. One circuit preserves the prescribed complete
real frame with error at most $`2^{-L}`$ in the full initialized-isometry
norm. Count is optimal at every precision; count and depth match through
$`b\le\sqrt{NL/F(n,L)}`$ on every nonempty interval above the
literal lower threshold. For $`L>N/n^2`$, the stronger depth has
additive L in place of F.

R72 proves the same bounds for the prescribed magnitude frame
$`D_\phi W_{\mathbb R}`$, with one clean flag and
$`b\ge17(L+n+8)`$. It charges the native one-clean literal diagonal
and composes half-error factors on the same workspace, including earlier
leakage and common phase. Its depth matching interval has the corresponding
complex lower threshold. Optimal large-width depth remains open. The
[real R70 theorem](depth/ALL_PRECISION_DEPTH.md) and
[complex R71 theorem](depth/COMPLEX_FRAME_DEPTH.md) retain their two-clean
statements and supporting proofs.

For real frames with two clean flags, the retained [compressed-column theorem](depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains the stronger logarithm $`\log_2(n+2)`$ for $`6\le L\le cn`$,
with absolute c. The [older uniform construction](../docs/UNIFORM_PRECISION_DEPTH.md)
retains its role in the two-clean precision and width splices.

For every fixed $`0\lt\eta\le1/64`$, R72 gives the one-clean all-n bounds $`T=O_\eta(\sqrt N+N/b)`$, $`G=O_\eta(N)`$, and

```math
D_T=O_\eta(N/b^2+\min\{n,\sqrt n\log_2(n+2)\}).
```

At sufficient square-root dirty width it gives $`o(n)`$ depth.
The uniform extension also gives sublinear depth for shrinking error
when $`L=o(n/\log(n+2))`$ and $`b=\Theta(\sqrt{NL})`$ is
sufficiently large. For $`L=O(\sqrt n)`$ this retains the
$`O(\sqrt n\log(n+2))`$ fixed-accuracy depth order. Optimal depth
remains open.

| Current proof | Role |
|---|---|
| [One-clean count and depth](depth/ONE_CLEAN_DEPTH.md) | Retain both families' simultaneous bounds using native one-clean layer words, protected terminal work, explicit width/precision splices, and literal phase composition |
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

The terminal proof reuses the compact binary source and shared-source
composition R12/R14, with the inherited geometric weights, coherent
failure history and LCU recovery F1/F10/F11. R20 supplies exact bounded
controls, and the [grouped terminal alternative](depth/GROUPED_TERMINAL_DEPTH.md)
supplies the retained whole-word coarse-prefetch interface. That alternative
keeps its proved $`L^{3/2}`$ contribution and column-forest residuals;
the selected terminal kernel uses sixteen one-layer matrix-unit types.
The early shared phase gradient and the terminal geometric source have
separate paid preparations and return proofs. No new external premise
or priority claim is introduced.

The [entry-indexed encoding bound](depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
proves sharp worst-case normalization $`\sqrt{n+1}`$ for that
architecture using the balanced Haar frame. At operator error epsilon,
the lower bound is $`(1-\epsilon)\sqrt{n+1}`$. Extra entry-local
channels and reweighting obey the same bound. It does not constrain
unrestricted block encodings or establish a native T-depth lower bound.

The [clean-aware depth obstruction](depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md),
R73, instead applies to arbitrary coherent Clifford+T circuits. For every
fixed clean budget $`a\ge0`$, a complete real frame requires at least four
T layers below some positive accuracy threshold $`\varepsilon_a`$,
uniformly over dimension and dirty width. The threshold is not evaluated.
An explicit one-clean witness separately requires at least three layers
for integer $`L\ge8`$. Both witnesses also lie in the phase-dressed
family with zero leaf phases. These are constant lower bounds, not bounds
growing with dimension or precision. Optimal depth remains an open problem.

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
