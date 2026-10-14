# Research results for complete-frame compilation

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Current gaps](../docs/OPEN_PROBLEM.md) · [Core proof map](../docs/README.md)

These notes give constructive reductions, promised-family compilers, and
limits of particular interfaces for the constant-clean count and T-depth
questions. They are separate from the proof chain for Results A–D.
Each result retains its mathematical hypotheses; a method-specific
obstruction is not a lower bound for every compiler.

The [state-based QBP supplement](../supplements/state_based_qbp/README.md)
uses a different decoder and has its own task-level guarantee.

## Endpoint results

### Constant-clean high-precision count

The target is the complete prescribed frame at the literal allocation in
the [endpoint statement](../docs/OPEN_PROBLEM.md#constant-clean-high-precision-count).
The [growing-source compiler](endpoint/GROWING_SOURCE_COMPILER.md) attains
$`T=O(N)`$ for every real frame with two clean flags and $`N+n+7`$
dirty qubits, matching the endpoint lower bound in order.
Compact classical descriptions, first-column agreement, or uncharged
program access do not meet that target.

#### Native primitives and promised families

| Study | Result | Scope |
|---|---|---|
| [Growing precision source](endpoint/GROWING_SOURCE_COMPILER.md) | Unrestricted real-frame endpoint $`T=\Theta(N)`$ at error $`2^{-N}`$, two clean flags and $`N+n+7`$ dirty qubits | Growing sources, sealed rejection modes, frozen output tags and actual dirty replay are charged; no additional one-clean, complex-frame or depth bound |
| [Disjoint native blocks](endpoint/DISJOINT_BLOCK_COMPILER.md) | Linear-count synthesis of disjoint SU blocks with a linear sum of squared dimensions; disjoint reflection corollary | Actual coarse mesh, fine sources, literal phase, and complete return are charged |
| [Rank-weighted boundary compiler](endpoint/RANK_WEIGHTED_BOUNDARY_COMPILER.md) | Linear native count under a support-times-rank budget; boundary-column and collective common-basis fusion | Two-clean complex and one-clean real domains; growing sequential overlap does not inherit the disjoint budget |
| [Borrowed phase and recovery sectors](endpoint/CONJUGATE_SIGNAL_RECOVERY.md) | Full conjugate-sector phase word, zero-clean real reflection banks, and one-clean real block refinements | The mode/recovery flag remains necessary for the stated accepted-block construction; exact native counterexample and balanced-sector bound |
| [Borrowed selected axes](endpoint/BORROWED_AXIS_COMPILER.md) | Full-input selected Pauli/source-axis rotations with two arbitrary helpers; jointly compiled raw rejection orientations at O(Q+q+ell) T-count | Literal query phases and shared-core width are charged; raw canonicalization does not compile transported logical prefixes |
| [Matching and leaf deflation](endpoint/MATCHING_LEAF_DEFLATION.md) | Canonical rank equals twice forest matching size; matching-budget compilers and an explicitly closed nonidentity deflation | Paid input correction retained; separate-level precision still costs O(nN); no unrestricted endpoint improvement |
| [Precision-gauge compiler](endpoint/PRECISION_GAUGE_COMPILER.md) | Linear endpoint count on regular charts, fixed-ratio exponent intervals, and the stated cheap/ultraweak family; exact mixed-height scheduling bound | Rescheduling unchanged merge banks retains repeated precision cost; the next construction changes those sources |
| [Heterogeneous precision source](endpoint/HETEROGENEOUS_PRECISION_SOURCE.md) | At every integer height, $`O(N+\mathcal B_K)`$ count after certified clipping; linear retained bottom trace gives linear count | Fixed common-source construction with mixed local gauges and actual dirty replay; includes the earlier exponent-envelope and mixed-height families |
| [Structured forests](endpoint/STRUCTURED_FOREST_COMPILERS.md) | Linear endpoint count under sub-binary degree-window growth, including fixed computational outdegree; over a cheap or fixed-ratio base, either bounded exceptional path depth or bounded exits per base component suffices | Bounded exits allow arbitrarily deep exceptional paths; original computational-label tree and fixed promise constants, with exact signed and corner tuples included |
| [Antichain changes](endpoint/ANTICHAIN_COMPILER.md) | **Proved special case:** linear endpoint count for a literal native-baseline antichain promise | Generic rounded frames need not satisfy that promise |
| [Sparse nested changes](endpoint/SPARSE_UPDATE_COMPILER.md) | **Proved special case:** linear endpoint count under the stated ancestor-closure condition, including a changed path | This is a restricted update family, not the unrestricted frame family |

#### Full-frame reductions and remaining native steps

| Study | Result | Scope |
|---|---|---|
| [Fixed-tree Cayley reduction](endpoint/TREE_CAYLEY_REDUCTION.md) | All-angle phase gauge, regular parameter chart, exact O(N)-T permutations, and a two-diagonal-algebra generator | The core has an O(nN)-T fine realization; O(N)-T synthesis remains open |
| [Boundary propagation](endpoint/BOUNDARY_PROPAGATION.md) | Sparse preconditioning, repeated-call precision within the exact width, and complete feedback cancellation to one prefix encoder | Repeated calls cost O(RN); the surviving fine encoder has a direct O(nN)-T certificate |
| [Coarse encoder](endpoint/COARSE_PREFIX_ENCODER.md) and [collective refinement](endpoint/COLLECTIVE_PRECISION_REFINEMENT.md) | Joint prefix compilation, exact geometric history, quadratic physical purification, and a cubic midpoint replacement at O(N) T-count | Fixed-order refinement of the terminal regular-core frame; an all-order linear-cost ledger remains open |
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

The selected theorem gives linear-in-n depth at fixed accuracy and
sufficient square-root dirty width. The following refinements isolate
particular costs; they do not remove both remaining linear allowances.

| Study | What survives | Limit of the result |
|---|---|---|
| [Restricted source depth](depth/SOURCE_T_DEPTH.md) | Exact source schedules and matching depth in a specified Majorana-linear architecture | No unrestricted full-frame depth lower bound |
| [Two-layer source obstruction](depth/SHALLOW_SOURCE_OBSTRUCTION.md) | Full-input source approximation gaps in the stated two-layer model | Does not rule out other source architectures or frame compilers |
| [Source leakage](depth/SOURCE_LEAKAGE_DIAGNOSTICS.md) | Explicit common-precision source families with coherent leakage accumulation | Ideal-angle stability cannot simply be substituted for the full output error; other circuits remain eligible |
| [Flag echoes](depth/HOPF_FLAG_ECHO.md) | Exact errors and a same-mask cancellation exception | Generic radial leakage survives these proposed echo words |
| [Radial filter](depth/HOPF_RADIAL_FILTER.md) | A phase-correct filter with quadratic radial error and a smaller precision cap | Charged extra calls leave the asymptotic complete-frame frontier unchanged |
| [Batched lookup](depth/BATCHED_DIRTY_LOOKUP.md) | A valid intermediate schedule and its distinct ledger | The selected proof uses the later self-contained amortized query; this schedule is retained as a predecessor |
| [Common-source rearrangement](depth/COMMON_SOURCE_REUSE.md) | Exact conjugation and counterexamples to a stale success monitor or one-use fresh bank | The conjugated reflection or another work-return mechanism remains necessary |
| [Nonuniform indicators](depth/NONUNIFORM_DIRTY_INDICATOR.md) | Linear count/work and logarithmic address T-depth; a sublinear eligible late tail | Early logical stages and program queries still cost linear depth |
| [Protected unary source](depth/PROTECTED_UNARY_SOURCE.md) | One source preparation/actual return across early groups, with a global error charge | Removes source-boundary overhead only |
| [Windowed predicates](depth/WINDOWED_GROUP_PREDICATES.md) | Consumed activity caches with exact inactive cancellation | Predicate overhead becomes sublinear; logical transport and queries remain |
| [Shared-prefix queries](depth/SHARED_PREFIX_QUERY_AUDIT.md) | Exact cache-capacity restrictions, charged refresh, and a fresh-query reduction | Special affine differences help only with an extra promise; generic correction is still charged |
| [Retained-source fusion](depth/RETAINED_SOURCE_FUSION.md) | Exact factors and a two-shift, logarithmic-depth stabilizer completion | The complementary ordered transport remains linear in group height; the original full-group schedule is still retained |

At fixed accuracy and sufficient width, logical transport and program
query/unload each retain an O(n) allowance. Improving only one is a
component milestone. A complete-frame improvement must reduce both or
bypass them with one proved circuit, with work and error charged.

## Attribution and reproducibility

The [research literature notes](RELATED_WORK.md) preserve comparisons for
these studies. Stable source/result identifiers remain in the
[source catalogue](../docs/reference/SOURCE_CATALOGUE.md); the
[finite-evidence catalogue](../docs/reference/VERIFICATION_CATALOGUE.md)
identifies the corresponding checks and their limitations.

Each proof links its executable evidence in the shared package and test
directories. The complete research corpus is included in mathematics,
link and presentation validation.
