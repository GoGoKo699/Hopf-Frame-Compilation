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
Compact classical descriptions, first-column agreement, or uncharged
program access do not meet that target.

| Study | Result | Scope |
|---|---|---|
| [Bounded diagonal factorization](endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md) | Certified factor search at the literal dirty budget and a stable nested-projector criterion for complete-column coverage | Bounded Haar-mask words are excluded; coverage by a different mixer family is a separate premise |
| [Fixed-tree Cayley reduction](endpoint/TREE_CAYLEY_REDUCTION.md) | All-angle phase gauge, regular parameter chart, exact O(N)-T permutations, and a two-diagonal-algebra generator | The remaining Cayley operator needs a charged native implementation |
| [Boundary propagation](endpoint/BOUNDARY_PROPAGATION.md) | Sparse preconditioning, repeated-call precision within the exact width, and complete feedback cancellation to one prefix encoder | Repeated calls cost O(RN); fine encoder synthesis remains the surviving operation |
| [Coarse prefix encoder](endpoint/COARSE_PREFIX_ENCODER.md) | One O(N)-T circuit for all prefix frames at error $`2^{-\lceil N/n\rceil}`$, with full-operator SU(2) return and small-conjugation error suppression | Coarse precision and one quadratic correction do not establish the endpoint |
| [Structural compilation limits](endpoint/STRUCTURAL_COMPILATION_LIMITS.md) | Sharp normalization; arbitrary-order Haar and fixed-alphabet restrictions; reflection-rank, diagonal-displacement and original-angle query bounds | Each theorem retains its representation assumptions; the stronger alternating-Haar bound also remains valid |
| [Source reuse](endpoint/SOURCE_REUSE_LIMITS.md) | Exact transformed-mask fusion, encoded-source and Majorana extraction restrictions, Spin covariance, single-Clifford spectral obstruction, and exact dirty-work advantage | Programmed source transformations retain their charge; the restrictions do not apply to every native compiler |
| [Tree transport](endpoint/ENDPOINT_TREE_TRANSPORT.md) | Sparse generators, weighted norms, complete transport columns, and native small-mode constructions | Representation size and classical conditioning do not price the joint native operation; unchanged local scattering retains its query limitation |
| [Weighted transport block](endpoint/WEIGHTED_TRANSPORT_BLOCK.md) | Valid complete dilation with exact/approximate dirty-return variants | The constructed native word retains the repeated precision cost |
| [Residual assembly](endpoint/RESIDUAL_ASSEMBLY.md) | Forward/reverse assembly, coupled repairs, branch-rank bound, independent-forest packing theorem and polar leakage alignment | The $`\Omega(Nn)`$ forest bound concerns independent coefficients; actual Hopf path products are correlated |
| [Canonical grouped scalar](endpoint/CANONICAL_SCALAR_COMPLETION.md) | A canonical scalar and direction-controlled inverse fit the existing group interface | Consolidating the scalar calls also helps the old construction; precision is still paid per group |
| [Antichain changes](endpoint/ANTICHAIN_COMPILER.md) | **Proved special case:** linear endpoint count for a literal native-baseline antichain promise | Generic rounded frames need not satisfy that promise |
| [Sparse nested changes](endpoint/SPARSE_UPDATE_COMPILER.md) | **Proved special case:** linear endpoint count under the stated ancestor-closure condition, including a changed path | This is a restricted update family, not the unrestricted frame family |

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
