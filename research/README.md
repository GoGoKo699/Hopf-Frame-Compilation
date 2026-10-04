# Research archive: the two remaining compiler gaps

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Current gaps](../docs/OPEN_PROBLEM.md) · [Core proof map](../docs/README.md)

This archive preserves attempts to improve the constant-clean count endpoint
and large-workspace T-depth. These notes are outside the proof chain selected
for Results A–D. Some contain proved special cases or useful components;
others identify why a specific proposal fails. Their status is stated below.
An obstruction to one interface is not a lower bound for every compiler.

The scientific scope is frozen after the
[bounded internal audit](../docs/CORE_CLAIM_AUDIT.md). Nothing in an archived
“next task” paragraph automatically reopens research. The
[current gap statement](../docs/OPEN_PROBLEM.md) gives the live conditions.
The separately completed [state-based QBP supplement](../supplements/state_based_qbp/README.md)
is a different decoder and is not failed endpoint work.

## What we tried

### Constant-clean high-precision count

The target is the complete prescribed frame at the literal allocation in
the [endpoint statement](../docs/OPEN_PROBLEM.md#constant-clean-high-precision-count).
Compact classical descriptions, first-column agreement, or uncharged
program access do not meet that target.

| Study | What survives | Why the general endpoint remains open |
|---|---|---|
| [Source reuse](endpoint/SOURCE_REUSE_LIMITS.md) | Encoded-source restrictions, exact source-width transitions, and small source-merging identities | Transformed masks, correlations, and rejected sectors still require charged work; a reusable precision state alone supplies no joint compiler |
| [Tree transport](endpoint/ENDPOINT_TREE_TRANSPORT.md) | Sparse generators, weighted norms, complete transport columns, and native small-mode constructions | Representation size and classical conditioning do not price the joint native operation; unchanged local scattering retains its query limitation |
| [Weighted transport block](endpoint/WEIGHTED_TRANSPORT_BLOCK.md) | Valid complete dilation with exact/approximate dirty-return variants | The constructed native word retains the repeated precision cost |
| [Residual assembly](endpoint/RESIDUAL_ASSEMBLY.md) | Complete forward/reverse assembly, actual inverses, and coupled/commutator repairs | The proved accounting retains the per-depth precision term |
| [Canonical grouped scalar](endpoint/CANONICAL_SCALAR_COMPLETION.md) | A canonical scalar and direction-controlled inverse fit the existing group interface | Consolidating the scalar calls also helps the old construction; precision is still paid per group |
| [Antichain changes](endpoint/ANTICHAIN_COMPILER.md) | **Proved special case:** linear endpoint count for a literal native-baseline antichain promise | Generic rounded frames need not satisfy that promise |
| [Sparse nested changes](endpoint/SPARSE_UPDATE_COMPILER.md) | **Proved special case:** linear endpoint count under the stated ancestor-closure condition, including a changed path | This is a restricted update family, not the unrestricted frame family |

The [retained route assessments](ROUTE_HISTORY.md) explain the earlier
selection decisions and sufficient recurrences. They are historical
assessments; the current rule is to require a new priced native mechanism
before another fixture pass.

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

Code and tests remain in their shared package and test directories, where
their imports and full validation still work. Each archived proof retains
its links to that evidence. Moving a note does not downgrade or remove its
checks, and the archive is not excluded from mathematics or link validation.

No attempt is resumed merely because another finite example can be run.
Reopen for a concrete qualifying rule, an unrestricted lower-bound idea,
a discovered defect, or an explicit change of scope. Manuscript writing
remains on hold.
