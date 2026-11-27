# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for complete
Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

This checkpoint owns the current question and next action. Replace it when
the question changes; proofs and route assessments stay in their linked homes.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

The current continuation started from the verified [R103 merge](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/284)
at [23cfd8b](https://github.com/GoGoKo699/Hopf-Frame-Compilation/commit/23cfd8bfdeefb05a31f87e56d17f3663aa205e5e).
The [returned-frame test](research/depth/ROBUST_DECISION_BIAS.md#complete-return-four-band-floor) raises the one-clean floor to five bands. The [rare cycle](research/depth/ROBUST_DECISION_BIAS.md#zero-clean-six-band-floor) gives six with zero clean work at small fixed error.

<a id="current-depth-frontier"></a>

## Destination and strongest established statements

[R99](research/depth/PROGRAMMED_ROTATION_DEPTH.md) gives
$`O(NL/b^2+\log(n+2)+\log^2(L+2))`$ depth through $`L\le cn`$
with the same simultaneous budgets as [all-precision R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md).
[R100](research/depth/GEOMETRIC_CHUNK_DEPTH.md) gives sublinear precision depth
for $`\log_2 L=o(n)`$. [R102](research/depth/COMPACT_PROGRAM_PRECISION.md) gives $`O(NL/b^2+(1+L/n)\log_2(n+2))`$ through $`L\le2^{cn}`$.
[R104](research/depth/EXACT_KEY_COVER_PRECISION.md) makes the additive precision term sublinear for $`2^{cn}\le L=o(N)`$ at the original widths. The fixed-accuracy gap stays five bands to logarithmic at square-root width.
[R103](research/depth/PRECISION_COUNT_DEPTH_TRADEOFF.md) matches R96 at $`L\ge N`$ with one clean and $`O(L)`$ T gates. [R105](research/depth/CLEAN_SPACE_PRECISION_DEPTH.md) comes within a logarithmic factor of R103 when the clean budget grows through $`\Theta(L^{3/4})`$.
[R79 remains selected Result D](manuscript/PUBLICATION_SCOPE.md).

## Main unfinished route and exact missing operation

[Prescribed XOR embeds in a complete real frame](research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md).
The [strictly positive version](research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md#near-uniform-positive-frames)
retains lookup in an arbitrarily narrow fixed near-uniform density window.
Lookup also suffices for that subclass with one paid balanced frame.

In the constructive direction, [R98](research/depth/QUERY_DEPTH_FRAME_TRANSFER.md)
accepts a complete **word** query. A separate [scalar transfer](research/depth/PRESCRIBED_SCALAR_TABLES.md)
uses prescribed sign tables across widths with an iterated-logarithm cost.
A sufficiently sublogarithmic scalar source improves full-frame depth
at fixed precision. The proof charges its terminal, inverse and work.
The [word-query accuracy result](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
and the [scalar precision wrapper](research/depth/PRESCRIBED_SCALAR_TABLES.md#pair-reflection-scalar-amplification)
have different source hypotheses. A fixed-error scalar source with exact
physical address preservation is one sufficient premise. More generally,
a [full-space signed gap](research/depth/PRESCRIBED_SCALAR_TABLES.md#full-space-signed-gap-source)
permits nonflat spectra and completed address motion, with logarithmic
amplification cost and no additional clean qubits. A
[one-clean coherent decision](research/depth/PRESCRIBED_SCALAR_TABLES.md#one-clean-decision-source)
weakens the source requirement further: amplify on its initialized code
before promoting it to complete return. The raw source may leave garbage,
and its rejected clean sector need not carry the table's prescribed sign.

**Resumption question.** Can distributed counters enforce a common odd
Fourier-key shift in $`o(\log p)`$ native T-depth on arbitrary dirty inputs,
or can a joint word enforce the same orthogonality without that shift?
The [shifted-sector promotion](research/depth/PARALLEL_CLEAN_PROMOTION.md#random-shift-clean-promotion)
pays complete return for a shareable exactly false-silent probe with
$`p^{o(1)}`$ cleans. The [distributed counter](research/depth/COHERENT_COUNTER_SOURCE.md#distributed-cyclic-counter)
pays local visits in $`O((\log\log(p+4))^3)`$ T-depth. Its common shift costs logarithmic depth.
[Low-degree monomial repairs](research/depth/COHERENT_COUNTER_SOURCE.md#monomial-clock-repair-bank) fail in the literal bank at fixed error.

A separate sufficient route gives the correct signed gap on diagonal blocks of
$`M=J^\dagger V^\dagger Z_cVJ`$, with bounded off-diagonal row sum
in each preserved work sector? The labels may include the address
and any dirty inputs. The paid
[quantum records](research/depth/PRESCRIBED_SCALAR_TABLES.md#expander-quantum-address-tags)
then complete the sign with no final clean-width penalty. The
[linear-tag route](research/depth/PRESCRIBED_SCALAR_TABLES.md#summable-address-motion)
uses CNOTs alone if the sum over all address differences is bounded.
Both charge angular error separately from clean leakage.
Pay V and its actual reverse with one clean input, T-count and dirty
width $`O(2^{a/2})`$, and Clifford count $`O(2^a)`$.
Depth $`O(\log a/\log\log a)`$ still suffices at
$`b=K_L\sqrt N`$. A source that proves the full compressed sign directly
can use the [smaller original amplification factor](research/depth/PRESCRIBED_SCALAR_TABLES.md#square-root-width-signed-gap-transfer).
Plain diagonal success cannot support a general short-tag repair.
The [native adversary](research/depth/PRESCRIBED_SCALAR_TABLES.md#adaptive-classical-tag-limit)
defeats short deterministic classical tags even when chosen from its
description. Its sectorwise row sum permits quantum recovery.

The [four-call shear reduction](research/depth/PRESCRIBED_SCALAR_TABLES.md#self-inverse-shear-tables)
still gives a concrete entry readout target for known self-inverse CNOT
shears. Amplify each decision before applying that reduction.

The [separated-block margin bound](research/depth/PRESCRIBED_SCALAR_TABLES.md#two-address-block-margin)
excludes endpoint repairs of the affine Gauss phase. The
[source route register](research/ROUTE_DECISIONS.md#complete-sources-and-promotion)
retains the boundary, additive-visit and paired-Pauli scopes.

The [unrestricted three-band bias theorem](research/depth/ROBUST_DECISION_BIAS.md)
forces exponentially small average signed bias for typical tables. Complete
return and the actual inverse let a sparse Pauli probe extend its sharp
error floor to four bands for a fixed clean count. At five bands a one-clean
compiler needs [large axis intersections](research/depth/ROBUST_DECISION_BIAS.md#complete-return-five-band-intersection)
at both ends. Remove that intersection hypothesis or strengthen the
[prefix and suffix spreading test](research/depth/ROBUST_DECISION_BIAS.md#complete-return-pauli-spreading).
The [four-band decision screens](research/depth/ROBUST_DECISION_BIAS.md#four-band-complementary-spans)
still concern raw decision bias. For
arbitrary depth at square-root T-count, [temporal cut and Gram
codes](research/depth/TEMPORAL_AXIS_BIAS.md) require square-root-scale cut
ranks and dense anticommutation under the stated sign promises. One
[orthogonal component](research/depth/AXIS_TRAJECTORY_BIAS.md#same-component-symplectic-strengthening)
must carry the dependence, shared span and noncommuting rank. For complete
return, [linear traces](research/depth/AXIS_TRAJECTORY_BIAS.md#linear-trace-cut-dependence)
also require dependence within the halves of every chronological cut.
The [raw-decision count bound](research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md#raw-decision-count)
is now tight for strict computational-basis bias: a [paid controlled Gauss word](research/depth/GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions)
represents every table with square-root T-count/width and constant native depth.
Every computational dirty input has the correct exponentially small bias.
Two product superpositions give opposite deterministic answers at every
address, so strict basis signs alone cannot prove growing depth.
An isolated marked address cannot force growing fixed-error depth:
[expander zero reflections](research/depth/EXPANDER_ZERO_REFLECTION.md) already have constant depth.
Exact [two-visit row readouts](research/depth/SCALAR_REPRESENTATION_PRODUCT.md#two-visit-marked-width)
need exponential width for sublinear clean count. Fixed-error dense-row fusion still needs simultaneous error control.

The [feedback-edge pencil](research/depth/PRESCRIBED_SCALAR_TABLES.md#thin-kernel-rank-obstruction)
still needs a changed signal support or baseline. The lower route has to
handle substantial total dependence within the halves of a native cut.

### Other retained results and alternatives

At $`b=\Theta(\sqrt N)`$, [prefix-shared rotations](research/depth/PREFIX_ROTATION_PROFILES.md)
give sublogarithmic depth on complete real frames with $`q=n^{o(1)}`$
arbitrary top levels and one rotation per lower depth within each resulting
subtree. Independent deepest-layer entries still include every Boolean table.
The [four-layer scalar gap](research/depth/FOUR_LAYER_SCALAR_GAP.md)
closes the zero-clean exactly silent single-mark route at seven or more
address bits. It concerns a different allocation from the general one-clean gap.

The [augmented field recoupler](research/depth/AUGMENTED_FIELD_RECOUPLER.md#7-a-logarithmic-depth-decoder-for-the-full-augmented-query)
uses the established logarithmic-depth indicator within the discovery budgets.
Substituting it back into the indicator still pays those query calls.

The complete indicator remains sufficient: on a p-bit address and arbitrary
$`m=2^p`$ bank, toggle the selected bit and return every helper with
zero clean inputs, fixed full-operator error, $`O(m)`$ T gates and dirty
qubits, $`O(m^2)`$ Cliffords and $`o(\log(p+2))`$ native T-depth.
[Sparse-bank preprocessing](research/depth/GATEWISE_READ_ONLY_EQUALITY.md#sparse-code-preprocessing)
cannot upgrade the sparse source: weight $`o(m)`$ needs at least $`m-o(m)`$
clean inputs under any uniform intermediate leakage below one.
[Indicator fusion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) · [Field decisions](research/ROUTE_DECISIONS.md#field-recoupling)

<a id="the-two-useful-source-interfaces"></a>

### Choose a complete source interface

The [interface register](research/ROUTE_DECISIONS.md#complete-interfaces)
retains [complete Pauli derivatives](research/depth/BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)
and [encoded dirty copies](research/depth/ENCODED_DIRTY_COPY_INTERFACE.md).
The latter needs a jointly paid transition and physical copy, but no physical
implementation of the abstract decoder. The [Gram interface](research/ROUTE_DECISIONS.md#gram-interfaces)
needs a normalized complete-input source; use its actual-accuracy
[one-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
or [joint two-clean](research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout)
recovery ledger. Paid wrappers do not construct their source.

Retained alternatives: [four-band front end](research/depth/CLIFFORD_INVOLUTION_SOURCE.md),
[weighted outer phases](research/depth/JOINT_ADDRESS_HELPER_SOURCE.md),
[reduced central gauge](research/depth/COHERENT_CENTRAL_TRANSPORT.md#reduced-carrier-gauge),
[minimal-tree probe](research/depth/ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source),
and [joint precision recovery](research/ROUTE_DECISIONS.md#joint-precision-route)
beyond the separately paid two-path pilots. Consult their route screens before resuming.

## Stopped routes, reopening conditions and superseded instructions

[Route decisions](research/ROUTE_DECISIONS.md) own the exact excluded class
and reopening condition for each mechanism. Proofs, counterexamples and unique
unfinished attempts remain in their linked homes. Historical prompts are
records, not current work orders. A scoped failure does not retire the
unrestricted destination; a finite pilot is not an all-size construction.

<a id="research-decision-and-stopping-rules"></a>
<a id="continuing-toward-the-general-answer"></a>
<a id="read-by-claim"></a>

## Working contract and reading map

Prefer analytic identities and a few exact fixtures. Preserve literal phases,
every logical column, coherent address and arbitrary dirty/reference inputs,
actual native adjoints, and all clean leakage and work return. No QRAM, free
loading or precision states, unpriced clean history, resets or postselection.
A suffix is clean only on its active sector; inactive action needs its own
proof. Separate general theorems, conditional interfaces, components and
representation limits. Internal audits and finite checks do not replace
analytic proof or external technical review.

| Information that changed | Update its authoritative home |
| --- | --- |
| Construction or proof | Existing proof home; [research catalogue](research/README.md) locates it. Change the index only for a new home or changed role. |
| Mechanism assessment or reopening condition | [Route decisions](research/ROUTE_DECISIONS.md), linked to the proof. |
| Bottleneck or next action | Replace this checkpoint; keep detailed outcomes in existing proof/route homes. |
| General theorem or global gap | [Resource frontier](docs/OPEN_PROBLEM.md). Selected A–D statements belong in [publication scope](manuscript/PUBLICATION_SCOPE.md); [proof map](docs/README.md) owns dependencies. |
| Evidence or attribution | Relevant tests and [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md); [source catalogue](docs/reference/SOURCE_CATALOGUE.md) owns attribution and stable IDs. |
| Validation procedure or upstream baseline | [Research validation](scripts/README.md#research-validation) or [SYNC](SYNC.md), respectively. |

<a id="refurnishing-policy"></a>

Update only affected homes. Keep summaries as links plus decision-relevant
scope. Preserve unique attempts, stable identifiers, inbound anchors, tests
and receipt hashes; change mechanisms when complete synthesis restores the old cost.

## Restore and verify

Check worktree status and remote, fetch main, and inspect changes since the
scientific baseline. Preserve unrelated work. Read the checkpoint and the
dependency branch needed for the next action; use catalogues for targeted
lookup rather than replaying session history.

[Research validation](scripts/README.md#research-validation) owns focused
iteration checks and the complete ten-command final gate. Complete focused
[visual checks](assets/README.md#rendering-checks) before the memory-heavy
numerical suite; full presentation CI covers the corpus. For proof moves,
compare mathematics and audit old and new links. Record exact-tree receipts
in the review record and merge after required checks pass.
