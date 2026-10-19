# Continuing workspace

The repository is organized around the
[claims selected for the manuscript](manuscript/PUBLICATION_SCOPE.md).
The [resource statement](docs/OPEN_PROBLEM.md) records the proved
one-clean count theorem in Result C, its optimal banked count, the
sharper two-clean complex endpoint, all-precision count/depth bounds,
precision-uniform subtree depth, and the fixed-accuracy general
logarithmic-depth theorem R90, the precision refinements R92–R96, and
R97's outer-overlap obstruction. The
[research index](research/README.md) groups native primitives and promised families, full-frame reductions, source programs,
and representation limits. The [evidence catalogue](docs/reference/VERIFICATION_CATALOGUE.md)
links their implementations and checks. Each result has one proof home;
a completed source isometry or compact representation does not by itself
supply a linear-count native recovery circuit.
Final manuscript drafting is a separate task.

## Current checkpoint: joint indicator mechanisms and native orbit closure

The [current revision](research/depth/JOINT_INDICATOR_SESSION_RESULT.md)
began at verified main `b598e4aa55c308f11f183a11ab92f0e4fc70463a`
(PR #195). R96 remains the strongest unconditional general upper theorem:
$`D_T=O(NL/b^2+L+\log_2(n+2))`$, optimal T-count, $`G=O(NL)`$,
one clean qubit and the original literal widths. R79 remains selected
Result D. Optimal general depth remains the destination.

The main constructive target is the complete approximate dirty indicator.
At every fixed requested tolerance it needs linear-in-bank-size T-count
and dirty width, quadratic Clifford count, and sublogarithmic address
T-depth in full operator norm. The current restored-address R98 and its
[constant-accuracy extension](research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
provide the paid transfer; temporary address changes are allowed inside
completed joint queries.

The [maximal source-gauge theorem](research/depth/WEAKER_INDICATOR_INTERFACES.md)
allows arbitrary bank-controlled address and work motion, with neither
returned at the source boundary. A two-call Clifford probe still extracts
the complete indicator. [Selected routing](research/depth/JOINT_ROUTING_INDICATOR_INTERFACES.md)
is likewise equivalent under a paid constant-call reduction. Its direct
switch tree keeps a sequential address cost; fixed Clifford re-encoding
cannot simultaneously diagonalize all bank translations. Neither result
is an unrestricted depth obstruction.

A [noncommuting dirty phase loop](research/depth/DIRTY_GEOMETRIC_PHASE_LOOPS.md)
uses a geometric scalar on both signal eigenstates, avoiding the earlier
rank-one-reflection fixed-space premise. Its full-bank accuracy and native
digital recurrence retain costs outside the target envelope. A collective
native word that avoids independent-row precision accumulation remains
an open constructive question.

The parallel [native-orbit theorem](research/depth/GROWING_DISPLACEMENT_NATIVE_ORBITS.md)
closes the growing-total-displacement case left by PR #195 in the actual
no-middle-pair reduction. The full cubic tensor forces actual pair
permutations plus one exceptional transvection. Keeping every permutation
phase inside an actual three-band residual gives sign bound
$`\max\{D_3+576,D_{\mathrm{ins}}(0,0,0,3)\}`$, independently of
active ranks and dirty width. Original clean-transition middle pairs and
higher-layer asymptotics remain open. These logical pair orbits are distinct
from the original physical middle pairs excluded by the reduction.

Every column, literal phase, dirty input, reference, actual inverse,
leakage and return requirement remains in force. The complete-frame
resource frontier is unchanged. The reviewed documentation-only
Hopf-QBP change is reconciled in [SYNC](SYNC.md); no scientific port was
needed. Final manuscript drafting is separate.

## Mandate

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge are authorized. Preserve the user's separation between scientific
work and final manuscript writing. Do not add an application-advantage
project, large simulations, QRAM, free supplied precision states, or
uncharged coherent evaluators.

The central question is whether a prescribed differential-frame completion
retains state-preparation resource tradeoffs. Exact and fault-tolerant
models price that same operator differently. Complete-input action, literal
phases, actual inverses, work return, and external references are essential.
A logical suffix is clean only on its active sector; inactive behavior
must be proved separately.

Optimal general depth remains the research destination. The user's
direction, reaffirmed on 16 October 2026 (Asia/Shanghai), is to pursue it
with confidence and curiosity, allowing pauses to consolidate what has
been learned. R90 reaches logarithmic general depth at fixed accuracy;
R92 extends a slow unbounded precision range, R94 reaches
$`O(\ell\log_2(\ell+2))`$, and R95 reaches $`O(\ell)`$ at
inverse-polynomial accuracy.
R96 extends the depth law $`O(NL/b^2+L+\ell)`$ to all precisions.
Optimal general depth and matching dimension-growing lower bounds remain
mathematical questions.
A failed route or an ended session does not retire the objective.

## Read by claim

| Need | Authoritative starting point |
|---|---|
| What we claim and under which assumptions | [Publication scope](manuscript/PUBLICATION_SCOPE.md) |
| How the argument works | [Claim-led narrative](REVIEW.md) |
| Which proof supports each claim | [Proof dependency map](docs/README.md) |
| Evidence and limits | [Verification](docs/VERIFICATION.md), [core audit](docs/CORE_CLAIM_AUDIT.md) |
| Imported premises and attribution | [Source map](docs/SOURCE_MAP.md), [related work](docs/RELATED_WORK.md) |
| One-clean count frontier and what remains unresolved | [Resource questions](docs/OPEN_PROBLEM.md) |
| Current depth costs, proof gap, and routes worth reopening | [Depth preparation map](research/depth/DEPTH_OBSTACLE_MAP.md) |
| Executed research design and outcome | [Five-hour depth prompt](research/depth/RESEARCH_SESSION_PROMPT.md), [session result](research/depth/RESEARCH_SESSION_RESULT.md) |
| Precision-uniform follow-up | [R79 proof](research/depth/UNIFORM_SUBTREE_DEPTH.md), [session record](research/depth/UNIFORM_PRECISION_SESSION_RESULT.md) |
| Joint decoder/terminal follow-up | [Common-flag depth](research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md), [recursive terminal](research/depth/RECURSIVE_TERMINAL_DEPTH.md), [bounded-density depth](research/depth/BOUNDED_DENSITY_DEPTH.md), [session outcome](research/depth/JOINT_DEPTH_SESSION_RESULT.md) |
| Branch-orientation follow-up and historical gap | [Affine orientations](research/depth/AFFINE_ORIENTATION_DEPTH.md), [exact deterministic frames](research/depth/DETERMINISTIC_FRAME_DEPTH.md), [session outcome](research/depth/BRANCH_ORIENTATION_SESSION_RESULT.md) |
| Bounded in-place conversion checkpoint | [Permutation session](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md): no qualifying general construction; scoped conversion limits and continuing direction |
| Joint orientation and recovery | [Executed design](research/depth/JOINT_ORIENTATION_RECOVERY_PROMPT.md), [session outcome](research/depth/JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md), and [geometric-band compiler](research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md) |
| General-magnitude recovery checkpoint | [Executed design](research/depth/GENERAL_MAGNITUDE_RECOVERY_PROMPT.md), [session outcome](research/depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md), [interval-dilation theorem R86](research/depth/INTERVAL_DILATION_DEPTH.md), and [selected dirty routing](research/depth/SELECTED_DIRTY_ROUTING.md) |
| Current fixed-accuracy general depth | [Dirty catalogue terminal R90](research/depth/DIRTY_CATALOGUE_TERMINAL.md) and [session outcome](research/depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md) |
| Executed precision-terminal session | [Design](research/depth/PRECISION_TERMINAL_PROMPT.md), [outcome](research/depth/PRECISION_TERMINAL_SESSION_RESULT.md): R91's parallel reader, R92's absolute precision envelope, and R93's grouped-terminal improvement |
| Executed catalogue bootstrap session | [Design](research/depth/CATALOGUE_BOOTSTRAP_PROMPT.md), [R94 proof](research/depth/PORTABLE_WORKSPACE_CATALOGUE.md), [outcome](research/depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md): carrying the full local workspace through selection permits a constant number of groups compiled directly to final accuracy |
| Executed full-flag and local-precision session | [Design](research/depth/FULL_FLAG_DEPTH_PROMPT.md), [R95 proof](research/depth/LOGARITHMIC_PRECISION_DEPTH.md), [outcome](research/depth/FULL_FLAG_DEPTH_SESSION_RESULT.md): coarse flags with sparse correction and a faster local compiler give uniform logarithmic depth through inverse-polynomial accuracy |
| Executed precision bridge and direct terminal | [Design](research/depth/BORROWED_PRECISION_BRIDGE_PROMPT.md), [R96 proof](research/depth/ALL_PRECISION_LINEAR_DEPTH.md), [outcome](research/depth/BORROWED_PRECISION_BRIDGE_SESSION_RESULT.md): a global borrowed-precision stage, fast references and directly programmed terminal give $`O(NL/b^2+L+\ell)`$ at every precision |
| Executed one-clean four-layer session | [Design](research/depth/ONE_CLEAN_FOUR_LAYER_PROMPT.md), [R97 proof](research/depth/OUTER_OVERLAP_RIGIDITY.md), [outcome](research/depth/ONE_CLEAN_FOUR_LAYER_SESSION_RESULT.md): bounded outer-axis misalignment has a positive error floor; alignment with either physical axis list gives the uniform R76 floor; the unrestricted gap remains open |
| Executed whole-terminal session | [Borrowed-work design](research/depth/TERMINAL_SOURCE_GAUGE_PROMPT.md): a fixed number of R89 reductions and dirty catalogue commutators replace the remaining recursive terminal |
| Retained causal compact-coin construction | [R89 proof](research/depth/CAUSAL_COMPACT_COIN_DEPTH.md) and [continuation record](research/depth/COIN_CONTINUATION_SESSION_RESULT.md); its coarse reduction supports R90 |
| Retained product-coin construction | [R87 proof](research/depth/COIN_DILATION_DEPTH.md), [executed factorization design](research/depth/INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md), and [session outcome](research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md) |
| Compact-coin and subset primitives | [Fixed-coin flags](research/depth/COMPACT_COIN_FLAGS.md) and [short-label subset routing](research/depth/SHORT_LABEL_SUBSET_ROUTING.md) |
| Terminal and routing refinements | [Prefix-selected terminal R88](research/depth/PREFIX_SELECTED_COIN_TERMINAL.md), [logarithmic-depth contiguous routing](research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals), and [continuation record](research/depth/COIN_CONTINUATION_SESSION_RESULT.md) |
| Endpoint compiler, reductions, and scoped limits | [Research results](research/README.md) |
| Completed changed-decoder work | [State-based QBP supplement](supplements/state_based_qbp/README.md) |

Results A–D are the principal resource theorems. Result C is the
precision-uniform one-clean count theorem and its optimal banked frontier.
The sharper two-clean complex endpoint has its own literal allocation. Retained
diagonal, multiplexor, complex-magnitude, and borrowed-workspace corollaries keep their own literal
reservations. Older proof components remain where a selected
construction actually uses them; an older headline bound is not itself a
reason to delete its lemmas.

Exploratory gap-closing notes and positive restricted cases are indexed in
`research/`. Their historical next-step passages are not current work
orders. Completed state-based work is separate from both those studies and
the prescribed-frame proof chain. Native code, tests, and exact receipts
remain in their shared directories and continue to be verified.

## Research decision and stopping rules

These rules distinguish completed claims and individual candidate stops
from the continuing general research objective. The selected count
package has reached its stopping point. The
[publication scope](manuscript/PUBLICATION_SCOPE.md) is the canonical
statement of A–D, literal workspace reservations, matching intervals,
and retained corollaries. In particular, C proves the one-clean count
frontier; D is the one-clean all-precision count/depth theorem R79 for
both frame families. The sharper two-clean complex count endpoint keeps
its distinct allocation. The older two-clean linear-precision depth
construction retains its supporting source and recovery proofs. The
[core audit](docs/CORE_CLAIM_AUDIT.md) reconstructs these proofs and their
minimal dependencies. Do not reopen resolved count questions or transfer
a supporting theorem's clean allocation to the selected construction.

### Current depth frontier

Use [the resource statement](docs/OPEN_PROBLEM.md#large-width-t-depth)
as the authoritative frontier. The selected general upper result is
[R79, uniform subtree depth](research/depth/UNIFORM_SUBTREE_DEPTH.md).
For all integer $`n\ge1,L\ge6`$, one clean qubit suffices for
both prescribed frame families, with absolute constants in

```math
D_T=O\!\left(NL/b^2+\min\{nL,(L+\log_2(n+2))\log_2(n+L+2)\}\right),
```

optimal $`T=O(\sqrt{NL}+L+NL/b)`$, and $`G=O(NL)`$ on the same
circuit. Literal widths remain $`b\ge17(L+n+7)`$ for real and
$`b\ge17(L+n+8)`$ for complex magnitude frames. The stronger
$`L\gt N/n^2`$ depth $`O(NL/b^2+L)`$ is retained.
R78 supplies the decoder geometry and full-mode error proof; R72 supplies
the one-clean terminal, fallback, and phase circuits. R79 explicitly
absorbs the exponential precision costs of the original per-call source.
At $`L=\Theta(\log n)`$ and $`b=\Theta(\sqrt{NL})`$, depth is
$`O(\log^2 n)`$ with inverse-polynomial error in n. The later fixed-accuracy
[R90 dirty catalogue theorem](research/depth/DIRTY_CATALOGUE_TERMINAL.md)
gives, for arbitrary node-dependent magnitudes and the same literal widths,

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell),\qquad \ell=\log_2(n+2).
```

It covers signed and zero coefficients and the prescribed complex leaf
phases in full initialized-isometry norm. Three R89 coarse reductions
leave a small terminal whose complete rounded program selects a finite
catalogue through arbitrary dirty registers. A commutator on two disjoint
borrowed arrays returns both arrays exactly. At fixed accuracy, optimal
square-root-width depth is now between a constant lower bound and
$`O_L(\ell)`$. R79 remains the selected manuscript theorem;
R90 does not claim constants uniform in growing L.

The later [R92 envelope](research/depth/UNIFORM_CATALOGUE_PRECISION.md)
makes the catalogue precision charge absolute, giving an additive
$`\ell+2^{C_*(L+1)^3}`$ alternative to R79 for all L.
[R93](research/depth/GROUPED_CATALOGUE_PRECISION.md) proves, uniformly on
$`6\le L\le2\ell`$,

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)
                        +L\ell/\log_2(\ell+2)\right),
```

with the same optimal T-count, Clifford budget, one clean qubit,
literal real/complex widths,
and full-isometry guarantee. Its square-root-width benchmark at
$`L=\Theta(\ell)`$ improves to $`O(\ell^2/\log\ell)`$.
The [R94 portable-workspace catalogue](research/depth/PORTABLE_WORKSPACE_CATALOGUE.md)
reaches the stronger bound on that same range:

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right).
```

It retains all those count, width and full-isometry guarantees with absolute
constants. R79 compiles each larger group directly to final accuracy;
the entire local workspace travels with the selected payload and its
borrowed slot action cancels exactly. A constant number of groups uses
the [R91 parallel word query and paired router](research/depth/PARALLEL_BILINEAR_WORD_QUERY.md).
R91 remains independent of R92 and R93, whose proofs do not use it.
The subsequent [R95 theorem](research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
improves that same range to

```math
D_T=O\!\left(NL/b^2+L+\ell\right).
```

It retains optimal T-count, $`G=O(NL)`$, one external clean bit,
the original real/complex literal widths, and full initialized-isometry
error with absolute constants. Its coarse-frame sparse correction removes
the opening's repeated spectral search. Its local compiler uses only
logarithmic initialized source scratch and much larger dirty resources;
R94's complete-payload catalogue absorbs those local costs.

The subsequent [R96 theorem](research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
proves that same bound for all integer $`n\ge1,L\ge6`$. Its
[global stage](research/depth/GLOBAL_BORROWED_PRECISION_STAGE.md) reserves
only $`O(\log L+\log^3m)`$ original logical zeros and prices every
long-address table against the original global prefix. A shorter exact
geometric source approximates the original reference with a full-space
comparison. Three reductions leave $`O(\log L)`$ layers where needed;
a [directly programmed terminal](research/depth/DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md)
uses fine masks and a paid finite native-gate schedule without enumerating
precision programs. The literal one-clean widths, optimal counts and full
error contract are unchanged. Optimal depth remains open. The selected
manuscript Result D remains R79.

The retained [R88 terminal theorem](research/depth/PREFIX_SELECTED_COIN_TERMINAL.md)
handles every real terminal height $`0\le M\le n`$, with a preserved
external prefix, one clean qubit, the literal real width, and the same
count budgets. Its depth is
$`O_L(N/b^2+\ell\log_2\log_2(M+4))`$ in full initialized-isometry norm.
For $`M=O(\ell)`$ this gives an additive
$`O_L(\ell\log_2\log_2(\ell+4))`$ contribution.
The [contiguous coin router](research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
has $`O(\ell)`$ T-depth. These component refinements initially left
R87's $`O_L(N/b^2+\ell\log_2(\ell+2))`$ headline unchanged.
R89 then batches the coarse decoder in $`O_L(\ell)`$ depth and
applies that mechanism recursively to a preserved external prefix.
Its stronger real-terminal bound is
$`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$.
Here $`\log_2^*x`$ counts repeated base-two logarithms until the value
is at most 2. R88 retains its independent terminal theorem and supplies
sharp-query and fallback interfaces to R89 and R90.

The subsequent [joint-depth session](research/depth/JOINT_DEPTH_SESSION_RESULT.md)
left the general frontier at R79 at that checkpoint. At fixed accuracy
and square-root width, R80 proves logarithmic depth for
depth-homogeneous angle magnitudes with arbitrary node signs and supplied
leaf phases; its signed balanced real case is exact. R81 gives a uniform
terminal-height bound, with a log-times-log-log contribution for a
logarithmic protected suffix. At that same accuracy and width scale, R82 gives a complete
sub-square-logarithmic circuit when
the target density has fixed positive bounds relative to a classically
specified positive product reference. Its constant-normalization
weighted correction and recursive source schedule are priced in full.
The early-decoder cost remained unresolved in that session, which records
exact failures of generic dirty-cache hoisting and row-wise fingerprint substitution, with explicit
limits on those arguments.

The later [branch-orientation session](research/depth/BRANCH_ORIENTATION_SESSION_RESULT.md)
proves R83's fixed-accuracy logarithmic depth when the unordered magnitude
pair depends only on depth and each branch choice is affine in the node
prefix. Its leaf and marker permutations have complete native circuits.
R84 gives exact $`O(\log(n+2)(1+\log^*n))`$ depth for arbitrary
node-dependent signed corner rotations. Both retain one clean qubit,
square-root T-count, and linear Clifford count at square-root dirty width.
R82 now also permits a supplied positive affine-oriented reference, with
the same fixed-density bounds and literal allocations. At that checkpoint,
arbitrary nonlinear orientation at interior bias remained unresolved:
the factorization alone did not price its in-place leaf and marker maps.
A parallel XOR lookup alone still does not supply those whole-register maps.

The later [joint orientation session](research/depth/JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
proves [R85](research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md): arbitrary
nonlinear orientations with shared unordered depth magnitudes admit
fixed-accuracy $`O_L(N/b^2+\ell\log_2(\ell+2))`$ T-depth, where
$`\ell=\log_2(n+2)`$, with the same one-clean count budgets and literal
real/complex widths. A coarse geometric band determines a subtree
unloader that exposes enough clean suffix for the local orientation
conversion. Only the coarse band then needs spectral search. The proof
includes arbitrary signs, zeros, supplied leaf phases, and full work
return. R87 subsequently attained the same depth order for arbitrary
magnitudes, and R89 improves the general bound further. R90 reaches
logarithmic depth for arbitrary magnitudes, including this promised
family. R85 retains its independent proof.

The [bounded permutation follow-up](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
did not pass its new-mechanism checkpoint and stopped construction early.
It retains a valid unary-bank relabeling identity whose recursive selected
permutation remains unpaid, and restrictions on four XOR shears, exact
encoded Clifford actions, and corner words with limited interleavers.
None is an unrestricted native lower bound; an explicit approximate
encoding also prevents extending the exact Clifford-order argument by
continuity. R85 subsequently bypasses that conversion interface through
joint recovery; it does not give a general whole-register permutation
compiler. The subsequent
[general-magnitude session](research/depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
proves R86 through a positive interval embedding, a balanced-Haar change
of basis, and paid recovery with exposed clean work. The later
[factorization session](research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
proves R87: a product of coin words and tensor Hadamards exposes the
marker depth directly, removing R86's interval-overlap amplification.
Independent technical review can strengthen the proved results alongside
that search. It does not replace the general objective. Final manuscript
drafting remains separate.

The latest obstruction result is
[R97, outer-overlap rigidity](research/depth/OUTER_OVERLAP_RIGIDITY.md).
For the explicit complete frame with $`\cos(2\theta)=1/3`$:

- R75 proves a zero-clean five-layer lower bound for integer $`L\ge3`$,
  at arbitrary dirty width. With one clean qubit, sufficiently small
  error forces exact terminal clean-sector preservation.
- R76 excludes Pauli and Hermitian-Clifford encoded midpoints under its
  stated thresholds. R77 also excludes non-Clifford midpoints when
  either two-layer encoder has bounded inner clean spread.
- R74 and R77 force both middle crossing rank and both encoder spreads
  to grow in any hypothetical arbitrarily accurate four-layer one-clean
  sequence.
- R97 gives a positive floor depending only on the smaller outer-axis
  misalignment dimension: the intersection dimension minus the number
  of individual axes from that list lying in the intersection. Bounded
  intersection is a consequence; alignment with either physical list
  gives R76's uniform floor at arbitrary intersection dimension.

The unrestricted one-clean five-layer lower bound remains open. The
remaining compressed-conjugation freezing problem must allow growing
inner spread and crossing rank, with diverging outer-axis misalignment
on both sides in every admissible native representation. Commutant lifting
and phase discreteness already hold; terminal return alone does not finish the proof. R77's exact
native counterexample rules out a universal common compatible sector flip,
not the remaining gap conjecture. Its occupied defect is a positive constant.
Its tensor extensions have growing displayed overlap and center mixing,
but an alternative Pauli-preserving three-layer realization; they do not
exclude a different native split.
The separate four-layer Hilbert–Schmidt identity-isolation question also
remains open; R75's operator-norm rigidity does not settle that weaker metric.

Optimal large-width depth is the broader unresolved question. The
[preparation map](research/depth/DEPTH_OBSTACLE_MAP.md) reconstructs the
pre-R78 recovery/query costs and the compressed freezing gap. The
[executed session](research/depth/RESEARCH_SESSION_RESULT.md) records how
native subtree reflections replaced the old overlap recovery. The
[uniform follow-up](research/depth/UNIFORM_PRECISION_SESSION_RESULT.md)
records the absolute-precision extension and its independent audits. The map also
explains why the spectator-padded witness has constant depth at fixed
accuracy: a dimension-growing lower bound needs a different witness. The
[entry-indexed bound](research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
restricts one encoding architecture, not arbitrary native circuits.
Success requires a concrete depth improvement or an invariant valid for
the unrestricted circuit model. Mechanism discovery, reformulation, and
resolving a named proof gap are legitimate steps toward that success;
more finite fixtures alone do not establish it. No practical crossover, efficient preprocessing, or
general native emitter is claimed. Internal audits do not replace external
technical review; final manuscript drafting requires the user's request.

### Continuing toward the general answer

The [whole-terminal session](research/depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
advances the complete general frame to R90's $`O_L(\ell)`$ depth
at square-root width. Three R89 coarse reductions shrink the remaining
frame enough to encode its complete rounded program in a short dirty word.
Exact XOR translations select catalogue slots, and a commutator using two
disjoint borrowed arrays cancels the unwanted array actions. The native
catalogue factors require no clean ancilla; their approximation error is
the maximum over the selected factors. The complete proof prices the
catalogue, queries, routing, actual inverses, and full work return.
R87 supplies the fine decoder and R88 the sharp queries and literal-width
fallback used through R89. The [precision-terminal follow-up](research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
proves R91–R93. Its grouped coarse operator and shared correction stream
reduce repeated long-prefix reads, with exact dirty return and an explicit
conditional-work ledger. The [catalogue bootstrap](research/depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md)
then proves R94 by carrying R79's complete local workspace through selected
transfer and compiling each group directly to final accuracy. Its constant
number of groups bypasses the old residual labels and repeated prefix reads.
The [executed full-flag session](research/depth/FULL_FLAG_DEPTH_SESSION_RESULT.md)
then improves both remaining costs and proves R95's
$`O(NL/b^2+L+\ell)`$ on $`6\le L\le2\ell`$.
The opening implements a coarse exact frame at final accuracy and adds a
signed sparse correction. The local word replaces precision-sized clean
source work by borrowed real rotations, uses exact compact reference
coins, and contracts twice before a constant number of high-precision
R79 groups. Its larger local dirty costs fit the complete-payload catalogue.
The [executed precision bridge](research/depth/BORROWED_PRECISION_BRIDGE_SESSION_RESULT.md)
then establishes [R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md),
the same depth law for every precision. Its first qualified bridge reaches
$`L\le n^{3/4}`$; direct prefix programming then removes the catalogue
restriction. Fine precision masks share a fixed local grammar, and only
the shorter coarse native word is read instruction by instruction.
All controls, original conditional zeros, dirty pools and actual inverses
are paid. Optimal general depth and matching dimension-growing lower
bounds remain open. The precision restriction in this upper law is now
resolved; do not repeat it or the earlier R94 opening/local questions.
The [executed four-layer design](research/depth/ONE_CLEAN_FOUR_LAYER_PROMPT.md)
produces [R97](research/depth/OUTER_OVERLAP_RIGIDITY.md): a gap at bounded
outer-axis misalignment and a uniform gap when the intersection is aligned
with either physical axis list. The
[session outcome](research/depth/ONE_CLEAN_FOUR_LAYER_SESSION_RESULT.md)
records the remaining unbounded misalignment on both sides and the exact
limits of proposed descent and center reductions. The unrestricted
five-layer obstruction remains open; it is a constant lower-bound target,
separate from dimension-growing depth.
No efficient preprocessing, practical crossover, or scalable native emitter is claimed.
Preserve one clean qubit,
square-root T-count, linear Clifford count, and the complete operator
contract throughout. The [whole-register reduction](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md#1-the-sufficient-target-and-its-scope)
remains sufficient, but is not required by R85's complete recovery word.

Separate discovery from certification. Begin a bounded exploration with
a question, a changed assumption or representation, and the uncertainty
it aims to resolve. An unfinished identity can be a useful research
lead when its missing step is explicit. Before claiming an improvement,
complete the native identity and resource, error, and return proof.
The joint session's successful change was to use workspace exposed on
the reached spectral band, and to prove return for the complete word.
R86 uses that same freedom for a complete local frame compiler, then
erases the fine label through its actual inverse and forward word.
Its normalization bound and native recovery were proved together.
The subsequent [factorization design](research/depth/INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md)
has been executed. Its weighted-chain route is retained separately in
the [session record](research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md);
R87 instead uses spectral readout on the embedded logical image.
The interval-overlap normalization is no longer a cost of the current
general circuit.

Use the [full-circuit contract](docs/OPEN_PROBLEM.md#full-circuit-requirement):
a native identity or encoding rule and symbolic resource, error, and cleanup
recurrences must precede new fixtures. At large width, price the complete
recovery, query, arithmetic, and terminal ledger in the preparation map.
The older grouped route must reduce both its logical-transport and query
allowances or bypass that interface; the current column construction
already bypasses it.
The one-clean count frontier retains literal dirty reservations and charges
every precision source and query. It supplies no new depth, efficient
classical preprocessing, or general emitter guarantee.

An unrestricted lower bound is also a legitimate resolution, but needs an
invariant valid for the actual circuit model. Existing interface failures
are not unrestricted lower bounds. Stop a candidate that retains its old
asymptotic cost, needs unpriced clean history, or cannot return its work.

At a checkpoint, decide whether to continue, reformulate, or pause that
candidate. Record what was learned, the exact unresolved step, and a
question from which work can resume. Exhausting the current list of
ideas is a session outcome, not evidence that the target is unattainable.
Use bounded sessions to learn and avoid repeating unchanged failures;
discovering the next idea remains part of the authorized research.
Final writing begins only when requested.

## Refurnishing policy

Each selected claim has one authoritative proof home. Front pages summarize
and link; they do not reproduce full proofs. Remove duplicated exposition
when its complete argument remains in that home. Preserve unique attempts,
negative results, positive special cases, and their evidence in the research
archive, with an outcome and reopening condition. Preserve stable source
identifiers and all required dependency links.

Repository paths and section links must survive moves. The full mathematics
and presentation checks include archived notes and the completed supplement.
Changes to documentation tests should reflect the new reading structure;
scientific tests and receipt hashes must not be weakened or regenerated.

## Restore and verify

Start with `git status --short --branch` and the relevant proof chapter.
For substantive theorem/circuit changes, use:

```bash
python -m compileall -q compiler_robust_hopf scripts tests validate.py
python scripts/reviewer_walkthrough.py
python scripts/coarse_frame_native_example.py
python scripts/complex_coarse_native_example.py
python scripts/residual_qbp_native_example.py --q 16
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
python scripts/unified_resource_ledger.py --n 12 --format json
python scripts/strict_zero_echo_ledger.py --n 12 --format json
```

The [presentation guide](assets/README.md#rendering-checks) gives the complete
rendered-document check. Record checks on the actual reviewed tree, then
merge only after required gates pass. Do not confuse finite validation
with an asymptotic proof or a scalable native emitter.
