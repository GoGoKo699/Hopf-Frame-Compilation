# Continuing workspace

The repository is organized around the
[claims selected for the manuscript](manuscript/PUBLICATION_SCOPE.md).
The [resource statement](docs/OPEN_PROBLEM.md) records the proved
one-clean count theorem in Result C, its optimal banked count, the
sharper two-clean complex endpoint, all-precision count/depth bounds,
and precision-uniform subtree depth. The
[research index](research/README.md) groups native primitives and promised families, full-frame reductions, source programs,
and representation limits. The [evidence catalogue](docs/reference/VERIFICATION_CATALOGUE.md)
links their implementations and checks. Each result has one proof home;
a completed source isometry or compact representation does not by itself
supply a linear-count native recovery circuit.
Final manuscript drafting is a separate task.

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

The general construction remains the research destination. The user's
direction, reaffirmed on 16 October 2026 (Asia/Shanghai), is to pursue it
with confidence and curiosity, allowing pauses to consolidate what has
been learned. The beauty of the question motivates the search; whether
the conjectured improvement is true remains a mathematical question.
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

Use [the resource statement](docs/OPEN_PROBLEM.md#constant-depth-frontier)
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
$`O(\log^2 n)`$ with inverse-polynomial error in n. At fixed accuracy,
optimal square-root-width depth remains between a constant lower bound
and $`O_\eta(\log^2 n)`$.

The subsequent [joint-depth session](research/depth/JOINT_DEPTH_SESSION_RESULT.md)
leaves that general frontier unchanged. At fixed accuracy and square-root
width, R80 proves logarithmic depth for
depth-homogeneous angle magnitudes with arbitrary node signs and supplied
leaf phases; its signed balanced real case is exact. R81 gives a uniform
terminal-height bound, with a log-times-log-log contribution for a
logarithmic protected suffix. At that same accuracy and width scale, R82 gives a complete
sub-square-logarithmic circuit when
the target density has fixed positive bounds relative to a classically
specified positive product reference. Its constant-normalization
weighted correction and recursive source schedule are priced in full.
A general improvement still needs a cheaper complete early
decoder or native query. The session records exact failures of generic
dirty-cache hoisting and row-wise fingerprint substitution, with explicit
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
return. Logarithmic depth for this family and a corresponding improvement
for arbitrary node-dependent magnitudes remain open.

The [bounded permutation follow-up](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
did not pass its new-mechanism checkpoint and stopped construction early.
It retains a valid unary-bank relabeling identity whose recursive selected
permutation remains unpaid, and restrictions on four XOR shears, exact
encoded Clifford actions, and corner words with limited interleavers.
None is an unrestricted native lower bound; an explicit approximate
encoding also prevents extending the exact Clifford-order argument by
continuity. R85 subsequently bypasses that conversion interface through
joint recovery; it does not give a general whole-register permutation
compiler. The unrestricted R79 frontier is unchanged.
Independent technical review can strengthen the proved results alongside
that search. It does not replace the general objective. Final manuscript
drafting remains separate.

The latest obstruction result is
[R77, inner clean-spread rigidity](research/depth/ONE_CLEAN_INNER_SPREAD.md).
For the explicit complete frame with $`\cos(2\theta)=1/3`$:

- R75 proves a zero-clean five-layer lower bound for integer $`L\ge3`$,
  at arbitrary dirty width. With one clean qubit, sufficiently small
  error forces exact terminal clean-sector preservation.
- R76 excludes Pauli and Hermitian-Clifford encoded midpoints under its
  stated thresholds. R77 also excludes non-Clifford midpoints when
  either two-layer encoder has bounded inner clean spread.
- R74 and R77 force both middle crossing rank and both encoder spreads
  to grow in any hypothetical arbitrarily accurate four-layer one-clean
  sequence. These necessary size conditions do not exclude such a sequence.

The unrestricted one-clean five-layer lower bound remains open. The
remaining compressed-conjugation freezing problem must allow growing
inner spread and crossing rank. Commutant lifting and phase discreteness
already hold; terminal return alone does not finish the proof. R77's exact
native counterexample rules out a universal common compatible sector flip,
not the remaining gap conjecture. Its occupied defect is a positive constant.
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

R85 reaches the sub-square-logarithmic target for arbitrary nonlinear
branch orientation with a shared unordered magnitude pair at each depth.
The next local question is whether its remaining coarse-band label can
be erased with a constant number of native source layers. This would
remove the repeated-comparison cost in early recovery; logarithmic total
depth would also require compatible logarithmic terminal and routing
schedules. Arbitrary node-dependent magnitudes and optimal
unrestricted depth remain broader questions. Preserve one clean qubit,
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
The next design should ask whether that same freedom can remove the
remaining coarse search or replace the local product reference for
unequal node-dependent magnitudes. Neither further saving is established.

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
