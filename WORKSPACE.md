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
| Branch-orientation follow-up and current constructive gap | [Affine orientations](research/depth/AFFINE_ORIENTATION_DEPTH.md), [exact deterministic frames](research/depth/DETERMINISTIC_FRAME_DEPTH.md), [session outcome](research/depth/BRANCH_ORIENTATION_SESSION_RESULT.md) |
| Bounded in-place conversion checkpoint | [Permutation session](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md): stopped without a qualifying general construction; scoped conversion limits and reopening condition |
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

The selected count package has reached its stopping point. The
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
leaf phases; its signed balanced real case is exact. R81 improves only
the deepest logarithmically many general layers to a log-times-log-log
contribution. At that same accuracy and width scale, R82 gives a complete
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
the same fixed-density bounds and literal allocations. The current
constructive gap is arbitrary nonlinear orientation at interior bias:
the factorization $`W=P W_0 Q`$ still needs priced in-place leaf and marker maps, or a
direct complete-word construction. A parallel XOR lookup alone does not
supply them. The branch session records the precise reopening condition.

The [bounded permutation follow-up](research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
did not pass its new-mechanism checkpoint and stopped construction early.
It retains a valid unary-bank relabeling identity whose recursive selected
permutation remains unpaid, and restrictions on four XOR shears, exact
encoded Clifford actions, and corner words with limited interleavers.
None is an unrestricted native lower bound; an explicit approximate
encoding also prevents extending the exact Clifford-order argument by
continuity. R79 and the arbitrary-orientation frontier are unchanged.
Without a new qualifying mechanism, the next phase is independent
technical review of the proved results, rather than another automatic
construction session. Final manuscript drafting remains separate.

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
Further research must supply a concrete depth improvement or an invariant
valid for the unrestricted circuit model. More finite fixtures alone are
not a next step. No practical crossover, efficient preprocessing, or
general native emitter is claimed. Internal audits do not replace external
technical review; final manuscript drafting requires the user's request.

### Reopening research requires a qualifying mechanism

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

Reopen for a concrete qualifying idea, a discovered defect, or an explicit
change of scope. Do not restart exploration merely to keep activity going.
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
