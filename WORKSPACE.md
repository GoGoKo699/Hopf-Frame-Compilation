# Continuing workspace

The repository is organized around the
[claims selected for the manuscript](manuscript/PUBLICATION_SCOPE.md).
The [resource statement](docs/OPEN_PROBLEM.md) records the proved
one-clean count theorem in Result C, its optimal banked count, the
sharper two-clean complex endpoint, and sublinear large-width T-depth. The
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

The selected scientific package has reached its stopping point.
[Result C](research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives $`T=O(N+L)`$ and $`G=O(NL)`$ for every integer $`L\ge6`$
and $`n\ge1`$, at base dirty reservation $`B=L+n+7`$ for real
frames or $`B=L+n+8`$ for complex magnitude frames. At $`b\ge2B`$,
the [banked theorem](research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
matches the full worst-case count frontier
$`\Theta(\sqrt{NL}+L+NL/b)`$ with one clean qubit.
The [sharp two-clean endpoint](research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
keeps $`b=N+n+7`$ for both families at $`L=N,n\ge3`$.
A, B, and D retain their separate resource statements.

The [current proof reconstruction](docs/CORE_CLAIM_AUDIT.md) records
the native contracts and minimal dependencies. The split at $`L=N/n^2`$
uses the direct layerwise baseline below the cutoff and growing sources
above it. Iterated grouping is a retained standalone result, not a
dependency of C's current count proof. The conditioning lemma in its
chapter remains necessary. Historical audits stay in the reference
catalogue; current claims are stated in the publication scope.

The user explicitly reopened large-width depth research. The resulting
[compressed-column theorem](research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
now has an absolute-constant precision-uniform refinement. For some absolute
$`c>0`$, integers $`6\le L\le cn`$, two clean flags, and every
$`b\ge17(L+n+7)`$, it gives full initialized-isometry error at most
$`2^{-L}`$ and simultaneous bounds

```math
T=O(\sqrt{NL}+NL/b),\qquad G=O(NL),\qquad
D_T=O\!\left(NL/b^2+\min\{nL,\sqrt n\log_2(n+2)+L^{3/2}\}\right).
```

The [shared-gradient proof](research/depth/SHARED_GRADIENT_RECOVERY.md)
pays native preparation and actual return once across early recovery.
All arithmetic/coefficient tables, complete-column recovery, exact inactive
guards, simultaneous conditional scratch, and the final tail are charged.
The [grouped terminal proof](research/depth/GROUPED_TERMINAL_DEPTH.md)
prefetches complete coarse programs and pays one fine residual per group;
its final capped layers preserve the stated count and Clifford budgets.
At sufficient
$`b=\Theta(\sqrt{NL})`$, shrinking error retains sublinear depth
when $`L=o(n^{2/3})`$. Count is optimal throughout the stated range;
depth also matches through $`b\le\sqrt{NL/F(n,L)}`$ whenever
that interval is nonempty, where

```math
F(n,L)=\min\{nL,\sqrt n\log_2(n+2)+L^{3/2}\}.
```

The retained all-precision D theorem covers
other precisions; the fixed-accuracy statement covers every n.

The [entry-indexed encoding bound](research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
fixes the scope of the remaining square-root normalization cost; it is
not an unrestricted native-depth lower bound. Optimal depth remains open.
The small absolute precision constant and dimension cutoff carry no
practical crossover claim. A general native emitter, practical crossovers
and further fixed fixtures are not default requirements. The internal
audit does not replace external technical review.

### Reopening research requires a qualifying mechanism

Use the [full-circuit contract](docs/OPEN_PROBLEM.md#full-circuit-requirement):
a native identity or encoding rule and symbolic resource, error, and cleanup
recurrences must precede new fixtures. At large width, a complete improvement
must reduce both the logical-transport and query allowances or bypass them.
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
