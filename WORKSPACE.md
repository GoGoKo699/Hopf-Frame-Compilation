# Continuing workspace

The repository is organized around the
[claims selected for the manuscript](manuscript/PUBLICATION_SCOPE.md).
The [resource statement](docs/OPEN_PROBLEM.md) records the proved
one-clean count theorem in Result C, its optimal banked count, the
sharper two-clean complex endpoint, and all-precision count/depth bounds
with sublinear large-width special cases. The
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

The selected count package has reached its stopping point.
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

The selected depth result is the
[one-clean all-precision theorem](research/depth/ONE_CLEAN_DEPTH.md), R72.
For every $`n\ge1`$, integer $`L\ge6`$, one clean flag, and
$`b\ge17(L+n+7)`$, one complete real-frame circuit has full
initialized-isometry error at most $`2^{-L}`$ and absolute-constant bounds

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)),
```

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
```

Count is optimal at every eligible width; depth also matches through
$`b\le\sqrt{NL/F(n,L)}`$ whenever that interval is nonempty.
For $`L>N/n^2`$, depth strengthens to $`O(NL/b^2+L)`$ with
the same count and Clifford budgets. R72 prices native one-clean words,
relocates the terminal rejection flag into the protected logical bank,
and supplies its own width and precision splices, including bounded
dimensions. Source preparation, coherent queries, actual inverses,
conditional work, and external-reference error are all charged. The earlier
[real theorem R70](research/depth/ALL_PRECISION_DEPTH.md) and
[complex theorem R71](research/depth/COMPLEX_FRAME_DEPTH.md) retain their
two-clean statements and remain supporting proofs.

R72 gives the same simultaneous bounds for $`D_\phi W_{\mathbb R}`$
at $`b\ge17(L+n+8)`$ with one clean flag. Its literal-diagonal
schedule has depth $`O(NL/b^2+L+\log_2(n+2))`$ and preserves
optimal count. Both half-error factors reuse the sole flag and arbitrary
dirty pool; the full-isometry hybrid includes prior leakage. The stronger
high-precision bound also holds for this complex magnitude family.
Independently prescribed phases include zero-amplitude leaves. This is
not arbitrary complex-unitary synthesis.

For some absolute $`c>0`$, the
[compressed-column theorem](research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains its two-clean assumption and stronger logarithm $`\log_2(n+2)`$ on
$`6\le L\le cn`$. At sufficient $`b=\Theta(\sqrt{NL})`$,
shrinking error retains sublinear depth for $`L=o(n/\log(n+2))`$;
the fixed-accuracy all-n statement remains valid. Its early gradient and
terminal binary source are separate resources. The grouped terminal
proof and its exact coarse-program interface remain at their original paths.

The one-clean real and phase-dressed all-precision extensions are complete at the
stated scope. Further
research must address a concrete depth improvement or a valid unrestricted
lower bound; producing more fixtures is not a default next step. The
[entry-indexed bound](research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
remains an architecture-specific normalization bound, not an unrestricted
native-depth lower bound. Optimal depth remains open. No practical
crossover, efficient preprocessing, or general native emitter is claimed.
Internal audits do not replace external technical review. Final manuscript
drafting remains a separate task requiring the user's request.

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
