# Precision-uniform follow-up: execution record

[R79 proof](UNIFORM_SUBTREE_DEPTH.md) · [Original session](RESEARCH_SESSION_RESULT.md) · [Current frontier](../../docs/OPEN_PROBLEM.md)

The user authorized this continuous research session on 2026-10-15,
with a five-hour ceiling. It started at 12:29:36 UTC from main
`88008c6bd95c250211212678550b33d98441065b`, the merged R78 result.
The task was to replace accuracy-dependent constants with an absolute
all-precision bound while retaining one clean qubit, literal widths,
optimal T-count, the Clifford budget, and full operator/return error.

## Outcome and decisive steps

**Local result R79** establishes the all-precision theorem in its
[single proof home](UNIFORM_SUBTREE_DEPTH.md#1-the-all-precision-theorem).
In particular, the $`O(\log^2 n)`$ large-width depth now covers
$`L=\Theta(\log n)`$, hence inverse-polynomial error in n. The
complete real and prescribed complex magnitude families keep their
original one-clean widths. The stronger high-precision clause and
worst-case count optimality are retained.

Three details make the precision extension work:

1. The absolute source ledger exposes $`2^L`$ and $`2^{2L}`$
   costs, including standalone overhead. A protected suffix of length
   $`M=\lceil A[L+\log_2(n+L+2)]\rceil`$ suppresses the large
   tables; the standalone terms fit for a fixed sufficiently small
   $`L/n`$. One absolute dimension cutoff handles the whole range.
2. The one-clean shared binary terminal stream handles most of the
   remaining layers. Only $`O(\log n)`$ deepest layers use direct
   precision sources, avoiding a quadratic cost in L.
3. R72 already supplies the desired order outside this low-precision,
   large-dimension range. Its complete circuits cover those cases,
   including unbounded precision at bounded dimension. A separate width
   splice preserves every literal eligible allocation.

The original per-call gradient is retained. A new shared-gradient
invariant is unnecessary. No new external premise, circuit identity,
finite fixture, or scalable native emitter is introduced.

## Independent review

Separate read-only audits checked the absolute source/anchor costs,
terminal/all-width composition, and adversarial precision/leakage issues.
All three then reviewed the written proof and found no blocking defect.
The review explicitly checked the standalone exponential terms, occupied
inactive work, full-mode reflections, actual inverses, one-clean reuse,
finite-dimension/unbounded-precision fallback, half-error reselection,
literal complex phases, and the matching interval. The terminal workspace
constant is chosen after the enlarged early-bank constant.

These are internal audits, not external technical review. The unchanged
projector and source/terminal fixtures test their finite algebra; the
uniform cutoff and resource estimates are proved analytically.

## Validation and remaining limits

Local validation passed **922 tests**, including the unchanged decoder
checks. The focused documentation/decoder run passed 45 checks. All four
fault-tolerant receipt suites reproduced their immutable fixtures.
Python compilation, the reviewer walkthrough, coarse real and complex
native examples, the residual-QBP example, offline provenance, and both
required resource ledgers passed. The integrating pull request records
the reviewed tree, rendered-document audit, and remote CI results.

The source catalogues and claim maps select R79 while retaining R78's
decoder proof and R72/R66's supporting constructions. Historical session
records remain separate from the current frontier.

The unrestricted optimal depth at large width remains open. Improving
fixed-accuracy $`O(\log^2 n)`$ to $`O(\log n)`$ would require
reducing both the sequential spectral-query cost and the terminal-layer
cost, or a construction that bypasses both. The current ledger does not
prove those costs necessary. No practical crossover, efficient classical
preprocessing, general emitter, or external priority claim is added.
