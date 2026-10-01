# Continuing research workspace

This is the entry point for resuming work when a previous conversation or
execution workspace is unavailable. The proofs and research decisions live
in the repository; a conversation summary is only a retrieval aid.

The 2026-10-01 source-carry pass continues from main commit
`e6e9cd5d993d78e7c78198169e82b257500d4968`. Check the current branch and
later commits before continuing. New source-width constructions and scoped
cost obstructions refine the research decision; the compiler frontier is
unchanged.

## Mandate and reading order

Continue the theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository write permission is limited to this repository. Keep manuscript
writing and release work outside the current research pass.

1. Read [research status and the next bounded pass](docs/OPEN_PROBLEM.md).
   This is the primary home for the current frontier and research decision.
2. Read the [overview](README.md), [technical narrative](REVIEW.md), and
   [publication scope](manuscript/PUBLICATION_SCOPE.md) for established claims.
3. For the active pass, read [conditional-suffix grouping](docs/CONDITIONAL_SUFFIX_COMPILER.md),
   especially its complete group contract and resource sums, followed by
   [source-reuse limits](docs/SOURCE_REUSE_LIMITS.md) and
   [the completed repair audit](docs/RESIDUAL_ASSEMBLY.md#9-a-repair-word-without-an-ill-conditioned-transported-basis).
4. Consult the [verification map](docs/VERIFICATION.md),
   [source map](docs/SOURCE_MAP.md), and [related work](docs/RELATED_WORK.md)
   before extending a claim or changing its attribution.

## Exact point of continuation

The exact all-workspace size/depth theorem and the sufficient-clean
fault-tolerant T-count theorem are established. The latest antichain and
sparse nested-update results give linear precision cost for their stated
native-baseline promises. Generic rounding need not satisfy either promise.
The complete-frame and fixed-parameter QBP guarantees already stand.

For the prescribed complete real frame, the selected unresolved endpoint is

```math
N=2^n,\qquad a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

```math
\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n)),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

The upper bound also works with one clean qubit. The two promised update
classes have not narrowed this generic gap. Optimal T-depth remains a
separate open question.

The branching passes have established a complete one-signal coupled boundary
and an explicit commutator repair supported on at most four transported
modes. It needs no division by a vanishing boundary norm. Its two child
calls cancel with the anchored merge, leaving one child call and the
original fine-precision local target wrappers. This cancellation survives
matched actual inverse substitutions on the whole dirty space. The small
repair norm does not discount the complete merge's local precision error.
A separate native common-conjugator word fails even after retuning its
masks. The generic resource gap is unchanged.

The revision replaces the previous two-depth task: a fixed fork already
costs $`O(L)`$, and the best compiler already shares precision across
much larger groups. Rechecking low rank, basis access, or child cancellation
does not target its remaining $`\ell_*(n)`$ factor. Keep those results
as completed infrastructure. The recent passes have not narrowed the
generic gap or established that its resolution is close.

## Latest result and next bounded pass

The [source-carry audit](docs/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation)
solves the width-only part of the previous task. Reordering the source
loader preserves its certified coefficient grid and gives literal native
bridges whose T costs telescope. Initial/final loaders and all monotone
width changes cost exactly $`2(m_{\max}-1)`$ T gates, with Clifford
cost within $`O(NL)`$. This is a complete-space identity, including
occupied flags and arbitrary dirty correlations.

The interior programs remain charged. A coefficient arising from a legal
group gives a transformed mask requiring
$`T\ge m-2\log_2s-9`$ at its stated fine accuracy when synthesized
separately. A second, flag-correlated source code also has cheap width
changes, but a legal scalar query leaves its code by constant norm.
Renewing the complete source syndrome costs $`\Omega(L)`$ at compilation
accuracy. These are scoped interface bounds, not additive frame lower
bounds. Neither candidate yet provides a cheap complete group program.

Keep the existing groups: their total table/coarse-program work is already
$`O(N)`$. The next bounded task is a **changed source-and-program word
or a logical-dependent encoded boundary** for two unequal groups, with a
recurrence closed through a third and arbitrarily many more. A successful
joint interior may spend $`O(L)`$ once; its total target is
$`O(\sum_gQ_g+L+\mathrm{poly}(n))`$. Dirty-only outer loaders do
not supply precision to otherwise precision-free independent group bodies.
Individual group action may be deferred instead of renewing each source
or exposing each transformed mask separately.

1. Derive the complete boundary action and native cost recurrence first.
   Explain where the internal precision charge occurs and why it is paid
   only once. No successful joint interior has yet been constructed.
2. Retain at most two clean flags and $`b=L+n+7`$ arbitrary dirty wires.
   Source tails can become selectors, never initialized work. Include
   changing suffix predicates, actual query unloading, literal inverses,
   and complete dirty/reference return. Keep the fixed deepest-layer tail
   with its existing $`O(N+L)`$ compiler.
3. Only after deriving a changed parametric word, test unequal widths,
   changed logical addresses, a third group, occupied ports, and a real
   target over a complex native baseline. These finite checks falsify a
   recurrence; its asymptotic bound needs a proof.

The [research decision](docs/OPEN_PROBLEM.md#revision-decision-and-next-bounded-pass)
gives the full ledger and stop conditions. Stop a candidate that recharges
precision per group, requires full syndrome renewal per query, expands the linear
table work, or assumes unproved clean work or return. Change representation
instead of repeating the now-completed width, mode-closure, or commutator
repair audits. No reset, free catalyst, or target oracle is supplied.
Other fully charged constructions remain allowed, including a linear-T
endpoint circuit with larger Clifford cost. The generic gap has not
narrowed, and these results do not establish that its resolution is close.

## Restore and verify

Run from the repository root with the dependencies in `requirements.txt`:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

The baseline before this pass had 251 tests and passed Python 3.11/3.13
CI, all four fault-tolerant receipt suites, and rendered presentation.
The source-carry pass adds 12 focused checks for native phases, width
bridges, complete scalar words, legal grouped coefficients, correlated
codes, and syndrome extraction. All 263 tests passed locally on Python
3.12, as did the reviewer walkthrough, offline synchronization check, and
protected-math rendering handoff. CI on the merged change records the
Python 3.11/3.13, exact-receipt, and full rendered-presentation results.
Its proof sections received separate internal mathematical checks; this
is not external peer review. Finite fixtures check identities and failure
controls; native resource bounds remain analytic arguments.

Update this entry point when the task changes, while keeping detailed
research conclusions in [the existing checkpoint](docs/OPEN_PROBLEM.md)
and proofs in their existing chapters. Keep the frontier and the unresolved
step explicit so continuation does not depend on access to an old chat.
