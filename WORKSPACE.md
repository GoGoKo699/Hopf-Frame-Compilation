# Continuing research workspace

This is the entry point for resuming work when a previous conversation or
execution workspace is unavailable. The proofs and research decisions live
in the repository; a conversation summary is only a retrieval aid.

The 2026-10-01 small-product pass continues from main commit
`efb56deda5660f6726f8d93c8f551c272c1b1567`. Check the current branch and
later commits before continuing. Small complete examples now guide the
next construction; the compiler frontier is unchanged.

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
   [small products](docs/SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis),
   [the Cayley residual](docs/ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation),
   and [the completed repair audit](docs/RESIDUAL_ASSEMBLY.md#9-a-repair-word-without-an-ill-conditioned-transported-basis).
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

Start from small complete examples and let them suggest a construction.
The positive [fixed-address example](docs/SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis)
shows that arbitrarily many noncommuting rotations on one target can be
multiplied into four quaternion coordinates and compiled once. This uses
the existing multiplexor theorem and its stated workspace allocation.
Changing the address or adding overlapping target pairs breaks that
fixed two-mode description; noncommutation alone is not the obstruction.

Four- and eight-mode branching examples led to a
[Cayley representation](docs/ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation)
of the complete residual $`\mathcal R=C^\dagger W`$. Its
skew-Hermitian generator $`K=(\mathcal R-I)(\mathcal R+I)^{-1}`$
has an exact tree recursion with two-dimensional local corrections. It
handles actual complex native coarse words, has $`O(N)`$ classical
recursive data, and avoids division by vanishing defects. The retained
uniform coarse approximation bounds its small inverses. Inverse Cayley
is stable in operator norm, but has no established cheap native circuit.

The earlier source-width result remains useful: outer loaders and monotone
bridges cost $`2(m_{\max}-1)`$ T gates. Separately exposed transformed
masks and repeated complete syndrome renewal still cost precision. That
pass solved transport, not the joint program. The new Cayley coordinates
also do not yet improve the generic T-count.

The next experiment is a **joint native program for the small Cayley
residual**, first on four logical modes and then eight. Include a complex
native baseline, zero defects, and nonzero changes in both branches. Try
to use the combined recursive data without synthesizing every local
update separately. A generic resolvent or polynomial conversion must count
every generator call; its bounded condition number is not a free circuit.

1. Small examples may come before a general recurrence. Compare complete
   matrices and actual source counts, including all rejected flag ports,
   dirty inputs, and literal inverses. A fixed-size $`O(L)`$ synthesis is
   already available and alone does not establish a scalable gain.
2. If a pattern survives, derive its recurrence as the tree or number of
   existing groups grows. The target remains $`O(N+L)`$ T gates with
   the current linear table/coarse-program work; allow one internal
   $`O(L)`$ precision charge and deferred individual group action.
3. Keep at most two clean flags and $`b=L+n+7`$ arbitrary dirty wires.
   Track actual unloading, changing suffix predicates, released source
   tails, and full dirty/reference return. Keep the fixed deepest-layer
   tail with its existing compiler. No reset, free catalyst, target-frame
   oracle, or newly initialized logical sector is supplied.

The [research decision](docs/OPEN_PROBLEM.md#revision-decision-and-next-bounded-pass)
gives the ledger and stop conditions. Stop a native candidate that merely
returns to fine local target wrappers, renews precision per group, expands
the linear table work, or assumes free transport/inversion. Retain useful
operator identities, then change the native candidate. Other fully charged
constructions remain allowed, including a linear-T endpoint circuit with
larger Clifford cost. These examples give a specific next experiment;
they do not establish that the generic gap is close to resolution.

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

The source-carry baseline has 263 tests and passed Python 3.11/3.13 CI,
all four exact-receipt suites, and rendered presentation. The small-product
pass adds focused checks for fixed-address product compression and the
complete Cayley residual; matrices have at most sixteen logical modes.
All 273 tests passed locally, including ten new small-example checks.
The reviewer walkthrough, offline synchronization, and documentation
checks also passed. CI on the merged change records the Python 3.11/3.13,
exact-receipt, and rendered-presentation results. Separate internal reviews
checked the algebra and conditioning; this is not external peer review.
Finite fixtures do not establish native asymptotic resource bounds.

Update this entry point when the task changes, while keeping detailed
research conclusions in [the existing checkpoint](docs/OPEN_PROBLEM.md)
and proofs in their existing chapters. Keep the frontier and the unresolved
step explicit so continuation does not depend on access to an old chat.
