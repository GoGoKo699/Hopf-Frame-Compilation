# Continuing research workspace

This is the entry point for resuming work when a previous conversation or
execution workspace is unavailable. The proofs and research decisions live
in the repository; a conversation summary is only a retrieval aid.

The 2026-10-01 revision reviewed main commit
`b0eab2af22287c9fa8384df3d811afaed21cda3b`, including the complete
transported-repair audit. Check the current branch and later commits before
continuing. This revision changes the research decision, not the theorems.

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

## Next bounded pass: reuse precision across existing groups

Keep the current groups and their linear total table/coarse-program work.
Seek a complete circuit interface that composes across arbitrarily many
adjacent groups, rather than another fixed-height improvement. The source
widths are $`m_g=L+\lfloor r_g/4\rfloor+8`$ and vary monotonically
in physical execution order. One sufficient, **unproved** mechanism would
pay $`O(\max_g m_g)`$ initially and finally, and only
$`O(|m_{g+1}-m_g|+\mathrm{poly}(n))`$ between groups. The widths'
total variation is $`O(n)`$. Together with the retained linear table
work this would give $`T=O(N+L)`$. Leave the fixed deepest-layer tail
with its existing $`O(N+L)`$ compiler.

1. First specify the shared boundary and a changed native transition word.
   Show its complete accepted/rejected action and why its output satisfies
   the same interface for the next, unequal-size group. State the actual
   precision-charge recurrence before adding more numerical fixtures.
2. Give the live-wire schedule at every transition. The source tail may need
   to become selector workspace; holding the largest source throughout is
   not automatically compatible with $`b=L+n+7`$. Include actual unloading,
   literal inverses, changing suffix predicates, arbitrary dirty/reference
   inputs, and at most two external clean flags. No reset, free catalyst,
   target oracle, or newly initialized logical sector is available.
3. Test the proposed native word on two unequal groups and then a third,
   with nonzero changes and a real target over a complex native baseline.
   These are falsification checks for the stated unbounded composition
   rule, not evidence of an asymptotic bound by themselves.

The detailed sufficient ledger and stop conditions are in
[the research decision](docs/OPEN_PROBLEM.md#revision-decision-and-next-bounded-pass).
Stop this mechanism if it charges another length-L preparation at each
boundary, shifts that cost into separately synthesized masks, expands the
linear tables, or assumes unavailable clean work or intermediate return.
If no closed interface can be specified, record that absence and change
the representation; do not repeat the completed canonical repair analysis.
This is a chosen sufficient route, not a lower bound or a necessity theorem.
A linear-T endpoint construction with larger fully charged Clifford cost
would also settle the T-count question. Other promised-family results are
not the priority without a reduction covering arbitrary updates.

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

The reviewed baseline has 251 tests. Both Python 3.11 and 3.13 CI jobs,
all four fault-tolerant receipt suites, and rendered presentation passed.
Local validation, the reviewer walkthrough, final repair tests, and the
offline synchronization check also passed during the repair pass. The
proof sections received separate internal mathematical checks; this is
not external peer review. Record fresh results after changes. Finite
fixtures check identities and failure controls; native resource bounds
remain analytic arguments.

Update this entry point when the task changes, while keeping detailed
research conclusions in [the existing checkpoint](docs/OPEN_PROBLEM.md)
and proofs in their existing chapters. Keep the frontier and the unresolved
step explicit so continuation does not depend on access to an old chat.
