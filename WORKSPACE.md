# Continuing research workspace

This is the entry point for resuming work when a previous conversation or
execution workspace is unavailable. The proofs and research decisions live
in the repository; a conversation summary is only a retrieval aid.

The recovery on 2026-10-01 started from main commit
`9baa1e1a1b3cf8c650e4b9febf90ec6b0b51dda8`, “Consolidate the endpoint research
checkpoint and next test.” Check the current branch and later commits before
continuing. This handoff adds no compiler theorem.

## Mandate and reading order

Continue the theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository write permission is limited to this repository. Keep manuscript
writing and release work outside the current research pass.

1. Read [research status and the next bounded pass](docs/OPEN_PROBLEM.md).
   This is the primary home for the current frontier and research decision.
2. Read the [overview](README.md), [technical narrative](REVIEW.md), and
   [publication scope](manuscript/PUBLICATION_SCOPE.md) for established claims.
3. For the active pass, read [residual assembly](docs/RESIDUAL_ASSEMBLY.md),
   [weighted transport](docs/WEIGHTED_TRANSPORT_BLOCK.md),
   [tree transport](docs/ENDPOINT_TREE_TRANSPORT.md), and
   [source-reuse limits](docs/SOURCE_REUSE_LIMITS.md).
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

The next deliverable is **joint precision synthesis across two adjacent
depths, with an explicit native source-sharing word and recursive ledger**. Read
[the coupled completion and repair](docs/RESIDUAL_ASSEMBLY.md#8-a-coupled-completion-and-its-native-cost)
and [the stopped source-sharing ansatz](docs/SOURCE_REUSE_LIMITS.md#5-a-shared-conjugator-does-not-close-a-branching-fork).
The four modes can now be accessed by the compiled child and its actual
inverse. Rechecking low rank, basis access, or this child cancellation
does not advance the remaining precision ledger.

## First research pass

1. Give a changed native word that shares precision work through two adjacent
   depths. Either jointly compile the reduced ordered target wrappers or
   supply another complete rejected action with a proved sharing identity.
   The stopped common-conjugator word needs a changed construction. Specify
   actual source calls and unloading after logical-mode changes, on arbitrary
   dirty/reference inputs with at most two external clean flags.
2. Audit independently changed parent and children using an actual complex
   native coarse frame even when the target is real. If the contract
   survives, apply the same word at height three and state the band invariant.
   Check native sharing, not only dense matrix closure. Do not introduce
   a larger logical register that silently supplies initialized work.
3. Prove an induction for the expanded word, including any multiplicative
   growth in recursive child calls, and charging native programs,
   source preparation, masks, queries and unloading, literal controls,
   actual inverses, routing, amplification, and complete-word error.
   Normalization two and the declared workspace must survive composition.

The desired uniform ledger is

```math
T=O\!\left(\sum_v p_v+L\right),\qquad
\sum_v p_v=O(N),\qquad G=O(NL).
```

Here each local cost $p_v$ includes its native program and routing. A valid
linear-T circuit at the selected endpoint also settles the T-count question
with a larger, fully charged Clifford count. Constant global source uses
and linear tables are sufficient route invariants, not necessary properties
of every solution.

Stop a proposed word if its claimed induction hides repeated precision
charges, changed-address query unloading, rejected returns, or unavailable
clean work. Record a scoped failure in its existing proof home; it is not
a lower bound on all compilers. Further special-family results are not the
priority without a proved reduction for arbitrary updates.

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

The recovered checkpoint has 228 tests. Its GitHub checks passed for
Python 3.11, Python 3.13, and rendered presentation. The 2026-10-01 recovery
also reproduced all 228 tests locally on Python 3.12.14 / NumPy 2.3.5,
the reviewer walkthrough, all four fault-tolerant receipt suites, and the
offline synchronization check. The handoff edit passed 38 documentation
and narrative checks. Record fresh results when resuming; these results
do not certify later edits. Finite tests support the proofs and do not
establish asymptotic bounds.

The subsequent branching pass added 15 focused tests; all 243 local tests,
the reviewer walkthrough, and the offline synchronization check passed.
The new proof sections also received separate internal mathematical checks.
The subsequent repair pass supplies basis-free mode access, exact recursive
child cancellation, and a complete-word precision audit. The independent
local-wrapper emission still pays precision per depth; a jointly synthesized
precision ledger remains open.
This pass adds eight focused checks. The local full validation, final
repair suite, reviewer walkthrough, and offline synchronization check passed.
The fixtures check the operator identities and failure controls; the
native resource accounting remains an analytic argument.

Update this entry point when the task changes, while keeping detailed
research conclusions in [the existing checkpoint](docs/OPEN_PROBLEM.md)
and proofs in their existing chapters. Keep the frontier and the unresolved
step explicit so continuation does not depend on access to an old chat.
