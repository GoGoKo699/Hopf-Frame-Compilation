# Continuing research workspace

This is the entry point for resuming work when a previous conversation or
execution workspace is unavailable. The proofs and research decisions live
in the repository; a conversation summary is only a retrieval aid.

The 2026-10-01 joint-source pass continues from main commit
`d5b0d19245d7949913ecfa5cf65f472794655003`. Check the current branch and
later commits before continuing. Complete eight- and sixteen-mode native
words now share a source body, but their fair precision-cost comparison
leaves the compiler frontier unchanged.

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
   [the shared source body](docs/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets),
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
is stable in operator norm; a general cheap native conversion remains open.

The [four-mode native benchmark](docs/ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark)
uses the standard magic basis to reduce a real four-mode target to two
one-qubit factors. At an unchanged k-bit address, the existing one-target
compiler gives $`T=O(2^k+L)`$ and $`G=O(2^kL)`$ with one clean flag
and $`b=L+n+7`$, where $`n=k+2`$. Each stage borrows the idle logical
target as its extra dirty helper; complete-isometry errors telescope
without resetting the flag. Both fixed basis changes cost fourteen
Clifford gates in the displayed word. A complex native coarse baseline
requires its separately charged actual inverse.

At eight modes the root becomes a controlled Bell-projector rotation.
It has four commuting Pauli factors, but it cannot be absorbed into
independent prefix and child factors: the exact product-gate distance is
$`2\sin(\theta/4)`$ for $`0\le\theta\le\pi/2`$. This restricts the
simple magic-basis extension, not joint controlled synthesis. Direct
Woodbury inversion of the recursive Cayley update reconstructs the
original fine target wrapper, so that implementation does not remove
the repeated precision charge.

The earlier source-width result remains useful: outer loaders and monotone
bridges cost $`2(m_{\max}-1)`$ T gates. Separately exposed transformed
masks and repeated complete syndrome renewal still cost precision. That
pass solved transport, not the joint program. The new Cayley coordinates
also do not yet improve the generic T-count.

The [joint native word](docs/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
uses the original tree basis. A fixed scalar conjugator independent of
the logical target permits exact cancellation across changing targets and
predicates. Two parity boundaries return an arbitrary borrowed signal,
with no intermediate reset. The declared source counts fall from 135 to
87 at eight modes and from 180 to 114 at sixteen. These are structural
word comparisons, not leading precision savings: after hoisting the
common loader, the fixed mask commutes with its tail, leaving only a
two-T seed. Both comparison words then have the same leading loader and
variable-mask term $`(40g+4)q`$ for g stages at source width q.

For local depth three or four and an unchanged k-bit prefix, the complete
ledger gives $`T=O(2^k+L)`$, $`G=O(2^kL)`$ at the selected two-clean
allocation. The clean wires are used only as borrowed signal and helper.
This fixed-depth construction is not a scalable gain: its ten programmed
mask appearances per layer still renew precision, and its constant-width
reservation does not extend to arbitrary depth without adjustment.

The next bounded experiment is **joint synthesis of the programmable
mask interiors**, starting with two changing-target stages of this
explicit word. Expand actual inverses and compare against the already
simplified two-T fixed-mask baseline. Seek a native identity or encoding
that changes the precision recurrence as more stages are added. Further
outer-loader or fixed-mask cancellation alone does not address that task.
Keep zero defects, independent branch angles, a complex native coarse
inverse, and complete signal/dirty ports in the small examples.

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

The native small-example baseline has 282 tests and passed Python
3.11/3.13 CI, all four exact-receipt suites, and rendered presentation.
This pass adds three joint-source checks. Their complete matrices have
dimensions 128 and 256, with deliberately coarse q=2 source precision.
They verify actual gate words, literal phases, changed predicates and
targets, occupied signal ports, helper return, and a complex coarse
inverse. Fully simplified emitted words enforce the same leading
precision charge in both comparisons. Fixed-mask cancellation is also
checked at q=2, 3, and 5.

Run the commands above for the current 285-test suite and the latest
verification result. CI records Python 3.11/3.13, exact-receipt, and
rendered-presentation results on the corresponding commit. Separate
internal reviews checked the algebra and workspace ledger; this is not
external peer review. Finite fixtures do not establish native asymptotic
resource bounds.

Update this entry point when the task changes, while keeping detailed
research conclusions in [the existing checkpoint](docs/OPEN_PROBLEM.md)
and proofs in their existing chapters. Keep the frontier and the unresolved
step explicit so continuation does not depend on access to an old chat.
