# Continuing research workspace

This is the entry point when a previous conversation or execution
workspace is unavailable. Proofs and research decisions live in the
repository; conversation summaries are retrieval aids.

The 2026-10-01 revision starts from main commit
`7856bfa8bf907dcb12d7e856387784bf6af6f27e`. Check the current branch
and later commits before continuing. The recent small-example pass is
complete. It supplied useful native identities but did not narrow the
generic endpoint gap. The next task is a specific compatibility test at
the best grouped compiler's interface, not further generic mask search.

## Mandate and reading order

Continue theorem-led research in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge are authorized for this repository.
Keep manuscript writing and release work outside this research pass.
Use small analytic examples and finite checks; no large simulations,
QRAM, resets, supplied catalysts, or hidden initialized work.

1. Read the [current assessment](docs/OPEN_PROBLEM.md#current-assessment-what-the-results-establish)
   and [revised next task](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule).
2. Read [conditional-suffix Sections 3--7](docs/CONDITIONAL_SUFFIX_COMPILER.md#3-a-small-coefficient-table-and-a-streamed-coarse-circuit)
   and its one-clean extension. These give the actual best compiler's
   scalar/atom interface, private work, error budget, and group sums.
3. Read the [shared source body](docs/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
   and [source-width analysis](docs/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation).
   Their different source families and scope must remain explicit.
4. Use [verification](docs/VERIFICATION.md), [attribution](docs/SOURCE_MAP.md),
   and [related work](docs/RELATED_WORK.md) before extending a claim.
   The [overview](README.md), [technical narrative](REVIEW.md), and
   [publication scope](manuscript/PUBLICATION_SCOPE.md) retain the main story.

## What is established

| Result | Capability and boundary |
|---|---|
| Exact compilation | Matching size and CNOT/depth tradeoffs for the prescribed complete frame, for every clean-work budget |
| Fault-tolerant frontier | Matching T-count under the stated sufficient-clean reservation; complete-frame and fixed-parameter QBP guarantees stand |
| Constant-clean compiler | One clean qubit gives $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ at $`b\ge L+n+7`$; T-depth is a separate open problem |
| Promised updates | Antichain and sparse ancestor-closed changes have $`O(N+L)`$ T-count under their distinct assumptions; generic rounding supplies neither promise |
| Compact residuals and transport | Linear classical Cayley/weighted data, complete coupled boundaries, and cheap source-width transitions; generic coherent conversion and interior programming remain charged |
| Small native programs | Fixed-address quaternion compression, four-mode magic-basis factors, and eight-/sixteen-mode changing-target words with full borrowed-signal return; fixed-size cost savings do not improve the generic bound |

For the selected complete real-frame endpoint,

```math
N=2^n,\qquad a=2,\qquad b=N+n+7,\qquad L=N,\qquad n\ge3,
```

```math
\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n)),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

The upper bound also uses only one clean qubit. Recent passes have not
narrowed this gap or established that its resolution is close. The
established publication results do not depend on closing it.

The latest shared word reduces declared source appearances from 135 to
87 at eight modes and from 180 to 114 at sixteen. Hoisting the loader
and commuting the fixed mask through its tail leaves both comparison
words the same leading precision term $`(40g+4)q`$, for g ordinary
stages at source width q. The remaining savings concern fixed scalar
work. Its precision charge still grows with the number of stages.

This matters strategically: the best compiler already replaces n ordinary
layers by $`R=O(\ell_*(n))`$ unequal groups while keeping total table
work linear. A better constant in a layerwise word does not remove that
remaining group factor. The paired-source real-Y circuit has not yet
been transferred to the grouped forward/reverse scalar SELECT, its
complex phases, distinct rejection flags, and private-suffix reflection.

## Next bounded task

First audit a concrete canonical scalar completion in the actual grouped
interface. For each retained coefficient $`c\in[0,1/4]`$, propose
$`R_c=\exp[-i\arccos(c)Y_\sigma]`$ on the existing scalar rejection
flag sigma. Its accepted entry is c. Keep the atom flag separate, retain
literal phases, and use actual inverses for reverse atoms.

Use the allowed two-clean layout with external h and sigma. A new native
rotation signal would be borrowed from the dirty pool and excluded from
the initialized-work reflection. Keep sigma out of temporary buffers.
The first deliverable is a symbolic forward/inverse atom and a register-
lifetime table that audits every intervening SELECT, atom, coarse, and
reflection word against $`Z_\sigma`$. This is a compatibility question,
not a certified replacement theorem or an elementary compiler project.

Use fully amplified inner rotations and derive a fresh full-operator
error budget, with the constant width increase $`q_g=m_g+O(1)`$
charged. A zero coefficient still gives a nonidentity rotation.
Name the actual paired source: cheap bridges for the retained one-tail
loader do not automatically transfer. The detailed
[decision](docs/OPEN_PROBLEM.md#next-bounded-task-and-stopping-rule)
records these checks and the live-work contract.

Compatibility alone is not an endpoint gain. Before proceeding to two
unequal groups, specify an actual fusion rule and its variable-R cost
hypothesis. For the retained sufficient route, seek

```math
T_{\rm joint}\le c_RL+O\!\left(\sum_gQ_g+RP(n)\right),
```

with $`c_R`$ bounded independently of R and P a fixed polynomial
independent of L. Table and coarse-program sizes must combine additively.
If a proposed rule closes for two groups, test that same rule on a third
and derive the recurrence. Intermediate group action may be deferred,
but final decoding, changed predicates, released dirty tails, literal
inverses, and all rejected/reference action remain charged.

If the expanded word still renews precision per group, close that
candidate. Do not continue with larger mask sweeps or another fixed-size
cancellation. Other complete representations remain allowed; this budget
is sufficient, not a lower bound or compulsory architecture. A focused
literature check should address a concrete alternative mechanism.

## Restore and verify

Run from the repository root with `requirements.txt` installed:

```bash
git status --short --branch
git rev-parse HEAD
python scripts/reviewer_walkthrough.py
python validate.py --quiet
python scripts/verify_fault_tolerant.py
python scripts/check_upstream_sync.py --offline
```

The merged native checkpoint has 285 passing tests, Python 3.11/3.13 CI,
all four exact-receipt suites, and rendered presentation. Its newest
complete matrices have dimensions 128 and 256 at deliberately coarse
q=2 precision. They check actual gate words, complete dirty/signal ports,
helper return, and fair simplified counts; they do not certify asymptotic
costs. Internal algebra and resource reviews are not external peer review.

This revision changes the research decision and navigation, not the
proofs or executable construction. Keep detailed conclusions in
[OPEN_PROBLEM.md](docs/OPEN_PROBLEM.md) and proofs in their existing
chapters so continuation does not depend on an old chat.
