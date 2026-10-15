# Joint orientation and recovery: session result

[Executed design](JOINT_ORIENTATION_RECOVERY_PROMPT.md) · [R85 proof](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) · [Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) · [Discovery notes](JOINT_ORIENTATION_DISCOVERY_NOTES.md)

**Outcome: the session's sub-square-logarithmic target is reached for
arbitrary nonlinear orientations of shared unordered depth magnitudes.**
The proof is [local result R85](GEOMETRIC_BAND_ORIENTATION_DEPTH.md).
The unrestricted arbitrary-magnitude frontier remains R79; logarithmic
depth for the orientation family remains a further question.

Execution began on 15 October 2026 at 22:04:30 UTC from main commit
`d6259c2`, on branch `research/joint-orientation-session`. The authorized
ceiling was five elapsed hours, ending at 03:04:30 UTC on 16 October.
The construction and independent proof reviews were completed within the
first hour. Repository validation and integration are recorded below;
nominal design allocations are not reported as elapsed work.

## 1. The new bound and its scope

Fix an integer $`L\ge6`$, let $`\eta=2^{-L}`$,
$`N=2^n`$, and $`\ell=\log_2(n+2)`$. R85 gives one coherent
Clifford+T circuit with one external clean qubit and simultaneous bounds

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\bigl(N/b^2+\ell\log_2(\ell+2)\bigr).
```

The literal dirty thresholds remain $`17(L+n+7)`$ for the real frame
and $`17(L+n+8)`$ for the prescribed phase-dressed complex magnitude
completion. Every node may choose either ordering of its depth's common
unordered magnitude pair. Orientations may be nonlinear; signs, zero
entries, balanced pairs, and supplied leaf phases are included.

At square-root dirty width this is log-times-log-log depth, strictly
below the previous square-logarithmic upper order for this family. The
norm includes all columns, literal phases, actual inverses, work leakage
and return, arbitrary dirty inputs, and reference entanglement.

This is a fixed-accuracy promised-family upper bound. It does not assert
optimal depth, a lower bound specific to this family, absolute constants
uniform in precision, efficient classical preprocessing, a general native
emitter, or a compiler for arbitrary node-dependent magnitude pairs.
The principal general A–D claims and their immutable receipts remain
unchanged. Final manuscript drafting was not part of this session.

## 2. The change that made the construction work

The previous conversion interface demanded a complete in-place leaf
permutation with no clean output word. The new construction asks where
the source image itself makes such a word available.

Partition the depth labels into geometrically sized intervals. There
are only $`O(\log n)`$ intervals, and each interval's width is at most
half the suffix space available after unloading its upper endpoint.
The old source tag already identifies its interval; computing that
coarse label is an ordinary short classical query.

On that entire spectral band, the selected descendant unloader maps
the data to a shorter frame and a zero suffix. Exchange part of that
suffix with a fixed arbitrary dirty output word. The output is now
clean, and its old contents are parked in positions ignored by the
local permutation tables. Forward XOR, swap, and inverse XOR implement
the local triangular orientation permutation; swapping back returns
both the external dirty word and the exposed suffix.

The local homogeneous product reference then changes the fine depth
label into a computational last-one label. XOR it into the old tag,
reverse the local transport, and reload the descendant state. The data
returns with its full marker phase; only the coarse interval label
remains. Its spectral binary search takes $`O(\log\log n)`$ source
layers rather than $`O(\log n)`$. Each layer costs $`O(\log n)`$
T-depth, and the recursive terminal has the same final depth order.

The crucial zero-space assertion holds on the ideal reached band.
It is never imposed on actual leaked states or the inactive sector.
Full unitary extensions and actual inverses propagate all earlier error.
The complete recovery word, rather than each speculative factor alone,
is the unit whose action and return are proved.

## 3. What the independent reviews checked

The proof was reconstructed separately for geometry, native resources,
error propagation, input certification, and the final combined theorem.
A fresh native review followed those first audits. The material checks
were:

- Unloading the entire band, including coherent mixtures of markers,
  exposes the claimed suffix; this is an isometry statement, not a
  collection of independent column error estimates.
- The two permutation queries read the fixed data and bin registers and
  target a disjoint dirty word. Their tables ignore the parked dirty
  tail and have $`2^h\mathrm{poly}(n)`$ rows, not $`2^{2h}`$.
- Parallel query outputs have private non-Clifford supports. Selected
  swap banks use exact dirty-involution echoes, and selected rotations
  use one parallel native product with its actual inverse.
- The local permutation returns the external dirty output even on an
  occupied logical tail. Only its intended permutation action needs
  that tail to be zero.
- Fine recovery never reads the tag it erases. Coarse recovery never
  reads the coarse tag it erases. Both use the same actual inverse words.
- The early error is bounded by
  $`\delta_0+2\delta_U+2\delta_C+8k\delta_R`$ on the complete source
  image. Coarse reflection comparisons include their zero-source-work
  projector, including on occupied source work.
- A protected bank $`M=A_L\ell`$ fits all records and suppresses the
  fixed polynomial count and dirty-width overhead. Original-accuracy
  narrow-width and bounded-dimension fallbacks retain the literal
  supplied-width thresholds.
- Finite certified common-angle rounding and arbitrary finite orientation
  tables cover balanced and zero pairs without an equality oracle.
  Existing sign and prescribed-phase wrappers retain literal phases.

The proof home contains the complete argument. These internal reviews
are not a substitute for future external technical review.

## 4. Other retained findings and the reformulation path

The opening whole-operator lane produced an exact two-level mixed-basis
fusion. Its Clifford change of basis sends the uniform anchor to a Bell
state, so the isolated factorization does not by itself close under
recursion. Two visible orientation levels are already affine. The
[discovery notes](JOINT_ORIENTATION_DISCOVERY_NOTES.md) preserve the
identity and its reopening condition without claiming a new family bound.

The recovery lane first considered simultaneous common-data reflections,
a common flag transporter, and weighted tree-shift chains. A transporter
may mix vectors within a depth band; this freedom is valid. Fixed
nonlocal descendant jumps, however, do not erase the chain-start depth
for free. A batch of commuting ideal reflections also needs a physical
shared-data schedule. Neither unfinished interface is used as a premise
in R85.

The decisive reformulation combined two ideas: coarse depth intervals
and clean space exposed after unloading a reached band. Equal-sized
intervals with a recursive fine decoder retain the old leading depth.
Geometric intervals instead make direct local orientation conversion
possible, and remove the fine decoder entirely.

An independent audit also enlarged [R81's terminal-height statement](RECURSIVE_TERMINAL_DEPTH.md)
to every integer $`0\le M\le n`$, uniformly in M, with depth
$`O_L(N/b^2+\ell\log_2(M+2))`$ and unchanged count/width guarantees.
This initially suggested larger protected histories. It remains a useful
proved extension, but the final R85 construction needs only the original
logarithmic-height terminal case. At M equal to n the extension still
gives square-logarithmic depth; it is not a separate unrestricted
improvement.

The bounded XOR-shear reserve found no generic constant-round conversion.
The discovery notes preserve the fully dirty five-shear restriction and
an exact six-shear example escaping its exponent. Neither limits the
successful joint recovery word or unrestricted native circuits.

## 5. Verification and repository state

The focused geometric-band tests reconstruct frame columns independently
from path products and exercise the new recovery identities. They include
nonlinear orientations, singular pairs, literal sign and leaf-phase
wrappers, dirty/reference correlations, and occupied inactive records.
These are bounded fixtures; native source/query synthesis and asymptotic
cutoffs are proved analytically. Existing scientific receipts are not
regenerated.

At the local validation checkpoint, 15 October 2026 at 22:54 UTC
(about fifty minutes after the start), all 945 tests passed. The 56-test
focused mechanism, narrative, mathematics, presentation, and provenance
suite also passed. Compile checks, the reviewer walkthrough, real and
complex native examples, the residual-QBP example, both resource ledgers,
and offline upstream synchronization all passed. The four exact
fault-tolerant receipt suites reproduced their preserved scientific
fixtures.

The integration record is the reviewed pull request: it records the
final 147-document desktop/narrow render check, required remote analytic
checks on Python 3.11 and 3.13, the remote presentation job, the merged
commit, and actual elapsed session time. These gates apply to the reviewed
tree; an initial narrow-screen equation overflow was corrected before
publication. Construction and proof review finished within the first
hour; the five-hour ceiling is not a claim of five hours spent. Internal
proof audits and bounded tests do not replace external technical review.

## 6. Precise resumption question

Can the remaining coarse-band label be erased with a constant number of
native source layers, while retaining one clean qubit, the count budgets,
and the full initialized-isometry contract? That would remove the
early recovery comparison factor. Logarithmic total depth also requires
compatible $`O_L(\ell)`$ terminal and routing schedules: the current R81
terminal contributes $`O_L(\ell\log(\ell+2))`$, and the current routing
bound also exceeds $`O_L(\ell)`$.

The broader magnitude question has a separate concrete gap: after
unloading a band, can its local flag be transported cheaply when the
magnitude pair depends on the node, so that a single tensor-product
reference is no longer available? R85 identifies exactly where its
promise is used; it does not settle that extension.
