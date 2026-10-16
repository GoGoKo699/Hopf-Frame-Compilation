# Interval-recovery factorization: session result

[Executed design](INTERVAL_RECOVERY_FACTORIZATION_PROMPT.md) ·
[R87 coin-dilation theorem](COIN_DILATION_DEPTH.md) ·
[Independent weighted-chain factorization](WEIGHTED_CHAIN_FACTORIZATION.md) ·
[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer)

**Outcome: the preferred general fixed-accuracy depth improvement is
proved.** Product-coin dilation makes the fine spectral depth a
classical label after Hadamards, removing the repeated interval-overlap
recovery. The complete compiler is
[R87](COIN_DILATION_DEPTH.md). An independent complete weighted-chain
factorization also resolves the design's central mathematical scheduling
question; its native implementation remains separately scoped.

Execution began on 16 October 2026 at 01:19:57 UTC from main commit
`8cfe25631398074b0563c0666bf5c79440c4253d`, on branch
`research/interval-recovery-factorization`. The authorized ceiling
was 06:19:57 UTC. Construction, internal review, validation, and integration
are distinct checkpoints. The five-hour allowance is not a claim that
five hours were spent.

## 1. The new general bound

For every fixed integer $`L\ge6`$, put $`N=2^n`$,
$`\eta=2^{-L}`$, and $`\ell=\log_2(n+2)`$. One circuit with
one external clean qubit has

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell\log_2(\ell+2)\right).
```

The literal dirty thresholds remain $`17(L+n+7)`$ for real frames
and $`17(L+n+8)`$ for prescribed complex magnitude completions.
Arbitrary node-dependent magnitudes, signs, zeros, and supplied leaf
phases are included. The full initialized-isometry norm covers all
columns, literal phase, leakage and work return, dirty inputs, and
external references.

At square-root dirty width this improves R86's general
$`O_L(\ell^{3/2})`$ bound to $`O_L(\ell\log(\ell+2))`$.
It reaches R85's depth order for the general family. R79 remains the
absolute all-precision theorem. Optimal general depth, a
dimension-growing unrestricted lower bound, and a precision-uniform
extension of this gain remain open. No efficient preprocessing,
scalable native emitter, or final manuscript is asserted.

## 2. The representation change that removes recovery

The original interval dilation assigns consecutive atom intervals to
leaves. Its balanced-Haar overlap is an isometry, but taking entrywise
absolute values loses interference and requires square-root
normalization. The opening admissible-spine witness confirms that
this loss cannot be removed merely by better positive weights for
the unchanged overlap.

The new representation assigns one independent tau-bit coin word to
each local tree level. At each node an integer threshold determines
the branch. A logical leaf corresponds to the product rectangle of
coin strings consistent with its path. These rectangles partition the
whole coin space. A generic paid family source prepares the normalized
rectangle while preserving the leaf; a paid classical coin-to-leaf
lookup then clears that old address.

Every scaling vector becomes a normalized prefix rectangle with uniform
later coins. Every depth-j marker becomes that prefix rectangle, a
normalized mean-zero threshold vector in coin j, and uniform later
coins. Tensor Hadamards send the later coins to zero and leave the
j-th coin with no zero component. The last nonzero coin block is
therefore exactly j, for every column and every superposition in its
depth band.

The complete local word is the embedding, Hadamards, a fine-tag XOR,
Hadamards, and the actual inverse embedding. It returns the original
data with literal phase +1. A single source-isometry error delta gives
full error at most $`2\delta`$, including rejected work. The image
is contained in the corresponding coin sector; it need not fill it.
No reflection about that reached image is supplied or assumed.

The local depth is $`O_L(\ell)`$. Its largest table address remains
$`a+q+r`$, with $`r=q\tau`$, so the existing thin bins, selected
dirty routing, protected count margin, coarse search, and terminal
compiler apply. Coarse recovery and terminal processing give the
remaining $`O_L(\ell\log(\ell+2))`$ term.

## 3. Independent factorization also completed

The design's whole-comb estimate was reconstructed through a finite
weighted maximal inequality, including the residual terminal leaf.
The absolute matrix norm is at most $`2+2\sqrt2\lt5`$,
independent of chain length and weights.

The global scheduling question admits a different answer from the
proposed boundary repairs. Decompose the complete unit-atom tree into
ordinary heavy paths, each continued to a unit atom. Sort every path's
side masses in descending order. A finite sorting lemma preserves the
lower child-fraction bound and ensures that each peel has at least its
own mass remaining. Greedy packets then have no unresolved trailing
subtree. One pass halves logarithmic imbalance up to a fixed additive
constant.

Complete path transitions fix their scaling coordinates and act on
disjoint marker-coordinate blocks, even when physical subtree supports
nest. Thus a whole pass has bounded absolute norm. Iteration takes
$`O(\log(\tau+2))`$ passes. Intermediate subtrees may be arbitrary
atom sets; one final depth-first atom permutation restores intervals
before the bounded-bias terminal comparison.

The [factorization proof](WEIGHTED_CHAIN_FACTORIZATION.md) is retained
as an independent exact representation result. Its square-factor source
wrappers require two atom indices and the corresponding native ledger.
It is not presented as a second unconditionally compiled theorem.
R87 uses the simpler direct coin readout.

## 4. What the original partial-transport reserve missed

The design fixed an initial reached frame $`V_0=B_rR`$ and an
interval-embedded final frame $`V_s=EW`$. A chain
$`T_j=V_jV_{j-1}^\dagger`$ necessarily has product
$`V_sV_0^\dagger`$. That endpoint can itself have much larger
absolute norm than the primary coefficient overlap.

For clarity, in the repository's marker ordering
$`R|x\rangle=|0^{r-q},x\rangle`$. The nonroot columns of
$`B_rR`$ are the q-level Haar markers on the first $`2^q`$
atom coordinates; its root is uniform on all $`2^r`$ coordinates.
Take a homogeneous local tree with left probability p and right
probability $`1-p`$. Each marker has a distinct near-corner leaf
whose coefficient magnitude is at least $`(1-p)^{q/2}`$:
for a marker at prefix v of length $`j-1`$, choose
$`x=v0\,1^{q-j}`$. Its last zero uniquely identifies the marker.

For marker i with balanced support size $`L_i`$, test that leaf row
of $`WR^\dagger B_r^\dagger`$ against the balanced marker's sign
vector, which is bounded by one and equals $`\sqrt{L_i}`$ times
the marker. Orthogonality extracts its W coefficient, so the row's
absolute sum on the first atom block is at least
$`\sqrt{L_i}(1-p)^{q/2}`$. Test the absolute matrix against the
uniform unit vector on that block and sum over these distinct rows.
Since $`\sum_iL_i/2^q=q`$,

```math
\bigl\|\,|EW R^\dagger B_r^\dagger|\,\bigr\|
\ge\sqrt q\,(1-p)^{q/2}.
```

Positive disjoint E preserves this absolute norm. When qp is bounded,
constant-absolute-norm factors for these fixed endpoints require
$`\Omega(\log q)`$ stages. This is a restriction on that exact
transport interface, not a native-depth lower bound. The coin route
changes the embedding and reads the needed flag directly, so it does
not use that factorization.

## 5. Independent reviews and retained safeguards

The opening comb and sharpness calculations, sorted heavy-path pass,
product-coin identity, native circuit, and complete outer composition
were reviewed independently. The primary written proof received
separate full-chapter mathematical and native-resource reconstruction.
These are internal technical reviews, not external peer review.

The written review made the active-sector domain explicit:
$`\|UJ-\overline E\|\le\delta_E`$ is an H-one source statement.
The H-zero word is proved separately as exact identity on occupied
registers. It also retained the following points.

- Sector containment is sufficient; no surjectivity onto all coin
  strings of a given last-block depth is claimed.
- Membership is a total classical function. Its actual inverse
  reconstructs the leaf before reversing the conditional source.
- Hadamards on the entire fixed coin bank cancel on unused arbitrary
  padding; no padded-zero test or extra clean bank is assumed.
- Shared query targets are disjoint from the union of relevant public
  addresses, with private simultaneous dirty pools fully charged.
- All approximation estimates use actual inverses and ideal-image
  hybrids. Leakage is propagated, never reset.
- The original-accuracy narrow-width and finite-dimension fallbacks
  preserve the literal real/complex reservations.

## 6. Validation and integration record

The new [seven bounded fixtures](../../tests/test_coin_dilation_depth.py)
cover heterogeneous rectangles, all scaling flags and literal markers,
exact last-block support, physical source and membership return,
fine-tag erasure with references, nonzero source-work leakage, and
inactive occupied-work cancellation. They pass locally. Their bounded
matrix completions are interface checks, not native circuit emitters.

Local validation passed:

- compileall and all 959 tests, including the seven new fixtures;
- the reviewer walkthrough and real, complex, and q=16 residual native
  examples;
- all four immutable fault-tolerant scientific receipt suites, unchanged;
- offline upstream provenance and both n=12 resource ledgers;
- 50 focused narrative, presentation, math, and literature-policy checks;
- all 155 rendered Markdown pages, with no reported failures, and
  visual inspection of the new proof pages at desktop and narrow widths.

The associated pull request records the reviewed tree, required CI
conclusions, integration commit, and final elapsed time. The merge gate
is successful analytic validation on Python 3.11 and 3.13 and the
rendered-presentation job on that reviewed tree. Internal review and
finite tests do not substitute for external technical review.

The next general-depth question is to improve the remaining complete
coarse-recovery, routing, and terminal costs together while preserving
the simultaneous count and one-clean return contract. The local
spectral-readout problem has a direct solution in the coin
representation; optimal general depth remains the destination.
