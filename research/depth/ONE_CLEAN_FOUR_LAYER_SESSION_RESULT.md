# One-clean four-layer session: outer misalignment and the remaining gap

[Executed design](ONE_CLEAN_FOUR_LAYER_PROMPT.md) ·
[Outer-overlap theorem R97](OUTER_OVERLAP_RIGIDITY.md) ·
[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer)

The session proves **R97**, a positive four-layer error floor controlled
by the misalignment of the two outer physical axis lists. It gives a
uniform floor for aligned intersections of arbitrary size. The
unrestricted one-clean five-layer obstruction remains open; no native
counterexample to that gap was found. Optimal general depth remains
the destination, and the all-precision upper theorem R96 is unchanged.

Execution began on **16 October 2026 at 13:02:47 UTC**, from main commit
`539ba460299c7ced94cf68752641bc3d2fcbb490`, tree
`5a2fa257d1ff32014d62a9d2a2fc924dd86c6ad4`, after design PR #163.
The planned qualification gate was 14:32:47 UTC and the five-hour ceiling
was 18:02:47 UTC, including integration. The work branch is
`research/one-clean-four-layer-session`. The research claim was frozen
early, at 13:50:31 UTC, after the strengthened theorem and the remaining
native-independence issue were independently identified. Five hours was
a ceiling, not time reported as spent. Final manuscript drafting remains
separate.

## 1. The complete new result

For an actual two-plus-two split with one initialized clean qubit, let
$`L_A,L_B`$ span the outer physical T axes and $`I=L_A\cap L_B`$.
Define

```math
\mu_A=\dim I-\#\{\text{individual A axes lying in }I\},
\qquad
\mu_B=\dim I-\#\{\text{individual B axes lying in }I\}.
```

The lists are independent, so the counts are dimensions. For
$`r=\min\{\mu_A,\mu_B\}`$, the complete-frame witness

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3)
```

has a positive initialized-isometry error floor depending only on r.
The [proof](OUTER_OVERLAP_RIGIDITY.md) gives its conservative constant
through $`k(r)=2r+(r+1)(2^r-1)`$, $`t(r)=8r+8k(r)`$, and R77's
finite-boundary constant at t(r). The floor permits arbitrary dirty
inputs and references, optimized common phase, unbounded inner clean
spreads, crossing rank and total intersection dimension.

When either physical list spans I using a subset of its individual
axes, r=0 and the constant is exactly R76's uniform Clifford-midpoint
floor. No matching list is needed on the other side. Equal outer spaces
with different bases are covered. Bounded total intersection is a weaker
sufficient condition.

An arbitrarily accurate four-layer sequence must have both misalignment
dimensions diverging in every admissible native representation. The
radical of the clean form $`\beta(i,j)=\langle i,J_Aj\rangle`$ must
also diverge. These add to R74's crossing-rank and R77's inner-spread
requirements. They do not exclude a sequence in which all grow.
This witness has constant depth at fixed error by compiling its one
active qubit. Even an unrestricted five-layer obstruction for it would
not be a matching logarithmic lower bound for R96.

## 2. What changed the proof

The designed joint inner-algebra description was completed: at most two
noncommuting qubits and 2048 joint blocks suffice, but an unbounded
commuting center remains. Exact native examples show that the data
operation can mix every character of that center. Bounded inner block
size therefore does not pass the uniformity gate by itself.

The useful constraint comes from the actual outer physical axes.
If a physical phase layer conjugates a Clifford to another Clifford,
that Clifford must signed-permute the individual independent axes.
The shared clean image restricts these permutations to one matching of
fixed points and disjoint pairs. At bounded overlap, all but finitely
many axes can be transferred or reduced to an invariant isotropic space.
Fixing boundedly many central characters makes the residual phase
Clifford, yielding the finite expansion needed by R77.

Simultaneous diagonal Clifford extraction then strengthens overlap to
misalignment. A phase-polynomial layer mapping a commuting Pauli family
to Paulis admits one diagonal Clifford with all the same signed actions.
Its remainder commutes with that family. This allows one outer layer
in the proof to be phase-polynomial while the other stays physical.
Peeling individually contained axes from either native list reduces the
intersection to its misalignment dimension. The analytical transformed
layers need not themselves be physical T-depth-one layers; the proof
uses R76's established phase-polynomial interface.

## 3. Exact scope checks

The [representation limits](FOUR_LAYER_REPRESENTATION_LIMITS.md) give
the analytic arguments, with one proof home for each fact.

| Proposed shortcut | Supported conclusion |
|---|---|
| Fix the shared center and apply a finite-block argument | The native R77 tensor family mixes all center characters. Data invariance does not follow from clean invariance. |
| Treat growing displayed ranks as intrinsic | That family has an exact three-layer rewrite with a Pauli clean midpoint. Its displayed complexity is removable. |
| Remove the clean qubit by exact same-space phase-polynomial descent | A three-layer four-wire circuit has occupied determinant -i, impossible for any finite Clifford/weighted-cubic word on the original three logical wires, even up to phase. One additional dirty wire is still insufficient. |
| Recover the outer axes from Paulis sent to Cliffords | A native non-Clifford midpoint in the third hierarchy level sends every Pauli to a Clifford, so that intrinsic set loses the relevant axes. |
| Substitute a one-sided Clifford gap | The one-sided gap is proved; general two-sided compression has an extra projection and remains outside its scope. |

None gives a native family with nonzero occupied defect tending to zero.
The determinant example has finite projective order and does not
challenge the approximation floor. Two additional dirty wires remove
its particular determinant obstruction, leaving width-enlarged descent
logically possible.

## 4. The precise unrestricted step

After individually contained axes are removed, the remaining intersection
can be noncoordinate in both lists. On its radical part, S transports a
common Pauli i to a Pauli p. The product-label $`h=i+p`$ commutes with
the outer layers and lies in the clean displacement images. Correlated
inner data cores fix this label; small occupied defect forces their sign
actions to agree, so the relative data action preserves this smaller
center. This does not make each data image preserve its character sectors,
and h can vanish while the unaligned intersection remains large.

Fixing center characters can identify independent axes, leaving **both**
outer descriptions as dependent phase-polynomial layers. The physical
normalizer and matching argument then cease to apply: four-element
orbits of dependent labels can sum to zero. R97 pays for this loss
through the misalignment dimension. No proof bounds its minimum over
all exact native rewrites.

The early qualification decision was to stop the original center-block
reduction as a proof of the unrestricted gap, retain R97, and preserve
two concrete reopening questions:

1. Can the correlated dependent phase-polynomial pair be controlled in
   operator norm without bounding its center or misalignment?
2. Can every shared-image native pair be rewritten with bounded
   misalignment while preserving the two-plus-two interface, perhaps
   with explicitly returned additional dirty work?

The R77 rewrite makes the second question concrete; it does not answer
it generally. A counterexample to a proposed reduction must still be
distinguished from an actual vanishing-defect native family.

## 5. The short constructive scout

The [Gauss-sum query note](GAUSS_SUM_QUERY_ENCODING.md) preserves a
different encoding. An arbitrary Q-row Boolean table is one inverse
entry of an affine binary matrix of dimension $`m=O(\sqrt Q)`$.
Its quadratic Gauss sum has a constant-T-depth native realization with
$`T=O(\sqrt Q),G=O(Q)`$ and $`O(\sqrt Q)`$ dirty width. The table
sign appears in a normalized trace and a matrix element of magnitude
$`2^{-m}`$, rather than as a scalar on arbitrary work.

Every address has both ±1 eigenspaces. A clean-only spectral extractor
cannot distinguish their multiplicities, giving full-operator error at
least $`\sqrt2`$ for a nonconstant table. An address-dependent
transformation of the Gaussian work is the unpaid step. The scout ended
within its 30-minute cap; no faster dirty query or complete frame
compiler qualified.

Two primary-source scope checks supplied no missing theorem:
[Bao et al., arXiv:2605.26983v1, §6.1](https://arxiv.org/html/2605.26983v1)
leaves a separation question for products of six third-level operations
beyond its proved inverse theorem;
[Sinha–Lovitz, arXiv:2610.06819v1](https://arxiv.org/html/2610.06819v1)
provides a Clifford rank/coefficient estimate at fixed actual rank, while
the unrestricted outer dressing has no such bound. Neither is a new
premise of R97.

## 6. Review, evidence and integration

Construction and skeptical lanes independently reconstructed the normalizer,
matching constraints, mixed-image count, invariant isotropic space,
tensor quotient and finite expansion. A fresh audit checked asymmetric
extraction, misalignment, ring phases and the witness consequence.
The complete canonical R97 proof passed a separate end-to-end review.
The data images are never assumed to preserve sectors used solely to
expand the midpoint.

The [nine exact fixtures](../../tests/test_outer_overlap_rigidity.py)
check native-axis compatibility and routing, a non-Clifford transverse
midpoint and its gauge, one-clean parity encoding, complete center mixing,
the alternative three-layer word with its scalar, and every occupied
entry, leakage and determinant in the descent example. They reuse
retained exact arithmetic. Infinite-width gaps and uniform constants
are analytic conclusions.

Local validation passed all **1018 tests** in 374.248 seconds, including
the nine new fixtures. All 50 focused documentation checks pass.
Compilation, the reviewer walkthrough, real and complex native examples,
the q=16 residual example, both n=12 resource ledgers and offline upstream
verification pass. All four retained scientific receipts match without
regeneration.

Rendered checks cover all 14 changed Markdown pages at desktop and narrow
widths, with SVG and native math; there are no layout failures or TeX
errors. All five explanatory diagrams pass at both widths. Representative
proof and session previews were visually inspected. The final whitespace
check passes. The associated pull request records the exact reviewed tree,
remote CI and merge receipt.

## 7. Continuing destination

R97 narrows a constant-layer problem and identifies the remaining
algebraic interface. It does not settle optimal general depth. R96
remains the all-precision upper theorem, R79 remains selected manuscript
Result D, and final manuscript drafting remains separate. The next
mathematical step should address the unbounded dependent-axis pair,
prove a useful native rewrite, or introduce a dimension-sensitive
witness and invariant. Repeating finite tests or expanding the same
center with a larger parameter does not answer those questions.
