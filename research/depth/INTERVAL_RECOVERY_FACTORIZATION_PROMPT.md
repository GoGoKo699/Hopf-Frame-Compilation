# Five-hour design: factor the interval recovery

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) · [R86 proof](INTERVAL_DILATION_DEPTH.md) · [Previous outcome](GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md) · [Selected routing](SELECTED_DIRTY_ROUTING.md)

**Status: executed; see the [session outcome](INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md).**
The design baseline is main commit
`4f0fc7347965991230e3aaec95183f597dc02223`, containing R86.
The next execution may use at most five elapsed hours, including review,
validation, waits, and integration. Its clock starts when execution
begins. This design's preliminary calculations were opening checks,
not themselves a new general-depth theorem. Execution produced the
separate [R87 coin-dilation proof](COIN_DILATION_DEPTH.md) and
[weighted-chain factorization](WEIGHTED_CHAIN_FACTORIZATION.md).

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge remain authorized. Preserve the distinction between research
and final manuscript drafting.

## 1. Target and the new bottleneck

Put $`N=2^n`$ and $`\ell=\log_2(n+2)`$. Start at fixed
$`\eta=1/64`$, with one external initialized qubit and
$`b=\lceil C_b\sqrt N\rceil`$ arbitrary dirty qubits for a sufficiently
large absolute constant. The preferred target is a general compiler with

```math
T=O(\sqrt N),\qquad G=O(N),\qquad
D_T=O(\ell\log(\ell+2)),
```

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta.
```

W is the prescribed complete real frame with arbitrary node-dependent
normalized pairs, signs and zeros, or its prescribed complex magnitude
completion with supplied leaf phases. A completed statement must cover
every column, literal phase, leakage and work return, arbitrary dirty
input, and external references. No QRAM, measurement, reset, free
coherent evaluator, supplied precision state, or extra external clean
word is allowed. Elementary Clifford counts remain charged. Every table
and stage choice must admit terminating certified classical preprocessing;
its running time need not be efficient.

Any strict general improvement $`D_T=o(\ell^{3/2})`$ with these
simultaneous budgets is a successful theorem outcome. Matching the
preferred bound is the primary aim. After the large-width proof, restore
the original literal real/complex width thresholds through a proved
splice. An absolute all-precision extension and optimal depth are
separate later questions; neither replaces this session's target.

R86 rounds once with $`\tau=\Theta_L(\ell)`$, embeds a q-level local
tree in $`r=q\tau`$ atom bits, and uses
$`K=B^\dagger EW`$. Its absolute-overlap normalization is
$`O(\sqrt\tau)`$. Each recovery round costs $`O_L(\ell)`$,
giving $`O_L(\ell^{3/2})`$ local depth. Coarse recovery and the
terminal already fit $`O_L(\ell\log(\ell+2))`$; boundary routing
also fits this target.

The repeated cost includes sources, balanced transforms, membership
operations, and anchor reflections. Saving only one query subroutine
does not save the complete round. The proposed change is to replace
one large-normalization conversion by a short sequence of paid
conversions, retaining interference between their entries.

## 2. Opening diagnostic: unchanged positive weights are insufficient

Reconstruct this admissible family before choosing a normalization route.
Take integers $`s\ge4`$, $`q\ge2^s-1`$, and $`\tau\ge s`$.
On the first $`m=2^s-1`$ nodes of the all-right spine, use left
probability $`p=2^{-s}`$ and right probability $`1-p`$.
Set every other node to a balanced pair. These probabilities are exactly
dyadic at denominator $`2^\tau`$, and the complete tree has positive
integer atom masses at $`r=q\tau`$.

The peeled intervals have lengths $`p(1-p)^t\in[p/4,p]`$ in the
unit-interval normalization, and cover $`[0,1/2]`$. A dyadic J
inside that half with length $`u\ge8p`$ has at least
$`u/(4p)`$ whole peels in either designated half of J. For peels
in its first half, direct integration of the two Haar wavelets gives
overlap magnitude at least $`\sqrt{p/(8u)}`$. Testing the absolute
overlap on the unit vector constant on the m spine-marker columns gives

```math
\|\,|B^\dagger EW|\,\|\ge\frac{\sqrt{s-3}}{16}.
```

Check the sign convention, complete-peel counting, and atom integrality.
For compatibility with R86, one can take $`n=2^{2s}`$:
$`\tau=4s+2L+12`$, and the largest thin bins eventually accommodate
the spine. Verify that placement rather than treating q as independent
of the outer allocation.

Two independent planning calculations agree on this witness. If its
opening reconstruction confirms it, do not spend the session optimizing
positive Schur weights for the unchanged K. The
[entry-indexed normalization lemma](ENTRY_INDEXED_ENCODING_LIMIT.md#1-the-encoding-contract)
also applies to rectangular matrices and prevents entry-local channels
from beating their absolute norm. This is an exact representation
obstruction, not an approximation-robust claim about every changed
encoding and not a native-depth lower bound. A different reference,
cross-entry interference, or a different reached-space word remains open.

## 3. First constructive lemma: a whole weighted chain

Ordinary tree rotations can require too many serial rounds. The proposed
primitive handles an entire comb at once.

Let positive weights $`\mu_0,\ldots,\mu_m`$ sum to one, including
the residual terminal leaf $`\mu_m`$. Put
$`M_i=\sum_{j\ge i}\mu_j`$. The comb frame H has root column
$`(\sqrt{\mu_j})_j`$ and marker columns, for $`i\lt m`$,

```math
w_i(i)=-\sqrt{\frac{M_{i+1}}{M_i}},\qquad
w_i(j)=\sqrt{\frac{\mu_i\mu_j}{M_iM_{i+1}}}\quad(j>i),
```

with zero entries for $`j\lt i`$. Derive its complete orthogonality,
root action, and literal marker signs. The proposed estimate is

```math
\|\,|H|\,\|\le2+2\sqrt2\lt 5.
```

There is a concrete proof to check. Split the absolute matrix into its
root column, a diagonal contraction, and a positive tail operator S.
Writing the latter in weighted tail-average coordinates gives

```math
\lambda_i=\frac{\mu_iM_{i+1}}{M_i},\qquad
\sum_{j=i}^{m-1}\lambda_j\le2M_{i+1}\quad(0\le i\lt m).
```

The finite weighted maximal/Carleson estimate would then give
$`\|S\|\le2\sqrt2`$. Prove that finite inequality with these exact
tail and endpoint conventions, or reconstruct the relevant finite
argument in [R82](BOUNDED_DENSITY_DEPTH.md#4-a-dimension-independent-entrywise-positive-correction-norm).
Do not treat an unverified named inequality as an imported oracle.

If this lemma passes, changing a comb into an outer comb of inner
packet-combs uses at most three constant-normalization factors.
The complete coefficient transition fixes the common scaling coordinate;
its auxiliary factors need not individually do so. Batch the complete
transitions on disjoint affected marker-coordinate blocks.
A conservative bound $`5^3`$ is enough:
constant optimization is not the research objective.

Compile the complete coefficient transition with the paid row/column
source and square-swap wrapper. Include a residual terminal leaf,
arbitrary packet count, identity on unused labels, literal signs,
and guarded inactive action. A source for one column is insufficient.
The fact that packet scaling vectors are abstract quantum vectors
does not supply an uncharged physical change of basis.

## 4. Primary lane: halve logarithmic imbalance

Complete the atom embedding's columns EW to a full atom interval-Haar
basis F. Within each original leaf interval, add the Haar wavelets from
integer midpoint bisection. Label the original columns so that

```math
FR=EW,\qquad R|x\rangle=|0^{r-q},x\rangle.
```

The added splits have child fractions at least one third. Thus the
complete interval tree starts with logarithmic imbalance

```math
\beta=\max_I\log_2\frac{|I|}{\min(|I_0|,|I_1|)}
\le\max\{\tau,\log_2 3\}.
```

The central unproved lemma is:

> Construct a new interval-Haar basis with
> $`\beta'\le\beta/2+C`$, using a bounded number of whole-comb
> reassociations and explicitly paid repairs, with constants independent
> of q, r, the number of nodes, and the masses.

One specific construction attempt sets
$`\epsilon=2^{-\lceil\beta/2\rceil}`$, follows maximal chains
whose continuing child has fraction greater than $`1-\epsilon`$,
and groups consecutive side intervals into packets of mass between
$`\epsilon P`$ and $`2\epsilon P`$, where P is the mass at
the packet's start. Reassociate into an outer comb of packet-combs.
Derive the actual child-ratio estimates; do not assume that local
packet estimates automatically control the resulting whole tree.

The first serious obstacle is a trailing packet of total mass below
$`\epsilon P`$. A tiny child followed by two large balanced children
already prevents a proof using only reassociation of whole children.
Try extracting a contiguous boundary chunk from the terminal heavy
subtree, along its boundary spine, through the same comb primitive.
Keep atom lengths integral; a proposed fractional boundary needs an
explicit rounding/repair bound, or an exact dyadic alternative.

The decisive global question is whether those boundary repairs can be
scheduled in a bounded number of disjoint marker-coordinate batches.
Nested repairs may share a boundary spine. Write their dependency and
conflict relation explicitly. Disjoint interval supports are sufficient
but not necessary; nested moves must prove that shared scaling vectors
are fixed. A local repair picture is not a global parallel schedule.

Let $`S(\beta)`$ count serial paid conversion layers after all actual
inverses and repairs. The preferred recurrence is

```math
S(\beta)\le S(\beta/2+C)+O(1),
\qquad S(\beta)=O(\log(\beta+2)).
```

Stop the imbalance-reduction passes at a fixed threshold
$`\beta_0>2C`$, then use the bounded-imbalance terminal conversion.
The recurrence is not iterated below that threshold.

A proved $`O(\log\beta)`$ cost per halving pass would instead give
$`S(\beta)=O(\log^2(\beta+2))`$, still a strict improvement over
$`\sqrt\beta`$. Retain such an intermediate general gain.
Two serial recursive children must be summed; a maximum requires a
single paid block-diagonal implementation with its simultaneous work
charged. Ordinary $`O(\log q)`$ rotation rounds give the wrong
parameter and can restore square-logarithmic depth.

For bases $`H_0=F,H_1,\ldots,H_s`$ ending at bounded imbalance,
write the actual coefficient transitions

```math
Q_j=H_{j+1}^{\dagger}H_j,\qquad Q_{\rm end}=B^\dagger H_s.
```

R86's norm argument gives constant normalization for $`Q_{\rm end}`$:
it uses interval supports, orthonormality, ancestor telescoping, and the
child-ratio bound, without requiring all leaves at one tree depth.
The exact required identity is

```math
W=E^\dagger BQ_{\rm end}Q_{s-1}\cdots Q_0R.
```

Every Q needs either a proved constant absolute norm and paid source
implementation, or an explicit cheaper structured native word.
The intermediate H matrices themselves are not free gates.
Within an all-constant-normalization chain, the opening witness forces
$`S=\Omega(\log\tau)`$ by submultiplicativity of absolute matrices.
That explains the chosen stage target within this model only.
An explicitly cheap factor with larger absolute norm is also allowed.

## 5. Reserve lane: control only the reached atom space

If full atom completion is the obstruction, weaken that requirement.
Seek explicitly specified isometries

```math
V_0=BR,\qquad V_s=EW,\qquad
T_j=V_jV_{j-1}^{\dagger}.
```

Then $`T_s\cdots T_1V_0=V_s`$. The physical order in this reserve
is R, B, $`T_1,\ldots,T_s`$, and the actual $`E^\dagger`$.
It differs from the full-basis order in Section 6. A useful transport must have a
certified constant absolute norm or an independently paid native word.
Naming a changing subspace is not a construction, and
$`V_j^\dagger V_{j-1}`$ is generally a contraction rather than the
needed transport.

Prove a reached-isometry amplification lemma before using this lane.
The nonzero singular values of each ideal $`T_j`$ are one.
Ordinary all-work-anchor amplification can therefore act correctly on
the reached input range, without a supplied reflection about that range.
Derive the literal phase and the full initialized-output error from the
actual source word. In the block repair, the rejection Gram matrix for
$`\zeta T_j`$ remains bounded below when $`\zeta\le1/2`$.
Check uniform perturbation, leakage propagation, guards, and every
tested work field; this is an obligation, not a quoted new lemma.

Moving positive disjoint interval embeddings alone cannot solve the
transport problem. If $`V_j=E_jW_j`$, then

```math
|T_j|=E_j|W_jW_{j-1}^{\dagger}|E_{j-1}^{\dagger}.
```

Its norm equals that of the middle absolute matrix. Dropping that middle
factor would drop the frame change. The reserve must identify a genuine
new signed transport or return identity.

Generic prefix caching is not a third lane. R86's row, atom, and source
indices change through recovery; only the outer prefix, bin, and enable
persist. The retained [dirty-cache failures](JOINT_DEPTH_SESSION_RESULT.md#5-a-generic-dirty-cache-correction-already-contains-a-lookup)
and [source-monitor failures](COMMON_SOURCE_REUSE.md#11-a-common-source-identity-and-the-remaining-reflection)
must not be renamed as savings. A different clean-history proposal is
admissible only after stating its stored information, live capacity,
dynamic update, reflection, and complete return word.

## 6. Native allocation and success gates

Full atom transitions now use addresses of length $`a+2r`$,
not R86's $`a+q+r`$. For the full-basis lane, a sufficient fixed
layout is the q-bit local input X, an r-bit atom bank Z, an r-bit
auxiliary bank Y, and separately reserved short source/flag metadata.
Move X into the designated R coordinates of zero Z. Each square
transition uses row/column sources and a final swap between the r-bit
banks. Apply B and the actual inverse embedding only at the ends.
The initial X is zero when that inverse writes the final local frame
back to X.

With $`v=h-d+m_0`$, $`d=a+q`$, and $`r=q\tau`$, use, for
example,

```math
q=\min\left\{d,\left\lfloor\frac{v}{4(\tau+2)}\right\rfloor\right\}.
```

Then $`2r\le v/2`$ and $`a+2r\le h+m_0`$. Fix
$`m_0=\lceil c_L\ell\rceil\ge8(\tau+2)`$ and all metadata before increasing the independent
protected margin $`M-m_0`$. Recompute rounding at the first bins,
the finite-dimension fallback, and any extra bank required by a repair.
This constant-factor thinning still permits $`O_L(\ell^2)`$ bins
and short coarse search.

Keep one selected physical layout around the entire local forward/inverse
erasure episode. Use active-only source and amplification anchors,
excluding parked arbitrary dirty values and inactive padding. Compile
the ragged family of actual addresses with the
[selected-query lemma](SELECTED_DIRTY_ROUTING.md#2-parallel-selected-tables-with-different-relevant-addresses).
Never replace its summed table volume by an unpriced maximum-width
union of physical banks.

For S serial stages, allocate complete local stage errors summing to the
chosen local tolerance; $`O(\eta/S)`$ per stage is sufficient.
Reinsert the resulting precision into source depth and metadata.
Actual banks can be reused by initialized-isometry hybrid arguments;
do not assume perfect intermediate return or silently reset them.

The desired full local-family envelope, for some fixed polynomial
degree chosen before the protected margin, is

```math
T_{\rm loc},w_{\rm loc}
\le P_L(n)\,2^{(h+m_0)/2},\qquad
G_{\rm loc}\le P_L(n)\,2^{h+m_0},
```

```math
D_{\rm loc}=O_L(\ell S(\tau)),\qquad
S(\tau)=O(\log(\tau+2))
```

for the preferred result. A larger but $`o(\sqrt\tau)`$ stage count
also improves R86. Add coarse search, terminal, boundary routing,
sign gauges, and prescribed phases before claiming a general result.
Restore singular inputs through one consistent global rounding, and
use the original-accuracy narrow-width and finite-dimension fallbacks.

## 7. Five-hour execution and decisions

Use two independent constructors: one for the whole-comb and global
packet schedule, one for the reached-space reserve and its amplification
contract. A native-resource auditor follows both. A skeptical reviewer
tries to break the opening witnesses, nested repairs, and recurrences.
Reserve a fresh reviewer for the surviving complete argument.

| Elapsed minutes | Work and decision |
|---|---|
| 0–25 | Verify the baseline and full operator contract; reconstruct the spine witness and comb bound, including all endpoints. |
| 25–85 | Derive explicit packet reassociation and boundary repair; independently develop the reached-space interface. Produce equations, live registers, and a dependency graph. |
| 85–125 | Exchange counterexamples. Test nested trailing repairs first. Choose the primary mechanism and one distinct reserve; change the representation if its actual recurrence is unchanged. |
| 125–205 | Resolve the central global scheduling or transport lemma, maintaining native counts, precision, and return throughout. |
| 205–250 | Reconstruct the complete general bound or isolate the exact remaining lemma. Begin independent proof review and the necessary long checks as soon as the statement stabilizes. |
| 250–300 | Complete validation and integration, or preserve the reviewed branch with explicit pending gates and one precise resumption question. |

These are ceilings and decision windows, not time to fill. A completed,
reviewed and integrated target may finish earlier. The five-hour clock
includes tool waits and CI.

Stop an individual candidate when it hides a serial height-dependent
repair, assumes the unknown intermediate frame as a gate, omits an
address register, needs unavailable initialized history, or fails
literal work return. Reopen it only after changing that exact premise.
An unresolved but narrowing construction can use the discovery budget.
Ending a candidate or a session does not retire optimal general depth
as the continuing objective.

Keep one compact checkpoint with the exact identity, chosen
representation, register map, recurrence, error/return argument, strongest
counterexample, and next derivation. Use only a few meaningful finite
fixtures after symbolic identities survive. No broad simulations,
finite-size asymptotic fits, manuscript drafting, or external outreach.
Search primary literature only for a named missing lemma, then translate
all access and initialization assumptions before adopting it.

A theorem outcome requires the complete simultaneous bound and a fresh
reconstruction of its pivotal lemma. A mechanism outcome must state its
remaining full-circuit obligations. A failed route should retain its
scoped counterexample and reopening condition. Add no result number for
a conjecture. R86 and the existing scope of R79 remain authoritative
unless a proved correction or improvement is integrated.

Record the eventual execution in a separate result note and mark this
design executed. Use the [workspace validation requirements](../../WORKSPACE.md#restore-and-verify)
for substantive proof or circuit changes; design-only changes need the
relevant narrative, link, mathematical-display, and rendered checks.
Merge only after required CI passes on the reviewed tree. Report actual
elapsed time and distinguish internal review from external review.
