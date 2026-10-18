# Native partial conjugation: cube relations and returned phase mechanisms

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Five-layer proof home](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Previous checkpoint](FIVE_LAYER_QUERY_SESSION_RESULT.md)

This round began at verified main
`d484d25d2a6d24f7f9594cb5a92a5594aba1f78f` (PR #169).
It addressed the missing bounded-degree action of a native two-T-dressed
Pauli on the Clifford comparisons it conjugates back into Cliffords.
The destination remains optimal general T-depth for the complete prescribed
frame. No unrestricted depth conclusion is inferred from the scoped results.

## 1. Established frontier and the unresolved question

The count package is settled under its hypotheses. R96 remains the strongest
unconditional complete-frame upper law, with one clean qubit, optimal
T-count, linear-in-NL elementary Clifford count and

```math
D_T=O(NL/b^2+L+\log_2(n+2))
```

at all precisions and the retained literal real/complex width thresholds.
R79 remains selected Result D; R97 keeps its outer-overlap hypotheses;
unrestricted one-clean four-layer rigidity remains open. R98 remains
conditional on a faster exact multi-output query satisfying its simultaneous
count, Clifford and workspace envelope.

For the read-only query interface, the zero-clean three-layer affine and
four-layer quadratic XOR bounds are established. The five-layer theorem
with a Clifford-valued middle still gives degree at most seven, or four
for a visible outer shear, at error below $`\sqrt2/12`$.
The **general five-layer case remains open**. A uniform degree-versus-depth
law, robust at growing depth and compatible with unrestricted address
mixing, has not been proved. Neither the earlier precision restriction in
R96 nor the settled count questions is reopened.

## 2. What this round establishes

**A robust quadratic subcase without a Clifford middle.**
[Section 5.1 of the five-layer home](READ_ONLY_FIVE_LAYER_REDUCTION.md#51-a-pauli-valued-comparison-robust-quadratic-bound)
proves degree at most two when the specified comparison family K is
Pauli-valued, at full operator error below $`\sqrt2/12`$. B may be
non-Clifford. The Pauli label of K is a quadratic baseline plus a fixed
nonzero label times the unknown Boolean function. An eight-factor cube
argument first gave a cubic bound at the smaller radius $`\sqrt2/48`$.
A native peeling step then freezes the projective cocycle at the stronger
radius, without rounding B to a Clifford. The quadratic baseline's
coefficients lie in the outer native-axis span; their conjugation is
controlled by the independent inner axes. A cubic target would require
contraction rank at least three where those axes permit rank at most one.
This excludes the cubic case. The extra Pauli-comparison hypothesis remains
explicit and is not inferred from return. No sharpness claim is made.

**Exact and robust matching of one native band.**
[The one-band theorem](ONE_BAND_AXIS_MATCHING.md) recovers independent T
axes from Pauli coefficient norms. If two one-T bands match through
Cliffords to error below $`\sqrt2/4`$, the intervening Clifford must map
the individual active axes bijectively, allowing signs and inactive wires.
Preserving just their span is insufficient. Conversely that matching
makes the band passage Clifford. This is a necessary native invariant,
not an upper construction.

**A cheaper complete scalar phase circuit.**
[The hybrid phase construction](READ_ONLY_PHASE_RELOCATION.md#3-parallel-phase-words-and-a-sharp-scalar-example)
retains four phases on an original dirty bit and relocates only four
others. It improves the previous four-layer scalar CCZ word from
105 T / 346 Cliffords / 28 dirty helpers to

```math
T=64,\qquad G=218,\qquad D_T=4,\qquad b=17,\qquad a=0.
```

All three CCZ arguments remain read-only address wires. Every T target
within a layer is distinct; every bank computation and its actual inverse
is charged. Error and leakage are zero, and all dirty work returns on every
input, including external references. The retained three-layer obstruction
makes the depth sharp within this interface. The cheaper existing AND
fixture remains the appropriate XOR baseline. For M R91 bilinear phase
terms, the alternative costs 64M T, 17M simultaneous dirty bits and
$`196M+22\sum\ell_{ij}`$ Cliffords in four layers. Both old and new phase
components already have constant depth; the logarithmic address-indicator
cost and complete-frame law do not improve.

**The diagonal-echo limitation includes arbitrary phases.**
[The authoritative echo proof](RETURNED_DIAGONAL_ECHO_DEGREE.md) now allows
any diagonal source with phases in the full unit circle, including arbitrary
precision and irrational angles. Exact scalar return itself forces a
Boolean phase, still of degree at most three. Signed derivative exchanges
and triangles of subset exchanges force two-torsion; the Boolean-source
proof becomes a special case. The fully charged cubic native witness
remains sharp. There is no approximate-return or non-diagonal extension.

## 3. Failed inferences and their diagnoses

The one-band theorem cannot simply be used twice: the intermediate
$`D^\dagger K D`$ need not be Clifford even when the complete
$`D A D^\dagger K D A^\dagger D^\dagger`$ is Clifford.
An exact one-qubit cancellation demonstrates this, with both its displayed
native word and simplified Clifford baseline charged in the proof home.
Thus the apparent reduction would move the original missing assumption to
an intermediate operator.

Discarding the native constraint is also unsuccessful. The
[full-tableau counterexample](ONE_BAND_AXIS_MATCHING.md#5-an-abstract-normalizer-can-have-arbitrarily-high-tableau-degree)
uses a controlled translation of a bank of arbitrary work wires. It is an
involutory normalizer of an abelian Clifford subgroup; a spectator X makes
it an abstract Pauli image. Its partial tableau action has unbounded degree.
Every source and inverse can be implemented natively with the displayed
growing cost. This is neither a two-T construction nor a lower bound against
all resyntheses. It proves that abstract subgroup normalization alone cannot
supply the desired bounded-degree map.

Pairwise Choi overlaps or relative spectra also lose essential information:
for any Boolean g, the Pauli family $`Z^x\otimes X^{g(x)}`$ has orthogonal
Choi vectors and the same balanced relative spectra at distinct addresses.
The new cube proof instead transports multi-operator relations. A general
bounded-length product-of-representations reformulation was explored but
not closed; no degree theorem follows from it in this round.

## 4. Concrete resumption question

Can the eight-factor cube argument be replaced by a **uniformly bounded
length relation for the non-Pauli Clifford comparisons**, cancelling the
two native one-T baseline families while detecting a high Boolean derivative
of the target? The relation must survive conjugation by the actual native
B and retain the five-layer address cocycle. It cannot assume either
intermediate band passage is Clifford. A length bound independent of dirty
width would also give a controlled error budget, as the preliminary cube
argument demonstrates. The stronger final subcase uses additional native
axis information; it does not establish such a relation in general.

This is the most useful next lower-bound question: it extends the relation and native-axis
mechanisms at precisely the point where the general five-layer
proof remains incomplete. A bounded-degree native partial-tableau map is
still an alternative sufficient route. An explicit non-diagonal family
violating such a relation would guide constructive work, but all source,
query, control, inverse and recovery costs would then have to fit R91's
simultaneous envelope before R98 improves general frame depth.

One closed subcase is not a uniform depth law. Dimension-growing unrestricted
lower bounds and an asymptotically faster complete-frame construction remain
open. The unsuccessful candidates have been reformulated; the general
objective is unchanged. Final manuscript drafting remains separate.

## 5. Evidence

Each result has one authoritative proof home, linked above, and bounded
checks for its native identities, phases, inverse order, error argument or
resource ledger. Independent mathematical audits checked all four result
packages. Finite enumeration supports the analytic arguments and does not
extrapolate in width. Validation and presentation results are recorded with
the merge on the reviewed tree.

## Subsequent checkpoint

The [returned-relations round](RETURNED_RELATIONS_SESSION_RESULT.md), starting
at main `fe45d5a`, extends the native-axis subcase, proves robust diagonal-echo
rigidity, and tests the proposed operator-word and scalar-factor shortcuts.
The next-step discussion above is the historical handoff to that round.
