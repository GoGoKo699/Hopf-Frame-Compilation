# Span-preserving comparisons and the cost of returned sources

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Five-layer proof home](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Previous checkpoint](RETURNED_RELATIONS_SESSION_RESULT.md)

This round began at verified main
`b8ecbfbe37856f4d67a1cc479b1b5ba9317d4426` (PR #171), on 18 October
2026. It addressed the previous checkpoint's missing right-action invariant
and tested whether compressing its constructive counterexample could save
native resources. Optimal general T-depth remains the destination.

## 1. Established frontier and the exact remaining question

The count package is settled under its stated hypotheses. R96 remains the
strongest unconditional general upper law for both prescribed complete
frame families, with one clean qubit and full error at most $`2^{-L}`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It holds for every $`n\ge1,L\ge6`$, at real dirty width
$`b\ge17(L+n+7)`$ or complex magnitude width
$`b\ge17(L+n+8)`$. R79 remains selected Result D. R97 retains its
outer-overlap hypotheses; unrestricted one-clean four-layer rigidity and
dimension-growing unrestricted depth lower bounds remain open.

The exact general question is the optimal worst-case depth with the same
operator, precision, clean/dirty allocation and simultaneous count budgets,
especially beyond the range where count divided by width matches depth.
R98 would improve general fixed-accuracy depth from a faster exact
multi-output dirty query meeting its simultaneous resource envelope.
R91 still supplies logarithmic additive address depth; no faster general
query is constructed here.

Previously resolved questions remain resolved: one-clean optimal count,
R96's all-precision general upper law, zero-clean read-only affine
three-layer and quadratic four-layer XOR rigidity, and the five-layer
Pauli-comparison and Clifford-middle subcases. This round additionally
resolves the **span-preserving five-layer branch** posed in PR #171.

## 2. A finite invariant for axis mixing

[Section 5.4 of the proof home](READ_ONLY_FIVE_LAYER_REDUCTION.md#54-span-preserving-comparisons-a-finite-degree-bound)
proves that a zero-clean read-only five-layer XOR query has Boolean degree
at most 112 if its specified comparison K preserves the span of the active
outer native T axes and full operator error is below $`\sqrt2/24`$.
Individual axes may map to products; K need not be Pauli and the middle B
need not be Clifford. The bound is conservative. It does not assert that
full return implies the extra span hypothesis.

The earlier cocycle reduction was exact projectively at this radius but
left a non-Clifford right factor after peeling. The new proof changes the
features. It uses the group algebra of the **actual comparison group**, so
actual middle conjugation and right multiplication act linearly without
requiring a normalizer of the larger ambient Clifford group.

For nonaffine output, span preservation separates the native baseline and
target shears. Their associative increment algebra satisfies $`J^4=0`$.
The resulting ambient projective Clifford 2-group has trivial eighth
Jennings dimension subgroup. A locally proved subgroup-extension lemma
then gives bounded-degree native coordinates for the intrinsic features.
The nonidentity target survives the eighth augmentation truncation; a
translation-space dimension argument yields 112. The classical Jennings
formula and weighted basis are explicitly attributed at their point of use.

This handles arbitrary finite dirty width without a width-dependent error
loss. The proof adds no physical gate, source, inverse or workspace. Literal
phases are retained through the full operator approximation before taking
projective algebraic features. All dirty inputs, output columns and external
references remain in the norm. The earlier bounds at the larger radius
$`\sqrt2/12`$ keep their separate scopes: two for Pauli comparisons,
seven/four for a Clifford middle, and 244/41 for individual-axis permutations.

## 3. Shared native prefixes give an unrestricted comparison invariant

A further reformulation tested whether the odd-order independent baselines
could actually occur as comparisons of a returned five-layer circuit.
[Section 5.5](READ_ONLY_FIVE_LAYER_REDUCTION.md#55-shared-prefix-comparisons-a-robust-quadratic-group-map)
answers this without a span-preservation hypothesis. At error below
$`\sqrt2/40`$, the actual projective comparison group has nilpotency class
at most two and exponent dividing four. The address map has a central
bilinear multiplication defect. Consequently the coordinates of
$`[K_x]`$ in $`\mathbb F_2[G]/I_G^3`$ have Boolean degree at most four,
and distinct elements of G remain distinct in this quotient.

This uses the shared prefixes explicitly: B times J is a fixed conjugate
of a Pauli family times a single-T-band Pauli-image family. Their cross
commutators are native-axis Paulis central in that inner group. Short word
relations then transfer to K by full operator error and the projective
Clifford gap. No long word accumulates a width-dependent error.

Thus arbitrary odd-order **actual comparison groups** are excluded at the
stated radius. Odd order in the larger independent-baseline group remains
possible. The distinction matters: low-degree intrinsic features now exist
without the span hypothesis, but extracting the target still compares K
with a baseline that need not lie in G. The general degree theorem remains
open. No assumption that the original B normalizes every Clifford or every
Pauli is introduced.

Combining the quadratic group normal form with the span-preserving target
detector gives the sharper conditional bound **21** below $`\sqrt2/40`$:
the eighth intrinsic augmentation quotient has feature degree at most 14,
and the ambient inverse baseline adds at most seven. The degree-112 result
retains its larger error radius. Neither numerical constant is claimed sharp.

## 4. A constructive candidate fails for a measurable reason

The [scalar-source proof](SCALAR_REPRESENTATION_PRODUCT.md#5-a-source-cost-barrier-survives-pauli-carrier-compression)
now goes beyond counting the displayed gate expansion. Any Clifford
translation source on Pauli carriers with commutation matrix
$`F_{uv}=f(u+v)`$ obeys

```math
w\ge\tfrac12\mathrm{rank}_{\mathbb F_2}F,\qquad
T(U)\ge\mathrm{rank}_{\mathbb F_2}\bigl(F(P_x+I)\bigr),
```

where w is carrier width and P_x is translation of the table labels.
These bounds permit arbitrary carrier embeddings and additional arbitrary
helpers returned as identity. For degree d they force
$`w,T(U)\ge2^{d-1}`$. Full degree on an a-bit address therefore costs
at least $`2^{a-1}`$ T gates per source. The original uncompressed bank
has the stronger bound $`T(U)\ge2^a`$, even if the final scalar is trivial.

The mechanism is operator-level: a Clifford's Pauli expansion is flat on a
coset whose size is its symplectic displacement rank exponent, while each
physical T gate can at most double a conjugated Pauli's support. Full-source
approximation below $`1/(2\sqrt2)`$ preserves the integer count bound.
Eight separately implemented source/actual-inverse calls must therefore
pay at least four times the table size for full-degree tables. Their errors
add forward; small completed-word error does not imply source accuracy.

At comparable dirty width, the fairly simplified R91 scalar baseline has
square-root T-count, linear elementary Clifford count and logarithmic
address depth, with both exact returned queries charged. Carrier compression
cannot place the source-call mechanism inside that envelope. Its count-to-width
depth bound is explicit but gives no growing-depth lower bound at large width.
No independent elementary Clifford lower bound is claimed.

This stops **separately paid Clifford translation sources** as the proposed
improvement. It does not stop fusion or a different synthesis of the entire
returned scalar, non-Clifford sources, or general query constructions.
The full existing source ledger and simultaneous allocation remain in the
proof home; no new compiler or complete-frame upper bound is claimed.

## 5. Precise resumption question

The most useful lower-bound continuation is now **target extraction from
the actual comparison features when K moves the outer axis span**. Write

```math
K_x=K_x^0 Q^{g(x)},\qquad
v(x)=[K_x]\in\mathbb F_2[G]/I_G^3.
```

The new shared-prefix theorem gives faithful degree-four features v. Can
one recover g from these features and the specified outer native baseline
with Boolean polynomial degree bounded independently of dirty width?
A concrete sufficient form is $`g\in\mathcal P_r\,\mathrm{span}(v)`$
for an absolute r. That would give a finite general five-layer degree bound.
The question must retain the actual two native one-T baseline factors and
Q, rather than allow an arbitrary address-dependent reference operator.

The missing step is a **target detector**, not bounded nilpotence of the
actual group. The span theorem supplies a detector by putting baseline and
target into one bounded-filtered ambient group. Outside span preservation,
the [independent-baseline odd-order witnesses](BASELINE_WORD_OBSTRUCTION.md)
rule out assigning that same filtration to an arbitrary ambient group.
They do not refute the shared-prefix theorem. A different relative feature
or a compatibility identity between the inner and outer descriptions is
needed. Removing this last comparison hypothesis would settle the retained
zero-clean read-only five-layer branch, while unrestricted address mixing,
one clean qubit and dimension-growing depth would still require further work.

The constructive alternative is equally specific: eliminate or fuse the
full Clifford translation source and exhibit a native returned word query
that improves R91's additive address depth while retaining its T-count,
elementary Clifford count and simultaneous dirty-width bounds. The source
lower bound cannot be applied after such a change of mechanism.

The original full-frame contract remains unchanged: every logical column,
literal signs and phases, singular and zero coefficients, arbitrary dirty
inputs and references, actual inverses, leakage and complete workspace
return. Final manuscript drafting remains separate. A closed subcase or a
failed candidate does not retire the general-depth objective.

## 6. Evidence

The authoritative arguments are the two linked proof homes. The
[augmentation checks](../../tests/test_span_preserving_augmentation.py)
test native shear nilpotence, subgroup filtration behavior and the degree
threshold. The [shared-prefix checks](../../tests/test_shared_prefix_comparison.py)
exercise the central defect, short relation lengths and intrinsic features. The [source checks](../../tests/test_scalar_representation_product.py)
retain full-column scalar/source/inverse tests and add commutation rank,
Clifford support and approximation checks. Independent audits checked the
analytic arguments and their scope; finite instances do not prove the
arbitrary-width conclusions. Validation and integration receipts are recorded
in the associated pull request.
