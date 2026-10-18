# Native defects close the five-layer read-only degree gap

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound) ·
[Previous checkpoint](SPAN_PRESERVING_SESSION_RESULT.md)

[Subsequent checkpoint](SIX_LAYER_TRANSPORT_SESSION_RESULT.md): the six-layer
Pauli-comparison case is now bounded; the general six-layer question remains open.

This round began on 18 October 2026 at verified main
`6e20eb47bc7d857f88f4a4da8577b3a23457419e` (PR #172).
It addressed target extraction outside the outer-axis span and tested a
whole-word constructive alternative. Optimal general T-depth remains the
destination; final manuscript drafting is separate.

## 1. Established frontier and exact general question

The count package remains settled under its stated hypotheses. R96 remains
the strongest unconditional upper law for both prescribed complete frame
families, for every $`n\ge1,L\ge6`$, at full error at most $`2^{-L}`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad N=2^n.
```

One clean qubit suffices, with real dirty width $`b\ge17(L+n+7)`$
or complex magnitude width $`b\ge17(L+n+8)`$. These bounds hold on
one circuit. Every logical column, literal sign and phase, singular or
zero coefficient, arbitrary dirty/reference input, actual inverse,
leakage and complete workspace return remains part of the contract.
R79 remains selected Result D. R97 retains its outer-overlap hypotheses.
Unrestricted one-clean four-layer rigidity remains open.

The general question is the optimal worst-case T-depth at the same precision,
allocation and simultaneous T-count and elementary Clifford budgets,
especially where count divided by width does not match R96. A growing
unrestricted native-query lower bound would transfer through the two-call
corner reduction. A faster exact read-only multi-output query satisfying
R98's count/Clifford/width envelope would improve the fixed-accuracy frame
upper bound. These are different hypotheses; neither transfer has its
missing premise supplied here.

Previously open questions already resolved include one-clean optimal count,
R96's all-precision upper law, zero-clean read-only affine three-layer and
quadratic four-layer XOR rigidity, and the five-layer conditional subcases.
This round resolves the remaining **width-independent five-layer Boolean
degree bound**, without a comparison-span or Clifford-middle assumption.
The sharp degree and optimal general depth remain open.

## 2. The supported new theorem and why the mechanism changes

[Section 5.6](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
proves that a zero-clean read-only five-layer XOR query at full operator
error $`\eta\lt\sqrt2/24`$ has Boolean degree at most **14**, independently
of finite dirty width. An invisible outer target shear improves this to
**11**. Completed Clifford interlayers preserve the original address
pointwise, and T bands act on work only. T-count and elementary Clifford
count are unrestricted in this necessary condition. The proof adds no
physical source, control, inverse or workspace. It does not give a compiler.

The preceding intrinsic comparison features did not directly detect the
output relative to a baseline outside their group. Instead of extending
that group, this round uses the ordered defect of the actual comparisons.
Its shared native prefix pulls back to a Pauli with a bilinear address label.
The first of the two T-band passages therefore gives a known Clifford.
Robust one-band matching applies to the last passage at error $`6\eta`$,
freezing the defect's projective action without any width-dependent loss.

For each fixed translation, the defect has symplectic degree at most five.
Its literal sign on a degree-two moving Pauli probe has degree at most ten.
The existing visible-target linear extractor and a translation-dimension
argument give 14; the invisible-target signed derivative identity gives 11.
Unlike the earlier failed partial-normalizer argument, the needed Clifford
intermediate is proved from a Pauli input. Projective scalars are discarded
only after the full-operator comparison, for conjugation actions; literal
Hermitian-Pauli signs are retained.

The span-preserving 112/21 bounds are superseded at their error radii; their
augmentation and shared-prefix structural results remain valid. The sharper
Pauli-comparison and Clifford-middle subcases, and the larger-radius
individual-axis bound, retain their respective hypotheses. Degree at least
15 now excludes five or fewer layers in this read-only zero-clean model at
the stated error, but this is still a constant-depth obstruction.

## 3. Failed alternatives and their useful diagnosis

The [returned affine example](BASELINE_WORD_OBSTRUCTION.md#5-a-returned-affine-query-can-have-an-odd-order-outer-baseline)
has an actual comparison group of order two while its formal outer baseline
has odd symplectic order $`2m+1`$. It is a literal fully returned query,
not just independent algebraic data. Its displayed native ledger is
$`4m`$ T gates, four nonempty T layers, $`6m-4+2ms+s`$ elementary
Cliffords, zero clean and $`m-1`$ dirty helpers, for an address parity
using s bits. Fair simplification gives s CNOTs and zero T depth. This
stops an unconditional ambient-2-group inference. It supplies neither a
high-degree query nor a resource gain, and does not contradict the theorem.

For the [constructive fusion attempt](SCALAR_REPRESENTATION_PRODUCT.md#6-whole-word-fusion-removes-carrier-motion-but-leaves-the-table-bit),
the four-factor scalar reduces to the commutator square of one transported
Pauli and a fixed probe. All irrelevant carrier motion cancels; the sole
surviving coefficient is exactly the table bit f(x). In a rectangular table
factorization, loading that coefficient still requires the retained dirty
indicators. R91 supplies $`O(\log_2(a+2))`$ additive depth with its
simultaneous square-root T-count and linear Clifford count. No native
identity replacing that step was found. The separately paid source-cost
barrier does not exclude a different fused or non-Clifford implementation.

## 4. Precise resumption question

The most useful next question is a **native degree-versus-depth recurrence**,
starting with six layers. This is preferable to optimizing the constant 14:
a recurrence could turn growing Boolean degree into growing read-only depth.
The [six-layer boundary calculation](READ_ONLY_FIVE_LAYER_REDUCTION.md#71-the-six-layer-native-boundary)
retains three shared T bands and a Clifford-valued inner family formed from
two differently oriented one-band Pauli images and a Pauli factor. Its
ordered defect need not be Pauli. A one-qubit native inner example has odd
symplectic order three; it is not a returned six-layer query counterexample.

Can full return constrain this particular three-band conjugation of the
native defect enough to bound both its tableau degree and its moving-probe
sign degree independently of dirty width? The next argument must prove any
intermediate Clifford property, or replace matching by a different native
invariant. Track word lengths and full-operator error radii as functions of
depth from the outset. A shrinking error threshold without suitable control
would not settle fixed-accuracy depth growth. Even a successful read-only
recurrence leaves initialized work and arbitrary address mixing to be
addressed before the unrestricted complete-frame lower-bound transfer.

The constructive alternative remains an exact reader with sublogarithmic
additive address depth inside R98's simultaneous resource envelope. The
five-layer theorem does not exclude six-or-more-layer readers. The next
candidate must replace the paid table-bit loading operation itself and
include all source, query, control, recovery and actual inverse costs.

## 5. Evidence and limits

The new [native-defect checks](../../tests/test_five_layer_native_defect.py)
exercise the bilinear Pauli pullback, actual inverses, matching and its
failed arbitrary-intermediate alternative, signed axis corrections,
Boolean moving-probe formulas, the error budget and translation dimensions.
The [returned-affine checks](../../tests/test_returned_affine_baseline.py)
cover every small address/work column, literal comparison extraction,
odd orders and the complete ledger. The
[six-layer boundary checks](../../tests/test_six_layer_defect_boundary.py)
test the next native factorization and its non-Pauli defect.
These bounded checks support the analytic arguments; they do not prove an
arbitrary-width statement by enumeration. Repository validation and rendered
presentation receipts are recorded with the integrating pull request.
