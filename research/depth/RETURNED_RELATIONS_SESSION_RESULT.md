# Native axes, scalar return and the missing relation

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Five-layer proof home](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Previous checkpoint](NATIVE_PARTIAL_ACTION_SESSION_RESULT.md)

This round began at verified main
`fe45d5ad04ccf10296142cf7feaffc8257adfd85` (PR #170).
It tested the proposed bounded operator relation for non-Pauli five-layer
comparisons, then reformulated after finding a native obstruction to
baseline-only cancellation. The destination remains optimal general
T-depth for the prescribed complete differential frame.

## 1. Established results and the exact general question

The count package remains settled under its stated hypotheses. R96 still
provides the strongest unconditional general complete-frame upper law:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)),
```

with one clean qubit, full error at most $`2^{-L}`$, all
$`n\ge1,L\ge6`$, and literal dirty thresholds
$`b\ge17(L+n+7)`$ for real frames and
$`b\ge17(L+n+8)`$ for complex magnitude frames. R79 remains selected
Result D. R97's outer-overlap hypotheses remain essential; unrestricted
one-clean four-layer rigidity is open.

The general question is to determine the optimal worst-case T-depth on the
same circuit and resource model, including the large-width range where
count divided by physical width does not match the additive depth above.
A robust dimension-growing unrestricted lower bound is missing. The
sufficient constructive route R98 still needs a faster exact multi-output
dirty query satisfying its simultaneous T-count, Clifford-count, depth and
workspace envelope. R91's logarithmic address-indicator cost remains.

Formerly open issues already resolved are not reopened: the one-clean
count frontier, the all-precision upper law, zero-clean read-only affine
three-layer and quadratic four-layer XOR bounds, and the five-layer
Pauli-comparison subcase. The latter remains degree at most two below
error $`\sqrt2/12`$. The separate Clifford-middle subcase remains degree
at most seven, or four with a visible outer target shear.

Every result here preserves its full operator contract: all logical
columns, literal phases, arbitrary output/dirty/reference inputs, actual
inverses, leakage and complete work return. No new complete-frame compiler
is claimed, so its treatment of signs, singular and zero coefficients and
simultaneous allocations is unchanged. Final manuscript drafting is separate.

## 2. New supported results

### A larger five-layer subcase

[Section 5.2 of the proof home](READ_ONLY_FIVE_LAYER_REDUCTION.md#52-an-outer-axis-normalizer-comparison-a-finite-degree-bound)
handles comparisons K that permute the individual active outer T axes,
with signs allowed. The middle B need not be Clifford and K need not be
Pauli. At the same full error radius $`\sqrt2/12`$, the Boolean degree is
at most 244; an invisible original target shear improves it to 41, or 31
when the peeled target is Pauli. These are conservative finite bounds,
not sharpness claims. The earlier stronger bounds still apply in their
own subcases.

Peeling the outer native band now leaves a Clifford cocycle. A finite set
of tableau features transforms linearly under fixed left/right Clifford
multiplication. The proof explicitly handles phase carries for a varying
Clifford tableau; it does not incorrectly divide polynomial coefficients
by two. Native degree bounds and a translation-space dimension argument
then bound the output degree independently of dirty width. An existing
fully returned native fixture satisfies this new hypothesis and neither
of the earlier two extra hypotheses.

[Section 5.3](READ_ONLY_FIVE_LAYER_REDUCTION.md#53-preserving-the-axis-span-freezes-a-third-level-cocycle)
extends the **cocycle reduction**, but not the degree theorem, to K that
preserves only the outer axis span. The peeled operators lie in a specific
subgroup of the third Clifford-hierarchy level. At error below
$`\sqrt2/24`$, the projective cocycle becomes exact. This uses closure of
that conjugated subgroup, not closure of the whole third level. Its
right factor need not be Clifford, which is the remaining feature-space
obstruction in this branch.

### A robust diagonal-source obstruction

The [returned diagonal echo theorem](RETURNED_DIAGONAL_ECHO_DEGREE.md#robust-return-at-a-width-independent-error-radius)
now forces Boolean degree at most three at full operator error below
$`2/275`$, for arbitrary unit-circle source phases. No phase-denominator,
work-width, independent-direction or exact-return assumption is needed.
A four-direction integer translation identity controls the error without
choosing phase logarithms or accumulating a width-dependent constant.

If an actual native source approximates the diagonal source in full
operator norm by delta, all four source/actual-inverse calls add at most
$`4\delta`$. Thus completed error epsilon still gives the cubic bound
when $`\varepsilon+4\delta\lt 2/275`$. The generic source ledger and all
simultaneous arbitrary work are explicit in the proof home. Relaxing
return slightly or increasing diagonal-source precision does not create
a high-degree query through this echo. This does not cover non-diagonal
sources or general native circuits.

### Two failed shortcuts now have explicit counterexamples

The [baseline-word obstruction](BASELINE_WORD_OBSTRUCTION.md) shows that
any fixed word cancelling every pair of independently chosen native
one-T baseline families also cancels a separate commuting Boolean target.
This remains true for universal projective error below $`1/\sqrt2`$.
An explicit native family has arbitrary odd symplectic order; a binary
incidence matrix then forces the net occurrence of every nonzero address
in the word to vanish. The witness costs 4m T gates, four layers and
$`4m-2+2ms`$ elementary Cliffords on m arbitrary work qubits, for an
s-bit address parity. Its actual inverse is charged. It is a source,
not a returned query. The theorem does not discard or constrain the
additional shared-source return cocycle of an actual five-layer circuit.

A second reformulation asked whether scalar return plus a bounded number
of abstract Pauli-representation factors could suffice. The
[four-factor construction](SCALAR_REPRESENTATION_PRODUCT.md) disproves that:
for any Boolean f with f(0)=0, honest representations A and B satisfy
$`A_xB_xA_xB_x=(-1)^{f(x)}I`$ with literal full return. Both use the same
canonical linear Pauli gauge and explicit fixed unitary conjugators.
There is no hidden address-dependent scalar.

The sources encode f in a bank of $`Q=2^a`$ arbitrary qubits. Writing
$`s=aQ/2`$ and $`m=Q\mathrm{wt}(f)/2`$, the expanded returned word
has the conservative native ledger

```math
T=56(s+m),\qquad G=112s+80m+4a+16,
\qquad D_T\le40(s+m),\qquad b=a+Q,
```

with zero clean work, zero error and every source/inverse charged.
The fair R91 scalar baseline at the same asymptotic width has
$`T=O(\sqrt Q)`$, $`G=O(Q)`$ and logarithmic address depth. The new
construction supplies no resource saving and no constant-depth native
counterexample. It shows exactly how abstract factor counting can move
all the difficulty into an unpriced source.

## 3. Concrete resumption question

The most useful next question is the **span-preserving but axis-mixing
five-layer branch**. Its projective cocycle is now exact at a fixed error
radius. Can one find a bounded-degree feature space for the peeled
operators that remains closed under their actual native conjugation and
right multiplication, and still detects the target Boolean function?

More concretely, write $`C_x=D^\dagger K_xD`$. When K preserves the
outer axis span, C is a commuting T-phase polynomial times a Clifford
normalizing that span. The native cocycle is

```math
C_{x+y}\doteq A'_y C_x A_y'^\dagger C_y.
```

Here A' is still the specified one-T Pauli-image Clifford family. Unlike
the individual-axis case, C_y need not be Clifford. The missing step is
control of this right action and its phase correction, not another
assumption of intermediate Cliffordness. A successful invariant must
retain these native factors; the two counterexamples above show why
baseline-only words and arbitrary representation counts cannot replace
them. A counterexample inside this stricter native branch would also be
useful, provided its actual circuit and all resources are exhibited.

This question extends the proved branch while keeping the exact obstacle
visible. It would not by itself prove a uniform degree-versus-depth law:
general axis mixing, growing-depth robustness and unrestricted address
mixing remain further obligations. Constructively, a non-diagonal query
mechanism must still beat R91's fairly simplified simultaneous envelope
before conditional R98 improves general frame depth.

The unsuccessful candidates have been stopped and reformulated. They do
not retire the continuing optimal-depth objective.

## 4. Evidence and integration

Each theorem has one authoritative proof home linked above. Bounded checks
verify native sources, actual inverses, all small operator columns, literal
phase carries, symplectic orders, translation identities, source ledgers
and the strict boundaries of the new subcases. Independent audits checked
the analytic arguments; finite enumerations are not used to extrapolate
in dimension. Validation and presentation results are recorded with the
reviewed merge. No manuscript theorem or production compiler is replaced.
