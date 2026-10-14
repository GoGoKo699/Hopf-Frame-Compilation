<a id="bounded-audit-of-the-selected-compiler-claims"></a>
<a id="reconstruction-of-the-one-clean-count-theorem"></a>

# Reconstruction of the count and depth theorems

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [Proof roadmap](../REVIEW.md) · [Verification](VERIFICATION.md) · [Source map](SOURCE_MAP.md)

**Decision, 14 October 2026:** the reconstructed proof supports Result C
with one clean qubit, its literal real/complex dirty reservations, and
the full initialized-isometry contract. No unresolved proof objection
remains from this bounded internal review. The reconstruction starts
from `753bb86472340ff4f92fdfe8edf263bf9e48dbf3` and includes the
shorter precision splice and small-dimension corollary below.

The authoritative statement is
[growing-source Section 11](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count),
with its banked proof in Section 12. The
[3 October audit](reference/CORE_CLAIM_AUDIT_2026_10_03.md)
preserves the earlier A–D decisions, imported-interface checks, correction,
baseline, and finite evidence. Its older Result C is a retained supporting
construction. A, B, and D retain their distinct resource statements.

A separate reconstruction below checks the
[one-clean depth theorem](../research/depth/ONE_CLEAN_DEPTH.md), R72,
which gives Result D optimal count at every eligible precision and width.
The compressed-column construction remains a stronger linear-precision
special case; fixed accuracy retains sublinear depth at sufficiently
large constant square-root dirty width.

## Claim-to-proof decisions

Write $`N=2^n`$, $`n\ge1`$, and integer $`L\ge6`$. Put
$`B=L+n+7`$ for real frames and $`B=L+n+8`$ for complex
magnitude frames. The claim fixes $`a=1`$ and error $`2^{-L}`$.

| Resource point | Reconstructed conclusion | Matching qualification |
|---|---|---|
| $`b\ge B`$ | $`T=O(N+L)`$, $`G=O(NL)`$ | At exactly $`b=B`$ and $`L\ge\max\{6,n\}`$, $`T^\star=\Theta(N+L)`$ |
| $`b\ge2B`$ | $`T^\star=\Theta(\sqrt{NL}+L+NL/b)`$, with constructive $`G=O(NL)`$ | Every $`n\ge1`$ and integer $`L\ge6`$, at the fixed one-clean allocation |

This does not declare every width between B and $`2B`$ optimal.
An arbitrary larger clean budget changes the lower-bound width to
$`n+a+b`$. B therefore retains its sufficient-clean regime,
including zero dirty work; D retains its simultaneous count/depth circuit.

### Fragile dependencies retained by the audit

| Obligation | Retained proof and reconstruction |
|---|---|
| Exact native programming | [Dirty interpreter](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter), [whole-word banks](OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost): full-input XOR action, disjoint live core/selectors/banks, literal Pauli order, and actual inverse |
| Conditioned native symbols | [Conditional-suffix Section 1](CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally): returned predicate and phase echoes price controlled symbols; controlled T is not a free primitive |
| Scalar and phase words | [One-clean Sections 2–5](ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks): native weight-two loaders, scalar anticommutator, five-call error, and exact inactive identity |
| Mixed source and full mode arena | [Heterogeneous sources](../research/endpoint/HETEROGENEOUS_PRECISION_SOURCE.md), [precision gauges](../research/endpoint/PRECISION_GAUGE_COMPILER.md), and [growing Sections 2–5](../research/endpoint/GROWING_SOURCE_COMPILER.md#2-a-polynomial-schedule-and-the-mixed-local-gauge): effective clipping, finite-shift estimate, normalized refinement, sealed rejection, same-core replay, and compensated overlap |
| Independent block/source interfaces | [Disjoint compiler](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#2-exact-routing-and-coarse-words): dirty routing, exact logical coarse words, local charts, and constant fine slots; these do not invoke the growing theorem |
| Borrowed signal and recovery | [Conjugate Sections 1–2](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md#1-the-actual-finite-phase-word): full-operator phase error, every native interlayer lifted, genuine initialized anchor retained, arbitrary signal excluded from the reflection |
| Uniform count and width | [Layerwise banked baseline](ONE_CLEAN_COMPILER.md#layerwise-banked-baseline), [banked interfaces](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces), and [growing Section 12](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier): rounded bank allocation and every fine, coarse, routing, predicate, inverse, and Clifford charge |
| Phase dressing and lower bounds | [One-clean Section 8](ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames), [complex diagonal subfamily](FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame), and [real lower bounds](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage) |

## Reconstruction findings

The phase source's Hermitian off-diagonal operator satisfies
$`B^2=(1-|z|^2)I`$. Its five-call error square is scalar on the
entire signal/core space. Thus the flag-zero estimate genuinely supplies
the full-operator bound on both arbitrary signal sectors. All S/T
interlayers and actual inverses are lifted with literal phases retained.

The borrowed signal remains outside every source core. Recovery retains
the original initialized suffix and its marker h, but excludes that
signal. Only the initialized final target must be real; the complementary
action may be complex. The retained two-sector counterexample explains
why the same argument cannot remove genuine recovery initialization.

Certified attenuation, state-chart, and replay data retain amplitude or
source-column error bounds. Only native fine-slot replacements require
full-operator bounds. Comparing independently initialized source columns
controls the overlap without assuming inverse-source accuracy on arbitrary
outputs. The recovery estimate includes the full rejected column.

Choose the fixed anchor height first, then one universal dimension
threshold. In the high regime $`L>N/n^2`$, rounded precision
$`K\in[L,2L)`$ dominates every required polynomial in n. Dyadic
integrality, multiplier bounds, finite-source error, self-borrowing, and
fixed helper slack therefore hold uniformly in L and the target data.
Coarse dirty work is $`o(N/n^2)`$; the fine core uses
$`L+O(1)`$ wires. Since $`\sqrt{NL}>N/n`$, the charged coarse
T-cost $`N^{3/4}\mathrm{poly}(n)`$ is absorbed. Coarse Clifford
cost $`N\mathrm{poly}(n)`$ is $`O(NL)`$. Signed-normalization
permutations use this same dirty routing ledger.

Below the cutoff, the layerwise banked bound is
$`O(\sqrt{NL}+nL+NL/b)`$, and $`nL\le\sqrt{NL}`$.
Bounded n uses that baseline for arbitrary L with one universal constant.
In particular, $`n=1,2`$ gives $`O(L)`$ at the original
literal width; the canonical root-angle witness gives $`\Omega(L)`$.
This removes iterated grouping from C's proof chain. The grouped theorem
remains a standalone result, and its native conditioning lemma remains
a necessary primitive.

The complex lower bound uses arbitrary leaf phases with all magnitude
angles zero. This diagonal subfamily lies inside the canonical complex
domain. It replaces the overly broad claim that zero leaf phases embed
the entire real-angle domain. The specific root-angle $`\pi/16`$
witness belongs to both families.

<a id="fixed-accuracy-large-width-depth-audit"></a>

<a id="precision-uniform-large-width-depth-audit"></a>

## All-precision count and depth audit

The [one-clean theorem](../research/depth/ONE_CLEAN_DEPTH.md), R72,
gives complete real frames a circuit with one initialized flag, error
at most $`2^{-L}`$, and absolute constants for every $`n\ge1`$,
integer $`L\ge6`$, and $`b\ge17(L+n+7)`$. Define

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
```

Then

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

The bounds hold simultaneously in the full initialized-isometry contract.
Count is optimal throughout; count and depth match through
$`b\le\sqrt{NL/F(n,L)}`$ when the interval is nonempty.
For $`L>N/n^2`$, the stronger depth bound is $`O(NL/b^2+L)`$.
Optimal depth beyond the stated matching range remains open.

R72 gives the same bounds for $`D_\phi W_{\mathbb R}`$ at
$`b\ge17(L+n+8)`$. Its literal-diagonal depth is proved on the
same optimal-count circuit. Both factors use precision $`L+1`$ and
reuse the same clean flag and dirty pool. The full-isometry hybrid propagates
the real factor's earlier leakage unitarily; no intermediate reset or exact
core return is assumed. The diagonal subclass with zero magnitude angles
supplies the complex family's count and width-dependent depth lower bounds.

| Reconstructed obligation | Finding and proof boundary |
|---|---|
| One-clean native baseline | [R72](../research/depth/ONE_CLEAN_DEPTH.md) schedules the existing one-clean source word with complete rectangular queries and returned predicate helpers. Its uncapped depth is $`O(NL/b^2+n(L+n))`$; the extra term is charged before the width split, not discarded. |
| Width and precision splices | R72 splits at $`b_c=\sqrt{NL/[n(L+n)]}`$: below it the width term absorbs the uncapped remainder; above it the protected native schedules fit after a fixed dimension cutoff. The proof checks optimal count, elementary Clifford cost, full error, and bounded dimensions on each selected complete circuit. |
| Terminal flag and recovery | The early shared-gradient construction uses one external activity flag. Terminal padding rejection moves into the protected logical bank, is included in the recovery reflection, and is returned before erasing the original predicate. The inactive full-space identity and source-return comparison are retained. |
| Literal diagonal and complex composition | R72 charges twenty completed mask queries and forty-five native source appearances for the one-clean phase word. Nonuniform indicators and the routed fallback respect the literal allocation. The source has precision $`q=\ell+5`$ and core width $`\ell+6`$; ordered masks, actual inverses, and literal phases retain the full operator. Half-error composition on the sole clean flag preserves the global F bound and high-precision clause. |
| High-precision branch | The unbounded-dimension [high-precision schedule](../research/depth/HIGH_PRECISION_DEPTH.md) already leaves its second clean flag unused. R72 supplies a one-clean bounded-dimension fallback and charges the completed-mask replacements and native remainder, giving $`O(NL/b^2+L)`$ for $`L>N/n^2`$. |
| Retained two-clean proofs | [R70](../research/depth/ALL_PRECISION_DEPTH.md) and [R71](../research/depth/COMPLEX_FRAME_DEPTH.md) retain their original statements. The [uniform baseline](UNIFORM_PRECISION_DEPTH.md), [whole-frame stream](../research/depth/WHOLE_FRAME_BINARY_DEPTH.md), and [high-precision schedule](../research/depth/HIGH_PRECISION_DEPTH.md) remain separately scoped supporting proofs; their two-clean assumptions are not transferred without the R72 arguments above. |
| Stronger linear-precision result | The retained two-clean [compressed-column proof](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) retains $`\log_2(n+2)`$ for $`6\le L\le cn`$. Its [descendant source](../research/depth/COHERENT_DESCENDANT_LOADER.md), [early shared gradient](../research/depth/SHARED_GRADIENT_RECOVERY.md), and separate [binary terminal stream](../research/depth/SHARED_BINARY_TERMINAL.md) keep their own exact inactive guards, complete-column recovery, failure history, and source-return comparisons. |
| Imported state source | GKW v3, Lemma 3.2 and Section 3.3, supplies the coarse Boolean-state premise used by the compressed-column construction. Its ancillary allocation is not imported as a two-clean implementation; local coherent source and native schedules remain necessary. |
| Scoped normalization barrier | The [entry-indexed bound](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md) concerns normalization under its explicit row/column factorization contract. It supplies no unrestricted native-depth lower bound. |

At sufficient $`b=\Theta(\sqrt{NL})`$, the linear-precision
special case gives sublinear depth for $`L=o(n/\log(n+2))`$.
R72 independently gives the one-clean fixed-accuracy all-n consequence
$`O_\eta(\sqrt n\log(n+2))`$ at sufficient square-root width.

Every selected branch includes source leakage, dirty/reference inputs,
literal phases, and actual inverse words in its full comparison.
Occupied programs are not silently identified with initialized programs.
Original activity predicates are erased only after the corresponding
protected work returns; unitary error composition never assumes a reset.
Absolute cutoffs enforce resource bounds without making practical
crossover claims. Complete arguments remain in the linked proof homes.

## Shared operator and operational contract

The separate [R73 obstruction](../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md)
uses that same initialized-isometry contract for a constant lower bound.
Its clean-projector sum contains $`2^a`$ full-space transfer entries, but
the witness error remains $`2\eta`$. A bounded-normalization three-layer
reduction and the cited polynomial-limit theorem yield a positive error
gap for each fixed clean budget. The elementary initialized-sector
argument in [R74, Section 5.1](../research/depth/FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector)
now makes that threshold explicit without the limit premise. At four
layers it gives minimum middle crossing rank $`\Omega_a(\sqrt L)`$.
The separate [operator-rigidity theorem R75](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md)
requires at least five layers for the zero-clean witness at integer
$`L\ge3`$, with arbitrary dirty width. Its one-clean consequence is
exact terminal clean-sector preservation only. These results do not
close the growing-depth gap.

The initialized-isometry norm covers every logical column, arbitrary
dirty/reference inputs, literal phases, and clean leakage. Actual inverses
inherit the same estimate. Sequential reuse carries previous errors
unitarily; approximate work is never assumed to reset between factors.
Exact bank/selector return remains distinct from approximate operator-core
return. QBP retains separate observable, interference-work, sampling, and
classical-output charges.

<a id="imported-premises-checked"></a>
<a id="correction-evidence-and-stopping-decision"></a>

## Evidence and stopping decision

The [source map](SOURCE_MAP.md) identifies the imported synthesis and
lower-bound premises; the [dated audit](reference/CORE_CLAIM_AUDIT_2026_10_03.md#imported-premises-checked)
retains their earlier interface checks. The count reconstruction imports
no additional external premise. The depth extension uses the cited GKW
state-preparation result with its local coherent and native realization;
no literature-priority claim follows from this audit. The [verification map](VERIFICATION.md)
distinguishes native words, exact identities, supplied-matrix diagnostics,
and analytic resource sums. The existing scientific tests and exact
receipt hashes are retained.

This is internal proof reconstruction, not external technical review or
formal proof certification. The selected count package is ready for
manuscript development when requested. Result D now has all-precision
optimal count and simultaneous depth bounds. Optimal depth beyond the
stated matching intervals remains a separate question.
