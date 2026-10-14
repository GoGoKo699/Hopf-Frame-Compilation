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
[fixed-accuracy large-width depth construction](../research/depth/SUBLINEAR_T_DEPTH.md).
It preserves C and the uniform D theorem while establishing a sublinear
depth upper bound at sufficiently large constant square-root dirty width,
with a width-sensitive extension at every eligible dirty allocation.

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

## Fixed-accuracy large-width depth audit

For every fixed $`0\lt\eta\le1/64`$, the
[complete construction](../research/depth/SUBLINEAR_T_DEPTH.md) gives
all prescribed real frames on $`n\ge1`$ qubits a circuit with two
initialized flags and $`b=\lceil C_\eta\sqrt N\rceil`$ dirty qubits:

```math
T=O_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(\sqrt n\,\log(n+2))=o(n).
```

These are simultaneous bounds with the full initialized-isometry error
at most eta. They settle the stated fixed-accuracy sublinear-depth
existence question; they do not determine optimal depth or replace D's
precision-uniform bounds. The fixed-accuracy width extension retains the
literal lower threshold $`b\ge17(L+n+7)`$ and gives
$`D_T=O_\eta(N/b^2+\min\{n,\sqrt n\log_2(n+2)\})`$,
$`T=O_\eta(\sqrt N+N/b)`$ and $`G=O_\eta(N)`$.
Its matching interval ends at
$`\sqrt N/[n^{1/4}\sqrt{\log_2(n+2)}]`$ when nonempty.

| Reconstructed obligation | Finding and proof boundary |
|---|---|
| Imported state source | [GKW v3, Lemma 3.2 and Section 3.3](https://arxiv.org/html/2411.04790v3#S3) supplies coarse Boolean-phase state approximation and LCU/amplification structure. Its ancillary allocation is not treated as a two-clean implementation. |
| Effective coherent preparation | The [descendant loader](../research/depth/COHERENT_DESCENDANT_LOADER.md) uses a local strict residual search, a common geometric coefficient state, and source-label-dependent attenuation. It terminates at zero residuals without equality decisions. Orthogonal tag/prefix sectors prevent a dimension factor in source-isometry error. |
| Native guards and live work | [Shallow dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md) expose private T-layer supports despite shared read-only addresses. The [loader](../research/depth/COHERENT_DESCENDANT_LOADER.md) prepares binary phase factors in parallel and uses exact signed-addition echoes. Its actual preparation inverse gives exact inactive identity on occupied source/program registers; active source leakage is charged in the full comparison. Fixed coefficient echoes have disjoint T supports. |
| Whole-operator recovery | The tag source and uniform row isometry encode the complete early frame with normalization $`\sqrt K`$, $`K=O(n)`$. A mathematical near-isometry completion controls actual forward and inverse source words in full operator norm, including leaked source scratch. Guarded amplitude amplification retains its literal active-sector sign. |
| Same-circuit resource budgets | A terminal block of $`M=\lceil\log_2^2(n+2)\rceil`$ supplies the simultaneous $`O_\eta(\log(n+2))`$ conditional clean reservation. Early tables have $`2^{n-M}\mathrm{poly}_\eta(n)`$ rows. The [nonuniform dirty indicator](../research/depth/NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence) keeps the final M layers within the count and Clifford budgets; they contribute only polylogarithmic depth at square-root width, with the additional $`N/b^2`$ term at general widths. |
| Literal-width splice | Below $`b=\sqrt{N/n}`$, the retained compiler absorbs n into $`N/b^2`$. Above it, the early width is $`o_\eta(\sqrt{N/n})`$ and the tail keeps its $`N/b`$, $`N/b^2`$ terms. The old compiler covers finite dimensions; the new indicator is not assumed to fit at every literal base allocation. |
| Scoped normalization barrier | The [entry-indexed bound](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md) uses an explicit factorization contract and the balanced Haar frame. Its exact and approximate bounds concern normalization only; they supply no unrestricted native-depth lower bound. |

An occupied active angle program represents its XOR-offset word, so the
audit does not identify initialized-program accuracy with a specified
full-space row rotation. The recovery comparison instead includes the
actual normalized source outputs and all scratch leakage. The original
terminal-zero predicate is reversed only after recovery; unitary
telescoping covers its leakage and the following tail without a reset.
Finite exceptional dimensions use the retained compiler with enlarged
fixed-accuracy constants.

## Shared operator and operational contract

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
manuscript development when requested. The large-width construction
establishes sublinear depth and extends the matching width range at fixed
accuracy. Optimal depth beyond that range and precision-uniform refinements
remain separate questions.
