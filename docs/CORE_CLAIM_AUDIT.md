# Bounded audit of the selected compiler claims

[Publication scope](../manuscript/PUBLICATION_SCOPE.md) · [Verification](VERIFICATION.md) · [Source map](SOURCE_MAP.md) · [Research decision](../WORKSPACE.md#research-decision-and-stopping-rules)

**Decision, 3 October 2026:** no unresolved claim-level blocker was found in
Results A–D and their retained compiler corollaries. Freeze the selected
scientific scope. One exact-workspace sentence is corrected below; the
theorem statements and resource frontiers are unchanged.

The reviewed baseline is `67a4cf25be30ec0440dd6e63a85181ca9c38121d`.
This bounded internal review read the proof arguments, checked the fragile
algebra and resource sums, and verified load-bearing source interfaces.
It is not formal proof certification, independent peer review, a new
priority survey, or a proof supplied by passing tests. The separate
state-based QBP theorem and later component refinements are not additional
principal results of this audit.

## Claim-to-proof decisions

Notation and full statements are fixed in the
[publication scope](../manuscript/PUBLICATION_SCOPE.md#four-principal-results).
Here $`N=2^n`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
$`h=1+\lceil\log_2(L+n+2)\rceil`$, $`q=n+a+b`$,
$`B_0=L+n+7`$, and $`\ell_*(n)=1+\log_2^*(n+2)`$.
Approximate results assume $`n\ge1`$ and $`0\lt\eta\le1/64`$.
Upper bounds apply to every admitted parameter tuple; matching lower bounds
are worst-case over the stated family.

| Result and decision | Proof homes | Checks and qualifications |
|---|---|---|
| **A: pass.** Exact size $`\Theta(N)`$, CNOT count $`\Theta(N)`$ for $`n\ge2`$, depth $`\Theta(n+N/(n+m))`$ | [Exact compiler](COMPILER_THEOREM.md), Sections 4–10 | Four-sector zero-work echo; decoder code-space preservation; copy erasure before flag reuse; coherent routing and actual unroute; maximal-cut splice for every $`m\ge0`$. Parameter count, CNOT slots with arbitrary ancillary width, and backward cones give the matching lower bounds. Literal complex phases compose as a complete multiplexor. CNOT count is zero at $`n=1`$. |
| **B: pass.** $`T^\star=\Theta(\sqrt{NL}+L+NL/q)`$, $`G=O(NL)`$ | [Fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md), Sections 2–10 | Exact dirty-bank cancellation; logarithmic-width finite precision source; complete sparse residuals; retained-source shifts; compressed failure history; full-output amplification; count and Clifford sums. Requires $`a\ge C(n+h)`$ for a sufficiently large absolute constant. Real and phase-dressed complex magnitude families are separately covered. |
| **C: pass.** One-clean real-frame $`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ | [One-clean compiler](ONE_CLEAN_COMPILER.md), Sections 2–6; [conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md), Sections 1–10 | Paired-source signs, anticommutator identity, literal five-call amplification, disjoint residual columns, separate scalar/atom flags, inactive identity, suffix relocation, iterated grouping, and finite-size fallback. The literal reservation is $`a=1,b\ge B_0`$; no intermediate reset or renewed clean-input assumption is used. |
| **D: pass.** Simultaneous uniform and low-precision T-count/T-depth bounds | [Uniform precision](UNIFORM_PRECISION_DEPTH.md), Sections 2–5; [blocked queries](BLOCKED_BILINEAR_LOOKUP.md); [unary groups](UNARY_PHASE_GRADIENT.md); [parallel lookup](PARALLEL_DIRTY_LOOKUP.md) | Rectangular power-of-two allocation, full-input bilinear cancellation, convergent weighted tail sums, uniform early cutoff, helper peaks, source preparation/return, and finite-size fallback all fit $`a=2,b\ge17B_0`$. Applies to complete real frames. The selected circuit meets its count and depth bounds simultaneously. |

### Fragile dependencies retained by the audit

For A, the decoder inverse remains valid because the intervening Givens
word preserves the complete one-hot code space. Routing preserves the
prefix needed to erase coherent copies before those wires become flags.
For the CNOT lower bound, at most twice the CNOT count ancillary wires
participate; untouched wires contribute only a removable common phase.
Fusing local gates leaves the stated finite-dimensional parameter bound.
The exact arbitrary-one-qubit model is essential here.

For B, the failure counter receives at most one increment per intermediate
factor and cannot wrap. Its final zero sector therefore selects the intended
product of accepted blocks. The independent normalization mode is excluded
from those tests. The cubic accepted-block estimate is upgraded locally to
a bound on the entire output, including rejected work. The clean reservation
leaves a constant fraction of total width available for the lookup bank.
The real final-layer subfamily contains arbitrary diagonals on one fewer
qubit with a stabilizer target; the single-qubit case is handled separately.
This permits the imported diagonal lower bound without confusing state and
operator synthesis.

For C, every completed inactive subroutine is exactly identity even when
the suffix is arbitrary. On the active sector the suffix supplies only its
declared conditional clean work. Group errors and the one-clean tail fit
the total error budget, and actual uncomputation carries existing leakage.
The banked bound is
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ at $`b\ge2B_0`$.
It matches the inherited lower bound when $`L\ell_*^2(n)\le N`$ or
$`b\le N/\ell_*(n)`$. Literal diagonals and complete U(2) multiplexors
retain their own reservations. Splitting real-frame and diagonal error
raises the complex magnitude thresholds to $`L+n+8`$ and $`2(L+n+8)`$.
The zero-clean corollary is only the separately proved layerwise real-frame
construction; it does not extend grouping or literal phase synthesis.

For D, the rectangular allocation bounds count and depth using the same
query circuit. The weighted tail sums converge uniformly in precision.
In the low-precision branch, the early cutoff has exponent
$`9\log_2(3)/16\lt1`$, so its costs and
helpers fit the claimed asymptotic allowances; the finitely many smaller
dimensions admit one absolute fallback constant. At general precision,

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+nL\right),\qquad G=O(NL).
```

For $`6\le L\le\log_2(n+2)/16`$, the count drops its $`nL`$ term
and the additive depth term becomes n. The matching intervals remain
$`17B_0\le b\le\sqrt{N/n}`$ at general precision and
$`17B_0\le b\le\sqrt{NL/n}`$ in that low-precision range, whenever
nonempty. Dividing the count lower bound by physical width supplies the
depth lower bound. Outside those intervals no unrestricted depth optimality
is inferred. The sufficient condition $`L\le N/n^2`$ extends count
optimality only. T-depth permits arbitrary Clifford interlayers and is
not total physical depth.

## Shared operator and operational contract

The [frame-safe proof](FRAME_SAFE_COMPILATION.md) and
[QBP approximation proof](QBP_APPROXIMATION.md) use

```math
\left\|VJ_a-J_a(W\otimes I_D)\right\|\le\eta.
```

This controls all logical inputs, arbitrary dirty work, external references,
literal phases, and clean-work leakage. The actual adjoint satisfies the
same initialized-isometry bound. Sequential reuse telescopes against ideal
initialized embeddings while unitarity carries previous leakage; it does
not assume that actual work becomes exactly zero after each approximate
factor. Exact dirty return is a separate stronger promise for B and exact
lookup/control components. One-clean operator-core return is approximate
within this norm.

Fixed-decoder frame necessity retains its quantifiers: exact prepared state,
regular nonzero coordinate amplitude, and every allowed observable force
the prescribed marker columns up to common phase. Singular coordinates or
one fixed observable do not. A common phase cancels between an actual word
and its actual inverse; the literal compiler convention is stronger.
The changed-decoder state-based result is a different interface.

The QBP mean and gradient perturbation bounds follow from the same
full-output contract. Dirty-bank reuse permits conditional means rather
than independent shots, with fresh declared initialized inputs for each
execution. The Hilbert martingale bound used for concentration allows
bounded increments without independence or conditional symmetry. Sampling,
observable work, initialized interference branch, supplied evaluators,
and classical preprocessing remain charged separately. These facts imply
no general end-to-end gradient advantage.

## Imported premises checked

These are source-interface checks, not a claim that all the resulting
complete-frame arguments were already in the literature. The
[source map](SOURCE_MAP.md) retains the full attribution.

| Primary source and locator | Admitted premise and boundary |
|---|---|
| [Yuan–Zhang, v3](https://arxiv.org/html/2202.11302v3), Theorem 2; Lemmas 5, 6, 9; Section 2 | Exact state-preparation comparison and complete-unitary toolbox with returned clean work. The state theorem alone does not prescribe the remaining columns. The ancillary-free MCX primitive uses arbitrary exact one-qubit phases; it is not imported into the Clifford+T ledger. |
| [Low–Kliuchnikov–Schaeffer, v2](https://arxiv.org/html/1812.00954v2), Eq. (8), Table 2, Fig. 1(d), Appendix C | Coherent dirty-bank cancellation and arbitrary-output XOR queries. Initialized selector/output work remains distinct from borrowed banks. |
| [Gosset–Kothari–Wu, v3](https://arxiv.org/html/2411.04790v3), Lemma 2.3, Theorem 4.2, Lemma B.1 and Corollary B.2 | Logarithmic determinant-one synthesis, adaptive ancillary diagonal lower bound, and separately phase-correct complete multiplexor synthesis. The repository supplies its real-frame reduction and operator-to-channel comparison. |
| [Khattar–Gidney, v1](https://arxiv.org/html/2407.17966v1), Sections 4, 5.4 and Fig. 6 | Logarithmic-depth MCX with exactly two arbitrary dirty helpers, returned coherently and disjoint from outer toggle controls. Exact-Toffoli substitution preserves the asymptotic bounds without measurement uncomputation. Section 5.4 becomes 5.5 in v2. |
| [Low–Wiebe](https://arxiv.org/pdf/1805.00675), Appendix B Lemma 13; [Fang–Lin–Tong, v2](https://arxiv.org/html/2208.06941v2), Section 2.4 Lemma 3 and Appendix D | Logarithmic coherent product history retains private block workspace. The repository proves its own failure-counter convention and charges that work. |
| [Berry et al.](https://arxiv.org/pdf/1412.4687), Eqs. (11)–(15) | Cubic oblivious amplification lineage. The complete initialized-isometry estimate, leakage, and literal signs are proved locally. |
| [Bausch](https://arxiv.org/pdf/2009.10709), Eqs. (4), (6), Section 2.3.3 | Geometric weights and digit-source lineage. The logarithmic peak initialized width is an additional local construction. |
| [Pinelis](https://arxiv.org/pdf/1208.2200), Theorem 3.5 | Hilbert-space bounded-increment martingale tail; the stated QBP increment bound gives denominator 32 in the exponent. The separate conditional-symmetry hypothesis of Theorem 3.6 is not needed. |

## Correction, evidence, and stopping decision

The exact compiler's Section 8 formerly called its displayed tail-work
expression the peak. It is a conservative **upper bound**: at
$`n=2,t=1,s=1`$, four work qubits suffice while the expression gives five.
The enclosing eight-qubit reservation remains valid. Changing the sentence
to “at most” fixes the overstatement without changing the theorem.
No other claim-level repair was identified in this bounded review.

The reviewed baseline passes 505 tests, four exact fault-tolerant receipt
suites, and rendered-presentation checks. Those finite checks support their
specified contracts; the decisions above rest on the proof arguments and
source premises. This pass adds no circuit, simulation, or new theorem.
The [verification map](VERIFICATION.md) continues to distinguish native
small words, exact identities, illustrative ledgers, and asymptotic proofs.
It does not promise a general grouped elementary emitter or useful numerical
crossover constants for existential grouping thresholds.

The scientific stopping gate is satisfied for the selected package.
Keep the high-precision constant-clean count gap and the large-width depth
gap explicit. Both remain optional research. A discovered defect reopens its
repair; a new construction requires the native rule and symbolic resource,
error, and cleanup recurrences in the
[workspace gate](../WORKSPACE.md#reopening-research-requires-a-qualifying-mechanism).
Manuscript writing remains on hold until requested. Further activity alone
is not a reason to enlarge the scope or repeat finite fixtures.
