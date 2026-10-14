# Source and dependency map

[← Verification](VERIFICATION.md) · [Complete narrative](../REVIEW.md) · [Related work →](RELATED_WORK.md)

This is the claim-to-source map for Results A–D, including C's
one-clean banked frontier and the fixed-accuracy large-width depth theorem. It separates
imported synthesis tools, inherited Hopf/QBP interfaces, and the operator
constructions established here. The [source catalogue](reference/SOURCE_CATALOGUE.md)
preserves the complete stable C/H/Q/F/R register and detailed contribution
boundaries. Placement in that register does not select a research study as
a premise of the core theorem.

The [publication scope](../manuscript/PUBLICATION_SCOPE.md) fixes the exact
parameter ranges and corollaries. The exact and Clifford+T theorems price
the same prescribed completion in different gate models; they do not
assert simultaneous optimality of one circuit in both models.

## 1. Imported state-preparation toolkit

The normative exact compiler source is P. Yuan and S. Zhang,
“Optimal (controlled) quantum state preparation and improved unitary synthesis
by quantum circuits with any number of ancillary qubits,” *Quantum* **7**,
956 (2023), doi:10.22331/q-2023-03-20-956. The published article corresponds
to `arXiv:2202.11302v2`; the imported statements were also checked in
`arXiv:2202.11302v3`.

| Stable ID | Precise imported premise | Role in A |
|---|---|---|
| C1 | Theorem 2: exact state-preparation size and all-workspace depth frontier | Comparison scale; the local theorem still proves the prescribed completion |
| C2 | Lemma 5: ancilla-free multi-controlled X, linear size/depth | Borrowed-suffix predicates and controlled branches |
| C3 | Lemma 6: complete UCG synthesis at arbitrary clean width | Narrow echo UCGs, subtree frames, and phase diagonals |
| C4 | Lemma 9: coherent CNOT copy–uncopy | Prefix decoder and coherent router |
| C5 | Theorem 1: generic controlled state preparation | Comparison with a generic all-column construction |

The earlier Sun–Tian–Yang–Yuan–Zhang paper is the historical predecessor.
The active proof uses the uniform Yuan–Zhang framework at every workspace
budget. [Full source locators and local consumers](reference/SOURCE_CATALOGUE.md#1-imported-state-preparation-toolkit)
retain the individual imported resource statements.

## 2. Inherited Hopf operator interface

The [Hopf interface](HOPF_INTERFACE.md) specifies the balanced binary chart,
state column, marker columns, and addressed layers (H1–H9). The oriented
incoming amplitude satisfies

```math
\partial_{\theta_j}|\psi\rangle=a_j|e_j\rangle,
\qquad g_{j,j}=a_j^2.
```

On the canonical domains, $a_j\ge0$. At a singular magnitude coordinate
the raw derivative vanishes while the parameter tuple selects a unit
marker-frame continuation. The complex magnitude frame is
$W_{\mathbb C,\mathrm{mag}}=D_{\mathrm{ph}}W_{\mathbb R}$; leaf-phase
derivatives use their own stream. These are inherited geometry/interface
facts, not consequences of the compiler's finite matrix checks.

The [H1–H9 register](reference/SOURCE_CATALOGUE.md#2-inherited-hopf-operator-interface)
links each fact to its local restatement and implementation.

## 3. Inherited QBP interface

The tracked `Hopf-QBP/main` baseline is

```text
a9885317cf998a7df87ca07ba86e3bd4f0f419ef
```

and is reconciled in [SYNC](../SYNC.md) and
[upstream provenance](../provenance/upstream.json).

Q1–Q7 fix the complete inverse-frame record, shared magnitude parities,
separate phase stream, and phase-calibrated controlled access to a Hermitian
unitary observable. The statistical target is simultaneous absolute accuracy
of the raw coordinate gradient at a fixed ideal parameter tuple. Complete-vector,
relative, normalized-frame, and natural-gradient targets have different
conditioning. A checkpoint needs its active interface, not merely one state
column. See the [QBP register](reference/SOURCE_CATALOGUE.md#3-inherited-qbp-interface)
and [matched-program consequence](QBP_CONSEQUENCE.md).

## 4. Results established in this repository

Write $N=2^n$, $q=n+a+b$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ for $0\lt\eta\le1/64$.
Upper bounds hold for every parameter tuple; matching lower bounds are
worst-case statements over the specified family.

| Selected result | Local proof and construction | Imported premises and boundary |
|---|---|---|
| **A. Exact complete frame**: size $\Theta(N)$ and depth $\Theta(n+N/(n+m))$ for every clean budget | [Compiler theorem](COMPILER_THEOREM.md): strict-zero echo, conditioned tree cut, decoder, explicit router, and local lower-bound reductions; R1–R11 | C1–C5 and H1–H9; real and phase-dressed complex magnitude frames; exact work return. CNOT count is $\Theta(N)$ for $n\ge2$ and zero for $n=1$. |
| **B. Matching finite-precision frontier**: $T^\star=\Theta(\sqrt{NL}+L+NL/q)$ | [Fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md): compact source, complete corrections, shared-source composition; R12–R16 | F1–F5, F10–F11; sufficient clean reservation $a\ge C(n+h)$, $h=1+\lceil\log_2(L+n+2)\rceil$; real and complex magnitude families; $O(NL)$ Clifford cost. |
| **C. One-clean count frontier**: $`T=O(N+L)`$ at base width; $`T^\star=\Theta(\sqrt{NL}+L+NL/b)`$ at doubled base width | [One-clean growing-source proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count) and [banked proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier); R64–R65, with R25's small-dimension corollary | Every $`n\ge1`$, $`L\ge6`$, fixed $`a=1`$; base $`B=L+n+7`$ real, $`B=L+n+8`$ complex magnitude; banked match at $`b\ge2B`$; $`G=O(NL)`$. Exact-base matching at $`L\ge\max\{6,n\}`$. Existing F2–F5, F8–F9, F11, F19–F20 and R17/R20/R25/R62/R63 interfaces; full isometry, literal phases and actual inverses. |
| **D. Same-circuit count and T-depth** | [Uniform-precision theorem](UNIFORM_PRECISION_DEPTH.md), with retained [unary source](UNARY_PHASE_GRADIENT.md) and dirty-query constructions; R43–R45 | Two clean flags and $b\ge17(L+n+7)$; complete real frames. Native sources, inverses, predicates, queries, and conditional work are charged. |
| **Fixed-accuracy width-sensitive depth**: $`D_T=O_\eta(N/b^2+\sqrt n\log(n+2))`$ | [Compressed-column construction](../research/depth/SUBLINEAR_T_DEPTH.md), [coherent descendant loader](../research/depth/COHERENT_DESCENDANT_LOADER.md), and [shallow dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md); R66 | Every fixed $`0\lt\eta\le1/64`$, two clean flags, and $`b\ge17(L+n+7)`$; complete real frames, same-circuit $`T=O_\eta(\sqrt N+N/b)`$ and $`G=O_\eta(N)`$. GKW's Boolean-state lemma and TTK's exact modular addition are imported; coherent loading, guarded phase-gradient calibration, native query scheduling, width splice and full-operator recovery are local. The linked theorem retains the stronger minimum with n. No optimal-depth claim beyond its matching interval or precision-uniform improvement. |
| **Endpoint count theorem and complex-magnitude corollary**: worst-case $`T^\star=\Theta(N)`$ | [Growing-source proof](../research/endpoint/GROWING_SOURCE_COMPILER.md) and [complex-magnitude composition](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint); R62 | Existing F2/F5/F8–F9/F11/F19–F20 ingredients through the retained native interfaces; $`a=2,b=N+n+7,\eta=2^{-N},L=N,n\ge3`$. Complete prescribed real and phase-dressed complex magnitude frames; literal phases and full work-return error. No additional Clifford-count or T-depth bound. |
| **Precision-uniform two-clean count**: $`T=O(N+L)`$ | [Uniform growing-source proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count); R63 | Existing R17/R20/R25/R62 native interfaces, with the precision-parametrized disjoint-block/state-source form; $`n\ge3`$, integer $`L\ge6`$, $`\eta=2^{-L}`$, $`a=2`$, $`b\ge L+n+7`$ for real frames and $`b\ge L+n+8`$ for phase-dressed complex magnitude frames. Universal constants, full initialized-isometry error, and actual inverses; no new Clifford-count or T-depth bound. |

The separate [entry-indexed encoding limit](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
(R67) gives a sharp square-root normalization barrier within its explicit
row/column factorization contract. It is not a native-depth lower bound.

For D, the uniform upper bounds are

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad
D_T=O\!\left(\frac{NL}{b^2}+nL\right),\qquad G=O(NL).
```

For $6\le L\le\log_2(n+2)/16$, replace the count's $nL$ term by zero
and the depth's $nL$ term by $n$. Simultaneous worst-case matching holds on
$17(L+n+7)\le b\le\sqrt{N/n}$ for every $L\ge6$, and on
$17(L+n+7)\le b\le\sqrt{NL/n}$ in the stated slow-precision range,
whenever the intervals are nonempty. Outside those ranges the upper bounds
retain their actual precision terms; no general T-depth optimum is claimed.

The fixed-accuracy large-width theorem has a separate dependency chain;
it does not change the uniform theorem's literal thresholds or matching
intervals:

| Obligation | Source and local proof boundary |
|---|---|
| Coarse real-state approximation | [GKW v3, Lemma 3.2 and Section 3.3](https://arxiv.org/html/2411.04790v3#S3): Boolean-phase state preparation and LCU/amplification architecture; its ancillary allocation is not imported |
| Coherent descendant source | [Loader Sections 2–5](../research/depth/COHERENT_DESCENDANT_LOADER.md#2-a-strictly-certifiable-geometric-expansion): local strict certified residual recursion, preserved tag/prefix sectors, calibrated real rotations, exact inactive identity, and source-isometry error including leakage |
| Native data access | [Shallow dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md): returned arbitrary-input XOR/phase queries, private non-Clifford supports under shared read-only addresses, and charged mask/guard words |
| Complete-column recovery | [Depth construction Sections 3–6](../research/depth/SUBLINEAR_T_DEPTH.md#3-compress-the-complete-columns): depth tags, normalization $`\sqrt K`$ with $`K=O(n)`$, full-unitary comparison, literal amplification signs, and polylogarithmic terminal cutoff |
| Count-efficient final layers | [Retained nonuniform indicator](../research/depth/NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence) and capped operator source: linear indicator count, separate simultaneous dirty reservations, and complete tail error |

The [stable R register](reference/SOURCE_CATALOGUE.md#4-results-established-in-this-repository)
also records supporting corollaries, separate state-based QBP theorems, and
research constructions. Their individual hypotheses are not suppressed by
the selected claims.

R62 uses the [mixed local gauges](../research/endpoint/HETEROGENEOUS_PRECISION_SOURCE.md),
[disjoint-block compilation of the fixed bottom levels](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks),
and [disjoint state extensions](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources),
together with R20's conditional-suffix interface and R25's actual dirty
scalar words. The local proof supplies the growing arenas, exact source
refinement, compensated overlap, and summable native charges. R17's packed
literal diagonal gives the complex-magnitude corollary on the same wires;
the [single-angle precision witness](FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
gives $`2N-1`$ for both families. This composition imports no new synthesis
premise or literature-priority claim.

R63 separates the dyadic source precision from the logical dimension:
the retained arenas depend on N, while the source charges sum to $`O(L)`$.
The [high-precision disjoint-block/state-source variant](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#precision-parametrized-high-precision-form)
supplies the same complete native interfaces for $`L>N/n^2`$ beyond a
universal dimension threshold. R17 supplies both the $`L\le N/n^2`$
splice and the bounded-n, arbitrary-L fallback. For complex magnitude
frames, precision $`L+1`$ in the real factor and R17's packed literal
diagonal give the stated $`L+n+8`$ reservation. No additional external
premise is imported, and R62 retains its tighter complex endpoint width.
For R63 at exactly $`b=L+n+7`$ (real) or $`b=L+n+8`$ (complex
magnitude), $`L\ge\max\{6,n\}`$ makes physical width $`q=\Theta(L)`$.
The existing F3/F4 lower bound $`\Omega(NL/q+L)`$ therefore gives
worst-case $`\Theta(N+L)`$ at these widths. The uniform upper bound
also allows larger dirty allocations; this matching corollary does not.

R64 uses the [whole-source conjugate lift](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md#2-lift-the-whole-high-level-word)
before outer recovery. The borrowed phase signal stays outside every core;
the initialized anchor marker and logical recovery coordinates remain.
The same real target is obtained in both conjugate sectors, with actual
dirty replay and leakage included. R25's one-clean literal diagonal then
composes the complex magnitude frame at precision $`L+1`$.

R65 uses the [banked high-precision interfaces](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces)
and the retained exact dirty word-bank query. Self-borrowed selectors and
returned banks implement the same complete masks. Their summed native
charges remove the repeated precision term, and the full Clifford ledger
is $`O(NL)`$. At the stated bank thresholds, $`q=n+1+b=\Theta(b)`$,
so F3/F4's existing $`\Omega(\sqrt{NL}+L+NL/q)`$ bound matches the
upper construction. C uses R25's layerwise banked baseline for
$`L\le N/n^2`$, where $`nL\le\sqrt{NL}`$, and for bounded n.
The cases $`n=1,2`$ have $`N\le4`$, so that baseline costs $`O(L)`$
and the root-angle witness supplies the matching $`\Omega(L)`$.
R20's exact conditional-work guards remain a dependency; its grouped
precision bound is a separate retained corollary. R64–R65 introduce no
new external synthesis or lower-bound premise, priority claim, or
T-depth guarantee. For complex magnitude frames, the arbitrary-diagonal
subfamily supplies the lower bound even on the canonical angle domain.

R17 includes the packed two-clean diagonal reservation. The endpoint
[bounded-factor bridge](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
and [regular tree reduction](../research/endpoint/TREE_CAYLEY_REDUCTION.md)
are R52–R53. Their factor-coverage and Cayley-synthesis hypotheses distinguish
these constructive reductions from the complete-frame upper bounds.

## 5. Fault-tolerant sources and contribution boundaries

| Premise | Primary locator | Imported role and local specialization |
|---|---|---|
| Dirty lookup and finite-width counting | LKS, Section 2, Eq. (8), Appendices B.2/C, Section 5; F2–F3 | Exact dirty-bank echo and counting method; local literal interpreters and frame reductions |
| Optimal T-count benchmarks | GKW, Theorems 1.1–1.2, 4.1–4.2 and Appendix B; F4–F5 | Unrestricted count, diagonal lower bounds, complete multiplexors, and error allocation; prescribed frame and small-clean contracts remain local |
| Boolean-phase real-state preparation | [GKW v3, Lemma 3.2 and Section 3.3](https://arxiv.org/html/2411.04790v3#S3) | Coarse real-state approximation and constant-query LCU structure; the coherent descendant extension, strictly certifiable selection, dirty implementation, and native depth ledger are proved locally |
| Precision source | Bausch, Eqs. (4), (6), Section 2.3.3; F1 | Geometric weights and bit oracle; local compact and dirty-operator source implementations charge peak workspace |
| Conditional work and predicates | Khattar–Gidney, Sections 3–4, 5.4, 7.3; F8 | Two returned dirty helpers for logarithmic-depth MCX; helper lifetimes are disjoint from live query storage |
| Full-space loaders | Kerenidis–Prakash, Definitions 4.4/4.6, Theorem 4.9; F9 | Standard anticommuting-loader representation; local geometric source and programmed block words |
| Product compression and amplification | Low–Wiebe Lemma 13; Fang–Lin–Tong Lemma 3/Appendix D; Berry et al. Eqs. (7)–(15); Brassard et al. Eq. (8); F10–F11/F20 | Failure history, LCU, and amplification; local complete-isometry error and literal-phase accounting |
| Dirty-query depth and unary reuse | F29–F34/F36–F37 | Controlled Clifford circuits, grouped selection, phase synthesis, arithmetic, reusable references, and bilinear recursion; local all-dirty return and width allocation |

The [complete F1–F42 register](reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries)
provides primary URLs, versions, exact locators, and comparison-only sources.
The [related-work map](RELATED_WORK.md) separates established ingredients
from the frame-specific construction; bounded source checking certifies
neither priority nor novelty by absence of a matching search result.

With $J_a$ appending initialized clean work, every finite-precision claim uses

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

This covers arbitrary logical and dirty/reference inputs, literal phases,
and all final clean leakage. The sufficient-clean construction returns dirty
work exactly. Operator-source and one-clean cores return approximately within
this norm; their lookup banks, selectors, and suffix-control subroutines
return exactly. Grouped conditional-suffix work and its clean predicate also
have final return error included in the norm. No measurements, resets,
postselection, or free supplied precision sources enter the stated model.

D retains [dirty-sum compression](DIRTY_SUM_COMPRESSION.md),
[chunked indicators](CHUNKED_DIRTY_INDICATOR.md), and
[parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md) as supporting constructions.
Standalone source-depth bounds and radial filtering are optional research
studies and supply no additional premise for its uniform theorem.

## 6. Evidence classification

Analytic identities, explicit reversible schedules, imported elementary
synthesis, and finite regression checks are distinct evidence levels.
The [verification map](VERIFICATION.md) connects the selected claims to their checks; the
[fixture catalogue](reference/VERIFICATION_CATALOGUE.md) retains every
individual input and software boundary. Matrix identities and receipt
arithmetic expose convention, phase, cleanup, and resource errors without
proving asymptotic optimality or supplying a general native frame emitter.

## 7. Provenance records

- [Upstream provenance](../provenance/upstream.json): commits, file lineage, and reconciliation.
- [Literature provenance](../provenance/literature.json): source versions and role assignments.
- [Synchronization policy](../SYNC.md): tracked baseline and reconciliation rules.

These records preserve attribution; they are not additional scientific premises.
The [borrowed-workspace appendix](BORROWED_WORKSPACE_COMPILER.md) contains
the full argument used here. Its historical sources,
`Hopf_Fault_Tolerant_Research.md` §§4.5–4.6, 4.8, 5.19 and
`Hopf_Small_Clean_Workspace.md`, identify lineage rather than extra proof
or runtime dependencies.
