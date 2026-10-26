# Source and dependency map

[← Verification](VERIFICATION.md) · [Complete narrative](../REVIEW.md) · [Related work →](RELATED_WORK.md)

This is the claim-to-source map for Results A–D, including C's
one-clean banked frontier and precision-uniform subtree depth bounds. It separates
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
bcc7333c94cd6bd5f1c4abeb559594b7777b6b0b
```

and is reconciled in [SYNC](../SYNC.md) and
[upstream provenance](../provenance/upstream.json).
The upstream [paper crosswalk](https://github.com/GoGoKo699/Hopf-QBP/blob/bcc7333c94cd6bd5f1c4abeb559594b7777b6b0b/review/README.md)
locates its complete scientific statements and executable checks.

Q1–Q7 fix the complete inverse-frame record, shared magnitude parities,
separate phase stream, and phase-calibrated controlled access to a Hermitian
unitary observable. The statistical target is simultaneous absolute accuracy
of the raw coordinate gradient at a fixed ideal parameter tuple. Complete-vector,
relative, normalized-frame, and natural-gradient targets have different
conditioning. A checkpoint needs its active interface, not merely one state
column. See the [QBP register](reference/SOURCE_CATALOGUE.md#3-inherited-qbp-interface)
and [matched-program consequence](QBP_CONSEQUENCE.md).

## 4. Results established in this repository

The [publication scope](../manuscript/PUBLICATION_SCOPE.md) owns the
selected formulas, precision ranges and literal reservations. The
[resource frontier](OPEN_PROBLEM.md) distinguishes those claims from
the strongest research upper bound, R96, and the remaining depth question.
R79 remains selected Result D. This page maps dependencies rather than
restating their resource ledgers.

| Claim | Local construction | Imported premises and boundary |
|---|---|---|
| A: exact complete frame, R1–R11 | [Compiler theorem](COMPILER_THEOREM.md): borrowed-suffix echo, conditioned tree cut, decoder, router and local lower-bound reductions | C1–C5 and H1–H9; prescribed real and phase-dressed complex magnitude frames, exact workspace return |
| B: sufficient-clean T-count, R12–R16 | [Fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md): compact source, complete residual corrections and shared-source composition | F1–F5 and F10–F11; sufficient clean allocation, complete targets and elementary Clifford accounting |
| C: one-clean count, R64–R65 | [Uniform growing-source theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count) and [banked theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier), with R25's small-dimension baseline | F2–F5, F8–F9, F11 and F19–F20 through R17/R20/R25/R62/R63; exact-base and doubled-base matching have different hypotheses; no T-depth consequence |
| D: one-clean simultaneous count/depth, R79 | [Uniform subtree theorem](../research/depth/UNIFORM_SUBTREE_DEPTH.md), using [R78's native spectral decoder](../research/depth/SUBTREE_PROJECTOR_DECODER.md) and R72's one-clean terminal, width/precision splice and literal diagonal | Existing descendant-source, native query, recovery and conditional-work premises; all per-call precision, protected storage, actual inverses and leakage are charged; no unrestricted depth optimum |
| Stronger research upper, R96 | [All-precision linear-depth theorem](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md): borrowed precision, fast references and directly programmed terminal | Existing exact bounded-control, source, recovery and query premises; same one-clean real/complex reservations and full operator contract; no new imported premise, efficient-preprocessing or general-emitter claim |

The [stable R register](reference/SOURCE_CATALOGUE.md#4-results-established-in-this-repository)
retains the full statements, separate state-based QBP results and earlier
resource regimes. Upper bounds hold for all admitted parameter tuples;
matching lower bounds are worst-case over the named family. The exact and
Clifford+T constructions do not assert simultaneous optimality of one
circuit in both gate models.

### Count construction dependencies

| Result | Dependency and local contribution |
|---|---|
| R62: two-clean high-precision endpoint | [Mixed local gauges](../research/endpoint/HETEROGENEOUS_PRECISION_SOURCE.md), [fixed-bottom disjoint blocks](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks), [disjoint state extensions](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources), R20's conditional guards and R25's dirty scalar words supply the native interfaces. The [growing-source proof](../research/endpoint/GROWING_SOURCE_COMPILER.md) supplies growing arenas, exact refinement, compensated overlap and summable charges. |
| R62: complex endpoint and lower bound | R17's packed literal diagonal gives the [same-wire complex corollary](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint); the [single-angle witness](FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness) applies to both families. This retains the sharper two-clean complex reservation. |
| R63: uniform precision | The [two-clean proof](../research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count) separates logical source size from accuracy. Its [high-precision block/state variant](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#precision-parametrized-high-precision-form) and R17's low-precision/bounded-dimension baseline cover complementary regimes. F3/F4 match the exact base width only under the stated precision condition. |
| R64: one-clean recovery | The [whole-source conjugate lift](../research/endpoint/CONJUGATE_SIGNAL_RECOVERY.md#2-lift-the-whole-high-level-word) borrows the phase signal outside every core while retaining the initialized recovery anchor. Both conjugate sectors give the same real target; actual dirty replay and leakage are included. R25 supplies literal complex phase composition. |
| R65: banked count | [Banked high-precision interfaces](../research/endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces), self-borrowed selectors and exact returned word banks preserve complete masks. R25 covers low precision and bounded dimension, including one and two logical qubits. F3/F4 supply matching count bounds; the arbitrary-diagonal complex subfamily retains the canonical-domain lower bound. |
| R52–R53: separate reductions | The [bounded-factor bridge](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md) and [regular tree reduction](../research/endpoint/TREE_CAYLEY_REDUCTION.md) retain their factor-coverage and Cayley-synthesis hypotheses; they are not substitutes for the complete-frame theorem. |

R62–R65 reuse existing synthesis and lower-bound premises. Their complete
proofs price guards, every source, actual inverse, dirty/reference action
and work return; neither the count theorem nor its banked extension alone
supplies a depth bound.

### Depth construction dependencies

R66 is a retained two-clean, linear-precision predecessor, not the current
general frontier. Its [full theorem](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
keeps its own range, source proofs and explicit resource ledger.

| Obligation | Source and local proof boundary |
|---|---|
| Coarse real-state approximation | [GKW v3, Lemma 3.2 and Section 3.3](https://arxiv.org/html/2411.04790v3#S3): Boolean-phase preparation and LCU/amplification architecture; its ancillary allocation is not imported |
| Coherent descendant source | [Loader Sections 2–5](../research/depth/COHERENT_DESCENDANT_LOADER.md#2-a-strictly-certifiable-geometric-expansion) and [complete loader](../research/depth/COHERENT_DESCENDANT_LOADER.md): strict certified recursion, preserved tag/prefix sectors, calibrated rotations, exact inactive identity and leakage-aware source error |
| Native access | [Shallow dirty queries](../research/depth/SHALLOW_DIRTY_QUERIES.md): arbitrary-input XOR/phase queries, private non-Clifford supports under shared read-only addresses, charged masks and guards |
| Early shared precision | [Shared-gradient recovery](../research/depth/SHARED_GRADIENT_RECOVERY.md): one paid preparation, compute/swap/uncompute arithmetic and tables, ideal-gradient invariant sector, and reflection excluding the live core; F32 supplies the retained per-call calibration adder |
| Complete-column recovery | [R66's column construction](../research/depth/SUBLINEAR_T_DEPTH.md#3-compress-the-complete-columns): depth tags, square-root normalization, full-unitary comparison, literal amplification signs and terminal cutoff |
| Terminal source and history | [Shared binary terminal](../research/depth/SHARED_BINARY_TERMINAL.md): R12/R14 compact source and accepted-shift composition; F1/F10/F11 geometric weights, coherent failure history and LCU recovery; a sixteen-type one-layer dictionary and protected suffix |
| Conditional controls and programs | [R20's bounded-control lemma](CONDITIONAL_SUFFIX_COMPILER.md#exact-conditioning-without-a-clifford-count-penalty) and [whole-word prefetch](../research/depth/GROUPED_TERMINAL_DEPTH.md#2-prefetch-the-exact-coarse-native-word): literal guarded phases, stable addresses and actual unloading |
| Deepest layers | [Rectangular nonuniform query](../research/depth/NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence), capped operator source and [uniform baseline](UNIFORM_PRECISION_DEPTH.md): simultaneous dirty reservations, complete tail error and fallback |
| Separate interfaces | [Grouped terminal](../research/depth/GROUPED_TERMINAL_DEPTH.md) retains its column-forest proof; the [entry-indexed limit R67](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md) concerns its stated row/column factorization, not unrestricted depth |

R79 uses the original per-call descendant loader and R78 spectral decoder
instead of the early shared-gradient/overlap-recovery route. Its terminal
uses R72's one-clean relocation and uncapped deepest-layer interface.
R45, R66, R70 and R71 retain their original two-clean statements; the
terminal binary source and early binary phase gradient are distinct
charged mechanisms, never supplied precision states.

### Lower-bound dependencies and scope

| Result | Imported or retained input | Local contribution and limit |
|---|---|---|
| [Clean-aware R73](../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md) | F35/F44 stabilizer facts; F43 in the retained scalar-closure proof only | Fixed-clean three-layer obstruction; [R74's initialized-sector proof](../research/depth/FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector) independently makes it explicit without F43 |
| [Crossing-rank R74](../research/depth/FOUR_LAYER_CROSSING_RANK.md) | F44 symplectic realization and stabilizer coordinates | Active-span normalization and explicit fixed-clean rank-dependent gap; unbounded crossing rank is not excluded |
| [Operator rigidity R75](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md) | F35/F44, R73 local phase gap and unit-circle arithmetic | Two-layer Pauli-image isolation and exact commutants yield the zero-clean five-layer lower bound; at one clean qubit, only terminal clean-Z preservation follows |
| [Midpoint R76](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md) | R75 phase-polynomial extension and R74 finite-boundary form | Compressed-commutant lift for Pauli/Hermitian-Clifford midpoints; separate local cyclotomic norm argument for finite-projective-order blocks; terminal return alone does not imply either hypothesis |
| [Inner spread R77](../research/depth/ONE_CLEAN_INNER_SPREAD.md) | R74 finite-boundary theorem and R76 lifting/phase discreteness | Flat Pauli expansion covers bounded inner spread; its exact common-sign-flip counterexample limits that certificate, not unrestricted rigidity |
| [Outer overlap R97](../research/depth/OUTER_OVERLAP_RIGIDITY.md) | R74 estimate, R76 lifting/closure and R77 flat expansion | Physical-axis normalizer, shared-image compatibility and diagonal-Clifford extraction give a bounded-misalignment gap; [representation limits](../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md) distinguish same-space descent from non-invariant center mixing |

The [frontier](OPEN_PROBLEM.md#constant-depth-frontier) owns the exact
thresholds and quantifiers. These statements do not prove unrestricted one-clean four-layer rigidity or
dimension-growing optimal depth.

## 5. Fault-tolerant sources and contribution boundaries

The following supporting-family constructions reuse existing premises.
Their local promises and ranges remain in their proof homes, not in a
second general resource ledger.

| Supporting result | Retained ingredients and local addition |
|---|---|
| [Common-flag R80](../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md) | One-qubit synthesis, conjunction/query interfaces, [quadrant-sign gauge](../research/endpoint/PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization), and R72/R79 fallbacks. Actual shared product factors give exact return; certified common-angle perturbation preserves the promise. |
| [Recursive terminal R81](../research/depth/RECURSIVE_TERMINAL_DEPTH.md) | R78 full-mode decoder and rectangular query with unchanged external prefix. The uniform-height resource sums are local; the early decoder must still be charged. |
| [Bounded-density R82](../research/depth/BOUNDED_DENSITY_DEPTH.md) | Precision gauge, native real-vector source, R80 common flag, sharp queries and [R66 isometry extension](../research/depth/SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison). Finite variance/Carleson bounds, sparse Schur overlaps and recursive density control are local. |
| [Positive-product R82 corollary](../research/depth/BOUNDED_DENSITY_DEPTH.md#12-positive-product-reference-corollary) | Effectively supplied strictly positive classical reference and fixed density-ratio bounds; constants are independent of biases. Reference zeros and efficient reference search are not covered. |
| [Affine orientations R83](../research/depth/AFFINE_ORIENTATION_DEPTH.md) and [deterministic frames R84](../research/depth/DETERMINISTIC_FRAME_DEPTH.md) | Existing sign gauge/queries plus paid triangular affine leaf maps and short-tag marker corrections, or exact selected involutions, signed strings and last-deviation inverse; no free in-place permutation. |
| [Geometric-band R85](../research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md) | GKW-based sources, R78 reflections, native queries/echoes, sign/phase composition and the logarithmic-height R81 terminal. Coarse bins, upper-frame unloading, local flag transport and a temporarily clean dirty output word are proved locally. |

None of these scoped families replaces the strongest general upper R96.
The detailed [source register](reference/SOURCE_CATALOGUE.md) distinguishes
local proof additions from imported synthesis or lower-bound premises.

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

The [complete F1–F54 register](reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries)
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

### Research-only group-algebra premise

The [span-preserving five-layer proof](../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#54-span-preserving-comparisons-a-finite-degree-bound)
uses the classical Jennings dimension-subgroup formula and weighted
augmentation basis for finite 2-groups. The precise statements are reproduced
at their point of use and verified against H. Usefi,
[Identifications in Modular Group Algebras](https://www.math.mun.ca/~usefi/papers/JPAA.pdf),
p. 2 and Corollary 3.4; Lemma 5.1 also verifies the subgroup-adapted basis
construction. The filtration-extension bound, native nilpotence and Boolean
degree argument are proved locally. These algebraic premises add no circuit
oracle or resource assumption and are not premises of selected Results A–D.

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
