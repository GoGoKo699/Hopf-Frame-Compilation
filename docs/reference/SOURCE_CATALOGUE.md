# Source catalogue and extended lineage

[Core source map](../SOURCE_MAP.md) · [Core related work](../RELATED_WORK.md) · [Verification catalogue](VERIFICATION_CATALOGUE.md)

This reference catalogue preserves the stable source and result identifiers
and the detailed contribution boundaries behind the concise core maps.
Its register includes selected Results A–D, supporting constructions,
separate state-based QBP results, and research studies. Inclusion here is
not selection as a principal result or a claim that every listed source is
a premise of the current compiler.

## Catalogue index

| Material | Location | Reading status |
|---|---|---|
| C1–C5: exact synthesis imports | [Imported toolkit](#1-imported-state-preparation-toolkit) | Core premises for A |
| H1–H9 and Q1–Q7: inherited interfaces | [Hopf](#2-inherited-hopf-operator-interface) and [QBP](#3-inherited-qbp-interface) | Shared operator and output contracts |
| R1–R11: exact compiler results | [Exact results](#4-results-established-in-this-repository) | Core and supporting results for A |
| F1–F42 and R12–R51: finite-precision register | [Fault-tolerant register](#5-fault-tolerant-sources-and-contribution-boundaries) | Mixed register; each row retains its actual scope |
| Detailed literature comparison | [Extended lineage](#extended-lineage-for-the-core-and-its-supporting-constructions) | Attribution and comparison detail |
| Developmental depth and endpoint studies | [Research literature](../../research/RELATED_WORK.md) | Research record; not additional principal results |
| State-based QBP literature | [Separate supplement](../../supplements/state_based_qbp/RELATED_WORK.md) | Separate task and decoder contracts |

The F/R identifiers are stable and are not renumbered by document placement.
The authoritative selection is [publication scope](../../manuscript/PUBLICATION_SCOPE.md).
In particular, standalone source-depth obstructions and radial-filter
variants are optional studies, not premises of Result D. Their register
entries preserve attribution and the restricted statements actually proved.

## 1. Imported state-preparation toolkit

The normative compiler citation is:

> P. Yuan and S. Zhang, “Optimal (controlled) quantum state preparation and
> improved unitary synthesis by quantum circuits with any number of ancillary
> qubits,” *Quantum* **7**, 956 (2023),
> doi:10.22331/q-2023-03-20-956.

The published article corresponds to `arXiv:2202.11302v2`.  The imported
statements below were also checked in `arXiv:2202.11302v3`; their conclusions and
exact arbitrary-one-qubit+CNOT model are unchanged for the uses made here.

| ID | Imported result | Exact use in this repository | Local consumer |
|---|---|---|---|
| C1 | Theorem 2: exact arbitrary state preparation has $\Theta(2^n)$ size and $`\Theta\!\left(n+\frac{2^n}{n+m}\right)`$ depth for every clean-workspace budget | comparison frontier and target lower-bound scale | [compiler theorem](../COMPILER_THEOREM.md), [`unified_compiler.py`](../../compiler_robust_hopf/unified_compiler.py) |
| C2 | Lemma 5: exact ancilla-free multi-controlled X has linear size and depth | strict-zero toggles, direct suffix flags, and branch predicates | [`strict_zero_echo.py`](../../compiler_robust_hopf/strict_zero_echo.py), [`router.py`](../../compiler_robust_hopf/router.py) |
| C3 | Lemma 6: a total-width-$`q`$ UCG with $w$ clean work qubits has $O(2^q)$ size and $`O\!\left(q+\frac{2^q}{q+w}\right)`$ depth | half-angle UCGs, flagged layers, subtree frames, and the phase diagonal | [`strict_zero_echo.py`](../../compiler_robust_hopf/strict_zero_echo.py), [`unified_compiler.py`](../../compiler_robust_hopf/unified_compiler.py) |
| C4 | Lemma 9: coherent CNOT-tree copying and exact uncopying | control fanout in the conditioned prefix and coherent router | [`tree_decoder.py`](../../compiler_robust_hopf/tree_decoder.py), [`router.py`](../../compiler_robust_hopf/router.py) |
| C5 | Theorem 1: generic controlled state preparation depends on combined index and target width | comparison showing why a generic all-column construction would have quadratic Hilbert-space scale | [`resource_bounds.py`](../../compiler_robust_hopf/resource_bounds.py), [related work](../RELATED_WORK.md) |

The preceding state-preparation paper by Sun, Tian, Yang, Yuan, and Zhang is
retained as the historical source of the earlier time–space landscape and of
selected primitives credited by the later article.  The active proof uses the
uniform all-workspace framework above rather than selecting between the two
papers by regime.

## 2. Inherited Hopf operator interface

The compiler construction does not depend on the complete optimization
framework.  It consumes the following facts.

| ID | Inherited fact | Source role | Local restatement or implementation |
|---|---|---|---|
| H1 | a balanced complete binary tree with $N=2^n$ leaves gives $N-1$ real magnitude coordinates | chart definition | [Hopf interface](../HOPF_INTERFACE.md), [`frames.py`](../../compiler_robust_hopf/frames.py) |
| H2 | the sine–cosine path map covers normalized real states; leaf phases give the complex chart | forward state map | [`real_tree_data`](../../compiler_robust_hopf/frames.py), [`complex_state`](../../compiler_robust_hopf/complex_analysis.py) |
| H3 | real nonfinal angles use $[0,\pi/2]$, the final real depth uses $[0,2\pi)$, and complex magnitudes use $[0,\pi/2]$ | canonical chart domains | [`canonical_magnitude_angle_mask`](../../compiler_robust_hopf/frames.py), [`test_frames.py`](../../tests/test_frames.py) |
| H4 | the oriented incoming amplitude obeys $\partial_{\theta_j}\lvert\psi\rangle=a_j\lvert e_j\rangle$ and $g_{j,j}=a_j^2$ | differential geometry | [`RealTreeData.incoming_amplitude`](../../compiler_robust_hopf/frames.py), [`test_complex_analysis.py`](../../tests/test_complex_analysis.py) |
| H5 | on the canonical domains $a_j\geq0$, so $a_j=\sqrt{g_{j,j}}$ | canonical-domain consequence | [`test_frames.py`](../../tests/test_frames.py) |
| H6 | if $g_{j,j}=0$, the raw differential vanishes while the parameter tuple still selects a unit marker-frame continuation | singular-coordinate boundary | [`regular_coordinate_mask`](../../compiler_robust_hopf/frames.py), singular tests in [`test_frames.py`](../../tests/test_frames.py) |
| H7 | the marker $\lambda(j)$ places the state and coordinate-frame directions in one known computational basis | differential-frame interface | [`conventions.py`](../../compiler_robust_hopf/conventions.py), [`frames.py`](../../compiler_robust_hopf/frames.py) |
| H8 | the complete real frame is a product of prefix-selected rotations restricted to the zero-suffix sector | addressed-frame structure | [Hopf interface](../HOPF_INTERFACE.md), [`direct_addressed_depth_layer`](../../compiler_robust_hopf/frames.py) |
| H9 | the phase-dressed complex magnitude frame is $W_{\mathbb C,\mathrm{mag}}=D_{\mathrm{ph}}W_{\mathbb R}$ | complex magnitude interface | [`complex_magnitude_frame_matrix`](../../compiler_robust_hopf/frames.py) |

The complete inverse map and optimization context are contained in the first
Hopf paper and repository.  The addressed frame, global magnitude record,
direct phase record, checkpoint interface, and statistical task boundaries are
contained in `Hopf-QBP`.

## 3. Inherited QBP interface

The tracked `Hopf-QBP/main` baseline is

```text
a9885317cf998a7df87ca07ba86e3bd4f0f419ef
```

and is reconciled in [`SYNC.md`](../../SYNC.md) and
[`provenance/upstream.json`](../../provenance/upstream.json).

| ID | QBP fact used here | Role in the compiler consequence | Local counterpart |
|---|---|---|---|
| Q1 | the global magnitude circuit applies the complete inverse differential frame | fixes the required compiler contract | [frame-safe compilation](../FRAME_SAFE_COMPILATION.md) |
| Q2 | one X-basis outcome contributes a parity record to every magnitude coordinate | shared-record mechanism | [`decoders.py`](../../compiler_robust_hopf/decoders.py), [`test_decoders.py`](../../tests/test_decoders.py) |
| Q3 | leaf-phase derivatives use a separate signed one-hot stream | separates the complex magnitude frame from the phase record | [`decoders.py`](../../compiler_robust_hopf/decoders.py), [`complex_analysis.py`](../../compiler_robust_hopf/complex_analysis.py) |
| Q4 | the primary finite-shot target is simultaneous absolute accuracy of the raw coordinate gradient | defines the $O(\log n)$ execution statement | [QBP consequence, Section 5](../QBP_CONSEQUENCE.md#5-statistical-target) |
| Q5 | complete-vector, relative, normalized-frame, and natural-gradient targets have different conditioning | prevents overextension of the raw-coordinate claim | [QBP consequence, Section 5](../QBP_CONSEQUENCE.md#5-statistical-target) |
| Q6 | the validated core assumes phase-calibrated controlled access to a Hermitian unitary observable | access-model premise | [QBP consequence, Section 2](../QBP_CONSEQUENCE.md#2-controlled-observable-interface) |
| Q7 | a checkpoint compiler needs equality on its active interface, not merely one state column | separate checkpoint boundary | [compiler boundaries](../COMPILER_BOUNDARIES.md) |

## 4. Results established in this repository

| ID | Result | Analytic proof | Implementation and finite evidence |
|---|---|---|---|
| R1 | complete frame safety implies substitution and, at regular points with exact prepared state, is necessary up to common phase for all observable-dependent fixed-decoder means | [frame-safe compilation](../FRAME_SAFE_COMPILATION.md) | actual probability-decoder tests with complex phases, workspace leakage, and singular markers |
| R2 | one correct state column need not preserve the decoded gradient | [two-qubit argument](../COMPILER_BOUNDARIES.md#2-two-qubit-global-state-column-counterexample) | [`test_compiler_boundaries.py`](../../tests/test_compiler_boundaries.py) |
| R3 | the borrowed-suffix echo implements one nonfinal addressed layer with no ancillary wire | [compiler theorem, Lemma Z](../COMPILER_THEOREM.md#lemma-z-borrowed-suffix-echo) | [`strict_zero_echo.py`](../../compiler_robust_hopf/strict_zero_echo.py), [`test_strict_zero_echo.py`](../../tests/test_strict_zero_echo.py) |
| R4 | the strict-zero frame has $\Theta(N)$ size and $\Theta(n+N/n)$ depth | [compiler theorem, Proposition Z](../COMPILER_THEOREM.md#proposition-z-strict-zero-resources) | [`strict_zero_audit.py`](../../compiler_robust_hopf/strict_zero_audit.py), exact-rational tests |
| R5 | the prefix is a conditioned smaller frame and the tail is a direct sum of subtree frames | [compiler theorem, Lemmas T1–T2](../COMPILER_THEOREM.md#6-exact-tree-cut) | [`tree_structure.py`](../../compiler_robust_hopf/tree_structure.py), cut tests |
| R6 | a clean binary–one-hot decoder has $3\cdot2^t-2-t$ workspace, $O(t)$ depth, and $O(2^t)$ size | [compiler theorem, Lemma P](../COMPILER_THEOREM.md#lemma-p-clean-binaryone-hot-decoder) | [`tree_decoder.py`](../../compiler_robust_hopf/tree_decoder.py), [`test_tree_decoder.py`](../../tests/test_tree_decoder.py) |
| R7 | an explicit coherent router realizes the tail direct sum and clears all data, token, copy, and flag work registers | [compiler theorem, Lemma R](../COMPILER_THEOREM.md#lemma-r-explicit-coherent-router) | [`router.py`](../../compiler_robust_hopf/router.py), [`test_router.py`](../../tests/test_router.py) |
| R8 | the maximal feasible cut attains $`O\!\left(n+\frac{N}{n+m}\right)`$ depth for every large workspace, including $s=1$ | [compiler theorem, Proposition R](../COMPILER_THEOREM.md#proposition-r-maximal-cut-depth) | [`resource_bounds.py`](../../compiler_robust_hopf/resource_bounds.py), broad-grid tests |
| R9 | parameter capacity and output light cones give matching size/depth bounds; fusing free one-qubit slots gives the separate CNOT lower bound for arbitrary clean workspace | [compiler theorem, Section 9](../COMPILER_THEOREM.md#9-matching-lower-bounds), using the standard parameter method of [Iten et al., Section III](https://arxiv.org/html/1501.06911v4#S3) | analytic dimension proof and resource diagnostics |
| R10 | the real frame and phase-dressed complex magnitude frame attain the all-workspace optimum | [compiler theorem, main theorem](../COMPILER_THEOREM.md#main-theorem-optimal-exact-hopf-frame-compilation) | [`unified_compiler.py`](../../compiler_robust_hopf/unified_compiler.py), both resource ledgers |
| R11 | frame-safe compilation preserves the global record and introduces no additional asymptotic depth factor in the matched program | [QBP consequence](../QBP_CONSEQUENCE.md) | decoder and boundary tests; reviewer walkthrough |

## 5. Fault-tolerant sources and contribution boundaries

The [fault-tolerant compiler](../FAULT_TOLERANT_COMPILER.md) uses a complete
initialized-isometry approximation, uniformly over logical and dirty inputs
and their references. In that sufficient-clean construction, dirty workspace
is restored exactly; initialized logical work may have residual error included
in the isometry norm. The operator-source constructions below, including their one-clean
refinement, instead include the dirty operator core's return error in that norm, while returning
lookup banks, selectors, and suffix-control work exactly. These are distinct
return guarantees within the complete-input approximation model, separate
from the exact clean-workspace size–depth theorem.

| ID | Source and locator | Imported fact or comparison | Boundary |
|---|---|---|---|
| F1 | [Bausch, arXiv:2009.10709v4](https://arxiv.org/pdf/2009.10709v4), Eqs. (4), (6), Section 2.3.3 | geometric precision weights, address-and-digit bit oracle, capped source | these ideas predate this project; the explicit exact preparation uses a unary temporary before binary conversion |
| F2 | [LKS, arXiv:1812.00954v2](https://arxiv.org/html/1812.00954v2), Section 2, Table 2, Figure 1(d), Eq. (8), Appendix B.2 and Appendix C Theorems 1–2 | dirty-assisted lookup, restored banks, swap routing, and parallel dirty indicators | general lookup and parallel-selector ideas are inherited; conjugating a single X by the exact bank router supplies the scratch-free indicator used here; partial batches reuse the dirty linear echo and bank decomposition, with a local smaller-workspace allocation, literal inverse/phase ledger, and Hopf composition |
| F3 | LKS, Section 5 | Clifford+T circuit counting at finite width | the frame packing and specialization of dirty inputs supply the local reduction; the counting method is imported |
| F4 | [GKW, published Quantum article](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/), Theorems 1.1–1.2 and 4.1–4.2 | unrestricted optimal T-count benchmarks and ancilla-independent lower bounds | generic state preparation does not prescribe the Hopf completion; the diagonal subfamily supplies the frame reduction |
| F5 | GKW, Lemmas 2.1, 2.3, B.1 and Corollary B.2 | Boolean synthesis, one-qubit approximation, complete multiplexors, geometric layer-error allocation | neither multiplexor synthesis nor geometric allocation is claimed as new; the direct composition retains a repeated precision charge |
| F6 | [Tan, published PRX Quantum article](https://journals.aps.org/prxquantum/pdf/10.1103/pxhd-9s9q), Definition I.2, Theorem I.1, Lemma IV.1 and Remark IV.2 | general-unitary comparison and shared Boolean instruction synthesis | the instruction registers in the displayed construction are initialized; this is not the prescribed small-clean/dirty tradeoff |
| F7 | [Li–Ou–Wang–Yao–Yuan–Zhang, arXiv:2607.28260v1](https://arxiv.org/html/2607.28260v1), Sections 3–4 | sparse QROM and sparse-state comparison | different input families; no general full-frame conclusion is imported |
| F8 | [Khattar–Gidney, arXiv:2407.17966v1](https://arxiv.org/html/2407.17966v1), Sections 3–4, 5.4 and 7.3 | conditional work, dirty selectors, and logarithmic-depth MCX with two dirty helpers | these are inherited techniques; frame predicates borrow idle query storage, while batch guards reserve two disjoint helpers throughout each live query; both use exact native Toffolis without measurements |
| F9 | [Kerenidis–Prakash, arXiv:2202.00054v2](https://arxiv.org/html/2202.00054v2), Definitions 4.4/4.6 and Theorem 4.9; [Chee et al., arXiv:2301.07477](https://arxiv.org/pdf/2301.07477), Appendix C; [Bravyi, arXiv:quant-ph/0404180](https://arxiv.org/pdf/quant-ph/0404180), Section II, Eqs. (2)–(5) | full-space Clifford loaders, scalar/antisymmetric product decomposition, and paired-Majorana rotations | the Pauli representation and overlap algebra predate this work; the operator-source notes supply native geometric specializations and dirty-programmed block words; paired Majoranas and their anticommutation are standard representations |
| F10 | [Low–Wiebe, arXiv:1805.00675v2](https://arxiv.org/html/1805.00675v2), Lemma 13; [Fang–Lin–Tong, arXiv:2208.06941v2](https://arxiv.org/html/2208.06941v2), Section 2.4, Lemma 3 and Appendix D | coherent product compression with logarithmic failure-history work | failure tracking is inherited; residual factorization, shared-source reuse and the charged frame schedule are the local specialization |
| F11 | [Berry–Childs–Cleve–Kothari–Somma, arXiv:1412.4687](https://arxiv.org/pdf/1412.4687), Eqs. (7)–(15) | prepared linear combinations of selected unitary blocks, normalization-two oblivious amplification, and its cubic accepted block | LCU and amplification are inherited; the local proofs explicitly price the selected branches and bound the complete initialized isometry, including rejected work retained coherently |
| F12 | [Yamazaki–Akibue–Sano, arXiv:2603.14202v3](https://arxiv.org/html/2603.14202v3), Theorem 1, Eqs. (21)–(24), (27), Section 3.1 | deterministic complete controlled SU(2) synthesis and probabilistic precision constants | precision-length clean instructions remain; the matching coefficient-three lower bound holds for fixed control width and independent Haar blocks |
| F13 | [Yuan–Zhang–Zi, arXiv:2608.17846v2](https://arxiv.org/html/2608.17846v2), Theorem I.1 and Definition II.1; [Fang–Heunen–Wang, arXiv:2607.12907v1](https://arxiv.org/html/2607.12907v1), Theorem 1.2 and Corollaries 3.8–3.9 | current generic-unitary and near-Clifford comparisons; the [endpoint substitution](../../research/ROUTE_HISTORY.md#selection-audit-after-canonical-completion) keeps precision factors explicit | different target families and clean allocations; the initialized-isometry contract already appears in general-unitary synthesis |
| F14 | [Vasconcelos–Gilyén, arXiv:2507.07900v2](https://arxiv.org/html/2507.07900v2), Sections 2–3 and Appendix B | block-work uncomputation, exact-history lower bounds and approximate product compression | original clean work is still needed during a query; the lower bound is for the defined coherent-measurement class, and approximate compression needs near-identity dilations |
| F15 | [Ma–Joven–Liu, arXiv:2609.11153v1](https://arxiv.org/html/2609.11153v1), Theorem 3.1 | clean ancilla compression for block-encoding counting bounds | at most $`n+2T`$ clean ancillas, with possible normalization change; no two-clean unitary conclusion |
| F16 | [Motlagh–Pocrnic, arXiv:2605.20334v1](https://arxiv.org/html/2605.20334v1), Section II | improved dirty-QROM constants | the displayed output word remains initialized; no asymptotic or constant-factor lookup improvement is claimed here |
| F17 | [Lai, arXiv:1411.5408v3](https://arxiv.org/pdf/1411.5408v3), Theorems 1.6 and 1.8 | finite-tree Carleson embedding and L2 maximal inequality | the inequalities are inherited; the tree-transport note supplies the local overlap-defect identity and its weighted residual application, not a free block-encoding circuit |
| F18 | [Boyd–Vandenberghe, Convex Optimization](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf), Appendix A.5.5 | quadratic elimination and Schur complements | standard linear algebra; the weighted-block and residual-assembly notes supply their scalar tree recursions, local rejection corrections, and mode allocations |
| F19 | [Barenco et al., arXiv:quant-ph/9503016v1](https://arxiv.org/pdf/quant-ph/9503016v1), Lemma 4.1 and Section 8 | Euler factors, two-level unitary decomposition, and Gray-code basis routing | inherited generic synthesis; the local pricing uses the retained borrowed-sector echo and does not claim an optimal native implementation |
| F20 | [Brassard–Høyer–Mosca–Tapp, arXiv:quant-ph/0005055](https://arxiv.org/pdf/quant-ph/0005055), Section 2, Eqs. (7)–(8) | repeated amplitude amplification and its sine/cosine law | the five-call one-clean amplification is a specialization at amplitude $`\sin(\pi/10)`$; the local proof supplies the oblivious complete-isometry error and literal-phase accounting |
| F21 | [Vatan–Williams, arXiv:quant-ph/0308006](https://arxiv.org/pdf/quant-ph/0308006), Section III Theorem 1 and Section IV | magic-basis equivalence between SO(4) and two SU(2) factors | the reduction and fixed-size real two-qubit synthesis are inherited; the local Cayley formulas fix a literal Clifford convention and reuse the existing one-target compiler with a borrowed logical spectator |
| F22 | [Beals–Buhrman–Cleve–Mosca–de Wolf, Quantum Lower Bounds by Polynomials](https://homepages.cwi.nl/~rdewolf/publ/qc/polynomials.pdf), Lemma 4.1 | query-by-query amplitude-degree growth | the [Hopf scattering audit](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer) derives its Laurent-degree specialization directly; it concerns unchanged local-angle queries, not native T counts |
| F23 | [Toulouse–Umrigar, arXiv:physics/0701039v2](https://arxiv.org/pdf/physics/0701039), Eq. (45); [Mitarai–Fujii, arXiv:1901.00015](https://arxiv.org/pdf/1901.00015), Fig. 1 | logarithmic wave-function derivatives and Hadamard-test interference | inherited estimator ingredients; the reference-state QBP note proves its own leaf scores, Hopf envelope, support bounds, and native resource comparison |
| F24 | [Möttönen–Vartiainen–Bergholm–Salomaa, arXiv:quant-ph/0407010v1](https://arxiv.org/pdf/quant-ph/0407010v1), Section III, Eqs. (4), (5), and (7) | uniformly controlled phase cascade and its residual arithmetic-mean common phase | inherited decomposition; the complex coarse proof fixes local signs and applies the retained exact-return dirty interpreter to determinant-one row words |

| F25 | [Ross–Selinger, arXiv:1403.2975v3](https://arxiv.org/abs/1403.2975v3), abstract and synthesis-runtime discussion | distinguishes short native words from efficient search; optimal synthesis uses a factoring oracle, and the efficient expected runtime without it is conditional | computational context only; the bounded-input proof uses F5 word-length existence and guarded exhaustive coarse enumeration, with no factoring oracle or runtime conjecture |
| F26 | [Baur–Strassen, *Theoretical Computer Science* 22(3), 317–330 (1983)](https://www.sciencedirect.com/science/article/pii/030439758390110X) | classical reverse differentiation background | the explicit Pauli baseline derives its Hopf reverse recurrence and dyadic error bound directly; reverse differentiation is not claimed as a new algorithmic principle |
| F27 | [Casas et al., arXiv:2602.05425v1](https://arxiv.org/html/2602.05425v1#S3.SS2.SSS2), Section III.2.2, Eq. (29) | exact denominator-exponent depth bound for Clifford-matchgate plus T layers | the [source-depth proof](../../research/depth/SOURCE_T_DEPTH.md) applies this existing method to two particular sources, permitting signed-permutation Clifford stages including the parity-odd extension; no unrestricted depth lower bound follows |
| F28 | [Vasconcelos, arXiv:2609.34659v1](https://arxiv.org/html/2609.34659v1#S3.SS4.SSS1), Theorem 8; [Kim, arXiv:2506.15147v3](https://arxiv.org/pdf/2506.15147), Sections 2.2, 2.4 and 3 | precision-depth and catalytic-rotation comparisons | the former uses precision-sized clean work; the latter uses a prepared catalyst, leaves the clean/dirty-only constant-depth question open, and records subsequent improvements in its final note; see F29 for the universal-catalyst refinement |
| F29 | [Kim–Laakkonen, arXiv:2512.24982v1](https://arxiv.org/html/2512.24982v1), Theorems 3, 5 and 6; Section 5.1 | constant-depth control of CNOT/Clifford circuits without ancillas, constant-factor depth overhead for controlled Clifford+T, and universal catalytic rotations | the controlled-shear lemma is a rank-sensitive specialization; the catalytic comparison uses classically compiled angle-dependent matrices and measurement-assisted dynamic preparation, not the unary source's coherent one-hot program and charged unitary boundary pair |
| F30 | [Boyd, arXiv:2312.00696v2](https://arxiv.org/html/2312.00696v2), Section III and Appendix A | commuting SELECT/QROM groups and Clifford changes of basis for parallel action | its address copies use initialized registers; the local all-dirty shear echo and workspace allocation are proved separately |
| F31 | [Selinger, arXiv:1210.0974v2](https://arxiv.org/pdf/1210.0974), Section 2, Eqs. (5)–(6), Theorem 4.1; Proposition 5.1, Eqs. (17)–(21) | parallel parity-phase synthesis with initialized work, and one-T-layer Pauli-conjugation normal form | parallel phase synthesis is established; the unary-source proof separately gives guarded trilinear phases on conditional work using the retained native Toffoli word; the rational-transfer specialization concerns full-input nonaffine permutations with returned dirty helpers, not initialized-clean isometries |
| F32 | [Takahashi–Tani–Kunihiro, arXiv:0910.2530v1](https://arxiv.org/pdf/0910.2530v1), Sections 2.1–2.3 | linear-size, linear-depth exact ripple-carry addition without initialized work | deleting the two gates targeting the arbitrary carry-output wire gives the modular adder used in the signed dirty increment and baseline sum tree; its literal CNOT/Toffoli word is emitted in the bounded checks; later clean-work/fanout constructions are not used |
| F33 | [Remaud–Vandaele, arXiv:2501.16802v2](https://arxiv.org/html/2501.16802v2), Lemmas 2/4, Algorithm 3, Theorem 2 | exact helper-free addition via shallow CNOT and Toffoli ladders | truncate the carry output at the abstract ladder level, then synthesize the shorter ladders; applies only to private counters; bounded checks audit the reduced macro, while the optimized ladder-depth bound is imported analytically |
| F34 | [Vandaele, arXiv:2603.12917v1](https://arxiv.org/html/2603.12917v1), Section 5, Theorem 4 and Corollary 7 | exact logarithmic-depth increment and controlled increment with one returned dirty helper | supplies the retained round-based compressor and the separate two-dirty-bit read-only increment; temporary control borrowing stays on private supports; the carry-pipeline refinement instead uses linear TTK arithmetic |
| F35 | [Aaronson–Gottesman, arXiv:quant-ph/0406196v5](https://arxiv.org/pdf/quant-ph/0406196v5), Section III; [Zhang–Zhang, arXiv:2409.13809v2](https://arxiv.org/html/2409.13809v2#S3.SS1), Theorem III.1, Eqs. (10)–(11) | stabilizer-overlap quantization and Pauli conjugation by one T layer into a Hermitian Clifford | the [two-layer source obstruction](../../research/depth/SHALLOW_SOURCE_OBSTRUCTION.md) derives a full-space transfer alphabet and robust source witnesses; initialized-clean isometries and growing frame-depth lower bounds are excluded |
| F36 | [Jones et al., arXiv:1204.0567](https://arxiv.org/pdf/1204.0567), Section 2.1, Eqs. (2)–(4); Section 4.1, Fig. 15 and Eqs. (19)–(20) | Fourier-state phase kickback, programmed shifts, and reusable phase references | the unary source changes the encoding and implements coherent one-hot-selected shifts in conditional Hopf work; neither phase kickback nor reference reuse is new |
| F37 | [Iggy van Hoof, arXiv:1910.02849v2](https://arxiv.org/html/1910.02849v2), Section 4.2; Section 4.3, Theorem 4.1 | standard three-product Karatsuba recursion and subquadratic binary-polynomial multiplication | the unary-source proof uses the bilinear rank recursion followed by linear cyclic reduction; it does not import constant depth from the space-efficient reversible multiplication schedule |
| F38 | [Wu et al., arXiv:2609.36574v1](https://arxiv.org/html/2609.36574v1), Sections II.C, III.B and IV.C–E; Theorem 1 and Proposition 1 | shared phase arithmetic, reversible encoding, and amortized reference preparation | comparison checked 3 October 2026; its binary ripple-adder construction permits measurement/feedforward and has linear-width depth; no coherent unary-program, two-clean Hopf, or constant-depth convolution interface is imported |
| F39 | [Parham, arXiv:2504.19966v1](https://arxiv.org/html/2504.19966v1), Proposition 1.8, Theorems 1.14–1.15 and Section 6 | magic-hierarchy and classical-complexity comparisons | comparison only: the state/Boolean reductions do not establish a blanket barrier for prescribed-unitary lower bounds; no compiler premise or growing frame-depth lower bound is imported |
| F40 | [Al-Ghattas–Gamarnik–Kiani, arXiv:2610.02166v1](https://arxiv.org/html/2610.02166v1), Corollary 1.7, Lemma 4.4 and Proposition 4.5 | state-complexity comparisons with unrestricted Clifford blocks | comparison only: the multiple-round theorem requires linear total width; arbitrary-width extensions concern specific one-round classes and do not supply a growing bound at $`b\asymp\sqrt N`$ |
| F41 | [Nielsen–Chuang, arXiv:quant-ph/9703032v1](https://arxiv.org/pdf/quant-ph/9703032), pp. 1–2, Eq. (3), Result and Eqs. (6)–(10) | exact programmable-unitary orthogonality | attribution for R49's local cache-capacity boundary; no native depth lower bound |
| F42 | [Zhang–Tan–Kothari–Gosset–Gidney, arXiv:2609.39092v1](https://arxiv.org/html/2609.39092v1), Theorems 1–3, Section 5 and Appendix B, Proposition 16 | adaptive constant-overhead compilation at inverse-polynomial error and unitary linear-cost phase-gradient preparation | Theorem 3 is measurement-free but uses precision-sized initialized work; the [source-reuse audit](../../research/endpoint/SOURCE_REUSE_LIMITS.md#9-dirty-programs-and-phase-gradient-sources) checks dirty substitutions and catalyst return without importing a two-clean compiler |

Standard Pauli linear combinations, reversible arithmetic, and oblivious
amplitude amplification are used with their actual preparations and adjoints.
The operator and workspace arguments that instantiate them are part of the
local construction, not additional oracle assumptions.

The exact operator-source lower bound uses
[Gosset–Kliuchnikov–Mosca–Russo, Sections 2.3 and 4](https://arxiv.org/html/1308.4134v1),
for the Pauli-transfer denominator method. The local proof supplies its
source-specific witnesses and normalization through returned helpers.
The correlated-shot extension uses
[Pinelis, Theorem 3.5](https://arxiv.org/pdf/1208.2200v2), with bounded
Hilbert-space martingale increments; the concentration inequality itself
is inherited.

| ID | Local fault-tolerant result | Contribution and support |
|---|---|---|
| R12 | exact compact capped geometric preparation | simultaneous linear T count and logarithmic peak clean width, with ordinary binary labels and exact temporary return; refines the implementation of the established geometric source |
| R13 | addressed SU(2) sampling and direct frame composition | complete operator contract with small clean work; combines fixed Pauli atoms, geometric bit sampling, and charged dirty lookup |
| R14 | shared-source composition of frame residuals | the uniform-precision construction combines accepted-branch shifts and a retained source with the inherited compression gadget F10 and amplification F11; see [theorem and proof](../FAULT_TOLERANT_COMPILER.md) |
| R15 | matching T-count in the stated workspace regimes | the sufficient-clean [theorem](../FAULT_TOLERANT_COMPILER.md) plus the F3/F4 reductions; the separate all-clean-budget [corollary and borrowed-workspace proof](../BORROWED_WORKSPACE_COMPILER.md) also uses the arbitrary-budget compiler under its additional condition |
| R16 | fixed-parameter bounded-score robustness | [QBP approximation](../QBP_APPROXIMATION.md): complete complex gradient, reflection sums, rounded weights, correlated dirty-bank reuse, and quantum/classical budgets; no derivative of a compiled word or gradient-query optimum |
| R17 | two-clean operator-source compiler and corollaries | [proof](../OPERATOR_SOURCE_COMPILER.md): full-frame $O(N+nL)$ bound and dirty-bank refinement; packed literal diagonals use $`b=n+1+\lceil(\ell+4)/2\rceil`$ for error $`2^{-\ell}`$ and $`T=O(N+\ell)`$, with matching banked count at twice that base reservation. General one-target U(2) multiplexors retain their separate reservation. F2/F9/F11 and R25's paired sign masks supply the source/query/amplification ingredients; core return error is included. |
| R18 | exact geometric operator-source costs | [source proof](../OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers): $2m-4$ uncontrolled and $2m-2$ controlled T gates, even with returned helpers; native words and transfer witnesses checked for small $m$ |
| R19 | sufficient-clean complex-frame extension | [Corollary 7](../FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame): literal diagonal SU(2) embedding, same-pool composition, and diagonal-subfamily lower bounds |
| R20 | conditional-suffix compiler | [proof](../CONDITIONAL_SUFFIX_COMPILER.md): $`O(N+L\ell_*(n))`$ T gates and $`O(NL)`$ Clifford gates for $`n\ge1`$, $`L\ge6`$, $`a=2`$, $`b\ge L+n+7`$; star residuals and streamed coarse programs use $`O(\log(s+1))`$ private clean work per group of $`s`$ depths |
| R21 | two-clean T-depth upper bound | [schedule](../T_DEPTH_COMPILER.md): $`D_T=O(NL/b+L\ell_*(n)+n^3)`$, $`T,G=O(NL)`$, at $`b\ge2(L+n+7)`$; combines inherited F2 routing and F8 predicates with the complete-frame source, coarse-program, and conditional-work ledgers; no matching depth claim |
| R22 | structured residual and scoped endpoint diagnostics | [source-reuse analysis](../../research/endpoint/SOURCE_REUSE_LIMITS.md) and [tree transport](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md): linear generators, normalized transport, weighted norms using F17, finite-order correction witnesses, and Pauli-routed leakage; the [Frobenius-coarse specialization](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#frobenius-small-residuals-also-fit-the-endpoint-budget) is a direct consequence of the retained borrowed compiler, not a new frontier; the inherited Pauli-transfer method prices one paired-source transformed mask; a literal shared-conjugator fork has a coefficient-ellipse obstruction even under mask retuning, derived directly from the R25 Majorana moments; these are scoped interface restrictions, not additive or unrestricted frame lower bounds |
| R23 | simultaneous two-clean count and depth | [parallel dirty lookup](../PARALLEL_DIRTY_LOOKUP.md): a single X conjugated by the F2 router gives a scratch-free indicator with linear address T-depth; F8 predicates and the full-frame ledger reduce the fixed-accuracy bound to $`O(n^2)`$; complete real-frame composition retains count-efficient banks under a sufficient square-root-scale width condition; no new general lookup tradeoff or optimal T-depth claim |
| R24 | forward weighted transport component | [weighted block](../../research/endpoint/WEIGHTED_TRANSPORT_BLOCK.md): exact weighted Gram witness, scalar recursion using F18, complete one-signal-flag dilation, and certified singular-case preprocessing; F19 and the borrowed reflection interpreter give $`O(L\sqrt N)`$ native T count with exact dirty return, while its borrowed-signal extension gives $`O(N+nL)`$ with one block signal and approximate dirty return; routing is charged, and neither improves the full-frame endpoint bound |
| R25 | one-clean operator-source compiler and corollaries | [proof](../ONE_CLEAN_COMPILER.md): conjugated scalar-source word and Pauli-routed anticommutator, native head-and-two-tail source, and fixed five-call amplification using F20; real grouped/banked bounds extend to $`a=1`$, as do phase-dressed complex magnitude frames by sequential composition at their separate threshold; X symmetry gives a zero-clean layerwise real-frame corollary; matched literal-diagonal and complete U(2) multiplexor bounds retain their own reservations; F2/F9 supply inherited lookup and Clifford-algebra ingredients, and no new T-depth claim is made |
| R26 | complete two-clean residual assembly | [proof](../../research/endpoint/RESIDUAL_ASSEMBLY.md): affine forward dilation using F18, native controls using F19 and the borrowed-signal compiler, and two-term selection/amplification using F11 give $`T=O(N+nL)`$, $`G=O(NL)`$ at $`a=2`$, $`b\ge L+n+7`$, including leakage and dirty return; the full-mode hierarchy and elementary coupled completion have complete recursive contracts; direct block algebra gives a rank-at-most-four commutator repair without a mode-gap assumption, exact child-call cancellation, and a weighted repair-error bound; the surviving local wrappers retain their precision cost, so no grouped endpoint improvement or priority claim for the elementary completion follows |
| R27 | native-baseline antichain correction | [proof](../../research/endpoint/ANTICHAIN_COMPILER.md): exact descendant-forest conjugation and prefix-pair packing reduce promised antichain-only changes to one addressed SU(2) table; the retained borrowed interpreter and three Euler factors, using F19 and the borrowed-signal refinement of R25, give $`T=O(N+L)`$, $`G=O(NL)`$ at $`a=0`$, $`b\ge L+n+7`$; the supplied baseline has determinant-one native local words of length $`O(n-d+1)`$ at depth d, and the target agrees with it literally outside the antichain; neither generic rounding, a general complex-magnitude extension, nor closure of the unrestricted endpoint is claimed |
| R28 | native-baseline sparse nested updates | [proof](../../research/endpoint/SPARSE_UPDATE_COMPILER.md): exact off-support forest factorization and charged basis packing reduce an ancestor-closed set S to $`m=\lvert S\rvert+1`$ modes; for $`s=\lceil\log_2m\rceil`$ and $`n\ge2s+32`$, a shared dense column-residual dictionary and the R17 scalar-source block, using the inherited lookup, routing, and normalization-two amplification ingredients F2/F11/F19, give $`T=O(N+L)`$, $`G=O(NL)`$ at $`a=1`$, $`b\ge L+n+7`$; the complete initialized-isometry bound includes approximate work return and arbitrary dirty references; all changes on one root-to-leaf path are covered with a finite small-n fallback, but the unrestricted frame endpoint remains open |
| R29 | native source-width boundaries | [source-width proof](../../research/endpoint/SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation): the reversed geometric loader retains the coefficient grid and gives $`O(L+n)`$ T gates and $`O(NL)`$ Clifford gates for the outer loaders and monotone width bridges; the transformed group bodies remain separately charged; inherited Pauli-transfer and algebraic-norm methods price the original bridge and one mask from a legal grouped coefficient; [native checks](../../tests/test_precision_carry.py) support the literal phases, bridge direction, and full-input identities, without improving the compiler frontier |
| R30 | correlated flag/source boundary | [correlated-boundary proof](../../research/endpoint/SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost): charged initialization and decoding, direct code identities, and native physical width changes include arbitrary dirty references and released tail wires; a legal grouped mask produces code leakage, and complete coherent two-syndrome renewal reduces to the separately priced returned-source operation; [code checks](../../tests/test_correlated_precision_carry.py) include actual inverses and approximate flag leakage; this restriction concerns the full renewal interface, not accepted-only blocks or arbitrary carried-source compilers |
| R31 | small-product residual coordinates | [product-first examples](../../research/endpoint/SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis) reuse the existing fixed-address one-target compiler; [Cayley recursion](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) derives constant-size local data for the complete complex-coarse residual using standard Cayley/low-rank matrix algebra, with uniform conditioning and stable inverse conversion; [small-product](../../tests/test_small_product_compilation.py) and [Cayley](../../tests/test_tree_cayley.py) fixtures check complete matrices; no free coherent evaluator, improved native T-count, or priority for the classical transform is claimed |
| R32 | native four/eight-mode benchmarks | [four-mode primitive](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) combines F21 with the retained one-target compiler, certified factor coordinates, and a complete borrowed-spectator ledger; [eight-mode coupling](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling) gives a charged four-Pauli construction and a scoped product-factor distance; [inverse audit](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#9-inverse-cayley-recovers-the-original-local-scattering-word) reconstructs the original target wrapper; no generic endpoint improvement or additive precision lower bound follows |
| R33 | joint changing-target source word | [shared body](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) reuses R25's scalar source and amplification, with an unconditional fixed conjugator, complete SU(2) signal error, and a borrowed-signal parity return; direct word cancellation gives $`27g+6`$ source appearances, while fixed-mask tail commutation leaves both comparisons the same leading precision term; eight-/sixteen-mode allocations and [native checks](../../tests/test_joint_source_body.py) include queries, helpers, literal inverses, and complex coarse baselines; no generic T-count improvement or optimized-word minimum is claimed |
| R34 | canonical scalar compatibility with grouped SELECT | [group audit](../../research/endpoint/CANONICAL_SCALAR_COMPLETION.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge) instantiates R33 on the existing scalar flag, with a distinct atom flag and borrowed signal, full-operator perturbation, and the retained group workspace slack; direction-controlled Z selects actual inverses for both the canonical and old scalar words; [finite checks](../../tests/test_grouped_scalar_completion.py) cover the complete SELECT and amplification interfaces; the displayed precision recurrence remains linear in group count, with no improved endpoint or new lower bound |
| R35 | packed Hopf scattering and its boundary distinction | [tree-scattering proof](../../research/endpoint/ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer) gives complete root/marker ports, standard feedback elimination, and a charged $`O(N+L)`$-T step using R25 and exact borrowed routing; F22's degree method tests the unchanged-coin conversion via the Hopf path amplitude; [small checks](../../tests/test_hopf_scattering.py) verify full transfers and native port permutations; no priced feedback implementation or stronger native T lower bound follows |
| R36 | two-clean state-only preparation | [state-only proof](../../supplements/state_based_qbp/STATE_ONLY_COMPILER.md) combines the actual borrowed coarse circuit with one globally computed residual table, R25's full-operator SU(2) synthesis, and F20's one-step amplification at amplitude one half; it gives $`O(N+L)`$ T count for $`L\ge\max\{6,n\}`$ with exact stated dirty allocation and approximate joint return; no other system columns are prescribed |
| R37 | raw gradients from reference-state interference | [decoder proof](../../supplements/state_based_qbp/REFERENCE_STATE_QBP.md) specializes F23 to all real Hopf derivatives and all allowed observables, gives a division-free derivative-envelope recurrence and a bounded-branch corollary, and charges R36 plus the observable; generic sampling overhead can grow as n, and no complete-frame or total-gradient optimality follows |
| R38 | all-angle coarse-frame gradient decoding | [coarse-frame proof](../../supplements/state_based_qbp/COARSE_FRAME_QBP.md) combines R36's common-coarse reference preparation, F23's X/Y interference, and Walsh spreading with exact classical correction from the actual native C; the depth-record norm is at most five, with $`O(N+L')`$ T-count per execution apart from the observable and $`O(S+Nn)`$ classical arithmetic after preprocessing; no fine full-frame or total-gradient optimality claim |
| R39 | gauge-fixed complex state-based QBP | [native interface](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md) applies the borrowed reflection interpreter to F24's determinant-one phase tables, giving an actual logical coarse C with exact dirty return and $`O(N)`$ T count; R36's complex residual uses the same two flags and $`b\ge L+n+7`$; [two-stream proof](../../supplements/state_based_qbp/COMPLEX_COARSE_QBP.md) supplies magnitude and phase gradients with the R38 depth bound and $`O(S+Nn)`$ reconstruction; the literal original common frame phase and fine complete frame are not compiled |
| R40 | banked state preparation and fair QBP cost comparison | [banked state corollary](../../supplements/state_based_qbp/COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) reuses R17's exact whole-word dirty-bank query, inherited from F2, in R36/R39's preparation: $`T=O(\sqrt{NL}+L+NL/b+n\sqrt N)`$, $`G=O(NL)`$, at two clean flags, $`L\ge\max\{6,n\}`$, and $`b\ge2(L+n+7)`$; [gauged borrowed baseline](../../supplements/state_based_qbp/QBP_COST_COMPARISON.md#2-a-gauged-complex-borrowed-frame-baseline) combines the retained interpreter and F24's phase cascade with weighted row precision to give $`T=O(NK/(n+b)+K\sqrt N)`$, $`G=O(NK)`$, exact dirty return, and zero compiler clean work for $`b\ge2n`$; the [task comparison](../../supplements/state_based_qbp/QBP_COST_COMPARISON.md) retains separate frame/state precisions, literal workspace thresholds, observable costs, sampling, and preprocessing; it compares constructive upper bounds without a literal common-phase-frame or end-to-end gradient optimality claim |
| R41 | bounded-input construction and classical comparison | [algebraic residual coefficients](../../supplements/state_based_qbp/RESIDUAL_TABLE_PREPROCESSING.md) use standard half-phase identities with one shared root, certified rational intervals, and a finite-radius cutoff to preserve R36/R39's error constants without Euler search; [bounded-input audit](../../supplements/state_based_qbp/BOUNDED_INPUT_QBP.md) combines F5 existence with coarse enumeration, explicit masks and instruction output to prove polynomial construction for the listed grouped/state alternatives; the banked small-system source uses $`P+13`$ dirty wires. Deterministic and term-sampled classical Pauli baselines are charged; no unconditional efficient fine-word search, generic Euler-runtime theorem, or end-to-end quantum advantage is claimed |
| R42 | state-based QBP T-depth composition | [depth proof](../../supplements/state_based_qbp/STATE_QBP_DEPTH.md) composes R21/R23's exact schedules, inherited from F2, with R36/R39's constant number of residual rotations and the actual exact-return coarse interpreter; with $`B_0=P+n+7`$, $`b\ge2B_0`$ gives $`D_T=O(NP/b+P+n^3)`$, $`T,G=O(NP)`$, while $`b\ge16(B_0+\sqrt{NP})`$ gives one circuit with $`T=O(\sqrt{NP}+P+n\sqrt N)`$, $`G=O(NP)`$, and $`D_T=O(P+n^3)`$; both real/complex task streams and oracle depth are charged; no new lookup primitive, depth optimality, total-runtime gain, or general emitter is claimed |
| R43 | fixed-accuracy linear T-depth complete real frame | [unary phase-source proof](../UNARY_PHASE_GRADIENT.md), 3 October 2026: coherent one-hot shifts, guarded bilinear work, unitary source preparation/return, and the Hopf angle-stability/group/query allocation give $`D_T=O_\eta(n)`$, $`T=O_\eta(\sqrt N)`$, $`G=O_\eta(N)`$ with two external clean flags and sufficient $`C_\eta\sqrt N`$ dirty work; F5/F8/F31/F36/F37 are attributed ingredients, F29/F38 are comparisons; no supplied catalyst, generic synthesis priority, matching depth lower bound, or high-precision endpoint follows |
| R44 | fixed-accuracy dirty-width/depth tradeoff | [blocked bilinear lookup](../BLOCKED_BILINEAR_LOOKUP.md), composed with R43's early groups: for fixed $`0\lt\eta\le1/64`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and $`b\ge17(L+n+7)`$, one two-clean complete real-frame circuit has $`T=O_\eta(\sqrt N+N/b)`$, $`G=O_\eta(N)`$ and $`D_T=O_\eta(N/b^2+n)`$; simultaneous worst-case orders are $`\Theta_\eta(N/b)`$ and $`\Theta_\eta(N/b^2)`$ when $`b\le\sqrt{N/n}`$; the full-input leaf and allocation use existing F2/F8/F29/F30/F31 ingredients, and F3/F4 supply the lower-bound lineage; large-width depth optimality and the high-precision endpoint remain open |
| R45 | uniform precision and dirty-width depth tradeoff | [uniform composition](../UNIFORM_PRECISION_DEPTH.md), 3 October 2026: at two clean flags, $`L\ge6`$ and $`b\ge17(L+n+7)`$, absolute constants give $`T=O(\sqrt{NL}+NL/b+nL)`$, $`G=O(NL)`$, $`D_T=O(NL/b^2+nL)`$; for $`L\le\log_2(n+2)/16`$, the sharper same-circuit bounds omit nL from T and replace it by n in depth. The matching width endpoints are respectively $`\sqrt{N/n}`$ and $`\sqrt{NL/n}`$ when eligible; this is rectangular allocation and uniform composition of R43/R44 with the existing hybrid, using the same external premises and F3/F4 lower bounds, with no new native query or high-precision endpoint claim |
| R46 | linear-size dirty indicator with logarithmic address T-depth | [nonuniform chunk schedule](../../research/depth/NONUNIFORM_DIRTY_INDICATOR.md), 3 October 2026: an exact s-bit indicator has $`T,G,w=O(2^s)`$ and $`D_T=O(\log(s+2))`$, including arbitrary dirty-reference return. Unequal chunks reuse the existing tree echo and read-only conjunction; F2/F8/F31 retain their indicator, dirty-work and phase-synthesis roles. Composition with the existing rectangular bilinear query improves the late-tail ledger; early unary groups retain their $`O(n)`$ depth contribution. No new external premise, synthesis-priority claim or complete-frame depth order follows |
| R47 | one protected source across the early unary groups | [protected-source proof](../../research/depth/PROTECTED_UNARY_SOURCE.md), 3 October 2026: two initialized flags record the initial zero sector of a fixed logical source bank and each group's nonbank activity. One unconditional native U and its actual inverse give a global $`2\delta`$ source-error bound including final bank/flag return; inactive initial-bank sectors cancel exactly. With the stated reservations, the source-boundary depth is $`O(L+\log(n+2))`$, while three other early $`O(n)`$ contributions remain. F36 supplies phase-reference reuse, F8 conditional work, and F28/F29 remain catalytic comparisons; no supplied catalyst, new external premise, priority claim or improved full-frame depth order |
| R48 | cached activity predicates across early-group windows | [windowed-predicate proof](../../research/depth/WINDOWED_GROUP_PREDICATES.md), 3 October 2026: in the uniform low-precision regime, activity T-depth is $`O(n\log\log(n+2)/\log(n+2)+\log(n+2))`$. Windows use $`2J+1`$ additional conditional logical cache bits, $`J=\lceil\log_2(n+2)\rceil`$, the existing two clean flags and the same two returned dirty predicate helpers. The proof combines F8's borrowed conjunctions with consume-before-change cleanup and R47's protected bank; exact inactive-cache return preserves the global $`2\delta`$ source bound. Logical stages and program queries retain separate $`O(n)`$ allowances; no new external premise, priority claim or complete-frame depth order |
| R49 | shared-prefix query interface audit | [local audit](../../research/depth/SHARED_PREFIX_QUERY_AUDIT.md), 3 October 2026: separates reusable dirty work from initialized cache capacity and charges the shared-prefix interface. F41 contextualizes the exact local rank boundary. No improved complete-frame bound, unrestricted T-depth obstruction, or priority claim |
| R50 | charged two-group program refresh and its correction boundary | [shared-prefix audit, Sections 6–8](../../research/depth/SHARED_PREFIX_QUERY_AUDIT.md#6-a-charged-two-group-program-refresh): exact retained-program/source identity with actual inverse return and a reused activity flag; the generic difference family contains fresh queries, while explicit affine factors admit a constant-depth dirty echo. Existing bilinear and literal-Toffoli premises suffice; no new external premise, generic depth lower bound, or improved complete-frame order |

The [retained-source continuation](../../research/depth/RETAINED_SOURCE_FUSION.md#9-retained-source-fusion-without-a-larger-modulus)
extends R32 and R43 with the following component result.

| ID | Result | Scope and dependencies |
|---|---|---|
| R51 | original-modulus fusion and a shallow stabilizer completion | [exact source-ring factors](../../research/depth/RETAINED_SOURCE_FUSION.md#9-retained-source-fusion-without-a-larger-modulus) retain literal opposite determinant phases and every source input. The [growing even-height completion](../../research/depth/RETAINED_SOURCE_FUSION.md#10-two-source-shifts-for-the-stabilizer-completion) caches invariant last-nonzero-pair sectors and uses two selected shifts, with $`D_T=O(\log(g+2))`$, $`T,w=O(q2^g+R)`$, and $`G=O(q2^g+qR)`$, excluding program queries, outer predicates, and source boundaries. F21 supplies the magic basis; existing F8/F31/F36/F37-based primitives supply conditional work and source shifts. The remaining complete transport is still height-linear; no new external premise, optimality, priority, or improved complete-frame order is claimed |
| R52 | bounded-factor endpoint bridge | [Proof](../../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md): R17's packed source makes any fixed number of diagonals affordable at the literal endpoint width, with finite-dimensional fallback. Rational-circle enumeration gives terminating certified selection from exact coverage alone; uniform factor coverage remains a hypothesis. |
| R53 | regular fixed-tree Cayley reduction | [Proof](../../research/endpoint/TREE_CAYLEY_REDUCTION.md): literal phase gauge, polynomial resolvent, certified compact chart, exact O(N)-T coarse permutations, and a Haar-diagonal generator identity. R17 supplies exterior diagonals; the fixed-tree Cayley core still needs native synthesis. |

The [Hopf error audit](../HOPF_ERROR_ACCUMULATION.md) derives a sharp ideal-angle
stability recurrence and an exact finite relative-spectrum recursion from
the inherited nested frame supports H8, using standard block-matrix
spectral algebra. A midpoint-grid corollary is restricted to independent
nearest angular rounding. The audit also applies the existing scalar
source and amplification identities to an explicit family with coherent
linear leakage on reused flags. It does not import a general composition
lower bound or claim an unrestricted frame-depth obstruction. Its finite
fixtures use one common reduced source algebra, not a new native compiler.

The [four-echo audit](../../research/depth/HOPF_FLAG_ECHO.md) derives exact compressed and
full-isometry errors from those same source identities. The positive
[radial filter](../../research/depth/HOPF_RADIAL_FILTER.md) imports
[Grover's fixed-point phase sequence](https://arxiv.org/abs/quant-ph/0503205),
Eq. (1) and Section 3, and F5's determinant-one single-qubit word-length
theorem. Fixed-point amplification and ancilla-free phase approximation
are inherited. The local work supplies the literal global correction,
the full error including the residual accepted phase, a two-flag native
phase ledger, and the complete-frame precision cap. Failure-probability
suppression is not identified with phase-sensitive operator error, and
no new asymptotic count/depth frontier or fine-search runtime is claimed.

The [amortized dirty-lookup theorem](../AMORTIZED_DIRTY_LOOKUP.md) combines
F2/F8's dirty traversal and echo with F29/F30's controlled-linear and
commuting-operator techniques. Its query interface is retained by the
[variable-width hybrid](../PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision).
For every $`L\ge6`$, two clean flags and $`b\ge17(L+n+7)`$ give
one complete real-frame circuit with

```math
T=O(\sqrt{NL}+NL/b+nL),\qquad G=O(NL),
\qquad D_T=O(NL/b^2+nL+n\chi(n)),
```

```math
\chi(t)=\log_2(t+2).
```

Both count and depth match when
$`b\le\sqrt{NL/(nL+n\chi(n))}`$ and the interval is nonempty.
The older routed matching intervals remain valid separately.

The [bilinear query](../PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction)
specializes F2's indicator/bilinear framework. The
[dirty-counter indicator](../PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters)
uses exact sum and cyclic-rotation echoes. Its retained signed increment
applies two F32 adders and Clifford gates. The two-dirty-bit refinement
uses modular negation as an involution and imports F34's completed
controlled increment: public literals touch only CNOTs, while the source
borrows only private controls. F33 supplies the retained baseline
sum-tree adder.
The [masked carry pipeline](../DIRTY_SUM_COMPRESSION.md) uses F32's linear
TTK arithmetic for its growing blocks and final readout. Deferred parity
forests expose the carry predicates as linear forms in stable raw bits.
An exact private-helper phase identity lets overlapping blocks read those
bits without changing them. A bit at distance d becomes final within
$`O(d)`$ T layers; the column launch schedule therefore finishes in
$`O(m)`$ depth. A geometric resource sum and balanced-forest support
bound charge all dirty work and Clifford gates. The helper-only offsets
cancel in the outer echo. The new pipeline does not depend on F33/F34's
optimized depths; those remain analytic imports in their separate uses.
A width-dependent early/late cutoff absorbs polynomial indicator overhead
without changing the clean reservation, query error, or T-count order.
Reversible addition, phase polynomials, routing, and echo cancellation are
established ingredients; no generic synthesis priority or unrestricted
matching large-width depth is claimed. The separate two-bit indicator
certificate uses F31 for its full-input depth lower bound.

The [exact source-depth certificate](../../research/depth/SOURCE_T_DEPTH.md) specializes F27
to the existing geometric and paired sources, with a parallel paired-tail
schedule. Its lower bounds are for the stated Majorana-layer architecture
and exact targets; it changes no asymptotic frame or state-QBP theorem.
The [precision-depth comparison](../../research/RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026)
records F28's distinct workspace contracts.
The [two-layer source obstruction](../../research/depth/SHALLOW_SOURCE_OBSTRUCTION.md)
combines F35's standard facts with exact geometric-source coefficients.
The [conditional geometric source](../CONDITIONAL_GEOMETRIC_SOURCE.md)
instead supplies an explicit reversible prefix preparation and its
active-suffix block-encoding contract, using the existing exact lookup,
native controlled-H, reflection, and amplification ingredients. By itself
it changes the precision component only; it is not all-dirty source
resynthesis.

The retained [grouped program construction](../GROUPED_PROGRAM_PREFETCH.md) combines
that source with exact conditional prefetch, read-only program copies,
private conjunction trees, and an internal-enable phase correction.
The [chunked dirty indicator](../CHUNKED_DIRTY_INDICATOR.md) uses the existing
read-only dirty counter in a top-down reversible tree, then conjugates
a root flip to obtain an exact selected-path translation. Its late-layer
chunk allocation and the adaptive early grouping prove, at fixed accuracy
and sufficient square-root-scale dirty width, simultaneous
$`T=O(\sqrt N)`$, $`G=O(N)`$, and
$`D_T=O(n\log\log(n+2))`$ for one complete real-frame circuit with two
external clean flags. These are local compositions of the attributed
counter, bilinear query, conditional-work, and amplification interfaces;
they require no new external premise. Neither a matching large-width
depth bound nor the high-precision endpoint follows.

The [incremental selectors](../GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
reuse this same conditional-work interface and exact native Toffolis.
Their local dependency schedule reduces selector/enable maintenance to
linear group depth within the original reservation. The
[common-source audit](../../research/depth/COMMON_SOURCE_REUSE.md#11-a-common-source-identity-and-the-remaining-reflection)
uses direct unitary conjugation and the existing normalization-two block:
boundary preparations cancel, but conjugated reflections and source-bank
return remain charged. The legal zero-angle witness excludes only the
stated stale-monitor and one-use-bank substitutions. Neither statement
imports a new synthesis premise or improves the global depth order.

The [unary phase-source construction](../UNARY_PHASE_GRADIENT.md), added
**3 October 2026**, uses F36's phase-kickback identity with a coherent
one-hot program, F37's Karatsuba rank bound, and a separately proved
guarded bilinear circuit in conditionally initialized suffix work.
Parallel non-Clifford synthesis has the precedent F31. The actual source
preparation and inverse are charged once per group; source error is
bounded for the entire group, including rejected components and work
return. The Hopf angle-stability bound and early/late query allocation
then give R43's $`O_\eta(n)`$ T-depth with optimal-order T count at
fixed accuracy. This improves the preceding geometric-source schedule
without giving a shallow implementation of its conjugated reflection.
F29 and F38 give nearby catalytic/shared-arithmetic interfaces, with
different program, preparation, and workspace contracts. No priority
claim is inferred from these comparisons, and the high-precision
endpoint and unrestricted T-depth optimality remain open.

The [blocked bilinear query](../BLOCKED_BILINEAR_LOOKUP.md), also added
**3 October 2026**, extends that fixed-accuracy conclusion to every
$`b\ge17(L+n+7)`$ through R44. Its controlled leaf polarizes a
bilinear Boolean phase with a returned dirty bit; the completed leaves
then use the existing two-pass dirty traversal and bilinear indicator
echo. F2/F8/F29/F30/F31 retain their established query, traversal,
controlled-linear and phase-synthesis roles. The local proof adds the
literal full-input leaf, charged chunk allocation and complete-frame
composition, without a new external premise or generic priority claim.
R43 remains a valid construction and a square-root-width corollary.
The previous variable-accuracy theorem remains valid separately; R44's
constants depend on fixed eta. Its enlarged matching window follows
from the existing F3/F4 count lower bound divided by physical width,
not a new depth lower-bound method.

The [uniform precision refinement](../UNIFORM_PRECISION_DEPTH.md), added
**3 October 2026**, chooses rectangular blocks in R44's existing query
and makes the unary cutoff and the split between precision regimes
uniform. It reuses the same native circuits, source-return certificate,
and F2/F8/F29/F30/F31/F36/F37 ingredients; no new external synthesis
premise is introduced. R45's general count is optimal in order within
its matching interval, or at all eligible widths when $`L\le N/n^2`$;
its low-precision count is optimal throughout its stated regime. Both
matching intervals use the existing F3/F4 lower bounds. The earlier
theorems remain valid, and neither generic priority, unrestricted
large-width depth optimality, nor the selected high-precision endpoint
is claimed.

The [consolidated state-based QBP theorem](../../supplements/state_based_qbp/STATE_BASED_QBP_THEOREM.md)
collects R36 and R38–R42 under one input, precision, workspace, and
sampling contract. It introduces no additional compiler bound or
optimality claim. R37 remains an optional reference-state predecessor.

The [native coarse-frame example](../../supplements/state_based_qbp/NATIVE_COARSE_QBP.md) is finite implementation
evidence for R38 and R36's small-register fallback. Its commutator distance
certificate, per-reflection dirty echo, and literal gate ledger are stated
explicitly. It does not add a new asymptotic compiler theorem or replace
the existing analytic proof with numerical scaling evidence.

The [native complex extension](../../supplements/state_based_qbp/NATIVE_COMPLEX_COARSE_QBP.md) supplies finite
implementation evidence for R39: elementary all-suffix phase selection,
complete coherent preparation, both gradient streams, and a literal cost
comparison. Its sharper K2 trace and rational K3 certificate specialize
the same elementary commutator algebra. They do not change the asymptotic
compiler frontier or establish a general native residual-table emitter.

The [certified native residual rows and tables](../../supplements/state_based_qbp/NATIVE_RESIDUAL_ROTATION.md) implement
R25's paired-source rotations and borrowed-signal identity with R41's
algebraic coefficient intervals. Exact rational programming, literal
elementary words, and full-input finite checks close one implementation
bridge. The two-row extension uses exact Clifford mask selection and
R25's source-center conditioning with one enable literal; its inactive
sector is exactly identity. This control arity needs no extra helper.
These implementations add no asymptotic theorem: general lookup tables,
larger predicates, and the complete fine state compiler remain work.

For polynomial accuracy-bit budgets, the direct sampler can already attain the
matching T count. That regime is not attributed to the later shared-source
composition. The latter removes the repeated precision cost uniformly over
precision. Write $`\ell_*(n)=1+\log_2^*(n+2)`$, where $`\log_2^*`$
counts base-two logarithms until the value is at most one. Conditional-suffix
grouping and its one-clean refinement give the constant-clean endpoint upper bound $`O(N\ell_*(n))`$.
The gap from $`\Omega(N)`$ remains at the
[selected high-precision endpoint](../OPEN_PROBLEM.md).
At $`b\ge2(L+n+7)`$, the grouped compiler's banked form gives
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ T gates and $`O(NL)`$ Clifford gates.

The two-clean proof uses the
[two-pass dirty lookup](../BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
and [borrowed predicate toggle](../BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
proved in the borrowed-workspace appendix, together with amplification and
the lower bounds stated in the fault-tolerant chapter.
The conditional-suffix refinement resolves the sparse grouped-residual
support into forward and reverse rank-one stars and a diagonal. Depth labels,
marker predicates, and variable suffix Hadamards replace full endpoint words;
coarse programs stream one symbol at a time. Its private work is initialized
only in an explicitly selected logical sector. The dirty operator source
supplies precision without an initialized precision register.

The source audit supports these precise dependencies and comparisons. It does
not certify priority or infer novelty from a bounded search finding no match.
The current general-unitary benchmark includes Yuan–Zhang–Zi rather than
treating Tan as the newest result. The literal-diagonal theorem and its one-clean refinement
are separately compared with GKW and Yamazaki–Akibue–Sano in
[related work, Section 12](#12-contemporary-comparisons-and-the-broader-compiler-contribution).

## Extended lineage for the core and its supporting constructions

The following detailed comparisons retain their original section headings.
The concise [core account](../RELATED_WORK.md) selects their role in Results A–D;
nearby endpoint attempts and software boundaries remain supporting context.

## 1. Exact state preparation and the target frontier

The state-preparation line develops exact circuits whose depth decreases as
clean ancillary workspace increases.  Its uniform all-workspace conclusion is

```math
S_{\mathrm{QSP}}(n,m)=\Theta(2^n),
```

```math
D_{\mathrm{QSP}}(n,m)
=\Theta\left(n+\frac{2^n}{n+m}\right).
```

The active compiler framework in this repository is P. Yuan and S. Zhang,
*Quantum* **7**, 956 (2023).  We use its exact results for:

- ancilla-free multi-controlled X;
- uniformly controlled one-qubit gates at arbitrary workspace;
- coherent CNOT copy–uncopy;
- the optimal state-preparation comparison frontier.

The published article corresponds to `arXiv:2202.11302v2`; the imported
statements were also checked in v3.

The earlier paper by Sun, Tian, Yang, Yuan, and Zhang established the preceding
time–space landscape and is retained as the historical predecessor and original
source of selected primitives.  For the present theorem, the later uniform
framework is sufficient for all workspace regimes.

The connection to the Hopf problem is constructive rather than automatic.  Once
the Hopf completion is written as addressed complete-operator layers, the
state-preparation primitives can be adapted while preserving its designated
columns.

## 2. Uniformly controlled gates and the Möttönen route

Möttönen and coauthors introduced the uniformly controlled rotation language
for state preparation.  Bergholm and coauthors developed the corresponding
general uniformly controlled one-qubit gates.

The earlier Hopf-QBP compiler follows this route: compute the shared lower-
suffix-zero predicate into one clean flag, then apply a smaller multiplexed
rotation selected by the prefix and flag.  This already shows that the complete
Hopf frame can be compiled with linear asymptotic gate count and constant
additional workspace.

At strict zero workspace, inserting identity blocks into one full-width
multiplexor is not enough.  Although the logical block table is sparse, the
standard Walsh/Gray-code angle transform generically repeats nonzero
coefficients across all suffix frequencies.  The physical angle table is not
sparse in the required sense.

The strict-zero schedule instead retains a narrow UCG.  One original suffix data
qubit carries the predicate temporarily, and two half-angle width-$`d+2`$ UCGs
produce the full addressed depth while restoring that qubit exactly.

## 3. Borrowed workspace and controlled-unitary roots

Several established techniques are close to the strict-zero circuit.

- Barenco and coauthors use roots and conjugation identities in controlled-
  unitary decompositions.
- Later borrowed-qubit constructions use unknown-state logical wires while
  restoring them exactly.
- Conditionally clean ancillas and toggle-detection patterns organize
  cancellation across different original values of a borrowed wire.
- Ancilla-free multi-controlled $\mathrm{SU}(2)$ constructions give related individual
  controlled-rotation primitives.

For one Hopf prefix, put

```math
C_p=R_y(\theta_p/2).
```

The identities

```math
C_p^2=R_y(\theta_p),
\qquad
XC_pX=C_p^{-1}
```

make the original zero-suffix sector square to the required rotation and make
the unwanted original value of the borrowed bit cancel.

The borrowed bit is not an extra dirty ancillary wire.  It is part of the
logical input, the desired operation depends on its original value, and the
proof must restore it on arbitrary entangled inputs.

### Narrow strict-zero contribution

The repository does not claim the invention of borrowed qubits, conditional
cleanliness, toggle detection, controlled-unitary roots, Pauli conjugation, or
UCGs.

The Hopf-specific statement is:

> One original suffix data qubit serves as a restored in-place predicate
> carrier, reducing every addressed Hopf depth to two total-width-$`d+2`$ UCGs
> and linear predicate toggles.  Summed over the tree, this gives the optimal
> strict-zero complete-frame frontier.

## 4. Sparse and restricted multiplexors

[Xu et al., “A Unified Framework for Optimizing Uniformly Controlled
Structures in Quantum Circuits”](https://arxiv.org/abs/2512.08675) study
restricted UCGs and synthesis savings from limited control participation.
The addressed Hopf layer is related, but its sparsity has a particular form:
one long all-zero suffix predicate is shared by every prefix-selected rotation.

The borrowed-suffix factorization converts this logical predicate into smaller
physical UCG width rather than attempting to prune the full-width Möttönen angle
table.  That is the structural distinction relevant to the strict-zero bound.

## 5. Coherent routing and workspace parallelism

For large workspace, the Hopf tail has the exact direct-sum form

```math
R_t^{(n)}
=\bigoplus_{r=0}^{2^t-1}W_s^{(r)}.
```

The compiler realizes this form by standard reversible ingredients:

- balanced CNOT fanout of the prefix bits;
- Fredkin routing of one suffix-token block;
- disjoint token-controlled subtree frames;
- exact uncomputation and inverse routing.

The state-preparation framework supplies the coherent-copy and UCG primitives.
The Hopf-specific part is the tree-cut identity, the register layout, the
route–operate–unroute schedule, and the shared workspace envelope for prefix and
tail.

The implementation gives the router explicitly and tests arbitrary complex
inputs in which the prefix and suffix are entangled.

## 6. Hopf coordinates and inverse-frame gradients

The first Hopf work develops the balanced binary chart, inverse map, diagonal
metric, and coordinate directions.  The second develops the addressed
state-and-marker frame, global magnitude record, direct phase record, and
checkpoint interface.

The compiler target is

```math
W_{\mathbb R}|0^n\rangle=|\psi\rangle,
\qquad
W_{\mathbb R}|\lambda(j)\rangle=|e_j\rangle.
```

For unrestricted angles,

```math
\partial_{\theta_j}|\psi\rangle
=a_j|e_j\rangle,
\qquad
g_{j,j}=a_j^2.
```

On the canonical domains, $a_j$ is nonnegative.  At a singular coordinate, the
raw derivative vanishes while the complete parameter tuple still selects an
orthogonal marker-frame continuation.

The phase-dressed complex magnitude frame is

```math
W_{\mathbb C,\mathrm{mag}}
=D_{\mathrm{ph}}W_{\mathbb R}.
```

The leaf-phase coordinates remain a separate direct stream.

The synthesis theorem proves that these prescribed magnitude-frame operators,
not only their first columns, attain the all-workspace state-preparation
frontier.

## 7. Comparison by mathematical role

| Line of work | Object or result | Role used here | Additional Hopf structure |
|---|---|---|---|
| Möttönen; Bergholm | multiplexed rotations and general UCGs | block-diagonal one-qubit-gate language | addressed depths are reduced to narrow UCGs rather than full-width sparse tables |
| exact QSP time–space tradeoffs | arbitrary state preparation with clean workspace | target frontier and compiler primitives | designated frame columns must survive the adaptation |
| controlled-unitary roots and borrowed workspace | exact conditioned operations with limited ancillas | cancellation and restored logical carriers | the predicate includes the borrowed suffix bit's original value |
| restricted UCGs | repeated or sparse block structure | neighboring structural context | one shared zero-suffix predicate is factored through a restored logical bit |
| reversible routing | CNOT fanout, Fredkin networks, uncomputation | coherent branch movement and parallel action | realizes the Hopf tail direct sum within the prefix workspace envelope |
| Hopf chart | state map, inverse, metric, coordinate directions | specifies the structured columns | yields the addressed tree layers |
| Hopf QBP | inverse-frame magnitude record and direct phase record | makes the completion operational | requires frame-safe rather than state-column compilation |

## 8. Contribution at theorem level

The exact compiler contribution is the following chain.

1. Identify complete frame safety as the compiler contract used by the global
   inverse-frame record.
2. Give an exact state-column counterexample.
3. Express the frame as zero-suffix-addressed tree layers.
4. Close the zero-workspace endpoint with the borrowed-suffix echo.
5. Build a clean binary–one-hot prefix decoder and explicit coherent router.
6. Match the state-preparation size–depth frontier for every $m\geq0$.
7. Transfer the result to the phase-dressed complex magnitude frame.
8. State the QBP consequence under matched program and output-task conventions.

The theorem is specific to the structured Hopf completion.  It does not imply
that an arbitrary family of prescribed unitary columns can be implemented at
state-preparation cost.

## 9. Fault-tolerant lineage and the precision register

Let $N=2^n$ and let $L$ denote the requested number of accuracy bits. Three
established results are particularly close to the fault-tolerant construction.

**Gosset–Kothari–Wu (GKW).** Their Theorems 1.1–1.2 give the optimal unrestricted
T-count $\Theta(\sqrt{NL}+L)$ for state preparation and diagonal synthesis.
Their Appendix B already treats complete single-qubit multiplexors and allocates
error geometrically across state-preparation layers, obtaining
$O(\sqrt{NL}+nL)$. These are imported benchmarks and methods. A state-preparation
theorem alone does not specify the remaining columns of the Hopf frame, and an
unrestricted-ancilla theorem does not supply a prescribed clean/dirty budget.
Their Theorem 4.2 supplies the ancilla-independent lower bound used through the
frame's diagonal subfamily.

**Low–Kliuchnikov–Schaeffer (LKS).** Section 2, Table 2, Figure 1(d), and Eq. (8)
provide dirty-assisted SelectSwap lookup with exact bank restoration. For an
$m$-bit table with $D$ entries, its T count is $O(D/\lambda+m\lambda)$ with
$m\lambda$ dirty bank qubits and $O(\log D+m)$ clean wires. Appendix C also
gives exact Boolean XOR constructions using dirty auxiliaries. The basis-state
identity in Eq. (8) extends to reference-entangled inputs by linearity. We do
not claim dirty lookup or its cancellation identity as a contribution. Its
dirty bank does not replace the clean output word needed by an instruction-word
interpreter. Section 5 supplies the finite-width counting method.

Appendix C also separates indicator parallelism from count-efficient bank
selection. The [parallel lookup proof](../PARALLEL_DIRTY_LOOKUP.md) gives a
literal exact realization by conjugating one X with the inherited bank
router, with scratch-free indicator work and linear address T-depth.
The predicate schedule uses Khattar–Gidney
[Section 5.4](https://arxiv.org/html/2407.17966v1#S5.SS4), allocating its two
dirty helpers only between completed queries. Its role here is to
compose the inherited lookup idea with complete two-clean real frames,
retaining their T-count while reducing T-depth. It does not claim a new
general lookup tradeoff or equally small Clifford depth.

The [partial-batch schedule](../../research/depth/BATCHED_DIRTY_LOOKUP.md) reuses these
indicator and linear-echo mechanisms with a smaller simultaneous pool.
LKS Appendix C, Theorem 2 states a count-efficient tradeoff with
$`O(\lambda\sqrt{Qm})`$ dirty workspace, $`\lambda\ge1`$, so its
displayed range begins at square-root-scale width. Partial batches
retain the count-efficient bank choice below that scale and supply the
intermediate-workspace T-depth ledger and full return contract for Hopf
composition. The underlying lookup mechanism is inherited; no priority
claim over general lookup tradeoffs is made. Motlagh–Pocrnic,
[arXiv:2605.20334v1, Sections II.2–II.4](https://arxiv.org/html/2605.20334v1),
improve Toffoli-count prefactors through SelectCopy and sequential
output-bit packets. Their displayed constructions retain initialized
outputs and do not state the simultaneous $`B^{-2}`$ T-depth dependence
audited here. The batch guards reserve F8's two dirty helpers during
the live query, separately from its occupied banks and indicators.

The [amortized refinement](../AMORTIZED_DIRTY_LOOKUP.md) combines the same
dirty echo with a multiplexed family of linear shears. The dirty tree
traversal and cancellation are established techniques; see
[Khattar–Gidney, arXiv:2407.17966v2, Sections 4, 7.1 and 7.4](https://arxiv.org/html/2407.17966v2).
[Boyd, arXiv:2312.00696v2, Section III and Appendix A](https://arxiv.org/html/2312.00696v2)
already groups commuting SELECT/QROM operators and changes Clifford
bases to expose parallel action, using clean address-fanout registers.
[Kim–Laakkonen, arXiv:2512.24982v1, Theorems 3 and 5](https://arxiv.org/pdf/2512.24982v1)
give constant-depth control of CNOT and Clifford circuits without
ancillas or measurements. The local rank-reduction lemma specializes
that established possibility with explicit rank-sensitive T-count and
Clifford-count bounds. The new composition moves the low-address
indicator outside the high-address traversal and supplies its all-dirty
return identity and complete Hopf resource ledger. It removes the batch
logarithm and proves matching count and depth in an explicit
precision/workspace range, including fixed and inverse-polynomial
accuracy in N. The variable-accuracy extension uses the same native
query and the existing layer error bounds.
No priority is claimed for dirty iteration, commuting-operator grouping,
or constant-depth controlled Clifford circuits.

**Bausch.** Equations (4) and (6) of *Fast Black-Box Quantum State Preparation*
already use a geometric precision register and a bit oracle addressed by the
data index and precision position. Section 2.3.3 includes a capped geometric
source. Its explicit exact preparation uses a linear-size unary temporary,
then binary conversion and routing cleanup. Consequently, logarithmic output
width there does not by itself establish logarithmic peak clean workspace.
The compact Gray construction used here supplies the specific simultaneous
resource guarantee: exact $O(L)$ T count, $O(\log L)$ peak clean workspace,
ordinary binary labels, and restored temporary work. The geometric weighting
and bit-oracle idea are credited to Bausch; no strict gate-count improvement
over his routing variants is asserted.

The addressed SU(2) compiler combines this source with a constant-size Pauli
linear combination and oblivious amplitude amplification. It is a
space-efficient application of these ingredients, with all preparations,
adjoints, reflections, and one-bit tables charged. It is not a new general
principle of dyadic sampling or amplitude amplification.

## 10. What the full-frame T theorem adds

The main additional construction is the uniform-in-precision composition of
complete frame corrections. A shared geometric source supplies approximate
shift amplitudes on accepted branches; the source is retained between stages.
The residual decomposition, failure tracking, and one final amplification
control the full initialized isometry while returning dirty workspace exactly.
The failure counter uses the established block-product compression gadget of
[Low–Wiebe, Lemma 13](https://arxiv.org/html/1805.00675v2), in the form explained
by [Fang–Lin–Tong, Section 2.4, Lemma 3 and Appendix D](https://arxiv.org/html/2208.06941v2).
It coherently retains failed histories while reusing block work. Neither
logarithmic history storage nor final amplification is a new general
mechanism here; the additional ingredients are the residual representation,
shared precision source, and charged frame-specific schedule.
Under the conditions in the [theorem](../FAULT_TOLERANT_COMPILER.md), this gives

```math
T=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{n+a+b}\right),
```

where $a$ and $b$ count clean and dirty ancillary qubits. The theorem retains
its sufficient clean-workspace reservation. A separate retained corollary for
all clean budgets, under an additional restriction on precision and dirty
width, is proved in the [borrowed-workspace appendix](../BORROWED_WORKSPACE_COMPILER.md).
The lower bounds follow from GKW and LKS with the stated reductions; the
contribution is the upper construction and the resulting match in those
regimes. The all-clean-budget extension also uses the earlier borrowed-work
compiler; it is a compiler splice, not an independent lower-bound technique.

The simpler direct sampler already costs
$O(\sqrt{NL}+nL+NL/(n+a+b))$ with logarithmic precision workspace. For example,
at $L=n^2$, $a=\Theta(n)$ with sufficient headroom, and
$b=\Theta(n\sqrt N)$, it already attains $\Theta(n\sqrt N)$ T gates.
This illustrates the small-clean sampler; it is not an additional improvement
due to shared-source composition. Eliminating repeated precision charges is
the reason for the stronger uniform-precision construction.

Tan's published Theorem I.1 and Lemma IV.1 treat general unitary synthesis and
batched controlled gates. Remark IV.2 uses zeroed precision-length instruction
registers; it does not directly give the clean/dirty contract above. The
Yuan–Zhang depth theorem uses arbitrary continuous one-qubit gates, so it is not
a T-count comparison. Recent sparse-QROM and sparse-state bounds address
different input families. These distinctions delimit the use of each result;
they do not establish priority by absence of a matching theorem.

At a fixed parameter value, a full-isometry error can bound the bias of the
[bounded-score QBP estimator](../QBP_APPROXIMATION.md).
This does not differentiate the compiled Clifford+T word as a function of the
parameters and does not prove optimal gradient-query complexity.

## 11. Operator sources and the one-clean refinement

The [operator-source compiler](../OPERATOR_SOURCE_COMPILER.md)
uses two initialized flags and an arbitrary dirty bank. Full-space
anticommuting loaders are established by
[Kerenidis–Prakash](https://arxiv.org/html/2202.00054v2), Definitions 4.4/4.6
and Theorem 4.9; their scalar-overlap product algebra also appears in
[Chee et al.](https://arxiv.org/pdf/2301.07477), Appendix C. Dirty SelectSwap
is credited to LKS, and normalization-two amplification to
[Berry et al.](https://arxiv.org/pdf/1412.4687), Eqs. (11)–(15).

The local construction specializes the loader to native geometric weights,
programs its signs through whole-word dirty queries, and combines a scalar
block, suffix echo, and full-output amplification into the prescribed Hopf
frame. With $`b\ge L+n+7`$ it gives $`T=O(N+nL)`$; with
$`b\ge2(L+n+7)`$ its banked form gives
$`T=O(\sqrt{NL}+nL+NL/b)`$. Both use $`O(NL)`$ Clifford gates.
The operator core returns approximately within the full-isometry bound;
lookup banks and selectors return exactly. Literal diagonal compilation
also gives the phase-dressed magnitude-frame corollary with its own error
allocation and dirty-width threshold.

The [conditional-suffix compiler](../CONDITIONAL_SUFFIX_COMPILER.md) gives
$`O(N+L\ell_*(n))`$ T gates and $`O(NL)`$ Clifford gates for real
frames at $`n\ge1`$, $`L\ge6`$, $`a=2`$, and $`b\ge L+n+7`$.
Here $`\ell_*(n)=1+\log_2^*(n+2)`$, with repeated base-two logarithms
stopping at a value at most one. Star residuals and streamed coarse programs
use logical suffix qubits as clean work only in the active sector. At
$`L=N`$, the upper bound is $`O(N\ell_*(n))`$; the $`\Omega(N)`$
lower bound remains unchanged.
With $`b\ge2(L+n+7)`$, the grouped banked bound is
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, still with $`O(NL)`$ Clifford gates.
These are explicitly charged constructions in their stated regimes;
loaders, dirty lookup, overlap algebra, and amplification retain their
earlier attribution.

Hierarchical block organization and in-place index updates are also
established. [Nguyen–Kiani–Lloyd](https://arxiv.org/pdf/2201.11329),
Section 4.1 and Appendix C.1, use a separate initialized local-vector
register; their low-rank construction in Section 4.3 uses factor-state
oracles. [Setty](https://arxiv.org/pdf/2508.21667), Theorem 1 and Section 2,
includes in-place logical shift/delete operations while retaining initialized
data-label/PREP work. The local ancestor blocks instead use the logical word
as the coefficient address, a dirty scalar precision source, and explicitly
charged streaming lookup. The contribution is this complete-workspace
realization for the prescribed residual structure, not a new generic LCU
principle or an assumed QRAM interface.

### One initialized flag

The [one-clean compiler](../ONE_CLEAN_COMPILER.md) retains the established
Clifford-loader and dirty-lookup ingredients. Its additional construction
conjugates one scalar-source block by another, routing their rejected
components through different logical Pauli operators. The Clifford
anticommutator reduces the accepted block to a real rotation with two
programmable scalar coefficients. A head coefficient and two geometric
tails supply those coefficients. The standard paired-Majorana algebra and
its rotations are recalled in
[Bravyi, Section II, Eqs. (2)–(5)](https://arxiv.org/pdf/quant-ph/0404180).
The local source packs both tails into essentially one dirty qubit per
accuracy bit. The
contribution is this native source and conjugated block word, with its
complete-input error and workspace accounting, not Majorana operators,
Clifford-algebra overlap identities, or source conjugation in general.

At the fixed amplitude $`\beta=\sin(\pi/10)`$, two amplification
iterations use five forward or inverse block calls. The amplification
mechanism is inherited from
[Brassard–Høyer–Mosca–Tapp, Section 2, Eqs. (7)–(8)](https://arxiv.org/pdf/quant-ph/0005055)
and the oblivious block setting above. The local proof verifies its literal
phase and complete-isometry error for the programmed source, including dirty
core return; it does not claim a new amplitude-amplification principle.

The real-frame baseline and grouped bounds therefore require only
$`a\ge1`$, with the same $`b\ge L+n+7`$ dirty reservation. Their
banked refinements retain the doubled dirty threshold. A flag-only Pauli
routing also gives one-clean literal diagonals and, through the retained
certified Euler decomposition, complete one-target U(2) multiplexors.
Sequential composition with the literal diagonal also gives the grouped
phase-dressed complex magnitude-frame corollary, at $`b\ge L+n+8`$ or
its doubled banked threshold. This is closure under the existing
complete-isometry contract, with no new phase-gradient columns.

The real-rotation word also commutes with X on its amplification signal.
The [full-space symmetry argument](../ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
therefore permits borrowing that signal, giving the zero-clean layerwise
bound $`T=O(N+nL)`$, $`G=O(NL)`$ at $`b\ge L+n+7`$.
The scalar-phase word does not share this symmetry, and the grouped proof
retains its external clean predicate. These corollaries do not lower the
proved full-frame endpoint T count or extend the separately proved
two-clean T-depth schedules.

The same primitive prices the
[weighted forward component](../../research/endpoint/WEIGHTED_TRANSPORT_BLOCK.md) at
$`O(N+nL)`$ T gates using only borrowed synthesis work alongside its
logical dilation signal. Approximate return of the core and borrowed
amplification signal is included in the full-operator error. The earlier
$`O(L\sqrt N)`$ route remains useful for its exact dirty-work return.
Neither component bound removes the full-frame endpoint gap.

The [complete residual assembly](../../research/endpoint/RESIDUAL_ASSEMBLY.md) incorporates the
diagonal and forward term into one affine tree map, then selects it or the
reverse dilation using two clean flags. The scalar elimination uses
standard Schur complements (F18 in the source map); the local contribution is the affine
recursion, complete tree-mode allocation, and native controlled-branch
accounting. The prepared two-term linear combination and normalization-two
amplification are inherited from
[Berry et al., Eqs. (7)–(15)](https://arxiv.org/pdf/1412.4687).
Borrowed synthesis work and full-operator error estimates make the selected
branches usable while both flags are occupied. This gives
$`T=O(N+nL)`$ and $`G=O(NL)`$ at $`a=2`$,
$`b\ge L+n+7`$; it is not a new generic LCU method or an improvement
over the grouped frame bound. The constant-call reduction applies to the
new affine and reverse dilations, not to an uncontrolled opaque forward
block alone.

The [source map](../SOURCE_MAP.md) gives exact theorem numbers and local consumers.
The comparisons identify dependencies and specific additional constructions;
they do not certify priority.

## 12. Contemporary comparisons and the broader compiler contribution

The following comparisons were checked against primary sources on
**22 September 2026**, with the F12 comparison updated to v3 on
**9 October 2026**. Precision, initialization, and error contracts
matter as much as T-count exponents.

**General unitary synthesis.** Tan is no longer the newest generic upper
bound. [Yuan–Zhang–Zi, arXiv:2608.17846v2, Theorem I.1 and Definition II.1](https://arxiv.org/html/2608.17846v2)
give, for $`n+L\le N`$,

```math
T=O\!\left(nN^{5/4}(n+L)^{5/8}\right),
\qquad a=O\!\left(N\sqrt{n+L}\right).
```

Their ancillary qubits are clean. Their error criterion is the complete
initialized-isometry norm, including leakage and literal phase, so that
criterion itself is not a novelty claim here. They treat arbitrary dense
unitaries; the Hopf family has only $`N-1`$ real parameters and has additional
tree structure. The comparison establishes the importance of using that
structure, not an improvement to arbitrary-unitary synthesis.
[Fang–Heunen–Wang, arXiv:2607.12907v1, Theorem 1.2 and Corollary 3.8](https://arxiv.org/html/2607.12907v1)
give an instance-dependent bound in terms of phase-optimized Frobenius
distance to the Clifford group. Their circuit uses initialized ancillas and
their detailed conclusion is in diamond distance. It supplies neither a
two-clean guarantee nor the prescribed frame's worst-case frontier.

**The closest precision and workspace comparison.**
[Yamazaki–Akibue–Sano, arXiv:2603.14202v3, Theorem 1 and Section 3.1](https://arxiv.org/html/2603.14202v3)
give the deterministic bound $`m+O(\sqrt{N(m+1)})`$ for every supplied
SU(2) block table at diamond distance $`\varepsilon`$. Here m is the
largest minimum even T count of an ancilla-free block approximation.
The circuit uses $`O(m+n+1)`$ clean and $`O(\sqrt{N(m+1)})`$ dirty
ancillas, exactly restored; its clean instruction word has
$`m+K+O(1)`$ bits with $`K\le\min\{N,m/2+1\}`$.
The leading precision coefficient three and its matching lower bound
are probabilistic statements for independent Haar blocks, with the
matching limit taken at fixed n. Equation (27) aligns determinant-one
block lifts and their relative signs. The decisive endpoint mismatch is
the $`O(m+n+1)`$ clean-work allocation.

GKW's diagonal proof likewise computes a precision-length instruction word
into an initialized register, then applies that word and uncomputes it. LKS
can replace its lookup scratch with dirty banks, but an arbitrary dirty
output word is not a known instruction word. The operator-source constructions
avoid that instruction register by using full-space operator identities on
the dirty core; the one-clean refinement changes the scalar-block word.

This distinction already gives a result beyond the Hopf family. For **every
literal diagonal** on $`n`$ qubits, the
[one-clean proof](../ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors) establishes

```math
a=1,\quad b\ge2(L+n+7),\qquad
T=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O(NL).
```

The upper bound includes approximate joint return of the dirty operator
core; lookup banks and selectors return exactly. Since $`n+1+b=\Theta(b)`$
in this range, GKW's diagonal lower bound and LKS-style finite-width counting
give the stated worst-case match. At $`L=N`$ and $`b=2(N+n+7)`$ this is
$`\Theta(N)`$ T count with one clean qubit. This diagonal endpoint is
settled by the retained proof, whereas the corresponding complete-frame
endpoint still has a gap between $`\Omega(N)`$ and
$`O(N\ell_*(n))`$. The broader scientific contribution
is therefore a precision/workspace compiler for a standard operator family,
together with the structured complete-frame extensions. It is not a claim
that the entire constant-clean frame frontier is matched.

The same refinement covers **arbitrary complete one-target U(2)
multiplexors**, with $`N=2^n`$ blocks, one initialized flag, and
$`b\ge2(L+n+9)`$ dirty qubits, at
$`\Theta(\sqrt{NL}+L+NL/b)`$ T gates. Four addressed Euler
factors suffice; their implementation reuses the flag and dirty core with
all return error included. The [one-clean corollary](../ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
uses the retained certified matrix-entry procedure through finite approximate
Euler search,
so no nonsingular-chart or exact-zero promise is introduced. This directly
addresses the complete multiplexor task in the GKW and Yamazaki–Akibue–Sano
comparisons, while supplying the one-clean guarantee. It does not
improve their leading constants, and its general classical coordinate search
is not claimed efficient. It also does not combine all Hopf tree depths into
one jointly charged precision source.

**Why general compression does not remove the distinction.**
[Vasconcelos–Gilyén, arXiv:2507.07900v2](https://arxiv.org/html/2507.07900v2)
give an ancilla-uncomputation procedure that still queries the original
block encoding on its initialized work. Thus it reduces retained work
between calls, not the initial peak-clean requirement to a constant number
of qubits. Their
exact-product logarithmic ancilla lower bound is for the specified multiple
coherent measurement circuit class. Their approximate compression theorem
requires near-identity supplied dilations, not merely near-identity accepted
blocks. Neither statement proves or resolves our constant-clean frame gap.
[Ma–Joven–Liu, arXiv:2609.11153v1, Theorem 3.1](https://arxiv.org/html/2609.11153v1)
compress a block encoding to at most $`n+2T`$ clean ancillas without raising
T count, possibly changing normalization. This supports counting lower
bounds; it is not a constant-clean unitary implementation theorem.

**Lookup constants.**
[Motlagh–Pocrnic, arXiv:2605.20334v1, Section II](https://arxiv.org/html/2605.20334v1)
improve dirty-QROM constants using SelectCopy and shared work across bit
packets. The output instruction register remains initialized in their
construction. These improvements preserve the asymptotic lookup tradeoff
used here and do not supply the one-clean diagonal compiler. We claim no
improvement to their constant factors.

## References highlighted here

- A. Barenco et al., “Elementary gates for quantum computation,” *Physical
  Review A* **52**, 3457–3467 (1995).
- M. Möttönen et al., “Transformation of quantum states using uniformly
  controlled rotations,” *Quantum Information and Computation* **5**, 467–473
  (2005).
- V. Bergholm et al., “Quantum circuits with uniformly controlled one-qubit
  gates,” *Physical Review A* **71**, 052330 (2005).
- X. Sun, G. Tian, S. Yang, P. Yuan, and S. Zhang, “Asymptotically Optimal
  Circuit Depth for Quantum State Preparation and General Unitary Synthesis,”
  *IEEE TCAD* **42**, 3301–3314 (2023).
- P. Yuan and S. Zhang, “Optimal (controlled) quantum state preparation and
  improved unitary synthesis by quantum circuits with any number of ancillary
  qubits,” *Quantum* **7**, 956 (2023).
- B. Claudon et al., “Polylogarithmic-depth controlled-NOT gates without
  ancilla qubits,” *Nature Communications* **15**, 5886 (2024).
- T. Khattar and C. Gidney, “Rise of conditionally clean ancillae for efficient
  quantum circuit constructions,” *Quantum* **9**, 1752 (2025).
- B. Zindorf and S. Bose, “Efficient implementation of multi-controlled quantum
  gates,” *Physical Review Applied* **24**, 044030 (2025).
- C. Xu et al., [“A Unified Framework for Optimizing Uniformly Controlled
  Structures in Quantum Circuits”](https://arxiv.org/abs/2512.08675),
  arXiv:2512.08675 (2025).
- J. Bausch, [“Fast Black-Box Quantum State Preparation”](https://arxiv.org/pdf/2009.10709v4),
  arXiv:2009.10709v4 (2022), Eqs. (4), (6), and Section 2.3.3.
- G. H. Low, V. Kliuchnikov, and L. Schaeffer,
  [“Trading T gates for dirty qubits in state preparation and unitary synthesis”](https://arxiv.org/html/1812.00954v2),
  *Quantum* **8**, 1375 (2024).
- G. Brassard, P. Høyer, M. Mosca, and A. Tapp,
  [“Quantum Amplitude Amplification and Estimation”](https://arxiv.org/pdf/quant-ph/0005055),
  arXiv:quant-ph/0005055, Section 2, Eqs. (7)–(8).
- D. W. Berry, A. M. Childs, R. Cleve, R. Kothari, and R. D. Somma,
  [“Simulating Hamiltonian dynamics with a truncated Taylor series”](https://arxiv.org/pdf/1412.4687),
  *Physical Review Letters* **114**, 090502 (2015), Eqs. (11)–(15).
- D. Gosset, R. Kothari, and K. Wu,
  [“Quantum state preparation with optimal T-count”](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/),
  *Quantum* **10**, 2168 (2026).
- X. Tan, [“Unitary Synthesis with Fewer T Gates”](https://journals.aps.org/prxquantum/pdf/10.1103/pxhd-9s9q),
  *PRX Quantum* **7**, 033058 (2026).
- T. Li, F. Ou, X. Wang, P. Yao, P. Yuan, and S. Zhang,
  [“Optimal T Counts under Sparsity: from QROM to State Preparation and Block Encoding”](https://arxiv.org/html/2607.28260v1),
  arXiv:2607.28260v1 (2026).
