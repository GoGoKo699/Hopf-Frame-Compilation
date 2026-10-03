# Source and dependency map

[← Verification](VERIFICATION.md) · [Complete narrative](../REVIEW.md) · [Related work →](RELATED_WORK.md)

This page separates the proof into three layers:

1. an exact state-preparation compiler toolkit imported from the literature;
2. a compact Hopf frame and QBP interface inherited from the earlier Hopf work;
3. complete-operator factorizations, schedules, and bounds proved in this
   repository.

The tables identify the precise fact used, rather than citing an entire paper as
a premise. Section 5 adds the separate Clifford+T toolkit. Yuan–Zhang remains
the normative framework for the exact arbitrary-one-qubit+CNOT theorem; the
fault-tolerant theorem has its own precision and workspace hypotheses.

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
| C1 | Theorem 2: exact arbitrary state preparation has $\Theta(2^n)$ size and $`\Theta\!\left(n+\frac{2^n}{n+m}\right)`$ depth for every clean-workspace budget | comparison frontier and target lower-bound scale | [compiler theorem](COMPILER_THEOREM.md), [`unified_compiler.py`](../compiler_robust_hopf/unified_compiler.py) |
| C2 | Lemma 5: exact ancilla-free multi-controlled X has linear size and depth | strict-zero toggles, direct suffix flags, and branch predicates | [`strict_zero_echo.py`](../compiler_robust_hopf/strict_zero_echo.py), [`router.py`](../compiler_robust_hopf/router.py) |
| C3 | Lemma 6: a total-width-$`q`$ UCG with $w$ clean work qubits has $O(2^q)$ size and $`O\!\left(q+\frac{2^q}{q+w}\right)`$ depth | half-angle UCGs, flagged layers, subtree frames, and the phase diagonal | [`strict_zero_echo.py`](../compiler_robust_hopf/strict_zero_echo.py), [`unified_compiler.py`](../compiler_robust_hopf/unified_compiler.py) |
| C4 | Lemma 9: coherent CNOT-tree copying and exact uncopying | control fanout in the conditioned prefix and coherent router | [`tree_decoder.py`](../compiler_robust_hopf/tree_decoder.py), [`router.py`](../compiler_robust_hopf/router.py) |
| C5 | Theorem 1: generic controlled state preparation depends on combined index and target width | comparison showing why a generic all-column construction would have quadratic Hilbert-space scale | [`resource_bounds.py`](../compiler_robust_hopf/resource_bounds.py), [related work](RELATED_WORK.md) |

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
| H1 | a balanced complete binary tree with $N=2^n$ leaves gives $N-1$ real magnitude coordinates | chart definition | [Hopf interface](HOPF_INTERFACE.md), [`frames.py`](../compiler_robust_hopf/frames.py) |
| H2 | the sine–cosine path map covers normalized real states; leaf phases give the complex chart | forward state map | [`real_tree_data`](../compiler_robust_hopf/frames.py), [`complex_state`](../compiler_robust_hopf/complex_analysis.py) |
| H3 | real nonfinal angles use $[0,\pi/2]$, the final real depth uses $[0,2\pi)$, and complex magnitudes use $[0,\pi/2]$ | canonical chart domains | [`canonical_magnitude_angle_mask`](../compiler_robust_hopf/frames.py), [`test_frames.py`](../tests/test_frames.py) |
| H4 | the oriented incoming amplitude obeys $\partial_{\theta_j}\lvert\psi\rangle=a_j\lvert e_j\rangle$ and $g_{j,j}=a_j^2$ | differential geometry | [`RealTreeData.incoming_amplitude`](../compiler_robust_hopf/frames.py), [`test_complex_analysis.py`](../tests/test_complex_analysis.py) |
| H5 | on the canonical domains $a_j\geq0$, so $a_j=\sqrt{g_{j,j}}$ | canonical-domain consequence | [`test_frames.py`](../tests/test_frames.py) |
| H6 | if $g_{j,j}=0$, the raw differential vanishes while the parameter tuple still selects a unit marker-frame continuation | singular-coordinate boundary | [`regular_coordinate_mask`](../compiler_robust_hopf/frames.py), singular tests in [`test_frames.py`](../tests/test_frames.py) |
| H7 | the marker $\lambda(j)$ places the state and coordinate-frame directions in one known computational basis | differential-frame interface | [`conventions.py`](../compiler_robust_hopf/conventions.py), [`frames.py`](../compiler_robust_hopf/frames.py) |
| H8 | the complete real frame is a product of prefix-selected rotations restricted to the zero-suffix sector | addressed-frame structure | [Hopf interface](HOPF_INTERFACE.md), [`direct_addressed_depth_layer`](../compiler_robust_hopf/frames.py) |
| H9 | the phase-dressed complex magnitude frame is $W_{\mathbb C,\mathrm{mag}}=D_{\mathrm{ph}}W_{\mathbb R}$ | complex magnitude interface | [`complex_magnitude_frame_matrix`](../compiler_robust_hopf/frames.py) |

The complete inverse map and optimization context are contained in the first
Hopf paper and repository.  The addressed frame, global magnitude record,
direct phase record, checkpoint interface, and statistical task boundaries are
contained in `Hopf-QBP`.

## 3. Inherited QBP interface

The tracked `Hopf-QBP/main` baseline is

```text
a9885317cf998a7df87ca07ba86e3bd4f0f419ef
```

and is reconciled in [`SYNC.md`](../SYNC.md) and
[`provenance/upstream.json`](../provenance/upstream.json).

| ID | QBP fact used here | Role in the compiler consequence | Local counterpart |
|---|---|---|---|
| Q1 | the global magnitude circuit applies the complete inverse differential frame | fixes the required compiler contract | [frame-safe compilation](FRAME_SAFE_COMPILATION.md) |
| Q2 | one X-basis outcome contributes a parity record to every magnitude coordinate | shared-record mechanism | [`decoders.py`](../compiler_robust_hopf/decoders.py), [`test_decoders.py`](../tests/test_decoders.py) |
| Q3 | leaf-phase derivatives use a separate signed one-hot stream | separates the complex magnitude frame from the phase record | [`decoders.py`](../compiler_robust_hopf/decoders.py), [`complex_analysis.py`](../compiler_robust_hopf/complex_analysis.py) |
| Q4 | the primary finite-shot target is simultaneous absolute accuracy of the raw coordinate gradient | defines the $O(\log n)$ execution statement | [QBP consequence, Section 5](QBP_CONSEQUENCE.md#5-statistical-target) |
| Q5 | complete-vector, relative, normalized-frame, and natural-gradient targets have different conditioning | prevents overextension of the raw-coordinate claim | [QBP consequence, Section 5](QBP_CONSEQUENCE.md#5-statistical-target) |
| Q6 | the validated core assumes phase-calibrated controlled access to a Hermitian unitary observable | access-model premise | [QBP consequence, Section 2](QBP_CONSEQUENCE.md#2-controlled-observable-interface) |
| Q7 | a checkpoint compiler needs equality on its active interface, not merely one state column | separate checkpoint boundary | [compiler boundaries](COMPILER_BOUNDARIES.md) |

## 4. Results established in this repository

| ID | Result | Analytic proof | Implementation and finite evidence |
|---|---|---|---|
| R1 | complete frame safety implies substitution and, at regular points with exact prepared state, is necessary up to common phase for all observable-dependent fixed-decoder means | [frame-safe compilation](FRAME_SAFE_COMPILATION.md) | actual probability-decoder tests with complex phases, workspace leakage, and singular markers |
| R2 | one correct state column need not preserve the decoded gradient | [two-qubit argument](../REVIEW.md#12-a-complete-two-qubit-obstruction) | [`test_compiler_boundaries.py`](../tests/test_compiler_boundaries.py) |
| R3 | the borrowed-suffix echo implements one nonfinal addressed layer with no ancillary wire | [compiler theorem, Lemma Z](COMPILER_THEOREM.md#lemma-z-borrowed-suffix-echo) | [`strict_zero_echo.py`](../compiler_robust_hopf/strict_zero_echo.py), [`test_strict_zero_echo.py`](../tests/test_strict_zero_echo.py) |
| R4 | the strict-zero frame has $\Theta(N)$ size and $\Theta(n+N/n)$ depth | [compiler theorem, Proposition Z](COMPILER_THEOREM.md#proposition-z-strict-zero-resources) | [`strict_zero_audit.py`](../compiler_robust_hopf/strict_zero_audit.py), exact-rational tests |
| R5 | the prefix is a conditioned smaller frame and the tail is a direct sum of subtree frames | [compiler theorem, Lemmas T1–T2](COMPILER_THEOREM.md#6-exact-tree-cut) | [`tree_structure.py`](../compiler_robust_hopf/tree_structure.py), cut tests |
| R6 | a clean binary–one-hot decoder has $3\cdot2^t-2-t$ workspace, $O(t)$ depth, and $O(2^t)$ size | [compiler theorem, Lemma P](COMPILER_THEOREM.md#lemma-p-clean-binaryone-hot-decoder) | [`tree_decoder.py`](../compiler_robust_hopf/tree_decoder.py), [`test_tree_decoder.py`](../tests/test_tree_decoder.py) |
| R7 | an explicit coherent router realizes the tail direct sum and clears all data, token, copy, and flag work registers | [compiler theorem, Lemma R](COMPILER_THEOREM.md#lemma-r-explicit-coherent-router) | [`router.py`](../compiler_robust_hopf/router.py), [`test_router.py`](../tests/test_router.py) |
| R8 | the maximal feasible cut attains $`O\!\left(n+\frac{N}{n+m}\right)`$ depth for every large workspace, including $s=1$ | [compiler theorem, Proposition R](COMPILER_THEOREM.md#proposition-r-maximal-cut-depth) | [`resource_bounds.py`](../compiler_robust_hopf/resource_bounds.py), broad-grid tests |
| R9 | parameter capacity and output light cones give matching size/depth bounds; fusing free one-qubit slots gives the separate CNOT lower bound for arbitrary clean workspace | [compiler theorem, Section 9](COMPILER_THEOREM.md#9-matching-lower-bounds), using the standard parameter method of [Iten et al., Section III](https://arxiv.org/html/1501.06911v4#S3) | analytic dimension proof and resource diagnostics |
| R10 | the real frame and phase-dressed complex magnitude frame attain the all-workspace optimum | [compiler theorem, main theorem](COMPILER_THEOREM.md#main-theorem-optimal-exact-hopf-frame-compilation) | [`unified_compiler.py`](../compiler_robust_hopf/unified_compiler.py), both resource ledgers |
| R11 | frame-safe compilation preserves the global record and introduces no additional asymptotic depth factor in the matched program | [QBP consequence](QBP_CONSEQUENCE.md) | decoder and boundary tests; reviewer walkthrough |

## 5. Fault-tolerant sources and contribution boundaries

The [fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md) uses a complete
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
| F12 | [Yamazaki–Akibue, arXiv:2603.14202v1](https://arxiv.org/html/2603.14202v1), Theorem 1, Section 3 and Theorem 4 | leading precision constants for complete controlled SU(2) synthesis | the dirty-assisted construction retains a precision-length clean instruction register; the ancilla-free construction has different scaling and a typical-target guarantee |
| F13 | [Yuan–Zhang–Zi, arXiv:2608.17846v2](https://arxiv.org/html/2608.17846v2), Theorem I.1 and Definition II.1; [Fang–Heunen–Wang, arXiv:2607.12907v1](https://arxiv.org/html/2607.12907v1), Theorem 1.2 and Corollaries 3.8–3.9 | current generic-unitary and near-Clifford comparisons; the [endpoint substitution](OPEN_PROBLEM.md#selection-audit-after-canonical-completion) keeps precision factors explicit | different target families and clean allocations; the initialized-isometry contract already appears in general-unitary synthesis |
| F14 | [Vasconcelos–Gilyén, arXiv:2507.07900v2](https://arxiv.org/html/2507.07900v2), Sections 2–3 and Appendix B | block-work uncomputation, exact-history lower bounds and approximate product compression | original clean work is still needed during a query; the lower bound is for the defined coherent-measurement class, and approximate compression needs near-identity dilations |
| F15 | [Ma–Joven–Liu, arXiv:2609.11153v1](https://arxiv.org/html/2609.11153v1), Theorem 3.1 | clean ancilla compression for block-encoding counting bounds | at most $`n+2T`$ clean ancillas, with possible normalization change; no two-clean unitary conclusion |
| F16 | [Motlagh–Pocrnic, arXiv:2605.20334v1](https://arxiv.org/html/2605.20334v1), Section II | improved dirty-QROM constants | the displayed output word remains initialized; no asymptotic or constant-factor lookup improvement is claimed here |
| F17 | [Lai, arXiv:1411.5408v3](https://arxiv.org/pdf/1411.5408v3), Theorems 1.6 and 1.8 | finite-tree Carleson embedding and L2 maximal inequality | the inequalities are inherited; the tree-transport note supplies the local overlap-defect identity and its weighted residual application, not a free block-encoding circuit |
| F18 | [Boyd–Vandenberghe, Convex Optimization](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf), Appendix A.5.5 | quadratic elimination and Schur complements | standard linear algebra; the weighted-block and residual-assembly notes supply their scalar tree recursions, local rejection corrections, and mode allocations |
| F19 | [Barenco et al., arXiv:quant-ph/9503016v1](https://arxiv.org/pdf/quant-ph/9503016v1), Lemma 4.1 and Section 8 | Euler factors, two-level unitary decomposition, and Gray-code basis routing | inherited generic synthesis; the local pricing uses the retained borrowed-sector echo and does not claim an optimal native implementation |
| F20 | [Brassard–Høyer–Mosca–Tapp, arXiv:quant-ph/0005055](https://arxiv.org/pdf/quant-ph/0005055), Section 2, Eqs. (7)–(8) | repeated amplitude amplification and its sine/cosine law | the five-call one-clean amplification is a specialization at amplitude $`\sin(\pi/10)`$; the local proof supplies the oblivious complete-isometry error and literal-phase accounting |
| F21 | [Vatan–Williams, arXiv:quant-ph/0308006](https://arxiv.org/pdf/quant-ph/0308006), Section III Theorem 1 and Section IV | magic-basis equivalence between SO(4) and two SU(2) factors | the reduction and fixed-size real two-qubit synthesis are inherited; the local Cayley formulas fix a literal Clifford convention and reuse the existing one-target compiler with a borrowed logical spectator |
| F22 | [Beals–Buhrman–Cleve–Mosca–de Wolf, Quantum Lower Bounds by Polynomials](https://homepages.cwi.nl/~rdewolf/publ/qc/polynomials.pdf), Lemma 4.1 | query-by-query amplitude-degree growth | the [Hopf scattering audit](ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer) derives its Laurent-degree specialization directly; it concerns unchanged local-angle queries, not native T counts |
| F23 | [Toulouse–Umrigar, arXiv:physics/0701039v2](https://arxiv.org/pdf/physics/0701039), Eq. (45); [Mitarai–Fujii, arXiv:1901.00015](https://arxiv.org/pdf/1901.00015), Fig. 1 | logarithmic wave-function derivatives and Hadamard-test interference | inherited estimator ingredients; the reference-state QBP note proves its own leaf scores, Hopf envelope, support bounds, and native resource comparison |
| F24 | [Möttönen–Vartiainen–Bergholm–Salomaa, arXiv:quant-ph/0407010v1](https://arxiv.org/pdf/quant-ph/0407010v1), Section III, Eqs. (4), (5), and (7) | uniformly controlled phase cascade and its residual arithmetic-mean common phase | inherited decomposition; the complex coarse proof fixes local signs and applies the retained exact-return dirty interpreter to determinant-one row words |

| F25 | [Ross–Selinger, arXiv:1403.2975v3](https://arxiv.org/abs/1403.2975v3), abstract and synthesis-runtime discussion | distinguishes short native words from efficient search; optimal synthesis uses a factoring oracle, and the efficient expected runtime without it is conditional | computational context only; the bounded-input proof uses F5 word-length existence and guarded exhaustive coarse enumeration, with no factoring oracle or runtime conjecture |
| F26 | [Baur–Strassen, *Theoretical Computer Science* 22(3), 317–330 (1983)](https://www.sciencedirect.com/science/article/pii/030439758390110X) | classical reverse differentiation background | the explicit Pauli baseline derives its Hopf reverse recurrence and dyadic error bound directly; reverse differentiation is not claimed as a new algorithmic principle |
| F27 | [Casas et al., arXiv:2602.05425v1](https://arxiv.org/html/2602.05425v1#S3.SS2.SSS2), Section III.2.2, Eq. (29) | exact denominator-exponent depth bound for Clifford-matchgate plus T layers | the [source-depth proof](SOURCE_T_DEPTH.md) applies this existing method to two particular sources, permitting signed-permutation Clifford stages including the parity-odd extension; no unrestricted depth lower bound follows |
| F28 | [Vasconcelos, arXiv:2609.34659v1](https://arxiv.org/html/2609.34659v1#S3.SS4.SSS1), Theorem 8; [Kim, arXiv:2506.15147v3](https://arxiv.org/pdf/2506.15147v3), Section 3 | current precision-depth comparisons | the former uses precision-sized clean work, the latter a prepared catalyst and leaves the clean/dirty-only constant-T-depth question open; neither is used as a theorem premise for the Hopf schedules |
| F29 | [Kim–Laakkonen, arXiv:2512.24982v1](https://arxiv.org/pdf/2512.24982v1), Theorems 3 and 5 | constant non-Clifford-depth control of CNOT and Clifford circuits without ancillas | the local controlled-shear lemma is a rank-sensitive specialization with a literal four-T-layer word and explicit Clifford ledger; constant-depth control is inherited |
| F30 | [Boyd, arXiv:2312.00696v2](https://arxiv.org/html/2312.00696v2), Section III and Appendix A | commuting SELECT/QROM groups and Clifford changes of basis for parallel action | its address copies use initialized registers; the local all-dirty shear echo and workspace allocation are proved separately |
| F31 | [Selinger, arXiv:1210.0974v2](https://arxiv.org/html/1210.0974v2), Proposition 5.1, Eqs. (17)–(21) | one-T-layer Pauli-conjugation normal form | the local rational-transfer specialization excludes T-depth one for full-input nonaffine classical permutations with arbitrary returned dirty helpers; initialized-clean-subspace implementations are outside that claim |
| F32 | [Takahashi–Tani–Kunihiro, arXiv:0910.2530v1](https://arxiv.org/pdf/0910.2530v1), Sections 2.1–2.3 | linear-size, linear-depth exact ripple-carry addition without initialized work | deleting the two gates targeting the arbitrary carry-output wire gives the modular adder used in the signed dirty increment and baseline sum tree; its literal CNOT/Toffoli word is emitted in the bounded checks; later clean-work/fanout constructions are not used |
| F33 | [Remaud–Vandaele, arXiv:2501.16802v2](https://arxiv.org/html/2501.16802v2), Lemmas 2/4, Algorithm 3, Theorem 2 | exact helper-free addition via shallow CNOT and Toffoli ladders | truncate the carry output at the abstract ladder level, then synthesize the shorter ladders; applies only to private counters; bounded checks audit the reduced macro, while the optimized ladder-depth bound is imported analytically |
| F34 | [Vandaele, arXiv:2603.12917v1](https://arxiv.org/html/2603.12917v1), Section 5, Theorem 4 and Corollary 7 | exact logarithmic-depth increment and controlled increment with one returned dirty helper | supplies the retained round-based compressor and the separate two-dirty-bit read-only increment; temporary control borrowing stays on private supports; the carry-pipeline refinement instead uses linear TTK arithmetic |
| F35 | [Aaronson–Gottesman, arXiv:quant-ph/0406196v5](https://arxiv.org/pdf/quant-ph/0406196v5), Section III; [Zhang–Zhang, arXiv:2409.13809v2](https://arxiv.org/html/2409.13809v2#S3.SS1), Theorem III.1, Eqs. (10)–(11) | stabilizer-overlap quantization and Pauli conjugation by one T layer into a Hermitian Clifford | the [two-layer source obstruction](SHALLOW_SOURCE_OBSTRUCTION.md) derives a full-space transfer alphabet and robust source witnesses; initialized-clean isometries and growing frame-depth lower bounds are excluded |

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
| R14 | shared-source composition of frame residuals | the uniform-precision construction combines accepted-branch shifts and a retained source with the inherited compression gadget F10 and amplification F11; see [theorem and proof](FAULT_TOLERANT_COMPILER.md) |
| R15 | matching T-count in the stated workspace regimes | the sufficient-clean [theorem](FAULT_TOLERANT_COMPILER.md) plus the F3/F4 reductions; the separate all-clean-budget [corollary and borrowed-workspace proof](BORROWED_WORKSPACE_COMPILER.md) also uses the arbitrary-budget compiler under its additional condition |
| R16 | fixed-parameter bounded-score robustness | [QBP approximation](QBP_APPROXIMATION.md): complete complex gradient, reflection sums, rounded weights, correlated dirty-bank reuse, and quantum/classical budgets; no derivative of a compiled word or gradient-query optimum |
| R17 | two-clean operator-source compiler and corollaries | [proof](OPERATOR_SOURCE_COMPILER.md): full-frame $O(N+nL)$ bound and dirty-bank refinement; matched literal diagonal and general one-target U(2) multiplexor frontiers; all have separate explicit reservations and include operator-core return error |
| R18 | exact geometric operator-source costs | [source proof](OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers): $2m-4$ uncontrolled and $2m-2$ controlled T gates, even with returned helpers; native words and transfer witnesses checked for small $m$ |
| R19 | sufficient-clean complex-frame extension | [Corollary 7](FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame): literal diagonal SU(2) embedding, same-pool composition, and diagonal-subfamily lower bounds |
| R20 | conditional-suffix compiler | [proof](CONDITIONAL_SUFFIX_COMPILER.md): $`O(N+L\ell_*(n))`$ T gates and $`O(NL)`$ Clifford gates for $`n\ge1`$, $`L\ge6`$, $`a=2`$, $`b\ge L+n+7`$; star residuals and streamed coarse programs use $`O(\log(s+1))`$ private clean work per group of $`s`$ depths |
| R21 | two-clean T-depth upper bound | [schedule](T_DEPTH_COMPILER.md): $`D_T=O(NL/b+L\ell_*(n)+n^3)`$, $`T,G=O(NL)`$, at $`b\ge2(L+n+7)`$; combines inherited F2 routing and F8 predicates with the complete-frame source, coarse-program, and conditional-work ledgers; no matching depth claim |
| R22 | structured residual and scoped endpoint diagnostics | [source-reuse analysis](SOURCE_REUSE_LIMITS.md) and [tree transport](ENDPOINT_TREE_TRANSPORT.md): linear generators, normalized transport, weighted norms using F17, finite-order correction witnesses, and Pauli-routed leakage; the [Frobenius-coarse specialization](ENDPOINT_TREE_TRANSPORT.md#frobenius-small-residuals-also-fit-the-endpoint-budget) is a direct consequence of the retained borrowed compiler, not a new frontier; the inherited Pauli-transfer method prices one paired-source transformed mask; a literal shared-conjugator fork has a coefficient-ellipse obstruction even under mask retuning, derived directly from the R25 Majorana moments; these are scoped interface restrictions, not additive or unrestricted frame lower bounds |
| R23 | simultaneous two-clean count and depth | [parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md): a single X conjugated by the F2 router gives a scratch-free indicator with linear address T-depth; F8 predicates and the full-frame ledger reduce the fixed-accuracy bound to $`O(n^2)`$; complete real-frame composition retains count-efficient banks under a sufficient square-root-scale width condition; no new general lookup tradeoff or optimal T-depth claim |
| R24 | forward weighted transport component | [weighted block](WEIGHTED_TRANSPORT_BLOCK.md): exact weighted Gram witness, scalar recursion using F18, complete one-signal-flag dilation, and certified singular-case preprocessing; F19 and the borrowed reflection interpreter give $`O(L\sqrt N)`$ native T count with exact dirty return, while its borrowed-signal extension gives $`O(N+nL)`$ with one block signal and approximate dirty return; routing is charged, and neither improves the full-frame endpoint bound |
| R25 | one-clean operator-source compiler and corollaries | [proof](ONE_CLEAN_COMPILER.md): conjugated scalar-source word and Pauli-routed anticommutator, native head-and-two-tail source, and fixed five-call amplification using F20; real grouped/banked bounds extend to $`a=1`$, as do phase-dressed complex magnitude frames by sequential composition at their separate threshold; X symmetry gives a zero-clean layerwise real-frame corollary; matched literal-diagonal and complete U(2) multiplexor bounds retain their own reservations; F2/F9 supply inherited lookup and Clifford-algebra ingredients, and no new T-depth claim is made |
| R26 | complete two-clean residual assembly | [proof](RESIDUAL_ASSEMBLY.md): affine forward dilation using F18, native controls using F19 and the borrowed-signal compiler, and two-term selection/amplification using F11 give $`T=O(N+nL)`$, $`G=O(NL)`$ at $`a=2`$, $`b\ge L+n+7`$, including leakage and dirty return; the full-mode hierarchy and elementary coupled completion have complete recursive contracts; direct block algebra gives a rank-at-most-four commutator repair without a mode-gap assumption, exact child-call cancellation, and a weighted repair-error bound; the surviving local wrappers retain their precision cost, so no grouped endpoint improvement or priority claim for the elementary completion follows |
| R27 | native-baseline antichain correction | [proof](ANTICHAIN_COMPILER.md): exact descendant-forest conjugation and prefix-pair packing reduce promised antichain-only changes to one addressed SU(2) table; the retained borrowed interpreter and three Euler factors, using F19 and the borrowed-signal refinement of R25, give $`T=O(N+L)`$, $`G=O(NL)`$ at $`a=0`$, $`b\ge L+n+7`$; the supplied baseline has determinant-one native local words of length $`O(n-d+1)`$ at depth d, and the target agrees with it literally outside the antichain; neither generic rounding, a general complex-magnitude extension, nor closure of the unrestricted endpoint is claimed |
| R28 | native-baseline sparse nested updates | [proof](SPARSE_UPDATE_COMPILER.md): exact off-support forest factorization and charged basis packing reduce an ancestor-closed set S to $`m=\lvert S\rvert+1`$ modes; for $`s=\lceil\log_2m\rceil`$ and $`n\ge2s+32`$, a shared dense column-residual dictionary and the R17 scalar-source block, using the inherited lookup, routing, and normalization-two amplification ingredients F2/F11/F19, give $`T=O(N+L)`$, $`G=O(NL)`$ at $`a=1`$, $`b\ge L+n+7`$; the complete initialized-isometry bound includes approximate work return and arbitrary dirty references; all changes on one root-to-leaf path are covered with a finite small-n fallback, but the unrestricted frame endpoint remains open |
| R29 | native source-width boundaries | [source-width proof](SOURCE_REUSE_LIMITS.md#6-changing-source-width-without-renewing-its-preparation): the reversed geometric loader retains the coefficient grid and gives $`O(L+n)`$ T gates and $`O(NL)`$ Clifford gates for the outer loaders and monotone width bridges; the transformed group bodies remain separately charged; inherited Pauli-transfer and algebraic-norm methods price the original bridge and one mask from a legal grouped coefficient; [native checks](../tests/test_precision_carry.py) support the literal phases, bridge direction, and full-input identities, without improving the compiler frontier |
| R30 | correlated flag/source boundary | [correlated-boundary proof](SOURCE_REUSE_LIMITS.md#7-a-flag-correlated-source-boundary-and-its-query-cost): charged initialization and decoding, direct code identities, and native physical width changes include arbitrary dirty references and released tail wires; a legal grouped mask produces code leakage, and complete coherent two-syndrome renewal reduces to the separately priced returned-source operation; [code checks](../tests/test_correlated_precision_carry.py) include actual inverses and approximate flag leakage; this restriction concerns the full renewal interface, not accepted-only blocks or arbitrary carried-source compilers |
| R31 | small-product residual coordinates | [product-first examples](SOURCE_REUSE_LIMITS.md#8-small-products-compress-before-synthesis) reuse the existing fixed-address one-target compiler; [Cayley recursion](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation) derives constant-size local data for the complete complex-coarse residual using standard Cayley/low-rank matrix algebra, with uniform conditioning and stable inverse conversion; [small-product](../tests/test_small_product_compilation.py) and [Cayley](../tests/test_tree_cayley.py) fixtures check complete matrices; no free coherent evaluator, improved native T-count, or priority for the classical transform is claimed |
| R32 | native four/eight-mode benchmarks | [four-mode primitive](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark) combines F21 with the retained one-target compiler, certified factor coordinates, and a complete borrowed-spectator ledger; [eight-mode coupling](ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling) gives a charged four-Pauli construction and a scoped product-factor distance; [inverse audit](ENDPOINT_TREE_TRANSPORT.md#9-inverse-cayley-recovers-the-original-local-scattering-word) reconstructs the original target wrapper; no generic endpoint improvement or additive precision lower bound follows |
| R33 | joint changing-target source word | [shared body](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets) reuses R25's scalar source and amplification, with an unconditional fixed conjugator, complete SU(2) signal error, and a borrowed-signal parity return; direct word cancellation gives $`27g+6`$ source appearances, while fixed-mask tail commutation leaves both comparisons the same leading precision term; eight-/sixteen-mode allocations and [native checks](../tests/test_joint_source_body.py) include queries, helpers, literal inverses, and complex coarse baselines; no generic T-count improvement or optimized-word minimum is claimed |
| R34 | canonical scalar compatibility with grouped SELECT | [group audit](CONDITIONAL_SUFFIX_COMPILER.md#11-a-canonical-scalar-fits-the-group-interface-but-retains-its-precision-charge) instantiates R33 on the existing scalar flag, with a distinct atom flag and borrowed signal, full-operator perturbation, and the retained group workspace slack; direction-controlled Z selects actual inverses for both the canonical and old scalar words; [finite checks](../tests/test_grouped_scalar_completion.py) cover the complete SELECT and amplification interfaces; the displayed precision recurrence remains linear in group count, with no improved endpoint or new lower bound |
| R35 | packed Hopf scattering and its boundary distinction | [tree-scattering proof](ENDPOINT_TREE_TRANSPORT.md#11-a-packed-hopf-scattering-step-and-its-boundary-transfer) gives complete root/marker ports, standard feedback elimination, and a charged $`O(N+L)`$-T step using R25 and exact borrowed routing; F22's degree method tests the unchanged-coin conversion via the Hopf path amplitude; [small checks](../tests/test_hopf_scattering.py) verify full transfers and native port permutations; no priced feedback implementation or stronger native T lower bound follows |
| R36 | two-clean state-only preparation | [state-only proof](STATE_ONLY_COMPILER.md) combines the actual borrowed coarse circuit with one globally computed residual table, R25's full-operator SU(2) synthesis, and F20's one-step amplification at amplitude one half; it gives $`O(N+L)`$ T count for $`L\ge\max\{6,n\}`$ with exact stated dirty allocation and approximate joint return; no other system columns are prescribed |
| R37 | raw gradients from reference-state interference | [decoder proof](REFERENCE_STATE_QBP.md) specializes F23 to all real Hopf derivatives and all allowed observables, gives a division-free derivative-envelope recurrence and a bounded-branch corollary, and charges R36 plus the observable; generic sampling overhead can grow as n, and no complete-frame or total-gradient optimality follows |
| R38 | all-angle coarse-frame gradient decoding | [coarse-frame proof](COARSE_FRAME_QBP.md) combines R36's common-coarse reference preparation, F23's X/Y interference, and Walsh spreading with exact classical correction from the actual native C; the depth-record norm is at most five, with $`O(N+L')`$ T-count per execution apart from the observable and $`O(S+Nn)`$ classical arithmetic after preprocessing; no fine full-frame or total-gradient optimality claim |
| R39 | gauge-fixed complex state-based QBP | [native interface](COMPLEX_COARSE_COMPILER.md) applies the borrowed reflection interpreter to F24's determinant-one phase tables, giving an actual logical coarse C with exact dirty return and $`O(N)`$ T count; R36's complex residual uses the same two flags and $`b\ge L+n+7`$; [two-stream proof](COMPLEX_COARSE_QBP.md) supplies magnitude and phase gradients with the R38 depth bound and $`O(S+Nn)`$ reconstruction; the literal original common frame phase and fine complete frame are not compiled |
| R40 | banked state preparation and fair QBP cost comparison | [banked state corollary](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation) reuses R17's exact whole-word dirty-bank query, inherited from F2, in R36/R39's preparation: $`T=O(\sqrt{NL}+L+NL/b+n\sqrt N)`$, $`G=O(NL)`$, at two clean flags, $`L\ge\max\{6,n\}`$, and $`b\ge2(L+n+7)`$; [gauged borrowed baseline](QBP_COST_COMPARISON.md#2-a-gauged-complex-borrowed-frame-baseline) combines the retained interpreter and F24's phase cascade with weighted row precision to give $`T=O(NK/(n+b)+K\sqrt N)`$, $`G=O(NK)`$, exact dirty return, and zero compiler clean work for $`b\ge2n`$; the [task comparison](QBP_COST_COMPARISON.md) retains separate frame/state precisions, literal workspace thresholds, observable costs, sampling, and preprocessing; it compares constructive upper bounds without a literal common-phase-frame or end-to-end gradient optimality claim |
| R41 | bounded-input construction and classical comparison | [algebraic residual coefficients](RESIDUAL_TABLE_PREPROCESSING.md) use standard half-phase identities with one shared root, certified rational intervals, and a finite-radius cutoff to preserve R36/R39's error constants without Euler search; [bounded-input audit](BOUNDED_INPUT_QBP.md) combines F5 existence with coarse enumeration, explicit masks and instruction output to prove polynomial construction for the listed grouped/state alternatives; the banked small-system source uses $`P+13`$ dirty wires. Deterministic and term-sampled classical Pauli baselines are charged; no unconditional efficient fine-word search, generic Euler-runtime theorem, or end-to-end quantum advantage is claimed |
| R42 | state-based QBP T-depth composition | [depth proof](STATE_QBP_DEPTH.md) composes R21/R23's exact schedules, inherited from F2, with R36/R39's constant number of residual rotations and the actual exact-return coarse interpreter; with $`B_0=P+n+7`$, $`b\ge2B_0`$ gives $`D_T=O(NP/b+P+n^3)`$, $`T,G=O(NP)`$, while $`b\ge16(B_0+\sqrt{NP})`$ gives one circuit with $`T=O(\sqrt{NP}+P+n\sqrt N)`$, $`G=O(NP)`$, and $`D_T=O(P+n^3)`$; both real/complex task streams and oracle depth are charged; no new lookup primitive, depth optimality, total-runtime gain, or general emitter is claimed |

The [Hopf error audit](HOPF_ERROR_ACCUMULATION.md) derives a sharp ideal-angle
stability recurrence and an exact finite relative-spectrum recursion from
the inherited nested frame supports H8, using standard block-matrix
spectral algebra. A midpoint-grid corollary is restricted to independent
nearest angular rounding. The audit also applies the existing scalar
source and amplification identities to an explicit family with coherent
linear leakage on reused flags. It does not import a general composition
lower bound or claim an unrestricted frame-depth obstruction. Its finite
fixtures use one common reduced source algebra, not a new native compiler.

The [four-echo audit](HOPF_FLAG_ECHO.md) derives exact compressed and
full-isometry errors from those same source identities. The positive
[radial filter](HOPF_RADIAL_FILTER.md) imports
[Grover's fixed-point phase sequence](https://arxiv.org/abs/quant-ph/0503205),
Eq. (1) and Section 3, and F5's determinant-one single-qubit word-length
theorem. Fixed-point amplification and ancilla-free phase approximation
are inherited. The local work supplies the literal global correction,
the full error including the residual accepted phase, a two-flag native
phase ledger, and the complete-frame precision cap. Failure-probability
suppression is not identified with phase-sensitive operator error, and
no new asymptotic count/depth frontier or fine-search runtime is claimed.

The [amortized dirty-lookup theorem](AMORTIZED_DIRTY_LOOKUP.md) combines
F2/F8's dirty traversal and echo with F29/F30's controlled-linear and
commuting-operator techniques. Its query interface is retained by the
[variable-width hybrid](PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision).
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

The [bilinear query](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction)
specializes F2's indicator/bilinear framework. The
[dirty-counter indicator](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters)
uses exact sum and cyclic-rotation echoes. Its retained signed increment
applies two F32 adders and Clifford gates. The two-dirty-bit refinement
uses modular negation as an involution and imports F34's completed
controlled increment: public literals touch only CNOTs, while the source
borrows only private controls. F33 supplies the retained baseline
sum-tree adder.
The [masked carry pipeline](DIRTY_SUM_COMPRESSION.md) uses F32's linear
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

The [exact source-depth certificate](SOURCE_T_DEPTH.md) specializes F27
to the existing geometric and paired sources, with a parallel paired-tail
schedule. Its lower bounds are for the stated Majorana-layer architecture
and exact targets; it changes no asymptotic frame or state-QBP theorem.
The [precision-depth comparison](RELATED_WORK.md#16-precision-depth-and-workspace-assumptions-2-october-2026)
records F28's distinct workspace contracts.
The [two-layer source obstruction](SHALLOW_SOURCE_OBSTRUCTION.md)
combines F35's standard facts with exact geometric-source coefficients.
The [conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md)
instead supplies an explicit reversible prefix preparation and its
active-suffix block-encoding contract, using the existing exact lookup,
native controlled-H, reflection, and amplification ingredients. By itself
it changes the precision component only; it is not all-dirty source
resynthesis.

The [grouped program construction](GROUPED_PROGRAM_PREFETCH.md) combines
that source with exact conditional prefetch, read-only program copies,
private conjunction trees, and an internal-enable phase correction.
The [chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md) uses the existing
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

The [incremental selectors](GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
reuse this same conditional-work interface and exact native Toffolis.
Their local dependency schedule reduces selector/enable maintenance to
linear group depth within the original reservation. The
[common-source audit](GROUPED_PROGRAM_PREFETCH.md#11-a-common-source-identity-and-the-remaining-reflection)
uses direct unitary conjugation and the existing normalization-two block:
boundary preparations cancel, but conjugated reflections and source-bank
return remain charged. The legal zero-angle witness excludes only the
stated stale-monitor and one-use-bank substitutions. Neither statement
imports a new synthesis premise or improves the global depth order.

The [consolidated state-based QBP theorem](STATE_BASED_QBP_THEOREM.md)
collects R36 and R38–R42 under one input, precision, workspace, and
sampling contract. It introduces no additional compiler bound or
optimality claim. R37 remains an optional reference-state predecessor.

The [native coarse-frame example](NATIVE_COARSE_QBP.md) is finite implementation
evidence for R38 and R36's small-register fallback. Its commutator distance
certificate, per-reflection dirty echo, and literal gate ledger are stated
explicitly. It does not add a new asymptotic compiler theorem or replace
the existing analytic proof with numerical scaling evidence.

The [native complex extension](NATIVE_COMPLEX_COARSE_QBP.md) supplies finite
implementation evidence for R39: elementary all-suffix phase selection,
complete coherent preparation, both gradient streams, and a literal cost
comparison. Its sharper K2 trace and rational K3 certificate specialize
the same elementary commutator algebra. They do not change the asymptotic
compiler frontier or establish a general native residual-table emitter.

The [certified native residual rows and tables](NATIVE_RESIDUAL_ROTATION.md) implement
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
[selected high-precision endpoint](OPEN_PROBLEM.md).
At $`b\ge2(L+n+7)`$, the grouped compiler's banked form gives
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ T gates and $`O(NL)`$ Clifford gates.

The two-clean proof uses the
[two-pass dirty lookup](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
and [borrowed predicate toggle](BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
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
are separately compared with GKW and Yamazaki–Akibue in
[related work, Section 12](RELATED_WORK.md#12-contemporary-comparisons-and-the-broader-compiler-contribution).

## 6. Evidence classification

The local evidence separates the following kinds of support.

| Level | Examples |
|---|---|
| analytic identity | addressed layers, echo sectors, tree cut, lower bounds |
| explicit reversible construction | binary–one-hot decoder and coherent router |
| imported exact synthesis | elementary UCG and multi-controlled-X circuits |
| finite regression evidence | matrix equality, entangled-input routing, cleanup, and resource ledgers |
| fault-tolerant analytic construction | compact source, complete sampler block, shared-source composition, and charged resource proof |
| fault-tolerant finite checks | exact small circuits and arithmetic certificates described in the [verification index](../verification/fault_tolerant/README.md) |

The finite checks are designed to expose convention, indexing, order, phase,
cleanup, and resource errors. They do not replace the asymptotic proofs. The
fault-tolerant checks are not a general emitted Clifford+T frame compiler.

## 7. Provenance records

- [`provenance/upstream.json`](../provenance/upstream.json): upstream commits,
  file lineage, and reconciliation;
- [`provenance/literature.json`](../provenance/literature.json): compiler source
  versions and role assignments;
- [`SYNC.md`](../SYNC.md): human-readable synchronization policy.

These records preserve provenance and source discipline; they are not additional
scientific assumptions.

The [borrowed-workspace appendix](BORROWED_WORKSPACE_COMPILER.md) condenses
the original research derivations in `Hopf_Fault_Tolerant_Research.md`,
Sections 4.5–4.6, 4.8 and 5.19, and `Hopf_Small_Clean_Workspace.md`.
Those filenames identify its historical source, not additional proof or runtime
dependencies. The appendix contains the complete real-frame argument used by
this publication, including the dirty lookup, exact sector echo, resource sum,
and restricted matching splice.

---

[← Verification](VERIFICATION.md) · [Complete narrative](../REVIEW.md) · [Related work →](RELATED_WORK.md)
