# Proof and executable correspondence

[← QBP consequence](QBP_CONSEQUENCE.md) · [Complete narrative](../REVIEW.md) · [Source map →](SOURCE_MAP.md)

The repository uses executable checks to test the conventions and constructions
that are easiest to get wrong: column labels, chronological order, relative
phases, coherent routing, work-register cleanup, and simultaneous resource
peaks.  The asymptotic theorem itself is established analytically.

## 1. Evidence levels

| Level | Meaning in this repository |
|---|---|
| analytic proof | a dimension-independent operator or resource argument |
| explicit construction | a complete reversible or logical-gate schedule is supplied |
| imported exact synthesis | an elementary compiler theorem is used under its stated model |
| finite regression check | independently constructed matrices, states, schedules, or ledgers are compared |

These levels are kept distinct.  A matrix test does not prove an asymptotic
bound, and an asymptotic bound does not certify an implementation's bit order.

## Claim-to-evidence index

The selected [Results A–D](../manuscript/PUBLICATION_SCOPE.md) are analytic
theorems. Their finite fixtures check conventions, literal phases, work
return, and resource accounting under the individual recorded contracts.

| Claim | Analytic proof | Principal executable evidence | Boundary |
|---|---|---|---|
| A: exact complete frame | [Compiler theorem](COMPILER_THEOREM.md) | Exact-frame, strict-zero, decoder, router, and resource suites in §§2–8 below | UCG and multi-controlled-X elementary synthesis are imported |
| B: sufficient-clean matching T count | [Fault-tolerant compiler](FAULT_TOLERANT_COMPILER.md) | [Four exact receipt suites](../verification/fault_tolerant/README.md), resource ledgers, and full-input approximation checks | Finite kernels and analytic resource proxies; no general native frame emitter |
| C: one-clean count frontier | [Growing-source theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count), [banked frontier](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier), [one-clean baseline](ONE_CLEAN_COMPILER.md) | [Operator-source checks](../tests/test_operator_source_compiler.py), native phase and query checks; bounded growing-source, conjugate-recovery and bank-allocation diagnostics | Finite emitted components and supplied diagnostic matrices; no general native frame emitter; full joint return error is included |
| D: one-clean optimal count and T-depth | [Precision-uniform subtree theorem, R79](../research/depth/UNIFORM_SUBTREE_DEPTH.md), using [R78's decoder](../research/depth/SUBTREE_PROJECTOR_DECODER.md) and [R72's one-clean interfaces](../research/depth/ONE_CLEAN_DEPTH.md) | [Exact subtree and bounded spectral-decoder checks](../tests/test_subtree_projector_decoder.py), [one-clean native words](../tests/test_one_clean_compiler.py), and retained shared-binary-source, conjugate-recovery, blocked-query, dirty-indicator, and uniform-allocation suites | Analytic absolute precision ledger, enlarged protected bank, terminal and width/precision splices, and full-isometry composition; no new circuit identity or fixture. Finite checks do not certify unbounded precision, asymptotic cutoffs, or emit every variable-size frame |
| All-precision general depth R96 | [All-precision depth theorem](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md), using [the global stage](../research/depth/GLOBAL_BORROWED_PRECISION_STAGE.md), [fast references](../research/depth/FAST_CAPPED_GEOMETRIC_REFERENCE.md), [polynomial-precision bridge](../research/depth/POLYNOMIAL_PRECISION_DEPTH.md), and [direct-prefix terminal](../research/depth/DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md) | [Four bounded capped-reference checks](../tests/test_fast_capped_geometric_reference.py); retained native grouped-terminal controls, growing-source/recovery, R91 query and R94 payload checks | The global address/reserve ledger, common fine grammar, separate precision/coarse support, instruction schedule, simultaneous widths, all-precision splice and full-isometry error are analytic. Dense comparison completions are verification objects; no new general emitter or scientific receipt |
| Logarithmic-precision general depth R95 | [Full-flag depth theorem](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md), using [sparse correction](../research/depth/COARSE_FLAG_SPARSE_CORRECTION.md), [small-workspace sources](../research/depth/SMALL_WORKSPACE_DESCENDANT_SOURCE.md), and [the relaxed-count local compiler](../research/depth/RELAXED_COUNT_LOCAL_DEPTH.md) | [Four bounded sparse-correction checks](../tests/test_coarse_flag_sparse_correction.py) for the live root, signed Schur sources, complete-work recovery, actual inverses and leakage; retained native/source, R91 query and R94 payload tests | Dense source completions are verification objects. Native source scheduling, the two-stage local construction, catalogue absorption, full-isometry composition and uniform resource bounds remain analytic; no general native emitter or new scientific receipt |
| Fixed-accuracy logarithmic general depth R90 | [Dirty catalogue terminal](../research/depth/DIRTY_CATALOGUE_TERMINAL.md), using R89 coarse reductions | [Bounded catalogue checks](../tests/test_dirty_catalogue_terminal.py) of dirty program selection, two-array commutator return, literal native factors, reference inputs, and the one-array failure, plus retained coarse/query evidence | Analytic fixed-number reduction, zero-clean local synthesis, selected-factor error, catalogue absorption and literal-width splice; fixed-L constants. No efficient preprocessing, depth optimality, or scalable native emitter is certified by these finite checks |
| Retained causal compact-coin theorem R89 | [Causal compact-coin proof](../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md), with compact coins and subset primitives | [Four bounded causal-batching fixtures](../tests/test_causal_compact_coin_depth.py), [seven compact-coin checks](../tests/test_compact_coin_flags.py), and [two subset-routing checks](../tests/test_short_label_subset_routing.py), plus retained source, query, routing, and terminal evidence | Analytic catalogue existence, complete coarse flags, native simultaneous ledger, prefix recursion, full-isometry error, and literal-width splice; fixed-L constants. Its coarse reductions support R90 |
| Retained product-coin theorem R87 | [Product-coin proof](../research/depth/COIN_DILATION_DEPTH.md), using retained rounding and selected routing | [Seven bounded coin-dilation fixtures](../tests/test_coin_dilation_depth.py), plus retained source, query, coarse-recovery, and terminal checks | Analytic product partition, spectral-band readout, paid source and actual inverse, full-isometry error, and original literal-width splice; fixed-L constants. Its fine decoder supports R89 |
| Retained terminal R88 and contiguous routing | [Prefix-selected terminal proof](../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md), [contiguous-layout specialization](../research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals) | [Three bounded selected-terminal fixtures](../tests/test_prefix_selected_coin_terminal.py) and [three logical routing fixtures](../tests/test_interval_barrel_routing.py) | Analytic recursive errors, native query/source and routing ledgers, literal-width splice, and arbitrary dirty/reference contract; sharp queries and fallback support R89's stronger terminal bound |
| Interval-dilation predecessor R86 | [Interval-dilation proof](../research/depth/INTERVAL_DILATION_DEPTH.md), with [selected dirty routing](../research/depth/SELECTED_DIRTY_ROUTING.md) | [Seven bounded interval-dilation fixtures](../tests/test_interval_dilation_depth.py), plus retained source, query, recovery, and terminal checks | Analytic absolute-norm Haar estimate, paid local compiler, global sign/zero rounding, full-isometry composition, and literal-width splice; fixed-L constants. Bounded matrices do not emit the general native compiler or prove asymptotics |
| Retained linear-precision depth predecessor | [Complete construction](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision), [shared-gradient recovery](../research/depth/SHARED_GRADIENT_RECOVERY.md), [shared binary terminal source](../research/depth/SHARED_BINARY_TERMINAL.md) | [Bounded all-column recovery](../tests/test_sublinear_depth_recovery.py), [native calibration](../tests/test_phase_gradient_calibration.py), [shared-gradient checks](../tests/test_shared_gradient_recovery.py), [shared binary terminal checks](../tests/test_shared_binary_terminal.py), [native coarse-program checks](../tests/test_grouped_terminal_depth.py), and retained query/indicator suites | Absolute constants in the stated linear-precision range; analytic source, invariant-sector recovery, uniform workspace splice and native schedule; finite checks are not a general compiler |
| Fixed-accuracy polylogarithmic corollary | [Subtree-projector decoder, R78](../research/depth/SUBTREE_PROJECTOR_DECODER.md) | [Exact subtree and bounded spectral-decoder checks](../tests/test_subtree_projector_decoder.py), plus retained source/query/one-clean suites | Analytic full-mode reflection, reversible search, native ledger and width splice; retained as the decoder proof and fixed-accuracy consequence of R79. Finite fixtures do not emit a native descendant loader or prove asymptotics |
| Entry-indexed normalization barrier | [Encoding limit](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md) | Exact Haar spectrum, explicit optimal factorization and nuclear-norm dual witness | Analytic proof for the stated factorization contract; no unrestricted native-depth or recovery-query lower bound |
| Zero-clean five-layer obstruction | [Operator rigidity, R75](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md) | [Exact image orbit, native split words, Clifford boundaries, and common-commutant audits](../tests/test_four_layer_operator_rigidity.py) | Uniform operator-norm theorem, arbitrary dirty width; unrestricted positive-clean five-layer and scalar-gap questions remain open |
| Fixed-clean four-layer rank requirement | [Initialized-sector rigidity, R74](../research/depth/FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector) | [Symplectic, weighted-phase, and initialized-projector audits](../tests/test_four_layer_crossing_rank.py) | Explicit analytic minimum-rank bound $`\Omega_a(\sqrt L)`$; no unrestricted four-layer exclusion with positive clean work |
| Restricted one-clean five-layer obstructions | [Midpoint rigidity and finite projective order, R76](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md) | [Phase descent, native midpoint, compressed-commutant, finite-boundary and cyclotomic checks](../tests/test_one_clean_midpoint_rigidity.py) | Pauli or Hermitian-Clifford encoded midpoint with unrestricted crossing rank; separate finite-projective-order occupied-block obstruction; neither condition follows from terminal return |
| Non-Clifford midpoint obstruction at bounded inner spread | [Inner clean-spread rigidity, R77](../research/depth/ONE_CLEAN_INNER_SPREAD.md) | [Flat expansions, exact sector moments, literal four-layer return, and finite-boundary audits](../tests/test_one_clean_inner_spread.py) | Explicit gap at each spread t; necessary minimum $`\Omega(\log\log L)`$ spread; exact common-flip counterexample has constant occupied defect and does not refute the unrestricted conjecture |
| One-clean obstruction at bounded outer-axis misalignment | [Outer-overlap rigidity, R97](../research/depth/OUTER_OVERLAP_RIGIDITY.md), with [representation limits](../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md) | [Nine exact fixtures](../tests/test_outer_overlap_rigidity.py): physical-axis signs and routing, transverse midpoint gauge, small parity families and center mixing, literal three-layer resplitting, and every occupied entry and determinant of the descent counterexample | The gap at bounded misalignment and uniform aligned threshold are analytic; finite fixtures check native identities and scope, not the unrestricted gap. Displayed complexity does not exclude another split; two added dirty wires remove only the demonstrated determinant obstruction |

<a id="endpoint-count-evidence"></a>

<a id="two-clean-count-evidence"></a>

### One-clean count evidence

Result C's [one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
and [banked frontier](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
cover every logical column, literal phase, returned clean work, arbitrary
dirty input and external reference. The [claim ledger](../manuscript/PUBLICATION_SCOPE.md)
states the exact thresholds, matching ranges and sharper two-clean complex endpoint.

The catalogue's [count-source checks](reference/VERIFICATION_CATALOGUE.md#count-source-checks)
record non-dyadic precision scales, phase dressing and conjugate recovery,
including supplied-matrix and native-emission limits. Those finite checks
do not establish the all-size compiler or its resource bound.

### Subtree-projector decoder evidence

The catalogue's [complete-frame depth checks](reference/VERIFICATION_CATALOGUE.md#complete-frame-depth-checks)
record R78–R96's decoder, orientation, dilation, routing, sparse-correction
and precision fixtures. Each retains its own negative controls, finite
range, return contract and analytic-versus-implemented boundary. Result D
uses R78's retained decoder; R79's all-precision composition adds no new
circuit identity or fixture.

### Sublinear-depth evidence

The catalogue's [shared-source and recovery checks](reference/VERIFICATION_CATALOGUE.md#shared-source-and-recovery-checks)
record the sublinear predecessor, native calibration, shared-gradient,
grouped-terminal and shared-binary-source fixtures. Their bounded words
and supplied table actions do not emit a complete native descendant
loader or certify asymptotic cutoffs and preprocessing.

### Detailed evidence records

The [verification catalogue](reference/VERIFICATION_CATALOGUE.md) preserves
the per-fixture descriptions, negative controls, and software boundaries for
the supporting constructions and earlier routes. Its
[detailed finite checks](reference/VERIFICATION_CATALOGUE.md#detailed-finite-checks)
cover packed diagonals, source borrowing, tree reductions, refinement and
forest compilers, operator stability, and representation restrictions. The
[endpoint source and representation map](reference/VERIFICATION_CATALOGUE.md#endpoint-source-and-representation-checks)
links the additional source, bank, recovery, encoding, and reflection
fixtures to their proof and module homes. These records retain their stated
interface restrictions; they are not additional claims of A–D or unrestricted
circuit lower bounds.

The separate [state-based QBP coverage](reference/VERIFICATION_CATALOGUE.md#state-based-qbp-coverage)
records its changed decoder and bounded native integration. Use the
[test index](../tests/README.md) to find executable suites and
[§9](#9-reproduce-the-checks) to reproduce the checks.

## 2. What is represented locally

| Component | Local representation | Principal executable check | Imported ingredient |
|---|---|---|---|
| Hopf states and frames | dense matrices from independent recursive and addressed-layer constructions | complete equality, orthogonality, state and marker columns, chart domains, singular coordinates | Hopf geometry inherited from the earlier work |
| strict-zero echo | dense logical gates and exact reversible permutations | all four sectors, every nonfinal depth, full frames, inverse, and complex composition | UCG and ancilla-free MCT resource theorems |
| binary–one-hot decoder | explicit X/CNOT/Toffoli layers | basis action, arbitrary-basis reversibility, disjoint layers, counts, and clean return | constant-cost elementary decompositions |
| coherent router | explicit CNOT-fanout and Fredkin layers with sparse complex-state simulation | basis routing, entangled inputs, tail direct sum, complete cut, and zero leakage | coherent copy and constant-cost controlled gates |
| controlled subtree frames | exact logical token/flag-controlled rotations | equality to the ideal block-diagonal tail and flag cleanup | UCG and MCT synthesis |
| leaf-phase diagonal | exact block-diagonal UCG matrix | complete diagonal, inverse, common phase, and complex magnitude composition | UCG synthesis |
| resource theorem | integer and exact-rational ledgers | workspace peaks, schedule dispatch, endpoints, geometric sums, and cut inequalities | asymptotic primitive bounds |
| gradient records | parity, signed-histogram, and fast Walsh–Hadamard decoders | agreement of routes, empirical means, and deterministic record norms | QBP record identities |

The repository does not reproduce the elementary UCG or multi-controlled-X
compiler.  Those exact results are imported from the all-workspace
state-preparation framework.  Toffoli, Fredkin, controlled one-qubit gates, and
fixed-width controlled Givens rotations have exact constant-size,
constant-depth decompositions in the declared arbitrary-one-qubit+CNOT model.

## 3. Geometry and compiler contract

The frame suite checks:

- equality of the recursive and addressed-layer frames;
- real-frame orthogonality;
- phase-dressed complex magnitude-frame unitarity;
- the breadth-first marker convention;
- $g_{j,j}=a_j^2$ for unrestricted angles;
- $a_j=\sqrt{g_{j,j}}$ on the canonical domains;
- tolerance-aware regularity at floating-point chart boundaries;
- zero raw derivative and a unit chart-selected marker continuation at a
  singular coordinate.

The compiler-boundary fixtures check:

- two unitaries with the same prepared state but different marker columns;
- the resulting gradient change
  ```math
  (2,0,0)\longmapsto(0,\sqrt2,0);
  ```
- a checkpoint suffix that preserves one state but changes a derivative;
- an active-interface-safe suffix that preserves designated means without
  preserving the complete distribution.

Additional boundary tests exercise the sharp worst-observable sensitivity,
common-phase cancellation, projected marker error versus physical leakage,
and singular-marker ambiguity across all two-qubit Pauli observables.

Files:

- [frame implementation](../compiler_robust_hopf/frames.py)
- [frame tests](../tests/test_frames.py)
- [complex geometry tests](../tests/test_complex_analysis.py)
- [compiler-boundary implementation](../compiler_robust_hopf/compiler_boundaries.py)
- [compiler-boundary tests](../tests/test_compiler_boundaries.py)

## 4. Strict-zero construction

The strict-zero suite checks:

1. $C^2=R_y(\theta)$ and $XCX=C^{-1}$;
2. all four $(h,b)$ sectors;
3. exact restoration of the borrowed logical suffix bit;
4. absence of hidden workspace;
5. every nonfinal addressed depth through $n=8$;
6. complete real frames through $n=8$;
7. inverse frames and phase-dressed complex magnitude frames;
8. the endpoints $n=1$, $d=0$, $d=n-2$, and the final depth.

The resource audit checks, using exact arithmetic,

```math
\sum_{q=2}^{n}\frac{2^q}{q}
\leq6\frac{2^n}{n}
```

and the absorption of the polynomial predicate terms into $O(2^n/n)$.

Files:

- [strict-zero construction](../compiler_robust_hopf/strict_zero_echo.py)
- [operator tests](../tests/test_strict_zero_echo.py)
- [exact-rational resource audit](../compiler_robust_hopf/strict_zero_audit.py)
- [resource-audit tests](../tests/test_strict_zero_audit.py)

## 5. Binary–one-hot prefix decoder

The decoder is tested as a complete reversible permutation, not only on its
intended clean input.  The suite verifies:

- $\lvert x\rangle\lvert0\rangle$ maps to $\lvert0\rangle\lvert e_x\rangle\lvert0\rangle$ and returns under the inverse;
- arbitrary computational-basis contents return after forward and inverse;
- every declared layer has disjoint wire support;
- the exact workspace formula
  ```math
  3\,2^t-2-t;
  ```
- the exact gate and depth formulas;
- the one-hot Givens network equals the complete prefix Hopf frame on the code.

Files:

- [decoder schedule](../compiler_robust_hopf/tree_decoder.py)
- [decoder tests](../tests/test_tree_decoder.py)

## 6. Coherent route–operate–unroute

The router is supplied as an explicit schedule.  The tests cover:

- branch-data, token, copy, and reusable-flag register allocation;
- balanced prefix fanout and exact uncopy;
- disjoint Fredkin layers;
- every clean basis input routed to its prefix-selected branch;
- forward route followed by inverse route on arbitrary complex inputs;
- prefix–suffix-entangled inputs;
- token-controlled action of all subtree frames;
- cleanup of every branch flag before inverse routing;
- equality to
  ```math
  \bigoplus_r W_s^{(r)};
  ```
- equality of the complete routed cut to the direct Hopf frame;
- zero probability outside the clean-workspace subspace.

The explicit schedule is cross-checked against

```math
\text{copy wires}=(2^t-1)(s+1)-t,
```

```math
\text{forward Fredkins}=(2^t-1)(s+1).
```

The simulator iterates over branches for convenience.  The circuit-depth claim
uses the declared parallel schedule: branch data, tokens, flags, and subtree
frames have disjoint physical support.

Files:

- [router and sparse-state simulator](../compiler_robust_hopf/router.py)
- [router operator tests](../tests/test_router.py)
- [tree identities](../compiler_robust_hopf/tree_structure.py)

## 7. Resource theorem

The ledgers use integer or exact-rational checks rather than fitted slopes.  They
cover:

- strict zero workspace;
- the one-clean-flag endpoint;
- the direct-to-routed transition;
- maximal feasible cuts over broad $n,m$ grids;
- copy-pool reuse as branch flags;
- the $s=1$ routed endpoint;
- arbitrarily large workspace;
- matching real-state parameter and light-cone lower bounds.

Representative checked implications are

```math
n^2
=O\left(\frac{2^n}{n+m}\right)
\qquad(1\leq m\lt 4n),
```

and

```math
\frac{2^s}{s}
=O\left(1+\frac{2^n}{n+m}\right)
```

for the maximal routed cut, including the saturated endpoint $`s=1`$.

Files:

- [all-workspace resource selection](../compiler_robust_hopf/unified_compiler.py)
- [resource diagnostics](../compiler_robust_hopf/resource_bounds.py)
- [general resource tests](../tests/test_resource_bounds.py)
- [unified compiler tests](../tests/test_unified_compiler.py)

## 8. Gradient records

The decoder suite compares:

- direct record-wise parity averages;
- dense Walsh transforms;
- the fast Walsh–Hadamard implementation;
- direct phase-stream signed one-hot records.

It verifies the deterministic norm-two record property used by the
coordinatewise concentration statement.

Files:

- [decoders](../compiler_robust_hopf/decoders.py)
- [decoder tests](../tests/test_decoders.py)

## 9. Reproduce the checks

Use Python 3.11 or 3.13. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/reviewer_walkthrough.py
python validate.py
python scripts/verify_fault_tolerant.py
```

The last two commands run the deterministic unittest collection and the
four separate exact receipt suites. Resource ledgers and offline provenance:

```bash
python scripts/unified_resource_ledger.py --n 12
python scripts/strict_zero_echo_ledger.py --n 12
python scripts/check_upstream_sync.py --offline
```

The [rendering guide](../assets/README.md#rendering-checks) covers diagrams,
Markdown, MathJax, and native MathML; [detailed coverage](reference/VERIFICATION_CATALOGUE.md#rendering-checks)
records its limits. Local previews do not reproduce GitHub's private renderer.

## 10. Evidence boundary

The repository does not use finite experiments to establish asymptotic
optimality.  It also does not test:

- device connectivity or routing overhead;
- a hardware-native gate set;
- a general elementary-gate emitter for the full asymptotic Clifford+T compiler;
- noisy execution or readout mitigation;
- application-specific controlled-observable implementations;
- arbitrary non-Hopf differential frames.

The proof should be assessed in four separate steps:

1. the Hopf operator identities;
2. the logical circuits and workspace cleanup;
3. the imported synthesis bounds and resource sums;
4. the matched QBP output and access conventions.

## 11. Fault-tolerant evidence and approximate QBP

The [fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md) separates the analytic
resource proof from the [focused finite checks](../verification/fault_tolerant/README.md).
Run `python scripts/verify_fault_tolerant.py` from the repository root. Four
stdlib-only suites reproduce the retained exact Gray-source, shift-kernel,
resource, and reflection receipts without changing their saved reference files.
The runner compares scientific receipt fields and separately verifies the
deterministic source-file hashes, whose values changed when the sources were
relocated into this repository.

| Standalone receipt suite | What it checks | Evidence boundary |
|---|---|---|
| Gray geometric source | Emitted logical source words, actual inverse, scratch return, reference witnesses, and gate-count recurrences | Finite clean-scratch fixtures; the full source is not emitted as an elementary Clifford+T circuit |
| Shift kernel | Complete finite kernel columns, actual inverses, source defects, failure tracking, and negative controls | Small fixed parameters, without the production residual construction or emitted outer amplification |
| Shift resource | Exact rational scales, coefficient rounding, scalar error budgets, and resource ledgers | Analytic proxies and inequalities; the rows are not measured gate counts |
| Geometric reflection | Small source/reflection matrices, exact Pauli-transfer arithmetic, and algebraic norm separation | Supporting source-cost evidence and a kernel-suite dependency; no additive lower bound follows from these fixtures |

The [receipt map](../verification/fault_tolerant/README.md) records the exact
fixture ranges, source dependencies, and provenance for all four suites.

The [approximation bridge](QBP_APPROXIMATION.md) has tests for full-input error,
actual-adjoint transfer, leakage, reference-entangled borrowed inputs, and
bounded estimator bias. These validate the finite examples behind the general
proof, not differentiation of a synthesized family.

The approximation suite also checks the complex phase stream directly,
reflection-term sampling, pathwise finite-weight error, and a two-shot
measurement instrument whose retained dirty state changes the next-shot law.
The last fixture tests the conditional-mean premise of the martingale proof;
it is not a fitted statistical scaling experiment.

---

[← QBP consequence](QBP_CONSEQUENCE.md) · [Complete narrative](../REVIEW.md) · [Source map →](SOURCE_MAP.md)
