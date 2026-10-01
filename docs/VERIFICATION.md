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

The [operator-source proof](OPERATOR_SOURCE_COMPILER.md) is accompanied by
[finite circuit checks](../tests/test_operator_source_compiler.py).
These test the native geometric operator, its programmable scalar block,
the two-flag rotation, coherent amplification, dirty word-bank routing, the
restored suffix echo, and literal phase-diagonal composition on arbitrary dirty
inputs.
The analytic proof establishes the improved upper bound; finite matrices test
its fragile identities and conventions.

The suite now includes a complete two-layer **native Clifford+T** fixture
with nine wires and exactly two initialized flags. All 128 logical/dirty
input columns are propagated through the actual gate words, including an
exact seven-T Toffoli decomposition, the literal amplification phase,
actual reversed-word inverses, and retained intermediate leakage. Its
deliberately coarse precision tests composition rather than claiming the
theorem's asymptotic accuracy. Separate fixtures verify general U(2)
multiplexor Euler order, address-dependent phases, and four-stage composition,
as well as the optimal exact source words and Pauli-transfer witnesses.

The [source-reuse limits](SOURCE_REUSE_LIMITS.md) have
[separate finite checks](../tests/test_source_reuse_limits.py) for nilpotent
contractions, arbitrary encoding bases, dirty-dimension independence, and
transformed-mask operator identities. Counterexamples test why nilpotence
and complete-output accuracy cannot be omitted. These support the scoped
analytic restrictions; they do not establish a full-frame impossibility.

The [conditional-suffix proof](CONDITIONAL_SUFFIX_COMPILER.md) has
[focused grouped-block checks](../tests/test_conditional_suffix_compiler.py).
They combine a native small operator source with separate scalar and
atom flags, retain every dirty input column through amplification,
and check literal phases, the disjoint ancestor-column support partition,
forward and reverse column blocks, inactive suffix sectors, and predicate
uncomputation after leakage. A separate integer ledger tests exponentially
growing groups and the dirty-width/error inequalities. Its sample
grouping constants do not certify the unspecified fixed constants of the
coarse synthesis primitive. These are finite interface checks, not a native
elementary circuit emitter for the whole asymptotic grouped compiler.

The [T-depth schedule](T_DEPTH_COMPILER.md) has
[native routing checks](../tests/test_t_depth.py) for literal four-layer
shared-control Fredkins, disjoint T-layer supports, actual inverses, and
whole dirty-bank queries. Clifford layers may have substantial depth;
the fixture does not treat T-depth as total execution depth.
The [parallel dirty-lookup checks](../tests/test_parallel_dirty_lookup.py)
audit exact bilinear cancellation, literal native phases, disjoint T layers,
recursive scratch reuse, and arbitrary-input return using symbolic Boolean
polynomials. Their [analytic composition](PARALLEL_DIRTY_LOOKUP.md) retains
the count bound while reducing T-depth under its sufficient dirty-width
condition. The linear table maps still have a charged Clifford-depth cost.
The [tree-residual checks](../tests/test_tree_residual_structure.py)
reconstruct complete small residuals from classical tree generators,
including complex coarse words and singular angles.
The [tree-transport checks](../tests/test_tree_transport.py) also reconstruct
the sparse resolvent, verify the exact transport Gram matrix and complete
three-mode unitary columns, and test subtree overlap-defect telescoping,
weighted norms, and omitted perturbation terms. These support the
[endpoint candidate analysis](ENDPOINT_TREE_TRANSPORT.md); they do not
construct a cheaper native joint block.
The [weighted-block checks](../tests/test_weighted_transport_block.py)
verify its nonorthogonal-column witness, scalar recursion against independent
subtree Schur solves, complete one-flag unitaries and actual inverses, and
weighted-map perturbation bounds at zero and near-zero defects. They also
check the batched allocation, literal level packing and gathering, the
physical marker permutation, and the selected completion's discontinuity.
The [batched native bound](WEIGHTED_TRANSPORT_BLOCK.md#5-a-native-implementation-by-depth-batching)
is proved analytically. These operator matrices have dimension at most 32;
the separate permutation checks enumerate basis labels without dense matrices.
They do not emit the elementary circuit or establish the linear T-count target.
The [residual-assembly checks](../tests/test_residual_assembly.py) compare the
affine forward recursion with independent subtree Schur solves, reconstruct
the complete affine and reverse dilations, and assemble and amplify their
two-flag selection on every logical input column. They check literal mode
packing and gathering, singular and zero-defect cases, actual inverses,
and perturbations that retain rejected-flag leakage and arbitrary dirty
inputs. Negative controls detect the wrong reverse branch, a changed
relative root phase, and replacement of an actual inverse by a forward call.
The [native assembly and resource bound](RESIDUAL_ASSEMBLY.md) are proved
analytically. These operator fixtures use matrices of dimension at most 64;
they do not emit the native compiler or demonstrate an improved endpoint
T-count.
The [affine-tree fusion checks](../tests/test_affine_tree_fusion.py) compare
the ten-mode and recursive full-input merges with literal local products,
actual inverses, and independent path maps. They test continuation ranks,
growing column support, eager entry counts, and a bottom-only perturbation
whose internal transport is undamped despite a strict external norm margin.
Native paired-source fixtures separately check transformed-mask correlations,
changed-address uncomputation, and rejected-space return when a scalar signal
is reused. Matrices have dimension at most 64. These support the
[fusion audit](RESIDUAL_ASSEMBLY.md#7-a-bounded-audit-of-fusion-across-tree-depths)
and [scoped source-reuse arguments](SOURCE_REUSE_LIMITS.md); the mask
T-count bound is analytic, and neither the fixtures nor the representation
witnesses establish a general gate lower bound or a cheaper native compiler.
The [coupled-merge checks](../tests/test_coupled_residual_merge.py) compare
the complete normalization-two recursion and its unfolded target wrappers
at a fork and height three, with close complex native coarse words and real
targets. They retain every signal port, actual inverse, and correlated
spectator input; they test the stability estimate, anchored mixed term,
and the complete rank-four repair. These small operator fixtures do not
price the transported repair modes or emit a faster native compiler.
The [transported-repair checks](../tests/test_transported_repair.py) test
the subsequent commutator circuit using only the supplied child and local
parent words. They cover intersecting and vanishing transported modes,
the weighted repair-error bound, exact cancellation with noncanonical
actual words on shared dirty work, and the failure of an unmatched ideal
inverse. The repaired full merge retains the local half-error bound;
these tests do not establish shared precision or a better T-count.
The [shared-conjugator checks](../tests/test_shared_conjugator_merge.py)
use literal Clifford+T source and routing words on all 128 core/logical/signal
columns. They test valid outer cancellation, the surviving two fork returns,
a direct half-unitary target failure, and the precision-independent
coefficient-ellipse obstruction to retuning the same word. A second flag
is an untouched spectator. The literal source ledger concerns that emitted
word; it is not a general gate lower bound.
The [antichain checks](../tests/test_antichain_compiler.py) compare the exact
strict-descendant forest factorization with complete complex tree words,
pack mixed-depth disjoint updates into one last-bit multiplexor, and check
native dirty-Fredkin echoes on every input, including inactive sectors.
They also test reflection-by-reflection predicate echoes and detect the
failure of the stated factorization when updates are comparable or the
forest is omitted. Matrices have dimension at most 64. The
[antichain compiler](ANTICHAIN_COMPILER.md) proves the native resource and
full-operator error bounds analytically; these fixtures neither emit its
fine-precision synthesis circuit nor establish the unrestricted endpoint.
The [sparse-update checks](../tests/test_sparse_update_compiler.py) verify
ancestor-closed support, exact off-support forest factorization, affine
basis-state transpositions, and the packed operator's identity complement.
They reconstruct the dense column dictionary and detect reversed atom/filter
order and an incorrectly shared rejection flag. Small component matrices
and batched input columns test actual inverses, complete amplification,
inactive-sector identity, and dirty-core perturbations including rejection.
The [sparse-update proof](SPARSE_UPDATE_COMPILER.md) supplies the native
resource and conditional-workspace bounds. These checks are finite interface
tests, not an elementary emitter for the precision-dependent compiler.
The [one-clean checks](../tests/test_one_clean_compiler.py) reconstruct the
paired-Majorana source and general Pauli masks from native gates, audit the
conjugated scalar word and five-call amplification on all dirty input
columns, and check addressed relative phases and exact inactive action.
They also test X symmetry, the full-operator borrowed-signal estimate,
reference stability, and a scalar-phase counterexample to that extension.
The fine-precision checks use an independent small Clifford-algebra
representation, rather than a large precision-core simulation. The
[one-clean theorem](ONE_CLEAN_COMPILER.md) supplies the analytic error and
workspace proof; these fixtures do not emit its asymptotic grouped circuit.
The [source-merge checks](../tests/test_source_merge.py) use native scalar
sources and noncommuting three-level logical operations to exhibit the
Pauli-routed cubic return on every dirty input. These support the scoped
claims in [source-reuse limits](SOURCE_REUSE_LIMITS.md); neither classical
compression nor a failed merge fixture settles the linear frame endpoint.
The direct single-flag merge is also checked against its actual native
branch word, including its constant error and identity-angle amplification
failure.

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

Use Python 3.11 or 3.13.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the short orientation:

```bash
python scripts/reviewer_walkthrough.py
```

Run the deterministic unittest collection and the separate exact
fault-tolerant receipt suites:

```bash
python validate.py
python scripts/verify_fault_tolerant.py
```

The first command discovers the files in `tests/`; the second reproduces the
four retained receipts in `verification/fault_tolerant/`.

Print the two resource ledgers:

```bash
python scripts/unified_resource_ledger.py --n 12
python scripts/strict_zero_echo_ledger.py --n 12
```

Check the recorded upstream versions without network access:

```bash
python scripts/check_upstream_sync.py --offline
```

For the optional browser audit, follow the
[rendering guide](../assets/README.md#rendering-checks). It checks all diagrams
and complete Markdown pages, including table mathematics, display equations,
missing expressions, and horizontal overflow. Equations are checked both with
MathJax SVG and the browser's native MathML layout. The native checks reject
unsupported numbered rows and detect the narrow vertical stacking that can
otherwise pass expression-count and overflow checks. Equation numbers are
printed as ordinary mathematical text so native rendering preserves them.
Comparison signs use HTML-safe TeX commands, guarded against literal
less-than characters that GitHub can misinterpret before math rendering.
Saved desktop and narrow-screen previews support human inspection. These
local checks cover both rendering paths but do not reproduce GitHub's private
renderer or guarantee identical layout across all browsers.

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
