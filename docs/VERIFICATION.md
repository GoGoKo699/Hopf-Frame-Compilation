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
| Retained linear-precision depth predecessor | [Complete construction](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision), [shared-gradient recovery](../research/depth/SHARED_GRADIENT_RECOVERY.md), [shared binary terminal source](../research/depth/SHARED_BINARY_TERMINAL.md) | [Bounded all-column recovery](../tests/test_sublinear_depth_recovery.py), [native calibration](../tests/test_phase_gradient_calibration.py), [shared-gradient checks](../tests/test_shared_gradient_recovery.py), [shared binary terminal checks](../tests/test_shared_binary_terminal.py), [native coarse-program checks](../tests/test_grouped_terminal_depth.py), and retained query/indicator suites | Absolute constants in the stated linear-precision range; analytic source, invariant-sector recovery, uniform workspace splice and native schedule; finite checks are not a general compiler |
| Fixed-accuracy polylogarithmic corollary | [Subtree-projector decoder, R78](../research/depth/SUBTREE_PROJECTOR_DECODER.md) | [Exact subtree and bounded spectral-decoder checks](../tests/test_subtree_projector_decoder.py), plus retained source/query/one-clean suites | Analytic full-mode reflection, reversible search, native ledger and width splice; retained as the decoder proof and fixed-accuracy consequence of R79. Finite fixtures do not emit a native descendant loader or prove asymptotics |
| Entry-indexed normalization barrier | [Encoding limit](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md) | Exact Haar spectrum, explicit optimal factorization and nuclear-norm dual witness | Analytic proof for the stated factorization contract; no unrestricted native-depth or recovery-query lower bound |
| Zero-clean five-layer obstruction | [Operator rigidity, R75](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md) | [Exact image orbit, native split words, Clifford boundaries, and common-commutant audits](../tests/test_four_layer_operator_rigidity.py) | Uniform operator-norm theorem, arbitrary dirty width; unrestricted positive-clean five-layer and scalar-gap questions remain open |
| Fixed-clean four-layer rank requirement | [Initialized-sector rigidity, R74](../research/depth/FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector) | [Symplectic, weighted-phase, and initialized-projector audits](../tests/test_four_layer_crossing_rank.py) | Explicit analytic minimum-rank bound $`\Omega_a(\sqrt L)`$; no unrestricted four-layer exclusion with positive clean work |
| Restricted one-clean five-layer obstructions | [Midpoint rigidity and finite projective order, R76](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md) | [Phase descent, native midpoint, compressed-commutant, finite-boundary and cyclotomic checks](../tests/test_one_clean_midpoint_rigidity.py) | Pauli or Hermitian-Clifford encoded midpoint with unrestricted crossing rank; separate finite-projective-order occupied-block obstruction; neither condition follows from terminal return |
| Non-Clifford midpoint obstruction at bounded inner spread | [Inner clean-spread rigidity, R77](../research/depth/ONE_CLEAN_INNER_SPREAD.md) | [Flat expansions, exact sector moments, literal four-layer return, and finite-boundary audits](../tests/test_one_clean_inner_spread.py) | Explicit gap at each spread t; necessary minimum $`\Omega(\log\log L)`$ spread; exact common-flip counterexample has constant occupied defect and does not refute the unrestricted conjecture |

<a id="endpoint-count-evidence"></a>

<a id="two-clean-count-evidence"></a>

### One-clean count evidence

Result C's [one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
and [banked frontier](../research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
give the analytic count and workspace guarantees for both complete-frame
families. Their complete initialized-isometry contract covers every logical
column, literal phase, returned clean work, arbitrary dirty input and
external reference. The [claim ledger](../manuscript/PUBLICATION_SCOPE.md)
states the exact base and doubled-base thresholds, matching ranges and
the separate sharper two-clean complex endpoint.

The all-size conclusions, including the small-dimension baseline, follow
from the proofs. The following fixtures check specific source, recovery
and composition identities; they do not certify an asymptotic compiler.

The real [growing-source diagnostics](../tests/test_growing_precision_source.py) test
normalized refinement, sealed occupied and entangled failures, retained
source tags, attenuation compensation, actual dirty replay, original-anchor
recovery, and linear full-output error. Exact integer checks preserve the
source-scale relation across refinements. The optional requested precision
uses the least power of two $`K\ge L`$ and checks $`q_dM_d=4K`$
for non-dyadic values below and above $`L=N`$, without changing the
logical mode schedule. Those symbolic checks also distinguish incompatible
source precisions and schedules with nonintegral q; they do not certify the
growing-source regime's hypotheses or dispatch its fallbacks. The native arena allocation,
summable charges, and unrestricted count bound come from the proof. These
finite diagnostics are neither a scalable native emitter nor independent
external review.

The same suite separately checks phase dressing of a diagnostic recovered
real factor with nonzero leakage by finite packed-diagonal native words.
The five-wire fixture reuses two flags and two arbitrary dirty-core bits,
checks every logical/dirty input column, literal phases, the actual inverse,
and reference-entangled inputs. Discarded leakage, reversed factor order,
and omitted common phase have counterchecks. Its addressed direct sum is
supplied; it emits neither the growing-source compiler nor a general
diagonal table. The precision and workspace bounds remain analytic.

The [conjugate-recovery suite](../tests/test_conjugate_recovery.py) adds a
separate 16-dimensional diagnostic retaining the original initialized
mode anchor and an arbitrary borrowed signal. It checks a lifted native
interlayer, both signal sectors, dirty/reference inputs, actual inverses,
and a full-operator perturbation that leaves nonzero recovery leakage.
Negative controls incorrectly initialize the borrowed signal or omit an
interlayer lift. Its source is a supplied matrix; it does not emit a
growing-source circuit or establish the all-size banked resource bound.

### Subtree-projector decoder evidence

The [structured common-flag theorem R80](../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md)
has four [bounded algebra checks](../tests/test_homogeneous_depth.py).
Independent guarded-layer and literal-column constructions test general
complex one-qubit factors. Actual native conjugators test the two separate
involution echoes on every small dirty/tag column, including padded tags.
The complete tag word checks exact dirty return and inactive occupied
tags. Exact rational tree products check the retained sign factorization,
including markers below zero incoming edges. These tests distinguish
general unitary geometry from the narrower native factorization used to
approximate real rotations. They do not establish the recursive resource
bound or emit a general compiler.

The [recursive terminal theorem R81](../research/depth/RECURSIVE_TERMINAL_DEPTH.md)
is an analytic extension of retained source/query interfaces to a preserved
prefix, with summable errors and an explicit width splice. Its uniform
terminal-height bound for all $`0\le M\le n`$ is analytic and adds no
separate native primitive or numerical test of asymptotic cutoffs.

The [bounded-density theorem R82](../research/depth/BOUNDED_DENSITY_DEPTH.md)
has six [bounded fixtures](../tests/test_weighted_haar_depth.py): independent
weighted-frame/gauge reconstruction with vanishing corrections; exact
rational subtree variance and Carleson bounds; single-leaf/depth row and
column sources with separate failure padding; a biased product-reference
gauge, weighted energy and overlap check; the complete three-branch LCU
and literal amplification with all-work return and inactive occupied
work; and a perturbed source using its actual inverse with nonzero
leakage retained. The product-reference fixture checks exact reference
geometry, not its native approximation. Source-unitary completions and
bounded matrices are verification fixtures, not supplied native oracles.
Native table synthesis, arbitrary dirty/reference guarantees, recursive
workspace, asymptotic resources, and bias-independent constants follow
from the proof rather than these finite checks. No scalable emitter or
new scientific receipt is supplied; R79 remains the general frontier.

The [affine-orientation theorem R83](../research/depth/AFFINE_ORIENTATION_DEPTH.md)
and [exact deterministic theorem R84](../research/depth/DETERMINISTIC_FRAME_DEPTH.md)
have five [focused checks](../tests/test_orientation_depth.py). Independent
path-product columns are compared with complete packed-tag circuits,
including actual complex native one-qubit words, affine marker signs,
outside-prefix dependence, signed deterministic continuations, arbitrary
dirty inputs, and occupied inactive tags. A guarded affine-map fixture
expands its Toffolis into literal Clifford+T words; a negative control
exhibits the changed-address dirty-update failure. Lookup actions in these
bounded checks are reversible truth-table fixtures. Native query synthesis,
parallel T supports, recursion, literal widths, and asymptotic bounds are
proved in the chapters. R82's affine-reference extension uses these
reference-core identities with its unchanged weighted-source geometry.
No additional scientific receipt or general emitter is introduced.

The [geometric-band orientation theorem R85](../research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md)
has six [bounded checks](../tests/test_geometric_band_orientation.py):
integer bin capacities; independently reconstructed truncated flags and
local bands; the physical dirty-buffer word and its actual inverse;
coherent fine-tag erasure/reload with nonlinear orientations, zero pairs,
signs, leaf phases, and a dirty reference; five-bin coarse search/XOR/reverse
with padded thresholds; and inactive fine-word identity with occupied
protected records. The tests use ideal matrices and exact basis
permutations. The coarse check reuses R78's literal spectral-search word.
The complete native schedule, source-error composition, literal widths,
and asymptotic resources remain analytic. These fixtures emit no general
native compiler and introduce no scientific receipt.

The [in-place permutation checkpoint](../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
has two [boundary fixtures](../tests/test_permutation_bridge_limits.py).
One exhausts all permitted clean/dirty inputs of the four-shear order-eight
example and checks the order-sixteen triangular witness. The other builds
explicit binary-linear Clifford permutations and dense spectral isometries
to check why an exact encoding-order restriction does not extend to fixed
error. Neither fixture supplies a native encoder, a general conversion,
or an asymptotic proof. The checkpoint records its scoped restrictions
and the frontier at that stage; scientific receipts are unchanged.

The [precision-uniform theorem R79](../research/depth/UNIFORM_SUBTREE_DEPTH.md)
uses the unchanged R78 decoder identities below. Its new work is analytic:
retain the original source's explicit precision factors, enlarge the
protected bank, prove one absolute workspace cutoff, and combine R72's
one-clean middle stream, logarithmic terminal and remaining precision
branches. No new circuit identity, imported primitive or native fixture
is introduced. Existing source, terminal and query suites support their
component contracts; finite checks do not establish the new all-precision
resource envelope or certify its absolute constants.

The [fixed-accuracy theorem R78](../research/depth/SUBTREE_PROJECTOR_DECODER.md)
has five [focused checks](../tests/test_subtree_projector_decoder.py).
Exact rational path products independently reconstruct signed, unequal and
singular frame columns and verify their nested projectors and marker bands.
Small physical Hadamard/reflection words check coherent search, padded
thresholds, every occupied old tag, and literal phases. An entire small
occupied-work space checks inactive cancellation and actual inverse use.
Leaked source columns with deliberately different unitary completions test
the full-space projector bound and the complete decoder error, retaining
nonzero returned-work leakage. Missing signs, incorrect padding, an unguarded
tag XOR, and a substituted inverse have explicit counterchecks.

These are finite algebraic fixtures. The source/query native schedules,
all-size polynomial absorption, literal one-clean allocation, terminal
composition, and width splice are proved analytically in R78. The suite
emits neither a native descendant loader nor a scalable frame compiler.

### Sublinear-depth evidence

The [compressed-column theorem](../research/depth/SUBLINEAR_T_DEPTH.md)
has separate source, exact dirty-query, complete-recovery and tail proofs.
Its five [bounded recovery checks](../tests/test_sublinear_depth_recovery.py)
use independently written two-target descendant vectors and small numerical
unitary completions. They check every frame column for unequal signed and
singular tuples, marker packing, padded-tag overlap, actual inverse use,
active-only amplification phase, and residual work/flag leakage. Negative
controls preserve the root while corrupting a marker, substitute a different
completion's inverse, or apply the active minus globally. The last error
check includes final guard erasure and exact identity on inactive columns.

Five [binary-gradient calibration checks](../tests/test_phase_gradient_calibration.py)
reuse the emitted helper-free TTK adder. Exhaustive small register inputs
cover occupied programs, gradient offsets, both target signs, inactivity,
and actual inverses. A two-bit Clifford+T matrix checks literal phases;
larger small kickback blocks check the real rotation. A coarse native
product source produces nonzero leakage and verifies its complete-output
error bound. Negative controls drop a guard or target sign, or substitute
an ideal preparation inverse. The source is deliberately coarse, not an
implementation of arbitrary-accuracy one-qubit synthesis.

The [shared-gradient checks](../tests/test_shared_gradient_recovery.py)
add a literal two-bit compute/swap/uncompute adder, occupied-work inactive
echoes, and a complete coefficient-bank phase table. Recovery-shaped
noncommuting words test one native preparation and its actual inverse at
the two boundaries, including nonzero final leakage and a preparation-error
bound independent of the number of middle calls. Counterchecks reject a
reflection that tests the live gradient, a substituted ideal inverse, and
a discarded coefficient-bank scalar phase. Programmable angle accuracy
requires a zero program; the fixed outer damping instead ignores program
contents. These checks use small exact permutations and matrices of
dimension at most 64. They test the local identities and source boundary
comparison, not a complete emitted amplitude-amplification compiler.

The [grouped-terminal checks](../tests/test_grouped_terminal_depth.py)
combine literal native H/T/T-dagger interpretation with complete-word
prefetch and actual unloading. They test occupied inactive programs,
predicates and dirty helpers, with negative controls for an omitted phase
correction, a mutable query address and a missing activity guard. Integer
checks cover rounded group sizes, execution order and the source/coarse
balance. Their small column batches and illustrative allocation constants
do not establish a scalable emitter or a numerical crossover.

The [shared binary terminal checks](../tests/test_shared_binary_terminal.py)
add an exact sparse two-stage, sixteen-type accepted product on every
source/logical column. Noncommuting coarse words produce coherent prefix
addresses for the next stage; failure amplitudes and literal complex phases
are retained. The source is unprepared only at the final boundary. Missing
history and intermediate source projection give different accepted blocks.
This algebra fixture supplies coefficient-table actions and folds the local
suffix predicate into atom rejection; it does not emit the full query or
production-precision circuit. Separate native shift/history words check all
occupied register inputs, dirty-helper return, literal phases and actual
inverses. A native source/body/inverse word checks inactive cancellation on
occupied coherent inputs. Wrong inverse ordering, a raw disabled source,
and preparation substituted for its inverse have negative controls. The
pinned generic source/shift/history receipts remain unchanged.

These fixtures do not emit the derived Boolean tables or a complete native
descendant loader. The existing native dirty-query, shared-control phase,
returned-helper predicate and nonuniform-indicator tests cover their stated
finite interfaces. Terminating source selection, polynomial and arithmetic-table overhead
absorption, and the same-circuit asymptotic theorem are proved analytically;
no numerical extrapolation or efficient preprocessing claim is made.

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
