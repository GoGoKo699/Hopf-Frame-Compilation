# Documentation map

[Landing page](../README.md) · [Read the argument](../REVIEW.md) · [Publication scope](../manuscript/PUBLICATION_SCOPE.md)

The selected compiler package has four principal results: exact compilation,
sufficient-clean matching T-count, one-clean compilation, and simultaneous
T-count/T-depth tradeoffs for complete real frames. Start with the landing
page and `REVIEW.md`, then follow the relevant claim below.

## Proof chapters

### Shared target and error contract

State preparation fixes one column. The compiler target fixes the state and
prescribed differential-frame columns, with the declared workspace returned.
Approximate compilation includes all logical inputs, arbitrary dirty inputs
and their references, literal phases, and clean-work leakage.

| Read | Purpose |
|---|---|
| [Hopf interface](HOPF_INTERFACE.md) | Coordinates, marker columns, singular charts, addressed layers, and the full-frame tangent metric |
| [Frame-safe compilation](FRAME_SAFE_COMPILATION.md) | Complete-input error, actual adjoints, fixed-decoder necessity, and singular exceptions |
| [Compiler boundaries](COMPILER_BOUNDARIES.md) | Concrete distinctions between state preparation, checkpoints, and prescribed complete frames |

### A. Exact compilation at every clean-workspace budget

The [exact compiler theorem](COMPILER_THEOREM.md) proves matching elementary
size and depth, and a separate CNOT-count bound, for every number of clean
work qubits. It covers real and phase-dressed complex magnitude frames in
the arbitrary-one-qubit/CNOT model.

Read its echo, decoder, router, and maximal-cut constructions together with
the lower bounds. The proof returns workspace exactly and preserves every
required column; a state-preparation theorem alone does not establish that
contract. The exact gate model is distinct from the Clifford+T model below.

### B. Matching T-count with sufficient clean workspace

The [fault-tolerant compiler theorem](FAULT_TOLERANT_COMPILER.md) gives the
matching worst-case precision/workspace frontier under its sufficient clean
reservation. Real and phase-dressed complex magnitude families have their
stated extensions. Follow the complete source, residual dictionaries,
compressed failure history, amplification, and resource ledger.

Section 10 gives the real-frame reduction to diagonal synthesis, the
fixed-width counting argument, and an explicit single-angle precision
witness valid at arbitrary width. These lower bounds also support C and D;
they are not obtained by adding costs of particular source calls.

### C. One-clean compilation and its corollaries

The [one-clean compiler](ONE_CLEAN_COMPILER.md) and
[conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) establish the
grouped real-frame bound and its banked refinement. The latter's Section 10
extends the complete group construction to one initialized compiler flag.
Conditional logical zeros are usable only on the active sector; the proofs
also prove inactive identity on arbitrary inputs and the full return-error estimate.

| Supporting proof | Why it remains in the package |
|---|---|
| [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) | Full-input two-clean source, packed diagonal precision/workspace bound, amplification/error certificate, layerwise frame baseline, and separate multiplexor reservation |
| [Conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) | Residual columns, active logical work, iterated grouping, and banked one-clean composition |
| [One-clean compiler](ONE_CLEAN_COMPILER.md) | Native one-clean primitive; literal diagonals, complete one-target U(2) multiplexors, phase-dressed magnitude frames, and the separate zero-clean real-frame corollary |
| [Borrowed-workspace compiler](BORROWED_WORKSPACE_COMPILER.md) | Arbitrary clean/dirty allocations, exact dirty traversal, and the restricted matching splice |

Keep each corollary's literal dirty reservation. Reducing the clean count
does not supersede every earlier resource point. The complex magnitude
extension is a real frame followed by supplied leaf phases, not arbitrary
complex-unitary synthesis.

### D. Simultaneous real-frame T-count and T-depth

Start with [uniform precision and depth](UNIFORM_PRECISION_DEPTH.md).
Its theorem gives absolute-constant bounds on one circuit, with two clean
flags and its literal dirty threshold. It proves simultaneous matching
count and depth on explicit intervals, plus a stronger slowly growing
precision range. At fixed accuracy and sufficient square-root dirty width,
optimal-order T-count accompanies an O(n) T-depth upper bound.

The proof has an early unary-group branch, an exact rectangular-query tail,
and an all-precision fallback. The following chapters are mathematical
dependencies, not a sequence of research tasks to repeat.

| Part of the argument | Proof dependency and role |
|---|---|
| Uniform composition | [Uniform precision and depth](UNIFORM_PRECISION_DEPTH.md): rectangular allocation, uniform cutoff, weighted sums, literal-width fallback, and matching inequalities |
| Selected table queries | [Blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md): controlled bilinear helper echo, selected-block dirty traversal, four-call cancellation, and fixed-accuracy composition |
| Charged early source | [Unary phase-source groups](UNARY_PHASE_GRADIENT.md): exact selected cyclic shifts, conditional work, source preparation/actual inverse, and complete group error |
| Program and selector lifetime | [Grouped program prefetch](GROUPED_PROGRAM_PREFETCH.md), especially Section 10: preserved programs, incremental prefix selectors, and consume-before-change suffix enables |
| Geometric predecessor | [Conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md): the complete conditional-source interface used by the retained geometric grouped construction |
| Tail indicators | [Chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md): full-input dirty-tree echo, live-work bound, and the summable chunk allocation |
| Shallow dirty arithmetic | [Masked dirty-sum compression](DIRTY_SUM_COMPRESSION.md): compressor invariant, weighted-count recurrence, stable-prefix pipeline, and exact offset cancellation |
| Literal-width fallback | [Amortized dirty lookup](AMORTIZED_DIRTY_LOOKUP.md): exact selected-shear query, capped source precision, and the same-circuit resource/error ledger |
| General-precision hybrid | [Parallel dirty lookup](PARALLEL_DIRTY_LOOKUP.md): bilinear queries, dirty-counter primitives, and the hybrid used outside the low-precision branch |
| Native scheduling and depth bounds | [T-depth schedule](T_DEPTH_COMPILER.md): exact shared-control phase/Fredkin batches, routing, two-dirty predicates, and count-to-depth lower-bound transfer |
| Angular approximation | [Hopf error accumulation](HOPF_ERROR_ACCUMULATION.md), Section 1: full-frame square-sum bound for ideal angle perturbations |

Older constructions remain where they supply a primitive, error certificate,
or fallback with a distinct reservation. The masked compressor and its count
recurrence are needed by the later pipeline; the grouped selector schedule
is needed by the unary source. Their proofs must not be replaced by a link
to a stronger headline bound.

D is a real complete-frame theorem. The complex extensions of A–C and the
state-based protocol have separate contracts. T-depth allows arbitrary
Clifford interlayers; their gate counts are charged, and their physical
depth is not declared zero. The large-width depth optimum remains open.

### Fixed-decoder QBP consequence

The [QBP consequence](QBP_CONSEQUENCE.md) explains why preserving the
additional columns matters for the original inverse-frame decoder.
[Approximate QBP](QBP_APPROXIMATION.md) proves fixed-parameter robustness,
complete complex gradients, observable-sum and classical-weight errors,
and concentration with correlated dirty-bank reuse.

These statements separately charge quantum executions, observable access,
initialized protocol work, and classical output. They do not establish
optimality among all gradient algorithms or an end-to-end advantage.

### Separate supplements and optional research

The completed [state-based QBP supplement](../supplements/state_based_qbp/README.md)
uses a changed decoder and an initialized state-preparation interface. Its
real/complex theorem, computational input audit, task costs, and bounded
native integrations are grouped there. It does not replace the complete-frame
contract or become an extra prerequisite for A–D.

The [research index](../research/README.md) maps constructive endpoint
reductions, promised-family compilers and scoped interface bounds. The
[bounded-factor bridge](../research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
uses the packed diagonal compiler; the
[regular tree Cayley reduction](../research/endpoint/TREE_CAYLEY_REDUCTION.md)
supplies an all-angle representation with charged exterior permutations.
Their remaining hypotheses are separate from Results A–D.

## Evidence and sources

| Read | Purpose |
|---|---|
| [Core-claim audit](CORE_CLAIM_AUDIT.md) | Recorded internal claim-to-proof review and its evidence limits |
| [Verification](VERIFICATION.md) | Analytic proofs, exact certificates, native finite checks, and missing general-emitter coverage |
| [Fault-tolerant receipts](../verification/fault_tolerant/README.md) | Four standalone exact source/kernel/resource checks |
| [Source map](SOURCE_MAP.md) | Imported premises, inherited interfaces, and local constructions |
| [Related work](RELATED_WORK.md) | Primary-source comparisons by target, gate model, precision, and workspace |
| [Provenance](../provenance/README.md) | Source versions and upstream lineage |

Finite checks test their specified contracts; they do not prove asymptotic
theorems, establish priority, or supply independent peer review. No general
elementary Clifford+T emitter for the asymptotic grouped compiler is claimed.

The [implementation](../compiler_robust_hopf/README.md),
[tests](../tests/README.md), [scripts](../scripts/README.md), and
[diagrams](../assets/README.md) provide reproduction maps.
The [manuscript guide](../manuscript/README.md) gives the selected argument's
section order.
