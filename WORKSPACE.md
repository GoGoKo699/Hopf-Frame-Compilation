# Continuing workspace

<a id="mandate"></a>

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
are authorized. The destination is optimal general native T-depth for complete
Hopf frames, with simultaneous count, width, precision and work return.
Manuscript drafting and the changed-decoder supplement are separate tasks.

This checkpoint owns the current question and next action. Use the proof homes for derivations and the route register for reopening conditions.

<a id="current-checkpoint-native-source-scope-and-complete-readout"></a>

## Current checkpoint

The verified scientific baseline is the [two-angle checkpoint](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/287)
at [4d018bd](https://github.com/GoGoKo699/Hopf-Frame-Compilation/commit/4d018bd34c262fc43f7e14d270f68e64e882d193).
The [shared-layer construction](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/286)
enlarges the allowed source budget. The two-angle reduction supplies one sufficient route.
The continuation prioritizes a joint fixed-error construction, keeps a bounded phase attempt,
and investigates growing lower bounds independently. This plan changes no proved frontier.

<a id="current-depth-frontier"></a>

## Destination and strongest established statements

At fixed precision $`L\ge6`$, with one clean input and $`b=\Theta(\sqrt N)`$ dirty qubits,
general complete frames have a [five-band lower bound](research/depth/ROBUST_DECISION_BIAS.md#complete-return-four-band-floor)
and an $`O(\log n)`$ native T-depth upper bound, where $`N=2^n`$.
[R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md) gives the all-precision upper with optimal-order T-count and
$`G=O(NL)`$ at its stated widths. [R79](manuscript/PUBLICATION_SCOPE.md) remains selected Result D.
The [resource frontier](docs/OPEN_PROBLEM.md#large-width-t-depth) gives the one-clean precision refinements,
the matching precision-dominated lower bound and the separate larger-clean tradeoff.

<a id="main-unfinished-route-and-exact-missing-operation"></a>

## Research plan: joint fixed-error construction first

The first constructive milestone is general fixed-accuracy depth $`o_L(\log n)`$
at one clean input and square-root dirty width, retaining $`T=O_L(\sqrt N)`$,
$`G=O_L(N)`$ and the complete return contract. Optimal depth remains the destination;
a growing lower bound is an independent way to narrow that gap.
The precision and larger-clean results retain their separate parameter regimes.

### 1. Main construction: control the whole bank error

Start with one joint complete-indicator candidate at fixed operator tolerance.
Use the [source contract](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#fixed-singly-exponential-source-overhead):
zero clean inputs, sublogarithmic depth, and fixed singly exponential resources
with a common cost exponent across fixed tolerances. One suitable fixed tolerance
also suffices under the completed-address hypotheses in the interface table below.

**First feasibility question:** can the errors of completed Boolean rows occupy
orthogonal work sectors, or cancel coherently, while their joint native word has
a uniform bank-error bound? Test completed rows and their actual inverses on
arbitrary dense banks and superpositions. Sector selection, transport and return
must have explicit native costs. A bound on the fraction of bad dirty inputs
does not establish a small operator error.

The [row-product bounds](research/depth/COHERENT_COUNTER_SOURCE.md#a-finite-angle-precision-target-suffices)
already constrain independent rows and ordered subset products on shared work.
Identify which restrictive premise changes, or prove that the shared-work
necessary condition is compatible with a uniform operator-norm bound.
Merely sharing a register or allocating nonuniform row precision does not establish an escape.
Stop a candidate when its complete synthesis restores the old logarithmic
lookup, carry or conjunction cost. If the complete-indicator requirement is the
obstacle, evaluate the [one-clean coherent decision](research/depth/PRESCRIBED_SCALAR_TABLES.md#one-clean-decision-source)
instead: prove its sign and uniform gap on the entire initialized code before recovery.

### 2. Bounded construction: test the public third-root phase

**Can the two fixed phase operations below be synthesized to error $`2^{-s}`$ in $`o(\log(s+2))`$ native T-depth?**
The angles are $`2\pi/3`$ and $`\arccos(1/\sqrt3)`$, with their negatives obtained by actual reversal.
Allow $`2^{O(s)}`$ T gates, elementary Cliffords and dirty workspace, and zero clean inputs.
The third-root word uses its public bit only as a CNOT control, so address and eigenlabel parities can drive simultaneous calls.
The second operation is an ordinary phase rotation on a private qubit.
Each primitive may carry one common phase. The row cancels it with actual inverses. Bound the full operator error on arbitrary dirty inputs. [Phase word and eigenlabels](research/depth/COHERENT_COUNTER_SOURCE.md#a-finite-angle-precision-target-suffices)

| Step toward complete frames | What the proof pays |
| --- | --- |
| Phases to one MOD3 row | The existing [dirty gauge](research/depth/COHERENT_COUNTER_SOURCE.md#private-gauge-mod3) returns the whole register to the requested tolerance and costs $`O(\log(s+2))`$ depth. Independent rows need a joint bank error bound. For a fixed bank tolerance $`\tau`$, row error $`O(\tau 2^{-p})`$ suffices on a p-bit address, so individual phase precision has $`s=p+O(\log(p+2)+\log(1/\tau))`$. |
| Rows to a complete indicator | An [odd MOD3 residue kernel](research/depth/COHERENT_COUNTER_SOURCE.md#odd-parity-row-kernel) is its own inverse under binary convolution. A paid CNOT basis change, the parallel row phases and their recovery toggle the selected bit of an arbitrary bank. |
| Indicator to general frames | [Shared-layer fusion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#fixed-singly-exponential-source-overhead) absorbs singly exponential source costs with a common exponent across fixed tolerances. The [complete word-query transfer](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#7-conditional-complete-frame-consequence-with-its-source-tolerance) then gives $`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$ and $`D_T=O_L(N/b^2)+o_L(\log n)`$ with one clean input at the original real and complex width thresholds. |

A new phase word should identify how it pays the non-Clifford controls in parallel and bound the full-space error. A cyclic gauge must also account for wraparound. The [component bound](research/depth/ROBUST_DECISION_BIAS.md#the-same-components-constrain-literal-phase-precision) applies to the native product of the public word specialized at one and its adjoint specialized at zero. Their common phase cancels, giving a zero-clean approximation to $`e^{2\pi i/3}I`$. For that product, $`d(r+1)=\Omega(s)`$ for depth d and maximum interaction-component half-rank r. A joint bank construction can instead target the complete indicator norm directly.

Begin with the public word, which has the tighter control convention. Require
an all-size native family, a full spectral or operator-error bound, and a
precision-depth recurrence that beats $`O(\log s)`$. Its coupled components must
satisfy the rank requirement above. Stop a family that only reconstructs the
existing increment-and-gradient gauge or imports an initialized precision state.
The private-angle word is the next obligation only after the public mechanism
passes this feasibility test. Failure here closes that candidate, not the frame problem.

### 3. Independent lower bound: iterate complete-return constraints

Start from the [returned-query probe identity](research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md#4-what-a-growing-lower-bound-would-now-need).
Try a one-band extension rule for significant complete-return probe descriptions,
with quantitative dependence on depth, fixed error and the simultaneous count budget.
This is a proposed induction task, not an established invariant.
It must allow arbitrary Clifford walls and large connected interaction components.
Complete return holds for the finished circuit. Keep both the forward prefix and
actual-adjoint suffix at each cut; do not assume that an intermediate prefix returns work.

Before extending the rule, test it against [exact cancellation at large crossing rank](research/depth/ROBUST_DECISION_BIAS.md#exact-cancellation-at-large-crossing-rank)
including the [single dense connected identity](research/depth/TEMPORAL_AXIS_BIAS.md#the-anticommutation-graph-and-its-temporal-cuts),
and [constant-depth Gauss basis bias](research/depth/GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions).
Raw rank, basis-only bias or a fixed source grammar cannot supply the needed bound.
An $`O(T^2)`$ description bound alone is insufficient at the allowed square-root
T-count; a bounded cut or component rank would restrict the circuit model.
Retain a new constant-band exclusion as a scaling route only if it exposes an
iterable argument; otherwise record its scoped use in screening a live candidate.
The [precision tradeoff](research/depth/PRECISION_COUNT_DEPTH_TRADEOFF.md) alone gives
no depth lower bound growing with n at fixed L and the allowed square-root T-count.

### Work cycle and decision points

1. Write one candidate's full operator, register allocation and proposed scaling mechanism. Check its reopening condition in the route register.
2. Derive an all-size identity or resource/error recurrence. For a recursive proposal, inspect two symbolic levels before assuming it iterates. Use small exact fixtures to test specific identities or counterexamples.
3. Audit the complete native costs, actual inverse and arbitrary dirty/reference action. A successful component must pass the existing full-frame transfer at its actual requested tolerance.
4. Choose explicitly: continue a qualifying mechanism, revise a concrete failed premise, or retire the candidate. Keep the joint construction and lower-bound investigations independent of the phase outcome.
5. Integrate only established results and decision-relevant failures in their existing homes, run the required checks, and replace this checkpoint with the next unpaid operation.

Judge progress by an improved unconditional bound, a qualifying scalable native
source, an iterable lower-bound mechanism, or a decisive test of a live candidate.
Another sufficient interface alone does not establish that the final gap has narrowed.

<a id="the-two-useful-source-interfaces"></a>
<a id="other-retained-results-and-alternatives"></a>

### Choose a complete source interface

| Alternative entry point | First operation to establish |
| --- | --- |
| [Complete indicator](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#shared-layer-indicator-fusion) | Toggle the addressed bank bit with zero clean inputs at every fixed full-space tolerance in sublogarithmic native depth and fixed singly exponential resources. One fixed error below $`1/\sqrt2`$ also suffices under [exact completed address preservation](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#fixed-tolerance-shared-layer-source) or [summable address motion](research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#fixed-tolerance-summable-address-source). |
| [One-clean scalar decision](research/depth/PRESCRIBED_SCALAR_TABLES.md#square-root-width-signed-gap-transfer) | On the entire initialized code, $`M=J^\dagger V^\dagger Z_c VJ`$ needs sign $`D_f\otimes I`$ and a uniform positive gap for any prescribed table f on a bits. Raw garbage is allowed. With $`T,b_{\rm dirty}=O(2^{a/2})`$ and $`G=O(2^a)`$, depth $`O(\log a/\log\log a)`$ suffices at $`b=\lceil K_L\sqrt N\rceil`$. The transfer pays for amplification, promotion and reversal. |
| [Direct positive-frame construction](research/depth/BOUNDED_DENSITY_DEPTH.md#bounded-density-scalar-transfer) | A paid reference and target-to-reference density in a fixed positive interval give constant scalar-call multiplicity. Extending to general frames must accommodate the [embedded Boolean table](research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md#near-uniform-positive-frames), already present for strictly positive near-uniform targets. |
| [Lower bound from complete return](research/depth/ROBUST_DECISION_BIAS.md#native-interaction-spectral-volume) | On the one-clean rare cycle, $`(d+1)(r+1)=\Omega(\log(p+2))`$. To prove a growing unrestricted lower bound, handle large connected components. The [zero-clean six-band floor](research/depth/ROBUST_DECISION_BIAS.md#zero-clean-six-band-floor) uses a different allocation. |

## Stopped routes, reopening conditions and superseded instructions

Look up the proposed operation in [route decisions](research/ROUTE_DECISIONS.md#find-the-decision-for-the-proposed-operation) before revisiting it.
In particular, check the [separated-block margin](research/depth/PRESCRIBED_SCALAR_TABLES.md#two-address-block-margin) for split-address sources,
the [Gauss example](research/depth/GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions) for basis-only decisions and
the [sparse-bank dimension bound](research/depth/GATEWISE_READ_ONLY_EQUALITY.md#sparse-code-preprocessing) for preprocessing.
The register records the excluded class and what would reopen it. Historical session prompts describe earlier work.

<a id="research-decision-and-stopping-rules"></a>
<a id="continuing-toward-the-general-answer"></a>
<a id="read-by-claim"></a>

## Working contract and reading map

Prefer analytic identities and a few exact fixtures. Preserve literal phases,
every logical column, coherent address and arbitrary dirty/reference inputs,
actual native adjoints, and all clean leakage and work return. No QRAM, free
loading or precision states, unpriced clean history, resets or postselection.
A suffix is clean only on its active sector; inactive action needs its own
proof. Separate general theorems, conditional interfaces, components and
representation limits. Internal audits and finite checks do not replace
analytic proof or external technical review.

| Information that changed | Update its authoritative home |
| --- | --- |
| Construction or proof | Existing proof home; [research catalogue](research/README.md) locates it. Change the index only for a new home or changed role. |
| Mechanism assessment or reopening condition | [Route decisions](research/ROUTE_DECISIONS.md), linked to the proof. |
| Bottleneck or next action | Replace this checkpoint; keep detailed outcomes in existing proof/route homes. |
| General theorem or global gap | [Resource frontier](docs/OPEN_PROBLEM.md). Selected A–D statements belong in [publication scope](manuscript/PUBLICATION_SCOPE.md); [proof map](docs/README.md) owns dependencies. |
| Evidence or attribution | Relevant tests and [verification catalogue](docs/reference/VERIFICATION_CATALOGUE.md); [source catalogue](docs/reference/SOURCE_CATALOGUE.md) owns attribution and stable IDs. |
| Validation procedure or upstream baseline | [Research validation](scripts/README.md#research-validation) or [SYNC](SYNC.md), respectively. |

<a id="refurnishing-policy"></a>

Update only affected homes. Keep summaries as links plus decision-relevant
scope. Preserve unique attempts, stable identifiers, inbound anchors, tests
and receipt hashes; change mechanisms when complete synthesis restores the old cost.

## Restore and verify

Check worktree status and remote, fetch main, and compare it to the scientific baseline above.
Preserve unrelated work. Follow the proof links for the chosen operation.

[Research validation](scripts/README.md#research-validation) owns focused
iteration checks and the complete ten-command final gate. Complete focused
[visual checks](assets/README.md#rendering-checks) before the memory-heavy
numerical suite; full presentation CI covers the corpus. For proof moves,
compare mathematics and audit old and new links. Record exact-tree receipts
in the review record and merge after required checks pass.
