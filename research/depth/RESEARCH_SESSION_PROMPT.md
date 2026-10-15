# Runnable prompt for a five-hour depth research session

[Preparation map](DEPTH_OBSTACLE_MAP.md) · [Current frontier](../../docs/OPEN_PROBLEM.md) · [Workspace mandate](../../WORKSPACE.md)

**Status: designed, not yet executed.** This prompt starts a new research
session when the user asks to run it. Its baseline is main commit
`321101d`, containing the independently reviewed preparation map. Check
for subsequent changes before execution; a later result may supersede a
route or its starting bound. This document proposes attacks and decision
rules. It establishes no new mathematical claim.

## 1. Mission and success criteria

Work in **GoGoKo699/Hopf-Frame-Compilation** for at most five hours of
elapsed session time. Pursue a genuine change of mechanism in the remaining
large-width T-depth problem. Use the available time autonomously within
the existing repository authorization. Progress updates do not require
stopping for user input. Five hours is a ceiling, not a requirement to
continue an exhausted route or to report a positive theorem.

The primary target is a complete-family construction with

```math
\eta_0=1/64,\quad N=2^n,\quad a=1,\quad
b=\lceil C_b\sqrt N\rceil,\quad
T=O(\sqrt N),\quad G=O(N),\quad
D_T=o(\sqrt n\log(n+2)).
```

Here $`C_b`$ is a fixed sufficient constant, not a function growing with n.
All constants in the gate bounds are independent of n and the angle tuple.
A polylogarithmic depth bound is the stretch target. Do not assume that
it is achievable, or discard a rigorous smaller asymptotic improvement.
Cover sufficiently large n uniformly and state the finite-dimension
fallback. The initial family is the complete real Hopf frame; extending
to shrinking error and complex phase dressing comes after the main proof.

A dimension-growing lower bound is an alternative objective, eligible
only after an explicit internal selection decision. State its quantifiers
and circuit class. A bound using the simultaneous T/G budgets must be
labelled as such; an unrestricted lower bound may not import those
restrictions. An additional constant-depth exclusion for the existing
spectator-padded witness is outside this session's main success criterion.

Rank the actual outcome honestly:

| Outcome | Required content |
|---|---|
| Complete theorem | A new asymptotic upper bound under the benchmark, or a selected dimension-growing lower bound with its precise circuit class |
| Complete native mechanism | A proved improvement to a dominant component, its full ledger, and the exact remaining composition obstacle; not yet a compiler theorem |
| Decisive route exclusion | A proved scoped obstruction or exact counterexample to a precise promising premise, with surviving alternatives stated |
| Sharper unresolved lemma | A materially narrower obligation and useful failure evidence; explicitly not a proved result |
| No qualifying advance | A concise record of the tested premises and why none survived; no expansion of the scientific claims |

A verified defect in an inherited theorem takes priority: isolate it,
obtain independent review, and repair or delimit the affected claim before
using it further.

## 2. Freeze the task, challenge the method

Retain the complete initialized-isometry guarantee

```math
\|VJ_1-J_1(W(\boldsymbol\theta)\otimes I_b)\|_{\mathrm{op}}
\le\eta_0
```

for every admitted angle tuple, including signed, unequal, and singular
tuples, every logical column, and arbitrary dirty/reference inputs.
Preserve literal phases, clean leakage accounting, actual inverses, and
work return. Charge every precision source, query, reflection, guard,
basis transformation, and classical-data interface. T-depth allows
arbitrary Clifford interlayers, with elementary Clifford count included
in G. No QRAM, reset, measurement, free evaluator, supplied precision
state, extra clean register, or state-column substitution is allowed.
Conditional zero suffixes require a saved original predicate and a proof
of identity on arbitrary occupied inactive work.

The existing tag basis, row/column support, amplification word, protected
split, query circuit, and lower-bound witness may all change. Ask first:
**which cost is forced by a theorem, and which belongs to our chosen
implementation?** The preparation map answers this for the current route;
it does not certify that any proposed replacement is cheap.

Read the [preparation map](DEPTH_OBSTACLE_MAP.md) and the relevant proof
homes, not the whole archive again. At fixed accuracy the starting ledger is

```math
D_T=O\!\left(\sqrt n\log n+\sqrt n\log\log n+
\log^2 n+N/b^2\right).
```

Its gradient is already prepared once and returned once. R67 bounds
entry-indexed normalization, not arbitrary native depth or all recovery
query complexity. R75–R77's current witness has fixed-accuracy constant
depth by the existing one-qubit compiler. These are starting facts, not
research tasks to rediscover.

## 3. Independent opening attacks

Start two constructive lanes independently. Each must propose a precise
new lemma, identify the old assumption it changes, and compute the whole
conditional ledger before requesting sustained attention. Keep only one
main candidate per lane. A short third lane challenges the upper-bound
framing by looking for dimension-sensitive hardness.

### A. Erase the depth label using nested subtree projectors

This is the preferred initial hypothesis, **not an established circuit**.
For the h-bit early frame, define normalized subtree scaling states
$`|\psi_{d,p}\rangle`$ directly from the descendant angles, without
dividing by an incoming amplitude that can vanish, and set

```math
A_d=\sum_{p\in\{0,1\}^d}
 |p\,\psi_{d,p}\rangle\langle p\,\psi_{d,p}|,
\qquad 0\le d\le h.
```

First derive, with the literal marker convention and every admitted tuple,

```math
A_d\le A_{d+1},\qquad A_h=I,\qquad
\Pi_0=A_0,\qquad \Pi_{d+1}=A_{d+1}-A_d.
```

The intended $`\Pi_0`$ is the state-column projector; $`\Pi_{d+1}`$
should be the span of all depth-d marker columns. The ideal tagged source in
[R66 §3](SUBLINEAR_T_DEPTH.md#3-compress-the-complete-columns)
has action $`W_h|j\rangle|t(j)\rangle`$; its native word approximates
that isometry. Try to implement the full spectral-label XOR

```math
E=\sum_{t=0}^{h}\Pi_t\otimes X_t,
\qquad X_t|z\rangle=|z\mathbin\oplus t\rangle.
```

If realized, E erases that tag and retains every logical column. The
formula only specifies the desired unitary; it does not price it.

**The first substantive gate is a native reflection lemma.** Use the
[descendant-state loader](COHERENT_DESCENDANT_LOADER.md) on computational
anchors $`|p,0^{h-d}\rangle`$. Let $`P_{\mathrm{anchor},d}`$ also
test the fresh initialized loader-work bits. Exclude the occupied old tag,
the variable-d/search record, and the saved activity flag from this zero
test. Establish that

```math
\widehat{\mathcal R}_d
=\widehat U_d(I-2P_{\mathrm{anchor},d})\widehat U_d^\dagger
```

obeys the initialized-isometry bound

```math
\|\widehat{\mathcal R}_dJ_{\mathrm{work}}
 -J_{\mathrm{work}}(I-2A_d)\|_{\mathrm{op}}\le\epsilon_d.
```

It must hold on arbitrary logical inputs and include leakage and
actual-inverse error; no full-space identity on occupied work is assumed.
The shared live gradient, if used, is treated in its invariant sector;
the reflection must not incorrectly test it for zero. Prove coherent
variable-d selection, leaving the old tag untouched and reading search
controls without overwriting them. Put the control and activity guard on
the central anchor phase, with the loader and its actual inverse outside
it. Prove exact control-off and activity-off identity on arbitrary occupied
protected work; do not discard a conditional scalar minus. Do not use W
or its inverse as a subroutine.

Next test whether coherent binary search among the nested projectors can
compute the sector label with $`O(\log h)`$ controlled reflections.
Give the reversible midpoint arithmetic, clean control bits, complete
history lifetime, and uncomputation. Compute the label, XOR it into the
old tag, and uncompute the new label without assuming access to the erased
tag. A classical binary-search description is not a quantum circuit.

Derive a new source-plus-erasure hybrid: the initial tagged-source error,
every controlled-reflection error, search compute/uncompute, and final
tag/work leakage must all be covered. In particular,
$`D_{\mathrm{early}}=D_{\mathrm{source}}+D_E`$ and
$`\epsilon_{\mathrm{early}}\le\delta_{\mathrm{source}}+\delta_E`$;
R66's old amplification error allocation does not automatically transfer.

The conditional target is $`O(\log h)`$ native reflections at
$`O(\log n)`$ depth each, hence $`O(\log^2 n)`$ early depth and the
same order after the retained tail. If there are q approximate calls,
derive an explicit allocation $`\delta=\Theta(\eta_0/q)`$ and its
precision cost. First try to refit the table and width ledgers with
$`M=\lceil A\log(n+2)\rceil`$ for a fixed sufficiently large A.
Polynomial overhead relative to the smaller early tables may fit the
exponential budgets; standalone arithmetic and live clean work require
separate estimates. Do not assume this fit before deriving the schedule.
A larger paid polylogarithmic suffix is eligible if needed for reversible
search, but then derive its terminal depth anew and state the resulting
bound; the old $`O(\log^2 n)`$ tail cannot simply be retained.

Immediate rejection tests are incorrect literal marker projectors,
an unpaid long clean index/history, uncontrolled reflection leakage,
failure of coherent variable-depth selection, or inability to return the
search record. Test balanced, unequal multilevel, and signed singular
families, then prove the all-tuple identity. Check
[joint history compression](../endpoint/JOINT_HISTORY_COMPRESSION.md)
and [common-source reuse](COMMON_SOURCE_REUSE.md): a named decoder or
transported reflection is not already compiled, and a static success
monitor may be invalid. Explain precisely why the realized decoder leaves
[R67's entry-indexed overlap architecture](ENTRY_INDEXED_ENCODING_LIMIT.md).

### B. Change the native cost of repeated large queries

Keep the current recovery initially and attack the exact large-address
operations that produce $`O(\log n)`$ per call. Begin with their actual
Boolean-phase, marker-extraction, calibration, and source-reflection
interfaces. Determine which inputs remain unchanged between uses and
which are altered by recovery. Prove any cache lifetime from that word.

The generic query contract is

```math
|x,y,w\rangle\longmapsto|x,y\mathbin\oplus F(x),w\rangle
```

for arbitrary y, w, and references. A phase or marker specialization is
eligible if it replaces every leading-cost occurrence needed by the
complete compiler. A faster promised table family is insufficient unless
all generated tables obey that promise by proof.

Explore whether the smaller early table volume permits a different
space/depth schedule, or whether a fused sequence avoids repeatedly
constructing the same nonlinear information. Recompute the protected
suffix size and the *simultaneous* T/G/width fit if the tradeoff changes.
Clifford address fanout is permitted; it does not supply initialized
copies, erase unknown masks, or prove that non-Clifford supports are
disjoint. Give the physical parallel schedule.

Use [shallow dirty queries](SHALLOW_DIRTY_QUERIES.md),
[nonuniform indicators](NONUNIFORM_DIRTY_INDICATOR.md),
[shared-prefix limits](SHARED_PREFIX_QUERY_AUDIT.md), and
[common-source reuse](COMMON_SOURCE_REUSE.md) as collision checks.
Reject stale monitors, unchanged-address assumptions contradicted by the
recovery word, and table differencing that reconstructs a generic lookup.

As a conditional decision test, reducing every leading query to
$`O(\log\log n)`$ while retaining the other allowances would give
$`O(\sqrt n\log\log n+\log^2 n)`$. This is a target calculation,
not a claimed primitive. Arithmetic, reflections, recovery multiplicity,
and the terminal tail must still be counted.

### C. A short challenge from dimension-sensitive hardness

For at most forty minutes, seek a witness family whose active structure
grows with n and a robust invariant that detects it under arbitrary
Clifford interlayers and one-clean dirty return. Two concrete candidates
are the balanced Haar frame and the deepest-layer Boolean multiplexor

```math
U_f=\sum_{x\in\{0,1\}^{n-1}}|x\rangle\langle x|
 \otimes R_y(\theta_{f(x)}),
```

obtained by setting all upper angles to zero and choosing two fixed
final-layer rotations separated by more than $`2\eta_0`$ in the relevant
operator metric. Compare generic truth tables with a
point-indicator table. Neither family is presumed hard; a product of
legal frames need not itself belong to the family.

Before developing a lower-bound proof, try to compile the proposed witness
cheaply. Test spectator reduction, Clifford conjugation, disjoint or
parallel decomposition, and use of the existing count/depth compiler.
Many parameters or many active qubits alone do not establish depth
hardness. Ordinary light cones, support size, and entanglement growth need
separate justification because Clifford interlayers can be arbitrarily deep.
In particular, $`H^{\otimes n}`$ has T-depth zero and entrywise absolute
norm $`\sqrt N`$; a phase depending on a parity can be Clifford-conjugate
to one local phase gate. Reject an invariant that incorrectly charges
either example. The existing word-count exponent $`O(q^2+qT)`$, with
physical width q and $`q,T=\Theta(\sqrt N)`$, is already $`O(N)`$;
reusing that count against a fixed-error frame packing gives no growing
depth bound without a new justified restriction.

For the multiplexor, test a structural lemma about the *joint* occupied
conjugation of the preserved address algebra and the target Pauli. Write
the claimed per-layer restriction before using Boolean degree or table
complexity. Approximate return cannot be promoted to exact return at
general depth using the repository's at-most-four-layer result.

Write the desired invariant bound for *every* eligible native circuit,
its quantitative stability under fixed operator error, and its dependence
on width. A claimed growing lower bound must give the quantifiers

```math
\begin{gathered}
\forall n\ge n_0\ \exists W_n\in\mathcal F_n\ \forall V\in\mathcal C_n:\\
D_T(V)\lt f(n)\ \Longrightarrow\
\|VJ_1-J_1(W_n\otimes I_b)\|_{\mathrm{op}}\gt\eta_0,\\
f(n)\longrightarrow\infty.
\end{gathered}
```

Specify $`\mathcal C_n`$, any fixed gate-budget constants, and whether
the witness is explicit.
If several observable errors are combined, budget their aggregation
uniformly in n; an n-term error sum is not $`O(\eta_0)`$ without a proof.
If clean compression, crossing rank, or inner spread enters,
show why the parameter is controlled; changing its name is not progress.
By the first selection gate, either present a credible dimension-growing
mechanism or return this worker to adversarial review of the upper route.
Do not fall back to the old single-qubit witness to obtain a consolation
constant-depth result.

## 4. People, clock, and selection gates

Use a coordinator and up to five workers: two constructors, the short
hardness challenger, an adversarial mathematician, and a native/resource
auditor. Keep one concurrency slot available for a narrow literature check
or fresh final reviewer. Do not allocate workers merely to reread all files.
After selection, the losing constructor becomes a checker or works on the
reserve candidate's single weakest lemma.

| Elapsed minutes | Work and required output |
|---|---|
| 0–20 | Verify the baseline and freeze the contract. Each constructor writes its changed assumption and proposed lemma independently. |
| 20–55 | Derive the two candidate mechanisms and conditional whole-circuit ledgers. Run the short hardness challenge in parallel. |
| 55–75 | Exchange the strongest objections. Select one primary mechanism and at most one qualified reserve. Explicitly record any switch to a lower-bound objective. |
| 75–155 | Concentrate on the pivotal lemma. Resource and counterexample audits run alongside symbolic work. |
| 155–175 | Recompute the complete bound. Continue only if the uncertainty is narrowing; otherwise switch once to the already-qualified reserve. |
| 175–235 | Finish or precisely delimit the selected argument. Begin fresh proof reconstruction as soon as a stable lemma exists. No new broad lane. |
| 235–275 | Complete independent review and relevant validation; integrate verified material. Allow CI to run while the final composition is checked. |
| 275–300 | Resolve review findings, complete permitted integration, and report the actual outcome. Save a resumable checkpoint if any gate remains unfinished. |

There are at most two mechanism-selection pivots: the first gate and the
reserve switch. Earlier decisive falsification permits an earlier gate,
not unlimited variants. If no candidate qualifies at the first gate,
allow one twenty-minute reframing pass using the *specific* failed premises;
count it as the reserve switch and shorten the proof phase accordingly.
If it also fails, end with an honest negative checkpoint. Do not consume
the remaining hours manufacturing conditional reductions.

A targeted literature check is justified only by a precise missing lemma
or a proposed transfer from another model. Use primary papers or official
sources, read their assumptions, and record the model translation: clean
ancillas, measurements, gates, precision, access model, and returned work.
An oracle-query theorem does not automatically imply native T-depth.
Stop the search once it settles that premise; avoid a general literature
survey and do not make unverified novelty claims.

## 5. The candidate must survive its strongest objection

At every gate, reverse roles. Ask:

1. Which line actually leaves the old architecture, and what replaces it?
2. If the proposed lemma were true, what is the largest term in the complete
   circuit after all multiplicities and precisions are substituted?
3. Does a hidden subroutine implement the original problem under another
   name? Draw an acyclic dependency graph separating inherited and new lemmas.
4. Which legal input is most likely to defeat the claim? Include coherent
   superpositions across active/inactive sectors and arbitrary dirty inputs.
5. For a lower bound, which legal native circuit lies outside the proposed
   argument? Is the claim about all representations or only a convenient one?

Do not score candidates by elegance or confidence. Prefer a larger provable
asymptotic gain with fewer unresolved interfaces; a small complete gain
outranks a dramatic claim whose decoder is unpaid. Stop a construction if
its full ledger retains the old leading order, unless it proves a distinct
useful mechanism whose remaining obstacle is explicitly identified.

For each surviving construction maintain

```math
D_{\mathrm{total}}=D_{\mathrm{boundary}}
 +\sum_{j=1}^{q}D_{\mathrm{call},j}+D_{\mathrm{tail}},
\qquad
\epsilon_{\mathrm{total}}\le
\epsilon_{\mathrm{boundary}}+\sum_j\epsilon_j+\epsilon_{\mathrm{tail}}.
```

Expand each call into queries, arithmetic, reflections, guards, and actual
inverses. Give corresponding T and G sums. Width is the peak live
reservation of the physical schedule: sequential pools may be reused,
concurrent pools may not. Distinguish exact returned primitives from
approximate initialized-isometry composition. Use the shared-gradient
comparison only within its invariant sector and account for exterior
preparation error separately.

Do not transfer ideal angular stability to physical rejection leakage,
replace actual inverses with ideal inverses, or assume intermediate clean
work from terminal return. A small normalized trace is not a small operator
norm. A hypothetical algebraic counterexample must satisfy native
realizability before it can refute a native-circuit claim.

## 6. Checkpoints, proof review, and repository integration

Keep one short current checkpoint per candidate, updated at the gates:
target/status; changed assumption; exact lemma; dependency graph;
native operation order and input sectors; T/G/depth/width/error ledger;
strongest stress case; and continue/pivot/proved/disproved/unresolved
decision with one precise reason. Preserve these during a long session so
another coordinator can resume without reconstructing the conversation.

Prove identities and recurrences before building fixtures. Use small exact
checks only to discriminate a stated premise or validate a new native word.
Do not run a broad numerical search or treat a finite fit as an asymptotic
proof. A fresh reviewer must reconstruct the pivotal lemma and the final
composition independently. Separate proved facts, imported premises,
unproved lemmas, and finite evidence in the resulting note.

Retain one proof home per result. If a genuine result survives, write its
scoped statement and proof, add only relevant evidence, and update the
frontier and index without duplicating proofs. If the outcome is negative
or unresolved, preserve only the unique mechanism, precise failure, and
reopening condition in a compact research note; leave selected claims
unchanged. Do not draft the manuscript or redesign the repository.

Use the [workspace validation instructions](../../WORKSPACE.md#restore-and-verify)
for substantive mathematical/circuit changes. Run checks on the reviewed
tree and merge only after the required gates pass. Existing authorization
allows routine repository integration; do not request it again. Never
weaken tests, regenerate exact receipts merely to pass, or mark pending
checks successful. The five-hour cap includes tool waits. If validation or
review cannot finish before it, leave the exact tree and pending work
identified in the checkpoint rather than bypassing the gate.

The final report must state the outcome tier, the strongest proved
statement and its scope, what changed in the obstacle map, verification
performed, the remaining precise obstruction, and the repository state.
Report elapsed time accurately. If no stronger bound was proved, say so
directly. Make the next decision from evidence, not from the wish to keep
the session moving.
