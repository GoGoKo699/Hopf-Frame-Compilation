# Five-hour design: joint orientation and recovery

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) · [Current frontier](../../docs/OPEN_PROBLEM.md) · [Permutation checkpoint](PERMUTATION_BRIDGE_SESSION_RESULT.md)

**Status: executed.** The [session outcome](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
records R85's geometric-band construction, its reviews, and verification.
The design baseline was main commit `ed5597e`; execution began from
`d6259c2` on 15 October 2026 at 22:04:30 UTC, after the design was merged.
This document preserves the executed instructions. The
[older five-hour prompt](RESEARCH_SESSION_PROMPT.md) records the separate
R78 search.

## 1. Destination and the change of method

Work only in **GoGoKo699/Hopf-Frame-Compilation**, for at most five hours
of elapsed session time, including review, tool waits, and integration.
Use the existing authorization to modify and merge the repository.
This session pursues a complete construction for arbitrary nonlinear
branch orientations. Finding the mechanism is part of the work.

The previous conversion search asked for a cheap in-place leaf
permutation. This time the unit of design is the **complete oriented
frame, or a source followed by its recovery**. Intermediate work may
remain correlated across those operations. Only a proved identity can
justify reusing it or declaring it returned.

The question to keep visible is: **which repeated costs are required by
the whole operator, and which arise because we implement its factors
separately?** Confidence and the appeal of the problem motivate sustained
exploration; they do not determine the mathematical outcome. A pause
must preserve what was learned and a precise question from which to resume.

## 2. Freeze the target, allow the representation to change

Let $`N=2^n`$, fix $`\eta=1/64`$, and supply one external clean qubit
and $`b=\lceil C_b\sqrt N\rceil`$ arbitrary dirty qubits, for a fixed
sufficient constant $`C_b`$. Seek, on the same native Clifford+T circuit,

```math
T=O(\sqrt N),\qquad G=O(N),\qquad
D_T=o(\log^2(n+2)),
```

```math
\|VJ_1-J_1(W_r\otimes I_b)\|_{\rm op}\le\eta.
```

Here $`J_1`$ appends the sole clean qubit in $`|0\rangle`$, and
$`I_b`$ is the identity on all b arbitrary dirty qubits.

Logarithmic depth is the stretch goal. A complete
$`O(\log n\log\log(n+4))`$ construction would already cross the
current target. State a finite-dimension fallback and constants uniform
in n and the admitted orientation and magnitude tables.

At depth j, the unordered magnitude pair $`\{a_j,b_j\}`$ is shared
by all nodes, with $`a_j^2+b_j^2=1`$. Each node may choose its ordering
by an arbitrary Boolean function of its prefix. Begin the symbolic
analysis with positive pairs and an arbitrary orientation table; a
fixed interior bias is a diagnostic instance. A theorem for the full
promised family must also treat zeros, node signs, and supplied leaf
phases through explicitly priced constructions. A fixed-bias theorem
alone must be labelled with that narrower scope. Arbitrary node-dependent
magnitude pairs and optimal unrestricted depth are broader questions.

Keep every logical column, literal marker sign, actual inverse,
initialized-isometry error, and dirty/reference input. T-depth permits
arbitrary Clifford interlayers, whose elementary gate count belongs in G.
There is no QRAM, measurement, reset, uncharged coherent evaluator,
extra initialized output word, or supplied precision state. Classical
table preprocessing may terminate without an efficiency guarantee.

Protected logical zero suffixes remain available only in their certified
active sector. Save its predicate before occupying the suffix, account
for arbitrary occupied inactive inputs, and restore the predicate by
the actual inverse. Inherited exact-return primitives keep their exact
contracts. The final approximation criterion already accounts for any
allowed residual leakage; do not impose exact final cleanup as an
additional restriction on every proposed approximate construction.

## 3. Established starting points

Read these interfaces, then work on the missing identity rather than
rereading the whole archive:

| Interface | What is already available |
|---|---|
| [Whole-register reduction](PERMUTATION_BRIDGE_SESSION_RESULT.md#1-the-sufficient-target-and-its-scope) | A full in-place compiler is sufficient, but is not necessary for a direct frame word |
| [R78 nested flag](SUBTREE_PROJECTOR_DECODER.md#2-the-literal-nested-subtree-spaces), [R79 ledger](UNIFORM_SUBTREE_DEPTH.md) | Orthogonal depth bands, a complete tagged source, and a paid spectral decoder; logarithmically many serial source/reflection calls remain |
| [R80 actual common-flag word](HOMOGENEOUS_MAGNITUDE_DEPTH.md#5-one-protected-stage-and-exact-tag-return) | A cheap simultaneous flag conjugator erases the tag without adaptive spectral comparisons in the homogeneous case |
| [R81 terminal](RECURSIVE_TERMINAL_DEPTH.md#1-fixed-accuracy-terminal-theorem) | The deepest logarithmically many general layers have a paid log-times-log-log depth bound |
| [R83](AFFINE_ORIENTATION_DEPTH.md), [R84](DETERMINISTIC_FRAME_DEPTH.md) | Affine orientations and deterministic corners are solved within their stated scopes; they are diagnostic limits, not new targets |
| [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) | Exact forward/inverse XOR queries and Boolean phases; they do not supply in-place permutations or a clean output word |

If a new early stage leaves only $`O(\log n)`$ layers, attempt to
compose it with R81 before designing another terminal. Verify its actual
height, preserved-prefix, width, and error premises. An early-stage
improvement is not a full theorem until that composition is paid.

## 4. Two independent opening lanes

### A. Compile the relative frame as one operator

Reserve B for the positive balanced frame. The established identities are

```math
W_r=P_F W_0Q,\qquad B=P_F BQ,
```

```math
K_r:=W_rB^\dagger
 =P_F(W_0B^\dagger)P_F^{-1},\qquad W_r=K_rB.
```

**Research hypothesis:** the conjugated operator $`K_r`$ may admit a
native construction cheaper than separately implementing its two leaf
permutations. The equation itself supplies no circuit or depth saving.

The first symbolic task is to derive and verify a concrete recursive
representation. For a positive subtree v of height m, write

```math
W_v=(W_{v0}\oplus W_{v1})R_v(\theta_v),
```

where $`R_v`$ rotates the two child zero anchors. Chronologically the
root rotation acts first, followed by the child frames. With
$`K_v=W_vB_v^\dagger`$, derive

```math
K_v=(K_{v0}\oplus K_{v1})J_v(\delta_v),
\qquad \delta_v=\theta_v-\pi/4.
```

Here $`J_v(\delta)`$ is identity off the span of
$`|0\rangle|+\rangle^{\otimes(m-1)}`$ and
$`|1\rangle|+\rangle^{\otimes(m-1)}`$, and acts there by the literal
real rotation through delta. In the full register it is controlled by
the computational prefix v. Orientation changes
$`\delta_v`$ between $`\delta_{|v|}`$ and
$`-\delta_{|v|}`$. This is an ordered word with computational prefix
controls and uniform-state suffix projectors. R80 does not directly
compile that mixed-basis interface.

Verify this identity on the whole two-level subtree algebra, including
every marker column. Compute the adjacent parent/child commutator before
attempting parallelism. Then seek an exact fusion or echo that combines
those factors while retaining and eventually returning their shared
work. Compare balanced, interior, affine, and corner instances. The
pivotal question is whether this ordered word can avoid an independent
full-address query at each depth. Merely deriving its recursion does
not establish a native saving.

The first lane checkpoint must contain:

1. An exact block identity or a precisely stated conjectured native word.
2. Its action on every live register, including arbitrary dirty inputs.
3. The specific separate-factor cost it could remove, with a provisional
   count/depth/width recurrence and all unknown terms named.
4. One falsifiable lemma about a cancellation, selection, or transport
   operation that would decide whether the approach advances.

A useful candidate may retain address information or correlated dirty
work until the final correction. Do not require its invented intermediate
macros to be independent returned-work primitives. Instead track their
joint state and prove the final cancellation. Reinterpreting a lookup
table without physically restoring its bank is insufficient.

Reject circular implementations that call $`P_F`$, $`W_r`$, or an
equivalent selected router without a new native schedule. A decomposition
into n sequential levels has not improved the target merely because
each level is cheap. If a Hamiltonian or polynomial representation is
proposed, price its encoding, normalization, degree, and complete recovery
before using it as the new circuit model.

### B. Recover the band label using the freedom inside each band

For an h-bit local early frame, let $`A_d^r`$ be its nested scaling
projectors and let

```math
\Pi_0^r=A_0^r,\qquad
\Pi_t^r=A_t^r-A_{t-1}^r\quad(1\le t\le h).
```

Let $`t(j)`$ be the last-one position of marker input j, with
$`t(0)=0`$. The ideal tagged source is the isometry

```math
S_r|j\rangle=|w_j^r\rangle|t(j)\rangle_T,
\qquad |w_j^r\rangle=W_r|j\rangle.
```

An abstract complete band-label XOR is

```math
E_r=\sum_{t=0}^h\Pi_t^r\otimes X_T^{\mathrm{bin}(t)},
\qquad E_rS_r|j\rangle=|w_j^r\rangle|0\rangle_T.
```

The XOR definition acts on every tag value, including padded ones. It
does not assert a cheap native implementation. A joint source/recovery
word needs the displayed action on the reached tagged-source subspace;
it need not implement this particular extension everywhere else.
It still must satisfy the entire active/inactive and dirty-input contract.

For a concrete way to expose the available freedom, define

```math
C_d=I_{2^d}\otimes|0^{h-d}\rangle\langle0^{h-d}|,
\quad \Delta C_0=C_0,\quad \Delta C_t=C_t-C_{t-1}.
```

If U transports every $`C_d`$ to $`A_d^r`$, then

```math
E_r=(U\otimes I_T)
 \left(\sum_{t=0}^h\Delta C_t\otimes X_T^{\mathrm{bin}(t)}\right)
 (U^\dagger\otimes I_T).
```

Replacing U by $`UD`$, for any unitary D acting separately within the
$`\Delta C_t`$ bands, leaves this identity unchanged. Thus the
intermediate transporter may mix vectors inside a band. The complete
source/recovery composition must still produce each prescribed column
with its literal phase. Choosing $`U=W_r`$ only restates the problem.
The within-band D is freedom to design the final native word, not an
additional free unitary that may be inserted into the circuit.

**Research hypothesis:** this freedom permits a joint band-label word
with fewer serial large queries than R79, even though the individual
marker permutation is expensive in the present decomposition.

Begin by deriving a joint encoding of two or more band-label bits, or a
common flag transport, together with its actual inverse. The exact bit
phases

```math
Z_{r,k}=\sum_{t=0}^h(-1)^{\mathrm{bit}_k(t)}\Pi_t^r
```

commute algebraically. Determine whether a shared native construction
can use that fact without a separate source bank for each bit or a new
sequence of adaptive full-address queries. A direct product of threshold
reflections has the known linear-boundary cost; changing their notation
does not evade it. A classical XOR computed only from the output leaf
cannot identify overlapping positive marker bands.
The first nontrivial character to examine is depth parity, k zero;
another single threshold reflection would repeat an available primitive.

An equivalent concrete target for the shared word is the joint character

```math
\mathcal Z_r=\sum_{z\in\{0,1\}^s}
 \left(\sum_{t=0}^h(-1)^{z\cdot\mathrm{bin}(t)}\Pi_t^r\right)
 \otimes|z\rangle\langle z|_T,
\qquad s=\lceil\log_2(h+1)\rceil.
```

Hadamards on the existing tag before and after this word give $`E_r`$;
no second clean label is required. Try to construct the character from
the tree data and a shared recovery word. Writing it as the unknown
frame conjugating a computational phase is circular. Its inactive
implementation must make the complete Hadamard/character/Hadamard word
identity on occupied tags; all character controls and work are charged.

The decisive candidate lemma must specify the source-correlated input
sector, the joint eraser or transporter, its live work, and the number
and native depth of its large table accesses. It may leave one named
construction step open at the discovery gate. A single descendant
loader's unspecified completion is not a transporter of the entire flag.

For error composition, if the actual source has isometry error delta and
a unitary recovery has error epsilon on the ideal source image, their
composition has error at most $`\delta+\epsilon`$ by the triangle
inequality. Establish those bounds uniformly on the initialized and
dirty/reference sectors; individual-column estimates are insufficient.
R79's transported all-zero-source-work projector must not be replaced
by a logical projector tensored with identity on arbitrary occupied work.

## 5. Discovery gates and adversarial review

Use a coordinator, two constructors, a native-resource auditor, and a
mathematical skeptic. Keep capacity for a narrow source check or a fresh
reviewer. The constructors begin independently, then exchange their
strongest objections and any compatible identities.

At the first gate, require a precise hypothesis, a plausible reason that
it changes the serial cost, and a named missing lemma. A finished general
compiler is not an admission requirement. If neither opening lane has a
promising representation, use the reserved reformulation block to change
the encoding, the permitted within-band action, the location of retained
history, or the separation between source and recovery. Do not simply
enumerate more versions of the same failed word.

Keep one primary candidate and at most one distinct reserve after that
block. Judge progress by what uncertainty has been resolved. A candidate
with a narrowing, explicit gap may continue without a finished theorem.
A proved failure ends that candidate under its assumptions; it does not
end the project. Do not substitute another solved subclass or an
unrestricted-impossibility claim for the constructive destination.

The auditor and skeptic must challenge:

- Hidden clean output, copied controls, selected wire permutations, or
  encoders whose implementation is the original problem.
- Unchanged leading depth after all calls, actual inverses, arithmetic,
  precision, recovery multiplicity, and the terminal are counted.
- Claimed cancellations between ideal factors that fail for the actual
  chosen native words, or for literal phases and occupied work.
- Global prefixes omitted from addresses, simultaneous helper pools
  counted as reusable, or shared non-Clifford supports called parallel.
- A state-column argument substituted for the complete frame, including
  normalized marker columns below very small incoming amplitudes.
- Algebraic commutativity substituted for a physical schedule, or exact
  encoded-order restrictions extended to approximate encodings.

The [previous session](PERMUTATION_BRIDGE_SESSION_RESULT.md) and
[joint-depth record](JOINT_DEPTH_SESSION_RESULT.md) delimit four-shear
conversion, masked-bank hoisting, leaf-only erasure, threshold products,
and restricted word/encoding arguments. They guide which assumption to
change; they do not bound unrestricted native circuits.

## 6. Five-hour allocation

| Elapsed minutes | Work and checkpoint |
|---|---|
| 0–25 | Verify the current baseline; freeze the full target and identify the repeated cost. Write one question per lane. |
| 25–95 | Derive the two representations independently. Produce exact algebra where available, provisional recurrences, and the pivotal unknown. |
| 95–135 | Exchange objections and deliberately reformulate failed premises. Select a primary candidate and, if useful, one genuinely different reserve. |
| 135–220 | Concentrate on the pivotal identity and its native schedule. Keep the resource and return audit alongside the derivation; switch to the reserve only for a stated mathematical reason. |
| 220–265 | Reconstruct the whole argument or delimit its remaining gap. Start independent proof reconstruction and minimal relevant validation as soon as a stable statement exists. |
| 265–300 | Finish review, validate the reviewed tree, integrate permitted changes, and preserve a precise resumption question. |

The clock is a ceiling, not a promised theorem or an obligation to spend
time on an unchanged dead end. Do not end discovery merely because the
first candidate list failed. Complete the reformulation pass and evaluate
its remaining question. If a candidate still has useful unresolved work,
use the available research block. Report actual elapsed time accurately.
Begin long validation early enough to fit; an unfinished required check
must remain explicitly pending rather than being bypassed at the deadline.

## 7. Evidence, outputs, and integration

Keep one concise working checkpoint with the target, current identity or
conjecture, changed assumption, live-register contract, resource/error
recurrences, strongest objection, and the exact next derivation. Save
expensive reasoning on the research branch so the work is resumable.

Follow the theorem-first workflow: write the symbolic identity or
encoding and resource/error/cleanup recurrences before adding fixtures.
Use a small exact check only for a specific surviving premise or native
word. No broad numerical search or asymptotic claims from finite fits.
A focused literature check is useful only for a named missing primitive;
read primary sources and translate their workspace, precision, access,
and return assumptions before importing a claim.

| Outcome | Honest retained result |
|---|---|
| Complete target | A full promised-family circuit below square-logarithmic depth with all simultaneous bounds and error/return proof |
| Proved mechanism | A native improvement to a named bottleneck, with its scope and remaining full-circuit obligation explicit |
| Unresolved mechanism | New algebra or a materially sharper conjectured interface, one precise missing lemma, and a concrete resumption question |
| Candidate failure | A scoped counterexample or obstruction and the assumption it changes for the next attempt; no general impossibility claim |

A fresh reviewer must reconstruct the pivotal lemma and final
composition before any theorem is promoted. A discovered inherited
defect takes priority over using that result. Distinguish proved facts,
imported premises, hypotheses, and finite evidence throughout.

Retain one proof home for new mathematics. Record the execution outcome
in a compact research note, update this prompt's status and its index
link, and update the frontier only if a proved bound changes. Existing
R79–R84 statements and immutable scientific receipts remain authoritative.
Do not create a new result label for a plan or an unresolved conjecture.
Final manuscript drafting and external reviewer outreach are separate.

Use the [workspace checks](../../WORKSPACE.md#restore-and-verify) for
substantive mathematical or circuit changes, and the relevant document
and presentation checks for a plan or unresolved note. Merge only after
required checks pass on the reviewed tree. If the five-hour ceiling
arrives first, preserve the exact branch, proof status, and pending
checks without claiming completion. The final report states the strongest
proved conclusion, its limits, what was learned, the next unresolved
question, verification, elapsed time, and repository state.
