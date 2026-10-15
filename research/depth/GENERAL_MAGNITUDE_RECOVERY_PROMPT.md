# Five-hour design: general magnitudes from exposed clean space

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) · [R85 proof](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) · [Previous session](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md) · [General frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth)

**Status: designed, not executed.** The design baseline is main commit
`d52efbc`, containing the reviewed R85 result. Record the actual execution
baseline and start time when this prompt is run. This document is a
research plan; none of its proposed primitives is an established theorem.

## 1. Priority and destination

Work only in **GoGoKo699/Hopf-Frame-Compilation**, for at most five elapsed
hours including discovery, review, validation, tool waits, and integration.
Repository modification and merge are already authorized. Preserve the
separation between scientific work and final manuscript drafting.

The primary target is to remove R85's shared-unordered-magnitude promise.
Its geometric unloading identity already holds for general frames. The
promise enters the local flag transporter: a triangular permutation of a
tensor-product reference works when the magnitude pair depends only on
depth. With node-dependent pairs that particular reference disappears.

The question for this session is: **can the clean workspace exposed on a
reached band support a cheap general flag-label eraser?** The changed
resource is a linear number of conditionally initialized local wires,
available together with the tagged source image. It is not an extra
external clean register supplied on every logical input.

Improving R85 further to logarithmic depth remains a secondary direction.
That requires a compatible improvement of coarse recovery, terminal
processing, and selected routing. Removing only its coarse search does
not suffice. Give priority to removing the structural promise rather
than beginning three separate optimizations of the solved family.

Confidence in the destination supports exploration; it does not decide
whether a proposed identity is true. Complete a deliberate reformulation
pass before ending discovery. A failed candidate leaves a scoped lesson
and a precise next question, not a conclusion that the general target
cannot be reached.

## 2. Freeze the full-circuit benchmark

Let $`N=2^n`$, fix $`\eta=1/64`$, and supply one external initialized
qubit and $`b=\lceil C_b\sqrt N\rceil`$ arbitrary dirty qubits, for
a sufficiently large constant $`C_b`$ independent of n and the frame.
Seek one native Clifford+T circuit V with

```math
T=O(\sqrt N),\qquad G=O(N),\qquad
D_T=o(\log^2(n+2)),
```

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta.
```

Here $`J_1`$ appends the sole external qubit in $`|0\rangle`$, and
$`I_b`$ is identity on all b supplied arbitrary dirty wires.
W is a prescribed complete real frame with arbitrary node-dependent
normalized pairs, including signs and zeros, or its prescribed complex
magnitude completion with leaf phases. Begin discovery with finite
positive pairs; a completed general claim must restore the full scope
through priced certified rounding, sign gauges, and phase composition.
An interior-only theorem is a restricted result, not the primary target.

Matching R85's log-times-log-log depth for general magnitudes is the
preferred target; any strict improvement below square-logarithmic depth
with these simultaneous budgets counts. First prove the large-width,
fixed-accuracy statement. Extension to all eligible widths or absolute
all-precision constants is a separate obligation, not an opening demand.
Give a finite-dimension fallback before promoting a theorem.

Keep every logical column, literal phase, actual inverse, work leakage
and return, arbitrary dirty input, and external reference. Charge elementary
Clifford gates even though arbitrary Clifford interlayers are allowed by
the T-depth model. There is no QRAM, measurement, reset, free coherent
evaluator, supplied precision state, or extra external initialized word.
Terminating classical preprocessing need not be efficient.

## 3. The exact missing local interface

Use the notation of R85: $`\ell=\log_2(n+2)`$,
$`M=\lceil A\ell\rceil`$, $`h=n-M`$, and $`Q=2^h`$ initially.
The old fine tag T selects a geometric bin with endpoints a and d and
width $`q=d-a`$. The retained bin register R identifies that bin.
After the selected descendant unloader, the ideal reached subspace has
a preserved computational prefix u of length a, a q-bit local frame
band, and a zero suffix. This statement holds on coherent sums of
columns and prefixes, including arbitrary magnitudes.

Let $`A_j^{u,q}`$ be the local scaling projectors and put

```math
\Pi_0^{u,q}=A_0^{u,q},\qquad
\Pi_j^{u,q}=A_j^{u,q}-A_{j-1}^{u,q}\quad(1\le j\le q).
```

The minimally sufficient new word E has, on the reached input,

```math
E\bigl(|\psi_j\rangle_X|a+j\rangle_T
       |0^r\rangle_C|\xi\rangle_D\bigr)
\approx
|\psi_j\rangle_X|0\rangle_T
       |0^r\rangle_C|\xi\rangle_D,
\qquad |\psi_j\rangle\in\mathrm{Ran}\,\Pi_j^{u,q}.
```

The first bin includes j zero; later bins have $`j\ge1`$. C denotes
only the exposed initialized wires reserved for this local construction;
D denotes arbitrary dirty work, possibly entangled with a reference.
The error must be a single operator-norm bound on the whole tagged
subspace, uniform in u and the selected bin. Separate column estimates
do not establish it. The complete word restores the original vector
and its phase; intermediate transformations may mix within a band.
The default composition also preserves u, the coarse label R, and the
saved activity bit H on its ideal output. A joint construction may change
or erase these records too, but then must supply a new composition and
return argument rather than inherit the displayed R85 ledger.

The exposed capacity before new allocations is
$`r_{\rm available}=h-d+m_0`$, with
$`m_0=\Theta(\ell)`$ and $`r_{\rm available}\ge2q`$ in the
original geometric schedule. Tags, source banks, thresholds, and every
live history record are disjoint from this reserve. A construction must
allocate its actual C, routing work, and retained history within it.

An approximate transporter of the full local flag is sufficient, but
not necessary. A convenient stronger goal is a band-label XOR on every
input tag value. Do not require either stronger extension if a direct
word on the reached tagged source image closes the full circuit.
Likewise, E may use the old tag before erasing it; any subsequent inverse
that depends on that tag needs an explicit coherent return identity.

Require a physical unitary extension on occupied work. On the inactive
saved-predicate sector, the completed early word must have literal
identity action on arbitrary logical and work inputs. Prove error by
comparing with the ideal reached source image and propagating leakage
through actual unitaries. Never apply an initialized-tail lemma to
actual leaked states without that comparison.

## 4. Two independent construction lanes

### A. Use the clean bank to contract the unresolved flag

Derive an exact split of a general local flag across a cut, tracking the
data, cut label, descendant source, and every retained record. Seek a
native step that reduces the unresolved flag height substantially while
preserving the remaining quantum information in explicitly tracked
registers. The data need not remain in its original basis between steps.

The first falsifiable question is whether a step using $`O(q)`$ clean
local wires can reduce the unresolved height from q to at most
$`q^\alpha`$, for a fixed $`0\lt\alpha\lt1`$, with only a
constant number of charged large-source layers and only one serial
recursive child. The intended depth recurrence is

```math
F(q,\epsilon)\le F(q^\alpha,\epsilon')
                  +O(\ell+\log(1/\epsilon')).
```

It would need $`O(\log\log q)`$ stages. This is a design hypothesis,
not an available primitive; choose the recursive errors so their full
sum fits epsilon and price each resulting precision. Prove the initialized-work and return bounds
for the entire retained history; along that proposed contraction,
$`\sum_r q_r=O(q)`$ and $`\sum_r\log(q_r+2)=O(\log q)`$
are useful capacity estimates, not supplied ancillas.
A different recurrence with sub-square-logarithmic total depth is
equally admissible.

Write the actual recurrence before claiming progress. In particular,
a halving recursion with $`O(\ell)`$ depth per stage still has
$`O(\ell\log q)`$ depth. For the largest geometric bins q is of
order n, so that reproduces the old square-logarithmic cost. A native
step must change the contraction, the per-stage cost, or the physical
parallel schedule. Replacing the local transporter by the unknown local
frame, or invoking R79/R81 on that largest local frame, does not do so.
Likewise, two serial children of height square root q give the recurrence
$`F(q)\le2F(\sqrt q)+O(\ell)`$, whose usual bound is still
$`O(\ell\log q)`$. A smaller child height alone proves no saving.
All recursive multiplicities, forward calls, and actual inverses belong
in the recurrence. A maximum in place of their sum requires a physical
parallel schedule with all simultaneous workspace charged.

The first lane deliverable is an exact cut identity, a live-register
table, and one explicit contraction or shared-work lemma whose proof
would change this recurrence. Include the eventual return of any
retained cut history. Independent subspaces on one quantum data register
are not independent physical copies.

### B. Use the exposed space for a richer joint label encoding

Do not restrict the new encoding to a single leaf address. From the
original geometric capacity and $`d=a+q`$,

```math
a+3q\le h+m_0.
```

Thus a selected classical table with a prefix and as many as three local
q-bit indices can have at most $`2^{h+m_0}`$ rows. With the default
logarithmic reserve, that is $`Q\,\mathrm{poly}(n)`$. This is table
capacity, not a circuit theorem: output length, private dirty pools,
normalization, invocation count, and coherent return remain to be paid.
The count margin must be explicit:

```math
2^{h+m_0}=N2^{-(M-m_0)}.
```

Fix the coefficient of $`m_0=c_0\ell`$ before increasing A, or reserve
$`M-m_0\ge A_{\rm slack}\ell`$ with a sufficiently large independent
count margin. The same wires cannot simultaneously enlarge the local
table and supply an unchanged count suppression. Extra protected wires
may remain unused; do not enlarge the donor bank automatically with M.

Start with a concrete algebraic opening: let $`C=U_q\cdots U_1`$ be
the ordinary prefix-controlled rotation word, where the jth rotation
angle is selected by the already updated prefix. Derive its row identity

```math
C[y,x]=\prod_{j=1}^q
 R\bigl(\theta_j(u y_{\lt j})\bigr)[y_j,x_j],
```

where R is the real two-by-two rotation. Thus the row state
$`C^\dagger|y\rangle`$ is a product state conditional on an unchanged
y register. Verify the full flag action and literal phases first, with
three independent angles at height two as the smallest heterogeneous
diagnostic. The sequential word for C is not the desired shallow circuit.

A selected row source would have the action

```math
|y\rangle|0^q\rangle\longmapsto
|y\rangle C^\dagger|y\rangle.
```

Its prefix/row/leaf address has a+2q bits. Check the existing family-source
interface and its small fresh source bank before inventing another state
preparer. The missing step is coherent removal of the old y and any
program history, possibly only on the reached tagged bands. Preparing
orthogonal rows does not remove that record. Ordinary query uncomputation
after changing its address is invalid. A positive-amplitude overlap
encoding must also pay its normalization; even the balanced q-bit
rotation tensor product has entrywise-magnitude norm $`2^{q/2}`$.

More generally, seek a joint band character, a larger source encoding,
or a reversible history using this capacity and the clean bank. One
stronger target, after shifting the absolute tag by a, is

```math
\mathcal Z=
\sum_{z\in\{0,1\}^s}|z\rangle\langle z|_T\otimes
\sum_{j=0}^q(-1)^{z\cdot\mathrm{bin}(j)}\Pi_j^{u,q},
\qquad s=\lceil\log_2(q+1)\rceil.
```

Hadamards around this character give a full local band-label XOR.
Its formula supplies no implementation. A useful first test is a
non-threshold character such as depth parity on a genuinely heterogeneous
tree, with the clean bank and its actual return made explicit. A direct
reached-tag eraser is acceptable even if it never builds this character.

The first lane deliverable is a native encoding or chronological word
with its table exponent, success amplitude or normalization, correction
depth, and full return contract. Count capacity cannot rescue exponential
amplitude-amplification or serial-repetition depth. A frame-conjugated
computational phase is circular unless the conjugation has a new paid
implementation. Do not copy an unknown quantum data state, assume a
fresh source for each character, or infer parallelism from commutativity.

## 5. Resource and composition gates

Price the **coherently selected family of local words**, including all
bins and prefixes. Per-bin estimates alone are insufficient. The retained
source, unloading/reloading, and coarse recovery contribute
$`O(\ell\log(\ell+2))`$ depth in R85's schedule. A sufficient new
local-family envelope, at the final allocated error, is

```math
T_{\rm loc},w_{\rm loc}
 \le\mathrm{poly}(n)(\sqrt Q+1),\qquad
G_{\rm loc}\le\mathrm{poly}(n)(Q+1),
```

```math
D_{T,\rm loc}=o(\ell^2),
```

with a fixed polynomial degree independent of the later choice of A,
and with all simultaneous work counted.
The preferred local depth is $`O(\ell\log(\ell+2))`$. Once the
degree and clean allocations are known, a sufficiently large constant
A suppresses those polynomial count and dirty-width factors through
$`Q=N2^{-M}`$. For a richer table of size $`2^{h+m_0}`$, use the
actual residual margin $`M-m_0`$ from Section 4 rather than absorbing
$`2^{m_0}`$ into a degree that silently grows with A.
Explicitly verify that the new initialized work fits
the exposed reserve; the polynomial envelope does not make it clean.

For a direct local eraser of error epsilon on the whole reached image,
the early comparison can take the form

```math
\epsilon_{\rm early}
\le\delta_0+2\delta_U+\epsilon+8k\delta_R,
\qquad k=\lceil\log_2p\rceil,
```

where p is the actual bin count. This uses the R85 source and coarse
recovery interfaces only if their controls, zero-work projectors, and
actual-inverse premises still hold. Derive any changed bound afresh.
Orthogonal preserved-prefix and bin sectors take a maximum; serial
stages require an error sum. Increase source precision accordingly and
reinsert it into counts, depth, and clean-bank requirements.

The default terminal is the existing R81 construction with
$`M=O(\ell)`$. A larger protected bank is a legitimate reserve:
R81 prices every M, and permits a sub-square-logarithmic terminal when
$`\log(M+2)=o(\ell)`$. If choosing such an M, recompute bin count,
all table exponents, routing, source work, and terminal depth together.
In particular, $`2^{m_0}`$ is no longer automatically polynomial in n
when the free reserve grows beyond logarithmic size.

Making every bin small can destroy the gain: if $`p=n^{\Omega(1)}`$,
binary search needs $`\Omega(\log n)`$ full-address comparisons.
Polylogarithmically many bins can still permit the desired saving. Requiring
superlinear clean work in q changes the geometric capacity and may
cause exactly that problem. Count the new bins and selected routing;
do not inherit R85's estimates unchanged. A single coarse-source-layer
improvement for the R85 family is a component result until the separate
terminal and routing costs are also removed.

## 6. Five-hour execution and reformulation

Use a coordinator, two independent constructors, a native-resource
auditor, and a skeptical reviewer; retain capacity for a fresh proof
review or a narrow primary-source check. Constructors exchange both
identities and strongest objections at the first gate.

| Elapsed minutes | Work and decision |
|---|---|
| 0–25 | Verify the baseline; reconstruct the general unloading premise and freeze the reached-tag contract and live registers. |
| 25–85 | Develop lanes A and B independently. Produce exact algebra, provisional recurrences, and the first falsifiable lemma. |
| 85–125 | Exchange objections; deliberately reformulate the encoding, the location of history, or the bin/workspace tradeoff. Choose one primary mechanism and at most one distinct reserve. |
| 125–215 | Concentrate on the decisive native identity, maintaining its simultaneous resource, error, and return ledger. Switch only for a stated mathematical reason. |
| 215–255 | Reconstruct the complete argument or isolate the exact remaining lemma. Start fresh review and necessary validation once the statement is stable. |
| 255–300 | Finish validation and integration within the ceiling; preserve a precise resumption question and honest theorem/check status. |

An unresolved but narrowing construction is allowed to use the discovery
budget. A finished compiler is not required at the first gate. Stop a
candidate whose complete recurrence is unchanged, whose clean bank is
unpaid, or whose return identity fails; use the reformulation block to
change that premise rather than terminating the general objective.

The following are specific audit questions:

- Where is the first operation that actually uses the newly exposed clean
  wires to improve depth? Would the proposal be identical without them?
- Does the selected table include the preserved prefix and every local
  row/column/history index? Which wires contain parked dirty data, and
  why are they absent from addresses that assume zeros?
- Are parallel query outputs using disjoint non-Clifford supports, and
  are all simultaneously live helper pools included in width?
- Does a recursive or encoded operator merely rename the original frame
  compiler? Is its claimed saving still present after all inverses and
  source precisions are counted?
- Does the full proof cover normalized marker columns below tiny incoming
  amplitudes, singular pairs, occupied inactive records, and references?
  State-preparation closeness or clipping rare branches is insufficient.

Only reopen the secondary logarithmic-depth route if it has a concrete
new complete schedule or a clearly scoped useful lemma. Do not let an
already solved subclass replace the general-magnitude destination.

## 7. Evidence, output, and integration

Keep one compact working checkpoint: target, exact identity or conjecture,
changed resource, register map, recurrence, error/return argument,
strongest objection, and the next derivation. Preserve expensive reasoning
on the research branch so a session boundary does not lose it.

Write symbolic identities and native ledgers before numerical fixtures.
Use a few exact checks only for a surviving mechanism; do not run broad
simulations or infer asymptotics from finite fits. Search literature only
for a named missing primitive, using primary sources and translating its
initialization, access, precision, normalization, and return assumptions.
No external theorem is imported by this design.

| Outcome | What to retain |
|---|---|
| Complete general target | The full fixed-accuracy circuit and simultaneous resource/error/return proof, with a finite fallback and independently reconstructed pivotal lemma |
| Proved mechanism | Its exact scope, new native saving, and remaining full-circuit obligations |
| Unresolved mechanism | New algebra or a materially sharper interface, one precise missing lemma, and a concrete resumption question |
| Candidate failure | A scoped counterexample or recurrence obstruction and the premise needed to reopen it |

A fresh reviewer must reconstruct the pivotal lemma before promoting a
theorem. An inherited defect takes priority over building on it. Existing
R79–R85 proofs remain authoritative until a proved correction or extension
is integrated; create no result label for a plan or conjecture. Keep one
proof home and update frontier pages only when a proved bound changes.

Record the execution outcome separately and mark this prompt executed.
Use the [workspace validation requirements](../../WORKSPACE.md#restore-and-verify)
for substantive mathematics or circuit changes; document-only outcomes
need the relevant narrative, link, mathematics, and rendered checks.
Merge only after required CI passes on the reviewed tree. The five-hour
clock includes waits and integration: start long checks early, and
preserve the branch with explicit pending gates if the ceiling arrives.
Report actual elapsed time. No final manuscript drafting or external
reviewer outreach is part of this session.
