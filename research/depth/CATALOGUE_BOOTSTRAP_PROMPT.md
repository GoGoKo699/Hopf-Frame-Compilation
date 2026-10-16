# Five-hour design: carry the workspace through the catalogue

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[Current grouped theorem R93](GROUPED_CATALOGUE_PRECISION.md) ·
[Previous session](PRECISION_TERMINAL_SESSION_RESULT.md) ·
[Local compiler R79](UNIFORM_SUBTREE_DEPTH.md)

**Status: executed.** The [session result](CATALOGUE_BOOTSTRAP_SESSION_RESULT.md)
records R94's portable-workspace construction, which reaches the primary
benchmark below. The session began on 16 October 2026 at **07:56:34 UTC**,
with a five-hour ceiling of **12:56:34 UTC**, including review and integration.
The design was prepared on 16 October 2026 after
[PR #156](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/156),
with scientific baseline main commit
`b07e5902cfca633012677b5a3ce78b29b781c052`.
Planning and preliminary algebraic reviews alone do not establish the
proposed theorem. The separate [R94 proof](PORTABLE_WORKSPACE_CATALOGUE.md)
and independent audit support the new result; older theorem statements
remain retained.

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge remain authorized. This is research, separate from final
manuscript drafting. Optimal general depth remains the destination.
The immediate question is whether the established one-clean compiler
can supply its own faster catalogue entries when their workspace travels
with the selected logical group.

## 1. Target and unchanged contract

Put $`N=2^n`$ and $`\ell=\log_2(n+2)`$. The primary benchmark is

```math
L=\max\{6,\lceil\ell\rceil\},\qquad
\eta=2^{-L},\qquad b=\lceil C_b\sqrt{NL}\rceil,
```

for a sufficiently large absolute constant. Seek a complete arbitrary
real-frame compiler with absolute constants and

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad
D_T=O\!\left(\ell\log_2(\ell+2)\right).
```

R93 gives $`O(\ell^2/\log_2(\ell+2))`$ at this benchmark.
The requested improvement is a mathematical target, not a proved
consequence of the plan. It would still not establish optimal depth.
The error is inverse polynomial in the logical qubit count n.

Require one external clean qubit and the full initialized-isometry bound

```math
\|\mathcal VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L}.
```

Every logical column, literal phase, sign, zero, clean-work leakage,
dirty input and external reference is included. Use actual inverse
native words. There is no measurement, reset, QRAM, supplied precision
state, extra initialized word, or uncharged coherent evaluator.
Elementary Clifford gates remain charged. Classical preprocessing must
terminate; efficiency is not required by the retained model.

Once the real benchmark closes, seek the same construction throughout
R93's range, including the paid complex phase completion:

```math
6\le L\le2\ell,\qquad
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
```

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right).
```

Preserve the literal widths $`b\ge17(L+n+7)`$ for real frames and
$`b\ge17(L+n+8)`$ for complex magnitude frames. Prove the necessary
real precision shifts, for example internally through $`L\le3\ell`$,
before invoking the complex result at L+1. A benchmark-only or real-only
result must be stated at that scope. Do not pursue a stronger additive
$`O(\ell)`$ target before this construction is complete.

## 2. The changed assumption

R93 uses generic zero-clean local factors, forcing short groups to
control their depth. Its fine correction then needs initialized endpoint
labels. This design changes both costs at once: compile each entire group
directly to final accuracy with R79, and carry the compiler's workspace
through the dirty catalogue transfer.

The local compiler requires a clean input only when acting on the
logical payload. Its copies in borrowed catalogue slots may receive
arbitrary occupied inputs. They need only be fixed physical unitaries
whose unwanted slot action cancels exactly. That cancellation must be
proved for the complete payload, including every local dirty wire.

This is an acyclic use of an established theorem: the new global circuit
calls the retained R79 compiler at small dimension. It must not invoke
the proposed improved theorem to compile its own entries.

**Compile directly to final accuracy.** Do not attach R93's old residual
source to an approximate enlarged payload. The actual coarse word can
leak into its selector and initialized flag and act on its dirty pool;
its residual is then not an operator on the h-bit group alone. The old
matrix-unit correction does not apply to that residual. The present
route removes this correction, its endpoint labels, and its history.

## 3. First gate: the complete transfer word

For a rounded h-bit real group frame $`W_v`$, use one selector c and form

```math
A_v=\mathrm{diag}(I,W_v),\qquad B=X_c.
```

First verify directly from the prescribed Hopf completion that $`A_v`$
is an (h+1)-bit frame: zero root angle, identity left subtree, and
$`W_v`$ in the right subtree. Do not add a generic control to a native
R79 circuit. Its Clifford interlayers could become costly controlled
operations, destroying the claimed T-depth.

Compile this controlled frame with the already-proved real R79 theorem.
Let d be its initialized flag, $`D_*`$ its arbitrary local dirty pool,
and $`\widehat A_v`$ the resulting actual native unitary. The entire
payload is

```math
S=(c,X,d,D_*),\qquad w=h+2+b_{\rm loc}.
```

Adapt R90 Section 1's selected swap to move **all** of S into a slot.
For the corresponding transfer word, prove the exact full-space identity

```math
K_v=\widehat A_v(S)\otimes
       \widehat A_v^\dagger(\mathcal A_v)
       \otimes I_{\rm other\ slots}.
```

Since B acts only on the logical selector, the proposed cancellation is

```math
K_vB_SK_v^\dagger B_S^\dagger
=\bigl(\widehat A_vB\widehat A_v^\dagger B^\dagger\bigr)_S
   \otimes I_{\mathcal A}.
```

One borrowed catalogue array should suffice for this word: B is a local
Clifford operation and is not transferred from another slot. This is a
different use from the invalid merger of R90's two transferred factors
onto one array. The auxiliary paired-router array is still charged.

State that rightmost factors act first, and give the chronological native
word separately. Verify the branch direction:

```math
A_vBA_v^\dagger B^\dagger
=\mathrm{diag}(W_v^\dagger,W_v).
```

Thus c=1 must implement $`W_v`$, with d initially zero. Derive the inverse
initialized-isometry estimate from the forward estimate, using the actual
inverse word. Because B preserves d=0, a two-term unitary hybrid should
bound the commutator error by twice the local compiler error on arbitrary
$`c,X,D_*`$ inputs with d=0. No assertion about clean flags in occupied
array slots is permitted. Their return must be exact on the whole space.

Use [R91's paired router](PARALLEL_BILINEAR_WORD_QUERY.md#6-a-faster-selected-slot-router-with-a-larger-dirty-allocation)
and its zero-clean word query. The paired translation temporarily
permutes its auxiliary array E; only the complete selected swap returns
E. Include both the dirty-program echo and the complete swap when
claiming work return. The reserved zero program must select the actual
identity native word, not merely an approximation to identity.

### Bounded falsification fixtures

After expanding the symbolic word, test noncommuting selected factors,
occupied slot flags, arbitrary dirty pools, and literal phases. Include
a four-mode frame where an earlier rotation changes a later selector.
Use a negative control that leaves $`D_*`$ stationary while swapping only
part of S: choose a unitary whose initialized action is correct but whose
occupied-flag sector couples to that pool. This should expose why the
whole payload must travel. Check the c-branch direction and the actual
inverse, rather than replacing the inverse by a second approximation.

Small operator fixtures can reject an identity. They cannot establish
the asymptotic theorem, the full native emitter, or the width bound.
Avoid large simulations and tests that only repeat implementation steps.

## 4. Final-accuracy groups and simultaneous native resources

Reuse R93 Sections 2–3: one uniform R89 reduction followed by the
prefix-selected R78/R79 early decoder leaves

```math
M=O(L+\log_2(\ell+2))=O(\ell)
```

through the intended precision range, with the already-priced
$`O(\ell\log_2(\ell+2))`$ opening depth. Do not replace the long
preserved prefix by a short address in any query estimate.

Choose one fixed absolute grouping constant, initially

```math
h_{\max}=\max\{1,\lfloor\ell/8\rfloor\}.
```

Leave a fixed bottom reserve Z and partition the remaining height into
groups of at most this size. Prove that the number R of groups is bounded
by one absolute constant. Empty middles and a final short group need
explicit treatment. At each group, use final precision

```math
\kappa=L+\lceil\log_2(R+1)\rceil+C_\epsilon,
```

with a fixed error margin covering angle rounding, the two local compiler
occurrences, the opening decoder, and the bottom. For its (h+1)-bit
controlled frame, choose a sufficiently large absolute local constant and
reserve

```math
b_{\rm loc}=\left\lceil\max\left\{17(\kappa+h+8),
                 C_{\rm loc}\sqrt{2^{h+1}\kappa}\right\}\right\rceil.
```

The required R79 specialization has separate native bounds

```math
t_{\rm loc}=O\!\left(\sqrt{2^{h+1}\kappa}+\kappa\right),
\qquad g_{\rm loc}=O(2^{h+1}\kappa),
```

```math
d_{\rm loc}=O\!\left((\kappa+\log_2(h+3))
                          \log_2(h+\kappa+3)\right)
            =O\!\left(\ell\log_2(\ell+2)\right).
```

Reconstruct R90's depth-weighted angle encoding at final accuracy. A
sufficient complete program bound is

```math
q_h\le C_q(\kappa+1)2^h=O(\ell(n+2)^{1/8})=o(n),
\qquad V_h=2^{q_h}.
```

Price a fixed local native word for every code, including unused codes,
and make code zero exactly identity. The crucial room is exponential in
n: the catalogue and router are $`2^{o(n)}`$, or $`N^{o(1)}`$.
Do not call these quantities subpolynomial in n.

**Separate T-count, Clifford count, and depth.** R91 Section 6's displayed
catalogue estimate uses one S to bound all three local costs. Reusing
that maximum here would erase the improvement. Write the refined ledger
from parallel circuits on disjoint slots. With Q prefix addresses and B
exclusive dirty query wires, sufficient bounds for one group are

```math
T_i=O\!\left(\sqrt{Qq_h}+Qq_h/B
                +V_ht_{\rm loc}+wV_h^2\right),
```

```math
G_i=O\!\left(Qq_h+V_hg_{\rm loc}+wV_h^2\right),
```

```math
D_i=O\!\left(\log_2(n+2)+Qq_h/B^2
                +d_{\rm loc}+\log_2(q_h+2)\right).
```

Every simultaneous pool must appear in an explicit wire table. A
conservative reservation may use

```math
q_h+3wV_h+C_{\rm router}wV_h^2+B+b_{\rm loc}
       +\text{predicate helpers}.
```

The query pool must be disjoint from the logical payload, all slots,
program register, router grids, public address, and protected logical
reserve. No local compiler may quietly borrow another slot's workspace.
Reusing returned pools between groups is allowed; overlap within a
native layer needs a proved schedule.

For s lower logical bits, prove

```math
Q\le C2^{n-s-h},\qquad Qq_h\le CN\kappa2^{-s}.
```

Sum query volumes and their square roots geometrically. Constantly many
local depths and long-prefix reads then meet the primary target if all
the other obligations hold. A short final group must satisfy the same
ledger; do not charge it as a fresh length-n endpoint table.

## 5. Logical guards, approximation, and the original width

Choose a fixed bottom reserve Z with distinct locations for c, d, and
one group guard G, plus any explicitly charged constant helpers. Save
$`H=[Z=0]`$ in the sole external clean flag before changing Z. Prepare
c=1 by a CNOT from H; d is initially zero on that sector.

For each group, compute G with a direct H control. Its lower-bit
predicate must exclude Z and every bit in the current payload X.
Select the rounded group code precisely on H=G=1 and select zero
otherwise. Require exact full-space identity for the zero code, including
occupied logical c,d and arbitrary dirty arrays. The inactive-sector
proof may not assume those logical bits are clean.

The prefix P is preserved during its group and is freshly read for the
next group. It may contain outputs changed by earlier groups. Do not
freeze its initial value in a retained program. G must uncompute exactly
because its predicate excludes every changed payload bit, even when the
local compiler leaks into c or d. Give this argument as a full operator
statement; an ideal-path diagram is insufficient.

Compare the composed actual word with an ideal sequence that returns
c,d and $`D_*`$ after every group. Telescope the initialized-isometry
errors against those ideal states. Never reset a leaked flag or assume
the next actual input is exactly initialized. Undo c preparation and the
saved-H predicate with their actual inverse words, retaining all leakage
in the norm bound. Finish the constant-height bottom using the retained
smaller-error R72 construction, with its long-prefix and native costs.

For the all-width extension, split at $`b=\sqrt N/n`$. In the wide branch,
prove one absolute dimension cutoff, uniform through the required
precision shifts, beyond which

```math
wV_h^2+V_hg_{\rm loc}+b_{\rm loc}+q_h
=2^{o(n)}\le2^{n/4}.
```

After including predicate work and all constants, leave a fixed fraction
of b for an exclusive global-query pool, and check its native minimum
$`B\ge C(n+2)`$. In the narrow branch, use the complete retained R79
circuit at the original precision and literal width:
$`NL/b^2\gt n^2L`$ absorbs its depth. Handle finite dimensions by one
absolute cutoff; L is bounded in this range, so no L-dependent maximum
may be hidden in a uniform constant.

Finally price signs, zeros, and the prescribed complex diagonal. Use the
proved shifted real instance and the paid diagonal at half error within
the original complex width. Do not insert a stricter-precision fallback
at the real minimum width. Preserve R79's other regimes, R90's fixed-L
result, and R92's absolute envelope regardless of this session's outcome.

## 6. Five-hour execution and decision gates

Record the actual UTC start and five-hour ceiling in the eventual session
record. The ceiling includes tool time, review, rendering, CI, and
integration; it is a budget, not a requirement to remain busy. Use two
cooperating constructive lanes and a separate skeptical reviewer. One
lane owns the full payload transfer and local native catalogue; the
other owns guards, composition, and simultaneous global resources.

| Elapsed time | Work and required checkpoint |
|---|---|
| 0–40 minutes | Re-read the retained interfaces. Expand the complete word, controlled-frame embedding, inverse estimate, and occupied-input cancellation. Reject a branch or workspace error immediately. |
| 40–100 minutes | Derive final-accuracy code volume and separate native count/depth ledgers. Produce the complete wire reservation and fixed-bottom guard design. By 100 minutes require a credible simultaneous path to the target, not isolated favorable estimates. |
| 100–190 minutes | Close the global initialized-isometry proof, geometric sums, uniform cutoff, literal-width fallback, and complex precision shift. Stop extending the scope when a missing premise cannot be resolved within the remaining proof budget. |
| 190–240 minutes | Independent audit of the same circuit and every resource. Add only bounded fixtures addressing concrete algebraic risks; resolve review findings. No new construction after minute 190. |
| 240–300 minutes | Integrate the supported outcome, update only justified claim/index entries, run required checks and desktop/narrow rendering, and complete the authorized repository workflow. Reserve time for CI. |

If the complete transfer fails the first gate, name the exact surviving
operator and its hypothesis. If the native or width ledger fails the
second gate, stop this theorem attempt and identify the charge that
prevents it. Do not spend the rest of the session tuning constants around
a failed asymptotic inequality or reattaching the invalid old residual.
The remaining time is for a precise component result or obstruction and
its integration. These gates stop a candidate, not the general objective.

## 7. Reviewable outcomes

The preferred outcome is a separate proof of the target with one complete
circuit satisfying all counts, depth, width, and operator requirements.
Only then promote its precise scope in the frontier and source map.
Use a new result identifier only for an established claim; this design
does not allocate one or change the selected manuscript Result D.

A narrower valid outcome is a portable initialized-workspace catalogue
lemma: exact full-space borrowed-array return, initialized logical-input
approximation, actual inverse, and separate native resource bounds.
State the missing global obligation explicitly. If even that lemma
fails, retain the exact counteridentity and the condition needed to
reopen it. No search failure is an unrestricted lower bound.

The final session report should say what was proved, what was disproved
at a stated interface, what remains open, which checks passed, and where
the reviewed work was integrated. Keep the distinction between a
successful identity, a native component, and a complete frame theorem.
The next destination remains the general depth question.
