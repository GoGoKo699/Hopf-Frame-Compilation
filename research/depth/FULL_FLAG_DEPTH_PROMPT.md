# Five-hour design: read the complete flag and pay precision once

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[Session baseline R94](PORTABLE_WORKSPACE_CATALOGUE.md) ·
[Previous session](CATALOGUE_BOOTSTRAP_SESSION_RESULT.md) ·
[Causal flag construction R89](CAUSAL_COMPACT_COIN_DEPTH.md)

**Status: executed; [outcome and validation](FULL_FLAG_DEPTH_SESSION_RESULT.md).**
Prepared on 16 October 2026 after
[PR #158](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/158),
from scientific baseline main commit
`6302418330eef6901b06154c5ccb0aeb9386a286`.
Execution began on 16 October 2026 at 09:39:11 UTC from main
`5833a19fc0faf58b74c7795664636d2cb3e28a77`, with ceiling 14:39:11 UTC.
The original design below is retained as the work order. Its supported
outcome is [R95](LOGARITHMIC_PRECISION_DEPTH.md); proposed routes below
are not independently promoted to theorems.

Work only in **GoGoKo699/Hopf-Frame-Compilation**.
Repository modification and merge remain authorized. Final manuscript
drafting is separate. Optimal general depth remains the destination;
even the preferred result below would be an upper bound, not a proof of
optimality.

## 1. The next target

Put $`N=2^n`$, $`\ell=\log_2(n+2)`$, and $`\eta=2^{-L}`$.
Start with the real-frame benchmark

```math
L=\max\{6,\lceil\ell\rceil\},\qquad
b=\lceil C_b\sqrt{NL}\rceil,
```

for a sufficiently large absolute constant and sufficiently large n.
Seek one complete native
circuit with

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad D_T=O(\ell).
```

R94 gives $`O(\ell\log_2(\ell+2))`$ here. The new question is
whether its extra logarithm can be removed by changing the decoder and
the local synthesis architecture. The accuracy is inverse polynomial
in the logical qubit count n, not in N.

If the benchmark closes, seek the uniform real and complex magnitude
statement on $`6\le L\le2\ell`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\ell).
```

Require one external initialized qubit, the original real width
$`b\ge17(L+n+7)`$ or complex magnitude width
$`b\ge17(L+n+8)`$, and

```math
\|\mathcal VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L}.
```

The error includes all columns, signs, zeros, literal phases, clean-work
leakage, arbitrary dirty inputs and external references. Use actual
native inverse words. Every source, query, arithmetic step, guard,
recovery, workspace wire and elementary Clifford gate is charged. There
is no measurement, reset, QRAM, supplied precision state, or extra
external initialized bank. Classical preprocessing must terminate;
efficient preprocessing and a scalable emitter are separate questions.

## 2. Two costs must change together

| Part of R94 | Present cost at the benchmark | Required change |
|---|---|---|
| Opening after the first R89 reduction | R93 Section 3 reduces $`m=O(\ell^3)`$ to $`O(\ell)`$ using $`O(\log\ell)`$ sequential search/source positions; each pays the long-prefix query depth $`O(\ell)`$ | A complete reduction to $`O(\ell)`$ remaining layers in $`O(\ell)`$ depth |
| Native word inside each catalogue slot | R79 at group height $`h=\Theta(\ell)`$ and final precision $`\kappa=\Theta(\ell)`$ costs $`O(\ell\log\ell)`$ | A portable one-clean word of depth $`O(h+\kappa)`$, with its full workspace charged |

Improving only one row leaves the complete R94 asymptotic order
unchanged. R91's parallel word query and paired router already fit
$`O(\ell)`$; another routing optimization does not address these costs.

The useful new freedom is local: **a catalogue entry does not need
locally optimal gate count.** Its much larger count and arbitrary dirty
workspace can still fit below the global budget. Its initialized-work
promise cannot be enlarged for free.

Run two constructive lanes and one independent skeptical lane. A single
new construction may replace both interfaces, but must demonstrate both
savings in the complete circuit.

## 3. Opening lane: expose every spectral label at once

The starting objects are R78's nested subtree spaces and their
orthogonal difference projectors $`\Pi_t`$. R89 already exposes coarse
flags through compact coins, causal matching and reference transforms.
It still needs fine recovery and $`O(u^3)`$ protected coordinates,
where $`u=l+\lceil\log_2(m+l+2)\rceil`$.

For a remaining positive rounded m-bit tree, reserve its bottom r
original logical bits and put $`a=m-r`$. The unchanged external
prefix P has $`n-m`$ bits. Seek a full flag dilation for the early
a-bit frame $`W_{P,a}`$, with bands $`0\le t\le a`$;
the protected r-bit suffix supplies conditional-zero work on its saved
active sector. It is not an additional initialized allocation.
A sufficient interface consists of a native embedding $`E_P`$,
a fully charged reference transform R, and
a classical register $`\lambda`$ satisfying, on the ideal initialized
sector,

```math
R E_P\,\mathrm{Ran}(\Pi_t)
\subseteq\mathrm{Ran}[\lambda=t]
\qquad\text{for every }t.
```

Initialization is implicit in this display and must be explicit in the
proof. The image may contain additional work. With a tag T outside the
embedding, the complete chronological word is

```math
E_P,\quad R,\quad
(T\mathrel{\oplus}=\lambda),\quad R^\dagger,\quad E_P^\dagger.
```

It should implement the whole spectral-label XOR and return its work.
Preserving the root state alone is insufficient: every difference band
must have its correct label, and the completed readout word must have
the prescribed literal phase. Intermediate band-dependent phases in
the embedding may cancel through its actual inverse. A source completion chosen
arbitrarily on the other columns does not establish this property.

The concrete first candidates are one joint compact-coin representation
or a compressed causal record preserving the complete nested flag. Audit
R89's metadata, rank windows, matching selectors, reference banks and
fine decoder separately. Distinguish information stored in genuinely
initialized coordinates from data that can be moved through arbitrary
dirty work in a completed reversible word. R94's slot cancellation is a
tool for moving a unitary; it does not initialize an occupied record.

The desired additional conditional-zero reserve is
$`O(l+\log_2(m+l+2))`$, with a constant number of sequential
P-dependent query positions. Any P-dependent reference data are queried
and charged as well. A useful sufficient table-address target is

```math
n-r+C(l+\log_2(m+l+2)),
```

where r is the protected suffix. Give output lengths, instance counts,
query minimum widths and standalone costs before choosing
$`r=A(l+\log_2(m+l+2))`$. Constants and polynomial degrees must
be fixed independently of A. A different ledger is admissible if it
proves the same complete budgets.

Apply the resulting early compiler after R94's first uniform R89 stage,
at $`m=O(\ell^3)`$ and $`l=L+O(1)`$. It must leave
$`O(L+\log\ell)=O(\ell)`$ layers, not another cubic reserve.
Include the initial tagged source and fine labels in the claim; a coarse
label improvement alone does not close this lane.

Approximate embeddings need an operator-isometry error bound on the
entire initialized input space. Establish the comparison word and actual
inverse estimate before using a hybrid. Exact comparison bands may
commute while their native approximations do not. Saved-H inactive
identity must hold on arbitrary occupied work; the final predicate
inverse must propagate leakage without reset.

## 4. Local lane: spend dirty space to save depth

For each rounded h-bit group $`W_v`$, retain R94's controlled Hopf
frame $`A_v=\mathrm{diag}(I,W_v)`$. Seek a physical native word
$`U_v`$ with one initialized flag d and an arbitrary pool $`D_*`$:

```math
\|U_vJ_d-J_d(A_v\otimes I_{D_*})\|\le2^{-\kappa},
\qquad d_{\rm loc}=O(h+\kappa).
```

The important regime is $`h=\Theta(\ell)`$,
$`\kappa=\Theta(\ell)`$. Any short group must also have a
valid bound. A sufficient, deliberately generous local allowance is

```math
t_{\rm loc},\ g_{\rm loc},\ w_{\rm loc}
\le 2^{P(h+\kappa+1)},
```

for one fixed polynomial P with absolute coefficients. This is an
allowance for arbitrary dirty work, not initialized ancillary bits.
Prove a terminating choice of a physical word for every rounded code.

Retain $`h\le\lfloor\ell/8\rfloor`$ beyond a fixed cutoff and
R94's complete code $`q=O(\kappa2^h)=o(n)`$. Then

```math
\log_2\bigl(w_{\rm loc}2^{2q}
       +2^q(t_{\rm loc}+g_{\rm loc})+1\bigr)=o(n)
```

is the absorption test for the whole catalogue and paired router. This
permits a radically different local compiler without weakening the
global optimal T-count or Clifford budget. Count and depth remain
separate ledgers. All slot operations run in parallel only because their
complete payloads and work pools are disjoint.

Inspect the [shared-gradient construction](SHARED_GRADIENT_RECOVERY.md)
as a possible precision source. Its invariant-sector comparison and
actual exit inverse are valid tools. Its bank, however, contains
$`3r+q+s+2=O(\kappa+\log(h+2))`$ initialized coordinates.
These are not supplied by one clean flag and an arbitrary dirty pool.
At this local relative precision R79 can use its R72 fallback; improving
only R79's early source leaves its deepest precision-paying layers.
A proposed shared source must close that bottom and its initialization
cycle explicitly.

Once a local word qualifies, use R94's complete payload
$`S=(c,X,d,D_*)`$, exact transfer, and commutator with $`X_c`$.
The unwanted slot action then cancels on its full occupied space. Retain
the factor-two local error estimate, exact identity code, guard excluding
changing payload bits, freshly read prefix, and full leakage telescope.
Do not substitute a partial payload or an old X-only residual correction.

## 5. Early falsification gates

Reject the following shortcuts unless their missing premise is proved:

- Applying R94 to its own local entries: the needed
  $`\kappa=\Theta(\ell)`$ exceeds the local theorem range
  $`O(\log(h+2))=O(\log\ell)`$. Any recursion needs valid
  domains and a decreasing dependency, not a circular appeal.
- Iterating unchanged R89 at final precision: its reserve
  $`A(l+\log m)^3=\Theta(\ell^3)`$ need not contract.
- Using the R79 high-precision clause by shrinking h: making
  $`\kappa\gtrsim2^h/h^2`$ affordable restores a growing number
  of groups and their long-prefix reads.
- Cataloguing the entire cubic-height remainder: the complete program
  has $`q=O(L2^m)`$ bits; its $`2^q`$ catalogue is not covered
  by the affordable $`q=o(n)`$ argument.
- Reading a whole program once and executing it from arbitrary dirty
  storage: an XOR lookup produces an unknown XOR mask. Noncommuting
  program execution does not automatically cancel that mask. Require a
  literal full-space identity with native costs.
- Parallel threshold tests or source calls merely because their ideal
  projectors commute. Shared logical targets and changing addresses need
  a physical schedule; the retained threshold-product obstruction in
  [the joint-depth record](JOINT_DEPTH_SESSION_RESULT.md#4-threshold-products-and-arbitrary-source-completions)
  must be checked before reusing that architecture.
- Moving metadata into dirty registers without proving occupied-input
  correctness, or inferring a complete flag map from a root-state source.

Use tiny noncommuting examples to reject candidate identities after
their symbolic words are written. Useful risks are two same-root
completions with different marker bands, occupied precision records,
changed selectors, and a dirty-mask echo with noncommuting instructions.
Bounded fixtures test these identities; they cannot prove asymptotic
resource bounds or an unrestricted depth lower bound.

## 6. Five-hour execution and decision gates

The ceiling includes review, validation, rendering, CI and integration.
It is a maximum budget; completion can be earlier.

| Elapsed time | Required work and decision |
|---|---|
| 0–35 minutes | Reconstruct the two contracts and register lifetimes. The skeptical lane checks local precision domains, initialized work and the old threshold/cache counterexamples. |
| 35–100 minutes | Opening lane develops a complete flag representation; local lane develops a depth-first compiler under the relaxed dirty budget. By minute 100 require explicit native words or a concrete synthesis recurrence for both, or one complete bypass. Favorable isolated estimates do not pass. |
| 100–180 minutes | Close surviving operator, error and resource arguments. If only one lane survives, finish a reusable component at its actual scope and identify the precise missing partner. Do not claim a global improvement. |
| 180–235 minutes | Freeze the construction. Independent reconstruction checks every source, inverse, guard, table address, simultaneous pool and error term. Add only bounded fixtures needed to resolve concrete risks. |
| 235–300 minutes | Integrate the supported result or scoped obstruction; run required checks and desktop/narrow rendering, publish the reviewed tree and complete authorized integration after CI passes. |

The full theorem must retain all R94 obligations: a fixed number of
groups, geometric query sums, exclusive query pools, one absolute cutoff,
and the original-width narrow/finite fallback. Prove the internal real
precision margin before applying the complex L+1 instance and paid
diagonal. Use the original-precision complete fallback when a stricter
tolerance does not fit its literal allocation. Do not give the complex
extension or all-width theorem merely on the strength of the benchmark.

If a principal route fails, record its exact surviving operator or
unpaid charge. Do not spend the remaining time tuning unchanged
asymptotics. A broader-precision corollary of R94 can be considered only
after the main decision gate and with a complete uniform ledger; it is
a distinct outcome, not evidence that the logarithmic-depth target was
reached. A lower-bound pivot requires a new invariant for the actual
Clifford+T model. Finishing a constant-layer obstruction alone would not
establish dimension-growing optimal depth.

## 7. Reviewable outcomes

The preferred outcome is a complete stronger theorem with all resources
on the same circuit. A useful narrower outcome is a proved full flag
decoder or a one-clean local depth lemma with its exact precision and
workspace range. A failed construction can yield a precise counteridentity
and reopening condition. None of these outcomes ends the general search.

Keep R94 and the retained R79/R90/R92 frontiers intact unless a new proof
supports a precisely scoped improvement. Selected manuscript Result D
remains R79. Assign a new result identifier only to an established claim;
this design allocates none. The session report should distinguish the
proved theorem, component results, interface failures, remaining gap,
independent reviews and finite checks.
