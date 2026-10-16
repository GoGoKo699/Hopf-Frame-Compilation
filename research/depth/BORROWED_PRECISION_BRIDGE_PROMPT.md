# Five-hour design: move precision into borrowed work

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[Scientific baseline R95](LOGARITHMIC_PRECISION_DEPTH.md) ·
[Previous outcome](FULL_FLAG_DEPTH_SESSION_RESULT.md) ·
[Small-workspace sources](SMALL_WORKSPACE_DESCENDANT_SOURCE.md)

**Status: designed, not executed.** Prepared on 16 October 2026 from
main commit `67b7d724854fc46310104c8404ae2b9fb895bd3e`, after
[PR #160](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/160).
The proposed statements below are research targets, not additional
theorems. No result identifier is assigned.

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
remain authorized. Final manuscript drafting is separate. The execution
clock begins when this design is run; record its actual starting commit
and a ceiling five hours later. The ceiling includes review, validation,
CI and integration and need not be exhausted after completion.

## 1. Main question and measurable target

Can the initialized precision reserve of a globally selected early
compiler be reduced from linear to logarithmic in precision, while its
complete native count remains optimal after assembly?

Put $`N=2^n`$, $`\ell=\log_2(n+2)`$ and $`\eta=2^{-L}`$.
Seek the following uniform extension of R95:

```math
6\le L\le n^{1/4},\qquad
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
```

```math
D_T=O(NL/b^2+L+\ell).
\qquad\text{Target, not established here.}
```

Require one external initialized qubit, real dirty width
$`b\ge17(L+n+7)`$ or complex magnitude width $`b\ge17(L+n+8)`$,
and the same full initialized-isometry guarantee

```math
\|\mathcal VJ_1-J_1(W\otimes I_b)\|_{\mathrm{op}}\le2^{-L}.
```

All constants must be absolute on this fixed range. Include every
logical column, literal phases, signed and zero coefficients, work
leakage, arbitrary dirty inputs and external references. Count native
sources, every query, elementary Clifford gate, simultaneous workspace,
guards, arithmetic, recovery and actual inverses. No measurement, reset,
QRAM, supplied precision state, extra external clean bank or free
coherent evaluator is permitted. Classical preprocessing must terminate.

The first diagnostic slice is $`L=\max\{6,\lfloor n^{1/4}\rfloor\}`$
at $`b=\lceil C_b\sqrt{NL}\rceil`$ for sufficiently large n.
R79 currently gives $`O(n^{1/4}\log n)`$ there; the target is
$`O(n^{1/4})`$. Its accuracy reaches stretched-exponential error in n,
$`2^{-\Theta(n^{1/4})}`$, beyond R95's inverse-polynomial range.
This remains an upper-bound target, not a claim of optimal depth.

## 2. Why R95 does not already prove this

| R95 dependency | What grows outside its stated range | Required response |
|---|---|---|
| First uniform R89 opening | Protected height $`O((L+\ell)^3)`$ | Keep its full uniform ledger; the chosen range leaves this height $`o(n)`$ |
| Shared-gradient coarse/sparse stage | Remaining height $`O(L)`$ | Apply a new globally priced small-workspace stage to shrink it further |
| Final-accuracy catalogue | Groups of height $`O(\ell)`$ would repeat an $`O(L)`$ local cost $`O(1+L/\ell)`$ times | Reach height $`O(\log L)`$ before grouping |
| Relaxed local lemma | Its allowance $`2^{P(h+\kappa+1)}`$ has an unspecified fixed polynomial P | Avoid relying on it at polynomial L; use R79's already proved high-precision clause on smaller groups |

Increasing the fixed coefficient in $`L\le C\ell`$ is not the main
milestone. Nor is a wider range with the unchanged repeated precision
cost. Every asymptotic estimate below is conditional on the new stage
contract being proved.

## 3. Primary lane: a globally priced small-workspace early stage

Start from the actual word in
[local compiler Sections 2–3](RELAXED_COUNT_LOCAL_DEPTH.md#2-an-early-stage-with-a-logarithmic-precision-register):
coarse mathematical frame, constant-normalization signed sparse correction,
small-workspace descendant sources, and exact compact reference coins.
Promote that word to an unchanged global prefix P of length n-m.
Its existing local $`2^{P(m+l+1)}`$ allowance does not establish the
required global resource bounds.

For stage precision l define

```math
K=m+l+2,\qquad
\chi=\lceil\log_2(l+2)\rceil+\lceil\log_2(m+2)\rceil^3,
\qquad r=\lceil A\chi\rceil.
```

Use the stage only when $`r\lt m/2`$. Its early logical operator is

```math
\sum_P|P\rangle\langle P|\otimes
\left(W_{P,m-r}\otimes\Pi_Z+
I\otimes(I-\Pi_Z)\right),\qquad
\Pi_Z=|0^r\rangle\langle0^r|.
```

Save the original predicate into the sole external clean bit H. The
bottom Z coordinates are conditional work on H one; the completed body
must be literal identity on H zero even when that work is occupied.
Restore the saved predicate through its actual inverse.

A sufficient new ledger would have a constant number of sequential
P-dependent table positions, at most a fixed polynomial in K simultaneous
instances and outputs, and

```math
\overline Q=2^{n-r+C\chi},\qquad
T_{\mathrm{query}}\le CK^d(\sqrt{\overline Q}+\overline Q/b+1),
```

```math
G_{\mathrm{query}}\le CK^d(\overline Q+1),\qquad
D_{T,\mathrm{stage}}\le
C\left[l+m+\log_2(n+\chi+2)+K^d\overline Q/b^2\right].
```

Price standalone T-count, Clifford count and dirty width separately by
$`2^{C(m+\chi)}`$ or another explicit bound that fits the application.
State the sufficient private-query minimum widths and every simultaneous
allocation. Fix C, d and all other degrees before choosing A.
The intended use has $`m=O(L)`$, $`l=L+O(1)`$.
There the exact-reference depth $`O(m)`$ is affordable.

Reconstruct, rather than assume, these points:

- Append P to every table that depends on it, including hash/shear and
  matching catalogues. Inspect all source-anchor and transported-reflection
  addresses; no hidden $`Cl`$ address term may survive.
- Use the borrowed real calibration for every ordinary source. Its
  precision core is dirty, excluded from initialized anchor tests.
  Apply source completion on the entire physical initialization domain,
  including arbitrary dirty inputs.
- Retain exact padded geometric references and rebuilt causal matchings.
  Include their private initialized scratch in the augmented nonzero tests.
  Use their unconditional inverse/controlled-label/forward sandwich for
  occupied-sector cancellation.
- Preserve the live root residual, separate failure coordinates,
  normalization-two recovery, literal sign and all-work reflection.
  Use actual native inverses and full-isometry error estimates.
- Give a physical schedule. Commuting ideal projectors or disjoint logical
  targets alone do not justify parallel source calls with shared controls.

**Qualification gate:** by minute 90, require the complete word,
initialized-work lifetimes, global address envelope, standalone costs,
minimum widths, and full error/return argument. An unpriced local
exponential allowance or an extra clean precision bank does not pass.
The separated dependence $`\chi=O(\log l+\log^3 m)`$ is essential:
the weaker $`O(\log^3(m+l))`$ reserve would stall the second
contraction at $`O(\log^3 L)`$. Derive hash, matching and router
metadata from the coarse geometry of height m; show that the
$`O(\log l)`$ initialized source index appears only constantly
many times. Do not inherit a generic log K bound for this purpose.
An alternative ledger is acceptable if it proves the complete target
budgets. Do not add broad simulations to compensate for a missing bound.

## 4. Assembly lane: an explicit finite precision bridge

Use R95 unchanged on $`L\le2\ell`$. On the overlapping new branch
$`\ell\le L\le n^{1/4}`$, investigate this finite sequence:

| Stage | Intended remaining logical height |
|---|---|
| Retained first uniform R89 stage | $`m_1=O((L+\ell)^3)`$ |
| R95 shared-gradient coarse/sparse stage | $`m_2=O(L)`$ |
| First new small-workspace stage | $`m_3=O(\log_2^3(L+2))`$ |
| Second new small-workspace stage | $`m_4=O(\log_2(L+2))`$ |

Reprice the first two native stages from R89's uniform stage ledger and
the coarse/sparse component's full parameter bounds. Their final R95
application statements cannot be invoked outside the proved range.
In particular, retain every exponential precision term in an address
or standalone allocation before using $`m_1=o(n)`$.

When a reserve test fails, omit that stage and prove the same next-height
bound directly from the failed inequality. Treat heights at most three
with the retained bottom compiler, without negative-height blocks.
Both new stages have input height $`O(L)`$, hence exact-reference depth
$`O(L)`$. There are only constantly many sequential stages.

At the final tail, use R94's full-payload catalogue on groups of height

```math
h\le g=\max\{1,\lfloor\tfrac14\log_2(L+2)\rfloor\}.
```

Its controlled entries have h+1 logical bits. At precision
$`\kappa=L+O(1)`$, verify the R79 high-precision condition
$`\kappa\gt2^{h+1}/(h+1)^2`$ beyond one absolute cutoff.
It gives local depth $`O(L)`$ with polynomial-in-L counts and work.
The tail height gives a fixed upper bound on the number of groups.
Choose error padding only after that bound and all recovery factors
are fixed.

Retain the actual final-accuracy code, all slots and both routing arrays.
The intended absorption estimates are

```math
m_1=O(n^{3/4})=o(n),\qquad
q_h=O(L2^h)=O(L^{5/4})=O(n^{5/16})=o(n).
```

These leave room for the earlier standalone costs and the complete
$`w2^{2q_h}`$ router bill at $`b\ge\sqrt N/n`$. Prove that all
simultaneous pools fit fixed fractions of b. Reconstruct R91's original
global-address query sums; no local height may replace an actual prefix.

For $`b\lt\sqrt N/n`$, verify the original-precision R79 fallback
against $`NL/b^2`$. For bounded n, use its original-width full compiler.
Prove an internal real range with enough slack, for example
$`L\le2n^{1/4}`$, before paying the complex L+1 call and literal
diagonal. The final theorem must concern one circuit satisfying all
resources together, at the original real and complex width thresholds.

## 5. Optional strengthening after the main gate

The exponent $`1/4`$ is chosen for explicit slack, not asserted optimal.
A justified range enlargement may be included only after the primary
construction's ledger is complete.

A separate short exploratory lane may seek faster references with the
same small initialized footprint: parallel small-workspace reference
sources, or a hybrid of the exact capped geometric source and an
approximate tensor-product geometric source. These might reduce the
reference allowance from $`O(m)`$ to $`O(l+\log m)`$ and bypass
the first R89 cubic reserve. Require a literal vector/error formula,
full-space transported-anchor comparison, shared-H schedule and actual
standalone table sizes. The primary bridge does not depend on this idea.

Reaching $`L=\Theta(n)`$ through this constant-group assembly needs
another change. R94's actual code has

```math
q_h=\sum_{j=0}^{h-1}2^j(\kappa+h-j+6)
\ge\kappa(2^h-1).
```

For groups of the proposed size $`h\approx\tfrac14\log_2L`$ and
$`L=\Theta(n)`$, this is $`\Theta(n^{5/4})`$ code bits, and the
full $`2^{q_h}`$ catalogue exceeds the resource budget. Constant-height
groups can reduce the code size but restore an unbounded number of
sequential precision charges. This is a limitation of this assembly,
not a general lower bound.
One distinct candidate is a prefix-programmable high-precision terminal:
retain R79's actual high-precision word but replace its local data access
by charged global-prefix queries, without enumerating all precision
programs or all global prefix slots. Another is a coarse catalogue with
a separately paid fine correction. Neither interface is established here.
Do not let an unfinished stretch replace a completed primary proof.

## 6. Execution, review and stopping rules

Run the construction and assembly lanes in parallel with an independent
skeptical lane. The skeptic should reconstruct registers and addresses
from the native words, not accept the proposed resource summary.

| Elapsed time | Required work and decision |
|---|---|
| 0–25 minutes | Confirm the clean R95 baseline, immutable receipts, source contracts and all precision-dependent costs. Record actual start and deadline. |
| 25–90 minutes | Derive the new global stage and its full ledger. Independently check the finite height sequence and catalogue budget. Apply the qualification gate above. |
| 90–180 minutes | Complete the uniform assembly, width/precision splices and errors. If the primary route qualifies early, test a precisely scoped strengthening; otherwise finish the valid component or identify its exact unpaid term. |
| 180–235 minutes | Freeze the scientific claim. Obtain independent reconstruction and fresh adversarial review. Add only bounded fixtures addressing concrete surviving risks. |
| 235–300 minutes | Integrate supported results, run required repository checks and desktop/narrow rendering, publish the reviewed tree and merge only after CI passes. |

The preferred outcome is the complete precision bridge. A useful narrower
outcome is the globally priced source/correction stage at an explicitly
proved range. If a candidate fails, preserve its exact counteridentity
or unpaid resource, with a reopening question. Failure of that candidate
does not retire the general depth objective.

Do not claim dimension-growing lower bounds from query-architecture
failures, light cones that ignore free Clifford interlayers, or another
constant-layer obstruction. Those need separate witnesses and quantifiers.
Optimal general depth remains open even if the primary target succeeds.
Preserve R95 and the selected manuscript Result D (R79) until a new
proved statement justifies a precisely scoped update. No efficient
preprocessing, practical crossover or scalable emitter claim follows.
