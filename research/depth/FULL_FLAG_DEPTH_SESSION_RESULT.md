# Complete flags and local precision: execution outcome

[Executed design](FULL_FLAG_DEPTH_PROMPT.md) ·
[R95 theorem](LOGARITHMIC_PRECISION_DEPTH.md) ·
[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer)

The session establishes **R95**, a uniform logarithmic-depth upper bound
at inverse-polynomial accuracy in the logical qubit count. It improves
both costs that remained in R94: the opening spectral search and the
local precision schedule. Optimal general depth remains open.

Execution began on **16 October 2026, 09:39:11 UTC** from main commit
`5833a19fc0faf58b74c7795664636d2cb3e28a77`, tree
`ea8eb50cf4cb10c8a0c2ea0949912a6d00a13cca`, after design PR #159.
The authorized ceiling was **14:39:11 UTC**, including review, validation,
CI and integration. The work branch is `research/full-flag-depth-session`.
Final manuscript drafting remains separate.

## 1. Established result

For $`N=2^n`$, $`\ell=\log_2(n+2)`$ and every integer
$`6\le L\le2\ell`$, R95 gives one circuit with

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\ell).
```

There is one external initialized qubit. Literal dirty thresholds remain
$`17(L+n+7)`$ for complete real frames and $`17(L+n+8)`$ for
complex magnitude completions. Error is at most $`2^{-L}`$ in full
initialized-isometry norm, with arbitrary logical inputs, literal phases,
zeros, all work leakage, arbitrary dirty states and external references.
Constants are absolute on the stated range. The inherited count lower
bounds still match; no depth-optimality claim is made.

At $`L=\Theta(\ell)`$ and square-root dirty width, the new depth
is $`O(\ell)`$, compared with R94's $`O(\ell\log_2(\ell+2))`$.
R79 remains the selected manuscript Result D and the absolute
all-precision statement. No efficient preprocessing, practical crossover,
or scalable native emitter is asserted.

## 2. What changed the problem

The useful freedom was to spend far more resources inside a small
catalogue entry, while retaining the global optimal count. This does not
create extra clean qubits. The construction removes the local demand
for them and transports every remaining physical workspace wire.

| Proof home | New role |
|---|---|
| [Coarse flags and sparse correction](COARSE_FLAG_SPARSE_CORRECTION.md) | Round the mathematical tree only to inverse-polynomial accuracy in its local height, implement that frame at final accuracy, then correct its signed sparse difference with constant-normalization recovery |
| [Small-workspace descendant source](SMALL_WORKSPACE_DESCENDANT_SOURCE.md) | Use the existing borrowed real-rotation word for calibration; retain only a logarithmic initialized expansion index and one attenuation bit |
| [Local depth with relaxed counts](RELAXED_COUNT_LOCAL_DEPTH.md) | Use exact compact geometric references, two explicit height reductions, and actual-prefix catalogue slots to reach a constant number of R79 high-precision groups |
| [Complete R95 assembly](LOGARITHMIC_PRECISION_DEPTH.md) | Insert those local words into R94's complete-payload catalogue and retain the global query sums, original-width fallbacks, error budget and complex phase composition |

For the opening, each row of $`D=W-W^0`$ has at most m+1 supported
ancestor/root entries. Coarse per-node angular error
$`1/[64(m+1)^2]`$ therefore gives $`\||D|\|\lt 1/16`$.
Signed Schur sources encode D, including its nonzero state column.
The normalization-two block $`(W^0+D)/2=W/2`$ recovers the full frame.
The physical coarse word is itself compiled to final accuracy; its
leakage is not mistaken for a sparse matrix correction.

The opening's simultaneous references share one paid gradient. Exact
integer sums of individually rounded whole-bank phases preserve tensor
factorization on the gradient invariant sector. This gives conditional
reserve $`O(l+\log^3(m+l+2))`$ and constant sequential global-prefix
query positions. After the retained first R89 stage it leaves $`O(\ell)`$
original logical bits.

The local construction instead uses exact compact geometric sources.
Padding each dyadic exponent to one less than a power of two permits the
retained recursive Gray-path word: linear gate count and logarithmic clean
width. Borrowed real calibration removes the final-precision clean bank
from every other source. Two decreasing local reductions leave
$`O(\log k)`$ layers at precision $`2^{-k}`$; a constant number of
small groups then lies in R79's proved high-precision regime. Selecting
one slot per actual prefix avoids paying for a universal rounded program
inside this local construction.

This is a finite decreasing dependency. It uses R94's transfer identity,
not R94's restricted-precision theorem on its own entries. Every nested
slot carries its full compiler workspace; occupied slot actions cancel
by the same exact commutator identity.

## 3. Independent reconstruction and corrections

Separate opening, local-source, alternative-construction, and skeptical
lanes reconstructed the operator and resource arguments. The native-model
audit independently checked the exact reference source and alias/matching
padding. A further review checked the complete global assembly.
A final adversarial review, started without the earlier audit conclusions,
independently reconstructed all four chapters and found no blocking gap.
Its separate borrowed-source subaudit checked the full-operator estimate
and source completion with arbitrary dirty cores.

The reviews required these explicit corrections before promotion:

- Round each complete local coefficient-bank eigenphase once before
  summing integer programs, avoiding an omitted index-length error factor.
- Count P-dependent hash and matching tables in the global address envelope.
- Apply source completion on the entire physical space, including arbitrary
  dirty cores. Approximate dirty return is not exact tensor factorization.
- Rebuild causal matchings for the padded exact reference coins; preserve
  their class-density bounds and conditional scratch allocation.
- Handle every residual height at most three directly, including arbitrarily
  large precision; omit invalid negative-height catalogue blocks.
- Fix component constants and group counts before the error-padding constant,
  then enlarge the common absolute cutoff. The local precision inflation
  remains inside a fixed-polynomial exponential resource allowance.

The complete real error budget is $`499\eta/2048\lt \eta/4`$ on the
large wide branch. The complex case pays L+1 within the internal real
range and composes the retained literal diagonal. Narrow and bounded-size
cases use the original-precision R79 circuit at its original literal
width. Independent review found no remaining substantive gap after these
corrections. Internal review is not a substitute for external technical
review.

## 4. Scoped routes that did not supply the construction

Root-state preparation alone does not determine the complete spectral
flag. For example, $`U_1=H\otimes H`$ and
$`U_2=U_1\mathrm{SWAP}`$ have the same zero column, but
$`U_2^\dagger U_1=\mathrm{SWAP}`$ exchanges the last-nonzero labels
of $`|10\rangle`$ and $`|01\rangle`$. The successful construction
retains R89's causal band geometry.

There is also a restricted obstruction to basis-permutation/support
compression of independent fine geometric coins. For a product of a
coins with probabilities $`1/2,1/4,\ldots,2^{-\tau},2^{-\tau}`$,
$`\tau\ge2`$, let B count coordinates outside their first atom.
Then $`B\sim\mathrm{Binomial}(a,1/2)`$ and
$`\Pr(B\lt a/4)\le4/a`$. Every tuple with $`B\ge a/4`$ has weight
at most $`2^{-5a/4}`$. Hence any basis support of size at most
$`2^{a+r}`$ has probability mass at most

```math
4/a+2^{r-a/4}.
```

For $`a=\Theta(\ell^3)`$ and $`r=O(\ell)`$ this tends to zero,
so that particular support-pruning route cannot even approximate the
root. This is not a lower bound for general coherent compression or
native circuits; the coarse/fine construction bypasses it.

A borrowed phase clock also fails to remove a noncommuting program's
unknown phase sector automatically. Let A be H on the first data bit,
B be CNOT, $`W=BA`$, and
$`V=|+\rangle\langle+|\otimes I+|-\rangle\langle-|\otimes W`$.
For $`P=Z`$ on the clock,

```math
PVPV^\dagger
=|+\rangle\langle+|\otimes W
 +|-\rangle\langle-|\otimes W^\dagger.
```

It has occupied-clock error $`\sqrt2`$ from $`I\otimes W`$.
A cheap family of genuine powers of W would be a different premise.
Neither this counterexample nor the older dirty-cache failures establish
an unrestricted depth lower bound.

## 5. Bounded evidence and integration

The [four new operator fixtures](../../tests/test_coarse_flag_sparse_correction.py)
check the sparse absolute norm through height four; signed Schur overlap
through height three with a live root and zero residual; normalization-two
recovery with all work retained and literal phase; and actual-inverse
composition under source leakage. Coherent activity sectors and entangled
dirty/reference states are included. Negative fixtures detect deleted root
corrections, discarded signs, incomplete anchors and stale inverses.
Their dense source completions are verification matrices, not supplied
native oracles. Asymptotic resources and cutoffs are proved analytically.

Local validation passed:

- `python -m compileall -q compiler_robust_hopf scripts tests validate.py`.
- `python validate.py --quiet`: **1005 tests passed**, including the four
  new fixtures, in 366.905 seconds with one BLAS/OpenMP thread.
- All 54 focused operator/documentation guards passed: the new fixture and
  reviewer narrative, literature policy, inline-math, presentation and
  typography suites.
- The reviewer walkthrough, real and complex coarse native examples,
  residual QBP example at q=16, fault-tolerant receipt verifier, offline
  upstream check, and both n=12 resource ledgers passed. All four immutable
  fault-tolerant receipt suites matched without regeneration.
- All changed Markdown pages were checked at desktop and narrow widths,
  with both SVG mathematics and native MathML. No TeX errors, missing
  images, page overflow or oversized desktop equations were found.
  Wide mobile equations and tables retain local scrolling.
- `git diff --check` passed.

Remote CI runs, the published tree and the merge receipt are recorded in
the associated pull request. Integration requires all analytic and
presentation checks to pass and the remote head to match the reviewed
tree. No immutable scientific receipt is regenerated.

## 6. Where the journey continues

The session's uniform logarithmic-depth target is reached. The general
answer still requires either a better unrestricted upper bound or a
matching lower bound, and the intermediate precision region remains
broader than R95's range. A next design should examine which of those
questions can be reduced to a new, reviewable mechanism. It should not
repeat R94's resolved opening or local-precision obstacles, or treat
the absence of a matching lower bound as evidence of optimality.
