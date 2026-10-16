# Borrowed precision and direct programming: execution outcome

[Executed design](BORROWED_PRECISION_BRIDGE_PROMPT.md) ·
[All-precision theorem R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer)

The session establishes **R96**, the general upper bound

```math
D_T=O(NL/b^2+L+\ell),\qquad
N=2^n,\quad \ell=\log_2(n+2),\quad n\ge1,\quad L\ge6.
```

The designed target was the same law only through $`L\le n^{1/4}`$.
The first completed bridge reaches $`L\le n^{3/4}`$. A separately
reviewed direct-prefix terminal then removes the precision restriction.
Optimal general depth remains open; this is an upper-bound theorem.

Execution began on **16 October 2026 at 11:02:53 UTC**, from main commit
`fd815938d7c4c73fb0460be22615a133f83103a3`, tree
`5f7926a90493137c4e1225f09be601eefd91de91`, after design PR #161.
The qualification checkpoint was 12:32:53 UTC; the authorized integration
ceiling was 16:02:53 UTC. The work branch is
`research/borrowed-precision-bridge-session`. Final manuscript drafting
remains separate.

## 1. Complete result and its scope

One native Clifford+T circuit simultaneously has

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\ell).
```

It uses one external initialized qubit and any supplied dirty allocation
$`b\ge17(L+n+7)`$ for complete real frames, or
$`b\ge17(L+n+8)`$ for the prescribed complex magnitude completion.
All constants are absolute, uniformly over every precision. Error is at
most $`2^{-L}`$ in full initialized-isometry norm: every logical column,
literal phase, signed and zero coefficient, work leakage, arbitrary dirty
input and external reference is included.

At sufficient width $`b=\Theta(\sqrt{NL}+L+n)`$, the additive depth is
$`O(L+\ell)`$. In particular the designed slice
$`L=\Theta(n^{1/4})`$ improves R79's $`O(n^{1/4}\log n)`$ to
$`O(n^{1/4})`$, and the same law now holds throughout the previously
uncovered intermediate precision region. The inherited worst-case count
lower bounds still match. No matching depth lower bound, efficient
preprocessing, practical crossover or scalable native emitter is asserted.
R79 remains the selected manuscript Result D.

## 2. What changed the construction

| Proof home | Role |
|---|---|
| [Global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md) | Separate the initialized source index from its dirty precision core; charge every long-address table to the actual global query envelope |
| [Fast capped reference](FAST_CAPPED_GEOMETRIC_REFERENCE.md) | Approximate the original geometric reference by a shorter exact cap, with a full-space comparison and all original flag coordinates retained |
| [Polynomial-precision bridge](POLYNOMIAL_PRECISION_DEPTH.md) | Independently complete three reductions and the existing full-payload catalogue through $`L\le n^{3/4}`$ |
| [Direct-prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md) | Retain a common fine-mask grammar and program only the shorter native coarse word; allocate no catalogue of full precision programs |
| [All-precision assembly](ALL_PRECISION_LINEAR_DEPTH.md) | Combine the two new mechanisms with original-width low-, narrow-, finite- and high-precision cases and literal complex phases |

The global stage has reserve

```math
\chi=O(\log(l+2)+\log^3(m+2)),\qquad r=\lceil A\chi\rceil.
```

Its ordinary source keeps only logarithmic initialized scratch. Its
precision-sized core is arbitrary dirty work, excluded from every source
anchor. Coarse geometry and per-prefix hash/matching selectors depend on
height m, not on final accuracy l. Every long local-data lookup is priced
with address at most $`n-r+C\chi`$, even when it happens not to depend
on the external prefix. Only short metadata tables and fixed native
words remain in the standalone $`2^{C\chi}`$ allowance. This is what
permits the first stage at m=n.

For original and shortened power-of-two caps $`R\lt M`$, the exact
identity is

```math
\|g_M-g_R\|^2=(4-2\sqrt2)2^{-R}.
```

A cap of length $`O(l+\log(p+1))`$, with p the number of bins, suffices
at final stage accuracy. The mathematical matching still uses the
original cap. A two-plane unitary comparison transports its entire
zero-anchor projector; the original high coin bits and all private
scratch stay in the nonzero tests. The actual unconditional reference
inverse and forward cancel on every occupied inactive input.

These two changes qualify the polynomial bridge. The remaining obstacle
was its catalogue of all final-accuracy programs: that catalogue grows
too large at linear precision. The direct terminal reads the required
program from the actual prefix instead. Its fine masks have a common
height-dependent grammar, and only the precision-independent coarse
native remainder is interpreted gate by gate. At group height
$`h\le\log_2(L+2)/32`$, even a deliberately generous
$`2^{3h}\mathrm{poly}(h)`$ instruction schedule fits the linear-L
depth and optimal count budgets.

The interpreter uses three original suffix bits on the active sector:
a local clean flag, a returned instruction selector, and a presence bit
for optional finite bins. All gates are explicitly guarded, and the whole
word is identity on inactive occupied work. Its precision base plus
separately reserved coarse work fits below b/4; an exclusive b/2 query
pool remains available at the original literal threshold.

For $`n^{1/4}\le L\le2^{cn}`$, with a sufficiently small fixed c chosen after all
native constants, the fast stages leave $`O(\log L)`$ layers and their
standalone/pool costs fit the wide branch. For larger L the original
height already is $`O(\log L)`$, so those stages are omitted. The
retained R79 high-precision clause handles $`L\gt N/n^2`$ directly.
A constant number of directly programmed groups suffices in the remaining
cases. This finite dependency does not invoke the new theorem on its own
components.

## 3. Independent reconstruction and corrections

Construction, assembly and skeptical lanes independently reconstructed the
new global stage and polynomial bridge. A fresh review, without the prior
verdicts, checked the same complete proof and found no blocker. The direct
terminal and all-precision assembly then received their own independent
native-word, width, common-grammar and full-error reviews. A second fresh
adversarial review traced the inherited interfaces without relying on the
earlier conclusions and found no remaining blocker.

The work made these points explicit before promotion:

- Keep the fine donor reserve independent of the later coefficient A.
  Test only used conditional work; the full saved suffix predicate is a
  separate paid operation and never enlarges every global table address.
- Classify every long-data lookup under the suppressed global query bill.
  Calling a table local does not justify an exponential standalone cost.
- Keep the original geometric matching when shortening its physical
  reference preparation, and transport all augmented anchor coordinates.
- Fix fine precision widths, physical arenas, type slots and cheap-row
  grammar uniformly across prefixes before computing the actual mixed
  gauge and its column replay. Preserve the compressed node/power address;
  padding to a whole-mode address would lose the query sums.
- Bound elementary Clifford gates as well as T gates before adding native
  controls. The universal coarse schedule excludes the long precision
  core; its genuine shared signal and metadata locations remain fixed.
- Use the proved dirty-helper controlled-phase words, including their
  residual phase corrections. No extra clean AND target is assumed.
- Give optional absent bins literal identity on occupied inputs by keeping
  their presence bit in both the query and central native-gate guard.
- Prove actual simultaneous allocation inequalities, including the
  precision core, separate coarse helpers and exclusive global query pool.
  A generic $`O(L)`$ width assertion would not establish a literal minimum.
- Retain the bottom-layer indexing $`j=1,2,3`$ in the inherited
  $`q_j=L+j+7`$ schedule; the explicit index set must agree with its
  quoted error budget.

On the new real branch the complete error is at most
$`249\eta/1024\lt\eta/4`$, where $`\eta=2^{-L}`$. Existing branches
retain their own full error guarantees. The complex theorem calls the
completed real result and literal diagonal at L+1, at exactly
$`17(L+n+8)`$ dirty width. All leakage propagates through actual subsequent
unitaries. Internal review does not replace external technical review.

## 4. Bounded evidence and integration

The [four new fixtures](../../tests/test_fast_capped_geometric_reference.py)
check the exact cap overlap and distance, the full-space two-plane
completion, every original two-bin band and input label, and the actual
inverse on occupied work with coherent dirty/reference entanglement.
Negative checks detect deletion of original high coin bits or private
scratch and substitution of a stale comparison inverse. Dense source
completions in these tests are verification objects, not native oracles.

The retained [native program fixtures](../../tests/test_grouped_terminal_depth.py)
check literal controlled alphabet gates, dirty phase corrections, actual
unquery, immutable addresses and occupied inactive programs. Retained
source, query, full-payload, sparse-correction and one-clean tests keep
their individual scopes. The new common grammar, global query sums,
uniform cutoffs and all-precision composition are analytic proofs.
No broad simulation or regenerated scientific receipt substitutes for
those arguments.

Local validation passed: all **1009 tests** in 366.983 seconds, including
the four new fixtures; the focused theorem/presentation set passed all
58 tests. Compilation, the reviewer walkthrough, real and complex native
examples, the q=16 residual example, both n=12 resource ledgers, and the
offline upstream check passed. All four fault-tolerant scientific receipts
matched their retained fixtures without regeneration.

Rendered validation covers all 16 changed Markdown pages at desktop and
narrow widths, with SVG and native math checks; it reports no failures or
TeX errors. The final whitespace check also passes. The associated pull
request records the exact reviewed tree, remote CI, and merge receipt.

## 5. Where the journey continues

The precision restriction in the upper law is closed. The remaining
question is whether this general depth law can be improved, or matched
by a lower bound under the same one-clean and arbitrary-dirty model.
At fixed accuracy, the known gap between constant lower bounds and
logarithmic upper depth remains. A new session should target that gap or
another explicit limitation of the complete native circuit model. It
should not reopen the now-priced global precision reserve or catalogue
restriction as though they were still missing from this construction.
