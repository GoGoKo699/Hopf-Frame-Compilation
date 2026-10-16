# Product-coin continuation: causal batching and recursive depth

[Current frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth) ·
[Product-coin theorem R87](COIN_DILATION_DEPTH.md) ·
[Coin terminal R88](PREFIX_SELECTED_COIN_TERMINAL.md) ·
[Causal compact coins R89](CAUSAL_COMPACT_COIN_DEPTH.md) ·
[Selected routing](SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)

This is the continuation requested after the R87 integration, on
16 October 2026. The repository baseline is main commit
`83f2220cd0e3fafa53d47766bf248db35a128c4c`, merged in
[PR #151](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/151).
The original five-hour ceiling remains 06:19:57 UTC, measured from
01:19:57 UTC; the continuation does not restart that allowance.

The goal is still optimal general depth for the prescribed complete
operator. Every candidate retains all columns, literal phases, actual
inverses, one external initialized qubit, the original dirty widths,
and full initialized-isometry error including work and references.
The continuation improves the fixed-accuracy general upper bound to

```math
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(n+2))\right),
\qquad N=2^n,\quad\ell=\log_2(n+2).
```

Here $`\log_2^*x`$ counts repeated binary logarithms until the
argument is at most two. The same circuit retains
$`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, one external clean qubit,
and the original real and complex dirty thresholds. R89 replaces the
serial coarse search and applies the new mechanism recursively.
The absolute all-precision theorem R79 is unchanged. No final manuscript,
new external resource assumption, or scalable native emitter is part
of this continuation.

## 1. Exact routing without the extra selector-depth factor

The R87 layout exchanges two contiguous logical intervals with fixed
external banks. It does not require a fully general selected permutation
network. The new [routing specialization](SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
loads the two endpoints and lengths into four short protected words
once, uses guarded barrel rotations and variable-prefix swap echoes,
and returns those parameters by the actual inverse lookup.

Each fixed rotation bit is two matchings. Its guard reads only H and
that parameter bit, so its dirty query and literal Fredkin echo have
constant T-depth. The logarithmic number of rotation bits gives total
layout-and-return depth $`O(\ell)`$. The prefix-swap predicates
are paid parallel-output tables. Arbitrary padding and former bank
contents return exactly, including correlations with external systems.
Each completed macro is identity on H zero even with occupied
parameters. The additional conditional clean metadata is
$`O(\ell)`$ and is reserved before choosing the independent
protected-bank margin.

Two independent reviews checked the full permutation and resource
interfaces. Three bounded fixtures exercise all valid small intervals,
dirty-control offsets, occupied padding, actual inverses, inactive
parameter values, and a nontrivial body between layout and unlayout.
These fixtures check the logical interfaces; the native proof supplies
the phase and asymptotic bounds.

## 2. A smaller prefix-selected terminal bound

[R88](PREFIX_SELECTED_COIN_TERMINAL.md) proves, for every
$`0\le M\le n`$, the terminal estimate

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell\log_2\log_2(M+4)\right).
```

It uses the original real width $`b\ge17(L+n+7)`$, one
external initialized qubit, and terminal error strictly below
$`2^{-L}/4`$. In the complete complex construction the supplied
leaf diagonal remains separately priced, as in R87.
For the protected height $`M=O_L(\ell)`$, the terminal overhead
is now $`O_L(\ell\log\log(\ell+4))`$.

The proof rederives the local product-coin stage with a preserved
global prefix. Source precision, masks, donor space, and polynomial
factors depend on the remaining local height m. The long prefix
contributes table rows and logarithmic query depth, but no polynomial
in global n to the count multiplier. Parallel one-bit sharp queries
avoid a serial membership-output term; a read-only replacement of
the high-address traversal makes their physical T-layer parallelism
valid. The explicit degree-twenty envelope is fixed before choosing
the suppression margin.

At each stage $`m\mapsto O_L(\log m)`$, only that stage's
early factor is rounded. The exact early/terminal factorization
justifies different denominators at later stages. Errors and
suppressed counts sum, and the retained smaller-error bottom preserves
the literal width. Signed inputs use one finite sign reduction and
two sign diagonals outside the recursion. Two independent full-proof
reviews found no blocker. Three bounded matrix fixtures check every
column of the selected factorization, occupied protected sectors,
different stage denominators, and the maximum-over-prefix error rule.

This result does not imply an $`O_L(\ell)`$ terminal. Even a
constant-call coarse decoder at each recursive stage would leave an
unbounded number of serial global-prefix source episodes. A new
schedule or paid return-preserving amortization is needed to remove
that cost.

## 3. The coarse-control obligation and its causal repair

While the transformed coin state is present, its last nonzero block
can determine the old coarse label. Erasing that label at this point
does not complete the circuit: the selected source inverse, layout
inverse, and descendant reload still need it. Copying the label to a
parameter register moves the same return obligation.

Membership can reconstruct the local leaf label without changing the
coin state. The subsequent conditional source inverse changes the
coin state. Replacing its preserved selector by a function of those
changing coins is therefore not automatically its actual inverse.
A valid replacement must prove the whole consumed-control isometry
and its inactive extension, rather than assume gatewise controls stay
valid. This is a gap in that proposed word, not an impossibility theorem
for other coarse decoders.

There is also an exact common-basis formulation. If C simultaneously
conjugates the simple computational suffix projectors to all required
coarse scaling projectors, the entire adaptive search can be replaced
by C, one classical coarse-label XOR, and its actual inverse. The
ordinary complete prefix-controlled state-preparation unitary supplies
such a C algebraically, but its serial rotation word is not the desired
native compiler. Independent one-column source completions do not
automatically share a common C. Bases inside each corresponding
coarse annulus may be chosen freely.

The successful change is the [causal compact-coin construction](CAUSAL_COMPACT_COIN_DEPTH.md).
Each coarse block has a fixed weighted coin with only a logarithmic
number of added qubits. Fiber-cardinality groups and staged shear
codes compress each original child label enough to expose its rank
workspace. All original prefixes are still available when the short
compression metadata is loaded. The joint retained code then
determines every original leaf, so one paid table clears all discarded
high fields and one source prepares all conditional ranks.

A local matching assigns each custom atom to an equal-weight fixed
coin atom. Shared short XOR-shift catalogues allow all these
prefix-dependent matchings to be implemented in three query episodes:
load their indices from the original custom tuple, apply the disjoint
shifts, and erase the indices from the final tuple. Every old metadata
field becomes fixed zero padding. The inverse table decodes the
original prefixes in finite classical preprocessing; it is not an
uncharged coherent causal decoder.

This locality is essential. An earlier proposed global matching
preserved atom weights but could exchange the two blocks of a balanced
example, moving a marker into the wrong coarse band. The exact
[counterexample](COMPACT_COIN_FLAGS.md#6-native-scope-and-bounded-evidence)
is retained. R89 chooses its matchings locally, producing new
conditional partitions while preserving every product-prefix rectangle.
Parallel fixed-coin inverse sources
then expose the coarse band through the last nonzero augmented block.
The test includes each source's work as well as its coin.
The existing nearby-completion lemma supplies the full operator error
comparison; the circuit uses only actual physical inverses.

The new [short-label routing lemmas](SHORT_LABEL_SUBSET_ROUTING.md)
prove the dense matching and sparse compression primitives with their
exact query words. They choose their own matching and code, and do not
claim an arbitrary prescribed permutation compiler.

## 4. Prefix-local recursion and the remaining logarithmic target

The first complete application, with protected height
$`M=O_L(\ell^3)`$ and the independent R88 terminal, already gives
$`O_L(N/b^2+\ell\log\log(\ell+4))`$.
The stronger [prefix-selected extension](CAUSAL_COMPACT_COIN_DEPTH.md#10-prefix-selected-recursion-and-iterated-logarithm-depth)
chooses catalogues separately for each immutable external prefix.
Only the changing local prefix needs a common catalogue. Thus all
initialized metadata and source work depend on the remaining height m,
with protected remainder $`r=O_L(\log^3(m+2))`$.

The proof rederives the component bounds uniformly in the stage
precision; it does not invoke a fixed-accuracy theorem at a silently
varying accuracy. Parallel sharp queries avoid a polynomial count
factor in the long external prefix. An independent protected margin
suppresses every addressed metadata and source-work pattern. Counts
and errors sum, and every two recursive steps decrease the remaining
height at least as far as one binary logarithm. This proves the
iterated-logarithm depth and the original literal widths.

Independent reviews checked the causal all-column identity, the
selected rank-source extension, global read/target separation, actual
inverse and occupied-inactive behavior, simultaneous work, and recursive
precision/count sums. Four constructor fixtures compare with an
independently built complete frame, including real hash collisions,
different prefix assignments, all 64 frame columns, arbitrary old tags,
and a coherent external reference. They also exhibit why source work
must be included in the nonzero-block test.

The remaining factor is an unbounded number of serial episodes that
read the long unchanged prefix. The present query interface charges
$`O(\ell)`$ to each. Removing that factor requires a paid way to
combine or amortize those episodes, or a construction that avoids
the recursion. It cannot be removed by declaring $`\log^*n`$
constant. General $`O_L(\ell)`$ depth, optimal depth, and a
dimension-growing unrestricted lower bound remain open.

The bounded follow-up tested extending the same selector batching
across terminal stages. It does not transfer directly: a terminal
factor changes the prefix used by the next factor into a quantum
superposition. The next source selector is therefore not a short
classical function of the original basis input. A joint branch source
would need a new way to consume its old-input and quantum-history
coordinates while preserving every frame column. This isolates the
missing interface in that proposal; it is not a native lower bound.

## 5. Validation and integration

The linked proof chapters are the authoritative claim homes. This
record separates their completed arguments from the failed global
matching and the remaining logarithmic target.

The complete deterministic suite passed **978 tests**, including 19 new
bounded fixtures: three interval-router tests, three prefix-terminal
tests, seven compact-coin tests, two subset-routing tests, and four causal
constructor tests. Independent proof reviews cover the asymptotic and
full-isometry claims that these finite fixtures cannot establish.

Python compilation, the ten-check reviewer walkthrough, both retained
native coarse examples, the residual-QBP example at q=16, the exact
fault-tolerant fixtures, and offline upstream synchronization passed.
Both resource-ledger commands completed at n=12. Their retained audit
status fields are unchanged; executing a ledger is not a new theorem
certificate. The final documentation gate covers all 160 Markdown pages,
the repository diagrams, MathJax SVG and native MathML, and desktop and
narrow layouts. A display-width failure in R88 was corrected by
line-breaking the unchanged equation.

Integration is gated on the reviewed tree and the required Python 3.11,
Python 3.13, and rendered-presentation CI jobs. The pull request records
the exact tree and gate receipts. No final manuscript is drafted, and
the bounded fixtures are not a scalable native emitter. The continuing
objective remains optimal general depth.
