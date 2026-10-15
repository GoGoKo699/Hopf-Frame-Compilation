# Joint decoder and terminal depth session

This session started from R79 on 15 October 2026, with an authorized
five-hour ceiling. Its target was a complete fixed-accuracy circuit below
the general square-logarithmic depth bound, retaining one clean qubit,
square-root T-count, linear Clifford count, arbitrary dirty/reference
inputs, and literal work return. It did not change the general frontier.

Three constructive results survived independent internal review:

- [R80, homogeneous-magnitude depth](HOMOGENEOUS_MAGNITUDE_DEPTH.md),
  gives logarithmic depth at fixed accuracy and square-root width for
  the complete family whose angle magnitudes
  depend only on tree depth. Node signs and supplied leaf phases may vary.
  The signed balanced real subfamily has an exact native version.
- [R81, recursive terminal compilation](RECURSIVE_TERMINAL_DEPTH.md),
  reduces the last logarithmically many layers of an arbitrary frame to
  a log-times-log-log depth contribution. Its early decoder remains
  square-logarithmic, so this supporting lemma does not improve R79's
  complete-frame bound.
- [R82, bounded-density compilation](BOUNDED_DENSITY_DEPTH.md), gives
  a sub-square-logarithmic complete circuit when the root probabilities
  stay within fixed positive factors of a supplied classical product
  distribution. Its constant-normalization early construction and
  recursive terminal have depth
  $`O(\log(n+2)(1+\log^*(n+2)))`$ at square-root width and fixed
  accuracy and density bounds.

The first result uses a common subspace conjugator that is cheaper than
the full prescribed frame. The second revisits the terminal as a smaller
frame with an unchanged quantum prefix. The third bounds a weighted-Haar
correction and supplies its row and column preparations natively. None supplies
a free evaluator, an initialized dirty register, or a general emitter.
The selected principal theorem remains
[R79](UNIFORM_SUBTREE_DEPTH.md); final manuscript drafting was not part
of this session.

## 1. The complete circuit that survives

For depth-homogeneous magnitudes, a product of one-qubit unitaries B
transports every computational suffix-zero subspace to the corresponding
subtree scaling subspace. A tagged column source can therefore erase
its tag by conjugating a classical last-one-position XOR with B.
There is no binary sequence of transported threshold reflections.

The native construction uses the same actual synthesized factors in the
source and B. Each factor is a product of two exact involutions, so its
selected implementation has exact dirty-mask return. In this real
common-flag core, tag erasure, conditional-suffix return, and the sole
external activity flag then return exactly. Only the logical frame
differs from the target rotations.
This avoids treating an approximately synthesized real rotation as an
exact real gate.

Recursive protected suffixes are logarithmic in the current local
dimension. After the existing sign gauge is removed, the remaining
terminal is independent of the earlier logical prefix. Its successive
local logarithmic depth charges have a geometric sum. R80 gives the
native word, physical T-layer schedule, certified rounding, error sum,
and all-width fallback. The bounded
[operator fixtures](../../tests/test_homogeneous_depth.py) test its
literal algebra and occupied-work behavior, not its asymptotic ledger.

## 2. Bounded density removes the repeated early comparisons

The retained precision-gauge identity expresses a positive frame as a
diagonal dressing of a reference frame minus a subtree-scaling correction.
For density bounded above and below relative to a positive product
reference, a finite martingale variance identity bounds the correction's
Carleson sum. The resulting operator norm bound also holds for its
entrywise absolute matrix. This is the feature missing from a generic
entry-indexed encoding of the full frame.

Positive Schur weights then give normalized row and column source states
with constant overlap normalization. Their common coordinates are one
leaf register and a short depth tag. The column source uses a preserved
subtree prefix; the row source preserves the leaf and prepares only the
short depth state. Orthogonal failure flags complete both isometries.
The native descendant-source interface prices these preparations and
their actual inverses.

A constant-size linear combination joins the correction to the reference
frame. Distinct product-factor banks retain rejected components, and
tuned constant-round amplification returns the complete prescribed
operator, including literal phase. A full-unitary comparison controls
leakage; a block-error estimate alone is not substituted for work return.
The outer reflection tests all protected work and retains the original
saved activity guard on occupied inactive sectors.

Each early stage now uses a constant number of global-prefix source
queries. Recursive local heights decrease to their logarithms, leaving
an iterated-logarithmic number of stages. Their counts have geometric
budgets, while each stage still pays the global logarithmic query depth.
Conditioned density bounds are always taken from the original global
ratio; they do not deteriorate at successive stages.

The reference is a classically specified positive product distribution,
not a supplied quantum state or a freely available frame oracle. R80's
native common-flag word prices its complete reference frame. Uniform
reference gives the bounded-uniform-density special case. Reference
biases may be arbitrarily strong, but the density ratio bounds are fixed.
An arbitrary target cannot be made eligible by choosing that target as
its own unpriced reference.

The [weighted operator fixtures](../../tests/test_weighted_haar_depth.py)
check the gauge, finite variance identities, sparse source overlap,
literal amplified return, and retained leakage. The dimension-independent
norm bound, native ledger and all-width recursion remain analytic.

## 3. A shorter general eraser saves calls, not an order

The original R78 eraser reconstructs a separate label, XORs the old tag,
and uncomputes that label. There is also a direct exact eraser using the
old tag T. This observation is optional; R79 retains its original proof
and constants.

Let the data lie in band t, with T initially equal to t. Process bits
j from low to high. The unprocessed higher bits still equal those of t.
Compute a threshold in a separate zero selector D:

```math
d_j(T)=\sum_{k>j}2^kT_k+2^j-1.
```

Apply a Hadamard to bit j, the R78 controlled reflection about
the nested subtree space at D, another Hadamard, and a CNOT from the
saved activity flag H to bit j. Uncompute D from the unchanged higher
bits. On the active band, the reflection toggles bit j by
the predicate that t is at most the threshold; the final CNOT
complements this value. Thus the old bit is cleared. On the inactive
sector the completed word is exactly identity even on occupied work.

With s tag bits this uses s transported reflections, hence two s
descendant-source appearances and s anchor phases. The separate search
register and guarded tag XOR disappear. The live conditional bank is
two s selector/tag bits plus the source bank, rather than three s plus
the source bank. The R78 full-mode comparison gives erasure error at
most four s delta, using the same linear hybrid estimate.

This halves the reflected-source work but leaves logarithmically many
serial lookups, each with logarithmic address T-depth. It is not a
sub-square-logarithmic general compiler. Coherent errors need not add
as the square root of the number of bits.

## 4. Threshold products and arbitrary source completions

For an oblivious product of threshold reflections, let m_d be the
multiplicity of the reflection at d. Its sign on band t obeys

```math
\sigma(t)=(-1)^{\sum_d m_d[t\le d]},\qquad
\frac{\sigma(t+1)}{\sigma(t)}=(-1)^{m_t}.
```

Every boundary at which the desired band-bit sign changes must
therefore occur with odd multiplicity. The lowest bit changes at every
boundary. Replacing adaptive binary search by that direct product
requires linearly many thresholds. This is a restriction on the stated
product architecture, not a lower bound for arbitrary coherent algorithms
or native circuits.

A single descendant source compresses its own scaling subspace to
computational anchors. Its arbitrary unitary completion is not a
prescribed frame on the complementary marker bands. The higher branch
also retains a coherent prefix in subsequent lookup addresses. Dropping
either feature does not establish a smaller worst-case query recurrence.

The ordinary top-down state-preparation circuit does transport the entire
nested flag, without fixing the orientation inside every band. This is
more freedom than synthesizing the prescribed frame itself. For arbitrary
node-dependent magnitudes, however, its directly priced schedule still
has one stage per tree depth. R80 exploits the case where that same flag
has an independently cheap product conjugator.

## 5. A generic dirty-cache correction already contains a lookup

Let E add the one-hot vector at a read-only address a to an arbitrary
dirty bank Y. For an arbitrary Boolean table f and a fixed bank position
j, consider two exact Clifford payloads on one logical target:

```math
B_0(Y)=Z^{f^{\mathsf T}Y},\qquad B_1(Y)=X^{Y_j}.
```

Differencing each payload separately gives, in chronological order,
the desired Z and then X actions. Combining them first and applying one
outer dirty difference instead gives

```math
\begin{aligned}
F(Y)&=X^{Y_j}Z^{f^{\mathsf T}Y},\\
F(Y+e_a)F(Y)^\dagger
&=(-1)^{f(a)Y_j}X^{[a=j]}Z^{f(a)}.
\end{aligned}
```

The extra phase depends on arbitrary dirty input and is forbidden by the
full operator contract. Returning the cache does not remove it. If a
general correction C_f implemented that phase, then

```math
C_fX_{Y_j}C_fX_{Y_j}=(-1)^{f(a)}I.
```

Thus a depth-d generic correction already supplies an arbitrary Boolean
phase query in depth at most two d. This is a reduction for this proposed
cache interface, not an unrestricted query-depth lower bound.

There is a separate elementary limit on a more restrictive cache-only
architecture. Suppose caching acts only on the prefix and arbitrary
dirty work; the interior processor never touches the prefix; and
uncaching again touches only prefix and dirty work. Both caching and
uncaching must preserve each computational prefix label: the prefix is
read-only, controlling a unitary on the dirty work. On a basis prefix,
choose the dirty work maximally mixed. Caching leaves it unchanged and
independent of that prefix. The interior processor therefore induces a
prefix-independent reduced channel on the payload, and uncaching cannot
alter that reduced channel. Such an architecture cannot implement
distinct prefix-dependent payload unitaries. Continued Clifford access
to the prefix, a cache that changes the prefix label, initialized cache
bits, or a cache that also touches the payload lie outside this argument.
In particular, swapping the prefix into the dirty cache and back is
excluded by the read-only hypothesis.

## 6. Row-wise fingerprints do not approximate an arbitrary table

Replacing every one-hot row equality by a small linear fingerprint
does not give a global approximate lookup. Let a full-rank linear map
send s bits to m bits, with m smaller than s, and put

```math
v_a(y)=[H(y)=H(a)].
```

For the bilinear table matrix equal to the identity, and equal input
addresses, the fingerprinted query returns

```math
v_a^{\mathsf T}v_a=2^{s-m}\bmod2=0,
```

where the exact one-hot query returns one. This failure is deterministic,
despite a possible small collision probability for each separately fixed
false equality. A useful global approximate indicator needs its own
full-operator error proof. An improved isolated zero test does not remove
the arbitrary phase tables in the descendant source.

## 7. What would reopen the general target

The general early decoder still has logarithmically many sequential
source calls at logarithmic address T-depth. A qualifying next mechanism
must change that complete product, for example by supplying:

- A native common flag conjugator for arbitrary node-dependent magnitudes,
  with its arbitrary completion, work return, and all resource costs paid.
- A global dirty lookup with smaller depth and a full-input error bound,
  preserving the square-root count and linear Clifford budgets.
- A different complete tag-erasure circuit that uses additional frame
  structure beyond the current threshold interface.

R82 combines a different early mechanism with recursion on a promised
family, but its constants depend on fixed density bounds. Extending that
normalization to arbitrary node-dependent magnitudes, including singular
branches, would require a new argument. The recursive terminal and the
constant-factor direct eraser do not close the general gap by themselves.
The successful comparison cases show that arbitrary signs and substantial
magnitude variation can coexist with improved depth; they do not bound
the worst case over all prescribed frames.

## 8. Review and verification

Separate internal reviews checked the common-flag circuit, the
prefix-selected terminal ledger, the weighted source geometry and
amplification, and the product-reference recursion. A final scope review
kept R79 as the unrestricted frontier and made the read-only-prefix
hypothesis in the cache-only argument explicit. These are internal
technical audits, not external peer review or a priority assessment.

Ten new bounded fixtures accompany the proofs. Four cover the R80
frame/tag identities and native involution echoes; six cover the R82
weighted gauge, variance, sparse sources, amplification and leakage.
Independent test review confirmed their stated boundaries: R82's
inactive-sector fixture checks disabled source macros and actual PREP
cancellation, while the internal native enable wiring is proved in the
chapter. None of these checks substitutes for the asymptotic resource
argument or emits a general compiler.

The full validation suite passed all 932 tests, including the ten new
fixtures. Compilation, the reviewer walkthrough, real and complex native
examples, the residual-QBP example, offline upstream reconciliation, and
both resource ledgers passed. All four immutable fault-tolerant receipt
suites matched their existing scientific fixtures; no receipt was
regenerated. Rendered-document and remote-CI outcomes are recorded with
the reviewed pull request.
