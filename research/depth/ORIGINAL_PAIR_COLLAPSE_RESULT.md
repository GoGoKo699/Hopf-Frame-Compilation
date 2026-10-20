# Revision: use Booleanity before resynthesizing negative singletons

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Pair-collapse proof](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md) ·
[Previous structural revision](STRUCTURAL_PAIR_REVISION_RESULT.md)

This continuation starts from main
`3c901fdffb6b48548fc67ac7056481db5ac4136e`, tree
`0de9c050745f9504bb78bbd4b34608a59d989f4b`, after PR #209.
That revision replaced the large all-active feature estimate by a
structural six/five sign/XOR degree bound. The new step removes the
remaining physical-pair hypotheses rather than adjusting that estimate.

## 1. Destination, theorems, components and representation limits

| Role | Current assessment |
| --- | --- |
| Destination | Optimal unrestricted native T-depth with simultaneous count, Clifford, dirty-width, full operator error and return budgets |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md), with selected [R79](UNIFORM_SUBTREE_DEPTH.md) and all literal allocations unchanged |
| New original-pair theorem | [Boolean pair collapse](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md): if an exact complete one-clean native-three Boolean word has any genuine original middle pair, it has exactly one; its difference sums every actual active axis at both endpoints, and no negative singleton remains after literal repair |
| Small degree consequence | The retained [all-active diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md) now gives sign/XOR degree at most six/five for every paired source, and five/four at common active rank at least four |
| Unrestricted three-layer consequence | Combining the paired theorem with the retained no-pair and Pauli-interior theorems gives an absolute, width-independent Boolean-degree ceiling for all at-most-one-clean native-three words. Its universal constant is still the crude no-pair constant, not six/five |
| Useful component | Literal cancellation is paid when a negative singleton equals an actual endpoint axis. Membership in an endpoint span is not a native cancellation identity |
| Changed mechanism | The scalar pair sector preserves the actual silent-list independence. Boolean reality removes silent couplings directly, so no product-axis resynthesis is needed |
| Unpaid upper operation | A complete sublogarithmic weak quantum-fiber source or whole extraction within the retained simultaneous budgets |
| Still open | A small quantitative universal three-layer degree bound, higher-layer degree growth, and dimension-growing unrestricted depth |

The new theorem applies to the actual unrestricted **one-clean,
three-native-layer, complete Boolean model**, not merely to one chosen
pair representation. It does not cover arbitrary clean count, general
complex diagonal targets, every fixed approximation tolerance, or
arbitrary depth. In particular, the retained many-pair complex-diagonal
examples are not Boolean counterexamples.

## 2. The first unsupported inference, and its replacement

The initial negative-singleton question had a real unpaid step. When
a middle singleton is literally an endpoint native axis, its two odd
T powers combine into an actual Clifford, with the signed-axis scalar
retained. The residual lists are literal subsets. But a singleton that
is only a product of endpoint axes is not a gate in that list. Putting
it into a new basis changes the phase polynomial; its cancellation
cannot be credited without a separate identity and resource charge.

The [new proof](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md) does not perform
that replacement. It uses only the original middle pairs and the actual
silent endpoint list:

- Fix all original pair parities to their scalar branch, after the
  retained single shared computational-center normalization.
- The silent list remains independent even after the physical clean
  and these parity restrictions. An actual active-axis witness excludes
  the clean Pauli; full endpoint-list independence excludes every
  combination of pair differences.
- On that full remaining-register Boolean block, a one-band Pauli
  transport has a real Boolean derivative. Native independence forces
  every obstructing silent coupling to vanish. This promotes individual
  Pauli commutation signs, not an arbitrary sector identity.
- Comparing the all-scalar sector with each single-quarter-turn sector
  forces every genuine pair difference to be the same complete active
  sum. Original middle independence permits at most one such pair.
  Reversal gives the output conclusion, and the actual clean comparison
  then excludes negative singletons.

Thus neither a common active span nor positive singletons need be
assumed in advance. The complete Boolean action supplies the missing
constraint. The proof adds no native band, initialized helper,
approximation error or assumption on the unused clean block.

## 3. What this changes in the lower frontier

The exhaustive exact three-layer split is now:

| Actual source branch | Sign-degree ceiling |
| --- | --- |
| Either interior clean image is Pauli | Retained $`2^{26}`$ bound |
| Both are non-Pauli, with no original middle pair | Retained $`D_{\rm np}`$ bound |
| At least one original middle pair | Six, or five at common active rank at least four |

Here, with the constants defined in the retained proof homes,

```math
D_{\rm np}=\max\{D_3+576,D_{\rm ins}(0,0,0,3)\}.
```

Its finite-insertion term already dominates the Pauli-interior and
paired constants. Consequently every at-most-one-clean exact native-three
Boolean sign has degree at most $`D_{\rm np}`$, and every nonconstant
complete Boolean XOR predicate has degree at most $`D_{\rm np}-1`$.
Zero-clean words also have the much sharper retained five/four bound.

This closes the qualitative unrestricted one-clean three-layer degree
question. It does **not** replace the old large number by a small
universal six/five theorem: the no-pair proof still contains a very lossy
finite-feature estimate. Leaving its constant symbolic makes that
limitation explicit. A structural quantitative sharpening remains a
useful next lower question, not a reopening of the already proved
qualitative no-pair closure.

The retained approximation-to-exactness theorem applies first to the
original three-layer circuit. Its conservative radius
$`2^{-530690}`$ therefore gives the corresponding higher-degree
obstruction, allowing the retained common phase. It is not a claim at
arbitrary fixed tolerance. The resulting depth lower bound is still
constant; no dimension-growing unrestricted lower law follows.

## 4. Complete upper recurrence and the next research question

No new source was compiled. For the
[output-aliased weak source](OUTPUT_ALIASED_CHARACTER_SOURCE.md), let
$`m=2^b`$ and let eta bound full-operator error from the ideal
permutation with arbitrary quantum private-helper fibers. The paid
readout remains

```math
\begin{aligned}
T_I&\le4T_V,&D_T(I)&\le4d_V,\\
G_I&\le4G_V+bm+2b+2m+6,&h_I&\le m+h_V,\\
\varepsilon_I&\le4\eta,&a_{\rm clean}&=0.
\end{aligned}
```

These are actual source/adjoint calls. Their ideal private actions
cancel, and the full error includes arbitrary-input helper return and
reference behavior. Root fusion uses 32 calls, error at most
$`32\eta`$, and depth
$`32d_V(b_0+1)+O_\tau(\log\log(p+4))`$. Request
$`\eta=\tau/32`$; fixed-L frame transfer retains
$`\eta_L=2^{-L}/(512K_*)`$ and its literal allocations.

The distinct one-signal amplitude-$`1/\sqrt2`$ interface retains
depth $`24d_U+32`$ and complete error
$`16(3\sqrt2\epsilon+2\epsilon^2)`$. Its borrowed signal and
toggle return only approximately unless exact flatness is proved
separately. The two-band relative-gap theorem remains source-specific.

The next constructive question is still a paid weak quantum-fiber
source, coherent correction, paired transform, or complete extraction
with sublogarithmic native depth. The known source retains logarithmic
depth, so inserting it into either wrapper leaves the full resource
frontier unchanged. On the lower side, the original-pair scout has now
completed its exact three-layer branch; the next quantitative test is
whether the no-pair matrix core can be removed structurally, rather than
counted through another large alphabet. Higher layers require a new
argument.

The pair-collapse proof and its scalar-sector, reversal and coverage
steps were independently audited internally. Focused fixtures support
literal native identities and necessary hypotheses, not an enumeration
of the theorem. These checks do not replace external technical review.
