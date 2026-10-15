# Branch-orientation depth session

This session began on 15 October 2026 from commit `887dc8c`, after the
joint decoder/terminal session. The authorized ceiling was five hours.
The target was a complete fixed-accuracy circuit below R79's general
square-logarithmic depth at square-root dirty width, retaining one
external clean qubit, square-root T-count, linear Clifford count, literal
phases, all marker columns, and arbitrary dirty/reference inputs.

The unrestricted frontier remains unchanged. Two constructions survived
independent internal review:

- [R83, affine branch orientations](AFFINE_ORIENTATION_DEPTH.md), gives
  logarithmic depth when each depth has a common unordered pair of
  rotation magnitudes and its branch choices are affine Boolean
  functions of the prefix. Node signs and supplied leaf phases may vary.
  A native affine marker correction makes the complete factorization
  usable without a free permutation or selected-Clifford oracle.
- [R84, deterministic frames](DETERMINISTIC_FRAME_DEPTH.md), gives an
  exact real compiler for arbitrary node-dependent corner pairs. At
  square-root width its depth is
  $`O(\log(n+2)(1+\log^*(n+2)))`$. A direct last-deviation query erases
  the source tag; no approximate loader or amplification is needed.

[R82](BOUNDED_DENSITY_DEPTH.md) also extends from positive product
references to positive affine-oriented references using R83's priced
conditional core. Its fixed-density depth order and literal allocations
are unchanged. The geometric weighted correction never required product
structure; the native reference circuit did.

These are supporting-family results. They do not replace the selected
general theorem R79, implement a scalable general emitter, or claim a new
unrestricted lower bound. Final manuscript drafting was outside the task.

## 1. The chosen test family and the exact reduction

Fix $`0\lt a,b\lt1`$, $`a^2+b^2=1`$, and $`a\ne b`$. At every node v
an arbitrary bit $`r_v`$ chooses either $`(a,b)`$ or $`(b,a)`$.
This retains one bias while allowing arbitrary nonlinear branch choices.
Depth-dependent base pairs and zeros are included in the exact algebra.

Define the tree automorphism recursively by

```math
F(\emptyset)=\emptyset,\qquad
F(yi)=F(y)(i\mathbin\oplus r_{F(y)}).
```

Let P permute leaves by F. For the prescribed marker
$`\lambda(v)=v1\,0^{n-|v|-1}`$, put

```math
Q|0\rangle=|0\rangle,\qquad
Q|\lambda(v)\rangle=(-1)^{r_v}|\lambda(F^{-1}(v))\rangle.
```

The [full-column argument](AFFINE_ORIENTATION_DEPTH.md#3-full-column-leaf-and-marker-factorization)
does not need affinity: P transports each reference scaling vector to
its oriented counterpart, and reversing the two children contributes
exactly the displayed marker sign. Thus

```math
W_r=P W_0 Q.
```

Here $`W_0`$ is the homogeneous reference with the unswapped pairs.
Both permutations must be priced, or their joint word must be compiled
directly. A state-column identity alone would miss Q and its literal
signs. No division by an incoming amplitude is used.

At the positive balanced frame B the branch choices disappear. Therefore

```math
B=P B Q,\qquad Q=B^\dagger P^{-1}B,
```

```math
W_rB^\dagger=P(W_0B^\dagger)P^{-1},\qquad
W_r-B=P(W_0-B)Q.
```

The ordinary norm of the last difference and the norm of its entrywise
absolute matrix are both orientation-independent. This is useful
geometry, not a circuit for P. In the corner case $`W_0=I`$, the cheap
deterministic frame is $`C=P Q`$, which still differs from P.

## 2. A cheap XOR oracle does not complete the permutation

Offline expansion gives a triangular map

```math
F(x)_j=x_j\mathbin\oplus f_j(x_{\lt j}).
```

Its inverse has the same form. The diagonal on two arbitrary n-bit
registers

```math
D_F|x,z\rangle=(-1)^{z\cdot F(x)}|x,z\rangle
```

has an exact native square-root-width implementation. The terms
$`z_jx_j`$ are Clifford CZ gates. Each remaining term
$`z_jf_j(x_{\lt j})`$ is a Boolean phase on j address bits. Give each
term private dirty query work and use the
[shallow returned-dirty phase interface](SHALLOW_DIRTY_QUERIES.md).
Shared address interactions are read-only Clifford controls, so the
physical non-Clifford supports run in parallel. Geometric table volumes
give

```math
T,w=O\!\left(\sum_{j=1}^n2^{j/2}\right)=O(\sqrt N),\qquad
G=O\!\left(\sum_{j=1}^n2^j\right)=O(N),
```

```math
D_T=O(\log(n+2)).
```

The public second word fits in the stated width above a fixed dimension
cutoff. All helpers return, including on entangled inputs. Conjugating
by Hadamards on z gives the complete XOR oracle
$`|x,z\rangle\mapsto|x,z\mathbin\oplus F(x)\rangle`$.

This does not give $`P\otimes I_{\rm dirty}`$. With arbitrary initial
y, the usual compute/swap/inverse-query word instead gives

```math
(x,y)\longmapsto
\bigl(y\mathbin\oplus F(x),\,
x\mathbin\oplus F^{-1}(y\mathbin\oplus F(x))\bigr).
```

The second word need not return. A simultaneous dirty-mask echo also
changes its own query addresses. For example, let
$`f_1=1,f_2(x_1)=x_1`$, start with $`x=d=00`$, query all f into d,
XOR d into x, repeat the query, and repeat the XOR. In first-bit-first
notation the resulting pair is $`(11,01)`$ instead of $`(10,00)`$.

An initialized output word would repair the compute/swap/uncompute
construction. Protecting h such bits inside a complete frame, however,
requires a remaining logical suffix of at least $`h+O(\log h)`$.
The resulting recurrence halves the local height. It retains
$`\Theta(\log n)`$ stages, each with a global-address logarithmic
query charge. Even eliminating that address charge alone would leave
the ordinary biased-reference synthesis charges
$`\sum_i\log m_i=\Theta(\log^2 n)`$ for this halving schedule.
This benchmark preserves the old order and was stopped.

Fourier overlap supplies no free repair: with initialized output, its
direct block amplitude for the permutation is $`N^{-1/2}`$. Coherent
amplification would require $`\Theta(\sqrt N)`$ oracle appearances,
already raising the T count to $`O(N)`$ for this construction. These
are failures of specified conversions, not a general impossibility
theorem for native triangular permutations.

## 3. What the two successful mechanisms change

For affine orientation, $`F(x)=Lx+v`$ is unitriangular affine. R83
constructs a constant-T-depth outer guard by expressing its linear part
as two involutions, each Clifford-conjugate to disjoint CNOT pairs.
It then corrects the marker tail after the full inverse map. The tail
correction is an affine shear whose coefficients depend on a short
last-one tag. Selected suffix Hadamards turn it into quadratic Boolean
phases with short addresses and private dirty pools. Long outside-prefix
dependence contributes polynomially many separate terms, not a long
lookup address. The returned tag is reused by the direct homogeneous
common-flag core. This is why the recursive depth charges sum to
$`O_L(\log n)`$.

For deterministic corners, each complete column is a signed basis leaf.
Given its marker tag, a selected-suffix XOR prepares the precomputed
favored continuation. A phase query between selected Hadamard banks
implements that many-output XOR without an initialized output word.
For a resulting leaf x define its deviations

```math
z_j(x)=x_j\mathbin\oplus r_{x_{\lt j}},\qquad
t_{\rm out}(x)=\max\{j:z_j(x)=1\},
```

with empty maximum zero. The old tag equals this single classical
function of the output. One returned-dirty word query therefore erases
it exactly. Protected logarithmic suffixes give an iterated-logarithmic
number of stages, each with constant global-query multiplicity. R84
prices the tag-expanded tables and supplies an exact all-width fallback.

The same output-only eraser cannot interpolate to positive bias. Suppose
a processor preserves computational x and acts only on tag/work. If two
marker columns have nonzero amplitude at the same x with distinct tags,
their tag/work states are orthogonal. A unitary cannot send both to the
same zero tag and returned-work state. More quantitatively, for two
orthogonal inputs at most one can have overlap squared greater than
one half with that required output. For $`n\ge2`$, last-layer and parent markers share
an x with amplitudes at least $`\min(a,b)`$ and
$`\min(a,b)^2`$. One of their column errors is consequently at least
$`\min(a,b)^2/\sqrt2`$. A basis change acting on the data, or a
different recovery mechanism, is necessary. This restriction concerns
tag/work-only erasure, not arbitrary complete circuits.

## 4. The balanced positive correction has growing normalization

A constant-normalization positive entry-overlap correction to B does
not cover the fixed biased family. Let

```math
\psi=(a|0\rangle+b|1\rangle)^{\otimes n},\qquad
c=(a+b)/\sqrt2\lt1.
```

For a marker at depth d with homogeneous prefix amplitude $`A_p`$,
the homogeneous and balanced columns have the same entry signs. Their
absolute overlaps satisfy

```math
\langle\psi,|e_0|\rangle=2ab A_p,\qquad
\langle\psi,|e_B|\rangle=c^{n-d}A_p.
```

The reverse triangle inequality and
$`\sum_{|p|=d}A_p^2=1`$ give

```math
\bigl\|\,|W_0-B|\,\bigr\|^2
\ge\sum_{m=1}^n\max(0,2ab-c^m)^2
=4a^2b^2n-O_{a,b}(1).
```

The omitted initial terms are a fixed finite set; the remaining error
terms are bounded geometric sums. For an upper bound, partition a
frame's columns into its state column and n marker-depth blocks. Within
one depth their absolute columns have disjoint supports and unit norm.
Cauchy-Schwarz therefore gives
$`\|\,|W|\,\|\le\sqrt{n+1}`$. Thus

```math
\bigl\|\,|W_r-B|\,\bigr\|=\Theta_{a,b}(\sqrt n)
```

for every orientation table, by Section 1's signed-permutation identity.
In the [entry-indexed support model](ENTRY_INDEXED_ENCODING_LIMIT.md),
taking absolute values of the row/column sources retains contractions
and forces normalization at least the norm of the absolute represented
matrix. Hence an exact positive single-overlap encoding of this
correction cannot have dimension-independent normalization. This says
nothing against interference architectures outside that entry-indexed
support model, other reference families, or unrestricted native circuits.

## 5. A short balanced-Haar word cannot supply the missing basis change

There was a concrete possible use of the cheap corner frame C. With
$`C=P Q`$ and B balanced, Section 1 gives

```math
P B P^{-1}=B C^{-1},\qquad
P K P^{-1}=C K C^\dagger
```

whenever K is diagonal by marker depth, since K commutes with Q. A short
word for $`W_0`$ in B, its inverse, and such diagonals would then have
a priced conjugate. The following bound rules out that proposed short
word, even at fixed error and fixed strictly interior bias.

Write $`q=|+\rangle`$. The entrywise absolute balanced frame has the
exact decomposition

```math
|B|=K_0+\sum_{t=1}^nK_t,\qquad
K_0=|q^n\rangle\langle0^n|,
```

```math
K_t=I^{\otimes(t-1)}\otimes
|q^{n-t+1}\rangle\langle1\,0^{n-t}|.
```

Each $`K_t`$ or its adjoint sends a tensor product of single-qubit
states from $`\{|0\rangle,|1\rangle,q\}`$ to another such product,
times a scalar in $`[0,1]`$, or to zero. Consider an n-logical-qubit
word

```math
U=D_0B^{s_1}D_1\cdots B^{s_k}D_k,\qquad s_i\in\{1,-1\},
```

where all D are arbitrary unitary computational diagonals. With $`e_0=|0^n\rangle`$,
entrywise absolute values remove these phases and give

```math
|Ue_0|\le |B^{s_1}|\cdots|B^{s_k}|e_0.
```

Expanding the right side gives at most $`(n+1)^k`$ of the
product states just described, each with coefficient at most one.

Put

```math
\rho=\max\{a,b,(a+b)/\sqrt2\}\lt1.
```

Every expanded state has overlap at most $`\rho^n`$ with the biased
root column $`\psi`$. If $`\|U-W_0\|\le\varepsilon\lt1`$, then

```math
1-\varepsilon
\le|\langle\psi|U|0^n\rangle|
\le(n+1)^k\rho^n.
```

Consequently

```math
k\ge\frac{n\log(1/\rho)+\log(1-\varepsilon)}{\log(n+1)}
=\Omega_{a,b,\varepsilon}(n/\log n).
```

The orientations of B and its inverse need not alternate. Arbitrary
diagonal phases are allowed, so depth-only phases are covered. This is
a word-alphabet restriction on the n logical qubits. Ancillary
operations, non-diagonal embeddings, enlarged-register mixers,
additional Clifford changes of basis, and amortized implementations
of growing words lie outside its hypotheses. It is not a T-depth lower
bound for the biased frame: R80 already compiles that frame in
logarithmic depth using other native operations.

There is a simpler separate restriction on bare corner words. Across
the first-bit cut, a corner frame or its inverse has off-diagonal block
rank at most one: its computational-label tree word uses the sole
crossing edge at most once. Off-diagonal rank is subadditive under
products. A first-bit flip has crossing rank $`N/2`$, so it needs at
least $`N/2`$ such bare factors in that fixed basis. The flip itself
is Clifford. This argument therefore expressly permits no conclusion
about corner words interleaved with Clifford basis changes.

## 6. Bounded product patches are not a general partition rule

Take $`a>b>0`$ and let the orientation at a nonroot node be its last
physical prefix bit, with root orientation zero. Put
$`p(0)=a^2,p(1)=b^2`$. Its root probabilities are

```math
\mu(x)=p(x_1)\prod_{j=2}^hp(x_j\mathbin\oplus x_{j-1}).
```

In any conditioned h-level block, compare the uniform measure on the
two constant strings with the uniform measure on the two alternating
strings. Both have the same one-bit marginals. Their expected log
density under any positive product q is therefore equal. The first
transition contributes equally to expected log mu, including a fixed
outside boundary bit. The remaining transitions give a difference
$`(h-1)\log(a^2/b^2)`$. Since both expectations of log(mu/q) lie
between its minimum and maximum,

```math
\frac{\max_x(\mu_x/q_x)}{\min_x(\mu_x/q_x)}
\ge(a^2/b^2)^{h-1}.
```

Fixed positive ratio bounds allow only bounded-height product-reference
patches in every conditional block. Merely concatenating such patches
leaves linearly many stages. This example is nevertheless easy for
R83: its leaf map is linear prefix XOR. It motivates the proved affine
reference extension of R82 and shows why failure of a reference
restriction is not intrinsic frame hardness.

Pruning tiny-probability branches does not repair the full-frame error
contract. Below an arbitrarily rare prefix, changing one orientation
changes the normalized marker by $`\sqrt2|a-b|`$, independently of
the incoming mass. The root-state change may be tiny while that marker
column has constant error.

## 7. Remaining mechanism and stopping point

The exact bipartite phase/XOR interface in Section 2 is already within
the desired count and depth budgets. What remains missing is a native
conversion to an arbitrary in-place triangular permutation with one
external clean qubit and returned dirty work, or a direct compiler for
$`P W_0 Q`$ that avoids those separate permutations. Such a proposal
must also price marker recovery, all actual inverses, physical T-layer
supports, and a full initialized-isometry error bound.

The corner result removes recovery only where the output determines its
tag. The affine result removes the nonlinear in-place obstacle only for
a class with a native guarded-linear-map schedule. The balanced positive
correction, short Haar-word factorization, dirty-output conversion, and
halving recurrence above do not bridge the arbitrary fixed-bias case.
A corner/Clifford construction outside the proved word restrictions
remains possible, but no qualifying identity was obtained in this
session. Repeating those failed interfaces without changing an assumption
is not a further research step.

## 8. Review and verification

Separate internal lanes derived and cross-checked the affine marker
correction, guarded unipotent maps, corner source and decoder, all-width
recurrences, and the scoped obstructions. The native proof precedes the
bounded fixtures in
[`tests/test_orientation_depth.py`](../../tests/test_orientation_depth.py).
Those checks exercise literal algebra and work return; they do not
establish asymptotic resource bounds or implement a general emitter.
Internal audits are not external peer review or a priority assessment.

The reviewed pull request records the complete validation, unchanged
fault-tolerant receipts, rendered-document checks, and remote CI results.
