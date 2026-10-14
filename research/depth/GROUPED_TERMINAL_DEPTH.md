# Grouped terminal layers with whole-word coarse programs

[Precision-uniform frame theorem](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) · [Conditional suffix groups](../../docs/CONDITIONAL_SUFFIX_COMPILER.md) · [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md)

**Retained preceding construction.** The current precision-uniform
theorem uses the [shared binary terminal source](SHARED_BINARY_TERMINAL.md),
which replaces repeated fine-source preparation by one source and retained
failure history. The three-halves-power bound below remains a standalone
result. Its exact whole-word coarse interpreter is also a component of
the current proof.

Several terminal layers can share one fine-precision scalar source. The
retained column-forest construction already provides that sharing. Its
coarse circuit becomes shallower if each short native word is fetched
once into a conditional-zero program buffer, interpreted in its original
order, and unloaded by the actual query inverse. Reserving enough lower
logical wires for that buffer reduces the terminal precision term from
quadratic to three-halves power.

The argument uses the individual-group proof in Sections 1–5 of the
[conditional-suffix compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md),
not its particular exponentially growing partition. All one-qubit
synthesis constants, exact query interfaces, and native predicate
constants below are absolute. Effectively specified frame angles and
certified classical coefficient evaluation are inherited assumptions.
Classical preprocessing is separate from the native circuit resources.

## 1. Terminal lemma and partition

Write

```math
N=2^n,\qquad \ell=\log_2(n+2),\qquad p=L+\ell,
\qquad \eta_L=2^{-L}.
```

Fix the absolute constants A and c of the
[protected-suffix cutoff](SUBLINEAR_T_DEPTH.md#102-a-protected-suffix-and-an-absolute-dimension-cutoff),
and let

```math
6\le L\le cn,\qquad
M=\left\lceil A\left[L+\log_2(n+L+2)\right]\right\rceil.
```

There are absolute positive constants $`A_0,A_1,n_0`$ such that, for
$`n\ge n_0`$, $`M\le n/2`$, and

```math
a=2,\qquad b\ge17(L+n+7),\qquad b\ge\sqrt{N/n},
```

the final M layers of any prescribed complete real Hopf frame have one
coherent Clifford+T circuit with complete-isometry error less than
$`\eta_L/2`$ and simultaneous bounds

```math
T_T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
G_T=O(NL),
```

```math
D_{T,T}=O\!\left(\frac{NL}{b^2}+p^{3/2}+p\ell\right).
```

Every constant is independent of n, L, b, and the frame. The error
includes both external flags, logical workspace, arbitrary dirty inputs,
and arbitrary external references. No initialized work beyond the two
external flags is supplied: other initialized bits are logical suffix
bits used only on their active zero sector. There are no measurements,
resets, supplied resource states, or uncharged quantum oracles.

For the construction, set

```math
s=\lceil\sqrt p\rceil,\qquad
R_0=\min\{M,\lceil A_0\ell+A_1s\rceil\}.
```

Leave the deepest $`R_0`$ original layers to the retained capped
layerwise compiler. Partition the other $`M-R_0`$ layers into
consecutive groups of heights $`1\le t\le s`$, taking full-height
groups except possibly one. Execute groups in the prescribed
shallow-to-deep order, followed by the deepest layers. Their number g
and heights satisfy

```math
g\le M/s+1,\qquad \sum_{\rm groups}t\le M=O(p).
```

A group ending at absolute depth e has $`r=n-e\ge R_0`$ untouched
lower logical bits. The first external clean qubit saves the predicate
$`H=[\text{these }r\text{ bits are zero}]`$. The second is the
scalar-source flag. Every completed group subroutine preserves H and
is exactly identity on H zero, on the whole space of its other inputs.

## 2. Prefetch the exact coarse native word

For a t-layer group, the retained proof chooses

```math
K=2^{\lceil\log_2(8t+12)\rceil}=O(t),\qquad
\delta_c=\frac1{4tK\,2^{t/2}}.
```

Phase-calibrated one-qubit synthesis at accuracy $`\delta_c`$ gives
native words of length $`w=O(\log(1/\delta_c))=O(t)`$ over a fixed
finite alphabet. Pad them with identity symbols to a common length.
Their exact addressed product is the coarse circuit C. Its literal
native phases are retained, including phases of complex coarse words.
The classical residual is evaluated from this actual C:

```math
E=C^\dagger W_g-I.
```

At one local depth, replace the old symbol-at-a-time stream by:

1. Compute the predicate g that the remaining **local logical bits**
   below this word's target are zero.
2. Query the entire $`O(t)`$-bit encoded native word into a zero program
   buffer. Include H in the address and assign literal zero rows to H
   zero.
3. Interpret the same native symbols in their original order. Each
   literal alphabet gate has direct H and g controls, as well as its
   bounded symbol controls.
4. Apply the actual full-word query inverse, then the actual predicate
   inverse.

The program address consists of the preserved external prefix and the
preceding local bits; it excludes the target. The local predicate also
excludes the target. Both addresses remain unchanged throughout the
one-qubit word. In particular, the predicate does **not** test the
protected lower r-bit suffix, which now holds the program and other
private work.

Every completed symbol preserves each program bit exactly; any basis
changes for negative symbol literals are reversed within that symbol.
Consequently the program and predicate return exactly, even on occupied inputs.
On their initialized active inputs, the resulting logical operation
is precisely the chosen coarse word. An occupied active program can
select a different word; no desired-word approximation on that input
is claimed or needed.

Direct H control on every interpreted symbol is essential. Zero query
rows alone would not suppress interpretation of an occupied inactive
program. With the direct control, the interpreter is identity when H
is zero, and its surrounding query and predicate words cancel exactly.
Thus the complete coarse operation is identity on the entire inactive
space, including occupied program and predicate registers.

The retained
[bounded-control native lemma](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#exact-conditioning-without-a-clifford-count-penalty)
implements these controls with constant exact cost per symbol and a
fixed pool of returned dirty helpers. In particular, controlled T and
controlled phases are literal native phase echoes with their residual
corrections; no projective gate replacement is made. Helpers are
disjoint from controls, targets, and live program bits.

The original residual uses at most
$`A_{\rm priv}\log_2(t+2)+B_{\rm priv}`$ initialized internal bits.
The complete program adds at most $`C_{\rm buf}t`$ bits. Choose the
absolute $`A_1,A_0`$ so that

```math
C_{\rm buf}t+A_{\rm priv}\log_2(t+2)+B_{\rm priv}
\le A_1s+A_0\ell\le r
```

whenever middle groups exist. The program, label, rejection flags,
and predicate scratch have disjoint live allocations. A constant
temporary and the returned dirty helper pool can be reused only after
their actual return.

## 3. The residual and its native depth

The retained column-forest decomposition has $`O(t)`$ term labels and
$`O(t2^e)`$ coefficient rows. Its forward and reverse atoms use
distinct scalar and atom rejection flags. A common SELECT uses exactly
two shared scalar blocks, one forward and one inverse, independently
of t. It does not invoke a separate precision source for each depth
label. Literal phase labels remain outside the forward or inverse
real-base word, as in the retained proof.

Use one uniform scalar-source width for all middle groups:

```math
m=L+\lceil\log_2(n+2)\rceil+10=O(p).
```

The source core is arbitrary dirty work. Its coefficient query is an
exact XOR into that core, with no initialized m-bit precision output.
Source masks use the retained exact all-dirty operator source and its
actual inverse, scheduled serially in $`O(m)`$ T-depth. Each large
coefficient query and inverse completes before a following atom changes
its local logical address. In the reverse branch the scalar word comes
first and the inverse atom follows, so its query is also unloaded
before the address changes.

The [shallow word-query lemma](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address)
implements every coefficient query and every whole coarse-word query
in $`O(\ell)`$ T-depth, including polynomially many output bits.
Its public address appears only as a CNOT control; all simultaneous
non-Clifford supports for different outputs are private. Counts and
width are charged in Section 5 below.

At each of the t coarse depths, loading and unloading cost
$`O(\ell)`$, the one local-suffix predicate and its inverse cost
$`O(\ell)`$, and interpreting the $`O(t)`$ literal native symbols
costs $`O(t)`$. Therefore

```math
D_{T,\mathrm{coarse}}=O(t^2+t\ell).
```

For the atom, use the explicit
[selected-type schedule](../../docs/T_DEPTH_COMPILER.md#group-schedule).
Unroll its $`O(t)`$ labels. For each label, decode a type into a
conditional-zero temporary, toggle the atom rejection flag using a
constant number of marker conjunctions on at most $`O(n)`$ controls,
and erase the predicate scratch **before** changing the local word.
Apply at most t bounded-controlled Hadamards and one controlled X
serially, then erase the type temporary from the unchanged term label.
Negated markers use a selected-type toggle followed by a selected-type
positive-marker toggle. Exact read-only conjunctions have
$`O(\ell)`$ T-depth and returned arbitrary dirty helpers. The reverse
atom is the actual reversed native word. Hence

```math
D_{T,\mathrm{atom}}=O(t^2+t\ell).
```

This bound uses no simultaneous shared-control Hadamard premise and
no unspecified polynomial in t. Group activity predicates and the
outer initialized-work reflections each add $`O(\ell)`$ depth.
All source, atom, query, and reflection multiplicities in the
normalization-two amplification are absolute constants. Thus

```math
D_{T,g}=O(p+t^2+t\ell).
```

Summing over groups, using $`t\le s`$ and $`M=O(p)`$, gives

```math
D_{T,\mathrm{middle}}
=O\!\left(p(M/s+1)+Ms+M\ell\right)
=O\!\left(p^2/s+ps+p\ell\right).
```

## 4. Complete-input error and workspace return

The unchanged group accepted-block proof gives scalar coefficient
error at most $`(5/4)2^{1-m}`$ and, after normalization-two
amplification, complete-isometry error at most

```math
10\,2^{-m}.
```

Neither the column length nor its classical square-root coefficient
scaling multiplies this error. The bound holds uniformly over the
preserved prefix direct sum and all arbitrary dirty inputs and
references. With at most n groups,

```math
\sum_{\rm groups}10\,2^{-m}
\le\frac{10n}{1024\,2^{\lceil\log_2(n+2)\rceil}}\eta_L
\lt\frac{10}{1024}\eta_L.
```

The initialized-work projector used by the group amplification includes
the new program buffer. Equivalently, retain the original reflection
testing the entire active suffix and the scalar flag. The buffer is
idle during the residual half-block and starts at zero in its
initialized-work embedding. The physical amplification keeps its
guarded sign $`Z_H`$; a global minus sign would not preserve identity
on H zero.

Apply the actual coarse C after the residual amplification and then
reverse the original computation of H. On ideal initialized columns C
returns its work exactly. Residual rejection leakage need not have
returned the entire suffix exactly; the actual coarse circuit and the
actual predicate inverse simply propagate that existing error by
unitaries. No reset, discarded rejection space, or occupied-program
accuracy assumption is used.

There is no separate final coarse-approximation error: the residual
was computed from the exact native C, and
$`C(C^\dagger W_g)=W_g`$. Coarse accuracy is used to bound the
residual coefficients, not as an additional approximation after their
correction. The complete group is exactly identity on its original
inactive sector. These properties permit the same two external flags
and dirty pool to be reused between groups by full-isometry
telescoping.

## 5. Absolute native counts and dirty width

First suppose middle groups exist, so that
$`R_0=\lceil A_0\ell+A_1s\rceil`$. In the stated precision range,
m, t, and every address length are $`O(n)`$ with absolute constants.
Every coefficient table has $`O(t2^{n-r})`$ rows; a coarse word table
at any local depth has at most $`O(2^{n-r})`$ rows. Padding, phase
labels, and activity guards give the common row bound

```math
Q_{\rm table}\le C(n+2)^2N2^{-R_0},
```

with at most $`C(n+2)`$ address and output bits. There are at most
$`C(n+2)`$ completed long-query occurrences: a constant number per
residual group and a constant number per coarse local depth, rather
than one per native symbol.

For a q-bit address, m outputs, and Q rows, the conservative query
lemma bounds T-count and returned dirty width by
$`O(m(q+2)^3\sqrt Q)`$, and Clifford count by
$`O(mQ+m(q+2)^3\sqrt Q)`$. Substitution and summation above fit
degree six in n. A deliberately larger common majorant is convenient:

```math
T_{\rm middle},w_{\rm middle}
\le C(n+2)^{10}\left(\sqrt N\,2^{-R_0/2}+1\right),
```

```math
G_{\rm middle}
\le C(n+2)^{10}\left(N2^{-R_0}+1\right).
```

The additive polynomial charges all local source words, native
controlled symbols, label decoding, activity and marker predicates,
reflections, and their returned dirty helpers. Indeed, there are at
most $`O(n^2)`$ local symbols or atom gates and at most $`O(n^2)`$
predicate occurrences, each with the retained safe cubic count and
width majorant. These fit well within degree ten. Peak dirty width is
reused after actual return; multiplying it by the query count above
only overestimates its requirement. No private initialized logical bit
is also allocated as an arbitrary dirty bank.

Choose $`A_0\ge24`$ in addition to the conditional-work inequalities.
Then

```math
(n+2)^{10}2^{-R_0/2}
\le(n+2)^{10-A_0/2}\le(n+2)^{-2}.
```

The analogous Clifford factor has exponent $`10-A_0`$. Increasing an
absolute dimension cutoff absorbs the additive degree-ten terms into
the exponential budgets, and gives

```math
w_{\rm middle}=o(\sqrt{N/n}),\qquad
T_{\rm middle}=O(\sqrt N),\qquad G_{\rm middle}=O(N).
```

These estimates are uniform over $`6\le L\le cn`$. In particular,
the returned query/helper pools fit at $`b\ge\sqrt{N/n}`$ above
one absolute cutoff; none is a precision-dependent hidden constant.

The original linear core/selector reservation can also be retained.
Writing $`B_0=L+n+7`$, a coefficient-query address has length
$`k=n-r+O(\log_2(t+2))`$, including H and the mode bits. Thus

```math
m+k+O(1)
\le L+\lceil\ell\rceil+10+n-r+O(\log_2(t+2))+O(1)
\le B_0
```

after increasing the same fixed $`A_0,A_1`$. The saved suffix supplies
initialized work; the m-bit core and k selectors are external arbitrary
dirty work. Additional shallow-query banks and predicate helpers are
paid in the returned pool above. Disjoint pools are reserved while
simultaneously live; all later reuse follows an exact completed return.

## 6. The deepest layers and the terminal bound

For the deepest $`R_0`$ layers retain the capped precision schedule
and the nonuniform rectangular-query allocation from the
[layerwise tail bound](SUBLINEAR_T_DEPTH.md#7-the-final-m-layers-retain-the-count-budgets).
If $`1\le j\le R_0`$ is a layer's height above the leaves, with
$`j-1`$ lower logical bits, this is

```math
m_j=L+6+\min\{j,\lceil\log_2(8n)\rceil\},\qquad
Q_j=4N2^{-j}.
```

Restricting its weighted sums to a subset of layers preserves

```math
T_{\rm bottom}=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),
\qquad G_{\rm bottom}=O(NL).
```

These are uniform bounds: $`m_j\le L+6+j`$ controls the convergent
weighted geometric sums, and the precision cap gives $`m_j=O(p)`$.
Hence

```math
\sum_{j=1}^{R_0}Q_jm_j=O(NL),\qquad
\sum_{j=1}^{R_0}\sqrt{Q_jm_j}=O(\sqrt{NL}),\qquad
\sum_{j=1}^{R_0}m_j=O(R_0p).
```

with absolute constants for $`L\ge6`$. The retained bottom base is
$`L+n+9`$, allowing its extra two error-budget bits. Above an absolute
cutoff the returned pool $`B=b-(L+n+9)`$ satisfies $`B\ge b/2`$
and $`B\ge16c_I(n+3)`$, where $`c_I`$ is the retained absolute
indicator constant, since $`L\le cn`$ and
$`b\ge\sqrt{N/n}`$. Polynomial native predicate and source counts
fit the same exponential budgets after increasing that cutoff. No
fixed-accuracy hidden constant is used in these weighted sums.

Its source widths are $`O(p)`$. The same query allocation and returned
predicate helpers give

```math
D_{T,\rm bottom}
=O\!\left(\frac{NL}{b^2}+R_0p\right)
=O\!\left(\frac{NL}{b^2}+p(s+\ell)\right).
```

The capped source error over all deepest layers is less than
$`\eta_L/4`$. Together with Section 4 this gives terminal error
less than $`\eta_L/2`$. The bottom compiler and the middle groups
run sequentially, so their maximum dirty allocation, rather than its
sum, is required.

If $`R_0=M`$, there are no middle groups and no middle query ledger
is needed. In that case
$`M\le\lceil A_0\ell+A_1s\rceil=O(\ell+s)`$, so the same
bottom depth bound still holds. Thus this case does not silently
reintroduce a quadratic precision term.

Combining the two portions proves

```math
D_{T,T}
=O\!\left(\frac{NL}{b^2}+p^2/s+ps+p\ell\right)
=O\!\left(\frac{NL}{b^2}+p^{3/2}+p\ell\right),
```

with the simultaneous T and G bounds in Section 1.

## 7. Consequence for the full precision-uniform frame

The [shared-gradient early circuit](SHARED_GRADIENT_RECOVERY.md) has
depth $`O(L+\sqrt n\,\ell)`$ and complete error at most
$`9\eta_L/64`$. Its dirty work is sequential with this terminal
circuit. Combining the two errors stays below $`\eta_L`$, while

```math
p^{3/2}=O(L^{3/2}+\ell^{3/2}),\qquad
L\ell\le L^{3/2}+\ell^3,\qquad
\ell^3=O(\sqrt n\,\ell)
```

reduces their additive depth to
$`O(\sqrt n\,\ell+L^{3/2})`$. The middle inequality follows by
splitting at $`L=\ell^2`$; the last has an absolute constant for
$`n\ge1`$.

The [uniform width and finite-size splice](SUBLINEAR_T_DEPTH.md#104-every-eligible-width-and-a-uniform-finite-fallback)
uses the retained compiler when $`b\lt\sqrt{N/n}`$, because then
$`nL\le NL/b^2`$, and below the enlarged absolute dimension cutoff.
Selecting also against that retained circuit supplies the global bound

```math
D_T=O\!\left(\frac{NL}{b^2}
 +\min\{nL,\sqrt n\log_2(n+2)+L^{3/2}\}\right)
```

with the same literal base $`b\ge17(L+n+7)`$, two flags, T-count,
and Clifford count. This is a constructive upper bound, not a lower
bound on possible faster source synthesis or other grouped circuits.
The $`L\le cn`$ restriction and the source/preprocessing assumptions
are unchanged.

The associated [bounded grouped-terminal fixtures](../../tests/test_grouped_terminal_depth.py)
check native program interpretation, literal phases, actual unloading,
occupied inactive work, and rounded group partitions. Their scope is
finite interface evidence; the absolute cutoffs, polynomial ledgers,
and uniform resource bounds are the analytic argument above.
