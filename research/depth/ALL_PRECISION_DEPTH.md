# Optimal count and a depth bound at every precision

[Resource statement](../../docs/OPEN_PROBLEM.md) · [Linear-precision compressed columns](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) · [Whole-frame binary stream](WHOLE_FRAME_BINARY_DEPTH.md) · [High-precision depth](HIGH_PRECISION_DEPTH.md)

This chapter combines the complete-frame constructions into one theorem
for every dimension, precision, and eligible dirty width. The whole-frame
binary stream and the high-precision construction need not have overlapping
precision ranges: in the intervening range, the retained layerwise depth
already fits the logarithmic envelope. The combination introduces no new
native primitive or source-state premise. Its stable result identifier is
**R70**.

## 1. The uniform theorem

Let $`n\ge1`$ and $`L\ge6`$ be integers, and put

```math
N=2^n,\qquad \eta=2^{-L},\qquad B_0=L+n+7,
```

```math
F(n,L)=\min\left\{nL,(\sqrt n+L)\log_2(n+L+2)\right\}.
```

For every integer $`b\ge17B_0`$ and every prescribed complete real
Hopf frame W, there is one coherent Clifford+T circuit V with two
external clean flags and at most b arbitrary dirty qubits such that

```math
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+F(n,L)\right).
```

Here $`J_2`$ initializes the two external flags and leaves the logical
and dirty inputs unrestricted. All constants are absolute and all three
resource bounds hold on the same selected circuit. The norm includes
every logical input, arbitrary dirty inputs and reference correlations,
literal phases, work return, and clean or logical-work leakage. Unused
supplied dirty wires remain untouched. G counts elementary Clifford
gates; arbitrary intervening Clifford circuits are permitted in T-depth
and remain charged in G.

Angles admit terminating certified evaluation, including singular and
unequal tuples. Classical preprocessing is terminating, without an
efficiency claim. Every source preparation, exact query, actual inverse,
guard, and recovery is included. The construction uses no measurement,
reset, supplied precision state, QRAM, or uncharged coherent evaluator.

There is a stronger bound throughout the high-precision range:

```math
L\gt N/n^2
\quad\Longrightarrow\quad
D_T=O\!\left(\frac{NL}{b^2}+L\right),
```

with the same count, Clifford-count, width, and error guarantees.
Equivalently, the depth theorem can use the sharper piecewise term

```math
H(n,L)=
\begin{cases}
F(n,L),&L\le N/n^2,\\
L,&L\gt N/n^2.
\end{cases}
```

The inherited worst-case T-count lower bound matches the displayed
count at every eligible width. The depth statement is an upper bound;
Section 4 gives its matching width intervals and their limits.

## 2. Proof dependencies

| Complete-circuit input | Precision range used here | Role in this proof |
|---|---|---|
| [Retained uniform compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary) | Every $`n\ge1,L\ge6`$ | At every $`b\ge17B_0`$, gives $`T=O(\sqrt{NL}+NL/b+nL)`$, $`G=O(NL)`$, and $`D_T=O(NL/b^2+nL)`$; supplies the intermediate-precision and finite-dimension branches |
| [Compressed-column theorem, R66](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) | $`6\le L\le cn`$, for an absolute $`c\gt0`$ | Optimal count and additive depth $`\min\{nL,(\sqrt n+L)\log_2(n+2)\}`$ |
| [Whole-frame binary stream, R68](WHOLE_FRAME_BINARY_DEPTH.md) | $`6\le L\le2^{\kappa n}`$, for an absolute $`0\lt\kappa\lt1`$ | Optimal count and additive depth $`\min\{nL,(n+L)\log_2(n+L+2)\}`$ |
| [High-precision theorem, R69](HIGH_PRECISION_DEPTH.md) | $`L\gt N/n^2`$ | Optimal count and additive depth L |
| [Inherited complete-frame lower bounds](../../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage) | Every dimension and precision at physical width $`\Theta(b)`$ | Count optimality and the matching depth interval |

The three refined constructions use the same two-clean full-isometry
contract and the same literal lower threshold $`b\ge17B_0`$.
Their depth improvements are not inferred from their count bounds.
Fix their absolute constants c and kappa and one absolute dimension
cutoff $`n_0`$ large enough for all their asymptotic branches.
The retained compiler handles smaller dimensions directly.

## 3. Select one complete circuit

First observe the uniform inequalities

```math
L\le F(n,L)\le nL,\qquad L\le H(n,L).
```

The first lower bound follows from $`n\ge1`$ and
$`\log_2(n+L+2)\gt1`$. The following branch order avoids any
assumption about the relative order of the precision cutoffs at small n.

### Bounded dimensions, with precision unrestricted

For $`1\le n\lt n_0`$, use the retained uniform compiler at the
supplied width. Its additional count and depth term satisfies

```math
nL\le n_0L\le n_0H(n,L).
```

Consequently its count is
$`O(\sqrt{NL}+L+NL/b)`$, its depth is
$`O(NL/b^2+H(n,L))`$, and its Clifford count is $`O(NL)`$.
The enlarged constant depends only on the single fixed absolute cutoff
$`n_0`$, not on L or b. This is the actual all-precision compiler
with its original literal dirty allocation. No finite-accuracy net is
being extended to unbounded precision, and no growing source arena is
assumed to fit in these dimensions.

Henceforth let $`n\ge n_0`$.

### High precision first

If $`L\gt N/n^2`$, use R69. It supplies the claimed optimal count,
Clifford count, and depth $`O(NL/b^2+L)`$ at the stated width.
Since $`L\le F(n,L)`$, it also implies the global F bound.

For all remaining branches $`L\le N/n^2`$, so

```math
nL\le\sqrt{NL}.
```

Thus the retained compiler has the required optimal count whenever it
is selected below.

### Linear precision

If $`L\le cn`$, use R66. Its additive depth is at most F because

```math
\log_2(n+2)\le\log_2(n+L+2).
```

Its simultaneous count and Clifford bounds are already within the
claimed bounds.

### Remaining precision covered by the whole-frame stream

If $`L\gt cn`$ and $`L\le2^{\kappa n}`$, use R68. Now

```math
n+L\le(1+c^{-1})L\le(1+c^{-1})(\sqrt n+L),
```

and therefore

```math
\min\{nL,(n+L)\log_2(n+L+2)\}
\le(1+c^{-1})F(n,L).
```

The multiplier is an absolute constant. R68's simultaneous count and
Clifford bounds are unchanged.

### The intervening precision range

In the remaining case $`L\gt2^{\kappa n}`$ while
$`L\le N/n^2`$. Use the retained uniform compiler. Its count
term $`nL`$ was already absorbed above, and

```math
nL\lt\kappa^{-1}L\log_2L
\le\kappa^{-1}(\sqrt n+L)\log_2(n+L+2).
```

Together with the identity bound $`nL\le nL`$, this gives
$`nL\le\kappa^{-1}F(n,L)`$. Thus its depth also satisfies
the global bound. No overlap between R68's range and R69's range is
required.

These branches exhaust every n and L. In each branch one complete
circuit supplies all three resource bounds, the full error estimate,
and the literal workspace allocation. The proof selects among circuits;
it does not combine a count bound from one circuit with a depth bound
from another. This proves the theorem, including the stronger
high-precision clause. ∎

## 4. Count optimality and matching depth intervals

At the supplied width the physical number of qubits is
$`n+b+2=\Theta(b)`$. The inherited lower bounds therefore give

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad D_T^\star=\Omega(NL/b^2).
```

The upper bound proves the full worst-case count frontier

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right)
```

at every eligible width. For depth, every nonempty interval

```math
17(L+n+7)\le b\le\sqrt{\frac{NL}{F(n,L)}}
```

has simultaneous matching

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

Indeed, the upper endpoint gives $`F(n,L)\le NL/b^2`$.
Since $`F(n,L)\ge L\ge6`$, it also gives
$`b\le\sqrt N\le N`$ and $`b\le\sqrt{NL}`$.
Thus both L and $`\sqrt{NL}`$ are bounded by $`NL/b`$.
For integer n and L, the interval is nonempty exactly when

```math
289(L+n+7)^2F(n,L)\le NL.
```

The stronger high-precision additive term L would instead give the
upper endpoint $`b\le\sqrt N`$. That interval is empty under
the present literal allocation when $`L\gt N/n^2`$. To see this,
$`n^2/\sqrt N\le9/2`$ for every positive integer n: its ratio
between consecutive n is below one from $`n=6`$ onward, and
the first six values have maximum $`9/2`$. Hence

```math
b\ge17(L+n+7)\gt17N/n^2\ge(34/9)\sqrt N.
```

The high-precision upper bound is nevertheless $`O(L)`$ at every
eligible width, since

```math
\frac{N}{b^2}\lt\frac{n^4}{289N}\le\frac{81}{1156}.
```

This simplification supplies no matching unrestricted depth lower bound.
In particular, neither F nor L is asserted necessary as an additive
depth term. Optimal large-width depth remains open.

## 5. Evidence and scope

This chapter is an analytic combination of complete-circuit theorems.
The native identities, source cleanup, literal phases, and bounded
implementation evidence remain in their linked proof homes. No new
finite fixture is needed to establish the branch inequalities or the
uniform finite-dimension fallback, and the existing source and kernel
receipts are unchanged. The result introduces no general native emitter,
practical crossover estimate, efficient preprocessing guarantee, or
unrestricted optimal-depth claim.
