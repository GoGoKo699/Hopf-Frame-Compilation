# Cached activity predicates across a window of early groups

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Protected unary source](PROTECTED_UNARY_SOURCE.md) · [Unary group interface](../../docs/UNARY_PHASE_GRADIENT.md) · [Uniform precision schedule](../../docs/UNIFORM_PRECISION_DEPTH.md) · [Current frontier](../../docs/OPEN_PROBLEM.md)

A window of consecutive early groups can share one predicate for its
unchanged outer suffix. Short cached block predicates and a suffix-product
chain supply each group's activity flag. The caches are consumed before
their logical controls change and returned exactly at the window boundary.

With $`J=\lceil\log_2(n+2)\rceil`$ groups per full window, this reduces
the aggregate early activity-predicate T-depth to

```math
O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}+\log(n+2)\right)
```

in the uniform low-precision range. It uses $`2J+1`$ additional protected
logical bits, the same two initialized flags, and the same two arbitrary
dirty predicate helpers. The logical group stages and program queries
retain separate $`O(n)`$ depth allowances. The complete-frame frontier
and its proved width thresholds therefore remain unchanged.

## 1. Registers and the exact window contract

Use $`n\ge1`$, $`N=2^n`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
The protected-source construction computes into an initialized flag H
whether a fixed terminal logical bank B was originally zero. Its other
initialized flag h returns to zero after every complete group. Enlarge B
by a disjoint cache slice of $`2J+1`$ bits. The unary preparation U acts
only on the source slice and leaves the cache slice untouched.

Consider a window of $`1\le t\le J`$ groups. Write its consecutive
logical target blocks as $`C_1,\ldots,C_t`$, of positive lengths
$`g_1,\ldots,g_t`$. Let Z be the logical suffix beyond the entire
window with B excluded. Reserve cache bits

```math
e_1,\ldots,e_t,\qquad u_1,\ldots,u_t,\qquad v
```

inside B; unused cache bits are untouched. On H equal to one these bits
start at zero. On H equal to zero their inputs may be arbitrary and
entangled. No cache bit belongs to the source, convolution work, program,
selected row, selector tree, or arbitrary dirty query pool.

Each complete group body $`G_i(h)`$ includes its program load, logical
stages, selector cleanup, and actual program inverse. It preserves h and
all cache bits. When h is zero it is literal identity on arbitrary source
and group-work inputs. When h is one with valid conditional-zero work,
it changes only its own logical block and the source core; all future
logical blocks and Z return exactly. This return holds for every source
core input and reference correlation, including source-preparation error.
These are the existing unary group's interfaces.

The window supplies the activity value

```math
h=H[Z=0]\prod_{j=i+1}^{t}[C_j=0]
```

around $`G_i`$. It returns h and every cache bit exactly. In the active
sector its completed action is the same ordered group word as individually
computing each outer predicate. In the sector $`Hh=00`$ its whole action
is identity on arbitrary cache, source, logical, dirty, and reference
inputs. The latter statement requires h initially zero.

## 2. Native cache schedule and consume-before-change cleanup

Use the established exact two-dirty-helper multi-control X for zero
predicates. For a length-a block it has $`O(a+1)`$ T and Clifford count
and $`O(\log(a+2))`$ T-depth, preserving its controls and returning
both arbitrary helpers. An empty block gives a Clifford X. Compute all
short predicates sequentially so that the same helper pair suffices
throughout the window.

First compute $`v\mathrel\oplus=[Z=0]`$ and then
$`e_i\mathrel\oplus=[C_i=0]`$ for every i. Compute the suffix chain
in descending order:

```math
u_t\mathrel\oplus=1,\qquad
u_i\mathrel\oplus=e_{i+1}u_{i+1}\quad(i=t-1,\ldots,1).
```

Every product toggle is a literal exact Toffoli. On zero cache input,
$`u_i=\prod_{j=i+1}^t[C_j=0]`$, including $`u_t=1`$.
For each i in ascending order:

1. Apply the actual inverse of the predicate that computed $`e_i`$.
   Its block $`C_i`$ has not yet been a logical target.
2. Toggle $`h\mathrel\oplus=Hvu_i`$, apply the complete body
   $`G_i(h)`$, and reverse that exact h-toggle word.
3. If $`i\lt t`$, undo the Toffoli with controls
   $`e_{i+1},u_{i+1}`$ and target $`u_i`$. For $`i=t`$, undo
   the X on $`u_t`$.

Finally apply the actual inverse of the Z-predicate on v. Its logical
controls have returned after every complete group.

The three-control h-toggle has an exact constant-size implementation
with one returned arbitrary dirty bit d from the idle helper pair.
Chronologically, apply

```math
\operatorname{CCX}(H,v;d),\quad
\operatorname{CCX}(d,u_i;h),\quad
\operatorname{CCX}(H,v;d),\quad
\operatorname{CCX}(d,u_i;h).
```

The two h toggles sum to $`u_i(d\oplus Hv)\oplus u_id=Hvu_i`$,
and d returns. Four exact seven-T Toffolis give constant T and Clifford
count and constant T-depth on arbitrary inputs. Use the actual reversed
native word for the second h-toggle. No relative-phase replacement or
additional initialized wire is used.

### Active cache correctness

On H equal to one the cache starts at zero. Immediately before group i,
its future block predicates and suffix-chain controls still describe
their original blocks. Every earlier completed group has restored these
logical blocks, even if it temporarily borrowed some as conditional work.
The cache value $`u_i`$ is therefore exactly the required suffix test.

Erasing $`e_i`$ before its group preserves the information needed for
activity because $`u_i`$ involves only later blocks. After the group,
the next block and its caches are still valid, so the chain gate erases
$`u_i`$ exactly. Nothing attempts to reconstruct an old predicate from
an already changed target. The active body's conditional work is valid:
$`v=u_i=1`$ certifies the entire nonbank outer suffix is zero. Source
and convolution work retain their separate protected-bank guarantees.
This proves the active contract by induction without resetting or
projecting the source.

### Arbitrary inactive-cache return

On $`Hh=00`$, the exact three-control toggle leaves h zero regardless
of its dirty helper or cache inputs. Every body is then identity. For a
computational-basis input, write the initial block cache as $`a_i`$,
the initial chain cache as $`w_i`$, and $`f_i=[C_i=0]`$. After
preparation their values are

```math
b_i=a_i\oplus f_i,\qquad
c_t=w_t\oplus1,\qquad
c_i=w_i\oplus b_{i+1}c_{i+1}.
```

Ascending cleanup first restores $`a_i`$ using the unchanged logical
block. It then restores $`w_i`$ using the still-intact future values
$`b_{i+1},c_{i+1}`$, or the matching X at the last block. The final
v inverse uses its unchanged logical controls. Thus every arbitrary
initial cache value and every dirty helper returns, with no phase.
Linearity extends the identity to coherent cache inputs and arbitrary
reference entanglement. The proof also covers $`t=1`$ and a partial
final window.

## 3. Simultaneous reservations and the revised cutoff

Let $`q=2^\ell`$ be the unary phase modulus and
$`\rho=\log_2 3`$. Choose the same absolute reservation constant
$`A\ge1`$ as in the protected-source construction. Define

```math
s_{\rm source}=\left\lceil A(q^\rho+q\ell+\ell+1)\right\rceil,
\qquad S=s_{\rm source}+2J+1,
\qquad K=16S.
```

The fixed terminal bank B has S logical bits. At remaining height
$`k\gt K`$ use

```math
g=\left\lfloor\log_2\frac{k}{16Aq}\right\rfloor.
```

The original estimates still give $`g\le\ell`$ and
$`g\ge\lfloor(\rho-1)\ell\rfloor`$. In particular,

```math
S\le k/16,\qquad Aq2^g\le k/16,\qquad
Ag\le s_{\rm source}\le k/16.
```

The bank plus dynamic group work $`A(q2^g+g)`$ uses at most
$`3k/16`$ logical bits. The current outer suffix has at least
$`15k/16`$ bits, so the dynamic program and selector work fit outside
B at every group. Partition the resulting early-group list into windows
of at most J groups. This does not require an additional simultaneous
bank of J programs: each complete group returns its program and dynamic
work before the next group.

The last early group leaves height $`K'\gt K-\ell\ge15S`$.
Every window target therefore remains outside B, including in a partial
last window. The actual source inverse and H erasure precede the tail.
When the early segment is empty, retain the established compiler.

Exactly two arbitrary dirty helpers suffice for all cache predicates,
the original B-zero predicate, and the constant-arity h-toggles. Their
uses are sequential and each returns before its next use. They are
disjoint from the logical bank, every live program, and the query pool.
With $`B_0=L+n+7`$ and tail precision $`L'=L+2`$, the same conservative
extra query pool remains

```math
B_{\rm extra}=b-(B_0+2)-2=b-B_0-4\ge b/2
\qquad\text{when }b\ge17B_0.
```

Use this refinement only where its bank and query reservations fit; the
existing finite-size and small-width fallbacks keep the literal
complete-frame threshold. No conditional logical cache bit is also counted as an
arbitrary dirty helper.

## 4. Predicate depth, counts, and the complete-frame scope

Let E be the number of early groups and $`W=\lceil E/J\rceil`$ the
number of nonempty windows. A window with t groups pays one long
Z-predicate pair, all short block-predicate pairs, the chain, and two
constant-arity h-toggles per group. Its activity cost is

```math
T,G=O\!\left(|Z|+1+\sum_{i=1}^t(g_i+1)\right),
\qquad
D_T=O\!\left(\log(n+2)+\sum_{i=1}^t\log(g_i+2)+t\right).
```

Both short-predicate preparation and consumed cleanup are charged
sequentially. Since $`g_i\le\ell`$, the aggregate depth is

```math
D_{T,\rm activity}
=O\!\left(W\log(n+2)+E\log(\ell+2)\right).
```

Throughout $`6\le L\le\log_2(n+2)/16`$, choose the least power of
two $`q\ge4\pi\sqrt{2n}/\eta`$. Then $`\ell=\Theta(\log(n+2))`$,
$`E=O(n/\log(n+2))`$, and $`W\le E/J+1`$, with absolute constants
for sufficiently large n. Consequently

```math
D_{T,\rm activity}
=O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}+\log(n+2)\right),
\qquad
T_{\rm activity},G_{\rm activity}
=O\!\left(\frac{n^2}{\log^2(n+2)}+n\right).
```

The count bound uses $`|Z|\le n`$ and
$`\sum_i g_i\le n`$ over the complete early segment. These polynomial
counts fit the retained $`O(\sqrt{NL})`$ T and $`O(NL)`$ Clifford
budgets. Adding $`2J+1=O(\log(n+2))`$ protected bits preserves

```math
K=O\!\left((n+2)^{9\rho/16}
 +(n+2)^{9/16}\log(n+2)\right),
\qquad K[L+\log(n+2)]=o(n).
```

The protected-source proof applies with the enlarged bank. The cache
starts and ends at zero in its active embedding and returns arbitrarily
on its inactive embedding. Hence the same global source comparison
charges only $`2\delta`$, with $`\delta=\eta/8`$; the exact window
words add no approximation error. The actual final source inverse and
H erasure still include any residual bank or flag leakage in that norm.

At fixed accuracy and sufficient square-root-scale dirty width, the
updated depth ledger is:

| Contribution | Total T-depth allowance |
|---|---:|
| Global protected-source and B-zero boundaries | $`O(\log(n+2))`$ |
| Windowed early activity predicates | $`O(n\log\log(n+2)/\log(n+2)+\log(n+2))`$ |
| Early logical stages and incremental selectors | $`O(n)`$ |
| Early program prefetch and unload | $`O(n)`$ |
| Improved late tail | $`o(n)`$ |

The two remaining linear allowances are upper-ledger costs, not lower
bounds. This local refinement gives no sublinear complete-frame theorem,
new matching interval, or high-precision endpoint claim.

## 5. Proof and evidence boundary

The arbitrary-cache identity, validity of consumed predicates, global
source-error comparison, and asymptotic reservation are analytic claims.
The five checks in
[test_windowed_group_predicates.py](../../tests/test_windowed_group_predicates.py)
audit exact active windows with one, two, and three unequal blocks;
noncommuting reduced group bodies that retain an arbitrary source and
temporarily borrow future logical work; arbitrary inactive cache and
dirty-helper inputs; the literal native four-Toffoli flag word and its
actual inverse; and negative cases with late block erasure, late chain
erasure, or a missing original-bank guard H.

The complete-window checks use exact rational reduced group operators.
Only the constant-arity flag gadget is emitted as a native Clifford+T
word. The fixture does not emit a logarithmic-depth long predicate,
loaded-program query, or scalable complete-frame compiler, and its
finite cases do not replace the symbolic proof.
