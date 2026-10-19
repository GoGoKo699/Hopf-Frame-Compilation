# Revision: keep the depth destination, relax the source count requirement

[Workspace](../../WORKSPACE.md) ·
[Full resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Preceding source-gap result](TWO_BAND_SOURCE_GAP_SESSION_RESULT.md)

The revision does not establish a faster general circuit. It identifies
an unsupported change of emphasis in the previous continuation and
replaces it with a paid constructive reduction. Adding a third source
band only leaves the hypotheses of the two-band obstruction; no actual
three-band source was supplied. The new useful fact is that polynomial
count and width overhead can be removed from a complete shallow
indicator by a root splitter and an exact short-address suffix.

## 1. Separate the destination, theorems, components and limits

| Role | Current conclusion |
| --- | --- |
| General destination | Optimal unrestricted native T-depth on the same complete-frame circuit as the count, width, precision and return guarantees |
| Strongest unconditional general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md): $`D_T=O(NL/B^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/B)`$, $`G=O(NL)`$, with one clean qubit and the original real/complex dirty reservations |
| Selected manuscript theorem | R79 remains Result D; this research revision does not change the selected statement |
| Unrestricted lower information | Indicator [Pauli nullity](INDICATOR_PAULI_NULLITY.md) gives linear T-count at constant full-operator error, not growing depth. The two-layer and zero-clean three-layer Boolean theorems retain their stated clean allocations. The general fixed-accuracy depth interval remains logarithmic versus constant |
| Paid useful components | Approximate R98, bounded-clean promotion, one-signal recovery, the new root fusion, and the branch-difference recurrence. None supplies the missing shallow source |
| Approach-specific limit | The [fixed Clifford-sum gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) excludes fixed-positive-amplitude, uniformly bounded-clean two-band sources at sufficiently small fixed error and sufficiently large address size, independent of count and width |
| Focused lower scout | Original physical pairs retain their profile/exceptional-module constraints. The new minimal-pair family has a complete zero-clean rewrite giving diagonal third-level closure and sign/XOR degree at most three/two. This is for a fixed input/middle family, not all original pairs or unrestricted depth |

Here $`N=2^n`$, $`n\ge1`$, $`L\ge6`$, and B is the complete
frame's dirty budget. R96 requires real $`B\ge17(L+n+7)`$ or
complex $`B\ge17(L+n+8)`$, and includes complete initialized-isometry
error at most $`2^{-L}`$. No allocation or operator promise is weakened.

The two-band theorem does not exclude every modest numerical tolerance,
three or more source bands, growing clean count without a uniform bound,
or unrestricted constant depth. Its finite-p dressed-family screens
retain their separate, often stronger constants. Failure of one
representation is not a lower bound on the full model.

## 2. The first unpaid operation in the previous candidate

The one-signal route asks for the complete compressed block
$`K\approx I_p/\sqrt2`$. Recovery, actual inversion and clean removal
have already been paid. For arbitrary source depth $`d_U`$, their
complete recurrence has

```math
D_T\le24d_U+32,\qquad
\epsilon_{\mathrm{ind}}
\le16(3\sqrt2\epsilon+2\epsilon^2).
```

The value 104 is just the specialization $`d_U=3`$. It is not an
implemented indicator. The first unpaid operation is the entire
logical transition, including the sparse address response across all
arbitrary output inputs. Growing probe rank or Clifford complexity
does not establish that response. Exactly three bands is stronger than
the actual sufficient goal of sublogarithmic depth.

## 3. Chosen change of mechanism and the proved recurrence

The [root-fusion theorem](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
allows a complete b-bit dirty source E with

```math
T_E,w_E=2^b\,\mathrm{poly}_\delta(b),\qquad
G_E=4^b\,\mathrm{poly}_\delta(b).
```

Choose $`b=p-r`$, where $`r=O_\tau(\log(p+2))`$, and use E
only for the first root split. Its controlled use is a two-source
bilinear echo, not a free controlled-E gate. Finish the r-bit suffix
with the exact retained tree. Write $`q=2^b`$, $`m=2^p`$ and
H for that exact suffix traversal. The complete recurrence is

```math
\begin{aligned}
T&\le8T_E+48q+4T_H,\\
G&\le8G_E+96q+4G_H+2m+32,\\
D_T&\le8d_E(b)+32+4D_T(H),\\
w&\le2m+1+\max\{w_E,w_H\},\\
\epsilon_{\mathrm{ind}}&\le8\delta.
\end{aligned}
```

Setting $`\delta=\tau/8`$ gives the original linear T/dirty,
quadratic Clifford target, with depth

```math
8d_{\tau/8}\bigl(p-O_\tau(\log(p+2))\bigr)
+O_\tau(\log\log(p+4)).
```

All source/inverse appearances, both old output banks when promotion
is used, live tree nodes, private pools, literal phases, and return
error are charged. The norm is on all address, output and dirty inputs
and their references. Only the first large source is approximate, so
its requested complete-bank accuracy stays constant with p.

The [differential refinement](DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md)
separately proves

```math
e_i\le e_{i-1}+\delta_i,\qquad
\epsilon_{\mathrm{ind}}\le2\sum_i\delta_i.
```

Only the selected parent's two-branch relative action enters. Common
right-side private action and errors at untouched parents cancel through
the actual inverse. Independent child errors under the selected parent
can still add. A simultaneous native schedule is required in addition
to this operator factorization.

## 4. Does the new missing operation change the full frontier?

Yes, if a complete source satisfying these relaxed budgets actually
has $`d_E(b)=o(\log(b+2))`$ at each fixed requested accuracy.
The paid R98 substitution preserves fixed-L count and literal widths,
with schematic full depth

```math
D_T=O_L\!\left(
N/B^2+\overline d_{\delta_L}(2n+10)+\log_2\log_2(n+4)
\right),
\qquad \delta_L=\frac{2^{-L}}{128K_*}.
```

The depth envelope and error allocation are detailed in the fusion
proof. Constant source depth would give a log-log additive frame term;
it would not establish optimality or a new all-precision L dependence.

The reduction alone changes no unconditional frontier. Substituting
the existing logarithmic-depth indicator keeps the old logarithmic
cost. The retained finite pulse word has excessive exponential count
overhead and depth, so it also does not qualify. Polynomial independent
row count is no longer by itself a terminal objection to a *smaller*
source, but a joint full-bank accuracy proof and faster native depth
remain necessary.

A precisely priced slowly growing clean bank is now another possible
representation. The fusion proof gives a conditional window with
$`k=\lfloor\beta\log_2\log_2(b+4)\rfloor`$ and source depth
$`O((\log(b+2))^\gamma)`$, where $`\beta+\gamma\lt1`$.
Generic approximate return requires source precision
$`\delta/(8\cdot2^k)`$; all source budgets must hold uniformly
along that varying precision schedule. No source meeting this window
has been constructed.

A [first test of that window](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#8-a-first-test-of-the-growing-clean-window)
examines weight-only phase counters. With q sequential phase accesses
and at most $`2^k`$ distinct frequencies in each whole access, a
pigeonhole collision forces a nonzero weight to imitate weight zero
when $`q2^k=O((\log b)^\sigma)`$ for fixed $`\sigma\lt1`$.
This stops that compressed-counter strategy. The spectral premise is
additional to the clean budget, and q is not automatically T-depth;
more general dirty-spectrum or address/bank interactions remain outside
the obstruction.

## 5. Keep the original physical-pair scout focused

The [actual minimal-pair source](MINIMAL_PHYSICAL_PAIR_COMPLETION.md)
has one original physical pair and exact occupied identity, with its
unused clean-sector action retained. The pair yields
$`T_B T_{-Z_aB}`$. Cancelling the first factor leaves its required
central partner. The native axis independence condition prevents
adding both cancellations to the same middle band.

For this fixed input and middle word, a spectator native band W that
moves B makes the required completion non-C3. Actual rank-two clean-image
matching shows every one-band output completion has C3 logical action,
including arbitrary endpoint Cliffords and unused-sector freedom. Thus
the fixed mechanism cannot produce occupied identity or another Clifford
target with that growing orbit. This is an actual physical-provenance
limit, stronger than rejecting just one displayed output word.

The same clean comparison also gives every complete target in this
family a literal zero-clean three-band representation with one middle
axis. The retained native-remainder theorem then makes every diagonal
target diagonal third-level: Boolean sign degree is at most three and
XOR predicate degree at most two. This closes non-Clifford Boolean
targets for the specified family as well.

It does not close altered source words, general fixed-pair degree,
growing original pair count or unrestricted depth.
The next lower question is whether the tied partner can be retained in
simultaneous native matching for those remaining actual sources. The
already closed no-original-pair branch is not reopened.

## 6. Verification boundary and next question

The new source fusion and its asymptotic normalization are analytic
compositions of retained native bilinear and exact-tree primitives.
Two new branch-difference fixtures retain every physical column and
actual inverse; the minimal-pair fixture checks its actual physical
source and unused sector. These finite checks do not construct the
missing shallow source or prove an unrestricted depth bound.

The next constructive question is now explicit: **find a complete
sublogarithmic-depth indicator with polynomial address overhead, or a
native joint branch difference meeting the same paid recurrence.**
Stop a substitution when its actual native depth remains logarithmic,
its complete-bank error is unpaid, or its clean-removal multiplier
restores the old cost. Keep the unrestricted destination open while
changing that mechanism.
