# A directly prefix-programmed high-precision terminal

[Fast geometric references](FAST_CAPPED_GEOMETRIC_REFERENCE.md) ·
[Global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md) ·
[Retained high-precision compiler](HIGH_PRECISION_DEPTH.md) ·
[Exact parallel queries](PARALLEL_BILINEAR_WORD_QUERY.md)

This chapter programs
short high-precision local compilers directly from their actual global
prefixes. A small conditional program bit selects a native gate; no
catalogue of all final-accuracy words or array of all prefix slots is
allocated. The complete global precision splice is a separate proof.

## 1. Conditional group interface

A global n-bit register has immutable prefix P (t bits), active group X
(h bits), and lower suffix Z (s bits), with $`n=t+h+s`$ and $`s\ge3`$. Let the
sole external initialized qubit H save $`[Z=0^s]`$. Choose three distinct
original Z bits d,z,w. On $`H=1`$ they are zero; on $`H=0`$ all may be occupied.
The group target is the direct sum of $`W_P`$ on X when $`Z=0`$
and identity otherwise.
Use d as the local one-clean R69/R72 flag, z as a reusable native-gate
program selector, and w only when a finite optional subroutine must be
present for some P and absent for others. All other local compiler wires
are arbitrary dirty, disjoint from P,H,z,w and the rest of Z.

Let $`\kappa\ge6`$ be local accuracy; the local large-height branch assumes
$`h\ge h_0`$ (the retained absolute R69 cutoff) and $`\kappa\gt2^h/h^2`$. The intended
choice is $`h\le\lfloor\log_2(\kappa+2)/32\rfloor`$, with $`\kappa=L+O(1)`$. The number of
groups will be constant. The direct-prefix terminal is needed only on
$`L\ge n^{1/4}`$, $`L\le N/n^2`$; other ranges have separate complete compilers.

The new quantum implementation must be EXACTLY the direct sum of chosen
actual local native words on $`H=1,z=w=0`$, including their dirty-core and
rejected-sector action. It must be literal identity on all $`H=0`$ inputs.
No extra source or synthesized controlled-unitary approximation is used.

## 2. A common local high-precision grammar

Retain the construction and explicit local ledger of
[R69 Sections 2–6](HIGH_PRECISION_DEPTH.md#2-exact-replacement-at-a-completed-mask-boundary), with n there replaced
by h and local dimension $`M=2^h`$. Its high-level fine-query positions are:

* At every $`d\gt d_0`$: constantly many growing-coefficient masks, with
  $`Q_d=O((M/2^d)P_d)`$, $`P_d\asymp(d+1)^2`$, and
  $`m_d=O(\kappa/P_d+P_d^3\log_2(P_d+2))`$;
* cheap-row masks and their ACTUAL column replays, $`Q_d=O(M/2^d)`$,
  $`m_d=O(\kappa/P_d+\log_2(P_d+2))`$;
* the bounded metadata/type masks, no larger row order;
* constantly many full fine/source/attenuation/tag/replay/bottom banks,
  $`Q=O(M)`$, $`m=O(\kappa)`$.

Refinement jumps are at the dimension-dependent positions where P_d
changes; the four binary substeps and their precision widths depend on
h and $`\kappa`$, not on the input data. Fix the padded precision/core
layout and cheap-scalar grammar at each height BEFORE computing any mixed
gauge. In the [growing-source construction](../endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count),
$`P_d,M_d,q_d,\Delta_d`$ and $`q_{\mathrm{ch},d}`$ are functions
only of height and requested precision. Every coefficient/type word is padded
to its dimension-dependent output length. Exponents remain offline mask
data, not quantum addresses. Cheap/expensive/zero nodes retain the same
allocated labels and their known complete guards. All invalid/padding
rows have their retained total definition.

The common skeleton is indexed by height, call type and one of a fixed
number of fine-source banks. The data-dependent permutations, mode
packing, finite coarse synthesis and chart-center words occur BETWEEN
these fine positions and belong to the coarse remainder below. They
need not have a common native gate list. The local mode layouts and
precision core bases are fixed before selecting P; the known packing
words map each P-dependent support arrangement into the corresponding
fixed physical fine-bank layout.

In particular, fix the same worst-case full arenas from
[GROWING_SOURCE_COMPILER Section 4](../endpoint/GROWING_SOURCE_COMPILER.md#4-full-unitary-mode-arenas-and-exact-guards)
for every P. Their complete dyadic reservations count all node/port,
source, rare/power, multiplier, flag and arithmetic values before data
are supplied; their aligned placement depends only on local height.
Zero, cheap and expensive nodes keep their allocated feature labels.
The node and power address coordinates of each coefficient mask can
therefore be fixed at that height. Prefix selection does not replace
these compressed addresses by the whole mode register. Data-dependent
live-cell transports remain in the coarse permutation words. This
preserves both the geometric row sum and its square-root sum below.

When a fixed bounded number of disjoint-block bins is not fully used,
pad its list with actual empty subroutines. To emit such a conditional
subroutine, query its P-dependent presence bit into w while $`H=1`$;
use $`H\wedge w`$ as its complete external guard, then reverse that presence
query. An absent bin is thus FULL physical identity, not a phase-zero
approximation. The coarse-gate selector z remains separate while w is
live. This adds only constant-many short P queries per original finite
interface. The compiler's fixed $`\kappa`$-dependent fine loaders themselves
are shared schedules; different P values do not select arbitrary
precision-length instruction strings.

The retained sums, INCLUDING actual inverses and recovery occurrences,
are

```math
\sum_jQ_jm_j\le CM\kappa,\qquad
\sum_j\sqrt{Q_jm_j}\le C\sqrt{M\kappa},
```

```math
\sum_jm_j\le C[\kappa+P(h)],\qquad
\#\{\mathrm{fine\ positions}\}\le P(h).
```

These use one fixed polynomial P and are the sums from
[R69 Section 5](HIGH_PRECISION_DEPTH.md#5-uniform-sums-of-the-rescheduled-queries),
not an inference from the final complete-compiler big-O statement.

## 3. Exact prefix selection of complete fine masks

For each completed fine XOR mask, append P and H (and w only while an
optional finite bank is selected) to its actual local address. The
inactive rows are zero. The local address itself and output/core retain
their original meaning. R91 then implements exactly the guarded direct
sum of the original completed masks. Z masks retain their output-H
conjugation and X/Z order; those elementary H banks are charged below.
The query's helpers are in one exclusive arbitrary pool B_Q and return
exactly. Neither the conjugate signal r nor any live core or program
wire is in B_Q.

The global row count is at most $`C2^tQ_j`$. Thus, writing $`N_s=2^{t+h}`$,

```math
\sum_jQ_j^{\mathrm{glob}}m_j\le CN_s\kappa,\qquad
\sum_j\sqrt{Q_j^{\mathrm{glob}}m_j}\le C\sqrt{N_s\kappa}.
```

Sequential R91 implementations have

```math
T_{\mathrm{fine}}\le C\left[\sqrt{N_s\kappa}+N_s\kappa/B_Q\right],
\qquad G_{\mathrm{fine}}\le CN_s\kappa,
```

```math
D_{T,\mathrm{fine}}\le
 C\left[N_s\kappa/B_Q^2+P(h)\log_2(n+2)\right].
```

This changes the complete operator only by the intended P/H selection,
not the actual scalar rejected operator. The actual row word/core is
fixed before computing its mixed gauge, and its column uses its actual
native replay/inverse. No $`2^t`$ factor enters approximation error.

## 4. Elementary length of precision-sized native words

All remaining precision-sized loader words have $`O(m_j)`$ ELEMENTARY
Clifford+T length, not only T-count. The Majorana split construction in
[the operator-source construction](../../docs/OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit) uses
$`2(m-1)`$ T and $`O(m)`$ Clifford gates; its two-tail construction
has $`O(q)`$ T and $`O(q)`$
Clifford gates. The one-clean source uses the same weight-at-most-two
split generators. There are constantly many such loader occurrences
per completed fine scalar call. Output H banks add $`O(m_j)`$.

The exact dyadic refinement source widths halve strictly, so their sum
is $`O(\kappa)`$; see the
[refinement sum](../endpoint/GROWING_SOURCE_COMPILER.md#7-summable-native-charges)
and its [precision replacement](../endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count).
Its rare-event/index arithmetic, multiplier fields and fixed quarter-turn
SELECTs operate on polynomial-in-h mode registers. Large exponents stay
offline. Their elementary gates, metadata, central predicates and all
other finite routing belong to P(h) or the coarse remainder. The
[conjugate lift](../endpoint/CONJUGATE_SIGNAL_RECOVERY.md#2-lift-the-whole-high-level-word) adds only a
constant number of native gates per original gate.

Consequently these common precision-sized elementary lists have total
length $`O(\kappa+P(h))`$. Apply the exact bounded-control native-gate
[construction for conditional suffixes](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#exact-conditioning-without-a-clifford-count-penalty)
to every such gate, with external guard H or $`H\wedge w`$. Its controlled-T construction
uses RETURNED DIRTY helpers, not an initialized AND target. It has
constant exact cost for the fixed number of added controls. Thus this
entire guarded portion has T, G and serial T-depth $`O(\kappa+P(h))`$.

The original local multi-literal arena guards remain part of the local
word; their original native gates are already counted. Adding the one
or two external literals does not replace their h-bit predicate by a
new n-bit predicate. P is used only in the paid queries.

## 5. Exact coarse programming with short selection queries

[R69 Section 6](HIGH_PRECISION_DEPTH.md#6-every-remaining-native-word-is-affordable-in-depth) independently bounds its coarse/routing/state-source
remainder by $`T\le CM^{3/4}P_*(h)`$, $`G\le CMP_*(h)`$. Therefore its TOTAL
native elementary length is at most $`S=CMP_1(h)`$. This remainder is
precision-independent in length and wire support even though its data
are evaluated from the requested frame and actual fine-word gauge.

Choose a common coarse wire set of size $`W_c\le CMP_2(h)`$, containing
all genuinely acted-on local mode/flag/signal/metadata wires and a
separate dimension-dependent coarse helper pool. Do not put the full
$`\kappa`$-bit precision core in this set. Original coarse helper layouts
may be renamed into the common coarse helper pool because their
complete exact logical words return those arbitrary helpers. Genuine
shared metadata and signal targets keep their prescribed locations.
Actual cheap-scalar cores and their replay remain in the fine class.
The one shared conjugate signal r is included at its actual fixed
location in the coarse support: lifted coarse T gates can target it.
This adds one wire, never the complete precision core.

Between successive fine positions, independently generate each P's
actual coarse word. Pad it to the uniform bound S; multiplying this
bound by the polynomial number of segments still has the form M P(h).
For every position, every elementary native gate type, and every ordered
one- or two-wire target tuple in W_c, let f(P) be its finite Boolean
indicator. Fixed target tuples are distinct from P,H,z,w. At most one
indicator is one for that P and position.

Chronologically emit

```math
Q_f:\ z\mathrel{\oplus}=Hf(P),\qquad
C_{H\wedge z}(G),\qquad Q_f^\dagger.
```

While an optional bank is present, replace the query enable by
$`H\wedge w`$ and the central gate control by $`H\wedge w\wedge z`$.

Q_f is an exact R91 query into the ORIGINAL conditional-zero z. On $`H=1`$
and $`z=0`$ it gives exactly $`G^{f(P)}`$ and returns z to zero. On $`H=0`$ the
query has zero rows and the central controlled gate is identity even
when z,d,w are occupied. The bounded-control native-gate lemma covers
H,S,T,T-dagger,CNOT with constant returned dirty helpers and literal
phases. No program catalogue or array of P slots is allocated.

Fine words never touch z. Thus z returns exactly after each use and
remains zero on the active sector even when prior source error occupies
d or disturbs arbitrary dirty cores. If w is live it is preserved until
its actual P-query inverse. Earlier errors propagate through complete
actual unitaries, not through an assumed clean reset.

The universal schedule has at most

```math
J\le C2^{3h}P_3(h)
```

serial Boolean-query/controlled-native positions: S times $`O(W_c^2)`$.
Every query has at most $`C2^t`$ rows, with only a fixed number of extra
H/w controls. Its output is one bit. Therefore

```math
T_{\mathrm{coarse}}\le CJ[2^{t/2}+2^t/B_Q+1],\qquad
G_{\mathrm{coarse}}\le CJ(2^t+1),
```

```math
D_{T,\mathrm{coarse}}\le
 CJ[\log_2(n+2)+2^t/B_Q^2+1].
```

Equivalently the exponential factors are

```math
J2^{t/2}\le C\sqrt{N_s}\,2^{5h/2}P_3(h),\qquad
J2^t\le CN_s2^{2h}P_3(h).
```

These are deliberately loose but explicit. P_3 and all constants are
fixed before choosing any global cutoff. A local Clifford-count bound
is essential here; a bound only on local T-count would not justify S.

## 6. Width at the original literal global threshold

Use the inherited local high-precision base $`\kappa+h+7`$, after stripping
its old fine-query pool. It includes the actual live precision core,
borrowed conjugate signal and original local guard fields. Reserve the
separate coarse-only support/helpers of at most $`W_c\le2^hP_4(h)`$, plus
a fixed number $`c_g`$ of returned dirty bounded-control helpers.
For $`h\le\log_2(\kappa+2)/32`$, beyond one absolute local cutoff,

```math
W_c\le\kappa.
```

Thus a sufficient nonquery base is at most $`2\kappa+h+7+c_g`$. At
$`\kappa=L+c_{\mathrm{err}}`$ and global $`b\ge17(L+n+7)`$, with $`h\le n`$ and a sufficiently
large absolute global cutoff, this base is below $`b/4`$. Reserve
$`B_Q=\lfloor b/2\rfloor`$ exclusively for R91 calls; the remaining quarter is slack.
The three original suffix bits d,z,w are not counted as dirty helpers.
None is an external initialized qubit.

The sufficient R91 minimum $`B_Q\ge C(n+2)`$ follows on the wide branch
$`b\ge\sqrt N/n`$ after an absolute n cutoff, independently of L. On large-L
minimal-width cases $`b\ge17L`$ only improves it. The fine output core and
all coarse targets are disjoint from B_Q. Coarse sources and fine
sources are sequential; the stated bound nevertheless reserves both.

If h is below the absolute R69 dimension cutoff, use the original
bounded-height one-clean native addressed-layer words. Their number
of layers is an absolute constant. Appending P to their masks gives
the same global volume orders; their common fixed loaders have
$`O(\kappa)`$ length. Their program-dependent coarse remainder is bounded,
and the same H guard/d initialization argument applies. No large-h
R69 asymptotic is invoked on a bounded-height group.

## 7. Error and saved-predicate return

The group begins with the exact saved predicate $`H=[Z=0^s]`$. A generic
R91 XOR query for this one Boolean function has at most $`2^s\le N`$ rows
and a separate H output, so a constant number of groups contributes
only $`O(\sqrt N+N/b)`$ T, $`O(N)`$ Clifford, and $`O(\log_2(n+2)+N/b^2)`$ depth.
This conservative price suffices; no unproved shallow MCX is required.

On $`H=1,z=w=0`$, the newly emitted word is the direct sum of actual local
words at precision $`\kappa`$. On $`H=0`$ it is exact identity on all occupied
work. Since the target local compiler has initialized d error at most
$`2^{-\kappa}`$, the group has the same full initialized-isometry error on
all X, P, dirty and reference inputs. All original signs, zeros,
source rejection, fine core return and actual inverse errors are
those of the retained local word. The global prefix sectors take a
maximum. No output scalar phase is discarded.

The actual saved-predicate inverse is applied last. The ideal group
returns d,z,w and the rest of Z, hence preserves the original suffix
predicate and returns H. Its actual leakage is transported unitarily
through that inverse. A constant number R of groups adds at most
$`R2^{-\kappa}`$ plus the separate bounded bottom error. Fix R first, then
choose the constant $`c_{\mathrm{err}}`$.

## 8. Absorption on the intended branch

For $`h\le\log_2(\kappa+2)/32`$ and $`\kappa=L+O(1)`$, every fixed polynomial P(h)
is below $`L^\epsilon`$ eventually. For example fix $`\epsilon=1/32`$ after all
degrees. Then

```math
J\le CL^{1/8},\qquad
2^{5h/2}P(h)\le CL^{7/64},\qquad
2^{2h}P(h)\le CL^{3/32}.
```

The first exponent uses $`3/32+1/32=1/8`$. Therefore coarse counts are bounded by
$`O(\sqrt{N_sL}+L+N_sL/b)`$, G by $`O(N_sL)`$, and depth by

```math
O\!\left(L+N_sL/b^2+L^{1/8}\log_2(n+2)\right).
```

On $`L\ge n^{1/4}`$, the last term is $`O(L)`$ after one absolute cutoff:
$`\log_2(n+2)=O(L^\epsilon)`$ for every fixed $`\epsilon\gt0`$ there. Fine query
$`P(h)\log n`$ and all common loaders satisfy the same bound. A bounded
number of groups therefore gives $`O(L+\log n+NL/b^2)`$ depth and the
retained count bounds. The terminal may be used whenever the remaining
height is $`O(\log L)`$, with an absolute constant; for $`L\ge2^{cn}`$, the
whole n-bit input already has this property. In both cases keep three
bottom bits for the retained bottom compiler.

This component supplies the terminal interface for
[the all-precision theorem R96](ALL_PRECISION_LINEAR_DEPTH.md). Its native
construction and resource bounds are independent of that global assembly;
it does not assert the full precision splice by itself.
