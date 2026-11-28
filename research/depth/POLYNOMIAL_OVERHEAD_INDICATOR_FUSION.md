# Complete indicator fusion with one source per tree layer

<a id="polynomial-source-overhead-can-be-absorbed-by-a-paid-indicator-fusion"></a>

[Retained exact dirty tree](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Branch-difference recurrence](DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md) ·
[Native bilinear batch and clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Complete approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Two-band source obstruction](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)

In [shared-layer fusion](#shared-layer-indicator-fusion), each tree layer
applies one complete zero-clean indicator and its actual inverse,
regardless of the number of parents. A constant number of balanced layers
absorbs fixed singly exponential source counts and dirty width when the
exponent is common across fixed tolerances. A sublogarithmic source depth then improves the
general fixed-accuracy frame bound through the [R98 transfer](#7-conditional-complete-frame-consequence-with-its-source-tolerance).

Sections 1–5 give the one-prefix construction. A large approximate
indicator and an exact short suffix absorb polynomial overhead even
when its exponent depends on tolerance. The completed echo has eight
source or actual-inverse occurrences. Section 6 also covers clean
sources after charging promotion.

## 1. The relaxed complete source contract

For b address bits put $`q=2^b`$. Let $`I_b`$ denote the literal
complete indicator

```math
I_b|u,T\rangle=|u,T\mathbin\oplus e_u\rangle,
```

where the q-bit bank T and every private work input W are arbitrary.
At each requested fixed $`0\lt\delta\lt1`$, suppose a terminating
native compiler supplies a zero-clean unitary $`E_{b,\delta}`$ with

```math
\|E_{b,\delta}-I_b\otimes I_W\|_{\mathrm{op}}\le\delta.
```

The identity on W requires all private work to return in the ideal
comparison. Every address, bank and work column is included, as are
arbitrary external references. Neither
exact intermediate return nor exact address preservation is assumed for
the approximate emitted word. Bounds below are common bounds for that
word and its actual reversed native inverse, including the inverse
Clifford cost in the positive-S alphabet.

Assume fixed nonnegative exponents $`c_T,c_w,c_G`$ and constants
$`C_\delta`$, uniform in b, such that

```math
\begin{aligned}
T_E(b,\delta)&\le C_\delta 2^b(b+2)^{c_T},\\
w_E(b,\delta)&\le C_\delta 2^b(b+2)^{c_w},\\
G_E(b,\delta)&\le C_\delta 4^b(b+2)^{c_G},\\
D_T(E_{b,\delta})&\le d_\delta(b).
\end{aligned}
```

Here $`w_E`$ excludes u and T. The exponents may depend on the
fixed requested tolerance, provided they do not depend on b. No
precision-uniform bound is inferred from these fixed-tolerance bounds.

**Reduction.** For every fixed requested indicator tolerance
$`0\lt\tau\lt1`$, take $`\delta=\tau/8`$ and choose any fixed

```math
a\gt\max\{c_T,c_w,c_G/2\}.
```

For all sufficiently large p, put

```math
r=\left\lceil a\log_2(p+2)\right\rceil,
\qquad b=p-r,\qquad m=2^p.
```

There is then a native zero-clean complete indicator
$`\widetilde I_{p,\tau}`$ with common forward/actual-inverse bounds

```math
\begin{aligned}
\|\widetilde I_{p,\tau}-I_p\otimes I_{\mathrm{borrowed}}\|_{\mathrm{op}}
&\le\tau,\\
T(\widetilde I_{p,\tau}),\ w(\widetilde I_{p,\tau})&=O_\tau(m),\\
G(\widetilde I_{p,\tau})&=O_\tau(m^2),\\
D_T(\widetilde I_{p,\tau})
&\le8d_{\tau/8}(p-r)+O_\tau(\log_2\log_2(p+4)).
\end{aligned}
```

All constants are uniform in p. In particular, a source satisfying
$`d_\delta(b)=o(\log_2(b+2))`$ at every fixed requested tolerance
would supply the still-missing normalized indicator interface.

## 2. Pay for a root-controlled complete bank transition

Split the p-bit address as $`x=(u,v)`$, with b prefix bits and r
suffix bits. Borrow a root z, a q-bit old source bank T, a distinct
q-bit first tree bank A, and the source work W. Every one of these
inputs is arbitrary. The source E acts on u,T,W and acts trivially
on v,z,A.

Define the exact Hermitian native batch

```math
B_q=\prod_{j=0}^{q-1}\mathrm{CCX}(z,T_j;A_j).
```

It adds $`zT_j`$ to each $`A_j`$. Apply Hadamards to the A bank
before and after the controlled bilinear CCZ batch from
[clean promotion, Section 3](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence).
This gives common actual-inverse bounds

```math
\begin{aligned}
T_B&\le6q+(q\bmod2),\\
G_B&\le12q+3,\\
D_T(B_q)&\le4,\qquad w_B=0.
\end{aligned}
```

The shared root appears only as a Clifford control between the four
physical T bands of that batch. No helper, controlled-E operation or
initialized control is assumed. For $`b\ge1`$, q is even and the
T bound is $`6q`$.

Use the splitter, in operator order,

```math
S=E^\dagger B_q E B_q.
```

Its chronological order is $`B_q,E,B_q,E^\dagger`$. Replacing E
by the ideal $`I_b`$ gives the exact basis identity

```math
S_0:\quad
(u,v,z,T,A,W)\longmapsto
(u,v,z,T,A\mathbin\oplus z e_u,W).
```

Indeed the first batch contributes $`zT`$, the source changes T
to $`T\oplus e_u`$, the second batch contributes
$`z(T\oplus e_u)`$, and the actual ideal inverse restores T.
The two old-T contributions cancel on every input. This is a literal
phase-free operator identity, not a promise on zero-filled banks.

Two-source telescoping gives

```math
\|S-S_0\|_{\mathrm{op}}\le2\delta.
```

No error is multiplied by q: E's premise already covers its complete
q-bit bank jointly. For $`b\ge1`$, both S and its actual inverse
have the common bounds

```math
\begin{aligned}
T_S&\le2T_E+12q,\\
G_S&\le2G_E+24q+6,\\
D_T(S)&\le2d_\delta(b)+8.
\end{aligned}
```

This is where the nonlinear prefix-to-all-output action is paid: it
occurs in the two E appearances. The additional multiplication by the
arbitrary root z is paid by two B batches. Giving a scalar equality
amplitude or one initialized output column would not establish this
splitter identity.

## 3. Complete the small address suffix with exact dirty trees

Starting from each $`A_j`$ as a dirty root, run the chronological
top-down word of the [retained nonuniform tree](NONUNIFORM_DIRTY_INDICATOR.md#1-the-exact-tree-interface-permits-unequal-chunks)
on the remaining r address bits. Use the same suffix chunk schedule in
all q copies. These are top-down tree words, before their own root/leaf
echoes; the outer echo below is performed only once.

The copies share suffix address controls. Their children and private
helpers are disjoint at each stage. The retained read-only conjunction
interface places all original address and parent controls solely in
Clifford interactions, so corresponding T layers of all copies run in
parallel. Their elementary Clifford interactions remain charged. The
same arbitrary private pool can be reused between consecutive exact
stages, because each such stage returns its helpers.

Call this combined exact suffix word H. For absolute constants
$`C_T,C_G,C_w,C_d`$, it and its actual inverse satisfy

```math
\begin{aligned}
T_H&\le C_T m,\qquad G_H\le C_G m,\\
w_H&\le C_w m,\qquad
D_T(H)\le C_d\log_2(r+2).
\end{aligned}
```

Here $`w_H`$ counts the largest private conjunction pool, excluding
every live tree node. For $`r=0`$, H is identity and all four quantities
can be set to zero. The constants follow directly from the geometric
stage budget in the retained proof; q copies multiply counts and width
by q, and do not multiply T-depth.

All boundary nodes stay live. If $`0=d_0\lt d_1\lt\cdots\lt d_t=r`$
are the suffix boundaries, their total number, including z and the
first bank A, is

```math
n_{\mathrm{tree}}
=1+q\sum_{i=0}^t2^{d_i}
\le1+q(2^{r+1}-1)
=2m-q+1.
```

The q-bit old source bank T is separate from these nodes and from the
specified output Y. Its inclusion gives
$`n_{\mathrm{tree}}+q\le2m+1`$. This is why omitting T from a
generic $`2m-1`$ tree bound would be an incorrect width audit.

The source-private and suffix-private pools can share a pool of size
$`\max\{w_E,w_H\}`$, still disjoint from T and every live node.
This reuse does not require E to return that pool exactly: H is an
exact returned word on every arbitrary private input. Alternatively,
separate pools give the safe bound $`2m+1+w_E+w_H`$. We use the
shared-pool bound

```math
w_{\mathrm{total}}\le2m+1+\max\{w_E,w_H\}.
```

Define $`F=HS`$ and $`F_0=HS_0`$, so the chronological top-down
word is the splitter followed by the exact suffix. Then

```math
\|F-F_0\|_{\mathrm{op}}\le2\delta.
```

For fixed x, $`F_0`$ is an exact linear transformation of the entire
arbitrary dirty tree word. A root difference propagates to one selected
node per boundary and to precisely the final leaf labelled x. T and
the private pool return exactly under $`F_0`$. No such exact
intermediate return assertion is made about F.

## 4. Complete resource, error and return recurrence

Let C copy each of the m final tree leaves into the corresponding
arbitrary output bit of Y. Its cost is m CNOTs. Set

```math
K=F X_z F^\dagger,
\qquad K_0=F_0 X_z F_0^\dagger.
```

The chronological outer echo is $`C,K,C,K^\dagger`$. With
$`F_0`$ in place of F, the root difference changes only the selected
leaf contribution to Y; the second K occurrence restores the entire
arbitrary tree. Thus the ideal echo is exactly
$`I_p\otimes I_{\mathrm{borrowed}}`$, including literal phase.

The emitted echo has four F/actual-inverse appearances, hence eight
E/actual-inverse appearances and eight B/actual-inverse appearances.
Unitary product telescoping gives

```math
\|\widetilde I_p-I_p\otimes I_{\mathrm{borrowed}}\|_{\mathrm{op}}
\le4\|F-F_0\|_{\mathrm{op}}
\le8\delta.
```

This is a full-operator bound. It includes address and output error,
all root, old-bank, tree and private-work return error, and arbitrary
reference entanglement. Approximate source motion or leakage is never
discarded, measured or reset. Each later actual word remains unitary
on that motion. Because the ideal target is a Hermitian involution,
the same bound holds for the actual emitted inverse.

For $`b\ge1`$, the complete common ledger is

```math
\begin{aligned}
T(\widetilde I_p)&\le8T_E+48q+4T_H,\\
G(\widetilde I_p)&\le8G_E+96q+4G_H+2m+32,\\
D_T(\widetilde I_p)&\le8d_\delta(b)+32+4D_T(H),\\
w(\widetilde I_p)&\le2m+1+\max\{w_E,w_H\},\\
\mathrm{error}(\widetilde I_p)&\le8\delta.
\end{aligned}
```

The $`2m+32`$ Clifford term includes two leaf-readout words, eight
B-batch additive allowances of three, and two root X words. Each X
is the four-Clifford native word HSSH. The actual inverse reverses all
native words; the common bounds never treat an inverse S gate as a
free single positive-S gate. No output bank or live helper is omitted
from the width formula, whose only exclusions are the specified
p-bit address and m-bit Y output.

## 5. Absorb polynomial overhead without changing source accuracy

Choose $`\delta=\tau/8`$ and a as in Section 1. For all sufficiently
large p the integer $`r=\lceil a\log_2(p+2)\rceil`$ obeys
$`1\le r\le p-1`$, so $`b=p-r\ge1`$. Then

```math
\begin{aligned}
2^b(b+2)^{c_T}&\le m(p+2)^{c_T-a}\le m,\\
2^b(b+2)^{c_w}&\le m(p+2)^{c_w-a}\le m,\\
4^b(b+2)^{c_G}&\le m^2(p+2)^{c_G-2a}\le m^2.
\end{aligned}
```

All source count and width contributions therefore fit the claimed
normalized budgets. The exact suffix uses $`r=O_\tau(\log(p+2))`$
address bits and contributes only
$`O_\tau(\log\log(p+4))`$ T-depth. This proves the reduction.
The non-strict choice $`a\ge\max\{c_T,c_w,c_G/2\}`$ also
suffices when a is positive; the strict choice avoids a degenerate
zero-suffix convention.

For every smaller p where the proposed first chunk is absent, use the
retained exact indicator directly, including its Clifford zero-address
case. At each fixed tolerance there are only finitely many such p;
their counts and depth are absorbed by the tolerance-dependent
constants. One may equivalently define
$`r(p)=\min\{p,\lceil a\log_2(p+2)\rceil\}`$ and select this
fallback whenever $`b=p-r(p)=0`$. There is no call to an absent
positive-address source in that branch.

No monotonicity assumption on $`d_\delta`$ is needed for an
individual p. For a subsequent construction using many smaller
addresses, define

```math
\overline d_\delta(B)
=\max\bigl(\{d_\delta(b):1\le b\le B\}\cup\{0\}\bigr).
```

If $`b_{\max}`$ is the largest first-chunk size actually used,
replace the source-depth term by $`8\overline d_\delta(b_{\max})`$.
When all calls take the exact fallback, $`b_{\max}=0`$ and this
source contribution is zero. If $`d_\delta(b)=o(\log(b+2))`$,
then $`\overline d_\delta(B)=o(\log(B+2))`$ as well: split its
maximum at a fixed threshold where the little-o bound holds; the
finite initial maximum divided by $`\log(B+2)`$ tends to zero.

Only the first large bank source is approximate. Its error tolerance
is fixed independently of p. The q exact suffix copies create no
q-fold approximation error and require no exponentially small
per-copy source tolerance. In contrast, proposing q approximate rows
as the source E still requires a proof of E's joint full-bank error.

The short exact suffix limits the one-prefix construction. A source
with T-count proportional to $`4^b`$ already needs $`r=\Theta(p)`$ to fit
$`O(2^p)`$ gates, bringing the retained exact suffix back to a
logarithmic depth bound. [Section 9](#shared-layer-indicator-fusion)
uses approximate sources on balanced blocks to absorb this larger count.
In both constructions a logarithmic source depth still gives a
logarithmic bound.

## 6. Optional clean-source windows and the remaining operation

A complete clean-source query V can first be converted to E using
[the paid clean-promotion recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence).
For clean count k, source error epsilon, and $`M=2^k`$, the bounds
needed here are

```math
\begin{aligned}
\delta_E&\le8M\epsilon,\\
T_E&\le M(4T_V+4T_C+24q+4),\\
G_E&\le M(4G_V+4G_C+40q+16)+2q,\\
d_E&\le M(4d_V+4d_C+16),\\
w_E&\le k+q+1+w_V+w_C,
\end{aligned}
```

where the returned zero-test has
$`T_C,G_C,w_C=O((k+3)^3)`$ and $`d_C=O(\log(k+3))`$.
For bounded k, setting $`\epsilon=\delta/(8M)`$ preserves fixed
polynomial overhead and changes source depth only by a constant
factor. Combining this with the fusion gives total error at most
$`64M\epsilon`$. The old output bank and clean-register cover
work are already included in $`w_E`$ before the fusion's separate
old E-output bank is added.

The following is the serial conditional growing-clean window. The
[paid parallel construction](PARALLEL_CLEAN_PROMOTION.md#51-paying-any-sufficiently-small-packing-permutation)
gives a wider sufficient window under its native address-sharing
hypothesis, and its Section 5.2 includes this fusion's eight-call error
and resource substitution. Neither construction supplies the source.
For constants
$`\beta\ge0`$ and $`\gamma\ge0`$ with
$`\beta+\gamma\lt1`$, suppose

```math
k(b)=\left\lfloor\beta\log_2\log_2(b+4)\right\rfloor,
\qquad
\epsilon_b=\frac{\delta}{8\cdot2^{k(b)}}.
```

Assume the source has polynomial address overhead in T-count and
width, at most quadratic-in-q times polynomial address overhead in
Clifford count, and depth $`O((\log(b+2))^\gamma)`$ **uniformly
along this precision schedule**. Fixed-epsilon bounds alone would
not justify this assumption. Promotion then has polynomial address
overhead and

```math
d_E(b)
=O\!\left(
 [\log(b+2)]^{\beta+\gamma}
 +[\log(b+2)]^\beta\log\log\log(b+16)
\right)
=o(\log(b+2)).
```

Fusion would consequently produce the desired normalized sublogarithmic
indicator. For k equal to zero use V directly. If the actual V preserves
its entire clean sector exactly, the stronger promotion error is
$`2\epsilon`$ rather than $`8M\epsilon`$; one can request
$`\epsilon=\delta/2`$. Exact sector preservation does not remove
the M sequential sectors from the emitted depth recurrence.

These are sufficient windows, not source constructions. They concern
complete initialized-isometry query promises. No recovery of a general
k-clean compressed amplitude block is supplied for free. In particular,
the fixed-clean, fixed-positive-amplitude two-band obstruction remains
in force even when polynomial overhead is permitted: its gap is
independent of source counts and width. Unbounded k leaves that
theorem's stated scope, but pays the quantitative promotion factors
above. No unrestricted depth lower bound follows in either case.

The [older frame hybrid](../../docs/PARALLEL_DIRTY_LOOKUP.md#hybrid-with-the-layerwise-frame-compiler)
already absorbed polynomial indicator overhead through an early/late
frame cutoff. The retained nonuniform tree already removed polynomial
overhead from exact read-only conjunction edges. Neither establishes
the present standalone approximate whole-bank reduction: its new
interface consists of the paid root splitter, a fixed-accuracy large
bank, and an exact small-address tail.

For the one-prefix construction the missing operation is the complete
sublogarithmic source in Section 1, possibly obtained from a clean
source within the windows above. An implementation assembled from
approximate rows needs a joint operator bound on the arbitrary output
bank and every borrowed input. The
[branch-difference recurrence](DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md#2-error-and-complete-return)
offers another source contract when internal errors cancel in selected
parent differences.

## 7. Conditional complete-frame consequence with its source tolerance

Apply the [approximate R98 transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md#3-the-bounded-r98-call-count-and-final-error)
to the one-prefix indicator of Sections 1–5. A completed word query
uses at most four complete indicator occurrences, including actual
inverses, so its added error is at most $`32\delta`$. The retained
complete frame has at most $`K_*`$ such word-query occurrences,
uniformly in dimension and dirty allocation. At fixed integer
$`L\ge6`$, put $`\eta=2^{-L}`$ and request the same relaxed
source tolerance at every length:

```math
\delta_L=\frac{\eta}{128K_*}.
```

The added full-frame error is at most
$`32K_*\delta_L=\eta/4`$. The retained real and complex error
budgets consequently give complete initialized-isometry error at most
$`2^{-L}`$, including literal phases, arbitrary dirty inputs and
references, and every work and sole-clean-flag return error. No source
accuracy depending on n, a table row count or an output-bank size is
required by this transfer.

With $`N=2^n`$, every $`n\ge1`$, and the same literal dirty
thresholds $`B\ge17(L+n+7)`$ for real frames and
$`B\ge17(L+n+8)`$ for complex magnitude completion, the conditional
bounds on one circuit are

```math
\begin{aligned}
T&=O_L(\sqrt N+N/B),\qquad G=O_L(N),\\
D_T&=O_L\!\left(
N/B^2+\overline d_{\delta_L}(2n+10)
+\log_2\log_2(n+4)\right).
\end{aligned}
```

There is one external clean qubit and at most B arbitrary dirty qubits.
The source envelope of Section 5 covers every smaller first-chunk
length; the multiplier eight is absorbed in the displayed depth
constant. R98's retained small-dimension and narrow-width branches
absorb the new tolerance-dependent constants without strengthening its
literal allocation thresholds. Its existing logarithm-of-logarithm
term absorbs the exact suffix depth introduced here.

This is a fixed-L conditional statement, not a uniform growing-L
improvement. If the relaxed complete source has sublogarithmic address
depth at every requested fixed tolerance, it improves the fixed-accuracy
frame law through this paid recurrence. Without that source, R96 and
the unrestricted depth gap remain unchanged; the original physical-pair
lower scout is independent of this sufficient upper-bound interface.

## 8. A first test of the growing-clean window

Consider a weight-only phase fingerprint: the address is preserved,
and its only influence is through $`w=|x|\in\{0,\ldots,b\}`$ in
h sequential phase calls. On the entire remaining physical register,
including arbitrary dirty helpers, suppose

```math
U(w)=C_hR_h(w)C_{h-1}\cdots C_1R_1(w)C_0,\qquad
R_t(w)=\sum_{\ell=1}^{d_t}e^{iw\theta_{t,\ell}}P_{t,\ell}.
```

All C_t and spectral projectors are w-independent; the projectors in
each call are orthogonal and sum to identity. The C_t may mix logical,
signal and dirty registers. Suppose $`d_t\le2^k`$. This is an
additional spectral assumption on the **whole phase call**, including
its dirty action; the clean allocation alone does not imply it.

Let R count phase angles after omitting zeros modulo 2π, so $`R\le h2^k`$.
For R≥1 put $`Q=\lfloor b^{1/R}\rfloor\ge1`$. There is some
nonzero weight w≤b for which

```math
\|U(w)-U(0)\|\le\frac{2\pi h}{Q}.
```

To prove this, partition the unit R-cube into $`Q^R`$ boxes and place
the $`Q^R+1`$ points $`j(\theta_1,\ldots,\theta_R)/(2\pi)`$
modulo one, for j from zero to $`Q^R`$. Two share a box. Their
difference gives a positive w≤b with every phase within 2π/Q of a
multiple of 2π. Thus $`\|R_t(w)-I\|\le2\pi/Q`$ for every call;
unitary telescoping proves the bound. The angles may depend on b.

Compression by any common initialized isometry J preserves this bound
for $`K(w)=J^\dagger U(w)J`$. If weight zero and every positive
weight must produce opposite literal phases with amplitude alpha>0
and error epsilon, even on one fixed unit input, then

```math
2\alpha-2\epsilon\le\|K(w)-K(0)\|
\le\frac{2\pi h}{Q},\qquad
\epsilon\ge\alpha-\frac{\pi h}{Q}.
```

If R=0, every call is identity and the error is at least alpha directly.
The bound may be vacuous for Q=1. A complete indicator phase contains
this witness: fix the bank with only its zero-address bit equal to one.
Its desired phase is negative at weight zero and positive otherwise.
The test applies when the implementation on this witness has the
displayed weight-only form; a general indicator need not have it.

There is also a quantitative access tradeoff. Put
$`\Delta=\alpha-\epsilon\gt0`$. For R at least one the displayed
floor gives $`Q\le\pi h/\Delta`$. Since
$`b^{1/R}\lt Q+1`$, every successful fingerprint must satisfy

```math
\begin{aligned}
b&\lt(1+\pi h/\Delta)^R,\\
R\log(1+\pi h/\Delta)&\gt\log b,\\
h2^k\log(1+\pi h/\Delta)&\gt\log b.
\end{aligned}
```

At fixed k and fixed positive Delta this requires
$`h=\Omega(\log b/\log\log b)`$. It leaves a near-logarithmic
but sublogarithmic access window; no construction in that window is
supplied. These are bounds on the specified sequential phase accesses.
Their native implementation, accuracy and simultaneous resource ledger
are additional obligations.

In the window $`2^k=O((\log b)^\beta)`$ and
$`h=O((\log b)^\gamma)`$, with $`\beta+\gamma\lt1`$, one has

```math
Q\ge\tfrac12\exp\!\left(c(\log b)^{1-\beta-\gamma}\right)
```

for some c>0 and sufficiently large b. Hence h/Q tends to zero and
the error is at least $`\alpha-o(1)`$, even granting exact free
phase calls. This stops this compressed weight-only fingerprint from
supplying the inverse-polylogarithmic source precision in Section 6.

Here h counts logical phase accesses. A native ledger is needed to
relate it to $`d_V`$. Additional frequency directions on dirty registers,
more phase accesses, or address and bank motion outside the weight-only
form can escape this test. The resulting source still needs a complete
bank and work-return bound.

<a id="shared-layer-indicator-fusion"></a>

## 9. Share a complete source across one tree layer

Let the p-bit address be divided into nonempty blocks of lengths
$`b_1,\ldots,b_J`$. Write $`s_0=0`$ and
$`s_i=\sum_{h\le i}b_h`$. At stage i the parent bank contains one
arbitrary bit $`z_a`$ for each prefix $`a\in\mathbb F_2^{s_{i-1}}`$.
The new tree bank $`A_{a,v}`$ is indexed by a and
$`v\in\mathbb F_2^{b_i}`$. Borrow **one** old source bank T with
$`2^{b_i}`$ bits. A source $`E_i=E_{b_i,\delta_i}`$ satisfies the
complete zero-clean assumption of Section 1 on the current address
block, T and private work, with error $`\delta_i`$. No completed
address preservation is needed.

Define the exact batch

```math
B_i=\prod_{a,v}\mathrm{CCX}(z_a,T_v;A_{a,v}).
```

Use Hadamards on the child bank around the corresponding CCZ batch.
The [private-target cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md#1-a-complete-diagonal-compiler)
uses one dirty phase target per triple. Its eight T bands act on
distinct targets; both shared controls appear only in Clifford
moves. The phase exponent at the boundary is $`4z_aT_v A_{a,v}`$
modulo eight, independently of the initial helper. Thus each batch has
at most eight T bands, $`8\cdot2^{s_i}`$ T gates,
$`10\cdot2^{s_i}`$ Cliffords, and $`2^{s_i}`$ returned dirty
helpers. These bounds also hold for its actual reversed word.

Use the stage $`S_i=E_i^\dagger B_i E_i B_i`$ in operator order.
With the ideal indicator in both source slots, the two batches add
$`z_a T_v`$ and $`z_a(T_v\oplus[v=x^{(i)}])`$ to each child.
The inverse restores T; the net addition is
$`z_a[v=x^{(i)}]`$. This identity holds for every parent bit
simultaneously. It is literal on the entire physical space. There
is only one E occurrence and its inverse, regardless of the number
of parents. The unitary hybrid gives
$`\|S_i-S_{i,0}\|\le2\delta_i`$ even when the actual E moves
the address or leaves residual work error.

Let $`F=S_J\cdots S_1`$ and $`K=F X_{\varnothing}F^\dagger`$,
where the Pauli flips the root. In the ideal circuit each stage adds
a linear function of its arbitrary parents to the child bank.
Conjugating the root flip therefore flips the root and precisely
the one selected child at each successive level. Every other
input-dependent contribution cancels. Let C copy the leaves into the
external output bank by CNOT. In chronological order
$`C,K,C,K^\dagger`$ adds the selected leaf bit difference to the
output and cancels the entire internal flip. Each correction uses the
actual reversed word. The net approximation satisfies

```math
\begin{aligned}
\|\widetilde I_p-I_p\otimes I_{\mathrm{borrowed}}\|
  &\le 8\sum_i\delta_i,\\
T(\widetilde I_p)&\le 8\sum_i T(E_i)+O(2^p),\\
G(\widetilde I_p)&\le 8\sum_i G(E_i)+O(2^p),\\
D_T(\widetilde I_p)&\le 8\sum_i D_T(E_i)+64J,\\
w(\widetilde I_p)&\le O(2^p)+\max_i\bigl(2^{b_i}+w(E_i)\bigr).
\end{aligned}
```

Here T, G and depth on the right are common source/actual-inverse
bounds, and $`w(E_i)`$ excludes its old output bank. There are four
traversals of F including actual reverses and two source slots per
stage. The non-source counts sum over the number of children,
$`\sum_i2^{s_i}\lt2^{p+1}`$. Tree banks persist across traversals.
The old source bank, private work and cubic flags can be reused. Their
ideal action returns them and the full-space hybrid bounds actual
leakage across reuse. The same bounds, including the error, hold for
the actual native inverse. The compiler has zero initialized inputs.

<a id="fixed-singly-exponential-source-overhead"></a>

### Fixed singly exponential source costs suffice

Suppose that for every fixed requested $`0\lt\delta\lt1`$ a
terminating native compiler supplies complete indicators of arity b
with zero clean inputs, independently of the final tree size, and
common forward/actual-inverse bounds

```math
T(E_{b,\delta}),\,w(E_{b,\delta})\le C_\delta 2^{cb},
\qquad
G(E_{b,\delta})\le C_\delta 4^{cb},
\qquad
D_T(E_{b,\delta})\le d_\delta(b).
```

The constant c is finite and common to all tolerances. For any fixed
target tolerance tau, choose a fixed integer $`J\gt\max(0,c)`$ and
a balanced J-block partition with $`b_i\le p/J+1`$. Request the
fixed tolerance $`\delta_i=\tau/(8J)`$ on each block. The above
ledger gives $`O_\tau(2^p)`$ T gates and borrowed width,
$`O_\tau(4^p)`$ Cliffords, and depth

```math
D_T\le 8J\max_{b\in\{\lfloor p/J\rfloor,\lceil p/J\rceil\}}
               d_{\tau/(8J)}(b)+64J.
```

A bounded set of arities $`p\lt J`$ uses the exact dirty tree.
Consequently $`d_\delta(b)=o(\log(b+2))`$ at each fixed tolerance
suffices. The uniform exponent lets J be chosen before the source
tolerance. If the exponent grows with required accuracy, the
displayed layer ledger still applies but this choice of J needs a
separate bound.

The resulting normalized indicator can replace the one in
[Section 7](#7-conditional-complete-frame-consequence-with-its-source-tolerance).
With $`\tau=2^{-L}/(16K_*)`$, the completed word queries add at
most $`4K_*\tau=2^{-L}/4`$ error. At fixed L, the same real and
complex dirty thresholds and one-clean frame allocation give
$`T=O_L(\sqrt N+N/B)`$, $`G=O_L(N)`$ and, for a sublogarithmic
source, depth $`O_L(N/B^2)+o_L(\log(n+2))`$. The complete frame
comparison includes every logical column and dirty input, references,
literal phase and work return to the stated tolerance.

<a id="fixed-tolerance-shared-layer-source"></a>

### A fixed-error address-preserving source also suffices

Suppose instead that a zero-clean E is supplied at one error
$`0\lt\epsilon_0\lt1/\sqrt2`$ with the bounds above, and that its completed
word exactly preserves each computational address sector on the full
physical space. Internal address motion is permitted. Give the
amplified source a new arbitrary bank Y and one extra dirty bit t.
With the old source bank T define

```math
P_t=\prod_v\mathrm{CCZ}(t,T_v,Y_v),\qquad
V=E^\dagger P_t E P_t.
```

The two occurrences of $`P_t`$ use the same paid shared-control
cubic schedule. At the ideal E they cancel the old bank's contribution
and yield precisely $`(-1)^{tY_u}`$. The actual V preserves each
computational $`(u,Y)`$ sector and has full error at most
$`2\epsilon_0`$.

Apply the [pair-reflection wrapper](PRESCRIBED_SCALAR_TABLES.md#pair-reflection-scalar-amplification)
within those sectors. Its physical signal is the same arbitrary t,
irrespective of the number of labels. Since $`V=I\oplus U_Y`$
in the t sectors, the product $`V^\dagger X_t V X_t`$ is exactly
$`U_Y\oplus U_Y^\dagger`$. It is within $`2\epsilon_0`$ of the sector
sign $`(-1)^{Y_u}`$, so its Hermitian part has the signed gap
$`1-2\epsilon_0^2`$. The wrapper supplies literal phase
$`(-1)^{Y_u}`$ to any fixed tolerance with a tolerance-dependent
number of V and actual-inverse visits and one-qubit gates. Hadamards
on Y before and after recover the new complete XOR indicator, with
the old T and all work included in its full-operator comparison.
The extra T, Clifford and width costs per visit are $`O(2^b)`$;
the extra depth per visit is constant. Taking the exponent at least
one absorbs these costs uniformly in the requested tolerance.

Exact preservation here pins the sign subspaces to the desired
address. Close agreement of E on the whole space alone gives no
such guarantee for the wrapper. For example, if E is the ideal
indicator preceded by a small address-only rotation R, the
t=1 sector of V has the conjugated sign
$`R^\dagger(-1)^{Y_u}R`$. Sign sharpening with the same pair
of reflections retains those rotated eigenspaces. The general
shared-layer theorem above uses full-error estimates directly and
allows such address motion.

<a id="fixed-tolerance-summable-address-source"></a>

### Summable completed address motion also permits sharpening

For a fixed-error source E with $`\epsilon_0\lt1/\sqrt2`$, let
$`E^{(a)}`$ retain its address blocks from u to $`u\oplus a`$.
Suppose a uniform bound $`\sum_a\|E^{(a)}\|\le B_0`$ holds on
the full space, with $`B_0`$ independent of the address length.
A fixed finite number of supported differences suffices. Keep
the same fixed exponential resource bounds and depth premise.

Use the same old bank T, new bank Y and exact word
$`V=E^\dagger P_t E P_t`$. V preserves t and every computational
Y sector, including when E moves the address. Its t-zero sector
is exactly identity, and its t-one sector is a unitary $`U_Y`$
with $`\|U_Y-D_Y\|\le2\epsilon_0`$, where
$`D_Y|u\rangle=(-1)^{Y_u}|u\rangle`$ acts trivially on other work.
Initialize t to zero and apply $`H_t V H_t`$. Its compressed
measurement of $`Z_t`$ gives $`M_Y=(U_Y+U_Y^\dagger)/2`$.
On an arbitrary unit vector in a fixed u,Y sector, expansion
of the squared error gives the signed diagonal bound
$`(-1)^{Y_u}(M_Y)_{uu}\succeq\gamma I`$ with
$`\gamma=1-2\epsilon_0^2\gt0`$. This calculation requires no
commutation of $`U_Y`$ with $`D_Y`$.

The sum of address-grade norms is submultiplicative, since
displacements add under multiplication. The t-one restriction
of $`P_t`$ has zero address displacement. Thus the address-grade
sum of the full compression M, including all Y sectors in the
same operator, is at most $`B_0^2`$. The
[linear recording construction](PRESCRIBED_SCALAR_TABLES.md#summable-address-motion)
gives k initialized tags that retain only a kernel of displacement
grades. Choose $`k=O_{B_0,\gamma}(1+\log(1/\tau))`$ and one
linear tag so that the recorded compression $`M'`$ differs
from its address diagonal by at most $`\delta=\gamma\tau/64`$,
simultaneously on all Y and arbitrary work. Its spectral sign is
within $`\delta/(\gamma-\delta)`$ of the required sign.

Use the same clean-code QSP recovery against the **actual** sign
of $`M'`$, with initialized-isometry error $`\varepsilon`$. After
Hadamards on Y its comparison unitary is within
$`\rho=\delta/(\gamma-\delta)`$ of the specified complete vector
XOR. Apply the [vector XOR promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#4-approximate-clean-return-still-gives-a-full-operator-query)
with k+1 clean bits and a fresh arbitrary output bank. Compare
first to the mathematical clean-sector extension with that
comparison unitary on its clean block. Its full-space discrepancy
costs $`8\cdot2^{k+1}\varepsilon`$ in the emitted promotion.
For the extension itself, the selective echoes cover disjoint exact
clean sectors. In each sector the two appearances of its restricted
unitary cost at most $`2\rho`$ in operator norm by replacement
with the desired vector XOR. Their direct sum has the same bound,
independent of the number of sectors. This argument uses neither
address preservation of the restricted unitary nor separate exact
return of the actual emitted QSP word. A common source phase
cancels inside the actual conjugations. Choosing
$`\varepsilon=\tau/(16\cdot2^{k+1})`$ gives total error at most tau.

At fixed tolerance, all tag and recovery parameters are bounded
independently of the address length. A linear tag uses only
Cliffords. The promotion and its actual inverse use constantly
many source visits, $`O_\tau(2^b)`$ additional gates and dirty
space, and constant additional T-depth beyond the charged visits.
Hence the zero-clean amplified compiler retains a uniform
singly exponential exponent and sublogarithmic depth. The
displacement premise is an additional source condition.
