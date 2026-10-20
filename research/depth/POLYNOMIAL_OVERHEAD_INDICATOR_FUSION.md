# Polynomial source overhead can be absorbed by a paid indicator fusion

[Retained exact dirty tree](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Branch-difference recurrence](DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md) ·
[Native bilinear batch and clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Complete approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Two-band source obstruction](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)

The complete approximate indicator need not first be discovered with
linear T-count and dirty width. Polynomial overhead in its address size
can be removed by using that indicator once near the root of a dirty
tree, then completing only a logarithmically sized address suffix with
the retained exact construction. A native controlled bilinear batch
pays for the root-controlled use of the source. The complete echo has
eight source/actual-inverse appearances and error at most eight times
the source's full-operator error.

This is a constructive reduction between source interfaces. It supplies
no new shallow source and does not change the unconditional depth
frontier. Its first unpaid operation is still a complete approximate
indicator, now with polynomial address overhead allowed in its counts
and width. Its source T-depth must still beat logarithmic address depth.

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

This absorbs polynomial overhead in b, not arbitrary count growth.
An emitted source T-count proportional to $`4^b`$ times a polynomial
cannot fit $`O(2^p)`$ with $`r=O(\log p)`$. Paying that cost by
shrinking b would require $`r=\Theta(p)`$, returning the exact tail
to logarithmic-in-p depth. Likewise substituting the retained exact
$`d_E(b)=O(\log(b+2))`$ indicator into this reduction leaves that
depth order unchanged. The reduction does not turn an existing
logarithmic-depth or exponentially over-budget source into a faster one.

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

There is also a conditional growing-clean window. For constants
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

The first unpaid operation is therefore a complete shallow source E
with the relaxed polynomial overhead of Section 1, or an explicitly
paid clean source satisfying one of the windows above. Requiring
exactly three physical T bands is unnecessary; any source depth that
is little-o of logarithmic address depth suffices. Independent-row
polynomial count overhead alone is no longer a terminal objection to
this relaxed source interface, but rowwise approximation, arbitrary
output action and complete work return still need a joint operator
proof. The existing filter and clean-history failures are not repaired
by declaring the source or its nonlinear bank action available.
The [branch-difference recurrence](DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md#2-error-and-complete-return)
gives a different possible source contract when internal errors cancel
in selected-parent differences; it does not remove the complete-bank
premise used in this reduction.

## 7. Conditional complete-frame consequence with its source tolerance

Apply the [approximate R98 transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md#3-the-bounded-r98-call-count-and-final-error)
to the normalized indicator from this chapter. A completed word query
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

Here h counts logical phase accesses, **not native T-depth**; relating
it to d_V needs a separate native ledger. This test does not exclude
all sources in the growing-clean window. Additional usable frequency
directions in the full dirty-register action, more phase accesses, or
address/bank motion outside the weight-only form can leave its premise.
Such a changed mechanism must still supply complete bank accuracy and
work return; polynomial count overhead is now permitted by the fusion.
