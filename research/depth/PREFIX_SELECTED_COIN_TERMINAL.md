# Prefix-selected terminal compilation with product coins

[Product-coin decoder](COIN_DILATION_DEPTH.md) ·
[Retained recursive terminal](RECURSIVE_TERMINAL_DEPTH.md) ·
[Sharp dirty queries](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence) ·
[Selected routing](SELECTED_DIRTY_ROUTING.md)

**Local result R88.** A preserved external prefix can select an arbitrary
terminal frame while the product-coin construction uses work and
polynomial overhead depending only on the remaining local height.
Recursive protected banks then improve the terminal depth in R81 from
$`O_L(N/b^2+\ell\log_2(M+2))`$ to
$`O_L(N/b^2+\ell\log_2\log_2(M+4))`$.
This chapter is the proof home for that terminal improvement. It does
not claim an $`O_L(\ell)`$ terminal or remove the separate coarse
recovery cost of the complete R87 construction.
The later [R89 prefix-selected recursion](CAUSAL_COMPACT_COIN_DEPTH.md#10-prefix-selected-recursion-and-iterated-logarithm-depth)
replaces that coarse search and attains
$`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$.
This chapter retains its independent terminal proof and supplies the
read-only sharp-query and literal-width interfaces used by R89.

## 1. Terminal theorem and scope

Fix an integer $`L\ge6`$, and put

```math
N=2^n,\qquad \eta=2^{-L},\qquad \ell=\log_2(n+2),\qquad
F(m)=\log_2\log_2(m+4).
```

For every integer $`n\ge1`$, every $`0\le M\le n`$, and every
prescribed product $`W_{\rm tail}`$ of the deepest M layers of a
complete real Hopf frame, there is one Clifford+T circuit V with one
external initialized qubit and at most b arbitrary dirty qubits whenever
the integer b satisfies
$`b\ge17(L+n+7)`$, such that

```math
\|VJ_1-J_1(W_{\rm tail}\otimes I_b)\|_{\rm op}\lt\eta/4,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell F(M)\right).
\qquad(1)
```

The norm includes every column, literal phase, source leakage and work
return, arbitrary dirty inputs, and external references. The constants
are independent of n, M, b, and the angle tuple; they may depend on the
fixed L. For M zero use identity. The unchanged upper $`n-M`$ bits
may select unrelated local frames, and need not be zero.

For $`M\le\min\{n,C_0\ell\}`$ with fixed $`C_0`$, (1) gives

```math
D_T=O_{L,C_0}\!\left(N/b^2+\ell\log_2\log_2(\ell+4)\right).
\qquad(2)
```

The angle interface is the existing terminating certified evaluation
interface. Signed and singular inputs are included. Every query, source,
gradient preparation, actual inverse, predicate, and routing operation
is charged. There is no QRAM, measurement, reset, supplied precision
state, or free coherent evaluator. No efficient classical preprocessing,
scalable native emitter, or precision-uniform strengthening is asserted.

## 2. Exact factorization and where rounding occurs

When m bottom layers remain, the first $`n-m`$ logical bits form a
preserved prefix P. Choose a bottom bank Z of r logical bits, and put
$`h=m-r`$. The next h layers have the exact target

```math
E_m=\sum_p|p\rangle\langle p|_P\otimes
\left(W_{p,h}\otimes\Pi_Z+I_X\otimes(I_Z-\Pi_Z)\right),
\qquad \Pi_Z=|0^r\rangle\langle0^r|.
\qquad(3)
```

The remaining terminal operator factors as $`W_m=W_rE_m`$,
where the enlarged preserved prefix for $`W_r`$ has length $`n-r`$.
This is R81's [prefix-selected factorization](RECURSIVE_TERMINAL_DEPTH.md#2-a-decoder-with-a-preserved-external-prefix).
Save $`H=[Z=0^r]`$ using the sole external initialized qubit before
using Z as work; reverse that actual predicate word last.

First suppose all prescribed pairs are positive magnitudes. At one
stage let l be its accuracy parameter and set

```math
\zeta=2^{-l},\qquad K=m+l+2,\qquad
\tau=\left\lceil2\log_2(m/\zeta)\right\rceil+12.
\qquad(4)
```

Apply R86's [certified finite rounding](INTERVAL_DILATION_DEPTH.md#2-finite-rounding-and-the-atom-embedding)
only to the pairs in $`W_{p,h}`$, uniformly over p. It produces
strictly interior probabilities $`k_v/2^\tau`$, with
$`1\le k_v\le2^\tau-1`$, and

```math
\|E_m-E'_m\|\le2h\,2^{-\tau/2}\le\zeta/32.
\qquad(5)
```

Every source, projector, and coin rectangle within this stage uses the
same rounded early tuple. The unprocessed factor $`W_r`$ retains its
prescribed coefficients until its own stage. Thus later stages may use
different rounding denominators. Equation (3) and unitary telescoping,
not a presumed globally common denominator, justify this choice.
The error is a maximum over orthogonal P sectors, not their sum.

The intended application to R87 already has a positive tuple. To cover
standalone signed or singular terminal inputs, apply the same finite
rounding/sign reduction once to the whole selected M-level family,
with error at most $`\eta/512`$. Its exact full-column sign gauge is

```math
W'_{\rm tail}=D_{\rm sign}W_{\rm pos}S_{\rm sign}.
\qquad(6)
```

The signs are finite preprocessing data, including below originally
zero incoming amplitudes; no exact sign or zero test on an original
analytic coefficient is required. Each diagonal is one exact Boolean
phase on the full n-bit logical address. Pay these two phases once.
Subsequent stage rounding applies only to positive magnitudes. Paying
an unsuppressed global sign query at every recursive stage would add
an unjustified count factor and is not this construction.

## 3. A prefix-selected product-coin stage

Reserve disjoint short metadata of size $`O(l+\log_2(m+2))`$ and
free donor bits

```math
m_{\rm free}=\left\lceil c_0[l+\log_2K]\right\rceil,
\qquad m_{\rm free}\ge4(\tau+2).
\qquad(7)
```

The absolute constant $`c_0`$ can ensure the last inequality. Metadata
contains the old fine tag, coarse label, search and threshold records,
source index/program/gradient, and other fixed interfaces. None belongs
to the free donor reserve. Choose r after these constants and the
resource degree below, so that both the metadata and (7) fit in Z.

Use the R87 bins within the local h-bit early register. Starting
backward at d=h, repeatedly take

```math
v=h-d+m_{\rm free},\qquad
q=\min\left\{d,\left\lfloor\frac{v}{2(\tau+2)}\right\rfloor\right\},
\qquad a=d-q,\qquad c=q\tau.
\qquad(8)
```

Here c is the number of coin bits, and r still denotes the protected
remainder. The local descendant inverse exposes $`h-d`$ zero bits,
and

```math
v\ge2(c+2q),\qquad a+q+c=d+c\le h+m_{\rm free}.
\qquad(9)
```

For each bin and prefix, the R87 rectangle source preserves the local
q-bit leaf label x and prepares the normalized indicator of its
product rectangle on c coin bits. The exact membership query then
clears x. If E denotes this isometry and $`\mathcal H`$ the tensor
Hadamard, the last nonzero complete coin block gives the local spectral
depth exactly. The physical chronological word is

```math
U_E,\quad\mathcal H,\quad
T\mathrel{\oplus}=H\bigl(a+\lambda(C)\bigr),\quad
\mathcal H,\quad U_E^\dagger.
\qquad(10)
```

The inverse in (10) is the actual reversed membership/source word.
R87's [initialized-isometry estimate](COIN_DILATION_DEPTH.md#4-paid-physical-embedding-and-actual-inverse)
gives error at most $`2\delta_E`$, for arbitrary old-tag input on
the active initialized sector. This includes all coin/source leakage
and return. The middle label is zeroed on H zero, so the full sandwich
is exactly identity there even on occupied work. Hadamards on inactive
fixed-bank padding commute through the selected middle action and
cancel. No source completion on occupied initialized work is assumed.

One completed selected layout surrounds (10); its actual inverse
restores all parked arbitrary bank values and their references. The
layout acts only on local logical bits and the free reserve. It excludes
P and all metadata. The retained selected router suffices: its boundary
depth is $`O(\log(m+2)\log\log\log(m+8))`$. This proof needs no
new faster routing lemma.

Including P, the relevant large data addresses are at most

| Query | Number of relevant data bits |
|---|---:|
| Rectangle family source | $`(n-m)+a+q+c`$ |
| Membership | $`(n-m)+a+c`$ |
| Fine-label readout | c |

By (9) all are bounded by $`n-r+m_{\rm free}`$.
Short control and expansion labels are charged separately below.
The membership output lies in the fixed external local-leaf bank,
disjoint from the union of every bin's original-prefix and coin address
wires. Fine-tag and source-program outputs have their own disjoint
targets. This full union exclusion permits the parallel selected-query
echo; selected-bin disjointness alone would not suffice.

Run the complete R87 early word: initial tagged source, coarse-label
computation, selected descendant inverse, (10), descendant reload,
endpoint uncomputation, and coarse-label erasure. The latter uses
binary search on transported scaling projectors, with the original
fine tag and coarse label excluded from its input. P is preserved and
is never tested for zero. Padded endpoint comparisons are always true;
the endpoint-h reflection is
$`I-2(A_h\otimes|0\rangle\langle0|_A)`$, hence minus identity
on initialized A, rather than an identity reflection.

Let p_bin be the number of bins and $`k=\lceil\log_2p_{\rm bin}\rceil`$.
The growth in (8) gives

```math
p_{\rm bin}=O(\tau\log_2(m+2)),\qquad
k=O(1+\log_2\tau+\log_2\log_2(m+2)).
\qquad(11)
```

Set the completed source tolerances to

```math
\delta_0=\delta_U=\delta_E=\zeta/128,\qquad
\delta_R=\frac{\zeta}{512(k+1)}.
\qquad(12)
```

R87's full-word hybrid gives implementation error
$`\delta_0+2\delta_U+2\delta_E+8k\delta_R\lt\zeta/8`$.
Together with (5), one stage differs from (3) by less than
$`\zeta/4`$. Its finest source precision is
$`O(l+\log(k+1))`$, including the reflection tolerances. All
conditionally initialized source work is $`O(l+\log(m+2))`$.
Every inactive completed word has literal identity action. The actual
saved-H inverse propagates leakage unitarily; reusing H at the next
stage requires no reset or assumption of exact actual return.

## 4. Parallel sharp queries with a read-only public address

For a one-bit table with Q rows and a disjoint arbitrary dirty pool B,
the sharp blocked query has

```math
\begin{aligned}
T_Q&=O(1+\sqrt Q+Q/B),\qquad G_Q=O(Q),\\
D_{T,Q}&=O(1+Q/B^2+\log_2(s+2)),\\
s&=\lceil\log_2Q\rceil.
\end{aligned}
\qquad(13)
```

provided $`B\ge C_Q(s+2)`$ for an absolute constant $`C_Q`$.
Here the address interface is strengthened to be read-only at every
non-Clifford boundary. The low nonuniform indicators already use
private returned conjunction circuits, whose public bits enter only
as CNOT controls. In the blocked query's
[high-address dirty traversal](../../docs/BLOCKED_BILINEAR_LOOKUP.md#2-select-a-bilinear-block-with-dirty-traversal),
replace each two-control selector update by the k=2 instance of the
[read-only conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences).
Its two inputs are the private preceding selector and a public positive
or negative address literal. It has constant count, constant T-depth,
and a constant-size private dirty pool, returned before reuse. Negative
literals use the same read-only literal interface. Root updates are
Clifford operations. Bilinear leaf actions touch private indicators,
selectors, outputs, and helpers. Enlarging $`C_Q`$ by an absolute
constant therefore proves (13) with public wires appearing only as
Clifford controls, without changing its asymptotic count or depth.

For each Boolean-query position in the common source grammar, run its
per-bin, per-output one-bit instances in parallel with private dirty
pools. Non-Clifford supports are disjoint. Shared-address Clifford
operations may serialize; their elementary gate counts are summed.
If Q denotes that parallel query word and C copies each private dirty
output into the shared final target by CNOT, execute
$`Q,C,Q^\dagger,C`$. All dirty offsets cancel and every helper
returns, as in the
[ragged selected-query lemma](SELECTED_DIRTY_ROUTING.md#2-parallel-selected-tables-with-different-relevant-addresses).
The union address/output exclusion in Section 3 is essential here.

Whole sources sharing quantum targets are not run in parallel. Their
common grammar selects active target prefixes, exact zero reflections,
and returned short masks. Only each Boolean-query position is
parallelized. In particular a long membership output does not inherit
the serial output-length term of a multioutput blocked query. Its
one-bit pools are simultaneously allocated and counted. Identity
padding aligns different native schedules.

## 5. Explicit local polynomial envelope

Put

```math
\overline Q=N2^{-(r-m_{\rm free})},\qquad K=m+l+2.
```

We now fix the degree before choosing the protected margin. Constants
in this paragraph are absolute; choosing the eventual fixed-L schedule
only strengthens the estimates. Equations (4), (11), and (12) give
$`\tau\le CK`$, $`p_{\rm bin}\le CK^2`$, $`k\le CK`$,
source precision $`p_{\rm src}\le CK`$, and expansion length
$`J_{\rm exp}\le CK`$. There are at most CK output bits at one
query position, hence at most $`I=CK^3`$ simultaneous one-bit
bin/output instances. In particular $`q_{\max}\le m`$ and
$`c_{\max}\le m+m_{\rm free}\le CK`$; the short source bank
and each calibration output also have at most CK bits.

The following deliberately loose ledger is sufficient.

| Factor | Bound |
|---|---:|
| Relevant data configurations in one bin | $`\overline Q`$ |
| Short bin, threshold/tag, expansion combinations | $`CK^4`$ |
| Full source-work configurations, if an anchor is tabulated | $`C4^lK^3`$ |
| Rows of any one-bit large-table instance | $`CK^8 4^l\overline Q`$ |
| Simultaneous one-bit instances | $`CK^3`$ |
| Sequential large-table positions, for count bounds | $`CK`$ |

The source-work entry follows from the retained source-bank estimate
$`2^u\le C\delta^{-2}p_{\rm src}`$, since (12) gives
$`\delta^{-1}\le CK2^l`$. This charges even a tabulated all-work
anchor; structural conjunction implementations are smaller. These
source-work patterns are address rows, not simultaneous output
instances. The factor $`4^l`$ therefore does not enter I.
The source has a constant number of table positions, since its
coherent expansion has an index and is not J_exp sequential queries.
The number of source appearances is $`O(k+1)`$, including every
actual inverse and the complete coarse search.

Allocate $`B_{\rm inst}=\lfloor b/(C'I)\rfloor`$ to each
simultaneous instance, reserving a fixed fraction of b for their total.
Once the width condition of Section 7 holds, (13) bounds one parallel
query position by

```math
T\le C\left[K^3+K^7 2^l\sqrt{\overline Q}
                  +K^{14}4^l\overline Q/b\right],
\qquad G\le CK^{11}4^l\overline Q,
```

```math
D_T\le C\left[\ell+l+\log_2(m+2)
                  +K^{14}4^l\overline Q/b^2\right].
\qquad(14)
```

The square-root term uses $`I\sqrt{Q_{\rm inst}}`$; the
linear term uses $`I^2Q_{\rm inst}/b`$. The parallel depth uses
$`I^2Q_{\rm inst}/b^2`$, with no extra factor I. Summing at most
CK sequential positions raises these powers to at most fifteen.
For the non-table depth, retain the sharper $`O(k+1)`$ position
count, rather than replacing it by CK.

The remaining local tables depend only on a short depth or bin label:
tag packing has at most $`2^h`$ rows and is included in the large
ledger; coarse-label/endpoint/threshold tables have at most $`CK^2`$
rows; selected masks have at most $`CK^2`$ rows and CK outputs.
They do not read all metadata jointly. Applying the conservative
parallel-output interface to these short fields, with all bins and
inverses included, costs at most $`CK^8`$ T gates, Clifford gates,
and dirty helpers. Structural predicates have at most CK literals;
their returned read-only helpers and the selected permutation networks
also fit this polynomial. Gradient/coefficient synthesis and modular
shifts have polynomial count in p_src and depth $`O(p_{\rm src})`$.
Thus $`CK^{10}4^l`$ safely bounds all standalone counts and local
dirty work. Their depth is at most
$`C(k+1)(l+\log(m+2))`$ plus the stated router depth.

Consequently **degree twenty** covers one complete stage:

```math
\begin{aligned}
T_i&\le C\left[K^{20}2^l\sqrt{N2^{-(r-m_{\rm free})}}
       +K^{20}4^lN2^{-(r-m_{\rm free})}/b+K^{20}4^l\right],\\
G_i&\le C\left[K^{20}4^lN2^{-(r-m_{\rm free})}+K^{20}4^l\right].
\end{aligned}
\qquad(15)
```

On the fixed-L schedule below, $`l=O_L(\log(m+2))`$ and
$`k=O_L(F(m))`$. Including the router, whose local depth is smaller,

```math
D_{T,i}\le C_L\left[
(\ell+l+\log_2(m+2))(1+F(m))
+K^{20}4^lN2^{-(r-m_{\rm free})}/b^2\right].
\qquad(16)
```

There is no polynomial in global n in (15). The global address appears
in (13) only through its logarithmic depth and its minimum dirty-pool
requirement. The elementary Clifford counts in the sharp query already
include every high-address traversal and private coordinate change.

## 6. Recursive banks, errors, and resource sums

Starting from $`m_0=M`$, choose at each nonfinal stage

```math
l_i=L+\lceil\log_2(m_i+2)\rceil+10,\qquad K_i=m_i+l_i+2,
```

```math
r_i=\left\lceil A[l_i+\log_2K_i]\right\rceil,
\qquad m_{i+1}=r_i.
\qquad(17)
```

Fix the absolute constant A after the metadata constant, c_0, and
degree twenty. Increase it enough that r_i contains all metadata and
free donors and that

```math
r_i-m_{{\rm free},i}\ge100[l_i+\log_2K_i]+O(1).
\qquad(18)
```

Choose a fixed stopping height $`m_*=m_*(L,A)`$ such that every
$`m>m_*`$ satisfies

```math
r(m)\le m/2,\qquad F(r(m))\le\tfrac12F(m).
\qquad(19)
```

Such a finite cutoff exists because $`r(m)=O_L(\log m)`$.
Below it, stop and use R81's retained smaller-error constant-height
bottom. No asymptotic assertion is made about the practical size of
this cutoff. Every nonfinal protected bank fits within the m remaining
logical bits; no external initialized bank is added.

The complete stage error, including (5), is less than

```math
2^{-l_i}/4\le\eta/[4096(m_i+2)].
```

The heights increase at least geometrically when read backward, so
$`\sum_i1/(m_i+2)\le2`$. Thus all recursive errors sum to at most
$`\eta/2048`$. For the remaining prescribed positive layers use
R72's smaller-error allocation $`q_j=L+j+7`$, whose error is strictly below
$`15\eta/64`$. Even including the optional single sign-rounding
error $`\eta/512`$, the total is strictly below $`\eta/4`$.
When (6) was used, these remaining layers belong to $`W_{\rm pos}`$,
not to the original signed tuple; they have not been rerounded at any
earlier recursive stage.
Unitary telescoping uses ideal initialized stage inputs for each local
estimate and propagates previous actual leakage unchanged. It does
not assert that an approximate preceding stage has returned exact zeros.

Equation (18) gives

```math
K_i^{20}2^{l_i}2^{-(r_i-m_{{\rm free},i})/2}\le K_i^{-10},
\qquad
K_i^{20}4^{l_i}2^{-(r_i-m_{{\rm free},i})}\le K_i^{-20}.
\qquad(20)
```

Both inverse-power sums are bounded by (19). Therefore (15)–(16)
sum to $`O_L(\sqrt N+N/b)`$ T count, $`O_L(N)`$ Clifford
count, and $`O_L(N/b^2)`$ table-depth penalty. Standalone terms
obey

```math
\sum_iK_i^{20}4^{l_i}\le C_L(M+2)^{22}
\le C_L(n+2)^{22}=O_L(\sqrt N),
\qquad(21)
```

and also fit $`O_L(N)`$ Cliffords. The unbounded number of recursive
stages has not been multiplied into an unsuppressed count term.

The second contraction in (19), with $`F(m)\ge1`$, gives

```math
\sum_i(1+F(m_i))=O_L(F(M)).
```

Since $`l_i+\log_2(m_i+2)=O_L(\ell)`$, the non-table terms in
(16) sum to $`O_L(\ell F(M))`$. The constant-height final bottom
has $`T=O_L(\sqrt N+N/b+n^2)`$, $`G=O_L(N)`$, and
$`D_T=O_L(N/b^2+\ell)`$. The optional two sign phases have those
same count/depth orders by (13). These bounds establish (1).

## 7. Literal width threshold and finite dimensions

For sufficiently large n, split at $`b_c=\sqrt N/n`$.
On $`b\ge b_c`$, reserve a fixed fraction of b for the sum of all
simultaneous one-bit query pools and another for local arbitrary dirty
work, external fixed banks, masks, and returned helpers. Every stage
has $`m_i\le M\le n`$, $`l_i=O_L(\log(n+2))`$, and
$`I\le C_L(n+2)^3`$. All its local dirty work is bounded by one
fixed polynomial in n. Its public query address has at most
$`n+O_L(\log(n+2))`$ bits. Exponential dominance provides one
cutoff $`n_0(L)`$, independent of M, b, and the angles, such that
all local work fits and every private query pool satisfies (13).
Query pools are reused between serial positions and stages; simultaneous
pools are summed. The final bottom's sharp-query pool also fits beyond
this cutoff. All these pools are arbitrary dirty work, disjoint from
live public addresses, outputs, and the protected logical bank.

For $`b\lt b_c`$, use the retained ungrouped smaller-error terminal
circuit for the deepest M layers, at the original prescribed angles.
Its reservation is

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad B'_0=L+n+9.
```

At $`b\ge17(L+n+7)`$ its remaining query pool satisfies
$`B=b-B'_0\ge16b/17-2\ge b/2`$. This is the explicitly proved
[R81 fallback](RECURSIVE_TERMINAL_DEPTH.md#6-all-dirty-widths-at-the-literal-threshold),
not a generic theorem invoked at a larger L and hence a larger width
threshold. It has error below $`\eta/4`$ and the same count bounds,
while

```math
D_T=O_L(N/b^2+nM+M^2)=O_L(N/b^2)
```

because $`M\le n`$ and $`N/b^2>n^2`$. For every n below the
fixed dimension cutoff, use that same fallback for every b. Its
bounded remaining depth is absorbed by $`O_L(\ell F(M))`$,
uniformly over those finite pairs $`(n,M)`$ and the angle tuple,
even for arbitrarily large supplied b. Unneeded wires remain untouched.

## 8. What remains for logarithmic depth

For a logarithmic terminal M, (2) removes one logarithm from R81's
overhead but leaves $`\log\log(\ell+4)`$. This is an upper bound
for the explicit construction, not a depth lower bound. Even if the
coarse search within every stage were replaced by a constant-call
method, the present recursion would still have an unbounded number of
serial global-prefix-dependent source episodes. The current query
interface charges $`O(\ell)`$ to each. Declaring their iterated-log
number constant would not prove an $`O_L(\ell)`$ terminal.

R89 subsequently supplies the constant-call causal coarse decoder and
obtains an iterated-logarithm stage bound. It still has an unbounded
number of serial global-prefix episodes. The logarithmic target
requires a different schedule, a paid amortization
of those global-prefix episodes with returned work, or a construction
that avoids the recursion. Neither a generic dirty cache nor a
fixed-accuracy theorem used at silently varying precision supplies it.
The separate full-frame coarse recovery and selected routing must also
be included in any improved complete-frame claim.

The [bounded checks](../../tests/test_prefix_selected_coin_terminal.py)
compare independently constructed selected-frame columns with ordered
native layer factors, test stage-specific rounding on every column,
and retain occupied protected sectors. They do not emit the large
sources or prove the asymptotic allocation; those are the arguments
above.
