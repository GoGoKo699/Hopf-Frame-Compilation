# General frame depth by causally batched compact coins

[Product-coin theorem](COIN_DILATION_DEPTH.md) ·
[Compact fixed coins](COMPACT_COIN_FLAGS.md) ·
[Short-label subset routing](SHORT_LABEL_SUBSET_ROUTING.md) ·
[Native family sources](COHERENT_DESCENDANT_LOADER.md) ·
[Coin terminal](PREFIX_SELECTED_COIN_TERMINAL.md)

**Local result R89.** At fixed accuracy, the general complete-frame
compiler has depth
$`O_L(N/b^2+\ell(1+\log_2^*(n+2)))`$, with the retained T-count,
Clifford-count, one-clean interface, and literal widths. Here
$`N=2^n`$, $`\ell=\log_2(n+2)`$, and $`\log_2^*`$ counts binary
logarithms until the argument is at most two.
The construction replaces R87's final serial coarse-label search by one
compact product-coin embedding, a parallel fixed reference transform,
and their actual inverses. The early coarse decoder has
$`O_L(\ell)`$ depth. A prefix-selected recursion leaves only a
polylogarithmic terminal height at each stage. Sections 1–9 first prove
the simpler $`O_L(N/b^2+\ell\log_2\log_2(\ell+4))`$ bound;
Section 10 proves the displayed iterated-logarithm refinement.

This is the proof home for the causal batching construction. Its
matching is local at each original prefix; an arbitrary global
weight-preserving matching would preserve the root but need not preserve
the coarse frame flags. No efficient classical preprocessing, scalable
native emitter, precision-uniform strengthening, unrestricted lower
bound, or optimal-depth assertion is made.

## 1. Theorem and retained early interface

Fix integers $`n\ge1,L\ge6`$, put $`\eta=2^{-L}`$, and let W be
any prescribed complete real frame. For every integer
$`b\ge17(L+n+7)`$, there is a Clifford+T circuit V with one external
initialized qubit and at most b arbitrary dirty qubits such that

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(n+2))\right).
\qquad(1)
```

The prescribed complex magnitude completion, including supplied leaf
phases, satisfies the same statement at $`b\ge17(L+n+8)`$.
All bounds hold on the same circuit. The norm includes every column,
literal phase, source leakage and work return, arbitrary dirty inputs,
and external references. There is no measurement, reset, supplied
precision state, QRAM, or free coherent evaluator.

Use R87's [single finite positive rounding](COIN_DILATION_DEPTH.md#1-theorem-and-inherited-interface),
with

```math
\tau=\left\lceil2\log_2(n/\eta)\right\rceil+12,
\qquad p_v=k_v/2^\tau,\quad1\le k_v\le2^\tau-1.
\qquad(2)
```

Its complete-operator error is at most $`\eta/32`$. The exact finite
sign gauges and the supplied complex leaf phases are restored in
Section 9. All early sources, probabilities, and coins below refer to
this same rounded tuple, including below originally zero amplitudes.

Reserve a protected logical suffix of length

```math
M=\lceil A_L\ell^3\rceil,\qquad h=n-M,
\qquad m_0=\lceil c_L\ell\rceil.
\qquad(3)
```

The constants needed for metadata, source accuracy, and table overhead
are fixed before choosing the independent suppression constant
$`A_L`$. Work beyond the resulting fixed L-dependent dimension cutoff;
Section 9 supplies the finite and narrow-width branches.

On the h early bits retain R87's geometric bins
$`0=d_0\lt d_1\lt\cdots\lt d_p=h`$, obtained backwards from
$`q=\min\{d,\lfloor(h-d+m_0)/(2(\tau+2))\rfloor\}`$.
They satisfy

```math
p=O_L(\ell^2),\qquad q_i=d_i-d_{i-1}\ge1.
\qquad(4)
```

Run the same tagged early source, coarse-bin computation, selected
endpoint inverse, local product-coin fine erasure, and endpoint reload
as [R87 Section 6](COIN_DILATION_DEPTH.md#6-complete-early-word-and-simultaneous-global-budget).
The fine tag is then zero and the data is restored. Only the old coarse
label remains. The first old bin contains depths $`0,\ldots,d_1`$;
old bin $`i-1`$ contains depths $`d_{i-1}+1,\ldots,d_i`$.
This chapter replaces its final binary-search cleanup. In the retained
fine episode, implement both local/donor interval exchanges by the
[logarithmic-depth contiguous router](SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals),
including its actual reverse and loaded interval parameters. Those
intervals already satisfy R87's disjointness and donor bounds; the
extra live parameter words have $`O_L(\ell)`$ bits. This makes the
retained fine episode $`O_L(\ell)`$ depth as used in Section 8.

The saved activity bit is $`H=[Z=0^M]`$, computed before protected Z is
used as work. Every new primitive is exactly identity when H is zero,
including on arbitrary occupied source, metadata, and dirty registers.

## 2. Compact local coins and positive atom classes

For block i and an original preceding prefix u of length $`d_{i-1}`$,
let $`p_{i,u,x}>0`$ be the conditional probability of its
$`M_i=2^{q_i}`$ coarse children. Its dyadic denominator is
$`K_i=2^{s_i}`$, where $`s_i=q_i\tau`$.

Use the [acyclic two-outcome construction](COMPACT_COIN_FLAGS.md#2-a-second-finite-construction-with-two-outcomes-per-column).
Its fixed coin is a uniform M_i-valued column times the geometric
binary coin with positive weights

```math
\frac12,\frac14,\ldots,2^{-s_i},2^{-s_i}.
```

Thus its number of positive atoms and physical coin length are

```math
D_i=M_i(s_i+1),\qquad
r_i=q_i+\lceil\log_2(s_i+1)\rceil=q_i+O_L(\ell).
\qquad(5)
```

Write its fixed positive atom weights as $`w_{i,z}`$ and its state as
$`\gamma_i=\sum_z\sqrt{w_{i,z}}|z\rangle`$, padding unused binary
positions with zero. At each prefix u the finite alias assignment
partitions these atoms into fibers $`\mathcal A_{i,u,x}`$ with total
weight $`p_{i,u,x}`$.

A weight class contains all atoms with a given geometric position,
merging the two final positions if desired. Every such class has at
least M_i atoms, independently of u and the conditional probabilities.
This cardinality, rather than a lower bound on individual branch
probabilities, supplies the density used in Section 6.

## 3. Shared hash catalogues and local codes

Put $`R_{i,u,x}=|\mathcal A_{i,u,x}|`$ and
$`g=\lceil\log_2R_{i,u,x}\rceil`$, with g zero for a singleton.
For fixed i,u,g let $`S_{i,u,g}`$ be the corresponding child set and
$`B=|S_{i,u,g}|`$. For a nonempty set use

```math
k=\min\{q_i,\lceil\log_2(4B)\rceil\}.
```

The fiber-cardinality groups satisfy

```math
B2^g\le2D_i,\qquad k+g\le q_i+O_L(\ell).
\qquad(6)
```

Apply the [staged shear code](SHORT_LABEL_SUBSET_ROUTING.md#3-a-fixed-width-staged-shear-compressor)
to each group. Physically define Low as $`Q_i[0:k]`$ and High as
$`Q_i[k:q_i]`$; these names specify the split and need not denote
binary significance. Write $`x=(H_x,L_x)`$ for the two parts in the
shear formula. Its code is the low hash together with the group, split,
stage, and selected hash seed.
The case $`k=q_i`$ uses the identity and a fixed null stage.

The hash catalogue can be common to every preceding prefix. To see
this, fix i,k and one stage, and condition on the earlier catalogues.
Every remaining prefix/group set is now fixed. A uniformly random
$`k\times(q_i-k)`$ binary shear matrix isolates at least half its
points with probability at least one half: the expected number of
collision victims is at most one quarter of its size. A catalogue of
$`O(d_{i-1}+\log(n+2))`$ independent matrices therefore contains a
good one for all the at most $`2^{d_{i-1}}\mathrm{poly}(n)`$
remaining families simultaneously. This follows by a union bound;
choose a successful finite catalogue and proceed to the next stage.
There are at most $`q_i+1`$ stages and splits. The total catalogue size
is polynomial in n, so every selected seed has $`O_L(\ell)`$ bits.

For an assigned child, store metadata
$`(g_i,k_i,j_i,l_i)`$ and low code

```math
z_i=L_i\oplus M_{i,k_i,j_i,l_i}H_i.
\qquad(7)
```

Given the original preceding prefix and this code, the child is unique.
The selected seed permits a physical shear without rereading that
prefix. It need not be invariant under the eventual matching.

## 4. One global compression word with disjoint query banks

For each block allocate a physical register Q_i consisting of its
original q_i logical bits, $`O_L(\ell)`$ extra zero slots, and
$`O_L(\ell)`$ metadata slots in fixed positions. Choose the constants
so that (6) fits a rank window after the low code. Let Q be their union.
Then

```math
m:=|Q|=h+O_L(p\ell)=h+O_L(\ell^3).
\qquad(8)
```

Activity H, the old coarse tag, every source bank, and every query
helper lie outside Q. The extra Q coordinates come from the protected
suffix and are zero on the active initialized sector.

First, one parallel-output query on the original whole leaf computes
all metadata in (7). Every original prefix is still available then.
For each block, route its low and high intervals into two fixed,
separate external arbitrary dirty banks. The old bank contents are
parked in the original intervals. The
[contiguous interval routers](SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
for different blocks have disjoint arrays and pools and run in parallel,
with $`O_L(\ell)`$ total depth.

Apply all local shears to the Low banks. Each reads only its short
catalogue index and its own High bank; its Low target bank is disjoint.
The shears use private query pools and run in parallel. The tuple of
Low codes and metadata is injective: decode block 1, use its original
child as the prefix for block 2, and continue. Hence one global finite
XOR table on that tuple supplies every original High part. XOR it into
the High banks, clearing their active portions. Reading the whole
Low-bank union costs at most h data bits; dirty padding is ignored.

Reverse the actual routers. The low codes and cleared highs return to
the original block positions, and every external bank is restored,
including its correlations with references. The query read/target
separation is global, not only true on a selected basis row:

| Operation | Read address | XOR target |
|---|---|---|
| Metadata | Original whole leaf | Fixed metadata positions |
| Local shear | Short seed and fixed High bank | Fixed Low bank |
| High erasure | All Low banks and metadata | All High banks |

The retained metadata is part of the compressed code and later part of
the custom atom address. It is not claimed to return at this stage.
Every row has a total finite extension, and the whole word has an
actual inverse. A zero inactive row in every query, together with the
guarded routers, makes the complete word identity on H zero even when
all these banks are arbitrarily occupied.

## 5. One joint source on metadata-selected rank windows

Inside Q_i, place its g_i rank bits in the cleared window
$`[k_i,k_i+g_i)`$, using extra zero slots when necessary. Fixed metadata
positions are excluded from all rank windows. For each decoded child,
enumerate its original local fiber in a fixed finite order. On the
joint rank-window union prepare amplitudes

```math
\prod_{i=1}^p
\sqrt{\frac{w_{i,z(i,u_i,x_i,t_i)}}{p_{i,u_i,x_i}}},
\qquad 0\le t_i\lt R_{i,u_i,x_i},
\qquad(9)
```

and zero on padded ranks. This is one normalized conditional source,
whose preserved label is the entire injective compressed-code tuple.
It is not a sequence of p source calls.

Here is the required extension of the retained
[Boolean-phase source grammar](COHERENT_DESCENDANT_LOADER.md#3-constant-query-flagging-and-literal-inner-amplification).
Replace its selected suffix by the metadata-selected rank union.
Every large Boolean phase may read the one physical m-bit Q register,
plus the ordinary source index and activity bit. Hadamard masks read
only the fixed preserved metadata and H; their possible target wires
exclude metadata positions entirely. Selected zero anchors test only
the rank union and the stated source fields. Calibration tables depend
only on the preserved code and metadata and ignore rank coordinates.
Program, attenuation, gradient, and phase targets remain separate
from Q. Nonselected padding is preserved and ignored by these tables.

The original proof uses suffix geometry only in those masks and
zero-anchor predicates. Its coefficient expansion, literal inner
amplification, source approximation, and actual inverse are unchanged.
In particular there is no full h-bit Low bank plus another full h-bit
rank bank in a source address: the address is Q itself, of size (8).
No new source oracle or clean output register is assumed. Formally,
fixed metadata $`\mu`$ selects a target set $`S(\mu)`$ disjoint
from metadata. The source labels are
$`(\mu,Q\setminus S(\mu))`$. Different metadata sectors are
orthogonal; at fixed metadata every macro preserves the complement.
This is an ordinary direct-sum family source on variable target sets.

Let $`U_E`$ include the exact compressor, this physical source, and
the exact matching below. Its only approximation is the joint source,
with uniform initialized-isometry error $`\delta_E`$. Orthogonal
preserved code sectors take a maximum, so this is the operator error
on all logical columns and external references, not an average under
the state-column distribution.

## 6. Three queries batch the causal local matchings

Each reached custom atom in Q_i consists of metadata, low code, and
rank. Its weight is inherited from its enumerated original atom.
For fixed original prefix u and weight class t, let $`S_{i,u,t}`$
be these custom locations. Let $`T_{i,t}`$ be the standard fixed-coin
locations of that weight, with all unused Q_i positions zero. They
have equal cardinality at least M_i. Since
$`|Q_i|=q_i+O_L(\ell)`$, their density obeys

```math
\alpha_{i,u,t}\ge(n+2)^{-C_L}.
\qquad(10)
```

There is no probability cutoff or rare-prefix exception. Q_i was fixed
using only rank/code/hash metadata and padding. The matching index is
stored outside Q, so its later catalogue size does not enlarge Q_i or
feed back into this density estimate.

The [dense XOR matching lemma](SHORT_LABEL_SUBSET_ROUTING.md#1-a-short-shift-catalogue-between-dense-sets)
can use one catalogue common to every u,t at this block. Union its
spectral and degree probability bounds over at most
$`2^{d_{i-1}}\mathrm{poly}(n)`$ families; a polynomial number
of shifts still suffices by (10). Choose a finite successful catalogue
$`v_{i,1},\ldots,v_{i,e_i}`$, where $`e_i=\mathrm{poly}(n)`$.
Its vectors depend only on i and their short index. The selected
matching $`f_{i,u,t}`$ may depend on u,t.

A local matching may change the original alias assignment. Define its
new membership $`\chi'_{i,u}`$ by decoding the preimage custom atom.
Equal weights preserve every child's total mass exactly. This choice
is made locally at each original prefix, so it remains causal.

From the full custom tuple Q, decode the original leaf from its codes,
independently of the ranks, and select each local matching index
$`c_i(Q)`$. A separate initialized register J of $`O_L(p\ell)`$
bits holds their concatenation. Execute

```math
\begin{aligned}
J&\mathrel{\oplus}=(c_1(Q),\ldots,c_p(Q)),\\
Q_i&\mathrel{\oplus}=v_{i,J_i}\quad\text{for every i in parallel},\\
J&\mathrel{\oplus}=(b_1(Q),\ldots,b_p(Q)).
\end{aligned}
\qquad(11)
```

The final table is defined on the fixed-coin image by causal decoding:
block 1's coin and weight class determine its inverse matching and
original child; that supplies the prefix for block 2, and so on.
It returns exactly the same indices, $`b_i=c_i`$. This causal decoding
is finite preprocessing for one paid table, not a serial coherent
subroutine. The first and third queries read Q and target J; the
middle queries read only disjoint short J_i and target Q_i. All native
query work returns, and J returns on the reached image.

The output has fixed coin coordinates and zero padding in every Q_i.
Thus old group, split, stage, seed, and low-code metadata are consumed
into the matched output; no additional leaf or bin label survives.
The word has a full-input unitary extension and its actual inverse.
Only reached-image clean return is asserted for J; H-zero rows are
identity on every input.

## 7. Whole coarse flags, reference work, and actual inverses

For the new local memberships, the ideal complete embedding is

```math
E|x\rangle=
\bigotimes_{i=1}^p
\left(\frac1{\sqrt{p_{i,u_i,x_i}}}
\sum_{\chi'_{i,u_i}(z)=x_i}\sqrt{w_{i,z}}|z\rangle\right),
\qquad u_i=x_1\cdots x_{i-1}.
\qquad(12)
```

Its columns have disjoint support and unit norm. At every coarse
boundary, coefficient cancellation maps a normalized logical scaling
vector to its normalized prefix-fiber state tensored with unrestricted
future fixed coins. Every original marker in band i maps to a prefix
state, a block-i vector orthogonal to $`\gamma_i`$, and unrestricted
future coins. This is the full
[compact-coin flag identity](COMPACT_COIN_FLAGS.md#3-products-of-compact-coins-for-consecutive-tree-blocks),
including all depths within a block and arbitrary superpositions.

Let $`U_i`$ be an actual fixed-coin source, with complete private
initialized source bank A_i and arbitrary dirty bank D_i, satisfying
the full initialized-isometry contract

```math
\left\|U_i\bigl(|0^{C_iA_i}\rangle\otimes I_{D_i}\bigr)
- |\gamma_i\rangle|0^{A_i}\rangle\otimes I_{D_i}\right\|_{\rm op}
\le\epsilon_i.
```

This is uniform over dirty inputs and external references. The banks
D_i of simultaneous reference sources are disjoint.

The [nearby-completion lemma](SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison)
gives an analysis-only unitary $`\widehat U_i`$ with exactly that
initialized isometry, including $`I_{D_i}`$, and
$`\|U_i-\widehat U_i\|\le6\epsilon_i`$.
No native gate implementing this completion is invoked.
Apply all actual $`U_i^\dagger`$ in parallel. The computational
last-nonzero test for a block includes its coin register **and every
qubit of A_i**. Testing the coin alone would not follow from a
one-column source guarantee. Local encoding padding remains zero on
the ideal matched image and is untouched by these reference sources.

Let $`\lambda`$ be the last nonzero augmented block, or zero if all
are zero. The ideal transforms $`\widehat U_i^\dagger`$ map the
root to all zero and band i to $`\lambda=i`$. To reproduce R87's
zero-based coarse bin, XOR

```math
\kappa(\lambda)=
\begin{cases}
0,&\lambda=0,\\
\lambda-1,&\lambda\ge1
\end{cases}
```

into the old coarse label, guarded by H. Follow this by all actual
$`U_i`$ and $`U_E^\dagger`$. The ideal middle XOR intertwines the
embedding with the desired coarse spectral XOR, with literal phase
$`+1`$. Therefore the complete physical word has error at most

```math
2\delta_E+12\sum_{i=1}^p\epsilon_i
\qquad(13)
```

on the full initialized input isometry, including every tag value,
logical column, arbitrary reference, leakage, and returned work. This
is the usual forward/middle/actual-inverse inequality; no postselection
or unproved off-anchor completion is used.

Take
$`\delta_E=\eta/512`$ and
$`\epsilon_i=\eta/[4096(p+1)]`$. Then (13) is below
$`\eta/128`$. With R87's original source, endpoint, and fine-decoder
tolerances $`\eta/128`$, the complete early error is below
$`\eta/8`$. The old coarse label is erased on the ideal source image;
actual earlier leakage propagates through the physical unitaries.
Every new component is guarded identity on H zero, and the saved
protected-zero predicate is reversed last as before.

## 8. Simultaneous native resource envelope

All per-block metadata, source rank padding, local fixed-coin work,
compression-router parameter words, and the concatenated matching
index use $`O_L(p\ell)=O_L(\ell^3)`$ conditional clean bits. The
router parameters occupy $`O_L(\ell)`$ bits per block and remain
live through route/body/unroute. The permanent Q extension, external
matching index J, reference banks A_i, short source banks, and live
router parameters are disjoint simultaneous reservations. The source
work needed for the fixed U_i has
precision $`O_L(\log(p+2))`$; its total over p parallel sources fits
this same reservation. These banks are disjoint from H, the old tags,
fixed source programs, and the free donor reserve used by the earlier
R87 fine episode. Serial episodes may reuse returned dirty pools.

Every large table reads at most $`h+C_L\ell^3`$ bits, including its
metadata, activity, and source indices. The first and last matching
tables read Q; the middle catalogue queries have only
$`O_L(\ell)`$ address bits per block. The global high-part erasure
reads h Low-bank bits plus metadata and targets separate High banks.
The joint source reads one Q bank and targets separate source fields.
The last-nonzero query reads fixed coins and their A_i and targets the
separate old coarse label. No query charges a smaller address than its
actual physical public union.

Using the charged source/query bounds, summing every simultaneous
private pool and every parallel local call, and absorbing fixed
polynomials into the exponent, gives constants independent of A_L:

```math
T_E,w_E\le C_L\left(2^{(h+C_L\ell^3)/2}+1\right),
\qquad
G_E\le C_L\left(2^{h+C_L\ell^3}+1\right),
\qquad D_E=O_L(\ell).
\qquad(14)
```

Here $`w_E`$ includes the external Low/High banks, all routing banks,
source pools, and reference-entangled dirty helpers at their actual
simultaneous peak. The compression has a constant number of global
query/router appearances; its block-local shears and routers run in
parallel. The joint rank preparation is one source. Equation (11) is
three query appearances. Fixed reference sources run on disjoint
blocks with private banks. Their precision is inverse-polynomial in p,
so their depth remains $`O_L(\ell)`$. All actual inverse appearances
are included in (14).

Fix these constants first, then choose A_L larger than the clean
reservation coefficient and the exponent coefficient in (14), with a
further fixed margin. Since h=n-M,

```math
2^{h+C_L\ell^3}
\le N\,2^{-(A_L-C_L)\ell^3+O(1)}.
```

It follows, beyond the fixed cutoff, that

```math
T_E=O_L(\sqrt N),\qquad G_E=O_L(N),\qquad
w_E=o(\sqrt N/n).
\qquad(15)
```

The independent protected margin suppresses the table volume; its full
M wires are not also appended to the public source address.

## 9. Terminal, literal widths, and scope

Use the wide construction when $`b\ge\sqrt N/n`$ and n exceeds
the chosen cutoff. Equation (15) fits the entire early dirty peak into
a fixed fraction of b. The rounded frame has the retained exact
early/terminal factorization. Apply the
[R88 terminal theorem](PREFIX_SELECTED_COIN_TERMINAL.md#1-terminal-theorem-and-scope)
to its remaining M layers, reusing the returned dirty pool, with error
$`\eta/4`$. It has the retained counts and

```math
D_{\rm terminal}
=O_L\!\left(N/b^2+\ell\log_2\log_2(M+4)\right)
=O_L\!\left(N/b^2+\ell\log_2\log_2(\ell+4)\right),
\qquad(16)
```

because $`M=O_L(\ell^3)`$. R88 invokes the prior R87 local stage,
not this theorem, so this use has no recursive proof dependency.

Restore the finite real sign gauges with their retained exact Boolean
phase queries. The prescribed complex leaf diagonal is compiled at
error $`\eta/4`$ with the retained one-clean diagonal theorem. Each
has depth $`O_L(N/b^2+\ell)`$ and the retained counts. The total
error is at most
$`\eta/32+\eta/8+\eta/4+\eta/4\lt\eta`$.
This proves the base bound (16) in the wide branch. Section 10
strengthens its terminal step and proves (1).

For $`b\lt\sqrt N/n`$, use the original-accuracy
[R72 compiler](ONE_CLEAN_DEPTH.md): its
$`O_L(N/b^2+nL)`$ depth is $`O_L(N/b^2)`$ because
$`N/b^2>n^2`$. The same compiler covers every dimension below the
fixed cutoff, absorbed in L-dependent constants. No stricter internal
tolerance changes its literal reservations. Thus the original real
width $`17(L+n+7)`$, complex width $`17(L+n+8)`$, and one external
clean qubit hold throughout.

R79 remains the absolute all-precision theorem. The base construction
already replaces serial coarse search by one causal coarse eraser.
The following prefix-selected extension improves its terminal further.

## 10. Prefix-selected recursion and iterated-logarithm depth

The local construction in Sections 2–8 extends to an unchanged external
prefix. This gives the following terminal statement, independently of
the R88 terminal used in Section 9. For $`0\le M\le n`$, the deepest
M layers of any prescribed real frame admit error below $`\eta/4`$
and

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(M+2))\right),
\qquad(17)
```

at the same literal real width and one-clean interface. For M zero use
identity. Taking M=n proves the real version of (1); adding the retained
complex leaf diagonal at error $`\eta/4`$ proves its complex version.
Alternatively, (17) can replace the terminal in Section 9. The extension
uses only the local causal stage, the retained sharp-query interface,
and the retained constant-height bottom. It does not invoke the whole
statement (1) in its proof.

### 10.1 A local stage and catalogues selected by an immutable prefix

When m layers remain, let P be the unchanged $`n-m`$-bit prefix.
Choose r bottom protected bits, put $`h=m-r`$, and use the exact
prefix-selected factorization

```math
W_m=W_r E_m,\qquad
E_m=\sum_P|P\rangle\langle P|\otimes
\left(W_{P,h}\otimes\Pi_Z+I\otimes(I-\Pi_Z)\right).
```

Save $`H=[Z=0^r]`$ in the sole external clean qubit. For stage accuracy
$`2^{-l}`$, define

```math
K=m+l+2,\qquad u=l+\lceil\log_2K\rceil.
\qquad(18)
```

Round only this stage's early pairs at
$`\tau=\lceil2\log_2(m/2^{-l})\rceil+12`$, with error at most
$`2^{-l}/32`$. Every stage-local source and coin uses that tuple.
The unprocessed $`W_r`$ retains its previous prescribed values until
its own stage. Unitary telescoping, rather than a common denominator
across recursion stages, controls the error.

The fine bins have $`p\le C\tau\log_2(m+2)\le C u^2`$.
Each compact block adds $`O(\log K)`$ custom coordinates and short
catalogue fields. Its fixed reference source needs
$`O(l+\log(p+1))`$ initialized work. All code padding, fine/source
records, matcher indices, router parameters, reference banks, and free
donor reserve consequently fit $`C u^3`$ additional protected bits,
after fixing the native constants. No such logical bit is also counted
as an arbitrary dirty query helper.

Choose catalogues **separately for each fixed P**, common across all
internal original prefixes. There are only
$`2^m\mathrm{poly}(K)`$ internal prefix/group families. Thus the
shared hash and matching catalogues have $`\mathrm{poly}(K)`$
entries, with $`O(\log K)`$-bit indices, uniformly in P. Classical
preprocessing may enumerate all P sectors; its efficiency is not an
assumption. The catalogue vectors and matrices may depend on P:
a shear reads P, its local High bank, and its short selector; a shift
reads P and its short selector. Their targets are separate local data
banks. P is never modified, routed, or tested for zero.

These simultaneous calls use the
[R88 read-only sharp-query interface](PREFIX_SELECTED_COIN_TERMINAL.md#4-parallel-sharp-queries-with-a-read-only-public-address),
with private dirty pools and output controls. Shared P enters only as
Clifford controls at non-Clifford boundaries. The full causal matching
and rank-source argument holds in each orthogonal P sector, so errors
take a maximum over P. Requiring one catalogue common across all global
P sectors would unnecessarily introduce global-n-sized labels; that
is not this construction.

### 10.2 One-stage addresses, accuracy, and resource bounds

The physical local Q bank has h data bits plus $`O(p\log K)`$
custom coordinates. Including P and every metadata/source field, each
large table reads at most

```math
n-r+C u^3
\qquad(19)
```

bits. This covers metadata load, global High erasure, the single joint
rank source, both matching-selector queries, the augmented-block
last-nonzero query, all anchors, and their actual inverses. Local
shears and catalogue shifts have smaller addresses. Fine recovery has
the retained $`n-r+O(u)`$ address bound.

An arbitrary table on all $`O(u^3)`$ metadata bits is charged within
(19); it is not described as polynomial in K. Only local
$`O(\log K)`$-bit masks/catalogue fields give polynomial-size short
tables. The number of parallel one-bit outputs and instances is at most
$`K^{d_0}`$ for a fixed $`d_0`$, with no polynomial multiplier in
global n.
The source changes only local Q and its local work; none of these
operations produces an n-bit output merely because it reads P.

Assign each constant-number nonreference source error at most
$`2^{-l}/1024`$, and each reference source error at most
$`2^{-l}/[4096(p+1)]`$. Their precisions are $`O(u)`$.
Equations (12)–(13), the retained fine-source hybrid, and the rounding
error give total stage error below $`2^{-l}/4`$. All source inverses
are actual physical inverses. H-zero rows and complete occupied-work
identity hold as in the base stage.

The constants here are uniform in the varying stage precision l;
we do not substitute l into the fixed-L asymptotic statement (14).
They follow directly from the components. Namely,
$`\tau=O(u)`$, $`p=O(u^2)`$, each custom block adds
$`O(\log K)`$ code coordinates, and the shared hash/matching
catalogues have $`\mathrm{poly}(K)`$ size with one fixed degree.
Each reference bank has $`O(l+\log(p+1))`$ bits. The normalized
source grammar has a constant number of large Boolean-query positions;
its precision-dependent row factors, including factors of size
$`4^{O(l)}`$, are included in the address allowance $`C u^3`$.
They are not mislabeled as polynomial in K. All output lengths, local
mask banks, routing arrays, coefficient programs, and simultaneous
instance counts are polynomial in K with fixed degree, since
$`l\le K`$ and $`u=O(K)`$. The sharp query contributes no
polynomial in global n: that address length enters only its logarithmic
depth and minimum pool width. This fixes absolute C and an instance/row degree $`d_0`$ before
the reserve choice, uniformly over l,m,P and the angle tuple.

Put $`\overline Q=2^{n-r+C u^3}`$. There are at most
$`K^{d_0}`$ simultaneous instances, each with at most
$`K^{d_0}\overline Q`$ rows. Give each one a private fraction
$`b/(C'K^{d_0})`$ of the dirty pool. A distinct fixed envelope
degree $`d\ge3d_0`$, increased once for standalone local terms,
then gives

```math
\begin{aligned}
T_{\rm stage}&\le C K^d(\sqrt{\overline Q}+\overline Q/b+1),\\
G_{\rm stage}&\le C K^d(\overline Q+1),\\
D_{\rm stage}&\le C\left[\log_2(n+C u^3+2)+u+
                         K^d\overline Q/b^2\right].
\end{aligned}
\qquad(20)
```

The coefficient and degree depend only on fixed native constructions,
not on the protected size chosen next. The logarithm is
$`O(\ell)`$ at each nonfinal stage below. Local interval routing has
$`O(\log(m+2))`$ depth and polynomial local count/width. The source
has a fixed number of query positions, while its reference sources
run in parallel; neither hash stages nor causal inverse decoding are
serial physical loops.

### 10.3 Recursive reserve, summability, and stage count

Start at $`m_0=M`$, and at a nonfinal stage choose

```math
\begin{aligned}
l_i&=L+\lceil\log_2(m_i+2)\rceil+10,\\
K_i&=m_i+l_i+2,\qquad u_i=l_i+\lceil\log_2K_i\rceil,\\
r_i&=\lceil A u_i^3\rceil,\qquad m_{i+1}=r_i.
\end{aligned}
\qquad(21)
```

Fix A after the constants in (20) and the clean reservation. Since
$`u_i\ge\log_2K_i`$, it can ensure

```math
r_i-Cu_i^3\ge(2d+80)u_i,
```

and hence

```math
K_i^d2^{-(r_i-Cu_i^3)/2}\le K_i^{-40},\qquad
K_i^d2^{-(r_i-Cu_i^3)}\le K_i^{-80}.
\qquad(22)
```

For fixed L the nondecreasing function $`r(m)`$ is
$`O_L(\log_2^3(m+2))`$. Choose a fixed stopping height
$`m_*=m_*(L,A)`$ so that for all $`m>m_*`$,

```math
r(m)\le m/2,\qquad r(r(m))\le\log_2 m.
\qquad(23)
```

Both inequalities eventually hold. Stop on entering $`m\le m_*`$
and use the retained
[smaller-error constant-height bottom](PREFIX_SELECTED_COIN_TERMINAL.md#6-recursive-banks-errors-and-resource-sums).
Its error is below $`15\eta/64`$. This is the existing terminal
interface, not a substitution of stricter accuracy into the original
literal-width formula.

Read backwards, nonfinal heights grow at least geometrically.
Therefore the inverse-power sums in (22) are bounded. The large-table
terms in (20) sum to
$`O_L(\sqrt N+N/b)`$, $`O_L(N)`$, and $`O_L(N/b^2)`$ respectively.
Also

```math
\sum_iK_i^d\le C_L(M+2)^d\le C_L(n+2)^d=O_L(\sqrt N)
```

beyond a fixed cutoff, so the standalone gate terms do not incur an
unpaid stage multiplier. The same polynomial resources fit the wide
dirty pool.

Every two nonfinal steps contract the height at least as much as one
binary logarithm. By monotonicity, while those stages exist,
$`m_{2j}\le\log_2^{(j)}M`$. The stage count is therefore
$`O_L(1+\log_2^*(M+2))`$. Since $`u_i=O_L(\ell)`$, the serial
non-table depths in (20) sum to the depth in (17). The constant-height
bottom adds $`O_L(N/b^2+\ell)`$.

The stage errors obey

```math
2^{-l_i}/4\le\eta/[4096(m_i+2)],
\qquad \sum_i2^{-l_i}/4\le\eta/2048.
```

For general signed or singular input, first make the same finite
whole-family magnitude/sign reduction as R88, at error
$`\eta/512`$, and pay its two exact Boolean sign phases once. The
resulting positive tuple is then factored and rounded stage by stage.
Including the smaller-error bottom, the total remains below
$`\eta/4`$. The bottom uses the remaining layers of that positive
tuple, rather than returning to the pre-gauge signed coefficients.

### 10.4 Literal widths and all terminal heights

On the sufficiently large branch $`b\ge\sqrt N/n`$, every number
of simultaneous private instances and every standalone routing/mask
bank is bounded by a fixed polynomial in n. The public address length
is $`O_L(n)`$. Exponential dominance ensures each allocated private
pool exceeds the sharp-query minimum width, uniformly over all
$`M\le n`$ and all tuples. Its addresses and targets are excluded
from the pool. The logical clean reserve fits by the stopping rule.
Returned pools are reused only between completed positions or stages.

For $`b\lt\sqrt N/n`$, use the retained
[ungrouped smaller-error terminal](PREFIX_SELECTED_COIN_TERMINAL.md#7-literal-width-threshold-and-finite-dimensions)
at the original prescribed angles. Its reservation is
$`B'_0=L+n+9`$, so $`b\ge17(L+n+7)`$ leaves
$`b-B'_0\ge16b/17-2\ge b/2`$ for its dirty query pool.
Its error is below $`\eta/4`$, its counts are retained, and

```math
D_T=O_L(N/b^2+nM+M^2)=O_L(N/b^2).
```

Use the same fallback below the fixed dimension cutoff, uniformly in
all M and arbitrarily large b. Its remaining bounded depth is absorbed
by $`O_L(\ell(1+\log_2^*(M+2)))`$. Thus (17) holds at the stated
literal width, and the final complete-frame theorem (1) follows.

This refinement still has an unbounded iterated-logarithm stage count.
An $`O_L(\ell)`$ general compiler and optimal depth remain open.

## 11. Bounded evidence

The [causal compact-coin tests](../../tests/test_causal_compact_coin_depth.py)
exercise an actual batched compressor with collisions and staged
assignment, its inverse and metadata return, prefix-dependent local
partitions, the common-catalogue global selector word, every marker
column, and the full tag word on a coherent reference input. A separate
case shows why the reference source bank must be included in the
last-nonzero block. These finite tests use explicit common shift unions;
they do not certify the asymptotic catalogue size, native resource
bounds, or recursive stage count.

The [subset-routing tests](../../tests/test_short_label_subset_routing.py)
and [compact-coin tests](../../tests/test_compact_coin_flags.py) check
their separate finite primitives. All asymptotic, initialized-work,
source-precision, and full-input inactive contracts are proved above
and in the cited native lemmas. No scalable circuit emitter is claimed.
