# Local frame depth when dirty workspace is inexpensive

[Coarse flags and sparse correction](COARSE_FLAG_SPARSE_CORRECTION.md) ·
[Small-workspace descendant sources](SMALL_WORKSPACE_DESCENDANT_SOURCE.md) ·
[Portable payloads](PORTABLE_WORKSPACE_CATALOGUE.md) ·
[Exact compact geometric source](../../docs/FAULT_TOLERANT_COMPILER.md#3-an-exact-geometric-source-of-logarithmic-width)

**Local component of [R95](LOGARITHMIC_PRECISION_DEPTH.md).** This chapter uses
a deliberately generous local count and dirty-space budget. It does not
assert locally optimal counts. The complete-frame application places its
physical word, including all its work, inside a portable catalogue slot.

## 1. Local contract

There are absolute constants C and a fixed polynomial P such that every
prescribed real h-bit Hopf frame W, with integers $`h\ge1`$ and
$`\kappa\ge6`$, has a physical Clifford+T word U with one initialized
qubit d and an arbitrary dirty pool D satisfying

```math
\|UJ_d-J_d(W\otimes I_D)\|_{\rm op}\le2^{-\kappa},
\qquad D_T(U)\le C(h+\kappa),
```

```math
T(U),\ G(U),\ |D|\le2^{P(h+\kappa+1)}.                 \qquad(1)
```

Unused supplied wires are spectators. All columns, signs, singular pairs,
literal phases, clean leakage and arbitrary dirty/reference inputs are
included. The word and its actual inverse exist by terminating certified
preprocessing; no efficient preprocessing or scalable emitter is claimed.

For bookkeeping put $`k=\max\{h,\kappa,6\}+c_*`$. Individual
tolerances below are fixed constant multiples of $`2^{-k}`$. Prove the
component bounds uniformly in the internal precision parameter k, then
fix A, the number of stages and the final group bound. Only afterward
choose the fixed integer $`c_*`$ to cover their error factors and enlarge
the absolute cutoff if needed. Thus $`k=O(h+\kappa)`$.
There is no invocation of the theorem being proved on its original
dimension. Section 5 gives an explicit two-stage dependency.

## 2. An early stage with a logarithmic precision register

Apply the algebra of the
[coarse-flag construction](COARSE_FLAG_SPARSE_CORRECTION.md#2-a-coarse-exact-tree-and-a-small-positive-residual-norm)
to a single m-bit frame, with no external prefix. Its coarse mathematical
tree has dyadic probability precision
$`\tau_0=O(\log_2(m+2))`$, independent of the final tolerance. The
sparse difference has absolute matrix norm below $`1/16`$. A constant
normalization overlap combines the coarse frame and this difference;
constant-length oblivious amplification recovers the desired frame.

All ordinary real sources, including the initial tagged source, joint
rank source, endpoint sources, and sparse row/column sources, now use the
[small-workspace source](SMALL_WORKSPACE_DESCENDANT_SOURCE.md). At tolerance
$`2^{-(k+O(1))}`$, its only initialized source fields are
$`O(\log_2(k+2))`$ coefficient-index bits and one attenuation bit.
The precision core and borrowed real-rotation signals are arbitrary dirty
work. Its source depth is $`O(k+\log_2(m+k+2))`$. The constant
number of source positions does not multiply this bound by the number of
coarse bins. For transported anchors and the sparse overlap, use the
source lemma's nearby completion on the ENTIRE physical initialized
input space, including arbitrary dirty work. Its ideal source isometry
tensors identity on that work. The actual native source need not return
its dirty precision core exactly, and no exact tensor-identity action
on its unused columns is assumed.

For the joint rank source retain R89's
[metadata-selected window extension](CAUSAL_COMPACT_COIN_DEPTH.md#5-one-joint-source-on-metadata-selected-rank-windows):
only its preserved labels, window masks and anchor predicates change.
The same logarithmic source-scratch bound is independent of the number
of labels; their addresses and native table sizes remain charged below.

The simultaneous reference sources are instead implemented EXACTLY as
described in Section 3. Their private initialized registers have size
$`O(\log_2(m+2))`$ per bin, with no final-precision term.
The geometric bins number $`O(\log_2^2(m+2))`$. Consequently all
live code, hash, matching, router, tag, failure and reference fields,
together with the separate coarse/correction banks and LCU branch bit,
fit within

```math
r(m,k)=\left\lceil A\left[
 \log_2(k+2)+\log_2^3(m+2)\right]\right\rceil            \qquad(2)
```

original bottom bits, after fixing a sufficiently large absolute A.
All constants and table degrees are fixed before A is chosen. Source
labels and tags are counted separately from source scratch. The two LCU
branches do not share simultaneous initialized banks.

Use the early stage only when $`r(m,k)\lt m/2`$. Save its bottom-zero
predicate into its single initialized external bit, use the original
bottom bits as conditional work, and invert the actual predicate word
last. Every complete body is exactly identity on the inactive occupied
sector. The source-isometry and full-unitary comparisons, all-work
reflection and actual-inverse telescope from the coarse-flag proof apply
unchanged. With fixed error margins, this stage implements the first
$`m-r(m,k)`$ layers with error $`O(2^{-k})`$.

With unrestricted local dirty allocations, each finite XOR or Boolean
phase table is implemented by R91 with a private pool at least a fixed
constant times the square root of its row/output product. Its T-depth is
logarithmic in the address length. Every address and output length is
bounded by a fixed polynomial in $`m+k+1`$, and there are only a
fixed-polynomial number of simultaneous local instances. All source
loaders, elementary Clifford gates, pools and actual inverses are
therefore included in the safe bounds

```math
D_{T,\rm early}=O(k+m),\qquad
T_{\rm early},G_{\rm early},w_{\rm early}
 \le2^{P_0(m+k+1)}                                    \qquad(3)
```

for one fixed polynomial $`P_0`$. This is a local statement. It does not
replace the separately proved global count ledger of the coarse-flag
construction by an unspecified large table.

## 3. Exact reference coins without a precision bank

For a bin of height $`q_i`$, put $`s_i=q_i\tau_0`$. Enlarge this
geometric exponent to

```math
s_i'=2^{\lceil\log_2(s_i+1)\rceil}-1,
\qquad s_i\le s_i'\le2s_i+1.                           \qquad(4)
```

Each original alias acceptance probability with denominator
$`2^{s_i}`$ also has denominator $`2^{s_i'}`$. Use that same
acceptance probability in the finite two-outcome alias construction.
Its fixed coin is the tensor product of the uniform $`2^{q_i}`$-column
state and the geometric state with probabilities

```math
\frac12,\frac14,\ldots,2^{-s_i'},2^{-s_i'}.
```

It has $`2^{q_i}(s_i'+1)`$ positive atoms, and still uses
$`q_i+O(\log_2(m+2))`$ coin bits. Each geometric weight class still
has at least $`2^{q_i}`$ atoms. The rank-window inequalities and
inverse-polynomial matching densities of R89 therefore retain their
same orders and absolute constants after enlargement. The causal matching
and hash catalogues are rebuilt for these enlarged exact coins; an old
matching is not silently reused.

The geometric part has the power-of-two dimension $`M_i=s_i'+1`$.
Apply the retained exact capped geometric source, whose recursive
Gray-path implementation has

```math
T_i,G_i=O(M_i),\qquad
w_{i,\rm clean}=O(\log_2 M_i),\qquad D_{T,i}\le T_i.
```

This is the complete recursive source, not a serial list of independently
synthesized multi-controlled Hadamards. Its private scratch returns
exactly on every coin input, and its actual inverse is used.

The geometric-bin rule is

```math
q_i\le\frac{a-d_i+m_{\rm free}}{2(\tau_0+2)},
\qquad a=m-r,
```

where $`d_i\ge0`$ is the bin endpoint and the free donor reserve is
contained in the protected bottom bank. Thus $`m_{\rm free}\le r`$
and $`s_i\le m/2`$. Equation (4) gives
$`\max_i M_i=O(m)`$. Every reference source has disjoint coin and
private registers, so the tensor reference word has T-depth $`O(m)`$.
The total private initialized work is

```math
\sum_iO(\log_2 M_i)=O(\log_2^3(m+2)).                 \qquad(5)
```

The uniform-column Hadamards are Clifford. No final-accuracy synthesis
appears in these fixed reference preparations.

Include every private reference scratch bit in its augmented
last-nonzero-block test. Each exact reference inverse maps the root coin
to its initialized column and maps its orthogonal complement away from
that column. R89's full-band argument therefore applies literally,
including the zero padding. Reference words may run unconditionally:
their complete inverse/controlled-label/forward sandwich is identity on
the entire inactive sector because the middle label XOR is identity
there. Their private banks can be arbitrarily occupied on that sector.
This avoids assuming independent initialized flags for simultaneously
enabled sources or omitting a shared-control depth charge.

## 4. Selecting complete local words by the actual prefix

A local tail may be selected by an unchanged prefix P. Local count
optimality is unnecessary, so use one slot for each ACTUAL prefix value,
rather than a universal catalogue of all rounded frame codes.
Reserve the bottom three original logical bits $`(c,d,G)`$, save their
all-zero predicate in the one external clean bit H, and apply
CNOT from H to c. For each prefix p compile the controlled Hopf frame
$`A_p=\mathrm{diag}(I,W_p)`$ on $`(c,X)`$ with its own
physical one-clean word $`U_p`$, initialized d and arbitrary pool
$`D_*`$. The complete payload is $`S=(c,X,d,D_*)`$.

The selected index is a distinct inactive code zero unless H and the
returned group guard G are both one. Active codes encode P with one
additional tag bit. For an a-bit prefix there are at most $`2^{a+1}`$
slots, and the index has at most $`a+1`$ bits. Code zero is the actual
empty native word. Every other padded code can also be identity.

Let $`F=\bigotimes_pU_p`$, and let $`S_p`$ swap the whole payload
with slot p. With rightmost factors acting first,

```math
K_p=S_pFS_pF^\dagger
 =U_p(S)\otimes U_p^\dagger(\mathrm{slot}_p),
```

```math
K_pX_cK_p^\dagger X_c
 =(U_pX_cU_p^\dagger X_c)(S)\otimes I_{\rm array}.      \qquad(6)
```

Use the exact R91 paired-array router and its complete returned-work
selected swap at each occurrence. Its selected index is obtained by a
paid XOR query into an arbitrary program; no initialized program is
assumed. If a slot has width w and V slots are allocated, this costs
$`O(wV^2)`$ gates and arbitrary dirty helpers, and logarithmic depth
in the index length. All four F appearances run their disjoint physical
slot words in parallel. Therefore if each slot has depth $`d_0`$ and
counts/width at most $`2^{P_1(h+k+1)}`$, the wrapper has

```math
D_{T,\rm wrap}=O(d_0+\log_2(h+k+2)),\qquad
T_{\rm wrap},G_{\rm wrap},w_{\rm wrap}
 \le2^{P_2(h+k+1)}.                                   \qquad(7)
```

All live arrays, query pools, predicates and payloads are disjoint.
The original local frame height bounds a. A fixed number of such
wrappers preserves a fixed-polynomial exponent in (7).

The R94 inverse-isometry estimate gives error at most twice the entry
error for (6) on c one and d zero. On c one the ideal commutator is
$`W_p`$. Direct sums over P take a maximum, not a sum over slots.
The group guard excludes c,d,X and every travelling dirty bit. It tests
only H and unchanged lower logical bits; its actual inverse restores
the original G value exactly. The full word is identity on H zero,
including all occupied slots and payloads. The ideal-embedding telescope
propagates c/d leakage through subsequent full unitaries. Finally undo
the c preparation and original H predicate, with their actual inverses.
These are precisely R94's complete-payload and guard interfaces.

## 5. A two-stage construction, followed by constant-size height groups

First suppose k is above one absolute cutoff. Starting at $`m=h\le k`$,
apply Section 2 if its protected reserve is less than m/2. Otherwise
omit that early stage. In either case the remaining height $`m_1`$
satisfies

```math
m_1\le C_1\log_2^3(k+2).                              \qquad(8)
```

When the stage runs, this follows from (2). When it is omitted,
$`m\le2r(m,k)`$ and $`m\le k`$ give the same inequality.

If $`m_1\le3`$, use the fixed-dimensional bottom construction below
and omit this wrapper. Otherwise apply the wrapper of Section 4 to the
WHOLE remaining $`m_1-3`$-bit
frame above its bottom three bits. Every controlled slot entry is a
Hopf frame on at most $`m_1-2`$ logical bits. Construct that entry by
one more Section 2 opening, now at its OWN local height and with no
external prefix. If its reserve is too large, omit that opening too.
Its residual height $`m_2`$ obeys

```math
m_2\le C_2\left[\log_2(k+2)
                  +\log_2^3(m_1+3)\right]
     \le C_3\log_2(k+2).                              \qquad(9)
```

The last inequality holds beyond one absolute cutoff because
$`(\log\log k)^3=o(\log k)`$. The failed reserve test gives
the same inequality with a factor two. This second opening is performed
inside independent disjoint prefix slots: its hash and matching
catalogues depend on its OWN logical height, not on the original long
prefix. No prefix-uniform metadata estimate is being substituted into
(9).

It remains to compile frames of height at most the right side of (9).
Heights at most three use the fixed-dimensional bottom construction
directly. At larger heights reserve bottom three bits and partition the
rest into groups of height
at most

```math
g=\max\{1,\lfloor\tfrac14\log_2(k+2)\rfloor\}.
```

There are at most a fixed absolute number $`R_*`$ of groups, including
a possible short first group. For every actual prefix and every group,
compile its controlled frame of height at most g+1 using R79 at
precision k. For sufficiently large k,

```math
k>\frac{2^{g+1}}{(g+1)^2}.
```

Thus the already proved high-precision clause gives entry depth O(k).
Take its dirty pool at least
$`17(k+g+8)`$ and a sufficiently large constant times
$`\sqrt{2^{g+1}k}`$, so its width penalty is bounded.
The entry count and width are polynomial in k. Use Section 4's full
payload wrapper for each group. The actual prefix length is at most
$`m_2`$; its catalogue, router and query costs satisfy (7), while
every group still has depth O(k).

The deepest three bits are a prefix-selected bank of fixed-dimensional
real frames. They can be compiled by a fixed-length decomposition into
addressed real rotations. The zero-clean borrowed-signal primitive,
with R91 exact mask replacements, gives depth O(k+log(h+2)), arbitrary
dirty work and error $`O(2^{-k})`$. This is also the bottom-three
compiler after the first wrapper. It requires no additional initialized
flag. The saved H predicate has already been inverted before this
bottom factor runs.

This is an explicit finite construction: first early opening, one
parallel bank of second-opening entries, and one parallel bank of
R79 high-precision entries inside each of a fixed number of bottom
groups. There is no appeal to the new lemma at the same dimension or
precision, and no use of R94 beyond its local precision range.

If k is below the fixed cutoff, then h and k are both bounded. The
retained one-clean layerwise compiler at precision k supplies
(1) after increasing the absolute constants. Its actual word can be
used without a finite-net or precision-dependent-cutoff argument.

## 6. Full error, counts and the intended use

All early stages, slot words and fixed-dimensional bottoms have the
full initialized-isometry or full-unitary guarantees stated above.
The number of sequential early stages and bottom groups is bounded by
one absolute constant. Parallel prefix banks contribute a maximum
error. Each complete transfer commutator contributes a factor at most
two, and the nesting depth of these commutators is fixed. Choose the
single constant $`c_*`$ after these finite factors so that their total
is at most $`2^{-\kappa}`$. All inverse appearances use actual
native words. The hybrid compares each bracket on its ideal initialized
input and propagates earlier leakage unitarily; no bank is reset.

Equations (3), (7), (8) and (9) give depth O(k) on this same circuit.
A fixed number of products, maxima, catalogue multiplicities
$`2^{O(h)}`$ and squared-array costs applied to fixed-polynomial
exponential bounds remains $`2^{P(h+k+1)}`$ for one fixed P.
The exponent and all cutoffs are independent of the frame and of
precision. Thus all Clifford gates and the simultaneous arbitrary
dirty reservation satisfy (1), completing the local lemma.

For the outer R94 application, take a group of height
$`h\le\ell/8`$ and precision $`\kappa=O(\ell)`$. Equation (1)
has logarithmic count/width bound polynomial in ell. Together with
R94's rounded code length $`q=O(\kappa2^h)=o(n)`$,

```math
\log_2\!\left(w_{\rm loc}2^{2q}
 +2^q(t_{\rm loc}+g_{\rm loc})+1\right)=o(n).
```

The complete local word and all its work therefore fit the original
outer catalogue absorption test. Its depth is now $`O(h+\kappa)`$.
This local lemma alone does not improve the global R94 bound: the
separate opening improvement and its simultaneous global ledger must
also be established.
