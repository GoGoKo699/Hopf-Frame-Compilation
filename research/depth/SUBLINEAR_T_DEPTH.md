# Sublinear T-depth at square-root dirty width

[Resource statement](../../docs/OPEN_PROBLEM.md) · [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) · [Coherent descendant loader](COHERENT_DESCENDANT_LOADER.md) · [Linear-size dirty indicators](NONUNIFORM_DIRTY_INDICATOR.md)

This construction bypasses the successive early Hopf layers by encoding
their complete columns into orthogonal depth-and-prefix sectors. A
coherent real-state loader implements that encoding; a Clifford row
isometry gives normalization $`O(\sqrt n)`$. Amplitude amplification
recovers the whole early frame. The last polylogarithmically many original
layers use the retained native layerwise compiler.

The construction uses the Boolean-phase state-preparation method of
Gosset–Kothari–Wu, with a separately proved terminating selection rule,
coherent descendant labels, and arbitrary-dirty query implementations.
It is not an invocation of an uncharged state-preparation oracle.

## 1. Fixed-accuracy theorem

For every fixed $`0\lt\eta\le1/64`$, there are constants
$`C_\eta,c_\eta`$ such that every prescribed complete real Hopf frame W
on $`n\ge1`$ qubits, $`N=2^n`$, has a coherent Clifford+T circuit V
with two external initialized flags and at most
$`b=\lceil C_\eta\sqrt N\rceil`$ arbitrary dirty qubits satisfying

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T\le c_\eta\sqrt N,\qquad
G\le c_\eta N,\qquad
D_T\le c_\eta\sqrt n\,[\log_2(n+2)]^2=o(n).
```

All bounds hold on the same circuit. G counts elementary Clifford gates.
T-depth permits arbitrary intervening Clifford circuits, whose count is
still charged. The norm includes every logical column, arbitrary dirty
inputs and reference correlations, literal phases, and all clean and
logical-work leakage. Unused supplied dirty wires remain untouched.

Angles admit terminating certified evaluation, including unequal and
singular tuples. Classical preprocessing is terminating but has no
claimed efficiency bound. Every source, query, actual inverse, selective
phase and recovery call is compiled and counted. There are no measurements,
resets, supplied precision states, QRAM, or uncharged coherent evaluators.

This answers the fixed-accuracy $`o(n)`$ existence question positively.
It does not establish optimal T-depth, a precision-uniform refinement of
Result D, or a general native emitter. The unrestricted growing-depth
lower-bound problem remains open.

## 2. Protect one terminal logical block

For all sufficiently large n, choose

```math
M=\lceil[\log_2(n+2)]^2\rceil,\qquad
h=n-M\ge1,\qquad Q=2^h=N2^{-M}.
```

Let X be the first h logical bits and Z the final M bits. The first h
addressed Hopf layers have the exact complete-operator product

```math
W_E=W_h\otimes\Pi_Z+I_X\otimes(I_Z-\Pi_Z),
\qquad \Pi_Z=|0^M\rangle\langle0^M|.
```

Here $`W_h`$ is the complete h-bit real frame using precisely the first h
levels of the original angle tree. Every early layer requires its full
lower suffix to be zero; on Z zero its remaining suffix condition is
exactly that of $`W_h`$. This proves the identity on every X input.
The final M original layers give $`W=W_TW_E`$.

Compute $`H=[Z=0^M]`$ into one initialized external flag by the retained
[exact logarithmic-depth predicate](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
Its returned dirty helpers are separate from live source/query registers.
Keep H unchanged until the early circuit ends. On H equal to one, Z
provides $`O_\eta(\log(n+2))`$ initialized working bits. On H equal to
zero, Z is arbitrary. Every completed word below is literally identity
on that sector, on the entire Hilbert space. Z is never retested while it
holds source work. The other external flag may remain unused until the tail.

The dirty query pool is disjoint from X, Z, and both external flags.
All Z reservations are simultaneous: tag, LCU index, inner attenuation,
reused angle program, and outer attenuation. Their sum is
$`O_\eta(\log(n+2))`$, which fits M for sufficiently large n.

## 3. Compress the complete columns

Set

```math
s=\lceil\log_2(h+1)\rceil,\qquad K=2^s<2(h+1).
```

The s-bit tag D lies in Z. Give the state column j equal to zero tag zero.
For a nonzero marker $`j=p1\,0^{h-d-1}`$ at depth d, give it tag
$`t(j)=d+1`$. The marker convention is the literal
[Hopf interface](../../docs/HOPF_INTERFACE.md#2-tree-and-marker-convention).

Define a source isometry R by

```math
R|j\rangle
=\sum_i(W_h)_{ij}|i\rangle_X|t(j)\rangle_D|0\rangle_{\rm other}.
```

At a fixed nonzero tag, columns have distinct prefixes p and disjoint
descendant supports. Different tags are orthogonal; the state column has
its own tag. Every column is normalized. Hence R is an isometry, and the
same remains true after replacing each descendant vector by any normalized
vector within its own tag-and-prefix sector.

Reach the source anchors using two reversible XOR operations:

```math
D\mathrel{\oplus}=H\,t(X),\qquad
X\mathrel{\oplus}=H\,m(D),
```

where $`m(0)=0`$, $`m(d+1)=2^{h-d-1}`$, and unused tags have zero mask.
On initialized D and active H this sends a marker to
$`|p0^{h-d}\rangle_X|d+1\rangle_D`$, and sends zero to the root anchor.
Both operations are permutations on their full input spaces; their actual
reverse is the decoder. The first has Q rows and s output bits. The second
has K rows and h output bits and may use the exact dirty-mask CNOT echo.
Neither requires a second initialized h-bit data register.

Apply the [coherent descendant loader](COHERENT_DESCENDANT_LOADER.md) to
each anchor. It prepares the signed root or marker vector in the zero
suffix, preserving its tag and upper prefix. Its derived Boolean tables
are indexed by tag, current X, and the small LCU index. The current X
already contains both the unchanged prefix and the varying suffix, so
the table size is $`Q\,\mathrm{poly}_\eta(n)`$, not $`Q^2`$.

Let $`\widehat R`$ denote the entire native compression-and-loader word,
and J initialize its tag and all source work. For a chosen tolerance
$`\delta`$, including all synthesized rotations and source leakage,

```math
\|\widehat R J-R\|\le\delta.
```

This is an operator-norm source bound. The error vectors for distinct
anchors remain in orthogonal tag-and-prefix sectors; thus their norm is
the maximum branch error. Every arbitrary dirty query/mask wire returns
exactly on all inputs. The source bound consequently tensors with dirty
systems and references. It makes no assertion that an occupied angle
program represents the intended target rotation.

Define the row isometry

```math
P|i\rangle=|i\rangle_XH^{\otimes s}|0^s\rangle_D
|0\rangle_{\rm other}.
```

Since every used tag has the same Hadamard overlap,

```math
P^\dagger R=W_h/\sqrt K.
```

The unused padded tags are harmless. No target-dependent row source,
Perron weights, or second long clean index is needed.

## 4. From a source error to an honest unitary comparison

The following dimension-independent lemma prevents a source promise from
being silently substituted for full operator accuracy.

Let A and B be isometries with $`\|A-B\|\le\delta\lt1/3`$. Put

```math
F=BA^\dagger+(I-BB^\dagger)(I-AA^\dagger).
```

Then $`FA=B`$, $`F^\dagger B=A`$, and

```math
\|F-I\|\le\|B-A\|+\|BB^\dagger-AA^\dagger\|\le3\delta.
```

F is invertible. Its polar unitary
$`C=F(F^\dagger F)^{-1/2}`$ obeys $`CA=B`$, because
$`F^\dagger FA=A`$. All singular values of F are within
$`3\delta`$ of one, giving

```math
\|C-I\|\le\|C-F\|+\|F-I\|\le6\delta.
```

Apply the lemma with $`A=\widehat RJ`$ and B equal to R.
The mathematical completion $`R_* = C\widehat R`$ satisfies

```math
R_*J=R,\qquad \|R_* -\widehat R\|\le6\delta.
```

C is used only in this proof and is never implemented. It may act on
source scratch that contains synthesis leakage. The comparison acts on
the mode space alone and is tensored with identity on exactly returned
dirty work. Extend it by identity on H equal to zero and on unused Z bits.
Actual inverse words obey the same full norm estimate. No intermediate
source work is assumed to reset.

## 5. Recover the full early unitary with its literal phase

Choose

```math
k=\lceil\pi\sqrt K/4\rceil,\qquad
\vartheta=\pi/(4k+2),\qquad
\gamma=\sqrt K\sin\vartheta\le1.
```

The inequality follows from $`\vartheta\le1/\sqrt K`$ and
$`\sin\vartheta\le\vartheta`$. On one extra zero bit z inside Z,
use the literal real rotation

```math
G_z|0\rangle=\gamma|0\rangle+\sqrt{1-\gamma^2}|1\rangle.
```

For the ideal comparison unitary
$`U_*=G_z\otimes[(H_D^{\otimes s}\otimes I)R_*]`$, enlarge J to
initialize z as well. Its accepted block is

```math
J^\dagger U_*J=\sin\vartheta\,W_h.
```

Write $`\Pi=JJ^\dagger`$ and $`S=I-2\Pi`$. The iterate

```math
A_*=-U_*SU_*^\dagger S
```

satisfies $`A_*^kU_*J=JW_h`$. Indeed, for every orthonormal input
column, $`U_*J`$ has good component $`\sin\vartheta\,JW_h`$
and a normalized bad component of norm $`\cos\vartheta`$. These
good/bad pairs are mutually orthogonal because $`W_h`$ is unitary.
On each pair the iterate has matrix

```math
\begin{pmatrix}
\cos2\vartheta&\sin2\vartheta\\
-\sin2\vartheta&\cos2\vartheta
\end{pmatrix}.
```

It sends the initial coordinates to
$`(\sin((2k+1)\vartheta),\cos((2k+1)\vartheta))=(1,0)`$.
This is an all-column unitary recovery, with positive literal phase.

Implement the source, row Hadamards and damping as H-guarded native words.
Use the exact selective phase $`S_H=I-2[H=1]\Pi`$ and the physical
iterate

```math
\widehat A_H=Z_H\widehat U_HS_H\widehat U_H^\dagger S_H.
```

Here $`\Pi`$ tests every initialized live Z bit: tag, source index,
inner attenuation, angle program, and outer attenuation. It excludes X,
the dirty pool, and H. Unused Z bits are spectators, and the comparison
extension above is identity on them. Compile this selective phase by
the exact returned-helper zero predicate, with its controls additionally
including H. Its polynomial-logarithmic count and
$`O_\eta(\log\log(n+2))`$ T-depth are charged below; the reflection
is not an oracle assumption.

$`Z_H`$ supplies precisely the active minus. On H equal to zero every
factor is identity on arbitrary source/program inputs, so the iterate
is identity too. A scalar global minus would give the wrong relative
phase between original Z sectors and must not replace $`Z_H`$.

The native word contains $`2k+1=O(\sqrt n)`$ actual forward/inverse
source calls and $`2k`$ selective phases. The fixed damping rotation is
compiled by the exact inactive echo in the loader's guarded-rotation
lemma; its active error is a full operator error.

Give the source and damping errors the common bound

```math
\delta=\frac{\eta}{64(2k+1)}.
```

The source comparison costs $`6\delta`$ and damping costs
$`\delta`$. Unitary telescoping, including actual inverse occurrences,
therefore bounds the recovered early error by

```math
(2k+1)(6\delta+\delta)=7\eta/64<\eta/2.
```

Finally reverse the initial exact H predicate. The ideal active output
has all Z work zero and returns H to zero; the original nonzero-Z sector
was exactly unchanged. Applying the same actual predicate inverse to
both compared outputs preserves the norm bound and includes any residual
tag, source, damping, or H leakage. This proves the complete early
initialized-isometry bound without an intermediate reset.

## 6. Native depth and exponentially smaller early tables

The [query lemma](SHALLOW_DIRTY_QUERIES.md) implements arbitrary m-bit
tables on an r-bit read-only address with no initialized helper, using

```math
T,w=O(m(r+2)^3\sqrt{2^r}),\qquad
G=O(m2^r+m(r+2)^3\sqrt{2^r}),\qquad
D_T=O(\log_2(r+2)).
```

It proves that non-Clifford supports are private even when address
controls are shared. All elementary Clifford interactions are counted.
Boolean phase and variable-suffix Hadamard echoes have literal full-space
return, and have exact inactive identity when their table rows vanish.

Here $`\delta^{-1}=\mathrm{poly}_\eta(n)`$. One descendant preparation
uses a constant number of large Boolean phase queries and actual inverses;
its source reflection is another such query. It has
$`O_\eta(\log\log(n+2))`$ fixed coefficient rotations, whose total
serial synthesis depth is $`O_\eta(\log(n+2)\log\log(n+2))`$.
There is only one source-label-dependent calibration rotation per
unamplified loader call. Its binary angle has $`O_\eta(\log(n+2))`$
bits; serial exact guarded half-angle echoes give
$`O_\eta(\log^2(n+2))`$ depth including synthesis. Constantly many
inner amplification calls preserve this bound. The tag packing, variable
Hadamards, row Hadamards, outer damping and selective phases also fit it.
Thus

```math
D_{T,E}=O_\eta(\sqrt n\,\log^2(n+2)).
```

All table row counts are at most $`Q\,\mathrm{poly}_\eta(n)`$;
their output widths, number of uses and all additional gate/dirty costs
are fixed polynomials in n. Consequently there are fixed exponents for
which the entire early segment has the bounds

```math
T_E\le\mathrm{poly}_\eta(n)\sqrt N\,2^{-M/2}
        +\mathrm{poly}_\eta(n)=o_\eta(\sqrt N),
```

```math
G_E\le\mathrm{poly}_\eta(n)N2^{-M}
        +\mathrm{poly}_\eta(n)=o_\eta(N),
```

```math
w_E\le\mathrm{poly}_\eta(n)\sqrt N\,2^{-M/2}
        +\mathrm{poly}_\eta(n)=o_\eta(\sqrt N).
```

The same dirty pool is reused after completed return, so peak width is
not summed over calls. These limits follow from
$`M/\log(n+2)\longrightarrow\infty`$. This suppression pays for
native data access; it cannot remove the number of sequential recovery
calls, which has been counted separately.

## 7. The final M layers retain the count budgets

Put $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and
$`L'=L+2`$. Use the existing
[capped full-input operator-source compiler](../../docs/AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
on the final M original layers. For height $`1\le j\le M`$ above
the leaves set

```math
Q_j=4N2^{-j},\qquad r_j=n-j+2,\qquad
m_j=L+6+\min\{j,\lceil\log_2(8n)\rceil\}.
```

Its complete-isometry error per layer is at most
$`10\sqrt2\,2^{-m_j}`$. The retained geometric-plus-capped error
sum over all n layers is strictly less than $`2^{-L'}\le\eta/4`$;
this subset has no larger error. Sources, suffix conditioning, actual
inverses and rejected-space error retain that full contract.

Keep the base dirty reservation $`B_0'=L+n+9`$ and the two external
clean flags. The source fits because $`m_j+r_j+1\le B_0'`$.
Choose $`C_\eta`$ sufficiently large that for large n the extra
query pool $`B=b-B_0'`$ satisfies
$`B\ge b/2`$ and $`B\ge16c_I(n+3)`$, where $`c_I`$ is the
fixed constant of the [nonuniform indicator](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).

That retained indicator has genuinely linear count and dirty width,
and logarithmic address T-depth. Its nonuniform chunk sum is bounded by
$`(67^3+2/65)2^r`$ before fixed native constants. No polynomial-in-n
indicator overhead is charged at the bottom tables of size comparable
to N. The resulting exact rectangular query has

```math
T_j=O(m_j+\sqrt{Q_jm_j}+Q_jm_j/B),\qquad G_j=O(Q_jm_j),
```

```math
D_j=O(m_j[1+Q_j/B^2]+\log_2(r_j+2)).
```

Each amplified layer has a constant number of completed queries and
actual inverses. The same pool is reused between them. Its returned
predicate helpers are borrowed only between completed queries.

Since $`m_j\le L+6+j`$ and $`m_j=O_\eta(\log(n+2))`$,
convergent geometric sums give

```math
\sum_{j=1}^M Q_jm_j=O_\eta(N),\qquad
\sum_{j=1}^M\sqrt{Q_jm_j}=O_\eta(\sqrt N),\qquad
\sum_{j=1}^M m_j=O_\eta(M\log(n+2)).
```

Native source words add $`O(\sum_jm_j)`$ count and depth. Exact
suffix predicates add $`O(M\log(n+2))`$ depth and polynomial count;
two-clean reflections cost a constant per layer. Thus the same tail has

```math
T_T=O_\eta(\sqrt N),\qquad G_T=O_\eta(N),\qquad
D_{T,T}=O_\eta(M\log(n+2)+N/B^2)=O_\eta(\log^3(n+2)).
```

Only the early segment uses the protected logical bank. The tail uses
the original arbitrary-dirty operator source and has no supplied unary
phase-state reservation.

## 8. Compose and cover every dimension

The early and tail dirty reservations are used sequentially. Enlarge
$`C_\eta`$ to cover their maximum. Their errors compose in the full
initialized-isometry norm: apply the actual tail unitary to the early
error, then compare the tail on the ideal returned-work embedding.
This gives total error less than $`\eta/2+\eta/4<\eta`$ and requires
no physical reset after the approximate early circuit.

The early depth dominates the tail because
$`\log(n+2)=O(\sqrt n)`$. Their T and elementary Clifford counts
give the stated simultaneous bounds. For the finitely many n before
$`h\ge1`$, the scratch fit, and the asymptotic width inequalities hold,
use the retained [two-clean full-frame compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary).
At fixed eta, enlarging $`C_\eta,c_\eta`$ covers all these dimensions
and their literal base reservations. This is an asymptotic theorem, not
a practical crossover estimate.

## 9. Scope and evidence

The proof uses three distinct ingredients: the new coherent descendant
source, the native shallow query schedule, and the retained linear-size
indicator for the late layers. The column encoding and robust recovery
above supply the missing complete-operator step. A state source alone
would not prove the theorem.

The result closes the previously stated $`o(n)`$ existence question at
fixed accuracy and sufficiently large constant square-root dirty width.
It leaves the optimal depth order and other width/precision regimes open.
The [five bounded recovery checks](../../tests/test_sublinear_depth_recovery.py)
test signed and singular complete columns, marker tags, actual inverses,
active-only phase, and final leaked-work return. Their small numerical
unitary completions are verification devices, not compiler primitives.
They are not evidence for the asymptotic count absorption or a claim that
the resulting exponentially sized circuits have been emitted.
