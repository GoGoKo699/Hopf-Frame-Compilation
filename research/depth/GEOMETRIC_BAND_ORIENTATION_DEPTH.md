# Geometric-band recovery for arbitrary branch orientations

[Nested subtree spaces](SUBTREE_PROJECTOR_DECODER.md) · [Homogeneous flag transport](HOMOGENEOUS_MAGNITUDE_DEPTH.md) · [Recursive terminal](RECURSIVE_TERMINAL_DEPTH.md) · [Joint session design](JOINT_ORIENTATION_RECOVERY_PROMPT.md)

**Local result R85.** The
unordered rotation magnitudes may depend on depth, and their ordering may
be chosen by an arbitrary Boolean function of the node prefix. A complete
frame in this family has fixed-accuracy T-depth
$`O_L(N/b^2+\log(n+2)\log\log(n+4))`$, retaining one external clean
qubit, the inherited T-count bound, and linear Clifford count.

The new resource is conditional clean space exposed by the source itself.
A coarse band label identifies a descendant loader whose inverse maps
that entire band into a shorter frame and a zero suffix. Geometric band
widths make that suffix large enough for a local clean output word.
This permits an exact local permutation, although no general dirty-only
in-place permutation compiler is supplied. The source's old fine tag
chooses the coarse band; only the much shorter coarse tag needs spectral
recovery afterward.

## 1. Promised family and theorem

Number logical positions from 1 to n. At depth j assume that one
nonnegative normalized pair $`(a_j,b_j)`$ satisfies, for every prefix p,

```math
(|c_{j,p}|,|s_{j,p}|)=
\begin{cases}
(a_j,b_j),&r_j(p)=0,\\
(b_j,a_j),&r_j(p)=1,
\end{cases}
\qquad a_j^2+b_j^2=1.
```

Each $`r_j`$ is an arbitrary Boolean function. Node signs are arbitrary;
zeros and balanced pairs are allowed. The promise is existence of these
pairs and finite Boolean tables. They need not be supplied, and there
is no efficiency claim for their classical discovery.

Fix an integer $`L\ge6`$, put $`\eta=2^{-L}`$,
$`N=2^n`$, and $`\ell=\log_2(n+2)`$. For every integer
$`n\ge1`$ and every indicated arbitrary dirty allocation, one external
initialized qubit suffices:

| Prescribed operator W | Required arbitrary dirty width |
|---|---|
| Complete real frame with the promised unordered pairs | $`b\ge17(L+n+7)`$ |
| Prescribed magnitude completion $`D_\phi W_{\mathbb R}`$ with that real-frame promise | $`b\ge17(L+n+8)`$ |

There is one Clifford+T circuit V satisfying simultaneously

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell\log_2(\ell+2)\right).
```

Here $`J_1`$ appends the sole initialized qubit. Constants depend only
on fixed L, uniformly over dimension, magnitudes, orientations, signs,
and eligible allocations. The norm includes every logical column,
literal phase, final work leakage and return, arbitrary dirty inputs,
and external references. Unused supplied wires remain untouched.
G counts elementary Clifford gates, including the arbitrary Clifford
interlayers permitted by the T-depth model. Every source, query, control,
preparation, actual inverse, and recovery is charged. There is no QRAM,
measurement, reset, supplied precision state, or free permutation oracle.
The complex claim concerns the prescribed magnitude completion only.

At $`b=\Theta(\sqrt N)`$ the depth is

```math
O_L(\log(n+2)\log\log(n+4))=o(\log^2(n+2)),
```

with square-root T-count and linear Clifford count. This is a supporting
family upper bound, not an optimal-depth claim or an all-precision bound
with absolute constants. Arbitrary node-dependent unordered magnitude
pairs remain outside this theorem.

## 2. Certified finite reduction, including zeros and signs

Use the finite rounding argument of
[R83, Section 2](AFFINE_ORIENTATION_DEPTH.md#2-finite-rounding-with-affine-orientations-and-arbitrary-signs),
replacing its finite list of affine masks by the finite list of all
Boolean tables at the given depth. Set $`\epsilon=\eta/(64n)`$.
Enumerate a finite angle grid in $`[0,\pi/2]`$ with mesh at most
$`\epsilon/2`$, and dovetail certified strict component tests for every
grid angle and candidate orientation table at that depth. A grid point
within $`\epsilon/4`$ of a true common angle, together with one true
orientation table, passes with a strict margin. Hence this finite search
terminates even when orientations are ambiguous at balanced pairs.

Choose signs from rational approximants of the signed coefficients at
absolute accuracy epsilon. A wrong sign implies that the original
coefficient has magnitude at most epsilon. The accepted finite frame
$`W'`$ has coordinate errors at most $`3\epsilon`$, local rotation
errors at most $`3\sqrt2\epsilon`$, and therefore

```math
\|W-W'\|\le3\sqrt2 n\epsilon\lt \eta/8.
```

Its unordered magnitudes are exactly common within each depth, and its
orientation and sign tables are finite. The
[sign gauge](HOMOGENEOUS_MAGNITUDE_DEPTH.md#2-separate-signs-without-a-sign-oracle)
gives the exact full-column identity

```math
W'=D_{\rm sign}W_{\rm mag}S_{\rm sign}.
```

The output sign is the product of edge signs along the leaf path. The
root input sign is one; a marker input sign is its ancestor-sign product
times both edge signs at its own node. This identity divides by no
amplitude and permits either sign on a zero edge. Both exterior signs
are priced exact Boolean phases. Sections 3–8 treat the positive finite
frame $`W_{\rm mag}`$; Section 9 composes these reductions and the
prescribed leaf phase.

## 3. Protected registers and geometric bins

For sufficiently large n choose

```math
M=\lceil A_L\ell\rceil\lt n,\qquad h=n-M,\qquad Q=2^h,
```

where the fixed constant $`A_L`$ will dominate all clean allocations
and polynomial resource envelopes. The finite remaining dimensions
use Section 9's fallback. Split the original logical register into its
first h bits X and last M bits Z. The exact early/terminal factorization
from R78 is

```math
W=W_TW_E,\qquad
W_E=W_h\otimes|0^M\rangle\langle0^M|_Z
     +I_X\otimes(I_Z-|0^M\rangle\langle0^M|_Z).
```

Before occupying Z, save $`H=[Z=0^M]`$ in the sole external clean
qubit. The physical word never retests Z while it is work. On H one,
reserve disjoint portions of Z for the following registers:

| Register | Purpose | Size |
|---|---|---|
| T | Old fine source tag | $`s=\lceil\log_2(h+1)\rceil`$ |
| R | Coarse bin, encoded from zero | $`k=\lceil\log_2 p\rceil`$ |
| Y | Fresh coarse-search record | k |
| D | Selected descendant threshold | s |
| A | Complete fresh initialized source bank | $`O(\log(1/\delta))`$ |
| Other work | Short exact arithmetic and fixed interface bits | $`O(s+k+1)`$ |
| $`Z_{\rm free}`$ | Untouched reserve available to local permutations | $`m_0=\Theta_L(\ell)`$ |

The bin count p and tolerance delta are specified below. Taking
$`A_L`$ sufficiently large accommodates this entire table with positive
$`m_0\ge2`$. No source bank, tag, threshold, or search record belongs
to $`Z_{\rm free}`$. All these protected registers may be occupied
arbitrarily on H zero.

Choose boundaries

```math
0=d_0\lt d_1\lt \cdots\lt d_p=h,
\qquad a_i=d_{i-1},\quad b_i=d_i,\quad q_i=b_i-a_i,
```

such that

```math
q_i\le\frac{h-b_i+m_0}{2}.
\qquad\text{(1)}
```

A concrete choice starts backward at $`b=h,r=m_0`$. Set
$`q=\min\{b,\lfloor r/2\rfloor\}`$, prepend the boundary
$`b-q`$, and replace $`(b,r)`$ by $`(b-q,r+q)`$ until b is zero.
Apart from the final truncated step, r grows by a fixed factor greater
than one; thus

```math
p=O(1+\log(1+h/m_0))=O(\ell),\qquad
k=O(\log(\ell+2)).
```

The first bin contains $`0\le t\le d_1`$; bin i greater than one
contains $`d_{i-1}\lt t\le d_i`$. Define its zero-based label by

```math
g(t)=\min\{i\in\{1,\ldots,p\}:t\le d_i\}-1.
```

Padded tag values receive any fixed valid bin label. For each i choose
$`q_i`$ distinct physical wires $`B_i`$ from
$`X_{b_i+1:h}`$ and $`Z_{\rm free}`$. Equation (1) gives twice the
required capacity. These wires exclude every preserved record and source
work bit. They are claimed clean only on the ideal reached bin-i sector.

A separate fixed word O of $`q_{\max}=\max_iq_i\le h`$ arbitrary
dirty qubits will serve as a physically disjoint query output. O is not
initialized. Its temporary exchange with $`B_i`$ is the only clean
output word used below.

## 4. Unloading a whole band exposes the required zeros

Use the literal nested scaling spaces of R78. At depth d, let
$`v_{d,u}`$ be the full h-bit subtree scaling vector with prefix u,
and define

```math
A_d=\sum_{|u|=d}|v_{d,u}\rangle\langle v_{d,u}|,
\qquad \Pi_0=A_0,\quad \Pi_t=A_t-A_{t-1}\ (t\ge1).
```

For a selected d the descendant loader has anchor and output isometries

```math
J_d|u\rangle=|u0^{h-d}\rangle_X|0\rangle_A,
\qquad B_d|u\rangle=v_{d,u}|0\rangle_A,
```

and, on H one, a native unitary $`\widehat U_d`$ satisfying

```math
\|\widehat U_dJ_d-B_d\|\le\delta_U,
\qquad
\|\widehat U_d^\dagger B_d-J_d\|\le\delta_U.
\qquad\text{(2)}
```

It preserves D and the prefix, returns its arbitrary dirty pool exactly,
and is exactly identity on the full H-zero space. Its dagger is its
actual reversed word. The old T and coarse R are spectators.

Let $`W^{[d]}`$ be the upper d-level frame, with the original
coefficients at depths 1 through d. Expanding the defining scaling and
marker recurrences shows, for every column whose band depth is at most d,

```math
|w_j\rangle_X|0\rangle_A
   =B_d\,|w^{[d]}_{j_{1:d}}\rangle.
\qquad\text{(3)}
```

Equivalently, $`B_d`$ transports the entire upper frame flag to the
full nested flag through depth d. This is an isometry identity on the
whole subspace $`\mathrm{Ran}\,A_d`$, so it gives an operator-norm
premise for (2), not just separate column estimates.

In bin i greater than one, the depth t is larger than $`a_i`$.
Consequently the upper frame column in (3) has an unchanged computational
prefix of length $`a_i`$; its next $`q_i`$ bits belong to the local
oriented frame band of depth $`t-a_i`$. Its remaining $`h-b_i`$
bits are zero after $`J_{b_i}`$. For the first bin $`a_1=0`$,
the same statement includes the root column and local depth zero.
Thus $`\widehat U_{b_i}^\dagger`$ exposes the zeros required by
(1), with only its source-isometry error. No assumption is made about
its completion on the other inputs.

## 5. A local flag transporter with a paid clean-output permutation

### Local flag identity

Fix i and an unchanged physical prefix u of length $`a_i`$. On the
next $`q_i`$ bits define the triangular permutation recursively by

```math
y_j=z_j\mathbin\oplus r_{a_i+j}(u y_{\lt j}),\qquad
F_{i,u}(z)=y,\qquad P_{i,u}|z\rangle=|F_{i,u}(z)\rangle.
\qquad\text{(4)}
```

It is bijective, since each next z bit is recovered from the corresponding
y bit and its already known prefix. It maps each homogeneous scaling
vector to its oriented scaling vector, with a permutation of prefix
labels at each depth. Therefore its local scaling projectors satisfy

```math
A^{\rm loc}_j=P_{i,u}A^{\rm hom}_jP_{i,u}^\dagger.
```

Let $`R(\theta_j)`$ have first column $`(a_j,b_j)^\mathsf T`$,
put

```math
U_{0,i}=\bigotimes_{a_i\lt j\le b_i}R(\theta_j),\qquad
C_i=P_{i,u}U_{0,i},
```

and write

```math
C^{\rm loc}_j=I_{2^j}\otimes|0^{q_i-j}\rangle\langle0^{q_i-j}|.
```

Then

```math
C_i C^{\rm loc}_j C_i^\dagger=A^{\rm loc}_j.
\qquad\text{(5)}
```

It follows that $`C_i^\dagger`$ takes local band depth j to the
computational last-one band of depth j. It may mix vectors within that
band. This does not affect the label XOR, and the subsequent actual
inverse restores the full vector and its literal phase. In particular,
no separate signed marker permutation is needed. Equation (5) uses no
division and includes zero or balanced magnitude pairs.

### Exact permutation word and all dirty inputs

On the initialized-tail sector exposed in Section 4, implement
$`P_{i,u}`$ by the following chronological word, selected by
$`H=1,R=i-1`$:

1. Swap $`O_{1:q_i}`$ with $`B_i`$. The output O now has a zero
   first word, and $`B_i`$ holds the original arbitrary dirty word.
2. XOR $`F_{i,u}(X_{a_i+1:b_i})`$ into $`O_{1:q_i}`$.
3. Swap $`X_{a_i+1:b_i}`$ with $`O_{1:q_i}`$.
4. XOR $`F_{i,u}^{-1}(X_{a_i+1:b_i})`$ into $`O_{1:q_i}`$.
5. Swap $`O_{1:q_i}`$ with $`B_i`$ again.

The local data becomes its F image, O returns its original dirty word,
and $`B_i`$ returns zero. This is a basis-permutation identity for
every initial dirty word and hence holds with arbitrary entanglement
and references. Its actual reverse implements the inverse permutation
on that same whole initialized-tail sector.

In fact the external word O returns on the entire input space. If the
local data, parked tail, and external dirty words initially have values
$`x,r,y`$, the output is

```math
x'=r\mathbin\oplus F_{i,u}(x),\qquad
r'=x\mathbin\oplus F_{i,u}^{-1}(r\mathbin\oplus F_{i,u}(x)),
\qquad y'=y.
```

Thus the initialized-tail premise is needed for the local permutation
action and zero-tail return, not for restoring the external dirty output.
Any actual logical-tail leakage remains part of the outer error bound.

Both XOR tables read the fixed physical register X and the preserved
R,H; their target is the physically disjoint O. On sector i they ignore
all X positions beyond $`b_i`$, so the temporarily moved dirty word
cannot affect an address. They output zero in positions beyond $`q_i`$
and when H is zero. Padded R values select identity. There are
$`O(p2^h)`$ rows and at most n output bits. Exact bit-parallel native
queries therefore implement both tables with polynomial overhead in n.
No clean h-bit index or supplied in-place permutation is used.

For completeness, each selected swap above has the following native
schedule. For each bin i, its $`q_i`$ physical pairs are disjoint.
Compute $`f_i(H,R)=H[R=i-1]`$ into a separate arbitrary dirty control
for each pair, apply the corresponding controlled SWAPs in parallel,
reverse the query, and apply the same controlled SWAPs again. For a
pair involution S and dirty control c, this is exactly

```math
Q_{f_i},\ S^c,\ Q_{f_i}^\dagger,\ S^c
\quad=\quad S^{f_i},
\qquad\text{(6)}
```

because the read-only address is unaffected. The dirty controls and
helpers return exactly. Within one bin all non-Clifford supports are
private; shared R,H occur only as read-only query controls. An elementary
controlled SWAP has constant T-depth and constant Clifford+T size.
Process the p bins serially, reusing these private pools. Since only one
bin is selected, their exact product is the desired selected matching.
The schedule is identity at H zero even on occupied protected work.

The short-address query bound gives, for all bins of one matching,

```math
T,w=O\!\left(n(k+3)^3\sqrt p\right),\qquad
G=O\!\left(np+n(k+3)^3\sqrt p\right),\qquad
D_T=O\!\left(p\log_2(k+3)\right).
\qquad\text{(7)}
```

Here the count uses $`\sum_iq_i=h`$; the width sum is a conservative
bound. Thus the selected routing costs at most
$`O(\ell\log(\ell+2))`$ depth and polynomial resources, without a
free selected Clifford network.

## 6. Parallel selected rotations and the fine-tag XOR

Select intervals without selecting different angles in a serial loop.
For each physical position j choose an ancilla-free native word
$`V_j`$ approximating $`R(\theta_j/2)`$ to error $`e_j`$.
Define $`V=\bigotimes_{j=1}^hV_j`$ and the exact Boolean phase

```math
Q_I|x,R=i-1,H\rangle
 =(-1)^{H\sum_{a_i\lt j\le b_i}x_j}|x,R=i-1,H\rangle,
```

where the exponent is a binary sum and padded bins select zero. It is
one paid phase table on X,R,H. Set

```math
\widehat U_0=VQ_IV^\dagger Q_I.
\qquad\text{(8)}
```

Chronological execution in (8) is $`Q_I,V^\dagger,Q_I,V`$.
On a selected position its local factor is
$`V_jZV_j^\dagger Z`$; on every other position its factor is exactly
identity. Since

```math
R(\theta/2)ZR(-\theta/2)Z=R(\theta),
```

we have $`\|\widehat U_0-U_{0,i}\|\le2\sum_je_j`$ on sector i.
All phases of the synthesized words cancel against their actual inverses.
Equation (8) is exactly identity on H zero, even on arbitrary data and
occupied protected records. Taking $`e_j=\delta_C/(2h)`$ gives
error at most $`\delta_C`$ and parallel one-qubit T-depth
$`O(\log(h/\delta_C+2))`$.

Combine this word with Section 5's exact selected permutation to obtain
$`\widehat C_i=P_{i,u}\widehat U_0`$. The equality denotes its
initialized-tail action, not a supplied oracle. Its physical word
preserves that entire anchor sector exactly, returns every arbitrary
dirty input exactly, and is within $`\delta_C`$ of $`C_i`$ there.
Its actual reverse has the same approximation property.

For local data x define

```math
f_i(x)=a_i+\mathrm{lastone}(x_{a_i+1:b_i}).
\qquad\text{(9)}
```

Last-one is zero on the all-zero word. In the first bin this gives
$`f_1(0)=0`$, as required for the root. In later bins local depth zero
is absent from the ideal reached sector; any output may be assigned
there. An exact bit-parallel XOR table, reading X,R,H and targeting T,
implements $`T\mathrel{\oplus}=Hf_i(X)`$. It never reads T.
Its $`O(\log n)`$ output bits, $`O(p2^h)`$ rows, and
$`O(\log n)`$ T-depth are charged by the same query ledger.

## 7. Complete recovery word, phases, and operator error

The packed source from R78 has ideal isometry

```math
S|j\rangle=|w_j\rangle_X|t(j)\rangle_T
            |0\rangle_{\text{other protected work}},
\qquad \|\widehat S-S\|\le\delta_0.
\qquad\text{(10)}
```

Its completed physical word is exactly identity at H zero. Starting
from (10), compute $`R\mathrel{\oplus}=Hg(T)`$ by an exact short
table. This uses the original source tag directly and leaves T unchanged.
Compute $`D=b_{R+1}`$ in its zero record. On each selected bin execute
chronologically

```math
\widehat U_{b_i}^\dagger,\quad
\widehat C_i^\dagger,\quad
T\mathrel{\oplus}=Hf_i(X),\quad
\widehat C_i,\quad
\widehat U_{b_i}.
\qquad\text{(11)}
```

Then uncompute D using the unchanged R. The loader source bank, every
selected operation, and both query addresses are disjoint from T. In
the ideal comparison, (3) exposes the zeros; (5) puts the local data in
its computational last-one band; and (9) clears T. Reversing the local
transport and reloading restores the original complete frame column
with its literal phase. All dynamic clean space and every dirty input
return. Only $`R=g(t(j))`$ remains.

More precisely, (2) bounds the first actual inverse-loader error on the
whole ideal source image by $`\delta_U`$. The actual inverse local
transport adds at most $`\delta_C`$ on the ideal anchor sector; the
exact XOR gives the correct common depth there. The forward local
transport and forward loader add at most the same respective errors.
Thus (11) differs from its ideal reached-sector action by at most

```math
2\delta_U+2\delta_C.
\qquad\text{(12)}
```

This is a uniform operator-norm estimate on the complete selected source
image. Orthogonal bin and preserved-prefix sectors take a maximum,
not a sum over columns. Actual earlier leakage is propagated through
unitaries; it is never assumed initialized when invoking the local
clean-tail lemma. The source loader's completion outside its image is
irrelevant. All daggers in (11) are actual reversed native words.

It remains to erase R. Coherently search for its zero-based label among
p ordered bins using $`k=\lceil\log_2p\rceil`$ comparisons. A
comparison of the unknown bin label with an integer m uses

```math
[g(t)\le m]=[t\le d_{\min\{m+1,p\}}].
\qquad\text{(13)}
```

Thus padded endpoints use $`d_p=h`$ and always return true. Each
threshold is an exact finite function of the previously recovered search
prefix and is computed into D and uncomputed. The R78 Hadamard/controlled
reflection/complement-bit word writes the coarse label into fresh Y,
with no scalar phase, while leaving every data vector in its original
band. The old fine tag T and the old coarse R are not queried.

Run this search $`\widehat Q`$, XOR HY into R, and reverse the actual
search. The loader reflections use their full initialized-source anchor,
including A. On H one, R78's un-controlled transported-projector
comparison holds on the full space, including occupied A:

```math
\widehat U_d(I-2J_dJ_d^\dagger)\widehat U_d^\dagger
\quad\hbox{is within }4\delta_R\hbox{ of}\quad
I-2(A_d\otimes|0\rangle\langle0|_A).
```

The comparison is identity off that transported zero-work subspace;
it is not $`(I-2A_d)\otimes I_A`$ on occupied A. Including H and the
current search bit in the actual central predicate makes both the
actual and ideal controlled reflections identity when either control
is zero, and retains the global $`4\delta_R`$ bound. Thus

```math
\|\widehat Q-Q_*\|\le4k\delta_R,
\qquad \|\widehat E_{\mathrm{coarse}}-E_{\mathrm{coarse},*}\|
       \le8k\delta_R.
\qquad\text{(14)}
```

The ideal word clears R and Y, returns every source and threshold record,
and preserves the frame data with phase +1. It does not read the erased T.
The exact comparison reflections nest; the approximate physical ones
need not nest or commute, because (14) uses unitary telescoping in full
operator norm.

Combining (10), (12), and (14) gives the complete early initialized error

```math
\epsilon_E\le\delta_0+2\delta_U+2\delta_C+8k\delta_R.
\qquad\text{(15)}
```

For example choose

```math
\delta_0=\delta_U=\delta_C=\eta/128,\qquad
\delta_R=\eta/[512(k+1)],
```

so $`\epsilon_E\lt \eta/8`$. The initial source and local unloaders
use fixed-L precision; the finest comparison source tolerance has
$`\log(1/\delta_R)=O_L(\log(k+2))`$.

On H zero, all local permutation and rotation words and the fine XOR are
exactly identity. Threshold preparation meets its actual inverse. The
coarse erasure is its actual search/inverse sandwich with disabled
middle XOR, so it is identity even on occupied records. Every completed
source is itself H-zero identity. Thus the full early circuit has
literal inactive phase +1. On its ideal active output all Z bits are
again zero. Apply the actual saved-predicate inverse; it returns H on
that ideal output and propagates (15) unchanged on the actual output.
This proves the full early initialized-isometry contract, including
arbitrary dirty systems and references.

## 8. Native resource ledger

Use the fresh per-call descendant loader, including its own gradient
preparation and actual inverse. No retained-gradient extension is needed.
The absolute ledger from R79 gives, at source tolerance delta,

```math
T,w\le C(h+\log(1/\delta)+2)^6\delta^{-1}(\sqrt Q+1),
```

```math
G\le C(h+\log(1/\delta)+2)^6\delta^{-2}(Q+1),\qquad
D_T\le C[\log(h+2)+\log(1/\delta)].
\qquad\text{(16)}
```

Its source bank has $`O(\log(1/\delta))`$ initialized bits; the
central reflection tests all that work. Its selected source labels and
preserved-prefix direct sums are already part of this interface. Adding
the short bin selector changes rows and counts by a polynomial factor.

The present word has one initial source, two unloader appearances in
(11), and $`O(k)`$ reflected-source calls in the final coarse erasure.
It has a constant number of full-X forward/inverse XOR and phase queries,
bit-parallel fine XOR, and the polynomial selected-swap schedule (7).
The fresh source tolerance in (16) is at least
$`\eta/[512(k+1)]`$. A conservative common envelope is therefore

```math
T_E,w_E\le C_L(n+2)^{10}(\sqrt Q+1),\qquad
G_E\le C_L(n+2)^{10}(Q+1),
```

```math
D_{T,E}\le C_L\ell\log_2(\ell+2).
\qquad\text{(17)}
```

For the exact local F and inverse-F queries alone, at most
$`4p2^h`$ rows and n outputs give
$`T,w\le Cn^4\sqrt p\sqrt Q`$ and
$`G\le CnpQ+Cn^4\sqrt p\sqrt Q`$, inside (17).
Shared read-only addresses appear only as Clifford controls; each
parallel output has its private dirty non-Clifford supports. This is the
bit-parallel query bound, not a serial output-word implementation.
The one-qubit words in (8) run in parallel and contribute
$`O_L(\ell)`$ depth. Exact saved predicates and short arithmetic
have polynomial resources and $`O(\ell)`$ depth. All are included
in (17).

Take $`A_L`$ sufficiently large, in particular larger than 40 after
increasing it for the simultaneous clean register constants. Since
$`Q=N2^{-M}`$, (17) has

```math
T_E=O_L(\sqrt N),\qquad G_E=O_L(N),\qquad
w_E=o(\sqrt N/n).
```

These estimates include the arbitrary dirty output O, every private
query/switch pool, and all actual inverses. They hold uniformly over
angles and orientation tables beyond one L-dependent dimension cutoff.
The live initialized work is precisely the protected allocation in
Section 3; no new external initialized qubit is introduced.

## 9. Terminal composition, literal widths, and finite fallback

First suppose $`b\ge\sqrt N/n`$ and n exceeds the fixed cutoff.
The early pools in (17) fit within a fixed fraction of b. The
[R81 terminal theorem](RECURSIVE_TERMINAL_DEPTH.md#1-fixed-accuracy-terminal-theorem)
applies to the remaining $`M=O_L(\ell)`$ levels, with their
preserved upper prefix, arbitrary coefficients, one returned external
clean qubit, and the same original dirty pool reused. At a fixed fraction
of the total accuracy it gives

```math
T_T=O_L(\sqrt N+N/b),\qquad G_T=O_L(N),\qquad
D_{T,T}=O_L(N/b^2+\ell\log_2(\ell+2)).
```

No extension of its terminal-height range is required here. The terminal
word composes with the early initialized-isometry bound by unitary error
propagation, including possible actual leakage of the returned external
qubit. Its exact target is $`W_TW_E`$.

The two exact sign gauges from Section 2 use the paid all-width Boolean
phase construction of R80, with
$`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and
$`D_T=O_L(N/b^2+\ell)`$. In the complex case use the retained R72
prescribed diagonal phase word at a fixed fraction of the total error;
its bounds have the same fixed-L orders. Reuse the sole initialized
signal and dirty pools between completed factors. The original complex
literal reservation is $`17(L+n+8)`$. All stricter internal
precision reservations fit in the wide regime beyond the same fixed
cutoff. Allocate fixed fractions of eta to finite rounding, the early
word, terminal word, and prescribed phase. Their initialized norm errors
sum to at most eta.

For $`b\lt \sqrt N/n`$, invoke the retained whole-frame R72 theorem
at the original total eta, not an unnecessarily stricter internal
precision. It has the original literal reservation, the required T and
G bounds, and

```math
D_T=O_L(N/b^2+nL)=O_L(N/b^2),
```

because $`N/b^2\gt n^2`$. This fallback applies to every frame, so its
use does not require the orientation promise. For the finitely many n
below the fixed cutoff, including any case where the protected choice
would leave no early block, use that same theorem at the original total
accuracy. Increasing an L-dependent constant absorbs their finite
additive depths uniformly over all angles and allocations.

These branches preserve the literal widths stated in Section 1 and
establish the theorem for all n and every eligible b. Unused supplied
wires remain untouched. The improvement comes from complete source and
recovery together; it does not provide a standalone arbitrary triangular
permutation compiler with only dirty output space.

## 10. Scope and verification

The proof uses the common unordered magnitude pair at each physical
depth only in the local tensor-product transporter (5). With arbitrary
node-dependent magnitude pairs that local factor is no longer a tensor
product of fixed one-qubit rotations. No analogous bound for that broader
family is asserted here. The fixed-L dependence is also essential to
the stated resource constants; no absolute all-precision theorem is
inferred by substituting a growing L. The later
[R86 interval-dilation theorem](INTERVAL_DILATION_DEPTH.md) gives
$`O_L(N/b^2+\ell^{3/2})`$ for arbitrary node-dependent magnitudes.
The present $`O_L(N/b^2+\ell\log_2(\ell+2))`$ bound remains
stronger within its shared-unordered-pair promise.

The six [bounded mechanism tests](../../tests/test_geometric_band_orientation.py)
check geometric capacity, the full truncated-band identity, local flag
transport, dirty-output restoration, complete fine-tag erasure/reload,
coarse-tag erasure with padded thresholds, and occupied inactive sectors.
They support the symbolic
proof and do not replace its native schedule, error bound, or asymptotic
resource ledger.
