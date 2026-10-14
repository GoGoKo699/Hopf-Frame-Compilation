<a id="sublinear-t-depth-at-square-root-dirty-width"></a>

<a id="sublinear-t-depth-and-a-fixed-accuracy-width-frontier"></a>

# Sublinear T-depth and a width and precision frontier

[Resource statement](../../docs/OPEN_PROBLEM.md) · [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) · [Coherent descendant loader](COHERENT_DESCENDANT_LOADER.md) · [Linear-size dirty indicators](NONUNIFORM_DIRTY_INDICATOR.md)

This construction bypasses the successive early Hopf layers by encoding
their complete columns into orthogonal depth-and-prefix sectors. A
coherent real-state loader implements that encoding; a Clifford row
isometry gives normalization $`O(\sqrt n)`$. Amplitude amplification
recovers the whole early frame. The terminal layers share one binary
geometric precision source with retained failure history, followed by
capped layerwise compilation at the bottom. This retains the terminal
width penalty.
Selecting the retained compiler at smaller widths gives an
all-eligible-width theorem.

Sections 1–9 retain the fixed-accuracy construction. Section 10 proves
an absolute-constant precision-uniform extension, using
[one charged gradient across recovery](SHARED_GRADIENT_RECOVERY.md)
and the [shared binary terminal source](SHARED_BINARY_TERMINAL.md).
The fixed-accuracy and uniform constructions have separate gradient
reflection rules and share the same complete-input contract and native
query lemmas.

The construction uses the Boolean-phase state-preparation method of
Gosset–Kothari–Wu, with a separately proved terminating selection rule,
coherent descendant labels, and arbitrary-dirty query implementations.
It is not an invocation of an uncharged state-preparation oracle.

## 1. Fixed-accuracy theorem

For every fixed $`0\lt\eta\le1/64`$, put

```math
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad
B_0=L+n+7,\qquad f(n)=\sqrt n\,\log_2(n+2).
```

There is a constant $`c_\eta`$ such that for every $`n\ge1`$,
$`N=2^n`$, and integer $`b\ge17B_0`$, every prescribed complete
real Hopf frame W has one coherent Clifford+T circuit V with two
external initialized flags and at most b arbitrary dirty qubits satisfying

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T\le c_\eta\left(\sqrt N+\frac Nb\right),\qquad
G\le c_\eta N,\qquad
D_T\le c_\eta\left(\frac N{b^2}+\min\{n,f(n)\}\right).
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

In particular, there is a constant $`C_\eta`$ such that at
$`b=\lceil C_\eta\sqrt N\rceil`$ the same guarantees become

```math
T=O_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(f(n))=o(n).
```

The inherited worst-case lower bounds match the T-count at every
eligible width. They also match the depth throughout every nonempty interval

```math
17B_0\le b\le\sqrt{\frac N{f(n)}}
=\frac{\sqrt N}{n^{1/4}\sqrt{\log_2(n+2)}}:
\qquad
T^\star=\Theta_\eta(N/b),\qquad
D_T^\star=\Theta_\eta(N/b^2).
```

These are fixed-accuracy statements; their constants may depend on eta.
They leave the precision-uniform constants and original matching intervals
of [Result D](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary)
unchanged. They answer the $`o(n)`$ existence question positively and
extend the fixed-accuracy interval of matching width-sensitive depth.
Optimal depth beyond the matching interval, an unrestricted growing-depth
lower bound at square-root width, and a general native emitter remain open.
Section 10 gives a separate precision-uniform extension through a
sufficiently small absolute linear-precision range.

## 2. Protect one terminal logical block

First consider sufficiently large n and
$`b\ge\sqrt{N/n}`$. Section 8 covers every other eligible allocation
using the retained compiler. Choose

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
reused angle program, gradient core, and outer attenuation. Their sum is
$`O_\eta(\log(n+2))`$, which fits M for sufficiently large n.

## 3. Compress the complete columns

Set

```math
s=\lceil\log_2(h+1)\rceil,\qquad K=2^s\lt2(h+1).
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
inner attenuation, angle program, gradient core, and outer attenuation.
It excludes X, the dirty pool, and H. Unused Z bits are spectators, and the comparison
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
(2k+1)(6\delta+\delta)=7\eta/64\lt\eta/2.
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
its source reflection is another such query. Its
$`O_\eta(\log\log(n+2))`$ fixed coefficient rotations act on distinct
index bits and run in parallel, with depth $`O_\eta(\log(n+2))`$.
Their shared activity control occurs only in Clifford CZ banks.

The source-label-dependent calibration uses the loader's
[signed modular-shift construction](COHERENT_DESCENDANT_LOADER.md#guarded-real-rotations-and-the-calibration-table).
Its angle program and separate gradient core each have
$`O_\eta(\log(n+2))`$ bits. Parallel native product preparation,
two exact helper-free modular adders, Clifford complement banks, and
the actual preparation inverse have total depth
$`O_\eta(\log(n+2))`$. Every calibration prepares and unprepares
its own gradient. The inactive sector is literally identity on arbitrary
program and core inputs; active gradient leakage enters the source error.
There is no supplied phase state, reset, or uncharged reuse assumption.

Constantly many inner amplification calls preserve this depth bound.
The tag packing, variable Hadamards, row Hadamards, outer damping and
selective phases also fit it. Thus

```math
D_{T,E}=O_\eta(\sqrt n\,\log(n+2))=O_\eta(f(n)).
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
        +\mathrm{poly}_\eta(n)=o_\eta(\sqrt{N/n}).
```

The same dirty pool is reused after completed return, so peak width is
not summed over calls. These limits follow from
$`M/\log(n+2)\longrightarrow\infty`$: even after multiplying
the normalized width by $`\sqrt n`$, the exponential suppression
beats every fixed polynomial. Thus all early dirty reservations fit
every $`b\ge\sqrt{N/n}`$ for sufficiently large n. This suppression pays for
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
For large n in this branch the extra query pool $`B=b-B_0'`$ satisfies
$`B\ge b/2`$ and $`B\ge16c_I(n+3)`$, where $`c_I`$ is the
fixed constant of the [nonuniform indicator](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
Both inequalities follow from $`b\ge\sqrt{N/n}`$ and the fact
that the base and helper reservations are only linear in n at fixed eta.
This does not assert that the new indicator fits the literal minimum
$`b=17B_0`$; Section 8 retains the old compiler where it does not fit.

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
T_T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad
G_T=O_\eta(N),\qquad
D_{T,T}=O_\eta\!\left(\log^3(n+2)+\frac N{b^2}\right).
```

Here $`B\ge b/2`$ converts the query penalties to b, and the
unweighted source and predicate counts are polynomial in n and fit the
exponential count budgets. At square-root width this specializes to the
previous $`O_\eta(\log^3(n+2))`$ tail depth.

Only the early segment uses the protected logical bank. The tail uses
the original arbitrary-dirty operator source and has no supplied unary
phase-state reservation.

## 8. Compose and cover every dimension

The early and tail dirty reservations are used sequentially. In the
large-n branch $`b\ge\sqrt{N/n}`$, each fits the supplied allocation,
so their maximum fits too. Their errors compose in the full
initialized-isometry norm: apply the actual tail unitary to the early
error, then compare the tail on the ideal returned-work embedding.
This gives total error less than $`\eta/2+\eta/4\lt\eta`$ and requires
no physical reset after the approximate early circuit.

The additive tail depth is absorbed because
$`\log^3(n+2)=O(f(n))`$. Summing the same circuit's resources gives

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(\frac N{b^2}+f(n)\right).
```

Choose a fixed threshold $`n_0(\eta)`$ after which $`h\ge1`$,
the conditional scratch fits, and every width inequality used above
holds whenever $`b\ge\sqrt{N/n}`$. This threshold is independent
of b and of the prescribed frame. For $`n\ge n_0(\eta)`$ and
$`17B_0\le b\lt\sqrt{N/n}`$, use the retained
[two-clean full-frame compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary).
At fixed eta its bounds are

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(\frac N{b^2}+n\right).
```

The original additive T-count term n is absorbed by $`\sqrt N`$.
In this width branch $`n\le N/b^2`$, so its depth is
$`O_\eta(N/b^2)`$ and satisfies the desired bound.

For the finitely many $`n\lt n_0(\eta)`$, use the same retained
compiler at every $`b\ge17B_0`$. The finite maximum of
$`n/f(n)`$ is absorbed into a constant depending only on eta,
uniformly for all those b. This preserves the literal base allocation:
no unproved claim about the new indicator at $`17B_0`$ is needed.
Consequently the displayed $`N/b^2+f(n)`$ upper bound holds for
every eligible width and dimension. Selecting between this construction
and the retained $`N/b^2+n`$ construction gives the stronger
$`N/b^2+\min\{n,f(n)\}`$ statement. Each selection is one circuit
with all the stated count, depth, error, and workspace guarantees.

Choose $`C_\eta`$ so that $`\lceil C_\eta\sqrt N\rceil\ge17B_0`$
for every n. Such a constant exists because $`(L+n+7)/2^{n/2}`$
is bounded at fixed eta. Substitution proves the square-root-width
corollary. These constants are asymptotic guarantees, not practical
crossover estimates.

Finally, at physical width $`n+b+2=\Theta(b)`$, the
[inherited lower bounds](../../docs/UNIFORM_PRECISION_DEPTH.md#5-complete-the-uniform-theorem-and-its-matching-statements)
give

```math
T^\star=\Omega_\eta\!\left(\sqrt N+\frac Nb\right),\qquad
D_T^\star=\Omega_\eta(N/b^2).
```

If $`b\le\sqrt{N/f(n)}`$, then $`f(n)\le N/b^2`$ and
$`b\le\sqrt N`$. Both upper bounds therefore match these lower
bounds, proving the interval in Section 1. The unrestricted optimum
beyond that interval is not determined by this argument.

## 9. Scope and evidence

The proof uses three distinct ingredients: the new coherent descendant
source, the native shallow query schedule, and the retained linear-size
indicator for the late layers. The column encoding and robust recovery
above supply the missing complete-operator step. A state source alone
would not prove the theorem.

The result closes the previously stated $`o(n)`$ existence question at
fixed accuracy and sufficiently large constant square-root dirty width,
and supplies an all-eligible-width fixed-accuracy bound with a larger
matching interval. It leaves depth beyond that interval and a
general high-precision depth refinement open. The absolute-constant
extension in Section 10 covers its stated range of growing precision.
The [five bounded recovery checks](../../tests/test_sublinear_depth_recovery.py)
test signed and singular complete columns, marker tags, actual inverses,
active-only phase, and final leaked-work return. Their small numerical
unitary completions are verification devices, not compiler primitives.
The [five bounded calibration checks](../../tests/test_phase_gradient_calibration.py)
exercise the emitted signed modular-shift echo, literal rotation phase,
inactive occupied inputs, actual native preparation inverses, and gradient
leakage. They support that local native interface, not a general frame emitter.
They are not evidence for the asymptotic count absorption or a claim that
the resulting exponentially sized circuits have been emitted.

## 10. A precision-uniform extension through linear precision

There are absolute constants $`c\gt0`$ and C such that for every
$`n\ge1`$, integer $`6\le L\le cn`$, and integer
$`b\ge17(L+n+7)`$, every prescribed complete real Hopf frame on
$`N=2^n`$ dimensions has one circuit with two external clean flags,
at most b arbitrary dirty qubits, and

```math
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-L}.
```

Put

```math
\ell=\log_2(n+2),\qquad p=L+\ell,\qquad
F(n,L)=\min\{nL,(\sqrt n+L)\ell\}.
```

The same circuit satisfies

```math
T\le C\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
G\le CNL,
```

```math
D_T\le C\left(\frac{NL}{b^2}+F(n,L)\right).
```

All constants here are independent of n, L, b, and the prescribed frame.
The full initialized-isometry contract, literal phases, arbitrary dirty
inputs and references, terminating certified preprocessing, and absence
of measurements or supplied precision states are the same as in Section 1.
For general eta, choosing its ceiling precision L gives
$`2^{-L}\le\eta`$ whenever this L lies in the displayed range.

The inherited lower bounds give simultaneous worst-case matching on
every nonempty interval

```math
17(L+n+7)\le b\le\sqrt{\frac{NL}{F(n,L)}}:
\qquad
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

The T-count also matches in order at every eligible width. For a
sufficient absolute $`C_b`$, the allocation
$`b=\lceil C_b\sqrt{NL}\rceil`$ is eligible in every stated dimension
and gives

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad
D_T=O((\sqrt n+L)\log_2(n+2)).
```

In particular, $`L=O(\sqrt n)`$ retains the
$`O(\sqrt n\log(n+2))`$ depth order, and
$`L=o(n/\log(n+2))`$ gives $`D_T=o(n)`$ at this allocation.
For example, $`L=\lceil\sqrt n\rceil`$ gives error at most
$`2^{-\sqrt n}`$ while retaining that
$`O(\sqrt n\log(n+2))`$ depth, for sufficiently large n.
These statements remain within $`L\le cn`$. This is a sufficient
precision range, not a claim that the logical scratch allocation works
at arbitrary L. The fixed-accuracy theorem and its proof above remain
valid separately.

### 10.1 An absolute source and query ledger

Use the same recovery with $`R=2k+1\le C_R\sqrt n`$ source calls,
where $`C_R`$ is absolute, but set

```math
\delta=\frac{2^{-L}}{64R},\qquad
\rho=1+\lceil\log_2(1/\delta)\rceil\le C_\rho p,
\qquad u=n+\rho+2.
```

The constants in the [descendant source](COHERENT_DESCENDANT_LOADER.md#6-native-resource-interface)
are absolute when written in h and rho. Its geometric expansion uses
the fixed numbers $`3/5`$ and $`4/5`$; its length is $`O(\rho)`$.
Parallel gradient synthesis and its actual inverse have
$`O(\rho^2)`$ T and Clifford counts and $`O(\rho)`$ T-depth.
The helper-free modular adders have $`O(\rho)`$ count and depth.
The fixed coefficient bank and outer damping fit these bounds using
the same uniform one-qubit synthesis premise. Thus no source constant
depends separately on the requested error.

All simultaneously initialized internal work fits $`C_s p`$ bits
for an absolute $`C_s`$. This includes tag, source index, both
attenuations, calibration program, gradient core, and the additional
returned addition register used by the shared version below. In the
independently prepared source, the outer selective phase tests every
live work bit as in Section 5. The shared version has the different
reflection specified below. Both reflections are native zero predicates
with dirty returned helpers, not lookup tables indexed by the data and
every precision-work bit.

For $`Q=2^h`$, every large table has at most
$`C(h+1)Q\rho\le Cu^2Q`$ rows, at most $`Cu`$ output bits,
and at most $`Cu`$ address bits. Padded rows and the activity guard
only change the absolute constant. The conservative shallow-query
lemma therefore gives, per completed query,

```math
T,w\le Cu^5\sqrt Q,\qquad G\le Cu^5Q,
\qquad D_T\le C\log_2(u+2).
```

There are a constant number of such queries per source call, including
their actual inverses and the source reflection. The marker-to-tag query
has Q rows and $`O(\log(n+2))`$ outputs and fits the same ledger.
Short-tag masks have $`O(n)`$ rows and $`O(n)`$ outputs; applying
the query lemma to them costs at most $`Cu^5`$. Native zero predicates
on at most $`Cu`$ controls have count and returned dirty work at most
$`Cu^3`$, by the retained read-only conjunction bound. This also
charges the initial and final protected-suffix predicates. None of
these operations adds a new initialized register outside the protected
logical block.

Multiplying counts by $`O(\sqrt n)`$ recovery calls and bounding all
additional local words gives the following convenient degree-six bound
for the independently prepared version:

```math
T_E,w_E\le C_Eu^6(\sqrt Q+1),\qquad
G_E\le C_Eu^6(Q+1),\qquad
D_{T,E}\le C_E\sqrt n\,p.
```

Peak width is reused after completed return and need not be multiplied
by the number of calls; the displayed larger bound is sufficient.
Every constant in this ledger is absolute. In particular, the earlier
$`\mathrm{poly}_\eta(n)`$ notation is not being used to assert a
precision-uniform bound without controlling its degree or coefficients.

To obtain the stronger theorem, use the
[shared-gradient early recovery](SHARED_GRADIENT_RECOVERY.md) instead.
Its arithmetic precision is

```math
r=L+\lceil\log_2R\rceil+16
\le L+\tfrac12\log_2(n+2)+20=O(p).
```

The three r-bit registers are the program B, retained gradient F, and
returned zero addition register C. Their total allocation, including
the tag, source index and two attenuations, is $`3r+q+s+2\le C_s p`$
after enlarging the same absolute constant. The new modular adder uses
two exact queries on 2r address bits with r outputs, and a Clifford
swap. It has

```math
T,w=O((r+1)^4 2^r),\qquad
G=O(r4^r+(r+1)^4 2^r),\qquad D_T=O(\log(r+2)).
```

Its literal addition contract holds on C zero for arbitrary program
and gradient inputs; its actual inverse gives subtraction there. The
completed activity echo is identity on every inactive input, including
occupied C. Coefficient banks use one short phase table, calibration
retains its original table, and fixed outer damping is hardwired into
arithmetic tables so that it has the required full-operator contract
even on occupied B.

The local proof charges one preparation and one actual inverse at the
exterior boundaries. Its recovery takes place inside the exactly
invariant sector $`|\Gamma_r\rangle_F|0^r\rangle_C`$. Source-angle
comparison first uses the invariant B-zero subspace; source-isometry
completion then applies on the larger mode space allowing occupied B.
The physical outer reflection omits F and includes C. The reflection
of Section 5, which tests F zero, must not be used for this shared
construction. No gradient reset or source-return promise is introduced.

The standalone arithmetic tables add genuine precision-dependent costs.
Including them, the early circuit has absolute bounds

```math
T_E,w_E\le C_Eu^6(\sqrt Q+2^r+1),\qquad
G_E\le C_Eu^6(Q+4^r+1),
```

```math
D_{T,E}\le C_E\left(L+\log(r+2)
       +\sqrt n[\ell+\log(r+2)]\right).
```

Here $`r=O(\rho)`$, so all per-call polynomial factors, short
tables, and the single preparation pair fit the same degree-six
majorant. The source data and angular errors telescope across the R
calls, whereas preparation error is paid only at the two boundaries.
The local proof gives complete early error at most
$`9\cdot2^{-L}/64\lt2^{-L}/4`$, including all work and flag
leakage. This is an entire early-unitary lemma, not merely a replacement
state-source promise.

### 10.2 A protected suffix and an absolute dimension cutoff

Fix an absolute integer

```math
A\ge\max\{16,C_s+1\},\qquad c=\frac1{8A},\qquad
M=\left\lceil A\left[L+\log_2(n+L+2)\right]\right\rceil.
```

This ensures $`M\ge C_s p`$. For all sufficiently large n,
uniformly over $`L\le cn`$, it also gives $`M\le n/2`$ and
$`h=n-M\ge n/2`$. For example, it suffices that $`n\ge8`$ and

```math
\log_2(2n+2)\le\frac{n}{8A}.
```

Indeed $`L\le n/(8A)`$ and $`n+L+2\le2n+2`$ imply
$`M\le n/4+1\le n/2`$. The early/tail operator decomposition
in Section 2 holds for this M with no change to its all-column proof.

There is an absolute $`C_u`$ with
$`u\le C_u(n+L+2)`$. Hence

```math
u^6 2^{-M/2}
\le C_u^6(n+L+2)^{6-A/2}2^{-AL/2}
\le C_u^6(n+L+2)^{-2}.
```

The analogous Clifford factor obeys
$`u^6 2^{-M}\le C_u^6(n+L+2)^{6-A}`$.
Since $`L\le cn\le n`$, one has $`u^6=O(n^6)`$ with an
absolute constant. Increase an absolute dimension cutoff until also

```math
cn+\tfrac12\log_2(n+2)+20\le n/4.
```

Then $`r\le n/4`$. The standalone arithmetic width and T-count
are at most $`O(n^6 2^{n/4})`$, and their Clifford count is at most
$`O(n^6 2^{n/2})`$. These costs are suppressed by their own precision
bound, not by $`2^{-M}`$. Thus, for absolute $`C_1,C_2`$,

```math
w_E\le C_1 n^{-2}\sqrt N+C_2n^6(2^{n/4}+1).
```

Choose once and for all an absolute integer $`n_0`$ after which both
preceding cutoff inequalities and the following inequalities hold:

```math
C_1n^{-2}\sqrt N+C_2n^6(2^{n/4}+1)\le\frac12\sqrt{N/n},
```

```math
\sqrt{N/n}\ge\max\{4n+18,32c_I(n+3)\}.
```

Include also the arithmetic Clifford condition

```math
C_2n^6(2^{n/2}+1)\le N,
```

enlarging the same absolute $`C_2`$ to cover its native constant.

Increase $`C_2`$ if necessary to cover the cubic returned-helper pools
of the tail predicates as well. The first inequality then also ensures
that those polynomial-size pools fit the returned query workspace.
Such an $`n_0`$ exists because $`c\le1/128`$, the logarithmic
terms are sublinear, and exponential growth dominates every fixed
polynomial; it is independent of L, b, and the frame.

For $`n\ge n_0`$ and $`b\ge\sqrt{N/n}`$, all early dirty
work fits b. The same bounds give $`T_E=O(\sqrt N)`$ and
$`G_E=O(N)`$, and hence the desired $`O(\sqrt{NL})`$ and
$`O(NL)`$ budgets. Conditional source work fits M by construction.
Since $`r\le n/4`$, its early depth is
$`O(L+\sqrt n\,\ell)`$.
No accuracy-dependent crossover or supplied phase-state reservation
is used.

<a id="103-the-uniform-tail-and-error-composition"></a>

### 10.3 Share the terminal precision source and compose the errors

Apply the [shared binary terminal lemma](SHARED_BINARY_TERMINAL.md)
to the final M layers. It reserves
$`R_0=\min\{M,\lceil A_0\ell\rceil\}`$ bottom logical bits for
a sufficiently large absolute $`A_0`$. On their active zero sector,
one logarithmic-width capped binary source serves all remaining terminal
single layers. Each layer has a fixed sixteen-term matrix-unit residual.
Its coefficient query reads preserved prefix and label bits, so it
unloads exactly after the selected target operation.

A nonwrapping failure counter inserts the intermediate private-zero
projections without projecting or preparing the source again. Uniform
rare-mode rejection gives the same constant attenuation on locally
inactive layers. Native preparation and actual inverse are paid once
per complete stream; global normalization-two recovery uses only three
stream appearances. Coarse programs and arithmetic have separate,
exactly returned scratch, including on failed branches. The completed
middle circuit is identity on the original inactive sector with arbitrary
occupied work: the unconditional preparations cancel with their actual
inverses around its guarded interior.

The lemma's absolute width and count ledger fits
$`b\ge\sqrt{N/n}`$ after enlarging the same absolute $`n_0`$
from Section 10.2. Its deepest $`R_0`$ layers retain the capped
arbitrary-dirty source and nonuniform-indicator query allocation.
Together they give one terminal circuit with

```math
T_T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad
G_T=O(NL),\qquad
D_{T,T}=O\!\left(\frac{NL}{b^2}+p\ell\right).
```

Its full-isometry error is less than $`2^{-L}/2`$, including
the bottom-only case $`R_0=M`$. Early and terminal dirty pools
are used sequentially, so their maximum is the required allocation.
The shared-gradient early error is at most $`9\cdot2^{-L}/64`$.
Full-isometry telescoping gives total error less than $`2^{-L}`$:
actual later unitaries propagate earlier leakage, with no reset
assumption and with all dirty/reference correlations retained.

Adding the early depth gives

```math
D_T=O\!\left(\frac{NL}{b^2}+L+\sqrt n\,\ell+p\ell\right)
=O\!\left(\frac{NL}{b^2}+(\sqrt n+L)\ell\right).
```

Here $`p\ell=L\ell+\ell^2`$, $`\ell\ge1`$,
and $`\ell^2=O(\sqrt n\,\ell)`$ with an absolute constant.
The earlier [grouped terminal proof](GROUPED_TERMINAL_DEPTH.md) remains
a standalone three-halves-power bound and supplies the exact coarse
word-prefetch interface reused by the new terminal construction.

### 10.4 Every eligible width and a uniform finite fallback

For $`n\ge n_0`$ and $`17(L+n+7)\le b\lt\sqrt{N/n}`$,
use the retained precision-uniform compiler. It has
$`D_T=O(NL/b^2+nL)`$, and here $`nL\le NL/b^2`$.
Its count term $`nL`$ is absorbed into $`\sqrt{NL}`$ throughout
$`L\le cn\le n`$, because

```math
\frac{nL}{\sqrt{NL}}
=\frac{n\sqrt L}{2^{n/2}}
\le\frac{n^{3/2}}{2^{n/2}}=O(1).
```

This proves all desired bounds on the smaller-width branch without
changing its literal base allocation.

For $`n\lt n_0`$, use that same retained compiler at every eligible b.
Its constants are already absolute. Moreover,

```math
nL\le n_0L\ell,
```

so its depth satisfies the new upper bound after increasing one absolute
constant. This fallback is uniform in L and b; it is not a separate
finite-size constant for each accuracy. It also avoids claiming that
the new indicator fits at the literal minimum $`17(L+n+7)`$.

Finally, the retained compiler itself satisfies the theorem's T and G
bounds throughout this precision range and has depth
$`O(NL/b^2+nL)`$. Selecting between the two complete circuits gives
the minimum in the theorem. Every choice retains the full return
contract and all three resource bounds on one circuit.

At every eligible width the physical width is $`\Theta(b)`$, so the
inherited lower bounds are

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
D_T^\star=\Omega(NL/b^2).
```

Here $`L=O(\sqrt{NL})`$, so the T-count matches throughout the
stated range. If $`b\le\sqrt{NL/F(n,L)}`$, then
$`F(n,L)\le NL/b^2`$ and $`b\le\sqrt{NL}`$, since
$`F(n,L)\ge1`$. This proves the displayed simultaneous matching
interval. Choosing a sufficient absolute $`C_b`$ makes
$`\lceil C_b\sqrt{NL}\rceil\ge17(L+n+7)`$ for every stated
n and L: use $`L\ge6`$, $`L\le n`$, and boundedness of
$`(2n+7)/2^{n/2}`$. Substitution proves the square-root-width
corollaries, including sublinear depth when
$`L=o(n/\log(n+2))`$.

The [all-precision theorem](ALL_PRECISION_DEPTH.md) now combines this
low-precision branch with a whole-frame binary stream and a separately
proved high-precision schedule. It retains optimal counts at every
precision; this chapter remains the stronger low-precision proof home.

The linear-precision restriction pays for the early simultaneously
initialized program and gradient work inside the protected logical suffix.
The terminal source and history use only logarithmic initialized work. This
proof supplies no uniform extension to arbitrary precision, no optimal
depth claim beyond inherited matching lower bounds, and no efficient
classical preprocessing or general native-emitter guarantee.
The [bounded shared-gradient checks](../../tests/test_shared_gradient_recovery.py)
exercise exact addition, occupied inactive work, literal coefficient
phases, preparation paid once, occupied-program damping, and
wrong-reflection/inverse controls. The
[grouped-terminal checks](../../tests/test_grouped_terminal_depth.py)
exercise native coarse program interpretation, actual unloading, occupied
inactive work, and rounded group partitions.
The [shared-binary-terminal checks](../../tests/test_shared_binary_terminal.py)
exercise the addressed matrix-unit kernels, retained failure history,
final-only source return, and occupied inactive preparation boundaries.
These checks support the native interfaces; the uniform degree, cutoff,
and resource statements above are analytic proofs, not numerical
extrapolations.
