# Linear-precision depth at every accuracy

[Polynomial-precision bridge](POLYNOMIAL_PRECISION_DEPTH.md) ·
[Direct prefix-programmed terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md) ·
[Global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md) ·
[Retained all-precision R79 theorem](UNIFORM_SUBTREE_DEPTH.md)

**Theorem R96.** The small-workspace
opening and a directly prefix-programmed high-precision terminal fit together
at every accuracy. The terminal queries the actual data of its local native
word instead of enumerating every fine-precision program. At sufficiently
large precision the original frame already has constantly many groups,
so no early contraction is required.

## 1. The all-precision statement

For integers $`n\ge1`$ and $`L\ge6`$, put
$`N=2^n`$, $`\ell=\log_2(n+2)`$ and $`\eta=2^{-L}`$.
Every prescribed complete real Hopf frame W and every integer
$`b\ge17(L+n+7)`$ admit one native Clifford+T circuit V with one
external initialized qubit and at most b arbitrary dirty qubits such that

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+L+\ell\right).              \qquad(1)
```

The prescribed complex magnitude completion satisfies the same bounds at
$`b\ge17(L+n+8)`$. All constants are absolute, independent of n,L,b
and the prescribed angle and phase data. These resources belong to one
actual circuit in every branch. G counts elementary Clifford gates.

The norm includes every logical column, literal signs and phases, singular
continuations, all work leakage and return, arbitrary dirty inputs, and
external references. Classical preprocessing terminates but need not be
efficient. There is no QRAM, measurement, reset, supplied precision state,
additional external clean word or free coherent evaluator. Every source,
query, guard, recovery, actual inverse and simultaneous workspace is paid.

The inherited lower bounds match the T-count. The depth bound is an upper
bound; optimal unrestricted depth remains open. The complex family is the
prescribed magnitude completion, not arbitrary complex-unitary synthesis.
No practical crossover or scalable emitter follows.

## 2. The two native component interfaces

### Fast early stage

Use the [global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md)
with its [fast geometric references](FAST_CAPPED_GEOMETRIC_REFERENCE.md).
For remaining height m, unchanged original prefix of length n-m and
stage precision l, put

```math
K=m+l+2,\qquad
\chi=\lceil\log_2(l+2)\rceil+\lceil\log_2(m+2)\rceil^3,\qquad
r=\lceil A\chi\rceil.                                  \qquad(2)
```

When $`r\lt m/2`$, this compiles the earlier m-r layers with full
initialized-isometry error at most $`2^{-l}/16`$, leaving the original
bottom r coordinates. Its saved-predicate word uses the sole external
clean bit; the completed body is literal identity on the inactive sector,
including occupied conditional work. Its sharp standalone count and dirty
width are at most $`S_{\rm loc}=K^{C_s}2^{C_s\chi}`$. With
$`\overline Q=2^{n-r+C_q\chi}`$, its native ledger is

```math
\begin{aligned}
T_E&\le CK^d(\sqrt{\overline Q}+\overline Q/b+1)+S_{\rm loc},\\
G_E&\le CK^d(\overline Q+1)+S_{\rm loc},\\
D_E&\le C\left[l+\log_2(m+2)+\log_2(n+\chi+2)
                         +K^d\overline Q/b^2\right].
\end{aligned}                                               \qquad(3)
```

The sufficient simultaneous allocation is

```math
b\ge8W_{\rm loc}+CK^{d_w}(n+\chi+2),\qquad
W_{\rm loc}\le S_{\rm loc}.                             \qquad(4)
```

All coefficients, degrees and bank constants are fixed before A.
Since $`\chi\ge\log_2K`$, choose A so that

```math
r-C_q\chi\ge(2d+80)\chi.                                \qquad(5)
```

The large-table terms in (3) then fit $`O(\sqrt N+N/b)`$ T gates,
$`O(N)`$ Clifford gates and $`O(N/b^2)`$ depth. Every long table
retains its actual global prefix, including data-dependent catalogues and
source anchors. Dirty precision cores and unused protected padding are
excluded from initialized-work anchor tests.

### Direct high-precision terminal

The [direct-prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md)
acts on h consecutive logical bits X, preserves the upper prefix P of
length $`n-s-h`$, and has at least three original lower suffix bits Z.
For each group it separately saves $`H=[Z=0^s]`$ in the sole external
clean bit. On that sector, three distinct original Z bits supply the
local one-clean flag d, native instruction selector z and optional
finite-bank presence selector w. All payload-gate targets exclude P,H,z,w;
w is preserved throughout its optional bank and actually unqueried.
Every completed group is literal identity on H zero, including occupied
d,z,w, and the actual predicate inverse is applied last. Its full
initialized-isometry error is at most $`2^{-\kappa}`$.

The terminal proves a common high-precision native grammar across P,
including identical live precision layouts and actual row/column replay.
It prefixes complete fine masks and programs only a separately bounded
coarse native remainder. The fixed precision-sized nonquery words have
elementary length $`O(\kappa+P(h))`$, including Clifford gates and
exact bounded controls. The programmable support excludes the long
precision core, and its universal instruction schedule has at most
$`2^{3h}P(h)`$ Boolean positions. Here P is one fixed polynomial,
not a precision-dependent choice.

Use this terminal only when $`L\ge n^{1/4}`$, $`\kappa=L+O(1)`$ and

```math
1\le h\le g=\max\{1,\lfloor\log_2(L+2)/32\rfloor\}.       \qquad(6)
```

Writing $`N_s=2^{n-s}`$, its actual global fine-mask product sums are
$`O(N_s\kappa)`$ and $`O(\sqrt{N_s\kappa})`$. The coarse schedule
has corresponding factors bounded by
$`C\sqrt{N_s}\,2^{5h/2}P(h)`$ and $`CN_s2^{2h}P(h)`$.
Equation (6) absorbs these into $`O(\sqrt{N_sL})`$ and $`O(N_sL)`$.
Its serial prefix-query depth is $`O(2^{3h}P(h)\ell)=O(L)`$:
$`L\ge n^{1/4}`$ gives $`\ell=O(\log_2(L+2))`$, and the remaining
fixed powers of logarithms are absorbed by the strict power slack.
R79's local high-precision condition holds there; bounded local heights
use directly prefixed R72 layerwise words.

The saved full-suffix predicate has an exact Boolean query with at most
$`2^s\le N`$ rows. Including it and its actual inverse gives

```math
\begin{aligned}
T_i&=O(\sqrt{N_sL}+L+N_sL/b+\sqrt N+N/b),\\
G_i&=O(N_sL+N),\\
D_i&=O(L+\ell+N_sL/b^2+N/b^2).
\end{aligned}                                               \qquad(7)
```

The terminal's physical nonquery peak is at most $`2\kappa+h+C`$
after an absolute cutoff: its live precision base has coefficient one
in kappa, and the additional coarse support and helpers are at most kappa.
With an exclusive $`B_Q=\lfloor b/2\rfloor`$ pool, sufficient conditions are

```math
2\kappa+h+C\le b/4,\qquad B_Q\ge C_Q(n+2).               \qquad(8)
```

The three original suffix bits are separate from this dirty allocation.
No full local $`17\kappa`$ pool is reserved before adding the global pool.
The predicate query uses the same returned pool outside the local body.

## 3. Fix constants before selecting dimensions

First fix the fast-stage constants, A, all instance degrees and all
standalone and minimum-pool coefficients in (2)–(5). Fix the terminal's
polynomial degrees, common native grammar and factor 1/32 in (6).
Then choose one absolute $`0\lt c\lt1/16`$ sufficiently small.

For $`L\le2^{cn}`$ and $`m\le n`$, fixed absolute coefficients give

```math
\log_2K\le cn+C_1\ell,\qquad
\chi\le cn+C_2\ell^3,\qquad
\log_2S_{\rm loc}\le C_3cn+C_3\ell^3,
```

```math
\log_2\!\left(CK^{d_w}(n+\chi+2)\right)
 \le C_4cn+C_4\ell^3.                                  \qquad(9)
```

Choose c so that $`Ac\lt1/8`$, $`C_3c\lt1/16`$ and
$`C_4c\lt1/16`$. The terms $`\ell^3=o(n)`$ are absorbed beyond
one absolute dimension cutoff. In this exponential-precision range,
fixed powers of K have not been silently described as polynomial in n.

The height estimates below determine a uniform group-count bound
$`R_*`$, depending only on these fixed constants and c. After fixing
that bound, choose

```math
\kappa=L+\lceil\log_2(R+1)\rceil+10.                    \qquad(10)
```

Finally choose one absolute $`n_0`$ covering all dominance estimates,
bounded local cases, precision padding and width conditions. It may
depend on c and $`R_*`$, but not on L,b or the target. This order
avoids a circular precision or reserve choice.

## 4. Moderate precision: three early stages and a short terminal

Work on

```math
n\ge n_0,\qquad b\ge\sqrt N/n,\qquad
n^{1/4}\le L\le2^{cn}.                                  \qquad(11)
```

Make the retained finite whole-family magnitude/sign reduction at error
$`\eta/512`$. Its two exact sign phases have full n-bit addresses and
cost $`O(\sqrt N+N/b)`$ T gates, $`O(N)`$ Clifford gates and
$`O(N/b^2+\ell)`$ depth. Singular continuations remain prescribed;
subsequent factors refer to the resulting positive tuple.

Starting at $`m_0=n`$, use at most three fast stages at $`l=L+5`$.
When $`r\lt m/2`$, execute and retain the original r-bit suffix;
otherwise omit and retain $`m\le2r`$. A height at most three goes
straight to the bottom. With $`a=\log_2(L+2)`$,

```math
m_1\le C(a+\ell^3),\qquad
m_2\le C[a+\log_2^3(a+\ell^3+2)],\qquad
m_3\le C[a+\log_2^3(m_2+2)]=O(a).                       \qquad(12)
```

Indeed (11) gives $`a\ge c_0\ell`$ beyond one absolute cutoff.
Hence $`m_1=O(a^3)`$, and $`\log_2^3(a+2)=O(a)`$ proves
$`m_2=O(a)`$ and then $`m_3=O(a)`$. Failed reserve tests give these
same bounds. The choice $`Ac\lt1/8`$ also makes the first reserve
proper beyond $`n_0`$.

The bounds (9) and the choice of c give uniformly

```math
S_{\rm loc}\le2^{n/8},\qquad
CK^{d_w}(n+\chi+2)\le2^{n/8}.                           \qquad(13)
```

Enlarge $`n_0`$ so $`W_{\rm loc}\le b/16`$ and the second term is
at most b/2. This proves the actual allocation (4). The reserve margin
(5) suppresses every large-table polynomial factor. The complete opening,
including sign phases and actual inverses, satisfies

```math
T_{\rm open}=O(\sqrt N+N/b),\qquad G_{\rm open}=O(N),\qquad
D_{\rm open}=O(N/b^2+L+\ell).                            \qquad(14)
```

Each stage contributes at most $`\eta/512`$ error. Returned allocations
are reused sequentially; previous leakage propagates through later full
unitaries and is not reset.

Write $`M=m_3=O(\log_2(L+2))`$. Keep the final three original bits for
the bottom and divide the remaining layers into groups of height at most
g, with one possible short first group. Their number has an absolute
bound $`R_*`$. Use (10). Every group saves and reverses its own entire
current lower-suffix predicate, and freshly queries its preserved prefix.

Here $`L\le N^c=o(\sqrt N/n)`$. Thus the terminal's precision base and
coarse support are small compared with the wide allocation, and (8) holds.
The query pool is disjoint from every live source bank, logical coordinate,
output, core and helper. This prices one actual simultaneous allocation.

## 5. Larger precision needs no early contraction

First use R79 directly whenever $`L\gt N/n^2`$. Its retained
high-precision clause already proves (1) at the original threshold.
The remaining larger-precision branch is

```math
n\ge n_0,\qquad b\ge\sqrt N/n,\qquad
2^{cn}\lt L\le N/n^2.                                  \qquad(15)
```

Make the same finite magnitude/sign reduction and divide all $`n-3`$
layers above the original three bottom bits directly into groups of height
at most g. No fast stage is used, so its polynomial-in-L standalone cost
is absent. Since $`\log_2L\gt cn`$,

```math
R\le1+n/g\le C/c,                                      \qquad(16)
```

an absolute bound. Also $`h\le g\le n/32`$ after the common cutoff.
The short first group keeps its actual height, including the terminal's
bounded-height fallback when needed.

The literal width pays for precision directly: $`b\ge17(L+n+7)`$.
The additional coarse support obeys $`2^hP(h)=o(L)`$ at (6).
Using $`\kappa=L+O(1)`$ and $`h\le n`$, (8) follows from

```math
2\kappa+h+C\le\frac{17}{4}(L+n+7)\le b/4               \qquad(17)
```

beyond the absolute cutoff. The original d,z,w bits do not count as dirty
helpers. The exclusive $`B_Q=\lfloor b/2\rfloor`$ pool is disjoint from
all live output, core, signal, metadata and coarse work. Its minimum
$`C_Q(n+2)`$ follows from the wide branch; the literal $`17L`$ term
also suffices on (15). The same returned pool implements the saved-suffix
query outside the local body.

This is a physical peak bound. It does not assume that precision work is
$`N^{o(1)}`$ in (15), and it does not add a separate local $`17\kappa`$
pool. Each group satisfies (7) at the original global allocation.

## 6. Counts, depth, bottom and full error

Both new branches have at most $`R_*`$ groups. For group i let h_i be
its actual height and $`s_i=3+(R-i)g`$ its lower suffix size. Its
preserved prefix has $`n-s_i-h_i`$ bits. The fine-mask products include
that entire factor, giving

```math
\sum_iN_{s_i}L=O(NL),\qquad
\sum_i\sqrt{N_{s_i}L}=O(\sqrt{NL}),\qquad
N_{s_i}=N2^{-s_i}.                                    \qquad(18)
```

The short first group's local $`2^{h_i}`$ factor cancels its prefix
factor. The separate full-suffix predicate queries in (7) cost at most
$`O(\sqrt N+N/b)`$, $`O(N)`$ and $`O(\ell+N/b^2)`$ per group.
Their constant-number sum fits (1), as does $`R\kappa=O(L)`$.
Equations (7), (14) and (18) therefore give the opening and middle bounds
on the same circuit.

The deepest $`\min(M,3)`$ layers after an opening, or the final three
layers in (15), use the retained R72 smaller-error bottom. It has
$`q_j=L+j+7`$, $`m_j=L+j+8`$ and actual global rows $`Q_j=N2^{-j}`$
for the relevant $`j\in\{1,2,3\}`$. Its base $`L+n+9`$ and a separate
$`\Theta(b)`$ pool fit these large wide branches. At most three native
sources and their long-prefix masks give

```math
T_{\rm bot}=O(\sqrt{NL}+L+NL/b),\qquad G_{\rm bot}=O(NL),\qquad
D_{\rm bot}=O(NL/b^2+L+\ell).                           \qquad(19)
```

The retained error is below $`15\eta/64`$. Its L depth term is kept
explicitly, rather than inherited as a smaller-range $`O(\ell)`$ shorthand.

By (10), the group errors sum to at most $`R2^{-\kappa}\le\eta/1024`$.
Every selected local program is emitted exactly, including its rejected
operator and actual inverse; programming introduces no approximation.
Coherent prefix sectors take a maximum, not a sum over their number.
The complete moderate-branch real error is at most

```math
\frac{\eta}{512}+\frac{3\eta}{512}
 +\frac{\eta}{1024}+\frac{15\eta}{64}
 =\frac{249\eta}{1024}\lt\frac\eta4.                     \qquad(20)
```

The larger branch omits the three early-stage errors. Each saved external
predicate is inverted on the ideal returned suffix. Its actual inverse
propagates prior leakage unitarily. Full initialized-isometry telescoping
justifies reuse of the external bit and returned pools, without assuming
exact return of approximate local cores.

## 7. Every dimension, width and complex phase

Use the independently proved [polynomial-precision theorem](POLYNOMIAL_PRECISION_DEPTH.md)
for $`6\le L\le n^{3/4}`$, and original-width R79 for $`L\gt N/n^2`$.
Above the common dimension cutoff, every other precision requiring a new
circuit exceeds $`n^{1/4}`$, so (11) and (15) cover it. Possible empty
ranges create no gap.

For $`b\lt\sqrt N/n`$, use whole R79 at the original precision and width:
its additive depth is absorbed because

```math
F_{79}(n,L)\le nL\le n^2L\lt NL/b^2.                    \qquad(21)
```

For $`n\lt n_0`$, use that same original-precision compiler at every
width and every unbounded L. Then $`nL\le n_0L=O(L+\ell)`$ with one
absolute constant. No stronger precision is imposed on a real fallback at
its literal minimum. This proves the real theorem at all parameters.

For complex magnitude frames, apply the real theorem at L+1 and then the
retained one-clean literal diagonal at L+1. The real threshold is exactly

```math
17((L+1)+n+7)=17(L+n+8),                                \qquad(22)
```

also sufficient for the diagonal. Both costs obey (1) because
$`L+1=\Theta(L)`$. Each error is at most $`\eta/2`$; initialized-isometry
telescoping gives at most eta while reusing the single initialized qubit.
Prescribed phases at zero-amplitude leaves remain literal. Narrow and
bounded-dimension complex cases may equivalently use whole original-
precision R79 without a separate diagonal.

Every branch selects one actual circuit with simultaneous T-count,
Clifford count, depth, width and full error. Sequential allocations are
reused only through their completed native words. This establishes (1).

## 8. Dependency and scope

The low-precision component is proved independently using R79 local words
inside R94's complete-workspace catalogue. The new terminal uses the actual
R69/R72 high-precision grammar, prefixes its paid masks and programs its
separately bounded coarse remainder. Its common fine skeleton, coefficient-
one precision base, exact bounded controls and full inactive behavior have
their own proof home. The global stage and fast references do not invoke
this theorem, so the dependency graph is acyclic.

The all-precision result changes an upper bound. It does not prove that
the additive $`L+\ell`$ term is unavoidable. Existing constant-layer
obstructions and architecture-specific failures retain their quantifiers;
optimal general depth remains a separate destination.
