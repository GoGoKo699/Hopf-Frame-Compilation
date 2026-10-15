# Precision-uniform depth from subtree-projector decoding

[Projector decoder](SUBTREE_PROJECTOR_DECODER.md) · [One-clean baseline](ONE_CLEAN_DEPTH.md) · [Resource statement](../../docs/OPEN_PROBLEM.md)

**Local result R79.** The native subtree-projector decoder extends to
variable precision with absolute constants. Its all-work anchor queries
have exponential precision costs. A larger protected logical suffix
absorbs those costs, and the one-clean binary terminal stream avoids
paying separately for every remaining precision source. The resulting
all-precision theorem improves R72's additive depth while preserving
its literal widths, optimal T-count, and linear Clifford budget.

## 1. The all-precision theorem

Let $`n\ge1`$ and $`L\ge6`$ be integers, and define

```math
N=2^n,\qquad \eta_L=2^{-L},\qquad
\ell=\log_2(n+2),\qquad \lambda=\log_2(n+L+2),
```

```math
F_{79}(n,L)=\min\{nL,(L+\ell)\lambda\}.
```

For either prescribed operator W below, every indicated integer dirty
allocation admits one coherent Clifford+T circuit V:

| Prescribed family | Required arbitrary dirty width |
|---|---|
| Complete real Hopf frame | $`b\ge17(L+n+7)`$ |
| Complex magnitude frame $`D_\phi W_{\mathbb R}`$ | $`b\ge17(L+n+8)`$ |

There is one external initialized qubit. With its initialization
isometry $`J_1`$, the same circuit satisfies

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+F_{79}(n,L)\right).
```

All constants are absolute: they are independent of n, L, b, and the
prescribed angles and phases. The norm includes every logical column,
literal phases, work leakage and return, arbitrary dirty inputs, and
external references. Unused supplied wires remain untouched. G counts
elementary Clifford gates, including the arbitrary Clifford interlayers
allowed by the T-depth model. The stronger inherited clause is

```math
L\gt N/n^2\quad\Longrightarrow\quad
D_T=O\!\left(\frac{NL}{b^2}+L\right).
```

The inputs have R72's terminating certified evaluation interface,
including signed, unequal, and singular angles and specified phases at
zero-amplitude leaves. Classical preprocessing has no efficiency claim.
Every query, preparation, guard, actual inverse, and recovery is charged.
There is no QRAM, measurement, reset, supplied precision state, or free
coherent evaluator. The complex family is the prescribed magnitude
completion, not arbitrary complex-unitary synthesis.

For fixed L this recovers R78's $`O(N/b^2+\log^2(n+2))`$ depth.
More generally, for $`L=\Theta(\log_2(n+2))`$ and sufficiently
large n, at $`b=\Theta(\sqrt{NL})`$ it gives

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad D_T=O(\log^2(n+2)).
```

Here the error is inverse polynomial in n, the number of logical qubits;
this does not assert the same depth at inverse-polynomial error in N.
Section 7 gives the matching intervals and remaining limits.

## 2. The unchanged early decoder and its precision budget

Use the exact geometry and native words of
[R78, Sections 2–5](SUBTREE_PROJECTOR_DECODER.md#2-the-literal-nested-subtree-spaces).
For an h-bit early frame, its nested subtree scaling projectors obey

```math
A_d=\sum_{t=0}^d\Pi_t,\qquad
(I-2A_d)x_t=(-1)^{[t\le d]}x_t.
```

The complete tagged source prepares the literal frame column together
with its band label t. With $`s=\lceil\log_2(h+1)\rceil`$, a
reversible binary search using these projectors recovers t, XORs it
into the old tag, and reverses the search. R78's full-mode reflection
comparison uses the original per-call descendant loader, including
fresh local gradient preparation and its actual inverse. Its error
is at most $`4\delta`$ per reflection, also on occupied source work.
Its comparison reflection is identity off the transported all-zero-work
subspace; it is not a logical reflection tensored with identity on
arbitrary source work. Approximate projectors need not commute or nest.

Save $`H=[Z=0^M]`$ before using the protected logical suffix Z.
Only H is externally initialized. The tag, search record, threshold
record, and reused source bank all occupy Z on H one. On H zero
the completed circuit is exactly identity on every occupied input,
with literal phase +1. The final actual predicate inverse restores H
on the ideal output and propagates all actual leakage unitarily.

Set

```math
\delta_0=\delta=\frac{2^{-L}}{64(s+1)},\qquad
p=1+\lceil\log_2(1/\delta)\rceil
 =L+7+\lceil\log_2(s+1)\rceil.
```

The complete early error from R78 is

```math
\|V_EJ_1-J_1(W_E\otimes I_b)\|
\le\delta_0+8s\delta\le\eta_L/8.
```

This estimate is already uniform in precision. The new work is to
bound its native costs uniformly, instead of substituting L into
R78's accuracy-dependent asymptotic constants.

## 3. Absolute early resource bounds

The absolute [per-call ledger in R78, Section 6](SUBTREE_PROJECTOR_DECODER.md#6-native-resource-ledger-and-protected-block-choice)
has source bank size $`u=O(p)`$ and $`2^u\le C\delta^{-2}p`$.
The central reflection tests this entire bank. For a source or
reflected-source macro, that ledger gives

```math
\begin{aligned}
T,w&\le C(h+p+2)^6\delta^{-1}(\sqrt Q+1),\\
G&\le C(h+p+2)^6\delta^{-2}(Q+1),\\
D_T&\le C[\log_2(h+2)+p],\qquad Q=2^h.
\end{aligned}
```

Here and below C denotes an absolute constant which may increase.
In particular, the anchor query has up to
$`C(h+1)Qp\delta^{-2}`$ rows. Its precision cost is explicitly
included, even though its output is just one bit. The initial tagged
source, all $`4s`$ descendant-source appearances, $`2s`$ central
phases, threshold arithmetic, guarded XOR, and two saved-predicate
boundaries are included in the following complete early envelope.
Put $`K=n+L+2`$. Since $`p=O(L+\log_2(n+2))=O(K)`$
and $`s=O(\log_2(n+2))`$, a safe bound is

```math
\begin{aligned}
T_E,w_E&\le C K^{10}2^L(\sqrt Q+1),\\
G_E&\le C K^{10}2^{2L}(Q+1),\\
D_{T,E}&\le C(L+\ell)\ell.
\end{aligned}
\qquad\text{(1)}
```

For the first two lines, substitute
$`\delta^{-1}=64(s+1)2^L`$ into the per-call bounds and
multiply by the number of appearances. Degree ten dominates these
polynomials and the predicate boundaries. Width need only be maximized
over the returned query pools; the displayed sum bound is conservative.
For depth, retain $`s=O(\ell)`$ and $`p=O(L+\ell)`$ instead.
No shared-gradient invariant or additional clean source is used.

### Protected-bank choice and the standalone costs

Fix an absolute $`A\ge40`$, large enough also for the simultaneous
conditional clean-work bound $`3s+u+O(1)`$. Set

```math
c=\frac1{8A},\qquad
6\le L\le cn,\qquad
M=\lceil A[L+\lambda]\rceil,\qquad h=n-M.
```

Then $`K=\Theta(n+2)`$ on this range, with absolute comparison
constants. Choose an absolute dimension cutoff large enough that
$`M\le n/2`$ and $`L\le N/n^2`$ uniformly on the range.
For example, $`AL\le n/8`$ and
$`A\log_2(n+L+2)+1=o(n)`$ establish the first condition.
The conditional bank fits because $`3s+u+O(1)\le C_{\rm bank}(L+\lambda)`$
and A was chosen after that native constant.

The large-table terms in (1) satisfy

```math
K^{10}2^L\sqrt Q
\le\sqrt N\,2^{(1-A/2)L}K^{10-A/2}
\le\sqrt N K^{-10},
```

```math
K^{10}2^{2L}Q
\le N\,2^{(2-A)L}K^{10-A}
\le NK^{-30}.
```

The standalone terms in (1) are not polynomial in precision. They
are bounded by $`CK^{10}2^{cn}`$ and $`CK^{10}2^{2cn}`$.
Since $`c\lt1/2`$ and $`2c\lt1`$, exponential dominance gives,
uniformly for all $`6\le L\le cn`$,

```math
w_E=o(\sqrt N/n),\qquad T_E=O(\sqrt N),\qquad G_E=O(N).
\qquad\text{(2)}
```

All cutoffs here depend only on fixed native constants. In particular,
there is no hidden cutoff depending on L or on the error tolerance.

## 4. The remaining layers with one clean qubit

The new protected suffix may have $`M=\Theta(L+\ell)`$ bits.
Use the actual one-clean terminal construction from
[R72, Sections 3–5](ONE_CLEAN_DEPTH.md#the-final-logical-layers-at-smaller-error),
including its relocated padding-rejection bit and its uncapped bottom
layers. The original two-clean terminal lemma alone would not justify
the clean allocation or the smaller width boundary below.

After fixing A, take the terminal construction's absolute workspace
constant $`C_Z`$, choose $`A_0\ge\max\{24,C_Z+3\}`$, and set

```math
R_0=\min\{M,\lceil A_0\ell\rceil\},\qquad R=M-R_0.
```

If R is positive, a shared binary source serves these R middle layers.
All its initialized source, padding, private failure, retained history,
arithmetic, and coarse-program fields, including a distinct relocated
padding-rejection coordinate, fit in the bottom $`R_0`$ logical bits
on the saved active sector. Only the same external H is used.

Enlarging A changes neither the stream's source interface nor its
native polynomial degree: its cap and stage count remain
$`J,R=O(L+\ell)=O(n)`$. The
[proved terminal ledger](SHARED_BINARY_TERMINAL.md#7-simultaneous-resources)
and R72's relocation give

```math
\begin{aligned}
T_M,w_M&\le C(n+2)^{10}(\sqrt N\,2^{-R_0/2}+1),\\
G_M&\le C(n+2)^{10}(N2^{-R_0}+1),\\
D_{T,M}&=O((L+\ell)\ell).
\end{aligned}
```

Since $`R\gt0`$ implies $`R_0=\lceil A_0\ell\rceil`$,
these bounds imply $`w_M=o(\sqrt N/n)`$, $`T_M=O(\sqrt N)`$,
and $`G_M=O(N)`$ after one absolute cutoff. The full stream error
is at most $`3\eta_L/256`$. Its actual preparation words cancel
on H zero for arbitrary occupied work. On H one its zero-work
reflection includes the relocated rejection bit; the actual recovery
and final predicate inverse propagate leakage without reset.

For the deepest $`R_0`$ layers, R72's smaller-error bottom has

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad Q_j=N2^{-j}.
```

Its dirty base is $`B'_0=L+n+9`$. The original supplied real
threshold leaves $`B=b-B'_0\ge b/2`$. Whenever
$`b/2\ge16c_I(n+2)`$, the shallow rectangular queries fit and

```math
\begin{aligned}
T_B&=O(\sqrt{NL}+NL/b+R_0L+n^2),\\
G_B&=O(NL),\\
D_{T,B}&=O(NL/b^2+R_0[L+R_0+\ell]).
\end{aligned}
```

Its error is less than $`15\eta_L/64`$. Since $`R_0=O(\ell)`$,
the additive depth is $`O((L+\ell)\ell)`$. The count overhead
is absorbed by $`R_0L\le nL\le\sqrt{NL}`$ and
$`n^2=O(\sqrt N)`$. If R is zero, omit the middle stream and
apply this same bottom estimate to $`R_0=M=O(\ell)`$.

For any two successive actual factors and their ideal logical/dirty
operators, use the full initialized-isometry identity

```math
V_2V_1J_1-J_1R_2R_1
=V_2(V_1J_1-J_1R_1)+(V_2J_1-J_1R_2)R_1.
```

Thus the same H can save the early M-bit predicate, then the stream's
$`R_0`$-bit predicate, and then serve as the bottom signal. Only
ideal comparison terms have exact clean return; later actual unitaries
propagate prior leakage. The complete error is at most

```math
\eta_L/8+3\eta_L/256+15\eta_L/64
=95\eta_L/256\lt\eta_L.
```

This argument includes arbitrary dirty/reference inputs and preserves
literal phases. Sequential dirty reservations are maximized, not added.

## 5. Every width and precision

Use R72's width boundary

```math
b_c=\sqrt{\frac{NL}{n(L+n)}},\qquad
b_c\ge\sqrt{6/7}\,\frac{\sqrt N}{n}.
```

The second inequality follows from $`nL/(L+n)\ge6/7`$ for
$`n\ge1,L\ge6`$. Fix a single absolute $`n_0`$ meeting all
previous cutoffs and, uniformly on $`6\le L\le cn`$,

```math
w_E,w_M\le b_c/4,\qquad
b_c/2\ge16c_I(n+2).
```

Equations (1)–(2) and the terminal bounds prove such a cutoff exists.
Returned polynomial helper pools are included in these reservations.
For $`n\ge n_0`$ and $`6\le L\le cn`$, select one complete
real circuit as follows:

- If $`b\ge b_c`$, use Sections 2–4. It has simultaneous
  $`T=O(\sqrt{NL}+NL/b)`$, $`G=O(NL)`$, and
  $`D_T=O(NL/b^2+(L+\ell)\ell)`$ at the literal threshold.
- If $`b\lt b_c`$, use
  [R72's all-layer one-clean fallback](ONE_CLEAN_DEPTH.md#3-a-one-clean-fallback-at-every-width).
  Its additive depth $`n(L+n)`$ is less than $`NL/b^2`$.
  Its extra $`nL`$ count is at most $`\sqrt{NL}`$ on this
  low-precision range. Its full error is less than $`15\eta_L/16`$.

The refined branch has the claimed minimum, since

```math
(L+\ell)\ell\le(L+\ell)\lambda,\qquad
\frac{(L+\ell)\ell}{nL}
=\frac\ell n+\frac{\ell^2}{nL}=O(1).
```

Both ratios are uniformly bounded for $`n\ge1,L\ge6`$.

For all other n and L, invoke the complete R72 theorem, with
$`F_{72}(n,L)=\min\{nL,(\sqrt n+L)\lambda\}`$.
If $`L\gt cn`$, then

```math
\sqrt n+L\le(1+c^{-1})L\le(1+c^{-1})(L+\ell).
```

If $`n\lt n_0`$, uniformly for unbounded $`L\ge6`$,

```math
\sqrt n+L\le(1+\sqrt{n_0}/6)L
\le(1+\sqrt{n_0}/6)(L+\ell).
```

In either case $`F_{72}=O(F_{79})`$, by
$`\min\{a,Cb\}\le C\min\{a,b\}`$ for $`C\ge1`$.
R72 selects one actual native circuit supplying all the stated bounds
at the original width, including its bounded-dimension all-precision
fallback. No finite-accuracy net or mixing of resource bounds from
different circuits is used. Whenever $`L\gt N/n^2`$, select
R72's stronger high-precision circuit; it also satisfies the general
claim because $`F_{79}\ge L`$. This proves the real theorem.

## 6. Literal complex phases

At $`b\ge17(L+n+8)`$, select the new real theorem at precision
$`k=L+1`$ and the
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at the same k. Both fit exactly this threshold. The full-isometry
hybrid reuses the same sole external signal and dirty pool; the two
errors add to at most $`2^{-L}`$. All leaf phases, including those
at zero-amplitude leaves, remain prescribed literally.

The diagonal has

```math
T_D=O(\sqrt{Nk}+k+Nk/b),\qquad G_D=O(Nk),
```

```math
D_{T,D}=O(Nk/b^2+k+\ell).
```

Uniformly in n and L,

```math
F_{79}(n,L+1)=O(F_{79}(n,L)),\qquad
F_{79}(n,L)\ge L,\qquad F_{79}(n,L)\ge\ell.
```

The first follows by bounding each entry of the minimum under the
unit precision shift. For the others each entry separately dominates
L and ell; in particular $`nL\ge\ell`$ for $`n\ge1,L\ge6`$.
Thus the half-error shift and diagonal overhead fit the theorem.
The real circuit is reselected at k, even if this crosses a branch
boundary. The stronger high-precision clause may be inherited directly
from R72's complex theorem at L, with all the same simultaneous bounds.

## 7. Matching intervals, evidence, and remaining limits

The [inherited complete-frame lower bounds](ONE_CLEAN_DEPTH.md#8-matching-intervals-and-evidence)
give

```math
T^\star=\Omega(\sqrt{NL}+L+NL/b),\qquad
D_T^\star=\Omega(NL/b^2),
```

because the physical width $`n+b+1=\Theta(b)`$ at the literal
thresholds. Therefore T is worst-case optimal at every eligible width.
For the family-specific $`B_0=L+n+7`$ or $`L+n+8`$, every
nonempty integer interval

```math
17B_0\le b\le\sqrt{NL/F_{79}(n,L)}
```

has $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$. Nonemptiness means
$`17B_0\le\lfloor\sqrt{NL/F_{79}}\rfloor`$. The depth term
is absorbed by the width term, and $`F_{79}\ge L`$ makes
$`NL/b`$ absorb both other count terms. The stronger high-precision
clause also retains R72's matching interval through $`b\le\sqrt N`$
when that interval is eligible.

This proof changes precision allocation, protected-bank size, and
circuit selection. It introduces no new circuit identity or external
premise. The failure-sensitive
[projector-decoder checks](../../tests/test_subtree_projector_decoder.py)
remain the evidence for the unchanged nested-space, phase, and
occupied-work algebra; the source and terminal evidence is linked in
their proof homes. The exponential-versus-polynomial estimates and
the all-precision splice are analytic proofs, not conclusions from
finite tests. Immutable native receipts are unchanged.

R78 remains the original fixed-accuracy construction and proof of the
decoder. R72 remains the all-precision baseline and supplies the
one-clean terminal, fallback, diagonal, and high-precision circuits.
R66's older two-clean linear-precision resource bound is superseded,
while its source and terminal proofs remain dependencies. No
unrestricted growing-depth lower bound follows at square-root width.
The remaining constant-layer obstructions, R67's overlap-encoding
limit, practical crossover, efficient preprocessing, and general native
emitter questions retain their previous scopes. Internal audits are
not external technical review or a priority claim.
