# An absolute precision envelope for dirty catalogue compilation

[Catalogue terminal R90](DIRTY_CATALOGUE_TERMINAL.md) ·
[Uniform coarse stage R89](CAUSAL_COMPACT_COIN_DEPTH.md#10-prefix-selected-recursion-and-iterated-logarithm-depth) ·
[All-precision theorem R79](UNIFORM_SUBTREE_DEPTH.md)

**Local result R92.** The R90 construction admits an explicit absolute
precision cost. Three coarse stages and its unchanged dirty catalogue
give logarithmic depth on an unbounded, slowly growing precision range.
This chapter removes the accuracy-dependent constants from that part of
the general theorem. It introduces no new quantum primitive and does not
improve the inverse-polynomial-error benchmark by itself.

## 1. The all-parameter statement

Let $`n\ge1,L\ge6`$ be integers and put

```math
N=2^n,\qquad \eta=2^{-L},\qquad \ell=\log_2(n+2),
```

```math
F_{79}(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

There is an absolute constant $`C_*`$ such that every prescribed operator
W in the following table admits a coherent Clifford+T circuit V at every
indicated integer dirty allocation.

| Prescribed family | Arbitrary dirty width |
|---|---|
| Complete real Hopf frame | $`b\ge17(L+n+7)`$ |
| Complex magnitude frame $`D_\phi W_{\mathbb R}`$ | $`b\ge17(L+n+8)`$ |

One external initialized qubit suffices. With its initialization isometry
$`J_1`$, the same circuit satisfies

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad G=O(NL),
```

```math
D_T=O\!\left(\frac{NL}{b^2}+
\min\{F_{79}(n,L),\ell+2^{C_*(L+1)^3}\}\right).       \qquad(1)
```

Every constant is absolute, independent of n, L, b and the prescribed
angles and phases. G includes elementary Clifford gates. The norm covers
every logical column, literal phases, work leakage and return, arbitrary
dirty inputs, and external references. Signed and singular data retain
the terminating certified evaluation interface of R72. The prescribed
complex completion includes phases at zero-amplitude leaves.

The stronger retained high-precision clause remains

```math
L\gt N/n^2\quad\Longrightarrow\quad
D_T=O(NL/b^2+L),                                      \qquad(2)
```

with the same simultaneous guarantees. No QRAM, measurement, reset,
supplied precision state, additional initialized word, or free evaluator
is used. Classical preprocessing terminates but need not be efficient.

The proof first establishes the alternative term in (1) for real frames
at all parameters, then takes its minimum with R79. A unit precision shift
supplies the literal complex extension. Section 9 gives the growing-
precision corollary.

## 2. A uniform stopping rule for the existing coarse stage

Use only the explicit local-stage ledger in
[R89 Sections 10.1–10.3](CAUSAL_COMPACT_COIN_DEPTH.md#101-a-local-stage-and-catalogues-selected-by-an-immutable-prefix),
whose constants are uniform in the local precision. Its fixed-L recursive
theorem and unspecified fixed-L stopping height are not needed here.
For remaining height m define

```math
\begin{aligned}
l&=L+\lceil\log_2(m+2)\rceil+10,\\
K&=m+l+2,\qquad u=l+\lceil\log_2K\rceil,\\
r_L(m)&=\lceil A_0u^3\rceil.
\end{aligned}                                         \qquad(3)
```

Fix the absolute reserve coefficient $`A_0`$ after the native constants
$`C_0,d_0,d`$ and conditional clean-work reservation in R89. It can ensure

```math
r_L(m)-C_0u^3\ge(2d+80)u.                             \qquad(4)
```

All logical source, mask, metadata, reference, and donor coordinates fit
inside this reserve. The R89 proof also bounds every additional local
dirty field, routing array and returned helper by $`CK^e`$ for one fixed
exponent e. Its exponential precision-dependent row factors are included
in the address allowance $`C_0u^3`$, rather than in that polynomial width.

Write $`p=L+1\ge7`$. Ceilings and
$`K\le C(m+2)(L+1)`$ give an absolute $`C_1`$ with

```math
r_L(m)\le C_1[p+\log_2(m+2)]^3                       \qquad(5)
```

for every $`m\ge1`$. Choose an absolute $`m_0`$ such that

```math
4C_1\log_2^3(m+2)\le m/4\qquad(m\ge m_0).
```

Choose an absolute integer $`B\ge16C_1`$ with $`B7^3\ge m_0`$,
and set

```math
H(L)=\lceil Bp^3\rceil.
```

For $`m\gt H(L)`$, the inequality
$`(x+y)^3\le4(x^3+y^3)`$ proves

```math
r_L(m)\le4C_1p^3+4C_1\log_2^3(m+2)\le m/2.           \qquad(6)
```

Starting at $`m_0'=n`$, execute the stage only while its current height
exceeds H(L), replacing m by $`r_L(m)`$. Stop after at most three stages
and finish with the R90 catalogue. Equation (6) supplies a legitimate
proper protected suffix at every executed stage and prevents a growing
reserve from invalidating the factorization. It is a uniform stopping
rule of order $`(L+1)^3`$.

If all three stages occur, (5) gives, with absolute constants,

```math
\begin{aligned}
m_1&\le C_2[p+\ell]^3,\\
m_2&\le C_3[p+\log_2(\ell+2)]^3,\\
m_3&\le C_4[p+\log_2\log_2(\ell+4)]^3.
\end{aligned}                                         \qquad(7)
```

For each substitution use
the inequality

```math
\log_2(p+x+2)\le\log_2(p+1)+\log_2(x+2)+C
```

and
$`\log_2(p+1)\le Cp`$. The fixed shifts keep every logarithm
positive. An early stop has $`m\le Bp^3+1`$, so increasing the
absolute $`C_4`$ once gives, for the final height in every case,

```math
m\le C_4(p+z)^3,\qquad z=\log_2\log_2(\ell+4).        \qquad(8)
```

## 3. Price the complete terminal before choosing a dimension branch

Use the unchanged rounded program of
[R90 Section 8](DIRTY_CATALOGUE_TERMINAL.md#8-a-finite-full-frame-program-with-linear-in-dimension-length).
Its q bits and preserved-address table satisfy

```math
q\le(L+12)2^m,\qquad Q=2^{n-m},\qquad Qq\le3NL.        \qquad(9)
```

For a fixed catalogue factor at error an adequately small absolute
multiple of eta, the unabridged native synthesis bound in R90 (6) is

```math
S\le C2^{3m}(L+m+1).                                 \qquad(10)
```

Here S bounds its T-count, elementary Clifford count and T-depth with
an absolute constant. No fixed-L shorthand is substituted for (10).
The factors use no ancillary work beyond their m target wires.

Put $`U=q+S+mq+m+2`$. Equations (8)–(10), including their polynomial
prefactors, imply

```math
U\le2^{C_5(p^3+z^3)}.
```

For this fixed $`C_5`$ there is an absolute $`C_6`$ with

```math
C_5[\log_2\log_2(\ell+4)]^3
\le\tfrac14\log_2(\ell+4)+C_6
```

for all $`\ell\ge\log_2 3`$. This follows from
$`(\log_2 y)^3/y\to0`$; the remaining bounded interval has an
absolute maximum independent of L. Absorbing $`C_6`$ into a coefficient
of $`p^3`$ and applying $`xy\le(x^2+y^2)/2`$ gives

```math
U\le2^{C_7p^3}(\ell+4)^{1/4}
 \le2^{C_8p^3}+\sqrt{\ell+4}.                          \qquad(11)
```

In particular the additive q charge in the word query and XOR router,
and the local synthesis depth S, are explicitly bounded.

Pay for every catalogue slot and both disjoint arrays. Up to absolute
coefficients their total gate and workspace bill outside the query is

```math
M_{\rm cat}=2^q(S+mq)+q+2m2^q.
```

Since $`M_{\rm cat}\le C2^qU`$, (11) implies

```math
\log_2(M_{\rm cat}+2)
\le C[2^{C_8p^3}+\sqrt{\ell+4}].                       \qquad(12)
```

The factor $`2^q`$ is retained in both gate counts and array width,
including entries never selected by a supplied prefix.

## 4. An absolute large-dimension allocation

Fix an absolute integer $`A_1\gt4C_8`$. Enlarge it below to cover
finitely many absolute cutoffs, after all native constants have been
fixed. The catalogue construction will be used only on the branch

```math
n\ge2^{A_1p^3},\qquad b\ge\sqrt N/n.                  \qquad(13)
```

Then $`2^{C_8p^3}\le n^{1/4}`$ and
$`\sqrt{\ell+4}=O(n^{1/4})`$. Equation (12) is $`O(n^{1/4})`$,
so beyond an absolute dimension cutoff,

```math
M_{\rm cat}\le2^{n/8}\le\sqrt N/(4n).                 \qquad(14)
```

Enlarging $`A_1`$ makes (13), with $`p\ge7`$, imply this cutoff.
The program and both arrays occupy at most $`b/4`$. Reserve
$`B_Q=\lfloor b/2\rfloor`$ additional arbitrary dirty wires for the
program query, disjoint from every live program, array, address and
logical-target wire. At the same type of absolute cutoff,

```math
B_Q\ge16c_I(n-m+2),\qquad B_Q=\Theta(b),               \qquad(15)
```

so the retained rectangular query's actual minimum pool fits. Its
complete operator, including all actual inverses, has

```math
\begin{aligned}
T_Q&=O(q+\sqrt{Qq}+Qq/B_Q),\\
G_Q&=O(Qq),\\
D_Q&=O(q+\log_2(n-m+2)+Qq/B_Q^2).
\end{aligned}
```

Substituting (9), (11), and (14)–(15) into the complete R90 local
catalogue ledger proves

```math
\begin{aligned}
T_{\rm term}&=O(\sqrt{NL}+NL/b+\sqrt N),\\
G_{\rm term}&=O(NL),\\
D_{\rm term}&=O(NL/b^2+\ell+2^{C_8p^3}).
\end{aligned}                                         \qquad(16)
```

This pays the program-output term rather than inferring query depth from
the number of query invocations. The fixed local circuits and routers
require no additional initialized bank.

## 5. Coarse-stage pools and standalone terms on the same branch

On (13), $`L=O(\ell^{1/3})`$. Every executed height is at most n,
and (3) therefore gives, uniformly over its at most three stages,

```math
l_i,u_i=O(\ell),\qquad K_i=O(n).
```

Including the polynomial row factors in R89, each query address has
length at most

```math
n-r_i+C_0u_i^3+d_0\log_2K_i+C\le Cn.
```

Here $`\ell^3=o(n)`$ is used only beyond an absolute cutoff. There
are at most $`CK_i^{d_0}`$ simultaneous one-bit instances. Give each
a private pool $`\lfloor b/(C'K_i^{d_0})\rfloor`$, taking $`C'`$
large enough that their sum occupies a fixed fraction of b. At the
smallest wide allocation, each pool is at least

```math
c\,2^{n/2}/n^{d_0+1}-1,
```

which exceeds the sharp-query minimum $`C_Q(Cn+2)`$ after one absolute
cutoff. The other local dirty fields occupy $`O(n^e)=o(\sqrt N/n)`$
and fit in another fixed fraction. Their conditionally clean logical
fields have already been reserved inside $`r_i`$. No live local or
logical wire is also counted as an arbitrary dirty query helper.

The uniform one-stage ledger, with
$`\overline Q_i=2^{n-r_i+C_0u_i^3}`$, is

```math
\begin{aligned}
T_i&\le CK_i^d(\sqrt{\overline Q_i}+\overline Q_i/b+1),\\
G_i&\le CK_i^d(\overline Q_i+1),\\
D_i&\le C[\log_2(n+C_0u_i^3+2)+u_i+
                         K_i^d\overline Q_i/b^2].
\end{aligned}
```

The protected margin (4) gives

```math
K_i^d2^{-(r_i-C_0u_i^3)/2}\le K_i^{-40},\qquad
K_i^d2^{-(r_i-C_0u_i^3)}\le K_i^{-80}.
```

Summing at most three stages consequently yields

```math
\begin{aligned}
T_{\rm coarse}&=O(\sqrt N+N/b+n^d),\\
G_{\rm coarse}&=O(N+n^d),\\
D_{\rm coarse}&=O(N/b^2+\ell).
\end{aligned}
```

The standalone polynomial terms satisfy $`n^d=O(\sqrt N)`$ after
an absolute cutoff. Enlarge $`A_1`$ once to include all these cutoffs.
Each choice depends only on fixed native constants, never on L.

Completed query pools are reused only after return. Coarse stages and
the terminal occur in series, so their dirty peaks are maximized rather
than added. The two retained exact sign phases also have absolute costs

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\ell),
```

and their one-bit query minima fit this same wide branch. Together with
(16), this proves the claimed alternative resource envelope on (13).

## 6. Complete fallbacks at original accuracy

The alternative envelope must hold at all parameters before it can be
minimized with R79. Two complete existing compilers supply its remaining
branches without changing a literal width.

If $`b\lt\sqrt N/n`$, use the complete
[R72 real theorem](ONE_CLEAN_DEPTH.md#1-the-one-clean-theorem) at the
original tolerance eta. Its absolute additive depth is at most nL, while

```math
NL/b^2\gt n^2L\ge nL.
```

It therefore has $`D_T=O(NL/b^2)`$ with the required optimal count and
original real width. This invokes the full R72 theorem, including its
count-efficient circuit selection, rather than only its layerwise
fallback with the separate extra nL count term.

If $`n\lt2^{A_1p^3}`$, at any supplied width use the complete R79
real theorem at original eta. Since $`A_1`$ is now fixed,

```math
\ell\le A_1p^3+2,\qquad
\log_2(n+L+2)\le\ell+\log_2(L+2)=O(p^3).
```

Consequently

```math
F_{79}(n,L)\le Cp^6\le2^{C_9p^3}                     \qquad(17)
```

for an absolute $`C_9`$. The last inequality holds for all $`p\ge7`$
after fixing that one constant. No accuracy-dependent maximum over
finite dimensions is used.

The remaining branch is (13). Taking
$`C_R\ge\max\{C_8,C_9\}`$ proves, at all n, L and original real
widths, one circuit with the desired counts and

```math
D_T=O(NL/b^2+\ell+2^{C_R(L+1)^3}).                    \qquad(18)
```

Select this complete circuit or R79 according to the smaller displayed
additive bound. Every candidate already satisfies all count, depth,
width and error bounds simultaneously. This proves the real minimum
in (1), without combining resource estimates from different circuits.

## 7. Full-space error and returned work

On the new wide branch keep R90's whole-family magnitude/sign reduction
at error $`\eta/512`$. Each executed coarse stage has error at most

```math
2^{-l_i}/4\le\eta/[4096(m_i+2)].
```

There are at most three stages, so their sum is at most
$`\eta/2048`$. R90's depth-dependent program rounding costs less
than $`\eta/256`$, and its certified native commutator catalogue
costs less than $`\eta/32`$. The real error is therefore below

```math
\frac{\eta}{512}+\frac{\eta}{2048}+
\frac{\eta}{256}+\frac{\eta}{32}
=\frac{77\eta}{2048}\lt\frac{\eta}{16}.
```

The two disjoint arrays return exactly through the complete commutator
word, even when its fixed factors are approximate. Actual inverse words
cancel the unselected factors and both selected array residues. The
logical error is a maximum over the catalogue labels, with no multiplier
equal to the number of slots.

The catalogue acts on the entire input space and does not assume that
the earlier coarse stage returned an exactly initialized flag. Ordinary
initialized-isometry telescoping propagates all prior leakage unitarily.
This preserves every logical column, literal sign and phase, arbitrary
dirty input, and external reference. The complementary whole-theorem
fallbacks retain their original full error at most eta.

## 8. Complex phases and the retained high-precision clause

At $`b\ge17(L+n+8)`$, apply the real alternative theorem (18) at
$`k=L+1`$ and the
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at the same k. Their required width is exactly
$`17(k+n+7)=17(L+n+8)`$, as in R79 Section 6. Reuse the same sole
initialized qubit and returned dirty pool. Their full-isometry errors
add to at most $`2^{-L}`$.

The diagonal has absolute bounds

```math
T_D=O(\sqrt{Nk}+k+Nk/b),\qquad G_D=O(Nk),
```

```math
D_{T,D}=O(Nk/b^2+k+\ell).
```

Since $`k=\Theta(L)`$, the counts are unchanged. For $`L\ge6`$,

```math
(L+2)^3\le(8/7)^3(L+1)^3\lt2(L+1)^3.
```

Choose the final absolute $`C_*\ge2C_R`$, increasing it if needed
to absorb k. The real precision shift and diagonal overhead then fit
the alternative complex term in (1). Select it or the complete R79
complex circuit to obtain the minimum at the original complex width.

For either family, when $`L\gt N/n^2`$, select the stronger retained
R79/R72 high-precision circuit. It also satisfies (1), since both
$`F_{79}(n,L)`$ and $`\ell+2^{C_*(L+1)^3}`$ dominate L. This proves
(2) without weakening any simultaneous resource guarantee.

## 9. A growing-precision logarithmic-depth corollary

Fix the final $`C_*`$ after complex completion. Choose an absolute
$`c\gt0`$ small enough that $`C_*(7c/6)^3\le1/2`$. Whenever

```math
6\le L\le c[\log_2(\ell+2)]^{1/3},                    \qquad(19)
```

the inequality $`L+1\le7L/6`$ gives

```math
2^{C_*(L+1)^3}\le(\ell+2)^{1/2}=O(\ell).
```

Thus both prescribed families have, with absolute constants,

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+\ell),                                   \qquad(20)
```

at their original widths and full operator accuracy. The range is
unbounded as n grows. No L-dependent dimension cutoff remains; for
small n the permitted integer interval in (19) can simply be empty.

This corollary gives error tending to zero, but only very slowly. At
$`L=\Theta(\ell)`$, the explicit exponential precision cost in (1)
is larger than the R79 alternative, so this construction alone retains
its $`O(\ell^2)`$ additive bound. The precision exponent is an upper
bound for this particular terminal construction, not a necessary cost
for arbitrary circuits. Optimal general depth remains open.

The quantum words are unchanged from R89 and R90. Their bounded operator
tests retain their stated roles; the new result rests on the uniform
stopping, catalogue, width and fallback inequalities proved here. No
scalable circuit emitter or practical crossover is asserted.
