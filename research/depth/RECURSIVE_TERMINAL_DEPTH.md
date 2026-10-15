# Recursive compilation of prefix-selected terminal frames

[Subtree decoder](SUBTREE_PROJECTOR_DECODER.md) · [Uniform frame theorem](UNIFORM_SUBTREE_DEPTH.md) · [One-clean primitives](ONE_CLEAN_DEPTH.md) · [Common-flag family](HOMOGENEOUS_MAGNITUDE_DEPTH.md)

**Local result R81.** Any number of final layers of an arbitrary real
Hopf frame admits a recursive subtree-decoder construction. A preserved
external prefix selects the local frame. Its table rows retain that
prefix, while its initialized work and other polynomial costs depend only
on the local frame height. The sharp dirty query removes the global
address polynomial that would otherwise prevent recursion.

This improves a terminal component. The generic complete-frame frontier
remains [R79](UNIFORM_SUBTREE_DEPTH.md): its first decoder still has
sequential spectral comparisons on the full early frame. The
[level-homogeneous construction](HOMOGENEOUS_MAGNITUDE_DEPTH.md) has a
different common-flag mechanism that also removes those comparisons;
that mechanism is not assumed here.

## 1. Fixed-accuracy terminal theorem

Fix an integer $`L\ge6`$. For every integer $`n\ge1`$, put

```math
N=2^n,\qquad \eta=2^{-L},\qquad \ell=\log_2(n+2).
```

Let $`W_{\mathrm{tail}}`$ be the prescribed product of the deepest M
layers of a complete real Hopf frame, where M is an integer and

```math
0\le M\le n.
```

For every integer $`b\ge17(L+n+7)`$, one coherent Clifford+T circuit
V uses one external initialized qubit and at most b arbitrary dirty
qubits, and satisfies

```math
\|VJ_1-J_1(W_{\mathrm{tail}}\otimes I_b)\|_{\mathrm{op}}
\lt\eta/4,
```

```math
T=O_L\!\left(\sqrt N+\frac Nb\right),\qquad
G=O_L(N),
```

```math
D_T=O_L\!\left(
\frac N{b^2}+\ell\log_2(M+2)\right).
```

All bounds hold on the same circuit. The norm includes every logical
column, literal phases, source leakage and work return, arbitrary dirty
inputs, and external references. Unused supplied wires remain untouched.
The constants may depend on the fixed L, but not on n, M, b, or the angle
tuple. This is not an absolute all-precision theorem. For $`M=0`$,
the target and circuit are identity.

The angle interface is the existing terminating certified evaluation
interface, including signed, unequal, and singular angles. Classical
preprocessing has no efficiency claim. Every query, preparation, actual
inverse, and guard is charged. No QRAM, measurement, reset, supplied
precision state, or free coherent evaluator is used.

**Logarithmic-height corollary.** For every fixed $`C_0>0`$, the
previous range $`M\le\min\{n,C_0\ell\}`$ retains

```math
D_T=O_{L,C_0}\!\left(N/b^2+\ell\log_2(\ell+2)\right)
```

with the same count, workspace, and error guarantees.

**Subpolynomial protected-bank corollary.** Let $`M=M(n)\le n`$
satisfy $`\log_2(M+2)=o(\ell)`$. At sufficient square-root dirty
width, the terminal circuit then has $`D_T=o(\ell^2)`$ while
retaining $`T=O_L(\sqrt N)`$ and $`G=O_L(N)`$. Examples are

| Protected logical suffix M | Terminal T-depth at square-root dirty width |
|---|---|
| $`\lceil C\ell^c\rceil`$, fixed $`C>0,c\ge1`$ | $`O_{L,C,c}(\ell\log_2(\ell+2))`$ |
| $`\lceil2^{\sqrt\ell}\rceil`$ | $`O_L(\ell^{3/2})`$ |
| $`\lceil n^{1/\log_2\ell}\rceil`$ | $`O_L(\ell^2/\log_2\ell)`$ |

Cap each displayed M by n in finite dimensions. This prices the tail
left by an early construction that uses M logical suffix bits as
conditional initialized work. That early construction must still prove
its own native action, inactive identity, resource bounds, and work return;
the corollary does not supply an early decoder.

## 2. A decoder with a preserved external prefix

Suppose m bottom layers remain. Their unchanged upper
$`a=n-m`$ logical bits form an external prefix P. In each P sector the
remaining operator is a complete m-bit Hopf frame, with independently
prescribed angles. Choose a protected bottom bank Z of r bits and put
$`h=m-r`$. The next h layers have the exact target

```math
\sum_p |p\rangle\langle p|_P\otimes
\left(W_{p,h}\otimes\Pi_Z
      +I_X\otimes(I_Z-\Pi_Z)\right),
\qquad \Pi_Z=|0^r\rangle\langle0^r|.
```

Save $`H=[Z=0^r]`$ using the sole external clean qubit. Apply the
native [R78 decoder](SUBTREE_PROJECTOR_DECODER.md#2-the-literal-nested-subtree-spaces)
within each P sector, then reverse the actual saved-predicate word.
The external prefix is preserved by every operation. Descendant phase
and calibration queries include P in their addresses; local masks,
tag packing, threshold arithmetic, and saved predicates need no P input.

The full-mode reflection proof is unchanged. Its transported anchor
tests the entire fresh source bank, including the gradient, and the
appropriate local data suffix. It never tests P for zero. Source errors
therefore take a maximum over the orthogonal P sectors, rather than a
sum over their number. The reflection comparator is identity outside
the transported initialized-work subspace, as in R78.

On H zero, each completed source, query, and search stage is exactly
identity on arbitrary occupied work. On H one, the ideal stage restores
Z. Its final actual predicate inverse propagates all source leakage
unitarily, including possible final H leakage. This full-isometry
contract permits the same H to be reused at the next stage without
a reset.

At local precision l, set

```math
s=\lceil\log_2(h+1)\rceil,\qquad
\delta=\frac{2^{-l}}{64(s+1)},\qquad
p=1+\lceil\log_2(1/\delta)\rceil.
```

R78 gives stage error at most $`2^{-l}/8`$. Its three s-bit records
and fresh source bank require

```math
3s+u+O(1)=O(l+\log_2(m+2)),\qquad
2^u\le C\delta^{-2}p.
```

All these initialized bits lie inside Z. The external prefix changes
the table sizes, but not this initialized-work requirement.

## 3. Native costs depend polynomially on the local height

Write

```math
\overline Q=2^{a+h}=N2^{-r},\qquad K=m+l+2.
```

Use the [sharp rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence)
for the large tables. For R rows, q output bits, and a disjoint extra
dirty pool B, its bounds are

```math
T_Q=O\!\left(q+\sqrt{Rq}+\frac{Rq}{B}\right),\qquad
G_Q=O(Rq),
```

```math
D_{T,Q}=O\!\left(q+\frac{Rq}{B^2}
                  +\log_2(\lceil\log_2R\rceil+2)\right).
```

Its sufficient width condition is checked in Section 6. The query is
an exact arbitrary-output XOR, including its dirty return and literal
inactive rows. The deliberately conservative address-cubic bound in
R78 is not used for these global-prefix tables.

The [source ledger](COHERENT_DESCENDANT_LOADER.md#6-native-resource-interface)
and the complete R78 query multiplicities give the following factors.
Here C is absolute and may increase between bounds.

| Operation | Rows and outputs | Consequence |
|---|---|---|
| Descendant phase or calibration | At most $`C(h+1)\overline Qp`$ rows and $`Cp`$ outputs | Volume at most $`CK^3\overline Q`$ |
| All-work anchor phase | At most $`C(h+1)\overline Qp\delta^{-2}`$ rows, one output | Volume at most $`CK^4 4^l\overline Q`$ |
| Initial tag packing | $`2^h`$ rows and s outputs | Volume at most $`CK\overline Q`$ |
| Local masks and tag unpacking | $`O(K)`$ rows and $`O(K)`$ outputs | $`O(K^2)`$ Clifford count per call |

The anchor does not need P in its address; including it in the displayed
upper bound is harmless. The inequalities use
$`h+1,s\le K`$, $`p\le CK`$, and
$`\delta^{-1}\le CK2^l`$. There are $`O(s)\le CK`$ large-table
appearances, counting all actual inverses and the fixed expansion inside
each source. Thus their summed volumes are at most
$`CK^5 4^l\overline Q`$, and their summed square-root terms are at
most $`CK^3 2^l\sqrt{\overline Q}`$.

The wide local masks use the [parallel-output query](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address).
They have polynomial local dirty width and $`O(\log_2(K+2))`$
depth. Their h output bits do not create a serial h-depth term.
Large calibration tables have only $`O(p)`$ outputs, so their
rectangular output-depth term is already part of the precision cost.
Local gradient and coefficient preparations, arithmetic, native guards,
tag XOR, and saved-predicate boundaries have a safe total polynomial
envelope $`CK^4`$. All query pools return before reuse.

With $`B=\Theta(b)`$, degree ten conservatively covers these explicit
factors. One complete stage has

```math
\begin{aligned}
T_i\le C\bigg[&K^{10}2^l\sqrt{N2^{-r}}
 +\frac{K^{10}4^lN2^{-r}}b+K^{10}4^l\bigg],\\
G_i\le C\big[&K^{10}4^lN2^{-r}+K^{10}4^l\big].
\end{aligned}
```

Retaining the sharper logarithmic bounds on s and p gives

```math
D_{T,i}\le C\left[
\log_2(m+2)\big(\ell+l+\log_2(m+2)\big)
+\frac{K^{10}4^lN2^{-r}}{b^2}\right].
```

For the fixed-L construction below, every address has at most
$`n+O_L(\log(m+2))`$ bits, so its query-depth logarithm is
$`O_L(\ell)`$. No polynomial in global n has entered the count
or local-work envelope. The standalone $`4^l`$ term is deliberately
loose and will also be charged.

## 4. Recursive banks and summable errors

Choose one absolute $`A\ge40`$, also large enough for the
initialized-work constant in Section 2. Starting at $`m_0=M`$, set
at each nonfinal stage

```math
l_i=L+\lceil\log_2(m_i+2)\rceil+4,\qquad
K_i=m_i+l_i+2,
```

```math
r_i=\left\lceil A\big[l_i+\log_2K_i\big]\right\rceil,
\qquad m_{i+1}=r_i.
```

Choose a fixed stopping height $`m_*=m_*(L,A)`$ so that for every
$`m>m_*`$ the resulting remainder obeys

```math
r(m)\le m/2,\qquad
\log_2(r(m)+2)\le\tfrac12\log_2(m+2).
```

Such a height exists because $`r(m)=O_L(\log m)`$; choose it once,
independently of n and M. Stop when
$`m_k\le m_*`$. At every preceding stage
$`h_i=m_i-r_i\ge m_i/2`$, and the protected bank fits all three
records and the fresh source work.
If $`M\le m_*`$, omit recursion and use the final bottom directly.

The stage errors satisfy

```math
\epsilon_i\le\frac{2^{-l_i}}8
\le\frac{\eta}{128(m_i+2)}.
```

Reading the decreasing heights backward gives a geometric sequence,
so $`\sum_i1/(m_i+2)\le2`$. The total recursive error is at most
$`\eta/64`$. For the final $`m_k`$ original layers, use the
[smaller-error one-clean bottom](ONE_CLEAN_DEPTH.md#the-final-logical-layers-at-smaller-error)
with $`q_j=L+j+7`$. Its full error is strictly below
$`15\eta/64`$. Full-isometry telescoping therefore gives total
error strictly below $`\eta/4`$. The number of prefixes and
recursion stages adds no further error factor.

## 5. Resource sums

The protected-bank choice gives, at every nonfinal stage,

```math
K_i^{10}2^{l_i}2^{-r_i/2}\le K_i^{-10},\qquad
K_i^{10}4^{l_i}2^{-r_i}\le K_i^{-30}.
```

The sums of these inverse powers are bounded, because the heights
increase at least geometrically when read backward. Hence the large
tables contribute

```math
T=O\!\left(\sqrt N+\frac Nb\right),\qquad
G=O(N),\qquad D_{T,\mathrm{tables}}=O(N/b^2).
```

The standalone terms are also finite-budget costs:
$`4^{l_i}\le C_L(m_i+2)^2`$ and
$`K_i=O_L(m_i+2)`$, so

```math
\sum_iK_i^{10}4^{l_i}\le C_L(M+2)^{12}
\le C_L(n+2)^{12}=O_L(\sqrt N).
```

The polynomial-versus-exponential comparison is uniform over every
$`M\le n`$. The same terms fit $`O_L(N)`$ Clifford gates. There
is no unpaid iterated-logarithm count factor.

The stronger logarithmic contraction gives

```math
\sum_i\log_2(m_i+2)=O_L(\log_2(M+2)),
```

```math
\sum_i\log_2^2(m_i+2)=O_L(\log_2^2(M+2)).
```

Since $`l_i\le L+\log_2(m_i+2)+5`$, the complete recursive depth
is bounded by

```math
O_L\!\left(N/b^2+\ell\log_2(M+2)+\log_2^2(M+2)\right).
```

The square-logarithmic term is included in this sum. It is absorbed by
$`\ell\log_2(M+2)`$ because $`M\le n`$. The final constant-height
bottom has

```math
\begin{aligned}
T_{\mathrm{bottom}}&=O_L(\sqrt N+N/b+n^2),\\
G_{\mathrm{bottom}}&=O_L(N),\\
D_{T,\mathrm{bottom}}&=O_L(N/b^2+\ell).
\end{aligned}
```

Use $`n^2=O(\sqrt N)`$ and $`\log_2(M+2)\ge1`$ to obtain
the theorem's simultaneous bounds, uniformly in M.

## 6. All dirty widths at the literal threshold

For sufficiently large n, depending only on the fixed L and native
constants, split at

```math
b_c=\sqrt N/n.
```

When $`b\ge b_c`$, reserve $`B=b/4`$ for each sequential large
query. A separate $`b/4`$ allocation suffices for every local
polynomial-in-M mask and helper pool. Uniformly over the recursion,
$`m_i\le M\le n`$, $`l_i=O_L(\log(n+2))`$, and
$`K_i=O_L(n+2)`$. Thus one fixed polynomial in n bounds all local
dirty work. The query addresses have
$`n+O_L(\log(M+2))=O_L(n)`$ bits. Exponential dominance gives
one dimension cutoff, independent of M, at which the local pool fits
and the sharp query condition
$`B\ge16c_I(\lceil\log_2R\rceil+2)`$ holds.
The final bottom also satisfies $`b/2\ge16c_I(n+2)`$, as required
for its shallow query schedule. The pools are arbitrary dirty work,
disjoint from live addresses, targets, and protected logical work.

When $`b\lt b_c`$, retain the original ungrouped one-clean circuit for
just the deepest M layers. Specifically, use R72's smaller-error
allocation

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad B'_0=L+n+9.
```

At the original supplied threshold its remaining query pool satisfies

```math
B=b-B'_0\ge16b/17-2\ge b/2.
```

Thus the routed query fits. This is the proved smaller-error
reservation; it is not an invocation of the generic theorem with L
replaced by $`L+2`$ and a larger dirty-width threshold. Its error is
below $`\eta/4`$. The weighted query sums over $`1\le j\le M`$
are bounded by their convergent infinite sums:

```math
\sum_{j=1}^M\sqrt{N2^{-j}(L+j+8)}=O_L(\sqrt N),\qquad
\sum_{j=1}^MN2^{-j}(L+j+8)=O_L(N).
```

The standalone source and predicate terms are $`O_L(nM+M^2)`$,
at most $`O_L(n^2)=O_L(\sqrt N)`$. Thus the same routed word
has $`T=O_L(\sqrt N+N/b)`$ and $`G=O_L(N)`$ uniformly in M.
Its depth is at most

```math
O_L(N/b^2+nM+M^2)=O_L(N/b^2),
```

because $`M\le n`$ and $`b\lt b_c`$ imply $`N/b^2>n^2`$.

For the finitely many n below the dimension cutoff, use the same routed
fallback for every b. The cutoff depends only on L and native constants,
so there are finitely many pairs $`(n,M)`$ with $`0\le M\le n`$
below it. The routed count inequalities still hold for every b; the
bounded remaining depth is absorbed by
$`O_L(\ell\log_2(M+2))`$, uniformly over those pairs and the
angle tuple, even as b grows arbitrarily large. Extra supplied wires may
remain unused. This completes the literal all-width statement.

## 7. Scope of the improvement

At $`M=\Theta(\log n)`$, the first recursive stage still has
$`\Theta(\log\log n)`$ sequential R78 comparisons. Its source
tables retain the global prefix, and their proved query T-depth is
$`O(\log n)`$. Their product gives the displayed terminal bound.
Sharper error sums do not remove this dependency.

This is a limit of the present construction, not an unrestricted
depth lower bound. Starting the same recursion at $`m=n`$ leaves
the generic $`O(\log^2 n)`$ complete-frame bound. At variable
precision, the original gradient and program also require
$`\Theta(L)`$ initialized bits, so this argument does not establish
an absolute $`O(L+\log^2 n)`$ theorem.

The proof is an analytic extension and resource allocation of retained
native primitives. Small-dimensional circuit checks do not establish
its asymptotic recursion or substitute for the operator and width
arguments above.
