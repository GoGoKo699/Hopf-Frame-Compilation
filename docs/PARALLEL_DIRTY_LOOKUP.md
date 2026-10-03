# Parallel dirty lookup while retaining the T-count bound

[T-depth model and routing](T_DEPTH_COMPILER.md) · [Grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md) · [Research status](OPEN_PROBLEM.md)

Additional dirty workspace can parallelize the selector computation without
increasing the count-efficient number of word banks. This gives a simultaneous
count and depth guarantee for the prescribed complete real Hopf frame.

**Routed baseline.** Write $`N=2^n`$, $`n\ge1`$,
$`0\lt\eta\le1/64`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
$`B_0=L+n+7`$, and $`\ell_*(n)=1+\log_2^*(n+2)`$.
There is a fixed sufficient constant C such that, for

```math
a\ge2,\qquad b\ge C\bigl(B_0+\sqrt{NL}\bigr),
```

one coherent Clifford+T circuit satisfies all of

```math
\|VJ_a-J_a(W\otimes I_b)\|\le\eta,
\qquad T=O\!\left(\sqrt{NL}+L\ell_*(n)\right),
\qquad G=O(NL),
```

```math
D_T=O\!\left(\min\{nL+n^2,\ L\ell_*(n)+n^3\}\right).
```

Only two external clean qubits are used. The complete-input error includes
clean leakage and arbitrary dirty inputs with references. The coefficient
tables remain classically specified and their quantum lookups are fully
charged. The sufficient constant C is not a practical crossover estimate.

**Hybrid theorem.** Define

```math
\chi(t)=\log_2(t+2).
```

For every $`b\ge17B_0`$, Sections 5–6 and the
[masked-sum refinement](DIRTY_SUM_COMPRESSION.md) give a two-clean complete
real-frame circuit with the same error contract and

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL+n\chi(n)\right).
```

Both count and depth match their worst-case lower bounds when
$`17B_0\le b\le\sqrt{NL/(nL+n\chi(n))}`$. At fixed accuracy,
sufficiently large $`b=\Theta(\sqrt N)`$ permits **one circuit** to have

```math
T=O(\sqrt N),\qquad D_T=O\!\left(\min\{n^2,n\chi(n)\}\right),
\qquad G=O(N).
```

The earlier routed circuit supplies the alternative quadratic bound.
The existing worst-case lower bound makes this T-count optimal in order.
The T-depth upper bound is not proved optimal, and does not bound total
elementary depth. The selected allocation $`L=N,b=N+n+7`$ is not covered
by the new sufficient-width hypothesis; its linear T-count endpoint remains
open.

At fixed accuracy, the [blocked bilinear refinement](BLOCKED_BILINEAR_LOOKUP.md)
uses conditional unary groups and width-constrained bilinear blocks to
give $`D_T=O(N/b^2+n)`$, retaining optimal-order T-count throughout
$`b\ge17B_0`$. Count and depth match through $`b\le\sqrt{N/n}`$
when nonempty. Its fixed-accuracy proof does not replace the uniform
precision statements above.

## 1. An exact dirty indicator by conjugated routing

All sequences in this chapter are chronological. For a k-bit address x,
put $`H=2^k`$ and let Y consist of H arbitrary dirty bits. Regard these
bits as one-bit banks and use the exact
[bank router](T_DEPTH_COMPILER.md#a-bank-router) $`\mathcal R_x`$.
It moves the original bit $`Y_x`$ to position zero, preserving x. Apply

```math
\mathcal R_x,\quad X_{Y_0},\quad\mathcal R_x^\dagger.
```

The resulting indicator $`I_k`$ has the literal action

```math
|x,Y\rangle\longmapsto|x,Y\oplus e_x\rangle.
```

To prove this, fix a basis address x. The router is a permutation of the
Y wires. Flipping its position zero and undoing that same permutation
flips exactly the original position x; every other bit is restored.
The actual inverse is essential: the router need not be an involution.
Linearity and the exact Fredkin phases extend the identity to arbitrary
superpositions of addresses, dirty outputs, and reference correlations.
No bit of Y is initialized and no additional scratch is used.

At level j the router uses address bit j to swap representatives of
adjacent blocks of size $`2^j`$. Thus it selects
$`x=\sum_{j=0}^{k-1}2^jx_j`$, consistent with little-endian table labels.
There are $`H-1`$ Fredkins in k shared-control batches. The literal
four-T-layer batch decomposition gives

```math
D_T(I_k)\le8k,\qquad T(I_k)\le14(H-1),\qquad G(I_k)=O(H).
```

For $`k=0`$ the word is just X. A serial elementary schedule costs
$`O(H)`$ depth. In particular, the T-depth bound does not assume
constant-depth Clifford fanout or routing. This direct conjugation
replaces the earlier recursive bilinear indicator, whose quadratic
address-depth and extra scratch were unnecessary for this T-depth task.

### 1.1. A full-input lower bound

An exact nonaffine classical permutation needs at least two T layers,
even with arbitrarily many returned dirty helpers and arbitrary Clifford
circuits between T layers. This is a full-Hilbert-space statement: no
helper is promised to start in a particular state.

Here is a short specialization of the one-layer Pauli-conjugation
argument in [Selinger, Proposition 5.1](https://arxiv.org/html/1210.0974v2).
For a permutation matrix U on w qubits, every Pauli-transfer entry

```math
2^{-w}\mathrm{Tr}(QUPU^\dagger)
```

is a dyadic rational. Indeed, the matrix entries in the trace sum are
Gaussian integers, and the trace is real for Hermitian Paulis P and Q.
Tensoring U with an identity on arbitrary dirty helpers preserves this
property. Suppose instead that an implementation has one nonempty T
layer, written $`V=C_2\tau C_1`$, where the C operators are arbitrary
Cliffords and tau is a tensor product of T, its inverse, and identities.
Choose an active wire j and $`P=C_1^\dagger X_jC_1`$. Then

```math
VPV^\dagger=(Q\pm R)/\sqrt2,
```

where Q and R are distinct Hermitian Paulis. Its two nonzero
Pauli-transfer coefficients are irrational, a contradiction. An empty
T layer leaves a Clifford; a Clifford that is a classical permutation
acts affinely on bit strings. Thus a nonaffine permutation cannot have
T-depth zero or one. Global phases do not affect this argument.

For $`k\ge2`$, a coordinate of $`I_k`$ contains the degree-k Boolean
monomial in the address bits, so its classical action is nonaffine.
Consequently $`D_T(I_k)\ge2`$, independently of returned dirty width.
This constant bound is not a growing lower bound and does not establish
optimality of the routed construction. It also does not apply to a
circuit required to work only on an initialized clean subspace; the
physical unitary outside that subspace need not be a permutation.

### 1.2. The two-bit indicator has optimal T-depth two

First implement a Toffoli with one arbitrary dirty helper. For four bits
a, b, c, d, the integer parity identity

```math
\sum_{s\in\{0,1\}^3}(-1)^{|s|}
 \bigl(d\oplus(s_1a\oplus s_2b\oplus s_3c)\bigr)
 =-4abc(-1)^d\equiv4abc\pmod8
```

implements the literal CCZ phase, independently of d. Split its eight
parities into two ordered bases over the binary field:

```math
\mathcal A=(d,\ d\oplus a,\ d\oplus b,\ d\oplus c),
```

```math
\mathcal B=(d\oplus a\oplus b,\ d\oplus a\oplus c,
 d\oplus b\oplus c,\ d\oplus a\oplus b\oplus c).
```

Both bases are invertible: differences with d in the first, or with the
last entry in the second, recover a, b, c and then d. For each basis,
compute it by CNOTs on the four physical wires, apply four simultaneous
T or inverse-T gates with the displayed parity signs, and undo the
actual CNOT word. The signs are $`(+,-,-,-)`$ for A and
$`(+,+,+,-)`$ for B. The result is exactly CCZ on a, b, c tensor the
identity on d, with eight T gates in two T layers. Hadamards on c turn
this into a Toffoli with the same resources. All-input equality also
returns a helper entangled with any external reference.

Now let a be the high address bit and b the low bit, and order the four
dirty outputs by $`00,01,10,11`$. The indicator increment is

```math
e_{2a+b}=(1\oplus a\oplus b,\ b,\ a,\ 0)
       \oplus ab(1,1,1,1).
```

The first vector uses only X and CNOT. For the second, let F fan out
$`Y_{11}`$ into the other three output wires. Chronologically apply
F, the above Toffoli with target $`Y_{11}`$, and F again. During the
Toffoli, borrow $`Y_{00}`$ as its arbitrary dirty helper. The Toffoli
returns that wire exactly before the final fanout, so the net change is
ab in all four outputs. This uses the existing six wires, no additional
helper, eight T gates, and T-depth two. Together with Section 1.1,

```math
D_T^\star(I_2)=2.
```

No optimal T-count is asserted. For comparison, $`I_0`$ and $`I_1`$
are Clifford and have T-depth zero. The
[native checks](../tests/test_dirty_indicator_depth.py) verify both
parity bases, literal phases, actual inverses, and the full six-wire
indicator matrix. This closes a bounded base case; Section 6 separately
supplies the general shallow indicator used in Section 5. The phase-polynomial technique
and one-layer Pauli normal form are established tools; no generic
Toffoli synthesis priority is claimed.

## 2. Indicator workspace and exact return

The H outputs are also the router's only target workspace. During the
indicator they are disjoint from the address, lookup word banks, and
query output. The intermediate permutation may disturb every indicator
bit; only the completed XOR action is used by the loader below.

The indicator is logically self-inverse, but use its actual reversed
native word when an inverse is required. In a completed loader, both
indicator calls finish before another register is reused. Therefore an
H-bit pool suffices, with no recursive scratch, clean copies, or assumed
zero sector. The same pool is reused between completed queries.

## 3. Count-efficient whole-word queries

Consider a Q-row, m-bit XOR table, padding Q to a power of two with zero
rows if needed. Let $`1\le\mu\le Q`$ be a power-of-two bank count and
put $`H=Q/\mu`$. The high address specifies the concatenation of mu
table words. Its length is $`\mu m`$.

Reserve the dirty bank word Z of length $`\mu m`$, an H-bit dirty
indicator Y, with no additional indicator scratch. Let A be the classical
binary matrix whose column for each high address is its concatenated
bank word. The exact linear map

```math
C_A:\ (Y,Z)\longmapsto(Y,Z\oplus AY)
```

uses at most $`H\mu m=Qm`$ CNOTs. The chronological echo

```math
C_A,\ I_k,\ C_A,\ I_k^\dagger
```

changes Z by $`Ae_x`$ and restores all indicator work. This is an
all-input high-address loader, including arbitrary initial banks.
It costs $`O(H)`$ T gates, $`O(Qm)`$ Clifford gates, and
$`O(1+\log H)`$ T-depth. A serial schedule of the CNOT matrix
has depth at most $`Qm`$; no stronger claim about that depth is needed.

Use this loader in the [completed SelectSwap query](T_DEPTH_COMPILER.md#2-exact-lookup-with-depth-optimized-banks): load, route, copy to the
arbitrary output, undo routing and loading, then route, copy, and undo
routing again. The original selected bank contents cancel. The exact
resources are bounded by

```math
T_{\rm query}=O(Q/\mu+\mu m),\qquad G_{\rm query}=O(Qm),
```

```math
D_{T,\rm query}=O\!\left(1+\log(Q/\mu)+\log\mu\right)=O(1+\log Q).
```

An elementary-depth upper schedule is
$`O(Qm+H+\mu m)`$: the linear maps can be serial, and
the routing uses only $`O(\mu m)`$ elementary gates. In particular,
low T-depth does not imply equally low Clifford or total depth.

| Live registers | Dirty width | Lifetime and return |
|---|---|---|
| Existing source, selectors, and helpers | $`B_0`$ | Reserved throughout; existing source return is charged in the frame error |
| Additional word banks | $`\mu m`$ | Live throughout the completed query; returned exactly |
| Indicator output | $`H`$ | Live during each loader; returned exactly |

Coefficient-query outputs already lie in the base reservation;
coarse-symbol queries reuse their existing private buffers. Neither
requires an additional initialized output. All new work fits when

```math
B_0+\mu m+Q/\mu\le b.
```

Zero completed table rows give the identity on inactive sectors, even
though the internal indicator and router may act there. There are no
new source calls, logical clean-bit assumptions, or supplied quantum
oracles in this replacement.

Choose mu as the largest power of two at most $`\sqrt{Q/m}`$ when
$`Q\ge m`$, and choose $`\mu=1`$ otherwise. For the former case,

```math
\mu m+Q/\mu\le3\sqrt{Qm};
```

for the latter it is at most $`2m`$. Thus a sufficient extra dirty pool
is $`O(m+\sqrt{Qm})`$. This retains the count-efficient bank scale;
the extra qubits hold returned selectors rather than more word banks.

## 4. Full-frame composition

Use the same coefficient tables, source precision, coarse programs,
amplification, and conditional suffix sectors as the existing grouped
compiler. Its weighted table estimates give

```math
\sum_gQ_gm_g=O(NL),\qquad
\sum_g\sqrt{Q_gm_g}=O(\sqrt{NL}),\qquad
\sum_gm_g=O(L\ell_*(n)+n).
```

Every source width is at most $`B_0`$. Every coefficient or coarse-symbol
query has $`Qm=O(NL)`$ and $`\log Q=O(n+1)`$, including the fixed
reserved tail. For example, a group ending above r suffix bits has
$`Q_g=O(s2^{n-r})`$ with $`s\le2^{r/C_1}`$ and $`C_1>2`$;
mode padding changes only fixed constants. A single pool of size
$`C(B_0+\sqrt{NL})`$ consequently fits every query, and is reused
between queries instead of summed over groups.

The preceding query lemma gives the established count-efficient T-count
for all tables. The coarse-symbol stream additionally uses the existing
estimate $`\sum_gs_g2^{e_g/2}=O(\sqrt N)`$, which includes its
repeated symbol queries. The scalar sources retain their total
$`O(L\ell_*(n)+n)`$ count and serial depth. The old
polynomial predicate and interpreter work fits the same T-count bound
because every fixed polynomial in n is $`O(\sqrt N)`$.
Hence the grouped schedule has

```math
T=O(\sqrt{NL}+L\ell_*(n)),\qquad G=O(NL).
```

Each coefficient query now has T-depth $`O(n+1)`$. There are at most
n groups and $`O(s_g^2)`$ coarse-symbol queries per group. Since
$`\sum_gs_g^2\le n^2`$, their total query depth is $`O(n^3)`$.
All other schedules are exactly those already charged in the
[depth proof](T_DEPTH_COMPILER.md#3-composition-with-the-grouped-full-frame-compiler).
This gives $`D_T=O(L\ell_*(n)+n^3)`$.

Applying the same query replacement to the older layerwise compiler gives

```math
T=O(\sqrt{NL}+nL),\qquad
D_T=O(nL+n^2),\qquad G=O(NL).
```

Select the schedule with the smaller displayed depth expression. If the
layerwise expression is smaller, then
$`nL\le L\ell_*(n)+n^3`$, and $`n^3=O(\sqrt N)`$ absorbs the
additional term in its T-count. Thus the selected circuit has the
simultaneous bounds in the theorem, not just two guarantees on different
circuits.

Every replaced query implements the same literal unitary as before.
The amplification, inverse, inactive-sector, and complete-input error
proofs therefore carry over unchanged. The QBP bias and dirty-reference
guarantees follow from the existing approximation theorem.

## 5. A bilinear query reduction

The following exact query removes the word-bank router by using two dirty
indicators. Its depth depends on the indicator implementation. The
general hybrid below isolates that dependence, and Section 6 supplies
the polylogarithmic-depth implementation used in the improved bound.

### Exact bilinear oracle

Split an r-bit address into a and b of lengths $`r_a+r_b=r`$.
Put $`H=2^{r_a}`$, $`J=2^{r_b}`$, and
$`Q=HJ`$. For an arbitrary m-bit table, let
$`D_\ell\in\mathbb F_2^{H\times J}`$ contain output bit ell. Introduce
dirty registers $`Y\in\mathbb F_2^H`$ and $`X\in\mathbb F_2^J`$,
disjoint from the address and arbitrary query output Z. Define

```math
\mathcal B:(Y,X,Z)\longmapsto
\left(Y,X,\left(Z_\ell\oplus Y^{\mathsf T}D_\ell X\right)_{\ell=1}^m\right).
```

This is an involution. For $`\rho_\ell=\mathrm{rank}(D_\ell)`$,
binary row and column elimination gives invertible $`P_\ell,Q_\ell`$ with
$`P_\ell D_\ell Q_\ell=J_{\rho_\ell}`$, where the rectangular matrix on the
right has rho diagonal ones. For this output bit, change coordinates by

```math
Y'=P_\ell^{-\mathsf T}Y,\qquad X'=Q_\ell^{-1}X.
```

Then $`Y^{\mathsf T}D_\ell X=Y'^{\mathsf T}J_{\rho_\ell}X'`$.
The middle operation consists of rho Toffolis with disjoint control pairs
$`(Y'_i,X'_i)`$ and common target $`Z_\ell`$. Conjugating that target
by a Hadamard makes them shared-control CCZ gates. The
[four-layer phase-polynomial schedule](T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
therefore implements them with T-depth at most four and T-count
$`6\rho_\ell+(\rho_\ell\bmod2)`$, preserving literal phase. Undo
both coordinate changes before proceeding to another output bit.

Elimination and its actual inverse use
$`O(\rho_\ell(H+J))=O(HJ)`$ elementary Clifford gates for each nonzero
matrix. No auxiliary qubit is required. Thus, writing
$`R=\sum_\ell\rho_\ell`$, a sequential schedule over the output bits has

```math
T(\mathcal B)\le6R+\sum_\ell(\rho_\ell\bmod2),\qquad
D_T(\mathcal B)\le4m,\qquad G(\mathcal B)=O(Qm).
```

In particular, $`T(\mathcal B)=O(m\min(H,J))`$. Every completed
output-bit operation restores Y and X, although its internal Clifford
changes need not preserve them individually.

### The two indicator echoes

Let $`I_a:Y\mapsto Y\oplus e_a`$ and
$`I_b:X\mapsto X\oplus e_b`$ be exact dirty indicators, preserving their
addresses and returning all their own helpers. Execute the chronological
sequence

```math
\mathcal B,\ I_a,\ \mathcal B^\dagger,\ I_b,\quad
\mathcal B,\ I_a^\dagger,\ \mathcal B^\dagger,\ I_b^\dagger.
```

The four bilinear calls add, for each output bit,

```math
Y^{\mathsf T}D_\ell X
\oplus(Y+e_a)^{\mathsf T}D_\ell X
\oplus(Y+e_a)^{\mathsf T}D_\ell(X+e_b)
\oplus Y^{\mathsf T}D_\ell(X+e_b)
=D_\ell[a,b].
```

All of X, Y, and the indicator helpers return exactly. Use the actual
reversed native word for every dagger. This computational-basis identity
retains literal phase and hence holds for arbitrary inputs and reference
correlations.

The query uses four bilinear oracles, two a-indicators, and two
b-indicators, counting actual inverses at the same cost. If $`w_a,w_b`$
are their additional dirty helper widths, the live dirty width, excluding
address and query output, is

```math
H+J+\max(w_a,w_b).
```

The common helper pool is reused only after a completed indicator. The
bilinear oracle requires none of it. For each resource
$`C\in\{T,G,D_T\}`$, the query has the explicit upper ledger
$`4C(\mathcal B)+2C(I_a)+2C(I_b)`$. In particular its T-depth is at
most $`16m+2D_T(I_a)+2D_T(I_b)`$. There is no remaining bank router.

The routed indicators already proved in Section 1 make this an
unconditional alternative query. For $`Q\ge m`$, choose H as the
largest power of two at most $`\sqrt{Q/m}`$ and put $`J=Q/H`$.
Then $`H\le J`$ and $`H+J+mH=O(\sqrt{Qm})`$, giving

```math
T=O(\sqrt{Qm}),\qquad G=O(Qm),\qquad
D_T=O(m+\log(Q+1)),\qquad w=O(\sqrt{Qm}).
```

For $`Q\lt m`$, taking $`H=1,J=Q`$ instead gives T-count
$`O(m+Q)=O(m)`$ and dirty width $`O(Q)`$. These bounds require no
shallow-indicator hypothesis. They recover count efficiency but do not
improve the existing full-frame depth order with the routed indicators.

### Hybrid with the layerwise frame compiler

Suppose that, for every k, an exact all-input dirty indicator on
$`S=2^k`$ outputs has T-count, elementary Clifford count, and additional
dirty width each at most $`O(SP(k+1))`$, for a fixed nondecreasing
polynomial $`P\ge1`$, and T-depth at most $`f(k)`$. The address must
return unchanged and every helper must return exactly, including on
reference-entangled inputs. The routed indicator of Section 1 satisfies
the count and width conditions with $`f(k)=O(k)`$. Section 6's counter
construction, with [masked compression](DIRTY_SUM_COMPRESSION.md), supplies
$`f(k)=O(\chi(k))`$, defined below. We first prove
the composition for a general f.

For $`r_a=\lfloor r/2\rfloor`$, $`r_b=\lceil r/2\rceil`$, the
exact query just proved then has

```math
T=O\!\left(\sqrt Q\,[m+P(r+1)]\right),\qquad
G=O\!\left(Qm+\sqrt Q\,P(r+1)\right),\qquad
w=O\!\left(\sqrt Q\,P(r+1)\right),
```

and depth $`O(m+f(r_a)+f(r_b))`$. The polynomial overhead is acceptable
on early frame layers even though it need not preserve the optimal
standalone lookup count.

### Every eligible width and precision

The hybrid does not require fixed accuracy or square-root-scale width.
For every $`L\ge6`$ and $`b\ge17B_0`$, put

```math
M_n=L+4+\lceil\log_2(8n)\rceil,\qquad
A_n=M_n+P(n+2),\qquad
w_b=\min\{b-B_0,\sqrt{NL}\}.
```

The [capped layer precisions](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision)
satisfy $`m_d\le M_n`$. Set $`k=n-d`$, so
$`Q_d=4N2^{-k}`$ and $`r_d=d+2\le n+1`$. Choose a sufficiently
large fixed constant kappa and the cutoff

```math
k_0=\min\!\left\{n,
 \left\lceil2\log_2\frac{\kappa\sqrt N\,A_n}{w_b}\right\rceil
 \right\}.
```

Use the bilinear query when $`k\gt k_0`$, and the existing amortized
query for the last $`k_0`$ layers. If the early set is nonempty,
geometric summation gives

```math
\sum_{k>k_0}\sqrt{Q_d}\,A_n
=O\!\left(\sqrt N\,A_n2^{-k_0/2}\right)=O(w_b).
```

Increasing kappa makes every early query's helper pool fit inside
$`b-B_0`$; the complete base reservation is retained. Thus the literal
threshold $`17B_0`$ is unchanged. An empty early set simply uses the
old circuit. Early query T-count is $`O(w_b)\le O(\sqrt{NL})`$.
Their indicator Clifford count has the same bound; bilinear Clifford
cost is bounded by $`\sum_dQ_dm_d=O(NL)`$. Late queries are a subset
of the old schedule and retain its count bounds. Sources, predicates,
and their complete error proof are unchanged. Therefore

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL).
```

Write $`F_n=\max_{0\le k\le n+1}f(k)`$. Early query depth is
$`O(n[M_n+F_n])`$. Late routing costs $`O(nk_0)`$, and late
chunk-depth terms sum to $`O(NL/b^2)`$. Sources and suffix predicates
cost $`O(nL+n\log(n+1))`$ depth. Since P is fixed polynomial and
$`b-B_0\ge16b/17`$,

```math
k_0=O\!\left(L+\log(n+2)
 +\log_+\frac{\sqrt{NL}}{b}\right),
\qquad \log_+x=\max\{0,\log_2x\}.
```

The apparent extra workspace logarithm is absorbed in the existing
terms. With $`x=NL/b^2`$, the inequality
$`\ln x\le\ln n+x/n`$ for positive x implies

```math
n\log_+x=O\!\left(x+n\log(n+2)\right).
```

Consequently the general indicator interface gives, on one two-clean
complete real-frame circuit at every eligible width and precision,

```math
D_T=O\!\left(\frac{NL}{b^2}+nL+n[F_n+\log(n+2)]\right).
```

Exact query replacement preserves the full-isometry error, dirty/reference
return, literal phases, and actual-inverse amplification. This argument
uses an adjustable cutoff, not additional initialized work.

### The improved matching range

The [compressed counter construction](DIRTY_SUM_COMPRESSION.md) supplies
$`f(k)=O(\chi(k))`$, where

```math
\chi(t)=\log_2(t+2).
```

It follows that, for all $`L\ge6`$ and $`b\ge17B_0`$,

```math
D_T=O\!\left(\frac{NL}{b^2}+nL+n\chi(n)\right)
```

with the preceding same-circuit T and Clifford counts. When nonempty,

```math
17B_0\le b\le\sqrt{\frac{NL}{nL+n\chi(n)}}
```

is a simultaneous matching range. Indeed the additive depth terms are
at most $`NL/b^2`$. Also $`b\le\sqrt{NL}`$ and $`b\ge1`$, so
$`\sqrt{NL}\le NL/b`$ and $`nL\le NL/b^2\le NL/b`$. The inherited worst-case
count and depth lower bounds therefore match:

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

At fixed L, the upper endpoint has order
$`\sqrt{N/(n\chi(n))}`$, asymptotically larger than the old
$`\sqrt N/n`$ interval; retain the old interval separately at small n.
At sufficient $`b=\Theta(\sqrt N)`$, the bound becomes
$`D_T=O(n\chi(n))`$ with optimal-order $`T=O(\sqrt N)`$ and
$`G=O(N)`$. More generally, $`b=\Theta(\sqrt{NL})`$ gives depth
$`O(nL+n\chi(n))`$ whenever the allocation is eligible; T-count is
optimal in order when $`L\le N/n^2`$. No high-precision endpoint
improvement or unrestricted matching large-width depth follows.

This is a specialization of the dirty-indicator and bilinear framework in
[Low–Kliuchnikov–Schaeffer, Appendix C](https://arxiv.org/html/1812.00954v2),
combined with the existing rank reduction and exact shared-control phase
schedule. The local result is an explicit resource and all-input interface
for the Hopf composition, not a new generic bilinear principle.
The [bounded bilinear checks](../tests/test_bilinear_dirty_lookup.py)
verify rectangular basis orientation, literal phases, actual inverses,
complete dirty return, and emitted resource counts using the existing
routed indicators. The counter construction has separate bounded checks.

## 6. A polylogarithmic-depth indicator using dirty counters

This section gives the binary-adder-tree baseline and the surrounding
counter echoes. [Pipelined masked compression](DIRTY_SUM_COMPRESSION.md)
replaces its sum operation, improving the indicator depth to $`O(\chi(k))`$
with polynomially larger per-row counts and width. The rest of this
section's exact circuit remains in use.

For $`k\ge1`$, put $`m=\lceil\log_2(k+1)\rceil`$ and
$`M=2^m\gt k`$. An exact indicator on $`H=2^k`$ arbitrary dirty
outputs can be implemented with

```math
T,G=O\!\left(H(k+1)\log(k+2)\log\log(k+4)\right),\qquad
w=O\!\left(H(k+1)\log(k+2)\right),
```

```math
D_T=O\!\left(\log(k+2)[\log\log(k+4)]^2\right).
```

Here w is additional dirty helper width; all helpers return exactly, on
all inputs. No clean qubit is used. The construction first computes a
read-only conjunction, then runs all equality rows in parallel with an
explicit shared-address schedule. It trades polynomially more dirty
work for depth; Section 5 absorbs that overhead on early frame layers.

### Exact modular adders

For addition $`\mathrm{ADD}_m(A;B):B\mapsto B+A\pmod M`$
preserving arbitrary A, use the exact ripple-carry circuit of
[Takahashi–Tani–Kunihiro, Sections 2.1–2.3](https://arxiv.org/pdf/0910.2530v1).
Their carry-output wire z is arbitrary and is targeted only by the
CNOT from $`A_{m-1}`$ in Step 2 and the final Toffoli of Step 3. Delete those two gates
and z; no remaining gate depends on z. The resulting two-register word
is exactly modular addition, on all inputs, with no helper or initialized
bit. For $`m\ge2`$ it uses $`2m-2`$ Toffolis and $`5m-6`$
CNOTs; for m equal to one it is CNOT. Thus $`T,G,D_T=O(m)`$.
The source's later clean-work and unbounded-fanout constructions are
not used. A may change during this adder and is restored at completion;
both registers are private in the application below. Every subtraction
uses the actual reversed native addition word.

A faster private adder follows from
[Remaud–Vandaele, Algorithm 3 and Lemmas 2/4](https://arxiv.org/html/2501.16802v2).
Remove the carry output at the abstract ladder level, before optimizing:
shorten Slice 2's CNOT ladder to $`(A_1,\ldots,A_{m-1})`$, and
Slice 3's inverse Toffoli ladder to

```math
(A_0,B_0,\ldots,A_{m-2},B_{m-2},A_{m-1}).
```

In the serial forms this removes precisely the two gates targeting the
carry output; all other slices remain. Apply the source's exact ladder
synthesis to these shortened lists. No carry wire is then present to be
borrowed internally. The result has no additional helper and

```math
T,G=O(m\log(m+2)),\qquad D_T=O(\log^2(m+2)).
```

Both operands are private, so the synthesis may temporarily borrow their
bits. The complete adder restores its first operand and adds it modulo
M into the second. The linear-size TTK adder supplies the retained
signed-increment baseline; the faster adder supplies the sum tree.

### A read-only controlled increment using two additions

Let U be an arbitrary m-bit counter and g an arbitrary m-bit helper.
Write $`g_0`$ for its low bit and $`e=(-1)^{g_0}`$. Let
$`Q=\mathrm{ADD}_m(g;U)`$, and let E XOR the literal ell into
$`g_0`$. The chronological word

```math
F=Q^\dagger,\ E,\ Q,\ E
```

returns g and adds $`\ell e`$ to U modulo M: flipping the low bit
changes the numerical value of g by $`\ell(1-2g_0)`$. Let $`C_g`$ be
CNOT fanout from $`g_0`$ to all bits of U. Its arithmetic action is
$`U\mapsto eU-g_0`$, so

```math
C_g,\ F,\ C_g:\qquad
U\longmapsto e(eU-g_0+\ell e)-g_0=U+\ell\pmod M.
```

The last equality uses $`e^2=1`$ and $`(e+1)g_0=0`$. Thus two TTK
adders and Clifford gates implement the controlled increment with
$`T,G,D_T=O(m)`$, using exactly m arbitrary returned helper bits.
For a positive address literal E is CNOT into $`g_0`$; for a negative
literal it is that CNOT followed by X on $`g_0`$. The address is never
a target and never participates in a non-Clifford gate. This is an
all-input identity, including coherent helpers and address controls.

### A logarithmic-depth increment using two dirty bits

The following involution echo improves that increment's private workspace
and depth. Let $`S(U)=-U-1\pmod M`$ be bitwise complement and
$`R(U)=-U\pmod M`$ be modular negation. Both are involutions, and
chronological S followed by R is increment. Reserve two arbitrary dirty
bits d,e, disjoint from U and the public literal ell.

Implement $`R^d`$ by a CNOT from d to each counter bit, followed by
the d-controlled increment from
[Vandaele, Section 5, Corollary 7](https://arxiv.org/html/2603.12917v1#S5).
That completed increment may borrow d internally, but d is private here;
it returns d and the separate dirty helper e exactly. Write E for XOR of
ell into d. Then the chronological word

```math
F=R^d,\ E,\ R^d,\ E
```

acts on U as $`R^{d\oplus\ell}R^d=R^\ell`$ and returns both dirty
bits. Therefore the desired controlled increment is

```math
S^\ell,\ R^d,\ E,\ R^d,\ E:
\qquad U\longmapsto U+\ell\pmod M.
```

The order matters: placing $`S^\ell`$ last gives decrement instead.
For a positive literal, its public gates are CNOTs into U and d. For a
negative literal, append X on each corresponding private target. The
public address is never changed and participates only as a CNOT control.
All non-Clifford gates belong to private row work. This implements

```math
T,G=O(m),\qquad D_T=O(\log(m+2)),\qquad w_{\mathrm{extra}}=2.
```

For m equal to one use a single literal-controlled X directly. For all
other m, the two helper bits are sufficient; no minimality is claimed.
The full permutation identity, including literal phase, proves the
contract on arbitrary coherent dirty inputs and reference systems.
The [bounded native checks](../tests/test_readonly_dirty_increment.py)
use an explicitly slower serial controlled increment; the optimized
depth above is inherited from the completed source theorem.

### Add the Hamming weight without initializing a counter

For one conjunction of k literals $`\ell_i`$, allocate private
arbitrary m-bit registers $`a_1,\ldots,a_k,c`$. Form a balanced
reversible addition tree on the a registers: add each left representative
into its disjoint right representative and carry an unpaired block to
the next level. Its root contains $`\sum_i a_i\pmod M`$.
Add that root into c and reverse the entire tree. Denote this word by L;
its exact action is

```math
L:\quad c\longmapsto c+\sum_i a_i\pmod M,
\qquad (a_1,\ldots,a_k)\longmapsto(a_1,\ldots,a_k).
```

Each tree level consists of disjoint additions. Thus L has
$`T,G=O(km\log(m+2))`$ and
$`D_T=O(\log(k+1)\log^2(m+2))`$. Two private helper bits per a
register suffice for its controlled increment. The older m-bit helper
reservation remains a valid conservative bound for m at least two;
the one-bit increment needs no helper. The modular adders themselves
need no helper.

Let J increment every $`a_i`$ by its literal $`\ell_i`$ in parallel.
The chronological translation echo

```math
A=L^\dagger,\ J,\ L,\ J^\dagger
```

has the exact action $`c\mapsto c+s\pmod M`$, where
$`s=\sum_i\ell_i`$, and returns every a register and increment helper.
Indeed, the two additions into c are minus the old sum of a registers
and plus their sum after the literal increments. Unknown initial offsets
cancel by arithmetic modulo M. No counter stores s by itself, and no
zero sector is assumed.

### Remove the remaining counter offset by cyclic routing

Allocate M arbitrary dirty selector bits V. Define P to move the old
bit $`V_j`$ to position $`j-1\pmod M`$, and let $`B_c=P^c`$
act on V while preserving c. For each bit $`c_j`$, apply the fixed
permutation $`P^{2^j}`$ controlled by that bit. Every fixed permutation
is a product of two involutions, each a matching of disjoint swaps:
on a cycle indexed by t, the reflections $`t\mapsto-t`$ and
$`t\mapsto1-t`$ compose to one step. Choose their order for the
desired orientation, and omit fixed points. The existing shared-control
Fredkin schedule implements each matching in at most four T layers.
Consequently $`B_c`$ has $`D_T\le8m`$, $`T,G=O(Mm)`$,
and no additional work.

The chronological word

```math
R=A,\ B_c,\ A^\dagger,\ B_c^\dagger
```

acts as $`P^s`$ on V and returns every counter and increment helper.
The first rotation sees $`c+s`$, the second sees the restored c, so
their product is $`P^{c+s}P^{-c}=P^s`$. Since $`P^s`$ moves old position
s to zero, the word

```math
K=R,\ X_{V_0},\ R^\dagger
```

flips exactly $`V_s`$, returning all other registers. Finally use

```math
\mathrm{CX}(V_k;y),\ K,\ \mathrm{CX}(V_k;y),\ K^\dagger.
```

This toggles the arbitrary output y by $`[s=k]`$ and returns every
bit of V. Because $`0\le s\le k\lt M`$, this is precisely the
conjunction of all literals, with no modular alias. Every identity is
literal on every computational-basis input; the native Toffoli and
Fredkin words introduce no residual phases. Linearity therefore gives
the complete arbitrary-input and reference-return statement.

### Parallel equality rows and the resource ledger

For row r, choose $`\ell_i=[x_i=r_i]`$, take y to be $`Y_r`$,
and allocate its own counters, selector, and increment helpers. All H rows
are private except for the original address. Its only occurrences are
in J or its actual inverse, as CNOT controls in E. Negative literals add
only X gates on private helper bits. During J, every adder and native T gate
acts within one private counter/helper block. Hence all non-Clifford
layers run on disjoint wires across rows and literal positions. Shared
CNOT fanout has a charged gate count and is not assumed to have constant
total depth. This proves $`D_T(J)=O(m)`$ for the whole family, without
an H or k serialization factor; the address may be coherent.

The live helper width per row is at most $`2(k+1)m+M`$: k input
counters, one accumulator, their private increment helpers, and the selector.
All are simultaneously reserved; none is borrowed from a live output
or another row. L uses $`2(k-1)+1`$ additions. A, R, K, and the
final output echo use only a fixed number of their constituent words
and actual inverses. Thus a single row has

```math
T,G=O(km\log(m+2)+Mm),\qquad
D_T=O\!\left(\log(k+1)\log^2(m+2)+m\right).
```

Multiplying counts and width by H, while retaining the batched depth,
proves the opening bounds. The case k equal to zero is X and needs no
helper. No total Clifford-depth bound is being inferred from T-depth.

### Complete-frame consequence

The baseline indicator satisfies the polynomial-overhead interface of
Section 5. Its [masked-sum refinement](DIRTY_SUM_COMPRESSION.md) uses
$`T,G,w=O(2^k(k+1)^3)`$ and $`f(k)=O(\chi(k))`$;
the same majorant $`P(t)=C(t+1)^3`$ remains sufficient. At fixed
accuracy and a sufficiently large $`b=\Theta(\sqrt N)`$, the refined
composition gives one two-clean complete real-frame circuit with

```math
T=O(\sqrt N),\qquad G=O(N),\qquad
D_T=O\!\left(n\chi(n)\right).
```

It retains the same full-isometry error bound and optimal-order
worst-case T-count. Choosing between this circuit and the earlier routed
one gives depth $`O(\min\{n^2,n\chi(n)\})`$ at the same
sufficient square-root-scale dirty width. This is an asymptotic
improvement, not a practical crossover estimate or an optimal-depth
theorem. The high-precision complete-frame endpoint remains open.

The [baseline counter checks](../tests/test_counter_dirty_indicator.py) audit
both modular-adder actions, signed increments, actual inverses, arbitrary
dirty offsets, cyclic orientation, and shared-address scheduling. Full
counter fixtures emit the linear TTK adder. The reduced RV macro is
checked separately; its optimized ladder depth is imported analytically,
not emitted or inferred from the serial macro fixtures. The general
resource bounds and Hopf composition are analytic arguments above. Separate
[compression checks](../tests/test_dirty_sum_interfaces.py) audit the
round-based sum, its helper offsets, actual inverses, and column schedule.
The [pipeline checks](../tests/test_pipelined_dirty_sum.py) audit the
deferred parity forests, doubling blocks, full sum echo, and release
schedule. This refinement uses linear TTK arithmetic inside each block
and for readout; its depth improvement comes from the pipeline proof.

## 7. Attribution and evidence

The [partial-batch extension](BATCHED_DIRTY_LOOKUP.md) retains the same
bank echo while reusing fewer indicator wires. At fixed accuracy it
gives $`T=O(\sqrt N+N/b)`$ and
$`D_T=O(N\log(b+2)/b^2+n^2)`$ for complete real frames with two
clean qubits and $`b\ge17(L+n+7)`$. Two additional dirty helpers are
reserved during each live batch guard; they are not borrowed from its
occupied bank or indicator storage.
The [amortized construction](AMORTIZED_DIRTY_LOOKUP.md) further reduces
that fixed-accuracy depth to $`O(N/b^2+n^2)`$ at the same threshold,
using returned dirty traversal selectors and one outer indicator echo.
Its variable-accuracy composition gives
$`D_T=O(NL/b^2+nL+n^2)`$ with
$`T=O(\sqrt{NL}+NL/b+nL)`$ at that threshold; both resources match
their lower bounds for $`b\le\sqrt{NL/(nL+n^2)}`$ when eligible.

Parallel dirty indicators and the separation of selector parallelism from
word-bank count have primary precedent in Low, Kliuchnikov, and Schaeffer,
*Trading T gates for dirty qubits in state preparation and unitary synthesis*,
[arXiv:1812.00954v2, Appendix C](https://arxiv.org/html/1812.00954v2).
The scratch-free indicator of Section 1 is a direct conjugation of that established routing
primitive; no new general lookup tradeoff is claimed. Literal Fredkin
words and the complete dirty-register ledger make it suitable for the
repository's full-input contract. The local result is the sharper
simultaneous count/depth composition for complete two-clean real Hopf
frames. Its predicates use the separately credited borrowed-MCX schedule
in the linked depth proof.

[Focused checks](../tests/test_parallel_dirty_lookup.py) test the selected
bit, actual inverse, literal native phases, scratch-free width, and
completed lookup return. Finite checks do not establish the asymptotic
theorem or optimal T-depth. A matching depth lower bound and the
high-precision linear endpoint remain separate open questions.
