# Parallel dirty lookup while retaining the T-count bound

[T-depth model and routing](T_DEPTH_COMPILER.md) · [Grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md) · [Research status](OPEN_PROBLEM.md)

Additional dirty workspace can parallelize the selector computation without
increasing the count-efficient number of word banks. This gives a simultaneous
count and depth guarantee for the prescribed complete real Hopf frame.

**Theorem.** Write $`N=2^n`$, $`n\ge1`$,
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

At fixed accuracy, a sufficiently large $`b=\Theta(\sqrt N)`$ therefore
permits **the same circuit** to have

```math
T=O(\sqrt N),\qquad D_T=O(n^2),\qquad G=O(N).
```

The existing worst-case lower bound makes this T-count optimal in order.
The T-depth upper bound is not proved optimal, and does not bound total
elementary depth. The selected allocation $`L=N,b=N+n+7`$ is not covered
by the new sufficient-width hypothesis; its linear T-count endpoint remains
open.

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

## 5. Attribution and evidence

Parallel dirty indicators and the separation of selector parallelism from
word-bank count have primary precedent in Low, Kliuchnikov, and Schaeffer,
*Trading T gates for dirty qubits in state preparation and unitary synthesis*,
[arXiv:1812.00954v2, Appendix C](https://arxiv.org/html/1812.00954v2).
The indicator above is a direct conjugation of that established routing
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
