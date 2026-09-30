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
D_T=O\!\left(\min\{nL+n^3,\ L\ell_*(n)+n^4\}\right).
```

Only two external clean qubits are used. The complete-input error includes
clean leakage and arbitrary dirty inputs with references. The coefficient
tables remain classically specified and their quantum lookups are fully
charged. The sufficient constant C is not a practical crossover estimate.

At fixed accuracy, a sufficiently large $`b=\Theta(\sqrt N)`$ therefore
permits **the same circuit** to have

```math
T=O(\sqrt N),\qquad D_T=O(n^3),\qquad G=O(N).
```

The existing worst-case lower bound makes this T-count optimal in order.
The T-depth upper bound is not proved optimal, and does not bound total
elementary depth. The selected allocation $`L=N,b=N+n+7`$ is not covered
by the new sufficient-width hypothesis; its linear T-count endpoint remains
open.

## 1. An exact dirty bilinear echo

All sequences in this chapter are chronological. Work in Boolean XOR and
AND. Suppose dirty vectors U and V can be toggled by fixed functions u and v
of an unchanged address. Let F XOR their outer product into distinct target
bits Y. The sequence

```math
F,\ U\mathrel{\oplus}=u,\ F,\ V\mathrel{\oplus}=v,\
F,\ U\mathrel{\oplus}=u,\ F,\ V\mathrel{\oplus}=v
```

returns U and V and changes Y by

```math
UV^{\mathsf T}\oplus(U\oplus u)V^{\mathsf T}
\oplus(U\oplus u)(V\oplus v)^{\mathsf T}
\oplus U(V\oplus v)^{\mathsf T}=uv^{\mathsf T}.
```

No initial value of U or V is assumed. These are exact permutation
identities on every basis input, so they extend to arbitrary superpositions
and reference correlations with no input-dependent phase.

### Parallel implementation of the outer product

Let U have p bits and V have q bits, and put $`H=pq`$. Borrow two dirty
bits $`A_{ij},B_{ij}`$ for every pair. Four batches of disjoint exact
Toffolis into $`Y_{ij}`$, interleaved with fanout toggles, give

| Batch | Contents used as the two controls |
|---|---|
| 1 | $`A_{ij},B_{ij}`$ |
| 2 | $`A_{ij}\oplus U_i,B_{ij}`$ |
| 3 | $`A_{ij}\oplus U_i,B_{ij}\oplus V_j`$ |
| 4 | $`A_{ij},B_{ij}\oplus V_j`$ |

Toggle A between batches 1–2 and 3–4; toggle B between batches 2–3 and
after batch 4. The same bilinear identity leaves precisely
$`Y_{ij}\mathrel{\oplus}=U_iV_j`$ and restores every A and B.
Each Toffoli is a literal exact Clifford+T gate word; no faulty-sign
decomposition is used. Within a batch the three-wire supports are disjoint.

A shared-control CNOT fanout on t arbitrary targets needs no clean copies.
Choose a binary CNOT tree P on those targets that maps the first unit
vector to the all-ones vector. The chronological word
$`P^\dagger,\mathrm{CNOT}_{c,1},P`$ toggles all targets by c,
using $`O(t)`$ CNOTs and $`O(1+\log t)`$ elementary depth. The separate
rows of A, or columns of B, have disjoint supports and run in parallel.
Consequently F uses

```math
4H\ \text{Toffolis},\qquad T,G=O(H),\qquad D_T=O(1),
\qquad D_{\rm elem}=O(1+\log H),
```

and exactly $`2H`$ additional dirty bits. Actual reversed-word inverses
are available with the same resources.

## 2. A literal dirty indicator

For a k-bit address x, let $`I_k`$ implement

```math
|x,Y,w\rangle\longmapsto|x,Y\oplus e_x,w\rangle,
\qquad H=2^k,
```

where $`e_x`$ is the length-H unit vector and both Y and w are arbitrary.
For $`k=0`$, use X on the sole output. For $`k=1`$, use X on output
zero and CNOT from x into each of the two outputs.

For $`k\ge2`$, split x into r and s bits, with
$`r=\lfloor k/2\rfloor`$, $`s=\lceil k/2\rceil`$.
Borrow vectors U and V of lengths $`p=2^r`$, $`q=2^s`$.
Apply the bilinear echo from Section 1, using two recursive calls to
$`I_r`$ to toggle U and two to $`I_s`$ to toggle V. The four F calls
target the H-bit output grid. The net change is
$`e_{x_{\rm hi}}e_{x_{\rm lo}}^{\mathsf T}=e_x`$.
Every recursive call and outer-product batch restores its borrowed work.
Use actual reversed-word inverses for the second recursive toggles;
their logical action equals the forward indicator XOR.

Let K(k) count dirty scratch beyond the H outputs. The $`2H`$ work bits
needed during F can be reused for recursive scratch, since F and the
recursive calls occur at different times. Thus

```math
K(k)=2^r+2^s+\max\{2H,K(r),K(s)\}\le3H.
```

The base cases need no scratch. The inequality follows inductively from
$`2^r+2^s\le H`$ and $`3\cdot2^s\le2H`$ for $`k\ge2`$.
Including Y, the indicator therefore uses at most $`4H`$ dirty bits.

If A(k) counts its Toffolis, the explicit construction gives

```math
A(k)=2A(r)+2A(s)+16\cdot2^k,
\qquad A(0)=A(1)=0.
```

Both T-count and Clifford count are $`O(2^k)`$. For sufficiently large
k, the child sizes total only $`O(2^{k/2})`$; their constant multiplicity
is absorbed by the top-level $`2^k`$ term. The finitely many smaller k
set the fixed constant. Elementary depth satisfies

```math
D(k)\le2D(r)+2D(s)+O(k)=O((k+1)^2).
```

The same upper bound holds for T-depth. In particular, the shared-control
fanouts are charged at their actual Clifford depth; they are not assumed
to be physical constant-depth operations.

## 3. Count-efficient whole-word queries

Consider a Q-row, m-bit XOR table, padding Q to a power of two with zero
rows if needed. Let $`1\le\mu\le Q`$ be a power-of-two bank count and
put $`H=Q/\mu`$. The high address specifies the concatenation of mu
table words. Its length is $`\mu m`$.

Reserve the dirty bank word Z of length $`\mu m`$, an H-bit dirty
indicator Y, and at most $`3H`$ indicator scratch. Let A be the classical
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
$`O((1+\log H)^2)`$ T-depth. A serial schedule of the CNOT matrix
has depth at most $`Qm`$; no stronger claim about that depth is needed.

Use this loader in the [completed SelectSwap query](T_DEPTH_COMPILER.md#2-exact-lookup-with-depth-optimized-banks): load, route, copy to the
arbitrary output, undo routing and loading, then route, copy, and undo
routing again. The original selected bank contents cancel. The exact
resources are bounded by

```math
T_{\rm query}=O(Q/\mu+\mu m),\qquad G_{\rm query}=O(Qm),
```

```math
D_{T,\rm query}=O\!\left((1+\log(Q/\mu))^2+\log\mu\right).
```

An elementary-depth upper schedule is
$`O(Qm+(1+\log H)^2+\mu m)`$: the linear maps can be serial, and
the routing uses only $`O(\mu m)`$ elementary gates. In particular,
low T-depth does not imply equally low Clifford or total depth.

| Live registers | Dirty width | Lifetime and return |
|---|---|---|
| Existing source, selectors, and helpers | $`B_0`$ | Reserved throughout; existing source return is charged in the frame error |
| Additional word banks | $`\mu m`$ | Live throughout the completed query; returned exactly |
| Indicator output | $`H`$ | Live during each loader; returned exactly |
| Indicator scratch | at most $`3H`$ | Reused between outer products and recursive calls; returned exactly |

Coefficient-query outputs already lie in the base reservation;
coarse-symbol queries reuse their existing private buffers. Neither
requires an additional initialized output. All new work fits when

```math
B_0+\mu m+4Q/\mu\le b.
```

Zero completed table rows give the identity on inactive sectors, even
though the internal indicator and router may act there. There are no
new source calls, logical clean-bit assumptions, or supplied quantum
oracles in this replacement.

Choose mu as the largest power of two at most $`\sqrt{Q/m}`$ when
$`Q\ge m`$, and choose $`\mu=1`$ otherwise. For the former case,

```math
\mu m+4Q/\mu\le9\sqrt{Qm};
```

for the latter it is at most $`5m`$. Thus a sufficient extra dirty pool
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

Each coefficient query now has T-depth $`O(n^2)`$. There are at most
n groups and $`O(s_g^2)`$ coarse-symbol queries per group. Since
$`\sum_gs_g^2\le n^2`$, their total query depth is $`O(n^4)`$.
All other schedules are exactly those already charged in the
[depth proof](T_DEPTH_COMPILER.md#3-composition-with-the-grouped-full-frame-compiler).
This gives $`D_T=O(L\ell_*(n)+n^4)`$.

Applying the same query replacement to the older layerwise compiler gives

```math
T=O(\sqrt{NL}+nL),\qquad
D_T=O(nL+n^3),\qquad G=O(NL).
```

Select the schedule with the smaller displayed depth expression. If the
layerwise expression is smaller, then
$`nL\le L\ell_*(n)+n^4`$, and $`n^4=O(\sqrt N)`$ absorbs the
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
Its proof allows faulty-sign Toffolis. The exact bilinear echoes above
give a self-contained literal implementation with a stated dirty-width
constant, sufficient for the repository's stricter complete-input contract.
They do not claim a new general lookup tradeoff or improve that source's
elementary-depth theorem. The local result is the simultaneous count/depth
composition for complete two-clean real Hopf frames.

[Focused checks](../tests/test_parallel_dirty_lookup.py) test the cancellation,
native phases, returned work, and schedule supports. Finite checks do not
establish the asymptotic theorem or optimal T-depth. A matching depth lower
bound and the high-precision linear endpoint remain separate open questions.
