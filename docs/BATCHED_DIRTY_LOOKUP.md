# A batched dirty-workspace tradeoff for fixed-accuracy Hopf frames

[Depth model and routing](T_DEPTH_COMPILER.md) · [Full indicator](PARALLEL_DIRTY_LOOKUP.md) · [Complete-frame source](OPERATOR_SOURCE_COMPILER.md)

A dirty indicator need not hold every high-address selector at once.
Processing bounded groups of selectors gives a depth tradeoff below the
workspace threshold of the full-indicator construction, while retaining
the available optimal-order T-count on the same circuit.

**Fixed-accuracy theorem.** Let $`n\ge1`$, $`N=2^n`$, and fix
$`0\lt\eta\le1/64`$ independently of n. Put
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and $`B_0=L+n+7`$.
For every $`b\ge17B_0`$, the prescribed complete real Hopf frame has
a coherent Clifford+T implementation using two initialized compiler flags
and at most b arbitrary dirty qubits, with

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T=O\!\left(\sqrt N+\frac Nb\right),\qquad G=O(N),
\qquad
D_T=O\!\left(\frac{N\log_2(b+2)}{b^2}+n^2\right).
```

The constants may depend on the fixed accuracy. All three bounds hold
for the same circuit. The complete-input norm includes flag leakage,
arbitrary dirty inputs, and their reference correlations. No resets,
measurements, initialized indicators, supplied catalysts, or quantum
memory oracles are used. The coefficient tables are classical; their
coherent queries and every inverse are implemented and charged.

The T-count is worst-case optimal in order under these workspace
conditions. The T-depth matches the available lower bound within a
logarithmic factor in the small-workspace regimes identified in Section 5.
The additive $`n^2`$ term leaves a larger gap at large workspace. This
is neither a general optimal T-depth theorem nor a total elementary-depth
bound. The threshold 17 is a literal sufficient reservation for this
proof, not a necessary threshold or a practical crossover estimate.

## 1. A guarded indicator on arbitrary dirty bits

Let an address have the form $`(c,u)`$, where u selects one of s
indicator bits Y and c selects one of K chunks. Assume s and K are
powers of two. For a fixed classical chunk index t, let
$`\mathcal R_u`$ be the existing exact router on the s one-bit banks Y.
It moves the original $`Y_u`$ to position zero. Apply, chronologically,

```math
\mathcal R_u,\quad
\mathrm{MCX}_{[c=t]\to Y_0},\quad
\mathcal R_u^\dagger.
```

Call this word $`I_t`$. For every basis input it implements

```math
|c,u,Y\rangle\longmapsto
|c,u,Y\oplus[c=t]e_u\rangle.
```

Indeed, on the selected chunk the central gate flips the routed bit;
undoing that same permutation flips the original position u. On other
chunks the central gate is identity and the router cancels. The actual
inverse router is required. The native Fredkin router and central MCX
are literal exact circuits, so the identity extends to superpositions,
arbitrary dirty indicator contents, and entangled references without
address-dependent phases.

Use the two-dirty-helper MCX construction of Khattar and Gidney,
[*Rise of conditionally clean ancillae for optimizing quantum circuits*,
Section 5.4](https://arxiv.org/html/2407.17966v1#S5.SS4), with literal
exact Toffoli decompositions. A conjunction of p literals then has
$`O(p)`$ T and Clifford count and $`O(\log(p+1))`$ T-depth for
$`p\ge1`$; bounded p is handled directly. Negative literals add X gates.
The two helpers are arbitrary, disjoint from every address, Y bit, and
word bank, and returned exactly before the inverse router. When K=1,
the central gate is an unconditional X and neither helper is used.

Put $`p=\log_2K`$. The guarded indicator therefore has

```math
T(I_t),G(I_t)=O(s+p),\qquad
D_T(I_t)=O\!\left(1+\log_2s+\log_2(p+1)\right).
```

Its only work beyond Y is the two returned MCX helpers. They are not
borrowed from a live word bank. The address is preserved throughout the
completed indicator, and no clean predicate flag is required.

## 2. Chunked high-address loading and the completed query

Consider a Q-row, m-bit XOR table $`f`$, with $`Q=2^r`$ and $`m\ge1`$.
A non-power-of-two table may first be padded by zero rows. Choose powers
of two $`1\le\mu\le Q`$ and $`1\le s\le H=Q/\mu`$, and set
$`K=H/s`$. Split the address as $`(c,u,v)`$, with c the chunk, u its
row within the chunk, and v the bank selector. These address fields are
disjoint and unchanged by the completed query.

Let Z consist of $`\mu m`$ arbitrary dirty bank bits. For each t,
form a classical binary matrix $`A_t`$ with s columns and $`\mu m`$
rows. Column u concatenates the table words $`f(t,u,v)`$ over all
bank choices v. The Clifford linear map

```math
C_{A_t}:\quad (Y,Z)\longmapsto(Y,Z\oplus A_tY)
```

uses at most $`s\mu m`$ CNOTs. The chronological echo

```math
C_{A_t},\quad I_t,\quad C_{A_t},\quad I_t^\dagger
```

changes Z by $`[c=t]A_te_u`$ and restores Y and both helpers.
For an initial indicator word y, the two bank contributions are
$`A_ty`$ and $`A_t(y\oplus[c=t]e_u)`$; their unwanted parts cancel.
This proof does not assume that Z or Y is zero. Concatenate these K
completed echoes to obtain a loader $`\mathcal L`$ for the complete
high address. The same Y and helper pool is reused between chunks.

Next use the established completed SelectSwap word. Let $`\mathcal R_v`$
route the selected m-bit bank to position zero, and let C copy that bank
by CNOTs into the separate, arbitrary m-bit query output. Apply

```math
\mathcal L,\ \mathcal R_v,\ C,\ \mathcal R_v^\dagger,\
\mathcal L^\dagger,\ \mathcal R_v,\ C,\ \mathcal R_v^\dagger.
```

The output receives the loaded selected bank and its original contents,
leaving exactly $`f(c,u,v)`$. Every bank, indicator bit, and helper
returns exactly. The loader inverse reverses the entire native loader,
including chunk order, and occurs only after routing is undone. It never
uses the separate output as work. All identities are on the complete
input space. Zero table rows remain exactly inactive even if internal
routing and echoes act on those address sectors.

The resulting exact query has

```math
\begin{aligned}
T_{\rm query}&=O\!\left(H+Kp+\mu m\right),\\
G_{\rm query}&=O\!\left(Qm+H+Kp+\mu m\right),\\
D_{T,\rm query}&=O\!\left(
K[1+\log_2s+\log_2(p+1)]+\log_2\mu\right).
\end{aligned}
```

Clifford linear maps can have substantial elementary depth; only their
gate count is used here. The resource bounds do not hide their execution
in a quantum table oracle.

| Additional live register | Dirty width | Lifetime and return |
|---|---|---|
| Word banks Z | $`\mu m`$ | Live throughout the completed query; returned exactly |
| Indicator Y | s | Reused between completed chunk echoes; returned exactly |
| Dedicated MCX helpers | 2 | Used only by the central conjunction gate; returned before inverse routing |

The m-bit query output, address, source core, and any previously reserved
flags or selectors are outside this table. The literal additional dirty
width is

```math
\mu m+s+2.
```

## 3. Allocation at a prescribed extra width

Let B be the available additional dirty width, excluding the query
output and all existing reservations. A convenient sufficient condition is

```math
B\ge16(m+r+1),\qquad r=\log_2Q.
```

For $`Q\ge m`$, choose mu as the largest power of two at most

```math
\min\!\left\{\sqrt{Q/m},\frac{B}{4m}\right\}.
```

This minimum is at least one. For $`Q\lt m`$, choose $`\mu=1`$.
With $`H=Q/\mu`$, choose s as the largest power of two at most
$`\min\{H,B/4\}`$. Both choices divide Q, and

```math
\mu m\le B/4,\qquad s\le B/4,\qquad
\mu m+s+2\le B.
```

The count cap on mu prevents excess workspace from forcing more word
banks than the count-efficient choice. Unused wires remain untouched.

If K=1, there is no high-chunk predicate. If $`K\gt1`$, rounding gives
$`s\gt B/8`$, and $`p=\log_2K\le r`$. Hence

```math
Kp=\frac Hs p\le\frac{8Hr}{B}=O(H).
```

For $`Q\ge m`$, power-of-two rounding also gives

```math
H=O\!\left(\sqrt{Qm}+\frac{Qm}{B}\right),
\qquad \mu m\le\sqrt{Qm}.
```

For $`Q\lt m`$, we instead have $`H=Q\lt m\le B/16`$,
so s=H and K=1. These cases prove

```math
T_{\rm query}=O\!\left(\sqrt{Qm}+m+\frac{Qm}{B}\right),
\qquad G_{\rm query}=O(Qm).
```

For the Clifford bound, $`H\le Q`$ and $`\mu m\le Qm`$;
the predicate count is absorbed as above. When Q=1, the sole row may
simply be emitted as an X word on the output, with zero T-count,
$`O(m)`$ Clifford gates, and no additional work.

To bound depth, first consider K=1. The indicator and bank routing
contribute $`O(1+\log_2H+\log_2\mu)=O(\log_2(Q+1))`$.
For $`K\gt1`$, the case $`Q\lt m`$ is impossible. If the mu cap is
$`B/(4m)`$, then $`\mu\gt B/(8m)`$, and
$`K\lt64Qm/B^2`$. If the mu cap is $`\sqrt{Q/m}`$, then
$`H\lt2\sqrt{Qm}`$. Since $`H\gt s\gt B/8`$, this implies
$`\sqrt{Qm}/B\gt1/16`$, so again $`K=O(Qm/B^2)`$.
Moreover $`s\le B/4`$ and $`p\le r\le B/16`$, so one chunk
costs $`O(\log_2(B+2))`$ T-depth. Thus the same exact query has

```math
D_{T,\rm query}=O\!\left(
\log_2(Q+1)+\frac{Qm}{B^2}\log_2(B+2)\right).
```

This query lemma allows arbitrary m. The Hopf theorem below is specialized
to fixed accuracy; no variable-precision complete-frame theorem is asserted
by that specialization.

## 4. Fixed-accuracy complete-frame composition

Use the existing two-flag layerwise complete-frame compiler and retain
its coefficient tables, source words, suffix echo, amplification, and
error allocation. At depth $`d=0,\ldots,n-1`$, it has a fixed number
of mask queries with

```math
Q_d=2^{d+2},\qquad r_d=d+2,\qquad m_d=L+n-d+4.
```

The old source core, selectors, and separate dirty suffix control occupy
exactly $`B_0=L+n+7`$ wires. Its two initialized flags are separate.
The crucial simultaneous reservation is

```math
m_d+r_d+1=B_0,\qquad
B=b-B_0\ge16B_0.
```

The allocation of Section 3 therefore fits every query, with all banks,
indicators, and both dedicated MCX helpers outside the occupied base.
These extra registers are reused between completed queries and layers.
Between queries, the suffix predicates use two returned extra-pool wires
as in the [existing depth proof](T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
The source and its actual controlled versions remain on their stated
base registers. No initialized flag is reused as dirty work.

Because every replaced query has the same literal all-input unitary,
the complete-input approximation proof is unchanged. In particular,
earlier leakage is propagated by actual unitaries, the suffix echo still
cancels on every input, and no good-block projection or reset is inserted
between amplification uses. The original layer errors telescope to
less than $`2^{-L}\le\eta`$ for the full prescribed frame.

The existing geometric sums give

```math
\sum_d Q_dm_d=O(NL),\qquad
\sum_d\sqrt{Q_dm_d}=O(\sqrt{NL}),\qquad
\sum_dm_d=O(nL+n^2).
```

At fixed L, the last expression is $`O(n^2)`$, which is
$`O(\sqrt N)`$. The fixed number of source calls per layer contributes
that same count and serial depth. Suffix predicates have $`O(n^2)`$
total T and Clifford count and $`O(n\log(n+1))`$ depth; the two-flag
reflections have constant cost per use. These polynomial counts are also
absorbed into $`O(\sqrt N)`$ and $`O(N)`$, respectively.

As $`B=\Theta(b)`$ and the number of queries per layer is fixed,
the query and source bounds give, on one schedule,

```math
T=O\!\left(\sqrt N+\frac Nb\right),\qquad G=O(N),
```

```math
\begin{aligned}
D_T
&=O\!\left(\sum_d\log_2(Q_d+1)
+\frac{\log_2(B+2)}{B^2}\sum_dQ_dm_d+n^2\right)\\
&=O\!\left(n^2+\frac{N\log_2(b+2)}{b^2}\right).
\end{aligned}
```

No cap on the available b is necessary: the per-query bank cap already
avoids needless T-count, and full indicators automatically appear when
the useful pool is large enough. The norm transfers to the existing
complete-frame Hopf-QBP approximation and inverse guarantees. Observable
cost and the number of executions remain separate charges.

## 5. Optimal-order count and the remaining depth gap

With two clean flags and $`b\ge17B_0`$, physical width
$`q=n+2+b`$ is $`\Theta(b)`$. At fixed accuracy the retained
worst-case complete-frame T-count lower bound therefore gives

```math
T^\star=\Omega\!\left(\sqrt N+\frac Nb\right).
```

The theorem attains this order while also satisfying its displayed
depth bound. Dividing the count lower bound by q gives

```math
D_T^\star=\Omega\!\left(\frac{\sqrt N}{b}+\frac N{b^2}\right).
```

There is also the existing constant worst-case lower bound. For
$`b\le\sqrt N`$, the term $`N/b^2`$ dominates the displayed
expression. In particular, when the interval is nonempty,

```math
17B_0\le b\le\frac{\sqrt N}{n}
```

implies $`n^2\le N/b^2`$. The achieved depth is then within
$`O(\log(b+2))`$ of this lower bound. At $`b=\Theta(n)`$ with a
sufficient fixed constant, the concrete asymptotic consequences are

```math
T^\star=\Theta(N/n),\qquad
\Omega(N/n^2)\le D_T^\star\le O(N\log(n+2)/n^2).
```

The upper count and depth still belong to the same construction.
At $`b=\Theta(\sqrt N)`$ the retained upper depth is $`O(n^2)`$
and the lower bound is only constant; no logarithmic matching claim is
made there. The exact source-depth certificate in its restricted
Majorana architecture is not added to these unrestricted lower bounds.
The high-precision constant-clean complete-frame endpoint is unchanged.

## 6. Attribution and finite evidence

Low, Kliuchnikov, and Schaeffer,
[*Trading T gates for dirty qubits in state preparation and unitary synthesis*](https://arxiv.org/html/1812.00954v2),
supply the SelectSwap and dirty-indicator framework. The
[routed indicator](PARALLEL_DIRTY_LOOKUP.md) supplies the literal router
conjugation used here; Khattar and Gidney supply the two-dirty-helper
conjunction implementation. The new local composition batches these
indicators under an explicit live-width budget and applies the resulting
schedule to prescribed complete real Hopf frames. No claim of a new
general lookup primitive or unrestricted optimal T-depth is made.

[Focused finite checks](../tests/test_batched_dirty_lookup.py) cover the
guarded router, multiple chunk echoes, actual inverse orientation,
literal native phases, and arbitrary dirty inputs. Their small parameter
choices need only satisfy the literal width ledger; the stronger
sufficient reservation in Section 3 is used for the uniform asymptotic
bounds. Finite examples do not establish those bounds or the lower-bound
comparison; the preceding proofs do.
