# A blocked bilinear query and the fixed-accuracy width frontier

[Amortized dirty traversal](AMORTIZED_DIRTY_LOOKUP.md) · [Bilinear queries](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction) · [Chunked indicators](CHUNKED_DIRTY_INDICATOR.md) · [Unary group compiler](UNARY_PHASE_GRADIENT.md)

A bilinear table query can process several address blocks without routing
an output bank. Its selected middle operation uses a dirty traversal of
the high address. A constant-depth controlled bilinear operation supplies
each leaf. The two low-address indicators are computed only a constant
number of times around that whole traversal. Their chunk sizes can then
vary with the layer's distance from the leaves.

Combining this exact query with the conditional unary group compiler
gives the following fixed-accuracy theorem. Fix $`0\lt\eta\le1/64`$,
set $`N=2^n`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
and $`B_0=L+n+7`$. For every $`n\ge1`$ and $`b\ge17B_0`$,
every prescribed complete real Hopf frame W has one coherent Clifford+T
circuit V with two external clean flags and at most b arbitrary dirty
qubits satisfying

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
```

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad
G=O_\eta(N),\qquad
D_T=O_\eta\!\left(\frac N{b^2}+n\right).
```

All three bounds describe the same circuit. The error includes every
work return and arbitrary reference correlations. Every query and
inverse below is a literal full-input native circuit. There are no
measurements, resets, QRAM, supplied phase states, or extra external
initialized work. T-depth allows arbitrary intervening Clifford
circuits, whose elementary gate count remains included in G.

The T-count has the established optimal worst-case order. When the
interval is nonempty, count and depth both match their existing lower
bounds for

```math
17B_0\le b\le\sqrt{N/n},\qquad
T^\star=\Theta_\eta(N/b),\qquad
D_T^\star=\Theta_\eta(N/b^2).
```

The additive n in the upper bound is not an unrestricted depth lower
bound. Depth optimality at larger widths and the variable-accuracy
complete-frame endpoint remain open. Constants may depend on the fixed
accuracy in this chapter. The subsequent
[uniform-precision refinement](UNIFORM_PRECISION_DEPTH.md) exposes that
dependence and uses a rectangular allocation of the same query to prove
absolute-constant depth O(NL/b²+nL), with a stronger slowly growing
precision corollary.

## 1. A controlled bilinear output with one arbitrary dirty helper

Let Y and X be disjoint registers of H and J arbitrary bits, let z be
one arbitrary output bit, and let h be a separate arbitrary control.
For a fixed binary matrix $`D\in\mathbb F_2^{H\times J}`$, we need

```math
(h,Y,X,z)\longmapsto(h,Y,X,z\oplus hY^{\mathsf T}DX).
```

Reserve one arbitrary dirty helper e, disjoint from these registers.
Binary elimination, with actual inverse coordinate changes, puts the
bilinear form into $`\sum_{i=1}^{\rho}Y_iX_i`$, where
$`\rho=\operatorname{rank}(D)`$. The orientation is the same as in
the [existing bilinear lemma](PARALLEL_DIRTY_LOOKUP.md#exact-bilinear-oracle):
if $`PDQ=J_\rho`$, use $`Y'=P^{-\mathsf T}Y`$ and
$`X'=Q^{-1}X`$. These basis changes and their inverses have
$`O(\rho(H+J))`$ Clifford gates and require no extra wire.

In these coordinates, apply H to z and define the exact diagonal word
and exact Toffoli

```math
D_e=\prod_{i=1}^{\rho}\mathrm{CCZ}(e,Y_i,X_i),\qquad
E=\mathrm{CCX}(h,z;e).
```

Execute chronologically

```math
D_e,\quad E,\quad D_e^\dagger,\quad E^\dagger,
```

then reverse the Hadamard and the coordinate changes. Every dagger is
the actual reversed native word. Writing
$`b(Y,X)=\bigoplus_iY_iX_i`$, the two diagonal actions have combined
phase

```math
(-1)^{e b(Y,X)+(e\oplus hz)b(Y,X)}
=(-1)^{hz b(Y,X)}.
```

The final E inverse restores e. Conjugating z by Hadamards therefore
implements the required controlled XOR, with every helper returned.
This proof permits arbitrary h, z, e, Y, and X. In particular, for h
equal to zero the full completed word is identity on arbitrary dirty
work. The computational-basis identity, including literal phase,
extends to arbitrary superpositions and reference correlations.

The [shared-control phase schedule](T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
implements $`D_e`$ with at most four T layers and
$`6\rho+(\rho\bmod2)`$ T gates. Its control e may be arbitrary;
no initialized copies are used. Each E has seven T gates and four
T layers. Thus, for $`\rho>0`$, the controlled bilinear output has

```math
T\le12\rho+2(\rho\bmod2)+14,\qquad
D_T\le16,\qquad G=O(\rho(H+J)).
```

For rho equal to zero, emit identity. For an m-bit output Z with
matrices $`D_1,\ldots,D_m`$, process the output bits sequentially,
restoring both coordinate changes after every output. Reuse the same
returned helper e. Writing $`R_D=\sum_j\operatorname{rank}(D_j)`$,
the complete controlled bilinear oracle has

```math
T=O(R_D),\qquad D_T=O(m),\qquad G=O(mHJ).
```

Completed output-bit words preserve Y, X, h, and e. Intermediate
coordinate changes need not preserve individual Y or X bits.

## 2. Select a bilinear block with dirty traversal

Split the r-bit table address into three disjoint fields $`(c,a,b)`$,
with $`K=2^p`$ high blocks, H low-a values, and J low-b values, so
$`Q=KHJ=2^r`$. For output bit j in block t, let
$`D_{t,j}\in\mathbb F_2^{H\times J}`$ be its classical table matrix.
Define the desired selected bilinear map

```math
F_c:(Y,X,Z)\longmapsto
\left(Y,X,\left(Z_j\oplus Y^{\mathsf T}D_{c,j}X\right)_{j=1}^m\right).
```

Reserve p arbitrary dirty traversal selectors and the separate helper
e from Section 1. Use the existing
[two-pass dirty traversal](AMORTIZED_DIRTY_LOOKUP.md#2-selecting-a-shear-with-dirty-unary-traversal),
with each leaf now the controlled bilinear oracle for that block.
Its terminal selector is the control h in Section 1.

For completeness, at leaf t write $`\lambda_i=[c_i=t_i]`$ and let
$`v_1=w_1\oplus\lambda_1`$,
$`v_i=w_i\oplus v_{i-1}\lambda_i`$. Its expanded terminal value is
$`v_p=[c=t]\oplus h_t(c,w)`$. One traversal inserts the leaf
operation controlled by $`v_p`$; the second omits the first selector
update, retaining the control $`h_t`$, and is applied by its actual
inverse. Each traversal restores all its selectors. For p equal to
one, the second traversal still visits its two leaves; only the root
updates are omitted.

Every completed leaf preserves Y, X, and the selectors. All completed
bilinear maps commute, since they preserve Y and X and XOR their
contributions into Z. The two traversals consequently cancel every
unknown dirty-selector contribution and retain only block c. Each
leaf also returns e. This proves F on the full input space, including
arbitrary selector and helper states. For p equal to zero, use the
existing unguarded bilinear oracle directly.

The complete depth-first traversals have $`O(K)`$ nodes. Each leaf
has at most $`m\min(H,J)`$ total rank. Therefore

```math
T(F)=O(Km\min(H,J)),\qquad
D_T(F)=O(Km),\qquad G(F)=O(KmHJ)=O(Qm).
```

The additional work is exactly the p selectors and one dirty helper;
Y and X are separately allocated indicator registers. The Clifford
ledger includes all coordinate changes for every block and every
output bit. No classical matrix transformation is treated as a free
quantum oracle.

## 3. Two indicator echoes complete the exact query

Let $`I_a:Y\mapsto Y\oplus e_a`$ and
$`I_b:X\mapsto X\oplus e_b`$ be exact dirty indicators, preserving
their addresses and returning all their own helpers. Execute

```math
F_c,\ I_a,\ F_c^\dagger,\ I_b,\quad
F_c,\ I_a^\dagger,\ F_c^\dagger,\ I_b^\dagger.
```

For every fixed c and output bit, its four bilinear contributions sum
over $`\mathbb F_2`$ to

```math
Y^{\mathsf T}D_cX
+(Y+e_a)^{\mathsf T}D_cX
+(Y+e_a)^{\mathsf T}D_c(X+e_b)
+Y^{\mathsf T}D_c(X+e_b)
=D_c[a,b].
```

Thus the output is XORed with exactly its selected table word. Both
indicators, their helpers, all traversal selectors, and e return
exactly. Every address is preserved. A zero table row gives exact
identity even if the internal words act nontrivially on that sector.
The query works for arbitrary initial outputs and dirty work, and
for coherent addresses and references.

If the indicators' respective resources are $`T_a,G_a,D_a,w_a`$
and $`T_b,G_b,D_b,w_b`$, where widths include their outputs, a
sufficient simultaneous query width beyond the original address and
output is

```math
w_a+w_b+p+1.
```

This conservative allocation keeps both private indicator pools
disjoint. Reusing returned helper pools could improve its constant but
is unnecessary. The exact resource ledger is at most four completed F
calls, two a-indicators, and two b-indicators, counting actual inverses
at the same cost. In particular,

```math
T=O(Km\min(H,J))+2T_a+2T_b,
```

```math
G=O(Qm)+2G_a+2G_b,
\qquad D_T=O(Km)+2D_a+2D_b.
```

There is no output-bank router. This is the existing bilinear echo
with a new, fully charged selected middle operation.

## 4. Allocate a truncated balanced query at a prescribed width

Use the [chunked dirty indicator](CHUNKED_DIRTY_INDICATOR.md) with a
positive integer chunk cap a. Put $`P_a=(a+2)^3`$. Fix an absolute
$`c_I\ge1`$ large enough that an s-output indicator, including its
output word, has T-count, Clifford count, and dirty width each at most
$`c_IsP_a`$. This holds even if its address length is below a, since
its actual chunks are truncated to that length. A zero-bit indicator
is simply X on its sole output.

Let B be additional dirty width, excluding the query address, output,
and every existing compiler reservation. The following sufficient
condition makes the allocation direct:

```math
B\ge16c_I(P_a+r+1).
```

For $`Q\ge2`$, choose s as the largest power of two at most

```math
\min\!\left\{2^{\lfloor r/2\rfloor},\frac{B}{8c_IP_a}\right\}.
```

Take $`H=J=s`$ and $`K=Q/s^2`$. The unassigned address bits are
the high-block field c. Both H and J divide Q, and K is an integer
power of two; if r is odd the balanced count cap leaves at least one
high-block bit. Each indicator gets its own private dirty pool. Their
combined width is at most B divided by four, and the traversal selectors
plus e fit the stated remaining allocation. Thus the whole query fits B.
For $`Q=1`$, emit its sole row directly as a Clifford X word.

Power-of-two rounding gives

```math
\frac1s=O\!\left(\frac1{\sqrt Q}+\frac{P_a}{B}\right),
\qquad
K=O\!\left(1+\frac{QP_a^2}{B^2}\right).
```

Since $`s\le\sqrt Q`$, Section 3 yields the uniform upper bounds

```math
T=O\!\left(m\sqrt Q+\frac{QmP_a}{B}+\sqrt QP_a\right),
```

```math
G=O(Qm+\sqrt QP_a),
```

```math
D_T=O\!\left(
 m\left[1+\frac{QP_a^2}{B^2}\right]
 +\left\lceil\frac{\log_2s}{\min\{a,\log_2s\}}\right\rceil
       \log_2(\min\{a,\log_2s\}+2)\right).
```

The indicator-depth term is defined as zero when s equals one. This
standalone query bound does not assert optimal dependence on arbitrary
word precision m. In the fixed-accuracy Hopf composition, its extra
precision factors occur in geometrically weighted layer sums.

## 5. Fixed-accuracy composition at every eligible width

Use the existing theorem's literal threshold $`b\ge17B_0`$ from
the opening statement. Separate a small-width branch from the new
construction. If

```math
b\le\sqrt N/n,
```

then the [amortized fixed-accuracy theorem](AMORTIZED_DIRTY_LOOKUP.md#fixed-accuracy-corollary)
already gives
$`D_T=O_\eta(N/b^2+n^2)=O_\eta(N/b^2)`$ with the required
same-circuit T and Clifford counts. No query or source is changed on
this branch.

For the remainder suppose $`b>\sqrt N/n`$ and n is sufficiently
large, with the threshold depending only on eta and fixed circuit
constants. Use the unary compiler's
[early groups and error allocation](UNARY_PHASE_GRADIENT.md#7-complete-frame-theorem-at-fixed-accuracy).
In detail, choose the least power of two
$`q\ge4\pi\sqrt{2n}/\eta`$, put $`\rho=\log_2 3`$ and
$`\delta=\eta/(8n)`$, and retain its sufficient local constant A
and cutoff

```math
K_*=\left\lceil16A(q^\rho+q\log_2q+\log_2q+1)\right\rceil,
\qquad C=16A.
```

At remaining height $`k>K_*`$, use
$`g=\lfloor\log_2(k/(Cq))\rfloor`$. Their program width is
$`w=q(2^g-1)\le k/C`$. These groups fit their conditionally zero
logical suffixes, number $`O_\eta(n/\log(n+2))`$, and have
$`O_\eta(n)`$ total internal, prefetch, and predicate depth.
Their local gate counts are polynomial in n. The cutoff satisfies

```math
K_*=O_\eta(n^{\rho/2}+\sqrt n\log(n+2)),
\qquad K_*\log(n+2)=o_\eta(n).
```

The early parallel query helper pools still fit the smaller available
dirty width. Their peak is bounded by

```math
O\!\left(\sqrt N(n+2)^3
                  \sum_{k>K_*}k2^{-k/2}\right)
=o_\eta(\sqrt N/n).
```

Their T-count is $`O_\eta(\sqrt N)`$ and their Clifford count is
$`O_\eta(N)`$, exactly as in the square-root-width proof. These
pools are separate from the conditional logical suffix and are returned
before every group begins.

Set

```math
L'=\max\{6,\lceil\log_2(4/\eta)\rceil\},
\qquad B'_0=L'+n+7,
\qquad B=b-B'_0-2.
```

Since $`L'\le L+2`$, the original literal threshold ensures
$`B\ge b/2`$. On the large-width branch,
$`B>\sqrt N/(2n)`$, so this pool is larger than any fixed polynomial
in n for sufficiently large n. The source base, two predicate helpers,
and two external clean flags remain separately reserved. The predicate
helpers are arbitrary dirty inputs and return exactly.

### A width-constrained query on the entire tail

After the last unary group, the remaining height is
$`K'\le K_*`$. For every $`1\le k\le K'`$, use the old
full-input operator source and literal amplification with its capped
precision

```math
m_k=L'+4+\min\{k,\lceil\log_2(8n)\rceil\}
\le m=O_\eta(\log(n+2)).
```

Its exact table has $`Q_k=4N2^{-k}`$ rows and
$`r_k=n-k+2\le n+1`$ address bits. Use Section 4 with

```math
a_k=\min\{n+2,2^{\lfloor k/12\rfloor}\},
\qquad P_k=(a_k+2)^3.
```

Both useful estimates hold:

```math
P_k\le(n+4)^3,\qquad P_k\le27\,2^{k/4}.
```

The first ensures the sufficient width condition
$`B\ge16c_I(P_k+r_k+1)`$ for every tail query on this branch,
including the smallest k. The second makes the following resource sums
converge. No source or predicate borrows a live query register.

The query's indicator depth is

```math
O\!\left(n(k+1)2^{-k/12}+\log(n+2)\right).
```

Indeed, when the untruncated exponential cap is below the indicator
address length, its reciprocal is within a fixed factor of
$`2^{-k/12}`$ and its logarithm is $`O(k+1)`$. Otherwise the
indicator has one chunk and logarithmic address-length depth. This
argument also covers s equal to one by the Clifford convention above.
Thus each completed query obeys

```math
T_k=O\!\left(m_k\sqrt{Q_k}
       +\frac{Q_km_kP_k}{B}+\sqrt{Q_k}P_k\right),
```

```math
G_k=O(Q_km_k+\sqrt{Q_k}P_k),
```

```math
D_{T,k}=O\!\left(m_k+
       \frac{Q_km_kP_k^2}{B^2}
       +n(k+1)2^{-k/12}+\log(n+2)\right).
```

The bounded number of queries and actual inverses per amplified stage
only changes constants. These are exact full-input replacements for
the old table words, so the source's accepted block, entire rejected
action, and full error certificate are unchanged.

### Summing the tail resources

Use $`m_k\le L'+4+k`$ for the weighted count and chunk-depth
terms. With $`Q_k=4N2^{-k}`$ and $`P_k\le27\,2^{k/4}`$,

```math
\sum_{k=1}^{K'}m_k\sqrt{Q_k}=O_\eta(\sqrt N),\qquad
\sum_{k=1}^{K'}\sqrt{Q_k}P_k=O(\sqrt N),
```

```math
\sum_{k=1}^{K'}Q_km_kP_k=O_\eta(N),\qquad
\sum_{k=1}^{K'}Q_km_kP_k^2=O_\eta(N),\qquad
\sum_{k=1}^{K'}Q_km_k=O_\eta(N).
```

For example, the P-squared sum is bounded by a constant times
$`N\sum_{k\ge1}(L'+4+k)2^{-k/2}`$. Hence the whole tail has

```math
T=O_\eta\!\left(\sqrt N+\frac NB\right),\qquad G=O_\eta(N).
```

For the unweighted terms retain the cap $`m_k\le m`$ rather
than the looser bound linear in k. The sum
$`\sum_{k\ge1}(k+1)2^{-k/12}`$ is finite, so

```math
\sum_{k=1}^{K'}D_{T,k}
=O_\eta\!\left(\frac N{B^2}+n+K_*\log(n+2)\right)
=O_\eta\!\left(\frac N{b^2}+n\right).
```

The old full-input source words, reflections, and suffix predicates add
$`O_\eta(K_*\log(n+2))=o_\eta(n)`$ depth and polynomial gate
counts. Every query fits the extra pool B by construction; this pool is
reused only after its exact return. The early-group resources already
fit the same b. Consequently all three advertised bounds hold on the
same circuit.

### Error, the literal threshold, and the matching interval

The unary group's unchanged error argument charges at most eta divided
by four for rounding all early angles, using the full-frame angular
square-sum bound; at most eta divided by four for their native source
preparations and actual inverses; and at most eta divided by four for
the remaining old source layers with accuracy parameter eta divided
by four. This gives total initialized-isometry error at most eta.
All work return and reference correlations are included. Exact table
replacement introduces no additional approximation error, and no
intermediate source is assumed reset.

The finitely many n below the large-n thresholds use the original
amortized theorem at accuracy eta and literal width $`b\ge17B_0`$.
Its depth $`O_\eta(N/b^2+n^2)`$ becomes
$`O_\eta(N/b^2+n)`$ over that finite set by increasing only the
eta-dependent constant. The same fallback already handles the
small-width branch for every n. Thus the global statement retains
$`17(L+n+7)`$ despite using $`L'`$ in the large-width construction;
it does not silently impose $`17(L'+n+7)`$ on the whole theorem.

Finally, the inherited complete-frame lower bounds give
$`T^\star=\Omega_\eta(\sqrt N+N/b)`$ and
$`D_T^\star=\Omega_\eta(N/b^2)`$ under this width threshold.
For $`b\le\sqrt{N/n}`$, the n term in the new depth bound is at
most $`N/b^2`$, and the square-root count term is at most $`N/b`$.
This proves the simultaneous matching interval in the opening
statement. It makes no lower-bound claim for the additive n outside
that interval.

## 6. Proof and evidence boundary

The controlled bilinear echo, high-block selection, and two-indicator
identity are exact algebraic circuit proofs. They use the established
shared-control native phase schedule, rank reduction, dirty traversal,
and chunked dirty indicators; none is an assumed quantum table oracle.
The global theorem combines these identities with the already proved
conditional unary group, angular error bound, and original full-input
source certificate. Its count and depth estimates follow from the
explicit reservations and geometric sums above, not finite-size fits.

The [five bounded checks](../tests/test_blocked_bilinear_lookup.py) verify
native controlled bilinear leaves on every input, one complete eight-wire
native selected-block word, and its actual inverse. Exact Boolean
polynomials check complete queries, including arbitrary indicator words,
selector stacks, outputs, and helpers. Emitted schedules check literal
T-count, T-depth, and disjoint T targets. Negative controls expose a
missing helper inverse, omitted traversal, and incorrect block labels.
The complete-query checks use a reduced controlled-bilinear gate and the
existing routed indicator; they do not emit asymptotic chunked indicators
or a variable-size complete-frame compiler. Their exact scope is recorded
in [verification](VERIFICATION.md).

The [source map](reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries)
and [related-work comparison](../research/RELATED_WORK.md#20-blocked-bilinear-queries-and-the-dirty-width-tradeoff-3-october-2026)
attribute the dirty traversal, phase synthesis, and bilinear framework.
The additional result is their charged selected-block interface and
fixed-accuracy width composition, not generic priority for phase
polarization or dirty cancellation. No unrestricted depth lower bound
is inferred from the finite checks.
