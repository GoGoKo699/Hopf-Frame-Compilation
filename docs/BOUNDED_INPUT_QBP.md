# Bounded-input costs for Hopf gradients

[Task theorem](STATE_BASED_QBP_THEOREM.md) · [Task-cost comparison](QBP_COST_COMPARISON.md) · [Complex gradient protocol](COMPLEX_COARSE_QBP.md) · [Classical source tables](OPERATOR_SOURCE_COMPILER.md#10-classical-table-construction)

This note fixes an explicit classical input and execution model. In this
model the observable is a supplied Pauli list, so a classical algorithm
can apply it to an amplitude vector. That access is stronger than an
uninterpreted controlled-observable oracle. A comparison must charge it
consistently on both sides.

## 1. Explicit inputs and the cost model

Let $`N=2^n`$, $`n\ge1`$. Supply the $`N-1`$ real Hopf angles and,
for a complex chart, N real leaf-phase representatives as dyadic numbers
in $`[-8,8]`$, each with at most $`B\ge1`$ fractional bits. Constantly many sign
and integer bits suffice for these bounded inputs. No amplitude-to-angle
conversion, phase unwrapping, or transcendental range reduction is included
in the input contract. A real chart omits the phase array.

Supply an explicit list

```math
H=\sum_{r=1}^{M}c_rQ_r,\qquad
Q_r\in\{I,X,Y,Z\}^{\otimes n},\qquad
c_r\in[-1,1]\cap2^{-B}\mathbb Z,
\qquad \Lambda=\sum_r|c_r|\gt0.
```

Coefficients with fewer fractional bits are padded to this common
denominator. A negative Pauli sign is absorbed into its coefficient.
Repeated terms and cancellations are allowed. The dyadic representation
here is ordinary fixed point, not a succinct exponent encoding whose
expanded denominator length is left uncharged. The list length M and input
precision B are unrestricted unless a later statement adds a promise.
Reading the parameter tuple and observable takes input length

```math
I_{\rm in}=O(NB+M(n+B)),
```

in addition to the supplied accuracy and confidence specifications. The
actual description length is charged when it is smaller or larger because
of a chosen encoding. Computing Lambda exactly uses integer sums with
$`B+\lceil\log_2(M+1)\rceil+O(1)`$ bits.

Request every original raw gradient coordinate within
$`0\lt\varepsilon_\infty\le\Lambda`$, and retain the comparison's

```math
K=\max\{6,\lceil\log_2(80\Lambda/\varepsilon_\infty)\rceil\},
\qquad P=\max\{n,K\}.
```

A probabilistic algorithm also receives $`0\lt\delta\lt1`$. The
deterministic baseline below has no failure probability and needs no
confidence amplification. Observable coefficients and the parameter tuple
are fixed while taking the derivatives; they are not differentiated.

Use a classical random-access bit-cost model with schoolbook integer
arithmetic. Addition on w-bit integers costs $`O(w)`$ and multiplication
or division costs $`O(w^2)`$. Accessing n-bit basis indices and reading
their Pauli descriptors is charged. Explicit classical tables and circuit
words must be constructed and stored; they are not free oracles. A fixed
compiled program may be reused across executions, but sequential emission
of its instructions and each executed elementary quantum gate remain
charged operations. T-count and Clifford count are still separately
reported; converting these ledgers into physical time requires explicit
gate-cost assumptions. No QRAM or coherent unit-cost access to the input
arrays is assumed. Ordinary classical arrays do not provide such quantum
access.

## 2. A deterministic classical all-gradient baseline

The signed-permutation Pauli action and the reverse chain rule are standard
classical tools. For the general arithmetic derivative-complexity framework,
see [Baur and Strassen, *The complexity of partial derivatives*, Theoretical
Computer Science 22 (1983), 317–330](https://doi.org/10.1016/0304-3975(83)90110-X).
The formulas below specialize that reverse calculation to the Hopf tree
and give a direct precision and input-cost proof; no priority is claimed
for reverse-mode differentiation or classical state-vector contraction.

### 2.1 Amplitudes, signed Pauli permutations, and one reverse pass

Write $`c_j=\cos\theta_j`$, $`s_j=\sin\theta_j`$ for the tree
factors; these subscripts label nodes and are distinct from the observable
coefficients above. Starting at $`a_{\rm root}=1`$, compute

```math
a_{j0}=a_jc_j,\qquad a_{j1}=a_js_j.
```

At a leaf x set $`r_x=a_x`$ and
$`\psi_x=e^{i\varphi_x}r_x`$. In a real chart use phase factor one.
There are $`O(N)`$ forward operations after trigonometric evaluation.
The original supplied phases can be used directly; no common-phase gauge
or native coarse circuit is needed for this classical calculation.

For Pauli string $`Q_r`$, let $`u_r`$ mark its X or Y positions and
$`z_r`$ mark its Z or Y positions. Let $`y_r`$ count its Y positions.
In the same bit order as the computational basis,

```math
Q_r|x\rangle=i^{y_r}(-1)^{z_r\cdot x}|x\oplus u_r\rangle,
```

where the dot product is modulo two. Thus compute $`v=H\psi`$ by
permuting amplitudes and multiplying by signs or powers of i, accumulating
the supplied exact coefficients. This takes $`O(MN)`$ scalar arithmetic
operations, plus at most $`O(MNn)`$ elementary bit work for naive index and
parity calculation. It does not construct an $`N\times N`$ matrix.

Initialize real leaf weights
$`\beta_x=\operatorname{Re}(e^{-i\varphi_x}v_x)`$. Traverse the
original real tree upward:

```math
\beta_j=c_j\beta_{j0}+s_j\beta_{j1},\qquad
g_{\theta_j}=2a_j(-s_j\beta_{j0}+c_j\beta_{j1}).
```

At the leaves also output

```math
g_{\varphi_x}=2\operatorname{Im}(\overline{\psi_x}v_x).
```

These identities are $`2\operatorname{Re}(J^\dagger H\psi)`$ and
the original direct phase derivative. The reverse pass holds v and the
phase factors fixed; it is the analytic Hopf contraction, not a derivative
through an approximate evaluator. It costs $`O(N)`$ operations and
outputs all $`N-1`$ magnitude coordinates and, when present, all N phase
coordinates. No dense Jacobian, division by an amplitude, metric inverse,
or exact-zero test is used. The formulas include arbitrary signs and
singular tuples. Zero-support phase derivatives vanish automatically.

### 2.2 Certified digits without unstable intermediate divisions

A deliberately conservative implementation makes the bit bound explicit:
round the trigonometric factors once to dyadics, then perform the forward
pass, Pauli accumulation, and reverse pass with exact rational arithmetic.
All denominators remain powers of two. Only the final answer is rounded.
This avoids assumptions about floating-point cancellation or conditioning.

Choose R large enough that every real sine/cosine factor and every complex
unit phase is certified to absolute error at most $`e=2^{-R}`$. Clamp
the real factors to $`[-1,1]`$; this cannot increase their error. A complex
phase can be rounded with a fixed guard margin in its real and imaginary
parts. For example it suffices to take

```math
R\ge\max\!\left\{6,
 \left\lceil\frac n2\right\rceil+
 \left\lceil\log_2
  \frac{12(n+1)\Lambda}{\varepsilon_\infty}\right\rceil\right\}.
```

Rounding to R fractional bits with a fixed additional guard changes this
choice only by an absolute constant. In particular
$`R=O(K+n+\log(n+1))`$.

To verify the error bound, each leaf amplitude and each nonzero raw
magnitude derivative is a signed product of n real factors and one unit
phase. Products of the clamped real factors have absolute value at most
one. Telescoping their perturbations gives component error at most
$`(n+1)e`$. Structural zeros in derivative vectors are retained exactly.
For the conceptual exact and approximate state and derivative vectors,

```math
\|\widetilde\psi-\psi\|_2\le E,\qquad
\|\widetilde d_j-d_j\|_2\le E,
\qquad E=(n+1)\sqrt N\,e.
```

These vectors are used only to prove the error, not stored as a Jacobian.
The exact state has norm one and every raw magnitude derivative has norm
at most one. The phase derivative is $`i|x\rangle\langle x|\psi`$,
so its approximation obeys the same bound. Since
$`\|H\|\le\Lambda`$, the exact arithmetic result satisfies, for every
magnitude or phase coordinate,

```math
\begin{aligned}
|\widetilde g_j-g_j|
&\le2\Lambda\left(E(1+E)+E\right)\\
&\le6\Lambda E\le\varepsilon_\infty/2.
\end{aligned}
```

Here the chosen R ensures $`E\le\varepsilon_\infty/(12\Lambda)\le1/12`$.
Final dyadic rounding to absolute error at most
$`\varepsilon_\infty/2`$ therefore yields the requested deterministic
coordinate guarantee. The same bound covers singular angles; no inverse
gap enters the required precision.

### 2.3 A conservative schoolbook bit bound

For bounded dyadic inputs in $`[-8,8]`$, the
[explicit Taylor procedure](OPERATOR_SOURCE_COMPILER.md#10-classical-table-construction)
computes a certified sine/cosine or complex phase to R bits in
$`O(R^4)`$ bit operations after reading the input, with fixed guard
margins. It uses exact integer arithmetic and a certified Taylor tail, not
an assumed elementary-function oracle. There are $`O(N)`$ such values.

Every final polynomial term contains at most $`2n+2`$ rounded
trigonometric factors and one observable coefficient. The tree partial
sums and Pauli accumulation have the same degree bound or a smaller one.
After aligning their power-of-two denominators, it suffices to use

```math
W_{\rm bit}=O\!\left(B+nR+\log(M+1)\right)
```

bits per rational intermediate, including sign and integer range.
For a coarse bound one may take
$`B+(2n+4)(R+1)+\lceil\log_2(M+1)\rceil+10`$ bits, enlarging the
absolute constant for the chosen fixed trigonometric guard margin.
This also covers exact computation of Lambda and final rounding. No
normalization by Lambda is needed in the arithmetic; its exact value
only sets the error budget.

Consequently a deterministic implementation has bit complexity

```math
\boxed{
O\!\left(I_{\rm in}+N R^4+
 MN[W_{\rm bit}^2+n]\right).
}
```

This includes forward/reverse arithmetic and output rounding since
$`M\ge1`$. Accuracy-input processing and its description length are
also charged. Auxiliary storage is
$`O(NW_{\rm bit}+R^2)`$ bits beyond the supplied input list; terms can
be processed sequentially. Writing the rounded output costs its actual
bit length, bounded by the displayed estimate. With exact coefficients
and singular-safe polynomial contractions, no Monte Carlo shots or
native synthesis preprocessing enter this baseline.

Identical Pauli labels can also be merged exactly before the contraction.
There are only $`4^n=N^2`$ labels, so the number of nonzero resulting
terms satisfies

```math
M_{\rm eff}\le\min\{M,N^2\}.
```

The common denominator remains $`2^B`$ and merged numerators need at most
$`B+\lceil\log_2(M+1)\rceil+O(1)`$ bits. Deterministic comparison
sorting and accumulation, for example, cost at most

```math
O\!\left(M\log(M+1)[n+B+\log(M+1)]\right)
```

bit operations, with its input and sorting storage charged. Afterwards
replace the factor M in the contraction count by $`M_{\rm eff}`$,
while retaining the coefficient-bit growth in $`W_{\rm bit}`$. This
gives at most $`O(N^3)`$ scalar contraction operations after aggregation,
with the same certified polynomial bit overhead. It is an available
baseline, not a claim of optimal classical complexity.

Aggregation is optional preprocessing available to both quantum and
classical routes. If it makes H zero, all gradients are zero. Otherwise
it may lower the coefficient norm: one may recompute
$`\Lambda_{\rm eff}`$, the common accuracy parameter K, and the sampling
distribution from that same merged representation for both routes.
Neither a cancellation-inflated unmerged norm nor its conservative shot
count should be used to claim an advantage. When applying a theorem stated
for $`\varepsilon_\infty\le\Lambda_{\rm eff}`$, a larger requested
error can simply be replaced by the smaller valid target
$`\min(\varepsilon_\infty,\Lambda_{\rm eff})`$. The unmerged formulas
above remain valid without this optional preprocessing.

### 2.4 A second classical baseline when the Pauli list is long

The explicit access model also permits classical term-only sampling.
Prepare the same amplitude data classically, sample r with probability
$`|c_r|/\Lambda`$, and consider the whole coordinate vector

```math
G_r=\Lambda\operatorname{sgn}(c_r)
       \nabla\langle\psi|Q_r|\psi\rangle.
```

Its mean is the requested gradient. At a magnitude depth, the derivative
columns are orthonormal marker columns multiplied by the real incoming
amplitudes. Projection onto those columns and $`|a_j|\le1`$ give
$`\|G_r^{(d)}\|_2\le2\Lambda`$. The entire phase block has norm at
most $`2\Lambda`$ as well, since it is
$`2\Lambda\operatorname{sgn}(c_r)\operatorname{Im}(\overline\psi\odot Q_r\psi)`$.
Thus the retained
Hilbert-space concentration bound gives the same
$`O(\Lambda^2\varepsilon_\infty^{-2}[1+\log((n+1)/\delta)])`$
sufficient sampling order. The common
conservative allocation in [the cost comparison](QBP_COST_COMPARISON.md#1-common-precision-workspace-and-execution-ledger)
covers both its sampling error and the deterministic rounding allocation.
There is no additional quantum branch or leaf measurement noise in this
classical comparator.

Do not compute and retain one gradient vector per draw. First collect
$`h_r=\operatorname{sgn}(c_r)\,\#\{\text{draws of }r\}`$ and form

```math
\widehat v=\frac\Lambda S\sum_{r:h_r\ne0}h_rQ_r\psi.
```

Then run one reverse tree and one phase contraction. If D distinct terms
were sampled, $`D\le\min(M,S)`$, the arithmetic count is
$`O(ND+S+M)`$ after coefficient and trigonometric evaluation, plus index
operations. This is the exact average of the sampled gradient vectors.
The sampled effective operator has norm at most Lambda for every record,
so the same uniform trigonometric error bound applies pathwise. Tighten R
by a fixed constant if allocating only $`\varepsilon_\infty/4`$ to this
rounding. Counts use $`O(1+\log S)`$ bits; multiply the integer-weighted
sum first and apply the common rational factor $`\Lambda/S`$ at the
end. All additional integer multiplication, division, and output bits are
charged.

An exact term sampler is also explicit. Align $`|c_r|`$ over $`2^B`$,
compute integer cumulative sums, draw a uniform integer below their total
by binary rejection, and locate its interval by binary search. The
integer width is $`O(B+\log(M+1))`$; rejection uses fewer than two
trials in expectation. Prefix sums, input reading, random bits, searches,
and count increments are not free. For example, their expected sampling
bit cost is bounded by

```math
O\!\left(
 M[B+\log(M+1)]+S\{[B+\log(M+1)]
 [1+\log(M+1)]+\log(S+1)\}\right).
```

The amplitude arithmetic can use the preceding exact-rational scheme with
an additional $`O(\log(S+1))`$ bits for counts and the final division.
This gives a fully charged polynomial-precision comparator for large M.
If evaluating all M terms is cheaper, use the deterministic algorithm
instead. These alternatives depend on the explicit Pauli list and do not
supply classical access to an arbitrary unknown quantum observable.

## 3. Which native searches are polynomial in the input parameters

### 3.1 Exhaustive search uses a length theorem, not a synthesis-time assumption

The imported [F5 one-qubit approximation lemma](SOURCE_MAP.md#5-fault-tolerant-sources-and-contribution-boundaries),
GKW Lemma 2.3, supplies a Clifford+T word of length
$`O(p+1)`$ approximating a determinant-one one-qubit target to
$`2^{-p}`$. Its Clifford count has the same order. Use this statement
only for existence and length; it is not an unconditional polynomial-time
algorithm for finding such a word.

There is nevertheless a deterministic, certified search costing

```math
2^{O(p)}\operatorname{poly}(B+p+1)
```

bit operations for the bounded-angle targets needed here. Enumerate
constant-alphabet words in increasing length. Evaluate each candidate and
the target with a fixed guard margin, and accept only when a rational
upper enclosure for their Frobenius distance is below the requested
error. A word promised at a sufficiently smaller constant error has a
strict acceptance margin, so no exact comparison at a rounding boundary
is required. An inconclusive candidate is skipped. The existence lemma
ensures acceptance by length $`O(p+1)`$ without requiring the search
algorithm to know its implicit constant.

Native matrix entries belong to a fixed algebraic number ring and are
certifiably evaluable in time polynomial in word length and requested
digits. Bounded target angles have the certified trigonometric evaluator
of Section 2. Thus each candidate test is polynomial-time, while the
number of tested words is $`2^{O(p)}`$. This preserves the logarithmic
quantum word length. The retained phase-cancelling wrappers and exact
reflection interpreters then produce their specified literal native
words with only constant-factor changes. Neither factoring nor a
number-theoretic runtime conjecture is used.

The distinction between coarse and fine p is decisive:

| Construction step | Native-search precision p | Consequence of this search bound |
|---|---|---|
| State route's actual coarse real rows and complex prefix-phase rows | $`O(n)`$ | $`N^{O(1)}`$ search time |
| Original grouped frame's coarse rows in a group of height s | $`O(s)\le O(n)`$ | $`N^{O(1)}`$ search time |
| Grouped residuals, the fixed deepest tail, and literal phase tables | No fine native-word search | Fine precision is supplied by explicit operator-source words and masks |
| Direct borrowed fine frame, real or gauge-fixed complex | $`O(K+n-d)`$ at depth d | In general $`2^{O(K+n)}`$ search time |
| The retained direct finite-size state fallback | $`O(P)`$ | In general $`2^{O(P)}`$ search time |

For the grouped coarse rows,
[the displayed coarse tolerance](CONDITIONAL_SUFFIX_COMPILER.md#3-a-small-coefficient-table-and-a-streamed-coarse-circuit)
is $`1/(4sK_g2^{s/2})`$ with $`K_g=O(s)`$, hence needs
$`O(s)`$ bits. Here $`K_g`$ is the group's padded type count, not
the accuracy parameter K. For state preparation, the actual coarse
distance is $`O(N^{-1/2})`$ and geometric layer allocation needs at
most $`O(n)`$ bits per row. Complex subtree means are rational averages
of the supplied bounded phases, and the prefix angles remain bounded.
The number of these coarse row targets is polynomial in N. Consequently
their combined exhaustive search is polynomial in N, even though its
fixed exponent need not be small.

### 3.2 Fine tables and explicit sources avoid a fine word search

For the original grouped compiler, construct the actual coarse words
first. Then evaluate the target group, its actual coarse matrix, and the
residual coefficients to certified precision. A conservative algorithm
may use dense matrices of dimension at most N for each group/prefix
piece. There are polynomially many such pieces. This is sufficient for
a polynomial classical bound; it is not an assertion that dense
processing is the preferred implementation.

All these matrices are finite products of explicitly supplied rotations
or recorded native gates. Computing entries to $`K+O(n)`$ guarded
fractional bits suffices after enlarging the constant for polynomially
many entrywise operations and the displayed coefficient normalizations.
The nonnegative phase split uses certified estimates and clipping, not
an exact sign decision. Rational mask rounding then follows
[the operator-source rule](OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients).
The weighted group ledger gives $`O(NK)`$ table bits. The fixed deepest
tail uses the same explicit sources, including when the entire tree is
below the fixed grouping threshold. Thus it introduces no exceptional
fine one-qubit search for small n. The
[one-clean literal complex diagonal](ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames)
also uses explicit sine/cosine masks, so it has the same conclusion.
Banking changes the exact query program, not these coefficient data.

For the state route, evaluate the target amplitudes, recorded coarse
blocks, and $`C^\dagger\psi`$ after the coarse search. A conservative
working precision is $`O(P+n+\log(n+1))`$ bits: the factor
$`\sqrt N`$ in the residual table and the completion's square-root
continuity are included by taking the coefficient enclosure to
$`2^{-2(P+10)-O(1)}`$. The
[direct algebraic residual-table construction](RESIDUAL_TABLE_PREPROCESSING.md#3-one-half-phase-used-twice)
uses certified sine/cosine pairs for three rotations;
it does not solve for a fine Euler angle or search a fine Euler grid.
Its rational arithmetic, square-root enclosures, clipping, and source
mask rounding have polynomial bit cost. Source width is $`O(P+n)`$,
and the complete state tables have $`O(NP)`$ bits because $`P\ge n`$.
The [rational certificates](RESIDUAL_TABLE_PREPROCESSING.md#4-fixed-precision-rational-certificates)
preserve the established rounding and full-isometry error allocation.

For $`n\ge6`$, these explicit fine source words fit the basic state
reservation $`b\ge P+n+7`$. For every $`n\ge1`$ they fit the
banked reservation $`b\ge2(P+n+7)`$: when $`n\le5`$, fix all
$`n+2`$ table-address bits into at most 128 sectors. There is one row
per sector and no free address bit. With $`q=P+10`$, the source core,
helper, and borrowed signal use $`P+13`$ dirty wires, which fit this
reservation. Its masks are known Pauli words, so no table-query bank is
needed in these sectors. Their number is constant, and the fine cost is
$`O(P)`$. This [explicit source implementation](RESIDUAL_TABLE_PREPROCESSING.md#6-small-systems-at-the-banked-reservation)
replaces the direct fine-word search at the stronger reservation. It does not change the
basic-budget fallback for $`n\le5`$.

## 4. A uniform construction and explicit-output bound

Use the input model of Section 1 and the literal workspace eligibility
conditions in [the task comparison](QBP_COST_COMPARISON.md#3-the-available-bounds-at-a-common-workspace).
There are deterministic classical algorithms that construct the
following native programs with their already proved quantum counts and
complete error/dirty-return contracts:

* the original grouped or eligible banked real frame, and its one-clean
  literal complex extension, for every $`n\ge1`$;
* the real or complex state program and coherent reference/target
  program for $`n\ge6`$ at the basic state reservation;
* those state programs for every $`n\ge1`$ at the banked reservation.

The spare-flag allocation for the one-clean complex frame and the
one-wire difference between the old and new bank thresholds remain
exactly as in that comparison. A preprocessing bound does not relax a
quantum workspace condition.

Put

```math
J=\left\lceil\log_2(n+b+K+3)\right\rceil.
```

Including input reading, coarse-word search, certified coefficient and
mask construction, native-program emission, and a library of the M
controlled Pauli terms, conservative bit-time bounds are

```math
\begin{aligned}
\mathcal P_{\rm grouped}
 &\le I_{\rm in}+N^{O(1)}\operatorname{poly}(B+K+J)
             +O(MnJ),\\
\mathcal P_{\rm state}
 &\le I_{\rm in}+N^{O(1)}\operatorname{poly}(B+P+J)
             +O(MnJ).
\end{aligned}
```

The constants and polynomial degrees are uniform and depend only on
the fixed native gate set and the retained compiler templates, not on
the supplied angles, K, B, M, or b. These intentionally conservative
bounds assert polynomial construction, not a practical synthesis time
or an explicit small exponent. Storage is polynomial as well; candidates
in the coarse search can be discarded sequentially.

To account for program output directly, the retained Clifford-plus-T
lengths are $`O(NK)`$ for the grouped/banked original words and
$`O(NP)`$ for the state words. Exact predicate, dirty-bank, source,
reflection, and inverse templates are finite constructive gate rules.
Expanding them gives, respectively,

```math
O(NKJ)\quad\hbox{and}\quad O(NPJ)
```

bits of elementary-gate instructions, including wire addresses. The
controlled-Pauli library adds $`O(MnJ)`$ bits. One may instead emit
only requested terms, charging their actual generation. The word is
classical output and can be reused; running it on each shot still costs
its executed gate count. Explicit static masks are compiled into these
gate lists, with no QRAM assumption. An unnecessarily huge declared
dirty pool can be left unused; retaining J also charges its binary
wire labels rather than treating their length as a constant.

This theorem does not make every quantum minimum in the comparison
polynomial-time constructible in K by the same proof. Direct fine
borrowed-frame synthesis, including the gauge-fixed complex corollary,
has the sufficient exhaustive-search bound
$`2^{O(K+n)}\operatorname{poly}(B+K+J)`$ after summing over rows.
It is polynomial in N when $`K=O(n)`$, but the bound can be exponential
in K outside that regime. The basic-budget state fallback at fixed
$`n\le5`$ has the analogous limitation. Faster synthesis algorithms
may improve these construction times; none is silently assumed here.
The explicit grouped and banked-state alternatives above already have
polynomial construction under their stated reservations.

## 5. What this establishes about complete costs

### 5.1 A constructive compiler statement, with its input restriction

The bounded-input construction closes a specific computational gap:
coefficient evaluation, the state completion, and program generation need
not be treated as unpriced searches at fine precision. The quantum T and
Clifford bounds remain those already proved, with their literal clean/dirty
reservations. Their classical construction is now polynomial for the
listed alternatives. The polynomial's coarse-search exponent is a fixed
gate-set constant, not a practical runtime estimate.

This statement does not apply to arbitrary supplied computable evaluators,
whose evaluation time can be unrestricted. It also does not turn the
controlled-observable interface into classical access to an unknown
quantum operation. The direct and term-sampled algorithms in Section 2
depend on the explicit Pauli description. A supplied general quantum
observable needs its own access and computational comparison.

Short native words and fast native-word search are distinct resources.
For context, [Ross–Selinger](https://arxiv.org/abs/1403.2975v3) give an
optimal algorithm with a factoring oracle; their efficient expected
runtime without that oracle is proved under a number-theoretic hypothesis.
No such hypothesis or oracle enters the coarse enumeration above. Nor do
we infer an unconditional polynomial-in-K search from a logarithmic
quantum word-length bound. The algebraic residual procedure is specific
to its displayed completion; it does not price generic SU(2) or U(2)
Euler preprocessing in the earlier complete-multiplexor theorem.

### 5.2 Uniform logical-gate costs remove the displayed T separation

Let $`\ell=1+\log_2^*(n+2)`$. In the high-precision banked comparison,
$`K\ge n^2`$ and $`P=K`$. The available per-execution T expressions
are

```math
Q=\sqrt{NK}+K+NK/b,\qquad
R=\sqrt{NK}+K\ell+NK/b.
```

Both retained Clifford expressions are $`O(NK)`$. Controlled Pauli
strings use $`O(n)`$ exact Clifford gates and no T gates, so their
execution cost is included at this order. If every elementary logical
gate has the same unit cost, adding the retained Clifford expressions
gives the same $`O(SNK)`$ quantum gate upper bound for both routes.
The T improvement alone therefore gives no separation between these
all-gate upper expressions. This is not a lower bound on Clifford count
or an optimal quantum gradient algorithm. Different physical prices for
T and Clifford operations require an explicit machine model.

Construction, term selection, classical records, and output bits are
additional charges. They consist of Section 4's preprocessing, Section
2.4's shared term sampler when it is used, and the applicable
[histogram/decoder bit ledger](QBP_COST_COMPARISON.md#5-classical-work-and-preprocessing).
Each sampled label and counter is processed; neither the $`S`$ records
nor the output coordinates disappear when T gates become cheaper.

### 5.3 High precision also has a deterministic classical comparator

Use the same Pauli representation and Lambda for both protocols. If
duplicates are merged, both may use the resulting representation and
recompute K from its Lambda before making this comparison. An identically
zero observable returns the zero gradient directly.

The diverging T-expression ratio requires
$`K\gg N/\ell^2`$ and $`b\gg N/\ell`$ along
$`\ell\to\infty`$. Since $`\ell^2=o(\sqrt N)`$, the first
condition implies $`K/\sqrt N\to\infty`$, hence eventually
$`N\le K^2`$. When M, B, and the accuracy-specification length are
polynomial in $`N+K`$, the deterministic bit bound in Section 2 is
consequently polynomial in K.
By contrast, the retained quantum protocols at the comparison's common
sufficient allocation execute

```math
S=\Theta\!\left(4^K[1+\ln((n+1)/\delta)]\right)
```

shots per stream, already requiring $`\Omega(S)`$ sequential executions.
Polynomial preprocessing does not convert that high-precision T saving
into an end-to-end advantage over this deterministic classical algorithm.
This concerns the displayed protocols and sufficient shot schedule in the
explicit Pauli model. It is not a sampling lower bound for every quantum
algorithm, and does not assert that the sufficient shot schedule is
optimal on every instance.

For an arbitrarily long supplied list, exact duplicate aggregation still
leaves at most $`N^2`$ distinct Paulis. Its ingestion, sorting, coefficient
growth, and any unusually long input precision remain charged on both
sides; they must not be dropped to obtain a stronger runtime comparison.
Classical term-only sampling is another available option, especially
before all distinct terms are worth applying. No universal advantage or
disadvantage is inferred for other observable-access models or accuracy
regimes.

### 5.4 Research boundary

The selected Hopf state-based construction has a complete quantum resource
ledger and now a uniform bounded-input construction. Its additional-bank
T improvement is a compiler result with a precise precision/workspace
range. The explicit Pauli case supplies useful classical comparators,
not an end-to-end speedup example.

The [consolidated state-based theorem](STATE_BASED_QBP_THEOREM.md) states
these contracts alongside their complete-frame boundary. A new algorithmic
advantage claim needs a concrete
observable-access model and an appropriate classical comparator; none is
selected by this audit. A full elementary fine-precision emitter remains
an implementation task. The separate constant-clean complete-frame
endpoint remains open, and neither question is settled by another special
finite fixture.
