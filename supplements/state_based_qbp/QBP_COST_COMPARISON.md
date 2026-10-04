# Comparing complete Hopf-gradient costs

[Task theorem](STATE_BASED_QBP_THEOREM.md) · [Original QBP ledger](../../docs/QBP_APPROXIMATION.md) · [Real coarse decoder](COARSE_FRAME_QBP.md) · [Complex coarse decoder](COMPLEX_COARSE_QBP.md) · [Research status](../../docs/OPEN_PROBLEM.md)

The state-only construction and the original frame protocol solve the same
raw-gradient task with different quantum programs and classical decoders.
A fair comparison uses the same accuracy, confidence, observable access,
and physical workspace. The state compiler's precision includes a dimension
floor; that floor does not apply to the original frame compiler.

## 1. Common precision, workspace, and execution ledger

Let $`N=2^n`$, $`n\ge1`$, and supply a Hermitian observable as a real
reflection sum $`H=\sum_r c_rO_r`$, with
$`\Lambda=\sum_r|c_r|\gt0`$. Set coordinate accuracy
$`0\lt\varepsilon_\infty\le\Lambda`$ and failure probability
$`0\lt\delta\lt1`$. Use the common precision parameters

```math
K=\max\!\left\{6,\left\lceil\log_2
                \frac{80\Lambda}{\varepsilon_\infty}\right\rceil\right\},
\qquad P=\max\{n,K\}.
```

An original frame approximation at error $`2^{-K}`$ and an oracle at
that error satisfy the retained accuracy allocation. The state-only
preparation uses error $`2^{-P}\le2^{-K}`$ because its theorem also
requires $`P\ge n`$. Applying that floor to the original frame precision
would overprice the original protocol. Conversely, pricing the state
program at K when $`K\lt n`$ would omit one of its hypotheses.

For a common sufficient shot allocation, take

```math
S=\left\lceil
\frac{800\Lambda^2[1+\ln(2(n+1)/\delta)]}
     {\varepsilon_\infty^2}
\right\rceil.
```

Use S executions of the magnitude stream and, for a complex chart, S
executions of its direct phase stream. This conservative common choice
covers both retained concentration bounds. Classical coefficient and
decoder rounding must still satisfy their stated error budgets. The
original protocol can use smaller sufficient constants; sharing S here
compares quantum construction bounds at one common guarantee rather than
claiming equal optimal sample counts. No derivative of a synthesis word
or observable implementation is taken.

Fix two initialized **compiler** flags and $`b\ge P+n+7`$ arbitrary
dirty qubits. Each protocol additionally reserves the n initialized
logical system qubits, its interference branch, and the controlled
observable's own work. These resources are the same physical allocation
for both choices. Dirty inputs may retain correlations with external
references or previous executions; the complete error and concentration
contracts already include this possibility.

An original compiler may leave flags unused or borrow a flag in its
known-zero state: a guarantee for arbitrary borrowed input includes that
special case. In particular, the retained one-clean complex compiler needs
$`K+n+8`$ borrowed bits. Its effective pool here is $`b+1`$, after
borrowing the spare compiler flag, and

```math
b+1\ge P+n+8\ge K+n+8.
```

This reallocation introduces no third initialized compiler flag. An
original zero-clean construction may similarly borrow both available
flags, although the corollary below already fits using only the b dirty
qubits. Additional banked constructions retain their own literal workspace
thresholds; they are not implied merely by writing $`b=\Theta(P+n)`$.

Let $`\overline t_O=\sum_r |c_r|t_{O_r}/\Lambda`$ be the expected
T-count of the same phase-calibrated controlled observable in either
protocol. Write $`t_{F,\mathbb R}(K)`$ and $`t_{F,\mathbb C}(K)`$
for the selected original frame words; the complex choice may be the
literal prescribed frame or the consistent gauge-fixed frame proved
below. Let $`t_{\rm pair}(P)`$ and $`t_{\rm state}(P)`$ denote the
applicable coherent reference/target preparation and single-state
preparation, and let $`t_C`$ be the actual coarse word's T-count.
The relevant preparation functions depend on whether the chart is real
or complex. Each preparation count includes its forward coarse word.

| Protocol | Expected T-count for the gradient streams |
|---|---|
| Original real chart | $`S(2t_{F,\mathbb R}(K)+\overline t_O)`$ |
| State-only real chart | $`S(t_{\rm pair}(P)+t_C+\overline t_O)`$ |
| Original complex chart | $`S(3t_{F,\mathbb C}(K)+2\overline t_O)`$ |
| State-only complex chart | $`S(t_{\rm pair}(P)+t_{\rm state}(P)+t_C+2\overline t_O)`$ |

The extra coarse occurrence is its actual inverse in the magnitude
readout. The direct phase stream has no inverse frame. Replacing the
expected oracle count by the largest supplied oracle count gives a
deterministic bound. Clifford counts follow the same composition, with
the branch/system readout gates included separately. Compilation,
classical preprocessing, arithmetic precision, histogram reconstruction,
and writing the gradient remain additional costs. Quantum gate bounds
alone do not establish an end-to-end speedup or optimality of gradient
estimation.

## 2. A gauged complex borrowed-frame baseline

The retained arbitrary-budget borrowed theorem is stated for real frames.
It must not be silently applied to a complex phase-dressed frame. For the
common workspace above, the following narrower corollary provides a valid
complex baseline using the already proved phase-prefix construction.

**Corollary.** Supply real phase representatives $`\varphi_x`$ and a
real Hopf frame $`W_{\mathbb R}`$. Define

```math
\mu=\frac1N\sum_x\varphi_x,\qquad
D_0=e^{-i\mu}\operatorname{diag}(e^{i\varphi_x}),\qquad
W'=D_0W_{\mathbb R}.
```

For $`0\lt\eta\le1/64`$, set
$`K=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. If $`b\ge2n`$,
there is an actual logical native unitary $`F_\eta`$ such that its
Clifford+T implementation is exactly $`F_\eta\otimes I_b`$ and

```math
\|F_\eta-W'\|\le\eta,\qquad
T(F_\eta)=O\!\left(\frac{NK}{n+b}+K\sqrt N\right),\qquad
G(F_\eta)=O(NK).
```

No initialized compiler work is used. Every borrowed wire returns
exactly on all logical and dirty inputs, including arbitrary reference
correlations. At $`\eta=2^{-K}`$ the statement uses the common K from
Section 1; $`b\ge P+n+7`$ implies $`b\ge2n`$.

*Proof.* At depth $`d=0,\ldots,n-1`$, let $`P_d`$ be the exact
prefix-addressed phase layer from
[the complex coarse construction](COMPLEX_COARSE_COMPILER.md#2-exact-prefix-phase-factors).
It has $`S_d=2^d`$ determinant-one $`R_z`$ rows, each repeated over
every logical suffix. Its subtree means telescope to
$`D_0=P_{n-1}\cdots P_0`$. The same construction works at fine
precision; its earlier coarse parameter is not needed here.

Allocate the phase-layer errors by

```math
\epsilon_d=\eta\,2^{d-n-1},\qquad
\sum_{d=0}^{n-1}\epsilon_d\lt\eta/2.
```

The [phase-cancelled row words](COMPLEX_COARSE_COMPILER.md#3-literal-determinant-one-row-words)
are literal determinant-one native unitaries of length

```math
w_d=O(1+\log_2(1/\epsilon_d))=O(K+n-d).
```

Their exact reflection rewrites preserve every row-dependent scalar
phase. In particular, the actual native rows need not remain diagonal.
They are applied in the prescribed layer order, and no synthesis error
is hidden in a change of basis or a dirty-work state.

Use the [borrowed reflection interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
with the optional extra control omitted. At depth d choose the largest
power of two

```math
\lambda_d\le\min\{\sqrt{S_d},b-d\}.
```

It uses $`d-\log_2\lambda_d`$ arbitrary external selectors and
$`\lambda_d`$ arbitrary external bank bits. The conservative reservation
$`d+\lambda_d\le b`$ suffices. Prefix bits are unchanged logical
addresses; the row target is distinct. No initialized bank, synthesis
signal, or extra predicate helper is required. Complete the dirty
load/use/unload echo for each reflection before advancing to the next
reflection. The resulting logical layer $`\widehat P_d`$ is implemented
exactly tensored with identity on the dirty pool, at cost

```math
T_d=O\!\left(w_d(S_d/\lambda_d+\lambda_d+1)\right),
\qquad G_d=O(S_dw_d).
```

Let $`B=b-n+1`$. Rounding the bank count to a power of two changes
these bounds by a constant factor, and

```math
T_d=O\!\left(w_d\left(\frac{2^d}{B}+2^{d/2}+1\right)\right).
```

The depth dependence of $`w_d`$ must be retained. With $`j=n-d`$,

```math
\begin{aligned}
\sum_{d=0}^{n-1}2^d(K+n-d)
 &=N\sum_{j=1}^{n}2^{-j}(K+j)=O(NK),\\
\sum_{d=0}^{n-1}2^{d/2}(K+n-d)
 &=\sqrt N\sum_{j=1}^{n}2^{-j/2}(K+j)=O(K\sqrt N),\\
\sum_{d=0}^{n-1}(K+n-d)
 &=O(nK+n^2)=O(K\sqrt N).
\end{aligned}
```

The last estimate uses $`K\ge6`$ and the uniform elementary bounds
$`n,n^2=O(2^{n/2})`$. Consequently the complete phase circuit costs

```math
T(\widehat P)=O(NK/B+K\sqrt N),\qquad G(\widehat P)=O(NK).
```

This weighted sum removes the unnecessary $`n\sqrt N`$ term obtained
by replacing every word length by $`O(K+n)`$. Since $`b\ge2n`$,
$`B=\Theta(n+b)`$.

Compile $`W_{\mathbb R}`$ at error $`\eta/2`$ with the
[real borrowed theorem](../../docs/BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements),
using the same pool. Its precision is $`K+O(1)`$, and its actual logical
word $`\widetilde W_{\mathbb R}`$ has the stated count and exact-return
properties. Set

```math
F_\eta=\widehat P_{n-1}\cdots\widehat P_0
          \widetilde W_{\mathbb R}.
```

Unitary telescoping gives
$`\|F_\eta-D_0W_{\mathbb R}\|\le\eta/2+\sum_d\epsilon_d\lt\eta`$.
Exact dirty return holds after every completed subroutine, so composition
and the actual inverse preserve it on the full input space. The counts
add to the claimed bounds. The depth-zero row, and $`n=1`$, can also be
emitted directly; the same estimates cover them. ∎

The original gradient decoder remains valid with this frame. Magnitude
derivatives obey $`d'_j=e^{-i\mu}d_j=a_jW'|\lambda(j)\rangle`$,
so the usual incoming-amplitude and marker scores are unchanged. The
common scalar cancels from $`W'^\dagger O W'`$. In the direct phase
stream, replacing psi by $`e^{-i\mu}\psi`$ leaves
$`2\operatorname{Im}(\overline{\psi_x}(O\psi)_x)`$ unchanged.
Equivalently, differentiating the selected gauge adds a common-phase
direction whose Hermitian-energy derivative vanishes. Thus both streams
estimate the original physical raw gradients, with the same complete
perturbation ledger and no assumption on nonzero leaf amplitudes.

This corollary approximates the gauge-fixed frame W', not the prescribed
literal $`D_\varphi W_{\mathbb R}`$. It does not assert synthesis of
the removed common phase, extend the real-only theorem to arbitrary
complex unitaries, or resolve the fine complete-frame endpoint. Phase
representatives, certified classical evaluation, native word synthesis,
and the coefficients needed for decoding remain separately charged.

## 3. The available bounds at a common workspace

Write $`\ell=1+\log_2^*(n+2)`$. The following expressions are scales
of **constructive upper bounds**, with independent hidden constants. They
are not measured gate counts, lower bounds, or formulas for optimal costs.
Taking a minimum means choosing a construction whose hypotheses hold.

For either real QBP or consistently gauged complex QBP, the original
protocol can use

```math
F_{\rm borrow}=\frac{NK}{n+b}+K\sqrt N,\qquad
F_{\rm group}=N+K\ell.
```

The borrowed choice follows from the retained real theorem or Section 2;
the grouped choice uses the spare-flag allocation in Section 1. When its
additional bank reservation holds, it also has

```math
F_{\rm bank}=\sqrt{NK}+K\ell+\frac{NK}{b}.
```

For the old literal complex banked compiler, the exact condition at two
available flags is $`b+1\ge2(K+n+8)`$. Use the slightly stronger common
condition $`b\ge2(P+n+8)`$ whenever comparing both banked choices below.
At $`P=K`$, the new bank minimum $`2(P+n+7)`$ alone is one wire short
of the old condition even after borrowing the spare flag. This constant
difference does not change asymptotic regimes, but cannot be omitted from
a literal allocation.

The state-only route always has the scale

```math
A_{\rm state}=N+P.
```

The [banked state corollary](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation)
now supplies, at $`b\ge2(P+n+7)`$,

```math
B_{\rm state}=\sqrt{NP}+P+\frac{NP}{b}+n\sqrt N.
```

It replaces the fine residual table's mask queries by exact whole-word
dirty-bank queries. The additional banks are disjoint from the occupied
core and the borrowed synthesis signal. Exact query equality preserves
the existing full-operator error proof and needs no new clean flag.
The last term is the charged exact-return coarse circuit; an approximately
returned synthesis core cannot be substituted for that logical circuit.
The extra coarse inverse in magnitude readout changes only a constant
factor in these expressions.

Thus the quantum T budgets from Section 1 have order at most

```math
S\bigl(\min F+\overline t_O(K)\bigr),\qquad
S\bigl(\min\{A_{\rm state},B_{\rm state}\}
                         +\overline t_O(K)\bigr),
```

where unavailable choices are omitted. Fixed factors for the two complex
streams remain explicit in Section 1. In particular the controlled
observable is compiled at K, **not P**. Oracle calls are exactly S for a
real chart and $`2S`$ for a complex chart in either protocol.

The corresponding Clifford bounds are

```math
G_{\rm original}=O\bigl(S(NK+\overline g_O(K)+n)\bigr),\qquad
G_{\rm state}=O\bigl(S(NP+\overline g_O(K)+n)\bigr).
```

There is no improved Clifford order here. If $`K\lt n`$, the state's
retained bound is larger. With sufficiently many initialized compiler
qubits, the original frame theorem already supplies
$`O(\sqrt{NK}+K+NK/(n+a+b))`$ T gates per frame. That construction
retains its sufficient-clean hypothesis; it is not available merely
because the same number of arbitrary dirty qubits exists.

## 4. Where the bound changes

All statements in this section compare the expressions above along
$`n\to\infty`$. Equal order does not locate a finite numerical crossover.
Even a diverging ratio improves the available upper bound, rather than
proving a separation between optimal gradient algorithms.

At a minimum-order pool $`b=\Theta(P+n)`$, with the stated literal
eligibility enforced, the comparison is:

| Accuracy-bit regime | Original bound expression | State bound expression |
|---|---|---|
| $`K=o(n)`$ | $`NK/n`$ | $`N`$ |
| $`K=\Omega(n)`$, $`K\ell=O(N)`$ | $`N`$ | $`N`$ |
| $`K\ell\gg N`$ | $`K\ell`$ | $`N+K`$ |

Indeed the borrowed term is $`\Theta(NK/n)`$ in the first row.
For $`K\ge n`$, its $`NK/(n+b)`$ term is of order N. If
$`K\ell\gg N`$, both its other term and every grouped/banked frame
expression are at least the scale $`K\ell`$, while grouping attains
that scale. Increasing a minimum-order pool by only a fixed factor does
not change these conclusions. Fixed-accuracy QBP belongs to the first
row; it does not require the high-precision endpoint $`K=N`$.

Additional banks produce a sharper comparison. Assume
$`K\ge n^2`$ and $`b\ge2(P+n+8)`$, so $`P=K`$. The coarse term
$`n\sqrt N`$ is then absorbed by $`\sqrt{NK}`$. The minima of the
available expressions are, up to constant factors,

```math
Q=\sqrt{NK}+K+\frac{NK}{b},\qquad
R=\sqrt{NK}+K\ell+\frac{NK}{b}
```

for state and original routes respectively. To see that no omitted
choice changes their order, use $`b\ge2K`$ and split at $`K=N`$:
Q never exceeds a constant times $`N+K`$, and R never exceeds a constant
times $`N+K\ell`$. The borrowed expression also bounds R below up to
a constant, since $`K\sqrt N`$ dominates both $`\sqrt{NK}`$ and
$`K\ell`$, and $`n+b=\Theta(b)`$.

The exact ratio of these two displayed scales is

```math
\frac RQ
=1+\frac{\ell-1}{\sqrt{N/K}+1+N/b}.
```

Consequently a diverging ratio occurs precisely when

```math
K\gg N/\ell^2,\qquad b\gg N/\ell,
```

along $`\ell\to\infty`$. These regimes lie beyond both sufficient
matching conditions retained for the one-clean banked frame
bound. The result improves a task-specific construction there; it does
not resolve that frame-synthesis gap.

For illustration, if $`b=\Omega(\sqrt{NK})`$ in addition to the common
bank reservation, then:

| Precision regime | Original expression | State expression |
|---|---|---|
| $`K=o(n^2)`$ | $`\sqrt{NK}`$ | $`n\sqrt N`$ |
| $`n^2\le K\le N/\ell^2`$ | $`\sqrt{NK}`$ | $`\sqrt{NK}`$ |
| $`N/\ell^2\ll K\ll N`$ | $`K\ell`$ | $`\sqrt{NK}`$ |
| $`K=\Theta(N)`$ | $`N\ell`$ | $`N`$ |
| $`K\gg N`$ | $`K\ell`$ | $`K`$ |

The first row uses $`P=\max(n,K)`$ and retains the coarse cost; it does
not extend the simplification $`P\ge n^2`$ below its domain. For example,
at fixed accuracy and $`b=\Theta(\sqrt N)`$, the original bound is
$`O(\sqrt N)`$, whereas this banked state construction gives
$`O(n\sqrt N)`$. In the high-precision range,
$`K=N/\ell^{3/2}`$ and a sufficiently large bank give a ratio of order
$`\ell^{1/4}`$ between the displayed expressions.

Oracle and sampling costs remain decisive. If $`\overline t_O=O(Q)`$,
a diverging R/Q ratio survives in the complete quantum T budgets.
If the oracle term is $`\Omega(R)`$, those budgets have the same order.
For supplied Pauli-string observables, their controlled implementations
can be Clifford circuits, so an added oracle T charge need not dominate.
No oracle-cost assumption is silently built into the comparison.

Moreover K measures the requested gradient accuracy. Under Section 1's
choice, $`\Lambda/\varepsilon_\infty=\Theta(2^K)`$, and hence

```math
S=\Theta\!\left(4^K[1+\ln((n+1)/\delta)]\right).
```

Removing a repeated precision term does not reduce this sampling cost.
The regime $`K=N`$ therefore describes extremely fine gradient accuracy,
not a practical fixed-accuracy sampling advantage. The bounded native
examples likewise establish implementation correctness; their literal
counts favor the original protocol on those exact targets.

## 5. Classical work and preprocessing

The original decoder already admits an integer histogram implementation.
Let A be its signed magnitude histogram, including observable-term signs,
and let $`F_N`$ be the unnormalized Walsh matrix. Then

```math
k=F_NA,\qquad \widehat g_j=\frac{2\Lambda}{S}a_j k_{\lambda(j)}.
```

Integer butterflies cost $`O(Nn)`$ arithmetic, and evaluating incoming
amplitudes and producing the output costs $`O(N)`$. With its record-wise
alternative, the retained original bound is
$`O(S+N\min\{S,n\})`$. Its separate phase histogram adds $`O(S+N)`$.
The corrected state decoder costs $`O(S+Nn)`$, including actual-C
application and its reverse tree traversal. Both use $`O(N)`$ stored
values after their coefficients exist. There is no improved classical
decoding order, and the original can be cheaper when $`S\lt n`$.

Counters need $`O(1+\log S)`$ exact bits. Original incoming weights
need $`O(K+\log(n+1))`$ fractional bits; corrected coefficients and
contractions use $`O(K+n+\log(n+1))`$. Arithmetic-operation counts must
therefore be multiplied by the costs of integer addition, certified
coefficient evaluation, and the applicable precision-dependent arithmetic.
Reading shot labels and writing all $`\Theta(N)`$ output coordinates
are also charged.

State compiler preprocessing first evaluates the amplitudes and the actual
native coarse coefficients, applies $`C^\dagger`$ to the state in
$`O(Nn)`$ scalar operations for a complex chart, and constructs the
residual SU(2) tables. A conservative sufficient coefficient precision is

```math
B\ge2(P+10)+\frac32 n+\log_2(n+1)+O(1).
```

Indeed rounded stable coarse updates have error bounded by
$`O(Nn\,2^{-B})`$; multiplying the tail by $`\sqrt N`$ gives
$`O(N^{3/2}n\,2^{-B})`$. This meets the retained table accuracy
$`2^{-2q-12}`$, with $`q=P+10`$. The completion square root obeys
the uniform bound $`e+\sqrt{2e}`$, so proximity of the head coefficient
to the unit circle requires no inverse-gap precision assumption. This
bounds the required digits, not the runtime of arbitrary supplied
computable evaluators or a native-word synthesis search.

The [bounded-input source-table result](../../docs/OPERATOR_SOURCE_COMPILER.md#10-classical-table-construction)
prices a particular table-generation procedure. It does not automatically
price all coarse-word synthesis, grouped compilation, or certified Euler
search for every fastest construction in this comparison. Keep separate
compiler-preprocessing costs $`\mathcal P_{\rm old}`$ and
$`\mathcal P_{\rm state}`$, together with the shared classical sampler
setup and term-selection costs. A complete resource account comprises
those preprocessing costs, the T and Clifford budgets in Section 3,
classical record processing, and output bit costs. They cannot be collapsed
into a single runtime without a machine and input model.

Preprocessing may be reused across shots and observable batches at the
same parameter tuple. A parameter update generally changes the state,
coarse word, residual tables, and corrected decoder; it requires the new
applicable preprocessing. Large supplied phase windings and certified
range reduction are not free. Neither a free amplitude oracle nor QRAM
is assumed.

The subsequent [bounded-input audit](BOUNDED_INPUT_QBP.md) makes these
costs constructive for explicit dyadic angles/phases and Pauli lists. Its
[algebraic residual procedure](RESIDUAL_TABLE_PREPROCESSING.md) removes
the Euler search for this particular U(z). Grouped/state construction then
uses only coarse native-word searches, giving a conservative polynomial
preprocessing bound with the stated small-register qualifications. A
direct fine borrowed-frame search retains its different runtime. These
are bounded-input results; arbitrary effective evaluators still have no
uniform running-time guarantee.

## 6. Research conclusion and the next boundary

The task-specific route removes the repeated grouped precision term for
state preparation and keeps the original logarithmic depth-block shot
order. Additional dirty banks strengthen its T bound, with two compiler
flags, in the high-precision regimes identified above. Fixed-accuracy
gradients retain a better original frame construction. Clifford order,
oracle query count, and classical decoding order do not improve here.

The bounded-input question is now answered in its
[own chapter](BOUNDED_INPUT_QBP.md). Explicit Pauli access admits a
deterministic division-free classical gradient calculation, as well as
term-only classical sampling. Under uniform elementary-gate costs, the
retained Clifford upper expressions also remove the quantum T separation.
The result remains a quantum T-resource improvement in stated regimes,
not an end-to-end gradient speedup. A new advantage claim needs its own
observable-access and execution model; another special native fixture
alone would not supply it. The fine complete-frame endpoint remains a
separate problem.

The [state-based depth composition](STATE_QBP_DEPTH.md) adds a separate
workspace-versus-T-depth result. Section 7 compares its complete execution
ledger with eligible original schedules, without identifying T-depth with
total circuit depth or elapsed runtime.

## 7. State-based T-depth comparison

Keep Section 1's K, $`P=\max(n,K)`$, S, physical allocation, and observable
access. Put $`B_0=P+n+7`$ and $`\ell=1+\log_2^*(n+2)`$.
For a chosen literal circuit let d denote its T-depth, and write
$`\overline d_O=\sum_r|c_r|d_{O_r}(K)/\Lambda`$. Each observable uses
precision K, not P. Serial composition gives the following expected
T-depth upper schedules; replacing the oracle average by its maximum
gives deterministic upper schedules.

| Protocol | Serial gradient T-depth upper schedule |
|---|---|
| Original real chart | $`S(2d_{F,\mathbb R}(K)+\overline d_O)`$ |
| State-based real chart | $`S(d_{\rm pair}(P)+d_C+\overline d_O)`$ |
| Original complex chart | $`S(3d_{F,\mathbb C}(K)+2\overline d_O)`$ |
| State-based complex chart | $`S(d_{\rm pair}(P)+d_{\rm single}(P)+d_C+2\overline d_O)`$ |

The pair and single-state depths include their forward coarse word.
Magnitude readout adds its actual inverse, of depth $`d_C`$; the direct
phase stream does not. Branch and system readout are Clifford and add
zero T-depth, but their elementary depth is not zero. Both choices retain
the same S magnitude executions and, for complex charts, S phase
executions. Parallel shots require additional systems, branches, compiler
flags, and work; the ledger does not grant them for free.

The [new composition theorem](STATE_QBP_DEPTH.md) supplies two schedules
for the state-based compiler work, including magnitude readout:

| Schedule and dirty reservation | T-depth | T-count | Clifford count |
|---|---|---|---|
| A: $`b\ge2B_0`$ | $`O(NP/b+P+n^3)`$ | $`O(NP)`$ | $`O(NP)`$ |
| B: $`b\ge16(B_0+\sqrt{NP})`$ | $`O(P+n^3)`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ | $`O(NP)`$ |

These bounds apply to real and gauge-fixed complex charts. The observable
is added by the preceding ledger. Each row's count and depth hold for the
same circuit. In particular, filling the pool with banks in schedule A
does not retain schedule B's smaller T-count.

The original **real-frame** schedules remain eligible competitors at
their own precision K. The [depth-optimized construction](../../docs/T_DEPTH_COMPILER.md)
gives, when $`b\ge2(K+n+7)`$,

```math
D_{A,\rm old}=O\!\left(\frac{NK}{b}
 +\min\{nK+n^2,\ K\ell+n^3\}\right),\qquad T,G=O(NK).
```

The [count-preserving construction](../../docs/PARALLEL_DIRTY_LOOKUP.md) has a fixed
sufficient constant $`c_{\rm old}`$ and, at
$`b\ge c_{\rm old}(K+n+7+\sqrt{NK})`$, gives

```math
D_{B,\rm old}=O\!\left(\min\{nK+n^2,\ K\ell+n^3\}\right),
\qquad T=O(\sqrt{NK}+K\ell),\qquad G=O(NK).
```

The [amortized real-frame schedule](../../docs/AMORTIZED_DIRTY_LOOKUP.md) adds another
eligible same-circuit choice, at $`b\ge17(K+n+7)`$:

```math
D_{\rm amortized,old}=O\!\left(\frac{NK}{b^2}+nK+n^2\right),
\qquad T=O\!\left(\sqrt{NK}+\frac{NK}{b}+nK\right),\qquad G=O(NK).
```

Use K here, not the state preparation floor P. In its nonempty matching
range $`b\le\sqrt{NK/(nK+n^2)}`$, frame count and depth are both
optimal in order. In particular, $`K=\Theta(n)`$ and sufficient
$`b=\Theta(n)`$ give per-frame $`T=\Theta(N)`$ and
$`D_T=\Theta(N/n)`$. Substitution into the first row's serial gradient
ledger still charges both frame calls, the observable, and all S shots.
This is a frame-compiler theorem, not a lower bound for gradient estimation
or a change to the state-based schedules. The A/B comparison below remains
a comparison of its stated expressions; all eligible schedules may be used.

The subsequent [uniform real-frame theorem](../../docs/UNIFORM_PRECISION_DEPTH.md)
is also eligible at $`b\ge17(K+n+7)`$. It retains the last
same-circuit T/Clifford bounds and improves depth to
$`O(NK/b^2+nK)`$. For $`6\le K\le\log_2(n+2)/16`$, it gives
$`T=O(\sqrt{NK}+NK/b)`$ and $`D_T=O(NK/b^2+n)`$.
The fixed-accuracy predecessor already gives $`O_K(n)`$ depth at
sufficient square-root dirty width. Thus the historical original-B
expression compared below is not the best current fixed-accuracy frame
bound. These improvements still use K, charge both frame calls and all
shots, and assert no optimal gradient-estimation or complex-frame depth.

The common A pool satisfies the old A threshold. For a common B comparison,
require literally

```math
c_* =\max\{16,c_{\rm old}\},\qquad
b\ge c_*(B_0+\sqrt{NP}).
```

Since $`P\ge K`$, this allocation supplies both B schedules. Every eligible
original count construction also gives $`D_T\le T`$, so keep the smaller
such upper expression if useful. For complex original frames this includes
the gauge-fixed borrowed baseline in Section 2 and the eligible literal
grouped/banked constructions in Section 3. The real depth-composition
theorems above are not themselves complex-frame depth theorems. In
particular, at $`P=K,b=2(P+n+7)`$, the old literal complex bank reservation
is still one wire short after borrowing the spare clean flag; it requires
$`b+1\ge2(K+n+8)`$. Its unbanked and gauged borrowed choices remain valid.

For a comparison of the two count-preserving **real** depth expressions,
define

```math
Q_D=P+n^3,\qquad R_D=\min\{nK+n^2,\ K\ell+n^3\}.
```

Along $`n\to\infty`$, a diverging ratio $`R_D/Q_D`$ occurs precisely
when $`K\ell\gg n^3`$. In that regime $`P=K`$ and the grouped old
expression is selected, so

```math
\frac{R_D}{Q_D}=\frac{K\ell+n^3}{K+n^3}
 =1+\frac{\ell-1}{1+n^3/K}.
```

If $`K\ell=O(n^3)`$, $`R_D\le K\ell+n^3=O(n^3)\le O(Q_D)`$,
so no divergent gain is implied. At fixed accuracy the original B
expression is $`\Theta(n^2)`$, while the new B expression is
$`\Theta(n^3)`$. These are comparisons of displayed upper expressions,
not lower bounds, optimal schedules, or a minimum over all new choices:
with additional banks schedule A can improve depth at its separate
T-count cost.

For example, take $`K=n^3`$ and the common B pool. Then $`P=K`$,
the new depth expression is $`\Theta(n^3)`$, and the original one is
$`\Theta(n^3\ell)`$. Both simultaneous T-count expressions have
leading order $`\sqrt{NK}=n^{3/2}\sqrt N`$, with Clifford bound
$`O(Nn^3)`$. Thus the available depth expression can improve even when
the T-count orders coincide. Independent hidden constants prevent this
asymptotic comparison from locating a finite crossover.

In the improving regime, the ratio survives the serial oracle ledger if
$`\overline d_O=O(Q_D)`$; if $`\overline d_O=\Omega(R_D)`$, those total
depth expressions have the same order. S still has the exponential-in-K
dependence in Section 4.
Clifford depth, measurement and classical latency, preprocessing, and
gradient reconstruction remain additional costs. Arbitrary Clifford
circuits between T layers can have substantial depth. This comparison
therefore establishes neither total-runtime advantage nor optimal
gradient-query complexity, and it does not remove the explicit Pauli
classical baselines in the bounded-input audit.
