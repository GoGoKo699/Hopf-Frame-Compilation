# Raw Hopf gradients in a coarse-frame measurement basis

[State-only compiler](STATE_ONLY_COMPILER.md) · [Reference-state decoder](REFERENCE_STATE_QBP.md) · [Original approximation contract](QBP_APPROXIMATION.md)

A coarse approximation to the complete Hopf frame supplies a measurement
basis with uniformly bounded raw-gradient records. Fine precision is needed
only for the prepared state and the controlled observable. The decoder uses
the actual coarse circuit, including its complex phases, and has a classical
histogram implementation with $`O(S+Nn)`$ arithmetic operations. It requires
neither a fine inverse frame nor a dense table of all derivative vectors.

For every real Hopf angle tuple and every allowed Hermitian-unitary observable,
the ideal means are the original raw magnitude derivatives. Their sufficient
sample order is $`O(\Lambda^2\varepsilon_\infty^{-2}[1+\log(n/\delta)])`$,
without the derivative-envelope factor in the general
[reference-state construction](REFERENCE_STATE_QBP.md#3-a-reference-for-every-real-angle-tuple).
With the stated two-clean state compiler, one execution costs
$`O(N+L')+t_O`$ T gates. This is a different decoder and measurement
distribution, not a fine compilation of the prescribed full frame.

The [complex extension](COMPLEX_COARSE_QBP.md) supplies the same magnitude
bound for phase-dressed states and adds the direct leaf-phase stream. Its
[compiler proof](COMPLEX_COARSE_COMPILER.md) fixes a consistent common
phase and records actual native prefix rows with exact dirty return.

The interference and change-of-basis identities use the standard ingredients
attributed in [the reference-state decoder](REFERENCE_STATE_QBP.md).
The statements below price their Hopf-specific combination. No priority claim
is made for Hadamard interference or classical reverse differentiation.

## 1. The actual coarse circuit and transformed derivatives

Fix $`N=2^n`$, $`n\ge1`$, the real complete frame $`W`$, its state
$`|\psi\rangle=W|0^n\rangle`$, and its raw derivatives

```math
|d_j\rangle=\partial_{\theta_j}|\psi\rangle
           =a_jW|\lambda(j)\rangle.
```

Here $`a_j`$ is the real oriented incoming amplitude, including its sign.
Choose the actual recorded native coarse word C with exact dirty return and

```math
\|C-W\|\le\delta_c,
\qquad \delta_c\le\min\{1/64,1/(4\sqrt N)\}.
```

The [coarse construction](STATE_ONLY_COMPILER.md#2-a-native-coarse-circuit-makes-the-tail-small)
gives $`T(C)=O(N)`$ and $`G(C)=O(Nn)`$ under the reservation below.
It is a logical product of recorded two-mode native words. C can be complex;
it is never replaced by an ideal real surrogate in this protocol.

Write $`H_n=H^{\otimes n}`$ and define

```math
|r\rangle=C|0^n\rangle,\qquad U=H_nC^\dagger,\qquad
|f_j\rangle=U|d_j\rangle.
```

The exact reference after U is uniform:
$`U|r\rangle=H_n|0^n\rangle=N^{-1/2}\sum_x|x\rangle`$.
The actual inverse of C is used, so there is no residual error in this
identity. Coarse error will bound record sizes; it is not added to the
gradient bias.

## 2. Universal means from randomized X/Y interference

Prepare the coherent joint state

```math
|\Omega\rangle
=\frac{|0\rangle_BC|0^n\rangle+|1\rangle_B|\psi\rangle}{\sqrt2},
```

apply a phase-calibrated controlled Hermitian-unitary O on branch one, and
then apply the same U to both branches. Choose X or Y measurement on B with
equal probability, independently of the previous data, and measure the
system in the computational basis. Denote the branch eigenvalue by
$`s=(-1)^b`$, the chosen basis by $`q\in\{X,Y\}`$, and the leaf by x.
For $`v=UO\psi`$, the signed probabilities conditioned on the chosen basis
obey

```math
\sum_s s\Pr(s,x\mid X)=\frac{\operatorname{Re}v_x}{\sqrt N},\qquad
\sum_s s\Pr(s,x\mid Y)=\frac{\operatorname{Im}v_x}{\sqrt N}.
```

Use the records

```math
Y_j(q,s,x)=4s\sqrt N\begin{cases}
\operatorname{Re}f_{jx},&q=X,\\
\operatorname{Im}f_{jx},&q=Y.
\end{cases}
```

The factor four includes the probability one half of each basis choice.
Consequently

```math
\begin{aligned}
\mathbb E Y_j
&=2\operatorname{Re}\sum_x\overline{f_{jx}}v_x\\
&=2\operatorname{Re}\langle d_j|O|\psi\rangle
=\partial_{\theta_j}\langle\psi|O|\psi\rangle.
\end{aligned}
```

The observable may be complex. The target state and its parameter derivatives
remain real, but f need not be real because the native C may be complex.
Omitting the Y records or replacing f by ideal Walsh coefficients generally
changes the means. No derivative of C, its native gate word, or the chosen
measurement basis is taken: the tuple and C are fixed throughout the protocol.
The same tuple, coarse word, and decoder are retained throughout each batch.
This is a real magnitude-gradient result; it makes no additional claim about
the separate leaf-phase stream.

For $`H=\sum_t c_tO_t`$, choose t independently with probability
$`|c_t|/\Lambda`$, where $`\Lambda=\sum_t|c_t|\gt0`$, and multiply each
record by $`\Lambda\operatorname{sgn}(c_t)`$. Coefficients and sampling
probabilities have the [existing exact classical contract](QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights).
The means then give the raw derivatives of $`\langle\psi|H|\psi\rangle`$.
Classically preparing and sampling this distribution remains charged.

## 3. Every depth record has norm at most five

Put $`E=C^\dagger W-I`$, so $`\|E\|\le\delta_c`$. At depth d, let
$`J_d`$ inject its marker labels and let $`A_d=\operatorname{diag}(a_j)`$
on those nodes. The matrix of transformed derivative columns at this depth is

```math
F_d=H_n(I+E)J_dA_d.
```

The squared incoming amplitudes sum to one at each depth. Thus every row of
$`H_nJ_dA_d`$ has Euclidean norm $`1/\sqrt N`$. Every row of the remaining
term has norm at most
$`\|E\|\|A_d\|\le\delta_c`$. Therefore, for every x,

```math
\|(f_{jx})_{j\text{ at depth }d}\|_2
\le\frac1{\sqrt N}+\delta_c
\le\frac5{4\sqrt N}.
```

Taking real or imaginary parts cannot increase the norm. Each sampled depth
record consequently satisfies

```math
\boxed{\|Y^{(d)}\|_2\le5\Lambda.}
```

No inverse incoming amplitude or metric coefficient appears. This bound
holds at singular tuples and for unrestricted real angles. It is a bound for
each depth and hence for raw-coordinate accuracy after a union bound, not a
dimension-free bound on the Euclidean norm of the concatenated gradient.

## 4. Preparation error, finite decoding, and reused dirty work

Require complete joint-state preparation error at most $`\eta_P`$, including
relative branch phase, both compiler flags, arbitrary dirty input, and any
external reference. Each controlled observable has full initialized-isometry
error at most $`\eta_O`$ on arbitrary branch-system inputs and retained work.
All work remains coherent until the final measurement. The exact native
$`C^\dagger`$ and H gates do not increase the preceding state error, which
is at most $`\eta_P+\eta_O`$.

For any real unit vector in one depth's coordinate space, its scalar decoded
observable has norm at most $`5\Lambda`$. The bounded-observable inequality
and Euclidean duality therefore give depth-block bias at most

```math
10\Lambda(\eta_P+\eta_O).
```

Define normalized decoder coefficients
$`\kappa^X_{jx}=\sqrt N\operatorname{Re}f_{jx}`$ and
$`\kappa^Y_{jx}=\sqrt N\operatorname{Im}f_{jx}`$.
Uniform coefficient error at most $`\tau`$ changes every empirical
coordinate average by at most $`4\Lambda\tau`$. Alternatively, the
histogram algorithm below may certify its final deterministic rounding error
directly, without constructing these coefficients.

For $`0\lt\varepsilon_\infty\le\Lambda`$ and $`0\lt\delta\lt1`$, choose

```math
\eta_P,\eta_O\le\frac{\varepsilon_\infty}{80\Lambda},\qquad
\tau\le\frac{\varepsilon_\infty}{16\Lambda},\qquad
S=\left\lceil\frac{800\Lambda^2[1+\ln(2n/\delta)]}{\varepsilon_\infty^2}\right\rceil.
```

These allocate at most $`\varepsilon_\infty/4`$ to bias and
$`\varepsilon_\infty/4`$ to decoding, leaving
$`\varepsilon_\infty/2`$ for sampling. If the dirty bank is reused,
conditional centered depth records have norm at most $`10\Lambda`$.
The [same Hilbert-space martingale bound](QBP_APPROXIMATION.md#101-reusing-a-dirty-bank-across-executions)
as in the original protocol gives

```math
\Pr\!\left(\left\|\frac1S\sum_{t=1}^S
  (Y_t^{(d)}-\mathbb E[Y_t^{(d)}\mid\mathcal F_{t-1}])
  \right\|_2\ge\varepsilon_{\rm stat}\right)
\le2\exp\!\left(-\frac{S\varepsilon_{\rm stat}^2}{200\Lambda^2}\right).
```

Every conditional bias obeys the same uniform estimate. A union bound over
the n depths proves raw-gradient $`\ell_\infty`$ error at most
$`\varepsilon_\infty`$ with probability at least $`1-\delta`$.
The system and all declared clean inputs are initialized between completed
executions. No intermediate reset or independence assumption about a reused
dirty bank is needed. The independently sampled observable label and X/Y
basis retain their specified distributions.
Both choices are drawn freshly and independently of the previous history in
each execution, and independently of each other.

## 5. A histogram decoder without a dense derivative table

Keep signed integer histograms

```math
A_x=\sum_{t:q_t=X,\,x_t=x}\operatorname{sgn}(c_{u_t})s_t,\qquad
B_x=\sum_{t:q_t=Y,\,x_t=x}\operatorname{sgn}(c_{u_t})s_t,
```

where $`u_t`$ is the classically selected observable term. Let
$`h=(4\Lambda\sqrt N/S)(A+iB)`$. For the real Hopf Jacobian J whose
columns are $`d_j`$, the exact empirical estimate is

```math
\widehat g
=\operatorname{Re}(J^\dagger U^\dagger h)
=J^{\mathsf T}\operatorname{Re}(C H_n h).
```

First apply $`H_n`$ by a fast Walsh transform. Then apply the recorded
logical two-mode words of C in their actual chronological order, including
their complex entries. This takes $`O(Nn)`$ and $`O(N)`$ scalar operations,
respectively, once those two-by-two coefficients are available. It does not
simulate the Hilbert space of the dirty implementation.

At the leaves initialize $`\beta_x=\operatorname{Re}(C H_n h)_x`$.
Using the original tuple's signed incoming amplitudes $`a_j`$, propagate

```math
\beta_j=c_j\beta_L+s_j\beta_R,\qquad
\widehat g_j=a_j(-s_j\beta_L+c_j\beta_R).
```

This final forward/reverse tree computation takes $`O(N)`$ operations.
It differentiates the original Hopf state with the measured leaf weights
held fixed. In particular, it does not differentiate C or the histograms.
Updating the histograms costs $`O(S)`$, and writing the output costs
$`\Theta(N)`$. The total is $`O(S+Nn)`$ arithmetic operations and
$`O(N)`$ storage, with no $`N\times(N-1)`$ derivative table.

### Exact integer Walsh transform and finite precision

The square-root normalization can be canceled before floating arithmetic.
Let $`F_N=\sqrt N H_n`$, the unnormalized sign-valued Walsh matrix, and set

```math
k=F_N(A+iB),\qquad u=k/S.
```

Compute the Walsh butterflies using exact integers. Since
$`\|A+iB\|_1\le S`$, each intermediate component has magnitude at most S;
its real and imaginary integer counters need $`O(1+\log S)`$ bits. Then
$`|u_x|\le1`$, $`\|u\|_2\le\sqrt N`$, and the same estimate is

```math
\widehat g=4\Lambda J^{\mathsf T}\operatorname{Re}(Cu).
```

Here is a sufficient conservative rounding analysis. Approximate each actual
two-by-two coefficient of C and each sine/cosine by certified absolute error
at most $`2^{-P}`$, clamp the real trigonometric factors to $`[-1,1]`$, and
round each scalar arithmetic operation with absolute error at most
$`2^{-P}`$. Use sufficient integer range. When $`N2^{-P}`$ is bounded by
a small constant, the rounded divisions $`u=k/S`$ contribute at most
$`O(\sqrt N2^{-P})`$ in Euclidean norm. The subsequent $`O(N)`$
near-unitary two-mode updates give the dominating bound

```math
\|\widetilde y-Cu\|_2=O(N^{3/2}2^{-P}).
```

Indeed, accumulated coefficient error is $`O(N2^{-P})`$ in operator norm,
the input norm is at most $`\sqrt N`$, and update rounding contributes
$`O(N2^{-P})`$. A raw derivative has norm at most one. Its trigonometric
products contain at most n factors; the division-free reverse recurrence
has only n levels and row norms $`1+O(2^{-P})`$. A conservative bound on
each final decoded coordinate is therefore

```math
|\widetilde g_j-\widehat g_j|
\le O\!\left(\Lambda[N^{3/2}+Nn]2^{-P}\right).
```

These computations are performed before the final known scaling by
$`4\Lambda`$. Evaluate that rescaling with absolute error
$`O(\Lambda2^{-P})`$; this preserves the displayed bound also when
$`\Lambda\lt1`$. The description length of the supplied coefficients,
the representation of that scale, and final output rounding remain charged.

Thus $`P=O(n+\log(n+1)+\log(\Lambda/\varepsilon_\infty))`$ fractional
bits, with suitable fixed constants and guard bits, suffice for deterministic
decoding error at most $`\varepsilon_\infty/4`$. This replaces the
coefficientwise $`\tau`$ condition above. Angle descriptions, certified
trigonometric evaluation, evaluation of the actual native coarse coefficients,
and arithmetic bit costs remain charged. In particular, these decoder bounds
do not establish a bit-complexity theorem for the state compiler's separate
fine-program preprocessing.

## 6. Native preparation, workspace, and scope

Set

```math
\ell=\max\{6,\lceil\log_2(80\Lambda/\varepsilon_\infty)\rceil\},
\qquad L'=\max\{n,\ell\}.
```

The [common-coarse state preparation](STATE_ONLY_COMPILER.md#7-a-complex-native-reference-from-the-same-coarse-word)
prepares the required branch pair coherently with two compiler clean flags,
$`b_F\ge L'+n+7`$ arbitrary dirty qubits, and
$`T=O(N+L')`$, $`G=O(NL')`$. Its reference is the actual complex
$`C|0^n\rangle`$. It uses the same C in both branches and prepares the
appropriate residual states before that shared C; this is not an invocation
of a two-real-state theorem on a potentially complex reference.

The protocol additionally applies the exact inverse C and n Hadamards for
measurement. There are therefore two charged coarse appearances, C in
preparation and $`C^\dagger`$ in readout. Their combined costs remain
$`O(N)`$ T gates and $`O(Nn)`$ Cliffords. There is no additive
$`\delta_c`$ bias: the reference, basis, and classical coefficients all use
the same literal C. No QRAM, free state oracle, reset, or supplied catalyst
is assumed.

Reserve the initialized interference branch and the observable's own work
in addition to the two compiler flags. The branch is preserved by the
coherent preparation. Dirty reservations use the actual $`L'`$, including
its constants; they are not replaced by an asymptotic $`\Theta(L')`$
allowance. With average selected-oracle costs $`\overline t_O,\overline g_O`$,

```math
\mathbb E\mathcal T_{\nabla,\rm coarse}
=O\!\left(\frac{\Lambda^2[1+\log(n/\delta)]}{\varepsilon_\infty^2}
 [N+L'+\overline t_O]\right),
```

and the expected Clifford count is
$`O(S[NL'+\overline g_O+n])`$. Replacing average oracle costs by their
maxima gives deterministic bounds. Classical preprocessing and the histogram
decoder above are additional end-to-end costs.

This removes the general reference decoder's growing envelope factor while
using the new state-only precision bound. It does not establish optimal
gradient-query complexity, a universal speedup, or a fine full-frame compiler.
At fixed accuracy the original protocol's available per-frame
$`\Theta(\sqrt N)`$ T-count can still be preferable. Reading a dense angle
tuple, compiling its tables, implementing the observable, and writing all
gradient coordinates remain part of any meaningful advantage comparison.

The [banked state refinement](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation)
also applies here: $`b\ge2(L'+n+7)`$ gives per-execution T-count
$`O(\sqrt{NL'}+L'+NL'/b+n\sqrt N)`$ apart from the observable,
including the charged coarse inverse. The
[task-cost comparison](QBP_COST_COMPARISON.md) uses common accuracy,
workspace, and oracle access, and includes the original histogram decoder
and both routes' preprocessing obligations.

## 7. Finite evidence

[Coarse-frame decoder tests](../tests/test_coarse_frame_qbp.py) check randomized
X/Y means, the complete depth-record bound, the histogram contraction, and
complex coarse phases on four and eight modes, including singular angles.
A concrete phase fixture shows that dropping Y or using ideal Walsh scores
can create a nonzero decoded gradient when the true gradient is zero.
These tests are finite algebraic diagnostics, not a native emitter or a
large-instance performance claim.

The subsequent [native integration example](NATIVE_COARSE_QBP.md) emits
the complete two-qubit protocol through the state compiler's exact
finite-size fallback. It checks an actively used arbitrary dirty helper,
literal branch phase, complex controlled observables, actual adjoints,
and corrected readout. Its finite gate counts are compared with the
original protocol on the same exact target; it does not emit the general
fine-precision residual table.
