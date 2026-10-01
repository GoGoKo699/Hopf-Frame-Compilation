# State-based Hopf QBP: consolidated theorem

[Real protocol](COARSE_FRAME_QBP.md) · [Complex protocol](COMPLEX_COARSE_QBP.md) · [Task-cost comparison](QBP_COST_COMPARISON.md) · [Bounded-input construction](BOUNDED_INPUT_QBP.md)

This note collects the proved state-preparation, decoding, and resource
statements into one corollary. The linked chapters contain the proofs;
no stronger compiler or optimality assertion is introduced here.

## 1. Input, access, and accuracy

Let $`N=2^n`$, $`n\ge1`$, and supply the $`N-1`$ real Hopf angles
theta defining $`\psi_{\mathbb R}=W_{\mathbb R}(\theta)|0^n\rangle`$.
For a complex chart, also supply N real phase representatives phi and set

```math
\psi=D_\varphi\psi_{\mathbb R},\qquad
D_\varphi=\operatorname{diag}(e^{i\varphi_x}).
```

For a real chart take $`D_\varphi=I`$ and request only magnitude
derivatives. Angles are unrestricted; singular tuples and zero amplitudes
are included. The data admit certified classical evaluation. Their input
length, evaluation time, and any large-phase range reduction are charged.

Supply a fixed Hermitian observable
$`H=\sum_{t=1}^M c_tO_t`$, where each $`O_t`$ is Hermitian unitary,
each $`c_t`$ is real, and $`\Lambda=\sum_t|c_t|\gt0`$. Classical
term selection has probabilities $`|c_t|/\Lambda`$ and retains the
coefficient sign. Its setup and sampling are charged. A controlled
$`O_t`$ is supplied or compiled with its literal branch phase calibrated
and complete initialized-isometry error at most $`2^{-K}`$ on arbitrary
branch-system inputs, including its declared work and external references.
Rounded weights require the additional allocation in the
[reflection-sum contract](QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights).

For exact term weights, fix $`0\lt\varepsilon_\infty\le\Lambda`$,
$`0\lt\delta\lt1`$, and

```math
K=\max\!\left\{6,\left\lceil\log_2
 \frac{80\Lambda}{\varepsilon_\infty}\right\rceil\right\},\qquad
P=\max\{n,K\},\qquad
S=\left\lceil\frac{800\Lambda^2[1+\ln(2(n+1)/\delta)]}
 {\varepsilon_\infty^2}\right\rceil.
```

The corollary estimates all original raw magnitude derivatives and, in a
complex chart, all N original raw leaf-phase derivatives of
$`\langle\psi|H|\psi\rangle`$ with simultaneous coordinate error
at most $`\varepsilon_\infty`$ and probability at least $`1-\delta`$.
The observable and its coefficients are held fixed when differentiating.
This is coordinatewise accuracy, not a bound on the full gradient's
Euclidean error or on metric-normalized derivatives.

## 2. Workspace, literal phase, and the preparation contract

Reserve two initialized compiler flags and
$`b\ge P+n+7`$ arbitrary dirty qubits. Separately reserve the n
initialized system qubits, one interference branch, and the controlled
observable's declared work. Thus “two flags” is not the total number of
initialized physical qubits. No further initialized compiler work or
supplied resource state is assumed.

For a complex chart use the supplied representatives to define

```math
\mu=N^{-1}\sum_x\varphi_x,\qquad D_0=e^{-i\mu}D_\varphi,\qquad
W'=D_0W_{\mathbb R},\qquad \psi'=W'|0^n\rangle.
```

For a real chart put $`\mu=0`$, $`D_0=I`$, and $`W'=W_{\mathbb R}`$.
The physical energy and requested derivatives are unchanged. The common
phase contribution to a leaf-phase energy derivative vanishes; all N
coordinates, including the redundant common-phase direction, are retained.
Keep this gauge, the parameter tuple, the compiled words, and the decoder
fixed throughout a batch. No synthesis word or decoder coefficient is
differentiated.

There is an actual recorded logical native unitary C with

```math
\|C-W'\|\le\min\{1/64,1/(4\sqrt N)\},\qquad
T(C)=O(N),\qquad G(C)=O(Nn).
```

Its physical implementation is exactly $`C\otimes I_b`$, and its
actual inverse is exactly $`C^\dagger\otimes I_b`$. This exact dirty
return holds on all logical inputs. In particular, C may be complex even
for a real target, and actual complex prefix rows cannot be replaced by
ideal diagonal phases in the decoder.

The fine preparation has error at most $`2^{-P}`$ as a complete
initialized isometry. For a single state its contract is

```math
\sup_{\|\xi\|=1}\left\|
 V(|0^n\rangle|00\rangle|\xi\rangle)
 -|\psi'\rangle|00\rangle|\xi\rangle\right\|\le2^{-P}.
```

The coherent version preserves an arbitrary protocol branch and prepares
$`C|0^n\rangle`$ on branch zero and psi' on branch one, with the
same isometry error and literal relative phase. On a plus branch its
ideal output is $`(|0\rangle C|0^n\rangle+|1\rangle\psi')/\sqrt2`$.
Both contracts include flag leakage, approximate dirty return, and
arbitrary external-reference correlations. They require initialized
logical input; they do not prescribe the action on other logical states.

## 3. Executions, records, and the guarantee

Use S magnitude executions. Prepare the coherent reference/target pair,
apply the selected controlled observable on branch one, and apply
$`U=H^{\otimes n}C^\dagger`$ to both branches. Independently choose
X or Y branch measurement with equal probability and measure the leaf x.
The exact reference becomes uniform because the actual inverse of C is
used. For $`d'_j=\partial_{\theta_j}\psi'`$, $`f_j=Ud'_j`$, and
branch eigenvalue s, the ideal magnitude record is

```math
Y_j=4\Lambda\operatorname{sgn}(c_t)s\sqrt N
\begin{cases}\operatorname{Re}f_{jx},&X,\\
\operatorname{Im}f_{jx},&Y.\end{cases}
```

Its mean is the original raw derivative. At every magnitude depth,
$`\|Y^{(d)}\|_2\le5\Lambda`$, uniformly in the supplied angles.
Coarse error bounds this record size; it adds no bias when the same actual
C is used in the reference, measurement, and correction.

A complex chart uses S additional phase executions: prepare psi', use a
plus branch and the controlled observable, then measure branch Y and the
computational leaf. The record
$`2\Lambda\operatorname{sgn}(c_t)s e_x`$ has the original phase
derivative as its mean and norm at most $`2\Lambda`$. This stream uses
no inverse coarse or inverse fine frame. There are exactly S controlled
observable calls in a real chart and $`2S`$ in a complex chart.

Choose fresh term labels independently of previous history; magnitude
X/Y choices are also fresh and independent of the term labels. With
preparation error $`2^{-P}`$ and oracle error $`2^{-K}`$, magnitude
depth bias is at most $`\varepsilon_\infty/4`$. Allocate at most
$`\varepsilon_\infty/4`$ to deterministic reconstruction rounding.
The retained martingale concentration and union bound over n depths and
the optional phase block prove the guarantee in Section 1.

The dirty pool may be reused and correlated with past executions. Initialize
the system and declared clean inputs between completed executions; there
is no intermediate reset or postselection. The guarantee uses each
conditional all-input error bound, without assuming independent dirty
states or claiming that total dirty disturbance stays below one-shot error.

## 4. Quantum and classical resource ledger

Both single-state and coherent-pair preparation have per-use bounds

```math
A_{\rm state}=N+P,\qquad T=O(A_{\rm state}),\qquad G=O(NP).
```

At the stronger literal reservation $`b\ge2(P+n+7)`$, either may use

```math
B_{\rm state}=\sqrt{NP}+P+\frac{NP}{b}+n\sqrt N,\qquad
T=O(B_{\rm state}),\qquad G=O(NP).
```

At this reservation the exact coarse bound is
$`T(C)=O(Nn/(n+b)+n\sqrt N)`$; its first term is absorbed by
$`NP/b`$. The final term therefore remains charged. It may be absorbed
when $`P\ge n^2`$; it is not generally omitted. The additional actual
coarse inverse in magnitude readout changes these orders only by a
constant factor. Choose the smaller eligible construction; these are
constructive upper bounds with independent hidden constants.

Let $`t_{\rm pair},t_{\rm single},t_C`$ be actual T-counts, including
each preparation's forward coarse call, and let
$`\overline t_O=\sum_t|c_t|t_{O_t}/\Lambda`$ at oracle precision K.

| Chart | Expected T-count for all gradient streams |
|---|---|
| Real | $`S(t_{\rm pair}+t_C+\overline t_O)`$ |
| Complex | $`S(t_{\rm pair}+t_{\rm single}+t_C+2\overline t_O)`$ |

The Clifford count is $`O(S(NP+\overline g_O+n))`$ in either case,
with the fixed two-stream factor included for the complex chart. Oracle
precision is K, not P. Replacing expected oracle counts by their maxima
gives a deterministic count bound.

Signed integer histograms, an exact integer Walsh transform, actual-C
application, and one division-free real-tree reverse traversal reconstruct
the empirical gradients in $`O(S+Nn)`$ arithmetic operations and
$`O(N)`$ stored values after coefficients exist. The phase histogram
adds $`O(S+N)`$. Counters use $`O(1+\log S)`$ bits; a separate guarded
fractional precision $`D_{\rm dec}=O(n+K+\log(n+1))`$ suffices for
final rounding, rather than simply P fractional bits. Input
evaluation, arithmetic bit costs, term sampling, and output bits remain
charged; ordinary floating-point utilities are not accuracy certificates.

### Depth schedules for the same task

The [state-based depth proof](STATE_QBP_DEPTH.md) gives two schedules,
with $`B_0=P+n+7`$. Each row describes one circuit, including forward
coarse preparation and the additional inverse in magnitude readout:

| Sufficient dirty work | Compiler T-count | Compiler T-depth |
|---|---|---|
| $`b\ge2B_0`$ | $`O(NP)`$ | $`O(NP/b+P+n^3)`$ |
| $`b\ge16(B_0+\sqrt{NP})`$ | $`O(\sqrt{NP}+P+n\sqrt N)`$ | $`O(P+n^4)`$ |

Both have $`G=O(NP)`$, two compiler flags, and the unchanged state/error
contract. The first row deliberately permits a higher T-count. Oracle
depth at precision K and all S or $`2S`$ executions remain charged in the
[serial task ledger](QBP_COST_COMPARISON.md#7-state-based-t-depth-comparison).
These are T-depth upper bounds; Clifford and total elementary depth are
separate, and no matching depth lower bound is asserted.

## 5. Construction, proof dependencies, and scope

For bounded dyadic angles/phases in $`[-8,8]`$ with at most B fractional
bits, and an explicit M-term Pauli list with bounded B-bit coefficients,
[the computational corollary](BOUNDED_INPUT_QBP.md#4-a-uniform-construction-and-explicit-output-bound)
gives deterministic preprocessing polynomial in N, B, P, M, and wire-label
length, including static tables and native-program output: for
$`n\ge6`$ at the basic reservation, or every $`n\ge1`$ at the banked
reservation. Native output has $`O(NPJ+MnJ)`$ bits for
$`J=\lceil\log_2(n+b+K+3)\rceil`$. The basic-budget fallback for
$`n\le5`$ retains a sufficient $`2^{O(P)}`$ search bound; the quantum
gate theorem still holds there. No polynomial evaluator bound is inferred
for arbitrary computable inputs, and no QRAM or factoring oracle is used.
Explicit controlled Pauli terms use $`O(n)`$ Clifford gates each;
their descriptions, generation, and execution are still charged.

The implemented [residual coefficient helper](RESIDUAL_TABLE_PREPROCESSING.md)
certifies classical rotation data. Complete native [real](NATIVE_COARSE_QBP.md)
and [complex](NATIVE_COMPLEX_COARSE_QBP.md) examples verify special two-qubit
targets. The general fine-precision native emitter is not implemented;
the construction bound above is analytic.

| Ingredient | Canonical proof |
|---|---|
| State isometry and literal common reference | [State compiler §§1–7](STATE_ONLY_COMPILER.md), especially [§7](STATE_ONLY_COMPILER.md#7-a-complex-native-reference-from-the-same-coarse-word) |
| Complex gauge, exact coarse return, and additional banks | [Complex compiler §§1–8](COMPLEX_COARSE_COMPILER.md) |
| Magnitude means, depth bound, and real histogram | [Real QBP §§1–5](COARSE_FRAME_QBP.md) |
| Original complex phase gradients, confidence, and certified decoder | [Complex QBP §§1–6](COMPLEX_COARSE_QBP.md) |
| Complete compiler T-depth schedules | [State-based depth proof](STATE_QBP_DEPTH.md) |
| Reused dirty-work concentration and observable access | [QBP approximation §10](QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights) |
| Algebraic residual tables and polynomial construction qualifiers | [Residual preprocessing](RESIDUAL_TABLE_PREPROCESSING.md) and [bounded-input §§3–4](BOUNDED_INPUT_QBP.md#3-which-native-searches-are-polynomial-in-the-input-parameters) |

The result prepares and measures the required states; it does not compile
the fine prescribed complete Hopf frame or close that frame's endpoint.
The [fair quantum comparison](QBP_COST_COMPARISON.md) retains original
compilers that can be cheaper in some regimes; no T-count or shot
optimality is asserted. Under explicit Pauli-list access, the
[classical all-gradient and term-sampled baselines](BOUNDED_INPUT_QBP.md#2-a-deterministic-classical-all-gradient-baseline)
must also be charged and compared. A favorable quantum T upper bound alone
is neither a total-runtime speedup nor a quantum lower-bound statement.
