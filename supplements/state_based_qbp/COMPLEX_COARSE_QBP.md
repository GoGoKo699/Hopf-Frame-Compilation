# Complex Hopf gradients with an actual coarse measurement frame

[Task theorem](STATE_BASED_QBP_THEOREM.md) · [Complex coarse compiler](COMPLEX_COARSE_COMPILER.md) · [Real coarse-frame decoder](COARSE_FRAME_QBP.md) · [Original complex-gradient contract](../../docs/QBP_APPROXIMATION.md#9-the-complete-complex-gradient)

The actual gauge-fixed complex coarse circuit supplies the same bounded
magnitude records as the real construction. A separate computational-basis
phase stream retains the original signed one-hot record. Together they
estimate every raw magnitude and leaf-phase derivative, including singular
tuples, with no derivative-envelope sampling penalty.

With two initialized compiler flags, the stated dirty reservation, and a
separately reserved protocol branch, each execution costs
$`O(N+L')+t_O`$ T gates. The classical reconstruction uses signed integer
histograms, the recorded coarse circuit, and a real Hopf reverse traversal.
This is a state-preparation and measurement result. It does not compile the
fine prescribed complete frame, establish a universal advantage, or supply
a general native emitter for this complex protocol. The separate
[bounded native example](NATIVE_COMPLEX_COARSE_QBP.md) now emits both
streams completely for a special exact two-qubit target.
The [native residual integration](NATIVE_RESIDUAL_QBP.md) also emits both
streams for a certified complex one-qubit target, using fine residual
preparation, charged common C and its actual inverse, and exact rational
histogram reconstruction.

## 1. Supplied phases, the fixed gauge, and the actual native interface

Fix $`N=2^n`$, $`n\ge1`$, real Hopf angles theta, and supplied real
leaf-phase representatives $`\varphi_x`$. Write

```math
\begin{gathered}
|\psi_{\mathbb R}\rangle=W_{\mathbb R}|0^n\rangle,
\qquad D_\varphi=\mathrm{diag}(e^{i\varphi_x}),
\qquad \mu=N^{-1}\sum_x\varphi_x,\\
D_0=e^{-i\mu}D_\varphi,\qquad W'=D_0W_{\mathbb R},
\qquad |\psi'\rangle=W'|0^n\rangle.
\end{gathered}
```

The energy of psi' equals that of the originally specified physical state
$`D_\varphi\psi_{\mathbb R}`$. The
[native coarse theorem](COMPLEX_COARSE_COMPILER.md#1-target-gauge-and-workspace)
constructs an actual recorded logical unitary C such that

```math
\|C-W'\|\le\delta_c,
\qquad \delta_c\le\min\{1/64,1/(4\sqrt N)\},
\qquad T(C)=O(N),\quad G(C)=O(Nn).
```

Its physical word implements $`C\otimes I_b`$ exactly on arbitrary dirty
inputs and external references; its actual inverse implements
$`C^\dagger\otimes I_b`$. The recorded logical action can be applied
classically in $`O(Nn)`$ two-mode operations. Synthesized phase-prefix
rows act on every suffix, and their actual native blocks need not remain
diagonal. They must not be replaced by ideal diagonal phases when applying
C or its inverse.

Keep the tuple, supplied phase representatives, chosen mu, actual coarse
word, and decoder fixed throughout each batch. No native synthesis decision
or measurement-basis coefficient is differentiated. A different winding
representative can change the chosen scalar gauge; it does not change the
physical energy, but every branch and classical coefficient in that batch
must use the same chosen gauge.

Let the fixed observable be $`H=\sum_t c_tO_t`$, with Hermitian-unitary
terms, real coefficients, and $`\Lambda=\sum_t|c_t|\gt0`$. Draw t with
probability $`|c_t|/\Lambda`$ and retain its sign. The controlled observable
must have its literal branch phase calibrated. The
[retained coefficient and oracle contracts](../../docs/QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights)
apply; preparing and sampling the classical term distribution is charged.

## 2. Magnitude means and the unchanged depth bound

Since mu and the leaf phases do not depend on the magnitude coordinates,

```math
|d'_j\rangle=\partial_{\theta_j}|\psi'\rangle
 =D_0|d_j\rangle=a_jW'|\lambda(j)\rangle.
```

Here $`a_j`$ is the signed real incoming amplitude. At each depth,
$`\sum_j a_j^2=1`$, including at singular tuples. Put

```math
|r\rangle=C|0^n\rangle,\qquad U=H_nC^\dagger,
\qquad H_n=H^{\otimes n},\qquad |f_j\rangle=U|d'_j\rangle.
```

Prepare $`(|0\rangle_Br+|1\rangle_B\psi')/\sqrt2`$, apply the
selected controlled observable on branch one, then U on both branches.
The reference becomes exactly $`H_n|0^n\rangle`$. Choose an independent
fair X or Y measurement on B and measure the system leaf x. If s is the
branch eigenvalue, the records are

```math
Y_j(q,s,x)=4\Lambda\,\mathrm{sgn}(c_t)s\sqrt N
\begin{cases}
\mathrm{Re}\,f_{jx},&q=X,\\
\mathrm{Im}\,f_{jx},&q=Y.
\end{cases}
```

For $`v=UO_t\psi'`$, the signed leaf probabilities conditioned on X
and Y are respectively $`\mathrm{Re}\,v_x/\sqrt N`$ and
$`\mathrm{Im}\,v_x/\sqrt N`$. Averaging the fair basis choice gives

```math
\mathbb E Y_j
=2\,\mathrm{Re}\langle d'_j|H|\psi'\rangle
=\partial_{\theta_j}\langle D_\varphi\psi_{\mathbb R}|
 H|D_\varphi\psi_{\mathbb R}\rangle.
```

For the norm bound put $`E=C^\dagger W'-I`$. If $`J_d`$ injects the
marker labels at depth d and $`A_d=\mathrm{diag}(a_j)`$, their
transformed derivative columns form

```math
F_d=H_n(I+E)J_dA_d.
```

Every row of the first term has norm $`1/\sqrt N`$. The error term has
row norm at most $`\|E\|\|A_d\|\le\delta_c`$. Taking real or imaginary
parts cannot increase the norm, so

```math
\|(f_{jx})_{j\text{ at depth }d}\|_2\le\frac5{4\sqrt N},
\qquad \boxed{\|Y^{(d)}\|_2\le5\Lambda.}
```

This is a bound for each depth, not for the concatenated magnitude vector.
There is no division by an incoming amplitude and no coarse-error bias:
the same actual C defines the reference, measurement basis, and scores.

## 3. The separate original phase stream

Prepare psi', then prepare B in plus, apply the selected controlled
observable, and measure B in Y and the system in the computational basis.
There is no inverse frame or inverse coarse circuit in this stream. Use

```math
P=2\Lambda\,\mathrm{sgn}(c_t)s\,e_x,
\qquad \|P\|_2=2\Lambda.
```

The standard Pauli-Y convention gives

```math
\mathbb E P_x
=2\,\mathrm{Im}\bigl(\overline{\psi'_x}(H\psi')_x\bigr)
=\partial_{\varphi_x}\langle D_\varphi\psi_{\mathbb R}|
 H|D_\varphi\psi_{\mathbb R}\rangle.
```

Indeed $`\partial_{\varphi_x}\psi'=i(|x\rangle\langle x|-I/N)\psi'`$.
The scalar term contributes zero to the energy derivative because
$`\langle\psi'|H|\psi'\rangle`$ is real. Equivalently, freeze the scalar
gauge at the tuple when forming the local differential records. This keeps
all N original phase coordinates, including the redundant common-phase
direction; no phase metric is inverted. In expectation their sum is zero.

Both the phase-stream distribution and its mean are invariant under a
common scalar phase. This does not license changing only one branch of the
magnitude interference state. That would change a measurable relative
phase and, in general, its gradient means.

## 4. Joint precision, confidence, and reused dirty inputs

Require complete preparation isometry error at most $`\eta_P`$ in either
stream and controlled-observable error at most $`\eta_O`$. The magnitude
preparation contract includes its literal relative plus sign, both compiler
flags, all dirty inputs, and external references. The phase preparation
requires the same complete work guarantee. Unitary measurement wrappers do
not increase the preceding error.

Bounded decoded observables and Euclidean duality give magnitude depth
bias at most $`10\Lambda(\eta_P+\eta_O)`$ and whole-phase-vector bias at
most $`4\Lambda(\eta_P+\eta_O)`$. For
$`0\lt\varepsilon_\infty\le\Lambda`$, $`0\lt\delta\lt1`$, choose

```math
\eta_P,\eta_O\le\frac{\varepsilon_\infty}{80\Lambda},\qquad
S=\left\lceil
\frac{800\Lambda^2[1+\ln(2(n+1)/\delta)]}
 {\varepsilon_\infty^2}\right\rceil
```

executions of **each** stream. Allocate at most
$`\varepsilon_\infty/4`$ to deterministic final reconstruction rounding
and $`\varepsilon_\infty/2`$ to sampling. The magnitude bias is at most
$`\varepsilon_\infty/4`$ and phase bias at most
$`\varepsilon_\infty/10`$.

Fresh term choices, and magnitude X/Y choices, are independent of the
previous history; term and basis choices are also mutually independent.
Conditional centered records have norm at most $`10\Lambda`$ for a
magnitude depth and $`4\Lambda`$ for the entire phase vector. The
[retained Hilbert-space martingale inequality](../../docs/QBP_APPROXIMATION.md#101-reusing-a-dirty-bank-across-executions)
therefore bounds each block's sampling tail by

```math
2\exp\!\left(-\frac{S\varepsilon_{\rm stat}^2}{200\Lambda^2}\right).
```

Union over the n magnitude depths and one phase block, with
$`\varepsilon_{\rm stat}=\varepsilon_\infty/2`$, proves simultaneous
raw-coordinate error at most $`\varepsilon_\infty`$ with probability at
least $`1-\delta`$. The all-input preparation and oracle contracts apply
to each conditional dirty state. Reusing the dirty pool consequently needs
no synthesis precision proportional to the number of executions. Initialize
the system and declared clean inputs between completed executions; no
intermediate reset, postselection, or independence of the dirty bank is
assumed. This does not bound its total accumulated disturbance by one
execution's error.

## 5. Integer histograms and the real-tree adjoint

Let A and B be the magnitude stream's signed integer leaf histograms for
X and Y, respectively, including the selected term sign. The explicit
number of executions is S; it cannot be recovered from signed histograms
after cancellations. Put

```math
h=\frac{4\Lambda\sqrt N}{S}(A+iB),\qquad
J'=D_0J_{\mathbb R}.
```

The exact empirical magnitude estimate is

```math
\widehat g_{\rm mag}
=\mathrm{Re}[(J')^\dagger C H_nh]
=J_{\mathbb R}^{\mathsf T}\,\mathrm{Re}(D_0^\dagger C H_nh).
```

Cancel square-root normalizations before numerical arithmetic. With the
integer Walsh matrix $`F_N=\sqrt N H_n`$, compute

```math
k=F_N(A+iB),\qquad u=k/S,\qquad
\boxed{\widehat g_{\rm mag}
=4\Lambda J_{\mathbb R}^{\mathsf T}\,\mathrm{Re}(D_0^\dagger Cu).}
```

Use exact integer butterflies separately on A and B. Since
$`\sum_x(|A_x|+|B_x|)\le S`$, every real or imaginary intermediate
counter is bounded in magnitude by S. Counters need $`O(1+\log S)`$
bits, $`|u_x|\le1`$, and $`\|u\|_2\le\sqrt N`$.

Apply the actual recorded coarse factors to u, respecting their order and
their action on every suffix. Multiply leaf x by
$`e^{i(\mu-\varphi_x)}`$ and take its real part to initialize
$`\beta_x`$. For $`c_j=\cos\theta_j`$, $`s_j=\sin\theta_j`$, compute
the signed incoming amplitudes forward from $`a_{\rm root}=1`$ and run

```math
\beta_j=c_j\beta_L+s_j\beta_R,\qquad
z_j=a_j(-s_j\beta_L+c_j\beta_R)
```

backward through the real tree. Output $`4\Lambda z_j`$. Every measured
leaf weight, C, and the phase factors are held fixed in this reverse
calculation. It does not differentiate their dependence on the tuple.

For the separate phase stream, let Q be its signed integer leaf histogram.
Its empirical estimate is simply $`\widehat g_{\rm phase}=2\Lambda Q/S`$.
Together the two histogram updates cost $`O(S)`$; Walsh and actual-C
application cost $`O(Nn)`$; phase multiplication, tree traversals, phase
output, and writing all coordinates cost $`O(N)`$. Total reconstruction
uses $`O(S+Nn)`$ arithmetic operations and $`O(N)`$ storage after the
recorded coefficients exist. There is no dense Jacobian or simulation of
the dirty Hilbert space.

The executable
[complex histogram utility](../../compiler_robust_hopf/complex_coarse_decoder.py)
exposes

```python
decode_complex_coarse_frame_histograms(
    theta, coarse_blocks, phase_blocks, leaf_phases,
    hist_x, hist_y, shots, *, coefficient_scale=1,
)
```

Here `coarse_blocks` records the actual coarse real-tree rows, and
`phase_blocks` records the actual prefix rows, in chronological heap order,
so that $`C=C_P C_R`$. The supplied `leaf_phases` fix the same arithmetic
mean gauge as the target. The result contains only the magnitude
coordinates; reconstruct the separate phase histogram as above. Walsh
counts use exact Python integers, while phase evaluation, native-block
application, and final contractions use ordinary NumPy floating point.
This reconstructs empirical records and does not implement the certified
arithmetic below, certify sampling accuracy, or emit a native circuit.
The [bounded complex checks](../../tests/test_complex_coarse_qbp.py) compare
both means and histogram contractions with independent analytic gradients.
Their local row blocks are native words; their assembled prefix selection
and residual amplification are matrix diagnostics.

## 6. Certified arithmetic and supplied-input costs

The preceding operation count needs a stable representation for its
rounding claim. The native coarse construction supplies
$`M=O(Nn)`$ logical two-mode updates, each unitary, including the actual
finite native row words. Certify each block coefficient, each sine/cosine,
and each unit phase in $`D_0^\dagger`$ to absolute error at most
$`e=2^{-P}`$, and round elementary scalar arithmetic to absolute error e.
Use sufficient integer range and clamp real trigonometric factors to
$`[-1,1]`$. For $`Me`$ smaller than a fixed constant, approximate block
norms have a bounded product.

Rounded divisions $`u=k/S`$ contribute $`O(\sqrt N e)`$ Euclidean error.
A telescoping product bound gives $`O(Me)`$ coarse operator error; acting
on a vector of norm at most $`\sqrt N`$ and including update rounding
therefore gives $`O(M\sqrt N e)`$ error. The final diagonal phase factors
add $`O(\sqrt N e)`$. Each exact raw derivative has norm at most one.
The division-free tree recurrence has n levels; its certified
trigonometric products and arithmetic are bounded conservatively by
$`O(Nn e)`$ per normalized output. Thus

```math
|\widetilde g_{{\rm mag},j}-\widehat g_{{\rm mag},j}|
\le O\!\left(\Lambda[N^{3/2}n+Nn]2^{-P}\right).
```

The phase histogram only needs certified division and final scaling. Do
the normalized computations before the final known factors $`4\Lambda`$
or $`2\Lambda`$, and evaluate those rescalings to absolute error
$`O(\Lambda2^{-P})`$. Consequently

```math
P=O\!\left(n+\log(n+1)+\log(\Lambda/\varepsilon_\infty)\right)
```

fractional bits, with suitable fixed constants and guard bits, suffice for
the deterministic coordinate rounding allocation in Section 4. This is a
certified-arithmetic existence bound, not a claim about ordinary NumPy
floating-point contractions. An arithmetic-operation count alone would
not establish it without the recorded stable factors.

Supply certified real phase representatives and angle evaluators. Their
description lengths, certified trigonometric evaluation, and range
reduction for large phase windings are charged. As functions of the supplied
phases, subtree means and mu have absolute condition number one in the
maximum norm; subsequent phase differences incur only a constant factor
of input error. Guard precision also covers arithmetic rounding through
the n levels of subtree means without changing the displayed order of P.
Evaluating their
unit phases needs no exact principal argument, phase unwrapping, or
division by a state amplitude. In particular zero-amplitude leaves create
no singular decoder condition. The cost of reading large integer parts is
not included for free in P.

Likewise charge certified evaluation of the actual native row coefficients,
term-distribution preprocessing, all arithmetic bit costs, and output
representation. A finite row word can be evaluated with guard precision
for its recorded length; none of these coefficients are a free supplied
quantum oracle. The [algebraic residual procedure](RESIDUAL_TABLE_PREPROCESSING.md)
programs the fine table directly from certified sine/cosine pairs, without
Euler search. Its square-root and coefficient bounds, and the
[bounded-input construction](BOUNDED_INPUT_QBP.md), separately price that
preprocessing. The histogram arithmetic count alone is not its bit bound.

## 7. Complete quantum ledger and the comparison boundary

Set

```math
L'=\max\{n,6,\lceil\log_2(80\Lambda/\varepsilon_\infty)\rceil\},
\qquad a=2,\qquad b\ge L'+n+7.
```

The [complex state and common-reference theorem](COMPLEX_COARSE_COMPILER.md#6-two-flag-complex-state-and-common-reference-preparation)
prepares psi' alone, or the coherent pair $`(C|0^n\rangle,\psi')`$,
with complete error at most $`\varepsilon_\infty/(80\Lambda)`$ and
$`T=O(N+L')`$, $`G=O(NL')`$. It uses the same actual C and complex
residual $`C^\dagger\psi'`$; it does not invoke a two-real-state theorem
on a complex target. The literal relative branch phase is included.

Reserve the initialized interference branch and the observable's own work
in addition to the two compiler flags. The streams run separately and can
reuse their workspace. If $`t_{\rm pair}`$ and $`t_{\rm state}`$ are their
respective preparation costs and bars denote term-sampling averages, then

```math
\begin{aligned}
\mathbb E\mathcal T_{\nabla,\mathbb C}
&\le S(t_{\rm pair}+t_C+t_{\rm state}+2\overline t_O)\\
&=O(S(N+L'))+2S\overline t_O,\\
\mathbb E\mathcal G_{\nabla,\mathbb C}
&=O(SNL')+2S\overline g_O.
\end{aligned}
```

There are three charged coarse appearances per pair of executions: C in
each preparation and $`C^\dagger`$ in magnitude readout. The phase stream
has no inverse C. Hadamards and branch readouts are included in the
Clifford bound. No QRAM, free prepared reference, supplied catalyst, or
intermediate reset is assumed.

The retained original complete-frame protocol has expected cost
$`S(3t_F+2\overline t_O)`$ at its applicable confidence and precision
allocation. The new result replaces that compilation task and decoder;
it does not establish uniform dominance, optimal gradient cost, or a
fine-frame endpoint. Its sufficient sampling constant also differs from
the original bound. With $`b\ge2(L'+n+7)`$, the
[banked preparation corollary](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation)
improves its per-execution T bound to
$`O(\sqrt{NL'}+L'+NL'/b+n\sqrt N)`$ apart from the observable.
The [complete task comparison](QBP_COST_COMPARISON.md) keeps the original
precision $`K=\max\{6,\lceil\log_2(80\Lambda/\varepsilon_\infty)\rceil\}`$
distinct from $`L'=\max(n,K)`$: the observable only needs precision K.
It identifies the high-precision improvement and the fixed-accuracy
limitation without claiming an end-to-end speedup.
Source attribution for Hadamard interference,
histogram differentiation, and dirty-input martingale concentration is
retained in the [real decoder](COARSE_FRAME_QBP.md),
[reference-state discussion](REFERENCE_STATE_QBP.md), and
[original approximation contract](../../docs/QBP_APPROXIMATION.md). The separate
native construction supplies the new complex interface used here.
