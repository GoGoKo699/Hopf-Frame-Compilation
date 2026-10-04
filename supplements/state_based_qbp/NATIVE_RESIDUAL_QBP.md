# Native residual preparation inside both QBP streams

[Coherent selector](NATIVE_RESIDUAL_BRANCH.md) · [Complex QBP decoder](COMPLEX_COARSE_QBP.md) · [Verification](../../docs/VERIFICATION.md)

The [bounded integration](../../compiler_robust_hopf/native_residual_qbp.py)
connects the precision-dependent residual emitter to an actual common
coarse circuit, a calibrated controlled observable, and both gradient
readout streams. The target has one system qubit and genuinely complex
residual coefficients. Every quantum gate is elementary Clifford+T;
the classical program and histogram decoders use exact rational arithmetic.
There is no intermediate reset, projection, or replacement of the actual
inverse by an ideal inverse.

The earlier [complex native example](NATIVE_COMPLEX_COARSE_QBP.md) uses
an exact finite-size preparation fallback. Here the preparation itself
uses the fine residual tables and their conditional error certificates.
This is a bounded executable integration of the existing theorem, not
a general native compiler or a new asymptotic result.

## 1. A target with an exact coarse certificate

Use the conventions

```math
R_y(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix},
\qquad R_z(\phi)=\operatorname{diag}(e^{-i\phi},e^{i\phi}).
```

Fix

```math
u=2\arctan(1/1024),\qquad v=2\arctan(1/256),
```

```math
C=R_z(\pi/4)R_y(\pi/4),\qquad
W'=R_z(\pi/4+v)R_y(\pi/4+u).
```

The real magnitude coordinate is $`\theta=\pi/4+u`$; the supplied
leaf phases are $`(-\phi,\phi)`$ with $`\phi=\pi/4+v`$.
Their mean gauge is exactly zero. These formulas specify the target;
no transcendental angle is evaluated when programming its native word.

For $`t=1/1024`$ and $`r=1/256`$, set

```math
c_u=\frac{1-t^2}{1+t^2},\quad s_u=\frac{2t}{1+t^2},\qquad
c_v=\frac{1-r^2}{1+r^2},\quad s_v=\frac{2r}{1+r^2}.
```

All four numbers are rational. The determinant-one relative matrix
has trace $`2c_uc_v`$, so its two eigenvalues give the exact certificate

```math
\|C-W'\|^2=2(1-c_uc_v)
=\frac{262144}{4042387697}\lt\frac1{4096}.
```

Thus C meets the actual 1/64 coarse threshold on the complete logical
space, not only on its prepared column. Chronological Z,H implements
$`R_y(\pi/4)`$; chronological X,TDG,X,T implements
$`R_z(\pi/4)`$ with literal phase. C costs two T/TDG gates and
acts only on the system, so its borrowed-work return is exact.

The true residual is

```math
C^\dagger W'|0\rangle=a|0\rangle+b|1\rangle,
\qquad a=c_vc_u+i s_vs_u,\quad b=c_vs_u+i s_vc_u.
```

Both complex pairs are rational and satisfy $`|a|^2+|b|^2=1`$
exactly. The state emitter uses $`w=\sqrt2 b`$, so
$`|a|^2+|w|^2/2=1`$ and $`|w|\lt1`$. The same squared
distance above places this residual inside the 1/64 neighborhood.

At precision q, certified rational endpoints enclose the single square
root of two. Midpoint multiplication and nearest-dyadic rounding give
the supplied a,w pairs with error below $`2^{-2q-20}`$ each.
Their exact error certificates discharge the coefficient promises for
this fixture. They do not substitute floating-point angles or a numerical
Euler solver. Fine-q checks construct words and certificates without
simulating their full Hilbert spaces.

## 2. Complete native stream words

The layout is the coherent selector's: core 0 through q, signal q+1,
flag t=q+2, system x=q+3, protocol branch c=q+4, and flag s=q+5.
Each execution initializes x,c,s,t to zero; the q+2 core/signal wires
are arbitrary and may be entangled with external references. The branch
is separately counted from the two clean compiler flags. Initialization
between completed executions does not permit resets inside a word.

Let $`\widehat A_{\rm pair}`$ select residuals $`|0\rangle`$ and
$`a|0\rangle+b|1\rangle`$ under c. Let
$`\widehat A_{\rm one}`$ prepare the latter residual, using the
unbranched emitter with its mode flag remapped to q+5. This leaves c
untouched and uses the same physical layout and dirty pool.
The fixed Hermitian-unitary observable is the system Hadamard,
$`O=H_x=(X_x+Z_x)/\sqrt2`$, with coefficient scale one.
Its branch-one controlled word is the exact two-T controlled Hadamard.
O is held fixed when taking every derivative.

| Chronological stage | Magnitude stream | Phase stream |
|---|---|---|
| Residual preparation | H on c, then the coherent pair emitter | Unbranched emitter; c remains zero |
| Common coarse circuit | C on x | C on x, then H on c |
| Observable | Controlled O on branch one | The same controlled O |
| System readout | Actual C inverse, then H on x | Computational basis |
| Branch readout | Fair independent choice of X or Y | Y |

X readout uses H on c; Y uses SDG followed by H. Both then measure c
in the computational basis, with sign $`s=(-1)^c`$. The system outcome
is the leaf $`x\in\{0,1\}`$. All flag and dirty outcomes are traced
over when forming these probabilities. No successful-flag postselection
or renormalization is allowed.

Both residual emitters retain their literal amplification sign and use
their actual Q inverse. Applying C, the exact observable, and the readout
wrappers preserves the corresponding preparation isometry error. The
phase stream uses no inverse coarse or inverse fine frame.

## 3. Both raw gradient decoders

Write $`a=A+iB`$ and $`b=D+iE`$, with rational A,B,D,E fixed above.
For $`|d_\theta\rangle=\partial_\theta W'|0\rangle`$, put
$`f=HC^\dagger d_\theta`$. Direct multiplication gives

```math
f=\frac1{\sqrt2}
 \begin{pmatrix}\overline a-\overline b\\-\overline a-\overline b\end{pmatrix}.
```

The two leaf-weight vectors for the fair X/Y magnitude experiment are
therefore exactly rational:

```math
m_X=(4(A-D),-4(A+D)),\qquad m_Y=(4(E-B),4(E+B)).
```

The record is $`s\,m_q(x)`$. Its ideal mean equals the magnitude
derivative; omitting Y generally changes that mean. For signed integer
leaf histograms $`h_X,h_Y`$ from S total magnitude executions,

```math
\widehat g_\theta=\frac{m_X\cdot h_X+m_Y\cdot h_Y}{S}.
```

S includes both independently chosen measurement settings; it is not
reconstructed from signed counts. The phase record is $`2s\,e_x`$.
For its separate S-execution signed histogram h,

```math
(\widehat g_{\varphi_0},\widehat g_{\varphi_1})=2h/S.
```

The two decoder functions return exact Fractions. They retain both raw
phase coordinates, including the redundant common-phase direction.
Their ideal sum is zero; a finite histogram is not projected onto that
constraint. Input checks require integral counts, a positive integral S,
and the necessary absolute-count budget.

Independent analytic differentiation gives

```math
g_\theta=\sqrt2\bigl[-(c_u^2-s_u^2)+4c_us_uc_vs_v\bigr],
```

```math
g_{\varphi_0}=\frac{(c_u^2-s_u^2)(c_v^2-s_v^2)}{\sqrt2},
\qquad g_{\varphi_1}=-g_{\varphi_0}.
```

All three are nonzero. The coefficient programming, quantum protocol,
and fixed histogram reconstruction refer to the same target, gauge,
observable, and actual coarse C.

## 4. Literal costs, errors, and evidence

Let $`T_{\rm pair}=3240q+3793+210(2k_z+k_y)`$ be the coherent
residual count and $`T_{\rm one}=3240q+3793`$ the unbranched count.
The complete emitted stream counts are

| Stream | Charged T/TDG gates |
|---|---|
| Magnitude, either measurement setting | T_pair + 2 for C + 2 for controlled O + 2 for actual C inverse |
| Phase | T_one + 2 for C + 2 for controlled O |

Every Clifford is also present in the word and included in its gate
inventory. The two magnitude settings are alternative executions,
not two successive circuits in one shot. With S executions of each
stream, their total T count is $`S(T_{\rm pair}+T_{\rm one}+10)`$.
Observable cost, forward C, and the additional magnitude inverse are
all charged. The dirty allocation remains L+12 at $`q=L+10`$,
four above the basic n=1 reservation. It fits the banked reservation
and does not replace the minimum-budget fallback.

For the respective certified preparation errors
$`\eta_{\rm mag},\eta_{\rm phase}\lt390\,2^{-q}`$, the existing
bounded-record argument gives

```math
|\mathbb E\widehat g_\theta-g_\theta|\le10\eta_{\rm mag},
\qquad
\|\mathbb E\widehat g_{\rm phase}-g_{\rm phase}\|_2
\le4\eta_{\rm phase}.
```

These bounds include arbitrary dirty/reference inputs and leakage.
The oracle and readout are exact. Sampling cost and confidence remain
those of the [two-stream theorem](COMPLEX_COARSE_QBP.md#4-joint-precision-confidence-and-reused-dirty-inputs);
the finite checks below evaluate means and do not claim sampled accuracy.

The [tests](../../tests/test_native_residual_qbp.py) compare emitted words
with independent logical matrices and source algebra, sum probabilities
over all work and flag outcomes, and check both decoded means.
Small q=5 native columns exercise coherent dirty/reference inputs;
fine precision uses exact certificates and gate ledgers. The q=5 analytic
error constants are loose, so independent native-oracle agreement and
nonvacuous numerical checks are kept separate from the fine-q proof.

The executable [ledger example](../../scripts/residual_qbp_native_example.py)
reports the exact fixture, resource, and bias certificates. No exact
native fine-frame circuit for W' is supplied, so this pass makes no
same-accuracy gate-count comparison with an original fine-frame word.
The general fine compiler, larger addressed schedules, and end-to-end
advantage remain separate questions.

```bash
python scripts/residual_qbp_native_example.py --q 16
python scripts/residual_qbp_native_example.py --q 16 --format json
```

The selected bounded residual-to-readout integration is complete. The
[coverage audit](../../docs/reference/VERIFICATION_CATALOGUE.md#state-based-qbp-coverage) records its
proof dependencies and the separate general software boundaries.
No larger fixed example is selected; a further implementation pass needs
a concrete Hopf-QBP use requirement identifying a missing interface.
