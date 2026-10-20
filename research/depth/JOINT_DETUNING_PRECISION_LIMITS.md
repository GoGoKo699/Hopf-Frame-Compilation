# Joint detuning and the cost of independent fine-angle replacement

[Paid fixed-angle detuning](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#52-joint-resynthesis-pays-one-special-angle-detuning-call) ·
[Polynomial-overhead discovery target](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Phase-signal interfaces](PHASE_SIGNAL_QUERY_LIMITS.md)

The fixed-angle bank at pi/4 has a paid joint native resynthesis. Replacing
that angle by a finer dyadic angle does not supply the same native primitive.
This chapter gives a complete arbitrary-angle fallback ledger and proves
a narrower stopping rule: independently synthesized fine rotations, when
certified as a complete tensor layer, retain at least logarithmic native
depth in the address length. Fixed finer angles have a stronger linear
depth cost in that representation.

This is not a bound on joint non-Clifford synthesis of the bank or the
complete selective pulse. It permits a common Clifford compute/uncompute,
but excludes additional synthesis helpers, cross-coordinate non-Clifford
gates and cancellation of errors across later pulse layers. The full
resource frontier is unchanged. The next operation is a complete selective
word with genuinely joint native synthesis and an error/return proof.

## 1. The arbitrary-angle component and its paid recurrence

Let b be the address length, $`m=2^b`$, and let x, Y and the m signal
wires $`s_u`$ be arbitrary, including reference entanglement. Define

```math
d_u=|x\mathbin\oplus u|+1-Y_u,
\qquad
R_\theta=\prod_u\exp(i\theta d_u Z_{s_u}).
```

The literal identity

```math
\begin{aligned}
R_\theta={}&
 \prod_u\exp\!\left(\frac{i\theta(b+1)}2Z_{s_u}\right)
          \exp\!\left(\frac{i\theta}2Z_{Y_u}Z_{s_u}\right)\\
&\quad\cdot
 \prod_{j,u}\exp\!\left(-\frac{i\theta(-1)^{u_j}}2
                         Z_{x_j}Z_{s_u}\right)
\end{aligned}                                                    \qquad(1)
```

has no omitted scalar. The local axes in the first line form one
independent set: compute $`Y_u\leftarrow Y_u\oplus s_u`$, apply the
two types of one-qubit rotation on their respective Y and s coordinates
in parallel, and uncompute. This costs 2m CNOTs. For each j, compute
$`s_u\leftarrow s_u\oplus x_j`$ on all signals, apply the indicated
rotations on their distinct signal coordinates, and uncompute. These b
rounds cost 2mb CNOTs. No initialized copy or additional helper is used.

Put $`q=m(b+2)`$. Let $`t_\epsilon,g_\epsilon,d_\epsilon`$ be
common forward/actual-inverse native bounds for literal ancilla-free
one-qubit rotations at operator error epsilon, uniformly over the angles
in (1). These are synthesis costs, not free rotation gates. Telescoping
the q replacements gives the simultaneous bounds

```math
\begin{aligned}
T_R&\le q\,t_\epsilon,\\
G_R&\le q\,g_\epsilon+2m(b+1),\\
D_T(R)&\le(b+1)d_\epsilon,\\
b_{\mathrm{signal}}&=m,\qquad
b_{\mathrm{additional}}=a_{\mathrm{clean}}=0,\\
\|\widetilde R_\theta-R_\theta\|&\le q\epsilon.
\end{aligned}                                                    \qquad(2)
```

The actual adjoint of the emitted word has the same error and common
resource bounds; positive-S costs for inverse Clifford gates belong in
$`g_\epsilon`$. Equation (2) is a full-operator statement. Although
the ideal component is diagonal, it can entangle arbitrary signals with
the data. It does not return every borrowed signal to its initial state.
Signal return is an obligation of the completed selective word.

For clarity, consider a word with h such detuning calls, tolerances
$`\epsilon_j`$, and interlayers whose common resources are
$`T_j,G_j,d_j`$ and operator errors $`\rho_j`$. Reserve one common
arbitrary helper space of width w for those interlayers; (2) requires
none beyond the m signals. On that same space,

```math
\begin{aligned}
T_W&\le q\sum_{j=1}^h t_{\epsilon_j}+\sum_jT_j,\\
G_W&\le q\sum_{j=1}^h g_{\epsilon_j}
          +2hm(b+1)+\sum_jG_j,\\
D_T(W)&\le(b+1)\sum_{j=1}^h d_{\epsilon_j}+\sum_jd_j,\\
b_W&\le m+w,\qquad a_{\mathrm{clean}}=0,\\
\|\widetilde W-W\|&\le
 \eta:=q\sum_{j=1}^h\epsilon_j+\sum_j\rho_j.
\end{aligned}                                                    \qquad(3)
```

The interlayers must themselves satisfy the stated zero-clean contract;
their completion and sharing are not supplied by this formula. No
intermediate return is inferred from the norm estimate.

If an ideal selective word separately proves

```math
\|Z_S W^\dagger Z_S W-D_b\otimes I_{S,w}\|\le\delta_0,
\qquad
Z_S=\prod_uZ_{s_u},\qquad
D_b|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle,
```

then using the actual emitted adjoint and Hadamards on Y yields a complete
indicator with

```math
\begin{aligned}
\epsilon_{\mathrm{ind}}&\le\delta_0+2\eta,\\
T_{\mathrm{ind}}&\le2T_W,\qquad
G_{\mathrm{ind}}\le2G_W+6m,\\
D_T(\mathrm{ind})&\le2D_T(W),\qquad
b_{\mathrm{ind}}\le m+w.
\end{aligned}                                                    \qquad(4)
```

Here the two $`Z_S`$ words cost 4m positive-S Cliffords and the two
output-Hadamard words cost 2m Cliffords. The norm includes return of all
signals and helpers, all logical columns and arbitrary references.
Equation (4) is conditional on the ideal selective estimate, which a
detuning component alone does not establish.

## 2. Finer dyadic banks are not exact native primitives

For $`\theta=\pi/2^r`$ with integer $`r\ge3`$, no finite
Clifford+T word implements $`R_\theta`$ exactly with returned standard
work, even allowing a common scalar. Fix x and all signals to zero and
compare two computational inputs differing only in one Y bit. The ratio
of their nonzero diagonal entries is

```math
e^{-i\theta}=e^{-i\pi/2^r}.
```

Every entry of a finite native word is in $`\mathbb Q(e^{i\pi/4})`$.
Returned work preserves this statement for the relevant matrix entries,
and their ratio also lies in that field. But the displayed primitive
$`2^{r+1}`$-th root has field degree $`2^r>4`$. A common scalar
cancels in the ratio. This proves the exact obstruction.

Approximation is a different question. This field argument supplies no
approximate depth bound for arbitrary joint native words. In particular,
it does not prevent exact joint resynthesis of a completed pulse whose
intermediate fine angles cancel.

## 3. Full tensor-layer accuracy forces summed local precision

Let $`P_\theta=\mathrm{diag}(1,e^{i\theta})`$. The following
lemma allows different local implementations, either sign of theta,
and an arbitrary common scalar. Suppose

```math
\left\|e^{i\gamma}\bigotimes_{i=1}^m V_i
              -\bigotimes_{i=1}^m P_{\sigma_i\theta}\right\|
\le\delta,\qquad
0\lt\delta\lt1,\quad \sigma_i\in\{1,-1\}.
```

Write

```math
\beta=2\arcsin(\delta/2)\lt\pi/3,
\qquad
\eta_i=\min_\phi\|V_i-e^{i\phi}P_{\sigma_i\theta}\|.
```

Then

```math
\sum_i\eta_i\le\beta.
                                                               \qquad(5)
```

**Phase-wrap proof.** Let the two eigenvalues of
$`P_{\sigma_i\theta}^\dagger V_i`$ be
$`\lambda_{i,0},\lambda_{i,1}`$, and let
$`\vartheta_i`$ be the principal argument of their ratio. Every
tensor eigenvalue, including the common scalar, has a unique argument
$`\phi_A\in[-\beta,\beta]`$, indexed by the subset A of factors
using their second eigenvalue. Along an edge of this subset cube,

```math
\phi_{A\cup\{i\}}-\phi_A\in[-2\beta,2\beta]
\subset(-\pi,\pi).
```

This difference must therefore equal the same real number
$`\vartheta_i`$, with no possible multiple of $`2\pi`$ added.
Consequently
$`\phi_A=\phi_\varnothing+\sum_{i\in A}\vartheta_i`$ as
real numbers, and the range over all subsets has width
$`\sum_i|\vartheta_i|\le2\beta`$. For a two-eigenvalue unitary,

```math
\eta_i=2\sin(|\vartheta_i|/4)
\le |\vartheta_i|/2.
```

Summation proves (5). In particular, some coordinate has
$`\eta_i\le\beta/m`$. This argument uses the complete tensor
operator promise, not accuracy on one initialized vector.

## 4. A one-qubit native algebraic gap

Consider a one-qubit Clifford+T word containing t T or inverse-T gates.
Its Bloch rotation matrix has entries

```math
x=\frac{a}{(\sqrt2)^t},\qquad a\in\mathbb Z[\sqrt2].
                                                               \qquad(6)
```

Arbitrary one-qubit Cliffords give signed permutation matrices. Under
the field conjugation $`\sqrt2\mapsto-\sqrt2`$, each T or
inverse-T Bloch matrix remains real orthogonal. Thus both embeddings
of every entry satisfy $`|x|,|x'|\le1`$.

Put $`c=\cos(\pi/8)`$, $`s=\sin(\pi/8)`$ and
$`L=\max(t,2)`$. Since c is not in $`\mathbb Q(\sqrt2)`$,

```math
z=2^L(x^2-c^2)\in\mathbb Z[\sqrt2]\setminus\{0\}.
```

The integrality uses $`c^2=(2+\sqrt2)/4`$. Its conjugate is
$`2^L((x')^2-s^2)`$, and the latter parenthesis has absolute value
at most one. The nonzero integer norm of z therefore gives

```math
|x^2-c^2|\ge2^{-2L},\qquad
|x-c|\ge2^{-2L-1}.
                                                               \qquad(7)
```

Projective operator error eta between one-qubit unitaries bounds each
Bloch-entry error by $`2\eta`$: choose the minimizing common phase
and telescope their adjoint actions on a Pauli. Applying (7) to the
cosine entry of a target $`P_{\pm\pi/8}`$ yields

```math
\eta\ge2^{-2\max(t,2)-2}\ge2^{-2t-6}.
                                                               \qquad(8)
```

For m independently compiled pi/8 rotations, (5), (8) and convexity
therefore imply the simultaneous necessary bounds

```math
\sum_i t_i\ge\frac m2\left(\log_2\frac m\beta-6\right),
\qquad
\max_i t_i\ge\frac12\log_2\frac m\beta-3.
                                                               \qquad(9)
```

On one physical qubit its T gates are sequential. The depth of the
independent layer is at least $`\max_i t_i`$, even though different
coordinates run together. The $`\Omega(mb)`$ count in (9) is compatible
with the relaxed polynomial-overhead discovery budget; the
$`\Omega(b)`$ depth does not improve that frontier.

## 5. Shrinking dyadic angles still retain logarithmic depth

Let $`m=2^b`$, $`b\ge6`$, and let
$`\theta=\pi/2^r`$ with integer $`3\le r\le b/2`$.
For an independent tensor layer satisfying the complete promise in
Section 3, some local word with t T gates obeys

```math
t\ge\max\left\{
 \frac{b-2r-\log_2\beta}{2^{r-2}},\quad
 r-\log_2(3\pi/2)
\right\}.
                                                               \qquad(10)
```

**Precision bound.** Put $`n=2^{r-3}`$ and let $`C_n`$ be the
degree-n Chebyshev polynomial, so
$`C_n(\cos\theta)=\cos(\pi/8)`$. For the actual cosine entry
x in (6), set $`y=C_n(x)`$. The polynomial has integer coefficients
and only powers with the same parity as n. Hence

```math
y=\frac{A}{(\sqrt2)^{nt}},\quad A\in\mathbb Z[\sqrt2],
\qquad |y|,|y'|\le1.
```

The last bounds follow because both x embeddings lie in
$`[-1,1]`$. The argument of (7), with denominator exponent nt,
gives

```math
|y-c|\ge2^{-2\max(nt,2)-1}\ge2^{-2nt-5}.
```

On $`[-1,1]`$, the Chebyshev derivative has absolute value at most
$`n^2`$. If the local projective error is eta, then

```math
|y-c|\le n^2|x-\cos\theta|\le2n^2\eta,
\qquad
\eta\ge\frac{2^{-2nt-6}}{n^2}.
```

Choose the coordinate from (5) with $`\eta\le\beta2^{-b}`$.
Taking logarithms and substituting $`n=2^{r-3}`$ gives the first
bound in (10).

**Angle-resolution bound.** For the same coordinate,
$`\eta\le\sin\theta/4`$. Indeed,
$`\sin\theta\ge2^{1-r}`$ and
$`\beta\lt\pi/3\lt2^{b-r-1}`$ in the stated range.
Its actual sine entry v is therefore nonzero and satisfies
$`|v|\le(3/2)\sin\theta`$. From (6),

```math
vv'\in2^{-t}\mathbb Z\setminus\{0\},\qquad |v'|\le1,
```

so $`|v|\ge2^{-t}`$. Since
$`\sin\theta\le\pi2^{-r}`$, this proves the second bound
in (10). The sign of the target angle does not change either argument.

For $`r\ge\tfrac12\log_2 b`$, the second bound is
$`\Omega(\log b)`$. Below that threshold the first bound is
$`\Omega(\sqrt b)`$ and is stronger. Thus independent replacement
has $`\Omega(\log b)`$ native depth uniformly over the stated
dyadic range, including every polynomially small angle
$`r=O(\log b)`$. At any fixed $`r\ge3`$, its first bound is
$`\Omega(b)`$. No claim here covers arbitrarily smaller angles
outside the range in (10).

## 6. Application and exact scope of the stopping rule

Clifford conjugation preserves operator norm. The tensor-layer theorem
therefore also applies after any common Clifford compute/uncompute that
maps independent Pauli axes to distinct coordinates. It applies, for
example, to independently replacing the fine rotations on the
$`Y_u\oplus s_u`$ coordinates of (1), when that complete local layer
is certified to the stated tolerance. If the other m local rotations
are present, fixing one eigenvalue of each of their relative-error
factors gives the same m-coordinate promise up to a common scalar.
An intact independent edge batch has the same interpretation.

The premise is an independent tensor implementation with a complete-layer
accuracy certificate. It does not follow merely from a final pulse being
accurate: intermediate approximation errors could cancel. It also does
not include synthesis helpers or joint non-Clifford operations between
coordinates. The pi/4 joint construction already demonstrates why an
unchanged list of rotations and a resynthesized product must be treated
separately.

Resolving only the scheduling of independent finer rotations therefore
does not reach the sublogarithmic source target. A possible improvement
must give a genuinely joint native bank operation or resynthesize the
complete selective pulse, with all its count, width, Clifford, depth,
error and return costs on the same circuit. Neither the exact field
obstruction nor the independent-layer lower bound supplies an
unrestricted depth obstruction. These are analytic statements; no new
numerical fixture is needed to establish them.
