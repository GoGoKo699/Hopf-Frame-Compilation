# Dirty geometric phase loops and the complete-bank cost

[Approximate query transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Phase-signal interfaces](PHASE_SIGNAL_QUERY_LIMITS.md) ·
[Complete-input filter limits](APPROXIMATE_INDICATOR_FILTER_LIMITS.md)

A two-level geometric loop can produce the same Boolean sign on both
states of an arbitrary dirty signal. Thus the fixed-subspace obstruction
for a rank-one signal reflection is not an obstruction to every dirty
signal mechanism. This chapter proves that statement with a uniform
finite-time operator estimate, then charges a complete-bank digital
implementation. The resulting explicit construction is much more
expensive than the retained indicator. Its independent signal errors
obey a sharp sum condition; improving only a single initialized signal
would not pay the complete operation. R96 and the frame frontier remain
unchanged.

All phases below are literal. There is no measurement, initialized
signal, adiabatic ground-state promise, or supplied analog gate in the
native conclusion. The continuous evolution is an intermediate identity
whose digital cost is charged in Section 4.

The [selective-inversion continuation](SELECTIVE_INVERSION_DIRTY_INDICATOR.md)
improves the analog bank-error profile. Its native row-count obstruction
shows why this improvement still does not give a linear-count indicator.

## 1. A scalar echo on both signal eigenstates

Use X, Y and Z for the three Hermitian Pauli matrices on one signal
qubit. For $`a\in\{1/2,3/2,5/2,\ldots\}`$ define

```math
H_a(\theta)=\sin\theta\,X+(a+\cos\theta)Z,
\qquad 0\le\theta\le2\pi.
```

Put $`r_a(\theta)=\sqrt{1+a^2+2a\cos\theta}`$. The eigenvalues
are plus and minus r, and $`r\ge1/2`$ for every such a. Let
$`\phi_a(\theta)`$ be the continuous argument of
$`a+\cos\theta+i\sin\theta`$, starting at zero, and write
$`R(\phi)=e^{-i\phi Y/2}`$. Then

```math
\begin{gathered}
H_a=r_aR(\phi_a)ZR(\phi_a)^\dagger,\\
\phi_a(2\pi)=2\pi\nu_a,\qquad
\nu_a=\mathbf1[a=1/2],\\
A_a:=\int_0^\pi r_a(\theta)\,d\theta
     =\int_\pi^{2\pi}r_a(\theta)\,d\theta.
\end{gathered}                                                   \qquad(1)
```

The equal dynamic integrals are important: the second half is the
reflection of the first in the X coordinate. For a time scale t, the
ideal adiabatic transport on an interval is

```math
U_a^{\mathrm{ad}}(\theta_2,\theta_1)
=R(\phi_a(\theta_2))
 e^{-itZ\int_{\theta_1}^{\theta_2}r_a(\theta)\,d\theta}
 R(\phi_a(\theta_1))^\dagger.
```

The literal two-half echo, with the rightmost factor acting first, is

```math
E_a^{\mathrm{ad}}
=Y U_a^{\mathrm{ad}}(2\pi,\pi)
 Y U_a^{\mathrm{ad}}(\pi,0)
=R(2\pi\nu_a)=(-1)^{\nu_a}I.                       \qquad(2)
```

Indeed, Y commutes with every R and anticommutes with Z, so the equal
dynamic phases cancel. The remaining spinor rotation is minus identity
for the winding loop and identity otherwise. Both eigenstates receive
the same sign. This is a half-loop echo inside one full loop; replacing
it by two full winding loops would square away the desired sign.

Equation (2) is an explicit matrix identity, not an appeal to an
initialized eigenstate. Geometric spin-echo constructions are an
established ingredient; for example, see
[Jones et al., *Geometric quantum computation with NMR*](https://arxiv.org/abs/quant-ph/9910052).
The particular half-integer family, complete dirty scalar contract,
finite-time estimate and bank accounting here are derived below.

## 2. Uniform finite-time accuracy

Let $`U_a(\theta_2,\theta_1)`$ be the actual propagator satisfying
$`i\partial_\theta U=tH_a(\theta)U`$, and use those actual propagators
in the echo E_a. Define the absolute constant

```math
C=(2+7\pi)(1+\pi).
```

For every positive half-integer a and every positive t,

```math
\|E_a-(-1)^{\nu_a}I\|_{\mathrm{op}}\le\frac{2C}{t}.
                                                               \qquad(3)
```

Here is a direct two-dimensional estimate, including both signal
eigenstates. On either half interval move to the instantaneous frame R
and then to the interaction frame of $`tr_aZ`$. The remaining generator
K is anti-Hermitian, has norm $`|g|`$ with $`g=\phi_a'/2`$, and
its off-diagonal entries have the form

```math
\pm g(\theta)e^{\pm2itF(\theta)},\qquad F'=r_a.
```

The vector defining H has first and second derivative of norm one.
Consequently

```math
r_a\ge\tfrac12,\quad |r_a'|\le1,\quad
|\phi_a'|\le2,\quad |\phi_a''|\le10,
\quad |g|\le1,\quad |g'|\le5.                         \qquad(4)
```

For example, differentiating the angular velocity gives the bounds
$`|\phi_a'|\le1/r_a`$ and
$`|\phi_a''|\le1/r_a+2/r_a^2`$. Integration by parts on any
subinterval of a half loop yields

```math
\left|\int g(\theta)e^{2itF(\theta)}\,d\theta\right|
\le\frac{2+7\pi}{t}.                                  \qquad(5)
```

The two boundary terms contribute at most $`2/t`$;
$`|(g/r_a)'|\le14`$ bounds the integral term by $`7\pi/t`$.
The primitive B of K has the same operator-norm bound. If W is the
interaction-frame propagator, integration by parts in
$`W-I=\int KW`$ gives

```math
W-I=BW-\int BKW,
\qquad
\|W-I\|\le\frac{(2+7\pi)(1+\pi)}{t}.
```

Thus each actual half-loop differs from its displayed adiabatic
transport by at most C/t. Product telescoping across the two halves
proves (3). This estimate is independent of a even though
$`\|H_a\|`$ grows with a. It is not a native simulation estimate.

## 3. Complete bank and the exact error-sum requirement

Let $`p\ge1`$, $`m=2^p`$, x be the address and Y the arbitrary
m-bit output bank. Borrow m arbitrary signal qubits, one per label u,
and set

```math
a_u(x,Y)=|x\oplus u|+\tfrac32-Y_u.
```

The unique winding condition is $`x=u`$ and $`Y_u=1`$. Running
the m ideal echoes in parallel therefore gives

```math
\bigotimes_u E_{a_u}^{\mathrm{ad}}
=(-1)^{Y_x}I_{\mathrm{signals}}.                       \qquad(6)
```

Address and bank are computationally preserved in this intermediate
construction; each signal undergoes noncommuting basis changes.
The arbitrary joint signal state, including entanglement among signals
and with references, returns with the same desired phase. Hadamards on
the bank before and after give the complete dirty indicator.

The finite-time bank norm is not the maximum row error. A useful exact
statement is the following. Suppose the target-relative row words are
$`V_u\in SU(2)`$ with eigenvalues $`e^{\pm i\delta_u}`$,
where $`0\le\delta_u\lt\pi/2`$. For
$`0\lt\tau\lt\sqrt2`$, put $`\alpha=2\arcsin(\tau/2)`$.
Then

```math
\left\|\bigotimes_u V_u-I\right\|\le\tau
\quad\Longleftrightarrow\quad
\sum_u\delta_u\le\alpha.                               \qquad(7)
```

For necessity, all signed sums of the delta values must lie modulo
$`2\pi`$ in $`[-\alpha,\alpha]`$. Flipping one sign first
shows $`\delta_u\le\alpha`$: the two phases are separated by
$`2\delta_u\lt\pi`$, whereas two points in the allowed arc have
circular distance at most $`2\alpha`$. Now move monotonically
from minus the sum to plus the sum by flipping the signs one at a time.
Each increment is at most $`2\alpha`$, smaller than the forbidden
gap $`2\pi-2\alpha`$. The path cannot move between distinct
allowed intervals on the real line. Its symmetric endpoints must both
belong to the interval centered at zero, proving (7). Sufficiency
follows by bounding every signed sum by the positive sum.

In particular, for all Y bits equal to one, the row angle depends only
on the Hamming distance w. In the stated convergence regime (7) requires

```math
\sum_{w=0}^p\binom pw\delta_{w+1/2}\le\alpha,
\qquad \delta_{3/2}\le\frac\alpha p.                    \qquad(8)
```

Thus even the p distance-one rows require improving accuracy as p
grows. This is an exact condition for independent two-level row words.
It does not impose a lower bound on time t or native depth: special
cancellation times, better row schedules and collective native
synthesis are not excluded. A common gauge near minus identity is
also outside the local $`\delta_u\lt\pi/2`$ convention unless it
is explicitly moved into the target scalar.

The simple sufficient estimate from (3), with no cancellation claim, is

```math
\left\|\bigotimes_u E_{a_u}-(-1)^{Y_x}I\right\|
\le\frac{2Cm}{t}.                                      \qquad(9)
```

The direct-sum structure in x,Y makes this a full-operator estimate,
not an average over addresses or signal inputs.

## 4. A paid native digital recurrence

The following deliberately elementary implementation exposes the cost;
it is an upper bound for this construction, not a lower bound on a
better compiler. Fix $`0\lt\tau\lt1`$ and choose

```math
\begin{gathered}
t=\frac{6Cm}{\tau},\\
J=\left\lceil\frac{3m\pi^2}{\tau}
       \left[t+2t^2(p+\tfrac52)\right]\right\rceil,
\qquad h=\frac\pi J,\\
M=2Jm(p+3),\qquad \epsilon=\frac{\tau}{3M}.
\end{gathered}                                                   \qquad(10)
```

On each half loop use J left-endpoint intervals and the first-order
split into the X part followed in matrix order by the Z part, or the
reverse fixed order. Freezing the Hamiltonian incurs at most
$`t\pi^2/(2J)`$ error per half and row since
$`\|H_a'\|=1`$. The split incurs at most
$`t^2\pi^2(a+1)/J`$. Across both halves and m rows the error is
therefore at most

```math
\frac{m\pi^2}{J}\left[t+2t^2(p+\tfrac52)\right]
\le\frac\tau3.                                         \qquad(11)
```

The first term in (10) makes the adiabatic error (9) at most tau/3.
The diagonal Hamiltonian of each time step is literally

```math
\begin{aligned}
H_Z(\theta)
={}&\sum_u(\tfrac p2+1+\cos\theta)Z_{s_u}
 +\tfrac12\sum_u Z_{Y_u}Z_{s_u}\\
&-\tfrac12\sum_{j=1}^p\sum_u(-1)^{u_j}Z_{x_j}Z_{s_u},\\
H_X(\theta)={}&\sin\theta\sum_u X_{s_u}.
\end{aligned}                                                   \qquad(12)
```

There are p+3 elementary Pauli rotations per row and time step.
Each two-body Z rotation is a CNOT, one single-qubit Z rotation on its
signal, and the actual inverse CNOT. For each fixed address bit j all
m such rotations can use distinct signal T targets in parallel;
their shared address occurs only in the intervening Clifford blocks.
The local X, local Z and bank-signal rounds add three more rotation
rounds. This is a legal schedule without initialized copies.

Let $`t_\epsilon,g_\epsilon,d_\epsilon`$ be common forward and
actual-inverse native bounds for a literal single-qubit rotation to
operator error epsilon, with no work. Such complete SU(2) synthesis,
including its scalar, is the usual ancilla-free rotation problem;
see [Selinger, arXiv:1212.6253v2, Sections 2 and 8](https://arxiv.org/pdf/1212.6253).
The recurrence does not assume that synthesis depth is free. It gives

```math
\begin{aligned}
T_{\mathrm{native}}&\le M t_\epsilon,\\
G_{\mathrm{native}}&\le M g_\epsilon+O(Jm(p+1)+m),\\
D_T&\le2J(p+3)d_\epsilon,\\
b_{\mathrm{dirty}}&=m,\qquad a_{\mathrm{clean}}=0.
\end{aligned}                                                   \qquad(13)
```

The O(m) term includes both literal Y echoes and the bank Hadamard
endpoints. Native approximation of all M rotations adds at most
$`M\epsilon=\tau/3`$ to the full operator norm. The three errors
therefore sum to at most tau, including every signal return error.
The actual reversed native word is the inverse; its norm error is
identical and its separately expanded Clifford count is included in
the common rotation bound. Tensoring with arbitrary reference systems
does not change any operator norm used above.

At fixed tau this elementary schedule has
$`t=O_\tau(m)`$, $`J=O_\tau(m^3(p+1))`$ and
$`M=O_\tau(m^4(p+1)^2)`$ before charging native synthesis.
The unmerged implementation actually emits $`2Jmp`$ address-signal
rotations, each with absolute angle $`\lambda=t\pi/(2J)`$.
The choices above give $`0\lt\lambda\lt\pi/8`$ and

```math
\frac\lambda\epsilon=\frac{3\pi tm(p+3)}\tau.
```

In particular epsilon is smaller than $`2\sin(\lambda/2)`$, the
distance to the nearest literal single-qubit Clifford, which here is
identity. Every separately synthesized address rotation therefore
requires at least one T gate. This explicit unmerged emitter exceeds
the target linear T-count even before its other rotations are charged.
Its large upper depth ledger does not prove a lower bound on the depth
of a different schedule. These statements concern the displayed
digital word; no optimality or limit on joint resynthesis is inferred. Higher-order
simulation or a better rowwise convergence bound might improve them
without supplying the desired simultaneous frontier.

## 5. What changed and the next mechanism question

The uncertainty resolved here is whether a dirty two-level signal can
receive an equality sign on its entire space through a noncommuting
loop. It can: equations (2)--(3) exhibit it without an initialized
eigenstate or a low-rank reflected signal projector. Consequently the
earlier fixed-subspace obstruction should not be extended to these
loops.

The first unpaid native step in the geometric proposal is replacing
the continuous transport by a jointly cheap complete bank word.
Independent row accuracy has the multiplicity charge (7), while the
simple diagonal interaction decomposition (12) carries p address
couplings for each of m signals at every time step. The explicit
recurrence pays both and retains a much larger cost than R91.

The justified next question is whether a collective native word can
implement the winding scalar while correlating or cancelling the
signal errors across the bank, and compressing the repeated address
couplings, on the entire arbitrary input space. Such a word may
temporarily mix the bank and address and need not factor into row
signals or preserve the Hamiltonian form (12). Success must supply
its own full-operator estimate and actual-inverse ledger; it cannot
inherit the rowwise adiabatic limit at no cost. This chapter supplies
neither that word nor a lower bound against it.

The most direct common-signal replacement already loses the selected
information. Identifying all signal Paulis in (12) and summing the
Hamiltonians gives

```math
H_{\mathrm{common}}(\theta)
=m\sin\theta\,X+
 \left[m(\tfrac p2+1+\cos\theta)
       +\tfrac12\sum_u Z_{Y_u}\right]Z.
```

For every address coordinate j,
$`\sum_u(-1)^{u_j}=0`$, so no address operator remains. Its echo
cannot distinguish two addresses selecting different bank values.
This excludes only that literal sum on a shared signal; interleaved
address rotations, nonlinear bank access and joint noncommuting words
leave the premise.

The [focused fixtures](../../tests/test_dirty_geometric_phase_loops.py)
check the literal native central minus sign, every physical column and
actual inverse of a small Clifford-angle digital word, the Hamiltonian
decomposition, and the tensor-error condition including its excluded
near-minus-identity regime. The coarse digital fixture is explicitly
not an accurate indicator. The uniform finite-time bound and all-size
resource recurrence are analytic results.
