# Selective inversion, dirty scalar echoes, and independent-row cost

[Geometric precursor](DIRTY_GEOMETRIC_PHASE_LOOPS.md) ·
[Complete approximate transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Joint routing interfaces](JOINT_ROUTING_INDICATOR_INTERFACES.md)

**Current scope clarification.** The independent-row count theorem below
still excludes a standalone linear-count bank. The later
[polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
can absorb polynomial overhead in a smaller first bank, so that count
objection alone no longer retires a source at the relaxed interface.
The explicit digital pulse below has exponential overhead and excessive
depth; it remains unsuitable. A replacement would need both the complete
bank error proof and genuinely sublogarithmic native depth.

A selective inversion pulse supplies a dirty scalar sign through an actual
adjoint, without an adiabatic or initialized-eigenstate premise. Its
Hamming-distance profile can control the complete bank with a logarithmic
pulse-width parameter. This improves the earlier uniform analog estimate,
but does not provide a faster native indicator. A separate native theorem
shows why independently compiled rows cannot meet the linear T-count:
every row needs at least p T gates at full-bank error below one half.
The theorem permits temporary address motion inside each completed row;
it does not constrain a jointly synthesized bank word. R96 remains the
strongest unconditional general upper bound.

## 1. Complete target and the independent-row stopping theorem

Let $`p\ge2`$, $`m=2^p`$, and consider the phase form of the
complete indicator,

```math
D_p|x,Y,w\rangle=(-1)^{Y_x}|x,Y,w\rangle,
\qquad I_p=H_YD_pH_Y.
```

All input systems and references are arbitrary. A literal full-operator
estimate is required; no input-dependent scalar is ignored.

The [native sharing theorem](NATIVE_INDICATOR_SHARING.md) proves that
an independently emitted row factorization

```math
U=\prod_{u\in\mathbb F_2^p}U_u,
\qquad U_u=\sum_x|x\rangle\langle x|\otimes V_{u,x}
                                                               \qquad(1)
```

cannot meet the linear T budget. Here $`V_{u,x}`$ is any unitary on
$`Y_u`$ and a private dirty register $`W_u`$, and those private
registers are disjoint. A row can temporarily mix its address, bank and
private work, provided its completed operator has (1). Neither its bank
qubit nor its private work must return at that boundary. If row u's
native word contains $`q_u`$ T gates, then

```math
\|U-D_p\otimes I_W\|\le\tau\lt \tfrac12
\quad\Longrightarrow\quad
q_u\ge p\ \hbox{for every }u,\qquad
\sum_uq_u\ge mp.                                  \qquad(2)
```

There is no separate row-accuracy premise. The source theorem isolates
each completed row through the Heisenberg image of $`X_{Y_u}`$,
then uses the affine Pauli-support dimension of a native word. The
remaining sections improve the analog pulse profile, but their native
emitter has precisely the independent structure covered by (2).
Shared native synthesis across rows, nonlocal bank interactions and
address mixing across incomplete row boundaries leave this premise.
Neither (2) nor the pulse-width bound below is an unrestricted T-depth
lower bound.

## 2. An actual-adjoint echo is exactly selective inversion

Write any $`U\in SU(2)`$ as

```math
U=\begin{pmatrix}a&b\\-\overline b&\overline a\end{pmatrix},
\qquad |a|^2+|b|^2=1.
```

Its actual-adjoint Z echo is

```math
\begin{aligned}
E(U)&=ZU^\dagger ZU\ =XU^{\mathsf T}XU,\\
E(U)&=\begin{pmatrix}
1-2|b|^2&2\overline a b\\
-2a\overline b&1-2|b|^2
\end{pmatrix},\\
\|E(U)-I\|&=2|b|,\qquad
\|E(U)+I\|=2|a|.                                    \qquad(3)
\end{aligned}
```

Direct multiplication proves the displayed matrix. Its eigenvalues are
$`e^{\pm i\delta}`$, where $`\delta=2\arcsin|b|`$, proving
the norm statements. A scalar multiplying U cancels in the echo, so the
same conclusion applies after choosing an SU(2) representative of a
one-qubit unitary. The implemented word always uses the actual adjoint;
there is no need to supply a transpose operation.

Perfect population inversion, $`a=0`$, gives literal minus identity
on both dirty signal states, irrespective of the inversion phase.
Perfect nontransition, $`b=0`$, gives identity, irrespective of the
diagonal dynamical phases. No signal initialization is needed. The
one-qubit signal also has no unused complementary subspace of the kind
in the earlier multi-qubit rank-one reflection obstruction.

## 3. A logarithmic-width selective pulse controls the full bank

The Rosen–Zener scattering formula is the only external analytic pulse
ingredient used here. In the conventions

```math
H_d(s)=dZ+\frac{1}{2\kappa}\,\mathrm{sech}(s/\kappa)X,
\qquad -\infty\lt s\lt \infty,
```

the scattering matrix has transition magnitude

```math
|b_d|=\mathrm{sech}(\pi\kappa d).                     \qquad(4)
```

This is equation (13b) of
[Torosov and Vitanov, *Smooth composite pulses for high-fidelity quantum information processing*, Phys. Rev. A 83, 053420](https://link.aps.org/accepted/10.1103/PhysRevA.83.053420),
with their $`\Omega_0\kappa=1`$ and
$`\Delta=2d`$. The scattering representation removes the free
asymptotic Z evolution. At d equal to zero the pulse area gives
$`U_0=-iX`$, so its echo is exactly minus identity. We use the
verified transition formula, not a conjecture about general pulse shapes.

Borrow one arbitrary signal per output label and put

```math
d_u(x,Y)=|x\oplus u|+1-Y_u\in\{0,1,\ldots,p+1\}.
                                                               \qquad(5)
```

There is a resonant row precisely when $`u=x`$ and $`Y_u=1`$.
Its echo gives the required minus sign exactly. Every other row's full
signal error from identity is $`2\,\mathrm{sech}(\pi\kappa d_u)`$.
For $`q=e^{-\pi\kappa}`$, tensor-product telescoping gives

```math
\left\|\bigotimes_u E(U_{d_u})-(-1)^{Y_x}I\right\|
\le4\left[(1+q)^p-1+q\right].                       \qquad(6)
```

For rows at distance $`w\ge1`$, the smallest detuning is w and
there are $`\binom pw`$ such rows. The distance-zero row is either
exactly resonant or has detuning one. Now use
$`\mathrm{sech}(z)\le2e^{-z}`$ and the binomial theorem.
This proof is uniform in x and every bank string; the preserved-sector
direct sum therefore gives full operator accuracy, including entangled
signals and references.

Fix $`0\lt \tau\lt 1/2`$ and choose

```math
\kappa=\frac1\pi\log\frac{32(p+1)}\tau,
\qquad q=\frac\tau{32(p+1)}.                         \qquad(7)
```

Since $`pq\lt 1`$ and $`e^{pq}-1\le2pq`$, (6) is at most tau/4.
Thus the pulse-width parameter is $`O_\tau(\log(p+2))`$,
improving the geometric precursor's uniform inverse-time estimate. The
bank error is still a sum; the exponentially decaying distance profile
pays that sum more efficiently.

There is also a matching necessary width order within this convergent
pulse family. Put $`\alpha=2\arcsin(\tau/2)`$ and suppose
$`2\arcsin(\mathrm{sech}(\pi\kappa))\lt \pi/2`$. For the
all-one bank, the [exact tensor-arc criterion](DIRTY_GEOMETRIC_PHASE_LOOPS.md#3-complete-bank-and-the-exact-error-sum-requirement)
implies

```math
2p\arcsin(\mathrm{sech}(\pi\kappa))\le\alpha,
\qquad
\kappa\ge\frac1\pi\mathrm{arcosh}
 \left(\frac1{\sin(\alpha/(2p))}\right).              \qquad(8)
```

This is $`\Omega_\tau(\log(p+2))`$ for the analog width
parameter, not a native T-depth lower bound. It assumes the displayed
pulse family and its local convergence branch.

## 4. Finite truncation does not grant free detuning phases

Use scaled time $`v=s/\kappa`$ and truncate to $`[-R,R]`$,
where

```math
R=\log\frac{16m}\tau.
```

The actual finite propagator $`U_{d,R}`$ solves

```math
i\partial_v U=
 \left[\kappa dZ+\tfrac12\mathrm{sech}(v)X\right]U.
```

In the interaction frame of $`\kappa dZ`$, the omitted drive has
integrated norm at most

```math
\int_{|v|>R}\tfrac12\mathrm{sech}(v)\,dv\le2e^{-R}.
```

Duhamel's formula bounds the scattering matrix's distance from
$`F U_{d,R}F`$ by that quantity, where
$`F=e^{i\kappa dRZ}`$. There is no need to compile F:

```math
E(FUF)=F^\dagger E(U)F.
```

The scalar target commutes with F, so its operator distance is unchanged.
Changing U by epsilon changes the actual-adjoint echo by at most twice
epsilon. Across all rows, truncation therefore contributes at most
$`4me^{-R}=\tau/4`$, including signal return error. The physical
pulse interval has length $`2\kappa R=O_\tau(p\log(p+2))`$.
Its length is not identified with a native depth.

## 5. A fully paid finite native word

Split the scaled interval into J equal pieces of length $`h=2R/J`$.
On each piece freeze $`a(v)=\tfrac12\mathrm{sech}(v)`$ at the
left endpoint and use one fixed-order product of the Z and X evolutions.
Since $`|a'|\le1/2`$, freezing contributes at most $`R^2/J`$
per row. The commutator has norm at most $`\kappa(p+1)`$, so
first-order splitting contributes at most
$`2R^2\kappa(p+1)/J`$. Thus one pulse has row error at most

```math
\frac{R^2}J\left[1+2\kappa(p+1)\right].
```

Choose

```math
\begin{gathered}
J=\left\lceil\frac{8mR^2}\tau
       \left[1+2\kappa(p+1)\right]\right\rceil,\\
M=2Jm(p+3),\qquad \epsilon=\frac\tau{4M}.
\end{gathered}                                                   \qquad(9)
```

Both actual-adjoint pulse orientations and all rows then contribute at
most tau/4 discretization error. The diagonal generator at one step is

```math
\kappa\sum_u\left[
 \frac{p+1}{2}Z_{s_u}
 +\frac12 Z_{Y_u}Z_{s_u}
 -\frac12\sum_{j=1}^p(-1)^{u_j}Z_{x_j}Z_{s_u}
\right].                                                       \qquad(10)
```

The X part is $`a(v)\sum_uX_{s_u}`$. There are p+3 Pauli
rotations per row and step. A ZZ rotation uses a CNOT, one Z rotation
on the signal, and the actual inverse CNOT. For a fixed j, all rows have
different signal targets and shared address bits appear only as Clifford
controls; their T layers can run together. The p address rounds, local
Z, bank-signal and local X rounds are all charged. No clean fanout is
used.

Let $`t_\epsilon,g_\epsilon,d_\epsilon`$ be common forward and
actual-inverse native bounds for a literal ancilla-free single-qubit
rotation at operator error epsilon. The finite word, including its
actual inverse, has

```math
\begin{aligned}
T_{\rm native}&\le M t_\epsilon,\\
G_{\rm native}&\le M g_\epsilon+O(Jm(p+1)+m),\\
D_T&\le2J(p+3)d_\epsilon,\\
b_{\rm dirty}&=m,\qquad a_{\rm clean}=0.
\end{aligned}                                                   \qquad(11)
```

The O(m) term includes literal Z echoes and the bank Hadamard endpoints.
All M rotation errors sum to at most tau/4. Equations (6)--(11) therefore
supply full error at most tau for the complete indicator, with the four
parts accounting for selectivity, truncation, discretization and native
rotation synthesis. Reversing the actual emitted native gates supplies
the inverse with the same error. This is an explicit, deliberately
unoptimized construction, not a free analog primitive.

At fixed tau, its parameters satisfy

```math
\kappa=\Theta_\tau(\log(p+2)),\quad
R=\Theta_\tau(p),\quad
J=\Theta_\tau(mp^3\log(p+2)),\quad
M=\Theta_\tau(m^2p^4\log(p+2)).                      \qquad(12)
```

This complete recurrence does not meet the requested native resource
envelope. More decisively, (2) shows that even an optimal independent
native compiler for each completed row still needs at least mp T gates
at the requested accuracy. Improving the simulation formula alone cannot
remove that count obligation. Joint resynthesis of all rows is outside
the independent-emitter theorem.

## 6. What changed and where to continue

Independent-row eigenphase accumulation does not force the earlier
exponentially large analog width. The selective-inversion echo and its
exponential detuning profile prove a logarithmic width is enough, with a
full-bank estimate and a fully paid finite implementation. However,
fixing that analytic precision problem does not fix the resource
frontier. The native independent-row theorem locates a separate linear
factor p in the count, irrespective of pulse design or row compilation.

The next useful mechanism must share native non-Clifford work across
rows, or mix address/bank systems across the boundaries at which (1)
would apply. Merely improving isolated pulse selectivity cannot achieve
$`T=O_\tau(m)`$. The general complete approximate indicator and
optimal general T-depth remain open. No obstruction here applies to the
unrestricted joint model, and no T-depth lower bound is inferred from
analog duration, a first-order simulation ledger, or a row T-count.

The [focused fixtures](../../tests/test_selective_inversion_dirty_indicator.py)
check actual native echo phases and inverses, a finite digital word on
all its physical columns, and the selective-inversion operator identities.
They do not replace the all-size proofs or certify a coarse fixture as
an accurate indicator.
