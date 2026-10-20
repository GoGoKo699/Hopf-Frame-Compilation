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

### 5.1. A complete distance certificate and the native scheduling gap

For a replacement source on $`b\ge1`$ address bits, put $`m=2^b`$
and retain the actual-adjoint echo $`E_d=ZU_d^\dagger ZU_d`$.
Define its full dirty-signal errors by

```math
e_0=\|E_0+I\|,\qquad e_d=\|E_d-I\|\quad(1\le d\le b+1).
```

Suppose the completed row words preserve x,Y and act on disjoint
arbitrary signals, with $`d=|x\oplus u|+1-Y_u`$. A sufficient
complete-bank certificate is

```math
\max(e_0,e_1)+
\sum_{w=1}^b\binom bw\max(e_w,e_{w+1})\le\delta.
```

The distance-zero row gives the first term. For each positive distance
w there are exactly $`\binom bw`$ rows, and either bank-bit value
selects one of the two displayed errors. Tensor-product telescoping
and the preserved-sector direct sum prove full operator error at most
delta, including every dirty-signal input, its return and references.
This applies to the actual emitted echoes, so their synthesis errors
must be included in the e_d. The pulse in Section 3 supplies one analog
profile; constant accuracy for one isolated row does not supply this
certificate.

There is a separate simultaneous-scheduling obligation. Within one row,
the axes $`Z_s,Z_YZ_s,Z_{x_j}Z_s`$ are independent, so a Clifford
coordinate change can place their rotations on distinct wires. Some
of those wires are shared address wires. This does not give simultaneous
depth for all rows.

In fact the unchanged address--signal rotation list has an exact packing
limit. Its mb binary Pauli labels are the edge-incidence vectors of the
connected complete bipartite graph $`K_{b,m}`$. Their rank over
$`\mathbb F_2`$ is $`m+b-1`$: every edge has even total parity,
and a spanning tree gives that many independent edges. A batch obtained
by one Clifford coordinate change to distinct rotation targets must be
an independent subset. Any partition into such batches therefore needs

```math
\left\lceil\frac{mb}{m+b-1}\right\rceil=b.
```

The equality follows from $`m=2^b\gt(b-1)^2`$ for every positive
integer b. Conversely, fixing one address index gives a star of m
independent axes, so the b address rounds in Section 5 attain this
bound. Adding idle dirty helpers and changing Clifford coordinates do
not change the rank of the unchanged axis list.

This is an exact limit on packing these intact rotations, not a lower
bound on unrestricted native T-depth. Each rotation batch still needs
its native synthesis ledger; special angles, cross-row resynthesis,
new helper-dependent axes or a different word can leave the packing
premise. The first unpaid operation in this proposed improvement is a
complete native source with the distance certificate above and a valid
sublogarithmic simultaneous schedule. A coordinate improvement for one
row or an analog pulse width does not provide it.

### 5.2. Joint resynthesis pays one special-angle detuning call

The intact-axis limit is not a depth bound even for the same diagonal
operator. For $`b\ge1`$, $`m=2^b`$, the whole-bank call

```math
R=\prod_u\exp\left(\frac{i\pi}{4}
       (|x\oplus u|+1-Y_u)Z_{s_u}\right)
```

has an exact native implementation with no additional helpers, no clean
work, and common forward/actual-inverse bounds

```math
D_T\le11,\qquad T\le5m+3b+5,\qquad
G=O(m^2+b^2+mb).
```

All x, Y and signal inputs are arbitrary. The statement retains the
literal scalar and returns every computational coordinate. This is a
paid detuning component, not a selective inversion or complete indicator.

Put $`\omega=e^{i\pi/4}`$ and
$`T_P=\exp(i\pi(1-P)/8)`$. Expanding the detuning gives exactly

```math
R=\omega^{m(b+2)/2}
 \left(\prod_uT_{Z_{s_u}}^{-(b+1)}T_{Z_{Y_u}Z_{s_u}}^{-1}\right)
 \left(\prod_{u,j}T_{Z_{x_j}Z_{s_u}}^{(-1)^{u_j}}\right).
```

The edge signs sum to zero for each address index; this fixes the
displayed scalar. Replacing each negative edge power by
$`T_P^{-1}=T_PS_P^\dagger`$, with $`S_P=T_P^2`$, costs only
Cliffords. Write the integer Hamming weights
$`S_X=\sum_jx_j`$ and $`S_S=\sum_us_u`$. The all-positive
edge product has computational phase

```math
\omega^{\sum_{u,j}(x_j\oplus s_u)}
 =\omega^{mS_X+bS_S-2S_XS_S}.
```

Combining its single-bit powers with the preceding local factors leaves
$`T_{Z_{s_u}}^{-1}`$ and $`T_{Z_{Y_u}Z_{s_u}}^{-1}`$ for
each u, together with $`T_{Z_{x_j}}^m`$, which is Clifford since
m is even. The remaining 2m local axes are jointly independent: compute
$`Y_u\mathbin\oplus s_u`$ into each existing Y wire, apply T-daggers
on all Y and signal wires in one band, then undo those CNOTs.

For either group A, let $`P_A=S_A\bmod2`$ and
$`Q_A=\binom{S_A}{2}\bmod2`$. Since
$`S_A\equiv P_A+2Q_A\pmod4`$, the remaining cross phase obeys

```math
-2S_XS_S\equiv-2P_XP_S+4(Q_XP_S+P_XQ_S)\pmod8.
```

Extract each parity into an existing wire by reversible linear
Cliffords. The first term is a CS-dagger between the parity wires:
T-daggers on both wires and a T on their XOR give two T bands and three
T gates. Each other term is a controlled quadratic sign. A reversible
linear change on its target group brings the quadratic into disjoint
pair products plus a linear form, by elimination of its alternating
bilinear form over $`\mathbb F_2`$. The linear part is controlled-Z
Clifford work. The disjoint pair products use the
[shared-control CCZ batch](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers):
for r pairs, four T bands and $`6r+(r\bmod2)`$ T gates,
with no helper. Each complete quadratic block uncomputes its coordinate
change before the next block. This remains valid if a parity wire in
the target group moves during that block.

The depth is therefore at most $`1+2+4+4=11`$. With at most
$`\lfloor b/2\rfloor`$ and $`\lfloor m/2\rfloor`$ pairs,
the T count is at most
$`2m+3+6\lfloor b/2\rfloor+6\lfloor m/2\rfloor+2`$,
which gives the displayed bound. Linear changes cost
$`O(b^2+m^2)`$ CNOTs; signed-edge Clifford corrections cost
$`O(mb)`$. The literal scalar is emitted using
$`(HS)^3=\omega I`$ with its exponent reduced modulo eight.
Reversing the actual native word supplies the inverse; replacing every
inverse S by three positive S gates preserves the stated common bounds.

There is an exact restriction on using only this angle in independent
rows. A row call $`R_d=\exp(i\pi dZ/4)`$ has period eight.
Thus a completed row word whose only address dependence enters through
these calls at $`d=|x\oplus u|+1-Y_u`$ and their inverses, with
address-independent interlayers on that row,
has $`V_{u,x}=F_u(|x\oplus u|)`$ and
$`F_u(w+8)=F_u(w)`$. Suppose the completed bank has the disjoint,
address-preserving row factorization (1). For $`b\ge8`$, the
Heisenberg image of $`X_{Y_u}`$ is identical at distances zero and
eight. The corresponding target images for $`D_b`$ are opposite.
The same row-isolation argument as in Section 1 therefore gives

```math
\|U-D_b\otimes I\|\le\tau
\quad\Longrightarrow\quad 2\le4\tau,
\qquad \tau\ge\tfrac12.
```

This guards completed independent-row words at the fixed angle. It does
not exclude joint-bank interlayers: the whole vector of detunings modulo
eight can distinguish addresses even when one row aliases. The paid
global call above remains available to such a different mechanism.
No claim about an arbitrary joint word follows from a single-row period.

## 6. What changed and where to continue

Independent-row eigenphase accumulation does not force the earlier
exponentially large analog width. The selective-inversion echo and its
exponential detuning profile prove a logarithmic width is enough, with a
full-bank estimate and a fully paid finite implementation. However,
fixing that analytic precision problem does not fix the resource
frontier. The native independent-row theorem locates a separate linear
factor p in the count, irrespective of pulse design or row compilation.

For a standalone linear-count bank, native work must be shared across
rows or the completed-row factorization must change. The later
[root-fusion reduction](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
allows polynomial address overhead in a smaller source, so the factor
p in (2) alone no longer stops that discovery route. The current digital
pulse still has exponential overhead and excessive depth. A replacement
must supply the complete distance certificate and a genuinely faster
native schedule; changing only the row count does not do so. Section 5.2
pays one whole-bank native call and shows why intact-axis packing must
not be mistaken for a depth obstruction. It does not yet supply the
selective word or its complete distance certificate. The general
complete approximate indicator and optimal general T-depth remain open.
No unrestricted depth lower bound is inferred from analog duration,
the intact-rotation packing limit, a simulation ledger or a row T-count.

The [focused fixtures](../../tests/test_selective_inversion_dirty_indicator.py)
check actual native echo phases and inverses, a finite digital word on
all its physical columns, and the selective-inversion operator identities.
They do not replace the all-size proofs or certify a coarse fixture as
an accurate indicator.
The [joint-detuning fixture](../../tests/test_joint_detuning_resynthesis.py)
checks a literal 11-band native word and its actual inverse on all
physical columns at b equal to two; the all-size claim is the analytic
resynthesis above.
