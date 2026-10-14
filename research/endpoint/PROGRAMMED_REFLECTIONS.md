# Programmed reflections: coverage, selection, and native cost

[Center-phase family](CENTERED_PHASE_REFLECTIONS.md) · [Paid encodings](PAID_LOGICAL_ENCODINGS.md) · [Diagonal compiler](BOUNDED_DIAGONAL_FACTORIZATION.md) · [Endpoint](../../docs/OPEN_PROBLEM.md)

This chapter separates phase-orbit membership, selection of reversing
reflections, and native implementation. It gives exact programming tests,
a certified rejected pair of complex mixers, a regular spectral chart,
and a query restriction for separately reconstructing its chosen factor.
A cheap alternate baseline is kept beside that restriction. These interfaces
do not themselves supply an unrestricted linear-T compiler. The separate
[growing-source compiler](GROWING_SOURCE_COMPILER.md) resolves the complete-frame
endpoint by a different construction.

Throughout, the original target has N=2^n logical modes, n at least three,
error eta=2^(-N), two initialized flags, and N+n+7 arbitrary dirty qubits.
Literal phases, every logical column, clean leakage, dirty return, and
external references belong to the initialized-isometry norm. Only T-count
is optimized. Classical spectral calculations are not supplied quantum gates.

## 1. One pair: reversibility and exact programming

Let J be a specified Hermitian involution, C an explicitly priced unitary,
and D(d) a diagonal with unit-modulus entries. Define

```math
B=C^\dagger JC,\qquad A(d)=CD(d)C^\dagger,
\qquad V(d)=A(d)JA(d)^\dagger J.
```

For every A, the identity $`JVJ=V^\dagger`$ holds. Conversely,
V is such a pair precisely when VJ is a Hermitian involution with J's
signature. For a diagonal sign center, this forces
$`|V_{xy}|=|V_{yx}|`$, independently of the center rank or the cost of A.

The legal ordered edges (0,N/2), then (N/2,3N/4), at angles alpha,beta,
have one directed entry equal to sin(alpha)sin(beta), with reverse entry
zero. Therefore every diagonal-center pair satisfies

```math
\|W-V\|\ge\tfrac12|\sin\alpha\sin\beta|.
```

If a physical diagonal center Z preserves an initialized embedding Gamma,
write $`Z\Gamma=\Gamma z`$. The compressed block of A Z A-dagger Z obeys
$`z(\Gamma^\dagger V\Gamma)z=(\Gamma^\dagger V\Gamma)^\dagger`$.
The same lower bound applies on any computational dirty input, even with
approximate return. A center swapping the initialized code is outside this
argument: A=diag(W,I) and a flag X give diag(W,W-dagger), but A then contains
the unpriced target. A conjugating exterior operation commuting with J
preserves reversibility; arbitrary extra multiplication need not do so.

### Exact nearest reversible matrix

For a unitary W with $`\|W-I\|\lt1`$, put

```math
Q=WJ,\qquad B_Q=(Q+Q^\dagger)/2,
\qquad\omega=\|(Q-Q^\dagger)/(2i)\|.
```

Let H_* act as +1 on the positive eigenspace of B_Q and as -1 on its
negative eigenspace. The Hermitian part B_Q is within distance one of J
and has the same inertia, so there is no zero eigenspace.
For each unit eigenvector of Q, the nearest possible eigenvalue of an
involution is one of -1 and 1. H_* makes those choices simultaneously.
Thus the unrestricted matrix optimum is

```math
\inf_A\|AJA^\dagger J-W\|
=\omega\sqrt{\frac{2}{1+\sqrt{1-\omega^2}}},
\qquad V_*=H_*J.
```

This prices neither the matrix sign nor its eigenbasis. The earlier
[rank-one center optimum](CENTERED_PHASE_REFLECTIONS.md) has its own,
more restrictive programming geometry.

### Phase-orbit criterion

For the specified C,J, define $`K_W=C^\dagger WJC`$. Then

```math
W=V(d)\quad\Longleftrightarrow\quad K_W=D(d)BD(d)^\dagger.
```

Exact membership requires the same diagonal, zero pattern, and entrywise
magnitudes, with K_W Hermitian. On the graph of nonzero off-diagonal entries
of B, the gains $`g_{ij}=(K_W)_{ij}/B_{ij}`$ must multiply to one around
every oriented cycle. These conditions suffice: set one phase per connected
component, propagate along a spanning tree, and use cycle consistency on
remaining edges. Isolated coordinates have arbitrary phases.

This is the unit-complex graph-gauge argument; see
[Cavaleri and Donno](https://doi.org/10.1515/spma-2022-0169).
For each cycle of length p there is also a division-free bound

```math
\|W-V(d)\|\ge\frac1p\left|
\prod_{r=1}^p(K_W)_{i_r i_{r+1}}-
\prod_{r=1}^p B_{i_r i_{r+1}}\right|,\qquad i_{p+1}=i_1.
```

Each entry of a unitary has magnitude at most one; a product telescope
proves the bound. Diagonal and entrywise-magnitude discrepancies are
separate lower bounds. Matching magnitudes alone is insufficient.
For non-Hermitian comparisons in later sections, opposite directed gains
must be checked separately rather than inferred by conjugation.

## 2. Two programs sharing a basis

Write $`R(u)=CD_uBD_u^\dagger C^\dagger`$, and seek W=R(v)R(u).
Set $`M=C^\dagger WC`$ and $`d_i=\overline v_i u_i`$. Cancellation gives

```math
M=D_v(BD_dBD_d^\dagger)D_v^\dagger.
```

Let $`P_{ij}=|B_{ij}|^2`$ and $`\mu_i=M_{ii}`$. P is real, symmetric,
nonnegative, and doubly stochastic. Taking the diagonal yields

```math
\boxed{(P-\mathrm{diag}(\mu))d=0,\qquad |d_i|=1.}
```

In a one-dimensional kernel the normalized vector must have every modulus
1/sqrt(N). A singular matrix without such a vector does not give a program.
For a valid d, apply the complete directed gauge test to
$`M`$ and $`BD_dBD_d^\dagger`$, recovering v and then u=v*d.
These steps are necessary and sufficient. Higher-dimensional kernels leave
a genuine unit-modulus feasibility problem.

For every candidate pair, its relative phase vector has norm sqrt(N).
Applying the target kernel to that vector is bounded by sqrt(N) times the
operator error. Therefore

```math
\inf_{u,v}\|W-R(v)R(u)\|
\ge\sigma_{\min}(P-\mathrm{diag}(\mathrm{diag}(C^\dagger WC))).
```

A floating singular value is a diagnostic, not a certified endpoint bound.

### The shared real-basis restriction

Suppose C is real up to one scalar, W is real orthogonal, and the support
graph of P is connected. Put $`g=1-\lambda_2(P)\gt0`$ and
$`a_i=1-M_{ii}=\|(M-I)e_i\|^2/2`$. The kernel is
$`P-I+\mathrm{diag}(a)`$. Its largest eigenvalue is at least the mean
of a; all others are at most $`-g+\max a_i`$. If
$`\|W-I\|^2\le g`$, then $`\max a_i\le g/2`$, and

```math
\boxed{\inf_{u,v}\|W-R(v)R(u)\|
\ge\frac{\|W-I\|_F^2}{2N}.}
```

For a disconnected P the family preserves its components; apply this
argument within them. Identity is isolated in the real image of each fixed
real-basis family. This is not uniform over target-dependent bases with
vanishing gaps. Genuinely complex bases, different orbits, and additional
factors are not excluded.

For example, on each bit choose $`C_1=(SH)T^\dagger(SH)^\dagger`$.
It is real up to scalar, $`C_1^\dagger ZC_1=H`$, and its tensor power
has T-count n, with $`P_{ij}=1/N`$ and g=1. The unequal-direction
three-mode path from Section 1, at equal angle theta, has floor
$`(1-\cos\theta)(3+\cos\theta)/N`$. With
$`s=\lceil\sqrt N\rceil+6`$, $`u=2^{-s}/32`$, and
$`\theta=2\arctan u`$, the floor exceeds $`2\eta`$ for n at least eight.
The target itself is an inexpensive special-unitary block.

Complex mixing is a substantive exception. If $`C=Q(SH)_0`$ with Q a
real native circuit, phase pairs $`d_{x\oplus1}=\overline d_x`$ make
$`CD_dC^\dagger`$ real. Full-rank real examples do not, however, establish
coverage of the prescribed family. The target tests below are stronger
than producing matrices already known to lie in a chosen orbit.

## 3. Use the target's reversing equations

For Hermitian involutions, W=R_2 R_1 is equivalent to

```math
H=H^\dagger,\quad H^2=I,\quad WH=HW^\dagger,
\qquad R_1=H,\quad R_2=WH,
```

with the required signatures. For a real W with distinct nonreal conjugate
pairs and no eigenvalues 1 or -1, a real rotation-plane basis gives

```math
W=Q\left(\bigoplus_jR(\theta_j)\right)Q^T,
\qquad
H=Q\left(\bigoplus_j
\begin{pmatrix}a_j&b_j\\b_j&-a_j\end{pmatrix}\right)Q^T,
\quad a_j^2+b_j^2=1.
```

The reversing equation preserves each plane. Anticommutation with its real
skew generator forces even a complex Hermitian candidate to have the real
form displayed. Degenerate eigenspaces have additional freedom and must not
be forced into this parameterization.

For different mixers C_1,C_2 and common center J, the necessary moments are

```math
\mathrm{diag}(C_1^\dagger HC_1)=\mathrm{diag}(C_1^\dagger JC_1),
\qquad
\mathrm{diag}(C_2^\dagger WHC_2)=\mathrm{diag}(C_2^\dagger JC_2).
```

They are linear in the plane coordinates. Involutivity and both complete
phase-orbit tests remain necessary. Without a spectral promise, use all
N^2 real coordinates of Hermitian H. Record real strict-upper and imaginary
upper entries of WH-HW^T, then both real moment lists. This defines a real
system A_W h=b with N^2+2N rows. For a candidate product V within epsilon
of real W, its H=R_1 satisfies reversing residual at most 2epsilon, exact
first moments, and second-moment residual at most epsilon.

### One exact, robust rejected pair

Let Q_1 be the chronological little-endian circuit

```text
H(1), CCX(2,1;0), CCX(0,1;2), H(1), H(2),
CCX(1,2;0), H(0), CX(1;2), CCX(0,2;1), H(2).
```

Let Q_2 reverse this list and relabel j as 2-j. Set
$`C_1=Q_1(SH)_0`$, $`C_2=Q_2(SH)_1`$, J=Z_0. Each fully expanded
mixer uses 68 gates, including 28 T or inverse-T gates, and no helpers.
For the true three-qubit Hopf target with all angles pi/2,

```math
W_*e_0=e_7,\qquad W_*e_j=-e_{j-1}\quad(1\le j\le7).
```

Its characteristic polynomial is lambda^8+1. Exact elimination of the
80-by-64 Hermitian system has rank 62. Every solution's first column is

```math
He_0=(4/3-x,\ -5x/3-y/3,\ 2x/3+y/3-1,\ 2x-2,
       -x,\ -y,\ -1/3,\ -x-2/3)^T.
```

Its squared norm minus one satisfies

```math
\frac{92x^2+14xy-96x+11y^2-6y+57}{9}
\ge\frac{3415}{963}.
```

Thus the linear equations already preclude an involution.

For a quantitative certificate put

```math
q=(259,-275,-215,-304,-169,20,-107,-383)^T/321,
\qquad Q_0=\|q\|^2=4378/963.
```

The [rational certificate](../../tests/fixtures/programmed_reflection_dual.json)
gives lambda with $`A_{W_*}^T\lambda=\ell_q`$ and
$`\lambda^Tb=Q_0`$, where ell_q is the real functional q^T H e_0.
Its error-weighted absolute sum is
$`2\|\lambda_D\|_1+\|\lambda_2\|_1=94889/321`$;
first-moment entries incur no error. Since $`\|q\|\lt15/7`$, every
candidate satisfies

```math
\boxed{\|R_2R_1-W_*\|\ge16201/1992669\gt1/128=2\eta.}
```

The exact test rebuilds all Hermitian coordinates and checks the dual
identities; no phase optimizer enters the certificate. A signed decrement
implements W_* with 14 T gates and no helpers. The restriction belongs to
these two mixers, not to the target. Extending them and W_* by spectator
identities preserves the gap for every n at least three, while the same
14-T target circuit still suffices. Arbitrary larger mixers and physical
embeddings remain outside this statement.

## 4. Baseline selection: a restriction and a paid escape

For even N, every SO(N) matrix is a product of two balanced real reflections.
This is the [Lim--Lu--Ye theorem](https://arxiv.org/abs/2501.18172), not a
new native compiler. A real rotation-plane decomposition also proves this
case directly: each plane rotation is the product of two trace-zero
reflections, and the 1 and -1 eigenvalues can be paired within their signs.
The conjugating spectral basis has not been priced.

There is nevertheless no continuous choice of two logical Hermitian
involutions R_1(theta),R_2(theta), continuous at the zero tuple, whose product
is every nearby complete Hopf frame. At zero their common baseline is R_0.
The reversing equation on each single-edge curve and continuity imply
$`K_eR_0+R_0K_e=0`$. On a connected tree with at least three vertices these
equations force R_0 to be diagonal with opposite values on adjacent vertices:

```math
R_0=\pm J_{\rm tree},\qquad
(J_{\rm tree})_{xx}=(-1)^{\mathrm{popcount}(x)}.
```

Activate three root edges at small nonzero angles. On the four active modes,
$`\det(I-W)=2\prod_{j=1}^3(1-\cos\theta_j)\gt0`$, and there is no
-1 eigenvalue near identity. A reversing reflection preserves those modes
and has signature (2,2) there. The forced baseline has signature (1,3) or
(3,1). Their operator distance is two, contradicting continuity. This is an
exact two-logical-reflection selection theorem, not a finite-error endpoint
bound. More factors, ancilla encodings, or discontinuous selection may avoid it.

The following fixed offset avoids the cancelling-baseline premise entirely.
Its implementation and its subsequent cost limitation are both part of the
statement; neither may be omitted.

## 5. A regular signed-cycle chart

Let M=2^m and define $`Pe_x=e_{x+1}`$ except $`Pe_{M-1}=-e_0`$.
Let $`Se_x=e_{M-1-x}`$. Then

```math
P^M=-I,\qquad SPS=P^T,\qquad P=(PS)S,
\qquad g=2\sin(\pi/M)\ge4/M.
```

S and PS are balanced real reflections. P's distinct eigenvalues are
$`\lambda_j=e^{i(2j+1)\pi/M}`$; let P_j be their projectors.
For E in SO(M), assume

```math
\delta=\|E-I\|\le g/(16\sqrt M),\qquad T=EP.
```

Use fixed circles of radius g/3 around the lambda_j, and let Q_j(T) be the
Riesz projectors of T on those circles. The resolvent norms are at most
3/g and 4/g respectively. The resolvent identity gives
$`\|Q_j-P_j\|\le4\delta/g`$. Define

```math
A=\sum_jQ_jP_j,\qquad
\alpha=4\sqrt M\delta/g\le1/4,\qquad
O=A(A^\dagger A)^{-1/2}.
```

Orthogonality of the P_j decomposition gives $`\|A-I\|\le\alpha`$.
The polar factor transports each P_j to Q_j. Conjugate pairing makes A and
O real, and $`\|O-I\|\le2\alpha`$. Consequently

```math
H_1(E)=OSO^T,\qquad H_2(E)=TH_1(E),\qquad
\boxed{E=H_2(E)H_1(E)P^T.}
```

Both factors are balanced real reflections, selected real-analytically on
the stated neighborhood. They satisfy

```math
H_1(I)=S,\quad H_2(I)=PS,\quad
\|H_1-S\|\le16\sqrt M\delta/g,\quad
\|H_2-PS\|\le\delta+16\sqrt M\delta/g.
```

The extra P^T is a paid circuit, not a phase convention. This is not the
two-cancelling-factor form excluded in Section 4.

### Certified classical evaluation

For rho=g/3 and L equally spaced unit-circle points zeta_k, use

```math
\widetilde Q_{j,L}=\frac1L\sum_{k=0}^{L-1}
\rho\zeta_k[(\lambda_j+\rho\zeta_k)I-T]^{-1}.
```

Diagonalizing the normal T evaluates this quadrature as a scalar geometric
sum on each eigenline. The inside ratio is at most 3delta/g; outside ratios
are at most 4/11. Hence

```math
\|\widetilde Q_{j,L}-Q_j\|
\le\frac{(4/11)^L}{1-(4/11)^L}.
```

All sample inverses have known norm bounds, and A-dagger A stays uniformly
positive. Certified finite arithmetic can therefore compute factor entries
to any requested accuracy without spectral sorting or an exact-zero oracle.
The floating eigensolver and polar SVD in the diagnostic code do not perform
this certification and do not emit native variable reflections.

### The offset and the original workspace

To emit P, apply a phase at the all-ones label and increment the data in
place. A k-control X uses 4k-8 Toffolis for k at least three, with k-2 arbitrary
chain bits. The echo toggles the target first with and then without the
leading conjunction, canceling every unknown chain term; all helpers return
exactly. Expanding the retained seven-T Toffoli gives, for m at least four,

```math
T(P)=7(2m^2-6m+1),\qquad b_P=m-3,\qquad a_P=0.
```

The m=3 count is 14 T gates. The four-data-bit word has 142 elementary gates,
63 T gates, and one arbitrary helper. S is bitwise complement. Actual
inverses preserve the signed wrap and literal phase.

The original coarse word may be complex. Reserve one of the allowed dirty
bits r, use b_0=N+n+6 remaining helpers, and set
$`s=\lceil\sqrt N\rceil+6`$. The
[borrowed coarse theorem](../../docs/BORROWED_WORKSPACE_COMPILER.md) gives
an exactly returned logical Q with error below 2^(-s) and O(N) T-count.
[Native realification](NATIVE_REALIFICATION.md) uses the arbitrary bit r,
at most doubles T-count, and gives a real Q_R on M=2N modes. Put

```math
E=(I_r\otimes W)Q_R^T,\qquad
I_r\otimes W=H_2(E)H_1(E)P^TQ_R.
```

This E lies within the chart because $`s\ge2+3(n+1)/2`$. The allocation is
$`(n+1)+2+(N+n+6)=N+2n+9`$. The bit r is unknown data; its final identity
action returns all its reference correlations. The coarse error belongs to
the exact residual and is not added again. Only Q_R and the fixed P are
natively priced here. The two selected fine reflections remain uncompiled.

## 6. Quantum reconstruction of the selected factor can be harder

Use M=2N and the actual last-depth-only Hopf family

```math
a_p=-\sin\frac{2\pi(2p+1)}N,\qquad
W_t=\bigoplus_{p=0}^{N/2-1}R_y(ta_p),\qquad
E_t=I_2\otimes W_t=e^{tK},\qquad\|K\|\le1.
```

This is a genuine residual family with the common coarse choice Q=I; no
claim is made that every fixed coarse algorithm selects identity. Set
$`T_t=E_tP`$, and $`S`$ as above. For any matrix H define
$`h=\|H-S\|`$ and $`\rho_t(H)=\|T_tH-HT_t^T\|`$. Put

```math
d_N=2\sin(\pi/N),\qquad
f_N=\cos(\pi/N)\cos(\pi/(2N)).
```

Then

```math
\boxed{(d_N+2|t|)h\ge f_N|t|-t^2-\rho_t(H).}
```

For a direct proof use centered signed-cycle eigenvectors

```math
f_j(x)=M^{-1/2}\exp[-i\omega_j(x-(M-1)/2)],\qquad
\omega_j=(2j+1)\pi/M.
```

They satisfy P f_j=exp(i omega_j)f_j and S f_j=conjugate(f_j).
For j=N/2-1, k=N/2+1, summing the matching entries gives

```math
f_j^\dagger Kf_k=-\tfrac12\cos(\pi/(2N)),\qquad
|f_j^\dagger(KPS+PSK)\overline f_k|=f_N.
```

The unitary Taylor remainder in E_t is at most t^2/2. Thus the corresponding
entry of T_t S-S T_t^T has magnitude at least f_N|t|-t^2. Write H=S+X.
In that entry, P X-X P^T has coefficient of magnitude d_N; replacing P by
T_t adds at most 2|t| times the norm of X. This proves the boxed bound.
It assumes neither differentiability nor the particular polar choice.

For an exact reverser and $`0\lt t\le\min(1/4,d_N/2)`$,

```math
\|H-S\|\ge t/(4d_N)\ge Nt/32.
```

Now consider a fixed-query-skeleton reconstruction with Q calls to T_t,
its adjoint, or controlled versions. All other gates and encodings may depend
on N but not t. Arbitrary clean/dirty widths and references are allowed.
Each query changes by at most |t| from t=0. If the requested factor has
H_0=S and errors epsilon_t,epsilon_0, the full-unitary hybrid gives

```math
\|H_t-S\|\le Q|t|+\epsilon_t+\epsilon_0.
```

At $`t_*=2^{-s}/32`$, $`s=\lceil\sqrt N\rceil+6`$, and both errors
at most eta, the preceding inequality implies

```math
\boxed{Q\ge N/64\quad(n\ge8).}
```

Indeed $`2\eta/t_*=2^{s+6-N}\le N/64`$. This is a query bound for an
individually reconstructed S-anchored factor. It is not a bound on classical
factor evaluation, direct synthesis with angle-dependent gates, or the joint
product without separately realizing its factors. Multiplying it by an
independent module price does not establish an unrestricted T lower bound.
The hybrid mechanism is the standard one from
[Bennett et al.](https://arxiv.org/abs/quant-ph/9701001); the forcing family
and finite inequality are proved above.

### A zero-T baseline avoids that penalty for the same family

Set $`J_*=P^{M/4}S`$. The quarter shift acts only on the two highest data
bits: chronologically CZ, CNOT from the lower to the higher, then X on the
lower. Together with bitwise complement this is an exact Clifford word with
no helpers. It obeys

```math
J_*PJ_*=P^T,\qquad J_*KJ_*=-K.
```

Hence $`\widetilde H_1=P^TJ_*`$ is constant and
$`\widetilde H_2(t)=E_tJ_*`$ gives

```math
E_t=\widetilde H_2(t)\widetilde H_1P^T,\qquad
\|\widetilde H_2(t)-J_*\|=\|E_t-I\|\le|t|.
```

Both factors are balanced real reflections. The first baseline is not S:
its distance is $`2\cos(\pi/(2M))`$. The same construction holds for
arbitrary final-layer angles satisfying
$`\theta_p=\theta_{(N/4-1-p)\bmod(N/2)}`$, leaving N/4 independently
variable, possibly signed or singular angles. One packed diagonal compiles
the complete matching target with O(N) T-count. Fixed P and J_* are already
priced, and may cancel before synthesis. This is an escape for a previously
inexpensive family, not generic branching coverage.

## 7. Native ledgers and evidence boundaries

Let exact mixers and centers act on the stated logical wires with their own
O(N) T counts. The packed diagonal theorem uses
$`q=\lceil(\ell+4)/2\rceil`$ dirty core bits, q+n+1 dirty helpers in
total, two clean flags, and at most $`1344N+36q-18`$ T gates per occurrence.
The following are sufficient ledgers only when the complete orbit tests pass.

| Representation | Fine slots | Precision exponent | Additional fine-slot T-count |
|---|---:|---|---|
| One programmed pair | 2 | N+3 | Twice the per-slot bound |
| Two reflections, shared mixer C | 3 after cancellation | N+4 | 4086N+378 |
| Two reflections, different C_1,C_2 | 4 | N+4 | 5448N+504 |

The shared word adds $`6T(C)+2T(J)`$; the different-basis word adds
$`4[T(C_1)+T(C_2)]+2T(J)`$. Their native diagonal errors are at most
eta/4, with another strict eta/4 available for certified phase selection.
Finite rational-circle search under exact coverage avoids exact real
support or rank decisions. The imported initialized-isometry telescope
propagates actual leakage through later unitary factors; no reset occurs.
On M=2N unknown data modes, four slots at exponent N+4 use at most
$`N/2+n+6\le b_0`$ dirty helpers and 10824N+504 T gates, before pricing
all mixers. This is a width check, not coverage of the spectral reflections.

[programmed_reflections.py](../../compiler_robust_hopf/programmed_reflections.py)
contains bounded numerical orbit, kernel, and moment diagnostics. Every
floating support/rank decision is labeled as such. It emits no arbitrary-angle
program or certified eigenbasis.
[spectral_reflection_chart.py](../../compiler_robust_hopf/spectral_reflection_chart.py)
reuses the existing native expansion, Toffoli, and inverse conventions. It
regenerates the two diagnostic mixers, the 14-T corner target, the 63-T
four-bit signed cycle, and the zero-T alternate baseline; gate lists are not
stored redundantly.

The [orbit tests](../../tests/test_programmed_reflections.py) reconstruct the
exact rational dual using the retained coefficient ring, check full native
corner/inverse matrices, and distinguish simple-kernel recovery from a
complete programming proof. The [chart tests](../../tests/test_spectral_reflection_chart.py)
check full native dirty-return and inverse action, the Fourier forcing
coefficient, the alternate baseline, quadrature error, and the physical width
and query inequalities. Finite matrices are diagnostics, not endpoint-precision
simulations. These routines emit no selected variable spectral reflection or
general complete-frame linear-T word. Their results constrain and implement
the specified interfaces; the separate growing-source construction does not
remove the unpaid native steps of this reflection route.
