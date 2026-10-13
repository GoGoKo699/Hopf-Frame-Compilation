# Source-program resynthesis and completed-group gauges

[Source reuse](SOURCE_REUSE_LIMITS.md) · [Grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md) · [Endpoint](../../docs/OPEN_PROBLEM.md)

This chapter distinguishes three tasks: reproduce a raw source program's
complete core action, replace a recovered group that returns its core, and
choose a common physical precision core for the existing endpoint compiler.
The first task can require a paid basis and address connectors. The second
admits independent commuting core gauges without physical basis changes or
connectors. The third removes width transitions at the selected endpoint but
retains the iterated-logarithmic number of precision charges.

None of the results below proves a uniform linear-T complete-frame compiler.
The target remains N=2^n, n at least three, two initialized external flags,
N+n+7 arbitrary dirty qubits, and error eta=2^(-N). Literal phase, all logical
columns, clean leakage and dirty/reference return belong to one initialized-
isometry norm. Only T-count is optimized; classical preprocessing may be
expensive but must terminate. There is no supplied source state, QRAM,
measurement, reset, postselection, or uncharged unitary evaluator.

## 1. The actual scalar word and exact versus approximate resynthesis

Use the [operator source](../../docs/OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit)
on m arbitrary core qubits. Its Hermitian anticommuting generators and weights are

```math
\Gamma_j=Z_0\cdots Z_{j-1}X_j,\qquad
M=\sum_{j=0}^{m-1}a_j\Gamma_j,
```

```math
a_j=2^{-(j+1)/2}\ (j\lt m-1),\qquad
 a_{m-1}=2^{-(m-1)/2}.
```

A mask f gives $`N_f=\sum_j a_j(-1)^{f_j}\Gamma_j`$. The actual one-flag
scalar word, with all rejected components retained, is

```math
\mathcal S_f=c_f I+D_f\otimes X_\sigma,\qquad
c_f=\sum_j a_j^2(-1)^{f_j},\qquad
D_f=(MN_f-N_fM)/2.
```

Here $`D_f^\dagger=-D_f`$ and $`D_f^2=-(1-c_f^2)I`$. These are identities
on the entire core space, not expectation values in a prepared source.
For one repeated mask, the zero-signal block of $`\mathcal S_f^R`$ is
$`T_R(c_f)I`$, not $`c_f^R I`$: already the second and third powers give
$`2c_f^2-1`$ and $`4c_f^3-3c_f`$. Rejected amplitudes can re-enter acceptance.

### Exact repetition can retain an unnecessary denominator cost

For m at least four choose f to flip only position m-2, leaving the final
encoding bit zero. Set

```math
a=2^{-(m-1)/2},\quad \Gamma=\Gamma_{m-2},\quad
c=1-2a^2,\quad \theta=2\arcsin a.
```

The Pauli $`P=Z_{m-2}`$ conjugates the scalar word to its adjoint. Its
normalized Pauli-transfer correlation for R repetitions is $`T_{2R}(c)`$.
Writing d=m-2 and p=2^d-1, the reduced rational Chebyshev values obey

```math
T_k(c)=\frac{q_k}{2^{k(d-1)+1}},\qquad
q_1=p,\quad q_2=p^2-2^{2d-1},
```

```math
q_k=pq_{k-1}-2^{2d-2}q_{k-2}.
```

Each numerator is odd. The retained
[returned-helper Pauli-transfer denominator argument](../../docs/OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers)
therefore proves

```math
T_{\mathrm{exact}}(\mathcal S_f^R)\ge4R(m-3)+2.
```

The displayed literal scalar implementation uses 6m-8 T gates per repetition.
The lower bound allows exactly returned stabilizer and arbitrary dirty helpers.
It does not apply to approximate resynthesis or the unrestricted Hopf target.

### One paid axis change removes the repetition count under approximation

Define

```math
Y_*={i(M\Gamma-aI)\over\sqrt{1-a^2}},\qquad
\mathcal S_f=e^{i\theta Y_*\otimes X_\sigma}.
```

The matrix $`Y_*`$ is a Hermitian involution anticommuting with Gamma. Put

```math
\beta=\tfrac12\mathrm{atan2}(\sqrt{1-2a^2},-a^2),\qquad
\phi=-\tfrac12\mathrm{atan2}(a,\sqrt{1-2a^2}),
```

```math
A_m=e^{i\phi\Gamma}e^{i\beta M},\qquad
A_m\Gamma A_m^\dagger=Y_*.
```

The identity follows in the two-generator algebra by setting Gamma=Z and
$`M=aZ+\sqrt{1-a^2}X`$; it holds on every physical copy of that algebra.
Consequently

```math
\mathcal S_f^R=(A_m\otimes I)
 e^{iR\theta\Gamma\otimes X_\sigma}(A_m^\dagger\otimes I).
```

The source rotation uses $`e^{i\beta M}=U_m e^{i\beta X_0}U_m^\dagger`$,
with the actual geometric loader. There are four loader appearances and
five Pauli-rotation appearances. The zero-helper one-qubit specialization of
[borrowed synthesis](../../docs/BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements),
followed by Clifford parity conjugations, gives

```math
T=O(m+\log(1/\epsilon)),\qquad \text{additional helpers}=0.
```

Using the actual inverse of the synthesized basis gives full error at most
$`2\epsilon_\phi+2\epsilon_\beta+\epsilon_{\mathrm{mid}}`$, independent
of R. Computing the accumulated angle needs additional classical precision,
not repeated quantum calls. Signal-only Clifford interlayers also compress
through this same basis into one two-qubit operation. All signal/core inputs
may be arbitrary. No intermediate core return is claimed.

## 2. A paid common basis for several masks

For h finite masks partition coordinates by their complete h-bit sign pattern.
Let C range over the q nonempty classes, with $`q\le\min(m,2^h)`$, and put

```math
w_C=\sum_{j\in C}a_j^2,\qquad
G_C=\sum_{j\in C}{a_j\over\sqrt{w_C}}\Gamma_j.
```

These are anticommuting Hermitian involutions. For p=min(C), the normalized
tail $`A_p=2^{p/2}\sum_{j\ge p}a_j\Gamma_j`$ has unit coefficient norm.
Sign-reverse its coordinates outside C to obtain B_p. Then

```math
\{A_p,B_p\}/2=c_C I,\qquad c_C=2^{p+1}w_C-1\in[0,1],
```

```math
G_C={A_p+B_p\over\sqrt{2+2c_C}},\qquad
{a_p\over\sqrt{w_C}}\ge1/\sqrt2.
```

The first coefficient contributes at least half the tail weight. Thus even an
exponentially small class has well-conditioned alignment overlaps; there is
no inverse-class-weight quantum cost.

### Explicit full and half moves

For Hermitian Clifford vectors A,B with $`\{A,B\}=2cI`$, $`0\le c\le1`$,
set $`t=\sqrt{(1+c)/2}`$. Choose

```math
r={1\over2t}\quad\text{for a full move},\qquad
r={1\over4t\sqrt{(1+t)/2}}\quad\text{for a half move}.
```

Define

```math
\beta=\arcsin r,\quad
\zeta=\mathrm{atan2}(cr,\sqrt{1-r^2}),\quad
\alpha=\pi/4-\zeta/2,\quad \gamma=-\pi/4-\zeta/2,
```

```math
U(A,B)=e^{i\alpha A}e^{i\beta B}e^{i\gamma A}.
```

The full move sends A to B; the half move sends A to
$`(A+B)/\sqrt{2+2c}`$. To verify this, set $`A=Z`$ and
$`B=\cos\vartheta Z+\sin\vartheta X`$. The product has opposite real
off-diagonal entries of magnitude $`r\sin\vartheta`$ and diagonal
$`\sqrt{1-r^2\sin^2\vartheta}`$. The formulas are continuous at c=1
and produce identity there.

Apply a half move to the signed tails and then a full move from Gamma_p to
G_C. The latter complete rotor lies inside the class span and fixes every
Gamma outside C, even though its implementation temporarily uses the tail.
It uses nine axis exponentials: two Gamma_p rotations and seven tail or
masked-tail rotations. Its fixed loader charge is at most 14(m-p-1) T gates.
The complete class basis therefore costs

```math
O(qm+q\log(q/\epsilon_B)).
```

Approximate class rotors need not commute: synthesize their fixed ordered
product and use its actual inverse.

### Compress the even algebra, retaining central parity

The source-vector dimension is

```math
d=\dim_{\mathbb R}\mathrm{span}\{M,N_{f_1},\ldots,N_{f_h}\}\le h+1.
```

Since every a_j is nonzero, this is the exact rational rank of the sign matrix
with row $`(1,(-1)^{f_{1j}},\ldots,(-1)^{f_{hj}})`$. No floating rank
threshold or unknown-real equality test is required. After class alignment,
signed coordinate exchanges and at most O(q^2) coefficient-plane rotations
put the vectors into the first d Clifford coordinates. Every plane rotation
is a physical Pauli rotation, with its parity work and precision charged.

The scalar word uses the even algebra, generated by the independent
anticommuting Paulis $`P_j=i\Gamma_0\Gamma_j`$, j=1,...,d-1. Binary symplectic
elimination puts it on

```math
s_* = \lfloor d/2\rfloor
```

core qubits, with identity spectators. For odd d-1 a central parity sector is
retained as an arbitrary qubit, not assigned an eigenvalue. The entire common
basis has cost $`O(qm+(q^2+q)\log((q+1)/\epsilon_B))`$.

For k interacting non-core qubits, any finite program using these sources,
actual inverses and Clifford interlayers becomes one operation on w=k+s_* qubits
in this basis. Multiply it classically before synthesis. A finite decomposition
into its known Clifford scalar phase and selected SU(2) rotations, with selectors
expanded into commuting Pauli terms, gives a conservative O(8^w) Pauli-rotation
count. Arbitrary non-Clifford scalar interlayers are not included in this
zero-additional-helper statement.
Thus

```math
T=O\left(qm+(q^2+q+8^w)
 \log{q^2+q+8^w\over\epsilon}\right).
```

For fixed h,k this is $`O_{h,k}(m+\log(1/\epsilon))`$, with no additional
helpers and no program-length factor. The guarantee reproduces the raw
program's whole core action. The construction does not replace a growing w
by a constant, or discard scalar phases from the compact product.

## 3. Growing source rank and the signal-only exception

### Actual grouped coefficients can force a growing common core

In the [actual grouped tables](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#3-a-small-coefficient-table-and-a-streamed-coarse-circuit),
let K be the padded term count and s the group height. Choose v at least
three with $`2^{-v}\le1/(8s)`$, and h+1 distinct external prefixes. Set

```math
c_*=2^{-v},\qquad c_j=2^{-v}-2^{-(v+j)}\quad(1\le j\le h).
```

At prefix x activate only the first local root edge with
$`\sin\theta_x=c_x/(K2^{s/2})`$ and positive cosine. Identity is a legal
actual coarse word: its error is below $`1/(4sK2^{s/2})`$. The displayed c's
are exactly zero-column forward rows. All other edges are identity.

At source width $`m\ge v+h+2`$, canonical rounding yields

```math
f_*={\bf1}_{\{1,\ldots,v\}},\qquad
f_j=f_*\mathbin{\oplus}e_{v+j}.
```

The entries lie at exact grid centers, so the permitted certified-estimate
error cannot change the rounded integer. There are h+2 sign classes and
h+2 independent source vectors. Other actual coefficient rows cannot reduce
that span. A fixed legal height/suffix leaves a fixed positive fraction of
N external prefixes, so h can grow proportionally to N.

This growth survives approximate simultaneous core-only compression. Put
p=h+1, $`\alpha_j=a_{v+j}`$ and
$`\beta=\sqrt{(1-c_*)/(1+c_*)}`$. The Hermitian generators

```math
K_0={iD_*\over\sqrt{1-c_*^2}},\qquad
K_j={-i(D_j-D_*)\over2\alpha_j}
```

have anticommutator Gram matrix

```math
G=\begin{pmatrix}1&\beta\alpha^T\\
\beta\alpha&I-\alpha\alpha^T\end{pmatrix},
\qquad \lambda_{\min}(G)\gt11/16.
```

Indeed $`\|\alpha\|^2\lt2^{-v-1}\le1/16`$. Whitening gives p
anticommuting Hermitian involutions whose 2^p ordered monomials are
orthonormal in the full normalized trace. They are even products of p+1
independent original Clifford directions; no chirality is selected.

If one core-only basis makes all scalar words epsilon-close to b active core
qubits plus the signal, extraction of their signal-odd coefficient is
contractive. The generator error is at most epsilon divided by
$`\alpha_{\min}=2^{-(v+h+1)/2}`$. Whitening increases it by at most
$`4\sqrt{p/11}`$. Rescaling to contractions and telescoping monomials gives

```math
4^b\ge2^p(1-\tau^2),\qquad
\tau\lt3p^{3/2}2^{(v+h+1)/2}\epsilon.
```

This is Bessel's inequality in the common matrix algebra. For tau below 1/2,
$`b\ge\lceil(h+1)/2\rceil`$. At epsilon=eta,
$`v=\lceil\sqrt N\rceil+6`$, h at most N/16 and n at least eight, the bound
on tau is below 1/2. This is not a T-count lower bound. It restricts one
prefix-independent, core-only representation of every full scalar word;
the logical target is itself an inexpensive disjoint matching.

### Signal-only programs have polynomial Lie size

A growing Hilbert representation need not require generic dense synthesis.
After the paid basis, write

```math
\mathcal S_f=c_f I-i\sum_{j=1}^p b_{f,j}P_j\otimes X,
\qquad c_f^2+\sum_jb_{f,j}^2=1.
```

For signal Pauli matrices sigma_1=X, sigma_2=Y, sigma_3=Z, define

```math
E_{ij}=P_iP_j\otimes I,\qquad
E_{i,p+a}=iP_i\otimes\sigma_a,\qquad
E_{p+a,p+b}=I\otimes\sigma_a\sigma_b.
```

Use antisymmetric indices. These independent anti-Hermitian matrices satisfy

```math
[E_{ij},E_{kl}]=2(\delta_{jk}E_{il}-\delta_{ik}E_{jl}
 -\delta_{jl}E_{ik}+\delta_{il}E_{jk}).
```

They represent so(p+3), of dimension $`\binom{p+3}{2}`$. Scalar-source
calls and signal SU(2) interlayers lie in its spin representation. Track the
known Clifford scalar phase separately. Multiplying one (p+3)-dimensional
real orthogonal matrix classically and applying Givens elimination gives
O(p^2) Pauli rotations. The lift's sign must be tracked or determined by
certified comparison with the original finite word; it is not discarded.
The classical coordinates are not additional qubits. Including the basis,

```math
T=O\left(qm+(q^2+p^2+q)
 \log{q+p+1\over\epsilon}\right).
```

For the displayed table family only two classes are nonsingletons, so the
loader term improves to O(m), with $`O((h+1)^2\log((h+2)/\epsilon))`$ fine
cost. Neither is linear at h proportional to N and epsilon=eta.

This closure does not survive arbitrary signal/address interlayers. With two
external qubits, distinct commuting Pauli labels can multiply without
cancelling. For example

```math
[iP_i\otimes X\otimes I,\ iP_j\otimes X\otimes Z]
=-2P_iP_j\otimes I\otimes Z\quad(i\ne j).
```

Closing the generators $`iP_i\otimes Q`$ and $`I\otimes iQ`$ over all
nonidentity two-qubit Q produces every noncentral label, of dimension
$`2^{p+4}-2^{p\bmod2}`$. One verifies generation inductively using the
symplectic commutator rule: nonorthogonal labels produce their binary sum;
the only excluded labels are the radical of the core form. For even p that
radical is zero; for odd p it has dimension one. This is a larger control
envelope, not a claim that the actual grouped word explores it. Exponential
Lie dimension is not an unrestricted native lower bound.

## 4. Conserved-prefix raw programs have a coherently priced basis

Let y be a conserved a-bit prefix, Q=2^a. At each y there are at most h masks
and k interacting non-core wires. The mask union over all y can be large.
For fixed h,k and $`L=\lceil\log_2(1/\epsilon)\rceil`$, there is a fixed
constant c_(h,k) such that

```math
m\ge\max\{L,2a\}+c_{h,k}
```

suffices for a whole raw-program compiler with two temporary initialized
helpers and

```math
T=O_{h,k}(m+Q+L).
```

The original signal is included in k and may be occupied. Both added helpers
are returned in the full isometry norm. Program length affects classical
multiplication and parameter evaluation, not the number of fine quantum calls.
The following identities price the prefix-dependent basis rather than treating
it as a coherent oracle.

### One full loader selects every tail

For p at least one the same full geometric loader obeys

```math
K_p=U_m\Gamma_pU_m^\dagger
={-\Gamma_{p-1}+A_p\over\sqrt2},\qquad
A_p=-K_p\Gamma_{p-1}K_p.
```

For p=0, $`A_0=U_m\Gamma_0U_m^\dagger`$. The identity is obtained by
tracking the first source rotation touching Gamma_p and then its tail.
A prefix-selected Gamma is an addressed X/Z word. Self-borrowed whole-word
queries cost O(Q) T and O(Qm) Clifford gates when m is at least 2a. They return
selectors before another operation changes the core. Identity entries handle
p=0 separately; a dirty sign echo retains the minus sign on the other rows.
There are six full-loader appearances in the displayed complete selection.
Therefore SELECT of arbitrary p(y) costs O(m+Q) T with no extra helpers.

Likewise

```math
e^{i\theta(y)A_{p(y)}}=
K_{p(y)}e^{-i\theta(y)\Gamma_{p(y)-1}}K_{p(y)}
```

on the p-positive rows. For a prefix-selected Hermitian Pauli P_y, introduce
an initialized eigenvalue helper t and use

```math
F_y=H_t C_t(P_y)H_t,\qquad
F_y e^{i\theta(y)Z_t}F_yJ_t
=J_t e^{i\theta(y)P_y}.
```

The controlled Pauli is a completed dirty query with its scalar retained.
The middle addressed diagonal uses a distinct initialized phase signal b and
the original core as arbitrary dirty precision work. Helper t is coherent
address data, not a second zero input during that call. A q-bit precision
choice gives error below $`30\,2^{-q}`$ and cost O(Q+q), by the
[one-clean diagonal primitive](../../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors).
All selectors can self-borrow when the stated fixed slack is available.
The primitive permits reference-entangled core inputs, so using the entangled
core as dirty work is valid.

### Prefix-dependent alignment and complete compact products

Apply Section 2's class alignment coherently, with p_C and signs supplied by
each prefix table. The selected-tail identity shares every geometric loader;
there is no O(mQ) T charge. Pack representatives and rotate the coefficient
span with prefix-selected Pauli rotations. Determine each source rank from
its exact finite sign matrix, pad to $`d_0=\min(m,h+1)`$, and retain
$`s_*=\lfloor d_0/2\rfloor`$ active core qubits including central parity.

The resulting $`\mathcal B=\sum_y|y\rangle\langle y|\otimes B_y`$
commutes with every prefix-preserving non-core interlayer. Thus the full raw
program equals one B, one prefix-selected compact product on w=k+s_* qubits,
and the actual B inverse. Its conservative count is

```math
T=O\bigl(qm+(q^2+q+8^w)(Q+\ell)\bigr),\qquad q\le\min(m,2^h).
```

Take a fixed rotation-count bound F depending on h,k and choose precision
at least $`L+\lceil\log_2(120F)\rceil`$, with enough arbitrary core spectators
for that phase source and self-borrowed queries. This is the meaning of
c_(h,k). Certified finite selection includes singular charts and literal
compact phases. Basis and compact errors contribute $`2\epsilon_B+\epsilon_V`$,
not a program-length multiple. For growing h or k these constants and the
spectator allocation must be kept explicitly.

A raw prefix change still has the exact connector

```math
\mathcal B^\dagger(|y'\rangle\langle y|\otimes I)\mathcal B
=|y'\rangle\langle y|\otimes B_{y'}^\dagger B_y.
```

Dropping this connector is not justified by cheap selection. The next section
explains the different rule after complete recovery.

## 5. Completed groups admit independent gauges without connectors

Let J initialize the designated external flags, F be the unitary target on all
unknown logical and dirty inputs, and G an ambient unitary with

```math
GJ=Jg,\qquad [g,F]=0.
```

For any actual unitary V, set $`K=G^\dagger VG`$. Direct substitution gives

```math
KJ-JF=G^\dagger(VJ-JF)g,\qquad
\lVert KJ-JF\rVert=\lVert VJ-JF\rVert.
```

This extends the existing dirty-only
[unitary-invariance observation](SOURCE_REUSE_LIMITS.md) to a conditional
prefix gauge. It uses no trace or channel conversion. With approximate
commutation, add $`\lVert g^\dagger Fg-F\rVert`$ to the bound.

The complete recovered Hopf group preserves its external prefix and targets

```math
F_g=\sum_y|y\rangle\langle y|\otimes W_{g,y}\otimes I_{\mathrm{dirty}}.
```

A gauge controlled by y and acting only on an arbitrary core subset commutes
with this target and preserves J. It does not act on the logical suffix or
external flags. On inactive suffix inputs, conjugation still preserves exact
identity. The full group V_g must include SELECT, rejected sectors, recovery,
actual coarse action, and predicate uncomputation; a scalar accepted block
alone does not satisfy this premise.

If $`\lVert V_gJ-JF_g\rVert\le\delta_g`$ and an emitted replacement obeys
$`\lVert(\widehat K_g-K_g)J\rVert\le\xi_g`$, the complete hybrid is

```math
\left\lVert\widehat K_R\cdots\widehat K_1J
-JF_R\cdots F_1\right\rVert
\le\sum_g(\delta_g+\xi_g).
```

The gauges may differ, fail to commute, use different dirty subsets, or depend
on a prefix changed later. No G_(g+1)-inverse G_g connector is executed.
This is not a claim that the two physical products are equal or related by
one common conjugation: each is a valid approximant to the same target.
What must be compiled is the transformed interior, not G itself.

### Exact counterexample to deleting a raw connector

On two core qubits set Gamma_0=X tensor I, Gamma_1=Z tensor X, and
$`M=(\Gamma_0+\Gamma_1)/\sqrt2`$. The mask P=Z tensor I gives c=0 and
$`D=(M(PMP)-(PMP)M)/2=\Gamma_0\Gamma_1`$. Then
$`S=X_\sigma\otimes D`$ has $`S^2=-I`$. Let G apply P only on prefix y=1.
It gives $`GSG=Z_yS`$, so the complete identities are

```math
SX_yS=-X_y,\qquad (GSG)X_y(GSG)=+X_y.
```

The norm error is two. Both accepted scalar coefficients are zero. This exact
Gaussian-integer witness separates raw source action from completed-group
return and retains the literal sign.

### Price the remaining compact interior, not an omitted basis

If the group gauge gives a prefix-selected interior on w active qubits,
including all genuinely used logical, label, signal, rejection and recovery
registers, synthesize that interior directly. An O(8^w) fixed-axis Pauli schedule
uses one additional temporary phase signal. A fixed Clifford sends each axis
to one active Z; unlike prefix-dependent axis selection, no eigenvalue helper
is needed. A sufficient explicit allocation for F rotation appearances is

```math
\ell\ge\max\{\lceil\log_2(1/\xi_g)\rceil+
\lceil\log_2(120F)\rceil,\ 2(a+2),\ 5\},\qquad
m-s_*\ge\ell+2.
```

Its cost is $`O(8^w(2^a+\ell))`$ plus the original polynomial predicates.
For conditional fixed-height groups, the temporary signal occupies an extra
active-suffix position with exact inactive conditioning. Used private flags
remain arbitrary inputs of the compact body; they are not reused as fresh
zero ancillas. Unused suffix-zero tests can be removed only after proving that
the body leaves them untouched. The allocation and O(8^w) factor are not
asserted uniform in growing height. No linear total fine cost follows yet.

### Independent completed-group accuracy forces repeated precision

For any fixed, angle-independent partition into R nonempty depth groups,
the additive objective $`\sum_gT(\widehat K_g)=O(N)`$ is impossible
when R grows and each completed group retains the endpoint error guarantee.
Choose one legal angle $`\pi/16`$ in every group and set all other angles
to zero. Each group is then one two-level rotation on a computational-tree
edge. Its two labels differ in one bit, so a stabilizer eigenvector has
literal eigenphase $`e^{\pm i\pi/16}`$. The
[single-angle precision bound](../../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
therefore gives

```math
T(\widehat K_g)\ge2N-1,\qquad
\sum_{g=1}^R T(\widehat K_g)\ge R(2N-1).
```

An admissible gauge preserves the group target and cannot remove this
witness. The hybrid budget $`\sum_g(\delta_g+\xi_g)\le\eta`$
implies that every separately emitted group has error at most eta.
This is an application of the existing
[modular-bank separation](ORDERED_BLOCK_BANKS.md#8-a-modular-bank-separation),
not an unrestricted lower bound.

The combined witness acts on at most $`2R\le2n`$ computational labels,
with identity on the complement. Its one nontrivial SU block has squared
dimension $`O(n^2)=O(N)`$, so the disjoint-block theorem compiles the
whole target in O(N) T gates. The remaining constructive objective must
therefore permit joint native synthesis across group boundaries or use
intermediate guarantees that do not require separately accurate returned
groups. The gauges above remain valid, but do not make their costs additive
at linear total count. A [global native-word formulation](JOINT_NATIVE_WORDS.md)
keeps this distinction explicit.

## 6. Endpoint-only normal form with one fixed physical source core

Use the original two-external-clean group layout and the original two-flag,
one-tail compiler for the fixed deepest tail. At L=N set, once,

```math
\mu=N+\lceil\log_2(64n)\rceil=N+6+\lceil\log_2 n\rceil.
```

Reserve these same mu physical dirty wires for every source call. The external
dirty remainder is $`n+1-\lceil\log_2n\rceil\ge2`$. A tail query has at
most n+1 address bits. A group query has
$`k\le n-r+\log_2K+O(1)`$, with $`\log_2K\le r/C_1+O(1)`$.
Enlarge the existing fixed sufficient suffix threshold once so k is at most
n+2. Since $`\mu\ge2(n+2)`$ for n at least three, all coefficient queries
use the retained self-borrowed whole-word construction. They may change which
core bits serve as selectors internally, but each complete query returns them.

Two outside dirty wires suffice for a suffix-toggle bit and a bounded-source-
control helper. A three-controlled source center has a returned one-dirty-bit
Toffoli echo. Long logical-suffix tests use the existing returned borrowed local
bit; coarse-symbol lookup may temporarily borrow core wires because its complete
word returns them exactly. The active suffix supplies its private initialized
work only on the active sector and is not also counted as source storage.
Thus the literal total width remains

```math
n+2+\mu+(N+n+7-\mu)=N+2n+9.
```

At width mu, group error is at most $`10\,2^{-\mu}`$ and tail-layer error
at most $`10\sqrt2\,2^{-\mu}`$. There are at most n groups-plus-tail-layers,
so

```math
\epsilon_{\mathrm{baseline}}
\lt15n2^{-\mu}\le(15/64)\eta.
```

The coefficient tables are re-rounded at the common width. Known inactive
coefficients retain their exact encodings. The actual coarse word is fixed
before its residual is evaluated; its discrepancy is not added again.
The complete-isometry hybrid propagates previous leakage without resets.

The old total table-row cost is still O(N); self-borrowing changes constants
only. There are $`O(1+\log_2^*(n+2))`$ groups and a fixed tail. Hence

```math
T=O\bigl(N+\mu[1+\log_2^*(n+2)]\bigr),\qquad G=O(N^2).
```

This preserves the known bound. It is specific to L=N and two available clean
flags; it does not replace the general-precision or one-clean theorem. A fixed
source permits exact exterior loader hoisting, and the returned-work target
permits omitting a global dirty-only exterior gauge. The surviving masks become
$`U_\mu^\dagger P_fU_\mu`$; their native cost is not zero. Width transitions
are no longer compulsory, but the growing fine interior remains unresolved.

## 7. Reproducibility and interpretation

The shared [source module](../../compiler_robust_hopf/source_resynthesis.py)
contains rational mask/rank/denominator checks, bounded common-basis matrices,
and exact binary Lie closures. The
[group module](../../compiler_robust_hopf/source_group_resynthesis.py) builds
actual rational coefficient tables and complete half-block/recovery words in
faithful compact Clifford coordinates. It checks both central sectors, not
2^m physical amplitudes at large m. The
[native module](../../compiler_robust_hopf/native_source_resynthesis.py)
regenerates selected-tail words using the retained exact Toffoli/phase helpers.
It emits no arbitrary-angle fine resynthesis compiler.

The [source tests](../../tests/test_source_resynthesis.py),
[group tests](../../tests/test_source_group_resynthesis.py), and
[native tests](../../tests/test_native_source_resynthesis.py) check the distinct
interfaces. Exact arithmetic is used for dyadic masks, source ranks, denominator
growth, the raw sign counterexample and the finite native selected tail. Floating
QR, SVD, common-basis and recovery checks are bounded diagnostics, not certified
endpoint simulations or proofs of uniform synthesis complexity.

The unrestricted comparison remains

```math
2N-1\le T^\star\le O(N[1+\log_2^*(n+2)]).
```

The usable changes are jointly resynthesizing raw
programs under explicit algebra promises, omitting eligible gauges only after
completed recovery, and using one fixed source width at the selected endpoint.
None licenses a free growing compact gate, uncharged quantum basis, omitted
rejected paths, or addition of primitive lower bounds across a fused circuit.
