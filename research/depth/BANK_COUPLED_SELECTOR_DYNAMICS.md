# Bank-coupled selector dynamics and complete acceptance

[Working checkpoint](../../WORKSPACE.md) ·
[Co-moving stages](CO_MOVING_LOOKUP_STAGE.md) ·
[Selector-source screens](ENTANGLING_SELECTOR_ACCEPTANCE.md) ·
[Complete Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Bank-dependent selector operations can leave the common-interlayer premise
of the earlier stage-rank screen. Their actual action, rather than this
change of premise, determines whether they improve acceptance. This note
pays two such layers, derives their complete bank/selector dynamics, and
identifies the fixed spaces that survive paired and unpaired use.
A Bell-cell witness strengthens the earlier finite transpose-zeta failure
to an all-size obstruction for those specified frames, even with added
bank-controlled swaps. The general native-depth question remains open.

## 1. Registers and the exact obligation

Use p address bits, $`m=2^p`$ arbitrary bank bits Y, and two arbitrary
selector/work registers s,t, each of length m. The second register is
**principal arbitrary work**, not initialized storage. Write bitwise AND
as $`\odot`$, complement as an overbar, and define

```math
h_Y(s)=\sum_uY_us_u\pmod2,\qquad
n_C(x,Y,s)=\sum_uY_us_u d_C(x,u),\qquad
d_C(x,u)=\sum_{c\in S_C}(x_c\oplus u_c).
```

The retained stage is $`R_C=i^{n_C}L_C^{h_Y(s)}`$, where L_C
relabels x,Y,s together and leaves t alone. Its complete native word and
actual adjoint have depth at most 116, $`T,b=O(pm)`$ and
$`G=O(m^2+pm)`$, with zero clean work and zero error.

For an initialized signal embedding J and an actual unitary U, the
accepted block $`K=J^\dagger UJ`$ must satisfy

```math
\left\|K^\dagger K-\frac{I+\Omega}{4}\otimes I_t\right\|\le\delta,
\qquad \Omega|x,Y,s\rangle=(-1)^{Y_x}|x,Y,s\rangle.       \qquad(1)
```

Every t input and every reference is retained. The accepted fiber may
move t, address, bank and selectors. The paid four-call Gram interface
enforces final return once (1) is established; intermediate t movement
is not itself a source failure. In contrast, the implementation helpers
below return after every component and may be reused. The unpaid step
is complete acceptance, not a new assumption of source-level t return.

<a id="paid-bank-couplings"></a>

## 2. Paid bank-controlled swaps and Hadamards

Define the exact involutions

```math
B_Y=\prod_u\mathrm{CSWAP}(Y_u;s_u,t_u),\qquad
D_Y=\prod_u H_{s_u}^{Y_u}.                              \qquad(2)
```

For q rows on disjoint triples, one common forward/actual-inverse ledger
in the positive-S native alphabet is:

| Operator | T-count | Clifford count | T-depth | Returned dirty helpers |
| --- | --- | --- | --- | --- |
| B | $`7q`$ | $`14q`$ | 4 | 0 |
| Signal-controlled B | $`26q+2(q\bmod2)`$ | $`\le48q+6`$ | 16 | q |
| D | $`2q`$ | $`15q`$ | 2 | 0 |
| Signal-controlled D | $`8q+(q\bmod2)`$ | $`\le22q+3`$ | 6 | 0 |

The full bank uses q=m. No additional clean qubits or approximation error
are required. Counting the external control as clean is a separate source
choice. All identities hold for arbitrary controls and coherent banks.

To obtain the swap rows, emit a Fredkin chronologically as

```text
CX(s,t), H(s), CCZ(Y,t,s), H(s), CX(s,t).
```

The disjoint CCZ layer has seven T gates and ten Cliffords per row in four
T layers. For a signal-controlled swap, leave the outer Clifford walls
unconditional and control only the middle CCZ product. When the signal
is zero the walls cancel. The
[controlled nonlinear phase](ENTANGLING_SELECTOR_ACCEPTANCE.md#nonclifford-selector-control)
costs $`26q+2(q\bmod2)`$ T gates, at most $`44q+6`$ Cliffords,
16 T layers and q returned dirty helpers. The four outer Cliffords per
row give the displayed bound, including actual reversal.

There is also a direct interpretation of the signal control. Borrow
arbitrary r_u, put
$`C=\prod_u\mathrm{CSWAP}(r_u;s_u,t_u)`$ and
$`A=\prod_u\mathrm{CCX}(a,Y_u;r_u)`$, and emit C,A,C,A.
The two swap exponents XOR to aY_u, and every r_u returns. This identity
does not initialize t; the middle-CCZ construction above has the better
Clifford bound.

For the Hadamard rows, use the retained literal
$`VZV^\dagger=H`$, with $`V=SHTHS^\dagger`$.
Ordinary controlled-H is V-adjoint, CZ, V. Signal control instead uses
V-adjoint, $`\prod_u\mathrm{CCZ}(a,Y_u,s_u)`$, V.
The shared-root CCZ batch costs $`6q+(q\bmod2)`$ T gates,
at most $`10q+3`$ Cliffords and four T layers. The V walls add
2q T gates, at most 12q Cliffords and two layers. Distinct Y_u controls
avoid the two-shared-control construction's extra borrowed rectangle.
All inverse S powers and literal phases are charged.

<a id="complete-swap-dynamics"></a>

## 3. Complete paired and unpaired swap dynamics

First consider $`B_YR_CB_Y`$, with the rightmost factor acting first.
Put

```math
a=h_Y(t),\qquad P=C^a,\qquad
t_B=\overline Y\odot t\oplus Y\odot s.
```

Its complete output and phase are

```math
\begin{aligned}
x'&=Px,&Y'&=PY,\\
s'&=\overline{PY}\odot Ps\oplus(PY)\odot t_B,\\
t'&=\overline{PY}\odot t_B\oplus(PY)\odot Pt,\\
\text{phase}&=i^{n_C(x,Y,t)}.
\end{aligned}                                                \qquad(3)
```

The first swap places t at every occupied selector row, so both the
stage trigger and its quarter phase read t there. Relabeling and the
second swap give (3). Both a and its integer phase exponent are invariant
under this output permutation. The actual inverse therefore has the same
basis permutation and opposite phase. In particular, a zero trigger can
leave a nontrivial t-dependent quarter phase even when every bit returns.

The dirty-register movement is exactly

```math
t'\oplus t=
\overline{PY}\odot Y\odot(s\oplus t)
\oplus(PY)\odot(Pt\oplus t).                              \qquad(4)
```

Only nonfixed rows in the union of Y and PY can change. A transposed pair
with two occupied rows swaps its t bits; with two vacant rows it swaps
only s; a singly occupied pair can transfer data between registers.
Every fixed physical row, including row zero, returns its t bit.
This is an exact description of principal work, not permission to
discard that register.

For a genuinely unpaired word $`U=R_DB_YR_C`$, define

```math
h=h_Y(s),\qquad P=C^h,\qquad
k=h_{PY}(t),\qquad Q=D^k.
```

The complete action is

```math
\begin{aligned}
x'&=QPx,&Y'&=QPY,\\
s'&=\overline{Y'}\odot QPs\oplus Y'\odot Qt,\\
t'&=\overline{PY}\odot t\oplus(PY)\odot Ps,\\
\text{phase}&=i^{n_C(x,Y,s)+n_D(Px,PY,t)}.
\end{aligned}                                                \qquad(5)
```

The second stage now reads t against the bank transported by the first.
Thus the operation creates bank-dependent memory. Its actual inverse is
$`R_C^\dagger B_YR_D^\dagger`$: reverse the matching order as
well as both phase signs. A same-order word with negated phases is not
generally the inverse.

At x=0,Y=e_0, however, (5) simplifies to a swap of s_0,t_0 and the
nonzero-selector permutation $`D^{t_0}C^{s_0}`$, with phase one.
Every word using only these stages, their adjoints and B_Y fixes the
entire subspace $`s_0=t_0=0`$. On the empty bank every such generator
is identity. Signal-only mixing and controlled uses cannot distinguish
these fixed inputs. If their common acceptance is a, then

```math
\delta\ge\max\{a,|a-1/2|\}\ge1/4.                       \qquad(6)
```

Adding D_Y alone replaces the witness rather than removing it. Choose
both s_0,t_0 in the positive Hadamard eigenstate
$`|h_+\rangle=\cos(\pi/8)|0\rangle+\sin(\pi/8)|1\rangle`$
and every nonzero selector/work wire in $`|+\rangle`$. The nonzero
selector tail is permutation-invariant; R_C, B_Y and D_Y all fix the
marked-origin vector. Their empty-bank restrictions are identity. The
same quarter floor holds for this generator grammar. Additional selector
operations can leave this common ray.

<a id="wrapper-and-pruning-limits"></a>

## 4. What wrappers and occupied-row tags can change

There is a full-space endpoint statement independent of any particular
witness. Let $`H_Y=(I+\Omega)/4`$, extend K by arbitrary dirty work W,
and let E_in,E_out be signal-independent unitaries. If E_in commutes
with $`H_Y\otimes I_W`$, then

```math
\begin{aligned}
\widetilde K&=E_{\rm out}(K\otimes I_W)E_{\rm in},\\
\|\widetilde K^\dagger\widetilde K-H_Y\otimes I_W\|
 &=\|K^\dagger K-H_Y\|.
\end{aligned}                                                \qquad(7)
```

The output unitary cancels; conjugation by E_in preserves norm and the
target. This includes arbitrary bank/address-preserving selector-work
wrappers. It does not apply to arbitrary interleaving through the source.

For stage-dependent whole-branch conjugations by B_Y, any old common
marked vector in $`s_0=0`$ survives after assigning every distinct
principal t_0 input zero. This is a permissible dirty input used to test
the complete norm, not a preparation assumption. Reusing a t register is
allowed. Unpaired swaps between already completed fixed branches also
preserve this vector.

A different proposed wrapper, replacing a selector frame J by
$`B_Y^\dagger(J\otimes I_t)B_Y`$, is more restrictive at the
origin: it moves J's physical s_0 wire to t_0, so the new frame does not
act on s_0 at all. It preserves the entire s_0=0 sector, and conjugating
R_C by such frames fixes that sector pointwise. These two wrapper
architectures should not be conflated.

Occupied-row tags also have an exact pruning identity. Let
$`T=\prod_u V_u(x)(s_u)`$ obey the literal rowwise covariance
$`V_{Cu}(Cx)=V_u(x)`$, and factor it into its occupied and vacant
parts, $`T=A I_{\rm vacant}`$. Global commutation up to compensating
row phases alone is not this premise. The vacant part commutes
with h_Y(s), n_C and L_C: it acts where Y_u=0 and follows the same
relabeling. It also commutes with A. Hence

```math
R_C^\dagger T^\dagger R_CT
 =R_C^\dagger A^\dagger R_CA.                             \qquad(8)
```

For example D_Y and global H_s have exactly the same completed echo.
The same holds for covariant mismatch-dependent row tags restricted to
occupied rows. The retained
[same-bank singlet proof](CO_MOVING_SELECTOR_READOUT.md#5-a-two-row-singlet-shield-for-arbitrary-angles-and-local-tags)
therefore still gives a one-half floor for the balanced difference source
and a quarter floor for its specified completed-echo scalar filters.

Likewise $`D_X(Y)=\prod_u\mathrm{CX}(Y_u;s_u)`$ replaces h by
$`h+|Y|\pmod2`$ and n_C by N_C-n_C, where
$`N_C=\sum_uY_ud_C(x,u)`$. It commutes with L_C, giving exactly
the earlier global-flip echo:

```math
R_C^\dagger D_XR_CD_X=i^{N_C-2n_C}L_C^{|Y|\bmod2}.          \qquad(9)
```

There can be no universal extension to all unpaired bank access. Put
$`D_Z(Y)=\prod_u\mathrm{CZ}(Y_u,s_u)`$ and let cF flip every
selector conditioned on signal z. The complete Clifford identity is

```math
(\mathrm cF)D_Z(\mathrm cF)D_Z
 =(-1)^{z|Y|}I_s.                                         \qquad(10)
```

It distinguishes the empty bank from an origin singleton while returning
all selectors. It reads bank parity, which is unchanged between addresses
of the same bank. Thus breaking that one comparison does not supply
addressed lookup. General unpaired words require a new complete acceptance
calculation; none of these limits is an unrestricted native-depth bound.

<a id="bell-cell-witness"></a>

## 5. An all-size witness for the specified nonlocal frames and swaps

Retain exactly the three local frames from the
[transpose-zeta scout](ENTANGLING_SELECTOR_ACCEPTANCE.md#nonlocal-selector-scout):
$`J_j=Q_jH_{A_j}P_jH_{A_j}Q_j^\dagger`$, with one frame repeated
in every four-row cell, $`Q_j=X^{u_j}Z^{v_j}`$, and
$`V_j=N_pJ_jN_p^\dagger`$. Here $`N_p=M_p^T`$ is the already
paid transpose-zeta involution. The low two row-label bits specify corners.

| Frame | Omitted CCZ corner | H mask | X mask | Z mask |
| --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 15 | 2 | 1 |
| 3 | 3 | 3 | 4 | 5 |

For every p at least two, define the normalized product vector

```math
\psi_p=\bigotimes_{a=0}^{m/4-1}|0,+,0,0\rangle_{4a,\ldots,4a+3}.
                                                               \qquad(11)
```

This uniformly sums selector strings supported on corner one. Write row
labels as (a,b), with b the low two bits. N_p maps such strings to
$`(z'_a,z'_a,0,0)`$ in each cell, where
$`z'=N_{p-2}z`$. The higher-coordinate transform is invertible,
so it preserves the uniform sum. Therefore

```math
N_p^\dagger\psi_p
 =\bigotimes_a\bigl(|\Phi^+\rangle_{4a,4a+1}|00\rangle_{4a+2,4a+3}\bigr).
                                                               \qquad(12)
```

Each J_j fixes one Bell cell, with literal phase +1. Frame 1's zero
corners 2,3 kill its CCZ. In frame 2, the X_1,Z_0 conjugation converts
the Bell pair to a singlet; H on both controls preserves its absence of
an 11 component. In frame 3, Z_0 makes the Bell pair Phi-minus and the
two Hadamards make Psi-plus, again with no 11 component on controls
0,1. Its X_2 does not change that fact. The reverse walls return every
cell, including their literal Pauli phases. Consequently

```math
V_j\psi_p=\psi_p\qquad(j=1,2,3).                           \qquad(13)
```

The physical s_0 bit is zero. Every marked-origin R_C and its actual
adjoint fixes psi_p, for any matching C. Choose t_0=0 for every distinct
principal t bank; every B_Y then fixes the extended vector as well.
On the empty bank, R_C and B_Y are identity and the same V_j fix psi_p.

Thus **arbitrary repetition, ordering, actual adjoints and signal controls
of these specified generators**, with common signal-only interlayers,
still have identical signal action on a marked and an unmarked input.
Equation (6) gives a quarter Gram floor at
every p at least two. It covers a swap placed inside a Fourier frame as
well as between completed branches. It strengthens the earlier finite
wrapped-frame failure; the old explicit witness and crossing-rank fixture
remain valid. It does not cover arbitrary nonlocal frames or a newly
added generator that moves this Bell-cell code.

If a physical source is within epsilon of one of the excluded initialized
isometries, contraction gives an accepted-Gram change at most
$`2\epsilon`$. Every corresponding quarter floor becomes
$`1/4-2\epsilon`$. This robustness statement includes leakage and
arbitrary external references.

<a id="inside-frame-bank-h-source"></a>

## 6. An inside-frame bank-H source and an exact three-bit Gram floor

An actual change of relative signal action is to put the bank couplings
inside a selector frame, after the completed stage. Consider

```math
A_{j,0}=R_{C_j},\qquad
A_{j,1}=V_j^\dagger D_YB_YR_{C_j}V_j.                    \qquad(14)
```

The one-signal multiplexor is emitted chronologically as

```text
controlled-V_j, R_Cj, controlled-B_Y, controlled-D_Y, controlled-V_j†.
```

The bank gates use the current Y after R; no controlled whole R is
assumed. The actual branch adjoint is
$`V_j^\dagger R_{C_j}^\dagger B_YD_YV_j`$, with the order
reversed. Omitting D_Y gives the swap-only attempt already excluded by
Section 5, even though that placement breaks the earlier two-term witness.
Adding D_Y moves the new Bell-cell witness and requires a new calculation.

For three such multiplexors with actual signal matrices S_0 through S_3,
the accepted block is

```math
K=\sum_{b_1,b_2,b_3\in\{0,1\}}c_bW_b,\qquad
W_b=A_{3,b_3}A_{2,b_2}A_{1,b_1},                         \qquad(15)
```

```math
c_b=(S_3)_{0b_3}(S_2)_{b_3b_2}(S_1)_{b_2b_1}(S_0)_{b_1\,0}.
```

Every empty-bank branch is identity. Arbitrary complex c_b are a useful
screening relaxation of these factored signal coefficients, not a native
linear-combination primitive.

The frame controls use the paid controlled cell-CCZ layer surrounded by
ordinary Clifford walls, including transpose zeta, Pauli conjugations
and partial Hadamards. The three stages, six frame controls, three swap
controls and three bank-H controls have simultaneous bounds

```math
D_T\le3(116)+6(16)+3(16)+3(6)=510,\qquad
T,b=O(pm),\qquad G=O(m^2+pm),                            \qquad(16)
```

before the charged signal synthesis. There is one clean signal, one
arbitrary principal t bank, and reusable returned implementation helpers.
Any signal approximation is charged in the complete initialized-isometry
error. This is a paid circuit family, not a successful lookup source.

At p=3, use C_01,C_12,C_20 and the specified frames of Section 5.
Restrict to x=0,Y=e_0. Every other t wire is an identity factor, but both
t_0 inputs remain: the effective work dimension is 512. Define the
complete marked trace Gram

```math
G_{bc}=\frac1{512}\mathrm{tr}(W_b^\dagger W_c).             \qquad(17)
```

Paths are ordered lexicographically by (b_1,b_2,b_3). Set
$`\Delta=\mathrm{diag}(1,2,2,4,2,4,4,8)`$.
Exact evaluation of all 512 columns, with literal stage factors and
actual adjoints, gives $`G=\Delta^{-1/2}A\Delta^{-1/2}/1024`$,
where the integer matrix is

```math
A=\begin{pmatrix}
1024&32&0&280&32&252&314&82\\
32&2048&556&0&318&64&61&628\\
0&556&2048&64&411&37&64&348\\
280&0&64&4096&22&822&323&128\\
32&318&411&22&2048&64&0&560\\
252&64&37&822&64&4096&1112&0\\
314&61&64&323&0&1112&4096&128\\
82&628&348&128&560&0&128&8192
\end{pmatrix}.                                             \qquad(18)
```

The square roots only encode the number of Hadamards along each path:
$`2^{|b|/2}W_b`$ has rational entries. No floating-point rank threshold
is used. All diagonal entries of G are one. Bounding
$`\sqrt2\le3/2`$ gives the following upper numerators for its eight
off-diagonal row sums, over the common denominator 131072:

```math
(64224,85104,80944,61760,71776,85936,75088,59232).
```

Each row sum is below two thirds. Symmetric diagonal dominance therefore
proves the simple certificate

```math
G\succeq I/3.                                             \qquad(19)
```

In particular, exact marked cancellation forces every c_b to vanish.
There is also a fixed complete error floor, without assuming exact
cancellation. Put $`a=|\sum_b c_b|^2`$, the empty-bank acceptance.
For the marked block, normalized trace is at most operator norm, so

```math
\begin{aligned}
\delta&\ge\|K_{\rm mark}^\dagger K_{\rm mark}\|
 \ge c^\dagger Gc\ge\tfrac13\|c\|_2^2\ge a/24,\\
\delta&\ge|a-1/2|.
\end{aligned}                                                \qquad(20)
```

These inequalities imply the conservative bound $`\delta\ge1/64`$.
It holds for every actual signal choice, and even for the larger set of
arbitrary path coefficients. A physical implementation within epsilon
of this source retains the floor $`1/64-2\epsilon`$.

This certificate concerns the specified p=3 source only. It does not
exclude a family patched at that size, other nonlocal frames, or arbitrary
bank-dependent interleavings. Unlike the all-size result in Section 5,
its scope comes from a complete finite operator calculation.

## 7. Consequence for the next source

Bank-dependent access is affordable and can change relative signal action.
The swap word creates memory, and an inside-frame bank-H moves the known
Bell-cell code. Neither fact establishes the target Gram. For the latter
candidate, the marked path operators instead become linearly independent:
there is no nonzero coefficient choice that cancels their complete action.
Tuning signal angles cannot repair that failure.

The next useful construction should start with an exact marked-sector
operator relation, or a uniform quantitative near-relation within the
recovery budget. It must then give normalized coverage of every unmarked
bank/address sector. Empty-versus-singleton distinction alone is weaker
than addressed lookup, as the Clifford parity calibration shows. Any new
relation must retain arbitrary principal t inputs, literal phases and
actual inverses; the existing recovery and promotion interfaces remain
available once the complete Gram is proved. No improvement to R96 or
unrestricted native-depth lower bound is supplied here.

Exact fixtures: [native bank controls](../../tests/test_bank_coupled_selector_controls.py)
check columns, dirty-helper return, actual positive-S reverses and layers;
[complete swap dynamics](../../tests/test_bank_coupled_selector_dynamics.py)
check every two-bit column and structured larger moving-bank columns;
[Bell-cell and finite Gram fixtures](../../tests/test_bank_coupled_selector_frames.py)
check the code transport, local literal phases, surviving all-size-witness
structure, and the complete three-bit matrix certificate. Finite fixtures
corroborate the analytic statements and establish the explicitly finite
certificate; they do not assert a general successful source.
