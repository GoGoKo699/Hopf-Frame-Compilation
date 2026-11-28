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
The signed-decision test in Section 8 uses an alternative nonflat acceptance
condition. Its full-bank trace bound excludes the proposed centered
recovery schedule for the specified three-stage family at every angle.
Section 9 adds an unmarked single-hole sector: the exact polar replacement
has a constant wrong-sign mean there, and arbitrary signal-only walls
cannot give this three-stage architecture a usable single-threshold gap.

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

For the flat Gram target, an exact marked-sector operator relation or a
uniform quantitative near-relation remains necessary. The
[one-clean signed-decision interface](PRESCRIBED_SCALAR_TABLES.md#one-clean-decision-source)
has an alternative sufficient premise: marked acceptance may be nonzero, provided
the complete marked and unmarked spectra lie on opposite sides of one half
with a usable gap. The finite floor above does not exclude that premise.
Sections 8 and 9 test it for a specified scalable word. Both routes retain
arbitrary principal t inputs, literal phases and actual inverses.
Empty-versus-singleton distinction alone is weaker than addressed lookup.
No improvement to R96 or unrestricted native-depth lower bound is supplied.

Exact fixtures: [native bank controls](../../tests/test_bank_coupled_selector_controls.py)
check columns, dirty-helper return, actual positive-S reverses and layers;
[complete swap dynamics](../../tests/test_bank_coupled_selector_dynamics.py)
check every two-bit column and structured larger moving-bank columns;
[Bell-cell and finite Gram fixtures](../../tests/test_bank_coupled_selector_frames.py)
check the code transport, local literal phases, surviving all-size-witness
structure, and the complete three-bit matrix certificate. Finite fixtures
corroborate the analytic statements and establish the explicitly finite
certificate; they do not assert a general successful source.

<a id="full-bank-signed-gap"></a>

## 8. A calibrated signed-decision word has exponentially small possible margin

### The explicit family and its alternative target

For $`p\ge3`$, use three address matchings, each with disjoint pairs,

```math
 C_0=\prod_{\substack{0\le j\le p-2\\j\ {\rm even}}}\mathrm{CX}(j,j+1),
 \qquad
 C_1=\prod_{\substack{0\le j\le p-2\\j\ {\rm odd}}}\mathrm{CX}(j,j+1),
 \qquad C_2=\mathrm{CX}(p-1,0).
```

Every address coordinate is a control in one of them. Use the specified
frames of Section 5 and the paid multiplexors of (14), with
$`A_{j,0}=R_{C_{j-1}}`$ for $`j=1,2,3`$. For one clean signal a, define

```math
 \begin{aligned}
 S_0&=H,&S_1&=R_z(\alpha)H,&S_2&=R_z(\beta)H,
 &S_3&=(S_2S_1S_0)^\dagger,\\
 \mathcal M_j&=\mathrm{diag}(A_{j,0},A_{j,1}),&
 U_{\alpha,\beta}&=S_3\mathcal M_3S_2\mathcal M_2S_1\mathcal M_1S_0 .
 \end{aligned}
```

The last signal word is the actual inverse of the first three signal
words. The empty-bank unitary is therefore exactly identity. All branch
operators commute with Omega: the completed stage relabels x,Y together,
and the frames and bank-controlled couplings leave x,Y unchanged.
For J inserting a=0, put

```math
 K=J^\dagger U_{\alpha,\beta}J,\qquad
 M=J^\dagger U_{\alpha,\beta}^\dagger Z_aU_{\alpha,\beta}J
   =2K^\dagger K-I .
```

The sufficient decision condition is $`[M,\Omega]=0`$ and
$`\Omega M\succeq\gamma I`$. It permits every marked Gram eigenvalue
below $`(1-\gamma)/2`$, rather than requiring marked cancellation.
At fixed margin and fixed requested error, the existing same-signal
sharpening, bank Hadamards and
[one-clean vector promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
would give the complete zero-clean indicator.

Before signal synthesis the source has the ledger (16):
$`D_T\le510`$, $`T,b=O(pm)`$, $`G=O(m^2+pm)`$, one clean signal,
and two principal arbitrary work banks. These singly exponential costs
fit the shared-layer discovery budget. They are not the square-root
source bounds for an independently prescribed classical table.
For a shrinking margin, the charged centered recovery has fixed-tolerance
depth $`O(\gamma^{-1}[d_U+1+\log(1/\gamma)])`$.

There is only one effective angle. Since a signal Z rotation commutes
with each multiplexor,

```math
 U_{\alpha,\beta}
 =R_z(-\alpha)H\mathcal M_3H\mathcal M_2R_z(\alpha)H\mathcal M_1H .
```

In the eight-path notation (15), with chronological bits $`b_1,b_2,b_3`$,

```math
 c_b=\tfrac14e^{i\alpha b_2}(-1)^{b_2(b_1+b_3)},\qquad
 \sum_b c_b=1,\qquad \sum_b|c_b|^2=\tfrac12 .
```

This is an exact cancellation of beta for these ideal Z-axis rotations.
The half-norm identity is specific to this signal family, not arbitrary
calibrated single-qubit walls.

### Complete tensor-swap action on a full bank

Set $`x=0,Y=\boldsymbol1`$. This entire sector is invariant under every
completed factor and has lookup charge minus one. Both s and t remain
arbitrary. Put $`D=2^m`$, $`k=m/4=2^{p-2}`$, let
$`r_j=R_{C_{j-1}}|_{x=0,Y=\boldsymbol1}`$, $`v_j=V_j`$,
$`h=H^{\otimes m}`$, and let F swap the two work banks.
Each $`r_j`$ is a literal monomial unitary, including its quarter phases.
The frames are real Hermitian involutions. Define
$`\ell_j=v_jh`$ and $`k_j=r_jv_j`$. In tensor order s,t,

```math
 a_{j,0}=r_j\otimes I,\qquad
 a_{j,1}=(\ell_j\otimes k_j)F .
```

Thus every complete path is $`W_b=(X_b\otimes Y_b)F^{|b|}`$:

| b | $`X_b`$ | $`Y_b`$ |
| --- | --- | --- |
| 000 | $`r_3r_2r_1`$ | I |
| 001 | $`\ell_3`$ | $`k_3r_2r_1`$ |
| 010 | $`r_3\ell_2`$ | $`k_2r_1`$ |
| 011 | $`\ell_3k_2r_1`$ | $`k_3\ell_2`$ |
| 100 | $`r_3r_2\ell_1`$ | $`k_1`$ |
| 101 | $`\ell_3k_1`$ | $`k_3r_2\ell_1`$ |
| 110 | $`r_3\ell_2k_1`$ | $`k_2\ell_1`$ |
| 111 | $`\ell_3k_2\ell_1`$ | $`k_3\ell_2k_1`$ |

Write $`\tau_D=D^{-1}\mathrm{Tr}`$ and
$`G_{bc}=D^{-2}\mathrm{Tr}(W_b^\dagger W_c)`$. Equal-parity paths obey

```math
 G_{bc}=\tau_D(X_b^\dagger X_c)\tau_D(Y_b^\dagger Y_c).
```

For opposite parity, $`\mathrm{Tr}((A\otimes B)F)=\mathrm{Tr}(AB)`$
and unitarity give $`|G_{bc}|\le D^{-1}`$. This retains every work column,
not an average over a promised subset of t inputs.

### Cell bounds for every equal-parity pair

The remaining estimates use the explicit frames, not arbitrary mixers.
Let $`P_{j,\pm}^{X,Z}`$ denote the indicated one-qubit Pauli projectors.
Direct conjugation of each literal cell CCZ gives

```math
 \begin{aligned}
 J_1&=I-2P_{0,+}^X P_{2,-}^Z P_{3,-}^Z,\\
 J_2&=I-2P_{0,+}^X P_{1,-}^X P_{3,-}^X,\\
 J_3&=I-2P_{0,+}^X P_{1,-}^X P_{2,+}^Z .
 \end{aligned}
```

In particular J1 is monomial. J2 has largest entry magnitude 3/4 and
entrywise row and column sum 5/2 on every row and column.
J3 has row sums 1 on eight rows and 2 on eight rows, with row maxima
respectively 1 and 1/2. Its column data are the same.
The products J1J2 and J1J3 retain these data because J1 is monomial.
Also

```math
 J_2J_3=I-2P_{0,+}^X P_{1,-}^X(I+X_3Z_2)/2,
```

so J2J3 has the same row and column bounds as J2. Consequently each
distinct pair product $`J_iJ_j`$ has average row maximum and average
column maximum 3/4, and maximum row and column sum at most 5/2.
For the pair 1,3 the latter bound improves to 2.

One further product is needed. Conjugating the first cell projector by
J2 gives

```math
 J_{212}:=J_2J_1J_2
 =I-2P_{0,+}^X P_{2,-}^Z(I-X_1Z_3)/2 .
```

It is identity on $`Z_2=+1`$. On $`Z_2=-1,Z_3=z`$ its other factor is
$`(I-X_0+zX_1+zX_0X_1)/2`$, with four entries of magnitude one half
per row. Hence J212 has the same squared-row-sum mean 5/2 as J3.
Their active Z2 sectors are disjoint, so J212J3 has four entries of
magnitude one half on every row and column.

These cell statements tensorize through
$`v_j=N_pJ_j^{\otimes k}N_p^\dagger`$; the zeta matrix is monomial.
For any monomial r, a normalized trace of $`rA`$ is bounded by the
average column maximum of A. A normalized trace against the flat h is
bounded by $`D^{-3/2}`$ times the entrywise sum of its other factor.
Thus distinct frame pairs give, respectively, bounds
$`(3/4)^k`$ and $`(5/8)^k`$; J1J3 or a monomial-surrounded J3 gives
$`(1/2)^k`$ in the latter estimate.

For completeness, put $`\rho=(3/4)^k`$ and $`\kappa=(5/8)^k`$.
Substitution in the path table, cyclicity of trace and the preceding
cell estimates give all twelve equal-parity off-diagonal bounds:

| Pairs b,c | Bound on $`\lvert G_{bc}\rvert`$ |
| --- | --- |
| 001,010; 001,100; 010,100; 011,101; 101,110 | $`\rho`$ |
| 000,011; 000,110; 001,111; 100,111 | $`\kappa`$ |
| 000,101; 011,110 | $`(1/2)^k`$ |
| 010,111 | $`\kappa`$ |

The last pair needs more than a maximum row-sum estimate. One of its
trace factors cycles to

```math
 \tau_D\!\left(h r_2(v_2v_1v_2)r_3^\dagger v_3\right).
```

The column sums of $`v_2v_1v_2`$ and row sums of v3 both have squared
sum $`D(5/2)^k`$. Monomial factors only permute these sums.
Cauchy--Schwarz bounds the entrywise sum of their product by
$`D(5/2)^k`$; flatness of h then gives $`\kappa`$.
For 011,110, one factor is, up to complex conjugation,
$`\tau_D(r_3^\dagger r_2v_2v_1v_2v_3)`$; the J212J3 row maxima give
the stated $`(1/2)^k`$ bound.

### Consequence for the possible margin and the paid recovery

Every row of $`G-I`$ has three equal-parity entries of magnitude at most
rho and four opposite-parity entries of magnitude at most $`D^{-1}`$.
Therefore, for every real alpha and beta in the specified signal family,

```math
 \begin{aligned}
 |\tau_{D^2}(M)|
 &=|2c^\dagger(G-I)c|
 \le 3(3/4)^k+4/D
 =:\eta_p
 \le7(3/4)^{2^{p-2}} .
 \end{aligned}
```

If the decision premise holds, M is negative definite on this full-bank
sector. It follows that $`\gamma\le-\tau_{D^2}(M)\le\eta_p`$.
The estimate does not assert that every angle has a wrong-sign eigenvalue:
an exponentially small negative margin is not excluded by a small trace.
For an emitted source within initialized-isometry error epsilon of this
ideal word, the compressed observable changes by at most $`2\epsilon`$,
giving $`\gamma\le\eta_p+2\epsilon`$. In particular a certified margin
with $`\epsilon\le\gamma/4`$ still satisfies $`\gamma\le2\eta_p`$.
A larger gap from a changed emitted signal word would need its own proof.
An actual final inverse still calibrates the empty bank, but beta cancels
exactly only for the stated Z-axis words.

Even a bespoke phase list cannot repair the ideal family through the
retained **centered odd-polynomial** recovery at singly exponential source
costs. Either M has a nonnegative eigenvalue and violates the sign premise,
or its negative spectrum contains $`\lambda\in[-\eta_p,0)`$.
On the corresponding two-reflection Jordan block choose coordinates
with physical signal phase axis Z and

```math
 W_\lambda=\lambda I+i\sqrt{1-\lambda^2}X,\qquad
 \|W_\lambda-W_0\|
 =2\sin\!\left(\tfrac12\arcsin|\lambda|\right)
 \le2|\lambda| \quad (|\lambda|\le1/2).
```

A fixed phase-list word with q occurrences of W or its actual inverse
differs between lambda and zero by at most $`2q|\lambda|`$, by telescoping.
The real polynomial of the centered odd filter satisfies $`f(0)=0`$.
Fixed-error approximation to the required minus sign needs
$`|f(\lambda)|\ge1/2`$, and hence, for all sufficiently large p,

```math
 q\ge(4\eta_p)^{-1}=2^{\Omega(m)} .
```

This violates the allowed $`2^{O(p)}`$ emitted-call ledger, as well as the
sublogarithmic depth target, for this recovery schedule. It is not a
native lower bound after arbitrary joint resynthesis of those calls.
This full-bank estimate alone leaves shifted or nonodd recovery outside
its conclusion: for example,
the abstract two-point spectrum $`+1,-\eta`$ is separated to fixed accuracy
by $`2x^2-1`$. Such a replacement here would need a new complete separation
argument for every bank, including all unmarked spectra and dirty inputs.
The general frame frontier and unrestricted joint construction remain open.

The next retained boundary test uses the same eight work paths:
take the exact native word $`W=SHTHS^\dagger`$ and set
$`S_0=S_1=S_2=W,\ S_3=W^{-3}`$. Up to its common phase, W is
$`R_y(\pi/4)`$; all common phases cancel in the complete source.
The actual three forward and three reversed W words add six T gates,
constant Clifford cost and no work. Empty-bank action is still exactly I.
For identical $`R_y(\theta)`$ walls and the actual inverse calibration,
the squared path-coefficient sum obeys
$`2\|\widetilde c\|_2^2-1=\cos^3\theta\cos(3\theta)`$:
insert a computational dephasing after each of the first three signal
walls and propagate the signal's Z component. At $`\theta=\pi/4`$ this
gives $`\|\widetilde c\|_2^2=3/8`$, and the same path Gram bound yields

```math
 \left|\tau_{D^2}(\widetilde M)+\tfrac14\right|
 \le\tfrac34\eta_p .
```

This fixed choice leaves the balanced-coefficient obstruction. It supplies
only a negative full-bank mean, not a signed operator gap.
[Section 9](#single-hole-signed-gap) resolves its complete decision test
with an unmarked dense-bank obstruction. No new selector mixer or recovery
interface is assumed.

The bounded fixtures check the exact cell identities and statistics,
complete small tensor-swap paths and actual inverses, and signal-angle
cancellation, including the exact one-qubit coefficient calibration above.
The tensorization, all-size trace bound and recovery obstruction are analytic.

<a id="single-hole-signed-gap"></a>

## 9. A single hole obstructs the polar source and arbitrary signal-wall tuning

Keep exactly the three matchings, frames and completed branches of
Section 8. The signal walls may now be any four bank-independent
single-qubit unitaries, including actual native words. Two invariant
sectors have opposite lookup charges but nearly the same mean compressed
observable. This restricts one ordered spectral threshold; it does not
assert identical spectra or identical complete source action.

### An invariant unmarked sector with every remaining work input retained

Set $`x=0,Y=\boldsymbol1-e_0`$ and restrict the physical selector bit
to $`X_{s_0}=-1`$. The lookup charge is plus one. Every matching fixes
row zero and permutes the other rows. Moreover, $`h_Y`$ omits $`s_0`$ and
$`d_C(0,0)=0`$. Thus each completed stage is
$`I_{s_0}\otimes r'_j`$, with $`r'_j`$ a literal monomial unitary on
the other selector bits. Both $`B_Y`$ and $`D_Y`$ ignore $`s_0,t_0`$.

The frames also preserve this sector. The transpose-zeta action is

```math
 (N_ps)_u=\bigoplus_{v\supseteq u}s_v,\qquad
 N_p(s_0,z)=(s_0\oplus\ell(z),N_{\rm tail}z),\qquad
 \ell(z)=\bigoplus_{v\ne0}z_v .
```

Here z comprises all nonzero-row selector coordinates. The tail matrix
is invertible, and the restriction to the negative X sector is the
real signed permutation

```math
 N_-|z\rangle=(-1)^{\ell(z)}|N_{\rm tail}z\rangle .
```

In particular $`N_p`$ commutes with physical $`X_{s_0}`$. Each local frame in
Section 8 commutes with its corner-zero X, and all its nonidentity terms
contain $`P_{0,+}^X`$. On the origin cell, restriction to $`X_{s_0}=-1`$
therefore leaves identity on the other three cell qubits. With
$`k=m/4`$, the restricted frame is exactly

```math
 v'_j=N_-\left(I_8\otimes J_j^{\otimes(k-1)}\right)N_-^\dagger .
```

Every $`t_0`$ input is an identity factor throughout the source. Remove that
spectator only from the notation, and write s',t' for the remaining
$`m-1`$ qubits in each bank. With

```math
 D'=2^{m-1}=8\,16^{k-1},\qquad
 h'=H^{\otimes(m-1)},
```

the partial bank swap is the full swap F' of s' and t'. Consequently
the complete restricted branches are

```math
 a'_{j,0}=r'_j\otimes I,\qquad
 a'_{j,1}=(v'_jh'\otimes r'_jv'_j)F' .
```

These are the same tensor-swap paths as in Section 8, with $`k-1`$ active
cells and an $`I_8`$ factor. No selector or work promise has been added to
the source contract: this is an invariant family of legal inputs on
which a complete decision source must also be correct.

### The single-hole Gram bound

The signed permutation $`N_-`$ preserves every entrywise bound used above.
The $`I_8`$ factor changes the mean row-maximum estimate to
$`\rho'=(3/4)^{k-1}`$. In a trace against the flat h', its dimension
contributes the additional factor $`1/\sqrt8`$, giving

```math
 \kappa'=\frac1{\sqrt8}(5/8)^{k-1}\le\rho' .
```

For the exceptional 010,111 pair, the squared row- and column-sum
totals are now $`D'(5/2)^{k-1}`$; the same Cauchy--Schwarz proof gives
$`\kappa'`$. The $`J_{212}J_3`$ bound becomes $`(1/2)^{k-1}\le\rho'`$.
All other equal-parity estimates follow from the same pair products.
Opposite-parity entries are at most $`1/D'`$. Hence the complete
single-hole path Gram, normalized over all s',t' columns, obeys

```math
 \|G_{\rm hole}^- -I\|
 \le 3(3/4)^{k-1}+4/D'
 =:\eta'_p .
```

The arbitrary spectator $`t_0`$ does not change this bound or any mean below.
The full-bank result of Section 8 remains
$`\|G_{\rm full}-I\|\le\eta_p`$ with
$`\eta_p=3(3/4)^k+4/2^m`$.

### A constant wrong sign for the exact polar word

For $`S_0=S_1=S_2=SHTHS^\dagger`$ and the actual inverse
$`S_3=(S_2S_1S_0)^\dagger`$, the squared coefficient norm is 3/8.
Let $`M_{\rm hole}^-`$ be the compressed observable on the preceding
unmarked sector. The new Gram estimate gives

```math
 \left|\tau_{D'^2}(M_{\rm hole}^-)+\tfrac14\right|
 \le\tfrac34\eta'_p .
```

For $`p\ge6`$, we have $`k-1\ge15`$,
$`(3/4)^{15}\lt1/64`$ and $`4/D'\lt1/64`$. Therefore

```math
 \eta'_p\lt1/16,\qquad
 \tau_{D'^2}(M_{\rm hole}^-)\lt-\frac{13}{64}.
```

This legal unmarked sector requires a positive observable. Its mean is
instead negative by a constant, so the polar source fails the
signed-decision premise at every p at least six. A witness density
operator uses $`|-\rangle\langle-|_{s_0}`$, maximally mixed s',t',
and any $`t_0`$ state; it can also be purified into an external reference.
This use of a particular input proves failure on the complete contract,
without giving the construction any initialized work.

### Arbitrary physical signal walls and one ordered threshold

For any four signal-only unitaries, let c be their factored coefficients
in (15) and put $`\nu=\sum_b|c_b|^2`$. Inserting computational dephasing
after each of the first three signal walls makes this sum exactly the
final signal-zero probability. Thus $`0\le\nu\le1`$; no empty-bank
calibration or ideal-angle assumption is needed. Since the walls are
bank-independent, the same c applies to both sectors. Write their
normalized observable means as $`\mu_{\rm full}`$ and
$`\mu_{\rm hole}`$. Then

```math
 \begin{aligned}
 |\mu_{\rm full}-(2\nu-1)|&\le2\nu\eta_p,\\
 |\mu_{\rm hole}-(2\nu-1)|&\le2\nu\eta'_p,\\
 |\mu_{\rm hole}-\mu_{\rm full}|
 &\le2\nu(\eta_p+\eta'_p).
 \end{aligned}
```

In particular, any centered decision margin must satisfy
$`\gamma\le\nu(\eta_p+\eta'_p)\le\eta_p+\eta'_p`$.
The same bound holds if all unmarked spectra lie at least gamma above
one scalar threshold shared across banks and all marked spectra lie at
least gamma below it. Reversing the ordering gives the same conclusion.

There is also an emitted-query obstruction for recovery through such an
ordered separation. In the displayed ordering choose

```math
 \lambda_+=\min\sigma(M_{\rm hole}^-),\qquad
 \lambda_-=\max\sigma(M_{\rm full}).
```

If the ordered separation holds, these opposite-charge eigenvalues obey

```math
 0\lt\delta:=\lambda_+-\lambda_-
 \le\mu_{\rm hole}-\mu_{\rm full}
 \le2\nu(\eta_p+\eta'_p).
```

On the two-reflection Jordan blocks, choose the common physical phase
axis Z and write
$`\mathcal W_\lambda=\lambda I+i\sqrt{1-\lambda^2}X`$.
The upper-semicircle chord estimate, including the endpoint limits, is

```math
 \|\mathcal W_\lambda-\mathcal W_\mu\|^2
 =2\left(1-\lambda\mu-\sqrt{1-\lambda^2}\sqrt{1-\mu^2}\right)
 \le2|\lambda-\mu|.
```

A fixed list of physical Z phases with q queries to this step or its
actual inverse therefore differs between the two eigenvalues by at most
$`q\sqrt{2\delta}`$. Initialized error at most one half from the
required opposite signs forces its two initialized matrix elements to
differ by at least one, even allowing one common global phase.
Consequently, whenever the ordered separation exists,

```math
 q\ge\frac1{2\sqrt{\nu(\eta_p+\eta'_p)}}
 \ge\frac1{2\sqrt{\eta_p+\eta'_p}}
 =2^{\Omega(m)} .
```

If $`\nu=0`$, both observables are already minus identity and the
ordered separation fails. This query bound does not assume an odd
polynomial or a threshold at zero. It concerns one phase list shared
across banks in the same-signal two-reflection recovery. Its emitted
source-call count exceeds the allowed $`2^{O(p)}`$ discovery budget.

### Native and scope consequences

No new gate or state preparation is assumed in this obstruction. The
original branch ledger (16) is unchanged, and the polar word uses the
already charged six extra T gates and actual signal inverses. Arbitrary
native single-qubit signal words satisfy the arbitrary-wall result
directly, including their literal phases; they need not approximate the
balanced or polar choices. If a further implementation is within
initialized-isometry error epsilon of a source considered here, each
compressed mean changes by at most $`2\epsilon`$. Returned implementation
helpers remain identity factors throughout this calculation.

The result stops signal-wall tuning of these specified three work stages
as a source for centered or single ordered-threshold recovery. It does
not prove coincident opposite-charge spectra, rule out a nonmonotone
spectral classification, or constrain arbitrary processing of the raw U.
Nor is it a native lower bound after arbitrary joint resynthesis of
recovery calls. Different work interlayers, frames, stage counts or
bank-dependent signal controls require a new complete calculation.

The [single-hole fixture](../../tests/test_bank_coupled_selector_frames.py)
checks the signed zeta restriction, all reduced three-bit stage columns and
actual inverses, and the complete tensor path traces with exact native polar
coefficients. Its unmarked mean is exactly negative. The all-size Gram bounds,
constant polar failure and arbitrary-wall recovery cost are analytic.
