# Native Pauli transport and nested bank coefficients

[MOD3 transport](MOD3_ZETA_TRANSPORT.md#pauli-transport-interface) ·
[Working checkpoint](../../WORKSPACE.md) ·
[Route decisions](../ROUTE_DECISIONS.md)

The complete transport contract permits a different source representative
when its relative error commutes with the input Pauli. This freedom gives
small native constructions. One controlled swap transports the specified
Pauli through any one signed MOD3 coordinate, but those words do not
compose into the complete address orbit. A separate three-layer word
below pays the complete two-address bank, including its joint sector,
with no initialized or borrowed work. Its final band commutes with the
residual lookup sign; retaining that condition prevents the same
three-band mechanism from reaching three address bits. A four-band word
also pays two nested bank coefficients, including a quartic term. A
separate two-band coefficient shear completes all eight coefficients in
five source layers, and gives a complete all-size source with linear
address depth. The missing operation is now a shallower joint assembly;
the complete source supplied by that serial recurrence retains the old cost.

<a id="one-coordinate-compression"></a>

## 1. A complete two-layer source for one coordinate

Use the normalized binary coordinates of the
[MOD3 reduction](MOD3_ZETA_TRANSPORT.md#order-three-blocks). There are
$`p\ge1`$ address bits, $`m=2^p`$ arbitrary bank bits in
$`q=m/2`$ pairs, and sign representatives $`\sigma_j\in\{+1,-1\}^p`$
whose first sign is positive. Put the all-positive representative first.
For coordinate $`i\in\{0,\ldots,p-1\}`$, define

```math
F_i=\bigoplus_{j=0}^{q-1}A^{\sigma_{j,i}},\qquad
v=\bigoplus_j(1,0),\qquad w_i=F_iv,\qquad P_0=X^v.
```

Thus each block of w_i is (0,1) at a positive sign and (1,1) at a
negative sign. The pivot block j=0 has positive sign for every i.
Let C_i be the following fixed CNOT word on the bank:

- fan the pivot first bit into every other first bit;
- fan the pivot second bit into every other second bit;
- additionally fan the pivot second bit into every negative-sign first bit.

All controls lie in the unchanged pivot pair and all targets lie outside
it. These CNOTs commute, $`C_i^2=I`$, and, as binary maps,

```math
C_i e_{0,0}=v,\qquad C_i e_{0,1}=w_i.
```

Let J swap the two pivot wires and let $`R_i=C_iJC_i`$. Emit C_i,
a swap of the pivot pair controlled by x_i, and C_i, in that chronological
order. The resulting full unitary E_i has bank action $`R_i^{x_i}`$.
In particular,

```math
E_iP_0E_i^\dagger
 =X^{F_i^{x_i}v}\otimes I_d.                         \qquad(1)
```

Here d is the arbitrary private helper used by the native swap below.
Equation (1) covers every coherent address, bank, helper and reference
input. Other address bits are arbitrary spectators. It agrees with the
complete MOD3 target on the address face where only x_i may be nonzero;
it is a partial transport statement for general p. For p=1 that face is
the entire address space, so it supplies the complete weak source.

The operator action can be written without a circuit diagram. For a
binary bank column z put $`\ell(z)=z_{0,0}\oplus z_{0,1}`$ and
$`u_i=v+w_i`$. Then

```math
R_i z=z+\ell(z)u_i,\qquad \ell(u_i)=0.                \qquad(2)
```

Thus R_i is an involution, its binary displacement rank is one, and it
exchanges v and w_i. The source generally changes bank columns; its
Pauli echo has precisely the full operator action in (1).

### Native schedule, precision and actual inverse

The retained [two-layer dirty-helper Toffoli](DIAGONAL_ENCODING_SESSION_RESULT.md#2-constructive-check-and-first-unpaid-cost)
uses eight T/T-dagger gates, sixteen CNOTs and two H gates. Its one
helper d is arbitrary and returns exactly, including literal phase.
For pivot bits a,b the chronological word

```text
CX(a,b), CCX(x_i,b;a), CX(a,b)
```

is controlled swap. Hence it costs eight T gates, twenty elementary
Cliffords and two physical T layers, with one arbitrary returned helper.
The Clifford walls contain
$`c_i=2(q-1)+n_i`$ CNOTs, where n_i is the number of negative signs.
Here $`n_0=0`$, and $`n_i=q/2`$ for i>0. The common forward and
actual-reverse source ledger is

```math
T_{E_i}=8,\qquad D_T(E_i)=2,\qquad
G_{E_i}=20+2c_i=4q+16+2n_i,\qquad w_{E_i}=1.           \qquad(3)
```

There are zero clean inputs and zero implementation error. The actual
reverse reverses every native gate and changes T to T-dagger and
conversely; no ideal involution is substituted for an emitted inverse.
Both physical T bands act on four distinct wires. Using
$`X=HS^2H`$ on the q fixed first bits, the emitted observable echo has
T-count 16, T-depth at most four, Clifford count $`2G_{E_i}+4q`$,
and the same one returned helper. Signal control of this observable
requires only q CNOTs in place of the middle Pauli; E_i remains
unconditional.

This construction does not contradict the
[order-three spectrum obstruction](MOD3_ZETA_TRANSPORT.md#scoped-mod3-limits).
That obstruction concerns conjugating a full order-three unitary into
copies of itself with identity helpers. The source here matches one
Pauli image for the two values x_i=0,1 and uses a different bank unitary.

<a id="composition-failure"></a>

## 2. Exact one-coordinate identities do not compose

Consider the candidate $`E=E_{p-1}\cdots E_0`$, with the same helper
reused after its exact return. Every source has the correct image of
P_0 individually. After the first source, however, the next source must
act on the already transported Pauli, rather than on P_0 again.
Equation (2) makes this failure explicit. All the R_i commute, since
every u_i is in the kernel of the same linear functional ell. Therefore

```math
R_{p-1}^{x_{p-1}}\cdots R_0^{x_0}v
 =v+\sum_i x_i(v+w_i).                                \qquad(4)
```

This mask is affine in the address. The required mask
$`F(x)v=\prod_iF_i^{x_i}v`$ need not be affine.
For every p at least two, set x_0=x_1=1 and all other address bits zero.
The two masks on a block j are then

| Sign $`\sigma_{j,1}`$ | Candidate mask | Required mask |
| --- | --- | --- |
| +1 | (1,0) | (1,1) |
| -1 | (0,0) | (1,0) |

Both sign classes are present. The candidate and target are distinct
commuting X Paulis, whose operator distance is two. Consequently

```math
\|EP_0E^\dagger-\mathcal B_pP_0\mathcal B_p^\dagger
                    \otimes I_d\|=2\qquad(p\ge2).      \qquad(5)
```

The equality follows from this address block and the general upper bound
two for the difference of unitaries. Private-helper return is already
exact. Additional precision cannot repair this exact sector mismatch.

<a id="prefix-orbit-rank"></a>

## 3. A prefix-orbit rank requirement for fixed coordinate factors

There is a stronger necessary condition than checking one Pauli at each
stage. Consider arbitrary fixed invertible binary maps V_0,...,V_{p-1}
on the bank and any auxiliary bits, with source action
$`V_{p-1}^{x_{p-1}}\cdots V_0^{x_0}`$. Active factors may move all
those bits. The address is the read-only selector of this stated
factorization, and a factor is identity when its selector bit is zero.
Assume the source transports v, with zero auxiliary mask, to F(x)v
with zero auxiliary mask for every address x. Then

```math
\mathrm{rank}(V_i-I)\ge2^i,\qquad
\sum_{i=0}^{p-1}\mathrm{rank}(V_i-I)\ge m-1.             \qquad(6)
```

To prove this, use $`F(x)v=Se_x`$ from the fixed binary reduction.
Setting every later address bit to zero forces the completed prefix
to send v to $`Se_y`$ for each $`0\le y\lt2^i`$. Turning on x_i
then requires

```math
(V_i-I)(Se_y)=S(e_y+e_{y+2^i})\qquad(0\le y\lt2^i),      \qquad(7)
```

with all vectors embedded by zero auxiliary coordinates. These $`2^i`$
differences are independent: their two-element supports are disjoint
before the invertible S wall. This proves (6). No commutation assumption
on the V_i is used. In the phase-free binary-permutation grammar, an
observable error strictly below two forces the same mask identity on
each address, so the same necessary rank bounds apply.

In particular, any binary conjugate of a single swap has displacement
rank one. Choosing a different C_i cannot make p such single-swap
representatives compose into the full transport when p is at least two:
the second coordinate already requires rank at least two. A stage
conjugate to k_i disjoint swaps requires $`k_i\ge2^i`$.

The label argument also extends to fixed completed controlled-Clifford
factors on the bank and arbitrary work. Use their full symplectic
matrices in (6) and embed each target X label in the symplectic space.
The same prefix differences are independent. An observable error
strictly below $`\sqrt2`$ forces equality of the Hermitian Pauli images,
since distinct signed Pauli operators have distance at least sqrt(2).
Literal phases remain an additional obligation. A physical swap has
symplectic displacement rank two, so this extension has a different
active-swap count than the binary restriction.

This is a rank requirement for the specified completed factorization.
It is not a native T-depth lower bound: the retained
[shared-root CCX batch](MOD3_ZETA_TRANSPORT.md#3-a-paid-linear-depth-native-baseline-and-optional-signal-control)
implements arbitrarily many controlled disjoint CNOTs in four T layers,
while the active bank map's displacement rank grows with their number.
The available width permits m directions. General quantum sources need
not have these fixed controlled-Clifford factors or preserve the address
at their boundaries.

<a id="joint-two-address-source"></a>

## 4. A three-layer source for the complete two-address bank

Let a,b be arbitrary address bits and d,t,v,w four arbitrary bank
bits. There is an exact source for the complete cubic normal form

```math
\begin{aligned}
f&=d+av+bw+abt\pmod2,\qquad E Z_d E^\dagger=(-1)^f,\\
T_E&=12,\quad D_T(E)=3,\quad G_E=68,\quad w_E=0.       \qquad(8)
\end{aligned}
```

No clean or dirty helper is supplied. The source may have
input-dependent phases; its actual Pauli echo cancels them. The
forward and actual-inverse bounds in (8) use adjacent S-power
reduction modulo four.

### 4.1. A complete cubic extractor

Put P=Z_d and define

```math
\begin{aligned}
F_{\rm in}&=\exp[-i\pi Y_d(1-Z_t)/8],\\
K_{\rm in}=F_{\rm in}PF_{\rm in}^\dagger
 &=|0\rangle\langle0|_t Z_d
   +|1\rangle\langle1|_t X_d.                       \qquad(9)
\end{aligned}
```

The two independent commuting axes Y_d and Y_d Z_t have opposite
T orientations. A literal chronological word is

```text
S_d^3, H_d, CX(d,t), [T_d, T_t-dagger], CX(d,t), H_d, S_d.
```

The bracket is one physical band: two T gates and eight Cliffords.
Opposite orientations cancel their common scalar phases. Define

```math
q=d-(b\mathbin\oplus d)
 +(a\mathbin\oplus b\mathbin\oplus d)
 -(a\mathbin\oplus d\mathbin\oplus t),\qquad
D=\mathrm{diag}\,e^{i\pi q/4}.                       \qquad(10)
```

The four parities are independent on a,b,d,t. Compute them with
CX(d,b), CX(b,a), CX(a,t), CX(b,t), CX(d,t), apply the
single band [T_d, T_b-dagger, T_a, T_t-dagger], and reverse the
five CNOTs. Thus D costs four T gates and ten CNOTs. Its native
basis temporarily moves the address wires.

At t=0, K_in commutes with D. At t=1 the phase difference on
flipping d is $`4ab+2(-1)^{a+d}`$ modulo eight. Consequently

```math
DK_{\rm in}D^\dagger=(-1)^{abt}K_0,\qquad
K_0=|0\rangle\langle0|_t Z_d
    +|1\rangle\langle1|_t Z_aY_d.                   \qquad(11)
```

For the cubic term alone,
$`F_0=\exp[i\pi X_d Z_a(1-Z_t)/8]`$ has
$`F_0PF_0^\dagger=K_0`$. The source
$`F_0^\dagger D F_{\rm in}`$ is therefore a relative-phase
three-control NOT on all sixteen inputs: it toggles d by abt,
using eight T gates in three bands and twenty-four Cliffords.
This is a construction, with no minimum-depth claim.

### 4.2. Both quadratic terms fit the same final band

Write $`z_j=(-1)^j`$ for a computational Z eigenvalue and put

```math
\begin{aligned}
s&=4av+4bw+4at-2t\pmod8\\
 &=2-2z_a-z_v-z_b-z_w+z_a z_v+z_b z_w+z_a z_t\pmod8,\\
F_{\rm out}&=\exp[-i\pi X_d s/8].                    \qquad(12)
\end{aligned}
```

Since $`e^{i\pi s/4}=(-1)^{av+bw}(-i z_a)^t`$, rotation
of Z_d around X_d gives

```math
K_{\rm out}=F_{\rm out}PF_{\rm out}^\dagger
          =(-1)^{av+bw}K_0.                         \qquad(13)
```

The two even coefficients are Clifford rotations. The six odd axes,
with their signs retained from (12), are

```math
X_dZ_v,\quad X_dZ_b,\quad X_dZ_w,\quad
X_dZ_aZ_v,\quad X_dZ_bZ_w,\quad X_dZ_aZ_t.            \qquad(14)
```

They commute and are physically independent. After H_d their parity
rows are $`d+v`$, $`d+b`$, $`d+w`$, $`d+a+v`$,
$`d+b+w`$ and $`d+a+t`$. From rows
ell_1,...,ell_6 recover $`d=\ell_2+\ell_3+\ell_5`$ and
$`a=\ell_1+\ell_4`$, then the other four coordinates. The
fixture emits their invertible basis by twenty CNOT row operations.
Use its actual reverse around the one six-target T band. No clean
parity wire is used.

The two even rotations have the literal word H_d, S_d, CX(a,d),
S_d^3, CX(a,d), H_d; surround the parity word by two further
H_d gates. F_out and its actual inverse each cost six T gates
and fifty Cliffords after adjacent-S reduction. Three T and three
T-dagger orientations cancel their common scalar phases; the two
Clifford rotations have opposite common scalars as well.

Set $`R=(-1)^{av+bw+abt}`$ and

```math
E=F_{\rm out}^\dagger D F_{\rm in}.                 \qquad(15)
```

F_out commutes with R, so (11)--(13) prove (8) as a full operator
identity. The bands have two, four and six distinct targets. Reversing
every gate and exchanging T with T-dagger pays the actual inverse.
The complete source preserves a,b,t,v,w. Its remaining action is the
toggle of d by av+bw+abt with possible literal relative phases.
No phase is discarded when the source is used in its Pauli echo.

### 4.3. The literal four-row indicator

Order the bank as Y_0,Y_1,Y_2,Y_3 for addresses a+2b. Start with
d=Y_0, v=Y_1, w=Y_2, t=Y_3. The fixed chronological word C

```text
CX(d,v), CX(w,t), CX(d,w), CX(v,t)
```

maps these to $`d=Y_0,v=Y_0+Y_1,w=Y_0+Y_2`$ and
$`t=Y_0+Y_1+Y_2+Y_3`$. It preserves Z_d, and f in (8)
is exactly $`Y_{a+2b}`$ in these coordinates. Put H_B equal
to Hadamard on all four bank wires. Then

```math
E_2=H_B C^\dagger E C H_B,\qquad
E_2 X_{Y_0} E_2^\dagger=I_2,\qquad
T_{E_2}=12,\quad D_T(E_2)=3,\quad G_{E_2}=84,
\quad w_{E_2}=0.                                    \qquad(16)
```

Every address sector, including a=b=1, and every coherent bank and
reference input are covered. Fixed S walls convert (16) to normalized
MOD3 Pauli coordinates if that interface is used. The direct native
echo costs twenty-four T gates, at most six T layers and 172
Cliffords, counting $`X=HS^2H`$. This ledger uses the actual
inverse and makes no unimplemented central-folding saving.

### 4.4. A parameterized joint two-coordinate source

The construction also pays a partial joint operation for every p at
least two. In the normalized MOD3 bank, the four masks

```math
v,\quad F_0v,\quad F_1v,\quad F_1F_0v
```

are independent: they are S applied to four different standard basis
vectors. Complete them, in that order, to a binary bank basis C_p.
For example, append standard coordinate vectors whenever they increase
the span. Binary Gaussian elimination emits C_p and its actual inverse
using at most $`m^2+3m`$ CNOTs each and no helper.

Apply E_2 on the first four bank coordinates and the two selected
address bits, with all other coordinates arbitrary spectators, and set

```math
\begin{aligned}
\widehat E_p&=U_{C_p}E_2U_{C_p}^\dagger,\\
\widehat E_p P_0\widehat E_p^\dagger
 &=\sum_{a,b}|ab\rangle\langle ab|
       \otimes X^{F_1^bF_0^av}\otimes I_{\rm spectator},\\
T_{\widehat E_p}&=12,\qquad D_T(\widehat E_p)=3,\qquad
G_{\widehat E_p}\le84+2(m^2+3m),\qquad w_{\widehat E_p}=0.
                                                               \qquad(16a)
\end{aligned}
```

The spectator identity in the observable refers to the other address
bits and any external work; the displayed Pauli acts on the entire
m-bit bank. This is a complete operator identity for that partial
operation on arbitrary inputs. It matches the full MOD3 target on
the face where every other address bit is zero, and does not claim
the full target when those bits are active. The fixed C_p words and
the native E_2 word give an actual parameterized family, including
its separately reversed source, not just an existence of Pauli labels.

<a id="final-band-commutation-boundary"></a>

## 5. What the final band must change at three address bits

Let P be a fixed Hermitian Pauli and R a Boolean diagonal involution
commuting with P. Suppose

```math
E=A^\dagger B C,\qquad EPE^\dagger=PR,\qquad [A,R]=0,  \qquad(17)
```

where each of A,B,C has at most one native T layer with arbitrary
Clifford boundaries. Then the Boolean degree of R is at most three.
Indeed $`K_C=CPC^\dagger`$ and $`K_A=APA^\dagger`$ are
literal Hermitian Cliffords, and (17) gives

```math
R=B K_C B^\dagger K_A.                              \qquad(18)
```

This is an actual circuit with at most two T layers. The retained
[exact two-layer diagonal rigidity theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md)
makes R diagonal C3, whose literal Boolean signs have degree at
most three. Monomiality of K_A,K_C and a width bound are not assumed.

For lookup use $`P=Z_{Y_0}`$ and $`R=(-1)^{Y_0+Y_x}`$.
Its degree is p+1: each bank variable has a nonzero full-address
monomial, which adding Y_0 does not cancel. Thus (17) is impossible
for p at least three. The retained
[three-layer rounding theorem](MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness)
also excludes error strictly below tau_* in this class, because that
complete approximation would already be exact.

The commutation condition is on the entire final one-band block A,
including its Clifford boundaries, in the stated Pauli coordinates.
This is not an unrestricted three-layer lower bound. Arbitrary
intermediate address and work movement are allowed; sources whose
final block does not commute with R are outside the argument. At
four or more source layers (18) no longer supplies the stated
two-layer theorem. A continuation must pay a noncommuting final
action, or change to a different approximate joint source.

<a id="four-band-nested-coefficient-source"></a>

## 6. A four-band source for two nested bank coefficients

A fourth source band pays a genuine quartic, bank-linear term. Let b,c,t
be three arbitrary address bits and d,a,v three arbitrary bank bits.
There is an exact native source satisfying

```math
E_4 Z_d E_4^\dagger=(-1)^{d+bt(a+cv)},\qquad
T_{E_4}=12,\quad D_T(E_4)=4,\quad G_{E_4}=34,
\quad w_{E_4}=a_{\rm clean}=0.                         \qquad(19)
```

All six inputs and their references are arbitrary. The exponent contains
the constant bank coefficient d, the bt coefficient a and the bct
coefficient v. It is not the complete three-address lookup polynomial.

### Construction and complete action

Use the cubic source $`E_{\rm cub}=F_0^\dagger D F_{\rm in}`$
from Section 4.1, acting on a,b,d,t. Its actual basis permutation toggles
d by abt. Define

```math
R=X_aZ_d,\qquad
G=\exp(i\pi Rcv/2),\qquad E_4=G E_{\rm cub}.           \qquad(20)
```

At cv=0, G is identity; at cv=1 it is iR. Thus G toggles a by cv,
with a retained relative phase, and does not change b,c,d,t,v. The
whole source maps a to a+cv and d to d+abt, where the latter a is
the input value. Its inverse therefore recovers the input d as
$`d+bt(a+cv)`$. Conjugating Z_d gives (19). Relative phases cancel
between the source and its actual adjoint on every input; no helper
sector or basis-state preparation is used.

The extra source is not a serial six-band composition. Put

```math
Q_1=X_dZ_a,\qquad Q_2=X_dZ_aZ_t.
```

The last band of the cubic source is
$`T(Q_1)T(Q_2)^\dagger`$, with
$`T(P)=e^{i\pi/8}e^{-i\pi P/8}`$. Every axis in
$`R,RZ_c,RZ_v,RZ_cZ_v`$ commutes with both Q_1 and Q_2: the
anticommutations at a and d cancel. Consequently its combined final
word is exactly

```math
G F_0^\dagger
=T(Q_1)T(Q_2)^\dagger T(R)^\dagger
 T(RZ_c)T(RZ_v)T(RZ_cZ_v)^\dagger.                    \qquad(21)
```

The three positive and three negative T orientations cancel their
common scalars. Five displayed axes are independent; the last is the
product of R, RZ_c and RZ_v. The chronological Clifford

```text
H_d, CX(a,d), H_a
```

sends Q_1 to Z_d and R to Z_a, while leaving Z_c,Z_t,Z_v fixed.
Follow it by CX(d,t), CX(a,c), CX(a,v). The first five rotations
then use the distinct physical targets d,t,a,c,v with orientations
T, T-dagger, T-dagger, T, T. Compute the parity of the current a,c,v
wires into v using CX(a,v), CX(c,v), apply T-dagger to v, and undo
those two CNOTs. Reverse the three parity CNOTs and the three-gate
Clifford. This pays (21) in two physical bands and sixteen Cliffords.

The retained first two cubic bands cost six T gates and eighteen
Cliffords. Hence the complete source has bands of sizes 2,4,5,1 and
the ledger (19), also for its actual native reverse. The complete
observable echo has twenty-four T gates, depth at most eight and
seventy Cliffords, counting $`Z_d=S_d^2`$. It returns every address
and bank input exactly, with zero implementation error.

### Parameterized partial family and its omitted coefficients

For every p at least three, the same word acts on three selected address
bits and three selected bank coefficients, with all other inputs arbitrary
spectators. In a physical $`m=2^p`$ lookup bank, the fixed binary
Möbius change to all multilinear coefficients costs pm/2 CNOTs. Its
constant coefficient is Y_0. Placing d,a,v at coefficients 0,5,7 for
addresses b,c,t, then conjugating the source by this change and its
actual reverse, gives a literal partial family with
$`T=12,D_T=4,G=34+pm`$ and no helpers. The four-band count is
independent of p; the transported polynomial still contains only the
three stated coefficients. Its phase echo has $`T=24,D_T\le8`$ and
$`G=70+2pm`$, with no helpers and zero complete error for that partial
polynomial.

Already at p=3, the independent coefficients of b,c,t,bc,ct remain
unpaid. For example, at address b=1,c=t=0 with only Y_1=1, the
partial echo is +1 whereas the required lookup sign is -1. Its full
lookup error is therefore two. Adding the omitted coefficient actions
requires a joint native schedule and a new complete identity; no
parallelism follows from the completed diagonal phases commuting.
The construction opens a four-band joint-interference mechanism but
does not move the complete-frame frontier.

The [native fixtures](../../tests/test_joint_lookup_pauli_source.py)
check all 64 source columns, its actual inverse and all 4096 echo entries,
the physical band and count ledger, and the fixed-bank partial lift.
These finite checks support the complete analytic identity above; they
are not an all-p lookup compiler or external mathematical review.

<a id="complete-coefficient-shear"></a>

### Complete coefficient assembly by a two-band shear

The omitted coefficients can be supplied by an actual native word, but
its serial extension does not improve asymptotic depth. Let c be an
arbitrary address bit and let A and B be two disjoint q-bit banks. Define

```math
\begin{aligned}
G_c^{(q)}
 &=\prod_{j=0}^{q-1}\exp(i\pi X_{A_j}cB_j/2),\\
G_c^{(q)}|c,A,B\rangle
 &=i^{c\sum_j B_j}|c,A+cB,B\rangle.
\end{aligned}                                                    \qquad(22)
```

The sum in the scalar exponent is an ordinary integer sum; the bank
addition is over the binary field. Thus (22) is a full monomial unitary,
with a retained relative phase, on arbitrary coherent inputs. It is not
an initialized copy of B, and no helper is used. For each j its four
Pauli rotations are

```math
T(X_{A_j})^\dagger T(X_{A_j}Z_{B_j})
T(X_{A_j}Z_c)T(X_{A_j}Z_cZ_{B_j})^\dagger.          \qquad(23)
```

Their common scalar phases cancel exactly. All the axes in (23), across
all j, commute. Their full list need not be independent, so one physical
band is not asserted. The first two axes for each j form one independent
list of 2q axes; the last two form another. In the second list, the
B_j label forces the coefficient of $`X_{A_j}Z_cZ_{B_j}`$ to vanish
in any relation, and then the A_j label forces the other coefficient
to vanish. Both lists therefore have a literal schedule on distinct
physical wires.

Here is that schedule. Apply H to every A_j. Compute each parity
$`B_j\mathrel{\oplus}=A_j`$, apply T-dagger to all A_j and T to
all B_j in one band, and undo the parity CNOTs. Next compute
$`A_j\mathrel{\oplus}=c`$ for every j, then
$`B_j\mathrel{\oplus}=A_j`$ using the current A_j. Apply T to all
A_j and T-dagger to all B_j in one band. Undo the latter two CNOT
lists in reverse order and apply H to all A_j. Each band has 2q
distinct targets. This emitted word and its actual native reverse have

```math
T_G=4q,\qquad D_T(G)=2,\qquad G_G=8q,\qquad
w_G=a_{\rm clean}=0.                               \qquad(24)
```

The reverse generally differs from the forward unitary by relative
phases, even though their basis permutations agree. Every echo below
uses the actual reverse of the emitted gates.

To complete three-address lookup, use addresses $`x_0,x_1,c`$ and
eight arbitrary ANF coefficients A_0,...,A_3,B_0,...,B_3. Let

```math
f_A=A_0+x_0A_1+x_1A_2+x_0x_1A_3,
\qquad f_B=B_0+x_0B_1+x_1B_2+x_0x_1B_3.            \qquad(25)
```

Use the Section 4 source E on the A bank, so that
$`EZ_{A_0}E^\dagger=(-1)^{f_A}`$, and append the chronological
word $`G_c^{(4)}`$. Its complete Pauli image is

```math
E_3=G_c^{(4)}E,\qquad
E_3Z_{A_0}E_3^\dagger=(-1)^{f_A+c f_B}.             \qquad(26)
```

Indeed conjugating a diagonal sign by (22) substitutes
$`A\mapsto A+cB`$ and cancels the relative phase, giving exactly
(26). All eight coefficients are independent. In particular the five
terms missing from (19) are now paid; no bank-coordinate promise is
being used. The coefficient-coordinate ledger is

```math
T_{E_3}=28,\qquad D_T(E_3)=5,\qquad
G_{E_3}=100,\qquad w_{E_3}=a_{\rm clean}=0.          \qquad(27)
```

The physical bands have sizes 2,4,6,8,8. No cancellation of adjacent
Clifford gates or unimplemented fusion is included in this ledger.
Conjugating by the three-bit Möbius bank word and its reverse adds
24 CNOTs. The raw-bank phase source therefore costs 124 Cliffords;
its complete $`Z_{Y_0}`$ echo costs 56 T gates, at most ten T layers
and 250 Cliffords. It implements $`(-1)^{Y_x}`$ with zero complete
error and returns every address and bank input. Adding bank Hadamards
on both sides gives the corresponding indicator source at 140
Cliffords; its actual $`X_{Y_0}=HS^2H`$ echo costs 284 Cliffords,
with the same T-count, depth bound and zero extra work.

### The all-size recurrence and why this completion is paused

The same construction gives a uniform complete family for every
$`p\ge2`$, with $`m=2^p`$. Work in the full ANF coefficient bank.
Apply the Section 4 source on the first four coefficients and the
first two address bits. For $`k=2,\ldots,p-1`$ in increasing order,
append $`G_{x_k}^{(2^k)}`$ on coefficient pairs
$`(A_j,A_{j+2^k})`$, $`0\le j\lt2^k`$. Every other coefficient
is an arbitrary spectator at that stage. Induction using (26) gives

```math
E_pZ_{A_0}E_p^\dagger
 =(-1)^{\sum_{S\subseteq\{0,\ldots,p-1\}}
              A_S\prod_{i\in S}x_i}.                         \qquad(28)
```

Equivalently, the inverse basis permutation first folds the upper
coefficient half into the lower half when $`x_{p-1}=1`$, then
recurses on that lower half. This accounts for every independent
coefficient, including mixed address sectors. The source may move the
coefficient bank, but its literal Pauli echo restores it. Tensoring
the full operator identity with any reference proves the coherent
return requirement without a basis-state restriction.

The source and its actual reverse have the simultaneous exact ledgers

```math
\begin{aligned}
T_{E_p}&=12+4\sum_{k=2}^{p-1}2^k=4m-4,\\
D_T(E_p)&\le3+2(p-2)=2p-1,\\
G_{E_p}&=68+8\sum_{k=2}^{p-1}2^k=8m+36,\\
w_{E_p}&=a_{\rm clean}=0.
\end{aligned}                                                   \qquad(29)
```

After the fixed Möbius wall and its reverse, the raw-bank phase-source
Clifford count is $`pm+8m+36`$. Its complete phase echo has
$`T=8m-8`$, depth at most $`4p-2`$ and
$`G=2pm+16m+74`$. The indicator source adds 2m bank Hadamards;
its complete indicator echo has
$`G=2pm+20m+76`$. These are simultaneous zero-error bounds with no
initialized or borrowed work. Every inverse is emitted and charged.

The displayed native schedule has two serial bands for each new
address bit, so it retains linear address depth. This is a cost of
this schedule, not a lower bound on every equivalent circuit.
Polynomial-overhead fusion does not turn this linear-depth source
into the required $`o(\log(p+2))`$ source. Merely adding further
coefficients through the same serial recurrence is therefore paused.
The next unpaid operation is a changed joint schedule that combines
growing coefficient blocks at sublogarithmic native depth, or a
different complete source interface. The five missing coefficients
at p=3 are no longer an open identity; reducing this finite example
alone would still not move the complete-frame frontier.

### A distinct all-bank word and the joint scheduling question

There is a different, more redundant source representation. For each
address coordinate i, apply (22) to every coefficient pair
$`(A_S,A_{S\cup\{i\}})`$ with i absent from S, using control x_i.
Call this full-bank relative-phase word $`\Gamma_i`$. Its binary
permutation is $`M_i=I+x_iN_i`$, where

```math
(N_iA)_S=
\begin{cases}
A_{S\cup\{i\}},&i\notin S,\\
0,&i\in S.
\end{cases}
```

Here $`N_i^2=0`$ and $`N_iN_j=N_jN_i`$, so the completed binary
maps M_i commute. The zeroth row of their product is the full polynomial
in (28). Thus any fixed chronological ordering of the actual
$`\Gamma_i`$ supplies another complete Pauli source. In coefficient
coordinates its explicit serial ledger is $`T=2pm`$,
$`G=4pm`$, at most 2p T layers and no helpers. The Möbius walls,
actual inverse and Pauli echo are paid exactly as above. Polynomial
count overhead is allowed by fusion; this schedule's depth still fails.

These are **all-bank** pairs, distinct from the economical prefix pairs
in (29); the prefix binary maps need not commute. Even the all-bank
native unitaries need not commute: at p=2 with both addresses one,
reversing their order changes the phase by
$`(-1)^{A_{\{0,1\}}(A_{\{0\}}+A_{\{1\}})}`$.
Each ordering nevertheless has the same required Pauli image, because
its own relative phases cancel in its actual echo. Reordering completed
binary maps therefore does not justify reordering their physical T gates.
This identifies a concrete remaining operation: a joint native
resynthesis of this whole Pauli image, or a different representative,
with sublogarithmic depth and a complete fixed-tolerance error bound.
It does not assert that such a resynthesis exists or that this
representation is necessary for an unrestricted source.

The [coefficient-shear fixtures](../../tests/test_complete_coefficient_shear.py)
check every native column and actual inverse for one and two pairs,
the two distinct-target bands and full source ledgers, all address
rows through six bits, and exact p=3 echo columns for every address,
every independent raw-bank direction and several dense banks. These
finite checks challenge the native phase and return identities;
(22)--(29) supply the all-size analytic proof. Neither is external
mathematical review.

## 7. Research consequence and evidence

The two-layer component pays weak transport of a single orbit direction.
The unrestricted two-band screen for the complete source, at p at least
three and the stated fixed error threshold, remains unchanged. The
first unsupported composition is replacing each coordinate action by
a representative certified only on the original Pauli. A successful
joint construction must transport the growing prefix orbit or leave this
completed-factor representation.

Pause the product of single-swap representatives. Reopen it only with
paid action on the whole required prefix orbit, or change to a source
whose non-Clifford interference supplies the joint observable directly.
The separate three-layer construction pays the complete two-address bank
and joint action of two selected coordinates at every bank size, with no
helper. Its final-block commutation condition explains why adding parity
terms within that mechanism does not give the whole address orbit.
The four-band construction adds a quartic term while retaining arbitrary
bank inputs. The separate two-band shear completes all three-address
coefficients in five source layers, and its uniform recurrence pays the
whole bank with linear address depth. This resolves the missing finite
coefficient identity but does not improve the known logarithmic-depth
complete indicator. Pause the serial recurrence as a frontier route;
a new growing-block schedule or a different interface is needed.
The scoped three-layer requirements remain representation limits.
The approximate scalar zero-reflection result still rules out
conjunction alone as a growing fixed-accuracy witness.

The [exact fixtures](../../tests/test_weak_pauli_orbit_compression.py)
check the native source, literal Pauli echo, private helper and actual
reverse on every column in small instances. Binary fixtures check all
addresses through p=7 and the independent prefix differences. They
corroborate the analytic identities; no general quantum compiler,
unrestricted depth lower bound or external mathematical review is claimed.
The [joint-source fixtures](../../tests/test_joint_lookup_pauli_source.py)
also check every column of the native cubic source and complete four-row
echo, their actual inverses, and all physical band and Clifford ledgers.
