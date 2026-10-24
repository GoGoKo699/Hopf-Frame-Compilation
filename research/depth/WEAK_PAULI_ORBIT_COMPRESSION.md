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
also pays two nested bank coefficients, including a quartic term, while
five independent coefficients of the three-address bank remain absent.
The general joint-transport operation remains unpaid.

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
f=d+av+bw+abt\pmod2,\qquad
E Z_d E^\dagger=(-1)^f,\qquad
T_E=12,\quad D_T(E)=3,\quad G_E=68,\quad w_E=0.       \qquad(8)
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
rows are $`d+v,d+b,d+w,d+a+v,d+b+w,d+a+t`$. From rows
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

### Parameterized partial family and first unpaid operation

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
bank inputs; its missing five coefficients require a joint completion.
Neither construction nor either scoped requirement moves the
complete-frame resource frontier. The known logarithmic-depth complete
indicator remains the baseline, and the approximate scalar zero-reflection
result still rules out conjunction alone as a growing fixed-accuracy witness.

The [exact fixtures](../../tests/test_weak_pauli_orbit_compression.py)
check the native source, literal Pauli echo, private helper and actual
reverse on every column in small instances. Binary fixtures check all
addresses through p=7 and the independent prefix differences. They
corroborate the analytic identities; no general quantum compiler,
unrestricted depth lower bound or external mathematical review is claimed.
The [joint-source fixtures](../../tests/test_joint_lookup_pauli_source.py)
also check every column of the native cubic source and complete four-row
echo, their actual inverses, and all physical band and Clifford ledgers.
