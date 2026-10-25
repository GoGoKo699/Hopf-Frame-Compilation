# Coupled order interference and anchored reflection trials

[Canonical paired target](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md) ·
[Paid relative phase](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md) ·
[Two-clean Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout)

This note gives a genuinely coupled two-signal source with one joint
acceptance. Its normalization is exact and its Gram is controlled by the
square of a product of two reflections. Two explicit non-Clifford,
bank-moving reflection pairs have constant native T-depth and fit the
polynomial intermediate budgets. Both pass the entire zero-address
sector, but fail a nonzero-address fixture. They are not qualifying
complete indicators. The all-size uncertainty is the squared-product
identity below, not extraction or clean-removal constants.

## 1. The joint order word and complete accepted operator

Let R and S be Hermitian involutions on the entire carrier, including any
address and bank registers. Their private implementation helpers return
on arbitrary inputs. Introduce two clean signal bits a,b. The following
list is chronological:

```math
\begin{gathered}
H_aH_b,\quad \mathrm C_a(R),\\
X_b,\ \mathrm{CX}(a,b),\ \mathrm C_b(S),\
\mathrm{CX}(a,b),\ X_b,\\
X_a,\ \mathrm C_a(R),\ X_a,\quad
\mathrm{CX}(a,b),\ H_a.
\end{gathered}                                                   \qquad(1)
```

Both signals start at zero and only their joint zero output is accepted.
Every displayed X means the literal four-Clifford word HSSH. Before the
last CNOT and H, the four signal paths carry the following words;
operator products act rightmost first:

| Signal path ab | Carrier word |
| --- | --- |
| 00 | RS |
| 01 | R |
| 10 | R |
| 11 | SR |

The middle S is controlled by $`1\mathbin\oplus a\mathbin\oplus b`$;
the final R is controlled by $`1\mathbin\oplus a`$. The final gates
decode the Bell bra $`(\langle00|+\langle11|)/\sqrt2`$. Thus, with
$`J|\psi\rangle=|\psi\rangle|00\rangle`$,

```math
K_+=J^\dagger UJ=\frac{RS+SR}{2\sqrt2},\qquad
K_+^\dagger K_+=\frac{I+\mathrm{Re}[(RS)^2]}4.
                                                               \qquad(2)
```

Inserting Z_a immediately before the Bell decoder instead gives
$`K_-=(RS-SR)/(2\sqrt2)`$ and Gram
$`[I-\mathrm{Re}((RS)^2)]/4`$. The emitted word couples the two signals
before and after the queries. Its accepted polynomial nevertheless
admits an alternative factorization: it is the one-signal order block
$`(RS+SR)/2`$ attenuated by $`1/\sqrt2`$. Thus physical signal coupling
does not establish an irreducible two-signal acceptance mechanism. The
useful change from the old spectral source is the squared-product
interface, not an obstruction to this alternative implementation.

If the common forward/actual-inverse bounds of controlled R and S are
$`(T_R,G_R,d_R,h_R)`$ and $`(T_S,G_S,d_S,h_S)`$, respectively, then

```math
\begin{aligned}
T_U&\le2T_R+T_S,&G_U&\le2G_R+G_S+22,\\
d_U&\le2d_R+d_S,&h_U&\le\max(h_R,h_S).
\end{aligned}                                                   \qquad(3)
```

The 22 signal Cliffords are sixteen from the four literal X words,
three H gates and three CNOTs. The minus variant adds two S gates.
The source has two initialized signals; they need not return before
acceptance. Every implementation helper does return. All identities
are operator identities, so include arbitrary coherent address, bank,
dirty-work and reference inputs; no basis-input promise is used.

<a id="squared-product-interface"></a>

## 2. The squared-product interface already has zero clean work

In the retained standard dual/carrier coordinates, write the desired
paired reflection as Q_x; the accepted target is $`(I+Q_x)/4`$. It is the
conjugate of $`\Pi_{\rm eq}/2`$ under the paid fixed basis envelope.
For any unitary W and Hermitian involution Q,

```math
\mathrm{Re}W=Q\quad\Longleftrightarrow\quad W=Q.
                                                               \qquad(4)
```

Indeed, A=Re W and B=Im W are commuting Hermitian operators and
$`A^2+B^2=I`$. If A=Q, then B squared is zero and B=0. Moreover,

```math
\left\|\frac{I+\mathrm{Re}W}4-\frac{I+Q}4\right\|\le\delta
\quad\Longrightarrow\quad
\|W-Q\|\le4\delta+\sqrt{8\delta}.
                                                               \qquad(5)
```

To prove (5), $`\|A-Q\|\le4\delta`$ and

```math
\|I-A^2\|\le\|Q(Q-A)+(Q-A)A\|\le8\delta,
```

since $`\|A\|\le1`$. Hence $`\|B\|\le\sqrt{8\delta}`$.
This proof does not assume that Q commutes with W.

Consequently the precise missing identity for (1) is

```math
(R_xS_x)^2\approx Q_x.                                      \qquad(6)
```

There is a useful necessary sector condition specific to this interface.
For W=(RS) squared, $`RWR=SWS=W^\dagger`$. Hence exact equality
W=Q implies $`RQR=SQS=Q`$ and $`SRS=RQ`$: the two reflections
commute on the positive Q sector and anticommute on the negative Q
sector. Conversely these relations imply W=Q. If
$`\|W-Q\|\le\varepsilon`$, then
$`\|[R,Q]\|,\|[S,Q]\|\le2\varepsilon`$, by inserting W in
$`RQR-Q`$ and its S analogue. Thus candidate reflections must preserve
the correct lookup sectors approximately, not merely have the correct
zero-address spectra. This is a screen for (6), not a general native
depth obstruction.

There is also a direct Gram test without the square-root loss in (5).
For $`G=(I+\mathrm{Re}W)/4`$, the same conjugation identities give
$`[G,R]=[G,S]=0`$. If $`\|G-(I+Q)/4\|\le\delta`$, then

```math
\|Q-RQR\|\le8\delta,\qquad \|Q-SQS\|\le8\delta.
```

Indeed compare the target with its R conjugate through G; each of
the two differences has norm at most delta. The S proof is identical.
This is a necessary condition on this order-interference Gram, not
a restriction on arbitrary normalized sources.

If solved with error controlled by (5), the four-reflection word in
(6) itself is a complete zero-clean implementation. The two-clean
Gram recovery theorem is available for the source (1), but is
unnecessary for this particular interface. The fixed physical-basis
envelope costs $`4c_n`$ CNOTs, where $`c_n=n^2+3n`$; the carrier XOR
wall costs another 2n CNOTs, and converting the m selected phase outputs
to XOR costs 2m H gates. Thus the direct complete indicator has the same
T-count and T-depth as (6), additional Clifford count at most
$`4c_n+2n+2m`$, and dirty width $`2n-m+h_W`$. Its actual inverse is
the reversed native word, not a new oracle. These are the same
[paired-query transfer operations](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md)
used by the Gram route.

For each fixed requested tolerance tau, a source Gram estimate
$`\delta=O(\tau^2)`$ suffices for direct complete accuracy by (5).
For example, when $`0\lt\tau\le1`$, $`\delta\le\tau^2/4096`$
gives direct error below $`\tau/8`$, suitable for the eight-call
[polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md).
Neither finite fixture below provides such a source estimate.

<a id="balanced-two-unitary-gram-readout"></a>

### A balanced two-unitary Gram source also gives a direct square query

This distinct reduction applies to any actual zero-clean unitary branches
$`V_0,V_1`$ on the complete address, bank and arbitrary-work space. They
need not be reflections, commute, preserve the address, or return their
principal selector registers separately. Any additional private implementation
helpers must return, so that these branches are complete unitaries on their
common principal registers. Suppose an accepted block has the balanced form

```math
K=(V_0-V_1)/2,\qquad Q=V_0^\dagger V_1,\qquad
\|K^\dagger K-P/2\|\le\delta,
```

where P is a projector, including identity on every work register that the
final query must return. Put $`A=\mathrm{Re}Q`$ and $`B=I-P`$.
The Gram is $`(I-A)/2`$, hence $`\|A-B\|\le2\delta`$. The exact identities

```math
\begin{aligned}
Q^2-(2B-I)&=2Q(A-B)+2(Q-I)B,\\
[(Q-I)B]^\dagger[(Q-I)B]&=2B(I-A)B
\end{aligned}
```

give $`\|(Q-I)B\|\le2\sqrt\delta`$ and therefore

```math
\|Q^2-(I-2P)\|\le4\delta+4\sqrt\delta.                 \qquad(6a)
```

No commutation between Q and P is used. For
$`P=(I+\Omega)/2`$, the literal word $`-Q^2`$ is a complete
zero-clean approximation to the lookup reflection $`\Omega`$. It uses
two actual $`V_0^\dagger`$ calls and two $`V_1`$ calls. The literal minus
is $`(XZ)^2=-I`$ on any existing principal wire: with X=HSSH and Z=SS,
it costs twelve elementary Cliffords and no T gates. If the branches have
common forward/inverse ledgers $`(T_j,G_j,d_j,h_j)`$, the direct word obeys

```math
\begin{aligned}
T_{\rm out}&\le2(T_0+T_1),&
G_{\rm out}&\le2(G_0+G_1)+12,\\
d_{\rm out}&\le2(d_0+d_1),&
h_{\rm out}&\le\max(h_0,h_1).
\end{aligned}
```

The private-pool maximum requires complete return of each private pool;
all principal selectors remain counted as arbitrary work. A direct
single-bank phase query becomes the XOR indicator by two m-Hadamard walls.
In paired coordinates the previously stated carrier/basis envelope must
instead be added explicitly. Every final return error is included in (6a).
For $`0\lt\tau\le1`$, requesting $`\delta\le\tau^2/4096`$ gives error
at most $`65\tau/1024\lt\tau/8`$, suitable for eight-call fusion.

The square-root dependence cannot generally be replaced by a linear bound:
for scalar P=0 and $`Q=e^{i\theta}`$, the Gram defect is
$`\delta=(1-\cos\theta)/2`$ and
$`|Q^2-1|=4\sqrt{\delta(1-\delta)}`$.
The squared-query condition alone is not a converse: Q=-I, P=0 has the
correct square and the wrong Gram. Exact equivalence additionally requires
$`Q|_{P^\perp}=I`$.

Thus a successful balanced two-branch source already solves a direct
zero-clean query; its nonunitary accepted block does not relax that
obligation. Different coefficients, attenuated blocks and genuine
multipath acceptance are outside this reduction. The
[coupled quarter-phase source](COUPLED_QUARTER_PHASE_SOURCE.md) is a
fully paid example whose complete Gram fails. Exact native and noncommuting
algebra fixtures are in [the order-source tests](../../tests/test_coupled_order_gram_source.py).

## 3. A paid anchor with the correct zero-address square

Retain p at least one, $`m=2^p`$, n=pm, and standard dual/carrier bits
z,t from the canonical source. Let $`P_0=Z_{z_0}Z_{t_0}`$. The trivial
anchors R=I, S=P_0 do not work in (2), since their product squares to I.
Instead encode $`l=z_0,c=z_0\mathbin\oplus t_0`$ with
$`C_0=\mathrm{CX}(z_0,t_0)`$ and use

```math
\begin{aligned}
r_0&=Z_l,\\
s_0&=|0\rangle\langle0|_c Z_l+|1\rangle\langle1|_c X_l
     =S_c^\dagger\,\mathrm{CY}(c,l)\,Z_l,\\
R_0&=Z_{z_0},\\
S_0&=\frac{I+P_0}2 Z_{z_0}
       +\frac{I-P_0}2 X_{z_0}X_{t_0}.
\end{aligned}                                                   \qquad(7)
```

Both are Hermitian Clifford involutions, and
$`(r_0s_0)^2=Z_c`$, hence $`(R_0S_0)^2=P_0`$ including its literal
scalar phase. The encoded s_0 has a ten-Clifford word: two S gates for
Z_l, three S_l gates, CX(c,l), S_l, then three S_c gates.

Controlling s_0 is also paid. Chronologically emit CZ(a,l), three S_l
gates, the private-phase CCX(a,c;l), S_l, then CS-dagger(a,c).
CZ is H-CNOT-H; the returned private-phase Toffoli has eight T gates,
ten Cliffords and T-depth at most eight; CS-dagger is the actual
inverse of the three-T, two-CNOT CS word. Thus the common bound is

```math
T_{\mathrm C(s_0)}\le11,\quad G_{\mathrm C(s_0)}\le19,\quad
d_{\mathrm C(s_0)}\le11,\quad h_{\mathrm C(s_0)}\le1.
                                                               \qquad(8)
```

The helper is arbitrary, not initialized. Controlling r_0 costs three
Cliffords. Conjugating either anchor by a bank Clifford conjugates its
controlled version with the same uncontrolled Clifford; it does not
require a controlled-Clifford synthesis oracle.

## 4. Uncontrolled relative phases suffice for the actual trials

Use the actual relative phase and mixer from the
[signal source](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#literal-controlled-t-wall):

```math
\begin{aligned}
D_x&=Q_PQ_{B_x}^\dagger P_tQ_{B_x}P_t^\dagger Q_P,\\
D_x|z,t\rangle&=\omega^{|t+Pz+B_xz|-|t+Pz|}|z,t\rangle,\\
E_x&=\mathcal R^\dagger D_x\mathcal R,\qquad D_0=E_0=I.
\end{aligned}                                                   \qquad(9)
```

Here $`\omega=e^{i\pi/4}`$, additions inside bit strings are XOR,
$`B_x=\bar N_x^T-I`$, P deletes binary coordinate zero, and P_t
is the n-fold T wall. Put $`r=(m-1)p^2`$ and
$`s\le(m-1)p^3`$. The literal bank shear has common bounds
$`T_B\le8r`$, $`d_B\le8`$, $`G_B\le2n+16r+2s`$,
and $`h_B\le r`$.
The two opposite T walls and two Q_P walls therefore give

```math
T_D\le16r+2n,\quad d_D\le18,\quad
G_D\le2G_B+2(n-1),\quad h_D\le r.
                                                               \qquad(10)
```

These are common bounds for the actual reverse as well. The mixer has
$`G_{\mathcal R}\le4n^2+8n`$ and zero T-count. Crucially,
$`\mathrm C(V^\dagger R_0V)=V^\dagger\mathrm C(R_0)V`$.
Only uncontrolled D calls occur below, not controlled T-wall calls.
All helpers can be reused serially. No precision approximation occurs
in these exact Clifford+T words.

<a id="minimal-bank-motion-trial"></a>

## 5. Minimal changed bank motion: a complete paid failure

Let j=n-p and
$`L=\mathrm{CX}(t_j,z_0)\mathrm{CX}(t_j,t_0)`$. Its two CNOTs
commute and it preserves P_0. Set

```math
\begin{aligned}
R'_0&=L^\dagger R_0L=Z_{z_0}Z_{t_j},\\
S'_0&=L^\dagger S_0L
 =\frac{I+P_0}2 Z_{z_0}Z_{t_j}
  +\frac{I-P_0}2 X_{z_0}X_{t_0},\\
R_x&=E_x^\dagger R'_0E_x,\qquad
S_x=D_x^\dagger S'_0D_x.
\end{aligned}                                                   \qquad(11)
```

The anchors use the encoder $`C_0L`$, costing three CNOTs. The pair
is Hermitian and involutory for every coherent address. At zero address
its square is P_0, so (2) is correct on the entire zero-address bank
space. Unlike words made only from E_x and $`P_0D_x^\dagger`$ with
signal interlayers, (11) changes the algebra seen by the last-block
bank bit: R'_0 anticommutes with $`X_{t_j}`$.

Controlled R_x uses two D calls, four mixer calls and nine anchor
Cliffords. Controlled S_x uses two D calls and the controlled anchor
with eleven T gates and 25 Cliffords. The physical source, including
its fixed basis envelope, consequently satisfies

```math
\begin{aligned}
T_U&\le6T_D+11=96r+12n+11,&d_U&\le119,\\
G_U&\le6G_D+8G_{\mathcal R}+4c_n+65,&h_U&\le r.
\end{aligned}                                                   \qquad(12)
```

The direct zero-clean word $`W_x=(R_xS_x)^2`$ has

```math
T_W\le8T_D,\quad d_W\le144,\quad h_W\le r,\quad
G_W\le8G_D+8G_{\mathcal R}+40.                              \qquad(13)
```

Here R'_0 costs four S gates and S'_0 costs sixteen Cliffords. Add
the complete-indicator Clifford envelope in Section 2 to (13).
Thus all intermediate resources have the polynomial overhead allowed
by fusion. The failure is the sector identity, not a hidden depth cost.

The exact p=1 fixture, with carrier ordering z_0,z_1,t_0,t_1, gives
at x=1 and $`|w\rangle=|0000\rangle`$:

```math
\langle w|K_+^\dagger K_+|w\rangle=\frac14,\qquad
\langle w|(I+Q_1)/4|w\rangle=\frac12,\qquad
\|(W_1-Q_1)|w\rangle\|^2=2.                               \qquad(14)
```

In this fixture $`Q_1=Z_{z_0}Z_{z_1}Z_{t_0}Z_{t_1}`$. Thus the
source Gram gap is at least one quarter at this size. This is a finite
failure of (11), not an all-p exclusion: a separately paid small-size
patch would still leave the general identity (6) to prove.

## 6. A broad-anchor, bank-swapped trial also fails its first test

The second trial changes more than the one tail coordinate. Define the
chronological encoder C by C_0, then all CX(z_0,j), then all CX(j,z_0),
where j ranges over the bank wires other than z_0,t_0. It costs 4n-3
CNOTs and obeys

```math
C^\dagger Z_lC=\prod_{j\ne t_0}Z_j,\qquad
C^\dagger X_lC=\prod_jX_j,\qquad C^\dagger Z_cC=P_0.
                                                               \qquad(15)
```

Let $`\bar R_0=C^\dagger r_0C`$ and
$`\bar S_0=C^\dagger s_0C`$. With H the complete 2n-qubit
Hadamard wall and B the exchange of the two n-bit banks, use

```math
V_x=HD_xH,\qquad A_x=BD_xB,\qquad
R_x=V_x^\dagger\bar R_0V_x,\qquad
S_x=A_x^\dagger\bar S_0A_x.                                \qquad(16)
```

B costs 3n CNOTs. Again the pair consists of exact Hermitian
involutions and the entire zero-address sector is correct. Its source
has the same T-count, T-depth and helper bounds as (12), but

```math
G_U\le6G_D+52n+4c_n+29.                                    \qquad(17)
```

For the direct word, the corresponding bound is
$`G_W\le8G_D+64n+4`$, with the T, depth and helper bounds of (13).
This uses the direct 4n-2-Clifford Pauli for $`\bar R_0`$ and the
8n+4-Clifford encoded $`\bar S_0`$.

At p=1,x=1 on the same w, exact arithmetic instead gives

```math
\langle w|K_+^\dagger K_+|w\rangle=\frac38,\qquad
\|(W_1-Q_1)|w\rangle\|^2=1.                               \qquad(18)
```

The Gram gap is at least one eighth at this size. This distinct trial
therefore also needs a changed identity, not just a paid extraction.

<a id="nonrectangular-four-path-acceptance"></a>
<a id="four-path-accepted-fibers"></a>

## 7. Nonrectangular four-path acceptance

The balanced two-unitary source in Section 2 does not exhaust four-path
acceptance. The following calculation tests an actual changed polynomial,
then supplies a separate paid diagnostic with nonpairwise cancellation.
Neither supplies lookup. The [coherent-feedback source assessment](COHERENT_FEEDBACK_SOURCES.md)
records the complementary attempt to mix the key after a nonlinear stage.

All operators act on the entire arbitrary data and private space, with
identity on any external reference. Let A and B be Hermitian involutions
and put $`\omega=e^{i\pi/4}`$. Two clean signals, uniform preparation and
uniform acceptance give the nonrectangular candidate

```math
K=\frac{\omega(I+A)+\bar\omega(B+BAB)}4.
                                                               \qquad(19)
```

The four branches are I, A, B and BAB; their literal SELECT is, chronologically,
$`\mathrm C_a(B),\mathrm C_b(A),\mathrm C_{ab}(B)`$. A subsequent
$`S_a^\dagger`$ supplies (19) up to the common scalar $`\bar\omega`$,
which has no effect on its Gram. This describes controls to be paid,
not free controls of an arbitrary source. At B=I the Gram is correctly
normalized as $`(I+A)/4`$. Writing S=I+A and T=B+BAB gives
$`S^2=2I+2A`$ and $`T^2=2I+AB+BA`$, hence the exact full Gram is

```math
G=K^\dagger K
 =\frac{4I+2A+\{A,B\}+i([B,A]+[BAB,A])}{16}.
                                                               \qquad(20)
```

This candidate has the sharp exact restriction

```math
G(G-I/2)=0\quad\Longleftrightarrow\quad B=I,
\qquad B=I\ \Longrightarrow\ G=(I+A)/4.                \qquad(21)
```

To prove it, the flat-spectrum condition is equivalent to $`Q_*=4G-I`$
being a Hermitian involution. A simultaneous one-dimensional eigenspace
with A=a and B=b gives $`Q_*=a(1+b)/2`$, so b must be plus one.
For a noncommuting two-dimensional block, choose
$`AB=\mathrm{diag}(e^{i\theta},e^{-i\theta})`$ and A=X.
With $`c=\cos\theta,s=\sin\theta`$ and $`|c|\lt1`$, (20) becomes

```math
Q_*=\frac c2I+\frac12X+\left(\frac s2+cs\right)Z.
                                                               \qquad(22)
```

Its trace is c. Two eigenvalues from plus or minus one would therefore
have to be opposite, forcing c=0; at c=0 they are instead plus and minus
$`1/\sqrt2`$. Such a block is impossible. These blocks exhaust the pair:
$`A(AB)A=(AB)^\dagger`$ pairs the nonreal eigenphases, while real
eigenphases give commuting eigenspaces. This proves (21), without a
prescribed rank for G. Exact lookup in (19) consequently requires A
already to be its complete sign.

There is also a conservative uniform approximate version, for any target
Hermitian involution Q, including targets mixing these two-dimensional blocks:

```math
\left\|G-\frac{I+Q}4\right\|\le\delta\lt\frac1{32}
\quad\Longrightarrow\quad
B=I,\qquad\|A-Q\|\le4\delta.                         \qquad(22a)
```

Indeed, an error below $`1/32`$ puts every eigenvalue of $`Q_*`$ within
$`1/8`$ of plus or minus one, by the spectral perturbation bound for
Hermitian operators. A commuting block with b equal to minus one has
$`Q_*=0`$ and is excluded. In a noncommuting block the eigenvalues are
$`c/2\pm R`$, where
$`R^2=1/4+(1-c^2)(1/2+c)^2`$. Both cannot be near the same sign,
since their sum c lies between minus one and one. Opposite signs imply,
with $`t=|c|`$, that $`t\lt1/4`$ and $`R\gt7/8+t/2`$.
However $`R^2\le1/2+t+t^2`$, whereas throughout $`0\le t\le1/4`$,

```math
\left(\frac78+\frac t2\right)^2
 -\left(\frac12+t+t^2\right)
 =\frac3{16}+\left(\frac14-t\right)
                  \left(\frac5{16}+\frac{3t}4\right)
 \ge\frac3{16}.
```

This contradiction proves (22a); the last inequality there follows from
$`G=(I+A)/4`$. The constant is not asserted optimal. This bounds the displayed
acceptance polynomial, not arbitrary multipath sources or native T-depth.

### A different four-path polynomial retains genuine cancellation

Now assume $`AB=-BA`$ and choose the physical branch phases

```math
(U_{00},U_{01},U_{10},U_{11})
 =(I,\omega A,iB,i\omega BAB),\qquad
K=\frac{I+\sqrt2 A+iB}{4}.
                                                               \qquad(23)
```

The same SELECT followed by $`S_aT_b`$ realizes these phases. Direct
multiplication gives

```math
\begin{aligned}
Q&=\frac{A+iAB}{\sqrt2},&Q'&=\frac{A-iAB}{\sqrt2},\\
K^\dagger K&=\frac{I+Q}{4},&KK^\dagger&=\frac{I+Q'}4,\\
K&=e^{i\pi B/4}\frac{\Pi_Q}{\sqrt2},
&\Pi_Q&=\frac{I+Q}{2}.
\end{aligned}                                                   \qquad(24)
```

Both signs are Hermitian involutions with equally sized sectors. Thus K
annihilates the entire Q-negative sector and has singular value
$`1/\sqrt2`$ on its entire positive sector. Its output sector is different:
the polar part in (24) moves it to the Q'-positive sector. There is no
thin-subspace acceptance deficit.

Every pair sum $`U_j+U_k`$ is invertible: the relative spectra are among
$`\{\pm\omega\},\{\pm i\},\{\pm\bar\omega\}`$ and the scalar
minus i, none containing minus one. The zero on the wrong sector therefore
cannot be explained by cancellation of any two paths. Nor does any
individual relative branch word equal Q. This disproves a universal
pairwise-collapse inference; it does not give an asymptotic advantage.
Here Q is itself a simple rotation of A. Producing cheap A,B for a specified
lookup Q remains unpaid.

For a completely paid diagnostic take $`A=Z_d,B=X_d`$ on one arbitrary
data wire. The complete source on two clean signals is chronological:

```math
H_aH_b,\quad\mathrm{CX}(a,d),\quad\mathrm{CZ}(b,d),\quad
\mathrm{CCX}(a,b,d),\quad S_aT_b,\quad H_aH_b.
                                                               \qquad(25)
```

Its accepted block has $`Q=(Z_d-Y_d)/\sqrt2`$ and
$`Q'=(Z_d+Y_d)/\sqrt2`$. A literal three-band CCZ uses the parity identity

```math
4abd=a+b+d-(a\mathbin\oplus b)-(a\mathbin\oplus d)
 -(b\mathbin\oplus d)+(a\mathbin\oplus b\mathbin\oplus d)
 \pmod8.                                                     \qquad(26)
```

Apply T on the three single wires first. Three CNOTs, chronologically
$`\mathrm{CX}(b,a),\mathrm{CX}(a,d),\mathrm{CX}(d,b)`$, expose
the independent parities $`a\oplus b,a\oplus d,a\oplus b\oplus d`$;
apply $`T^\dagger,T^\dagger,T`$ and reverse those CNOTs. Finally expose
$`b\oplus d`$, apply $`T^\dagger`$, and uncompute. Hadamards on d turn
this literal CCZ into CCX. Thus (25), including the CZ Hadamard walls,
has 8 T gates, at most four T bands, 27 forward elementary gates, two clean
signal inputs and no private helpers. Its actual reverse uses 29 elementary
gates when $`S^\dagger`$ is expanded into three positive S gates; hence 29
is a common gate bound, with the same T-count and depth bounds.
Acceptance denotes the Gram interface; no measurement or postselection
is inserted into the eventual complete recovery. This three-qubit fixture
is not a lookup source at arbitrary address length.

<a id="7-scope-evidence-and-reopening-condition"></a>

## 8. Scope, evidence and reopening condition

Equations (1)--(13), (15)--(17) are analytic identities and native
resource bounds. The exact guards in
[`test_coupled_order_gram_source.py`](../../tests/test_coupled_order_gram_source.py)
check the literal controlled anchor, its arbitrary-helper return and
actual inverse; the joint source and complete reverse on a small
noncommuting pair; the actual coherent-address p=1 relative phase; and
both distinct p=1 sector failures (14), (18). The finite fixtures do
not prove an asymptotic exclusion. These are internal analytic and
finite checks, not external mathematical review.

Equations (19)--(26) are separate analytic identities. The exact
[`test_four_path_accepted_sources.py`](../../tests/test_four_path_accepted_sources.py)
checks both four-path Grams, the literal native source and actual inverse,
its arbitrary spectator/reference extension, the moved polar sector,
and nonpairwise cancellation. These checks do not supply inexpensive
general A,B or their native controls.

The complete-frame frontier is unchanged. The source normalization and
joint interference are now paid; the first unsupported inference is
(6) on the correct lookup sectors for a different parameterized pair.
Continuing to reuse only the old E/F signal algebra is separately
blocked by its all-query blind-bank symmetry. The present trials break
that premise but do not repair the desired sector action.

Reopen this mechanism with an actual bank/address-moving pair and a
derivation of its squared product, or a concrete unresolved algebraic
identity that can supply the desired Q_x with uniform error. A p=1
patch alone is not such an identity. More layers or more complicated
readout do not by themselves repair either displayed pair. The subsequent
[affine orbit construction](AFFINE_ORBIT_REFLECTION_SOURCE.md) pays a
fixed-point-free, charge-preserving S and its orbit coordinates. Its local
and scalar-zero choices of R fail; a general second reflection is still unpaid.
The paid
approximate scalar zero reflection remains compatible with these
scoped failures; conjunction is not a dimension-growing fixed-accuracy
depth witness.

The distinct four-path question is to supply an actual native tuple
$`U_0,U_1,U_2,U_3`$ whose six cross terms satisfy

```math
\sum_{j\lt k}(U_j^\dagger U_k+U_k^\dagger U_j)=4Q_x.
```

For equal-amplitude acceptance this is exactly the required complete Gram;
an operator error at most $`16\delta`$ in this cross-term identity gives
Gram error at most delta. Section 7 shows why the zero wrong-sector condition
does not justify imposing pairwise cancellation or a fixed accepted polar
part. A candidate must still pay its SELECT controls and actual inverses,
retain all private inputs, and prove the full uniform identity before it can
use the existing complete recovery.
