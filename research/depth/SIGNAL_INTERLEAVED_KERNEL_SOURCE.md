# Paid signal interleaving and its complete Gram obstructions

[Canonical pencil and kernel](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md) ·
[Relative XOR-carry phase](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-xor-carry-source) ·
[Complete four-call recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Two explicit native families change the acceptance polynomial instead
of conjugating a fixed half-projector. The first uses two clean signals;
the second removes the independent attenuation and uses one clean signal.
Both have constant native T-depth and the permitted polynomial source
budgets. Their full accepted blocks are derived below. Neither supplies
the lookup Gram: each has a uniform error floor at every nonzero fixed
address. These are obstructions to the displayed words, not general
depth bounds. The native controlled phase is a useful completed component.
The changed final coupling in Section 7 is also resolved: a grading
identity excludes every final conjugation preserving the distinguished
X_z0, including conjugations that move the address.

## 1. Registers and the actual phase

Use p address bits, m=2^p, n=pm and the standard dual carrier coordinates
z,t of the pencil. Both n-bit carriers are arbitrary. The physical changes
are $`C_E=U_{E^{-\mathsf T}}\otimes U_{E^{-\mathsf T}}`$; they cost at
most 4c_n Cliffords around a source, with $`c_n=n^2+3n`$.
Let $`B_x=\bar N_x^{\mathsf T}-I`$, and let P delete only the first
binary coordinate e_0. The corank-one matrix
$`K_x=(I+B_x)P`$ has image $`\ker(v_x^{\mathsf T})`$.
The p-1 other binary directions in the first field block remain present.

For a binary matrix L, write
$`Q_L|z,t\rangle=|z,t+Lz\rangle`$ and
$`D(L)=Q_L^\dagger P_tQ_LP_t^\dagger`$, where
$`P_t=\bigotimes_jT_{t_j}`$ and $`\omega=e^{i\pi/4}`$.
This note uses the **relative** kernel phase

```math
 \begin{aligned}
 D_x&=D(K_x)D(P)^\dagger=Q_PD(B_x)Q_P,\\
 D_x|z,t\rangle&=\omega^{h_x(z,t)}|z,t\rangle,\\
 h_x(z,t)&=|t+Pz+B_xz|-|t+Pz|\pmod8.
 \end{aligned}                                                    \qquad(1)
```

The additions inside the weights are binary. These are literal phases,
not integer replacements for XORs. In particular D_0=I, and D_x is
independent of z_0. Every query helper is returned on arbitrary input.
The paid standard-coordinate query has

```math
 T_B\le8r,\quad d_B\le8,\quad h_B\le r,\quad
 G_B\le2n+16r+2s,\qquad
 r=(m-1)p^2,\quad s\le(m-1)p^3.                         \qquad(2)
```

Q_P uses n-1 fixed CNOTs. The fixed bank mixer used below is
$`\mathcal R=H^{\otimes2n}U_M`$, $`M(z,t)=(Ez+t,z)`$, with
$`G_{\mathcal R}\le4n^2+8n`$. Define

```math
 E_x=\mathcal R^\dagger D_x\mathcal R,\qquad
 F_x=P_0D_x^\dagger,\qquad P_0=Z_{z_0}Z_{t_0}.           \qquad(3)
```

The accepted lookup Gram in these coordinates is
$`G_*(x)=(I+Q_x)/4`$, where $`Q_x=Z_z(v_x)Z_t(v_x)`$.
No principal carrier, phase register or coefficient bank is initialized.

<a id="literal-controlled-t-wall"></a>

## 2. A literal controlled T wall and controlled relative phase

For a unitary V, the notation

```math
 \mathrm C_a(V)=|0\rangle\langle0|_a\otimes I+
 |1\rangle\langle1|_a\otimes V
```

denotes the full controlled operator.
Its native implementation is paid here, rather than supplied as an oracle.

First implement CT(a,b), whose phase is $`\omega^{ab}`$.
Borrow four arbitrary bits d,e,f,g. The returned private-phase Toffoli
$`T_{ab;d}`$ uses eight T/T-dagger gates, eight CNOTs and two Hadamards,
with g as its phase helper. The chronological word

```math
 T_{ab;d},\quad T_d,\quad T_{ab;d},\quad T_d^\dagger
```

has phase $`\omega^{ab-2abd}`$ and returns d,g. Follow it by the
[52-T CCS word](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-half-bilinear-source)
on logical controls a,b,d, using e,f,g as its three arbitrary helpers.
Its phase $`i^{abd}=\omega^{2abd}`$ cancels the entire dirty-dependent
residual. The resulting CT and its actual reverse therefore have

```math
 T_{\mathrm{CT}}=70,\qquad G_{\mathrm{CT}}=88,\qquad
 d_{\mathrm{CT}}\le70,\qquad h_{\mathrm{CT}}=4.          \qquad(4)
```

There is no discarded global phase. Every physical T target is one of
d,e,f,g; the logical a,b wires occur only as Clifford controls. Give
each target t_j its own four helpers to emit a full controlled P_t
wall in the same 70 T layers, despite the shared signal a. Every helper
returns before the next wall or query call.

Controlled conjugation then gives the exact complete identity

```math
 \mathrm C_a(D_x)
 =Q_PQ_{B_x}^\dagger\,\mathrm C_a(P_t)\,
   Q_{B_x}\,\mathrm C_a(P_t)^\dagger Q_P.              \qquad(5)
```

On a=0, the uncontrolled queries and fixed shears cancel literally;
on a=1, (5) is (1). Thus no controlled quadratic query is needed.
The helpers used by the query and the two phase walls can be reused
after exact return. Safe common forward/actual-reverse bounds are

```math
 \begin{aligned}
 T_C&\le16r+140n,&d_C&\le156,\\
 G_C&\le2G_B+178n,&h_C&\le\max(r,4n).
 \end{aligned}                                                   \qquad(6)
```

These identities hold with coherent signal and address inputs and with
arbitrary carrier, helper and reference states. They also pay
$`\mathrm C_a(E_x)=\mathcal R^\dagger\mathrm C_a(D_x)\mathcal R`$
and $`\mathrm C_a(F_x)=\mathrm C_a(P_0)\mathrm C_a(D_x)^\dagger`$.
The former uses two ordinary mixer calls; the latter adds two literal
CZ gates, costing six elementary Cliffords. Multiplication of E or F
by i is paid by S_a, not an unimplemented scalar on the controlled branch.

<a id="attenuated-signal-polynomial"></a>

## 3. First interleaving: complete polynomial and an all-size failure

At each address choose $`A=E_xX_{z_0}`$ and $`B=F_x`$, and emit on signal a

```math
 W=H_a\,\mathrm C_a(B)\,H_a\,\mathrm C_a(A)\,H_a.
```

A second independent clean signal b receives H_b. With both signals
inserted and accepted at zero, direct block multiplication gives

```math
 K_2=\frac{I+A+B-BA}{4},\qquad
 G_2=\frac1{16}(I+A+B-BA)^\dagger(I+A+B-BA).             \qquad(7)
```

This is the full accepted action on both arbitrary banks, not a trace
promise or a conjugated fixed projector. At x=0, A=X_z0 and B=P_0
anticommute, and (7) gives $`G_2=(I+P_0)/4`$, the correct base sector.
The source and its actual reverse have

```math
 \begin{aligned}
 T_2&\le32r+280n,&d_2&\le312,&h_2&\le\max(r,4n),\\
 G_2^{\rm native}&\le4G_B+356n+2G_{\mathcal R}+4c_n+12.
 \end{aligned}                                                   \qquad(8)
```

The conservative Clifford constant also covers the one tested change
$`X_{z_0}\mapsto X_{z_0}X_{z_p}`$. A controlled X product is just
one or two signal-to-bank CNOTs. This source has two clean inputs;
there is no claim that the one-clean recovery theorem applies to it.

The independent H_b attenuation also means that $`G_2=G_W/2`$.
Matching the desired half-projector would require $`G_W=\Pi_{\rm eq}`$:
the one-signal core would have to accept the equal sector with probability
one and the unequal sector with probability zero. Thus this particular
attenuation imposes a stronger intermediate discrimination task.

More decisively, (7) fails for every choice of the unitary A with the
stated B. If Bv=v for a unit vector v, then

```math
 \langle v|K_2|v\rangle=\frac12,\qquad
 \|K_2v\|^2\ge\frac14.                                 \qquad(9)
```

For any nonzero x, choose j different from zero with (v_x)_j=1, and
put $`v=|z=0,t=e_j\rangle`$. Such a j exists in the second field block.
Equation (1) gives D_xv=v, and P_0v=v, hence Bv=v. But Q_xv=-v, so
G_*v=0. Therefore

```math
 \|G_2-G_*(x)\|\ge1/4\qquad(x\ne0),                    \qquad(10)
```

independently of the mixer, the chosen X coupling or any other change
to A. It is the noncancellable I+B branch, not a normalization constant
or a small-size anomaly, that stops this polynomial with this final B.

<a id="one-clean-signal-polynomial"></a>

## 4. A direct one-clean interleaving and its full Gram

Remove b entirely and change both controlled branches. Put
$`A=iE_x`$, $`B=iF_x`$, retain the three-H word W, and use

```math
 U_1=T(Y_a)\,W\,S_a\,T(Y_a),\qquad
 T(Y)=e^{i\pi/8}R_y(\pi/4).                             \qquad(11)
```

The actual chronological word for T(Y) is S,S,S,H,T,H,S. Its actual
inverse has the same one-T, six-Clifford bound after reducing adjacent
S words. There is no independent attenuation signal in (11).
With only a inserted and accepted at zero, define

```math
 \alpha=1+i(\sqrt2-1),\qquad \beta=1+\sqrt2+i.
```

Multiplying the two-by-two signal blocks, retaining their literal
phases, gives the complete accepted operator

```math
 K_1=\frac{1+i}{8}
       \left(\alpha I+\overline\alpha A+\beta B-
                         \overline\beta BA\right).     \qquad(12)
```

For an expanded full Gram, put $`M_F=(F+F^\dagger)/2`$ and
$`N_F=(F-F^\dagger)/(2i)`$. Then

```math
 \begin{aligned}
 G_1={}&\frac I2+\frac{\sqrt2}{8}
                 (E^\dagger N_FE-N_F)\\
       &-\frac{(1+i)E+(1-i)E^\dagger}{8}\\
       &+\frac{(1-i)M_FE+(1+i)E^\dagger M_F}{8}.
 \end{aligned}                                                   \qquad(13)
```

No commutation of E and F is assumed. One way to check (13) is to
write (12) as $`(1+i)(X+YE)/8`$, with
$`X=\alpha I+i\beta F`$, $`Y=i\overline\alpha I+\overline\beta F`$.
Then $`X^\dagger X=8I-4\sqrt2N_F`$,
$`Y^\dagger Y=8I+4\sqrt2N_F`$ and
$`X^\dagger Y=-4(1+i)I+4(1-i)M_F`$.

At x=0, E=I and F=P_0. On a P_0 eigensector q, the literal core
obeys $`W_qS_a=\omega X_a^{(1-q)/2}`$. Consequently (11) is
$`iR_y(\pi/2)`$ on q=+1 and iX on q=-1, on either signal input.
In particular

```math
 K_1(0)=i\Pi_0/\sqrt2,\qquad G_1(0)=\Pi_0/2.           \qquad(14)
```

The physical source is C_E-dagger U_1 C_E. Its common forward and
actual-reverse ledger is

```math
 \begin{aligned}
 T_1&\le32r+280n+2,&d_1&\le314,&h_1&\le\max(r,4n),\\
 G_1^{\rm native}&\le4G_B+356n+2G_{\mathcal R}+4c_n+30,
 &\varepsilon_1&=0.
 \end{aligned}                                                   \qquad(15)
```

The constant 30 includes both T(Y) Clifford words, three Hadamards,
two CZs and three signal S occurrences charged at three Cliffords each
for a common inverse bound. The two carriers are principal inputs;
embedding them around the m output wires adds (2p-1)m arbitrary wires
besides h_1. These source counts and widths fit the polynomial fusion
envelope. Clean return is a desired source/readout consequence, not an
assumption used to derive (12); private implementation helpers return
exactly at each completed primitive.

<a id="flat-phase-one-clean-obstruction"></a>

## 5. The flat mixer still has the wrong accepted sectors

Use the same wrong-sector vector $`v=|0,e_j\rangle`$ from Section 3.
Now Fv=v and Bv=iv. The coefficient identities
$`\alpha+i\beta=2i\sqrt2`$ and
$`\overline\alpha-i\overline\beta=-2i\sqrt2`$ give

```math
 \langle v|K_1|v\rangle
   =\frac{(i-1)\sqrt2}{4}
        \left(1-i\langle v|E|v\rangle\right).           \qquad(16)
```

Every computational column of the actual mixer
$`\mathcal R=H^{\otimes2n}U_M`$ has uniform modulus. Thus
$`\langle v|E_x|v\rangle=2^{-2n}\mathrm{Tr}(D_x)=\mu_x`$.
Conjugation by Q_P does not change this trace. For fixed z, averaging
the bare XOR-carry phase over t contributes
$`(1/\sqrt2)^{|B_xz|}`$, so

```math
 \mu_x=2^{-n}\sum_z(1/\sqrt2)^{|B_xz|}
       =\left(\frac{2+\sqrt2}{4}\right)^{n-p}\gt0
       \qquad(x\ne0).                                  \qquad(17)
```

The last equality uses that B_x maps onto all first n-p binary
coordinates when x is nonzero; its last field block is zero. Since the
required Gram annihilates v, (16) gives the all-size complete bound

```math
 \|G_1-G_*(x)\|
   \ge\|K_1v\|^2
   \ge\frac{1+\mu_x^2}{4}\ge\frac14
       \qquad(x\ne0).                                  \qquad(18)
```

This includes arbitrary spectator identities. No physical trace
estimation is used: the trace only evaluates a diagonal amplitude of
the complete native word. The proof also covers any replacement mixer
whose relevant column is uniformly flat over the phase basis. A full
address Hadamard together with such a bank mixer merely averages the
real quantities mu_x, and does not remove this witness.

There is an exact missing sector condition behind (18). On **any**
B=i eigenvector v, (12) has zero acceptance if and only if

```math
 E v=-iv.                                                \qquad(19)
```

Indeed zero acceptance forces its diagonal amplitude to vanish, hence
$`\langle v|A|v\rangle=1`$. Since A is unitary, Av=v, or Ev=-iv.
Conversely those eigenvalue identities cancel all four terms in (12).
The flat mixer supplies a real diagonal expectation instead. This
stops the displayed one-clean phase coupling; it does not exclude
all signal interleavings or all one-clean normalized sources.

## 6. Exact fixtures for the first two words

The [exact fixtures](../../tests/test_signal_interleaved_kernel_source.py)
check CT on all logical and dirty columns and its actual reverse, the
controlled XOR phase including coherent address and every dirty input,
the signal-block polynomial and base normalization, and the following
finite unequal-sector receipts at x=1,z=e_p,t=0:

| Source | p=1 accepted norm squared | p=2 accepted norm squared |
| --- | --- | --- |
| Two-clean word, X_z0 coupling | $`1/4+\sqrt2/8`$ | $`7991/32768+(1391/16384)\sqrt2`$ |
| Single changed coupling X_z0 X_zp | $`3/8+\sqrt2/16`$ | $`9305/32768+\sqrt2/16`$ |
| Direct one-clean word | $`1/2`$ | $`1197/2048-(155/8192)\sqrt2`$ |

The p=2 calculation uses integer coefficients in Z[omega], with Python
integers before products. These fixtures guard the complete analytic
identities; the all-size conclusions are (10) and (18), not an
extrapolation from the table. They and the independent internal audit
are not external mathematical review. No resource frontier moves.

<a id="partial-h-final-coupling"></a>

## 7. The changed final coupling and its complete intersection

The next specified coupling changes the eigenspaces of the **final**
controlled phase. Let S apply H to z_j and t_j for every j at least
one, leaving z_0,t_0 unchanged, and replace B in (11) by

```math
 B'_x=iS^\dagger F_xS.                                   \qquad(20)
```

S commutes with P_0, so the exact x=0 action (14) is retained.
The literal controlled operation is
$`S^\dagger\mathrm C_a(B)S`$, costing 4(n-1) extra Cliffords
and no extra T gates, width, clean inputs or precision. Equation (12)
still gives the actual accepted block with B replaced by B'. Put
$`F'_x=S^\dagger F_xS`$, $`M'=(F'+F'^\dagger)/2`$ and
$`N'=(F'-F'^\dagger)/(2i)`$. The complete accepted Gram is

```math
 \begin{aligned}
 G'={}&\frac I2+\frac{\sqrt2}{8}(E^\dagger N'E-N')\\
      &-\frac{(1+i)E+(1-i)E^\dagger}{8}\\
      &+\frac{(1-i)M'E+(1+i)E^\dagger M'}{8}.
 \end{aligned}                                                   \qquad(21)
```

There is no discarded off-diagonal action in this formula. The native
controlled conjugation and every returned helper from Section 2 remain
valid on arbitrary coherent carrier, address and reference inputs.

The intersection suggested by (19) can be calculated exactly. At a
fixed nonzero address set $`w=v_x+e_0`$, so w_0=0 and w is nonzero.
In S coordinates the target charge becomes

```math
 S Q_xS=P_0X_z(w)X_t(w).                                  \qquad(22)
```

For $`y=(z,t)`$, put $`y'=(z+w,t+w)`$ and
$`b=z_0\mathbin\oplus t_0`$. These are disjoint two-element orbits,
and b is unchanged on each orbit. Since F is diagonal, the whole
intersection
$`S^\dagger\ker(F_x-I)\cap\mathrm{ran}\Pi_{\rm neq}(x)`$
is spanned by exactly the vectors

```math
 S^\dagger\frac{|y\rangle-(-1)^b|y'\rangle}{\sqrt2}
 \quad\text{for orbits with}\quad
 h_x(y)\equiv h_x(y')\equiv4b\pmod8.                    \qquad(23)
```

To prove completeness, membership in ker(F-I) permits only basis
labels with F eigenvalue one; the wrong-charge equation from (22)
then pairs both coefficients on each orbit with the displayed sign.
The variable a=t+Pz is unchanged when y is replaced by y', whereas
B_xz changes to B_xz+B_xw. Therefore

```math
 h_x(y')-h_x(y)\equiv|B_xw|\pmod2.                       \qquad(24)
```

An odd weight makes (23) empty. In particular, at x=1 the matrix B_1
shifts field blocks left by one, and v_1 consists of m copies of the
binary unit representing the field element one. The vector w deletes
only its first copy. Thus B_1w has m-1 one-bits, an odd number, and
the intersection is zero for every p. The earlier B=i eigenvector
witness does not survive this coupling. This does not establish the
Gram contract: the independent grading argument below still excludes it.

There is also a direct symmetry check specific to this S. The operator
E commutes with every bank X because it is a reversible linear
conjugate of the X-diagonal operator H_all D H_all. The operator F'
commutes with X_zj and X_tj for j at least one because S maps them to
diagonal Z operators. Hence G' commutes with them. At a nonzero x,
some such X_zj flips Q_x, already forcing Gram error at least 1/4.
The next section gives a stronger obstruction with a wider scope.

<a id="distinguished-x-grading-obstruction"></a>

## 8. The distinguished-coordinate grading excludes a larger family

Keep E and the signal polynomial (12), and replace F by
$`\widetilde F=V^\dagger FV`$ for any unitary V commuting with
$`X=X_{z_0}`$. This includes S, without requiring V to be Clifford,
local or diagonal. First work at a fixed nonzero address. The actual
first branch obeys [E,X]=0 as above. Since D is independent of z_0,
F=P_0D-dagger anticommutes with X, and so does its conjugate:

```math
 XEX=E,\qquad X\widetilde FX=-\widetilde F.              \qquad(25)
```

Write the Gram as (21) using the real and imaginary parts of
tilde F. Averaging it under X kills every term involving those parts.
The averaging map $`\mathcal A_X(L)=(L+XLX)/2`$ is a contraction
in operator norm. Since (v_x)_0=1, the target charge anticommutes
with X, so

```math
 \begin{aligned}
 \mathcal A_X(\widetilde G)
   &=\frac I2-\frac{(1+i)E+(1-i)E^\dagger}{8},\\
 \mathcal A_X(G_*)&=\frac I4,\\
 \|\widetilde G-G_*\|
   &\ge\left\|\frac I4-
                 \frac{(1+i)E+(1-i)E^\dagger}{8}\right\|.
 \end{aligned}                                                   \qquad(26)
```

This is an exact spectral obstruction, independent of the chosen
final conjugation. Since E is unitary and conjugate to D, its right
side is

```math
 \Delta_x=\max_{\lambda\in\mathrm{spec}(D_x)}
       \left|\frac14-\frac{\mathrm{Re}((1+i)\lambda)}4\right|.
                                                                    \qquad(27)
```

For nonzero x the image of B_x is all first n-p binary coordinates.
For every integer k between zero and n-p choose z with |B_xz|=k and
**set t=Pz**. Then a=t+Pz=0 in (1), so h_x=k and omega^k is an
actual eigenvalue of D_x. For p=1, its phases are exactly
$`\{1,\omega,\omega^{-1}\}`$, and lambda=omega gives Delta=1/4.
For p at least two, n-p is at least three; lambda=omega^3 gives
$`\mathrm{Re}((1+i)\omega^3)=-\sqrt2`$. This attains the maximum
possible value in (27) over the entire unit circle. Consequently

```math
 \|\widetilde G-G_*(x)\|\ge\Delta_x
 =\begin{cases}
     1/4,&p=1,\\
     (1+\sqrt2)/4,&p\ge2,
   \end{cases}
 \qquad x\ne0.                                           \qquad(28)
```

The same proof applies to the full coherent address space when V
moves addresses but still commutes with X_z0: take E,F and the target
as their address-controlled direct sums in (25)--(27). The resulting
bound is global; a per-address bound requires address-preserving V.
Arbitrary spectator identities do not change either statement. No
native cost is asserted for a general V; the specified S has the
paid ledger above.

Even before using the actual XOR spectrum, exact matching in (26)
would force

```math
 \mathrm{spec}(E)\subseteq\{1,-i\}.                      \qquad(29)
```

Indeed the real part in (27) equals one on the unit circle only at
those two points. Thus a different final coupling preserving this
grading cannot compensate for the current first-branch spectrum.
These are obstructions to this signal polynomial and its indicated
grading, not unrestricted bounds on one-clean Gram sources.

<a id="local-distinguished-rotation-obstruction"></a>

## 9. A paid distinguished-coordinate rotation and its trace obstruction

A local rotation of z_0 breaks the grading in Section 8. This is an
actual alternative native family, rather than an unpriced change of
basis. Let $`V_x=I`$ at x=0 and
$`V_x=R_y(\theta)_{z_0}`$ otherwise, and replace the final F by
$`F_x^{(\theta)}=V_x^\dagger F_xV_x`$. For theta equal to pi/4,

```math
 F_x^{(\pi/4)}
   =\frac{Z_{z_0}-X_{z_0}}{\sqrt2}\,Z_{t_0}D_x^\dagger
       \quad(x\ne0).                                      \qquad(30)
```

The complete accepted block is still (12), with B equal to
$`iF_x^{(\theta)}`$, and its complete Gram is (13) with this changed
F. In particular the base normalization (14) is retained.

To emit the address-dependent V, use the
[complete scalar zero reflection](EXPANDER_ZERO_REFLECTION.md#7-arbitrary-lengths-and-the-complete-zero-clean-ledger)
on the p+1 specified bits (x,z_0). Let its actual word be
$`\widetilde R_0`$, with full-operator error at most zeta, including
identity on every arbitrary helper. Put
$`\widetilde N=H_{z_0}\widetilde R_0H_{z_0}`$.
Its ideal comparison N is -X_z0 at x=0 and I otherwise.
Let $`\widetilde A`$ be a finite native approximation to
$`A=R_y(\theta/2)_{z_0}`$ within sigma, retaining the literal
single-qubit phase. Emit

```math
 \widetilde V=\widetilde A\widetilde N
               \widetilde A\widetilde N^\dagger.           \qquad(31)
```

At x=0, its ideal word is A X A X=I; otherwise it is A squared,
which is the required rotation. The two minus signs at x=0 cancel.
The actual inverse in (31) returns the same arbitrary helper pool
in the comparison, and telescoping gives
$`\|\widetilde V-V\otimes I\|\le2\zeta+2\sigma`$.
For the fixed theta=pi/4, its half-angle is synthesized at the fixed
requested precision; it is not treated as a native fine-angle gate.

Conjugate the paid controlled F by the actual V and its actual
inverse. If $`T_0,G_0,d_0,b_0`$ are common costs of the scalar word
and its actual inverse and $`t_A,g_A,d_A`$ are common costs for A,
the resulting source has the following additions to (15):

```math
 \begin{aligned}
 \Delta T&\le4T_0+4t_A,&
 \Delta G&\le4G_0+4g_A+8,&
 \Delta d&\le4d_0+4d_A,\\
 h&\le\max(h_1,b_0),&
 \varepsilon_U&\le4\zeta+4\sigma,&
 \varepsilon_G&\le8\zeta+8\sigma.
 \end{aligned}                                                   \qquad(32)
```

Here h_1 is the private width in (15); neither the carrier z_0 nor
the coherent address is an initialized helper. Sequential calls may
reuse the returned pools. The emitted source's actual reverse has the
same bounds. At fixed requested precision these additions have
constant T-depth and O(p) T-count and dirty width, so they fit the
permitted polynomial source envelope. The estimates are complete
operator and accepted-Gram estimates on arbitrary coherent address,
carrier, dirty-work and reference inputs.

Nevertheless no address-controlled single-qubit unitary V_x on z_0
can repair this acceptance polynomial with the present first branch E.
This stronger statement does not assume that V_x commutes with X_z0.
Fix a nonzero address and let the bank dimension be $`d=2^{2n}`$.
Write $`F=Z_{z_0}L`$, where $`L=Z_{t_0}D_x^\dagger`$ is independent
of z_0. For every local V,
$`V^\dagger Z_{z_0}V=aX_{z_0}+bY_{z_0}+cZ_{z_0}`$ with real
a,b,c. The operator E is X-diagonal: all its computational diagonal
entries and each fixed-displacement matrix entry are constant.
The z_0 signs therefore cancel the Y and Z contributions to
$`\mathrm{Tr}(V^\dagger FV E)`$ and
$`\mathrm{Tr}((V^\dagger FV)^\dagger E)`$.

The X contribution vanishes as well. The canonical binary change of
basis has $`Ee_0=e_0`$, because its first column is beta_0 v_0=e_0.
For the linear mixer M in Section 1 this gives
$`M(e_0,0)=(e_0,e_0)`$. Consequently the constant displacement
coefficient of the operator E_x is

```math
 \kappa_x=\frac1d\mathrm{Tr}(E_xX_{z_0})
  =\frac1d\mathrm{Tr}(D_xZ_{z_0}Z_{t_0})=0,             \qquad(33)
```

where the last equality uses the exact independence of D_x from z_0.
Thus, if $`M_V=(V^\dagger FV+(V^\dagger FV)^\dagger)/2`$,
then $`\mathrm{Tr}(M_VE)=0`$. Taking the trace of (13), the
N terms cancel by cyclicity. Equation (17) gives the remaining
real trace of E, so

```math
 \begin{aligned}
 \frac1d\mathrm{Tr}(G_V)&=\frac12-\frac{\mu_x}{4},\\
 \frac1d\mathrm{Tr}(G_*(x))&=\frac14,\\
 \|G_V-G_*(x)\|
   &\ge\frac{1-\mu_x}{4}
    \ge\frac{2-\sqrt2}{16}\gt0,\qquad x\ne0.
 \end{aligned}                                                   \qquad(34)
```

The first lower bound tends to 1/4 with p. It applies to every
address-controlled single-qubit V_x, including arbitrary coherent
address inputs through the direct-sum operator. It does not cover
couplings of z_0 to another carrier or address-changing V. For the
emitted approximate family (31), subtract the Gram error in (32)
from the bound. At p=1,x=1,theta=pi/4 the exact unequal-sector input
z=e_1,t=0 has accepted norm squared 7/16. This finite stronger
receipt is not extrapolated to larger p.

For a different final phase F, even matching the target trace with
this first branch and signal polynomial requires the mixed moment

```math
 \mathrm{Re}\!\left[(1-i)\frac{\mathrm{Tr}(M_FE_x)}d\right]
      =\mu_x-1.                                             \qquad(35)
```

The local family forces its left side to zero, whereas the required
value tends to -1. A final coupling outside this family must generate
that mixed moment as well as satisfy the full operator equation (13).

<a id="branch-moment-obstruction"></a>

## 10. Branch moments exclude every conjugation repair for p at least three

The mixed-moment requirement (35) gives a stronger obstruction that
does not require the actual mixer, any Pauli grading or a local final
coupling. Keep the signal polynomial (12), but allow independent
arbitrary unitary conjugates of both branch operators:

```math
 E'=R^\dagger D_xR,\qquad
 F'=V^\dagger(P_0D_x^\dagger)V.
```

At first fix a nonzero address. Any arbitrary spectator space is
included by tensoring the displayed original operators with identity;
R and V may mix those spectators with the carriers. Both act only on
the source's nonsignal registers: mixing the initialized signal would
change the polynomial premise. External references always tensor identity.
Write t for the
normalized trace on this complete space and put

```math
 \ell=n-p,\qquad a={2+\sqrt2\over4},\qquad
 \mu=a^\ell,\qquad \rho=2^{-\ell},\qquad
 H_{E'}={(1-i)E'+(1+i)E'^\dagger\over2}.
```

Here H is Hermitian. The two moments of the original diagonal phase
are exact:

```math
 \mathrm t(D_x)=\mu,\qquad \mathrm t(D_x^2)=\rho.
                                                               \qquad(36)
```

Indeed B_xz is uniform on the first ell binary coordinates, and
u=t+Pz is independent and uniform when taking the complete trace.
Each coordinate contributes h=0 with probability one half and h=+1
or -1 with probability one quarter each. Its first phase moment is
a and its second is one half. Independence proves (36). The carrier
t in this explanation is a bit string, distinct from normalized
trace t in the formulas.

Unitary conjugation preserves both moments. Also P_0 commutes with D,
so the square of the original F is D-dagger squared. For
$`M'=(F'+F'^\dagger)/2`$, these facts give

```math
 \mathrm t(M'^2)={1+\rho\over2},\qquad
 H_{E'}^2=I+\mathrm{Im}(E'^2),\qquad
 \mathrm t(H_{E'}^2)=1.                                  \qquad(37)
```

Take the trace of the complete Gram (13). The N terms cancel by
cyclicity, and the remaining cross term is a Hermitian inner product:

```math
 \mathrm t(G')={1\over2}-{\mu\over4}
                 +{\mathrm t(M'H_{E'})\over4},\qquad
 \mathrm t(G_*)={1\over4}.                               \qquad(38)
```

Cauchy--Schwarz for the normalized Hilbert--Schmidt inner product
therefore proves, for every such independent pair of conjugations,

```math
 \|G'-G_*(x)\|\ge
 {1-\mu-\sqrt{(1+\rho)/2}\over4}.                        \qquad(39)
```

The right side is positive for every p at least three. For an exact
rational certificate, ell is then at least 21, a is less than 7/8,
and $`(7/8)^{21}\lt1/16`$, while
$`\sqrt{(1+2^{-21})/2}\lt3/4`$. Thus (39) is strictly greater
than 3/64 throughout that range. Its exact expression tends to
$`(1-1/\sqrt2)/4`$ as p grows. A nonpositive right side at p=1 or
p=2 is not an exclusion of those cases by this argument; the earlier
grading and local-coupling results have their own smaller-p scopes.

### Coherent address motion in both conjugations

The same proof applies when R and V act on the full coherent address
space and need not preserve any address fiber. Replace D_x and
P_0D_x-dagger by their address-controlled direct sums, and use the
normalized trace on all address, carrier and spectator inputs.
At x=0, D_0=I, whereas at each other address (36) holds. Hence put

```math
 \overline\mu={1\over m}+\left(1-{1\over m}\right)\mu,
 \qquad
 \overline\rho={1\over m}+\left(1-{1\over m}\right)\rho.
```

Both branch conjugations preserve these global moments. The complete
target still has normalized trace one quarter, giving the global bound

```math
 \|G'-G_*\|\ge
 {1-\overline\mu-\sqrt{(1+\overline\rho)/2}\over4}
       \gt {7\over512}\qquad(p\ge3).                    \qquad(40)
```

For the last strict inequality, m is at least eight,
$`\overline\mu\lt23/128`$ and
$`\sqrt{(1+\overline\rho)/2}\lt49/64`$, using the same
rational estimates above. This is a global bound, not a claim about
fixed-address blocks after address motion. Arbitrary dirty dimensions
and reference identities do not change the normalized moments.

Thus changing either mixer, including genuinely entangling couplings
or coherent address motion, cannot repair this signal polynomial
while retaining these two branch spectra. The theorem needs only the
moments (36)--(37), so it also applies to any other branches satisfying
those same moment identities. If an emitted source is within full
operator error epsilon of a displayed ideal source, its Gram error
floor decreases by at most 2 epsilon. No native schedule for arbitrary
R or V is assumed: the obstruction already permits them at no cost.
This is not an obstruction to general higher-band sources or other
acceptance polynomials.

## 11. Verification scope and the next unpaid operation

The focused exact fixtures multiply the full p=1 signal operator
with the specified S, verify its accepted polynomial and grading
average, and check that the intersection in (23) is empty. They also
check the attained phase and spectral-gap arithmetic. At p=1,x=1,
the unequal basis input z=e_1,t=0 has acceptance probability 1/2;
that sharper finite receipt is not claimed as an all-p bound.
The all-size results are the orbit proof (23)--(24), the grading
proof (25)--(29), the local-rotation trace proof (33)--(34), and
the branch-moment obstruction (36)--(40) in its p-at-least-three scope.
A separate exact fixture verifies the full pi/4 signal word, the trace
cancellation, and the 7/16 finite unequal-sector receipt. Exact phase
counts also verify both moments and the rational positivity certificate.
Internal independent audit and exact fixtures are
distinct from external mathematical review.

This final-coupling route is paused with its native components
retained. A next source must change the actual branch moments or the
acceptance polynomial. Changing only the first mixer, adding a joint
final coupling, or moving addresses cannot overcome (40). The earlier
suggestion to escape solely by mixing the distinguished coordinate is
therefore superseded. The first unpaid operation remains a normalized
complete accepted Gram on the correct lookup sectors; none of these
necessary changes supplies that operation by itself.
