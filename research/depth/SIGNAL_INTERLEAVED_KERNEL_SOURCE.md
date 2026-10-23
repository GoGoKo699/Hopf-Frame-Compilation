# Paid signal interleaving and its complete Gram obstructions

[Canonical pencil and kernel](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md) ·
[Relative XOR-carry phase](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-xor-carry-source) ·
[Complete four-call recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Explicit native families change the acceptance polynomial instead
of conjugating a fixed half-projector. The first uses two clean signals;
the second removes the independent attenuation and uses one clean signal.
Both have constant native T-depth and the permitted polynomial source
budgets. Their full accepted blocks are derived below. Neither supplies
the lookup Gram: each has a uniform error floor at every nonzero fixed
address. These are obstructions to the displayed words, not general
depth bounds. The native controlled phase is a useful completed component.
The changed final coupling in Section 7 is also resolved: a grading
identity excludes every final conjugation preserving the distinguished
X_z0, including conjugations that move the address. Section 11 changes
the polynomial again with three controlled branches, including a native
signal seed; its complete action gives a separate scoped obstruction.
Section 12 makes the outer branches unequal. A complete rejection
condition excludes every left unitary and scalar signal seed while
retaining the specified right and middle branches; a moving-address
right-branch trial instead fails its base normalization.
Section 13 tests an actual two-signal coupling and identifies a common
bank Pauli that excludes every word built from the retained branches,
regardless of signal count or call count. A changed bank Hadamard leaves
that grammar but its displayed three-query word fails a finite sector test.

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

<a id="three-branch-seeded-source"></a>

## 11. Three controlled branches and a paid signal seed

Conjugating a controlled E by two controlled F branches is an actual
change of acceptance algebra. Keep the literal rotation
$`R_a=T(Y_a)`$ from (11), and put

```math
 \begin{aligned}
 P&=(I+P_0)/2,&
 B_F&=H_a\mathrm C_a(F)H_a,\\
 U_0&=R_a^3\mathrm C_a(P_0)R_a^\dagger,&
 W_\theta&=B_F R_y(\theta)_a^\dagger
           \mathrm C_a(E)R_y(\theta)_aB_F^\dagger,\\
 U_\theta&=U_0W_\theta.
 \end{aligned}                                                   \qquad(41)
```

The two outer branches are actual inverses. Since E_0=I, W_theta is
identity at address zero for every theta. The base source has accepted
row $`\omega[P/\sqrt2,\ I-P-P/\sqrt2]`$; hence its accepted
block is $`(1+i)P/2`$ and its Gram is P/2. Thus (41) retains
the exact required base normalization without another clean signal.

Here is the complete accepted polynomial. Write
$`s=\sin(\theta/2)`$, $`c=\cos(\theta/2)`$,
$`u=(s+c)/2`$, $`v=(s-c)/2`$, and define bank operators

```math
 L=uI+vF,\qquad H=uI-vF,\qquad
 C=\sqrt2\,vPF+(I-P)H,\qquad \Delta=E-I.
```

The letters H and C in this display are bank operators, not a
Hadamard or a controlled-gate abbreviation. Multiplying the signal
blocks gives

```math
 \begin{aligned}
 (W_\theta)_{00}&=I+L\Delta L^\dagger,&
 (W_\theta)_{10}&=H\Delta L^\dagger,\\
 K_\theta&=\omega\left(P/\sqrt2+C\Delta L^\dagger\right),\\
 G_\theta&=P/2+
     {PC\Delta L^\dagger+L\Delta^\dagger C^\dagger P\over\sqrt2}
       +L\Delta^\dagger C^\dagger C\Delta L^\dagger.
 \end{aligned}                                                   \qquad(42)
```

This is a full operator identity; E and F need not commute. Every
address, carrier and private-work input can be coherent or entangled
with a reference. Only the one signal is initialized and accepted.

### Native implementation and the protected unseeded sector

At theta equal to zero, (41) uses three paid controlled phase calls.
At theta equal to pi/4, use R_a and its actual inverse for the seed:
their scalar phases cancel exactly. Together with the physical C_E
envelope, safe common source and actual-inverse bounds are

```math
 \begin{gathered}
 \begin{array}{c|c|c}
  &\theta=0&\theta=\pi/4\\ \hline
 T&48r+420n+2&48r+420n+4\\
 d_T&470&472\\
 G&6G_B+534n+2G_{\mathcal R}+4c_n+37
   &6G_B+534n+2G_{\mathcal R}+4c_n+49
 \end{array}\\
 h\le\max(r,4n),\quad a_{\rm clean}=1,
 \quad\varepsilon_U=0
 \end{gathered}                                          \qquad(43)
```

The common Clifford constants allow 15 gates for the two orientations
of R cubed and R-dagger, six for controlled P_0, twelve for the P_0
factors in controlled F and its inverse, and four signal Hadamards.
The seeded word adds two six-Clifford one-T rotations. Replacing that
seed by the complex-axis native word R_a S_a adds four more Cliffords
for S and its actual inverse, with all other bounds unchanged. The returned
private pool is reused between the three controlled phase calls.
All helpers return exactly; no clean return of the source signal is
assumed before the paid recovery. The budgets fit polynomial-overhead
fusion. An arbitrary real theta below is an algebraic comparison,
not an assertion of exact native synthesis for every real angle.

The unseeded word has a protected sector for **every** unitary E.
If Fw=w then B_F and its inverse fix $`|0\rangle_aw`$;
the middle controlled E also fixes that vector. For any nonzero
address choose j different from zero with (v_x)_j=1 and set
$`w=|z=0,t=e_j\rangle`$. Then Fw=P_0w=w but Q_xw=-w.
Consequently

```math
 K_0w=\omega w/\sqrt2,\qquad
 \|G_0-G_*(x)\|\ge1/2\quad(x\ne0).                      \qquad(44)
```

Changing E, its spectrum or its mixer cannot repair this unseeded
sandwich while retaining that common fixed sector of F and P_0.
This explains why merely adding an inverse controlled branch does
not create a lookup source.

### Every real-axis signal seed retains a trace excess

The seed removes the premise that the middle control remains zero,
so (44) alone does not settle it. A separate calculation does. Fix a
nonzero address, use normalized trace t, and put

```math
 \begin{aligned}
 \mu&=\mathrm t(E)=a^{n-p},& a&=(2+\sqrt2)/4,\\
 \rho&=2^{-(n-p)},& A_2&=\mathrm t(M^2)=(1+\rho)/2,\\
 M&=(F+F^\dagger)/2,&
 Z&=\mathrm t(EME^\dagger M),\\
 J&=\mathrm t\!\left(M(E^\dagger-I)P_0(E-I)\right),&
 q&=\sin\theta,\quad b=\cos\theta.
 \end{aligned}
```

The symbol J here is a scalar moment, distinct from clean insertion
in the general Gram interface. The actual E has constant computational
diagonal mu and commutes with X_z0. Both M and P_0 anticommute with
X_z0, while $`\mathrm t(MP_0)=\mu`$. These facts and (42) give

```math
 \mathrm t(G_\theta)
 =\frac14+\frac{1-\mu}{4}
   +\frac{qbJ-b^2[Z+(1-2\mu)A_2]}8.                    \qquad(45)
```

For details of the cancellation, $`L^\dagger L=(I-bM)/2`$ and
$`C^\dagger C=(I-qP+b(I-P)M)/2`$. The normalized trace of
the linear term in (42) is
$`(\mu-1)v(u\mu+v)`$. In the quadratic term, the three
contractions against I, P and (I-P)M respectively are

```math
 \begin{aligned}
 \mathrm t((I-bM)\Delta^\dagger\Delta)&=2(1-\mu),\\
 \mathrm t((I-bM)\Delta^\dagger P\Delta)&=1-\mu-bJ/2,\\
 \mathrm t((I-bM)\Delta^\dagger(I-P)M\Delta)
   &=-\mu(1-\mu)-\frac b2[Z+(1-2\mu)A_2].
 \end{aligned}
```

These identities prove (45) without assuming any extra commutation.
Cyclicity also gives
$`J=\mathrm t(EME^\dagger P_0)+(1-2\mu)\mu`$.
Cauchy--Schwarz therefore implies
$`|Z|\le A_2`$ and
$`|J|\le\sqrt{A_2}+|1-2\mu|\mu`$.
For p at least three, mu is less than 1/16, A_2 is less than 9/16
and its square root is less than 3/4. Using $`|qb|\le1/2`$
and $`b^2\le1`$ in (45) proves the uniform bound

```math
 \begin{aligned}
 \|G_\theta-G_*(x)\|
 &\ge\mathrm t(G_\theta)-\frac14\\
 &\ge\frac{(1-\mu)(1-A_2)}4
           -\frac{\sqrt{A_2}+\mu}{16}
   \gt\frac{53}{1024},\qquad p\ge3,
       \quad x\ne0,
       \quad\theta\in\mathbb R.
 \end{aligned}                                                   \qquad(46)
```

The strict rational certificate is

```math
\frac{(15/16)(7/16)}4-\frac{3/4+1/16}{16}=\frac{53}{1024}.
```

### Complex-axis seeds obey the same bound

In fact every one-qubit unitary seed V in place of R_y(theta) is
excluded in the same range. Write
$`V^\dagger|1\rangle=\alpha|0\rangle+\beta|1\rangle`$,
and put

```math
 a=2\mathrm{Re}(\overline\alpha\beta),\qquad
 b=2\mathrm{Im}(\overline\alpha\beta),\qquad
 c=|\beta|^2-|\alpha|^2,\qquad a^2+b^2+c^2=1.
```

The polynomial (42) remains valid with
$`u=(\alpha+\beta)/2`$, $`v=(\alpha-\beta)/2`$.
Let $`N=(F-F^\dagger)/(2i)`$ and $`A=cM-bN`$.
Then $`L^\dagger L=(I-A)/2`$ and
$`C^\dagger C=(I-aP+(I-P)A)/2`$. Since
$`\mathrm t(MN)=0`$ and $`\mathrm t(NP_0)=0`$, define

```math
 \begin{aligned}
 A_{2,V}&=\mathrm t(A^2)
      =c^2 A_2+b^2(1-\rho)/2\le(1-a^2)A_2,\\
 Z_V&=\mathrm t(EAE^\dagger A),\\
 J_V&=\mathrm t(A\Delta^\dagger P_0\Delta)
      =\mathrm t(EAE^\dagger P_0)+(1-2\mu)c\mu.
 \end{aligned}
```

The same three contractions now give the complete normalized trace

```math
 \mathrm t(G_V)=\frac14+\frac{1-\mu}{4}
   +\frac{aJ_V-Z_V-(1-2\mu)A_{2,V}}8.                   \qquad(47)
```

For p at least three, Cauchy--Schwarz gives
$`|Z_V|\le A_{2,V}`$ and
$`|aJ_V|\le(\sqrt{A_2}+\mu)/2`$, because
$`|a|\sqrt{1-a^2}\le1/2`$ and $`|ac|\le1/2`$.
Also $`Z_V+(1-2\mu)A_{2,V}\le2(1-\mu)A_2`$.
Substitution in (47) proves the same strict bound 53/1024 in (46)
for **every one-qubit seed**, including complex axes.

These theorems use the actual flat, X-diagonal first branch and the
actual relative phase F. They do not allow arbitrary conjugations
of those operators, unequal outer controlled branches or different
acceptance words. They are not general three-band or higher-band
depth bounds. Address-controlled signal seeds are allowed fiber by
fiber; genuine address motion is not asserted for this trace proof.
Spectator identities and external reference identities preserve all
estimates. If an emitted source is within full error epsilon of
(41), subtract at most 2 epsilon from either Gram lower bound.

Exact p=1 fixtures multiply the complete signal word, compare (42)
on every carrier column, and check the following normalized traces
and wrong-sector receipts at x=1,z=e_1,t=0:

| Seed | Normalized trace of the Gram | Wrong-sector accepted norm squared |
| --- | --- | --- |
| theta=0 | $`9/32-\sqrt2/64`$ | $`1/2`$ |
| theta=pi/4 | $`11/32-5\sqrt2/128`$ | $`7/16+\sqrt2/16`$ |
| Native complex seed R_a S_a | $`41/128-\sqrt2/32`$ | $`1/2`$ |

The p=1 fixtures do not extend (46) to p=1 or p=2. The unseeded
bound (44) already covers those sizes analytically.

<a id="unequal-outer-branch-source"></a>

## 12. Unequal outer branches and their complete rejection condition

Changing only the left outer branch is a further actual source family.
Keep U_0 and B_F from (41), and use the paid partial Hadamard S of
Section 7 to put $`G=S^\dagger FS`$. For any signal seed V, define

```math
 B_G=H_a\mathrm C_a(G)H_a,\qquad
 U_{G,V}=U_0 B_G V_a^\dagger\mathrm C_a(E)V_aB_F^\dagger.
                                                               \qquad(48)
```

Here G_0=F_0=P_0, so the exact base Gram is still P/2. The
left branch and the right inverse are now different. For the actual
partial-H G, its controlled word adds two S calls, or 4(n-1)
Cliffords, to the ledger (43). In particular V=I has

```math
 \begin{aligned}
 T&\le48r+420n+2,&d_T&\le470,\\
 G_{\rm native}&\le6G_B+534n+2G_{\mathcal R}+4c_n+37+4(n-1).
 \end{aligned}
```

The native V=R_a seed adds two T, two T layers and twelve Cliffords;
V=R_a S_a adds four further Cliffords in the common inverse bound.
The private width remains at most max(r,4n), with one clean signal.
All private work returns exactly, actual inverses are emitted, and
the exact polynomial intermediate budgets still fit fusion. A general
unitary V or G in the theorem below is algebraic freedom, not a free
native implementation.

The following complete polynomial also applies to any unitary left
branch G. Write $`V^\dagger|1\rangle=\alpha|0\rangle+\beta|1\rangle`$
and put

```math
 \begin{aligned}
 u&=(\alpha+\beta)/2,&v&=(\alpha-\beta)/2,&
 L_F&=uI+vF,&\Delta&=E-I,\\
 C&=I-P,&D&=\sqrt2P-C,\\
 X&=I/2+u\Delta L_F^\dagger,&
 Y&=F^\dagger/2+v\Delta L_F^\dagger.
 \end{aligned}
```

The symbols X,Y,D in this subsection are block coefficients, not
Paulis or the original phase D_x. Multiplication of the two signal
blocks gives the entire accepted operator and Gram:

```math
 \begin{aligned}
 K_{G,V}&=\omega(CX+DGY),\\
 K_{G,V}^\dagger K_{G,V}
 &=X^\dagger CX-X^\dagger CGY-Y^\dagger G^\dagger CX
       +Y^\dagger G^\dagger(I+P)GY,\\
 X^\dagger X+Y^\dagger Y&=I/2.
 \end{aligned}                                                   \qquad(49)
```

The last identity follows from the unit norm of the first block
column before B_G. It also follows directly from the definitions.
These equations retain all off-diagonal action, coherent addresses,
arbitrary bank and dirty inputs, and reference identities.

With V=I, any unit vector satisfying Fw=Pw=w has

```math
 K_{G,I}w=\omega DGw/2,\qquad
 \|K_{G,I}w\|^2=(1+\|PGw\|^2)/4\ge1/4.              \qquad(50)
```

The wrong-charge basis vector in (44) has these properties. Thus
the unseeded unequal word fails for every E and every unitary G,
even if G moves addresses or mixes arbitrary spectators. For the
specified partial-H G, which commutes with P, its acceptance remains
one half. An unequal left branch alone does not remove the protected
input.

### A seed-independent rejection obstruction with the actual right branch

Allow any unitary seed V and any unitary G, but retain the actual E,F
of (3). For every fixed nonzero address and every p at least two,

```math
 \left\|K_{G,V}^\dagger K_{G,V}-G_*(x)\right\|\gt2^{-22}.
                                                               \qquad(51)
```

This permits arbitrary changes of the left branch, beyond conjugation
of F. It is not a theorem about every three-branch source: its actual
right branch, middle branch, placement of the single seed and U_0
readout remain hypotheses.

To prove (51), exact rejection of a vector w would require
GYw=CXw. More quantitatively, if w has wrong lookup charge and the
complete Gram error is at most delta, then (49) implies

```math
 \begin{aligned}
 \|K_{G,V}w\|^2
 &=2\|PGYw\|^2+\|C(X-GY)w\|^2\le\delta,\\
 d(w)&=\|Yw\|^2-\|CXw\|^2,\qquad
 |d(w)|\le\sqrt\delta.
 \end{aligned}                                                   \qquad(52)
```

For the last bound, the distance from GYw to CXw is at most
$`\sqrt\delta`$, while their two norms sum to at most one by the last
line of (49). This removes the left branch from a necessary condition.
Since the argument holds for every wrong-charge unit vector, it also
gives the operator condition

```math
 \|\Pi_{\rm neq}(Y^\dagger Y-X^\dagger CX)
                    \Pi_{\rm neq}\|\le\sqrt\delta.
```

Two computational inputs suffice. First take the earlier
w_+=|z=0,t=e_j>, with j different from zero and (v_x)_j=1.
It obeys Fw_+=w_+, P_0w_+=w_+, and Q_xw_+=-w_+.
For w_-, choose z with B_xz of weight four, set z_0=1 and let t be
z with coordinate zero set to zero.
This is possible because B_x is independent of z_0 and has the full
first n-p coordinates as image, with n-p at least six when p is at
least two. Then h_x=4, P_0=-1, F=+1, and
v_x dot (z+t)=v_x dot e_0=1, so Q_x=-1 as well.

The actual E is X-diagonal with constant computational diagonal
mu from (36). Its autocorrelation kappa is also common to every
computational basis column w:

```math
 \langle w|E^\dagger P_0E|w\rangle
 =\sigma\kappa\quad(P_0w=\sigma w),\qquad
 0\le\kappa\le1.                                      \qquad(53)
```

This assertion concerns computational basis vectors, not arbitrary
P_0 eigenvectors. Indeed $`E^\dagger P_0 E P_0`$ is X-diagonal,
so every computational
diagonal equals its normalized trace. Conjugating that trace by the
linear Fourier mixer turns P_0 into a fixed bit translation. In the
independent phase coordinates $`a=t+P_{\rm del}z`$ and b=B_xz,
where P_del deletes coordinate zero, the exponent is
the sum of (-1) raised to a_j times b_j over n-p active coordinates.
The contribution to its phase autocorrelation from one coordinate is
1 for translation (0,0), 1/2 for (1,0), and 1/sqrt(2) for either
translation with the second bit one. The total autocorrelation is
their product, proving (53). No cancellation or independence of the
chosen witness is assumed beyond this uniform-coordinate trace.

Put $`a=2\mathrm{Re}(\overline\alpha\beta)`$,
$`c=|\beta|^2-|\alpha|^2`$, q=1-mu, and b=(1-c)(1+a).
Here $`a^2+c^2\le1`$, including complex-axis seeds.
Using Fw_+=w_+ and Fw_-=w_-
in the definitions of X,Y gives the two exact defects

```math
 \begin{aligned}
 d_+&=1/4+acq/4-b(1-\kappa)/16,\\
 d_-&=acq/2+b(1-\kappa)/16.
 \end{aligned}                                                   \qquad(54)
```

For p=2, $`\mu=((2+\sqrt2)/4)^6\gt3/8`$. Hence
d_++d_-=(1+3acq)/4 is greater than 1/64, using ac at least
-1/2. Equation (52) already gives $`\delta\gt2^{-14}`$ in this case.

For p at least three, mu is less than 1/24: the single-coordinate
moment is less than 6/7, n-p is at least 21, and (6/7)^21 is less
than 1/24. Suppose $`\delta\le2^{-22}`$ and write
$`\epsilon=\sqrt\delta`$, so epsilon is at most 1/2048. Equations
(52)--(54), together with 1-kappa at most one, imply

```math
 \begin{aligned}
 ac&\ge-8/23-64\epsilon/23,\\
 b&\ge8/3-16\epsilon,\\
 a-c=b-1+ac&\ge91/69-432\epsilon/23,\\
 a^2+c^2
 &\ge(91/69-432\epsilon/23)^2-16/23-128\epsilon/23
 \gt1.
 \end{aligned}                                                   \qquad(55)
```

The last expression decreases throughout this epsilon interval and
is still greater than one at its right endpoint, as the exact rational
fixture checks. This contradicts the seed normalization and proves
(51). The constant is a conservative certificate, not an optimized
extraction or resource constant.

The same two witness inputs may be embedded into the full coherent
address space and tensored with any spectator state. Thus (51) is
also a global bound if the left G moves addresses or mixes spectators,
while E,F and the scalar signal seed retain the specified form.
An address-controlled seed is covered fiber by fiber. A general
seed acting on bank or dirty inputs is not included. Full source
approximation error eta decreases the bound by at most 2 eta.

The actual partial-H unequal source and its seeds are therefore
paused. Reopening requires changing the right branch or middle
dynamics, changing the signal/readout architecture, or another source
interface. The theorem does not supply a general depth lower bound.

### A concrete moving-address right branch does not preserve the base

One changed right branch was also tested on its complete operator.
Let A apply H to all p address bits and replace the right F in (48)
by $`F^{\rm addr}=A^\dagger F A`$, with left G=F, the actual E,
and native seed V=R_a. All operators here act on the full coherent
address space. The right controlled word is implemented by two A
calls around the paid controlled F inverse, adding only 2p Cliffords
to the seeded ledger (43). Its count, depth, width, inverse and private
return budgets remain qualified intermediate budgets.

This changes a hypothesis of (51), but it does not inherit the old
base normalization. At p=1 the complete 64-dimensional source has
the following exact receipt:

```math
 \|K|x=0,z=e_1,t=e_0\rangle\|^2
   =7/32-13\sqrt2/128\gt0,\qquad
 G_*|x=0,z=e_1,t=e_0\rangle=0.                          \qquad(56)
```

The fixture checks the full accepted polynomial (49), with this new
right branch, rather than restricting an intermediate address fiber.
This is a failure of the specified moving-address word, not a general
address-motion obstruction. A replacement architecture must explicitly
pay and verify its normalization instead of assuming that a fixed
address-zero action survives the Hadamards.

<a id="tail-bank-commutant-obstruction"></a>

## 13. Joint signals, a common bank Pauli, and a changed bank coupling

Two signals can interact physically without supplying new lookup sectors.
For $`A=E_xX_{z_0}`$, consider the actual two-clean word

```math
U_{\rm CZ}=H_aH_b\,\mathrm{CZ}_{ab}S_a\,
       \mathrm C_b(F_x)\mathrm C_a(A)\,H_aH_b.
                                                               \qquad(57)
```

Products act rightmost first. Both signals are inserted and accepted at
zero, after their CZ coupling. The four signal paths give the complete
accepted block and Gram

```math
\begin{aligned}
K_{\rm CZ}&=\frac{I+iA+F-iFA}{4},\\
K_{\rm CZ}^\dagger K_{\rm CZ}
 &=\frac1{16}(I+iA+F-iFA)^\dagger(I+iA+F-iFA).
\end{aligned}                                                   \qquad(58)
```

At x=0, A=X_z0 anticommutes with F=P_0, so
$`K_{\rm CZ}=(I+iX_{z_0})(I+P_0)/4`$ and the Gram is exactly
$`(I+P_0)/4`$. The signal coupling has therefore retained the base
normalization. Nevertheless the accepted polynomial is still an
independently attenuated one-signal core: its complete coefficients,
not the presence of a physical CZ, determine that fact. Explicitly,
$`K_{\rm CZ}=\langle0|W_a|0\rangle/\sqrt2`$ for
$`W_a=H_a\mathrm C_a(F)H_a\mathrm C_a(A)S_aH_a`$.
On any
wrong-charge F=+1 vector w, its accepted diagonal entry is one half,
since the term $`(I-F)A`$ has zero matrix element against w. Hence
the error is at least one quarter for this word, independently of A.

The controlled A uses a signal-to-z_0 CNOT followed by the paid
controlled E. The signal-only skeleton has eight Cliffords forward
and ten in its actual reverse, where S-dagger is three S gates.
Including controlled X and the six Cliffords in controlled P_0,
common forward/actual-inverse bounds are

```math
\begin{aligned}
T&\le32r+280n,&d_T&\le312,\\
G&\le4G_B+356n+2G_{\mathcal R}+4c_n+17,&
h&\le\max(r,4n),\qquad a_{\rm clean}=2.
\end{aligned}                                                   \qquad(59)
```

All private helpers return exactly on arbitrary inputs. The two source
signals need not return before acceptance; complete recovery is a separate
operation. Every branch, its literal phase and its actual inverse are
native, with zero synthesis error here. The intermediate polynomial
budgets fit the retained fusion allowance.

### The actual E,F branch algebra has a charge-flipping commutant

There is a simpler obstruction that does not count signal paths or
native bands. Put $`j_*=n-p`$, the first binary coordinate of the last
field block in t, and let $`W=X_{t_{j_*}}`$. The image of B_x occupies
only the first n-p binary target coordinates. The last p terms in the
weight difference (1) therefore cancel, making D_x independent of all
last-block t bits. Since j_* is positive, P_0 also acts trivially on
this bit. The actual E is X-diagonal. Consequently, at every address,

```math
[W,D_x]=[W,E_x]=[W,F_x]=[W,P_0]=[W,X_{z_0}]=0.
                                                               \qquad(60)
```

For a nonzero field element x, the last field entry of v_x is
$`x^{m-1}=1`$. In the retained binary basis, field one has first
coordinate one. Thus $`(v_x)_{j_*}=1`$, and
$`WQ_xW=-Q_x`$ at every nonzero address.

Consider any word of the completed controlled branches in (60), their
actual inverses, and interlayers acting as identity on this principal
bank qubit. The signal bank, all other work and the address may be
arbitrarily large; the interlayers may move the address coherently.
More generally any additional completed branch commuting with W is
allowed. Signal preparation, helper insertion and joint acceptance must
also act as identity on this bank qubit. The whole source and its
compression then commute with W. In particular its accepted Gram G
satisfies WGW=G, whereas the full target obeys

```math
\|G_*-WG_*W\|=1/2,\qquad
\|G-G_*\|\ge1/4\quad(p\ge1).                           \qquad(61)
```

The first norm is attained in any nonzero-address fiber; at x=0 the
target difference vanishes. The triangle inequality and unitary
invariance prove the second norm without restricting the source to
fixed-address fibers. Arbitrary spectator and reference identities do
not change the conclusion. Full source approximation error eta can
reduce the bound by at most $`2\eta`$.

This is an all-call-count limitation of the displayed branch algebra,
not a lower bound on native T-depth. Individual gates inside a paid
controlled-D word need not commute with W: it is the completed branch
identity, with returned helpers, that enters (60). Additional clean
signals, scalar signal seeds, higher-degree words in E,F or address-only
motion cannot remove this commutant. A branch or coupling acting on the
blind bank degree of freedom must change. General left unitaries in
Section 12 need not commute with W, so (61) does not replace that
section's differently scoped theorem.

### A tail-Hadamard change leaves the grammar but does not prove the Gram

Let H_tail apply H to every last-block t bit and put
$`E'=H_{\rm tail}E H_{\rm tail}`$. The x=0 branch remains identity.
This changes the actual bank algebra, not only signal processing; at
p=1 the exact fixture confirms that E' no longer commutes with W.
To change the accepted polynomial as well, insert a controlled E-prime
inverse after controlled F:

```math
U_{\rm tail}=H_aH_b\,\mathrm{CZ}_{ab}S_a\,
 \mathrm C_a(E'^\dagger)\mathrm C_b(F)
 \mathrm C_a(E'X_{z_0})\,H_aH_b.
                                                               \qquad(62)
```

Its four path operators are $`I,X_{z_0},F,E'^\dagger FE'X_{z_0}`$.
With $`F'=E'^\dagger FE'`$, the complete accepted block is

```math
K_{\rm tail}=\frac{I+F+iX_{z_0}-iF'X_{z_0}}4.            \qquad(63)
```

No commutation of F with E-prime is used. Both F and F-prime
anticommute with X_z0. Writing $`M=F+F^\dagger`$ and
$`M'=F'+F'^\dagger`$, the full Gram is

```math
K_{\rm tail}^\dagger K_{\rm tail}
 =\frac{4I+M+M'+i(M-M'+F'^\dagger F-F^\dagger F')X_{z_0}}{16}.
```

This expression retains every off-diagonal entry.
At x=0, (63) again equals $`(I+iX_{z_0})(I+P_0)/4`$, so its
normalization is correct. Relative to (59), the additional controlled
E-prime inverse and the two tail-Hadamard conjugations give

```math
\begin{aligned}
T&\le48r+420n,&d_T&\le468,\\
G&\le6G_B+534n+4G_{\mathcal R}+4c_n+17+4p,&
h&\le\max(r,4n),\qquad a_{\rm clean}=2.
\end{aligned}                                                   \qquad(64)
```

These are exact common forward/reverse bounds, with the same complete
private return contract. A full p=1 fixture nevertheless gives

```math
\|K_{\rm tail}|x=1,z=0,t=e_1\rangle\|^2=1/2,
\qquad Q_1|0,e_1\rangle=-|0,e_1\rangle.                 \qquad(65)
```

Thus breaking (60) is necessary for this branch route, not sufficient
for lookup. The p=1 fixture also verifies that the accepted block and
its Gram no longer commute with W, so this is not a disguised
application of (61). Equation (65) is a finite failure of this specified word;
it is not an all-size exclusion of the changed bank algebra, nor an
obstruction to correcting finitely many small sizes. A further source
must establish the full sector identity, rather than infer it from the
absence of the old commutant.

## 14. Verification scope and the next unpaid operation

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
acceptance polynomial. The three-branch sandwich (41) is also paused:
its unseeded control has a protected wrong sector, and every one-qubit
seed has the scoped trace excess (46)--(47). The unequal left branch
in (48) changes those moments but is independently stopped by (50)--(55),
even with an arbitrary one-qubit seed and arbitrary left unitary.
A changed right branch or different branch algebra remains unpaid.
Changing only the first mixer, adding a joint
final coupling, or moving addresses cannot overcome (40). The earlier
suggestion to escape solely by mixing the distinguished coordinate is
therefore superseded. The first unpaid operation remains a normalized
complete accepted Gram on the correct lookup sectors; none of these
necessary changes supplies that operation by itself.

The joint-signal words (57) and (62) now have complete polynomials,
native ledgers and exact base tests. The former retains independent
attenuation; the latter breaks the blind tail-X commutant but fails its
displayed p=1 lookup test. The all-call-count result (61) means that
changing signal processing alone cannot repair the retained E,F algebra.
The next source must change actual bank dynamics and establish the
complete lookup-sector identity; merely removing one old witness does
not establish that identity. Recovery and bounded-clean promotion must
be paid together with any successful source. None of these results
changes the complete-frame frontier.
