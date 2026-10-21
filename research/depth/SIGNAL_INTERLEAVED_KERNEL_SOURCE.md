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

## 6. Exact fixtures and the next changed coupling

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

A concrete next coupling changes the eigenspaces of the **final**
controlled phase. Let S apply H to z_j and t_j for every j at least
one, leaving z_0,t_0 unchanged, and replace B in (11) by

```math
 B'_x=iS^\dagger F_xS.                                   \qquad(20)
```

S commutes with P_0, so the exact x=0 action (14) is retained.
The literal controlled operation is
$`S^\dagger\mathrm C_a(B)S`$, costing 4(n-1) extra Cliffords
and no extra T gates, width, clean inputs or precision. Equation (12)
still gives the actual accepted block with B replaced by B'; it makes
no lookup-sector claim. The first concrete check is the intersection
$`S^\dagger\ker(F_x-I)\cap\mathrm{ran}\Pi_{\rm neq}(x)`$ and
whether its vectors can satisfy (19). If that test leaves a candidate,
its whole Gram, including all off-diagonal entries and every return,
must be derived. This particular final-coupling change is untested here.
It supersedes additional sign or fixed-mixer scans of the rejected words.
