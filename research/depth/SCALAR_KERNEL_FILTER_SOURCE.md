# A normalized scalar-zero source on a thin part of the kernel charge

[Canonical pencil and paired Gram](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) ·
[Paid scalar zero reflection](EXPANDER_ZERO_REFLECTION.md) ·
[Changed signal-interleaving question](CROSSED_CUBIC_KERNEL_SOURCE.md#124-charge-preservation-requires-a-different-acceptance-polynomial)

A scalar zero test gives an actual non-Gaussian, one-clean source with
a paid full-operator approximation to a normalized ideal Gram. That
ideal Gram rejects the entire wrong-charge sector but accepts only a
thin subset of the correct sector. Its exact zero eigenspace
contains both charges, so ordinary spectral processing of its accepted
Gram cannot complete the source. This is a paid partial construction,
not a complete indicator or a frame-depth improvement.

## 1. The returned kernel query and the scalar XOR

Use the canonical standard dual coordinates with $`p\ge1`$,
$`m=2^p`$, and $`n=mp`$. Let P delete only the first binary coordinate,
and put

```math
B_x=\bar N_x^{\mathsf T}-I,\qquad
K_x=(I+B_x)P=P+B_x.
                                                               \qquad(1)
```

The second equality uses $`B_xP=B_x`$. The
[corank-one query](CROSSED_CUBIC_KERNEL_SOURCE.md#121-the-actual-phase-and-native-schedule)
is

```math
Q_K|x,z,t\rangle=|x,z,t+K_xz\rangle,
\qquad
\ker K_x=\mathrm{span}\{e_0\},\qquad
\mathrm{im}K_x=\ker(v_x^{\mathsf T}).                    \qquad(2)
```

Both n-bit banks are arbitrary. In particular z is a coherent dirty
translation input, not a prepared choice of kernel element. Equation
(2) deletes one binary coordinate, not the first field block. The
query is the paid nonlinear shear $`Q_{B_x}`$ and the fixed CNOT shear
$`Q_P`$. With the existing notation

```math
r=(m-1)p^2,\qquad s\le(m-1)p^3,
```

safe common forward/actual-inverse bounds are

```math
T_K\le8r,\qquad d_K\le8,\qquad h_K\le r,\qquad
G_K\le3n+16r+2s.                                      \qquad(3)
```

Every query helper returns exactly, on every coherent address and bank
input. All its nonlinear T targets are private helpers.

Let $`R_0^{(n+1)}`$ be the reflection with sign minus precisely when
all n target bits t and a signal bit a are zero. The scalar theorem
supplies a zero-clean native $`\widetilde R_0`$ with

```math
\|\widetilde R_0-R_0^{(n+1)}\otimes I_W\|\le\epsilon,
                                                               \qquad(4)
```

including every arbitrary private work input and reference. Denote its
common forward/actual-reverse ledger by $`T_0,G_0,d_0,h_0`$. At fixed
epsilon these are $`T_0,G_0,h_0=O_\epsilon(n)`$ and
$`d_0=O_\epsilon(1)`$. There is no fine-angle or prepared work promise.
Define the emitted complete XOR approximation

```math
\widetilde Q_0=H_aX_a\widetilde R_0X_aH_a.              \qquad(5)
```

Its ideal operator toggles a exactly when t is zero, with no additional
phase. On that target sector $`R_0^{(n+1)}=-Z_a`$;
$`X_a(-Z_a)X_a=Z_a`$ and $`H_aZ_aH_a=X_a`$. On every other target
sector the ideal operator is identity. Thus the full error in (5) is
at most epsilon, with ten additional elementary Cliffords: two H gates
and two literal $`X=HSSH`$ words.

<a id="normalized-scalar-kernel-source"></a>

## 2. An actual signal-interleaved source and its complete Gram

Use the following literal products, with rightmost factors acting
first:

```math
\begin{aligned}
\widetilde Q_f&=Q_K^\dagger\widetilde Q_0Q_K\widetilde Q_0,\\
\widetilde M&=X_a\widetilde Q_f,\qquad R=T_Y(a),\\
\widetilde U&=R^3\widetilde M R^\dagger\widetilde M^\dagger.
\end{aligned}                                                   \qquad(6)
```

The last factor is the **actual adjoint** of the emitted M. Approximate
zero queries need not be involutions and are not replaced by themselves
in an inverse occurrence. No controlled source or controlled kernel
query is being assumed. The signal is directly the scalar XOR output.
The Pauli-axis T is the fixed Clifford conjugate taking Z to Y.

In the ideal word, the query toggles the Boolean function

```math
f_x(z,t)=[t=0]\mathbin\oplus[t+K_xz=0],                 \qquad(7)
```

so M acts as $`X_a^{1\oplus f_x(z,t)}`$ while preserving all other
inputs. Let $`\Pi_f`$ denote the complete diagonal projector onto
$`f_x(z,t)=1`$, with identity on every private helper. With
$`\omega=e^{i\pi/4}`$, the full ideal signal operator is

```math
U=\Pi_f\otimes R^2+(I-\Pi_f)\otimes\omega^{-1}Y_a.
                                                               \qquad(8)
```

Indeed the M-identity branch gives $`R^2`$, and
$`R^3X_aR^\dagger X_a=\omega^{-1}Y_a`$ on the other branch.
The scalar $`\omega^{-1}`$ is retained. For J inserting only
$`a=0`$, the complete accepted block and Gram are therefore

```math
A=J^\dagger UJ=\frac{\omega}{\sqrt2}\Pi_f,\qquad
A^\dagger A=\frac12\Pi_f.                              \qquad(9)
```

This is an operator identity on both arbitrary banks, the full coherent
address, all private dirty inputs, and their reference extensions. It
is not an average over the translation bank z or a scalar trace.

Two zero-query errors give
$`\|\widetilde M-M\|\le2\epsilon`$. Using the actual adjoint
in (6), unitary telescoping then gives

```math
\|\widetilde U-U\|\le4\epsilon,\qquad
\|\widetilde A^\dagger\widetilde A-\tfrac12\Pi_f\|
\le8\epsilon.                                        \qquad(10)
```

The second bound uses that both accepted blocks are contractions.
All approximate dirty return is included in the first norm. These are
accuracy statements for the partial projector in (9), not for the
complete lookup projector.

## 3. Paired physical lookup and one simultaneous native ledger

Let C,D be the arbitrary physical carriers in the pencil coordinates,
and let E be its fixed Vandermonde binary matrix. Put

```math
S|C,D\rangle=|C,C+D\rangle,\qquad
C_E=U_{E^{-\mathsf T}}\otimes U_{E^{-\mathsf T}},\qquad
B=C_ES.
```

Thus the standard coordinates in (6) are
$`z=E^{-\mathsf T}C`$ and $`t=E^{-\mathsf T}(C+D)`$.
The charge on t is exactly the paired physical charge:

```math
(-1)^{v_x^{\mathsf T}t}
 =\Omega_{x,C}\Omega_{x,D},\qquad
\Pi_{\rm eq}=(I+\Omega_{x,C}\Omega_{x,D})/2.
```

Emit $`\widetilde U_{\rm phys}=B^\dagger\widetilde U B`$.
Its accepted Gram is within $`8\epsilon`$ of
$`B^\dagger\Pi_f B/2`$. Both B orientations are paid: the basis
changes use at most $`4c_n`$ CNOTs, with $`c_n=n^2+3n`$, and the two
carrier-XOR walls use $`2n`$ more. They are not initialized-bank copies.

The source and its actual inverse have the simultaneous safe bounds

```math
\begin{aligned}
T_U&\le4T_K+4T_0+2,\\
D_T(U)&\le4d_K+4d_0+2,\\
G_U&\le4G_K+4G_0+63+4c_n+2n,\\
h_U&\le\max\{r,h_0\},\qquad a_{\rm clean}=1.
\end{aligned}                                                   \qquad(11)
```

There are four kernel-query appearances and four scalar-zero appearances.
Their physical T bands are serial. The scalar-to-XOR walls contribute
40 Cliffords, the two X factors in M and its reverse contribute eight,
and fifteen is a safe common forward/reverse allowance for the two
signal rotations. Each rotation needs one T or T-dagger, including the
$`R^3`$ factor, whose remaining power is Clifford.

The private pool can be reused: every primitive is certified on its
entire arbitrary-work space, so an earlier approximate return does not
invalidate a later bound. The final error still follows from full
operator telescoping, not from an assumption of exact intermediate
return. The two principal carriers additionally allocate
$`(2p-1)m`$ wires beyond the specified m-wire output bank. At each fixed
epsilon, (11) has constant native T-depth and polynomial address
overhead in count and dirty width. These budgets fit the
[polynomial-overhead fusion interface](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md),
but its sector contract does not.

<a id="scalar-kernel-sector-failure"></a>

## 4. Correct charge, missing support and the exact gap

Every vector $`K_xz`$ lies in the even-charge hyperplane. Hence
$`f_x(z,t)=1`$ implies $`v_x^{\mathsf T}t=0`$: either t is zero
or it equals $`K_xz`$. The wrong-charge sector is rejected exactly
in the ideal source. Nevertheless the accepted support is too small.
For each of the $`2^n-2`$ coefficient vectors outside the two-element
kernel of K, precisely the two targets $`t=0,K_xz`$ are accepted.
For its two kernel vectors no target is accepted. Thus, at each fixed
address and excluding spectator factors,

```math
\mathrm{rank}\Pi_f=2^{n+1}-4,\qquad
\mathrm{rank}\Pi_{\rm eq}=2^{2n-1},\qquad
B^\dagger\Pi_fB\lt\Pi_{\rm eq}.                       \qquad(12)
```

The strict inequality also has a direct witness: $`z=t=0`$ has even
charge but f=0. Consequently

```math
\left\|\tfrac12 B^\dagger\Pi_fB-\tfrac12\Pi_{\rm eq}\right\|
 =\frac12,\qquad
\left\|\widetilde A_{\rm phys}^\dagger\widetilde A_{\rm phys}
             -\tfrac12\Pi_{\rm eq}\right\|
 \ge\frac12-8\epsilon.                                \qquad(13)
```

This failure holds for every address size, not merely for a small
fixture. Both an odd-charge vector and the displayed missing-even
vector have Gram eigenvalue zero. Any scalar spectral function of the
ideal Gram assigns the same value to both, and therefore has error at
least one quarter against $`\Pi_{\rm eq}/2`$. If it keeps the odd
sector exactly zero, the missing-even error remains one half. This
statement does not exclude a changed source using additional
interleaving or information in its rejected block.

The general rank bounds for short words of thin reflections are already
proved in the [scalar reflection home](EXPANDER_ZERO_REFLECTION.md#11-a-paid-graph-reflection-still-needs-complete-readout).
They are not new here. The new constructive content is (6)--(11): a
normalized, paid, correctly charged partial source, together with its
specific support failure. Enlarging that support is an unpaid operation;
precision or four-call recovery alone cannot enlarge it.

<a id="fixed-diagonal-kernel-phase-rigidity"></a>

## 5. A fixed diagonal broad phase cannot keep the exact rejection property

One tempting repair is to replace the scalar zero phase by a broader
fixed phase. In physical carrier coordinates, let F be any
address-independent diagonal unitary, with entries $`F(t)`$. For a
coherent kernel translation h define its diagonal derivative

```math
D_h(t)=F(t+h)\overline{F(t)},\qquad h\in H_x,
\qquad H_x=\{h:h_{Y_x}=0\}.                            \qquad(14)
```

Suppose this derivative is identity for every arbitrary translation
input on every wrong-charge target: $`D_h(t)=1`$ whenever
$`t_{Y_x}=1`$, for every x and every $`h\in H_x`$. Then F is constant
on the whole affine halfspace $`t_{Y_x}=1`$. These halfspaces overlap,
so their constants agree. It follows that there is a unit scalar c such
that

```math
F(t)=c\quad\hbox{whenever any specified bank bit is one},
\qquad
\mathrm{rank}(F-cI)\le2^{n-m}.                         \qquad(15)
```

The unspecified carrier bits are arbitrary in every halfspace. They
do not provide an exception to (15). Only the subspace where all m
specified bank bits vanish can carry a different phase.

The same argument is robust. If
$`|D_h(t)-1|\le\eta`$ uniformly on those wrong-charge inputs, choose
one carrier string with all specified bank bits one. Every t outside
the all-zero-bank subspace lies in a common halfspace with that string.
Thus F is within eta, in operator norm, of a diagonal unitary that is
constant outside that subspace. The exceptional rank bound is unchanged.

This is a restriction on the explicit fixed-diagonal derivative
interface, not a no-go theorem for arbitrary non-Gaussian or
signal-interleaved sources. An address-dependent phase, a non-diagonal
insertion, or a different acceptance protocol must be priced and
analyzed separately. Within this interface, retaining exact rejection
of all wrong-charge inputs cannot also turn the scalar phase into a
constant-fraction insertion.

The route is paused. Its first missing operation is a paid interference
mechanism that accepts the entire correct-charge sector while
continuing to reject the other sector. The
[exact fixtures](../../tests/test_scalar_kernel_filter_source.py) check
all columns of a small literal native source and its actual reverse,
the corank-one support count, and the finite fixed-diagonal rigidity
condition. The small source uses an exact two-control zero XOR with a
returned arbitrary phase helper; it checks the ideal identities and
literal native phases, not the asymptotic approximate expander emitter.
Its seven-wire word has 66 T gates, 177 forward Cliffords and 179 for
the actual inverse, with one initialized signal and one arbitrary helper.
They support the identities; the all-size claims are analytic
proofs, internally audited rather than externally reviewed.
