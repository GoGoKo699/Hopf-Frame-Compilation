# Actual native center witness and safe-axis derivative pruning

This note gives a stronger sufficient condition for two-layer Boolean
derivatives, then an actual one-clean, three-native-band source where it
does not cover the remaining center direction. The target is CCZ, of
Boolean degree three. The witness limits these cancellation arguments;
it supplies no growing-degree family or unrestricted lower bound.

Throughout, $`T_P=\Pi_P^++\omega\Pi_P^-`$ for a signed Hermitian
Pauli P, $`\omega=e^{i\pi/4}`$, and products of operators act from right
to left. Native word lists below are chronological. All scalar phases
and the actual reversed inverse are retained.

## 1. Activating individually safe endpoint axes still gives a shallow derivative

Let the exact full logical operator be

```math
G=K_C^\dagger E_C^\dagger D E_AK_A=(-1)^{g(x)},
```

where each endpoint E is an independent native Pauli band, the K's are
actual Cliffords, and D is a commuting Pauli C3 operation. All logical
and dirty inputs remain arbitrary. Call an endpoint axis **safe** when
it commutes with D as an operator; let $`L_A^{\mathrm{bad}}`$ and
$`L_C^{\mathrm{bad}}`$ span the respective remaining axes.

Suppose a Pauli R activates only safe axes of E_C. The literal operator

```math
K_R=R E_CR E_C^\dagger
```

is a Clifford product of signed S rotations on just those axes, with
its scalar included. Each factor commutes with D. Since
$`E_CR E_C^\dagger=R K_R`$, it follows that

```math
D^\dagger(E_CR E_C^\dagger)D=(D^\dagger RD)K_R
\quad\text{is an actual Clifford}.                         \qquad(1)
```

For $`Q=K_C^\dagger RK_C`$, the full conjugate $`G^\dagger QG`$
therefore has an actual two-native-band representation using E_A and
its actual inverse. Multiplying by the Pauli Q extracts the Boolean
derivative in Q's computational X direction. The retained zero-clean
two-layer theorem bounds its degree by three. This constructed
derivative uses at most twice E_A's T-count, depth at most two, no
initialized helper, and paid Clifford words of size $`O(N^2)`$ on N
logical wires. The reverse statement uses E_C instead. The Clifford
in (1) is synthesized as a Clifford; D is not credited as a free gate.

Consequently the derivative argument of
[Computational phase derivatives §6](COMPUTATIONAL_PHASE_DERIVATIVES.md#6-a-broader-derivative-theorem-for-a-commuting-c3-sandwich)
improves to

```math
g=q_4+h,\qquad \deg q_4\le4,\qquad
h\text{ depends only on }
T_{\mathrm{bad}}=M\cap K_A^\dagger L_A^{\mathrm{bad}}K_A
                  \cap K_C^\dagger L_C^{\mathrm{bad}}K_C .         \qquad(2)
```

Here M is the computational Z parity space. For a target independent
of the dirty inputs, take the intersection with the designated-input
parity space as well. As in the earlier theorem, (2) is Boolean degree
accounting; removing q_4 is not a free native operation. In particular,
vanishing T_bad closes this subcase at sign degree four, while a
nonzero space need not support any phase dependence at all.

### 1.1. The retained invariant center can survive this pruning

In the genuine quadratic-control normal form, write

```math
D=H^{b(z)}E_B,\qquad S=R+L_B+\langle h\rangle,
```

where b has a nonzero quadratic part, R is the actual shared control
center, and E_B is the original independent native middle band. Its
axes commute with H and R; every endpoint axis commutes with R. Then
an endpoint axis a is safe exactly when it commutes with H and every
individual B axis, or equivalently when $`a\in S^\perp`$.

Indeed put $`X=E_BaE_B^\dagger`$ and
$`e=\langle a,h\rangle`$. Every Pauli term of X has this same
H-pairing, so $`DaD^\dagger=(-1)^{e b}X`$. If e=1 and this equals
a, then X is a times the non-Pauli Clifford $`(-1)^b`$. Its intrinsic
Clifford displacement lies in R. Native independence identifies this
displacement with the span of the B axes activated by a. Every such
axis would lie in R, contradicting the fact that a centralizes R.
Empty activation would make $`(-1)^b=I`$, also impossible. If e=0,
the equality X=a forces empty activation directly. The converse is
literal commutation.

Let R_0 be the simultaneous displacement center from
[Quadratic control native matching §5.3](QUADRATIC_CONTROL_NATIVE_MATCHING.md#53-a-larger-actual-center-gives-monomial-clifford-endpoint-images).
For each computational probe u and s in S, the native displacement
$`(\Phi_{J_u}-I)s`$ is a sum of A axes pairing nontrivially with s.
Only bad axes can occur. The same assertion holds at C because the
two probe actions agree on S. Thus

```math
R_0\subseteq L_A^{\mathrm{bad}}\cap L_C^{\mathrm{bad}}.                \qquad(3)
```

After the paid shared-center normalization, its computational parity
image survives in T_bad. Exact identity modes and dirty-independent
target restrictions still matter: (3) alone does not establish a
remaining relevant variable. The following actual source makes the
nonzero surviving direction concrete.

## 2. An actual four-wire one-clean source

Use physical wires u,v,r,s, with u the initialized input. Put
$`P_A=X_uX_v`$ and choose the following independent native lists:

```math
\begin{aligned}
E_{A,a}&=T_{Z_u}T_{Z_uZ_r}, &
E_{A,s}&=T_{X_v}T_{X_vZ_rZ_s}^\dagger,\\
E_B&=T_{Z_uZ_rZ_s}T_{Z_uZ_s}
      T_{Z_vZ_rZ_s}T_{Z_vZ_r},\\
E_{C,a}&=T_{Z_v}T_{Z_vZ_s}, &
E_{C,s}&=T_{X_u}T_{X_uZ_rZ_s}^\dagger,\\
E_A&=E_{A,a}E_{A,s}, & E_C&=E_{C,a}E_{C,s}.
\end{aligned}                                                    \qquad(4)
```

Each full endpoint list has four independent commuting labels, as
does B. Define actual Cliffords

```math
\begin{aligned}
C_0&=Z_uZ_rS_s\,CZ_{u,r}CZ_{u,s}CZ_{v,r}CZ_{r,s},\\
K_A&=CX_{u,v}H_u,\qquad K_C=C_0K_A,\\
V&=K_C^\dagger E_C^\dagger E_BE_AK_A.
\end{aligned}                                                    \qquad(5)
```

The diagonal parity-phase identity is literal on all sixteen columns:

```math
F_0:=E_{C,a}^\dagger E_BE_{A,a}
 =C_0(-1)^{(u+v)rs},\qquad
F_0P_AF_0^\dagger=C_0P_AC_0^\dagger=-P_AZ_s=:P_C.              \qquad(6)
```

The sum u+v in the Boolean exponent is modulo two. For example, the
Clifford exponent of C_0 modulo eight is
$`4u+4r+2s+4ur+4us+4vr+4rs`$, which fixes its scalar convention.
Conjugation by C_0 sends the two C-silent axes to the negatives of
one another. The two scalar factors from
$`T_{-P}=\omega T_P^\dagger`$ cancel, giving
$`C_0E_{C,s}C_0^\dagger=E_{C,s}`$ exactly.

Both endpoint clean-active ranks are two. The two clean images
$`J=E_AP_AE_A^\dagger`$ and $`S=E_CP_CE_C^\dagger`$ are
non-Pauli, each with four nonzero Pauli coefficients. Every B axis
is a negative singleton for J, so the middle pair count is
$`m_{\mathrm{pair}}=0`$. Equations (5)--(6) give the actual clean return.

On the initialized sector, let t=v denote the remaining core and put

```math
E=T_{X_t}T_{X_tZ_rZ_s}^\dagger,
\qquad F=Z_t^{rs}=(-1)^{trs}.
```

The occupied word is $`E^\dagger FE=F`$. In fact, E is I when
$`r+s=0`$ and $`\omega^{-1}S_{X_t}`$ when $`r+s=1`$; F differs
from I only at r=s=1. This proves the occupied CCZ target without
selecting or initializing either center bit.

## 3. Full columns, no shared Pauli flip, and native accounting

Writing c for the computational input bit on u, the full operator is

```math
V|c,t,r,s\rangle
 =(-i)^{c(r+s)}(-1)^{trs}
   |c,t+c(r+s),r,s\rangle.                                  \qquad(7)
```

Every addition inside a bit is modulo two. Equivalently
$`V=\mathrm{diag}(G,GL)`$ in the clean bit, where
$`G=CCZ_{t,r,s}`$ and L is I on r+s=0 and $`-iX_t`$ on r+s=1.
The inverse has the same permutation and the conjugate phases. This
is a fixed full-space extension, not a choice of unused action after
the occupied calculation.

There is no shared Pauli clean flip. Any clean-odd Pauli has an
off-diagonal logical Pauli block P. One block of its conjugate by V
contains $`G^\dagger PG L`$. The first factor has P's constant XOR
permutation; L adds the nonconstant translation $`t\mapsto t+r+s`$.
It cannot be a Pauli permutation. This excludes every Pauli with a
Pauli clean-odd pullback through V, and therefore excludes the shared
flip required by the two original encoders as well.

For an explicit native synthesis, let a change word W be applied
chronologically before its listed simultaneous phases and followed
by its reversed inverse. It implements $`W^\dagger T W`$.

| Band | Chronological change W | Simultaneous phases |
| --- | --- | --- |
| A | H v; CX v,s; CX r,s; CX u,r | T u; T r; T v; T-dagger s |
| B | CX s,v; CX s,u; CX r,s; CX r,u; CX v,s; CX v,r; CX r,v; CX u,v | T u; T v; T r; T s |
| C | H u; CX u,r; CX s,r; CX v,s | T v; T s; T u; T-dagger r |

Use the chronological word H u; CX u,v for K_A. Expand each displayed
CZ in C_0 as H on its target, CX, H; expand Z as S S. K_C is K_A
followed by C_0. The actual word is K_A, A, B, C-inverse, K_C-inverse.
Its separately reversed word retains the actual inverse of every gate.

| Resource | Forward word | Actual inverse |
| --- | --- | --- |
| Physical wires / clean / extra helpers | 4 / 1 / 0 | 4 / 1 / 0 |
| Native bands and T-count | 3 and 12 | 3 and 12 |
| T / T-dagger | 8 / 4 | 4 / 8 |
| H / CX | 14 / 34 | 14 / 34 |
| S / S-dagger, before positive-S expansion | 0 / 5 | 5 / 0 |
| Elementary Cliffords using H,S,CX | 63 | 53 |

The last row expands each $`S^\dagger=S^3`$ individually, without
post-expansion cancellations, in the forward and independently reversed
schedules. No hidden helper or phase oracle occurs. Tensoring untouched dirty wires gives the
same literal operator tensored with identity, including arbitrary
superpositions and external references. The different unused sector
in (7) is retained rather than replaced by a zero-clean extension.

## 4. The surviving center and encoded algebra

In the occupied form above, both endpoint lists are
$`\{X_t,X_tZ_rZ_s\}`$, the middle B-positive band is empty, and
$`R=\langle Z_r,Z_s\rangle`$, $`H=Z_t`$, $`b=rs`$. Both axes
are individually bad for F, despite their full band commuting with F.
Directly,

```math
T_{\mathrm{bad}}=R_0=\langle Z_rZ_s\rangle.                       \qquad(8)
```

For the computational Z probe, $`J_t=EZ_tE^\dagger`$ satisfies
$`J_tZ_tJ_t=Z_tZ_rZ_s`$ and fixes R. Its displacement of h is the
label Z_rZ_s; the other computational Z probes add no displacement.
The center is already saturated. This is the actual varying
translation r+s discussed in
[Native center saturation limits](NATIVE_CENTER_SATURATION_LIMITS.md).

Separately, the derivative-direction probe X_r gives the Clifford
$`J_r=EX_rE^\dagger`$ with

```math
J_rZ_tJ_r=Y_tZ_rZ_s.                                         \qquad(9)
```

This X probe is used for the diagnostic below, not for the
computational-Z family defining R_0.

There is also no fixed Clifford decoder of this encoded algebra into
a rank-three Pauli algebra. The eight-dimensional abelian algebra
$`E^\dagger\mathrm{Alg}(Z_t,Z_r,Z_s)E`$ has twelve distinct Pauli
support labels. In particular,

```math
E^\dagger Z_tE
 =\tfrac12\bigl(Z_t+Z_tZ_rZ_s+Y_t-Y_tZ_rZ_s\bigr).
```

Multiplication by the four center Paulis gives eight core support
labels, in addition to the four center labels. A fixed Clifford
permutes signed Pauli labels and preserves this support cardinality;
a rank-three Pauli algebra has just eight. This is a restriction on
that decoder, not an obstruction to other circuits or paid non-Clifford
decoding.

Finally, even replacing individual safety by the condition that the
whole endpoint image passes through F as a Clifford does not cover
the dual X_r direction here. Equation (9) gives

```math
(FJ_rF)Z_t(FJ_rF)=Y_tZ_rZ_s\,CZ_{r,s},                       \qquad(10)
```

which is not a Pauli. Hence FJ_rF is not Clifford. The **completed**
target derivative is nevertheless exactly $`CZ_{t,s}`$. Failure of
this intermediate condition must not be read as failure of all
two-layer derivative implementations. The next inference would need
actual cancellation involving both endpoint bands, or another paid
global mechanism.

The [two exact fixtures](../../tests/test_actual_native_center_witness.py)
check every original physical column, actual reverse words, the
native lists and counts, all clean-odd Pauli permutations, and the
encoded-algebra support and diagnostic (10). The witness has polar
rank two and is compatible with the retained bounded-polar-rank
theorem. It changes neither the complete-frame frontier nor the open
unrestricted one-clean three-layer question.
