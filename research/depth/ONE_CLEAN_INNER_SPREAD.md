# One-clean inner-spread rigidity

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Midpoint rigidity](ONE_CLEAN_MIDPOINT_RIGIDITY.md) · [Operator rigidity](FOUR_LAYER_OPERATOR_RIGIDITY.md) · [Crossing rank](FOUR_LAYER_CROSSING_RANK.md) · [Research index](../README.md)

This result, R77, gives a width-independent one-clean obstruction for
non-Clifford midpoints when either two-layer encoder has bounded inner
clean spread. This spread counts the inner T-layer axes anticommuting
with the encoded input clean Pauli. It is the binary symplectic
displacement rank of the corresponding Clifford core, and can grow
with physical width.

For the complete-frame witness

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),
\qquad \theta=\tfrac12\arccos(1/3),
```

a native four-layer realization admitting a split of inner clean spread
at most t has initialized-isometry error at least

```math
\eta_t=2^{-\mathcal M_t/2-2},
```

where

```math
K_t=2^{t+3},\qquad a_t=2(K_t^2-1),
```

```math
c_t=7K_t+5t+17,
\qquad
\mathcal M_t=c_t\frac{a_t^{K_t}-1}{a_t-1}.
```

The constants are conservative. The bound allows arbitrary dirty inputs,
unrestricted middle crossing rank, and a freely optimized common phase.
It does not assume a Clifford midpoint. An arbitrarily accurate
four-layer sequence would need unbounded inner clean spread at both
ends, as well as the growing middle crossing rank required by R74.
R76's Pauli and Clifford midpoint floors remain independently applicable,
without bounding this spread.
The unrestricted one-clean case and optimal depth remain open. The
compiler and its count/depth upper bounds are unchanged.

## 1. Inner clean spread

Let a native two-layer encoder be

```math
U=C_2D_2C_1D_1C_0.
```

Each $`D_j`$ is a physical T/T-dagger layer and each $`C_j`$ is Clifford. For
the physical clean Pauli $`Z_C`$, put $`P_0=C_0Z_CC_0^\dagger`$.
Among the active axes of $`D_1`$, let $`a_1,\ldots,a_t`$ be those
anticommuting with $`P_0`$. These axes are independent and isotropic.

Write the resulting clean image as

```math
S=UZ_CU^\dagger=EJE^\dagger,
\qquad E=C_2D_2C_2^\dagger,
```

```math
J=C_2C_1D_1P_0D_1^\dagger C_1^\dagger C_2^\dagger.
```

The core J is a Hermitian Clifford. Expanding the t nontrivial T
conjugations gives exactly $`2^t`$ distinct Hermitian Paulis:

```math
J=2^{-t/2}\sum_{j=1}^{2^t}\sigma_jP_j,
\qquad \sigma_j\in\{1,-1\}.
```

The silent axes cancel. Independence of the essential axes makes the
Pauli labels distinct; intervening Clifford conjugations retain the
flat coefficients.

The integer t is intrinsic to the specified Clifford core. In the
inner-layer coordinates its symplectic action satisfies

```math
F_J(v)=v+\sum_i\langle a_i,v\rangle a_i.
```

If M has the independent columns $`a_i`$ and $`\Omega`$ is the
symplectic form, then $`F_J-I=M M^{\mathsf T}\Omega`$.
The first map M is injective and $`M^{\mathsf T}\Omega`$ is
surjective onto $`\mathbb F_2^t`$. Therefore

```math
t=\mathrm{rank}_{\mathbb F_2}(F_J-I).
```

Clifford coordinate changes and T versus T-dagger signs do not change
this rank. Its value can depend on the native representation and
split; it is not asserted to be an invariant of S alone. These are
axes acting on the encoded clean observable, rather than necessarily
gates on the physical clean wire.

## 2. Compressed-image theorem

Let A,B,S be Hermitian involutions, with S traceless and
$`[A,S]=[B,S]=0`$. Suppose native two-layer encoders give the
correlated representations

```math
A=E_AK_AE_A^\dagger,
\qquad S=E_AJ_AE_A^\dagger,
```

```math
B=E_BK_BE_B^\dagger,
\qquad S=E_BJ_BE_B^\dagger.
```

All four cores are Hermitian Cliffords, and $`E_A`$,$`E_B`$ are Clifford
conjugates of physical T layers. Let $`Q=(I+S)/2`$ and let t be
the inner clean spread of the A encoder. With the constants above,

```math
\kappa_t=2^{-\mathcal M_t/2-1},
```

```math
\|(A-B)Q\|\lt\kappa_t
\quad\Longrightarrow\quad (A-B)Q=0.
```

The result is symmetric in A and B, so the smaller encoder spread
may be used. The proof also applies to phase-polynomial outer layers
when the stated flat Pauli expansion of the inner core is supplied.

### A finite Clifford expansion of the midpoint

Let $`L_A`$ be the isotropic axis space of $`E_A`$ and take any Hermitian
Pauli $`R\in L_A^\perp`$. Put

```math
G=ARA,\qquad H=BRB.
```

Since R commutes with $`E_A`$, the Pauli $`P_R=K_ARK_A`$ gives

```math
G=E_AP_RE_A^\dagger.
```

Thus G is a Hermitian Clifford. The same is true of SGS:

```math
SGS=E_A(J_AP_RJ_A)E_A^\dagger,
```

because $`J_AP_RJ_A`$ is a Pauli. In particular, SGS does not need
to be expanded quadratically in the number of terms of S.

The flat expansion from Section 1 gives

```math
S=2^{-t/2}\sum_{j=1}^{2^t}\sigma_jC_j,
\qquad C_j=E_AP_jE_A^\dagger.
```

Every $`C_j`$ is a Hermitian Clifford. Independently,

```math
H=E_BK_RE_B^\dagger,
```

```math
K_R=K_B(E_B^\dagger R E_B)K_B
```

is a Hermitian Clifford core. Consequently the compressed discrepancy
is

```math
\begin{aligned}
Q(G-H)Q={}&\tfrac14(G+SGS)\\
&+2^{-(t+4)/2}\sum_j\sigma_j(C_jG+GC_j)\\
&-\tfrac14E_B\mathcal B_RE_B^\dagger,
\end{aligned}
```

where

```math
\mathcal B_R=
K_R+J_BK_R+K_RJ_B+J_BK_RJ_B.
```

All left boundaries are bare Cliffords. The four right boundaries
share the single outer phase $`E_B`$. There are six terms with coefficient
magnitude $`1/4`$ and $`2^{t+1}`$ terms with coefficient magnitude
$`2^{-(t+4)/2}`$. Non-Hermitian products are allowed: they are
fixed ring-valued Clifford representatives with eighth-root common
phases.

### Equal coefficients and the finite-boundary gap

Use the common coefficient magnitude $`2^{-(t+4)/2}`$. Each of
the six larger coefficients has multiplier $`2^{t/2}`$. For even
t this is an integer and can be represented by repetitions. For odd
t use

```math
\sqrt2=\zeta_8+\zeta_8^{-1},
\qquad \zeta_8=e^{i\pi/4},
```

absorbing these eighth roots into the Clifford representatives.
Each larger term then contributes $`2^{\lceil t/2\rceil}`$
terms, for the total

```math
N_t=6\,2^{\lceil t/2\rceil}+2^{t+1}
\le K_t=2^{t+3}.
```

Both counts are even. Pad to $`K_t`$ with cancelling pairs if needed.
After Clifford-diagonalizing $`E_B`$, the normalized Choi discrepancy is
exactly the [finite-boundary form of R74](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
with

```math
K=K_t,\qquad b=t+4,\qquad m=0.
```

No crossing rank remains in this expression. The finite-boundary
proof fixes each summand's weighted-cubic phase separately; one common
phase polynomial is not required. Its recurrence has

```math
a_t=2(K_t^2-1),
```

```math
c_t=7K_t+(t+4)+4\log_2K_t+1
    =7K_t+5t+17.
```

Thus a nonzero discrepancy satisfies

```math
\|Q(G-H)Q\|_{2,\mathrm{tr}}^2
\ge2^{-\mathcal M_t}.
```

The norm uses normalized trace on the full physical space. For $`t=0`$,
the parameters are $`K_0=8`$, $`a_0=126`$, $`c_0=73`$, and
the exponent is exactly the eight-boundary exponent of R76. This
recovers its constant on the $`t=0`$ subfamily: a Clifford midpoint can
also have positive inner spread in a displayed representation.

## 3. Freezing and occupied phases

Put $`\epsilon=\|(A-B)Q\|`$. The restrictions of A and B to
$`Q\mathcal H`$ are unitaries, and $`\|QRQ\|\le1`$. Therefore

```math
\|Q(ARA-BRB)Q\|\le2\epsilon.
```

If $`\epsilon\lt\kappa_t`$, the finite-boundary gap forces
this discrepancy to vanish for every $`R\in L_A^\perp`$.
There is no union bound over the tested Paulis.

The [compressed commutant lift of R76](ONE_CLEAN_MIDPOINT_RIGIDITY.md#lifting-the-compressed-commutant)
does not require S to be Clifford. Applied to $`Y=BA|_{Q\mathcal H}`$,
it gives a unitary $`D_0`$ in the abelian Pauli algebra on $`L_A`$ with

```math
BAQ=D_0Q=QD_0.
```

In Clifford coordinates diagonalizing $`L_A`$, $`D_0`$ acts by unit scalars
$`\lambda_z`$ on its occupied blocks. The same coordinates make
$`E_A`$ diagonal with eighth-root entries up to an irrelevant common
phase. Every nonzero entry of $`S=E_AJ_AE_A^\dagger`$ is a
flat Clifford entry times an eighth root. Hence it is a unit of
$`\mathcal R=\mathbb Z[i,1/\sqrt2]`$.

Every nonzero row of Q contains a unit of this ring. If the row of S
has a nonzero off-diagonal entry, use $`Q_{ij}=S_{ij}/2`$.
Otherwise unitarity and Hermiticity give $`S_{ii}=\pm1`$,
and an occupied row has $`Q_{ii}=1`$. Since BAQ has entries in
the same ring, division in such an entry puts every occupied
$`\lambda_z`$ in $`\mathcal R`$. R73's unit-circle lemma
then makes it an eighth root.

A nonidentity occupied scalar would imply

```math
\|BAQ-Q\|\ge\sqrt{2-\sqrt2},
```

contradicting $`\epsilon\lt\kappa_t`$. Thus $`BAQ=Q`$, proving
the compressed-image theorem. The freezing step is where the bound
on inner clean spread is used.

## 4. Precision and the four layer sizes

Let J append the clean zero and suppose

```math
\|VJ-e^{i\phi}JU\|\le\eta
```

for the complete-frame witness tensored with dirty identity. Below
$`\tau=\sin(\pi/8)/(2\sqrt2)`$, R75 forces exact terminal
clean return. Every two-plus-two split $`V=V_2V_1`$ then has

```math
S=V_1Z_CV_1^\dagger=V_2^\dagger Z_CV_2.
```

For a target-preserved data Pauli P, take

```math
A=V_1PV_1^\dagger,
\qquad B=V_2^\dagger PV_2.
```

The initialized-isometry estimate gives $`\|(A-B)Q\|\le2\eta`$.
If either encoder has inner spread t, then
$`\eta\lt\eta_t=\kappa_t/2`$ freezes the required occupied
Paulis. As in R76, spectator and dirty X,Z and the active Y force
the occupied block to be a scalar times a discrete rotation. Its
distance from the witness is at least

```math
d_{\mathrm{rot}}
=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}
\gt\tau\gt\eta_t.
```

This contradiction proves the stated error floor. No assertion about
the unused clean sector's logical action is needed.

Define $`t_{\min}(V)`$ as the minimum of $`\min(t_A,t_B)`$
over native representations of V with at most four layers and their
compatible two-plus-two splits. Empty layers may be inserted. This
minimum exists because the possible values form a nonempty subset
of the nonnegative integers; no efficient method for finding it is
asserted.

An approximation of error at most $`2^{-L}`$, with integer
$`L\ge3`$, necessarily satisfies

```math
\boxed{\quad
\mathcal M_{t_{\min}(V)}\ge2L-4.
\quad}
```

The same inequality holds separately for both encoder spreads in
every representation. To express its scaling, the elementary bounds
$`c_t\le10K_t`$ and $`a_t\lt2K_t^2`$ give

```math
\log_2\mathcal M_t
\le16(t+4)2^t.
```

Consequently

```math
\log_2(2L-4)
\le16(t_{\min}+4)2^{t_{\min}}.
```

In particular,

```math
t_{\min}
\ge\log_2\log_2 L
-O(\log_2\log_2\log_2 L).
```

In a displayed four-layer native circuit, $`t_A`$ counts clean-active
axes in the first layer. The encoder $`V_2^\dagger`$ reads the
second half backwards, so $`t_B`$ counts them in the last layer. Together
with R74 this gives the following necessary growth in a hypothetical
arbitrarily accurate sequence:

| Native layer | Necessary number of active T/T-dagger gates |
|---|---|
| First | $`\Omega(\log\log L)`$ clean-active axes |
| Second | $`\Omega(\sqrt L)`$ from middle crossing rank |
| Third | $`\Omega(\sqrt L)`$ from middle crossing rank |
| Fourth | $`\Omega(\log\log L)`$ clean-active axes |

These are necessary resource conditions, not a construction and not
an unrestricted five-layer lower bound.

The theorem includes genuinely non-Clifford midpoints. For example,
on one qubit $`U=THT`$ sends X to

```math
UXU^\dagger=\frac Z{\sqrt2}+\frac{X-Y}{2}.
```

This is not a Hermitian Clifford, and the inner layer has $`t=1`$.
An input Clifford can relabel the clean Z as X; identity spectators
retain the same t at arbitrary width.

## 5. A common compatible flip need not exist

R76 gives a sufficient certificate using a Hermitian unitary C that
exchanges the S sectors and has simultaneous Pauli-sign actions on
A and B. The following exact native example shows that such a common
compatible flip does not always exist. It does not refute the open
occupied-commutator gap.

Let a be the clean qubit and b an arbitrary-input data qubit. Pauli
strings are ordered a,b. Set

```math
R_Q=e^{-i\pi Q/8}.
```

```math
C_A=(SH)_a\otimes H_b,
\qquad C_B=(SH)_a\otimes I_b,
```

where the single-qubit S in $`C_A`$,$`C_B`$ is the phase Clifford. Define

```math
U_A=(R_{IY}R_{ZY}^\dagger)
    (R_{YY}R_{ZX})C_A,
```

```math
U_B=(R_{IZ}^\dagger R_{ZZ}^\dagger)
    R_{ZX}C_B.
```

Operators act right to left. Each displayed pair has commuting,
independent axes. Thus $`U_A`$ has two native T layers and four T-type
gates; $`U_B`$ has two native layers and three T-type gates. Physical
routing is given below.

For the input data Pauli $`P=IZ`$, exact multiplication gives the
shared clean image

```math
\begin{aligned}
S&=U_A(ZI)U_A^\dagger\\
 &=U_B(ZI)U_B^\dagger\\
 &=\frac{-XX+XY+XZ+YI}{2}.
\end{aligned}
```

and data images

```math
A=\frac{IX-IZ+ZX+ZZ}{2\sqrt2}
  -\frac{YX+YZ}{2},
```

```math
B=\frac{-IX+IY-ZX-ZY}{2\sqrt2}
  +\frac{IZ}{\sqrt2}.
```

These are Hermitian involutions with $`[A,S]=[B,S]=0`$ and
$`\mathrm{Tr}(S)=0`$. The midpoint is not Clifford, since

```math
S(IX)S=\frac{-IY-IZ-ZY+ZZ}{2}
```

is not a Pauli. It does have a sector flip, namely $`\{ZI,S\}=0`$.
The obstruction concerns compatibility with both data images.

The actual circuit $`V=U_B^\dagger U_A`$ has at most four
physical T layers and seven T-type gates. Equality of the midpoint
images gives $`[V,ZI]=0`$, so it returns the initialized clean
qubit exactly for every data input.

### Trace obstruction to a common flip

Put $`Q_\pm=(I\pm S)/2`$. Exact products give

```math
\frac{\mathrm{Tr}(AB)}4=-\frac12,
\qquad
\frac{\mathrm{Tr}(SAB)}4=-\frac1{2\sqrt2}.
```

Both sectors have dimension two. Their normalized product traces are

```math
t_+=\frac{\mathrm{Tr}(Q_+AB)}2
    =-\frac12-\frac1{2\sqrt2},
```

```math
t_-=\frac{\mathrm{Tr}(Q_-AB)}2
    =-\frac12+\frac1{2\sqrt2}.
```

Suppose a Hermitian unitary C obeyed

```math
\{C,S\}=0,\qquad CAC=\alpha A,
\qquad CBC=\beta B,
```

with $`\alpha,\beta\in\{1,-1\}`$. Conjugation by C exchanges
the two sectors and sends AB to $`\alpha\beta AB`$. Trace
invariance would force $`t_-=\alpha\beta t_+`$, contrary to
their different absolute values. Therefore no such C exists, even
if it is allowed to be an arbitrary unitary. The same argument
excludes a non-Hermitian sector-exchanging unitary with the corresponding
conjugation sign actions.

### Occupied and full-space norms

Each A and B restricts to a traceless Hermitian involution on either
two-dimensional sector. The sector traces therefore give

```math
\|(A-B)Q_\pm\|
=\sqrt{3\pm1/\sqrt2}.
```

For a variant with a smaller occupied defect, replace $`U_B`$ by
$`U_B'=U_B(IX)`$. This only changes its initial Clifford and
gives $`B'=-B`$, with the same S and native layer counts. The
exact-return circuit $`V'=U_B'^\dagger U_A`$ satisfies

```math
\begin{aligned}
\|[V',\Pi P]\|
&=\|(A-B')Q_+\|\\
&=\sqrt{1-1/\sqrt2},
\end{aligned}
```

where $`\Pi=(I+ZI)/2`$. In contrast,

```math
\begin{aligned}
\|(A-B')Q_-\|
&=\|A-B'\|\\
&=\sqrt{1+1/\sqrt2}.
\end{aligned}
```

```math
\|A-SB'\|=\sqrt{3-1/\sqrt2}.
```

Neither multiplication choice $`S^kB'`$, $`k\in\{0,1\}`$,
has full-space difference norm equal to the occupied difference
norm. Thus an unconditional version of the common-flip norm-transport
identity is false. The occupied defect here is a fixed positive
constant, approximately 0.541196; spectator padding does not decrease
it.

### Physical-layer routing

The following Clifford routings specify literal physical T layers
on the same two wires. Gate lists are in application order.

| Routing | Application-order Clifford gates | Image of $`Z_a`$ | Image of $`Z_b`$ |
|---|---|---|---|
| $`B_{A1}`$ | H(a), CX(a,b), S(a), S-dagger(b), H(b) | YY | ZX |
| $`B_{A2}`$ | CX(b,a), H(b), S(b) | ZY | IY |
| $`B_{B1}`$ | CX(a,b), H(b) | ZI | ZX |
| $`B_{B2}`$ | CX(b,a) | ZZ | IZ |

Up to irrelevant common phases, the four layers are

```math
R_{YY}R_{ZX}=B_{A1}(T_aT_b)B_{A1}^\dagger,
```

```math
R_{IY}R_{ZY}^\dagger
=B_{A2}(T_a^\dagger T_b)B_{A2}^\dagger,
```

```math
R_{ZX}=B_{B1}T_bB_{B1}^\dagger,
```

```math
R_{IZ}^\dagger R_{ZZ}^\dagger
=B_{B2}(T_a^\dagger T_b^\dagger)B_{B2}^\dagger.
```

For each conjugated phase layer, apply the inverse routing list,
the indicated simultaneous T/T-dagger gates, then the routing list.
The input Clifford $`C_A`$ has list H(a), S(a), H(b), and $`C_B`$ has
H(a), S(a). The signed variant prepends X(b) to the $`C_B`$ list. No
additional clean or dirty qubits are used.

## 6. Verification and scope

The [exact algebra checks](../../tests/test_one_clean_inner_spread.py)
verify the finite-boundary bookkeeping, native core spread, the
non-Clifford midpoint examples, the shared midpoint, and the trace
and norm certificates for the common-flip counterexample. These
checks do not infer an asymptotic theorem from finite examples.

The proof uses R74's existing equal-coefficient finite-boundary gap,
R75's terminal clean-return rigidity, and R76's compressed commutant
lift and entry-ring argument. It uses no scalar polynomial-limit
premise. Unbounded inner spread remains allowed by the unrestricted
four-layer model, so the general one-clean five-layer obstruction
remains open.
