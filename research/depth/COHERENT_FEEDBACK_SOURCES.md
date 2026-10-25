# Local HRT operations and complete coherent feedback

[Residual-probe interface](DIRTY_CHOI_OBSERVABLE_SOURCE.md) ·
[HRT permutation](FULL_BANK_CHARACTER_TRANSPORT.md) ·
[Shared-control native batches](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers) ·
[Exact fixtures](../../tests/test_coherent_feedback_sources.py)

Every individual native one-qubit gate in the HRT coordinate frame has a
direct native implementation. In particular, an encoded Hadamard costs at
most sixteen T layers and needs no helper qubits. It mixes the physical key
after a nonlinear stage without implementing the full HRT permutation.

A complete retained-outcome feedback word built from this component cancels
its arbitrary borrowed register exactly. Its probe nevertheless loses the
required key response. A zero-key calibration leaves operator error two at
every size with at least two algebra generators. This is a failure of the
specified word, not a lower bound for general sources. No faster complete
indicator or improvement to the general frame frontier follows.

## 1. An exact local change of HRT coordinates

Use the square-zero algebra and complete permutation from the linked HRT
proof. There are $`w=2^r-1`$ arbitrary physical coefficient wires, labelled by
nonempty subsets of $`[r]`$. Write $`F|c\rangle=|f(c)\rangle`$ and
$`k=f^{-1}`$. For a nonempty S, put $`s=|S|`$ and define

```math
q_S(y)=\sum_{\substack{\{I,J\}: I,J\ne\varnothing\\I\sqcup J=S}}y_Iy_J,
\qquad k_S(y)=y_S+q_S(y).
```

Let $`B_S`$ add $`q_S(y)`$ into $`y_S`$, and let

```math
M_S=\prod_{\substack{T\ne\varnothing\\T\cap S=\varnothing}}
       \mathrm{CX}(y_T;y_{S\cup T}),
\qquad
C_S=\prod_{\substack{T\ne\varnothing\\T\cap S=\varnothing}}
       \mathrm{CCX}(y_S,y_T;y_{S\cup T}).
```

The pairs in $`M_S`$ are disjoint, and $`C_S`$ is its $`y_S`$-controlled
version. All coordinates in $`q_S`$ are strict subsets of S; neither their
values nor $`y_S`$ are changed by $`C_S`$. Both shears are involutions.
With $`D_{q_S}|y\rangle=(-1)^{q_S(y)}|y\rangle`$, the retained HRT identities give

```math
FX_SF^\dagger=X_SM_S,\qquad FZ_SF^\dagger=Z_SD_{q_S}.
```

The first equality follows by multiplying $`1+f(c)`$ by $`1+e_S`$;
the second reads the inverse coordinate $`k_S(y)`$. Now put
$`U_S=C_SB_S`$, in rightmost-first operator order. Direct conjugation gives

```math
U_SX_SU_S^\dagger=X_SM_S,\qquad
U_SZ_SU_S^\dagger=Z_SD_{q_S}.
```

The X and Z images determine the whole one-qubit matrix algebra, including
Y with its literal phase. Thus for every one-qubit operator A on coordinate S,

```math
FA_SF^\dagger=U_SA_SU_S^\dagger.
```

In particular, the encoded Hadamard has the exact factorization

```math
E_S:=FH_SF^\dagger=C_SB_SH_SB_SC_S.                 \qquad(1)
```

This palindrome can also be read chronologically. It is an equality on the
complete physical space; there is no initialized code or scalar gauge.

## 2. Sixteen native T layers without helpers

The numbers of quadratic terms and matching edges are

```math
a_s=2^{s-1}-1,\qquad b_s=2^{r-s}-1.
```

In $`B_S`$, each strict nonempty I is paired uniquely with its complement
$`S\setminus I`$. Hadamards on $`y_S`$ turn this shared-target Toffoli batch
into a shared-control CCZ batch on disjoint pairs. In $`C_S`$, Hadamards
on its $`b_s`$ targets give the same kind of CCZ batch. The linked exact
phase-polynomial schedule supplies four physical T bands per nonempty batch,
with no initialized or arbitrary helper.

For a batch of q pairs the schedule has $`6q+(q\bmod2)`$ T/T-dagger gates
and at most $`10q+3`$ Clifford gates before the Hadamard envelope. Reduce the
common root phase $`T^q`$ and its inverse modulo eight; its positive-S
remainder uses at most three S gates. This reduction is also made in the
actual reversed native word. Consequently valid common bounds for (1) are

```math
\begin{aligned}
T(E_S)&\le12(a_s+b_s)+2\bigl((a_s\bmod2)+(b_s\bmod2)\bigr),\\
D_T(E_S)&\le16,\\
G_{\rm Cliff}(E_S)&\le20a_s+24b_s+17,\\
G_{\rm all}(E_S)&\le T(E_S)+20a_s+24b_s+17,\\
h(E_S)&=a_{\rm clean}(E_S)=0.
\end{aligned}                                             \qquad(2)
```

Empty shears are omitted. At $`r=1`$ the source is simply H. Singleton and
top coordinates need only one shear type and at most eight T bands. For a
different native one-qubit gate A, replace the central H by its actual word
and add its resources, subtracting the one counted H where appropriate.
Every count remains $`O(2^r)`$ for a fixed-size A.

This supplies genuine later key mixing:

```math
E_SZ_SE_S^\dagger=X_SM_SD_{q_S}.                           \qquad(3)
```

The right side flips $`y_S`$ and is not diagonal. E is generally non-Clifford;
its four nonlinear shears must be charged. It is not a free Clifford wall
inside the three-block theorem. Distinct encoded Hadamards commute as
completed operators, but their native schedules share changing wires. No
parallel implementation of their product follows from that commutation.

## 3. A retained-outcome feedback word

Let $`\Omega=[r]`$, $`d=2^w`$, and let C,A,B be three arbitrary w-qubit
registers. The complete dirty Choi word from the linked proof has action

```math
\mathcal T|c,a,b\rangle
 =d^{-1}\sum_{t,y}(-1)^{a\cdot y+t\cdot b}
       |t,c+b+k(y),y\rangle.
                                                               \qquad(4)
```

Call its output registers t,u,y. Use the paid top encoded Hadamard
$`E_\Omega(t)`$ on the retained Fourier outcome t, then feed t into a
second quadratic query on u:

```math
J=Q_k(t;u)E_\Omega(t)\mathcal T.
```

This is a specified nonmonomial feedback word, with genuine key mixing
between nonlinear queries. Its complete probe action is

```math
J^\dagger Z_{u,\Omega}J
 =Z_{c,\Omega}\bigl(H_A D_{k_\Omega}(A)H_A\bigr)\otimes I_B.
                                                               \qquad(5)
```

To prove it, first use
$`Q_k^\dagger Z_{u,\Omega}Q_k=Z_{u,\Omega}D_{k_\Omega}(t)`$ and
$`E_\Omega D_{k_\Omega}E_\Omega=X_{t,\Omega}`$. In (4), the latter
translation contributes $`(-1)^{b_\Omega}`$. The u probe contributes
$`(-1)^{c_\Omega+b_\Omega+k_\Omega(y)}`$. The two old-b phases cancel,
leaving (5) after the actual A Fourier transform. Thus every B input and
reference correlation is included, rather than discarded as a Bell outcome.

The target is the residual phase

```math
\Xi_r(c,a)=(-1)^{\lambda(f(c+k(a)))}.
```

At $`c=0`$ it is exactly $`Z_{a,\Omega}`$. The A factor in (5) is
$`X_{a,\Omega}`$ times an operator on the other A coordinates, so it
anticommutes with that target. The uncalibrated derivative therefore has
operator error at least $`\sqrt2`$.

### Exact calibration does not recover the missing key response

Put $`L_A=H_AB_\Omega(A)`$. Since
$`D_{k_\Omega}=B_\Omega Z_\Omega B_\Omega`$, (5) gives

```math
(JL_A)^\dagger Z_{u,\Omega}(JL_A)
       =Z_{c,\Omega}Z_{a,\Omega}\otimes I_B.              \qquad(6)
```

This is exact on the full zero-key sector. For $`r\ge2`$, choose
$`c=e_i`$ and $`a=e_{\Omega\setminus\{i\}}`$. Then
$`f(c)=c`$, $`k(a)=a`$, and

```math
\lambda(f(c+k(a)))=\lambda(c+a+ca)=1,
\qquad c_\Omega+a_\Omega=0.
```

The calibrated derivative and target have opposite signs. Both are unitaries,
so the complete operator error is exactly two, for every B input. The same
witness persists as r grows. At $`r=1`$, (6) does equal the residual phase.

The source $`W=(JL_A)^\dagger`$ need not return B separately; its completed
derivative does return B and all private work. The failure is the lost
dependence on the non-top c coordinates. Reopening this word requires a
changed feedback or readout retaining that dependence, not a more accurate
implementation of its existing pieces.

### Complete native cost of the calibrated trial

Let $`t_k,g_k,d_k,h_k`$ be common forward/actual-inverse bounds for the paid
quadratic query on w coefficients, with $`g_k`$ counting every native gate.
Put $`a=2^{r-1}-1`$ and $`\pi_a=a\bmod2`$. Use the reduced Choi word,
which has one quadratic query and at most $`9w`$ additional Clifford gates.
The later query contributes another copy. The top mixer and calibration use
three $`B_\Omega`$ shears, one physical top H, and one full A-Hadamard layer.
Thus valid common bounds are

```math
\begin{aligned}
T(W)&\le2t_k+18a+3\pi_a,\\
D_T(W)&\le2d_k+12\le92,\\
G_{\rm all}(W)&\le2g_k+10w+48a+3\pi_a+16,\\
h_W&\le w+h_k,\qquad a_{\rm clean}=0.
\end{aligned}                                               \qquad(7)
```

Here $`h_W`$ excludes the principal C,A inputs and includes the arbitrary B
register. The private query pool can be reused after each complete query;
the local mixer and calibration need no private work. The word is exact,
and its native inverse is the actual reverse with the stated phase reductions.
It meets the relaxed count and width budgets but fails the operator target
by (6); pricing the circuit therefore does not create a qualifying source.

### Feeding both retained labels changes the work-return problem

A concrete repair adds u into t immediately before the last query:

```math
J_+=Q_k(t;u)\mathrm{CX}(u;t)E_\Omega(t)\mathcal T.
```

Its query consumes $`t+u`$. Write $`q=q_\Omega`$ and define the fixed
linear complement map on augmentation coefficients by

```math
(R(v))_S=v_{\Omega\setminus S}
\quad(\varnothing\ne S\ne\Omega),\qquad (R(v))_\Omega=0.
```

Quadratic polarization gives
$`q(t+u)=q(t)+q(u)+\lambda(tu)`$ and
$`\lambda(tu)=t\cdot R(u)`$. The output probe, pulled back through the
query and the encoded mixer, becomes exactly

```math
E_\Omega(t)\bigl[Z_{u,\Omega}D_{k_\Omega}(t+u)\bigr]E_\Omega(t)
       =X_{t,\Omega}(-1)^{q(u)+\lambda(tu)}.
```

The top mixer fixes the other t coordinates, so neither remaining phase
changes in this step. Undoing the t Fourier transform converts the bilinear
phase into $`b\mapsto b+R(u)`$. Undoing the Choi query and CNOT then shifts
c by the same vector. Explicitly, for $`v=c+b+k(a)`$, let

```math
\mathcal K|c,a,b\rangle
 =(-1)^{b_\Omega+q(v)}
       |c+R(v),a,b+R(v)\rangle.                            \qquad(8)
```

This is a Hermitian signed permutation: v and $`b_\Omega`$ are invariant
under its simultaneous shifts. The exact pullbacks are

```math
J_+^\dagger Z_{u,\Omega}J_+=H_A\mathcal K H_A,
\qquad
(J_+L_A)^\dagger Z_{u,\Omega}(J_+L_A)=\mathcal K.           \qquad(9)
```

For the second identity, the remaining two $`B_\Omega(A)`$ factors cancel:
they change only the top coordinate of $`k(a)`$, which neither q nor R uses.
Every original C,A,B input remains in (8); no outcome has been initialized
or discarded. The changed word restores dependence on the full key through
an explicit motion of the original dirty B register.

For every $`r\ge1`$, take $`c=b=0`$ and $`a=e_\Omega`$. Then
$`v=e_\Omega`$ and $`q(v)=R(v)=0`$, so (8) has sign +1, whereas
$`\Xi_r(0,e_\Omega)=-1`$. The complete calibrated error is again exactly
two, now already on the zero-key sector. This repair adds w CNOTs to (7)
and changes none of its T or width bounds.

The first word returns its dirty bank but loses the desired key response;
the repair introduces key-dependent bank motion without the target phase.
Jointly cancelling that motion while preserving the required residual phase
is unpaid. These two explicit failures do not exclude other feedback words.

## 4. Evidence and remaining source question

The finite fixtures compare complete encoded-H columns with the independent
HRT permutation, compile the literal four-band batches and actual inverses,
and check the feedback identity and calibrated witness. They support the
conventions and resource schedules; the displayed arguments establish the
all-size statements.

An individual local mixer is now paid. A joint source must still satisfy the
complete residual-probe identity, or a different sufficient accepted-Gram
contract, with sublogarithmic depth and the simultaneous count, width and
accuracy bounds. Serially applying all local mixers is not a faster source.
Neither the component nor this failed feedback word changes R96.
