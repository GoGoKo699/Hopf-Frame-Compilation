# Complete bank response and the limits of row isolation

[Current checkpoint](../../WORKSPACE.md) ·
[Coupled sources](COUPLED_QUARTER_PHASE_SOURCE.md) ·
[Matching stages](QUARTER_PHASE_MATCHING_LIMITS.md)

Composing two matching stages before taking a row-selective difference
does not preserve the single-stage cancellation of background rows.
This page derives the complete surviving conjugation and an exact
two-row Gram defect of three quarters. Echoing each stage first removes
that conjugation, but leaves a fixed-row component, not a complete lookup.
The aggregation statements below identify the missing native operation.
None is an unrestricted native-depth lower bound or an improvement to R96.

## 1. The response test and what two rows can establish

Put $`m=2^p`$, $`P_v=|v\rangle\langle v|\otimes I_{\rm work}`$,
$`D_Y=\sum_vY_vP_v`$, and $`H_Y=(I-D_Y)/2=(I+\Omega_Y)/4`$.
For a bank-preserving candidate, write $`G_Y=K(Y)^\dagger K(Y)`$.
The exact conditions

```math
G_0=I/2,\qquad
G_{Y\oplus e_v}-G_Y=-\frac{(-1)^{Y_v}}2P_v                 \qquad(1)
```

for every row and background telescope to $`G_Y=H_Y`$. Uniform fiber
identities include arbitrary coherent banks and external references.
The general source model may move the bank; (1) is a diagnostic for
bank-preserving candidates, not a restriction on that general model.

For a background S not containing u or v, the informative two-row test is

```math
G_{S\cup\{u,v\}}-G_{S\cup\{u\}}-G_{S\cup\{v\}}+G_S=0.       \qquad(2)
```

If (2) holds for every background, exact empty-bank and singleton values
imply (1): adding a different bit cannot change a row's increment.
Checking (2) only at the empty bank misses cubic and higher bank terms.
Commutation of the two difference operations, by contrast, is automatic
and supplies no test. A mixed-response witness of norm c forces target
error at least c/4 on one of its four banks, since the target has zero
mixed response.

A per-flip residual bounded by epsilon gives only
$`\|G_Y-H_Y\|\le |Y|\varepsilon`$ when the empty value is exact.
This is sharp even for positive contraction Grams: for
$`0\lt\varepsilon\le1/(2m)`$, the family

```math
G_Y=H_Y+\varepsilon|Y|I                                  \qquad(3)
```

has exact empty value, zero mixed response, and per-flip residual norm
epsilon, but final error $`|Y|\varepsilon`$. These abstract Grams make
no native implementation claim. Dimension-independent error control
needs additional structure, such as correctly localized residuals;
positivity and vanishing mixed response alone do not provide it.

## 2. Compose first: one background survives

Let M,N be disjoint CNOT matchings, with the
[paid dressed stages](QUARTER_PHASE_MATCHING_LIMITS.md#1-a-complete-matching-stage-with-constant-native-depth)

```math
R_M=\prod_u a_u^{s_uY_u},\qquad R_N=\prod_u b_u^{s_uY_u},
\qquad W=R_NR_M,\qquad F_v=X_{s_v}.
```

Factors commute within each stage, not between stages. Use the same
arbitrary selectors in both stages. With one clean signal, emit
chronologically $`H_a,CX(a;s_v),W,CX(a;s_v),Z_a,H_a`$, with Z
implemented by two S gates. The accepted block is

```math
K_v=(W-F_vWF_v)/2.                                      \qquad(4)
```

Write $`A=\prod_{u\ne v}a_u^{s_uY_u}`$ and
$`B=\prod_{u\ne v}b_u^{s_uY_u}`$. When $`Y_v=0`$, K vanishes.
When $`Y_v=1`$, the two selector-v fibers of W are
$`W_0=BA`$ and $`W_1=Bb_vAa_v`$. Their differences have opposite
signs, so the entire Gram is independent of s_v and equals

```math
G_v=\frac{Y_v}{2}
 \left[I-\mathrm{Re}(A^\dagger b_vAa_v)\right]
 \otimes I_{s_v}.                                       \qquad(5)
```

All second-stage background factors B cancel. The first-stage
background conjugates b_v and does not generally cancel. Formula (5)
is an operator identity on every other selector fiber, not a basis
probability assertion. It therefore includes their coherent inputs.

For stage ledgers $`(T_M,G_M,d_M,b_M)`$ and similarly for N, both the
source and its actual inverse have

```math
T\le T_M+T_N,\quad G\le G_M+G_N+6,\quad D_T\le248,
\qquad b_{\rm dirty}\le m+\max(b_M,b_N),\quad a_{\rm clean}=1.
                                                               \qquad(6)
```

There is zero implementation error. The additional stage helpers return
exactly and their pool is reused. The arbitrary selectors are included
in dirty width; accepting the signal does not promise their return.
For two single-CNOT stages, (6) gives
$`T\le196m+32`$, $`G\le222m+46`$ and $`b_{\rm dirty}\le4m`$.

### Incidence already destroys the marked kernel

Suppose M contains c-to-t and t is a control of N. Choose active rows
v and $`u=v\oplus e_c`$, with every other bank bit zero. At s_u=0,
$`A=I`$ and both a_v,b_v fix v, so (5) has zero marked expectation.
At s_u=1, $`A=a_u`$ and $`a_u|v\rangle=i|v\oplus e_t\rangle`$.
The latter basis vector is moved by b_v. Hence
$`\langle v|a_u^\dagger b_va_ua_v|v\rangle=0`$, and the marked
acceptance is one half.

If the two matching control sets cover every coordinate, they must be
disjoint perfect half-sets; every first-stage target is then a
second-stage control. Thus covering all coordinates with two such
stages cannot preserve this row echo's common marked kernel on every
background. This is a failure of (4), not of arbitrary coherent repairs.

## 3. Exact two-row defect, including coherent inputs

Take p=2, $`a=S_0CX(0;1)`$, $`b=S_1CX(1;0)`$, v=0, u=1,
and $`Y=e_0+e_1`$. Label addresses by $`x=x_0+2x_1`$.
For s_u=0,1, let $`Q_0=b_0a_0`$ and
$`Q_1=a_1^\dagger b_0a_1a_0`$. Their complete address actions are

| Input | Q_0 output | Q_1 output |
| --- | --- | --- |
| $`\lvert0\rangle`$ | $`\lvert0\rangle`$ | $`-\lvert3\rangle`$ |
| $`\lvert1\rangle`$ | $`-\lvert2\rangle`$ | $`i\lvert0\rangle`$ |
| $`\lvert2\rangle`$ | $`i\lvert3\rangle`$ | $`\lvert2\rangle`$ |
| $`\lvert3\rangle`$ | $`i\lvert1\rangle`$ | $`i\lvert1\rangle`$ |

Put $`G_j=(I-\mathrm{Re}Q_j)/2`$. Direct multiplication gives
$`G_j^2=3G_j/4`$ and $`\mathrm{tr}G_j=3/2`$: each is three
quarters of a rank-two projector. Their difference is

```math
D=G_1-G_0=\frac14
\begin{pmatrix}
2&-i&0&1\\ i&0&-1&0\\0&-1&-2&-i\\1&0&i&0
\end{pmatrix},\qquad
D^3=D/2,\quad\mathrm{tr}D=0,\quad\mathrm{tr}D^2=1.
                                                               \qquad(7)
```

Thus $`\|D\|=1/\sqrt2`$. Any selector-independent operator is at
distance at least $`1/(2\sqrt2)`$ from these two Gram fibers. Their
midpoint attains that minimax distance. At fixed s_u=1, (7) is also
the empty-background mixed response (2) for the partial row source,
which vanishes whenever $`Y_v=0`$.

For the actual two-row lookup target, the defect is larger and sharp:

```math
\psi=(|0\rangle+i|1\rangle)/\sqrt2,\qquad
H_Y\psi=0,\qquad G_1\psi=3\psi/4,
\qquad\|G_v-H_Y\|=3/4.                                  \qquad(8)
```

The upper bound follows from $`0\preceq G_j\preceq3I/4`$ and
$`0\preceq H_Y\preceq I/2`$. The witness is entirely inside the
marked subspace. Embedding these two coordinates and rows gives the
same defect for every p at least two, with spectator address bits
fixed to zero. This is an exact two-row-bank result, not an assertion
that all other banks have the same defect.

Scalar processing of this relative word cannot repair it. At the same
bank, Q_0 fixes the marked basis vector 0 and Q_1 fixes the unmarked
basis vector 2. Signal/auxiliary interlayers acting trivially on all
address, bank and dirty-selector work see the same eigenvalue one.
Their Gram error is therefore at least one quarter. Processing the
stages separately, changing the dirty work, or changing acceptance is
outside this spectral-only statement. A positive sum of accepted-channel
Grams also cannot cancel marked leakage: weight lambda on (8) contributes
at least $`3\lambda/4`$ there. Coherent recombination need not be positive
term by term and is not excluded.

## 4. Echo first: exact fixed-row isolation is possible

For each individual stage, the opposite order of operations works:

```math
\begin{aligned}
E_{M,v}&=R_M^\dagger F_vR_MF_v
       =a_v^{(1-2s_v)Y_v},\\
\left(\frac{R_M-F_vR_MF_v}{2}\right)^\dagger
\left(\frac{R_M-F_vR_MF_v}{2}\right)
 &=\frac{Y_v}{2}(I-\mathrm{Re}a_v)\otimes I_s.
\end{aligned}                                                   \qquad(9)
```

Every other row cancels before different generators are combined.
The unitary echo has bounds $`2T_M,2G_M+8,2d_M`$, reuses the
stage's private pool, and has zero clean inputs and exact error.
The two literal selector X gates each cost four Cliffords.

For the fourteen [dressed expander generators](QUARTER_PHASE_MATCHING_LIMITS.md#dressed-expander-gap)
on their admitted lengths $`p=3s`$, translate their Hamiltonian to v,
calling it H_v and retaining its proved positive constant gamma. Then

```math
\frac1{14}\sum_j K_{j,v}^\dagger K_{j,v}=Y_vH_v/2,
\qquad
\frac\gamma2Y_v(I-|v\rangle\langle v|)
\preceq Y_vH_v/2\preceq Y_v(I-|v\rangle\langle v|).       \qquad(10)
```

This statement is uniform over every bank and all work, but it remains
fixed-row. An average of Grams is not automatically the Gram of the
coherently averaged accepted blocks.

A square compression of the unitary echoes can be paid explicitly.
Four clean label bits select fourteen echoes and two identities.
For each j, implement the selected echo as
$`R_j^\dagger C_{\lambda=j}(F_v)R_jC_{\lambda=j}(F_v)`$.
Inactive labels cancel by the actual inverse; no controlled whole stage
is assumed. A four-control XOR is the fixed-size returned dirty ladder
in [the scalar construction](EXPANDER_ZERO_REFLECTION.md#5-literal-native-implementation-and-source-costs).
Use label Hadamards before and after SELECT. Its complete compression is

```math
A_v=\begin{cases}
I,&Y_v=0,\\
X_vAX_v,&Y_v=1,\ s_v=0,\\
X_vA^\dagger X_v,&Y_v=1,\ s_v=1,
\end{cases}\qquad A=(2I+\sum_jr_j)/16.                    \qquad(11)
```

Here X_v is the address translation, not a supplied register. With
$`P_v^{\rm fix}=(1-Y_v)I+Y_v|v\rangle\langle v|`$, (11) is identity
on that projector and contracts its complement by
$`c=\sqrt{1-7\gamma/32}\lt1`$. All other bank and selector bits
cancel. The complete SELECT and actual inverse have constant T-depth,
$`T,G,b_{\rm dirty}=O(pm)`$, four clean labels, and zero synthesis
error. This is a fixed-row compression, not the target lookup source.
The already-paid scalar zero test also supplies fixed-row predicates;
(11) changes background cancellation, not the general frontier.

## 5. Aggregation is an operation, not a normalization convention

Making the fixed label v coherent would require

```math
F_V=\sum_v|v\rangle\langle v|_V\otimes X_{s_v}.           \qquad(12)
```

This is exactly the complete dirty indicator on the selector bank.
It cannot be silently used to solve that same missing primitive.
Granting a uniform p-clean row label instead yields an average, for example

```math
\frac1m\sum_v P_v^{\rm fix}
 =\left(1-\frac{|Y|}{m}\right)I+\frac{D_Y}{m},\qquad
\frac1m\sum_vY_v|v\rangle\langle v|=D_Y/m.                \qquad(13)
```

The normalization and p initialized bits are real costs.

There is no universal clean-history obstruction to combining ideal
orthogonal row phases. In fact

```math
W_Y=\prod_v\exp(i\pi Y_v|v\rangle\langle v|/2)
   =I+(i-1)D_Y,\qquad
K=(I+iW_Y)/2,\quad K^\dagger K=H_Y,\quad W_Y^2=\Omega_Y.
                                                               \qquad(14)
```

The exact identity needs no growing failure history. Its joint native
implementation is unpaid and would already supply the complete phase
query by squaring. Orthogonality of the ideal row supports does not
localize approximation errors: multiplying approximate row unitaries
with only per-row norm error epsilon guarantees at most m epsilon.
Errors supported inside the respective orthogonal row blocks would give
a maximum-error bound instead; that stronger physical property needs proof.

### A global selector flip loses a whole affine cube

Replacing every F_v by $`F_{\rm all}=\prod_vX_{s_v}`$ is native but
does not implement (12). Let Y mark any three-dimensional affine cube
of eight row labels, with p at least three, and set all selectors to zero.
For every dressed disjoint CNOT matching M,

```math
E_{M,\rm all}=R_M^\dagger F_{\rm all}R_MF_{\rm all}
 =\prod_u r_{M,u}^{(1-2s_u)Y_u}=I
 \quad\hbox{on the selector-zero fiber}.                 \qquad(15)
```

For each control coordinate, the mismatch count across the cube is
0,4 or 8. Its contribution is $`(iX_t)^{\rm count}=I`$.
This holds at every address and for every choice of matching, including
changed control partitions. Each completed echo preserves selectors.

Consequently any word of these completed echoes and their inverses,
with arbitrary common address/clean-signal/private-work interlayers
that act trivially on the dirty selector and bank registers, has the
same accepted action at this cube bank as at the empty bank on that
fiber. Weight phases $`i^{|Y|}`$ also coincide. The two required Gram
expectations differ by one half at each cube address, giving worst-case
error at least one quarter. Dirty-selector mixing or other bank access
leaves this restricted grammar; there is no unrestricted depth conclusion.

## 6. Decision and finite evidence

Pause compose-then-row-echo acceptance, scalar processing of its relative
word, and the unmixed global-flip batching repair. Fixed-row isolation
and its uniform spectral statement remain valid components. A complete
repair must pay collective acceptance without assuming (12), losing the
normalization in (13), or restoring the old cost through m-fold precision.
The [co-moving stage](CO_MOVING_LOOKUP_STAGE.md) changes the mechanism
by moving the physical bank together with the address; its own accepted
operator, not (1) on fixed banks, must then be assessed.

The [exact fixtures](../../tests/test_bank_bit_response.py) check complete
small-bank operator identities, the coherent two-row witnesses, and the
restricted cube collision. A seven-wire native fixture fixes two rows
active and checks all 128 columns and the actual reverse, including both
arbitrary private helpers. It is not a scalable arbitrary-bank native
emitter. These finite checks support analytic identities, not asymptotic
inference or external mathematical review.
