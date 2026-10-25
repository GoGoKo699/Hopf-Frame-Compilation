# Bank mixing and the complete lookup phase

[Complete-indicator contract](../../WORKSPACE.md#main-unfinished-route-and-exact-missing-operation) ·
[Native relative-block source](COUPLED_FINAL_ENCODER.md#relative-probe-involution-repair) ·
[Source-specific route decisions](../ROUTE_DECISIONS.md#source-specific-representation-limits)

This chapter separates complete-action restrictions on two mixing grammars
from the exact phase/indicator reduction. Affine/Fourier dressings fail in
the repaired source word. Two Clifford mixers retain a rank-one exception,
but its required derivative has the same resource target as the complete
indicator. None of these statements is an unrestricted native-depth bound.

The bank B has $`m=2^p`$ bits indexed by binary p-bit strings a. Let H be the
binary syndrome matrix whose column a is a, let $`\pi`$ denote bank parity,
and put $`h_t=H^{\mathsf T}t`$. The co-moving address s and complete target are

```math
\begin{aligned}
s&=x+Hy,& \pi(y)&=\sum_a y_a,\\
q_s(y)&=\pi(y)+y_{s+Hy},\\
Q_s|y\rangle&=(-1)^{q_s(y)}|y\rangle.
\end{aligned}
```

All additions in these Boolean expressions are binary. Every bank, work
and reference input is arbitrary. The native ledger uses the
[co-moving Clifford W](COUPLED_FINAL_ENCODER.md#81-literal-native-words-and-their-simultaneous-ledger)
and [corrected middle](COUPLED_FINAL_ENCODER.md#relative-probe-involution-repair)
from the relative-block source. Analytic freedom in a phase or mixer does
not supply a native implementation.

- [Affine/Fourier transition gap and exact monomial screen](#affine-fourier-complete-transition-gap).
- [Two-Clifford mixing grammar and its rank-one exception](#two-clifford-mixer-screen).
- [Complete phase/indicator reductions and the full commutant freedom](#rank-one-indicator-equivalence).

<a id="affine-fourier-complete-transition-gap"></a>

## 1. Affine bank mixing around a full Fourier layer still misses lookup

A full bank Hadamard layer can be combined with an arbitrary invertible
linear change of bank coordinates, instead of changing only a Hadamard
mask. This leaves the diagonal-dressing premise of [Section 9.2](COUPLED_FINAL_ENCODER.md#relative-repair-diagonal-modulus) and makes
nonlocal bank transitions. Nevertheless the [corrected source (62)](COUPLED_FINAL_ENCODER.md#relative-probe-involution-repair) has a
complete error floor of square root two for both orders of those layers.
The obstruction persists under arbitrary X-basis phase functions in T
and arbitrary computational-diagonal middle phases. It concerns the
displayed source word; it is not a lower bound for arbitrary native
circuits or for two genuinely nonmonomial mixing stages.

### 1.1. Complete transition support for both affine/Fourier orders

Let $`p\ge2`$, $`m=2^p`$, and let L be any phase-free affine bank
permutation

```math
L|y\rangle=|Ay+b\rangle,\qquad A\in\mathrm{GL}(m,2).
```

Take either $`C=H_B L`$ or $`C=L H_B`$, where
$`H_B=H^{\otimes m}`$. The affine offset b is arbitrary; no initialization
of bank bits is allowed. For this subsection allow T to be any unitary
diagonal in the joint X basis, and $`D_s`$ to be any computational-diagonal
unitary. Retain the same source form
$`E_s=C^\dagger T^\dagger D_s^\dagger C T D_s`$ and
$`F_s=E_sZ_\pi E_s^\dagger`$. This analytic generality supplies no free
native implementation of T or $`D_s`$. Write $`\mathbf1=\boldsymbol1_m`$
and $`H_0=T Z_\pi T^\dagger`$. Conjugating actual and target by TC gives

```math
\mathcal A_s=D_s^\dagger C H_0 C^\dagger D_s,
\qquad
\mathcal R_s=T C Q_s C^\dagger T^\dagger.                 \qquad(71)
```

Because T commutes with every X-type Pauli while $`Z_\pi`$
anticommutes with the odd-weight ones,

```math
H_0 X^u H_0=(-1)^{\mathbf1\cdot u}X^u.
```

The affine L conjugates $`X^u`$ to $`X^{Au}`$ without a phase.
It follows that $`C H_0 C^\dagger`$ is computational-monomial,
with a single fixed displacement v in every column:

| C | Displacement v |
|---|---|
| $`H_B L`$ | $`A^{-\mathsf T}\mathbf1`$ |
| $`L H_B`$ | $`A\mathbf1`$ |

For the first row, conjugation of each $`Z^u`$ by $`C H_0 C^\dagger`$
returns $`(-1)^{(A^{-\mathsf T}\mathbf1)\cdot u}Z^u`$; the
joint computational eigenspaces are one dimensional. The second row
follows by first using $`H_B`$ and then applying L. The affine offset can
change the phase of a column but not its displacement. The diagonal
$`D_s`$ in (71) likewise changes only the column phase. Thus for every y,

```math
\mathcal A_s|y\rangle=\gamma_{s,y}|y+v\rangle,
\qquad |\gamma_{s,y}|=1.                                  \qquad(72)
```

The [Walsh identity (65)](COUPLED_FINAL_ENCODER.md#relative-repair-diagonal-modulus), which remains valid even when its masks are not
distinct at p=2, has support among

```math
\mathbf1+h_t+e_a,\qquad a,t\in\mathbb F_2^p.
```

None of these masks equals $`\mathbf1`$: $`h_0`$ is zero, and each
nonzero $`h_t`$ has weight $`m/2\ge2`$, whereas $`e_a`$ has weight one.
For $`C=H_B L`$, the corresponding X-displacement masks of
$`C Q_s C^\dagger`$ are $`A^{-\mathsf T}(\mathbf1+h_t+e_a)`$;
for $`C=L H_B`$ they are $`A(\mathbf1+h_t+e_a)`$. Affine
phases affect coefficients only. Neither support contains v.

Both conjugated targets are sums of X-type Paulis, so they commute with
the X-diagonal T. Therefore $`\mathcal R_s=C Q_s C^\dagger`$
and

```math
\langle y+v|\mathcal R_s|y\rangle=0
\quad\text{for every }s,y.
```

Equations (71)--(72) give normalized orthogonal output columns. Hence

```math
\|F_s(C)-Q_s\|\ge\sqrt2
\qquad\text{for every }s,\ p\ge2.                         \qquad(73)
```

Before the TC change of basis, an explicit witness is
$`C^\dagger T^\dagger|y\rangle`$. Embedding it in the minus-Y
co-moving probe fiber gives the same complete native Pauli-image error.
These witnesses test permitted arbitrary inputs; no state preparation
resource is supplied to the source.

### 1.2. The native costs are paid, but do not repair the transition

Suppose L is emitted using $`g_L`$ elementary Cliffords in a common
forward/actual-inverse bound, with CNOTs, swaps expanded into CNOTs,
and X expanded as HSSH. Both displayed choices of C have a physical
implementation $`C_{\rm phys}=W^\dagger C W`$ with

```math
G(C_{\rm phys}),\ G(C_{\rm phys}^\dagger)
\le pm+3m+g_L.
```

It uses no T gates or extra helpers. For the native product angles and
corrected middle of [Section 9.1](COUPLED_FINAL_ENCODER.md#relative-probe-involution-repair), the original X-axis source bands
therefore give the complete candidate ledger

```math
T=4m,\qquad D_T\le4,\qquad
G\le6pm+26m-12+2g_L,
\qquad b_{\rm extra}=1,\quad a_{\rm clean}=0.              \qquad(74)
```

The extra bit is the arbitrary physical probe. The actual Pauli echo has
T-count 8m, depth at most eight, and Clifford count at most twice the
bound in (74) plus two. Gaussian elimination supplies an emitted
$`g_L=O(m^2)`$ for any specified affine map, so these candidates fit the
permitted Clifford budget. The failure in (73) is their complete
operator action, not an unpaid control or uncounted native gate.

### 1.3. No monomial dressing gives exact lookup at the native angles

Return here to the native product $`T=\prod_j e^{-i\pi X_j/8}`$
and the corrected middle $`D_s`$ of [Section 9.1](COUPLED_FINAL_ENCODER.md#relative-probe-involution-repair). A separate exact statement
covers every computational-monomial C on the bank, including nonlinear basis permutations and arbitrary diagonal
phases, without granting it a free implementation. Such a C preserves
the fact that every entry of $`H_0`$ has modulus $`2^{-m/2}`$. Thus the
actual $`\mathcal A_s`$ in (71) has that same modulus everywhere.
The transformed target is
$`\mathcal R_s=T\mathrm{diag}(f)T^\dagger`$, where
$`\mathrm{diag}(f)=C Q_s C^\dagger`$ and f takes values
in $`\{+1,-1\}`$.

Let $`\widehat f(S)`$ be its normalized Boolean Walsh coefficients.
A direct product of the one-qubit transition probabilities of T gives

```math
\langle y|\mathcal R_s|y\rangle
 =\sum_{S\subseteq[m]}\widehat f(S)\,2^{-|S|/2}(-1)^{S\cdot y}.
```

Parseval consequently yields

```math
2^{-m}\sum_y|\langle y|\mathcal R_s|y\rangle|^2
 =\sum_S\widehat f(S)^2 2^{-|S|}\ge2^{-m}.                \qquad(75)
```

Here $`\sum_S\widehat f(S)^2=1`$, and equality in (75) is
possible only when all Fourier mass is at $`S=[m]`$. If
$`F_s(C)=Q_s`$ exactly, the flat actual matrix forces equality, hence

```math
C Q_s C^\dagger=\pm Z_\pi.                              \qquad(76)
```

Thus a monomial dressing that completed this exact four-band word would
already transport the full lookup sign to parity. The remaining flat
matrix identity is stronger: it forces the dressing's basis permutation
to be affine, so no monomial dressing completes the lookup for p at least
two.

To prove this, write $`S_B=\prod_j S_{B_j}`$. The native product has
the literal identity

```math
H_0=S_B^\dagger H_B S_B.
```

Let P be the basis permutation underlying C, with
$`P|y\rangle=|\sigma(y)\rangle`$. Substituting (76) into (71)
gives $`D_s^\dagger C H_0 C^\dagger D_s=\pm H_0`$. Absorbing
all diagonal factors on both sides, including the arbitrary phases of C,
therefore gives a diagonal unitary J such that

```math
J P H_B P^\dagger J^\dagger=\pm H_B.                    \qquad(77)
```

Take the multiplicative rectangle ratio of entries at
$`(x,y),(x,0),(0,y),(0,0)`$. Every entry is nonzero; the four
row/column phases of J and the common sign cancel. The Walsh kernel
$`(H_B)_{xy}=2^{-m/2}(-1)^{x\cdot y}`$ then gives, over the
binary field,

```math
g(x)\cdot g(y)=x\cdot y,
\qquad g(x)=\sigma^{-1}(x)+\sigma^{-1}(0).              \qquad(78)
```

The map g is bijective and fixes zero. For all x,x',y, (78) implies
$`[g(x+x')+g(x)+g(x')]\cdot g(y)=0`$. Surjectivity in y makes
the bracket zero. Hence g is linear and orthogonal, and sigma is affine.
But (76) would then make the exponent of $`Q_s`$ affine as well. This
contradicts its degree p+1 from [Section 10.4](COUPLED_FINAL_ENCODER.md#completed-coordinate-round-boundary) whenever p is at least two.
Consequently **no computational-monomial C**, even with a nonlinear
permutation and arbitrary phases, supplies exact complete lookup in this
source word at the fixed native T angles. No synthesis-cost assumption
on C was used.

This exact screen does **not** give a dimension-independent approximation
gap for arbitrary monomial dressings. The diagonal restriction only gives

```math
\|F_s(C)-Q_s\|\ge
\left(\sum_S\widehat f(S)^2 2^{-|S|}\right)^{1/2}-2^{-m/2}.
```

This bound can be small; it must not be promoted to a fixed-accuracy exclusion.
The robust result (73) has the separate affine/Fourier hypotheses and
already excludes angle tuning or a changed diagonal middle in that word.
Neither result excludes a genuinely different nonmonomial interlayer,
a nondiagonal changed middle or an accepted-Gram source. R96 remains
unchanged.

The [affine/Fourier fixtures](../../tests/test_affine_fourier_repair.py)
check both layer orders with a nontrivial CNOT basis and affine offset,
fixed transition displacement, exact target support exclusion, actual
p=3 source inverses and native echoes on the coherent witnesses, and
the native ledger through p=6. A separate changed-word fixture uses
correlated X-basis cubic phases and a new diagonal middle. The monomial
fixture includes cubic diagonal phases and checks the exact noise
identity and flat actual columns. A nonlinear permutation that does
map the lookup sign to parity is rejected by the exact Walsh-rectangle
identity, testing the second step of the monomial proof. These finite checks challenge the
analytic statements; they are not a general synthesis or depth proof.

<a id="two-clifford-mixer-screen"></a>

## 2. Two Clifford mixing stages: a transition screen and its exact exception

This screen concerns a specified complete bank-source grammar. It does not
bound arbitrary native T-depth, the earlier four-band X-rotation family,
or sources using additional active work registers or intermediate address
motion. Fix $`p\ge3`$, $`m=2^p`$, and a co-moving address $`s`$. All $`m`$ bank
bits and any passive reference are arbitrary. Put

```math
P=Z_\pi,\qquad Q_s|y\rangle=(-1)^{q_s(y)}|y\rangle,
\qquad q_s(y)=\pi(y)+y_{s+Hy}.
```

Allow the actual bank word

```math
E=D_2 C_2 D_1 C_1 D_0,\qquad F=EPE^\dagger,             \qquad(T1)
```

where $`C_1,C_2`$ are arbitrary bank Cliffords and the $`D_i`$ are arbitrary
bank-diagonal unitaries. No native implementation of a general $`D_i`$ is
granted. The screen is therefore valid even after giving this diagonal
freedom for free. The rightmost operation acts first; the inverse is the
actual adjoint. Diagonal boundaries $`D_0,D_2`$ may have literal phases.

### 2.1. A Fourier-coset identity for the complete diagonal

Write $`C_1PC_1^\dagger`$ as a Hermitian Pauli with binary X label $`a`$. Since
$`D_1`$ is diagonal, $`D_1(C_1PC_1^\dagger)D_1^\dagger`$ has matrix entries only
between $`y`$ and $`y+a`$. Its Pauli expansion is contained in one affine Pauli
coset with fixed X label $`a`$. Let

```math
\mathcal A=C_2\langle Z_{B_0},\ldots,Z_{B_{m-1}}\rangle C_2^\dagger,
\qquad V=\{u:Z^u\text{ belongs to }\mathcal A
                         \text{ up to scalar}\}.
```

Define $`h=m-\dim V`$; equivalently $`h`$ is the rank of the binary X-label
projection of this maximal commuting Pauli group. After $`C_2`$, the Pauli
support is in one affine coset of $`\mathcal A`$. Its intersection with the diagonal
Paulis is either empty or an affine binary coset $`u_0+V`$. $`D_0`$ commutes with
$`P`$, while $`D_2`$ leaves every diagonal entry unchanged. Consequently the
complete diagonal $`f(y)=\langle y|F|y\rangle`$, which is real, has Fourier support in
$`u_0+V`$, or vanishes identically. In the nonempty case,

```math
f(y+v)=(-1)^{u_0\cdot v}f(y)
\qquad\text{for every }v\in V^\perp.                    \qquad(T2)
```

This is an identity of the actual source, not a proposed classical
approximation of its transitions. If $`\|F-Q_s\|\lt1`$, then
$`|f(y)-(-1)^{q_s(y)}|\lt1`$ on every basis column. Thus $`f`$ is nonzero and has
the target's sign everywhere. Equation (T2) forces

```math
q_s(y+v)+q_s(y)=u_0\cdot v
\qquad(v\in V^\perp),                                  \qquad(T3)
```

uniformly in $`y`$. The complete lookup sign has only one nonzero linear
structure: $`v=\boldsymbol1_m`$, with derivative one. This follows from its distinct
Walsh masks $`\boldsymbol1_m+h_t+e_a`$: constant derivative requires $`v_a`$ to be independent
of $`a`$, and then $`Hv=0`$. Since $`m`$ is even, $`v=\boldsymbol1_m`$ has odd pairing with every
one of these masks. Hence (T3) requires $`h\le1`$.

The $`h=0`$ case also fails. Here $`C_2`$ is a monomial Clifford. If $`a\ne0`$,
the conjugated operator still has zero diagonal, giving norm distance at
least one from $`Q_s`$. If $`a=0`$, $`D_1`$ commutes with the first Pauli image and $`F`$
is a diagonal Pauli; $`Q_s`$ is not a Pauli and their norm distance is two.
Therefore a necessary condition for error strictly below one in (T1) is

```math
h=1,\qquad V^\perp=\mathrm{span}\{\boldsymbol1_m\},
\qquad u_0\cdot\boldsymbol1_m=1.                                 \qquad(T4)
```

Every other case has complete operator error at least one. The statement
includes arbitrary diagonal phases, arbitrary bank Clifford choices and
all bank columns, but it retains the source grammar (T1). In particular,
it is not a general restriction on the number of T layers.

<a id="common-affine-two-mask-screen"></a>

### 2.2. Two Fourier masks in one affine bank frame cannot meet the exception

Let $`L`$ be any invertible binary $`m\times m`$ matrix, let $`t`$ be any fixed bank
translation, and let $`U_{L,t}|z\rangle=|Lz+t\rangle`$. In that one common affine frame,
consider

```math
E=U_{L,t}D_2 H_A D_1 H_B D_0U_{L,t}^\dagger,            \qquad(T5)
```

where $`A`$ and $`B`$ are arbitrary sets of bank coordinates. The two masks need
not agree, and the diagonal middle and boundaries are unrestricted.
Every such word has $`\|EPE^\dagger-Q_s\|\ge1`$.

To prove this, the input parity in $`z`$ coordinates is a scalar times
$`Z^\ell`$ with $`\ell=L^{\mathsf T}\boldsymbol1_m`$. The unique nonzero target linear structure is
$`v_*=L^{-1}\boldsymbol1_m`$, still with derivative one. For $`C_2=H_A`$ the space $`V^\perp`$
is the span of the coordinate vectors indexed by $`A`$. The empty $`A`$ case is
the excluded $`h=0`$ case. Otherwise (T3) forces $`A=\{j\}`$ and $`v_*=e_j`$. But then
$`\ell_j=\boldsymbol1_m^{\mathsf T}Le_j=m=0\pmod2`$. The first mask $`H_B`$ therefore creates no
X component on $`j`$. Any X component it creates outside $`j`$ survives the
last $`H_j`$ and gives zero diagonal. If there is no such component, $`D_1`$
commutes with the first Pauli image and the result is a Pauli. Both cases
fail as above. A change of the two masks, arbitrary diagonal phases, and
a shared invertible affine exposure cannot produce the missing lookup.

The fixed affine walls are Clifford and can be emitted with $`O(m^2)`$
CNOTs plus $`O(m)`$ X gates and their actual inverses. This pays those walls,
not an arbitrary middle $`D_1`$. The exclusion remains valid regardless of
that middle's cost. Different affine frames between the mixers, active
extra registers, non-Clifford mixing stages and general source words are
outside (T5).

<a id="rank-one-mixer-lookup-derivative"></a>

### 2.3. The aligned rank-one case is real and exposes the missing phase

Condition (T4) must not be promoted to an impossibility claim. Define the
fixed involutory linear map $`L`$ by

```math
y_0=z_0,\qquad y_j=z_j+z_0\quad(j\ne0).
```

Its first column is $`\boldsymbol1_m`$. Since the target flips under complete bank
complement, there is a Boolean function $`g_s`$ of $`z_1,\ldots,z_{m-1}`$ with

```math
q_s(Lz)=z_0+g_s(z_1,\ldots,z_{m-1}).                   \qquad(T6)
```

Let $`R`$ be the involutory CNOT gather $`z_0\mapsto z_0+\sum_{j\ne0}z_j`$ and set

```math
\begin{aligned}
C_1&=H_{B_0}R,& C_2&=U_LH_{B_0},\\
D_s|z\rangle&=(-1)^{z_0g_s(z_{\ne0})}|z\rangle,\\
E_s&=C_2D_sC_1.
\end{aligned}                                                \qquad(T7)
```

The gather maps $`P`$ to $`Z_0`$, so $`C_1PC_1^\dagger=X_0`$. The exact complete
diagonal derivative is

```math
D_sX_0D_s^\dagger=(-1)^{g_s(z_{\ne0})}X_0.
```

The last $`H_{B_0}`$ and $`U_L`$ give $`E_sPE_s^\dagger=Q_s`$ on every bank column.
There is no initialized bank or omitted inverse. This realizes precisely
the rank-one exception. The Clifford envelope uses two Hadamards and
$`2(m-1)`$ CNOTs in each direction, zero T gates and no additional work.

Its missing operation is $`D_s`$. Formula (T7) merely writes the complete
lookup information into a controlled diagonal phase. Since $`q_s`$ has
Boolean degree $`p+1`$ and an invertible linear substitution preserves degree,
(T6) has $`\deg g_s=p+1`$. The displayed sign phase $`z_0g_s`$ has degree $`p+2`$.
No native sublogarithmic compiler or simultaneous resource ledger for that
middle is supplied. A general complex diagonal with the same derivative
can multiply (T7) by an arbitrary phase independent of $`z_0`$; the required
lookup-dependent derivative is unchanged.

Thus two genuinely different Clifford mixing stages are not excluded.
The complete transition calculation identifies exactly where this route
stores the unpaid operation: in a high-degree diagonal derivative.
Replacing that phase by a charged shallow construction would be new
progress; declaring its truth table, normalized spectrum, or formal
diagonal unitary does not improve the complete-indicator or frame frontier.

The [exact two-mixer tests](../../tests/test_two_mixer_lookup_screen.py) check the $`p=3`$ Walsh
support and every bank translation, a full Clifford-coset source fixture,
the same-frame slope mismatch, and all columns plus actual inverses of
the aligned rank-one construction. Its middle is evaluated as a truth-table
diagonal, explicitly not as a scalable native implementation. These finite
checks support the identities and boundaries; they do not replace the
all-size proof or supply the missing phase compiler.

<a id="rank-one-indicator-equivalence"></a>

### 2.4. The phase derivative is a complete-indicator problem

The rank-one exception admits a quantitative converse. This prevents the
phase in (T7) from being treated as an easier primitive merely because its
truth table is diagonal. Throughout this subsection $`p\ge3`$, $`m=2^p`$,
and every address, bank, work and reference input is arbitrary.
Let $`c=z_0`$ and introduce the *formal* string $`w`$ with
$`w_0=0`$ and $`w_j=z_j`$ for $`j\ne0`$. No physical zero input is introduced.
Since $`H\boldsymbol1_m=0`$ and $`m`$ is even, (T6) explicitly reads

```math
g_s(z_{\ne0})=\sum_{j\ne0}z_j+w_{s+Hw},\qquad
\Gamma_s=(-1)^{g_s(z_{\ne0})}.                         \qquad(T8)
```

**From a derivative to the full indicator.** Allow a native unitary V on
the address, bank and arbitrary work, with no diagonal or fixed-address
restriction. Suppose its complete relative action satisfies

```math
\|V X_0 V^\dagger-(\Gamma_s X_0)\otimes I_{\rm work}\|
 \le\eta.                                             \qquad(T9)
```

The source V need not return its work separately. Its actual inverse in
(T9) must be emitted, and the displayed complete operator promise includes
all work return. Let K be the Clifford shear
$`|x,y\rangle\mapsto|x+Hy,y\rangle`$, let $`H_B`$ denote the bank Hadamards,
and retain $`C_2=U_LH_{B_0}`$. Then the literal word

```math
\widetilde I_p=
 H_B Z_\pi K C_2 V X_0 V^\dagger C_2^\dagger K H_B       \qquad(T10)
```

approximates the complete indicator $`I_p`$, tensored with the identity on
work, within eta. Indeed $`C_2(\Gamma_s X_0)C_2^\dagger=Q_s`$,
$`KQ_sK=Z_\pi\Omega_x`$, and $`H_B\Omega_xH_B=I_p`$.
These are operator identities, so neither coherent address motion inside V
nor entanglement with a reference changes the estimate.
If t, d and G are common bounds for the emitted V and its actual native
inverse, with b work wires, one simultaneous ledger for (T10) is

```math
T\le2t,\qquad D_T\le2d,\qquad
G_{I}\le2G+pm+6m+4,\qquad b_I=b.                      \qquad(T11)
```

Here G counts **all** elementary native gates. Each K has pm/2 CNOTs;
the two $`C_2`$ walls cost 2m gates, the bank Hadamards cost 2m,
$`Z_\pi=S_B^2`$ costs 2m, and $`X_0=H_{B_0}S_{B_0}^2H_{B_0}`$ costs four.
No extra clean or dirty input is added by these walls.

In particular, an implementation within delta of $`D_s\otimes I`$ gives
$`\eta\le2\delta`$. The actual derivative promise (T9) is weaker than that
standalone approximation. At zero error its complete solution set is

```math
V=(D_s\otimes I)J,\qquad [J,X_0]=0.                   \qquad(T12)
```

This follows by multiplying the exact conjugation identity by
$`D_s^\dagger`$ on the left and D on the right. J may be nondiagonal, may
move the address and work, and need not return them on its own. Its action
cancels in the actual relative word. A phase independent of c is only the
diagonal part of this larger freedom. None of this makes J free to synthesize.

The [complete-indicator nullity bound](INDICATOR_PAULI_NULLITY.md) now gives
an unrestricted **count** consequence for every V satisfying (T9):
if $`\eta\lt1/\sqrt2`$, then $`2t\ge m-1`$, hence $`t\ge m/2`$.
For a standalone delta approximation of D, it is enough that
$`\delta\lt1/(2\sqrt2)`$. This uses the full indicator produced by (T10),
not a degree argument about a selected diagonal representation. It does
not imply a growing T-depth lower bound when dirty width may grow.

**From the indicator to the displayed phase.** Conversely, suppose a
native $`\widetilde I_p`$ approximates $`I_p\otimes I`$ within delta, with
common forward/actual-inverse bounds $`T_I,d_I,G_I,b_I`$ and zero clean
inputs. Add an arbitrary m-bit
bank R. Apply the Clifford shear $`s\mapsto x=s+Hw`$ and define

```math
A_c=\prod_{j\ne0}{\rm CCX}(c,z_j;R_j),\qquad
\Omega_x(R)=(-1)^{R_x}.                                \qquad(T13)
```

A leaves both $`R_0`$ and x unchanged. On every complete column,

```math
A_c\Omega_x(R) A_c\Omega_x(R)
       =(-1)^{c w_x}I_R.                               \qquad(T14)
```

In particular, the arbitrary value of $`R_0`$ cancels even when $`x=0`$.
Realize each Omega by bank Hadamards around the indicator. Use the actual
query adjoint in one of the two calls, giving
$`A_c\widetilde\Omega^\dagger A_c\widetilde\Omega`$.
The operator error is at most $`2\delta`$. Multiply by the Clifford phase
$`\prod_{j\ne0}{\rm CZ}(c,z_j)`$ and undo the address shear. Equations
(T8) and (T14) then give $`D_s\otimes I_{R,\rm work}`$ within $`2\delta`$.
No bank row was initialized, no postselection was used, and R and all query
work return in the same complete operator estimate.

The shared-control Toffoli bank A has a literal four-T-band construction
without any helper. Put $`k=m-1`$. Its Hadamard-conjugated sign phase is
$`(-1)^{c\sum_{j\ne0}z_jR_j}`$. The integer phase-polynomial identity is

```math
\begin{aligned}
4c\sum_{j\ne0}z_jR_j
&\equiv kc+\sum_{j\ne0}\bigl[z_j+R_j-(c\mathbin\oplus z_j)
 -(c\mathbin\oplus R_j)\\
&\hspace{42mm}-(z_j\mathbin\oplus R_j)
 +(c\mathbin\oplus z_j\mathbin\oplus R_j)\bigr]\pmod8.
\end{aligned}                                          \qquad(T15)
```

Since $`k\equiv7\pmod8`$, emit the c term as one T-dagger. In the first
band apply that gate together with T on every $`z_j,R_j`$. In the second,
compute both families $`c\oplus z_j`$ and $`c\oplus R_j`$ into their
respective wires, apply T-dagger to all of them, and undo the CNOTs.
The third band computes $`z_j\oplus R_j`$ into $`R_j`$, applies T-dagger,
and uncomputes. The fourth similarly computes
$`c\oplus z_j\oplus R_j`$ into $`R_j`$ and applies T. Every parity map
is invertible and restored; every physical T band uses distinct wires.
Shared CNOT controls incur Clifford scheduling, not extra T layers.
The phase core uses $`6k+1`$ T gates and 10k CNOTs. Hadamards on the k
R targets on both sides yield A with
$`T_A=6k+1`$, $`D_T(A)\le4`$, and $`G_A=18k+1`$.
The actual inverse has the same ledger.

Using these words in (T14), with a three-gate native implementation of
each CZ, gives the simultaneous reverse-reduction bounds

```math
\begin{aligned}
T_D&\le2T_I+12m-10,& D_T(D)&\le2d_I+8,\\
G_D&\le2G_I+pm+43m-37,& b_D&\le m+b_I,\\
\|\widetilde D-D_s\otimes I\|&\le2\delta.&
\end{aligned}                                          \qquad(T16)
```

The additive gate count includes pm CNOTs for the two address shears,
4m query Hadamards, two A words, and $`3(m-1)`$ gates for the parity phase.
All additional workspace is dirty. Thus, at fixed requested accuracy,
the displayed phase and the complete indicator have equivalent
sublogarithmic-depth targets up to constant call factors and additive
constant T-depth, with the same linear count/width and quadratic gate
orders. Equation (T12) leaves a genuinely different, potentially useful
source representative open, but paying for that representative is still
the complete transport problem. No general depth frontier changes here.

[Exact reduction fixtures](../../tests/test_rank_one_indicator_equivalence.py)
check the physical four-band phase polynomial and actual reverse,
all auxiliary-bank columns including row zero, the complete forward
operator, and a nondiagonal address/work-moving gauge. They verify these
reductions, not a new native implementation of D or the indicator.
