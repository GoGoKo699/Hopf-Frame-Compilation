# Joint transfer constraints for unrestricted two-layer queries

[Full-space transfer alphabet](SHALLOW_SOURCE_OBSTRUCTION.md#1-a-discrete-pauli-transfer-alphabet-at-t-depth-at-most-two) ·
[Initialized transfer and error](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#1-initialized-work-and-a-robust-transfer-witness) ·
[Read-only comparison](READ_ONLY_TWO_LAYER_RIGIDITY.md) ·
[Complete-frame reductions](BOOLEAN_QUERY_FRAME_REDUCTION.md)

With one exactly returned clean qubit, every occupied Pauli-transfer
coefficient of a native two-T-layer circuit lies in the same discrete
alphabet as a full-space two-layer coefficient. A common active axis of
the two one-band encoders forces one of the two apparent trace summands
to vanish. This excludes three-control AND and every three-address
Boolean cubic, and strengthens the exact affine-distance condition to
one quarter. R75's terminal-return rigidity makes the AND3 error floor
$`(3-2\sqrt2)/8`$ robust to leakage and approximate return.

The proof allows arbitrary Clifford interlayers, temporary address
changes and arbitrarily many dirty qubits. It is a constant-layer result,
not a growing depth lower bound or an improvement to R96. Sections 1–3
retain the weaker joint-row argument that led to the correlated-axis
theorem in Section 4; the later theorem supersedes their exact and
robust query bounds.

## 1. Model and the joint row constraint

Let S contain all logical inputs and outputs, D contain b arbitrary dirty
qubits, and c be one clean qubit. Put $`d=|S|+b`$, $`q=d+1`$,
and $`J|\psi\rangle=|\psi\rangle|0\rangle_c`$. The native circuit is

```math
V=C_2D_2C_1D_1C_0,
```

where each C is any Clifford word and each D is a physical layer of
T/T-dagger gates on distinct wires. No intermediate address restriction
is imposed. There are no measurements, resets or supplied states.
All T gates, Clifford gates, data wires and helpers belong to V; no
count or width bound beyond the two T layers and one clean qubit is
needed for the obstruction.

For Hermitian logical Paulis P,Q, write $`\widetilde P=P\otimes I_D`$
and $`\widetilde Q=Q\otimes I_D`$. Define the two full-space entries

```math
u=2^{-q}\mathrm{Tr}[(\widetilde P\otimes I_c)
                       V(\widetilde Q\otimes I_c)V^\dagger],
\qquad
v=2^{-q}\mathrm{Tr}[(\widetilde P\otimes I_c)
                       V(\widetilde Q\otimes Z_c)V^\dagger].
```

Expanding the clean-zero projector gives the initialized transfer
coefficient $`x_V(P,Q)=u+v`$, with no extra factor. The retained
two-layer overlap lemma gives

```math
u,v\in\mathcal A=\{0\}\cup\{\pm2^{-j/2}:j\ge0\}.
```

The new point is that u and v occupy distinct columns of the **same**
orthogonal Pauli-transfer row. Hilbert–Schmidt norm preservation gives

```math
u^2+v^2\le1.                                         \qquad(1)
```

In particular, a coefficient of magnitude one forces the other to
vanish. This does not assume clean return, exact query action, or any
classification of the Clifford interlayers. Also $`|u+v|\le1`$
because an initialized transfer coefficient is a normalized trace of
an operator of norm at most one. Define the necessary alphabet

```math
\mathcal S=\{u+v:u,v\in\mathcal A,\ u^2+v^2\le1,\ |u+v|\le1\}.
                                                               \qquad(2)
```

Membership in this set is necessary, not a native synthesis theorem.
The earlier independent-summand test $`\mathcal A+\mathcal A`$
misses (1): it accepts $`1-2^{-k}`$ via summands whose squared norm
exceeds one.

If the actual circuit obeys the full initialized-isometry guarantee

```math
\|VJ-J(U\otimes I_D)\|\le\eta,
```

the existing hybrid estimate, including clean leakage and dirty/reference
return, gives

```math
\left|2^{-|S|}\mathrm{Tr}(PUQU^\dagger)-(u+v)\right|\le2\eta.
                                                               \qquad(3)
```

Consequently any target entry t gives the valid floor
$`\eta\ge\tfrac12\mathrm{dist}(t,\mathcal S)`$. This proof witness
averages dirty inputs only after using the full operator contract. It
does not grant a maximally mixed input promise. Scalar phases do not
affect this necessary condition; the compiler still owes literal phases.

## 2. An exact affine-distance condition

This intermediate one-eighth condition is strengthened to one quarter
by the exact-return theorem in Section 4.

**Rational gap.** No element of $`\mathcal S`$ lies in
$`\mathbb Q\cap(3/4,1)`$ or its negative.

For a nonzero rational sum u+v, odd-exponent terms are irrational
rational multiples of $`\sqrt2`$. Two such terms can give a rational
sum only by canceling to zero. Thus the nonzero rational case uses
signed dyadic powers. If either term has magnitude one, (1) forces the
other to be zero. Otherwise each positive term is at most 1/2. Their
sum is at most 3/4 unless both are 1/2, which gives exactly one.
Opposite signs cannot increase the positive maximum. Negation gives
the other interval. This proves the gap without a finite-width cutoff.

Let $`Q_f`$ be an exact returned XOR query on a address bits and an
arbitrary output z:

```math
Q_f|x,z\rangle=|x,z\mathbin\oplus f(x)\rangle.
```

Taking output Z and its product with an address Z character as the
two logical probes gives every normalized Walsh coefficient

```math
W_f(\ell)=2^{-a}\sum_x(-1)^{f(x)+\ell\cdot x}.
```

These coefficients are rational. Therefore an exact one-clean native
two-layer query satisfies

```math
W_f(\ell)\notin (3/4,1)\cup(-1,-3/4)
\quad\hbox{for every }\ell.                           \qquad(4)
```

**Corollary.** If f is nonaffine, its relative Hamming distance from
every affine Boolean function is at least 1/8. Indeed the affine
distance is $`(1-(-1)^c W_f(\ell))/2`$; distance zero would make
f itself affine, and (4) excludes every distance strictly between zero
and 1/8. This is a necessary condition, not a converse and not a
degree bound. Balanced high-degree functions can satisfy it.

In particular, any nonzero predicate supported on less than one eighth
of the addresses is excluded. For
$`f_a(x)=\prod_{j=1}^a x_j`$, its bias is

```math
t_a=W_{f_a}(0)=1-2^{1-a}.
```

Every exact one-clean implementation with $`a\ge4`$ therefore needs
at least three physical T layers, independently of dirty width or
internal address mixing. Four-control AND here means four address bits
and a separate arbitrary XOR output. The cubic case $`a=3`$ has bias
3/4 and is not excluded by this argument.

## 3. A robust unrestricted four-control obstruction

These intermediate joint-row floors require no exact-return reduction.
Section 4 combines the native encoder correlation with R75 to improve
the AND4 floor to $`1/16`$ and also exclude AND3 robustly.

For four-control AND the exact distance is

```math
\mathrm{dist}(7/8,\mathcal S)=\frac{5\sqrt2-7}{8}.
                                                               \qquad(5)
```

To check the infinite alphabet, a summand of magnitude one gives only
the admissible sums $`\pm1`$. Otherwise order the magnitudes. If the
larger is at most $`1/(2\sqrt2)`$, the sum is at most $`1/\sqrt2`$
and too far below 7/8. The remaining possible positive larger terms
are $`1/\sqrt2`$ and 1/2; a negative smaller term is again too small.
For larger term 1/2, the nearest allowed smaller term to 3/8 is
$`1/(2\sqrt2)`$, giving $`(2+\sqrt2)/4`$. For larger term
$`1/\sqrt2`$, the nearest smaller term to $`7/8-1/\sqrt2`$
is $`1/(4\sqrt2)`$, giving $`5\sqrt2/8`$. Comparing these and
one proves (5); the attaining pair has squared norm $`17/32`$.

Thus every native circuit of T-depth at most two with at most one
clean qubit has full initialized-isometry error at least

```math
\eta_4=\frac{5\sqrt2-7}{16}>2^{-8}                    \qquad(6)
```

from four-control AND with all dirty work returned. Arbitrary controls,
Clifford address mixing, external references and approximate clean return
are included. Fewer clean qubits are covered by an untouched clean
spectator. Neither T-count nor Clifford count is restricted. At precision
$`L\ge8`$, this query consequently needs at least three T layers.

For completeness, the same averaged witness yields the following floors:

| AND arity a | Half the distance of $`t_a`$ to $`\mathcal S`$ |
|---|---:|
| 4 | $`(5\sqrt2-7)/16`$ |
| 5 | $`(8\sqrt2-11)/32`$ |
| 6 | $`(23-16\sqrt2)/64`$ |
| $`a\ge7`$ | $`2^{-a}`$ |

The same two largest nonunit summands suffice to check this table:
above $`(2+\sqrt2)/4`$, the only relevant candidates up to one are
$`5\sqrt2/8`$, $`(1+2\sqrt2)/4`$, and one. The floor decreases
with a in the last row. It is not a fixed-error, dimension-growing
depth lower bound. R96 and R97 are unchanged; the two-call frame
reduction does not turn this constant query exclusion into a stronger
complete-frame frontier.

## 4. A shared active axis restores the single transfer alphabet

**Theorem (one clean, two native T layers).** Suppose V returns its one
clean qubit exactly on every initialized logical and dirty input. Let
$`V_0`$ be its occupied unitary. Every Pauli-transfer coefficient of
$`V_0`$, including probes on dirty wires, belongs to

```math
\mathcal A=\{0\}\cup\{\pm2^{-k/2}:k\ge0\}.
```

No rationality, Booleanity or intermediate address-preservation premise
is imposed on $`V_0`$.

### Exact return and the common essential span

Exact initialized return makes the $`Z_c=+1`$ halfspace invariant.
Unitarity makes its orthogonal complement invariant as well, so
$`[V,Z_c]=0`$ on the full space. This implication uses one clean
qubit: preservation of a smaller multi-clean code need not preserve
each individual clean Pauli.

Factor the same physical circuit as

```math
V=U_B^\dagger U_A,
\qquad U_A=C_1D_1C_0,\qquad U_B=D_2^\dagger C_2^\dagger.
```

Each encoder has at most one physical T layer, including actual inverse
gates. Write $`U_j=E_jK_j`$, with $`K_j`$ Clifford and $`E_j`$
a Clifford conjugate of one physical band. Its signed Pauli axes are
independent and mutually commuting. Exact clean return gives the common
Hermitian Clifford

```math
S=U_AZ_cU_A^\dagger=U_BZ_cU_B^\dagger.
```

Put $`R_j=K_jZ_cK_j^\dagger`$. The clean-active axes are those
anticommuting with $`R_j`$. The
[intrinsic-axis lemma](ONE_CLEAN_MIDPOINT_RIGIDITY.md#the-essential-axes-are-intrinsic)
gives, for their binary labels $`a_i`$,

```math
F_S(v)=v+\sum_i\langle a_i,v\rangle a_i,
\qquad L:=\mathrm{im}(F_S-I)=\mathrm{span}\{a_i\}.
```

The equality of the image and span uses physical independence; signs
and T versus T-dagger do not change it. Thus both encoders have the same
active span L. Every axis in either full band commutes with every Pauli
label in L, including the axes silent on its clean core.

### The two traces cannot both survive

For Hermitian Paulis P,Q on the remaining wires, lift them trivially
over c and define

```math
A=U_A Q U_A^\dagger,\qquad
B=U_B P U_B^\dagger,\qquad G=BA.
```

A and B are Hermitian Cliffords and each commutes with S. With
$`d_{\rm phys}=2^q`$, the occupied transfer coefficient is

```math
x_{V_0}(P,Q)
=\frac{2}{d_{\rm phys}}\mathrm{Tr}\!\left[\frac{I+S}{2}BA\right]
=\frac{\mathrm{Tr}G+\mathrm{Tr}(GS)}{d_{\rm phys}}.
```

Both traces are real by Hermiticity and cyclicity. Suppose L is nonzero
and choose a signed clean-active axis H of the first encoder. Its label
belongs to both active spans, so H commutes with both $`E_A`$ and
$`E_B`$. Since it anticommutes with $`R_A`$, it obeys
$`SHS=-H`$. Conjugating H by either A or B only applies the
commutation sign of that encoder's Pauli core:

```math
AHA=\epsilon_A H,\qquad BHB=\epsilon_B H,
\qquad \epsilon_A,\epsilon_B\in\{\pm1\}.
```

For $`\epsilon=\epsilon_A\epsilon_B`$ this gives

```math
HGH=\epsilon G,\qquad HGSH=-\epsilon GS.
```

Trace invariance forces $`\mathrm{Tr}G=0`$ when
$`\epsilon=-1`$ and $`\mathrm{Tr}(GS)=0`$ when
$`\epsilon=1`$. At most one summand survives. Its normalized trace
is a real Clifford trace, whose magnitude is zero or $`2^{-k/2}`$.

If L vanishes, S is a signed Pauli. A and B preserve its +1 stabilizer
code; after a Clifford identification of that code with the remaining
logical space, their restrictions are Hermitian logical Cliffords.
The occupied coefficient is again a real normalized Clifford trace.
This proves the theorem in both cases without an additional physical
construction or any extra initialized work.

The correlation of the two native encoders is essential. An arbitrary
Hermitian Clifford midpoint and an arbitrary Clifford commuting with it
do not supply the common-axis sign action used here.

### Exact Boolean consequences

For an exact returned XOR query, all Walsh coefficients are rational.
The rational part of $`\mathcal A`$ is
$`\{0\}\cup\{\pm2^{-k}:k\ge0\}`$. If f is nonaffine, no
Walsh coefficient has magnitude one, so all have magnitude at most 1/2.
Therefore its relative Hamming distance from every affine Boolean
function is at least **one quarter**. This strengthens Section 2 and
excludes every nonzero predicate supported on less than one quarter
of the addresses. It is a necessary condition, not a converse or an
unrestricted degree classification.

In particular, AND3 has bias 3/4 and is excluded. More generally, any
three-address cubic has odd truth-table weight. Its eight Walsh
coefficients are odd multiples of 1/4; Parseval forces exactly one
coefficient of magnitude 3/4 and seven of magnitude 1/4. Every such
cubic is excluded by the same theorem.

### Robust floors include leakage and work return

[R75, Section 6](FOUR_LAYER_OPERATOR_RIGIDITY.md#6-one-clean-qubit-terminal-sectors-only)
shows that initialized-isometry error $`\eta\lt\tau`$ against a
clean-returning target forces the **same physical circuit** to preserve
$`Z_c`$ exactly whenever its T-depth is at most four, where

```math
\tau=\frac{\sin(\pi/8)}{2\sqrt2}.
```

Its argument bounds the leakage block by eta, the reverse block by the
same norm through unitarity, and hence $`\|[V,Z_c]\|\le2\eta`$.
No replacement circuit or increase in depth is involved. It applies
in particular to the present two-layer circuit and allows a freely
optimized common phase.

Below tau, exact return lets us use the theorem. Occupied unitaries
within eta up to common phase have Pauli images within $`2\eta`$,
so every target transfer entry t obeys

```math
\eta\ge\min\!\left\{\tau,\frac12\mathrm{dist}(t,\mathcal A)\right\}.
```

For AND3, the nearest alphabet value to 3/4 is $`1/\sqrt2`$.
For every AND arity $`a\ge4`$, the nearest value to
$`1-2^{1-a}`$ is one. Consequently the unconditional
initialized-isometry floors are

```math
\eta\ge
\begin{cases}
(3-2\sqrt2)/8\approx0.0214466,&a=3,\\
2^{-a},&a\ge4.
\end{cases}
```

All these values lie below tau. The AND3 floor also holds for every
three-address cubic, using its Walsh coefficient of magnitude 3/4.
Thus accuracy $`2^{-L}`$ with integer $`L\ge6`$ requires at least
three layers for these cubic queries; AND4 needs three layers for
$`L\ge5`$. Dirty width, T-count and elementary Clifford count are
unrestricted. Literal phases, all occupied columns, arbitrary external
references and complete workspace return remain in the model. These
fixed-depth obstructions do not give a dimension-growing depth law.

## 5. Exact rational return has a native Galois constraint

There is a second restriction that retains the physical independence of
the two T bands. Let $`\Pi=(I+Z_c)/2`$ and suppose V returns the one
clean qubit exactly, with occupied unitary $`U_0`$ on all remaining
logical and dirty wires. Assume the entries of $`U_0`$ lie in
$`\mathbb Q(i)`$. Every literal XOR query with returned dirty work
satisfies this premise. No address-preservation condition is imposed.

Let sigma be the automorphism of $`\mathbb Q(\sqrt2,i)`$ defined by

```math
\sigma(\sqrt2)=-\sqrt2,\qquad \sigma(i)=i.
```

**Galois return theorem.** Under these exact hypotheses,

```math
V^\dagger\sigma(V)\in\{I,Z_c\}.                       \qquad(7)
```

The statement concerns the actual circuit matrix, including its scalar
phase. It holds on all physical columns, although the target is specified
only on the clean-zero code.

### Signed relative operator

The native gates satisfy the exact identities

```math
\begin{aligned}
\sigma(H)&=-H,&\sigma(S)&=S,&\sigma(\mathrm{CX})&=\mathrm{CX},\\
\sigma(T)&=ZT,&\sigma(T^\dagger)&=ZT^\dagger.
\end{aligned}
```

Thus every actual Clifford word obeys
$`\sigma(C_j)=\epsilon_j C_j`$, with its literal
$`\epsilon_j\in\{\pm1\}`$. Let $`P_j`$ be the product of Z over
the active physical wires of $`D_j`$, irrespective of T orientation,
and set $`\epsilon=\epsilon_0\epsilon_1\epsilon_2`$. Then
$`\sigma(D_j)=P_jD_j`$. For the signed Hermitian Pauli
$`Q=C_1^\dagger P_2C_1`$, actual adjacent cancellations give

```math
G:=V^\dagger\sigma(V)
 =\epsilon C_0^\dagger(D_1^\dagger QD_1)P_1C_0.       \qquad(8)
```

A one-band Pauli image is a Hermitian Clifford, so G is a Clifford
with its actual scalar. Sigma fixes the occupied target and commutes
with adjoint. Consequently $`G\Pi=\Pi`$ literally.

Unitarity now gives $`G=\mathrm{diag}(I,W)`$ in the two clean
sectors. Since $`GX_cG^\dagger`$ is a signed Hermitian Pauli
anticommuting with $`Z_c`$, its off-diagonal blocks force

```math
G=\begin{pmatrix}I&0\\0&\lambda P\end{pmatrix},
\qquad \lambda\in\{\pm1,\pm i\},                    \qquad(9)
```

for a Hermitian Pauli P on the remaining wires, possibly identity.
Every entry of (9) is fixed by sigma. On the other hand, its definition
gives $`\sigma(G)=\sigma(V)^\dagger V=G^\dagger`$. Therefore G
is Hermitian and $`\lambda=\pm1`$. This is a controlled signed
Pauli; no action on the unused code has been assumed in deriving it.

### Independent native axes remove the controlled-Pauli alternative

For a binary symplectic Clifford map F define

```math
q_F(v)=\langle v,(F-I)v\rangle.
```

Under symplectic coordinate change S,
$`q_{S^{-1}FS}(v)=q_F(Sv)`$. Scalar phases and Paulis do not change F.
For (9), let c denote the clean Z label and p the logical P label.
The binary action is

```math
F_G(v)=v+\langle c,v\rangle p+\langle p,v\rangle c.
```

The two contributions to $`q_{F_G}`$ cancel over the binary field,
so this form vanishes identically, including when P is identity.

Now let $`a_1,\ldots,a_h`$ be exactly the first physical band's Z
axes anticommuting with Q. They are independent and mutually commuting.
The symplectic action of $`D_1^\dagger QD_1`$ is their transvection
product, independently of T versus T-dagger orientations:

```math
\begin{aligned}
F(v)&=v+\sum_{j=1}^h\langle a_j,v\rangle a_j,\\
q_F(v)&=\sum_{j=1}^h\langle a_j,v\rangle^2
       =\left\langle\sum_{j=1}^h a_j,v\right\rangle.
\end{aligned}                                                     \qquad(10)
```

Equation (8) changes this action only by the fixed $`C_0`$ coordinates.
The vanishing form and nondegeneracy of the symplectic pairing give
$`\sum_j a_j=0`$. Physical independence forces h to vanish. Thus Q
commutes with $`D_1`$, and (8) is a Pauli including its scalar. The
only Paulis equal to identity on every clean-zero column are I and
$`Z_c`$, proving (7).

The independence step cannot be replaced by dependent phase-polynomial
axes. The commuting triangle a,b,a+b has zero sum but can have a
nonidentity transvection product. For distinct Z axes its S-phase
product is a signed CZ, whose displacement form is zero. This explains
why a descent to arbitrary dependent phase-polynomial layers loses a
constraint needed here.

### Rational columns and the total-axis condition

For every physical Pauli R commuting with $`Z_c`$, equation (7) gives

```math
\sigma(VRV^\dagger)=VRV^\dagger.
```

Every transfer coefficient in these columns is real and sigma-fixed,
hence rational. Combining this with the full two-layer alphabet yields

```math
c_V(P,R)\in\{0\}\cup\{\pm2^{-k}:k\ge0\},
\qquad [R,Z_c]=0.                                  \qquad(11)
```

This includes both columns entering the occupied logical coefficient in
Section 1. It excludes the odd square-root powers for exact rational
return. In native coordinates it also gives the signed total-axis
constraint

```math
[Q,a_j]=0\ \hbox{for every first-band axis},\qquad
\epsilon C_0^\dagger QP_1C_0\in\{I,Z_c\}.           \qquad(12)
```

The pulled-back total second-band label therefore equals the total
first-band label or differs by the incoming clean Pauli label. Equation
(12) retains the actual signs.

The theorem is exact-only. Exact clean return alone does not make the
occupied target rational, and approximation to a rational target does
not commute with sigma. Thus (7) and (11) provide no additional robust
rounding radius beyond the independent argument in Section 4.

## 6. Scope fixtures for the Galois theorem

**Code-flipping columns may be irrational.** A controlled-H with c as
control is identity on clean zero and satisfies
$`V^\dagger\sigma(V)=Z_c`$. Its target-Pauli columns have rational
half coefficients, but

```math
VX_cV^\dagger=X_c\otimes H
 =\frac{X_cX_t+X_cZ_t}{\sqrt2}.
```

This has an actual two-T word: put $`A=SHTHS^\dagger`$ and execute
chronologically $`A^\dagger,\mathrm{CZ},A`$, with A on the target.
The identity $`AZA^\dagger=H`$ fixes the literal phase. Expanding
S-dagger and CZ into the fixed H/S/CNOT alphabet gives two T/T-dagger
gates and fifteen elementary Cliffords. Its actual reversed word also
has two T layers; expansion of its eight S-daggers gives thirty-one
elementary Cliffords. There are one clean control, one arbitrary target,
and no helpers. It is not valid to extend (11) to code-flipping columns.

**Rationality does not imply C3.** On two arbitrary wires a,b define

```math
C=\mathrm{CX}_{a\to b}H_a\mathrm{CX}_{a\to b},\qquad
D=T_aT_b,\qquad U=CDCD.
```

Its chronological word is D,C,D,C. It uses four T gates, six elementary
Cliffords and two physical T layers, without helpers. Its actual inverse
is $`D^\dagger CD^\dagger C`$, with the same ledger. Exact matrix
arithmetic gives $`\sigma(U)=U`$. Nevertheless, with
$`P=X_a`$ and $`W=UPU^\dagger`$,

```math
\begin{aligned}
W&=\tfrac12(X_a+Y_a+Z_aX_b+Z_aY_b),\\
WPW^\dagger&=\tfrac12(-X_a+Y_a+Z_aX_b+Z_aY_b).
\end{aligned}
```

The second line is not a Pauli, so U is not C3. An idle clean qubit
supplies exact rational occupied return with G equal to identity.
U is nonmonomial and is not a Boolean-query counterexample. The example
rules out inferring source C3 or a Boolean degree bound from (7) alone.

**Two initialized clean qubits change the conclusion.** Put
$`W=THT`$ and let $`B=\mathrm{CX}_{0\to1}(H_0\otimes I_1)`$
prepare the Bell state. Then

```math
V=B^\dagger(W\otimes W^*)B
```

fixes $`|00\rangle`$ with literal phase one, because
$`W\otimes W^*`$ fixes the Bell state for every unitary W. Its
chronological native word is

```math
H_0,\ \mathrm{CX}_{0\to1},\ (T_0T_1^\dagger),\
H_0H_1,\ (T_0T_1^\dagger),\ \mathrm{CX}_{0\to1},\ H_0.
```

It has four T/T-dagger gates, four H and two CNOT gates in two T layers.
The actual inverse has the same counts and returns the same initialized
pair. Arbitrary additional data and references are untouched. Here
$`V^\dagger\sigma(V)`$ is a non-Pauli Clifford. Moreover, the full
transfer coefficient with input $`Z_0`$ and output $`X_0X_1`$ is
$`\sqrt2/4`$, although the input commutes with both clean Z checks.
Thus neither (7) nor its rational-column corollary extends to two clean
qubits. The fair target baseline for this two-clean identity fixture is
the empty circuit; it is a scope test, not a synthesis improvement.

## 7. A cubic boundary for the full scalar alphabet and finite checks

Three-control AND is excluded by Section 4. A genuinely cubic query on
more address bits can nevertheless pass the **entire** scalar transfer
alphabet. Consider

```math
f(x,y)=x\cdot y+x_1x_2x_3,
\qquad x,y\in\mathbb F_2^3.
```

This is cubic and bent: summing over y gives
$`|W_f(w)|=1/8`$ for every w. To check every Pauli-transfer entry
of its XOR query, let a be the input Pauli's address translation and t
its target-Z label. Nonzero output entries retain both labels. Let w be
the difference of address-Z labels and r the difference of target
translation labels. Up to the fixed sign from Hermitian Pauli
conventions, the coefficient is

```math
\frac12\left[W_{tf}(w)+(-1)^r W_{tf+\Delta_a f}(w)\right],
\qquad \Delta_a f(z)=f(z)+f(z+a).
```

For $`t=1`$, the second function is $`f(z+a)`$, whose Walsh
coefficient is $`(-1)^{w\cdot a}W_f(w)`$. The resulting coefficient
is zero or has magnitude 1/8. For $`t=0`$ and nonzero
$`a=(u,v)`$,

```math
\Delta_a f(x,y)=u\cdot y+v\cdot x+u\cdot v
                 +\Delta_u(x_1x_2x_3).
```

This derivative is balanced. When u vanishes it is nonzero linear;
otherwise its quadratic polar matrix on x is

```math
\begin{pmatrix}
0&u_3&u_2\\u_3&0&u_1\\u_2&u_1&0
\end{pmatrix},
```

which has rank two. Thus the nonzero derivative Walsh magnitudes are
respectively one or 1/2. The displayed transfer formula, including
$`a=0`$, proves that every transfer entry has magnitude in

```math
\{0,1,1/2,1/4,1/8\}\subset\mathcal A.
```

This is an exact boundary certificate, **not** a native two-layer
construction. Applying the scalar alphabet separately to every entry
cannot exclude this cubic candidate. The concrete resumption question
is whether its full XOR permutation admits correlated one-band encoders
with the same returned clean code, or whether a simultaneous encoder
invariant excludes it. A resolution must retain physical independence
and the full occupied action; rationality alone does not imply C3, and
arbitrary dependent phase-polynomial descent loses native information.

The [joint-transfer tests](../../tests/test_joint_transfer_two_layer.py)
check the coefficient normalization, same-row norm restriction, bounded
alphabet/error witnesses, exact returned native encoder fixtures with
non-Pauli midpoints and different active bases, and the cubic bent
query's full scalar-alphabet boundary. The seven
[exact Galois tests](../../tests/test_two_layer_galois_return.py) use
$`\mathbb Q(\sqrt2,i)`$ arithmetic for the native automorphism,
controlled-H and its actual inverse, the two-layer dirty Toffoli with
and without an unused clean phase, independent and dependent axes,
all signed controlled-Pauli extensions, the rational non-C3 word, and
the two-clean Bell-code failure. Every fixture includes the relevant
occupied and unused columns. These checks support the identities and
scope boundaries; they are not an asymptotic compiler or an unrestricted
degree classification.
