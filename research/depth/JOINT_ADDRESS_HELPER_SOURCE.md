# Joint address, carrier and helper mixing on the first quartic face

[Ordered matching source](ORDERED_MATCHING_PROBE_SOURCE.md) ·
[Two outer cubic conjugations](CLIFFORD_INVOLUTION_SOURCE.md#two-outer-conjugations) ·
[Returned cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md)

A paid twenty-band trial mixes address, carrier and arbitrary helper
coordinates through Bell changes and swaps. It creates the required quartic
matching phase, but also moves carrier bits. Two anticommuting Pauli
coefficients prove that no final Clifford wall completes this trial.
The actual missing carrier correction is paid below, giving a complete
finite-face phase circuit. The [nonmonomial Hadamard-pair screen](#nonmonomial-hadamard-pair-screen)
below excludes a distinct first coupling for arbitrary later Clifford walls
and cubic-sign phases; an explicit weighted-phase word escapes its coefficient
screen. None supplies an all-size shallow lookup source or changes R96.

## 1. The complete quartic face and the preceding rank failure

Use four arbitrary address bits and four arbitrary carrier bits
$`r,c_1,c_2,d`$. Write

```math
a=x_1x_2,\qquad b=x_3x_4,\qquad
Q_\square=(-1)^{r+a c_1+b c_2+ab d}.                    \qquad(1)
```

This is the coefficient-lookup face with only the empty, 12, 34 and
1234 coefficients retained. In the full bank, setting all other coefficients
to zero gives a legitimate test subspace, not an initialized construction
resource. The finite target (1) is required on the entire eight-qubit input space,
tensored with identity on every arbitrary helper and reference.

A preceding trial used $`K_0=Z_r\mathrm{CX}(d;t)`$ and
$`D_0=\mathrm{CCZ}(t,u,v)`$, with t,u,v arbitrary helper bits. For any
unitary W on the address and u,v registers, the probe after these two
conjugations and a final $`H_t`$ is

```math
A_0=Z_r\mathrm{CZ}(d,t)\,
             W(-1)^{duv}W^\dagger.                    \qquad(2)
```

Here W holds r,d,t and both coefficient bits $`c_1,c_2`$; it may move
the address and both helpers jointly, or act on further arbitrary work.
The cubic identity follows by toggling t by d inside $`tuv`$. On the
invariant sector $`r=t=c_1=c_2=0,d=1`$, the actual negative eigenspace
has one quarter of the address/helper dimension. The target has one
sixteenth, since $`ab=1`$ at only one of the sixteen addresses. A larger
negative eigenspace has a nonzero intersection with the target positive
eigenspace, so their operator distance is two. No address/helper unitary W
repairs (2) while those carrier sectors remain fixed. This is a statement
about that coupling, not a restriction on general joint sources.

For the more specific $`W=BQ_{a,b}B^\dagger`$, where
$`Q_{a,b}:(u,v)\mapsto(u+a,v+b)`$ and B acts only on the address,
the residual phase is $`d(u+a)(v+b)`$ in B's address coordinates.
At $`u=1,v=0,d=1`$ its negative address rank is three, versus the
target rank one. Adding an address Fourier wall therefore already fails
on a single dirty sector. The dimension argument above also permits
W to mix u and v, showing why the next trial moves carrier coordinates.

<a id="phased-swap-joint-trial"></a>

## 2. A changed Bell-and-swap word creates the quartic phase

Retain (1), add two arbitrary helpers $`u_1,u_2`$, and let

```math
K=Z_r\,\mathrm{SWAP}(u_1,c_1)\mathrm{SWAP}(u_2,c_2).
```

K is a traceless Hermitian Clifford. For a pair u,c define the Bell
Clifford $`B_{uc}=\mathrm{CX}(u;c)H_u`$, in rightmost-first order.
The source

```math
U=B_{u_1c_1}B_{u_2c_2}\,
       \mathrm{CCX}(u_1,c_1;r)\mathrm{CCX}(u_2,c_2;r)
\quad\text{satisfies}\quad UZ_rU^\dagger=K.             \qquad(3)
```

Indeed each shared-target CCX adds its pair CZ to the Z probe, and
$`B_{uc}\mathrm{CZ}(u,c)B_{uc}^\dagger=\mathrm{SWAP}(u,c)`$.
This is an explicit four-band instance of the balanced-Clifford source.

Use two outer cubic-sign batches and fixed Clifford walls:

```math
\begin{aligned}
D_2&=\mathrm{CCZ}(u_1,x_1,x_2)\mathrm{CCZ}(u_2,x_3,x_4),\\
C_2&=H_{x_2}H_{x_4}B_{u_1c_1}^\dagger B_{u_2c_2}^\dagger,\\
D_3&=\mathrm{CCZ}(d,x_2,x_4),\\
C_3&=H_{c_1}H_{c_2}\mathrm{SWAP}(x_2,c_1)
                         \mathrm{SWAP}(x_4,c_2),\\
V&=C_3D_3C_2D_2U,\qquad A=VZ_rV^\dagger.
\end{aligned}                                                   \qquad(4)
```

These walls genuinely mix the address and carrier, and the Bell changes
mix the helpers and carrier. Every input remains arbitrary. There is no
projection onto a Bell outcome or a zero helper sector.

The complete resulting observable is

```math
A=Z_r\mathrm{CZ}(u_1,x_2)\mathrm{CZ}(u_2,x_4)
       (-1)^{dab}Z_{c_1}^{a}Z_{c_2}^{b}
                    X_{c_1}^{db}X_{c_2}^{da}.           \qquad(5)
```

To derive it, first conjugate one phased swap. A phase $`Z_u^a`$
around $`\mathrm{SWAP}(u,c)`$ gives
$`(Z_uZ_c)^a\mathrm{SWAP}(u,c)`$. The inverse Bell transform maps
the swap to CZ and $`Z_uZ_c`$ to $`Z_c`$. The two address Hadamards
therefore turn the result into

```math
C_2D_2KD_2C_2^\dagger
=Z_r\mathrm{CZ}(u_1,c_1)\mathrm{CZ}(u_2,c_2)
       \mathrm{CCX}(c_1,x_1;x_2)\mathrm{CCX}(c_2,x_3;x_4).
```

Conjugating the last two shears by $`D_3`$ adds the phase
$`d(c_1x_1x_4+c_2x_3x_2+c_1c_2x_1x_3)`$. The swaps in $`C_3`$
change the shears to $`c_1\mapsto c_1+a,c_2\mapsto c_2+b`$ and that
phase to $`d(a c_2+b c_1+ab)`$. Its phase commutes with those two
translations. The final carrier Hadamards give (5), including its literal
$`(-1)^{dab}`$ sign.

Equivalently, on every computational input the output has the target
phase from (1), an extra sign $`(-1)^{u_1x_2+u_2x_4}`$, and carrier
motion $`(c_1,c_2)\mapsto(c_1+db,c_2+da)`$. The correct basis phase
alone is therefore insufficient for the complete observable.

### A common forward and actual-inverse native ledger

The two shared-root CCZ gates inside (3), enclosed by $`H_r`$, cost
12 T gates, at most four bands, and at most 25 Clifford gates in a
common forward/inverse allowance. The two Bell words add four Cliffords.
Thus U costs at most 41 native gates. The returned-dirty cubic compiler
pays $`D_2`$ at 16 T gates and 32 total gates, using two arbitrary
private masks; $`D_3`$ costs eight T gates and sixteen total gates,
reusing one mask after return. $`C_2,C_3`$ cost six and eight Cliffords.
Consequently the literal candidate and its actual inverse obey

```math
T_V\le36,\qquad D_T(V)\le20,\qquad G_V\le103,
\qquad w_{\rm extra}=4,\qquad a_{\rm clean}=0.          \qquad(6)
```

Extra width consists of $`u_1,u_2`$ and the two returned masks, beyond
the eight logical inputs. Native emission error is zero. Equation (6)
does not assert small error relative to the target: the next section
excludes it. The actual inverse returns the full source input, and the
Gray masks return at every batch. The remaining helper dependence in
(5) is retained as part of the failed complete action.

<a id="anticommuting-support-gap"></a>

## 3. No final Clifford can remove the unwanted carrier action

For a Hermitian Pauli P, write
$`\widehat A(P)=2^{-n}\mathrm{Tr}(PA)`$, where n counts every physical
wire. Two coefficients of (5) are

```math
P=Z_rZ_{c_1}Z_{u_1},\qquad
R=Z_rX_{c_1}Z_{c_2}Z_{u_2},\qquad
\widehat A(P)=\widehat A(R)=1/16,\qquad PR=-RP.         \qquad(7)
```

For P, the carrier Pauli in (5) is $`Z_{c_1}`$ precisely when
$`a=1,b=0,d=0`$. Averaging the spectator CZ phases against $`Z_{u_1}`$
then selects $`x_2=1,x_4=0`$. Thus $`x_1=1`$, $`x_3`$ is free,
and d is zero, giving two of the 32 address/d choices. For R, the
carrier Pauli is $`X_{c_1}Z_{c_2}`$ when $`a=0,b=1,d=1`$.
The helper average selects $`x_2=0,x_4=1`$, leaving $`x_1`$ free
and forcing $`x_3=1`$, again two of 32. No literal minus sign occurs
in either selected branch.

Every computational-diagonal operator has only mutually commuting
Paulis in its expansion. Clifford conjugation preserves both coefficient
magnitudes and Pauli commutation. Therefore, for every final Clifford F,
at least one of $`FPF^\dagger,FRF^\dagger`$ is absent from the diagonal
target's Pauli support. Normalized trace against that Pauli gives

```math
\boxed{\ \|FAF^\dagger-Q_\square\otimes I\|\ge1/16.\ } \qquad(8)
```

The bound permits F to mix the address, carrier and every helper, and to
use further arbitrary wires. Identity extensions preserve (7). Replacing
the particular final wall in (4) by any other Clifford is included by
absorbing the old wall's inverse into F. This is a uniform error floor
for the specified two-cubic trial, not a theorem about all two-outer-
conjugation choices. In particular, changing the earlier cubic or Clifford
walls leaves its hypotheses.

## 4. The exact missing correction is paid on this finite face

Let $`J_0=\mathrm{CZ}(u_1,x_2)\mathrm{CZ}(u_2,x_4)`$, and let G
toggle $`c_1`$ by db and $`c_2`$ by da. From (5),

```math
A=J_0Q_\square G,\qquad [G,Q_\square]=0,\qquad
GJ_0A=Q_\square\otimes I_{\rm work}.                  \qquad(9)
```

The commutation follows because the phase change is
$`a(db)+b(da)=0`$. The two helper CZ gates cancel the actual helper
phase before their wires are reused. They cost six Cliffords.

Use those now returned $`u_1,u_2`$ as arbitrary outer masks. Define

```math
U_0=\mathrm{CCX}(x_3,x_4;u_1)\mathrm{CCX}(x_1,x_2;u_2),
\qquad
W_0=\mathrm{CCX}(d,u_1;c_1)\mathrm{CCX}(d,u_2;c_2).
```

The chronological word $`U_0,W_0,U_0,W_0`$ implements G and returns
both outer masks: the two occurrences of each old helper value cancel.
Each two-CCX batch costs at most 16 T gates, eight bands and 36 total
gates, using the same two separate returned Gray masks. Thus G costs
at most 64 T gates, 32 bands and 144 total gates.

Execute the actual probe word $`V^\dagger,Z_r,V`$, then $`J_0`$, then
this correction. Its exact complete finite-face ledger is

```math
T\le136,\qquad D_T\le72,\qquad G_{\rm all}\le358,
\qquad w_{\rm extra}=4,\qquad a_{\rm clean}=0,\qquad\epsilon=0.
                                                               \qquad(10)
```

Z is paid as two S gates. All helpers and arbitrary reference correlations
return in (9); reversing the actual word gives the same common bounds.
This pays the repair rather than calling the unwanted translation free.
It does not meet the three-block discovery target, implement the other
coefficient terms of the full p=4 bank, or supply an all-size recurrence.
The new design question is to cancel these carrier translations within
a genuinely changed joint word while retaining the quartic sign.

The [exact joint-source fixtures](../../tests/test_joint_address_helper_source.py)
check the complete operator identity, both noncommuting Pauli coefficients,
literal native pieces, actual inverses and the returned finite-face repair.
These fixed checks support the displayed analytic identities, not an
asymptotic compiler.

<a id="nonmonomial-hadamard-pair-screen"></a>

## 5. A nonmonomial Hadamard-pair source still needs a changed coupling

The Clifford-symmetry screen need not constrain every computational Z
probe of a nonmonomial first image. On u,v, the Hermitian Clifford
$`H_uH_v`$ exchanges the Z and X algebras, whose intersection consists
only of the identity. This is a different first image from the paired
swaps in Section 2. Its complete action nevertheless has a separate
obstruction for the following natural outer coupling.

Keep all eight logical wires of (1), two arbitrary helpers u,v and any
further arbitrary work. Define

```math
K_H=Z_rH_uH_v,\qquad
D_2=\mathrm{CCZ}(u,x_1,x_2)\mathrm{CCZ}(v,x_3,x_4),
\qquad R_H=D_2K_HD_2.
                                                               \qquad(11)
```

The first image $`K_H`$ is traceless, Hermitian and Clifford. Let
$`C_2,C_3`$ be arbitrary Cliffords on the whole physical register, with
no address, carrier or helper preservation requirement. Let
$`D_3|z\rangle=(-1)^{g(z)}|z\rangle`$, where g is any Boolean
polynomial of degree at most three on all those wires. Thus D_3 may be
an arbitrary cubic-sign batch, not just the single CCZ in Section 2.
For

```math
B_H=C_3D_3C_2R_HC_2^\dagger D_3^\dagger C_3^\dagger
```

the complete error obeys

```math
\boxed{\ \|B_H-Q_\square\otimes I_{\rm work}\|\ge1/16.\ }
                                                               \qquad(12)
```

The claim is about the specified K_H and D_2 with arbitrary later
Clifford mixing and cubic-sign D_3. It does not exclude general weighted
cubic phases, other first images or other first couplings. No source is
silently controlled, initialized or replaced by a projective inverse.

### Two commuting support classes anticommute across the partition

Conjugation in (11) gives the exact Hermitian operator

```math
R_H=\frac{Z_r}{2}\left[
 Z_uZ_v+Z_uX_v\mathrm{CZ}(x_3,x_4)
 +X_uZ_v\mathrm{CZ}(x_1,x_2)
 +X_uX_v\mathrm{CZ}(x_1,x_2)\mathrm{CZ}(x_3,x_4)\right].
                                                               \qquad(13)
```

Use $`\mathrm{CZ}(i,j)=(I+Z_i+Z_j-Z_iZ_j)/2`$. Its complete
Pauli support consists of 25 distinct labels in two classes:

| Class | Labels and coefficient magnitudes |
| --- | --- |
| E | One label $`e=Z_rZ_uZ_v`$ of magnitude 1/2, and the sixteen labels $`Z_rX_uX_v\prod_{i=1}^4Z_{x_i}^{s_i}`$, each of magnitude 1/8 |
| O | Four labels $`Z_rX_uZ_vZ_{x_1}^{s_1}Z_{x_2}^{s_2}`$ and four labels $`Z_rZ_uX_vZ_{x_3}^{s_3}Z_{x_4}^{s_4}`$, each of magnitude 1/4 |

Each class is mutually commuting, while every E label anticommutes
with every O label. The sixteen-label part of E is an affine four-flat
L in binary Pauli-label space, and $`e\notin L`$. Identity extensions
on arbitrary work preserve these statements, as does simultaneous
Clifford conjugation. Coefficient signs from the CZ expansions are
retained in (13); only the upper estimates below discard their signs.

### A cubic-sign Pauli image has a commuting affine support

For a Hermitian Pauli P with computational displacement a,
$`D_3PD_3`$ is P times the Boolean derivative sign
$`(-1)^{g(z+a)+g(z)}`$. Its exponent q is quadratic and satisfies
$`q(z+a)=q(z)`$. A quadratic Boolean sign has Walsh support an affine
space of size $`2^{2h}=4^h`$, with every nonzero coefficient of magnitude
$`2^{-h}`$ for some integer h at least zero.

For completeness, let B be q's alternating polar form. Squaring its
normalized Walsh coefficient gives

```math
\widehat{(-1)^q}(s)^2
 =2^{-w}\sum_{t\in\ker B}(-1)^{q(t)+q(0)+s\cdot t}.
                                                               \qquad(14)
```

The restriction $`q(t)+q(0)`$ is linear on $`\ker B`$. Thus the
sum is nonzero on one affine coset of $`(\ker B)^\perp`$, where its
value is $`2^{-\mathrm{rank}B}`$. Alternating binary forms
have even rank $`2h`$, proving the claim. Invariance under a additionally
puts every supported Fourier label s in $`a^\perp`$. Consequently
all Pauli terms $`PZ^s`$ commute, including their literal Hermitian
phases. Their binary labels form the same affine-size support.

It follows that, for any Hermitian Pauli P, the operator
$`C_2^\dagger D_3C_3^\dagger PC_3D_3C_2`$ has a commuting affine
support F of size $`4^h`$, with coefficient magnitudes $`2^{-h}`$.
This statement includes arbitrary mixing by both Clifford walls.

### Every resulting Pauli coefficient is at most one half

Normalized trace against P expresses its coefficient in B_H as the
Pauli inner product of R_H with the preceding operator. Since F is
commuting, it cannot intersect both E and O. If it intersects O, the
absolute inner product is at most

```math
2^{-h}\frac{\min(4^h,8)}4\le\frac12.                  \qquad(15)
```

For E, h=0 gives at most 1/2. If $`h=1,2`$ and $`e\in F`$,
the intersection $`F\cap L`$ is a proper affine subspace of F,
because $`e\notin L`$. It has size at most $`4^h/2`$, and hence

```math
2^{-h}\left(\frac12+\frac{4^h}{16}\right)=\frac38.
                                                               \qquad(16)
```

If $`e\notin F`$, the bound instead is
$`2^{-h}\min(4^h,16)/8\le1/2`$. For $`h\ge3`$, even the
entire E support gives at most $`(1/2+16/8)2^{-h}\le5/16`$.
These cases prove the uniform bound

```math
\max_P
 \left|2^{-w}\mathrm{Tr}(PB_H)\right|\le1/2.    \qquad(17)
```

The coefficient ceiling itself is sharp: D_3=I retains the coefficient
1/2 of e, and a Clifford can move e to any specified nonidentity Pauli.
No claim is made that the complete error bound (12) is sharp.

For the target, averaging over the independent c_1,c_2,d inputs in (1)
leaves only addresses with $`a=b=0`$. Therefore

```math
2^{-w}\mathrm{Tr}\bigl(Z_r(Q_\square\otimes I)\bigr)
 =\Pr[a=b=0]=(3/4)^2=9/16.                            \qquad(18)
```

The normalized Pauli trace of an operator is bounded by its operator
norm. Equations (17)--(18) prove (12), on every input and after tensoring
any external reference. Extra arbitrary wires cannot dilute either
normalized coefficient.

### Weighted outer phases leave this screen

The cubic-sign restriction is essential. Set

```math
C_2=(HSH)_u(HSH)_v,\qquad D_3=T_uT_v,\qquad C_3=I.
```

With literal adjoints,
$`(HSH)^\dagger T^\dagger X T(HSH)=H`$. Thus, for
$`P=Z_rX_uX_v`$, its pullback through the two later factors is K_H.
Equation (13) now gives the actual coefficient

```math
2^{-w}\mathrm{Tr}(PB_H)
 =\frac12\left(\frac12+\frac14+\frac14+\frac18\right)
 =9/16.                                              \qquad(19)
```

These are allowed weighted third-level outer phases, but they are not
Boolean-sign phases. Their Pauli images can have noncommuting support,
so (17) cannot be extended to them. The complete action of this particular
word is nevertheless still wrong:

```math
B_H=Z_r[(1-a)X_u-aY_u][(1-b)X_v-bY_v].                \qquad(20)
```

Each helper axis is X at its zero address product and minus Y at its
one product. The coefficients at $`Z_rX_uX_v`$ and $`Z_rY_uX_v`$
are respectively $`9/16`$ and $`-3/16`$, since the two address pairs
are independent. These Paulis anticommute. As in Section 3, at least one
remains absent from any computational-diagonal target after every final
Clifford F, including arbitrary additional work. Consequently this
specific weighted word satisfies
$`\|FB_HF^\dagger-Q_\square\otimes I\|\ge3/16`$.

Thus (19) escapes one screen, while (20) diagnoses its complete residual
action. Neither gives lookup or an all-size source. For sign-only
discovery, changing K_H or its first coupling is necessary. General
weighted discovery remains open and requires a changed full operator.

The [exact nonmonomial fixtures](../../tests/test_nonmonomial_source_screen.py)
check (13) on every active column, its two support classes, quadratic
Walsh support and commutation, the target coefficient (18), and the
literal weighted-phase escape (19). The unrestricted later-Clifford
and all-helper scope of (12) is analytic, not inferred from finite tests.

<a id="fourier-swap-fixed-coefficient-face"></a>

## 6. A Fourier-dressed swap completes the fixed-coefficient face

A changed nonmonomial Clifford gives a positive complete identity. It
implements the seven-wire target

```math
Q_1=(-1)^{r+a v+b u+ab},\qquad a=x_1x_2,\quad b=x_3x_4. \qquad(21)
```

Here u,v are arbitrary logical coefficient bits, identified with
$`c_2,c_1`$ respectively. Equation (21) is the d=1 restriction of (1),
stated as a separate complete seven-wire operator. No initialized d is
supplied, and tensoring (21) with identity on an arbitrary d does not
implement (1).

Let $`B=\mathrm{CX}(u;v)H_u`$ and define

```math
\begin{aligned}
K_F&=Z_r H_v\mathrm{SWAP}(u,v)H_v,\\
U_F&=H_vB\,\mathrm{CCX}(u,v;r),\\
D_2&=\mathrm{CCZ}(u,x_1,x_2)\mathrm{CCZ}(v,x_3,x_4),\\
C_2&=H_rB^\dagger H_v,\qquad
D_3=\mathrm{CCZ}(r,u,v),\qquad C_3=H_r,\\
V_F&=C_3D_3C_2D_2U_F.
\end{aligned}                                                   \qquad(22)
```

K_F is a traceless Hermitian Clifford, genuinely nonmonomial on u,v.
Its preserved computational Pauli intersection on that pair is trivial:
it sends $`Z_u,Z_v`$ to $`X_v,X_u`$. The first source satisfies
$`U_F Z_rU_F^\dagger=K_F`$ because B maps CZ to SWAP.

Put $`R_a=Z_uX_v`$ and $`R_b=X_uZ_v`$. They commute with each other
and with the dressed swap. Direct Pauli multiplication gives the literal
cross sign

```math
D_2K_FD_2=(-1)^{ab}R_a^aR_b^bK_F.                     \qquad(23)
```

The Clifford $`B^\dagger H_v`$ maps $`R_a,R_b`$ to $`Z_v,Z_u`$
and the dressed swap to CZ. Thus the observable after C_2 is
$`X_r\mathrm{CZ}(u,v)(-1)^{a v+b u+ab}`$. Conjugation by D_3 adds
another CZ to the X_r factor and cancels that background exactly; the
last H_r restores Z_r. Consequently

```math
V_F Z_rV_F^\dagger=Q_1\otimes I_{\rm work}.             \qquad(24)
```

There is no residual carrier translation or helper phase. This success
changes the earlier coupling and therefore does not contradict Section 3.

### Literal native cost

The one-pair shared-root CCZ compiler inside U_F costs seven T gates in
four bands and ten Cliffords. Its two root Hadamards, B and H_v add five
Cliffords, giving $`G_{\rm all}(U_F)=22`$. The two returned-dirty cubic batches
in (22) cost respectively 16 and eight T gates, eight bands each, and
32 and 16 total gates. C_2 and C_3 add four and one Cliffords. The
same bounds hold for the actual inverse:

```math
T_{V_F}=31,\qquad D_T(V_F)\le20,\qquad G_{\rm all}(V_F)=75,
\qquad w_{\rm extra}=2,\qquad a_{\rm clean}=0.          \qquad(25)
```

The two extra wires are arbitrary Gray masks returned after each batch;
one is reused by D_3. The literal phase word $`V_F^\dagger,Z_r,V_F`$
therefore costs 62 T gates, at most 40 bands and 152 total gates,
including two S gates for Z_r. Its error is zero, with every logical,
mask and reference input allowed. This finite count is not a comparison
with the full variable-coefficient face or an all-size lookup bound.

The word also has an ordinary classical normal form. Put
$`M=\mathrm{CCX}(u,v;r)`$ and
$`T_{ab}=\mathrm{CCX}(x_1,x_2;u)\mathrm{CCX}(x_3,x_4;v)`$.
Then the same logical source and completed phase satisfy

```math
V_F=M T_{ab} M,\qquad
Q_1=Z_r T_{ab}\mathrm{CZ}(u,v)T_{ab}\mathrm{CZ}(u,v).
```

Indeed $`B^\dagger H_vD_2H_vB=X_u^aX_v^b`$. The source translates
u,v by a,b and r by $`av+bu+ab`$; its completed probe returns every
carrier. Each helper-free Toffoli costs seven T gates, four bands and
19 total gates. The middle two occupy disjoint triples, so the direct
source costs 28 T gates, at most 12 bands and $`G_{\rm all}=76`$, with no
extra wires. The direct phase costs 28 T gates, at most eight bands and
$`G_{\rm all}=84`$, also with no extra wires. These are literal finite
simplifications of (24), not a new all-size mechanism: the successful
nonmonomial pilot reduces to three ordinary Toffoli stages and a
quadratic finite difference.

<a id="common-readout-four-corner-gap"></a>

### Why a common later readout cannot supply arbitrary d

The direct extension of (22) has a uniform obstruction, even when the
common later readout is an arbitrary unitary. Let P,R be commuting
Hermitian Paulis, K any Hermitian Clifford, and

```math
O_{ab}=P^aR^b K R^bP^a,\qquad a,b\in\{0,1\}.
```

All other registers, including d and any work, remain arbitrary. The
same unitary W must complete all four address fibers. It may mix those
other registers freely, but it holds the address and is independent of
a,b. Since $`KPK`$ is a Hermitian Pauli,

```math
O_{00}O_{10}O_{11}O_{01}
 =(KPK)R(KPK)R\in\{+I,-I\}.                            \qquad(26)
```

The same four-corner product after W is still that scalar. In contrast,
the desired fibers $`Q_{ab}=(-1)^{r+a c_1+b c_2+ab d}`$ obey
$`Q_{00}Q_{10}Q_{11}Q_{01}=Z_d`$. Their product has distance two from
either scalar on the full arbitrary-d space. Telescoping four unitary
factors proves

```math
\max_{a,b}\|W O_{ab}W^\dagger-Q_{ab}\otimes I\|\ge1/2.  \qquad(27)
```

For the concrete address products in (21), the four fibers can be taken
at $`x_2=x_4=1`$, with $`x_1=a,x_3=b`$. The theorem applies only when
the later word is the same held-address operation on those four fibers.
It includes arbitrary later weighted phases on the other registers and
arbitrary helper width. Later address-dependent operations, genuine
address motion, or a changed first coupling need a fresh analysis.
The scalar minus sign allowed by (26) is exactly what makes the d=1
positive identity possible; it cannot encode an arbitrary coefficient.

The [Fourier-swap fixtures](../../tests/test_fourier_swap_face.py) check
every seven-wire column, the nonmonomial Clifford and its first source,
both actual native cubic inverses on arbitrary masks, the common ledger,
and examples of both signs in (26). Equations (23)--(27), rather than
those finite checks, establish the complete identities and scoped gap.
The next construction must change the later address participation or
leave this two-Pauli coupling grammar; repeating the fixed-coefficient
pilot does not supply the missing uniform source.
