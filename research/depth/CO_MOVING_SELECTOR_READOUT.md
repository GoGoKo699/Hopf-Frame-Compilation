# Paid selector tags and a co-moving readout obstruction

[Working checkpoint](../../WORKSPACE.md) ·
[Co-moving stage](CO_MOVING_LOOKUP_STAGE.md) ·
[Complete recovery interface](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Address-relative selector Hadamards have a constant native-depth
implementation, including their clean-controlled versions and actual
inverses. They genuinely mix selectors, but the resulting one-echo
readout still fails on the complete bank space. The exact echo has a
sixth-root collapse, and a stronger two-row singlet obstruction survives
arbitrary phase angles and identical local rotations on co-moving row
pairs. These are paid components and scoped negative results, not a
complete lookup source or a new general depth theorem.

## 1. Conventions and the complete target

Use the notation of the [stage proof](CO_MOVING_LOOKUP_STAGE.md#1-a-complete-charge-preserving-stage).
There are p address bits, $`m=2^p`$ bank bits and m arbitrary selectors.
C is a nonempty matching of q disjoint address CNOTs, with control set S.
All bank, selector and private dirty inputs may be entangled with a
reference. Define

```math
\begin{aligned}
L|x,Y,s\rangle&=|Cx,CY,Cs\rangle,&
h&=\sum_uY_us_u\pmod2,\\
\beta_u(x)&=(x_c\oplus u_c)_{c\in S},&
d_u(x)&=|\beta_u(x)|,\\
n&=\sum_uY_us_ud_u(x),&N&=\sum_uY_ud_u(x),\\
R&=i^nL^h,&\Omega|x,Y,s\rangle&=(-1)^{Y_x}|x,Y,s\rangle.
\end{aligned}                                                    \qquad(1)
```

The simultaneous relabeling preserves h, n and N. It also preserves each
control mismatch, in the sense $`\beta_{Cu}(Cx)=\beta_u(x)`$.
The already paid stage R and its actual inverse have common bounds
$`T_R,G_R,d_R`$ given in that proof, with $`d_R\le116`$,
at most 3qm returned private dirty helpers, no clean input and zero error.
The m selectors are principal arbitrary inputs, not initialized helpers.

The required accepted Gram is

```math
H_Y=(I+\Omega)/4=\tfrac12\Pi_+,
\qquad \Pi_\pm=(I\pm\Omega)/2.                              \qquad(2)
```

Thus a marked address has target zero and an unmarked address has target
one half. Every error below is the complete operator-norm Gram error,
not a computational-basis diagonal check.

## 2. Paid affine address-relative Hadamards

For a nonzero binary linear form ell of support size k and a constant
$`\alpha\in\{0,1\}`$, define the selector unitary

```math
B_{\ell,\alpha}=\prod_u H_{s_u}^{\ell(x\oplus u)\oplus\alpha}.
                                                                    \qquad(3)
```

This is not treated as a free address-controlled gate. Let

```math
V=SHTHS^\dagger,\qquad VZV^\dagger=H.                       \qquad(4)
```

The equality includes the literal scalar phase. Emit controlled-H on
control z and selector s chronologically as

```text
V†(s), H(s), CX(z,s), H(s), V(s).
```

Expand each $`S^\dagger`$ as three S gates and reduce adjacent S powers
modulo four, including when reversing the actual native word. The word
and its actual inverse each cost two T-type gates and 15 Cliffords. Only the
selector receives T gates. Across all rows, the first V-adjoint T band,
the shared-control Clifford CZ gates, and the final V T band give
T-depth two; no concurrency of gates sharing z is assumed inside a
Clifford band.

Compute $`z=\ell(x)`$ into an existing address pivot using k-1 CNOTs.
For row u add a fixed H if $`\ell(u)\oplus\alpha=1`$, and undo the
coordinate change afterwards. Exactly half the rows need this H. Hence
both (3) and its actual native inverse have

```math
T_B\le2m,\qquad G_B\le31m/2+2(k-1),\qquad d_B\le2.          \qquad(5)
```

There are no helpers, clean inputs or implementation error. The zero
linear form is simply identity or a global selector Hadamard. Adjacent
tags of the same H type collapse exactly:

```math
B_{\ell,\alpha}B_{r,\gamma}
 =B_{\ell\oplus r,\alpha\oplus\gamma}.                     \qquad(6)
```

Many uninterrupted affine labels therefore do not create additional
selector dynamics. The implementation (5) allows every affine label;
the sector identities below require the specified covariance under C.

## 3. Paying the signal control and the whole difference source

Let a be an arbitrary principal control. Controlling (3) requires the
row action $`H_s^{a(z\oplus j_u)}`$, with
$`j_u=\ell(u)\oplus\alpha`$. Between V-adjoint and V, implement
$`\mathrm{CCZ}(a,z,s)`$ and, when $`j_u=1`$, $`\mathrm{CZ}(a,s)`$.
To keep T gates off both shared controls, give each row two private
arbitrary dirty bits $`b_1,b_2`$ and use the chronology

```text
CCZ(b1,b2,s), CX(a,b1), CCZ(b1,b2,s), CX(z,b2),
CCZ(b1,b2,s), CX(a,b1), CCZ(b1,b2,s), CX(z,b2).
```

Both dirty bits return exactly. The accumulated binary phase exponent is

```math
b_1b_2s+(b_1\oplus a)b_2s
 +(b_1\oplus a)(b_2\oplus z)s+b_1(b_2\oplus z)s=azs.         \qquad(7)
```

Each CCZ uses the ordinary exact seven-T, six-CNOT parity decomposition.
The four corners cost 28 T gates and 28 Cliffords, including the four
dirty-bit flips. With V-adjoint and V this becomes 30 T and 40 Cliffords;
the optional CZ costs three more Cliffords. The triples $`(b_1,b_2,s)`$
are disjoint across rows, so four sequential seven-layer CCZ bands and
the two V bands give a legitimate depth-30 schedule. The shared a and z
receive only Clifford gates. Thus

```math
T_{\mathrm cB}\le30m,\qquad
G_{\mathrm cB}\le83m/2+2(k-1),\qquad
d_{\mathrm cB}\le30,\qquad b_{\mathrm cB}=2m.                \qquad(8)
```

This is a common bound for the actual forward and inverse native words,
with zero clean requirements and zero error. Counting a as a clean
signal is a separate choice made by the following source.

Initialize a to zero, accept a equal to zero, and emit

```text
H(a), controlled-B(a), R, controlled-B†(a), Z(a), H(a).
```

No controlled R is used. Its exact accepted block is

```math
K=(R-B^\dagger RB)/2=R(I-E)/2,\qquad
E=R^\dagger B^\dagger RB,\qquad
K^\dagger K=(I-\mathrm{Re}E)/2.                            \qquad(9)
```

The accepted fiber may move address, bank and selectors. With
$`Z=S^2`$, one simultaneous source/actual-adjoint ledger is

```math
\begin{aligned}
T_U&\le T_R+60m,&G_U&\le G_R+83m+4(k-1)+4,\\
D_T(U)&\le176,&b_{\rm private}&\le\max\{3qm,2m\}=3qm.
\end{aligned}                                                   \qquad(10)
```

There is one clean signal and zero implementation error. All private
work returns between components, so its pool is reusable; the m arbitrary
selectors remain separate. This pays a normalized source, not the Gram
identity (2), recovery, or a lookup implementation.

## 4. Complete echo and the sixth-root collapse

Choose $`\ell(x)=\bigoplus_{c\in S}x_c`$, $`\alpha=0`$, and call
the tag $`B_S`$. Its row exponent is $`a_u=d_u\pmod2`$.
It commutes with L. Put $`P_\sigma=(I+\sigma L)/2`$ for
$`\sigma\in\{+1,-1\}`$ and

```math
g_\eta=S^{-\eta}HS^\eta H,\qquad \eta\in\{+1,-1\}.
                                                                    \qquad(11)
```

On the L-eigenspace sigma, the stage R is a product of selector phases
$`S^{(d_u+2\epsilon_\sigma)Y_u}`$, where
$`\epsilon_\sigma=(1-\sigma)/2`$. Matched even-distance factors are
untouched by the tag. For odd d the exponent modulo four is
$`\eta_{\sigma,d}=\sigma(-1)^{(d-1)/2}`$. Therefore the full-space
relative operator, on arbitrary coherent banks, is

```math
E_H=R^\dagger B_SR B_S
 =\sum_{\sigma=\pm1}P_\sigma
   \prod_{u:\ d_u\ {m odd}}
    g_{\sigma(-1)^{(d_u-1)/2}}(s_u)^{Y_u}.                  \qquad(12)
```

The products in (12) are diagonal in address and bank, not in selectors.
They are covariant under simultaneous relabeling and hence commute with
the projectors $`P_\sigma`$. This is why multiplying or cubing the
sector expression does not discard cross-sector terms improperly.
Direct one-qubit algebra gives

```math
\det g_\eta=1,\qquad \mathrm{tr}\,g_\eta=1,\qquad
g_\eta^2-g_\eta+I=0,\qquad g_\eta^3=-I.                     \qquad(13)
```

Consequently, on the complete space,

```math
E_H^3=(-1)^N I,\qquad E_H^6=I.                              \qquad(14)
```

The cube is only a low-moment Clifford phase. Writing
$`w=\sum_uY_u`$ as an integer, its binary exponent is

```math
N\equiv w\sum_{c\in S}x_c
 +\sum_uY_u\sum_{c\in S}u_c\pmod2.                         \qquad(15)
```

Computing the address parity in a pivot, applying CZ
between it and every bank bit, applying the indicated fixed bank Z mask,
and undoing the pivot change realizes (15) with no T gates.

For a single CNOT and a singleton bank $`Y=e_u`$, (12) is identity
when $`x_c=u_c`$ and has $`\mathrm{Re}E_H=I/2`$ when
$`x_c\ne u_c`$. Thus (9) has Gram zero or one quarter according to
one control mismatch. This exact singleton contrast does not establish
the full-address target and does not survive arbitrary banks as a gap.

## 5. A two-row singlet shield for arbitrary angles and local tags

The following obstruction is not special to the roots in (14). Allow
any real angle theta and any one-qubit unitary $`V_\beta`$ for each
control-mismatch vector beta, and define

```math
R_\theta=e^{i\theta n}L^h,\qquad
B_V=\prod_u V_{\beta_u(x)}(s_u),\qquad
E_V=R_\theta^\dagger B_V^\dagger R_\theta B_V.               \qquad(16)
```

This is an algebraic extension, not a native synthesis claim for an
arbitrary angle or unitary. The tag still commutes with L. With

```math
Q_{\sigma,d}=\mathrm{diag}(1,\sigma e^{i\theta d}),\qquad
g_{\sigma,\beta}=Q_{\sigma,|\beta|}^\dagger
 V_\beta^\dagger Q_{\sigma,|\beta|}V_\beta,
```

the same sector calculation gives

```math
E_V=\sum_{\sigma=\pm1}P_\sigma
       \prod_u g_{\sigma,\beta_u(x)}(s_u)^{Y_u},\qquad
\det g_{\sigma,\beta}=1.                                   \qquad(17)
```

Choose any moved row u and let $`v=Cu\ne u`$. Such a row exists
because the matching is nonempty. Fix the same two-row bank
$`Y=e_u+e_v`$ and put its two active selectors in the singlet

```math
|\psi^-\rangle_{u,v}=(|01\rangle-|10\rangle)/\sqrt2.         \qquad(18)
```

Every matching control coordinate is unchanged by C, so
$`\beta_u(x)=\beta_v(x)`$ for every address x. The product in each
sector of (17) acts on the active pair as g tensor g and acts identically
on inactive selectors. For any two-by-two matrix g,

```math
(g\otimes g)|\psi^-\rangle=(\det g)|\psi^-\rangle.           \qquad(19)
```

Define the complete input and its exact echo action by

```math
\begin{aligned}
|\Psi_{x,\chi}\rangle
 &=|x,e_u+e_v\rangle|\psi^-\rangle_{u,v}|\chi\rangle_{\rm inactive},\\
E_V|\Psi_{x,\chi}\rangle&=|\Psi_{x,\chi}\rangle.
\end{aligned}
                                                                    \qquad(20)
```

The inactive state is arbitrary, including entanglement with a reference.
No L-eigenstate assumption is needed: each sector product fixes the
input, and $`P_++P_-=I`$. In particular (20) holds at both marked
addresses u,v and at every unmarked address, which exists for $`p\ge2`$.
If $`V_0=I`$, the echo at the marked addresses is additionally identity
on all selector states; the singlet collision does not need this premise.

For the balanced difference source in (9), the unmarked state (20) has
accepted norm zero instead of squared norm one half. Its complete Gram
error is therefore at least one half. More generally, any common
signal/auxiliary preparation, interleaving and acceptance using only
$`E_V,E_V^\dagger`$ and acting trivially on address, bank and selectors
between calls has the same accepted work action on these marked and
unmarked eigenvalue-one inputs. The target expectations differ by one
half, so its worst-case Gram error is at least one quarter. This includes
scalar Laurent-polynomial filtering and signal operations depending
only on bank weight: both witnesses use the same weight-two bank.
Arbitrary bank-dependent operations are not included, and this is not
a bound on arbitrary interleaved co-moving sources.

## 6. Boundary of the result and the next complete operation

The paid H tags escape the selector-zero premise of unmixed stages,
but replace it with the complete orbit-pair collision above. Changing
only the phase angle or a common local mismatch rotation cannot remove
that collision. Applying many adjacent same-H affine tags only invokes
(6), and filtering the completed echo only invokes the one-quarter floor.

Unequal local labels within a C-pair, entangling selector operations, or
genuinely interleaved different matchings are outside the hypotheses.
They are possible mechanism changes, not established successful sources.
Their first required identity is still a complete accepted block K with
$`K^\dagger K\approx H_Y`$, retaining every bank and selector input,
marked-sector cancellation, off-diagonal terms, actual inverses and all
signal controls. No coherent row selection or controlled stage is free.

[Exact selector fixtures](../../tests/test_co_moving_selector_readout.py)
check (4), the controlled-H word, the two-dirty
four-corner construction and their actual reverses on every small input
column. Logical fixtures check (12) on every two-bit address/bank/selector
column in both matching orientations, (13), and (20) with H and a
noninvolutive SH tag, arbitrary inactive-selector basis states, and a
two-control matching example. These finite checks corroborate the
identities; the all-size proofs and common resource bounds are above.

## 7. Affine row motion does not repair the unchanged four-path word

A separate relative-coordinate screen applies even if physical row zero
moves. This section is an exact logical identity, not a native ledger for
affine stages. Let $`g(x)=Ax\oplus b`$, where A is a disjoint CNOT
matching and $`Ab=b`$, so g is an involution. Let S be a subset of coordinates
unchanged by A, and define

```math
R_{A,b,S}=i^{\sum_us_uY_ud_S(x,u)}L_g^h,\qquad
d_S(x,u)=\sum_{c\in S}(x_c\oplus u_c).
```

Write $`T_v`$ for selector-row translation $`u\mapsto u\oplus v`$.
The map from relative address z and arbitrary selector bank q,

```math
J|z,q\rangle=m^{-1/2}\sum_v|z\oplus v,e_v,T_vq\rangle,       \qquad(21)
```

is an isometry: different v have orthogonal bank labels, and for each v
the other two changes of coordinates are permutations. On every summand,
$`h=q_0`$ and $`n=q_0d_S(z,0)`$. If $`h=1`$, set
$`v'=Av\oplus b`$. The output address is $`Az\oplus v'`$,
the output bank is $`e_{v'}`$, and the selector bank is $`T_{v'}Aq`$.
Reindexing the uniform sum gives the literal intertwining identity

```math
R_{A,b,S}J=J\widehat R_{A,S},\qquad
\widehat R_{A,S}|z,q\rangle
 =i^{q_0d_S(z,0)}|A^{q_0}z,A^{q_0}q\rangle.                 \qquad(22)
```

The reduced stage is independent of b. Both sides are unitary on the
invariant range of J, so the actual adjoints intertwine as well.
A selector-only M commuting with every $`T_v`$ obeys $`MJ=JM`$;
this includes global selector X and H, their signal-controlled forms,
and arbitrary signal-only interlayers. A bank-weight phase $`f(w)`$
reduces to the scalar $`f(1)`$. Consequently every word in these
operations and the affine stages restricts to its offset-free word.
The target also intertwines:

```math
H_YJ=J\left[\tfrac12(I-|0\rangle\langle0|_z)\otimes I_q\right].
                                                                    \qquad(23)
```

The [two-bit four-path marked witness](CO_MOVING_LOOKUP_STAGE.md#4-the-earlier-four-path-acceptance-also-fails-after-transplantation)
therefore lifts through J. Set $`z=0`$, put selectors (1,2,3) in
$`(|100\rangle-|001\rangle)/\sqrt2`$, and leave selector zero
arbitrary. Its acceptance is still three quarters for every admissible
pair of affine offsets, although the target is zero. The coherent bank
input in (21) is allowed by the complete contract.

This excludes repairing that unchanged word only by adding affine
offsets. Even relative-row product tags $`\prod_uV_{f(x\oplus u)}(s_u)`$
intertwine J with the same formula on (z,q), whether or not they commute
with A; mere row-nonuniformity need not break this screen. Fixed physical
row labels or address-only interlayers may break it and need separate
analysis. Changed acceptance remains outside the unchanged four-path
failure even when J is preserved.
The [exact affine fixture](../../tests/test_affine_singleton_transport.py)
checks (21)-(23), the allowed global interlayers, actual adjoints and all
two-bit offset-pair marked witnesses. Removing a physical fixed row is
therefore not by itself a complete readout mechanism.
