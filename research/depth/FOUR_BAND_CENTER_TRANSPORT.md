# Four-band center transport and the first unmatched native axis

[Original physical-pair scout](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md) ·
[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md) ·
[Actual native matching](ONE_BAND_AXIS_MATCHING.md)

These results revise the first unpaid inference in the fixed-endpoint
four-band scout. They concern actual physical pairs and independent
native lists. None is an unrestricted four-layer Boolean-degree theorem
or a change to the complete approximate dirty-indicator frontier.

| Result | Exact scope and consequence |
| --- | --- |
| Section A: one pair-orientation mismatch | With the same original endpoint probes, one original B pair and the preserved center, D has the same signed J action; Section F then gives paid two-band resynthesis |
| Section B: lost initial matching | A nontrivial complete Boolean controlled-SWAP source has a non-Clifford intermediate probe, while its paid two-band, eight-T realization remains available |
| Section C: one anticommuting unmatched orbit | With the fixed endpoints and preserved center, exact return forces an actual inverse rotation and removes two T gates |
| Section D: automatic center preservation | Complete occupied Booleanity forces every D axis to preserve the original center for these fixed endpoints, without complete B matching |
| Section E: a general original-pair derivative | Arbitrary original endpoints are allowed; the center's logical X direction has a complete two-band Boolean derivative of degree at most three |
| Section F: paid whole-word resynthesis | With an all-active common endpoint and one support-two matched pair per middle band, the full word has an actual two-band realization with its Clifford cost paid |

Section D therefore supplies the center premise used in A and C when
the complete Boolean contract holds. The original equal endpoints and
their computational center transport remain essential hypotheses.

## Fixed physical data

On wires x1,x2,z,d,w, put

```math
a_1=Y_1X_2,\quad a_2=X_1Y_2,\quad r=a_1a_2=Z_1Z_2,
\qquad p=X_1X_2,\quad h=Z_d.
```

Use the actual endpoints

```math
E_A=E_C=T_{a_1}T_{a_2}^{\dagger}T_{Z_w},\qquad
K_A=K_C=H_w\,CNOT_{w\to1},
```

with the chronological order for K being CNOT then H. Their clean
reflection is

```math
J=J_{\rm core}J_w,\qquad
J_{\rm core}=SWAP_{12}CZ_{12},\qquad
J_w=(X_w+Y_w)/\sqrt2.
```

Additional logical and dirty wires are allowed. These endpoint words
act identically on them; the actual middle lists may act on every
wire. The explicit fixture below needs only the five displayed wires.

The complete source is

```math
V=K^\dagger E_A^\dagger E_D E_B E_AK.
```

There is one initialized clean wire x1, arbitrary logical inputs and
arbitrary returned dirty inputs/references. Complete occupied Boolean
diagonality retains one literal common scalar. The unused clean block
is unrestricted unless a full-column identity is proved. Clean return
alone gives the full reflection identity

```math
E_D(E_BJE_B^\dagger)E_D^\dagger=J.
```

Each native list consists of independent commuting signed Hermitian
Pauli axes, each with exponent plus or minus one. Signed scalar phases
are part of the actual gates. Except where explicitly relaxed below,
every actual D axis commutes with the original physical center r.

## A. One B pair-orientation mismatch also restores common J action

Assume B has exactly one genuine original pair under J, of difference
r, and any independent negative-singleton list A. All other orbits are
positive singletons. Suppose the genuine pair is orientation-mismatched
relative to its literal J-aligned orientation. The aligned case is
already covered by [negative-span restoration](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md#15-negative-singleton-restoration-and-an-actual-gate-removal).
Let b be the actual
pair axis used in its S repair; its partner has binary label b+r.
The original endpoint probes, before absorbing that repair, are both J.

Let $`S=E_BJE_B^\dagger`$ and write its binary symplectic action as I+N.
Let Omega be the symplectic form. The literal one-pair repair and
negative-singleton update give

```math
\begin{aligned}
N_J-N&=(br^{\mathsf T}+rb^{\mathsf T}+AA^{\mathsf T})\Omega,\\
Nb&=r,\qquad Nr=NA=0.
\end{aligned}
```

S is Clifford, so matching of the actual D axes under S is justified
by the actual one-band theorem. Let C contain one chosen repair axis
from each mismatched D pair, R=NC its pair differences, and L the
negative D singleton columns. Literal D repair gives

```math
N_J-N=(CR^{\mathsf T}+RC^{\mathsf T}+LL^{\mathsf T})\Omega.
```

Left multiplication by N yields

```math
RR^{\mathsf T}=rr^{\mathsf T}.
```

Actual D independence makes R full column rank. Consequently there is
exactly one mismatched D pair, its difference is r, and its chosen
repair axis c obeys Nc=r. Put $`h_0=b+c`$. The remaining equation is

```math
LL^{\mathsf T}=AA^{\mathsf T}
                         +h_0r^{\mathsf T}+rh_0^{\mathsf T}.
```

Native D independence forbids r in span(L), because the genuine D pair
already spans r. The columns of [A,r] are independent by original B
independence. If $`h_0`$ were outside span(A,r), the right side would
factor through the independent columns $`[A,r,h_0]`$ with a nonsingular
middle matrix: identity on A and the exchange matrix on the last two
coordinates. Its image would contain r, a contradiction. Thus
$`h_0=At+\gamma r`$. In the basis [A,r], the remaining middle matrix is

```math
\begin{pmatrix}I&t\\t^{\mathsf T}&0\end{pmatrix}.
```

Its image contains r exactly when $`t^{\mathsf T}t=1`$. Therefore t has
even parity and

```math
\mathrm{span}(L)=\{Ax+r(t^{\mathsf T}x):x\in\mathbf F_2^{|A|}\}.
```

In particular At lies in span(L), so
$`b=c+At+\gamma r`$ lies in the actual D span, as does every column of A.
Every D axis therefore commutes with the actual B repair axis b,
the pair center r, and all B negative-singleton axes.

The image under J of a D axis also commutes with these axes: J sends
b to its signed partner b+r and negates every negative singleton.
Hence conjugation by the B repair and negative rotations does not
change that image. S and J have the same literal signed action on
every actual D axis.

The actual D mismatched pair already supplies one r-pair. For these
fixed endpoints, any J-pair in a native D list preserving r must have
difference r: write $`N_Jd=t_1a_1+t_2a_2+t_3Z_w`$. Center preservation
gives $`t_1=t_2`$; commutation of the two actual pair axes gives
$`0=\langle d,N_Jd\rangle=t_1+t_2+t_3`$, hence $`t_3=0`$.
Native independence permits at most
one such pair. Thus the original, unmodified fixed-endpoint source
falls under the retained
[cubic theorem](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md#14-native-cubic-closure-of-the-fixed-endpoint-one-pair-class).
Section F below strengthens this branch to an actual two-band
realization with its Clifford cost paid.

This does not assert a free same-budget replacement, nor a theorem for
genuinely unequal original endpoint probes. In the special case A is
empty, $`h_0`$ is zero or r, so the actual B/D pair axes coincide and the
[shared-pair removal](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md#8-remove-common-actual-pairs-before-studying-the-remaining-interaction)
applies literally.

## B. A nontrivial Boolean source with a genuinely non-Clifford middle probe

The following repartition retains all eight literal commuting factors
of the [controlled-SWAP example](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md#11-a-nondiagonal-middle-defeats-endpoint-commutation-but-retains-two-bands):

| Band | Literal signed Hermitian axes | T exponents |
| --- | --- | --- |
| B | hp, hpr, hpZ_z, h | +1, -1, -1, -1 |
| D | hprZ_z, hZ_z, hr, hrZ_z | +1, +1, +1, -1 |

Here pr is the literal product XX times ZZ, namely -YY. Those minus
signs must not be silently dropped from the axes.

In the diagonal parity coordinates (x1,x2,z,d), the row masks are

```math
B=[9,11,13,8],\qquad D=[15,12,10,14].
```

Both matrices have rank four. Subtract B's row 8 to expose x1,x2,z.
For D, 14 xor 12 gives x2, 14 xor 10 gives z, 15 xor 14 gives x1,
and then d is recovered. Thus these are actual independent native
bands, with no clean parity workspace.

A concrete B parity encoder is the chronological sequence
CX(d->x1), CX(x1->x2), CX(x1->z). A concrete D parity encoder first
adds x2,z,d into x1, then swaps x2 and z, then applies
CX(d->x2), CX(d->z), CX(x2->d), CX(z->d). The swap uses three CNOTs,
so these encoders use 3 and 10 CNOTs, respectively. Apply the actual
inverse encoder after each native diagonal band. Use the same fixed
Q from Section 11 (chronological CX(x1->x2), H(x1)) to carry the
diagonal axes to the displayed physical axes. This gives a fully
paid five-wire native implementation; the middle encoders and their
inverses cost 30 elementary H/CNOT gates after canceling the adjacent
Q-inverse/Q boundary.

The B list retains the genuine negative signed original pair
$`hp\leftrightarrow-hpr`$. But J sends $`hpZ_z`$ to $`-hprZ_z`$,
outside the B list. Therefore $`E_BJE_B^\dagger`$ is not Clifford by
the exact one-band
matching theorem. This is not a logical-pair substitution.

A direct literal certificate is also available. Put $`A_0=hpZ_z`$ and
$`A_1=hprZ_z`$, so $`JA_0J=-A_1`$. The remaining aligned pair factor
and positive singleton factor commute with J and both axes, giving

```math
E_BJE_B^\dagger=T_{A_0}^{\dagger}JT_{A_0}
                    =\omega T_{A_0}^{\dagger}T_{A_1}^{\dagger}J.
```

$`A_0`$ and $`A_1`$ are commuting independent Pauli axes. Their two-T
band is not Clifford: a Pauli anticommuting only $`A_0`$ acquires two nonzero Pauli
terms. Multiplication by Clifford J cannot change that conclusion.

Nevertheless, reordering commuting actual factors proves the complete
operator identity

```math
E_DE_B=CSWAP_{z;1,2}\otimes I_{d,w}.
```

All literal scalar phases are identical to that example. Consequently
the full physical source is

```math
V=(-1)^{(x_1\mathbin\oplus w)(1+x_2)z}.
```

Its occupied action is the nontrivial Boolean phase
$`(-1)^{w(1+x_2)z}`$, with exact dirty and clean return on every reference.
The displayed four native band sizes are [3,4,4,3], totaling 14 T
gates. The already proved exact whole-word resynthesis still has two
native bands and 8 T gates. There are no inverse-padded middle axes:
the two lists are disjoint. This witnesses the failed initial matching
inference without proving an increased optimal depth or improving the
approximate-source frontier.

## C. Any sole anticommuting unmatched orbit cancels in this physical class

Retain the fixed endpoints and actual preserved r. Suppose B has one
unmatched axis u with exponent epsilon, an aligned original pair and
positive J singleton factors. More generally, require that the actual
remaining native product V_B commutes J literally. Native commutation
also makes it commute u. Thus

```math
\begin{aligned}
E_B&=T_u^\epsilon V_B,\qquad [V_B,J]=[V_B,u]=0,\\
E_BJE_B^\dagger&=T_u^\epsilon J T_u^{-\epsilon}.
\end{aligned}
```

Assume [u,r]=0 and put v=JuJ. The new case is {u,v}=0. The actual D
band has independent commuting axes, each preserving r, and exact
return gives

```math
E_D(T_u^\epsilon JT_u^{-\epsilon})E_D^\dagger=J.
```

No assumption that the intermediate image is Clifford is made.

Define the encoded Pauli triple and the commuting factor

```math
w'=-iuv,\qquad H'=(u+v)/\sqrt2,\qquad F=JH'.
```

H' and F are Hermitian Cliffords. F commutes u,v,w' and r. The source
and target reflections factor as

```math
\begin{aligned}
T_u^\epsilon JT_u^{-\epsilon}
 &=F\left(\frac{u}{\sqrt2}+\frac{v+\epsilon w'}2\right),\\
J&=F\frac{u+v}{\sqrt2}.
\end{aligned}
```

The physical role of the preserved center is essential. On its actual
r=+1 and r=-1 codes, J_core becomes the Pauli P_s=Z_q or X_q in the
fixed CNOT pair coordinates. If $`t_i=\langle u,a_i\rangle`$ and
$`t_3=\langle u,Z_w\rangle`$, then $`[u,r]=0`$ gives
$`t_1=t_2`$, while $`\{u,JuJ\}=0`$ gives
$`t_1+t_2+t_3=1`$. Thus u has w type X or Y.
Write $`u_s=A_sR_w`$, where R is X or Y, and let $`\kappa_s`$
be the commutation sign of A_s with P_s. Then the restricted F_s is
a genuine Hermitian Pauli:

```math
F_s=\begin{cases}
P_sA_s,&\kappa_s=+1,\\
-iP_sA_sZ_w,&\kappa_s=-1,\ R=X,\\
+iP_sA_sZ_w,&\kappa_s=-1,\ R=Y.
\end{cases}
```

It commutes the restricted encoded Pauli plane generated by u_s,v_s.

Restrict D to either actual r code. The projected list can lose
independence, and no independence claim is needed here. Its binary
span L_s remains isotropic. Conjugation of any Pauli by these
commuting rotations has support in its input label plus L_s.

If L_s had zero intersection with the encoded plane, the three source
cosets F_s+u_s+L_s, F_s+v_s+L_s and F_s+w'_s+L_s would be disjoint.
The nonzero third source component could not disappear into the two
target cosets. Thus the intersection contains a nonzero encoded Pauli
R_s. Isotropy makes this intersection one-dimensional, and D_s fixes
R_s because R_s belongs to its commuting axis span.

Apply the anticommutator with R_s, divided by two, to the return
identity. Since F_s commutes the encoded plane, it gives

```math
a_R D_sF_sD_s^\dagger=b_R F_s,
```

where source coefficients are (1/sqrt(2),1/2,epsilon/2) and target
coefficients are (1/sqrt(2),1/sqrt(2),0) in the u,v,w' coordinates.
Unitary norm excludes R=v (ratio sqrt(2)) and R=w' (ratio zero).
Therefore R=u and D_s fixes F_s literally. Both physical r sectors
reassemble to global fixation of u and F by E_D.

Let $`Q=E_DvE_D^\dagger`$. Global fixation of u gives
$`E_Dw'E_D^\dagger=-iuQ`$. Canceling F from the exact return identity
therefore yields

```math
\begin{aligned}
Q+\epsilon(-iuQ)&=\sqrt2 v,\\
E_DvE_D^\dagger&=(v-\epsilon w')/\sqrt2.
\end{aligned}
```

Now use independence of the original full physical D list. Its
conjugation of v has exactly 2^k nonzero equal-magnitude Pauli terms,
where k is the number of actual axes anticommuting with v. The last
display has precisely two terms, so k=1. Its second label differs
from v by u. Hence the unique actual active D axis is a signed
representative sigma*u, and its effective T orientation is -epsilon.

This axis commutes every other actual D factor, and u commutes every
other actual B factor. The two native gates therefore delete literally:

```math
T_{\sigma u}^{\delta}T_u^\epsilon
 =\omega^{\delta(1-\sigma)/2}I,
\qquad \sigma\delta=-\epsilon.
```

The reduction removes two actual T gates, introduces no helper or
error, and retains all physical columns, phases, references and the
actual inverse. It need not remove a whole band. A scalar may be
carried literally or emitted by a constant-size Clifford scalar word;
for example, $`(HS)^3=\omega I`$ on an existing wire.

The proof covers core-coupled axes with N_Ju=r+Z_w as well as local
X_w. It does not assert this cancellation for an arbitrary Clifford
involution or without preserved r: the decisive Pauli property of
F_s was derived from these actual endpoints and actual sector codes.
It also does not cover additional B factors that fail to commute J.

## D. Complete Boolean return forces the actual D band to preserve r

The center premise itself can be removed in the same fixed-endpoint
class. Assume only that B contains the actual original pair with center
r, so every actual B axis commutes r. Do not assume complete B matching
or that any actual D axis preserves r. Retain the full occupied Boolean
action and exact clean return.

Put $`L=E_AK`$, $`U=E_DE_B`$ and $`\Pi=(I+J)/2`$.
The actual center has computational
logical pullback

```math
L^\dagger rL=K^\dagger rK=Z_1Z_2Z_w.
```

It therefore commutes the complete Boolean target on the initialized
code. Together with exact clean return, this gives

```math
Q\Pi=r\Pi,\qquad
Q=UrU^\dagger=E_DrE_D^\dagger.
```

The last equality uses actual B commutation with r. The operator Q is
an actual full-space Hermitian Clifford because it is the image of a
Pauli under one independent native band. Consequently the Clifford
C=rQ fixes Pi pointwise.

Apply the retained
[pointwise stabilizer theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#1-the-exact-pointwise-clifford-stabilizer). The clean
image J has intrinsic active rank three, so C is literally I or J:

```math
Q=r\quad\hbox{or}\quad Q=rJ.
```

The second possibility is incompatible with actual D provenance.
Multiplication by Pauli r does not change binary symplectic action,
so Q=rJ would have displacement image

```math
L_A=\mathrm{span}\{a_1,a_2,Z_w\}.
```

The intrinsic-axis theorem for the actual native image $`E_DrE_D^\dagger`$
says its active D axes span exactly L_A. But r=a1+a2 lies in the
isotropic space L_A. Every axis in that space commutes r, contradicting
the claimed activity. Hence Q=r. Independent-axis flat expansion then
forces every actual D axis to commute r.

This proof does not promote an arbitrary occupied Clifford restriction
to a full Clifford: Q is independently known to be an actual one-band
Pauli image before clean-code rigidity is used. It does not impose a
Boolean unused block. It also does not use B's full signed-orbit
matching. The center is forced by the same fixed endpoints, the actual
original B center and its computationally diagonal logical pullback.

Combining this lemma with Section C removes the separate D-center
premise from that sole-anticommuting-unmatched-axis cancellation theorem
whenever the complete occupied Boolean contract is available.

## E. A paid derivative for a general original physical pair

This section allows arbitrary original endpoint bands and Cliffords.
Use the complete actual four-band source

```math
V=U_C^\dagger E_DE_BU_A,\qquad
U_A=E_AK_A,\qquad U_C=E_CK_C,
```

with exact occupied Boolean action
$`G=\alpha\,\mathrm{diag}((-1)^{g(x)})`$ and the same full clean,
dirty and reference return contract. Suppose two actual commuting
B axes form a genuine pair under $`J_A=E_AP_AE_A^\dagger`$.
Let r be their signed Hermitian product, including the pair sign.
No matching premise is imposed on the other B axes, and no center
preservation premise is imposed on D or the output endpoint.

Actual pair provenance makes $`J_ArJ_A=r`$. Its binary difference
lies in the actual A-clean-active span, so $`[E_A,r]=0`$. Therefore
$`[P_A,r]=0`$, and

```math
R_A=U_A^\dagger rU_A=K_A^\dagger rK_A,
\qquad [R_A,Z_c]=0.
```

Thus R_A is a legitimate physical Pauli preserving the initialized
clean code. This clean-evenness uses the actual pair relation; mere
membership in the active span would not suffice. The pair can have
any even support in the original A-active coordinates.

Every B axis commutes r. Put $`Q=E_DrE_D^\dagger`$, an actual
full-space Clifford obtained by conjugating a Pauli through one
independent native band. The literal full-space identity is

```math
U_C^\dagger Q U_C=VR_AV^\dagger.
```

Consequently the physically implemented circuit

```math
\mathcal D_r
 =R_AK_C^\dagger E_C^\dagger Q E_CK_C
 =R_AVR_AV^\dagger
```

has occupied action
$`(-1)^{g(x)+g(x+v)}`$, where v is the logical X mask of R_A.
All Pauli Y/Z phases and the common target scalar cancel literally.
Dirty directions in v contribute nothing when g is dirty-independent.
The complete derivative returns the clean and arbitrary dirty inputs
exactly, including external references; its unused clean block remains
unrestricted.

This is a paid construction using the two actual C bands. Synthesize
Q as its actual Clifford, with its literal scalar, rather than crediting
it as a free native emitter. On q physical wires, let T_C count the
C band T gates and G_C its existing Clifford exposure words. Then

```math
\begin{aligned}
T(\mathcal D_r)&\le2T_C,\qquad d_T(\mathcal D_r)\le2,\\
G(\mathcal D_r)&\le2G_C+2G(K_C)+O(q^2).
\end{aligned}
```

The last term pays for Q and R_A. No additional clean or dirty helper
is used, the error is zero, and the actual reversed word supplies
the inverse. The retained
[two-layer real-diagonal theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#3-the-pauli-midpoint-real-diagonal-rigidity-with-one-clean-qubit)
therefore proves

```math
\deg\bigl(g(x)+g(x+v)\bigr)\le3.
```

This constrains only the derivative in the center's logical X
direction. It does not imply a total-degree bound of four for g,
particularly when v is zero, and does not prove that D preserves r.
The endpoint differences and all original native pair signs have been
retained throughout.

The derivative also gives a robust restriction. Let the original source
have complete initialized-isometry error epsilon against G, allowing
clean leakage and dirty-return error. Pair provenance still makes R_A
clean-even. The actual inverse has the same error epsilon on the
initialized code: choose any code-preserving unitary extension of G,
write $`\Pi_0`$ for the initialized-code projector, and use unitarity in

```math
(V^\dagger-G^\dagger)\Pi_0
 =V^\dagger(G-V)G^\dagger\Pi_0.
```

Telescoping the two source calls in the exact identity for
$`\mathcal D_r`$ therefore gives complete derivative error at most
$`2\epsilon`$. This includes leakage, dirty return and references;
it does not promote approximate return to exact return. The same
physical two-C-band circuit and resource ledger apply.

Put $`\tau=\sin(\pi/8)/(2\sqrt2)`$. The proof of the retained
[robust two-layer obstruction](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction)
applies to every Boolean diagonal phase of degree greater than three,
not only to its displayed XOR specialization. Below tau, it first
rounds the clean return and all preserved Z operators, then uses the
native eighth-root phases to force the exact target up to one common
scalar. The real-diagonal theorem then forces degree at most three.
Consequently the directional degree restriction above remains valid
whenever $`\epsilon\lt\tau/2`$.

For the exact indicator sign $`g(x,Y)=Y_x`$ on p address bits, this
directional conclusion gives a concrete restriction on every original
B pair center. Write its logical X mask as (a,b), where a shifts the
address and b is a fixed bank mask; arbitrary dirty directions can be
omitted because the target is dirty-independent. With
$`\delta_t(x)=[x=t]`$, the derivative is

```math
\Delta_{(a,b)}g
 =\sum_tY_t\bigl(\delta_t(x)+\delta_{t+a}(x)\bigr)+b_{x+a}.
```

If a is nonzero, a linear address-coordinate change makes it the first
coordinate vector. Each coefficient of Y_t then becomes the indicator
of the other p-1 coordinates and has degree exactly p-1. The bank-free
term cannot cancel these Y-linear monomials, so the whole derivative
has degree exactly p. If a is zero, the derivative is simply b_x.
Consequently, for p greater than three, every such center has zero
address X mask, and its bank mask is the truth table of a Boolean
polynomial of degree at most three, that is, a word of RM(3,p).
The same restriction holds for the span of these projected center
masks because RM(3,p) is linear. This constrains only their X projection:
native independence does not bound the number of original pairs from
it, since independent physical centers can all have zero X mask.
The same mask restriction holds at complete error
$`\epsilon\lt\sin(\pi/8)/(4\sqrt2)`$ by the preceding paid derivative
and rounding argument.

## F. A paid two-band resynthesis of the fully matched mechanism

The fixed endpoint above has only clean-active A axes. More generally,
suppose every axis in the actual independent list

```math
E_A=\prod_jT_{a_j}^{\epsilon_j},\qquad
J=E_APE_A^\dagger
```

anticommutes the Pauli P. Let an actual middle band E_B be fully
matched under J, with exactly one genuine pair and any positive or
negative singletons. Suppose the pair difference uses exactly two
original A axes, a_i and a_j. These are physical native hypotheses,
not hypotheses about a freely chosen logical pair.

Write the pair as b and $`b'=JbJ`$. If its displayed second axis is c
with $`JbJ=\sigma c`$, use the literal representative conversion

```math
T_c^{\theta_0}
 =\omega^{\theta_0(1-\sigma)/2}T_{b'}^{\sigma\theta_0},
\qquad \omega=e^{i\pi/4},
```

and retain its scalar. Let eta and theta be the resulting effective
orientations. Put $`PbP=\kappa b`$. The signed pair formula is

```math
b'=-\kappa\epsilon_i\epsilon_j a_i a_j b.
```

Only a_i and a_j are active on b, so

```math
E_A^\dagger bE_A
 =\frac{b+i\epsilon_i a_i b+i\epsilon_j a_j b
              -\epsilon_i\epsilon_j a_i a_j b}{2}.
```

Multiply this expression by
$`\eta-\theta\kappa\epsilon_i\epsilon_j a_i a_j`$.
If $`\kappa\theta/\eta=1`$, the weighted sum
$`\eta b+\theta b'`$ is unchanged. If
$`\kappa\theta/\eta=-1`$, its image is

```math
i\eta(\epsilon_i a_i b+\epsilon_j a_j b)
 =F(\eta b+\theta b')F^\dagger,
\qquad F=S_{a_i}^{-\epsilon_i}.
```

Every singleton commutes every a_j: its J-fixed binary label has zero
active displacement, and the original A axes are independent.
It therefore also commutes F. Applying the preceding equality to the
commuting pair exponential, while retaining its scalar
$`e^{i\pi(\eta+\theta)/8}`$, proves the literal whole-band identity

```math
E_A^\dagger E_BE_A=F_BE_BF_B^\dagger,
\qquad F_B\in\{I,S_{a_i}^{-\epsilon_i}\}.
```

The right side is an actual native band exposed by a Clifford. Its
axes retain the original independence and commutation. This is stronger
than identifying a C3 operator or a convenient logical restriction.

If E_D meets the same hypotheses and the original endpoint bands obey
$`E_C=E_A`$, apply the identity separately to both middle bands:

```math
V=K_C^\dagger F_DE_DF_D^\dagger F_BE_BF_B^\dagger K_A.
```

This is a full-space two-band realization, including the unused clean
block. The input and output Cliffords need not be equal. If T_B,T_D
count the original middle T gates and G_B,G_D count their Clifford
exposures, its paid ledger is

```math
\begin{aligned}
T(V)&=T_B+T_D,\\
d_T(V)&\le2,\\
G(V)&\le G_B+G_D+G(K_A)+G(K_C^\dagger)\\
    &\quad+\sum_{j=B,D}\bigl(G(F_j)+G(F_j^\dagger)\bigr).
\end{aligned}
```

Each nontrivial F_j is an S rotation about one actual Pauli. Local
H/S basis changes, a CNOT parity computation onto an existing wire,
one S or S-dagger, and the actual inverse exposure implement it in
$`O(q)`$ elementary Cliffords with no helper. A signed axis reverses
the effective S orientation; a common scalar of the chosen F_j
representative cancels with its actual inverse. In the positive native
alphabet S-dagger costs three S gates.

The inverse uses the reversed literal word, with the inverse B/D
exposures and $`G(K_C)+G(K_A^\dagger)`$ in its Clifford ledger.
Its T count and depth are the same. Added error is zero. Every exact
or approximate clean, dirty and reference return contract of the
original source is preserved. If its occupied action is Boolean, the
retained two-layer real-diagonal theorem gives degree at most three.

The all-active condition is essential to this whole-word reduction.
For a general endpoint, J-singletons commute only its clean-active
axes; noncommuting passive A factors may retain the old cost. A pair
with a larger active difference also lies outside this proof. The
[Pauli-domain scout](HIGHER_SUPPORT_PAIR_BOOLEANITY.md) gives an exact
local obstruction in those cases and tests the complete Boolean
condition separately. No unrestricted four-layer collapse follows.

## G. A sole unmatched anticommuting factor must merge natively

Let J be a Hermitian Clifford involution and
$`E=B_0T_u^\epsilon`$ an actual independent commuting native band.
Put $`v=JuJ`$. Assume $`\{u,v\}=0`$ and that v commutes with
every actual $`B_0`$ axis. The whole $`B_0`$ phase need not
commute with J. Let D be another actual independent commuting
native band, with orientations $`\delta_j\in\{\pm1\}`$, and require
the actual transported probe $`J'=DEJ E^\dagger D^\dagger`$ to
be Clifford.

Then an actual D axis has label u. Writing $`d_j=\sigma u`$ gives
the full physical merger

```math
\begin{aligned}
DE&=\alpha D_{\rm rest}B_0S_u^k,\\
\alpha&=\omega^{\delta_j(1-\sigma)/2},\qquad
k=\frac{\epsilon+\sigma\delta_j}{2}\in\{-1,0,1\}.
\end{aligned}
```

If $`J'=J`$, necessarily $`k=0`$: the factors cancel. This includes
a sole factor $`E=T_u^\epsilon`$ and the earlier hypothesis
$`[B_0,J]=0`$. In the latter case $`B_0`$ fixes u and therefore v;
its independent native expansion then makes each actual axis
commute v.

For the proof, set $`U=DE`$. E fixes u, and only its u factor
acts on v. Thus

```math
\begin{aligned}
UuU^\dagger&=DuD^\dagger=C_u,\\
UvU^\dagger&=DHD^\dagger=J'C_uJ'^\dagger,\\
H&=T_u^\epsilon vT_u^{-\epsilon}.
\end{aligned}
```

Both $`C_u`$ and H are actual one-band Pauli images, hence
Clifford. The second equality independently establishes that
$`DHD^\dagger`$ is Clifford before invoking the
[one-band matching theorem](ONE_BAND_AXIS_MATCHING.md).
H has the binary transvection

```math
p\longmapsto p+\langle u,p\rangle u.
```

Any actual D axis d anticommuting with u would therefore have its
mate $`d+u`$ in the same D list. Those two labels anticommute,
which is impossible for an actual native band. Every D axis
commutes u, so now $`C_u=u`$ is Pauli. The actual Clifford J'
then makes $`DHD^\dagger=J'uJ'^\dagger`$ Pauli. This stronger
conclusion does not follow merely from its Clifford membership.

Compare H's two terms with the reverse D-band expansion. Native
independence forces exactly one active D axis, whose label is
their difference u. That common axis commutes both complete
lists, proving the displayed signed merger.

For the same-probe strengthening, U commutes J and fixes u, so it
also fixes v. After the merger,

```math
D_{\rm rest}S_u^k vS_u^{-k}D_{\rm rest}^\dagger=v.
```

The argument of $`D_{\rm rest}`$ is Pauli. An independent native
band can send a Pauli to a Pauli only by fixing it literally.
For $`k=\pm1`$, the S rotation instead gives label $`u+v`$.
Hence $`k=0`$, without assuming $`[B_0,J]=0`$.

This deletes two actual T rotations and adds at most one Pauli S
rotation, paid with $`O(q)`$ elementary basis/parity Cliffords.
The scalar is retained or emitted by a fixed Clifford scalar word.
There is no helper or added error; all physical columns, actual
inverse, dirty/reference return and remaining native independence
are preserved. No merging of the two remaining bands is inferred.

With one clean qubit, complete exact return between endpoint codes
implies full probe transport: their positive eigenspaces are half
dimensional and the orthogonal complements are transported too.
Return on a smaller code, or approximate return, does not supply
this exact hypothesis.

For involutive J, the binary matrix
$`C_{ik}=\langle a_i,Ja_kJ\rangle`$ on actual E axes is symmetric.
This result removes a row whose only nonzero entry is its diagonal.
A surviving anticommuting orbit must therefore couple to another
actual E axis through a nonzero off-diagonal entry. Commuting
unmatched orbits remain a separate possibility.

## H. Active support matching and its paid rank-four exception

For arbitrary actual independent bands $`U=DE`$, a Pauli image
$`URU^\dagger=S`$ requires $`ERE^\dagger=D^\dagger SD`$.
The flat supports have equal dimension and active span. If A and B
are their active axis matrices, equality of their Clifford
symplectic actions gives

```math
AA^{\mathsf T}=BB^{\mathsf T},\qquad
B=AO,\qquad OO^{\mathsf T}=I.
```

Signs and orientations do not change this binary displacement.
The necessary basis condition does not replace coefficient matching.
In dimensions at most three, an orthonormal binary basis is a
permutation. Thus every active E axis occurs as an actual D axis.
Remove common axes by literal Clifford mergers and write $`U=CU_0`$.
If either residual list has at most three axes, then

```math
M_{U_0}=\langle E_{\rm rest},D_{\rm rest}\rangle^\perp,
```

fixed pointwise by $`U_0`$. The original U has the same domain but
acts there by C, which cannot be dropped.

At rank four the only additional basis type, up to permutation, is
the four complementary triples $`b_i=\sum_{j\ne i}a_j`$.
A weight-one and a weight-three vector cannot coexist in a full
orthonormal basis. If this match occurs inside larger lists, its
equal active spans make all eight block rotations commute with both
complete lists. They factor literally from the whole word.

The positive four-wire block is

```math
E=\prod_{i=1}^4T_{Z_i},\qquad
D=\prod_{i=1}^4T_{\prod_{j\ne i}Z_j},\qquad
DE=(-1)^{e_1+e_2+e_3}.
```

Here $`e_j`$ are elementary symmetric Boolean polynomials. The sign
is positive at weights zero and four and negative at weights one,
two and three. The disjoint actual lists are independent. This
complete component costs eight T gates in two bands and 18 CNOTs
for the complementary-parity encoder and its actual inverse.
Arbitrary orientations change it only by paid Clifford rotations
and scalars; its eight-T bill is not erased.

For actual endpoints put $`R=X_1X_2X_3X_4`$,
$`a=Y_1X_2X_3X_4`$, $`A=C=T_a`$ and $`P=Z_1`$. Then

```math
J=APA^\dagger=\frac{Z_1+R}{\sqrt2},
\qquad JZ_iJ=Z_1Z_iR.
```

Every E axis anticommutes with its J image, which is outside E.
Nevertheless DE commutes a, and $`A^\dagger DEA=DE`$ on the full
space. This is a complete Boolean representation with native sizes
$`[1,4,4,1]`$ and ten T gates; the explicit two-band realization is
cheaper. Endpoint exposures and their actual inverses remain paid.

This rank-one endpoint supplies no original physical pair. Tensoring
it with the earlier scout's Section 11 full-space Boolean source,
and encoding the sole clean probe as the product of the old clean
probes by a CNOT and its inverse, gives a concise provenance guard.
On nine wires it has one clean, the original returned dirty wire,
sizes $`[4,8,8,4]`$ and 24 T gates. The genuine pair and four
unmatched anticommuting axes survive. The two known full-space
reductions give a paid 16-T, two-band baseline. No frontier gain
or free parallel repacking is claimed.

## I. An independent native union excludes anticommuting orbits

Section I.3 gives the strongest result: the cross-commutation rank
is bounded by native overlap on the affected E coordinates. The support
arguments below retain useful independent-union and dependent-code
components without supplying a free native replacement.

Suppose the combined actual E/D label list is independent, a
stronger condition than having no common individual axis. Write
$`W=\langle a_i,b_j\rangle`$ and
$`M_{ji}=\langle b_j,a_i\rangle`$. For a physical Pauli R use
$`e_i=\langle a_i,R\rangle`$ and
$`d_j=\langle b_j,R\rangle`$. No full expansion branches collide:

```math
N(e,d)=\sum_{s\subseteq\mathrm{supp}(e)}
                       2^{\mathrm{wt}(d+Ms)}
```

is the exact Pauli support size, for all signs and orientations.
Size one occurs exactly at $`e=d=0`$, so $`M_U=W^\perp`$ is
fixed pointwise. Size two occurs only at
$`e=0,\mathrm{wt}(d)=1`$, or when e selects one zero M column
and $`d=0`$.

Let $`A_c`$ span the actual E axes with nonzero M columns.
The span of all inputs with support size at most two is
$`A_c^\perp`$: the physical commutation map onto $`(e,d)`$ is
surjective, and these inputs span precisely the vectors whose
coupled e coordinates vanish, including its kernel $`W^\perp`$.
This also handles empty lists.

For actual Clifford transport $`J'=UJU^\dagger`$, intertwining
preserves every physical support-size class. Hence J preserves W,
its radical, and $`A_c`$. A coupled E axis maps into isotropic
$`A_c`$; a silent E axis belongs to the radical and maps there.
Therefore

```math
[a_i,Ja_kJ^\dagger]=0\qquad\text{for all }i,k.
```

The actual inverse gives the symmetric D/J' conclusion. Arbitrary
physical spectator mixing is covered. This is full-operator
eligibility, not a native resynthesis or a consequence of one
occupied matrix element.

### I.1. One actual pair on both sides permits a justified restriction

Suppose instead that the actual spans satisfy
$`A\cap B=\langle r\rangle`$. Let
$`r=a_1a_2`$ be a genuine original E-pair center, fixed positively
by J, and suppose two actual D axes obey
$`b_1b_2=\tau r`$, $`\tau\in\{\pm1\}`$.
All middle factors commute r, so the actual output probe fixes it
too. On its character $`r=s`$, the literal pair factors are

```math
\begin{aligned}
F_{E,s}
 &=\omega^{\epsilon_2(1-s)/2}
   S_{a_1}^{(\epsilon_1+s\epsilon_2)/2},\\
F_{D,s}
 &=\omega^{\delta_2(1-\tau s)/2}
   S_{b_1}^{(\delta_1+\tau s\delta_2)/2}.
\end{aligned}
```

They are Clifford factors, with their scalar phases retained.
After omitting the two actual pairs, write

```math
U_s=F_{D,s}D_{{\rm rest},s}E_{{\rm rest},s}F_{E,s}.
```

Absorb the outer factors into the restricted Clifford probes.
The combined retained axis list is independent on this code:
a relation modulo r would yield a vector in $`A\cap B`$, whereas
native independence prevents either retained list from summing to r.
Thus the preceding eligibility theorem applies legitimately.
Since $`F_{E,s}`$ fixes every retained E axis, its cross-commutation
conclusion transfers unchanged to J. Genuine pair matching supplies
the omitted cross terms. Every E axis again commutes with every
J-image of an E axis.

No D-pair matching premise is required. Positive center fixation,
the stated one-dimensional intersection and full probe transport
are required. The sector Clifford factors and analytical code
coordinates are not a free native replacement of the physical word.

### I.2. One shared center with higher native support

Retain $`A\cap B=\langle r\rangle`$, the genuine original E pair
with center r, positive fixation by J, and actual full Clifford
probe transport. Let k be the number of actual D axes in the
unique expression for r. Their signed product may be $`\tau r`$,
$`\tau\in\{\pm1\}`$. This sign is retained on both center characters.

For $`k=1`$, the D factor restricts to its literal scalar
$`\omega^{\delta(1-\tau s)/2}`$ on $`r=s`$. Remove that scalar and
the E-pair Clifford exactly as in I.1. The remaining combined list
is independent on the code, so Section I excludes anticommuting
E orbits. The $`k=2`$ case is I.1. Neither restriction is a free
physical resynthesis.

For $`k\ge3`$, remove only the E pair as its scalar-times-Clifford
factor $`F_E`$. Write the restricted D list as $`D_0D_{\rm rel}`$,
where $`D_{\rm rel}`$ contains the k axes expressing r. The
remaining E span and the restricted D span intersect trivially.
The combined lists have exactly one relation, internal to
$`D_{\rm rel}`$. These statements follow from the original
intersection and actual-list independence; they are not assumed
after the quotient.

For a code Pauli, let c be its commutation vector against
$`D_{\rm rel}`$. It has even parity. Every proper active subset is
independent, so its Pauli image has $`2^{\mathrm{wt}(c)}`$ terms.
If all k axes are active, k is even. Absorbing the actual signs and
orientations into $`C_j`$, the Clifford product in that image is,
up to its retained common phase,

```math
2^{-k/2}\prod_{j=1}^k(I-iC_j),
\qquad \prod_{j=1}^kC_j=\nu I,\quad \nu\in\{\pm1\}.
```

A subset S and its complement share a code label. Their combined
coefficient contains
$`1+\nu(-1)^{k/2+|S|}`$. Exactly one parity class survives.
There are $`2^{k-2}`$ surviving labels, all with coefficient
magnitude $`2^{1-k/2}`$. Thus the exact support function is

```math
F(c)=
\begin{cases}
2^{\mathrm{wt}(c)},&c\ne\mathbf1,\\
2^{k-2},&c=\mathbf1\quad(k\text{ even}).
\end{cases}
```

In particular $`F(0)=1`$, every nonzero allowed c has
$`F(c)\ge4`$, and weight-two c has $`F(c)=4`$. For $`k=3`$ only
weights zero and two occur. For $`k=4`$, full activity also leaves
four terms, each of magnitude $`1/2`$. The signed parity coset can
change, but these counts cannot.

Let e and d be commutation vectors against the retained E list and
$`D_0`$, and let $`M_0,M_c`$ be their cross matrices with E.
Distinct E branches lie in distinct cosets of the D span. Hence
the complete restricted support size is exactly

```math
N(e,d,c)=
\sum_{s'\subseteq\mathrm{supp}(e)}
  2^{\mathrm{wt}(d+M_0s')}F(c+M_cs').
```

Size one occurs precisely at $`e=d=c=0`$, fixed pointwise by the
residual middle. Define $`A_{\rm heavy}`$ as the span of actual
retained E axes i for which

```math
(M_c)_{:,i}\ne0
\quad\text{or}\quad
\mathrm{wt}((M_0)_{:,i})\ge3.
```

The span of all code inputs with $`N\le4`$ is exactly
$`A_{\rm heavy}^{\perp}`$. To see this, at $`e=0`$ the d units and
weight-two c vectors span every allowed D commutation coordinate.
For one active e coordinate, a nonzero $`M_c`$ column makes
$`N\ge5`$; if that column is zero but $`c\ne0`$, then $`N\ge8`$.
With $`c=0`$, two powers of two sum to at most four only for an
$`M_0`$ column of weight zero, one or two. Each such e coordinate
does occur: choose $`d=0`$ for weights zero and one, and one
supported d unit for weight two. Two active e coordinates require
both columns to be wholly silent and $`d=c=0`$; more give at least
eight terms. The commutation map onto $`(e,d,c)`$ is surjective
subject only to c's even parity, and its kernel already has
support one. This proves equality of spans, including empty lists.

Clifford intertwining preserves support-size classes, so the
adjusted input probe preserves $`A_{\rm heavy}`$. This space is
isotropic. Every heavy E axis therefore commutes with every probe
image of the retained E list. The conclusion transfers back through
$`F_E`$ because it fixes that entire list, and then to physical
commutation because adding r does not change it. Actual pair
matching supplies the omitted cross terms. Fully silent retained
axes also cannot have anticommuting orbits: they lie in the radical
of the residual span, preserved by the adjusted probe.

At the level of this filtration, a remaining anticommuting E orbit
would have to commute with every actual $`D_{\rm rel}`$ axis and
couple to exactly one or two
actual $`D_0`$ axes. If $`D_0`$ is empty, none remains, for any k.
The next native-defect argument also closes those remaining
couplings in the stated one-center class.

An invariant heavy span does not make its actual native phase
commute with J, nor permit that phase to pass through $`D_0`$.
There is no free heavy-block deletion. The actual one-native output
endpoint, complete Boolean action on all dirty/reference columns,
and a better whole-word schedule remain unpaid. No count, depth,
error, inverse or return recurrence has changed.


### I.3. Native defects close the one-center cross-commutation branch

The same defect used in G strengthens the preceding filtration.
Let $`U=DE`$ be the actual two-band middle, J a Hermitian Clifford
involution, and $`J'=UJU^\dagger`$ the actual Clifford output
probe. For the actual E axes define

```math
\begin{aligned}
C_{ji}&=\langle a_j,Ja_iJ\rangle,\\
S&=\{j:C_{j,:}\ne0\},\qquad
A_S=\langle a_j:j\in S\rangle,\\
H&=\langle E\rangle\cap\langle D\rangle.
\end{aligned}
```

C is symmetric because J has involutive symplectic action. The
quantitative necessary condition is

```math
\mathrm{rank}\,C\le\dim(H\cap A_S)\le\dim H.
```

For each i, both $`H_i=E(Ja_iJ)E^\dagger`$ and
$`D a_iD^\dagger`$ are actual one-band Pauli images, hence
Clifford. The identity

```math
DH_iD^\dagger=J'(Da_iD^\dagger)J'^\dagger
```

independently proves Clifford output before applying one-band
matching. Thus $`H_i`$ permutes individual actual D-axis labels.
Its binary symplectic action is

```math
\begin{aligned}
F_i&=I+N_i,\\
N_i&=\sum_j C_{ji}a_ja_j^{\mathsf T}\Omega,
\qquad N_iN_\ell=0.
\end{aligned}
```

The last identity follows from E isotropy. These binary actions
commute; no literal commutation of the Clifford defects is asserted.
For an actual D axis d, write $`m_j=\langle a_j,d\rangle`$.
Adding an E vector does not change m, so its full binary orbit is

```math
\begin{aligned}
d+L_d,\qquad
L_d&=\left\langle
       \sum_jm_jC_{ji}a_j:\text{all }i
     \right\rangle .
\end{aligned}
```

Every orbit point is an actual D-axis label. If $`\dim L_d\ge2`$,
four distinct such labels
$`d,d+h_1,d+h_2,d+h_1+h_2`$ sum to zero, contradicting native
D independence. Every defect-group orbit is therefore a singleton
or a pair, and

```math
\mathrm{rank}\bigl(\mathrm{diag}(m)C\bigr)\le1.
```

A pair's nonzero displacement belongs to H. Displacements from
different orbits are independent, since their relation would be a
relation among disjoint pairs of actual D axes. These are orbits
of the defect group, not automatically original J-matched pairs.

Over the binary field, all nonzero C rows selected by m must be
identical. If R is that row class, the pair displacement is

```math
h=\sum_{\substack{j\in R\\m_j=1}}a_j.
```

The two actual D axes commute, so
$`0=\langle d,h\rangle=\sum_{j\in R}m_j`$. Each displacement
has even support, at least two actual E axes, within one row
class. Symmetry makes identical rows identical columns, so its
E-coordinate vector also lies in the kernel of C. All statements
retain the actual signed pair center and gate orientations; no
native pair phase has disappeared.

Call an E axis silent if it commutes every actual D axis and
coupled otherwise. For a silent $`a_i`$, U fixes $`a_i`$
literally. Clifford transport makes $`U(Ja_iJ)U^\dagger`$
Pauli. Comparing the flat active supports in

```math
E(Ja_iJ)E^\dagger=D^\dagger(\text{Pauli})D
```

puts every C-neighbor $`a_j`$ of i in H: the affine support's
difference span on either side is its actual active-axis span.
Every such neighbor is silent, since H lies in the isotropic D
span. Thus C has no silent-to-coupled edge. If i is also active,
choose one neighbor j. Symmetry and the same argument at the
silent j put $`a_i`$ itself in H.

Let s count active silent E axes and t count distinct nonzero
C row classes consisting of coupled axes. A nonzero row class
cannot mix silent and coupled vertices: a silent row has a
nonzero silent neighbor, and an equal coupled row would violate
the absence of mixed edges. Every coupled row class is touched
by an actual D axis, so it supplies a nonzero h in H supported
entirely on that class. Choose one such h per class.

These t displacements and the s active silent axes have disjoint
actual E-coordinate supports, hence are independent in
$`H\cap A_S`$. Their s+t row types span C. Consequently

```math
\mathrm{rank}\,C\le s+t\le\dim(H\cap A_S)\le\dim H,
```

as claimed. In particular, $`H\cap A_S=\{0\}`$ forces C=0.

Now suppose H is exactly the span of the genuine original pair
center $`r=a_1+a_2`$. Actual J matching of those two axes makes
both corresponding C rows zero. Native E independence then gives
$`\langle r\rangle\cap A_S=\{0\}`$. Consequently

```math
\langle a_i,Ja_jJ\rangle=0\qquad\text{for every }i,j.
```

This closes the light-coupling escape left by I.2 for the complete
one-center class, regardless of the length of r's actual D
expression or the number of outside D axes.

Cross-commutation still does not make J permute actual E axes.
The defect pairs become scalar-times-Clifford only on their center
characters; simultaneous restriction may introduce relations in
other actual lists, and J need not preserve those characters.
There is no free global Clifford or independent native quotient.
The rank bound allows arbitrary width and spectators, and its
overlap dimension can grow with width. It is not an unrestricted
depth bound.

Commuting unmatched images can remain, and no heavy native phase
has been deleted or moved through another band. The actual output
endpoint, complete Boolean action, and better whole-word bill
remain separate obligations. No resource, error, actual-inverse
or return recurrence changes.

### I.4. A nonsingular silent block has a paid Clifford removal

Retain I.3's actual-band and Clifford-probe hypotheses. Let
$`S_0`$ index the C-active E axes silent on D, put
$`s=|S_0|`$ and $`A_0=\langle a_i:i\in S_0\rangle`$, and
let $`C_0`$ be C restricted to these indices. If $`C_0`$ is
nonsingular, a literal full-space Clifford factor removes
$`2s`$ T gates from the actual two-band middle.

For each $`i\in S_0`$, U fixes $`a_i`$, and

```math
S_i=U(Ja_iJ)U^\dagger=J'a_iJ'^\dagger
```

is Pauli. Flat active-support matching in

```math
E(Ja_iJ)E^\dagger=D^\dagger S_iD
```

puts every active D axis individually in the span of i's actual
E neighbors, hence in $`A_0`$. Those neighbor spans together
span $`A_0`$: I.3 excludes edges to coupled axes, and $`C_0`$
has no zero rows. The actual D axes in $`A_0`$ therefore form
exactly one independent basis of size s. Call their subband
$`D_0`$, and call the corresponding E subband $`E_0`$.
They commute every residual native factor, giving the literal
factorization

```math
U=F U_{\rm rest}=U_{\rm rest}F,\qquad F=D_0E_0.
```

The original F word has $`2s`$ T gates in two actual bands.
The same support identity gives
$`F(Ja_iJ)F^\dagger=S_i`$: the residual E axes fix $`Ja_iJ`$,
and all D axes active on $`S_i`$ are in $`D_0`$. F also
fixes $`A_0^\perp`$ pointwise. Nonsingularity of $`C_0`$
makes the labels $`Ja_iJ`$ a basis of the dual quotient by
$`A_0^\perp`$. Together with $`A_0^\perp`$, they span the full physical
Pauli space, so F is a genuine Clifford.

Synthesize that exact Clifford, including its literal scalar, in
$`O(q^2)`$ Clifford gates and no helper. To retain the scalar,
expose $`A_0`$ by a Clifford and normalize the resulting diagonal
Clifford to phase one at its all-positive character. The original
commuting T word has an $`\omega^k`$ phase there; emit this
remaining scalar with $`(HS)^3=\omega I`$ if needed. Existing
band-exposure Cliffords may all be retained after deleting the
selected rotations. Thus the complete-word ledger satisfies

```math
\begin{aligned}
T_{\rm new}&=T_{\rm old}-2s,\qquad
d_{\rm new}\le d_{\rm old},\\
G_{\rm new}&\le G_{\rm old}+O(q^2).
\end{aligned}
```

The actual inverse uses its own previous Clifford count plus
$`O(q^2)`$. This is a full-column identity, including arbitrary
dirty inputs and references; workspace, approximation and return
errors are unchanged.

If $`C_0`$ is singular, this proof does not remove F. The
rank-four commuting cubic phase in Section H is the retained guard
against calling active-support matching a free Clifford
replacement. Even the nonsingular removal is a source-specific
component, not an unrestricted depth bound or a change to the
approximate-indicator frontier.

<a id="g-next-structural-question-and-the-unchanged-frontier"></a>

## J. Next structural question and the unchanged frontier

The sole-factor merger is paid, while the first multi-factor
matching exception is a commuting cubic component retaining the
old two-band cost. Extending that component does not change the
resource frontier.

For a genuinely noncommuting unmatched branch, the next necessary
ingredient is shared native overlap on the C-active E coordinates
of dimension at least the cross-commutation rank, not only the
center of a matched physical pair. When that center is already
shared, another independent overlap is needed; changing
the original center's transport is another route. Commuting
unmatched completions remain open as well. Any proposed completion
must supply actual endpoints and full Boolean return, retain all
native phases and inverses, and compare its entire cost with paid
resynthesis. Span invariance does not pay the heavy native phase.

Retain original physical-pair provenance, use the paid directional
derivative rather than an unsupported total-degree bound, and
compare any completed word with whole-word resynthesis. The complete
approximate dirty-indicator resource, error and return recurrences
remain unchanged.

The original center-forcing statement remains scoped to its fixed
endpoint Boolean class. General endpoint motion and a nondiagonal
logical center are outside that statement. The
[Pauli-domain criterion](HIGHER_SUPPORT_PAIR_BOOLEANITY.md#1-a-necessary-pauli-domain-condition-for-diagonal-completion)
still tests full Boolean completion separately.

The [focused exact fixture](../../tests/test_four_band_original_pair_scout.py)
checks the new signed repartition, its actual native encoders, all 32
physical columns, the reversed inverse, and the shorter full-word
realization. The structural conclusions in A, C and D follow from the
proofs, not from a finite search. The main approximate-indicator
resource, error and return recurrences remain unchanged.
