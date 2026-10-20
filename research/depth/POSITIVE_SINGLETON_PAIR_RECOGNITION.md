# Boolean recognition of an all-active pair from positive singletons

[Literal original-pair repair](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[General occupied computational center](OCCUPIED_BOOLEAN_PAIR_CENTER.md#4-computational-center-and-a-finite-feature-description) ·
[All-active diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)

**Later extension:** [complete Boolean pair collapse](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md)
removes both the positive-singleton and one-pair premises. Every paired
complete Boolean source has exactly one all-active pair and no negative
singletons after repair. The proof below remains a useful component;
its final negative-singleton reopening question is now closed.

For one genuine original middle pair, the absence of negative singletons
already forces the all-active support condition under the complete
Boolean contract. The proof uses the pair's two literal parity branches
and a global commutation identity for Boolean Pauli derivatives. It does
not use the all-active normal form to establish that support condition.

## 1. Original physical hypotheses

Use the actual one-clean three-native-band word

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,\qquad
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger.
```

Each phase list is an independent commuting list of signed Hermitian
Paulis, with its literal T/T-dagger orientations. Assume exact clean
return and complete occupied action

```math
G=e^{i\theta}\,\mathrm{diag}((-1)^{g(x)})
```

on the entire logical and dirty register, including arbitrary reference
inputs. The unused clean block is unrestricted. Suppose the original
middle list has exactly one genuine pair B,D under
$`J=E_AP_AE_A^\dagger`$, with $`JBJ=sD`$ and
$`r=b+d\ne0`$. After the literal orientation repair below, require
every remaining original singleton U to be **positive**:
$`JUJ=U`$. No bound is imposed on active ranks, singleton or silent
lists, or dirty width.

Then r is the sum of every actual clean-active axis at both endpoints.
In particular, the [all-active diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)
applies: sign/XOR degree is at most six/five uniformly, and at most
five/four when the common active rank is at least four. A proper-support
Boolean source in this one-pair setting must therefore have a negative
singleton. This conclusion concerns actual physical lists, not a choice
of basis in their spans.

## 2. Equal clean images and the missing active direction

For the original pair orientations epsilon,delta, use the exact repair

```math
F=T_B^\epsilon T_D^{s\epsilon},\qquad
C_{\rm rep}=T_D^{\delta-s\epsilon},\qquad
E_B=C_{\rm rep}FE.
```

The actual Clifford $`C_{\rm rep}`$ is absorbed into the output
encoder and phase list as in the original-pair proof. It fixes B,D and r.
Retain corrected C-side notation. The singleton band E has its original
axes and orientations; positivity means that it commutes J literally.
The repaired pair F also commutes J. Hence the actual clean comparison is

```math
E_CP_CE_C^\dagger=J=E_AP_AE_A^\dagger.
```

Let $`L_A,L_C`$ be the endpoint clean-active spans. The intrinsic
displacement images of these independent native lists give
$`L_A=L_C=L`$. Their complete nonzero Pauli supports are the affine
spaces $`p_A+L`$ and $`p_C+L`$, so $`p_A+p_C\in L`$.
These are the same equal-image facts used in the earlier all-active
proof, obtained here directly from the assumed positive singletons.

Define two linear functionals on the actual isotropic span L:

```math
\mu(l)=\langle p_A,l\rangle,\qquad
\lambda(l)=\langle b,l\rangle.
```

Both are nonzero. Every actual active A axis has mu equal to one;
the genuine pair has
$`r=(F_J-I)b=\sum_{a\in A_{\rm active}}\langle a,b\rangle a\ne0`$,
so lambda cannot vanish. If $`\lambda\ne\mu`$, the two nonzero
binary functionals are independent. Choose an actual Pauli label t with

```math
t\in L,\qquad \langle p_A,t\rangle=0,\qquad
\langle b,t\rangle=1.
```

Choose a signed Hermitian representative $`P_t`$. It commutes both
clean Paulis, since $`p_A+p_C\in L`$ and L is isotropic. It also
commutes **every actual endpoint phase axis**, active or silent: both
complete endpoint lists commute their common active span L. Every
singleton commutes L as well. Indeed, its binary fixation by J and
independence of the A-active axes force each of its pairings with those
axes to vanish. Thus the full singleton band E commutes $`P_t`$.

The contradiction below excludes this t. It will follow that
$`\lambda=\mu`$.

## 3. Restrict the clean wire only after the full branch identities

Set $`\Gamma=sBD`$, $`\Pi_\pm=(I\pm\Gamma)/2`$,
$`\eta=(1-s)/2`$ and $`\omega=e^{i\pi/4}`$. The literal native
identity, including its scalar, is

```math
F=\omega^{-\epsilon\eta}
 \left(\Pi_+S_B^\epsilon+\Pi_-\omega^\epsilon I\right).
```

The general [occupied-pair center theorem](OCCUPIED_BOOLEAN_PAIR_CENTER.md#2-a-genuine-center-at-bounded-original-pair-count)
applies before any all-active conclusion. It makes Gamma a genuine
transported signed Pauli commuting both clean Paulis and every original
phase band. Its label r lies in L. It is independent of either clean
Pauli: an actual active axis commutes r and anticommutes that endpoint's
clean Pauli. Therefore its logical restriction is nontrivial and both
logical Gamma eigenspaces are nonempty.

Use the proved shared computational-center normalization from that
general theorem. Its extracted quadratic Boolean phase and fixed
Cliffords retain a complete Boolean diagonal target. The endpoint
images of Gamma become computational-Z Paulis. Their actual transport
through a diagonal target then makes them the **same signed** Pauli,
denoted $`\Gamma_0`$. This uses no separate frame on different parity
sectors. There is no restriction by the entire center.

Retain the resulting endpoint notation and write the actual logical
Paulis as

```math
t_A=K_A^\dagger P_tK_A,\qquad
t_C=K_C^\dagger P_tK_C,
\qquad T=Gt_AG^\dagger.
```

Here and below endpoint expressions denote their legitimate positive
clean restrictions. They are actual Hermitian logical Paulis because
$`P_t`$ commutes both physical clean Paulis. Both commute
$`\Gamma_0`$, and so does T.

The needed comparison follows first on the full physical space: commute
$`P_t`$ through E_A and E, conjugate it by the actual F, and then
through the output band. Restrict this complete clean-preserving
identity, not an isolated non-Pauli code operation. Since E_C commutes
$`P_t`$, the scalar minus branch and the quarter-turn plus branch give

```math
T\Pi_{0,-}=t_C\Pi_{0,-},\qquad
(Tt_C+t_CT)\Pi_{0,+}=0,\qquad
\Pi_{0,\pm}=(I\pm\Gamma_0)/2.
```

For the second identity, $`P_t`$ anticommutes with B, so
$`S_B^\epsilon P_tS_B^{-\epsilon}`$ anticommutes with
$`P_t`$. All surrounding factors in this comparison commute
$`P_t`$. The literal scalar cancels against the actual adjoint in
each conjugation; it has not been discarded from the source word.

## 4. Boolean derivatives cannot switch their commutation sign

Conjugation by any diagonal unitary preserves the computational X mask
of a Pauli. Because $`\Gamma_0`$ is computational Z and its minus
eigenspace is nonempty, the first branch identity forces t_A and t_C
to have the same X mask u: their actions on any basis vector in that
sector have nonzero amplitudes at the same output basis vector.

For a Boolean diagonal G, define its literal diagonal derivative

```math
\Delta_u(x)=(-1)^{g(x)+g(x+u)}.
```

The common scalar of G cancels, and

```math
T=t_A\Delta_u,\qquad \Delta_u(x+u)=\Delta_u(x).
```

Thus Delta_u commutes with every Pauli of X mask u, including t_C.
If sigma is the global Pauli commutation sign of t_A and t_C, then

```math
Tt_C=\sigma t_CT
```

on the **entire** logical space. The minus-branch equality forces
$`\sigma=1`$: both products restrict there to identity. Hence T and
t_C commute globally. Their plus-branch anticommutation would then
give $`2Tt_C\Pi_{0,+}=0`$, impossible because Tt_C is unitary and
$`\Pi_{0,+}\ne0`$. This contradiction uses the complete Boolean
action; clean return or a general complex diagonal target alone would
not give the derivative invariance used here.

Consequently $`\lambda=\mu`$. Every actual active A axis pairs
with b by one, so the independent-native displacement formula gives
$`r=\sum_{a\in A_{\rm active}}a`$. The same holds for the actual
C list: its clean image is J, and its clean functional agrees with mu
on L because $`p_A+p_C\in L`$. Thus
$`r=\sum_{c\in C_{\rm active}}c`$ too. Finally undo the literal
orientation repair. Its Clifford fixes r and transports the original
C-active list bijectively, so this is also the all-active condition on
the original output list, not only the repaired one.
The common active rank is now even and at least two, since the genuine
pair has $`0=\langle b,r\rangle=k\pmod2`$. No even-rank premise
was used in the contradiction.

## 5. Consequence and remaining physical question

The all-active theorem now supplies the stated six/five uniform bounds,
with five/four when the common active rank is at least four. Its
original-circuit rounding consequence at error below $`2^{-530690}`$
also applies when the resulting actual physical comparison satisfies
the present one-pair, positive-singleton hypotheses. No arbitrary
fixed-accuracy or unrestricted three-layer conclusion follows.

The proof preserves literal phases, actual inverses, arbitrary logical,
dirty and reference inputs, exact clean return and unrestricted unused
clean action. It adds no initialized helper or approximation error.
Fixed coordinate Cliffords have the retained $`O((Q+1)^2)`$ costs on
the same register if emitted. This is an analytic necessary condition,
not a new compiler or a conclusion from finite enumeration.

The remaining one-pair problem must contain negative singletons. They
give $`J_C=E_-J_AE_-^\dagger`$ rather than equal clean images, so
the common-L argument and the Pauli commuting both endpoint bands are
not supplied. Distinguishing negative singletons that coincide with
actual endpoint axes from product axes is a possible next question;
removing the latter must not silently replace a dependent phase list
by a native band. That question, several original pairs, unrestricted
one-clean three-layer degree and dimension-growing depth remain open.
