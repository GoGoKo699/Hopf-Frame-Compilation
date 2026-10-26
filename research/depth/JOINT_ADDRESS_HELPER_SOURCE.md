# Joint address, carrier and helper mixing on the first quartic face

[Ordered matching source](ORDERED_MATCHING_PROBE_SOURCE.md) ·
[Two outer cubic conjugations](CLIFFORD_INVOLUTION_SOURCE.md#two-outer-conjugations) ·
[Returned cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md)

A paid twenty-band trial mixes address, carrier and arbitrary helper
coordinates through Bell changes and swaps. It creates the required quartic
matching phase, but also moves carrier bits. Two anticommuting Pauli
coefficients prove that no final Clifford wall completes this trial.
The actual missing carrier correction is paid below, giving a complete
finite-face phase circuit. Neither result supplies an all-size shallow
lookup source or changes R96.

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
