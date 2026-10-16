# Four-layer representations: exact scope checks

[Session outcome](ONE_CLEAN_FOUR_LAYER_SESSION_RESULT.md) ·
[Outer-overlap rigidity R97](OUTER_OVERLAP_RIGIDITY.md) ·
[Retained R77 example](ONE_CLEAN_INNER_SPREAD.md#5-a-common-compatible-flip-need-not-exist)

These analytic examples delimit reductions considered in the 16 October
2026 session. They neither disprove the unrestricted occupied gap nor
give a dimension-growing depth lower bound. The
[exact fixtures](../../tests/test_outer_overlap_rigidity.py) check native
axis compatibility and the signs, occupied matrices and identities in
Sections 2–4. The other statements below are analytic arguments.

## 1. A finite inner alphabet leaves an outer coupling problem

Let C,P be the commuting encoded clean and data Paulis before an inner
physical phase layer D, whose independent isotropic axis span is L.
The symplectic form on $`L+\mathrm{span}\{C,P\}`$ has rank at most
four. One Clifford coordinate change therefore identifies its Pauli
algebra with at most two noncommuting qubits, an unbounded commuting
center and an untouched multiplicity factor.

If there are r noncommuting qubits, $`r\le2`$, the inner axes can be
written $`Z^u\otimes Z^{\lambda_i}`$. On a center character z, define

```math
k_u(z)=\sum_{i:u_i=u}\epsilon_i(-1)^{\lambda_i\cdot z},
\qquad
D_z=\prod_{u\ne0}\exp(-i\pi k_u(z)Z^u/8).
```

The omitted u=0 factor cancels in conjugation. The joint clean/data
block depends only on the two Pauli signs and the three possible
nonzero coefficients modulo eight. Thus it has at most
$`4\cdot8^{2^r-1}\le2048`$ possible values. Its selector projectors
have finite Fourier expressions: for each coefficient use

```math
\frac18\sum_{j=0}^7\zeta_8^{-j\kappa_u}
\exp\!\left(\frac{i\pi j}{4}
              \sum_{i:u_i=u}\epsilon_i Z^{\lambda_i}\right),
```

and multiply by the two sign projectors. These central exponential
factors are Cliffords, but their symplectic ranks need not be bounded.
An alphabet bound is not a bound on center multiplicity or its mixing
through the two outer layers.

The exact typed correlations also survive. Partition the inner axes
by $`(\langle a,C\rangle,\langle a,P\rangle)`$ and put
$`G_{uv}=D_{uv}^2`$. Then

```math
J_C=G_{10}G_{11}C,\qquad
J_P=G_{01}G_{11}P,\qquad
J_{CP}=G_{10}G_{01}CP.
```

In binary coordinates the two displacement maps have images in L
and their products vanish in either order. The common 11 class cancels
in the product; it is not an independent boundary. These identities
alone do not control the occupied discrepancy after the outer layers.

## 2. One initialized parity can have unbounded center mixing

Use the actual two-qubit R77 encoders $`U_A,U_B`$, with
$`U'_B=U_B(IX)`$. Copy their wires as $`a_i,b_i`$, for
$`1\le i\le r`$. Only $`a_1`$ is initialized. Let the input Clifford
$`C_r`$ be the fanins $`\mathrm{CX}(a_i,a_1)`$ and
$`\mathrm{CX}(b_i,b_1)`$, for $`i\ge2`$. They send the physical
clean and tested data Z to their respective parity products. Set

```math
U_{A,r}=(\bigotimes_i U_{A,i})C_r,\qquad
U_{B,r}=(\bigotimes_i U'_{B,i})C_r.
```

Both encoders have two physical T layers. All other wires are arbitrary
dirty inputs. Their common clean image and data images are

```math
S_r=\bigotimes_i S_i,\qquad A_r=\bigotimes_i A_i,
\qquad B_r=\bigotimes_i(-B_i).
```

The occupied condition is $`\prod_i s_i=+1`$, rather than r initialized
local sectors. Both displayed inner clean spreads, the outer crossing
rank and the outer-intersection dimension equal r. The local outer
lists are IY,ZY on A and IZ,ZZ on B; their crossing matrix is the
all-ones two-by-two matrix.

For $`d=ZI`$, exact Pauli multiplication gives

```math
SdS=-d,\qquad BdB=d,\qquad AdA=(XY-YI)/\sqrt2.
```

Thus the clean image normalizes the common radical center, but A does
not preserve it. With $`\Pi_z=\prod_i(I+z_i d_i)/2`$, one obtains

```math
\Pi_w A_r\Pi_z A_r\Pi_w=2^{-r}\Pi_w
\qquad\text{for every }w,z.
```

Indeed each nonempty product of the $`A_id_iA_i`$ is off-diagonal in
at least one center character. All $`2^r`$ characters mix, even though
the clean action has character orbits of size two. This refutes
automatic data invariance or bounded mixing of that center.

There is no finite alphabet of occupied relative eigenvalues either.
On the all-local-plus subspace, the eigenvalues of $`A_rB_r`$ include

```math
e^{i(r-2j)\alpha},\quad 0\le j\le r,\qquad
\cos\alpha=\frac{1+1/\sqrt2}{2}.
```

These are distinct: if $`\alpha/\pi`$ were rational,
$`2\cos\alpha=1+1/\sqrt2`$ would be an algebraic integer, contradicting
its norm 1/2. However the occupied defect stays bounded away from zero.
For odd r it is at least $`\sqrt{1-1/\sqrt2}`$; for even
$`r\ge2`$ use the eigenvalue $`e^{2i\alpha}`$. The maximum distance
from one tends to two along either parity. This is a stress test of a
representation claim, not a small-norm counterexample.

## 3. The same family has a simpler native representation

The displayed growing parameters in Section 2 are removable. Write
$`R_P=\exp(-i\pi P/8)`$ and use the R77 Pauli-rotation words,
before converting them to physical T gates. Define

```math
G=R_{IY}R_{ZY},\qquad
C'_B=Z_b\,\mathrm{CZ}_{ab}H_a,\qquad C'_A=C'_B S_bX_b.
```

The following identities include their literal phases:

```math
GU_A=\zeta_8(R_{YY}R_{ZZ})C'_A,\qquad
GU_B=i(R_{IZ}R_{ZZ})R_{ZX}C'_B.
```

For A, the common gauge gives $`R_{IY}^2R_{YY}R_{ZX}C_A`$;
the Clifford $`R_{IY}^2`$ fixes YY and sends ZX to -ZZ. The remaining
Clifford identity is
$`R_{ZZ}^{-2}S_a=\zeta_8Z_b\mathrm{CZ}_{ab}S_b`$.
For B, put $`M=C'_BC_B^\dagger=Z_b\mathrm{CZ}_{ab}S_a^\dagger`$.
The required identity is block diagonal in $`Z_a`$. Its -1 block is
immediate from $`M=-iI`$; its +1 block is

```math
R_Y^2R_Z^{-2}R_X=iR_Z^2R_XZ=R_Z^{-2}R_X^\dagger,
```

which follows from
$`R_Z^2R_Y^2R_Z^{-2}=R_X^{-2}`$ and $`iR_Z^2Z=R_Z^{-2}`$.
The shared midpoint after this gauge is the Hermitian Clifford
$`(YI+XZ)/\sqrt2`$.

Canceling the two ZZ rotations and conjugating by $`C'_B`$ yields

```math
U_B^\dagger U_A
=\zeta_8^{-1}R_{IX}R_{IZ}^\dagger R_{ZX}^\dagger(I\otimes SX).
```

All three rotation axes commute with the physical clean ZI. The signed
B variant adds only an IX Clifford on the left. Tensoring these words
and using the fanins gives a three-layer realization of every member
of Section 2, with a Pauli clean midpoint and minimum inner spread zero.
The earlier center-mixing identity remains true for its specified
representation. It cannot obstruct an existential change of split.

## 4. Exact same-space phase-polynomial descent is false

Use three data wires 0,1,2 and one clean wire a=3. Let $`C_0`$ apply
$`H_a`$ followed by $`\mathrm{CX}(a,j)`$ for j=0,1,2, so

```math
C_0|z,0\rangle=(|z,0\rangle+|\bar z,1\rangle)/\sqrt2.
```

Put $`D=T^{\otimes4}`$, $`E=R_{XXXX}R_{YYYY}`$ and
$`V=C_0^\dagger D^\dagger EDC_0`$. E is one physical layer after
Clifford routing: in application order use

```
CX(2,1), CX(3,1), CX(0,1), H(0), CX(0,1), CX(0,2), CX(0,3).
```

This B sends Z0,Z1 to XXXX,YYYY, so
$`B(T_0T_1)B^\dagger=\zeta_8 E`$. The full word has three native
T layers. The clean core
$`J=D X^{\otimes4}D^\dagger=((X+Y)/\sqrt2)^{\otimes4}`$
exchanges XXXX and YYYY, hence commutes with E. Both physical clean
sectors return exactly.

On the occupied cat state for a string of weight t, XXXX·YYYY has
eigenvalue $`(-1)^t`$. For even t, $`D^\dagger XXXX D`$ has
eigenvalue $`\zeta_8^{2t-4}`$; for odd t the two conjugated axes cancel.
The occupied action is therefore diagonal with phase

```math
\zeta_8^{f(z)},\qquad
f(z)=1-z_0-z_1-z_2+2z_0z_1z_2\pmod8.
```

Its determinant is $`\zeta_8^{1-3}=-i`$. Substituting the literal
physical-T routing multiplies this eight-dimensional block by
$`\zeta_8`$ and leaves that determinant unchanged.

Every three-qubit weighted-cubic phase polynomial has form

```math
p(x)=c+\sum_i l_ix_i+2\sum_{i\lt j}q_{ij}x_ix_j+4t x_0x_1x_2.
```

Its total exponent is $`8c+4\sum_i l_i+4\sum_{i\lt j}q_{ij}+4t`$,
so its determinant is ±1. Standard three-qubit Clifford representatives
also have determinant ±1. Any finite word in these operations has that
same determinant restriction. Equality up to common phase cannot evade
it: the standard word would be diagonal with entries in
$`\mathbb Z[i,1/\sqrt2]`$ of modulus one. The retained unit-circle
lemma makes them eighth roots, hence the common phase an eighth root
whose eighth power is one.

Thus exact descent to the original three logical wires is impossible,
even with unlimited phase-polynomial layers. One appended arbitrary
dirty wire gives determinant -1, while every standard four-qubit
Clifford/weighted-cubic word has determinant +1, so it is still excluded.
Two appended dirty wires give determinant +1 and remove this particular
obstruction. A width-enlarged descent theorem remains possible. The
example has finite projective order and does not refute an approximation
floor.

## 5. An intrinsic Clifford-image centralizer can be trivial

For a Hermitian involution S, consider the binary subgroup

```math
\mathcal G(S)=\{r:S R S\text{ is Clifford}\}.
```

If $`S=EJE^\dagger`$ with E one native phase layer and J Clifford,
then $`L_E^\perp\subseteq\mathcal G(S)`$. Thus
$`\mathcal G(S)^\perp\subseteq L_E`$ for every such representation.
This intrinsic subspace nevertheless need not detect any non-Clifford
outer action.

Use three wires a,b,c, with a clean, and define

```math
D=R_{YZI}^\dagger R_{YIX}^\dagger R_{YZX},\qquad
U=T_cDH_a.
```

The three axes of D commute and are independent. A Clifford B applying

```
CX(c,b), CX(a,c), CX(b,a), H(a), S(a), H(c)
```

sends $`Z_a,Z_b,Z_c`$ to YZI,YIX,YZX, respectively. Therefore
$`B(T_a^\dagger T_b^\dagger T_c)B^\dagger=\zeta_8^{-1}D`$,
and U has two actual native T layers up to a common phase. Exact Pauli
multiplication gives

```math
D(XII)D^\dagger
=H_a\otimes\mathrm{CX}(b,c)
=\frac{XII+XIX+XZI-XZX+ZII+ZIX+ZZI-ZZX}{2\sqrt2}.
```

Consequently its clean image is

```math
S=UZ_aU^\dagger=H_a\otimes M,\qquad
M=T_c\mathrm{CX}(b,c)T_c^\dagger
=\mathrm{CX}(b,c)\mathrm{diag}(\zeta_8^{b-2bc}).
```

The diagonal factor is a weighted-cubic phase polynomial, so S belongs
to the third level of the Clifford hierarchy and maps every Pauli to a
Clifford. Hence $`\mathcal G(S)`$ is the entire Pauli label space and
$`\mathcal G(S)^\perp=0`$. But S is not Clifford, since

```math
S(IXI)S=IX(X+Y)/\sqrt2
```

is not a Pauli. The same U may be used on both sides of a shared-image
pair, giving zero discrepancy. This example refutes only the proposed
intrinsic-axis extraction; it is not a nonzero small-defect family.

## 6. A uniform one-sided Clifford comparison

There is a useful positive lemma that does not solve the general
two-sided compression. Let $`S=EJE^\dagger`$ be a traceless Hermitian
involution, with J a Hermitian Clifford and E a phase-polynomial layer,
and put $`Q=(I+S)/2`$. For fixed ring-valued Clifford representatives
G,H, either $`(G-H)Q=0`$ or its operator norm is at least an absolute
positive constant. One conservative choice is

```math
A_*=131070,\qquad C_*=1834,\qquad
\gamma=2^{-\frac{C_*}{2}\frac{A_*^{256}-1}{A_*-1}}.
```

To prove this, put $`W=H^\dagger G`$ and use normalized full trace.
For $`\varepsilon=\|(W-I)Q\|`$,

```math
0\le1-\Re\mathrm{tr}(W)-\Re\mathrm{tr}(WS)
=\|(W-I)Q\|_{2,\mathrm{tr}}^2\le\varepsilon^2.
```

The one-phase near-one argument of R75 applies to
$`\Re\mathrm{tr}(WEJE^\dagger)`$: its phase-polynomial and cosine
steps require a fixed eighth-root Clifford scalar, rather than
Hermiticity of W. Put $`\delta=(1-1/\sqrt2)/2`$. If
$`\Re\mathrm{tr}(W)\le\delta/2`$ and
$`\varepsilon\lt\sqrt{\delta/2}`$, then
$`\Re\mathrm{tr}(WS)\gt1-\delta`$. The gap forces WS=I, whence
W=S and $`(W-I)Q=0`$.

Otherwise $`|\mathrm{tr}(W)|\gt\delta/2`$. For
$`s=\mathrm{rank}(F_W-I)`$, Clifford trace orthogonality gives either
zero trace or $`|\mathrm{tr}(W)|=2^{-s/2}`$, so $`s\le7`$.
Indeed $`|\mathrm{Tr}(W)|^2`$ is the trace of conjugation by W on
the Pauli basis; its diagonal signs form a character of
$`\ker(F_W-I)`$ and sum to zero or $`2^{2q-s}`$ on q wires.
Applying the same identity to PW shows that W has exactly
$`2^s\le128`$ nonzero Pauli coefficients, of magnitude
$`2^{-s/2}`$. Their phases are eighth roots by the entry ring and
unit-circle lemma.

Expand W in

```math
(W-I)Q=(W-I+WS-S)/2.
```

Each product PS is
$`E(E^\dagger PE)JE^\dagger`$, with a Clifford boundary. W and I
are bare Clifford boundaries. Express every coefficient at the common
scale $`2^{-9/2}`$ using eighth-root phases. The three terms W,-I,-S
cost 16 summands each. The WS terms cost at most 128 summands; when a
half-integer power of two occurs, use
$`\sqrt2=\zeta_8+\zeta_8^{-1}`$. The total is even and at most 176.
Canceling pairs pad it to 256. R74's finite-boundary theorem with
$`K=256,b=9,m=0`$ gives the stated nonzero norm floor.

If both G,H preserve Q, their two-sided compressed difference equals
this one-sided difference. In general $`Q(G-H)Q`$ need not equal
$`(G-H)Q`$. No uniform gap for that unrestricted compression, or for
the needed four-Clifford sum, follows from this lemma.

## 7. Reopening conditions

One smaller center is genuinely preserved by the relative data action.
In the native correlated setting, suppose $`SiS=p`$ is a Pauli and both
i,p commute with both outer spaces. Put $`h=i+p`$ at the binary-label
level. This label lies in the image of each inner clean displacement
map. The correlated inner data maps satisfy
$`(F_K-I)(F_J-I)=0`$, because both have images in the same inner
isotropic axis span. They therefore fix h, and the outer layers commute
with h. Thus the actual data images obey

```math
AhA=\epsilon_A(h)h,\qquad BhB=\epsilon_B(h)h.
```

The Hermitian product representatives of these labels commute with S.
If the two signs differed, applying A and B to an occupied h-eigenvector
would give orthogonal unit vectors, with difference norm $`\sqrt2`$.
Below that threshold their sign characters agree, and BA preserves the
generated center. This does not imply that A and B separately preserve
its sectors. Its binary rank may also be zero: in the R77 example
$`S(ZI)S=-ZI`$, so $`h=i+p=0`$ while the data still mixes the original
center. Tensor copies make that displayed intersection arbitrarily large.

The common-center route needs a uniform estimate that handles the actual
data mixing, or an exact change of native representation that removes
it. A finite inner alphabet or clean-invariant center does not supply
that estimate. Section 3 makes representation optimization a concrete
possibility rather than an assumed invariance.

Any proposed clean-removal theorem must accommodate the determinant
obstruction in Section 4, specify additional arbitrary work if used,
and preserve the claimed number of phase layers. No native family with
nonzero occupied defect tending to zero, and no four-layer approximation
to the complete witness, was obtained from these examples.
