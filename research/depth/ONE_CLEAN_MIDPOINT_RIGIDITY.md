# One-clean midpoint rigidity

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Operator rigidity](FOUR_LAYER_OPERATOR_RIGIDITY.md) · [Crossing rank](FOUR_LAYER_CROSSING_RANK.md) · [Research index](../README.md)

This result, R76, gives one-clean obstructions with an explicit condition
on the encoded clean observable at a two-layer midpoint. Dirty width and
middle crossing rank are unrestricted. The midpoint condition is not a
consequence of terminal clean return.

Let a native circuit of T-depth at most four be split as
$`V=V_2V_1`$, with each half of T-depth at most two. Write

```math
\Pi=(I+Z_C)/2,\qquad S=V_1Z_CV_1^\dagger.
```

For the complete-frame witness

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad n\ge1,
```

the following initialized-isometry error floors hold, including arbitrary
dirty inputs and a freely optimized common phase.

| Hypothesis on a four-layer realization | Error floor |
|---|---|
| It admits a split with S a Pauli | $`\tau=\sin(\pi/8)/(2\sqrt2)\approx0.135299`$ |
| It admits a split with S a Hermitian Clifford | $`\eta_C=2^{-\mathcal M/2-2}`$, where $`\mathcal M=73(126^8-1)/125`$ |
| Its occupied block has finite projective order after exact clean return | $`1/196608`$ |

The first row requires at least five layers for integer $`L\ge3`$.
The second gives a positive width-independent threshold; its explicit
constant is very conservative. The third row is a separate spectral
obstruction, and applies to a finite-projective-order full circuit at
$`L\ge18`$, because R75 first supplies exact clean return. It does not
require a Clifford midpoint.

An unrestricted four-layer one-clean approximation could have a
non-Clifford midpoint and an occupied block of infinite projective order.
These theorems do not exclude that possibility.
Optimal depth remains an open problem.

## 1. Phase-polynomial layers and operator rigidity

A phase-polynomial layer is a Clifford conjugate of a diagonal operation
$`\mathrm{diag}(\zeta_8^{p(x)})`$, where, modulo eight,

```math
p(x)=c+\sum_i\ell_i x_i
+2\sum_{i\lt j}q_{ij}x_ix_j
+4\sum_{i\lt j\lt k}t_{ijk}x_ix_jx_k.
```

These are diagonal third-level Clifford-hierarchy operations in the chosen
coordinates. Equivalently, they are products of commuting Pauli phase
rotations with integer multiples of the T angle. The axes may be
dependent. This definition makes no claim that the operation has physical
T-depth one.

The R75 operator-rigidity proof holds unchanged with these layers:

* The finite difference of p along a binary translation is a diagonal
  Clifford phase. Hence a layer sends every Pauli to a Hermitian Clifford.
  A two-layer Pauli image therefore has the form $`EKE^\dagger`$ with
  K a Hermitian Clifford.
* E lies in an isotropic Pauli algebra L, so every Pauli in
  $`L^\perp`$ commutes with E. A full diagonal Lagrangian may also be
  used. Independence of a displayed list of phase axes is unnecessary.
* The doubled Choi phase is still weighted cubic, including after affine
  restriction. R75's one-phase gap remains
  $`\delta=(1-1/\sqrt2)/2`$.
* The close-conjugation argument puts the relative unitary into a common
  isotropic Pauli algebra. Its exact Clifford diagonalization has entries
  in $`\mathcal R=\mathbb Z[i,1/\sqrt2]`$. R73's unit-circle lemma
  makes those entries eighth roots, exactly as in R75.

Consequently distinct Pauli images under two phase-polynomial layers
have operator-norm distance at least

```math
\rho=\sqrt{\delta/2}=\frac{\sin(\pi/8)}{\sqrt2}.
```

A circuit with at most four such layers satisfies
$`\|[U,P]\|\lt\rho\Rightarrow[U,P]=0`$ for every Hermitian
Pauli P. Approximation to a Pauli-preserving target below
$`\tau=\rho/2`$ therefore forces exact preservation.

## 2. Two native layers with exact clean return

**Descent theorem.** A native circuit W of T-depth at most two that
preserves a physical Pauli P induces, on either eigenspace of P, at most
two phase-polynomial layers with Clifford blocks, after a Clifford
identification of that eigenspace with the logical space.

### Diagonal Clifford extraction

Suppose D is a phase-polynomial operation on an isotropic space L and
$`DRD^\dagger=R'`$ for Hermitian Paulis R,R'. There is a diagonal
Clifford Q on L with the same action. To prove this, make L Z-type.
The two Paulis have the same translation v on these coordinates. If
v=0, take Q=I. Otherwise change binary coordinates so $`v=e_1`$.
The required diagonal derivative is a Pauli character

```math
i^t(-1)^{\sum_j b_jx_j},\qquad t\bmod2=b_1.
```

The parity condition follows by applying the translation twice. For the
derivative convention $`q(x+e_1)-q(x)`$, the Clifford polynomial

```math
q(x)=2t x_1+4x_1\sum_{j\gt1}b_jx_j\pmod8
```

has exactly this derivative. The opposite convention uses its inverse.
Thus Q exists, and $`D_0=DQ^\dagger`$ commutes with R'. If R'
translates a coordinate of L, this commutation makes D0 independent of
that coordinate. Its remaining phase support commutes with R'. Fixing
a diagonal stabilizer or removing this translated coordinate preserves
weighted cubic degree. The resulting operation on the stabilizer
quotient is again a phase-polynomial layer.

### The essential axes are intrinsic

For a conjugated physical T/T-dagger layer E and a Pauli R, split
$`E=E_aE_s`$ according to axes anticommuting or commuting with R.
The silent part Es commutes with R, and
$`M=ERE^\dagger=E_aRE_a^\dagger`$ is a Hermitian Clifford.
If the essential axis labels are $`a_i`$, its binary symplectic map is

```math
F_M(v)=v+\sum_i\langle a_i,v\rangle a_i,\qquad
\mathrm{im}(F_M-I)=\mathrm{span}\{a_i\}.
```

The last equality uses independence of the physical layer's axes.
Signs and T inverses leave the binary transvections unchanged. Hence
the essential span is intrinsic to M. The silent layer commutes with
every axis in that span.

### Descent proof

Split $`W=U_2^\dagger U_1`$ with each Uj having at most one native
T layer, and write

```math
U_j=E_jC_j,\qquad R_j=C_jPC_j^\dagger.
```

Preservation of P gives the same Hermitian Clifford image
$`M=U_1PU_1^\dagger=U_2PU_2^\dagger`$. Both essential spans are
therefore the same isotropic space $`L=\mathrm{im}(F_M-I)`$.
The operation

```math
D=E_{2,a}^\dagger E_{1,a}
```

is a phase-polynomial layer on L and maps R1 to R2. Extract a diagonal
Clifford Q with that action and put $`D_0=DQ^\dagger`$. Both silent
layers commute with L, Q and D0. The first silent layer also commutes
with R2 because it commutes with R1 and Q. Reordering these commuting
factors gives

```math
W=C_2^\dagger(D_0E_{2,s}^\dagger)E_{1,s}QC_1.
```

Both phase factors preserve R2. The first is a single phase-polynomial
operation: D0 is supported on the R2-commuting part of L, and the
silent axes commute with L. Each factor descends to a phase-polynomial
layer on the R2 code. Finally $`C_2^\dagger QC_1`$ is a Clifford
preserving P, so its restriction is a logical Clifford. This proves
the theorem. The case L=0 is included. Dependencies created by the
stabilizer quotient are allowed in the descended phase layers.

### Pauli midpoint consequence

If $`V=V_2V_1`$ preserves the terminal clean Pauli and
$`S=V_1Z_CV_1^\dagger`$ is Pauli, choose a Clifford B sending S
to $`Z_C`$. Both $`BV_1`$ and $`V_2B^\dagger`$ preserve
$`Z_C`$ and retain depth at most two. Their occupied restrictions
have at most four phase-polynomial layers in total. Section 1 therefore
gives

```math
\|[V,\Pi P]\|\lt\rho\quad\Longrightarrow\quad[V,\Pi P]=0
```

for every data Pauli P. The clean qubit need not be a separate wire at
the original midpoint; a Clifford encoding of it is allowed.

## 3. A Clifford midpoint with unrestricted crossing rank

**Compressed-image theorem.** Let A,B,S be Hermitian involutions, with
S traceless and Clifford and $`[A,S]=[B,S]=0`$. Suppose

```math
A=E_AK_AE_A^\dagger,\quad S=E_AJ_AE_A^\dagger,\qquad
B=E_BK_BE_B^\dagger,\quad S=E_BJ_BE_B^\dagger,
```

where the four cores are Hermitian Cliffords and EA,EB are
phase-polynomial layers. Put $`Q=(I+S)/2`$ and

```math
\mathcal M=73\frac{126^8-1}{125}=37100281773115459,\qquad
\kappa=2^{-\mathcal M/2-1}.
```

Then $`\|(A-B)Q\|\lt\kappa`$ implies $`(A-B)Q=0`$.
Native two-layer encoders of the same clean Pauli supply these correlated
representations. No crossing-rank bound is imposed.

### Eight Clifford boundaries

Let R be any Pauli commuting with the diagonal axis space of EA.
Then $`G=ARA`$ is a Hermitian Clifford. Independently,
$`BRB=E_BK_RE_B^\dagger`$ has a Hermitian Clifford core, since
$`E_B^\dagger RE_B`$ is Clifford. The compressed discrepancy is

```math
Q(ARA-BRB)Q
=\tfrac14(G+SG+GS+SGS)
-\tfrac14 E_B(K_R+J_BK_R+K_RJ_B+J_BK_RJ_B)E_B^\dagger.
```

All eight boundaries are Cliffords. R need not commute with S, and the
individual products need not be Hermitian. Their fixed ring-valued
Clifford representatives have flat affine Choi supports with eighth-root
common phases. After diagonalizing EB, the right boundaries acquire one
weighted-cubic phase.

The [finite-boundary form of R74](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
applies with K=8 summands, coefficient magnitude
$`1/4=2^{-b/2}`$ for b=4, and m=0 active Hadamard coordinates.
Its recurrence is

```math
t_0=0,\qquad t_{j+1}=126t_j+73\quad(0\le j\lt8),\qquad
t_8=\mathcal M.
```

Thus a nonzero discrepancy has normalized Hilbert--Schmidt squared norm
at least $`2^{-\mathcal M}`$. This is a gap for a squared vector
norm, not an asserted gap at an interior scalar transfer value.

If $`\epsilon=\|(A-B)Q\|`$, the discrepancy has operator norm
at most $`2\epsilon`$. Indeed, A and B restrict to unitaries on QH,
and $`\|QRQ\|\le1`$. Below kappa the discrepancy must vanish
for every tested R. No union bound is needed.

### Lifting the compressed commutant

Let Y be the unitary BA restricted to QH. The frozen identities imply
that Y commutes with all compressed operators QRQ for
$`R\in L_A^\perp`$. These Paulis span an algebra
$`\mathcal N`$ whose commutant is the abelian Pauli algebra on LA.
After a Clifford change of coordinates,

```math
\mathcal H=\bigoplus_z\mathcal H_z,\qquad
\mathcal N=\bigoplus_z\mathcal B(\mathcal H_z).
```

Let J embed QH into H, let Pz project onto Hz, and put $`J_z=P_zJ`$.
The operators $`J_z^\dagger TJ_z`$, as T ranges over
$`\mathcal B(\mathcal H_z)`$, span all matrices on
$`\mathcal K_z=\mathrm{range}(J_z^\dagger)`$, extended by zero
on its orthogonal complement. This follows by inverting the restriction
of Jz to Kz onto its range. Commutation therefore makes Y a scalar
$`\lambda_z`$ on each occupied Kz. These spaces span QH and
$`|\lambda_z|=1`$.

Define $`D_0=\bigoplus_z\lambda_z I_{\mathcal H_z}`$, choosing
one on unused blocks. The scalar actions give
$`D_0J=JY`$ and $`J^\dagger D_0=YJ^\dagger`$. Consequently

```math
X:=JYJ^\dagger=BAQ=D_0Q=QD_0.
```

This elementary argument lifts the compressed commutant. It does not
infer that X commutes with the original uncompressed test Paulis.

### Discrete occupied phases

In these Clifford coordinates X,Q have entries in
$`\mathcal R=\mathbb Z[i,1/\sqrt2]`$. For an occupied block z,
choose a row i with a nonzero row of Q. Then
$`X_{ij}=\lambda_zQ_{ij}`$. If that row of S has a nonzero
off-diagonal entry, $`Q_{ij}=S_{ij}/2`$ is a unit of R: every
nonzero Clifford entry is an eighth root times an integer power of
$`1/\sqrt2`$. Hence $`\lambda_z\in\mathcal R`$. Otherwise
the occupied row of Q is diagonal and $`Q_{ii}=1`$, with the same
conclusion.

The unit-circle lemma makes every occupied lambda an eighth root.
A nonidentity occupied scalar contributes at least
$`\sqrt{2-\sqrt2}`$ to $`\|X-Q\|=\epsilon`$, contradicting
$`\epsilon\lt\kappa`$. Thus X=Q and the theorem follows.

## 4. Initialized error and the frame witness

Let J append the clean zero and suppose
$`\|VJ-e^{i\phi}JU\|\le\eta`$ for the target frame tensored
with dirty identity. R75 forces $`[V,Z_C]=0`$ whenever
$`\eta\lt\tau`$. For any two-plus-two split this gives the shared
midpoint

```math
S=V_1Z_CV_1^\dagger=V_2^\dagger Z_CV_2.
```

For a target-preserved data Pauli P, set
$`A=V_1PV_1^\dagger`$ and $`B=V_2^\dagger PV_2`$.
The initialized-isometry estimate gives

```math
\|(A-B)(I+S)/2\|
=\|V\Pi PV^\dagger-\Pi P\|\le2\eta.
```

If S is Clifford and $`\eta\lt\kappa/2`$, Section 3 freezes every
required initialized Pauli observable. If S is Pauli, Section 2 gives
the stronger threshold tau. Spectator and dirty X,Z and the active Y
force the occupied block to be
$`e^{i\alpha_0}R_y(\alpha)\otimes I`$. Its transfer entries
$`\cos(2\alpha),\sin(2\alpha)`$ lie in
$`\mathbb Z[1/\sqrt2]`$, so the unit-circle lemma forces
$`\alpha\in(\pi/8)\mathbb Z`$ modulo common phase. The nearest
such rotation to the witness has distance

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\gt\tau.
```

This proves the first two error floors. Every decomposition with the
stated midpoint type is excluded below its threshold. Accordingly, an
approximating four-layer sequence at arbitrarily high precision must
eventually have a non-Clifford clean midpoint in every two-plus-two
decomposition. This is a necessary condition, not a construction.

## 5. A complementary finite-projective-order obstruction

The following argument does not use a midpoint condition or a T-depth
bound. Suppose an exactly returned occupied unitary A has entries in
$`K=\mathbb Q(\zeta_8)`$ and finite projective order: some positive
power of A is scalar. Then

```math
\inf_\phi\|A-e^{i\phi}(R_y(\theta)\otimes I)\|
\ge\frac1{196608}.
```

To prove this, put $`B=A\otimes A^\dagger`$. B has finite order,
entries in K, and distance at most $`2\eta`$ from
$`(R_y(\theta)\otimes I)\otimes(R_y(\theta)^\dagger\otimes I)`$
if the claimed error is eta. Common phase cancels. Every eigenvalue mu
of B is a root of unity within $`2\eta`$ of
$`\{1,e^{2i\theta},e^{-2i\theta}\}`$. This spectral inclusion
follows by applying the target resolvent to an eigenvector of B.

Use the integer polynomial

```math
p(z)=(z-1)(3z^2-2z+3)=3z^3-5z^2+5z-3.
```

On the unit circle, $`|p(z)|\le16`$ and
$`|p(z)-p(w)|\le24|z-w|`$. Hence $`|p(\mu)|\le48\eta`$
for each eigenvalue and each member of its Galois orbit over K: the
characteristic polynomial of B has coefficients in K.

For $`F=\mathbb Q(\mu)`$, the K-orbit occupies at least one
quarter of the rational conjugates, because
$`[F\cap K:\mathbb Q]\le[K:\mathbb Q]=4`$. If mu is not one,
p(mu) is a nonzero algebraic integer. The quadratic roots of p are not
roots of unity: their irreducible monic quadratic over the rationals has
trace 2/3, so they are not algebraic integers. The nonzero integer norm
of p(mu), using the bounds on its conjugates, now gives

```math
48\eta\ge16^{-3},\qquad \eta\ge\frac1{48\cdot16^3}.
```

If all eigenvalues of B are one, its unitarity gives B=I and A is
scalar. Its distance from the target is at least
$`2\sin(\theta/2)`$, which is larger still.

This excludes, for example, exactly returning conjugated-Clifford
families $`V=U^\dagger CU`$ with C Clifford, regardless of the
complexity of their intermediate code. For an at-most-four-layer full
circuit of finite projective order, approximation below tau first forces
clean return; the occupied block then inherits finite projective order.
Since $`2^{-18}\lt1/196608`$, this family is excluded for
$`L\ge18`$. No assertion of finite projective order is made for a
generic Clifford+T circuit.

## 6. The remaining one-clean question

The sufficient unrestricted target is a width-independent c such that

```math
D_T(V)\le4,\quad[V,\Pi]=0,\quad\|[V,\Pi P]\|\lt c
\quad\Longrightarrow\quad[V,\Pi P]=0.
```

Equivalently, set $`J_P=I-\Pi+\Pi P`$, a controlled-Pauli
Clifford involution. Its two split images have the same projector image,
and their difference is exactly the compressed Pauli-image difference.
Jp has symplectic displacement rank at most two. This does not by itself
bound the crossing of the two outer phase layers.

The lifting and phase-discreteness steps in Section 3 already extend to
a general native two-layer midpoint S: in coordinates diagonalizing EA,
$`S=E_AJ_AE_A^\dagger`$ has nonzero entries equal to flat Clifford
entries times eighth roots. The remaining step is freezing the compressed
conjugations. When S is non-Clifford, the left terms SG and GS in the
eight-term discrepancy are no longer Clifford boundaries. Their second
phase layer can have unbounded crossing rank, so the zero-active-coordinate
finite-boundary argument does not apply. Pairwise R75 separation does
not establish a gap for this correlated sum.

The [inner clean-spread extension R77](ONE_CLEAN_INNER_SPREAD.md)
does freeze these conjugations when either encoder's first layer has
bounded clean-active spread t. The inner clean Clifford has a flat
$`2^t`$-term Pauli expansion, turning SG and GS into finite sums of
bare Clifford boundaries; G and SGS are already Clifford. This gives
an explicit positive gap even when S is non-Clifford. The unrestricted
step remains open when that spread grows in both encoders.

Terminal clean return alone supplies no Pauli or Clifford decoder for
that intermediate code. Scalar overlaps approaching one, or a rare
sector with an order-one defect, also do not refute the operator-norm
conjecture. A counterexample must be an actual native four-layer family
with exact clean return and a nonzero occupied commutator tending to zero.

There is an alternative exact certificate for an individual compressed
comparison. Suppose a Hermitian unitary C anticommutes with S and has
Pauli sign actions $`CAC=\alpha A`$, $`CBC=\beta B`$, where
$`\alpha,\beta\in\{1,-1\}`$. Choose $`k\in\{0,1\}`$ with
$`(-1)^k=\alpha\beta`$. Then C transports the two clean-sector
blocks of $`A-S^kB`$ into one another with a common sign, so

```math
\|A-S^kB\|=\|(A-B)(I+S)/2\|.
```

Both A and $`S^kB`$ are native two-layer Pauli images. R75 therefore
gives exact occupied equality below rho. Each encoder of a non-Clifford
S separately supplies a Hermitian Clifford flip: in the inner phase
coordinates the input clean Pauli cannot lie in the active isotropic
span, since that would leave only a one-layer Pauli image. Symplectic
duality supplies an anticommuting Pauli in that span's centralizer.
Its inner layer cancels, and its outer image is Clifford. A common flip
with both sign actions need not exist:
[R77's exact two-qubit example](ONE_CLEAN_INNER_SPREAD.md#5-a-common-compatible-flip-need-not-exist)
has unequal absolute product traces on the two sectors, excluding even
an arbitrary sector-exchanging unitary with those actions. The sufficient
certificate remains valid when a common flip is available.

## 7. Verification and dependency boundary

The proofs use R73's ring arithmetic, R74's elementary finite-boundary
argument, and R75's one-phase gap. They do not use the polynomial-limit
premise retained in R73/R74's alternative scalar proofs. The compressed
commutant lift and the finite-projective-order obstruction are proved
directly above.

[Finite algebra audits](../../tests/test_one_clean_midpoint_rigidity.py)
check phase derivatives, native clean-return examples, the distinction
between compressed and full commutants, and the cluster constants. They
do not infer a width-independent theorem from finite instances. The
main one-clean compiler and its count/depth upper bounds are unchanged.
