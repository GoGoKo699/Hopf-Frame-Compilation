# Four-layer operator rigidity and a zero-clean five-layer obstruction

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Three-layer obstruction](CLEAN_AWARE_DEPTH_OBSTRUCTION.md) · [Crossing-rank obstruction](FOUR_LAYER_CROSSING_RANK.md) · [Research index](../README.md)

This result, R75, gives a width-independent obstruction at four T
layers in the zero-clean, arbitrary-dirty model. Its premise is operator
norm closeness. It does not require a uniform gap for scalar four-layer
Pauli-transfer coefficients, and does not assert such a scalar gap.

Put

```math
\delta=\frac{1-1/\sqrt2}{2},\qquad
\rho=\sqrt{\delta/2}
=\frac{\sin(\pi/8)}{\sqrt2},\qquad
\tau=\frac\rho2
=\frac{\sin(\pi/8)}{2\sqrt2}\approx0.135299.
```

**Operator-rigidity theorem.** If A and B are Hermitian Pauli images
under circuits with at most two T layers, on the same physical space,
then

```math
A=B\quad\text{or}\quad\|A-B\|\ge\rho.
```

Consequently an at-most-four-layer circuit V satisfying
$`\|VPV^\dagger-P\|\lt\rho`$ for a Hermitian Pauli P must preserve
P exactly. In particular, if a target U preserves P and
$`\|V-e^{i\phi}U\|\lt\tau`$ for some common phase, then V
preserves P exactly.

For the complete-frame witness

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad n\ge1,
```

at most four T layers have full-space operator error at least tau,
even with arbitrary returned dirty width and a freely optimized common
phase. Hence every integer accuracy exponent $`L\ge3`$ requires at
least five T layers in the zero-clean model. The same real witness is
contained in the complex magnitude family by choosing zero leaf phases.
This is a lower bound on T-depth, not an optimal-depth theorem.
Optimal depth remains an open problem.

## 1. Two-layer Pauli images and commuting axis spaces

A native circuit has arbitrary Clifford blocks and physical layers of
T or T-dagger gates on distinct wires. Identity layers may be inserted.
Every two-layer Pauli image has the form

```math
A=E_AK_AE_A^\dagger,
```

where $`K_A`$ is a Hermitian Clifford and $`E_A`$ is a Clifford
conjugate of one physical T layer. Indeed, in
$`U=C_2D_2C_1D_1C_0`$, the operator
$`C_1D_1(C_0PC_0^\dagger)D_1^\dagger C_1^\dagger`$ is a
Hermitian Clifford: a single T layer sends each Pauli to a tensor product
of Hermitian one-qubit Cliffords. Conjugating this operator and D2 by C2
gives the stated form. The same holds for B.

Let $`L_A`$ be the binary span of the active signed Pauli axes of
$`E_A`$, and define $`L_B`$ in the same way. Both spaces are isotropic.
Every Hermitian Pauli R whose label belongs to $`L_A^\perp`$ commutes
with $`E_A`$. The use of these full Pauli centralizers requires no
Lagrangian extension, additional T gate, or initialized qubit.

Two elementary consequences will be used:

* For $`R\in L_A^\perp`$, the operator
  $`ARA=E_A(K_ARK_A^\dagger)E_A^\dagger`$ is a Hermitian Clifford,
  since the inner operator is a Pauli.
* For every Hermitian Pauli R, without a restriction to $`L_B^\perp`$,

  ```math
  BRB=E_BK'_RE_B^\dagger,\qquad
  K'_R=K_B(E_B^\dagger R E_B)K_B^\dagger
  ```

  has a Hermitian Clifford core $`K'_R`$. Here a Clifford-conjugated
  physical T layer still sends every Pauli to a Hermitian Clifford.

## 2. The one-phase gap with arbitrary Hermitian Clifford boundaries

**Lemma.** For arbitrary Hermitian Clifford unitaries G,H and a
Clifford-conjugated physical T layer E, the real number

```math
x=2^{-q}\mathrm{Tr}(G E H E^\dagger)
```

satisfies

```math
x=1\quad\text{or}\quad x\le1-\delta.
```

The boundaries need not themselves be images of Paulis under one T
layer. This slight generality is important in the application above.

The proof is precisely the weighted-phase argument in
[the three-layer obstruction, Section 4.1](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#41-a-sharp-gap-below-a-preserved-pauli).
For completeness, change Clifford coordinates to make E a physical
T layer and use normalized Choi vectors for G,H. They are stabilizer
states. Their computational supports are affine spaces of dimensions
a,b. If their supports are disjoint, x=0. Otherwise their intersection
has dimension k, and their flat amplitudes give

```math
x=2^{k-(a+b)/2}\mathbb E_u\zeta_8^{p(u)}
=2^{-s/2}\mathbb E_u\zeta_8^{p(u)},\qquad
s=a+b-2k\ge0.
```

Hermiticity fixes each boundary's otherwise arbitrary common phase to
an eighth root of unity. To see this directly, write all nonzero matrix
entries as a common phase times fourth roots of unity times the common
magnitude. A nonzero entry and its transposed conjugate then force the
square of that common phase to be a fourth root. The doubled phase word
$`E\otimes\overline E`$ contributes a mod-eight parity polynomial.
Together with stabilizer phases and affine restriction, it gives

```math
p(u)=c+\sum_i\ell_i u_i
+2\sum_{i\lt j}q_{ij}u_i u_j
+4\sum_{i\lt j\lt k}t_{ijk}u_i u_j u_k\pmod8.
```

If $`s\ge1`$, then $`x\le1/\sqrt2`$. If s=0 and p modulo two
is nonzero, it is affine and takes the value one on at least half the
inputs. Each odd residue contributes at least $`1-1/\sqrt2`$ to
$`1-\cos(\pi p/4)`$, giving loss at least delta. If that parity
vanishes, divide p by two: a nonzero next parity has degree at most two,
occurs on at least a quarter of inputs, and contributes loss at least
one there. If it too vanishes, the final parity p/4 has degree at most
three and, if nonzero, contributes loss two on at least one eighth of
inputs. These latter cases give loss at least 1/4. If all parities
vanish, p is identically zero and x=1. The elementary Boolean support
bound used here is the same one proved in the cited section.

## 3. Exact commutation from close conjugation actions

Assume $`\varepsilon=\|A-B\|\lt\rho`$. For every Hermitian
Pauli R,

```math
\|ARA-BRB\|\le2\varepsilon.
```

For $`R\in L_A^\perp`$, Section 1 puts the two operators in the class of
Section 2. Their normalized trace overlap y therefore obeys the gap
alternative. Since both are Hermitian unitaries,

```math
1-y
=\tfrac12\|ARA-BRB\|_{2,\mathrm{tr}}^2
\le2\varepsilon^2
\lt\delta.
```

Thus y=1 and $`ARA=BRB`$. The same argument with A,B interchanged
works for every $`R\in L_B^\perp`$. Put $`W=BA`$. Equality of these
conjugation actions means that W commutes with every Pauli in
$`L_A^\perp+L_B^\perp`$.

Expand W in the orthogonal Pauli basis. A coefficient on a label v can
survive these commutation conditions only if v is symplectically
orthogonal to $`L_A^\perp+L_B^\perp`$. Consequently

```math
W\in\mathrm{span}_{\mathbb C}\{P_v:v\in I\},\qquad
I=(L_A^\perp+L_B^\perp)^\perp=L_A\cap L_B.
```

The space I is isotropic. A single Clifford change of coordinates sends
it to a Z-type subspace, so W becomes diagonal. This exact commutant
conclusion has no dependence on the width or the crossing rank.

## 4. Discrete phases in the common commutant

Pauli images under exact Clifford+T circuits have matrix entries in

```math
\mathcal R=\mathbb Z[i,1/\sqrt2].
```

Their common circuit phases cancel in conjugation. Products and
Clifford conjugation preserve this entry ring, using the usual
H,S,CNOT representatives for Clifford operators. Therefore every
diagonal entry of the Clifford-diagonalized W belongs to this ring.
Unitarity makes every such entry have modulus one.

The unit-circle classification proved in
[the three-layer obstruction, Section 4.2](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#42-exact-dirty-return-and-the-explicit-zero-clean-floor)
states that the only modulus-one elements of this ring are the eighth
roots of unity. Its hypotheses apply entrywise here. This is not an
assertion that arbitrary eigenvalues of a Clifford+T matrix lie in the
entry ring: the exact Clifford diagonalization established in Section 3
is what makes the assertion available.

If W is not the identity, at least one diagonal entry is a nontrivial
eighth root. Hence

```math
\|W-I\|\ge|\zeta_8-1|=\sqrt{2-\sqrt2}.
```

But $`\|W-I\|=\|A-B\|\lt\rho\lt\sqrt{2-\sqrt2}`$.
Therefore W=I and A=B, proving operator rigidity.

## 5. Four-layer preservation and the zero-clean witness

Split any at-most-four-layer circuit as $`V=V_2V_1`$, with each half
having at most two T layers. For a Hermitian Pauli P, set

```math
A=V_2^\dagger P V_2,\qquad
B=V_1 P V_1^\dagger.
```

Both are two-layer Pauli images, including when either half is inverted.
Moreover

```math
\|A-B\|=\|P-VPV^\dagger\|.
```

If the latter norm is below rho, operator rigidity gives A=B and hence
$`VPV^\dagger=P`$. If U preserves P and
$`\|V-e^{i\phi}U\|\le\eta`$, the latter norm is at most 2 eta.
Thus $`\eta\lt\tau`$ forces exact preservation, simultaneously
for any number of target-preserved Paulis and without a union bound.

With no initialized clean work, approximation to
$`W_n\otimes I_{\rm dirty}`$ is a full-space operator condition.
The target preserves X and Z on every spectator logical and dirty wire,
and Y on the final logical wire. Below tau, V preserves all these
Paulis exactly. Commutation with their algebra gives

```math
V=e^{i\varphi}R_y(\alpha)\otimes I_{\rm spectators}.
```

Its transfer entries $`\cos(2\alpha),\sin(2\alpha)`$ lie in
$`\mathbb Z[1/\sqrt2]`$. The same unit-circle classification forces
$`\alpha\in(\pi/8)\mathbb Z`$ modulo a common phase. As calculated
in the cited Section 4.2, the nearest allowed rotation to this witness
has phase-optimized operator distance

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\gt\tau.
```

This contradicts any assumed error below tau. The result is therefore
the lower bound tau, not a claim that the exact four-layer optimum is
c0. Since $`2^{-L}\le1/8\lt\tau`$ for every integer $`L\ge3`$,
at least five T layers are required at those accuracies.

## 6. One clean qubit: terminal sectors only

For one initialized clean qubit, write

```math
V=\begin{pmatrix}A_0&B_0\\C_0&D_0\end{pmatrix}.
```

An initialized-isometry error at most eta against a clean-returning
target gives $`\|C_0\|\le\eta`$. Unitarity gives
$`\|B_0\|=\|C_0\|`$, because
$`B_0B_0^\dagger=I-A_0A_0^\dagger`$ and
$`C_0^\dagger C_0=I-A_0^\dagger A_0`$ have the same eigenvalues.
Consequently

```math
\|VZ_CV^\dagger-Z_C\|=\|[V,Z_C]\|\le2\eta.
```

At at most four T layers and $`\eta\lt\tau`$, the preservation
theorem forces $`[V,Z_C]=0`$: both terminal clean sectors are preserved
exactly. This does not force exact preservation of logical or dirty
Paulis on the occupied sector. The unoccupied sector has no specified
logical action, and the preceding zero-clean argument cannot be
applied to it. No one-clean five-layer lower bound follows here.

The [midpoint theorem R76](ONE_CLEAN_MIDPOINT_RIGIDITY.md) proves such
a bound under a further hypothesis on a two-plus-two split. A Pauli
encoded clean midpoint gives the same error floor tau; a Hermitian
Clifford midpoint gives a smaller explicit positive threshold, with no
crossing-rank bound. R76 also shows that the operator-rigidity proof
extends to commuting phase-polynomial layers with dependent axes. The
unrestricted one-clean case still permits a non-Clifford midpoint.

## 7. What remains open for scalar overlaps

Let $`x=2^{-q}\mathrm{Tr}(AB)`$ for two two-layer Hermitian
Pauli images. The existence of a dimension-independent gap below one
for all unequal A,B remains open. The operator-rigidity theorem does
not establish that gap.

There is a useful exact reduction. Sections 1–3 work also with
$`\varepsilon=\|A-B\|_{2,\mathrm{tr}}`$, because multiplication by
unitaries preserves this norm. Thus whenever

```math
x>1-\delta/4,
```

W=BA already belongs to the common Pauli algebra of
$`I=L_A\cap L_B`$ and is Clifford-diagonal with eighth-root entries.
If $`d_I=\dim I`$, its joint eigensectors all have the same fraction
$`2^{-d_I}`$ of the physical dimension. Therefore

```math
x=1\quad\text{or}\quad
1-x\ge(1-1/\sqrt2)\,2^{-d_I}
```

within this near-one regime. A nontrivial phase confined to rare joint
eigensectors could have a small normalized-trace deficit, while its
operator distance from identity is still at least
$`\sqrt{2-\sqrt2}`$. Whether such increasingly rare central phases
actually arise from native pairs of two-layer Pauli images is not
settled by this reduction. It is a remaining realization question,
not a counterexample asserted here.

This explains why full-space operator approximation admits the
five-layer obstruction even while the scalar four-layer closure
question and positive-clean frontier remain unresolved.
