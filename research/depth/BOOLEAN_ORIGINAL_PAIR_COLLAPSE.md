# Complete Booleanity collapses the original physical pair count

[Occupied transported center](OCCUPIED_BOOLEAN_PAIR_CENTER.md) ·
[Literal pair repair](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[All-active diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md) ·
[Earlier positive-singleton recognition](POSITIVE_SINGLETON_PAIR_RECOGNITION.md)

The scalar branch of the original pairs supplies a full Boolean
one-band comparison on a particularly safe stabilizer code. Unlike an
arbitrary central fiber, this code preserves independence of the actual
output-silent list. Boolean reality then removes its entire dressing
of a missing active direction. Comparing the other pair branches on
the full logical space forces every pair difference to equal the sum
of all actual active axes. Original middle independence permits only
one such pair.

## 1. Exact statement and operator contract

Use an actual one-clean three-native-band circuit

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,\qquad
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger.
```

Every native band is an independent commuting list of signed Hermitian
Pauli T/T-dagger axes, with its literal orientations and scalar phases.
Assume exact clean return and complete occupied action

```math
G=e^{i\theta}\mathrm{diag}((-1)^{g(x)})
```

on the entire logical and arbitrary dirty register, including arbitrary
reference inputs. The common scalar is retained. The unused clean
block is unrestricted. Let m count the genuine pairs in the original
physical middle list under $`J=E_AP_AE_A^\dagger`$, not logical
orbits produced by a later restriction or change of basis.

**Original-pair collapse theorem.** If $`m\ge1`$, then:

1. $`m=1`$.
2. Its difference $`r=b+d`$ is the sum of every actual clean-active
   axis at each original endpoint.
3. After literal orientation repair, every original middle singleton
   is positive. The corrected clean images agree and have one common
   active rank k, which is even and at least two.

Consequently the existing all-active diagonal reduction gives

| Paired branch | Boolean sign degree | Complete Boolean XOR degree |
| --- | --- | --- |
| Arbitrary active ranks before recognition | $`\deg g\le6`$ | $`\deg f\le5`$ |
| Resulting common active rank $`k\ge4`$ | $`\deg g\le5`$ | $`\deg f\le4`$ |

There is no bound on dirty width, original singleton count, endpoint
silent lists or their coupling ranks. The no-original-pair branch keeps
its separate [no-pair theorem](GROWING_DISPLACEMENT_NATIVE_ORBITS.md);
it is not assigned the six/five bound by this result.

## 2. Literal repair and the genuine pair parities

Write the actual pair relations and orientations as

```math
JB_jJ=s_jD_j,\qquad JD_jJ=s_jB_j,\qquad
\epsilon_j,\delta_j\in\{1,-1\},\quad s_j\in\{1,-1\}.
```

Keep the exact repair from the original insertion identity:

```math
\begin{aligned}
F_j&=T_{B_j}^{\epsilon_j}T_{D_j}^{s_j\epsilon_j},&
C_j&=T_{D_j}^{\delta_j-s_j\epsilon_j},\\
F&=\prod_jF_j,& C_{\rm rep}&=\prod_jC_j,\\
E_B&=C_{\rm rep}F E_-E_+.
\end{aligned}
```

Here $`E_-`$ and $`E_+`$ retain the original negative and positive
singleton axes and their original orientations. Each $`C_j`$ is an
actual Clifford, and every displayed middle factor commutes with the
others. The signed-axis identity $`T_{-P}^a=\omega^aT_P^{-a}`$,
where $`\omega=e^{i\pi/4}`$, proves $`[F_j,J]=0`$ literally;
the scalar factors for a negative pair sign cancel.

Absorb only $`C_{\rm rep}`$ into the output encoder and conjugate
the actual output list. Use corrected C-side notation until Section 7.
The complete physical operator is unchanged, and

```math
V=K_C^\dagger E_C^\dagger F E_-E_+E_AK_A,\qquad
S:=E_CP_CE_C^\dagger=E_-JE_-^\dagger.                 \qquad(1)
```

Both endpoint active sets are nonempty when m is positive. Denote their
actual spans by $`L_A,L_C`$; equality of these spans is not assumed.
Every singleton commutes $`L_A`$: its fixed binary label under J
and independence of the A-active list force each active-axis pairing
to vanish.

The [general occupied-center theorem](OCCUPIED_BOOLEAN_PAIR_CENTER.md#2-a-genuine-center-at-bounded-original-pair-count)
applies before any all-active conclusion. Define

```math
\Gamma_j=s_jB_jD_j,\quad r_j=b_j+d_j,\quad
\Pi_{j,\pm}=(I\pm\Gamma_j)/2,\quad
R_{\rm pair}=\mathrm{span}\{r_1,\ldots,r_m\}.
```

The $`r_j`$ are independent, belong to $`L_A\cap L_C`$,
and commute both clean Paulis and every actual native axis. The signed
operators $`\Gamma_j`$ are transported exactly between endpoints.
The repair fixes each one pointwise. Their literal pair formulas are

```math
F_j=\omega^{-\epsilon_j\eta_j}
 \left(\Pi_{j,+}S_{B_j}^{\epsilon_j}
                 +\Pi_{j,-}\omega^{\epsilon_j}I\right),
\qquad \eta_j=(1-s_j)/2.                              \qquad(2)
```

In particular all pair factors are scalar on their joint minus sector.
No pair factor is deleted from the complete word.

## 3. One-band Boolean reality on a full logical space

We use the following explicit form of the
[Boolean-reality lemma](COMPUTATIONAL_PHASE_DERIVATIVES.md#7-boolean-reality-gives-a-shared-computational-center).
Let D be a complete Boolean diagonal unitary up to one common scalar,
let a,c be Hermitian Paulis, and let B be an actual independent
commuting native T/T-dagger band on that same full register. If

```math
DaD^\dagger=B^\dagger cB,                             \qquad(3)
```

then every B axis commutes c, and the right side equals c exactly.

Indeed the right side is a Hermitian Clifford. The left side is a
Pauli times a diagonal Boolean derivative, so its Clifford displacement
image lies in the computational Z space. For a one-band conjugation,
the displacement image is precisely the span of its axes active on c:
their independence and the nondegeneracy of the ambient symplectic
pairing make the coefficient pairing map surjective. Thus every active
axis is a computational Z Pauli.

The two Pauli cores a,c now have the same computational X mask. Removing
a from (3) expresses the literal Boolean derivative as a scalar, a
computational Z Pauli, and a product of odd S powers on the active
axes. Squaring gives identity. The binary Pauli label of that square
is the sum of those active axis labels; the scalar and the Z core
cannot change it. Independence forces the active list to be empty.
This uses Boolean reality, not merely monomiality of the Clifford.

## 4. The scalar-sector code retains native independence

Take any physical label

```math
t\in L_A\cap p_A^\perp
```

and a signed Hermitian representative $`P_t`$. It commutes the
complete A list, every singleton and J. Equation (1) therefore gives
$`[P_t,S]=0`$, without requiring $`[P_t,F]=0`$.

For the corrected independent output-active list, binary fixation
under S gives

```math
0=(F_S-I)t
 =\sum_{c\in C_{\rm active}}\langle c,t\rangle c.
```

Independence forces $`\langle c,t\rangle=0`$ for each actual
C-active axis. Their phase rotations commute $`P_t`$ literally;
the remaining signed commutation with S is therefore its commutation
with $`P_C`$. Hence $`[P_t,P_C]=0`$ too. Only output-silent
axes could still dress $`P_t`$.

Use the existing shared computational-center normalization for the
genuine transported center. Its extracted quadratic Boolean phase
and fixed Clifford coordinate changes retain a complete Boolean
diagonal target. No separate frame is chosen on different sectors.
The two endpoint images of every $`\Gamma_j`$ are then the same
signed computational-Z Pauli, denoted $`\Gamma_{0,j}`$: their
actual transport through the diagonal target forces signed equality.

Let $`t_A,t_C`$ be the legitimate positive-clean restrictions of
$`K_A^\dagger P_tK_A,K_C^\dagger P_tK_C`$, respectively.
They are Hermitian logical Paulis and commute all $`\Gamma_{0,j}`$.
Write $`C_s`$ for the actual C-silent band pulled through the output
encoder and restricted to the clean code. Define the joint projector

```math
\Pi_{\boldsymbol -}=\prod_{j=1}^m(I-\Gamma_{0,j})/2.
```

First transport $`P_t`$ through the actual full physical word.
It commutes A and the singleton band; (2) makes every pair scalar
on the joint minus sector, and C-active commutes $`P_t`$.
Only then restrict the resulting clean-preserving equality. This gives

```math
(Gt_AG^\dagger)\Pi_{\boldsymbol -}
 =(C_s^\dagger t_CC_s)\Pi_{\boldsymbol -}.             \qquad(4)
```

The necessary native quotient is injective. If a sum of actual
C-silent labels satisfies

```math
\sum_{c\in I\subseteq C_{\rm silent}}c
      =\alpha p_C+r,\qquad r\in R_{\rm pair},
```

pairing with any actual C-active axis gives $`\alpha=0`$:
all C axes and r commute that axis, whereas $`p_C`$ anticommutes
with it. Now $`r\in L_C`$, and independence of the full actual
C list makes both the silent sum and r zero. Thus

```math
\mathrm{span}(C_{\rm silent})
 \cap\bigl(\langle p_C\rangle+R_{\rm pair}\bigr)=\{0\}.
                                                               \qquad(5)
```

Consequently restricting the C-silent list to the positive clean code
and the joint minus pair sector creates neither an identity axis nor
any relation. Its actual signed T/T-dagger factors remain one native
band, including all scalar phases.

Every joint pair sector is nonempty on the positive clean code.
The $`r_j`$ are independent; an actual active axis commutes all of
them and anticommutes the clean Pauli, so the clean label is independent
of their span as well. This is a family of m+1 independent commuting
stabilizer constraints, with every sign assignment allowed.

Use one affine computational coordinate change to turn the independent
signed $`\Gamma_{0,j}`$ into coordinate Z operators, then fix their
minus eigenvalues in (4). The restricted G is a Boolean diagonal
operator on the entire remaining register, not merely on some selected
vectors. Both restricted endpoint operators are Paulis, and (5) gives
an independent native band on the same register. Section 3 applies to
this full-register equality and makes every projected C-silent axis
commute the projected $`t_C`$.

All these original physical operators preserve the clean and pair
constraints. Their commutation signs therefore survive restriction;
an anticommuting pair of unitaries cannot become commuting on a nonzero
invariant code. Hence $`P_t`$ commutes **every actual C-silent axis
on the full physical space**. Together with the preceding calculation,
it commutes the entire output band.

This is a specific safe quotient, not a claim that restriction by the
entire transported center preserves native independence. The restriction
is a proof step and supplies no initialized helper or implementation.

## 5. Boolean derivatives cannot change commutation between pair sectors

The previous section applies to every $`t\in L_A\cap p_A^\perp`$.
Put $`T=Gt_AG^\dagger`$. Since $`P_t`$ now commutes both
complete endpoint bands and every singleton, the literal full-space
branch formulas give

```math
T\Pi_{\boldsymbol -}=t_C\Pi_{\boldsymbol -}.          \qquad(6)
```

For a fixed j, let $`\Pi^{(j)}`$ denote the sector with
$`\Gamma_{0,j}=+1`$ and every other pair parity minus. If
$`\langle b_j,t\rangle=1`$, the only non-scalar pair action
there is a quarter-turn about $`B_j`$, which sends $`P_t`$ to
an operator anticommuting with $`P_t`$. Exact conjugation through
the surrounding commuting bands and legitimate clean restriction give

```math
(Tt_C+t_CT)\Pi^{(j)}=0.                              \qquad(7)
```

All literal pair scalars cancel against their actual adjoints in these
conjugations; they remain present in the source word. Both projectors
in (6)--(7) are nonzero by Section 4.

Equation (6) forces $`t_A,t_C`$ to have the same full computational
X mask u. Acting on any computational basis vector in that nonempty
affine sector, each operator has one nonzero output basis vector.
Different masks cannot yield the same output vector; multiplying by
a computational-Z stabilizer changes a Z part, not an X mask.

For the complete Boolean target, the literal derivative satisfies

```math
T=t_A\Delta_u,\qquad
\Delta_u(x)=(-1)^{g(x)+g(x+u)},\qquad
\Delta_u(x+u)=\Delta_u(x).
```

The derivative therefore commutes with every Pauli of X mask u,
including $`t_C`$. If sigma is the Pauli commutation sign of
$`t_A,t_C`$, then

```math
Tt_C=\sigma t_CT
```

on the entire logical space. Equation (6) forces $`\sigma=1`$,
since both products restrict to identity on that sector. Equation (7)
would contradict this: $`Tt_C`$ is unitary and
$`\Pi^{(j)}\ne0`$. Therefore

```math
\langle b_j,t\rangle=0
\quad\text{for every }t\in L_A\cap p_A^\perp
\text{ and every }j.                                  \qquad(8)
```

No pair-parity promise replaces the complete target: the contradiction
compares its globally defined Boolean derivative across two actual
nonempty sectors.

## 6. All pair differences coincide, so at most one survives

On $`L_A`$ define the nonzero linear functional
$`\mu(l)=\langle p_A,l\rangle`$ and, for each original pair,
$`\lambda_j(l)=\langle b_j,l\rangle`$. Each lambda is nonzero:
the native displacement formula gives

```math
r_j=(F_J-I)b_j
   =\sum_{a\in A_{\rm active}}\langle a,b_j\rangle a\ne0.
```

Equation (8) says $`\ker\mu\subseteq\ker\lambda_j`$.
Over $`\mathbb F_2`$, the two nonzero functionals must therefore
be equal. Every actual A-active axis has mu equal to one, so

```math
r_j=\sum_{a\in A_{\rm active}}a
\qquad\text{for every }j.                             \qquad(9)
```

Distinct original pair differences are independent: a relation among
them would be the corresponding relation among their disjoint actual
middle axes. Equation (9) therefore proves $`m\le1`$. Under the
standing positive-pair premise, $`m=1`$ and input all-active support
has been proved without any assumption on singleton signs.

Apply the same argument to the reversed repaired circuit, with actual
middle band $`\widetilde E_B=F E_-E_+`$. This band has the same
independent signed physical axes, with the explicitly repaired native
orientations. Its clean image at the reversed input is the corrected
output image $`S=\widetilde E_BJ\widetilde E_B^\dagger`$.
Each original middle axis commutes $`\widetilde E_B`$, so

```math
SB_jS^\dagger=\widetilde E_B(JB_jJ)\widetilde E_B^\dagger=s_jD_j,
\qquad SUS^\dagger=\pm U
```

with the same actual pairs and singleton signs. The reversed occupied
target is again Boolean up to the conjugate common scalar, and the
actual reversed bands are still independent native lists. Thus the
same proof gives output all-active support in the corrected physical
list. Undo the repair: $`C_{\rm rep}`$ fixes r and transports
the actual output-active list bijectively, preserving its sum.

## 7. Negative singletons vanish and the existing sharp reduction applies

Return to the corrected one-pair comparison (1). Every negative
singleton commutes $`L_A`$ and anticommutes $`P_A`$. Its
literal T conjugation contributes an odd S rotation to J. Thus the
binary displacement difference is

```math
N_S=N_J+\sum_{u\in B_-}u\langle u,\,\cdot\,\rangle.
```

There are no cross terms because singleton axes commute the active
span. For either independent endpoint active list, its diagonal
displacement expression is

```math
\langle v,N_Jv\rangle
 =\left\langle\sum_{a\in A_{\rm active}}a,v\right\rangle,
\qquad
\langle v,N_Sv\rangle
 =\left\langle\sum_{c\in C_{\rm active}}c,v\right\rangle.
```

Both actual sums are now r. Nondegeneracy gives
$`\sum_{u\in B_-}u=0`$, and original middle independence
forces the negative list to be empty. This is the retained
[all-active clean comparison](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md#2-orientation-repair-removes-every-negative-singleton),
now reached without a source-class all-active assumption.

Consequently $`S=J`$ literally. Their intrinsic displacement images
give $`L_A=L_C=L`$, and their complete nonzero Pauli supports give
$`p_A+p_C\in L`$. The common rank k is even and at least two:
$`\langle b,a\rangle=1`$ for each actual active axis, while
$`0=\langle b,r\rangle=k\pmod2`$ because the paired middle
axes commute. The [fixed full-space diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)
now applies with all its original native provenance intact, proving
the bounds in Section 1.

## 8. Scope, resources and the next distinction

This is an exact structural theorem about actual one-clean native
three-band Boolean circuits. It uses the full logical, dirty and
reference contract, retains the common phase and the unrestricted
unused clean block, and adds no clean qubit or approximation error.
All native gates and their actual inverses are retained. The literal
repair and fixed coordinate Cliffords have their ordinary
$`O((Q+1)^2)`$ elementary-Clifford charges on the same register if
emitted. No scalar-sector restriction is an emitted operation, and
the proof is not a new source compiler or a free native-band deletion.

Booleanity is essential. The
[clean-return family with growing physical pair count](CLEAN_CODE_CLIFFORD_RIGIDITY.md#5-arbitrarily-large-rank-still-need-not-supply-a-shared-clean-flip)
has a complex diagonal T-phase occupied target, not a Boolean sign.
Its existence does not contradict the scalar-sector reality step.

The retained [original-circuit rounding theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
applies at initialized operator error strictly below $`2^{-530690}`$.
It first gives exact clean return and exact occupied Boolean action
for that same native circuit, up to its retained common phase when
phase optimization is permitted. If its actual original comparison
contains a pair, the collapse and degree bounds above then apply.
No pair hypothesis is imposed on an approximately defined clean image
before this exact rounding step. This conservative tiny-error
consequence does not establish the same result at arbitrary fixed
accuracy or replace the separate no-pair branch.

The paired branch is now exhausted at native depth three; the
no-original-pair branch keeps its separate established bounds. This
does not turn a limitation of a two-band source, a retained native
representation, or a particular extraction mechanism into a
dimension-growing lower bound for unrestricted T-depth. The complete
approximate dirty-indicator source and the general resource frontier
remain separate questions.

## 9. Complete native-three corollary, with the no-pair constant retained

The remaining case split is exhaustive. For an exact complete Boolean
sign circuit of T-depth at most three, on the full logical and dirty
register with at most one initialized clean qubit, retain the bound

```math
D_{\rm np}=D_3+576.
```

Here $`D_3`$ is the constant of the
[partial-diagonal feature theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants),
and the later [one-parity reduction](CENTRAL_PARITY_QUADRATIC_REDUCTION.md)
replaces the old rank-two finite-insertion application by 2049.
The high-rank value displayed here exceeds both 2049 and $`2^{26}`$.
Neither generic feature theorem is changed: the new rank-two argument
avoids using the finite-insertion theorem in this case. The original
pair-collapse proof and its six/five conclusion are unchanged.

| Actual circuit case | Established sign-degree bound |
| --- | --- |
| Zero clean qubits | $`5`$, by [full-space native-three rigidity](THREE_LAYER_NATIVE_STRUCTURE.md) |
| One clean qubit and either interior clean image Pauli | $`2^{26}`$, by the [Pauli-interior theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#7-one-pauli-interior-clean-image-suffices) |
| Both interior clean images non-Pauli, with no original middle pair | $`D_{\rm np}`$, by the [complete no-pair closure](GROWING_DISPLACEMENT_NATIVE_ORBITS.md#1-actual-source-hypotheses-and-conclusion) |
| At least one original middle pair | $`6`$, by Sections 1--7 above |

Therefore every such exact at-most-one-clean native-three circuit obeys

```math
\deg g\le D_{\rm np},\qquad
\deg f\le D_{\rm np}-1
\quad\text{for a complete Boolean XOR query}.                    \qquad(10)
```

For the XOR statement, conjugate the actual output by Hadamards to
obtain the sign $`(-1)^{y f(x)}`$. Those are Clifford endpoint changes
and preserve native depth, the clean/dirty/reference contract and
operator error. A nonzero predicate has sign degree one greater than
its own; the zero predicate is immediate. Empty native bands cover
smaller depth without changing any hypothesis.

The [R74-based original-circuit rounding argument](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
now gives the corresponding unrestricted-at-native-depth-three tiny-error
consequence: initialized operator error strictly below
$`2^{-530690}`$ forces the **same original circuit** to have exact
clean return and exact occupied Boolean action, up to the retained
common phase when allowed. Apply the exhaustive case split only after
that rounding. Thus targets above (10) cannot be approximated below
that threshold by any at-most-one-clean depth-three circuit, with
arbitrary dirty width. No original-pair premise is needed for this
combined corollary.

This establishes an absolute width-independent ceiling at native depth
three. The universal ceiling is the retained crude no-pair constant,
**not six/five**. A small universal constant remains open. Higher
native depths, arbitrary fixed-accuracy bounds, and dimension-growing
unrestricted-depth lower bounds are not covered; the full resource
frontier is unchanged.
