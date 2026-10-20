# One central parity removes the rank-two insertion

[Actual quadratic-control word](QUADRATIC_CONTROL_NATIVE_MATCHING.md) ·
[Pauli-boundary theorem](PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md) ·
[Complete no-pair closure](GROWING_DISPLACEMENT_NATIVE_ORBITS.md) ·
[Original-pair collapse](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md)

The rank-two quadratic insertion does not require its generic finite
matrix core. Fixing one genuine computational center parity makes that
insertion an actual Clifford. Each sector is then covered by the retained
Pauli-boundary three-band theorem, giving sign degree at most **2049**
after reconstructing the original Boolean function. This replaces the
dominant finite-insertion term in the previous universal bound. The
remaining universal ceiling is the retained high-rank constant
$`D_3+576`$; the asymptotic resource frontier is unchanged.

## 1. Full-space hypotheses and the paid replacement

Retain the literal no-original-pair normal form

```math
V=K_C^\dagger E_C^\dagger H^{b(z)}E_BE_AK_A
  =e^{i\theta}D_g.                                      \qquad(1)
```

It acts on the entire logical register, including arbitrary dirty inputs
and references. Each E is an actual independent commuting signed Pauli
T/T-dagger list. The signed Hermitian Pauli H commutes with every middle
axis and a genuine computational Z center R. Every band and both
normalized endpoint Cliffords fix R pointwise. The real quadratic
normalization and central cubic sign are retained separately, exactly as
in the [source proof](QUADRATIC_CONTROL_NATIVE_MATCHING.md#1-actual-input-controls-and-the-full-space-word).

Suppose the polar rank of b is two. Choose one fixed invertible linear
change of the actual computational center coordinates so that

```math
b(z)=z_1z_2+\ell(z),\qquad \deg\ell\le1.            \qquad(2)
```

Extend this to one fixed computational coordinate system on the entire
logical register. Such a change preserves Boolean degree and conjugates
all actual factors together. It does not depend on a selected sector.

For a bit s define the full-space operators

```math
\begin{gathered}
\Pi_s=\frac{I+(-1)^sZ_1}{2},\qquad
F_s=H^{s z_2+\ell(z)},\\
V^{(s)}=K_C^\dagger E_C^\dagger F_sE_BE_AK_A.
\end{gathered}                                        \qquad(3)
```

The exponent of F_s is affine in the actual center bits, including any
z_1 term in ell. Therefore F_s is an actual Clifford fixing R pointwise.
For a signed center parity P selecting a linear control, its literal
controlled-Pauli formula is

```math
\frac{I+P}{2}+\frac{I-P}{2}H.                       \qquad(4)
```

When P and H have independent commuting binary labels, a fixed Clifford
coordinate change identifies (4) with CZ. Coincident labels reduce it
to a Pauli or identity; an identity control or target gives the actual
Pauli or scalar factor. Negative representatives retain their signs in
(4). Products for affine controls remain Clifford because the controls
and H commute and $`H^2=I`$. No input-dependent scalar is discarded.
The retained affine-control construction uses $`O(Q^2)`$ elementary
Cliffords on Q logical wires, no T gates and no initialized helper.

On the selected full sector, the exact equality is

```math
F_s\Pi_s=H^{b(z)}\Pi_s,\qquad
V^{(s)}\Pi_s=V\Pi_s.                               \qquad(5)
```

All surrounding factors commute with the sector projector. Equation (5)
holds on every remaining logical column and after tensoring any reference
identity. It follows from the displayed full-space replacement; it is
not an assumed extension of an occupied identity.

## 2. Each sector has an actual three-band Pauli-code witness

Absorb the Clifford F_s by the exact definitions

```math
E_{C,s}=F_s^\dagger E_CF_s,\qquad
K_{C,s}=F_s^\dagger K_C.
```

Then

```math
V^{(s)}=K_{C,s}^\dagger E_{C,s}^\dagger E_BE_AK_A.   \qquad(6)
```

The three lists in (6) are actual independent native lists, with their
original sizes and literal orientations. Their axes and both endpoints
fix Z_1. Regard z_1 as the one initialized coordinate of this proof
witness. For s=1, conjugate the entire witness by X_1 to identify its
occupied sector with clean zero. Every native-band boundary has a Pauli
clean image: the whole bands preserve Z_1, and the Clifford encoders
used to expose their native rotations map Pauli codes to Pauli codes.

By (5), the occupied action is exactly

```math
e^{i\theta}(-1)^{g_s(x)},\qquad g_s(x)=g(s,x),      \qquad(7)
```

on the complete remaining register x. The common phase theta is the
same inherited phase in both sectors. The
[Pauli-boundary theorem](PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md#1-statement-and-physical-hypothesis)
therefore gives

```math
\deg g_s\le2048\qquad(s=0,1).                       \qquad(8)
```

One may also apply its full-logical-space proof directly. Quotienting
an independent commuting list by the one stabilizer direction creates
at most one binary relation. A projected identity retains its literal
scalar, and repeated axes retain the Clifford power and scalar produced
by combining them. The restricted endpoint operators are full logical
Cliffords, because their signed Pauli actions descend to the symplectic
quotient. The retained argument after its equation (2) uses only these
three relation bounds and complete Boolean diagonality. The
[bounded-relation theorem](BOUNDED_RELATION_NATIVE_MATCHING.md#1-statement-and-deliberately-loose-constants)
also explicitly retains the stronger 2048 bound in this case.

This quotient is a proof operation. It does not make a dependent
projected phase list a free native physical band. The full-space witness
(6) is an independent way to discharge that distinction.

## 3. Reconstruct the same Boolean function

The two sectors in (7) use the same remaining computational coordinates.
The Boolean identity

```math
g(z_1,x)=g_0(x)+z_1\bigl(g_0(x)+g_1(x)\bigr)        \qquad(9)
```

therefore proves

```math
\deg g\le2049.                                     \qquad(10)
```

The sums in (9) are over the binary field. No sector-dependent Clifford
change of the target variables is used in this reconstruction. Returning
through the fixed computational coordinate change preserves degree.
Restoring the retained central cubic and real quadratic signs adds
Boolean polynomials of degree at most three and two, so it does not
increase (10).

For a complete Boolean XOR query, apply this statement to its actual
output-Hadamard sign representation and verify the rank-two normal-form
hypotheses there. The sign is $`(-1)^{y f(x)}`$, so the corresponding
nonzero predicate has degree at most 2048. This is a conditional bound
for that branch, not a universal XOR ceiling of 2048.

## 4. A general polar-rank and relation-budget component

At polar rank $`2r`$, choose fixed actual center coordinates with

```math
b(z)=\ell(z)+\sum_{j=1}^r z_{2j-1}z_{2j}.          \qquad(11)
```

Fix the r odd-indexed coordinates to a string s. Replace the insertion
on the full register by

```math
F_s=H^{\ell(z)+\sum_{j=1}^r s_jz_{2j}}.            \qquad(12)
```

It is an actual Clifford and agrees literally with the old insertion
on that entire sector. All factors preserve the r selected stabilizers.
Their restricted endpoint Cliffords are actual logical Cliffords. The
kernel of each list's quotient map lies in the span of those r
stabilizers, so each restricted list has relation dimension at most r.
The full-space bounded-relation theorem consequently gives

```math
\deg g_s\le D(r,r,r),
```

where D is its explicit retained function. The fixed-coordinate selector
for s is $`\prod_{j=1}^r(1+z_{2j-1}+s_j)`$. Summing its product with
g_s over all s reconstructs g, proving

```math
\deg g\le D(r,r,r)+r.                              \qquad(13)
```

For r=1 use the stronger 2048 sector theorem, giving (10). For r=0
absorb the whole insertion globally and retain the sharper zero-clean
three-band sign bound five. These are alternative degree proofs to the
[finite-insertion bound](QUADRATIC_CONTROL_NATIVE_MATCHING.md#2-fixed-polar-rank-gives-a-uniform-degree-bound);
neither claims an unchanged-cost physical synthesis from a dependent
phase expression. The r-bit restriction cannot supply a uniform bound
when r grows, and no bounded-rank premise has been inferred from the
actual source.

## 5. The new universal ceiling and its precise scope

The retained
[complete no-pair case split](GROWING_DISPLACEMENT_NATIVE_ORBITS.md#1-actual-source-hypotheses-and-conclusion)
had the ceiling
$`\max\{D_3+576,D_{\mathrm{ins}}(0,0,0,3)\}`$.
Equation (10) replaces its rank-two finite-insertion term. Its rank-zero,
constant-flag and shared-clean-flip cases already obey the retained
smaller bounds. Thus the improved no-pair ceiling is

```math
D_{\mathrm{np,new}}=D_3+576.                        \qquad(14)
```

Here D_3 is the retained
[partial-diagonal constant](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants).
Its definition gives $`D_3\ge480\,2^{43}`$, so (14) exceeds both
2049 and $`2^{26}`$. Together with the
[original-pair collapse](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md#9-complete-native-three-corollary-with-the-no-pair-constant-retained),
the exhaustive actual-circuit cases are now:

| Actual circuit case | Sign-degree ceiling |
| --- | --- |
| Zero clean qubits | 5 |
| One clean qubit and either interior clean image Pauli | $`2^{26}`$ |
| Both interior clean images non-Pauli, no original middle pair | $`D_3+576`$ |
| At least one original middle pair | 6 |

Every exact complete Boolean circuit of native T-depth at most three
with at most one initialized clean qubit therefore satisfies

```math
\deg g\le D_3+576,\qquad
\deg f\le D_3+575
\quad\text{for a complete Boolean XOR query}.       \qquad(15)
```

Empty native bands cover smaller depth. The statement permits arbitrary
dirty width, intermediate Clifford mixing and unrestricted unused-clean
action, with complete clean/dirty return on the occupied code. Actual
adjoints preserve every identity: reverse the displayed words and retain
the conjugated scalars. No inverse is independently selected.

For approximation, apply the
[original-circuit rounding theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
first. Initialized operator error strictly below $`2^{-530690}`$
forces exact clean return and occupied Boolean action for that same
original circuit, up to the retained common phase when permitted. Only
then use the exact normal form, sector argument and case split. This
gives the corresponding higher-degree obstruction at that conservative
error threshold, including arbitrary references. It does not give an
obstruction at every fixed accuracy or at growing native depth.

The resource frontier remains R96. This reduction neither supplies the
complete sublogarithmic dirty-indicator source nor improves its count,
depth, error or return recurrence. The quantitative improvement removes
one expensive proof representation; the high-rank branch still controls
the universal maximum.

## 6. The next restriction must retain its native relation

The high-rank continuation can leave two endpoint translation masks
selected by one central parity. Fixing that parity may make the mask
constant, but can also create a relation in an actual endpoint list.
The independent-axis diagonal-peeling argument cannot then be applied
without a replacement proof.

Here is an exact one-band guard. On a quantum wire t and four center
wires, take

```math
E=\prod_{i=1}^4T_{X_tZ_i},\qquad H=Z_t,
\qquad P=Z_1Z_2Z_3Z_4.                             \qquad(16)
```

The four full-space axes commute and are independent. On the positive
P sector choose fixed computational coordinates
$`z_4=z_1+z_2+z_3`$. Their four projected axes have exactly one relation,
the sum of all four labels. Every one still has X_t support. The literal
sector operator is

| Remaining center string | E on the entire t fiber |
| --- | --- |
| 000 | $`X_t`$ |
| 111 | $`-X_t`$ |
| Each other string | $`iI`$ |

It follows that

```math
E_+^\dagger Z_tE_+=(-1)^{q(z)}Z_t,
\qquad
q(z)=1+z_1+z_2+z_3+z_1z_2+z_1z_3+z_2z_3.           \qquad(17)
```

This Hermitian Clifford has one constant computational X mask, namely
zero. Nevertheless none of its four activated projected axes is
computationally diagonal. E_+ itself is not Clifford: (17) is not a
Pauli because q has a nonzero quadratic part. The scalar i in the table
is retained, as are the signs in its two exceptional fibers.

This is a native one-band restriction counterexample to the proposed
dependent-axis inference. It is not a completed actual no-pair Boolean
source, and therefore does not rule out a stronger theorem using the
full source's matching and Boolean identities. The next structural
question is whether those additional identities eliminate the exceptional
relation after the parity restriction, or permit a direct reduction
that keeps and controls it. Merely obtaining a constant translation mask
does not yet pay that operation. The
[exact regression fixtures](../../tests/test_central_parity_quadratic_reduction.py)
check the literal relation guard; the degree conclusions follow from
the analytic arguments above.
