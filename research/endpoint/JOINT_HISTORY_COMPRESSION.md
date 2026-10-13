# Joint history compression and fine-prefix extraction

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Packed histories](PACKED_REJECTION_BANKS.md) · [Research index](../README.md)

A joint compressor may mix rejection banks and their logical coordinates.
It need not implement the separate inverse prefixes of the bankwise
construction. Nevertheless, three calls to any such compressor and its
actual inverse recover a complete fine prefix. At the depthwise endpoint
schedule, adding a constant number of existing layers completes the target.
Thus a linear-count joint compressor would itself suffice for the endpoint,
without executing the raw cascade. No such compressor is constructed here.
Mixing banks while preserving their decoded logical labels has a separate
constant-error obstruction; that restricted obstruction is not an
unrestricted T-count lower bound.

## 1. The normalized history is an isometry

Use the ideal depthwise schedule in
[Section 12 of the packed-history proof](PACKED_REJECTION_BANKS.md#12-depthwise-balance-fits-the-literal-precision-width).
There are D retained depths, three stages per depth, and nested supports
$`P_d`$. Let $`J_C`$ insert a logical input into the cap, and let
$`J_{d,a}`$ identify $`P_d S`$ with its private stage bank, extended by
zero on $`P_d^\perp`$, so $`J_{d,a}^\dagger J_{d,a}=P_d`$.
Write $`W_{d,a}`$ for the complete chronological prefix before that stage,
including the actual separate coarse operations. It preserves $`P_d`$.
For each logical coordinate x, let e(x) be its first active depth, or D
if none. The cap and stage part of the uncompressed initialized output is

```math
\mathsf F
=J_C\,\mathrm{diag}(c_{e(x)})
 +\sum_{d=0}^{D-1}\sum_{a=0}^2
 c_d J_{d,a}W_{d,a}P_d,
\qquad c_d=2^{d-D}/8.
\qquad\text{(1)}
```

All maps are tensored with the identity on every dirty input. Orthogonal
banks and the scalar weight identity give

```math
\mathsf F^\dagger\mathsf F=I/64,
\qquad X_0=8\mathsf F,
\qquad X_0^\dagger X_0=I.
\qquad\text{(2)}
```

In particular, $`X_0`$ is defined on every logical and dirty column, not just
the frame's first column. A full unitary H is an epsilon-accurate joint
compressor here when

```math
\|H\mathsf F-J_C/8\|\le\epsilon.
\qquad\text{(3)}
```

H may act arbitrarily on the occupied complement, even on the success
and preweight banks. Only (3), unitarity, and use of its actual inverse
are required. Set $`X=H^\dagger J_C`$. Unitary invariance immediately gives
$`\|X-X_0\|\le8\epsilon`$. Reversing the compressor therefore prepares
the normalized history on the entire input space. This statement prices
neither H nor the abstract history in (1).

## 2. Three calls recover the last fine prefix

Choose the last stage $`(d,a)=(D-1,2)`$, put $`P=P_{D-1}`$, and abbreviate
its embedding as $`J_F`$ and its prefix as $`W_*`$. Extend $`W_*`$ by the
identity on $`P^\perp`$. Since $`c_{D-1}=1/16`$ and an input outside P
never enters a retained stage,

```math
J_F^\dagger X_0=\tfrac12 W_*P,
\qquad X_0(I-P)=J_C(I-P).
\qquad\text{(4)}
```

Let $`R_C=I-2J_CJ_C^\dagger`$ and
$`R_F=I-2J_FJ_F^\dagger`$. The complete three-call word is

```math
\mathcal Q_H=-H^\dagger R_C H R_F H^\dagger,
\qquad
\mathcal Q_HJ_C=-(I-2XX^\dagger)R_FX.
\qquad\text{(5)}
```

For the ideal column $`X_0`$, write $`A=J_F^\dagger X_0`$. Direct multiplication,
without any choice of an ideal occupied-input completion, yields

```math
-(I-2X_0X_0^\dagger)R_FX_0
=X_0+2J_FA-4X_0A^\dagger A
=J_FW_*P+J_C(I-P).
\qquad\text{(6)}
```

Here $`A^\dagger A=P/4`$. Thus both the active and inactive logical
columns have the required literal phase. Let L be the full swap of
$`J_C(P S)`$ with the stage bank, fixing every other mode. Then (6) becomes
$`J_CW_*`$ after L.

For isometries X and Y of equal finite domain dimension,

```math
\|XX^\dagger-YY^\dagger\|\le\|X-Y\|.
```

Indeed the projector difference norm is the larger of the two
off-subspace norms, and each is at most the column difference norm.
Apply this fact with $`Y=X_0`$ to (5). One changed column and one changed
reflection give the full-output estimate

```math
\|L\mathcal Q_HJ_C-J_CW_*\|
\le3\|X-X_0\|\le24\epsilon.
\qquad\text{(7)}
```

No small-epsilon assumption or bound on H's other columns is needed.
The estimate includes leakage, arbitrary dirty inputs and external
references. It does not replace the actual inverse with an ideal inverse.

### Native routing and the endpoint consequence

Prepare and remove the cap flag pattern with the same fixed flag X.
The resulting initialized-success circuit is
$`X_f L\mathcal Q_H X_f`$. The leading minus in (5) is the exact Clifford
scalar $`-I=(HS)^{12}`$ on a retained wire. The cap reflection tests only
the two flags. The stage reflection tests its actual flags, type bits,
sentinel and fixed zeros. The swap L is the
[complete cap-bank route](PACKED_REJECTION_BANKS.md#9-physical-allocation-and-exact-native-scalar-compression)
with X as its pair gate. These are full-input permutations and reflections;
all occupied and inactive sectors are specified.

The route and predicates cost $`O(n^2)`$ T gates and use one arbitrary
returned helper. This helper is borrowed sequentially from the same dirty
pool, including when that pool is entangled after H. There is no additional
initialized flag or clean suffix. Consequently the prefix circuit has

```math
T_{\rm prefix}\le3T(H)+O(n^2),
\qquad \epsilon_{\rm prefix}\le24\epsilon.
\qquad\text{(8)}
```

For $`D=n-8`$, append the final residual axis, the final exactly logical
coarse word, and the last eight target depth layers. This follows the
residual-then-coarse ordering of the packed construction. The nine fine
layers use the [disjoint SU(2)-block theorem](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks),
whose fixed-volume, fixed-error-budget guarantee includes the literal endpoint width;
their total row volume is O(N). Each has cost $`O(N+\ell)`$ at
$`\ell=N+O(1)`$, and the final coarse word costs O(N). Give each fine
layer error at most $`\eta/64`$ and take $`\epsilon\le\eta/64`$ in (3).
The complete target error is then at most $`33\eta/64\lt\eta`$, with

```math
T_{\rm target}\le3T(H)+O(N),
\qquad a=2,\qquad b=N+n+7.
\qquad\text{(9)}
```

The constant-precision-offset reservations of those nine existing layers
fit this literal width for sufficiently large n; the same finite small-n
fallback as in the packed proof applies. H itself must fit that reservation.
No raw cascade, rejection-orientation bank, or preweight circuit is called
in this reduction. Equation (9) is a conditional compiler theorem, not a
construction of a short H or a lower bound against one. It shows why a
joint compressor with this contract is already a sufficient fine-synthesis
primitive.

## 3. Preserving logical labels leaves a constant error floor

A narrower joint mixer can mix all copies of a decoded logical coordinate
x, including arbitrary dirty work, while preserving x. Let $`\Pi_x`$
project onto its cap and stage copies, tensored with the full dirty space.
The remaining modes can be left fixed, or assigned to orthogonal sectors
that complete this decomposition. Assume $`[H,\Pi_x]=0`$ for every x.
H may otherwise depend on the target and on x.

For an input label j and any unit dirty vector z, (1) gives

```math
\|\Pi_x\mathsf F|j,z\rangle\|^2
=c_{e(j)}^2\,\delta_{xj}
 +\sum_{d,a}c_d^2
   |(W_{d,a}P_d)_{xj}|^2.
\qquad\text{(10)}
```

Missing bank coordinates contribute zero. A unitary H commuting with
$`\Pi_x`$ preserves this probability exactly. The desired scalar cap
column has no support at x different from j, so

```math
\|H\mathsf F-J_C/8\|
\ge\left(\sum_{x\ne j}\sum_{d,a}
c_d^2|(W_{d,a}P_d)_{xj}|^2\right)^{1/2}.
\qquad\text{(11)}
```

The probabilities in different banks add before mixing. Unitary
interference within a fixed logical label cannot erase their total norm.

For a native-compatible constant-gap example, take $`D\ge2`$, all ideal
residual axes equal to identity, and all coarse coins identity except the
root coin $`HZ\in SU(2)`$, applied after its three residual stages.
This is an exact Clifford coin. The three last-depth prefixes all contain
it. For either root-pair input, the other root coordinate has squared
amplitude 1/2 in each of those three prefixes. Since the last rejection
weight is 1/16, (11) gives

```math
\|H\mathsf F-J_C/8\|\ge\sqrt{3/512}.
\qquad\text{(12)}
```

The zero residuals satisfy the near-unit source chart. Earlier banks can
only increase this lower bound. Arbitrarily large dirty space and
target-dependent mixers do not alter it. Logical-label mixing, including
the fine-prefix bank, is outside this restriction.

## 4. Evidence and remaining question

The [joint-compression tests](../../tests/test_joint_history_compression.py)
check normalized histories, the complete extraction word, its inactive
columns and literal phase, arbitrary occupied-input completions, entangled
dirty inputs, the perturbation estimate, and the conserved-label example.
They are bounded matrix diagnostics, not an emitted compressor or an
asymptotic gate-count certificate.

Scalar or logical-label-preserving bank mixing is stopped by (11).
General mixing remains possible, but by (9) a linear-count implementation
would itself complete the endpoint. Reopening that route requires an
actual native synthesis rule for H, with its inverse, literal workspace
and full column error priced. The existing endpoint bounds are unchanged.
