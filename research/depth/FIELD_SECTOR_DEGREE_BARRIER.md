# A degree barrier for Clifford sectors of augmented field decoders

The [moving-helper recoupler](AUGMENTED_FIELD_RECOUPLER.md) makes every recursive stage a Clifford in a
fixed helper-X sector. This gives a useful compilation interface, but does
not allow all stages to be collected into a bounded number of batches whose
sector symplectic matrices depend affinely on the helper characters. Their
required degree grows with the address length.

The result concerns the complete query-intertwining equation and a specified
helper-sector model. It is not a lower bound on unrestricted native T depth,
on the weaker complete-probe contract, or on a compiler whose intermediate
operations leave that model.

## 1. Complete contract and conclusions

Let $`G=\mathbb F_2^p`$, $`p\ge2`$, and $`m=2^p`$. The principal core consists
of an m-wire bank and q further arbitrary qubits. Let $`Q_x`$ be Clifford
queries on this whole core whose binary symplectic matrices have the form

```math
F_x=I+N_x,\qquad N_{x+y}=N_x+N_y,\qquad N_xN_y=0
\quad(x,y\in G).
```

The field shear, with any additional identity or Pauli-dressed core helpers,
has this square-zero label action. Literal query phases do not enter the
argument. In particular $`F_x^{-1}=F_x`$.

Let $`W_x`$ translate the bank wires by x. The target core query is
$`T_x=W_x\otimes R_x`$, where each $`R_x`$ is any fixed Clifford on the q
remaining core qubits, independent of helper-sector coordinates. Its
symplectic matrix is $`V_x`$. Write $`\Pi_x e_z=e_{z+x}`$ for the binary
wire-index permutation, so its bank block is $`\Pi_x\oplus\Pi_x`$. No restriction to identity decoded core helpers
is needed.

There is also an h-qubit sector-helper register. In its common X basis,
indexed by $`r\in\mathbb F_2^h`$, the helper query acts as

```math
P_x|r_X\rangle=\gamma_x(r)|(r+\Lambda x)_X\rangle,
\qquad |\gamma_x(r)|=1,
```

for one fixed binary linear map $`\Lambda:G\to\mathbb F_2^h`$. The retained
augmented construction uses a Pauli Z query, so $`\gamma_x(r)=1`$; the proof
also allows the scalar phases of a general Pauli query.

Assume the candidate decoder fixes every helper-X sector and has
Clifford core blocks:

```math
E=\sum_r |r_X\rangle\langle r_X|\otimes E_r,
```

where each $`E_r`$ is a Clifford up to an arbitrary scalar phase. Let $`M(r)`$
be its binary symplectic matrix. Suppose the complete operator error obeys

```math
\eta=\max_x\left\|
 E(Q_x\otimes P_x)E^\dagger-T_x\otimes P_x
\right\|_{\rm op}\lt \frac1{\sqrt2}.
```

All inputs, including the sector helpers and reference systems, are arbitrary.
Then:

1. $`\mathrm{rank}\Lambda\ge p-1`$, so at least p-1 independent helper
   shift coordinates are necessary in this model.
2. The maximum Boolean algebraic degree among the entries of $`M(r)`$ is at
   least p-1.
3. If $`\mathrm{rank}\Lambda=p-1`$, the restriction to **every** active
   helper coset has degree exactly p-1. Arbitrary dependence on the invariant
   helper coordinates is allowed.
4. A factorization into k stages with helper-affine core symplectic matrices,
   with arbitrary constant core Clifford walls, requires $`k\ge p-1`$.

Here degree means the unique square-free algebraic normal form over
$`\mathbb F_2`$. It concerns symplectic labels, not the scalar or Pauli-sign
functions of a Clifford tableau. The radius is sufficient and is not claimed
to be optimal.

## 2. Constant-error rounding gives an exact symplectic cocycle

At helper input r, both sides of the comparison shift the helper to
$`r+\Lambda x`$ with the same scalar $`\gamma_x(r)`$. Thus the full operator
norm is the maximum norm over these blocks, and

```math
\|E_{r+\Lambda x}Q_xE_r^\dagger-T_x\|_{\rm op}\le\eta.
```

Both unitaries in this comparison are Cliffords up to scalar phase. For any
Hermitian Pauli A, their conjugated Pauli images differ by at most
$`2\eta\lt \sqrt2`$. Distinct signed Hermitian Paulis have operator distance
at least $`\sqrt2`$. Their Pauli actions therefore agree exactly. The scalar
phases of the sector blocks need not agree, and are not discarded from the
operator-error premise; they simply cancel in this conjugation argument.

Taking binary labels gives the exact identity

```math
M(r+\Lambda x)F_x=V_xM(r),\qquad
M(r+\Lambda x)=V_xM(r)F_x.                         \qquad(1)
```

No approximate invariance or continuity of the discrete Clifford labels is
assumed. The discrete Pauli gap supplies it.

## 3. The helper shift needs rank at least p-1

If $`\ker\Lambda`$ contained independent b,c, equation (1) in any fixed
sector would imply

```math
M(r)N_b=(V_b-I)M(r),\qquad
M(r)N_c=(V_c-I)M(r).
```

Since $`N_bN_c=0`$ and $`M(r)`$ is invertible, it would follow that
$`(V_b-I)(V_c-I)=0`$. Its bank block is nonzero: on a coordinate vector
it produces the sum of the four distinct coordinates indexed by
$`z,z+b,z+c,z+b+c`$. Thus $`\dim\ker\Lambda\le1`$, proving the rank bound.
Any decoded residual Clifford block is irrelevant to this nonzero bank block.

## 4. A one-dimensional kernel forces a top-degree coefficient

First suppose $`\mathrm{rank}\Lambda=p-1`$. Let b span its kernel,
and choose a complementary subspace $`H\subset G`$ of dimension p-1. For
arbitrary base sector r, use coordinates on H to parametrize the active
helper coset $`r+\Lambda H`$. Its top algebraic-normal-form coefficient is

```math
C_r=\sum_{t\in H}M(r+\Lambda t).
```

All sums and products in this section are binary. The square-zero source
action gives $`F_tN_b=N_b`$, and equation (1) for b gives
$`M(r)N_b=(V_b-I)M(r)`$. Therefore

```math
C_rN_b
 =\sum_{t\in H}V_tM(r)F_tN_b
 =\left[\sum_{t\in H}V_t(V_b-I)\right]M(r).        \qquad(2)
```

The bank block inside brackets is

```math
\sum_{t\in H}\Pi_t(\Pi_b-I)=\sum_{x\in G}\Pi_x=J_m,
```

on both X and Z labels, where $`J_m`$ is the nonzero all-ones matrix. Hence
the bracketed matrix is nonzero, even with arbitrary decoded residual core
helpers. Right multiplication by the invertible $`M(r)`$ cannot make it zero.
Thus $`C_rN_b\ne0`$, and in particular $`C_r\ne0`$.

There are p-1 active variables, so this proves degree exactly p-1 on every
active coset. It is a nonzero matrix coefficient, not a conclusion drawn
from a few entries or a sampled set of helper characters.

For the identity-residual target, (2) is the particularly simple witness

```math
C_rN_b=(J_m\oplus J_m\oplus0)M(r).                \qquad(3)
```

## 5. Giving every address bit a helper shift does not remove the floor

The only other possible rank is $`\mathrm{rank}\Lambda=p`$. Restrict
to a full p-dimensional active helper coset. Rearranging (1) gives

```math
V_x=M(r+\Lambda x)F_xM(r)^{-1}.                   \qquad(4)
```

An entry of the regular bank-translation matrix is the indicator of one
address value, whose Boolean degree is exactly p. Consequently
$`\deg V_x=p`$. The source matrix is affine in x and the rightmost matrix
in (4) is constant. If the restricted degree of M were d, the right side
would have degree at most d+1. Thus $`d\ge p-1`$.

This case supplies a lower bound, not an assertion that every such decoder
has degree exactly p-1. Extra invariant helper coordinates and arbitrary
commutant choices cannot evade either case of the argument.

## 6. Consequence for jointly compiled affine stages

Suppose a proposed compilation factors its sector symplectic matrices as

```math
M(r)=C_k A_k(r)C_{k-1}\cdots A_1(r)C_0,
\qquad \deg A_j\le1,
```

with constant invertible symplectic matrices $`C_j`$. Every entry of the
product has degree at most k. The preceding theorem therefore forces
$`k\ge p-1`$.

Fixed affine reparametrizations of helper characters preserve Boolean degree.
They can also be inserted between stages: pull the fixed helper permutations
through the factorization and substitute the corresponding affine coordinates
in the following stages. Constant invertible core Clifford walls preserve
the maximum matrix degree as well. Arbitrary joint core/helper Clifford walls
need not preserve the chosen helper-sector algebra and are outside this model.

For example, one batch of phases that are helper-linear and quadratic in
core bits, in one fixed Clifford frame, has a helper-affine sector symplectic
matrix. It cannot implement this complete decoder from p=3 onward. A bounded
number of such stages cannot implement the family at growing p.

This is a stage-count statement for this compilation interface. A native
T band need not preserve the helper-X algebra or be a Clifford in every
sector. Conversely, a stage with a high-degree helper predicate need not
have an affine symplectic matrix. Neither possibility is excluded here.
The theorem gives no unrestricted native T-depth lower bound.

An explicit escape from an overbroad stage count is a reversible helper
Toffoli $`(r_0,r_1,r_2)\mapsto(r_0,r_1,r_2+r_0r_1)`$, followed by one core
Clifford shear controlled by the current third helper bit, followed by the
same helper Toffoli. The helpers return, but the completed sector matrix is
$`I+(r_2+r_0r_1)N`$, of degree two, despite only one affine core stage.
The nonlinear helper motion carries the additional degree and must be paid.
This does not construct a decoder for the present query family.

More generally, if the current helper coordinates at core stage i have
Boolean degree at most $`\Delta_i`$ and that stage has matrix degree at most
$`a_i`$ in those coordinates, then the final matrix degree is at most
$`\sum_i a_i\Delta_i`$. The coordinate degree is bounded by the product of
the degrees of preceding helper permutations. Thus bounded-degree nonlinear
helper stages can change the stage-count law, but cannot be counted as free
fixed Clifford walls. This bookkeeping remains within a classical helper
permutation and Clifford-core-stage model.

## 7. The augmented construction reaches the degree floor

In one recursion cell, let its helper character be $`\rho\in F`$. The
complete cell identity holds for every lower field parameter, not just the
chosen address subspace. At upper address bit zero it gives

```math
R_\rho=A_\rho R_0Q_\rho,
\qquad A_\rho=Q_{p-1}(\rho)\otimes Q_{p-1}(\rho). \qquad(5)
```

Use trace-F coordinates for both X and Z labels. The two field shears have
symplectic matrices $`I+\rho A`$ and $`I+\rho B`$ for constant F-linear
matrices A,B. The zero-character recoupler has an F-linear symplectic matrix
$`M_0`$: its coordinate changes and fixed trace-Fourier transforms are all
F-linear on these Pauli coordinates. Thus

```math
M_\rho
 =M_0+\rho(AM_0+M_0B)+\rho^2 AM_0B.               \qquad(6)
```

Both $`\rho\mapsto\rho`$ and $`\rho\mapsto\rho^2`$ are binary-linear.
Therefore this complete cell's binary symplectic matrix is helper-affine.
Parallel cells sharing the same level helper give a block-diagonal affine
matrix. There are p-1 levels, so the whole augmented decoder has degree at
most p-1. Its helper-shift rank is p-1; the lower bound proves degree exactly
p-1, uniformly.

Equation (6) concerns Pauli labels. Literal Pauli signs and overall sector
phases remain owned by the full augmented operator identity and its actual
inverse proof. The positive native ledger still pays sixteen T bands per
level; the degree theorem does not prove that this constant is optimal.

## 8. Verification and next operation

The [degree fixtures](../../tests/test_field_sector_degree.py) check the square-zero source against regular
translation derivatives, the robust signed-Pauli gap and scalar-phase caveat,
the active-coset coefficient identity, and the helper-affine local stages.
The [symbolic tableau fixtures](../../tests/test_field_sector_tableaux.py) separately retain signed conjugations,
actual inverses and cocycles at bounded sizes.

A successful improvement must change the compilation interface: for example,
use higher-degree sector stages with a genuinely cheaper complete native
implementation, move the helper algebra during synthesis, or allow non-Clifford
intermediate sectors. Each remains an open operation, with the entire native
word, actual adjoint and arbitrary-helper return still to be paid. The allowed
fusion overhead remains $`T,w=2^p\mathrm{poly}(p)`$ and
$`G=4^p\mathrm{poly}(p)`$, rather than arbitrary higher powers of $`2^p`$.
No improved general frame bound follows from this representation barrier.
