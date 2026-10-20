# A complete indicator from fermionic bank motion

[Current revision](BANK_MOTION_CONTINUATION_RESULT.md) ·
[Polynomial-overhead source interface](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Earlier phase-signal interfaces](PHASE_SIGNAL_QUERY_LIMITS.md)

An actual bank-moving word can implement the complete indicator without
initializing any bank state. A fixed fermionic Walsh transform turns
address-controlled translation into Clifford phases; conjugating one
bank Z by that translation gives exactly the selected bank Z. Fermionic
permutation signs cancel on every occupation sector, not only on a
single-excitation input.

The resulting native construction has polynomial address overhead in
count and no additional work, but its displayed T-depth is linear in the
address length. It therefore does not improve the retained logarithmic
source depth or R96. The first unpaid operation is a sufficiently shallow
synthesis of the fixed transform or of the whole routing echo. The
Gaussian representation limits below do not constrain unrestricted
Clifford interlayers.

## 1. The literal full-bank identity

Let $`b\ge1`$, $`m=2^b`$, and index the arbitrary m-qubit bank Y
by $`u\in\mathbb F_2^b`$. Identify its computational basis with
the occupation basis of m fermionic modes in this fixed order. For a
single-particle unitary V define

```math
\Gamma(V)=\bigoplus_{r=0}^m\bigwedge^r V.
```

In particular Gamma fixes the vacuum and is a literal homomorphism:
$`\Gamma(VW)=\Gamma(V)\Gamma(W)`$. Its action on the fully
occupied sector is the determinant, which is retained throughout.
Writing $`n_u=(I-Z_u)/2`$, a mode permutation P obeys
$`\Gamma(P)^\dagger n_u\Gamma(P)=n_{P^{-1}u}`$, including
all occupation sectors.

Let $`W=H^{\otimes b}`$ on the m-dimensional **mode-index space**,
and put $`F=\Gamma(W)`$. This notation does not identify F with b
physical Hadamards. For each arbitrary address x define

```math
\begin{aligned}
\Lambda_x&=\mathrm{diag}_{u}((-1)^{x\cdot u}),
&P_x|u\rangle&=|u\mathbin\oplus x\rangle,\\
C&=\prod_{i=1}^b\prod_{u:u_i=1}\mathrm{CZ}_{x_i,Y_u},
&C_x&=\Gamma(\Lambda_x).
\end{aligned}                                                   \qquad(1)
```

The ordinary Walsh identity gives $`P_x=W\Lambda_xW`$. Hence the
complete controlled word $`R=F^\dagger C F`$ has branches
$`R_x=\Gamma(P_x)`$, and

```math
\begin{aligned}
R_x^\dagger Z_0R_x&=Z_x,\\
D_b&=R^\dagger Z_0R
   =\sum_x|x\rangle\langle x|\otimes Z_x,\\
I_b&=H_YD_bH_Y,\qquad
I_b|x,Y\rangle=|x,Y\mathbin\oplus e_x\rangle .
\end{aligned}                                                   \qquad(2)
```

Here $`H_Y`$ applies H to every bank wire. Translation can produce
occupation-dependent fermionic signs, but the actual inverse in (2)
cancels them. Equation (2) is a full-operator identity, so arbitrary
address and bank superpositions, correlations and external references
are included. No helper or initialized single-excitation state is used.

## 2. Native butterfly stages, including their phases

For modes $`i\lt j`$, let $`H_{ij}`$ act as H on their
single-particle coordinates and as identity on the other coordinates.
Put

```math
S_{ij}=\prod_{i\lt r\lt j}Z_r,\qquad
P=X_iS_{ij}Y_j,\qquad Q=Y_iS_{ij}X_j.
```

The exact occupation-space identity is

```math
\Gamma(H_{ij})
=e^{-i\pi P/8}e^{i\pi Q/8}Z_j
=T(P)T(Q)^\dagger Z_j,
\qquad T(P)=e^{i\pi/8}e^{-i\pi P/8}.
                                                               \qquad(3)
```

The two T scalars cancel literally. After fermionic reordering, on the
two-mode factor's vacuum, one-particle and double-occupation sectors,
(3) has respective actions 1, H and minus one; the last sign is the
determinant of H. For a matching of mode
pairs, the resulting m axes are disjoint Majorana bilinears. They
commute and are independent, so they admit one native T layer.

An explicit elementary Clifford construction also pays the layer.
First bring all matched mode pairs next to one another by a fermionic
mode permutation. An adjacent fermionic swap is SWAP followed by CZ,
costing six elementary Cliffords: three CNOTs for SWAP and H,CNOT,H
for CZ. At most $`m(m-1)/2`$ adjacent swaps suffice for any desired
mode order. The complete routing and its actual inverse therefore cost
at most $`6m(m-1)`$ Cliffords.

For an adjacent pair, now with $`S_{ij}=I`$, use

```math
A=S_j\,\mathrm{CX}_{i,j}\,H_i\,\mathrm{CX}_{i,j},
\qquad
AZ_iA^\dagger=P,\qquad AZ_jA^\dagger=Q.
```

Thus (3) is $`A(T_iT_j^\dagger)A^\dagger Z_j`$.
A costs four Cliffords, its actual inverse costs six after expanding
$`S^\dagger=S^3`$, and Z costs two. The m/2 pairs run their
T gates on distinct physical wires. A complete matching stage has
common forward/actual-inverse bounds

```math
T_{\rm stage}=m,\qquad D_T(\mathrm{stage})=1,\qquad
G_{\rm stage}\le6m(m-1)+6m=6m^2,\qquad w=0.
                                                               \qquad(4)
```

In the inverse, reduce the whole adjacent Z inverse power
$`S^{-2}=S^2`$; no reduction crosses a noncommuting gate. This
keeps the same bound while using the actual inverse operator. All
permutation routing is also actually reversed.

The b index-coordinate Hadamards partition into b such matching stages.
Their exterior powers multiply to F, including their determinant
phases. Consequently

```math
T_F=mb,\qquad D_T(F)\le b,\qquad
G_F\le6bm^2,\qquad w_F=0,\qquad\epsilon_F=0.
                                                               \qquad(5)
```

Logical commutativity of the index-coordinate stages does not merge
their overlapping rotation axes into a native independent-axis band.
Equation (5) pays the displayed circuit and its actual adjoint; it is
not a lower bound on unrestricted synthesis of F.

## 3. Complete resource, error and return recurrence

Suppose a native implementation $`\widetilde F`$, possibly using
arbitrary private work, satisfies

```math
\|\widetilde F-F\otimes I_w\|\le\eta.
```

Let $`T_F,G_F,d_F,w_F`$ be common forward/actual-inverse bounds.
Use its actual adjoint in
$`\widetilde R=\widetilde F^\dagger C\widetilde F`$ and in
$`\widetilde I_b=H_Y\widetilde R^\dagger Z_0\widetilde R H_Y`$.
There are exactly four transform occurrences. Both C words together
contain bm controlled-Z gates. Therefore

```math
\begin{aligned}
T_E&\le4T_F,&G_E&\le4G_F+3bm+2m+2,\\
d_E&\le4d_F,&w_E&\le w_F,\\
\|\widetilde I_b-I_b\otimes I_w\|&\le4\eta,
&a_{\rm clean}&=0.
\end{aligned}                                                   \qquad(6)
```

Unitary telescoping proves the error bound without assuming intermediate
work return. It includes every return error, all bank values and all
reference correlations. The wrapper is an exact involution as an emitted
unitary even when the transform is approximate; this fact alone does
not give exact private-work return. With the exact transform (5), the
native indicator instead has

```math
T_E\le4bm,\qquad d_E\le4b,\qquad
G_E\le24bm^2+3bm+2m+2,\qquad w_E=0,
\qquad\epsilon_E=0.                                    \qquad(7)
```

All counts fit the relaxed polynomial-overhead source interface, but
the displayed depth does not. A transform satisfying (6) at every
fixed requested eta with polynomial count/width overhead and
$`d_F=o(\log(b+2))`$ would supply the missing complete source.
The paid root fusion would use eight source/actual-inverse occurrences,
at total error at most $`32\eta`$, and its retained exact suffix.
Its complete resource recurrence remains in the linked fusion proof.
Substituting (5) does not improve the frontier.

## 4. What the Gaussian representation still costs

For a fixed transform, restrict Clifford boundaries to preserve the
linear span of the standard $`2m`$ Majorana operators. Such a
Clifford signed-permutes those operators: a nontrivial linear combination
of distinct Majoranas cannot itself be a Pauli. Restrict each native
T layer to commuting quadratic Majorana axes. Distinct such axes commute
only when their Majorana pairs are disjoint, so one layer at most
doubles the support of a propagated Majorana. F sends each Majorana
to a linear combination of m distinct Majoranas. Thus this representation
requires $`d_F\ge b`$, matching (5).

There is also a robust obstruction for a directly synthesized mode
transform inside the echo. Replace W in (1)--(2) by any
number-preserving single-particle unitary V. Write

```math
a=V|0\rangle,\qquad s=|\mathrm{supp}(a)|,\qquad
r_x=V^\dagger\Lambda_xV|0\rangle.
```

Here $`n_r`$ denotes occupation of the normalized mode r.
Its phase echo is $`I-2n_{r_x}`$. On the two-mode span of
$`r_x`$ and $`|x\rangle`$, the number-operator difference has
eigenvalues plus and minus their angle's sine on the one-particle
sector, and zero on vacuum and double occupation. Hence its complete
Fock-space error is exactly

```math
\|(I-2n_{r_x})-Z_x\|
=2\sqrt{1-|\langle x|r_x\rangle|^2}.                 \qquad(8)
```

All r_x lie in the same s-dimensional subspace
$`L=V^\dagger\mathrm{span}(\mathrm{supp}(a))`$.
If the complete indicator error is at most $`\tau\lt2`$,
then (8) gives
$`\|\Pi_L|x\rangle\|^2\ge1-\tau^2/4`$ for every x.
Summing over x proves

```math
s\ge m(1-\tau^2/4).                                  \qquad(9)
```

If V is implemented by d disjoint two-mode mixing rounds, interspersed
with arbitrary mode phases and permutations, then $`s\le2^d`$.
This direct echo representation consequently needs

```math
d\ge b+\log_2(1-\tau^2/4),
```

and integer $`d\ge b`$ when $`\tau\lt\sqrt2`$.
This does not assume that V separately approximates W. It also does
not turn mode-mixing rounds into an unrestricted native-depth bound.
The stronger [whole-word Majorana theorem](WHOLE_WORD_MAJORANA_LIMITS.md)
now also covers fusion of the whole echo when **every actual Clifford
boundary** preserves the read-only address and the standard Majorana
space on all bank and arbitrary dirty-work wires. Its full-operator
argument does not require separate accuracy certificates for F or its
inverse. Arbitrary non-Gaussian Clifford boundaries and address motion
remain outside that theorem. A changed mechanism must still pay the
complete identity, native operations, simultaneous resources and actual
return.

The [focused native fixture](../../tests/test_bank_motion_walsh_indicator.py)
checks the two-address-bit construction on all 64 columns of its six
physical wires, including every bank occupation sector and the actual
inverse. Its emitted forward word has $`T=32,D_T\le8,G=258`$
after local inverse-S reductions, with zero clean and zero additional
dirty work. This supports the literal identity and gate ledger, not
the universal or asymptotic claims above.
