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

<a id="fermionic-spectral-error-contract"></a>

## 5. Complete accuracy of a fermionic lift is a spectral phase condition

The [approximate scalar zero reflection](EXPANDER_ZERO_REFLECTION.md)
suggests a different one-particle input to Gamma. Its constant-depth
full-operator approximation on the **index register** does not by itself
certify the induced operation on all occupations. The exact distinction
is as follows.

Let V and R be any unitaries on an M-dimensional mode-index space, and
let the eigenvalues of $`E=R^\dagger V`$ be
$`e^{i\theta_1},\ldots,e^{i\theta_M}`$, with principal phases
$`\theta_j\in(-\pi,\pi]`$. Multiplicativity of Gamma and the
exterior-power eigenvalues give the exact identity

```math
\|\Gamma(V)-\Gamma(R)\|_{\rm op}
 =\max_{S\subseteq[M]}
       \left|\exp\!\left(i\sum_{j\in S}\theta_j\right)-1\right|.
                                                               \qquad(10)
```

This includes every occupation sector and every reference extension.
The empty subset is the vacuum, whose phase is fixed. A common phase
of the **one-particle** V becomes a particle-number-dependent phase
under Gamma and cannot be discarded.

For $`0\lt\tau\lt1`$, put $`a=2\arcsin(\tau/2)`$. Then (10)
is at most tau if and only if

```math
\sum_{\theta_j\gt0}\theta_j\le a,
\qquad
\sum_{\theta_j\lt0}|\theta_j|\le a.                       \qquad(11)
```

Indeed, singleton subsets first force $`|\theta_j|\le a`$.
If the positive total exceeds a, the first positive partial sum
exceeding a lies in $`(a,2a]`$. Since $`a\lt\pi/3`$, that sum
lies in the forbidden arc $`(a,2\pi-a)`$, contradicting (10).
The negative phases have the same argument. Conversely, (11) places
every subset sum in $`[-a,a]`$, which proves sufficiency.

Thus the complete error controls total positive and negative phase
mass, rather than just the largest individual phase. A convenient
sufficient certificate is

```math
\|\Gamma(V)-\Gamma(R)\|
 \le\sum_{j=1}^M|e^{i\theta_j}-1|
 =\|E-I\|_1
 \le M\|V-R\|.                                         \qquad(12)
```

The first inequality follows by telescoping each subset product. If
$`\mathrm{rank}(E-I)\le r`$, its final bound improves to
$`r\|V-R\|`$. This is a sufficient certificate, not a necessary
rank restriction on successful mechanisms.

The loss can occur with determinant one: for even M, take M/2 phases
$`\pi/M`$ and M/2 phases $`-\pi/M`$. Then the one-particle
error tends to zero, but occupying only the positive-phase modes gives
complete error at least $`\sqrt2`$. A rank-one reflection has a
special protection. If $`R=I-2|x\rangle\langle x|`$ and
$`V=I-2|\psi\rangle\langle\psi|`$, with
$`\alpha=\arccos|\langle x|\psi\rangle|`$, the only nontrivial
relative eigenvalues are $`e^{2i\alpha}`$ and
$`e^{-2i\alpha}`$. Equation (10) gives $`2\sin\alpha`$,
independent of M, as in (8). Merely being an approximate scalar
zero-test does not establish this relative-rank property.

<a id="expander-fermionic-lift"></a>

## 6. An actual expander lift, its constant auxiliary dimension and its unpaid steps

There is an explicit native lift of the scalar source with only a
constant factor more modes. It is useful to state the family before
assessing its accuracy and addressing cost. In this section p is the
external address length and $`m=2^p`$.

First define the logical scalar word of
[the zero-reflection proof](EXPANDER_ZERO_REFLECTION.md#3-fixed-point-cubing-with-the-same-four-clean-bits)
using its exact specified permutations, rather than their dirty-helper
native decompositions. SELECT is a signed computational-basis
permutation on the tested bits and four selectors. Its CNOT matchings,
selected Z, and selector masks are explicit. The four-selector query
$`C_4`$ is the specified classical permutation; the bounded-clean
promotion's zero-test C is likewise the specified four-literal
permutation. Its K is the specified CCZ phase. All these operators are
exact, with their actual inverses.

The only logical auxiliary index bits needed are the four selectors,
the old query output, the new query output, and the promotion toggle:
seven bits. For $`n=3\lceil p/3\rceil=p+r`$, keep also the
$`r\le2`$ padding bits and perform the complete padding-sector cover.
Consequently the logical source has

```math
k=7+r\le9,\qquad d=2^k\le512,\qquad M=dm
                                                               \qquad(13)
```

mode coordinates. The logical reflection word, denoted
$`V_{p,\eta}`$, obeys the scalar proof's complete estimate

```math
\|V_{p,\eta}-(R_0\otimes I_d)\|\le\eta,
\qquad R_0=I_m-2|0^p\rangle\langle0^p|.                \qquad(14)
```

Here every logical auxiliary input is arbitrary. The four initialized
selectors occur only in the source proof before promotion; (14) has no
initialized input. This construction does not delete helpers from an
emitted native scalar circuit. It defines the same logical operators
without introducing their implementation helpers, then implements their
Gamma images directly on the M physical occupation wires.

That implementation is literal:

* A signed index permutation lifts to a fixed fermionic mode permutation
  and single-mode Z gates. Adjacent fermionic swaps are the six-Clifford
  words from Section 2; at most $`M(M-1)/2`$ such swaps suffice.
* An index H lifts to the matching stage (3)--(4), with M T gates,
  one physical T band, and at most $`6M^2`$ Cliffords. Index S or T
  lifts to S or T on each mode whose corresponding index bit is one.
  Index CNOTs and parity exposures are fixed index permutations.
* The determinant-one selector rotations use the actual finite
  single-qubit words from the scalar proof. Lift every H, S, T and
  T-dagger occurrence, retaining its literal phase and using actual
  reversals. No continuously variable occupation phase is supplied for
  free.

For each fixed one-particle tolerance eta the logical word has only
$`O_\eta(1)`$ signed-permutation blocks and single-qubit factors:
SELECT is one such permutation block, regardless of p. Therefore its
native lift $`L_{p,\eta}=\Gamma(V_{p,\eta})`$ has simultaneous
common forward/actual-inverse bounds

```math
T_L=O_\eta(M),\qquad G_L=O_\eta(M^2),\qquad
D_T(L)=O_\eta(1),\qquad a_{\rm clean}=0.                \qquad(15)
```

There is no work beyond its M occupation wires. These estimates are
native circuit bounds at fixed **one-particle** accuracy. They are not
complete-bank accuracy bounds.

The ideal lift marks all d auxiliary modes in row zero:

```math
\Gamma(R_0\otimes I_d)=\prod_{a=1}^d Z_{Y_{0,a}}.
```

Let $`P_x|u,a\rangle=|u+x,a\rangle`$ on mode coordinates and
$`\mathcal P=\sum_x|x\rangle\langle x|\otimes\Gamma(P_x)`$.
Define the fixed Clifford

```math
C=\prod_{u\in\mathbb F_2^p}\prod_{a=2}^d
          \mathrm{CX}_{Y_{u,a},Y_{u,1}}.
```

The complete emitted candidate is

```math
\widetilde D_{p,\eta}
 =C\mathcal P^\dagger L_{p,\eta}\mathcal P C^\dagger.
                                                               \qquad(16)
```

With the ideal lift in place of L, (16) is exactly
$`\sum_x|x\rangle\langle x|\otimes Z_{Y_{x,1}}`$
tensored with identity on all other modes. Indeed, mode translation
moves the marked row to x, and C maps that row's d-mode parity to its
first wire. The $`(d-1)m`$ helper modes need no initialization;
the ideal return is exact. The two C words cost $`2m(d-1)`$ CNOTs.
The norm error in (16) equals
$`\|\Gamma(V_{p,\eta})-\Gamma(R_0\otimes I_d)\|`$,
so (10)--(11) are precisely its full bank/helper/reference contract.
Hadamards on the m specified bank wires turn this phase into the
complete bit indicator.

There are two independent unpaid improvements.

**Precision.** The scalar estimate (14) supplies only the last bound in
(12). Taking $`\eta=\tau/M`$ is sufficient for full tolerance tau.
The logical factor count is $`O(N\ell)=O(L^2)`$, where
$`L=1+\log(1/\eta)`$, $`N=O(L)`$, and the synthesized selector
rotation word has length $`\ell=O(L+\log L)`$. Lifting those
factors serially therefore gives depth $`O_\tau(p^2)`$. This is the
cost of the displayed precision allocation, not a lower bound on the
actual expander lift. A direct proof of (11), a bounded-relative-rank
construction, or another total-phase estimate could improve this
allocation in a changed family. Section 7 shows that merely keeping
the present source's branching length fixed cannot do so, even with
arbitrarily precise auxiliary rotations.
The resulting count and width overheads are polynomial in p and hence
fit the paid fusion interface; its depth certificate does not.

**Coherent centering.** Equation (16) must implement
$`\mathcal P`$, not separately choose a classical circuit for each x.
An actual implementation uses p successive matchings: at stage j,
$`x_j`$ controls the $`M/2`$ fermionic swaps pairing
$`(u,a)`$ with $`(u+e_j,a)`$. Route each matching to adjacent modes
by fixed fermionic swaps, apply the controlled matching, and actually
unroute. A controlled adjacent fermionic swap is a Fredkin followed by
CCZ on its control and two modes. The
[shared-control batch schedules](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
and the [CCZ batch schedule](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
give at most eight T bands and $`6M+2`$ T gates per stage, with
$`O(M^2)`$ Cliffords and no clean or extra dirty work. Both orientations
of $`\mathcal P`$ thus add at most $`16p`$ T bands and
$`12pM+4p`$ T gates. These are genuine shared-control native
schedules, but the displayed depth is linear in p.

Second quantization does not automatically lift a circuit that moves
the **external physical address** through its intermediate operations.
It is a homomorphism on the mode-index operator, or fiberwise on an
explicit address-preserving family. Centering the scalar expander at x
conjugates its index generators by $`P_x`$; it does not remove the
coherent implementation obligation. Conversely, moving address and
bank labels together enters the
[separate joint-relabeling interface](EXPANDER_ZERO_REFLECTION.md#9-the-uniform-address-obstruction-to-joint-relabeling).

Thus this is a complete specified partial family, not a new shallow
indicator. The route is paused pending both a changed source with a
complete spectral-phase certificate or another way to avoid accumulation, and a
paid centering or replacement interface of the required depth. The
[polynomial-overhead fusion theorem](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
absorbs the displayed polynomial count overhead, not either growing
depth cost. This does not restrict general non-Gaussian native
interfaces, and it does not restore conjunction as a dimension-growing
fixed-accuracy depth witness.

The [spectral fixtures](../../tests/test_fermionic_lift_spectrum.py)
check exact determinant-one phase accumulation, the tilted rank-two
relative spectrum, the native two-mode lift's determinant sign, and
full-sector auxiliary parity decoding. They are finite identity checks;
(10)--(16), the native family and the resource conclusions are analytic.

<a id="sparse-expander-lift-obstruction"></a>

## 7. Fixed branching cannot give a complete expander lift

For unpadded lengths p divisible by three, the displayed scalar word has an
additional exact property: it is minus identity when the tested data
is zero, on **every** logical auxiliary input. The selected generators
fix zero, so $`S_j`$ is exactly identity there even with approximate
selector rotations and their actual inverses. The scalar query and
bounded-clean promotion are consequently exact on that data sector;
the final phase echo is minus identity. This argument is made only for
the unpadded lengths. The approximate padding-sector cover does not
automatically preserve an exact zero-data eigenvector.

Here is a general sparse-column consequence. Let
$`f=m^{-1/2}\sum_u|u\rangle`$, let U be a unitary on
$`\mathbb C^m\otimes\mathbb C^d`$ satisfying
$`U(f\otimes|a\rangle)=-f\otimes|a\rangle`$ for every a,
and suppose each column of U has at most s nonzero entries. Put

```math
V=(I_m-2|f\rangle\langle f|)\otimes I_d,
\qquad A=U-I.
```

For the column $`j=(u,a)`$, the sum of A's entries in rows with
auxiliary label a is exactly minus two. This follows from the left
eigenvector relation, which follows in turn from unitarity and the
specified right eigenvector. The column has at most $`s+1`$ nonzero
entries, so Cauchy--Schwarz gives

```math
\|Ae_j\|^2\ge\frac4{s+1},\qquad
\|(U-V)e_j\|^2=\|Ae_j\|^2-\frac4m.
```

The second equality retains the exact cross term from that same
minus-two row sum. Summing over all dm columns and using
$`\|U-V\|\le2`$ gives

```math
\|U-V\|_1
 \ge\frac12\|U-V\|_{\mathrm F}^2
 \ge2d\left(\frac m{s+1}-1\right).                    \qquad(17)
```

If $`\|\Gamma(U)-\Gamma(V)\|\le\tau\lt1`$, (11) implies
$`\|U-V\|_1\le2a`$, where $`a=2\arcsin(\tau/2)`$.
Thus a necessary condition is

```math
m\le(s+1)(1+a/d).                                     \qquad(18)
```

Apply this after conjugating the logical scalar word by Fourier
Hadamards on its p **tested data** bits only. This is a mathematical
basis change, used only in the proof; no native Fourier transform is
being supplied. Every data CNOT becomes the reversed CNOT, and the
selected data Z becomes a selected X. Thus SELECT and all other
non-auxiliary blocks remain signed permutations in this basis. The
remaining mixing blocks act only on the d-dimensional auxiliary
register. With b such blocks, even allowing each to be an arbitrary
unitary of unlimited precision, every column has support at most

```math
s\le d^b.
```

The exact zero-data sector becomes precisely the prescribed
$`f\otimes\mathbb C^d`$ eigenspace. Equations (17)--(18) therefore
apply to the actual scalar word and its literal fermionic lift. They
require

```math
b\log_2d\ge p-O_\tau(1).                              \qquad(19)
```

Even allowing $`d\le(p+2)^c`$ for fixed c requires
$`b=\Omega_\tau(p/\log(p+2))`$ within this block grammar.
For the unpadded source $`d=128`$ is constant, while the number of
auxiliary mixing blocks is $`O(N)`$, with N the fixed-point source
length. It follows that complete fixed-tolerance accuracy requires
$`N=\Omega_\tau(p)`$ in this family. Increasing only the precision
of its selector rotations cannot repair a fixed N. The explicitly
serial matching compiler in Section 6 consequently has growing depth;
its constant-length version fails even before coherent centering.

This is a structural obstruction for the specified sparse index word,
not an unrestricted native T-depth theorem. The variable b counts
auxiliary mixing blocks, not arbitrary native bands. A new family
could violate the sparse-column premise or use a different physical
interface. The unbounded subsequence of unpadded lengths already
excludes this fixed-branching family as a general complete indicator.
The spectral fixture also checks (17)'s exact column cross term on a
small signed-permutation example. No finite fixture proves (18)--(19).
