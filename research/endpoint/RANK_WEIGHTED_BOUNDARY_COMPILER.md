# Rank-weighted blocks and collective boundary fusion

[Disjoint-block primitive](DISJOINT_BLOCK_COMPILER.md) · [Signal and recovery contracts](CONJUGATE_SIGNAL_RECOVERY.md) · [Canonical forests](MATCHING_LEAF_DEFLATION.md) · [Endpoint](../../docs/OPEN_PROBLEM.md)

A support-times-rank budget replaces the squared-size hypothesis for disjoint
blocks. A rank-preserving coarse circuit makes the residual small enough for
a constant-normalization sum of reflections. Branches occupy mode coordinates,
not a newly initialized logarithmic register. The same theorem prices a whole
sequence of boundary corrections when its collective basis and core satisfy
the weighted budgets below. It does not assert those budgets for every Hopf
frame, or make independently compiled overlapping factors free to compose.

## 1. Contract and theorem

Let $`N=2^n`$, $`n\ge3`$, $`\eta=2^{-N}`$. Fix $`\kappa,c\gt0`$,
independent of n, target entries and ranks. Let U be a computable direct sum
of $`U_j\in SU(m_j)`$ on known disjoint subsets of N computational labels,
with identity elsewhere. Suppose supplied or proved integer bounds satisfy

```math
\mathrm{rank}(U_j-I)\le r_j\le m_j,\qquad
\sum_jm_jr_j\le\kappa N.
```

There is a terminating coherent Clifford+T compiler with

```math
T=O_{\kappa,c}(N),\qquad G=O_{\kappa,c}(N^2),\qquad
\|\widehat UJ_2-J_2(U\otimes I_b)\|\le\eta/c,
\qquad b\le N+n+7.
```

All logical/dirty inputs, reference systems, literal phases and clean leakage
are included. Only T-count is optimized. Certified finite classical search
is allowed; no exact-rank/eigenvector oracle, supplied source, measurement,
reset or uncharged coherent operation is used. For real final blocks the
[conjugate-signal replacement](CONJUGATE_SIGNAL_RECOVERY.md#one-clean-real-block-corollary)
reduces the initialized budget to one, with the same T and dirty bounds.
It does not remove the mode/recovery flag.

A special unitary with $`r_j\le1`$ is identity; discard these using the
promise, not an exact test on a computable matrix. The theorem includes
$`r_j=m_j`$ and permits, for example, $`m=N^{3/4},r=N^{1/4}`$.
The disjoint-support condition is not replaced by a sum of ranks of
arbitrary overlapping sequential operations.

## 2. Exactly logical coarse synthesis preserving the rank budget

The following calls use the [dirty interpreter and sector echo](../../docs/BORROWED_WORKSPACE_COMPILER.md),
[full-mode reflection banks](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources),
and [literal fine phases](../../docs/ONE_CLEAN_COMPILER.md#literal-diagonal-phases-use-the-same-flag).
All their preparation, controls, routing, inverses and return are included.

First handle supports $`m_j\gt N/128`$. Their total rank is below
$`128\kappa`$. Spectrally each target is a product of selective phases
$`I+(e^{i\theta_a}-1)v_av_a^\dagger`$. A complex state extension, center
phase and actual inverse implement one selective phase with at most thirteen
fine diagonal slots. Only a fixed number is required for these large supports.
Phase ambiguity of a vector cancels in its projector. Certified searches over
orthonormal-frame/phase meshes replace exact spectral selection, with a strict
final-product tolerance; no condition-number promise is imposed.

For the remaining supports put $`\rho=\max_jr_j`$ and
$`\delta_0=1/(512\rho)`$. Ideally

```math
U_j=B_j(D_j\oplus I)B_j^\dagger,\qquad D_j\in SU(r_j)
```

with diagonal D. The first $`r_j`$ columns of B have an extension using at
most $`r_j`$ complex Householder reflections, up to column phases commuting
with D. At each position replace the ideal reflection by the actual word

```math
\widehat P(I-2|e\rangle\langle e|)\widehat P^\dagger.
```

Here $`\widehat P`$ is an exactly logical native state preparation on the
specified block, with exact dirty return and disabled padding. Its actual
reflection still has rank one relative to identity. Magnitudes use block
Hopf preparation; complex phases use adjacent determinant-one Rz differences,
with a common vector phase immaterial to the projector. Previous pivots need
not be fixed exactly by the approximate basis: only unitarity and block
support are required for the next rank bound.

Approximate D by a determinant-one native word on the chosen $`r_j`$ pivot
coordinates. It need not remain diagonal. The resulting coarse block is

```math
U_{0,j}=B_{0,j}(D_{0,j}\oplus I)B_{0,j}^\dagger,
\qquad \mathrm{rank}(U_{0,j}-I)\le r_j,
\qquad \|U_{0,j}-U_j\|\lt\delta_0.
```

Physically the full word is exactly $`U_0\otimes I_{\rm dirty}`$.
Accuracy $`O(\delta_0/\rho)`$ per frame step is inverse-polynomial in N,
so single-qubit word length is a fixed polynomial in n.

### Compressed-region resource sum

For a dyadic rank class $`R/2\lt r_j\le R`$ its total support is
$`S_R\le2\kappa N/R`$, and $`R\le2\sqrt{\kappa N}`$. Group by support
size as well, padding and packing each class in compressed regions. There
are at most R reflection positions. At one position, disjoint magnitude
preparations, phases, center signs and pivot exchanges use tables on
$`O(S_R)`$ rows. The fixed-prefix echo does not expand them to N rows.
The cost per position is $`\mathrm{poly}(n)\sqrt{S_R}`$. One-time
pack/unpack costs over the polynomially many classes are
$`\mathrm{poly}(n)\sqrt N`$; later permutations stay inside the compressed
regions. Therefore

```math
T(U_0)\le\mathrm{poly}(n)\sum_R R\sqrt{S_R}
 +\mathrm{poly}(n)\sqrt N
 =O_\kappa(\mathrm{poly}(n)N^{3/4})=O_\kappa(N).
```

The polynomial degree is fixed, not rank-dependent. Clifford cost is
$`O_\kappa(\mathrm{poly}(n)N)=O_\kappa(N^2)`$, and dirty work is
$`O_\kappa(\sqrt N+\mathrm{poly}(n))`$. Paying a new arbitrary N-label
permutation at every position would not prove this estimate.

For effective selection, enumerate the bounded-length native frame and
SU-core words with their syntactically guaranteed support/determinant/rank
properties. Refine matrix comparisons until a candidate is strictly within
$`\delta_0`$. The approximation argument supplies one within half that
tolerance, so this search terminates even at degeneracies. Classical time
is not optimized.

## 3. Small residual, constant coefficient normalization

Define $`E_j=U_{0,j}^\dagger U_j`$ from the actual coarse word. Then

```math
\|E_j-I\|\lt\delta_0,\qquad
\mathrm{rank}(E_j-I)\le q_j:=\min(2r_j,m_j).
```

The factor two is necessary when the coarse and target active spaces differ.
No coarse error is added again: the exact factorization is $`U=U_0E`$.
Select exact-unitary spectral candidates approximating E with strict final
error margin and norm distance below $`2\delta_0`$ from identity. Rational
complex unitary-frame meshes and rational unit-circle phases are dense;
unused slots have phase zero. Only rational candidate zeros are tested.
This supplies finite computable parameters without exact spectral selection.

For a spectral piece write

```math
E=I+\sum_{a=1}^{q}z_av_av_a^\dagger,
\quad z_a=e^{i\theta_a}-1,\quad H_a=I-2v_av_a^\dagger.
```

Its exact linear-combination identity is

```math
E=a_0I+\sum_{a=1}^{q}a_aH_a,\qquad
 a_0=1+\tfrac12\sum_a z_a,\quad a_a=-z_a/2.
```

With the preceding rank and norm bounds,

```math
\lambda=|a_0|+\sum_a|a_a|
\le1+q\|E-I\|\lt5/4.
```

Add cancelling identity coefficients $`\pm(2-\lambda)/2`$. Their sum is
zero but the sum of absolute coefficients becomes exactly two. This is
standard LCU/oblivious amplification; the restricted-workspace implementation
below, not a free SELECT/PREPARE oracle, is required here. See Berry et al.,
[PRL 114, 090502](https://doi.org/10.1103/PhysRevLett.114.090502).

### Branch copies and bins

For an m-dimensional spectral piece, use $`q+3`$ copies of its coordinates.
Branch zero is the original input support. The other $`m(q+2)`$ coordinates
lie in the mode-flag-one sector. For every input x prepare

```math
P|x\rangle=\sum_{\ell=0}^{q+2}\sqrt{|a_\ell|/2}\,|\ell,x\rangle.
```

These supports are disjoint across x. The full preparation is a real
state-extension bank. SELECT applies $`H_a`$ in its own branch copy; those
complex reflection supports are disjoint. A literal diagonal D supplies
coefficient phases. Thus $`\mathcal U=P^\dagger D\,\mathrm{SELECT}\,P`$
has initialized block E/2. Inactive original labels have block I. Eigenvectors
appear as disjoint reflection normals, not as a free complete eigenbasis.

For $`m_j\le N/128`$, split the spectral terms into chunks of at most

```math
L_m=\lfloor N/(8m)\rfloor-3\ge13.
```

Each chunk has weight $`m(q+3)\le N/8`$. The number $`t_j`$ of chunks
obeys $`t_j\le32m_jq_j/(3N)+1`$, since $`L_m\ge3N/(32m)`$. Then

```math
\sum_jm_j(q_j+3t_j)
\le\tfrac54\sum_jm_jq_j+3\sum_jm_j\le4\kappa N,
\qquad \max_jt_j\le64\kappa/3+1.
```

Process equal-index chunks in rounds; each round retains disjoint original
supports. Greedy bins of weight at most N/4 have every nonlast weight above
N/8. At most
$`H_\kappa=32\kappa+\lceil64\kappa/3\rceil+2`$ bins suffice.
Within a bin, extra modes are below N/4; padded preparation and SELECT
supports are each at most N/2. They can be independently packed and
unpacked inside 2N modes, with full action fixing borrowed padding.

## 4. Native recovery and uniform accounting

One overlap uses at most fifteen fine diagonal slots: four for preparation,
four for its inverse, six for complex SELECT, and one for coefficient phase.
The [full recovery identity](CONJUGATE_SIGNAL_RECOVERY.md#3-recovery-must-be-checked-on-both-sectors)
for $`-\mathcal UZ_f\mathcal U^\dagger Z_f\mathcal U`$ gives active output
E with zero ideal rejection. Inactive inputs acquire a known minus sign,
removed by the exact Boolean phase primitive. A bin uses at most 45 fine
slots, not a number proportional to its rank.

Including the large-support cases, a conservative slot bound is
$`F_\kappa=10000\lceil\kappa+1\rceil`$. Choose

```math
q=N+\lceil\log_2(240cF_\kappa)\rceil.
```

The native hybrid is below $`\eta/(8c)`$. Assign less than
$`\eta/(4c)`$ separately to spectral-product selection and certified
source/phase evaluation; retain strict final slack. Actual inverse words
satisfy the inverse-target isometry bound. A hybrid compares ideal initialized
prefixes and propagates prior actual leakage by unitarity, without reset.

Only the mode flag and fine signal are initialized. Fine address width is
at most $`n+1`$; the shared dirty core has $`q+1`$ qubits and self-borrowing
requires $`q+1\ge2(n+1)`$. For fixed promises this and
$`q+1\le N+n+7`$ hold for sufficiently large n. Coarse/routing work reuses
the pool. The existing direct literal SU-synthesis fallback handles the
finite smaller dimensions. This enlarges only a uniform constant. No
conditioning-dependent quantum count or unspecified coherent call remains.

For real final blocks, apply the conjugate-sector replacement to this whole
word before phase compilation. It uses one initialized mode flag and
$`q+2`$ arbitrary dirty qubits, with the real-target isometry error unchanged
apart from the allocated fine hybrid. Arbitrary complex targets retain the
two-clean theorem. Neither statement resolves a sequential overlapping word.

## 5. Boundary-column reduction and collective closure

Let $`e_i`$ be anchors with disjoint interiors $`Y_i`$. Every lower stage
V preserves $`\mathrm{span}\{e_i\}\oplus Y_i`$ separately; U acts only
on anchor groups. Choose a full extension S with $`Se_i=Ve_i`$ in each
child. Then $`S^TV`$ fixes the anchors and commutes with U, proving

```math
\mathscr K_U(V):=UVU^TV^T=USU^TS^T.
```

For real unit $`v_i=Ve_i`$, choose a sign $`\sigma_i`$ with
$`1+\sigma_i(e_i^Tv_i)\gt1/2`$. The normalized vector
$`(e_i+\sigma_iv_i)/\sqrt{2+2\sigma_i e_i^Tv_i}`$, its Householder
reflection, and a center sign supply S, including $`v_i=\pm e_i`$.
Disjoint state-extension banks therefore price this correction without
implementing the lower block interiors. The alternative representative-pair
basis yields a core of dimension at most $`2r_g`$ in upper group g, with
squared-size budget $`4\sum_g r_g^2`$. Boundary columns, not arbitrary
state-preparation replacements for V, justify the cancellation.

Transported and untransported products are distinct:

```math
\mathscr K_U(AB)=\mathscr K_U(A)A\mathscr K_U(B)A^T.
```

Thus a correctly transported sweep through a fixed child decomposition
reduces to the final boundary columns. An ordinary product of corrections
needs the collective spaces below.

For a finite sequence $`K_t=U_tV_tU_t^TV_t^T`$, suppose the same anchor
groups are preserved at every step. Let $`d_i\le m_i:=\dim Y_i`$ bound

```math
\mathcal F_i=\mathrm{span}_t\{V_te_i-(e_i^TV_te_i)e_i\}\subseteq Y_i.
```

Choose one common basis B fixing anchors and mapping $`d_i`$ representative
coordinates to a padded frame for $`\mathcal F_i`$. It can be chosen with
$`B_i\in SO(m_i)`$ and $`\mathrm{rank}(B_i-I)\le\min(2d_i,m_i)`$,
using at most $`d_i`$ pairs of reflections; frame signs are free. The full
ordered product has the exact form

```math
K_R\cdots K_1=B\left(\bigoplus_g C_g\oplus I\right)B^T,
\qquad D_g=\dim C_g=r_g+\sum_{i\in g}d_i.
```

Indeed $`V_tU_{t,g}^TV_t^T=I+F_{t,g}(A_{t,g}^T-I)F_{t,g}^T`$,
where F contains the boundary columns. Each correction acts only on the
anchors and their collective transverse spaces. All intermediate B bases
cancel literally. Compute the whole $`C_g`$ product classically before its
single native compilation.

If known bounds $`q_g\le D_g`$ satisfy $`\mathrm{rank}(C_g-I)\le q_g`$
and

```math
\sum_i m_i\min(2d_i,m_i)\le\kappa_BN,\qquad
\sum_gD_gq_g\le\kappa_CN,
```

the rank-weighted theorem prices B, the complete core, and its actual inverse
in $`O_{\kappa_B,\kappa_C}(N)`$ T gates, independent of R and
$`\max_i d_i`$. Real targets need one clean flag. Certified SO-frame/core
search with strict product margin avoids exact rank selection. Allocating
$`\eta/4`$ to selection, $`\eta/8`$ to each B appearance and
$`\eta/4`$ to the core gives total at most $`3\eta/4`$.
The fixed-channel theorem is the corollary with bounded $`d_i`$ and
$`q_g=D_g`$, not a competing stronger claim.

## 6. Scope tests for repetition

For one shared anchor, the commutator is the product $`H_uH_vH_uH_v`$
of two reflection normals $`e-Ue`$, $`e-Ve`$. Its rank relative to identity
is at most two, and its norm is
$`d_Ud_V\sqrt{1-d_U^2d_V^2/16}`$, where $`d_U=\|Ue-e\|`$ and similarly
for V. For an r-dimensional computational intersection, expanding the
commutator through its intersection projector gives rank at most $`2r`$.
This rank can grow under ordinary multiplication; two noncommuting
rank-two corrections on a five-label path can have rank-four product.

In particular, for $`U=A\oplus I`$ and a matching lower rotation with
nonzero sine s and cosine c, the correction is

```math
K=\begin{pmatrix}
c^2I+s^2A&cs(I-A)\\cs(A^T-I)&c^2I+s^2A^T
\end{pmatrix},\qquad
\mathrm{rank}(K-I)=2\mathrm{rank}(A-I).
```

A full common core can have quadratic support-times-rank weight. This is
not an intrinsic cost lower bound: A can itself be a cheap signed permutation.
Likewise a path's successive boundary columns can span a growing space.
For equal positive rotations their transverse matrix is
$`\mathrm{diag}(s,s^2,\ldots,s^h)T`$, where T has diagonal one and every
entry above it equal to c. Its inverse has upper entries
$`-c(1-c)^{j-i-1}`$, giving minimum singular value at least $`s^h/2`$.
Intermediate directions cannot simply be declared irrelevant at endpoint
precision; a transported sweep may instead telescope them away.

For approximate boundary columns, an isometry replacement obeys
$`\|K_t-\widetilde K_t\|\le2\|U_t-I\|\|F_t-\widetilde F_t\|`$.
Small interactions can reduce the required boundary accuracy, but this
inequality supplies neither a new source nor a summed linear T ledger.

[Focused tests](../../tests/test_boundary_refinement.py) check residual-rank
doubling, reflection-LCU normalization and complete recovery, bin bounds,
boundary cancellation, transported and common-basis products, and rank-growth
witnesses. Their finite matrix checks are not native compiler emissions or
proofs of the asymptotic bound; the charged construction above supplies that
bound. Householder synthesis is standard, for example Malvetti et al.,
[Quantum 5, 412](https://doi.org/10.22331/q-2021-03-15-412). The unrestricted
Hopf hierarchy still requires a global cost certificate for its actual
changing boundaries or another complete representation.
