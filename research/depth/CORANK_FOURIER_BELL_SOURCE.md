# Corank-one Fourier kernels and complete Bell acceptance

[Kernel translation query](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-xor-carry-source) ·
[Paired Gram contract](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) ·
[Four-call recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

The paid corank-one translation matrix has an exact singular Fourier
operator whose left Gram is the desired central-character projector.
Its complete two-bank unitary dilation is also paid. This does not
supply the required lookup source: selecting a Bell fiber uses an entire
auxiliary bank, and the paired shared-outcome construction has success
weight exponential in minus the carrier width. A concrete coherent-SWAP
repair instead has a fully normalized one-clean Gram and maximal lookup
error at every size. Sections 7--8 pay a normalized source for the thin
projector using the scalar zero reflection, and then a two-reflection
source that succeeds completely at n=2 but has an explicit missing
even-charge subspace at every n at least three. Section 9 pays a broad four-path swap/relative-T source, whose coherent
auxiliary sectors also have exact Gram gap one half. These are statements
about the displayed words, not general native-depth obstructions. Section 10
transfers the complete charge-readout error to auxiliary acceptance. The separate [Clifford projector-Gram theorem](CLIFFORD_PROJECTOR_GRAM.md)
covers fully interleaved Clifford fibers and gives a bounded-clean
affine-query count obstruction. The complete-frame frontier is
unchanged. The [current research question](../../WORKSPACE.md#main-unfinished-route-and-exact-missing-operation)
organizes continuation; the source-specific reopening conditions below retain
their stated hypotheses.

## 1. The singular Fourier operator and its unpaid dilation

Use the standard dual coordinates of the
[nilpotent-pencil source](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#3-the-full-quadratic-xor-query-is-paid).
There are p address bits, m=2^p and n=pm bits in each carrier. Let
q=2^n be the carrier Hilbert-space dimension; q is not m. Put

```math
 B_x=\bar N_x^{\mathsf T}-I,\qquad
 K_x=(I+B_x)P=B_x+P,
 \qquad P=I+e_0e_0^{\mathsf T}.
```

All matrices here are binary. P deletes only the first binary
coordinate. With the retained binary basis for the field powers,

```math
 \ker K_x=\langle e_0\rangle,\qquad
 \ker K_x^{\mathsf T}=\langle v_x\rangle,\qquad
 v_x=(1,x,\ldots,x^{m-1})^{\mathsf T},\qquad
 \mathrm{im}K_x=\ker(v_x^{\mathsf T}).                    \qquad(1)
```

The field entries of v_x are expanded as binary coordinates. In
particular v_x is nonzero even at x=0. Define the q-by-q matrix

```math
 (F_x)_{b,z}=q^{-1/2}(-1)^{b^{\mathsf T}K_xz}.
```

Summing the binary characters over b or z gives, respectively,

```math
 F_x^\dagger F_x=I+X(e_0),\qquad
 F_xF_x^\dagger=I+X(v_x).                              \qquad(2)
```

Thus F_x divided by sqrt(2) is a partial isometry. The abstract accepted
operator F_x-dagger divided by two has Gram
$`(I+X(v_x))/4`$. The rank defect is exactly one binary direction, but
(2) is not a native synthesis of its one-clean unitary dilation.

For the paired contract, write
$`Q_x=X_C(v_x)X_D(v_x)`$ and $`\Pi_{\rm eq}=(I+Q_x)/2`$.
The fixed envelope consisting of $`E^{-\mathsf T}`$ on each carrier,
then carrier Hadamards, takes the physical paired lookup charge to
this X-basis charge. Its forward and reverse use at most
$`4c_n+4n`$ Cliffords in total, with $`c_n=n^2+3n`$.
No address-controlled inverse is included in this envelope.

The usual Pauli corrections already locate the difficulty. For every
binary a,

```math
 F_xX(a)=Z(K_xa)F_x.
```

In the other direction, $`F_xZ(s)=X(c)F_x`$ requires
$`K_x^{\mathsf T}c=s`$. Solvability is the known condition s_0=0;
selecting a solution c still requires the inverse linear solve modulo
the unknown vector v_x. Neither the constant corank nor the first
covariance identity supplies that correction.

## 2. An actual unitary Bell kernel on arbitrary banks

Let A be the n-bit data bank and B an n-bit auxiliary bank, both
arbitrary. The paid XOR query and complete Clifford walls give

```math
 Q_{K_x}|z,y\rangle=|z,y+K_xz\rangle,\qquad
 W_x=\mathrm{SWAP}_{AB}
       (H_A^{\otimes n}\otimes H_B^{\otimes n})Q_{K_x}.
                                                               \qquad(3)
```

Operator products are rightmost first. If the output data and
auxiliary labels are b and a, direct application of the two Hadamard
walls gives the complete matrix element

```math
 \langle b,a|W_x|z,y\rangle
   =q^{-1}(-1)^{a\cdot z+b\cdot(y+K_xz)},\qquad
 \langle a|_B W_x|y\rangle_B
   =q^{-1/2}Z(y)F_xZ(a).                                \qquad(4)
```

Every a,y block is retained. This proves the action on arbitrary
coherent address, data, auxiliary, native-helper and reference inputs.
In particular (4) is not a promise that B starts or ends at zero.

For the reversed unitary define
$`T_{y,a}=\langle y|_B W_x^\dagger|a\rangle_B`$. Equations (2)
and (4) give

```math
 \begin{aligned}
 T_{y,a}&=q^{-1/2}Z(a)F_x^\dagger Z(y),\\
 T_{y,a}^\dagger T_{y,a}
   &=q^{-1}\bigl(I+(-1)^{v_x\cdot y}X(v_x)\bigr).
 \end{aligned}                                             \qquad(5)
```

The formula holds for every fixed a; it does not delete the off-diagonal
blocks between different a inputs. Initializing a and selecting one y
would yield $`2^{1-n}\Pi_{v_x,+}`$ for even v_x dot y, not the required
$`\Pi_{v_x,+}/2`$ except at n=2. It also uses n initialized/selected
auxiliary bits. Accepting all y gives identity, as unitarity requires.

For diagonal acceptance weights $`0\le w(y)\le1`$, attaining
$`\Pi_{v_x,+}/2`$ on a fixed-a input requires

```math
 \sum_y w(y)=q/4,\qquad
 \sum_y(-1)^{v_x\cdot y}w(y)=q/4.                       \qquad(6)
```

Nonnegativity forces w(y)=0 on every odd v_x dot y. Thus this
acceptance rule must identify the missing kernel parity. Selecting
even ell dot y for one fixed nonzero linear label ell instead gives
I/2 unless ell equals v_x. This last statement concerns the fixed-a diagonal
block, not the full arbitrary-auxiliary operator.

## 3. Paired shared outcomes preserve the charge but dilute the norm

Consider two copies of (5) with the same selected outcome y and any
fixed auxiliary inputs a_C,a_D. The accepted map for this outcome is
$`T_{y,a_C}\otimes T_{y,a_D}`$. Write
$`s_y=(-1)^{v_x\cdot y}`$. Its Gram is

```math
 q^{-2}(I+s_yX_C(v_x))(I+s_yX_D(v_x)).
```

Keeping distinct outcomes orthogonal and summing over their common y
removes the individual characters. Since v_x is nonzero, the complete
sum is exactly

```math
 \sum_y
 (T_{y,a_C}\otimes T_{y,a_D})^\dagger
 (T_{y,a_C}\otimes T_{y,a_D})
   =q^{-1}(I+Q_x)=2^{1-n}\Pi_{\rm eq}.                 \qquad(7)
```

Shared byproducts do preserve the paired sector. Their normalization
is the problem: equality of two n-bit outcomes retains only a
q-inverse fraction of their possible pairs. At p=1, n=2, (7) happens
to have the desired coefficient one half; at p=2, n=8, it is one
over 128. The general coefficient is exponential in minus pm, so
inverse-success amplification is outside the allowed polynomial
overhead in m. Standard square-root amplification, if separately
available, would already use $`\Omega(2^{n/2})`$ calls. Its required
state preparation and success reflections are not supplied here.

Independent acceptance of both outcomes sums to identity. Accepting
only $`v_x\cdot(y_C+y_D)=0`$ would instead select the correct paired
parity, but this is the missing lookup-character readout. Coherently
erasing distinct outcome labels is a different operation from (7)
and requires its own unitary and complete arbitrary-input analysis.
Equation (7) is a fixed-auxiliary diagonal-block obstruction; no
claim is made that a prepared Bell resource or an outcome projection
returns arbitrary dirty auxiliary input.

The full arbitrary-auxiliary operator is also explicit. Let
$`P_{=,YZ}`$ project onto equal computational labels of the two
output auxiliary banks, and let $`W_{2,x}=W_{x,CY}W_{x,DZ}`$.
For the reverse Bell word its complete accepted Gram is

```math
 \begin{aligned}
 A_x&=W_{2,x}P_{=,YZ}W_{2,x}^\dagger\\
 &=q^{-1}\sum_{\ell\in\mathbb F_2^n}
 X_C(\ell)X_D(\ell)
 X_Y(K_x^{\mathsf T}\ell)X_Z(K_x^{\mathsf T}\ell).
 \end{aligned}                                           \qquad(7a)
```

Indeed $`P_{=,YZ}=q^{-1}\sum_\ell Z_Y(\ell)Z_Z(\ell)`$ and
$`W_xZ_B(\ell)W_x^\dagger=X_A(\ell)X_B(K_x^{\mathsf T}\ell)`$.
The q summands are a faithful commuting translation group, so A_x
is a projector of rank $`q^3`$ on the $`q^4`$-dimensional four-bank
space. The summand ell=v_x is exactly Q_x. Thus
$`Q_xA_x=A_x`$: this projector imposes the desired paired charge
and n minus one extra independent stabilizer conditions. Its range
is smaller than the desired half-space by a factor $`2^{n-1}`$.
Taking a fixed computational auxiliary diagonal in (7a) retains
only ell=0 and ell=v_x, which recovers (7). The other summands move
the arbitrary auxiliaries coherently and cannot be discarded. No
one-clean implementation of the equality projector is assumed.
Section 7 supplies a native approximate source for $`A_x/2`$. Its ideal Gram distance from
$`\Pi_{\rm eq}/2`$ would be exactly one half: A_x is a proper
subprojector for every $`n\ge2`$. In particular the coincident p=1
fixed-auxiliary coefficient in (7) is not a complete-input exception.

## 4. A normalized one-clean coherent-SWAP repair fails completely

There is a paid coherent alternative to testing equality of all outcome
bits. Let $`W_{2,x}=W_{x,CY}W_{x,DZ}`$, where Y,Z are arbitrary
n-bit auxiliary banks, and conjugate their whole-bank swap:

```math
 S_x'=W_{2,x}\mathrm{SWAP}_{YZ}W_{2,x}^\dagger.
```

Multiplying the three binary shear matrices, or summing (4), gives
its complete action

```math
 S_x'|c,d;y,z\rangle
   =|d,c;\ y+K_x^{\mathsf T}(c+d),\
                 z+K_x^{\mathsf T}(c+d)\rangle.         \qquad(8)
```

The transpose matters. This is a Hermitian involution with no scalar
phase. The auxiliary shifts cannot be replaced by their trace or by
their action on an initialized plus state.

Use one clean signal a, and define the literal native rotation
$`R_a=(S_aH_a)T_a(S_aH_a)^\dagger`$. Write the controlled involution as

```math
 \mathcal C_a(S_x')=|0\rangle\langle0|_a\otimes I
       +|1\rangle\langle1|_a\otimes S_x'.
```

Then emit

```math
 U_x=R_a^3\mathcal C_a(S_x')R_a^\dagger,\qquad
 K_{{\rm acc},x}=\langle0|_aU_x|0\rangle_a
             ={1+i\over4}(I+S_x'),\qquad
 G_x=K_{{\rm acc},x}^\dagger K_{{\rm acc},x}
             ={I+S_x'\over4}.                          \qquad(9)
```

On the plus eigenspace of S_x' the signal word is R_a squared,
with accepted scalar (1+i)/2. On the minus eigenspace it has zero
accepted amplitude. Thus the one-half projector normalization is
exact, including the native global phase. This is a complete source
on arbitrary Y,Z, not a normalized trace.

For any auxiliary basis labels y,z, the normalized vector

```math
 |\psi_x\rangle=2^{-1/2}
   (|0,0\rangle-|v_x,v_x\rangle)_{CD}\otimes|y,z\rangle_{YZ}
```

has $`S_x'|\psi_x\rangle=|\psi_x\rangle`$ by (8), while
$`Q_x|\psi_x\rangle=-|\psi_x\rangle`$. Consequently

```math
 \left\|G_x-\Pi_{\rm eq}/2\otimes I_{YZ}\right\|=1/2
 \qquad\text{for every p and every x}.                  \qquad(10)
```

The upper bound follows because both operators are one half times
orthogonal projectors. The same witness also defeats multiplying
S_x' by a fixed data-bank swap: c=d on both components. The result
stops this concrete coherent-SWAP repair, not arbitrary coherent
interference or signal-interleaved acceptance.

## 5. The native schedule, actual inverses and precision

Put $`r=(m-1)p^2`$ and $`s\le(m-1)p^3`$. The returned bilinear
query has the paid common bounds

```math
 T_B\le8r,\qquad d_B\le8,\qquad h_B\le r,\qquad
 G_B\le2n+16r+2s.
```

Since $`Q_{K_x}=Q_{B_x}Q_P`$, adding n minus one fixed CNOTs,
2n Hadamards and n literal three-CNOT swaps proves

```math
 T_W\le8r,\qquad d_W\le8,\qquad h_W\le r,\qquad
 G_W\le G_B+6n-1.                                     \qquad(11)
```

The extra B bank is a principal arbitrary input of W, not one of its
returned r query helpers. Every query helper returns exactly. The
reverse is the actual reversed native word and has the same bounds.
No fine-angle gate or approximation occurs in (3).

In (9), implement the controlled S_x' using two W calls, a
signal-controlled swap of Y,Z, and two actual reversed W calls.
Compile each Fredkin as two outer CNOTs and an exact returned
eight-T Toffoli. A private phase mask gives the Toffoli as eight
T/T-dagger gates, eight parity CNOTs and two target Hadamards. The
n masks are distinct; shared signal a occurs only as a Clifford
control, so all n Toffolis use eight T layers together. Masks may
reuse the query pool between sequential calls. The two signal
rotations use two T gates and at most 15 Cliffords in the common
forward/actual-inverse bound, using T cubed equals ST and positive
S representatives of actual adjoints.

Including the physical-to-X-dual envelope from Section 1, the source
and its actual inverse therefore obey

```math
 \begin{aligned}
 T_U&\le32r+8n+2,&d_U&\le42,\\
 G_U&\le4G_W+12n+15+4c_n+4n,\\
 b_U&\le2n+\max(r,n),&a_{\rm clean}&=1.
 \end{aligned}                                           \qquad(12)
```

Here the two n-bit principal carriers and the p-bit address are
excluded from b_U; Y,Z and every private phase mask are included.
Physical T targets in each query are private returned masks. The
query calls, Fredkin word and signal rotations are sequential, so
no shared T target is scheduled twice in a layer. Address is
preserved, and every native helper returns on arbitrary inputs;
Y,Z have the explicitly displayed action (8).

The whole family is exact in the native alphabet. If each emitted W
call has full operator error at most eta, with actual inverses used
in reverse positions and other primitives exact, telescoping gives
$`\|\widetilde U-U\|\le4\eta`$ and Gram error at most 8 eta
relative to (9). This does not repair its target gap (10). All
errors and identities tensor with arbitrary references. The budgets
in (12) fit the relaxed polynomial source envelopes, but the wrong
projector prevents invoking the completed recovery/fusion theorems.

## 6. Evidence and the next unpaid operation

The proofs above are analytic. The
[exact fixtures](../../tests/test_corank_fourier_bell_source.py)
check the p=1 literal native word on every address, bank and helper
column; the p=1,2 field matrices and full singular Fourier Grams;
the p=2 shared-outcome coefficient using actual integer Walsh
propagation; the complete p=1 equality projector including off-diagonal
auxiliary blocks; the coherent-SWAP transpose; and the one-signal native
phase and actual inverse. These finite checks challenge conventions
and normalization. They are internal verification, not external
mathematical review or an all-circuit lower bound.

The singular Fourier partial isometry remains a precise target for
a different acceptance circuit. The later construction below supplies
a normalized native thin-projector source and a complete n=2 repair.
At growing size the missing operation is still a paid extraction of
the central character on its entire half-space, with arbitrary
auxiliary inputs. Shared outcomes, the coherent swap and the two
reflections below have different complete operators, but none supplies
that growing-size extraction. Reopening requires a changed acceptance
word whose Gram covers the missing even-charge subspace.

<a id="paid-thin-bell-source"></a>

## 7. The scalar zero reflection pays the thin half-projector source

There is now a complete implementation of the previously hypothetical
thin source. This does not require a clean Bell bank. Let
$`Q_{K_x^{\mathsf T}}(C;Y)`$ add $`K_x^{\mathsf T}c`$ into Y,
and define the returned shear word

```math
 \mathcal L_x=Q_{K_x^{\mathsf T}}(C;Z)
              Q_{K_x^{\mathsf T}}(C;Y)\,
              \mathrm{CX}_{C\to D}.
```

The n CNOTs in the last factor act coordinatewise. These three ideal
shears commute, preserve C and the address, and are involutions.
Conjugating X_C(ell) gives the summand in (7a), so the complete identity is

```math
 A_x=\mathcal L_x
       (|+^n\rangle\langle+^n|_C\otimes I_{DYZ})
       \mathcal L_x^\dagger.                            \qquad(13)
```

The plus vector in (13) identifies a projector, not an initialization.
Put $`R_{A,x}=I-2A_x`$ and let $`R_0^{(n+1)}`$ be the complete
zero reflection on the existing signal a and carrier C. The actual
controlled reflection has the ideal full-space word

```math
 \begin{aligned}
 \mathcal C_{A,x}
   &=\mathcal L_x H_C X_a R_0^{(n+1)}X_a H_C\mathcal L_x^\dagger\\
   &=|0\rangle\langle0|_a\otimes I
     +|1\rangle\langle1|_a\otimes R_{A,x}.
 \end{aligned}                                           \qquad(14)
```

Here H_C means all n carrier Hadamards. Testing zero after X_a tests
a=1 and C=0; consequently (14) uses a zero reflection on n+1 tested
bits, not a supplied controlled version of an approximate unitary.
The [complete scalar zero-reflection construction](EXPANDER_ZERO_REFLECTION.md#7-arbitrary-lengths-and-the-complete-zero-clean-ledger)
provides an emitted $`\widetilde R_0^{(n+1)}`$ with no initialized
work and full error at most eta. Denote its common native
forward/actual-inverse resources by $`T_0,G_0,d_0,h_0`$.
All its borrowed inputs and reference correlations are arbitrary.
Its return is approximate within eta, not assumed exact.

Keeping the literal R_a from (9), emit the one-clean source

```math
 U_{A,x}=R_a^3 Z_a\mathcal C_{A,x}R_a^\dagger,
 \qquad
 K_{A,x}={1+i\over2}A_x,
 \qquad K_{A,x}^\dagger K_{A,x}=A_x/2.                  \qquad(15)
```

On A_x the middle signal operation is identity and the accepted
amplitude of R_a squared is (1+i)/2; on its complement the middle
operation is Z_a and the accepted amplitude vanishes. Equation (15)
holds on the whole auxiliary space, including the moving auxiliary
blocks in (7a). The emitted source satisfies

```math
 \|\widetilde U_A-U_A\|\le\eta,\qquad
 \|\widetilde K_A^\dagger\widetilde K_A-A_x/2\|\le2\eta.
                                                               \qquad(16)
```

The actual reversed word uses the actual inverse of the emitted zero
reflection. An approximate zero reflection need not be exactly
Hermitian or involutory.

For the simultaneous ledger, $`PB_x^{\mathsf T}=B_x^{\mathsf T}`$,
so $`K_x^{\mathsf T}=P+B_x^{\mathsf T}`$. The transposed bilinear
query therefore has the same r-mask, eight-layer bound, plus n minus
one fixed CNOTs. With $`g_K=G_B+n-1`$,

```math
 \begin{aligned}
 T_{\mathcal L}&\le16r,&d_{\mathcal L}&\le16,
 &G_{\mathcal L}&\le2g_K+n,\\
 T_{\mathcal C}&\le32r+T_0,&d_{\mathcal C}&\le32+d_0,
 &G_{\mathcal C}&\le4g_K+4n+8+G_0.
 \end{aligned}                                           \qquad(17)
```

The last eight Cliffords are two literal X_a gates. Query calls are
sequential, so one mask pool suffices. Reuse the pool with the zero
reflection's borrowed work; every primitive has a complete arbitrary-input
guarantee, and its error is counted in (16). Including the physical
X-dual envelope, common source and actual-inverse bounds are

```math
 \begin{aligned}
 T_{U_A}&\le32r+T_0+2,&d_{U_A}&\le34+d_0,\\
 G_{U_A}&\le G_{\mathcal C}+17+4c_n+4n,&
 b_{U_A}&\le2n+\max(r,h_0),\qquad a_{\rm clean}=1.
 \end{aligned}                                           \qquad(18)
```

The 17 Cliffords include the signal Z and both rotation orientations.
The two principal carriers and address are excluded from b; Y,Z and
every private helper are included. At fixed eta this is constant
native T-depth with polynomial count/width overhead in m. It is a
paid correctly charged partial source, but A_x is a proper subprojector
of the required half-space. Recovery cannot enlarge its support.

<a id="two-reflection-bell-source"></a>

## 8. Two reflections give a complete small case and an exact scaling failure

Change the acceptance circuit using a Hadamard on only Y:

```math
 A_x'=H_YA_xH_Y,\qquad
 V_x=(I-2A_x')(I-2A_x),\qquad
 U_{2,x}=H_a Z_a(H_Y\mathcal C_{A,x}H_Y)\mathcal C_{A,x}H_a.
                                                               \qquad(19)
```

Only a starts clean. All four banks and the private helpers remain
arbitrary. The ideal accepted block is the anti-Hadamard-test block

```math
 K_{2,x}=(I-V_x)/2.                                     \qquad(20)
```

There is a complete fiber calculation of its support. Fix X-basis
eigenlabels alpha, beta, gamma on C,D,Z. The summand of A_x on Y is
$`(-1)^{\ell\cdot s}X_Y(K_x^{\mathsf T}\ell)`$, where
$`s=\alpha+\beta+K_x\gamma`$. If
$`v_x\cdot(\alpha+\beta)=1`$, the terms ell and ell+v_x cancel,
so A_x and A_x' both vanish on that entire fiber.

Otherwise s belongs to the image of K_x. Its restriction to the fixed
domain $`e_0^\perp`$ is bijective onto that image, so there is a
unique mathematical label u with $`K_xu=s`$ and u_0=0. This solve
labels the proof's fibers; no inverse query is emitted. On the last
n minus one Y bits, put $`d=2^{n-1}`$ and

```math
 |\chi_u\rangle=d^{-1/2}
       \sum_{y\in\mathbb F_2^{n-1}}(-1)^{u\cdot y}|y\rangle.
```

In this even-charge fiber the actual two projectors are

```math
 A_x=|\chi_u\rangle\langle\chi_u|\otimes I_{Y_0},
 \qquad
 A_x'=|u\rangle\langle u|\otimes I_{Y_0},
 \qquad |\langle u|\chi_u\rangle|^2=\epsilon:=2/q.
                                                               \qquad(21)
```

Thus for every n at least two,
$`A_xA_x'A_x=\epsilon A_x`$. The two rays are independent. Their
span projector, on the complete four-bank space, is

```math
 \begin{aligned}
 P_{{\rm span},x}
   &={A_x+A_x'-A_xA_x'-A_x'A_x\over1-\epsilon},\\
 K_{2,x}^\dagger K_{2,x}
   &={2I-V_x-V_x^\dagger\over4}
     =(1-\epsilon)P_{{\rm span},x},\\
 P_{{\rm span},x}&\le\Pi_{\rm eq},\qquad
 \mathrm{rank}(P_{{\rm span},x})=2q^3.
 \end{aligned}                                           \qquad(22)
```

At n=2, d=2 and the two rays span the whole Y tail. Hence
$`P_{{\rm span},x}=\Pi_{\rm eq}`$ and (22) is exactly
$`\Pi_{\rm eq}/2`$. This is a genuine complete one-clean source
for p=1, including every auxiliary input and coherent address, not
the misleading fixed-auxiliary coincidence in (7).

For every n at least three there is an explicit missing even-charge
input, independent of x. Put C,D,Z in their all-plus states, leave
Y_0 arbitrary, and use

```math
 |\psi_{\rm tail}\rangle
   ={\,|e_1\rangle-|e_2\rangle\over\sqrt2}
```

on the Y tail, where e_1,e_2 are two different nonzero standard basis
labels. In this fiber u=0. The displayed vector is orthogonal to
both the zero vector and the uniform plus vector, so both A_x and
A_x' annihilate it, whereas the desired paired charge is plus.
Consequently

```math
 \left\|K_{2,x}^\dagger K_{2,x}-\Pi_{\rm eq}/2\right\|=1/2
 \qquad(n\ge3).                                        \qquad(23)
```

The upper bound follows from (22): on the span the discrepancy is
$`1/2-\epsilon`$, on its even-charge complement it is minus one
half, and on odd charge it is zero. The witness remains valid for
an arbitrary coherent address and arbitrary references. It also
shows why more signal operations using only these same two ideal
reflections cannot distinguish the missing even-charge vectors from
odd-charge vectors: both reflections are identity on both spaces.
Conjugating both Y and Z by Hadamards instead makes the stabilizer
families commute and supplies only additional thin intersections.

### Native precision, complete return and recovery

Use (14) twice, with actual emitted scalar zero words of error eta.
The two additional H_Y walls and the signal H,Z,H are Clifford.
With the physical envelope included, the common ledger is

```math
 \begin{aligned}
 T_{U_2}&\le64r+2T_0,&d_{U_2}&\le64+2d_0,\\
 G_{U_2}&\le2G_{\mathcal C}+6n+4+4c_n,&
 b_{U_2}&\le2n+\max(r,h_0),\qquad a_{\rm clean}=1,\\
 \|\widetilde U_2-U_2\|&\le2\eta,&
 \|\widetilde G_2-G_2\|&\le4\eta.
 \end{aligned}                                           \qquad(24)
```

Every native query helper returns exactly. All return errors of the
approximate scalar zero words, including any movement of their tested
signal/carrier wires between boundaries, are included in (24).
The actual source inverse reverses both emitted scalar zero words;
neither is replaced by a second forward approximation. On the
missing-even witness its accepted norm is at most 2 eta, so the
emitted lookup Gram error is at least $`1/2-4\eta^2`$.

For p=1, the scalar zero test has only three tested bits. Use an exact
eight-T returned CCZ with negative C literals, rather than an
approximation. Each controlled reflection then uses four returned
bilinear query words and that CCZ. The explicit fixture emits the
whole source with 80 T gates and 168 Cliffords in dual coordinates;
the common positive-S actual inverse uses the same bounds after
adjacent S powers are reduced modulo four. With the safe fixed basis
envelope this gives

```math
 T_{U_2}\le80,\quad G_{U_2}\le216,\quad d_{U_2}\le80,
 \quad b_{U_2}\le5,
 \qquad(p=1,\ n=m=2).                                  \qquad(25)
```

The source Gram error is zero. The already proved
[four-call recovery and serial dirty promotion](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
therefore give a complete exact zero-clean p=1 indicator, for example
with the simultaneous safe bounds
$`T\le2664`$, $`G\le7276`$, $`D_T\le2592`$ and dirty width
at most 11. All auxiliary return is part of that theorem. If using
the approximate scalar primitive instead, eta at most
$`\tau^2/16384`$ makes the p=1 Gram error at most
$`\tau^2/4096`$ and the same recovery has final error at most tau.
For larger p no precision choice repairs (23).

The [short-word thin-reflection rank bounds](EXPANDER_ZERO_REFLECTION.md#11-a-paid-graph-reflection-still-needs-complete-readout)
already cover the relevant baseline-dependent call obstruction;
they are not a new general depth bound here. Their exact controlled-call
specialization is useful. Let D be the dimension of the banks, address
and arbitrary helpers, excluding the signal. Each ideal controlled
$`R_{A,x}`$ call differs from identity on rank D/q in the full
2D-dimensional space. For k such calls with arbitrary unitary
interlayers, replace every call by identity to define the baseline
accepted block K_0. Telescoping the rank defects and compressing to
the signal-zero subspace gives

```math
 \mathrm{rank}(K-K_0)\le kD/q.
                                                               \qquad(26)
```

If K_0=0 and k is strictly less than q/2, the kernel of K intersects
the desired rank-D/2 even-charge space. Hence the ideal Gram error is
at least one half. If the emitted whole word is within full operator
error nu of this ideal word, its accepted norm on that kernel vector
is at most nu, giving Gram error at least $`1/2-\nu^2`$.
Thus a sufficiently accurate repair within this zero-accepted-baseline
controlled-call grammar requires at least $`q/2=2^{n-1}`$ ideal
calls. This constant counts controlled reflections: an uncontrolled
reflection on the two-signal-sector space has twice the rank defect.
The baseline assumption is essential, and this is not a native
T-depth bound.

The constructive progress is (14)--(25): controlled and normalized thin sources are paid, a real
complete small case is identified, and the next unpaid operation is
removing the common-kernel blindness on the full even-charge sector.
Reopen with a changed broad phase or signal coupling that acts on that
space, not another precision choice or scalar function of this walk.

The additional exact fixtures check the full p=1 Gram, the literal
source and actual inverse with both arbitrary private-mask values,
and the p=2 common-kernel witness. These support the analytic proof;
they are internal verification, not external mathematical review.

<a id="four-path-bell-relative-phase"></a>

## 9. A broad four-path Bell source has the wrong coherent sectors

A changed trial coherently combines a swap and a wide non-Clifford phase,
without a thin outcome projector. Retain the complete two-copy Bell word
$`W_{2,x}=W_{x,CY}W_{x,DZ}`$, with all four n-bit banks arbitrary.
On the two auxiliary banks define

```math
S=\mathrm{SWAP}_{YZ},\qquad
B=T_Y^{\otimes n}(T_Z^\dagger)^{\otimes n},\qquad
B|y,z\rangle=\omega^{|y|-|z|}|y,z\rangle,
\qquad SBS=B^\dagger.                                  \qquad(27)
```

Here T_Y denotes literal T gates on bank Y, not a Pauli-Y-axis rotation.
The weight difference is an ordinary integer. Use two clean signals a,b,
and emit chronologically

```math
W_{2,x}^\dagger,\ H_a,\ \mathcal C_a(S),\ H_a,
\mathcal C_a(B),\ H_a,\ W_{2,x},\ H_b.                 \qquad(28)
```

Both signals are accepted jointly at zero. The second signal gives a
paid attenuation, while the first supplies four interfering paths. Every
inverse in (28) is the actual native reverse. Direct signal multiplication
and (27) give the **complete** accepted block and Gram:

```math
\begin{aligned}
K_{\rm acc,x}&=\frac14W_{2,x}(I+S+B-BS)W_{2,x}^\dagger,\\
G_x&=W_{2,x}\left[\frac I4+\frac{(B^\dagger-B)S}{8}\right]
       W_{2,x}^\dagger,\qquad 0\preceq G_x\preceq I/2.
                                                               \qquad(29)
\end{aligned}
```

For the Gram calculation, write the numerator as
$`(I+B)+(I-B)S`$. Its two squared terms sum to 4I, and each cross
term equals $`(B^\dagger-B)S`$. This last operator is Hermitian,
because S exchanges B and its adjoint. No auxiliary label is measured,
initialized, averaged or erased. Formula (29) includes the off-diagonal
coherences between swapped outcomes.

### Native parallel controls and the simultaneous ledger

The controlled swap is the retained private-phase Fredkin batch from
Section 5: 8n T gates, 12n Cliffords, eight T layers and n returned
private helpers. Shared signal a is only a Clifford control; each row's
physical T target is private, and its bank targets are disjoint.

The [literal controlled-T word](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#literal-controlled-t-wall)
has 70 T gates, 88 Cliffords, at most 70 T layers and four arbitrary
returned helpers. **All its physical T targets are private.** Its two
logical controls, including a, occur only as Clifford controls. Give each
of the 2n controlled T or T-dagger factors in (27) its own four-helper
pool. The retained word and its actual inverse therefore align in the
same 70 T layers, with 140n T gates and 176n Cliffords. This proves a
native shared-control schedule; no standard controlled-T decomposition
with a shared signal T target is being parallelized.

The two Bell copies also run in parallel: they have disjoint carriers
and private pools, and their native query words use the common address
only as Clifford controls. Each outer two-copy wall has depth at most
eight and uses 2r returned query helpers. The two outer walls, Fredkin
batch and controlled-phase batch are serial. Reuse their private pools
only after their certified complete return.

With the notation and physical X-dual envelope of Sections 1 and 5,
one common forward/actual-inverse ledger for (28) is

```math
\begin{aligned}
T_U&\le32r+148n,&D_T(U)&\le94,\\
G_U&\le4G_W+192n+4c_n+4,
&b_U&\le2n+\max\{2r,8n\},\qquad a_{\rm clean}=2.
                                                               \qquad(30)
\end{aligned}
```

The 2n term in b_U is the arbitrary Y,Z auxiliary allocation; the two
principal carriers and address are excluded. The four signal Hadamards
are included. The physical envelope contributes $`4c_n+4n`$
Cliffords, and the two Bell walls contain four W calls in total.
All components are exact native words. If the four W occurrences
instead have full error at most eta, retaining actual inverse partners,
source error is at most 4eta and Gram error at most 8eta. These are
full-space errors, including every private return and reference.
The simultaneous costs fit the relaxed polynomial source budgets;
no unpaid normalization or simultaneous-control step remains.

### The entire target sector, not a single outcome, is required

The Bell kernel fixes the target's orientation exactly. Equation (4),
or direct Pauli propagation through (3), gives

```math
W_x^\dagger X_A(\ell)W_x
 =Z_A(K_x^{\mathsf T}\ell)Z_B(\ell),\qquad
W_{2,x}^\dagger Q_xW_{2,x}=Z_Y(v_x)Z_Z(v_x).            \qquad(31)
```

The second identity uses $`K_x^{\mathsf T}v_x=0`$. Thus in the
coordinates inside the W walls, the desired Gram is
$`(I+Z_Y(v_x)Z_Z(v_x))/4`$, with identity on both data carriers.

Equal-weight auxiliary outcomes have B=I and actual Gram I/4. The
vacuum y=z=0 is in the correct charge sector, where the desired value
is I/2. This already gives error at least one quarter for every size.
There is a sharper coherent witness for every canonical $`n=pm\ge2`$.
Choose distinct auxiliary coordinates i,j, put $`y=e_i+e_j`$, z=0,
and use the ordered pair of states $`|y,0\rangle,|0,y\rangle`$.
On their invariant span the exact Gram from (29) is

```math
\frac14\begin{pmatrix}1&-i\\i&1\end{pmatrix},       \qquad(32)
```

with eigenvalues zero and one half. Both states have the same target
charge $`(-1)^{v_x\cdot y}`$. The desired Gram on their whole span
is therefore either zero or I/2; it cannot assign the two different
eigenvalues in (32). Since both complete Grams lie between zero and
I/2, this proves the exact distance

```math
\left\|G_x-\Pi_{\rm eq}/2\otimes I_{YZ}\right\|=\frac12
\qquad(p\ge1),                                        \qquad(33)
```

and at least $`1/2-8\eta`$ for the stated approximate implementation.
The witness is coherent in the auxiliary banks. Discarding the off-diagonal
entries of (32) would lose the full error. Restriction to any fixed address
is enough for the operator-norm bound on coherent addresses and references.

Pause this particular broad four-path source. Its normalization, native
controls, actual inverses and arbitrary-input action are paid, but the
accepted eigenvectors distinguish auxiliary swap coherence instead of the
required complete lookup charge. More precision or the existing recovery
cannot change (33). A changed auxiliary coupling must first supply that
charge-dependent action; arbitrary multipath sources and different
entangling interlayers are outside this calculation. No complete-indicator
or frame-depth frontier changes.

The [exact multipath fixtures](../../tests/test_multipath_bell_source.py)
check a literal one-cell control word with 148 T gates and 192 Cliffords,
all principal/signal columns at two dirty-helper assignments, the actual
inverse, and the accepted four-path coefficients. Separate exact checks
verify the complete two-cell Gram, its coherent eigenvalue witness, and
small native Bell conjugations fixing the transpose in (31). They are
finite internal checks, not a general native emitter or external
mathematical review; the all-size operator and scheduling claims above
are analytic.

<a id="bell-acceptance-charge-transfer"></a>

## 10. A complete Bell envelope transfers the same charge to acceptance

The previous failures are not repaired merely by allowing arbitrary broad,
non-diagonal auxiliary acceptance after the same Bell envelope. The charge identity (31) gives
an exact consequence for arbitrary acceptance effects. It does not exclude an
address-dependent acceptance that actually computes the missing charge,
or a source with genuinely interleaved data and auxiliary dynamics.

Keep the notation of Sections 1--3 and put

```math
W_{2,x}=W_{x,CY}W_{x,DZ},\qquad
Q_x=X_C(v_x)X_D(v_x),\qquad
R_x=Z_Y(v_x)Z_Z(v_x).
```

Here Y,Z are entire arbitrary auxiliary banks. Equation (31) gives
$`W_{2,x}^\dagger Q_xW_{2,x}=R_x`$ using the left kernel $`v_x`$,
not the known right kernel $`e_0`$. It includes identity on every other
principal and implementation register.

Let an arbitrary acceptance mechanism on x,Y,Z and additional work have
effect E, a positive contraction. Neither diagonality nor a chosen
auxiliary input is assumed. Extend E by identity on C,D and let
$`W_2=\sum_x|x\rangle\langle x|\otimes W_{2,x}`$. Suppose the
complete source Gram has the envelope form

```math
G=W_2(I_{CD}\otimes E)W_2^\dagger.
```

The notation $`I_{CD}\otimes E`$ permits E to act on or even move the
address; it only requires that acceptance itself act trivially on C,D.
Define the coherent direct sums Q and R from $`Q_x`$ and $`R_x`$. Unitary
invariance of operator norm and (31) give the exact equality

```math
\left\|G-\frac{I+Q}{4}\right\|
 =\left\|I_{CD}\otimes E-\frac{I+R}{4}\right\|.      \qquad(34)
```

All target factors are tensored with identity on arbitrary extra work.
Thus the same complete lookup charge must be resolved on the auxiliaries,
with exactly the same error. Coherent off-diagonal auxiliary effects and
arbitrary accepted unitary fibers do not change (34). It is an equality
of full operator errors, not a fixed-input probability or a normalized
trace comparison. Equation (34) remains valid on an external reference.

There is a sharp consequence when E is address independent. For distinct
addresses $`x,x'`$, the Vandermonde columns $`v_x,v_{x'}`$ are distinct. Choose a
binary auxiliary difference y+z on which their two characters disagree.
Then

```math
\left\|\frac{I+R_x}{4}-\frac{I+R_{x'}}{4}\right\|=\frac12.
```

The triangle inequality forces

```math
\max_x\left\|E-\frac{I+R_x}{4}\right\|\ge\frac14.
                                                               \qquad(35)
```

The positive contraction $`E=I/4`$ attains one quarter simultaneously
at every address, so the constant is sharp for this effect class. This
holds already at p=1. If acceptance is address dependent, (34), rather
than (35), is the retained statement: implementing its required effect
is still an unpaid source operation.

No new native gate is assumed. $`W_2`$ is the already paid two-copy Bell
envelope; E and any mechanism realizing it retain their own count,
depth, width, precision and return obligations. In particular an
abstract effect satisfying (34) is not a source construction. The
equivalence does not cover a different envelope whose data and auxiliary
couplings are interleaved with acceptance.

<a id="hyperplane-permutation-acceptance-gap"></a>

### Pure charge-preserving permutation expansion has two fixed coset vectors

A separate proposed use of the same kernel is to expand within
$`H_x=\mathrm{im}K_x=\ker(v_x^{\mathsf T})`$. This requires a
different action from merely preserving charge. Let $`e=e_0`$, so
$`v_x^{\mathsf T}e=1`$, and define the two normalized vectors

```math
|u_{x,\sigma}\rangle
 =2^{-(n-1)/2}\sum_{y:v_x^{\mathsf T}y=\sigma}|y\rangle,
\qquad \sigma\in\{0,1\}.
```

Every computational-basis permutation preserving both cosets fixes
both vectors exactly. This includes affine maps y to Ay+c with
$`v_x^{\mathsf T}A=v_x^{\mathsf T}`$ and
$`v_x^{\mathsf T}c=0`$. In particular it includes each transvection

```math
A=I+(K_xa)b^{\mathsf T},\qquad b^{\mathsf T}K_xa=0.
```

The last condition makes A an involution, and the left-kernel identity
makes it charge preserving. These assertions do not assume a native
implementation of arbitrary choices of a,b.

Consider any complete word of such permutations controlled by signals
or other work, their actual inverses, and arbitrary interlayers acting
as identity on this tested bank. The address is held fixed between
these stages; interlayers may be address controlled but may not mix
different x values. Helpers and acceptance have the same insertion and
projection for both bank states. On either uniform coset vector every
controlled permutation acts as identity, so the induced signal/work
word is exactly the same for sigma=0 and sigma=1. Its accepted Gram
has equal expectations on these two bank states and any common work
input.

The desired $`(I+Z(v_x))/4`$ instead has expectations one half and
zero. Therefore the complete accepted Gram error is at least one
quarter, regardless of how many permutation calls or selectors are
used, or what spectral gap an average of those permutations has. For
a complete unitary phase implementation the corresponding error is
at least $`\sqrt2`$: its common induced work unitary V obeys

```math
\|(V-I)\psi\|^2+\|(V+I)\psi\|^2=4
```

for every normalized work input psi. At least one of the two required
charge sectors therefore has error at least square root two.

This is a restriction on the specified permutation grammar. Bank phases,
other bank dynamics, or intermediate address motion can leave it. A
constant gap for motion within the cosets, by itself, would select
invariant orbit vectors and would not supply the full desired charge
sector. The next operation must distinguish the sectors on all bank
columns; charge preservation and expansion do not establish that readout.

The exact [acceptance-transfer fixtures](../../tests/test_bell_acceptance_transfer.py)
use literal p=1 Bell matrices on every data and auxiliary column, a
non-diagonal coherent auxiliary effect, p=1,2 distinct-address charge
comparisons, and finite hyperplane transvections. They check the native
orientation and complete-space identities; the all-size reductions and
scope statements above are analytic, not inferred from finite fixtures.

## Clifford projector-Gram results

<a id="interleaved-clifford-gram-rigidity"></a>
<a id="11-interleaved-clifford-fibers-still-require-a-growing-query-word"></a>

[Interleaved Clifford projector-Gram theorem](CLIFFORD_PROJECTOR_GRAM.md#interleaved-clifford-gram-rigidity).

<a id="stabilizer-compression-proves-exact-gram-rigidity"></a>

[Stabilizer-compression proof](CLIFFORD_PROJECTOR_GRAM.md#stabilizer-compression-proves-exact-gram-rigidity).

<a id="a-query-count-consequence-for-bounded-clean-allocation"></a>

[Bounded-clean affine-query limit](CLIFFORD_PROJECTOR_GRAM.md#a-query-count-consequence-for-bounded-clean-allocation).
