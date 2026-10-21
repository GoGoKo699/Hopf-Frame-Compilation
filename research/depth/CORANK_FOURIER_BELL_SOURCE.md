# Corank-one Fourier kernels and complete Bell acceptance

[Kernel translation query](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-xor-carry-source) ·
[Paired Gram contract](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) ·
[Four-call recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

The paid corank-one translation matrix has an exact singular Fourier
operator whose left Gram is the desired central-character projector.
Its complete two-bank unitary dilation is also paid. This does not
supply a one-clean source: selecting a Bell fiber uses an entire
auxiliary bank, and the paired shared-outcome construction has success
weight exponential in minus the carrier width. A concrete coherent-SWAP
repair instead has a fully normalized one-clean Gram and maximal lookup
error at every size. These are statements about the displayed words,
not general native-depth obstructions. The complete-frame frontier is
unchanged.

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
Even if a source for $`A_x/2`$ were supplied, its Gram distance from
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
a different acceptance circuit. The missing operation is a paid
coherent extraction of its left central character with constant
accepted amplitude and arbitrary auxiliary inputs. Shared Pauli
byproducts alone give (7), while the concrete coherent swap gives
(8), not that character. Reopening this route requires a changed
unitary acceptance word, with its complete Gram derived; another
claim that corank one supplies a free one-clean dilation does not
change the source status.
