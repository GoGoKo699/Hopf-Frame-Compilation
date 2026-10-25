# Clifford interlayers and complete lookup

[Weighted pulses](COORDINATE_WEIGHTED_PULSE_LIMITS.md) ·
[Complete derivative](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) ·
[Central transport](COHERENT_CENTRAL_TRANSPORT.md) ·
[Dirty cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md)

A growing fixed Clifford network can destroy the coordinate-signature witness.
Nevertheless two nonlinear blocks remain insufficient in a broader, basis-free
sense: the first block may be any third-level Clifford-hierarchy unitary, and
the second may be any diagonal third-level unitary surrounded by arbitrary
Cliffords. The complete lookup derivative has a uniform error floor in this
class, including arbitrary work and address movement.

This is a two-block obstruction. It does not give a growing lower bound on
unrestricted native T-depth, exclude three such blocks, or change the general
frame frontier.

## 1. Complete model and claim

Let p be at least three and let the coefficient bank have $`m=2^p`$ independent
bits $`c_S`$, one for every subset S of the address labels. Define

```math
f(x,c)=\sum_{S\subseteq[p]}c_S\prod_{i\in S}x_i\pmod2,
\qquad Q|x,c,h\rangle=(-1)^{f(x,c)}|x,c,h\rangle.          \qquad(1)
```

The work h consists of arbitrarily many arbitrary input qubits. Q acts as
identity on that work. The fixed Möbius bank Clifford converts this coefficient
lookup into ordinary addressed lookup, so the statement is unchanged by that
already paid coordinate choice. All errors are complete operator norms,
including arbitrary external references; no work input is initialized.

Write $`\mathcal C_3`$ for the third level of the Clifford hierarchy: U belongs
to this level when $`UPU^\dagger`$ is a Clifford for every Pauli P. Allow sources

```math
V=C D U,\qquad C\text{ Clifford},\quad
D\text{ diagonal and in }\mathcal C_3,\quad U\in\mathcal C_3.
                                                               \qquad(2)
```

For any Hermitian Pauli probe P, the complete derivative obeys

```math
\boxed{\ \|VPV^\dagger-Q\|\ge\sin(\pi/8).\ }            \qquad(3)
```

Any two weighted quarter-angle parity pulses with arbitrary initial, middle
and final Clifford networks fall within (2). Each pulse has the form
$`H_JD_jH_J`$ with a diagonal weighted cubic $`D_j`$, and Clifford multiplication
on either side preserves $`\mathcal C_3`$. The class also contains every source of native
T-depth at most two. Conversely, a $`\mathcal C_3`$ block need not have native depth one;
(3) therefore concerns a broader block grammar, not an identification of
block count with native depth. Arbitrary clean initialized inputs are outside
the complete-space premise.

## 2. Two elementary Clifford thresholds

Put

```math
\delta=2\sin(\pi/8)=\sqrt{2-\sqrt2}.
```

**Diagonal-sign threshold.** If a computational-diagonal sign unitary R is
within strict distance delta of a Clifford K, even allowing an arbitrary
common scalar on K, then R is itself a Clifford up to a common scalar.

A Clifford column is a stabilizer state, flat on an affine support. If its
support has more than one element, its amplitude at the target basis state
has modulus at most $`1/\sqrt2`$, giving column distance at least delta. A wrong
monomial permutation gives an orthogonal output column and distance $`\sqrt2`$.
Thus a closer K must be diagonal. Its entries have the form $`\lambda i^{q(z)}`$,
where q is quadratic modulo four with even cross coefficients. Dividing by
the target signs gives values in lambda times the fourth roots of unity.
The strict condition $`|r-1|`$ less than delta places every ratio in the open
angular interval $`(-\pi/4,\pi/4)`$. That interval cannot contain two distinct
rotated fourth roots. All ratios are therefore one common scalar, proving
the threshold.

**Hermitian monomial quantization.** Every nonzero entry of a Hermitian
monomial Clifford belongs literally to the eighth roots of unity. Write
its action as $`\lambda i^{q(z)}`$ times an affine basis permutation pi. Its
square is identity, so

```math
\lambda^2 i^{q(z)+q(\pi(z))}=1.
```

Hence lambda squared is a fourth root, lambda is an eighth root, and so is
every nonzero entry. The minimum distance between distinct eighth roots is
delta. For a diagonal D in $`\mathcal C_3`$, the operator $`DX^aD^\dagger`$ is itself a Hermitian
monomial Clifford. Its nonzero entries are exactly the relative phases
of D at pairs of computational inputs separated by a. The quantization
just proved therefore makes every such relative phase an eighth root,
without a choice of global phase for D. Consequently diagonal $`\mathcal C_3`$
conjugation of a Hermitian monomial Clifford still has literal eighth-root
entries.

## 3. A good Pauli cannot translate the address

For a Hermitian unitary R, call a Pauli label good if RPR is a Clifford.
If P has a nonzero address-X translation a, the complete lookup gives

```math
QPQP^\dagger|x,c,h\rangle
 =(-1)^{f(x,c)+f(x+a,c+b)}|x,c,h\rangle,                 \qquad(4)
```

where b is P's bank translation. Its Z factors, scalar phase and helper
translation cancel from (4). In the coefficient of the independent bank
bit $`c_{[p]}`$, the exponent contains

```math
\prod_{i=1}^p x_i+\prod_{i=1}^p(x_i+a_i).
```

For any j with $`a_j`$ equal to one, the coefficient of the degree-(p-1)
monomial over all i other than j is one. Multiplication by $`c_{[p]}`$ therefore
produces a nonzero degree-p Boolean term. No other bank bit or bank-translation
term can cancel it. Since p is at least three, (4) is not a diagonal Clifford:
a diagonal sign Clifford has Boolean degree at most two. Thus every good
Pauli for Q has zero address-X translation and commutes with all address Z's.

This conclusion uses the complete programmable bank. It is not a degree
lower bound for arbitrary native circuits or for a single scalar Boolean
function computed with arbitrary dirty work.

## 4. Approximation forces an address-preserving final frame

Let $`A=VPV^\dagger`$ and suppose, for contradiction, that

```math
\|A-Q\|\lt\delta/2.
```

With $`K=UPU^\dagger`$, K is a Hermitian Clifford and
$`A=CDKD^\dagger C^\dagger`$. Consider the maximal commuting Pauli label space

```math
\mathcal L=C\langle Z_1,\ldots,Z_n\rangle C^\dagger,
```

where n counts every address, bank and work wire. For every Hermitian Pauli
R in this space, $`C^\dagger RC`$ is Z type and commutes with D. Therefore

```math
ARA=C D\bigl(K(C^\dagger RC)K^\dagger\bigr)D^\dagger C^\dagger
```

is a Clifford: the parenthesized operator is a Pauli, and D belongs to $`\mathcal C_3`$.
Unitary telescoping gives

```math
\|QRQ-ARA\|\le2\|A-Q\|\lt\delta.
```

Multiplication by $`R^\dagger`$ compares the diagonal sign $`QRQR^\dagger`$ with a
Clifford. The diagonal-sign threshold makes QRQ a Clifford exactly. Thus
all of L consists of good labels for Q. By Section 3 every address Z commutes
with L. A maximal commuting Pauli label space equals its own symplectic
centralizer, so every address Z belongs to L.

Equivalently, the p independent Paulis $`C^\dagger Z_{x_i}C`$ are all diagonal.
A fixed monomial Clifford B can map those independent Z labels to the
standard address labels, including correcting their signs with X gates:

```math
B^\dagger C^\dagger Z_{x_i}CB=Z_{x_i}\quad\text{for all }i.
```

This is only a change of coordinates in the proof, not a free native resource.
Put $`C'=CB`$, $`D'=B^\dagger DB`$ and $`U'=B^\dagger U`$. Then

```math
V=C'D'U',\qquad C'^\dagger Z_{x_i}C'=Z_{x_i}.
```

The monomial Clifford B preserves the diagonal $`\mathcal C_3`$ phase class by affine
substitution. Clifford multiplication preserves $`\mathcal C_3`$, so D' is diagonal $`\mathcal C_3`$
and U' is $`\mathcal C_3`$. The final Clifford C' now preserves every address value.

## 5. Monomial rounding closes the contradiction

Define

```math
Q'=C'^\dagger Q C',\qquad K'=U'PU'^\dagger,\qquad
A'=D'K'D'^\dagger.
```

The error remains less than delta/2. Because C' holds the address, Q' is a
bank/work Pauli in every fixed-address fiber. One can see this without an
assumption about intermediate fibers: for each fixed bank Z(v), the global
Clifford conjugate $`C'^\dagger Z(v)C'`$ has no address-X component, so restriction
to an address gives a Hermitian Pauli on bank/work. Apply this to the bank
label selected by that address in (1). Therefore Q' is computational-monomial
and all its nonzero column phases belong to the fourth roots of unity.

Diagonal D' cannot change a column support or modulus. If K' has a
nonmonomial column, its overlap with the single target output has magnitude
at most $`1/\sqrt2`$, giving distance at least delta. If K' has a wrong
monomial permutation, the distance is $`\sqrt2`$ on that column. Neither is
consistent with error less than delta/2. Hence K' is monomial with exactly
the same permutation as Q'.

By Hermitian monomial quantization, every nonzero entry of A' is an eighth
root. Every nonzero entry of Q' is a fourth root. Distinct such entries
are at distance at least delta. The strict error bound therefore forces
literal equality of every entry:

```math
Q'=D'K'D'^\dagger.                                      \qquad(5)
```

The right side is itself $`\mathcal C_3`$. Indeed K' is a monomial Clifford, so
$`K'D'^\dagger K'^\dagger`$ is diagonal $`\mathcal C_3`$. Diagonal $`\mathcal C_3`$ phases form a group:
for each binary translation a, the multiplicative difference of a product
is the product of the two multiplicative differences, which are diagonal
Cliffords. Thus $`D'(K'D'^\dagger K'^\dagger)K'`$ belongs to $`\mathcal C_3`$. Clifford conjugation
preserves $`\mathcal C_3`$, and (5) would make Q belong to $`\mathcal C_3`$ as well.

This is impossible. The diagonal sign Q in (1) has the independent term
$`c_{[p]}`$ times the product of all p address bits, of degree p+1 at least four.
A diagonal $`\mathcal C_3`$ sign has Boolean degree at most three: its difference in
each individual input bit must be a diagonal Clifford sign of degree at
most two. This contradiction
proves (3).

If an actual native circuit approximates a source in (2), tensored with its
promised returned work, within epsilon in complete operator norm, its lookup
derivative error is at least $`\sin(\pi/8)-2\varepsilon`$. The complete-source
assumption includes all work-return errors. The theorem itself already allows
arbitrary work and intermediate address motion; only initialized-work promises
fall outside it.

## 6. A long physical Clifford frame and its calibrated repair

The following candidate gives every coordinate its own place in a growing
Clifford network. It therefore tests a different mechanism from choosing a
fixed small list of coordinate masks. Define the quarter-angle pulses

```math
\begin{aligned}
E&=\prod_{|S|\ {\rm even},\ i\notin S}
      e^{i\pi X_{A_S}n_{S+i}/4},&
E_x&=\prod_{|S|\ {\rm even},\ i\notin S}
      e^{-i\pi Z_{x_i}X_{A_S}n_{S+i}/4},\\
O&=\prod_{|S|\ {\rm odd},\ i\notin S}
      e^{i\pi X_{A_S}n_{S+i}/4},&
O_x&=\prod_{|S|\ {\rm odd},\ i\notin S}
      e^{-i\pi Z_{x_i}X_{A_S}n_{S+i}/4}.
\end{aligned}                                                   \qquad(6)
```

Here $`S+i=S\cup\{i\}`$ and $`n_T=(I-Z_{A_T})/2`$. For each
coordinate, take the **physical** full-coordinate Clifford

```math
C_i=\prod_{S:i\notin S}e^{i\pi X_{A_S}n_{S+i}/2},\qquad
C=C_{p-1}\cdots C_0,\qquad V_x=E_x C O_x C^\dagger E.
                                                               \qquad(7)
```

Each $`C_i`$ includes both tail parities. Its elementary edges are the
literal products $`S_{S+i}\,\mathrm{CX}_{S+i\to S}`$; their phases remain present. The binary action
of C is the complete fixed subset-zeta map, which can lower coefficient
labels through arbitrarily many ranks at zero T-depth. This is the
intended use of a growing Clifford network. On the face with only the
empty and singleton coefficients occupied, the middle odd pulse is
identity and C commutes with the even pulses; the source gives the
correct affine coefficient action there.

Nevertheless, for every $`p\ge2`$,

```math
\|V_x Z_{A_\varnothing}V_x^\dagger-Q\|\ge\sqrt3.         \qquad(8)
```

### The first failure is zero-address calibration on a rank-two face

Restrict to any coordinate pair $`i\lt j`$, with bank bits
$`a=A_\varnothing,b=A_i,c=A_j,d=A_{ij}`$. Set every other coefficient
and every address bit to zero. This is an allowed invariant input
subspace: every exterior edge has a zero head. All exterior coordinate
Cliffords act identically, and the surviving ordered frame is
$`C_jC_i`$. Its exact action is

```math
C|a,b,c,d\rangle
 =i^{b+(c\mathbin\oplus d)+2d}
 |a\mathbin\oplus b\mathbin\oplus c\mathbin\oplus d,
   b\mathbin\oplus d,c\mathbin\oplus d,d\rangle.
                                                               \qquad(9)
```

The exponent is an integer exponent, not a binary phase. Fix d=1,
which remains invariant. Every factor commutes with $`X_a`$. In its
eigenspace $`X_a=\lambda\in\{+1,-1\}`$, put
$`D_\lambda=\mathrm{diag}(1,i\lambda)`$. Equation (9) gives

```math
\begin{aligned}
C_\lambda&=-i\lambda
 (XD_\lambda)_b\otimes(XD_{-\lambda})_c,\\
C_\lambda X_b C_\lambda^\dagger&=-\lambda Y_b,&
C_\lambda X_c C_\lambda^\dagger&=+\lambda Y_c.
\end{aligned}                                                   \qquad(10)
```

At zero address, $`E_x=E^\dagger`$ and $`O_x=O^\dagger`$.
Writing $`F_\lambda=e^{i\pi\lambda(n_b+n_c)/4}`$ gives

```math
\begin{aligned}
(V_0)_\lambda
 &=F_\lambda^\dagger
       e^{i\lambda\pi(Y_b-Y_c)/4}F_\lambda
   =B_\lambda\otimes B_\lambda^\dagger,\\
B_\lambda&=\frac{I}{\sqrt2}+\frac{i}{2}(X+\lambda Y).
\end{aligned}                                                   \qquad(11)
```

In this basis $`Z_a`$ exchanges the two lambda sectors, and the target
at zero address is $`Z_a`$. The derivative error on this block is

```math
\|R-I\|,\qquad
R=(B_+B_-^\dagger)\otimes(B_+^\dagger B_-).
```

Each two-by-two factor is special unitary with trace one, hence has
eigenvalues $`e^{\pm i\pi/3}`$. Thus
$`\mathrm{spec}(R)=\{1,1,e^{2\pi i/3},e^{-2\pi i/3}\}`$,
and the block error is square root three. The invariant embedding proves
(8) for every p. The failure already occurs at zero address; the correct
action on affine coefficient faces does not supply complete calibration.

### Paying the zero-address repair still leaves maximal lookup error

Use the actual inverse of the zero-address source to calibrate it. The
completed product simplifies literally to

```math
W_x=V_xV_0^\dagger=E_x C O_x O C^\dagger E,
\qquad W_0=I.                                        \qquad(12)
```

This repair pays the missing source factors; it is not a free calibration.
Now choose the same invariant two-coordinate face with
$`x_i=x_j=1`$ and d=1. On this block,

```math
O_xO=-X_bX_c,\qquad C_\lambda(O_xO)C_\lambda^\dagger=Y_bY_c.
```

Also $`E_x=E`$. Since $`Y_bY_c`$ sends $`n_b+n_c`$ to
$`2-n_b-n_c`$, (10) gives

```math
(W_x)_\lambda=F_\lambda Y_bY_c F_\lambda
 =i\lambda Y_bY_c,
\qquad W_x=iX_aY_bY_c.                               \qquad(13)
```

Its transported probe is $`-Z_a`$. The required target on this block
is $`-Z_aZ_bZ_c`$, which equals $`+Z_a`$ on b=1,c=0. Therefore

```math
\|W_x Z_{A_\varnothing}W_x^\dagger-Q\|=2
\qquad(p\ge2).                                      \qquad(14)
```

The upper bound two holds for the difference of unitaries; the displayed
allowed input supplies the matching lower bound. Equations (9)--(13)
retain all physical Clifford phases. Replacing C by its binary zeta
permutation would be a different candidate.

### Angle tuning does not repair this calibrated family

Keep the same physical C and replace the even quarter angle by a real
alpha and the odd quarter angle by a real beta:

```math
W_x(\alpha,\beta)
 =E_x(\alpha)C O_x(\beta)O(\beta)C^\dagger E(\alpha).
```

This family is still identity at zero address. On the same block
$`x_i=x_j=1,d=1`$, put $`c=\cos(2\beta)`$. Its two root-X blocks
have the exact form

```math
\begin{aligned}
W_\lambda&=e^{2i\alpha\lambda}A_\lambda\otimes B_\lambda,\\
A_\lambda&=c\cos\alpha\,I
 -i\lambda\bigl(\sin(2\beta)Y+c\sin\alpha\,Z\bigr),\\
B_\lambda&=c\cos\alpha\,I
 +i\lambda\sin(2\beta)Y-i\lambda c\sin\alpha\,Z.
\end{aligned}
```

Both A and B are special unitary and their minus-lambda blocks are
their adjoints. Put $`u=2\alpha`$ and choose
$`\gamma=\arcsin[-c^2\sin u]`$. The matrices $`ZA_+^2`$ and
$`ZB_+^2`$ each have determinant minus one and trace
$`-2ic^2\sin u`$. Comparing the transported root probe with the
target $`-Z_aZ_bZ_c`$ therefore gives the relative spectrum

```math
\{e^{2iu},e^{2iu},-e^{i(2u+2\gamma)},-e^{i(2u-2\gamma)}\}.
```

The block error is
$`\max\{2|\sin u|,2|\cos(u+\gamma)|,2|\cos(u-\gamma)|\}`$.
It is at least square root two. Indeed, if $`y=\sin^2u\ge1/2`$,
the first term suffices. Otherwise $`z=\sin^2\gamma\le y`$ and
the mean of the two squared cosine terms is
$`1-y-z+2yz\ge1-2y+2y^2\ge1/2`$. Thus for every pair of real
angles and every $`p\ge2`$,

```math
\|W_x(\alpha,\beta)Z_{A_\varnothing}
       W_x(\alpha,\beta)^\dagger-Q\|\ge\sqrt2.
```

This is an analytic statement about the specified calibrated family.
No native synthesis cost is asserted for arbitrary real angles; the
ledgers below concern the original quarter angles only.

### Complete native ledgers

Put $`r=pm/4`$. The frame C has $`pm/2`$ edges, each costing one S
and one CX, so its forward Clifford count is pm. Its actual inverse
costs at most 2pm with $`S^\dagger=S^3`$. The pair of occurrences costs
at most $`3pm=12r`$ in either orientation. Compiling the three quarter
pulses separately gives common forward/actual-inverse bounds

```math
T(V)\le3m+96r,\quad D_T(V)\le39,\quad
w(V)\le5r,\quad G_{\rm all}(V)\le18m+225r.             \qquad(15)
```

For the calibrated word, the completed odd block can be fused explicitly.
Put $`\omega=e^{i\pi/4}`$. After Hadamards on its odd targets, an edge of $`O_xO`$ has phase
$`\omega^{2x_i b-4x_iab}`$. Thus the dirty cubic compiler uses r
quadratic occurrences, r cubic occurrences and no linear term. Its
middle-block bounds are

```math
T\le20r,\qquad D_T\le12,\qquad w\le3r,\qquad
G_{\rm all}\le m+40r.                                \qquad(16)
```

The m counts the two Hadamard walls. Combining (16), the two endpoint
pulses and the actual C pair gives

```math
T(W)\le2m+84r,\quad D_T(W)\le38,\quad
w(W)\le5r,\quad G_{\rm all}(W)\le13m+194r.             \qquad(17)
```

All bounds have zero clean inputs and zero synthesis error. Each pulse
returns its arbitrary helpers, so the same pool is reused. All phases,
actual inverses and arbitrary external references are included. A
whole-source implementation error epsilon weakens (8) to
$`\sqrt3-2\varepsilon`$ and (14) to $`2-2\varepsilon`$.

The paid long Clifford frame supplies substantial coefficient transport,
but its first three-pulse word fails calibration and its calibrated repair
has the wrong complete address sector. These are two specified-word
failures. They do not exclude arbitrary three- or four-pulse sources,
other Clifford networks, different nonlinear bases or address motion.

The [exact fixtures](../../tests/test_clifford_interleaved_lookup.py) check
phase thresholds, complete commutators, address/helper normalization, native
source words, actual inverses and invariant-sector witnesses. The all-size
claims follow from the proofs above, not from finite fixtures. The next source
must change the failed joint transition or implement the unpaid
[coherent central transport](COHERENT_CENTRAL_TRANSPORT.md).
