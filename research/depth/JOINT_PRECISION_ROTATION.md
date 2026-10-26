# A jointly compiled one-clean precision rotation

[One-clean source](../../docs/ONE_CLEAN_COMPILER.md) ·
[Paired native tails](SOURCE_T_DEPTH.md#3-two-tails-in-parallel) ·
[General depth frontier](../../docs/OPEN_PROBLEM.md)

An isolated prescribed rotation admits a smaller version of the existing
paired-source circuit. A finite Clifford chart raises its accepted
normalization to one half, so three source calls suffice for recovery.
A common dirty conjugation can then be omitted from the entire recovered
word. The only growing operations left are programmed masked-tail
conjugations inside a fixed signal/target circuit.

This is a local component and an explicit simplification of that source
architecture. Its depth remains linear in precision. It does not improve
R96, supply a coherently addressed chart selector, or assert an optimal
isolated-rotation width: ordinary single-qubit approximation already gives
a different construction without a dirty bank.

## 1. Complete local statement

Use $`R_y(\theta)=e^{-i\theta Y}`$. Let J append one zero flag a,
leaving the target t and every dirty input arbitrary. For every integer
$`q\ge5`$ and effectively prescribed sine/cosine pair, the construction
below emits a literal native Clifford+T word F with

```math
\|FJ-J(R_y(\theta)\otimes I_{q+1})\|_{\rm op}
\le\frac{5\sqrt{481}}{16}\,2^{-q}\lt7\,2^{-q}.
```

It uses one clean qubit, q+1 arbitrary dirty qubits, no additional helpers,
and the simultaneous conservative bounds

```math
T\le24q+84,\qquad D_T\le12q+96,\qquad G=O(q).
```

Here G counts elementary Cliffords, including actual inverse phase words.
For integer $`L\ge6`$, take $`q=L+3`$: the complete error is below
$`2^{-L}`$ with L+4 dirty qubits. The norm includes clean leakage,
dirty return, all logical columns, literal phases and arbitrary references.
There is no supplied precision state, measurement or reset.

## 2. A finite Clifford chart gives half normalization

The literal Clifford $`HZ=R_y(\pi/4)`$. Choose
$`k\in\{0,\ldots,7\}`$ so that the residual angle
$`\alpha=\theta-k\pi/4`$, modulo $`2\pi`$, lies in
$`(-\pi/6,\pi/6)`$. The nearest center has distance at most
$`\pi/8`$, so this larger open chart has a fixed overlap margin.
Certified interval evaluation can test
$`\cos(\theta-k\pi/4)\gt\sqrt3/2`$ for each k; at least one
strict test eventually succeeds. No exact boundary or sign decision is
required. At the end apply $`C_k=(HZ)^k`$, including its literal sign
when k is at least four.

The paired source on m=q+1 dirty qubits uses the Hermitian Majoranas

```math
\Gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\Gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j,
\qquad\{\Gamma_i,\Gamma_j\}=2\delta_{ij}I.
```

Its squared weights are one quarter on the head, one half on the odd
tail, and one quarter on the even tail. Each tail has normalized weights
$`w_j=2^{-j-1}`$ for $`0\le j\lt q-1`$ and
$`w_{q-1}=2^{1-q}`$. Let M denote this source, and let
$`N_f=P_fMP_f^\dagger`$ flip exactly the even tail. Then
$`\{M,N_f\}/2=cI`$ with $`c=1/2`$.

For a programmed Pauli mask $`P_g`$, put
$`N_g=P_gMP_g^\dagger`$. Write its head sign as +1 and its signed
tail means as u,v. The scalar moments are

```math
p_+=\frac14+\frac u2,\qquad p_-=\frac v4,
\qquad s=p_++p_-,\qquad\tau=p_+-p_-.
```

Define the actual scalar word, with rightmost factors acting first, by

```math
\mathcal S(M,N)=H_a C_{a=0}(M)N C_{a=1}(M)H_a.
```

If $`\{M,N\}/2=sI`$, this equals
$`sI+X_a(MN-NM)/2`$. With the literal routing Cliffords, set

```math
A=\mathrm{CZ}_{a,t}\mathcal S(M,N_f)\mathrm{CZ}_{a,t},
\qquad
B_g=\mathrm{CX}_{a,t}\mathcal S(M,N_g)\mathrm{CX}_{a,t},
\qquad Q=A^\dagger B_gA.
```

The [full scalar-word identity](../../docs/ONE_CLEAN_COMPILER.md#3-conjugating-scalar-blocks-produces-a-rotation)
gives, on every dirty input,

```math
K=J^\dagger QJ=sI+tXZ_t,\qquad t=s/2-\tau.
```

The rejected operators are retained in Q. In particular, this identity
does not multiply accepted blocks in isolation.

To target $`K_*=R_y(\alpha)/2`$, choose

```math
s_*=\frac12\cos\alpha,\qquad
\tau_*=\frac{s_*}{2}-\frac12\sin\alpha,
```

```math
u_*=\frac34\cos\alpha-\frac12\sin\alpha-\frac12,
\qquad
v_*=\frac12\cos\alpha+\sin\alpha.
```

On the closed larger chart $`|\alpha|\le\pi/6`$,

```math
\frac{3\sqrt3-6}{8}\le u_*\le\frac{3\sqrt3-2}{8},
\qquad
\frac{\sqrt3-2}{4}\le v_*\le\frac{\sqrt3+2}{4}.
```

Both intervals lie strictly inside [-1,1]. Apply the existing certified
signed-dyadic rule to the two means. Its errors in the two coefficient
sums obey

```math
|a|=|\Delta p_+|\le\frac54\,2^{-q},\qquad
|b|=|\Delta p_-|\le\frac58\,2^{-q}.
```

Consequently the actual scalar block has error

```math
\begin{aligned}
\delta^2=\|K-K_*\|^2
&=(a+b)^2+(-a/2+3b/2)^2\\
&=\frac54a^2+\frac12ab+\frac{13}{4}b^2
\le\frac{925}{256}\,2^{-2q}.
\end{aligned}
```

The maximum of this last quadratic on the indicated rectangle is at
equal-sign extreme corners. Thus
$`\delta\le(5\sqrt{37}/16)2^{-q}`$, which is below 1/8 for q at
least five. The full source and masks are exact native circuits; only
these two coefficient means are rounded.

## 3. Three-call recovery and its full error

Let $`Z=Z_a`$. Since $`I-2JJ^\dagger=-Z_a`$, normalization-two
recovery is the literal word

```math
F_0=-QZQ^\dagger ZQ.
```

The leading minus is necessary and is a Clifford scalar, for example
$`XZXZ=-I`$. The middle call is the actual inverse of Q.
The [normalization-two identity](../../docs/FAULT_TOLERANT_COMPILER.md#4-from-a-half-unitary-block-to-the-complete-output)
gives

```math
J^\dagger F_0J=3K-4KK^\dagger K=(3-4r^2)K,
\qquad r^2=s^2+t^2.
```

Here the real two-coefficient block permits a direct full-output estimate.
Put $`W=R_y(\alpha)`$ and
$`h=s\cos\alpha+t\sin\alpha`$. Then
$`\delta^2=r^2+1/4-h`$ and $`W^\dagger K+K^\dagger W=2hI`$.
Both $`F_0J`$ and JW are isometries, so

```math
\begin{aligned}
\|F_0J-JW\|^2
&=2-2(3-4r^2)h\\
&=8(r^2-1/4)^2+2(3-4r^2)\delta^2.
\end{aligned}
```

This equality includes all rejected output. Since
$`|r-1/2|\le\delta\le1/8`$, the last expression is at most

```math
\bigl[8(r+1/2)^2+6-8r^2\bigr]\delta^2
=(8r+8)\delta^2\le13\delta^2.
```

Combining this with the coefficient bound proves
$`\|F_0J-JW\|\le(5\sqrt{481}/16)2^{-q}\lt7\,2^{-q}`$.
The final target-only Clifford C_k restores the original angle without
changing this error.

## 4. Remove a common dirty conjugation from the whole word

Write the paired loader as a fixed two-rotation seed followed by its
two tails. With

```math
R_{ij}=\exp\!\left(\frac{i\pi}{8}(i\Gamma_i\Gamma_j)\right),
\qquad U_{\rm seed}=R_{02}R_{01},
```

the fixed source is

```math
M_0=U_{\rm seed}\Gamma_0U_{\rm seed}^\dagger
=\frac12\Gamma_0+\frac1{\sqrt2}\Gamma_1+\frac12\Gamma_2.
```

Let V contain the q-1 paired tail rounds, chronologically applying
$`R_{2j+1,2j+3}`$ and $`R_{2j+2,2j+4}`$ for
$`j=0,\ldots,q-2`$. Then $`M=VM_0V^\dagger`$.
Each tail plane joins two Majoranas with the same fixed f sign, so
$`[P_f,V]=0`$. Therefore

```math
N_f=VN_{f0}V^\dagger,\qquad
N_{f0}=\frac12\Gamma_0+\frac1{\sqrt2}\Gamma_1-\frac12\Gamma_2.
```

Both seed sources act only on the first two dirty qubits; explicitly
$`N_{f0}=Z_1M_0Z_1`$. Define

```math
\widetilde P_g=V^\dagger P_gV,\qquad
\overline N_g=\widetilde P_gM_0\widetilde P_g^\dagger,
```

and obtain $`A_0,\overline B_g`$ from the displayed scalar/routing
words by using $`(M_0,N_{f0})`$ and $`(M_0,\overline N_g)`$.
All operations outside $`\widetilde P_g`$ and its actual inverse have
bounded support. Extending V by identity on the flag and target gives

```math
Q=V\overline QV^\dagger,\qquad
\overline Q=A_0^\dagger\overline B_gA_0.
```

Because V acts only on arbitrary dirty inputs, it intertwines J with
$`I_t\otimes V`$, which commutes with the desired logical target. For
any physical word F and logical W,

```math
\|FJ-J(W\otimes I)\|
=\|(V^\dagger FV)J-J(W\otimes I)\|.
```

This is unitary invariance of the norm; exact return by F is not assumed.
Since V also commutes with $`Z_a`$, the emitted recovery may be

```math
\overline F_0=-\overline QZ_a\overline Q^\dagger Z_a\overline Q,
```

with both exterior dirty conjugations omitted. The target-only correction
C_k also commutes with V. Thus the same complete error holds for the
smaller emitted word. This omission does not supply a cheap arbitrary
change of basis for the three source vectors; the programmed masked-tail
operation remains in the circuit.

## 5. Native ledger and the remaining precision problem

The [paired-tail schedule](SOURCE_T_DEPTH.md#3-two-tails-in-parallel)
implements V and its actual inverse with common bounds

```math
T(V)=2(q-1),\qquad D_T(V)=q-1,\qquad G(V)=O(q),
```

without an additional helper. Each seed source, including either
flag-controlled version, uses $`U_{\rm seed}`$, a central X or
controlled X, and the actual seed inverse: at most four T gates and
four T layers, with bounded Clifford count. Native representatives of
the two seed rotations can carry common scalars; those cancel in every
actual source conjugation. A fixed Pauli gives N_f0 with the same bound.

A_0 contains three bounded seed-source occurrences. The programmable
$`\overline B_g`$ contains three and four V/V-adjoint occurrences,
from its two masked-tail words. Hence $`\overline Q`$ contains nine
bounded sources and four tails. Three-call recovery contains 27 bounded
sources and twelve tails, giving

```math
T\le108+24(q-1),\qquad
D_T\le108+12(q-1),\qquad G=O(q).
```

The Pauli programs, routing Cliffords, actual inverses, literal minus,
and C_k are included. No query or online chart selector is used: the
isolated angle and its chart are classically prescribed. All q+1 core
wires remain arbitrary dirty inputs; only a starts initialized.

The remaining growing word is $`V^\dagger P_gV`$ inside this fixed
signal/target interaction. A joint resynthesis may leave the
Majorana-preserving class. Independently balancing the existing tail
loader does not: its restricted denominator bound remains applicable.
A sublinear-depth complete masked-tail/recovery word, or a different
source with the same full scalar block, is still required.

General coherent addresses and active guards also need their own paid
interfaces. In particular, three-call recovery with a global leading
minus would act as -I where Q is inactive identity; an active-only sign
must be supplied there. None of these additional interfaces follows from
the isolated construction, and R96's general depth bound is unchanged.

## 6. A third-level reflection cannot supply the calibration

A simpler hierarchy-constrained replacement for M has an explicit
limitation. Suppose M is a Hermitian involution in the third Clifford
hierarchy level and P is a Hermitian Pauli. Put $`N=PMP`$,
$`C=MPM`$ and $`K=MN=CP`$. Then C is a Hermitian Clifford, so

```math
K^2=CPCP=P'P,\qquad P'=CPC
```

is a product of Hermitian Paulis. Its fourth power is identity, hence
$`K^8=I`$. Consequently the spectrum of
$`(MN+NM)/2=(K+K^\dagger)/2`$ is contained in

```math
\mathcal E=\{0,\pm1/\sqrt2,\pm1\}.
```

For any real c, spectral calculus gives

```math
\left\|\frac{MN+NM}{2}-cI\right\|
\ge\mathrm{dist}(c,\mathcal E).
```

In particular the required fixed calibration $`c=1/2`$ has error at
least $`1/\sqrt2-1/2`$, regardless of physical width. This is a
restriction on this reflection/Pauli-mask interface, not a lower bound
on the error of an arbitrary complete rotation circuit.

The calibration gap is sharp. On two qubits let

```math
M=X_u\mathrm{CX}_{u,v}(T_u^\dagger T_v)\mathrm{CX}_{u,v},
\qquad P=X_v.
```

The diagonal part has phase $`\omega^{v-2uv}`$,
$`\omega=e^{i\pi/4}`$, so M is a Hermitian third-level involution.
The relative word K has phase
$`\omega^{(1-2u)(1-2v)}`$ and therefore real part $`I/\sqrt2`$.
This uses two T gates in one native layer and no helper.

More generally, if $`M\in\mathcal C_k`$ is Hermitian, repeated
squaring of MP descends the hierarchy: inductively
$`(MP)^{2^j}=C_jP`$, where C_j is Hermitian and belongs to
$`\mathcal C_{k-j}`$. At j=k-1 the right side is a phased Pauli.
Thus $`(MN)^{2^k}=I`$. Exact calibration one half would instead force
the eigenvalues of MN to be $`e^{\pm i\pi/3}`$, which have order six.
The fixed seed M_0 already realizes this calibration with $`P=Z_1`$.
It therefore lies outside every finite hierarchy level despite its
constant native implementation. Hierarchy level cannot replace native
T-depth in this precision problem.

## 7. Universal Pauli-mask closure fixes the source algebra

There is also an exact converse behind the scalar-overlap representation.
Let M be any Hermitian involution on m qubits and expand it in the
Hermitian Pauli basis:

```math
M=\sum_u a_uP_u,\qquad a_u\in\mathbb R,\qquad\sum_u a_u^2=1.
```

Suppose $`\{M,P_vMP_v\}/2`$ is scalar for every Pauli label v.
Then the distinct Pauli terms with nonzero coefficients anticommute
pairwise. Conversely, that pairwise anticommutation implies the
scalar-overlap property for every Pauli mask.

To prove the forward direction, label Paulis by
$`\mathbb F_2^{2m}`$ and write
$`\chi_u(v)=(-1)^{[u,v]}`$ for their symplectic characters. Put
$`M_v=P_vMP_v`$. Conjugating the assumed scalar identity by another
Pauli shows

```math
\frac{\{M_x,M_y\}}2=c_{x+y}I
\quad\hbox{for every }x,y.
```

Take the normalized Fourier coefficient with characters
$`\chi_u(x)\chi_v(y)`$. The left side is
$`a_ua_v\{P_u,P_v\}/2`$. The right side vanishes when u and v
are distinct, because its scalar depends only on x+y. Thus every
distinct pair of nonzero terms anticommutes. Conversely, those cross
terms cancel directly, leaving

```math
\frac{\{M,M_v\}}2=\sum_u a_u^2\chi_u(v)I.
```

If the identity coefficient is nonzero, M is simply $`\pm I`$.
Otherwise there are at most 2m+1 nonzero Pauli terms: the symplectic
Gram matrix of r pairwise anticommuting labels has zero diagonal and
one off the diagonal. Its binary rank is r for even r and r-1 for
odd r, and cannot exceed 2m.

Thus retaining exact scalar overlap for every physical Pauli mask
already forces an anticommuting Pauli-vector source. This does not
force a Majorana-preserving native synthesis or a growing depth bound.
A finite chosen program family, approximate scalar identities, or a
complete accepted block without universal mask closure can avoid the
premise. The current joint word remains a candidate for unrestricted
resynthesis even though its original source has this closure.

The [bounded native checks](../../tests/test_joint_precision_rotation.py)
verify the chart margins, coefficient constants, actual inverses, source
and recovery conjugations, literal phases, and native ledgers. Small
matrices also distinguish the scalar seed from a non-scalar Clifford
mask example. These fixtures support the identities; the all-precision
error and the source-class statements use the analytic arguments above.
