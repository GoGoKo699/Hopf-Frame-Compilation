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

## 8. Exactly which accepted blocks three-call recovery returns

The half-unitary scalar block is sufficient but not necessary. Let Q now
be any unitary, set $`K=J^\dagger QJ`$, and let W be the desired unitary
on the complete input space. With $`Z=I-2JJ^\dagger`$ and
$`F=-QZQ^\dagger ZQ`$, the exact characterization is

```math
FJ=JW
\quad\Longleftrightarrow\quad
K=W\left(\frac32P-I\right)
\quad\hbox{for some orthogonal projector }P.
```

No commutation between P and W is required. Indeed the accepted block is
$`3K-4KK^\dagger K`$. If it equals W, then
$`H=W^\dagger K=3K^\dagger K-4(K^\dagger K)^2`$ is Hermitian.
Substitution of K=WH gives $`3H-4H^3=I`$. Since K is a contraction,
the spectrum of H lies in [-1,1]; factoring the scalar equation gives

```math
1-3h+4h^3=(h+1)(2h-1)^2,
\qquad
H=\frac32P-I.
```

Conversely this H gives accepted block W. Since FJ and JW are
isometries, an accepted unitary block leaves no rejected component,
proving the complete returned-work equality.

There is also an exact stability identity for every Hermitian contraction
H when K=WH. Writing $`E=FJ-JW`$ gives

```math
E^\dagger E
=2I-2(3H-4H^3)
=8(H+I)(H-I/2)^2,
\qquad
\|E\|^2=\left\|8(H+I)(H-I/2)^2\right\|.
```

This includes all leakage and dirty return. The two allowed amplitudes
have different stability: perturbing h=1/2 gives error proportional to
the perturbation, whereas $`h=-1+\epsilon`$ gives squared error
$`8\epsilon(3/2-\epsilon)^2`$. The projector-amplitude option therefore
does not automatically inherit the linear error estimate of Section 3.
Producing either admissible accepted block with the required native
resources remains a separate construction problem.

## 9. A balanced scalar product and its paid native cost

A concrete alternative joins scalar words without requiring a
Majorana-preserving loader. Suppose the unitary U has one distinguished signal a and
otherwise arbitrary inputs, with
$`(U+U^\dagger)/2=cI`$ and $`Z_aUZ_a=U^\dagger`$.
On two identical copies let $`C=\mathrm{CZ}_{a,b}`$ and let
$`\Sigma`$ swap the complete copies. Define

```math
R=U_A C U_B C,\qquad
\mathcal G=Z_aZ_bC\Sigma.
```

Writing U=cI+A, its skew part anticommutes with the signal Z.
C changes the second skew part to $`Z_aA_B`$, which anticommutes
with the first. Hence $`(R+R^\dagger)/2=c^2I`$ on every input.
$`\mathcal G`$ is a Hermitian Clifford involution: it exchanges the two twisted
child words while reversing each, so $`\mathcal G R\mathcal G=R^\dagger`$.
The scalar-product identity permits unequal c values, but this physical
swap grading uses identical children.

The grading can be converted to a Pauli with a paid native word.
For children of width n, the signal factor is
$`g=Z_aZ_b\mathrm{CZ}_{a,b}\mathrm{SWAP}_{a,b}`$.
The chronological word

```math
\mathrm{CX}_{a,b},\quad
\mathrm{CH}_{b,a},\quad
\mathrm{CX}_{b,a}
```

conjugates g to $`Z_a`$. A literal controlled H uses two T gates:
$`\mathrm{CH}=V\mathrm{CZ}V^\dagger`$ with
$`V=SHTHS^\dagger`$ and $`VZV^\dagger=H`$.
Inverse Bell transforms on the n-1 core pairs change their SWAPs to CZs.
The resulting grading is $`Z_a\prod_j\mathrm{CZ}_{x_j,y_j}`$.
Toffolis with these control pairs and common target a remove that product.

For m=n-1 pairs, the combined CCZ phase exponent, in units of pi/4, is

```math
ma+\sum_j\left[
x_j+y_j-(a\mathbin\oplus x_j)-(a\mathbin\oplus y_j)
-(x_j\mathbin\oplus y_j)+(a\mathbin\oplus x_j\mathbin\oplus y_j)
\right].
```

Aggregate the common target phase modulo eight and compute each of the
four parity families on existing core wires, applying its phase and
uncomputing. Together with H on the common target before and after,
this gives the Toffoli product with $`6m+(m\bmod2)`$ T gates,
at most five T layers and O(m) elementary Cliffords. No helper is
initialized. Thus a normalizer N with $`N\mathcal G N^\dagger=Z_a`$ costs
two additional T gates and at most seven T layers when n is at least two.
For n=1 only the two-layer controlled H is needed.

Set $`U_{\rm new}=NRN^\dagger`$ and use actual inverses throughout.
Starting from the literal Clifford $`U_1=\exp(i\pi X/4)`$, including
its Clifford scalar, iteration at dyadic width w gives

```math
\frac{U_w+U_w^\dagger}{2}=2^{-w/2}I,\qquad
Z_aU_wZ_a=U_w^\dagger,\qquad
J^\dagger U_wJ=2^{-w/2}I.
```

This is an exact one-clean scalar block on w-1 arbitrary dirty inputs.
It is not a returned-work rotation: its rejected norm is
$`\sqrt{1-2^{-w}}`$. It also supplies only this sequence of coefficients,
not an arbitrary-angle program. For dyadic w at least two, the emitted
word has

```math
T=6w\log_2w-7w+6,\qquad
G=O(w\log w),\qquad
D_T\le9w-14.
```

The T recurrence is T(2)=4 and $`T(w)=2T(w/2)+6w-6`$ for w at
least four; the depth recurrence pays both child words and both
normalizers. Native commutation checks of this particular word, including
equal-axis cancellation and Clifford absorption, still found growing
dependencies at bounded widths. They do not prove an all-size or
unrestricted depth lower bound. Thus this paid attempt currently has
linear depth and superlinear count, and does not improve R96.

The [native fixtures](../../tests/test_joint_precision_rotation.py)
check literal small joins, their actual inverses and ledgers, and the
projector-amplitude recovery with noncommuting W and P. The remaining
precision task is to produce a scalar or projector-amplitude accepted
block, or directly its complete recovery, with linear count and
sublinear precision depth while paying every initialization and return.

<a id="two-path-recovery-bypass"></a>

## 10. Balanced two-path recovery has a direct native bypass

The projector-amplitude option suggests a concrete pilot: interfere two
native unitaries U and V on the complete target/dirty space, using one
Hadamard-prepared signal, to obtain $`K=(U+V)/2`$. This pilot has an
exact simplification. Its recovery can be replaced by either of two
ordinary native words, with at most a square-root-two error factor and
no initialized signal. This conclusion also holds for any other unitary
completion Q with that same accepted block.

Set $`A=-UV^\dagger U`$ and $`B=-VU^\dagger V`$, retaining both
literal minus signs. For an arbitrary complete unitary target W and
$`E=FJ-JW`$, the full-error identity is

```math
E^\dagger E
=\frac12\left[(A-W)^\dagger(A-W)+(B-W)^\dagger(B-W)\right].
```

In particular, if the complete recovered error is at most epsilon, then

```math
\|A-W\|,\ \|B-W\|\le\sqrt2\,\epsilon.
```

These are full operator bounds, including every dirty input and reference;
neither branch is assumed to return work before the argument. To verify
the identity, put $`R=U^\dagger V`$. Direct multiplication gives

```math
3K-4KK^\dagger K
=-\frac12U(R^\dagger+R^2)=\frac{A+B}{2}.
```

For every unitary completion, isometry of FJ and JW gives

```math
E^\dagger E=2I-W^\dagger(3K-4KK^\dagger K)
-(3K-4KK^\dagger K)^\dagger W,
```

which proves the claim.
For the canonical source
$`Q=H_a\,\mathrm{diag}(U,V)\,H_a`$, there is also the literal identity

```math
F=H_a\,\mathrm{diag}(A,B)\,H_a.
```

Thus its rejected block is (A-B)/2. Discarding that block would lose the
condition which makes both direct words accurate. The factor square root
of two is sharp: take U=I, V=T and W=-T-dagger on one qubit. Then A=W,
whereas B=-T-squared differs from W, and the error identity is exactly
half the squared B error. All phases in this witness are literal.

Exact recovery occurs precisely when A=B=W. Equivalently,

```math
R^3=I,\qquad W=-UR^\dagger,
\qquad P=I-\frac{I+R+R^2}{3},
\qquad K=W\left(\frac32P-I\right).
```

Here P is the spectral projector onto the two nontrivial cube-root
eigenspaces of R. This provides genuine noncommuting W/P examples, but
does not synthesize W for free: the direct word already implements it.
For example, let R be the order-three two-CNOT permutation on target t
and dirty wire d, choose the literal $`W=H_tT_tH_t`$, and put
$`U=-WR`$, $`V=-WR^2`$. The associated P does not commute with W.
The canonical Q has an explicit native implementation: apply R,
signal-controlled R, and -W between the signal Hadamards. The two
controlled CNOTs are two seven-T Toffolis. This implementation has
15 T gates, and its three-call recovery has 45; either bypass word has
three before further simplification. This is a bounded demonstration of
the reduction, not a new arbitrary-angle family.

More generally, let $`(t_U,d_U,g_U)`$ and $`(t_V,d_V,g_V)`$ bound the
T count, T-depth and elementary Clifford count of the literal branches
and their actual inverses. The two alternatives have simultaneous bounds

```math
\begin{array}{c|ccc}
 &T&D_T&G\\ \hline
A&2t_U+t_V&2d_U+d_V&2g_U+g_V+O(1)\\
B&t_U+2t_V&d_U+2d_V&g_U+2g_V+O(1).
\end{array}
```

Choose one complete row according to the desired resource tradeoff; minima
from different rows are not a simultaneous ledger. The branches use the
same physical target/dirty space, their actual inverses reuse it, and the
signal is absent. The literal -I costs zero T gates and at most twelve
H/S Cliffords via $`XZXZ=-I`$, with $`X=HS^2H`$ and $`Z=S^2`$.
Requesting complete recovered error $`2^{-L}/\sqrt2`$ is sufficient for direct error
$`2^{-L}`$; one additional integer accuracy bit suffices.

This reduces the pilot built from separately paid branches to direct
alternating-word synthesis. A cheap joint Q does not by itself supply
cheap independent branch words; the resource comparison assumes those
words are available. The result does not exclude such a joint source,
a joint synthesis of the alternating word, unequal or additional paths,
other signal phases, or accepted blocks without this two-unitary form.
R96 and the isolated general-angle bound above are unchanged.
The [three exact fixtures](../../tests/test_two_path_precision_bypass.py)
check all native columns and actual inverses, a noncommuting W/P source,
the completion-independent error identity, and its sharp leakage witness.

<a id="actual-clifford-signal-precision"></a>

## 11. A Clifford image of the actual signal leaves no recovery gain

The [four-band Clifford-involution compiler](CLIFFORD_INVOLUTION_SOURCE.md)
does not by itself supply a precision source. This follows by applying the
existing [one-clean Clifford-compression lemma](CROSSED_CUBIC_KERNEL_SOURCE.md#clifford-signal-gram-gap)
to the actual accepted signal, rather than to an unrelated transported probe.
Let Q be any source with one clean signal a, put

```math
A=Q^\dagger Z_aQ,\qquad M=J^\dagger AJ,
\qquad K=J^\dagger QJ,\qquad G=K^\dagger K=(I+M)/2,
```

and assume A is Clifford. It is also Hermitian, so the compression lemma
makes every eigenvalue of M belong to $`\{0,\pm1/\sqrt2,\pm1\}`$.
Consequently

```math
\mathrm{spec}(G)\subseteq
\left\{0,\frac12,1,\frac{2-\sqrt2}{4},\frac{2+\sqrt2}{4}\right\}.
```

For the three-call word F, its accepted block is
$`B=K(3I-4G)`$. A Gram eigenvalue g gives accepted singular value
$`\sqrt g\,|3-4g|`$. Among the displayed possibilities, g=1 gives one;
every other g gives at most $`\cos(\pi/8)`$. If G differs from identity,
choose a unit right singular vector v for such a smaller value. For any
complete unitary target W, the full isometry identity gives

```math
\|(FJ-JW)v\|^2
=2-2\mathrm{Re}\langle Wv,Bv\rangle
\ge2-2\cos(\pi/8).
```

Thus complete recovered error strictly below $`2\sin(\pi/16)`$ forces
G=I. In that case $`QJ=JK`$ already returns the signal exactly and
$`FJ=-JK`$: the direct source with a paid literal minus has exactly the
same complete error to W, including all dirty-return error. In particular,
exact recovery $`K=W(3P/2-I)`$ forces P=0, since its nonzero P sector
would require the unavailable Gram eigenvalue one quarter.

For the compiler's source grammar $`Q=CDB`$, with Clifford C,B and diagonal
third-level D, the actual image A is Clifford for every signal choice:
negating D's weighted cubic phase keeps $`D^\dagger`$ diagonal third-level,
and the Clifford boundaries preserve that level in
$`Q^\dagger=B^\dagger D^\dagger C^\dagger`$.
Arbitrary Clifford boundary mixing therefore does not change this conclusion.
Knowing only that some other specified Pauli has a Clifford image would not
establish the hypothesis on A. A genuinely composed source whose actual
signal image is non-Clifford remains eligible. This is an application of
the retained compression lemma, not a new general precision-depth bound.

<a id="biased-two-path-precision"></a>

## 12. A paid unequal-path pilot also reduces to a direct branch

An unequal signal seed leaves the balanced two-path grammar in Section 10.
The literal chronological one-qubit word H,T,H,T,H has two T gates and
zero-state probability three quarters. Call its unitary S and use

```math
Q=S_a^\dagger\,\mathrm{diag}(U,V)\,S_a,
\qquad K=\frac{3U+V}{4}=\frac14U(3I+R),
\qquad R=U^\dagger V.
```

Here U and V are separately paid native words on the entire target/dirty
space. Their controlled implementation is an additional cost of Q; it is
not granted for free by this formula. No precision-dependent signal state
is supplied. Put $`T_R=(R+R^\dagger)/2`$. Direct multiplication gives

```math
\begin{aligned}
G&=\frac{5I+3T_R}{8},\\
B=3K-4KK^\dagger K
 &=\frac1{16}U(3I-7R-9R^\dagger-3R^2),\\
\Lambda=I-B^\dagger B
 &=\frac{27}{32}(I-T_R)(I+T_R)^2.
\end{aligned}
```

The positive operator Lambda is the complete rejected Gram of FJ, for
every unitary completion with this K. Exact recovery requires Lambda=0.
Since R is unitary, this is equivalent to $`R^2=I`$, and then
$`B=-UR=-V`$. Thus exact recovery implements an already supplied branch.
The projector-amplitude form is retained, with $`P=(I-R)/2`$ and
$`W=-V`$; no commutation of W with P is required.

There is also a uniform approximate reduction. Write D=B+V. Functional
calculus of the unitary R gives

```math
D^\dagger D
=\frac9{32}(5I-3T_R)(I-T_R)(I+T_R),
\qquad (D^\dagger D)^2\preceq12\Lambda.
```

For the inequality, on $`-1\le t\le1`$ the ratio of the two sides before
the constant 12 is $`(3/32)(5-3t)^2(1-t)\le12`$; endpoints follow by
continuity. If the complete recovered error to W is at most epsilon,
then $`\|B-W\|\le\epsilon`$ and $`\|\Lambda\|\le\epsilon^2`$.
Therefore

```math
\|-V-W\|\le\epsilon+12^{1/4}\sqrt\epsilon
\le3\sqrt\epsilon\qquad(0\le\epsilon\le1).
```

Requesting complete recovered error at most $`2^{-2L}/9`$ suffices for
direct error $`2^{-L}`$; $`2L+4`$ integer accuracy bits suffice. The direct
word -V uses V's T count, T-depth and dirty allocation, with at most twelve
additional elementary Cliffords for the literal minus. Thus a uniform
family whose separately paid branches have linear precision count/work and
sublinear precision depth already provides a direct branch with those same
asymptotic bounds after the constant-factor accuracy rescaling. This
comparison requires the separate branch words; a cheaply synthesized
joint completion need not supply them or their resource bounds.

The pilot does genuinely leave the Clifford-signal hypothesis. Take U=I
and let V=R be the two-CNOT order-three permutation on two arbitrary wires.
Its controlled version is two seven-T Toffolis, so the source above has
18 T gates, at most fourteen T bands and 22 elementary Cliffords; its
literal three-call recovery has 54 T gates, at most 42 bands and 82
Cliffords. These are common actual-inverse bounds with one clean signal,
two arbitrary data/work wires and no additional helpers. Let P project onto
the two nontrivial cube-root eigenspaces of R. Then

```math
G=I-\frac9{16}P,\qquad M=2G-I=I-\frac98P.
```

The eigenvalue minus one eighth in M is excluded by Section 11, proving
that the actual signal image is non-Clifford. Nevertheless the accepted
recovery singular value on P is $`5\sqrt7/16`$. Every complete unitary
target consequently has recovered error at least
$`\sqrt{2-5\sqrt7/8}`$. This exact native trial changes the signal image
but does not give a precision-scaling construction.

The [signal and biased-path fixtures](../../tests/test_precision_signal_and_bias.py)
retain all native columns, actual inverses, rejected Grams and the exact
non-Clifford-signal witness. The single unequal-path pilot stops at its
direct-branch reduction; other weights, additional paths and jointly
compiled sources remain open. R96 and the general precision frontier
are unchanged.
