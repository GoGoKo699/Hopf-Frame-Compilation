# Joint Fourier words and native source limits

[Gaussian scalar-loop limit](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md) ·
[Native sharing requirement](NATIVE_INDICATOR_SHARING.md) ·
[Paid removal of a bounded clean bank](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)

Joint noncommuting words remain a live route to the complete approximate
dirty indicator. The tests below concern particular words that use
nonquadratic phases, so they lie outside the earlier Gaussian-round
argument. They identify where these candidates stop: a repeated
Fourier/Boolean-phase loop has a robust discrete closure condition at
bounded work degree, and a single Boolean-phase block between arbitrary
Cliffords has the wrong singular values for normalization-two recovery.

None of these screens is an unrestricted T-depth bound. Different phase rounds,
address motion between rounds, and more general noncommuting native
words remain outside their hypotheses. No complete-frame resource
frontier changes here; R96 remains the strongest unconditional general
upper theorem.

## 1. Exact fourfold closure retains only a simple sign

Let f be any Boolean function on n work bits, with n at least one, and
put

```math
F(w)=(-1)^{f(w)},\qquad D_f=\mathrm{diag}F,\qquad
H=H^{\otimes n},\qquad C=HD_fH.
```

Use the normalized Walsh coefficients
$`\widehat F(t)=2^{-n}\sum_w(-1)^{t\cdot w}F(w)`$.
Then $`C_{w+t,w}=\widehat F(t)`$. For $`h\in\{0,1\}`$,
the literal scalar-return condition is equivalent to

```math
(HD_f)^4=(-1)^hI
\quad\Longleftrightarrow\quad
CD_f=(-1)^hD_fC.
\qquad(1)
```

Indeed, $`(HD_f)^2=CD_f`$ and $`(CD_f)^{-1}=D_fC`$.
Comparing the indicated matrix entries proves that (1) holds exactly
when every t in the Fourier support satisfies

```math
f(w+t)+f(w)=h\quad\hbox{for every }w.
\qquad(2)
```

Fourier transforming this translation identity gives
$`s\cdot t=h`$ for every pair s,t in the same Fourier support.
Set s equal to t. Every supported frequency has parity h, and hence
$`F(w+\mathbf1)=(-1)^hF(w)`$. In particular,

```math
h=f(\mathbf1)+f(0).
\qquad(3)
```

This does not force f to be quadratic. For example,
$`f(w)=(w_1+w_2)(w_3+w_4)(w_5+w_6)`$ has a cubic phase and
Fourier support in the span of the three disjoint pair masks. That span
is self-orthogonal, so (2) holds with h zero. What is restricted is the
scalar extracted by this repeated word, not the possible internal phase.

## 2. Bounded work degree makes the closure test robust

Suppose $`\deg f\le d`$, where $`d\ge1`$. For either sign s,

```math
(HD_f)^4\ne sI
\quad\Longrightarrow\quad
\|(HD_f)^4-sI\|_{\rm op}\ge2^{2-d}.
\qquad(4)
```

The bound is independent of the number of dirty work bits. In
particular, a cubic phase whose fourth word is within error strictly
less than one half of a scalar must satisfy that scalar identity exactly.

To prove (4), write $`s=(-1)^h`$ and define

```math
\mu(t)=|\widehat F(t)|^2,\qquad
q_h(t)=\Pr_w[f(w+t)+f(w)\ne h],\qquad
G_h=\{t:q_h(t)=0\}.
```

Parseval gives $`\sum_t\mu(t)=1`$. Translation directions on
which f has constant derivative form a linear subspace, and their
constant derivative is a linear functional there. Thus G_h is an
affine subspace or is empty. Every nonzero q_h(t) is at least
$`a=2^{1-d}`$: its mismatch polynomial is nonzero and has degree
at most d-1. The elementary Boolean polynomial weight bound used here
follows by induction: if a degree-r polynomial depends on the last
variable, its two restrictions differ by a nonzero polynomial of degree
at most r-1, and on each point of that difference at least one
restriction is one. If it does not depend on that variable, delete it.

If Fourier mass is not entirely in a nonempty G_h, there is a Fourier
hyperplane $`v\cdot t=b`$ containing G_h but excluding some
positive-mass frequency. Its outside mass obeys

```math
\sum_{t:\,v\cdot t\ne b}\mu(t)
 =\Pr_w[f(w+v)+f(w)\ne b]\ge a.
\qquad(5)
```

The equality is the Fourier autocorrelation identity; the last
inequality uses the same nonzero derivative weight bound. Therefore
$`\mu(G_h^c)\ge a`$ whenever it is nonzero. If G_h is empty,
its outside mass is one and the same bound holds. Finally, the
normalized Hilbert--Schmidt norm gives

```math
\begin{aligned}
\|(HD_f)^4-sI\|_{\rm op}^2
 &=\|CD_f-sD_fC\|_{\rm op}^2\\
 &\ge2^{-n}\|CD_f-sD_fC\|_{\rm HS}^2\\
 &=4\sum_t\mu(t)q_h(t)\ge4a^2.
\end{aligned}
\qquad(6)
```

If the outside mass is zero instead, (2) holds and the word closes
exactly. This proves the asserted alternative.

Now allow a preserved address x and phases f_x of uniformly bounded
work degree d. A full-operator approximation to
$`(-1)^{h(x)}I`$ with error below $`2^{2-d}`$ forces (3)
at every x. If every coefficient of f_x as a work polynomial is affine
in x, h is affine in x as well. Thus replacing quadratic work phases
by cubic phases does not overcome the scalar-sign obstruction in this
particular fourfold word. The statement permits arbitrary correlated
work and references, but does not permit address motion within the word.
If d grows with the address length, the threshold in (4) shrinks; fixed
requested accuracy then need not force exact closure.

The bounded-degree hypothesis is essential. Take
$`D=I-2|0^n\rangle\langle0^n|`$ and $`c=2^{-n/2}`$.
D and HDH are reflections about vectors with overlap c. On their
two-dimensional span the product has eigenvalues
$`e^{\pm2i\arccos c}`$, and it is identity on the orthogonal
complement. Consequently

```math
\|(HD)^4-I\|_{\rm op}=4c\sqrt{1-c^2}\longrightarrow0,
```

although $`f(\mathbf1)+f(0)=1`$. This is a mathematical
counterexample to a degree-free robust extension, not a cheap native
source: the complete zero-test phase still has to be compiled.

## 3. One Boolean phase cannot supply a half-unitary block

The paid clean-removal interface allows a different constructive target:
a native source U with one signal a whose complete compression
$`A=\langle0_a|U|0_a\rangle`$ approximates half of the desired
unitary on all remaining logical and dirty registers. The signal may
have arbitrary rejected action. Here is a screen for the simplest such
source. If

```math
U=C_{\rm out}D_fC_{\rm in},
\qquad C_{\rm in},C_{\rm out}\text{ Clifford},
```

with any Boolean diagonal D_f on the entire register, then

```math
\mathrm{spec}(A^\dagger A)\subseteq\{0,1/2,1\},\qquad
\|A-V/2\|_{\rm op}\ge1/\sqrt2-1/2
\qquad(7)
```

for every unitary V on the other registers. No degree or count bound on
f is needed. The endpoint Cliffords may mix the address, output, signal,
and work arbitrarily.

Let $`Q=C_{\rm in}Z_aC_{\rm in}^\dagger`$,
$`R=C_{\rm out}^\dagger Z_aC_{\rm out}`$, and $`M=D_fRD_f`$.
Q and M are Hermitian monomial involutions with fixed binary XOR
translations. A Hermitian Pauli has matrix entries consisting of one
constant scalar in $`\{1,i,-1,-i\}`$ times real signs on its
translation; multiplication by D_f only changes those signs. Thus
$`W=QM`$ is monomial with a fixed translation and a constant scalar
in that set. Its square is diagonal with real signs, proving
$`W^4=I`$ with its literal phase retained.

It follows that the Hermitian operator
$`B=(QM+MQ)/2=(W+W^\dagger)/2`$ has spectrum contained in
$`\{-1,0,1\}`$. It commutes with Q. On the Q-positive sector,
the compression of M equals the restriction of B. Under the isometry
from the original signal-zero sector to this sector, A-dagger A is
therefore the restriction of $`(I+B)/2`$ to the Q-positive sector.
This proves the spectral claim in (7).
All singular values of V/2 are one half; the singular-value perturbation
inequality proves the stated distance, since the closest allowed
singular value is $`1/\sqrt2`$.

There is also a native one-T-layer version. If
$`U=C_{\rm out}D C_{\rm in}`$, where D is one physical layer
of T and T-dagger gates, set $`M=D^\dagger R D`$ instead.
Write the diagonal entries of D as
$`\exp(i\pi g(w)/4)`$, with g an integer linear combination
of individual computational bits. The square of $`W=QM`$ is
diagonal. Apart from a real Pauli sign, its phase exponent is

```math
g(w)-g(w+r)+g(w+q+r)-g(w+q),
```

where q,r are the XOR translations of Q,R. This exponent is even
because g modulo two is affine. Hence $`W^8=I`$, and the same
compression argument gives

```math
\begin{aligned}
\mathrm{sing}(A)&\subseteq
 \{0,\sin(\pi/8),1/\sqrt2,\cos(\pi/8),1\},\\
\|A-V/2\|_{\rm op}&\ge1/2-\sin(\pi/8).
\end{aligned}
\qquad(8)
```

Thus one physical T band with arbitrary Clifford endpoints cannot
supply the half-unitary interface at sufficiently small constant error.
This native statement imposes no read-only-address restriction. Its
constant gap alone does not exclude the allowed normalization
$`1/\sqrt2`$. Section 4 uses its flatness to resolve that one-band
case. Two or more noncommuting T bands remain outside these screens.

Consequently one cannot supply a half-unitary block merely by inserting
a high-degree or bent Boolean phase between two Clifford changes of
basis. This does not exclude other constant normalizations: the allowed
singular value $`1/\sqrt2`$ has its own paid three-call recovery
using signal-S phases in the
[clean-removal interface](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md).
Indeed $`U=H_a\otimes D_f`$ has accepted block
$`D_f/\sqrt2`$ for any Boolean f, although compiling D_f may
retain the cost one sought to remove. A source with multiple distinct
noncommuting phases is also outside the screen.

The missing operation is still a native constant-normalization block
with the complete indicator action, linear T-count and dirty width,
quadratic Clifford count, and the requested smaller T-depth. Recovery
and clean removal are paid once that source exists; neither test supplies
it or rules it out in the unrestricted model.

## 4. A flat one-band block reduces to two logical native bands

Keep one signal, a nonempty remaining register, and one physical T layer
between arbitrary Clifford endpoints. If its accepted block is flat,
$`A^\dagger A=I/2`$, the normalized logical unitary
$`L=\sqrt2A`$ has a zero-clean native realization with at most
two T layers on all remaining logical and dirty registers. In particular,
if L is a literal Boolean diagonal phase, its degree is at most three;
if L is a complete Boolean XOR query, its degree is at most two.

Here is the reduction, retaining the original independent native axes.
Use Q,R,M from Section 3. Flatness says the compression of M to the
Q-positive sector is zero. Since M is a Hermitian unitary and both Q
sectors have equal dimension, its off-diagonal block is unitary and its
other diagonal block is zero. Thus $`\{Q,M\}=0`$ on the
whole space. The Pauli expansion of $`M=D^\dagger RD`$ has
nonzero support exactly
$`R+\mathrm{span}\{P_j:P_j\text{ anticommutes with }R\}`$:
the physical band axes P_j are independent, so its subset terms cannot
cancel. Conjugating by Q now shows that Q anticommutes with R and
commutes with every axis active on R. No band axis is active on both
clean Paulis.

A fixed Clifford change of coordinates sends Q to $`Z_a`$ and
R to $`X_a`$. Every band axis then belongs to one of three groups,
with its actual Hermitian sign retained:

```math
X_aP_j,\qquad Z_aQ_j,\qquad I_aL_j.
```

All the original axes commute. In the compression
$`\langle+_a|D|0_a\rangle`$, the X-group acts on the left as
its logical Pauli rotations, the Z-group acts on the right as its
logical rotations, and the silent group can be assigned to either.
The overlap $`\langle+|0\rangle=1/\sqrt2`$ is literal.
Returning to the original endpoint coordinates only adds logical
Cliffords. Equivalently, the exact Clifford
$`(I+RQ)/\sqrt2`$ maps the Q-positive code to the R-positive
code, and their fixed-code overlap is that Clifford divided by
$`\sqrt2`$.

If both active groups are nonempty, the projected X-plus-silent list
is independent. Otherwise a relation would lift to $`X_a`$ in
the original list, which is impossible: that operator anticommutes
with an original Z-group axis, while every original axis commutes
with that axis. The projected Z list is independent by the same argument
using any X-group axis. Each projected list therefore emits as one
physical T layer after a logical Clifford change of coordinates.

If one active group is empty, all projected axes commute and the
projection has kernel of dimension at most one. Delete one axis if a
relation exists. The remaining independent list emits in one band, and
the deleted nonidentity Pauli rotation emits in a second band. A
projected identity axis contributes only an eighth-root scalar, paid by
a constant Clifford word; for example $`(HS)^3=e^{i\pi/4}I`$.
Thus the two-band statement includes these cases as well. The Boolean
degree assertions now follow from the retained
[full-space two-layer diagonal theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#1-model-and-full-space-diagonal-theorem).

This reduction is robust at fixed requested error. Let V be a unitary
target on the entire nonsignal space. By (8),

```math
\|A-V/\sqrt2\|\le\epsilon
 \lt\gamma:=\cos(\pi/8)-1/\sqrt2
\quad\Longrightarrow\quad A^\dagger A=I/2.
\qquad(9)
```

Indeed, every singular value must be within epsilon of
$`1/\sqrt2`$, and gamma is the distance to the nearest other
allowed value. The normalized L is then an exact two-layer unitary and
$`\|L-V\|\le\sqrt2\epsilon`$. The retained
[two-layer Boolean-query floor](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction)
is $`\delta=\sin(\pi/8)/(2\sqrt2)`$, even after optimizing
a common scalar. Therefore every nonquadratic complete XOR target obeys

```math
\|A-V/\sqrt2\|\ge\frac{\sin(\pi/8)}4.
\qquad(10)
```

The asserted floor is smaller than gamma, so (9) justifies its use.
For the complete indicator, conjugating the entire output bank by
Hadamards gives $`(-1)^{Y_x}`$, a Boolean phase of degree p+1.
Thus (10) applies when $`p\ge3`$. It includes arbitrary dirty
width, arbitrary endpoint address/output mixing, and all reference
inputs. This closes the one-physical-band source at both tested
normalizations, not the normalized-block route. The next source question
allows at least two genuinely noncommuting native bands and still owes
its full count, width, Clifford and depth recurrence.

## 5. A precise next interface for two source bands

A two-band source can be written $`U=U_B^\dagger U_A`$, where
U_A,U_B are actual one-band encoders with Clifford endpoints. With
$`J|\psi\rangle=|0\rangle_a|\psi\rangle`$, put

```math
K=J^\dagger U_B^\dagger U_AJ,\qquad
A=U_AZ_aU_A^\dagger,\qquad B=U_BZ_aU_B^\dagger.
```

A and B are Hermitian Clifford involutions. Exact flatness
$`K^\dagger K=I/2`$ is equivalent to $`\{A,B\}=0`$,
by the equal-sector block argument in Section 4. A constant-accuracy
version is available even though the discrete one-band spectrum (8)
need not apply to U.

Suppose a unitary target V on the complete nonsignal space satisfies
$`\|K-V/\sqrt2\|\le\epsilon`$. Then

```math
\|\{A,B\}\|=4\|K^\dagger K-I/2\|
 \le4(\sqrt2\epsilon+\epsilon^2).                 \qquad(11)
```

For the equality, express B in the two equal-dimensional A sectors.
Its diagonal blocks D,F have the same norm: unitarity gives
$`D^2=I-EE^\dagger`$ and $`F^2=I-E^\dagger E`$ for the
square off-diagonal block E. The anticommutator has diagonal blocks
2D and minus 2F, while the compressed Gram defect is D/2 in the
actual U_A coordinates.

Two distinct Hermitian Cliffords have operator distance at least
$`\sqrt{2-\sqrt2}`$. If they differ only by scalar, Hermiticity
makes that scalar minus one and their distance is two. Otherwise their
normalized Choi states are distinct stabilizer states, whose overlap
has magnitude at most $`1/\sqrt2`$; the normalized
Hilbert--Schmidt distance already gives the stated bound. Apply this
to the actual Hermitian Cliffords ABA and minus B. Equation (11) shows

```math
4(\sqrt2\epsilon+\epsilon^2)\lt\sqrt{2-\sqrt2}
\quad\Longrightarrow\quad \{A,B\}=0.             \qquad(12)
```

In particular, $`\epsilon\le1/10`$ suffices. Thus a sufficiently
accurate two-band source at this amplitude must already be exactly
flat; its normalized logical unitary $`\sqrt2K`$ is within
$`\sqrt2\epsilon`$ of V. This is not yet a degree bound on
that logical unitary.

Nor is aligning the flat sectors a new unpaid operation. For an actual
anticommuting pair define $`R_P=e^{-i\pi P/4}`$. Then

```math
W=\frac{I+BA}{\sqrt2}=R_B^\dagger R_AR_B,
\qquad WAW^\dagger=B.                             \qquad(13)
```

Each R uses its actual one-band encoder, a signal-S or S-dagger
Clifford, and the actual inverse. Its constant eighth-root scalar can
be emitted by a Clifford word as above. W uses six one-band encoder
occurrences. If the encoders have linear T-count and quadratic
Clifford count, W preserves those orders, with constant depth and no
additional clean work. With $`\omega=e^{i\pi/4}`$, cancelling
adjacent actual inverses gives the literal completed word

```math
U_B^\dagger WU_A=\omega^{-1}S_a^\dagger U S_aU^\dagger S_aU.
```

In promotion, the leading signal S-dagger commutes with the controlled
bilinear probe and therefore cancels in each source conjugation; the
fixed scalar cancels as well. What remains is exactly the three-call
raw signal-S recovery already paid in the clean-removal theorem.
Thus its subsequent zero-clean promotion retains that theorem's complete
forward/actual-inverse ledger.

The next concrete question is therefore to construct actual native
anticommuting reflection pairs whose complete normalized logical
transition approximates the indicator, with linear source count and
width. Anticommutation guarantees flat success; it does not determine
the logical transition. Neither an arbitrary supplied Clifford pair nor
a formal alignment substitutes for the required original one-band
encoders. Address/output mixing is still allowed throughout.

### 5.1. Two source bands at half amplitude reduce to two clean bands

This continuation after PR #198 closes the alternative half-amplitude
normalization for two physical source bands. The argument received an
independent analytic audit; it uses the retained two-layer theorem and
adds no new numerical or native-circuit fixture. Its proof is distinct
from the still-open inverse-square-root-two transition in Section 5.
A constant-normalization block remains a paid sufficient interface when
more source bands are allowed. The complete-frame frontier is unchanged.

Let a be the sole initialized signal and let J insert a equal to zero.
The nonsignal space includes the entire arbitrary address, bank and
dirty-work input; all identities extend by identity on external
references. Factor an actual two-band native source as

```math
U=U_B^\dagger U_A,\qquad
K=J^\dagger U_B^\dagger U_AJ,\qquad
A=U_AZ_aU_A^\dagger,\quad B=U_BZ_aU_B^\dagger.
```

Each encoder has at most one physical T/T-dagger layer, with arbitrary Clifford endpoints. Thus A and B are actual Hermitian Clifford involutions. Both eigenspaces of each have the same dimension. No intermediate address-preservation or output-basis assumption is made.

For any unitary V on the complete nonsignal space, suppose

```math
\|K-V/2\|\le\epsilon,\qquad \epsilon\lt1/2.
\qquad\text{(H1)}
```

The rounding condition is

```math
8(\epsilon+\epsilon^2)\lt\sqrt{2-\sqrt2}.
\qquad\text{(H2)}
```

Under these conditions the derivation gives exact flatness K-dagger K equal to I/4, and a native one-clean circuit with at most two physical T bands whose exact occupied unitary is 2K. Its error against V is at most 2 epsilon, with exact signal return and no new clean resource.

Use U_A coordinates, where A is Z on the signal, and write

```math
B=\begin{pmatrix}D&E\\E^\dagger&F\end{pmatrix},\qquad
K^\dagger K=(I+D)/2.
```

The singular values of K lie in the interval [1/2 minus epsilon, 1/2 plus epsilon]. By (H1) this interval lies strictly inside (0,1). Hence the spectrum of D lies in (-1,1), so E is invertible: unitarity of B gives

```math
EE^\dagger=I-D^2,\qquad DE+EF=0.
```

In the left polar decomposition E equal to R times V_E, R is the positive square root of I minus D squared and commutes with D. The second identity therefore gives F equal to minus V_E-dagger D V_E. Consequently the two diagonal blocks of the following defect are unitarily equivalent:

```math
\{A,B\}+I=
\begin{pmatrix}2D+I&0\\0&I-2F\end{pmatrix}.
```

It follows that

```math
\|\{A,B\}+I\|
 =4\|K^\dagger K-I/4\|
 \le4(\epsilon+\epsilon^2).
\qquad\text{(H3)}
```

The invertibility step is essential. Equal sector dimensions alone do not justify replacing this shifted norm by a single compressed block.

#### Rounding to the exact order-three identity

Set W equal to AB. This is an actual Clifford, W-dagger equals BA, and A W A equals W-dagger. Put E_0 equal to W plus W-dagger plus I. Then

```math
W^2+W+I=WE_0,\qquad
W^3-I=(W-I)WE_0,
```

so (H3) implies

```math
\|W^3-I\|\le8(\epsilon+\epsilon^2).
\qquad\text{(H4)}
```

Any Clifford that is not a scalar has distance at least square root of (2 minus square root of 2) from I: its normalized Choi stabilizer state has overlap of magnitude at most 1/square root of 2 with that of I, and normalized Hilbert--Schmidt distance is bounded by operator distance. Thus (H2) and (H4) force W cubed to be scalar. Conjugating that scalar by A and using A W cubed A equal to (W cubed)-dagger forces it to be real, hence plus or minus I. Equation (H4) is strictly below two, excluding minus I. Therefore W cubed equals I.

Its possible eigenvalues are 1 and the two primitive cube roots of one. E_0 consequently has possible eigenvalues 3 and 0. Its norm in (H3) is strictly below 3 under (H2), so the eigenvalue 1 of W is absent. This proves

```math
W^2+W+I=0,\qquad \{A,B\}=-I,\qquad K^\dagger K=I/4.
\qquad\text{(H5)}
```

This rounding uses actual Clifford structure; it is not a statement about arbitrary pairs of reflections or arbitrary multiband encoders.

#### The literal polar aligner is already Clifford

Equation (H5) gives the unusually simple identity

```math
R=I+BA=-AB.
\qquad\text{(H6)}
```

Thus R is a literal unitary Clifford, including its displayed minus sign. It satisfies R A R-dagger equal to B. More directly, for a vector v in the A-positive sector,

```math
Rv=(I+B)v=2P_Bv,\qquad P_B=(I+B)/2.
```

The native completed word

```math
V_{\mathrm{clean}}=U_B^\dagger R U_A
\qquad\text{(H7)}
```

therefore has the exact complete-code identity

```math
V_{\mathrm{clean}}J=J(2K).
\qquad\text{(H8)}
```

There is no rejected amplitude and no residual input-dependent phase in (H8). It includes all arbitrary dirty-work inputs. The actual adjoint of (H7) supplies the inverse. Neither controlled encoders nor postselection are emitted.

Let T_A,G_A and T_B,G_B be common forward/actual-inverse source
bounds and let $`n_{\rm phys}`$ be the total physical width. Compile
R as a Clifford from its tableau, rather than executing its T-containing
source expressions. This costs $`O(n_{\rm phys}^2)`$ elementary
Cliffords and no T gates; its literal minus sign costs a constant
Clifford word. The same applies to its actual inverse. Equation (H7)
therefore has common bounds

```math
T\le T_A+T_B,\qquad D_T\le2,\qquad
G\le G_A+G_B+O(n_{\rm phys}^2).
```

It uses the same sole clean signal and all original dirty inputs.
When $`n_{\rm phys}=O(m)`$, this fits the quadratic Clifford
budget. The tableau is explicit from the source; it is not a supplied
free operation.

Generic promotion with one clean bit would subsequently use eight copies of this recovered source or its actual inverse and its already-paid bilinear tests. Its complete promoted error is at most 32 epsilon. The crucial depth-screening conclusion, however, already follows from (H7), before promotion.

#### Consequence for the complete indicator

Let

```math
\delta=\frac{\sin(\pi/8)}{2\sqrt2}.
```

The retained one-clean two-layer Boolean theorem has full initialized-isometry error floor delta for a nonquadratic Boolean XOR query. The same proof applies to vector outputs: Hadamards on the entire bank turn the complete indicator into the Boolean phase minus one to the power Y_x, of degree p plus one. For p at least three, that phase has degree greater than three.

If a two-band source had

```math
\|K-I_p/2\|\lt\frac\delta2
 =\frac{\sin(\pi/8)}{4\sqrt2}\approx0.0676495,
\qquad\text{(H9)}
```

then (H2) holds, the exact recovery above is available, and (H8) has error strictly below delta. This contradicts the retained theorem. Thus the right side of (H9) is a constant lower bound for this particular half-normalized two-band source. A freely optimized common phase does not evade the retained floor. It permits arbitrary dirty width, arbitrary intermediate address/bank mixing and references. It is not a growing-depth lower bound or a restriction on three or more source bands.


### 5.2. Paused discovery and the next complete transition

The continuation stopped for a fresh-workspace handoff, not because the
inverse-square-root-two mechanism was resolved. Its first unpaid
operation is a family of actual one-band encoders whose anticommuting
clean-image reflections have a complete normalized transition
approximating the indicator, with linear T-count/dirty width and
quadratic Clifford count. Exact flatness and a paid aligner alone do
not select that transition. All address, bank, dirty and reference
inputs, the literal phase and the actual inverse remain part of it.

One separated construction is a dead end within the existing one-band
screen. On disjoint logical banks R,S let C_R,C_S copy computational
parities r,s into the signal, and let V_R,V_S be one-band words supported
on their respective banks. With

```math
U_A=V_RC_R,\qquad U_B=V_SH_aC_S,
```

the clean-image reflections anticommute, but direct contraction gives

```math
\sqrt2\langle0_a|U_B^\dagger U_A|0_a\rangle
   =V_R(-1)^{rs}V_S^\dagger.
```

More decisively, the source itself is
$`C_S^\dagger H_a(V_S^\dagger V_R)C_R`$: the disjoint
T bands merge into one physical band between Cliffords. The Section 4
screen therefore applies. This construction can be reopened only with
genuinely overlapping logical dynamics that prevent that merge; high
reflection rank by itself is insufficient.

A structural lead, not a degree theorem, is to write the symplectic
maps of the two reflections as $`I+N_A,I+N_B`$. Each N squares
to zero, and anticommutation implies
$`N_AN_B=N_BN_A`$. The product image lies in the intersection
of their actual active isotropic axis spans. Nonzero product or
nontrivial shared-axis geometry is a candidate source of overlapping
dynamics, not a proved necessary condition for every successful source.
The next task is to compute a complete transition from such actual
encoders and expose its resource recurrence before drawing a depth
conclusion. More than two source bands also remain allowed.

## 6. Bounded receipts

[Exact fixtures](../../tests/test_joint_fourier_word_limits.py) check
the Fourier/derivative identity on every three-bit Boolean function,
literal small native fourth words and actual inverses, and the
one-signal compression polynomial on native examples with nonquadratic
phase and noncommuting Clifford endpoints. These finite checks do not
replace the dimension-independent proofs or establish a general native
depth bound.
