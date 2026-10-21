# Clifford averages cannot select the complete lookup eigenspace

[Complete observable contract](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Scalar expander and its earlier screens](EXPANDER_ZERO_REFLECTION.md) ·
[Current revision](READOUT_AND_NATIVE_SOURCE_REVISION_RESULT.md)

The scalar zero test has a pointwise fixed space that is easy for
Clifford generators to preserve. The complete lookup has a different
geometry. No nontrivial Clifford can fix its entire positive eigenspace
pointwise. More strongly, even an arbitrarily wide Clifford has small
normalized trace correlation with the lookup observable. Thus a
normalized Clifford average cannot uniformly spectrally separate its
two whole eigenspaces by a constant gap.

This removes one proposed escape from the earlier joint-relabeling and
bounded operator-Schmidt-rank screens: **wide Clifford coupling alone**
does not supply a constant-size Clifford-average filter for the complete
bank. Wide non-Clifford coupling, signed combinations with a growing
coefficient budget, and general native circuits remain outside the
conclusion. In particular, this is not a native T-depth lower bound.

## 1. Complete space and normalized quantities

Let p be at least two, let m equal two to the p, and put

```math
\Omega_p=\sum_{x\in\mathbb F_2^p}|x\rangle\langle x|
                 \otimes Z_{Y_x}\otimes I_H,
\qquad P_\pm=(I\pm\Omega_p)/2.
```

Y has m arbitrary output bits. H contains any number of arbitrary
helper bits. Both eigenspaces have half the full dimension D. Write
$`\tau(A)=D^{-1}\mathrm{tr}A`$ and
$`\tau_\pm(A)=2\tau(P_\pm A)`$. Literal phases are retained.
The target and every estimate below include identity on H, rather than
an initialized helper promise. No physical trace operation is supplied:
normalized traces are only witnesses for full-operator claims.

Define the sharp correlation bound

```math
b_p=\frac{p+1}{2^p}.
```

For **every** Clifford C on this complete space, with arbitrary
address, bank and helper mixing and any common scalar phase,

```math
|\tau(\Omega_p C)|\le b_p.                              \qquad(1)
```

The proof is elementary stabilizer amplitude algebra and is included
below. No restriction to tensor-product branch operators, bounded
Schmidt rank, read-only addresses, or a fixed number of helper wires
is imposed.

## 2. The quadratic-phase correlation lemma

First consider a modulo-four quadratic phase in all bits $`x,z`$,
where $`z=(Y,H)`$, with arbitrary linear coefficients and even cross
coefficients. Decompose it as

```math
q(x,z)=q_A(x)+q_B(z)+2(Lx)\cdot z\pmod4.
```

Absorb the constant into q_A, so q_B is zero at zero. If
$`q_B(z)=\sum_i a_i z_i+2\sum_{i\lt j}b_{ij}z_i z_j`$, let B be
the symmetric binary matrix with diagonal $`a_i\bmod2`$ and
off-diagonal entries b_ij. Bitwise addition gives the literal identity

```math
q_B(z+t)-q_B(z)=q_B(t)+2z^{\mathsf T}Bt\pmod4.
```

For t in the radical $`\ker B`$, q_B(t) is even because its parity
is $`t^{\mathsf T}Bt`$. The same identity makes
$`i^{q_B(t)}`$ a real character of that radical. If B has binary rank
r, expanding the squared Walsh coefficient and summing first over z
therefore gives

```math
\left|\mathbb E_z i^{q_B(z)}(-1)^{\lambda\cdot z}\right|^2
 =2^{-\dim z}\sum_{t\in\ker B}
                       i^{q_B(t)}(-1)^{\lambda\cdot t}.
```

This is zero off one affine r-dimensional subspace S and equals
$`2^{-r}`$ on S. Thus the nonzero Walsh magnitudes are
$`2^{-r/2}`$. Here r may be odd; no Boolean-sign-only premise is used.

The bank label of $`Z_{Y_x}`$ is the unit vector e_x, extended by
zeros on H. Only addresses in

```math
\mathcal A=\{x:e_x+Lx\in S\}
```

can contribute to $`\mathbb E_{x,z}(-1)^{Y_x}i^{q(x,z)}`$. If this
set is nonempty, fix $`x_0\in\mathcal A`$. Its distinct vectors
$`e_x+e_{x_0}`$, for $`x\ne x_0`$, are linearly independent, and
all belong to

```math
\mathrm{im}L+\mathrm{dir}S.
```

Consequently $`|\mathcal A|-1\le p+r`$. Arbitrary helper
coordinates do not change that independence or dimension bound. Thus

```math
\begin{aligned}
\left|\mathbb E_{x,z}(-1)^{Y_x}i^{q(x,z)}\right|
 &\le \frac{\min(m,p+r+1)}{m}\,2^{-r/2}\\
 &\le \frac{p+1}{m}.                                    \qquad(2)
\end{aligned}
```

For the last inequality, the ratio between successive values of
$`(p+r+1)2^{-r/2}`$ is
$`(p+r+2)/(\sqrt2(p+r+1))`$, which is at most one for p at least
two and r at least zero. The address-only phase q_A can only cancel
contributions and was not omitted from the operator.

## 3. From quadratic phases to every Clifford diagonal

A diagonal Clifford phase has the form

```math
\gamma\,i^{q(v)},\qquad |\gamma|=1,\qquad
q(v)=\sum_i a_i v_i+2\sum_{i\lt j}b_{ij}v_i v_j\pmod4.       \qquad(3)
```

where v denotes all computational input bits. Equation (2) therefore
proves (1) for diagonal Cliffords, without a coefficient or phase loss.

The same argument covers a general Clifford. The affine stabilizer
amplitude fact used here is also the retained
[source-catalogue F44](../../docs/reference/SOURCE_CATALOGUE.md),
Dehaene--De Moor, Section V, Theorem 5(ii). Its normalized Choi
vector is a stabilizer state. Computational amplitudes of a stabilizer
state have affine support, constant magnitude on that support, and a
phase that is linear modulo four plus twice a Boolean quadratic.
These facts follow by putting its Pauli stabilizers into row-echelon
form: diagonal stabilizers specify the affine support, and each
remaining X translation fixes the amplitude ratio by a linear Pauli
phase. Restricting the Choi support to equal input and output labels
therefore gives either an identically zero diagonal or

```math
\langle v|C|v\rangle
 =a\,\mathbf1_{v\in D_0}\,i^{q(v)},
 \qquad |a|\le1,                                       \qquad(4)
```

for an affine subspace D_0. The phase is initially defined on D_0.
Choose affine coordinates on D_0 and an affine retraction onto those
coordinates to extend it to the whole bit space. Substitution of
Boolean affine forms preserves the displayed modulo-four quadratic
phase class; in particular, cross terms remain even. No circuit
implementing this extension is being assumed.

An affine subspace indicator has a character expansion

```math
\mathbf1_{v\in D_0}
 =2^{-s}\sum_{t\in(\mathrm{dir}D_0)^\perp}
                         (-1)^{t\cdot(v-v_0)},
 \qquad s=\mathrm{codim}D_0,                     \qquad(5)
```

whose total coefficient modulus is one. Multiplication by any of
these affine characters merely adds twice an affine Boolean function
to q in (3), still a permitted modulo-four quadratic. Applying (2)
term by term to (3)--(5), and using $`|a|\le1`$,
proves (1). This proof retains all diagonal phases, affine constraints
and arbitrary helper variables. It does not replace C by a supplied
classical permutation or by an unimplemented diagonal circuit.

The bound is attained. Let e_i now denote the i-th unit **address**,
and take the actual diagonal Clifford whose fixed-address action is

```math
C_x=Z_{Y_0}\prod_{i=1}^p(Z_{Y_0}Z_{Y_{e_i}})^{x_i}.
```

It consists of one Z and two address-to-bank CZ gates per i. Its bank
Pauli equals $`Z_{Y_x}`$ exactly at the p+1 addresses zero and e_i.
At every address of weight at least two its label is supported on
zero and the unit addresses, whereas $`Z_{Y_x}`$ has a different
label. Those fibers are trace-orthogonal. Consequently
$`\tau(\Omega_p C)=(p+1)/m`$, with identity on all helpers.

## 4. Exact pointwise fixed-space rigidity

For p at least two, a Clifford satisfying

```math
C P_+=P_+                                               \qquad(6)
```

is the identity, including its scalar phase. This exact statement is
stronger than what (1) alone gives at small p.

Fixing one computational basis vector makes a Clifford monomial:
each conjugated computational Z stabilizer must be diagonal, since a
Pauli with a nonzero X label cannot fix that vector. Hence its
computational permutation is affine. The good basis strings
$`Y_x=0`$ affinely span the full input-bit space. Indeed all
$`(x,0,h)`$ are good, and for each output coordinate u there is
an address x distinct from u for which both $`(x,0,h)`$ and
$`(x,e_u,h)`$ are good. The affine permutation fixing all these
strings is therefore the identity.

The remaining phase is diagonal Clifford. Setting Y to zero in (6)
makes every address-only and helper-only phase vanish. The coefficient
of a bank bit Y_u has the form
$`a_u+2\sum_i b_{ui}x_i`$ modulo four, after bank-helper terms
are removed by varying h at any x distinct from u. This affine phase
coefficient vanishes at every address other than u. A two-dimensional
address face through u has three other vertices; its affine
parallelogram identity forces the coefficient at u to vanish as well.
Every bank-bank quadratic coefficient vanishes by choosing an address
distinct from its two labels, which is possible since m is at least
four. Thus the entire phase is identically one.

There is a useful clean-code qualification. If a Clifford instead
satisfies $`CJ P_+=J P_+`$, where J inserts any number of
computational-zero signal bits, the same argument on that affine code
proves $`CJ=J`$ on the **entire** logical code. C can still act
away from that code. This does not supply a full-space identity for
arbitrary signal inputs and is not used as one below.

## 5. Positive or lazy Clifford averages

Let

```math
A=\sum_j w_j C_j,\qquad w_j\ge0,\qquad \sum_j w_j=1,
\qquad H=I-\mathrm{Re}A.
```

H is positive semidefinite. Every branch may be a different, arbitrarily
wide Clifford; coherent address motion and common phases are permitted.
Equation (1) gives

```math
|\tau_-(H)-\tau_+(H)|
 =2|\mathrm{Re}\tau(\Omega_p A)|\le2b_p.           \qquad(7)
```

In particular, if all vectors in the positive sector have H energy at
most eta and all vectors in the negative sector have H energy at least
g, then

```math
g-\eta\le2b_p.                                         \qquad(8)
```

No invariance of either sector is needed for this necessary condition:
the two quadratic-form inequalities suffice. A normalized average
therefore cannot have a constant uniform spectral separation between
the whole sectors as p grows. A nonzero common baseline energy does
not remove this conclusion. The same contrast bound applies to a
signed or complex combination with coefficient l1 at most one.

An additional relative-energy estimate does not require normalizing
the positive weights. For a non-scalar Clifford, its normalized trace
has magnitude at most $`1/\sqrt2`$: its Choi vector and the
identity Choi vector are distinct projective stabilizer states, whose
overlap has at most that magnitude. Therefore

```math
\tau_+(I-\mathrm{Re}C)
 \ge\kappa_p:=1-1/\sqrt2-b_p\ge1/8\qquad(p\ge6),
\qquad
\tau_-(I-\mathrm{Re}C)\le2.                        \qquad(9)
```

A scalar Clifford contributes the same nonnegative energy in both
sectors, even if arbitrary common phases are allowed. Summing (9)
with arbitrary nonnegative weights gives

```math
\tau_-(H)\le\frac{2}{\kappa_p}\tau_+(H)\le16\tau_+(H),
\qquad H=\sum_jw_j(I-\mathrm{Re}C_j).              \qquad(10)
```

Thus the same low/high energy premises require
$`g\le2\eta/\kappa_p\le16\eta`$. For an exactly
frustration-free positive sum, (6) gives the sharper conclusion:
every positive-weight branch is identity, so the sum is zero.

## 6. The coefficient budget of a complete Clifford filter

For a literal decomposition $`E=\sum_j a_j C_j`$, put
$`\Lambda=\sum_j|a_j|`$. If its complete operator satisfies

```math
\|E-\alpha\Omega_p\|_{\rm op}\le\epsilon\lt|\alpha|,
```

then normalized trace duality and (1) imply

```math
\Lambda\ge\frac{|\alpha|-\epsilon}{b_p}.               \qquad(11)
```

This charges a coefficient budget, not a free LCU implementation. It
allows an arbitrary number of branches and arbitrary cancellations.
For a k-clean computational compression
$`K=J^\dagger E J`$, the same argument expands the clean projector
as $`2^{-k}\sum_tZ_A(t)`$. Each $`Z_A(t)C_j`$ is a Clifford,
and normalized trace conversion gives

```math
\|K-\alpha\Omega_p\|\le\epsilon
\quad\Longrightarrow\quad
\Lambda\ge\frac{|\alpha|-\epsilon}{2^k b_p}.           \qquad(12)
```

Every non-clean input remains arbitrary in (12). No return promise is
deduced merely from this compression.

The same witness applies to a partial-isometry Gram source. Let
$`K=J^\dagger E J`$, let E have Clifford coefficient budget Lambda,
and let $`\Pi_{\rm eq}=(I+\Omega_p Z_s)/2`$ for an arbitrary
signal bit s. If
$`\|K^\dagger K-\Pi_{\rm eq}/2\|\le\epsilon`$, form the
algebraic full-space expression $`F=4E^\dagger P_0E-I`$, where
$`P_0=JJ^\dagger`$. It has budget at most
$`4\Lambda^2+1`$ and its clean compression differs from
$`\Omega_p Z_s`$ by at most four epsilon. Multiplying the witness
by the Clifford Z_s reduces it to (12), so necessarily

```math
(4\Lambda^2+1)2^k b_p\ge1-4\epsilon.                   \qquad(12a)
```

This is a constraint on the supplied source representation, not a
physical implementation of F or a substitute for the source's actual
inverse and return proof.

Clifford coefficient budgets are submultiplicative under actual
products and unchanged by actual adjoints. A completed s-branch
Clifford SELECT on a computational selector has a decomposition with
budget at most s, by expanding each selector projector into Paulis.
A determinant-one selector phase with relative phase
$`e^{i\pi/3}`$ has budget at most two, since it is a common scalar
times $`I+(e^{i\pi/3}-1)P`$, and a computational projector P
has coefficient budget one. Consequently the literal cubing word

```math
S_{j+1}=S_j R S_j^\dagger R^\dagger S_j,\qquad N=3^j,
```

starting from a word of budget at most s, has the explicit budget

```math
\Lambda(S_j)\le s^N2^{N-1}.                             \qquad(13)
```

This is an algebraic bound for that actual word, not an alternate
implementation. Fixed-width selector operations and fixed numbers of
source/actual-inverse calls multiply it only by their corresponding
bounded factors. Equations (11)--(13) show why replacing the scalar
expander branches by arbitrary wide Cliffords does not yield a
constant-call, constant-amplitude complete lookup filter. With a
fixed branch budget per letter, (12) requires a number of such letters
of order p minus k minus log p. The actual number of calls is not an
unrestricted native-depth lower bound: jointly resynthesized words and
one wide T layer can have much larger Clifford coefficient budgets.

### A complete bank needs growing numbers of bare scalar zero tests

An exact zero reflection on any chosen set of wires has the algebraic
form $`R=I-2P_0`$. The computational-zero projector P_0 has a Pauli
expansion with total coefficient modulus one, so R has Clifford
coefficient budget at most three, independently of the tested width.
Identity on arbitrary other inputs, literal Clifford conjugation and
actual inversion do not enlarge this bound.

Consider an emitted word containing r actual approximate zero-reflection
calls, each within full-operator error zeta of its stipulated exact
reflection, interleaved with arbitrary wide Cliffords. The calls and
interlayers may mix the address, bank and helpers; no address-sharing
or fixed-target assumption is imposed. Replace just the r approximate
calls by their exact counterparts to obtain W_0. Telescoping on the
complete space gives $`\|W-W_0\|\le r\zeta`$, and
$`\Lambda(W_0)\le3^r`$. Consequently

```math
\|J^\dagger WJ-\alpha\Omega_p\|\le\epsilon
\quad\Longrightarrow\quad
3^r2^k b_p\ge|\alpha|-\epsilon-r\zeta,                 \qquad(14)
```

where J inserts k computational-zero signals and every other input is
arbitrary. Actual inverse calls have the same error and coefficient
budgets and each counts among the r occurrences. If the right side is
bounded below by a positive constant, then

```math
r\log_2 3\ge p-k-\log_2(p+1)-O(1).
```

Thus a constant number of paid scalar zero tests plus arbitrary wide
Clifford processing cannot supply a constant-amplitude complete bank
with a small clean signal. This is an actual-call count in the stated
representation, not native T-depth: many disjoint calls may run in one
layer, and a jointly resynthesized native word need not retain the
displayed call factorization. Wide non-Clifford interlayers are outside
the premise. The bound neither prices those interlayers for free nor
supplies the still-missing simultaneous bank schedule and return proof.

## 7. Revision and next paid operation

The scalar approximate zero reflection remains a valid paid component;
its sparse marked subspace is not the complete lookup half-space.
The new screen rules out a direct complete-bank replacement based on
a low-budget Clifford average, including wide Clifford branches that
escape the earlier Schmidt-rank screen. It neither invalidates the
scalar result nor changes the general complete-frame frontier.

The next constructive mechanism must supply more than wide Clifford
coupling. A genuinely wide non-Clifford address-bank interaction can
already have the required coefficient budget within the relaxed
polynomial overhead. The paid quadratic-query kernels are examples
of such available interactions, but their constant-fraction complete
readout is still unpaid. The remaining target is an actual whole-bank
word or compression, not a further bare-Clifford expander: simultaneous
native count, Clifford count, dirty width, depth, literal phase,
actual inverse, complete error and arbitrary-work return must all be
proved for that word. These representation screens do not exclude
address-moving overlapping two-band sources with growing clean count
at fixed error, higher-band mechanisms, or unrestricted sublogarithmic
native depth. The separate read-only two-band screen retains its own
stronger source-specific clean-count restriction.
