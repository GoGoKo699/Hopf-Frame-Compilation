# Involutory Gaussian loops and their unsquared scalar limit

[Complete observable source](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Closed Gaussian round limits](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md) ·
[Paid inverse extraction](GAUSSIAN_INVERSE_EXTRACTION.md) ·
[Growing-clean recovery](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag)

A symmetric involutory binary pencil gives a literal sixfold Fourier
loop that returns every arbitrary carrier input. Its scalar is a cubic
sign for affine coefficients. Retaining the unsquared Gauss phase does
not supply an asymptotic escape: if the matrix entries have Boolean
degree at most D on the full external input cube, any constant-root
normalization of that phase to a sign has degree at most $`9D`$,
independently of carrier width. This conclusion holds even if the
in-place matrix permutation is granted for free.

These are limits of this involutory scalar mechanism. They do not give
an unrestricted native T-depth lower bound, rule out all Gaussian
completions, or construct a source in the growing-clean window.

## 1. Conventions and the literal third power

Let A be an n-by-n symmetric binary matrix with $`A^2=I`$ over
$`\mathbb F_2`$. The n carrier wires are arbitrary. External registers
are read-only in this calculation and may be in arbitrary superpositions
or entangled with the carrier and a reference. Define

```math
q_A(z)=\sum_i A_{ii}z_i+2\sum_{i\lt j}A_{ij}z_i z_j\pmod4,
\qquad D_A|z\rangle=i^{q_A(z)}|z\rangle,
\qquad F_A=H^{\otimes n}D_A.
                                                               \qquad(1)
```

Write $`P_A|z\rangle=|Az\rangle`$ for the distinct in-place
permutation, and put

```math
G_A(s)=\sum_z i^{q_A(z)}(-1)^{s\cdot z},
\qquad \gamma_A=2^{-n/2}G_A(0).
```

The refinement identity is

```math
q_A(z+t)=q_A(z)+q_A(t)+2z^{\mathsf T}At\pmod4.
```

Shifting by As and using character orthogonality gives

```math
G_A(s)=i^{-q_A(As)}G_A(0),\qquad |\gamma_A|=1.          \qquad(2)
```

The two refinements $`q_A(Az)`$ and $`q_A(z)`$ have the same polar
form and the same reduction modulo two. Their difference is therefore
twice a linear Boolean function:

```math
q_A(Az)-q_A(z)=2\ell\cdot z,
\qquad
\ell_i=\frac{q_A(Ae_i)-A_{ii}}2\pmod2.
                                                               \qquad(3)
```

The numerator is even, and its value modulo four determines the quotient.
Applying (3) twice gives $`A\ell=\ell`$. Direct evaluation of the kernels
then gives, with every carrier column retained,

```math
\begin{aligned}
\langle y|F_A^2|x\rangle
 &=2^{-n/2}\gamma_A i^{q_A(x)-q_A(A(x+y))},\\
\langle w|F_A^3|x\rangle
 &=\gamma_A i^{q_A(x)-q_A(Ax)}
   2^{-n}\sum_y(-1)^{(w+Ax+\ell)\cdot y}\\
 &=\gamma_A(-1)^{\ell\cdot x}\,\mathbf1_{w=Ax+\ell}.
\end{aligned}
```

Thus the exact full-space identity, including its Pauli correction and
literal scalar, is

```math
\boxed{F_A^3=\gamma_A X^\ell P_A Z^\ell.}             \qquad(4)
```

No initialized impulse, postselection, or stationary-phase assumption
enters this identity. The expression for $`\ell`$ is not a free native
implementation of either parameter-dependent Pauli correction.

## 2. The sixth power squares away the Arf sign

Binary orthogonality gives $`A\mathbf1=\mathbf1`$: every row has odd
self-inner-product. Evaluating (3) at $`z=\mathbf1`$ shows that
$`\ell\cdot\mathbf1=0`$. Squaring (4) therefore yields

```math
F_A^6=\gamma_A^2(-1)^{\ell\cdot\ell}I=\gamma_A^2I.
                                                               \qquad(5)
```

Put $`d=\mathrm{diag}(A)`$. Complex conjugation of the Gauss sum gives,
for any nonsingular symmetric binary A,

```math
\overline{G_A(0)}=G_A(d)
 =i^{-q_A(A^{-1}d)}G_A(0),
\qquad \gamma_A^2=i^{q_A(A^{-1}d)}.
```

In the involutory case $`Ad=d`$. Indeed, symmetry gives
$`\mathrm{diag}(A^3)=A\,\mathrm{diag}(A)`$ over $`\mathbb F_2`$,
while $`A^3=A`$. Hence

```math
\boxed{(H^{\otimes n}D_A)^6=i^{q_A(d)}I,\qquad
q_A(d)=\sum_i d_i+2\sum_{i\lt j}A_{ij}d_i d_j\pmod4.}
                                                               \qquad(6)
```

If the entries of A have Boolean degree at most D in external variables
$`\xi`$, the exponent in (6) has a multilinear polynomial over
$`\mathbb Z/4\mathbb Z`$ of degree at most $`3D`$. A degree-D Boolean
bit has an integer lift modulo four of degree at most $`2D`$; in the
second sum only the Boolean cubic products matter because of their
factor two. Also $`\mathrm{tr}(A)=n\pmod2`$, since $`(A+I)^2=0`$.
Consequently $`q_A(d)-n`$ is always even, and coefficientwise division
in the unique multilinear expansion gives

```math
i^{q_A(d)}=i^n(-1)^{h(\xi)},\qquad \deg h\le3D.       \qquad(7)
```

The fixed factor $`i^n`$ remains part of the literal circuit identity.
Removing it here is a mathematical comparison of sign functions, not
an uncharged native phase correction. For affine coefficients in all
of $`(x,Y)`$, this sign cannot be the complete addressed bit $`Y_x`$
when $`p\ge3`$, because that bit has degree $`p+1`$.

## 3. A doubled-block family computes trace rather than rank

For arbitrary symmetric r-by-r binary B, the family

```math
A=I+\begin{pmatrix}B&B\\B&B\end{pmatrix}              \qquad(8)
```

is symmetric and involutory. Write $`z=(u,v)`$, $`w=u+v`$, and
$`b=\mathrm{diag}(B)`$. As an identity modulo four,

```math
q_A(u,v)=|w|+q_B(w)+2b\cdot w
          +2\sum_i(1-w_i)v_i.
```

Summing over v kills every w except $`\mathbf1`$, so

```math
\gamma_A=i^{r+q_B(\mathbf1)+2\mathrm{tr}(B)},\qquad
F_A^6=(-1)^{r+\mathrm{tr}(B)}I.                     \qquad(9)
```

For $`r=2`$, both $`B=0`$ and the rank-one all-ones B give $`+I`$
in the sixth power, whereas $`B=\mathrm{diag}(1,0)`$ gives $`-I`$.
Rank parity is not the computed scalar. This family also has a fixed
basis shear implementation of its in-place map. Let E consist of r
CNOTs from v to u, so $`E(u,v)=(u+v,v)`$. Then

```math
P_A=E\,Q_B\,E,\qquad Q_B(w,v)=(w,v+Bw).              \qquad(9a)
```

Thus an actual returned implementation of $`Q_B`$ supplies $`P_A`$
with its T-count, T-depth, dirty reservation and full error unchanged,
and with $`2r`$ additional Cliffords. The actual inverse reverses this
word. This paid component still leaves the low-degree phase (9).

## 4. A characteristic-vector formula for the unsquared phase

The following identity fixes all sign and determinant conventions.
Let M be an integer symmetric matrix with odd determinant, and choose
an integer characteristic vector w satisfying
$`Mw=\mathrm{diag}(M)\pmod2`$. For every odd signed integer a, define

```math
\chi_2(a)=(-1)^{(a^2-1)/8},\qquad \omega=e^{\pi i/4}.
```

Thus $`\chi_2(a)=+1`$ for $`a=\pm1\pmod8`$, and it is $`-1`$
for $`a=\pm3\pmod8`$, including negative determinants. Then

```math
2^{-n/2}\sum_{z\in\mathbb F_2^n}i^{z^{\mathsf T}Mz}
 =\chi_2(\det M)\,\omega^{w^{\mathsf T}Mw}.          \qquad(10)
```

Here and below, integer products are used when an exponent is reduced
modulo four or eight. A direct proof is short. Congruence modulo eight
by an odd-determinant integer matrix reduces M to odd one-dimensional
blocks and even-diagonal two-dimensional blocks. An odd diagonal pivot
can be split by elimination modulo eight. If all diagonals are even,
nonsingularity supplies an odd off-diagonal entry, whose invertible
two-dimensional block can be split in the same way. These steps finish
because the complementary block remains nonsingular modulo two.

For an odd scalar block $`[a]`$,

```math
\frac{1+i^a}{\sqrt2}=\chi_2(a)\omega^a.
```

For the block

```math
\begin{pmatrix}2a&b\\b&2c\end{pmatrix}
```

with b odd,
the normalized sum is $`(-1)^{ac}`$, the characteristic vector is zero
modulo two, and the determinant is $`4ac-1\pmod8`$. Thus (10) holds
for both block types and for their orthogonal sums.

Both sides are invariant under the congruence used above. The Gauss sum
is unchanged by the induced invertible binary variable change: replacing
an integer vector by its binary lift changes its quadratic value by a
multiple of four. The determinant acquires an odd square, whose
$`\chi_2`$ is one. Finally, replacing a characteristic lift w by
$`w+2v`$ changes its norm by

```math
4(v^{\mathsf T}Mw+v^{\mathsf T}Mv)=0\pmod8.
```

The same calculation compares characteristic vectors across a
congruence. This proves (10) without choosing a parameter-dependent
native basis. The congruence is only a proof device.

For primary-source context, [Taylor, *Gauss Sums in Algebra and
Topology*, Theorem 2.8](https://arxiv.org/pdf/2208.06319) gives the
corresponding mod-eight Gauss formula and reduced-block calculation.
Its $`\sigma_2(M)`$ is $`w^{\mathsf T}Mw\pmod8`$ on the blocks
above, and its determinant correction is $`\chi_2(\det M)`$.
The elementary argument here supplies the specific form needed below.

Apply (10) to the integer zero-one lift $`\widetilde A`$ of the
involutory binary matrix A. Since $`Ad=d`$, one may take w=d. Therefore

```math
\boxed{\gamma_A=\omega^{d^{\mathsf T}\widetilde A d}
 (-1)^{(\det(\widetilde A^{\,2})-1)/8}.}             \qquad(11)
```

The square in (11) is an integer matrix product. It is not replaced by
the binary identity before taking its determinant.

## 5. Low-degree involutory coefficients cannot hide a high-degree sign

Assume A depends on external bits $`\xi`$, is symmetric and involutory
at **every point of their full Boolean cube**, and every entry has
Boolean degree at most D. Suppose there is a single constant
$`c\in\mathbb Z/8\mathbb Z`$ such that

```math
\gamma_{A(\xi)}=\omega^c(-1)^{f(\xi)}.
```

The constant root is independent of every external bit. It is not
assumed to be $`i^n`$, nor may a parameter-dependent choice of root be
discarded. Then

```math
\boxed{\deg f\le9D.}                               \qquad(12)
```

First, an integer monomial in at most r Boolean entry functions, each
of degree at most D, has a multilinear lift modulo $`2^s`$ of degree
at most $`(r+s-1)D`$. To see this, express each entry as an XOR of its
ANF monomials and use

```math
u_1\mathbin\oplus\cdots\mathbin\oplus u_N
 =\sum_{\varnothing\ne S\subseteq[N]}
   (-2)^{|S|-1}\prod_{j\in S}u_j.
```

In a product of r such expansions, a surviving term modulo $`2^s`$
contains at most $`r+s-1`$ ANF factors. Multilinear reduction can only
decrease its degree. Summing integer monomials preserves the bound.

Set $`C=\widetilde A^{\,2}-I`$. Its entries are pointwise even.
Expanding the determinant by principal minors therefore gives

```math
\det(\widetilde A^{\,2})-1
 \equiv\sum_{1\le|S|\le3}\det C[S]\pmod{16}.         \qquad(13)
```

Every larger minor is divisible by sixteen. The right side of (13)
is an integer polynomial of degree at most six in the zero-one entries
of A. The lifting bound with $`s=4`$ gives degree at most $`9D`$
in the external bits. Its values are divisible by eight, because the
square of any odd integer is one modulo eight. The unique multilinear
coefficients modulo sixteen are also divisible by eight: each is an
integer alternating sum of function values. Dividing those coefficients
by eight proves that

```math
s_A(\xi)=\frac{\det(\widetilde A^{\,2})-1}{8}\pmod2
\quad\text{has degree at most }9D.                  \qquad(14)
```

Similarly, $`t_A=d^{\mathsf T}\widetilde A d\pmod8`$ has a
multilinear lift of degree at most $`5D`$: its integer expression
has degree at most three in matrix-entry bits. Equation (11) and the
constant-root hypothesis imply $`t_A-c=0\pmod4`$ everywhere.
Coefficientwise division modulo eight shows that
$`(t_A-c)/4\pmod2`$ has degree at most $`5D`$. Adding it to
$`s_A`$ proves (12).

The full-cube condition is essential to these divisibility arguments.
An involution identity valid only on a promised one-hot output bank,
an initialized carrier, or selected external inputs does not satisfy
this theorem or the complete arbitrary-bank source contract.

For p address bits and an arbitrary $`2^p`$-bit output bank, the required
scalar bit is

```math
Y_x=\sum_{u\in\mathbb F_2^p}Y_u
 \prod_{j:u_j=1}x_j\prod_{j:u_j=0}(1+x_j)
\quad\text{over }\mathbb F_2,
\qquad \deg Y_x=p+1.                                \qquad(15)
```

Thus this unsquared scalar mechanism requires $`p+1\le9D`$, regardless
of carrier width. In particular fixed-degree coefficients cannot supply
the growing-address source, even if $`P_A`$ and the corrections in (4)
were free. The coefficient degree D is not native T-depth.

There is also a fixed-accuracy version for the ideal scalar. Distinct
eighth roots have distance at least $`2\sin(\pi/8)`$. If a common
phase $`e^{i\theta}`$ satisfies

```math
\sup_\xi|\gamma_{A(\xi)}-e^{i\theta}(-1)^{f(\xi)}|
 \lt\sin(\pi/8),                                    \qquad(16)
```

each of the two opposite disks contains at most one eighth root;
when f takes both values, those roots must be opposite. Hence the
constant-root form in (12) holds exactly. A constant f is already
covered. For an approximate emitted loop, its complete implementation
error must first be added to the desired scalar error before applying
(16). This does not infer quantization for an arbitrary approximate
source outside this Gaussian family.

## 6. Native obligations, actual inverses, and complete ledgers

An XOR query does not by itself pay the single-bank map $`P_A`$.
For $`Q_A:(c,t)\mapsto(c,t+Ac)`$, the chronological word

```math
(c,t)\longmapsto(c,t+Ac)\longmapsto(At,t+Ac)
 \longmapsto(At,Ac)\longmapsto(Ac,At)                 \qquad(17)
```

uses forward, reverse-direction, forward XOR queries and a final bank
SWAP. It pays $`P_A\otimes P_A`$ on two arbitrary banks. Treating
the second bank as zero would change the operator contract.

Assume each supplied primitive j has common bounds
$`(T_j,G_j,d_j,h_j)`$ for its actual forward circuit and actual adjoint,
returns its private arbitrary h_j-bit pool ideally, and has complete
operator error at most $`\eta_j`$. The primitive labels are D for
$`D_A`$, P for $`P_A`$, X and Z for the parameter-dependent Paulis
in (4), and Q for either declared query direction in (17). These are
assumptions about actual native implementations, including their
coefficient access; no loader or inverse is supplied by matrix notation.

Let $`R_A=X^\ell P_AZ^\ell`$. Exact cancellation gives
$`R_A^\dagger F_A^3=\gamma_AI`$. Its implementation uses the actual
adjoint of the supplied correction word. Sufficient bounds are

| Complete word | T-count | Clifford count | T-depth | Full error |
| --- | --- | --- | --- | --- |
| Sixth power | $`6T_D`$ | $`6G_D+6n`$ | $`6d_D`$ | $`6\eta_D`$ |
| Paired permutation (17) | $`3T_Q`$ | $`3G_Q+3n`$ | $`3d_Q`$ | $`3\eta_Q`$ |
| Unsquared cancellation | $`3T_D+T_P+T_X+T_Z`$ | $`3G_D+G_P+G_X+G_Z+3n`$ | $`3d_D+d_P+d_X+d_Z`$ | $`3\eta_D+\eta_P+\eta_X+\eta_Z`$ |

The sixth power needs one arbitrary n-bit principal carrier and h_D
private bits. The paired permutation needs two arbitrary n-bit principal
banks and h_Q private bits. A sufficient separate-pool allocation for
unsquared cancellation is n principal bits plus
$`h_D+h_P+h_X+h_Z`$ private bits. Each pool is reused sequentially.
If carrier wires themselves are allocated as private work, include their
n or 2n bits in the dirty-width budget.

One can also avoid a supplied single-bank P by pairing two third powers.
On each bank, $`X^\ell F_A^3Z^\ell=\gamma_AP_A`$. Cancel their
paired permutation using the actual inverse of (17). A conservative
serial ledger is

```math
\begin{aligned}
T&\le6T_D+2T_X+2T_Z+3T_Q,\\
G&\le6G_D+2G_X+2G_Z+3G_Q+9n,\\
d&\le6d_D+2d_X+2d_Z+3d_Q,\\
\epsilon&\le6\eta_D+2\eta_X+2\eta_Z+3\eta_Q.
\end{aligned}                                                    \qquad(18)
```

Two arbitrary principal banks and
$`h_D+h_X+h_Z+h_Q`$ separate private bits suffice. This word returns
the banks with scalar $`\gamma_A^2`$, losing the unsquared sign.
The ledger claims no parallel depth reduction by overlapping helper
pools.

All tables are zero-clean conditional ledgers. If a primitive instead
needs initialized helpers, its clean reservation and initialized-isometry
error must be declared separately and passed to the
[paid recovery and promotion recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md).
Nothing here provides those initialized wires for free. Dense matrix
entries also need not fit the simultaneous source count and width
envelopes merely because a phase or XOR query has a constant-depth
implementation.

The actual inverse reverses each emitted word, including every supplied
adjoint; it has the displayed common bounds. Its ideal scalar phase is
the complex conjugate of the forward phase. With exact primitives,
all final carrier and private-work returns are exact up to the stated
external scalar. With approximate primitives, unitary telescoping gives
the displayed full errors, including leakage, arbitrary references and
all return errors. No intermediate carrier reset is assumed.

## 7. Assessment and bounded evidence

The [bounded exact fixture](../../tests/test_involutory_gauss_source_limits.py)
retains literal phases and every carrier column for two fixed two-bit
and two fixed four-bit examples. It checks the third power, sixth power,
actual inverse, characteristic determinant formula, and doubled-block
trace value. It includes a nonzero Pauli correction and a negative
determinant. This is not a parameter-wide search or an asymptotic proof.

The first unsupported inference in the proposed unsquared route was that
an affine involutory pencil could hide an arbitrarily high-degree addressed
sign in its Gauss phase. Section 5 rules that out before any in-place
compilation question. Resolving the native $`P_A`$ obligation alone would
not improve the full source frontier through this mechanism.

General non-involutory forms with a separately paid inverse or completion,
quantum motion of external address registers, and other non-Clifford
interference remain distinct research questions. The earlier
[closed Gaussian round theorem](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
has its own actual-round hypotheses and return argument; the present
coefficient theorem does not replace it with a native-depth claim.
Neither result changes the unrestricted frontier or the full operator
contract of the intended approximate dirty indicator.
