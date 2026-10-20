# A whole-word native-depth limit for Majorana-preserving indicator circuits

[Complete bank-moving indicator](BANK_MOTION_WALSH_INDICATOR.md) ·
[Earlier approximate Gaussian source limit](APPROXIMATE_QUERY_FRAME_TRANSFER.md#4-constant-error-does-not-repair-the-gaussian-conjugator) ·
[Fixed Majorana-Walsh source](MAJORANA_WALSH_QUERY.md) ·
[Paid polynomial-overhead source fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Full resource frontier](../../docs/OPEN_PROBLEM.md)

A complete approximate indicator needs at least b native T bands at
sufficiently small fixed error if every actual Clifford boundary preserves
the address and the standard Majorana span of the arbitrary bank and
private work. This bounds the entire emitted word. No separately
implemented Walsh transform, routing unitary or source conjugator is
assumed. Arbitrarily many dirty helpers, intermediate mixing with those
helpers, and whole-word cancellation are included within this boundary
grammar.

The restriction on the actual Clifford boundaries is essential. This is
not an unrestricted Clifford+T depth lower bound. It shows that fusing the
bank-moving echo while retaining this Gaussian representation cannot
supply the shallow source required by the paid fusion theorem. R96 and
the unrestricted frontier remain unchanged. The next constructive
question is a complete source using non-Gaussian Clifford boundaries or
address motion across native bands, with its native costs, actual inverse
and complete return error paid.

## 1. Full operator contract and actual native boundaries

Let $`b\ge2`$, $`m=2^b`$, and let Y be an arbitrary bank of m
physical qubits indexed by $`u\in\mathbb F_2^b`$. Let W contain
any finite number w of arbitrary private qubits, and put $`Q=m+w`$.
The complete phase indicator is

```math
D_b=\sum_{x\in\mathbb F_2^b}|x\rangle\langle x|\otimes Z_{Y_x}.
                                                               \qquad(1)
```

Fix one physical ordering of all Q bank and work qubits. Its standard
Majoranas, acting trivially on the address, are

```math
\gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j,
\qquad 0\le j\lt Q.
                                                               \qquad(2)
```

They satisfy $`\{\gamma_j,\gamma_k\}=2\delta_{jk}I`$ and
generate the full matrix algebra on the bank and W. This fixed physical
Majorana frame is used at every native boundary.

Consider an actual emitted native word

```math
U=C_dB_dC_{d-1}\cdots B_1C_0,                            \qquad(3)
```

with the following properties.

* Each B is one physical T/T-dagger band: a product of T or T-dagger
  gates on distinct physical qubits. Address T gates are allowed.
* Each C is an actual global Clifford that preserves every computational
  address sector pointwise, so
  $`C=\sum_x|x\rangle\langle x|\otimes C_x`$.
* For every x, conjugation by $`C_x`$ preserves the real linear span
  of the Majoranas (2).

The restriction concerns the complete Clifford interlayer between
successive native T bands, including the two endpoint blocks. Its
elementary Clifford decomposition need not preserve the span at every
individual gate. Conversely, merely having a Gaussian endpoint U is
insufficient: every actual interlayer in (3) must satisfy the hypothesis.
There is no charge-free change to a different spectral or logical frame.

Assume that a single global phase phi gives

```math
\|U-e^{i\phi}(D_b\otimes I_W)\|_{\rm op}\le\epsilon\lt 1,
\qquad a_{\rm clean}=0.                                \qquad(4)
```

The comparison includes all address, bank and work inputs and their
arbitrary external references. There is no initialized occupation,
one-excitation, fixed-parity or clean-helper promise. Identity on W is
part of the target, so (4) includes its complete return error. Work may
mix arbitrarily with the bank at intermediate steps within the stated
grammar; intermediate return is not assumed. The phase in (4) is common
to the full operator, not a freely chosen phase for each input sector.

**Theorem.** Every word (3) satisfying (4) obeys

```math
2^d\ge 2^b(1-\epsilon^2/2)^2,
\qquad
d\ge\left\lceil b+2\log_2(1-\epsilon^2/2)\right\rceil.
                                                               \qquad(5)
```

In particular,

```math
\begin{aligned}
\epsilon\lt \sqrt{2-\sqrt2}&\quad\Longrightarrow\quad d\ge b,\\
\epsilon\lt 1&\quad\Longrightarrow\quad d\ge b-1.
\end{aligned}                                                   \qquad(6)
```

This d counts actual native T bands of the entire word, not logical
rounds or calls to a separately implemented transform. There is no bound
on the number of dirty helpers or on the intervening Clifford count.

The bit indicator is $`I_b=H_YD_bH_Y`$, where H is applied to
every bank qubit. The same theorem applies to an implementation of
$`I_b`$ when its actual $`H_Y`$-conjugated word for $`D_b`$
has the boundary grammar (3). These endpoint Clifford gates do not change
T-depth, but their preservation of the fixed Majorana frame must not be
assumed automatically.

## 2. Actual Clifford boundaries give address characters

Fix an interlayer C and one work Majorana $`\gamma_j`$.
Global Clifford conjugation sends it to a single Hermitian global Pauli.
Since C preserves each address sector, it commutes with every address Z.
The conjugated Pauli therefore also commutes with every address Z and
contains only I or Z on each address wire. Consequently

```math
C\gamma_jC^\dagger
=(-1)^{c_j}Z_x^{a_j}\otimes P_j                     \qquad(7)
```

for a fixed work Pauli $`P_j`$ and fixed binary data a_j,c_j. On
each address sector, (7) is
$`(-1)^{a_j\cdot x+c_j}P_j`$. The span-preservation hypothesis
forces $`P_j`$ to be one of the standard Majoranas up to sign:
distinct Pauli words are linearly independent, so a single Pauli in this
span cannot be a nontrivial linear combination of its basis elements.
Absorb its sign into c_j. Orthogonality of conjugation makes the resulting
assignment a permutation. Thus

```math
C_x\gamma_jC_x^\dagger
=(-1)^{a_j\cdot x+c_j}\gamma_{\pi(j)},               \qquad(8)
```

where pi is independent of x. In particular, the matrix coefficients of
each Clifford interlayer are either zero or one Walsh character.
An x-dependent permutation is not available from this actual global
Clifford. Scalar phases of $`C_x`$ cancel from (8).

A physical work T or T-dagger gate rotates the two Majoranas belonging
to that physical qubit by angle $`\pm\pi/4`$ and leaves the other
Majoranas unchanged. A whole physical band therefore consists of
disjoint two-coordinate rotations, with at most two nonzero entries in
each row and column. Address T gates are scalar on each x sector and
have no effect on this adjoint action.

Let $`O_U(x)`$ be the real orthogonal matrix of the resulting
Majorana action. Begin with the signed-permutation character matrix for
$`C_0`$. At a T band, each propagation path has at most two
continuations. At a Clifford boundary it has one continuation and gains
one character factor. A product of characters is again a character.
After d bands, every individual matrix coefficient is therefore a sum
of at most $`2^d`$ character terms. Cancellations can reduce its
support but cannot increase it. With normalized Walsh coefficients,

```math
\left|\mathrm{supp}\widehat{(O_U)_{ij}}\right|
\le K,\qquad K=2^d.                                  \qquad(9)
```

This is a path count for each entry, not a count of all paths in the
whole matrix. It is independent of Q and includes all intermediate
bank/helper mixing.

## 3. Full operator norm controls the entire Gaussian matrix

We need a bound independent of Q. Applying a separate adjoint norm
estimate to every Majorana would instead lose a factor proportional to
the width. The full operator contract supplies a stronger statement.

Since (3) preserves x, (4) holds on each branch. Define

```math
V_x=e^{-i\phi}(Z_{Y_x}\otimes I_W)^\dagger U_x,
\qquad \|V_x-I\|\le\epsilon,
\qquad \beta=2\arcsin(\epsilon/2)\lt \pi/3.              \qquad(10)
```

Both factors in $`V_x`$ preserve the Majorana span. Let O be its
orthogonal Majorana matrix. Write Pi for the work parity operator, the
Hermitian unitary proportional to the product of all $`2Q`$
Majoranas. Conjugation gives
$`V_x\Pi V_x^\dagger=\det(O)\Pi`$. Negative determinant
would imply

```math
2=\|V_x\Pi V_x^\dagger-\Pi\|
\le2\|V_x-I\|\le2\epsilon\lt 2,
```

a contradiction. Hence $`O\in SO(2Q)`$.

### Full-spectrum lemma

Choose a real orthogonal normal form for O with Q planar rotation
angles $`\theta_\ell\in[0,\pi]`$. This is an analytic choice
for the proof, not a circuit operation. The rotated Majoranas eta still
satisfy the same anticommutation relations. The Q commuting involutions
$`-i\eta_{2\ell-1}\eta_{2\ell}`$ have every independent
choice of signs as a joint eigenvalue. One way to see this is to expand
their joint spectral projectors: every nonempty Majorana monomial has
zero trace, so each projector has trace and rank one.

The product of the corresponding planar Gaussian rotations implements
O and has spectrum

```math
\exp\left(\frac i2\sum_{\ell=1}^Q s_\ell\theta_\ell\right),
\qquad (s_1,\ldots,s_Q)\in\{\pm1\}^Q.
```

Any two unitaries implementing the same action on all Majoranas differ
by a scalar: their quotient commutes with the full matrix algebra
generated by (2). Thus, for some scalar lambda, the complete spectrum
of $`V_x`$ is

```math
\lambda\exp\left(\frac i2\sum_{\ell=1}^Q
s_\ell\theta_\ell\right),
\qquad (s_1,\ldots,s_Q)\in\{\pm1\}^Q.                \qquad(11)
```

Every sign choice is present because the operator norm in (10) is on
the full arbitrary-work Hilbert space. Restricting to a fixed fermion
parity sector would remove sign choices and would not justify the
following argument.

### The spectral arc cannot wrap

For a unitary within epsilon of I, every eigenphase is in
$`[-\beta,\beta]`$. Compare the all-minus eigenvalue in (11)
with the eigenvalue obtained by changing any subset of signs. Their
ratio shows that every positive subset sum of the angles, modulo
$`2\pi`$, belongs to $`[-2\beta,2\beta]`$.

In particular, each $`\theta_\ell\le2\beta`$: its canonical
range is $`[0,\pi]`$, whereas the other allowed circular arc starts
at $`2\pi-2\beta\gt \pi`$. If the total angle exceeded
$`2\beta`$, its first partial sum crossing that value would lie in

```math
(2\beta,4\beta]\subset(2\beta,2\pi-2\beta),
```

where the inclusion uses $`\beta\lt \pi/3`$. That partial sum is a
forbidden subset sum. Therefore

```math
\sum_{\ell=1}^Q\theta_\ell\le2\beta.                \qquad(12)
```

The scalar lambda in (11) need not be known or independent of x. It
cancels in spectral ratios. This does not replace the common-phase
operator contract (4) by a sectorwise gauge convention.

### Frobenius estimate

The planar normal form and (12) give

```math
\|O-I\|_F^2=8\sum_{\ell=1}^Q\sin^2(\theta_\ell/2)
\le8\sin^2\beta.                                     \qquad(13)
```

For the last inequality, put $`u_\ell=\theta_\ell/2`$ and use
$`\sum_\ell u_\ell\le\beta\lt \pi/2`$. Repeatedly combine
terms using

```math
\sin^2(a+b)-\sin^2a-\sin^2b
=2\sin a\sin b\cos(a+b)\ge0
\quad\text{when }a,b\ge0,\ a+b\le\pi/2.
```

Orthogonal invariance of Frobenius norm converts (13) to the uniform
full-matrix estimate

```math
\|O_U(x)-O_D(x)\|_F^2\le8\sin^2\beta
\qquad\text{for every }x.                            \qquad(14)
```

Its right-hand side does not depend on the number of dirty helpers.

## 4. The selected bank response has full Walsh support

Conjugation by $`Z_{Y_u}`$ changes the signs of precisely the two
standard Majoranas belonging to that physical bank qubit. Thus the two
corresponding diagonal entries of $`O_D(x)`$ both equal

```math
h_u(x)=1-2[x=u].                                      \qquad(15)
```

For $`\chi_a(x)=(-1)^{a\cdot x}`$, use the normalized Fourier
transform $`\widehat f(a)=\mathbb E_x f(x)\chi_a(x)`$.
The coefficients of (15) are

```math
\widehat h_u(0)=1-2/m,\qquad
\widehat h_u(a)=-(2/m)(-1)^{a\cdot u}\quad(a\ne0).
                                                               \qquad(16)
```

Since $`b\ge2`$, all m coefficients have modulus at least
$`2/m`$. If an actual diagonal entry f has support at most K,
at least $`(m-K)_+`$ of these coefficients are absent from f.
Parseval then implies

```math
\mathbb E_x|f(x)-h_u(x)|^2
\ge\frac{4(m-K)_+}{m^2}.                              \qquad(17)
```

There are two such target diagonal entries for each of the m bank
qubits. Summing only these $`2m`$ entries, and then using (14),
gives

```math
8(1-K/m)_+
\le\mathbb E_x\|O_U(x)-O_D(x)\|_F^2
\le8\sin^2\beta.                                     \qquad(18)
```

If $`K\ge m`$, (5) already holds. Otherwise (18) gives

```math
K\ge m\cos^2\beta=m(1-\epsilon^2/2)^2,
```

proving (5). If $`\epsilon\lt \sqrt{2-\sqrt2}`$, the last
factor is strictly larger than $`1/2`$, so the integer band count
satisfies $`d\ge b`$. If $`\epsilon\lt 1`$, it is strictly
larger than $`1/4`$, so $`d\ge b-1`$. The restriction
$`b\ge2`$ is essential to this Fourier argument: for $`b=1`$,
the constant coefficient in (16) vanishes and each target sign function
is a single character.

## 5. What is closed and what remains unpaid

The earlier [robust Gaussian source limit](APPROXIMATE_QUERY_FRAME_TRANSFER.md)
requires a separately implemented Majorana-preserving conjugator. The
[bank-motion construction](BANK_MOTION_WALSH_INDICATOR.md) pays a complete
indicator and shows why a direct fixed transform retains linear depth
in b. The theorem here removes the separate-conjugator premise: any
whole-word synthesis satisfying (3) retains that linear native-depth
cost at sufficiently small constant full-operator error.

Actual inverse words are included. Reversing (3) preserves the same
boundary conditions and number of T bands, and (4) transfers to the
actual adjoint with conjugate global phase because the target is an
involution. All private work remains arbitrary, and its return is
already included in the full norm. No projector, initialized fermion
state, ideal inverse or free spectral basis change is used.

This does not bound all circuits with Gaussian final branches, all
Clifford+T implementations of the indicator, or unrestricted native
T-depth. Non-Gaussian actual Clifford boundaries and address motion
across native T bands leave the hypotheses. They are alternative
mechanisms requiring a new paid construction, not operations ruled out
by (5).

The first unsupported extension to the unrestricted model is concrete.
For the first two work qubits, the allowed native Clifford
$`\mathrm{CX}_{0,1}`$ sends $`\gamma_0=X_0`$ to
$`X_0X_1`$, outside the standard Majorana span (2). Thus arbitrary
Clifford interlayers do not satisfy (8). A composite Clifford block can
nevertheless preserve the span even if some elementary gates in its
decomposition do not. The theorem tests the complete interlayer between
T bands, not each constituent gate.

In particular, the paid phase-indicator construction in the
[bank-motion proof](BANK_MOTION_WALSH_INDICATOR.md#2-native-butterfly-stages-including-their-phases)
does fit the grammar. Its fermionic routing permutations preserve the
Majoranas, and its address-controlled bank Z masks have precisely the
character action (8). For an adjacent pair, the alignment
$`A=S_j\mathrm{CX}_{i,j}H_i\mathrm{CX}_{i,j}`$ acts on the
four local Majoranas, with their common preceding Z string included, as

```math
(\gamma_{2i},\gamma_{2i+1},\gamma_{2j},\gamma_{2j+1})
\longmapsto
(\gamma_{2j+1},-\gamma_{2i+1},\gamma_{2i},-\gamma_{2j}),
\qquad j=i+1.
```

It fixes the other Majoranas. Hence the complete alignment, its actual
inverse, and all resulting Clifford interlayers preserve the required
span, even though the displayed elementary CNOTs and Hadamard need not
do so individually. The bound applies after removing the outer bank
Hadamards to compare with the phase indicator (1).

The [source-fusion recurrence](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
still accepts polynomial count and width overhead but needs a complete
source with $`d_\delta(b)=o(\log(b+2))`$ at each fixed requested
accuracy. The present representation instead costs $`\Omega(b)`$
even with arbitrary finite dirty width. It cannot improve that full
resource frontier through fusion alone. R96 remains unchanged; the
focused original physical-pair lower scout remains a separate question.

The result is analytic and holds at all sizes $`b\ge2`$. No new
finite fixture is needed for the proof, and no empirical screen is used
as an asymptotic claim.
