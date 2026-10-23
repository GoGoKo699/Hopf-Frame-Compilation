# A paid constant-depth approximate zero reflection

[Bounded-clean query promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Returned dirty conjunctions](../../docs/DIRTY_SUM_COMPRESSION.md) ·
[Parallel promotion and its native-sharing requirement](PARALLEL_CLEAN_PROMOTION.md)

For every fixed tolerance $`0\lt \tau\lt 1`$, an approximate reflection about
the all-zero string has a literal zero-clean Clifford+T implementation
with

```math
T,G,b_{\mathrm{dirty}}=O_\tau(p),\qquad
D_T=O_\tau(1).
```

The error is in full operator norm, including every arbitrary borrowed
input and reference. The construction uses a four-clean complete query,
a fixed-point compression filter, and the retained bounded-clean
promotion. All inverses are actual native reversals. The constants are
large because the explicit expander gap is small.

This is an unconditional single-predicate component. Its native
boundaries move the tested register, so it supplies no simultaneous
schedule for the exponentially many rows of a complete indicator.
The [fermionic lift assessment](BANK_MOTION_WALSH_INDICATOR.md#expander-fermionic-lift)
gives an explicit constant-auxiliary lift, its exact full-occupation
spectral error contract, a fixed-branching obstruction, and the
still-unpaid coherent centering operation.

## 1. The supplied CNOT expander and its representation

Anshu, Balasubramanian, Haferkamp, Harrow and Tan,
[Depth-1 expanders on the unitary group and applications](https://arxiv.org/html/2609.01605v1),
Theorem 2.3 and Section 2.2, supply the following group-expander input.
The theorem applies to every finite-dimensional unitary representation,
which is stronger than a gap only in the conjugation representation.

For $`s\ge1`$ and $`n=3s`$, there are fourteen involutions
$`C_1,\ldots,C_{14}\in\mathrm{GL}(n,2)`$ such that

```math
\begin{aligned}
\left\|\frac1{14}\sum_{i=1}^{14}\rho(C_i)
       -\Pi_{\mathrm{inv}}\right\|&\le1-\delta,\\
\delta&=\frac1{48\cdot10^5}.
\end{aligned}                                                   \qquad(1)
```

Here $`\Pi_{\mathrm{inv}}`$ projects onto the invariant vectors of
the representation. Four generators are individual CNOTs. The other
ten are products of s disjoint CNOTs, with disjoint controls and
targets within each product. The displayed generating set comprises
four elementary CNOTs, six within-triple matchings and four cyclically
shifted matchings. Its total CNOT count is $`4+10s`$; coincident
generators at small s cause no difficulty.

Use the computational-basis permutation representation on n qubits.
Its two orbits are the zero string and all nonzero strings. Thus

```math
\begin{aligned}
\mathrm{ran}(\Pi_{\mathrm{inv}})
  &=\mathrm{span}\{|0^n\rangle,|u\rangle\},\\
|u\rangle&=\frac1{\sqrt{2^n-1}}
                 \sum_{x\ne0}|x\rangle.
\end{aligned}                                                   \qquad(2)
```

These vectors identify mathematical subspaces; neither is supplied as
an initialized physical register.

## 2. Selecting zero and obtaining a contracting compression

Put

```math
H_0=I-\frac1{14}\sum_{i=1}^{14}C_i,\qquad
P_1=\frac{I-Z_1}{2}.
```

Both operators annihilate $`|0^n\rangle`$, and
$`H_0\ge\delta(I-\Pi_{\mathrm{inv}})`$. Moreover

```math
a=\langle u|P_1|u\rangle
  =\frac{2^{n-1}}{2^n-1}\ge\frac12.
```

On the orthogonal complement of zero, the nontrivial two-dimensional
block of $`\delta(I-|u\rangle\langle u|)+P_1`$ has trace
$`1+\delta`$ and determinant $`\delta a`$. Every other block has
gap at least $`\min(\delta,1)`$. The smaller eigenvalue is at least
the determinant divided by the trace, giving

```math
H_0+P_1\ge
 \frac{\delta}{2(1+\delta)}
 (I-|0^n\rangle\langle0^n|).                                  \qquad(3)
```

For a directly paid compression, use four clean selector bits. SELECT
applies, uniformly over its sixteen labels, the fourteen CNOT
generators, $`Z_1`$, and I. Write

```math
\begin{aligned}
S_0&=H^{\otimes4}\,\mathrm{SELECT}\,H^{\otimes4},\\
P&=|0000\rangle\langle0000|\otimes I,\\
A&=PS_0P\big|_{\mathrm{ran}(P)}\\
 &=\frac{\sum_i C_i+Z_1+I}{16}
   =I-\frac78H_0-\frac18P_1.
\end{aligned}                                                   \qquad(4)
```

The compression fixes zero. Applying the same two-dimensional argument
to $`7\delta(I-|u\rangle\langle u|)+P_1`$ bounds its largest
eigenvalue on zero-perpendicular by $`1-\gamma`$, where

```math
\gamma=\frac{7\delta}{16(1+7\delta)}.
```

Its smallest eigenvalue is at least $`-7/8`$: each of the fifteen
nonidentity terms in the numerator of (4) is bounded below by -I.
Since $`\gamma\lt 1/8`$,

```math
A|0^n\rangle=|0^n\rangle,\qquad
\|A|_{\{0^n\}^{\perp}}\|\le c:=1-\gamma\lt 1.                  \qquad(5)
```

The emitted $`S_0`$ is an exact Hermitian unitary. It acts as identity
on the data-zero subspace for every selector state and every arbitrary
implementation helper, because each selected generator fixes zero.

## 3. Fixed-point cubing with the same four clean bits

The phase choice is the $`\pi/3`$ fixed-point mechanism of Grover,
[A different kind of quantum search](https://arxiv.org/pdf/quant-ph/0503205),
Sections 2--4 and equation (3). The following block calculation and
complete-return implementation specify its use here.

Let $`P_{\mathrm{sel}}=|0000\rangle\langle0000|`$ and use the
determinant-one selector rotation

```math
\begin{aligned}
R&=\exp\!\left[\frac{i\pi}{3}
       \left(P_{\mathrm{sel}}-\frac I{16}\right)\right]\\
 &=\prod_{\varnothing\ne B\subseteq[4]}
        \exp(i\pi Z_B/48).
\end{aligned}                                                   \qquad(6)
```

Its relative phase between P and its complement is
$`\omega=e^{i\pi/3}`$. Its common scalar is retained and cancels
against the actual inverse in

```math
S_{j+1}=S_j R S_j^\dagger R^\dagger S_j.                       \qquad(7)
```

For any unitary S, put $`A_S=PSP`$ on the range of P. Replacing R
by $`I+(\omega-1)P`$ does not change (7), because the removed scalar
and its inverse cancel. Direct multiplication gives

```math
\begin{aligned}
P(SRS^\dagger R^\dagger S)P
 &=(\omega+\bar\omega-1)A_S\\
 &\quad+|\omega-1|^2 A_S A_S^\dagger A_S\\
 &=A_S A_S^\dagger A_S.
\end{aligned}                                                   \qquad(8)
```

This identity holds for the entire compressed block; no decomposition
into supplied two-dimensional signal circuits is assumed. Since the
initial A is Hermitian, after j steps, with $`N=3^j`$,

```math
PS_jP\big|_{\mathrm{ran}(P)}=A^N.
```

The word contains N occurrences of $`S_0`$ or its actual inverse and
$`N-1`$ occurrences of R or its actual inverse. It still uses only
four clean selectors. On the data-zero subspace it remains exactly
identity, even if R is replaced by a native approximant and every
inverse occurrence uses that approximant's actual inverse.

To pay (6), synthesize its fifteen displayed single-qubit Pauli
rotations to full-norm error at most $`\nu/15`$ each. Computing and
uncomputing their parities uses 34 CNOTs in total. No extra register is
required. If $`t_r,g_r,d_r`$ bound one such single-qubit synthesis,
safe common forward/inverse costs are

```math
T_R=15t_r,\qquad
G_R=15g_r+34,\qquad d_R\le15d_r.                              \qquad(9)
```

The finite native words depend on the chosen tolerance. No fine-angle
gate or unimplemented global phase is treated as free.

## 4. A complete XOR query with four clean selectors

Let u be an arbitrary output bit. Let $`C_4`$ be the exact returned
query toggling u precisely when the four selectors are zero. Emit

```math
V=S_j^\dagger C_4S_j.                                         \qquad(10)
```

Let J initialize only the four selectors and let $`Q_0`$ toggle u
iff the data x is zero. On the data-zero subspace (10) is exact. On
its orthogonal complement,

```math
\|(V-I)J\|\le2\|PS_jJ\|\le2c^N.
```

The data-zero subspace and its orthogonal complement are invariant
throughout. Hence this is a complete initialized-isometry bound for
every data, output, dirty-work and reference input, with no measurement
or postselection. It needs neither a controlled $`S_j`$ nor a
controlled source.

If the emitted rotation approximates R within $`\nu`$, telescoping
its $`N-1`$ occurrences and using actual inverses yields

```math
\|\widetilde VJ-JQ_0\|
 \le\epsilon_V:=2c^N+2(N-1)\nu.                              \qquad(11)
```

The bound includes selector return. The implementation's dirty helpers
return exactly at this stage: $`S_0`$ and $`C_4`$ have exact
whole-space dirty-return identities, and R acts only on the selectors.
The selectors themselves return only approximately on general inputs.

## 5. Literal native implementation and source costs

Use the [returned private-phase Toffoli](../../docs/DIRTY_SUM_COMPRESSION.md#read-only-affine-controls-with-private-phase-work):
eight T/T-dagger gates, at most twelve elementary Cliffords, T-depth at
most eight, and one returned private dirty phase bit. Its controls are
never physical T targets. The two-ladder construction for a q-control
XOR, $`q\ge3`$, uses $`4q-8`$ such Toffolis, $`q-2`$ dirty
chain bits and one reused phase bit.

A selected CNOT has four selector literals plus its original CNOT
control. Its five-control XOR therefore has bounds

```math
T\le96,\qquad G\le144,\qquad d\le96,
\qquad b_{\mathrm{dirty}}\le4.
```

Within one matching, these words run in the same 96 T layers. Their
selector controls are shared only through Clifford controls, and all
physical T targets are private. The selected $`Z_1`$ is a five-wire
phase, with bounds $`T\le64`$, $`G\le98`$, $`d\le64`$, and
three private dirty helpers.

Apply each branch's negative-selector masks once around that branch.
The fifteen nontrivial SELECT branches need at most 480 elementary
Cliffords for these masks. Including both four-Hadamard envelopes,
safe common forward/actual-inverse bounds for $`S_0`$ are

```math
\begin{aligned}
T_0&\le96(4+10s)+64,\\
G_0&\le144(4+10s)+586,\\
d_0&\le14\cdot96+64=1408,\\
b_0&\le4s.
\end{aligned}                                                   \qquad(12)
```

For the negative four-literal query $`C_4`$, including its masks,
take

```math
T_4=64,\qquad G_4=128,\qquad d_4=64,\qquad b_4=3.
```

The dirty pool can be reused after each complete return, so
$`b_V\le4s`$. The source and its actual inverse obey

```math
\begin{aligned}
T_V&\le2NT_0+2(N-1)T_R+64,\\
G_V&\le2NG_0+2(N-1)G_R+128,\\
d_V&\le2Nd_0+2(N-1)d_R+64.
\end{aligned}                                                   \qquad(13)
```

The same emitted word satisfies (11)--(13). Specified data, output and
four source-clean selectors are excluded from $`b_V`$. There is no
unpaid fanout or simultaneous reuse of a physical T target.

## 6. Removing the clean selectors and extracting the phase

Apply [bounded-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
with four clean bits and one specified query output. Its
[approximate block-extension argument](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#4-approximate-clean-return-still-gives-a-full-operator-query)
allows the source's arbitrary off-code action and does not require
address-preserving native boundaries. It emits a zero-clean query
$`\widetilde Q`$ with full-operator error at most
$`128\epsilon_V`$ and common bounds

```math
\begin{aligned}
T_Q&\le64T_V+64T_4+448,\\
G_Q&\le64G_V+64G_4+898,\\
d_Q&\le64d_V+64d_4+256,\\
b_Q&\le6+b_V+b_4\le4s+9.
\end{aligned}                                                   \qquad(14)
```

The former clean selectors are now arbitrary borrowed bits. Their
return and every other helper return are included in the full norm
estimate. The specified data and specified query output are excluded
from $`b_Q`$.

Borrow an arbitrary echo bit u and emit

```math
\widetilde R_n=\widetilde Q^\dagger Z_u\widetilde QZ_u.
```

Its ideal full-space operator, including the literal phase, is

```math
R_n=I-2|0^n\rangle\langle0^n|
```

tensored with identity on u and all other helpers. Thus

```math
\|\widetilde R_n-R_n\otimes I\|
 \le256\epsilon_V.                                           \qquad(15)
```

The T-count and depth are at most $`2T_Q`$ and $`2d_Q`$; the
Clifford count is at most $`2G_Q+4`$; and the borrowed width is at
most $`4s+10`$.

For fixed $`\tau`$, choose a power of three $`N>1`$ such that

```math
2c^N\le\frac{\tau}{2048},\qquad
\nu\le\frac{\tau}{4096(N-1)}.                               \qquad(16)
```

Then $`\epsilon_V\le\tau/1024`$, which also covers the padding
below. N and all rotation-synthesis tolerances depend only on
$`\tau`$ and the fixed gap $`\delta`$, never on the data length.
The clean count stays four throughout this precision amplification.

## 7. Arbitrary lengths and the complete zero-clean ledger

For $`p\ge1`$, set $`n=3\lceil p/3\rceil=p+r`$, where
$`r\in\{0,1,2\}`$, and borrow r arbitrary padding bits z.
Let $`X_b`$ add b only to z. The exact full-space identity is

```math
\prod_{b\in\mathbb F_2^r}X_bR_nX_b
 =(I-2|0^p\rangle\langle0^p|)\otimes I_z.                   \qquad(17)
```

For every padding label, exactly one factor is negative when x is
zero, and none is negative otherwise. This proves (17) for arbitrary
padding superpositions and reference correlations, without a padding
promise. There are $`h=2^r\le4`$ calls. Replacing them by actual
approximate calls costs at most $`256h\epsilon_V\le\tau`$ in
full operator norm. Approximate padding return is included in this
bound, rather than claimed exact.

A cyclic Gray order uses at most h padding X gates. With $`n=3s`$,
one simultaneous final ledger is

```math
\begin{aligned}
T&\le h(128T_V+128T_4+896),\\
G&\le h(128G_V+128G_4+1800)+4h,\\
d&\le h(128d_V+128d_4+512),\\
b_{\mathrm{dirty}}&\le r+4s+10.
\end{aligned}                                                   \qquad(18)
```

When $`r=0`$, omit the final $`4h`$ Clifford allowance. All actual
inverse bounds are the same. Conjugating by X on all p data inputs
changes zero to all ones, adding $`8p`$ elementary Cliffords and no
T gates. Equations (11)--(18) prove the stated
$`T,G,b_{\mathrm{dirty}}=O_\tau(p)`$ and $`d=O_\tau(1)`$
bounds for one complete approximate zero phase, or all-ones
multi-controlled Z.

## 8. Uniform precision and complete-indicator scope

The same ledger also makes the accuracy dependence explicit. Put
$`L=1+\log_2(1/\tau)`$. The gap $`\delta`$ is an absolute
constant, so (16) permits $`N=O(L)`$ and
$`\log_2(1/\nu)=O(L)`$. The retained determinant-one
single-qubit word-length bound
([F5 in the source catalogue](../../docs/reference/SOURCE_CATALOGUE.md),
GKW Lemma 2.3) therefore gives
$`t_r,g_r,d_r=O(L)`$, without extra workspace. Substitution in
(9), (12)--(18) yields

```math
\begin{aligned}
T,G&=O(pL+L^2),\\
b_{\mathrm{dirty}}&=O(p),\qquad d=O(L^2).
\end{aligned}                                                   \qquad(19)
```

The constants are absolute and include the stated small expander gap.
These are simultaneous bounds for the same emitted full-norm word
and its actual inverse, with no initialized work. The precision-depth
term is quadratic in L; it is not a linear-precision depth claim.

The single-predicate primitive is paid, including its precision,
arbitrary-helper return and native inverse. Its CNOT generators move
the tested data between native boundaries. Independent row copies
therefore cannot share that data by the native parallel-source
schedule. Logical commutation of the completed row reflections does
not change this dependency.

Serial constant-accuracy row calls also lack a dimension-independent
full-bank error bound. A complete indicator still requires a paid
simultaneous bank schedule and its whole-operator error analysis. This
proof establishes neither and does not change the general
complete-frame upper frontier.

## 9. The uniform-address obstruction to joint relabeling

A direct whole-bank extension would let the expander act jointly on
the address and on the labels of the arbitrary output bank. This
extension fails even after arbitrary address-independent quantum
operations on the bank, selectors and helpers. The obstruction has
no gate-count, depth or helper-width hypothesis.

Let the address set have size $`m=2^p`$, with $`p\ge1`$. A
permutation g acts jointly by

```math
U_g|x,Y\rangle=|gx,gY\rangle,\qquad
(gY)_{gx}=Y_x.
```

Allow these joint permutations, including selector- or helper-controlled
versions, and all unitaries that act only on the bank and helper
registers. The latter class includes arbitrary bank-changing and
nonmonomial operations, every bank-diagonal phase, all fixed bank Z
probes, and arbitrary quantum selector dynamics. Every allowed
operation and its inverse commutes with

```math
\begin{aligned}
P_+&=|+_{\mathcal D}\rangle\langle+_{\mathcal D}|
       \otimes I_{Y,\mathrm{helpers}},\\
|+_{\mathcal D}\rangle&=\frac1{\sqrt m}\sum_x|x\rangle.
\end{aligned}                                                   \qquad(20)
```

For joint relabeling this follows because every permutation fixes the
uniform address vector; its simultaneous bank permutation does not
change the identity factor. Address-independent operations commute
with (20) directly. The same invariant therefore survives any complete
word, arbitrary common helper preparation, and final helper
compression or readout.

The target is the complete addressed phase

```math
\Omega|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle.
```

Choose any balanced computational bank Y, with exactly m/2 ones.
Then

```math
\langle+_{\mathcal D},Y|
 \Omega|+_{\mathcal D},Y\rangle
 =\frac1m\sum_x(-1)^{Y_x}=0.                                  \qquad(21)
```

Thus $`\Omega|+_{\mathcal D},Y\rangle`$ lies in the orthogonal
complement of the uniform-address subspace. The input is a legitimate
state of the requested arbitrary-bank domain, not an extra promise.

Let K be a compression of an allowed word, using any common
address-independent helper insertion and readout. Its output on this
input remains in the uniform-address subspace. For every complex
scalar $`\alpha`$, orthogonality gives

```math
\|K-\alpha\Omega\|\ge|\alpha|.                             \qquad(22)
```

In particular the amplitude-$`1/\sqrt2`$ contract cannot be
approximated below absolute error $`1/\sqrt2`$ by this grammar.

For a complete initialized-unitary implementation with helpers
returned, both the actual and desired output vectors are normalized
and orthogonal. Consequently

```math
\|UJ-J\Omega\|\ge\sqrt2,                                   \qquad(23)
```

where J inserts the stipulated helper state. Leakage from that state,
arbitrary quantum helper mixing, a common external scalar phase and
reference systems do not alter the bound: actual outputs still have
uniform address, while the required output has orthogonal address.

A word $`J_0^\dagger Z_{Y_0}J_0`$ also preserves (20) whenever
$`J_0`$ does, so the same $`\sqrt2`$ gap applies directly to its
complete addressed-observable transport. Moreover the bank Hadamard
wall commutes with (20) and conjugates the complete XOR indicator to
$`\Omega`$. The obstruction therefore covers that query contract as
well. A lookup XOR into an extra arbitrary bit similarly reduces to
the phase contract on input $`|{-}\rangle`$.

The representation-wide Cayley gap used for the single zero test does
not remove this common invariant in the joint construction. A new
whole-bank mechanism must include an address-sensitive operation that
breaks it; bank-changing operations alone are insufficient. Such an
operation can be cheap, and this result places no native-depth lower
bound on it. What remains unpaid is a complete word outside this
algebra with a uniform whole-bank error estimate, arbitrary-helper
return and a simultaneous native resource ledger.

## 10. Bounded-choice joint SELECTs need growing product length

Local address phases can leave the preceding invariant subspace, but
they do not by themselves supply the missing full-bank compression.
There is a separate rank obstruction across the cut between the p-bit
address and the bank together with **all** selectors and private helpers.

Use normalized Hilbert--Schmidt inner products on both sides of this
cut. The target has the operator-Schmidt decomposition

```math
\Omega\otimes I_h
 =\frac1{\sqrt m}\sum_{x=0}^{m-1}
   \bigl(\sqrt m\,|x\rangle\langle x|\bigr)
      \otimes\bigl(Z_{Y_x}\otimes I_h\bigr),
 \qquad m=2^p.                                               \qquad(24)
```

Each displayed family is orthonormal. Thus all m Schmidt coefficients
equal $`1/\sqrt m`$, independently of the arbitrary dirty width.
If K has operator-Schmidt rank at most R and
$`\|K-\alpha\Omega\otimes I_h\|_{\mathrm{op}}\le\epsilon`$,
normalized Hilbert--Schmidt error is also at most epsilon. The best
rank-R approximation to (24) discards at least m-R equal squared
coefficients. Consequently, for nonzero alpha and epsilon below its
modulus,

```math
\epsilon^2\ge|\alpha|^2\left(1-\frac Rm\right),
\qquad
R\ge m\left(1-\frac{\epsilon^2}{|\alpha|^2}\right).
                                                               \qquad(25)
```

A completed s-choice joint relabeling SELECT, counting any identity
padding branch among its choices, has the form

```math
\sum_{j=1}^{s}A_j\otimes
       \bigl(|j\rangle\langle j|_{\mathrm{sel}}\otimes B_j\bigr),
                                                               \qquad(26)
```

where $`A_j`$ relabels the address and $`B_j`$ relabels the bank.
It has Schmidt rank at most s, as does its actual inverse. Arbitrary
address-local operations, including Z phases and Hadamards, have rank
one. Every bank/selector/helper-local operation also has rank one,
including bank-changing operations. A fixed-arity interaction across
the cut has a fixed rank bound; a single cross-cut CNOT has rank two.

Schmidt rank is submultiplicative under products. Helper insertion and final
compression act only on the bank side of this cut, so they cannot
increase the rank. Hence any clean compression made from completed
SELECTs of sizes $`s_j`$ and additional cross-cut probes of rank bounds
$`r_\ell`$ satisfies

```math
\prod_j s_j\prod_\ell r_\ell
 \ge m\left(1-\frac{\epsilon^2}{|\alpha|^2}\right)             \qquad(27)
```

whenever it meets the target error. In particular, with a uniform
constant bound greater than one on these ranks, the number of such
completed factors is $`\Omega(p)`$ at fixed relative error below one.
A constant-size LCU followed by a constant number of filter calls
therefore cannot supply the required nonzero-amplitude complete source
merely by adding local address phases or local bank basis changes.

This is a restriction on the displayed product representation, **not**
a native T-depth bound. For example, p CNOTs from the address into p
distinct dirty helper bits have Schmidt rank $`2^p`$ and fit in one
Clifford boundary. Such a wide coupling leaves the bounded-rank grammar
at no T cost. The existing private quadratic-query constructions also
contain wide couplings and are not excluded here. A next source may use
that freedom, but must still prove its full compression, simultaneous
resources, literal phase and arbitrary-helper return.

## 11. A paid graph reflection still needs complete readout

The [nilpotent-pencil query](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#3-the-full-quadratic-xor-query-is-paid)
supplies a wide coupling outside the preceding product-rank restriction.
Put $`w=p2^p`$, $`d=2^w`$ and $`A_x=N_x^{\mathsf T}`$, using
that note's binary transpose convention. On two arbitrary w-bit
principal banks z,t, its exact returned query is
$`Q_A:(z,t)\mapsto(z,t+A_xz)`$. The new zero reflection gives

```math
\begin{aligned}
R_{\Gamma_x}&=Q_A^\dagger
       (I-2I_z\otimes|0^w\rangle\langle0^w|_t)Q_A
       =I-2P_{\Gamma_x},\\
P_{\Gamma_x}&=\sum_z|z,A_xz\rangle\langle z,A_xz|.
\end{aligned}                                                   \qquad(28)
```

This is a complete graph-predicate reflection, with no initialized z
or t. If $`R_t`$ is the emitted zero reflection at full error zeta,
safe common bounds for the emitted graph word are

```math
\begin{aligned}
T_\Gamma&\le2T_Q+T_{R_t},&G_\Gamma&\le2G_Q+G_{R_t},\\
d_\Gamma&\le16+d_{R_t},&
b_\Gamma&\le\max(b_Q,b_{R_t}),\\
\epsilon_\Gamma&\le\zeta.&
\end{aligned}                                                   \qquad(29)
```

The private pool is reused after complete query return. Both w-bit
banks remain principal inputs, excluded from $`b_\Gamma`$. Embedded
around the specified output bank, they contribute
$`(2p-1)2^p`$ additional arbitrary wires. At fixed zeta, (29) has
polynomial overhead in the relaxed source envelopes and constant
native T-depth. Every inverse and return error is included.

### Exact overlaps and their normalization

For one fixed address, write A for the invertible binary matrix
$`A_x`$ and introduce the mathematical code embeddings

```math
J_\Gamma|z\rangle=|z,Az\rangle,\qquad
J_L|z\rangle=|z,+\rangle,\qquad
J_R|t\rangle=|+,t\rangle.
```

Here each plus state is uniform on a whole w-bit bank; these embeddings
are not supplied physical initializations. Direct matrix multiplication
gives

```math
\begin{aligned}
J_\Gamma^\dagger J_L&=d^{-1/2}I,\\
J_\Gamma^\dagger J_R&=d^{-1/2}U_{A^{-1}},\\
J_L^\dagger P_\Gamma J_R&=d^{-1}U_{A^{-1}},
\qquad U_A|z\rangle=|Az\rangle.
\end{aligned}                                                   \qquad(30)
```

All principal cosines between the graph and either plus-bank code are
$`d^{-1/2}`$, independent of A. The usual two-projector rotation
therefore needs order $`\sqrt d`$ calls to obtain constant transfer
amplitude from that small overlap. This describes that amplification
mechanism, not an unrestricted inversion lower bound. The initial
plus-bank code itself would reserve w clean bits. Against the
computational code $`I_z\otimes|0\rangle\langle0|_t`$, the graph
instead has one principal cosine equal to one and all others zero.
Its intersection does not supply the inverse map on arbitrary inputs.

The usual forward Hadamard postselection has exactly the same missing
normalization:

```math
(\langle+|_z\otimes I_t)Q_A
       (I_z\otimes|0\rangle_t)=d^{-1/2}U_A.                 \qquad(31)
```

This block is flat, but at amplitude $`2^{-w/2}`$, with an entire
initialized bank. Paying the reflections does not turn (30) or (31)
into a flat amplitude-$`1/\sqrt2`$ source with a bounded clean bank.

For completeness, the fixed translation-Bell code has basis
$`|\Phi_u\rangle=d^{-1/2}\sum_z|z,z+u\rangle`$. If
$`\rho=\mathrm{rank}(A+I)`$, its nonzero principal cosines with
the graph are $`2^{-\rho/2}`$, with multiplicity $`2^\rho`$;
the remaining $`d-2^\rho`$ are zero. An apparent
$`1/\sqrt2`$ at rho equal to one therefore covers only rank two.
For the nilpotent pencil at nonzero x, $`\rho=w-p`$, so this
comparison also retains small amplitudes and a nontrivial kernel.

### Thin reflections cannot repair a separated full-space baseline

There is a broader but explicitly baseline-dependent check. On a
physical space of dimension D containing the two arbitrary banks,
an ideal graph reflection differs from identity on rank $`D/d`$.
Arbitrary dirty factors preserve this fraction; signal controls can
only decrease it. Consider a unitary word W containing r such calls
and arbitrary unitary interlayers. Let $`W_0`$ be the same word with
those calls replaced by identity. Telescoping the differences gives

```math
\mathrm{rank}(W-W_0)\le rD/d,\qquad \|W-W_0\|_{\mathrm{op}}\le2.
```

Compress any k clean bits, with logical dimension $`D/2^k`$, and
call the resulting operators K and $`K_0`$. In normalized
Hilbert--Schmidt norm,

```math
\|K-K_0\|_{2,\mathrm{norm}}
       \le2\sqrt{r2^k/d}.                                  \qquad(32)
```

The initial and final clean embeddings may differ. No principal bank
has been silently initialized; its w-bit initialization, if chosen,
must count toward k. The same rank statement allows coherent address
motion in the interlayers, since every graph fiber has the same rank.

Suppose each actual graph call has full-space operator error at most
zeta, and the emitted compression approximates a desired complete
operator $`K_*`$ to operator error epsilon. If the erased baseline has
normalized Hilbert--Schmidt gap
$`\Delta=\|K_0-K_*\|_{2,\mathrm{norm}}`$, then

```math
\epsilon\ge\Delta-r\zeta-2\sqrt{r2^k/d}.                    \qquad(33)
```

For observable transport $`W^\dagger P W`$, replace r by $`2r`$
in both correction terms. If completed interlayers preserve x and
the erased fixed-x baseline is Clifford, use the same computational-zero
embedding for the k clean qubits on both sides of the observable compression.
Its compressed Pauli is then either zero or a signed logical Pauli. A wrong label then has gap
$`\sqrt2`$ from the desired unit Pauli, and a zero compression has
gap one. Small numbers of thin graph calls cannot correct such a
sector with a small clean bank and small total approximation error.

No baseline gap is presumed for arbitrary interlayers: they might
already implement the desired inverse-column transport. Thus (33)
does not bound unrestricted inversion or native depth. It identifies
the unpaid operation in this candidate: complete extraction from the
graph, or a different baseline carrying the required information,
rather than the graph predicate, its normalization, or its native
reflection cost alone.

### Coset labels preserve the missing extraction cost

For fixed A, the graph cosets
$`P_c=\sum_z|z,Az+c\rangle\langle z,Az+c|`$ each have
rank d and sum to identity. Adding an arbitrary w-bit label c gives
$`\sum_c|c\rangle\langle c|\otimes P_c`$, of rank
$`d^2`$ in dimension $`d^3`$: the fraction remains $`1/d`$.
An arbitrary label does not promise that it equals the actual syndrome.
Summing all cosets removes the dependence on A. Weighted sums instead give

```math
\begin{aligned}
D_f&=\sum_c f(c)P_c
     =Q_A^\dagger(I\otimes\mathrm{diag}(f))Q_A,\\
D_f|z,t\rangle&=f(t+Az)|z,t\rangle.
\end{aligned}                                                   \qquad(34)
```

A fixed character $`f(c)=(-1)^{\ell^{\mathsf T}c}`$ yields
the forward Pauli $`Z_z(A^{\mathsf T}\ell)Z_t(\ell)`$.
To obtain $`Z_t(A^{-\mathsf T}\ell)`$ after the known correction
$`Z_z(\ell)`$, (34) requires the weights
$`f_x(c)=(-1)^{\ell^{\mathsf T}A_x^{-1}c}`$ already; setting
z to zero proves necessity. The inverse-label readout has moved into
the weights. Likewise, the paid graph-Bell basis change $`Q_AH_z`$
has a reverse-coordinate description with labels
$`A^{-\mathsf T}s,A^{-1}c`$ and scalar
$`(-1)^{s^{\mathsf T}A^{-1}c}`$. That description is not a
free implementation of the reverse labels.

A constant-fraction image under a partial trace also retains its
normalization. For nonzero ell, let $`P_\ell`$ project onto
$`\ell^{\mathsf T}z=0`$ and let $`P_{\rm inv}`$ project onto
$`\ell^{\mathsf T}A^{-1}t=0`$. Then

```math
\begin{aligned}
\mathrm{Tr}_z[(P_\ell\otimes I)P_\Gamma]&=P_{\rm inv},\\
J_R^\dagger(P_\ell\otimes I)P_\Gamma J_R&=d^{-1}P_{\rm inv},\\
J_R^\dagger Z_z(\ell)P_\Gamma J_R&=d^{-1}Z_t(A^{-\mathsf T}\ell).
\end{aligned}                                                   \qquad(35)
```

The first line is an unnormalized partial trace, not a returned
compression. The half-graph has nonzero principal cosines
$`d^{-1/2}`$ against the plus-bank code, with multiplicity d/2.
Adding a balancing arbitrary bit merely replaces the identity
baseline by its Z Pauli; the thin perturbation and (32) remain.

### Constant-fraction syndrome predicates have an affine-label limit

There is a separate obstruction after replacing the thin graph by
a fixed rank-k syndrome constraint. For a fixed binary matrix M of
rank k, the ideal reflection is

```math
R_M=Q_A^\dagger(I-2\Pi_{Mt=0})Q_A
   =I-2^{1-k}\sum_{s\in\mathbb F_2^k}
          Z_z(A_x^{\mathsf T}M^{\mathsf T}s)Z_t(M^{\mathsf T}s).
                                                               \qquad(36)
```

Its Pauli support lies in a subgroup of size $`2^k`$, whose labels
are affine in x because $`A_x`$ is affine. Any phase function of
these same k syndrome bits has the same property. A constant-fraction
predicate therefore avoids the rank fraction in (32), but not the
following support restriction.

Let W use r such syndrome phases, with ranks $`k_1,\ldots,k_r`$,
interspersed with fixed bank/helper Cliffords independent of x.
All completed factors preserve the address. For any fixed Pauli P,
the transported observable $`W^\dagger P W`$ is a sum of at most
$`2^{\sum_i k_i}`$ affine Pauli-label branches; coefficients may
depend on x. To see the count, conjugation of a Pauli by one syndrome
phase stays in that Pauli's coset of the same syndrome subgroup.
The two Fourier indices combine inside the subgroup, so the support
multiplier is $`2^{k_i}`$, not $`2^{2k_i}`$. Fixed Clifford
conjugation preserves affine labels, giving the claim by induction.

Project each label onto the Z-label of the specified m-bit output
bank, where $`m=2^p`$. One affine map $`v(x)`$ can equal
$`e_x`$ at at most p+1 addresses: its evaluation matrix factors
through $`(1,x_1,\ldots,x_p)`$, whereas those matching columns
are distinct standard basis vectors. Thus if

```math
2^{\sum_i k_i}(p+1)\lt 2^p,                                   \qquad(37)
```

some address has no desired $`Z_{Y_x}`$ component. Standard clean
compression can delete or rescale helper Pauli factors but cannot
create that missing bank label. Orthogonality in normalized
Hilbert--Schmidt norm then gives operator error at least
$`|\alpha|`$ to $`\alpha\Omega`$, including arbitrary dirty
helpers and any fixed number of computational-zero clean bits.
If each actual syndrome-phase call has full-space error at most zeta,
the lower bound is $`|\alpha|-2r\zeta`$ for the emitted observable:
it physically uses both W and its actual inverse. The support count
is not doubled, but the call and error ledger is.

For uniform k and $`\epsilon+2r\zeta\lt |\alpha|`$, the bound therefore requires
$`kr\ge p-\log_2(p+1)`$. In particular, constant-k predicates
and $`O(\log\log p)`$ such rounds do not complete this readout.
This is an observable-transport statement; a direct source
compression with an arbitrary net Clifford requires its own support
analysis. Predicates on all w syndrome bits, address motion, or
other couplings can leave the stated grammar. Nor does (37) bound
native T-depth: it counts completed factors in this representation,
which might admit a different joint synthesis. The next mechanism
must pay complete readout while using more than a bounded collection
of these affine-label branches.

<a id="synchronized-anchored-generator-stage"></a>

## 12. Shared address-generator stages and their bank-moment limit

The data motion in the scalar filter motivates a different native
schedule. A single translated generator can in fact be shared across
all rows without copying or initializing the address. The exact
identity below pays that stage. Its simple reuse with the physical
bank as controls, however, retains only p+1 bank parities and cannot
give the complete lookup. Neither fact asserts a general depth bound.

Let C be a product of q CNOTs on p address bits, with all control and
target wires pairwise distinct. Write its pairs as
$`(c_j,t_j)`$, so $`C^2=I`$. For each row label
$`u\in\mathbb F_2^p`$, define the anchored affine permutation

```math
g_u(x)=C(x+u)+u=Cx+(I+C)u.
```

The translations $`(I+C)u`$ are fixed by C. Thus every g_u is an
involution, and all these maps for the same C commute. Give each row
an arbitrary coherent control a_u and put

```math
M_0(a)=\bigoplus_u a_u,\qquad
M_1(a)=\sum_u a_u u\in\mathbb F_2^p.
```

The complete product, in any row order, is

```math
\prod_u g_u^{a_u}(x)
   =C^{M_0(a)}x+(I+C)M_1(a).                         \qquad(38)
```

This identity includes coherent controls: both sides are the same
permutation of every joint computational column, with scalar +1.
It does not replace the controls by sampled classical values.

### Native word and simultaneous ledger

Borrow an arbitrary parity bit h and one private phase helper for
each CNOT pair. First XOR every a_u into h, apply C controlled by h,
undo those parity CNOTs, and apply C controlled by h again. The two
controlled powers multiply to $`C^{M_0(a)}`$, independently of the
old h, which returns. Each controlled matching is a wall of q exact
returned private-phase Toffolis. Their physical T targets are the q
distinct phase helpers; h is only a Clifford control. Each wall has
eight T layers, 8q T gates and 10q elementary Cliffords, including
the target Hadamards.

Finally, for each pair $`(c_j,t_j)`$ and every label with
$`u_{c_j}=1`$, apply CNOT from a_u to the address target t_j. These
$`qm/2`$ CNOTs add the translation in (38), where $`m=2^p`$.
Together with the two m-CNOT parity words, the forward word and its
actual native reverse have the common complete ledger

```math
\begin{aligned}
T&=16q,&D_T&\le16,\\
G&=20q+2m+qm/2,&b_{\rm dirty}&=q+1,\\
a_{\rm clean}&=0,&\varepsilon&=0.
\end{aligned}                                                   \qquad(39)
```

The p address bits and m row controls are principal arbitrary inputs,
outside the private width in (39). Every private helper returns
exactly; arbitrary references tensor identity. In particular this
stage has T-count O(p), Clifford count O(mp) and constant T-depth.
If a_u is instead a nontrivial function of a row selector and a bank
bit, computing that function is a separate operation not included in
(39).

Different generator matchings need not commute. Equation (38) permits
reordering rows only within one fixed generator stage. It does not
justify moving all equal-generator factors together inside a product
of complete row filters. For example, take two address bits,
$`C=\mathrm{CX}(0;1)`$, $`D=\mathrm{CX}(1;0)`$, and active
labels 0 and e_0. Starting at address zero, applying C then D within
each row gives $`e_0+e_1`$; applying both C rows and then both D
rows gives e_1. Nor does the shared stage provide the multi-selector
compression or clean removal for a new synchronized filter word.

<a id="shared-stage-bank-moment-gap"></a>

### Reusing the same physical bank as controls cannot complete lookup

There is an exact limitation on the simplest proposed joint word.
At every completed stage (38), set $`a_u=Y_u`$ using the same
physical m-bit bank Y. Allow any number of stages, any CNOT matchings,
their actual inverses, and arbitrary intervening unitaries on the
address, signals and private work. Those interlayers act as identity
on Y, and all dependence on Y must occur through the displayed stages.
The address can move coherently. Helper insertion and final
compression are common to all Y and act outside Y.

By (38), each stage depends on Y only through $`M_0(Y),M_1(Y)`$.
The complete word, its accepted block K and its Gram therefore have
identical fibers whenever these p+1 parities agree. This remains true
after any number of layers, without a depth or private-width bound.

For p at least two, choose linearly independent row labels a,b. Let
Delta be the bank string supported on the four labels
$`0,a,b,a+b`$. Then

```math
M_0(\Delta)=0,\qquad M_1(\Delta)=0,\qquad \Delta_0=1.
                                                               \qquad(40)
```

Hence Y and $`Y+\Delta`$ have identical complete source fibers, but
the desired lookup signs at address zero are opposite. For
$`\Omega|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle`$, the triangle
inequality on these two identified fibers gives

```math
\|K-\alpha\Omega\|\ge|\alpha|.                       \qquad(41)
```

In particular a complete phase implementation in this grammar has
error at least one. For the half-projector Gram target, the desired
expectations on the same address/work input are one half and zero,
while the two actual Gram expectations are equal. Consequently

```math
\left\|K^\dagger K-\frac{I+\Omega}{4}\right\|
   \ge\frac14.                                         \qquad(42)
```

This is also the paired lookup Gram in its single-bank form after
the fixed bitwise CNOT change of coordinates. Approximation of a
displayed ideal word by full operator error epsilon decreases (41)
by at most epsilon and (42) by at most $`2\epsilon`$.

The premise is the complete bank dependence, not the number of
stages. Independent row selectors multiplying Y_u, other bank-sensitive
interactions, or intervening operations on Y can leave this grammar;
their full action and return must be paid separately. The next missing
operation is a joint source that retains enough bank information and
has bounded-clean complete return. Sharing one generator is now paid,
but alone supplies neither that operation nor a better frame frontier.

<a id="independent-row-selector-echo"></a>

### Independent row selectors: paid motion and a collapsing dirty echo

An explicit changed-input trial assigns an arbitrary selector s_u to
every row and replaces the control a_u by s_u Y_u. For one matching C,
define

```math
S_C(s,Y)=\prod_u g_u^{s_uY_u},\qquad
S_C(s,Y)x=C^{M_0(sY)}x+(I+C)M_1(sY).                    \qquad(43)
```

Here sY means coordinatewise binary multiplication, not a supplied
product bank. The operators preserve s and Y. Equation (43) is a
complete permutation identity with literal phase +1, not a statement
only about classical selector inputs.

There is a direct native realization without initializing the selectors.
Borrow one arbitrary parity bit h and emit chronologically:

1. The returned bilinear query B toggling h by the sum of s_u Y_u.
2. The matching C controlled by h.
3. The actual native inverse of B.
4. The actual native inverse of the controlled matching.
5. For each matching pair (c_j,t_j), toggle address bit t_j by
   the sum of s_u Y_u over labels with u at coordinate c_j equal to one.

The first four operations give the power C to M_0(sY), independently
of the original h, and return h. The last operation gives the translation
in (43). Each bilinear wall is a collection of exact private-phase
Toffolis: enclose their commuting diagonal words by one Hadamard pair
per output, and place all T gates on distinct private helpers. Shared
selector, bank and output wires occur only as Clifford controls inside
the phase words. Thus each of the five walls has at most eight T layers.
For q at least one, safe common forward/actual-inverse bounds are

```math
\begin{aligned}
T_S&\le16m+16q+4qm,&d_S&\le40,\\
G_S&\le20m+20q+5qm,&
b_S&\le1+\max(m,q,qm/2).
\end{aligned}                                                   \qquad(44)
```

The p address bits, m bank bits and m row selectors are principal
arbitrary inputs. All additional helpers return exactly. If borrowed
as source work, the m selectors must also be included in its width.
There is no precision error. Formula (44) fits polynomial intermediate
budgets, but supplies neither a complete indicator nor clean removal.

The natural zero-clean echo removes precisely the information it added.
Let X_s flip every row selector, using its literal four-Clifford X word.
Since all anchored maps for this one C commute and square to identity,

```math
S_C X_s S_C^\dagger X_s
 =\prod_u g_u^{s_uY_u+(s_u+1)Y_u}
 =\prod_u g_u^{Y_u}.                                    \qquad(45)
```

It returns every selector and has bounds 2T_S, 2d_S and
2G_S+8m, reusing the returned private pool. Its full action is exactly
the old stage (38), so any word assembled only from these echoes and
the interlayers allowed above retains the Gram gap one quarter.
Without an echo or internal selector changes, the s=0 fiber is
bank-independent and likewise cannot be a complete lookup.

Even granting m initialized selectors does not repair the unmixed
averaging version. Put Hadamards on all selectors only at entrance and
exit, use any number of stages (43), and allow common interlayers on
address and other work but not on s or Y. For fixed Y, every selector
path depends only on M(sY). Uniform s makes this moment vector uniform
on

```math
L(Y)=\mathrm{span}_{\mathbb F_2}\{(1,u):Y_u=1\}.             \qquad(46)
```

Indeed the linear map from s onto L(Y) has equally sized fibers.
The accepted block therefore depends only on L(Y), even when the
different matching stages do not commute. For p at least two, the full
bank and the full bank missing row zero have the same L: for independent
a,b, the three nonzero labels a,b,a+b already sum to (1,0) after
augmentation. Their lookup signs at address zero are opposite. The
same two-fiber comparison proves Gram gap at least one quarter for
this accepted word. Its m clean inputs are additionally outside the
bounded-clean source budget; the gap holds even if that cost is granted.

These are limits of the displayed dirty echo and unmixed acceptance,
not of general independent-selector words. Internal selector mixing,
bank interactions or different-generator commutators can leave their
premises. Such a candidate still owes its full operator and selector
return; the failed row/stage reordering above cannot supply them.

The proof is analytic. The exact native and bank-collision fixtures in
[the synchronized-stage tests](../../tests/test_synchronized_address_stage.py)
are finite internal checks, not external mathematical review.

<a id="mixed-selector-commutator"></a>

## 13. Internal selector mixing: a paid commutator with unresolved return

The next literal trial leaves the unmixed-selector premise. For two
matchings C,D, take the exact stages (43) and put, in operator order,

```math
A=S_C(s,Y),\qquad B=H_s^{\otimes m}S_D(s,Y)H_s^{\otimes m},
\qquad W=ABA^\dagger B^\dagger.                         \qquad(47)
```

Use the actual native reversals in (47), not a promise that each
emitted word is syntactically its inverse. The logical A,B are
Hermitian. With the common stage bounds (44), W has zero clean work,
error zero, T-count at most 2T_C+2T_D, T-depth at most 160, Clifford
count at most 2G_C+2G_D+4m, and private pool at most
max(b_C,b_D). The m arbitrary selectors are additional borrowed work.
Every private stage helper returns; the selectors do not generally do so.
These are simultaneous polynomial-overhead budgets, not a complete source.

Set the bank to one active row u. All selectors except s_u decouple.
Write g=g_{C,u}, h=g_{D,u} and $`P_h^\pm=(I\pm h)/2`$.
In the s_u computational basis the exact surviving factors are

```math
A=\begin{pmatrix}I&0\\0&g\end{pmatrix},\qquad
B=\begin{pmatrix}P_h^+&P_h^-\\P_h^-&P_h^+\end{pmatrix},
\qquad W_{01}=P_h^-gP_h^+.
                                                               \qquad(48)
```

Consequently $`\|W_{01}\|=\|[g,h]\|/2`$. Indeed the commutator
has off-diagonal blocks twice those of g in the h eigenspaces.
Selector return on this complete one-row sector holds if and only if
g and h commute. Necessity follows from (48) and Hermiticity of g;
if they commute the entire word is

```math
W=I_{s_u}\otimes\bigl(I-2P_g^-P_h^-\bigr).             \qquad(49)
```

For C=D, all anchored row maps commute with one another, so the full
m-selector echo returns every selector and is exactly
$`\prod_u g_{C,u}^{Y_u}`$: it recovers the old moment stage (38).
That returned special case does not retain new lookup information.

There is a dimension-independent leakage witness for different stages.
For p at least two choose C=CNOT(0,1), D=CNOT(1,0), u=0, and fix
the other address bits to zero. The two nontrivial address reflections
have negative vectors
$`|g_-\rangle=(|1\rangle-|3\rangle)/\sqrt2`$ and
$`|h_-\rangle=(|2\rangle-|3\rangle)/\sqrt2`$.
Their product has order three, giving

```math
\|[g,h]\|=\sqrt3,\qquad \|W_{01}\|=\sqrt3/2.          \qquad(50)
```

Thus W is at distance at least sqrt(3)/2 from every operator that
returns s_u and acts only on the address. This is an exact full-space
return failure, not a numerical failure on an initialized selector.

Taking a fixed power does not repair exact return in this fixture.
A and B are rank-one Householder reflections on eight dimensions;
their negative vectors are $`|1\rangle|g_-\rangle`$ and
$`|{-}\rangle|h_-\rangle`$, with overlap magnitude
$`1/(2\sqrt2)`$. W fixes their six-dimensional common complement.
Its other two eigenvalues satisfy $`\lambda+\bar\lambda=1/4`$.
They are not roots of unity: otherwise this rational trace would be
an integer, since it would be an algebraic integer. For every positive
integer k, the range of $`W^k-I`$ is therefore the same two-dimensional
plane. This plane is not a tensor product of the whole selector with
one address line, since g_- and h_- are independent. But a nonidentity
$`I_s\otimes V`$ of rank-two difference from identity has exactly
that range. No positive power W^k returns the selector exactly.
This does not give an approximation lower bound for arbitrary powers;
powers approaching identity are not excluded.

The opposite-CNOT trial also fixes both the marked address and the
uniform vector on its other three addresses. Every word in these
controlled permutations and selector-only interlayers acts identically
on those two address states. A common helper insertion/acceptance
therefore has Gram gap at least one quarter against the lookup
half-projector: the two required expectations differ by one half.
This is the same common-invariant principle as Section 9, now on the
one-row fiber, not a new unrestricted depth bound. A bank scalar
cannot distinguish those states either.

Pause this undressed commutator. A phase-dressed generator or another
address-sensitive insertion may break that invariant, but then its
complete accepted Gram and arbitrary-selector return must be derived
again. The six-wire native fixture checks all 64 columns with arbitrary
bank bit and both private helpers, 128 T gates and 156 Cliffords in
either orientation. A second exact fixture checks (48), the leakage
eigenvalue, the power polynomial and the failure of the three single
selector-Pauli echoes. The latter are finite checks, not an exclusion
of general echo words. The all-size scope above is analytic and
internally audited; it supplies no external mathematical review.

<a id="quarter-phase-selector-stage"></a>

### A paid quarter-phase stage for a changed acceptance question

For a single CNOT C with control c and target t, define
$`r_u=i^{x_c\oplus u_c}g_{C,u}`$ and the complete stage

```math
R_C=D_CS_C,\qquad
D_C|x,s,Y\rangle
 =i^{\sum_u s_uY_u(x_c\oplus u_c)}|x,s,Y\rangle,
\qquad R_C=\prod_u r_u^{s_uY_u}.                       \qquad(51)
```

The control bit is unchanged by every g_{C,u}, so the factors commute
as required. All x,s,Y inputs remain arbitrary. A private controlled-S
phase has an exact 12-T, 12-Clifford word: the four-T, four-CNOT
second difference returns its arbitrary dirty bit d and has exponent

```math
d-(d\oplus a)+(d\oplus a\oplus b)-(d\oplus b)
=-2ab+4dab\pmod8.
```

Its actual inverse followed by the private CCZ(d,a,b) gives
$`i^{ab}`$, returning both private bits.

For $`i^{abc}`$, use chronologically CS(d,c)-dagger,
Toffoli(a,b;d), CS(d,c), Toffoli(a,b;d), then C3Z(a,b,c,d).
The first four factors give $`i^{abc}(-1)^{abcd}`$; the last
removes that sign. The retained private C3X word, without its outer
target Hadamards, pays C3Z. Thus this controlled-controlled-S has
72 T gates, 80 Cliffords, at most 72 T layers and three returned
arbitrary private bits. All physical T targets are private.
For u_c=1, use CS(s_u,Y_u) times the actual inverse of this cubic
quarter phase; for u_c=0 use the cubic phase itself. Parallel private
row words and (44) give common forward/actual-reverse bounds

```math
T_{R_C}\le T_C+78m,\quad G_{R_C}\le G_C+86m,\quad
d_{R_C}\le124,\quad b_{R_C}\le\max(b_C,3m).             \qquad(52)
```

There are zero clean inputs and zero error, with all additional helpers
returned on the complete space; the m selectors remain principal inputs.
At one active row, $`(r_u+r_u^\dagger)/2=P_{x_c=u_c}`$.
For opposite CNOT directions at p=2, the two common +1 eigenspaces
therefore intersect only in the marked address. Their symmetrized
average is $`(P_{x_0=u_0}+P_{x_1=u_1})/2`$, with gap one half
below its marked eigenvalue one. This removes the preceding common
fixed-vector obstruction in that fixture. It does not make the old
commutator a source, or prove a uniform full-bank spectral gap.
The [flip-difference source below](#quarter-phase-flip-difference-source)
pays one nonunitary accepted polynomial of these stages. Its full lookup
Gram still fails. The exact tests check the private phases and the
p=2 projector identity; no complete source is supplied by this component.

<a id="quarter-phase-flip-difference-source"></a>

### A one-clean flip-difference source and its complete accepted Gram

For p at least two, the paid single-CNOT R_C gives an actual nonunitary
accepted block without a controlled-R_C implementation. Let F flip all
m selectors. With one initialized signal a, emit chronologically

```math
H_a,\quad \prod_u\mathrm{CX}(a,s_u),\quad R_C,\quad
\prod_u\mathrm{CX}(a,s_u),\quad Z_a,\quad H_a.
                                                               \qquad(53)
```

The Z is the literal two-S word. The complete native unitary U has
$`T_U=T_{R_C}`$, $`d_U\le124`$, $`G_U\le G_{R_C}+2m+4`$,
the same private pool, one clean signal and zero implementation error.
Counting the selectors, its dirty width outside the address and row bank
is $`b_U=m+b_{R_C}`$.
Its actual reverse has the same bounds. The selectors are arbitrary
principal work; only the stage's additional helpers return exactly.
For J inserting a=0, direct propagation of both signal branches gives

```math
K_-=J^\dagger UJ=(R_C-FR_CF)/2,\qquad
K_-^\dagger K_-=\tfrac12\bigl(I-\mathrm{Re}(R_C^\dagger FR_CF)\bigr).
                                                               \qquad(54)
```

Every address, row bank, selector, private-helper and reference input
is included. Nontrivial accepted fibers are allowed by the
[paid Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout);
U itself need not return its signal or selectors. That recovery would
give complete return only if the full target Gram were supplied.

For one active row u, (54) reduces to the literal identities

```math
K_-=P_{x_c\ne u_c}\,Z_{s_u}(I-iX_t)/2,\qquad
K_-^\dagger K_-=\tfrac12P_{x_c\ne u_c}\otimes I_s.       \qquad(55)
```

Thus the source has correct half normalization and an arbitrary-selector
unitary accepted fiber on the mismatch halfspace. It rejects the marked
address, but for p at least two it also misses the unmarked addresses
sharing that control bit. This is a complete partial-projector identity,
not a source for the whole lookup sector.

The general-bank action identifies the next failure precisely. At each
fixed control bit define the ordinary integer counts
$`N=\sum_{u:u_c\ne x_c}Y_u`$ and
$`N_s=\sum_{u:u_c\ne x_c}s_uY_u`$. All matching rows are identity
and each selected mismatching row contributes iX_t. Therefore

```math
R_C=(iX_t)^{N_s},\qquad R_C^\dagger FR_CF=(iX_t)^{N-2N_s},
\qquad K_-^\dagger K_-=
\begin{cases}
I/2,&N\text{ odd},\\
\bigl[1-(-1)^{N/2-N_s}\bigr]I/2,&N\text{ even}.
\end{cases}                                                   \qquad(56)
```

The even case depends on the arbitrary selector parity and takes values
zero or one, not one half. For p at least two choose a marked row x
and two further active rows with opposite control bit, and set s=0.
Then N=2 and N_s=0, so the actual Gram is one while the desired
$`(I+\Omega)/4`$ is zero. The full Gram distance is exactly one,
since both are positive contractions. Even an empty bank gives actual
zero versus desired one half. For the canonical paired $`\Pi_{\rm eq}/2`$,
put Y in the specified physical C rows and conjugate by the n bitwise
D-to-C CNOTs. Both orientations add 2n Cliffords, no T-depth, and
2n-m additional arbitrary carrier inputs. Total dirty width beyond the
specified m outputs is then $`2n+b_{R_C}`$. The gap persists; no
Vandermonde basis change or initialized carrier is needed.

Pause this single-stage difference. The [coupled-source audit](COUPLED_QUARTER_PHASE_SOURCE.md)
now pays two-stage and four-path acceptance, but exposes coherent Gram
defects. The [matching-family audit](QUARTER_PHASE_MATCHING_LIMITS.md)
proves a one-row spectral escape. The remaining question is complete
signal-interleaved acceptance that removes coherent interference and
retains arbitrary-bank normalization. The six-wire, 20-T, 28-Clifford native fixture fixes
one row active and checks every address/selector/helper/signal column
and its actual reverse; it is not an arbitrary-bank native emitter.
The separate algebra fixture checks all p=2 row-bank/selector fibers
of (54)--(56). These finite checks support
the analytic all-size action and sharp failure; they do not supply the
missing source or external mathematical review.
