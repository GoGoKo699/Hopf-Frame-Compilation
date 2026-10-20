# Spectral clean-source promotion with paid paving frames

[Complete clean-source contract](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Parallel promotion and its scalar filter](PARALLEL_CLEAN_PROMOTION.md) ·
[Native sharing audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md)

A capped selective echo and a finite paving construction pay the global
frame gap needed by the retained spectral filter. The frames are general
Clifford+T words, not necessarily Clifford frames. At fixed accuracy this
gives a zero-clean promotion with native depth
$`O_\tau(d_V+k)`$ and explicitly bounded count and dirty width. Exact
preservation of the source's clean sector permits constant source
accuracy, independent of k. This enlarges the sufficient-source window
to $`k=o(\log p)`$ under the resource assumptions in Section 7.
An actual qualifying shallow clean source remains absent.

Every norm below is a complete operator or initialized-isometry norm,
including arbitrary source work, borrowed registers and references.
Every inverse is the actual reversed emitted word. No classical
precomputation time bound is claimed for selecting the paving; a finite
exact terminating selection algorithm is given.

## 1. Source contract and synchronized native stages

Let V be the complete source query of
[clean promotion Section 1](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#1-the-source-is-a-complete-query-not-a-prepared-state).
It has k clean bits A, q arbitrary old output bits T, unchanged address
x, and arbitrary private work. Write J for insertion of A=0 and P=JJ†.
The exact-sector approximate contract is

```math
 [V,P]=0,\qquad \|VJ-JQ_f\|\le\eta.                  \qquad(1)
```

Arbitrary mixing within the nonzero A sector is allowed. The alternative
contract assumes only initialized-isometry error epsilon, with possible
leakage. Common forward/reverse resource bounds are T_V,G_V,d_V,b_V,
excluding x,T,A. The case k=0 is already a zero-clean source and needs
no promotion; assume k≥1 below.

Require V's actual completed Clifford boundaries to fix x and its
physical T targets to be private. Set $`M=h=2^k`$. Copy j has private
registers A_j,T_j, work, weight t_j, toggle e_j and quartic-phase work.
Only x and the specified q-bit output Z are shared. Put

```math
\begin{aligned}
 F(x,Z)&=f(x)\mathbin{\cdot}Z,\\
 q_j&=[A_j=0^k],\qquad K(A)=\sum_jq_j,\\
 \chi(A)&=[1\le K(A)\le2],\\
 g_\chi(A,t)&=\chi(A)\bigoplus_jq_jt_j.
\end{aligned}                                                       \qquad(2)
```

The source stage V_all is the parallel product of the V_j. During that
stage each A_j is private and may change. The selective stage defined
next reads every A_j and therefore runs only after the complete source
stage has ended. Independent local echoes are not interleaved across
these changing shared key registers.

## 2. A complete capped selective echo

For an arbitrary bit u and an additional stable control a, define

```math
\begin{aligned}
 K_{\rm all}&=\prod_{j,\ell}
  \mathrm C^3\mathrm Z(e_j,t_j,T_{j,\ell},Z_\ell),\\
 C_{\rm all}:\quad e_j&\longmapsto
  e_j\oplus a u\chi(A)q_j,\\
 \mathcal E(C,D)&=D^\dagger C^\dagger DC.
\end{aligned}                                                       \qquad(3)
```

Emit the synchronized whole word

```math
 \mathcal E(C_{\rm all},
       V_{\rm all}^\dagger K_{\rm all}V_{\rm all})
 \mathcal E(C_{\rm all},K_{\rm all}).                  \qquad(4)
```

For an exact source this equals $`(-1)^{auF g_\chi}`$ on the entire
space. To check both cancellation and accuracy, exact sector
preservation implies that V_all commutes with every q_j and with chi.
The e-toggle derivative of K_all is

```math
 B=\prod_{j:au\chi q_j=1}
       (-1)^{t_j T_j\cdot Z}.
```

Thus (4) is algebraically V_all† B V_all B. On chi=0, or au=0, it is
exactly identity. In a sector with chi=1 only one or two B_j factors
are active. Every inactive, potentially arbitrary off-code V_j cancels.
Each active clean-sector conjugation has error at most 2eta, giving

```math
 e_g\le
 \begin{cases}
 4\eta,&\text{contract (1)},\\
 8h\epsilon,&\text{initialized error with leakage}.
 \end{cases}                                         \qquad(5)
```

The first bound is a maximum over orthogonal joint q-sectors, not a
sum over all h copies. The second follows by replacing each V_j with
the retained exact block extension and telescoping the four literal
source appearances. Neither estimate drops arbitrary source-work or
reference errors. The physical receipt retains all four source stages;
the algebraic cancellation is used only to prove the operator contract.

Here is an explicit vector-query payment for C_all. Use the
[retained exact dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md), with
the [read-only native schedule](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences).
For an r-bit address its common forward/reverse T, G and private-width
bounds are $`2^{O(r)}\mathrm{poly}(r)`$, and its depth is
O(log(r+2)). Let $`m=hk+2`$
count all key bits and the stable controls a,u. A complete m-address
dirty indicator I_m writes equality toggles into an arbitrary private
$`2^m`$-bit bank Y. For each output e_j, XOR its specified truth-table
mask of Y into e_j. Apply I_m† and repeat every mask. All Y and indicator
work return, and every e_j has exactly the toggle in (3). If b_I excludes
the indicator's principal output bank, safe bounds are

```math
\begin{aligned}
 T_D&\le2T_I(m),\\
 G_D&\le2G_I(m)+2h2^m,\\
 d_D&\le2d_I(m),\\
 b_D&\le2^m+b_I(m).
\end{aligned}                                                       \qquad(6)
```

All masks are Clifford CNOT words. The indicator's only shared-key
interactions are read-only Clifford controls; all T targets are private.
The whole selective stage therefore has native depth O(log(m+2)),
with its complete exponential bank charged in (6).

Use the retained literal quartic-phase receipt from
[parallel promotion Section 2](PARALLEL_CLEAN_PROMOTION.md#2-a-paid-selective-phase-with-an-arbitrary-weight-bit):
one term costs 32 T gates, 44 Cliffords, depth 32 and two dirty bits.
Equation (4), including actual adjoints, has

```math
\begin{aligned}
 T_g&\le4hT_V+4T_D+128hq,\\
 G_g&\le4hG_V+4G_D+176hq,\\
 d_g&\le4d_V+4d_D+128,\\
 w_g&\le h(k+3q+2+b_V)+b_D.
\end{aligned}                                                       \qquad(7)
```

w_g excludes x,Z,u and the two later clean signals a,b. All work in
(7) is arbitrary. The word is controlled-a G_{uF}, without a controlled
source V. Off au=0 it is literally identity, including for approximate V.

## 3. A finite paving gives a whole-space gap

The key-and-weight register R has $`n=h(k+1)`$ bits. Let P_0 project
onto the computational labels with $`g_\chi=0`$. Under uniform basis
counting, K is binomial with h=M and success probability 1/M. For M≥2,

```math
\begin{aligned}
 \Pr(K=1)&=(1-1/M)^{M-1}>e^{-1},\\
 \Pr(K=2)&=\tfrac12\Pr(K=1),\\
 \delta:=2^{-n}\mathrm{rank}(P_0)
  &=1-\tfrac12\Pr(1\le K\le2)<\tfrac34.
\end{aligned}                                                       \qquad(8)
```

The last equality holds because, in every active key sector, the nonzero
parity in the arbitrary weight bits is balanced. Put W=H on every R bit
and Q=WP_0W. Every diagonal entry of Q is delta.

[Marcus, Spielman and Srivastava, *Interlacing Families II* (2015),
Corollary 1.5](https://arxiv.org/pdf/1306.3969) partitions an isotropic
family of vectors of squared norm at most delta into r parts whose
frame operators have norm at most
$`(1/\sqrt r+\sqrt\delta)^2`$. Apply it to the vectors Q|v> in the
range of Q, with r=64. It gives a coloring c(v) with six bits such that
each coordinate-color projection Pi_j satisfies

```math
\begin{aligned}
 \|\Pi_jQ\Pi_j\|&\le\alpha,\\
 \alpha&=(1/8+\sqrt\delta)^2
   <49/64+\sqrt3/8<63/64.
\end{aligned}                                                       \qquad(9)
```

The frame theorem first bounds Q Pi_j Q on the range of Q; the two
positive products have the same nonzero eigenvalues. No normalized
trace or typical-state estimate replaces the operator norm in (9).

This is a terminating exact classical construction. Enumerate the finite
six-bit colorings and retain the first for which every rational matrix
$`(63/64)I-\Pi_jQ\Pi_j`$ is positive semidefinite. Q has dyadic rational
entries, so exact rational matrix tests decide each candidate. The strict
slack in (9) guarantees termination. The emitted word's resource ledger
does not claim that this exhaustive classical search takes polynomial
time.

For b in the six-dimensional binary cube define

```math
 D_b|v\rangle=(-1)^{b\cdot c(v)}|v\rangle,
 \qquad C_b=WD_bW.
```

Character orthogonality removes cross-color entries, so

```math
\begin{aligned}
 \frac1{64}\sum_bD_bQD_b
   &=\sum_j\Pi_jQ\Pi_j\le\frac{63}{64}I,\\
 H_1:=\sum_b\frac{I-C_b(-1)^{g_\chi}C_b^\dagger}{2}
   &=64I-\sum_bC_bP_0C_b^\dagger,\\
 I&\le H_1\le64I.
\end{aligned}                                                       \qquad(10)
```

For uF=0 the corresponding H is zero. This proves the complete spectral
gap, including every arbitrary R input. It supplies general native
frames; it does not prove that Clifford-only frames suffice.

## 4. Literal native frames and their entire borrowed bank

For one coloring character f_b(v)=b·c(v), use a complete n-address dirty
indicator I_n with its own arbitrary $`2^n`$-bit output bank Y. Put
$`Z_{f_b}=\prod_{v:f_b(v)=1}Z_{Y_v}`$. Then

```math
 I_n^\dagger Z_{f_b}I_nZ_{f_b}
   =D_b\otimes I_Y\otimes I_{\rm work}.               \qquad(11)
```

Every old Y phase cancels. With Z implemented as SS, common forward and
actual-inverse frame bounds are

```math
\begin{aligned}
 T_C&\le2T_I(n),\\
 G_C&\le2G_I(n)+4\,2^n+2n,\\
 d_C&\le2d_I(n),\\
 b_C&\le2^n+b_I(n).
\end{aligned}                                                       \qquad(12)
```

The last 2n Cliffords are W and W†. The same returned pool is reused
across all 64 frames. Thus T_C,G_C,b_C are $`2^{O(n)}\mathrm{poly}(n)`$
and d_C=O(log(n+2)). The truth table, indicator output bank and native
word are all paid. No QRAM or prepared sign state is supplied.

The controlled framed operation is C_b†, controlled G_{uF}, C_b in
chronological order. The frame is not controlled. On control zero its
two actual inverse occurrences cancel. Its error remains e_g.

## 5. The two-clean scalar filter and complete tolerance choice

Use the literal construction of
[parallel promotion Section 6](PARALLEL_CLEAN_PROMOTION.md#6-alternative-a-constant-spectral-gap-permits-a-paid-scalar-filter)
with s=64 and gamma=1 from (10). Only its two signals a,b are clean.
For an ideal framed involution G_i, H_a ctrl_a(G_i) H_a computes its
eigenvalue bit. A literal two-bit controlled phase, followed by the
actual reverse, simulates a controlled exponential of
$`P_i=(I-G_i)/2`$. Its three determinant-one R_z rotations retain the
common scalar phase; the scalar cancels in the final conjugation.

Put $`t=\pi/(2s)`$, choose L first-order product-formula steps and let
B=sL. On a=0 the ideal extension simulates controlled exp(−itH).
On a=1 it instead simulates controlled exp(−it(sI−H)), with the same
retained scalar. It is not an identity extension on arbitrary a. The
same Trotter bound holds on both blocks. If each one-qubit exponential
has error at most mu_z, the Hadamard-test word A obeys

```math
 \delta_A\le\frac{\pi^2}{8L}+B(2e_g+3\mu_z).          \qquad(13)
```

On the ideal a=0 sector, A is scalar identity for uF=0. For uF=1 its
b=0 failure amplitude is at most $`c=\cos(\pi/256)\lt1`$.

The retained determinant-one fixed-point recurrence is

```math
\begin{aligned}
 R&=\mathrm{diag}(e^{i\pi/6},e^{-i\pi/6})\text{ on }b,\\
 S_0&=A,\qquad S_{j+1}=S_jRS_j^\dagger R^\dagger S_j.
\end{aligned}                                                       \qquad(14)
```

Its failure amplitude cubes. This is the explicitly phased adaptation
of [Grover's pi/3 step](https://arxiv.org/pdf/quant-ph/0503205) proved in
the retained filter; it uses the same two clean bits at every level.
For N=3^j, the complete initialized query
$`V_F=H_uS_j^\dagger Z_bS_jH_u`$ has error

```math
 \epsilon_F\le2c^N+2N\delta_A+2(N-1)\mu_R,             \qquad(15)
```

where mu_R is the full-norm synthesis error of R. Both clean returns and
every arbitrary borrowed-input error are included. All approximants use
their actual adjoints. No approximate G_i is assumed Hermitian: (13)–(15)
compare the actual word with the exact Hermitian ideal by telescoping.

For a desired post-fusion tolerance tau≤1, set sigma=tau/8. Choose

```math
\begin{aligned}
 2c^N&\le\sigma/128,\qquad N=3^j,\\
 L&\ge192\pi^2N/\sigma,\qquad B=64L,\\
 e_g&\le\sigma/(3072NB),\\
 \mu_z&\le\sigma/(4608NB),\\
 \mu_R&\le\sigma/(512N).
\end{aligned}                                                       \qquad(16)
```

Then delta_A≤sigma/(512N), epsilon_F≤sigma/64, the zero-clean phase
below has error at most sigma, and root fusion multiplies that error by
at most eight. All constants N,L and rotation accuracies depend only on
tau. Sufficient source tolerances are

```math
 \eta\le\frac{\sigma}{12288NB}
 \quad\text{in (1)},\qquad
 \epsilon\le\frac{\sigma}{24576hNB}
 \quad\text{with leakage}.                           \qquad(17)
```

Thus exact clean-sector preservation retains constant source accuracy.
The general leakage contract still needs a uniform shrinking-accuracy
source ledger.

## 6. All resources refer to the same emitted word

Let T_z,G_z,d_z and T_R,G_R,d_R bound the finite native single-qubit
rotation approximants in (13)–(16), also for their actual inverses.
Charging the non-Clifford frames in every simulation term gives

```math
\begin{aligned}
 T_A&\le B(2T_g+4T_C+3T_z),\\
 G_A&\le B(2G_g+4G_C+6+3G_z)+2,\\
 d_A&\le B(2d_g+4d_C+3d_z),\\
 T_F&\le2NT_A+2(N-1)T_R,\\
 G_F&\le2NG_A+2(N-1)G_R+4,\\
 d_F&\le2Nd_A+2(N-1)d_R.
\end{aligned}                                                       \qquad(18)
```

The four final Cliffords pay Z_b=SS and the two H_u gates. Apply the
retained complete promotion at two clean bits and one output, then the
two-query Z_u echo and final H_Z envelope. With C_2 the separate returned
two-bit zero-test, the complete pre-fusion query has zero clean work and

```math
\begin{aligned}
 T_{\rm out}&\le32T_F+32T_{C_2}+224,\\
 G_{\rm out}&\le32G_F+32G_{C_2}+456+2q,\\
 d_{\rm out}&\le32d_F+32d_{C_2}+128,\\
 w_{\rm out}&\le w_g+b_C+b_{C_2}+5,\\
 \delta_{\rm out}&\le64\epsilon_F.
\end{aligned}                                                       \qquad(19)
```

Width excludes only the specified address x and output Z. It includes
all h source-clean banks now treated as dirty, both exponential
indicator banks in (6) and (12), source and phase work, the two promoted
signals, and the arbitrary final echo bit. The full-norm estimate
includes return of every borrowed register, even under general source
leakage. The actual inverse has the same error and resource bounds.

Equations (6), (7), (12), (16), (18) and (19) are simultaneous count,
depth, width and precision recurrences; no cost is attached to a
different or idealized implementation. In particular, the 4T_C and
4d_C terms cannot be omitted by treating the paving frames as Clifford.

## 7. Compression sources and the larger sufficient window

An input source U with k_0 clean bits and compression amplitude
$`1/\sqrt2`$ can first use the retained
[complete compression recovery](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md).
Its recovered V has k=k_0+1 clean bits and

```math
\begin{aligned}
 T_V&\le3T_U+4T_{C_{k_0}},\\
 G_V&\le3G_U+4G_{C_{k_0}}+42,\\
 d_V&\le3d_U+4d_{C_{k_0}},\\
 b_V&\le b_U+b_{C_{k_0}}.
\end{aligned}                                                       \qquad(20)
```

The recovery flag and its zero-test work are private. The completed
address-fixing boundary and private-T-target requirements therefore
persist. If its compression K satisfies the actual flat identity
$`K^\dagger K=I/2`$, V preserves its clean sector exactly and
$`\eta\le\sqrt2\epsilon_U`$. Without flatness the initialized V error
is at most $`3\sqrt2\epsilon_U+2\epsilon_U^2`$. Substitute these bounds
into (17), and (20) into the same resource recurrences above. No growing
clean allocation is omitted: the actual recovered width is k_0+1 before
the h parallel copies are made arbitrary.

For fixed tau, the resulting source promotion has

```math
\begin{aligned}
 d_{\rm out}&=O_\tau(d_V+k),\\
 T_{\rm out}&=O_\tau(h(T_V+q)+\Lambda),\\
 G_{\rm out}&=O_\tau(h(G_V+q)+\Lambda),\\
 w_{\rm out}&=O_\tau(h(b_V+q+k)+\Lambda),\\
 \Lambda&=2^{O(h(k+1))}\mathrm{poly}(h,k).
\end{aligned}                                                       \qquad(21)
```

The separate exact ledgers above govern each resource. If $`k=o(\log p)`$, then
$`h(k+1)=p^{o(1)}`$. The extra exponential banks are consequently
$`2^{p^{o(1)}}`$, and the multiplicative h factor is polynomially
absorbable in the retained $`2^p\mathrm{poly}(p)`$ T/width and
$`4^p\mathrm{poly}(p)`$ Clifford allowances. A source with
$`d_V=o(\log p)`$ yields that same native-depth target. Under the flat
contract, its budgets are needed at one fixed accuracy depending only
on the final tolerance. The leakage version instead requires the
explicit accuracy in (17).

Finally apply the [paid root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
at source address size $`b=p-O_\tau(\log(p+2))`$. Its eight actual
source/adjoint appearances and charged exact tail give

```math
\begin{aligned}
 T_{\rm root},w_{\rm root}&=O_\tau(2^p),\\
 G_{\rm root}&=O_\tau(4^p),\\
 d_{\rm root}&\le8d_{\rm out}(b)
   +O_\tau(\log\log(p+4)),\\
 \delta_{\rm root}&\le8\delta_{\rm out}\le\tau.
\end{aligned}                                                       \qquad(22)
```

The choice sigma=tau/8 in (16) already pays this last error factor.
Equation (22) is conditional on the stated complete source budgets,
not only on its initialized-state success probability.

This pays the former global-frame hypothesis in the general native
class and bypasses the radix-packing obligation. It does not implement
the monomial packing W, establish a Clifford-only frame theorem, or
construct a qualifying shallow source. The complete unrestricted
frontier therefore remains open.
