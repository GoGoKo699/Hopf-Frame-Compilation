# Complete-input limits of two approximate equality filters

[Approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Retained complete indicator](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Native bilinear phase word](PARALLEL_BILINEAR_WORD_QUERY.md#1-a-phase-derivative-on-an-arbitrary-dirty-wire)

Constant-in-address source accuracy is enough for the conditional frame
transfer, but pointwise equality accuracy does not provide that source.
This chapter checks two concrete candidates. An additive real polynomial
phase filter has a sharp complete-operator error floor of square root
two. A native coherent hash filter approximates one controlled phase on
initialized work, but has error two on arbitrary dirty seeds and on a
complete indicator bank even when the seed is initialized. Neither
statement constrains general nonlinear access to the output bank or
arbitrary Clifford/Hadamard interference. The frontier is unchanged.

## 1. The literal additive phase-filter contract

Let $`p\ge1`$, $`m=2^p`$, let x be a p-bit address and let Y be an
arbitrary m-bit bank. For a real polynomial P consider the formal unitary

```math
U_P|x,Y\rangle
 =\exp\!\left(i\pi\sum_{u\in\mathbb F_2^p}
                      Y_uP(|x\oplus u|)\right)|x,Y\rangle.
```

The desired phase is $`D|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle`$.
Hadamards on all Y wires before and after these operators turn D into
the complete dirty indicator and preserve their operator distance.
The metric is literal: no Y-dependent phase, discarded work or
sectorwise phase identification is allowed. The polynomial coefficients
are arbitrary real numbers; no native realization of their rotations is
assumed. Thus failure of this formal word precedes its synthesis cost.

**Sharp filter bound.** If $`\deg P\lt p`$, then

```math
\|U_P-D\|_{\mathrm{op}}\ge\sqrt2.                         \qquad(1)
```

This concerns the displayed additive dependence on Y, not all circuits
with low-degree intermediate expressions.

### Proof, including phase wraparound

Suppose the norm were at most $`\tau\lt\sqrt2`$, and put
$`\beta=2\arcsin(\tau/2)/\pi\lt1/2`$. Fix x equal to zero and
write $`r_u=P(|u|)-\mathbf1[u=0]`$. For every subset S of labels,

```math
\sum_{u\in S}r_u\pmod2\ \in[-\beta,\beta].             \qquad(2)
```

Singleton subsets supply representatives $`e_u=r_u-2k_u`$ in that
interval, with integer k_u. All subset sums of e obey the same condition.
If the sum of positive e_u exceeded beta, the first partial sum crossing
beta would lie in $`(\beta,2\beta]`$. This interval lies strictly
between the allowed modulo-two arcs, since $`2\beta\lt2-\beta`$.
Consequently the positive total is at most beta. The negative total has
absolute value at most beta by the same argument. Therefore

```math
\sum_u|e_u|\le2\beta\lt1.                              \qquad(3)
```

The pth finite difference of a degree-less-than-p polynomial vanishes:

```math
\sum_u(-1)^{|u|}P(|u|)
 =\sum_{j=0}^p(-1)^j\binom pj P(j)=0.
```

Subtracting the delta and the even integers gives

```math
\sum_u(-1)^{|u|}e_u=-1-2\sum_u(-1)^{|u|}k_u.
```

This is an odd integer, whose absolute value is at least one, contradicting
(3). The argument permits arbitrary coefficient wraparound; it does not
assume that the original polynomial errors are already small real numbers.

### Sharpness despite pointwise error one over m

Prescribe the p+1 values

```math
P_*(j)=\mathbf1[j=0]-\frac{(-1)^j}{m},
\qquad 0\le j\le p.                                  \qquad(4)
```

Their pth finite difference is zero, so their interpolation polynomial
has degree at most p-1. Equivalently, (4) deletes just the full-parity
Fourier term from equality. Every coefficient error has magnitude 1/m.
There are m/2 positive and m/2 negative errors. All subset sums therefore
lie in $`[-1/2,1/2]`$, and choosing one parity class attains an endpoint.
It follows that $`\|U_{P_*}-D\|=\sqrt2`$. Pointwise error as small
as 1/m still does not imply useful complete-bank accuracy in this word.

## 2. A changed mechanism: coherent parity hashing

The next candidate uses only native hash arithmetic and a small exact
predicate, without a supplied fine-angle rotation. Its desired single
phase is

```math
M_p|x,z\rangle=(-1)^{z\mathbf1[x=0]}|x,z\rangle.
```

Choose k positive integer parity tests. A k-by-p binary matrix A is a
seed register, and t is a k-bit hash register. Initialize these pk+k
qubits to zero, put A in its uniform superposition with pk Hadamards,
compute $`t\leftarrow t\oplus Ax`$, apply the exact phase
$`(-1)^{z\mathbf1[t=0]}`$, actually uncompute the hash, and apply
the seed Hadamards again.

For nonzero x, exactly the fraction $`2^{-k}`$ of all matrices A
satisfy Ax equal to zero. For x equal to zero the desired phase is exact.
At nonzero x and z equal to one, precisely those colliding seed components
have the wrong sign. Orthogonality of the seed states gives the exact
complete initialized-isometry error

```math
\|VJ-J(M_p\otimes I)\|_{\mathrm{op}}=2^{1-k/2}.        \qquad(5)
```

J initializes A and t only; private dirty work is arbitrary. The norm
includes seed return error, all address superpositions, arbitrary z and
references. The hash register and private dirty work return exactly.
This is not a classical success probability asserted to be an operator
bound: (5) follows by taking the norm of the wrong-sign components in
each preserved x,z block.

### The hash arithmetic and actual inverse are paid

Use one private dirty bit $`d_{lj}`$ for each product $`A_{lj}x_j`$.
The R91 phase word has controls $`y=A_{lj},z=t_l`$ and toggles d
by $`s=x_j`$. Its exact phase is $`(-1)^{A_{lj}x_jt_l}`$,
with d returned. All pk instances have distinct T targets, while their
shared inputs occur only as Clifford controls. Hadamards on t give the
complete shear B computing Ax, even on arbitrary initial t:

```math
T_B\le8pk,\qquad G_B\le18pk+2k,\qquad D_T(B)\le8,
\qquad b_B=pk.                                           \qquad(6)
```

The Clifford upper bound uses literal parity compute/uncompute for each
term, rather than an optimized Gray walk. Let $`T_0,G_0,d_0,b_0`$
be common forward/actual-inverse bounds for the exact negative-control
k-bit zero-test phase. The retained dirty conjunction supplies polynomial
count and width and logarithmic depth in k, with no further clean work.
The entire single-phase candidate has common bounds

```math
\begin{aligned}
T&\le16pk+T_0,&G&\le38pk+4k+G_0,\\
D_T&\le16+d_0,&b_{\mathrm{dirty}}&\le\max(pk,b_0),\\
a_{\mathrm{clean}}&=pk+k.
\end{aligned}                                                   \qquad(7)
```

The seed Hadamards, both hash orientations and the zero test are included.
The private dirty pool is reused only after completed return. For fixed
desired error, k can be constant in (5), but the initialized allocation
in (7) is essential to that conclusion.

For an arbitrary dirty seed, fix any nonzero x and any colliding matrix
$`A_0`$, for example the zero matrix. Choose the initial seed state
$`H^{\otimes pk}|A_0\rangle`$, and take t zero and z one. The first
Hadamards produce that bad seed exactly. The candidate has phase minus
one while the target has phase plus one, and the seed returns to its
initial state. Thus its complete arbitrary-seed operator error is two.
More dirty storage does not initialize the seed of this particular word.

## 3. The complete bank fails even with the initialized seed

Reuse the hash to test every output label u. With t equal to Ax, the
formal completed phase for a fixed A is

```math
(-1)^{\sum_uY_u\mathbf1[Au=Ax]}.
                                                               \qquad(8)
```

This whole-bank candidate can also be charged natively. Its row controls
are the stable affine literals
$`1+t_l+\sum_j A_{lj}u_j`$. The retained
[read-only affine-control ladder](../../docs/DIRTY_SUM_COMPRESSION.md#read-only-affine-controls-with-private-phase-work)
implements each row's phase with O(k) T gates, private dirty work and
T-depth, and $`O(k(p+1))`$ Cliffords. Different rows use distinct
Y targets and helper pools; the shared A,t wires are only CNOT controls.
The constant affine terms are Clifford-paid. Consequently all rows may
run in parallel. Including hash computation,
its actual inverse, seed Hadamards and indicator endpoints gives

```math
\begin{aligned}
T&=O(mk+pk),&G&=O(mk(p+1)+pk),\\
D_T&=O(k+1),&b_{\mathrm{dirty}}&=O(mk+pk),\\
a_{\mathrm{clean}}&=pk+k.
\end{aligned}                                                   \qquad(9)
```

For fixed k these counts have the desired linear T/dirty and quadratic
Clifford shape, and constant native depth, if the clean work were granted.
The operator is nevertheless wrong by the maximum possible amount.

If $`k\lt p`$, every fiber $`x+\ker A`$ has even size, including
rank-deficient seeds. On the allowed phase-basis input with every Y_u
equal to one, (8) is plus one for every seed A. The desired phase is
minus one. This proves error exactly two, even for an initialized
coherent seed and even if all native components are implemented exactly.
The unwanted bucket toggles are not a small error in the complete
indicator contract. Hadamard endpoints preserve this distance.

## 4. The precise limit of nonlinear summary decoding

A nonlinear decoder of the same parity summaries does not automatically
repair the bucket construction. Suppose its only access to the
computational Y values is through

```math
\ell_f(Y)=\bigoplus_u f(u)Y_u,
```

where, for each fixed address and seed sector, f is a Boolean polynomial
over $`\mathbb F_2`$ of degree less than p. Every such f has even
truth-table weight: each monomial omits an index and therefore sums to
zero over the cube. Hence
$`\ell_f(Y\oplus\mathbf1)=\ell_f(Y)`$. Hash-bucket indicators
with fewer than p parity tests are examples.

Any reversible, nonlinear or coherent processing of these summaries
has the same work action for Y and its complement, whereas the required
sign $`(-1)^{Y_x}`$ changes. For any common resulting unit vector v
and ideal returned vector a, the two squared distances satisfy

```math
\|v-a\|^2+\|v+a\|^2=4.
```

At least one is therefore at least square root two. This argument also
permits a coherent seed: complementation fixes the summary in every
seed sector. It assumes Y is computationally preserved during the
summary program, apart from the final common Hadamard endpoints.
Additional direct Y access, interleaved changes of its basis, or features
outside this specified family are not covered. A raw Y correction can
break the displayed complement symmetry; its complete selection and
return costs would then need a new proof. This is an interface diagnosis,
not a lower bound on unrestricted nonlinear-Y readout.

## 5. What a successful replacement must supply

The initialized single-predicate bound (5) does not supply a zero-clean
complete indicator. Nor does independently replacing dirty tree nodes
by bounded-error predicates give a maximum-over-branches estimate:
unknown dirty internal flags do not make their errors mutually exclusive.
Any such proposal still owes its full-word error and return proof.

A source that actually implements the complete indicator at requested
fixed error, with $`O_\tau(m)`$ T and dirty width,
$`O_\tau(m^2)`$ Cliffords and $`o(\log(p+2))`$ T-depth,
would enter the paid approximate R91/R98 transfer. A whole-word
nonlinear or interference mechanism could meet that contract without
satisfying either filter premise above. Neither candidate here provides
that missing native operation.

The [focused fixtures](../../tests/test_approximate_indicator_filter_limits.py)
check bounded exact polynomial and native hash identities. The
all-size norm bounds and resource accounting are analytic; finite
checks do not construct a faster complete source.
