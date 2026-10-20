# Parallel clean-source promotion: packing and spectral alternatives

[Complete clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Native sharing audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

This is an exact conditional reduction. Independent dirty source banks
permit parallel source use even when the source genuinely uses its clean
bank between native boundaries. Section 5.1 pays the fixed packing by
exact permutation synthesis when its register is sufficiently small.
An actual qualifying clean source remains absent.

The source and error conventions are those of
[clean promotion Sections 1–4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#1-the-source-is-a-complete-query-not-a-prepared-state).
The native-sharing distinction is the one in the
[clean-sector parallelism audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md).
Section 6 retains a conditional extraction through a constant
whole-space gap for Clifford frames. The later
[MSS construction](MSS_SPECTRAL_PROMOTION.md) pays a different global
frame family with native non-Clifford Boolean phases, including their
T cost, and widens the sufficient window to clean count sublogarithmic
in p. Neither reduction constructs the missing shallow source itself.

## 1. Source contract and native sharing

Let V act on an address x, a q-bit old output T, k source-clean bits A,
and arbitrary private work. Let J insert A=0 and put P=JJ†. The exact
source contract is VJ=JQ_f, where Q_f adds f(x) to T and returns every
arbitrary private input. A fixed scalar multiplying the clean-code
action is harmless here because every V occurs inside a conjugation.

The first approximate variant assumes exact invariance of P and

```math
\|VJ-JQ_f\|\le\eta.
```

The second assumes only the initialized-isometry error epsilon, allowing
leakage from P. All norms are full operator norms, stable under tensoring
an arbitrary reference. Common forward/actual-inverse resource bounds
are T_V, G_V, d_V, b_V; b_V excludes x,T,A.

For the parallel native schedule, require that every actual completed
Clifford boundary of V fix the computational address x, and that its
physical T/T-dagger targets be private. No boundary need fix A: each
copy receives its own arbitrary A bank. Source words with genuine
address motion are outside this schedule. The final phase-output Z is
also read-only during the parallel batch.

Put M=2^k and take h=M source copies. Copy j has private registers
A_j,T_j, source work, a zero-test toggle e_j, an additional arbitrary
weight bit t_j, and the helper pools below. Only x and the specified
q-bit output Z are shared. None of these private registers is initialized.

## 2. A paid selective phase with an arbitrary weight bit

Write F(x,Z)=f(x)·Z over F_2. For one copy define

```math
\begin{aligned}
B_j&=\prod_{\ell=1}^q\mathrm{CCZ}(t_j,T_{j,\ell},Z_\ell),\\
K_j&=\prod_{\ell=1}^q\mathrm{C}^3\mathrm{Z}
 (e_j,t_j,T_{j,\ell},Z_\ell),\\
D_j&=V_j^\dagger K_jV_j,\\
C_j&:\ e_j\longmapsto e_j+[A_j=0^k].
\end{aligned}
```

C_j is the exact
[returned dirty zero-test](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences),
with common bounds
T_C,G_C,b_C≤c(k+3)^3 and d_C≤c log_2(k+3). Its private pool is disjoint
from the source pool. For an emitted unitary R put

```math
\mathcal E(C,R)=R^\dagger C^\dagger RC.
```

The [selective-echo proof](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#2-a-selective-native-echo-removes-the-unknown-off-code-action)
applies with t_j incorporated into
the phase. With actual adjoints in every occurrence,

```math
S_j=\mathcal E(C_j,D_j)\mathcal E(C_j,K_j)
   =(-1)^{F(x,Z)t_j[A_j=0^k]}.                    \qquad(1)
```

The equality is on the entire physical space. Arbitrary off-code source
behavior cancels; no controlled source is used. One S_j has four source,
four zero-test and four K_j occurrences, including their actual inverses.
For an exactly sector-preserving approximate source its error from (1)
is at most $`2\eta`$. Under only initialized-isometry error epsilon,
the [exact block-extension comparison](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#4-approximate-clean-return-still-gives-a-full-operator-query)
gives at most $`8\epsilon`$.

Here is a literal all-private-T-target implementation of K_j. For a
term abcd borrow arbitrary bits rho and w. The standard returned cubic
phase uses only w as a T target and implements CCZ(u,v,s) with
T=8, G=10, d≤8 and one returned dirty helper. Explicitly, the four-T
Gray word on w has exponent

```math
p(w;u,v)=w-(w\oplus u)-(w\oplus v)+(w\oplus u\oplus v)
        =(4w-2)uv.
```

The [literal private-phase receipt](BALANCED_ALGEBRA_PERMUTATION.md#4-literal-native-bilinear-word-and-exact-counts)
gives an echo under $`w\mapsto w+s`$ with phase $`(-1)^{uvs}`$,
eight T/T-dagger gates
and ten CNOTs. Conjugating its target by H gives the returned Toffoli
Q:rho→rho+ab, with T=8,G=12,d≤8. Let P_3=CCZ(rho,c,d), using the
same private w. The chronological sequence Q,P_3,Q,P_3 has phase
(-1)^{abcd} and returns both rho and w. Distinct terms use distinct
rho,w pairs. Their shared variables occur only as Clifford controls,
so each of the four stages runs across all q terms in parallel. Thus

```math
T_K\le32q,\qquad G_K\le44q,\qquad d_K\le32,
\qquad b_K=2q.                                    \qquad(2)
```

These bounds also hold for the literal reverse. In particular the
four-input controlled phase is not assumed to have a free ancilla-free
Clifford+T implementation.

Set

```math
\begin{aligned}
T_S&=4T_V+4T_C+128q,\\
G_S&=4G_V+4G_C+176q,\\
d_S&=4d_V+4d_C+128.
\end{aligned}                                      \qquad(3)
```

At a fixed shared x,Z sector, actual Clifford boundaries from different
copies act on disjoint private registers, times commuting scalar phases.
The T bands have disjoint physical targets. Combining corresponding
boundaries and T bands therefore emits all h selective words in depth
at most d_S, with T≤hT_S and G≤hG_S. This proves a native schedule;
commutativity of the completed operators alone would not prove it.

Let R consist of the h pairs (A_j,t_j), and define

```math
g(R)=\bigoplus_{j=1}^h t_j[A_j=0^k],\qquad
G_F=\prod_{j=1}^h S_j=(-1)^{F(x,Z)g(R)}.             \qquad(4)
```

For approximate sources, the corresponding emitted batch has error at
most 2hη or 8h epsilon. There is no direct-sum improvement: many private
A_j banks may simultaneously lie in their source-clean sectors.

## 3. Exact packing and the complete echo

The register R has d=h(k+1) bits and dimension L=2^d. The truth-table
weight r of g satisfies

```math
r=\frac{L-[2(M-1)]^h}{2}
 =\frac L2\left(1-(1-1/M)^h\right).                 \qquad(5)
```

For h=M, (1−1/M)^M<1/2, so L/4≤r≤L/2. Borrow an additional arbitrary
bit z and choose one existing bit b of R. On R,z define the subsets

```math
\begin{aligned}
E&=\{(R,z):g(R)=1\},\\
E'&=\{z=1,b=1\}
 \ \cup\ \{z=0,\ \mathrm{rank}(R)\lt2r-L/2\}.
\end{aligned}                                      \qquad(6)
```

Both have cardinality 2r. Here rank is the integer represented by a
fixed ordering of R. A concrete mathematical permutation W maps the
lexicographically ordered elements of E to those of E', and likewise
their complements. This specifies W without reference to f,x,Z.
That basis ordering alone is not a native implementation of W;
Section 5.1 pays its synthesis in the small-register window.

More generally any unitary W with WΠ_EW†=Π_E' suffices. An implementation
with arbitrary returned helpers must satisfy this identity after
tensoring both projectors by the identity on those helpers; a subspace
or initialized-helper promise does not suffice.

Put C_0=CNOT(z→b). On z=0 it fixes every basis vector. On z=1 it
complements membership in E'. Consequently

```math
E'\mathbin\triangle C_0(E')=\{z=1\}.
```

With P_F'=WG_FW†, form the complete words

```math
\begin{aligned}
E_F&=P_F'{}^\dagger C_0^\dagger P_F'C_0=(-1)^{Fz},\\
D_F&=E_F^\dagger X_zE_FX_z=(-1)^F I_{\rm private}.
\end{aligned}                                      \qquad(7)
```

Every equality holds for every R,z and every other arbitrary input.
The final H_Z envelope converts D_F to the complete XOR Q_f on the
specified output Z. Thus the identity neither creates clean copies nor
uses postselection, an initialized impulse, or a free uniform-state
reflection. For approximate emitted G_F or W, the expression uses their
actual adjoints, not an assumed involution of an approximate word.

## 4. Simultaneous conditional ledger and precision

Suppose the exact W above has a zero-clean implementation with common
forward/actual-inverse bounds T_W,G_W,d_W,b_W. Four G_F/actual-adjoint
occurrences and eight W/actual-adjoint occurrences appear in D_F. It
also has four C_0 gates and two X_z gates. In the positive-S native
alphabet X=HSSH, giving the common complete-query bounds

```math
\begin{aligned}
T_{\rm final}&\le4hT_S+8T_W,\\
G_{\rm final}&\le4hG_S+8G_W+12+2q,\\
d_{\rm final}&\le4d_S+8d_W,\\
b_{\rm final}&\le h(k+3q+2+b_V+b_C)+1+b_W.
\end{aligned}                                      \qquad(8)
```

The width excludes x and the specified Z output. In each copy it
includes A_j,T_j,e_j,t_j, source work, zero-test work and the 2q
private phase helpers. The extra one is z. Disjoint pools are charged;
no helper is assumed initialized. Exact sources and exact W return
every borrowed input literally. For approximate sources, all possible
address, output, helper and reference errors are included in

```math
\begin{aligned}
\|Q_{\rm final}-Q_f\otimes I_{\rm private}\|
 &\le8h\eta &&\text{(exact source-sector invariance)},\\
 &\le32h\epsilon &&\text{(general initialized-isometry promise)}.
\end{aligned}                                      \qquad(9)
```

If an emitted W approximates the stipulated exact W in full operator
norm by mu, add 8mu to either bound; its actual inverse has the same
error. This requires a literal full-space approximation, including
phase and arbitrary-helper return. The actual inverse of the final
query has the same complete error.

For exactly flat amplitude-$`1/\sqrt2`$ compression with error epsilon,
the [retained general-k recovery](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag)
gives an exactly invariant recovered
clean sector with η=sqrt(2)epsilon. Its additional recovery flag must
be included in the k used here. Equation (9) then becomes
8sqrt(2)h epsilon. No new recovery recurrence is assumed.

For a general unitary packing, target complete error tau therefore
requires η_b≤tau/(8·2^k), or epsilon_b≤tau/(32·2^k), with an additional
allowance if W is approximate. Fixed-accuracy polynomial source bounds
alone do not establish that varying-accuracy claim. The next subsection
removes this accuracy multiplier for exactly sector-preserving sources
by prescribing a stronger monomial packing.

The remaining packing obligation is a W with the required full-space
identity and simultaneous count, width and depth budgets. Section 5.1
pays it by generic exact synthesis for a restricted register-size window.
The constant-degree phase gates and source parallel schedule have already
been paid in (2)–(3).

### 4.1 A refined monomial packing retains constant source accuracy

Assume exact source-sector preservation, while allowing the clean-code
error η. Set a_j(R)=t_j[A_j=0] and let N(R) be their ordinary integer
sum. The fractions of R inputs with N=0 and N=1 are respectively

```math
\begin{aligned}
r_0/L&=(1-1/(2M))^M\ge1/2,\\
r_1/L&=\tfrac12(1-1/(2M))^{M-1}>1/4.
\end{aligned}                                      \qquad(10)
```

The inequalities follow directly from Bernoulli's inequality. Including
the arbitrary incoming bit z doubles both cardinalities. There are
therefore enough N=1 states to fill the target set {z=1,b=1}, of size
L/2, and enough N=0 states to fill {z=1,b=0}, also of size L/2.

Refine W as follows. Choose L/2 basis states from each of those two
input classes, map them bijectively to the respective target sets, and
map all remaining odd/even g classes to the remaining E'/complement
sets inside z=0. Ordering each set lexicographically specifies a
permutation. This is compatible with (6) and has the same resource
ledger. Section 5.1 supplies its exact native implementation, with
exponential cost in the packing register's size.

Exact sector preservation implies more than a norm estimate: S_j
preserves every A_j,t_j basis value and is literally identity if
a_j(R)=0. On a_j(R)=1 it has a private unitary fiber within 2η of
(-1)^F. Thus the actual batch is block diagonal in R, and a monomial W
keeps the packed word block diagonal in R,z. On z=1,b=0 its fiber is
I; on z=1,b=1 it is exactly one selective source fiber. On z=0 its
fiber may involve many copies, but C_0 is identity there, so E_F
cancels that fiber with its actual adjoint exactly.

Consequently E_F is exactly I on z=0. On z=1 its two b fibers are
one selective source fiber and its adjoint, each within 2η of (-1)^F.
The X_z echo encounters only one nonidentity E_F fiber on each z input.
It follows that the complete final query satisfies

```math
\|Q_{\rm final}-Q_f\otimes I_{\rm private}\|
 \le2\eta.                                         \qquad(11)
```

R,z return exactly for an exact monomial W; all possible old-output and
source-private errors are included in (11). Full-norm approximation of
this W by mu adds 8mu and includes its return errors. The general
unitary packing of Section 3 does not inherit (11): the direct-sum
argument uses the specified monomial preimages. General source leakage
still has the bound 32h epsilon from (9).

For exactly flat compression, the retained recovery gives
η=sqrt(2)epsilon and exact source-sector invariance. The refined bound
is therefore 2sqrt(2)epsilon, with the recovery flag included in k.
It needs only constant source accuracy at each fixed final tolerance.

At k=Theta(log b), h=2^k is polynomial in the source address length b.
Paying this refined W with polynomial simultaneous count/width overhead
and d_W=polylog(k) would remove the exponential-in-k depth multiplier
without a varying-accuracy source requirement in the exact-sector case.
The native-sharing hypothesis and simultaneous source budgets still
apply. General leakage instead retains the precision schedule above.
The [root-fusion error multiplier](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
must be included in either application.

## 5. Lexical degree and exact small-register synthesis

For the particular E' in (6), put

```math
K=2r-L/2=2^{d-1}-2^h(M-1)^h.
```

For k≥1 and h=M, the power of two dividing K is exactly 2^h. Hence
H(R)=[rank(R)<K] ignores the lowest h rank bits. On its remaining
d−h=hk bits it has odd truth-table weight K/2^h, so its full-degree
algebraic-normal-form coefficient is nonzero. Its Boolean degree is
exactly hk. The E' membership predicate is

```math
e'(R,z)=zb+(1+z)H(R),\qquad \deg e'=hk+1.            \qquad(12)
```

By contrast deg g=k+1. If W is monomial, its phases cancel when
conjugating Π_E, and g composed with its inverse basis permutation
equals e'. Thus some inverse coordinate of W has degree at least
(hk+1)/(k+1). If W is composed of ell monomial rounds, each with
inverse coordinate degree at most c>1, then

```math
\ell\ge\log_c\!\left(\frac{hk+1}{k+1}\right)
      =\Omega(\log h)=\Omega(k).                    \qquad(13)
```

Arbitrary additional dirty helpers do not change this composition
bound when the final predicate is helper-independent. This stops a
bounded-degree monomial implementation of this particular lexical
packing from furnishing polylog(k) rounds. It is not an unrestricted
native T-depth lower bound. A different E', a different echo, or
genuinely nonmonomial packing is outside this argument. Generic
truth-table permutation synthesis is also not a polynomial-count
receipt in d=h(k+1) itself. Its exponential cost can nevertheless be
affordable relative to the source address length, as follows.

### 5.1 Paying any sufficiently small packing permutation

Every permutation of N bits is a composition of at most 2N−1
single-target Boolean gates. Such a gate flips one bit by an arbitrary
Boolean function of the other N−1 bits. The target sequence can always
be chosen as N,N−1,...,1,...,N−1,N.

This bound is standard: Soeken, Abdessaied and Drechsler,
[Section 3, Theorem 1 and Corollary 1](https://arxiv.org/pdf/1407.5878),
give the Young-subgroup formulation, crediting De Vos and Van Rentergem.
The following pair-graph proof makes the recursive controls explicit;
the dirty-bank native payment is separate.

Here is a constructive induction. Write input and output as (z,u) and
(w,v), with u,v the final bit. Form a bipartite multigraph with one
left vertex per z and one right vertex per w, and an edge for each
permutation input. Every vertex has degree two. Alternately color the
edges of each even cycle, including a two-edge cycle, with colors 0,1.
The edge color has the form u+a(z) at its input and v+b(w) at its
output. Let A flip u by a(z), and B flip v by b(w). Then B pi A
preserves the final bit: the edge entering pi after A has the original
u as its color, and B makes its output bit equal that color. Apply
induction to the two (N−1)-bit restrictions of B pi A.
Their common target sequence permits combining each pair of Boolean
controls into one function also controlled by the final bit. The
one-bit base case is I or X. This proves the stated count and supplies
a terminating classical construction from the permutation table.
The refined W table itself is obtained by enumerating its 2^N basis
labels and ordering the specified classes. No rank/unrank circuit or
initialized table register is assumed in the emitted quantum word.

Each single-target gate has an exact zero-clean implementation. Put
m=N−1 and use the
[retained exact dirty indicator](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
I_m on the other data bits x, with an arbitrary 2^m-bit output bank B
and arbitrary private work. Thus B becomes B+e_x. For the desired
Boolean function f define the Clifford parity operation

```math
P_f:\ t\longmapsto t+\sum_{z:f(z)=1}B_z.
```

The chronological word I_m,P_f,I_m†,P_f sends t to t+f(x) and returns
B and every private input literally. The two unknown bank parities
cancel. This is a full operator identity with arbitrary old target,
bank, work and reference; it does not compute f into an initialized
register. Each P_f costs at most 2^m CNOTs.

Let L_N=2N−1 and let b_I exclude the indicator's specified x,B registers.
One returned pool can be reused for all single-target gates. Common
forward/actual-inverse bounds are

```math
\begin{aligned}
T_W&\le2L_NT_I,& d_W&\le2L_Nd_I,\\
G_W&\le2L_NG_I+2L_N2^m,& b_W&\le2^m+b_I.
\end{aligned}                                      \qquad(13a)
```

The retained indicator is exact and has safe majorants
T_I,G_I,b_I≤C 2^m(m+2)^3 and d_I≤C log_2(m+2). Consequently

```math
\begin{aligned}
T_W,G_W&=O\!\left(N2^N(N+2)^3\right),\\
b_W&=O\!\left(2^N(N+2)^3\right),\\
d_W&=O\!\left(N\log(N+1)\right).
\end{aligned}                                      \qquad(13b)
```

This synthesis is phase-exact. Its actual inverse reverses the emitted
native word and has the same bounds. It applies in particular to the
refined monomial W of Section 4.1, so its complete error remains 2η.
There is no W-approximation contribution.

For this packing N=h(k+1)+1 and h=2^k. Substitution in (8) gives

```math
d_{\rm final}=O\!\left(d_V+2^k(k+1)^2\right).       \qquad(13c)
```

Let p now denote the source address length, with q=2^p. Assume the
source has the relaxed simultaneous bounds T_V,b_V≤2^p poly(p),
G_V≤4^p poly(p), and the native-sharing contract of Section 1. If
2^k(k+1)^2=o(log(p+2)), then N log(N+1)=o(log(p+2)); the added
packing count and width are p^o(1), and the h private source copies
retain polynomial address overhead. Any d_V=o(log(p+2)) is then
sufficient for sublogarithmic promoted depth.

For example, take k=floor(beta log_2 log_2(p+4)) with beta<1 and
d_V=O((log(p+2))^gamma) with gamma<1. These two inequalities are
independent; their sum need not be less than one. This improves the
sufficient source window of the
[serial promotion recurrence](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#6-optional-clean-source-windows-and-the-remaining-operation).
For exact source-sector preservation, the refined packing needs only
fixed source accuracy at fixed requested tolerance. General leakage
still requires the uniform accuracy schedule in (9).

This is a paid wrapper improvement, conditional on an actual source
with the stated native sharing, full return, counts and precision.
The [root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) removes
polynomial address overhead and multiplies the complete source error
by eight. It does not provide the missing clean source.

### 5.2 Complete substitution for a compressed source

Suppose U has k_0 clean bits and compression K within epsilon≤1/8 of
Q_f/sqrt(2). All old output and other work inputs remain arbitrary.
Use the
[retained recovery with one extra flag](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag),
including its literal scalar correction. It gives k=k_0+1 clean bits,
h=2^k, and

```math
\begin{aligned}
T_V&\le3T_U+4T_{C_{k_0}},\\
G_V&\le3G_U+4G_{C_{k_0}}+42,\\
d_V&\le3d_U+4d_{C_{k_0}},\\
b_V&\le b_U+b_{C_{k_0}}.
\end{aligned}
```

Here C_j is the separate returned j-bit zero-test. Require U's completed
native Clifford boundaries to fix the address and its T targets to be
private to each copy. The recovery's extra flag and zero-test pool are
private, so the recovered V retains the sharing contract.

Combining this recovery, the parallel selective phases, the exact W of
(13a)–(13b), and the final echoes gives one zero-clean emitted word with

```math
\begin{aligned}
T_{\rm final}\le{}&h\bigl(48T_U+64T_{C_{k_0}}+16T_{C_k}
                         +512q\bigr)+8T_W,\\
G_{\rm final}\le{}&h\bigl(48G_U+64G_{C_{k_0}}+16G_{C_k}
                         +704q+672\bigr)\\
                  &+8G_W+12+2q,\\
d_{\rm final}\le{}&48d_U+64d_{C_{k_0}}+16d_{C_k}+512+8d_W,\\
b_{\rm final}\le{}&h\bigl(k_0+3q+3+b_U+b_{C_{k_0}}+b_{C_k}\bigr)
                    +1+b_W.
\end{aligned}                                      \qquad(13d)
```

Counts, width, error and depth refer to this same native word and its
actual reverse. Width excludes the specified address and output. All
arbitrary borrowed-register and reference errors are included in

```math
\delta_{\rm final}\le
\begin{cases}
2\sqrt2\epsilon,&K^\dagger K=I/2,\\
32h(3\sqrt2\epsilon+2\epsilon^2),&\text{general compression}.
\end{cases}                                         \qquad(13e)
```

Exact flatness makes the recovered clean sector invariant and permits
the refined constant-accuracy cancellation; it is a separate identity
to prove. General compression instead requires uniform source budgets
at the displayed shrinking tolerance. After root fusion the two errors
are at most 16sqrt(2)epsilon and 256h(3sqrt(2)epsilon+2epsilon^2).
For a final tolerance tau≤1, sufficient respective choices are
epsilon≤tau/(16sqrt(2)) and epsilon≤tau/(1152h), also with epsilon≤1/8.
The flat case returns its packed R,z registers exactly; possible
logical and source-private errors remain in the full norm. The general
case makes no stronger arbitrary-input return claim than (13e).

## 6. Alternative: a constant spectral gap permits a paid scalar filter

This is an alternative to packing W, not a way to implement it. Suppose
there are s fixed global Clifford frames C_i on R, with s independent
of k, for which

```math
\begin{aligned}
H_1&=\sum_{i=1}^s\frac{I-C_i(-1)^gC_i^\dagger}{2},\\
\mathrm{spec}(H_1)&\subseteq[\gamma,s],\qquad\gamma>0.
\end{aligned}                                          \qquad(14)
```

The lower bound is on the entire arbitrary R space. No initialized
state, average-state gap or restricted spectral support suffices.
The following construction pays extraction from (14); it does not
provide frames satisfying (14).

### 6.1 Controlled phase batches and constant-clean simulation

Introduce an arbitrary query-output bit u and two clean signal bits
a,b. The desired intermediate complete query adds F(x,Z) to u. For
controlled G_{uF}, replace each selective zero-test by

```math
C_{a,u}:\ e_j\longmapsto e_j+a u[A_j=0^k].
```

When au=0 both echoes are identity; when au=1 they are the original
echoes. Thus the emitted word is literally controlled-a G_{uF}; it
never controls V. The
[returned-conjunction compiler](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
on k+2
literals fixes a,u at every actual native boundary and has only private
T targets. The quartic K receipt (2) is unchanged. If T_D,G_D,d_D,b_D
are zero-test bounds at k+2 literals, safe batch bounds are

```math
\begin{aligned}
T_g&\le h(4T_V+4T_D+128q),\\
G_g&\le h(4G_V+4G_D+176q),\\
d_g&\le4d_V+4d_D+128,\\
w_g&\le h(k+3q+2+b_V+b_D).
\end{aligned}                                      \qquad(15)
```

Here a,b are the only clean wires, u is arbitrary, and w_g counts all
other borrowed wires outside x,Z. The controlled-batch error e_g is at
most 2hη for exact source-sector preservation, or 8h epsilon otherwise.
Conjugating by the uncontrolled C_i and C_i† gives controlled G_{uF,i}.

Let P_i=(I−G_{uF,i})/2. For an exact Hermitian G_{uF,i}, the word

```math
E_i=H_a\,\mathrm{ctrl}_a(G_{uF,i})\,H_a
```

computes its eigenvalue bit into a initially zero. Applying
CP_{b,a}(−theta), then the actual inverse E_i†, gives controlled-b
exp(−i theta P_i) and returns a. A literal CP implementation uses
three R_z rotations and two CNOTs. With
R_z(phi)=exp(−i phi Z/2), those rotations are at angles
−theta/2,−theta/2,+theta/2 on b,a,b xor a. The emitted operation is
e^{i theta/4} CP(−theta). This fixed scalar is retained; it cancels
in the final conjugation below. No scalar correction or extra clean
phase wire is assumed.

Take t=pi/(2s). On a=0, a first-order product formula with L steps and
B=sL terms simulates controlled exp(−itH_{uF}). Its literal full
extension on a=1 instead simulates controlled exp(−it(sI−H_{uF})),
because the computed eigenvalue bit is complemented there. Include the
common retained CP scalar on both blocks. This is not an identity
extension on arbitrary a. The same Trotter bound holds on both blocks,
since complementing every P_i leaves all commutator norms unchanged.
A safe full-norm error against this block extension is
delta_T≤t^2 s^2/(2L). If each of the three single-qubit rotations per
term is approximated in full norm by mu_z, the resulting one-test word
A, with H_b endpoints, has full-norm error at most

```math
\delta_A\le\delta_T+B(2e_g+3\mu_z).                 \qquad(16)
```

All inverses are actual native reverses. On the ideal a=0 code, A
preserves each H_{uF} eigenspace and acts on b as a two-dimensional
unitary. At uF=0 it is identity up to the retained fixed scalar. At
uF=1 its b=0 failure amplitude obeys

```math
 |A_{00}|\le c=\cos\!\left(\frac{\pi\gamma}{4s}\right)\lt1.
```

### 6.2 Fixed-point cubing with the same two clean signals

The pi/3 fixed-point step is the one introduced by
[Grover, *A different kind of quantum search* (2005), Sections 2–4,
equation (3)](https://arxiv.org/pdf/quant-ph/0503205). Here it is used
with determinant-one signal rotations, explicit arbitrary-bank return,
and a complete native precision ledger.

Put R=diag(exp(i pi/6),exp(−i pi/6)) on b and define

```math
S_0=A,\qquad S_{j+1}=S_j R S_j^\dagger R^\dagger S_j.
                                                               \qquad(17)
```

For any two-dimensional unitary S with failure amplitude a=S_00, the
new failure amplitude is a|a|^2. Indeed let omega=exp(i pi/3),
R_0=diag(omega,1) and R_1=diag(1,omega). The 00 entry of
S R_0 S† R_1 S is

```math
a\left[1+(\omega-1)\omega
       +(\omega-1)(1-\omega)|a|^2\right]
 =\omega a|a|^2.
```

Since R_0=omega^(1/2)R and R_1=omega^(1/2)R†, (17) has the claimed
entry, without an unimplemented scalar phase. The uF=0 block remains
scalar identity, and its scalar cancels in S_j† Z_b S_j.

For N=3^j, S_j has N A/actual-adjoint appearances and N−1 R/actual-
adjoint appearances, on the same two clean wires a,b. Therefore
V_F=H_u S_j† Z_b S_j H_u is a clean-code query for F with error

```math
\epsilon_F\le2c^N+2N\delta_A+2(N-1)\mu_R,           \qquad(18)
```

where mu_R bounds full-norm synthesis error of R. The error includes
both clean returns and every arbitrary bank/reference input. The
ideal-error term follows from
norm((S_j†Z_bS_j+I)|0_b>)=2|(S_j)_00| in the uF=1 block.

Unlike using j independent test qubits, the clean count in (18) is
always two. The [retained promotion recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
at k=2 consequently gives a zero-clean
query for F with error at most 32epsilon_F. Two such queries in the
Z_u echo give phase (−1)^F with error at most 64epsilon_F. The final
H_Z envelope gives the requested complete XOR query into Z. Every
borrowed register is covered by this full-operator return estimate.

For any fixed final tolerance tau, choose the power of three N large
enough that 2c^N≤tau/128, then choose delta_A≤tau/(512N) and
mu_R≤tau/(512N). These choices make (18) at most tau/64. Choose L,
rotation synthesis accuracy and source accuracy in (16) accordingly.
There is no growing-clean accuracy obstruction: N,L and all rotation
budgets depend only on s,gamma,tau. Source error still requires
e_g=O_tau(1), hence η=O_tau(1/h) or epsilon=O_tau(1/h), in contrast
to the refined monomial packing's constant-accuracy bound (11).

### 6.3 Native ledger and the remaining operation

Let C_* bound the elementary Clifford count of every C_i and its actual
inverse. Let (T_z,G_z,d_z) and (T_R,G_R,d_R) be common native word
bounds for the chosen single-qubit rotation approximants and actual
inverses. These are finite fixed-precision Clifford+T words, not free
rotations. Safe common bounds are

```math
\begin{aligned}
T_A&\le B(2T_g+3T_z),\\
G_A&\le B(2G_g+4C_*+6+3G_z)+2,\\
d_A&\le B(2d_g+3d_z),\\
T_F&\le2NT_A+2(N-1)T_R,\\
G_F&\le2NG_A+2(N-1)G_R+4,\\
d_F&\le2Nd_A+2(N-1)d_R.
\end{aligned}                                      \qquad(19)
```

The final four Cliffords in G_F pay Z_b=SS and the two H_u gates.
The two signal bits are excluded from w_g. Let C_2 be the retained
two-bit zero-test used by the outer promotion. Substituting k=2,q=1
in the retained promotion ledger, then paying the two-query Z_u echo
and final H_Z envelope, gives

```math
\begin{aligned}
T_{\rm spectral}&\le32T_F+32T_{C_2}+224,\\
G_{\rm spectral}&\le32G_F+32G_{C_2}+456+2q,\\
d_{\rm spectral}&\le32d_F+32d_{C_2}+128,\\
w_{\rm spectral}&\le w_g+b_{C_2}+5.
\end{aligned}                                      \qquad(20)
```

The helper count excludes x and specified Z, and includes the arbitrary
final echo bit u. The word has zero initialized work. Error is at most
64epsilon_F, also for the actual inverse. Counts, width, accuracy and
depth refer to this same word.

Thus a constant number of explicitly paid global frames satisfying
(14) would replace W by a constant-clean filter and its paid promotion.
No such Clifford-frame family is supplied here. The later
[MSS construction](MSS_SPECTRAL_PROMOTION.md) instead pays global native
frames and adds their T cost to this filter's ledger. Local product frames do
not establish (14), and no average-state or random-state argument can
replace its whole-space gap. Both this spectral alternative and the
packing alternative also require an actual qualifying clean source V,
which is still absent. Paying a wrapper alone does not improve the
complete unrestricted frontier.

[Structured spectral-frame limits](SPECTRAL_PROMOTION_FRAME_LIMITS.md)
give three separate obstructions using this tensor reflection: frames
sharing a Lagrangian Pauli space, bounded physical block cones, and
explicit Bell orbits for global graph frames. The block-cone statement
concerns ordinary block depth, not native T-depth. Dense frames with
different noncommuting global Pauli bases remain open.

## 7. Structured packing: concentration is paid, radix compression is not

The generic permutation compiler in Section 5.1 is sufficient below its
stated window. A structured alternative separates the useful component
from its remaining native operation. Let $`k\ge1`$, $`h=2^k`$,
$`n=hk`$, and retain arbitrary k-bit keys A_j and weights t_j. Put

```math
q_j=[A_j=0],\qquad g=\bigoplus_jq_jt_j,\qquad
Q=\bigvee_jq_j.
```

A balanced tree concentrates g into one weight t_*. Each subtree keeps
one of its actual weight bits. At a node, let u,v be the nonemptiness
predicates of the unchanged left/right key registers. Swap the two
representatives when $`(u,v)=(0,1)`$, and add the right representative
to the left when $`(u,v)=(1,1)`$. These disjoint conditions transform
the old contribution $`ut_L+vt_R`$ into $`(u\lor v)t_L`$.
For each fixed key input this is an invertible linear map on arbitrary
weights. The complete tree U therefore satisfies

```math
U(-1)^{Fg}U^\dagger=(-1)^{FQt_*}.
```

F is any external Boolean phase label. No predicate or weight is
initialized. A controlled SWAP is three single-target Boolean XORs;
the conditional CNOT is one. At a node with s leaves, each is paid by
the exact dirty-indicator parity echo on at most ks+1 controls, in
$`O(\log(ks+2))`$ T-depth and exponential count/dirty width. Distinct
nodes at one level have disjoint principal inputs and private pools.
Thus the same literal word and actual inverse obey

```math
d_U=O(k\log(n+2)),\qquad
T_U,G_U,b_U\le2^{O(n)}\mathrm{poly}(n).
```

The all-nonzero key set B has size
$`(2^k-1)^h\lt2^{n-1}`$. Suppose a zero-clean key permutation C
maps B into the halfspace a=0 for one designated actual key bit a.
With $`Q'=Q\circ C^{-1}`$, this would imply $`aQ'=a`$. Define

```math
\begin{aligned}
D&=(CU)(-1)^{Fg}(CU)^\dagger=(-1)^{FQ't_*},\\
L&=\mathrm{CNOT}(a,t_*),\\
E&=D^\dagger L^\dagger DL=(-1)^{Fa}.
\end{aligned}
```

Then $`E^\dagger X_aEX_a=(-1)^F`$ on the entire space. This uses
four original phase batches and eight CU/actual-inverse occurrences,
returning all arbitrary keys, weights and helpers. It pays extraction
conditional on C. It does not implement C.

This concentration also does not retain Section 4's one-active-fiber
error cancellation: a Q=1 sector can contain many active source copies.
The complete batch error and four-use telescoping bound remain necessary.
The refined $`2\eta`$ bound cannot simply be transferred to this word.

### 7.1 An elementary shuffle and the exact remaining operation

Let $`L=2^r`$ and $`0\le b\le L`$. The stable exchange of adjacent
integer intervals [b,L) and [L,L+b) is exactly three additions on an
(r+1)-bit register: add minus b modulo 2L; if the high bit is zero,
add b modulo L to the low bits; then add b modulo 2L. The middle
rotation, translated by the outer additions, rotates [b,L+b) and fixes
its complement. It maps the first interval forward by b and the second
backward by L minus b.

[Remaud and Vandaele, Sections 5--7](https://arxiv.org/html/2501.16802v2)
pay exact reversible addition without initialized work and constant
addition in $`O(\log^3(r+2))`$ Toffoli depth with at most one dirty
helper. A stable outer control is added to their constant-arity gates
using returned private helpers. Hence the elementary interval shuffle
has polylogarithmic native depth and polynomial count/width. Neither
measurement nor clean carry history is used.

For an m-bit row x and r-bit column y, the desired multi-bit stable
radix shuffle is the total permutation

```math
F_b(x,y)=
\begin{cases}
xb+y,&y\lt b,\\
2^m b+x(L-b)+(y-b),&y\ge b.
\end{cases}
```

Its first rectangle occupies the first $`2^m b`$ integers; the other
rectangle fills the complement. Successive row-bit shuffles implement
it, but use m dependent stages. At the top of a balanced packing of
the nonzero-radix digits, m is proportional to n. This retains linear
depth in n, up to logarithmic factors.

The first unpaid operation is a sublinear-depth full-space implementation
of this multi-bit radix permutation, or another actual key compression C.
Out-of-place XOR queries for a map and its inverse do not provide its
in-place action on an arbitrary target bank. The cited adder pays each
elementary shuffle, not removal of their serial dependence. The new
[approximate zero reflection](EXPANDER_ZERO_REFLECTION.md) likewise
does not supply this permutation. Neither composition changes the paid
packing window or the unrestricted source frontier.

The [MSS spectral wrapper](MSS_SPECTRAL_PROMOTION.md) subsequently
bypasses this permutation and pays a larger conditional source window.
The first unpaid radix operation above remains a limitation of the
structured packing mechanism, rather than a barrier to promotion by
other native dynamics.
