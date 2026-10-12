# Ordered block banks and finite histories

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Residual sources](ORIGINAL_TREE_RESIDUAL.md) · [Disjoint-block primitive](DISJOINT_BLOCK_COMPILER.md) · [Research map](../README.md)

Every prescribed complete real Hopf frame has an exact consecutive-band
factorization with linear total squared block size and logarithmically many
ordered factors. Two native banks store all these derived blocks at linear
T-count. A nested-support history construction supplies a dimensionally legal
source pair and an exact native right readout. Coarse prefix banks and fine
local SU(2) banks have separately charged native implementations.

None of these statements implements the ordered fine action at linear total
T-count. The fine left source, or an equivalent joint native word, remains
unproved. The modular separation below concerns unchanged, independently
compiled banks; it is not an unrestricted T-count lower bound.

## 1. Complete-operator contract and exact depth bands

Let $`N=2^n`$, $`n\ge3`$, and $`\eta=2^{-N}`$. Use the
[prescribed operator](../../docs/HOPF_INTERFACE.md), with
$`R_y(\theta)=e^{-i\theta Y}`$, every supplied signed or singular tuple,
two clean flags and at most $`N+n+7`$ arbitrary dirty qubits. Every claimed
native approximation uses the full initialized-isometry norm, literal phase,
actual inverses and dirty/reference return. Only T-count is optimized.

For $`0\le a<b\le n`$, define $`F_{a,b}=L_{b-1}\cdots L_a`$.
This is a direct sum of $`2^a`$ special-unitary blocks of dimension
$`2^{b-a}`$, acting on

```math
\{(p2^{b-a}+u)2^{n-b}:0\le u<2^{b-a}\},\qquad 0\le p<2^a,
```

with identity on all other labels. It follows directly by retaining the
prefix and zero suffix in each defining layer. This is an all-column
identity, not a zero-input promise. The active-mode and squared-size sums
are

```math
M_{a,b}=2^b,\qquad B_{a,b}=2^{2b-a}. \tag{1}
```

Any consecutive depth partition gives W exactly, in its original order.
Write $`\mathcal M=\sum M_{a,b}`$ and $`\mathcal V=\sum B_{a,b}`$.
Neither quantity is a T-count or a description-length lower bound.

## 2. Linear volume with logarithmically many factors

Choose band heights from the bottom upward as
$`3,2,3,5,9,17,\ldots`$, truncating only the final band at the root.
For $`j\ge2`$, the height is $`2^{j-2}+1`$, and the untruncated
cumulative height after j bands is $`t_j=2^{j-1}+j+1`$.
The bottom volume is $`8N`$; band j has volume at most
$`N2^{1-j}`$ for $`j\ge2`$. Hence

```math
\mathcal V<9N,\qquad K\le1+\lceil\log_2n\rceil. \tag{2}
```

The first earlier factor occupies at most $`N/8`$ modes; each subsequent
one occupies at most a quarter of the previous. Thus earlier factors together
have active-mode sum at most $`N/6`$, and volume below N.
In particular $`\mathcal M\le7N/6`$.

### A nearly optimal volume schedule

Keep the final p depths as individual factors, and use (2) for the remaining
$`n-p`$ depths. For $`0\le p<n`$, the singletons contribute
$`4N(1-2^{-p})`$ in volume and $`2N(1-2^{-p})`$ in modes. The
remaining problem has size $`N'=N2^{-p}`$. Consequently

```math
\mathcal V<(4+5\,2^{-p})N,\qquad
\mathcal M\le2N-\tfrac56N2^{-p},\qquad
K\le p+1+\lceil\log_2(n-p)\rceil. \tag{3}
```

At p=n the exact sums are $`4N-4`$ and $`2N-2`$. For fixed
$`\epsilon>0`$, fixed p gives volume below $`(4+\epsilon)N`$ and
$`K=\log_2n+O_\epsilon(1)`$. With
$`p=\lceil\log_2\lceil\log_2n\rceil\rceil`$, the representation
has $`\mathcal V=(4+o(1))N`$ and
$`K=(1+o(1))\log_2n`$. The lower bounds next justify both equalities.

### Sharpness in the literal block/monomial model

Consider K direct-sum block stages separated only by permutations and
unitary diagonals. Blocks, supports and interlayers may depend on the target.
Count the squared dimensions of nontrivial blocks as $`B_j`$ at stage j.
No free nonmonomial mixer is allowed in this representation model.

Starting from one basis mode, let $`s_j`$ be the number of reachable modes.
A touched m-block adds at most $`m-1`$ modes and
$`m^2\ge4(m-1)`$. Thus $`\mathcal V\ge4(s_K-1)`$.
The balanced Hopf root column is flat. At endpoint error, omitting even one
coordinate costs $`1/\sqrt N>\eta`$, so

```math
\mathcal V\ge4(N-1). \tag{4}
```

For consecutive bands this also follows by telescoping
$`2^{a+2h}-4(2^{a+h}-2^a)=2^a(2^h-2)^2\ge0`$.
Height-one layers attain equality.

If s previous modes touch u blocks of sizes $`d_1,\ldots,d_u`$ and
q unchanged singletons, then $`q+u\le s`$. Cauchy--Schwarz gives

```math
s_j^2\le s_{j-1}(B_j+s_{j-1}),\qquad
s_{j-1}\le1+\sum_{i<j}B_i,\qquad
s_K\le(\mathcal V+1)^{1-2^{-K}}. \tag{5}
```

At error $`\epsilon<1`$ the flat column needs at least
$`(1-\epsilon^2)N`$ reachable outputs; therefore

```math
\mathcal V\ge[(1-\epsilon^2)N]^{2^K/(2^K-1)}-1. \tag{6}
```

A linear volume budget implies $`K\ge\log_2n-O(1)`$. Conversely,
for fixed K choose endpoints
$`r_j=\lceil n(1-2^{-j})/(1-2^{-K})\rceil`$, omitting repeated
endpoints. Each band exponent $`2r_j-r_{j-1}`$ is at most
$`n2^K/(2^K-1)+2`$, giving
$`\mathcal V\le4K N^{2^K/(2^K-1)}`$ for the entire frame.

These bounds do not restrict arbitrary native circuits. A Hadamard tensor
spreads a basis state at zero T-count and is outside the monomial interlayers.
They neither price nor forbid a joint source using other mixing.

## 3. Two stored native banks, not two chronological calls

In (2), the final factor occupies N modes with volume $`8N`$. All earlier
factors together occupy at most $`N/6`$ modes with volume below N.
Pack their blocks into another N-mode direct sum with identity padding.
For (3) the earlier-mode sum is still below N and total volume below
$`9N`$, so a fixed linear bound applies to both banks.

Each is covered by the [disjoint SU-block theorem](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks):
for fixed $`\kappa,c`$, volume at most $`\kappa N`$ permits
$`O_{\kappa,c}(N)`$ T-count at error $`\eta/c`$ with precisely the
endpoint workspace. Exact packing, controls, fine source calls, inverses,
and return are included in that theorem. A growing c must be substituted
into its explicit precision/width proof, not hidden as a fixed constant.

The packed direct sum does not apply its factors chronologically to the
original input. An ordered evaluator still needs a legal selection,
transport and return mechanism. The bank theorem guarantees initialized
isometry accuracy, not unspecified occupied-flag action. Cached classical
matrices or gate lists are not a single quantum execution.

## 4. One nested-support history lemma

Let $`P_1\le\cdots\le P_K=I`$ be known computational projections.
Let unitary factors $`F_j,G_j`$ act as identity on $`P_j^\perp`$ and
preserve $`P_j`$. They need not preserve every earlier projection.
Write $`\Delta_j=F_j-G_j`$ and suppose
$`\sum_j\|\Delta_j\|\le\delta<1/4`$. In application order, expand

```math
\mathcal W(z)=\prod_{j=K}^{1}(G_j+z\Delta_j)
=\sum_{k=0}^{K}z^k C_k,\qquad
C_0=Q,\quad\sum_k C_k=W. \tag{7}
```

An input outside $`P_j`$ is unchanged in its first j factors. Its formal
column has degree at most $`K-j`$, proving

```math
C_k=C_kP_{K-k+1},\qquad
\|C_k\|\le\delta^k/k!. \tag{8}
```

The norm bound sums products of local defects and uses unit norm of the
intervening G factors; the elementary symmetric sum is bounded by
$`(\sum_j\|\Delta_j\|)^k/k!`$. The formal factors at general z
are not asserted unitary or natively implemented.

For consecutive bands $`[r_{j-1},r_j)`$, $`P_j`$ selects labels
multiple of $`2^{n-r_j}`$, of dimension $`2^{r_j}`$. Group the actual
coarse layers of [the residual corollary](ORIGINAL_TREE_RESIDUAL.md#3-an-exactly-logical-native-coarse-word)
in exactly the same bands. Telescoping inside disjoint blocks gives the
assumed sum of defect norms. Complex native coarse blocks are retained.

For the ungrouped n layers, label i has at most
$`h(0)=n`$, $`h(i)=\min(n,1+v_2(i))`$ active factors; equivalently
$`C_k`$ has right support on multiples of $`2^{k-1}`$. Their dimensions
sum to $`2N-2`$. The finer geometric local bound in that case is

```math
\|C_k\|\le\frac{\delta^k}{\prod_{r=1}^k(2^r-1)}. \tag{9}
```

To obtain (9), bound the finite elementary symmetric sum by that of
$`2^{-1},2^{-2},\ldots`$. Splitting its smallest index gives
$`e_k=e_{k-1}/(2^k-1)`$, $`e_0=1`$. It also follows by the
convergent product $`\prod_{j\ge1}(1+z2^{-j})`$.

### A normalized right source and one failure bank

Let h(i) count projections $`P_j`$ containing i. Use shared rows
$`(i,k)`$, $`0\le k\le h(i)`$, with the first N rows $`(i,0)`$.
Set

```math
\beta_{ki}=\begin{cases}2^{-(k+1)/2},&k<h(i),\\2^{-h(i)/2},&k=h(i),\end{cases}
\qquad S_B^{\rm sh}|i\rangle=\sum_{k=0}^{h(i)}\beta_{ki}|i,k\rangle.
\tag{10}
```

The last two weights agree, and their squared sum is one. Different input
labels have disjoint rows, so this right source is already an isometry.
It depends only on the support schedule, not the fine target coefficients.
There are $`m=N+\sum_j\dim P_j`$ shared modes.

Let the leading left rows be $`L_0=Q^\dagger/\sqrt2`$. For each
$`k\ge1`$ and eligible i use the row
$`C_k[:,i]^\dagger/(2\beta_{ki})`$. Then (8) gives

```math
L^\dagger S_B^{\rm sh}=W/2,\qquad
L^\dagger L=I/2+H,\qquad
0\preceq H,\quad\|H\|\le\sum_{k\ge1}2^{k-1}\|C_k\|^2
\le\frac{\delta^2}{1-2\delta^2}<2\delta^2. \tag{11}
```

A single failure bank completes the sources:

```math
S_A=\begin{pmatrix}L\\Q^\dagger(I-L^\dagger L)^{1/2}\end{pmatrix},
\qquad S_B=\begin{pmatrix}S_B^{\rm sh}\\0\end{pmatrix},
\quad S_A^\dagger S_A=S_B^\dagger S_B=I,
\quad S_A^\dagger S_B=W/2. \tag{12}
```

The nonpadding mode count is $`2N+\sum_j\dim P_j`$: $`4N-2`$
ungrouped, at most $`19N/6`$ under (2), and at most
$`4N-(5/6)N2^{-p}`$ under (3). Pad to the available 4N modes.
These are basis modes, not newly initialized dirty qubits.

The left deficit has eigenvalues between $`1/2-2\delta^2`$ and 1/2,
so finite classical coefficient evaluation is well conditioned. Its quantum
square root and the coefficient rows are not free. At zero defect the left
source is $`Q^\dagger/\sqrt2`$ in the initial and failure banks and
zero elsewhere. The actual Q inverse followed by a fixed matching prepares
that coarse source at linear T-count, and its distance to (12) is at most
$`2\delta`$. This is not the required fine source accuracy.

## 5. Exact native right readout and common routing

Start label i in row $`(i,0)`$. At round k apply
$`R_y(\pi/4)=HZ`$ to each pair $`((i,k),(i,k+1))`$ with h(i)>k.
Each round is a matching; successive equal splits give (10), including its
last unsplit weight. The word is a full unitary on all 4N modes.

Two-level HZ operations on arbitrary pairs are not native primitives. Use
the [exact dirty routing construction](DISJOINT_BLOCK_COMPILER.md#dirty-only-routing-and-boolean-phases):
a permutation of $`M=2^w`$ labels has $`2w-1`$ Beneš Boolean-X stages.
Alternately color the degree-two input-pair/output-pair cycles, route each
color through a half-size network, and recurse. Each Boolean stage has the
[dirty reflection interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
with $`\lambda=2^{\lfloor(w-1)/2\rfloor}`$ banks and
$`\lceil(w-1)/2\rceil`$ selectors. Its T-count is $`O(\sqrt M)`$
and Clifford count $`O(M)`$, with exact full-space dirty return.

Pack one matching into adjacent pairs, apply a selected Z then a selected H,
and invert the permutation. The selectors and routers return before the next
reflection. Controlled H uses two T gates; Toffoli uses seven. These are
literal matrix decompositions, not gates modulo a scalar phase. Thus K rounds
cost

```math
T(B)=O(Kn\sqrt N)=O(N),\qquad G(B)=O(KnN),\qquad
b_B=2^{\lfloor(n+1)/2\rfloor}+\lceil(n+1)/2\rceil\le N+n+7.
\tag{13}
```

Here $`K\le n`$ and $`n^2\sqrt N=O(N)`$ with an absolute
constant. Both flags are ordinary occupied mode wires; no extra clean
helper is needed. The exact native implementation fixes all dirty work on
its full input space, including arbitrary reference correlations.

The shared [native emitter](../../compiler_robust_hopf/native_ordered_banks.py)
contains one Beneš/loader/router/inverse implementation for both this readout
and the prefix bank below. The ungrouped n=3 readout emits 18,492 T/T-dagger
gates and 44,940 elementary gates on five mode wires and six dirty helpers.
This is an unoptimized finite example, not a practical gate-count claim.

## 6. An exactly returned bank of coarse prefixes

Let $`V_d=Q_{d-1}\cdots Q_0`$. Its relevant invariant subspace consists
of multiples of $`2^{n-d-1}`$, of dimension $`2^{d+1}`$; denote the
restriction by $`\bar V_d`$. Pack the restrictions as

```math
\mathcal P_Q=\left(\bigoplus_{d=0}^{n-1}\bar V_d\right)\oplus I_2,
\qquad\dim\mathcal P_Q=2N. \tag{14}
```

At coarse depth $`\ell`$, $`2^\ell`$ local words appear in every
history d>$`\ell`$. Their copies form a matching of size
$`R_\ell=(n-1-\ell)2^\ell`$, with
$`\sum_\ell R_\ell=N-n-1`$.
The words have length $`O(s+n-\ell)`$. Compressed matching tables and
exact packing give

```math
\sum_{\ell=0}^{n-2}(s+n-\ell)\sqrt{R_\ell}
=\sqrt N\sum_{h=2}^{n}(s+h)\sqrt{h-1}\,2^{-h/2}
=O((s+1)\sqrt N).
```

Include the packing/unpacking permutations and borrowed predicates:

```math
\widehat{\mathcal P}_Q=\mathcal P_Q\otimes I_b,\qquad
T=O((s+1)\sqrt N+n^2\sqrt N+n^3)=O(N),\quad G=O(N^{3/2}),
\quad b\le2^{\lfloor n/2\rfloor}+\lceil n/2\rceil+1. \tag{15}
```

The equality is exact. A compressed-prefix echo requires actual local words
$`D=(XP^\dagger XP)^2`$, satisfying $`XDX=D^\dagger`$. The emitter
accepts syntactically certified `EchoProgram` seeds rather than assuming this
identity for arbitrary SU(2) inputs. In high-prefix sectors its chronological
word is $`Y_zG_z^\dagger T_hG_zT_hY_z`$, with
$`Y_z=C_zX`$ and the borrowed toggle $`T_h:z\mapsto z\oplus h`$.
The middle exponent is $`(z\oplus h)-z`$; conjugation by $`X^z`$
removes its sign. Inactive action and borrowed return are exact. A seed word
for P is not itself a supplied arbitrary rotation or uncharged oracle.

## 7. A native fine local bank on occupied signals

Put $`B_{d,p}=Q_{d,p}^\dagger R_y(\theta_{d,p})\in SU(2)`$ and
pack the N−1 blocks with one idle pair on the 2N-mode register:

```math
\mathcal D=\left(\bigoplus_{d,p}B_{d,p}\right)\oplus I_2.
```

One original flag is part of the data-mode register. The other is a borrowed
signal, not assumed initialized. This is legal specifically because the
[real-rotation flag symmetry](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
extends that primitive to the complete signal input space. Fixed Clifford
axis conjugations cover three determinant-one Euler factors. The literal
scalar-phase primitive does **not** have that borrowed-signal guarantee.

For $`n\ge9`$, choose $`q=N+15`$. Self-borrowed whole-word queries
need $`q+1`$ core qubits and no additional selectors. Their required
address condition holds, and $`q+1=N+16\le N+n+7`$. Certified Euler
selection to error $`\eta/128`$ and three rotation calls give

```math
\|\widehat{\mathcal D}-\mathcal D\otimes I_{\rm signal,dirty}\|
<129\,2^{-N-15}+\eta/128<\eta/64,
\qquad T=O(N),\quad G=O(N^2). \tag{16}
```

For $`3\le n\le8`$, direct dirty-table synthesis gives a finite
$`O(N\sqrt N)`$ fallback with $`\sqrt N\le16`$. Its constant
is uniform over targets. Equation (16) is full-unitary approximate return;
it is not the exact-return identity in (15).

With $`A_d=V_d^\dagger Q_d^\dagger L_dV_d`$,
$`\mathcal P_Q^\dagger\mathcal D\mathcal P_Q`$ stores the
nontrivial restrictions of these operators. Exact cancellation gives

```math
Q^\dagger W=A_{n-1}\cdots A_0,\qquad
Q^\dagger C_k=\sum_{d_k>\cdots>d_1}(A_{d_k}-I)\cdots(A_{d_1}-I).
\tag{17}
```

This is a left-residual identity, not a transfer of the original right-cut
rank theorem. It exposes the correlated fine data in a single priced bank;
it does not perform the ordered products in (17).

## 8. A modular-bank separation

Restrict the processor so that its sign-dependent data enter only through
the unchanged direct sums of the derived band matrices. All other encodings
and interlayers are sign-independent. Queries may be controlled, inverted,
or fixed-conjugated; tensor occurrences count separately. Grant arbitrary
workspace and even one oracle containing every band's blocks.

Choose one edge in each band at increasing depth $`d_j`$, forming a path
$`v_0=0`$, $`v_j=v_{j-1}+2^{n-d_j-1}`$. All other angles are zero.
Set

```math
s=\lceil\sqrt N\rceil+6,\quad\delta=2^{-s},\quad
u=2^{-(s+n+2)},\qquad \theta_j=2\arctan(\sigma_j u),\quad\sigma_j=\pm1.
```

The common coarse word Q=I is valid: local error is at most 2u and
$`\|W-I\|\le2Ku<\delta/2`$. This is a permitted coarse choice,
not a statement that a particular synthesis search always selects identity.
Each bank entry has Boolean degree at most one. After R queries every
amplitude has degree at most R, whereas the target path amplitude is

```math
\langle v_K|W_\sigma|0\rangle
=\left(\prod_j\sigma_j\right)a,\qquad
 a=\left(\frac{2u}{1+u^2}\right)^K\ge u^K. \tag{18}
```

For R<K its degree-K Fourier coefficient cannot agree; averaging against
$`\prod_j\sigma_j`$ proves worst-case error at least a.
For both the p=0 and the p(n) schedules, when $`n\ge8`$,
$`K(s+n+2)\le N-n`$, hence $`a\ge2^n\eta`$.
At n=8 use respectively K=3 and K=5, s=22. For n≥9, K≤n and
$`s\le2^{n/2}+7`$ suffice because
$`n2^{n/2}+n^2+10n<2^n`$: it holds at nine and both positive
terms divided by $`2^n`$ decrease thereafter. Thus R≥K is necessary.

If each module has error at most $`\eta/64`$ on its actually used
subspace, an honest hybrid gives ideal-query error at most
$`\eta+R\eta/64<2^n\eta`$ for R<K. The restriction survives.
With initialized-isometry module guarantees the ideal call must meet that
initialization premise; no occupied-flag extension is inferred.

Each nonempty bank contains a two-level rotation with stabilizer eigenvector
$`(|a\rangle+i|b\rangle)/\sqrt2`$ and eigenphase
$`e^{\pm i\theta}`$. Its distance to the nearest eighth root, 1, is
$`2u/\sqrt{1+u^2}>u>\eta`$. The
[positive Pauli-failure gap and scalar lemma](../../docs/FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
then force at least $`2N-1`$ T gates in any separately compiled bank
accurate to $`\eta/64`$. Otherwise all stabilizer failures vanish,
forcing an eighth-root eigenphase. Clifford relabelling identifies the
active pair without changing that argument.

Consequently **additive standalone module accounting** gives

```math
T_{\rm modular}\ge K(2N-1)=\Omega(N\log n). \tag{19}
```

The same complete witness acts only on K+1 labels. Since
$`(K+1)^2\le(n+1)^2\le2N`$, it is one eligible SU(K+1) block
and has an O(N)-T compiler by the disjoint-block theorem. This separates
independent high-precision modules from joint synthesis on that family.
It does not add primitive lower bounds against fused gate words,
parameter-dependent interlayers, reparameterized banks, or different
intermediate precision guarantees. Nor does it make the independently
synthesized strategy an improvement over the existing grouped compiler.

## 9. Related restrictions and exact applicability boundaries

### A generic Cayley polynomial still needs growing degree

For $`f(x)=(1+ix)/(1-ix)`$, $`|x|\le r`$, and any polynomial of
degree p, put $`m=p+1`$, $`h=2r/m`$, $`x_j=-r+jh`$. The m-th
finite difference kills the polynomial, while

```math
\Delta_h^m f(-r)=\frac{2(ih)^m m!}{\prod_{j=0}^m(1-ix_j)}.
```

The absolute sum of its coefficients is $`2^m`$. Uniform error ε
therefore obeys

```math
\epsilon\ge\frac{2m!(r/m)^m}{(1+r^2)^{(m+1)/2}}
\ge\frac{2(r/e)^m}{(1+r^2)^{(m+1)/2}}. \tag{20}
```

At $`r=2^{-s-2}`$, ε≤$`2^{-N}`$ implies
$`p\ge\lfloor N/(s+4)\rfloor`$. Conversely, geometric truncation
of $`(I+Z)(I-Z)^{-1}=I+2\sum_{k\ge1}Z^k`$ has error at most
$`4\,2^{-s(p+1)}`$ for $`\|Z\|\le2^{-s}`$.
Thus its polynomial degree scale is $`\Theta(\sqrt N)`$. This is
not a T-count lower bound, not an admissible polynomial-circuit construction,
and not a restriction on the formal layer polynomial (7), target-dependent
reparameterization, or priced rational methods.

### Small unchanged local coins retain all-order information

Activating every edge on one all-right path with
$`t_d=2^{-(s+n-d+2)}`$, $`\theta_d=2\arctan(\sigma_dt_d)`$,
and Q=I gives the same Fourier argument for the local bank in Section 7.
The path amplitude is at least
$`2^{-E_n}`$, $`E_n=ns+n(n+5)/2\le N-n`$ for n≥8. Hence an
unchanged local bank requires n queries even in the coarse neighborhood.
The n=8 inequality is direct; for n≥9 it follows from the preceding
exponential comparison. The target is itself an eligible SU(n+1) block.
This is the small-angle version of the repository's
[original-angle query restriction](STRUCTURAL_COMPILATION_LIMITS.md#5-original-angle-parallel-coins-need-at-least-n-queries),
not an unrestricted impossibility statement.

### Averaging discards ordering that a scalar filter cannot restore

Let $`U=R_{01}(\theta)`$, $`V=R_{12}(\theta)`$, $`W=VU`$ on
three labels and $`B=(U+V-I)/2`$. For
$`v=(e_0+e_2)/\sqrt2`$, $`w=(e_0-e_2)/\sqrt2`$,

```math
Bv=B^\dagger v=(\cos\theta/2)v,\qquad
w^\dagger Wv=-\sin^2\theta/2.
```

Every element of the algebra generated by B and its adjoint preserves the
span of v. Its error to W is therefore at least $`\sin^2\theta/2`$.
This includes the accepted cubic amplification polynomial. With
$`\theta=\delta/(4N)`$, the lower bound exceeds η for n≥6: use
$`\sin\theta\ge\theta/2`$ and check
$`N>2s+2n+7`$, directly at six and monotonically thereafter.
The example obeys the local small-error budget. It excludes that averaged
accepted-block algebra, not a richer rejected-space encoding retaining order.

## 10. Remaining native theorem and reproducibility

A sufficient new result would implement the left source (12) with
$`T(\widehat A)=O(N)`$, the literal workspace, and
$`\|\widehat AJ_2-S_A\otimes I_b\|\le\eta/8`$.
The exactly returned right readout gives

```math
\widehat U=\widehat A^\dagger\widehat B,\qquad
\|2J_2^\dagger\widehat UJ_2-W\otimes I_b\|\le\eta/4.
```

The [robust normalization-two lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
then recovers the complete frame with two forward overlap calls and one
actual inverse. Only the stated initialized-column approximation of A is
needed; its unused extension is not prescribed. A directly fused overlap or
another complete native representation is equally eligible.

A shrinking input rank does not imply cheap descendant data. With Q=I,
$`C_n=\Delta_{n-1}\cdots\Delta_0`$ has at most two nonzero input
columns, but a root column can contain independent dependence on all N/2
bottom angles in disjoint output pairs. Likewise a constant normalization,
computable square root, or a bank prepared once does not remove the ordered
fine program. No fine left-source emitter or endpoint circuit is supplied.

[Combinatorial and matrix routines](../../compiler_robust_hopf/ordered_block_banks.py)
are separate from the [native readout/prefix emitters](../../compiler_robust_hopf/native_ordered_banks.py).
The [tests](../../tests/test_ordered_block_banks.py) include exact integer
schedule/volume checks, exact rational Fourier and root-series witnesses,
finite complex source identities, literal gate-phase tests, full-matrix
small dirty interpreters, and bounded occupied/entangled-vector readout tests.
They do not replace the all-dimension proofs. NumPy is the only nonstandard
runtime dependency. Regenerate compact native summaries, or optional local
gate lists, with

```sh
python scripts/ordered_banks_native_example.py --n 3
python -m unittest tests.test_ordered_block_banks
```

Large gate lists and duplicate scratch simulators are not maintained as
parallel proof sources. Results A–D and the unrestricted endpoint bounds
are unchanged.
