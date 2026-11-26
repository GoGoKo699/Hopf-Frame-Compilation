# Prescribed scalar signs as a sufficient source for general frame depth

A scalar sign on $`a`$ address qubits is the literal operator
$`D_g|x\rangle=(-1)^{g(x)}|x\rangle`$, tensored with identity on all
arbitrary dirty work. This chapter gives a conditional unrestricted
compiler, an average-case criterion for its source, and two table encodings:
an affine matrix pencil and a fixed self-inverse binary shear. The latter
reduces the source question to scalar entry readout with four calls and no
additional clean work. A rank argument screens a thin-kernel readout of the
former. All use the complete frame convention, including arbitrary references
and actual native inverses. The [unrestricted upper theorem R96](ALL_PRECISION_LINEAR_DEPTH.md)
still gives $`O_L(N/b^2+\log(n+2))`$ without the source hypothesis.

## 1. Scalar source and frame bound

Suppose that, for each arity $`a`$, table $`g`$ and additional dirty
allocation $`B\ge C(a+2)`$, a terminating classical compiler emits an
exact native scalar word and its reverse, using at most $`B`$ dirty
wires and a fixed number of initialized qubits. For nondecreasing
$`d\ge1`$ assume the same word satisfies

```math
 T=O(2^{a/2}+2^a/B),\qquad G=O(2^a),\qquad
 D_T=O(d(a)+2^a/B^2).
```

Then for fixed integer $`L\ge6`$, $`N=2^n`$, and the
[R98](QUERY_DEPTH_FRAME_TRANSFER.md) allocations
$`b\ge17(L+n+7)`$ for real frames or $`b\ge17(L+n+8)`$ for the
prescribed complex magnitude completion, there is a one-clean circuit
of complete initialized-isometry error at most $`2^{-L}`$ with

```math
 T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
 D_T=O_L\!\left(N/b^2+\log\log(n+4)
 +(1+\log^*(n+2))d(n+O_L(1))\right).
```

Thus a constant additive scalar depth gives a doubly logarithmic frame
depth at square-root dirty width. A doubly logarithmic additive scalar
depth gives $`o(\log n)`$ there. For a general $`d=o(\log n)`$ the
product with the number of stages matters. Conversely an arbitrary
$`p`$-bit sign occurs as a complete real frame on $`p+1`$ bits by
setting all earlier angles to zero and the deepest angle to $`\pi g(x)`$.
Both final columns acquire the prescribed sign, and the last logical
wire is returned for every input. [Canonical angle domains](../../docs/HOPF_INTERFACE.md)

**Eliminating the source's initialized work.** Compile the lifted table
$`t g(x)`$ and conjugate by Hadamards on $`t`$ to obtain an XOR.
[Complete XOR promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) removes
its bounded clean supply with constant overhead. For the promoted XOR
$`Q`$, the echo $`Q^\dagger Z_t QZ_t`$ gives $`D_g`$ even
when $`t`$ is occupied. One extra address bit only changes constants.

### Why the output bits fit the scalar budget

The Walsh transform of a complete word XOR is the sign
$`(-1)^{y\cdot F(x)}`$. Every relevant output character therefore
counts as a scalar input. In an early stage of remaining height $`m`$
let $`P`$ be the preserved original prefix of length $`n-m`$, let
$`\ell`$ denote local precision, and put
$`t=\lceil\log_2(m+2)\rceil`$,
$`\chi=t^3+\lceil\log_2(\ell+2)\rceil`$,
$`r=\lceil A\chi\rceil`$, $`u=m-r`$. A short preserved selector
chooses the relevant address sites and the output-character sites in
each row. Route both into the same carrier, call the scalar, and undo
the actual route before the final Hadamards. The paired permutation
wrapper for a selector of length $`s`$ and physical union of $`w`$
sites uses $`O(2^s w)`$ gates and dirty wires at depth
$`O(\log(s+2))`$. No row-by-row chronology is introduced. This
argument uses a single completed word even when its native implementation
temporarily moves the address. [Paired permutation routing, §5](PAIRED_PERMUTATION_ROUTING.md)

The free donor reserve in the fine bins is $`m_{\rm free}=O(t)`$.

| Long position | Address plus output characters outside $`P`$ |
| --- | --- |
| Tagged endpoints, joint rank, residuals, anchors, calibration and signs | One $`u+O(\chi)`$ bank and short metadata. Append the actual $`O(\ell+1)`$ calibration characters. |
| Fine rectangle, membership and labels | In selected bin $`i`$, membership reads $`\alpha_i+q_i\tau_0`$ bits and writes the active $`q_i`$-bit slice. Their sum is the rectangle width $`\alpha_i+q_i+q_i\tau_0\le u+m_{\rm free}`$, with the same short overhead as the fine tag. |
| Compact metadata, first/last matching and augmented coarse tag | The long bank has size $`u+O(\chi)`$ and its indices and metadata add $`O(\chi)`$. The combined reference test already occupies that bank. |
| Local shear and global high erasure | In block $`i`$ the active low and high lengths are $`k_i`$ and $`q_i-k_i`$. The two variable slices total $`q_i`$ in either direction. Summing gives $`u`$ plus short split/seed fields. |
| Matching shifts | Only prefix and short indices are read; the target is one $`u+O(\chi)`$ bank. |

For high erasure, decode the high part from the active low codes and
metadata on reached rows. Choose zero if no preimage exists and on
inactive slices or invalid splits. This specifies a total XOR on
arbitrary work while retaining the literal physical routers and their
inverses on malformed metadata. The compact embedding constrains only
the reached labels and inactive behavior. These choices leave its
closures exact. The row bounds use [R98, §2](QUERY_DEPTH_FRAME_TRANSFER.md),
[global stage, §§2–6](GLOBAL_BORROWED_PRECISION_STAGE.md),
[causal matching, §§3–7](CAUSAL_COMPACT_COIN_DEPTH.md), and
[coin dilation, (15)](COIN_DILATION_DEPTH.md).

**Recursion and budget.** Each stage uses constantly many chronological
long calls, of arity at most $`n-r+C'\chi`$. Set
$`h(m)=\max(1,\log^*(m+2))`$ and
$`\ell(m)=L+C_0+2\lceil\log_2(1+h(m))\rceil`$. Continue below the
usual four-stage stopping point until a fixed $`L`$-dependent height.
Above that cutoff $`r(m)\lt m/2`$, $`r(m)=O_L(\log^3(m+2))`$, and
$`r(r(m))\le\log_2 m`$. There are $`O(1+\log^*n)`$ stages with
bounded multiplicity per $`h`$. The errors $`2^{-\ell(m)}/16`$
sum within the prescribed tolerance. For the short depth,
$`\sum_m O(\ell(m)+\log(\chi(m)+2))`$ is $`O_L(\log\log(n+4))`$
because the log-size terms decrease geometrically above a fixed cutoff, and
the precision sum is $`O_L((1+\log^*n)\log(2+\log^*n))`$.

On $`b\ge\sqrt N/n`$ the available long-query pool is a constant
fraction of $`b`$. Increasing $`A`$ absorbs the polynomial stage
factors and leaves a summable factor $`2^{-4\chi(m)}`$ in each
long volume/count contribution. The auxiliary gates and dirty pools
cost at most $`n2^{O_L(\log^3(n+2))}=N^{o(1)}`$. At bounded
terminal height the entire angle program has $`O_L(1)`$ characters,
so prefix plus program fits $`n+O_L(1)`$ scalar bits. The R90 terminal
and the endpoint and complex phase words fit the same ledger. If
$`b\lt \sqrt N/n`$, $`N/b^2`$ covers R96's logarithmic term; finite
dimensions use R96 at their literal allocations. Stopping after only
four stages could leave a terminal with $`q`$ program bits and
$`2^{n-m+q}`$ scalar entries even though its original word table
had just $`q2^{n-m}`$ bits. The bounded-height terminal is essential
to this proof.

## 2. Density of exact easy signs

Fix one allocation. Let $`\mathcal A`$ be a family of exactly
implementable $`a`$-bit tables at the source budgets above and set
$`M=2^a`$. It suffices that a uniform truth table $`U`$ has expected
relative Hamming distance $`\delta\le M^{-\beta}`$ from
$`\mathcal A`$ for fixed $`\beta>0`$. In particular it suffices that

```math
 \alpha=|\mathcal A|/2^M\ \ge\ 2^{-M^{1-\eta}}
```

for any fixed $`\eta>0`$. For $`X=\mathrm{dist}_H(U,\mathcal A)`$,
changing one truth entry changes $`X`$ by at most one. The bounded
differences inequality gives
$`\alpha=\Pr[X=0]\le\exp(-2(\mathbb E X)^2/M)`$. Hence
$`\delta\le\sqrt{\ln(1/\alpha)/(2M)}=O(M^{-\eta/2})`$.

To program a $`p=a-1`$ bit sign $`f`$, lift it to
$`F(x,t)=t f(x)`$. Approximate a uniform mask table $`U`$ by
shifting its address by a random XOR, taking a nearest easy table,
and undoing the address shift with X gates. The shifted table and
shift are independent, so every individual entry is wrong with
probability $`\delta`$. Do the same for $`U+F`$. For the two
translated answers $`A_i,B_i`$ let $`e_i=A_i+B_i`$. At each
$`x`$, its derivative $`e_i(x,0)+e_i(x,1)`$ differs from $`f(x)`$
with probability at most $`4\delta`$. For $`K`$ independent pairs,
$`K`$ odd, the probability that *some* address gets a wrong majority
is at most $`2^p2^K(4\delta)^{(K+1)/2}`$. A fixed
$`K=O(1+1/\beta)`$ makes this less than one for all sufficiently large
arities. Finite search finds a list correct everywhere; small arities
use exact direct synthesis. Equivalently one can enumerate bounded
native words and verify exact scalar tables by their algebraic entries.

For each pair, conjugate the two easy signs $`D_{A_i}D_{B_i}`$ by
Hadamards on a fresh zero vote bit. This records the derivative and
the scalar $`(-1)^{e_i(x,0)}`$. Phase the majority of the votes
and reverse the actual recording word. All ancillary signs cancel on
every coherent input. This costs $`4K`$ source calls including
adjoints, constant extra clean qubits, $`O_\beta(p+1)`$ auxiliary
gates and constant additional T-depth. Apply the bounded-clean
promotion above for the zero-clean interface. The arity increment
multiplies each volume term by a constant. The easy family can vary
with the allocation. The denominator is the number of all truth tables. Applying the theorem
at sublogarithmic depth requires an easy family meeting that density.

A square raw table can be viewed as a binary crossing matrix. Uniform
invertible matrices have cube density $`\prod_{i=1}^{2^{a/2}}(1-2^{-i})`$,
bounded below by a positive constant. If mean distortion to easy signs
on those matrices is $`e`$, Markov gives constant cube density within
relative distance $`2e`$. Applying the same concentration bound to that
neighborhood gives uniform-table distortion $`\le2e+O(M^{-1/2})`$.
An unused input bit handles odd arities. This is a criterion on truth
tables actually realized by complete circuits.

For the restricted self-inverse matrices in §5, a direct four-call reduction
replaces this density criterion. It works for every table, including singular
crossing matrices, and also has a complete approximate-source guarantee.

## 3. Precision and completed address return

For sources at physical isometry error $`\epsilon`$, including
arbitrary leakage, the same density correction has error at most
$`4K\epsilon`$ by replacing its actual calls and reversals along the
ideal trajectory. This does not by itself provide the decreasing error
budget of the iterated transfer. [The original approximate word
transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) has only constantly
many completed long calls and a different terminal hypothesis.

One sufficient fixed-accuracy alternative requires exact preservation
of every computational address sector at each completed scalar
boundary on all physical clean/dirty sectors. With this premise, the
two-reflection construction below amplifies fixed full-operator accuracy
using logarithmically many actual source calls, one arbitrary signal
qubit and no initialized votes. Write $`H_n=1+\log^*(n+2)`$. At the
same fixed-L allocations and error as §1, its conditional frame bound is

```math
\begin{aligned}
 T&=O_L(\sqrt N+N/b),\qquad G=O_L(N),\\
 D_T&=O_L\!\left(N/b^2+\log\log(n+4)
       +H_n\log(H_n+2)d(n+O_L(1))\right).
\end{aligned}
```

It retains one final clean qubit. A shallow fixed-accuracy scalar source
is still a premise; the unconditional frontier is unchanged.

<a id="pair-reflection-scalar-amplification"></a>

### Amplification using two reflections

**Promote before amplifying.** Let the original scalar source use a fixed
$`c_0`$ clean qubits. For a p-bit table $`f(x)`$, compile the lifted
table $`t f(x)`$. Hadamards on t give an approximate XOR of f into
that old target. Apply [complete XOR promotion, §§3–4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
to a distinct arbitrary output z, then conjugate the promoted query by
Hadamards on z. Rename z as t. This gives a zero-clean physical unitary V
approximating $`D_{tf}`$ on its entire address/signal/work space. The old
target, former clean bank and every promotion helper are arbitrary dirty
work. The promotion error is at most $`8\cdot2^{c_0}`$ times the
original source error; thus a fixed original tolerance suffices. Its
count, width and depth overhead is constant because $`c_0`$ is fixed.
Only one lifted scalar address bit is used in the original call.

The original source preserves x on all physical sectors, as do its
actual inverse and all other promotion gates. Hence V exactly preserves
x even off the original initialized code. It need not preserve the new
signal t. This establishes the physical premise used below, rather than
inferring it from an initialized-isometry estimate. A zero-clean source
already satisfying this premise can be used directly.

Let $`0\lt\epsilon_0\le10^{-3}`$ be the fixed promoted error. Allow
a common source phase $`e^{i\theta_f}`$, independent of x, t and work,
and suppose

```math
 \|V-e^{i\theta_f}(D_{tf}\otimes I_{\mathrm{work}})\|_{\mathrm{op}}
 \le\epsilon_0.
```

Set

```math
 S=X_t,\qquad R=V^\dagger S V,\qquad W=RS,\qquad
 \mathcal H=\tfrac12(W+W^\dagger),\qquad c=1-2\epsilon_0^2.
```

The operators S and R are exact Hermitian involutions, and
$`SWS=W^\dagger`$. At fixed x, put $`\sigma_x=(-1)^{f(x)}`$.
The ideal commutator is $`\sigma_x I`$, so the two actual source
appearances give

```math
 \|W_x-\sigma_x I\|\le2\epsilon_0,\qquad
 \sigma_x\mathcal H_x\succeq cI.
```

The second inequality follows by expanding
$`(W_x-\sigma_x I)^\dagger(W_x-\sigma_x I)=2I-2\sigma_x\mathcal H_x`$.
The common source phase cancels exactly. The proof only needs a fixed
$`\epsilon_0\lt1/\sqrt2`$, so that $`c\gt0`$; the displayed
smaller tolerance also covers the retained vote alternative below.

Jordan's decomposition for two reflections applies separately inside
each preserved x-sector. On a two-dimensional block, a mathematical
basis change puts S at Z and W at
$`\lambda I+i\sqrt{1-\lambda^2}X`$, with
$`\lambda\in[-1,1]`$ and $`\sigma_x\lambda\ge c`$.
Each reflection has eigenvalues plus and minus one on such a block,
so W and every phase rotation have determinant one there.
This basis change is used only in the proof. The physical phase rotation
is always $`e^{i\phi S}=e^{i\phi X_t}`$ on the same arbitrary signal.

For odd $`k\ge1`$, define the explicit majority polynomial

```math
 p_k(\lambda)=
 2\sum_{j=(k+1)/2}^{k}\binom{k}{j}
   \left(\frac{1+\lambda}{2}\right)^j
   \left(\frac{1-\lambda}{2}\right)^{k-j}-1.
```

It is a real odd degree-k polynomial, bounded by one in absolute value
on $`[-1,1]`$, with $`p_k(\pm1)=\pm1`$. Viewing the sum as a
binomial majority and applying the bounded-sum tail estimate gives

```math
 1-\sigma p_k(\lambda)\le2e^{-kc^2/2}
 \quad\hbox{whenever}\quad \sigma\lambda\ge c,
 \qquad \sigma\in\{-1,1\}.
```

The real-polynomial quantum signal processing (QSP) theorem of
Gilyén–Su–Low–Wiebe, [F55](../../docs/reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries),
v1 Theorems 3–5 and Corollary 10, pp. 8–14, supplies a common list
$`\phi_0,\ldots,\phi_k`$ for the operator-order word

```math
 U_\Phi=e^{i\phi_0S}
          \prod_{j=1}^{k}\left(W e^{i\phi_jS}\right).
```

The rightmost factor acts first. On each two-dimensional block this is
an SU(2) matrix with upper-left polynomial P satisfying
$`\mathrm{Re}\,P=p_k`$. Every such matrix has Hermitian part
$`(\mathrm{Re}\,P)I`$, so the same physical word satisfies

```math
 \tfrac12(U_\Phi+U_\Phi^\dagger)=p_k(\mathcal H).
```

No circuit extracts a real part or projects onto a chosen block: this
identity controls the whole unitary on every block. The phase list is
common to all addresses and all dirty eigenvectors. On a one-dimensional
block write $`S=s`$, $`R=r`$, where $`s,r\in\{-1,1\}`$, and
$`\lambda=rs\in\{-1,1\}`$. The physical word is
$`\lambda^k e^{is\sum_j\phi_j}`$. QSP at the endpoints gives
$`P(\lambda)=\lambda^k e^{i\sum_j\phi_j}`$, so its Hermitian
part is again $`p_k(\lambda)`$. Since k is odd and
$`p_k(\pm1)=\pm1`$, exact phase factorization gives the literal
endpoint signs on both choices of s.

For a certified finite phase list, it suffices that its realized
polynomial $`\widehat p=\mathrm{Re}\,P`$ obey
$`\|\widehat p-p_k\|_{[-1,1]}\le\delta`$. The constructive
factorization in the cited theorem and corollary permits any positive
delta. An arbitrarily fine certified phase list is sufficient; no
classical runtime bound is needed here. The exact unitary associated
with that list still has Hermitian part $`\widehat p(\mathcal H)`$, including
the one-dimensional blocks. Therefore

```math
 \|U_{\Phi,x}-\sigma_x I\|^2
 =2\|I-\sigma_x\widehat p(\mathcal H_x)\|
 \le4e^{-kc^2/2}+2\delta.
```

For requested $`0\lt\tau\le1/2`$, take
$`\delta=\tau^2/16`$ and the least odd
$`k\ge(4/c^2)\ln(8/\tau)`$. The last bound is at most
$`3\tau^2/16`$, so phase factorization and polynomial approximation
together contribute less than $`\tau/2`$ in operator norm.

**Native phase rotations and complete ledger.** Compile each of the
$`k+1`$ rotations $`e^{i\phi_jX_t}`$ at literal operator error
$`\tau/[2(k+1)]`$ using the phase-calibrated one-qubit synthesis
primitive [F5](../../docs/reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries),
Lemma 2.3. The summed rotation error is at most $`\tau/2`$. Thus the
emitted word has full error at most tau against
$`D_f\otimes I_{t,\mathrm{work}}`$, including every dirty input and
external reference. All source appearances use the same actual V or
its actual adjoint; their original constant error establishes the
spectral gap and is not summed afresh over the QSP repetitions.

Let $`T_V,G_V,d_V`$ be common bounds for V and its actual inverse.
There are exactly $`2k`$ such appearances and $`2k`$ Pauli X gates.
The actual native inverse of the full word is charged in the same
positive-S alphabet. Common forward/inverse bounds are

```math
\begin{aligned}
 T_\tau&\le2kT_V+O\!\left(k\log\frac{k+1}{\tau}\right),\\
 G_\tau&\le2kG_V+O\!\left(k\log\frac{k+1}{\tau}\right),\\
 D_\tau&\le2kd_V+O\!\left(k\log\frac{k+1}{\tau}\right),\\
 w_\tau&\le w_V+1,\qquad a_{\mathrm{clean}}=0.
\end{aligned}
```

Here $`w_V`$ excludes x and the signal t but includes all promoted
work; t becomes the one additional dirty qubit of the final scalar.
One pool is reused sequentially throughout. Since
$`k=O(1+\log(1/\tau))`$, the local synthesis term is
$`O((1+\log(1/\tau))^2)`$. No controlled V, initialized signal,
measurement or postselection is used. This is one degree-k QSP word
factorized from the original pair of reflections; it does not assume
that an amplified cubic word still obeys $`SUS=U^\dagger`$.

### Substitution in the frame recurrence

For each early stage of remaining height m, request
$`\tau_m=c_L/(1+h(m))^2`$, with the function h of §1 and a
sufficiently small fixed $`c_L=O(2^{-L})`$. There are constantly many
chronological source slots per stage and bounded multiplicity per h,
so the added errors sum within a fixed share of the final tolerance.
The amplified degree is $`k_m=O_L(\log(h(m)+2))`$.

The early-stage arity margin $`r-C'\chi`$ in §1 absorbs these
logarithmic call factors in all volume terms after increasing its
fixed constant A. In particular, the inherited factor
$`2^{-4\chi(m)}`$ still bounds summable T-count, Clifford-count
and volume-depth contributions after absorbing $`k_m`$.
The one lifted input and the fixed promotion reservation only change
allocation constants. All other local gates total at most
$`O_L(H_n\log^2(H_n+2))=O_L(\log\log(n+4))`$, with one shared
dirty pool and no new clean supply.

The initial full-address sign words, the bounded-height terminal and
the complex phase words each have constantly many source occurrences.
Assign them fixed L-dependent tolerances, independently of n. Their
amplification therefore has constant overhead, preserving the
$`O_L(\sqrt N+N/b)`$ T-count and $`O_L(N)`$ Clifford bound at
full scalar volume. Shrinking every later tolerance according only to
its chronological stage number would fail this accounting at the
full-volume terminal.

Finally,

```math
 \sum_m k_m d(n+O_L(1))
 =O_L\!\left(H_n\log(H_n+2)d(n+O_L(1))\right).
```

Together with §1's short-stage sum this proves the displayed bound.
The original small-dimension and narrow-width R96 fallbacks retain the
literal dirty thresholds and the sole final clean qubit. A source with
a sufficiently strong explicit precision-uniform ledger can still be
called directly; the fixed-error amplification uses its separate exact
address-sector premise.

### Retained vote alternative and the address condition

The earlier vote construction remains an elementary, less efficient
proof under the same exact address-preservation premise. Convert a
lifted sign to an XOR, record odd k votes in fresh zero bits with the
same dirty core, phase their majority and reverse the actual recording.
At fixed x each wrong-record Kraus operator has norm at most
$`\epsilon_0`$; orthogonality of records bounds wrong-majority
squared mass by $`2^k\epsilon_0^k`$. Its clean-code error is at most
$`2(2\epsilon_0)^{k/2}`$. Applying this to $`u f(x)`$, promoting
the votes at error factor $`8\cdot2^k`$ and taking the Z echo gives
full error $`32(\sqrt{8\epsilon_0})^k`$. This uses at most two
lifted address bits and polynomial overhead in $`1/\tau`$; the QSP
construction supersedes that precision cost in the frame recurrence.

Generic coherent address motion does not satisfy either amplification
premise. For the vote example, tilt the $`|1\rangle`$ projector by
angle theta and control every vote XOR on that same tilted projector.
Each call has error $`2|\sin\theta|`$, but every odd majority returns
the same tilted phase. For QSP, the separation into the prescribed
$`\sigma_x`$ sectors similarly requires their exact preservation;
spectral rounding alone does not correct rotated address eigenspaces.

## 4. A paid matrix-vector product and the missing readout

For $`p=r+s`$ let $`R=2^r`$, $`S=2^s`$. Index $`R`$
coordinates by subsets of the first input half. Let $`L_u=I+J_u`$
be binary unitriangular, with zero empty row for $`J_u`$ and

```math
 (J_u z)_A=u_{\max A}z_{A\setminus\{\max A\}}\quad(A\ne\varnothing).
```

Writing $`a,b`$ for the two empty-subset vectors gives
$`L_u^{-1}a=m(u)`$, the vector of all squarefree monomials including
one. The same construction gives $`L_v^{-1}b=m(v)`$. Compute the
joint Boolean algebraic-normal-form coefficient matrix $`C`$ offline,
so $`f(u,v)=m(u)^TCm(v)`$. Define over the binary field

```math
 P=\begin{pmatrix}L_u^T&C\\0&L_v\end{pmatrix},\quad
 e=(a,0),\quad h=(0,b),\quad K=P+he^T .
```

Here $`e^TP^{-1}h=f(u,v)`$. The kernel equation forces
$`z=P^{-1}h(e^Tz)`$, so $`K`$ has nullity exactly $`f(u,v)`$.
The reversible operation $`w\gets w+Kz`$ needs one CNOT per constant
entry and one address-controlled Toffoli per varying entry, of which
there are $`R+S-2`$. Hadamard-transform the output and use disjoint
private phase masks for the CCZ terms; shared public wires serve as
Clifford controls. With a balanced split the complete matvec costs
$`O(\sqrt{2^p})`$ T gates and dirty work, $`O(2^p)`$ elementary
Cliffords and $`O(1)`$ T-depth without clean inputs.
[Native parallel phase masks](PARALLEL_BILINEAR_WORD_QUERY.md)

An equivalent raw-table graph has a right binary tree directed away
from its root, a left tree directed toward its root, all programmed
right-leaf to left-leaf crossings, and one marked edge from the left
root back to the right root. The right address enables the agreeing
child edge, and the left address enables the agreeing incoming edge
at each left parent. Exactly one right leaf can be reached from its
root and one left leaf can reach its root. The unique possible cycle
exists precisely when the selected table entry is one. Its length
is $`j=p+2`$. For adjacency $`G`$ the root-root entry of
$`G^j`$ equals $`f(u,v)`$. If $`Q`$ implements the shear
$`w\gets w+G^jz`$, set $`V=Q^\dagger Z_{w,\rm root}QZ_{w,\rm root}`$.
Then $`V X_{z,\rm root}V^\dagger X_{z,\rm root}=(-1)^fI`$
on arbitrary work. Four uses of the shear or its adjoint multiply its
complete approximation error by at most four.

The paid shear uses $`K`$ or $`I+G`$. A shallow shear by the power
$`G^j`$, or a complete readout of the selected corank bit, is the
remaining operation. For example, the uniform average of the matvec
phase $`(-1)^{y^TKz}`$ is $`2^{-\mathrm{rank}K}`$, an
exponentially small distinction in the carrier dimension. A shallow
matvec or an invertible crossing alone cannot be substituted for the
scalar sign in §1. [Related observable-pencil source](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md)

<a id="self-inverse-shear-tables"></a>

## 5. Self-inverse shear tables suffice

Let an arbitrary $`p=2r`$ bit function have the raw truth matrix
$`C\in\mathbb F_2^{m\times m}`$, $`m=2^r`$, with
$`C_{uv}=f(u,v)`$. No invertibility hypothesis on $`C`$ is needed. Set

```math
 A_C=\begin{pmatrix}I_m&C\\0&I_m\end{pmatrix},\qquad
 A_C^2=I_{2m},\qquad (A_C-I_{2m})^2=0.
```

As an operation on two physical banks, $`(Y,Z)\mapsto(Y+CZ,Z)`$ is
a fixed Clifford word with at most $`m^2`$ CNOTs and is its own inverse.
Its scalar **entry** table instead takes binary indices
$`(\alpha,u)`$ and $`(\beta,v)`$, where the two block bits are arbitrary:

```math
 h_C(\alpha,u;\beta,v)
 =\delta_{\alpha,\beta}\delta_{u,v}+(1+\alpha)\beta C_{uv}
 \quad\text{in }\mathbb F_2.
```

For every incoming $`\alpha,\beta,u,v`$,

```math
 \sum_{s,t\in\mathbb F_2}
 h_C(\alpha+s,u;\beta+t,v)=C_{uv}.
```

The diagonal contribution occurs twice and the upper-right contribution
once. Thus a complete scalar entry source for just these shear matrices
suffices for arbitrary prescribed scalar tables.

### Actual inverses, literal phase and leakage

Suppose a physical unitary $`V_C`$ obeys, for its initialized embedding $`J`$,

```math
 \|V_CJ-e^{i\theta_C}J(D_{h_C}\otimes I_{\rm work})\|\le\epsilon.
```

The phase $`\theta_C`$ may depend on the compiled table but is common to
all address and work inputs. There is no promise on the source outside
the initialized code. Unitarity gives the same bound for the actual
adjoint, with phase $`e^{-i\theta_C}`$: multiply the source error on the
left by $`V_C^\dagger`$ and on the right by
$`D_{h_C}^\dagger\otimes I_{\rm work}`$.

Apply the following operations in chronological order:

```math
 V_C,\ X_\alpha,\ V_C^\dagger,\ X_\beta,\
 V_C,\ X_\alpha,\ V_C^\dagger,\ X_\beta.
```

On the ideal trajectory the four calls visit all four block-bit settings.
Both block bits return, the fixed source phase cancels, and the mixed
difference above gives $`D_f\otimes I`$ on every borrowed-work input.
Replacing the four calls one at a time along that trajectory gives
initialized-isometry error at most $`4\epsilon`$. Each ideal prefix
preserves the initialized code; the remaining actual suffix is unitary.
This hybrid therefore includes intermediate address motion, clean leakage
and dirty work entangled with arbitrary references. It does not reset
approximate work or replace an actual inverse with a second forward call.

### All-width resource ledger

Suppose the source and its actual reverse have common bounds
$`T_s,G_s,d_s`$ and private dirty allocation $`B_s`$. The extraction uses
at most $`4T_s`$ T gates, $`4d_s`$ T-depth and $`4G_s+O(1)`$
elementary Cliffords, with the same initialized bank and private dirty pool
reused sequentially. The two block bits are borrowed, not initialized.
The four X gates and their literal native reverses have constant Clifford
cost and zero T-depth.

For odd $`p`$, first append one unused arbitrary bit to $`f`$ and split
the resulting even address into equal halves. This duplicates truth-table
rows or columns; the desired sign is independent of the padding bit and
returns it. The source arity is consequently $`a=p+2`$ for even $`p`$
and $`a=p+3`$ for odd $`p`$. Its table volume is at most $`8\cdot2^p`$.
Reserve those two or three bits from the requested total dirty allocation
$`B`$, leaving $`B_s=B-(a-p)`$. Increasing the absolute lower-width
constant ensures $`B_s\ge C(a+2)`$ and $`B_s=\Theta(B)`$ whenever
$`B\ge C'(p+2)`$.

It follows that the source budgets of §1, assumed only for scalar entry
tables of the matrices $`A_C`$, give the same budgets for every $`f`$,
with additive depth $`O(d(p+3))`$. Exact sources transfer exactly; at
requested error $`\tau`$, ask for source error $`\tau/4`$. A bounded
initialized supply remains bounded. If a fixed-accuracy source exactly
preserves all computational address sectors on all physical sectors, the
completed wrapper does too, so the amplification hypothesis of §3 is
retained. Generic initialized-isometry approximation alone still does not
establish that stronger hypothesis.

**The unpaid operation.** Applying the Clifford shear to $`2m`$ physical
bank bits does not implement its scalar entry table on $`2\log_2(2m)`$
binary address bits. For fixed indices, conjugating $`X_{Z,v}`$ by the
shear and taking its Pauli commutator with $`Z_{Y,u}`$ gives
$`(-1)^{C_{uv}}I`$. Coherently selecting those two physical Paulis is
the missing operation in that construction. A separately synthesized
dirty indicator retains the old logarithmic address depth. A useful next
construction must pay the joint indexed commutator or another complete
entry readout. Invertibility and an inexpensive inverse are already
available; this reduction does not improve the unconditional frame bound.

<a id="thin-kernel-rank-obstruction"></a>

## 6. A rank obstruction for thin-kernel filters

The marked pencil in §4 gives a particularly small invariant subspace.
This permits a direct screen for attempts to turn its invariant projector
into a scalar sign by a short word of reflections.

Let $`V=B_qR_qB_{q-1}\cdots R_1B_0`$ be a unitary word, let
$`V_0=B_qB_{q-1}\cdots B_0`$ be the word with those insertions removed,
and assume $`\mathrm{rank}(R_j-I)\le s_j`$. For an isometry
$`J:\mathbb C^R\to\mathcal H`$ with $`V_0J=J`$,

```math
 \sum_j s_j\lt R\quad\Longrightarrow\quad\|VJ+J\|_{\rm op}=2.
```

Indeed, telescoping the two words gives
$`\mathrm{rank}(V-V_0)\le\sum_j s_j`$. Hence
$`(V-V_0)J`$ has a nonzero kernel vector. On its normalization the full
word still acts as identity, at distance two from the required literal minus sign.
The upper bound two follows because both maps are isometries. Arbitrary
unitary interleavings and conjugations are allowed in this lemma.

For a fixed address, write $`m`$ for the binary dimension of $`K`$ in
§4. Averaging translations by $`\mathrm{im}K`$ on an $`m`$-qubit
carrier gives the invariant projector $`\Pi_K`$. Its rank equals the
number of translation cosets, $`2^{m-\mathrm{rank}K}=2^{f(x)}`$.
A unitary insertion $`I+(e^{i\phi}-1)\Pi_K`$ therefore differs from
identity in rank at most two. With $`h`$ spectator qubits the bound is
$`2^{h+1}`$. If $`k`$ of the $`m+h`$ work qubits are initialized, the
complete initialized code has dimension $`2^{m+h-k}`$.

Consequently, at an address with $`f(x)=1`$, fewer than
$`2^{m-k-1}`$ such insertions cannot approximate the scalar minus sign
with error below two, provided the erased word fixes that code. Arbitrary
dirty padding cancels from the bound. To use these fixed-address ranks,
the ideal insertions and interleavings must preserve computational address
sectors. The full-space rank lemma does not itself require that restriction.

If an emitted word $`\widetilde V`$ satisfies
$`\|(\widetilde V-V)J\|\le\Delta`$, the same argument gives
$`\|\widetilde VJ+J\|\ge2-\Delta`$. Full-operator primitive error
bounds imply such a $`\Delta`$ by telescoping. Initialized-only primitive
bounds need their own valid hybrid; they cannot be summed through arbitrary
intermediate leaked work without one.

This excludes the stated thin-filter extraction with bounded clean work
and polynomially many insertions in $`m`$ at target error $`\tau`$ when
$`\tau+\Delta\lt2`$. It does not
exclude a representation with a larger signal subspace, a baseline already
carrying the answer, or unrestricted scalar synthesis. It is an architecture
screen, not a native T-depth lower bound. No physical projection or
postselection is assumed in the argument.

## Finite evidence

The [finite checks](../../tests/test_prescribed_scalar_tables.py) recover
every two-bit target from shifted even-parity tables, verify the
derivative vote as a coherent unitary, and test selected Walsh characters,
coranks and marked-edge powers. They also check the shear completion and
four-call extraction on every two-bit table, selected larger tables and an
odd-arity table with arbitrary padding. A perturbed source fixture includes
address motion, nonzero clean leakage, arbitrary dirty inputs and a fixed
phase, and checks the actual-adjoint error bound. It is not a native emitter.
The asymptotic packing, density, rank and frame bounds above are analytic.
