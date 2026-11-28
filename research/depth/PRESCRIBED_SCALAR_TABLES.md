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

<a id="3-precision-and-completed-address-return"></a>
<a id="precision-and-completed-address-return"></a>

## 3. Precision and full-space sign return

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

It retains one final clean qubit. The full-space signed-gap condition
below is a weaker sufficient premise for the same amplification: it
allows nonflat spectra and completed address motion. A qualifying native
source is still a premise; the unconditional frontier is unchanged.

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

<a id="full-space-signed-gap-source"></a>

### A full-space signed gap suffices

The amplification does not require V to approximate $`D_{tf}`$, or
to preserve every address. Let V be an actual zero-clean native word
on an a-bit address, an arbitrary signal t and arbitrary work. Define
S, R, W and $`\mathcal H`$ as above. For
$`D=D_f\otimes I_{t,\mathrm{work}}`$, it suffices that

```math
 \mathrm{sgn}(\mathcal H)=D,\qquad
 |\mathcal H|\succeq\gamma I,\qquad 0\lt\gamma\le1.
```

These are full-physical-space statements. Equivalently, D commutes
with $`\mathcal H`$ and $`D\mathcal H\succeq\gamma I`$.
The positive and negative eigenspaces must be the prescribed Boolean
sign sectors, but magnitudes within them can vary with address and work.
Neither V nor its actual inverse has to preserve those sectors. A gap
only on an initialized code does not establish this interface; the
earlier bounded-clean XOR promotion does not automatically promote a
raw signed-gap promise. The [one-clean decision route](#one-clean-decision-source)
below instead amplifies on the initialized code before promoting the
resulting complete query.

**Whole-unitary proof.** The reflection identities give
$`[\mathcal H,S]=[\mathcal H,R]=0`$. Consequently D, as a
spectral function of $`\mathcal H`$, commutes with S, R and every
QSP word $`U_\Phi`$. Jordan's two-reflection decomposition and the
one-dimensional endpoint argument above apply on the entire physical
space, without selecting computational address sectors. Thus

```math
 \tfrac12(U_\Phi+U_\Phi^\dagger)=\widehat p(\mathcal H),
 \qquad
 (U_\Phi-D)^\dagger(U_\Phi-D)
   =2\bigl(I-D\widehat p(\mathcal H)\bigr).
```

For $`0\lt\tau\le1/2`$, the bounded odd sign polynomial of
[F55](../../docs/reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries),
v1 Lemma 25, has degree
$`q=O(\gamma^{-1}\log(1/\tau))`$, is bounded by one on
$`[-1,1]`$, and approximates the sign function to
$`\tau^2/16`$ on $`[-1,-\gamma]\cup[\gamma,1]`$.
Certify its QSP factorization to uniform error $`\tau^2/16`$.
The displayed identity then bounds the whole-unitary error by
$`\tau/2`$. Endpoint values need only have that accuracy; no
postselection or scalar real-part extraction is performed.
Compile the physical rotations $`e^{i\phi X_t}`$ with total
literal error at most $`\tau/2`$, as above.

The same forward/inverse ledger holds with k replaced by q:
$`2q`$ actual V/V-adjoint appearances, local synthesis cost
$`O(q\log((q+1)/\tau))`$ in each gate/depth bound, one reused
arbitrary signal and no clean work. The explicit majority polynomial
also works, with $`O(\gamma^{-2}\log(1/\tau))`$ degree; Lemma 25
gives the sharper gap dependence. For fixed positive gamma, both give
the logarithmic accuracy cost used in the frame recurrence.

This is a more flexible algebraic source target, not a construction of
it. Once supplied, a constant gap yields a complete fixed-accuracy scalar
with constant source-call overhead. If gamma decreases with arity, the
factor $`\gamma^{-1}`$ must be charged even at fixed final accuracy.
In particular, multiplying raw $`O(2^{a/2})`$ T and $`O(2^a)`$
Clifford budgets by an unbounded factor does not preserve the desired
full-volume counts without additional source savings. The
[two-source-band screen](FOUR_LAYER_SCALAR_GAP.md#two-band-signed-gap-sources)
permits arbitrary width and address motion, and forces an exponentially
small gap for a single marked table and its Pauli-balanced lift.

<a id="one-clean-decision-source"></a>

### Amplify a one-clean decision before imposing complete return

An actual source can instead use one initialized decision signal c.
Let J insert $`c=0`$, let V act on c, an a-bit address and arbitrary
work, and put

```math
 S=Z_c,\qquad R=V^\dagger S V,\qquad
 A=J^\dagger VJ,\qquad M=J^\dagger RJ=2A^\dagger A-I.
```

For $`D=D_f\otimes I_{\mathrm{work}}`$, require

```math
 \mathrm{sgn}(M)=D,\qquad |M|\succeq\gamma I,
 \qquad 0\lt\gamma\le1.
```

Equivalently, $`[M,D]=0`$ and $`DM\succeq\gamma I`$. The accepted
Gram $`A^\dagger A`$ can have any spectrum bounded below by
$`(1+\gamma)/2`$ on the f=0 subspace and bounded above by
$`(1-\gamma)/2`$ on the f=1 subspace. Its accepted polar factor is
unrestricted. V may move the address, leave arbitrary quantum garbage
and fail to return c. Its rejected-sector sign need not have the
prescribed address labeling. All bounds concern the actual V and its
actual inverse, with no source-return premise.

This is a coherent decision condition on every input in each target
sign subspace, including arbitrary dirty states and references. It is
stronger than a success probability checked only on computational
dirty inputs, or a small signed bias averaged over dirty inputs.
The exact sign subspaces also matter: a spectral gap with a rotated
positive/negative decomposition does not meet this condition.

**Two checks on candidate decision sources.** Suppose a candidate
preserves every address and has fibers
$`V_x=G_xV_0G_x^\dagger`$, where $`G_x`$ preserves the clean-zero
code. Then $`G_x`$ commutes with S, its restriction $`g_x=J^\dagger G_xJ`$
is unitary, and

```math
 M_x=g_xM_0g_x^\dagger.
```

Every address therefore has the same compressed spectrum. Such a
family cannot be positive definite at one address and negative definite
at another, so it cannot decide a nonconstant table with the required
gap. This check does not apply when the gauge moves the clean code or
when the source mixes addresses.

A separate rank check applies to a general accepted projector Pi.
For any specified accepting input subspace with isometry $`J_E`$ and
dimension $`R_E`$,

```math
 \mathrm{rank}(J_E^\dagger V^\dagger\Pi VJ_E)
 \le\mathrm{rank}(\Pi)=r.
```

If $`r\lt R_E`$, its compression has a kernel and
$`2J_E^\dagger V^\dagger\Pi VJ_E-I`$ has a minus-one eigenvalue.
It cannot have a positive gap on that entire accepting subspace.
For a read-only source, this can be applied to all dirty inputs at one
accepting address. With address motion, use the actual global accepted
projector and accepting input subspace; a rank counted only in one
address fiber is not a bound on the coupled operation. A kernel on a
prescribed rejecting subspace is allowed by the decision premise.

**Recovery on the same clean signal.** Set
$`W=RS`$ and $`\mathcal H=(RS+SR)/2`$. Reflection identities give

```math
 [\mathcal H,S]=[\mathcal H,R]=0,\qquad
 \mathcal HJ=JM.
```

Use the same real-polynomial QSP word, now with physical phases
$`e^{i\phi Z_c}`$ on the existing clean signal. For its exact
phase-list unitary U and realized polynomial $`\widehat p`$,

```math
 \tfrac12(U+U^\dagger)=\widehat p(\mathcal H),\qquad
 [U,\mathcal H]=0.
```

Writing $`B=J^\dagger UJ`$, the intertwining identity implies
$`[B,M]=0`$, hence $`[B,D]=0`$. Therefore the complete initialized
error, including leakage out of c=0, satisfies the exact identity

```math
 (UJ-JD)^\dagger(UJ-JD)
 =2I-B^\dagger D-DB
 =2\bigl(I-D\widehat p(M)\bigr).
```

For requested $`0\lt\epsilon\le1/2`$, use the bounded odd sign
polynomial from F55 at degree
$`q=O(\gamma^{-1}\log(1/\epsilon))`$, with approximation error
$`\epsilon^2/16`$ outside the gap. Certify the phase factorization
to the same polynomial error. The last identity gives initialized
error at most $`\epsilon/2`$. Compile the q+1 physical single-qubit
phases with total literal operator error at most $`\epsilon/2`$,
using F5 as above. The emitted word then obeys
$`\|UJ-JD\|\le\epsilon`$. The same bound holds for its actual
inverse by the initialized-isometry inverse identity; it does not
require exact clean return by the approximant.

There are exactly 2q appearances of V or its actual adjoint. For
common source bounds $`T_V,G_V,d_V,w_V`$, where w excludes the
address and c, the amplified scalar has common forward/inverse bounds

```math
 \begin{aligned}
 T_U&\le2qT_V+O\!\left(q\log\frac{q+1}{\epsilon}\right),\\
 G_U&\le2qG_V+O\!\left(q\log\frac{q+1}{\epsilon}\right),\\
 d_U&\le2qd_V+O\!\left(q\log\frac{q+1}{\epsilon}\right),\\
 w_U&\le w_V,\qquad a_{\mathrm{clean}}=1.
 \end{aligned}
```

No fresh phase-estimation clock, controlled V, reset or postselection
is used. The original decision bias is not a per-call error summed
over the filter: it supplies the spectral premise for that one word.

**Promotion after amplification.** To produce a zero-clean scalar for
f, separately compile the explicit lifted table $`g(x,t)=t f(x)`$.
This is an additional table invocation, not a free controlled source.
Amplify its one-clean decision to initialized error
$`\epsilon=\tau/32`$. Conjugation by $`H_t`$ makes a one-clean XOR
of f into arbitrary t. The [serial one-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
uses eight appearances of that amplified query or its actual inverse,
has full error at most $`16\epsilon`$, and permits arbitrary address
motion and off-code behavior. Its final Z echo uses two promoted
query appearances and gives the scalar at full operator error at most
$`32\epsilon=\tau`$.

Thus the zero-clean scalar uses exactly 32q appearances of the raw
lifted-table V or its actual adjoint, with common bounds

```math
 \begin{aligned}
 T_\tau&\le32qT_V+O\!\left(q\log\frac{q+1}{\tau}\right)+O(1),\\
 G_\tau&\le32qG_V+O\!\left(q\log\frac{q+1}{\tau}\right)+O(1),\\
 D_\tau&\le32qd_V+O\!\left(q\log\frac{q+1}{\tau}\right)+O(1),\\
 w_\tau&\le w_V+4,\qquad a_{\mathrm{clean}}=0.
 \end{aligned}
```

The four extra arbitrary wires are the lifted old target, former clean
signal, promotion toggle and final XOR target. Each pool is reused
sequentially. Fixed positive gamma therefore supplies the same
logarithmic precision cost as the full-space source, with one additional
source-address bit and constant allocation factors.

**A native example separating the premises.** On address bit x and
clean signal c, define actual Clifford words and one T gate by

```math
 C=\mathrm{CX}_{x\to c}\mathrm{CX}_{c\to x},\qquad
 K=\mathrm{CX}_{x\to c}H_x,\qquad V=CKT_xK^\dagger.
```

Rightmost factors act first. The native T conjugate is, up to its
common scalar, $`e^{-i\pi X_xX_c/8}`$; that scalar cancels in R.
Since $`C^\dagger Z_cC=Z_x`$, direct Pauli conjugation gives

```math
 R=\frac{Z_x+Y_xX_c}{\sqrt2},\qquad
 M=\frac{Z_x}{\sqrt2},\qquad
 \mathcal H=\frac{Z_xZ_c}{\sqrt2}.
```

This one-T-band decision has gap $`1/\sqrt2`$ but its full-space
sign is $`Z_xZ_c`$, not $`Z_x\otimes I_c`$. V is nonmonomial and
does not return c. A Clifford address change extends the example to
every nonconstant affine parity at every address size. It separates
the source hypotheses; it supplies no arbitrary-table depth bound.

The [projector-Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
already permits arbitrary accepted partial isometries, but uses a
prescribed approximately flat Gram $`\Pi/2`$ and a fixed polynomial.
The present decision premise uses nonflat spectra separated across one
half. The [sampled-vote source](RETURNED_MASK_AND_ZERO_TEST.md#4-sampling-with-arbitrary-dirty-origins-and-normalizing-size)
already obtains a compressed sign gap for its specific sampler. Here
QSP on the existing clean signal supplies logarithmic precision
recovery without a fresh clock. None of these interfaces constructs
the required general native source.

**The two-band screen still applies.** In the c=0 and c=1 blocks write

```math
R=\begin{pmatrix}M&B_0\\B_0^\dagger&N_0\end{pmatrix}.
```

Then $`\mathcal H=\mathrm{diag}(M,-N_0)`$. From R squared equal to
identity, $`MB_0+B_0N_0=0`$ and
$`B_0B_0^\dagger=I-M^2`$. These identities pair every eigenvalue of
M strictly between minus one and one with the same eigenvalue of
$`-N_0`$, with its multiplicity. The two blocks have equal size and
R and S are both balanced involutions, so their remaining endpoint
multiplicities also agree. Hence the two block spectra coincide.
For a table with smaller nonempty sign fraction epsilon, the full
$`\mathcal H`$ has that same fraction. The
[two-band spectral bound](FOUR_LAYER_SCALAR_GAP.md#two-band-signed-gap-sources)
therefore still excludes a uniform constant-gap decision source with
at most two raw T bands for all single-mark tables. The clean decision
interface does not circumvent that result or give a growing-depth
lower bound.

<a id="recorded-address-decision-source"></a>

### Record a short address label before reading a moving source

Let $`M=J^\dagger V^\dagger Z_cVJ`$ for a raw one-clean source V.
Write $`M_{xy}`$ for its computational address blocks on the entire
arbitrary work space. Suppose the diagonal blocks have the prescribed
sign and gap, in the sense that $`(-1)^{f(x)}M_{xx}\succeq\gamma I`$.
Let h be a k-bit label for which $`h(x)\ne h(y)`$ whenever $`x\ne y`$
and $`M_{xy}\ne0`$. To use this promise the compiler must supply h
and pay its full returned XOR into the tag register. For example, if
all nonzero differences $`x\oplus y`$ in the off-diagonal support lie
in a supplied set $`\mathcal S`$, a linear h exists for any k with
$`2^k\gt|\mathcal S|`$. A uniformly chosen binary k-by-a matrix kills
each prescribed nonzero difference with probability $`2^{-k}`$.
Finite search gives a valid matrix, and at most ak CNOTs compute its
label with no T layer or private work.

Append k clean tags, record h coherently, and apply V on its original
registers. Measuring the same output $`Z_c`$ as in the original source
now gives the compressed observable

```math
M'_{xy}=[h(x)=h(y)]M_{xy}.
```

This is a matrix-entry identity on the entire initialized code. In
particular the dirty input need not factor across address blocks or
from a reference. It makes M' block diagonal with sign $`D_f`$ and gap
$`\gamma`$. V can move the address and the record can remain in the
raw output. The actual adjoint reverses both V and the record. If
instead the principal blocks of M on the two prescribed sign sectors
already have the uniform signed gap, h only needs to separate the
off-diagonal blocks between *opposite* signs. Within either sign sector,
pinching to the h blocks preserves the positive lower bound.

For fixed k and fixed positive gamma this is another sufficient source
for the [square-root-width transfer](#square-root-width-signed-gap-transfer).
In the same-signal recovery replace S by $`2JJ^\dagger-I`$, with J now
inserting all k+1 initialized bits. Keep R as the conjugate of the
single measured $`Z_c`$. The identities $`\mathcal HJ=JM'`$ and
$`[R,\mathcal H]=[S,\mathcal H]=0`$ are unchanged. The QSP phases
$`e^{i\phi S}`$ use the $`2^{k+1}-1`$ nonidentity Z products in the
expansion of S and exact Clifford parity exposure. Omitting its identity
term gives only one common scalar for the entire source, independent of
the input. This scalar cancels in the subsequent
[clean-code XOR promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#1-the-source-is-a-complete-query-not-a-prepared-state).
For the lifted table, amplify first and then promote k+1 clean bits
with their actual reversals. With k fixed the promotion costs only a
constant factor and the final scalar has literal phase and complete
return. The degree and per-phase accuracy remain as in the one-clean
recovery. A linear tag adds $`O(ak)`$ Cliffords per raw source call, and
the tags become a constant number of borrowed wires after promotion.
The source's T-count, Clifford count, width and sublogarithmic-depth
criterion are therefore preserved in the frame transfer. Any growing
tag budget needs a new resource and accuracy ledger.

<a id="summable-address-motion"></a>

**Small total coupling permits many differences.** Let $`M^{(u)}`$
keep precisely the blocks $`M_{x,x\oplus u}`$ of M. It has norm
$`\max_x\|M_{x,x\oplus u}\|`$, since it is a block diagonal operator
times an address translation. Suppose the diagonal signed gap is
$`\gamma\gt0`$ and the source has the uniform bound

```math
\sum_{u\ne0}\|M^{(u)}\|\le B
```

for a fixed B, also for every lifted table used in promotion. For a
linear k-bit tag, the expected sum of the norms of surviving nonzero
differences is at most $`B2^{-k}`$. Use a certified upper bound on B
and a positive lower bound on gamma with strict slack in the search.
The finitely many block norms of a native word admit certified
evaluation. For fixed $`0\lt\tau\le1/2`$ take
$`k=O_{B,\gamma}(1+\log(1/\tau))`$ so the tagged compression M' is
within $`\delta=\gamma\tau/32`$ of the diagonal compression
$`M^{(0)}`$. Both spectra avoid the interval of radius
$`\gamma-\delta`$ about zero. The integral representation of the
sign by the sum of resolvents at $`\pm it`$, for $`t\ge0`$, gives

```math
\|\mathrm{sgn}(M')-D_f\otimes I\|
\le \frac{2\delta}{\pi}\int_0^\infty
       \frac{dt}{(\gamma-\delta)^2+t^2}
=\frac{\delta}{\gamma-\delta}.
```

Here the difference of each resolvent is the product of the two
resolvents with the perturbation inserted. The bound is in full
operator norm, including noncommuting compression errors.

To promote at tolerance tau, amplify the lifted decision to the
*actual* sign of M' with initialized-isometry error
$`\epsilon=\tau/(64\cdot2^{k+1})`$, up to its common scalar.
After the target Hadamards, compare with the exact clean-sector
extension of that unitary, as in
[promotion §4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#4-approximate-clean-return-still-gives-a-full-operator-query).
Only the discrepancy from this extension costs $`8\cdot2^{k+1}\epsilon`$.
The extension's clean-sector unitary is within
$`\delta/(\gamma-\delta)`$ of the ideal XOR. The disjoint-sector
estimate there charges this second error only twice. After the final
Z echo the sum of these two terms is still below tau. Applying
the leakage multiplier directly to the address mixing would lose
the gain from the tag.

At fixed B and gamma, all native overheads of this recovery are
polynomial in $`1/\tau`$. The QSP degree is
$`O_\gamma(1+\log(1/\tau)+k)`$. Each phase has at most
$`2^{k+1}-1`$ parity terms synthesized to a precision with logarithm
$`O_{B,\gamma}(1+\log(1/\tau))`$. Linear recording costs ak CNOTs.
Serial promotion visits $`2^{k+1}`$ sectors and reuses its returned
zero-test pool. Thus there is a fixed exponent C for which a raw
source with the square-root T-count and width and linear-volume
Clifford count yields a complete zero-clean scalar with

```math
\begin{aligned}
T_\tau&=O_{B,\gamma}(\tau^{-C}2^{a/2}),&
G_\tau&=O_{B,\gamma}(\tau^{-C}2^a),\\
D_\tau&=O_{B,\gamma}(\tau^{-C}(d(a)+1)),&
w_\tau&=w_V+O_{B,\gamma}(\tau^{-C}).
\end{aligned}
```

These bounds include the actual reverse. The lifted arity adds one
to a and the width uses that source's work bound. In the
[square-root-width recurrence](#square-root-width-signed-gap-transfer),
early-stage inverse tolerances are polynomial in $`h(m)+2`$.
Increasing the arity margin absorbs the resulting polynomial factors
in counts and width. The full-volume sources still have fixed
L-dependent tolerances. The new frame depth at fixed L is

```math
O_{L,B,\gamma}\!\left(\log\log(n+4)+H_n^{O(1)}(d(n+C_L)+1)\right)
```

with the same simultaneous frame counts and width, using the usual
nondecreasing depth envelope d. In particular
$`d(a)=O(\log a/\log\log a)`$ still suffices. This conditional
route requires a paid bound on the total off-diagonal coupling.
Correct diagonal blocks by themselves supply no such bound.

**A uniform recording rule needs long tags for unrestricted motion.**
Fix a sign table and let S be a set of m addresses of a common sign
$`s\in\{+1,-1\}`$, with m at least half the address count. Consider
any read-only recording isometry that maps an arbitrary R-dimensional
helper, at each address x, into an environment of dimension $`2^kR`$.
This includes arbitrary quantum records, not just computational tags.
The recording acts independently of the unknown source observable and
does not act on its original work B. Allow B to have dimension at least
mR, padding to a power of two if needed.

Take the maximally entangled input between the mR-dimensional address
and helper subspace supported on S and an mR-dimensional subspace of B.
After recording and tracing out its environment, the density matrix
$`\rho`$ on address and B has rank at most $`2^kR`$. Its address
diagonal blocks obey $`\rho_{xx}\preceq I_B/(mR)`$ by isometry in
each address fiber. Let P be its spectral projector on eigenvalues
at least $`t/(mR)`$ for a fixed $`t\gt2`$. The projector is supported
on S and satisfies

```math
P_{xx}\preceq I_B/t,\qquad
\mathrm{Tr}(P\rho)\ge 1-t2^k/m .
```

Indeed $`P\preceq(mR/t)\rho`$, and at most $`2^kR`$ nonzero
eigenvalues are omitted. Set $`M=D_f\otimes I_B-2sP`$. This is
a Hermitian involution since P lies in the s sector. Its diagonal
blocks satisfy the prescribed gap $`1-2/t`$. It is a legitimate
abstract one-clean decision observable: on its negative eigenspace
flip the clean output and on its positive eigenspace leave it fixed.
On the entangled input above the recorded observable has signed
expectation at most $`-1+2t2^k/m`$. Thus a uniform read-only recorder
requires $`k\ge\log_2 m-O_t(1)`$ merely to guarantee a positive signed
gap for every such M. The argument allows arbitrary dirty padding in
the recorder. It makes no native cost claim for this adversarial M
and does not bound quantum records chosen from a particular source's
description.

<a id="adaptive-classical-tag-limit"></a>

**A shallow decision defeats every small classical recording.**

The address tag in the preceding construction may depend on the emitted
source. Even allowing this choice, uniformly signed address diagonals
alone do not yield a bounded-size deterministic classical recording.
The following source has one clean signal, $`O_\epsilon(a)`$ T gates,
elementary Cliffords, and dirty wires, and $`O_\epsilon(1)`$ native
T-depth for any fixed $`0\lt\epsilon\lt1/4`$. The same bounds include
its actual inverse. Every address diagonal has signed gap at least
$`1/2-\epsilon`$. For every address function h with fewer than
$`2^{a-2}/3`$ tags, even if chosen knowing the source, the recorded
compression has a signed expectation at most $`-1+\epsilon`$.

Let $`a\ge4`$, let $`f(x)=x_1`$, and reserve four fixed anchor
addresses $`B_0,\ldots,B_3`$ with prefix $`01`$. Set the last two
bits of $`B_0`$ to zero and obtain the other anchors by varying just
those bits. Four arbitrary dirty lists of $`a-2`$ bits encode
addresses $`A_0,\ldots,A_3`$ with prefix $`00`$. Write z for all four
list words. On address X, let $`U_z`$ first apply Hadamards to the last
two bits and then perform the four transpositions $`(A_j,B_j)`$ in a
fixed order. Define the projector

```math
P=\sum_z |z\rangle\langle z|\otimes U_z|B_0\rangle\langle B_0|U_z^\dagger,
```

tensored with identity on private helpers. For every z the four
images of the anchors under the transposition product are distinct
because that product is a permutation. The normalized state defining
P in that sector is uniform on the four images, all of which have
first bit zero. If the four $`A_j`$ are distinct, the four images are
exactly those $`A_j`$. This description gives a full operator on the
entire arbitrary list register, including repetitions and coherences.

The ideal source toggles c by $`x_1`$ and on the range of P.
These actions commute because P lies in the $`x_1=0`$ sector. Its
clean-code observable is

```math
M_0=D_f\otimes I-2P,
\qquad
(-1)^{f(x)}(M_0)_{xx}\succeq\tfrac12 I.
```

Indeed $`P_{xx}\preceq I/4`$, with no condition on the dirty list or
the private work. To implement the second toggle, conjugate the
complete XOR at the single marked address $`B_0`$ by $`U_z`$ and its
actual inverse.

**Native transposition.** For any one of the address pairs A,B, let
K toggle an arbitrary shared helper w by
$`[x=A]\oplus[x=B]`$. It is two complete zero-test XORs with the
address shifted by A or B and shifted back by CNOTs and X gates.
The [zero-reflection proof](EXPANDER_ZERO_REFLECTION.md) supplies a
fixed-error complete zero-test XOR with $`O_\epsilon(a)`$ T-count,
Clifford count and dirty work at $`O_\epsilon(1)`$ native T-depth.
Its padding-sector product for the phase also works for the XOR,
since at the zero unpadded address exactly one padding shift triggers
the ideal XOR. The same full-operator hybrid pays the at most four
padded calls. In particular the zero-test word may move the address
internally.

Let L add $`w(A\oplus B)`$ into X, with all other wires returned.
This is one parallel batch of $`O(a)`$ Toffolis and CNOTs. Each
Toffoli uses its own private phase helper, so w and the list input
bits enter only as shared Clifford controls. The predicate in K is
invariant under translation by $`A\oplus B`$. The chronological
echo $`K,L,K^\dagger,L^\dagger`$ therefore adds
$`(A\oplus B)([x=A]\oplus[x=B])`$ to X and returns w in the exact
comparison. This is the desired transposition on every input, with
the list values preserved. When A and B coincide the same identity
gives the identity permutation. All approximate inverse calls in
the emitted echo use the actual reversals.

There are four transpositions in $`U_z`$ and four in its actual
inverse. Each has at most four appearances of the complete zero-test
word or its inverse. The marked-address XOR for $`B_0`$ adds one
call. With each complete zero-test XOR accurate to $`\epsilon/66`$,
the 33-call hybrid makes the source within $`\epsilon/2`$ of its
ideal unitary. Thus its actual compressed observable M obeys
$`\|M-M_0\|\le\epsilon`$. This proves the diagonal gap and all the
stated resource bounds. All zero-test and translation pools can be
reused at ideal returned boundaries by a full-operator hybrid.

**Failure of adaptive classical tagging.** If h has fewer than
$`2^{a-2}/3`$ values, four distinct addresses with prefix $`00`$ share
one value of h. Set the dirty list to these four labels and take X
uniform on them. This vector is in the positive f sector and is a
minus-one eigenvector of $`M_0`$. Copying h into initialized tags has
the same tag on every component of the vector, so its recorded
expectation remains minus one. The recorded observable of the actual
source differs by at most $`\epsilon`$, by compression of the same
unitary conjugation. The input is allowed even though other dirty
lists may contain duplicates. For this source such a classical tag
repair therefore needs at least $`a-3`$ tag bits. The counterexample
concerns the stated deterministic classical recording with tags left
untouched by the source. Other source repairs require a separate argument. Its table f itself
has a Clifford decision.

**Quantum tags from a single-band Pauli orbit.** The same example
also screens one family of source-dependent quantum records. Suppose
the range projector of the record isometry at x is
$`Q_x=D_xQ_0D_x^\dagger`$, where each $`D_x`$ is a Pauli and $`Q_0`$
is the image of k initialized bits under a common one-band native word.
The word and Pauli labels may be chosen knowing the source. Undoing the
last Clifford leaves $`Q_0=EPE^\dagger`$, where E is a diagonal T band
and P a stabilizer projector of codimension k in qubit count. Let
L be the span of the stabilizer X labels, of dimension $`r\le k`$.
Every Q keeps the computational cosets of L invariant and has rank
zero or one on each. For P the diagonal stabilizer sum gives constant
diagonal $`2^{-r}`$ on an active coset and zero on an inactive one.
E preserves each coset and the Pauli translates permute them.

Among the $`2^{a-2}`$ prefix-00 addresses, some coset is active on
at least $`2^{a-2+r-k}`$ by averaging their ranks. The unit sphere
on that coset has a one-eighth net of cardinality at most
$`17^{2^{r+1}}`$. For $`k=o(\log a)`$, four of its active address
lines have normalized vectors within one eighth of the same net
point. Use those four addresses as the dirty list in the native
example. In the recorded code take equal amplitudes on the four
addresses with their corresponding helper vectors. The bad-projector
expectation is the squared norm of the mean of these four vectors,
at least $`(7/8)^2`$. The signed observable therefore has expectation
at most $`1-2(7/8)^2+\epsilon\lt0`$. A common orbit encoder with sublogarithmic
clean width needs to leave this one-band representation to repair
the native example. Two noncommuting encoder bands can change the
coset structure, and their full-space behavior remains open.

<a id="expander-quantum-address-tags"></a>

**Quantum records suppress couplings within each work sector.**
Partition the principal input into p computational label wires and
an arbitrary work register. The label x can include both the address
and any designated dirty input bits. Assume each label-diagonal
block has the prescribed sign and gap $`\gamma\gt0`$, with the sign
depending only on the address. Allow a common orthogonal decomposition
of the remaining work into sectors z preserved by all $`M_{xy}`$. A
single sector is allowed. Suppose that

```math
\max_{x,z}\sum_{y\ne x}\|M_{xy}|_z\|\le B
```

for a uniform constant B, also for all lifted tables. This is a proof
premise that requires no circuit to find z. It allows a dense union of
coupling graphs over dirty inputs.

The [read-only field expander](GATEWISE_READ_ONLY_EQUALITY.md) has
selector branches $`G_j(x)`$ representing $`I+B_0(x)N_j`$, where
$`B_0`$ is the linear field embedding, $`N_j^2=0`$, and padded branches
have $`N_j=0`$. Thus $`G_j(x)^\dagger G_j(y)=G_j(x\oplus y)`$
literally. Let $`J_0`$ prepare the two field bits in states $`|1\rangle`$
and $`|{-}\rangle`$, leaving all other work arbitrary. Its range is
orthogonal to the invariant space at every nonzero parameter. Apply
Hadamards on the eleven clean selector bits and the paid SELECT, then
keep the selector and field outputs as the record. SELECT returns all
private masks exactly, so their identity factors are implicit:

```math
E_x=2^{-11/2}\sum_{j=0}^{2^{11}-1}|j\rangle\otimes G_j(x)J_0,
\qquad E_x^\dagger E_y=J_0^\dagger A_{x\oplus y}J_0.
```

Here A is the lazy field average, so the overlap is I at x=y and has
operator norm at most an absolute $`c\lt1`$ otherwise. Use L copies
with private selectors, field words and masks. The common label is
read-only at every record gate, so all copies align their native T layers.
With the [linear-rank field schedule](GATEWISE_READ_ONLY_EQUALITY.md#linear-rank-field-select),
the record and actual reverse have $`13L`$ clean inputs,
$`T,b_{\rm dirty}=O(L(p+1))`$, $`G=O(L(p+1)^2)`$ and constant
T-depth. Applying the source only on its original wires gives

```math
M'_{xy}=M_{xy}\otimes(E_x^\dagger E_y)^{\otimes L},
\qquad \|M'-M^{(0)}\|\le Bc^L.
```

For the bound, replace off-diagonal blocks within each work sector
by their norms. The resulting real symmetric matrix has row sums at
most $`Bc^L`$, which bounds its operator norm. The sector sum takes
its maximum, including arbitrary work superpositions and references.

Take $`L=O_{B,\gamma}(1+\log(1/\tau))`$ to make $`Bc^L\le\gamma\tau/32`$.
Apply the preceding sign perturbation and exact clean-sector comparison
with $`k=13L`$. The sign error avoids the promotion multiplier, and the
complete scalar costs polynomial overhead in $`1/\tau`$. Here p may be as large as the total source width. At fixed tau the
record uses at most square-root count and width and linear-volume
Clifford count even if every principal dirty wire is tagged. For the
variable-tolerance early stages the extra overhead is polynomial in
$`1/\tau`$, including the multiplicative width factor. The arity
margin and fixed-tolerance full-volume calls give the same frame
consequence as for the linear tags.

For the [adaptive classical counterexample](#adaptive-classical-tag-limit),
the exact four-list sectors have row sum at most $`3/2`$. Its approximate
native source can be chosen within $`O(\gamma\tau)`$ of the ideal
compressed observable using the retained zero-test precision ledger.
Quantum recording contracts the comparison error. Allocate both
that error and $`Bc^L`$ within $`\gamma\tau/32`$ and the same
complete repair follows. For arbitrary tables the signed diagonal
blocks and their row-sum bound remain to be realized.

<a id="two-address-block-margin"></a>

### Two separated address blocks have a width-independent bias bound

Independent input and output preparations do not by themselves supply a
decision source. There is a quantitative obstruction when the entire
read-only source factors into just two address blocks. Let the prescribed
sign table be $`S\in\{\pm1\}^{h\times\ell}`$ and suppose

```math
V_f(u,v)=B_{f,v}A_{f,u}.
```

The two families act on the same signal/work space. They may be arbitrary
unitaries, depend on the entire table f, and leave arbitrary garbage.
Their respective address dependence is only on u and only on v; the
physical address registers remain unchanged. Table-dependent but
address-independent input, middle and output unitaries on signal/work
can be absorbed into these two factors. In particular, the endpoint
unitaries need not be each other's inverses. A word with further
alternating address dependence, a gate depending jointly on u and v,
or physical address motion need not have this factorization.

Let J insert one clean signal $`c=0`$, let D be the full signal/work
dimension, and use $`\tau(X)=\mathrm{Tr}(X)/D`$. Suppress f in the
factors and put

```math
M_{uv}=J^\dagger V_f(u,v)^\dagger Z_cV_f(u,v)J,
\qquad b_{uv}=\frac2D\mathrm{Tr}(M_{uv}).
```

If $`S_{uv}M_{uv}\succeq\gamma I`$ for every address, then
$`S_{uv}b_{uv}\ge\gamma`$. Even the weaker promise that every
computational dirty input has signed bias at least gamma implies this
trace inequality. For either premise,

```math
\gamma\le\frac{\|S\|_{\rm op}}{\sqrt{h\ell}}.
```

**Proof.** Since $`JJ^\dagger=(I+Z_c)/2`$ and the output observable
is traceless, cyclicity gives the exact normalization

```math
b_{uv}=\tau\!\left[(A_uZ_cA_u^\dagger)
                         (B_v^\dagger Z_cB_v)\right].
```

Both factors are Hermitian involutions and therefore unit vectors in
the real normalized Hilbert--Schmidt space. Write their coordinates in
any orthonormal real basis as row vectors $`r_u,t_v`$. The ordinary
matrix operator norm and Frobenius Cauchy--Schwarz inequality give

```math
\begin{aligned}
\gamma h\ell
&\le\sum_{u,v}S_{uv}\langle r_u,t_v\rangle\\
&\le\|S\|_{\rm op}
      \left(\sum_u\|r_u\|^2\right)^{1/2}
      \left(\sum_v\|t_v\|^2\right)^{1/2}\\
&=\|S\|_{\rm op}\sqrt{h\ell}.
\end{aligned}
```

Neither the number of dirty wires nor any native gate count occurs in
this argument. It is a bound on a quantitative margin, not merely on
which strict sign patterns can appear.

For $`a=2r`$, $`m=2^r`$ and the Walsh table
$`S_{uv}=(-1)^{u\cdot v}`$, orthogonality gives
$`SS^{\mathsf T}=mI`$. Consequently

```math
\gamma\le m^{-1/2}=2^{-a/4}.
```

Row or column permutations and sign changes preserve the matrix norm.
Moreover the obstruction concerns far more than this single table.
For a uniformly random $`m\times m`$ sign matrix and any fixed
$`\gamma\gt0`$, the fraction permitting the stated one-clean bias is
at most

```math
\Pr\{\|S\|_{\rm op}\ge\gamma m\}
\le 2\exp\!\left(2m\ln9-\frac{\gamma^2m^2}{8}\right).
```

To see the estimate directly, take a one-quarter net of the real unit
sphere with at most $`9^m`$ points. Maximizing on the two nets loses
at most a factor two in the operator norm. For any two fixed unit
vectors x,y, the independent signs in $`x^{\mathsf T}Sy`$ have
squared coefficients summing to one. The exponential-moment bound
$`\Pr\{|x^{\mathsf T}Sy|\ge t\}\le2e^{-t^2/2}`$ and a union
bound prove the display. Thus only an exponentially small fraction of
tables can have a fixed positive margin in this two-block class. This
probabilistic counting is an analytic statement, not a sampled test.

**Clean inputs and approximate factorization.** With k clean inputs,
write $`P=JJ^\dagger`$, of rank $`D/2^k`$, and retain a fixed
traceless Pauli readout Z. Now

```math
\begin{aligned}
b_{uv}&=\frac{2^k}{D}\mathrm{Tr}(J^\dagger V_f^\dagger ZV_fJ)\\
&=\tau\!\left[A_u(2^kP-I)A_u^\dagger B_v^\dagger ZB_v\right],\\
\tau((2^kP-I)^2)&=2^k-1,\\
\gamma&\le\sqrt{2^k-1}\,
                 \frac{\|S\|_{\rm op}}{\sqrt{h\ell}}.
\end{aligned}
```

The same proof applies because only the first vector norm changes.
For the Walsh table, a constant margin therefore needs
$`k\ge\log_2(1+\gamma^2m)`$ clean inputs in this factorization;
arbitrary additional dirty width does not replace them. The output
measurement is still the one fixed traceless Pauli, not joint
acceptance of all k initialized bits.

If an actual source $`\widetilde V`$ is only approximately factored,
with $`\|(\widetilde V_f(u,v)-B_vA_u)J\|\le\eta`$ uniformly,
the two compressed observables differ by at most $`2\eta`$.
The corresponding bounds therefore hold with $`\gamma-2\eta`$
on the left. This uses initialized-isometry error and includes garbage
and leakage; no full-space approximation outside the initialized code
is needed. The estimate follows by expanding the two products of
isometries around the common readout Z. More generally the actual
source may move addresses if its complete initialized isometry is
within eta of the displayed read-only direct-sum source; the same
compression estimate still applies.

**Application to a single affine-address Boolean phase.** Suppose a
completed diagonal work word has phase $`(-1)^{q_f(x,w)}`$, where
$`q_f`$ is affine over the binary field in the read-only address x,
with arbitrary work-dependent coefficients. Splitting the address into
u and v gives commuting factors
$`D_f(u,v)=D_{f,0}D_{f,L}(u)D_{f,R}(v)`$.
Arbitrary table-dependent signal/work input and output unitaries can
then be absorbed into two factors of the required form. The bound
therefore excludes a universal constant-bias repair by those endpoints
alone. The [paid controlled Gauss phase](GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions)
has exactly this affine-address exponent; its private helpers return
identically and may also be included in the arbitrary endpoints.
This application permits changed, independent boundary codes. It does
not exclude multiple noncommuting phases with table/address operations
between them, or endpoints that move the physical address.

**The anticommutator mechanism attains the one-clean bound.** Let
$`\Gamma_1,\ldots,\Gamma_m`$ be pairwise anticommuting Hermitian
involutions on arbitrary dirty work. For a Walsh table define

```math
R_u=\frac1{\sqrt m}\sum_vS_{uv}\Gamma_v,
\qquad Q_v=\Gamma_v,
\qquad \frac12\{R_u,Q_v\}=\frac{S_{uv}}{\sqrt m}I.
```

The anticommutation relations imply $`R_u^2=Q_v^2=I`$. In signal
block coordinates choose

```math
A_u=\frac1{\sqrt2}
 \begin{pmatrix}I&-R_u\\R_u&I\end{pmatrix},
\qquad
B_v=\frac1{\sqrt2}
 \begin{pmatrix}I&Q_v\\-Q_v&I\end{pmatrix}.
```

These are unitaries, and their complete clean-code compression is

```math
J^\dagger(B_vA_u)^\dagger Z_c(B_vA_u)J
=\frac12\{R_u,Q_v\}
=\frac{S_{uv}}{\sqrt m}I.
```

Thus scalar anticommutators remove dependence on every dirty input but
retain exactly the vanishing normalization. This is an abstract
sharpness construction, not a native synthesis of the $`R_u`$ family.
The preceding bound applies to arbitrary two-block unitaries, not only
to this Majorana realization.

**Three address blocks already escape, with a paid native word.** For
r-bit halves u,v and an arbitrary r-bit bank w, define

```math
A_u|c,w\rangle=|c\mathbin\oplus(u\cdot w),w\rangle,
\qquad B_v|c,w\rangle=|c,w\mathbin\oplus v\rangle.
```

Direct action on every basis column gives

```math
A_uB_vA_u=X_c^{u\cdot v}X_w^v,
\qquad A_uB_vA_uB_v=X_c^{u\cdot v}.
```

The three-block word has full compressed sign $`S_{uv}I`$ and gap
one, with permitted dirty motion. The fourth block returns w. Both
identities preserve all coherent addresses and arbitrary references.
For an actual native realization write
$`A_u=H_c(\prod_i\mathrm{CCZ}(c,u_i,w_i))H_c`$ and implement
each CCZ by the retained eight-T, eight-CNOT Gray phase cycle on its
own arbitrary helper. At each of the eight steps the phase gates
have distinct private targets; shared controls cost no extra T band.
The two A calls use two such batches. The $`B_v`$ calls are r CNOTs
each. Commuting the signal Hadamards through B cancels the middle pair,
so the complete four-block word and its actual native inverse satisfy

```math
T=16r,\qquad D_T\le16,\qquad
G_{\rm Cliff}=18r+2,\qquad w_{\rm dirty}=2r.
```

For the three-block raw decision the Clifford count is $`17r+2`$;
all other displayed bounds are the same. The private helpers always
return exactly, and no initialized bank is assumed. The
[exact fixtures](../../tests/test_separated_decision_margin.py) check
the anticommutator compression and the literal native escape on all
logical and dirty columns. The all-size margin and resource statements
follow from the proofs above.

This closes two separated address blocks as a universal constant-bias
source. It does not close the already alternating interior-Clifford
word below, arbitrary independent boundaries around such a word,
joint table/address interleaving, or general native depth. The paid
Walsh example makes that distinction explicit. The general frame
frontier is unchanged.

<a id="single-projective-group-insertion"></a>

### One projective group insertion has an exact row-bias characterization

Three alternating address blocks allow the Walsh construction above. Consider the form

```math
V_f(u,v)=A'_{f,u}B_{f,v}A_{f,u},
\qquad v\in G=\mathbb F_2^r,
```

where the two outer signal/work unitaries are independent, and the middle family is a projective representation of G. The factors may all depend on the entire table. The address is read-only. Scalar phases of the family cancel in conjugation, so the result includes affine Pauli-frame banks with fixed signal/work changes of basis absorbed into the outer unitaries. One initialized signal c, arbitrary dirty width, and unrestricted garbage are allowed.

Let $`b_u(v)`$ be the normalized trace of the compressed observable $`J^\dagger V_f(u,v)^\dagger Z_cV_f(u,v)J`$ on the initialized code. For a function g of v write $`\widehat g(\chi)=|G|^{-1}\sum_v g(v)\chi(v)`$ for its normalized Fourier coefficients, with $`\chi(v)\in\{\pm1\}`$. Every such program satisfies

```math
\sum_\chi |\widehat b_u(\chi)|\le1\quad\text{for each }u.
```

To prove the bound, let $`\tau(X)=\mathrm{Tr}(X)/D`$ on the D-dimensional signal/work space and suppress f. As above the normalized clean bias is

```math
b_u(v)=\tau(Q_u B_v R_u B_v^\dagger),
\quad R_u=A_u Z_c A_u^\dagger,
\quad Q_u=(A'_u)^\dagger Z_c A'_u.
```

Both observables have unit normalized Hilbert--Schmidt norm. The conjugations by $`B_v`$ are an ordinary representation of G on that Hilbert space. Let $`\Pi_\chi`$ be its orthogonal projection on character $`\chi`$. Fourier transformation and Cauchy--Schwarz give

```math
\widehat b_u(\chi)=\langle \Pi_\chi Q_u,\Pi_\chi R_u\rangle_\tau,
\qquad
\sum_\chi |\widehat b_u(\chi)|\le
\Bigl(\sum_\chi\|\Pi_\chi Q_u\|_\tau^2\Bigr)^{1/2}
\Bigl(\sum_\chi\|\Pi_\chi R_u\|_\tau^2\Bigr)^{1/2}=1.
```

**Complete converse.** Every real row family with Fourier norm at most one can be realized with the exact compressed observable $`b_u(v)I`$, on the entire arbitrary dirty space, using one common projective Pauli family $`B_v`$. For $`r\ge1`$ use $`m+2`$ pairwise anticommuting Hermitian Pauli Majoranas, $`m=|G|`$, one for each character and two dummy axes. These exist on $`(m+2)/2`$ work qubits by the Jordan--Wigner Pauli list. For row u put $`\beta_\chi=\widehat b_u(\chi)`$, $`L=\sum_\chi|\beta_\chi|`$, and take

```math
\begin{aligned}
R&=\sum_\chi\sqrt{|\beta_\chi|}\,\Gamma_\chi+\sqrt{1-L}\,\Gamma_*,\\
Q&=\sum_\chi\mathrm{sgn}(\beta_\chi)\sqrt{|\beta_\chi|}\,\Gamma_\chi+\sqrt{1-L}\,\Gamma_{**}.
\end{aligned}
```

Both are Hermitian involutions. Choose $`B_v`$ so that $`B_v\Gamma_\chi B_v^\dagger=\chi(v)\Gamma_\chi`$ and both dummy axes are fixed. On an even-sized Majorana frame a product over an even subset flips exactly that subset under conjugation. To flip an odd subset, use its odd complement. Apply this rule to the character signs of each standard generator of G, making a Hermitian Pauli representative $`P_i`$ for its action. Set $`B_v=P_1^{v_1}\cdots P_r^{v_r}`$. Pauli multiplication makes this a projective representation. Each controlled tensor Pauli is a Clifford circuit of $`O(m)`$ elementary gates, using a control-Z for a negative sign and Clifford basis changes for its nonidentity local factors. The shared middle word and its actual inverse therefore have zero T gates and $`O(rm)`$ elementary Cliffords. They use no additional helper beyond the Majorana work bank.

On c and this work choose the outer blocks

```math
A_u=\frac1{\sqrt2}\begin{pmatrix}I&-R\\R&I\end{pmatrix},\qquad
A'_u=\frac1{\sqrt2}\begin{pmatrix}I&Q\\-Q&I\end{pmatrix}.
```

Multiplying their two clean-signal blocks gives the exact compression $`\{B_v^\dagger Q B_v,R\}/2=b_u(v)I`$. The two distinct dummy axes contribute zero to the anticommutator. The complete converse is an abstract source characterization with a paid middle word. The outer blocks need their own native synthesis when applied to a proposed compiler.

**Margin.** For a prescribed row $`s_u(v)\in\{\pm1\}`$, the best margin in this architecture is

```math
\max_{\|\beta\|_1\le1}\min_v s_u(v)\sum_\chi\beta_\chi\chi(v)
=\min_{\mu\in\mathcal P(G)}\max_\chi
     \left|\sum_v\mu(v)s_u(v)\chi(v)\right|..
```

The equality is finite bilinear minimax over the real unit $`\ell_1`$ ball and the probability simplex, followed by $`\ell_1`$/$`\ell_\infty`$ duality. Necessity uses normalized trace after a full signed operator gap, and sufficiency uses the complete converse. The same middle word works for every row, so the table margin is the minimum of the displayed value over rows. Uniform $`\mu`$ gives $`\gamma\le\max_\chi |\widehat s_u(\chi)|`$. At even r the bent row $`s_u(v)=(-1)^{v_1v_2+\cdots+v_{r-1}v_r}`$ has every Fourier magnitude $`2^{-r/2}`$, by multiplying the two-bit Fourier transforms. This bound is attained abstractly by $`b_u=2^{-r/2}s_u`$, whose Fourier $`\ell_1`$ norm is one. In contrast each Walsh row is a character, and the previously paid Walsh circuit attains gap one. One group-like visit to v thus does not suffice for arbitrary tables. Test a non-group-like family or separated v visits to leave this restriction.

<a id="paired-pauli-one-commuting-wall"></a>

### Paired affine Pauli visits around a commuting T wall have a quadratic margin bound

Consider the changed read-only source

```math
V_f(u,v)=A'_{f,u} B_v^\dagger W_f B_v A_{f,u},
```

with arbitrary independent u-boundaries, B a projective Pauli representation of G, and W a fixed product of commuting Pauli T or T-adjoint rotations. It may depend on the whole table. Allow any number of commuting axes, even dependent ones, in W. A fixed Clifford frame for all those axes is included by changing B and the outer boundaries. The architecture has two v visits with a non-Clifford interlayer, so the single-insertion theorem need not apply.

For an axis $`P_i`$ and sign $`\epsilon_i\in\{\pm1\}`$ write the literal rotation as $`T(P_i)^{\epsilon_i}=e^{i\epsilon_i\pi/8}e^{-i\epsilon_i\pi P_i/8}`$, and let $`\ell_i(v)`$ be the commutation bit of $`B_v`$ with $`P_i`$. It is affine linear as a Boolean function, and linear when the representation is normalized at zero. Since all axes commute,

```math
E_v=B_v^\dagger W_f B_v W_f^\dagger
=\prod_i\exp\bigl(i\epsilon_i\pi\ell_i(v)P_i/4\bigr).
```

Absorb the fixed W into the initial boundary. In one common simultaneous axis eigenbasis, the adjoint action of E on $`|s\rangle\langle t|`$ is multiplication by $`i^{q(v)}`$, where

```math
q(v)=\sum_i \epsilon_i\,\frac{p_i(s)-p_i(t)}2\,\ell_i(v)\pmod4
   =q_0+\sum_j a_j v_j+2\sum_{j\lt k}b_{jk}v_jv_k\pmod4.
```

Here $`p_i(s),p_i(t)\in\{\pm1\}`$ are the axis eigenvalues, $`a_j\in\mathbb Z_4`$, $`b_{jk}\in\mathbb F_2`$, and $`q_0\in\mathbb Z_4`$. The second equality follows by expanding each Boolean parity as an integer modulo four. Group eigenoperators by their common phase function into mutually orthogonal spaces. The normalized trace proof above now gives $`b_u(v)=\sum_q\beta_{u,q}i^{q(v)}`$ with $`\sum_q|\beta_{u,q}|\le1`$, possibly complex coefficients. Consequently a full signed margin on a row forces

```math
\gamma\le \max_q\bigl|\mathbb E_v s_u(v)i^{q(v)}\bigr|.
```

This condition is independent of width, number of commuting axes and the boundary unitaries. A uniform random sign row has correlation at most $`O(r/2^{r/2})`$ with all these phases with probability tending to one. More explicitly the probability that the maximum exceeds $`t`$ is at most $`4\cdot4^{r+1}\,2^{\binom r2}\exp(-2^r t^2/4)`$, by applying the independent-sign exponential bound to the real and imaginary parts and counting the possible q.

There is also an explicit adversarial row. Set
$`s(v)=(-1)^{\sum_j v_{3j-2}v_{3j-1}v_{3j}}`$ using disjoint triples.
Every phase $`i^{q(v)}`$ has third multiplicative derivative one.
Binary addition modulo four proves this for each linear term. Twice
each quadratic term has this property because its third Boolean
derivative vanishes. Three Cauchy--Schwarz applications therefore
bound the correlation by $`\|s\|_{U^3}`$, where

```math
\|s\|_{U^3}^8=\mathbb E_{v,h_1,h_2,h_3}
 \prod_{\omega\in\{0,1\}^3}\mathcal C^{|\omega|}
 s\left(v+\sum_j\omega_j h_j\right)
```

and $`\mathcal C`$ is complex conjugation. For each triple, the cubic
derivative is the determinant over $`\mathbb F_2`$ of the three
direction columns. Exactly $`(8-1)(8-2)(8-4)=168`$ of the $`2^9`$
direction matrices are nonsingular. Independence across triples gives

```math
\gamma\le(1-2\cdot168/512)^{\lfloor r/3\rfloor/8}
       =(11/32)^{\lfloor r/3\rfloor/8}.
```

**Arbitrary middle with one-band endpoints.** The same quadratic margin bound holds for arbitrary fixed W, even with many noncommuting layers, if each outer u-boundary has at most one native T band on signal/work for each fixed u. Fixed Clifford frames can be included without changing this conclusion. For one such endpoint write $`R=D P D^\dagger`$, with P a Pauli and D a commuting product of Pauli T rotations. An outer Clifford conjugates both P and the rotation axes and is therefore already included in this expression. A direct multiplication on the axes that anticommute with P gives

```math
B_vR B_v^\dagger=\chi_P(v)G_v R,\qquad
G_v=\prod_{i:\{P_i,P\}=0}(i\epsilon_iP_i)^{\ell_i(v)}.
```

Here $`\chi_P`$ is the character of the Pauli P under conjugation by $`B_v`$. The $`G_v`$ commute for all v because their axes belong to the same commuting band. Likewise write $`B_v Q B_v^\dagger=\chi_{P'}(v)H_vQ`$, with the family $`H_v`$ commuting in its own basis. The bias with the arbitrary middle W is therefore

```math
\begin{aligned}
b_u(v)&=\tau\bigl((B_vQB_v^\dagger)W(B_vRB_v^\dagger)W^\dagger\bigr)\\
 &=\chi_P(v)\chi_{P'}(v)\tau(RW^\dagger H_v Q W G_v)\\
 &=\chi_P(v)\chi_{P'}(v)\langle WR,H_v(QW)G_v\rangle_\tau.
\end{aligned}
```

Left multiplication by H and right multiplication by G commute even when their bases differ. On $`|h_s\rangle\langle g_t|`$, the simultaneous eigenphases are products of powers $`i^{\pm\ell_i(v)}`$. Including the two characters, they are exactly of the generalized quadratic form above. Both fixed Hilbert--Schmidt vectors have norm one, so the same orthogonal-space argument and cubic row apply. This proof places no native count or depth restriction on W.

In fact, for arbitrary middle W, a single boundary that is a character vector under $`\mathrm{Ad} B_v`$ already forces the Fourier $`\ell_1`$ bound above. Factor its character out of the first line above and use the one remaining group matrix coefficient. This includes one Clifford endpoint paired with any other endpoint. To escape both bounds while retaining the same paired Pauli address frames, both probes must have more than one character and at least one must leave the one-band architecture. The remaining source still needs a complete compressed sign after those changes.

**A paid escape for the cubic row.** For the explicit disjoint-triple cubic $`f(v)`$ above, its polarization has the form $`f(k+z)=f(k)+f(z)+\theta(k)\cdot\phi(z)`$. For each triple, use its six mixed terms, writing each as a product of one monomial in k and one in z, both of degree at most two. Borrow arbitrary registers k and t, with t one bit per mixed term. Write $`D_{c g}`$ for the diagonal phase $`(-1)^{c g}`$. Use

```math
\begin{aligned}
B_v&=X_k^v,\qquad A=D_{c[f(k)+t\cdot\phi(k)]}H_c,\\
A'&=H_cD_{c[t\cdot\phi(k)]},\\
W|c,k,t\rangle&=(-1)^{c f(k)}|c,k,t+\theta(k)\rangle.
\end{aligned}
```

The two probes are $`R=X_c(-1)^{f(k)+t\cdot\phi(k)}`$ and
$`Q=X_c(-1)^{t\cdot\phi(k)}`$. Conjugating $`B_vRB_v^\dagger`$ by W
adds $`f(k)`$ from the signal phase and replaces t by $`t+\theta(k)`$.
The resulting exponent satisfies

```math
 f(k+v)+f(k)+(t+\theta(k))\cdot\phi(k+v)=f(v)+t\cdot\phi(k+v).
```

Thus

```math
W(B_vRB_v^\dagger)W^\dagger=(-1)^{f(v)}B_vQB_v^\dagger,
\qquad
(A'B_v^\dagger W B_v A)^\dagger Z_c (A'B_v^\dagger W B_v A)
 =(-1)^{f(v)}Z_c.
```

This is a full physical operator identity with only the eventual readout signal initialized. The $`t+\theta(k)`$ update uses affine gates and Toffolis on distinct target bits. The other table-dependent terms are diagonal phases of arity at most four. Use the retained read-only private-wire phase batches on disjoint arbitrary helpers. They leave all principal controls available in Clifford interactions, including the shared c, so the T layers align across triples. The forward word and actual reverse use $`O(r)`$ T gates, elementary Cliffords and dirty work at constant native T-depth with zero implementation error. The raw source may move its dirty t bank. This restricted construction escapes the preceding architecture bounds using a changed middle and changed endpoints. General polarization incurs the large-degree monomial-source cost and gives no general frame-depth improvement.

<a id="fixed-boundary-decision-source"></a>

**Fixed boundaries cannot decode every table at square-root width.**
The code-preserving gauge check above leaves code-moving boundaries
open. A separate information bound applies when all table dependence
is confined to one unitary between fixed boundaries. Suppose the source
preserves the address and its fibers on w physical signal/work qubits
are

```math
 V_f(x)=F^\dagger G_x W_f G_x^\dagger F,
 \qquad x\in\{0,1\}^a.
```

F and every $`G_x`$ are independent of f. They may be arbitrary
unitaries and may move the clean-zero code. The table-dependent
$`W_f`$ is also arbitrary. Let $`D=2^w`$, let J insert the one
clean signal, and define

```math
 \begin{aligned}
 S&=Z_c,\qquad S_x=G_x^\dagger F S F^\dagger G_x,
 \qquad P_x=\frac{I+S_x}{2},\\
 M_f(x)&=J^\dagger V_f(x)^\dagger S V_f(x)J,
 \qquad \tau_w(A)=\frac{\mathrm{Tr}(A)}D.
 \end{aligned}
```

Assume that every table has the coherent decision promise
$`(-1)^{f(x)}M_f(x)\succeq\gamma I`$ at every address, for a
common positive gap gamma. The normalized compressed bias satisfies

```math
 \begin{aligned}
 b_f(x)&=\frac2D\mathrm{Tr}(M_f(x))\\
 &=\frac2D\mathrm{Tr}(P_x W_f^\dagger S_x W_f)\\
 &=\tau_w(W_f^\dagger S_x W_f S_x),
 \qquad (-1)^{f(x)}b_f(x)\ge\gamma.
 \end{aligned}
```

The last equality uses $`\mathrm{Tr}(S_x)=0`$. It does not need
the balanced-spectrum argument or a prescribed rejected-sector sign.
For the normalized Choi vector

```math
 |\Phi_D\rangle=\frac1{\sqrt D}\sum_{j=0}^{D-1}|j\rangle|j\rangle,
 \qquad |W_f\rangle=(I\otimes W_f)|\Phi_D\rangle,
```

the fixed Hermitian involution $`S_x^{\mathsf T}\otimes S_x`$
has expectation

```math
 \langle W_f|S_x^{\mathsf T}\otimes S_x|W_f\rangle=b_f(x).
```

Its sign outcome recovers bit $`f(x)`$ with error at most
$`(1-\gamma)/2`$. Thus the $`2^a`$ independently prescribed table
bits have a quantum random-access encoding in $`2w`$ qubits.
[F54](../../docs/reference/SOURCE_CATALOGUE.md),
[Nayak's Theorem 2.3](https://arxiv.org/abs/quant-ph/9904093), gives

```math
 2w\ge\alpha(\gamma)2^a,
 \qquad \alpha(\gamma)=1-H_2\!\left(\frac{1-\gamma}{2}\right),
 \qquad
 \gamma\le\sqrt{\frac{4\ln(2)w}{2^a}}.
```

The last inequality follows from
$`\alpha(\gamma)\ge\gamma^2/(2\ln2)`$. A constant positive gap
therefore requires $`w=\Omega(2^a)`$, irrespective of the gate count
or depth of the three constituent unitaries. At
$`w=O(2^{a/2})`$ the common gap is at most $`O(2^{-a/4})`$.
The Choi state and its measurement are proof objects, not assumed
physical resources of the source. The
[non-Clifford seed-spread screen](NONCLIFFORD_SEED_SPREAD.md#3-a-robust-bound-from-the-independent-commuting-half)
uses the same imported information bound on a compressed Pauli algebra;
the present application encodes the table in the program unitary and
does not assume Pauli label degree or a small support algebra.

**A table-dependent outer Clifford can still have a short descriptor.**
Suppose instead that $`G_x=\prod_{j=1}^a Q_j^{x_j}`$ for fixed
commuting Pauli generators and $`F_f=L_fF_0`$, where $`F_0`$ is
fixed and $`L_f`$ is an arbitrary table-dependent Clifford. Set

```math
 Q'_{f,j}=L_f^\dagger Q_jL_f,\qquad
 G'_{f,x}=\prod_{j=1}^a(Q'_{f,j})^{x_j},\qquad
 W'_f=L_f^\dagger W_fL_f.
```

Then $`V_f(x)=F_0^\dagger G'_{f,x}W'_fG_{f,x}^{\prime\dagger}F_0`$.
The signed Pauli images have a description of at most
$`a(2w+2)`$ classical bits, including their phases. Append this
description as a classical register to $`|W'_f\rangle`$. A fixed
decoder reads that register and uses the corresponding
$`S_{f,x}^{\mathsf T}\otimes S_{f,x}`$ measurement. The descriptor
is counted in the encoding, rather than supplied to a table-dependent
decoder for free. The same theorem therefore gives

```math
 2w+a(2w+2)\ge\alpha(\gamma)2^a.
```

At constant gap this implies $`w=\Omega(2^a/(a+1))`$ and still
excludes square-root width for this specified outer-Clifford family.

**Application to a paid twelve-band table word.** Let $`a=2r`$,
$`m=2^r`$, and index two dirty banks A and B by
$`i,j\in\{0,1\}^r`$. Write
$`Z_u^i=\prod_kZ_{u_k}^{i_k}`$ and similarly for v. For an arbitrary
binary table C, define

```math
 \begin{gathered}
 A_u=\prod_i e^{-i\pi X_cX_{A_i}Z_u^i/8},\qquad
 B_v=\prod_j e^{-i\pi Y_cX_{B_j}Z_v^j/8},\\
 E_C=\prod_{i,j}\mathrm{CZ}_{A_i,B_j}^{C_{ij}},\qquad
 F_C=L_CA_0B_0,\\
 K_C=A_u^\dagger(E_CB_v^\dagger E_C)A_u(E_CB_vE_C),\qquad
 K_0=A_u^\dagger B_v^\dagger A_uB_v,\\
 V_C=F_C^\dagger K_CK_0^\dagger F_C.
 \end{gathered}
```

Each pulse has m commuting independent Pauli axes and is one native T
band after Clifford transport. Actual conjugated T gates introduce a
common scalar which cancels with the actual inverses in the displayed
word. For Clifford $`L_C`$ the source uses $`12m`$ T gates, at most
twelve T bands, two m-qubit dirty banks, and one clean signal. Its
Clifford count is $`O(m^2+G_L)`$, where $`G_L`$ bounds either
native direction of $`L_C`$; thus it meets the proposed
$`T,w=O(m)`$, $`G=O(m^2)`$ budgets when $`G_L=O(m^2)`$.

Nevertheless the bank-Z gauge

```math
 G_{uv}=\prod_i Z_{A_i}^{i\cdot u}\prod_j Z_{B_j}^{j\cdot v}
```

conjugates $`A_0,B_0`$ to $`A_u,B_v`$ and commutes with
$`E_C`$. Hence $`K_CK_0^\dagger=G_{uv}W_CG_{uv}^\dagger`$
for an address-independent $`W_C`$. The gauge has a fixed commuting
Pauli generator for each of the a address bits. The fixed-boundary
bound applies when $`L_C`$ is table independent, and its counted
descriptor corollary applies to an arbitrary outer Clifford
$`L_C`$. Neither case supplies a uniform constant-gap decision for
every C, despite the paid native gate budgets.

These bounds concern the displayed read-only factorization. More
general table-dependent boundaries not covered by the short-descriptor
corollary, a table-dependent Clifford interleaved inside non-Gaussian
code preparation, further table/address interleaving, and physical
address motion remain outside their scope. They give no unrestricted
query-width or T-depth lower bound.

<a id="interior-clifford-decision-source"></a>

### An interior Clifford can give table-sized decoder information

Keep the twelve-band middle word above, but replace its preparation by
$`F_C=A_0L_CB_0`$, with the rightmost operator acting first. This
retains the displayed $`12m`$ T-count, twelve-band upper bound and
$`O(m^2+G_L)`$ Clifford ledger. All native pulse scalars cancel with
their actual inverses. Put $`w=2m+1`$ and define

```math
 \begin{aligned}
 \Gamma_{C;u,v}&=L_C^\dagger A_0^\dagger G_{uv}A_0L_C,\\
 W'_C&=L_C^\dagger A_0^\dagger W_CA_0L_C,
 \qquad D_0=B_0Z_cB_0^\dagger.
 \end{aligned}
```

The exact source and transformed Choi signal are

```math
 V_C=B_0^\dagger\Gamma_{C;u,v}W'_C\Gamma_{C;u,v}^\dagger B_0,
 \qquad S_C(u,v)=\Gamma_{C;u,v}^\dagger D_0\Gamma_{C;u,v}.
```

Thus $`|W'_C\rangle`$ is decoded by
$`S_C(u,v)^{\mathsf T}\otimes S_C(u,v)`$. For
$`p_i^0=X_cX_{A_i}`$ and $`U(u)=\{i:i\cdot u=1\}`$,

```math
 A_0^\dagger G_{uv}A_0
   =G_{uv}\prod_{i\in U(u)}e^{-i\pi p_i^0/4}.
```

The address generators have become Clifford involutions, generally
not Paulis. A sufficient common descriptor records the signed images
under $`L_C`$ of the $`m-1`$ axes with $`i\ne0`$ and the a
Pauli generators of $`G_{uv}`$. It costs at most
$`(m-1+a)(2w+2)`$ classical bits. The zero-index axis is never
selected. At $`w=O(m)`$ this descriptor can already have table-sized
$`O(m^2)`$ information, so the earlier short-descriptor bound does
not directly exclude this preparation.

**The larger dependence survives in the actual measurement family.**
For $`L_C=E_C`$, set

```math
 p_i=E_Cp_i^0E_C=X_cX_{A_i}Z_B^{C_{i,*}},\qquad
 q_j=Y_cX_{B_j},\qquad d_{ij}=1+C_{ij}\pmod2.
```

At $`v=0`$, abbreviate $`G_u=G_{u0}`$. Then
$`\Gamma_{C;u,0}=G_u\prod_{i\in U(u)}e^{-i\pi p_i/4}`$
commutes with every $`Z_{B_j}`$. Direct conjugation gives

```math
 D_0Z_{B_j}D_0=-i q_jZ_{B_j}=-Y_cY_{B_j}.
```

The $`p_i`$ commute with each other, while
$`p_iq_j=(-1)^{d_{ij}}q_jp_i`$ and $`G_u`$ commutes with
$`q_j`$. With the ordinary integer
$`s_j(u)=\sum_{i\in U(u)}d_{ij}`$, the literal Pauli identity is

```math
 S_C(u,0)Z_{B_j}S_C(u,0)
   =-i\,i^{s_j(u)}
      \left(\prod_{\substack{i\in U(u)\\d_{ij}=1}}p_i\right)
      q_jZ_{B_j}.
```

Its A-bank X support is exactly
$`\{i:i\cdot u=1,\ C_{ij}=0\}`$: no other factor acts on A,
and the distinct $`X_{A_i}`$ factors cannot cancel. Choosing
$`u=e_\ell`$ with $`i_\ell=1`$ recovers every entry of every
row $`i\ne0`$. Row zero enters none of these gauges, including
those with nonzero v, and changes none of the measurements. Therefore
the normalized measurement families have exactly
$`2^{m(m-1)}`$ distinct values. Their exact fixed-length information
content is $`m(m-1)`$ bits. This distinction also holds for the
two-register Choi measurements: they identify S up to a scalar sign,
and $`SZ_{B_j}S`$ is unchanged by that sign.

This concerns exact descriptions of this specified decoder family.
It is not a bound on approximate descriptions, source width, or all
possible Choi normalizations. It establishes a genuine failure of
the previous short-descriptor argument in this normal form; it does
not establish a coherent decision gap. The next two tests address
that separate requirement.

<a id="paired-row-decision-obstruction"></a>

### Paired rows leave an exact wrong-sign subspace

Let $`r\ge1`$ and suppose C has identical paired rows,

```math
 C_{(k,0),j}=C_{(k,1),j}=D_{kj}.
```

Put one singlet on each A pair, leaving the entire B bank arbitrary:

```math
 |\psi^-\rangle=\frac{|01\rangle-|10\rangle}{\sqrt2},\qquad
 |\omega_A\rangle=\bigotimes_k|\psi^-\rangle_{A_{k,0}A_{k,1}},
 \qquad
 \mathcal D_A=\mathcal H_c\otimes
     \mathrm{span}\{|\omega_A\rangle\}\otimes\mathcal H_B.
```

The identity $`(X_1+X_2)|\psi^-\rangle=0`$ implies that
$`A_u`$ acts as identity on this subspace whenever $`u=(p,0)`$:
the two character signs in each pair coincide. This also applies
to $`A_0`$. Every computational component of a singlet has odd
pair parity, so $`E_C`$ preserves the subspace and restricts to
$`Z_B^\kappa`$, where $`\kappa_j=\sum_kD_{kj}\pmod2`$.
Every B pulse preserves it as well. Consequently, for every v,

```math
 (K_CK_0^\dagger)|_{\mathcal D_A}=I
 \qquad\text{when }u=(p,0).
```

Indeed, the A pulse is identity there and the remaining factors
in each commutator are actual inverses.

The subspace $`\mathcal D_A`$ is precisely the simultaneous
minus-one eigenspace of the swaps exchanging the two wires in each
A pair. Any boundary commuting with those swaps preserves it.
The fixed $`A_0,B_0`$ commute with the swaps, so the interior
preparation has this property whenever $`L_C`$ does. Hence, for
every B-bank vector beta, its clean compression obeys

```math
 M_C((p,0),v)(|\omega_A\rangle|\beta\rangle)
       =|\omega_A\rangle|\beta\rangle,
 \qquad M_C=J^\dagger V_C^\dagger Z_cV_CJ.
```

For any nonzero D choose $`D_{pj}=1`$ and set
$`u=(p,0),v=j`$. The prescribed sign is negative, whereas the
displayed positive eigenspace has dirty dimension $`2^m`$.
Thus every positive coherent gap fails. This entangled dirty input
is an allowed witness, not an initialized circuit resource.

The condition holds for $`L_C=E_C`$ and both CNOT table shears

```math
 \Sigma_C^{A\to B}:|a,b\rangle\longmapsto|a,b+C^{\mathsf T}a\rangle,
 \qquad
 \Sigma_C^{B\to A}:|a,b\rangle\longmapsto|a+Cb,b\rangle.
```

Identical rows make the first shear invariant under exchanging the
two A controls; the second adds the same bit to both A targets.
Both commute with every paired swap. No column duplication is needed.
The column version uses B singlets and last-bit-zero v, with arbitrary u.

For the symmetric case $`C=D\otimes J_2`$, B singlets may also
be imposed. The graph phase and both shears then have literal scalar
$`(-1)^{\sum_{k,l}D_{kl}}`$ on the two-dimensional signal subspace:
each pair has odd Z parity, and $`X\otimes X`$ has singlet
eigenvalue minus one. That scalar cancels with its actual inverse.
The native T-pulse scalars cancel in the same forward/inverse pairs,
so the resulting identity is literal. These preserved-subspace
arguments do not exclude an arbitrary boundary that moves the code
away from the exhibited subspaces.

<a id="explicit-symmetry-breaker-failure"></a>

### A paid symmetry breaker avoids those subspaces but still fails

Use the first pair in each bank and define the constant-size Clifford
$`L_*=U_1U_0`$, again with the rightmost gate first:

```math
 \begin{aligned}
 U_0&=\mathrm{CX}_{c\to B_0}H_{B_0}H_cS_c^\dagger,\\
 U_1&=Z_{A_0}Z_{B_0}\mathrm{CX}_{A_0\to A_1}
       \mathrm{CX}_{B_0\to B_1}H_{A_0}H_{B_0}
       \mathrm{SWAP}_{c,A_0}.
 \end{aligned}
```

An index inside a gate or Pauli subscript denotes a bank qubit;
the standalone $`B_0`$ still denotes the zero-address pulse.
Z and SWAP have their usual fixed H,S,CNOT expansions, so this
boundary adds no T gates, clean qubits or growing Clifford cost.
Literal Pauli transport gives

```math
 L_*Y_cL_*^\dagger=-X_{A_0}X_{A_1},\qquad
 L_*(Y_cX_{B_0})L_*^\dagger=-X_{B_0}X_{B_1}.
```

Both input Paulis anticommute with $`D_0=B_0Z_cB_0^\dagger`$.
Their images stabilize the corresponding A and B singlet subspaces.
Thus the positive eigenspace of $`L_*D_0L_*^\dagger`$ has zero
intersection with either subspace. The preceding $`A_0`$ preserves
both. For repeated-block C, $`E_C`$ also commutes with the two
displayed stabilizers. Consequently the preparation with
$`L_C=E_CL_*`$ has a clean code with zero intersection with
both previously exhibited subspaces.

Nevertheless, take $`m=2,C=J_2,u=v=0`$. Propagating the Pauli
observable through the specified native word gives exactly

```math
 V_C^\dagger Z_cV_C
  =Z_c\frac{3I-X_{B_0}}4
    +X_c\frac{I+X_{B_0}}4(Q_1+Q_2-Q_3),
```

where

```math
 \begin{aligned}
 Q_1&=Z_{A_0}Y_{A_1}\frac{Z_{B_1}+Y_{B_1}}{\sqrt2},\\
 Q_2&=X_{A_0}X_{A_1}\frac{Z_{B_1}-Y_{B_1}}{\sqrt2},\\
 Q_3&=Y_{A_0}Z_{A_1}X_{B_1}.
 \end{aligned}
```

The three Q operators are pairwise anticommuting Hermitian involutions;
this also checks that the displayed reflection squares to identity.
Its twelve Pauli terms comprise two with $`Z_c`$ and ten with
$`X_c`$. Clean-zero compression therefore yields

```math
 M_C(0,0)=\frac{3I-X_{B_0}}4.
```

Its eigenvalues are one half and one. Every computational dirty basis
input to this five-qubit word has positive bias three quarters,
although the selected table entry is one. The whole active compression
has the wrong sign, beyond the previously identified singlet witnesses.

For every larger power-of-two m, put $`J_2`$ in the first two rows
and columns of C and set all other entries to zero. Put the additional
A and B pairs in singlets. At address zero their X sums vanish;
$`E_C`$ and $`L_*`$ act only on the active first pairs. This is an
invariant subspace on which the entire source restricts to the
five-qubit word above. Its compression and wrong sign persist at
every size; for $`m\ge4`$ these are nonconstant tables. Spectator
native pulse scalars cancel with the actual inverses as before.

The [interior-decision fixtures](../../tests/test_interior_decision.py)
check the finite operator identities; the explicit recovery and
spectator arguments above supply their arbitrary-size scope.
This paid repair shows that avoiding the known singlet subspaces
is insufficient, including on the embedded nonconstant tables.
By itself it does not exclude every other boundary. A separate
spectral test below treats all conjugating boundaries for a different,
constant-table family.

<a id="boundary-independent-spectral-obstruction"></a>

### An all-ones table defeats every conjugating boundary of this middle word

The literal-table prescription $`V_C=F_C^\dagger K_CK_0^\dagger F_C`$
cannot give a correct one-clean decision for every C, even if its
conjugating boundary is an arbitrary unitary. The obstruction is the
all-ones table $`C=J_m`$ at every $`m=2^r`$, $`r\ge1`$.
This concerns that prescribed middle word and its mutually inverse
boundaries; separate treatment or preprocessing of the table is not
being excluded.

**A spectral arc and a half-dimensional code.** Let Q be a unitary
on an even-dimensional space of dimension d. Suppose a spectral
subspace of dimension greater than $`d/2`$ lies in an arc
$`\{e^{it}:\phi\le t\le\phi+\omega\}`$, where
$`0\le\omega\le\pi`$. Every rank-$`d/2`$ projector P has a
unit vector psi in the intersection of its range with that subspace.
Rotating the arc to its midpoint gives

```math
 |\langle\psi,Q\psi\rangle|\ge\cos(\omega/2),
 \qquad
 \|PQ\psi\|\ge|\langle\psi,Q\psi\rangle|.
```

Consequently

```math
 \langle\psi|Q^\dagger(2P-I)Q|\psi\rangle
   =2\|PQ\psi\|^2-1\ge\cos\omega.
```

For one clean signal, $`JJ^\dagger=(I+Z_c)/2`$ has rank
$`d/2`$. For any boundary F take $`P=FJJ^\dagger F^\dagger`$;
the compressed observable of $`V=F^\dagger QF`$ is the compression
of the preceding operator to that code, expressed through the
isometry $`FJ`$. Hence

```math
 \lambda_{\max}(J^\dagger V^\dagger Z_cVJ)\ge\cos\omega.
```

**The all-ones middle word has such an arc.** For $`C=J_m`$,
every $`p_i=X_cX_{A_i}`$ commutes with every transported axis

```math
 E_C(Y_cX_{B_j})E_C=Y_cX_{B_j}Z_A^{\rm all}.
```

The two Pauli anticommutations cancel. Address-dependent character
signs do not change this relation, so $`K_C=I`$ literally and
$`Q=K_0^\dagger`$.

At address zero, simultaneous X-bank eigenvalues
$`\epsilon_i,\eta_j\in\{\pm1\}`$ reduce the two pulses to
$`A=e^{-i\alpha X_c}`$ and $`B=e^{-i\beta Y_c}`$, where

```math
 \alpha=\frac\pi8\sum_i\epsilon_i,\qquad
 \beta=\frac\pi8\sum_j\eta_j.
```

The signal-qubit commutator $`A^\dagger B^\dagger AB`$ has
eigenvalues $`e^{\pm i\theta}`$, with

```math
 \cos\theta=1-2\sin^2\alpha\sin^2\beta.
```

For a uniform bank eigenassignment let $`s=\sum_i\epsilon_i`$
and $`p=s/2`$, an integer because m is even. Its residue probabilities
are

```math
 \begin{aligned}
 q_0&=\Pr(p\equiv0\bmod4)=\frac14+\delta,\\
 q_2&=\Pr(p\equiv2\bmod4)=\frac14-\delta,\\
 \Pr(p\text{ odd})&=\frac12,\qquad \delta=2^{-m/2-1}.
 \end{aligned}
```

For example, these follow from
$`\mathbb E e^{i\pi p/2}=2^{-m/2}`$ and
$`\mathbb E(-1)^p=0`$, together with the conjugate Fourier moment.
The two banks are independent in this spectral count. The fraction
with eigenvalue one is $`2q_0-q_0^2`$. If both residues are odd,
the two eigenvalues are $`e^{\pm i\pi/3}`$; these blocks have
fraction one quarter, contributing one eighth to each eigenvalue.
Taking the adjoint preserves those fractions. Therefore Q has a
spectral subspace with eigenphases in $`[0,\pi/3]`$ and fractional dimension

```math
 2q_0-q_0^2+\frac18
   =\frac9{16}+\frac32\delta-\delta^2>\frac12.
```

The lemma yields, for every conjugating boundary F,

```math
 \lambda_{\max}(J^\dagger V^\dagger Z_cVJ)\ge\frac12.
```

The all-ones table requires a negative definite compression at any
positive decision gap, so no such boundary succeeds. Actual native
pulse scalars cancel within each commutator, and the boundary and
its actual inverse cancel their common scalar; no projective change
of the spectrum was made.

Extra arbitrary dirty registers on which the middle word is identity
preserve all spectral fractions. At any other address, the Walsh
character signs merely permute the independent bank eigenassignments,
so every address block has the same spectrum. The lemma thus also
applies to the complete address-register direct sum and permits F
to move the address: the required table sign is negative identity
on that entire register. The one-clean, rank-half premise and the
unchanged middle word remain essential to this statement.

The [spectral fixtures](../../tests/test_interior_decision.py) check
the small native middle word and residue-count identities; the arc
and intersection arguments give the arbitrary-boundary, all-size
conclusion. The earlier paired-row and explicit-repair results retain
their separate nonconstant counterfamilies. This stronger obstruction
excludes a universal prescription using the literal table in this
fixed conjugated middle word. It does not exclude separate compilation
of constant tables, table complementation or other preprocessing with
changed surroundings, independent input/output boundaries, or a
changed middle word. Other asymmetric boundaries for the remaining
nonconstant tables remain open. No unrestricted query-width or
growing T-depth lower bound follows.

<a id="calibrated-asymmetric-decision-source"></a>

### Independent code boundaries still have a repeated-row response

A concrete asymmetric alternative keeps the paid middle
$`Q_C=K_CK_0^\dagger`$ but uses different input and output preparations:

```math
 F_{\mathrm{in},C}=A_0E_CB_0,\qquad
 F_{\mathrm{out},C}=B_0E_CA_0,\qquad
 U_C=F_{\mathrm{out},C}^\dagger Q_C F_{\mathrm{in},C}
         F_{\mathrm{in},0}^\dagger F_{\mathrm{out},0} .
```

The rightmost calibration makes the zero-table source exactly identity.
It does not make the two table-dependent boundaries actual inverses, so
the preceding conjugation-only spectral-arc argument does not apply.
Cancelling the adjacent B0 and its actual inverse gives

```math
 U_C=A_0^\dagger E_CB_0^\dagger Q_C
          A_0E_CA_0^\dagger B_0A_0 .
```

Both native directions have at most $`14m`$ T gates and fourteen T bands,
$`O(m^2)`$ elementary Cliffords, two m-bit dirty banks and one clean
signal, outside the address. The native pulse scalars cancel. This pays
the proposed source, not its coherent decision gap.

It has a common response on a large nonconstant-table family. Suppose
each pair of A rows is identical,
$`C_{2i,j}=C_{2i+1,j}=c_{ij}`$, and put
$`z_j=\bigoplus_i c_{ij}`$. Choose every A pair in its singlet and
leave the **entire B bank arbitrary**. At u=0 the completed pulse A0 is
identity on this subspace, since each paired X sum vanishes. The graph
word preserves the subspace and restricts to $`Z_B^z`$: each A pair has
odd computational parity. Every Bv pulse also preserves it. Hence
$`Q_C=I`$ there for every v, and the complete source restricts to

```math
 U_C\big|_{\rm singlets}
 =Z_B^zB_0^\dagger Z_B^zB_0
 =\exp\!\left(-\frac{i\pi}{4}Y_c
                    \sum_{j:z_j=1}X_{B_j}\right),
```

independently of v. These are permissible dirty test inputs, not
initialized source resources. The exact clean-zero compression on that
subspace is

```math
 M_C(0,v)\big|_{\rm singlets}
 =\begin{cases}
 0,&|z|\text{ odd},\\
 (-1)^{|z|/2}\displaystyle\prod_{j:z_j=1}X_{B_j},&|z|\text{ even}.
 \end{cases}
```

Whenever row zero of C is nonconstant, two addresses v require opposite
decision signs but have the same complete source action on this subspace.
No positive uniform coherent gap is possible. The family includes dense
tables at every larger size; it is not only a sparse embedded exception.
Signal-only processing of the source and its actual inverse retains the
common response. New address-dependent work operations, other independent
boundaries or a changed middle word require their own analysis. The
[interior-decision fixtures](../../tests/test_interior_decision.py) check
the literal asymmetric word and its inverse on the full small singlet
sector, with arbitrary signal and B inputs, including zero and nonzero
compressed responses. This failure supplies no unrestricted native-depth
lower bound.

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
called directly. The fixed-error scalar premise above, a zero-clean
full-space signed gap, or the one-clean decision source with fixed
positive gamma supplies the same logarithmic amplification factor.

<a id="square-root-width-signed-gap-transfer"></a>

### The first improvement needs only square-root source width

An all-width source is unnecessary for a first general improvement at
square-root dirty width. Suppose that for every prescribed a-bit table
a terminating classical compiler supplies either the zero-clean
full-space signed-gap source or the one-clean decision source above,
with a uniform constant $`\gamma\gt0`$ and simultaneous bounds

```math
 T_V,w_V=O(2^{a/2}),\qquad G_V=O(2^a),\qquad
 d_V=O(d(a)),
```

where $`d\ge1`$ is nondecreasing and the bounds include the actual
inverse. Then, for every fixed $`L\ge6`$, a sufficiently large fixed
$`K_L`$ gives a one-clean general frame compiler at
$`b=\lceil K_L\sqrt N\rceil`$, complete error $`2^{-L}`$, and

```math
 \begin{aligned}
 T&=O_L(\sqrt N),\qquad G=O_L(N),\\
 D_T&=O_L\!\left(\log\log(n+4)
       +H_n\log(H_n+2)d(n+C_L)\right),
 \qquad H_n=1+\log^*(n+2).
 \end{aligned}
```

Here $`C_L`$ is a fixed arity allowance. In the early stages of §1,
the source arity is at most $`n-(A-C')\chi(m)+O_L(1)`$ and
$`q_m=O_L(\log(h(m)+2))`$. Increasing A makes both
$`\sum_m q_m2^{a_m/2}=O_L(\sqrt N)`$ and
$`\sum_m q_m2^{a_m}=O_L(N)`$. The full-volume initial, terminal
and complex-phase calls have fixed L-dependent tolerances and arity at
most $`n+C_L`$, so there are only constantly many of them at constant
amplification overhead. The largest source fits a fixed multiple of
$`\sqrt N`$ dirty wires. The concurrently live routing and metadata
banks cost $`N^{o(1)}`$; choose $`K_L`$ to hold both and reuse the
source pool sequentially. Finite dimensions use the R96 fallback.
This proves the stated width and count bounds without postulating a
source at narrower allocations. It does not establish an all-width
interpolation theorem.

The precise sufficient depth target is

```math
 (1+\log^*(a+2))\log(3+\log^*(a+2))\,d(a)=o(\log a).
```

For example, $`d(a)=O(\log a/\log\log a)`$ already suffices;
a doubly logarithmic source is a stronger goal than necessary.
An unspecified $`d=o(\log a)`$ does not by itself control the
iterated-logarithm multiplier. Even a constant-depth source leaves
the current transfer's $`O_L(\log\log n)`$ routing/recursion term;
this first improvement would not settle optimal frame depth.

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

Generic coherent address motion need not satisfy the signed-gap
premise. For the vote example, tilt the $`|1\rangle`$ projector by
angle theta and control every vote XOR on that same tilted projector.
Each call has error $`2|\sin\theta|`$, but every odd majority returns
the same tilted phase. For QSP, spectral rounding similarly leaves
the rotated sign eigenspaces unchanged. The exact condition is the
prescribed sign of $`\mathcal H`$, rather than completed preservation
of every address by V. Motion within a correct sign sector, or motion
that cancels from the two-reflection signal, is allowed.

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

Alternatively, a zero-clean full-space signed-gap source for each
$`h_C`$ can first be amplified to a complete entry scalar at error
$`\tau/4`$, then used in the four-call reduction. Constant source
gap preserves the square-root-width corollary of §3 as well. This
argument does not assume that the product of four raw gapped signals
retains the required sign or a gap.

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
The QSP fixtures check the full dirty-space polynomial, literal endpoints
and the failure to correct a tilted sign basis. A nonflat gapped fixture
also moves a completed address across sign sectors while its two-reflection
operator has the correct sign and mixes addresses within each sector.
Its rotations are analytic matrices; the theorem separately pays native
synthesis and every actual inverse.
The asymptotic packing, density, rank and frame bounds above are analytic.
