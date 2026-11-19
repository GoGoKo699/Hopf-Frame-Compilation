# Prescribed scalar signs as a sufficient source for general frame depth

A scalar sign on $`a`$ address qubits is the literal operator
$`D_g|x\rangle=(-1)^{g(x)}|x\rangle`$, tensored with identity on all
arbitrary dirty work. This chapter gives a conditional unrestricted
compiler, an average-case criterion for its source, and an explicit
affine matrix pencil encoding each source table. All three use the
complete frame convention, including arbitrary references and actual
native inverses. The [unrestricted upper theorem R96](ALL_PRECISION_LINEAR_DEPTH.md)
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
boundary on all physical clean/dirty sectors. After bounded-clean
promotion suppose the fixed full error is at most
$`\epsilon_0\le10^{-3}`$. Form the XOR of $`g(x)`$ by Hadamards
on the extra input of $`t g(x)`$. Record $`k`$ votes into fresh
zero bits with the same dirty core. At fixed $`x`$, the wrong-record
Kraus operator at each step has norm at most $`\epsilon_0`$ and later
calls leave earlier votes alone. Orthogonality of records bounds the
squared mass of wrong majorities by $`2^k\epsilon_0^k`$ for odd
$`k`$. Majority phase and actual uncomputation give clean-code error
$`2(2\epsilon_0)^{k/2}`$. Apply this to $`u g(x)`$, conjugate by
Hadamards on $`u`$, promote the $`k`$ votes at error factor
$`8\cdot2^k`$, and extract a scalar by the Z echo. The resulting
full error is at most $`32(\sqrt{8\epsilon_0})^k`$. The promotion
and the exact majority cost $`\mathrm{poly}(1/\tau)`$ for requested
error $`\tau`$, with at most two extra scalar address bits. [Promotion,
§§3–4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)

Request $`\tau_m=c2^{-L}/(1+h(m))^2`$. The reserve absorbs each
polynomial amplification cost in the volume terms, and the achieved
frame depth is

```math
 O_L\!\left(N/b^2+\log\log(n+4)
 +(1+\log^*(n+2))^{C_1}d(n+O_L(1))\right)
```

for a fixed constant $`C_1`$. The other resource bounds from §1 hold.
A source with an explicit precision-uniform ledger of the required
polynomial form can also be called directly. To see why generic
constant coherent leakage is insufficient for this particular vote,
tilt the $`|1\rangle`$ projector by angle $`\theta`$ and control
every vote XOR on that same tilted projector. Each call has error
$`2|\sin\theta|`$. The only records are all-zero and all-one on the
two tilted sectors; every odd majority returns the same tilted phase.

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

The [finite checks](../../tests/test_prescribed_scalar_tables.py) recover
every two-bit target from shifted even-parity tables, verify the
derivative vote as a coherent unitary, and test selected Walsh characters
and the coranks and marked-edge power on small tables.
The asymptotic packing, density and frame bounds above are analytic.
