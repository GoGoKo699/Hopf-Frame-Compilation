# Returned mask and moving-address zero-test reductions

[Complete indicator fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Bounded-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)

Put $`m=2^k`$. The target adds $`e_x`$ to an arbitrary $`m`$-bit bank,
returns the address and every borrowed wire, and has full operator error at
most the requested fixed $`\tau`$, including external references. Each
theorem below assumes a **physical source** with a terminating compilation
and its actual inverse. Let $`d^*(k)`$ be a common nondecreasing depth
bound for both at the fixed tolerances used. None of these reductions
constructs a sublogarithmic source.

## 1. A shifted odd mask on the same hidden sector

Let $`P_\sigma`$ be orthogonal helper projectors summing to identity and
$`h_\sigma`$ a binary word of odd weight for each sector. A callable on any
coherent address $`q`$ must add $`f_q(i)=h_\sigma(i+q)`$ into an arbitrary
bank on its whole sector cylinder, preserving the helper. All calls reuse
the **same** helper. Its basis and the masks can be unknown.

Borrow banks $`A,B`$ and address $`r`$. Write $`C`$ for binary group
convolution added into the problem bank, $`J`$ for the mask addition at $`r`$
into $`A`$ and $`K`$ for the mask addition at $`r+x`$ into $`B`$. CNOTs add
and remove the address offset. In chronological order apply
$`C,J,C,K,C,J^\dagger,C,K^\dagger`$. In each sector the four convolutions
add the mixed difference $`f_r\star f_{r+x}`$ and all borrowed inputs
return. In the binary group algebra $`h\star h=(\sum_i h_i)e_0`$:
distinct-index pairs cancel; every equal pair lands at zero. Translation
therefore gives exactly $`e_x`$. The proof extends over superpositions
of sectors. Supply the convolution using the
[ranked product and zeta changes](ALGEBRA_QUERY_EXTRACTION.md).
Four copies incur at most 160 additional T layers in the original
implementation, and $`O(mk^2)`$ additional T gates, Cliffords and dirty
width. Four full-operator source errors of size $`\delta`$ add to at most
$`4\delta`$. The counts also include four source visits and the retained
helper. Replicating unrelated sector variables would replace the square
by a product of two different masks.

A raw callable may move the address and helper and create quantum garbage
if it gives the specified *definite* updated bank on the whole cylinder.
Copy the old bank to a separate arbitrary target, run the raw callable,
copy the new bank and actually reverse. The deterministic second control
leaves its entire quantum output intact. This gives the strict
translation above with eight total raw visits in the mixed difference.
At $`m\,\mathrm{poly}(k)`$ T gates and dirty width and
$`m^2\,\mathrm{poly}(k)`$ elementary Cliffords, the linked indicator fusion
removes the polynomial factors with an additive doubly logarithmic depth.

## 2. Definite displacements on a positive fraction of inputs

A zero-clean unitary has designated output bank $`P(x,z)`$ definite on
every computational input; $`z`$ denotes all the dirty input, including the
old bank $`B(z)`$. The remaining outputs may be quantum. Suppose a known
binary linear map $`L_0`$ satisfies

```math
\Pr_{x,z}[F(x,z)=e_x]\ge\rho,\qquad
F(x,z)=L_0(P(x,z)+B(z)),
```

for uniform independent classical $`x,z`$ and fixed $`\rho\gt0`$. Assume
the relaxed resources at the end of §1, counting physical bank width.
Then for every fixed $`\tau`$ there is a full returned indicator with those
relaxed resources and depth $`O_{\rho,\tau}(d^*(k)+(\log\log(k+4))^C)`$
for an absolute $`C`$. A small full-operator error about such an ideal
unitary can be added per actual visit. Average success probabilities of
arbitrary quantum output banks do not imply the definite-bank promise.

**Classical predictions.** The same old-bank/new-bank copying echo,
with $`L_0`$ in the CNOTs, adds $`F`$ to a dirty carrier. Choose a fixed list
of $`O_\rho(k+1)`$ XOR shifts such that at every $`x`$ the mean success
over the shifted addresses and dirty words is at least $`\rho/2`$.
Independent sampling, concentration and a union bound over addresses
prove existence. Draw a shift $`s`$ from that list and add a uniform
offset to $`z`$. Validate the resulting proposal $`v`$ for $`u=x+s`$ by
requiring odd parity and testing whether

```math
\sum_i a_i v_i \langle q,i+u\rangle=0
```

for independent uniform $`a,q`$ of lengths $`m,k`$. A correct unit vector
always passes. An incorrect odd vector has a nonzero coordinate away
from $`u`$, so this test rejects with probability at least $`1/4`$. With
$`t`$ trials and $`\ell`$ validators per trial the chance of missing
all correct proposals or accepting a wrong one is at most
$`(1-\rho/2)^t+t(3/4)^\ell`$. This is one scalar bad event for *all*
output Fourier words $`y`$: return the parity of the first passing
proposal against $`y`$ shifted by its $`s`$. Both parameters are constants
at fixed error. All primitive parities are obtained by sign echoes on
the dirty carrier. Fixed Boolean combinations of returned bits use
borrowed mixed differences. Dense $`L_0`$ costs at most the allowed
quadratic Clifford count. Hardwired shift selectors of logarithmic
length use read-only conjunctions of doubly logarithmic T-depth.

**Coset smoothing.** Let $`W=O_\rho(m\,\mathrm{poly}(k))`$ count all
prediction coins. Dirty origins translate their coordinates. A fixed
matrix with $`W`$ rows and $`h=O((k+1)\log(W+2))`$ columns suffices:
every translate of its image has small bad fraction at every address.
For a bad set of density $`\delta`$, successively average its indicator
$`f_j`$ with its XOR shift in an independent random direction. The
exponential potential $`Z_j=\mathbb E_b\exp(\lambda 2^j f_j(b))`$ obeys
$`\mathbb E[Z_{j+1}\mid f_j]=Z_j^2`$. With probability at least $`3/4`$,
all steps obey $`Z_{j+1}\le8\cdot2^j Z_j^2`$ by Markov and summation.
Thus $`2^{-J}\ln Z_J\le\ln Z_0+\ln16`$. The maximum exceeds the log-mean
bound by at most $`W\ln2/(\lambda2^J)`$ in units of $`f_J`$. Choose
$`\lambda`$ large, $`\delta\le e^{-\lambda}`$, and $`2^J`$ a sufficiently
large constant times $`W`$. Independently concatenate $`k+1`$ groups of $`J`$
directions. A union bound gives a good group for each $`x`$. Averaging
over the sum of the groups preserves its bound on *every* coset.
The final short labels can be added through literal binary shears.
The scalar sampling argument in §4 with two ordinary-table stages gives
the stated fixed power; a sharper joint-density estimate would need its
shared bad-event and resource ledger checked at the proposed cost.

For a returned translation by an address-only vector of binary row rank
at least $`\rho m`$, pick independent pivot rows and map them linearly
to their labeled unit vectors. This gives the promise above. In this
case a direct short-list version gives $`O_{\rho,\tau}(d^*(k)+\log\log(k+4))`$.
Draw a separate list of trials at each of $`D=\mathrm{poly}(k)`$ labels
in each of the three mutually unbiased Pauli product bases. At each
address the bad-row subsets are independent Bernoulli restrictions in
different bases. For their flat transition matrices the random
restriction estimate gives $`\|P_iP_j\|\le\mu`$ simultaneously at all
addresses: choose moment order $`p=k+10`$, $`D`$ of order $`(4p)^4`$, and
failure density small enough that $`32\delta^{1/4}\lt\mu`$. The
tail at each of the $`3m`$ address/pair choices is at most $`2^{-2p}`$.
Finite search selects lists offline. Record the three votes successively
and use their majority. A history with two adjacent wrong votes contains
a product of two bad projectors; the separated wrong votes contain
$`P_3(I-P_2)P_1`$ of norm at most $`2\mu`$. Summing squared norms of wrong
histories costs at most $`7\mu^2`$. The majority phase and actual reverse
therefore cost at most $`2\sqrt7\mu`$ error before three-bit promotion.
Trial formulas use borrowed differences, so the number of initialized
votes is independent of $`\mu`$. Each label selects a fixed list, encoded
by polynomially many terms with logarithmically many controls. Dense
linear maps are paid as Clifford shears. [Tropp, Lemma 18 and Remark 10](https://tropp.caltech.edu/papers/Tro08-Linear-Independence-preprint.pdf)

## 3. A raw target-bit zero test with arbitrary address motion

The raw unitary $`R`$ uses zero clean wires. For both input target
values $`b`$, every nonzero address and any private state, the output
target is exactly $`b`$. On the zero address the operator amplitude of
the wrong target $`b`$ is at most a sufficiently small absolute
$`\delta`$ on all private inputs. Other outputs can be arbitrary.
Suppose the width and elementary gate counts fit $`2^{O(k)}`$ with an
arity exponent independent of subsequent fixed precision requests.

Copy the input bit to a second dirty target, apply $`R`$, copy the
output and actually reverse. The plus Fourier block on the new target
is identity. Its minus block is $`R^\dagger Z R Z`$. Off zero it is
identity on the whole input space, hence leaves the zero-address
subspace invariant. On zero, the same-bit error blocks have
orthogonal input and output target labels, so
$`\|(RZ+ZR)P_0\|\le2\delta`$. Thus the completed toggle and its
adjoint are within $`2\delta`$ of the desired flip there.

For candidate $`a`$, borrow $`m/4`$ disjoint $`k`$-bit origins $`w_t`$ and
private source lanes. In a prediction labeled by two field elements
$`A,B`$, take $`v_t=w_t+At+B`$ in $`\mathbb F_{2^k}`$. Add and remove all
affine offsets by CNOT, including the genuine logical $`x`$. In each lane
compute the product of the zero bits at $`v_t`$ and $`v_t+x+a`$ on two
different dirty carriers. Visit a Toffoli at the four successive
corners generated by the first toggle, the second toggle, and their
actual inverses. Exact displacements add their product with both
carriers returned. If either operation is identity the remaining
corners cancel exactly even for imperfect quantum behavior in the
other. For $`a\ne x`$ their supports are disjoint, giving exact silence
on the whole space. For $`a=x`$ the ideal parity is the parity of the
number $`H`$ of zero origins. In a fixed origin sector the perturbation
is at most $`C\delta H`$, including a parity echo to collect lanes.
Every completed single-address toggle preserves its address sectors.

For uniform $`A,B`$ the shifted origins are pairwise independent at
each fixed dirty column, so $`\mathbb EH=1/4`$ and
$`\mathbb E[H(H-1)]\le1/16`$. Hence $`\Pr[H=1]\ge3/16`$. Use a fixed
number of independent groups and OR their parities by borrowed
Boolean formulas. At the true address the miss probability decays
geometrically; false addresses remain exactly silent. Markov also
bounds total hits by a fixed cutoff outside a small constant event.
On the resulting common good event the Fourier forecast is uniformly
accurate for every remaining private quantum state. Choose the raw
$`\delta`$ after fixing the cutoff and the fixed formula lengths.
The sampler below supplies full, zero-clean error $`\tau`$ with
depth $`O_\tau(d^*(k)+\log\log(k+4))`$ in singly exponential volume.

## 4. Sampling with arbitrary dirty origins and normalizing size

The simulation requires an operator inequality for a **scalar bad
event** on the initial dirty coordinates. It never supplies independent
pure random states. Here are the ingredients for a label of
$`h=\mathrm{poly}(k)`$ bits.

* On $`r=O(h+1)`$ dirty sampler bits, coherently select a fixed large
  constant number of mutually unbiased Clifford bases. Tensor products
  of a fixed-size block suffice. For $`K=2^{\lfloor r/3\rfloor}`$,
  the projectors on the $`K`$ most probable outcomes in distinct bases
  have pair product norm at most $`K2^{-r/2}`$. The norm of their sum
  is at most one plus the sum of these pair bounds in a row of the
  block Gram matrix. Thus the *mean* mass above the $`1/K`$ cap is small.
  Capping and redistributing it yields nearby high-min-entropy laws.
* A strong linear extractor at fixed error has seed length
  $`u=O(\log^3(h+2))`$ and gives the requested output after increasing
  $`r`$ by a constant factor. Its selected binary matrix acts by a shear
  into an arbitrary label word. Toggle a dirty program bank by the
  seed indicator before and after a bilinear phase echo. The literal
  shear costs $`2^{O(u)}\,\mathrm{poly}(h)`$ volume and
  $`O(\log(u+2))`$ depth. An initial output mask translates a uniform law.
* Supply the $`u`$-bit seed with successive short ordinary table
  extractors and the same constant-basis preparation. At output
  length $`s`$, $`O(s+1)`$ independent random tables on $`O(s+1)`$ inputs
  suffice. For a fixed flat $`K`$-source and output subset the
  concentration exponent is proportional to $`K`$ times the number of
  tables. With constant slack it beats the logarithm of the count
  of source/subset pairs. An exhaustive offline search thus gives
  $`O(\log(s+2))`$ next seed bits. Each selected table XOR costs
  $`2^{O(s+1)}`$ volume and $`O(\log(s+2))`$ depth by the dirty indicator.
  With three tables the last seed occupies
  $`O(\log\log\log(u+16))`$ bits. A constant family of generic
  metric-uncertainty bases on that many bits plus a constant
  marginal suffices. Their brute-force native controlled synthesis
  at fixed error has depth $`2^{O(\log\log\log(u+16))}`$, within
  $`O(\log\log(k+4))`$. Two tables instead yield the fixed-power bound
  in §2. The discarded marginal is mathematical; every physical basis
  change is reversed. [Strong linear extractor, Lemma 4](https://drops.dagstuhl.de/storage/00lipics/lipics-vol124-itcs2019/LIPIcs.ITCS.2019.53/LIPIcs.ITCS.2019.53.pdf) · [Metric uncertainty, Theorem 2.5](https://arxiv.org/abs/1010.3007)

Condition on the outer coordinates and arbitrary original seed masks.
The inequality for each scalar event on the innermost sampler is
an operator inequality tensored with identity on other data. A direct
sum over the outer outcomes and conjugation by the chosen bases
preserve it. Repeating through the levels bounds the bad probability
on *any* input, including entangled dirty work. A Hermitian vote is
obtained by recording a prediction, phasing it and actually reversing.
Its compression to the initialized basis selectors has a fixed
sign gap. A fixed-gap phase-estimation filter and its actual reversal
give any fixed final accuracy. Increasing its clock length yields
exponentially decreasing filter error with logarithmically increasing
clock workspace, sufficient for the clean-count factor in
[clean-to-dirty XOR promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md).
The basis-selector count is constant. Comparisons at fixed precision
charge every call. Generic average quantum overlaps would not furnish
the scalar event used in this argument.

The direct zero-test construction can use $`2^{ck}`$ T gates, dirty
wires and elementary Cliffords for a uniform $`c`$. To obtain a
$`p`$-bit indicator, split $`p`$ into a fixed number $`s\gt c`$ of near-equal
blocks. On Fourier word $`y`$ phase by the sum of
$`y_j\prod_i B_{i,j_i}`$ over all full addresses $`j`$, where $`B_i`$ are
arbitrary short banks. Take the $`s`$-fold mixed difference with respect
to the translations $`B_i\gets B_i+e_{x_i}`$. Only $`y_x`$ survives.
Each corner can be a separate conjugation using the actual inverses,
so at most $`O_s(1)`$ small-source errors accumulate. The $`2^p`$
fixed-degree phases have private dirty implementations and constant
T-depth for fixed $`s`$; Clifford parities account for shared controls.
Choose the requested small-source error after fixing $`s`$. This gives
$`O_\tau(2^p)`$ T gates and dirty width, $`O_\tau(4^p)`$ elementary
Cliffords, and retains the depth order. The exponent uniformity
prevents a circular choice. This argument admits singly exponential
discovery volume in contrast to polynomial address overhead in §1.

## 5. Sparse support and the full-frame consequence

Let $`S`$ be known and nonempty. The raw bit is exactly silent off $`S`$.
On a known $`G\subseteq S`$ with $`|G|\ge\rho|S|`$ assume the *joint*
bound $`\|(ZR+RZ)P_G\|\le2\delta`$. For a fixed $`\alpha\gt0`$ assume
$`|S+S|\le m^{1-\alpha}`$. No constraint is placed on $`S\setminus G`$.
The completed echo is exactly identity off $`S`$ and preserves this
support, even when it can move an address within it. Use fixed
$`d\gt1/\alpha`$ nonzero field scalings $`\lambda_j`$ whose inverse
dilates of $`S+S`$ intersect only at zero. Independent uniform
nonzero scalings and a union bound over nonzero differences prove
existence at sufficiently large $`k`$; bounded exceptions can be compiled
directly.

For each dilation, use $`M=\max(1,\lfloor m/(4|S|)\rfloor)`$ disjoint
origins. At false addresses some dilation makes the two supports
disjoint, so its four-corner word is exactly identity, including on
quantum bad inputs. At the true address, pairwise independence bounds
the chance of a unique $`G`$ hit and no other $`S`$ hit below by

```math
\frac{M|G|}{m}\left(1-\frac{(M-1)|S|}{m}\right).
```

This is bounded away from zero. On this event the group incurs
$`O(\delta)`$ full error. Sectors with different occupied lane indices
are orthogonal and invariant, because each echo preserves membership
in $`S`$. Use independent affine seed pairs per dilation. Combine groups
by a balanced commutator tree: for distinct indices the transvections
$`I+pE_{ij}`$ and $`I+qE_{jl}`$ have commutator $`I+pqE_{il}`$. Implement
the local transvection by conditioning each four-corner parity phase on
the public source bit in the Fourier basis of its target. Any identity
child cancels the entire ancestor word. Thus every false address is
exactly silent. On the good intersection at the true address all child
transvections approximate their intended actions with $`O_d(\delta)`$
total error. Conjugate a Pauli X on the public source bit and cancel that
X to recover a returned flip.

Record a Hermitian Fourier vote $`V`$. If $`y_x=0`$ its result is positive
on every state. If $`y_x=1`$, for the favorable diagonal projector $`P`$
on initial origins, $`\|(V+I)P\|\le\epsilon`$ implies
$`V\le I-2P+3\epsilon I`$ by its two-by-two block decomposition. The
scalar sampler bounds the mean of $`P`$ in operator order even when
unfavorable branches scramble their work. The gap distinguishes the
two labels using a threshold just below one, fixed by $`\rho,\alpha`$;
choose $`\delta`$ accordingly. Two short tables, fixed-gap filtering,
bounded-clean promotion and the Cartesian normalization prove the
same counts as §4 with conservative added depth
$`(\log\log(k+4))^{O(1)}`$. A per-row error estimate on $`G`$ alone can
fail when wrong output components from several addresses interfere.

For all four sources, if $`d^*(k)=o(\log(k+2))`$ at every required
fixed accuracy, their normalized indicators and the
[approximate word-query transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md)
give $`o(\log n)`$ native depth at fixed precision $`L`$ and square-root
dirty width. At general admissible $`b`$ the depth is
$`O_L(N/b^2+d^*(2n+10)+a(2n+10))`$ with $`N=2^n`$, where $`a(k)`$ is
doubly logarithmic in §§1 and 3 and for the pivot-rank specialization.
Use the displayed fixed-power correction for joint density and sparse
support. The simultaneous counts are $`O_L(\sqrt N+N/b)`$ T gates and
$`O_L(N)`$ elementary Cliffords, with one external clean qubit and
the retained real/complex dirty thresholds. The native cost of a
qualifying raw source remains open. The finite tests exercise the
cancellation, affine hit probabilities and exact sparse support used above.
[Finite checks](../../tests/test_returned_mask_and_zero_test.py)
