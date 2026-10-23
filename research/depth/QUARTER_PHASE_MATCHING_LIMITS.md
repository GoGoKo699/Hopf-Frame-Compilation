# Quarter-phase matching stages: complete invariants and a spectral escape

[Resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Quarter-phase source](EXPANDER_ZERO_REFLECTION.md#quarter-phase-selector-stage) ·
[Route decisions](../ROUTE_DECISIONS.md)

The quarter-phase stage extends to a whole disjoint CNOT matching at
constant native T-depth. This pays simultaneous access to many address
coordinates. It does not by itself produce a normalized complete lookup
source. Two different obstructions below concern specified matching
grammars; neither is an unrestricted native T-depth lower bound.

Write $`m=2^p`$ and
$`\Omega|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle`$. The single-bank
Gram target is $`H_Y=(I+\Omega)/4`$. The fixed paired-bank coordinate
change in the [source interface](EXPANDER_ZERO_REFLECTION.md#quarter-phase-flip-difference-source)
transfers these statements to $`\Pi_{\rm eq}/2`$.

## 1. A complete matching stage with constant native depth

Let M contain q disjoint ordered pairs (c,t), with $`q\ge1`$.
All 2q address wires are distinct. Define

```math
r_{M,u}=\prod_{(c,t)\in M}
 \bigl(i^{x_c\oplus u_c}g_{(c,t),u}\bigr),\qquad
R_M=\prod_u r_{M,u}^{s_uY_u}.
```

Each control coordinate is unchanged by every factor. Thus both the
phase factors and anchored permutations commute across rows and pairs.
The address, bank Y and selectors s are arbitrary coherent inputs.
The matched permutation stage has the ledger (44) in the
[original proof](EXPANDER_ZERO_REFLECTION.md#independent-row-selector-echo).
For each pair and row, add its returned cubic quarter phase, with the
extra controlled-S and actual reversal for labels having $`u_c=1`$.
All physical T targets are private, so the qm row/pair phase words
share their Clifford controls and run in parallel. The same private
pool can be reused after the permutation stage has returned it.

The resulting simultaneous forward and actual-inverse bounds are

```math
\begin{aligned}
T_M&\le16m+16q+82qm,& d_M&\le124,\\
G_M&\le20m+20q+91qm,&
b_M&\le\max\{1+\max(m,q,qm/2),3qm\}.
\end{aligned}
```

There are zero clean inputs and zero implementation error. Every
additional private helper returns exactly on the complete space;
the m selectors are principal inputs and add m to source dirty width.
This is a composition of the already certified native primitives,
not a fresh emitted circuit for every qm. Two opposite perfect
matchings can therefore cover every address coordinate in constant
stage depth and polynomial intermediate resources. The next sections
show why coordinate coverage alone does not settle acceptance.

<a id="row-label-covariance"></a>

## 2. A row-label covariance screen

Consider any number of matching stages whose control coordinates lie
in a set S. They may have different targets and actual inverses.
Allow arbitrary address-only interlayers. Require all selector/signal
interlayers, clean insertions and acceptances to commute with swapping
two selector rows whose labels agree on S; all bank dependence occurs
through the displayed stages. In particular, uniform selector
Hadamards and global selector flips have this covariance.

If $`|S|\lt p`$, choose distinct labels u,v agreeing on S. The anchored
factor of every stage depends on its row label only through the
specified control coordinates. Swapping bank rows u,v and the
corresponding selectors therefore conjugates the whole source to
itself. Equivalently, its accepted blocks in the singleton banks
$`Y=e_u,e_v`$ are conjugate by that selector swap.

Use address u and any selector input invariant under the swap, with
the same other work. The two Gram expectations are equal, whereas
the desired expectations are respectively zero and one half. Hence

```math
\|K^\dagger K-H_Y\|\ge\frac14.
```

The analogous complete phase approximation error is at least one;
for target amplitude alpha it is at least $`|\alpha|`$. The proof
allows arbitrary stage count, clean count and address motion. It
does not apply to row-dependent selector operations that break the
covariance, nor to matching families covering all control coordinates.
The latter extension is paid in Section 1 and requires a different test.

<a id="matching-signature-singlet"></a>

## 3. Repeated local signatures leave an unmarked invariant

Fix a partition of $`p=2q`$ address bits into two-bit cells. On one
cell put

```math
r=S_0\mathrm{CX}(0;1),\qquad
h=S_1\mathrm{CX}(1;0),\qquad
D=S_0S_1=\mathrm{diag}(1,i,i,-1).
```

Both r and h fix $`|00\rangle`$ literally. On the nonzero states
1,2,3, their D-conjugates are real quarter rotations in two coordinate
planes. They generate the orientation-preserving signed-permutation
group of order 24. Thus the exact local group $`\mathcal G=\langle r,h\rangle`$
has 24 elements, and every $`L\in\mathcal G`$ fixes zero and has
a real orthogonal D-conjugate on the other three states.

The two-cell vector

```math
|\chi\rangle=
 \frac{-|1,1\rangle-|2,2\rangle+|3,3\rangle}{\sqrt3}
```

is fixed by $`L\otimes L`$ for every $`L\in\mathcal G`$.
Indeed the unnormalized maximally correlated vector is invariant
under $`O\otimes O`$ for real orthogonal O; applying $`D\otimes D`$
gives the displayed signs. This argument includes all local words
and their actual inverses, not only one r or h.

Now allow k anchored stage types. On cell j, stage a uses an arbitrary
local element $`L_{a,j}\in\mathcal G`$, conjugated by that cell's
two bits of the row label u. Stages may be controlled by the arbitrary
row selectors as above and used in any word, with arbitrary interleaved
operations on selectors, signals and other work. The latter operations
act as identity on the address. No other address operation or bank
dependence is included in this grammar. Count every distinct allowed
cellwise address layer in k.
These generalized stages are granted as logical operators for this
obstruction; only the actual matched form in Section 1 has its displayed
native ledger. The proof needs no reordering of different row factors.

Each cell has a signature
$`(L_{1,j},\ldots,L_{k,j})\in\mathcal G^k`$. If $`q>24^k`$,
two cells have the same signature. In the single-active-row bank
$`Y=e_u`$, translate chi on those cells by the corresponding bits of u
and put every other cell at its marked value. Both this unmarked
vector and the marked basis vector $`|u\rangle`$ are fixed by every
anchored stage, for every selector input. All inactive-row factors
are identity.

The complete source consequently acts identically on their work fibers.
Any common clean insertion and acceptance has equal Gram expectations
on the two vectors. Their target expectations differ by one half, so
the complete Gram error is at least one quarter. Thus error strictly
below one quarter in this grammar requires

```math
q\le24^k,\qquad k\ge\lceil\log_{24}q\rceil.
```

This is a signature-count obstruction for a fixed cell partition and
local group. It does not lower-bound unrestricted native T-depth;
native joint compilation is not identified with a count of logical
stages. Changed cell partitions, additional address interactions or
local elements outside this conjugate-real group can leave the premise.
Their complete accepted action still needs its own proof.

<a id="fixed-bipartition-invariant"></a>

## 4. Changing pairings inside a fixed bipartition is insufficient

A separate obstruction allows different cell partitions. Let
$`p=2q`$ with positive even q, and fix equal disjoint coordinate sets
A,B. Allow every perfect matching directed entirely from A to B or
entirely from B to A, together with their inverses. Put
$`z=i^q\in\{1,-1\}`$ and

```math
|\psi\rangle=
 \frac{|1_A0_B\rangle+|0_A1_B\rangle+z|1_A1_B\rangle}{\sqrt3}.
```

Every dressed A-to-B matching exchanges the first and third basis
vectors with phase z, and fixes the second. Every B-to-A matching
exchanges the second and third with phase z, and fixes the first.
Therefore all these operations and their inverses fix psi, independently
of the chosen pairing. They also fix the zero address.

On the single-active-row bank $`Y=e_u`$, translate psi by u. Arbitrary
selector/signal processing cannot distinguish its address fiber from
the marked fiber, by the same argument as Section 3. The complete
half-Gram error is at least one quarter for every word in this grammar,
regardless of the number of matchings. Permuting coordinates within
A and B also preserves the witness.

Changing the actual control partition can break this invariant. At
p=4 with $`A=\{0,2\}`$ and $`B=\{1,3\}`$, the unnormalized
vector $`|5\rangle+|10\rangle-|15\rangle`$ is sent by matching
$`(0;2),(1;3)`$ to
$`i|1\rangle+i|2\rangle+|3\rangle`$. The supports are disjoint.
This is an explicit escape from the invariant, not proof of a
complete lookup Gram or a uniform spectral gap.

For all preceding screens, approximating the ideal source unitary in
operator norm by epsilon changes its accepted Gram by at most
$`2\varepsilon`$. The respective one-quarter floors decrease by
at most that amount. Arbitrary reference inputs remain included.

<a id="dressed-expander-gap"></a>

## 5. A changed matching family has a uniform one-row spectral gap

Use the fourteen nonempty disjoint CNOT matching involutions
$`C_1,\ldots,C_{14}`$ on the unpadded admissible lengths $`p=3s`$
from [the supplied expander](EXPANDER_ZERO_REFLECTION.md#1-the-supplied-cnot-expander-and-its-representation).
Let $`\delta=1/(48\cdot10^5)`$ denote its positive constant gap on the
nonzero-address permutation representation, for every $`s\ge1`$.
This is the inherited expander premise;
no extension to other lengths is asserted here. Dress every matching by

```math
r_j=\left(\prod_{c\text{ a control of }C_j}S_c\right)C_j,
\qquad
H=I-\frac1k\sum_{j=1}^{k}\mathrm{Re}(r_j),\qquad k=14.
```

The controls and targets in each matching are disjoint, so
$`r_j^2=Z_{v_j}`$, where $`v_j\ne0`$ is its control mask.
Every r_j fixes zero. For any unit vector psi orthogonal to zero,
write a for its entrywise absolute values and u for the uniform unit
vector on nonzero addresses. With $`E=\langle\psi,H\psi\rangle`$,
the reverse triangle inequality and the undressed expander gap give

```math
E=\frac1{2k}\sum_j\|(I-r_j)\psi\|^2
 \ge\delta\bigl(1-\langle u,a\rangle^2\bigr),\qquad
\|a-u\|\le\sqrt{2E/\delta}.
```

For any j, let $`P_j=(I-Z_{v_j})/2`$. Factoring $`I-r_j^2`$
and using $`\|I+r_j\|\le2`$ gives
$`\|P_j\psi\|\le\|(I-r_j)\psi\|\le\sqrt{2kE}`$.
Moreover $`\|P_ja\|=\|P_j\psi\|`$ and
$`\|P_ju\|^2=2^{p-1}/(2^p-1)\ge1/2`$. Consequently

```math
\frac1{\sqrt2}\le\sqrt{2kE}+\sqrt{2E/\delta},\qquad
H\ge\gamma(I-|0\rangle\langle0|),\qquad
\gamma=\frac1{4(\sqrt{k}+1/\sqrt\delta)^2}.
```

This proves a uniform gap and removes every additional common fixed
vector for this changed family. There is also a contracting
non-Hermitian average

```math
A=\frac{2I+\sum_{j=1}^{14}r_j}{16},\qquad
\left\|A\big|_{\{|0\rangle\}^{\perp}}\right\|
 \le\sqrt{1-7\gamma/32}.
```

The convex-unitary variance identity, retaining just the identity-to-r_j
pairs, gives
$`1-\|A\psi\|^2\ge\sum_j\|(I-r_j)\psi\|^2/128=7E/32`$.
Conjugating every r_j by the same row translation proves the same
single-row statement at any marked address.

Individual matching stages have the complete native ledger of Section 1.
This spectral result supplies no new aggregate coherent SELECT or
accepted-source recovery ledger. In particular it does not establish
the arbitrary-bank source: products controlled by all row bits and
dirty selectors have a different Gram operator. The scalar zero
reflection is already paid in its original proof. The remaining
constructive task is to retain this changed address algebra while
deriving a normalized complete arbitrary-bank acceptance polynomial.

The [row-isolation continuation](BANK_BIT_RESPONSE_AND_ROW_ISOLATION.md)
now pays a fixed-row SELECT after cancelling each stage's background.
It does not supply coherent row selection or normalized full-bank acceptance.

The [exact fixtures](../../tests/test_quarter_phase_matching_limits.py)
check the local group, tensor invariant, heterogeneous signatures,
bipartition witnesses and an actual changed-partition action. They do
not certify the inherited asymptotic expander theorem or supply a
scalable native emitter. These are analytic proofs and internal finite
checks, not external mathematical review.
