# Masked sum compression with arbitrary dirty work

[Counter and query composition](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters) · [Amortized frame tradeoff](AMORTIZED_DIRTY_LOOKUP.md)

The dirty-counter indicator needs a reversible sum whose unknown helper
offset is independent of its input counters. Compressing columns into two
words supplies that interface. A schedule with stable output prefixes
improves the indicator T-depth to

```math
O(\chi(k)),\qquad
\chi(t)=\log_2(t+2),
```

with polynomial dirty width and gate count. All identities hold on the
full input space, with literal phase and arbitrary reference correlations.
There are no initialized helpers, measurements, resets, or catalysts.

## 1. The sum interface and imported primitives

Let $`k,m\ge1`$, $`M=2^m`$, and let $`a_1,\ldots,a_k,c`$
be arbitrary m-bit registers. We construct an exact reversible word

```math
L:(a,c,h)\longmapsto
\left(a,c+\sum_{i=1}^k a_i+g(h)\pmod M,h\right),
```

where h is arbitrary dirty work and g is a fixed function of its initial
bits only. The function g need not vanish. The construction has

```math
T=O\!\left(km+k(3/2)^m\right),\qquad D_T=O(m),
```

```math
G=O\!\left(km+k(3/2)^m
 [1+(\log(k+1)+m)\log(m+2)]\right).
```

Its full live width, including a and c, is
$`O(km+k(3/2)^m)`$. Every dirty wire returns at completion.
The depth is T-depth: arbitrary Clifford circuits between T layers have
their elementary gate count charged in G, without a total-depth claim.

The round-based baseline in Sections 2–3 imports the controlled increment of
[Vandaele, arXiv:2603.12917v1, Section 5, Theorem 4 and Corollary 7](https://arxiv.org/html/2603.12917v1#S5).
For two controls and an r-bit target, it uses $`O(r+1)`$
$`\{\mathrm{CCX},\mathrm{CX},X\}`$ gates, total depth
$`O(\log(r+2))`$, and one arbitrary returned dirty helper. Controls
may change internally but return exactly. Negative controls use X
conjugation. The completed corollary, rather than an internal promise
gate, supplies the all-input contract. Exact native Toffoli substitution
preserves these gate-count and depth orders.

The faster schedule proved in Section 4 instead needs only the exact
linear TTK modular adder and the
[TTK signed increment](PARALLEL_DIRTY_LOOKUP.md#a-read-only-controlled-increment-using-two-additions).
On an r-bit target they have $`T,G,D_T=O(r)`$; the adder needs no
helper, and the signed increment uses r arbitrary returned helper bits.
Both are already given as native words in the linked chapter. They may
borrow their private operands internally. Every subtraction or
uncomputation below uses the actual reversed native word. The optimized
increment and RV-adder depths are not needed for the new pipeline bound.

## 2. A local masked compressor

Take three active bits a,b,c of weight $`2^j`$, with
$`0\le j\lt m`$. For $`j\lt m-1`$, reserve a fresh arbitrary
dirty word D of length $`r=m-j-1`$, whose low bit has weight
$`2^{j+1}`$, and one separate returned dirty increment helper.
Execute CNOT from a to b, then from the updated b to c. The three
wires now hold

```math
a,\qquad p=a\oplus b,\qquad s=a\oplus b\oplus c.
```

The original integer carry has the disjoint-predicate expression

```math
\mu=\mathrm{Maj}(a,b,c)=a(1-p)+p(1-s),\qquad
a+b+c=s+2\mu.
```

Increment D once controlled on $`(a,\neg p)`$, and once controlled
on $`(p,\neg s)`$. The predicates cannot both hold. Each completed
increment preserves a,p,s and returns its helper, which is reused for
the second call. The result is $`D'=D+\mu\pmod{2^r}`$, so

```math
2^j(a+b+c)+2^{j+1}D
=2^j s+2^{j+1}D'\pmod M.
```

Keep s active at column j and every bit of D' active at its own higher
column. Retire a,p as inactive garbage; do not erase or reuse them.
The map is reversible: subtract the two predicate increments and undo
the two CNOTs. If $`j=m-1`$, D is empty: just perform the CNOTs and
discard the carry modulo M from the active-sum ledger. Nothing is
discarded from the physical circuit.

The compressor uses $`O(m-j)`$ gates and fresh wires and has depth
$`O(\log(m-j+1))`$. Its carry is encoded in a dirty upper word,
not in a designated bit promised to contain the unmasked majority.

## 3. The round-based baseline and weighted count

Initially column j contains the k bits at that position in the input
counters. In each round, partition every column's current active bits
into disjoint triples and apply the local compressor to all triples in
parallel. Reserve separate fresh D words and increment helpers for all
compressors, across every round. Newly produced bits enter the next
round. All gate supports within a round are disjoint, including controls;
the imported increment's temporary borrowing therefore causes no overlap.

Let $`h_j(t)`$ be the active height before round t and put
$`q_j(t)=\lfloor h_j(t)/3\rfloor`$. Every lower-column compressor
introduces one bit at column j, giving the exact recurrence

```math
h_j(t+1)=h_j(t)-2q_j(t)+\sum_{i\lt j}q_i(t),\qquad h_j(0)=k.
```

The schedule and all fresh reservations depend only on these integer
heights, never on quantum data. Put $`e_j=\max\{h_j-2,0\}`$.
Then $`q_j\le e_j`$ and
$`\max\{h_j-2q_j-2,0\}\le e_j/3`$, so

```math
e_j(t+1)\le\frac{e_j(t)}3+\sum_{i\lt j}q_i(t).
```

For $`\Phi(t)=\sum_{j=0}^{m-1}4^{-j}e_j(t)`$, summing the
geometric tails gives

```math
\Phi(t+1)\le\frac{\Phi(t)}3+
\sum_i q_i(t)\sum_{j\gt i}4^{-j}
\le\frac23\Phi(t).
```

Initially $`\Phi(0)\le4k/3`$; any nonzero value is at least
$`4^{-(m-1)}`$. Hence after $`O(\log(k+1)+m)`$ rounds every
column has at most two active bits. Each round has depth
$`O(\log(m+2))`$.

For the total count, let C_j be the number of compressors ever applied
at column j. Nonnegative final heights imply

```math
2C_j\le k+\sum_{i\lt j}C_i,\qquad
C_j\le\frac{k}{2}(3/2)^j.
```

Charging each compressor its actual upper-word length avoids an extra
factor m:

```math
\sum_{j=0}^{m-1}C_j(m-j)
\le\frac{k}{2}(3/2)^m\sum_{u=1}^m u(2/3)^u
=O\!\left(k(3/2)^m\right).
```

This pays for all permanently reserved fresh masks and increment helpers,
as well as their gates. The original counters remain separately charged.

## 4. Stable prefixes remove the round-depth factor

The new schedule processes the same masked carries, but does not wait for
a complete upper-word increment before using its finished low bits.
All controls shared between concurrent operations remain computationally
unchanged throughout those operations, not merely restored at their end.

### Static forests with deferred parity changes

For column j, its raw inputs are the k original counter bits and the
column-j bit of every fresh upper word allocated to a lower-column
compressor. Thus $`H_j=k+\sum_{i\lt j}C_i`$. Build a balanced ternary
forest on these wires, grouping triples at each level until at most two
survivors remain. Each internal node is a local compressor, so
$`2C_j\le H_j`$ and the preceding count bounds still apply.
All allocation and topology are determined before execution.

Defer every physical forest CNOT until all upper-word increments finish.
Record each node's a,p,s as explicit linear forms in its column's raw
wires. Its carry is still $`\mu=a(1-p)+p(1-s)`$. These forms can
overlap between nodes; no initialized copy is introduced. All raw wires
remain unchanged once they are available. A balanced forest has height
$`O(\log(H_j+2))`$, and the sum of the support sizes of its node
forms is $`O(H_j\log(H_j+2))`$: each leaf occurs only a constant
number of times per ancestor.

### Read-only affine controls with private phase work

For Boolean affine forms u,v and a private target t, implement their
Toffoli using H on t and a CCZ phase. One private dirty phase helper d
suffices for the exact identity

```math
\sum_{z\in\mathbb F_2^3}(-1)^{|z|}
 (d\oplus z_1u\oplus z_2v\oplus z_3t)
\equiv4uvt\pmod8.
```

For each of the eight terms, compute that parity into d, apply T or its
inverse with the displayed sign, and undo the parity word. Affine
constants use X on d. T gates act only on d; all raw-form wires occur
only as CNOT controls. The helper returns exactly and the phase is
literal, including for correlated or overlapping forms. Independent
targets with separate helpers therefore batch in eight T layers even
when their controls overlap. The Clifford count includes every parity
support; there is no constant-depth fanout assumption.

A conjunction of q stable forms can similarly toggle a private target
in $`O(q)`$ T-depth, count, and dirty width. For $`q\ge3`$, reserve
$`q-2`$ dirty chain bits. The full ladder toggles its first bit by the
first two controls, then each successive bit by the preceding chain bit
and the next control. Toggle the target using the last chain bit and
last control, then reverse the ladder. Repeat this word with the first
ladder update omitted. The unknown chain terms cancel, leaving exactly
the product of all q controls; all chain bits return. For q equal to
two use the preceding Toffoli directly. Apply the private-phase schedule
to every ladder gate. The original forms are never targets, including
when a literal is negative. Each form occurs only a constant number of
times in a completed conjunction word.

### A block increment with a permanently settled prefix

For a node's r-bit fresh upper word D, use consecutive blocks with
start and length $`(0,1),(1,1),(2,2),(4,4),\ldots`$, truncating
the last block. Before a block starting at ell, all lower bits have
already been incremented and will never change again. Its carry is

```math
f=\mu\prod_{i\lt\ell}(1-D_i).
```

Indeed, an enabled increment carries past the lower ell bits precisely
when their updated value is zero. This is also valid for ell equal to
zero, with the empty product one.

Allocate a private dirty bit g. Let Q increment just the current block
controlled on g, using the TTK signed increment and private helper bits.
Let E toggle g by f. It is the product of two conjunction words with
controls $`(a,\neg p,\neg D_0,\ldots,\neg D_{\ell-1})`$ and
$`(p,\neg s,\neg D_0,\ldots,\neg D_{\ell-1})`$; their enabling
predicates are disjoint. Let $`C_g`$ complement every current-block bit
controlled on g. The chronological word

```math
C_g,\ Q^\dagger,\ E,\ Q,\ E,\ C_g
```

increments this block by f and returns all helpers. To check it, the
middle four operations translate by $`f(1-2g)`$; complement has action
$`x\mapsto(1-2g)x-g`$, which conjugates that signed translation
to $`x\mapsto x+f`$. The enabling f depends on neither g nor the
current block, so this calculation holds on every input.

For block length u this costs $`O(\ell+u+1)`$ T gates, dirty wires,
and T-depth. Reserve separate helpers for each block and node. The sum
of ell and u over the doubling blocks is $`O(r)`$. More importantly,
each bit at offset v is permanently settled after at most
$`A(v+1)`$ T layers from its node's start, for one fixed constant A.
All later blocks read that prefix only through the phase-localized
conjunctions. Q and its inverse touch only their current block and
private helpers; they cannot borrow a settled bit.

### Column deadlines, exactness, and the refined ledger

Launch all nodes of column j at time Aj, measured in T layers. A raw
column-j bit originating at column i has offset $`v=j-i-1`$ in
its upper word. It settles by $`Ai+A(v+1)=Aj`$. Original counter
bits are available at time zero. Thus induction makes every launch
valid, and all upper-word increments finish by $`O(m)`$ T layers,
independently of k. Concurrent native gates have separate private
targets and work; any shared wire is a stable CNOT control.

Now emit all deferred forest CNOTs. This is the same exact compression
as the column-ordered sequence of local macros: within a column, the
recorded forms are its conceptual intermediate parities; higher-column
increments wait until every lower-column bit that they read is settled.
Postponing the parity changes prevents interference with still-running
prefix controls. No information is discarded. Call the resulting word U.

The T count and live width per node are $`O(r+1)`$, so the weighted
sum from Section 3 still bounds them by $`O(k(3/2)^m)`$.
Each node reads its forms $`O(\log(r+2))`$ times across blocks.
Since $`H_j\le k(3/2)^j`$, the additional Clifford count is at most

```math
O\!\left(\sum_j H_j\log(H_j+2)\log(m+2)\right)
=O\!\left(k(3/2)^m(\log(k+1)+m)\log(m+2)\right).
```

This accounts for shared parity fanout explicitly. The original counters,
final additions, and padding supply the remaining $`O(km)`$ terms.

## 5. Returning the work and cancelling its offset

Pad each final column to exactly two active bits with fresh arbitrary
dirty bits, at most 2m in total. Arrange the active bits into two m-bit
words $`z_0,z_1`$; choosing wire lists does not require moving data.
Let U be the complete compression word. Its local invariants give

```math
z_0+z_1=\sum_i a_i+g(h)\pmod M.
```

Explicitly, g sums each fresh upper word's initial value with its weight
$`2^{j+1}`$, plus the weighted initial padding bits. These words are
untouched before their assigned compressor; their initial offsets are
independent of a. Returned increment helpers have coefficient zero.
Retired garbage is retained physically throughout U.

Execute chronologically

```math
U,\quad\mathrm{ADD}_m(z_0;c),\quad
\mathrm{ADD}_m(z_1;c),\quad U^\dagger.
```

The two TTK adders preserve both words. Thus the actual inverse restores
all input counters, fresh masks, padding, and increment helpers, proving
L's stated action. The additions cost $`O(m)`$ gates and depth, so
L has T-depth $`O(m)`$. Its inverse has exactly the same resource
bounds and restores the work by reversing every native gate.

Now let J increment $`a_i`$ by a read-only literal $`\ell_i`$, using
the existing signed-increment circuit and a separate helper pool. J
returns that pool and does not touch h. The chronological echo

```math
L^\dagger,\ J,\ L,\ J^\dagger
```

adds $`\sum_i\ell_i\pmod M`$ into c: both the original counter
offsets and the identical $`g(h)`$ cancel. This identity needs neither
clean masks nor a separately computed carry bit. Literal permutations
and exact native decompositions extend it to arbitrary coherent inputs
and their reference systems.

## 6. Indicator and complete-frame consequences

Set $`m=\lceil\log_2(k+1)\rceil`$, so $`M\gt k`$, and retain
the [cyclic routing and output echoes](PARALLEL_DIRTY_LOOKUP.md#remove-the-remaining-counter-offset-by-cyclic-routing).
They extract the conjunction exactly and return every counter and mask.
The signed increments, cyclic router, and their work cost $`O(km)`$
gates and width and $`O(m)`$ T-depth. With $`\alpha=\log_2 3`$,
the complete row has T count and width $`O((k+1)^\alpha)`$,
Clifford count $`O((k+1)^\alpha\log(k+2)\log\log(k+4))`$,
and depth $`O(\chi(k))`$. For k equal to zero use X directly.

All $`H=2^k`$ equality rows have separate private work. Their only
shared-address interactions are Clifford CNOTs in J; all non-Clifford
layers have disjoint supports across rows. Consequently the exact dirty
indicator has

```math
T,w=O\!\left(H(k+1)^\alpha\right),\qquad
G=O\!\left(H(k+1)^\alpha\log(k+2)\log\log(k+4)\right),
\qquad D_T=O(\chi(k)).
```

No constant-depth Clifford fanout is assumed. The polynomial-overhead
interface in the [general hybrid proof](PARALLEL_DIRTY_LOOKUP.md#every-eligible-width-and-precision)
still admits $`P(t)=C(t+1)^3`$. Thus, for $`n\ge1`$,
$`N=2^n`$, every $`L\ge6`$, and $`b\ge17(L+n+7)`$, the same
two-clean prescribed complete real-frame circuit satisfies

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL+n\chi(n)\right).
```

The precision allocation, full-isometry error, literal width threshold,
and actual-inverse amplification are unchanged. Count and depth match
their worst-case lower bounds on the nonempty range
$`17(L+n+7)\le b\le\sqrt{NL/(nL+n\chi(n))}`$.
At fixed L and sufficient $`b=\Theta(\sqrt N)`$, this gives
optimal-order $`T=\Theta(\sqrt N)`$ and depth $`O(n\chi(n))`$.
It does not prove optimal large-width T-depth, improve the selected
high-precision endpoint, or extend the separate complex/state schedules.

## 7. Scope and evidence

The following interface obstructions explain the masks without giving
depth lower bounds. They allow arbitrary unitary action inside each
specified output-value sector, including changed helpers and garbage.

An exact unmasked parity/majority pair is already impossible on full
inputs: its weight-one input sector has dimension $`3\cdot2^q`$,
exceeding the required output-pair sector's $`2^{q+1}`$, even with
q dirty helpers. The masked modular output has a different contract.

**Canonical modular offset requirement.** Let $`M=2^{r+1}`$ and
suppose instead that the designated output word is $`s+2D'`$, with
r upper bits and $`s=a\oplus b\oplus c`$, and that it satisfies
$`s+2D'=a+b+c+g(h)\pmod M`$ on every input. Then the distribution
of g on uniform dirty inputs is uniform over all even residues modulo
M. Consequently $`g\bmod M`$ has exactly r bits of entropy and
requires at least r arbitrary input helper bits.

The complete output word is uniform on the maximally mixed state. With
$`z=\exp(2\pi i t/M)`$, its nontrivial Fourier coefficients vanish,
whereas the input coefficient is
$`((1+z)/2)^3\mathbb E[z^{g(h)}]`$. The first factor vanishes
only at $`t=M/2`$. Thus the offset's Fourier coefficients vanish
except at zero and $`M/2`$. The parity identity forces g even, so
both remaining coefficients equal one. Fourier inversion gives
$`\Pr[g(h)=2u]=2^{-r}`$ for every u. A function of q uniform bits
has entropy at most q, proving $`q\ge r`$.

The fresh r-bit upper word saturates this offset requirement only;
additional implementation helpers remain charged. This is not a lower
bound on circuit depth or on noncanonical modular representations.

There is also a restriction on a single classical additive-potential echo.
If $`L_h:(u,c)\mapsto(u,c+h(u)\bmod M)`$ and
$`J:(u,c)\mapsto(j(u),c)`$, the chronological word
$`L_h^\dagger,J,L_h,J^\dagger`$ adds $`h(j(u))-h(u)`$ to c.
If this equals a constant delta for every u, telescoping along any
length-r cycle of j gives $`r\delta=0\pmod M`$. For delta equal
to one every cycle has length divisible by M; at least
$`\lceil\log_2 M\rceil`$ modified bits are required, and an
involution cannot give a unit increment when $`M\gt2`$.
This concerns that single echo form only. The existing signed increment
uses an offset-dependent inner shift and Clifford sign correction, so
it is not excluded. Neither observation is an unrestricted depth lower
bound or an obstruction to the masked construction above.

The [baseline interface checks](../tests/test_dirty_sum_interfaces.py)
audit the local macro, round schedule, helper offsets, complete sum echo,
and actual inverses. The [pipeline checks](../tests/test_pipelined_dirty_sum.py)
separately audit deferred parity forests, signed block increments,
updated-prefix carries, full sum echoes, and static readiness deadlines.
The new pipeline uses the explicit TTK arithmetic and phase words above;
it does not infer a fast depth from a slower arithmetic fixture. Finite
checks support these interfaces; the uniform schedule and Hopf resource
bounds follow from the proofs above.
