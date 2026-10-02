# Masked sum compression with arbitrary dirty work

[Counter and query composition](PARALLEL_DIRTY_LOOKUP.md#6-a-polylogarithmic-depth-indicator-using-dirty-counters) · [Amortized frame tradeoff](AMORTIZED_DIRTY_LOOKUP.md)

The dirty-counter indicator needs a reversible sum whose unknown helper
offset is independent of its input counters. Compressing columns into two
words supplies that interface and improves the indicator T-depth to

```math
O(\chi(k)),\qquad
\chi(t)=\log_2(t+2)\,\log_2\log_2(t+4),
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
T,G=O\!\left(km+k(3/2)^m+m\log(m+2)\right),
\qquad D_T=O\!\left((\log(k+1)+m)\log(m+2)\right).
```

Its full live width, including a and c, is
$`O(km+k(3/2)^m+m)`$. Every dirty wire returns at completion.

We import the exact controlled increment of
[Vandaele, arXiv:2603.12917v1, Section 5, Theorem 4 and Corollary 7](https://arxiv.org/html/2603.12917v1#S5).
For two controls and an r-bit target, it uses $`O(r+1)`$
$`\{\mathrm{CCX},\mathrm{CX},X\}`$ gates, total depth
$`O(\log(r+2))`$, and one arbitrary returned dirty helper. Controls
may change internally but return exactly. Negative controls use X
conjugation. The completed corollary, rather than an internal promise
gate, supplies the all-input contract. Exact native Toffoli substitution
preserves these gate-count and depth orders.

For the final two modular additions, use the shortened-register
Remaud–Vandaele construction proved in the
[existing adder interface](PARALLEL_DIRTY_LOOKUP.md#exact-modular-adders).
It has no additional helper, preserves its first arbitrary operand,
and costs $`T,G=O(m\log(m+2))`$ and
$`D_T=O(\log^2(m+2))`$. Both operands are private; the source may
borrow them internally. Every subtraction or uncomputation below uses
the actual reversed native word.

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

## 3. A static parallel column network

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
h_j(t+1)=h_j(t)-2q_j(t)+\sum_{i<j}q_i(t),\qquad h_j(0)=k.
```

The schedule and all fresh reservations depend only on these integer
heights, never on quantum data. Put $`e_j=\max\{h_j-2,0\}`$.
Then $`q_j\le e_j`$ and
$`\max\{h_j-2q_j-2,0\}\le e_j/3`$, so

```math
e_j(t+1)\le\frac{e_j(t)}3+\sum_{i<j}q_i(t).
```

For $`\Phi(t)=\sum_{j=0}^{m-1}4^{-j}e_j(t)`$, summing the
geometric tails gives

```math
\Phi(t+1)\le\frac{\Phi(t)}3+
\sum_i q_i(t)\sum_{j>i}4^{-j}
\le\frac23\Phi(t).
```

Initially $`\Phi(0)\le4k/3`$; any nonzero value is at least
$`4^{-(m-1)}`$. Hence after $`O(\log(k+1)+m)`$ rounds every
column has at most two active bits. Each round has depth
$`O(\log(m+2))`$.

For the total count, let C_j be the number of compressors ever applied
at column j. Nonnegative final heights imply

```math
2C_j\le k+\sum_{i<j}C_i,\qquad
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

## 4. Returning the work and cancelling its offset

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

The adders preserve both words. Thus the actual inverse restores all
input counters, fresh masks, padding, and increment helpers, proving L's
stated action. Two final additions supply the $`m\log(m+2)`$ count
term; their $`O(\log^2(m+2))`$ depth is absorbed by the displayed
bound. Adding the original counters gives the complete width ledger.

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

## 5. Indicator and complete-frame consequences

Set $`m=\lceil\log_2(k+1)\rceil`$, so $`M\gt k`$, and retain
the [cyclic routing and output echoes](PARALLEL_DIRTY_LOOKUP.md#remove-the-remaining-counter-offset-by-cyclic-routing).
They extract the conjunction exactly and return every counter and mask.
The signed increments, cyclic router, and their work cost $`O(km)`$
gates and width and $`O(m)`$ T-depth. With $`\alpha=\log_2 3`$,
the complete row therefore has gates and width $`O((k+1)^\alpha)`$
and depth $`O(\chi(k))`$. For k equal to zero use X directly.

All $`H=2^k`$ equality rows have separate private work. Their only
shared-address interactions are Clifford CNOTs in J; all non-Clifford
layers have disjoint supports across rows. Consequently the exact dirty
indicator has

```math
T,G,w=O\!\left(H(k+1)^\alpha\right),\qquad D_T=O(\chi(k)).
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

## 6. Scope and evidence

Two narrow obstructions explain the use of masks and signed increments.
An exact unitary cannot send every three-bit input to designated bits
$`s=a\oplus b\oplus c`$, $`t=\mathrm{Maj}(a,b,c)`$ while allowing
only arbitrary dirty extra work. With q such helper bits, the three
weight-one inputs span an input sector of dimension $`3\cdot2^q`$,
but all require the output pair $`(s,t)=(1,0)`$, whose whole output
sector has dimension $`2^{q+1}`$. Modifying helpers or allowing
superpositions within that sector cannot make an isometry. The present
compressor instead outputs a masked upper word and retains a,p.

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

The [bounded interface checks](../tests/test_dirty_sum_interfaces.py)
audit the local macro, column schedule, helper offsets, complete sum
echo, and actual inverses. They use serial controlled-increment words;
Vandaele's logarithmic-depth synthesis and the final RV adder depth are
analytic imports, not depth inferred from those fixtures. The uniform
network bounds and Hopf composition are the proofs above.
