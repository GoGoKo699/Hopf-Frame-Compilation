# Table-specific probes and the remaining joint-word freedom

[Universal source gauges](WEAKER_INDICATOR_INTERFACES.md) ·
[Complete bilinear queries R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Restored-address frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md)

A fixed table can permit a source to change the dirty indicator bank,
as well as its address and private work. The exact freedom is determined
by a binary syndrome. This chapter identifies that freedom, gives its
paid equivalence to a smaller complete row query, and prices a safe
two-source composition. For full row rank the freedom collapses to the
universal indicator interface. Low rank is not guaranteed by the general
frame grammar. No faster source or complete-frame improvement follows;
R96 remains the strongest unconditional general upper bound.

This is a restriction on the specified universal bilinear probe at an
individual source boundary. It is not a necessity theorem for a whole
table query. In particular, a joint word can mix both address fields and
both borrowed banks without supplying either individual probe.

## 1. The table sees a syndrome, not necessarily the whole bank

Let the left address x have p bits, put $`H=2^p`$, and let A be an
arbitrary H-bit bank. Fix a binary matrix $`F\in\mathbb F_2^{H\times J}`$
of rank r. The right probe bank B has J arbitrary bits. Choose full-row-rank
matrices M,N giving a rank factorization

```math
F=M^{\mathsf T}N,
\qquad M\in\mathbb F_2^{r\times H},\quad
N\in\mathbb F_2^{r\times J}.
```

This is terminating classical preprocessing, with no efficiency claim.
The literal Clifford probe and desired selected phase are

```math
\begin{aligned}
P_F|A,B\rangle&=(-1)^{A^{\mathsf T}FB}|A,B\rangle,\\
D_F|x,B\rangle&=(-1)^{e_x^{\mathsf T}FB}|x,B\rangle.
\end{aligned}                                                    \qquad(1)
```

All statements allow arbitrary private work W and external references.
A source U acts only on x,A,W. It is not required to return any of them
at its own boundary. Its exact probe contract is

```math
U P_F U^\dagger P_F=D_F\otimes I_{A,W}.             \qquad(2)
```

The ideal indicator I on x,A satisfies (2). Consequently, with
$`R=I^\dagger U`$, equation (2) is equivalent to

```math
[R,Z_A^{Fb}]=0\quad\hbox{for every }b\in\mathbb F_2^J.
```

Since N is surjective, the operators tested are precisely
$`Z_A^{M^{\mathsf T}v}`$ for all r-bit v. Their simultaneous
eigenspaces are the syndrome fibers $`MA=s`$. Thus (2) holds if and
only if

```math
U=I G,\qquad G=\bigoplus_{s\in\mathbb F_2^r}G_s,    \qquad(3)
```

where each $`G_s`$ is an arbitrary unitary on the entire fiber
$`\{A:MA=s\}`$ together with x,W. In particular, G may mix
different physical A values, change the address, change private work,
and use noncommuting changes of basis. It need not be monomial.
The dimension of an A fiber is $`2^{H-r}`$. When r is zero the
probe and target are identity and the contract places no restriction
on U; the following row query is empty.

For full row rank $`r=H`$, the syndrome determines A, and (3) is
exactly the universal bank-diagonal gauge of the linked chapter.
This conclusion uses rank, not sparsity or a claimed difficulty of F.

## 2. Paid equivalence to an r-bit row query

Complete the rows of M to an invertible H-by-H binary matrix C. The
reversible linear coordinate change

```math
C_A|A\rangle=|CA\rangle=|S,K\rangle,
\qquad S=MA,
```

has r syndrome wires S and H-r remaining wires K. Gaussian elimination
implements it by CNOTs and swaps. For either orientation the bound
$`G_C\le4H^2`$ is sufficient: at most $`H(H-1)`$ row additions
and H swaps, each swap using three CNOTs. It uses no helper, clean bit,
or T gate. This coordinate change is paid when physically used.

Define the prescribed r-bit row and its complete XOR query by

```math
a(x)=Me_x,\qquad
Q_a|x,S\rangle=|x,S+a(x)\rangle.                    \qquad(4)
```

The canonical source

```math
J_F=C_A^\dagger(Q_a\otimes I_K)C_A                  \qquad(5)
```

has the same syndrome change as I. Therefore $`I^\dagger J_F`$
preserves every syndrome fiber, and (3) can equivalently be written
$`U=J_F G`$ with G preserving those fibers. In C coordinates, U
is precisely a complete row query followed on its input side by an
arbitrary unitary preserving S. No row lookup is supplied by knowing
the classical matrices M,F.

To extract the row query, introduce an arbitrary r-bit output V and
use the Clifford probe

```math
P_M=\prod_{j=1}^r\prod_{i:M_{ji}=1}\mathrm{CZ}(A_i,V_j).
```

The gauge commutes with $`P_M`$, so

```math
U P_M U^\dagger P_M=(-1)^{a(x)\cdot V}I_{A,W}.       \qquad(6)
```

Hadamards on V before and after give the complete row query, returning
x,A,W on every input. A chronological emitted word is

```math
H_V,\ P_M,\ U^\dagger,\ P_M,\ U,\ H_V.              \qquad(7)
```

Let $`\nu(M)`$ be the number of nonzero entries of M. If common
forward/actual-inverse source bounds are $`T_U,G_U,d_U,b_U`$,
the row-query wrapper satisfies

```math
\begin{aligned}
T_a&\le2T_U,&D_T(a)&\le2d_U,\\
G_a&\le2G_U+6\nu(M)+2r,&b_a&\le H+b_U.
\end{aligned}                                                    \qquad(8)
```

Each CZ is H,CNOT,H. The entire old bank A is borrowed and charged;
V is the specified output. The source's actual inverse is essential.
The direct parity probe avoids emitting C in this direction, and its
Clifford overhead is at most $`6rH+2r`$.

Conversely, a complete native implementation of (4) on arbitrary S,
with common bounds $`T_a,G_a,d_a,b_a`$, supplies (5) with

```math
T_U\le T_a,\quad D_T(U)\le d_a,\quad
G_U\le G_a+2G_C,\quad b_U\le b_a.                   \qquad(9)
```

The specified A bank supplies S and K, and is excluded from private
source width; K is untouched by the row query. Equations (8)--(9)
are a resource equivalence with the displayed H-bit borrowed bank and
quadratic Clifford overhead. They do not claim that an r-bit query
needs only O(r) work when extracted from an H-bit source.

## 3. Uniform approximation, including gauge and return errors

Suppose $`\widetilde U`$ is within delta in full operator norm of
some $`J_FG`$ from (3)--(5). Its actual adjoint has the same error.
The wrapper (7) has error at most $`2\delta`$, with exactly the
common ledger (8). All address, bank, private-work and reference errors
are included, without assuming exact intermediate return.

There is also a direct probe formulation. If $`\epsilon\lt1`$ and

```math
\|\widetilde U P_F\widetilde U^\dagger P_F-D_F\otimes I\|
\le\epsilon,                                                   \qquad(10)
```

then there exists a syndrome-preserving G such that

```math
\|\widetilde U-J_FG\|\le2\epsilon.                  \qquad(11)
```

Indeed, $`R=J_F^\dagger\widetilde U`$ approximately commutes,
with error at most epsilon, with every $`Z_A^{M^{\mathsf T}v}`$.
The mathematical average

```math
A_0=2^{-r}\sum_vZ_A^{M^{\mathsf T}v}RZ_A^{M^{\mathsf T}v}
```

preserves the syndrome fibers and is within epsilon of R. It is a
contraction whose singular values lie in $`[1-\epsilon,1]`$.
Its polar unitary G is fiber-preserving and within epsilon of A_0.
This proves (11), without a rank-dependent loss. Neither the average
nor the polar factor is a compiled operation. In fact the direct
extraction has the sharper error epsilon: decomposing (10) over B
tests every v=Nb, and N is surjective. The operator error of the
$`P_M`$ probe in (6) is therefore exactly the same as that of the
$`P_F`$ probe, not a sum over generators. A supplied row-query
approximation of error delta gives a source approximation of error
delta by (9), and its original probe error is at most $`2\delta`$.

For full row rank, M is invertible. Conjugating the extracted row
query by the paid linear coordinate change M on V turns (4) into
$`V\mapsto V+e_x`$. This adds at most $`8H^2`$ Clifford gates,
zero T depth and no initialized work. Consequently a full-row-rank
table probe with the requested linear T/dirty and quadratic Clifford
envelopes already supplies the complete indicator, with constant
source-call factors and the same direct probe error epsilon.

## 4. A safe whole-table composition uses completed probes

The gauge in (3) can change x. Replacing a returned-address indicator
inside an existing interleaving by U is therefore not justified by
(2). A complete probe can be composed safely as follows.

For a single-output table, split the original address into x,y and
take its table matrix F. Keep arbitrary A,B and an arbitrary phase
output t. Let

```math
P=(-1)^{t(MA)\cdot(NB)},\qquad
a(x)=Me_x,\quad b(y)=Ne_y.
```

Let U act only on x,A and its private pool, with ideal source
$`J_FG`$ preserving MA on the input side. Let V act only on y,B
and a disjoint private pool, with the corresponding ideal right
source shifting NB by b(y), up to its own syndrome-preserving gauge.
Define the completed left probe and the whole word by

```math
E=UPU^\dagger P,
\qquad W=VEV^\dagger E^\dagger.                    \qquad(12)
```

For the ideal sources,

```math
E=(-1)^{t a(x)\cdot NB}I_{A,W_L},\qquad
W=(-1)^{t a(x)\cdot b(y)}I_{A,B,W_L,W_R}.           \qquad(13)
```

The first completed word returns every left register. The right
gauge commutes with its ideal phase because it preserves NB and acts
on neither x nor t. This proves the second equality even though
neither individual source need preserve its address. Hadamards on t
give the complete one-bit XOR query. Every occurrence of an inverse
in (12) is the actual emitted adjoint, including the actual E adjoint.

R91's private-helper bilinear phase implements P with
$`T_P\le8H`$, $`D_T(P)\le8`$,
$`G_P\le16H+2\nu(F)`$, and H arbitrary private wires. The
Clifford estimate uses literal compute/uncompute of each phase parity,
without relying on Gray-walk cancellation. Hence a safe common ledger
for (12), including the output Hadamards, is

```math
\begin{aligned}
T_W&\le4T_U+2T_V+32H,\\
G_W&\le4G_U+2G_V+64H+8\nu(F)+2,\\
D_T(W)&\le4d_U+2d_V+32,\\
b_W&\le2H+J+b_U+b_V.
\end{aligned}                                                    \qquad(14)
```

The width includes both borrowed banks, disjoint source pools and
the H private bilinear wires. No clean bit is introduced. The right
source pool can remain changed while E runs, so it is not silently
reused as the left source pool. If U,V have full source-gauge errors
$`\delta_U,\delta_V`$, replacing the six source appearances gives
whole-query error at most $`4\delta_U+2\delta_V`$. This hybrid
also includes return and leakage errors on every later reuse.

The same ordered construction permits R91's selected, batched
multi-output P. Choose the left syndrome from the span of all columns
of every actual matrix $`F_{c,j}`$, and the right syndrome from
the span of all their rows. Both sources exclude the high selector c
and the output bank. Each completed left probe then precedes the
right-source conjugation. There are four left-source appearances,
two right-source appearances and four complete bilinear-word
appearances, outside all output batches and traversals. Thus the
source counts are constant, as required by the current R98 transfer;
the native bilinear count, Clifford and width charges retain the R91
allocation orders.

### The direct branch and the complete R98 transfer

R91 also has a direct branch when the q-bit table has Q rows,
$`q\ge Q`$, and its pool fits a Q-bit source bank. For this
branch take F to be the Q-by-q matrix of output rows. Use its source
and the original $`P_F`$ probe with the actual q-bit output as the
probe bank, with output Hadamards at both ends. This uses two source
calls, at most $`6\nu(F)+2q`$ further Clifford gates, and Q plus
the source's private width in borrowed wires. Source-gauge error delta
gives complete query error at most $`2\delta`$. No initialized
r-bit syndrome or extra row-lookup output is needed.

For an explicit conditional hypothesis, suppose that for every fixed
requested delta and every actual table family at a query position,
each length-H left source and length-J right source just described
has full-operator source-gauge error at most delta, zero clean work,
actual adjoints, and common bounds

```math
\begin{aligned}
T_U,b_U&=O_\delta(H),&G_U&=O_\delta(H^2),&
d_U&=O_\delta(d_\delta(\log_2H)),\\
T_V,b_V&=O_\delta(J),&G_V&=O_\delta(J^2),&
d_V&=O_\delta(d_\delta(\log_2J)).
\end{aligned}                                                    \qquad(15)
```

The nondecreasing function $`d_\delta\ge1`$, compiler constants
and pool constants are uniform in the tables, their selected blocks
and their output lengths. The direct branch has the analogous
length-Q source bound. This hypothesis concerns the source for the
joint syndrome of all matrices used by the completed bilinear word;
separately chosen gauges for individual batches are insufficient.

The nested word uses at most six approximate source appearances per
query, outside all batches and traversals. Its four complete bilinear
occurrences are a constant number, just as in R91. Enlarge the R91
allocation constant to reserve both live source pools. The retained
rounding and allocation proof then gives

```math
\begin{aligned}
T&=O_\delta(\sqrt{Qq}+Qq/B),&G&=O_\delta(Qq),\\
D_T&=O_\delta(d_\delta(a)+Qq/B^2),&
\|\widetilde U_f-U_f\otimes I\|&\le6\delta,
\end{aligned}                                                    \qquad(16)
```

for $`Q=2^a`$ and the same form of sufficient minimum pool
$`B\ge C_\delta(a+2)`$. The two source Clifford terms are
absorbed by $`H^2+J^2=O(Qq)`$ in the retained bilinear branches.
In the direct branch $`Q^2\le Qq`$ and
$`\nu(F)\le Qq`$. The zero-row-address query is a fixed
Clifford word. All pools fit the chosen physical allocation B after
the constant enlargement; neither source changes the other source's
original address or its bank.

Use only these completed long-query replacements in the current
R98 construction, leaving its short tables, paired routers and other
native components at their retained implementations. The
[approximate-transfer call audit](APPROXIMATE_QUERY_FRAME_TRANSFER.md#3-the-bounded-r98-call-count-and-final-error)
gives $`K_*=4C_E+56`$ completed long-query appearances, including
all actual inverses. At fixed final precision set

```math
\eta=2^{-L},\qquad \delta=\eta/(24K_*).
```

The added error is at most $`6K_*\delta=\eta/4`$. The complete
real and complex initialized-isometry errors are respectively at most
$`149\eta/512`$ and $`789\eta/1024`$. All dirty-return,
address-return and sole-clean-flag errors are included by full-unitary
telescoping; later sources are not assumed to start with freshly
returned work. The resulting conditional depth is

```math
D_T=O_L\!\left(N/b^2+d_{\delta_L}(2n+10)
                         +\log_2\log_2(n+4)\right).              \qquad(17)
```

The same circuit retains $`T=O_L(\sqrt N+N/b)`$,
$`G=O_L(N)`$, one clean qubit, and the original family-specific
literal dirty thresholds. R96 supplies the unchanged small-dimension
and narrow-width fallbacks; delta and all new constants depend on
fixed L, not on dimension or the target table. No growing-L improvement
is asserted. Equations (15)--(17) are a paid transfer of a supplied
table-specific source, not a compiler for that source.

## 5. What this changes, and what it leaves unpaid

The original one-output table can already have full row rank in a
balanced split: take $`H=J`$ and $`F=I_H`$. It describes equality
of the two half-addresses. A complete alternative adds x into y by
CNOTs, applies the retained read-only zero-conjunction to y and the
arbitrary output, and actually undoes the CNOTs. For
$`H=2^p`$ this has polynomial-in-p T, Clifford and dirty costs and
$`O(\log(p+2))`$ T depth. Thus F is a rank witness, not a
hard-query example or a lower bound on the complete word. It shows
why a small number of final outputs does not imply a small per-source
syndrome.
For several matrices the relevant rank is the joint column span, so
separate low-rank claims would not suffice either.

The actual R98 long-query positions admit prescribed signs, masks,
metadata and matching tables, including arbitrary original sign
tables. Their general theorem imposes no low-rank promise on these
balanced matrices. Consequently (8)--(14) provide no uniform
sublogarithmic source for those positions. The constant-accuracy
approximate transfer remains available, but its faster source premise
is still unpaid. The complete-frame resource frontier does not move.

The remaining constructive freedom is sharper than simply allowing
address or bank motion: synthesize the complete two-bank table word
without requiring its left source to answer every arbitrary B probe
on its own. The selected right address y is then part of the joint
native problem, rather than being replaced by a universal arbitrary
bank at an intermediate boundary. A successful proposal must still
price the complete word, actual inverse and all dirty-return errors;
the classification above neither constructs nor obstructs it.

The [exact fixture](../../tests/test_table_specific_indicator_interfaces.py)
uses a rank-two matrix with four left rows, a nontrivial paid linear
coordinate change, a nonlinear two-bit row query, and a source gauge
that mixes address and kernel-bank wires. It checks every column of
the native extraction and its actual inverse, literal phases and
gate counts. The same source fails an unpromised individual-bank
probe. This verifies the interface distinction, not its all-size
approximation proof or a new asymptotic compiler.
