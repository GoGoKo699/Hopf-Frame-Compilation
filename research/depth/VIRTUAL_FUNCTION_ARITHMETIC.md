# Exact arithmetic with virtual Boolean functions and dirty work

[Exact dirty indicators](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences) ·
[Parallel clean-source promotion](PARALLEL_CLEAN_PROMOTION.md) ·
[Native sharing audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md)

Arbitrary Boolean functions of unchanged inputs can control a balanced
addition recursion without being written into initialized registers.
This pays functional addition, odd multiplication modulo a power of two,
interval reversal, and modular addition in polylogarithmic native
T-depth. Their exponential count and dirty width are explicit. These
components do not yet pay the mixed-radix packing required at the
critical clean-source size.

Every circuit below acts on arbitrary principal inputs and arbitrary
borrowed work, including correlations with a reference. Its complete
operator is the stated permutation tensored with the identity on work.
There are zero clean qubits, zero approximation error, and no measurement.
An inverse is the actual reversed native word, with the same resource
bounds. Circuit descriptions and truth tables are fixed before execution.

## 1. The exact Boolean query and its native schedule

Let x consist of at most j unchanged bits, let t be a separate arbitrary
target, and let f be any Boolean function. The retained dirty indicator
I_j adds the equality vector of x into a private arbitrary bank Y. Put

```math
 P_f:\quad t\longmapsto
 t\oplus\bigoplus_{v:f(v)=1}Y_v.
```

The chronological word below is an exact Boolean query:

```math
 I_j,\quad P_f,\quad I_j^\dagger,\quad P_f
 \quad:\quad t\longmapsto t\oplus f(x).                \qquad(1)
```

The old Y contribution cancels and all indicator work returns. The
[indicator ledger](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
therefore supplies universal constants such that

```math
\begin{aligned}
 B(j)&=C2^j(j+3)^3,\\
 T_{\mathrm B}(j),G_{\mathrm B}(j),b_{\mathrm B}(j)&\le B(j),\\
 d_{\mathrm B}(j)&\le C\log_2(j+3).
\end{aligned}                                                       \qquad(2)
```

Here b includes the private indicator-output bank as well as its work.
The bound includes both parity CNOT words and is valid for the actual
inverse. A constant predicate is handled directly.

The raw x wires are unchanged at completed Clifford boundaries; their
shared interactions are CNOT controls. All physical non-Clifford targets
are private. Thus a logical layer of separate Boolean queries may share
unchanged input words when its principal targets and private pools are
disjoint. The Clifford gates can be serialized within their Clifford
blocks and are counted in (2); no constant Clifford-depth assumption or
initialized fanout copy is used. A principal target cannot be a changing
control of another concurrent query.

## 2. Replacing a formal control by a Boolean function

[Remaud and Vandaele, *Ancilla-free Quantum Adder with Sublinear
Depth*](https://arxiv.org/html/2501.16802v2), Section 7.1, give an exact
incrementer with one dirty helper, total X/CNOT/Toffoli depth
$`O(\log^2(s+2))`$ and count $`O(s\log(s+2))`$ on an s-bit word.
The helper returns for arbitrary inputs. The total logical depth is
needed here, since adding a control also makes Clifford gates nonlinear.

First formally control every elementary gate of this incrementer by a
read-only bit c. This includes every X and CNOT, not only its Toffolis.
For a stable Boolean function c(x), replace each resulting elementary
gate by (1), with predicate c(x) times the original gate's at most two
local controls. The target is always outside x. This is a substitution
at the classical logical-gate level; no compiled T gate is controlled.

Let N bound the total number of stable and active principal bits in the
application. Private pools for the substituted queries give

```math
\begin{aligned}
 L_{\mathrm{inc}}(s)&\le C(s+1)\log(s+2),\\
 T_{\mathrm{inc},c},G_{\mathrm{inc},c},b_{\mathrm{inc},c}
   &\le C L_{\mathrm{inc}}(s)B(N+4),\\
 d_{\mathrm{inc},c}
   &\le C\log^2(s+2)\log(N+6).
\end{aligned}                                                       \qquad(3)
```

Allocating a separate pool to every logical gate is a conservative width
bound; reuse can only decrease it. At one original logical layer, local
supports are disjoint. Only the unchanged x word is newly shared, so
Section 1 gives the claimed native parallel schedule. The old incrementer
helper and every new query pool return exactly. No physical wire is
initialized to c(x).

## 3. Functional addition into an arbitrary register

Let x be unchanged and let y be an arbitrary s-bit word. For any function
$`f(x)\in\mathbb Z/2^s\mathbb Z`$, the desired total permutation is

```math
 \mathrm{ADD}_f:\quad(x,y)\longmapsto
 (x,y+f(x)\bmod2^s).                                  \qquad(4)
```

Split s=t+u into halves differing by at most one, and write
$`y=l+2^t h`$ and $`f=f_l+2^t f_h`$. First increment h controlled by

```math
 c(x,l)=[l+f_l(x)\ge2^t].                              \qquad(5)
```

Both x and the old low word l are unchanged during this increment, so
(3) applies. Then run ADD of f_l into l and ADD of f_h into h in
parallel. Their output is

```math
\begin{aligned}
 l'&=l+f_l(x)\pmod{2^t},\\
 h'&=h+f_h(x)+c(x,l)\pmod{2^u}.
\end{aligned}                                                       \qquad(6)
```

This is (4). A one-bit leaf is (1). Parallel child calls share only the
unchanged external x word; neither uses the other child's changing
target as a control. The recursive arithmetic helper pools are private.

Let L_A(s) count substituted logical gates. With all principal input
sizes bounded by N, the complete recurrences are

```math
\begin{aligned}
 L_A(s)&\le L_A(t)+L_A(u)+Cs\log(s+2),\\
 d_A(s,N)&\le\max\{d_A(t,N),d_A(u,N)\}\\
 &\quad+C\log^2(s+2)\log(N+6),\\
 L_A(1)&\le C,\qquad d_A(1,N)\le C\log(N+6).
\end{aligned}                                                       \qquad(7)
```

Consequently one actual word and its actual inverse satisfy

```math
\begin{aligned}
 L_A(s)&=O(s\log^2(s+2)),\\
 T_A,G_A,b_A&\le C(s+1)\log^2(s+2)B(N+4),\\
 d_A(s,N)&=O(\log^3(s+2)\log(N+6)).
\end{aligned}                                                       \qquad(8)
```

The width bound includes original incrementer helpers and allocates all
Boolean-query pools independently. Every intermediate truth table is
a classical function of its stated unchanged controls. This pays its
evaluation by (1), without assuming a free arithmetic evaluator.

A stable outer Boolean control is incorporated by replacing f(x) by
c(x)f(x). The same recurrence then implements the controlled addition
with the same asymptotic bounds, including exact work return in its
inactive sector.

## 4. Odd multiplication modulo a power of two

Let a be odd and let an arbitrary n-bit word be
$`z=l+2^r h`$, where $`r=\lfloor n/2\rfloor`$ and $`u=n-r`$.
First apply the functional shear

```math
 h\longmapsto h+
 (a^{-1}\bmod2^u)
 \left\lfloor\frac{al}{2^r}\right\rfloor
 \pmod{2^u}.                                         \qquad(9)
```

Then multiply l by a modulo $`2^r`$ and h by a modulo $`2^u`$,
recursively and in parallel. The two output words are

```math
\begin{aligned}
 l'&=al\pmod{2^r},\\
 h'&=ah+\left\lfloor\frac{al}{2^r}\right\rfloor
       \pmod{2^u}.
\end{aligned}                                                       \qquad(10)
```

Thus their concatenation is exactly $`az\bmod2^n`$. Oddness is used
only to define the inverse in (9); a one-bit multiplication is identity.
No old copy of z is retained. The same construction allows a stable
external word x to specify any coefficient a(x) that is odd at every x.
Each inverse and carry in (9) is paid as part of an arbitrary Boolean
truth table, rather than computed into clean work.

For total principal size N, let L_M(n) count the logical Boolean gates.
Equations (7)–(9) give

```math
\begin{aligned}
 L_M(n)&\le L_M(r)+L_M(u)+C n\log^2(n+2),\\
 d_M(n,N)&\le\max\{d_M(r,N),d_M(u,N)\}\\
 &\quad+C\log^3(n+2)\log(N+6).
\end{aligned}                                                       \qquad(11)
```

Hence

```math
\begin{aligned}
 T_M,G_M,b_M&\le C(n+1)\log^3(n+2)B(N+4),\\
 d_M(n,N)&=O(\log^4(n+2)\log(N+6)).
\end{aligned}                                                       \qquad(12)
```

For a fixed coefficient and N=n this is $`O(\log^5(n+2))`$ native
T-depth. These are full-space exact in-place multipliers, with exponential
count and dirty width, rather than initialized-output multiply-adders.

## 5. Reversing a specified interval with one arbitrary toggle

Let stable x specify endpoints
$`0\le a(x)\le b(x)\le2^s`$. Put

```math
\begin{aligned}
 e(x,y)&=[a(x)\le y\lt b(x)],\\
 J_x(y)&=a(x)+b(x)-1-y\pmod{2^s}.
\end{aligned}                                                       \qquad(13)
```

J_x is a global involution. It maps the indicated interval to itself
by ordinary reversal and therefore preserves e on its complement too.
This includes empty intervals and endpoints equal to $`2^s`$.

Let Q_e toggle a separate arbitrary bit q by e. The chronological word

```math
 Q_e,\quad\mathrm{ctrl}_q J_x,\quad
 Q_e,\quad\mathrm{ctrl}_q J_x                         \qquad(14)
```

returns q and applies $`J_x^{e(x,y)}`$ to y. The two controlled
exponents are $`q\oplus e`$ and q, and e is unchanged between them.
Because J_x is an involution, their sum modulo two is e. This is a
literal dirty-bit cancellation, not a comparison flag initialized to
zero.

The controlled J_x is paid directly: conditionally complement y,
add q to obtain conditional negation, and then add
$`q(a(x)+b(x)-1)`$ modulo $`2^s`$. The additions are instances of
(4) with stable x,q. Every logical gate is accounted for before Boolean
query compilation; no control is attached to compiled T gates.

For total principal size at most N, (14) uses two Boolean queries and
two such controlled involutions. Safe common forward/reverse bounds are

```math
\begin{aligned}
 T_R,G_R,b_R&\le C(s+1)\log^2(s+2)B(N+6),\\
 d_R(s,N)&=O(\log^3(s+2)\log(N+8)).
\end{aligned}                                                       \qquad(15)
```

The extra two bits allowed in B include the arbitrary toggle and local
increment work. All of them return. The interval's endpoints need not
be computed into physical registers.

## 6. Total modular addition for an arbitrary modulus

Fix $`1\le N_0\le2^s`$ and let stable x specify
$`0\le c(x)\lt N_0`$. Apply the interval reversals in this order:

```math
 [0,N_0),\qquad[0,c(x)),\qquad[c(x),N_0).              \qquad(16)
```

For y below N_0, the first reversal sends it to N_0−1−y.
The appropriate one of the last two reversals then gives
$`y+c(x)\bmod N_0`$. Every y at least N_0 is fixed by all three.
Thus (16) is the total permutation

```math
 (x,y)\longmapsto
 \begin{cases}
 (x,y+c(x)\bmod N_0),&y\lt N_0,\\
 (x,y),&y\ge N_0.
 \end{cases}                                         \qquad(17)
```

Its T count, Clifford count and depth are at most three times (15),
while the private pools can be reused. Its error is exactly zero and
its actual inverse has the same complete return contract. A varying
modulus N_0(x) is also allowed when its endpoints and summand satisfy
the displayed inequalities for every stable x.

## 7. Packing scale and the remaining conversion

At the critical packing scale
$`k\asymp\log\log p`$, $`h=2^k`$, and $`n=hk`$, the operations above
have $`2^{O(n)}\mathrm{poly}(n)`$ count and dirty width. For
$`k=\log_2\log_2 p+O(1)`$, these are
$`2^{O(\log p\log\log p)}`$, well below the exponential-in-p allowance
of the retained polynomial-overhead source interface. Their native
T-depth is polylogarithmic in n and therefore $`o(\log p)`$.
This does not imply that a whole packing circuit has that depth: the
number and dependence of its arithmetic calls must also be bounded.

Nor may exponential truth-table cost be hidden when applying this
component to the full addressed-output source. If the controls include
an entire $`2^p`$-bit bank, (2) is doubly exponential in p and is outside
the allowed source budgets. Even a construction on only p address bits
still owes the complete bank action, native sharing and full-bank error.
These arithmetic components do not supply that source.

For $`m,r\ge1`$, an m-bit row x, an r-bit column y, and
$`0\le b\le L=2^r`$, the stable rectangle shuffle needed by the
packing approach is

```math
 F_b(x,y)=
 \begin{cases}
 bx+y,&y\lt b,\\
 2^m b+(L-b)x+y-b,&y\ge b.
 \end{cases}                                         \qquad(18)
```

This is a total permutation of $`2^{m+r}`$ labels. In the nonzero-digit
application, $`b=(2^k-1)^j`$ is positive and odd, as is L−b. On the first rectangle,
write $`z=2^m y+x`$. Its transpose satisfies

```math
 bx+y\equiv bz\pmod{2^m b-1},                         \qquad(19)
```

with the endpoint fixed. The modulus in (19) is not a power of two;
Section 4 does not implement it. The second rectangle has a separate
range and must be handled too.

Three paid modular-add shears and an exchange can give paired modular
multiplication $`(v,w)\mapsto(av,a^{-1}w)`$ on valid modular-register
pairs. A reciprocal bank may move inside a packing word if the packing
condition is uniform over it and the outer inverse restores it.
However, the valid-pair identity does not cover arbitrary out-of-range
labels of a separate dirty bank. Applied to two key subblocks already
known valid, it only permutes their valid rectangle and does not convert
that rectangle into a contiguous mixed-radix interval.

The first unpaid operation remains the complete rectangle-to-rank
conversion (18), or another literal key packing with a sufficiently
shallow total schedule. Out-of-place queries for it and its inverse do
not pay this in-place obligation. No larger sufficient packing window
or qualifying shallow clean source follows from the components alone.
