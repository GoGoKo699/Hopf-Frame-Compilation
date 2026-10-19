# Which incomplete indicator interfaces still contain the full operation

[Approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Complete bilinear queries](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Coherent-frame interface](PHASE_SIGNAL_QUERY_LIMITS.md)

A source boundary need not return the address or private work to supply
a complete indicator. This chapter identifies a large exact gauge
freedom and proves its full-operator approximate transfer. The freedom
is real, but does not give an asymptotically easier source: a paid
two-call Clifford probe recovers the complete indicator. Conversely the
complete indicator supplies the interface. This is an interface
classification, not a faster native compiler or a lower bound on
unrestricted synthesis. R96 and the complete-frame frontier are unchanged.

## 1. A source may mix the address and all private work

Let x have p bits, let Y have $`m=2^p`$ bits, and let W be an arbitrary
private work system. Define the target indicator on x,Y by

```math
I_p|x,Y\rangle=|x,Y\mathbin\oplus e_x\rangle.
```

For every computational Y value, choose an arbitrary unitary $`G_Y`$
on x,W. No efficient description or compiler for these unitaries is
assumed. Put

```math
G=\sum_Y|Y\rangle\langle Y|\otimes G_Y,
\qquad U=(I_p\otimes I_W)G.                         \qquad(1)
```

U may mix different address values, entangle the address with W, change
W, and carry arbitrary phases. It need not implement a classical
permutation or return either x or W at its boundary. This is more
freedom than a relative phase on the usual returned-address indicator.
The only premise is the completed factorization (1) on the full space;
there is no promise restricted to an initialized work state.

Introduce a separate arbitrary m-bit bank V. The Clifford bilinear
phase and desired selected phase are

```math
B=\prod_{j=0}^{m-1}\mathrm{CZ}(Y_j,V_j),
\qquad D_{x,V}|x,V\rangle=(-1)^{V_x}|x,V\rangle.
```

G commutes with B because it preserves every computational Y sector
and does not act on V. Therefore, in rightmost-first operator order,

```math
U B U^\dagger B
  =I_p B I_p^\dagger B
  =D_{x,V}\otimes I_{Y,W}.                           \qquad(2)
```

The second equality follows on every computational x,Y,V value from
$`(Y+e_x)\cdot V+Y\cdot V=V_x`$. It is an operator identity; the
first equality does not require U to preserve any of those sectors.
The actual U adjoint returns the possibly changed work in the complete
word. Hadamards on V before and after (2) give the literal indicator
into V, with x,Y,W all returned, including arbitrary reference inputs.
A chronological native implementation is

```math
H_V,\ B,\ U^\dagger,\ B,\ U,\ H_V.                 \qquad(3)
```

## 2. This is the maximal exact gauge for the universal probe

Suppose an arbitrary unitary U on x,Y,W satisfies (2) for the universal
coherent probe bank V. Restricting V to each computational value v gives

```math
U Z_Y^v U^\dagger
 =I_p Z_Y^v I_p^\dagger
\qquad(v\in\mathbb F_2^m).                          \qquad(4)
```

Thus $`R=(I_p^\dagger\otimes I_W)U`$ commutes with every product
$`Z_Y^v`$. The common eigenspaces of these products are exactly the
computational Y sectors, with the entire x,W system as multiplicity.
Consequently R has precisely the form G in (1). The characterization
allows arbitrary Y-controlled changes of the address and private work;
it imposes no unnecessary return condition on those systems.

An interface that specifies only selected probes can have a larger
commutant. That does not establish the universal identity (2). In
particular, correctness of separately tested generators with individual
error epsilon does not give an address-independent error for all their
products: a telescoping bound can grow with the number of generators.
The coherent-bank operator norm below tests all products at once.

## 3. Approximation, native resources and the actual inverse

There are two useful quantitative formulations. First, suppose a native
unitary $`\widetilde U`$ obeys, for some G of the form (1),

```math
\|\widetilde U-(I_p\otimes I_W)G\|_{\mathrm{op}}
 \le\delta.                                          \qquad(5)
```

No description or implementation of G is needed in (3). Two-call unitary
telescoping, using the actual $`\widetilde U^\dagger`$, gives complete
indicator error at most $`2\delta`$. This includes the return errors
of x,Y,W; no intermediate exact return is inferred or used.

Alternatively, a compiler can directly certify the universal probe
error

```math
\left\|\widetilde U B\widetilde U^\dagger B
       -D_{x,V}\otimes I_{Y,W}\right\|_{\mathrm{op}}
 \le\epsilon.                                       \qquad(6)
```

Then (3) has error at most epsilon. This statement is useful as an audit
of a joint word, but its certification must include the full coherent
bank V. Merely establishing a basis-state transition probability does
not prove (6).

The two notions are uniformly equivalent up to constants. If
$`\epsilon\lt1`$ in (6), put $`R=(I_p^\dagger\otimes I_W)\widetilde U`$.
Equation (6), decomposed over V, implies

```math
\|RZ_Y^v-Z_Y^vR\|_{\mathrm{op}}\le\epsilon
\qquad\hbox{for every }v.
```

Apply the diagonal twirl, as a mathematical average rather than a
compiled circuit:

```math
A=2^{-m}\sum_v Z_Y^vRZ_Y^v.
```

A preserves each Y sector, $`\|A\|\le1`$, and
$`\|A-R\|\le\epsilon`$. Since R is unitary, every singular value
of A lies in $`[1-\epsilon,1]`$. Its polar unitary
$`G=A(A^\dagger A)^{-1/2}`$ is therefore well defined and preserves
every Y sector. Moreover $`\|G-A\|\le\epsilon`$, so

```math
\|\widetilde U-(I_p\otimes I_W)G\|_{\mathrm{op}}
 \le2\epsilon.                                       \qquad(7)
```

The twirl has no dimension-dependent loss. Neither its exponentially
many formal terms nor the polar unitary is emitted as a native operation.
They prove the existence of the gauge in (7), not a synthesis algorithm.

Let common forward/actual-inverse source bounds be
$`T_U,G_U,d_U,b_U`$, where $`b_U`$ counts W and excludes x,Y.
Compile each controlled Z as H,CNOT,H. The two B occurrences and two
Hadamard endpoints in (3) give the simultaneous bounds

```math
\begin{aligned}
T_{\mathrm{ind}}&\le2T_U,&
G_{\mathrm{ind}}&\le2G_U+8m,\\
D_T(\mathrm{ind})&\le2d_U,&
b_{\mathrm{ind}}&\le m+b_U.
\end{aligned}                                        \qquad(8)
```

The borrowed old Y bank is charged in the last line; V is the specified
output bank. No clean qubit is added. The actual inverse of the complete
word satisfies the same error estimate, and the common Clifford source
bound must cover its inverse in the chosen native alphabet. Subsequent
simplification of adjacent Clifford inverse gates is optional and does
not improve the stated asymptotics.

## 4. Paid entry into the current R98 transfer

Suppose that for every fixed requested delta a zero-clean compiler of
(5) supplies

```math
T_U=O_\delta(m),\quad b_U=O_\delta(m),\quad
G_U=O_\delta(m^2),\quad d_U=O_\delta(d_\delta(p)),
                                                               \qquad(9)
```

uniformly over p, with $`d_\delta\ge1`$ nondecreasing. Equation (8)
then gives the complete indicator required by the R91 assembly, at
error $`2\delta`$. That assembly uses at most four completed indicator
appearances per word query, outside all output batches and traversals.
Its complete word-query error is at most $`8\delta`$, with the same
R91 count and width envelopes after changing constants. The source's
quadratic Clifford term is absorbed by the retained allocation proof:
$`H^2+J^2=O(Qq)`$ in the bilinear branch, and $`Q^2\le Qq`$ in
the direct branch. The indicator extraction is completed before any
query-level reuse; internal source changes of the address require no
parallel interleaving assumption.

Use the current restored-address R98 construction and the call count
$`K_*=4C_E+56`$ from the approximate-transfer chapter. At final fixed
precision $`\eta=2^{-L}`$, request

```math
\delta=\eta/(32K_*).                                 \qquad(10)
```

The added frame error is at most $`8K_*\delta=\eta/4`$. The
retained real and complex totals are respectively at most
$`149\eta/512`$ and $`789\eta/1024`$. All approximate source calls
and their actual inverses are included, without assuming that earlier
leakage vanishes at a pool-reuse boundary. The resulting conditional
complete-frame ledger is

```math
\begin{aligned}
T&=O_L(\sqrt N+N/b),\qquad G=O_L(N),\\
D_T&=O_L\!\left(N/b^2+d_{\delta_L}(2n+10)
                      +\log_2\log_2(n+4)\right),
\qquad N=2^n.
\end{aligned}                                        \qquad(11)
```

It retains R98's one-clean allocation, real threshold
$`b\ge17(L+n+7)`$, complex threshold $`b\ge17(L+n+8)`$, and
full initialized-isometry error at most $`2^{-L}`$. The existing R96
fallbacks cover narrow widths and bounded dimensions. No uniform
improvement for growing L is claimed. Under the direct probe hypothesis
(6), request $`\epsilon=\eta/(16K_*)`$ instead; the same proof uses
at most $`4K_*\epsilon`$ additional error.

Conversely a complete indicator already supplies (1) with G equal to
identity and supplies (6) through its two-call probe. Hence (9) with
$`d_\delta(p)=o(\log(p+2))`$ would improve the frontier, but this
source condition is equivalent, up to constant calls and linear extra
Clifford/dirty resources, to the original complete-indicator task.
Allowing these boundary gauges changes what a candidate circuit may
look like; it does not itself resolve the missing asymptotic operation.

## 5. What a genuinely weaker frame interface would have to establish

A table-specific probe can have fewer constraints than (6). To use that
fact, a proposal must name the actual R98 query positions it replaces,
prove the complete ordered replacement identity at those positions, and
bound its simultaneous counts and operator error. R91's general bilinear
assembly allows arbitrary matrices and arbitrary output words. A proof
for one selected probe, for initialized outputs, or for separately
correct generator conjugations cannot be substituted for that assembly
without a new audit. The present theorem establishes neither impossibility
nor a faster compiler for such a table-specific approach.

The [bounded native fixture](../../tests/test_weaker_indicator_interfaces.py)
uses a source that changes and entangles the address and its arbitrary
helper. It verifies the full two-call word and actual inverse on every
physical column, including the literal phase. The maximal commutant,
approximate polar argument and R98 transfer above are analytic results.
