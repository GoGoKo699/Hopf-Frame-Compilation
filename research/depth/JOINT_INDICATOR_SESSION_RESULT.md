# Revision: complete joint indicators and growing native displacement

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous checkpoint](EXTERNAL_CENTER_SESSION_RESULT.md) ·
[Current R98 transfer](QUERY_DEPTH_FRAME_TRANSFER.md) ·
[Constant-accuracy extension](APPROXIMATE_QUERY_FRAME_TRANSFER.md)

The [subsequent native-sharing revision](NATIVE_SHARING_SESSION_RESULT.md)
retains these results and continues the approximate indicator search from
PR #196. Its row-count and table-specific interface results have separate
proof homes; the present chapter remains the PR #196 checkpoint.

This session began on 19 October 2026 by fetching and checking main against
the reviewed PR #195 commit
`b598e4aa55c308f11f183a11ab92f0e4fc70463a`. They matched, and the working
tree was clean. The main constructive work returned to the complete
approximate dirty indicator. A focused parallel lower-bound scout resolved
the growing-total-displacement question in the actual no-middle-pair
reduction. The complete-frame resource frontier did not move.

## 1. Destination, strongest results and component scope

| Role | Statement |
| --- | --- |
| Destination | Optimal general T-depth for the prescribed complete frame, with all resource and operator requirements on the same circuit |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ |
| Literal domain | $`n\ge1,L\ge6`$; one clean; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error at most $`2^{-L}`$ |
| Selected theorem | R79 remains Result D; settled count frontiers and allocations are unchanged |
| Conditional route | Current R98 permits temporary address motion inside completed joint queries. Its approximate extension needs only constant source accuracy at fixed frame precision |
| New lower component | The actual no-middle-pair three-layer reduction has a uniform Boolean-degree bound, including growing total displacement and active ranks |
| New upper components | Maximal universal-probe gauge equivalence, paid selected-swap equivalence, and a noncommuting dirty geometric scalar with its complete-bank precision obligation |
| Remaining scope | Original physical middle pairs, higher-layer asymptotics, and a faster complete approximate indicator remain open |

At fixed accuracy and square-root dirty width, the general depth interval
is still a constant lower bound versus a logarithmic upper bound. A
three-layer component theorem does not close that asymptotic interval.
R74's exact-rounding floor remains unchanged. No supplied precision state,
QRAM, uncharged evaluator, measurement or reset is introduced.

## 2. The complete constructive target and paid transfer

For p address bits and $`m=2^p`$, the target is

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle.
```

At every fixed requested tolerance tau, the source must have zero clean
work, $`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$, and
$`D_T=o(\log(p+2))`$. Accuracy is full operator norm on arbitrary
address, bank, private dirty work and references. The actual inverse and
all leakage and return errors are charged. A single-predicate success
probability or a maximum over separately approximated bank rows is not
this contract.

The earlier additive-filter, Gaussian-conjugator, phase-reflection and
coherent-frame limits were read before choosing mechanisms. The present
search allows noncommuting bases and temporary address/bank motion beyond
those hypotheses. The current restored-address R98 assembly compiles
each shared-address query position as a complete joint word; it does not
interleave arbitrary returned-address sources gate by gate.

If a source with the required envelope and depth function d exists, the
retained constant-call transfer gives, at fixed L,

```math
\begin{aligned}
T&=O_L(\sqrt N+N/b),\qquad G=O_L(N),\\
D_T&=O_L\!\left(N/b^2+d_{\tau_L}(2n+10)
                  +\log_2\log_2(n+4)\right).
\end{aligned}
```

The literal real/complex widths and one-clean budget stay unchanged.
No such faster source is supplied by this session.

## 3. Mechanisms and their first unpaid operations

### A universal probe permits a larger boundary gauge

The [gauge theorem](WEAKER_INDICATOR_INTERFACES.md) allows a source
$`U=I_pG`$, where G is an arbitrary computational-bank-controlled unitary
on the address and private work. Neither the address nor private work
need return at this source boundary. A two-call Clifford bilinear probe
extracts the complete literal indicator and returns all systems.

This is exactly the maximal gauge for that universal probe. Approximate
probe error epsilon below one implies distance at most twice epsilon
from the gauge class, with no dimension factor. Gauge-class error delta
gives indicator error at most twice delta. The paid extraction has

```math
T\le2T_U,\quad G\le2G_U+8m,\quad
D_T\le2d_U,\quad b\le m+b_U.
```

R91 uses at most four complete indicators per query, outside all output
batches. Taking $`\delta=2^{-L}/(32K_*)`$ with the retained
$`K_*=4C_E+56`$ charges at most one quarter of the final error budget.
This freedom changes admissible candidate words but is constant-call
equivalent to the original task. Merely leaving a larger gauge unpaid
does not improve the frontier.

### A selected swap is also a complete sufficient interface

The [routing theorem](JOINT_ROUTING_INDICATOR_INTERFACES.md) derives the
indicator from two complete selected-swap calls. Conversely, native
gather/scatter/gather uses six indicators and noncommuting target bases,
with

```math
\begin{aligned}
T_S&\le6T_I+36m,&D_T(S)&\le6d_I+24,\\
G_S&\le6G_I+64m+26,&b_S&\le m+b_I.
\end{aligned}
```

Its error is at most six times the source error, independent of bank
size. A direct pruned switch tree has linear count and no helper, but
retains a depth proportional to p. A fixed Clifford cannot make every
full bank translation simultaneously diagonal: their symplectic
differences have nonzero products, whereas diagonal-Clifford differences
have pairwise zero products. That is a representation obstruction, not
a native-depth lower bound.

There is also a precise limit on classically closed recurrences. With r
monomial reversible rounds of coordinate degree at most k and arbitrary
affine address/bank/work mixing between them, full error below square
root two forces $`k^r\ge p`$. Thus quadratic classical rounds still
need logarithmically many rounds. This counts those specified rounds,
not native T layers; nonmonomial interference across the boundaries is
outside the premise.

The first unpaid operation is a joint synthesis that replaces the
successive selection decisions. Grouping them under a selected-router
name keeps the same missing operation unless its complete native
recurrence improves.

### A dirty geometric scalar changes the signal mechanism

The [geometric-loop study](DIRTY_GEOMETRIC_PHASE_LOOPS.md) uses a
noncommuting two-level loop and an echo that cancels dynamical phases.
Its ideal winding sign is the same on both signal eigenstates, so the
earlier rank-one-reflection fixed-space obstruction does not apply.
The finite-time comparison is on the entire dirty signal, rather than
on one prepared eigenvector.

Independent row loops nevertheless retain a full-bank precision cost.
For fixed address and output values their target-relative SU(2)
eigenphase errors must sum to a small arc, not merely be individually
small. The native digital construction also pays every controlled
rotation, approximation, interval and actual inverse. Its explicit
recurrence misses the target envelope. No limitation on a collective
joint word is inferred from this independent-row calculation.

## 4. The lower uncertainty that changed

The [native-orbit proof](GROWING_DISPLACEMENT_NATIVE_ORBITS.md) starts
from the actual one-clean source and its literal global Boolean phase.
In the remaining high-rank no-middle-pair case, let U be the total middle
displacement, K its orthogonal complement in the original native
coordinates, and v the common preimage supplied by PR #195. The full
cubic identity implies

```math
x*y+\beta(x,y)v\in K\qquad(x,y\in K).
```

Associativity makes $`K+\langle v\rangle`$ a coordinatewise algebra.
Its coisotropy forces actual coordinate blocks of size one or two.
Genuine-center fixation and the full cubic tensor then eliminate
cross-couplings between the pair blocks. The remaining action consists
of actual native pair permutations and one common exceptional
transvection.

The growing permutation factors are not discarded. They form the middle
band of an actual three-native-band residual. Normalizing the exceptional
fixed family leaves at most three interacting computational directions.
The retained partial-diagonal feature theorem gives

```math
\deg g\le D_3+576
```

in the high-rank case, and the other control ranks give the combined bound

```math
\deg g\le\max\{D_3+576,D_{\rm ins}(0,0,0,3)\}.
```

This is uniform in active ranks, control support, total displacement and
dirty width. XOR degree is one less only after checking the actual
output-Hadamard sign representation. The proof excludes growing Boolean
degree in this source class; it does not bound the displayed total
displacement itself. Its logical pair orbits differ from the original
physical clean-transition middle pairs that remain outside the reduction.

The proof retains literal phases, original native lists, full dirty and
reference inputs, clean return and actual inverses. A Gram-only model
does not satisfy its hypotheses. The exact component fixture includes
one such spurious cross-coupling that the full tensor excludes; it is
not represented as a completed Boolean source family.

## 5. Decision and next concrete question

The main constructive question remains a collective native word for the
complete approximate indicator. The new gauge characterization gives a
specific relaxed design target: synthesize only the required universal
bank-Pauli transport while allowing arbitrary bank-controlled address/work
motion, then use the paid two-call extraction. This is an equivalent
target with more implementation freedom, not a claimed shortcut.

The geometric candidate should reopen only with a collective mechanism
that avoids independent-row eigenphase accumulation and supplies the
native controlled rotations within the linear count budget. Repeating
the same independent loops at greater precision retains the diagnosed
cost. A genuinely weaker, table-specific frame interface must instead
identify actual R98 positions and prove their complete ordered
replacement; universal gauges and complete selected swaps do not provide
that weakening.

The next lower question is the original middle-pair branch at growing
active ranks: can its actual paired native factors support an analogous
finite residual argument while retaining their signs and relations?
The higher-layer dimension-growing question remains active independently.
Neither component closure nor a stopped mechanism retires optimal
general depth as the destination.

The reviewed Hopf-QBP documentation-only change is reconciled in
[SYNC](../../SYNC.md), advancing the tracked main to `99e5ee5` while
preserving both earlier reconciliations and all scientific interfaces.
No upstream source code or scientific correction was imported.

## 6. Verification boundary

The degree and source-interface statements have analytic proofs and
independent review. Bounded exact fixtures check native phases, whole
physical columns, actual inverses, returned arbitrary banks and the
specific algebra identities. They do not establish asymptotics by
extrapolation. Required repository and rendered-presentation gates are
recorded on the reviewed pull-request tree. Final manuscript drafting
remains separate.
