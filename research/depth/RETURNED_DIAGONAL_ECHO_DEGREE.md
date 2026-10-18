# Returned diagonal source words have bounded Boolean degree

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Read-only layer results](READ_ONLY_THREE_LAYER_RIGIDITY.md)

A fixed number of arbitrary diagonal work sources cannot expose an
arbitrarily high-degree Boolean scalar phase through affine transport of
fully arbitrary work. The theorem permits distinct sources, arbitrary
unit-circle phases, fixed invertible affine work mixing, and all literal
monomial-Clifford phases between source calls. It has a positive full-operator
error radius independent of address length, dirty width and source precision.
The earlier double-translation echo is the four-source special case: its
cubic degree bound now holds below error $`2\sin(\pi/32)`$.

The number of source occurrences is not native T-depth. Every source and
actual inverse must still be synthesized and charged, and a fused native
word may leave this architecture. This result gives no unrestricted native
T-depth recurrence. The explicit cubic native example below is a scope
witness, not a resource improvement.

## General diagonal source words

Let the address $`x\in\mathbb F_2^a`$ be read-only and let all
$`w`$ physical work bits be arbitrary, with zero initialized helpers.
Consider a word with $`q\ge0`$ occurrences of fixed work-diagonal sources

```math
D_j|z\rangle=e^{2\pi i h_j(z)}|z\rangle,
\qquad h_j:\mathbb F_2^w\longrightarrow\mathbb R/\mathbb Z,
\qquad 1\le j\le q.
```

Sources may differ, and may ignore some work bits. Their phases may be
irrational. They do not depend directly on the address: address dependence
enters through the intervening affine transformations of dirty work.
An actual inverse occurrence has its own position in q, with phase
$`-h_j`$ if it is the inverse of an earlier source.

Before, between and after the sources, allow completed read-only monomial
Clifford intervals. On computational columns each interval has the form

```math
|x,z\rangle\longmapsto
 e^{2\pi i c_j(x,z)}|x,M_jz+L_jx+b_j\rangle,
\qquad M_j\in\mathrm{GL}(w,\mathbb F_2).
```

Here the matrices and offsets are fixed. Literal monomial-Clifford phases,
including S phases, CZ phases and address-work cross phases, are retained.
A completed interval may be synthesized in any exact native way, but its
net action must have this form. Nonmonomial completed intervals and nonlinear
work permutations are outside this statement.

Write $`W_x`$ for the resulting work operator. For a prescribed Boolean
function f, put

```math
r=\max\{q,3\},\qquad
\rho_r=2\sin\!\left(\frac{\pi}{2^{r+1}}\right).
```

**Theorem.** If

```math
\max_x\|W_x-(-1)^{f(x)}I_W\|_{\mathrm{op}}\le\varepsilon,
\qquad \varepsilon\lt\rho_r,
```

then the completed affine work permutation is exactly the identity and

```math
\deg_{\mathbb F_2}f\le r-1=\max\{2,q-1\}.
```

Exact full scalar return to a prescribed Boolean phase is included by taking
$`\varepsilon=0`$. In particular, a target of degree $`d\ge3`$ requires at
least $`q\ge d+1`$ source occurrences in this architecture. The degree-two
allowance is necessary: an address CZ is already an allowed monomial
Clifford with no diagonal work source. The theorem does not assert that an
arbitrary returned scalar phase must be Boolean; that additional fact holds
for the specific echo below.

The norm is on the complete physical space, with every dirty input retained.
It includes coherent address superpositions and tensoring with arbitrary
external reference systems. There is no initialized-subspace promise,
measurement, reset, leakage allowance or intermediate return assumption.

## Proof by joint address-work translations

Following every computational column through the word gives

```math
W_x|z\rangle=e^{2\pi i\Phi(x,z)}|R(x,z)\rangle,
\qquad
\Phi(x,z)=c_{\mathrm{Cl}}(x,z)
 +\sum_{j=1}^q h_j(A_jz+B_jx+b_j)
 \quad\text{in }\mathbb R/\mathbb Z.                 (1)
```

Each $`A_j`$ is the product of the preceding invertible work blocks and is
therefore invertible. The offsets are unrestricted. The final map R is
likewise an affine permutation of work for each address.

The accumulated Clifford phase $`c_{\mathrm{Cl}}`$ has zero third additive
differences in the joint variables $`(x,z)`$. To see this directly, a
monomial Clifford phase is a sum of a constant, quarter-integer multiples
of affine Boolean parities, and half-integer multiples of products of two
affine parities. An affine parity divided by four has zero third differences
modulo one: any nonzero difference changes its sign with a factor of two
at each subsequent difference, so the third difference is integral.
A product of two affine parities divided by two has zero third differences
because its Boolean degree is at most two. These facts include repeated,
zero and dependent directions. Substitution of any affine prefix map leaves
this description unchanged, so it applies to the accumulated phase in (1).

Choose arbitrary address directions $`u_1,\ldots,u_r`$. Assign to each
source its own joint direction

```math
v_j=(u_j,A_j^{-1}B_ju_j)\quad(1\le j\le q),
\qquad v_j=(u_j,0)\quad(q\lt j\le r).
```

Translation by $`v_j`$ leaves the argument of source j unchanged. With
$`\Delta_v\Phi(y)=\Phi(y+v)-\Phi(y)`$, the additive differences commute.
Each source term is killed by its assigned factor, and the Clifford phase
is killed because $`r\ge3`$. Thus at every joint base point,

```math
\Delta_{v_1}\cdots\Delta_{v_r}\Phi(x,z)=0
                       \quad\text{in }\mathbb R/\mathbb Z.       (2)
```

This is a joint-cube identity on all work inputs. It does not fix work to a
prepared state or assume that any source is inexpensive.

Since $`\rho_r\lt\sqrt2`$, the operator error hypothesis first forces
$`R(x,z)=z`$ for every column. Otherwise the actual and target output basis
vectors are orthogonal and their difference has norm exactly $`\sqrt2`$,
regardless of their phases. The completed word is therefore diagonal.
Define its error phase relative to the literal Boolean target by

```math
E(x,z)=e^{2\pi i\Phi(x,z)}(-1)^{f(x)},\qquad
|E(x,z)-1|\le\varepsilon.
```

Its principal argument $`\theta(x,z)`$ satisfies

```math
|\theta(x,z)|\le2\arcsin(\varepsilon/2).
```

Only the error ratio is given a local principal argument; no real lift of
an arbitrary source phase is chosen. Exponentiating (2), the signed product
of the actual word phases on this r-cube equals one. If the Boolean
r-fold derivative of f in the address directions were one, the corresponding
target phase product would equal minus one. The signed product of the error
phases would then equal minus one, so the sum of their signed principal
arguments would be an odd multiple of pi. Consequently

```math
\pi\le2^r\,2\arcsin(\varepsilon/2),
\qquad
\varepsilon\ge2\sin\!\left(\frac{\pi}{2^{r+1}}\right),
```

contrary to the hypothesis. Every r-fold Boolean derivative of f vanishes
at every base address. This is exactly the algebraic-degree condition
$`\deg f\le r-1`$. The proof is independent of width and all source phase
denominators. The radius is a sufficient analytic bound; no optimality of
that radius is asserted.

## A translation echo has Boolean scalar phase

Let $`h:\mathbb F_2^w\to\mathbb R/\mathbb Z`$ be arbitrary and let
$`L,M:\mathbb F_2^a\to\mathbb F_2^w`$ be linear. Write

```math
D_h|z\rangle=e^{2\pi i h(z)}|z\rangle,\qquad
A_x=X^{Lx},\qquad B_x=X^{Mx}.
```

Define the literal commutator

```math
U_x=A_x(D_hB_xD_h^\dagger)A_x^\dagger
                     (D_hB_x^\dagger D_h^\dagger).
```

Direct action on every basis column gives

```math
U_x|z\rangle
=e^{-2\pi i\Delta_{Lx}\Delta_{Mx}h(z)}|z\rangle.
```

The minus sign follows from the actual right-to-left order of the displayed
commutator; the two source calls have positive phases and both actual inverse
calls have negative phases. If $`U_x=\lambda(x)I_W`$, then
$`\Delta_{Lx}\Delta_{Mx}h(z)=c(x)`$ is independent of z.
Translation by Lx negates this difference, since Lx has order at most two.
Thus $`2c(x)=0`$ in $`\mathbb R/\mathbb Z`$, including the case Lx is zero,
and uniquely $`c(x)=f(x)/2`$ with Boolean f. Hence
$`\lambda(x)=(-1)^{f(x)}`$, and $`f(0)=0`$ because $`U_0=I`$.

The echo has four diagonal source occurrences, including both actual
inverses. The general theorem with $`q=r=4`$ now gives $`\deg f\le3`$.
This includes overlapping or dependent translation images, arbitrary
superpositions, and external references. No scalar convention or
intermediate work return is used.

## Robust return at a width-independent error radius

For a prescribed Boolean f, the same echo satisfies

```math
\max_x\|U_x-(-1)^{f(x)}I_W\|_{\mathrm{op}}\le\eta,
\qquad \eta\lt2\sin(\pi/32)
\quad\Longrightarrow\quad f(0)=0,\quad\deg f\le3.
```

The radius is approximately $`0.196034`$. The degree conclusion is the
four-source case of the joint-cube proof, and the value at zero follows
from $`U_0=I`$. This improves the earlier sufficient radius $`2/275`$ and
allows a wider source-word architecture, while preserving the full operator
contract. The exact cubic example below shows that the degree conclusion
cannot be lowered for four sources; it does not assert that the error
radius is sharp.

### Approximate native sources remain charged

For the general word, suppose the actual native occurrence
$`\widetilde D_j`$ approximates $`D_j`$ within $`\delta_j`$ in full operator
norm on all physical work. Any extra source helpers are arbitrary, included
in that norm and in the simultaneous allocation, and have ideal identity
action. Initialized-only accuracy is insufficient. Every actual inverse
occurrence is included separately; the inverse of a source word has the
same full-operator error as the forward word.

Unitary telescoping between the ideal and compiled completed words gives

```math
\|\widetilde W-W\|_{\mathrm{op}}\le\sum_{j=1}^q\delta_j.
```

Thus completed native target error epsilon gives the same degree conclusion
whenever

```math
\varepsilon+\sum_{j=1}^q\delta_j\lt\rho_r.
```

This is a forward error ledger from explicit source accuracy; small error
of a fused completed word alone does not imply accurate individual sources.
All controls, phases, returns and reference systems remain included.

For native occurrence counts $`T_j,G_j`$ and T-depths $`d_j`$, the displayed
sequential source word has the charged upper ledger

```math
T\le\sum_{j=1}^qT_j,\qquad
D_T\le\sum_{j=1}^qd_j,\qquad
G\le\sum_{j=1}^qG_j+G_{\mathrm{intervals}}.
```

Here the monomial Clifford intervals are implemented as exact Clifford
words and $`G_{\mathrm{intervals}}`$ charges every elementary gate in them.
All simultaneously present work wires, including any source helpers, are
arbitrary and counted; there are zero clean helpers. This ledger does not
claim that a general arbitrary-phase source has an exact finite native
implementation. When it requires approximation, the displayed error budget
applies. Sources are sequential occurrences in this architecture; this is
not a lower bound on a different fused implementation.

In particular, if one native echo source $`\widetilde D`$ approximates
$`D_h\otimes I`$ within delta and its actual inverse is used at both dagger
positions, then

```math
\|\widetilde U-U\|_{\mathrm{op}}\le4\delta,
\qquad
\varepsilon+4\delta\lt2\sin(\pi/32)
\quad\Longrightarrow\quad\deg f\le3.
```

If that source has T-count $`T_D`$, T-depth $`d_D`$, total physical
arbitrary-work allocation W, and the actual forward and inverse words have
elementary Clifford counts $`G_D`$ and $`G_{D^\dagger}`$, the literal echo has

```math
T\le4T_D,\qquad D_T\le4d_D,\qquad
G\le2G_D+2G_{D^\dagger}+2\mathrm{nnz}(L)+2\mathrm{nnz}(M).
```

The actual inverse has the same T/T-dagger count and reversed T-layer
schedule, but its elementary Clifford count can differ in the native
alphabet: for example, an S-dagger may require three S gates. Both counts
are therefore charged separately. Each nonzero translation-matrix entry
costs one address-controlled CNOT. Both translation families occur twice. All four source calls, including
actual inverses, use the same W arbitrary wires; the simultaneous allocation
is W arbitrary work and zero clean work, excluding the read-only address.
There is no reset, measurement or uncharged recovery.

## Sharp cubic example, and a literal native ledger

Use six arbitrary work bits $`p_1,p_2,p_3,q_1,q_2,q_3`$, and three read-only
address bits $`x_1,x_2,x_3`$. Take

```math
h_2(p,q)=p_1q_2q_3+p_2q_1q_3+p_3q_1q_2\quad\text{in }\mathbb F_2,
\qquad h=h_2/2,\qquad Lx=(x,0),\qquad Mx=(0,x).
```

Here $`h=h_2/2`$ uses the embedding of $`\mathbb F_2`$ into the two-torsion
subgroup of $`\mathbb R/\mathbb Z`$. With Boolean differences,

```math
\Delta_{Lx}\Delta_{Mx}h_2
 =\sum_{\mathrm{cyc}(i,j,k)}
 x_i(x_jq_k+x_kq_j+x_jx_k)
 =x_1x_2x_3.
```

The sum ranges over the three cyclic permutations of (1,2,3).
Every term linear in q occurs twice. The three cubic terms add to one.
Hence $`U_x`$ is the literal scalar $`(-1)^{x_1x_2x_3}`$ on all six work bits.
The proof checks every physical input, rather than only a prepared code.

For a completely explicit finite native implementation, let $`D_h`$ be the
product of the three CCZ gates on $`(p_i,q_j,q_k)`$. One CCZ can be built from

```math
4abc\equiv a+b+c-(a\oplus b)-(a\oplus c)-(b\oplus c)
                         +(a\oplus b\oplus c)\pmod8.
```

Apply T to a,b,c in parallel. Implement each negative pair parity by a
CNOT into one member, T-dagger there, and the inverse CNOT. Implement the
positive triple parity by two CNOTs into c, T there, and those two CNOTs
in reverse order. This exact word has 7 T/T-dagger gates, 10 elementary
CNOTs, five declared physical T layers, no helper and no residual phase.
It realizes the displayed mod-8 identity using $`T=\mathrm{diag}(1,e^{i\pi/4})`$.

Execute the three CCZ words successively. $`D_h`$ has 21 T gates, 30 elementary
Cliffords and a schedule of at most 15 physical T layers. Use the actual
reversed inverse word for every $`D_h^\dagger`$ call. The displayed commutator
has two $`D_h`$ calls and two actual $`D_h^\dagger`$ calls, plus two A and two B
words of three address-controlled CNOTs each. Its literal ledger is

| Resource | Displayed native word |
|---|---:|
| T/T-dagger count | 84 |
| Elementary Clifford count | 132 CNOTs |
| T-depth | at most 60 |
| Clean work | 0 |
| Simultaneous arbitrary work | 6 bits |
| Error / leakage | 0 / 0 |
| Work return | Exact identity on all six bits |

There are no additional outputs or helpers hidden in that allocation.
All T gates are on work; every address interaction is a CNOT control.
The interlayer address-preservation requirement is therefore literal.
The counts describe this deliberately simple explicit word, not an optimized
compiler. No saving compared with existing native scalar/query baselines is
claimed. Its role is to certify that the cubic upper degree is attained within
the full native/operator contract. The abstract four-representation form
must not be misreported as four physical T layers.

## Scope and precise unresolved step

This closes a broader version of the diagonal-source fusion candidate:
hiding arbitrary table complexity in a fixed number of diagonal work sources
and exposing it by affine dirty-work transport cannot produce Boolean degree
above $`\max\{2,q-1\}`$, even when the source costs are temporarily ignored.
For the double-translation echo, four source occurrences can expose cubic
phase, but cannot expose higher degree below the stated error. Native
synthesis charges can only add to the displayed word's costs.

The theorem does not constrain nonmonomial completed Clifford intervals,
nonlinear work permutations, a source that directly accepts the address, or
arbitrary products of Pauli representations. In particular, it does not
identify q with native T-depth, prove a six-layer query obstruction, or
improve the R91/R96 complete-circuit upper laws. The sufficient radius
$`\rho_r`$ decreases with q, so the robust statement is not a fixed-error
bound for an unbounded number of source calls. A successful constructive
continuation must implement a mechanism outside this affine diagonal-source
architecture while paying the table-loading cost, all sources and actual
inverses, and full dirty-work recovery. The concrete unresolved question is
whether nonmonomial interference or a returned nonlinear intermediate can
remove the R91 address-indicator depth without moving that cost into source
synthesis or recovery.

The [general source-word checks](../../tests/test_diagonal_source_word.py)
verify joint-cube annihilation with distinct rational sources and invertible
affine work mixing, literal Clifford phases including cross terms and
repeated directions, and the operator-return and circular-error thresholds.
The [earlier echo checks](../../tests/test_returned_diagonal_echo_degree.py)
retain the complete native cubic word and actual inverse on all physical
columns, dependent and overlapping translations, non-Boolean rational source
phases, exact full-return elimination for all Boolean sources through four
address bits, and a perturbed source with positive full-input error. These
checks are bounded evidence, not the proof of the theorem.

### Historical four-direction certificate

The earlier proof expanded the fourth echo coefficient in the integer
translation group ring. Its overlapping terms cost at most 120 phase-error
units, seven paired disjoint terms cost at most 20 each, and comparison to
the Boolean fourth difference cost 15, giving $`120+7(20)+15=275`$ and the
older sufficient radius $`2/275`$. Its exact polynomial identities, including
signed additive differences and dependent translations, remain checked in
the earlier test file. A five-bit kernel check there also verifies that
four-direction tests through zero exclude fifth degree. This evidence is
retained for reproducibility; the current theorem and stronger radius use
the single joint-translation proof above.
