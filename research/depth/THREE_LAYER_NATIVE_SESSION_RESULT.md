# Three-layer revision: a native degree invariant and a paid cubic query

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof](THREE_LAYER_NATIVE_STRUCTURE.md) ·
[Projected-band component](ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation) ·
[Previous checkpoint](DIAGONAL_ENCODING_SESSION_RESULT.md)

This continuation began from verified main
`7cae235a80cda459602b537a0d72e1319fc7e9da` (PR #184), following
the 18 October 2026 checkpoint. It addresses the first new native gap
at three layers. The general complete-frame frontier is unchanged.

## 1. Revision and the selected question

The destination remains optimal general T-depth for the prescribed
complete Hopf differential frame, with simultaneous count, elementary
Clifford, precision and workspace budgets. Every column, literal phase,
singular and zero coefficient, arbitrary dirty input and reference,
actual inverse, leakage and complete work return is part of the contract.
Every source, control, query, recovery and inverse is charged.

| Category | Current statement |
| --- | --- |
| Strongest unconditional general upper law | R96: for all integers n at least one and L at least six, depth $`O(NL/b^2+L+\log_2(n+2))`$, T-count $`O(\sqrt{NL}+L+NL/b)`$ and elementary Clifford count $`O(NL)`$, on the same circuit |
| Literal allocation and precision | One clean qubit; dirty width $`b\ge17(L+n+7)`$ for real frames and $`b\ge17(L+n+8)`$ for complex magnitude frames; full error at most $`2^{-L}`$ |
| Selected publication claims | Count package settled under its hypotheses; R79 remains selected Result D; final manuscript drafting is separate |
| Useful lower components | R97 bounds outer overlap; R75 gives full-space Pauli rounding; the unrestricted at-most-one-clean two-layer Boolean degree bound is now proved |
| Useful constructive interface | R98 accepts exact completed address return, including internal address mixing, but still needs a faster arbitrary-table native query within every simultaneous budget |
| Unresolved destination | Improve the additive large-width depth while keeping the full budgets, or prove an unrestricted matching lower bound; dimension-growing lower bounds and unrestricted one-clean four-layer rigidity remain open |

At fixed accuracy and square-root dirty width, R96 has logarithmic
additive depth in n. A constant-layer query obstruction does not explain
whether this growth is necessary. The question selected here was whether
the joint native invariant survives the first asymmetric split: one
band on one side and two on the other. Arbitrary C3 source structure
cannot be inferred from a preserved maximal Pauli algebra, and the
previous one-band clean-flip argument cannot simply be iterated.

Both directions were pursued. A lower invariant valid for arbitrary
width could become a base for a depth-dependent bound. A construction
that pays the previously missing phase could expose a mechanism that
such a bound must admit. Neither route is judged by an abstract identity
alone.

## 2. What changed in the lower analysis

The [new full-space theorem](THREE_LAYER_NATIVE_STRUCTURE.md#1-exact-statement-and-conventions)
proves that a zero-clean native three-layer diagonal circuit has an
ordinary mod-eight phase polynomial of degree at most five in one
common affine coordinate system. Therefore literal Boolean signs have
degree at most five and literal XOR predicates degree at most four.
Arbitrary intermediate Clifford mixing and arbitrary returned dirty
width are allowed. R75 rounding and native root discreteness give the
same uniform error floor

```math
\tau=\frac{\sin(\pi/8)}{2\sqrt2}
```

for any XOR predicate of degree greater than four at this depth, even
after optimizing one common phase. Thus integer precision exponents
at least three require at least four T layers for those zero-clean
queries. This is an unrestricted zero-clean statement, not a read-only
architecture bound or a one-clean result.

The key step retains the entire preserved Pauli algebra. The middle
physical axes have only singleton and paired permutation orbits. A
joint constraint from the opposite native band correlates each pair's
sign and swap characters. Exact orientation changes then make each pair
phase commute with the whole algebra. After peeling these diagonal
factors, the remaining source lies in one commuting Pauli space.
Its symplectic displacement identifies the opposite band's active axes;
silent axes peel on the output side. This proves the residual C3
statement that was initially missing, instead of assuming it.

Every paired factor restricts a Hermitian Clifford to a parity sector,
giving a quadratic Boolean branch phase. The independent pair differences
can all be made coordinates at once; separately changing coordinates
for each pair would not justify the stated integer-degree bound.

The separate [projected-band theorem](ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation)
extends axis recovery through one Pauli stabilizer quotient. Only a
four-axis relation leaves an ambiguity, confined to one identified
three-dimensional span. Its necessary matching condition is robust
below $`\sqrt2/4`$ on the already identified quotient. It neither
supplies that quotient for a non-Pauli intermediate code nor proves
the one-clean three-layer degree bound.

## 3. The constructive phase is now paid

The earlier two-layer occupied operation was $`(iZ_t)^{q(x)}`$.
Its input-dependent phase prevented a literal higher-degree query.
The first missing operation was the correction of both its imaginary
phase and its dirty-dependent sign. A third independent physical band,
with all parity computations and Clifford corrections included, now
cancels both and yields the exact sign $`(-1)^{x_1x_2x_3x_4}`$.
Output Hadamards give the complete three-control AND query.

The [native construction](THREE_LAYER_NATIVE_STRUCTURE.md#8-paying-the-third-band-an-exact-cubic-xor-query)
uses **18 T gates, three T layers and 99 elementary Cliffords**, with
one clean qubit and one arbitrary dirty helper. Its actual inverse uses
18 T gates, three layers and 109 elementary Cliffords after adjacent-S
reduction. Error and leakage are zero; every logical input and dirty
reference is covered. The unused physical block is not Hermitian, so
repeating the forward word is not its full actual inverse.

A fairly simplified relative-phase compute/uncompute baseline uses
**16 T gates, 34 Cliffords and at most six layers**, with exactly the
same allocation. The new construction improves demonstrated depth,
not count. The previous unrestricted two-layer lower theorem makes
depth three optimal for this finite query and allocation.

The same mechanism gives a growing-arity native three-layer family
$`(-1)^{e_4(x)}`$ with linear counts and the same one-clean/one-dirty
allocation. Its degree remains four. Beyond four inputs, turning one
bit into an XOR output leaves the address sign $`(-1)^{e_4}`$ on the
other bits. This is the next unpaid operation in that attempted
high-arity query conversion. The family is a useful native stress test,
not a faster restored-address arbitrary-table compiler.

## 4. Scope failure, historical comparison and exact evidence

An exact two-qubit example disproves the proposed inference that an
arbitrary one-band/two-band split with a non-Pauli clean image has a
shared Pauli clean flip with Pauli preimages. It has a cheaper complete
one-layer realization, and its occupied target is just T-dagger. The
failure belongs to that representation and proof mechanism, not to the
unrestricted resource model. It cannot be promoted to an obstruction
against alternative native representations.

The known growing-degree C3 permutation family in
[native transport](NATIVE_CLIFFORD_TRANSPORT.md#7-third-level-membership-alone-does-not-bound-tableau-degree)
remains another warning against replacing native structure by abstract
hierarchy membership. Its paid realization has growing depth; it is
not a three-band counterexample. The source certificates and read-only
degree bounds keep their original hypotheses. Historical next-step
passages do not reopen the resolved count package, two-layer bent-query
exclusion or all-precision upper law.

The [native mechanism fixtures](../../tests/test_three_layer_native_mechanisms.py)
check literal occupied columns, three physical bands, inverse matrices,
clean and dirty return, and the simplified baseline. The
[projected-band fixtures](../../tests/test_projected_band_matching.py)
use exact cyclotomic coefficients for dependent Pauli images, Fourier
axis recovery and a native stabilizer lift. Universal degree and error
statements are analytic, not inferred from bounded enumeration.

## 5. Precise resumption question

**Can the joint native constraints on a non-Pauli intermediate clean
code bound the completed Boolean phase degree at three layers,
independently of address arity and dirty width?** The invariant must
use the entire logical Pauli algebra and tolerate the absence of a
shared Pauli clean flip. The projected-band result is available whenever
a genuine Pauli-code restriction can be justified; it is not a global
premise. Alternatively, a full native unbounded-degree query family
would resolve this question constructively, but its address phase and
actual recovery must be paid.

This question tests whether the new pair decomposition can support a
depth-dependent invariant rather than a sequence of isolated predicate
exclusions. A sharper zero-clean degree bound is secondary: abstract
quadratic branch phases alone can already produce quintic signs, so
sharpening needs further joint native compatibility. The parallel upper
target remains a faster restored-address arbitrary-table query within
R98's simultaneous counts, allocation, inverse and return budgets.

No complete-frame resource bound improves in this round. The two-call
query/frame reduction loses a factor two in depth, so the new
constant-layer query obstruction does not strengthen the retained frame
lower frontier. R96, selected R79 and scoped R97 remain unchanged. The
general objective remains active; this checkpoint records a proved
zero-clean step, an optimal finite query and the precise one-clean gap.
