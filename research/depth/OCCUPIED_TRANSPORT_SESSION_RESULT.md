# Revision: logical transport, native algebra costs, and the remaining depth question

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Occupied transport](OCCUPIED_LOGICAL_TRANSPORT.md) ·
[Returned algebra queries](ALGEBRA_QUERY_EXTRACTION.md) ·
[Previous checkpoint](CLEAN_CODE_ALGEBRA_SESSION_RESULT.md)

This revision began from verified main `b84b072` (PR #186), on
19 October 2026 in Asia/Shanghai. It replaces two unsupported shortcuts
with exact native reductions and complete resource accounts. The general
complete-frame frontier remains unchanged.

## 1. What is the destination, and what is already proved?

The destination is optimal general T-depth for the prescribed complete
Hopf differential frame. Every logical column, literal phase and sign,
zero and singular coefficient, dirty input and reference, actual inverse,
work return and leakage remains in the contract. All resources must be
charged on the same circuit.

| Role | Current statement |
| --- | --- |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, for all integer $`n\ge1,L\ge6`$ |
| Literal allocation | One clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error $`2^{-L}`$ |
| Selected claims | R79 remains Result D; the count package is settled under its stated hypotheses; manuscript drafting is separate |
| Unrestricted component lower results | At most one clean and two layers imply XOR degree at most two; zero clean and three layers imply degree at most four, with retained uniform error floors |
| Useful constructive components | Optimal-depth quadratic queries; complete balanced HRT recovery; conditional R98 transfers a faster arbitrary-table word query satisfying all its budgets |
| Scoped limitations | R97's overlap hypotheses, failed clean-flip or Clifford-extension arguments, and arithmetic-stage degree bounds do not exclude the unrestricted model |
| Open destination | Optimal general depth and dimension-growing matching lower bounds; exact one-clean three-layer degree and unrestricted one-clean four-layer rigidity remain unresolved |

At fixed accuracy and square-root dirty width, the general upper depth
is logarithmic in n and the available unrestricted lower bound remains
constant. A new component does not close that gap merely by having a
compact algebraic description or a better isolated gate count.

## 2. Lower route: use full Pauli transport instead of extending a code

The former first unsupported inference was that a two-band logical image
agreeing with a Clifford on the occupied code must itself be Clifford.
The new [complete-family witness](OCCUPIED_LOGICAL_TRANSPORT.md#5-complete-family-extension-witnesses)
refutes this for arbitrarily large clean-active rank, for every logical
generator simultaneously, with no middle pairs and no shared Pauli clean
flip. Its occupied target is identity and its full word has a paid
one-layer simplification. It is a failure of that proof inference,
not an unrestricted target obstruction.

The replacement is an actual identity. Every clean-swapped middle pair
has a difference Pauli commuting all three bands. It therefore transports
one physical Pauli to another on both clean sectors. For a Boolean sign
g, the active X masks of these transports span a space U on which all
Boolean derivatives are affine. In common linear coordinates,

```math
g=q+h,\qquad \deg q\le2,\qquad
\deg g\le\max\{2,n-\dim U\}.
```

Paid endpoint Cliffords remove q, align all transports, and diagonalize
their dirty parts. The same width, T-count and T-depth then suffice for
a representation in which every pair difference is a globally preserved
computational-Z parity. Forward and actual inverse Clifford overheads
are each $`O(Q^2)`$ at total physical width Q. This is a necessary
structural reduction; its overhead is not being inserted into the
optimal-count frame compiler.

For full-degree query predicates, including AND, U is already zero:
the reduction preserves the original query sign without mixing its active
coordinates. A further pair-free subcase closes to occupied C3 when the
clean-silent endpoint factors are occupied-diagonal. That additional
hypothesis matters: the retained three-layer quartic sign is pair-free
but not C3.

The first remaining unsupported step is now explicit. Sectorwise middle
Cliffords cannot be removed and followed by a one-relation theorem for
free: conditioning the conserved parities creates many projected-axis
relations, and the joint phase dependence still needs control. Pair-free
cases also retain clean-silent bands and non-Pauli unused-sector freedom.
Neither a growing relation count nor a failed extension proves hardness.

The [zero-clean phase proof](THREE_LAYER_NATIVE_STRUCTURE.md#6-exact-pair-phases-and-the-common-coordinate-bound)
was also sharpened: its special affine-parity times Boolean-quadratic
pair phases have integer degree at most five modulo eight in the original
coordinates. Arbitrary degree-five integer phase polynomials do not have
that affine invariance. The existing zero-clean XOR degree-four theorem
is unchanged.

## 3. Constructive route: pay the table extraction and change multiplication

The old HRT coordinate query could not treat nonsingleton coefficients as
free zeros. The new [two-commutator extraction](ALGEBRA_QUERY_EXTRACTION.md#1-an-arbitrary-table-scalar-on-the-entire-coefficient-space)
returns every coefficient input while producing any prescribed Boolean
table. It uses eight actual source or inverse occurrences; no eigenstate,
clean coefficient, measurement or postselection is supplied.

The elementary product count $`3^r`$ was then replaced by the established
ranked subset-convolution factorization, compiled through arbitrary dirty
linear temporaries and an eight-corner native CCZ phase cancellation.
All copies, transforms, inverse words, masks and output Hadamards are
charged. This gives a second complete balanced source with

```math
T_F,G_F=O(r^2 2^r\log r),\qquad
B_F=O(r^2 2^r),\qquad D_F=O(\log r),
```

zero clean work, exact return, zero leakage and the same bounds for its
actual inverse. The earlier $`O(3^r)`$ source retains its smaller depth
constant and separate literal ledger; neither variant is asserted to
dominate on small instances.

This does not supply R98. On an a-bit table, the extraction has
$`Q=2^a`$ rows but count $`O(a^2Q\log(a+2))`$, dirty width
$`O(a^2Q)`$, and logarithmic depth. R98 needs
$`T=O(\sqrt{Qq}+Qq/B)`$, $`G=O(Qq)`$, and
$`D_T=O(d(a)+Qq/B^2)`$ for all its word lengths and allocations;
sublogarithmic d would improve the fixed-accuracy frame frontier.
The initialization gap is closed, but those simultaneous costs remain.

Two exact diagnostics explain which further shortcuts fail. Deferring
arbitrary helper garbage cannot remove the inverse degree of an HRT
permutation, so bounded-degree completed classical stages still need
logarithmically many stages. A full-degree table functional has a
nondegenerate product pairing on the coefficient algebra, so quotienting
that algebra cannot shrink its one-qubit-per-coefficient carrier to
square-root width. Both statements are representation-specific.
Neither bounds unrestricted native T-depth or alternative encodings.

## 4. Continuing question and evidence

The next lower question is the completed Boolean phase in the paid
diagonal pair-center normal form, including its joint control dependence,
silent endpoint factors and unused-sector freedom. For a constructive
escape, a mechanism must improve the full arbitrary-table ledger, or
compile the complete frame through another paid route; another balancing
of bounded-degree classical recovery stages retains the old depth order.

The [transport fixtures](../../tests/test_occupied_logical_transport.py),
[algebra-query fixtures](../../tests/test_algebra_query_extraction.py),
and [three-layer phase fixtures](../../tests/test_three_layer_native_mechanisms.py)
retain native physical words, literal phases, all relevant small input
columns, arbitrary work, full unused blocks, and actual inverse costs.
Independent proof review covers the normal form, native fast transforms,
source allocation, extraction and the scopes of both limitations.
Finite checks support these identities; the arbitrary-width statements
follow from their analytic proofs. R96 and the general lower frontier
have not been promoted.
