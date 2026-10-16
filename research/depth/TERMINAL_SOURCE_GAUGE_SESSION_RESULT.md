# Whole-terminal session: logarithmic general depth

[Current frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth) ·
[Executed design](TERMINAL_SOURCE_GAUGE_PROMPT.md) ·
[R90 proof](DIRTY_CATALOGUE_TERMINAL.md) ·
[Previous R89 outcome](COIN_CONTINUATION_SESSION_RESULT.md)

This session executes the user's “Now go!” on 16 October 2026. Its
baseline is main commit
`53318f59e896a7f3750d35a397709860eaf57522`, containing the design
merged in [PR #153](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/153).
Execution began at **04:44:19 UTC**; the five-hour ceiling is
**09:44:19 UTC**, including review, validation and integration.
The scientific baseline was R89, merged in PR #152.

The session proves [R90](DIRTY_CATALOGUE_TERMINAL.md). For fixed
$`L\ge6`$, $`N=2^n`$ and $`\ell=\log_2(n+2)`$, every prescribed
complete real frame and phase-dressed complex magnitude frame admits

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell).
```

The same circuit retains one external clean qubit, the original real
width $`b\ge17(L+n+7)`$ and complex width $`b\ge17(L+n+8)`$,
and full initialized-isometry error at most $`2^{-L}`$. This includes
all columns, signs and zero coefficients, prescribed complex leaf
phases, actual inverses, arbitrary dirty inputs and external references.
R90 removes R89's unbounded iterated-logarithm multiplier. It does not
prove logarithmic depth optimal. R79 remains the absolute all-precision
theorem.

## 1. The opening ledger located the remaining cost

For R89's $`u(m)=l(m)+\lceil\log_2(m+l(m)+2)\rceil`$, with
$`l(m)=L+\lceil\log_2(m+2)\rceil+10`$, put
$`h(m)=\log_2(m+2)`$. At fixed L there is a constant $`B_L`$ with

```math
2h(m)\le u(m)\le2h(m)+B_L.
```

For $`r(m)=\lceil A u(m)^3\rceil`$, this implies
$`u(r(m))\le6\log_2u(m)+O_L(1)`$. Increasing the fixed stopping
cutoff therefore gives $`u_{i+1}\le u_i/2`$ and
$`\sum_i u_i=O_L(\ell)`$. The same argument makes every fixed-power
local logarithmic cost summable.

The extra factor in the existing circuit came from repeatedly paying
for its long global address. This was an upper-bound charge, not an
unrestricted per-query lower bound. R90 changes the schedule so only a
fixed number of those queries remains.

## 2. What the initial routes taught us

The source-gauge route first tracked two different diagonal Boolean
phase tables through a live Hadamard and an actual inverse source.
Here the masks M_1 and M_2 are diagonal sign involutions. If
$`U=F_2HF_1`$ and
$`\widetilde U=M_2F_2HM_1F_1`$, its output residue is

```math
\widetilde U U^\dagger
=M_2F_2(HM_1H)F_2^\dagger.
```

A diagonal dirty mask generally leaves the diagonal class at the
Hadamard. For a diagonal middle operation R, the complete sandwich
has the exact input-gauge relation

```math
\widetilde U^\dagger R\widetilde U
=M_1(U^\dagger R U)M_1.
```

A common input gauge can telescope across several such sandwiches,
but the actual source families do not automatically have the same
gauge. Coefficient preparation and the accepted-index projector add
another commutation obligation. A full-prefix one-hot bank would also
exceed the available width; balanced readouts retain their own
prefix-dependent correction. This route was not promoted to a native
compiler.

The independent compact-coin completion route reconstructed the
one-quarter/three-quarter alias example: in its eight-dimensional
padded atom space, the two positive fibers have sizes one and five.
Four orthogonal dirty-label images cannot fit in the first fiber even
if both padding atoms are added. Weight-class-preserving gauges also
fail because that fiber's maximum-weight class has only two atoms.
These statements concern that embedding and those gauge classes.
They are not obstructions to a general returned-dirty construction.

Directly reordering the recursive early factors also fails: later
rotations change the earlier suffix-zero predicates. For two bits,
the orders differ in a root-to-$`|11\rangle`$ amplitude. An abstract
direct sum of the resulting partial isometries still needs a native
implementation. Compact node-space and chain-profile representations
were explored but were not used as a proof premise.

These diagnoses suggested changing what must be returned. The useful
unit of cancellation became a complete selected unitary, rather than
an independently corrected phase table.

## 3. The successful change: cancel two dirty residues against themselves

A selected catalogue factor can be transferred onto the logical data
by swapping that data into its selected dirty slot, applying a fixed
catalogue circuit, swapping back, and using the actual catalogue
inverse. This leaves the wanted factor on the data and its inverse on
the selected dirty slot.

Use a separate array for a second factor. Every special unitary is a
commutator $`ABA^\dagger B^\dagger`$, and the unwanted A and B
residues now act on disjoint arrays. In the complete commutator each
dirty residue cancels with its own inverse. Every dirty input returns,
including coherent and reference-correlated inputs. The full operator
identity and a constructive factorization are proved in
[R90 Sections 1–2](DIRTY_CATALOGUE_TERMINAL.md#1-a-selected-small-unitary-and-two-borrowed-arrays).

Selection also has a complete dirty-program word. An XOR lookup
loads $`f(P)`$ into an arbitrary program R. Translate **all** array
slots by that current XOR label, reverse the lookup, and reverse the
translation using the restored R. The two translations commute, so
the unknown initial program cancels and the selected translation
remains. This is an explicitly priced translation network, not a
clean program or a free selected swap.

The approximate catalogue uses its actual native inverse words.
Its dirty arrays still return exactly. Only the selected logical
commutator is approximate, so the error takes a maximum over the
catalogue rather than a sum over all entries. This avoids requiring
catalogue-size precision.

## 4. Why the catalogue fits the original budgets

An entire height-m frame can be programmed with $`q=O_L(2^m)`$
bits by assigning finer angle precision near its root and coarser
precision near its many leaves. With $`2^{n-m}`$ prefix rows, the
whole program table has $`O_L(N)`$ bits. The existing rectangular
multibit query therefore keeps the square-root T-count and linear
Clifford count. Separate one-bit queries would not give that ledger.

At most three R89 stages leave
$`m=O_L((\log\log\ell)^3)`$, or a fixed smaller height. Then
$`q=\ell^{o(1)}`$. Each fixed catalogue factor has an elementary
ancilla-free synthesis of size $`2^{O_L(m)}`$, and all slots run
in parallel. Both arrays together and the total local catalogue gates
have size $`n^{o(1)}`$. They fit simultaneously before reserving
half the dirty pool for lookup in the wide branch.

The proof explicitly retains the original-accuracy complete real or
complex fallback on narrow widths and finite dimensions. Only the
new wide branch appends its separately budgeted complex leaf diagonal.
No stricter tolerance is inserted into a fallback's literal width
formula. Approximate early-source return is handled by unitary
telescoping; the terminal itself needs no clean input.

## 5. Verification and integration record

Independent reviews checked the operator identity, determinant and
literal phase, native no-clean factor synthesis, complete multibit
lookup, geometric program precision, three-stage cutoff, simultaneous
allocation and full-isometry error splice. A fresh full-proof review
then reconstructed the identities and the cited local-stage, query and
fallback interfaces independently. No mathematical blocker was found.
These internal reviews do not replace external technical review.

The [six new bounded tests](../../tests/test_dirty_catalogue_terminal.py)
pass. They include all 128 columns of a noncommuting catalogue word
with native H/S/T factors and exact routing permutations, a dirty Bell
pair with an explicit reference, nonlinear
two-bit program selection on all 512 basis inputs, a single-array
negative control, complete four-column real-frame factors, and the
literal shared-control phase polynomial.

The complete deterministic suite passed **984 tests in 426.242 seconds**
with Python 3.12.14 and NumPy 2.3.5. All 56 focused documentation,
mathematics, presentation-policy and catalogue checks passed. Python
compilation, the reviewer walkthrough, both retained native coarse
examples, the residual-QBP example at q=16, and offline synchronization
checks passed. Both resource ledgers completed at n=12; their retained
audit-status fields are unchanged. All four immutable fault-tolerant
receipt suites match without regeneration.

The final document gate covers all **163 Markdown pages and five
diagrams**, MathJax SVG and native MathML, and desktop and narrow
layouts. After a single-browser memory interruption, the same checks
ran in fresh browser batches. One long inline formula in R90 overflowed
the narrow layout and was moved unchanged into a display equation;
the corrected proof and final session record were rendered again.

Integration is gated on the reviewed tree and the required Python 3.11,
Python 3.13, and rendered-presentation CI jobs. The pull request records
the exact tree, gate receipts and merge commit. The proof home is R90;
bounded fixtures do not certify its asymptotic estimates. No final
manuscript drafting or scalable circuit emitter is part of this session.

## 6. The next mathematical questions

The fixed-accuracy logarithmic upper target is reached. The general
large-width optimum is still between the retained constant lower bound
and this logarithmic upper bound. A dimension-growing unrestricted lower
bound, or a still shallower general construction, would decide more.
A precision-uniform version of the new gain is also open: the catalogue
and stopping cutoffs here are allowed to depend on L.

The continuing destination is optimal depth under the same complete
operator and workspace contract. This session changes the proved
frontier; it does not change that destination.
