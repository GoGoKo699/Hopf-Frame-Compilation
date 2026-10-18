# Four-layer query closure: follow-up outcome

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof](READ_ONLY_THREE_LAYER_RIGIDITY.md#8-four-layers-with-zero-clean-qubits-degree-at-most-two) ·
[Three-layer round](THREE_LAYER_QUERY_SESSION_RESULT.md)

This follow-up began from verified main
`3b7c603a6688417e48b349dbed8facbc91ae07a0`, the merge of PR #167.
It tests the first case of that checkpoint's degree-growth question.
The original read-only interface and full arbitrary-work operator
contract are retained. No complete-frame upper bound changes.

## 1. The obstruction that prompted reformulation

The preceding round proved robust degree at most four with one clean
qubit and three physical T layers, and affine XOR output with zero
clean qubits at three layers. A phase-sensitive exact argument also
excluded nonaffine one-clean two-layer queries. Native AND fixtures
attain four zero-clean layers and three one-clean layers, with all
phases, inverses and work return explicit.

At four layers, the direct normal form has three dressed Pauli
representations. The product of the last two need not be a projective
representation. Thus the old cocycle cannot simply be iterated.
This follow-up changes the target representation: conjugating away the
outer physical layer leaves two projective representations, while the
comparison target becomes a **Clifford**, rather than a Pauli. This
preserves the actual-inverse estimate and makes a new cocycle available.

## 2. Supported result

Section 8 of the [same proof home](READ_ONLY_THREE_LAYER_RIGIDITY.md)
proves that a zero-clean, four-layer read-only XOR query satisfying

```math
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}
\le\eta\lt\frac{\sqrt2}{12}
```

has $`\deg_{\mathbb F_2}f\le2`$, independently of address length,
finite dirty width, T-count and elementary Clifford count. The target
and every dirty helper are arbitrary, including correlations with
external references. Every completed Clifford interlayer preserves the
address pointwise and physical T gates act only on work. Approximate
return is included in the full operator norm. No initialized helper
or intermediate workspace return is assumed.

A literal zero-clean Boolean scalar phase at four layers has degree
at most three, by the charged two-Hadamard derivative conversion.
This is a necessary condition; scalar attainability is not claimed.
The XOR degree-two bound is sharp, since the retained native AND query
has degree two and four layers. Its exact ledger remains eight T gates,
18 elementary Cliffords, zero clean qubits, one additional dirty helper
and zero error. The target is separate from that helper allocation.

The proof first treats the case where the target bit changes the outer
layer's symplectic action. Comparing highest-degree coefficients in
the Clifford cocycle makes every Boolean derivative affine.
If that symplectic contribution vanishes, a second invariant is needed;
the vanished parameter is not itself a lower-bound proof.

In that case the baseline Clifford family factors into a linear Pauli
family and a one-T-conjugated Pauli representation. The auxiliary
commutator with a fixed Clifford has a Pauli-valued cocycle. Independence
of the **physical** T axes forces its label to be linear. The remaining
defect has degree at most two, initially giving degree at most three
for the query. Its quadratic coefficients lie in the image of a native
independent-axis shear. A subspace of directions of codimension at most
one avoids the target label in that image, forcing affine derivatives
there and excluding the remaining cubic case.

Distinct projective Clifford actions have an operator separation at
least $`1/\sqrt2`$. Together with the retained $`6\eta`$
comparison, this gives the stated width-independent approximation
radius. Pauli-label comparisons also fit that error budget. The proof
retains scalar phases where they affect the target and uses projective
equalities only for necessary label invariants.

## 3. What is unchanged and what to ask next

The count package remains settled under its stated hypotheses. For
$`n\ge1,L\ge6`$, R96 still gives the strongest unconditional
general upper law

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)),
```

with one clean qubit, error $`2^{-L}`$, and literal dirty thresholds
$`17(L+n+7)`$ for real or $`17(L+n+8)`$ for complex magnitude
frames. Every logical column, literal sign and phase, singular and
zero coefficient, dirty/reference input, source, query, control,
recovery, actual inverse, leakage and work return remains charged.
R79 is still selected Result D. R97's four-layer obstruction keeps its
outer-misalignment hypotheses; unrestricted one-clean four-layer
rigidity remains open. The present zero-clean query theorem does not
settle that problem.

The most useful next question remains a **uniform degree-growth law**
for exact zero-clean read-only queries: can one prove
$`\deg f\le2^{Cd}`$ with an absolute C for d physical T layers,
or construct a literal native family escaping all width-independent
degree bounds at some fixed d? Intermediate growth rates are possible.
The four-layer test is now resolved; the first unclosed test is five.

The concrete missing step has changed. The four-layer proof handles
Clifford comparison targets by separating symplectic and Pauli defects.
At five layers, peeling another physical layer can produce a
non-Clifford comparison target. Neither the finite Clifford separation
nor the auxiliary Pauli cocycle then follows. A general argument must
control this new conjugation without assuming closure of the Clifford
hierarchy under products, discarding literal phases, or assigning
independent physical axes to a dependent transported representation.
This is a resumption question, not an asserted recurrence.

An exponential degree law would force $`\Omega(\log a)`$ depth
for an a-address AND query, obstructing a sublogarithmic improvement
through R98's existing exact zero-clean interface. A constructive escape
must additionally meet R91's multi-output T-count, elementary Clifford
and simultaneous-width envelope before it improves the complete circuit.
The bilinear echo already has constant depth; its address-indicator
cost must actually fall.

For a general frame lower bound, fixed-error robustness across growing
depth and unrestricted native address mixing remain separate obligations.
The retained unrestricted two-layer dirty Toffoli still falsifies a
naive enlargement of the present hypotheses. The two-call query/frame
reduction accepts unrestricted queries, while conditional R98 uses its
distinct read-only upper interface. Optimal general depth remains open.

## 4. Evidence and pause

The three bounded [four-layer checks](../../tests/test_read_only_four_layer_rigidity.py)
support the literal native normalization, actual inverse, auxiliary
defect and image identities, and codimension-one derivative step.
The degree and approximation bounds are analytic and were independently
reviewed, including the vanished-symplectic case and the final
codimension-one argument. The previous round's complete native fixtures,
resource ledgers and scope falsifiers remain in force. Required
validation and presentation results are recorded with the follow-up merge.

This pause records a closed four-layer case and a precise unclosed
five-layer mechanism. It does not turn a constant-layer theorem into a
dimension-growing lower bound or retire the general-depth objective.
