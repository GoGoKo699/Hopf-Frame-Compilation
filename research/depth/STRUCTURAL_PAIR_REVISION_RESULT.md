# Revision: remove the matrix core before counting its features

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous coherent-source checkpoint](COHERENT_SOURCE_AND_ONE_PAIR_RESULT.md) ·
[New proof](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)

This revision starts from main
`a513c534a81d9c257b5cfa4efc1da83db4affc5b`, tree
`9f1f6a23e72c8c64cfff3cf1b1a04c322a8f25c1`, after PR #208.
The earlier sign/XOR degree bound 37,748,865/37,748,864 was a generic
feature-count certificate, not a physical degree scale. The actual
one-pair structure supports a substantially sharper conclusion.

## 1. Destination, established results, components and limits

| Role | Revised assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth with simultaneous T-count, Clifford count, dirty width, complete operator error and return |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md), with selected [R79](UNIFORM_SUBTREE_DEPTH.md) and all literal allocations unchanged |
| Strongest general lower components | [Two-layer Boolean rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md), [zero-clean three-layer rigidity](THREE_LAYER_NATIVE_STRUCTURE.md), and the unrestricted indicator T-count bound remain distinct from a growing-depth theorem |
| Revised focused lower theorem | [Fixed full-space diagonal reduction](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md): in the all-active genuine original one-pair class, sign/XOR degree is at most **6/5 uniformly**, and **5/4 when the common active rank is at least four** |
| New physical recognition | [Positive-singleton recognition](POSITIVE_SINGLETON_PAIR_RECOGNITION.md): with exactly one genuine original pair and no negative singletons after literal orientation repair, complete Booleanity forces the all-active support condition. The same small degree bounds therefore apply without assuming that support |
| Retained component | The old proof's scalar-branch matching and one-X-mask conclusion remain essential. Its generic matrix-feature estimate is no longer needed for this class |
| Upper component and missing operation | Four-call extraction still permits arbitrary quantum private-helper fibers. A complete sublogarithmic native weak source or whole extraction remains unpaid |
| Scope of the obstruction | The direct degree proof assumes all-active support; the recognition theorem derives that support in the positive-singleton one-pair class. Negative singletons, arbitrary pair count and unrestricted one-clean three-layer circuits remain outside |

The revised lower theorem improves the known degree frontier inside its
specified physical class. It does not improve the full upper resource
frontier or the unrestricted logarithmic-versus-constant depth interval.
The huge old number came from

```math
37748865=48\cdot6\cdot2^{17}+4\cdot23+1+36.
```

The first term enumerated possible one-qubit matrix types, followed by
arbitrary Boolean decoding of their counts. That was a valid but very
lossy last step. A larger catalogue or a smaller enumeration constant
would leave the same structural information unused.

## 2. The replacement mechanism

The [proof home](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md) retains the
literal all-active reduction and then changes its final mechanism.

1. Scalar-branch native matching gives a dressed Hermitian Clifford
   with one computational X mask. The exact complete-target identity
   transports the **same** Pauli to the input side, giving the same mask
   there as well.
2. Independence of each actual native endpoint list forces every axis
   that dresses this Pauli to pull back to a computational Z. Peel these
   gates as literal diagonal C3 factors. Every retained band is an actual
   sublist, with no projected-list resynthesis.
3. A fixed diagonal Clifford aligns the remaining genuine endpoint
   Paulis. A fixed monomial Clifford makes their common axis X on one
   coordinate. The diagonal target commutes this X, so its phase is
   independent of that coordinate.
4. A Hadamard on that coordinate leaves the entire target unchanged and
   makes the pair insertion diagonal. The residual actual three-native-
   band word is therefore diagonal on its **full** logical input space.
   The sharp zero-clean degree-five theorem applies directly.

There is no restriction to another Pauli eigensector. In particular,
the reduction does not trade physical native independence for a
dependent projected phase list, and it does not add a clean promise.
The literal insertion phase has ordinary degree at most five modulo
eight. Restore the peeled phases in the same fixed coordinates and only
then use affine invariance of Boolean degree. This yields sign/XOR
degree five/four for active rank at least four. The retained rank-two
theorem supplies six/five for the remaining case.

The proof retains all occupied logical columns, arbitrary dirty and
reference inputs, exact clean return, the common scalar, and an
unrestricted unused clean block. It introduces no approximation error
or initialized helper. It is a necessary condition on an existing
circuit, not an uncharged new upper compiler. The conservative
original-circuit rounding floor remains $`2^{-530690}`$ under the
same actual original-list hypotheses.

## 3. What remains unpaid on the upper route

The [four-call source recurrence](OUTPUT_ALIASED_CHARACTER_SOURCE.md)
is unchanged. For $`m=2^b`$ and a weak quantum-fiber source V with
full-operator error eta, direct readout gives

```math
\begin{aligned}
T_I&\le4T_V,&D_T(I)&\le4d_V,\\
G_I&\le4G_V+bm+2b+2m+6,&h_I&\le m+h_V,\\
\varepsilon_I&\le4\eta,&a_{\rm clean}&=0.
\end{aligned}
```

Actual adjoints cancel the ideal private-helper fibers. The displayed
complete error includes return on arbitrary inputs and references.
Root fusion uses 32 source calls and error at most $`32\eta`$;
request $`\eta=\tau/32`$. Fixed-L frame transfer still requests
$`\eta_L=2^{-L}/(512K_*)`$ with the retained literal allocations.
The known source retains logarithmic depth.

The distinct amplitude-$`1/\sqrt2`$ interface still has depth
$`24d_U+32`$ and complete error
$`16(3\sqrt2\epsilon+2\epsilon^2)`$. Signal and toggle return
are approximate unless exact flatness is separately proved. The
bounded-clean two-band relative-gap theorem is a source-specific
obstruction; it does not supply an unrestricted native-depth bound.
Changing the lower proof does not pay a new upper source.

## 4. The next scout changes mechanism again

The first attempted proper-support extension did retain an unpaid step.
Writing its insertion as $`G_0=D F_h`$ would not supply the same-Pauli
transport when D does not commute P_h. Likewise, one X-mask coset modulo
an extra direction does not make the dressing gates diagonal. A
Hadamard on that extra direction need not preserve the complete target.
Another finite matrix catalogue would hide this issue.

The [recognition proof](POSITIVE_SINGLETON_PAIR_RECOGNITION.md) instead
returns to the actual original word. With no negative singletons, the
endpoint clean images agree and give a common active span. Proper
support would supply an actual Pauli t commuting with both endpoint
bands and every singleton, but anticommuting with the paired axis.
On the literal scalar pair branch its transported image equals the
output t Pauli; on the quarter-turn branch those two operators
anticommute.

The shared computational pair parity and scalar branch force the two
endpoint t Paulis to have the same X mask. A Boolean diagonal target's
phase derivative is invariant under translation by that mask, so its
transported Pauli must have one global commutation sign with the output
t Pauli. The two nonempty parity branches contradict this. Proper
support is therefore impossible in the positive-singleton one-pair
class: complete Booleanity forces all-active support at both actual
endpoints. Combining with Section 2 proves the same uniform six/five
bound throughout that class. No two-core feature estimate is needed.

The next unresolved original-one-pair question requires **negative
singletons**. They can give different endpoint clean images, so the
common-active-span inference used here is unavailable. A concrete first
question is whether negative singletons equal to actual endpoint axes
can be combined literally with their neighboring T phases, and what
remains when they are only products of endpoint axes. Product-axis
rewriting must not introduce dependent native lists. This is a research
question, not an additional closure theorem. General pair count,
unrestricted one-clean three-layer degree and dimension-growing depth
remain open.

The revision is analytic. The retained
[all-active fixture](../../tests/test_all_active_singleton_closure.py)
checks the central-controlled quarter-turn, literal scalar branch,
signed conventions, all physical columns and actual inverse. It does
not assert the new degree theorem by enumeration. Independent audits
checked both-end native peeling, full-space Pauli alignment, the phase
degree and restoration, and the positive-singleton Boolean commutation
contradiction. These internal checks do not replace external technical
review.
