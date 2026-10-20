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
| Retained component | The old proof's scalar-branch matching and one-X-mask conclusion remain essential. Its generic matrix-feature estimate is no longer needed for this class |
| Upper component and missing operation | Four-call extraction still permits arbitrary quantum private-helper fibers. A complete sublogarithmic native weak source or whole extraction remains unpaid |
| Scope of the obstruction | The new degree bound assumes the actual original pair difference sums every clean-active axis at both endpoints. It does not cover arbitrary original pairs or unrestricted one-clean three-layer circuits |

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

## 4. Next research question and evidence

Continue with an actual original pair whose difference uses a proper
subset of the endpoint-active axes. Determine whether the same complete
Boolean condition forces a removable Pauli symmetry, or whether an
additional retained transition prevents the fixed full-space
diagonalization. The first task is an original-list operator identity
and native matching condition, not another finite matrix catalogue.
Negative original singletons may give different endpoint clean images;
the all-active common-span reduction cannot simply be reused there.
General pair count, unrestricted one-clean three-layer degree and
dimension-growing depth remain open.

The next proposed proper-support extension has a precise audit point.
If its insertion is written $`G_0=D F_h`$, with F_h commuting P_h
but D not known to commute P_h, the exact transport is dressed by D;
the same-Pauli identity used above is unavailable. Likewise, one X-mask
coset modulo an extra direction does not force the dressing gates to
be diagonal. A Hadamard on that extra direction need not preserve the
complete target. The missing inference is a full-space Pauli symmetry
or an exact replacement cancellation, not a finite-feature estimate.
A focused next question is one genuine original pair with four-axis
support at both endpoints, at least two additional active axes, and
only positive singletons after literal orientation repair. Does complete
Booleanity force such a reduction? This extension is not a theorem of
the present checkpoint.

The revision is analytic. The retained
[all-active fixture](../../tests/test_all_active_singleton_closure.py)
checks the central-controlled quarter-turn, literal scalar branch,
signed conventions, all physical columns and actual inverse. It does
not assert the new degree theorem by enumeration. Independent audits
checked both-end native peeling, full-space Pauli alignment, the phase
degree and restoration. These internal checks do not replace external
technical review.
