# Revision: parallel access is paid; selective dynamics remain missing

[Workspace](../../WORKSPACE.md) ·
[Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous continuation](JOINT_SELECTIVE_CONTINUATION_RESULT.md)

This continuation starts from main
`61ce55be85c2f167b2be32fcde5ce946cc08961d`, tree
`17382985677e43cfd35e4e110763e52db128dc98`, after PR #204,
on 20 October 2026. It revises the candidate before extending it and
keeps the focused original physical-pair lower scout active. The new
results pay parallel affine access, delimit its completed Pauli echo,
repair a multioutput proof gap, and identify the exact cost of a larger
physical-pair sector representation. They do not improve the unrestricted
complete-frame frontier.

## 1. Destination, theorems, components and approach limits

| Role | Current assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth with simultaneous count, dirty width, complete error and return bounds on the same circuit |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md) remains the all-precision theorem; selected R79 remains Result D. At fixed accuracy and square-root dirty width, the logarithmic-versus-constant depth gap remains open |
| Unrestricted lower information | The retained Boolean layer obstructions and [linear indicator T-count bound](INDICATOR_PAULI_NULLITY.md); count alone does not force growing depth |
| Useful constructive reductions | [Root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) absorbs polynomial source overhead; [clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) pays the bounded-clean or amplitude-source wrapper; exact joint detuning has an 11-layer linear-count implementation |
| New paid component | [Parallel affine detuning](PARALLEL_AFFINE_DETUNING.md) implements polynomially many affine presentations in four T layers, at polynomial address overhead in T-count and arbitrary dirty width |
| New approach limits | One commuting batch has error at least one quarter. A completed Pauli echo with read-only address/bank and actual Clifford interlayers needs at least the ceiling of b/2 detuning rounds. Neither is an unrestricted native-depth lower bound |
| New lower-track components | The word-valued four-layer read-only degree-two theorem is proved. Larger physical-pair sector codes can need two endpoint bands; paying them does not recover a three-band comparison |

For clarity, R96 retains

```math
D_T=O(NL/B^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
```

Here $`N=2^n`$, $`n\ge1,L\ge6`$. The real/complex dirty
reservations remain respectively $`B\ge17(L+n+7)`$ and
$`B\ge17(L+n+8)`$, with one clean qubit and complete
initialized-isometry error at most $`2^{-L}`$. Discovery retains
arbitrary address, output, dirty-work and reference inputs, literal phase,
actual inverses, leakage and return.

The original overlapping two-band route at amplitude $`1/\sqrt2`$
remains stopped at the every-fixed-relative-accuracy scope for uniformly
bounded clean count. The [amplitude-free source defect](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
is uniform over positive amplitudes and dirty width. Additional source
bands remove a hypothesis; they do not provide a transition or improve
the complete resource frontier.

## 2. First unpaid operation after the scheduling repair

Let $`m=2^b`$. For q affine presentations $`z_j=L_jx+t_j`$,
fixed row permutations and signed forward/inverse calls, the commuting
full-detuning product is now paid by the retained four-layer dirty phase
relocation identity. There is no initialized affine address copy. The
proof includes the literal scalar, all inputs, actual inverses and exact
return of added helpers. For polynomial q its costs are

```math
D_T\le4,\qquad T,b_{\rm dirty}=O(qmb),\qquad G=O(qmb^2).
```

This fits root fusion's relaxed source budgets. It removes an access
scheduling question but does not return the original signal banks.
The first unpaid operation is still the complete selective identity

```math
\|P W^\dagger P W-D_b\otimes I\|\le\delta_0,
\qquad D_b|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle,
```

or a different complete source interface, together with qualifying
native depth and simultaneous resource bounds.

Two analytic tests decide what to try next.

* A single batch between arbitrary fixed bank-unitary endpoints has an
  exact vanishing third address difference. The target has a nontrivial
  eight-output Pauli difference, giving full-operator error at least
  $`1/4`$ for $`b\ge3`$, even up to a common phase. This holds
  independently of how many affine calls run in the batch.
* For an echo with a fixed Hermitian Pauli P on work only, h such rounds
  and actual Clifford interlayers preserving x and Y pointwise, fixing
  a valid one-hot Y table makes every branch Clifford. Ordered Pauli
  transport gives label degree at most h and literal modulo-four phase
  degree at most $`\max\{2h,1\}`$. Complete error below
  $`\sqrt2`$ forces scalar return, hence
  $`b\le\max\{2h,1\}`$. For $`b\ge2`$ this requires
  $`h\ge\lceil b/2\rceil`$.

The second proof bounds the echo, not the global scalar of an arbitrary
Clifford-valued word. It allows arbitrary joint work mixing by the actual
read-only Cliffords. It excludes bank/address motion, non-Clifford
interlayers, and other extraction interfaces. A round count is not a
lower bound on native T-depth after arbitrary whole-word resynthesis.
The four-layer upper schedule cannot be reversed into such an inference.

The next constructive question is therefore a complete selective word
whose intermediate dynamics leave this read-only Clifford-branch Pauli
echo, or a whole-word synthesis that avoids its sequential round cost.
A concrete candidate must exhibit its selected action and full return
before claiming a depth improvement. Merely accessing more affine forms
in parallel no longer addresses the missing operation.

## 3. Complete resource, error and return recurrence

The [batch proof](PARALLEL_AFFINE_DETUNING.md#1-four-native-layers-for-a-complete-affine-batch)
defines $`\chi=1`$ for even b and zero otherwise, and

```math
\begin{aligned}
M_\ell&=q_\ell m(b+1+\chi),\\
H_\ell&=q_\ell m(b+2+\chi)
             +m\sum_{j,a}\mathrm{wt}(L_{\ell,j,a}),\\
C_\ell&=m\sum_{j,a}t_{\ell,j,a}.
\end{aligned}
```

On n_S arbitrary signal wires, use a reusable private pool of
$`4\max_\ell M_\ell`$ arbitrary wires. Let the interlayers have
common forward/actual-inverse costs T_inter, G_inter, D_inter, additional
arbitrary width w and total full-operator error eta. The same emitted
word and actual adjoint have

```math
\begin{aligned}
T_W&\le15\sum_\ell M_\ell+T_{\rm inter},\\
G_W&\le\sum_\ell(22M_\ell+16H_\ell+64C_\ell
                         +3q_\ell m+84)+G_{\rm inter},\\
D_T(W)&\le4h+D_{\rm inter},\\
b_W&\le n_S+4\max_\ell M_\ell+w,\qquad
\|\widetilde W-W\|\le\eta,\qquad a_{\rm clean}=0.
\end{aligned}
```

For $`P=\prod_a Z_{s_a}`$, a separately proved ideal selective
echo at error delta-zero yields, after the paid output Hadamards,

```math
\begin{aligned}
\epsilon_{\rm ind}&\le\delta_0+2\eta,\\
T_{\rm ind}&\le2T_W,\qquad D_T(\mathrm{ind})\le2D_T(W),\\
G_{\rm ind}&\le2G_W+4n_S+2m,\\
b_{\rm ind}&\le n_S+4\max_\ell M_\ell+w.
\end{aligned}
```

The added relocation pool returns exactly after each batch. Original
signals and interlayer work return only through the hypothesized complete
echo; their error is included above. Arbitrary correlations and external
references are included. The earlier arbitrary-angle and one-signal
amplitude-$`1/\sqrt2`$ recurrences remain in the
[previous continuation](JOINT_SELECTIVE_CONTINUATION_RESULT.md#3-complete-resource-error-and-return-recurrence).
The latter still costs 24 source/actual-inverse calls, depth at most
$`24d_U+32`$ and error at most
$`16(3\sqrt2\epsilon+2\epsilon^2)`$; exact flatness separately
sharpens error and signal/toggle return. No qualifying source is supplied
by either ledger.

## 4. Repair the multioutput inference; retain the original-pair scout

The [word-valued four-band proof](READ_ONLY_THREE_LAYER_RIGIDITY.md#87-word-valued-queries-the-outer-shear-is-additive)
now proves that every component of a simultaneous read-only XOR query
has degree at most two below error $`\sqrt2/12`$, with arbitrary
output/work inputs and zero clean qubits. It uses native outer-axis
coordinate blocks, the vector cocycle and the inner-axis quotient; it
does not initialize the other outputs. Consequently the original full-call
transfer now excludes at most four calls for $`b\ge6`$ at that
floor. The at-most-three-call result for $`b\ge4`$ retains its
stronger floor $`\sqrt2/6`$. Neither is an unrestricted-address
result.

The [larger physical-support analysis](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#8-larger-physical-support-a-two-band-sector-encoder-and-a-one-band-obstruction)
identifies why simply fixing an original pair parity does not generalize
the previous two-axis proof. An endpoint quotient has at most one axis
relation, but this dependent list is not one native band. For a retained
all-active even-support family, the interior clean Clifford has nonzero
alternating symplectic displacement. Every nontrivial one-native-band
clean image instead has nonzero characteristic functional. The sector
code therefore has no one-band encoder, including arbitrary Clifford
endpoints and returned dirty spectators, with complete-isometry error
floor $`1/(4\sqrt2)`$.

There is a paid repair: emit a maximal independent part as one band and
the remaining parity phase as a second band, with unchanged T count,
no additional work, literal phases and the actual inverse. General
one-pair sectors then have at most five bands when middle singletons
survive, or four in the pair-only witness. This does not place them under
the retained three-band no-pair degree theorem. The obstruction belongs
to the sector representation; the witness is not a high-degree Boolean
family and supplies no unrestricted depth bound.

The next focused lower question is finite-insertion stability of the
complete occupied Boolean condition, or a different native matching
argument on the original physical pair. Deleting the insertion does not
preserve Booleanity automatically. Keep actual lists, signs, the central
partner and return throughout; do not replace this missing implication
by a basis change or free sector encoder. General original pairs,
growing pair count and unrestricted depth remain open.

## 5. Verification boundary

The affine schedule uses an existing exact native primitive, and both
new detuning obstructions are analytic. No new numerical search is
needed. The four-band extension adds one seven-wire exact simultaneous
quadratic-query receipt, all 128 columns and the actual inverse, with
$`T=16,D_T=4,G=36`$ and zero clean work. The physical-pair proof
extends its retained small witness with the paid sector encoder and
literal inverse checks. These fixtures support identities, not universal
or asymptotic claims. Independent audits checked the new proofs and
scopes. Required full validation, render and exact-tree publication
receipts are recorded on the pull request. The selected manuscript and
complete-frame frontier remain unchanged.
