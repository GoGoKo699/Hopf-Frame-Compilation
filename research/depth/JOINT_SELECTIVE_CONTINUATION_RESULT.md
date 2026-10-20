# Revision: the complete selective word remains unpaid

[Workspace](../../WORKSPACE.md) ·
[Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous continuation](AMPLITUDE_UNIFORM_SOURCE_RESULT.md)

This continuation starts from main
`9ecaa37d0c8be1a2c460de57425a04e21697dd0a`, tree
`2c9309c4f810d6bbfd63df9161d11ddd64f61f96`, after PR #203,
on 20 October 2026. It revises the joint-detuning candidate before
continuing it, and keeps the original physical-pair lower scout active.
No complete shallow indicator or stronger unrestricted depth bound is
claimed.

## 1. Separate the destination from the current mechanisms

| Role | Current assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth, with count, width, error and return on the same complete circuit |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md); selected R79 remains Result D. The logarithmic-versus-constant fixed-accuracy depth gap remains open |
| Unrestricted lower information | Retained Boolean layer obstructions and the [linear indicator T-count bound](INDICATOR_PAULI_NULLITY.md). The latter does not imply growing depth |
| Useful constructive components | [Root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) normalizes polynomial source overhead; [clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) pays the amplitude-source wrapper; the exact fixed-angle bank costs at most 11 T layers |
| Already closed source interface | Uniformly bounded-clean two-band sources fail at sufficiently small fixed relative error for all sufficiently large addresses, uniformly over every positive amplitude. This includes overlapping sources at amplitude one over square root two |
| New representation limits | Joint Clifford decoding, the specified read-only-bank pulse grammar, and independent fine-angle synthesis each stop under their own hypotheses. None is an unrestricted depth lower bound |
| New focused lower results | With an original middle pair, either endpoint active rank two gives sign/XOR bounds 197/196. A separate actual two-axis support condition closes a growing-rank one-pair family |

R96 retains

```math
D_T=O(NL/B^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
```

Here N is two to the n, n is at least one, L is at least six, and the
real/complex dirty reservations remain respectively
$`B\ge17(L+n+7)`$ and $`B\ge17(L+n+8)`$, with one clean
qubit and complete initialized-isometry error at most $`2^{-L}`$.
The full operator contract remains mandatory throughout the discovery
route: arbitrary address, output, dirty-work and reference inputs,
literal phase, actual inverses, leakage and return.

## 2. Locate the first unpaid operation

For b address bits and $`m=2^b`$, the paid component is

```math
R=\prod_u\exp\!\left(\frac{i\pi}{4}
 (|x\mathbin\oplus u|+1-Y_u)Z_{s_u}\right).
```

Its [joint native resynthesis](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#52-joint-resynthesis-pays-one-special-angle-detuning-call)
has $`T\le5m+3b+5`$, $`D_T\le11`$,
$`G=O(m^2+mb+b^2)`$, no clean qubit and no helper beyond the
m arbitrary signals. It can entangle those signals; it is not a return
certificate. Its complete vector distinguishes addresses, but that
observation supplies no native decoder.

The first unpaid operation is a complete selective word W whose echo
implements the selected output phase on all inputs and returns all work.
The following analytic tests now locate several unsupported shortcuts.

| Proposed shortcut | Proven limit and reopening condition |
| --- | --- |
| Normalize by the address-zero call, then use joint Clifford decoding | [Four-point identity](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#53-joint-clifford-decoding-of-the-relative-call-still-stops) gives error at least $`1/(2\sqrt2)`$ for b at least four, with any number of normalized calls and arbitrary fixed unitary endpoints. Internal non-Clifford interactions or another access representation remain outside |
| Insert just one internal one-band encoder into that normalized word | An eight-point identity gives error at least $`1/(4\sqrt2)`$ for b at least six, still allowing fixed arbitrary endpoints. Additional internal interactions are outside, but are not a construction |
| Use affine address presentations of normalized calls | If r parity features and s quadratic features suffice, accuracy below the four-point floor requires $`r+s\ge b-1`$. This is an access-feature constraint, not native T-depth |
| Use at most three original full calls with read-only address and joint bank Clifford interlayers | The actual restricted word falls under the retained word-valued read-only three-band theorem: for b at least four its error is at least $`\sqrt2/6`$. The four-band theorem is stated for one output; a componentwise extension is unproved |
| Keep all computational bank values read-only throughout a full-call pulse | [Degree-four closure](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#54-keeping-every-bank-value-read-only-still-stops-full-call-words) requires $`r\ge b-4`$ at error below $`1/(2\sqrt2)`$, where r is the rank of accessed parity forms. Repeated original calls have r equal to one and fail for b at least six, regardless of call count. Bank motion across call boundaries and other native dynamics remain outside |
| Replace the angle by finer independently compiled rotations | [Precision theorem](JOINT_DETUNING_PRECISION_LIMITS.md) retains $`\Omega(\log b)`$ depth for the stated dyadic range, and $`\Omega(b)`$ for every fixed finer angle. Joint non-Clifford synthesis, additional helpers and whole-pulse error cancellation are outside its premise |

The normalized identity is

```math
C_x=R_xR_0^\dagger
 =S_{\rm all}^{-P(x)}Z^{A(x)+Q(x)\mathbf1},
\quad A(x)_u=x\cdot u,\quad Q(x)=\binom{|x|}{2}\bmod2.
```

Here P is address parity. The relative symplectic action only sees P;
the remaining Pauli labels have degree two. The read-only-bank theorem
uses a different restriction: a valid fixed one-hot Y input makes every
full call a Clifford on arbitrary signals. Its phase on a common parity
fiber has degree at most four. This argument is about the specified
grammar, not arbitrary intermediate address/output motion.

The finer-angle exact field obstruction likewise supplies no approximate
depth theorem for a joint circuit. The independent-layer proof instead
uses its full tensor-operator promise: projective local errors must sum
to a constant, so some coordinate needs error of order one over m.
That precision cost is absent from a bare logical phase-access count.

## 3. Complete resource, error and return recurrence

The [arbitrary-angle ledger](JOINT_DETUNING_PRECISION_LIMITS.md#1-the-arbitrary-angle-component-and-its-paid-recurrence)
sets $`q=m(b+2)`$. For an h-call word with local rotation synthesis
costs $`t_{\epsilon_j},g_{\epsilon_j},d_{\epsilon_j}`$ and
charged interlayers, it proves

```math
\begin{aligned}
T_W&\le q\sum_jt_{\epsilon_j}+T_{\rm inter},\\
G_W&\le q\sum_jg_{\epsilon_j}+2hm(b+1)+G_{\rm inter},\\
D_T(W)&\le(b+1)\sum_jd_{\epsilon_j}+D_{\rm inter},\\
\eta&\le q\sum_j\epsilon_j+\sum_j\rho_j,\qquad
b_W\le m+w,\quad a_{\rm clean}=0.
\end{aligned}
```

The same majorants charge the actual adjoint. If the ideal word separately
proves the complete selective echo at error delta-zero, its native
replacement gives

```math
\begin{aligned}
\epsilon_{\rm ind}&\le\delta_0+2\eta,\\
T_{\rm ind}&\le2T_W,\qquad G_{\rm ind}\le2G_W+6m,\\
D_T(\mathrm{ind})&\le2D_T(W),\qquad b_{\rm ind}\le m+w.
\end{aligned}
```

The error includes every signal and helper return. The hypothesis on the
ideal echo is the missing selective operation, not a conclusion of this
ledger. Fixed-angle exact calls use their sharper 11-layer costs instead
of the arbitrary-angle fallback.

The alternative [one-signal source recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#62-the-next-source-band-retains-the-general-leakage-recurrence)
retains, for a complete block approximating $`I_p/\sqrt2`$ to epsilon,

```math
\begin{aligned}
T_{\rm ind}&\le24T_U+48m+8,\\
G_{\rm ind}&\le24G_U+82m+104,\\
D_T(\mathrm{ind})&\le24d_U+32,\qquad
b_{\rm ind}\le m+b_U+2,\\
\epsilon_{\rm ind}&\le16(3\sqrt2\epsilon+2\epsilon^2).
\end{aligned}
```

This is zero-clean promotion on arbitrary inputs. Return remains
approximate unless exact flatness is separately proved; then the retained
sharper $`2\sqrt2\epsilon`$ bound and exact signal/toggle return
apply. The general overlapping two-band source cannot meet this interface
at every fixed relative tolerance by the amplitude-uniform source gap.
Adding a band only removes a theorem hypothesis; it pays no transition.

## 4. Physical-pair scout: growing ranks with actual provenance

The [one-sided rank-two theorem](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#6-one-rank-two-endpoint-permits-the-other-active-rank-to-grow)
allows either endpoint active rank to equal two and the opposite rank to
grow without bound, with an original physical middle pair. When the
opposite rank exceeds two, negative singletons remain. On each common
pair-parity sector, an actual
negative singleton gives the shared clean flip, and the three comparison
phase lists have relation dimensions zero, zero, one. The retained
bounded-relation theorem gives degree 196 per sector, hence Boolean
sign degree at most 197 and XOR degree at most 196. The dependent output
list is not credited as one free native band. Both-rank-two six/five
remains sharper in its narrower class.

The [two-axis support extension](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#7-two-axis-pair-differences-with-both-endpoint-ranks-growing)
permits both ranks to grow when there is exactly one original pair and
its difference is the sum of exactly two axes in each actual endpoint
active list. The parity sectors become actual no-original-pair sources,
giving sign/XOR bounds $`D_{\rm np}+1`$/$`D_{\rm np}`$ from the
retained no-pair constant. Support in the physical list is essential;
choosing another basis of the same span does not meet the hypothesis.
Both extensions retain the conservative higher-degree floor
$`2^{-530690}`$ by rounding the original native three-band circuit
before taking any sector comparison.

The next focused lower question is an original one-pair source with both
active ranks growing and a pair difference of larger physical support at
one or both endpoints. Preserve actual native matching and the central
partner; no bound on a silent-axis span is assumed. General fixed pair
count, growing pair count and unrestricted depth remain open.

## 5. Research decision and verification boundary

Continue the complete approximate source problem with joint native bank
dynamics, address/output motion, or whole-word resynthesis that leaves
the stopped representations. First exhibit the selective operator and
return identity; then prove simultaneous sublogarithmic depth and the
polynomial count/width budgets permitted by root fusion. Resolving only
angle compilation, another band count or a decoder outside a proven
obstruction does not change the full resource frontier.

One new five-wire physical-pair fixture checks ranks two and five, three
negative singletons, a four-axis quotient relation, all original physical
columns, both parity sectors, the actual inverse and shared-flip return.
It is not a Boolean-target counterfamily. The other new statements are
analytic; no large numerical search was added. Independent internal
audits removed an unsupported four-band multi-output inference. Required
validation, rendering and exact reviewed-tree publication evidence are
recorded on the pull request. The selected manuscript is unchanged.
