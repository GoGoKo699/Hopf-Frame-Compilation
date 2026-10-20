# Continuation: amplitude-uniform source limits and a paid joint detuning

[Workspace](../../WORKSPACE.md) ·
[Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous revision](INDICATOR_FUSION_REVISION_RESULT.md)

This continuation starts from main
`98dff4157d4015fdeee79decfdef52f037e3e6b6`, after PR #202.
It closes decreasing amplitude as an escape from the bounded-clean
two-band source theorem, pays one joint detuning operation by a different
native mechanism, and strengthens the focused original-pair scout.
It does not construct a complete shallow indicator or improve the
unconditional general depth frontier.

## 1. Destination, established frontier, components and limits

| Role | Assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth with simultaneous count, width, error and return guarantees on the same complete circuit |
| Strongest unconditional upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md), with one clean qubit, the retained real/complex dirty reservations and complete initialized-isometry error; selected R79 remains Result D |
| Unrestricted lower information | The retained constant-depth Boolean obstructions and linear indicator T-count bound; no dimension-growing unrestricted depth bound follows here |
| Useful constructive component | [Root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) absorbs polynomial count/width overhead of a complete sublogarithmic-depth source. The new joint fixed-angle detuning below is an exact native component within those relaxed budgets |
| Source-specific limit | [Amplitude-free source defect](FIXED_CLIFFORD_SUM_SOURCE_GAP.md#7-an-amplitude-free-defect-gives-a-uniform-relative-gap): bounded-clean two-band sources have a positive relative-error gap uniformly over every positive amplitude at sufficiently large address size |
| Focused physical-pair result | [Paired rank-two endpoints](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#5-all-paired-sources-with-two-clean-active-axes-at-both-endpoints) have Boolean sign/XOR degree at most six/five. The fixed minimal core has the sharper three/two bound, including all admissible silent input extensions |

R96 still gives, for $`N=2^n`$ and precision L,

```math
\begin{aligned}
D_T&=O(NL/B^2+L+\log_2(n+2)),\\
T&=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
\end{aligned}
```

Here $`n\ge1,L\ge6`$, with real
$`B\ge17(L+n+7)`$ or complex $`B\ge17(L+n+8)`$.
The operator contract includes arbitrary dirty inputs and references.

## 2. The complete inverse-square-root-two recurrence stays paid

For $`m=2^p`$, the intended complete indicator is
$`I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle`$.
Let U be an actual native source with one initialized signal, insertion
J and complete block

```math
\|J^\dagger UJ-I_p/\sqrt2\|\le\epsilon\le1/8.
```

The [retained recovery and promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#6-a-one-signal-block-is-another-sufficient-source-interface)
give the simultaneous common forward/actual-inverse ledger

```math
\begin{aligned}
T_{\mathrm{ind}}&\le24T_U+48m+8,\\
G_{\mathrm{ind}}&\le24G_U+82m+104,\\
D_T(\mathrm{ind})&\le24d_U+32,\\
b_{\mathrm{ind}}&\le m+b_U+2,\qquad a_{\mathrm{clean}}=0,\\
\epsilon_{\mathrm{ind}}&\le16(3\sqrt2\epsilon+2\epsilon^2).
\end{aligned}
```

The final eight T gates can be omitted for even m. The error is on the
entire operator, including return of the borrowed signal and toggle,
all original work and arbitrary references. There are 24 source or
actual-inverse calls. Positive-S inverse costs and literal source phases
are included. Under the additional exact-flatness hypothesis, promotion
instead has error at most $`2\sqrt2\epsilon`$ and exact borrowed
signal/toggle return. That stronger return cannot be assumed for a
general deeper source.

The first unpaid operation remains the specified logical transition U.
Putting $`d_U=3`$ gives a conditional depth bound 104, not a native
construction. Changing its amplitude was another possible escape;
the new invariant closes that escape when the clean count stays bounded.

## 3. What the amplitude revision establishes

For actual one-band encoders, write
$`K=J^\dagger U_B^\dagger U_AJ`$. The new physical defect restricts
exactly to

```math
K^\dagger P K-PK^\dagger K.
```

Its expansion contains at most $`2^{3k+1}`$ literal Cliffords with
coefficients $`\pm2^{-3k}`$, independent of amplitude. The fixed-sum
gap therefore supplies $`\tau_0(k)>0`$ and $`m_0(k)`$ such that

```math
m\ge m_0(k),\quad 0\lt\alpha\le1
\quad\Longrightarrow\quad
\|K-\alpha I_p\|>\tau_0(k)\alpha.
```

The proof permits arbitrary source width/count, common phase and rejected
action. No flatness or clean-return assumption is added. For one signal
there is also the sharp finite-size alternative

```math
K=0\quad\hbox{or}\quad \|K\|\ge\sin(\pi/8).
```

Neither statement is an amplitude-independent absolute error floor or a
bound for unbounded clean count, additional source bands or unrestricted
depth. The constants are conservative; not every modest relative
tolerance is excluded. The one-signal norm gap does not bound each
singular value.

## 4. Change mechanism when intact rotations keep the old cost

The [native detuning audit](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#51-a-complete-distance-certificate-and-the-native-scheduling-gap)
finds exactly b independent-axis batches for the unchanged address/signal
rotation list on b address bits. The reason is its incidence rank
$`m+b-1`$, where $`m=2^b`$. This is a limit of keeping the rotations
intact. It is not a lower bound for their jointly resynthesized product.

The [joint fixed-angle construction](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#52-joint-resynthesis-pays-one-special-angle-detuning-call)
changes the phase-polynomial representation. It implements the whole bank

```math
R=\prod_u\exp(i\pi d_u Z_{s_u}/4),\qquad
d_u=|x\oplus u|+1-Y_u,
```

exactly with T-depth at most 11, T-count at most $`5m+3b+5`$, Clifford count
$`O(m^2+mb+b^2)`$, m arbitrary signal wires, and no additional
work or clean qubit. Its actual inverse and literal common phase are
paid. The diagonal call can entangle the signal states with the data;
return of borrowed signals is an obligation of the eventual complete
source, not a property of this component. The proof uses parity and quadratic phase resynthesis; it does not
assume free parallelization of the old dependent axes.

This is one detuning angle, not a complete selective-inversion pulse.
Its period-eight dependence cannot by itself distinguish zero from
weight eight in a completed independent weight-only row word. The
periodicity argument does not exclude joint bank interlayers that access
the full detuning vector. The old continuous pulse also retains its
large finite-simulation cost; the new component does not retroactively
repair that recurrence.

The next constructive question is therefore concrete: can this paid joint
bank operation, together with explicitly charged native interlayers,
produce the complete selected-output sign with sublogarithmic depth and
polynomial overhead? A positive answer must give a full-bank error and
return proof. For an independent distance-profile design, the retained
weighted sum certificate is required; an arbitrary joint design instead
needs its own full-operator certificate. Logical phase-access bounds
retain their spectral premises and are not automatically T-depth bounds.

## 5. Keep the physical-pair scout scoped and active

The [new proof](MINIMAL_PHYSICAL_PAIR_COMPLETION.md) uses actual native
axis independence and the common original-pair parity. When both endpoint
clean-active ranks are two and a genuine middle physical pair is present,
there is exactly one such pair. Fixing its common parity gives a literal
zero-clean three-band word in each sector. The retained sign-degree-five
theorem in each sector gives total sign degree at most six and XOR degree
at most five. Sector restriction is a proof device, not a free coherent
compiler. The existing approximation-to-exact theorem applies to the
same circuit and supplies its conservative scoped error floor.

This strengthens a previously bounded fixed-rank branch; it does not solve
the original-pair problem at growing endpoint rank. The no-pair branch
retains its own theorem and is not assigned the new six/five bound.
The next lower question is native matching and propagation of the actual
exceptional module beyond rank-two endpoints, preserving original-pair
provenance and all silent axes.

## 6. Verification boundary

Two new exact source fixtures check the amplitude-free factor order and
sharp compression boundary. One additional physical-pair fixture checks
the coupled silent input and both parity sectors. The joint detuning
fixture checks a small complete native word, literal phases and actual
inverse. These are bounded identity checks; the uniform theorems are
analytic. Required validation and publication results are recorded on the
reviewed pull request. No selected manuscript theorem was changed.
