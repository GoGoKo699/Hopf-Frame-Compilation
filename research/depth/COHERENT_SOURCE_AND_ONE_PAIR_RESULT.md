# Revision: weaker coherent sources and uniform all-active one-pair closure

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous character-source checkpoint](CHARACTER_SOURCE_REVISION_RESULT.md)

This is a historical checkpoint. The later
[structural pair revision](STRUCTURAL_PAIR_REVISION_RESULT.md) replaces
its large all-active lower bound with sign/XOR degree six/five uniformly
and five/four when the common active rank is at least four; see the
[full-space diagonal proof](ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md).
The coherent-source results below are retained. This revision does not
change the unrestricted resource frontier or settle a general original pair.

This continuation starts from main
`e6a1223c46dd3186ce3be2ab194ff0a70a6c6414`, tree
`67978e842513a1dfdd8c93d5ef13bcd7cee53c33`, after PR #207.
The revision removes an unnecessary extraction echo and an unnecessary
source-helper return promise. It also closes the retained all-active
genuine one-pair lower class without bounded singleton coupling.
Neither result changes the unrestricted complete-frame resource frontier.

## 1. Destination, strongest theorems, components and scoped limits

| Role | Current assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth with simultaneous T-count, Clifford count, dirty width, complete operator error and return |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md); selected [R79](UNIFORM_SUBTREE_DEPTH.md) and their literal allocations are unchanged |
| Useful source reduction | [Direct output aliasing](OUTPUT_ALIASED_CHARACTER_SOURCE.md) supplies the complete dirty indicator through four source/actual-adjoint occurrences and cancels arbitrary quantum private-helper fibers |
| Useful shallow components | The complete quadratic-inverse XOR and algebra-key group operation are paid. A paired transform and coherent correction are equivalent sufficient source targets up to paid constant depth |
| First unpaid upper operation | A sublogarithmic native weak fiber source, paired transform, coherent correction, or complete extraction within the polynomial-overhead source budgets |
| Fourier approach limitation | The shallow kernel exposes a trace character, not a returned scalar. Its particular Pauli-Z scalarization needs exponentially many masks and a growing clean label bank |
| Stronger focused lower theorem | [All-active singleton closure](ALL_ACTIVE_SINGLETON_CLOSURE.md) gives sign/XOR degree at most 37,748,865/37,748,864 for one genuine physical pair with all-active difference, allowing arbitrary original singletons, silent lists, endpoint ranks and dirty width |
| Remaining lower question | General original pair differences, then fixed pair count. The all-active proof does not cover a pair difference omitting actual endpoint-active axes |

For $`N=2^n`$, integer $`n\ge1,L\ge6`$, R96 retains

```math
D_T=O(NL/B^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
```

There is one clean qubit, real/complex dirty reservations respectively
$`B\ge17(L+n+7)`$ and $`B\ge17(L+n+8)`$, and complete
initialized-isometry error at most $`2^{-L}`$, including arbitrary
dirty/reference inputs. At fixed accuracy and square-root dirty width,
the unrestricted upper/lower depth interval remains logarithmic versus
constant. Internal analytic audits and finite fixtures support the new
claims; they do not replace external technical review.

## 2. Four-call extraction with quantum helper fibers

Put $`m=2^b`$, $`r=b+1`$ and $`w=2m-1`$. The retained HRT
permutation has principal action $`c\mapsto f(c)`$ on w arbitrary
coefficient bits. A sufficient source now has the weaker full-space form

```math
V_0=\sum_c|f(c)\rangle\langle c|\otimes S_c,
```

with arbitrary unitaries $`S_c`$ on private dirty helpers. These may
depend on c, carry phases, and mix computational basis states. They need
not return before the next source occurrence. A full-operator source
approximation $`\|V-V_0\|\le\eta`$ is required; columnwise success
probabilities do not suffice.

Compute the subset-zeta transform on the actual arbitrary output bank Y.
Its nonconstant coordinates are precisely the carrier coordinates used
by the retained character. Alias these physical wires rather than copy
them into a disjoint carrier. Only m additional arbitrary carrier wires
are needed. The character and its actual inverse restore all aliased
wires and cancel the helper action. Reverse zeta, correct the constant
character, and conjugate by output Hadamards. This gives the complete
indicator on arbitrary address/output/dirty/reference inputs:

```math
\begin{aligned}
T_I&\le4T_V,& D_T(I)&\le4d_V,\\
G_I&\le4G_V+bm+2b+2m+6,& h_I&\le m+h_V,\\
\varepsilon_I&\le4\eta,&a_{\mathrm{clean}}&=0.
\end{aligned}
```

The same readout accepts the coherent correction with principal action
$`R(z,y)=(z,y+f(z)+yf(z))`$ and arbitrary private unitary fibers.
Its arbitrary-key echo then gives

```math
\begin{aligned}
T_I&\le4T_R,&D_T(I)&\le4d_R,\\
G_I&\le4G_R+bm+4b+2m+6,&h_I&\le3m-1+h_R,\\
\varepsilon_I&\le4\eta,&a_{\mathrm{clean}}&=0.
\end{aligned}
```

All inverses are actual adjoints. Exact fiber sources give exact final
return; approximate sources give the displayed complete error, including
return and reference correlations. No intermediate exact return is
assumed. The [proof home](OUTPUT_ALIASED_CHARACTER_SOURCE.md) carries
the native Clifford counts and full principal identities.

Root fusion uses eight completed indicators, hence 32 V or R occurrences
and error at most $`32\eta`$. Its depth is at most

```math
32d_V(b_0+1)+O_\tau(\log\log(p+4)),\qquad
b_0=p-\lceil a\log_2(p+2)\rceil,
```

or the corresponding d_R expression. Request $`\eta=\tau/32`$.
For fixed-L frame transfer, request
$`\eta_L=2^{-L}/(512K_*)`$, with the retained dimension-independent
query count $`K_*`$. This pays the same full frame error and literal
allocations. The retained source still has logarithmic depth, so the
improved constants and source freedom do not by themselves improve R96.

## 3. Where the next mechanism stopped

The paid quadratic XOR $`Q_k`$ gives the exact principal identity
$`P=F\otimes F^\dagger=WQ_kRQ_k`$, where W swaps the two banks.
It also gives the converse $`R=Q_kWPQ_k`$. Each conversion adds
at most $`2T_{Q_k}`$, $`2G_{Q_k}+3w`$ and 80 T layers, preserves
full-operator error, and permits arbitrary private fibers. P is a
sufficient source because its second arbitrary bank can be counted as
private work and is restored by the completed character. This is not a
free standalone in-place F construction.

The algebra-key operation $`(a,y)\mapsto(a,a+y+ay)`$ is paid in
constant depth with polynomial address overhead. Its key a is not the
coherent ordinary key z: the missing conversion is $`a=f(z)`$.
Ordinary Hadamards merely yield another quadratic permutation and a
bilinear phase. A paid phase query for the group Fourier kernel likewise
does not supply the corresponding single-register Fourier transform.

Pursuing that kernel produces a native involution on two w-bit banks
whose every diagonal entry is $`p/2^w`$, with
$`p=(-1)^{\alpha\cdot f(\beta)}`$. Its full-operator distance
from that scalar is $`1+2^{-w}`$. Thus trace access is the first
unsupported inference in treating it as a returned scalar source.

The explicit Pauli-Z twirl repairs the scalarization, but has $`2w`$
clean labels and amplitude $`2^{-w}`$. For relative operator error
$`\epsilon\lt1`$, any linear combination of K of these mask conjugates
needs $`K\ge(1-\epsilon^2)2^{2w}`$. Standard promotion or
amplification retains exponential cost. This stops that representation,
not arbitrary native interference or another completed extraction word.
The classical inverse-degree argument also concerns completed classical
stages, even with classical helper garbage; it does not cover the newly
permitted quantum fiber action.

The next upper question is therefore a native weak fiber source or whole
extraction using genuinely nonclassical interference within the complete
source budgets. Another classical conversion or trace-only character
would retain the cost just identified.

The original amplitude-$`1/\sqrt2`$ route remains a distinct
conditional interface: 24 source occurrences, depth $`24d_U+32`$,
error $`16(3\sqrt2\epsilon+2\epsilon^2)`$, with approximate
signal/toggle return unless exact flatness is separately established.
The actual two-band, bounded-clean relative-gap theorem still applies
to that source class. No third-band source or unrestricted depth
obstruction has been supplied by the present revision.

## 4. The lower scout changes mechanism and closes growing coupling

The previous center proof had residual dimension at most $`2h+4`$
and needed fixed singleton-to-silent coupling rank h. Bounding h from
one pair remained unsupported. The new proof does not assert such a
bound. It instead uses the original pair's literal scalar parity branch
to obtain native matching modulo that single parity. On the other branch,
a central quadratic sign has rank k minus two, with k the actual common
active rank. Its nonaffinity forces the actual dressed Clifford to have
one computational X mask on the entire clean-reduced logical space.

A shared matrix-feature ledger with one core bit and nine extra features
then gives $`\deg g\le D_1+36=37748865`$ and XOR degree one less.
This allows unbounded h, both coupled silent lists, all original positive
singletons and arbitrary dirty/reference inputs. The physical clean
coordinate, literal phases and unrestricted unused clean block are
retained. The sharper rank-two six/five and one-sided 2428/2427 bounds
remain available in their subclasses. The retained original-circuit
rounding floor $`2^{-530690}`$ also applies when its rounded actual
physical lists meet this class's hypotheses.

The actual difference label must still be the sum of every clean-active
axis at both endpoints. General differences need not supply the same
common active span or central quadratic rank. That is the next focused
physical-pair question; multiple pairs, unrestricted one-clean three-layer
degree and dimension-growing depth remain open.

## 5. Focused evidence

The [output-alias fixture](../../tests/test_output_aliased_character.py)
checks every one of 2,048 physical columns and the actual reverse of a
two-address-bit phase indicator with a genuinely quantum conditional
helper. It checks the literal selected-output sign, helper return and
the four-call/Clifford ledger. The outer Hadamard conversion is analytic.

The [pair fixture](../../tests/test_all_active_singleton_closure.py)
checks the four-active-axis central-controlled quarter-turn on all 32
physical columns, two signed conventions, the literal scalar branch,
both clean sectors and actual inverse. It guards the structural identity;
it is not a Boolean query or an asymptotic degree witness. The uniform
theorems are analytic and do not follow from these finite tests.
