# Revision: pay parallel promotion, then test the complete source

This continuation starts from main commit
`af3a52e67c62746ae1d7dc137158367b9ef0ff22`, tree
`e81bcec94f14f93830a0e6ef5dbbaa10ab1400a7`, merged as PR #214 on
20 October 2026. The [preceding revision](GAUSSIAN_AND_CENTER_REVISION_RESULT.md)
separates small-clean Gaussian and involutory scalar restrictions from
the unrestricted model. The present results do not improve the general
resource frontier. They identify which new reductions are paid, which
candidate words fail, and which operation must change next.

## 1. Destination and present evidence

| Role | Assessment |
| --- | --- |
| General destination | Optimal complete-frame native T-depth under simultaneous count, width, Clifford, error and return budgets |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md) is unchanged; selected [R79](UNIFORM_SUBTREE_DEPTH.md) remains Result D |
| Complete constructive target | At every fixed tolerance, a zero-clean dirty indicator with linear exponential T-count/width, quadratic exponential Clifford count, and depth sublogarithmic in address size |
| Retained lower frontier | Constant-layer Boolean rigidity and the original three-band paired sign/XOR bounds six/five, or five/four at high active rank. No dimension-growing unrestricted depth lower bound is obtained here |
| New useful upper component | Exact native packing pays parallel promotion in a larger growing-clean window; a separate constant-clean spectral extraction remains conditional on a whole-space gap |
| New exact source identity | An odd-key counter uses one clean reference bit with arbitrary dirty lower bits. It supplies a flat compression identity, while its shallow joint arithmetic and full-bank schedule remain unpaid |
| New source assessment | Exact crossed-kernel collision equations expose complete arbitrary-key obstructions. A paid constant-depth word can fail the required observable by a constant norm gap |
| New lower component | Actual support-two paired factors admit paid one-band resynthesis. A Pauli-preserving-domain criterion rejects a specific higher-support center-moving Boolean completion |

The complete contract retains arbitrary coherent address, output-bank,
dirty-work and reference inputs, literal phase, actual inverse and all
leakage/return error. None of the new constructions is justified by a
favorable column, a trace, or an initialized private bank.

## 2. Parallel promotion: two distinct sufficient mechanisms

The [parallel-promotion proof](PARALLEL_CLEAN_PROMOTION.md) starts from
a k-clean complete query V. Every actual Clifford boundary must fix the
shared computational address, and every T target must be private to a
source copy. This is a native scheduling hypothesis, not an inference
from commutativity of completed source words. The source-clean bank may
move inside each copy.

Use h equal to two to the k private copies on arbitrary inputs. A paid
selective echo applies the desired phase only when a copy's clean bank
is zero and its arbitrary weight bit is one. Its quartic phase is emitted
with 32 T gates per output coordinate, private T targets and two returned
dirty helpers. All source calls, zero-tests and phase gates are charged.

An exact basis permutation W packs the odd-parity eigenspace so that
two further echoes extract the desired phase and return every borrowed
input. The following construction now pays W in a stated clean-count
window. With the notation of the proof home, the same completed circuit obeys

```math
\begin{aligned}
T_{\rm final}&\le4hT_S+8T_W,\\
G_{\rm final}&\le4hG_S+8G_W+12+2q,\\
d_{\rm final}&\le4d_S+8d_W,\\
b_{\rm final}&\le h(k+3q+2+b_V+b_C)+1+b_W.
\end{aligned}
```

For exact source-sector preservation, a refined monomial packing sends
only zero-active or one-active source configurations into the sector
used by the final echo. The other sector cancels with its actual inverse.
The full error is therefore at most twice the source error, independent
of h. For exactly flat amplitude-$`1/\sqrt2`$ compression, the retained
recovery gives complete error $`2\sqrt2\epsilon`$; its extra recovery
flag is included in k. General leakage has the separate bound
$`32h\epsilon`$. Approximation of W in full norm by mu adds eight mu.
These statements include every private return and the actual inverse.

The packing is now paid by exact native synthesis. On
$`N=2^k(k+1)+1`$ packing wires, every permutation decomposes into at
most 2N minus one single-target Boolean gates. Each such gate is a
returned parity echo around two calls to the retained exact dirty
indicator on N minus one controls. This returns its entire arbitrary
output bank and helper pool and has no source phase. The complete W
has depth $`O(N\log(N+2))`$ and count/width $`2^N\mathrm{poly}(N)`$.
Its exact inverse is charged, and no packing error is introduced.

Consequently

```math
2^k(k+1)^2=o(\log(p+2))
```

is a sufficient paid packing window. Its count and width overhead are
$`p^{o(1)}`$. A native-shareable source with
$`k\sim\beta\log_2\log_2 p`$ and depth
$`O((\log p)^\gamma)`$ now needs only separate conditions
$`\beta\lt1`$ and $`\gamma\lt1`$, instead of the serial sufficient
condition $`\beta+\gamma\lt1`$. Flat or exactly sector-preserving
sources retain constant requested accuracy. General leakage still
needs the stated shrinking accuracy envelope. The proof home gives
the fully substituted amplitude-recovery, promotion and fusion ledger.

This changes the conditional resource window, not the unrestricted
frontier: a qualifying source is still absent. At k of order log p,
this generic packing has excessive count and depth. Its lexical
bounded-degree monomial-round obstruction also remains valid; it is
not a lower bound for arbitrary native circuits or other packings.

There is a separate alternative. A constant number of global Clifford
conjugates of the selective reflection would suffice if their positive
projectors summed to an operator with a constant gap on the entire
arbitrary private space. The proof pays controlled queries without
controlling V, projector exponentials, actual inverse, rotation
synthesis and fixed-point filtering. The filter cubes its failure
amplitude while reusing the same two clean signals. Constant-clean
promotion and a final echo give a complete zero-clean query with
error at most 64 times the filter error.

The spectral precision schedule closes, but needs source error of order
one over h. No global-frame gap theorem is supplied. Local-frame or
average-state arguments do not prove the stipulated whole-space gap.
The [frame audit](SPECTRAL_PROMOTION_FRAME_LIMITS.md) proves separate
failures for common-Lagrangian frames, bounded physical support cones,
and specific dense graph/trace frames. Ordinary block depth is not
native T-depth, and dense noncommuting global frames remain open.
Neither route supplies a qualifying clean source. At k of order log p,
a better extraction would still be needed.

## 3. An exact counter source without a precise phase state

The [coherent-counter identity](COHERENT_COUNTER_SOURCE.md) uses
$`k=\lceil\log_2(p+2)\rceil`$ clean counter wires and one clean bit
in a second k-bit modular bank. Preparing that bit in minus selects
the direct sum of all odd Fourier keys; every other bit is arbitrary.
Each odd key is invertible modulo $`2^k`$, so the same coherent zero
test works in every occupied key sector and returns their coherences.
A row echo tests $`|x\oplus u|+1-Y_u=0`$ exactly. The complete product
returns all work and applies $`(-1)^{Y_x}`$. A final Hadamard on the
returned reference bit gives literal flat compression at amplitude
$`1/\sqrt2`$, without adding another clean signal.

This is a complete source identity with its own ledger, not a new
shallow source or a two-band claim. The modular operation must add
counter value times Hamming weight into an arbitrary accumulator.
Its present serial arithmetic schedule retains the old cost, and
reusing the clean counter sequentially across all rows is too deep.
Separate initialized row counters violate the clean allocation.
The full product's exact Hamming-shell spectrum explains its logical
commutation but does not provide a native parallel schedule.

The new odd-key mechanism pays reference preparation exactly. It leaves
a concrete joint arithmetic/whole-bank operation to construct. Its
clean count is of order log p, outside the generic packing window above,
so combining the two components does not yet improve the frontier.
The optional truncated-gradient lemma and the scope of recent
rotation compilers are retained separately in the proof home.

## 4. Crossed non-Gaussian kernels: the unpaid inference is readout

The [crossed-kernel proof](CROSSED_CUBIC_KERNEL_SOURCE.md) gives the full
matrix kernel and an exact collision criterion for a principal
permutation with arbitrary unitary helper fibers. An inverse map in a
Fourier phase is not such a permutation. The addressed-observable
contract is weaker, but its wrong-parity blocks must still vanish in
operator norm on all private inputs.

The explicit cubic-chirped trial has a paid native word with depth at
most 336 and polynomial overhead in the permitted exponential width.
Its temporary query work returns exactly. Its complete observable
nevertheless has gap at least square root of two: the candidate
preserves a dirty-key sector which the target flips. This is a full
operator failure, not a missing cleanup estimate.

The proof distinguishes several scopes. Computationally diagonal
middle words retain that sector regardless of their phase. A symmetric
plain-XOR middle has an exact uniform-key witness with gap two. More
general address-read-only translations with arbitrary phases retain a
square-root-two gap when the outer probe derivative is unbalanced. A
particular cheap balanced repair moves the obstruction to another
physical bit. These are restrictions on the stated kernels and middle
actions. They do not exclude all balanced chirps, arbitrary two-bank
motion, other outer frames, or unrestricted native circuits.

The useful next question is an actual whole-word observable identity
outside those preserved-sector hypotheses. Merely making the Fourier
kernel more nonlinear, or paying the already-shallow forward query
again, does not supply the missing readout.

## 5. Original physical pairs: remove the old cost before counting features

The [center proof](FOUR_BAND_CENTER_TRANSPORT.md) extends its paid
directional derivative to approximate sources: an error-epsilon source
gives an error-two-epsilon derivative, including return. Below the
retained robust two-band threshold, the indicator imposes the same
cubic Reed–Muller restriction on the projected physical-center mask.
It still does not impose a total-degree ceiling on invariant terms.

For an all-active endpoint and actual fully matched middle lists whose
single original pairs have support two, the conjugated middle factors
each resynthesize to one actual native band with paid Clifford
exposure. The complete four-band word reduces to two bands, preserving
literal phase, unused-block action, actual inverse and return. This
mechanism retains the old cost. Reframing it cannot create a new depth
obstruction.

The [higher-support continuation](HIGHER_SUPPORT_PAIR_BOOLEANITY.md)
first proves where that local resynthesis stops: a support-four pair,
or a support-two pair with relevant passive axes, has a Pauli image
whose support excludes a Clifford image. That is a local C3 test, not
a lower bound for the completed word.

A five-wire trial then has exact clean return and genuinely moves the
original physical center. It nevertheless cannot be made Boolean by
independent endpoint Cliffords. Its four-logical-qubit Pauli-preserving
domain has only two dimensions, whereas any diagonal unitary preserves
a four-dimensional logical Pauli Lagrangian. The deficit survives
arbitrary spectator mixing. The analytic proof retains the inactive
sector phases; the claim concerns this literal candidate, not every
higher-support pair.

For this unchanged middle on the same five wires, the proof also
excludes commuting D-phase corrections and every one-band output
endpoint refactorization. The next scout must change the actual middle
dynamics or add new native quantum coupling, and prove complete Booleanity. Clean return and center motion
alone are insufficient. A candidate's Pauli-preserving domain provides
a useful necessary test before trying to choose Boolean endpoint
coordinates.

## 6. Research decision

Keep the complete shallow source as the primary constructive question.
The paid packing enlarges the sufficient clean-count/depth window.
The odd-key code pays exact reference preparation. These improvements
change the full frontier only with an actual source satisfying the simultaneous
native resource, precision and return envelopes. Keep the focused
original physical-pair scout active outside the resynthesized and
explicitly rejected classes. The unrestricted logarithmic-versus-constant
depth interval remains open.
