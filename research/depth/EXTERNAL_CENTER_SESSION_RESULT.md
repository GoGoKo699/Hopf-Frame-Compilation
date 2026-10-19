# Revision: retain the native middle and bound its probe action

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[External-center correction](EXTERNAL_CENTER_MIDDLE_CORRECTION.md) ·
[Phase-signal and frame interfaces](PHASE_SIGNAL_QUERY_LIMITS.md) ·
[Previous checkpoint](MIDDLE_DEFECT_SESSION_RESULT.md)

This continuation began at verified main `981ebc6` (PR #194) on
19 October 2026. The revision separates a useful change of representation
from a paid native simplification. The surviving center displacement can
be removed from the controlled Pauli simultaneously for every probe.
Its phase then enters the middle band, where a concrete family cannot
be realized by one native T layer, even with arbitrary Clifford endpoints.
A changed argument keeps the original middle band and closes the case
of bounded total probe-displacement dimension, even with growing active
ranks and quadratic support. This is a new degree component; the general
resource frontier is unchanged.

## 1. Destination and the present frontier

| Role | Retained conclusion |
| --- | --- |
| Destination | Optimal general depth for the prescribed complete frame, with simultaneous count, width, full-operator accuracy and return guarantees |
| Strongest general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ |
| Literal domain | $`n\ge1,L\ge6`$; one clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; initialized-isometry error at most $`2^{-L}`$ |
| Selected result | R79 remains Result D, with its settled count frontier at the stated allocations |
| Unrestricted components | At most one clean and two native layers give XOR degree at most two; zero clean and three native layers give sign/XOR bounds five/four |
| One-clean components | Fixed endpoint clean-active ranks close, including arbitrary middle pairs. In the actual no-pair branch, fixed quadratic-control rank, fixed exceptional overlap, and fixed original active-intersection dimension close uniformly in other ranks |
| Current candidate | A growing even active intersection carries the entire quadratic support; nonconstant flags leave the original active algebra and enter the enlarged shared center |
| New structural result | One actual middle Pauli corrects all these flags at once. The residual is a native band times a specified Boolean cubic, with Clifford probe defects |
| New degree component | In the retained actual no-pair branch, fixed dimension of the total middle probe-displacement span gives a uniform sign-degree bound, even as both active ranks, quadratic support and dirty width grow |
| Remaining limits | A one-band replacement of that residual is not automatic; clean phase signals and coherent Pauli frames do not supply a free complete dirty query |

At fixed accuracy and square-root dirty width, the complete-frame depth
interval remains a constant lower bound versus a logarithmic upper bound.
Even closing the full no-pair three-layer case would still leave growing
middle pairs and higher depths. The exact-rounding floor and every
full-input, reference, literal-phase and actual-inverse requirement remain
unchanged.

## 2. What the native matching now supplies

Retain the actual no-pair normal form and its original clean origin.
Let beta be the quadratic-control symbol and S its support. In the
remaining high-rank case, S is the whole common clean-active intersection.
It lies in the positive middle span and is disjoint from each full silent
endpoint span. This disjointness is stronger than merely fixing S.

Use the original independent middle axes as coordinates. The endpoint
probe actions have square-zero displacements N_u, and their control-phase
flags are p_u. The full cubic identity, restricted to their common fixed
space, gives a common preimage:

```math
N_uv=p_u\quad\hbox{for every }u.
```

Here v is one fixed product of actual middle axes, not a collection of
fiber-dependent coordinate choices. Consequently the literal Pauli
$`H'=HP_v`$ has only constant probe translations, and

```math
H^{b(z)}=(H')^{b(z)}P_v^{b(z)}.
```

The new middle operation is $`E_B P_v^{b(z)}`$. It is a commuting
C3 phase and its probe defects are Clifford. Its Boolean cubic correction
is retained exactly. Applying the three-native-band theorem at this point
would be the first unsupported inference: that theorem requires an actual
independent native middle band.

This distinction is substantive. For an explicit growing orthogonal
transvection family, every correcting v makes the resulting diagonal
middle phase impossible to realize by any one-T-layer circuit with
arbitrary Clifford endpoints and returned dirty helpers. The proof uses
the literal characteristic-parity phase relation, not a restriction to
one chosen Clifford basis. It is an obstruction for this specific residual
operation; it neither excludes its realization with more layers nor
proves that every original source produces this family.

A rank-one middle action also admits a useful exact quotient description:
after fixing the quadratic support, its activated axes reduce to two
nonzero labels, each with odd multiplicity. The global phases and their
dependence on the fixed controls remain. A small fiber is not yet a paid
global simplification.

The successful change of mechanism keeps the original native middle.
Write

```math
U=\sum_u\mathrm{im}N_u,\qquad d=\dim U,\qquad
k=\frac{d(d+1)}2+d+1.
```

Its entire common fixed Pauli space can be made a genuine computational
center after extracting actual diagonal endpoint sublists. The residual
target may be complex at this stage; a literal diagonal Clifford aligns
the two transported Pauli families. Only exact identity modes undergo
Hadamards. All removed phases are retained in the original affine
coordinates.

The endpoint axes that interact with the middle use a common input-mask
space of dimension at most d(d+1)/2. The middle needs at most d further
directions modulo the fixed center, and H adds at most one. Keeping
these factors intact leaves an actual two-native-band word preserving
a computational-Z subgroup of codimension at most k. The existing
partial-diagonal finite-feature theorem therefore applies. At control
polar rank at least four, the resulting sign bound is

```math
\deg g\le D_k+72\,2^k.
```

Here D_k is the explicit retained partial-diagonal bound. Rank zero
already closes at five, and rank two retains finite insertion; a bound
covering all ranks takes the maximum with
$`D_{\mathrm{ins}}(0,0,0,3)`$. This proves uniformity for fixed d
without declaring the cubic-corrected middle native or independently
synthesizing its fibers.

## 3. Changed upper mechanisms and their first missing operation

The previous polynomial and short-hash failures concern their particular
interfaces. They do not forbid interference. This revision therefore
examines coherent phase signals and coherent Pauli-frame transport.

A phase signal can distinguish zero Hamming weight from positive weight
on an initialized signal register. Replacing that register by arbitrary
dirty work and using two signal reflections leaves a common fixed
subspace. On that subspace the intended negative scalar is missing and
the literal operator error is two. More generally, a product of q
conjugated rank-one reflections on a d-dimensional signal space has a
fixed subspace whenever q is smaller than d. This is a rank bound on
that return mechanism, not an unrestricted approximate-query bound.
The signal preparation, native angles, inverse and full return must all
be charged before importing a clean-input depth guarantee.

A coherent Pauli-frame route encounters a different boundary. Full frame
transport for the complete indicator already contains a complete indicator:
two transport calls, one of them the actual inverse, and Clifford work
isolate the selected dirty-bank phase while returning every arbitrary
frame key. The full-operator error grows by at most a factor of two.
Thus that correction interface retains the query being sought. A cheaper
description of a frame is not itself a cheaper coherent correction.

Neither route constructs a zero-clean indicator with linear T-count and
dirty width, quadratic Clifford count and sublogarithmic address depth.
The constant-accuracy transfer to the complete frame remains available
once such a source is supplied.

## 4. Decision and reopening questions

The remaining lower question concerns **growing total displacement
dimension**, alongside the separate arbitrary-middle-pair branch.
The corrected-middle description still supplies exact Clifford defects,
but no one-band replacement may be assumed. The fixed-d proof instead
retains the original native lists and bounds the part of their action
that is not centralized. Extending that mechanism requires controlling
the completed target as d grows, or supplying a full Boolean native
family that explains why such control fails. A raw model or a complex
clean-returning source does not answer that question.

The upper question remains a jointly accurate complete dirty indicator.
A changed mechanism must act correctly on the entire signal/work space
or supply its initialization and return within the same budgets. Reopening
the phase-signal route requires overcoming the common fixed subspace;
reopening the frame route requires a genuinely cheaper restricted
correction interface together with a proof that the complete query uses
only that restriction.

The new statements have analytic proofs and bounded exact fixtures.
Those checks distinguish native phase identities from Gram-only models
and initialized signals from arbitrary dirty work. They do not substitute
for an asymptotic proof or change the established resource frontier.
Final manuscript drafting remains separate.
