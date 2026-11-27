# Exact key covers give sublinear depth at intermediate precision

[Geometric chunk depth](GEOMETRIC_CHUNK_DEPTH.md) ·
[Compact program precision](COMPACT_PROGRAM_PRECISION.md) ·
[Precision count–depth tradeoff](PRECISION_COUNT_DEPTH_TRADEOFF.md)

**Theorem R104.** For an absolute $`c\gt0`$, let $`n\ge1`$ and $`L\ge6`$
be integers with $`2^{cn}\le L\le N=2^n`$. Keep the real and complex literal
dirty-width thresholds of R96, one external clean input, and its
full initialized-isometry error $`2^{-L}`$. Put

```math
 S=\sqrt{N/L},\qquad
 t=\max\{1,\lfloor c_0\log_2(2+S)\rfloor\}.
```

One circuit has simultaneous bounds

```math
 T=O(\sqrt{NL}+L+NL/b),\qquad
 G=O(NL),\qquad
 D_T=O(NL/b^2+L\log_2(t+2)/t).
```

There is one sufficiently small absolute choice of $`c_0>0`$.
The dirty peak fits the original literal real threshold
$`b\ge17(L+n+7)`$ and complex threshold $`b\ge17(L+n+8)`$.
The actual reverse obeys the same error and bounds. If
$`L=o(N)`$ in this range, the additive precision depth is
$`o(L)`$ at every eligible allocation. At
$`L=N/n^u`$ for fixed $`u>0`$, it is
$`O(L\log\log n/\log n)`$. At $`b=\Theta(\sqrt{NL})`$,
the width term is O(1). If $`S`$ is bounded we use R96 and
the displayed additive term stays linear in L, compatible with
the R103 lower bound.

For every fixed exponent $`\alpha\in[c,1)`$, taking
$`L=\lceil N^\alpha\rceil`$ gives additive depth $`O(L\log n/n)`$,
with constants allowed to depend on $`\alpha`$ in this specialization.

<a id="exact-key-cover"></a>

## 1. Exact selection and a finite basis permutation

Use two disjoint subgroups of $`q=2^k`$ slots each, with arbitrary
k-bit keys $`A_i`$ and small arbitrary payload/worker banks.
The initialized clean key in the local comparison will have
$`k=O(t+1)`$ bits. For one subgroup put
$`a_i=[A_i=0]`$ and $`j=\bigoplus_i i a_i`$.
Declare its sole candidate valid if $`a_j=1`$.
All keys are unchanged by this classifier, even on arbitrary
inputs with many zero keys. In the uniform counting measure, the
event that exactly one key is zero has probability

```math
(1-1/q)^{q-1}\ge e^{-1}.
```

It ensures success. Choose the first successful subgroup, if
any. Their raw key registers are disjoint. Thus the bad fraction
for the complete classifier is at most $`(1-e^{-1})^2\lt1/2`$.
It need not agree with the event that there are no zero keys.
Use this exact classifier in both the permutation and the native
selected centers.

The [R100 parity-index predictor](GEOMETRIC_CHUNK_DEPTH.md#1-exact-priority-on-arbitrary-occupied-keys)
gives each candidate's read-only exact phase and its returned
XOR with $`O(\log(k+2))`$ T-depth and $`2^{O(k+1)}`$ count
and width. Concretely, its residual vector contains the k index
parities and the one validity literal. The two shifted one-hot
words on this short vector and the four Clifford corners test
all residuals zero without initialized work. Each subgroup
success is the parity of its mutually exclusive valid predicates.
Combine the second subgroup's candidates with the first group's
failure, using the same returned short threshold oracle. This
combines only constantly many group Booleans. The global good
flag g and every final selected-slot predicate are exact returned
XOR or phase visits on a common arbitrary key bank. Extra fixed
control literals, payload or Fourier bits enter as read-only
controls. With each T support private, all parallel overlaps on
the shared key and payload controls are Clifford interactions.
The count and width stay $`2^{O(k+1)}`$ and the depth stays
$`O(\log(k+2))`$.

### Exact guarded interval exchanges

On an m-bit integer z, the global modular reflection
$`z\mapsto a+b-1-z\pmod{2^m}`$ is an involution. It reverses the
interval $`[a,b)`$ and preserves its complement. For a stable
external predicate p, the conjunction $`\chi=p[z\in[a,b)]`$ is
therefore invariant. Compute its exact returned XOR E into an
arbitrary bit d. Expand the fixed thresholds $`z\lt a`$ and
$`z\lt b`$ into disjoint tests of the first differing bit, then
include p by a bounded mixed-difference echo. The complete tests
preserve z and use polynomial m resources plus constantly many
visits to the paid XOR of p and its actual reverse.

The d-controlled modular reflection F is a complement and a modular
addition of a fixed integer, both controlled by d. Pay the addition
by descending bit toggles for each fixed increment, with the carries
implemented by exact returned conjunctions. Its native cost and work
are polynomial in m. Chronologically execute E, F, its predicate
reverse $`E^\dagger`$, then $`F^\dagger`$. The second predicate visit
sees the same $`\chi`$. The two flips return d and the net reflection
exponent is $`\chi`$. Thus the restricted reflection is exactly
guarded even though its interval test reads its own targets.
Reverse $`[a,b)`$, $`[b,c)`$, then $`[a,c)`$ to exchange the two
adjacent intervals while preserving the order within each. All
comparators, guards and arithmetic use arbitrary returned work.

### Pack the two-subgroup failure set

Flatten the $`2q`$ raw key slots in a fixed order. At each step the
prefix slots to the left of the active word are still raw. In each
group g use the prefix XOR summary $`s_g`$ of zero-key indices and
the final whole-group candidate $`r_g`$ on the original input.
For an original suffix, its descriptor is $`(r_1,u_1,r_2,u_2)`$.
Status $`u_g`$ is unread if the candidate lies in the raw prefix,
and otherwise records whether its original key was zero. For each
pair of prefix summaries there are at most $`(3q)^2`$ descriptors.
Assign each a contiguous interval of length equal to its exact
class cardinality, using a fixed descriptor order.

Start with the empty suffix, whose sole descriptor has $`r_g=s_g`$
and both statuses unread. Prepend raw slots from right to left. For
fixed new prefix summaries and a fixed raw value of the slot being
prepended, the shorter suffix has known intervals. Its old prefix
summary differs from the new one by the slot's index if that slot
is zero. The final candidates are unchanged, and only the status
of the candidate at the newly prepended slot may become known.
That status comes from the prepended raw key. Thus the new suffix
has at most $`O(q^3)`$ contiguous intervals with known new labels.
Group them by descriptor using $`O(q^6)`$ adjacent exchanges.
The interval lengths give the next class sizes by summation.

Compile this rearrangement for each of the $`q^2`$ summary pairs
with an external guard computed from the raw prefix. That prefix
remains unchanged throughout every exchange. Its returned guard
uses the same atomic zero tests and short parity-vector equality
as the selection above, with polynomial q and k resources. The
arithmetic only changes the active suffix and restores all helpers.
At the root both statuses are known. Put the classes with two
failed candidates first. This gives an exact basis permutation K
packing the bad set into $`[0,B)`$ of the $`2^{2kq}`$ key labels.

The density bound gives $`2B\lt2^{2kq}`$. Exchange $`[0,B)`$ and
$`[B,2B)`$ and undo K by its actual reverse. For B zero, take the
identity. The resulting $`S_{\rm key}`$ is an involution and
$`P_2=S_{\rm key}^\dagger P_1S_{\rm key}`$ is orthogonal to the
original bad projector $`P_1`$. Both projectors are computational.
The $`2q`$ stages, $`q^2`$ summary pairs and polynomially many
interval exchanges on at most $`2kq`$ bits give
$`\mathrm{poly}(q,k)=2^{O(k+1)}`$ gates, dirty peak and T-depth for
the permutation and its actual reverse. Only the key bank moves.
The exact basis action and literal returned helpers include all
coherent keys and arbitrary dirty inputs with references.
The permutation lies outside the two fixed Clifford frames covered
by [the slot bound](GEOMETRIC_CHUNK_DEPTH.md#9-slot-lower-bounds-for-the-two-basis-recovery).

## 2. One exact key frame transports a long source in short pieces

Retain the exact geometric chunk of
[R100 Section 3](GEOMETRIC_CHUNK_DEPTH.md#3-exact-geometric-chunks-and-the-powered-source-loader).
A split chunk has entrance axis $`\Gamma_j`$, at most t
consecutive splits on one path, and desired exact normalized
involution $`M_j`$ on those axes. Extend the old core's list
for this argument to its full set of paired Hermitian Majoranas,
using all original wires and no extra initialized state.
Their chirality is proportional to the all-wire Z Pauli.
Two spare anticommuting axes $`\Gamma_a,\Gamma_{aa}`$
are that chirality tensored with X and Y on one extra arbitrary
line. They anticommute with all original axes, including the
unpaired original list or the inactive last partner of a paired
two-tail source. Define

```math
  G'_j=(\Gamma_a+M_j)/\sqrt2,\qquad
  G''_j=(\Gamma_{aa}+G'_j)/\sqrt2.
```

Both are the exact padded geometric lists of R100, with one
or two additional leading splits. Its dyadic signed-Clifford
menu and exact half-amplitude oblivious amplification implement
$`-G''_j`$ with a local $`O(t+1)`$-bit initialized
selector. The code action returns the entire selector and
arbitrary private scratch exactly, so the actual word preserves
the projector of the initialized key on the full space.
Use one global padded size $`k\le C_k(t+1)`$ between chunks,
and keep the same physical key bank and dirties throughout.

### Compress only a short physical interval

For a source on m original qubits, a chunk occupies a contiguous
interval J of $`O(t+1)`$ physical qubits. For a paired path,
consecutive axes on either tail advance by one physical position.
All the original participating axes have the same
preceding Z string outside J. Both spares have the same exterior
chirality Z string. Thus the exterior labels span at most two
pure-Z parities. Choose at most two independent pivot columns
and clear all remaining exterior columns using at most two
CNOT eliminations per column, entirely outside J. This supplies
a fixed exact $`O(m)`$-Clifford and its inverse, confining the
participating axes to the interval, at most two old exterior wires
and the one added spare. Avoid any core-wide quadratic
symplectic elimination per chunk. Within that active support the
exact menu costs $`2^{O(t+1)}`$ elementary gates and private width
and depth $`O(\log(t+2))`$. The change of coordinates fixes
all keys and guard controls.

### Apply an exact chunk conditioned on the one frame's success

Work first in the unpermuted key frame. Apply its native
[selected-center word](GEOMETRIC_CHUNK_DEPTH.md#2-a-selected-center-and-the-two-seed-bases)
on the SMALL active payload. Swap that payload into the unique
selected good slot with the three returned XOR sweeps. Run its
local $`-G''_j`$ code word in all private slots. A selected
$`e^{i\pi\Gamma_{aa}/2}`$ center, then the actual local reverse
and selected swap reverse, acts as $`iG'_j`$ on the primary
small payload on success, and exactly I otherwise. This uses
$`G''_j\Gamma_{aa}G''_j=G'_j`$. The unchanged key-zero
projectors exactly support the center and the actual swap.
All candidate payloads, scratch and keys return, including
arbitrary entangled inputs. Pay the literal phase $`(-i)^g`$
to obtain the exact partial $`G'_j`$.

The remaining two factors are $`\Gamma_a`$ and the plane
Clifford $`C_j=\exp(\pi\Gamma_a\Gamma_j/4)`$.
Indeed $`C_j\Gamma_j C_j^\dagger=\Gamma_a`$ and
$`G'_j\Gamma_a C_j`$ conjugates $`\Gamma_j`$ to $`M_j`$.
Every outside original axis is fixed because $`C_j`$ commutes
with it and each of $`\Gamma_a,G'_j`$ negates it under
conjugation. To guard these two factors by g, use their
literal Pauli rotations. Diagonalize the fixed Pauli with a
Clifford on the small compressed payload. A guarded
$`e^{i\phi Z}`$ at $`\phi=\pm\pi/4,\pi/2`$ is a product
of eighth-root phases with predicates g and g times the data
bit. Use R100's exact returned derivative on arbitrary dirty
targets. For the Pauli itself correct the good-sector i by the
literal $`(-i)^g`$ phase. Unconditional conjugating Cliffords
cancel with their actual inverses on failure. The same parity
compression surrounds the complete guarded chunk.
All controls can be tested in $`O(\log(t+2))`$ T-depth
and $`2^{O(t+1)}`$ cost on the small support. No T cost scales
with m inside a chunk.

For the initial constant seed splits of the paired source,
guard the actual native Clifford-conjugated T/T-dagger
representatives. Their literal phase has eigenvalues 1 and an
eighth root after diagonalization, so only the corresponding
guarded eighth-root phase is needed. A projective
$`e^{i\pi P/8}`$ representative with an uncompensated
sixteenth-root scalar is unnecessary. Its native common scalar
on success cancels in the final source conjugation.

Compose the guarded seed and chunk factors in one key frame before
changing coordinates. Every factor returns the whole key bank and
all borrowed work. The exact complete word has operator

```math
 U_1=P_1\otimes I+(I-P_1)\otimes U_*,
 \qquad U_*X_0U_*^\dagger=M_{\rm old}\otimes I_{\rm spare},
```

where $`M_{\rm old}`$ is the original unpaired or paired
source reflection. There may be a fixed common scalar in $`U_*`$,
which cancels in its conjugation. Outside-axis transport proves
the long geometric coefficients by multiplying each terminal
coefficient by the next normalized list. Keep both paths for
the two-tail source. Key selection and the small private banks
are reused at every chunk; no full m-bit payload is copied.

## 3. Exact complete source and simultaneous resource count

For the first frame, conjugate the success-guarded
$`e^{i\phi X_0}`$ center by $`U_1`$ and its actual reverse.
For the second frame, surround the same complete word by
$`S_{\rm key}`$ and its actual reverse. These are exact
$`A_i(\phi)=P_i\otimes I+(I-P_i)\otimes e^{i\phi M_{\rm old}}`$
with all extra work returned. Use the unchanged R100 mirror
$`C_1=P_1\otimes I+(I-P_1)\otimes\Gamma_{aa}`$ on the primary
payload. The single long fixed Pauli in this mirror costs O(m)
Cliffords for axis diagonalization and its inverse, and only one
guarded short Pauli phase. It anticommutes with $`M_{\rm old}`$.
On the two commuting bad projectors with $`P_1P_2=0`$,

```math
 A_1(\pi/2) A_2(\pi/4) C_1 A_2(\pi/4) C_1^\dagger
   = iM_{\rm old}\otimes I .
```

On first success the second rotation and its mirrored version
cancel; on first failure both quarter-angle second rotations act.
The hardwired key permutation occurs only at the boundaries
of complete frame-2 visits, never at a chunk boundary.

For a complete original controlled source with stable F,
add F only at the whole success-guarded $`X_0`$ centers.
When F is false each actual source sandwich and inverse is
literal identity, and the two mirrors cancel. On F true correct
the known scalar with $`(-i)^F`$ using the exact R100 predicate
derivative. A conjunction F on p stable external literals is
one returned XOR oracle in the short predicate grammar. Its
depth is $`O(\log(p+2))`$ and it uses a separate returned
dirty lane, all paid. Its controls are outside the physical
core and the key/payload banks, as at the original complete
Majorana center. The old conjugate-signal convention can use
the same exact paired Pauli sign-mask Clifford as R100.
Every inverse is the actual inverse of the emitted complete word.
The effect is exactly the original completed operator tensored
with identity, on all arbitrary data, work and reference inputs.

There are $`O(1+m/t)`$ chunks. One physical short payload bank,
one key bank, and all their private work are reused at every
chunk. The constituent T count, Clifford count, native T-depth,
and extra dirty peak for one completed source, including actual
reverse, are bounded by

```math
\begin{aligned}
 T_{\rm src}&=O((1+m/t)2^{C t}+P(n)),\\
 G_{\rm src}&=O((1+m/t)(m+2^{C t})+P(n)),\\
 D_{\rm src}&=O((1+m/t)\log_2(t+2)+2^{C t}+P(n)),\\
 w_{\rm src}&=O(2^{C t}+P(n)).
\end{aligned}
```

The constants absorb the fixed number of physical centers,
mirrors and actual reverses. Here P denotes fixed polynomials
for the inherited stable guard reservations; the original
macro core is excluded from the last bound. The key permutation and
its reverse are included. No width factor depends on L through
tensor repetition or output accuracy. This is essential at
$`b=\Theta(L)`$.

## 4. Splice into the direct-prefix high-precision compiler

Choose an absolute c that overlaps R102, then a fixed sufficiently
small positive group-height fraction $`c_h\le1/64`$. For
$`L\ge2^{cn}`$ use the R96 large-precision assembly without an
opening contraction, grouping all depths above the original
three bottom bits at heights
$`h\le\max\{1,\lfloor c_h\log_2(L+2)\rfloor\}`$.
There are a constant number of groups because
$`n=O(\log(L+2))`$. The
[direct-prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md#8-absorption-on-the-intended-branch)
explicitly admits any remaining height $`O(\log L)`$.
Its local high-precision condition
$`\kappa>2^h/h^2`$ and its literal width proof apply here,
also when $`L>N/n^2`$ where the earlier global assembly
could use R79 instead. For bounded local heights use its
inherited directly prefixed R72 fallback. Choose the same
absolute accuracy padding of R96, with the constant number
of groups included.

Keep the finite magnitude/sign reduction, original global
addresses, separately saved H, original conditional d,z,w
from the protected suffix, the same actual mixed gauges and
local coarse program, and the retained R72 smaller-error
bottom. Every actual fine mask and its inverse uses the
unchanged complete R91 query. In particular the R72 bottom and
final complex R72 diagonal also use R91 for their complete
packed X/Z masks; older R72 query schedules with additive
output-width depth are not needed at $`b=\Omega(L)`$ here.
Use the present exact substitution only at completed geometric
source macros. Their full operator equivalence preserves every
rejected-sector action, replay, lower-suffix guard and inherited
error budget, without a new precision slack for promotion.
An exact center plus actual reverses suffices even after earlier
leakage.

The original direct-prefix skeleton and the bounded bottom and
phase have at most $`P_0(n)`$ completed source positions, including
actual inverses, for a fixed polynomial. Their total core
widths sum to $`O(L+P_0(n))=O(L)`$ in this branch.
Therefore their width squares sum to $`O(L^2)`$.
The fine-mask volume and square-root sums are still
$`O(NL)`$ and $`O(\sqrt{NL})`$ over the retained original
prefixes. The old nonquery coarse peak is below $`b/4`$ beyond
one cutoff at the literal threshold, and an exclusive R91
pool of $`\lfloor b/2\rfloor`$ stays disjoint. The unchanged
coarse, guard, saved predicates and indicator depth, apart from
the charged $`NL/b^2`$ volume term, is at most
$`O(L^{1/8}\log(n+2)+P_0(n))`$. Choosing small $`c_h`$
leaves a strict exponent reserve after every fixed polynomial
overhead. These estimates use existing elementary count bounds
and do not require the old global R96 cutoff at $`N/n^2`$.

Fix all native exponents before choosing $`c_0`$. Since
$`L\ge2^{cn}`$, taking $`c_0`$ sufficiently small gives
$`2^{Ct}=O(\min\{S,L^{1/4}\})`$ for all parameters, with
an absolute factor for bounded t. The extra work is $`o(L)`$.
Reuse it across fixed macros and reserve its peak below $`b/4`$,
next to the old precision base and exclusive query pool.

The added T count is $`O(L2^{Ct}+P_0(n)^2)=O(LS+L)`$.
Parity coordinate changes cost $`O(L^2+P_0(n)^2)`$ Cliffords,
and the remaining source Cliffords fit $`O(LS+L)`$. Both fit
$`O(NL)`$. The transported chunks use $`O(L\log(t+2)/t)`$
depth. The additive key permutations, guards and coarse terms cost
$`O(P_0(n)L^{1/4}+P_0(n)^2+L^{1/8}\log(n+2))`$.
This is $`o(L/t)`$ since $`L\ge2^{cn}`$ and $`t=O(n)`$.
The unchanged R91 volume term supplies $`O(NL/b^2)`$ separately.
All bounds apply on the same emitted word.

For bounded n or bounded t use the original R96 circuit with
its literal allocation and full error. In these cases
$`L\log(t+2)/t=\Theta(L)`$. For the complex frame use the real
part and phase at L+1 as in R96. When L=N the shift would cross
the stated window, and R96 at L covers that bounded-S case.
Otherwise L+1 is still in the admitted range and comparable to L.
The complex literal threshold matches the shifted real threshold.
The full error follows from the original initialized-isometry
telescoping, including arbitrary actual leakage and every
physical reverse. Reversal in the positive-S alphabet changes
the Clifford count by at most a constant.

The [bounded key-cover fixture](../../tests/test_tensor_sector_geometric_chunks.py)
checks suffix packing, both key frames, two geometric chunks, the literal guard
in both sectors, and the actual adjoint on the full represented space. Its
test of the frame algebra complements the separate selected-slot fixture.
The gate schedules and simultaneous asymptotic bounds follow from the
analytic construction above.
