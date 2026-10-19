# Revision: native sharing and table-specific indicator probes

[Next revision: clean-source promotion](CLEAN_PROMOTION_SESSION_RESULT.md) ·
[Resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous checkpoint](JOINT_INDICATOR_SESSION_RESULT.md) ·
[Current R98](QUERY_DEPTH_FRAME_TRANSFER.md) ·
[Constant-accuracy extension](APPROXIMATE_QUERY_FRAME_TRANSFER.md)

This continuation fetched main and verified the clean PR #196 commit
`55ae14ea5da2d307194ba95184846b69f7783998` before changing the
repository. The main research remained the complete approximate dirty
indicator. A focused parallel scout examined the original physical
middle pairs left outside the retained no-pair theorem. The session
establishes a native sharing requirement, a paid table-specific interface,
and a sharper diagnosis of two candidate mechanisms. The complete-frame
resource frontier does not move.

## 1. Destination, theorems, components and representation limits

| Role | Current statement |
| --- | --- |
| Destination | Optimal general T-depth for complete prescribed Hopf frames, with all resource and error requirements on the same native circuit |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ |
| Literal domain | $`n\ge1,L\ge6`$, one clean; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error at most $`2^{-L}`$ |
| Selected result | R79 remains Result D; the settled count frontier and its allocations are unchanged |
| Constructive transfer | Current restored-address R98 and its approximate extension need only constant source accuracy at each fixed final precision |
| New unrestricted necessary condition | Every bank observable of a zero-clean approximate indicator requires at least p active native T occurrences at full-operator error below one half |
| New restricted count consequence | Separately emitted b completed disjoint bank blocks cost at least bp T gates; independent rows cost at least mp |
| Useful positive components | Exact selective-inversion dirty echo, logarithmic analog pulse width, syndrome-gauge classification and paid row-query/whole-query transfer |
| Remaining lower scope | Original physical middle pairs, including the new occupied-code insertion question, and higher-layer asymptotics |

The new incidence condition is valid for arbitrary native Clifford and
address/bank mixing. Converting it to a large count lower bound requires
an additional bound on how many observables one T occurrence affects.
No such bound is imposed on the unrestricted model. The completed-block
theorem instead has an explicit representation premise. Neither is a
growing T-depth lower bound.

## 2. The full primitive and the changed mathematical uncertainty

For p address bits and $`m=2^p`$, the target remains

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle.
```

At every fixed requested tolerance tau, seek zero clean work,
$`T,b_{\mathrm{dirty}}=O_\tau(m)`$,
$`G=O_\tau(m^2)`$, and $`D_T=o(\log(p+2))`$.
Accuracy is literal full operator norm on arbitrary address, output,
dirty work and reference inputs. The actual inverse and all return and
leakage errors are part of the contract.

The earlier filter, Gaussian-conjugator, phase-signal and coherent-frame
limits were read before discovery. Their hypotheses do not prohibit
noncommuting joint address/bank words. The geometric continuation
suggested selective inversion as a way to return a scalar on both dirty
signal states. That mechanism does improve the analytic bank-error
profile: exponential detuning suppression pays the binomial row
multiplicity with logarithmic pulse width. Independent-row precision
alone therefore does not force the earlier exponentially large analog
width.

The more decisive uncertainty changed at native compilation. The
[sharing theorem](NATIVE_INDICATOR_SHARING.md) shows that an optimally
compiled complete row still needs at least p T gates when the product
of rows has full-bank error below one half. It assumes no individual
row-accuracy bound and allows internal address motion and arbitrary
row work. Improving the isolated pulse compiler cannot meet linear
bank count. The next mechanism must share native action across growing
bank sets, rather than merely sharpen each independent pulse.

## 3. Strongest new result and its exact scope

Propagate $`Z_{Y_u}`$ through the actual native indicator word. Each
active T occurrence can add at most one direction to its affine binary
Pauli support; a Clifford maps that support linearly. The target image
is

```math
Z_{Y_u}(-1)^{[x=u]}.
```

For $`p\ge2`$, this target is neither even nor odd under any nonzero
address translation. A support of dimension below p would force one
such signed translation symmetry. Its distance two from the target
contradicts full source error below one half. Hence each of the m
observables needs at least p active occurrences, and

```math
\sum_{\text{T occurrences }g}\iota_g\ge mp,
```

where $`\iota_g`$ counts the bank observables active at g. Thus
$`T\le C_\tau m`$ requires average incidence at least
$`p/C_\tau`$. One native gate may have large incidence; this is why
the result does not lower-bound unrestricted depth.

If the emitted word factors into b complete address-preserving blocks
on disjoint bank sets and private work, each block can be isolated by
one final bank observable. Its own native word then costs at least p T
gates. This proves $`T_{\mathrm{emitted}}\ge bp`$. With maximum
block size k, the lower bound is $`p\lceil m/k\rceil`$; even
growing blocks with $`k=o(p)`$ still fail linear count. Joint resynthesis
and address/bank motion across these boundaries are outside this
factorization premise. An exact seven-T joint p=2 fixture already lies
below the eight-T separate-row bound, illustrating the distinction.

## 4. Complete candidate recurrences and first unpaid operations

### Selective inversion

The [pulse proof](SELECTIVE_INVERSION_DIRTY_INDICATOR.md) establishes
for every SU(2) matrix with entries a,b the actual-adjoint echo

```math
E(U)=ZU^\dagger ZU,
\qquad \|E(U)-I\|=2|b|,
\qquad \|E(U)+I\|=2|a|.
```

Thus inversion gives minus identity on the entire dirty signal. The
verified Rosen–Zener amplitude supplies exponential detuning suppression
for $`d_u=|x\oplus u|+1-Y_u`$. The full bank has an analytic
constant-error construction with pulse-width parameter
$`\kappa=\Theta_\tau(\log(p+2))`$.

Finite truncation and actual digital synthesis are also paid. The
explicit first-order construction has

```math
J=\Theta_\tau(mp^3\log(p+2)),\qquad
M=\Theta_\tau(m^2p^4\log(p+2)),
```

with M single-qubit rotation occurrences, requested precision
$`\epsilon=\tau/(4M)`$, T-count at most $`M t_\epsilon`$,
Clifford count at most $`M g_\epsilon+O(Jm(p+1)+m)`$,
depth at most $`2J(p+3)d_\epsilon`$, and m dirty signals.
All four error sources and the actual inverse are included. These are
construction bounds, not lower bounds on joint native simulation.
The independent-row theorem supplies the separate architecture stop.
The first missing operation is cross-row native synthesis with the
linear count envelope, not an uncharged analog evolution.

### Table-specific gauges

The [table-specific theorem](TABLE_SPECIFIC_INDICATOR_INTERFACES.md)
revisits whether the standalone indicator is stronger than a frame
position needs. For a fixed binary table matrix
$`F=M^{\mathsf T}N`$ of rank r, the bilinear probe sees only the
syndrome MA. Its maximal source gauge can mix the bank inside each
syndrome fiber and arbitrarily mix address and private work there.
This is larger than the universal gauge when r is deficient.

The exact paid replacement is the complete r-bit row query
$`a(x)=Me_x`$. Two source calls extract it, with all old bank wires
counted as dirty work; a paid linear bank coordinate change gives the
converse. Direct probe accuracy transfers without a rank factor.
For full row rank, the same reduction recovers the entire indicator.

A safe nested two-source word uses four left-source calls, two
right-source calls and four complete bilinear words. It remains correct
when the source gauges temporarily change their addresses, unlike an
unjustified gate-by-gate substitution. The multi-output construction
uses the joint syndrome across all actual matrices and keeps the source
calls outside all batches and traversals. The proof home charges its
conditional R98 transfer.

Low rank is not a general hypothesis of the actual sign and program
tables. Even a one-output table can have an invertible balanced matrix;
the identity matrix is a rank witness, not a hard-query example. Thus
table specificity identifies a genuine restricted weakening, but supplies
no uniform faster source. The first unpaid operation remains a native
compiler for these row queries, or a complete joint table word which
does not expose either universal one-bank probe.

## 5. The focused original-middle-pair scout

The retained PR #196 degree theorem closes growing total displacement
in the actual no-middle-pair reduction. Its logical pair permutations
are not the original physical clean-transition pairs.

The [new insertion note](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) starts
from those actual signed physical pairs. A literal Clifford correction
reorients their factors into an insertion F commuting with the clean
image and every remaining middle axis. The residual clean comparison
has no physical pairs, and every endpoint list remains independent.
No sector quotient or assumed independence after restriction is used.

For s original pairs, F has original Pauli support span at most 2s.
Pulling it through an endpoint retains a sum of at most $`4^s`$
actual Hermitian Cliffords, but not a bounded Pauli support span. A
native one-pair family with growing endpoint rank demonstrates that
support growth. It is explicitly not a Boolean-target counterexample.

Deleting F preserves clean return but has not been shown to preserve
the required Boolean target. Applying the no-pair degree theorem after
that deletion is therefore unsupported. The existing full-space finite
insertion theorem also does not supply the needed occupied-code
comparison. The next lower question is stability under this actual
finite physical insertion at fixed s. Growing original pair count and
the higher-layer asymptotic question remain separate open cases.
No new one-clean degree closure is claimed in this continuation.

## 6. Frontier decision and next research question

R96 remains the strongest unconditional general upper bound. At fixed
accuracy and square-root dirty width, the constant-versus-logarithmic
general depth interval is unchanged. The exact-rounding floor, selected
R79, all literal allocations and the constant-accuracy R98 transfer
retain their scopes.

Stop independent row synthesis as a route to the requested linear T
count. Continue with a complete two-bank table word which shares native
operations across growing bank sets and need not supply each individual
universal bilinear probe. The concrete question is whether such a
noncommuting joint word can replace the complete R91 four-corner query
with sublogarithmic address T-depth while retaining its count, width,
Clifford and error envelopes. Its first proposed identity must expose
the actual native shared operations and the complete recurrence.

The lower scout should retain the finite Clifford sum of the actual
pair insertion and ask whether literal occupied Booleanity forces a
bounded-degree consequence at fixed physical-pair count. It should not
discard the insertion, quotient away its signs, or infer full-space
probe identities from an initialized-sector comparison.

The documentation-only Hopf-QBP difference was already reconciled in
PR #196 at tracked upstream `99e5ee5`; [SYNC](../../SYNC.md) and all
scientific provenance remain unchanged. This session adds a precise
primary-source attribution for the Rosen–Zener amplitude.

## 7. Verification boundary

The native support, pulse, source-gauge and physical-insertion arguments
received independent analytic review. Exact bounded fixtures check
literal phases, all physical columns, actual inverses, arbitrary work,
the signed-translation gap, and the distinction between separate and
joint bank synthesis. Their finite sizes are not evidence for an
asymptotic compiler or growing-depth lower bound. Required repository
and rendered-presentation checks are recorded on the reviewed tree.
Final manuscript drafting remains separate.
