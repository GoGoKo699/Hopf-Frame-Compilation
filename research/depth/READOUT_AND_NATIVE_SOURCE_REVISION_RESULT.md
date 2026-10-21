# Revision: native source scope and complete projector readout

This continuation fetched main commit
`3ac6ee925f1257aafaefc40cef201427767aea8d`, tree
`6ebf13b038c041643505d867ebbdcc4a3c42e7b5`, merged as PR #216.
The [preceding revision](EXPANDER_AND_SOURCE_SCOPE_REVISION_RESULT.md)
paid spectral promotion and a constant-depth approximate scalar zero
reflection. This revision checks whether those components actually
compose into the complete indicator, then changes the readout target.
It does not claim an improved unrestricted resource frontier.

## 1. Destination, theorems, components and approach limitations

| Role | Assessment after revision |
| --- | --- |
| Destination | Optimal general native T-depth for complete frames, with simultaneous count, width, precision and return budgets |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md) is unchanged; [R79](UNIFORM_SUBTREE_DEPTH.md) remains selected Result D |
| Complete constructive target | A zero-clean p-bit indicator at every fixed tolerance, with linear T-count and dirty width, quadratic Clifford count, and depth sublogarithmic in p |
| Paid components | Scalar expander zero tests, spectral clean-source promotion, polynomial-overhead root fusion, quadratic forward kernels, and actual-inverse readout recurrences |
| New source-specific theorem | Accurate complete native read-only two-band compression requires linearly many clean inputs, even without flatness or exact return |
| New representation theorem | The exact maximum normalized Clifford correlation with the addressed observable is at most the address length plus one divided by the bank size; normalized Clifford averages cannot uniformly separate its whole eigenspaces |
| New positive conditional interface | A one-clean source whose accepted Gram operator encodes the selected half-space gives a complete indicator through six-call recovery and paid serial promotion; the accepted unitary fiber need not be phase-aligned |
| New Gram-source screen | For p at least three, one-clean two-band sources cannot meet Gram error at most $`1/256`$, even with unrestricted intermediate address motion. At least three source bands or another clean allocation is required |
| Focused physical-pair component | Singular affected-silent blocks reduce to original native sublists bounded by pairing nullity; small nullity gives a full Clifford removal |
| Remaining limitations | No normalized shallow complete source has been constructed. Moving-address sources, genuinely wide non-Clifford readout and unrestricted native depth remain outside the new obstructions |

The full target is

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle,
\qquad |Y|=m=2^p.
```

Every statement retains arbitrary coherent address, output, dirty-work
and reference inputs. Literal phase, actual inverses and all return
error are part of the operator contract. The general R96 law remains

```math
\begin{aligned}
D_T&=O(NL/b^2+L+\log_2(n+2)),\\
T&=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
\end{aligned}
```

with its retained one-clean allocation and literal dirty thresholds.
At fixed accuracy and square-root dirty width, logarithmic general
upper depth versus constant lower depth remains the unresolved interval.

## 2. The first unsupported composition is now explicit

The [new read-only source theorem](READ_ONLY_TWO_BAND_SOURCE_GAP.md)
uses the actual three-Pauli-control normal form of a two-band word.
If each completed Clifford boundary fixes x, restriction to the common
zero clean-X syndrome leaves an address subspace of dimension at least
p minus twice the clean count k. A matrix coefficient on that subspace
has Fourier nuclear-norm sum at most one. The complete bank's independent
output labels force a sum growing as the square root of the subspace
size. Projective cocycles and all scalar phases are retained.

For amplitude alpha and complete compression error epsilon, the proof
gives

```math
1\ge\frac{\alpha M}{\sqrt{3M-2}}-\varepsilon\sqrt M,
\qquad M\ge2^{p-2k}.
```

In particular,

```math
\alpha=1/\sqrt2,\quad\varepsilon/\alpha\le1/8
\quad\Longrightarrow\quad
k\ge\left\lceil\frac{p-3}{2}\right\rceil.
```

No exact Gram flatness, exact clean return or dirty-width bound is
needed. This is much stronger than the general two-band rank--norm
screen, but only under its additional actual-boundary hypotheses.
Diagonal phases on an unchanged address do not escape the proof.

Consequently the current spectral wrapper cannot receive such a
two-band source in its sublogarithmic-clean window. For that emitted
wrapper, the private source-output banks alone would have width at
least two to the three-p-over-two minus a constant. Root fusion cannot
absorb this exponential factor into polynomial overhead.

The unsupported inference was that completed source address return
permits all private source copies to share a read-only native address.
It does not. An address-moving source needs its own paid schedule.
Neither this finding nor the wrapper's upper resource recurrence is
an unrestricted depth lower bound.

## 3. Wider Clifford branches still do not pay complete filtering

The [Clifford lookup theorem](CLIFFORD_LOOKUP_FIXED_SPACE.md) permits
arbitrary address, bank and helper mixing. For p at least two and any
full-space Clifford C, it gives the sharp bound

```math
|\tau(\Omega_p C)|\le b_p:=\frac{p+1}{2^p},
\qquad
\Omega_p=\sum_x|x\rangle\langle x|\otimes Z_{Y_x}\otimes I.
```

Here tau is normalized trace, used only as a witness. A wide
address-controlled Pauli Clifford attains the bound. The proof uses
flat affine Fourier support of quadratic phases and the independent
bank labels; arbitrary affine support constraints and dirty helpers
are included.

Thus a normalized Clifford average whose positive lookup sector has
energy at most eta and whose negative sector has energy at least g
must satisfy

```math
g-\eta\le2b_p.
```

Even the exact fixed-space statement is rigid: a Clifford fixing the
entire positive lookup sector pointwise is identity. A k-clean
Clifford-sum compression of the lookup has a corresponding coefficient
budget lower bound proportional to two to the p-minus-k divided by p.
These results stop a bare-Clifford replacement of the scalar expander,
including wide branches that leave the earlier Schmidt-rank screen.
They do not bound native T-depth: a single wide T layer can already
have a large Clifford coefficient budget.

## 4. Change the readout target, not merely its normalization claim

The [paired central readout](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md)
uses two arbitrary carrier banks of n equal to pm bits. The paired-Pauli
sum is a half-rank partial isometry

```math
B_x=\Pi_{\rm eq}(x)\,\mathrm{SWAP},\qquad
B_x^\dagger B_x=\Pi_{\rm eq}(x)
 =\frac{I+\Omega_x\otimes\Omega_x}{2}.
```

This avoids the thin Bell-code rank obstruction. It does not normalize
the source for free: ordinary Pauli SELECT supplies B divided by
$`2^{n-1}`$, rather than the required constant-amplitude block.
The normalization remains the first unpaid native operation.

A more useful sufficient target is the **Gram operator**. Let U have
one clean signal, J its zero embedding and K its complete accepted
block. It suffices to construct

```math
K^\dagger K=\Pi_{\rm eq}/2.
```

The accepted unitary fiber is unrestricted. Phase-S recovery followed
by a signal-Z echo reads this projector with six U/actual-adjoint calls,
returns the signal exactly, and cancels that fiber. Fixed CNOT and
output-Hadamard conjugations give the complete indicator. The serial
one-clean promotion accepts arbitrary intermediate address motion;
the spectral wrapper's native-sharing premise is not imposed here.

For the stronger approximate block promise
$`\|K-B_x/\sqrt2\|\le\eta\le1/8`$, put
$`e=3\sqrt2\eta+2\eta^2`$. The fully substituted zero-clean ledger is

```math
\begin{aligned}
T_{\rm out}&\le48T_U+48m+8,\\
G_{\rm out}&\le48G_U+16n+98m+152,\\
D_T({\rm out})&\le48d_U+32,\\
b_{\rm out}&\le2n+b_U+2,\\
\varepsilon_{\rm out}&\le32e.
\end{aligned}
```

The complete readout word is exactly Hermitian for every actual U;
its same-word inverse is therefore the actual adjoint, including phase.
At tolerance tau, requesting eta at most tau divided by 144 suffices.
All source quantities are evaluated at that requested precision.
The weaker approximate Gram promise with error delta instead gives
complete error at most $`32\sqrt{24\delta}`$ with the same resources.
It therefore suffices to request $`\delta\le\tau^2/24576`$.
Both estimates include all output and arbitrary-work return errors.

With a qualifying source, the extra carrier width is only polynomial
overhead, so the retained fusion reduction would normalize the counts.
No such source is supplied here. This is a paid alternative interface,
not a claim that restating the projector already changes the frontier.

The [principal-angle audit](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-band-projector-gram-gap)
also stops the most immediate attempt to realize it in two source bands.
The two one-band clean images are Hermitian Cliffords. Their product C
has principal-angle spectrum determined by the accepted Gram operator.
At Gram error at most $`1/256`$, C to the fourth is within one half
of identity; Clifford Pauli discreteness and literal reflection symmetry
force it to be exactly identity. A paid square then supplies a two-band
exactly returned clean action within $`1/64`$ of the lookup sign.
The retained Boolean rigidity floor excludes this for p at least three.

This argument permits arbitrary native address motion and dirty width.
It is specific to one clean signal and two source bands. It does not
follow from the unitary-compression theorem in Section 2, and does not
exclude a three-band or other sublogarithmic Gram source. The next
candidate should leave this closed two-band case, not merely change
its intermediate address convention.

## 5. Keep the original physical-pair scout focused

[Section I.5 of center transport](FOUR_BAND_CENTER_TRANSPORT.md#i5-singular-silent-blocks-have-a-paid-kernel-size-reduction)
removes an artificial nonsingularity restriction. For the actual
affected-silent block, put k equal to its axis count s minus its
restricted pairing rank. Literal common-axis mergers leave at most
$`\min\{s,\lfloor4k/3\rfloor\}`$ axes in each original band.
The retained axes are actual sublists, not an assumed native quotient.

For k at most two the entire block is Clifford. At k equal to three,
only the known rank-four complementary-triple cubic residue can remain.
The general reduction removes twice the number of merged pairs in T
count, never increases native depth, adds quadratic Clifford synthesis
cost and no workspace, and preserves every physical column exactly.

The small guard now explicitly compares eight T gates in two bands
with seven T gates in three bands, including actual encoders and
inverses. The latter saves count but not depth. Neither alternative
supplies the missing Boolean endpoint or a general four-band theorem.
The scout continues with genuinely coupled center transport outside
this silent factorization.

## 6. Next research question and certification boundary

The next constructive question is a normalized full-space Gram source
with at least three native bands, or another complete wide non-Clifford
readout, allowing genuine address motion while paying its native
schedule. A low-budget Clifford average
cannot supply it. A direct native read-only two-band unitary compression
cannot fill the current promotion window. These are reasons to change
the mechanism, not reasons to abandon optimal general depth.

Proofs are analytic, with independent internal audits and small exact
phase/native fixtures. Those audits are not external technical review.
The repository's ten-command validation and rendered-document gate
must pass on the reviewed tree before merge; their execution receipts
belong to the corresponding review record. No manuscript was changed.
