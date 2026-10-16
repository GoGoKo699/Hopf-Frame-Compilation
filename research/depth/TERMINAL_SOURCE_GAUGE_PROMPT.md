# Five-hour design: compose the terminal with returned borrowed work

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[R89 proof](CAUSAL_COMPACT_COIN_DEPTH.md) ·
[Previous outcome](COIN_CONTINUATION_SESSION_RESULT.md) ·
[Shared-prefix audit](SHARED_PREFIX_QUERY_AUDIT.md)

**Status: designed; not executed.** The scientific baseline is main commit
`eb4bde054ba310dec225e1d42881a5ab3d9d6a01`, integrated in
[PR #152](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/152).
The next execution may use at most five elapsed hours, including review,
validation, tool waits, and integration. Record its start and deadline
when execution begins. Designing this session does not establish any
new compiler. R89 remains the current fixed-accuracy frontier.

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge remain authorized. Preserve the separation between research
and final manuscript drafting.

## 1. Destination and the exact remaining cost

Put $`N=2^n`$ and $`\ell=\log_2(n+2)`$. Begin at
$`\eta=1/64`$ and $`b=\lceil C_b\sqrt N\rceil`$, for a fixed
sufficient constant. Seek a general complete-frame circuit with

```math
T=O(\sqrt N),\qquad G=O(N),\qquad D_T=O(\ell),
```

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta.
```

Here W includes arbitrary node-dependent real pairs, with signs and
zeros, and subsequently the prescribed complex magnitude completion.
The norm includes all logical columns, literal phases, source leakage,
clean-work return, arbitrary dirty inputs, and external references.
Keep one external initialized qubit. No QRAM, measurement, reset,
free coherent evaluator, supplied precision state, or extra clean word
is available. Elementary Clifford counts remain charged. Classical
preprocessing must terminate, but need not be efficient.

After completing the large-width argument, seek the fixed-L extension

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell),
```

at the original literal real and complex dirty thresholds. Do not
silently change those thresholds by substituting a stricter accuracy.
R79 remains the absolute all-precision theorem. Logarithmic general
depth is this session's target; it is not asserted to be optimal.

R89 Section 10 has heights $`m_i`$, local precision sizes $`u_i`$,
and $`m_{i+1}=O_L(\log^3(m_i+2))`$. Its paid stage ledger is

```math
D_i\le C\left[\log_2(n+Cu_i^3+2)+u_i+R_i/b^2\right],
\qquad R_i=K_i^d\overline Q_i.
```

The protected margin already makes the table terms and gate counts
summable. After increasing the fixed stopping cutoff, the local sizes
satisfy $`u_{i+1}\le u_i/2`$, hence $`\sum_i u_i=O_L(\ell)`$.
Reconstruct this elementary inequality at the start. Local routing and
precision costs therefore do not intrinsically carry the iterated
logarithm. The remaining repeated term is the global-address setup
$`\log_2(n+Cu_i^3+2)`$. This is an upper-bound charge in the current
implementation, not a per-query lower bound.

Pay the first R89 stage normally. Its remainder has height
$`M=O_L(\ell^3)`$. Throughout the remaining terminal, the first
$`n-M`$ logical bits form one genuinely immutable quantum prefix P.
Later appended prefix bits belong to the evolving local register Y;
they cannot be treated as classical functions of the initial input.
The M-bit remainder is now logical data, not a permanently clean cache.

## 2. Primary lemma: transport the borrowed-work error through a whole word

The unit of design is the complete ordered terminal, including its
sources, noncommuting local Hadamards, calibration, reflections,
coarse readout, and actual inverses. Intermediate borrowed work may
remain correlated. The final operator must return it.

Write the exact target as
$`W_{\rm tail}=E_s\cdots E_1`$. Seek explicitly represented physical
unitaries $`G_i`$ and implementable stage words $`\widetilde E_i`$
on the same logical, sole clean, and dirty registers, satisfying

```math
\|\widetilde E_iG_{i-1}J_1
 -G_iJ_1(E_i\otimes I_b)\|_{\rm op}\le\epsilon_i.
```

Then the physical word

```math
G_s^\dagger\widetilde E_s\cdots\widetilde E_1G_0
```

has the desired initialized action up to the summed errors. The
representations of every G and transformed stage must be concrete.
Writing a conjugated oracle, choosing an expensive abstract completion,
or defining G to contain the unknown target circuit does not prove a
native construction. Charge both endpoint words and their errors.

The retained source grammar is useful structure: a constant number
of Boolean-phase table positions separated by selected local
Hadamards, together with calibrated anchors and actual inverse
sandwiches. Its source completions and local matchings can be chosen
jointly if their required column/flag action is preserved. Begin with
this grammar rather than demand a universal reusable Boolean oracle.
Continued Clifford access to P is allowed and counted. Every changed
local address remains quantum.

A sufficient depth target for the remaining terminal is

```math
D_{\rm tail}\le C_L\left[
\ell+\sum_i\bigl(\log_2(m_i+2)\bigr)^c+\sum_iR_i/b^2\right]
```

for any fixed $`c\ge1`$, together with the original summed T/G/width budgets.
For $`M=O_L(\ell^3)`$, the local sum is
$`O_L((\log_2(\ell+2))^c)=O_L(\ell)`$ beyond a fixed cutoff.
Thus a polylogarithmic local correction is acceptable. The task does
not require every local correction to have logarithmic depth.
Persistent banks, output banks, traversal helpers, and source work
must fit simultaneously; stagewise maxima cannot be added selectively.

### Opening tests that can reject a proposed identity

First reconstruct the retained
[noncommuting dirty-cache witness](JOINT_DEPTH_SESSION_RESULT.md#5-a-generic-dirty-cache-correction-already-contains-a-lookup)
and [changed-address residue](SHARED_PREFIX_QUERY_AUDIT.md#1-the-exact-residue-of-a-retained-prefix-short-echo).
Plain retained one-hot tables and a single outer difference already
fail these tests. Their correction can itself contain an arbitrary
lookup. Reopening that interface requires a new complete identity.

For a proposed transported error or common output gauge, use two
different phase tables separated by a local Hadamard, a changed local
address, and the actual inverse source. Check every logical and dirty
basis input with literal phase, then a coherent external-reference
input. Extend the surviving identity to the anchor and reflection
positions. A commuting XOR example or root-column match is insufficient.

The [prefix-blind cache obstruction](SHARED_PREFIX_QUERY_AUDIT.md#2-a-dirty-only-cache-cannot-carry-an-exact-prefix-blind-reader)
does not cover every allowed representation here. State precisely
which assumption a candidate changes: continuing access to P,
encoding into the payload, or another explicitly priced interface.
Avoid introducing a clean copy of the arbitrary dirty baseline.

The first concrete ansatz is a private linear change of dirty-array
coordinates together with quadratic dirty-array phases whose
coefficients depend on short local words. Track the four bilinear
phase corners through the Hadamard and source pair. In particular, a
pair of retained prefix transformations with matrices $`A_P,B_P`$
changes a bilinear readout matrix to $`A_P^{\mathsf T}DB_P`$.
Its P-dependent coefficients
are not free. Either represent their action within the proposed native
budget or reject this ansatz. Closure cannot be assumed after a
noncommuting local operation.

## 3. Independent route: a compact-coin completion for every dirty label

Explore a full-basis version of the compact embedding. For a fixed P,
let E embed a D-dimensional local input into a DK-dimensional physical
atom space, where K is the dimension of an actual borrowed register.
Seek native-compatible unitaries $`G_y`$ satisfying

```math
E^\dagger G_y^\dagger G_zE=\delta_{yz}I_D,\qquad
\sum_yG_yEE^\dagger G_y^\dagger=I_{DK}.
```

The maps

```math
B|x,y\rangle=G_yE|x\rangle,\qquad
A|j,y\rangle=G_yEW_{\rm tail}|j\rangle
```

would then be full unitaries with
$`B^\dagger A=W_{\rm tail}\otimes I_K`$, including coherent dirty
labels. This encodes the unknown borrowed label in orthogonal sectors
and consumes it on return. It is a sufficient algebraic contract, not
a construction: implementing A cheaply is part of the missing lemma.
Additional helpers, the sole initialized bit, and prefix dependence
remain charged.

The first gate is dimension and support. A gauge restricted to each
original fiber needs enough room for K orthogonal images in that
fiber. Small or uneven fibers may immediately reject this subclass.
Padding, redistribution, or cross-fiber mixing must be priced and
must preserve or replace the required complete source/flag identity.
An arbitrary abstract orthogonal completion provides no native bound.

The small binary alias example with probabilities one quarter and
three quarters already has positive fiber sizes one and five in an
eight-dimensional padded space. With $`D=2,K=4`$, a fiber-preserving
completion fails the support test. Reconstruct that example first.
A candidate may use signed weight-class-preserving permutations that
cross fibers, provided both source words remain explicitly compilable.
An arbitrary lookup over every syndrome value can add an unpaid
factor K to the table volume.

Test two nonidentical blocks with interior nonbalanced probabilities
and prefix-dependent partitions. Check cross-syndrome orthogonality,
the entire marker family, literal signs, and the full relative word,
not only each fixed-y root state. A successful small fixture motivates
the symbolic general identity; it does not establish its catalogue size.
Construct A independently from the proposed source and marker formulas;
defining it as BW and checking the product would only test a tautology.

Keep this route independent of the first. If a whole-basis completion
works, it may remove the need for recursive clean suffixes altogether.
If it fails, record the exact subclass that fails, not a general
impossibility claim.

## 4. Proof and resource gates

Before promoting a surviving route to a theorem, it must supply:

1. A literal full-word identity or controlled approximation, including
   changed quantum selectors and actual inverses.
2. A register map with initialized, dirty, logical, and inactive roles
   at every live boundary. A protected suffix is clean only where its
   saved predicate proves it.
3. Separate recurrences for T-count, elementary Clifford count,
   T-depth, peak simultaneous width, and full-isometry error.
4. Charged global table addresses, short labels, source precision,
   both endpoint corrections, and all returns. A joint source cannot
   retain an unpaid old-input/history word or double its address exponent.
5. A complete terminal composition and the original-width finite/narrow
   fallback. Add signed/zero data and complex leaf phases through their
   proved paid interfaces.

Do not require exact cleanup from an approximate candidate beyond its
full initialized-isometry error contract. Exact borrowed primitives
keep their exact guarantees. Fixed-L constants and dimension cutoffs
must be stated; $`\log^*n`$ is never treated as constant.

Small exploratory checks may reject or clarify a specific identity
before these theorem gates are complete. Broad finite validation follows
a surviving native identity and symbolic resource ledger.

## 5. Five-hour execution and decision rules

| Elapsed time | Work and required decision |
|---|---|
| 0–25 minutes | Confirm baseline, reconstruct the residual ledger, identify the immutable prefix, and reproduce the old failure algebra |
| 25–120 minutes | Run the transported-gauge and full-basis-completion lanes independently; a skeptical reviewer checks the two-stage identities and allocation |
| 120–180 minutes | Select the strongest surviving mechanism; derive its native word and summed resources; stop or change a candidate that merely renames an expensive query |
| 180–235 minutes | Complete the full-frame splice, precision/error proof, literal widths, and targeted finite checks |
| 235–300 minutes | Independent proof audit, repository integration, required validation and rendered-document checks, CI, and merge if the gates pass |

These are checkpoints, not reasons to halt a promising derivation at an
arbitrary minute. Reserve the final validation window and respect the
overall ceiling. Use parallel derivation and skeptical review; assign
separate ownership before concurrent edits.

A full success is the logarithmic general theorem with all simultaneous
budgets. A strict general improvement below R89 is also a valid theorem
outcome if its asymptotics are proved. A narrower native composition lemma
is useful only with its actual hypotheses and an explicit remaining
composition gap. Exact failures should be preserved compactly when they
prevent repeating a route.

More fixtures, an unpriced source identity, or a changed constant in the
same bound is not a new general theorem. If neither mechanism closes,
leave R89 unchanged, record the precise failed interfaces and the best
remaining question, and decide whether to reformulate or pause the
candidate. The continuing destination is optimal general depth.
