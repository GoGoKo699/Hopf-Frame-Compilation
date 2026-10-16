# General-magnitude recovery: session result

[Executed design](GENERAL_MAGNITUDE_RECOVERY_PROMPT.md) · [R86 proof](INTERVAL_DILATION_DEPTH.md) · [Selected routing](SELECTED_DIRTY_ROUTING.md) · [Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer)

**Outcome: a complete fixed-accuracy construction for arbitrary
node-dependent magnitudes with sub-square-logarithmic T-depth.**
The proof is [R86, interval dilation](INTERVAL_DILATION_DEPTH.md).
It removes R85's shared unordered magnitude promise. Optimal depth
and an absolute all-precision improvement remain open.

Execution began on 15 October 2026 at 23:41:06 UTC from main commit
`97c18da15c8ec578562c25f95987277fa31db1b4`, on branch
`research/general-magnitude-session`. The authorized five-hour ceiling
was 04:41:06 UTC on 16 October. Construction, independent audits,
validation, and integration are distinct checkpoints; the ceiling is
not a claim of five hours spent.

## 1. The new bound

For every fixed integer $`L\ge6`$, put $`N=2^n`$,
$`\eta=2^{-L}`$, and $`\ell=\log_2(n+2)`$. The same one-clean
circuit has simultaneous bounds

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad
G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell^{3/2}).
```

The original dirty thresholds remain $`17(L+n+7)`$ for real frames
and $`17(L+n+8)`$ for prescribed complex magnitude completions.
Arbitrary node-dependent pairs, signs, zeros, and supplied leaf phases
are included. The norm covers all columns, literal phases, leakage and
work return, arbitrary dirty inputs, and external references.

At square-root dirty width, the general fixed-accuracy upper bound
improves from square-logarithmic depth to $`O_L(\ell^{3/2})`$.
R79 remains the absolute all-precision theorem. R85's smaller
$`O_L(\ell\log(\ell+2))`$ bound still applies to its promised family.
No optimal-depth claim, efficient classical preprocessing, scalable
native emitter, or final manuscript is asserted.

## 2. The decisive change of representation

The opening design tried to undo a general local flag using product
rows or a faster recursive contraction. Preparing a row did not erase
its old address, and shrinking one fine interval retained a coarse
recovery obligation. The successful reformulation changed the local
space before doing that conversion.

Round the entire tree once to strictly positive dyadic branch
probabilities with denominator $`2^\tau`$,
$`\tau=O_L(\ell)`$. A full-frame layer bound controls the rounding
error, including markers below originally zero branches. A q-level
local tree then has integer leaf masses in a uniform atom space of
$`r=q\tau`$ bits. Assign consecutive atom intervals to leaves.

The isometry E sends a leaf to the normalized uniform vector on its
interval. The atom itself uniquely identifies its leaf, so an ordinary
classical membership XOR erases the old leaf address after preparation.
Every weighted scaling vector becomes a normalized interval indicator;
every marker becomes an interval Haar wavelet.

Let B be the ordinary balanced Haar frame on the atom space. The
rectangular overlap $`K=B^\dagger EW`$ is an isometry. The key new
estimate is

```math
\bigl\|\,|K|\,\bigr\|\le C+\sqrt{3(\tau+2)}.
```

Its proof separates fine, middle, and coarse dyadic scales. The two
tails have geometrically decaying row energy and bounded column degree.
The middle has at most $`O(\tau)`$ entries per column and row energy
at most one. This gives square-root normalization without a density
promise.

Certified Schur weights turn K into a paid row/column source overlap.
Separate rejection registers, followed by B and the actual inverse of
E, yield a square block for the original local frame. Calibrated
amplification takes $`O(\sqrt\tau)`$ source rounds. The error comparison
is uniform as the success amplitude decreases; every rejected register
belongs to the tested anchor.

Thinner geometric bins expose $`O(q\tau)`$ clean logical space while
keeping only $`O_L(\ell^2)`$ bins. Exact selected routing moves that
space into fixed banks once around the entire local recovery word.
Each large source round costs $`O_L(\ell)`$ depth, giving the new
$`O_L(\ell^{3/2})`$ total. Coarse recovery and the terminal remain
$`O_L(\ell\log(\ell+2))`$.

## 3. Independent audits and corrections

The interval norm, physical square compression, uniform amplification
error, signs and zeros, selected routing, and outer theorem were reviewed
separately. The final proof retains the following material safeguards.

- Finite rounding uses rational approximants and deterministic rational
  ties. It does not test an original computable input for exact zero,
  sign, or equality at a rounding boundary.
- The fine-scale telescoping argument permits omitted positive terms;
  qualifying minimum-child lengths need not form an initial path segment.
- Coarse intervals with nonzero overlap lie near a dyadic wavelet jump.
  Localized Bessel bounds, rather than a false coarse-column Schur sum,
  control that tail.
- The square wrapper uses a separate output register. Rejected old
  inputs and source flags cannot re-enter the accepted anchor through
  the later inverse embedding.
- The amplification error constant is independent of its growing
  repetition count; the source tolerance pays that count explicitly.
- Variable-size source and balanced words preserve inactive dirty
  padding. The success predicate tests only selected initialized fields,
  never parked dirty values or the maximum physical bank.
- Per-bin tables read their own relevant inputs. Their volumes sum;
  a fully padded physical address is not assigned the smaller exponent.
  Parallel outputs have private non-Clifford supports.
- The free donor coefficient is fixed before enlarging the protected
  suffix. The residual $`M-m_0`$ suppresses every polynomial count
  and dirty-width overhead.
- Narrow-width and finite-dimension fallbacks use the original total
  accuracy, preserving the literal one-clean reservations.

These are internal mathematical audits. They do not replace future
external technical review.

## 4. Other explored routes and their scope

For the ordinary prefix-controlled rotation transporter C, its rows
factor as

```math
C[y,x]=\prod_{j=1}^q R(\theta_j(u y_{\lt j}))[y_j,x_j].
```

Thus a selected product-row source is valid, but it leaves the old y.
Changing that address before reversing its source is not a valid
uncomputation. In the Y eigenbasis, $`C(SH)^{\otimes q}`$ is a
complex Hadamard matrix with triangular phases. Standard Fourier
cleanup leaves corrections $`CY^zC^\dagger`$. Their root rotations
cancel, but general relative transports of height q minus one remain.
This is a useful normal form, not a fast recursive compiler.

A high-arity cut can expose a much shorter active local interval with
one descendant inverse. Its coarse label still has many values. Keeping
the old threshold search therefore retains square-logarithmic depth.
Copying computational prefixes of a scaling vector before unloading
entangles those records with the very state to be unprepared; they
cannot be silently treated as fixed clean history.

The shallow-query lane distinguished a fast single predicate from a
complete exponentially large dirty indicator. A useful scoped check is
the following. Suppose a coherent seed chooses a bank-toggle vector
$`v_s(x)`$, followed by seed unpreparation, and must approximate the
exact one-hot XOR $`e_x`$ with norm error epsilon on every arbitrary
dirty input. Testing dirty Fourier characters and averaging over them
gives

```math
\Pr_s[v_s(x)\ne e_x]\le\epsilon^2/2.
```

Indeed a Fourier character indexed by y has squared error
$`4\Pr_s[y\cdot(v_s(x)\oplus e_x)=1]`$; averaging y gives twice
the probability of any wrong toggle. Per-row fingerprint success is
therefore insufficient for this architecture. This is not an
unrestricted lower bound. No new imported oracle theorem is used by R86.

The new [selected-routing lemma](SELECTED_DIRTY_ROUTING.md) is retained
as a reusable supporting component. Its proof configures a Beneš network
classically and implements every switch layer by exact arbitrary-dirty
involution echoes. It needs no initialized routing bits.

## 5. Verification and integration

The seven [new bounded fixtures](../../tests/test_interval_dilation_depth.py)
cover heterogeneous integer intervals and all frame columns, the
rectangular overlap, membership erasure and actual inverse, separate
rejection histories, literal calibrated amplification, leakage with
actual inverses, reached-band erasure, and occupied inactive padding.
They do not establish the asymptotic bounds or constitute a native
emitter. Existing scientific receipts are preserved.

Two independent complete-proof reviews found no remaining substantive
gap after the corrections above. Local validation passed all 952 tests,
including the seven new fixtures. The reviewer walkthrough, native real
and complex examples, residual-QBP example, all four unchanged
fault-tolerant receipt suites, offline upstream check, compile check,
and both resource ledgers also passed. The presentation check rendered
all 151 Markdown pages at desktop and narrow widths, including MathJax
and native MathML, without failures.

The reviewed pull request records the final tree's presentation and
required remote checks, merge state, and actual elapsed session time.
These checks support the implementation and reading interfaces; the
asymptotic theorem rests on the analytic proof and its stated premises.

## 6. Precise next question

Can a general frame match the $`O_L(\ell\log(\ell+2))`$ depth already
proved for R85's promised family, while retaining the same clean,
count, and complete-operator contract? In R86 the new leading cost is
the $`O(\sqrt\tau)`$ amplification of the interval overlap. A
qualifying next mechanism must reduce that normalization cost, process
the scale decomposition more cheaply, or change its source interface
with a complete return proof.

Logarithmic or optimal unrestricted depth is a further question. The
general fixed-accuracy sub-square-logarithmic milestone is now separate
from those remaining goals.
