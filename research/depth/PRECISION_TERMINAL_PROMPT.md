# Five-hour design: precision without repeated long-prefix work

[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[R90 proof](DIRTY_CATALOGUE_TERMINAL.md) ·
[R79 uniform theorem](UNIFORM_SUBTREE_DEPTH.md) ·
[Previous outcome](TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)

**Status: executed.** The [session result](PRECISION_TERMINAL_SESSION_RESULT.md)
records R91's exact reader, R92's absolute precision envelope, and R93's
grouped-terminal improvement. Execution began on 16 October 2026 at
06:19:12 UTC, with an 11:19:12 UTC ceiling including integration.
The stronger primary benchmark below remains open. This design's scientific baseline
is main commit `aca1e24289743e06d0580996400f03f232d3b283`, integrated in
[PR #154](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/154).
No proposed target below is claimed merely by this design. The proved
results have separate proof homes linked from the session record.

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Repository modification
and merge remain authorized. Preserve the separation between research
and final manuscript drafting. Optimal general depth remains the
destination; this session tests whether precision can be added without
repeatedly paying for the long quantum address.

## 1. Main target and the contract that must survive

Put $`N=2^n`$ and $`\ell=\log_2(n+2)`$. The primary benchmark is

```math
L=\max\{6,\lceil\ell\rceil\},\qquad
\eta=2^{-L},\qquad b=\lceil C_b\sqrt{NL}\rceil,
```

with a sufficiently large absolute constant $`C_b`$. Seek a complete
arbitrary real-frame compiler with absolute constants and

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad
D_T=O\!\left(\ell\log_2(\ell+2)\right).
```

The stretch target is $`D_T=O(\ell)`$. R79 currently gives
$`O(\ell^2)`$ at this benchmark. The error is inverse polynomial in
the number n of logical qubits, not in the dimension N. R90 gives
logarithmic depth at fixed L; its constants cannot simply be reused
when L grows. Improving this benchmark does not prove optimal depth.

Require one external initialized qubit and

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta.
```

This includes every logical column, literal phases, signs and zeros,
source leakage, clean-work return, arbitrary dirty inputs, and external
references. Use actual inverse circuit words. Exact borrowed primitives
retain their exact return guarantees; approximate constructions must
meet the full initialized-isometry bound. No QRAM, measurement, reset,
supplied precision state, free coherent evaluator, or extra clean word
is available. Elementary Clifford gates remain charged. Classical
preprocessing must terminate but need not be efficient.

After the real benchmark closes, price the prescribed complex phase
completion, including its accuracy shift. Extend the width tradeoff only
if the same resource argument supports it:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
```

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right).
```

The original thresholds are $`b\ge17(L+n+7)`$ for real frames and
$`b\ge17(L+n+8)`$ for complex magnitude frames. A result only at
square-root width, or only for real frames, must be stated at that scope.
Do not enlarge a literal threshold by silently replacing L with L+1.
For composition, establish the real benchmark with uniform constants
under the required additive precision shifts, for example at each
$`k\in\{L,L+1,L+2\}`$. The exact slice $`k=\lceil\ell\rceil`$
alone does not supply its shifted instance.

## 2. Opening audit: expose all four precision costs

Read R90 Sections 5–10, R89's
[uniform local-stage ledger](CAUSAL_COMPACT_COIN_DEPTH.md#10-prefix-selected-recursion-and-iterated-logarithm-depth),
R79 Section 4, and the
[shared binary terminal](SHARED_BINARY_TERMINAL.md).
Reconstruct these costs before proposing a new interface:

| Current mechanism | Charge that a replacement must address |
|---|---|
| R89 protected suffix | A reserve of order $`(L+\log(m+L+2))^3`$ leaves a precision-dependent terminal |
| R90 complete program and catalogue | $`q\le(L+12)2^m`$ program bits and $`2^q`$ catalogue slots |
| R90 native local factors and program reader | Per-slot synthesis $`S=O(2^{3m}(L+m+1))`$; both the rectangular query and XOR slot router have an additive q depth charge |
| R79 shared binary middle | In the relevant nonempty middle regime, $`\Theta(L)`$ ordered stages retain serial global-prefix queries even though source preparation is shared |

For the complete-program query, $`Q=2^{n-m}`$ and $`Qq=O(NL)`$.
With B private dirty query wires, the retained depth bound is

```math
D_{\rm query}=O\!\left(q+\log_2(n+2)+Qq/B^2\right).
```

One query invocation therefore does not imply one logarithmic-depth
operation. Removing catalogue volume alone does not remove q or S.
Listing precision digits separately does not make successive logical
rotations commute.

After one coarse reduction, the remaining height is plausibly bounded
uniformly by $`m=O((L+\ell)^3)`$; prove the precise stopping version
when used. The first n-m logical bits are an immutable quantum prefix P.
The m-bit target X continues to evolve. A prefix-selected version of
R79's early decoder could have
$`O((L+\ell)\log_2(m+2))`$ depth, but that alone does not complete
the frame: its protected middle still has $`\Theta(L)`$ sequential
queries reading P, costing $`O(L\ell)`$ by the existing schedule.
Replacing n by m in that query charge is invalid. This is an obstruction
to that composition, not an unrestricted circuit lower bound.

For every proposed reader, list its immutable address, changing address,
output width, serial appearances, simultaneous work, and inverse. This
ledger is the first decision gate.

## 3. Main experiment: a structured returned program for the whole middle

The question is whether R90's exact cancellation can serve a structured
ordered stream without enumerating every precision-labelled terminal.
Its two disjoint borrowed arrays already transfer a selected commutator
to the target and return both arrays. It does not yet supply a compact
interpreter or a free decoded program.

Seek a complete native word for the precision-dependent middle, allowing
correlations with borrowed work between stages. The word must amortize
long-prefix selection and implement the genuinely ordered target,
including changing local selectors. A constant number of long reads
is a sufficient objective; $`O(\log(L+2))`$ rounds may also meet the
primary benchmark if all other costs fit. An equivalent changed
representation is welcome, but it must name the charge it removes.

Use two constructive lanes with a separate skeptical reviewer:

1. **Factor the catalogue while retaining the full-word cancellation.**
   Start with two disjoint banks and two short program blocks. Replace
   one full-terminal factor by an ordered product of short factors;
   expand the complete borrowed-work word before combining lookups.
   Determine whether internal selected swaps or translations cancel
   across the product, or whether noncommuting factors recreate an
   expensive mask-dependent correction. A scalable representation must
   price every catalogue, shared selector, program bit, and local factor.
2. **Change the reader and program schedule jointly.** Start from the
   exact XOR word query, retaining arbitrary dirty output offsets. Test
   a block or hierarchical reader with an explicit native output
   schedule, paired with the surviving interpreter from the first lane.
   Charge private pools and shared controls at simultaneous width. A
   faster reader without a valid ordered interpreter is a component,
   not the terminal theorem; an interpreter that still serializes q
   output bits does not pass the depth gate.

These lanes cooperate on one mechanism. They are not invitations to
start unrelated lower-bound, practical-emitter, or manuscript projects.
The existing prefix-blind dirty-cache obstruction remains applicable to
its stated interface. Continued Clifford access to P is allowed and
counted; a changed representation must identify which old hypothesis
it avoids. Never treat $`Y\mathbin{\oplus}f(P)`$ as a clean f(P).

### First falsifiable identity

Before a large resource derivation, use a preserved prefix selecting at
least four programs and two noncommuting native gates. Then include a
four-mode real frame whose earlier rotation changes a logical bit used
by a later node selector. Require arbitrary initial programs and both
dirty arrays, actual inverse words, all logical columns, and literal
phases. A full operator identity implies the external-reference contract;
a few computational inputs do not.

Reconstruct the retained
[noncommuting-cache failure](JOINT_DEPTH_SESSION_RESULT.md#5-a-generic-dirty-cache-correction-already-contains-a-lookup)
and [changed-address residue](SHARED_PREFIX_QUERY_AUDIT.md#1-the-exact-residue-of-a-retained-prefix-short-echo)
when the proposed word resembles them. Use the one-bank commutator as
a negative control if the new construction merges bank roles. A
commuting example, a clean program, or a root-column match is insufficient.

Small exact or numerical checks may reject a specific algebraic proposal.
A passing fixture earns a symbolic proof, not an asymptotic theorem.
Immediately after that proof, substitute the benchmark into the q, S,
catalogue-volume, query-round, and simultaneous-width bounds. Stop a
candidate whose improvement disappears there.

## 4. Bounded baseline: quantify R90's existing precision window

Cap the opening baseline audit at 45 minutes. An independently useful
target, using the retained construction, is an absolute c and one
absolute dimension cutoff for

```math
6\le L\le c\bigl[\log_2(\ell+2)\bigr]^{1/3},
\qquad D_T=O(NL/b^2+\ell),
```

with absolute count constants, original widths, and the full real/complex
contract. This is an unbounded but very slowly growing precision range,
not the inverse-polynomial benchmark. It is a proof target, not a result
already supplied by the fixed-L statement.

Prove a uniform stopping height $`O((L+1)^3)`$ and, after at most three
coarse stages, a bound of the form

```math
m\le C\left[L+\log_2\log_2(\ell+4)+1\right]^3.
```

Choose c after the native constants, with strict margins so that q and
S fit below order ell. Check the full catalogue count and both arrays,
not merely local synthesis depth. The query volume is $`O(NL)`$.
The coarse local tolerances include $`\lceil\log_2(m_i+2)\rceil`$;
they are not all equal to L. Use one uniform dimension cutoff, verify
every minimum query pool, and price the complex tolerance shift. In
the narrow branch, use the whole original-accuracy fallback at the
original literal width. No L-dependent finite maximum may be hidden
in an allegedly absolute constant.

If these inequalities close, record the result for independent audit.
If they do not close within the opening cap, record the missing
inequality and move to the main experiment. Complete the baseline
later if the main mechanism fails its gate. A possible stronger output
is an explicit absolute precision cost, such as

```math
D_T=O\!\left(NL/b^2+
\min\{F_{79}(n,L),\ell+2^{C(L+1)^3}\}\right),
```

where $`F_{79}=\min\{nL,(L+\ell)\log_2(n+L+2)\}`$.
This too requires proof, including its large-dimension catalogue branch
and complementary R79 fallback. Do not turn it into a third research
lane or replace a sharper retained high-precision R79 bound with it.

## 5. Five-hour schedule and evidence required to finish

| Elapsed time | Work and decision |
|---|---|
| 0–45 minutes | Freeze the baseline, audit its absolute precision window, and enumerate every serial long-query boundary |
| 45–120 minutes | Run the two constructive lanes and skeptical review; require a full-space identity and a credible simultaneous resource ledger by the checkpoint |
| 120–210 minutes | Develop the surviving mechanism into the benchmark theorem; if neither lane passes, complete the bounded precision theorem and preserve the exact failed interface |
| 210–255 minutes | Independent proof audit of the same complete circuit: counts, depth, peak width, all-column error, signs, complex scope, and original-width splices |
| 255–300 minutes | Integrate only established claims; run required validation, rendered-document checks and CI, and merge when all gates pass |

Stop new construction by the audit window. Parallel workers must have
separate ownership before editing. A failed 120-minute gate means
changing or pausing that candidate; it does not retire the general
question. Preserve the five-hour ceiling and report an unfinished gate
honestly rather than weakening a check to meet the clock.

A full primary success is the improved inverse-polynomial-error
benchmark on one circuit with all simultaneous budgets. A complete
absolute growing-precision theorem is a narrower success. A new native
lemma is useful at its proved scope with its remaining composition gap
identified. If neither theorem closes, retain R79/R90 unchanged and
write the smallest explicit failed identity or unpaid cost, plus the
condition under which that route should be reopened.

Before promoting any result, require a register map at every live
boundary; a full-word identity or norm proof; separate T, Clifford,
depth, peak-width, and error ledgers; and an independently checked
composition. Broad finite validation follows a surviving proof and
resource argument. Preserve immutable scientific receipts and normal
repository gates. A new constant, a fixture collection, or an abstract
unitary completion is not a new general-depth theorem.
