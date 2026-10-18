# Six-layer Pauli comparisons and diagonal-source bounds

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Six-layer proof](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Diagonal-source proof](RETURNED_DIAGONAL_ECHO_DEGREE.md) ·
[Previous checkpoint](NATIVE_DEFECT_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`67d886b1e79b6c63c7531860edfbaf8378d7dcb4` (PR #173).
It tested native six-layer defect transport and a constructive diagonal-source
alternative. Two scoped bounds and an actual-return counterexample are
supported. Optimal general T-depth remains the destination; final manuscript
drafting is separate.

## 1. Compact revision of the scientific state

The count package is settled under its stated hypotheses. R96 remains the
strongest unconditional complete-frame upper law for every
$`n\ge1,L\ge6`$, with full error at most $`2^{-L}`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad N=2^n.
```

These are simultaneous bounds on one circuit, with one clean qubit and
real dirty width $`b\ge17(L+n+7)`$ or complex magnitude width
$`b\ge17(L+n+8)`$. Every logical column, literal sign and phase,
singular and zero coefficient, arbitrary dirty/reference input, actual
inverse, leakage and complete workspace return remains required.
R79 remains selected Result D. R97 retains its outer-overlap hypotheses;
unrestricted one-clean four-layer rigidity remains open.

The exact general question is the optimal worst-case T-depth for the
prescribed complete frame at the same precision, allocation and simultaneous
T-count and elementary Clifford budgets. A faster exact multi-output reader
satisfying R98's native interface and resource envelope would improve the
fixed-accuracy upper bound. A growing unrestricted native-query lower bound
would transfer through the two-call corner reduction. Neither missing
premise is supplied in this round.

Formerly open questions already resolved include optimal one-clean count,
R96's all-precision upper law, zero-clean read-only three-layer affine and
four-layer quadratic XOR rigidity, and the general five-layer degree bound
14 below $`\sqrt2/24`$ (11 for invisible outer target shear). The latter
requires no extra comparison hypothesis. This round additionally resolves
a width-independent six-layer bound under the specified Pauli-comparison
hypothesis, and the fixed diagonal-source/affine-transport degree question.
Sharp degrees and optimal general depth remain open.

## 2. A six-layer subcase without intermediate Clifford assumptions

The [six-layer proof](READ_ONLY_SIX_LAYER_REDUCTION.md) establishes
$`\deg f\le51`$ below full operator error $`\sqrt2/64`$, independently
of dirty width, when its specified target comparison K is Pauli-valued.
The model has zero initialized work qubits, original address preserved by
every completed Clifford interlayer, and physical T bands acting on work.
T-count, Clifford count and finite dirty width are unrestricted in this
necessary condition; the analysis adds no physical operation or allocation.

The preceding five-layer mechanism starts with a Pauli defect. At six layers
the native inner family is instead a product of two one-band Pauli images
and a Pauli. Its Clifford-valued defect is conjugated through three bands.
Repeated one-band matching cannot assume the intervening operators Clifford.
The new subcase avoids that missing step: K has a quadratic Pauli-label
baseline plus a fixed target label times the output bit. An eight-factor
cube cancels the baseline. Its exact native Clifford pullback has symplectic
degree 16 and literal fixed-probe sign degree 48. A single frozen reference,
with the full $`16\eta/32\eta`$ error accounting, detects every third
difference. A non-Pauli reference improves the branch bound to 19.
The proof does not assert that the transported target itself is Clifford
under positive error, or that either conservative bound is attained.

## 3. Full return does not remove odd-order defects

[Section 8](READ_ONLY_SIX_LAYER_REDUCTION.md#8-an-actual-returned-comparison-with-an-odd-order-defect)
embeds the formerly inner-only diagnostic into an exact fully returned
identity query. Its actual specified comparison has projective order three;
its ordered defect has literal order three. Repeated multiplicative
differences remain non-Pauli indefinitely. The proposed universal 2-group
or repeated-Pauli-return recurrence therefore fails even with full return.

The displayed construction has four T gates in four nonempty slots of the
six-slot decomposition, fifteen elementary Cliffords, zero clean or extra
helper qubits, one arbitrary work/output qubit and two read-only address
bits. Its actual adjoint has four T gates and thirteen elementary Cliffords.
A literal address phase is essential. Both complete words simplify to
identity, whose fair baseline has zero gates and depth on the same logical
registers. Merging empty slots changes the defined comparison; this does not
contradict the five-layer theorem or improve a compiler.

## 4. A broader diagonal-source obstruction and larger echo radius

The [diagonal-source theorem](RETURNED_DIAGONAL_ECHO_DEGREE.md) permits q
occurrences of arbitrary diagonal work sources, including irrational
phases and actual inverses, joined by read-only monomial Clifford intervals.
All work is arbitrary. Invertible affine work mixing and literal Clifford
phases, including address/work cross phases, are included. Put
$`r=\max(q,3)`$. Full Boolean scalar return at error

```math
\varepsilon\lt2\sin\!\left(\frac{\pi}{2^{r+1}}\right)
```

forces degree at most $`r-1=\max(2,q-1)`$. Joint address/work difference
directions annihilate each source separately, while three differences kill
the intervening Clifford phase. A circular error argument needs no real lift
of the source phases. The four-call echo remains sharply cubic, now below
$`2\sin(\pi/32)`$ instead of the historical $`2/275`$ radius.

Every inverse is a separate source occurrence. Native source errors add,
forward and inverse elementary Clifford counts are distinguished, and every
simultaneously used helper is charged. The radius decreases with q; occurrence
count is not native T-depth. Nonmonomial interference, nonlinear work
permutations and unrestricted fused synthesis remain outside the theorem.
Thus diagonal-source fusion has not removed R91's logarithmic dirty-indicator
cost. No faster fully charged reader was found.

## 5. Precise unresolved step and resumption question

The most useful next question is bounded native transport of the general
six-layer ordered defect, allowing its odd-order modes. Before the common
three-band prefix, its symplectic degree is at most four and its literal
sign degree on degree-two moving probes is at most eighteen. Its conjugate
is within $`6\eta`$ of the actual Clifford comparison defect.

[Section 6](READ_ONLY_SIX_LAYER_REDUCTION.md#6-a-conditional-defect-detector-for-the-general-comparison)
proves the sufficient algebraic interface: if the actual defect has
symplectic degree at most d and moving-probe sign degree at most s, independent
of the unknown output degree, then the output degree is at most $`2d+4`$
for visible outer target shear and $`\max(s,6)+1`$ for invisible shear.
The missing step is a width-independent native bound on d and s, with a
controlled error radius, through that specified three-band prefix.
Intermediate Clifford membership cannot be assumed. The actual odd-order
fixture blocks nilpotence, but does not block this signed-coordinate route.

**Resumption question:** can the full-return compatibility and the two
native one-band factors bound this three-band transport, or can an actual
returned family make its degree amplification grow with dirty width?
Success on either side would resolve a genuinely open six-layer invariant
and guide a depth recurrence. The recurrence must track error as depth grows;
one-clean work and unrestricted address mixing remain separate requirements
before a dimension-growing general-frame lower bound follows. Constructively,
a nonmonomial or nonlinear returned mechanism is worth reopening only with
its loading, controls, recovery and actual inverse fully charged against R91.

## 6. Evidence and integration

The six-layer proof has one authoritative home. It retains the
[three native-boundary checks](../../tests/test_six_layer_defect_boundary.py)
and adds [six Pauli-comparison checks](../../tests/test_six_layer_pauli_comparison.py)
and [four actual-return checks](../../tests/test_returned_six_layer_odd_group.py).
The source theorem extends its existing home with
[four source-word checks](../../tests/test_diagonal_source_word.py), retaining
the [eight echo checks](../../tests/test_returned_diagonal_echo_degree.py).
Finite checks support the displayed identities and ledgers; the analytic
arguments establish the arbitrary-width bounds. The checkpoint and frontdoors
state their hypotheses without promoting them to an unrestricted depth law.
