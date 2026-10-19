# Two-band source obstruction and the continuing physical-pair scout

[Workspace](../../WORKSPACE.md) · [Resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous session](OVERLAPPING_SOURCE_SESSION_RESULT.md)

This continuation fetched main
`eb368d2fcccfac2c1dc22d66d0e8c8d979891384`, tree
`715b96aa427daa89ce1fffd15a024b4494959740`, on 19 October 2026.
It resolves the current two-band source question at the required
every-fixed-accuracy scope. It does not construct a faster indicator or
establish an unrestricted depth lower bound. R96 remains the strongest
unconditional general upper theorem, and R79 remains selected Result D.

## 1. The complete source interface is stopped

The target remains the whole indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle,
\qquad m=2^p,
```

on arbitrary address, output and dirty-work inputs, also with references.
At every fixed requested tolerance tau, the desired zero-clean word has
$`T,b=O_\tau(m)`$, $`G=O_\tau(m^2)`$, and depth
$`o(\log(p+2))`$. Actual inverses, phases, all return errors and leakage
remain charged.

The [fixed Clifford-sum theorem](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
proves: for any fixed clean count k and fixed positive amplitude alpha,
there are positive constants epsilon_0 and m_0 such that no actual
one-band encoder pair U_A,U_B satisfies

```math
\|J^\dagger U_B^\dagger U_AJ-\alpha I_p\|
\le\epsilon_0\qquad(m\ge m_0).
```

Here J inserts only the specified clean zeros. All other columns are
retained and the rejected blocks are unrestricted. The result allows
arbitrary width and native counts and does not need exact flatness.
Uniformly bounded clean counts are covered by padding to k_0.

The proof uses a uniform normalized Hilbert–Schmidt gap for nonzero
fixed-length, fixed-coefficient sums of literal native Cliffords.
A projected output-probe defect has exactly
$`2^{3k}+2^k`$ such terms. The indicator's sparse address action drives
its normalized Hilbert–Schmidt norm below the gap, so it must vanish exactly.
That contradicts the same probe's full operator action. No Clifford
sum, trace, or projector is supplied as a free native operation.

The proof includes a constructive but deliberately conservative gap.
For one clean signal and amplitude $`1/\sqrt2`$, it permits
$`g=2^{-1280499912}`$, $`\epsilon_0=g/4`$, and
$`p_0=2560999827`$. These constants illustrate existence; they are
not practical thresholds and do not rule out every modest numerical
tolerance. They do rule out a source meeting every fixed requested
accuracy for all address sizes. Unbounded clean count or changing
amplitude requires a new quantitative argument. The theorem does not
exclude a different constant-depth mechanism or close unrestricted
dimension-growing depth.

## 2. Complete transitions, repairs, and finite-size screens

The continuation first computes a paid Clifford repair of the earlier
growing signal-mixing family. Its old wrong-square test could not rule
out dressing: an explicit native Clifford makes the transition a
Hermitian involution. Its forward and actual-inverse counts, source
phase, recovered signal return and all logical columns are retained.

The [quarter-turn coefficient screen](OVERLAPPING_TWO_BAND_SOURCES.md#6-a-clifford-invariant-coefficient-screen-for-quarter-turn-transitions)
then excludes arbitrary independent logical Clifford pre/post dressing
of that repaired family. More generally, a quarter turn about a Hermitian
Clifford involution has either a Pauli output probe or every Pauli
coefficient at most $`1/\sqrt2`$. For p at least three, comparison with
the indicator's coefficient $`1-2/m`$ gives the source floor

```math
\frac{1-2/m-1/\sqrt2}{2\sqrt2}.
```

The specific growing family has the stronger coefficient cap one half
and floor at least $`1/(8\sqrt2)`$. These finite-p results retain value
despite the subsequent general theorem's much smaller explicit constant.

The [simultaneous parity analysis](SIMULTANEOUS_OUTPUT_PARITY_LIMITS.md)
shows that sufficiently accurate sources need output-mask support rank
at least $`m-\sum_{j=0}^3\binom pj`$. This is a linear native support
condition, compatible with the desired count budget. A paid two-band
example saturates both directional parity conditions while its complete
transition is identity. Hence rank saturation and nonscalar probes alone
do not identify the required indicator action. The fixed-sum theorem
supersedes this as the general two-band resumption question.

## 3. First unpaid operation and complete continuation recurrence

The next concrete interface is a one-signal source U with at least three
physical T bands and complete accepted block near $`I_p/\sqrt2`$.
The first unpaid operation is its **complete logical transition** with
uniform linear T-count/dirty width and quadratic Clifford count, including
the actual inverse. A formal unitary dilation or flat success probability
does not supply it. With a third source band the projected defect is no
longer automatically a fixed-length Clifford sum; this removes a proof
hypothesis but does not supply a source.

The [authoritative recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#62-the-next-source-band-retains-the-general-leakage-recurrence)
retains twenty-four source/actual-inverse calls. For a three-band source it gives

```math
\begin{aligned}
T_{\mathrm{ind}}&\le24T_U+48m+8,\\
G_{\mathrm{ind}}&\le24G_U+82m+104,\\
b_{\mathrm{ind}}&\le m+b_U+2,\\
D_T(\mathrm{ind})&\le24d_U+32\le104,\\
\epsilon_{\mathrm{ind}}&\le16(3\sqrt2\epsilon+2\epsilon^2).
\end{aligned}
```

Even m removes the last eight T gates. Requesting
$`\epsilon\le\tau/72`$ suffices for $`0\lt\tau\le1`$.
Do not carry over the two-band exact-flatness rounding. Signal and
toggle return are approximate and already included in this complete
operator error. Independently proved exact flatness would restore the
sharper $`2\sqrt2\epsilon`$ bound and exact signal/toggle return.

For R98, request $`\epsilon_L=2^{-L}/(1152K_*)`$. Four completed
indicators per query and at most K_* long-query calls then add at most
one quarter of the final error budget. The retained real and complex
totals are $`149\cdot2^{-L}/512`$ and $`789\cdot2^{-L}/1024`$.
The resulting conditional frame depth would be
$`O_L(N/b^2+\log_2\log_2(n+4))`$, with the retained one-clean
allocation, dirty thresholds and R96 fallback branches. It is conditional
on the unpaid source transition and does not change the current frontier.

## 4. The original physical-pair scout remains active

The earlier fixed-(s,t) degree theorem and clean-preserving Clifford-sum
reduction remain in the [occupied-center home](OCCUPIED_BOOLEAN_PAIR_CENTER.md).
The new [central-module profile theorem](CENTRAL_MODULE_NATIVE_PROFILES.md)
uses the actual original physical-pair reduction and genuine common
center R, then restricts only already clean-preserving factors.

It proves three corrected native profile distortions, bounded by
$`8c,4c,4c`$ for $`c\le3s+1`$, without bounding the coupled silent
span t. For each element of the actual Boolean matching family, the
escaped middle-axis span modulo R is at most $`8c`$. The profiles retain
central relation forms and native orientation signs; raw active rank
cannot replace them.

This does not bound the union of escaped axes, transformations within
the original middle span plus R, or t from s alone. The next precise
question is simultaneous native matching using the actual shared
square-zero displacement family and all three profiles. Fixed positive
s alone, growing pair count, and the unrestricted one-clean three-layer
degree question remain open. The closed no-middle-pair branch is not
reopened.

## 5. Verification boundary

Independent analytic audits checked the fixed-sum gap, native phases,
projected-defect normalization, explicit constants, source/return
recurrence, and central-module norms. Small exact fixtures cover the
repaired transition and actual inverse, quarter-turn coefficient cases,
both parity saturation examples, nonflat one-/two-clean projected defects,
literal ten-/sixty-eight-term expansions, sparse indicator probes, and
the corrected module profiles. These checks support the identities;
they do not substitute for their uniform proofs or emit a faster source.

Repository validation and complete rendered-document gates are recorded
on the final reviewed pull-request tree. No receipt hash, selected
manuscript theorem, attribution, or tracked upstream revision is changed.
This record and WORKSPACE supersede earlier prompts to continue general
two-band source discovery. The unrestricted depth destination and focused
original physical-pair scout remain active.
