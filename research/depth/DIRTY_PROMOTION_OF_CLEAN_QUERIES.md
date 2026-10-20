# Removing a bounded clean bank from a complete XOR query

[Current native-sharing revision](NATIVE_SHARING_SESSION_RESULT.md) ·
[Complete approximate frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Table-specific interfaces](TABLE_SPECIFIC_INDICATOR_INTERFACES.md) ·
[Native sharing requirement](NATIVE_INDICATOR_SHARING.md)

A complete XOR query with a bounded number of clean source qubits can be
promoted to a zero-clean query at constant overhead. The source may have
arbitrary nonmonomial behavior outside its clean sector. The construction
uses its actual inverse, a paid controlled bilinear phase, a returned
zero-test, and a cover of all clean-register sectors. It does not control
the source circuit or assume identity on its unpromised inputs.

This changes the sufficient source question: a uniformly bounded clean
bank is permitted if the whole approximate indicator, including every
arbitrary output and dirty-work input, meets the desired resource bounds.
It supplies no such faster source. The emitted cover has exponential
depth in the number of source-clean bits; a logarithmic clean bank does
not immediately give the requested sublogarithmic depth. R96 remains the
strongest unconditional general upper theorem, and the complete-frame
frontier is unchanged.

## 1. The source is a complete query, not a prepared state

Let x be the address, T a q-bit output bank with $`q\ge1`$, A a
k-bit clean source bank with $`k\ge1`$, and W arbitrary private
source work. For any Boolean vector $`f(x)\in\mathbb F_2^q`$,
write

```math
Q_f|x,T,W\rangle=|x,T+f(x),W\rangle.
```

Let J insert $`A=0^k`$ and let $`P=JJ^\dagger`$. A native
unitary V has the exact complete clean-code contract when

```math
VJ=JQ_f.                                                     \qquad(1)
```

Every x,T,W input, every superposition, and every external reference is
included. Since Q_f is unitary on that whole sector, (1) implies
$`[V,P]=0`$. A promise only for $`T=0`$, or a state-preparation
promise, would not imply this invariant-sector statement. Off the
sector P, V can have any unitary action, including address/output mixing
and noncommuting basis changes.

The approximate version is the initialized-isometry bound

```math
\|VJ-JQ_f\|_{\mathrm{op}}\le\epsilon\lt1.                    \qquad(2)
```

This includes all address, output, source-clean return, dirty-work return,
and reference errors. Exact preservation of P is not assumed in (2).
Common bounds $`T_V,G_V,d_V,b_V`$ charge the emitted V and its
actual native inverse; b_V excludes the specified x,T,A registers.
For zero clean bits no promotion is needed.

The same statements allow a fixed source phase:
$`VJ=e^{i\theta}JQ_f`$, or the corresponding error in (2).
Replace V by $`e^{-i\theta}V`$ only in the mathematical
comparison. Every emitted occurrence is inside a conjugation
$`V^\dagger K V`$, so the phase cancels exactly. No native
global-phase correction is supplied for free, and input-dependent
source phases are not admitted by this observation.

## 2. A selective native echo removes the unknown off-code action

Borrow A, a q-bit old output T, one toggle bit d, and the source's private
work W in arbitrary states. Let Z be the specified arbitrary q-bit output
of the new query. V acts only on x,T,A,W. Define Hermitian involutions

```math
B=\prod_{j=1}^q\mathrm{CZ}(T_j,Z_j),\qquad
K=\prod_{j=1}^q\mathrm{CCZ}(d,T_j,Z_j),\qquad
D=V^\dagger K V.                                            \qquad(3)
```

K is the d-controlled B. The word for D is chronologically V, K,
then the actual inverse of V. D is a Hermitian involution for every
unitary V, including an approximate source that leaks out of P.

Let C be an exact zero-clean returned query

```math
C:|A,d,w_C\rangle\longmapsto
       |A,d+[A=0^k],w_C\rangle.                              \qquad(4)
```

Its private dirty pool is disjoint from W; it has no action on x,T,Z.
For any unitary R define the emitted echo by

```math
\mathcal E(C,R)=R^\dagger C^\dagger R C.
```

The chronological order is C, R, the actual inverse of C, then the
actual inverse of R. When R is the d-controlled F, with
$`F^2=I`$ and $`[F,P]=0`$, C changes the d control precisely on
P. The two controlled F occurrences cancel on the complement and leave
one F on P, independently of the arbitrary incoming d. Thus

```math
\mathcal E(C,R)=P F+(I-P)I.                                 \qquad(5)
```

Under (1), $`F=V^\dagger B V`$ is an involution preserving P,
and its restriction there is $`Q_f^\dagger BQ_f`$. Since the
source is a complete XOR query,

```math
Q_f^\dagger BQ_f=(-1)^{f(x)\cdot Z}B.                       \qquad(6)
```

Emit the K echo first and the D echo second. In operator order the
result is

```math
S_0=\mathcal E(C,D)\mathcal E(C,K)
   =P(-1)^{f(x)\cdot Z}+(I-P)I.                             \qquad(7)
```

Every borrowed register returns. This is an exact full-operator
identity, including the literal scalar phase, with no promise on the
off-code V action. One selective word has four V/actual-inverse
appearances, four C/actual-inverse appearances, and four K/actual-inverse
appearances. No controlled-V is used. Controlling the Clifford B is
paid by K's native CCZ word; it is not treated as a Clifford operation.

## 3. Covering all sectors and paying the complete recurrence

For $`b\in\mathbb F_2^k`$, let X_b add b to A. The conjugate
$`S_b=X_bS_0X_b`$ implements the phase in (7) precisely on
$`A=b`$. The product of all these complete words is

```math
\prod_{b\in\mathbb F_2^k} S_b
      =(-1)^{f(x)\cdot Z}\otimes I_{A,T,d,W,w_C}.             \qquad(8)
```

Hadamards on Z before and after (8) give the desired complete XOR
query into its arbitrary output Z, with zero initialized work. The old
source output T is borrowed and returned. The Hadamards apply to the
whole Z bank, rather than being independently repeated inside a row
compiler.

Put $`M=2^k`$. A cyclic binary Gray order starts at zero, visits all
M masks, and returns to zero using M one-bit changes. Adjacent X_b
words therefore reduce to M individual X gates. Each X is the literal
native word HSSH, with four Cliffords. This saves the otherwise
unnecessary factor k in the translation count. It does not parallelize
the M selective words.

The controlled bilinear phase K has the native bounds

```math
T_K\le6q+(q\bmod2)\le6q+1,\qquad
G_K\le10q+3,\qquad D_T(K)\le4.                              \qquad(9)
```

For completeness, the ordinary exact CCZ parity polynomial has three
positive singleton terms, three negative pair terms, and one positive
triple term. Across the q triples, the common singleton phase on d
combines to $`T_d^q`$, requiring at most one T or T-dagger and
three positive-S Cliffords, also for its actual inverse. The remaining
six terms per triple give six T gates. Compute and uncompute the
two d/pair parities, the T/Z parity, and the triple parity with ten
CNOTs per pair. Their phases fit four physical T bands: the singleton
band, the two d/pair terms on their distinct T and Z wires, the T/Z
term, and the triple term. Shared d occurs only as a Clifford control
between these bands. No clean or helper register is needed. This is the
same shared-control scheduling principle used by the
[native batch compiler](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers).

The [returned dirty conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
supplies C, including its negative literals. Safe common majorants are

```math
T_C,G_C,b_C\le c(k+3)^3,\qquad d_C\le c\log_2(k+3).          \qquad(10)
```

All its helpers return on arbitrary inputs. For $`k=1`$, C is just
X_d followed by $`\mathrm{CX}(A,d)`$, so
$`T_C=d_C=b_C=0`$ and $`G_C\le5`$. The general recurrence is

```math
\begin{aligned}
T_{\mathrm{prom}}&\le M(4T_V+4T_C+24q+4),\\
G_{\mathrm{prom}}&\le M(4G_V+4G_C+40q+16)+2q,\\
D_T(\mathrm{prom})&\le M(4d_V+4d_C+16),\\
b_{\mathrm{prom}}&\le k+q+1+b_V+b_C.
\end{aligned}                                               \qquad(11)
```

These are common forward/actual-inverse bounds in the positive-S native
alphabet. The specified Z output and x address are excluded from the
helper count; the old q-bit T bank, k-bit A bank, d, and both private
pools are included. The same pools are reused sequentially, so width
does not receive the factor M. The helper count does not presume that
any reused register is clean.

With one clean source bit the concrete majorants become

```math
\begin{aligned}
T_{\mathrm{prom}}&\le8T_V+48q+8,\\
G_{\mathrm{prom}}&\le8G_V+82q+72,\\
D_T(\mathrm{prom})&\le8d_V+32,\qquad
b_{\mathrm{prom}}\le q+b_V+2.
\end{aligned}                                               \qquad(12)
```

If q is even, the final eight in the T bound may be omitted. This
includes indicators with $`q=2^p`$ and $`p\ge1`$.

## 4. Approximate clean return still gives a full-operator query

Under (2), a dimension-independent block-extension argument controls
the entirely arbitrary inputs of the promoted word. In the P and
complement decomposition write

```math
V=\begin{pmatrix} A_0&B_0\\ C_0&D_0\end{pmatrix}.
```

Equation (2) gives
$`\|A_0-Q_f\|\le\epsilon`$ and $`\|C_0\|\le\epsilon`$.
Unitarity and the square block A_0 imply
$`\|B_0\|=\|C_0\|`$: the spectra of $`A_0A_0^\dagger`$
and $`A_0^\dagger A_0`$ coincide, and the complementary defects
are $`B_0B_0^\dagger`$ and $`C_0^\dagger C_0`$.
Moreover,

```math
D_0^\dagger D_0=I-B_0^\dagger B_0\succeq(1-\epsilon^2)I.
```

The polar unitary W_0 of D_0 exists and satisfies

```math
\|D_0-W_0\|\le1-\sqrt{1-\epsilon^2}\le\epsilon.
```

Define the mathematical comparison extension
$`V_0=Q_f\oplus W_0`$. Its diagonal-block difference from V has
norm at most epsilon, and its off-diagonal-block difference has norm
$`\max\{\|B_0\|,\|C_0\|\}\le\epsilon`$. Consequently

```math
\|V-V_0\|_{\mathrm{op}}\le2\epsilon.                       \qquad(13)
```

V_0 is not an unpaid compiler operation: it is used only in the proof.
It meets (1), so replacing all source occurrences by V_0 in the word
of Sections 2–3 gives the exact target. There are four actual source
or inverse occurrences per selective sector. Unitary telescoping,
including the actual adjoints, therefore gives

```math
\|\widetilde Q_f-Q_f\otimes I_{\mathrm{borrowed}}\|_{\mathrm{op}}
\le8\cdot2^k\epsilon.                                      \qquad(14)
```

No independent-error or fresh-work assumption enters this estimate.
Every later actual word has norm one even on previously leaked
components. The distance is on the entire physical space, so includes
all output, address, borrowed-return and external-reference errors.
For a prescribed k, the sufficient source tolerance is
$`\epsilon=\tau/(8\cdot2^k)`$.

There is a stronger estimate under the distinct hypothesis that the
actual emitted approximate V preserves its clean sector exactly:

```math
[V,P]=0,\qquad R=J^\dagger VJ,\qquad \|R-Q_f\|\le\epsilon.
                                                               \qquad(15)
```

The two selective echoes still satisfy (5). S_0 is exactly identity
off P, and on P it is $`R^\dagger B R B`$, in that order. It
need not be a scalar phase. The X_b conjugates have disjoint exact
A-sector supports, with precisely the same restricted word on each.
Their product is therefore
$`I_A\otimes R^\dagger B R B`$, rather than an accumulation of
M error terms on any input sector. Two-source telescoping yields

```math
\|\widetilde Q_f-Q_f\otimes I_{\mathrm{borrowed}}\|_{\mathrm{op}}
\le2\epsilon.                                               \qquad(16)
```

A, d and C's private pool now return exactly. All possible residual
errors on x,T,Z,W are included in (16). The direct-sum estimate cannot
be inferred merely from small approximate clean-return error in (2).

## 5. A paid sufficient source for R98

Suppose a terminating native compiler supplies complete approximate
indicators with $`q=m=2^p`$, a clean count bounded by a fixed k_0
uniformly over p and all requested source tolerances, and, at each fixed
$`\epsilon>0`$, common forward/actual-inverse bounds

```math
\begin{aligned}
T_V,b_V&=O_\epsilon(m),\qquad G_V=O_\epsilon(m^2),\\
d_V&=O_\epsilon(d_\epsilon(p)),\qquad d_\epsilon(p)\ge1.
\end{aligned}                                               \qquad(17)
```

Padding the ignored clean bank to k_0, request
$`\epsilon=\tau/(8\cdot2^{k_0})`$. Equations (11) and (14)
then give a complete zero-clean indicator of tolerance tau with
$`T,b=O_\tau(m)`$, $`G=O_\tau(m^2)`$, and depth
$`O_\tau(d_\epsilon(p))`$. In particular,
$`d_\epsilon(p)=o(\log(p+2))`$ preserves the desired
sublogarithmic depth. The new q-bit intermediate output is absorbed
because $`q=m`$.

For a general query with many more output bits than rows, the explicit
q terms in (11) need not be absorbed by its original count or width
bound. No general R91 word-query frontier is asserted directly from
(11). The frame consequence here goes through the promoted complete
indicator and then the retained R91 allocation.

That allocation uses at most four completed indicator occurrences per
word query, including all actual inverses. Hence (14) adds at most
$`32\cdot2^{k_0}\epsilon`$ per complete query. With the fixed
R98 long-query count K_* and $`\eta=2^{-L}`$, set

```math
\epsilon_L=\frac{\eta}{128\cdot2^{k_0}K_*}.                 \qquad(18)
```

All added query error is at most $`\eta/4`$. The retained
[constant-accuracy frame audit](APPROXIMATE_QUERY_FRAME_TRANSFER.md#3-the-bounded-r98-call-count-and-final-error)
therefore gives total errors at most $`149\eta/512`$ for real
frames and $`789\eta/1024`$ for complex frames. The bilinear
branch pays both complete indicator pairs with its actual restored
addresses; $`H^2+J^2=O(Qq)`$ and the direct branch
$`Q^2\le Qq`$ absorb the added quadratic Clifford terms as in
that audit. All source constants and dirty-pool constants must be
uniform in address length and table; there is no new supplied table
operation.

For nondecreasing depth envelopes, the conditional result is

```math
\begin{aligned}
T&=O_L(\sqrt N+N/b),\qquad G=O_L(N),\\
D_T&=O_L\!\left(N/b^2+d_{\epsilon_L}(2n+10)
                         +\log_2\log_2(n+4)\right).
\end{aligned}                                               \qquad(19)
```

It retains one clean frame flag and the original dirty thresholds
$`b\ge17(L+n+7)`$ for real frames and
$`b\ge17(L+n+8)`$ for complex frames, using R96 for the small
or narrow fallback branches. The source's clean allocation has actually
been removed, rather than being charged to that already occupied flag.
No new unconditional general upper bound follows without the faster
source in (17).

If the clean count k_epsilon is only finite at each requested epsilon
but grows arbitrarily as epsilon decreases, (14) alone does not prove
that every target tolerance is attainable: one must solve
$`8\cdot2^{k_\epsilon}\epsilon\le\tau`$. Uniformly bounded
k_0 is a sufficient hypothesis. Under the stronger exact-invariance
condition (15), this circular accuracy issue disappears. Set
$`\epsilon=\tau/2`$ in (16); a finite k_epsilon independent
of p then gives a tolerance-dependent constant resource overhead.
For the corresponding frame statement use
$`\epsilon_L=\eta/(32K_*)`$. That is a separate sufficient
source interface, not a consequence of approximate return alone.

## 6. A one-signal block is another sufficient source interface

The existing
[calibrated amplification](INTERVAL_DILATION_DEPTH.md#6-uniform-calibrated-amplification)
and [normalization-two recovery](COARSE_FLAG_SPARSE_CORRECTION.md#6-native-sparse-residual-sources-and-normalization-two-recovery)
can be composed with this clean-removal word. The amplification method
itself is retained; the additional conclusion is a completely paid
zero-clean interface from a single signal qubit.

Let U be native on x,T,W and one signal a, let J insert
$`a=0`$, and require its entire compressed block to satisfy

```math
A_U=J^\dagger UJ,\qquad
\|A_U-Q_f/2\|_{\mathrm{op}}\le\varepsilon\le1/8.             \qquad(20)
```

Every logical address/output/dirty-work input is included. No promise
on the other blocks or on returned signal amplitude is made. Set
$`R=I-2JJ^\dagger=-Z_a`$ and
$`F=-URU^\dagger RU`$. Direct expansion gives the exact
full-isometry identity

```math
FJ=2JA_U+UJ(I-4A_U^\dagger A_U).                            \qquad(21)
```

Writing $`A_U=Q_f/2+E`$ with $`\|E\|\le\varepsilon`$
shows

```math
\|I-4A_U^\dagger A_U\|\le4\varepsilon+4\varepsilon^2,
\qquad
\|FJ-JQ_f\|\le6\varepsilon+4\varepsilon^2.                 \qquad(22)
```

In particular, (22) includes the rejected amplitude and all work-return
errors; a square-root leakage estimate is unnecessary. Emit the actual
word $`V=UZ_aU^\dagger Z_aU=-F`$. Its clean-code target has the
fixed phase minus one, which cancels in Section 2's conjugations. Two
literal Z words cost four positive-S Cliffords, so common bounds are
$`T_V\le3T_U`$, $`G_V\le3G_U+4`$,
$`d_V\le3d_U`$ and $`b_V=b_U`$. The inverse uses the actual
U adjoints and the same literal Z words.

The one-clean promotion consequently gives

```math
\begin{aligned}
T_{\mathrm{prom}}&\le24T_U+48q+8,\\
G_{\mathrm{prom}}&\le24G_U+82q+104,\\
D_T(\mathrm{prom})&\le24d_U+32,\qquad
b_{\mathrm{prom}}\le q+b_U+2,\\
\|\widetilde Q_f-Q_f\otimes I\|_{\mathrm{op}}
&\le16(6\varepsilon+4\varepsilon^2).
\end{aligned}                                               \qquad(23)
```

There are twenty-four U/actual-inverse occurrences. For even q the
final eight T gates may be omitted. For a fixed requested tolerance
$`0\lt\tau\le1`$, $`\varepsilon=\tau/112`$ is sufficient.
With $`q=m`$, a native U of linear T-count and dirty width,
quadratic Clifford count and sublogarithmic T-depth would therefore
supply the original indicator target and Section 5's frame transfer.

Half amplitude is only one sufficient normalization. A source with
$`\|A_U-Q_f/\sqrt2\|\le\varepsilon\le1/8`$ has an equally
paid three-call recovery using phase-S reflections. Put
$`\delta=-1-i`$ and
$`R_S=I+\delta JJ^\dagger=-iS_a`$. Direct expansion now gives

```math
(UR_SU^\dagger R_SU)J
 =\delta JA_U+UJ(-iI+\delta^2 A_U^\dagger A_U).
```

For the ideal block, the coefficient of UJ vanishes. The raw native
word $`V_S=US_aU^\dagger S_aU`$ has target
$`e^{i\pi/4}Q_f`$ on its clean code, and the complete clean-isometry
error obeys

```math
\|V_SJ-e^{i\pi/4}JQ_f\|
 \le3\sqrt2\,\varepsilon+2\varepsilon^2.
```

Indeed, $`|\delta|=\sqrt2`$, $`|\delta^2|=2`$, and
$`\|A_U^\dagger A_U-I/2\|\le\sqrt2\,\varepsilon+\varepsilon^2`$.
The fixed phase cancels in the promotion, just as in the half-amplitude
case. The forward V_S has at most $`3G_U+2`$ Cliffords; its
actual inverse has at most $`3G_U+6`$, since both S inverses are
three positive-S gates. The one-clean promotion has four forward and
four inverse source occurrences, so these additional Cliffords sum to
thirty-two. All four resource bounds in (23) therefore remain valid,
while its error improves to
$`16(3\sqrt2\,\varepsilon+2\varepsilon^2)`$.
For $`0\lt\tau\le1`$, requesting $`\varepsilon=\tau/72`$
is sufficient. No negative-S gate is counted as one elementary
Clifford in this ledger.

These interfaces ask for complete constant-amplitude blocks, which may
be easier to discover than complete returned indicators. Excluding one
amplitude within a source grammar does not exclude the other or a
different paid recovery. Constructing a suitable native block is still
unpaid; a formal contraction or unspecified unitary dilation supplies
neither U nor its inverse ledger.

### 6.1. Exact flatness sharpens the complete error and return recurrence

For an actual two-band source at inverse-square-root-two amplitude,
[Clifford rounding](JOINT_FOURIER_WORD_LIMITS.md#5-a-precise-next-interface-for-two-source-bands)
makes the accepted block exactly flat whenever its error epsilon against
the complete target is at most one tenth. Write

```math
K=J^\dagger UJ,\qquad K^\dagger K=I/2,\qquad
L=\sqrt2K,\qquad \|L-Q_f\|\le\sqrt2\varepsilon.
```

L is unitary on the **whole** nonsignal space. Substitution in the exact
phase-S identity above gives, with $`\omega=e^{i\pi/4}`$,

```math
V_SJ=(1+i)JK=\omega JL,\qquad
V_S^\dagger J=\omega^{-1}JL^\dagger,\qquad [V_S,JJ^\dagger]=0.
```

The returned signal is exact, even though L may only approximate Q_f.
Apply Section 4's exact-invariance estimate, not its leakage estimate.
The promoted full operator is literally

```math
\widetilde Q_f
 =I_{a,d}\otimes\bigl(H_Z L^\dagger B L B H_Z\bigr),
\qquad B=\prod_{j=1}^q\mathrm{CZ}(T_j,Z_j),
\qquad H_Z=H^{\otimes q}.
```

Here L acts on the address, old arbitrary output T and source work W;
identity factors on unused registers are understood. The old signal a
and toggle d are now arbitrary borrowed inputs and return **exactly**.
Errors on the address, old output, source work and actual output Z are
all included in the following full-operator bound, also with references:

```math
\begin{aligned}
\|\widetilde Q_f-Q_f\otimes I_{\rm borrowed}\|
 &\le2\|L-Q_f\|\le2\sqrt2\varepsilon,\\
T_{\rm prom}&\le24T_U+48q+8,\\
G_{\rm prom}&\le24G_U+82q+104,\\
D_T(\mathrm{prom})&\le24d_U+32,\qquad
b_{\rm prom}\le q+b_U+2.                         \qquad(23a)
\end{aligned}
```

These are common forward/actual-inverse bounds; q even removes the
last eight T gates. There are four forward and four inverse recoveries,
each using three U/actual-inverse occurrences. Their positive-S overhead
is $`4\cdot2+4\cdot6=32`$ Cliffords. No multiplication of the error
by the number of clean sectors is needed because those sectors are
exactly invariant. The actual adjoint of the entire promoted word obeys
the same full error bound against the inverse query. No independent
return promise or reset is used at any later call.

For $`q=m=2^p`$ and $`0\lt\tau\le1`$, it is enough to request
$`\varepsilon=\min\{1/10,\tau/(2\sqrt2)\}`$. A genuine two-band
source would give depth at most 80 with the displayed count, width and
Clifford budgets, provided its complete logical transition and uniform
source budgets exist. This is a conditional recurrence, not a compiler
for that transition. The subsequent
[fixed Clifford-sum theorem](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) excludes
this two-band interface at sufficiently small fixed error and large p.
Its arbitrary-accuracy indicator premise is therefore unavailable.

The retained restored-address R98 allocation uses at most four completed
indicators per query and at most K_* long-query calls. At final precision
$`\eta=2^{-L}`$, choose

```math
\varepsilon_L=\min\left\{\frac1{10},
                    \frac{\eta}{32\sqrt2K_*}\right\}.
```

Unitary telescoping, using actual adjoints on arbitrary leaked or returned
inputs, adds at most $`8\sqrt2K_*\varepsilon_L\le\eta/4`$.
The retained final real/complex errors are therefore again at most
$`149\eta/512`$ and $`789\eta/1024`$. All source constants must
be uniform in p at this fixed requested accuracy. Equation (19), its
one-clean frame allocation, its small/narrow R96 fallbacks and its
$`\log_2\log_2(n+4)`$ contribution remain unchanged. R96 is still
the strongest unconditional general upper bound.

### 6.2. The next source band retains the general leakage recurrence

The fixed Clifford-sum obstruction applies to two physical source bands
with uniformly bounded clean count and fixed positive block amplitude.
A one-signal source U with **three** physical T bands is a conditional
specialization of the recurrence below, not a constructed next mechanism.
The [current revision](JOINT_SELECTIVE_CONTINUATION_RESULT.md) permits
any qualifying sublogarithmic source depth. Its complete block must satisfy

```math
\|J^\dagger UJ-I_p/\sqrt2\|\le\varepsilon,
\qquad m=2^p.
```

The first unpaid operation is the actual complete logical transition,
with common forward/actual-inverse bounds
$`T_U,b_U=O_\varepsilon(m)`$ and $`G_U=O_\varepsilon(m^2)`$.
An unspecified block dilation does not provide these bounds. The third
band removes the automatic fixed-length Clifford expansion used in the
source-gap proof: a two-band encoder's Pauli images need not be Clifford.
This is a reason that proof no longer applies, not a construction or
evidence that such a transition exists.

Exact flatness must not be imported from Section 6.1. The general
signal-S recovery and the approximate sector cover instead give

```math
\begin{aligned}
\delta_{V_S}&\le3\sqrt2\varepsilon+2\varepsilon^2,\\
\epsilon_{\mathrm{ind}}&\le16\delta_{V_S},\\
T_{\mathrm{ind}}&\le24T_U+48m+8,\\
G_{\mathrm{ind}}&\le24G_U+82m+104,\\
b_{\mathrm{ind}}&\le m+b_U+2,\\
D_T(\mathrm{ind})&\le24d_U+32\le104.
\end{aligned}                                               \qquad(23b)
```

The twenty-four source/actual-inverse calls, positive-S overhead and
even-m removal of eight T gates are unchanged. For
$`0\lt\tau\le1`$, requesting $`\varepsilon\le\tau/72`$
gives $`\epsilon_{\mathrm{ind}}\le\tau`$. The norm is on the
complete arbitrary-input operator, including dirty work and references.
The borrowed signal and toggle need only return approximately here;
their return errors and signal leakage are included in (23b). No reset
or fresh clean input is assumed between calls. The actual inverse has
the same full error bound. If exact flatness is separately proved, the
stronger error and exact return conclusions of (23a) apply instead.

For restored-address R98 set $`\eta=2^{-L}`$ and request

```math
\varepsilon_L=\frac{\eta}{1152K_*}.
```

Four completed indicators per query and at most K_* long-query calls
then add at most $`4K_*\cdot72\varepsilon_L=\eta/4`$.
The retained real and complex total errors are
$`149\eta/512`$ and $`789\eta/1024`$. With uniform source budgets,
the conditional frame depth in (19) would be
$`O_L(N/b^2+\log_2\log_2(n+4))`$, with the same one-clean
allocation, dirty thresholds and R96 fallback branches. This remains
conditional on the unpaid complete transition. A source-class failure
does not establish an unrestricted depth bound, and R96 is unchanged.

### 6.3. Recovering a general clean-bank compression with one extra flag

A growing clean bank requires a paid conversion from a compressed block
to the initialized-isometry promise of Section 1. That conversion is
available at inverse-square-root-two amplitude. It does not construct
the compressed source, or make the later sector cover parallel.
The [paid parallel construction](PARALLEL_CLEAN_PROMOTION.md#51-paying-any-sufficiently-small-packing-permutation)
separately widens the sufficient growing-clean window when its actual
native address-sharing hypothesis holds; Section 5.2 there substitutes
this amplitude-recovery ledger into that complete circuit.

Let A contain k clean wires, with $`k\ge1`$, and let J insert
$`A=0^k`$. All other source inputs, including the complete q-bit
output and private work W, are arbitrary. Suppose an actual native
unitary U satisfies

```math
K=J^\dagger UJ,\qquad
\|K-Q_f/\sqrt2\|\le\varepsilon\le1/8.
```

Introduce one additional clean flag c, untouched by U. Let C_k be the
returned zero-test (4), now with target c and its own arbitrary private
pool, and emit

```math
\widehat R=C_k^\dagger S_c C_k.
```

On the complete c=0 sector this is exactly
$`R=I+(i-1)P`$, where $`P=JJ^\dagger`$ acts on A. Both c and
the zero-test pool return exactly there, for every A, logical and work
input, including superpositions. On c=1 the action is different and is
not discarded. Thus the construction uses k+1 clean wires before
promotion; c is not an implicitly initialized dirty helper.

Write J' for insertion of both A=0 and c=0. With tensor identities
understood, define the actual word

```math
V_{\rm raw}=U\widehat R U^\dagger\widehat R U.
```

Its rightmost factor acts first. Put $`a=i-1`$. The projector algebra,
on the entire initialized input space, gives

```math
V_{\rm raw}J'
 =aJ'K+UJ'\bigl(iI-2iK^\dagger K\bigr).
```

Indeed the more general coefficient of UJ' is
$`1+a+a^2K^\dagger K`$, and $`a^2=-2i`$. Since
$`a/\sqrt2=\omega^3`$ for $`\omega=e^{i\pi/4}`$, the full
initialized-isometry error is

```math
\begin{aligned}
\|K^\dagger K-I/2\|&\le\sqrt2\varepsilon+\varepsilon^2,\\
\|V_{\rm raw}J'-\omega^3J'Q_f\|
 &\le3\sqrt2\varepsilon+2\varepsilon^2
 =:\eta.
\end{aligned}
```

This bound includes rejected A amplitude, address/output error and
private-work return, with arbitrary external references. The flag and
the zero-test pool return exactly during this initialized recovery.
The inverse is the actual reversed native word, and has the same
initialized error against $`\omega^{-3}Q_f^\dagger`$. No claim of
exact A return follows without an additional flatness identity.

The fixed phase already cancels in the promotion of Sections 2--4.
For a standalone literal query, it can instead be corrected explicitly:
on any existing wire,

```math
F=(S^\dagger H)^9=\omega^{-3}I,\qquad V=F V_{\rm raw}.
```

The native identity $`(HS)^3=\omega I`$ fixes this phase. Expanding
each S dagger as three positive S gates costs 36 Cliffords for F;
its actual inverse reduces to $`(HS)^9`$, costing 18. This reduction
only combines adjacent powers of S. With common forward/actual-inverse
source and zero-test bounds, safe common bounds for the corrected V are

```math
\begin{aligned}
T_V&\le3T_U+4T_{C_k},\\
G_V&\le3G_U+4G_{C_k}+42,\\
d_V&\le3d_U+4d_{C_k},\\
b_V&\le b_U+b_{C_k},\qquad k_V=k+1,\\
\|VJ'-J'Q_f\|&\le\eta.
\end{aligned}
```

There are three U/actual-adjoint calls and four C_k/actual-adjoint
calls. The additive 42 allows six positive-S gates for the two inverse
phase gates and 36 for the optional phase correction. The raw word has
the smaller common additive allowance six. All zero-test helpers are
arbitrary and returned, and the displayed width uses disjoint private
pools. Sequential reuse does not multiply width by the call count.

For clarity, keep the explicitly corrected V in the following complete
ledger. Let $`M'=2^{k+1}`$ and let C_{k+1} be the promotion's separate
returned zero-test on the joint A,c bank. Substitution into (11) and
(14) gives a zero-clean full-operator query E with

```math
\begin{aligned}
T_E&\le M'\bigl(12T_U+16T_{C_k}+4T_{C_{k+1}}+24q+4\bigr),\\
G_E&\le M'\bigl(12G_U+16G_{C_k}+4G_{C_{k+1}}+40q+184\bigr)+2q,\\
d_E&\le M'\bigl(12d_U+16d_{C_k}+4d_{C_{k+1}}+16\bigr),\\
b_E&\le k+q+2+b_U+b_{C_k}+b_{C_{k+1}},\\
\|E-Q_f\otimes I_{\rm borrowed}\|
 &\le8M'\eta
 =16\cdot2^k(3\sqrt2\varepsilon+2\varepsilon^2).
\end{aligned}
```

Omitting F from every recovery is legitimate because its phase cancels
in each source conjugation; it replaces 184 by 40 in the Clifford
bound. The displayed conservative ledger nevertheless charges a literal
source correction. Every final borrowed-register return error is
included in the full norm. In particular, the flag's exact return on
the initialized recovery sector is not promoted to a claim of exact
return on arbitrary borrowed superpositions. The actual adjoint of E
obeys the same complete error bound against the inverse query.

If exact flatness $`K^\dagger K=I/2`$ is separately proved, then
$`L=\sqrt2K`$ is unitary and $`VJ'=J'L`$. The joint clean sector
is exactly invariant, so (16) instead gives error at most
$`2\sqrt2\varepsilon`$, with exact return of A, c and the promotion
toggle. This is an additional premise. Section 6.1's one-signal
two-band rounding does not establish it for a growing clean bank.

For $`0\lt\delta\le1`$, the sufficient compression accuracy

```math
\varepsilon\le\frac{\delta}{72\cdot2^k}
```

gives E error at most delta, since
$`3\sqrt2+2\varepsilon\lt9/2`$ in the stated range. This is a
compression accuracy, not an already returned-query accuracy.
For a b-bit indicator source, take

```math
k(b)=\max\{1,\lfloor\beta\log_2\log_2(b+4)\rfloor\},
\qquad \varepsilon_b=\frac{\delta}{72\cdot2^{k(b)}},
\qquad \beta,\gamma\ge0,\quad \beta+\gamma\lt1.
```

If U has depth $`O((\log(b+2))^\gamma)`$ uniformly along this
accuracy schedule, polynomial address overhead in T-count and width,
and at most quadratic-bank Clifford count times polynomial address
overhead, the zero-tests in (10) give the bound below. The simultaneous
count and width bounds must also hold uniformly along the same accuracy
schedule, with fixed polynomial exponents; fixed-epsilon constants alone
do not suffice.

```math
d_E(b)=O\!\left(
 [\log(b+2)]^{\beta+\gamma}
 +[\log(b+2)]^\beta\log\log\log(b+16)
\right)=o(\log(b+2)).
```

The [paid root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
then absorbs the polynomial overhead and multiplies E's complete error
by eight. To obtain final indicator tolerance tau, use
$`\delta=\tau/8`$ above. This pays recovery of a general-k amplitude
block inside the existing sufficient growing-clean window. It supplies
neither such a block nor a faster unrestricted circuit. The first
unpaid operation remains an actual source U meeting its full compressed
operator promise and simultaneous native budgets uniformly along
$`\varepsilon_b`$. Bounded-clean two-band sources retain the
[fixed-sum obstruction](FIXED_CLIFFORD_SUM_SOURCE_GAP.md); that theorem
does not give a gap uniform along this growing-k schedule.

## 7. Where the growing-clean candidate still pays depth

If $`k=O(\log(p+2))`$, the factor $`M=2^k`$ is polynomial
in p. When q and the source count and width bounds are polynomial in p
at the needed tolerance, this gives a polynomial-count, polynomial-width
zero-clean promotion. It is useful where a polynomial allowance is
available for a rare large conjunction. Its complete depth recurrence
remains (11), including M sequential selective words. It is not the
desired sublogarithmic-depth compiler. If k also depends on source
accuracy, both that dependence and (14) must be solved before asserting
any polynomial recurrence.

A tempting shortcut is to regroup the translated cover gate by gate.
For a native elementary gate g, consider the whole translated block
$`\prod_b X_b gX_b`$. Every Clifford g gives a Clifford block.
A T or T-dagger on a wire outside A is repeated M times and is Clifford
because M is even. On a wire of A, half the translations exchange its
two eigenvalues, and

```math
XTX=\omega T^\dagger,\qquad \omega=e^{i\pi/4}.
```

The complete translated block is then the scalar
$`\omega^{M/2}I`$ or its inverse, again a literal Clifford
scalar. Thus regrouping every elementary gate in this way produces
only a Clifford word. It cannot preserve a non-Clifford promoted
query such as a Toffoli query or an indicator with at least two address
bits. The exact sector-cover identity multiplies complete translated
words; it does not license this noncommuting reordering.

This is a failure of that specific regrouping, not an unrestricted
depth lower bound. A different joint native synthesis may share
non-Clifford work across sector words, use noncommuting basis changes,
or avoid the cover. The first unpaid operation for a growing clean
bank is precisely such a complete, fast sector-cover implementation.
For a constant clean bank the next concrete source question is instead
whether a complete approximate indicator with that clean allowance
admits linear T-count and sublogarithmic depth.

Applying the promotion independently to all one-bit indicator rows does
not evade the [native sharing theorem](NATIVE_INDICATOR_SHARING.md).
The useful sufficient interface is a whole-indicator or whole-query
source that shares nonlinear work across its growing output bank.

### An opaque source still requires exponentially many calls in k

The source-call overhead cannot be removed by a generic black-box
wrapper whose nonoracle interleavings are independent of f. This
statement allows arbitrary unitary interleavings and arbitrary mixing
of all address/output/dirty registers; it concerns sequential source
calls, not the native T-depth of a resynthesized circuit.

Consider two exact native source families on the same registers:
$`V_0=I`$ and
$`V_1=(I-P)+PX_T`$, with a single source output bit T.
The latter toggles T precisely on $`A=0^k`$, and can be emitted
by the retained native returned conjunction. Including arbitrary dirty
helpers, let D be the full physical dimension. Then

```math
\mathrm{rank}(V_1-V_0)=\frac{D}{2\cdot2^k}.                 \qquad(24)
```

A fixed wrapper with s appearances of V or its actual inverse has
completed words U_0,U_1 on these two sources. Rank subadditivity
applied to product telescoping gives

```math
\mathrm{rank}(U_1-U_0)\le\frac{sD}{2\cdot2^k}.              \qquad(25)
```

Suppose it promotes both sources with full-operator error at most
$`\tau\lt1`$: the targets are I and X_Z, each with identity
on every borrowed register. The difference U_1 minus U_0 is within
$`2\tau\lt2`$ of $`X_Z-I`$. The latter has D/2 singular
values equal to two, so a perturbation of norm less than two retains
at least D/2 nonzero singular values. Combining with (25) yields

```math
s\ge2^k.                                                   \qquad(26)
```

Both sources have actual inverses equal to themselves, so adjoint
access does not change this argument. Additional arbitrary dirty
workspace multiplies all dimensions equally. Initialized ancillary
subspaces would change the comparison and are outside this zero-clean
full-operator statement.

The bound is specific to source access through a fixed f-independent
wrapper. A compiler that inspects the native source, knows f, and
resynthesizes a joint word is outside its hypotheses. Even on the two
test sources, f is a constant and a compiler given its value can emit
I or X directly. Thus (26) is not a hardness claim for those queries,
an unrestricted T-count bound, or an unrestricted T-depth bound. It
does locate why merely parallelizing an opaque growing-clean cover
cannot be assumed free.

## 8. Exact native receipt and scope

The bounded fixture in
[`test_dirty_promotion_of_clean_queries.py`](../../tests/test_dirty_promotion_of_clean_queries.py)
starts with a one-clean complete AND-XOR query. Its source first applies
a Toffoli, then a controlled-H to an address bit controlled by the clean
wire. The latter is identity on the promised sector but genuinely
nonmonomial outside it. The source has nine T gates and twenty-seven
Cliffords, as does its actual inverse.

The six-wire promoted word uses 128 T gates and 346 Cliffords. Exact
native evaluation checks every column of the selective sector word,
the complete query, and its actual inverse, retaining arbitrary old
target, clean-bank-now-dirty, toggle, and actual output inputs. Literal
phases are compared in $`\mathbb Q(\sqrt2,i)`$. The test verifies
the operator mechanism and a concrete ledger; it does not construct
an asymptotically faster source or validate the approximate theorem by
sampling. Equations (13)–(16) provide the full-operator error proof.

A second three-wire exact native fixture takes
$`U=H_a\otimes\mathrm{CX}(x,T)`$. It checks the complete
$`Q_f/\sqrt2`$ compressed block and every column of the phase-S
recovery. The raw word is literally
$`e^{i\pi/4}S_a^\dagger\otimes Q_f`$, giving the stated clean-code
phase and exact signal return. Its forward eight-Clifford word and
actual twelve-Clifford inverse verify the asymmetric positive-S
accounting used in Section 6. The already checked sector-cover identity
then removes the fixed source phase; no extra phase correction is
inserted into this receipt.

A third five-wire fixture checks Section 6.3 with two initialized source
wires and the additional recovery flag. Its actual two-control zero-test
costs seven T gates and 28 Cliffords. All 32 columns of the raw recovery
and its actual inverse are checked, including a nonmonomial unpromised
sector. The literal phase correction costs 36 forward and 18 inverse
Cliffords; the corrected recovery costs 28 T gates and 156 forward or
142 inverse Cliffords. Its complete initialized action is exactly the
logical CNOT query. This is a bounded native identity receipt, not a
growing-clean source construction or an asymptotic emitter.

The [native sector-parallelism audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md)
separates completed-sector commutativity from an actual parallel schedule.
A clean register fixed at every native Clifford boundary can already be
removed directly. A universal zero-clean black-box wrapper needs at least
one source occurrence per clean sector; that count is not a native-depth
lower bound.
