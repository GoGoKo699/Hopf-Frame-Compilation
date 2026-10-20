# Revision: one algebra source serves the whole indicator bank

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous bank-motion checkpoint](BANK_MOTION_CONTINUATION_RESULT.md)

This continuation starts from main
`9282645cdfa33d37a09a29809306c23498d01b0e`, tree
`af1734f9660d08081ae8a345f6a4a7b746ba40ea`, after PR #206.
The revision changes the candidate mechanism without changing the general
destination. The retained HRT algebra permutation can supply the complete
indicator bank on one shared arbitrary coefficient carrier. Its count and
width now fit the polynomial-overhead source window. Its known native
depth remains logarithmic, so a faster source is still unpaid.

## 1. Separate the destination from the current results

| Role | Current assessment |
| --- | --- |
| General destination | Optimal unrestricted native T-depth with simultaneous T-count, Clifford count, dirty width, complete error and return |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md), with selected [R79](UNIFORM_SUBTREE_DEPTH.md) unchanged |
| Paid transfer | [Root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) absorbs polynomial address overhead using eight complete source/actual-inverse calls and an exact short suffix |
| Reopened constructive component | [Joint HRT extraction](FULL_BANK_CHARACTER_TRANSPORT.md) uses one coefficient carrier, eight algebra-source calls and a paid Clifford zeta readout for every arbitrary output simultaneously |
| Remaining constructive cost | The retained source has logarithmic depth. Its classical bounded-degree-stage limit does not constrain arbitrary native interference or synthesis of the complete eight-call word |
| Clean-promotion assessment | [Sector-parallelism audit](CLEAN_SECTOR_PARALLELISM_AUDIT.md): the simple native read-only condition that permits parallel covers also makes those clean bits dispensable. General black-box promotion requires at least one source occurrence per clean sector; this is a query-count statement |
| Focused lower scout | [Original singleton coupling](ORIGINAL_PAIR_SINGLETON_COUPLING.md) retains all positive middle singletons and derives a transported center at bounded actual singleton-coupling rank. No rank bound from one physical pair alone is inferred |

For $`N=2^n`$, integer $`n\ge1,L\ge6`$, R96 still gives

```math
D_T=O(NL/B^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
```

There is one clean qubit, real/complex dirty reservations respectively
$`B\ge17(L+n+7)`$ and $`B\ge17(L+n+8)`$, and complete
initialized-isometry error at most $`2^{-L}`$. The fixed-accuracy,
square-root-width unrestricted depth interval remains logarithmic versus
constant. None of the representation tests below closes that interval.

## 2. What the full-bank algebra extraction now pays

Put $`m=2^b`$, $`r=b+1`$, and use the square-zero algebra

```math
R_r=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2).
```

Every nonempty coefficient is a physical arbitrary bit. The retained
HRT permutation F maps coefficient c to

```math
f(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1.
```

Its carrier has $`2m-1`$ bits. The retained scalar extraction used a
fixed table-dependent carrier translation. Replace that translation by
one linear function of the entire arbitrary output bank Y. The coefficient
for the complementary carrier label $`[r]\setminus S`$, for each
nonempty $`S\subseteq[b]`$, is

```math
\sum_{u\subseteq S}Y_u.
```

These are precisely the algebraic-normal-form coefficients of the selected
bank bit $`Y_x`$. Compute the binary zeta transform in place on Y,
copy its nonempty coordinates to the carrier, and actually reverse the
transform. This exact Clifford translation costs $`bm+m-1`$ CNOTs;
Y returns before the next operation. No copy of an initialized output,
clean coefficient or separate source per row is used.

The two retained commutators, the constant-output correction and the
output Hadamards give the complete indicator. For common forward and
actual-inverse source bounds and full-operator source error eta,

```math
\begin{aligned}
T_E&\le8T_F,&D_T(E)&\le8d_F,\\
G_E&\le8G_F+2bm+4b+4m+10,&
w_E&\le2m-1+w_F,\\
\epsilon_E&\le8\eta,&a_{\rm clean}&=0.
\end{aligned}
```

The full operator comparison includes every address, output, coefficient,
private-work and reference input. Actual emitted inverses are used.
Approximate calls need not return their work exactly between occurrences;
the eight-call telescoping estimate pays every such error.

The [existing fast source](ALGEBRA_QUERY_EXTRACTION.md#5-a-smaller-count-complete-hrt-source)
has $`T_F,G_F=O(r^2 2^r\log(r+2))`$,
$`w_F=O(r^2 2^r)`$ and $`d_F=O(\log(r+2))`$, exactly.
Thus this complete bank already satisfies the relaxed source count and
width budgets. The old square-root-width objection for a single table
does not stop this full-bank source. The known depth still fails to improve
the retained frontier.

## 3. The first unpaid operation and the next constructive question

The concrete next question is a full-input native implementation of this
specific HRT permutation F, or its actual inverse, with the same polynomial
overhead and depth $`o(\log(r+2))`$ at every fixed requested full-operator
tolerance. A complete eight-call synthesis can also succeed without
separately approximating F. The separate-source recurrence is sufficient,
not a necessary form for a solution.

There is now a paid shallow component for its quadratic inverse k.
The complete query $`(a,z)\mapsto(a,z\mathbin\oplus k(a))`$ has at
most 40 T layers, $`O(r^3 2^r)`$ T/Clifford count and dirty width,
with every helper returned. The source and output banks are arbitrary.
This does not implement k in place: z still contains its original unknown
value, and a is still present. Swapping the registers does not erase either
obligation. A full-input conversion or a different coherent mechanism must
be supplied before this component can improve the frontier.

The same query pays a complete Choi-basis transform on two arbitrary
registers and the coherent Z correction. The next specific missing gate
is the coherent X correction $`R_z=F_rX^zF_r^\dagger`$. It acts as
$`y\mapsto y+f(z)+yf(z)`$ on an arbitrary carrier. A fixed key gives
a Clifford, but that does not synthesize the coherent family. Its
arbitrary-key echo gives another complete indicator recurrence:

```math
\begin{aligned}
T_E&\le8T_R,&d_E&\le8d_R,\\
G_E&\le8G_R+2bm+8b+4m+10,&w_E&\le4m-2+w_R,\\
\epsilon_E&\le8\eta,&a_{\rm clean}&=0.
\end{aligned}
```

This reduction pays the key return with actual adjoints and includes
all arbitrary inputs and references. A clean-Bell teleportation calculation
does not supply the missing native correction or a returned-work word
for arbitrary Bell labels. The full R source, or a different complete
fused mechanism, is now the concrete next question.

The current balanced arithmetic compiler does not answer that question.
F has degree r, while its inverse is quadratic. An exact full-input
implementation of that inverse whose completed classical inverse stages
have degree at most delta
therefore needs at least $`\lceil\log_\delta r\rceil`$ such stages,
even if dirty garbage remains live. These are completed classical stages,
not native T layers. Nonmonomial intermediate interference and a jointly
synthesized extraction word lie outside that argument. Third-level
Clifford-hierarchy membership also does not supply a native schedule.

Combining a qualifying F source with root fusion gives

```math
\epsilon_{\rm normalized}\le64\eta,\qquad
D_T\le64d_F(b_0+1)
       +O_\tau(\log\log(p+4)).
```

Here $`b_0=p-s`$ and $`s=\lceil a\log_2(p+2)\rceil`$ are the
actual large-prefix and exact-suffix lengths chosen by root fusion,
with fixed a large enough to absorb the source overhead.

Requesting $`\eta=\tau/64`$ is enough for the normalized indicator.
For the fixed-L frame transfer, request
$`\eta_L=2^{-L}/(1024K_*)`$, where $`K_*`$ is the retained
bounded query count. The frame error and literal one-clean allocations
are then those of the paid transfer. Substituting the known logarithmic
F source keeps the logarithmic term; it does not change R96.

## 4. Why the simple clean-cover parallelization is not a new source

If every actual Clifford boundary preserves each computational sector of
the clean A register, and no T gate targets A, restriction to A equal to
zero gives an actual Clifford+T word on the remaining arbitrary registers.
Its T-count and depth are unchanged. The restriction of a sector-preserving
Clifford is a Clifford, including its literal scalar, so the apparently
useful clean source already has a zero-clean implementation with paid
Clifford resynthesis.

Keeping A coherent between boundaries removes that argument, but also
removes the proposed parallel schedule: endpoint sector identities do not
make overlapping native T targets disjoint. This is the first unsupported
inference in that proposed shortcut.

A separate full-space rank argument shows that a universal black-box
promotion must use at least $`2^k`$ source or inverse occurrences for
k clean source bits, at output errors whose sum is below two. It permits
arbitrary source-independent interlayers and dirty workspace. It counts
occurrences, including parallel calls, rather than native depth. Inspecting
and resynthesizing a known source is outside this oracle statement.

## 5. Lower scout and verification boundary

The original all-active genuine one-pair scout now retains positive
middle singletons. Its new parameter h is the actual binary pairing rank
between the singleton span and either endpoint silent span; the two ranks
are proved equal. The physical endpoint transport identifies an enlarged
center with residual dimension at most $`2h+4`$. The shared-clean feature
proof then gives sign degree at most $`D_{2h+4}+4A_{2h+4}`$, and XOR
degree one less. The original higher-degree rounding floor remains
$`2^{-530690}`$ under this theorem's structural hypotheses.

This bounds a genuinely transported quotient, not an axis span inferred
from pairing rank alone. All actual middle gates, the unused clean block,
literal phases and the physical clean coordinate are retained. The sharper
one-sided bounds and the recovered pair-only constant remain available
in their own classes.
Unbounded h and general pair differences remain open; no unrestricted
one-clean three-layer degree theorem is asserted.

The new focused check verifies the joint zeta translation and phase
extraction on all small address, bank and coefficient columns, including
the actual inverse. It is a composition witness using the retained native
source proof, not a new large native emitter or an asymptotic proof.
The clean-sector and singleton-coupling arguments are analytic. Required
validation and exact-tree publication receipts belong on the pull request.
