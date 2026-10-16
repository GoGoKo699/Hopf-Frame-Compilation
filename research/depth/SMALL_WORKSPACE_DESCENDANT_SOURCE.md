# Descendant sources with logarithmic initialized scratch

[Signed descendant source](COHERENT_DESCENDANT_LOADER.md) ·
[Borrowed real rotations](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) ·
[Exact parallel queries](PARALLEL_BILINEAR_WORD_QUERY.md)

The signed descendant source can move its precision-sized program and
gradient into arbitrary dirty work. Its only initialized source scratch
is a logarithmic-size expansion index and one attenuation bit. The
replacement uses a full-unitary borrowed-signal rotation estimate; it
does not initialize an occupied phase register or infer a full-space
contract from a single prepared column.

This is a source component. A complete frame, spectral-label decoder,
and any accompanying depth improvement require their own composition,
protected logical reserve, and resource arguments.

## 1. Source statement

Retain the source labels of the
[signed descendant construction](COHERENT_DESCENDANT_LOADER.md#1-source-contract-and-disjoint-descendant-sectors).
An h-bit logical register X has at most h+1 valid tags k. Within tag k,
its first d(k) bits are an unchanged node prefix, and its remaining
s(k)=h-d(k) bits start at zero. For each such label a, prescribe a
computably evaluable real unit vector $`\psi_a`$ on those remaining
bits, with its literal sign. An optional immutable external prefix P
may also select this entire family.

Let $`0\lt\delta\le1/64`$ and put

```math
p=\lceil\log_2(1/\delta)\rceil,
\qquad
J=\min\{2^q:2(4/5)^{2^q}\le\delta/2\}.
```

Then $`J=O(p+1)`$ and $`q=\log_2J=O(\log(p+2))`$, with
absolute constants. There is a native source word U using exactly q
initialized expansion-index bits and one initialized attenuation bit,
besides the local zero input supplied by each source anchor. All other
additional wires may be arbitrary dirty inputs. Writing A for this
(q+1)-bit initialized scratch and D for the complete arbitrary dirty
allocation, the source has

```math
\|UJ_{\rm src}-\Psi_{\rm src}\|\le\delta,
```

```math
\Psi_{\rm src}:
|P,k,u\rangle|0^{s(k)}\rangle|0^A\rangle|\xi\rangle_D
\longmapsto
|P,k,u\rangle|\psi_{P,k,u}\rangle|0^A\rangle|\xi\rangle_D.
\qquad(1)
```

The map $`J_{\rm src}`$ embeds arbitrary superpositions of these
anchors and arbitrary dirty and reference inputs. Tags, P, and the
within-tag prefix u are preserved exactly by U. The norm includes every
source-work and dirty-work error; actual dirty core return is not assumed
exact.

If a preserved activity bit z is supplied, U can be made exactly
identity on its entire z-zero sector, including arbitrary occupied A,
local target, and D inputs. Only the z-one source anchors need A zero.
The activity bit is an unchanged control, not an extra source
initialization. In a complete-frame application it must be supplied by
the separately charged saved-predicate construction.

For the standard descendant family without P, one fixed absolute
constant C gives the sufficient standalone resource bounds

```math
\begin{aligned}
T,w_D&\le C(h+p+3)^3 2^{h/2},\\
G&\le C(h+p+3)^3 2^h,\\
D_T&\le C[p+\log_2(h+p+3)].
\end{aligned}                                                    \qquad(2)
```

Here $`w_D`$ counts all arbitrary dirty signals, precision cores, mask
outputs, and query helpers; A and the original logical/tag registers are
separate. Elementary Clifford gates are counted. The constant is fixed
by the retained native primitives and is independent of the target data
and requested accuracy. These deliberately generous bounds are not
claimed count optimal. Section 4 gives the address ledger when an
external prefix P is present.

## 2. Replace calibration by a borrowed real rotation

The existing geometric expansion and its literal five-call inner
amplification are unchanged. They use $`\gamma=3/5`$ and
$`\beta=4/5`$, with strictly certified Boolean-phase states
$`\phi_{a,j}`$, and satisfy

```math
\left\|\psi_a-\gamma\sum_{j=0}^{J-1}\beta^j\phi_{a,j}\right\|
\lt\beta^J.
```

Its normalized finite sum $`\Phi_a`$ therefore obeys
$`\|\Phi_a-\psi_a\|\le2\beta^J\le\delta/2`$.
The coefficient-index preparation is the real product

```math
C_\chi=\bigotimes_{j=0}^{q-1}R_y(\theta_j),
\qquad \theta_j=\arctan(\beta^{2^j/2}).
```

Keep the existing parallel native one-qubit approximation and guarded
half-angle echo for this bank. Allocate total operator error at most
$`\delta/64`$ to each complete coefficient-bank appearance, by
using at most $`\delta/(128q)`$ error in each half-angle word.
Its count is $`O(q[p+\log(q+2)])`$ and its T-depth is
$`O(p+\log(q+2))`$. It needs no ancillary initialized bits beyond
the q targets that are already part of A. Its actual inverse has the
same full-operator error. The z-zero sector is exact identity because
the native preparation and its actual inverse cancel there.

The calibration target t is the sole attenuation bit. Its ideal
label-dependent real rotation has zero-input coefficient

```math
\cos\vartheta_a=\sin(\pi/10)/\xi_a,
\qquad 1/3\le\xi_a\le1.
```

These are the same computable calibration data as in the original
source. The strict margin $`\sin(\pi/10)\lt1/3`$ makes the
calibration uniformly well conditioned. No exact zero or analytic sign
test is required.

Use the
[borrowed-signal real-rotation word](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
at scalar precision

```math
q_{\rm cal}=p+12.
```

Its borrowed signal r and its $`q_{\rm cal}+1`$-bit Majorana core
belong to D. Neither is initialized. The unchanged logical label a is
the table address, t is the target, and z is the optional one-literal
activity predicate. The complete actual word $`\widehat C`$ has

```math
\|\widehat C-C_{\rm cal}\otimes I_D\|
\lt43\,2^{-q_{\rm cal}}
\le\frac{43}{4096}\delta
\lt\frac{\delta}{64}.                                  \qquad(3)
```

Here $`C_{\rm cal}`$ is the addressed real rotation on the active
sector and identity on the inactive sector. Equation (3) is a
full-unitary estimate on all attenuation, borrowed-signal, core, and
helper inputs. In particular it does not require a clean phase compiler
flag. The exact flag symmetry of the real word is essential; the
one-clean literal-phase word by itself would not give this estimate.
The actual inverse obeys the same bound for
$`C_{\rm cal}^\dagger\otimes I_D`$.

Replace each completed mask in that scalar word by the exact R91 XOR
query, with the Z masks obtained by full-output Hadamard conjugation.
Keep the literal order of its X and Z masks and use actual native
inverses. Each replacement is the same complete operator tensored with
identity on its disjoint returned query pool. Thus it preserves (3),
the rejected operator, the exact flag symmetry, and exact z-zero
identity. The source loaders remain unconditional; only their central
X gates receive the one-literal guard. Their total T-depth is
$`O(p)`$, and each guarded center has constant native depth.

The former initialized calibration-program and gradient registers B,F
are absent. Their replacement is a scalar precision core with arbitrary
input, governed by the full-space estimate (3).

## 3. The complete source error and the larger completion space

Boolean phases, variable-suffix Hadamards, and inner reflections retain
their exact returned-work constructions. The inner reflection tests only
the source index, attenuation bit, and the appropriate local zero input.
It does not test the borrowed signal or precision core. That is exactly
the ideal source reflection: these dirty systems are identity factors in
both desired coefficient and calibration operations.

Compare the actual word with the ideal geometric-source word tensored
with $`I_D`$, replacing each complete coefficient bank and each complete
calibration by its ideal operation. The expanded two-round inner
amplification has five flagging-word appearances, each with two
coefficient banks and one calibration. Every actual inverse is the
reversed native word. The full-unitary estimates give, by a fifteen-term
hybrid,

```math
\|U-U_{\rm geom}\otimes I_D\|
\le15\delta/64.
```

This comparison does not require actual restoration of D between calls.
Earlier errors propagate through later unitaries. On initialized source
anchors, the ideal word returns A and prepares $`\Phi_a`$ with its
positive literal amplification phase. Consequently

```math
\|UJ_{\rm src}-\Psi_{\rm src}\|
\le\delta/2+15\delta/64
=47\delta/64\lt\delta.                                \qquad(4)
```

Preserved labels give a direct-sum maximum, not a sum over the number of
source columns. Tensoring arbitrary references preserves all estimates.
On z zero, each completed native macro is exactly identity, including
the correctly guarded reflection phases, so the complete source is
exactly identity on arbitrary occupied source and dirty work.

There is one change to the use of the
[dimension-free source-completion lemma](SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison).
The actual source need not factor as a mode-space unitary tensored with
$`I_D`$, since (3) permits small core disturbance. Apply the lemma
instead on the whole physical space, to the initialization isometry
$`J_{\rm src}`$ that already contains $`I_D`$. It supplies a
mathematical unitary correction C with

```math
\|C-I\|\le6\delta,\qquad U_*=CU,\qquad
U_*J_{\rm src}=\Psi_{\rm src},\qquad
\|U-U_*\|\le6\delta.                                  \qquad(5)
```

Construct this correction separately within the exactly preserved
label sectors and extend it by identity on z zero and unused sectors.
It therefore retains the preserved labels and inactive identity. It is
never physically implemented. Actual inverse words have the same
full-space comparison error.

In particular, let $`\Pi_{\rm src}=J_{\rm src}J_{\rm src}^\dagger`$.
Its physical source-anchor predicate includes validity of the padded tag;
invalid tags are excluded. All query tables retain a fixed total definition
on invalid tags, as in the original descendant loader, so no unspecified
padding sector is included in this projector.
Equation (5) gives the exact comparison-projector identity

```math
U_*\Pi_{\rm src}U_*^\dagger
=\Psi_{\rm src}\Psi_{\rm src}^\dagger.                 \qquad(6)
```

The right side has the prescribed descendant projectors and is exactly
identity on D. This is the property needed when a later argument
transports source-anchor reflections. It does not assume a chosen
complete action of $`U_*`$ on the complementary columns, and it does
not promote the actual source's approximate dirty return to an exact
one. Any later reference/flag construction must use (5)–(6) on this
larger space and charge its own finite number of source occurrences.

## 4. Native table and workspace ledger

First omit P. Let $`s=\lceil\log_2(h+1)\rceil`$ be the padded
tag length. A Boolean source-phase address contains the h logical bits,
s tag bits, q expansion-index bits, and the optional activity bit.
Thus a common padded row envelope is

```math
R\le C_0(h+1)(p+1)2^h,
\qquad a_{\max}=\lceil\log_2R\rceil
\le h+C_1\log_2(h+p+3).                                \qquad(7)
```

The calibration uses the tag and current logical register as a padded
address; its table ignores every descendant-target bit beyond the
preserved prefix. Its scalar mask output has
$`q_{\rm cal}+1=p+13`$ bits. Boolean phases have one output bit.
Variable-suffix Hadamard masks have at most $`2^{s+1}`$ rows and h
outputs. Invalid tags can be assigned identity rows. The constant number
of table positions includes the two Boolean source phases, calibration
masks, variable Hadamard masks, inner reflections, and every actual
inverse in the five-call expansion.

Choose one fixed dirty query-pool size

```math
B=\left\lceil C_2\left[a_{\max}+h+p+3+
              \sqrt{R(p+13)+h(h+1)}\right]\right\rceil, \qquad(8)
```

where $`C_2`$ exceeds the absolute R91 minimum-width constant. The
pool is disjoint from the live precision core, borrowed signal,
logical/label bits, initialized A, and the separate query outputs. It
may be reused only after each completed exact query returns it. Two
returned predicate helpers, when used instead of a reflection table,
are reserved separately or borrowed from this returned pool between
queries.

For every large query, R91 gives

```math
T_Q=O(\sqrt{Rm}),\qquad G_Q=O(Rm),\qquad
D_{T,Q}=O(\log_2(a_{\max}+2)+1),
```

where $`m\le p+13`$. Equation (8) also covers the short h-output
Hadamard-mask queries. The one-literal guarded scalar centers have
constant depth. Source loaders contribute $`O(p)`$ count and depth;
the parallel coefficient banks contribute
$`O(q[p+\log(q+2)])`$ count and
$`O(p+\log(q+2))`$ depth. All fixed routing, target basis changes,
and native scalar reflections are charged. Equations (7)–(8) then give
(2), after one absolute enlargement of C. The displayed cubic majorant
covers the precision core, signal, short masks, and standalone
polynomial terms as well as the query pool.

If an immutable P of length t is present, append it only where the data
depend on it. A sufficient large-table envelope is

```math
R_P\le C_0\,2^t(h+1)(p+1)2^h,
\qquad
a_{\max,P}\le t+h+C_1\log_2(h+p+3).                  \qquad(9)
```

Replace R and $`a_{\max}`$ by these quantities in (8) and in the
query ledger. Fixed coefficient banks and short Hadamard masks do not
acquire P addresses. The source-isometry error is still (4), with no
$`2^t`$ factor. In particular, at a sufficient private pool the native
depth is

```math
O\bigl(p+\log_2(t+h+p+3)\bigr).
```

No global count or width absorption is inferred merely from this local
statement: every application must retain the actual row factors in (9),
its exclusive query pools, and all simultaneously live source banks.
All constants and polynomial degrees in this lemma are fixed before a
caller chooses its protected logical reserve.

## 5. Scope

This construction prepares signed descendant columns with arbitrary
dirty inputs and logarithmic initialized scratch. It also supplies the
full-space comparison (5) needed for transported projectors. It does
not state that source columns determine an arbitrary completion, that a
prepared root alone determines an entire nested flag, or that source
words on shared targets can run in parallel. Repeated noncommuting
source calls still have to be scheduled and charged.

The proof changes one precision interface in an existing source. Its
geometric residual recursion, actual inner amplification, literal
inactive phases, borrowed-signal symmetry, and complete exact query
identities retain their cited proof homes. No QRAM, measurement, reset,
supplied precision state, efficient-preprocessing claim, or scalable
emitter is introduced.
