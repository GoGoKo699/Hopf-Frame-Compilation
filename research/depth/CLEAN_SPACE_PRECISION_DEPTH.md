# A clean-space tradeoff for large precision

[All-precision bound R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Exact key-cover bound R104](EXACT_KEY_COVER_PRECISION.md) ·
[Precision count–depth lower bound R103](PRECISION_COUNT_DEPTH_TRADEOFF.md)

**Theorem R105.** Fix a sufficiently small absolute $`c\gt0`$.
For integers $`n\ge1`$, $`L\ge6`$ and $`L\ge2^{cn}`$, put
$`N=2^n`$ and choose an integer $`1\le t\le L`$.
Every prescribed real complete Hopf frame and prescribed complex magnitude
completion has a native Clifford+T word and its actual reverse with
initialized-isometry error $`2^{-L}`$, at most $`1+C_{\rm cl}(t+1)`$
external clean inputs, and the literal R96 dirty thresholds
$`b\ge17(L+n+7)`$ and $`b\ge17(L+n+8)`$, respectively.
For a fixed polynomial $`P_{\rm c}`$, the simultaneous bounds are

```math
\begin{aligned}
T&=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\\
D_T&=O\!\left(NL/b^2+(L/t+P_{\rm c}(n))\log_2(t+2)\right).
\end{aligned}
```

For $`1\le t\le L^{3/4}`$, the depth is
$`O(NL/b^2+L\log_2(t+2)/t)`$. At $`t=1`$ R96 gives this
bound with one clean input. All additional clean inputs return
exactly. The error includes the original clean output, all logical
columns, arbitrary dirty inputs and references.

For $`L\ge N`$ at the stated width, $`T=O(L)`$ and $`NL/b^2=O(1)`$.
R103 gives worst-case depth $`\Omega(L/(k+1))`$ with $`k`$ external
clean inputs at that T count. For any sufficiently large fixed u,
taking $`t=\max\{1,\lfloor L/\log_2^u(L+2)\rfloor\}`$ uses
$`k=O(L/\log_2^u(L+2))`$ clean inputs and
$`D_T=O(\log_2^{u+1}(L+2))`$. This is within a logarithmic factor
of the lower bound at that clean allocation.

## 1. A linear-size menu with a shared comb selector

Let

```math
G_p=\sum_{j\lt p-1}2^{-(j+1)/2}\Gamma_j+
 2^{-(p-1)/2}\Gamma_{p-1},
```

where $`p\ge2`$ and the Hermitian Pauli axes anticommute. The R100
half-amplitude menu can use $`O(p)`$ atoms and a geometric prefix
selector. Write $`p=2r`$ or $`p=2r+1`$. The first uniform selector
bit chooses two equally weighted groups.

* The Pauli group gives every even-exponent coefficient its dyadic
  weight on its Pauli. For $`p=2r`$ these weights are
  $`2^{-1},\ldots,2^{-(r-1)}`$. Pad the remaining $`2^{-(r-1)}`$
  by opposite identity atoms of weight $`2^{-r}`$ each. For
  $`p=2r+1`$ the weights are $`2^{-1},\ldots,2^{-r},2^{-r}`$.
* In the pair group, each odd-exponent coefficient other than
  the first is $`z_j/\sqrt2`$. Put weight $`z_j`$ on
  $`(\Gamma_0+\Gamma_j)/\sqrt2`$. For $`p=2r`$ these weights
  are $`2^{-1},\ldots,2^{-(r-1)},2^{-(r-1)}`$, with only
  one weight-one atom for $`p=2`$. For $`p=2r+1`$, the
  weights $`2^{-1},\ldots,2^{-(r-1)}`$ leave $`2^{-(r-1)}`$.
  Divide this between $`(\Gamma_0\pm\Gamma_1)/\sqrt2`$.

The pair group's coefficient of $`\Gamma_0`$ is $`1/\sqrt2`$.
Every other coefficient is supplied once apart from the cancelling
padding. Including the unbiased group bit, the weighted operator is
$`G_p/2`$ with literal phases.

A comb of depth w uses disjoint prefixes
$`0,10,110,\ldots,1^{w-1}0,1^w`$ of a shared uniform tail.
Use the longer tail and ignore the surplus bits in the shorter
comb. At depth zero the only prefix is empty. Compute all prefix
conjunctions $`u_i=x_1\cdots x_i`$ with a balanced product tree and
a downward prefix pass on initialized wires. The left child inherits
the parent's exclusive prefix, while the right child's prefix is
that prefix times the stored product of the left child. Initialize
the root prefix to one and pad to a power of two with ones. This
uses $`O(p)`$ Toffoli and Clifford gates, $`O(p)`$ clean wires, and
$`O(\log_2(p+2))`$ Toffoli depth. Within each tree level the Toffoli
supports are disjoint. Each atom flag is the selected group, the
corresponding prefix and, if needed, its next zero literal. Use a
constant number of additional gates per atom and clean copies of any
shared inputs by CNOT and its actual reverse. The flags are one-hot
for every selector value. Reverse all flag computation after SELECT.

## 2. Multiplexed rotations with one-hot scratch

A pair atom is $`R_j\Gamma_0R_j^\dagger`$ with
$`R_j=e^{\pm i\pi P_j/8}`$ and $`P_j=i\Gamma_0\Gamma_j`$.
The plus exponential gives the plus pair. Let $`g_j`$ be the flag
and $`b_j`$ one initialized scratch bit per pair. Define E by
the following chronological layers.

1. Apply $`H_{b_j}`$ controlled by $`g_j`$ for every pair j.
2. Apply the Clifford product of $`P_j`$ controlled by $`b_j`$.
3. Repeat the controlled H batch.

Each controlled H uses two T gates by conjugating controlled Z with
a Clifford-conjugated native T or T-dagger on its target. The scalar
cancels and the supports in a batch are disjoint. On a selected
pair branch only $`b_j`$ enters step 2 in $`|+\rangle`$.
The third step writes the spectral projectors of $`P_j`$ into
that scratch, as $`|0\rangle\Pi_+(P_j)+|1\rangle\Pi_-(P_j)`$.
A parallel physical T/T-dagger layer on the scratch followed by
the actual reverse of E applies
$`\omega^{\pm(I-P_j)/2}`$ and returns all $`b_j`$.
On nonpair branches the scratch stays zero and the action is
identity. The full identity includes coherent selector and arbitrary
data and reference.

Let B denote the resulting rotation, choosing T-dagger for the plus
exponential. SELECT is B times the selected Pauli center times
$`B^\dagger`$. A pair branch selects $`\Gamma_0`$ and the
literal rotation scalar cancels. Other branches select their original
Pauli by a controlled Pauli Clifford or signed identity by a sign on
the atom flag. All scratch returns for arbitrary selector and data
after reversing the flags. SELECT uses $`O(p)`$ T gates, $`O(p)`$
initialized work and $`O(\log_2(p+2))`$ T-depth. Controlled Paulis
cost $`O(p^2)`$ elementary Cliffords on an $`O(p)`$-line core by
direct expansion. Section 4 has linear cost for the actual geometry.

## 3. Exact source replacement

Conjugate SELECT by Hadamards on the uniform selector to get U.
Its compression to zero selector is $`A=G_p/2`$, and the helpers
return on every selector value. For $`R=2\Pi_0-I`$ the zero-selector
reflection, the compression of $`URU^\dagger RU`$ is
$`4AA^\dagger A-3A=-G_p`$. Its unitarity ensures exact return on the
full initialized code. A clean balanced AND tree implements the sign
of R with $`O(p)`$ T and Cliffords and $`O(\log_2(p+2))`$ T-depth.
Flipping both reflection signs leaves the product unchanged.

Complete the underlying paired Majorana list in the original core.
Use its chirality Pauli times $`X_{\rm spare}`$ as $`\Gamma_a`$
on one arbitrary spare line. For a chunk of at most t splits,
let $`\Gamma_j`$ be its entrance axis and $`M_j`$ its output.
Then $`G'=(\Gamma_a+M_j)/\sqrt2`$ has the required geometric list,
with one leading axis. Apply $`-G'`$ on the clean code with the
fixed $`\Gamma_a`$ and plane Clifford
$`C_j=\exp(\pi\Gamma_a\Gamma_j/4)`$. The data action
$`-G'\Gamma_a C_j`$ conjugates $`\Gamma_j`$ to $`M_j`$.
It fixes all outside original axes since the plane commutes
with them and the other factors both negate them.

Repeat on both tails as needed, retaining the original constant seed
splits. Reuse the clean workspace. The comparison V transports
$`X_0`$ to the old reflection M on the enlarged arbitrary core.
Substitute $`V C_F(X_0)V^\dagger`$ at each complete original source,
with the same stable external guard F in the center. The substitute
agrees on the full added clean code, including arbitrary occupied
original work and references. On F false the actual sandwich is
literal identity even with occupied extra work. Retain the exact
controlled Pauli sign mask for the conjugate-signal lift. This
preserves the original rejection, mixed gauges, guards, actual
inverses and error budget.

## 4. Linear elementary Cliffords for the actual geometry

Let m be the old physical core size and $`J=[j,h)`$ the chunk's
contiguous interval of $`O(t+1)`$ old wires. If the prefix is
nonempty, use physical 0 as an outside pivot and apply CNOTs
from $`1\le r\lt j`$ into it. For a nonempty suffix use
physical $`m-1`$, with CNOTs from $`h\le r\lt m-1`$.
These take $`Z_{\lt j}`$ and $`Z_{\ge h}`$ to their respective
single Z pivots. The CNOT sets commute, as pivot targets
never serve as sources. Between successive chunks on one tail,
only the symmetric difference of the sets is needed. Its endpoints
are monotone, so the initial frame, all updates and the final
unframing cost $`O(m)`$ Cliffords in total. Resetting between two
tails costs another $`O(m)`$. This rebrackets exact frame
conjugations with actual inverses.

An original axis at physical r has form $`Q_{\rm pre}Z_{[j,r)}B_r`$,
where $`Q_{\rm pre}`$ is the prefix pivot and $`B_r`$ its on-site
X or Y. The spare axis is both outside pivots times $`Z_J`$
and the spare X, up to its known sign. Because the first axis of
$`G'`$ is this spare, each pair-rotation Pauli $`i\Gamma_a\Gamma_j`$
has a common outside-and-spare factor, a strict suffix in J,
and its signed on-site complementary X or Y.

For zero or one-hot computational controls, implement each of
these controlled lists as follows. Apply the local controlled
Paulis with their known signs. Apply the common outside factor
controlled by the parity of the controls, computed into and
uncomputed from a clean bit. For strict suffixes sweep J left
to right with a clean parity bit p. At each wire r apply CZ
from p to r, then toggle p by every control located at r.
Undo only the p toggles at the end, restoring zero. Sweep
right to left for strict prefixes. On each allowed label these
are the desired literal single Pauli's factors on disjoint
supports. The word is Clifford with $`O(|J|+\#\mathrm{controls})`$
elementary gates. It preserves the computational controls, so
its actual inverse agrees on the same code.

At step 2 of extraction, the b bank is supported on zero and one-hot
labels on every initialized selector branch. Thus the suffix sweep
has the required action. Direct Pauli atoms use the prefix sweep
with their one-hot flags. The common center $`\Gamma_a`$ is controlled
by the parity of the pair flags with $`O(|J|+1)`$ Cliffords.
The signed identity flags take their literal phases. All helpers
return. The fixed plane Clifford and spare Pauli have length
$`O(|J|+1)`$ in the compressed frame. A complete source with its
actual reverse costs $`O(m+1)`$ T and elementary Cliffords and
$`O((1+m/t)\log_2(t+2))`$ T-depth. The extra peak is $`O(t+1)`$
initialized wires and one arbitrary spare line beyond the
inherited center guard.

## 5. All-precision splice and the clean-depth regime

Use the R104 constant-group direct-prefix assembly, including
the unchanged R91 packed queries, R72 bottom, and R72 final
complex diagonal. At $`L\ge2^{cn}`$ the input height is
$`O(\log_2(L+2))`$. Choose the same small fixed group height
fraction as R104, with group heights bounded by
$`\max\{1,\lfloor c_h\log_2(L+2)\rfloor\}`$. The number of groups
is bounded. The local high-precision condition, literal width,
and retained query and coarse schedules hold for $`L\gt N`$ as well.
The direct-prefix component only requires the local height bound.
Use bounded-height R72 where required, including $`n\le3`$.
The complete source positions number at most $`P(n)`$ for a fixed
polynomial P. Their retained widths sum to $`O(L+P(n))=O(L)`$,
including actual inverse visits.

R104's residual native depth outside these sources is
$`O(NL/b^2+L^{1/8}\log_2(n+2)+P(n))`$. The old precision base,
coarse work and exclusive R91 pool have room for one spare line
at the literal threshold beyond one absolute cutoff. Cases of
bounded L have bounded n and can use R96. Every mask, query,
coarse and guard count stays as in the same direct-prefix
splice. The replacements cost $`O(L)`$ total T and elementary
Cliffords, keeping $`G=O(NL)`$ for $`L\gg N`$ as well. Their
native T-depth is

```math
O((L/t+P(n))\log_2(t+2)).
```

For $`t\le L^{3/4}`$, the $`L/t`$ term absorbs every fixed
polynomial in n and the coarse depth. Section 6 reschedules the
coarse work for larger t. The real part and final
complex phase use the original shift to $`L+1`$ at the same t.
The shifted literal real width matches the complex threshold.
The exact full-code substitutions preserve the original
telescoping, literal phases and complete work return. Unitarity
gives the same initialized error for the actual reverse.
Positive-S normalization changes the Clifford count by at
most a constant factor.

The [geometric chunk fixtures](../../tests/test_tensor_sector_geometric_chunks.py)
check the dyadic menu weights and propagate the actual native one-hot
axis-extraction word and its reverse on every tested initialized column.
The controlled-Pauli sweeps, rolling coordinate changes and all-size
frame splice are analytic.

<a id="structured-coarse-prefix"></a>

## 6. Structured coarse schedules for larger t

The serial coarse gate program in the
[direct-prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md#5-exact-coarse-programming-with-short-selection-queries)
can be replaced at complete native boundaries. Keep the original
coarse matrices, gates within each matching, mode arenas and
conditioned prefix. Their structure allows polynomially many
R91 queries with local mode bits in the address.

### Exact prefixed operations

Retain P, X, H, d, z and w from the direct-prefix terminal.
The local mode registers use at most $`h+O(1)`$ bits, including
any original mode flags. The optional presence bit w, when used,
is preserved for the entire visit. Let $`e=H`$, or $`e=H\wedge w`$
for a presence-guarded visit.

**Permutations and signs.** A family of phase-free permutations
$`\pi_P`$ of the local mode labels has a common Beneš stride
schedule by [the exact routing construction](../endpoint/DISJOINT_BLOCK_COMPILER.md#dirty-only-routing-and-boolean-phases).
Each stage toggles one physical mode bit by a Boolean table on
P and all other mode bits. Emit its complete R91 XOR directly
into that bit. The table has zero rows on $`e=0`$ and excludes
the target. All local modes, including occupied padding, have
the specified permutation. Where signed normalization is a
monomial, retain its separate literal Boolean sign in the
prescribed order.

The original two-bit clock implements the known center phase
$`i^f`$, with f including e. Use
$`Z_4 X_4^f Z_4^\dagger X_4^{-f}`$ with its rightmost factor
first. The high-bit XOR of $`X_4^f`$ reads the unchanged low
clock bit, and the low-bit XOR excludes its own target. These
are complete exact R91 queries. The identity
$`Z_4 X_4^f Z_4^\dagger X_4^{-f}=i^f I`$ includes arbitrary
clock inputs and their references. Squaring the clock phase
gives $`(-1)^f`$. Its conjugate lift on arbitrary r is the
same $`i^f`$ times the guarded sign $`(-1)^{rf}`$, with r
unchanged. On inactive rows the increment is identity, so the
two actual clock phases cancel.

**Addressed native symbols.** For a fixed one- or two-qubit gate
G, the relevant exact table f excludes the gate's changing
operands. Toggle $`z`$ by $`ef`$, apply $`C_{e\wedge z}(G)`$
with the inherited exact bounded control, then actually reverse
the query. This acts on $`e=1,z=0`$ as the selected G and
returns z. It is identity on $`e=0`$ for arbitrary z. All fine
words leave the physical z untouched, including on their
rejected columns. For genuinely small residual words, use the
same old tuple-and-type selectors with f depending only on P.
The original data, flags and local metadata remain physical
operands of the selected gates. They are excluded from those
program addresses. A size-polynomial residual word and
polynomial common support require only polynomially many such
queries. Returned helper positions can be reused exactly as in
the original program.

For a disjoint matching of oriented two-level gates, use an
exact permutation $`\rho_P`$ mapping every ordered edge to a
pair differing on one fixed physical bit. Complete this
permutation arbitrarily on unused labels, and use identity
words on unused pairs. Keep the original exact native
single-qubit word at every pair, including its phase, and pad
their lists with identities. With the other mode bits as the
pair address, emit the bounded symbol alphabet at each
chronological position by the addressed-symbol operation.
Undo $`\rho_P`$ by its actual reverse. The selected physical
bit never enters the query address. This gives the same
logical matching on the whole mode space and restores all
helpers. Each layer needs $`O(h+1)`$ permutation queries
plus the maximum native length of its single-qubit words
times a fixed alphabet size.

**Walsh and conjugate lift.** A block Walsh on aligned dyadic
blocks is one addressed H per local bit, enabled on the blocks
in which that bit lies in the varying suffix. Its enabling
predicate excludes its target. Use an exact permutation and
reverse to handle a nonaligned block arrangement. For a
complex native coarse symbol in the conjugate lift, include
unchanged r as one further selector bit and emit the actual
conjugate symbol on its second branch. Real permutations and
Walsh gates are unchanged. A fixed finite native phase targeting
r instead keeps its inherited support and excludes r from
the programming table. Equality is of the complete lifted
operators, so internal physical T locations need not agree.

### The complete coarse ledger

The [R69 coarse remainder](HIGH_PRECISION_DEPTH.md#6-every-remaining-native-word-is-affordable-in-depth)
uses the following operators in the retained complete frame assembly.

* The [growing-source arenas](../endpoint/GROWING_SOURCE_COMPILER.md#4-full-unitary-mode-arenas-and-exact-guards)
  use polynomially many exact mode permutations for packs,
  fresh-cell injections, signed normalization and routing.
  These move original local mode coordinates and return all
  external work. The fixed native splitters and bounded
  multiplier, type, rare-label, power and guard arithmetic
  have a polynomial-size native ledger. They retain the
  original stable physical controls in the P-only symbol
  program.
* The initial source and final tag use
  [disjoint state extensions](../endpoint/DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources).
  Their coarse state preparers are block Hopf words with at
  most $`h+O(1)`$ addressed two-level matching layers.
  The inverse-polynomial-in-local-dimension accuracy of
  each coarse rotation gives a native single-qubit word
  of length polynomial in h. Retain the actual
  determinant-one words and all literal phases. Block Walsh
  and known center phases are the exact operations above.
  The chart still uses the same coarse operator and its
  actual inverse.
* The only generic adjacent-row SU mesh in this assembly
  is for the original fixed-height bottom. Its
  [blocks](../endpoint/GROWING_SOURCE_COMPILER.md#1-statement-clipping-and-the-fixed-anchor)
  have dimension $`2^{d_0}`$ for an absolute constant
  $`d_0`$. The [disjoint block mesh](../endpoint/DISJOINT_BLOCK_COMPILER.md#3-coarse-synthesis-with-a-linear-block-entry-budget)
  consequently has a constant number of matching layers.
  Its native single-qubit words have the same coarse accuracy
  and polynomial length. The initial and final state
  extensions use the block Hopf preparers above. The
  growing coarse matrices in the overlap calculation are
  realized by the original local kernel stream. Retain its
  existing fine sources, coefficient positions and guards.

Thus the constant-group assembly has at most $`P_1(n)`$ coarse
query and controlled-gate visits for a fixed polynomial $`P_1`$.
A full logical table uses only P, the $`h+O(1)`$ mode or signal
bits, and a fixed number of guard or clock bits. The small
residual gate tables read P and the fixed enabling flags.
These are at most $`C N`$ rows and one output bit per query.
Neither metadata helper banks nor precision words join a
query address. R91 uses the same exclusive $`B_Q=\lfloor b/2\rfloor`$
pool. The total new coarse T and Clifford counts are

```math
O(P_1(n)(\sqrt N+N/b+1)),\qquad O(P_1(n)(N+1)).
```

Its native T-depth is
$`O(P_1(n)(\log_2(n+2)+N/b^2))`$.
Because $`L\ge2^{cn}`$, the fixed polynomial is absorbed in
the simultaneous count and volume envelopes. The additive
non-source native depth is polynomial in n. Returned clock
bits and bounded gate-control helpers fit the inherited room
outside R91 and the live precision core. On the active program
code the replacement has the old action on all occupied modes
and dirty work. The code remains valid after every old fine
word because z is disjoint. When the macro is inactive it is
literal identity for arbitrary occupied program and mode bits.
Keep every actual coarse inverse and fine replay at its old
boundary.

For bounded local height and the retained R72 bottom, the
fallback uses a constant number of native addressed-layer
words at precision $`L+O(1)`$. These have the same exact
unpaired or paired geometric sources and complete masks.
Use the unchanged R91 masks and replace the sources by
Sections 1–4. For $`n\le3`$ the same fallback applies to the
entire input. The final R72 complex diagonal similarly has
a constant number of packed masks and geometric-source
visits at the shifted precision. Cases with bounded L also
have bounded n and can use R96.

The Section 5 source estimate is valid for $`t\le L`$, since
each retained core may use its own shorter final chunk.
Adding the polynomial coarse depth proves the extended
bound with a larger fixed $`P_{\rm c}`$. At every source
boundary the new initialized wires return on the full
original space, even if an earlier approximation has
occupied an old signal or dirty work. The coarse scheduling
never uses the new clean bank. All exact substitutions and
actual inverses therefore preserve the initialized-isometry
error, literal phases and complete dirty/reference contract.
