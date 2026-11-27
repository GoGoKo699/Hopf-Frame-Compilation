# Precision depth from short geometric source chunks

[Programmed rotations R99](PROGRAMMED_ROTATION_DEPTH.md) ·
[All-precision bound R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Direct prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md)

**Theorem R100.** There is an absolute $`c\gt 0`$ such that, for integers
$`n\ge1`$ and $`L\ge6`$ with $`n^{1/4}\le L\le 2^{cn}`$, every prescribed
complete real Hopf frame and prescribed phase-dressed complex magnitude
completion admit a native Clifford+T compiler with one external initialized
qubit, arbitrary dirty width b at the original literal thresholds
$`b\ge17(L+n+7)`$ real and $`b\ge17(L+n+8)`$ complex, and
initialized-isometry error at most $`2^{-L}`$. Put $`N=2^n`$,
$`\ell=\log_2(n+2)`$ and $`a=\log_2(L+2)`$. Simultaneously,

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL).
```

The compiler may choose either of the following native T-depth bounds at
the same count, width and error:

```math
\begin{aligned}
D_{\rm layers}&=O\bigl(NL/b^2+a\ell+(1+L/n)a^2\bigr),\\
D_{\rm groups}&=O\bigl(NL/b^2+(1+L/n)a+L^{1/8}\ell\bigr).
\end{aligned}
```

The norm includes every column, all logical and dirty superpositions,
reference entanglement, work leakage and return, and literal prescribed
signs and phases. The bounds hold for the actual native reverse.
At $`b=\Theta(\sqrt{NL})`$, the first schedule gives $`O(\log^2 n)`$
when $`L=\Theta(n)`$. The second is $`o(L)`$ for
$`n^{1/4}\le L`$ and $`\log_2 L=o(n)`$, including
$`L=2^{n^\alpha}`$ for fixed $`0\lt \alpha\lt 1`$. R96 still supplies
the all-precision bound, R99 covers low precision, and the fixed-accuracy
depth gap is unchanged.

The construction splits fixed Majorana source loaders into short exact
geometric operators and gives each a zero-clean implementation with
exponentially small operator error. The original long exact mask queries
retain their actual addresses and multiplicities. None of the new accuracy
work enumerates the long mask words.

## 1. Exact priority on arbitrary occupied keys

Suppose a native word V has code action U on a payload and a clean
k-bit key A, with exact return of that key and private scratch. Private
scratch is arbitrary. The clean-code statement is over its entire
cylinder, so V preserves both sectors of $`[A=0]`$ on the whole
input space. The case k zero requires no promotion.

Use $`h=2^k`$ physical slots, with arbitrary private payloads and
keys $`A_i`$ indexed by k-bit strings i. For
$`a_i=[A_i=0]`$, define

```math
j=\bigoplus_i i a_i,\qquad v_i=[j=i]a_i,\qquad
g=\bigoplus_i v_i .
```

Every occupation pattern has at most one asserted $`v_i`$, including
patterns with many zero keys. Exactly one zero key implies g is one.
Uniform counting of the key words therefore gives success at least
$`(1-1/h)^{h-1}\ge e^{-1}`$.

For a fixed i the residual coordinates for $`v_i`$ are the k parities
$`j_b\oplus i_b`$ and the one bit $`a_i\oplus1`$.
Each is an exact returned XOR into an arbitrary target. For a parity,
apply the atomic read-only conjunction phases on O(k) key literals
in private banks, with the Fourier character of the XOR target read
only. Their T bands run in parallel. Apply the
[two displaced one-hot banks](PROGRAMMED_ROTATION_DEPTH.md#exact-phases-for-products-of-hash-tests)
to the residual XOR vector. This implements the literal sign of its
zero test, with the predicate inputs unchanged. Extra stable literal
controls cost one residual coordinate each. A Hadamard around the
phase target turns the sign into a literal XOR of $`v_i`$.
For g use all such phases in parallel, again with private T targets
and a shared read-only Fourier character. These complete XOR words
have depth $`O(\log(k+2))`$ and count and width $`2^{O(k+1)}`$.
Returned helpers for all these words are kept disjoint from the source
slots' own scratch.

Combine $`c_k=\lceil C(k+1)\rceil`$ disjoint subgroups of these
slots into a base group of w key bits, where
$`w=c_k h k`$. The fraction of configurations where all subgroups
fail is some $`\rho\le(1-e^{-1})^{c_k}`$.
For a power of two $`D\ge C'(w+1)/\rho`$, choose D translation
masks for the w key bits. Uniform independent masks make the number
of failures for each fixed key word binomial with mean $`D\rho`$.
The binomial bound $`\Pr[X\ge2D\rho]\le e^{-D\rho/3}`$ and a
union bound over the $`2^w`$ words give a list with at most
$`2D\rho`$ failing labels for every word. For example the bound
follows by applying Markov to $`2^X`$ and using
$`1+\rho\le e^\rho`$. Enumerate lists and check the integer
counts to choose one deterministically. Use the least eligible power
of two. Here $`\log_2 D=O(k+1)`$, since the all-fail fraction is
also at least the probability of no zero key in any subgroup, which
is exponential in $`c_k`$ with absolute constants. The exact
controlled translation uses one row-indicator phase on the dirty seed
label per row and a Clifford parity for its changed key coordinates.
Supply private lanes for all rows. Its count and width are
$`2^{O(k+1)}`$, and its native depth is $`O(\log(k+2))`$.

In basis 1 key and seed labels are computational. In basis 2 conjugate
by Hadamards on the seed alone. For each fixed key word both bad sets
have size at most $`2D\rho`$. Bounding the corresponding restricted
Walsh matrix by its Frobenius norm gives an overlap at most
$`\mu=2\rho\sqrt D`$. Choosing C large makes $`\mu\le1/4`$
uniformly for k at least one, since $`D=O((w+1)/\rho)`$
and $`w=O((k+1)^2 2^k)`$. Taking the maximum over key words
gives the whole-space operator bound.

For precision parameter s take $`R=O(s+1)`$ independent physical
base groups. In each seed basis select the first successful subgroup
across all groups, then its unique valid slot. The total failure projector
in basis t, denoted $`P_t`$, is the tensor product of the R base
failure projectors in that basis. Hence $`\|P_1P_2\|\le\mu^R`$
on every input. For each slot i in subgroup u the complete selection
predicate is

```math
b_{u,i}=v_{u,i}\prod_{u'\lt u}(1-g_{u'}).
```

Use the [returned dirty sum and conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
on $`O((k+1)(s+1))`$ Boolean inputs to implement each selected
XOR, including any further stable literal controls.
In the initial parallel controlled increments, use the
[two-addition construction](../../docs/PARALLEL_DIRTY_LOOKUP.md#a-read-only-controlled-increment-using-two-additions).
Its chronological six calls are
$`C_g,Q^\dagger,E,Q,E^\dagger,C_g`$ in the notation of that
construction. Replace the literal toggle E into the private low
control bit and its reverse by the complete read-only XOR for the
desired input and its actual reverse. These are the only visits to
that input in the increment. The arithmetic and cyclic test use that
complete increment and its inverse. Distinct instances have private
T targets. All resulting predicates are exact at depth
$`O(\log(k+2)+\log(s+2))`$ and with count and width
$`2^{O(k+1)}\mathrm{poly}(s+1)`$, including the disjoint
arithmetic and source control scratch for all copies.

## 2. A selected center and the two seed bases

Suppose $`B,\bar B`$ are Hermitian payload Paulis and
$`M=U^\dagger B U`$ anticommutes with $`\bar B`$.
In one seed basis swap the primary payload into the selected slot
and use identity if no slot is selected. Implement the swap by the
three XOR sweeps: toggle each slot bit by the corresponding primary
bit and its select predicate, toggle each primary bit by the parity
of selected slot bits, then repeat the first sweep. Each is an exact
returned Fourier phase batch with private T targets. Counts and width
grow polynomially with payload width, and depth is unchanged. Use
the actual inverse for the reverse swap.

Apply V to all slots, apply the center $`e^{i\phi B}`$
ONLY in the selected slot, and reverse all copies of V. The center
uses independent returned scratch and its flags depend only on the
preserved zero/nonzero indicators and the read-only seed. On every
joint indicator sector, all unselected V calls cancel exactly. On
the selected sector, its key is zero and the code action gives
$`e^{i\phi M}`$. Reverse the payload swap and the basis/shift
operations. The result acts exactly as

```math
A_t(\phi)=P_t\otimes I+(I-P_t)\otimes e^{i\phi M}.
```

All auxiliary payloads and scratch return. This identity includes
coherent superpositions of occupied slots.

The only center angles needed are $`\phi=\pi/2,\pi/4`$.
For literal implementation put $`\omega=e^{i\pi/4}`$.
For any of the exact stable predicates b, the chronological word
$`T_r^{-m}, E_b,T_r^m,E_b^\dagger`$ on an arbitrary bit r has
phase $`\omega^{mb-2mrb}`$. Correct the second term by the same
word at coefficient $`2m`$ for the predicate rb on a separate
borrowed bit, then at coefficient $`4m`$ for the newly introduced
residual predicate with a further private bit. The final residual
is a multiple of eight. Thus the literal phase $`\omega^{mb}`$
is paid by constantly many exact predicate visits. Diagonalize B to
Z on its slot by a fixed Clifford. A guarded $`e^{i\phi Z}`$
factors into the phases of b and the conjunction of b with that data
bit, including the known scalar on the positive sector. It uses the
same private parallel schedule. Use the parity of the exclusive
selection flags to implement the mirror
$`C_t=P_t\otimes I+(I-P_t)\otimes\bar B`$ on the primary
payload by the same literal sign method.

The zero-clean operator is

```math
W=A_1(\pi/2)A_2(\pi/4)C_1 A_2(\pi/4)C_1^\dagger .
```

For orthogonal bad projectors, their good indicators commute and at
least one is true. On the first good sector the mirror reverses the
angle of the second rotation, leaving the first angle of $`\pi/2`$.
Otherwise both copies of the second rotation act, with no mirror.
Thus $`W=iM\otimes I`$. For the actual projectors, the projector
onto $`(I-P_1)\mathrm{ran}P_2`$ has distance at most
$`\mu^R`$ from $`P_2`$ by principal angles. Replacing it
in the two $`A_2`$ factors changes the whole word by at most
$`4\mu^R`$. Choose R so this is below $`2^{-s}`$.
The construction returns the full arbitrary work space to this error,
with no supplied seed or success bit. Its depth is
$`O(D_T(V)+\log(k+2)+\log(s+2))`$, and count and width are
$`2^{O(k+1)}`$ times a fixed polynomial in $`s+1`$, payload
width and V's native resources. The same error and resource orders
hold for the physical reverse. All Clifford gates used to implement
the predicates, read-only interactions and swaps are charged.

## 3. Exact geometric chunks and the powered source loader

For any independent anticommuting Hermitian Pauli list
$`\Gamma_0,\ldots,\Gamma_{m-1}`$, $`m\ge2`$, define

```math
G_m=\sum_{j\lt m-1}2^{-(j+1)/2}\Gamma_j+2^{-(m-1)/2}\Gamma_{m-1}.
```

It is an involution by the geometric sum of squared coefficients. The
coefficients are rational dyadics or dyadics divided by $`\sqrt2`$.
The dyadic part sums to at most one. For the second part write the
numerators as $`z_j`$, with $`z_0=1`$ and
$`\sum z_j\le2`$. These inequalities are the two geometric
sums of alternating indices, counting the last coefficient at its
actual repeated exponent. For every $`j\ne0`$ in the second
part, use weight $`z_j`$ on the Hermitian Clifford
$`(\Gamma_0+\Gamma_j)/\sqrt2`$. The remaining nonnegative
$`z_0-\sum_{j\ne0}z_j`$ is represented by equal weights on the
plus and minus such Clifford for any other axis. Use the original
weights on the Pauli atoms of the dyadic part. The total atom weight
is at most two, with denominator exponent O(m). Divide by two
and pad to one with cancelling $`+I,-I`$ atoms using one extra bit.
A uniform O(m)-bit initialized binary selector, fixed intervals of
these dyadic weights and literal SELECT give compression $`G_m/2`$.

For a paid SELECT, allocate a private dirty payload slot per full
selector label, swap into the unique slot and reverse by the exact
three sweeps. Pauli atoms have an exact central Pauli controlled by
the read-only label equality. Pair atoms are a Pauli conjugated by one
Pauli-axis T rotation and its actual reverse, with only the center
guarded. For $`-I`$ use the literal equality phase. All label
tests use exact private conjunction phases. Thus even this direct
allocation has $`2^{O(m)}`$ T, Clifford and work cost on an
O(m)-line core, with $`O(\log(m+2))`$ T-depth. The
[exact half-amplitude amplification](PROGRAMMED_ROTATION_DEPTH.md#an-exact-invariant-clean-source-with-a-short-program)
has code action $`-G_m`$, with all arbitrary scratch returned.
Since that compression is unitary, the amplified word preserves the
whole initialized selector sector exactly.

Take a chunk of t consecutive splits on the chain axes
$`\Gamma_j,\ldots,\Gamma_{j+t}`$. The correct normalized
output for input axis $`\Gamma_j`$ is the geometric involution
$`M_{\rm chunk}`$ on that list. Add two spare axes
$`\Gamma_a,\Gamma_{aa}`$ anticommuting with all original axes.
They exist on one additional dirty qubit after extending the old
family to an independent even list. Explicitly, for an even list of
length m use its Hermitian chirality Pauli
$`\Omega=i^{m/2}\Gamma_0\cdots\Gamma_{m-1}`$ and new
axes $`\Omega\otimes X,\Omega\otimes Y`$. Both retained
source families use the standard Majorana lists, so any missing
last partner in the old family can be included before taking chirality.

Put $`G'=(\Gamma_a+M_{\rm chunk})/\sqrt2`$ and
$`G''=(\Gamma_{aa}+G')/\sqrt2`$. Both are geometric lists
with one additional leading split at each step. Use the exact amplified
word for $`-G''`$ as V in Section 2, with
$`B=\bar B=\Gamma_{aa}`$. Then
$`G''\Gamma_{aa}G''=G'`$ and $`G'`$ anticommutes
with the mirror. This implements $`iG'`$ on the whole space
with error $`2^{-s}`$.

For the selected slot circuits the O(t) participating axes can be
mapped to O(t) physical core lines by a fixed Clifford. Their
symplectic pairings agree with the standard list. To extend the map,
complete a symplectic basis of the span, pair its radical basis with
dual vectors in the ambient space and complete a symplectic basis
of the remaining complement. A Pauli corrects the signs of the
independent images. Ordinary Clifford elimination realizes the
resulting frame on an old core of m lines with $`O(m^2)`$
elementary Cliffords (eliminate a canonical X,Z pair at a time
by O(m) gates). The exact inverse of this frame map returns all
unused core coordinates. This charges long Pauli supports by a fixed
polynomial in the old core size.

Choose a Clifford plane rotation C carrying $`\Gamma_j`$
to $`\Gamma_a`$ and fixing every outside original axis.
The unitary $`G'\Gamma_a C`$ sends $`\Gamma_j`$ to
$`M_{\rm chunk}`$ by conjugation. Every outside original
axis is negated by each of the two factors $`\Gamma_a,G'`$
and hence fixed. This implements the exact required transformation of
the source axes up to the counted operator error, regardless of its
extension outside their span. The known scalar i cancels in loader
conjugations. Across successive chunks, the coefficient of the
last split axis multiplies the next normalized list. Keep the
original two seed splits for a paired source, then perform the
chunk replacements along both tails. The outside-axis identity
also fixes the other tail. For an original precision-size O(m)
loader, the comparison unitary $`U_*`$ has
$`U_* X_0 U_*^\dagger=M_{\rm old}\otimes I`$ on the
enlarged arbitrary core. Absorb the known chunk scalars into
$`U_*`$, which leaves that conjugation invariant. Its actual
approximation $`\widetilde U`$ uses $`O(1+m/t)`$ short chunks, with depth
$`O((1+m/t)[\log(t+2)+\log(s+2)])`$ and per-chunk
count and width $`2^{O(t+1)}\mathrm{poly}(m,s+1)`$.

Replace each completed controlled source by

```math
U_{\rm old} C_F(X_0) U_{\rm old}^\dagger
\quad\longmapsto\quad
\widetilde U C_F(X_0)\widetilde U^\dagger .
```

The center and its predicate keep their original exact implementation
and returned helpers, with fresh banks for the loader. The full
operator difference is at most twice the forward-loader error.
On F false, the actual new word is literal identity for every occupied
work state. An additional external guard on the old complete source
can likewise be applied solely in the center, since that guard is
disjoint from the precision core and the new auxiliary space.
For an uncontrolled source simply drop the center control.

The original conjugate-signal form can be replaced at the same complete
boundary. The unpaired axes are real. On the paired axes
$`\Gamma_{2j}=Z^{\lt j}X_j`$, $`\Gamma_{2j+1}=Z^{\lt j}Y_j`$,
the Pauli word $`Q=X^{\text{all}}Z^{(j\bmod 2)_j}`$ negates
exactly the imaginary axes. Conjugate the new completed word by
$`C_r(Q)`$ and its actual inverse as prescribed by the old lift.
It is an ordered product of CNOTs and CZs, so adds no T-depth even
with a shared r. It gives the intended real and conjugate source in
the two r sectors. It leaves the exact X/Z mask lists, their order
and their actual inverses intact.

## 4. Parallel capped geometric reference

Let M be a power of two with $`M\ge2`$. The normalized
uncapped column $`\psi`$ of amplitudes proportional to
$`2^{-j/2}`$, $`0\le j\lt M`$, is a product on its binary
digits. Its real prep uses one prescribed $`R_y`$ rotation per
digit in parallel. For the last atom b, put
$`p=|\langle b|\psi\rangle|^2=1/(2^M-1)`$.
The capped last probability is $`q=2p/(1+p)=2^{-(M-1)}`$.
Write $`\psi=\sqrt{1-p}g+\sqrt p\,b`$ with g supported
off b. Phase b by $`e^{i\phi}`$ and apply
$`I-2|\psi\rangle\!\langle\psi|`$. The last coefficient
becomes

```math
\sqrt p\,[(1-2p)e^{i\phi}-2(1-p)].
```

Its magnitude ranges from $`\sqrt p`$ to
$`(3-4p)\sqrt p\gt \sqrt q`$ since $`p\le1/3`$.
Choose the angle by the certified cosine formula, then correct the
relative phase by one last-atom phase and the common phase by one
unconditional phase. The component off b is a single nonzero scalar
multiple of g. These operations give the exact positive capped
column on the ideal initialized input.

For the phase corrections, implement the last-atom phase as
$`e^{i\phi[\text{coin}=b] Z_d}`$ on an additional arbitrary
dirty port d. Do the same for the common phase with predicate one.
On d's negative eigenspace the complete operation is complex
conjugated, so both sectors produce the same real capped column and
return d. This gives its precise initialized-cylinder projector with
identity on d and all new dirty work. The ideal prepare word may
depend on d on other columns. Use the
[borrowed-signal rotation](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
and a returned controlled Pauli center in a half-angle echo. The
exact comparison center of the reflection about $`\psi`$
is a zero-coin test. Use the actual product prep and its actual
reverse there. The phase guards are exact conjunctions on the coin
digits. Section 3 accelerates the fixed source loaders of the
borrowed-signal rotations. The parallel digit banks and the
constant number of phase/reflection stages have independent
private work. All comparisons are full-operator hybrids up to the
exact ideal cylinder map above. M one is the identity.

## 5. Opening and shared physical allocation

Build the [R96 moderate-precision comparison](ALL_PRECISION_LINEAR_DEPTH.md#4-moderate-precision-three-early-stages-and-a-short-terminal)
with its original finite magnitude and sign reduction and order of
factors. In an active stage of height m and precision l equal to L
plus an absolute constant, retain

```math
K=m+l+2,qquad
\chi=\lceil\log_2(l+2)\rceil+\lceil\log_2(m+2)\rceil^3,
\qquad r=\lceil A\chi\rceil .
```

The original global prefix has length n minus m. The stage uses its
saved-activity word and retains the original bottom r coordinates.
If the reserve is too large, skip the stage as before. No new fields
enter the initialized-work anchors or the long mask addresses.

Replace the three serial precision contributions in the
[original ledger](GLOBAL_BORROWED_PRECISION_STAGE.md#6-native-count-width-and-depth-ledger).

* **Calibration.** Keep the borrowed-signal real rotation and its
  original preserved source label and global prefix. Replace each
  complete fixed paired source by Section 3. Its O(l)-bit masks
  enter the original arbitrary core through only constantly many
  complete X and Z masks and their actual reverses. Retain the exact
  central Pauli guards of the existing borrowed-signal calibration.
  On the inactive sector the complete source and mask pairs cancel.
  Give any constant accuracy padding to the original table outputs.
* **Fixed coefficient rotations.** Use the same borrowed-signal
  rotation with fixed Clifford Pauli masks, the old angle data, and
  the accelerated fixed loaders. Pay the existing half-angle guard
  with a returned Pauli center. Use disjoint physical instances for
  the simultaneous coefficient sites. Allocate the old bank error
  divided by its number of sites, which needs
  $`l+O(\log(L+2))`$ local precision. Conditional work remains
  arbitrary off its initialized sector.
* **Geometric references.** Use Section 4 on each bin's original
  geometric size, with independent dirty work and no additional
  initialized field. Use a constant number of parallel angle banks
  and exact zero or last-atom tests on O(log(L+2)) literals.
  The comparison preparer needs the same initialized cylinder
  for the full-space transported tests, by the
  [reference projector identity](FAST_CAPPED_GEOMETRIC_REFERENCE.md#2-a-full-space-comparison-for-transported-zero-tests).
  Old initialized coin scratch can remain untouched. Its tests
  exclude all new dirty work. On the inactive sector the new actual
  preparations and their actual reverses cancel as before.

The first two replacements compare in full operator norm to the old
completed words, or to the same ideal borrowed-port rotation with its
original error allowance. For the last use the exact ideal cylinder
of Section 4 and a full-unitary hybrid of its own calls. The old
[source and complete-stage error estimates](GLOBAL_BORROWED_PRECISION_STAGE.md#7-complete-error-and-return)
then use the same ideal anchors. Only formerly initialized fields
enter those anchors. All leakage propagates through the full circuit.

Work first on $`b\ge\sqrt N/n`$ above one absolute cutoff. Set the maximum
chunk length $`t=\max(1,\lfloor\delta n\rfloor)`$ with fixed small
$`\delta\gt 0`$. The total number of chunk occurrences in either schedule
is bounded by $`\mathrm{poly}(n,L)2^{O(\ell^3)}`$. This counts
the fixed polynomial number of sources and scratch fields in the
opening layout, polynomially many fine macro positions in the group
height, and at most O(L) chunks per source. It also includes the
coefficient and reference banks, individual layers, and the final
phase. The old full-volume mask queries are unchanged and introduce no
substitution error. Use precision

```math
s=L+C_1(\log_2(n+L+2)+\ell^3)+C_2
```

with fixed sufficiently large constants. Then s is O(L) and the sum
of the full-operator loader hybrid errors is below $`2^{-L}/64`$,
including the actual reverses. For the individual angles in parallel
banks divide the original error share by their number. All those
angle sources still have widths O(L).

The old standalone bound is a fixed power of K times $`2^{O(\chi)}`$.
Additional short query pools, fixed transports, exact guards, core
work and all new primitive costs together are bounded by

```math
2^{C_3t}\mathrm{poly}(n,L)\,2^{C_4(\log_2(L+2)+\ell^3)} .
```

This includes all independent tensor groups, every exact priority
copy, every Clifford parity and frame change, every actual inverse,
and all the simultaneous instances within each chronological bank.
Fix the circuit exponents, constants and the opening reserve A first,
then choose $`\delta`$, then c. Choose $`C_3\delta\lt 1/16`$ and
c so the fixed coefficient of cn is less than 1/16, with
$`2Ac\lt 1/8`$. Allow the internal range $`L\le2^{2cn}`$ for the complex
shift in these estimates. Since $`\ell^3=o(n)`$, the standalone bill
and nonquery width are at most $`2^{n/4}=o(\sqrt N/n)`$ above one
fixed cutoff. The original long query pools keep a fixed fraction of
b, with physical T targets disjoint from new work.

Every original long opening table still has address length at most
$`n-r+C\chi`$, with the same fixed polynomial multiplicities. The
original reserve suppression applies after fixed coefficient enlargements.
Each shared address is read by Clifford interactions in independent
banks with all non-Clifford targets private. With at most three
openings and the two exact sign phases, the aggregate ledger is

```math
T_{\rm open}=O(\sqrt N+N/b),\quad G_{\rm open}=O(N),\quad
D_{\rm open}=O\bigl(N/b^2+\ell+(1+L/n)a\bigr).
```

Other opening depths remain logarithmic in the corresponding
addresses or short fields. An accelerated O(L)-bit source has
O(1+L/n) chunks of depth O(a). The inherited contractions leave
M=O(a) original bottom layers. Here $`a=\Omega(\ell)`$,
$`m_1=O(a+\ell^3)`$, and a further contraction, even after a
failed reserve test, bounds the remaining height by
$`O(a+\log_2^3(a+2))=O(a)`$.

## 6. Individual original tail layers

Use the [R72 addressed real rotation](../../docs/ONE_CLEAN_COMPILER.md#1-contract-and-addressed-primitive)
at every remaining layer. Index a layer by j, with j minus one
lower bits and n minus j upper address bits. Its full O(q)-bit masks
have $`Q_j=N2^{-j}`$ rows, and the lower bits enter only a separate
exact central conjunction. Use $`q_j=L+j+C`$ with an absolute C to
sum the rotation errors with a geometric allowance. The conjunction
has O(log(j+2)) depth by the read-only exact conjunction of Section 1
and its own returned helpers.

Implement each complete packed X or Z mask using the
[sharp R91 whole-word query](PARALLEL_BILINEAR_WORD_QUERY.md#4-an-allocation-at-every-sufficient-dirty-width),
with an exclusive pool of order b. This exact query removes the additive
q depth of the older layerwise schedule and preserves the ordered Pauli
word and its inverse. The suffix test stays outside the full prefix.
Use the borrowed flag or the external flag with its corresponding
full-operator or code contract. Section 3 replaces each fixed source
and the actual inverse appearances at the specified hybrid error.

The mask products and square-root products sum geometrically,

```math
\sum_j Q_jq_j=O(NL),\qquad
\sum_j\sqrt{Q_jq_j}=O(\sqrt{NL}) .
```

All standalone terms fit Section 5. The other depth on each layer is
$`O(\ell+(1+L/n)a)`$, since j is O(a). Summing with the original
width term gives the first depth bound. The term $`a\ell`$ is
absorbed into $`(1+L/n)a^2`$ on this precision range.

## 7. Groups with the original prefix grammar

Alternatively use the [direct-prefix terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md#2-a-common-local-high-precision-grammar)
on groups of height at most $`g=\max(1,\lfloor a/512\rfloor)`$
and retain up to three individual bottom layers. The number of nonempty
groups is bounded by an absolute constant. Each group uses its actual
original prefix and full-suffix flag with exact saved-predicate inverse.

First construct the original comparison circuit at $`\kappa=L+O(1)`$.
Fix its coarse instruction data, fine-mask rows, mixed gauges and row
and column words as in the direct-prefix proof. Replace only completed
fixed M or controlled-M source macros. The fine-query lists stay
unchanged. When a modified word is replayed in reverse, use its actual
reversed native word. The full-operator hybrid compares its rejected
action to the whole old word, so no gauge needs to be recomputed.

The terminal's [precision-sized elementary words](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md#4-elementary-length-of-precision-sized-native-words)
consist of the Majorana loaders, output Hadamard banks and strictly
halving [dyadic refinements](../endpoint/GROWING_SOURCE_COMPILER.md#3-exact-source-refinement-and-carried-kernels).
The long scalar in a refinement uses the raw Majorana source and a
fixed dyadic mask. Its rare prep, interleaving, controlled X,
multiplier and private mode fields use the original local polynomial
gate count. The completed source positions have total precision width
$`O(L+P(h))`$ and count at most $`P(h)`$ for a fixed polynomial
P and actual height h. Give each source the Section 3 replacement.
Put the external activity and optional bin-presence controls solely
into the exact Pauli center. Where the old grammar lifts a source
to the conjugate signal, use the completed source lift of Section 3.
A source already inside a complete conjugate-signal gadget needs
only its ordinary replacement. This preserves literal inactive
identity for occupied d, z and w. Full masks and their actual
inverses retain the original saved-predicate boundary.

For the original suffix size s, write $`N_s=N2^{-s}`$.
The original fine and coarse queries retain the sums of row products
$`O(N_sL)`$ and square-root products $`O(\sqrt{N_sL})`$, with
separate saved-predicate costs as in R96. The fine positions number at
most P(h) and the coarse positions at most $`2^{3h}P(h)`$, both
$`O(L^{1/8})`$ by the group size and fixed polynomial.
Their total depth is $`O(N_sL/b^2+L^{1/8}\ell)`$.
The fixed sources contribute

```math
O\bigl((1+L/n+P(h))a\bigr),
```

since the sum of chunk counts is $`O(1+L/n+P(h))`$.
The $`P(h)a`$ term is absorbed by $`L^{1/8}\ell`$.
Constantly many groups, their saved flags and up to three bottom
layers at their original addresses give the second native depth bound.
All original prefix-selected masks use the same exclusive returned
query pool. The simultaneous standalone work fits the remaining
allocation from Section 5.

## 8. Error, literal complex phase and inverse

For the real frame, retain the sign and stage errors of R96.
The arbitrary dirty reference helpers are included in its ideal
augmented projector. In the layer schedule use the geometric error
sum from Section 6. The grouped comparison keeps the original group
and bottom error allocation of R96. The substitutions are full-operator
hybrids with error at most the fixed fraction in Section 5. Padding
by absolute constants as needed gives real error below $`2^{-L}/2`$.
The ideal prefix keeps the initialized cylinder and the actual suffix
acts unitarily on all leaked components.

For the complex magnitude completion, compile the real part at L+1,
then use the [one-clean literal diagonal](../../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
at error $`2^{-(L+1)}`$, including the substitution error.
This uses the same external flag with its original initialized
phase contract. Its constant number of packed O(L)-bit masks have N
rows and use R91 with their original full addresses. Replace every
fixed source and actual reverse by Section 3. The added depth is
$`O(NL/b^2+\ell+(1+L/n)a)`$, with the same count bounds and
reused pool. The literal threshold for the real stage at L+1 is the
stated complex threshold. This handles phases on zero-amplitude leaves.

For $`b\lt \sqrt N/n`$, use R96 at its original precision and
minimum width, since $`NL/b^2\gt n^2L`$. Bounded dimensions in the
admitted precision range also use R96 at an absolute depth constant.
Every branch uses one external clean bit and only the priced dirty
work. For the ideal $`U=W\otimes I_b`$, multiplying by native
$`V^\dagger`$ on the left and $`U^\dagger`$ on the right gives
$`\|V^\dagger J_1-J_1U^\dagger\|=\|VJ_1-J_1U\|`$.
The actual inverse has the same T-depth and T count, and at most a
constant factor in Clifford count for the positive-S alphabet.

At eligible widths where the additive term of either bound is at most
$`NL/b^2`$, the [one-clean T-count lower bound](../../docs/OPEN_PROBLEM.md#optimal-banked-count)
divided by the O(b) physical targets per native layer gives matching
worst-case depth. The new accuracy work acts on fixed sources.
A sublogarithmic fixed-accuracy source for an arbitrary prescribed
Boolean table remains open.

The [five bounded fixtures](../../tests/test_tensor_sector_geometric_chunks.py)
check exact dyadic menu weights, complete chunk transport, occupied
priority selections, tensor bad overlap and both dirty-port sectors
of the capped reference. The all-size native schedule and frame splice
use the analytic arguments above.
