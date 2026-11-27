# A clean-space tradeoff for large precision

[All-precision bound R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Exact key-cover bound R104](EXACT_KEY_COVER_PRECISION.md) ·
[Precision count–depth lower bound R103](PRECISION_COUNT_DEPTH_TRADEOFF.md)

**Theorem R105.** Fix a sufficiently small absolute $`c\gt0`$.
For integers $`n\ge1`$, $`L\ge6`$ and $`L\ge2^{cn}`$, put
$`N=2^n`$ and choose an integer $`1\le t\le L^{3/4}`$.
Every prescribed real complete Hopf frame and prescribed complex magnitude
completion has a native Clifford+T word and its actual reverse with
initialized-isometry error $`2^{-L}`$, at most $`1+C_{\rm cl}(t+1)`$
external clean inputs, and the literal R96 dirty thresholds
$`b\ge17(L+n+7)`$ and $`b\ge17(L+n+8)`$, respectively.
The simultaneous resources are

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L\log_2(t+2)/t).
```

All additional clean inputs return exactly. The norm includes the
original clean output, all logical columns, arbitrary dirty inputs and
references. For $`L\ge N`$, R103 gives the lower bound
$`\Omega(L/(k+1))`$ on worst-case depth at $`O(L)`$ T gates with
$`k`$ external clean inputs. Taking $`k=1+C_{\rm cl}(t+1)`$ puts the
upper bound within a factor $`O(\log_2(k+2))`$. Here $`NL/b^2=O(1)`$.
At $`t=1`$ R96 gives the bound with one clean input.

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

Since $`t\le L^{3/4}`$, the $`L/t`$ term absorbs every fixed
polynomial in n and the coarse depth. The real part and final
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
