# Grouped programmed precision for complete Hopf frames

[Programmed rotations R99](PROGRAMMED_ROTATION_DEPTH.md) ·
[Parallel word queries R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Capped geometric references](FAST_CAPPED_GEOMETRIC_REFERENCE.md) ·
[All-precision baseline R96](ALL_PRECISION_LINEAR_DEPTH.md)

Group neighboring bottom layers into geometrically increasing blocks. A
short coarse program gives an exact dyadic unitary on each block. A
near-identity dyadic graph reflection carries the fine correction. The
selector examines correction coordinates at only one invariant candidate
per group, so the approximation error is paid for polynomially many
slots even though the construction has exponentially many dirty slots.

**Theorem R101.** There is an absolute $`\xi\gt 0`$ such that, for
$`N=2^n`$ and integers $`\max(6,\sqrt n)\le L\le\xi n`$, every
prescribed complete real Hopf frame admits a native word V with one
external initialized qubit and any dirty allocation $`b\ge17(L+n+7)`$.
For the prescribed phase-dressed complex magnitude completion the
literal minimum is $`b\ge17(L+n+8)`$. Write $`\ell=\log_2(n+2)`$
and $`\lambda=1+\log_2\log_2(L+4)`$. The same word obeys

```math
\begin{aligned}
\|VJ_1-J_1(W\otimes I_b)\|_{\mathrm{op}}&\le2^{-L},\\
T&=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\\
D_T&=O(NL/b^2+\ell\lambda).
\end{aligned}
```

All constants are absolute. The norm covers every logical column,
literal phase, arbitrary dirty input and reference, and all work return.
The bounds also hold for the actual native reverse. At
$`L=\Theta(n)`$ within this window and $`b=\Theta(\sqrt{NL})`$, the
depth is $`O(\log n\log\log n)`$, improving R99's $`O(\log^2 n)`$.
[Section 8](#8-small-clean-menus-through-quadratic-precision) permits
$`\sqrt n\le L\le n^2`$ with additive depth
$`O((1+L/n)\ell\log_2(\ell+2))`$, where $`\ell=\log_2(n+2)`$.
The fixed-accuracy frontier is unchanged. Classical preprocessing
terminates, with no running-time claim.

## 1. Block geometry and two borrowed zero bits

Keep R99's opening and its finite sign reduction, with the parallel cap
of Section 6. After at most three contractions it leaves
$`M=O(\log(L+2))`$ original bottom levels. In the
[complete addressed order](../../docs/HOPF_INTERFACE.md#4-complete-addressed-depth-operator),
a consecutive block of h levels with j original lower bits has
preserved prefix x of length $`n-j-h`$. On the sector where all j lower
bits are zero it applies the complete h-qubit real block $`U_x`$;
on every other suffix it is identity. Internal suffix controls remain
part of $`U_x`$.

Use R99's individual rotations for a fixed number of deepest levels.
For the other blocks, compute the exact suffix-zero predicate into the
external clean bit. On its active sector use two of the j lower bits
as $`c_1,c_2`$, both initialized conditionally. The physical block will
be literal identity on an inactive flag, including arbitrary occupied
$`c_1,c_2`$. Apply the actual suffix-predicate inverse last.

Design the block boundaries from the bottom, then execute them in
the original addressed order. Above a fixed bottom cutoff set
$`h=\min\{\lfloor j/C\rfloor,\lfloor c'\log_2(L+2)\rfloor,M-j\}`$, with fixed large C and fixed $`0\lt c'\lt 1/2`$.
Increase the dimension cutoff to ensure the cap is positive. Until
the cap is reached the bottom distance grows geometrically. The cap is a fixed fraction of $`M=O(\log(L+2))`$, so only a bounded
number of capped blocks are needed. There are
$`O(1+\log\log(L+4))`$ blocks. Every block below uses
$`2^{O(h)}(L+O(h+\log(n+2)))`$ program bits per addressed row.
Choose C after the exponent so its addressed volume is
$`O(NL2^{-j/2})`$.

## 2. Exact coarse block

Apply [R99's five-axis construction](PROGRAMMED_ROTATION_DEPTH.md#1-an-exact-invariant-clean-source-with-a-short-program)
to each whole addressed rotation at accuracy $`\epsilon_0/h`$, with
$`\epsilon_0=2^{-C_1(h+1)}`$ and $`C_1`$ a fixed large constant. On an
extra arbitrary rebit, realify the two imaginary axes $`iX_0X_1`$ and
$`iX_0Z_1`$ by replacing their factor i with $`iY`$ on the rebit.
The five resulting matrices are real, pairwise anticommuting Pauli
sections squaring to $`-I`$. The target real rotation tensored with
identity is unchanged. The product $`U_{0,x}`$ is therefore exactly real
dyadic orthogonal on the h logical bits and three arbitrary helpers,
and $`\|U_{0,x}-U_x\otimes I\|\le\epsilon_0`$.

The clipped cutoffs for the $`O(2^h)`$ local tree nodes occupy
$`O(2^h(h+1))=O(L)`$ bits in a slot program P. Broadcast the program
into arbitrary P fields from a central dirty R91 query output by
the four-call Clifford-copy/query/Clifford-copy/actual-query-inverse
echo of R99. Each slot additionally has an $`O(h+1)`$-bit local
selector and cutoff scratch A, and a direction bit z.

In each addressed layer, an exact short XOR into A loads the cutoff
of the current local prefix. Require the layer's lower suffix and z
to be zero to use that entry. Otherwise load the valid identity
code with $`a=2^K`$ and all other coefficients zero in R99's menu.
For fixed controls this lookup is affine in P, and its individual
read-only conjunction phases have private T targets. Apply the
R99 clean amplified rotation with its literal accepted minus sign
corrected by a Clifford $`-I`$. Unload the cutoff by its actual
inverse from the same preserved address and z. Finally reverse the
program broadcast. At original $`A=P=0`$ each layer returns its
local A; the completed slot is exactly $`U_{0,x}`$ at z zero and
identity at z one, preserving z. The projector $`[A=P=0]`$ is thus
invariant on the full slot space by unitarity. Slots factor at fixed
x. The local depth is $`O(h\log(h+2))`$, and the global program
uses constantly many visits.

## 3. An exact dyadic graph for the fine correction

On dimension $`D=2^{h+3}`$ put $`V_x=(U_x\otimes I)U_{0,x}^\top`$
and $`R_x=V_x^{1/4}`$, taking the principal real orthogonal root.
The binomial series at I gives certified classical approximations
since $`\|V_x-I\|\le\epsilon_0`$. Set
$`s=L+O(h+\log(L+2)+\log(n+2))`$ with fixed slack specified below.

Choose $`\nu=2^{-2s-C_2}`$. Round $`(1-\nu)R_x`$ entrywise
to a real dyadic A with denominator $`2^K`$, where
$`K=2s+O(\log(D+2))`$. Take K large enough that the row-sum
norm of $`A^\top A-(1-\nu)^2I`$ is less than $`\nu/2`$.
For example the entrywise rounding error $`O(2^{-K})`$
gives spectral error $`O(D2^{-K})`$ and row-sum error
$`O(D^2 2^{-K})`$ in this Gram matrix. Then
$`M=I-A^\top A`$ is diagonally dominant with positive diagonal
and $`\|M\|=O(\nu)`$.

For each off-diagonal entry $`M_{ij}`$, use Lagrange's four-square
theorem on $`2^{2K}|M_{ij}|`$. Four dyadic rows proportional to
$`e_i+\mathrm{sgn}(M_{ij})e_j`$ realize $`M_{ij}`$ and add
$`|M_{ij}|`$ to the corresponding diagonals. Cover the remaining
nonnegative diagonals by four squares each. The resulting B has
$`O(D^2)`$ rows, denominator $`2^K`$, $`B^\top B=M`$ and
$`\|B\|=O(2^{-s})`$. Pad vertically to dimension $`D'=O(D^2)`$,
a power of two. The matrix $`\mathcal K=[I;A;B;0]`$ has
$`\mathcal K^\top\mathcal K=2I`$, so
$`J'=\mathcal K\mathcal K^\top-I`$ is an exact reflection.
Replace $`(A,B)`$ by $`(R_x,0)`$ to define $`J_R`$, with
$`\|J'-J_R\|=O(2^{-s})`$. The two D-dimensional lanes
are indexed by helper label q zero and direction bit z.
For $`J_I`$ the same ideal reflection at $`R_x=I`$,
$`W'=J'J_I`$ is real dyadic orthogonal and approximates

```math
W_{\rm id}|_{q=0}=\mathrm{diag}(R_x^\top,R_x),
\qquad W_{\rm id}|_{q\ne 0}=I
```

in full operator norm at $`O(2^{-s})`$.

### A flat exact block encoding with one shared programme

Use the $`J=D'^2`$ literal real Pauli sections $`P_\alpha=X^uZ^v`$.
The coefficients $`a_\alpha=D'^{-1}\mathrm{tr}(P_\alpha^\top W')`$
are dyadic with $`O(s+h)`$ bits. The identity coefficient has
$`|a_0|\le1`$; every other one is at most
$`\|W'-I\|=O(\epsilon_0+2^{-s})\le1/J`$ in absolute value.
Fix the constants in the coarse and fine errors to ensure that bound.

Of two equally likely menu branches, the first emits $`\pm I`$
with mean $`a_0`$ using a dyadic threshold, and the second chooses
uniform $`\alpha`$. That branch emits the signed $`P_\alpha`$ with
probability $`J|a_\alpha|`$ at nonidentity labels. Otherwise it
emits independent uniform signed I, including the $`\alpha=0`$
label. Thus the uniform menu average is exactly $`W'/2`$.
Thresholds use $`O(s+h)`$ bits and an exact short R91 comparator
with $`O(\log(s+2))`$ T-depth and $`2^{O(s)}`$ private resources.
Ordered controlled X and Z batches implement the literal sections;
shared controls use the exact constant-depth
[shear schedule](PARALLEL_BILINEAR_WORD_QUERY.md#2-parallel-bilinear-outputs-and-their-controlled-version).
Include a separate selector bit for the uniform sign.

R91 displaces one shared arbitrary bank by the entire
$`O(J(s+h))`$-bit threshold table. A private lookup F for each slot's
coherent $`\alpha`$ reads only the corresponding word, by XORs of
its bank bits controlled on the label. These are exact short
conjunction phases. In chronological order $`F,Q,F,Q^\dagger`$
adds the true threshold into each slot's scratch and returns the
central bank. The F calls across slots share only stable bank reads
and have disjoint T targets. Unload with the actual reverse of the
same sequence. This loads one coefficient per slot without requiring
a clean copy of the full table. Each bank echo factors at fixed x on
the full inputs, even with coherent selectors.

Hadamards on the menu selector and their reverse give clean
compression $`W'/2`$. R99's exact oblivious amplification gives
$`-W'`$ with the entire selector, comparator and threshold bank
A of $`k=O(s+h)`$ bits returned to zero. Correct its minus sign
by a literal Clifford phase. The actual fine slot word W therefore
preserves its $`A=0`$ projector on the full space. The complete
slot array factors on x and uses constantly many long queries,
including the reverses in amplification.

## 4. Candidate selection and the local norm estimate

Use common dirty shift labels for the selection and the corresponding
masked execution. In one seed basis, group the slots into
$`g=\Theta(2^{|q|}(s+1))`$ groups, each with $`h_0=2^k`$ slots.
For each group form the XOR j of the binary indices i of all its A=0
slots. Its candidate is j if j itself has A=0. A qualifying candidate
also has q=0. Use the first qualifying group, so at most one slot is
chosen. The list of candidates before q qualification depends only
on invariant A-projector statuses, and its count is at most g.
If precisely one slot has A=0, its index gives a valid candidate.
This is the parity-index selection from the tensor-sector construction.

Let B_sel be the ideal operation applying Pauli $`X_z`$ to
the chosen slot orientation, identity if none. It tests q
only on those g candidate locations. Thus in $`T=W^\dagger B_{\mathrm{sel}}W`$, the local exact W's on every other slot
cancel entirely, including their arbitrary off-code actions,
conditional on any fixed A-status string. On candidates W
acts exact W'. Replace just those by W_id, with total operator
error $`O(g2^{-s})`$. The same bound holds on the coherent
direct sum over A statuses, x and seed labels. The native B_sel
below is exact.

Form $`F=T B_{\mathrm{sel}}`$. Its ideal comparison acts on
the selected q=0 slot as $`\mathrm{diag} (V_x^{1/2},V_x^{-1/2})`$ on z and the target payload, and
identity on all other slots. It preserves all original
masks exactly in the comparison. The norm bound involves only the g candidate slots.

### Exact native selection

Let $`a_{u,i}=[A_{u,i}=0]`$. Each group candidate predicate
tests $`a_{u,i}`$ and equality of i to the XOR of all indices
weighted by $`a_{u,i}`$. To qualify, also test $`q_{u,i}=0`$.
Use R99's two shifted indicators and four Clifford dot corners
on these O(s) bits, using a private dirty XOR query to form the
index parities with direct read-only conjunctions. One
candidate predicate or its exact Fourier-controlled XOR phase
has depth $`O(\log(s+2))`$ and cost/width $`2^{O(s)}`$.
These private calls can align their actual T layers, sharing
only raw key, q and Fourier character reads. Group success
$`g_u`$ is exactly the parity of these predicates, since at most
one qualifies. This prices its returned XOR into any arbitrary
target and the actual inverse.

For global priority take the candidate predicate and conjoin it
with all earlier groups' failure predicates, at most
$`g=\mathrm{poly}(L)`$ terms. Use the exact dirty counter
of DIRTY_SUM_COMPRESSION §§5–6. In each parallel private
controlled-increment word, replace the CNOT with control
$`g_u`$ by its returned XOR oracle above, including the actual
inverse where reversed. Neither the toggled line nor counter
needs initialization. Each control query has disjoint private
T targets and raw masks as read-only controls. The remaining
counter sum, offset cancellation and cyclic test has depth
$`O(\log(g+2))`$ and polynomial cost/width in g. The added fixed
number of XOR-oracle stages gives $`O(\log(s+2))`$ depth.
Include a fixed number of extra payload, Fourier character or
guarded-mode literals in the same conjunction. For each
requested sweep or success flag, perform all of its Fourier phase
terms with distinct private counters and oracles, align their
nonlinear layers, then undo the Fourier change on outputs.
This exactly realizes the B_sel label XOR, the complete success
XOR (parity of mutually exclusive final select predicates) and
the three sweeps of S. Each and its actual reverse have
$`2^{O(s)}`$ count/width and $`O(\log(s+2))`$ native depth,
without initialized work. Coarse mask selection uses this same
construction with q absent, g of order s, and A enlarged to
include P.

## 5. Extract and cover an arbitrary masked unitary

In the residual case use the actual
$`\mathcal E=S^\dagger F S`$, where S swaps
the main target+coarse-helper payload and $`c_1`$ into
the selected slot's corresponding payload and orientation
z. It never swaps A or q. In the ideal good sector,
the main sees $`V^{1/2}`$ at $`c_1=0`$ and $`V^{-1/2}`$
at $`c_1=1`$; the old slot returns exactly. With no qualifying slot
$`\mathcal E=I`$. The actual echo
```math
 (X_{c_1}\mathcal E X_{c_1})^\dagger\mathcal E
```
therefore compares to one full V on its main c1=0
sector on success, inverse on c1=1, all work returned.
Use only constant local and global words. The hybrid
is full norm since comparison of F was uniform across
all original A,q and the swap is exact.

For coarse $`U_0`$, use instead
$`\mathcal E=F^\dagger S^\dagger F S`$, with F equal
to its exact slot-controlled-local word (active when slot
z=0). S selects a first mask-valid slot and swaps same
payload+c1. F preserves A=P=0 each slot, so other slots
cancel and the chosen old slot inverse is independent
of c1. Main sees U0 at c1=0 and I at c1=1. The same
outer echo returns the old slot and yields desired masked
$`U_0`$ at c1=0. The exact selection is the above two-level
with q absent.

For each fixed original A,q, independently random shifts
on all A,q bits give each group a constant probability of
exactly one zero A-key. Its parity-index candidate then
qualifies with independent probability $`2^{-|q|}`$.
Thus total success
failure probability at most $`e^{-\Omega(g2^{-|q|})}`$.
For the coarse case q empty. Choose two independent
shift lists and a power-of-two dirty seed dimension
via exactly the R99 Chernoff/union proof, with failure
density $`\delta=e^{-C s}`$ and list size at least a
constant times the total raw-mask count divided by
$`\delta`$. It is $`2^{O(s)}`$. Opposite Z and X seed
bases yield bad projection overlap at most
$`2\delta\sqrt{\mathrm{dim}}`$, exponentially small
at fixed large C. R99's principal-angle correction
orthogonalizes the second bad space within this error,
uniformly on original mask sectors and arbitrary
references.

Let E1,E2 be actual masked orientation words for one
of the two target unitaries $`U_0,V`$, on the two seed
bases. On $`c_1=c_2=0`$ apply E1, then compute first-basis
failure as a complete XOR into $`c_2`$, apply E2 only
on $`c_2=1`$, and actual-uncompute the failure bit.
In the comparison the two bad projectors are orthogonal, so the
enabled branch applies the desired unitary exactly once and
returns both conditional bits. The actual error telescopes on
the complete physical space.

In the residual case implement the stable second-basis and outer
suffix controls on just the two $`B_{\rm sel}`$ centers in F,
by adding their literals to the exact local phase compiler.
When disabled $`F=W^\dagger W`$ and the S sandwich cancels
exactly. For coarse control only the local standalone S, also
by adding literals to its exact sweep predicates. When S=I,
the actual F†F cancels identically, and so does the c1 echo.
The shared original-address QROMs never need controls. On the
outer inactive branch E1,E2 are identities and the failure XOR
and its exact inverse cancel. Active branch uses the two known-zero
original suffix wires $`c_1,c_2`$ inside the external retained clean
guard. Their ideal return is followed by the actual inverse of
that external predicate.

## 6. Parallel capped geometric references

Use the [R100 capped reference](GEOMETRIC_CHUNK_DEPTH.md#4-parallel-capped-geometric-reference)
in each R99 bin, at its required truncated power-of-two size.
The cap uses a constant sequence of parallel digit preparations,
a reflected cap correction and the literal phase adjustments.
Use R99 programmed rotations at the original synthesis tolerance
for the angle calls in this range. The two dirty-port sectors
produce the same positive cap in the ideal initialized cylinder,
including arbitrary references. The transported projector is
therefore the original one tensored with identity on the borrowed
ports. The original analytic ray transport from the padded cap
to the full reference remains valid at the same truncation error.
Keep the exact original anchor tests and actual inverses. The
depth of these preparations is $`O(\log_2(L+2))`$ with separate
private banks for all simultaneous angle calls.

## 7. Simultaneous resources and complete return

For $`L\ge\sqrt n`$, the precision slack
$`s=L+O(h+\log(L+2)+\log(n+2))`$ is $`O(L)`$. Choose its constants
so each block's complete error is at most $`2^{-L}/(64(n+1))`$,
including its covered coarse and fine words. The ideal products apply
$`U_{0,x}`$ and then $`V_x`$, giving exactly $`U_x\otimes I`$.
The sign reduction, opening and bounded bottom use their retained
allocations with fixed precision slack. The initialized-isometry
hybrid applies to every logical column: the ideal prefix returns all
required work and each later actual word is unitary even on leakage.

Use $`b\ge\sqrt N/n`$ and dimensions above an absolute cutoff first.
Every standalone bank, swap, indicator, counter, seed lookup, and
comparator uses $`2^{O(s)}\mathrm{poly}(n)`$ total gates and peak
width, including all block occurrences and actual inverses. Fix
their exponents and the original opening constants before $`\xi`$.
These terms then fit below $`2^{n/4}`$ on $`L\le2\xi n`$ for small
enough $`\xi`$ and large enough n. Reserve a disjoint
$`\Theta(b)`$ R91 query pool. Local and query pools can be reused
across completed blocks, excluding all live logical and program
coordinates. No conditional clean bits are carried between blocks.

Each block makes a fixed number of long R91 calls, including
actual reverses, on x only. With volume $`Q_i=O(NL2^{-j_i/2})`$,
the geometric sums of $`Q_i`$ and $`\sqrt{Q_i}`$ give the
$`\sqrt{NL}`$ and $`NL/b`$ count terms, the $`NL`$ Clifford
term, and $`NL/b^2`$ depth. Query depth outside that volume
term is $`O(\log(n+2))`$ per block. Local coarse depth sums to
$`O(\sum_i h_i\log(h_i+2))`$, which is $`O(\ell\lambda)`$.
Local fine words, shifts and selection take $`O(\log(L+2))`$
per block, including suffix controls. The cap of Section 6
removes the serial reference angle positions in R99's opening.
The other opening words then contribute
$`O(N/b^2+\log(n+2)+\log(L+2))`$, with R99's original count,
guard and reserve geometry. The bounded deepest layers use R99
individually. This proves the stated simultaneous depth and counts.

For the prescribed complex magnitude completion, use the real branch
at $`L+1`$ with its corresponding literal allocation, followed by
R99's zero-clean addressed $`R_Z(\phi_x)`$ on the external clean bit.
On the ideal zero flag this supplies the prescribed literal phase,
including at zero-magnitude leaves. Use its full-operator error at
most $`2^{-(L+1)}`$; the two tolerances add to $`2^{-L}`$. The one
phase table has volume $`O(NL)`$ and a constant number of charged
queries. The internal allowance $`L+1\le 2\xi n`$ is why the
construction above used that wider interval. At $`b\lt \sqrt N/n`$,
R96 at the literal precision and width has its additive depth
absorbed by $`NL/b^2`$. Use the same fallback below the absolute
dimension cutoff.

Finally, writing $`U=W\otimes I_b`$, the norm of
$`V^\dagger J_1-J_1U^\dagger`$ equals that of $`VJ_1-J_1U`$
by multiplication by actual unitaries. Actual reversal preserves
T-depth and counts up to fixed alphabet factors, including replacement
of $`S^\dagger`$ by $`S^3`$. Both directions therefore obey the
same return and accuracy guarantee.

## 8. Small-clean menus through quadratic precision

The grouped construction also applies when $`\sqrt n\le L\le n^2`$.
With the same accuracy, work return and simultaneous counts, its depth is

```math
D_T=O\bigl(NL/b^2+(1+L/n)\ell\log_2(\ell+2)\bigr),
\qquad \ell=\log_2(n+2).
```

In particular this gives $`O(\log n\log\log n)`$ at $`L=\Theta(n)`$
and $`b=\Theta(\sqrt{NL})`$, with any fixed linear precision constant.
The real and complex circuits retain the literal R96 dirty thresholds.

### A dirty scalar in the residual menu

Keep the exact dyadic graph correction $`W'`$ and its literal real Pauli
coefficients $`a_\alpha`$. The first half of a uniform binary menu uses
$`a_0 I`$. The second selects a uniform label $`\alpha`$ and uses
$`J a_\alpha P_\alpha`$ for nonidentity labels, with a zero scalar on
its identity entry. All these scalars are dyadics in $`[-1,1]`$. Use the [exact grid and raw scalar block](../../docs/OPERATOR_SOURCE_COMPILER.md#the-scalar-block)
with an unpaired Majorana chain of length $`O(L+h+\log(n+L+2))`$.
The scalar is the compression on one local initialized flag, for all
arbitrary precision-core inputs. Apply the ordered $`P_\alpha=X^uZ^v`$
on its own payload with the menu label as stable controls. This makes the
compression of the complete menu exactly $`W'/2`$. Its
[R99 amplification](PROGRAMMED_ROTATION_DEPTH.md#an-exact-invariant-clean-source-with-a-short-program)
has action $`-W'`$ and preserves the whole initialized key sector.
The minus cancels in the selected graph conjugations. Only $`O(h+1)`$
local key bits are initialized, including the menu and scalar flags.

Tabulate the exact scalar sign words into a single program of width
$`O(J[L+h+\log(n+L+2)])`$ for each original outer prefix. Use the
complete R91 XOR into one arbitrary bus B. For each slot, a local
linear lookup $`F`$ reads the appropriate piece of B at that slot's
own coherent menu label and XORs it into its own arbitrary precision
core. Implement these local shears by read-only equality phases with
private T targets. Their joint depth is $`O(\log(h+2))`$ and their
local count is polynomial in the number of slots, $`2^h`$ and L.
The chronological echo $`F,Q_F,F,Q_F^\dagger`$ cancels the old bus
and gives the intended simultaneous displacements. Full-output
Hadamards make each completed displacement the required Z mask.
The complete menu uses only constantly many actual full-volume query
visits even when its private slot labels differ. Return and literal
phases are exact on arbitrary bus, core and selector inputs.

First form the complete exact comparison circuit using the old fixed
Majorana words. Then apply [R100 chunk substitution](GEOMETRIC_CHUNK_DEPTH.md#3-exact-geometric-chunks-and-the-powered-source-loader)
at completed fixed-source locations. This compares their full operators,
including the rejected action and all actual reverses. The old scalar
masks and every bus query stay unchanged. The new native word need not
preserve its initialized key exactly, since the full-operator hybrid
compares to the exact comparison word before using its slot identities.
In the inactive outer or second-basis branch disable the corresponding
actual selected center or selected swap as in the preceding construction.
Its surrounding actual words cancel literally on occupied work.

### A tensor cover at high accuracy

Use the [exact parity index and tensor priority](GEOMETRIC_CHUNK_DEPTH.md#1-exact-priority-on-arbitrary-occupied-keys).
For a subgroup of $`2^k`$ slots let $`a_i=[A_i=0]`$ and choose
$`j=\bigoplus_i i a_i`$. Qualify this sole candidate when $`a_j=1`$
and its helper label q is zero. In a coarse slot the helper label is
empty and the key includes the complete low-precision program. Uniform
translations give success probability at least $`2^{-(|q|+2)}`$,
since exactly one zero key suffices and its q label remains uniform.
Take $`c_{k,q}=\lceil C2^{|q|}(k+|q|+1)\rceil`$ subgroups per base
group. Use the explicit upper failure bound

```math
\rho_0=(1-2^{-(|q|+2)})^{c_{k,q}}.
```

For w raw key and helper bits, take the least power-of-two shift-list
size D above a sufficiently large absolute multiple of $`(w+1)/\rho_0`$. The count of failed labels for each
fixed input is stochastically bounded by a binomial of mean
$`D\rho_0`$. The Chernoff and union bound of R100 choose a list for
which every bad set has size at most $`2D\rho_0`$. Enumerating lists
and testing these integer counts terminates. The two seed bases use
the same physical seed register with Hadamards on that register alone.
Their bad overlap is at most $`2\rho_0\sqrt D\le1/4`$ for sufficiently
large fixed C. Here $`w,D\le 2^{O(k+|q|+1)}`$. No reciprocal of the
actual failure density is used.

Take $`O(L+\log(n+2))`$ independent base groups. The bad projectors
are tensor products of their respective group projectors, so the pair
overlap is exponentially small in the requested accuracy parameter.
The first-qualified-subgroup rule uses exact returned conjunctions of
complete parity-index flags, with one extra stable q-zero predicate.
Its local depth is $`O(\log(L+2)+|q|+\log(k+2))`$, and local count and
width are $`2^{O(k+|q|+1)}\mathrm{poly}(n,L)`$.
The graph comparison examines q only at the one key-valid candidate
per subgroup, as before. Its error is bounded by the number of subgroups
times the graph approximation error, irrespective of the number of
occupied slots. Choose the graph precision $`L+O(h+\log(n+L+2))`$
with sufficient fixed margins. All covers, mask reverses, failure
uncomputation and orientations use the same two returned conditional
suffix bits. The residual slots have $`k,|q|=O(h+1)`$.

### Complete-frame ledger

Use the R100 opening and bounded bottom at the original global prefix
addresses. Its remaining height is $`O(\log(L+2))`$. Use growing groups
as above, with height $`h\le j/C_0`$ at original suffix distance j and
$`h\le c_0\log_2(n+2)`$ for small fixed $`c_0\gt 0`$. The coarse program
key has size $`O(2^h(h+1))=o(n)`$. The residual slots have polynomial
count in n and L. The old exact coarse schedule has depth
$`O(h\log(h+2))`$ on disjoint private slots. There are
$`O(\log_2(\ell+2))`$ groups and their heights sum to $`O(\ell)`$.

Every program uses its full prefix of length n minus j minus h.
The same geometric envelope $`O(NL2^{-j/2})`$ holds for the complete
word volume of each block, since the program factor is
$`2^{O(h)}\mathrm{poly}(h+1)L`$. Its volume and square-root
sums pay all queries. The additional private bill including the coarse
slots, tensor covers and accelerated residual loaders is bounded by
$`2^{o(n)+O(t)}\mathrm{poly}(n,L)`$ for chunk length
$`t=\lfloor\delta n\rfloor`$ in the wide branch. Choose $`\delta`$
small as in R100. This fits next to the original exclusive query pool
when $`b\ge\sqrt N/n`$. Use R96 on narrower widths and bounded n.
Each residual source has $`O(1+L/n)`$ sequential chunks of logarithmic
native depth, at precision $`L+O(\log(n+L+2))`$ with enough slack for
the polynomial number of appearances. The mask depths add one $`\ell`$
per block plus the width term. The coarse depths sum to
$`O(\ell\log_2(\ell+2))`$. This gives the displayed depth with the
original T count, Clifford count and literal width. All accountings
include actual reverses. Use a slightly enlarged internal polynomial
range for the complex $`L+1`$ shift and the R100 literal phase word.

The [three bounded fixtures](../../tests/test_grouped_programmed_precision.py)
check exact dyadic graph orthogonality on every auxiliary lane,
the selected residual echo and inactive branch with actual adjoints,
and coherent coefficient displacement from an arbitrary shared bank.
The all-size gate schedule, error and width statements are analytic.
