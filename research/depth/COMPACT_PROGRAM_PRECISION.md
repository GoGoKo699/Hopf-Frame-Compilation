# Compact arithmetic programmes for complete precision tails

[Grouped programmed precision R101](GROUPED_PROGRAMMED_PRECISION.md) ·
[Geometric chunks R100](GEOMETRIC_CHUNK_DEPTH.md) ·
[Programmed rotations R99](PROGRAMMED_ROTATION_DEPTH.md)

**Theorem R102.** For an absolute $`c>0`$ and every integer precision
$`6\le L\le2^{cn}`$,
the two prescribed complete Hopf frame families admit one external
initialized qubit, the literal dirty thresholds
$`b\ge17(L+n+7)`$ real and $`b\ge17(L+n+8)`$ phase-dressed
complex, and full initialized-isometry error at most $`2^{-L}`$.
With $`N=2^n`$ and $`\ell=\log_2(n+2)`$, the simultaneous resources are

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O\bigl(NL/b^2+(1+L/n)\ell\bigr).
```

The actual reverse satisfies the same bounds for the inverse frame.
The error includes every logical column, coherent addresses, arbitrary
dirty inputs and references, and all work leakage and return.
In particular, $`b=\Theta(\sqrt{NL})`$ gives $`O(\log n)`$
native depth uniformly through linear precision. At higher precision
the additive term grows as $`(1+L/n)\ell`$, including at
$`L=2^{\Theta(n)}`$ inside the stated range. R96 retains coverage at
all precisions, and the fixed-accuracy gap persists.

A whole tail block reads only its original angles and a coarse
field program. A finite-field eigenstate permits constant native
T-depth per coarse layer after parallel preparation. The fine
correction calculates its Pauli coefficients from the same compact
row by returned dirty arithmetic of logarithmic native depth.
A fixed number of global word calls serves all arithmetic lanes,
slots and coefficient labels. This permits a constant number of
groups at the original table volume through $`L\le n^2`$.
For larger L, Section 6 promotes the entire fixed geometric source
after reuse of its short initialized key.

## 1. Constantly many blocks with the original address geometry

First work on $`L\le n^2`$. For small L with $`L\le2\ell`$, the
[R96 bound](ALL_PRECISION_LINEAR_DEPTH.md) has the desired depth. If dimensions
are below an absolute cutoff, the original R96 theorem gives the claim
because L is bounded by their squared maximum. For
$`b\lt\sqrt N/n`$ use R96 at the original width and precision,
since its additive depth is absorbed by $`NL/b^2\gt n^2L`$.
The rest of the proof uses the wide allocation above a fixed cutoff.

For $`L\gt2\ell`$ use the same sign reduction and three opening
contractions as R99 or R100. Use the parallel capped reference
of R100 with R99 rotational synthesis for
$`L\le c_*n/2`$, where the fixed $`c_*`$ is small enough for the
R99 exponential local bill including constant precision shifts.
At larger L use R100's accelerated fixed loaders in the opening.
The two branches overlap after increasing the absolute cutoff.
The opening has native depth
$`O(N/b^2+\ell+(1+L/n)\log_2(L+2))`$
and leaves $`M=O(\log_2(L+2))`$ original bottom layers.
Keep a fixed number of lowest levels for individual rotations,
as in R99 and R100. Split the other levels into blocks of
height at most $`\max(1,\lfloor c_0\min(\log_2(L+2),\ell)\rfloor)`$
with fixed small positive $`c_0`$. There are constantly many
blocks because $`L\le n^2`$.

For a block of h levels, j original lower bits stay outside its
address of length n minus j minus h. Keep that entire address.
Compute the exact lower-zero predicate into the external clean bit H.
There are at least two such lower bits for each nontrivial group.
Borrow two of them, $`c_1,c_2`$, as initialized bits only on
H one. The orientation and failed-sector selection use the
[R101 echo and two-basis cover](GROUPED_PROGRAMMED_PRECISION.md#5-extract-and-cover-an-arbitrary-masked-unitary).
Their actual control will be implemented at the slot swap, so the
body is exactly inactive on H zero with all old work occupied.
The separate suffix predicate and its actual inverse retain
their original literal allocation.

## 2. Coarse field phases in O(h) local layers

Fix a block height h and let $`m=C_0(h+1)`$ for a large absolute
integer constant. In $`\mathbb F_{2^m}`$ choose a primitive element g,
put $`Q=2^m-1`$ and $`\zeta=e^{2\pi i/Q}`$.
Round each local real angle at each preserved outer prefix to a
multiple $`2\pi d/Q`$ and store $`\beta=g^d,\beta^{-1}`$
in the compact classical row. Certified bounded-error rounding
with one extra slack avoids tie decisions. The ideal complete grid
block $`U_0`$ differs from the target block by $`O(h/Q)`$ in
operator norm. Its field programme has $`O((h+1)2^h)`$ bits
per row. These field elements are data, never new global addresses.

The reference for a local node is the initialized state

```math
|\chi\rangle=\frac{1}{\sqrt Q}\sum_{j=0}^{Q-1}\zeta^{-j}|g^j\rangle .
```

Multiplication by $`g^d`$ has its literal eigenvalue
$`\zeta^d`$, and its inverse has the conjugate eigenvalue.
In a Clifford basis where the target rotation is about Z, apply
the multiplier or its inverse according to the target eigenbit,
with signs chosen for the convention $`R_y(\theta)=e^{-i\theta Y}`$.
The completed ideal gate returns the reference exactly.
The field elements are nonzero on the valid program. A temporary
m-bit zero register suffices for multiplication. XOR
$`\beta r`$ into the temporary, XOR the inverse multiplier of
its output into r, then swap.
All three operations have the same stable guard, including the swap.
Implement the latter by three guarded XOR sweeps. Each field XOR is bilinear in
the stable data and the supplied multiplier. In the output
Fourier basis, each coefficient contributes a sign conditioned on
two binary inputs and the Fourier character. Include the constant
number of stable guard literals in this predicate. The exact
short conjunction compiler uses private lanes and returns all
helpers. Thus all predicates of a shear can run in constant T-depth
with a fixed polynomial number of Cliffords, T targets and dirty
wires per field multiplication. The same schedule runs across
distinct nodes with the common target eigenbit read-only.

Here is a phase reference preparation with the required full-space
comparison. Start with the uniform column $`\psi`$ on all
$`2^m`$ binary indices. Let t be the last index and
$`p=2^{-m}`$. Phase t by $`e^{i\phi}`$ and phase the
ray $`\psi`$ by $`e^{i\theta}`$. Put
$`w=1-p+pe^{i\phi}`$ and choose

```math
\cos\phi=1-\frac1{2(1-p)},\qquad
e^{i\theta}=1-\frac{e^{i\phi}}{w}.
```

The first equation gives $`|w-e^{i\phi}|=|w|`$ with w nonzero.
The last coefficient is then
$`\sqrt p[e^{i\phi}+(e^{i\theta}-1)w]=0`$.
All others have the same magnitude and phase. Implement the two
phases with a private arbitrary dirty port, using its two Z
eigenspaces for conjugate phase signs as in the R100 cap.
Correct the common phase by a prescribed Z rotation on that port.
The ideal prepared cylinder is therefore the positive uniform
column supported on indices below Q, tensored with the identity on
all arbitrary dirty work. Its ideal extension can act differently
off that cylinder. The predicate for the uniform ray is the
last or zero computational test after fixed Hadamards, as appropriate.

Next phase each binary index by $`\zeta^{-j}`$ using the
parallel product of single-bit Z rotations. Its known common scalar
is corrected by the literal minus identity, since the binary
place values sum to Q. Apply the fixed permutation $`j\mapsto g^j`$
for j below Q, and $`Q\mapsto0`$. One temporary initialized
m-bit register and the exact forward and inverse short word
queries implement it cleanly on this code, with its full arbitrary
off-code action retained in the comparison. The short queries and
conjunctions cost $`2^{O(m+1)}`$ gates and width and
$`O(\log_2(m+2))`$ T-depth. All angle calls are a constant
number of parallel banks of accelerated fixed rotations from R99
or R100, at the error specified below.

Supply one local reference per tree node in a coarse virtual slot.
Prepare them in parallel with independent cores, and reverse the
actual preparations after the whole coarse product. In the
real-rotation order, bits already rotated are preserved by later
layers. Compute unary prefix flags forward along the current
observed prefix, two children per node using constant-arity
conjunctions once per layer. They can be reversed from the final
logical prefix, after the last layer, using the actual reverse.
Precompute each lower-zero flag in parallel before any layer.
After a layer and before the next logical bit changes, erase its
lower-zero flag using the exact conjunction of the next lower-zero
flag and the next logical bit being zero. At the bottom erase
the empty-test flag by X. These erasures are literal XORs and
coincide with the intended reversals on every initialized column.
There are O(h) such chronological positions, each of constant
T-depth, disjoint from the parallel reference banks.
Use the private orientation z only as an extra guard at the
field multipliers. If z is one, the ideal code action is identity.
The prefix and suffix work still returns on that ideal code.
The whole physical word stays unitary on arbitrary occupied
slots and remains eligible for an actual reverse.

The complete coarse comparison on the initialized slot code is
$`\mathrm{diag}(U_0,I)`$ with its key returned. The initialized
key comprises the original zero field programme, the reference
fields and their permutation and multiplication temporaries,
and the prefix and suffix flags. It has $`k_c=O((h+1)2^h)`$
bits. All precision cores and private returned helpers are
arbitrary and excluded from that key. An exact shared XOR
displacement loads the field row into each originally zero
program via the R99 Clifford copy echo. Its actual reverse
unloads it after the comparison block. The actual slot replaces
only the fixed reference angle words by their native
approximations, so a full-space hybrid compares the actual and
ideal slot including all occupied keys. No assumption of actual
exact code preservation enters that hybrid.

## 3. Returned threshold arithmetic on the shared row

Let stable public registers carry a computational label $`w`$. For each
Boolean function $`f_i(w)`$, assume an exact complete XOR oracle into
any separate arbitrary target. Each oracle returns all private work on coherent public inputs and
arbitrary references. Let $`H=\sum_{i=1}^{p} f_i(w)`$ and consider the threshold
$`[H\ge t]`$. Constant thresholds and $`p=0`$ use fixed gates or identity.
In the following $`p\ge1`$ and $`0\le t\le p+1`$.

Use the output Fourier basis, with its bit $`\chi`$ unchanged throughout
the diagonal computation. Append $`p+1`$ occurrences of the literal
$`\chi`$ to the weight inputs. Occurrences refer to the same read-only
wire and do not require copies. Thus the new sum is
$`H'=H+(p+1)\chi`$, with $`0\le H'\le 2p+1`$, and

```math
[H'\ge t+p+1]=\chi[H\ge t].
```

Set $`m=\lceil\log_2(2p+2)\rceil`$ and $`M=2^m\gt 2p+1`$. We implement
a complete shift $`Q_{H'}:C\mapsto C+H'\pmod M`$ of an arbitrary m-bit
counter C. In the read-only two-addition increment of
[Parallel dirty lookup §6](../../docs/PARALLEL_DIRTY_LOOKUP.md#a-read-only-controlled-increment-using-two-additions),
let $`D`$ add a private arbitrary m-bit word g into a private
arbitrary m-bit addend a, and let $`E`$ XOR the desired input
$`f_i(w)`$ into the low bit $`g_0`$. Use the exact chronological word
```math
 C_g,\ D^\dagger,\ E,\ D,\ E^\dagger,\ C_g^\dagger ,
```
where $`C_g`$ complements a under $`g_0`$. The middle four operations
translate a by $`f_i(w)(1-2g_0)`$ and restore g. Conjugating by the
conditional complement cancels that sign and translates by $`f_i(w)`$.
Each actual E is the full child XOR, with private work distinct from a
and g and all public inputs preserved. The g bits for distinct
increments, counter outputs, and helpers are disjoint. All original
$`\chi`$ occurrences enter as Clifford controls. This is the
two-addition increment with no other use of the public literal.
The separate two-helper increment in that chapter has a different
control schedule.

Let J perform these parallel increments. Let S be the exact returned
sum of the initial addends into C, with the input-independent work offset
of [Dirty sum compression §§5–6](../../docs/DIRTY_SUM_COMPRESSION.md#5-returning-the-work-and-cancelling-its-offset).
Apply $`S^\dagger,J,S,J^\dagger`$ chronologically. The sum offset and
the arbitrary initial addends cancel, giving $`Q_{H'}`$ exactly and
restoring every addend and helper. Outside the child XORs this costs
$`O(m)`$ native T-depth and polynomial cost and dirty width in p.
There are four chronological batches of child XORs counting actual
adjoints.

Allocate separate arbitrary M-bit words U and V. Use the
[exact dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md)
to define A adding $`e_C`$ to U, and B adding $`e_{C+H'}`$ to V.
For B, apply $`Q_{H'}`$, the indicator, and its actual shift
inverse. Both completed shifts restore C and their private work.
Let $`\mathsf C`$ apply the fixed Clifford phase
```math
(-1)^{U^\top K V},\qquad K_{ij}=[(j-i)\bmod M\ge t+p+1].
```
This is a list of at most $`M^2`$ elementary CZs. The chronological
word
```math
 \mathsf C,A,\mathsf C,B,\mathsf C,A^\dagger,\mathsf C,B^\dagger
```
returns U,V and adds only $`(-1)^{e_C^\top K e_{C+H'}}`$.
Since $`H'\lt M`$, this is exactly
$`(-1)^{[H'\ge t+p+1]}=(-1)^{\chi[H\ge t]}`$, with no dependence
on the arbitrary old C. Output Hadamards therefore turn this phase
into the required complete XOR. The two B occurrences use at most
sixteen batches of child XORs or actual adjoints. Literal permutation
and phase identities give the operator equality on the full physical
space. No intermediate counter is assumed zero.

### Fixed-depth composition with one shared table bus

An a-bit address x selects a prescribed q-bit raw word r(x).
Take one shared arbitrary q-bit bus B, a disjoint exclusive R91 work
pool of $`B_Q\ge C_Q(a+2)`$ bits, and disjoint private outputs and
workers for every requested logical copy of a Boolean circuit.
In addition to r(x), each such circuit can read stable local inputs
on its own slot and stable public controls. Its output is a complete
bit XOR, possibly into arbitrary, entangled target bits.

Every leaf input of the unfolded circuit is a fixed position in
the raw row, a stable non-table bit, its complement, or a constant.
Simultaneous table leaf visits are Clifford shears F from the
corresponding bus coordinates into their private arbitrary targets.
Chronologically $`F,Q_r,F,Q_r^\dagger`$ XORs the actual raw row bits
into these targets and restores B and the query pool. The other
leaves are exact Clifford words. This paid leaf interface never
initializes a copy of the raw data. The programs may differ across
slots, provided all leaf indices are specified in their circuits.
A variable selection is expressed as Boolean gates above the leaves.

For a circuit of unweighted thresholds at depth d, recursively replace
each child XOR in the threshold construction above by the same construction. Unfold shared
Boolean gates as needed. Run all independent threshold words as one
batch, allotting separate private counters, indicators, and phase
targets. Stable public inputs, including any Fourier characters inherited
as literals, are read-only. Align their fixed sixteen child-query
positions, padding local arithmetic to the worst length. A lower-depth
threshold circuit can be padded with single-input identity thresholds.
Child actual inverses have the reverse schedule of the corresponding
forward words. All local phases and permutations use private T targets.
Any overlapping address, selector, or character interactions occur as
Clifford controls, whose elementary gates all count.

For definiteness let P majorize q, the number of stable local input
bits, the number of output gates, the number of gates, and maximum
fanin in one circuit. For R independent simultaneous copies, the
threshold construction unfolds to at most
$`R(P+2)^{C(d+1)}`$ local gates and private width for an absolute C.
To see the resource bound, each gate uses a polynomial in its fanin
of dirty wires and elementary gates from the counter and addend arrays, two
indicators, dense Clifford form, and sum-compression work. It calls its
children in a constant number of batches, with at most $`2P+1`$
simultaneous child instances per parent, including the repeated
Fourier control literals. The number of complete leaf batches is at
most $`16^d`$, independent of R and P, with at most two actual raw
queries per leaf batch, including inverses. At each level the local
depth outside child words is $`O(\log_2(P+2))`$. Thus, for fixed d,
the simultaneous resources obey
```math
\begin{aligned}
T&=O_d\!\left(\sqrt{2^a q}+2^a q/B_Q\right)
       +R(P+2)^{O(d+1)},\\
G&=O_d(2^a q)+R(P+2)^{O(d+1)},\\
D_T&=O_d\!\left(\log_2(a+2)+2^a q/B_Q^2+\log_2(P+2)\right),\\
w&\le B_Q+q+R(P+2)^{O(d+1)} .
\end{aligned}
```
The width excludes original public registers and the requested XOR
targets. The statement allows simultaneous exact adjoints through
reversing the emitted native word. It uses no initialized evaluator
work. It uses the R91 convention that Clifford interlayers need
no T-depth, with their full elementary gate count charged.

## 4. Fine correction from the raw angle row

Take one local h-bit complete real block at outer prefix x, and put $`D=2^h`$ . Let $`U_x`$ be
the prescribed local block and $`C_x`$ the exact ideal coarse block from the cyclic-field
angle grid. The coarse approximation satisfies $`\|C_x-U_x\|\le2^{-C(h+1)}`$ for a
sufficiently large absolute C. The desired exact residual is $`V_x=U_x C_x^\dagger`$ . With
slot orientation z the target becomes
$`W_x=|0\rangle\langle0|_z\otimes V_x+|1\rangle\langle1|_z\otimes I`$ . It is real orthogonal
and near I on the entire payload. For the real literal Pauli sections $`P_\alpha=X^u Z^v`$ ,
write $`W_x=\sum_\alpha a_\alpha P_\alpha`$ , with J labels. We have $`|a_0|\le1`$ and
$`J|a_\alpha|\le1/2`$ for nonidentity labels, because $`J=2^{O(h+1)}`$ and the coarse error
constant C is chosen first.

Tabulate certified dyadic sine and cosine approximations for both the original and coarse
angles, plus the coarse field multipliers and inverses. Use precision
$`S=s+C_S(h+\log_2(h+2))+C_S'`$ , where $`s=2L+C_s\log_2(n+L+2)+C_s'`$ . Each chronological
addressed rotation is then a dyadic matrix within $`O(2^{-S})`$ of its exact unitary in
operator norm, by direct sum across its controls. Leave inactive entries literal. Its norm is
at most $`1+O(2^{-S})`$ . Taking products of O(h) such factors and their transposes gives a
dyadic comparison matrix $`\widehat W_x`$ with $`\|\widehat W_x-W_x\|\le 2^{-s}/(10J)`$ for
the fixed constants above. This uses only accuracy of the dyadic data; the rounded matrices
need not be unitary and are never applied as gates. Each Pauli coefficient differs by at most
the same norm. The entire raw row has $`O(DL)`$ bits for $`L\gt 2\ell`$ , because the local
complete block has $`D-1`$ angle nodes and $`S=O(L)`$ .

A menu uses c equal to $`a_0`$ on its first half and $`J a_\alpha`$ on its second, with zero
on the latter's identity entry. Compute the corresponding dyadic $`\widehat c`$ from
$`\widehat W_x`$ , and set

```math
K_c=\min(2^s,\max(-2^s,\lfloor 2^s\widehat c+1/2\rfloor)),\qquad \widetilde c=2^{-s}K_c .
```

Projection to the interval cannot enlarge the error since the exact c lies in $`[-1,1]`$ . The
total error is at most $`2^{-s}`$ , even at endpoints. An unpaired raw Majorana scalar grid of
length $`s+2`$ supports every one of these dyadics. Encode its literal sign mask by binary
$`1-\widetilde c`$ on that signed grid using
[the separate all-minus endpoint](../../docs/OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients)
. Each sign is computed exactly.

For the online arithmetic, fix logical input and output basis labels. In a single
chronological block each logical target changes only once; its logical path is fixed by the
two endpoints because processed bits already have their output labels and unprocessed bits
retain their input labels. Each dyadic block matrix entry is therefore one product of at most
h prescribed local matrix entries, with its local active suffix determined by the fixed
labels. The entry of the product of two blocks sums over the D intermediate labels. Trace
against a Pauli adds at most 2D signed entries because Pauli sections are monomial. Every
coefficient is a signed sum of at most $`O(D^2)`$ products of $`O(h+1)`$ exact dyadic inputs.
Computing all coefficients and selecting with the stable coherent Pauli label needs only
polynomial size. All endpoint paths and their input field indices are offline constants. There
is no loop over them in the native word.

With common fixed denominators, products and sums are exact signed binary integer arithmetic
on $`O(hS+h)`$ bits.
[Hesse, Allender and Barrington](https://people.cs.rutgers.edu/~allender/papers/division.pdf)
give polynomial-size constant-depth threshold circuits for iterated binary multiplication and
addition. Sign handling, shifts by fixed powers of two, clipping, fixed-width comparisons,
selection and the final binary coding are constant-depth polynomial-size threshold arithmetic
as well. The whole formula has polynomial size in n,L in the block regime. Use Section 3 to
compute exactly the complete sign-mask XOR on each private scalar core, with only constantly
many full-volume queries of the entire raw row. Exact core Hadamards turn that XOR into the
literal Z phase. Uniform compilation here means the classical process terminates. It is not a
bound on preprocessing runtime for the native synthesis.

Form the exact raw-scalar menu using the ordered Pauli $`P_\alpha`$ in its second half and the
known half-selector Pauli controls. Its complete block differs from $`W_x/2`$ by at most
$`2^{-s}`$ in operator norm, by summing the *averaged* per-menu errors. No multiplicity factor
J is lost. Apply R99's actual amplification and the literal minus correction. The initialized
compression then differs from $`W_x`$ by $`O(2^{-s})`$ from its cubic singular-value
polynomial. Because the whole emitted word is unitary, its off-code leakage is at most
$`O(2^{-s/2})`$ . This already suffices at the chosen s; the linear normalization-two estimate
may also be used. By the
[complete block extension lemma](SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison)
, the actual amplified word is within $`O(2^{-s/2})`$ in full operator norm of an
analysis-only unitary whose exact whole-code action is $`W_x`$ and whose initialized key
sector is invariant. The result holds for arbitrary cores and references. Its residual key has
only $`O(h+1)`$ clean bits, including menu label and raw scalar flag, since all threshold
arithmetic and its outputs are the complete dirty phase masks. There is no graph-helper label
or graph midpoint. Original masks and separate retired helpers remain arbitrary.

## 5. Actual orientations, error and the complete ledger

Use the same two orientation and failure bits from the original lower suffix on the active
external guard. In either family the ideal slot word F preserves its initialized key A, and on
A=0 applies $`\mathrm{diag}(B_x,I)`$ in the slot orientation, with $`B_x=C_x`$ or $`V_x`$ .
Choose a valid slot by R100's exact parity index in each subgroup and exact first-success over
the tensor groups, with no helper q test now. Use both seed bases on the same physical seed
bank. With $`O(L+\log(n+2))`$ independent groups the bad projector overlap is exponentially
small. Swap the main logical payload and the initialized orientation $`c_1`$ into the selected
slot. As in grouped extraction, use $`E=F^\dagger S^\dagger F S`$ and
$`(X_{c_1} E X_{c_1})^\dagger E`$ . On success this is precisely $`B_x`$ on $`c_1=0`$ , with
every slot, A and borrowed helper returned. The old-slot factor cancels because it is
independent of the external $`c_1`$ . On failure it is identity. Use the second conditional
suffix bit $`c_2`$ for first-basis failure and second-basis control, then actually undo its
computation. Orthogonalize the second bad projector analytically as in R100 before asserting
ideal return. At each actual inactive branch disable all selected swap sweeps, so the actual F
and reverse cancel regardless of approximation, occupancy or possible source leakage. Neither
original prefix queries nor catalyst preparation need a shared external control.

For a coarse key of length k, all physical V and actual-reverse appearances across the slots
and tensor groups number $`2^{O(k+1)}\mathrm{poly}(n,L)`$ . Use catalyst synthesis precision
$`L+C(k+\log(n+L+2))+C'`$ with sufficiently large fixed constants, including parallel digits,
literal phase corrections and actual inverses. The whole coarse bank hybrid then costs only a
fixed small fraction of $`2^{-L}`$ , even including occupied off-code inputs. The conditioned
clean action is exact only in the comparison. The residual bank has $`k=O(h+1)`$ , so
appearances are only polynomial. Its block-extension errors $`O(2^{-s/2})`$ sum below the
accuracy slack, and its fixed-M replacement errors are assigned the same ledger, with
per-chunk precision $`L+O(\log(n+L+2))`$ . The complete bus, arithmetic and masks are exact
even when later work has leaked. The comparison residual and coarse blocks compose to
$`V_x C_x=U_x`$ exactly. All use the exact same two conditional suffix bits, returned in the
ideal block before the actual suffix predicate inverse. Earlier approximation error propagates
through full actual unitaries.

A group with lower suffix j uses exactly the original n minus j minus h prefix bits, with
constantly many optional guard literals. Its compact raw row has $`O(DL)`$ bits, so the entire
row volume is $`O(N L2^{-j})`$ . The coarse bulk programme uses only the multiplier portion at
width $`O(Dh)`$ , still bounded by the same envelope. Constantly many actual global queries
per group, including every returned threshold leaf echo, sum to $`O(\sqrt{NL}+NL/b)`$ T,
$`O(NL)`$ Clifford and $`O(NL/b^2+\ell)`$ query depth. The exact external saved-suffix
predicate only reads j original lower coordinates and returns with its actual inverse. The
constant bottom and literal complex phase use the accelerated R100 words with their unchanged
full masks.

Coarse k is $`O(2^h(h+1))=o(L)`$ and o(n), since h is capped by both logs. Coarse synthesis
including its all-slot accuracy charge is O(L), and residual scalar length and error precision
are O(L). Use chunks up to $`t=\max(1,\lfloor\delta n\rfloor)`$ , giving per-slot loader count
$`O(1+L/n)`$ . Standalone counts and dirty peak over all added stages, including the raw bus,
every private threshold oracle, field preparation, tensor groups and actual inverses, are
$`2^{o(n)+O(t)}\mathrm{poly}(n,L)`$ . Fix all exponents and choose small absolute $`\delta`$
so this fits under $`2^{n/4}`$ after a cutoff, disjoint from the original order-b query pool.
The local threshold delay is $`O(\log(n+L))=O(\ell)`$ . The coarse local stage uses $`O(h+1)`$
native T-depth plus its parallel prep and actual reverse. Both the latter and the residual
loaders have depth $`O((1+L/n)\log(L+2))`$ . The same bound covers exact tensor priorities and
payload returns. With constantly many groups this proves the displayed depth at the original
simultaneous counts and literal thresholds. Allow the internal polynomial range up to $`2n^2`$
for the real $`L+1`$ complex stage. Use the original one clean work flag for the final literal
phase and charge the actual native inverse by the same full-isometry norm and positive-S
constant factor as R100.

## 6. Promote the complete geometric word once

### Reuse the clean selector before removing it

Use the exact clean-code menus of [R100 Section 3](GEOMETRIC_CHUNK_DEPTH.md#3-exact-geometric-chunks-and-the-powered-source-loader).
For any short chunk that splits input Majorana $`\Gamma_j`$ along
at most t successive axes, set $`G'=(\Gamma_a+M_{\rm chunk})/\sqrt2`$
using a spare anticommuting axis. The exact clean-code action of the
R100 selector and amplitude amplification is $`-G'`$, with all
selector and arbitrary work returned. Follow a literal Clifford
plane rotation C sending $`\Gamma_j`$ to $`\Gamma_a`$ by
$`\Gamma_a`$ and then that exact code word. On the initialized
code the unitary is $`-G'\Gamma_a C`$. Conjugation carries
$`\Gamma_j`$ to $`M_{\rm chunk}`$ and fixes each original
outside axis by the R100 identity. The fixed minus is irrelevant
in this conjugation. The original core frame changes are fixed
Cliffords with the same elementary count as in R100.

Compose all these exact code words **before** using tensor promotion.
Their selectors can use the same padded zero key of length $`k=O(t+1)`$
because every exact completed menu returns its entire initialized key
on the whole code. The resulting loader acts exactly on the entire
initialized key cylinder, on arbitrary payload and private work. Keep
the original constant seed splits and both tails for the paired
source. For a source using O(m) splits the composed exact-code
loader U has

```math
U J=J(U_*\otimes I),\qquad U_*X_0U_*^\dagger=M_{\rm old}\otimes I_{\rm spare},
\qquad D_T(U)=O((1+m/t)\log_2(t+2)).
```

Here identities include the new arbitrary spectator axis and all
returned scratch; the selected payload is the entire O(m)-wire
old core plus that spectator. Elementary counts and width are
$`2^{O(t+1)}\mathrm{poly}(m+1)`$. This comparison concerns
the completed source action, without identifying the old whole
loader's other columns. Its actual reverse is a code-invariant
word with code action $`U_*^\dagger`$.

Apply [R100's selected-center and tensor cover](GEOMETRIC_CHUNK_DEPTH.md#2-a-selected-center-and-the-two-seed-bases)
once with $`V=U^\dagger`$, center Pauli $`B=X_0`$, and spare
mirror Pauli that anticommutes with every original source axis.
Thus the target is exactly $`M_{\rm old}`$ on every arbitrary core.
The physical word approximates $`iM_{\rm old}`$ to $`2^{-s}`$
with native depth

```math
O\bigl((1+m/t)\log_2(t+2)+\log_2(s+2)\bigr),
```

and count and width $`2^{O(t+1)}\mathrm{poly}(m,s+1)`$.
The repeated tensor bases no longer occur inside each short chunk.
Every physical key in this promoted word is arbitrary.

### Exact guarded source and complete spacetime ledger

For the original controlled source $`C_F(M_{\rm old})`$, place its
unchanged predicate F only in the selected Pauli-center predicates.
Its raw controls are stable and lie outside the precision core and slots.
Keep the mirror unchanged. On F false all actual selected rotations
are identities, so the mirrors and all surrounding actual unitary
words cancel. On F true correct $`iM_{\rm old}`$ by the literal
phase $`(-i)^F`$, using the constant-length native derivative for
stable Boolean scalar phases in R100 Section 2. If F has many
literals, add it as one complete exact returned predicate oracle in
the dirty-conjunction construction. Do not enumerate its input
strings among the private keys. The old arena/suffix conjunctions
cost only polynomial resources and logarithmic native T-depth.
The resulting controlled source has exact false-sector identity and
full-operator error $`2^{-s}`$ on the whole space. A fixed conjugate
signal uses the same elementary Clifford Q mask of R100 Section 3,
with its actual reverse. Use the actual reverse of the new promoted
source at inverse positions. The
original X/Z coefficient queries stay in
their exact order, with their full unmodified volume and addresses.

Use R102 through $`L\le n^2`$. Above that boundary retain the
original [R100 opening and direct-prefix groups](GEOMETRIC_CHUNK_DEPTH.md#7-groups-with-the-original-prefix-grammar)
including their fixed mixed gauges, masked sources, literal phase,
actual reversals and error budgets. Substitute the new complete
guarded source at the same fixed-loader boundaries in that
comparison word. Choose $`t=\max(1,\lfloor\delta n\rfloor)`$
with the same small absolute chunk fraction. The extra banks are
still $`2^{O(t)+O(\ell^3)}\mathrm{poly}(n,L)`$ across all
independent live instances. Choosing a sufficiently small c and
$`\delta`$ fits them in the same wide branch, disjoint from all
full-volume query pools. Count the coarse/fine queries verbatim,
and use the same literal-threshold R96 fallback below
$`b=\sqrt N/n`$ and at bounded n. The all-source hybrid uses
$`s=L+O(\log(n+L+2)+\ell^3)`$ with sufficient fixed margins;
this error allocation also covers every actual inverse.

The sum of fixed-source widths in a group remains $`O(L+P(h))`$
with at most $`P(h)`$ appearances for a fixed polynomial P. The
new loader term is therefore $`O((1+L/n)\ell+P(h)\log(L+2))`$
plus the logarithmic original guards. The opening and the final
bounded bottom and phase have the same improved loader term.
For the original R100 small group height, the coarse instruction
schedule and all polynomial-in-h factors cost at most
$`O(L^{1/8}\ell)`$. The latter is absorbed into $`O((L/n)\ell)`$
when $`L\gt n^2`$. The original query volume and square-root sums
remain unchanged, giving the claimed simultaneous resources.
For phase-dressed frames take the real precision shift inside a
slightly larger allowed exponential range and retain the actual
literal complex phase and its initialized flag from R100. The
forward full-isometry estimate also holds for the actual reverse.

The [four bounded fixtures](../../tests/test_compact_program_precision.py)
check the literal field character and multiplication with its actual reverse,
the prefix and suffix flags over every small logical column, and the returned
circulant threshold with arbitrary counters and both Fourier sectors.
The all-size arithmetic and frame error bounds are analytic.
