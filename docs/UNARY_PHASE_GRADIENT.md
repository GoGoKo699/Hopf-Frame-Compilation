# A conditional unary phase source for an entire Hopf group

[Grouped program interface](GROUPED_PROGRAM_PREFETCH.md) · [Exact angle stability](HOPF_ERROR_ACCUMULATION.md) · [Native depth convention](T_DEPTH_COMPILER.md)

This is a different local construction from geometric-source amplification.
A prepared unary phase state is an eigenvector of every cyclic wire shift.
A bilinear circuit implements a program-selected shift in constant T-depth,
using conditionally initialized work. The same source serves every height
of a group and is unprepared only at its boundary. Its preparation error
is charged twice for the entire group.

Sections 1–6 prove the local interface. Section 7 supplies the complete-frame
precision, group sizes, early/tail split, and query schedule, giving
fixed-accuracy T-depth $`O_\eta(n)`$ with optimal-order T-count and
square-root-scale dirty width. This construction does not replace the
previous conjugated geometric reflection by a shallow implementation.

The [blocked bilinear extension](BLOCKED_BILINEAR_LOOKUP.md) retains this
source construction and proves fixed-accuracy depth $`O(N/b^2+n)`$ with
optimal-order T-count throughout the original sufficient dirty-width
range. The square-root-width theorem below remains valid.

## 1. Local contract

Use the group registers of [the grouped interface](GROUPED_PROGRAM_PREFETCH.md#1-group-program-and-local-statement):
a preserved external prefix x, local bits $`t_0,\ldots,t_{g-1}`$,
and an outer suffix of length r. The external flag h records that this
outer suffix was zero. Put

```math
q=2^\ell,\qquad \ell\ge1,\qquad 1\le g\le\ell,\qquad
R=3^\ell,\qquad 0\lt\delta\le1/4,
\qquad \tau=1+\log_2((\ell+1)/\delta).
```

For each local height d and prefix $`p\in\{0,1\}^d`$, the loaded
program contains a q-bit one-hot word $`e_{a_{d,p}(x)}`$, where
$`a_{d,p}(x)\in\mathbb Z_q`$. Thus the program width is
$`w=q(2^g-1)`$. It is loaded by an exact arbitrary-input XOR lookup
with output $`hF(x)`$, using returned query helpers disjoint from the
conditional suffix. Loading and its actual inverse are separately priced.

Let $`G_a`$ be the exact local Hopf group with angles
$`\theta_{d,p}=2\pi a_{d,p}/q`$. There is an absolute C such that

```math
r\ge C\bigl(q2^g+R+q\log_2 q\bigr)
```

is sufficient to implement this group with initialized-isometry error
at most $`2\delta`$, including return of the source and the outer flag.
Its additional resources, beyond the lookup pair and outer-predicate pair,
are

```math
T=O\!\left(q2^g+gR+q+\ell\tau\right),
```

```math
G=O\!\left(q2^g+gqR+q\ell+\ell\tau\right),\qquad
D_T=O(g+\ell+\tau).
```

Here G counts elementary Clifford gates; long Clifford networks are not
assigned constant elementary depth. Two external clean flags suffice;
this local construction only uses h. Every other initialized input is
explicitly reserved in the active zero suffix. The construction has exact
identity on the full h-zero sector, including arbitrary suffix helpers
and their reference entanglement. There are no intermediate measurements,
resets, supplied phase states, or uncharged reflections.

## 2. A guarded bilinear XOR in sixteen T layers

Consider a fixed bilinear map over $`\mathbb F_2`$ with a rank-R
decomposition

```math
B(p,z)=\bigoplus_{\rho=1}^{R}
\gamma_\rho\,\alpha_\rho(p)\beta_\rho(z),
```

where the alpha and beta functions are linear forms and each
$`\gamma_\rho`$ is an output word. Inputs p and z and output y are
disjoint. Define $`\Gamma_\rho(y)=\gamma_\rho\cdot y`$ over
$`\mathbb F_2`$.

Apply Hadamards to every y wire. For each rho, use three private zero
leaves and CNOTs to compute
$`\alpha_\rho(p),\beta_\rho(z),\Gamma_\rho(y)`$.
With two more private zero wires compute their three-way AND:
first the AND of the first two leaves, then its AND with the third.
All first-round Toffolis have disjoint triples, as do all second-round
Toffolis. Apply $`\mathrm{CZ}(h,q_\rho)`$ to every final AND wire.
Reverse both AND rounds and all linear-form computations, then reverse
the output Hadamards. Every inverse is the actual gate inverse.

On h equal to one and zero private work this has the Fourier-basis phase

```math
(-1)^{\sum_\rho\alpha_\rho(p)\beta_\rho(z)\Gamma_\rho(y)}
=(-1)^{B(p,z)\cdot y}.
```

It is therefore exactly
$`(p,z,y)\mapsto(p,z,y\oplus B(p,z))`$, with all private work
returned. This is a basis identity and extends to arbitrary quantum
inputs and references, including arbitrary y.

On h equal to zero, the central CZ gates are all identity on the entire
Hilbert space. Everything around them cancels as its actual inverse.
This remains true when the private leaves and AND targets are arbitrary;
no copied control is interpreted as clean on that sector. In particular,
the original h is the only control of the central gates.

The circuit has exactly four Toffolis per rank term. Using the literal
seven-T, four-layer native word from the existing
[batch lemma](T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers),

```math
w_{\rm private}=5R,\qquad T\le28R,\qquad D_T\le16.
```

For q-bit operands and output, direct CNOT evaluation of every linear
form gives the conservative Clifford bound $`G=O(qR)`$. The private
copies ensure that no simultaneous native Toffolis share controls.

## 3. Karatsuba convolution and an in-place cyclic shift

Represent q-bit words as polynomials of degree below q over
$`\mathbb F_2`$, and let B be multiplication modulo $`X^q-1`$.
Its output coefficients are cyclic convolution. A self-contained
rank bound is $`R=3^\ell`$: split both polynomials into lower and
upper halves,

```math
A=A_0+X^{q/2}A_1,\qquad B=B_0+X^{q/2}B_1.
```

Compute the three half-size products
$`P_0=A_0B_0`$, $`P_2=A_1B_1`$, and
$`P_1=(A_0+A_1)(B_0+B_1)`$. Their reconstruction is

```math
AB=P_0+X^{q/2}(P_1+P_0+P_2)+X^qP_2.
```

All additions are linear. Recurse to scalar products and finally reduce
the polynomial modulo $`X^q-1`$, also linearly. Hence every scalar
product is a product of one linear form in each input, and the preceding
bilinear decomposition has $`3^\ell`$ terms. No field extension or
division is needed.

For the one-hot program $`p=e_a`$, convolution is the wire permutation

```math
B(e_a,z)=S_a z,\qquad (S_a z)_j=z_{j-a\bmod q}.
```

Let $`p^-_j=p_{-j\bmod q}`$, obtained by merely relabeling the input
wires. Reserve a zero q-bit temporary word y. Apply, chronologically,

1. The guarded bilinear XOR $`y\leftarrow y\oplus B(p,z)`$.
2. The guarded bilinear XOR $`z\leftarrow z\oplus B(p^-,y)`$.
3. Swap each pair $`(z_j,y_j)`$ controlled by the original h.

On h equal to one, one-hot p, and y equal to zero, the intermediate
words are $`(z,S_a z)`$ and $`(0,S_a z)`$, followed by
$`(S_a z,0)`$. The identity holds for every computational z, so the
same circuit permutes arbitrary quantum source inputs and returns y
exactly. It also holds coherently for superpositions of one-hot programs.

On h equal to zero, each bilinear XOR is exactly identity on arbitrary
work, and every controlled swap is identity. Thus the whole shift has
the required full inactive identity even for arbitrary, non-one-hot p.
The one-hot promise is needed only on the active initialized sector.

The existing shared-control Fredkin batch implements the last step
without clean helpers, with at most four T layers and
$`6q+(q\bmod2)`$ T gates, including its literal phase. Reusing the
same $`5R`$ private pool for the two XORs gives

```math
T_{\rm shift}\le56R+6q+1,\qquad
D_{T,\rm shift}\le36,\qquad
G_{\rm shift}=O(qR).
```

The temporary word y and the source word z are separately counted.

## 4. Prepare and decode one phase source

Write $`\omega_q=e^{2\pi i/q}`$ and define the unary source

```math
|\Phi_q\rangle=q^{-1/2}\sum_{j=0}^{q-1}\omega_q^j|e_j\rangle.
```

The exact cyclic permutation above satisfies the literal eigenvalue
identity

```math
S_a|\Phi_q\rangle=\omega_q^{-a}|\Phi_q\rangle.
```

It suffices to prepare this source to norm error delta once per group.
Use an ell-bit index register. Hadamards followed by the one-qubit
phases $`\operatorname{diag}(1,\omega_q^{2^k})`$ prepare the binary
phase state. Replace each phase, up to its fixed scalar, by its
determinant-one diagonal rotation and an actual Clifford+T approximation
of error at most $`\delta/\ell`$. The repository's existing
[single-qubit word-length premise](SOURCE_MAP.md#5-fault-tolerant-sources-and-contribution-boundaries)
supplies length $`O(\tau)`$ for each of these ell words. They act on
disjoint qubits and run in parallel. This imports a word-length bound,
not a new efficient classical synthesis-runtime guarantee.

Decode the binary index as follows. Compute a balanced binary prefix
tree of its literal predicates into private zero wires, using private
control copies at each level. Copy the q leaf predicates into the zero
unary source word. Reverse the entire prefix computation while the
binary index is unchanged. Finally erase each binary index bit by the
CNOT parity of the unary positions whose label has that bit set.

For every binary basis label j, this maps
$`|j\rangle|0^q\rangle|0^{\rm work}\rangle`$ to
$`|0^\ell\rangle|e_j\rangle|0^{\rm work}\rangle`$.
The prefix computation has $`O(q)`$ Toffolis, $`O(\log q)`$
T-depth, and $`O(q)`$ initialized work. The final parity erasure has
$`O(q\ell)`$ CNOTs. All of these maps are exact circuits.

Call the resulting approximate native preparation U. Its initialized
state lies in the unary subspace exactly, even if its synthesized
single-qubit words are not diagonal. By telescoping the ell approximants,
there is a fixed scalar $`\lambda`$, independent of all logical inputs,
such that

```math
\|U|0\rangle-\lambda|\Phi_q\rangle|0^{\rm work}\rangle\|
\le\delta.
```

The boundary inverse will be the actual $`U^\dagger`$. The scalar
lambda cancels in the resulting logical action; it is not discarded as
an input-dependent phase. The preparation and its inverse have

```math
T=O(q+\ell\tau),\qquad
G=O(q\ell+\ell\tau),\qquad
D_T=O(\ell+\tau).
```

## 5. Signed row selection and the literal rotation

At height d retain the exact prefix selectors
$`\pi_p=[t_0\cdots t_{d-1}=p]`$ and enable
$`u_d=[t_{d+1}\cdots t_{g-1}=0]`$. Use the
[incremental selector and consume-before-change enable schedule](GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables),
which has $`O(g)`$ total T-depth and $`O(2^g+g)`$ count and width.
Its validity requires only that a completed stage changes its current
target and preserves every other local bit, including on source leakage.
The stage below has precisely that property.

Put $`D=SH`$, so $`DZD^\dagger=Y`$. Apply $`D^\dagger`$ to
the current target, and call its computational value z during this
stage. In a zero q-bit selected-program word p compute

```math
p_j=[j=0]\oplus u_d[j=0]\oplus
\bigoplus_{v\in\{0,1\}^d}\pi_v u_d
\bigl((1-z)W_{d,v,j}\oplus zW_{d,v,-j}\bigr).
```

On the active sector, this is exactly $`e_0`$ if the inner enable
is false, and $`e_{a_{d,v}}`$ or $`e_{-a_{d,v}}`$ for the selected
row when z is zero or one. Each degree-four product is computed with
four private control copies and three private AND targets; its output
is CNOTed into p, then the computation is reversed. All terms run in
parallel within each of their three Toffoli rounds. Open controls use
literal X conjugations on private copies. The first two terms use an X
and a CNOT on $`p_0`$.

This selected-word computation and its actual inverse have constant
T-depth, $`O(q2^d)`$ T/Clifford count, and $`O(q2^d)`$ private
conditional-zero work. Their Clifford fanout gates may share controls
or targets; all simultaneous non-Clifford gates use private triples.

Apply the guarded cyclic shift from Section 3 to the common source, then
reverse the selected-word computation, and apply D to the target. The
shift preserves every selected-word control and p itself. Consequently
p and its private work erase exactly, including on an arbitrary source
state. On the ideal phase source, the target diagonal in the z basis is

```math
\operatorname{diag}(\omega_q^{-a},\omega_q^{a})
=e^{-i(2\pi a/q)Z}.
```

The complete target word is therefore exactly
$`D e^{-i(2\pi a/q)Z}D^\dagger=R_y(2\pi a/q)`$.
There is no half-angle factor or additional branch phase. An inner
disabled stage selects the identity shift, so its two target basis
changes cancel on every source input in the active zero-work sector.

On h equal to zero the guarded shift is identity on its entire input
space, including arbitrary selected programs and private work. The
selected-word compute/inverse and target basis changes then cancel
exactly. Thus every completed stage has full inactive identity without
assuming clean selector copies there.

## 6. Group return, error, and resources

Load the program, prepare the source with U, and perform the g stages
with the retained-selector schedule. Erase the selectors and source
with their stated actual inverses, unload the original program with
the original lookup inverse, and finally reverse the outer predicate.
The original program and external prefix are never modified.

On h equal to one, every stage is exact on the ideal phase source,
returns its selected program and shift work, and preserves that source.
This holds simultaneously for all logical columns. Let V denote the
exact interior group circuit, including selector computation and cleanup.
Let $`J_F`$ append the valid program $`F(x)`$ and the zero interior
work. Since the group preserves x, this programmed subspace is invariant.
Its ideal-source identity is

```math
V\bigl(|\Phi_q\rangle\otimes J_F|v\rangle\bigr)
=|\Phi_q\rangle\otimes J_FG_a|v\rangle.
```

All stages are unitary on rejected source components; they are neither
discarded nor separately approximated. Comparing the actual prepared
state to $`\lambda\Phi_q`$ before V costs delta. Comparing its
unpreparation to zero after V costs another delta, because the actual
inverse is used. On this valid-program embedding the bound is

```math
\|(U^\dagger VU)J_F-J_FG_a\|\le2\delta,
```

where the source-zero input is implicit in $`J_F`$ in this display.
The original lookup inverse erases $`F(x)`$ exactly and converts this
to the same bound with every program, source, and helper input zero.
The bound is uniform over logical inputs and arbitrary external
references. In particular it is not $`2g\delta`$.

On h equal to zero, every completed stage is exactly identity. The
selector/enable words cancel by the same preserved-control argument as
the incremental schedule, and U cancels with its actual inverse on
arbitrary source and preparation work. Loading the zero inactive row is
identity as well. Thus the inactive identity is exact. The final actual
outer-predicate inverse preserves the combined active/inactive norm
bound, including approximate return of h.

The program occupies $`q(2^g-1)`$ wires; the largest selected-word
computation uses $`O(q2^{g-1})`$ private wires. Source, temporary shift
word, selected program, and source preparation use $`O(q\ell)`$
wires conservatively. The convolution pool occupies $`5R`$ and is
reused after each completed XOR. Selectors and enables use
$`O(2^g+g)`$ more wires. Their disjoint union fits the Section 1
reservation for a fixed absolute C; no suffix wire is simultaneously
counted as a dirty query helper.

Summing selected-word costs gives $`O(q2^g)`$. The g shifts contribute
$`O(gR)`$ T gates and $`O(gqR)`$ Clifford gates at $`O(g)`$
T-depth. The boundary source and retained selectors supply the remaining
terms in the local ledger. The one-hot table queries, outer predicate,
and any angle-grid approximation retain their separate costs.

## 7. Complete-frame theorem at fixed accuracy

**Theorem.** Fix $`0\lt\eta\le1/64`$. There is a constant
$`C_\eta`$ such that, for every $`n\ge1`$, $`N=2^n`$, and
$`b\ge C_\eta\sqrt N`$, every prescribed complete real Hopf frame W
has one coherent Clifford+T circuit V using two external clean flags and
at most b arbitrary dirty qubits, with

```math
\|VJ_2-J_2(W\otimes I_b)\|\le\eta,
\qquad
T=O_\eta(\sqrt N),\quad G=O_\eta(N),\quad D_T=O_\eta(n).
```

The error is the full initialized-isometry norm, including returned
work and arbitrary references. The circuit has no measurements, resets,
supplied phase states, QRAM, or uncharged quantum oracles. Its phase state
is prepared and unprepared by charged native words. As elsewhere in this
repository, T-depth permits arbitrary intervening Clifford circuits;
their elementary gate count is included in G, and no matching bound on
total elementary depth is claimed.

The existing worst-case T-count lower bound is
$`\Omega_\eta(\sqrt N)`$, so the count has optimal order. The available
unrestricted depth lower bound at this width remains $`\Omega(1)`$.
This theorem improves the depth upper bound, not depth optimality or the
high-precision endpoint. Its constants may depend on the fixed accuracy;
it does not replace the all-precision, all-width resource theorem.

### Early groups and their conditional suffix reservation

Write $`\rho=\log_2 3`$. Choose q to be the least power of two with

```math
q\ge \frac{4\pi\sqrt{2n}}\eta,
\qquad
\delta=\frac\eta{8n}.
```

Thus $`q=\Theta_\eta(\sqrt n)`$. Let A be an absolute sufficient
constant from the local unary-group construction, so its simultaneously
live conditional suffix registers fit within

```math
A\left(q2^g+q^\rho+q\log_2 q+g\right)
```

bits. This includes the program, phase source, multiplication buffer,
all bilinear work and copies, and selector/enable maintenance. Increase
A to at least one if necessary, and set

```math
K=\left\lceil16A\bigl(q^\rho+q\log_2q+\log_2q+1\bigr)\right\rceil,
\qquad C=16A.
```

This cutoff directly reserves the fixed source and convolution pools.
It satisfies

```math
K=O_\eta(n^{\rho/2}+\sqrt n\log(n+2)),
\qquad K\log(n+2)=o_\eta(n).
```

For sufficiently large n, depending only on eta and the fixed circuit
constants, $`K\lt n`$. At a group start
$`d=n-k`$, while $`k>K`$, take

```math
g=\left\lfloor\log_2\frac{k}{Cq}\right\rfloor,
\qquad
w=q(2^g-1),
\qquad
r=k-g.
```

Each row stores the q-bit one-hot translation program; w is the full
program width for the group. Since $`q^2\ge n`$ and $`C\ge1`$,
all early groups have $`g\le\log_2q`$, as required by the local lemma.
The cutoff also gives

```math
g\ge\left\lfloor(\rho-1)\log_2q\right\rfloor
=\Omega_\eta(\log(n+2)).
```

For sufficiently large n this is positive. Moreover,

```math
Aq2^g\le k/16,
\qquad
A(q^\rho+q\log_2q+g)\le k/16,
\qquad g\le k/16.
```

These inequalities follow directly from the group-size choice and the
cutoff, using $`g\le\log_2q`$. The sufficient reservation is at most
$`k/8`$, while $`r=k-g\ge15k/16`$. These are logical suffix bits, zero only on the group's active
sector; they are not counted again as arbitrary dirty query work. The
external dirty query pools and two returned outer-predicate helpers are
reserved separately. All inactive-sector statements use the local
construction's literal guarded identities on arbitrary work.

There are $`O_\eta(n/\log(n+2))`$ early groups. The last group may
cross the cutoff by fewer than g layers; its remaining height
$`K'\le K`$ is the start of the individual-layer tail. No rounded early
angle is used by that tail.

### Full error, including the phase-source return

Round every early angle to its nearest $`2\pi/q`$ grid point, choosing
a real representative of the difference of magnitude at most
$`\pi/q`$. Leave all other ideal angles unchanged. The
[complete-frame angular bound](HOPF_ERROR_ACCUMULATION.md#1-a-sharp-angle-error-bound-on-the-full-frame)
gives

```math
\|W_{\rm grid}-W\|
\le\frac{\pi\sqrt{2n}}q\le\eta/4.
```

With the exact unary Fourier eigenstate, every early group is exactly
its prescribed grid-angle group and returns that state. On the active
sector, let $`J_F`$ append the valid loaded program $`F(x)`$ and zero
source, preparation, selector, and shift work. Let P be an ideal source
preparation, with its fixed scalar chosen as in Section 4, and
$`\widetilde P`$ its native replacement. Both act only on source and
preparation work, so Section 4 gives
$`\|\widetilde PJ_F-PJ_F\|\le\delta`$ uniformly over logical
inputs. If U is the exact guarded middle word of the complete group,
including selector computation and cleanup, then

```math
UPJ_F=PJ_FW_{\rm group}.
```

Consequently the charged prepare-use-unprepare word obeys

```math
\|\widetilde P^\dagger U\widetilde PJ_F-J_FW_{\rm group}\|
\le2\delta.
```

The first delta bounds the input-preparation difference. For the second,
use the exact return identity and

```math
(\widetilde P^\dagger-P^\dagger)PJ_F
=\widetilde P^\dagger(P-\widetilde P)J_F.
```

The grid-angle group preserves the external prefix x, so its action
also preserves the valid-program embedding $`J_F`$.
Thus this estimate requires neither an operator approximation on arbitrary
source inputs nor a fresh source at each height. It includes the final
source and preparation-work leakage. Exact program and selector return
remain valid through the middle word. On the outer inactive sector the
whole completed group is exactly identity even on arbitrary suffix work,
by the local guarded-cancellation proof. The original lookup inverse
erases $`F(x)`$ exactly, converting the active estimate to the contract
with zero program input and output. Exact load/unload and the actual
outer-predicate inverse preserve the norm estimate.

Prepare the binary Fourier product state with its
$`\log_2q`$ one-qubit factors approximated to error at most
$`\delta/\log_2q`$, then use the exact unary conversion. The established
single-qubit native approximation bound permits these factors to run on
disjoint wires in parallel. Their T-depth is
$`O(\log(\log(q)/\delta))=O_\eta(\log(n+2))`$; their count and the
charged exact conversion fit the local ledger. The actual circuit
inverse is used on exit. Any consistent global phase in the prepared
state cancels against that actual inverse.

There are at most n early groups, so unitary telescoping charges at most
$`2n\delta\le\eta/4`$ for their native phase preparations. For the
remaining layers, use the original capped additive source certificate
with accuracy parameter eta divided by four. Restricting that
nonnegative error sum to the tail costs at most $`\eta/4`$. Thus the
total error is at most $`3\eta/4\le\eta`$. Telescoping compares each
stage on its ideal initialized input and uses unitarity on preceding
leakage; it never resets actual work between groups. The estimate is
therefore a full-frame initialized-isometry bound, not a state-preparation
or accepted-block-only estimate.

### Early prefetches, native counts, and depth

Use the existing parallel single-output dirty-counter bilinear queries
to prefetch the w program bits from the $`d+1`$ address bits
$`(h,x)`$. Their dirty pools are separate, and the shared address enters
through Clifford controls. With $`Q=2^{d+1}`$,

```math
T_{\rm load},w_{\rm dirty}
=O\!\left(w\sqrt Q(n+2)^3\right),
\qquad
D_{T,\rm load}=O(\log(n+2)),
```

```math
G_{\rm load}
=O\!\left(wQ+w\sqrt Q(n+2)^3\right).
```

The $`wQ`$ Clifford cost of arbitrary table coordinate changes is
retained. The program outputs are the separately reserved logical
suffix wires. Actual unload has the same costs. Since $`w\le k/C`$
and early starts have distinct remaining heights k,

```math
\sum_{\rm early}T_{\rm load}
\le O\!\left(\sqrt N(n+2)^3
                 \sum_{k>K}k2^{-k/2}\right)
=O_\eta(\sqrt N),
```

```math
\sum_{\rm early}G_{\rm load}
\le O\!\left(N\sum_{k>K}k2^{-k}\right)+O_\eta(\sqrt N)
=O_\eta(N).
```

The first bound also bounds every peak private dirty pool. The large
cutoff absorbs these polynomial factors with ample exponential margin.
The charged local selection, translation, source preparation, and outer
predicate costs are polynomial in n. The local translation and selection
cost $`O(q2^g+gq^\rho)`$ T gates and the conservative
$`O(q2^g+gq^{\rho+1})`$ Clifford gates. Source preparation, conversion,
and predicates also have polynomial cost. Their sum over at most n
groups therefore fits both stated exponential gate-count bounds. No
per-group linear-in-k count claim is needed at this cutoff.

Each group's internal T-depth is

```math
O(g+\log(q)+\log(\log(q)/\delta))=O_\eta(\log(n+2)).
```

Its prefetch, unload, and outer predicate/inverse also cost
$`O(\log(n+2))`$ depth. Multiplying by the number of groups gives
$`O_\eta(n)`$ total early depth.

### One chunked schedule for the entire remaining tail

Set

```math
L'=\max\{6,\lceil\log_2(4/\eta)\rceil\},
\qquad h_n=\lceil\log_2(8n)\rceil,
\qquad m=L'+4+h_n.
```

On every remaining layer $`1\le k\le K'\le K`$, use the original
full-input operator source and literal amplification at the original
capped width

```math
m_k=L'+4+\min\{k,h_n\}\le m=O_\eta(\log(n+2)).
```

Use the already proved
[chunked-query schedule](CHUNKED_DIRTY_INDICATOR.md#5-a-summable-budget-for-the-final-hopf-layers)
for this entire tail. The coefficient query has
$`Q_k=4N2^{-k}`$ rows and $`n-k+2`$ address bits. Split that full
address into balanced halves and use chunk cap
$`a_k=2^{\lfloor k/12\rfloor}`$, truncated to each half's length.
The zero-bit half is a Clifford X. The established resource bounds give

```math
T_k,w_k
=O\!\left(\sqrt N\,2^{-k/2}[m_k+2^{k/4}]\right),
```

```math
G_k=O\!\left(N2^{-k}m_k+\sqrt N\,2^{-k/4}\right),
```

```math
D_{T,k}
=O\!\left(m_k+n(k+1)2^{-k/12}+\log(n+2)\right).
```

Here the sequential output-bit loop of the bilinear oracle costs
$`O(m_k)`$ depth; the two indicators are shared by those output bits.
The last bound follows directly from the truncated chunk sizes: if a
chunk cap is below its address-half length, its indicator depth is
$`O(n(k+1)2^{-k/12})`$; otherwise it is $`O(\log(n+2))`$.
The bounded number of queries and actual inverses in each amplified
layer changes only constants.

The count and workspace sums converge uniformly for every tail length
at most n, since $`m_k\le L'+4+k`$:

```math
\sum_{k=1}^{K'}T_k=O_\eta(\sqrt N),
\qquad
\sum_{k=1}^{K'}G_k=O_\eta(N),
\qquad
\max_{1\le k\le K'}w_k=O_\eta(\sqrt N).
```

For depth, retain the sharper cap $`m_k\le m`$, rather than the
loose bound $`L'+4+k`$. Because
$`\sum_{k\ge1}(k+1)2^{-k/12}\lt\infty`$,

```math
\sum_{k=1}^{K'}D_{T,k}
=O_\eta\!\left(n+K\log(n+2)\right)=O_\eta(n).
```

The original full-input sources have $`O(m_k)`$ serial T-depth.
Their reflections and the separately priced suffix predicates fit
$`O_\eta(\log(n+2))`$ additional depth per layer, with polynomial
counts. Thus all non-query tail depth is
$`O_\eta(K\log(n+2))=o_\eta(n)`$. This construction retains the
full n-bit external address and every complete-frame marker column; it
does not reinterpret the tail as a smaller independent frame.

Keep the original dirty base $`B_0=L'+n+7`$ and the separately
reserved predicate helpers. Reuse each query pool only after its exact
return. Every peak fits $`C_\eta\sqrt N`$ for a sufficiently large
constant. The early and tail ledgers now give all three bounds on the
same circuit. The finitely many n below the stated asymptotic thresholds
are covered by the previously proved fixed-accuracy compiler, increasing
only the eta-dependent constants.

## 8. Attribution and bounded evidence

Fourier-state phase kickback and reuse of a prepared phase reference are
established ingredients, as are the three-product Karatsuba identity,
parallel parity-phase synthesis with initialized ancillas, reversible
binary-to-unary conversion, and the single-qubit approximation bound.
The [source map](SOURCE_MAP.md#5-fault-tolerant-sources-and-contribution-boundaries)
records these dependencies, including Jones et al., Iggy van Hoof,
Selinger, and the existing native synthesis toolkit. The Karatsuba
reference supplies the multiplication identity; its reversible schedule
is not the source of this chapter's constant T-depth bound.

The additional interface proved here is the coherently loaded one-hot
translation with charged conditional work, literal inactive identity,
and exact temporary return on every active source state. Its combination
with a charged source boundary pair, retained local selectors, Hopf angle
stability, and the early/late query allocation gives the complete-frame
theorem. [Related work, Section 19](RELATED_WORK.md#19-unary-phase-source-reuse-and-linear-t-depth-3-october-2026)
compares the catalyst and shared-arithmetic constructions of
Kim–Laakkonen and Wu et al. Their results remain separate from the
specific coherent table and initialized-isometry contract here. These
comparisons identify inherited ingredients and resource distinctions;
they do not establish priority or a generic constant-depth rotation
theorem without supplied resources.

The [bounded checks](../tests/test_unary_phase_gradient.py) cover exact
small Karatsuba identities and source-bit permutations, a literal native
guarded phase gadget, small native source preparation and its actual
inverse, unequal three-stage row programs, and the retained-source
$`2\delta`$ error bound. The [verification index](VERIFICATION.md)
states their exact ranges and negative controls. The group fixtures use
reduced operators on the invariant unary source space; they do not emit
the complete native bilinear XOR, selected-program network, asymptotic
decoder, or full-frame lookup circuit. Uniform return, workspace, gate
counts, depth, and the global error bound follow from the analytic
construction above, not an extrapolation of those finite checks.
