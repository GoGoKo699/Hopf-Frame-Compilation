# Pauli-axis overlap across time and robust decision bias

[Three-band bias theorem](ROBUST_DECISION_BIAS.md) ·
[Axis relations and components](AXIS_TRAJECTORY_BIAS.md) ·
[Decision interfaces](PRESCRIBED_SCALAR_TABLES.md)

A circuit with square-root T-count that decides a typical prescribed
table at constant bias needs a square-root-dimensional shared Pauli
span across some cut between its actual native bands. The count uses
the spans of the pulled-back axes, with its cut taken in the actual
chronological order. Address motion, Clifford interlayers, dirty width
and the number of bands are unrestricted. A separate conditioning
argument requires both dense anticommutation and a large pairing
rank across an actual temporal cut. The given width allows all
these quantities to be large, and R96 remains the general upper.

## 1. From a native word to an exact cut contraction

Fix m computational inputs and average over all other inputs. Write
$`b(x)`$ for the expectation of a fixed output Pauli at the fixed input
$`x\in\mathbb F_2^m`$ and maximally mixed remaining wires. Use
the normalized Hilbert–Schmidt inner product
$`\langle A,B\rangle=2^{-w}\mathrm{Tr}(A^\dagger B)`$ on the
$`w`$ physical qubits. Pull the actual Clifford walls towards the output
and absorb their product into the measured Pauli. This leaves
$`d`$ chronological bands in a common Pauli coordinate frame. Band
$`j`$ has $`t_j`$ signed rotations and binary axis span $`L_j`$.
Its axes are isotropic and linearly independent because they are the
images of distinct single-qubit Z labels under one common Clifford.
Discard empty bands and put $`T=\sum_j t_j`$.

Write $`\sigma_q=X^{q_X}Z^{q_Z}`$ for the literal Pauli section.
If the input probe has label $`p(z)`$ and the measured Pauli label is
$`s`$, the Walsh coefficient $`\widehat b(z)`$ is the corresponding
matrix coefficient of the adjoint word on the orthonormal
$`\sigma_q`$ basis. A physical computational input probe is
$`Z^z`$. The same argument works for any affine input label $`p(z)`$
whose literal phase modulo four has binary digits of degree at most
two, as occurs when projecting a Clifford-coordinate probe onto a
tensor factor. The output phase is fixed. Pulling out scalar factors
of T and its inverse is exact in the adjoint.

Let $`L_{\le j}=\sum_{i\le j}L_i`$ and $`L_{>j}=\sum_{i>j}L_i`$.
At cut $`j`$ only labels lying in the intersection below can connect
the two endpoints. Write

```math
\begin{aligned}
E_j(z)&=(p(z)+L_{\le j})\cap(s+L_{>j}),\\
B_j&=L_{\le j}\cap L_{>j},&
b_j&=\dim B_j,& b_*&=\max_j b_j .
\end{aligned}
```

All the sets $`E_j(z)`$ are nonempty exactly on the same affine
domain $`A=\{z:p(z)+s\in\sum_i L_i\}`$. On A choose one affine
solution $`u^0(z)`$ of the full endpoint equation in the concatenated
axis lists. Its partial turns give compatible affine origins
$`q_j^0(z)`$ in every $`E_j(z)`$. Fix any basis for each $`B_j`$,
so cut labels are parametrized by $`y_j\in\mathbb F_2^{b_j}`$.
The initial and final cuts have dimension zero. Project the adjoint
word onto these label sets at the cuts. Paths that leave a cut set
can never reach the measured label, by the two span conditions.
Therefore the Fourier coefficient, up to the literal endpoint
phases, is the scalar product of the projected band matrices
$`K_d(z)\cdots K_1(z)`$. Each projected matrix has operator norm at
most one. This projection is an operation in the proof on an
orthonormal operator basis.

**A literal band entry.** If $`q`$ and $`q'`$ are predecessor and
successor labels, the turn vector $`u`$ is uniquely determined by
$`q+q'=\sum_i u_i p_i`$ in the band's independent axis list
$`(p_i)`$. If there is no such vector the entry is zero. Set
$`\alpha_i=\omega(p_i,q)`$ for the binary commutation form. Since
the axes commute within the band, every activation $`\alpha_i`$
is determined by the incoming label. A turn with $`u_i=1`$ is
allowed only when $`\alpha_i=1`$. At any allowed entry the exact
coefficient has absolute value $`2^{-|\alpha|/2}`$. Its phase is
$`i^\Phi`$, with

```math
\Phi=\sum_i c_i u_i
 +2\sum_{i\gt h}(p_{i,Z}\cdot p_{h,X})u_i u_h
 +2\sum_i(p_{i,Z}\cdot q_X)u_i\pmod4,
```

where $`c_i`$ includes the signed Hermitian Pauli convention and
the orientation of its rotation. Indeed an anticommuting
conjugation maps a Pauli to the sum of itself and its left product
by a signed $`iP_i`$, divided by $`\sqrt2`$. Successively multiplying
these left factors gives the displayed ordered cocycle from
$`\sigma_a\sigma_b=(-1)^{a_Z\cdot b_X}\sigma_{a+b}`$.
On the affine set where the endpoint difference is in $`L_j`$,
both the turns and the activations are affine. The two bits of
$`\Phi`$ then have degree at most one and two, respectively.
The linear sum modulo four has only a pairwise binary carry, and
every other term is twice a binary quadratic.

## 2. An exact description at a T-count budget

For caps T and $`b_*`$, the possible bias functions admit an exact
binary description with

```math
O\bigl(m^2+(m+b_*+1)(T+1)\bigr)
```

bits. The constant is independent of the physical qubit count.
To prove this, store the affine input domain A and endpoint
phases, and choose the cut bases successively as below. Counts,
orders and cut dimensions can be recorded within the bound.
In particular a positive composition of the T-count needs only
O(T) bits to delimit bands.

Put $`n_0=m+b_*+1`$. If $`t_j\ge n_0`$, use the full set of triples
of input z and predecessor and successor coordinates. Its affine
subset on which the endpoint difference is in $`L_j`$ takes
$`O(n_0^2)`$ bits. Store on this domain the affine activation and
turn maps, costing $`O(t_j n_0)`$ bits, and the quadratic phase
data, costing $`O(n_0^2)`$ bits. The decoder evaluates all entries
using Section 1 and contracts the matrices exactly. The
$`O(n_0^2)`$ terms fit $`O(t_j n_0)`$.

The small-band case needs a different encoding to keep the same
per-band bound. Suppose $`1\le t_j\lt n_0`$. Define

```math
K=L_{\le j-1}\cap L_{>j}=B_{j-1}\cap B_j.
```

The codimension of K in each neighboring cut span is at most
$`t_j`$. For the predecessor use
$`B_{j-1}\subseteq L_j+L_{>j}`$ and take the quotient by
$`L_{>j}`$. For the successor use
$`B_j\subseteq L_{\le j-1}+L_j`$ and take the quotient by
$`L_{\le j-1}`$. The quotient bounds charge the newly used band,
regardless of the ambient size of either cut span.

Choose a full-rank equation matrix for K in the current input
basis of $`B_{j-1}`$. Its codimension is $`k\le t_j`$ and its
description costs $`O(t_j b_{j-1})`$ bits. A deterministic choice
of pivot columns aligns the input coordinates with K plus at
most $`k`$ residual coordinates. In detail, each nonpivot basis
vector is corrected by the unique linear combination of the
pivot vectors that places it in the kernel. The decoder keeps the
nonpivot coordinates and updates the residual coordinates using
the equation matrix. Choose the outgoing basis to extend this
common K basis, so its contribution passes through unchanged.

For a candidate turn $`u`$ from an incoming label, subtract the
known origin turn $`u^0_j(z)`$. The remaining displacement
relative to the outgoing origin is the common K component plus
a combination of at most $`k`$ residual input directions and
the $`t_j`$ axis directions. Its membership in $`B_j`$ is determined
by reducing that combination modulo $`L_{>j}`$. The condition
uses at most $`t_j`$ equations on those $`k+t_j`$ coordinates:
all the directions already lie in $`L_j+L_{>j}`$ and inside
$`L_{\le j}`$. On the kernel of that condition, the map of the
combination to outgoing cut coordinates is linear. Store any
linear extension of it, using $`O(b_jt_j)`$ bits. The decoder adds
the unchanged K coordinates. The affine origin turn, activations,
signed axes' pairings and the affine factors $`p_{i,Z}\cdot q_X`$
in the phase together cost $`O(t_jn_0)`$ bits. Here
$`t_j^2\le t_jn_0`$. This allows exact enumeration of the allowed
turns with their phases and destination labels. It avoids storing
a dense change of basis on a large common span at a small band.

Summing over bands proves the exact description bound. The decoder
may take exponential classical time. No gate implementation of a
coordinate change or projected band is assumed.

## 3. A width-free approximate description

There is also a normalized $`L^2`$ net at tolerance
$`0\lt\epsilon\le1/2`$ with logarithmic cardinality
$`O(d n_*^2\log(n_*+2))`$, where
$`n_*=m+b_*+1+\log(d/\epsilon)`$ and $`d\ge1`$. There is no
T-count cap in this statement. Treat a Clifford-only word as an
additional trivial case.

Set $`K_0=\lceil 2b_*+m+2\log_2(4d/\epsilon)\rceil`$. For a
band retain valid entries with at most $`K_0`$ activations. Store
an affine basis chosen from the retained triples and the affine
hull it defines. The hull stays in the shift-matching domain.
On each of its at most $`m+b_{j-1}+b_j+1`$ basis points list the
nonzero activation and turn coordinates in a relabeled common
universe. Each list has at most $`K_0`$ entries, and the union
has size at most $`(m+b_{j-1}+b_j+1)K_0`$. These sparse lists
determine the full affine maps on the hull. The phase data are
quadratic there by Section 1. The hull, lists and phases therefore
fit $`O(n_*^2\log(n_*+2))`$ bits per band.

Evaluate the exact valid entries within the hull and set the rest
to zero. The additional hull points may have more than $`K_0`$
activations, which the decoder evaluates without approximation.
Invalid turns in the hull still yield zero. Every omitted valid
entry has magnitude at most $`2^{-K_0/2}`$. For each input the
Frobenius, and hence operator, error of the band is at most

```math
2^{(b_{j-1}+b_j-K_0)/2}\le
 \frac{\epsilon}{4d}\,2^{-m/2}.
```

Expanding the product around the exact matrices bounds its
error by the product of one plus the individual errors, minus
one. This is at most $`\epsilon 2^{-m/2}`$, including the
possible increases of the approximate matrix norms. Parseval
gives the normalized $`L^2`$ bound on b, taking the real part
of the approximate reconstruction.

**Random tables.** For any one fixed exact bias on the a-bit
address after setting the other $`m-a`$ bits, the Rademacher tail
for its sign-table correlation at threshold $`\gamma`$ is at
most $`2\exp(-2^a\gamma^2/2)`$ because the bias is bounded
by one. Under $`T\le C2^{a/2}`$, $`m=\mathrm{poly}(a)`$ and
fixed positive $`\gamma`$, the exact code and a union bound
imply that typical tables need $`b_*=\Omega_{C,\gamma}(2^{a/2})`$.
The bound holds for each successful native decomposition under
the count cap, with any number of bands. Both sides of a cut
attaining this bound must contain square-root-scale T-count.
For strict pointwise sign, regardless of the margin size, each
exact description realizes at most one Boolean table. This gives
the same order with constants independent of the margin.
For the width-free net, request ambient accuracy
$`\min(1/2,2^{-m/2})`$ so restriction to the configured input bits still
gives error at most $`2^{-a/2}`$. The same counting argument
requires $`d n_*^2\log(n_*+2)=\Omega_\gamma(2^a)`$ for typical
constant-gap tables. These results concern quantitative bias and
also apply if the stronger uniform operator contract implies it.

## 4. Conditioning a commuting radical

There is a separate budgeted screen for a signed gap on every
input of the specified address and clean code. It counts true
noncommuting Pauli pairs even when the total axis span has a
large commuting radical. Include in a binary span L the pulled-back
T axes, the fixed input Z generators and the output Pauli label.
Let $`2r`$ be the rank of the symplectic form on L.

Choose one Clifford coordinate change that puts L on $`r`$
active symplectic pairs and commuting Z labels on separate
radical wires. The simultaneous Pauli coordinate change is
analytical. The fixed input Z generators remain independent
and commute. On the kernel of their projection to the active
pairs they have purely radical labels, and those labels are
independent for an independent kernel basis. Their literal signs
on the kernel form a character because they are commuting Z
products in the new frame. Thus for every prescribed address and
clean bits the kernel constraints on the radical computational
label are consistent. Choose an affine solution $`y(x)`$ across
all addresses. After that choice the other fixed input constraints
give commuting independent active Paulis with affine signed
eigenvalues. Use the maximally mixed state on the common eigenspace
and on any spectators, with the chosen radical eigenlabel.

On the chosen inputs every T-axis orientation and the measured
Pauli sign depend affinely on x. All noncommuting evolution occurs
on $`r`$ active qubits. The bias for this selected input family
is described by the active Pauli labels, the affine orientation
bits, the input active stabilizers with their affine signs and
the measured Pauli. This costs

```math
O\bigl((T+m+1)(r+a+1)\bigr)
```

bits. The selected family is permitted to depend on the circuit.
The full diagonal operator gap on each address holds on it,
so counting these restricted functions still forces
$`r=\Omega_C(2^{a/2})`$ for typical tables at square-root
T-count and polynomially many fixed bits. Removing the input
generators and measured label from L reduces its number of
symplectic pairs by at most $`m+1`$. The same order therefore
holds for the span of the T axes alone under this contract.

There is a variant for a margin promised on each computational
dirty word. Enlarge L by all computational Z labels. Its radical
is then computational Z-diagonal in the original input frame,
since it commutes with every input Z. Conditioning on the
chosen radical label and fixed input bits gives a diagonal mixed
state in that original frame, so the pointwise promise still
applies. The enlarged span's number of symplectic pairs is
exactly the rank of the computational X masks of the T axes
and measured Pauli. This variant constrains that mask rank.
The smaller symplectic span statement uses the full operator
gap, which allows testing coherent dirty inputs.

### The anticommutation graph and its temporal cuts

The active labels of the whole ordered list of T axes, input
generators and output Pauli can instead be recovered from the
binary Gram matrix of their pairings. The map from formal list
coefficients to the nondegenerate quotient of L has kernel exactly
the Gram matrix's kernel. The quotient with its induced form
therefore recovers the listed labels up to a common symplectic
isomorphism. Fix deterministic Pauli coordinates for this
quotient. The Pauli signs from that change of coordinates
are fixed and are included in the stored affine sign bits. The
decoder can evaluate the chosen-input trace exactly. Pairs involving
input and measurement labels are stored entrywise.

Use chronological order for the T occurrences. If E unordered
pairs of them anticommute, their Gram matrix takes
$`O(\log(T+2)+E\log(e(T+1)^2/(E+1)))`$ bits, by describing the
set of E pairs. The remaining data take
$`O((T+m+1)(m+a+1))`$ bits. With $`T\le C2^{a/2}`$,
polynomially many fixed inputs and a strict coherent sign at every
address, a typical table requires $`E=\Omega_C(2^a)`$.
Each exact description realizes at most one strict sign table.
The sparse-graph codes with $`E\le c\,2^a`$ take fewer than
$`2^{a-1}`$ bits for sufficiently small c depending on C. This
forces square-root-scale degrees on a square-root-scale set of
occurrences. Each occurrence is counted at its actual position.

There is also a code in terms of the maximum prefix-suffix
pairing rank across actual band cuts. Let $`G^{(T)}`$ be the
T-axis portion of the Gram matrix. Write $`C_j`$ for the row
space of its block from bands at most j to bands greater than j.
Put $`c_j=\dim C_j`$ and $`c_*=\max_j c_j`$. The Gram matrix
then takes $`O(T(c_*+1))`$ bits. Keep coordinates for the
prefix rows in a basis of $`C_{j-1}`$. At band j store the
linear map that reads those basis vectors on the new
$`t_j`$ columns. This recovers all pairings involving those
columns and earlier rows. Deleting the columns projects to
a subspace of $`C_j`$. The kernel has dimension at most $`t_j`$
and is exactly the subspace supported on the deleted columns.
Store a kernel basis in the old coordinates in $`O(c_*t_j)`$ bits.
Extend it by coordinate vectors chosen by a fixed rule and use
the projected complement as the initial part of the new basis.
The new band's rows span the remaining quotient, of dimension
at most $`t_j`$. Store their full coordinates in the extended
basis in $`O(c_*t_j)`$ bits. The decoder projects each old row
by the encoded kernel and appends the new rows. Pairings inside
the band vanish because its axes commute. The endpoint row
spaces are zero and the code telescopes with the band counts.

The exact sign count now also forces $`c_*=\Omega_C(2^{a/2})`$
for typical coherent tables. Set intersection and pairing rank
are different conditions. For instance, distinct independent
Z axes and X axes can have zero shared span with full pairing
rank, while two identical Z bands have the reverse behavior.
The two requirements can select different cuts.

One large orthogonal component can meet all the rank and
density scales in a four-band identity. For even $`q\ge4`$, take
A on the single-qubit Z axes and B on X axes whose masks are
the columns of $`I+J_q`$, with $`J_q`$ the all-ones matrix.
This matrix squares to identity over $`\mathbb F_2`$. Its
incidence graph is connected and dense. Each independent
commuting list occurs in one actual band after a Clifford
basis change. Follow A and B by their actual inverses in
reverse order. At the middle cut the two sides have the same
span of dimension $`2q`$ and their pairing has rank $`2q`$.
The internal nullity is $`2q`$ and the symplectic pair count is q.

The fixed-error corner-frame readout from the
[three-band proof](ROBUST_DECISION_BIAS.md#a-direct-complete-frame-error-floor)
converts a complete frame approximation below its threshold to
such a uniform constant gap by measuring the output on the
fixed last input and clean sector. These bounds allow arbitrary
address motion and work return in that conversion. The remaining
depth question is how to control relative interference across
large shared spans and dense noncommuting couplings without
counting exactly reversed operations as useful bias.

<a id="correlated-packet-prefix-obstruction"></a>

## 5. Correlated packets cannot be compressed separately at a cut

A native band has an exact packet update, but a net for each side of
a cut in a fixed physical Pauli frame can already require
quadratic-in-T description size.
The obstruction below holds for two-band prefixes that extend to exact
full return. It does not exclude a joint description conditioned on the
completed word, or prove a growing frame-depth lower bound.

### The exact one-band update retains shared histories

Write a Pauli packet as a sum with labels $`\ell(x,z)`$ affine in the
probe parameters x and its retained history z. For a new physical band
with independent commuting axes $`a_1,\ldots,a_h`$, put
$`\alpha_i=\omega(a_i,\ell(x,z))`$. Each old term branches over
$`v\in\mathbb F_2^h`$ with $`v_i\le\alpha_i`$, receiving magnitude
$`2^{-|\alpha|/2}`$ and label

```math
 \ell(x,z)+\sum_i v_i a_i .
```

The literal phase is the quadratic modulo-four Pauli cocycle from
Section 1, including the actual T orientations. Equal final labels
must be summed coherently before any significance test. The activation
map contains every cross pairing of the new axes with the old history
directions. Section 2, with $`b_*\le T`$, already gives an exact
description of the completed probe contraction of size
$`O(T^2+mT+m^2)`$, independently of physical dirty width. This joint
description in common Pauli coordinates does not identify that width
with the number m of configured input bits. The following packing
shows why a uniformly accurate prefix-operator net in the fixed
physical frame cannot generally remove the quadratic term.

### A two-band packing within the simultaneous budgets

Use two q-qubit banks, $`q\ge1`$. For each
$`M\in\mathrm{GL}(q,2)`$ take the native axes

```math
 A_i=Z_i,\qquad B_j=X^{Me_j}\otimes Z_j,\qquad
 P=X^{\boldsymbol1}\otimes I .
```

The first-bank axes A and the axes B each form an independent commuting
list. Let $`E_A`$ and $`E_B(M)`$ be their literal native T bands, and
put $`U_M=E_B(M)E_A`$ and $`O_M=U_MPU_M^\dagger`$.
Conjugating the second-bank T targets by first-to-second CNOTs specified
by M, then first-bank Hadamards, realizes B. Thus $`U_M`$ has 2q T gates,
two native bands and $`O(q^2)`$ elementary Cliffords on 2q wires.

In the literal section $`\sigma(x,z)=X^xZ^z`$, its full Pauli expansion is

```math
 O_M=\sum_{u,v\in\mathbb F_2^q}
  2^{-(q+|M^{\mathsf T}u|)/2}\,
  i^{|u|-|v|}\,\mathbf1_{v\le M^{\mathsf T}u}\,
  \sigma(\boldsymbol1+Mv,u)\otimes Z^v .
```

The inequality in the indicator is coordinatewise. The first band
contributes $`i^{|u|}2^{-q/2}`$. The second band has activation
$`M^{\mathsf T}u`$ and contributes $`(-i)^{|v|}`$; its X-only
first-bank factors introduce no additional section phase. Every label
in this expansion has a unique pair (u,v).

Even the first band alone gives $`2^q`$ coefficients of magnitude
$`2^{-q/2}`$, although its overlap with the identical backward packet
is one. Thus coefficientwise significance can discard all of an
order-one cut contraction.

The squared coefficient distribution has u uniform and, conditional on
u, v uniform on the coordinates where $`M^{\mathsf T}u`$ is one.
Since M is invertible, the marginal bits of v are independent Bernoulli
$`1/4`$. The supports for M and $`M'`$ can coincide only where
$`(M-M')v=0`$. A rank-r binary system accepts an independent Bernoulli
$`1/4`$ vector with probability at most $`(3/4)^r`$: condition on
nonpivot coordinates and use the largest single-bit probability at
each pivot. Cauchy--Schwarz on the shared support gives

```math
 |\tau(O_M O_{M'})|\le(3/4)^{\mathrm{rank}(M-M')}.
```

Here $`\tau=2^{-2q}\mathrm{Tr}`$ and every $`O_M`$ has normalized
Hilbert--Schmidt norm one. Distinct matrices therefore give distance
at least $`1/\sqrt2`$. For $`\epsilon\lt1/(2\sqrt2)`$, any epsilon-net of
these prefix packets needs at least
$`|\mathrm{GL}(q,2)|=2^{q^2}\prod_{j=1}^q(1-2^{-j})\ge2^{q^2-2}`$
centers. The last bound follows by keeping the first factor and using
$`\prod_{j=2}^q(1-2^{-j})\ge1-\sum_{j=2}^q2^{-j}\ge1/2`$.
Likewise, a prefix synopsis
that predicts its contraction with every two-band packet to error below
$`1/8`$ must distinguish all M: the test against $`O_M`$ gives one for
M and at most $`3/4`$ for any other matrix.

A fixed CNOT fanout F satisfies $`FX_1F^\dagger=P`$, so the prefixes
$`U_MF`$ transport one logical X to these same packets. Their actual
reverses $`F^\dagger U_M^\dagger`$ complete identity on every logical,
dirty and reference input; one clean spectator may be included. With
$`q=\Theta(\sqrt N)`$, this four-band completed word has
$`T,b_{\rm dirty}=O(\sqrt N)`$ and $`G=O(N)`$. The prefixes are not
required to return work. Thus the packing applies even to prefixes
extendable to a circuit in the full permitted resource class. It
concerns a separately reusable prefix or suffix description, not a
description allowed to use its matched partner and terminal return.

### A relative matrix removes the common large description

For this same family define
$`K=(M')^{\mathsf T}M^{-\mathsf T}-I`$. At a matched label the literal
phases cancel. Changing variables to $`\alpha=M^{\mathsf T}u`$ gives
the exact positive overlap

```math
 \tau(O_M O_{M'})=2^{-q}\sum_{\alpha\in\mathbb F_2^q}
  2^{-(|\alpha|+|(I+K)\alpha|)/2}
  \left|\ker K^{\mathsf T}\cap
  \mathbb F_2^{\mathrm{supp}\alpha\cap
    \mathrm{supp}((I+K)\alpha)}\right| .
```

The coordinate space in this expression consists of vectors supported
on the indicated intersection. The formula depends on K alone. If the
overlap is at least $`\delta>0`$, the preceding bound forces
$`\mathrm{rank}K\le\log(1/\delta)/\log(4/3)`$. A rank
factorization then describes K in $`O(q(1+\log(1/\delta)))`$ bits,
instead of recording M and $`M'`$ separately.

This is a transparent special case of the relative-comparison strategy
in [the paired-half analysis](RELATIVE_FOUR_BAND_BIAS.md), not an
extension of its general inverse theorem. It describes this scalar
contraction. Closure of a comparably short joint description after an
additional arbitrary band, including coherent label collisions, remains
unproved. The [bounded native fixture](../../tests/test_native_layer_count_distribution.py)
compares the literal expansion and relative formula with all full-space
native words for $`\mathrm{GL}(2,2)`$ and their actual reverses. The
packing and all-size resource statements above are analytic.

<a id="actual-band-relative-contraction"></a>

## 6. Contract an inserted band with its actual suffix

The preceding relative matrix describes one completed scalar contraction.
It is not a state that can answer arbitrary later insertions. For a word
$`W=REL`$, retain the actual suffix and write
$`L_z=LP_zL^\dagger`$ and $`R_S=R^\dagger SR`$. If their literal Pauli
histories are $`c_L(z,u),\ell(z,u)`$ and $`c_R(v),r(v)`$, an inserted
band E with independent commuting axes $`a_i`$ gives the exact coefficient

```math
 \begin{aligned}
 \widehat b_W(z)&=\sum_{u,v,t}
 \overline{c_R(v)}c_L(z,u)\,
 2^{-|\alpha(z,u)|/2}i^{\Phi(z,u,t)}\\
 &\qquad\cdot\mathbf1_{t\le\alpha(z,u)}
 \mathbf1_{r(v)=\ell(z,u)+At},\\
 \alpha_i(z,u)&=\omega(a_i,\ell(z,u)),\qquad
 At=\sum_i t_i a_i .
 \end{aligned}
```

Here Phi is the literal band cocycle of Section 1. The formula follows
by applying that band entry to each left history and matching its final
label with the right history. It retains all collisions coherently.
Arbitrary Clifford walls act on both the labels and their literal phases;
no prefix-return premise or bound on physical dirty width is used.
For a normalized L2 profile error epsilon, a sufficient uniform coefficient
error is $`\epsilon 2^{-m/2}`$, by Parseval. Testing significance before
including E and its actual suffix does not have this guarantee.

### A matched pair does not determine later insertions

Use the packets $`O_M`$ of Section 5 and put
$`A_w=Z^w\otimes I`$ and $`E_w=e^{-i\pi A_w/8}`$. A literal native
T on this axis has the same conjugation. With $`\gamma=1/\sqrt2`$,

```math
 g_M(w):=\tau(O_M E_wO_ME_w^\dagger)
 =\frac{1+\gamma}{2}
  +\frac{1-\gamma}{2}(-1)^{w\cdot\boldsymbol1}
       2^{-|M^{\mathsf T}w|} .
```

Indeed, $`O_M^2=I`$ makes the commutator term in the conjugation have
zero trace. The remaining terms are
$`\cos^2(\pi/8)+\sin^2(\pi/8)\tau(O_MA_wO_MA_w)`$.
The Pauli marginal from Section 5 evaluates the last trace as
$`(-1)^{w\cdot\boldsymbol1}2^{-|M^{\mathsf T}w|}`$.
Thus pairs with old $`K=0`$ and overlap one can respond differently
to the same next band.

More precisely, any one stored description that answers every such
insertion to error below $`(1-\gamma)/16`$ needs
$`\Omega(q^2)`$ bits. Matrices related by a column permutation have
the same response. Otherwise
$`N=(M')^{\mathsf T}M^{-\mathsf T}`$ sends some unit vector to a
vector of weight at least two. Taking $`w=M^{-\mathsf T}e_i`$
separates the responses by at least $`(1-\gamma)/8`$.
There are $`|\mathrm{GL}(q,2)|/q!\ge2^{q^2-2-q\log_2q}`$ classes.
The original matched words extend to identity, but the hypothetical
insertions in this query family need not return work. This excludes
that reusable description requirement, not a description recomputed
from the actual complete word under its terminal promise.

### One actual band can remove a full-rank relative defect

Choose invertible M and $`M'`$ with $`M-M'`$ invertible. All axes in
the combined lists $`B(M),B(M')`$ commute. A dependence with turn
vectors u,v has second-bank label $`u+v=0`$ and first-bank label
$`Mu+M'v=0`$, so u=v=0. Thus the signed product

```math
 E=E_B(M')E_B(M)^\dagger
```

is one actual native band on 2q independent axes, with 2q T gates
and $`O(q^2)`$ Clifford cost. A CNOT change of basis in the common
commuting span, with first-bank Hadamards, realizes that band.
Its literal action replaces the left packet:

```math
 E E_B(M)E_A=E_B(M')E_A,\qquad
 E_A^\dagger E_B(M')^\dagger E E_B(M)E_A=I .
```

The completed five-band word therefore returns every input exactly,
including arbitrary dirty and reference inputs. It has 6q T gates,
2q physical wires and $`O(q^2)`$ elementary Cliffords. A fixed fanout
exposes the single logical probe as in Section 5; unused logical wires
and one clean spectator may be added. The old relative K has rank q,
and the old overlap is at most $`(3/4)^q`$, but the completed overlap
is one. This is an actual full-return obstruction to retaining the
old low-rank-K conclusion after insertion unchanged.

Such pairs exist for every $`q\ge2`$. Take $`M'=MC`$, where C is a
direct sum of companion blocks for $`t^2+t+1`$ and $`t^3+t+1`$.
Every $`q\ge2`$ is a sum of twos and threes, and neither block has
eigenvalue zero or one. M itself can be any invertible matrix.

This replacement is a positive control uniformly for all input probes:
it is an operator identity, not just the displayed overlap. Successive
replacements telescope to the final two-band packet. Against a fixed
right matrix $`M_R`$, the new relative matrix is simply
$`M_R^{\mathsf T}M_{\rm new}^{-\mathsf T}-I`$.
The existing [paired-half profile theorem](RELATIVE_FOUR_BAND_BIAS.md)
then applies with its original endpoint hypotheses and no accumulated
approximation error. This closed subclass always reduces back to a
four-band comparison; it does not establish closure for arbitrary bands.

Even the simpler inverse insertion revives a small contraction. For
$`Q=E_APE_A^\dagger`$, the unique matching turns give

```math
 \tau\bigl(QE_B(M)QE_B(M)^\dagger\bigr)
   =\left(\frac{1+1/\sqrt2}{2}\right)^q,
 \qquad
 E_B(M)^\dagger O_M E_B(M)=Q .
```

The first equality averages $`2^{-|M^{\mathsf T}u|/2}`$ over uniform u.
The second turns that contraction into one and completes the literal
four-band identity. Small partial contractions cannot be discarded
without accounting for the actual remaining inverse.

### A scalar calculation in the common commuting span

There is a further exact reduction when the inserted band's axes lie in
$`\mathcal H=\{X^a\otimes Z^b\}`$. Put
$`F=E_B(M')^\dagger E E_B(M)`$ and $`\omega_8=e^{i\pi/4}`$.
In the simultaneous X basis on the first bank and Z basis on the second,
F is diagonal with eigenphase $`\omega_8^{p(x,z)}`$, where the actual
signed parity polynomial is

```math
 \begin{aligned}
 p(x,z)&=c+|M^{\mathsf T}x+z|-|(M')^{\mathsf T}x+z|\\
 &\quad+\sum_{\ell=1}^h\varepsilon_\ell
 [a_\ell\cdot x+b_\ell\cdot z]_2\pmod8 .
 \end{aligned}
```

Here $`\varepsilon_\ell\in\{1,-1\}`$ records the actual orientation;
any constant from a signed axis is included in c. Both exterior B
words are part of p. For the fixed probe $`P=X^{\boldsymbol1}\otimes I`$,

```math
 \tau(O_{M'} E O_M E^\dagger)
   =\mathbb E_z\left|\mathbb E_x\omega_8^{p(x,z)}\right|^2 .
```

To prove this, write the contraction as $`\tau(QFQF^\dagger)`$.
After the first-bank Hadamards, Q is
$`((Z-Y)/\sqrt2)^{\otimes q}\otimes I`$. It preserves z, and
every first-bank entry has squared magnitude $`2^{-q}`$.
Expanding the trace therefore gives
$`2^{-3q}\sum_{z,x,x'}\omega_8^{p(x',z)-p(x,z)}`$, which is the
displayed formula. An arbitrary affine input probe need not have these
flat entries, so this calculation does not give a uniform whole-profile
formula without additional weights.

The replacement band above makes p identically zero even when the old
K has full rank. More generally, scalar value one only forces p to be
independent of x for each z. A nonconstant remaining dirty phase in z
is invisible to this scalar and need not satisfy terminal work return.
Recording the whole polynomial does not yet give a short description.
Thus the exact kernel, replacement identities and this controlled
calculation leave the required closure question open: compress the
actual completed contraction uniformly over the full affine probe
family, using its terminal return promise, with useful dependence on
depth and T rather than on the large physical dirty width.

The [bounded native fixture](../../tests/test_native_layer_count_distribution.py)
checks the single-axis response, a nonzero shared-span phase polynomial,
inverse revival and an odd-size replacement with $`q=3`$. Its common
Clifford frame emits the replacement as one signed band, and exact
native arithmetic verifies every column of the completed identity.

<a id="terminal-quotient-shared-span"></a>

## 7. Terminal return compresses a transverse commuting middle

The whole commuting family of Section 6 has a stronger conclusion than
its fixed-probe phase average. The relative geometry matters. Fix w
physical wires, with a logical wires, b arbitrary dirty wires and
$`c\in\{0,1\}`$ clean wires, so $`w=a+b+c`$. Let

```math
 \begin{aligned}
 A&=\prod_{j=1}^q T_j^{\epsilon_j},\qquad
 H_Q=\prod_{j=1}^q H_j,\qquad B=H_QAC,\\
 V&=B^\dagger D_pB
   =C^\dagger A^\dagger F A C,\qquad
 F=H_QD_pH_Q,\\
 D_p|y\rangle&=\omega_8^{p(y)}|y\rangle,\qquad
 \omega_8=e^{i\pi/4}.
 \end{aligned}
```

Here C is any common Clifford embedding, the q outer T axes are
distinct, and $`\epsilon_j\in\{1,-1\}`$. The actual middle F is a
product of any number of native bands whose axes lie in the fixed
transverse span $`\{X^u\otimes Z^v\}`$, with u on the first q wires
and v on the remaining wires. Its signed parity polynomial p includes
all those bands, orientations and literal scalar phases. Every such p
has degree at most three modulo eight: linear coefficients are arbitrary,
quadratic coefficients are even, and cubic coefficients are multiples
of four. Restriction along an affine binary map preserves this form.

This geometry includes all the shared-commuting insertions of Section 6
and arbitrary common Clifford conjugates. Intrinsically, let the actual
independent outer axes span $`L_A`$ and let a diagonal Lagrangian
$`L_F`$ contain the middle axes. Full crossing rank q, equivalently
$`L_A\cap L_F=0`$, supplies these coordinates without replacing the
individual outer axes. Choose $`f_i\in L_F`$ dual to the actual
$`A_i`$, so $`\omega(A_j,f_i)=\delta_{ij}`$. Complete them by a
basis of $`L_F\cap L_A^\perp`$ and extend to a symplectic basis.
A Clifford sends each actual $`A_i`$ to $`Z_i`$, its dual to $`X_i`$,
and the remaining diagonal basis to the remaining Z coordinates.
Axis signs are absorbed into the recorded orientations.

Dependent crossing does not justify this reduction while keeping the
same individual native T axes. Independent unmatched exterior
Cliffords and noncommuting middle insertions are also outside the
following result.

Let J initialize the c clean wires. Suppose that the actual completed
word satisfies, for some logical unitary U,

```math
 \|VJ-J(U\otimes I_b)\|_{\rm op}\le\eta
 \lt\kappa:=\frac{1-1/\sqrt2}{2}.
```

A common scalar may be absorbed into U. No prefix is assumed to return
work. For any affine Pauli probe family $`P_z`$ on the logical and dirty
input, with $`z\in\mathbb F_2^r`$ and literal phase digits of degree at
most two, and any fixed output Pauli S, the entire initialized profile

```math
 \beta(z)=2^{-(w-c)}\mathrm{Tr}
 \left(JP_zJ^\dagger V^\dagger S V\right)
```

has an exact description using

```math
 O\bigl(w(a+r+1)+(a+1)^3+r^2\bigr)
```

bits. The r parameters label the whole probe family; they are not the
physical width. A mathematical input slice, such as target input zero
or a dirty Pauli eigenspace, is handled by the corresponding constant
number of Pauli-probe sums. The decoder may take exponential time.

### A large monomial dirty subgroup supplies translations

Represent the dirty Pauli space modulo scalar phases as a binary space
of dimension 2b. Keep its subspace $`\mathcal D_0`$ whose conjugate
under C has zero X mask on the first q wires. Its dimension is at least
$`2b-q`$. These Paulis commute with A, so for every $`D\in\mathcal D_0`$,

```math
 Q_D:=BDB^\dagger=H_QCDC^\dagger H_Q
```

is a monomial Pauli in the y basis. Its translation vector consists of
the first-q Z mask and the remaining X mask of $`CDC^\dagger`$.
The kernel of this translation map on $`\mathcal D_0`$ has dimension
at most $`w-q`$: only the remaining Z coordinates can survive there.
The image H consequently satisfies

```math
 \dim H\ge2b-w=w-2(a+c),\qquad
 \mathrm{codim}\,H\le2(a+c).
```

This permits dense Clifford embeddings, large connected components and
large dirty width. Only the eventual small quotient will be recorded;
the full dirty-Pauli basis and C's quadratic-size tableau are not.

### Initialized return forces exact phase periods

Put $`\Pi=BJJ^\dagger B^\dagger`$. Every dirty D preserves the
initialized code and commutes with the ideal logical action. Applying
the complete error bound twice therefore gives

```math
 \|[D_p,Q_D]BJ\|_{\rm op}\le2\eta .
```

This does not apply a clean-code bound after a leaked prefix. It is a
commutator estimate for the completed word on its initialized input.
Also $`[Q_D,\Pi]=0`$. Write
$`Q_D|y\rangle=\zeta_D(y)|y+s_D\rangle`$ and
$`d_y=\omega_8^{p(y)}`$. The norm of row y of the preceding commutator
on the code is exactly

```math
 \left\|\langle y|[D_p,Q_D]BJ\right\|
   =|d_y-d_{y+s_D}|\sqrt{\Pi_{yy}}
   \le2\eta .
```

For c=0, every $`\Pi_{yy}=1`$. For c=1, conjugate the clean Z Pauli
through C, A and $`H_Q`$. Each first-q X or Y factor contributes a
diagonal coefficient of magnitude $`1/\sqrt2`$; a mandatory flip
makes the whole diagonal zero. Thus the transformed clean Pauli either
is a diagonal Pauli, or its diagonal entries have magnitude at most
$`1/\sqrt2`$. In the second case every $`\Pi_{yy}\ge\kappa`$.
In the first case Pi is exactly the projector onto one affine
hyperplane $`\mathcal S`$ in the y coordinates. Dirty translations
preserve this hyperplane, and $`\Pi_{yy}=1`$ on it. Use
$`\mathcal S=\mathbb F_2^w`$ in the full-support cases.

Unequal eighth roots have distance at least
$`2\sin(\pi/8)=2\sqrt\kappa`$. The strict error threshold therefore
forces $`d_y=d_{y+s_D}`$ for every y in the relevant support and every
$`D\in\mathcal D_0`$. Equivalently, p is exactly constant modulo eight
on the H cosets in $`\mathcal S`$. The conclusion is exact even though
the completed initialized return was approximate. There is no sum of
errors over generators: every translation in H comes from a dirty Pauli
to which the same two-error bound applies.

### A fixed quotient section gives a uniform decoder and update

Let $`k=\dim(\mathcal S/H)`$. Then $`k\le2(a+c)`$; in the diagonal
one-clean case the stronger bound is $`k\le2a+1`$. Choose a linear
quotient map R and a section $`y=y_0+Lt`$ such that
$`RL=I_k`$ and the section meets every H coset in $`\mathcal S`$ once.
On an affine clean slice, extend R to the ambient space and choose
$`y_0\in\mathcal S`$ and L within its direction space. Define

```math
 \bar p(t)=p(y_0+Lt),\qquad
 p'(y)=\bar p\bigl(R(y-y_0)\bigr).
```

The quotient polynomial has $`O((k+1)^3)`$ bits: its constant and linear
coefficients are modulo eight, its quadratic coefficients modulo four
after factoring out two, and its cubic coefficients modulo two after
factoring out four. Phase constancy gives $`D_{p'}BJ=D_pBJ`$ exactly.
In the full-support case it gives equality on every physical input.
In the diagonal clean-slice case it makes no assertion about clean-one
action, which is absent from the initialized profile.

Here is a constructive description, including the embedding costs.
Lift $`P_z`$ by identity on the clean wires, writing the result as
$`\widetilde P_z`$. Record q and its orientations, the affine labels of
$`C\widetilde P_zC^\dagger`$, the fixed Pauli $`CSC^\dagger`$, and
the clean Pauli $`CZ_cC^\dagger`$ when present. For the input phase,
record its original two degree-two binary digits and the Clifford
cocycle separately, and add their values modulo four when decoding.
The Clifford cocycle has affine low bit and quadratic high bit on the
affine labels, by ordered Pauli multiplication. Their sum need not
retain quadratic digits because of
the binary carry, but these separate records still cost $`O(r^2)`$.
The labels and phases together cost $`O(w(r+1)+r^2)`$ bits.
Record R using $`O(wk)`$ bits and the
small polynomial $`\bar p`$. The clean-slice equation, if any, is
known from the recorded clean Pauli. A fixed Gaussian-elimination rule
can choose $`y_0,L`$ from R and that equation, so a second large
embedding is unnecessary. The encoder uses this same section when
forming $`\bar p`$.

The decoder conjugates the recorded Paulis by the known $`H_QA`$,
reconstructs Pi, and evaluates the exact trace with $`D_{p'}`$.
Writing $`\mathcal P_z=B\widetilde P_zB^\dagger`$
and $`\mathcal S_0=BSB^\dagger`$, the trace is

```math
 \beta(z)=2^c\,2^{-w}\mathrm{Tr}
 \left(\Pi\mathcal P_zD_{p'}^\dagger\mathcal S_0D_{p'}\right).
```

This evaluates every z with the same recorded data and keeps the actual
probe-dependent matrix entries. No flat-entry assumption from the
single-probe average is used. It proves the asserted description bound
without recording a full Clifford tableau. If the clean projector has
full support, the same record also describes the unconditional affine
Pauli profile. For a diagonal clean slice, only the initialized profile
is claimed; arbitrary action on the rejected clean sector is not encoded.

The quotient section depends on the fixed outer embedding and register
allocation, not on the middle phase polynomial. An additional commuting
native band with signed parity phase e updates the record by

```math
 \bar p(t)\longmapsto\bar p(t)+e(y_0+Lt)\pmod8 .
```

Each axis restricts to one affine parity in k variables, so this is an
explicit update of the constant, linear, quadratic and cubic coefficients.
It needs no collection of new pairings with all earlier axes. The section
record is defined and updates exactly even when intermediate words fail
return and have different off-section actions. Only the actual final
return promise licenses reconstructing its completed profile from the
record. It is not a summary answering arbitrary hypothetical future
insertions without that promise.

For the counting consequence, take Walsh probes on r tested input bits
and write $`N=2^r`$ for the number of entries in their inverse Walsh
profile. At $`w=O(\sqrt N)`$ and $`a=O(\log N)`$, the description
length is $`o(N)`$. Under the displayed complete-return promise, this
gives only $`2^{o(N)}`$ initialized bias profiles, including a fixed
number of the input slices above. Thus this entire transverse family
cannot realize the robust signs of typical prescribed tables with
$`\Theta(N)`$ entries. This holds even with arbitrarily many commuting
middle bands. It does not count raw decisions that leave arbitrary
garbage. Adding more commuting bands in the same geometry does not
escape the returned-profile count.

The dependent-crossing case has a concrete missing premise. On two
wires take a diagonal basis generated by $`Z_1,Z_2`$ and the one native
band $`A=T(X_1)T(X_1Z_2)`$. Its two axes are independent and commute,
but have crossing rank one. Then

```math
 AZ_1A^\dagger
   =-Y_1\frac{I+Z_2}{2}+Z_1\frac{I-Z_2}{2}.
```

The transformed clean-zero projector is non-diagonal but has a zero
diagonal entry; its support is not an affine hyperplane. Thus the
uniform kappa-versus-affine-slice argument cannot simply be reused.
[Section 8](#terminal-quotient-crossing-defect) handles this geometry
using a pointwise Clifford gap and support-invariant cosets. Independent
unmatched exterior Cliffords and an arbitrary noncommuting next band
remain outside that extension.

The [exact native fixture](../../tests/test_native_layer_count_distribution.py)
uses three wires to compare every coefficient of an eight-member affine
probe family with the quotient decoder. It covers both a full-support
clean projector with diagonal weights $`(1\pm1/\sqrt2)/2`$ and an
affine clean slice whose rejected-sector action differs. It retains a
nonreturning intermediate dirty phase, its actual cancellation, and a
Clifford phase carry with cubic high bit. A two-wire native band checks
the deficient-crossing projector above. The all-size rank bound,
error threshold and description count are analytic.

<a id="terminal-quotient-crossing-defect"></a>

## 8. Terminal profiles at growing crossing defect

The affine-support dichotomy in Section 7 is unnecessary for one clean
input. A pointwise gap survives for every clean image of an actual
native T band. It permits a quotient description with explicit cost
for deficient crossing, including non-affine clean support.

Retain the logical, dirty and clean allocation, the completed
initialized-isometry error below kappa, and the affine probe family
of Section 7. Consider a completed word

```math
 V=C_0^\dagger A_0^\dagger F_0 A_0 C_0,
```

where $`C_0`$ is a common exterior Clifford, $`A_0`$ is one actual
native band of independent commuting signed Pauli T rotations, and
all actual middle axes belong to one isotropic binary Pauli span M.
The middle may contain arbitrarily many native bands in that span.
It may also retain an actual diagonal Clifford phase in the same
commuting Pauli algebra. Write $`T_F`$ for the number of actual
emitted middle T or inverse-T gates; Clifford factors are not charged
to that count.
For one completed word, M may be the span of all its actual middle
axes, including those generating a retained Clifford phase. For an
insertion sequence, fix the allowed M in advance.

First cancel every individual outer T factor whose axis lies in
$`M^\perp`$. It commutes with every middle gate and all other outer
factors, so its actual inverse and forward occurrence cancel across
$`F_0`$. This leaves the complete word exactly unchanged. For an
insertion sequence, cancel only factors commuting with all of the
fixed allowed M. Such a cancellation then remains valid after every
allowed insertion. This removes harmless outer-axis padding before
measuring the defect.

Let $`L_A`$ be the span of the retained actual outer axes and define

```math
 q=\dim L_A,\qquad
 \delta=\dim(L_A\cap M),\qquad s=q-\delta .
```

Under the same complete-return promise
$`\eta\lt\kappa=(1-1/\sqrt2)/2`$, the entire initialized profile has
an exact description using

```math
 \begin{aligned}
 u&=a+\delta+1,\\
 L_{\rm prof}&=O\bigl(w(a+r+1)+s\delta+r^2+
                  \min\{u^3,T_Fu+u^2\}\bigr)
 \end{aligned}
```

bits. The cubic choice is independent of the number of middle bands;
the sparse choice retains their actual emitted T count. The
description has a fixed-section update for every actual added band
in M, without a prefix-return assumption. Its decoder may take
exponential time. For a pure product of T bands, the sparse phase
term can be sharpened from $`T_Fu+u^2`$ to $`T_Fu+1`$.

### Choose the completion and retain each actual outer gate

The defect above is intrinsic to the retained outer axes and the
specified middle span; it need not depend on an arbitrary Lagrangian
completion. In the symplectic quotient $`M^\perp/M`$, the image

```math
 \overline E=\bigl((L_A\cap M^\perp)+M\bigr)/M
```

is isotropic. Extend it to a Lagrangian and choose a complementary
Lagrangian. The lift of this complement is a Lagrangian $`L_F`$
containing M, with

```math
 L_A\cap L_F=L_A\cap M.
```

Indeed any vector in $`L_A\cap L_F`$ lies in $`M^\perp`$ and
projects into both complementary spaces, so its image is zero and
the vector lies in M. Thus the crossing rank against $`L_F`$ is s.

Select s actual outer axes whose crossing functionals on $`L_F`$
are independent. Label these pivots $`P_1,\ldots,P_s`$, and label
the remaining actual axes $`P_{s+1},\ldots,P_q`$. There is a binary
matrix $`\Gamma\in\mathbb F_2^{\delta\times s}`$ for which the
literal commuting products

```math
 E_j=P_{s+j}\prod_{i=1}^sP_i^{\Gamma_{ji}},
 \qquad 1\le j\le\delta,
```

form an independent basis of $`L_A\cap L_F`$. Their Pauli signs
are included in these actual products. Choose $`f_i\in L_F`$ dual
to the pivot axes, complete them by the $`E_j`$ and the rest of
$`L_F\cap L_A^\perp`$, and extend to a symplectic basis. A Clifford
sends the pivots to $`Z_i`$, the $`f_i`$ to $`X_i`$, the $`E_j`$
to $`Z_{s+j}`$, and the remaining Lagrangian basis to the remaining
Z coordinates. Each nonpivot actual axis then becomes
$`Z_{s+j}\prod_iZ_i^{\Gamma_{ji}}`$.

Consequently, after absorbing this coordinate change into the common
Clifford C, the actual outer word is

```math
 A_\Gamma=\prod_{i=1}^sT(Z_i)^{\epsilon_i}
 \prod_{j=1}^{\delta}
 T\!\left(Z_{s+j}\prod_{i=1}^sZ_i^{\Gamma_{ji}}\right)^{\epsilon_{s+j}}.
```

No actual T factor has been replaced by rotations about a different
basis of its span. Axis signs are recorded in the orientations;
any resulting common outer scalar cancels between the actual inverse
and forward word. With $`H_s=\prod_{i=1}^sH_i`$, the middle and
complete word become

```math
 F=H_sD_pH_s,\qquad B=H_sA_\Gamma C,\qquad V=B^\dagger D_pB.
```

The signed parity polynomial p retains the literal middle scalar and
has the same weighted cubic form modulo eight as in Section 7.
The normal-form data $`s,\delta,\Gamma`$ and the orientations cost
$`O(s\delta+w)`$ bits.

### Every supported clean row has a uniform weight

For c=1, write
$`\Pi=BJJ^\dagger B^\dagger=(I+Q_c)/2`$, where
$`Q_c=BZ_cB^\dagger`$. A physical T band sends a Pauli to a literal
Hermitian Clifford, as in the
[native normal form](FOUR_LAYER_SCALAR_GAP.md#two-band-signed-gap-sources).
The conjugated band $`A_\Gamma`$ retains this property, and the
surrounding Cliffords preserve it. Thus $`Q_c`$ is a Hermitian Clifford.

For any fixed Clifford, all nonzero computational matrix entries have
one modulus $`2^{-h/2}`$ for an integer h. To see this, its zero
column is a stabilizer state, uniform in modulus on an affine support;
every other column is a Pauli translate of that state. Hermiticity
makes each diagonal entry of $`Q_c`$ real. If h=0, the projector
diagonal lies in $`\{0,1/2,1\}`$. If h is positive, every diagonal
entry of $`Q_c`$ has magnitude at most $`1/\sqrt2`$. Therefore

```math
 \Pi_{yy}\in\{0\}\cup[\kappa,1],
 \qquad \kappa=\frac{1-1/\sqrt2}{2}.
```

For c=0 the same assertion holds with $`\Pi=I`$. This is a
pointwise statement; the nonzero support may be non-affine, as in
the two-qubit example of Section 7.

### Dirty translations preserve that support

Keep the dirty-Pauli subspace $`\mathcal D_0`$ whose C-conjugates
commute with all q retained actual outer axes. It has dimension at
least $`2b-q`$. In the displayed coordinates their X masks vanish
on the first q wires, since those actual axes still span the first
q Z labels. Hence

```math
 Q_D=BDB^\dagger=H_sCDC^\dagger H_s
```

is a monomial Pauli. Its translation kernel has dimension at most
$`w-s`$: after imposing zero translation, only the Z labels after
the first s positions remain. Its translation image H satisfies

```math
 \dim H\ge(2b-q)-(w-s)=w-2(a+c)-\delta,
 \qquad k:=\mathrm{codim}\,H\le2(a+c)+\delta.
```

Every such $`Q_D`$ commutes with Pi because a dirty Pauli preserves
the initialized input code. Thus $`\Pi_{y+h,y+h}=\Pi_{yy}`$ for
every h in H. The set
$`\mathcal S=\{y:\Pi_{yy}\ne0\}`$ is a union of whole H cosets.
Its shape need not be affine.

Applying the completed return estimate twice, exactly as in Section 7,
gives

```math
 \|[D_p,Q_D]BJ\|_{\rm op}\le2\eta,
 \qquad
 |d_y-d_{y+h}|\sqrt{\Pi_{yy}}\le2\eta,
 \qquad d_y=\omega_8^{p(y)}.
```

On every supported row, distinct eighth roots would make the left
side at least $`2\sqrt\kappa\sqrt\kappa=2\kappa`$.
The strict error bound therefore forces exact equality
$`p(y+h)=p(y)\pmod8`$ on every supported H coset. Every h is the
translation of an actual dirty Pauli, so errors are not accumulated
along a chosen generating path. No initialized estimate is applied
after an unreturned prefix.

### Use one ambient quotient section

The translations in H have zero entries on the delta defect
coordinates $`s+1,\ldots,q`$: the first-q X masks were zero, and
$`H_s`$ leaves those defect positions unchanged. Thus H lies in
the coordinate space $`\mathcal E=\{h:h_{s+1}=\cdots=h_q=0\}`$.
Write

```math
 k_0=\mathrm{codim}_{\mathcal E}H=k-\delta\le2(a+c).
```

With the defect and remaining coordinates in their fixed order,
choose a rank-$`k_0`$ map $`R_0`$ on the remaining $`w-\delta`$
coordinates with kernel H, and its fixed Gaussian-elimination right
inverse $`L_0`$. The ambient quotient and section are

```math
 R(y)=(y_{\rm def},R_0y_{\rm rem}),\qquad
 L(t_{\rm def},t_0)=(t_{\rm def},L_0t_0).
```

All delta defect coordinates remain in the quotient; only their
identity embedding is implicit. R has rank k and kernel H. Define

```math
 \bar p(t)=p(Lt),\qquad p'(y)=\bar p(Ry),\qquad RL=I_k.
```

If y is supported, $`LRy`$ lies in its H coset and is also supported.
Exact phase constancy then gives $`p'(y)=p(y)\pmod8`$. Unsupported
rows of BJ are zero. Consequently

```math
 D_{p'}BJ=D_pBJ
```

exactly, even when p has different values off the initialized support.
A global linear section is sufficient; no list of the supported cosets
or separate support embedding is stored.

Record $`\Gamma`$, the outer orientations, the affine Pauli labels
of $`C\widetilde P_zC^\dagger`$, the fixed output image $`CSC^\dagger`$,
and the clean image $`CZ_cC^\dagger`$ when present. As in Section 7,
keep the original input phase digits separately from the Clifford
cocycle and add them modulo four during decoding; their carry need
not have quadratic digits. These probe data cost
$`O(w(r+1)+r^2)`$ bits. Record only $`R_0`$, using
$`O((w-\delta)k_0)=O(w(a+1))`$ bits, and the weighted cubic
$`\bar p`$ using $`O((k+1)^3)`$ bits. L follows from $`R_0`$ by
the fixed rule. This dense phase representation gives

```math
 \begin{aligned}
 O\bigl(s\delta+w(r+k_0+1)+r^2+(k+1)^3\bigr)
 =O\bigl(w(a+r+1)+s\delta+r^2+u^3\bigr).
 \end{aligned}
```

There is also a sparse phase representation at the actual count.
Each emitted middle T or inverse-T gate contributes one signed affine
parity to p, including the constant bit from its signed Pauli axis.
Restriction to $`y=Lt`$ leaves one signed affine parity in k variables,
which costs $`O(k+1)`$ bits. Thus the actual list costs
$`O(T_F(k+1))`$ bits. A separately retained diagonal Clifford phase
has even linear coefficients and quadratic coefficients divisible by
four modulo eight. Affine restriction preserves this form, so its
record costs $`O((k+1)^2)`$ bits. That term can be omitted for a pure
product of T bands, apart from a constant-size literal scalar record.
No arbitrary phase coefficient or uncharged T factor is supplied.

One mode bit selects the shorter phase representation. Its cost is

```math
 O\!\left(\min\{(k+1)^3,
                   T_F(k+1)+(k+1)^2\}\right).
```

Together with the embedding and probe data, this proves the stated
bound for $`L_{\rm prof}`$. The sparse list refers to the actual
emitted rotations; cancellation of their restricted functions does
not reduce the physical $`T_F`$ used in the resource comparison.

The decoder obtains Pi and the transported probes directly by
conjugating the recorded Pauli images through the known
$`H_sA_\Gamma`$, then uses the initialized trace formula in Section 7
with $`D_{p'}`$. It records neither a full Clifford tableau nor a
quadratic-size dirty-Pauli basis. The equality needed by that trace
holds on every initialized column. Arbitrary rejected-clean-sector
action is not part of this profile claim.

For a fixed allowed M, the cancelled factors, normal form, H, R and L
depend only on the outer word and allocation. An actual commuting band
with signed parity phase e updates

```math
 \bar p(t)\longmapsto\bar p(t)+e(Lt)\pmod8.
```

This updates the same constant, linear, quadratic and cubic coefficient
list in dense mode, or appends the actual restricted signed parities
in sparse mode. A fixed count cap can choose the mode in advance.
Alternatively, begin with the sparse list and convert to cubic
coefficients once that representation becomes shorter. The emitted
T count is nondecreasing under append, so no reverse conversion is
needed. A retained diagonal Clifford phase updates its quadratic
coefficient record. Intermediate words may fail return and differ
off the section. Only the final complete-return promise licenses
reconstructing their completed initialized profile. An insertion
outside M is outside this update rule.

### A quantitative obstruction when the defect grows

For the same Walsh-table reduction as in Section 7, take
$`N=2^r`$, $`w=O(\sqrt N)`$ and $`a=O(\log N)`$, with fixed
constants in those caps. With unrestricted middle count, the dense
record is $`o(N)`$ whenever $`\delta=o(N^{1/3})`$. With the paid
count $`T_F=O(\sqrt N)`$ and $`s\le w`$, the sparse record instead gives

```math
 L_{\rm prof}
   =O\bigl(\sqrt N(\log N+\delta+1)
            +(\log N+\delta+1)^2\bigr).
```

It is $`o(N)`$ whenever $`\delta=o(\sqrt N)`$. In either regime,
only $`2^{o(N)}`$ initialized profiles are possible, including a
fixed number of the mathematical input slices used in Section 7.
They cannot realize the robust signs of typical prescribed N-entry
tables under the complete-return promise.

The constants in these description bounds are uniform. For the paid
count, take $`w\le C_w\sqrt N`$, $`T_F\le C_T\sqrt N`$ and cap
the defect by $`\varepsilon\sqrt N`$. The leading code length is
bounded by a constant times
$`((C_w+C_T)\varepsilon+\varepsilon^2)N`$, plus $`o(N)`$.
A sufficiently small fixed epsilon makes this a fixed fraction of N.
Only an exponentially small fraction of all sign tables can then
occur. The dense count gives the same conclusion with a sufficiently
small constant multiple of $`N^{1/3}`$. Thus a typical qualifying
table requires, after cancellation of individually commuting outer
factors,

```math
 \begin{aligned}
 \delta&=\Omega(N^{1/3})
       &&\text{with unrestricted middle count},\\
 \delta&=\Omega(\sqrt N)
       &&\text{when }T_F=O(\sqrt N).
 \end{aligned}
```

These are necessary intersection sizes in every admissible
outer/middle representation. The defect is not a unitary invariant
under arbitrary resynthesis, and increasing it by harmless middle
padding does not establish a construction. The paid result is not a
growing native-depth lower bound: a wide outer band can already have
extensive defect. That remaining regime at square-root count,
unmatched exterior Cliffords, and noncommuting middle insertions
still require further arguments. The theorem also does not count
raw decisions that leave arbitrary garbage or supply an efficient
native compiler.

The [four-wire exact fixture](../../tests/test_native_layer_count_distribution.py)
retains a one-clean projector with non-affine support, a middle phase
periodic only on that support, and a quotient agreeing on every
initialized column while differing on the rejected sector. It compares
the entire eight-member affine probe profile with the decoder and
retains the literal phase carry. The pointwise Clifford gap,
completion lemma, growing-defect count and actual-band update are
analytic.
