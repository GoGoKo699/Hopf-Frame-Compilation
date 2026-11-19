# T-count outside the two widest native layers

Put $`n\ge1`$ for the number of logical qubits, $`a`$ for the number of initialized
qubits, $`h=n+a+1`$ and $`N=2^n`$. Allow arbitrary dirty width and arbitrary
Clifford interlayers. For each circuit let $`t`$ be its T-count and $`s`$ the
number left after deleting its two widest physical T layers (all layers if
there are fewer than two). The two layers may be separated. For a compiler,
use common upper bounds $`t,s`$ with the choice of layers made separately in
each circuit. The complete initialized-isometry contract includes dirty
return and external references as in [joint native words](../endpoint/JOINT_NATIVE_WORDS.md).

For all complete real frames at error $`2^{-L}`$, $`L\ge6`$, set
$`E=(N-1)(L-\log_2 6)`$. For all scalar sign tables on $`n`$ bits at
error below $`2\sin(\pi/16)`$, set $`E=N-1`$. In either case

```math
E\le O\!\left(\min\{t(s+h)+h^2,\,(s+h)^3\}\right).
```

At fixed precision, bounded clean supply and total count $`O(\sqrt N)`$,
some targets require $`\Omega(\sqrt N)`$ T gates outside **every** pair
of layers in any compilation meeting that count ceiling. Even with
unlimited sizes of the two layers, $`h=o(E^{1/3})`$ forces
$`s\ge(24E)^{1/3}(1-o(1))`$. Neither conclusion gives growing
depth for an unspecified constant in the total-count ceiling.

## 1. Full return bounds the physical Pauli span

Move the Cliffords to the output, retaining the signed Pauli rotation
axes in physical coordinates. Let $`v`$ be the binary dimension of their
span and let $`c=\cos(\pi/8)`$. For any input stabilizer state, including
references, let $`j`$ be the rank of the axes' commutation syndromes with
its full stabilizer. The overlap of the actual output with any stabilizer
is at most $`c^{2j-t}`$.

To see this, move the input to the computational zero state by a Clifford.
The syndrome becomes the X label. Keep the first chronological pivots
of those labels. For each of the remaining $`t-j`$ rotations use

```math
e^{\pm i\pi P/8}=(I+e^{\pm i\pi P/4})/(2c).
```

In each branch move the resulting quarter-turn Cliffords to the output.
Their transvections only add earlier labels to later ones, hence the
remaining $`j`$ pivot X labels stay independent. Expanding these rotations
in Pauli powers produces distinct basis labels indexed by the chosen pivot
subset. A binary-linear Clifford identifies these with $`|u,0\rangle`$.
The individual signs and ordering of Paulis contribute fourth-root linear
phases and binary quadratic signs, removable by a diagonal Clifford.
The remaining nontrivial state is $`(c|0\rangle+\sin(\pi/8)|1\rangle)^{\otimes j}`$.

Projecting spectator bits of a stabilizer onto zero gives a possibly
subnormalized stabilizer. Its normalized amplitudes are flat on an
affine space of dimension $`d`$. Choose $`d`$ information-set coordinates.
For $`z=\tan(\pi/8)=\sqrt2-1`$, the magnitude of the overlap is at most
$`2^{-d/2}c^j(1+z)^d=c^j`$: the chosen coordinates range over all
$`d`$-bit words and deleting other coordinates only increases $`z^{|u|}`$.
The branch coefficient sum is $`c^{-(t-j)}`$, proving the overlap bound.

Now put the logical input in a basis state, initialize the clean qubits,
and Bell-pair every dirty qubit with its own external reference.
Only logical/clean Z labels can have zero syndrome on physical axes.
Thus $`v\le j+n+a`$. For a target sign, the overlap with this same
stabilizer is at least $`1-\epsilon^2/2\gt c`$ at the stated error.
This gives $`v\le\lfloor t/2\rfloor+n+a`$. For a general frame the
ideal output belongs to the span of the $`N`$ stabilizers obtained by
varying the output logical basis word alone. The actual probability on
this subspace is at least $`1-\epsilon^2`$, so some such overlap is at
least $`\sqrt{(1-\epsilon^2)/N}`$. Consequently $`v\le t/2+O(h)`$.

## 2. A description that does not pay for dirty width

Write $`2r`$ for the symplectic rank on the span of the T axes. A physically
independent commuting layer has size at most $`v-r`$: at most the radical
dimension $`v-2r`$ plus a maximal isotropic dimension $`r`$ in the quotient.
The wider of the selected layers has size at least $`(t-s)/2`$. Applying
the preceding bound with the *actual* circuit lengths gives
$`r\le s/2+O(h)`$. When there are no layers the rank is zero.

Augment the span by the input logical and clean X/Z generators and the
logical output X/Z generators pulled back through the final Clifford.
There are at most $`4n+2a`$ added labels. Each increases the symplectic
rank by at most two. The resulting quotient by the radical has dimension
$`2q`$ with $`q\le s/2+O(h)`$. The radical commutes with every input
logical/clean Pauli, so it lives only on dirty wires. Fix a joint
character of independent commuting Hermitian representatives and use the
normalized projector onto its eigenspace as the dirty input. Symplectic
basis completion represents each listed axis by a signed Pauli on
$`q`$ effective qubits, times identity on any multiplicity space.

These signed axes and boundary labels determine the whole reduced logical
channel on this witness. Expand each input logical density matrix and the
clean-zero projector in the listed input Paulis, and expand any output
measurement in the listed output Paulis. Its expectation is given by the
normalized trace in the selected sector, with the known logical/clean
normalization. The multiplicity space cancels from that trace. The
dirty witness may depend on the circuit, since every dirty input is
covered by the accuracy guarantee.

An explicit ordered list of $`t`$ signed effective axes and $`O(h)`$
boundary labels costs $`O(t(s+h)+h^2)`$ bits, including the lengths.
Alternatively, describe each of the two deleted layers after diagonalizing
its projected commuting axes by a Clifford costing $`O(q^2+1)`$ bits
to specify. Its phase is a signed sum of affine parities modulo eight.
The integer coefficient of a degree-$`d`$ monomial in a parity is
$`(-2)^{d-1}`$ when supported on that parity and zero otherwise. Thus
the diagonal word, even with dependent, repeated or identity projected
axes, takes at most $`\binom q3+2\binom q2+3q+O(1)`$ bits. Keep the other
$`s`$ axes explicitly and mark both bulk-layer locations in this list.
This gives $`2\binom q3+O((s+h)^2)`$ bits, independent of both deleted
widths. It is a channel description and need not itself be a native
circuit. Substitution of $`q\le s/2+O(h)`$ gives the claimed constant 24.

For signs, quotient out a common sign. Each pair of distinct projective
sign tables maps a suitable two-basis-state superposition to orthogonal
outputs. One witness channel cannot approximate both to the stated error,
which is below one half. For frames, take first columns
$`(\sqrt{1-\|x\|^2},x)`$ over the radius-one-half ball in $`N-1`$
dimensions. Inner products are positive and the pure-state trace distance
is at least $`\|x-y\|/\sqrt2`$. A maximal $`3\epsilon`$-separated
packing has at least $`(1/(6\epsilon))^{N-1}`$ points by volume.
The prescribed real family realizes these columns. The reduced output
on the witness has trace distance at most $`\epsilon`$ from its target,
so at least $`2^E`$ descriptions are needed.

## 3. The coefficient for an entire layer-size profile

For layer sizes $`k_j`$ summing to $`t`$, the log-count of witness-channel
descriptions is also at most

```math
\frac{3t^2}{8}-\frac12\sum_j k_j^2+O(ht+h^2).
```

In the quotient of the T-axis span alone, count ordered commuting
$`k`$-tuples by the dimension $`u\le\min(k,r)`$ of their span. Choosing
the space by successively choosing independent orthogonal basis vectors,
dividing by its ordered bases, and then encoding the tuple gives

```math
2^{ku}\prod_{i=0}^{u-1}\frac{2^{2r-i}-2^i}{2^u-2^i}.
```

Its logarithm is at most $`(2r+k)u-3u^2/2+O(u+1)`$; the maximum on
this interval is at $`\min(k,r)`$. Summing dimensions costs $`O(k+1)`$
bits. The leading layer cost is $`2rk-g(k,r)`$, with
$`g(k,r)=(k^2-(k-r)_+^2)/2`$. The concatenated tuple spans, hence
the symplectic group acts freely. Dividing by its order removes
$`2r^2+r+O(1)`$ bits. These orbits encode the Gram matrix of T axes.
Adding boundary rows, signed labels, the length and layer partition costs
$`O(ht+h^2)`$: record all boundary pairings, then choose signs in a
symplectic realization of the enlarged Gram matrix on a radical witness.

Here $`r\le t/4+O(h)`$ and $`r\le t/2-\max_j k_j+O(h)`$ by Section 1.
Clip $`r`$ at $`t/4`$ at cost $`O(ht+h^2)`$ and put
$`\Delta=t/4-r`$. Then $`(k_j-r)_+\le2\Delta+O(h)`$, giving
$`\sum_j(k_j-r)_+^2/2\le\Delta t+O(ht)`$. Since
$`3t^2/8-(2rt-2r^2)=\Delta t+2\Delta^2`$, the profile bound follows.

Alternatively retain only the two largest penalties. Convexity gives
a combined penalty at least $`2g((t-s)/2,r)`$. The remaining radial
bound $`2rt-2r^2-2g((t-s)/2,r)`$ increases up to $`t/4`$. Evaluating
at $`\min(s/2,t/4)`$, up to $`O(ht+h^2)`$, yields the alternative
leading upper bound

```math
F(t,s)=
\begin{cases}
st/2+s^2/4,&s\le t/2,\\
3t^2/8-(t-s)^2/4,&s\ge t/2 .
\end{cases}
```

For nonuniform families take the maximum over allowed profiles before
comparing to $`E`$. The two-parameter function increases in both arguments
on $`0\le s\le t`$. For example, at a hypothetical count ceiling
$`t=(\sqrt{8/3}+o(1))\sqrt E`$, $`h=o(t)`$, the fraction of sign tables
attained with the two largest layers totaling at least $`\delta t`$ is
at most $`2^{-(\delta^2/4-o(1))t^2}`$ for any fixed $`\delta\gt0`$.
