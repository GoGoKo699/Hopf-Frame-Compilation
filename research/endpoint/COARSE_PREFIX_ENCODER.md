# A native coarse prefix encoder

The complete encoder in the
[boundary-feedback word](BOUNDARY_PROPAGATION.md#5-complete-feedback-and-the-surviving-prefix-encoder)
has an $`O(N)`$-T approximation at error
$`(43/64)2^{-\lceil N/n\rceil}`$, uniformly on the entire physical
space. This precision is weaker than the selected endpoint
$`\eta=2^{-N}`$. A separate small-conjugation lemma states exactly when
coarse encoder error is suppressed in a complete native word.

Throughout, $`N=2^n`$, $`n\ge3`$, the alphabet is literal
H,S,CNOT,T,T-dagger, and the workspace is two supplied clean qubits and
at most $`N+n+7`$ arbitrary dirty qubits. Operator estimates include
references and do not minimize over scalar phase.

## 1. All prefix frames in one mode register

Place a unitary coin $`Q_v`$ at each internal node of the full tree.
Its first column defines $`\mathcal B|v\rangle`$ on the two children;
its second column is the orthogonal wandering injection. The root and
these wandering columns, propagated to each arrival depth, form the
orthonormal chain basis. At depth d they are precisely the columns of
the complete depth-d zero-suffix frame $`W_d`$: earlier layers leave a
marker untouched, its own layer injects the second coin column, and
later layers propagate the first columns. Thus, in heap-node order,

```math
T_{\rm pre}=|0\rangle\langle0|\ \oplus\bigoplus_{d=0}^n W_d.
```

This is $`T_n`$ from the boundary recurrence, including its entire
depth-n frame. It acts on the n logical qubits and one supplied mode
bit a, which is treated as arbitrary input throughout the encoder.

Use the delimiter permutation

```math
P(0)=0,\qquad
P(2^d+p)=(2p+1)2^{n-d},\quad 0\le d\le n,\quad 0\le p\lt 2^d.
```

The image string is $`p\,1\,0^{n-d}`$. At stage
$`j=0,\ldots,n-1`$, the preceding j most-significant bits are the
unchanged address, bit j is the target, and the lower suffix must have
Hamming weight one. For a delimiter at depth $`d>j`$, that condition
is exactly that the remaining path bits are zero. For $`d\le j`$,
the suffix is zero and the stage is inactive. The delimiter never moves.
The n stage banks therefore implement exactly

```math
PT_{\rm pre}P^\dagger,
\qquad \sum_{j=0}^{n-1}2^j=N-1
```

table rows, simultaneously for every prefix frame and the padded zero.
This layout identity holds for arbitrary U(2) coins; the native estimate
below specializes to the actual regular-core SU(2) coins.

## 2. Exact routing and weight-one controls

To implement P, reverse the $`n+1`$ bits, then, for each possible
least-significant-one position k, reverse the bits above k conditioned
on the unchanged lower pattern $`0\cdots01`$. These sectors are
disjoint; all-zero input is fixed. Each controlled swap is two CNOTs
around an X controlled by the other swapped bit and the sector predicate.
The [exact borrowed-helper construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
gives $`O(n^2)`$ T per swap, with one arbitrary helper returned exactly.
There are $`O(n^2)`$ swaps, so P and its actual inverse cost
$`O(n^4)=O(N)`$ T. They reuse a source-pool wire while the source is
inactive.

For a suffix of length $`\ell=n-j`$, its weight-one predicate is the
XOR of its $`\ell`$ disjoint one-hot minterms. Conditioning a source
$`M=UX_0U^\dagger`$ changes only its central X:

```math
C_h(M)=U C_h(X_0)U^\dagger.
```

Replace that central toggle by the product of the minterm-controlled
toggles, including any scalar-signal literal. This is exact on every
occupied input. It uses one surrounding loader/inverse-loader pair,
with all minterms sharing that source. Each toggle has $`O(\ell^2+1)`$ T cost;
the whole predicate costs $`O(\ell^3+1)`$ and returns a separate dirty
helper exactly. Controls exclude the core target/helper and remain
unchanged during the scalar word. Thus the conditioned primitive is
exactly identity on the inactive sector. Summed predicate cost is
$`O(n^4)=O(N)`$.

## 3. The canonical SU(2) native circuit

For the regular core, use

```math
Q_v=\begin{pmatrix}ih_v&-\bar g_v\\g_v&-i\bar h_v\end{pmatrix}
 =R_y(\pi/2)R_x(\alpha_v)R_z(\alpha_v),
\qquad \alpha_v=\arctan t_v.
```

The two variable factors are fixed Clifford conjugates of the
[full-operator borrowed-signal real bank](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations).
Certified sine and cosine evaluations use
$`t_v/\sqrt{1+t_v^2}`$ and $`1/\sqrt{1+t_v^2}`$, so zero and
singular original charts require no equality oracle. Use the second
supplied wire c as the primitive's arbitrary signal. The guarantee holds
on both its input values. The fixed literal factor $`R_y(\pi/2)=XZ`$
is a predicate-controlled Z followed in time by predicate-controlled X.
Target Hadamards give the former from the same exact controlled-X word.

Choose

```math
s=\lceil N/n\rceil,\qquad q_j=s+n-j+7.
```

At stage j the external dirty reservation is

```math
b=(q_j+1)+j+1=s+n+9\le N+n+7,
```

for the precision core, j selectors, and one predicate helper. The mode
bit a and arbitrary signal c are already the two supplied wires. The
source-pool partition can change between stages without initialization.
Each variable bank has full-operator error below $`43\,2^{-q_j}`$;
ordinary unitary telescoping yields

```math
\boxed{\|\widehat T-T_{\rm pre}\otimes I_{c,b}\|
 \lt 86\sum_{j=0}^{n-1}2^{-(s+n-j+7)}
 \lt \frac{43}{64}\,2^{-s}.}
```

No intermediate exact work return is assumed. All columns, both supplied
wire values, dirty states, and references are included.

The source accounting per variable real bank is fixed by the actual
amplified word:

| Operation | Number per bank | Native cost per operation |
|---|---:|---|
| Source or inverse-source call | 45 | two loaders and a central toggle |
| Actual loader or inverse | 90 | $`2q_j`$ Pauli T rotations |
| Whole-word X/Z query, including inverses | 20 | $`O(2^j)`$ T; $`O(2^jq_j)`$ Clifford |
| Weight-one central toggle | 45 | $`O((n-j)^3+1)`$ T |
| Amplification signal reflection | 4 | Clifford |

The loader row expands the source row; these are not additive duplicate
charges. Two variable factors per coin multiply this ledger by two.
Fixed masks are Pauli gates, and fixed XZ conditioning has the already
priced predicate cost. Hence

```math
T=O\!\left(N+\sum_jq_j+n^4\right)
 =O(N+ns+n^4)=O(N),\qquad
G=O(N(s+n)+n^4).
```

The constants are absolute; finite certified table construction is
required, with no uniform classical runtime promise. This circuit
implements the complete correlated encoder at coarse precision, not a
fine-precision feedback completion.

## 4. Full-output error under a small conjugated gate

Let J initialize a scalar flag, and let E be a unitary on the remaining
registers, including identity on ideal returned work. Suppose the actual
native unitary V obeys
$`\|VJ-JE\|\le\delta`$. Let A preserve that initialized subspace,
$`AJ=JD`$, and suppose $`\|A-I\|\le\kappa`$ on the entire physical
space. Then

```math
\boxed{\|V^\dagger A VJ-JE^\dagger DE\|\le2\kappa\delta.}
```

Multiplying the hypothesis by actual inverse unitaries also gives
$`\|V^\dagger J-JE^\dagger\|\le\delta`$. Cancellation of the
identity part gives the exact decomposition

```math
\begin{aligned}
V^\dagger A VJ-JE^\dagger DE
={}&V^\dagger(A-I)(VJ-JE)\\
&+(V^\dagger J-JE^\dagger)(D-I)E.
\end{aligned}
```

Each summand is bounded by $`\kappa\delta`$. Thus flag leakage and
dirty disturbance are included. The extension
$`A=JDJ^\dagger+(I-JJ^\dagger)`$ is sufficient; smallness of an
accepted block alone is not. If an actual middle word has
$`\|\widehat A-A\|\le\epsilon`$ on the full space, then

```math
\|V^\dagger\widehat A VJ-JE^\dagger DE\|
\le\epsilon+2\kappa\delta.
```

### A native fine-precision middle bank

Let D be an addressed real-rotation bank on the $`n+1`$ node wires,
with every angle of magnitude at most $`\kappa`$. Condition it on
the encoder signal being zero, leaving the other signal sector identity.
Since $`\|R_y(\theta)-I\|\le|\theta|`$, its full extension A has
the required smallness. Fixed Clifford target conjugations also permit
the other rotation axes.

Use the full-operator real primitive at $`q=N+7`$. Of at most n
address bits, put at most three into fixed sector predicates, leaving
$`k\le n-3`$ free bits and at most eight sectors. The dirty allocation is

```math
(q+1)+k+1+1\le(N+8)+(n-3)+1+1=N+n+7,
```

including a separate helper and arbitrary synthesis signal. Neither
supplied wire is assumed available as that signal. Exact inactivity on
unchanged sectors makes their errors combine by a maximum. Therefore
$`\epsilon\lt 43\,2^{-(N+7)}\lt \eta/2`$ with $`O(N)`$ T count. A
weight-one suffix predicate uses the same central-toggle construction,
adding $`O(n^4)=O(N)`$ T. This contract applies to real rotations and
their fixed target conjugates, not arbitrary literal scalar-phase banks.

The [collective precision construction](COLLECTIVE_PRECISION_REFINEMENT.md)
uses this lemma with exact geometric histories to obtain quadratic and
cubic accuracy for the terminal regular-core frame, including the actual
physical work disturbance. The non-small reversal in the feedback word
does not satisfy this lemma's hypothesis.

## 5. Verification

The [coarse encoder tests](../../tests/test_coarse_prefix_encoder.py)
compare the simultaneous layout to independent complete depth-prefix
frames, including general complex unitary coins. They check every
delimiter-permutation input, the exact one-hot predicate truth table, the
geometric precision/workspace ledger, and a full-output conjugation
example with actual flag leakage. These are finite diagnostics of the
layout and algebra; the native cost and reference-safe error statements
follow from the exact primitive contracts and proofs above.
