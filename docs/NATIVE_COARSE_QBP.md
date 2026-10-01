# A small native integration of the coarse-frame decoder

[Coarse-frame proof](COARSE_FRAME_QBP.md) · [State-only preparation](STATE_ONLY_COMPILER.md#7-a-complex-native-reference-from-the-same-coarse-word) · [Verification map](VERIFICATION.md)

This example connects an emitted Clifford+T preparation, a complex controlled
observable, the actual coarse inverse, and the corrected X/Y decoder. The
logical system has two qubits. Every controlled operation is expanded into
elementary gates, and an arbitrary dirty helper is used and returned.

The target frame is exactly native at this special angle tuple. The example
therefore exercises the state compiler's finite-size fallback, rather than
its addressed residual table and amplification construction. The general
fine-precision resource theorem remains analytic. This example is an
integration and phase check, not evidence of a small-instance advantage.

## 1. An explicit complex coarse word

Take the two-qubit Hopf tuple with all three angles equal to $`\pi/4`$.
Each local real rotation is exactly $`R_y(\pi/4)=HZ`$, including its
literal phase. Let W be the prescribed complete frame.

Define the one-qubit commutator words as matrix products, with the rightmost
factor acting first:

```math
K_0=[T,HTH],\qquad
K_{r+1}=[K_r,SK_rS^\dagger],\qquad
[A,B]=ABA^\dagger B^\dagger.
```

These are determinant-one native words. Their inverses are obtained by
reversing the actual gate word and taking each gate's adjoint. Direct
two-by-two algebra gives

```math
\|K_0-I\|=1-\frac1{\sqrt2}.
```

For example, with $`\omega=e^{i\pi/4}`$, the diagonal and off-diagonal
entries of HTH are $`(1+\omega)/2`$ and $`(1-\omega)/2`$.
They give $`\operatorname{tr}K_0=2-(1-1/\sqrt2)^2`$; for a
determinant-one two-by-two unitary, the displayed distance follows.

For unitaries A and B,

```math
\|[A,B]-I\|=\|AB-BA\|
\le2\|A-I\|\|B-I\|.
```

Consequently

```math
\|K_3-I\|
\le128\left(1-\frac1{\sqrt2}\right)^8
\lt\frac1{64}.
```

This last comparison has a rational certificate: $`\sqrt2\gt7/5`$
implies $`(1-1/\sqrt2)^2\lt1/10`$, so the bound is smaller than
$`8/625\lt1/64`$.

The unexpanded one-qubit K3 word contains 256 T or T-adjoint gates.
This number is distinct from the larger count after adding address controls
and dirty-helper echoes. No synthesis search or precision extrapolation is
used to select the word.

Embed K3 on the root pair $`(|00\rangle,|10\rangle)`$, acting identically
on the other pair, and call that operator E. Set $`C=WE`$. This changes
only the root's logical two-mode block, so C retains a recorded tree of three
two-by-two blocks. It is complex and satisfies the complete-input bound
$`\|C-W\|=\|E-I\|\lt1/64`$. This is the coarse promise needed at
$`N=4`$, not merely closeness of the prepared column.

## 2. Coherent preparation and the dirty echo

The six named wires, numbered from the least significant bit, are:

| Wires | Role and input |
|---|---|
| 0, 1 | Two logical system qubits, initially zero for preparation |
| 2 | Protocol branch; arbitrary for the preparation contract, plus for execution |
| 3, 4 | Reserved compiler flags; untouched in this finite-size fallback |
| 5 | Arbitrary dirty helper; its state is never assumed to be zero |

Prepare by selecting E on branch zero and then applying the common W. On
the initialized system this gives the exact coherent pair

```math
|0\rangle_BC|00\rangle+|1\rangle_BW|00\rangle,
```

with each branch multiplied by its input amplitude. Initializing B in plus
therefore supplies the required normalized reference/target superposition.
The compiler flags are unused, rather than supplying hidden initialized
scratch. Additional unused dirty slots can be appended identically; one
active helper here is not a new asymptotic workspace theorem.

The selected root word is rewritten in the inherited reflection alphabet
$`\{T^jHT^{-j}\}_{j=0}^7\cup\{Z\}`$. Let f be the low-system-bit
predicate selecting the root pair, let beta be the branch-zero predicate,
and let z be the unknown helper bit. For each reflection R, the circuit
loads f into z, applies R controlled on beta and z, unloads f, and applies
that same controlled R again. The resulting exponent is

```math
\beta(z\oplus f)\oplus\beta z=\beta f.
```

Thus the helper returns exactly for both of its basis values and all coherent
superpositions. Complete each reflection's echo before the next reflection;
echoing a whole noncommuting word would be incorrect.

Conjugating a controlled Z by the literal word
$`V=SHTHS^\dagger`$, for which $`VZV^\dagger=H`$, implements controlled
H. Target T conjugations supply the other reflections. Twice-controlled Z
is an H-conjugated Toffoli, and every Toffoli is expanded into its exact
seven-T native word. No controlled-T gate or opaque rotation is admitted.

## 3. A charged complex observable and two readouts

Use the two concrete Hermitian-unitary observables

```math
O_+=W(TXT^\dagger\otimes I)W^\dagger,\qquad
O_Y=W(Y\otimes I)W^\dagger.
```

The first has raw gradient $`(\sqrt2,0,0)`$ at the balanced tuple; the
second has zero real-state gradient. Their controlled implementations use
uncontrolled basis conjugations around one branch-controlled X. This
preserves literal phase and charges both W wrappers; there is no supplied
observable-response vector in the execution circuit.
The W appearing inside each observable is fixed at the selected tuple.
Only the state is differentiated; varying the observable together with the
state would be a different gradient problem.

After the controlled observable, apply the actual inverse of C and a
Hadamard on each system bit. The reference is then uniform. Measure the
branch in X or Y with equal probability, and use the actual-C scores from
[the decoder proof](COARSE_FRAME_QBP.md#2-universal-means-from-randomized-xy-interference).
The Y contributions are nonzero even though the target is real. Omitting
them, or using the ideal Walsh coefficients in place of the corrected
coefficients, yields a biased result in these examples.

The original frame protocol is emitted for the same W and exactly the same
controlled observable. Both protocols have exact quantum words in this
fixture and use the same declared wire layout. Their counts are literal
counts of the expanded words, without claiming optimal circuit reduction.
The original protocol is cheaper here: its target already has a short exact
frame implementation. The long coarse perturbation deliberately exercises
complex corrections and borrowed-work return.

For the balanced tuple and $`O_+`$, the emitted counts are:

| Circuit word | T and T-adjoint | Clifford gates |
|---|---:|---:|
| W | 6 | 50 |
| C | 774 | 3,794 |
| Coherent preparation, before branch initialization | 5,126 | 12,659 |
| Controlled observable | 14 | 101 |
| Complete corrected X execution | 5,914 | 16,557 |
| Complete corrected Y execution | 5,914 | 16,558 |
| Complete original execution | 26 | 205 |

These are unsimplified counts of the supplied words. In particular, the
controlled-observable wrappers are retained in both executions; this table
is not a minimal compilation or an asymptotic comparison.

## 4. What the finite verification establishes

The initialized preparation embedding has four columns, one for every
branch/helper basis input, in a 64-dimensional circuit space. Comparing
the entire output isometry with the prescribed branch pair retains every
possible leakage row. The norm of the resulting four-column error covers
arbitrary branch/dirty/reference correlations, including relative phase.
The selected root word is also checked on all 16 active system/branch/helper
basis inputs, in four-column batches. Combining their errors checks the
complete operation, including every inactive low-bit sector.

For readout, two columns retain arbitrary dirty input after preparing the
branch in plus. Sandwiching each decoded score between these output columns
gives a two-by-two observable on the original dirty input. Its equality to
the analytic gradient times identity checks coherent dirty-input independence,
as well as the mean for each classical dirty basis state.

The histogram implementation uses exact Python-integer Walsh butterflies,
then applies the three actual logical blocks and a division-free reverse
Hopf traversal. It avoids a dense Jacobian. Its final contractions use
ordinary NumPy floating-point arithmetic; the separately proved arbitrary
precision rounding guarantee is not an implemented certified arithmetic
backend.

These are bounded deterministic calculations, with no shot Monte Carlo or
large simulation. Analytic gate identities establish exactness; finite
floating-point residuals check their implementation. They do not certify
the asymptotic theorem numerically or supply a general native frame/state
compiler.

## 5. Reproduce and inspect

Run from the repository root:

```bash
python scripts/coarse_frame_native_example.py
python scripts/coarse_frame_native_example.py --case singular --format json
python -m unittest tests.test_native_coarse_qbp tests.test_coarse_frame_decoder -q
```

The [native emitter](../compiler_robust_hopf/native_coarse_fixture.py) fixes
the register layout, word order, controls, and finite counts. The
[histogram utility](../compiler_robust_hopf/coarse_frame_decoder.py) accepts
the ordered actual two-by-two coarse blocks and signed X/Y histograms.
The shot count is supplied separately because signed counters lose canceled
records. It never constructs the full frame or derivative matrix.

The report's largest propagated batch has 64 rows and four columns. Its
gate identities and distance bound are analytic; the reported residuals
use ordinary floating point and a $`3\times10^{-10}`$ finite-check
tolerance. The script also compares a fixed signed-record histogram with
direct scores; those fixture counters are explicitly not sampled data.
