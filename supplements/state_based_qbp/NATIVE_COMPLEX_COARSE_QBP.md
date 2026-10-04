# A small native integration of both complex Hopf gradient streams

[Complex coarse compiler](COMPLEX_COARSE_COMPILER.md) · [Complex-gradient decoder](COMPLEX_COARSE_QBP.md) · [Real native integration](NATIVE_COARSE_QBP.md) · [Verification map](../../docs/VERIFICATION.md)

This two-qubit example emits the preparation, controlled observable, actual
coarse inverse, and readout for the complex magnitude stream, together with
the separate phase stream. Every operation is an elementary Clifford+T
gate. A dirty helper is actively used and returned; no controlled rotation
or supplied observable-response vector replaces a circuit operation.

The chosen target is exactly native. This is a finite-size integration of
the state compiler's fallback, including literal complex branch phases. It
does not emit a general addressed residual table or amplification circuit,
close the complete-frame endpoint, or demonstrate a compilation advantage.

## 1. An unwrapped phase table and its literal gauge

Let $`W_{\mathbb R}`$ be the prescribed two-qubit real Hopf frame. The
example uses balanced $`(1,1,1)`$, signed $`(1,-1,0)`$, and singular
$`(0,1,-1)`$ real angle tuples, in units of $`\pi/4`$. Each observable
and decoder is fixed at the selected tuple when taking derivatives.

Use arithmetic subtree mean $`\mu=.37+4\pi`$ and the three heap-ordered
phase-prefix angles

```math
(\gamma_1,\gamma_2,\gamma_3)=(9,-5,7)\pi/4.
```

The supplied, unwrapped leaf phases are

```math
(\varphi_{00},\varphi_{01},\varphi_{10},\varphi_{11})
=\mu+( -4,-14,2,16)\pi/4.
```

With $`R_z(\gamma)=\operatorname{diag}(e^{-i\gamma},e^{i\gamma})`$,
apply the root row on the high bit, then the two high-bit-selected child
rows on the low bit. Every row acts on all its suffix copies. Their product
is exactly

```math
D_0=e^{-i\mu}D_\varphi=\operatorname{diag}(-1,i,i,1),
\qquad V=D_0W_{\mathbb R},\qquad |\psi'\rangle=V|00\rangle.
```

In particular, $`D_0=-(S^\dagger\otimes S^\dagger)`$; the emitted
prefix word retains that literal minus sign. A scalar cannot be discarded
from only one component of a controlled or branch-selected word.

The stored winding representatives and mean are retained even though this
particular diagonal has a short expression. The circuit prepares the
gauge-fixed state psi', without emitting the scalar $`e^{i\mu}`$.
Its energy and physical phase derivatives agree with those of the supplied
state $`D_\varphi W_{\mathbb R}|00\rangle`$. Both branches and all
magnitude coefficients use this same gauge.

## 2. Actual complex coarse rows

Use the determinant-one commutator word from the
[real native fixture](NATIVE_COARSE_QBP.md#1-an-explicit-complex-coarse-word):

```math
K_0=[T,HTH],\qquad K_{r+1}=[K_r,SK_rS^\dagger],
\qquad [A,B]=ABA^\dagger B^\dagger.
```

Matrix products act from right to left. Every inverse is the literal
reversed, adjointed elementary word. The unexpanded one-qubit K3 contains
256 T or T-adjoint gates; address controls and dirty echoes add separately
charged gates.

Embed K3 on the real root pair $`(|00\rangle,|10\rangle)`$, with
identity on the other pair, and call this E. Replace each ideal phase row
by the actual native block

```math
\widehat P_j=R_z(\gamma_j)K_3.
```

Apply these blocks to every suffix copy in heap order to obtain
$`\widehat P`$, and set $`C=\widehat P W_{\mathbb R}E`$. These
phase rows are generally nondiagonal. Neither the quantum inverse nor the
classical decoder substitutes ideal diagonal phases for their actual
blocks. The recorded real root block is $`R_y(\theta_1)K_3`$; the two
real child blocks remain their exact native rotations.

Writing $`\delta_K=\|K_3-I\|`$, the two phase layers contribute at
most $`2\delta_K`$: the two child errors act on disjoint sectors, so
their layer norm is the maximum, not their sum. Including E gives the
complete-input estimate

```math
\|C-V\|\le3\delta_K.
```

This is an operator bound on the full logical space, not only a prepared
column. To certify the required constant, direct two-by-two multiplication
in $`\mathbb Q(\sqrt2,i)`$ gives

```math
\operatorname{tr}K_0=\frac12+\sqrt2,\qquad
\operatorname{tr}K_1=1+\frac{11\sqrt2}{16},\qquad
\operatorname{tr}K_2=\frac{1+2893\sqrt2}{2048}.
```

For a determinant-one two-by-two unitary, $`\|K-I\|^2=2-\operatorname{tr}K`$.
Since $`\sqrt2\gt41/29`$ (because $`1681\lt1682`$),

```math
\|K_2-I\|^2
=\frac{4095-2893\sqrt2}{2048}
\lt\frac{71}{29696}\lt\frac1{400}.
```

The unitary commutator inequality
$`\|[A,B]-I\|\le2\|A-I\|\|B-I\|`$ therefore yields

```math
\delta_K\lt\frac1{200},\qquad
\|C-V\|\lt\frac3{200}\lt\frac1{64}.
```

At $`N=4`$ this also meets the $`1/(4\sqrt N)`$ coarse promise. The
sharper K2 calculation matters: tripling the earlier fixture's coarser K3
bound alone would not certify $`1/64`$. The
[exact certificate regression](../../tests/test_complex_coarse_certificate.py)
checks these field identities with rational coefficients, independently of
the emitter and floating-point arithmetic.

## 3. Coherent preparation with a returned dirty helper

The six wires are numbered from the least significant bit:

| Wires | Role and input |
|---|---|
| 0, 1 | Logical system, initially zero for preparation |
| 2 | Protocol branch; arbitrary for preparation, plus for execution |
| 3, 4 | Reserved compiler flags; untouched in this finite-size fallback |
| 5 | Arbitrary dirty helper, actively used and returned |

First select E on branch zero and apply the common real frame. At each
phase row, in heap order, select its K3 on branch zero and then apply the
exact prefix Rz to both branches. Thus branch zero receives the actual
coarse row, while branch one receives the ideal row, without duplicating
the common phase rotation. For arbitrary branch amplitudes this prepares exactly

```math
a|0\rangle_B C|00\rangle+b|1\rangle_B V|00\rangle.
```

Initializing the branch in plus yields the required relative plus sign.
The two reserved flags supply no hidden scratch in this circuit. Arbitrary
external references and additional untouched dirty slots can be appended;
this finite fixture is not a new workspace bound.

The native implementation uses the inherited reflection alphabet
$`\{T^jHT^{-j}\}_{j=0}^7\cup\{Z\}`$. For each reflection R, load
its prefix predicate f into the unknown helper z, apply R controlled on
the selected branch beta and helper, unload f, and repeat the same
controlled reflection. Its net exponent is

```math
\beta(z\oplus f)\oplus\beta z=\beta f.
```

For an empty prefix, the load is an X on the helper. The echo completes
before the next reflection; grouping a noncommuting word into a single
echo would not implement this identity. All controls are expanded into
elementary gates, including the exact seven-T Toffoli. No helper input is
assumed to be zero.

The literal Clifford $`B=HS^\dagger`$ satisfies $`BYB^\dagger=Z`$.
It turns selected $`R_y(\pi/4)=HZ`$ reflection words into selected
$`R_z(\pi/4)`$ words, preserving their phases. Together with reflection
conjugations this avoids any opaque controlled-T primitive.

## 4. Two fixed observables and both measurement streams

Use the Hermitian-unitary observables

```math
O_+=W_{\mathbb R}(TXT^\dagger\otimes I)W_{\mathbb R}^\dagger,
\qquad
O_H=W_{\mathbb R}(I\otimes H)W_{\mathbb R}^\dagger.
```

Their controlled words include both real-frame wrappers. There are no
phase-frame wrappers around these observables. Consequently the complex
phase table affects the responses, and the retained cases exercise
nonzero magnitude and phase derivatives. The observable itself remains
fixed when differentiating the state.
For example, at the singular tuple $`(0,\pi/4,-\pi/4)`$, $`O_H`$ has
the raw gradients

```math
\nabla_\theta E=(0,\sqrt2,0),\qquad
\nabla_\varphi E=(-1/\sqrt2,1/\sqrt2,0,0).
```

Both streams therefore carry a nonzero signal even though half the target
leaves have zero amplitude.

For the magnitude stream, apply controlled O to branch one, the actual
$`C^\dagger`$ to both branches, and Hadamards on both system bits.
Choose X or Y branch measurement with equal probability. Decode using
the actual complex coefficients from
[the magnitude mean identity](COMPLEX_COARSE_QBP.md#2-magnitude-means-and-the-unchanged-depth-bound).
In particular, keep both the Y contributions and the correction by D0
and the actual C.

For the phase stream, prepare psi' on both branch components, apply the
same controlled O, then measure the branch in Y and the computational
leaf x. The record is $`2s e_x`$ for a single coefficient-one observable.
There is no inverse frame in this stream. Its signed means are all four
raw phase derivatives, including the redundant common-phase direction.

The original magnitude protocol is emitted with the same exact V and
the same controlled O. The original phase execution is the same phase
word just described.

## 5. Literal gate counts

For the balanced real tuple and $`O_+`$, the expanded words contain:

| Circuit word | T and T-adjoint | Clifford gates |
|---|---:|---:|
| Exact gauge-fixed V | 26 | 246 |
| Actual coarse C | 5,402 | 25,942 |
| Coherent preparation, before branch initialization | 20,506 | 48,630 |
| Controlled observable | 14 | 101 |
| Complete corrected magnitude X execution | 25,922 | 74,677 |
| Complete corrected magnitude Y execution | 25,922 | 74,678 |
| Complete original magnitude execution | 66 | 597 |
| Complete phase execution, either protocol | 40 | 350 |

Replacing $`O_+`$ by $`O_H`$ adds ten Clifford gates to its controlled
observable and each complete execution; the T counts are unchanged.
The complete executions include branch initialization and readout. Counts
retain all observable wrappers and make no cancellation or optimality
claim. The original magnitude word is much cheaper in this fixture,
whose target has a short exact implementation. The long coarse perturbation
serves to exercise complex correction and borrowed-work return.

## 6. Bounded verification and reproduction

The initialized preparation embedding has four columns, covering every
branch/helper basis input in a 64-dimensional circuit space. Comparing
the whole isometry keeps every output row, including leakage and relative
branch phase. Readout uses two columns to retain arbitrary dirty input.
Sandwiching decoded scores between these columns checks the entire
two-by-two operator on the dirty input against the analytic gradient
times identity, rather than only checking classical helper probabilities.
The fine and coarse phase tables, with either selected branch or no branch
selection, are also checked on all 16 active system/branch/helper inputs
in four-column batches. These checks include inactive branches and every
logical suffix sector.

The complex histogram utility contracts signed integer X/Y counts with
the actual real and phase blocks and a division-free Hopf reverse pass.
The phase record uses its signed computational-leaf histogram directly.
The final numerical contractions use ordinary floating point; the
separately proved rounding guarantee is not a certified arithmetic backend.

These are bounded deterministic checks, without shot Monte Carlo or large
simulation. Analytic gate identities establish the exact circuit contracts;
numerical residuals check their implementation. The report's largest
propagated batch has 64 rows and four columns, with a $`10^{-9}`$
floating-point check tolerance. General fine-precision resource claims
remain the analytic compiler's responsibility.

Run from the repository root:

```bash
python scripts/complex_coarse_native_example.py
python scripts/complex_coarse_native_example.py --case singular --format json
python -m unittest tests.test_native_complex_coarse_qbp tests.test_complex_coarse_certificate -q
```

See the [native emitter](../../compiler_robust_hopf/native_complex_coarse_fixture.py),
[complex histogram utility](../../compiler_robust_hopf/complex_coarse_decoder.py),
and [finite regression checks](../../tests/test_native_complex_coarse_qbp.py).
