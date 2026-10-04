# Native residual state preparation with two flags

[State compiler](STATE_ONLY_COMPILER.md) · [Native tables](NATIVE_RESIDUAL_ROTATION.md) · [Verification](../../docs/VERIFICATION.md)

The [one-qubit emitter](../../compiler_robust_hopf/native_residual_state.py)
composes two certified tables into a one-system-qubit state-preparation
word. Both clean compiler flags are used. The amplification calls the
actual inverse of the emitted half-amplitude word and retains all
intermediate leakage. This implements the state identity already proved
in [the state compiler, Sections 3–4](STATE_ONLY_COMPILER.md#3-two-flags-give-an-exact-half-amplitude-state).
The [two-qubit extension](#5-two-system-qubit-preparation-and-a-returned-core-helper)
uses four-row tables and an enlarged initial reflection with a returned
arbitrary core helper. Neither adds a new asymptotic claim.

## 1. Coefficients and initialized inputs

Write the normalized residual target as

```math
|\phi\rangle=a|0\rangle+\frac{w}{\sqrt2}|1\rangle,
\qquad |a|^2+\frac{|w|^2}{2}=1,\qquad |w|\le1.
```

The two complex inputs to `emit_residual_state(a, w, q)` are pairs of
exact dyadic rational approximations. For each pair the caller supplies
the same evaluation promise as the residual-table emitter:

```math
|a_0-a|\le e,\qquad |w_0-w|\le e,
\qquad e=2^{-2q-20},\qquad q\ge5.
```

The input must describe this same normalized state. The code rejects
floats, malformed inputs, invalid table coefficients, and an obviously
inconsistent normalization. In particular, with
$`D=|a_0|^2+|w_0|^2/2`$, compatibility requires

```math
|D-1|\le3e+\frac32e^2.
```

Indeed, each squared modulus changes by at most $`2e+e^2`$.
Passing this necessary check does not certify the caller's external
evaluation or normalization promise.

Qubit zero is the low integer bit. The physical allocation is:

| Wires | Role | Required input |
|---|---|---|
| 0 through q | Precision core | Arbitrary borrowed state |
| q+1 | Synthesis signal | Arbitrary borrowed state |
| q+2 | Flag t and table target | Zero |
| q+3 | Logical system x and table address | Zero |
| q+4 | Flag s and table enable | Zero |

There are two clean compiler flags, one initialized system qubit, and
$`q+2`$ borrowed wires. The latter can be entangled with each other and
an external reference. Neither the initial nor the good reflection tests
the borrowed wires. No extra helper or reset is used.

The coefficient condition is sufficient for this bounded state identity;
it does not require a coarse circuit. To connect it to the small residual
of the state theorem, additionally supply the actual coarse C and its
closeness promise, then apply C after this word. The emitter here ends
at phi and does not emit C or a gradient protocol.

## 2. The two tables and their half-amplitude column

Use the existing completion

```math
U(z)=\begin{pmatrix}z&-\sqrt{1-|z|^2}\\
\sqrt{1-|z|^2}&\overline z\end{pmatrix}.
```

Both tables have target t, address x, and enable s:

| Table | Enabled sector | Row x=0 | Row x=1 |
|---|---|---|---|
| M0 | s=0 | U(a) | U(a) |
| M1 | s=1 | U(0) | U(w) |

Each emitted table preserves x and s and is exactly identity in its
disabled sector, on the complete borrowed-work space. Thus their product
M has the required four rows. In particular, active U(0) retains its
rejected amplitude; it is not an identity row.

Put $`K=C_{s=1}(H_x)`$ and

```math
Q=H_s M K H_s,\qquad
P_{\rm good}=I_x\otimes|00\rangle\langle00|_{s,t}.
```

The rightmost operation acts first. The emitted chronological word is
`H_s, K, M0, M1, H_s`. For the ideal tables,

```math
P_{\rm good}Q|0_x00_{s,t}\rangle
=\frac12\left(a|0\rangle+\frac{w}{\sqrt2}|1\rangle\right)|00\rangle.
```

Let $`\delta_j`$ denote each native table's rational operator-error
certificate. Because the enabled sectors are disjoint and the inactive
actions are exactly identity,

```math
\|\widehat M-M\|\le\max(\delta_0,\delta_1),\qquad
\|\widehat Q-Q\|\le\delta_Q:=\max(\delta_0,\delta_1)
\lt130\,2^{-q}.
```

This maximum is valid before the final Hadamard changes s. Conjugation
by the exact surrounding gates preserves the bound. No independent
phase adjustment of the table sectors is permitted.

## 3. Literal amplification and charged reflections

Let J append $`|0_x00_{s,t}\rangle`$ to arbitrary borrowed input.
The exact reflections are

```math
R_{\rm init}=I-2JJ^\dagger,\qquad
R_{\rm good}=I-2P_{\rm good}.
```

The returned chronological stages are

```text
Q, R_good, actual_inverse(Q), R_init, Q, minus_identity
```

Equivalently, their operator is

```math
\widehat A=-\widehat Q R_{\rm init}\widehat Q^\dagger
R_{\rm good}\widehat Q.
```

The inverse reverses every emitted elementary gate and replaces T/S by
their adjoints. It is not a separately synthesized inverse target. The
ideal half-amplitude identity and one amplitude-amplification step give
$`AJ=|\phi,00\rangle\otimes I_{\rm dirty}`$. This is the standard
[Brassard–Høyer–Mosca–Tapp amplification law](https://arxiv.org/abs/quant-ph/0005055),
with the literal phase and dirty-input contract derived in the
[state proof](STATE_ONLY_COMPILER.md#4-one-state-amplification-step-with-actual-inverses).
Unitary telescoping therefore gives

```math
\left\|\widehat A J-|\phi,00\rangle\otimes I_{\rm dirty}\right\|
\le3\delta_Q\lt390\,2^{-q}.
```

The norm includes returned flags, borrowed work, and arbitrary reference
correlations. The circuit's action on other initial logical or flag
states is not the promised target isometry.

Every exact operation is expanded into elementary gates:

| Operation | Native implementation | T/TDG count |
|---|---|---|
| Controlled H | $`B\,CZ\,B^\dagger`$, with $`B=SHTHS^\dagger`$ and $`BZB^\dagger=H`$ | 2 |
| Initial reflection | X-conjugated CCZ on x,s,t, using the literal seven-T Toffoli | 7 |
| Good reflection | X-conjugated CZ on s,t | 0 |
| Leading minus sign | Literal X,Z,X,Z word on t | 0 |

The emitted unsimplified counts are therefore

```math
T(\widehat Q)=2(540q+630)+2=1080q+1262,
```

```math
T(\widehat A)=3T(\widehat Q)+7=3240q+3793.
```

Gate storage is $`O(q)`$. With $`q=L+10`$, the error is below
$`2^{-L}`$ and this particular word uses $`L+12`$ borrowed wires.
That exceeds the theorem's basic small-system reservation
$`L+n+7=L+8`$ at n=1; it fits the banked reservation
$`2(L+n+7)`$. This integration does not replace the theorem's
separately retained minimum-budget small-system fallback.

## 4. Evidence and the next boundary

The [tests](../../tests/test_native_residual_state.py) use certified dyadic
approximations to exactly normalized rational coefficient pairs. Small
native table blocks retain their literal sector phases and all borrowed
inputs. The complete preparation isometry then retains every intermediate
flag and borrowed-work component through the actual inverse and both
reflections. Independent elementary propagation also checks the flattened
word, while exact high-precision ledgers require no large statevector.

These checks complement the analytic error proof; coarse q=5 matrices
alone do not certify fine-precision accuracy. The literal minus sign,
the system bit in the initial reflection, and the active-zero table row
are essential to the ideal identity.

The one-system-qubit integration is complete. The following extension
uses the [four-row lookup](NATIVE_RESIDUAL_ROTATION.md#6-four-rows-with-two-address-bits-and-no-additional-helper)
and charges the enlarged reflection explicitly.

## 5. Two-system-qubit preparation and a returned core helper

The [two-qubit emitter](../../compiler_robust_hopf/native_two_qubit_residual_state.py)
exports `emit_two_qubit_residual_state(a, tails, q)`. The three entries of
`tails` approximate the scaled complex coefficients w1,w2,w3, in the
same exact-dyadic pair format as a. The target is

```math
|\phi\rangle=a|0\rangle+\frac12\sum_{j=1}^3w_j|j\rangle,
\qquad |a|^2+\frac14\sum_{j=1}^3|w_j|^2=1,\qquad |w_j|\le1.
```

Every supplied pair approximates its true coefficient within
$`e=2^{-2q-20}`$, with $`q\ge5`$. These evaluation and normalization
promises remain external. The necessary normalization sanity check is

```math
\left|\,|a_0|^2+\frac14\sum_{j=1}^3|w_{j,0}|^2-1\right|
\le\frac72e+\frac74e^2.
```

It sums the squared-modulus error $`2e+e^2`$ with weights totaling
$`7/4`$; passing it does not establish the external promises.

| Wires | Role | Required input |
|---|---|---|
| 0 through q | Precision core | Arbitrary borrowed state |
| q+1 | Synthesis signal | Arbitrary borrowed state |
| q+2 | Flag t and table target | Zero |
| q+3, q+4 | System bits x0,x1 and table addresses | Zero |
| q+5 | Flag s and table enable | Zero |

The row index is $`x=x_0+2x_1`$. There are two initialized system
qubits, two clean compiler flags, and q+2 borrowed wires. The allocation
is q+6 physical wires. Core wire zero is reused only during the initial
reflection below; it is not an additional reservation.

Use four-row table M0 with rows (a,a,a,a), enabled at s=0, and M1 with
rows (0,w1,w2,w3), enabled at s=1. Put

```math
K=C_{s=1}(H_{x_0}H_{x_1}),\qquad Q=H_sM_1M_0KH_s.
```

The two controlled Hadamards cost four T/TDG gates in total. Their active
branch is uniform on all four system labels. Thus, with good flags s=t=0,

```math
P_{\rm good}Q|00_x00_{s,t}\rangle
=\frac12\left(a|0\rangle+\frac12\sum_{j=1}^3w_j|j\rangle\right)|00\rangle.
```

The same chronological amplification as Section 3 uses the actual inverse
of this Q. The good reflection and literal minus sign are unchanged.
The initial reflection now tests all four logical zeros, excluding the
entire borrowed pool.

### Exact initial reflection on arbitrary leaked inputs

Let h be core wire zero and define exact native Toffoli words

```math
F=\operatorname{CCX}(x_0,x_1;h),\qquad
G=\operatorname{CCX}(h,s;t).
```

The chronological echo F,G,F,G has the Boolean action

```math
h_{\rm out}=h,\qquad
t_{\rm out}=t\oplus s(h\oplus x_0x_1)\oplus sh
=t\oplus sx_0x_1.
```

All other wires are unchanged. The Toffolis have literal phase one, so
linearity establishes the same identity on arbitrary superpositions and
reference correlations. No assumption that h is zero is used.
Conjugating this word by H on t gives the four-body controlled Z.
Conjugating that by X on x0,x1,s,t gives exactly

```math
R_{\rm init}=\left(I-2|0000\rangle\langle0000|_{x,s,t}\right)
\otimes I_{\rm dirty}.
```

This reflection uses four seven-T Toffolis, hence 28 T/TDG gates; the
other gates are Clifford. It returns h exactly even when Q or its actual
inverse has entangled it with the logical registers. Omitting the last G
would retain an unwanted dependence on the incoming helper bit.
The source core is reused between complete table subroutines, with no
reset, projection, or assumed cleanup.

### Error, literal count, and verification

If the two table certificates are $`\delta_0,\delta_1`$, their disjoint
enabled sectors give $`\delta_Q=\max(\delta_0,\delta_1)`$. The exact
reflections and actual inverse retain the complete initialized-isometry
bound

```math
\left\|\widehat A J-|\phi,00\rangle\otimes I_{\rm dirty}\right\|
\le3\delta_Q\lt390\,2^{-q}.
```

J now appends two system zeros and two flag zeros to arbitrary borrowed
input. All work return and flag leakage are included. No other logical
input action or supplied coarse circuit is promised.

M0 has identical rows, so both of its quadratic supports vanish. Let
$`k_z,k_y\in\{0,1,2\}`$ count the nonempty quadratic axes of M1's
outer and middle rotation tables. Including both tables, both controlled
Hadamards, all three Q calls, and the enlarged reflection gives

```math
T(\widehat Q)=1080q+1264+70(2k_z+k_y),
```

```math
T(\widehat A)=3240q+3820+210(2k_z+k_y)\le3240q+5080.
```

Gate storage remains $`O(q)`$. At $`q=L+10`$, error is below
$`2^{-L}`$ with L+12 borrowed wires. This is three more than the
minimum n=2 allocation L+9 in the state theorem, and fits its banked
pool. It does not replace the minimum-budget small-system fallback.

The [two-qubit tests](../../tests/test_native_two_qubit_residual_state.py)
check the reflection on every logical/helper input and the full
2048-by-128 preparation isometry at q=5 through native table sectors.
The actual inverse acts on all intermediate leakage; an independent
source-algebra oracle and direct flattened-word propagation retain
literal phases. Exact coefficient fixtures include a complex normalized
state inside the actual 1/64 coarse radius, with all three tails nonzero;
this certifies its residual neighborhood without supplying a coarse C.
Fine-q checks use rational certificates and emitted counts, not large
statevectors. The q=5 analytic constant is loose and is not the sole
numerical correctness oracle.

One- and two-system-qubit residual preparation are now implemented. A
[coherent one-system-qubit selector](NATIVE_RESIDUAL_BRANCH.md) also
preserves relative phase under a separate arbitrary protocol branch,
which its initial reflection excludes. The
[bounded native residual QBP example](NATIVE_RESIDUAL_QBP.md) integrates
that selector with charged common-coarse QBP and both decoder streams.
General lookup, the full fine state compiler, and its complete QBP
integration remain implementation tasks. The complete-frame endpoint
and end-to-end advantage questions are unchanged.
