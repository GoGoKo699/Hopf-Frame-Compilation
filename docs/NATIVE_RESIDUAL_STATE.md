# Native residual state preparation with two flags

[State compiler](STATE_ONLY_COMPILER.md) · [Native tables](NATIVE_RESIDUAL_ROTATION.md) · [Verification](VERIFICATION.md)

The [bounded emitter](../compiler_robust_hopf/native_residual_state.py)
composes two certified tables into a one-system-qubit state-preparation
word. Both clean compiler flags are used. The amplification calls the
actual inverse of the emitted half-amplitude word and retains all
intermediate leakage. This implements the state identity already proved
in [the state compiler, Sections 3–4](STATE_ONLY_COMPILER.md#3-two-flags-give-an-exact-half-amplitude-state).
It adds no new asymptotic claim.

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

The [tests](../tests/test_native_residual_state.py) use certified dyadic
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

This completes the selected one-system-qubit amplification integration.
General table sizing and dirty lookup remain unimplemented. The next
bounded step is an additional unchanged address bit: emit and verify its
lookup and predicate with every occupied flag and borrowed helper
charged before composing a larger state compiler. The coherent QBP
reference/target branch and the full fine state-preparation schedule
remain separate integration tasks. The complete-frame endpoint and
claims of end-to-end advantage are unchanged.
