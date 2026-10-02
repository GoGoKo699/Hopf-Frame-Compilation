# Coherent native selection of two residual states

[State preparation](NATIVE_RESIDUAL_STATE.md) · [Common-coarse reference](STATE_ONLY_COMPILER.md#7-a-complex-native-reference-from-the-same-coarse-word) · [Verification](VERIFICATION.md)

The [bounded branch emitter](../compiler_robust_hopf/native_branched_residual_state.py)
prepares either of two one-system-qubit residual states under an arbitrary
protocol branch. It preserves their literal relative phase and includes
all borrowed-work return and flag leakage in one coherent isometry bound.
It implements the branch extension already proved in
[the state compiler, Sections 6–7](STATE_ONLY_COMPILER.md#6-two-states-under-one-protocol-branch).
The separate branch is not initialized by this component and is not a
third clean compiler flag.

## 1. Coefficients, domain, and allocation

`emit_branched_residual_state(branches, q)` takes two entries `(a_j, w_j)`.
Each coefficient is an exact dyadic real/imaginary pair approximating a
true complex coefficient within $`e=2^{-2q-20}`$, with $`q\ge5`$.
For each branch $`j\in\{0,1\}`$, the external promises are

```math
|\phi_j\rangle=a_j|0\rangle+\frac{w_j}{\sqrt2}|1\rangle,
\qquad |a_j|^2+\frac{|w_j|^2}{2}=1,\qquad |w_j|\le1.
```

The emitter checks exact types, dyadic denominators, unit-disk
consistency, and the necessary normalization discrepancy
$`3e+3e^2/2`$ separately in each branch. Passing these checks does not
establish the external approximation or normalization promises.
Coefficients are neither reconstructed nor renormalized.

| Wires | Role | Required input |
|---|---|---|
| 0 through q | Precision core | Arbitrary borrowed state |
| q+1 | Synthesis signal | Arbitrary borrowed state |
| q+2 | Flag t and table target | Zero |
| q+3 | System x and low table address | Zero |
| q+4 | Protocol branch c and high table address | Arbitrary |
| q+5 | Flag s and table enable | Zero |

There are q+6 physical wires: q+2 dirty wires, one initialized system,
two clean compiler flags, and one arbitrary protocol branch. The branch
may be entangled with both the dirty pool and an external reference.
It is not counted as borrowed work either: the target map acts on it.

Let $`J_c`$ append the system zero and the two flag zeros to arbitrary
branch-and-dirty input, and define

```math
G_c=\sum_{j=0}^1|j\rangle\langle j|_c
 \otimes|\phi_j,00\rangle_{x,s,t}\otimes I_{\rm dirty}.
```

The promised error compares $`\widehat A J_c`$ directly with $`G_c`$.
It does not optimize separate global phases in the two branches or
promise a prescribed action on other system/flag inputs.

## 2. Four rows and a branch-independent reflection

Use the [four-row lookup](NATIVE_RESIDUAL_ROTATION.md#6-four-rows-with-two-address-bits-and-no-additional-helper)
with address index $`x+2c`$. Its two tables are

| Table | Enable | Rows in address order |
|---|---|---|
| M0 | s=0 | (a0,a0,a1,a1) |
| M1 | s=1 | (0,w0,0,w1) |

Only x receives the controlled Hadamard. In operator order, set

```math
K=C_{s=1}(H_x),\qquad Q=H_sM_1M_0KH_s.
```

Every complete table leaves the branch label unchanged and retains the
literal phase of its addressed action. Its elementary Toffoli expansion can temporarily target the high address;
preservation is a complete-word identity, not a gate-by-gate claim.
The active zero entries of M1 use the actual completion U(0), retaining
the rejected amplitude. In each branch, the accepted amplitude is
exactly one half of its normalized target. Therefore

```math
P_{\rm good}QJ_c=\tfrac12G_c,
\qquad P_{\rm good}=I_{c,x}\otimes|00\rangle\langle00|_{s,t}
 \otimes I_{\rm dirty}.
```

The exact initial reflection excludes c:

```math
R_{\rm init}=I_c\otimes
 \left(I-2|000\rangle\langle000|_{x,s,t}\right)
 \otimes I_{\rm dirty}=I-2J_cJ_c^\dagger.
```

It is the existing negative-control seven-T CCZ on x,s,t, and requires
no helper. The good reflection tests only s,t. With the literal common
minus sign, one amplification step gives

```math
A=-Q R_{\rm init}Q^\dagger R_{\rm good}Q,
\qquad AJ_c=G_c.
```

This is the same isometry calculation as unbranched preparation, now
with the branch in the input domain of J. It holds on superpositions
and external references by linearity. Including c in the initial
zero test changes the projector and invalidates the amplification on
one branch. Replacing a branch's literal phase by an arbitrary phase
also changes interference, even if its isolated output ray agrees.

## 3. Error and exact resource ledger

The native circuit uses the actual reversed and daggered Q word.
Intermediate leakage remains present through its inverse and both
reflections. Disjoint enabled sectors give the full-operator bound
$`\delta_Q=\max(\delta_0,\delta_1)`$ from the two table certificates.
Unitary telescoping then gives

```math
\|\widehat A J_c-G_c\|\le3\delta_Q\lt390\,2^{-q}.
```

This is one operator norm on the complete branch-and-dirty domain, so
there is no extra sum over branches or restriction to classical branch
inputs. The bound also includes correlations with an arbitrary reference.

M0 depends only on c. Its root rows are repeated in pairs, so its
quadratic mask supports vanish. Let $`k_z,k_y\in\{0,1,2\}`$ count
M1's nonempty quadratic axes for its outer and middle rotations. The
single controlled Hadamard costs two T/TDG gates. Thus

```math
T(\widehat Q)=1080q+1262+70(2k_z+k_y),
```

```math
T(\widehat A)=3240q+3793+210(2k_z+k_y)\le3240q+5053.
```

The latter charges three Q/Q-inverse occurrences and the seven-T
initial reflection. Good reflection and common minus are Clifford.
Storage is $`O(q)`$. At $`q=L+10`$ the error is below $`2^{-L}`$,
using L+12 dirty wires. This exceeds the basic n=1 reservation L+8 by
four wires and fits the banked pool. It does not replace the
minimum-budget small-system fallback.

## 4. Common reference and evidence boundary

Set $`a_0=1,w_0=0`$ for the common-coarse reference case. The residual
output on a plus branch is then

```math
\frac{|0\rangle_c|0\rangle_x+|1\rangle_c|\phi_1\rangle_x}{\sqrt2}
 \otimes|00\rangle_{s,t}\otimes|\xi\rangle_{\rm dirty}
```

within the stated error for every dirty input, and with the same
reference-safe bound. Applying one common exact-return native C would
send these states to $`C|0\rangle`$ and $`C|\phi_1\rangle`$ without
changing that bound. C and its gates are not supplied by this emitter.
The existing [common-coarse proof](STATE_ONLY_COMPILER.md#7-a-complex-native-reference-from-the-same-coarse-word)
charges that additional circuit and the subsequent protocol separately.

The [native tests](../tests/test_native_branched_residual_state.py)
check the complete 2048-by-256 initialized isometry at q=5, including
both branch values and every dirty-input column. Native table sectors
are compared with independent source algebra; coherent reference inputs
are also propagated through the flattened word and actual inverse.
Ideal counterexamples detect branch phase changes and a reflection
that incorrectly tests the branch. Small exact normalized complex
fixtures lie inside the actual 1/64 residual neighborhood; separate
phase and boundary witnesses need not. Fine-precision checks use exact
rational certificates and literal gate counts instead of large
statevectors. Finite numerical evidence supplements the analytic bound.

This completes the bounded coherent residual selector for one system
qubit. The next task is a bounded native common-coarse QBP integration,
charging forward C, the controlled observable, and inverse C in the
magnitude stream while retaining both decoder streams. General lookup
and the full fine state compiler remain implementation tasks. No new
asymptotic result, end-to-end advantage, or complete-frame closure is
claimed here.
