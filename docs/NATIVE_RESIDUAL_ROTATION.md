# Certified native residual rotations

[State-based QBP theorem](STATE_BASED_QBP_THEOREM.md) · [Residual coefficients](RESIDUAL_TABLE_PREPROCESSING.md) · [Verification](VERIFICATION.md)

The production modules now connect certified residual coefficients to one
literal, unaddressed Clifford+T completion. This implements the paired-source
identities of [the one-clean compiler](ONE_CLEAN_COMPILER.md), including its
borrowed-signal extension, and the algebraic preprocessing of
[the residual table](RESIDUAL_TABLE_PREPROCESSING.md). It adds executable
evidence for those constructions, with no new asymptotic compiler claim.

## 1. Input and output contract

```python
from fractions import Fraction
from compiler_robust_hopf.native_residual_rotation import emit_residual_row

row = emit_residual_row(Fraction(1, 2), Fraction(1, 4), q=80)
assert row.rotations[0] is row.rotations[2]
assert row.operator_error_bound < Fraction(130, 1 << 80)
```

The caller supplies dyadic real and imaginary parts of a coefficient
approximation, together with the promise

```math
|z|\le1,\qquad |z_0-z|\le2^{-2q-20},\qquad q\ge5.
```

Input checks reject floats and obviously inconsistent approximations;
they cannot establish the caller's external evaluation promise. The target is

```math
U(z)=\begin{pmatrix}z&-\sqrt{1-|z|^2}\\
\sqrt{1-|z|^2}&\overline z\end{pmatrix}.
```

`row.gates` is an immutable chronological tuple of `(name, wires)` pairs.
The elementary gate names are:

```text
X Z H S SDG T TDG CX
```

Qubit zero is the low integer bit. The layout is explicit:

| Wires | Role | Input contract |
|---|---|---|
| 0 through q | Precision core | Arbitrary borrowed state |
| q+1 | Synthesis signal | Arbitrary borrowed state |
| q+2 | Rotation target | Arbitrary data state |

Thus the row has q+3 total wires and q+2 borrowed wires excluding its
target. It needs no clean qubit, predicate helper, prepared source state,
measurement, or reset. The work return is approximate and included in
the full-operator error, which also covers entanglement with any external
reference. At q=P+10 this component uses P+12 dirty wires. An addressed
construction separately reserves its predicate helper; this does not
reduce the theorem's full workspace reservation.

## 2. Exact coefficient programming

[`rotation_programming.py`](../compiler_robust_hopf/rotation_programming.py)
accepts two `CoefficientInterval` objects through `program_rotation`.
Each interval has exact rational endpoints and width at most
$`2^{-q-20}`$. Their rectangle must intersect the unit circle. The caller
still promises that they enclose the intended common cosine/sine pair.
There is no inverse-angle reconstruction.

The programmer encloses $`\beta=(\sqrt5-1)/4`$ by integer square root and
forms rational intervals for

```math
p_+=\beta(3\cos\theta/4-\sin\theta/2),\qquad
p_-=\beta(\cos\theta/4+\sin\theta/2).
```

The sign of the midpoint of the first interval selects the head, with a
positive tie. It then estimates the tail means
$`u=2p_+-\sigma/2`$ and $`v=4p_-`$, clips to the known range, and rounds
to the geometric sign grid. Rounding ties choose the larger integer
index, hence the lower mean. The repeated last geometric weight encodes
the negative endpoint exactly. No floating-point rounding or uncertain
sign test is used.

The returned immutable program records the Majorana flip bits, exact
dyadic moments s and p, the input intervals, and rational error bounds.
The compressed operator is $`sI+pXZ`$. An integer-square-root enclosure
of the squared coefficient error certifies

```math
\left\|sI+pXZ-\beta R_y(\theta)\right\|
\le\delta_{\rm block}\lt\frac52\,2^{-q}.
```

The programmer checks the amplification regime as well. The finite-data
head decision and tail errors fit the budgets in the existing proof:
an incorrect head sign can only occur near zero, where the tail remains
in range. `emit_rotation(program, axis="y")` also recomputes the program
to reject altered signs or certificates. The supported alternative axis
is `"z"`.

## 3. Literal native word and error

[`native_residual_rotation.py`](../compiler_robust_hopf/native_residual_rotation.py)
builds the fixed paired source, literal Pauli masks, and their actual
reversed adjoints. The source appears through conjugated Pauli operators;
it is never supplied as an initialized state. The five-call amplification
uses four physical signal Z gates. Replacing the four ideal reflections
by these Z gates contributes no remaining scalar phase.

The existing amplification bound gives initialized-signal error less
than $`12\delta_{\rm block}`$. Exact commutation with signal X extends
this to arbitrary signal input with a factor $`\sqrt2`$. Therefore each
emitted rotation has full-operator error less than $`43\,2^{-q}`$.
For Rz, the chronological target gates are H, S, the Ry word, SDG, H;
they implement the literal Clifford conjugation with the required sign.

The residual helper supplies a shared outer half-phase and a middle real
rotation. The emitter uses the same outer program and word object twice
in chronological Rz, Ry, Rz order. Together with the helper's completion
bound, the row certificate is

```math
\|\widetilde U-U(z)\otimes I_{\rm work}\|
\le(129+1/16)2^{-q}\lt130\,2^{-q}.
```

No global or branch-relative phase is removed when forming or testing
the word. Production construction uses rational arithmetic and O(q)
gate storage, without dense matrices. Each source loader has 2q T/TDG
gates. The literal unsimplified amplification has 180q T/TDG gates per
rotation, hence 540q per residual row, with O(q) Clifford gates too.
These constants describe this implementation; they are not an optimized
practical synthesis estimate.

## 4. Evidence and next boundary

The [programming tests](../tests/test_rotation_programming.py) check exact
moments, rational certificates beyond floating-point precision, digit
and head boundaries, interval contracts, and endpoint encoding. The
[native tests](../tests/test_native_residual_rotation.py) inspect emitted
elementary gates on all input columns at small precision, including
literal phases, borrowed-signal symmetry, inverses, and the residual
composition. Fine-precision amplification is also checked in the small
Clifford-algebra representation; that check does not simulate a large
precision register.

One unaddressed row is now implemented. The general addressed residual
table, its dirty lookup/predicate schedules, and the complete fine state
compiler remain implementation work. A next bounded pass should expose
one small addressed two-row table with exact inactive-sector behavior
and a complete workspace/error contract before integrating state
amplification. The present row does not establish the full theorem's
count or depth schedules by itself.
