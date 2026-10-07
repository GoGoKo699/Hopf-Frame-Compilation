# Certified native residual rows and tables

[State-based QBP theorem](STATE_BASED_QBP_THEOREM.md) · [Residual coefficients](RESIDUAL_TABLE_PREPROCESSING.md) · [Verification](../../docs/VERIFICATION.md)

The production modules connect certified residual coefficients to literal
Clifford+T completions: one unaddressed row and tables of two or four rows
selected by one or two address bits and one enable literal. They implement the paired-source
identities of [the one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md), including its
borrowed-signal extension, and the algebraic preprocessing of
[the residual table](RESIDUAL_TABLE_PREPROCESSING.md). This adds executable
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
reference. At q=P+10 this component uses P+12 dirty wires. General addressed
constructions separately reserve lookup and predicate work; this does not
reduce the theorem's full workspace reservation. Section 5's bounded
one-address, one-enable case needs no additional helper.

## 2. Exact coefficient programming

[`rotation_programming.py`](../../compiler_robust_hopf/rotation_programming.py)
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

[`native_residual_rotation.py`](../../compiler_robust_hopf/native_residual_rotation.py)
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

The [programming tests](../../tests/test_rotation_programming.py) check exact
moments, rational certificates beyond floating-point precision, digit
and head boundaries, interval contracts, and endpoint encoding. The
[native tests](../../tests/test_native_residual_rotation.py) inspect emitted
elementary gates on all input columns at small precision, including
literal phases, borrowed-signal symmetry, inverses, and the residual
composition. Fine-precision amplification is also checked in the small
Clifford-algebra representation; that check does not simulate a large
precision register.

The unaddressed row and Sections 5–6's bounded two- and four-row extensions
are implemented. General addressed tables, larger dirty lookup/predicate
schedules, and the complete fine state compiler remain implementation
work. These components do not establish the full theorem's count or
depth schedules by themselves.

## 5. Two-row addressed tables with one enable literal

[`native_residual_table.py`](../../compiler_robust_hopf/native_residual_table.py)
emits exactly two residual rows at a common q, with one unchanged address
and an enable literal whose active value is zero or one. The per-row
dyadic input promise is the same as Section 1. For example:

```python
from fractions import Fraction
from compiler_robust_hopf.native_residual_table import emit_residual_table

table = emit_residual_table(((0, 0), (1, 0)), q=80, enable_value=0)
assert table.rotations[0] is table.rotations[2]
assert table.operator_error_bound < Fraction(130, 1 << 80)
```

The logical target, address, and enable are arbitrary occupied data,
including coherent superpositions and reference correlations. There is
no initialized predicate or helper:

| Wires | Role | Return contract |
|---|---|---|
| 0 through q | Precision core | Approximate identity, included in the error |
| q+1 | Synthesis signal | Approximate identity, included in the error |
| q+2 | Rotation target | Selected residual action |
| q+3 | Address | Preserved exactly |
| q+4 | Enable | Preserved exactly |

This is q+5 total wires, of which q+2 are borrowed work. For active enable
value v, the ideal action in register order enable, address, target, work is

```math
\mathcal T=
\Pi_v\otimes\left(\sum_{x=0}^1|x\rangle\langle x|\otimes U(z_x)\right)
\otimes I_{\rm work}+(I-\Pi_v)\otimes I,
\qquad \Pi_v=|v\rangle\langle v|.
```

### Exact address selection and inactive sectors

Each Pauli mask table consists of the row-zero mask and address-controlled
differences between the two rows, with all X factors before all Z factors.
Controlled X is CNOT; controlled Z is its target-Hadamard conjugation.
The address is never a gate target. These Clifford words equal the
literal mask in each address sector, including its scalar phase. Every
inverse reverses the actual emitted word.

The enable control is added only to the three Pauli centers inside each
scalar source word. The first center is a CNOT controlled by enable;
the other two are literal seven-T Toffolis controlled by enable and
signal. The loaders, masks, and routing remain
unconditional. When enable is inactive, each center is identity; loader
and mask inverse pairs cancel, followed by the signal Hadamard pair.
Both routed scalar words are then identity, and the four amplification
Z gates multiply to identity. The Rz basis changes cancel as well.
Thus the inactive sector is exactly identity on every work input, not
merely close to identity or correct on one prepared state. An outer
enable-X pair supplies the active-zero convention without a new flag.

This is the exact source-center conditioning of
[the one-clean proof §5](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger).
For one enable literal the largest center has only two controls, so the
Toffoli decomposition needs no helper. Larger predicates still need
their separately charged construction; no general predicate compiler
is implemented here.

### Error, cost, and verification

`emit_rotation_table` revalidates the two programs and requires the same
precision. On each active address sector it equals the corresponding
unaddressed emitted rotation, with full-operator error below
$`43\,2^{-q}`$. The address direct sum takes the maximum of row errors;
the inactive error is exactly zero. The residual table reuses one outer
Rz table object twice, with each row's consistently chosen half-phase.
Consequently

```math
\|\widetilde{\mathcal T}-\mathcal T\|
\le\max_{x=0,1}\left(\epsilon_{{\rm completion},x}
+\sum_{j=1}^3\epsilon_{x,j}\right)
=(129+1/16)2^{-q}\lt130\,2^{-q}.
```

There is no additional address approximation and no sum over the two
rows. The literal word has 180q+210 T/TDG gates per rotation table and
540q+630 per residual table. The extra 210 per rotation comes from
30 seven-T Toffolis in the unchanged five-call construction. Mask-table
gates are Clifford. Gate storage and total gate count remain O(q) for
this fixed table size; these are unsimplified counts, not optima.

The [table tests](../../tests/test_native_residual_table.py) check the literal
Toffoli, address-mask phases, active row words, and inactive identity.
At q=5 they propagate every core/signal/target input in each fixed address
and enable sector, retaining phases on the fixed control wires. This
checks the full ten-wire operator through 256-dimensional blocks without
a 1024-dimensional dense simulation. It also tests coherent relative
phases, active-zero selection, actual inverses, shared outer tables,
input contracts, and exact rational high-precision resource certificates.

The [bounded residual state emitter](NATIVE_RESIDUAL_STATE.md) now composes
two such tables into a one-system-qubit preparation, with two explicitly
initialized compiler flags and the actual inverse used in amplification.
General table sizing and the full count/depth schedules remain distinct
implementation tasks.

## 6. Four rows with two address bits and no additional helper

[`native_residual_lookup.py`](../../compiler_robust_hopf/native_residual_lookup.py)
adds `emit_four_row_rotation_table` and `emit_four_row_residual_table`.
Their input promises, elementary gate alphabet, active-zero convention,
and full-operator error contract are those of Section 5. The row index is
$`x=x_0+2x_1`$, with the low address bit listed first. Existing two-row
and one-qubit state APIs are unchanged.

```python
from fractions import Fraction
from compiler_robust_hopf.native_residual_lookup import emit_four_row_residual_table

table = emit_four_row_residual_table(
    ((0, 0), (1, 0), (0, 1), (Fraction(-3, 4), Fraction(1, 2))), q=80)
assert table.rotations[0] is table.rotations[2]
assert table.operator_error_bound < Fraction(130, 1 << 80)
```

| Wires | Role | Return contract |
|---|---|---|
| 0 through q | Precision core | Approximate identity, included in the error |
| q+1 | Synthesis signal | Approximate identity, included in the error |
| q+2 | Rotation target | Selected residual action |
| q+3, q+4 | Address x0, x1 | Preserved exactly by the complete word |
| q+5 | Enable | Preserved exactly by the complete word |

The total is q+6 wires, with the same q+2 borrowed work wires. No clean
input, lookup bank, or extra borrowed helper is required. Target and
controls may be occupied and arbitrarily entangled with the work or a
reference. The ideal action is the direct sum in Section 5 with four
active rows instead of two.

### Exact quadratic mask selection

For each core wire and each axis X or Z, let $`f_j`$ be that axis's
mask bit in row j. Its Boolean polynomial is

```math
f(x_0,x_1)=f_0\oplus(f_0\oplus f_1)x_0
\oplus(f_0\oplus f_2)x_1
\oplus(f_0\oplus f_1\oplus f_2\oplus f_3)x_0x_1.
```

The constant and linear terms use the same Pauli and controlled-Pauli
gates as Section 5. For a fixed axis, collect all wires whose quadratic
coefficient is one into a support S. All terms on that axis commute.
If S is nonempty, choose its smallest core wire p as a pivot and put

```math
C_X=\prod_{j\in S\setminus\{p\}}\mathrm{CX}_{p\to j},
\qquad
C_Z=\prod_{j\in S\setminus\{p\}}\mathrm{CX}_{j\to p}.
```

Each fanout or fanin is its own inverse and obeys

```math
C_X X_p C_X=\prod_{j\in S}X_j,\qquad
C_Z Z_p C_Z=\prod_{j\in S}Z_j.
```

Conjugating a single literal seven-T Toffoli on controls x0,x1 and
target p by $`C_X`$ therefore selects the entire X support. For Z,
conjugate the Toffoli target by H to obtain CCZ, then use $`C_Z`$.
This costs seven T/TDG gates per nonempty axis support, irrespective of
its size; all remaining fanout/fanin gates are charged Cliffords.
The pivot is existing arbitrary core data, not a new initialized wire.

Emit every X factor before every Z factor. This matches the chronological
order of `_mask` in each row, including any Pauli-product scalar phase.
Swapping the two axis blocks is not phase-neutral when their supports
overlap oddly. Each inverse reverses the actual emitted word.

The native seven-T Toffoli temporarily changes its second address control
through a CNOT between controls. Address preservation is an identity of
the complete Toffoli, mask, and table; it is not a claim that every
elementary gate leaves both address bits untouched. The mask is exactly
the required core Pauli in each complete address sector, with identity
on target and signal. It can therefore replace the two-row mask inside
the existing scalar source and amplified rotation without changing
their proofs.

Enable still conditions only the source centers. Inactive centers make
the literal loader and mask inverse pairs cancel on all inputs, including
all the quadratic gates. Hence every inactive sector is exactly identity;
active zero coefficients still implement U(0). The extra address bit
does not enlarge the enable predicate or require a hidden clean flag.

### Exact ledger and evidence boundary

Let k be the number of nonempty quadratic supports among X and Z in a
rotation's four-row mask. Then $`k\in\{0,1,2\}`$, and the mask uses
7k T/TDG gates. It or its actual inverse appears ten times in that
rotation word, so

```math
T_{\rm rotation}=180q+210+70k.
```

For a residual completion, let $`k_z`$ and $`k_y`$ refer to its outer
and middle rotation tables. The identical outer object is used twice:

```math
T_{\rm residual}=540q+630+70(2k_z+k_y)\le540q+1050.
```

All gates, including mask selection and source-center controls, are
counted. Storage and total gate count are $`O(q)`$ for this fixed table
size. These are literal unsimplified counts, not optimality claims.
The exact lookup adds no approximation error. Address direct sums still
take maximum row errors, followed by the three Euler-factor error sum,
giving full-operator error below $`130\,2^{-q}`$ on arbitrary work and
reference inputs.

The [lookup tests](../../tests/test_native_residual_lookup.py) compare small
native masks with independently specified literal Pauli rows and verify
all eight address/enable sectors at q=5. Each sector retains all 256
target/signal/core input columns. The sector evaluator tracks temporary
classical address changes and their diagonal phases; it rejects mixing
with free quantum inputs or failing to restore the controls. Independent
Majorana matrices check the active rotation words, and inactive words
are compared with identity without phase alignment. Fine-precision
certificates and emitted resource counts require no large statevector.

This completes the second-address lookup component. General table sizes,
larger predicates, and the complete fine state-preparation schedule
remain implementation work. The [two-system-qubit residual composition](NATIVE_RESIDUAL_STATE.md#5-two-system-qubit-preparation-and-a-returned-core-helper)
now uses these tables with two clean compiler flags, an enlarged charged
initial reflection, and amplification with the actual inverse. The
[coherent residual selector](NATIVE_RESIDUAL_BRANCH.md) instead uses the
second address as an arbitrary protocol branch, retaining its relative
phase and excluding it from the initial reflection. No complete-frame
endpoint or end-to-end advantage follows.
