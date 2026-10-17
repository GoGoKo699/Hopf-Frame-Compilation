# Gaussian inverse extraction: an exact identity and its paid boundary

[Gauss-sum encoding](GAUSS_SUM_QUERY_ENCODING.md) ·
[Retained query R91](PARALLEL_BILINEAR_WORD_QUERY.md)

The affine Gaussian encoding admits exact scalar extraction on arbitrary
dirty work if an **in-place affine linear map** is supplied. This avoids
the small normalized trace, but computing that map is a separate native
task. Within the affine-symplectic round architecture defined below,
its worst-case round complexity is exactly the address length. The
round obstruction applies even with initialized helpers that must return;
it is not a T-depth lower bound for general Clifford+T circuits, and it
does not obstruct the scalar table query itself.

## 1. Exact dirty Weyl extraction

Use the encoding's notation: an a-bit address x specifies an invertible
binary matrix $`A(x)=I+N(x)`$ of order $`m=O(\sqrt Q)`$, where
$`Q=2^a`$, and fixed vectors $`b=e_s,c=e_t`$ satisfy

```math
c^{\mathsf T}A(x)^{-1}b=f(x).
```

On an arbitrary m-bit work word define the in-place permutation
$`C_A|z\rangle=|A(x)z\rangle`$. If its native implementation
preserves x and returns its additional work, then

```math
M_b=C_A^\dagger X(b)C_A=X(A(x)^{-1}b).
```

Here $`X(v)|z\rangle=|z\oplus v\rangle`$ and
$`Z(c)|z\rangle=(-1)^{c^{\mathsf T}z}|z\rangle`$.
The chronological word

```math
M_b,\quad Z(c),\quad M_b^\dagger,\quad Z(c)
```

therefore equals $`(-1)^{f(x)}I`$ on the entire work space.
Its phase is literal, independent of the initial work word. It preserves
coherent addresses and arbitrary external references. In each M call,
the chronological circuit is $`C_A,X(b),C_A^\dagger`$; the inverse
is the actual reversed native word, not a new ideal oracle.

If one C implementation has simultaneous resources
$`(T_C,G_C,D_C,w_C)`$, excluding the address and m-bit work word,
the extraction has

```math
T\le4T_C,\qquad
G\le4G_C+2\mathrm{wt}(b)+4\mathrm{wt}(c),\qquad
D_T\le4D_C.
```

Its simultaneous work is m plus $`w_C`$, with the same clean
reservation as C. The four calls are sequential and reuse only a
returned pool. To control the scalar by a separate arbitrary bit h,
replace each $`Z(c)`$ by its product of $`\mathrm{wt}(c)`$
CZ gates controlled by h. The displayed scalar-query ledger expands
each Z as $`S^2`$; for the controlled query expand each CZ as H-CNOT-H,
replacing the last Clifford contribution by $`6\mathrm{wt}(c)`$.
No non-Clifford control of C is required.
This gives exact identity on the entire h-zero sector. Hadamards on an
arbitrary output h convert the controlled phase into the complete XOR
query $`h\mapsto h\oplus f(x)`$, adding two Clifford gates.

If the supplied C instead has full initialized-isometry error at most
$`\delta`$ relative to $`C_A`$ and returned work, its actual inverse
has the same bound on the ideal returned subspace. A four-term unitary
telescoping estimate gives extraction error at most $`4\delta`$.
All intervening ideal Pauli gates preserve that subspace. This estimate
includes leakage and reference correlations; no trace-norm or average
work-state guarantee is substituted for operator error.

## 2. The cheap Gaussian shear does not supply the inverse map

The bilinear phase $`(-1)^{y^{\mathsf T}A(x)z}`$, conjugated by
Hadamards on y, gives the doubled shear

```math
S_A:(y,z)\longmapsto(y\oplus A(x)z,z).
```

It has the encoding's constant native T-depth, $`O(m)`$ T-count,
$`O(Q)`$ elementary Clifford count, and $`O(m)`$ returned dirty
helpers. It is not $`C_A`$. Conjugating $`X_z(b)`$ by this shear
gives $`X_y(A b)X_z(b)`$. A Pauli commutator with $`Z_y(c)`$
therefore extracts $`(-1)^{c^{\mathsf T}A b}`$, rather than the
required inverse entry.

The same distinction appears in completing the square. A translation
$`(y,z)\mapsto(y\oplus u,z\oplus v)`$ changes the bilinear
form by

```math
y^{\mathsf T}A v+u^{\mathsf T}A z+u^{\mathsf T}A v.
```

Removing the encoding's two linear terms by such a translation
uniquely requires $`A v=b`$ and $`A^{\mathsf T}u=c`$.
The resulting constant is precisely $`c^{\mathsf T}A^{-1}b`$.
Thus a shift-only square-completion repair also owes the inverse
displacement. This uniqueness statement concerns translations; it does
not classify all possible linear, Clifford, or non-Clifford changes of
coordinates.

## 3. A fully native baseline for the owed map

The directed graph contains a left prefix tree of depth
$`p=\lfloor a/2\rfloor`$, the fixed table cross edges, and a reversed
right prefix tree of depth $`q=\lceil a/2\rceil`$.
In topological vertex order, apply every edge update

```math
w_v\longmapsto w_v\oplus N(x)_{v,u}w_u.
```

All incoming edges of a vertex are applied after its predecessors have
their final values. Consequently the final word satisfies
$`w'=w+Nw'`$, so this edge word is $`C_{A^{-1}}`$; its actual
reverse is $`C_A`$. This assertion holds on arbitrary work words,
not just a source impulse.

The constant cross-edge block uses at most Q CNOTs. Each variable tree
edge uses a Toffoli with one address literal and two work wires; a
complemented literal adds a CNOT. At one tree level all edges run
together in the eight-T-layer private-helper phase schedule of R91
Section 1, conjugated by target Hadamards. Multiple edges into a shared
target commute, and that target receives just the enclosing Hadamards;
all phase T gates are on distinct private dirty helpers. Sources and
targets belong to different levels. There are $`p+q=a`$ such levels,
and $`O(m)`$ variable edges in total. Every helper returns after its
completed level; the address only enters as a Clifford control in this
substitution. Thus both C maps, and the extraction above, have

```math
T=O(\sqrt Q),\qquad G=O(Q),\qquad
b=O(\sqrt Q),\qquad D_T=O(a),
```

with no initialized work. The live m-bit graph word and largest private
edge-helper pool are disjoint and included simultaneously. This is a
valid native recovery of the encoded sign, but its depth is worse than
R91's $`O(\log_2(a+2))`$ at the same count and width orders.

## 4. A tight affine-symplectic round obstruction

An **affine-symplectic round** is an address-preserving operation whose
restriction to each fixed computational address x is a Clifford on all
physical work wires, and whose binary symplectic matrix has entries of
algebraic degree at most one in x. Between rounds allow arbitrary
address-independent Clifford circuits on the physical work, leaving the
address untouched. The matrices and the number of
wires may depend on the prescribed table. Pauli signs and scalar phases
are unrestricted by this definition.

The Gaussian phase/shear blocks above are such rounds. Parallel blocks
on disjoint work spaces are allowed in a round; their combined matrix
is a direct sum of affine matrices. A collection of overlapping blocks
counts as one round only when its *completed* symplectic matrix is
affine. Arbitrary products of affine matrices are not silently declared
one round.

After r rounds, every binary symplectic entry has algebraic degree at
most r. Indeed the total matrix is a product of r affine matrices and
constant matrices. Expanding an entry gives products of at most r
affine polynomials; reducing $`x_i^2=x_i`$ cannot increase degree.
The inverse action of the target in-place permutation is visible in its
Pauli conjugations:

```math
C_A X(v) C_A^\dagger=X(A v),\qquad
C_A Z(w) C_A^\dagger=Z(A^{-\mathsf T}w).
```

Choose the table to be one only at the all-ones left and right leaves.
Then the prescribed inverse entry is
$`f(x)=\prod_{i=1}^a x_i`$, whose algebraic normal form has degree a.
It is an entry of $`A^{-\mathsf T}`$ with row and column swapped.
Therefore any exact r-round realization of $`C_A`$, with arbitrary
additional returned dirty work, satisfies $`r\ge a`$.

The bound is stable under full operator approximation. Distinct Pauli
labels have operator-norm separation at least $`\sqrt2`$. If two
Clifford unitaries are within error $`\varepsilon`$, their conjugations
of any unit Pauli are within $`2\varepsilon`$. For
$`\varepsilon\lt1/\sqrt2`$, the symplectic labels must therefore agree.
The same degree argument applies. A common or address-dependent scalar
phase in the comparison does not change the conjugations.

### Initialized helpers do not evade this round obstruction

The preceding approximation statement also permits any finite number
of clean work qubits, initially and finally zero. Let J embed arbitrary
nonclean work with those zeros. Suppose each fixed-address physical
circuit V is Clifford and

```math
\|VJ-J(C_A\otimes I_{\rm dirty})\|\le\varepsilon\lt1/\sqrt2.
```

For every output clean qubit j, the expectation of $`Z_j`$ on the
output obtained from maximally mixed nonclean input is at least
$`1-2\varepsilon^2\gt0`$. But $`V^\dagger Z_jV`$ is a signed
Pauli. Its expectation on zero clean qubits and maximally mixed
nonclean work belongs to $`\{0,+1,-1\}`$. It is therefore +1.
Every output clean stabilizer is exactly satisfied, so V preserves the
entire clean code exactly.

The induced unitary on the nonclean work is Clifford. In the physical
image of any nonclean Pauli, all clean factors are I or Z; deleting
those factors gives its logical Pauli label. Hence its logical
symplectic matrix is the corresponding row/column submatrix of the
physical one, whose entries still have degree at most r. The same
$`2\varepsilon\lt\sqrt2`$ argument fixes this logical matrix to that of
$`C_A\otimes I_{\rm dirty}`$. Thus $`r\ge a`$ persists with
initialized helpers and the full initialized-isometry error criterion.

The native edge baseline uses exactly a variable-affine rounds,
interspersed with fixed Clifford cross edges. Thus the worst-case bound
is tight for this round model. It does **not** give an unrestricted
native T-depth lower bound: rounds may have other syntheses, and a
general Clifford+T circuit need not be Clifford at fixed x or admit
affine round matrices. Nor can this symplectic argument lower-bound
the desired scalar $`(-1)^{f(x)}I`$, whose symplectic matrix is the
identity for every f. A different scalar extractor that avoids the
in-place inverse map remains possible.

The subsequent [scalar-loop obstruction](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
addresses that missing central phase with a different invariant. For
actual affine quadratic phase rounds and fixed Clifford interlayers,
zero-column return already forces $`\deg f\le4r+2`$. It includes
initialized helpers and a fixed error threshold. It does not hold merely
from specifying affine symplectic matrices with arbitrary scalar lifts,
and it remains distinct from a general native T-depth lower bound.

## 5. Evidence and resumption question

[Exact bounded fixtures](../../tests/test_gaussian_inverse_extraction.py)
check the full DAG map and its actual reverse on every work basis word
for every two-bit table and address, the literal Weyl phase on those
same arbitrary inputs, a failure of substituting the cheap shear,
and the maximal-degree inverse entry for increasing address sizes.
They support the logical edge-word interfaces. The native private-helper
substitution uses the retained [R91 finite native checks](../../tests/test_parallel_bilinear_word_query.py).
The symbolic arguments above prove the arbitrary-size resource and round
statements; neither finite suite supplies a general circuit lower bound.

The continuing question is whether the scalar table phase can be
extracted by a native mechanism that avoids implementing the in-place
inverse map, at $`O(\sqrt Q)`$ T-count, $`O(Q)`$ elementary Clifford
count, $`O(\sqrt Q)`$ dirty width and
$`o(\log_2(a+2))`$ T-depth, with all work returned. The Gaussian
round obstruction closes this particular inverse-based repair, not
that query question or the general Hopf depth objective.
