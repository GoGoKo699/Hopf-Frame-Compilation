# Revision: non-Pauli clean codes and the cost of changing representations

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Partial-diagonal theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md) ·
[Code decoding](NONPAULI_CODE_DECODING.md) ·
[Module query barrier](SYMPLECTIC_MODULE_QUERY_BARRIER.md) ·
[Previous checkpoint](OCCUPIED_TRANSPORT_SESSION_RESULT.md)

This continuation began from verified main `55ce8b0` (PR #187), on
19 October 2026 in Asia/Shanghai. It closes two further three-layer
subcases and identifies the native operation that prevents a general
extension. The general complete-frame resource frontier is unchanged.

## 1. Destination, strongest results, components, and scoped limits

The destination remains optimal general T-depth for the prescribed
complete Hopf differential frame. All logical columns, literal phases,
dirty inputs and references, clean leakage, work return and actual
inverses belong to the same resource contract.

| Role | Current statement |
| --- | --- |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, and $`G=O(NL)`$ for every integer $`n\ge1,L\ge6`$ |
| Literal allocation and error | One clean qubit, real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$, full initialized-isometry error at most $`2^{-L}`$ |
| Selected manuscript claims | R79 remains Result D; the count package is settled under its hypotheses; final manuscript drafting is separate |
| Unrestricted component lower results | At most one clean and two layers force XOR degree at most two; zero clean and three layers force XOR degree at most four, with their retained uniform error floors |
| New restricted three-layer result | Either actual interior clean image being Pauli forces sign degree at most $`2^{26}`$ and XOR degree at most $`2^{26}-1`$. Both Pauli images retain the stronger 2048/2047 bounds |
| Further restricted closure | Two non-Pauli interior images, each of clean-active rank one, admit a three-band Pauli-boundary representation and hence the stronger 2048/2047 bounds |
| Useful components | Full-space partial diagonality gives a finite low-degree matrix-feature description; a one-band non-Pauli code has an explicit paid decoder; a scoped zero-pair defect family has a common central Clifford determinant coboundary |
| Representation limits | Cheap classical decoder controls do not imply cheap native conjugation. A single elementary-abelian symplectic module has an exponential degree/carrier and degree/source-count barrier, but a two-family construction escapes that hypothesis |
| Still open | Uniform unrestricted one-clean three-layer degree, unrestricted one-clean four-layer rigidity, optimal general depth and dimension-growing matching lower bounds |

At fixed accuracy and square-root dirty width, the known general upper
depth remains logarithmic in n, while the unrestricted lower bound
remains constant. None of the new components alone changes that interval.

## 2. First unpaid operation: decoding a non-Pauli clean code

A native one-band clean image can be put into the form

```math
S=DQ_eD^\dagger,\qquad
D=\sum_z|z\rangle\langle z|\otimes S_t^{r(z)},\qquad
Q_e=T_t^eX_tT_t^{-e}.
```

Only four equatorial directions occur, controlled by two Boolean features
of degrees one and two. This compact description is a component result,
not a free Clifford change of basis. The
[literal decoder and inverse](NONPAULI_CODE_DECODING.md)
are explicitly charged. For k active axes, the controlled operation D
has at most two native bands, $`2k-1-e`$ T gates and $`2(k-1)`$
CNOTs in its canonical coordinates, with no added workspace; coordinate
Cliffords have separate costs. The complete code decoder has a different
central phase and cannot silently replace D in an isometry identity.

The first unsupported inference would be to decode the code and apply
native Pauli-axis matching while keeping the original three bands.
Conjugating a native Pauli by this decoder can give a Clifford with
unbounded symplectic displacement rank. The explicit conditional-S
family has only degree-one/two controls but displacement rank m or
$`m+1`$ on m controls. Low-degree controls therefore do not by themselves
give a bounded quantum core. This refutes that proposed reduction; it
does not exhibit a hard Boolean target or obstruct unrestricted circuits.

## 3. What can be closed without paying for a new decoder?

The [partial-diagonal theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md)
starts with an actual full-space three-band unitary preserving a
computational-Z subgroup of codimension k. Every exact k-qubit fiber is
determined by a bounded number of bounded-degree Boolean features,
depending only on k. The proof retains the native middle pairs, builds
an explicit Clifford correction on each active Pauli eigensector, and
counts the resulting finite matrix types modulo eight. The factors need
only commute on simultaneously attained fibers; their entire abstract
matrix alphabet need not commute. Silent endpoint factors and the
residual signed Clifford phase are included.

For k=1 the sign-degree bound is 37,748,829. One additional Pauli T
factor commuting with the preserved subgroup adds one affine feature,
still below $`2^{26}`$.

If either of the two actual interior clean images is Pauli, the retained
two-band clean-removal identity supplies a full-space extension of the
other transition. The final projected native band has at most one axis
relation, so it splits into one independent band and one explicit Pauli
T factor. The three-band prefix now fixes a computational-Z subgroup of
codimension at most one. This is a proved full-space property, rather
than an inference from clean-zero return, and the new theorem applies.
The reverse case uses the actual inverse word. The existing R74
rounding result transfers the bound to approximation error below
$`2^{-530690}`$ within this subcase.

A second closure covers two non-Pauli interior images with active ranks
one and one. Native matching forces the active middle axes to cancel
the endpoint active rotations into a literal Clifford. Repacking the
remaining three silent bands gives Pauli intermediate clean images, so
the earlier 2048/2047 bounds apply. The exact repacking does not assume
a shared clean-flip Pauli in the original split.

Thus an unbounded-degree counterfamily must have both interior clean
images non-Pauli, with at least one active rank exceeding one. In that
remaining branch, a scoped zero-pair, zero-probe-character reduction
also shows that the unused-sector defects commute and their phases
share one determinant coboundary from the actual middle band. Those
phases are central diagonal Clifford phases in intrinsic coordinates.
Identifying those coordinates with the target computational variables,
or implementing their gauge for free, remains unsupported. The
exceptional rank-two stabilizer and nonzero probe characters are not
covered by this residual reduction.

## 4. Constructive revision: a broader barrier, and its explicit escape

The earlier algebra-quotient barrier concerned a particular coefficient
representation. The [symplectic-module theorem](SYMPLECTIC_MODULE_QUERY_BARRIER.md)
now treats any elementary-abelian representation on a binary symplectic
space of dimension 2w. For a degree-d readout

```math
h(x)=\langle u,(S_x-I)v\rangle,
```

it proves $`w\ge2^{d-1}`$. If the represented Clifford family is
actually supplied by $`C_x=FP_xF^\dagger`$ for a linear physical
Pauli family, it also forces $`T(F)\ge2^{d-1}`$, including arbitrary
returned helpers. The source-count statement does not assume that F
lies in the third Clifford-hierarchy level. Its full-operator
approximation extension retains its separately stated error hypothesis.

For a full-degree a-bit table with $`Q=2^a`$ rows, eight separately
paid source occurrences therefore use at least $`4Q`$ T gates.
This architecture cannot meet R98's square-root count budget. The bound
does not apply to a fully fused synthesis, a different readout, or an
arbitrary native circuit.

There is an explicit escape from the module hypothesis: use two
generally noncommuting families on two half-address banks. A table is
programmed by a literal CNOT network of at most Q gates, and a complete
commutator extracts its XOR while returning every bank and helper.
Balanced carrier width is $`O(\sqrt Q)`$ even for full-degree tables.
This demonstrates why the module lower bound is not unrestricted.

It does not improve R91. The complete selector transports and their
actual inverses still cost depth $`O(a)`$ in the direct implementation;
the fair retained R91 interface already has square-root T count,
linear Clifford count and logarithmic address depth. The first missing
native operation is a sublogarithmic-depth complete selector or probe
transport within the same count and width budgets, followed by the
joint word-query schedule required by R98. Merely changing the carrier
or rebalancing the old recovery retains the old depth cost.

## 5. Next question and verification

The next lower question is whether joint native probe transport bounds
the completed Boolean phase when both interior clean images are
non-Pauli, after the paid pair-center reduction. The useful new residual
is a common central Clifford determinant phase; the first gap is its
connection to the actual computational input and the remaining native
bands. The constructive alternative must improve the complete selector
interface or change the full-frame mechanism while retaining every
resource budget. Neither route is abandoned because one representation
retains its old cost.

Ten new exact finite tests cover the
[code decoder](../../tests/test_nonpauli_code_decoding.py),
[partial-diagonal pair identities](../../tests/test_partial_diagonal_three_layer_features.py),
and [module barrier and escaped query](../../tests/test_symplectic_module_query_barrier.py).
They retain all small physical columns, native phases, full unused
blocks, returned arbitrary work and actual inverse words. Independent
proof review covers the arbitrary-width statements and the scope of
each obstruction. Finite fixtures support the analytic proofs; they do
not establish the uniform degree bounds by enumeration.
