# Returned-query invariants: execution outcome

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer)

This continuation began from verified main
`ac3038452451bf7a08e89ad5b6af8fb277a12eae`, the merge of PR #165.
It used the separate branch `research/returned-query-invariant-round`
and preserved earlier workspaces. Selected manuscript claims are unchanged.

## 1. Revision and exact cost addressed

The count package is settled under its stated allocations. R96 remains
the strongest unconditional general research upper bound:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

This holds for every $`n\ge1,L\ge6`$, one clean qubit, real dirty
width $`b\ge17(L+n+7)`$ or complex magnitude width
$`b\ge17(L+n+8)`$, and full initialized-isometry error
$`2^{-L}`$. Literal phases, all logical columns, singular and zero
coefficients, arbitrary dirty inputs, external references, actual
inverses, leakage and complete return remain part of the contract.
R79 remains selected Result D. R97 excludes bounded outer-axis
misalignment in four-layer one-clean circuits; unrestricted rigidity
and dimension-growing native lower bounds remain open.

General-magnitude recovery, recursive terminal overhead and the old
precision-window restrictions have already been resolved. At fixed
accuracy and square-root width, the general gap remains constant
lower depth versus logarithmic upper depth. Conditional R98 would
improve the complete circuit if R91's exact zero-clean, multi-output
query admitted sublogarithmic address depth with its original count,
Clifford, width and read-only interface.

The paid cost targeted here was R91's logarithmic address-indicator
depth. Its central bilinear echo already has a constant eight T layers.
The previous Gaussian round stopped at an in-place inverse map: its
symplectic invariant could not see a scalar sign. This round therefore
tested direct closed scalar loops, more finely split indicator banks,
and native query invariants with initialized work. A faster encoding or
an uncharged change of basis would not improve the baseline.

## 2. A scalar loop exists, but its Fourier map carries the missing cost

The [Gaussian scalar-loop proof](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
is the authoritative home for this round's main result. For a quadratic
form q with nonsingular polar B and zero constant term, its generalized
Fourier transform satisfies the literal identity

```math
(F_BD_q)^3=(-1)^{\mathrm{Arf}(q)}I.
```

This returns every work input. But $`F_B=H^{\otimes m}C_B`$ requires
the in-place map $`C_B|v\rangle=|Bv\rangle`$. Three map calls,
three phase calls, all Clifford gates, simultaneous helper allocation
and six-factor approximation error are charged in the proof home.
A cheap doubled shear does not supply that map. The identity is an
exact reformulation, not a faster query compiler.

The stronger obstruction addresses closed loops themselves. In a
circuit of r actual affine-address quadratic phase rounds, interleaved
with fixed work Cliffords and address-preserving Clifford controls,
literal zero-column return implies

```math
\deg_{\mathbb F_2} f\le4r+2.
```

The same conclusion holds at error below $`\sqrt{2-\sqrt2}`$,
with any number of initialized or dirty helpers. Full operator return
implies the zero-column premise. Thus a degree-a table needs linearly
many such rounds. Arbitrary scalar lifts of affine symplectic matrices
are not eligible primitives.

The missing invariant is supplied by inserting diagnostic Pauli
markers in a quadratic path sum. Complete zero-column return makes its
polar matrix nonsingular. Mixed marker phases recover its inverse
polar, whose entries have degree at most $`2r`$. A basis-independent
formula for the Arf sign then gives the displayed bound. Individual
Pauli-lift signs cancel from the mixed marker coefficients. Complex
Cliffords and the approximation threshold are handled explicitly.

This closes the previously untested bounded-round Gaussian scalar
escape. It does not prove a native T-depth lower bound: arbitrary
resynthesis may change the round architecture, and a native T layer
can give non-Clifford work dynamics at a fixed address.

## 3. Native two-layer queries with returned initialized work

The independent [read-only rigidity proof](READ_ONLY_TWO_LAYER_RIGIDITY.md)
uses physical T layers. Every completed Clifford interlayer must
preserve the address pointwise, and all T gates act on work. At full
error below $`\sqrt2/6`$, a two-layer XOR query with c clean qubits
obeys:

| Initialized work | Necessary Boolean structure |
|---|---|
| None | f is affine; a separate proof permits the larger threshold $`\sqrt{2-\sqrt2}/2`$ |
| One clean qubit | Degree at most two and quadratic polar rank at most two |
| Any $`c\ge1`$ | Degree at most $`2c`$ |

Dirty width is arbitrary. The proof retains actual inverses and clean
leakage. Endpoint clean syndromes expose a large address subspace on
which the logical Pauli labels form a projective representation; their
coset slopes give the degree bound. These are necessary conditions,
not synthesis or tightness claims.

A retained exact dirty Toffoli has two T layers when transformed
address wires may be T targets. Its complete native fixture explicitly
falsifies enlarging the zero-clean statement to arbitrary circuits.
The new theorem constrains the R98-style interface; it cannot be
inserted into the unrestricted query-to-frame lower-bound reduction.

## 4. Constructive alternatives and why they were stopped

[Multiway query cost barrier](MULTIWAY_QUERY_COST_BARRIER.md) is the
single proof home for the alternatives audited here.

* Splitting the address into k indicator fields shortens each bank,
  but representing arbitrary tables by r separated products requires
  $`r\sum_i2^{a_i}\ge2^a`$ for a fixed split. The termwise
  square-root-count envelope therefore leaves an $`\Omega(a)`$-bit
  field. Table-dependent splits and coherent approximate Boolean
  branches retain the asymptotic obstruction. Joint native synthesis
  is outside this representation-count argument.
* Replacing one-hot generation by the new paired router is circular:
  that router invokes indicators, and a full address-size phase grid
  also exceeds the query budget.
* A rank-deficient linear-hash bucket toggle is orthogonal to the
  required one-hot output on zero work. Coherent seeding and actual
  unpreparation do not remove its full-query error.
* A cyclic borrowed-slot commutator leaves adjacent program ratios on
  every slot. Those residual operators cannot be dropped; enforcing
  their cancellation restores the original prefix-product work.

These diagnoses stop the specific substitutions. They leave joint
native interference and other representations available.

## 5. Precise unresolved step and resumption question

**Subsequent resolution:** the
[three-layer continuation](THREE_LAYER_QUERY_SESSION_RESULT.md) proves
degree at most four with one clean qubit at error below $`\sqrt2/6`$,
within precisely the interface asked about below. It also sharpens the
zero-clean conclusion to affineness. The following question records this
round's historical handoff; the later checkpoint gives the current one.

The most useful immediate question is whether **three native T layers**
with preserved address boundaries and one clean qubit can implement
returned Boolean queries of unbounded degree as dirty width grows.
Either an explicit family or a width-independent degree bound would
test the first case outside this round's native invariant. It targets
the R98 query cost while allowing non-Clifford work interference,
which the Gaussian theorem excludes by hypothesis.
An obstruction with one clean qubit also covers zero clean qubits.
A positive one-clean construction would need a separately paid
clean-removal or borrowing step before it could meet R98's zero-clean
query interface; the frame compiler's sole flag is already reserved.

The missing step is concrete: for two layers, three Pauli-control
bands reduce to one projective representation on a clean-syndrome
kernel. With more bands, the intervening dressed families need not
preserve the initialized code or act as logical Paulis. Iterating the
two-layer argument is not a proof. A proposed higher-depth invariant
must survive their coherent cancellations, actual return and fixed
error; a proposed construction must price every nonlinear control,
source, Clifford gate and simultaneously occupied helper.

A positive high-degree example alone would not compile arbitrary
tables. It would identify an escape mechanism to test against R91's
full multi-output count envelope. A degree bound in the read-only
interface would constrain that upper route; obtaining an unrestricted
frame lower bound additionally requires an invariant stable under
address mixing and every native representation. The
[two-call reduction](BOOLEAN_QUERY_FRAME_REDUCTION.md) remains available
once that stronger query lower bound exists.

No unconditional complete-frame upper bound or dimension-growing
native T-depth lower bound changed in this round. R98 remains
conditional, and optimal general depth remains the continuing objective.

## 6. Evidence

Exact finite fixtures check the Arf identity, the literal Fourier cube,
marked path sums, read-only Clifford phases, initialized inverse hybrids
with nonzero leakage, and the unrestricted two-layer Toffoli scope
counterexample. The width-independent and approximate theorems are
analytic. Independent reviews checked the phase primitive, marker
degrees, realification, initialized-sector argument and the distinction
between native depth and completed Gaussian rounds.

All twelve new bounded tests passed. The reviewer walkthrough, native
real and complex examples, q=16 residual example, retained
fault-tolerant receipt suites, both n=12 resource ledgers, offline sync
check and Python compilation passed. Ten changed Markdown pages and
five retained diagrams passed desktop and narrow rendering checks.
The full-suite result is recorded with the merge; no retained scientific
receipt or manuscript claim is changed.
