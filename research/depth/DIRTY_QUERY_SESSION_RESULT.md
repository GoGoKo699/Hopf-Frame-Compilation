# Dirty-query reformulation: execution outcome

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer)

This round began from verified main commit
`6134f3296f3c8c731f20f662aad453b02288be93` and its R96/R97 checkpoint.
The repository had not advanced beyond that commit. Work took place on
18 October 2026 (Asia/Shanghai), on the separate branch
`research/dirty-scalar-query-round`; earlier workspaces were preserved.
Final manuscript drafting and the selected Results A–D were not changed.

## 1. Revision and choice of question

The count package is settled under its stated allocations. R96 proves
the all-precision upper law

```math
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad N=2^n,
```

with one clean qubit, optimal T-count, elementary Clifford count
$`G=O(NL)`$, the original literal dirty thresholds, and the full
initialized-isometry contract. R79 remains selected Result D.
R97 proves a four-layer floor at bounded outer-axis misalignment; it
does not prove unrestricted one-clean rigidity. Neither result supplies
a dimension-growing large-width depth lower bound.

The formerly open general-magnitude recovery, recursive terminal
overhead and precision-window restrictions have already been resolved.
The present fixed-accuracy, square-root-width gap is a constant lower
bound versus logarithmic upper depth. The spectator-padded R97 witness
cannot close that dimension gap: its only active qubit has a
fixed-accuracy constant-depth compiler.

The chosen question was whether an arbitrary Boolean table can become
a literal scalar phase on every dirty-work input in sublogarithmic
address depth, with the existing square-root T-count and linear
Clifford budget. This tests a common cost in the upper construction and
gives a dimension-sensitive subfamily for the lower-bound direction.
The previous Gauss-sum encoding was revisited because it already has
constant depth and the desired counts, but exposes its table sign only
through a small matrix element or trace.

## 2. Exact query and frame reductions

[Boolean query/frame reduction](BOOLEAN_QUERY_FRAME_REDUCTION.md) is the
proof home for a two-call equivalence between arbitrary Boolean XOR
queries and the deepest corner-frame subfamily. It preserves the clean
and dirty allocations when at least one dirty wire is available, doubles
T-count and depth at most, and adds two elementary Clifford gates.
The actual inverse is essential in the frame-to-query direction.
Approximation error increases by at most a factor of two, including
leakage and arbitrary external references. The precision increment is
charged at the original literal width threshold.

This reduction quantifies over arbitrary native representations. It
therefore transfers a robust dirty-query lower bound to the full frame
family, if such a lower bound is proved. It is not itself a growing
lower bound, and its two-call converse concerns the deepest subfamily.

## 3. A conditional route for the general construction

[Conditional R98](QUERY_DEPTH_FRAME_TRANSFER.md) gives a complete
fixed-accuracy construction under an exact dirty word-query hypothesis.
If the R91 interface has additive address depth d(a), with the same
T-count, elementary Clifford count, minimum dirty pool and native
read-only address controls, then both complete frame families admit

```math
D_T=O_L\!\left(N/b^2+d(2n+10)+\log_2\log_2(n+4)\right)
```

at the original one-clean dirty thresholds, retaining optimal T-count,
$`G=O_L(N)`$ and full initialized-isometry error $`2^{-L}`$.
The theorem is conditional; the known query has logarithmic d.

The residual cost audit matters. Substituting a faster query alone
leaves the old logarithmic barrel router. The new
[paired permutation wrapper](PAIRED_PERMUTATION_ROUTING.md) moves a
selected permutation onto the payload and its inverse onto a separate
dirty array. The complete body and actual inverse layout return that
array exactly. Its short-label phase-grid construction trades polynomial
dirty volume for doubly logarithmic routing depth. Four global early
stages then leave a fully enumerated terminal whose native factors,
selection queries and two-array commutator fit the displayed residual.
No source, anchor, original prefix, precision cost or simultaneous pool
is removed from the ledger.

The two transfers have different hypotheses. The lower-bound transfer
in Section 2 accepts arbitrary native one-bit queries, including a clean
flag. R98 requires exact, zero-clean, multi-output queries with read-only
native address controls so that private banks can run in parallel.
This round does not prove those two compiler interfaces equivalent.

## 4. The Gaussian candidate and its stopping point

[Gaussian inverse extraction](GAUSSIAN_INVERSE_EXTRACTION.md) replaces
the trace-extraction idea with an exact Pauli commutator. A paid
in-place map for the affine matrix A and its actual inverse expose
$`(-1)^{c^{\mathsf T}A^{-1}b}`$ as a scalar while returning all work.
Four map/inverse calls are charged. The cheap doubled-register shear
instead exposes an entry of A; substituting it would omit the inverse
computation.

The same proof home gives a robust limitation of that repair. A product
of r address-preserving affine symplectic work rounds has symplectic
entries of Boolean degree at most r. An explicit table gives an inverse
entry of degree equal to the number of address bits. Thus that in-place
conversion cannot be obtained from constantly many such rounds, even
with returned initialized helpers. The exact serial construction
matches the round bound, but its native depth allowance is worse than
the retained R91 query.

This stops the affine-round inverse repair in its present form. The
bound is on that work-conversion architecture. It does not lower-bound
arbitrary native T-depth, and the scalar phase itself is invisible to
the symplectic matrix. A direct scalar loop or a representation outside
the affine, fiberwise-Clifford interface remains eligible.

## 5. Resumption question

Can an arbitrary Boolean word query be implemented exactly with
returned dirty work, the R91 count envelope and read-only native address
controls, but with sublogarithmic address T-depth? On the lower-bound
side, can full dirty return constrain
$`Q_f^\dagger Z_tQ_f=(-1)^{f(x)}Z_t`$ in a way that is stable under
fixed operator error and every alternative native representation?

The Gaussian work-conversion result does not answer either question.
No unrestricted depth lower bound or improved unconditional full-frame
upper law is claimed by this round. R96 and the full resource and
operator contract remain the baseline and destination.

## 6. Evidence and integration checks

Twelve new bounded tests cover the two-call query/frame identities,
literal phases, the actual inverse with initialized leakage, the
Gaussian edge maps and their inverse-degree witness, and the paired
router's complete arbitrary-input action and selector-preservation
requirement. These tests support the stated interfaces; the proofs,
conditional assembly and asymptotic resource bounds are analytic.

The reviewer walkthrough, bounded native real and complex examples,
q=16 residual example, four retained fault-tolerant receipt suites,
both n=12 resource ledgers, offline synchronization check and Python
compilation passed. No retained scientific receipt was regenerated.
The ten changed Markdown pages and five existing diagrams passed the
desktop and narrow rendering checks, including protected inline math
and native MathML. Independent technical reviews checked the two
different query hypotheses, simultaneous work reservations, actual
inverse return, four-stage terminal asymptotics and scoped obstruction.
