# Where the selected compiler claims stop

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) · [Research archive](../research/README.md)

## Current decision

Results A–D and their retained corollaries pass the
[bounded internal claim-to-proof audit](CORE_CLAIM_AUDIT.md). The selected
scientific scope is frozen. This is an internal readiness decision, not
external peer review or formal proof certification. Manuscript writing
remains on hold.

Two resource questions remain open. They are distinct from the proved
claims and are not prerequisites for the selected paper. The
[research archive](../research/README.md#what-we-tried) records what was
tried, the useful results obtained, and why each route stopped.

Write $`N=2^n`$, $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
and $`\ell_*(n)=1+\log_2^*(n+2)`$. The coherent Clifford+T model
has a clean flags and b arbitrary dirty qubits, no QRAM, no supplied
precision source, and no intermediate measurement or reset. The full
initialized-isometry contract includes dirty/reference return and leakage.

## Constant-clean high-precision count

At the literal endpoint

```math
a=2,\qquad L=N,\qquad b=N+n+7,\qquad n\ge3,
```

the current worst-case real-frame bounds are

```math
\Omega(N)\le T^\star\le O(N\ell_*(n)).
```

The upper bound follows from the [one-clean grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit),
which may leave the second available flag unused. The lower bound is
inherited from [the full-frame counting argument](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage).
The sufficient-clean optimum in Result B uses a growing clean reservation
and therefore does not close this endpoint. Result D requires
$`b\ge17(L+n+7)`$ and does not apply at this literal width either.

A successful T-only resolution could have a larger fully charged Clifford
cost. What is missing is a complete native construction with one global
precision cost, or a stronger unrestricted lower bound. There is no
selected qualifying construction and no basis for claiming that closure
is close.

The [endpoint studies](../research/README.md#constant-clean-high-precision-count)
retain source reuse, tree transport, weighted blocks, residual assembly,
and canonical scalar attempts. Antichain and sparse updates are proved
restricted successes. They do not supply the promise for generic frames.

## Large-width T-depth

At fixed $`0\lt\eta\le1/64`$, two clean flags, and sufficient
$`b=\Theta_\eta(\sqrt N)`$, the proved construction has

```math
T=\Theta_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(n).
```

The unrestricted depth lower bound at that width is only $`\Omega(1)`$.
The explicit matching intervals at smaller widths remain those of
[Result D](UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary).
T-depth allows arbitrary Clifford interlayers; total physical depth is a
separate resource.

The smallest substantial new construction target is $`D_T=o(n)`$ with
those same clean, dirty, T-count, and Clifford budgets. The
[depth studies](../research/README.md#large-workspace-t-depth) reduce source,
predicate, and late-query overheads, but logical transport and program
query/unload each retain an O(n) allowance. These are costs of existing
constructions, not lower bounds on all possible circuits.

The exact retained-source factors and shallow stabilizer completion do
not remove the ordered transport. A source rewrite alone does not settle
the full-frame question.

## What would justify reopening

Before another construction pass, supply an explicit native identity or
encoding rule with symbolic count, depth, width, error, and cleanup bounds.
Then use small unequal examples to test that rule. A compact matrix or
classical recurrence alone is insufficient.

- For the endpoint, fit $`b=N+n+7`$ literally. Charge derived tables,
  basis changes, actual queries/unloads, and every precision-dependent call.
  A grouped route needs a precision coefficient independent of its group
  count; another full-frame route may bypass grouping.
- For depth, reduce both remaining linear allowances or give a complete
  circuit that bypasses them. A component proposal must state the other
  allowance that remains.
- An exact source rewrite must hold on every source input. A construction
  agreeing only on an ideal source needs its own full initialized-isometry
  proof, including inactive inputs and rejected work.
- A lower-bound proposal must survive the permitted Clifford interlayers,
  initialized flags, and returned dirty work. Failure of one query or source
  interface is not such a lower bound.

Stop a candidate whose expanded schedule retains the same asymptotic cost,
uses unpriced initialized history, or lacks the return proof. Do not repeat
larger instances of the same failed rule. A discovered defect in an existing
claim reopens that repair directly.

## Separate completed work

The [state-based QBP supplement](../supplements/state_based_qbp/README.md)
estimates the original raw gradients with a changed decoder. Its state
compiler, charged coarse inverse, classical correction, and bounded native
integration are complete at their stated scope. It does not implement the
fine prescribed frame, so it neither closes nor invalidates either gap.
Application-level advantage is outside the selected project scope.
