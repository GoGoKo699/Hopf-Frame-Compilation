# Exact and Fault-Tolerant Compilation of Hopf Differential Frames

Exact state preparation normally specifies one initialized input:

```math
|0^n\rangle|0^m\rangle
\longmapsto
|\psi\rangle|0^m\rangle.
```

The original inverse-frame Hopf decoder uses a prescribed unitary completion: the state in its
first column and coordinate-frame directions in designated other columns.
Its inverse resolves a common objective response.

| Synthesis task | Required action |
|---|---|
| Exact state preparation | fix $U\lvert 0^n\rangle$ |
| Hopf differential-frame compilation | fix $W\lvert x\rangle$ for every system basis state $\lvert x\rangle$ and restore workspace according to its input promise |

This repository asks whether the prescribed Hopf completion retains the
state-preparation size–depth frontier, and what finite precision costs.
The [publication scope](manuscript/PUBLICATION_SCOPE.md) fixes the full-frame claims.

<p align="center">
  <img src="assets/state-vs-frame.svg" width="900" alt="State preparation fixes one column, whereas Hopf differential-frame compilation fixes the state and designated frame columns." />
</p>

The [LLM reading guide](llms.txt) maps proofs and assumptions for complete-frame
synthesis, ancilla–depth and T-count–precision tradeoffs, and Hopf QBP robustness.
Resume research from the [workspace checkpoint](WORKSPACE.md).

## Two resource models

Both models target the complete frame. Their bounds use different
constructions, without simultaneous optimality of every cost.

| Model | Resource question | Scope of the matching theorem |
|---|---|---|
| Arbitrary one-qubit gates and CNOTs | exact size, CNOT count, and depth versus clean workspace | real and phase-dressed complex magnitude frames, every clean budget |
| Clifford+T | T-count versus accuracy and clean/dirty workspace | real and phase-dressed complex magnitude frames, every precision, under the sufficient clean reservation below |

### Exact logical theorem

Let

```math
N=2^n
```

and let $m\geq0$ be the number of clean ancillary qubits.  In the exact all-to-all
logical model with arbitrary one-qubit gates and CNOTs,

```math
\boxed{
S_{\mathbb R}(n,m)
=S_{\mathbb C,\mathrm{mag}}(n,m)
=\Theta(N)
}
```

and

```math
\boxed{
D_{\mathbb R}(n,m)
=D_{\mathbb C,\mathrm{mag}}(n,m)
=\Theta\left(n+\frac{N}{n+m}\right).
}
```

The upper bounds hold for every parameter tuple. The matching lower bounds
hold in the worst case over the Hopf-frame family, uniformly in the
clean-workspace budget. They do not impose the same cost on every individual
frame.

The worst-case **CNOT count alone is $\Theta(N)$ for $n\ge2$**,
with free one-qubit gates and arbitrary clean workspace; it is zero for $n=1$.
See the [separate lower-bound proof](docs/COMPILER_THEOREM.md#cnot-count-with-unrestricted-clean-workspace).

The real result concerns the complete Hopf differential frame.  The complex
result concerns the phase-dressed magnitude frame

```math
W_{\mathbb C,\mathrm{mag}}
=D_{\mathrm{ph}}W_{\mathbb R}.
```

The leaf-phase derivatives form a separate direct measurement stream; they are
not additional columns of the same $N$-dimensional unitary.

The construction uses at most the requested $m$ clean ancillary qubits and
returns them exactly to zero.  Thus a prescribed Hopf completion reaches the
optimal arbitrary-state-preparation frontier for every clean-workspace budget.

### Fault-tolerant theorem

Let $a$ be the number of initialized clean qubits and $b$ the number of borrowed
qubits in an arbitrary unknown state. Write

```math
q=n+a+b,\qquad
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad
h=1+\lceil\log_2(L+n+2)\rceil.
```

For $0\lt \eta\leq1/64$ and a sufficiently large fixed constant $C$, both the
real frame and phase-dressed complex magnitude frame have matching worst-case T-count

```math
\boxed{
T^\star(n,a,b,\eta)
=\Theta\left(\sqrt{NL}+L+\frac{NL}{n+a+b}\right),
\qquad a\geq C(n+h).
}
```

The upper construction uses $O(NL)$ Clifford gates. Every supplied precision
source is prepared and charged. Borrowed qubits are restored jointly with any
external reference; complete-input approximation includes clean-work leakage.
There are no measurements, resets, or free supplied catalysts in this model.

**Shared precision.** The [proof](docs/FAULT_TOLERANT_COMPILER.md) corrects a
coarse frame with one retained source and final amplification; all work is charged.

At $L=N$, sufficient $a=\Theta(n)$ and $b=\Theta(N)$ give $T^\star=\Theta(N)$.
With **one clean qubit** and $`b\ge L+n+7`$, the
[one-clean compiler](docs/ONE_CLEAN_COMPILER.md) gives
$`T=O(N+L\ell_*(n))`$, $`G=O(NL)`$ for real frames, where
$`\ell_*(n)=1+\log_2^*(n+2)`$. The iterated logarithm counts base-two logs
to at most one. At $`L=N`$, the
[open endpoint](docs/OPEN_PROBLEM.md) is $`\Omega(N)\le T^\star\le O(N\ell_*(n))`$.
Smaller allocations retain the
[borrowed-workspace bound](docs/BORROWED_WORKSPACE_COMPILER.md).
At $`b\ge L+n+7`$, [zero-clean real-frame synthesis](docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
gives $`T=O(N+nL)`$, $`G=O(NL)`$; the grouped improvement still uses one clean qubit.

Composing with a literal diagonal of separately supplied certified phases extends the
[one-clean bound](docs/ONE_CLEAN_COMPILER.md#8-phase-dressed-complex-magnitude-frames) to phase-dressed complex
magnitude frames at $`b\ge L+n+8`$. With $`b\ge2(L+n+8)`$, it gives
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$ T gates, still $`G=O(NL)`$.
Leaf-phase derivatives retain a separate QBP stream.

At fixed accuracy, [parallel lookup](docs/PARALLEL_DIRTY_LOOKUP.md)
gives $`T=O(\sqrt N)`$ and $`D_T=O(n^3)`$ in one real-frame circuit
with two clean and sufficiently large $`\Theta(\sqrt N)`$ dirty workspace.

Beyond Hopf frames, **literal diagonals and general one-target U(2)
multiplexors** attain $`\Theta(\sqrt{NL}+L+NL/b)`$ with one clean
qubit, respectively at $`b\ge2(L+n+7)`$ and $`b\ge2(L+n+9)`$.
For the multiplexor, $n$ counts address qubits and $N$ counts its blocks.
These matched one-stage frontiers do not close the multi-layer frame gap.

Two promised update classes attain $`T=O(N+L)`$, $`G=O(NL)`$ at
$`b\ge L+n+7`$. Both require literal agreement with a supplied baseline
outside the changed nodes; its native determinant-one local words have
length $`O(n-d+1)`$ at depth d.

- [Antichain changes](docs/ANTICHAIN_COMPILER.md) use **zero clean qubits**;
  no changed node is an ancestor of another.
- [Sparse nested changes](docs/SPARSE_UPDATE_COMPILER.md) use **one clean qubit**.
  If S is the ancestor closure of the changed nodes, a sufficient condition is

```math
n\ge 2\lceil\log_2(|S|+1)\rceil+32.
```

The sparse result includes arbitrary changes along one root-to-leaf path,
for all n using a finite small-n fallback. Their real Hopf specializations preserve the complete frame and QBP interface. Generic rounding
supplies neither promise; the unrestricted linear endpoint remains open.

## Exact construction at a glance

The Hopf frame is a product of addressed tree layers.  At depth $d$, the prefix
selects one rotation angle and the complete lower suffix supplies a shared
all-zero predicate.  Three schedules exploit this structure.

| Workspace | Schedule | Mechanism |
|---:|---|---|
| $m=0$ | borrowed-suffix echo | one original suffix data qubit carries the predicate temporarily and is restored exactly |
| $1\leq m\lt 4n$ | direct flagged UCG | one reusable clean flag stores the suffix-zero predicate |
| larger $m$ | routed parallel subframes | a tree cut turns the tail into a direct sum; the suffix is routed coherently and the subtree frames run in parallel |

<p align="center">
  <img src="assets/literature-lineage.svg" width="940" alt="The all-workspace state-preparation line and the Hopf differential-frame line meet in the optimal complete-frame compiler." />
</p>

## Read the repository in three passes

| Time | Route | Purpose |
|---:|---|---|
| 5 minutes | this page | problem, theorem, and construction map |
| 30–40 minutes | **[complete technical narrative](REVIEW.md)** | shared contract, the two resource models, and the QBP consequence |
| full audit | **[exact theorem](docs/COMPILER_THEOREM.md)**, **[T-count theorem](docs/FAULT_TOLERANT_COMPILER.md)**, **[one-clean theorem](docs/ONE_CLEAN_COMPILER.md)**, and **[verification map](docs/VERIFICATION.md)** | proofs, register schedules, evidence limits, and source dependencies |

The [documentation map](docs/README.md) includes the
[QBP consequence](docs/QBP_CONSEQUENCE.md), [related work](docs/RELATED_WORK.md),
and [Hopf geometry](docs/HOPF_INTERFACE.md).

## Why the completion matters

The complete frame satisfies

```math
W_{\mathbb R}|0^n\rangle=|\psi\rangle,
\qquad
W_{\mathbb R}|\lambda(j)\rangle=|e_j\rangle.
```

At a regular coordinate,

```math
\partial_{\theta_j}|\psi\rangle
=a_j|e_j\rangle,
\qquad
g_{j,j}=a_j^2,
```

where $a_j$ is the oriented amplitude entering the corresponding tree node.  On
the canonical Hopf domains, $a_j\geq0$ and equals the principal metric square
root.  At zero metric weight the raw derivative vanishes, while the marker
column remains the chart-selected orthogonal continuation determined by the
complete parameter tuple.

A state-preparation-equivalent completion may move these marker columns.  The
repository contains an exact two-qubit example in which the state is unchanged
but the decoded gradient changes from

```math
(2,0,0)
\quad\longmapsto\quad
(0,\sqrt2,0).
```

For this fixed decoder, the relevant compiler contract is

```math
\widetilde W
\bigl(|\varphi\rangle|0^m\rangle\bigr)
=(W|\varphi\rangle)|0^m\rangle
```

for every system input $\lvert\varphi\rangle$, not only the forward preparation input.
At regular points, universal preservation of the fixed decoder's means
also forces this contract up to common phase; see the
[necessity theorem](docs/FRAME_SAFE_COMPILATION.md#necessity-for-all-observable-dependent-gradient-means).

## Relation to the all-workspace state-preparation framework

The proof uses four exact results from the all-workspace state-preparation
framework.

| Imported result | Role here |
|---|---|
| optimal QSP frontier | benchmark and matching comparison |
| ancilla-free multi-controlled X | zero-suffix predicates and toggles |
| all-workspace UCG synthesis | prefix-selected rotations, subtree frames, and the phase diagonal |
| coherent CNOT copy–uncopy | control fanout for the decoder and router |

The Hopf-specific contributions are the complete-operator factorization and
three workspace schedules, which adapt these primitives while preserving
marker columns.

The normative compiler citation is:

> P. Yuan and S. Zhang, “Optimal (controlled) quantum state preparation and
> improved unitary synthesis by quantum circuits with any number of ancillary
> qubits,” *Quantum* **7**, 956 (2023).

## Consequence for quantum backpropagation

Frame-safe compilation preserves the complete inverse-frame measurement
distribution.  Under the primary finite-shot target of simultaneous absolute
accuracy for the raw Hopf-coordinate gradient, the global magnitude stream uses

```math
O(\log n)=O(\log\log M)
```

independent executions at fixed accuracy and confidence, where
$M=\Theta(2^n)$ is the number of coordinates.

This is a matched-program statement: scalar and gradient programs use the same
forward preparation family and controlled observable, while the gradient
program adds one inverse frame of the same asymptotic logical depth as optimal
state preparation.  Classical materialization of an $M$-entry output is not
included in that quantum-depth ratio.

At finite precision, the [approximation bridge](docs/QBP_APPROXIMATION.md)
controls the bias of the same fixed-parameter, bounded-score estimator. It uses
the actual compiled circuit and its actual adjoint. It does not differentiate a
discontinuous family of compiled words. Frame-synthesis lower bounds alone do
not establish optimal gradient-query or training complexity.
The same proof now covers the complete complex gradient, reflection-sum
observables, rounded classical weights, and correlated dirty-bank reuse
between executions with fresh declared clean inputs.

The separate [state-based QBP theorem](docs/STATE_BASED_QBP_THEOREM.md)
changes the decoder and avoids compiling the fine complete frame.
For real and phase-dressed complex Hopf states, observable coefficient norm $`\Lambda`$,
and error $`\epsilon`$, put
$`K=\max\{6,\lceil\log_2(80\Lambda/\epsilon)\rceil\}`$ and $`P=\max\{n,K\}`$.
Two initialized compiler flags and $`b\ge P+n+7`$ dirty wires give
$`T=O(N+P)`$, $`G=O(NP)`$ per execution apart from the observable;
the protocol branch is separate. The actual coarse inverse and classical
correction retain the shot scaling. [Depth schedules](docs/STATE_QBP_DEPTH.md)
and bounded-input costs are explicit.
[Bounded residual QBP](docs/NATIVE_RESIDUAL_QBP.md) implements both streams.
[General software](docs/VERIFICATION.md#state-based-qbp-coverage) remains
optional; the complete-frame endpoint is open.
No end-to-end speedup is claimed.

## Verification boundary

Checks cover complete logical identities, explicit decoder/router gates,
arbitrary-entangled inputs, exact resource ledgers, and gradient decoders.

UCG and multi-controlled-X decompositions are imported. Finite checks expose
indexing, phase, order, cleanup, and resource errors; the asymptotic claims
rest on proofs.

The fault-tolerant evidence adds exact finite source/kernel identities and
rational resource checks. It does not yet provide a general elementary emitter
for the complete asymptotic shared-source compiler. See the
[focused reproduction guide](verification/fault_tolerant/README.md).

Device connectivity, physical noise thresholds, optimal T-depth, arbitrary
non-Hopf charts, and application-independent observable costs remain outside
the claims. Clifford work, T work, quantum executions, and classical output
costs are reported separately.

## Reproduce the checks

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/reviewer_walkthrough.py
python validate.py
python scripts/verify_fault_tolerant.py
python scripts/unified_resource_ledger.py --n 12
python scripts/strict_zero_echo_ledger.py --n 12
```

## Relationship to the Hopf repositories

| Repository | Role |
|---|---|
| [`Hopf-ansatz`](https://github.com/GoGoKo699/Hopf-ansatz) | coordinate chart, inverse map, metric, tangent preparation, and optimization interface |
| [`Hopf-QBP`](https://github.com/GoGoKo699/Hopf-QBP) | global, direct-phase, and checkpoint gradient records, including the earlier Möttönen-style robustness result |
| **This repository** | complete-frame contracts, exact and fault-tolerant compilation, and their scoped QBP consequences |

## Status and license

[Evidence](docs/VERIFICATION.md) and [sources](docs/SOURCE_MAP.md) delimit the claims.
The repository uses the [MIT license](LICENSE).
