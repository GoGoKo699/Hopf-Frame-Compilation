# Exact and Fault-Tolerant Compilation of Hopf Differential Frames

State preparation fixes one column of a unitary. The inverse-frame Hopf
gradient decoder also uses designated coordinate-frame columns. This
repository compiles that prescribed **complete operator**, with explicit
costs for precision, initialized workspace, and borrowed workspace.

The principal results are A–D below. Their hypotheses and retained
corollaries are collected in the [publication scope](manuscript/PUBLICATION_SCOPE.md).

<p align="center">
  <img src="assets/state-vs-frame.svg" width="900" alt="State preparation fixes one column; complete-frame compilation also preserves the designated coordinate-frame columns." />
</p>

| Read next | Purpose |
|---|---|
| [Proof roadmap](REVIEW.md) | Follow the mechanisms connecting the four claims |
| [Documentation map](docs/README.md) | Find authoritative proofs and their dependencies |
| [Verification](docs/VERIFICATION.md) | Inspect executable evidence and its limits |
| [Research results and open gaps](research/README.md) | Find restricted successes and limits of specific compilation methods |

## Target and resource contract

Write $`N=2^n`$, $`n\ge1`$. The real frame has state column
$`W_{\mathbb R}|0^n\rangle=|\psi\rangle`$ and designated marker columns
$`W_{\mathbb R}|\lambda(j)\rangle=|e_j\rangle`$.
At regular coordinates, $`\partial_{\theta_j}|\psi\rangle=a_j|e_j\rangle`$.
The supplied parameter tuple also fixes the frame at singular charts.

The complex magnitude frame is $`D_{\rm ph}W_{\mathbb R}`$, with separately
supplied leaf phases. Phase derivatives use a separate measurement stream;
they are not extra columns of this unitary.

The exact model uses arbitrary one-qubit gates and CNOTs, all-to-all
connectivity, and $`m`$ clean ancillas returned exactly. The coherent
Clifford+T model uses $`a`$ clean and $`b`$ arbitrary dirty qubits. If
$`J_a`$ appends the clean zero state, its approximation contract is

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

This includes every system input, dirty work entangled with any reference,
literal phases, and clean-work leakage. Return error is included in this
same norm; Result B returns its dirty registers exactly. Source preparation,
queries, and actual inverses are charged. There are no measurements, resets,
or free supplied precision sources. Certified classical coefficient and
table generation are accounted for separately from quantum gate counts.

## Four principal results

The upper bounds hold for every parameter tuple under the stated hypotheses.
Matching lower bounds hold in the worst case over the Hopf-frame family.
For Result A they hold uniformly in the clean-workspace budget. Different
constructions attain the different resource guarantees.

**A. Optimal exact compilation, at every clean budget.** For real and
phase-dressed complex magnitude frames, every integer $`m\ge0`$ admits

```math
S=\Theta(N),\qquad D=\Theta\!\left(n+\frac{N}{n+m}\right).
```

The worst-case CNOT count alone is $`\Theta(N)`$ for $`n\ge2`$, even
with free one-qubit gates and unrestricted clean workspace; it is zero for
$`n=1`$. The [exact proof](docs/COMPILER_THEOREM.md) combines a borrowed-suffix
echo, a direct flagged schedule, and coherently routed parallel subframes.

For the approximate results, put

```math
0\lt\eta\le1/64,\quad L=\max\{6,\lceil\log_2(1/\eta)\rceil\},
\quad h=1+\lceil\log_2(L+n+2)\rceil,\quad q=n+a+b.
```

**B. Matching T-count with sufficient clean work.** For both frame families,
when $`a\ge C(n+h)`$ for a sufficiently large absolute constant $`C`$,

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{q}\right).
```

The [fault-tolerant proof](docs/FAULT_TOLERANT_COMPILER.md) shares one
charged precision source across the complete residual-correction stream.
The clean reservation is sufficient, not a clean-space lower bound.

**C. One-clean compilation.** For real frames, $`a=1`$ and
$`b\ge L+n+7`$ suffice for

```math
T=O(N+L\ell_*(n)),\qquad \ell_*(n)=1+\log_2^*(n+2).
```

Here the iterated logarithm counts base-two logs to at most one.
The [one-clean proof](docs/ONE_CLEAN_COMPILER.md) and
[conditional-suffix grouping](docs/CONDITIONAL_SUFFIX_COMPILER.md) also
give a banked tradeoff. The separately proved complex magnitude extension
requires $`b\ge L+n+8`$. Literal diagonals and one-target multiplexors
have their own matching one-clean corollaries and reservations.
For literal diagonals, the
[packed two-clean compiler](docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames)
uses only $`b\ge n+1+\lceil(L+4)/2\rceil`$ dirty qubits at
$`T=O(N+L)`$ and $`G=O(NL)`$.

**D. Simultaneous T-count and T-depth.** For complete real frames,
$`a=2`$ and $`b\ge17(L+n+7)`$ give one circuit with absolute constants:

```math
T=O(\sqrt{NL}+NL/b+nL),\qquad D_T=O(NL/b^2+nL).
```

For $`6\le L\le\log_2(n+2)/16`$, the same theorem strengthens this to

```math
T=O(\sqrt{NL}+NL/b),\qquad D_T=O(NL/b^2+n).
```

Both resources match simultaneously on the following intervals, when nonempty:

| Precision regime | Matching dirty-width interval |
|---|---|
| Every $`L\ge6`$ | $`17(L+n+7)\le b\le\sqrt{N/n}`$ |
| $`6\le L\le\log_2(n+2)/16`$ | $`17(L+n+7)\le b\le\sqrt{NL/n}`$ |

There, $`T^\star=\Theta(NL/b)`$ and $`D_T^\star=\Theta(NL/b^2)`$.
The [uniform-precision proof](docs/UNIFORM_PRECISION_DEPTH.md) charges
source preparation, queries, and return on this same circuit. T-depth allows
arbitrary Clifford interlayers and is distinct from total depth. No automatic
complex extension or unrestricted depth optimality is claimed.

All upper constructions in B–D use $`G=O(NL)`$ elementary Clifford gates.
The [full statements](manuscript/PUBLICATION_SCOPE.md) retain the additional
regimes, complex extensions, and older schedules useful at smaller reservations.

## Two open resource gaps

| Question | Current bounds |
|---|---|
| High precision, $`\eta=2^{-N},a=2,b=N+n+7,n\ge3`$ | $`2N-1\le T^\star\le O(N\ell_*(n))`$ |
| Fixed accuracy, $`a=2`$, sufficient $`b=\Theta_\eta(\sqrt N)`$ | Optimal count $`T=\Theta_\eta(\sqrt N)`$ is attained with $`D_T=O_\eta(n)`$; unrestricted depth lower bound remains $`\Omega(1)`$ |

An [explicit single-angle witness](docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
gives $`T\ge2N-1`$ at the count endpoint, even with unrestricted
Clifford interlayers and workspace, under the literal-phase contract.
Two constructive reductions isolate the count question:
[any fixed number of diagonal factors](research/endpoint/BOUNDED_DIAGONAL_FACTORIZATION.md)
fits the exact workspace with certified factor selection, and
[every frame reduces to a regular fixed-tree Cayley family](research/endpoint/TREE_CAYLEY_REDUCTION.md)
using two diagonals and exact O(N)-T permutations. Its
[propagation and feedback proof](research/endpoint/BOUNDARY_PROPAGATION.md)
gives a precision-feasible repeated-call block and an exact cancellation
leaving one complete prefix encoder. The
[joint coarse encoder](research/endpoint/COARSE_PREFIX_ENCODER.md) and
[collective refinement](research/endpoint/COLLECTIVE_PRECISION_REFINEMENT.md)
give a native cubic replacement for the regular-core terminal frame:
error $`40\,2^{-3s}+2^{-L}`$ at O(N) T-count, where
$`s=\lceil N/n\rceil`$, $`n\ge16`$, and $`3s\le L\le N`$.
This fixed-order correction includes dirty/reference return and leaves
the endpoint bounds unchanged. The
[structural restrictions](research/endpoint/STRUCTURAL_COMPILATION_LIMITS.md)
separate excluded Haar and dense-mixer families from unrestricted synthesis.
The [research index](research/README.md) collects these reductions and
method-specific limits with their precise hypotheses.

## Operational consequence

A state-equivalent completion can change the decoded gradient. An
[exact two-qubit example](docs/COMPILER_BOUNDARIES.md#2-two-qubit-global-state-column-counterexample)
changes $`(2,0,0)`$ to $`(0,\sqrt2,0)`$ without changing the prepared state.
Among completions preparing that same state exactly, universal preservation
of the fixed decoder's means at regular points forces the full frame up to
common phase.

Exact substitution preserves its complete measurement record. With
phase-calibrated controlled access to the specified Hermitian-unitary
observable, simultaneous absolute accuracy of the raw coordinate gradient
uses $`O(\log n)`$ independent executions at fixed accuracy and confidence.
The [QBP consequence](docs/QBP_CONSEQUENCE.md) and
[approximation bridge](docs/QBP_APPROXIMATION.md) separately charge executions,
circuit resources, and classical output; they use the actual compiled
circuit and its actual adjoint, not derivatives of compiled gate words.

The completed [state-based QBP branch](supplements/state_based_qbp/README.md)
changes the decoder and has a separate task-level guarantee. Full-frame
necessity here concerns the fixed decoder. No optimality among all gradient
algorithms or end-to-end training speedup is claimed.

## Evidence and reproduction

This is a theoretical construction with executable finite checks. The
[verification map](docs/VERIFICATION.md) distinguishes complete logical
identities, explicit decoder/router gates, native finite primitives, exact
resource ledgers, and analytic asymptotic arguments. Yuan–Zhang's exact UCG and multi-controlled-X synthesis supply the
imported toolkit; the [related-work comparison](docs/RELATED_WORK.md)
separates those premises from the prescribed-frame arguments. There is no
general elementary emitter for the complete asymptotic shared-source compiler.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reviewer_walkthrough.py
python validate.py
python scripts/verify_fault_tolerant.py
```

The [focused reproduction guide](verification/fault_tolerant/README.md)
and [source map](docs/SOURCE_MAP.md) identify evidence and attribution.
Hardware connectivity and physical noise thresholds are outside the model.

Hopf geometry and gradient records originate in
[Hopf-ansatz](https://github.com/GoGoKo699/Hopf-ansatz) and
[Hopf-QBP](https://github.com/GoGoKo699/Hopf-QBP); this repository supplies
the complete-frame compiler claims and their scoped consequences.
Contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
Citation metadata is in [CITATION.cff](CITATION.cff); code is under the
[MIT license](LICENSE).
