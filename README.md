# Exact and Fault-Tolerant Compilation of Hopf Differential Frames

State preparation fixes one column of a unitary. The inverse-frame Hopf
gradient decoder also uses designated coordinate-frame columns. This
repository compiles that prescribed **complete operator**, with explicit
costs for precision, initialized workspace, and borrowed workspace.

The principal results are A–D below. They include optimal exact compilation,
matching T-count frontiers with sufficient clean work or one clean qubit,
and simultaneous count/depth guarantees. Their hypotheses and retained
corollaries are collected in the [publication scope](manuscript/PUBLICATION_SCOPE.md).

<p align="center">
  <img src="assets/state-vs-frame.svg" width="900" alt="State preparation fixes one column; complete-frame compilation also preserves the designated coordinate-frame columns." />
</p>

| Read next | Purpose |
|---|---|
| [Proof roadmap](REVIEW.md) | Follow the four claims and their proof mechanisms |
| [Documentation map](docs/README.md) | Find authoritative proofs and their dependencies |
| [Verification](docs/VERIFICATION.md) | Inspect executable evidence and its limits |
| [Current depth frontier](docs/OPEN_PROBLEM.md) | Locate the proved resource bounds and the remaining one-clean question |
| [Research catalogue](research/README.md) | Find the endpoint proof dependencies, earlier constructions, and scoped limits |
| [Working checkpoint](WORKSPACE.md) | Resume the active research question and its full operator/resource contract |
| [LLM guide](llms.txt) | Follow the source-of-truth, scope and evidence rules |

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
postselection, QRAM or free supplied precision sources. Certified classical
coefficient and table generation are separate from quantum gate counts.

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
$`n=1`$. The [exact proof](docs/COMPILER_THEOREM.md) covers every clean budget.

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

**C. Optimal one-clean count.** For every $`n\ge1`$ and $`L\ge6`$,
one clean qubit suffices. Define the base dirty reservation B by

| Family | Base dirty reservation B |
|---|---|
| Complete real frame | $`L+n+7`$ |
| Phase-dressed complex magnitude frame | $`L+n+8`$ |

The [one-clean theorem](research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives $`T=O(N+L)`$ and $`G=O(NL)`$ at every $`b\ge B`$.
At every $`b\ge2B`$, it attains the full worst-case count frontier:

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right).
```

At exactly $`b=B`$ and $`L\ge\max\{6,n\}`$, the base-width count
also matches: $`T^\star=\Theta(N+L)`$. No general matching claim is made
at base width with $`L\lt n`$ or throughout $`B\lt b\lt2B`$.
Result C adds no T-depth bound.

**D. Optimal T-count with all-precision T-depth bounds.** The selected
[uniform subtree-depth theorem, R79](research/depth/UNIFORM_SUBTREE_DEPTH.md),
uses one clean qubit at $`b\ge17B`$, with the family-specific B above,
for every integer $`n\ge1,L\ge6`$. Put

```math
\ell=\log_2(n+2),\qquad
F(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

One circuit, with absolute constants and full error at most $`2^{-L}`$, gives

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

T-count is worst-case optimal at every eligible width. Count and depth
match simultaneously through $`b\le\sqrt{NL/F(n,L)}`$, whenever
the interval above the literal lower threshold is nonempty. In the
high-precision regime $`L>N/n^2`$, the stronger depth bound is
$`O(NL/b^2+L)`$, with the same count and Clifford budgets.

T-depth allows arbitrary Clifford interlayers and differs from total depth.
The all-precision constructions in B–D use $`G=O(NL)`$ elementary Clifford
gates. The [proof map](docs/README.md#d-simultaneous-frame-t-count-and-t-depth)
connects R79 to its native decoder, terminal schedules and
[uniform baseline](docs/UNIFORM_PRECISION_DEPTH.md).

<a id="two-open-resource-gaps"></a>
<a id="linear-count-endpoint-and-remaining-depth-gap"></a>

<a id="endpoint-count-claim-and-remaining-depth-gap"></a>

<a id="two-clean-count-theorem-and-remaining-depth-gap"></a>

<a id="one-clean-count-frontier-and-remaining-depth-gap"></a>

## Additional resource points and remaining depth gap

The strongest unconditional research upper, [R96](research/depth/ALL_PRECISION_LINEAR_DEPTH.md),
improves the general depth bound to

```math
D_T=O(NL/b^2+L+\ell),\qquad
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL).
```

It retains one clean qubit, $`b\ge17B`$, every integer $`n\ge1,L\ge6`$,
and full initialized-isometry error at most $`2^{-L}`$. Selected Result D
remains R79. At sufficient $`b=\Theta(\sqrt{NL}+L+n)`$, R96 gives
$`D_T=O(L+\ell)`$. **Optimal general T-depth remains open:** at fixed
accuracy and square-root dirty width, the gap is $`\Omega(1)`$ versus
$`O(\log n)`$.

The [resource frontier](docs/OPEN_PROBLEM.md) owns the complete bounds and
scoped lower results; [WORKSPACE](WORKSPACE.md) identifies the first unpaid
operation. The [publication scope](manuscript/PUBLICATION_SCOPE.md) retains
the distinct [two-clean endpoint](research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint),
[diagonal and multiplexor corollaries](docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors),
and smaller-workspace allocations. Earlier constructions and useful components
remain in the [research catalogue](research/README.md#large-workspace-t-depth).

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
