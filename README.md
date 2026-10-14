# Exact and Fault-Tolerant Compilation of Hopf Differential Frames

State preparation fixes one column of a unitary. The inverse-frame Hopf
gradient decoder also uses designated coordinate-frame columns. This
repository compiles that prescribed **complete operator**, with explicit
costs for precision, initialized workspace, and borrowed workspace.

The principal results are A–D below, accompanied by a precision-uniform
one-clean count theorem, optimal banked count, and a sharper two-clean
high-precision endpoint. Their
hypotheses and retained corollaries are collected in the
[publication scope](manuscript/PUBLICATION_SCOPE.md).

<p align="center">
  <img src="assets/state-vs-frame.svg" width="900" alt="State preparation fixes one column; complete-frame compilation also preserves the designated coordinate-frame columns." />
</p>

| Read next | Purpose |
|---|---|
| [Proof roadmap](REVIEW.md) | Follow A–D and the selected one-clean count frontier |
| [Documentation map](docs/README.md) | Find authoritative proofs and their dependencies |
| [Verification](docs/VERIFICATION.md) | Inspect executable evidence and its limits |
| [Count theorem and depth question](docs/OPEN_PROBLEM.md) | Read the one-clean count frontier, sharp endpoint, and remaining large-width depth target |
| [Research catalogue](research/README.md) | Find the endpoint proof dependencies, earlier constructions, and scoped limits |

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

<a id="two-open-resource-gaps"></a>
<a id="linear-count-endpoint-and-remaining-depth-gap"></a>

<a id="endpoint-count-claim-and-remaining-depth-gap"></a>

<a id="two-clean-count-theorem-and-remaining-depth-gap"></a>

## One-clean count frontier and remaining depth gap

For $`n\ge3`$ and every integer $`L\ge6`$, the selected
[precision-uniform theorem](research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
compiles both frame families to error $`2^{-L}`$ with one clean qubit
and universal constants:

| Family | Base dirty reservation B | T-count at $`b\ge B`$ |
|---|---|---|
| Complete real frame | $`L+n+7`$ | $`O(N+L)`$ |
| Phase-dressed complex magnitude frame | $`L+n+8`$ | $`O(N+L)`$ |

Both constructions have $`G=O(NL)`$ elementary Clifford gates.
At exactly these base reservations and $`L\ge\max\{6,n\}`$,
the count is worst-case optimal: $`T^\star=\Theta(N+L)`$.
For every $`L\ge6`$ and $`b\ge2B`$, the
[banked theorem](research/endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier)
attains the complete worst-case count frontier with the same one clean qubit:

```math
T^\star=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O(NL).
```

The growing-source proof separates precision from logical mode allocation.
A borrowed signal carries conjugate source sectors; their recovery agrees
on the real target. Exact word-bank queries reduce lookup cost while
returning arbitrary dirty banks. Literal leaf phases compose afterward
using the sole clean flag. Every column, phase, actual inverse, and
workspace return is covered by the common full-isometry contract.

| Separate resource point | Current guarantee |
|---|---|
| Sharper complex endpoint, $`\eta=2^{-N},a=2,b=N+n+7,n\ge3`$ | $`2N-1\le T^\star\le O(N)`$ for both complete real and phase-dressed complex magnitude frames; saves one dirty qubit relative to the uniform complex allocation |
| Fixed accuracy, $`a=2`$, sufficient $`b=\Theta_\eta(\sqrt N)`$ | For complete real frames, optimal count $`T=\Theta_\eta(\sqrt N)`$ is attained with $`D_T=O_\eta(n)`$; unrestricted depth lower bound remains $`\Omega(1)`$ |

The [earlier two-clean theorem](research/endpoint/GROWING_SOURCE_COMPILER.md#10-precision-uniform-two-clean-count)
and [sharper endpoint](research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
retain their own allocations. These selected count results supplement the
unchanged A–D statements. They establish no new T-depth guarantee or
general native emitter, and the base-width upper bound is not an
optimality claim at every smaller precision and dirty allocation.

The [endpoint statement](docs/OPEN_PROBLEM.md) records the remaining
large-width depth question. The [research index](research/README.md)
preserves the earlier promised-family compilers and the limits of specific
source, query, and factorization interfaces. Those restrictions retain
their stated hypotheses.

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
