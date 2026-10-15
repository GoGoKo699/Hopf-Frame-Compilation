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
| [Count and large-width depth](docs/OPEN_PROBLEM.md) | Read the one-clean count frontier, sharp endpoint, and sublinear-depth theorem |
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
also matches: $`T^\star=\Theta(N+L)`$. Its proof charges growing
precision sources, borrowed-signal recovery, exact dirty-bank queries, and
the layerwise baseline covering low precision and bounded dimension.
The [diagonal and multiplexor corollaries](docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
retain their own workspace reservations.

**D. Optimal T-count with all-precision T-depth bounds.** For every
$`n\ge1`$ and integer $`L\ge6`$, [one clean qubit suffices](research/depth/ONE_CLEAN_DEPTH.md) at these
family-specific dirty reservations:

| Target | Dirty reservation | Proof |
|---|---|---|
| Complete real frame | $`b\ge17(L+n+7)`$ | [One-clean theorem](research/depth/ONE_CLEAN_DEPTH.md) |
| Phase-dressed complex magnitude frame | $`b\ge17(L+n+8)`$ | [One-clean phase dressing](research/depth/ONE_CLEAN_DEPTH.md) |

Each theorem gives one circuit with absolute constants. Define

```math
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
```

Then

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

T-count is worst-case optimal at every eligible width. Count and depth
match simultaneously through $`b\le\sqrt{NL/F(n,L)}`$, whenever
the interval above the literal lower threshold is nonempty. In the
high-precision regime $`L>N/n^2`$, the stronger depth bound is
$`O(NL/b^2+L)`$, with the same count and Clifford budgets.

The error is at most $`2^{-L}`$ in the full initialized-isometry norm.
Source preparation, queries, actual inverses, and work return are charged.
T-depth allows arbitrary Clifford interlayers and is distinct from total
depth. The complex extension charges its literal leaf-phase diagonal and
composes both factors on the same workspace. Unrestricted depth optimality
remains open. R72 charges native one-clean words, moves terminal rejection
work into the protected logical bank, and proves its own width and precision
splices. The earlier [real](research/depth/ALL_PRECISION_DEPTH.md) and
[complex](research/depth/COMPLEX_FRAME_DEPTH.md) two-clean theorems retain
their stated proofs and allocations; their
[uniform baseline](docs/UNIFORM_PRECISION_DEPTH.md) keeps its separate
precision bounds and matching intervals.

All upper constructions in B–D use $`G=O(NL)`$ elementary Clifford gates.
The [full statements](manuscript/PUBLICATION_SCOPE.md) retain the additional
regimes, complex extensions, and older schedules useful at smaller reservations.

<a id="two-open-resource-gaps"></a>
<a id="linear-count-endpoint-and-remaining-depth-gap"></a>

<a id="endpoint-count-claim-and-remaining-depth-gap"></a>

<a id="two-clean-count-theorem-and-remaining-depth-gap"></a>

<a id="one-clean-count-frontier-and-remaining-depth-gap"></a>

## Additional resource points and remaining depth gap

Result C includes every column, literal phase, actual inverse, and workspace
return in the common full-isometry contract. Its matching statements retain
the displayed banked threshold and exact-base precision condition.
The following distinct allocations remain useful:

| Resource point | Current guarantee |
|---|---|
| Sharper complex endpoint, $`\eta=2^{-N},a=2,b=N+n+7,n\ge3`$ | $`2N-1\le T^\star\le O(N)`$ for both complete real and phase-dressed complex magnitude frames; saves one dirty qubit relative to the one-clean complex allocation |
| Fixed accuracy, $`a=1`$, sufficient $`b=\Theta_\eta(\sqrt N)`$ | For complete real frames, worst-case optimal count is attained with $`D_T=O_\eta(\sqrt n\log(n+2))=o(n)`$ and $`G=O_\eta(N)`$; optimal depth remains open |

The [sharper two-clean endpoint](research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
uses a packed literal diagonal. Result B covers allocations with sufficient
clean work even when Result C's dirty reservation is unavailable. Result D
supplies simultaneous count/depth guarantees on its own circuit; C adds no
T-depth bound. The retained two-clean real-frame [compressed-column special case](research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
replaces the logarithm above by $`\log_2(n+2)`$: for some absolute $`c>0`$,
$`6\le L\le cn`$ and every $`b\ge17(L+n+7)`$, two clean
flags suffice for

```math
T=O(\sqrt{NL}+NL/b),\qquad G=O(NL),\qquad
D_T=O\!\left(NL/b^2+\min\{nL,(\sqrt n+L)\log_2(n+2)\}\right).
```

At sufficient $`b=\Theta(\sqrt{NL})`$, this retains sublinear
T-depth while error $`2^{-L}`$ shrinks, whenever $`L=o(n/\log(n+2))`$.
The [resource statement](docs/OPEN_PROBLEM.md) gives its matching interval
and the remaining optimal-depth question; the [research index](research/README.md)
retains supporting constructions and the limits of specific interfaces.

The [constant-depth obstruction](research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md),
R73, gives an exact zero-clean error floor $`c_0\approx0.169714`$ for an
explicit complete frame at T-depth at most three, regardless of dirty
width. A Clifford attains this floor; integer $`L\ge3`$ requires at least
four layers. The [four-layer operator-rigidity theorem](research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md),
R75, strengthens the zero-clean conclusion: four layers have error at
least $`\rho/2\approx0.135299`$, where
$`\rho=\sin(\pi/8)/\sqrt2`$. Thus integer $`L\ge3`$ requires
at least five layers, with arbitrary dirty width and Clifford interlayers.

The [four-layer crossing-rank theorem](research/depth/FOUR_LAYER_CROSSING_RANK.md),
R74, shows that, at each fixed clean budget, four-layer approximations
to the same witness at error $`2^{-L}`$ require minimum middle active
symplectic crossing rank $`\Omega_a(\sqrt L)`$. Its explicit elementary
bound also supplies a fixed-clean three-layer accuracy threshold.
With one clean qubit, R75 forces exact terminal clean-sector preservation
below $`\rho/2`$. The [one-clean midpoint theorem](research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md),
R76, gives the same error floor when a four-layer circuit admits a
two-plus-two split whose encoded clean observable is Pauli, hence five
layers are required at $`L\ge3`$ in that subclass. A Hermitian-Clifford
midpoint gives a positive explicit threshold without a crossing-rank
bound. The unrestricted non-Clifford midpoint case remains open.
Optimal depth remains an open problem.

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
