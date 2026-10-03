# Related work and contribution boundary

[← Source map](SOURCE_MAP.md) · [Complete narrative](../REVIEW.md) · [Landing page →](../README.md)

The paper joins three technical lines around the same prescribed frame:

1. exact state preparation with arbitrary clean workspace;
2. Hopf coordinates and inverse-frame gradient readout;
3. precision-aware Clifford+T synthesis with clean and dirty workspace.

Borrowed-workspace and controlled-unitary constructions explain the strict-zero
schedule.  This page places the ingredients by mathematical role and keeps the
contribution narrow. The exact size–depth and approximate T-count results use
different circuit models; they do not assert simultaneous optimality of one
circuit. The [fault-tolerant theorem](FAULT_TOLERANT_COMPILER.md) states the
second model and its workspace conditions separately.

<p align="center">
  <img src="../assets/literature-lineage.svg" width="940" alt="The all-workspace state-preparation line and the Hopf differential-frame line meet in the prescribed-completion compiler." />
</p>

## 1. Exact state preparation and the target frontier

The state-preparation line develops exact circuits whose depth decreases as
clean ancillary workspace increases.  Its uniform all-workspace conclusion is

```math
S_{\mathrm{QSP}}(n,m)=\Theta(2^n),
```

```math
D_{\mathrm{QSP}}(n,m)
=\Theta\left(n+\frac{2^n}{n+m}\right).
```

The active compiler framework in this repository is P. Yuan and S. Zhang,
*Quantum* **7**, 956 (2023).  We use its exact results for:

- ancilla-free multi-controlled X;
- uniformly controlled one-qubit gates at arbitrary workspace;
- coherent CNOT copy–uncopy;
- the optimal state-preparation comparison frontier.

The published article corresponds to `arXiv:2202.11302v2`; the imported
statements were also checked in v3.

The earlier paper by Sun, Tian, Yang, Yuan, and Zhang established the preceding
time–space landscape and is retained as the historical predecessor and original
source of selected primitives.  For the present theorem, the later uniform
framework is sufficient for all workspace regimes.

The connection to the Hopf problem is constructive rather than automatic.  Once
the Hopf completion is written as addressed complete-operator layers, the
state-preparation primitives can be adapted while preserving its designated
columns.

## 2. Uniformly controlled gates and the Möttönen route

Möttönen and coauthors introduced the uniformly controlled rotation language
for state preparation.  Bergholm and coauthors developed the corresponding
general uniformly controlled one-qubit gates.

The earlier Hopf-QBP compiler follows this route: compute the shared lower-
suffix-zero predicate into one clean flag, then apply a smaller multiplexed
rotation selected by the prefix and flag.  This already shows that the complete
Hopf frame can be compiled with linear asymptotic gate count and constant
additional workspace.

At strict zero workspace, inserting identity blocks into one full-width
multiplexor is not enough.  Although the logical block table is sparse, the
standard Walsh/Gray-code angle transform generically repeats nonzero
coefficients across all suffix frequencies.  The physical angle table is not
sparse in the required sense.

The strict-zero schedule instead retains a narrow UCG.  One original suffix data
qubit carries the predicate temporarily, and two half-angle width-$`d+2`$ UCGs
produce the full addressed depth while restoring that qubit exactly.

## 3. Borrowed workspace and controlled-unitary roots

Several established techniques are close to the strict-zero circuit.

- Barenco and coauthors use roots and conjugation identities in controlled-
  unitary decompositions.
- Later borrowed-qubit constructions use unknown-state logical wires while
  restoring them exactly.
- Conditionally clean ancillas and toggle-detection patterns organize
  cancellation across different original values of a borrowed wire.
- Ancilla-free multi-controlled $\mathrm{SU}(2)$ constructions give related individual
  controlled-rotation primitives.

For one Hopf prefix, put

```math
C_p=R_y(\theta_p/2).
```

The identities

```math
C_p^2=R_y(\theta_p),
\qquad
XC_pX=C_p^{-1}
```

make the original zero-suffix sector square to the required rotation and make
the unwanted original value of the borrowed bit cancel.

The borrowed bit is not an extra dirty ancillary wire.  It is part of the
logical input, the desired operation depends on its original value, and the
proof must restore it on arbitrary entangled inputs.

### Narrow strict-zero contribution

The repository does not claim the invention of borrowed qubits, conditional
cleanliness, toggle detection, controlled-unitary roots, Pauli conjugation, or
UCGs.

The Hopf-specific statement is:

> One original suffix data qubit serves as a restored in-place predicate
> carrier, reducing every addressed Hopf depth to two total-width-$`d+2`$ UCGs
> and linear predicate toggles.  Summed over the tree, this gives the optimal
> strict-zero complete-frame frontier.

## 4. Sparse and restricted multiplexors

[Xu et al., “A Unified Framework for Optimizing Uniformly Controlled
Structures in Quantum Circuits”](https://arxiv.org/abs/2512.08675) study
restricted UCGs and synthesis savings from limited control participation.
The addressed Hopf layer is related, but its sparsity has a particular form:
one long all-zero suffix predicate is shared by every prefix-selected rotation.

The borrowed-suffix factorization converts this logical predicate into smaller
physical UCG width rather than attempting to prune the full-width Möttönen angle
table.  That is the structural distinction relevant to the strict-zero bound.

## 5. Coherent routing and workspace parallelism

For large workspace, the Hopf tail has the exact direct-sum form

```math
R_t^{(n)}
=\bigoplus_{r=0}^{2^t-1}W_s^{(r)}.
```

The compiler realizes this form by standard reversible ingredients:

- balanced CNOT fanout of the prefix bits;
- Fredkin routing of one suffix-token block;
- disjoint token-controlled subtree frames;
- exact uncomputation and inverse routing.

The state-preparation framework supplies the coherent-copy and UCG primitives.
The Hopf-specific part is the tree-cut identity, the register layout, the
route–operate–unroute schedule, and the shared workspace envelope for prefix and
tail.

The implementation gives the router explicitly and tests arbitrary complex
inputs in which the prefix and suffix are entangled.

## 6. Hopf coordinates and inverse-frame gradients

The first Hopf work develops the balanced binary chart, inverse map, diagonal
metric, and coordinate directions.  The second develops the addressed
state-and-marker frame, global magnitude record, direct phase record, and
checkpoint interface.

The compiler target is

```math
W_{\mathbb R}|0^n\rangle=|\psi\rangle,
\qquad
W_{\mathbb R}|\lambda(j)\rangle=|e_j\rangle.
```

For unrestricted angles,

```math
\partial_{\theta_j}|\psi\rangle
=a_j|e_j\rangle,
\qquad
g_{j,j}=a_j^2.
```

On the canonical domains, $a_j$ is nonnegative.  At a singular coordinate, the
raw derivative vanishes while the complete parameter tuple still selects an
orthogonal marker-frame continuation.

The phase-dressed complex magnitude frame is

```math
W_{\mathbb C,\mathrm{mag}}
=D_{\mathrm{ph}}W_{\mathbb R}.
```

The leaf-phase coordinates remain a separate direct stream.

The synthesis theorem proves that these prescribed magnitude-frame operators,
not only their first columns, attain the all-workspace state-preparation
frontier.

## 7. Comparison by mathematical role

| Line of work | Object or result | Role used here | Additional Hopf structure |
|---|---|---|---|
| Möttönen; Bergholm | multiplexed rotations and general UCGs | block-diagonal one-qubit-gate language | addressed depths are reduced to narrow UCGs rather than full-width sparse tables |
| exact QSP time–space tradeoffs | arbitrary state preparation with clean workspace | target frontier and compiler primitives | designated frame columns must survive the adaptation |
| controlled-unitary roots and borrowed workspace | exact conditioned operations with limited ancillas | cancellation and restored logical carriers | the predicate includes the borrowed suffix bit's original value |
| restricted UCGs | repeated or sparse block structure | neighboring structural context | one shared zero-suffix predicate is factored through a restored logical bit |
| reversible routing | CNOT fanout, Fredkin networks, uncomputation | coherent branch movement and parallel action | realizes the Hopf tail direct sum within the prefix workspace envelope |
| Hopf chart | state map, inverse, metric, coordinate directions | specifies the structured columns | yields the addressed tree layers |
| Hopf QBP | inverse-frame magnitude record and direct phase record | makes the completion operational | requires frame-safe rather than state-column compilation |

## 8. Contribution at theorem level

The exact compiler contribution is the following chain.

1. Identify complete frame safety as the compiler contract used by the global
   inverse-frame record.
2. Give an exact state-column counterexample.
3. Express the frame as zero-suffix-addressed tree layers.
4. Close the zero-workspace endpoint with the borrowed-suffix echo.
5. Build a clean binary–one-hot prefix decoder and explicit coherent router.
6. Match the state-preparation size–depth frontier for every $m\geq0$.
7. Transfer the result to the phase-dressed complex magnitude frame.
8. State the QBP consequence under matched program and output-task conventions.

The theorem is specific to the structured Hopf completion.  It does not imply
that an arbitrary family of prescribed unitary columns can be implemented at
state-preparation cost.

## 9. Fault-tolerant lineage and the precision register

Let $N=2^n$ and let $L$ denote the requested number of accuracy bits. Three
established results are particularly close to the fault-tolerant construction.

**Gosset–Kothari–Wu (GKW).** Their Theorems 1.1–1.2 give the optimal unrestricted
T-count $\Theta(\sqrt{NL}+L)$ for state preparation and diagonal synthesis.
Their Appendix B already treats complete single-qubit multiplexors and allocates
error geometrically across state-preparation layers, obtaining
$O(\sqrt{NL}+nL)$. These are imported benchmarks and methods. A state-preparation
theorem alone does not specify the remaining columns of the Hopf frame, and an
unrestricted-ancilla theorem does not supply a prescribed clean/dirty budget.
Their Theorem 4.2 supplies the ancilla-independent lower bound used through the
frame's diagonal subfamily.

**Low–Kliuchnikov–Schaeffer (LKS).** Section 2, Table 2, Figure 1(d), and Eq. (8)
provide dirty-assisted SelectSwap lookup with exact bank restoration. For an
$m$-bit table with $D$ entries, its T count is $O(D/\lambda+m\lambda)$ with
$m\lambda$ dirty bank qubits and $O(\log D+m)$ clean wires. Appendix C also
gives exact Boolean XOR constructions using dirty auxiliaries. The basis-state
identity in Eq. (8) extends to reference-entangled inputs by linearity. We do
not claim dirty lookup or its cancellation identity as a contribution. Its
dirty bank does not replace the clean output word needed by an instruction-word
interpreter. Section 5 supplies the finite-width counting method.

Appendix C also separates indicator parallelism from count-efficient bank
selection. The [parallel lookup proof](PARALLEL_DIRTY_LOOKUP.md) gives a
literal exact realization by conjugating one X with the inherited bank
router, with scratch-free indicator work and linear address T-depth.
The predicate schedule uses Khattar–Gidney
[Section 5.4](https://arxiv.org/html/2407.17966v1#S5.SS4), allocating its two
dirty helpers only between completed queries. Its role here is to
compose the inherited lookup idea with complete two-clean real frames,
retaining their T-count while reducing T-depth. It does not claim a new
general lookup tradeoff or equally small Clifford depth.

The [partial-batch schedule](BATCHED_DIRTY_LOOKUP.md) reuses these
indicator and linear-echo mechanisms with a smaller simultaneous pool.
LKS Appendix C, Theorem 2 states a count-efficient tradeoff with
$`O(\lambda\sqrt{Qm})`$ dirty workspace, $`\lambda\ge1`$, so its
displayed range begins at square-root-scale width. Partial batches
retain the count-efficient bank choice below that scale and supply the
intermediate-workspace T-depth ledger and full return contract for Hopf
composition. The underlying lookup mechanism is inherited; no priority
claim over general lookup tradeoffs is made. Motlagh–Pocrnic,
[arXiv:2605.20334v1, Sections II.2–II.4](https://arxiv.org/html/2605.20334v1),
improve Toffoli-count prefactors through SelectCopy and sequential
output-bit packets. Their displayed constructions retain initialized
outputs and do not state the simultaneous $`B^{-2}`$ T-depth dependence
audited here. The batch guards reserve F8's two dirty helpers during
the live query, separately from its occupied banks and indicators.

The [amortized refinement](AMORTIZED_DIRTY_LOOKUP.md) combines the same
dirty echo with a multiplexed family of linear shears. The dirty tree
traversal and cancellation are established techniques; see
[Khattar–Gidney, arXiv:2407.17966v2, Sections 4, 7.1 and 7.4](https://arxiv.org/html/2407.17966v2).
[Boyd, arXiv:2312.00696v2, Section III and Appendix A](https://arxiv.org/html/2312.00696v2)
already groups commuting SELECT/QROM operators and changes Clifford
bases to expose parallel action, using clean address-fanout registers.
[Kim–Laakkonen, arXiv:2512.24982v1, Theorems 3 and 5](https://arxiv.org/pdf/2512.24982v1)
give constant-depth control of CNOT and Clifford circuits without
ancillas or measurements. The local rank-reduction lemma specializes
that established possibility with explicit rank-sensitive T-count and
Clifford-count bounds. The new composition moves the low-address
indicator outside the high-address traversal and supplies its all-dirty
return identity and complete Hopf resource ledger. It removes the batch
logarithm and proves matching count and depth in an explicit
precision/workspace range, including fixed and inverse-polynomial
accuracy in N. The variable-accuracy extension uses the same native
query and the existing layer error bounds.
No priority is claimed for dirty iteration, commuting-operator grouping,
or constant-depth controlled Clifford circuits.

**Bausch.** Equations (4) and (6) of *Fast Black-Box Quantum State Preparation*
already use a geometric precision register and a bit oracle addressed by the
data index and precision position. Section 2.3.3 includes a capped geometric
source. Its explicit exact preparation uses a linear-size unary temporary,
then binary conversion and routing cleanup. Consequently, logarithmic output
width there does not by itself establish logarithmic peak clean workspace.
The compact Gray construction used here supplies the specific simultaneous
resource guarantee: exact $O(L)$ T count, $O(\log L)$ peak clean workspace,
ordinary binary labels, and restored temporary work. The geometric weighting
and bit-oracle idea are credited to Bausch; no strict gate-count improvement
over his routing variants is asserted.

The addressed SU(2) compiler combines this source with a constant-size Pauli
linear combination and oblivious amplitude amplification. It is a
space-efficient application of these ingredients, with all preparations,
adjoints, reflections, and one-bit tables charged. It is not a new general
principle of dyadic sampling or amplitude amplification.

## 10. What the full-frame T theorem adds

The main additional construction is the uniform-in-precision composition of
complete frame corrections. A shared geometric source supplies approximate
shift amplitudes on accepted branches; the source is retained between stages.
The residual decomposition, failure tracking, and one final amplification
control the full initialized isometry while returning dirty workspace exactly.
The failure counter uses the established block-product compression gadget of
[Low–Wiebe, Lemma 13](https://arxiv.org/html/1805.00675v2), in the form explained
by [Fang–Lin–Tong, Section 2.4, Lemma 3 and Appendix D](https://arxiv.org/html/2208.06941v2).
It coherently retains failed histories while reusing block work. Neither
logarithmic history storage nor final amplification is a new general
mechanism here; the additional ingredients are the residual representation,
shared precision source, and charged frame-specific schedule.
Under the conditions in the [theorem](FAULT_TOLERANT_COMPILER.md), this gives

```math
T=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{n+a+b}\right),
```

where $a$ and $b$ count clean and dirty ancillary qubits. The theorem retains
its sufficient clean-workspace reservation. A separate retained corollary for
all clean budgets, under an additional restriction on precision and dirty
width, is proved in the [borrowed-workspace appendix](BORROWED_WORKSPACE_COMPILER.md).
The lower bounds follow from GKW and LKS with the stated reductions; the
contribution is the upper construction and the resulting match in those
regimes. The all-clean-budget extension also uses the earlier borrowed-work
compiler; it is a compiler splice, not an independent lower-bound technique.

The simpler direct sampler already costs
$O(\sqrt{NL}+nL+NL/(n+a+b))$ with logarithmic precision workspace. For example,
at $L=n^2$, $a=\Theta(n)$ with sufficient headroom, and
$b=\Theta(n\sqrt N)$, it already attains $\Theta(n\sqrt N)$ T gates.
This illustrates the small-clean sampler; it is not an additional improvement
due to shared-source composition. Eliminating repeated precision charges is
the reason for the stronger uniform-precision construction.

Tan's published Theorem I.1 and Lemma IV.1 treat general unitary synthesis and
batched controlled gates. Remark IV.2 uses zeroed precision-length instruction
registers; it does not directly give the clean/dirty contract above. The
Yuan–Zhang depth theorem uses arbitrary continuous one-qubit gates, so it is not
a T-count comparison. Recent sparse-QROM and sparse-state bounds address
different input families. These distinctions delimit the use of each result;
they do not establish priority by absence of a matching theorem.

At a fixed parameter value, a full-isometry error can bound the bias of the
[bounded-score QBP estimator](QBP_APPROXIMATION.md).
This does not differentiate the compiled Clifford+T word as a function of the
parameters and does not prove optimal gradient-query complexity.

## 11. Operator sources and the one-clean refinement

The [operator-source compiler](OPERATOR_SOURCE_COMPILER.md)
uses two initialized flags and an arbitrary dirty bank. Full-space
anticommuting loaders are established by
[Kerenidis–Prakash](https://arxiv.org/html/2202.00054v2), Definitions 4.4/4.6
and Theorem 4.9; their scalar-overlap product algebra also appears in
[Chee et al.](https://arxiv.org/pdf/2301.07477), Appendix C. Dirty SelectSwap
is credited to LKS, and normalization-two amplification to
[Berry et al.](https://arxiv.org/pdf/1412.4687), Eqs. (11)–(15).

The local construction specializes the loader to native geometric weights,
programs its signs through whole-word dirty queries, and combines a scalar
block, suffix echo, and full-output amplification into the prescribed Hopf
frame. With $`b\ge L+n+7`$ it gives $`T=O(N+nL)`$; with
$`b\ge2(L+n+7)`$ its banked form gives
$`T=O(\sqrt{NL}+nL+NL/b)`$. Both use $`O(NL)`$ Clifford gates.
The operator core returns approximately within the full-isometry bound;
lookup banks and selectors return exactly. Literal diagonal compilation
also gives the phase-dressed magnitude-frame corollary with its own error
allocation and dirty-width threshold.

The [conditional-suffix compiler](CONDITIONAL_SUFFIX_COMPILER.md) gives
$`O(N+L\ell_*(n))`$ T gates and $`O(NL)`$ Clifford gates for real
frames at $`n\ge1`$, $`L\ge6`$, $`a=2`$, and $`b\ge L+n+7`$.
Here $`\ell_*(n)=1+\log_2^*(n+2)`$, with repeated base-two logarithms
stopping at a value at most one. Star residuals and streamed coarse programs
use logical suffix qubits as clean work only in the active sector. At
$`L=N`$, the upper bound is $`O(N\ell_*(n))`$; the $`\Omega(N)`$
lower bound remains unchanged.
With $`b\ge2(L+n+7)`$, the grouped banked bound is
$`O(\sqrt{NL}+L\ell_*(n)+NL/b)`$, still with $`O(NL)`$ Clifford gates.
These are explicitly charged constructions in their stated regimes;
loaders, dirty lookup, overlap algebra, and amplification retain their
earlier attribution.

Hierarchical block organization and in-place index updates are also
established. [Nguyen–Kiani–Lloyd](https://arxiv.org/pdf/2201.11329),
Section 4.1 and Appendix C.1, use a separate initialized local-vector
register; their low-rank construction in Section 4.3 uses factor-state
oracles. [Setty](https://arxiv.org/pdf/2508.21667), Theorem 1 and Section 2,
includes in-place logical shift/delete operations while retaining initialized
data-label/PREP work. The local ancestor blocks instead use the logical word
as the coefficient address, a dirty scalar precision source, and explicitly
charged streaming lookup. The contribution is this complete-workspace
realization for the prescribed residual structure, not a new generic LCU
principle or an assumed QRAM interface.

### One initialized flag

The [one-clean compiler](ONE_CLEAN_COMPILER.md) retains the established
Clifford-loader and dirty-lookup ingredients. Its additional construction
conjugates one scalar-source block by another, routing their rejected
components through different logical Pauli operators. The Clifford
anticommutator reduces the accepted block to a real rotation with two
programmable scalar coefficients. A head coefficient and two geometric
tails supply those coefficients. The standard paired-Majorana algebra and
its rotations are recalled in
[Bravyi, Section II, Eqs. (2)–(5)](https://arxiv.org/pdf/quant-ph/0404180).
The local source packs both tails into essentially one dirty qubit per
accuracy bit. The
contribution is this native source and conjugated block word, with its
complete-input error and workspace accounting, not Majorana operators,
Clifford-algebra overlap identities, or source conjugation in general.

At the fixed amplitude $`\beta=\sin(\pi/10)`$, two amplification
iterations use five forward or inverse block calls. The amplification
mechanism is inherited from
[Brassard–Høyer–Mosca–Tapp, Section 2, Eqs. (7)–(8)](https://arxiv.org/pdf/quant-ph/0005055)
and the oblivious block setting above. The local proof verifies its literal
phase and complete-isometry error for the programmed source, including dirty
core return; it does not claim a new amplitude-amplification principle.

The real-frame baseline and grouped bounds therefore require only
$`a\ge1`$, with the same $`b\ge L+n+7`$ dirty reservation. Their
banked refinements retain the doubled dirty threshold. A flag-only Pauli
routing also gives one-clean literal diagonals and, through the retained
certified Euler decomposition, complete one-target U(2) multiplexors.
Sequential composition with the literal diagonal also gives the grouped
phase-dressed complex magnitude-frame corollary, at $`b\ge L+n+8`$ or
its doubled banked threshold. This is closure under the existing
complete-isometry contract, with no new phase-gradient columns.

The real-rotation word also commutes with X on its amplification signal.
The [full-space symmetry argument](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
therefore permits borrowing that signal, giving the zero-clean layerwise
bound $`T=O(N+nL)`$, $`G=O(NL)`$ at $`b\ge L+n+7`$.
The scalar-phase word does not share this symmetry, and the grouped proof
retains its external clean predicate. These corollaries do not lower the
proved full-frame endpoint T count or extend the separately proved
two-clean T-depth schedules.

The same primitive prices the
[weighted forward component](WEIGHTED_TRANSPORT_BLOCK.md) at
$`O(N+nL)`$ T gates using only borrowed synthesis work alongside its
logical dilation signal. Approximate return of the core and borrowed
amplification signal is included in the full-operator error. The earlier
$`O(L\sqrt N)`$ route remains useful for its exact dirty-work return.
Neither component bound removes the full-frame endpoint gap.

The [complete residual assembly](RESIDUAL_ASSEMBLY.md) incorporates the
diagonal and forward term into one affine tree map, then selects it or the
reverse dilation using two clean flags. The scalar elimination uses
standard Schur complements (F18 in the source map); the local contribution is the affine
recursion, complete tree-mode allocation, and native controlled-branch
accounting. The prepared two-term linear combination and normalization-two
amplification are inherited from
[Berry et al., Eqs. (7)–(15)](https://arxiv.org/pdf/1412.4687).
Borrowed synthesis work and full-operator error estimates make the selected
branches usable while both flags are occupied. This gives
$`T=O(N+nL)`$ and $`G=O(NL)`$ at $`a=2`$,
$`b\ge L+n+7`$; it is not a new generic LCU method or an improvement
over the grouped frame bound. The constant-call reduction applies to the
new affine and reverse dilations, not to an uncontrolled opaque forward
block alone.

The [source map](SOURCE_MAP.md) gives exact theorem numbers and local consumers.
The comparisons identify dependencies and specific additional constructions;
they do not certify priority.

## 12. Contemporary comparisons and the broader compiler contribution

The following comparisons were checked against primary sources on
**22 September 2026**. Their precision, initialization, and error contracts
matter as much as their T-count exponents.

**General unitary synthesis.** Tan is no longer the newest generic upper
bound. [Yuan–Zhang–Zi, arXiv:2608.17846v2, Theorem I.1 and Definition II.1](https://arxiv.org/html/2608.17846v2)
give, for $`n+L\le N`$,

```math
T=O\!\left(nN^{5/4}(n+L)^{5/8}\right),
\qquad a=O\!\left(N\sqrt{n+L}\right).
```

Their ancillary qubits are clean. Their error criterion is the complete
initialized-isometry norm, including leakage and literal phase, so that
criterion itself is not a novelty claim here. They treat arbitrary dense
unitaries; the Hopf family has only $`N-1`$ real parameters and has additional
tree structure. The comparison establishes the importance of using that
structure, not an improvement to arbitrary-unitary synthesis.
[Fang–Heunen–Wang, arXiv:2607.12907v1, Theorem 1.2 and Corollary 3.8](https://arxiv.org/html/2607.12907v1)
give an instance-dependent bound in terms of phase-optimized Frobenius
distance to the Clifford group. Their circuit uses initialized ancillas and
their detailed conclusion is in diamond distance. It supplies neither a
two-clean guarantee nor the prescribed frame's worst-case frontier.

**The closest precision and workspace comparison.**
[Yamazaki–Akibue, arXiv:2603.14202v1, Theorem 1 and Section 3](https://arxiv.org/html/2603.14202v1)
sharpen the leading precision coefficient for complete multiplexed SU(2)
gates to $`3L+O(\sqrt{NL})+o(L)`$ for most targets. Their displayed
construction explicitly reserves $`3L+O(n)+o(L)`$ clean ancillas; its
remaining ancillas can be dirty. Theorem 4 removes ancillas with a cost
of order $`NL+Nn`$, under its stated typical-target guarantee. These are
channel-distance results with additional restrictions in the matching
leading-constant lower bound. They do not establish a literal-phase,
constant-clean implementation for every supplied table.

GKW's diagonal proof likewise computes a precision-length instruction word
into an initialized register, then applies that word and uncomputes it. LKS
can replace its lookup scratch with dirty banks, but an arbitrary dirty
output word is not a known instruction word. The operator-source constructions
avoid that instruction register by using full-space operator identities on
the dirty core; the one-clean refinement changes the scalar-block word.

This distinction already gives a result beyond the Hopf family. For **every
literal diagonal** on $`n`$ qubits, the
[one-clean proof](ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors) establishes

```math
a=1,\quad b\ge2(L+n+7),\qquad
T=\Theta\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O(NL).
```

The upper bound includes approximate joint return of the dirty operator
core; lookup banks and selectors return exactly. Since $`n+1+b=\Theta(b)`$
in this range, GKW's diagonal lower bound and LKS-style finite-width counting
give the stated worst-case match. At $`L=N`$ and $`b=2(N+n+7)`$ this is
$`\Theta(N)`$ T count with one clean qubit. This diagonal endpoint is
settled by the retained proof, whereas the corresponding complete-frame
endpoint still has a gap between $`\Omega(N)`$ and
$`O(N\ell_*(n))`$. The broader scientific contribution
is therefore a precision/workspace compiler for a standard operator family,
together with the structured complete-frame extensions. It is not a claim
that the entire constant-clean frame frontier is matched.

The same refinement covers **arbitrary complete one-target U(2)
multiplexors**, with $`N=2^n`$ blocks, one initialized flag, and
$`b\ge2(L+n+9)`$ dirty qubits, at
$`\Theta(\sqrt{NL}+L+NL/b)`$ T gates. Four addressed Euler
factors suffice; their implementation reuses the flag and dirty core with
all return error included. The [one-clean corollary](ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
uses the retained certified matrix-entry procedure through finite approximate
Euler search,
so no nonsingular-chart or exact-zero promise is introduced. This directly
addresses the complete multiplexor task in the GKW and Yamazaki–Akibue
comparisons, while supplying the one-clean guarantee. It does not
improve their leading constants, and its general classical coordinate search
is not claimed efficient. It also does not combine all Hopf tree depths into
one jointly charged precision source.

**Why general compression does not remove the distinction.**
[Vasconcelos–Gilyén, arXiv:2507.07900v2](https://arxiv.org/html/2507.07900v2)
give an ancilla-uncomputation procedure that still queries the original
block encoding on its initialized work. Thus it reduces retained work
between calls, not the initial peak-clean requirement to a constant number
of qubits. Their
exact-product logarithmic ancilla lower bound is for the specified multiple
coherent measurement circuit class. Their approximate compression theorem
requires near-identity supplied dilations, not merely near-identity accepted
blocks. Neither statement proves or resolves our constant-clean frame gap.
[Ma–Joven–Liu, arXiv:2609.11153v1, Theorem 3.1](https://arxiv.org/html/2609.11153v1)
compress a block encoding to at most $`n+2T`$ clean ancillas without raising
T count, possibly changing normalization. This supports counting lower
bounds; it is not a constant-clean unitary implementation theorem.

**Lookup constants.**
[Motlagh–Pocrnic, arXiv:2605.20334v1, Section II](https://arxiv.org/html/2605.20334v1)
improve dirty-QROM constants using SelectCopy and shared work across bit
packets. The output instruction register remains initialized in their
construction. These improvements preserve the asymptotic lookup tradeoff
used here and do not supply the one-clean diagonal compiler. We claim no
improvement to their constant factors.

## References highlighted here

- A. Barenco et al., “Elementary gates for quantum computation,” *Physical
  Review A* **52**, 3457–3467 (1995).
- M. Möttönen et al., “Transformation of quantum states using uniformly
  controlled rotations,” *Quantum Information and Computation* **5**, 467–473
  (2005).
- V. Bergholm et al., “Quantum circuits with uniformly controlled one-qubit
  gates,” *Physical Review A* **71**, 052330 (2005).
- X. Sun, G. Tian, S. Yang, P. Yuan, and S. Zhang, “Asymptotically Optimal
  Circuit Depth for Quantum State Preparation and General Unitary Synthesis,”
  *IEEE TCAD* **42**, 3301–3314 (2023).
- P. Yuan and S. Zhang, “Optimal (controlled) quantum state preparation and
  improved unitary synthesis by quantum circuits with any number of ancillary
  qubits,” *Quantum* **7**, 956 (2023).
- B. Claudon et al., “Polylogarithmic-depth controlled-NOT gates without
  ancilla qubits,” *Nature Communications* **15**, 5886 (2024).
- T. Khattar and C. Gidney, “Rise of conditionally clean ancillae for efficient
  quantum circuit constructions,” *Quantum* **9**, 1752 (2025).
- B. Zindorf and S. Bose, “Efficient implementation of multi-controlled quantum
  gates,” *Physical Review Applied* **24**, 044030 (2025).
- C. Xu et al., [“A Unified Framework for Optimizing Uniformly Controlled
  Structures in Quantum Circuits”](https://arxiv.org/abs/2512.08675),
  arXiv:2512.08675 (2025).
- J. Bausch, [“Fast Black-Box Quantum State Preparation”](https://arxiv.org/pdf/2009.10709v4),
  arXiv:2009.10709v4 (2022), Eqs. (4), (6), and Section 2.3.3.
- G. H. Low, V. Kliuchnikov, and L. Schaeffer,
  [“Trading T gates for dirty qubits in state preparation and unitary synthesis”](https://arxiv.org/html/1812.00954v2),
  *Quantum* **8**, 1375 (2024).
- G. Brassard, P. Høyer, M. Mosca, and A. Tapp,
  [“Quantum Amplitude Amplification and Estimation”](https://arxiv.org/pdf/quant-ph/0005055),
  arXiv:quant-ph/0005055, Section 2, Eqs. (7)–(8).
- D. W. Berry, A. M. Childs, R. Cleve, R. Kothari, and R. D. Somma,
  [“Simulating Hamiltonian dynamics with a truncated Taylor series”](https://arxiv.org/pdf/1412.4687),
  *Physical Review Letters* **114**, 090502 (2015), Eqs. (11)–(15).
- D. Gosset, R. Kothari, and K. Wu,
  [“Quantum state preparation with optimal T-count”](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/),
  *Quantum* **10**, 2168 (2026).
- X. Tan, [“Unitary Synthesis with Fewer T Gates”](https://journals.aps.org/prxquantum/pdf/10.1103/pxhd-9s9q),
  *PRX Quantum* **7**, 033058 (2026).
- T. Li, F. Ou, X. Wang, P. Yao, P. Yuan, and S. Zhang,
  [“Optimal T Counts under Sparsity: from QROM to State Preparation and Block Encoding”](https://arxiv.org/html/2607.28260v1),
  arXiv:2607.28260v1 (2026).

---

[← Source map](SOURCE_MAP.md) · [Complete narrative](../REVIEW.md) · [Landing page →](../README.md)


## 13. State-only preparation and a separate gradient decoder

The [state-only construction](STATE_ONLY_COMPILER.md) uses the retained
borrowed compiler for a coarse circuit, the full-operator rotation primitive
for a single correction table, and standard amplitude amplification.
[Brassard–Høyer–Mosca–Tapp, Section 2, Eq. (8)](https://arxiv.org/pdf/quant-ph/0005055)
gives the one-step success law at accepted amplitude one half. The local
proof supplies the actual coarse residual, two-flag word, coherent branch
extension, and literal dirty-space and error ledgers. At $`L=N`$ its linear
state-preparation T count agrees with the order of
[GKW's unrestricted-ancilla benchmark, Theorem 1.1](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/);
the additional local statement is the two-clean allocation and arbitrary
dirty/reference return. No complete-frame conclusion follows.

The [reference-state QBP protocol](REFERENCE_STATE_QBP.md) changes the
measurement and classical scores. Logarithmic wave-function derivatives
already appear in variational Monte Carlo, for example
[Toulouse–Umrigar, Eq. (45)](https://arxiv.org/pdf/physics/0701039).
Hadamard-test interference is also standard; see
[Mitarai–Fujii, Fig. 1](https://arxiv.org/pdf/1901.00015).
The general overlap identity is not a novelty claim. The local analysis
combines Hopf's one-node-per-depth leaf support with a positive derivative
envelope and the charged state-only compiler. A fixed lower bound on local
branch probabilities yields constant score overhead; arbitrary angles can
require an n-dependent factor in this sufficient sampling bound. This is a
specific alternative for raw gradients, not a generic training speedup or a
priority claim over all importance-sampling protocols.

The subsequent [coarse-frame decoder](COARSE_FRAME_QBP.md) combines the
same standard X/Y interference identity with a coarse full-frame word and
the original Walsh spreading step. Its local contribution is the
constant depth-record bound from coarse marker agreement, the coherent
reference matched to the actual native word, and the charged quantum and
classical reconstruction. Exact correction in the decoder removes coarse
bias without synthesizing a fine inverse frame. No claim of priority for
basis-dependent overlap estimation or classical adjoint accumulation is
made. The result covers all real angle tuples, including singular ones;
the complex phase-gradient stream remains separate.

## 14. Gauge-fixed complex state and gradient extension

The prefix phase cascade and its residual arithmetic-mean global phase
are inherited from [Möttönen–Vartiainen–Bergholm–Salomaa, Section III,
Eqs. (4), (5), and (7)](https://arxiv.org/pdf/quant-ph/0407010v1).
The [complex coarse compiler](COMPLEX_COARSE_COMPILER.md) fixes the sign
and half-angle convention, cancels the common phase for a state-based
task, and applies the retained literal reflection interpreter to each
determinant-one table. Its local result is the actual logical coarse C
with exact dirty return, linear T count, recorded possibly nondiagonal
rows, and the unchanged two-flag fine state-preparation reservation.

The [two-stream decoder](COMPLEX_COARSE_QBP.md) combines that interface
with the existing X/Y coarse-frame reconstruction and direct leaf-phase
records. It supplies the consistent gauge correction, arbitrary-angle
gradient identities, and complete quantum/classical precision ledger.
Neither the standard phase cascade nor global-phase invariance is a new
claim. The result does not compile the prescribed literal common phase of
the complete frame, prove universal gradient-cost improvement, or emit a
general native complex state-preparation package.

## 15. Bounded-input construction and the classical comparator

The [residual-table procedure](RESIDUAL_TABLE_PREPROCESSING.md) uses
ordinary half-angle algebra and the standard three-rotation SU(2)
factorization. Its local statement is a certified coefficient construction
for the particular residual completion: one shared half-phase preserves
the literal branch phase, finite rational decisions cover zero and boundary
cases, and direct sine/cosine programming retains the existing native error
and workspace contracts. No new general Euler decomposition is claimed.

The [bounded-input audit](BOUNDED_INPUT_QBP.md) uses GKW Lemma 2.3 only for
short-word existence. Exhaustive search at coarse precision gives a
conservative polynomial-in-N construction. This must be distinguished
from efficient fine synthesis: [Ross–Selinger](https://arxiv.org/abs/1403.2975v3)
give optimal synthesis with a factoring oracle and prove the efficient
expected runtime without that oracle under a number-theoretic hypothesis.
Neither is assumed in the bounded-input construction. Its residual fine
precision instead comes from explicit source masks.

For explicit Pauli inputs, applying signed permutations and propagating
derivatives backward are standard classical operations. Reverse
differentiation has a much broader established theory; see
[Baur–Strassen, *The complexity of partial derivatives*](https://www.sciencedirect.com/science/article/pii/030439758390110X),
*Theoretical Computer Science* **22**(3), 317–330 (1983).
The local audit provides the division-free Hopf recurrence, certified
dyadic rounding, duplicate-term accounting, and a term-sampling comparator.
These baselines limit the interpretation of the quantum T-count result;
they are not a claim to have invented reverse differentiation or classical
importance sampling. Unknown controlled observables retain their distinct
access model.

## 16. Precision depth and workspace assumptions (2 October 2026)

The [source-depth audit](SOURCE_T_DEPTH.md) separates a restriction of the
geometric source implementation from the unrestricted
[Hopf T-depth problem](T_DEPTH_COMPILER.md#4-lower-bounds-and-the-remaining-depth-gap).
The relevant denominator technique is already present in
[Casas et al., *Matchgate synthesis via Clifford matchgates and T gates*,
arXiv:2602.05425v1, Section III.2.2, Eq. (29)](https://arxiv.org/html/2602.05425v1#S3.SS2.SSS2).
In their Majorana representation, Clifford matchgates permute signed modes,
and one layer of disjoint non-Clifford mode rotations increases the least
square-root-of-two denominator exponent by at most one. Their resulting
depth lower bound concerns exact synthesis within this restricted gate
set. Applying that mechanism to the geometric source supplies a local
obstruction to parallelizing that source in the same representation; it
does not introduce a new general lower-bound method. Arbitrary Clifford
interlayers need not preserve the Majorana span, and the frame compiler
need not implement this source exactly, or use it at all.

There is also an established alternative when precision-sized **clean**
workspace is available. [Vasconcelos, *Depth-Optimal Quantum Compilation*,
arXiv:2609.34659v1, Theorem 8](https://arxiv.org/html/2609.34659v1#S3.SS4.SSS1)
constructs a fully unitary single-qubit rotation approximation with
$`O(L)`$ initialized ancillas, $`O(L)`$ gates, and $`O(\log L)`$
elementary depth, for $`L=\Theta(\log(1/\varepsilon))`$. Its error
includes ancilla return and a stated common phase. Fixed exact Toffoli
decompositions preserve these orders over Clifford+T. This demonstrates
that linear precision depth is not intrinsic to rotation approximation.
It does not give the same construction with two clean qubits and arbitrary
dirty work. Its bounded-arity elementary-depth lower bound also does not
bound T-depth when unrestricted Clifford circuits between T layers are free.

[Kim, *Catalytic z-rotations in constant T-depth*, arXiv:2506.15147v3,
Section 3](https://arxiv.org/pdf/2506.15147), published in *Quantum*
**10**, 2191 (2026), explicitly leaves constant-T-depth rotation using
only clean or dirty ancillas open. The depth-three construction assumes
a prepared nonstabilizer catalyst; the supplied-catalyst resource cannot
be replaced by arbitrary borrowed qubits. Charging catalyst preparation
and its initialized workspace is necessary before composition with the
present compiler. Its final note records subsequent depth-two and
measurement-assisted depth-one refinements. The universal-catalyst
comparison in Section 19 below includes Kim–Laakkonen's later construction.

These comparisons identify usable techniques and their workspace
conditions. The 2 October source-restriction pass yielded no asymptotic
improvement to the unrestricted complete-frame depth bound, and no
matching lower bound. Section 19 records the later unary-source compiler.
The exact source restriction must therefore remain separate from both
optimal approximate Hopf T-depth and the constant-clean T-count endpoint.

The [two-layer obstruction](SHALLOW_SOURCE_OBSTRUCTION.md) uses a
different, full-input argument. [Aaronson–Gottesman,
arXiv:quant-ph/0406196v5, Section III](https://arxiv.org/pdf/quant-ph/0406196v5)
gives the discrete magnitudes of stabilizer overlaps.
[Zhang–Zhang, arXiv:2409.13809v2, Theorem III.1,
Eqs. (10)–(11)](https://arxiv.org/html/2409.13809v2#S3.SS1)
explicitly uses the fact that a T layer conjugates Paulis to Hermitian
Cliffords. Splitting a two-layer transfer coefficient at the middle
Clifford reduces it to a normalized Clifford trace. The local proof
then derives a transfer alphabet and constant approximation gaps for
the geometric sources with any dirty width. These standard ingredients
are attributed; no growing unrestricted depth lower bound is inferred.

The [conditional geometric source](CONDITIONAL_GEOMETRIC_SOURCE.md)
supplies a different positive interface. It prepares the geometric
one-hot state using reversible prefix ORs and native controlled
Hadamards, clears noninvariant tree scratch before that batch, and
uses conditional suffix workspace rather than an externally initialized
precision register. Its preparation, selective reflection, actual
inverse, inactive-sector identity, and return error are charged locally.
This is an explicit construction from standard reversible and
block-encoding ingredients; it does not import a catalyst or claim
generic fast rotation synthesis with arbitrary dirty helpers. The
complete-frame query and suffix-predicate bottlenecks require a separate
composition, supplied in Section 18 below.

## 17. Fixed-point filtering and literal phase (2 October 2026)

[Grover, *Fixed-point quantum search*, PRL **95**, 150501 (2005)](https://arxiv.org/abs/quant-ph/0503205),
Eq. (1) and Section 3, gives the three-call pi-over-three phase sequence
that cubes the failure probability. The
[radial-filter chapter](HOPF_RADIAL_FILTER.md) applies this established
sequence to the actual amplified Hopf source. It retains the complex
accepted scalar: cubic suppression of rejected amplitude leaves a
quadratic full operator error because an accepted phase remains.

The selective phase is not a free exact Clifford+T gate. Its
determinant-one representative factors into three one-qubit Pauli-Z
rotations on the existing two flags, and a literal Clifford global
correction calibrates the complete word. The standard phase-sensitive
one-qubit approximation theorem, [GKW Lemma 2.3](https://arxiv.org/html/2411.04790v3#S2),
supplies the six ancilla-free native phase words at their charged
precision. The source-width cap can then use a half-logarithmic
coefficient, while the added calls and phase words leave the established
asymptotic T-count and T-depth unchanged. The fine-word search distinction
discussed above still applies.

The [short-echo audit](HOPF_FLAG_ECHO.md) is separate: its exact error
formulas rule out uniform radial cancellation for four specific Pauli
interpositions, and include a full-space equal-mask exception. They are
not a no-go theorem for fixed-point amplification or general composites.

## 18. Grouped programs and chunked dirty queries (2 October 2026)

The [grouped-program theorem](GROUPED_PROGRAM_PREFETCH.md#8-complete-frame-theorem-at-fixed-accuracy)
combines the conditional source with exact prefetch into an active zero
suffix. The program stays read-only through source leakage, allowing its
actual inverse to erase it exactly. Private predicate trees and literal
enable-dependent amplification phases complete the local interface.
The [chunked indicator](CHUNKED_DIRTY_INDICATOR.md) interpolates between
the existing routed and read-only counter constructions. A dirty linear
tree conjugates a root flip into a selected-path translation, and the
outer echo removes every unknown dirty offset.

Conditional workspace, reversible conjunctions, linear conjugation,
bilinear lookup, and block-encoding amplification are established tools.
The local contribution is their charged complete-frame composition:
an adaptive early prefetch schedule and a distinct late chunk allocation
give fixed-accuracy $`O(n\log\log(n+2))`$ T-depth with optimal-order
$`O(\sqrt N)`$ T-count, $`O(N)`$ Clifford count, two external clean
flags, and sufficient $`C_\eta\sqrt N`$ dirty workspace. No generic
priority claim is inferred. The variable-precision theorem keeps its
previous bounds. This geometric-source schedule remains a valid
predecessor to the linear fixed-accuracy bound below; unrestricted
large-width depth optimality and the constant-clean high-precision
endpoint remain open.

## 19. Unary phase-source reuse and linear T-depth (3 October 2026)

The [unary phase-source compiler](UNARY_PHASE_GRADIENT.md) uses established
phase kickback. [Jones et al., arXiv:1204.0567, Section 2.1,
Eqs. (2)–(4), and Section 4.1, Fig. 15](https://arxiv.org/pdf/1204.0567)
describe Fourier eigenstates of modular shifts, programmable phases,
and repeated reference reuse. The local source uses a unary encoding,
coherently selects its cyclic shift from a loaded one-hot program, and
prepares/unprepares it inside the active zero suffix. Neither phase
kickback nor shared phase-reference preparation is claimed as new.

The cyclic convolution uses the standard three-product Karatsuba
identity; see [Iggy van Hoof, arXiv:1910.02849v2, Section 4.2](https://arxiv.org/html/1910.02849v2).
Recursive scalar products give the bilinear rank bound, and reduction
modulo $`X^q-1`$ is linear. The constant T-depth does not come from
van Hoof's space-efficient reversible multiplication schedule. It comes
from the local guarded trilinear-phase construction with private
conditional work and the retained native Toffoli word. Parallel
parity-phase synthesis with initialized ancillas is already explicit in
[Selinger, arXiv:1210.0974v2, Section 2, Eqs. (5)–(6), and
Theorem 4.1](https://arxiv.org/pdf/1210.0974). The local proof must
additionally establish identity on arbitrary inactive inputs and exact
temporary return on every active source state.

[Kim–Laakkonen, arXiv:2512.24982v1, Theorems 3, 5 and 6,
and Section 5.1](https://arxiv.org/html/2512.24982v1) already give
constant-depth controlled CNOT/Clifford circuits and a universal
logarithmic-size catalyst for rotations. Their catalytic rotation
chooses an angle-dependent CNOT matrix classically; its depth-one
implementation uses measurement-assisted uncomputation. Their charged
preparation uses measured phase estimation, expected repetitions, and
dynamic compilation after selecting the catalyst eigenvalue. This does
not directly supply the unary compiler's coherently loaded angle table,
unitary source boundary pair, or conditional-work return contract.
The present bound does not settle Kim's generic clean/dirty-only
constant-T-depth rotation question: source preparation is charged and
the complete frame has linear, rather than constant, T-depth.

A recent comparison is [Wu et al., *Shared Phase Arithmetic for Parallel
Quantum Rotations*, arXiv:2609.36574v1, Sections II.C, III.B and
IV.C–E](https://arxiv.org/html/2609.36574v1), submitted 29 September 2026
and checked here on 3 October. Their Theorem 1 expresses reversible
phase-function evaluation, one binary addition, and decoding; Proposition 1
gives Clifford-only encoding for disjoint binary supports. The displayed
ripple adder permits measurement/feedforward and has depth linear in
the phase-register width. Preparation is charged separately, and parallel
batches require separate resources. These are useful shared-arithmetic
precedents, not the constant-depth unary-program interface used here.

The additional result is the complete charged composition: one unitary
source boundary pair per group, coherent one-hot shift selection,
inactive-sector identity, the Hopf angle-stability bound, and an
early/late workspace and query allocation. For every fixed accuracy
$`\eta`$, it gives one complete real-frame circuit with

```math
D_T=O_\eta(n),\qquad T=O_\eta(\sqrt N),\qquad G=O_\eta(N),
```

two external clean flags, and sufficient $`C_\eta\sqrt N`$ dirty
workspace. The initialized-isometry error includes all returned work
and arbitrary dirty-reference entanglement. It uses no supplied phase
state or intermediate measurement. This improves the fixed-accuracy
depth upper bound; it proves neither an unrestricted matching depth
lower bound nor the constant-clean high-precision endpoint. The sources
above identify inherited ingredients and interface distinctions, not
priority for the composite construction. This square-root-width result
is also a corollary of the broader fixed-accuracy tradeoff below.

## 20. Blocked bilinear queries and the dirty-width tradeoff (3 October 2026)

The [blocked bilinear lookup](BLOCKED_BILINEAR_LOOKUP.md) combines the
existing [two-pass dirty traversal](AMORTIZED_DIRTY_LOOKUP.md#2-selecting-a-shear-with-dirty-unary-traversal)
with the [bilinear indicator echo](PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction).
Their lineage remains [LKS, Appendix C](https://arxiv.org/html/1812.00954v2),
[Khattar–Gidney, Sections 4 and 7](https://arxiv.org/html/2407.17966v2),
and the controlled-linear and commuting-basis constructions of
Kim–Laakkonen and Boyd discussed in Section 9. Standard Boolean phase
polarization supplies the controlled leaf: for a bilinear phase P,
apply $`(-1)^{dP}`$, toggle the dirty bit d by hz, apply the phase
again, and undo the toggle. This full four-step word leaves exactly
$`(-1)^{hzP}`$ and returns d.
Neither that algebra nor generic dirty cancellation is a novelty claim.

The local proof establishes a rank-sensitive native leaf with preserved
controls and exact helper return, a selected-block traversal, and a chunk
allocation whose depth sums across the final Hopf layers. Together with
the unary early groups, it gives, for each fixed $`0\lt\eta\le1/64`$,
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and
$`b\ge17(L+n+7)`$, one complete real-frame circuit with

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad G=O_\eta(N),
\qquad D_T=O_\eta\!\left(\frac N{b^2}+n\right).
```

Two external clean flags suffice. Full initialized-isometry error,
literal phase and arbitrary dirty-reference return retain their existing
contracts. Sources, actual inverses and all queries are charged.

When the interval is nonempty, the inherited count lower bound and the
physical width $`n+2+b=\Theta(b)`$ give simultaneous worst-case matches:

```math
17(L+n+7)\le b\le\sqrt{N/n},\qquad
T^\star=\Theta_\eta(N/b),\quad D_T^\star=\Theta_\eta(N/b^2).
```

This extends the fixed-accuracy matching window; it introduces no new
external synthesis premise or depth lower-bound method. The preceding
variable-accuracy bounds remain valid; Section 21 gives a uniform
precision refinement. The unary
square-root-width schedule remains a valid predecessor and corollary;
unrestricted large-width depth optimality and the constant-clean
high-precision endpoint remain open. No generic lookup priority is claimed.

## 21. Uniform precision and rectangular query allocation (3 October 2026)

The [uniform precision theorem](UNIFORM_PRECISION_DEPTH.md) uses the same
native selected-block query, dirty traversal and bilinear echo as Section
20. Rectangular block dimensions balance word precision against indicator
cost. The proof makes the unary-source cutoff uniform at low precision
and uses the existing hybrid when its source term absorbs the logarithmic
depth contribution. These are allocation and composition results; the
external ingredients in Sections 19–20 are unchanged, and no new native
query or synthesis premise is assumed.

With $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$, two clean flags
and $`b\ge17(L+n+7)`$, one complete real-frame circuit has absolute,
precision-independent constants in

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

For $`6\le L\le\log_2(n+2)/16`$, the sharper bounds are

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+n\right).
```

The existing count lower bounds divided by physical width give
simultaneous worst-case $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$ through $`b\le\sqrt{N/n}`$
generally, and through $`b\le\sqrt{NL/n}`$ in the low-precision
regime, above the literal threshold and when the intervals are nonempty.
The low-precision count is optimal at every eligible width; the general
count retains nL and is asserted optimal outside its matching interval
only under a sufficient condition such as $`L\le N/n^2`$.
All preparation, queries, actual inverses and work return remain charged.
The selected high-precision endpoint and unrestricted large-width depth
optimality remain open. No priority claim follows from this composition.

## 22. Scope of recent depth lower bounds (3 October 2026)

These comparisons guide further research; neither is an imported compiler
premise. [Parham, arXiv:2504.19966v1](https://arxiv.org/html/2504.19966v1),
Proposition 1.8, relates T-depth to alternations of unrestricted Clifford
and shallow circuits. Theorems 1.14–1.15 connect sufficiently strong
explicit-state and Boolean-function lower bounds to classical threshold
circuit lower bounds, with polynomial clean workspace in the model.
Section 6 explicitly says that no analogous reduction is known for
general prescribed-unitary implementation. This is therefore not a
blanket complexity barrier to complete-frame operator lower bounds.

[Al-Ghattas–Gamarnik–Kiani, arXiv:2610.02166v1](https://arxiv.org/html/2610.02166v1),
submitted 1 October 2026, treats arbitrary Clifford blocks. Corollary
1.7(iii) and Section 4.2 cover every fixed number of shallow blocks at
total width $`M=O(n)`$; arbitrary-width extensions concern specific
one-round classes. Lemma 4.4 retains an $`O(kM^2)`$ entropy term, so
this theorem does not cover the current $`b\asymp\sqrt N`$ allocation.
No applicable growing depth bound was identified in this comparison;
that audit outcome does not prove that such a bound is impossible.
