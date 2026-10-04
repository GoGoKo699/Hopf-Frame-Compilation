# Proof roadmap: exact and fault-tolerant Hopf-frame compilation

[Claims at a glance](README.md) · [Full statements](manuscript/PUBLICATION_SCOPE.md) · [Proof index](docs/README.md) · [Evidence](docs/VERIFICATION.md)

The four principal results price the same prescribed operator in two gate
models. Result A reaches the exact all-workspace size/depth frontier.
Result B shares precision with sufficient clean work. Result C reduces the
clean reservation to one qubit. Result D schedules count and T-depth on one
complete real-frame circuit. This roadmap explains how the proofs connect;
the linked chapters contain the constructions and estimates.

The [bounded claim-to-proof audit](docs/CORE_CLAIM_AUDIT.md) found no
unresolved blocker in this selected package. It is an internal mathematical
review, not formal certification or independent peer review. Final manuscript
writing remains on hold.

## 1. Why the complete operator is the target

For $`N=2^n`$, the real Hopf frame specifies

```math
W|0^n\rangle=|\psi\rangle,\qquad
W|\lambda(j)\rangle=|e_j\rangle,\qquad
\partial_{\theta_j}|\psi\rangle=a_j|e_j\rangle.
```

The last identity uses the oriented amplitude entering a tree node. On the
canonical chart it is the nonnegative metric square root. At zero weight the
derivative vanishes, while the supplied parameter tuple still specifies the
marker column. The [Hopf interface](docs/HOPF_INTERFACE.md) fixes the angle,
marker, orientation, and singular-chart conventions consumed by every proof.

With $`R_y(\theta)=e^{-i\theta Y}`$, the compiler sees addressed layers

```math
W=L_{n-1}\cdots L_0,\qquad
L_d=I+\sum_{p=0}^{2^d-1}|p\rangle\langle p|
\otimes(R_y(\theta_{d,p})-I)
\otimes|0^{n-d-1}\rangle\langle0^{n-d-1}|.
```

These are identities on the full system space. In the inverse-frame gradient
protocol, a different completion can preserve the state and change the
answer: the [two-qubit counterexample](docs/COMPILER_BOUNDARIES.md#2-two-qubit-global-state-column-counterexample)
changes the decoded gradient from $`(2,0,0)`$ to $`(0,\sqrt2,0)`$.
The [necessity theorem](docs/FRAME_SAFE_COMPILATION.md#necessity-for-all-observable-dependent-gradient-means)
assumes exact preparation of the same state. At regular points, preserving
all observable-dependent means with the same compiler and its adjoint then
forces the full frame up to common phase. This does not assert necessity
for every gradient algorithm or for singular directions.

The exact compiler therefore satisfies $`VJ_m=J_mW`$, where $`J_m`$
appends clean zeros. Its actual inverse has the same clean-input guarantee.
In the approximate model, with $`a`$ clean and $`b`$ dirty qubits,

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

The norm covers all logical inputs, arbitrary dirty states and external
references, literal phases, and leakage outside the clean subspace.
The [actual-adjoint and composition proofs](docs/QBP_APPROXIMATION.md#2-the-actual-adjoint-has-the-same-error)
carry earlier leakage through later unitaries; they never assume that
approximate work has become exactly clean between factors.

The phase-dressed complex magnitude target is $`D_{\rm ph}W`$, with
separately supplied certified phases. Its phase derivatives use a separate
record. The weaker checkpoint active-interface contract and its distinct
mean-preservation guarantee are proved in [compiler boundaries](docs/COMPILER_BOUNDARIES.md).

## 2. Result A: exact schedules cover every clean budget

The upper bounds hold for every parameter tuple. The lower bounds hold in the
worst case over the Hopf-frame family, uniformly in the clean-workspace budget.

The exact model has arbitrary one-qubit gates, CNOTs, and all-to-all
connectivity. For every integer $`m\ge0`$, the target size and total depth are

```math
S=\Theta(N),\qquad D=\Theta\!\left(n+\frac{N}{n+m}\right).
```

The [exact theorem](docs/COMPILER_THEOREM.md) proves three schedules and
splices them without a workspace gap. Their common task is to implement
a prefix-selected rotation only when the complete lower suffix is zero.

| Clean budget | Mechanism | Authoritative argument |
|---|---|---|
| $`m=0`$ | Borrow one original suffix data qubit; a toggle-and-rotation echo removes dependence on its unknown input and restores it | [Strict-zero echo](docs/COMPILER_THEOREM.md#4-strict-zero-workspace) |
| $`1\le m\lt4n`$ | Compute the suffix predicate into one reusable clean flag, apply a uniformly controlled gate, then uncompute | [Direct flagged schedule](docs/COMPILER_THEOREM.md#5-small-positive-workspace) |
| $`m\ge4n`$ | Cut the tree, compile its conditioned prefix, route the remaining data coherently, and run subframes in parallel | [Tree cut](docs/COMPILER_THEOREM.md#6-exact-tree-cut) and [routed tail](docs/COMPILER_THEOREM.md#8-routed-parallel-tail) |

The zero-work proof checks all borrowed-bit and suffix-predicate sectors,
not only the preparation input. For small positive work, a direct schedule
has an additive $`n^2`$ depth term, which is absorbed in the claimed
bound throughout that budget range.

For larger work, the tree cut is a full-operator factorization into a
suffix-conditioned prefix and a direct sum of tail frames. The prefix uses
an explicit [binary–one-hot decoder](docs/COMPILER_THEOREM.md#7-conditioned-prefix-construction).
The intervening Givens word preserves the entire one-hot code space, so
its actual decoder inverse returns the work for arbitrary logical inputs.

The router preserves the prefix while moving the suffix to its selected
branch. It erases coherent control copies before reusing their wires as
flags, then reconstructs them for the actual reverse routing. Inactive
branches are identity. A cut with $`B=2^t`$ branches and suffix size
$`s=n-t`$ fits the conservative work envelope $`2B(s+1)`$.
The maximal feasible cut, including the saturated $`s=1`$ case, yields
the depth bound even when much supplied workspace is unused.

The [lower bounds](docs/COMPILER_THEOREM.md#9-matching-lower-bounds)
use family dimension for size and width-dependent depth, and backward
light cones for the additive $`\Omega(n)`$. CNOT count needs a separate
argument when arbitrary one-qubit gates are free. With $`K`$ CNOTs,
at most $`2K`$ ancillary wires participate; inactive clean wires supply
only a removable common phase. Fusing local gates leaves at most
$`n+4K`$ relevant one-qubit slots. This gives worst-case CNOT count
$`\Theta(N)`$ for $`n\ge2`$; the one-qubit case has no CNOTs.

A [literal phase diagonal](docs/COMPILER_THEOREM.md#10-phase-dressed-complex-magnitude-frame)
extends the result to the complex magnitude frame using the same pool
sequentially. Yuan–Zhang's UCG, multi-controlled-X, and coherent copy tools
are imported; the operator factorization and workspace schedules adapt them
to the prescribed completion. The [source map](docs/SOURCE_MAP.md) separates
these dependencies from local contributions.

## 3. Result B: pay fine precision once with sufficient clean work

For $`0\lt\eta\le1/64`$, define
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$,
$`h=1+\lceil\log_2(L+n+2)\rceil`$, and $`q=n+a+b`$.
The [fault-tolerant theorem](docs/FAULT_TOLERANT_COMPILER.md#1-target-resources-and-theorem)
proves, for sufficiently large absolute $`C`$ and $`a\ge C(n+h)`$,

```math
T^\star=\Theta(\sqrt{NL}+L+NL/q),\qquad G=O(NL).
```

The clean reservation holds the source, labels, arithmetic, and failure
record while leaving a constant fraction of available width for lookup.
An exact dirty-bank query cancels the borrowed input on every basis sector;
linearity gives joint return with arbitrary references. Dirty wires are
never treated as initialized program storage.

At low precision a direct addressed construction fits this bound. At high
precision, [grouped residual dictionaries](docs/FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame)
separate short coarse circuits from fine corrections. Their sizes $`Q_i`$
satisfy $`\sum_iQ_i=O(N)`$ and $`\sum_i\sqrt{Q_i}=O(\sqrt N)`$.
Coarse programs have length $`O(n+h)`$; fine coefficients are queried
one bit at a time using a separate digit address.

One [geometric source and local kernel](docs/FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel)
supply the weighted digits throughout the correction stream. A retained
failure counter prevents rejected branches from re-entering the final
accepted block. This uses the inherited block-product compression gadget,
with every history wire and update charged. Accepted blocks are a proof
device: the physical circuit performs no measurement or postselection.

The [full-output argument](docs/FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization)
controls accumulated source error, reverses the actual preparation once,
and performs final normalization-two amplification. Its reflection includes
all potentially disturbed initialized work. The source costs $`O(L)`$
once; summed lookup costs supply the other T terms. The construction returns
dirty banks exactly, while clean approximation work obeys the global norm.

The [matching lower bounds](docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
reduce a real-frame subfamily to diagonal synthesis and use fixed-width
coherent counting, retaining the published synthesis/counting premises.
The [complex extension](docs/FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame)
has its own literal-diagonal proof. Neither this clean reservation nor
optimal T-depth follows as a necessity from Result B.

## 4. Result C: operator precision and conditional suffix work

The [one-clean construction](docs/ONE_CLEAN_COMPILER.md) encodes binary
weights in a Pauli operator on arbitrary dirty work. Paired Majorana
sources and sign queries produce scalar overlaps proportional to the
identity on that work. A conjugated scalar word combines two overlaps into
a rotation; literal five-call amplification uses one clean signal flag.
The source, programs, actual inverses, and rejected-output error are charged.

Layerwise use gives an $`O(N+nL)`$ baseline. The
[conditional-suffix argument](docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit)
reduces repeated precision charges. A group is active only on its declared
zero-suffix sector. There that suffix supplies conditional clean workspace;
on the complementary sector the completed subroutine is identity.
Residual column forests become rank-one stars, and coarse programs stream
one symbol at a time. Their small private work permits rapidly growing
groups, with only $`O(\ell_*(n))`$ precision charges.

The complete-input error includes final suffix-scratch and predicate return.
It does not promote conditional zeros to globally initialized ancillas.
For $`\ell_*(n)=1+\log_2^*(n+2)`$, the resulting real-frame bound is

```math
a=1,\quad b\ge L+n+7:\qquad T=O(N+L\ell_*(n)),\quad G=O(NL).
```

At $`b\ge2(L+n+7)`$, banked lookup gives
$`T=O(\sqrt{NL}+L\ell_*(n)+NL/b)`$. This is an upper bound, with
matching subregimes specified in the [publication scope](manuscript/PUBLICATION_SCOPE.md).
A literal diagonal composes with the grouped frame using error $`\eta/2`$
per factor; the separately proved complex magnitude thresholds become
$`L+n+8`$ and $`2(L+n+8)`$.

The [diagonal and U(2) multiplexor corollaries](docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
are complete-operator capabilities with their own reservations. The
[zero-clean corollary](docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
is only the layerwise real-frame construction. Its borrowed signal does
not make the grouped or literal-phase constructions zero-clean.

## 5. Result D: count and depth on the same circuit

The [uniform-precision theorem](docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary)
uses two external clean flags and $`b\ge17(L+n+7)`$ dirty qubits:

```math
T=O(\sqrt{NL}+NL/b+nL),\quad D_T=O(NL/b^2+nL),\quad G=O(NL).
```

Its late groups use blocked bilinear queries. Unequal indicator widths
preserve square-root precision dependence in count while retaining the
query depth bound. The proof sums the capped-precision tails and fits the
simultaneous helper peaks inside the literal reservation.

For early groups, a charged unary-source schedule works through
$`6\le L\le\log_2(n+2)/16`$. The explicit cutoff is uniformly
sublinear in $`n`$, giving

```math
T=O(\sqrt{NL}+NL/b),\qquad D_T=O(NL/b^2+n),\qquad G=O(NL).
```

At higher precision the existing hybrid's additional logarithmic term is
absorbed by $`nL`$. A finite-size fallback keeps all constants absolute.
This is an asymptotic sufficient range, not a practical crossover estimate.

The inherited count lower bound divided by physical width gives matching
count and depth on $`17(L+n+7)\le b\le\sqrt{N/n}`$, extended to
$`17(L+n+7)\le b\le\sqrt{NL/n}`$ in the low-precision range,
whenever the relevant interval is nonempty. On those intervals,
$`T^\star=\Theta(NL/b)`$ and $`D_T^\star=\Theta(NL/b^2)`$.
The sufficient condition $`L\le N/n^2`$ also absorbs $`nL`$ in
count at every eligible width; it does not establish depth optimality there.

This is a real-frame theorem, with no automatic complex extension.
T-depth permits arbitrary Clifford interlayers; their gate count remains
charged, but this is not total circuit depth. Earlier schedules remain
useful at some smaller reservations or higher precisions. Later component
refinements do not change these complete-frame orders.

## 6. Operational consequence, open gaps, and evidence

[Exact substitution](docs/QBP_CONSEQUENCE.md#4-frame-safe-substitution)
preserves the global inverse-frame measurement record. For simultaneous
absolute accuracy of raw Hopf-coordinate gradients, the magnitude stream
uses $`O(\log n)`$ executions at fixed accuracy and confidence. The
comparison assumes phase-calibrated controlled access to the same specified
Hermitian-unitary observable in scalar and gradient programs.

The [approximation bridge](docs/QBP_APPROXIMATION.md) bounds fixed-parameter
bias using the actual compiled circuit and its actual adjoint. It also
covers the separate phase stream, reflection sums, rounded classical
weights, and correlated dirty-bank reuse with fresh declared clean inputs.
Quantum executions, per-execution resources, and classical output remain
separate costs. No derivative of discontinuously synthesized words or
optimality among all gradient algorithms is inferred.

The completed [state-based QBP supplement](supplements/state_based_qbp/README.md)
changes the decoder and avoids the fine complete frame. It is a distinct
task-level result, not a closure of either compiler gap:

| Unresolved question | Current guarantee |
|---|---|
| $`a=2,b=N+n+7,L=N,n\ge3`$ | $`\Omega(N)\le T^\star\le O(N\ell_*(n))`$ |
| Fixed accuracy, $`a=2`$, sufficient $`b=\Theta_\eta(\sqrt N)`$ | $`T=\Theta_\eta(\sqrt N)`$ with $`D_T=O_\eta(n)`$; unrestricted depth lower bound $`\Omega(1)`$ |

The [research index](research/README.md) preserves the literal open targets,
promised-update results, and scoped attempts. Interface obstructions do not
supply unrestricted lower bounds. Neither gap blocks the selected package.

Use the [verification map](docs/VERIFICATION.md) for proof-to-code coverage,
exact receipts, and reproduction commands. Finite checks test indexing,
phases, inverses, work return, and ledgers; they do not prove asymptotic
claims. Native finite components are distinguished from analytic
constructions, and a general elementary shared-source emitter remains
absent. The [source map](docs/SOURCE_MAP.md) records inherited primitives
and reductions. These boundaries are part of the claims, not additional
research requirements before final writing.
