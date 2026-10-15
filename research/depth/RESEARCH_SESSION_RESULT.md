# Fixed-accuracy depth research session: execution record

[Session design](RESEARCH_SESSION_PROMPT.md) · [Preparation map](DEPTH_OBSTACLE_MAP.md) · [R78 proof](SUBTREE_PROJECTOR_DECODER.md) · [Current resource statement](../../docs/OPEN_PROBLEM.md)

This records the original fixed-accuracy session. Its later
[precision-uniform follow-up](UNIFORM_PRECISION_SESSION_RESULT.md)
proved R79; the original decisions and evidence below are retained.

**Started:** 2026-10-15 10:44:44 UTC. **Execution baseline:** main
`aab766b`. The authorized budget was at most five hours; it was a ceiling,
not a requirement to continue after a qualified result. The design's earlier
baseline `321101d` is retained in the linked prompt.

The primary lane produced **local result R78**, a fixed-accuracy
polylogarithmic-depth compiler. The reserve lanes supplied two scoped
counterexamples and a cache reduction, preserved below. No growing depth
lower bound or external novelty claim was established. This record does
not replace the separate theorem and proof.

## 1. Result and exact scope

For fixed $`0\lt\eta\le1/64`$, let
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and $`N=2^n`$.
R78 gives one coherent circuit with one external initialized qubit and
arbitrary dirty width b, at the original literal thresholds:

| Family | Required width |
|---|---|
| Complete real Hopf frame | $`b\ge17(L+n+7)`$ |
| Complex magnitude frame $`D_\phi W_{\mathbb R}`$ | $`b\ge17(L+n+8)`$ |

The same circuit satisfies

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\mathrm{op}}\le\eta,
```

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(N/b^2+\log_2^2(n+2)\right).
```

This includes every logical column, literal phase, signed and singular
angles, initialized-work leakage and return, and arbitrary dirty inputs
and references. G counts elementary Clifford gates. Every coherent query,
source preparation, inverse, guard, and reflection is charged.

At square-root dirty width the depth is $`O_\eta(\log^2(n+2))`$,
meeting the session's stretch target with optimal worst-case T-count.
The inherited depth lower bound matches the new upper bound throughout
each nonempty interval from the applicable literal threshold to
$`\sqrt N/\log_2(n+2)`$. Optimal depth beyond that interval remains open.
The constants depend on fixed accuracy: R78 is not a precision-uniform
replacement for R72. Neither a scalable native emitter nor efficient
classical preprocessing is supplied.

## 2. Primary lane: change the recovery mechanism

The earlier compressed-column construction paid for
$`\Theta(\sqrt n)`$ outer recovery calls. The successful replacement
uses the complete frame's nested subtree scaling projectors
$`A_0\le\cdots\le A_h=I`$. Their differences are exactly the
orthogonal state/marker-depth bands. A coherent binary search identifies
the band in $`O(\log h)`$ projector queries. Computing its depth label,
XORing it into the old source tag, and reversing the search erases that
tag while retaining each literal frame column.

The projectors are compiled from the original descendant loader. Its
source-isometry guarantee transports an anchor projector with a global
operator-norm error bound; no implementation of the target frame or its
inverse is supplied to the decoder. All fresh source work, including the
locally prepared gradient, belongs to the anchor test. Controls occur on
the central phase, allowing the actual exterior source words to cancel
exactly on inactive occupied inputs.

The protected logical suffix now needs only $`M=O(\log n)`$ bits.
Its exponential reduction in early-table size absorbs the polynomial
query overhead; the direct one-clean bottom construction handles its
remaining layers. An explicit width split and the existing native
fallback cover every eligible allocation and finite dimension. The
one-clean diagonal supplies complex phase dressing by a full-isometry
hybrid, without reset or an assumption of exact intermediate return.

Independent reviews checked full-mode reflection error, literal search
signs, occupied inactive work, source multiplicities, elementary Clifford
costs, width reuse, the all-width splice, and inherited lower-bound
quantifiers. The details and error/resource ledger are in the
[R78 proof](SUBTREE_PROJECTOR_DECODER.md).

## 3. Reserve lane B: transporting a dirty cache retains the indicator problem

The alternative proposal retained a one-hot XOR cache across source
Hadamards to avoid repeatedly computing nonlinear address indicators.
It produced an exact reduction, not an improved query construction.

Let $`E_s|x,Y\rangle=|x,Y\oplus e_x\rangle`$ be the phase-free
indicator on an s-bit address and an arbitrary $`2^s`$-bit dirty bank.
For address coordinate i define $`C_i=E_sH_iE_s`$. Since $`E_s`$
commutes with $`Z_i`$ and is its own inverse,

```math
C_iZ_iC_i=E_sX_iE_s=X_iD_i.
```

Here $`D_i`$ flips both bank bits $`Y_{(a,0)},Y_{(a,1)}`$, where
a is x without coordinate i. Their unordered pair is independent of
$`x_i`$. Let $`\mathcal L`$ apply
$`\mathrm{CNOT}(Y_{(a,0)};Y_{(a,1)})`$ to every pair. It maps the
XOR shift (1,1) to (1,0), so

```math
\mathcal L X_i C_i Z_i C_i\mathcal L^\dagger
```

is exactly an $`(s-1)`$-bit one-hot XOR indicator into the first wire
of each pair, identity on the other bank and the address. This is a
full-space identity, including arbitrary dirty/reference inputs.

A full-input transported-Hadamard implementation of depth d, T-count t,
and Clifford count g would therefore give the shorter indicator with
depth at most $`2d`$, count at most $`2t`$, and Clifford count
$`2g+O(2^s)`$, with the same returned helpers. For a full-operator
approximation, using the actual inverse for the second occurrence gives
at most twice the implementation error by unitary telescoping.

Thus the natural cache update supplies almost the original indicator as
a subroutine. It does not establish a universal indicator lower bound.
A source-specific promised encoded action could evade this full-input
reduction, but would need a new complete recovery/return proof. None was
obtained. Lane B did not qualify as the primary positive mechanism.

## 4. Reserve lane C: two lower-bound falsifiers

No dimension-growing invariant qualified within the early scout. The
worker returned to independent review of the constructive lane. Two
derived counterexamples rule out specific tempting measures of hardness.

### Nonlinear multiplexors with growing derivative rank and constant depth

For $`n=2r+1`$, put $`f(a,b)=\bigoplus_{i=1}^r a_ib_i`$ and
$`D_f|a,b,t\rangle=(-1)^{tf(a,b)}|a,b,t\rangle`$. Every address
variable participates; the Boolean derivatives span dimension $`2r`$
and the quadratic polar form has rank $`2r`$. Nevertheless,

```math
4abt=a+b+t-(a\oplus b)-(a\oplus t)-(b\oplus t)
       +(a\oplus b\oplus t)\pmod8
```

packs $`D_f`$ into three physical T layers. The parity groups are
$`\sum_i(a_i+b_i)+rt`$,
$`-\sum_i[(a_i\oplus b_i)+(b_i\oplus t)]`$, and
$`\sum_i[-(a_i\oplus t)+(a_i\oplus b_i\oplus t)]`$.
The latter two use the invertible Clifford coordinate maps
$`(a_i,b_i,t)\mapsto(a_i\oplus b_i,b_i\oplus t,t)`$ and
$`(a_i,b_i,t)\mapsto(a_i\oplus t,a_i\oplus b_i\oplus t,t)`$.
The first uses $`T_t^r=S_t^{\lfloor r/2\rfloor}T_t^{r\bmod2}`$.
Thus $`D_f`$ needs no ancilla, at most three T layers,
$`6r+(r\bmod2)`$ T gates, and $`8r+O(1)`$ Clifford gates.

Set every upper Hopf angle to zero and the deepest angles to
$`\theta_{f(x)}`$. The complete frame is the multiplexor
$`U_f=\sum_x|x\rangle\langle x|\otimes R_y(\theta_{f(x)})`$.
For $`\theta_0=0,\theta_1=\pi/2`$,
$`U_f=H_tD_fH_tD_f`$, literally, since $`XZ=R_y(\pi/2)`$.
Its exact T-depth is at most six. More generally, with
$`c=(\theta_0+\theta_1)/2`$ and $`d=(\theta_0-\theta_1)/2`$,

```math
U_f=R_y(c)_tD_fR_y(d)_tD_f.
```

The [existing one-qubit one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md#6-complete-real-frames-and-the-grouped-extension)
handles these two rotations
with $`O_\eta(1)`$ depth and dirty work at fixed error. The resulting
complete-input circuit has $`D_T=O_\eta(1)`$ and $`T,G=O(n)+O_\eta(1)`$.
It fits the session budgets with unused wires untouched. Active-variable
count and derivative/polar rank alone therefore cannot force growing
depth. This says nothing about arbitrary truth tables or an AND witness.

### The balanced nested flag is already realized by a Clifford

For the frame with every angle $`\pi/4`$, let

```math
A_d=I_{2^d}\otimes|+\rangle\langle+|^{\otimes(n-d)},\qquad
B_d=I_{2^d}\otimes|0\rangle\langle0|^{\otimes(n-d)}.
```

The complete balanced frame and a free Clifford transport the entire
cumulative input-marker flag identically:

```math
W_{\rm bal}B_dW_{\rm bal}^\dagger
=H^{\otimes n}B_dH^{\otimes n}=A_d\quad\text{for every }d.
```

Consequently a lower-bound invariant depending only on these transported
projectors cannot distinguish this target from a Clifford. This does not
make $`W_{\rm bal}`$ Clifford: orientation inside each marker band is
additional information. Its ideal spectral label eraser also satisfies
$`E_{\rm spec}=(H^{\otimes n}\otimes I)E_{\rm marker}(H^{\otimes n}\otimes I)`$.
This is a balanced-case check, not a free generic decoder.

## 5. Validation and remaining limits

Five new tests in
[`test_subtree_projector_decoder.py`](../../tests/test_subtree_projector_decoder.py)
passed. They check finite algebra, phases, work leakage, and arbitrary
unitary completions; their completions are verification fixtures, not
native source implementations or asymptotic proofs. The scout separately
checked 2,728 exact mod-eight assignments, small nonsingular multiplexors,
and balanced projector transport through five qubits. The all-dimension
claims above follow from the displayed algebra, not those finite checks.

Local validation completed with **922 tests passing**, including the five
new checks. The 40 documentation checks passed; all 133 documents passed
the rendered-presentation audit. The four fault-tolerant receipt suites
matched their unchanged fixtures. Compilation, reviewer walkthrough,
coarse real/complex and residual-QBP examples, offline provenance, and
both required resource ledgers also passed. The integrating pull request
records the exact reviewed tree and remote CI gates.

No external literature theorem is imported by this execution record.
At the end of this session, precision-uniform improvements beyond R72
remained open; the later R79 follow-up resolves that upper-bound step.
The unrestricted growing-depth lower bound at square-root width,
practical crossover, and a general native emitter remain unresolved.
