# Raw Hopf gradients from state interference

[← Original approximation contract](../../docs/QBP_APPROXIMATION.md) · [State-only compiler](STATE_ONLY_COMPILER.md) · [Frame-safe interface](../../docs/FRAME_SAFE_COMPILATION.md)

A different decoder estimates every real Hopf magnitude derivative using
state preparation, one controlled observable, and leaf measurements. It does
not apply an inverse frame. For branch probabilities bounded away from zero
and one, the original raw-coordinate sampling order is retained. A positive
reference state extends the construction to arbitrary real angles, including
singular chart endpoints, with at most an additional factor of $`n+1`$ in the
sufficient sample count.

The subsequent [coarse-frame decoder](COARSE_FRAME_QBP.md) removes that
angle-dependent sampling factor by adding a charged coarse inverse and
using corrected complex scores. The leaf-only construction below remains
a simpler option under its bounded-branch promise.

This is a different measurement protocol with the same ideal gradient means.
It does not preserve the original Walsh-record distribution or compile the
prescribed frame. The universal-observable necessity result in
[frame safety](../../docs/FRAME_SAFE_COMPILATION.md#necessity-for-all-observable-dependent-gradient-means)
concerns the original fixed decoder and remains unchanged. The discussion
here is restricted to the real magnitude chart; the existing
[complex phase stream](../../docs/QBP_APPROXIMATION.md#9-the-complete-complex-gradient)
is a separate result.

Logarithmic wave-function derivatives are standard in variational Monte Carlo;
see Toulouse and Umrigar, [Eq. (45)](https://arxiv.org/pdf/physics/0701039).
Ancilla interference is the standard Hadamard test, for example
Mitarai and Fujii, [Fig. 1](https://arxiv.org/pdf/1901.00015).
Those ingredients are inherited. The Hopf support bounds, reference-state
construction, and native resource substitution are the claims developed below;
no priority claim is made for the general estimator identity.

## 1. An exact universal-observable identity

Fix the angle tuple, let $`N=2^n`$, and write the real state and its raw
derivatives as

```math
\psi_x=\langle x|\psi\rangle,\qquad
d_{jx}=\partial_{\theta_j}\psi_x
       =a_j\langle x|e_j\rangle.
```

Choose a known normalized real reference $`r`$ with $`r_x\ne0`$ wherever
some $`d_{jx}\ne0`$. Prepare

```math
|\Omega\rangle
=\frac{|0\rangle_B|r\rangle+|1\rangle_B|\psi\rangle}{\sqrt2},
```

apply the phase-calibrated controlled Hermitian-unitary observable $`O`$ to
branch one, and measure $`B`$ in the X basis and the system in the leaf basis.
For X eigenvalue $`s=(-1)^b`$ and $`v=O\psi`$,

```math
\Pr(s,x)=\frac14|r_x+s v_x|^2.
```

The record

```math
Y_j(s,x)=2s\kappa_{jx},\qquad \kappa_{jx}=d_{jx}/r_x
```

therefore satisfies

```math
\mathbb E Y_j
=2\operatorname{Re}\sum_x d_{jx}(O\psi)_x
=2\operatorname{Re}\langle\partial_j\psi|O|\psi\rangle
=\partial_j\langle\psi|O|\psi\rangle.
```

This permits complex Hermitian-unitary $`O`$; only $`\psi,r,d_j`$ are real.
If all derivatives vanish at a leaf, a zero-reference version may define its
decoder to be zero there. The positive construction below avoids that case
entirely. No classical evaluation of $`(O\psi)_x`$ is supplied.

At each depth, a leaf belongs to exactly one node subtree. Consequently only
one coordinate of that depth's record can be nonzero. Any bound

```math
|\kappa_{jx}|\le\sqrt K
```

gives the deterministic depth-record bound $`\|Y^{(d)}\|_2\le2\sqrt K`$.
The final guarantee is raw-coordinate $`\ell_\infty`$ accuracy, obtained from
the depth blocks; this is not a dimension-free full-gradient $`\ell_2`$ claim.

## 2. Balanced branches need only the original state

Suppose every node obeys the explicit promise

```math
p\le\sin^2\theta_j,\cos^2\theta_j\le1-p,
\qquad 0\lt p\le\tfrac12.
```

Take $`r=\psi`$, including its real signs. Every leaf amplitude is nonzero.
Only one preparation of $`\psi`$ is required: prepare it, append the branch
in $`|+\rangle`$, and apply the controlled observable. For a leaf below node j,

```math
\frac{d_{jx}}{\psi_x}
=\begin{cases}
-\tan\theta_j,&x\text{ is in the left child subtree},\\
\cot\theta_j,&x\text{ is in the right child subtree},
\end{cases}
```

and the ratio is zero outside that subtree. Hence Section 1 applies with

```math
K=K_p:=\frac{1-p}{p}.
```

For fixed p, this has the same asymptotic sample order as the original
fixed-tuple protocol. It uses the standard logarithmic derivative directly;
there is no reference-state table or controlled two-state preparation.
The branch promise is sufficient, not an assumption about generic Hopf
parameters. It neither restricts the observable to a classically easy family
nor supplies classical access to its response amplitudes.

## 3. A reference for every real angle tuple

Define the derivative envelope and its total mass by

```math
w_x=\max_j d_{jx}^2,\qquad Z=\sum_x w_x.
```

The root derivative has norm one. Also
$`\sum_{j\text{ at depth }d}\|d_j\|_2^2=\sum_{j\text{ at depth }d}a_j^2=1`$.
Thus

```math
1\le Z\le\sum_j\|d_j\|_2^2=n.
```

Both extremes occur. When all angles are $`\pi/4`$, every nonzero derivative
entry has magnitude $`1/\sqrt N`$, so $`Z=1`$. When all angles are zero,
the n left-spine derivatives are distinct computational basis vectors and
all other derivatives vanish, so $`Z=n`$. Regular angle tuples can approach
the latter value. Under the branch promise of Section 2, one also has
$`w_x\le K_p\psi_x^2`$ and hence $`Z\le K_p`$.

Use the stable positive reference

```math
q_x=w_x+\frac1N,\qquad K=Z+1,\qquad
r_x=\sqrt{q_x/K}.
```

It is normalized and satisfies

```math
r_x\ge\frac1{\sqrt{N(n+1)}},\qquad
\frac{|d_{jx}|}{r_x}\le\sqrt K,\qquad 2\le K\le n+1.
```

The floor prevents tiny or zero denominators without requiring a decision
about whether an exact amplitude vanishes. It adds only one to the envelope
mass. The derivative identity remains exact at zero incoming amplitudes and
at angles zero or $`\pi/2`$; no inverse metric factor occurs.

### Linear-size classical construction

At each tree node maintain $`A`$, the squared incoming state amplitude, and
$`M`$, the largest squared derivative of that incoming amplitude with
respect to an ancestor angle. Initialize $`(A,M)=(1,0)`$ at the root. For
$`c=\cos\theta_j`$ and $`s=\sin\theta_j`$, propagate

```math
(A_L,M_L)=(Ac^2,\max\{Mc^2,As^2\}),\qquad
(A_R,M_R)=(As^2,\max\{Ms^2,Ac^2\}).
```

At a leaf, $`M=w_x`$. This division-free traversal, summation of Z, and
construction of q use $`O(N)`$ scalar operations and storage. Positive
subtree sums of q give a certified positive Hopf chart for r: each node's
squared cosine is its left mass divided by its total mass. Every child mass
is positive, so obtaining these reference angles requires no singular
zero-test convention. The original state is supplied by its existing Hopf
angle tuple.

## 4. Bias, finite coefficients, and sample count

Allow $`H=\sum_t c_tO_t`$ with known exact real coefficients and
$`\Lambda=\sum_t|c_t|\gt0`$. Choose t independently with probability
$`|c_t|/\Lambda`$ and multiply the record by
$`\Lambda\operatorname{sgn}(c_t)`$. This classical sampling and its
preprocessing have the same contract as
[the original reflection-sum protocol](../../docs/QBP_APPROXIMATION.md#10-reflection-sums-and-finite-classical-weights).
One reflection is the case $`\Lambda=1`$.

Require the preparation to approximate the complete joint state, including
all declared clean work and every arbitrary dirty input/reference, with error
$`\eta_P`$. In the two-state version this must preserve the relative branch
phase; two unrelated state vectors specified only up to phase do not suffice.
Let each actual controlled observable have full initialized-isometry error
$`\eta_O`$ on arbitrary branch-system inputs and reference correlations.
Keep all work until the final measurement. The complete output-state error
is at most $`\eta_P+\eta_O`$.

For every unit vector in one depth's coordinate space, the scalar decoded
observable has norm at most $`2\Lambda\sqrt K`$. The bounded-observable
inequality and Euclidean duality give depth-block bias at most

```math
4\Lambda\sqrt K(\eta_P+\eta_O).
```

If finite decoder coefficients obey
$`|\widetilde\kappa_{jx}-\kappa_{jx}|\le\tau`$, each empirical coordinate
changes by at most $`2\Lambda\tau`$, independently of the observed data.
Rounded coefficients are emitted only for nodes on the observed leaf's
ancestor path; structural zeros outside that path remain exactly zero.
It suffices, for $`0\lt\varepsilon_\infty\le\Lambda`$, to choose

```math
\eta_P,\eta_O\le\frac{\varepsilon_\infty}{32\Lambda\sqrt K},\qquad
\tau\le\frac{\varepsilon_\infty}{8\Lambda},\qquad
S\ge\frac{128\Lambda^2K\,[1+\ln(2n/\delta)]}
{\varepsilon_\infty^2}.
```

The bias, finite-coefficient, and sampling budgets are respectively at most
$`\varepsilon_\infty/4,\varepsilon_\infty/4,\varepsilon_\infty/2`$.
The resulting estimate obeys
$`\|\widehat{\nabla E_H}-\nabla E_H\|_\infty\le\varepsilon_\infty`$
with probability at least $`1-\delta`$.

For independent executions this follows from the depth-record norm bound.
It also holds with a reused dirty bank: conditional centered records are
Hilbert-space martingale differences of norm at most $`4\Lambda\sqrt K`$.
The [same martingale bound already used for the original protocol](../../docs/QBP_APPROXIMATION.md#101-reusing-a-dirty-bank-across-executions)
gives tail probability at most
$`2\exp[-S\varepsilon_{\rm stat}^2/(32\Lambda^2K)]`$ per depth.
Each conditional bias obeys the same uniform estimate. The system and all
declared clean inputs are initialized between completed executions; no
intermediate reset, measurement, postselection, or supplied catalyst is used.

For the floor construction, a certified upper bound on K can replace K in
all allocations. For Section 2 use the known $`K_p`$ directly. As in the
original protocol, these are fixed-tuple estimates of ideal chart
derivatives, not derivatives of a synthesized gate word.

## 5. Charged quantum preparation and complete-gradient costs

Put

```math
\ell=\max\{6,\lceil\log_2(32\Lambda\sqrt K/\varepsilon_\infty)\rceil\},
\qquad L'=\max\{n,\ell\}.
```

The [state-only compiler](STATE_ONLY_COMPILER.md#1-state-only-contract) prepares a known real Hopf
state to complete state-isometry error $`2^{-L'}`$ with two compiler clean
flags, $`b_F\ge L'+n+7`$ arbitrary dirty qubits, and

```math
t_P=O(N+L'),\qquad g_P=O(NL').
```

Its [coherent two-state corollary](STATE_ONLY_COMPILER.md#6-two-states-under-one-protocol-branch) prepares r or $`\psi`$ under the existing
unchanged branch, with the same asymptotic counts and the same stated dirty
threshold. Thus the positive reference is an explicit additional classical
table and a charged native preparation, not a state oracle. The coherent
corollary is what permits the shared branch and avoids silently allocating
another clean system register. Section 2 uses only the ordinary single-state
version.

Reserve the protocol's initialized branch and the observable's own work in
addition to the two compiler flags. Exactly returned scratch may be reused
only where its all-input contract permits. The observable need only meet
the error allocated by $`\ell`$; the extra $`L'\ge n`$ condition belongs to
this particular state compiler. All workspace inequalities use the actual
$`L'`$, including precision constants.

With $`\overline t_O=\sum_t|c_t|t_{O_t}/\Lambda`$, the expected total is

```math
\mathbb E\mathcal T_{\nabla,\rm ref}
\le S(t_P+\overline t_O)
=O\!\left(\frac{\Lambda^2K[1+\log(n/\delta)]}{\varepsilon_\infty^2}
 [N+L'+\overline t_O]\right).
```

Use the maximum oracle cost for a deterministic upper bound. The Clifford
count is $`O(S[NL'+\overline g_O+n])`$ in expectation. For fixed balanced
branch promise p, K is constant; the fine-precision native cost is linear in
$`N+L'`$ without an asymptotic sampling penalty relative to the original
bound. For the generic floor construction, K can grow linearly in n.

This does not establish general dominance. At fixed accuracy the original
protocol already has per-frame $`\Theta(\sqrt N)`$ T-count with zero
compiler clean qubits and $`\Theta(\sqrt N)`$ dirty work; that can be better
than this $`O(N+L')`$ preparation. The new construction uses a different
decoder, different records, and its own workspace reservation. It does not
resolve the complete-frame $`L=N`$ endpoint or prove optimal total gradient
cost. An end-to-end quantum advantage also requires accounting for the
controlled observable and the classical costs below; it is not inferred
from a favorable state-preparation T-count alone.

## 6. Classical preprocessing and output

The supplied angle descriptions must be read. At a fixed tuple the envelope
and reference chart require $`O(N)`$ scalar evaluations and arithmetic
operations, plus the state compiler's own table-generation cost. These are
not unit-cost claims about arbitrary-precision real numbers. The positive
floor bounds the conditioning of divisions and square roots by powers of
$`N(n+1)`$. Certified arithmetic with
$`P=O(n+L'+\log(n/\tau))`$ bits suffices for the reference and decoder
rounding budgets, with the costs of trigonometric evaluation, arithmetic,
square roots, and angle recovery at P bits charged. The balanced decoder
only needs its bounded tan/cot tables; it does not require division by a
computed exponentially small leaf amplitude.
These precision estimates concern the reference and decoder arithmetic.
They do not establish a bit-complexity theorem for the native compiler's
separate certified table generation; that preprocessing cost remains charged
under its own stated scope.

Record-wise decoding computes the n derivative entries on the sampled path
using prefix and suffix products and updates those n coordinates. Its
arithmetic cost is $`O(Sn+N)`$, including output materialization. There is
also a histogram route with $`O(S+N)`$ arithmetic operations. Let

```math
h_x=\frac{2\Lambda}{S r_x}
       \sum_{t:\,x_t=x}\operatorname{sgn}(c_{u_t})s_t.
```

The estimate is $`\widehat g_j=\sum_x h_xd_{jx}`$. Hold the measured h
fixed and differentiate the linear functional $`\sum_x h_x\psi_x`$ by one
classical reverse traversal of the tree. Explicitly, initialize $`\beta_x=h_x`$
at leaves and compute at node j

```math
\beta_j=c_j\beta_L+s_j\beta_R,\qquad
\widehat g_j=a_j(-s_j\beta_L+c_j\beta_R).
```

This contracts the known Jacobian
without constructing it as an $`N\times(N-1)`$ matrix. It does not
differentiate r, the sampled data, or a compiled quantum circuit. In the
balanced version, equivalent left/right subtree histogram sums multiply
the node's bounded tan/cot coefficients.

Histogram counters need $`O(1+\log S)`$ bits. Finite multiplication,
division, and reverse accumulation must meet the same deterministic
decoder-error budget; the arithmetic-operation count does not waive those
bit costs. Reading the tuple and writing all $`N-1`$ coordinates already
cost $`\Omega(N)`$ classical work. Preprocessing can be reused over shots
at one fixed tuple; changing the tuple requires a new applicable preparation
and new classical weights.

## 7. Finite verification

[Reference-state tests](../../tests/test_reference_state_qbp.py) compare the
envelope traversal with explicit raw derivatives and check the complete
interference probabilities for complex Hermitian-unitary observables. They
include generic, negative, balanced, and endpoint angles, the extremal Z
values, and the deterministic finite-coefficient bound. The same small file
checks the state compiler's half-amplitude and amplification identities.
These are finite algebraic checks, not a native implementation or an
end-to-end advantage experiment.
