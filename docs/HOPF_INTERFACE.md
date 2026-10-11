# Hopf interface used by the compiler

[← Landing page](../README.md) · [Complete narrative](../REVIEW.md) · [Formal compiler theorem →](COMPILER_THEOREM.md)

The compiler theorem needs only a compact interface from the Hopf construction.
The inverse coordinate map, optimization experiments, and the complete QBP
protocol remain in the two earlier Hopf repositories.

## 1. Four facts consumed by the synthesis proof

Fix $n$ system qubits and write

```math
N=2^n.
```

The real Hopf construction supplies a unitary $W_{\mathbb R}(\boldsymbol\theta)$ with the following
properties.

### State column

```math
\boxed{
W_{\mathbb R}|0^n\rangle
=|\psi_{\mathbb R}(\boldsymbol\theta)\rangle.
}
```

### Marker columns

For every internal tree node $j$, a nonzero computational marker $\lambda(j)$ is
assigned to a unit frame direction $\lvert e_j\rangle$:

```math
\boxed{
W_{\mathbb R}|\lambda(j)\rangle
=|e_j(\boldsymbol\theta)\rangle.
}
```

### Coordinate differential

Let $a_j(\boldsymbol\theta)$ be the oriented amplitude entering node $j$.  Then

```math
\boxed{
\partial_{\theta_j}|\psi\rangle
=a_j|e_j\rangle,
\qquad
g_{j,j}=a_j^2.
}
```

### Addressed-layer product

The frame is a product of $n$ complete-operator layers:

```math
\boxed{
W_{\mathbb R}^{(n)}
=L_{n-1}^{(n)}\cdots L_1^{(n)}L_0^{(n)},
}
```

with $L_0$ acting first.

These four statements define the synthesis target.  In particular, the target
is not specified only by the state column.

## 2. Tree and marker convention

The $N-1$ magnitude coordinates are indexed by the internal nodes of a complete
binary tree.  Nodes are numbered breadth first, beginning at one.  If node $j$
has depth $d$ and position $r$, then

```math
j=2^d+r,
\qquad
0\leq r\lt 2^d.
```

Its marker is

```math
\boxed{
\lambda(j)=(2r+1)2^{n-d-1}.
}
```

The corresponding bit string consists of:

1. the $d$-bit prefix locating the node;
2. a one at the node's target position;
3. zeros in every lower position.

The state occupies marker zero.  The remaining markers place the frame
directions in a known computational basis, which is what makes inverse-frame
readout possible.

### Notation

| Symbol | Meaning |
|---|---|
| $n$ | number of system qubits |
| $N=2^n$ | Hilbert-space dimension |
| $m$ | clean ancillary qubits supplied to the compiler |
| $j$ | breadth-first internal-node index |
| $d$ | depth of node $j$ |
| $r$ or $p$ | position or prefix at depth $d$ |
| $\theta_{d,p}$ | magnitude angle at that node |
| $a_j$ | oriented incoming amplitude |
| $g_{j,j}=a_j^2$ | diagonal metric weight |
| $\lvert e_j\rangle$ | unit marker-frame direction |
| $\lambda(j)$ | computational marker assigned to $\lvert e_j\rangle$ |
| $W_{\mathbb R}$ | real Hopf differential frame |
| $D_{\mathrm{ph}}$ | diagonal leaf-phase layer |
| $W_{\mathbb C,\mathrm{mag}}=D_{\mathrm{ph}}W_{\mathbb R}$ | phase-dressed complex magnitude frame |

Basis labels are ordered from the most significant tree decision to the least
significant one.  The implementation uses the corresponding nonnegative integer
labels.

## 3. Incoming amplitude, metric, and chart boundary

At each internal node, the local split is

```math
|v_j\rangle
=\cos\theta_j|v_{2j}\rangle
+\sin\theta_j|v_{2j+1}\rangle.
```

Differentiating the local split produces the unit complement

```math
|e_j\rangle
=-\sin\theta_j|v_{2j}\rangle
+\cos\theta_j|v_{2j+1}\rangle.
```

The derivative of the complete state carries the amplitude that has already
reached node $j$.  This amplitude is the product of the ancestor sine and cosine
factors selected by the path to $j$:

```math
\partial_{\theta_j}|\psi\rangle
=a_j|e_j\rangle.
```

For unrestricted real angles, $a_j$ may be negative.  The metric entry and its
principal square root are

```math
g_{j,j}=a_j^2,
\qquad
\sqrt{g_{j,j}}=|a_j|.
```

### Canonical domains

For the real chart:

- depths $0,\ldots,n-2$ use
  ```math
  [0,\pi/2];
  ```
- the final magnitude depth uses
  ```math
  [0,2\pi),
  ```
  which carries the leaf signs.

For the complex chart, every magnitude angle uses

```math
[0,\pi/2],
```

and the leaf phases carry the complex arguments.

Only ancestor angles enter an internal node's incoming amplitude.  The final
real depth therefore does not affect the sign of any $a_j$.  On both canonical
magnitude domains,

```math
a_j\geq0,
\qquad
a_j=\sqrt{g_{j,j}}.
```

### Singular coordinates

If

```math
g_{j,j}=0,
```

then

```math
\partial_{\theta_j}|\psi\rangle=0.
```

There is no derivative-normalized tangent at that parameter value.  The unit
vector in marker column $\lambda(j)$ remains the chart-selected orthogonal
continuation determined by the complete parameter tuple.  The compiler target
is therefore still a well-defined unitary, while the corresponding raw
coordinate derivative has zero weight.

The public implementation keeps these notions separate:

- `incoming_amplitude`: the oriented value valid for unrestricted angles;
- `metric`: `incoming_amplitude**2`;
- `sqrt_metric`: the principal nonnegative square root;
- `regular_coordinate_mask(atol=...)`: a tolerance-aware numerical regularity
  classification;
- `in_canonical_magnitude_domain`: the declared chart-domain check.

## 4. Complete addressed depth operator

At depth $d$, split a computational-basis label as

```math
|p\rangle_P|x\rangle_T|z\rangle_Z,
```

where:

- $p$ is the $d$-bit prefix;
- $x$ is the next qubit and is the rotation target;
- $z$ is the lower suffix of length
  ```math
  s=n-d-1.
  ```

The complete depth operator is

```math
\boxed{
L_d^{(n)}
=I+
\sum_{p=0}^{2^d-1}
|p\rangle\!\langle p|
\otimes
\bigl(R_y(\theta_{d,p})-I\bigr)
\otimes
|0^s\rangle\!\langle0^s|.
}
```

The rotation convention is

```math
R_y(\alpha)=e^{-i\alpha Y}
=
\begin{pmatrix}
\cos\alpha&-\sin\alpha\\
\sin\alpha&\cos\alpha
\end{pmatrix}.
```

Thus the prefix chooses the angle, the next qubit is the target, and the
complete lower suffix supplies one shared all-zero predicate.  Every nonzero
suffix sector is fixed.

The zero-suffix projector is essential.  It prevents a deeper tree rotation
from changing frame columns established by earlier depths.  It is also the
structural feature exploited by all three compiler schedules.

> **Proof checkpoint.** $L_d$ is a full operator, not a rule restricted to the
> preparation path.  Every proposed compiler must preserve its identity action
> on all nonzero-suffix sectors.

### 4.1 Complete-frame parameter geometry

The state metric $`g_{j,j}=a_j^2`$ can vanish, but the complete operator
retains every angle as an independent local coordinate. With the
unnormalized Frobenius inner product, its parameter metric is exactly

```math
\boxed{
\mathrm{Tr}\!\left[
\left(\partial_{\theta_j}W_{\mathbb R}\right)^T
\partial_{\theta_k}W_{\mathbb R}
\right]=2\delta_{j,k}.
}
```

This identity holds for every real angle tuple, including every singular
state chart. Moreover, $`\boldsymbol\theta\mapsto W_{\mathbb R}`$ is
injective modulo the individual $`2\pi`$ angle periods.

**Proof of the metric identity.** Order the individual rotations breadth
first, as in the marker convention. For node $`j=2^d+r`$, write

```math
\mu_j=(2r)2^{n-d-1},\qquad b_j=\lambda(j),\qquad
K_j=|b_j\rangle\langle\mu_j|-|\mu_j\rangle\langle b_j|,
```

```math
R_j=e^{\theta_jK_j},\qquad
V_j=R_j\cdots R_1,\qquad V_0=I.
```

The marker $`b_j`$ has not been touched by any earlier rotation. The other
endpoint $`\mu_j`$ is zero or an earlier marker. Thus $`V_{j-1}`$
preserves the span of zero and the earlier markers, fixes $`|b_j\rangle`$,
and sends $`|\mu_j\rangle`$ under its transpose to a unit vector
$`|u_j\rangle=V_{j-1}^T|\mu_j\rangle`$ in that earlier span. Cancelling
the later rotations in the derivative gives

```math
A_j:=W_{\mathbb R}^T\partial_{\theta_j}W_{\mathbb R}
=V_{j-1}^TK_jV_{j-1}
=|b_j\rangle\langle u_j|-|u_j\rangle\langle b_j|.
```

Every nonzero entry of $`A_j`$ has one index equal to the fresh marker
$`b_j`$. For $`k\lt j`$, the support of $`A_k`$ lies entirely in the earlier
span. Hence $`\mathrm{Tr}(A_j^TA_k)=0`$ for distinct indices,
while $`\|A_j\|_F^2=2`$. Left multiplication by the orthogonal matrix
$`W_{\mathbb R}`$ preserves these inner products, proving the formula.

**Proof of injectivity.** A final-depth marker column is

```math
W_{\mathbb R}|b_j\rangle
=-\sin\theta_j|\mu_j\rangle+\cos\theta_j|b_j\rangle
\qquad(d=n-1).
```

Earlier layers fix this input basis vector, and the final-depth rotations
act on disjoint pairs. These columns therefore determine all final-depth
sines and cosines, hence their angles modulo $`2\pi`$. Remove that
determined layer by left multiplication by its transpose. The same
argument recovers the preceding layer, and induction recovers the whole
tuple. This is a uniqueness proof; the compiler need not extract angles
or decide whether an incoming amplitude vanishes. ∎

### 4.2 Complete-frame operator stability

The Frobenius metric above is distinct from the operator norm needed by
the compiler. The nested marker supports give two complementary bounds
on that norm, valid at every real angle tuple, including singular state
charts. Write $`W=W_{\mathbb R}`$ and consider a differentiable angle path.
Define the old and new subspaces and the preceding layer product by

```math
\mathcal S_d=\mathrm{span}\{|p0^{n-d}\rangle:|p|=d\},\qquad
\mathcal M_d=\mathcal S_{d+1}\ominus\mathcal S_d,\qquad
V_d=L_{d-1}\cdots L_0,\quad V_0=I.
```

Earlier layers preserve $`\mathcal S_d`$ and fix its orthogonal
complement pointwise. Put

```math
E_d=\sum_p\dot\theta_{d,p}
 |p1\,0^{n-d-1}\rangle\langle p0\,0^{n-d-1}|,
\qquad R_d=E_dV_d,\qquad R=\sum_dR_d.
```

Then $`R_d:\mathcal S_d\to\mathcal M_d`$, and the complete
logarithmic derivative and its positive semidefinite bound are

```math
\boxed{
W^T\dot W=R-R^T,\qquad
Q=R^TR=\sum_dR_d^TR_d,\qquad
\|\dot W\|\le2\sqrt{\|Q\|}.
}
```

In particular, $`Q`$ retains the locations and magnitudes of the
individual angle errors, rather than only a maximum at each depth:

```math
Q=\sum_d V_d^T
 \left(\sum_p\dot\theta_{d,p}^{\,2}
 |p0\,0^{n-d-1}\rangle\langle p0\,0^{n-d-1}|\right)V_d.
```

**Proof.** Disjoint rotations in one layer give
$`L_d^T\dot L_d=E_d-E_d^T`$. Since $`V_d`$ fixes
$`\mathcal M_d`$ pointwise, $`V_d^TE_dV_d=E_dV_d=R_d`$.
Cancelling the later layers in the product derivative therefore gives
$`W^T\dot W=\sum_d(R_d-R_d^T)`$. The spaces
$`\mathcal M_d`$ are mutually orthogonal, so
$`R_d^TR_e=0`$ for distinct depths and $`R^TR=Q`$.
Finally, $`\|R-R^T\|\le2\|R\|=2\sqrt{\|Q\|}`$. ∎

**The factor two is sharp.** At the zero angle tuple, set the derivative
to one only at prefix $`p=1^d`$ at each depth, and to zero elsewhere.
These edges form a path through the successive vertices
$`0^n,10^{n-1},110^{n-2},\ldots,1^n`$. On that path,
$`R`$ is a forward shift. Consequently

```math
Q^2=Q,\qquad \|Q\|=1,\qquad
\|R-R^T\|=2\cos\!\left(\frac{\pi}{n+2}\right)\longrightarrow2.
```

No smaller universal factor can replace two in this PSD bound.

#### Depth caps and finite perturbations

The complementary bound uses only the largest angle speed at each depth:

```math
e_d=\max_p|\dot\theta_{d,p}|,\qquad
\Lambda_0=0,\qquad
\Lambda_{d+1}=\frac{\Lambda_d+\sqrt{\Lambda_d^2+4e_d^2}}2.
```

Then $`\|\dot W\|\le\Lambda_n`$. If every derivative at depth
$`d`$ has magnitude $`e_d`$, equality holds at every base angle
tuple, with arbitrary signs at individual nodes and with zero caps
allowed. The [canonical stability proof](HOPF_ERROR_ACCUMULATION.md#nested-support-proves-the-recurrence)
derives the block recurrence

```math
K_{d+1}=\begin{pmatrix}K_d&-B_d^T\\B_d&0\end{pmatrix},
\qquad \|B_d\|=e_d.
```

For equal magnitudes within a depth, $`B_d^TB_d=e_d^2I`$.
Each eigenvalue $`\mu`$ of $`iK_d`$ therefore gives two eigenvalues

```math
\frac{\mu\pm\sqrt{\mu^2+4e_d^2}}2
```

of $`iK_{d+1}`$. This explains the exact norm recurrence and its
independence from base angles and derivative signs; the linked chapter
also proves the stronger finite relative-spectrum statement. Neither
the PSD bound nor the depth-cap bound dominates the other in every
instance: the selected-path example has PSD bound two, while uniform
depth magnitudes attain $`\Lambda_n`$.

For a straight path $`\theta(t)=\theta+t\Delta\theta`$, use the
constant caps $`e_d=\max_p|\Delta\theta_{d,p}|`$. Integration gives

```math
\|W(\theta+\Delta\theta)-W(\theta)\|
\le\min\left\{2,\Lambda_n,
 2\int_0^1\sqrt{\|Q(t)\|}\,dt\right\},
\qquad
\Lambda_n\le\sqrt{e_0^2+2\sum_{d=1}^{n-1}e_d^2}.
```

These are literal-phase operator-norm bounds; tensoring identity on
untouched dirty work or a reference preserves them. They do not certify
the full initialized output of an approximate source circuit. That
distinction and the sharp square-root dependence on depth are treated
in [ideal-angle error accumulation](HOPF_ERROR_ACCUMULATION.md).

**An isolated angle still needs absolute precision.** If only one angle
changes by $`\delta`$, all other factors cancel under unitary
invariance of the norm:

```math
\|W(\theta+\delta\mathbf e_j)-W(\theta)\|
=\|R_y(\theta_j+\delta)-R_y(\theta_j)\|
=2|\sin(\delta/2)|.
```

Here $`\mathbf e_j`$ denotes a coordinate unit vector, not a marker
column. The identity is independent of the incoming state amplitude,
even when that amplitude is zero. Thus the PSD bound improves the
accumulation of actual angular errors; it does not make a small branch
insensitive, regularize a relative-error source gauge, or supply a
precision-sharing compiler. The endpoint resource bounds are unchanged.

## 5. Two-qubit orientation

For two qubits, the three magnitude angles give

```math
|\psi\rangle
=
\cos\theta_1
\bigl(\cos\theta_2|00\rangle+
      \sin\theta_2|01\rangle\bigr)
+
\sin\theta_1
\bigl(\cos\theta_3|10\rangle+
      \sin\theta_3|11\rangle\bigr).
```

The unit frame directions are

```math
|e_1\rangle
=-\sin\theta_1
 \bigl(\cos\theta_2|00\rangle+
       \sin\theta_2|01\rangle\bigr)
+\cos\theta_1
 \bigl(\cos\theta_3|10\rangle+
       \sin\theta_3|11\rangle\bigr),
```

```math
|e_2\rangle
=-\sin\theta_2|00\rangle+
 \cos\theta_2|01\rangle,
```

```math
|e_3\rangle
=-\sin\theta_3|10\rangle+
 \cos\theta_3|11\rangle.
```

The differential weights are

```math
\partial_{\theta_1}|\psi\rangle=|e_1\rangle,
```

```math
\partial_{\theta_2}|\psi\rangle
=\cos\theta_1|e_2\rangle,
\qquad
\partial_{\theta_3}|\psi\rangle
=\sin\theta_1|e_3\rangle.
```

The marker assignment is

```math
\lambda(1)=10_2,
\qquad
\lambda(2)=01_2,
\qquad
\lambda(3)=11_2.
```

Hence, in computational column order,

```math
W_{\mathbb R}
=
\begin{pmatrix}
|&|&|&|\\
\psi&e_2&e_1&e_3\\
|&|&|&|
\end{pmatrix}.
```

This column order is the entire reason that a state-equivalent completion need
not be a valid differential-frame completion.  The numerical obstruction is
given in [the complete counterexample](COMPILER_BOUNDARIES.md#2-two-qubit-global-state-column-counterexample).

## 6. Phase-dressed complex magnitude frame

Attach one phase $\phi_\ell$ to every leaf and define

```math
D_{\mathrm{ph}}
=\sum_{\ell=0}^{N-1}
e^{i\phi_\ell}|\ell\rangle\!\langle\ell|.
```

The unitary used by the complex magnitude stream is

```math
\boxed{
W_{\mathbb C,\mathrm{mag}}
=D_{\mathrm{ph}}W_{\mathbb R}.
}
```

It contains the complex state and the $N-1$ phase-dressed magnitude-frame
directions.  The $N$ leaf-phase differentials are localized in the
computational basis and use a separate signed one-hot measurement record.  They
are not additional columns of this $N$-dimensional unitary.

Writing a basis label as $x=zb$, with the final bit as target,

```math
D_{\mathrm{ph}}
=
\sum_z |z\rangle\!\langle z|
\otimes
\begin{pmatrix}
e^{i\phi_{z0}}&0\\
0&e^{i\phi_{z1}}
\end{pmatrix}.
```

The phase dressing is therefore one exact $n$-qubit UCG.  Once the real frame
is compiled, the complex magnitude theorem is a short sequential-composition
corollary.

## 7. What the compiler theorem does not need

The following parts of the broader Hopf program are not premises of the
compiler construction:

- the inverse coordinate map;
- optimization convergence;
- the four-qubit pedagogical example;
- application-specific observables;
- a hardware-native gate set;
- the complete checkpoint protocol.

They remain important for the broader research program, but the synthesis
question is fully specified by the state column, marker columns, differential
weights, chart domains, and addressed-layer product above.

### Executable counterpart

- [Tree geometry and frame matrices](../compiler_robust_hopf/frames.py)
- [Marker convention](../compiler_robust_hopf/conventions.py)
- [Complex magnitude and phase derivatives](../compiler_robust_hopf/complex_analysis.py)
- [Frame and chart-domain tests](../tests/test_frames.py)
- [Full-operator differential and stability tests](../tests/test_frame_operator_stability.py)
- [Complex analysis tests](../tests/test_complex_analysis.py)

---

[← Landing page](../README.md) · [Complete narrative](../REVIEW.md) · [Formal compiler theorem →](COMPILER_THEOREM.md)
