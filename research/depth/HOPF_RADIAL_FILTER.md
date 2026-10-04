# A phase-calibrated radial filter for complete Hopf frames

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Error geometry](../../docs/HOPF_ERROR_ACCUMULATION.md) · [Short flag echoes](HOPF_FLAG_ECHO.md) · [Operator source](../../docs/OPERATOR_SOURCE_COMPILER.md) · [Current frontier](../../docs/OPEN_PROBLEM.md)

A three-call fixed-point filter changes the physical source circuit so
that its error relative to the encoded polar rotation is quadratic in
the radial coefficient defect. Its selective phases can be synthesized
on the same two flags, with their precision and literal global phase
charged. The resulting complete frame combines ideal-angle square-sum
stability with an additive quadratic remainder. It admits a smaller
source-precision cap, while retaining the current asymptotic count and
depth bounds.

The fixed-point phase sequence is inherited from
[Grover, *Physical Review Letters* **95**, 150501 (2005)](https://arxiv.org/abs/quant-ph/0503205),
Eq. (1) and Section 3. Its cubic suppression of failure probability is
not a new amplification method. The statements below supply the literal
phase-sensitive operator error, two-flag native phase implementation,
and complete-Hopf-frame precision allocation needed here.

## 1. The source error splits into angle and radius

On each addressed row, use the existing actual circuit

```math
K=XZ=-iY,\qquad B=J^\dagger QJ=(cI+sK)/2,
\qquad A=-QRQ^\dagger RQ,\qquad R=I-2JJ^\dagger.
```

J initializes only the two compiler flags. Every other wire may be
arbitrary and entangled with a reference. For
$`r=\sqrt{c^2+s^2}\gt0`$, put

```math
U=(cI+sK)/r=R_y(\phi),\qquad
\delta=1-r^2,\qquad a=\frac{r(3-r^2)}2.
```

The accepted amplified block is exactly $`J^\dagger AJ=aU`$.
In the regime $`|\delta|\le1`$, take $`r\gt0`$ so that
$`0\lt a\le1`$, and define the rejected probability

```math
t=1-a^2=\frac{\delta^2(3+\delta)}4.
```

Thus $`0\le t\le\delta^2`$. U contains the angular programming
error; t measures the remaining departure from that polar unitary.
The formulas apply separately on each coherent row. They do not discard
the rejected flag sectors between calls.

## 2. Exact fixed-point filtering, including its phase

Set $`z=e^{i\pi/3}`$ and

```math
R_z=I+(z-1)JJ^\dagger,\qquad
F=z^{-2}A R_z A^\dagger R_z A.
```

The inverse is the actual inverse of the same A. Since
$`(z-1)^2=-z`$ and $`z-1=z^2`$, expansion gives

```math
FJ=z^{-1}t\,AJ+a\,JU.
```

In particular,

```math
J^\dagger FJ=a(1+z^{-1}t)U,\qquad
\|(I-JJ^\dagger)FJ\|=t^{3/2}.
```

The accepted scalar is generally complex. Its phase must not be dropped:
the leakage amplitude is cubic in the original leakage amplitude, but
the full phase-sensitive distance is quadratic. Unitarity yields the
exact expression

```math
\|FJ-JU\|^2
=2-2\sqrt{1-t}(1+t/2)
=\frac{t^2(3+t)}{2[1+\sqrt{1-t}(1+t/2)]}
\le2t^2\le2\delta^4.
```

These are full initialized-isometry identities, including all dirty
inputs and reference extensions. Across a table, take the maximum over
the rows. The residual phase is included in this bound, with no
additional angle-dependent correction.

For comparison, consider the same three-call word with selective phases
$`\alpha=e^{i\varphi}`$ and $`\beta=e^{i\psi}`$. Its rejected
multiplier before global calibration is

```math
f=\alpha+(\alpha-1)(\beta-1)a^2.
```

First-order rejected cancellation as $`a\to1`$ requires
$`1+\beta(\alpha-1)=0`$. Unit modulus forces
$`\alpha=\beta=e^{i\pi/3}`$ or their common conjugate. This is a
statement about this three-call selective-phase form, not all composite
circuits.

## 3. Native selective phases use the same two flags

An exact pi-over-three selective phase is not treated as a native gate.
Instead use its determinant-one representative on flags a,b:

```math
\widetilde R=e^{-i\pi/12}R_z
=e^{i\pi Z_a/12}e^{i\pi Z_b/12}e^{i\pi Z_aZ_b/12}.
```

The equality follows from
$`JJ^\dagger=(I+Z_a)(I+Z_b)/4`$ on the flag register.
There is a literal Clifford global correction:

```math
F=(-i)A\widetilde R A^\dagger\widetilde R A.
```

The scalar $`-iI`$ is Clifford: the literal Pauli word $`YXZ=-iI`$
acts on any existing wire. Omitting this scalar would change
the complete-frame phase by $`\pi/2`$ even at zero radial defect.

Each factor in $`\widetilde R`$ is a determinant-one one-qubit Z
rotation, with the ZZ factor conjugated by two CNOTs. Approximate each
of the six one-qubit rotation occurrences in the complete word within
$`\tau/6`$ in literal operator norm. The phase-calibrated one-qubit
synthesis primitive, [GKW Lemma 2.3](https://arxiv.org/html/2411.04790v3#S2),
already used in the [fault-tolerant compiler](../../docs/FAULT_TOLERANT_COMPILER.md),
then costs $`O(\log(1/\tau))`$ native gates per occurrence and
uses no ancilla. Telescoping the six replacements gives a native circuit
$`\widetilde F`$ with

```math
\|\widetilde F-F\|\le\tau,\qquad
\|\widetilde FJ-JU\|\le\sqrt2\,\delta^2+\tau.
```

This is a word-length and circuit-resource statement. It does not add
an unconditional efficient fine-word search algorithm. The existing
distinction between native word existence and classical synthesis time
remains in force.

The filter contains two forward A calls and one actual inverse, hence
nine forward/inverse Q calls in total. Its fixed number of extra
reflections, scalar Cliffords, and CNOTs is charged. With
$`\tau=\Delta^2`$ and $`\Delta=2^{2-m}`$, the six phase
approximants cost $`O(m)`$ T gates, Clifford gates, and T-depth.
They touch only the existing flags. No clean work, measurement, reset,
or supplied phase state is added.

The lookup selectors and suffix work continue to return exactly after
each actual Q or inverse. The source core and flags return within the
full approximation bound. Approximate selective phases need not leave
inactive rows exactly unchanged; their error is included in tau.

## 4. A physical square-sum-plus-quadratic frame bound

For a target angle theta, certified coefficient programming gives

```math
\zeta=\|(c,s)-(\cos\theta,\sin\theta)\|_2
\le\frac{5\sqrt2}{8}\Delta\lt\Delta.
```

Assume $`\Delta\le1/4`$. Then
$`|\delta|\le2\zeta+\zeta^2\le3\Delta`$ and
$`|\delta|\lt1`$. Choosing the polar angle's real lift nearest theta
gives $`|\phi-\theta|\le2\Delta`$. One way to see the angle
bound is that a radius-Delta ball around a unit vector subtends angle
at most $`\arcsin\Delta\le2\Delta`$.

With $`\tau=\Delta^2`$, each filtered source layer therefore obeys

```math
\|\widetilde F_dJ-JL_d(\phi)\|
\le(9\sqrt2+1)\Delta_d^2\lt15\Delta_d^2.
```

Tensor identity on all unused workspace. Telescope these errors against
the exact polar-angle frame, then apply the
[sharp angle theorem](../../docs/HOPF_ERROR_ACCUMULATION.md#1-a-sharp-angle-error-bound-on-the-full-frame)
to that frame. For $`S=\sum_d2^{-2m_d}`$, this proves

```math
\|\widetilde VJ-JW(\theta)\|
\le2\sqrt{2\sum_d\Delta_d^2}+15\sum_d\Delta_d^2
=8\sqrt2\sqrt S+240S.
```

This bound holds for the modified physical circuit, including leakage
and arbitrary dirty/reference inputs. It does not assign a square-sum
bound to the unmodified stages ruled out by the earlier counterfamily.

## 5. A smaller source-precision cap

Let $`L\ge6`$, $`\eta=2^{-L}`$, and n be positive. Define

```math
h=\left\lceil\frac12\log_2(8n)\right\rceil,\qquad
m_d=L+4+\min\{n-d,h\},\qquad d=0,\ldots,n-1.
```

Equivalently, h is the least integer with $`4^h\ge8n`$.
Every assigned m is at least $`L+5`$, and is at most its value under
the earlier additive cap. Consequently all original layerwise dirty
reservations and weighted query upper bounds remain available.

The new squared-precision sum satisfies

```math
S\le4^{-(L+4)}\left(\sum_{k=1}^{h-1}4^{-k}+n4^{-h}\right)
\le\frac{11}{6144}\eta^2.
```

The physical full-frame error is therefore at most

```math
\sqrt{\frac{11}{48}}\,\eta+\frac{55}{128}\eta^2\lt\eta.
```

Writing $`\bar h=\min(n,h)`$, the assigned source precision totals

```math
\sum_dm_d=n(L+4)+n\bar h-\frac{\bar h(\bar h-1)}2
=nL+\frac12n\log_2n+O(n).
```

The logarithmic coefficient is halved for the source widths. This is
not a factor-two gate saving: the filter triples the amplified-stage
calls and adds the separately priced phase words. Those words have
precision of order $`2m_d`$. The total native cost still has order
$`nL+n\log(n+2)`$ for this part of the construction.

## 6. Resource frontier and evidence boundary

Apply the same old-query and dirty-counter hybrid schedules to the
smaller source tables, reserve the same base workspace, and multiply
their source/query calls by the fixed filter overhead. The extra phase
words add only $`O(\sum_dm_d)`$ to T count, Clifford count, and
serial T-depth. This fits the already retained orders. Thus, with two
clean flags and $`b\ge17(L+n+7)`$, the modified complete real-frame
circuit retains

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL+n\log_2(n+2)\right).
```

The same previously proved matching workspace interval remains valid.
The high-precision complete-frame endpoint and unrestricted large-width
depth remain open. Repeating a fixed number of radial filters would
still leave angular programming and the charged source depth to address.

The [bounded filter checks](../../tests/test_radial_filter.py) verify the
exact source/filter algebra, literal phase calibration, actual inverses,
phase-approximation error accounting, and the finite precision budget.
The phase fixtures use exact mathematical rotations and bounded
perturbations; they do not emit the general native approximating words.
Their word-length and complete-frame resource guarantees are analytic.
The separate [echo checks](../../tests/test_flag_echo.py) preserve the failed
two-call candidates and their exact equal-mask exception.
