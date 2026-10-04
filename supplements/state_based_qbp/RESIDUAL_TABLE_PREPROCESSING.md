# Certified preprocessing of the residual rotation table

This note replaces the finite Euler search in the state-only residual
table by direct algebraic programming of three rotations. It concerns
classical table construction after the actual coarse native word has
been supplied. It does not supply a polynomial-time algorithm for
arbitrary single-qubit native-word synthesis.

## 1. Input and output contract

Use the conventions

```math
U(z)=\begin{pmatrix}z&-\sqrt{1-|z|^2}\\
\sqrt{1-|z|^2}&\overline z\end{pmatrix},\qquad |z|\le1,
```

and

```math
R_y(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix},\qquad
R_z(\theta)=\operatorname{diag}(e^{-i\theta},e^{i\theta}).
```

In the residual construction, $`z`$ is either the root amplitude or
$`\sqrt N`$ times a nonroot amplitude. Its unit-disk promise follows
from the supplied coarse word's distance guarantee.

Let $`q\ge5`$ and supply dyadic $`z_0=a_0+ib_0`$ with the certificate

```math
|z_0-z|\le\delta,\qquad \delta=2^{-2q-20}.
```

The output consists of certified cosine/sine intervals for a chronological
$`R_z,R_y,R_z`$ word. It approximates $`U(z)`$ within $`2^{-q-4}`$
before native synthesis. No angle, inverse trigonometric function, or
test of an unknown exact zero is required. The helper
[residual_table_preprocessing.py](../../compiler_robust_hopf/residual_table_preprocessing.py)
implements precisely this interval interface; it emits no quantum gates.
Its consistency check on $`|z_0|`$ cannot replace the caller's error
certificate.

## 2. Shrink and a small-radius branch

Put $`\zeta=(1-2\delta)z_0=a+ib`$. Then

```math
|\zeta|\le(1-2\delta)(1+\delta)\lt1,
\qquad |\zeta-z|\le4\delta.
```

For any unit-disk entries at distance $`e`$, the completion obeys

```math
\|U(z)-U(\zeta)\|
\le e+\sqrt{2e}.
```

Indeed the diagonal difference has norm $`e`$, and the off-diagonal
difference has norm at most the square root of
$`\bigl||z|^2-|\zeta|^2\bigr|\le2e`$. Thus the shrink costs less
than $`2^{-q-8}`$: use $`\sqrt8\lt3`$ and
$`4\delta\le2^{-q-10}`$.

Set $`\tau=2^{-q-6}`$. If the finite rational comparison
$`a^2+b^2\le\tau^2`$ holds, emit $`U(0)=R_y(\pi/2)`$ with identity
outer rotations. For $`r=|\zeta|\le\tau`$,

```math
\|U(\zeta)-U(0)\|
\le r+1-\sqrt{1-r^2}\le2\tau.
```

The total classical completion error is therefore less than
$`(1/256+1/32)2^{-q}\lt2^{-q-4}`$. The comparison is on supplied
dyadic data, so a tiny nonzero original coefficient needs no special
oracle or separation promise.

## 3. One half-phase used twice

Otherwise $`r=\sqrt{a^2+b^2}\gt\tau`$. Set
$`s=\sqrt{1-r^2}`$ and choose $`w=u+iv`$ as follows:

```math
\begin{array}{ll}
a\ge0:&u=\sqrt{\frac12+\frac{a}{2r}},\quad v=\frac{b}{2ru};\\[2mm]
a\lt0:&v=\operatorname{sgn}_+(b)
 \sqrt{\frac12+\frac{|a|}{2r}},\quad u=\frac{b}{2rv}.
\end{array}
```

Here $`\operatorname{sgn}_+(0)=1`$ refers only to the finite dyadic
$`b`$. The primary component has magnitude at least $`1/\sqrt2`$;
the other division consequently has a known denominator bound.
The identities $`|w|=1`$ and $`w^2=\zeta/r`$ follow directly.
With $`D(w)=\operatorname{diag}(w,\overline w)`$,

```math
D(w)\begin{pmatrix}r&-s\\s&r\end{pmatrix}D(w)
=\begin{pmatrix}rw^2&-s\\s&r\overline w^2\end{pmatrix}
=U(\zeta).
```

Thus the required cosine/sine pairs are $`(u,-v),(r,s),(u,-v)`$.
The same half-phase must be used twice. Replacing both copies of $`w`$
by $`-w`$ changes neither the product nor its literal phase; replacing
only one would introduce a branch-relative minus sign. This also explains
why a discontinuous half-phase convention across the negative real axis
does not create a discontinuity in the completed target.

Feed these certified coefficients directly into the
[paired-source digit rule](../../docs/ONE_CLEAN_COMPILER.md#3-conjugating-scalar-blocks-produces-a-rotation).
The $`R_z`$ factors are fixed Clifford conjugates of the real primitive.
Its proof depends on the sine/cosine pair, not on a stored numerical angle.
Intervals of width $`2^{-(q+20)}`$, with the same precision for fixed
algebraic constants, more than suffice for the rule's required mean
error $`2^{-q-1}`$ and head-sign estimate error $`1/16`$.
The existing $`43\,2^{-q}`$ full-operator native error per rotation is
unchanged. Consequently

```math
\|\widehat U-U(z)\|\lt(129+1/16)2^{-q}\lt130\,2^{-q}.
```

Direct sums take the maximum row error. The actual-inverse three-call
state amplification therefore retains its bound
$`390\,2^{-q}`$. At $`q=P+10`$ this is less than $`2^{-P}`$,
including arbitrary dirty inputs and reference correlations, with the
same quantum word and workspace as the retained state construction.

## 4. Fixed-precision rational certificates

For a nonnegative rational $`x`$, let

```math
j=\left\lfloor\sqrt{\left\lfloor2^{2B}x\right\rfloor}\right\rfloor.
```

Integer square root encloses $`\sqrt x`$ between $`j2^{-B}`$ and
$`(j+1)2^{-B}`$; if the lower endpoint squares to $`x`$, both endpoints
may be that point. These are exact finite rational certificates.

For requested interval width $`2^{-t}`$, use
$`B=t+4q+40`$ and $`h=2^{-B}`$. This is a fixed precision choice,
not a refinement loop depending on the coefficient's distance to a
singularity. Enclose $`r`$ and $`s`$ separately from the exact dyadic
$`a^2+b^2`$. Intersect their intervals with $`[0,1]`$.
In the nonsmall branch the lower endpoint for $`r`$ is at least
$`\tau/2`$. Interval arithmetic gives the following conservative bounds:

| Quantity | Guaranteed range or interval width |
|---|---|
| $`r`$ | Width at most $`h`$; lower endpoint at least $`\tau/2`$ |
| $`1/2+\lvert a\rvert/(2r)`$ | Intersect with $`[1/2,1]`$; width at most $`h/\tau^2`$ |
| Primary half-phase magnitude | Intersect with $`[1/2,1]`$; width at most $`3h/\tau^2`$ |
| Denominator $`2r`$ times the primary component | Absolute value at least $`\tau/2`$; width at most $`8h/\tau^2`$ |
| Secondary half-phase component | Width at most $`32h/\tau^4`$ |

Clipping the final components to $`[-1,1]`$ preserves enclosure. Since
$`32h/\tau^4=2^{-t-11}`$, every returned interval has the required
width. All divisions have certified denominators. Near the unit circle,
$`s`$ is enclosed by direct square root rather than by a condition-number
estimate for inverse trigonometry. Square roots and finite rational
comparisons suffice even at the original endpoints $`z=0,\pm1,\pm i`$.

## 5. Preprocessing cost after the coarse word is supplied

Here is a conservative computational version for the bounded-input model.
The real Hopf angles and supplied real leaf-phase representatives lie in
$`[-8,8]`$, with at most $`B_{\rm in}`$ input bits each. The gauge is
the arithmetic mean of those representatives, as in the
[complex coarse compiler](COMPLEX_COARSE_COMPILER.md#1-target-gauge-and-workspace).
There is no phase unwrapping or free reduction of large angle descriptions.
The target state and the coarse word use that same literal gauge.

Supply the actual coarse word $`C`$, its exact dirty-return contract, and
its coarse distance guarantee. Let $`J`$ be the total length of its
recorded local native words, $`\ell`$ their maximum length, and $`M`$
the number of two-mode updates used to apply those blocks to an amplitude
vector. Repeated phase-prefix rows count every suffix in $`M`$.
For the retained complex coarse construction, $`M=O(Nn)`$.
Generating or verifying the supplied coarse native words is a separate
computational task.

Evaluate the Hopf state, the bounded phase exponentials, and the recorded
native blocks; then apply $`C^\dagger`$ and form the root and scaled-tail
entries. One sufficient fixed-point accuracy for these vector operations is

```math
R=2q+\lceil n/2\rceil+\lceil\log_2(N+M+1)\rceil+O(1).
```

Use another $`\lceil\log_2(\ell+1)\rceil+O(1)`$ guard bits while
evaluating an individual block. A rounded two-mode unitary update costs
$`O(2^{-R})`$ vector error at bounded intermediate norm; summing the
$`M`$ updates and initial state construction costs
$`O((N+M)2^{-R})`$. Multiplication of tail entries by $`\sqrt N`$
then still leaves the coefficient error at most $`\delta`$, after fixing
the absolute guard constant. The same reasoning covers the coherent
common-reference table, whose other residual is exactly $`e_0`$.

The bounded Taylor method of
[Operator Source Compiler Section 10](../../docs/OPERATOR_SOURCE_COMPILER.md#10-classical-table-construction)
evaluates each trigonometric input to $`D`$ bits in $`O(D^4)`$ bit
operations. Native Clifford+T coefficients and all square-root intervals
are also computable in this bound by elementary integer arithmetic.
Choose

```math
D=O\!\left(P+n+\log_2(M+1)+\log_2(\ell+1)\right),\qquad q=P+10.
```

Fixed-precision complex arithmetic for the recorded blocks and the
$`M`$ updates, the rational interval operations above, and writing the
$`O(Nq)`$ source bits give the conservative deterministic bound

```math
O\!\left(NB_{\rm in}+(N+J+M)D^4\right)
```

bit operations. This is polynomial in the specified input and output
lengths. It includes no search over an Euler grid, no exact-zero oracle,
and no inverse coefficient-gap dependence. A supplied coefficient
evaluator with no time bound, or a separately selected native-word
synthesizer, does not acquire a polynomial bound from this statement.

## 6. Small systems at the banked reservation

For $`n\ge6`$, the existing constant-sector construction applies directly:
with two occupied compiler flags and the separate protocol branch, fix
eight address bits, leave $`k=n-6`$ free, and use dirty core
$`q+1=P+11`$, $`k`$ selectors, one helper, and one borrowed signal.
This totals $`P+n+7`$ arbitrary dirty wires. Source-bank work, when used,
is disjoint from these live registers as in
[Complex Coarse Compiler Section 8](COMPLEX_COARSE_COMPILER.md#8-additional-dirty-banks-improve-fine-state-preparation).

At $`n\le5`$ and the banked reservation $`b\ge2(P+n+7)`$, fix all
$`n+2`$ address bits instead. There are at most $`128`$ sectors and
no free address: $`k=0`$. The core, helper, and borrowed signal need only

```math
(P+11)+1+1=P+13\le2(P+n+7).
```

Each one-row mask is emitted directly as its known Pauli word, costing
$`O(P)`$ Clifford gates; it requires neither table selectors nor query
banks. The source width and the full-operator error are unchanged.
Sector errors take a maximum, not a sum, and the bounded predicates use
the already counted arbitrary helper. Thus the same algebraic table
construction and two-flag amplification give polynomial preprocessing
for all small systems at the banked budget, without invoking a fine
single-qubit synthesis routine. This does not change the separately
retained finite-size fallback at the minimum unbanked reservation.

## 7. Bounded certificates

[test_residual_table_preprocessing.py](../../tests/test_residual_table_preprocessing.py)
checks integer-square-root enclosures, exact rational phase identities,
small-radius and unit-circle boundaries, both sides of the half-phase
cut, and the stated finite error and workspace inequalities. These are
classical arithmetic certificates. The separate
[native residual bridge](NATIVE_RESIDUAL_ROTATION.md) now consumes these
intervals to emit an unaddressed row or an enabled two-row table with a
full-operator certificate. General tables and the state compiler remain implementation
work; neither check proves arbitrary native-word synthesis complexity.
