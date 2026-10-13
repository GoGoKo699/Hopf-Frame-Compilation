# Center-phase reflections: exact image and approximation limits

[Paid encodings](PAID_LOGICAL_ENCODINGS.md) · [Full reflection compiler](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources) · [Endpoint](../../docs/OPEN_PROBLEM.md)

The center-phase word below has an explicitly solvable phase-selection
problem. Its optimum ranges over every phase table, not only a local
chart. Its exactly real outputs form a finite set of directions with one
continuous rotation angle. Generic real directions retain a quadratic
error floor, and an actual two-edge Hopf word has a separate floor against
the larger single-reflection family. These are restrictions on the stated
word, not on all non-Clifford baselines or the linear-T endpoint.

Throughout, $`N=2^n`$ and the original tolerance is $`\eta=2^{-N}`$.
Only T-count is optimized. Literal phases and full logical/work/reference
outputs remain part of the contract. The two diagonal occurrences here
are one programmed table and its actual inverse, not two independent slots.

## 1. The whole word and its global optimum

Let e be the zeroth basis vector, F the real normalized Walsh transform,
and define

```math
R=I-2ee^T,\qquad S=I+(i-1)ee^T,\qquad
Q(d)=F\mathrm{diag}(d)F,\quad |d_k|=1,
\qquad V(d)=S Q(d)R Q(d)^\dagger R S^\dagger.
```

Every table is the exponential of a real XOR-circulant Hermitian matrix.
Subtracting its mean eigenphase changes only a common scalar, which cancels.
Thus arbitrary tables include the original centered exponential candidate.
For $`H_v=I-2vv^\dagger`$ with a unit normal v, direct conjugation gives

```math
w(d)=SFd/\sqrt N,\qquad V(d)=H_{w(d)}R.
```

For the target $`W_v=H_vR`$, the rank-one projector norm identity implies

```math
\|V(d)-W_v\|=2\sqrt{1-|v^\dagger w(d)|^2}.
```

Put $`y=FS^\dagger v`$. The largest overlap with a unit-modulus table is
$`\|y\|_1/\sqrt N`$: triangle inequality proves the upper bound and
$`d_k=y_k/|y_k|`$ attains it at every nonzero coordinate. Zero coordinates
may receive any phase. Therefore

```math
\boxed{\inf_{|d_k|=1}\|V(d)-W_v\|
=2\sqrt{1-\|FS^\dagger v\|_1^2/N}.}
```

The infimum is attained. It allows nonsmooth selection, phases far from
identity, and complex target normals. Approximation of a general normal
can use certified finite phase search with a strict margin; the displayed
formula is not an assumed exact-zero oracle.

## 2. Exact real image and a regular optimizing table

For real unit b with $`b_0=0`$, put

```math
c=\cos(\gamma/2),\qquad s=\sin(\gamma/2),\qquad v=ce+sb,
\qquad W_b(\gamma)=H_vR=\exp\{\gamma(be^T-eb^T)\}.
```

Let $`\lambda_k=\sqrt N(Fb)_k`$. Their mean is zero and their mean square
is one. The optimal overlap is

```math
\mu(b,\gamma)=\frac1N\sum_k\sqrt{c^2+s^2\lambda_k^2}.
```

After a harmless common phase, an optimizing table is

```math
d_k^*=\frac{c+i s\lambda_k}{\sqrt{c^2+s^2\lambda_k^2}}.
```

For $`|\gamma|\le\pi/2`$ the denominators are at least $`1/\sqrt2`$.
Certified arithmetic therefore computes these phases without a singular
selection or nonlinear solver. Table error epsilon changes Q by at most
epsilon and the complete word by at most $`2\epsilon`$.

Strict concavity of the square root shows that, for nonzero s, exact
agreement occurs precisely when all $`\lambda_k^2=1`$. Hence

```math
b=F\sigma/\sqrt N,\qquad \sigma_k\in\{-1,1\},\qquad
\sum_k\sigma_k=0.
```

Conversely every balanced sign table supplies a real unit direction
orthogonal to e, and $`d_k=e^{i\gamma\sigma_k/2}`$ gives the exact
complete plane rotation for every gamma. If an arbitrary V(d) is real,
then its reflection projector is real; a common phase makes its normal
real. The same argument applies. Apart from identity, these flat-Walsh
directions exhaust the real image. At fixed N the direction set is finite,
not a neighborhood of all real directions. Complex intermediate operations
remain permitted and are used by this construction.

The earlier first-order calculation remains a corollary, not an exact
coverage claim. For a real vector a with $`a_0=0`$, define
$`A(a)_{xy}=a_{x\oplus y}`$. Differentiation gives

```math
\left.\frac{d}{dt}\right|_{t=0}
 S e^{itA(a)/2}R e^{-itA(a)/2}R S^\dagger
=ae^T-ea^T.
```

Every real star tangent occurs, but the curve can leave the real target
manifold at second order. The following optimum shows that better phase
fitting does not remove that error in general.

## 3. A sharp quadratic loss after all phase optimization

Choose distinct nonzero labels a,b and direction
$`b_*=(e_a+e_b)/\sqrt2`$. Its squared Walsh coefficients are zero on half
of the spectrum and two on the other half. With $`t=\sin^2(\gamma/2)`$,

```math
\boxed{\inf_d\|V(d)-W_{b_*}(\gamma)\|
=\sqrt{2-2\sqrt{1-t^2}},\qquad
 t\le\inf_d\|V(d)-W_{b_*}(\gamma)\|\le\sqrt2\,t.}
```

A stable expression is $`\sqrt2\,t/\sqrt{1+\sqrt{1-t^2}}`$. The small-angle
loss is $`\gamma^2/4+O(\gamma^4)`$. This is a global lower bound over the
table, not merely a Taylor upper estimate. More generally a fixed direction
with nonzero variance of $`\lambda^2`$ has leading loss
$`\sin^2(\gamma/2)\sqrt{\mathrm{Var}(\lambda^2)}`$; its higher-order
constant may depend on the direction and N. The explicit two-coordinate
formula supplies the uniform witness. This plane target is not asserted to
be an original ordered Hopf word.

## 4. The separate floor for a true ordered Hopf target

Let the (0,2) rotation act first and the (0,1) rotation second, both at
angle theta, with all other original angles zero. In order (0,2,1),

```math
W_\theta=\begin{pmatrix}
c^2&-cs&-s\\s&c&0\\sc&-s^2&c
\end{pmatrix},\qquad c=\cos\theta,\quad s=\sin\theta.
```

For $`|\theta|\le\pi/2`$, even the larger family with an arbitrary complex
unit normal satisfies

```math
\boxed{\inf_{\|w\|=1}\|H_wR-W_\theta\|=1-\cos\theta.}
```

Indeed, $`M=W_\theta R`$ is orthogonal with determinant -1 and eigenvalues
$`-1,e^{i\beta},e^{-i\beta}`$ on the active subspace, where
$`\cos\beta=1-(1-c)^2/2`$ and $`0\le\beta\le\pi/3`$. For an eigenvector
of M with eigenvalue $`e^{i\beta}`$, its distance from any Hermitian
involution is at least the distance to the set of eigenvalues -1 and 1,
namely $`2\sin(\beta/2)=1-c`$. A reflection in the -1 eigenvector attains
it. That normal is proportional to $`e_0+u e_2+u e_1`$ with
$`u=\tan(\theta/2)`$, regular at zero. Additional identity coordinates and
complex normals do not remove the spectral bound.

This also excludes merely replacing F by a different unitary inside the
same single-reflection/fixed-center geometry. It does not exclude multiple
centers, products of reflections, independently programmed slots, or another
whole-word construction. The three-label target is already cheaply compiled
by the disjoint-block theorem.

For an endpoint-scale comparison set

```math
\ell=\lceil\sqrt N\rceil+6,\quad \delta=2^{-\ell},\quad
u=\delta/32,\quad \theta=2\arctan u.
```

The true ordered word has $`\|W_\theta-I\|\le4u=\delta/8`$ and loss
$`1-\cos\theta=2u^2/(1+u^2)\ge u^2`$. The plane witness has loss at least
$`u^2/2`$. For every $`n\ge8`$, $`2\ell+12\lt N`$, so the common floor
$`2^{-2\ell-11}`$ exceeds $`2\eta`$. Thus an eta-accurate implementation
of an ideal member of this family cannot also approximate these targets to
eta. A different physical correction not accurately implementing that ideal
family is outside the statement.

## 5. Native exact subfamily and resource accounting

For any balanced sign table, the exact-real family has only two phase
values. A fixed number of existing literal diagonal and Boolean center-phase
calls gives $`O(N)`$ T-count at the original two-clean and $`N+n+7`$-dirty
allocation. All actual inverses are charged. This implements V(d); where
the global optimum is nonzero its representation error must be added to
native synthesis error. No smaller arbitrary-angle workspace claim follows.

[centered_phase_reflections.py](../../compiler_robust_hopf/centered_phase_reflections.py)
emits a nontrivial exact member using data bits 0,1,2 and clean work 3,4:

```math
\sigma(k)=(-1)^{k_0\oplus k_1k_2},\qquad
b=(0,1,0,1,0,1,0,-1)^T/2,\qquad\gamma=\pi/2.
```

Compute the Boolean sign into one clean bit, apply the literal chronological
word T-dagger, X, T, X for $`e^{i\pi Z/4}`$, and uncompute. Center S/Z
operations compute and uncompute the zero-data predicate using both work
bits. No live conjunction is mistaken for a zero state. Reusing the retained
seven-T Toffoli and expanding X/Z gives 438 elementary H/S/T/CNOT gates,
including 144 T/T-dagger appearances. The count is unoptimized.

The ideal initialized-input target is the real rational matrix

```math
I-ee^T-bb^T+be^T-eb^T.
```

Both work bits return exactly on every logical input. Occupied-work inputs
have a specified unitary action, not necessarily the same target. The
[tests](../../tests/test_centered_phase_reflections.py) use the existing
$`\mathbb Q(\sqrt2,i)`$ ring to verify all initialized columns and the full
actual inverse exactly. Bounded floating checks separately test the global
optimum, attained floors, real family, and first derivative. They do not
certify endpoint-precision general synthesis; only this finite member is
emitted. The large JSON gate list is regenerated from the function, not
copied into the repository.

Householder reflections are established background; see
[Ivanov, Kyoseva and Vitanov](https://doi.org/10.1103/PhysRevA.74.022323).
The projector and phase-alignment arguments used here are explicit above.
The existing [four-fine-slot reflection construction](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)
is a positive alternative for arbitrary real normals, using a different
chart and a priced coarse basis. This restricted optimum does not invalidate
it or establish an unrestricted compiler lower bound.
