# A binary shared-source stream through the whole upper frame

[Binary terminal circuit](SHARED_BINARY_TERMINAL.md) · [Precision-uniform compressed columns](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) · [Retained uniform compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary)

**Local result R68.** One logarithmic-width binary source can serve all
upper frame layers.
Choosing the protected bottom bank from $`\log(n+L+2)`$ extends the
native shared-source construction to an absolute exponential precision
range. The proof keeps source cardinality J explicit: only its binary
address is polynomial in n, while the paid source preparation is $`O(J)`$.
The exact native kernel and recovery are reused from the binary terminal
proof; the new work here is the uniform resource and width analysis.

The later [one-clean theorem R72](ONE_CLEAN_DEPTH.md) relocates the
padding-rejection bit into this protected logical bank and proves a
one-clean bottom schedule and width splice. This chapter keeps its
original two-clean theorem.

## 1. The claimed range and complete-circuit target

There are absolute constants $`A_0,C,n_0`$ and

```math
\kappa=\frac1{4A_0}\gt0
```

such that, for every integer $`n\ge1`$, integer
$`6\le L\le2^{\kappa n}`$, $`N=2^n`$, and $`b\ge17(L+n+7)`$,
every prescribed complete real Hopf frame W admits one coherent native
circuit with two external clean flags, at most b arbitrary dirty qubits,
and

```math
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-L},
```

```math
T\le C\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad G\le CNL,
```

```math
D_T\le C\left(\frac{NL}{b^2}
 +\min\{nL,(n+L)\log_2(n+L+2)\}\right).
```

All constants are independent of n, L, b, and the prescribed frame.
The complete-isometry norm includes every logical column, initialized
work return, arbitrary dirty inputs and references, and literal phase.
All tables are finite compiled data, not supplied quantum evaluators.
Certified classical coefficient computation has no efficiency claim.
There are no measurements, resets, or supplied precision states.

The theorem is deliberately restricted to the displayed exponential
precision range. The all-precision baseline retains an additive $`nL`$
T-count outside ranges where it can be absorbed; no unrestricted
optimal-count splice is inferred here.

Put

```math
\eta_L=2^{-L},\qquad
\ell=\log_2(n+2),\qquad x=n+L+2,\qquad \lambda=\log_2x,
```

```math
J=2^{\lceil\log_2(L+\lceil\ell\rceil+20)\rceil}.
```

Thus $`J=\Theta(L+\ell)`$ with absolute constants and $`J\le C_Jx`$.
Choose

```math
R_0=\min\{n,\lceil A_0\lambda\rceil\},\qquad R=n-R_0.
```

The R upper layers form the shared-source stream, in their original
shallow-to-deep order. The final $`R_0`$ layers use the retained capped
layerwise source. There is no compressed-column early segment and no
binary phase-gradient arithmetic in this construction.

## 2. Reuse the native stream with the new precision

Only the parameter and resource choices change. The exact circuit is
the [sixteen-type residual construction](SHARED_BINARY_TERMINAL.md#2-sixteen-residual-types-with-no-endpoint-register),
[paid source and local kernel](SHARED_BINARY_TERMINAL.md#3-paid-geometric-source-and-the-complete-local-kernel),
[history compression](SHARED_BINARY_TERMINAL.md#4-local-and-stream-accepted-column-estimates),
[global padding and inactive cancellation](SHARED_BINARY_TERMINAL.md#5-global-padding-and-full-occupied-space-inactive-identity),
and [complete-isometry recovery](SHARED_BINARY_TERMINAL.md#6-full-isometry-recovery-and-return)
from Sections 2–6 of the binary terminal proof. Their native identities
and error estimates require $`R\ge1`$ and the displayed source caps;
they do not require the earlier $`R=O(L+\log n)`$ or $`L=O(n)`$.
Sections 7–8 of that note are replaced by the resource argument below.

For clarity, the reused interface has the following literal features.
Save $`H=[Z=0^{R_0}]`$ in the first external clean flag. Each local
suffix predicate $`G_i`$ tests only logical bits below its target and
above Z, excluding the target and the protected bank. Its computation
and inverse use separate returned scratch and direct H control. Coarse
native words act on the target, conditional on its preserved prefix and
$`G_i`$; their residual dictionary has $`K=16`$ labels
$`(u,v,\omega)`$ with $`u,v\in\{0,1\}`$ and
$`\omega\in\{1,-1,i,-i\}`$. There is no full endpoint register.

Choose exactly

```math
\pi=2^{-\lceil\log_2(8R)\rceil},\qquad
\gamma=(1-\pi)^R,\qquad
\varepsilon=\frac{\pi}{(1-\pi)KJ},\qquad
\delta_c=\varepsilon/4,\qquad M_s=2J.
```

Then $`R\pi\le1/8`$ and $`7/8\le\gamma\le1`$. Actual,
literal-phase coarse words at accuracy $`\delta_c`$ have length

```math
w=O(\log(R+1)+\log J)=O(\lambda).
```

Compute the residual from those actual words. Certified nonnegative
lower rounding gives its J binary digits, with coefficient error
less than $`2\varepsilon2^{1-J}`$. Local inactive rows are literal
zero. Coarse approximation is therefore not a separate final error.
Prefetch the entire coarse word into conditional-zero work; every
consumed symbol keeps direct H and local-suffix controls. Its actual
unquery returns its program and helpers on every logical, source, and
failure input.

The [exact capped binary source](../../docs/FAULT_TOLERANT_COMPILER.md#3-an-exact-geometric-source-of-logarithmic-width)
has cap $`M_s=2J`$, paid count and conservative depth $`O(J)`$, and only
$`O(\log J)`$ initialized source/scratch width. The same forward
downward shift $`S_{2k}`$ is used for every matrix-unit orientation.
The digit query reads only the preserved tuple
$`(H,\text{prefix},G_i,u,v,\omega,k)`$, so its actual unquery remains
valid on all rejection branches.

Write F for retained private failures/labels, and $`g`$ for the
prepared source column. The accepted local contraction is exactly

```math
T_i=(1-\pi)I+\frac{\pi}{KJ}\sum_{u,v,\omega,k}
S_{2k}\otimes f_{i,u,v,\omega,k}\,\omega|u\rangle\langle v|.
```

On $`G_i=0`$ this is $`(1-\pi)I`$, because the zero-digit rare
branch still rejects. This common attenuation is essential. A
$`\lceil\log_2R\rceil`$-bit history counter records nonzero F between
stages, with at most $`R-1`$ increments and no wrap. For $`R=1`$ omit that
counter. Final F-zero/history-zero projection inserts the intermediate
F-zero projectors without projecting the source.

Source-preparation scratch and the coarse/predicate bank remain
separate from F/history. Rare-event and F-nonzero predicates use the
charged returned-helper conjunctions. No source preparation borrows
retained failures, and no many-controlled gate is free.

## 3. The same complete error and inactive-input contract

Let $`E_g=|g\rangle\otimes I`$ and
$`D_i=(1-\pi)C_i^\dagger W_i`$. The linked source-column estimate,
with $`e_s\lt2\,2^{-J}`$, is

```math
\|T_iE_g-E_gD_i\|
\le\pi e_s+\frac{2\pi}{J}2^{1-J}.
```

Contractive telescoping over all R stages gives

```math
\|V_{\rm acc}E_g-E_g\gamma W_{\rm upper}\|\le\xi,\qquad
\xi\le R\pi\left(e_s+\frac2J2^{1-J}\right)
\le\frac38\,2^{-J}\le\eta_L/1024.
```

The coherent direct sum over prefixes introduces no row-count factor.
At the final source boundary use the actual preparation inverse.
The independent padding sampler has power-of-two cap
$`L+12\le M_p\lt2(L+12)`$, targets the rational probability
$`1/(2\gamma)`$, and has probability error
$`\epsilon_0\le\eta_L/2048`$. Its native count/depth is $`O(L)`$ and
initialized width is $`O(\log L)`$.

Use exactly the nested chronological preparation and actual inverse
order of terminal Section 5. On H zero the stream, history and padding
operations are identity on arbitrary occupied inputs, so those
unconditional preparation words cancel in reverse order. The completed
recovery is exactly identity on the entire H-zero space; no cleanliness
claim is made for an individually disabled source preparation.

The physical zero-work reflection includes the protected bank, F,
source, history, padding, approximation work and second external flag,
after the actual source inversions; it excludes saved H and all
arbitrary dirty inputs. Thus the terminal Section 6 normalization-two
word gives complete upper error

```math
8(\xi+\epsilon_0)\le3\eta_L/256\lt\eta_L/16.
```

Every inverse is the actual inverse of the finite native circuit.
Finally reverse the original H predicate. Its unitary action propagates
the existing norm error and returns H on ideal outputs, with no reset.
These are full initialized-isometry estimates, including all logical
columns, literal phase, returned work, arbitrary dirty states and
external references.

## 4. Simultaneous clean width and stage depth

All live initialized fields occupy at most $`C_Z\lambda`$ bits
for an absolute $`C_Z`$: both binary sources and their separate preparation
scratch; rare, term, and digit labels; flags; $`\log R`$ history; global
padding; exact arithmetic/history scratch; and the separate coarse
program/predicate bank. Every size is $`O(\log J+\log(R+1))=O(\lambda)`$.
The first external flag saves H and the second may hold padding
rejection. No other externally initialized input is used.

Choose an absolute integer

```math
A_0\ge\max\{24,C_Z+2\},\qquad \kappa=\frac1{4A_0}.
```

When $`R\ge1`$, $`R_0=\lceil A_0\lambda\rceil`$ supplies all
this conditional-zero work. No source register has J initialized
qubits: its register and preparation scratch have logarithmic width.

A local stage has a constant number of shallow queries, each of depth
$`O(\log(n+\log J+2))`$; controlled comparison/subtraction and history cost
$`O(\log J+\log(R+1))`$; coarse interpretation costs $`O(w)`$; and the original
logical-suffix predicate costs $`O(\log(n+2))`$. Rare and F-nonzero events
are decoded by charged returned-helper conjunctions, not treated as
free many-controlled gates. Thus each full stage has depth $`O(\lambda)`$.

The two paid source preparations and inverses cost $`O(J)`$ per complete
stream. Outer recovery uses only three stream appearances. Therefore

```math
D_{T,\rm upper}=O(J+R\lambda)
=O(L+n\lambda).
```

## 5. Precision-aware table and native-work envelopes

Define $`u=n+\lceil\log_2J\rceil+20`$. At absolute layer depth d,
the digit query has at most $`CJ2^d`$ rows and constant output width;
the coarse-word query has at most $`C2^d`$ rows and $`O(w)`$ outputs.
The maximal address length is $`O(u)`$, with no target or source value
appended. The depths run consecutively from zero through $`R-1`$ and
$`R=n-R_0`$.

Apply the [conservative shallow-query bound](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address), including its private
non-Clifford supports for parallel coarse-word outputs. The sum
$`\sum_{d=0}^{R-1}2^{d/2}=O(2^{R/2})`$ and the analogous
geometric row sum give the safe common query envelopes

```math
T_{\rm query},w_{\rm query}
\le C u^6\sqrt J\,\sqrt N\,2^{-R_0/2},
\qquad
G_{\rm query}\le C u^6J\,N2^{-R_0}.
```

The coarse output factor w and every constant call multiplicity fit
$`u^6`$; this is a deliberately conservative bound. In particular J is
not silently replaced by a polynomial in n.

All remaining native source/padding words cost $`O(J)`$. Other local
operations have total count at most $`Cu^6`$ and peak returned dirty work
at most $`Cu^6`$: there are $`O(n)`$ stages, $`O(nw)`$ coarse symbols, and a
polynomial number of conjunctions/ripple gates, each with its charged
safe cubic helper/count bound. Source-preparation work is initialized
logarithmic work already counted in Z, not J dirty helpers. Thus

```math
T_{\rm upper}\le C u^6
 \left(\sqrt J\,\sqrt N\,2^{-R_0/2}+1\right)+CJ,
```

```math
w_{\rm upper}\le C u^6
 \left(\sqrt J\,\sqrt N\,2^{-R_0/2}+1\right),
```

```math
G_{\rm upper}\le C u^6\left(JN2^{-R_0}+1\right)+CJ.
```

Every returned pool is reused only after actual return and is disjoint
from simultaneously live query output/address and initialized work.
The only unweighted precision count is the paid $`O(J)`$ source/sampler
cost, not a polynomial-in-L factor multiplying the whole circuit.

## 6. One absolute dimension cutoff

In the range $`L\le2^{\kappa n}`$, $`u=O(n)`$ uniformly.
For example $`J\le C(L+n+2)`$ implies
$`\log_2J\le C+\kappa n+\log_2(n+2)=O(n)`$.
Also $`u\le Cx`$ and $`J\le Cx`$.

Choose an absolute $`n_0`$ large enough that, for all $`n\ge n_0`$,

```math
n+2\le2^{\kappa n},\qquad n\ge4(A_0+1).
```

Then $`x\le2^{\kappa n+1}`$, so

```math
R_0\le A_0(\kappa n+1)+1
\le n/4+A_0+1\le n/2.
```

In particular the nontrivial stream exists uniformly at these dimensions.
Since $`R_0=\lceil A_0\lambda\rceil`$,

```math
u^6\sqrt J\,2^{-R_0/2}
\le Cx^{13/2-A_0/2}\le Cx^{-11/2},
```

```math
u^6J\,2^{-R_0}\le Cx^{7-A_0}\le Cx^{-17}.
```

Thus the query width is at most $`C\sqrt N(n+2)^{-11/2}`$.
Increase $`n_0`$ so that

```math
C\sqrt N(n+2)^{-11/2}+Cn^6
\le\tfrac12\sqrt{N/n},\qquad
\sqrt{N/n}\ge32c_I(n+3).
```

Here C is fixed after the native constants, and $`c_I`$ is the absolute
nonuniform-indicator constant. These inequalities eventually hold
and do not depend on L, b, or the prescribed frame.
The upper query/helper pool fits every $`b\ge\sqrt{N/n}`$.

For counts, $`J=O(L+\ell)`$, $`L\le N`$, and exponential
growth imply, uniformly in this range,

```math
J+n^6=O(\sqrt{NL}),\qquad J+n^6=O(NL).
```

For the first bound use $`L\le\sqrt{NL}`$ and boundedness of
$`n^6/2^{n/2}`$, with $`L\ge6`$. The Clifford bound is weaker.
Consequently $`T_{\rm upper}=O(\sqrt{NL})`$ and
$`G_{\rm upper}=O(NL)`$ with absolute constants.

## 7. The bottom layers retain optimal counts

Use the [retained capped layerwise allocation](SHARED_BINARY_TERMINAL.md#8-bottom-and-terminal-composition) at bottom height j:

```math
m_j=L+6+\min\{j,\lceil\log_2(8n)\rceil\},\qquad
Q_j=4N2^{-j},\qquad 1\le j\le R_0.
```

The cap gives $`m_j=O(L+\ell)`$ and the geometric bound
$`m_j\le L+6+j`$ gives

```math
\sum_jQ_jm_j=O(NL),\qquad
\sum_j\sqrt{Q_jm_j}=O(\sqrt{NL}).
```

Its actual base is $`B_0'=L+n+9`$. With
$`B=b-B_0'`$, the literal supplied threshold gives

```math
B\ge b-b/17-2\ge b/2,
```

since $`b\ge17(L+n+7)\ge238`$.
On the large-width branch the cutoff in Section 6 then ensures
$`B\ge16c_I(n+3)`$. No assumption $`L\le n`$ is used here.
The [count-efficient rectangular queries](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence) therefore fit and
have their original $`Q_jm_j/B`$ count and $`Q_jm_j/B^2`$ depth
penalties.

Unweighted source and query costs sum to
$`O(R_0(L+\ell))`$. The important uniform absorption is

```math
\frac{nL}{\sqrt{NL}}
=\frac{n\sqrt L}{2^{n/2}}
\le n2^{-(1-\kappa)n/2}=O(1).
```

Since $`R_0\le n`$, this absorbs $`LR_0`$. The remaining
$`R_0\ell`$ and polynomial predicate counts are absorbed into
$`\sqrt{NL}`$ by exponential growth. Thus

```math
T_{\rm bottom}=O(\sqrt{NL}+NL/b),\qquad G_{\rm bottom}=O(NL).
```

Because $`R_0\le\lceil A_0\lambda\rceil`$ and $`\ell\le n+2`$,

```math
D_{T,\rm bottom}
=O(NL/b^2+R_0(L+\ell))
=O(NL/b^2+(n+L)\lambda).
```

Its complete-input error is below $`\eta_L/4`$. Upper and bottom
dirty pools are sequential, so their maximum is the required width.
Full-isometry telescoping gives total error less than
$`\eta_L/16+\eta_L/4\lt\eta_L`$, including all source/history,
flag, logical-work, and dirty/reference return. The entire frame,
including its continuation columns, is the prescribed product of
these upper and bottom logical layers.

Combining the same circuit's resources establishes the large-width,
large-dimension target

```math
T=O(\sqrt{NL}+NL/b),\quad G=O(NL),\quad
D_T=O(NL/b^2+(n+L)\lambda).
```

## 8. All widths, no-middle cases, and finite dimensions

The retained
[uniform two-clean compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary)
has absolute simultaneous bounds

```math
T=O(\sqrt{NL}+NL/b+nL),\quad G=O(NL),\quad
D_T=O(NL/b^2+nL)
```

and the required complete-isometry contract at the same literal
allocation. The Section 7 inequality absorbs $`nL`$ into $`\sqrt{NL}`$
throughout the present exponential range, for every n, not just
asymptotically.

For $`n\ge n_0`$ and $`b\lt\sqrt{N/n}`$, use that retained
compiler. Then $`nL\le NL/b^2`$, so its depth meets the new
large-width envelope as well as its own $`nL`$ bound.

If $`R_0=n`$ there is no source stream. The definition implies
$`n\le\lceil A_0\lambda\rceil`$ and hence
$`nL\le(A_0+1)L\lambda`$. Thus the same retained compiler
already has the new depth envelope and optimal counts in this case.

For any remaining $`n\lt n_0`$, use that compiler at every eligible
width. Its absolute constants and
$`nL\le n_0L\le n_0(n+L)\lambda`$ give the desired depth after
one absolute enlargement. Its $`nL`$ count is absorbed by the same
uniform inequality, not by a precision-dependent finite exception.

Finally the retained compiler itself supplies depth
$`O(NL/b^2+nL)`$ with the desired counts throughout the range.
Selecting between these two complete circuits gives the stated minimum.
No allocation is assembled by mixing the count of one circuit with
the depth of another.

The inherited count lower bound contains
$`\sqrt{NL}+L+NL/b`$. Since $`L\le N`$ here, its L term is
absorbed into $`\sqrt{NL}`$, so the stated count is worst-case optimal
throughout this range. The depth statement is an upper bound;
no unrestricted depth optimum is claimed.

## 9. Relation to the low-precision theorem

The [existing compressed-column theorem](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) remains separately stronger
for $`6\le L\le cn`$:
$`D_T=O(NL/b^2+\min\{nL,(\sqrt n+L)\ell\})`$.
This whole-upper-frame stream does not remove its $`\sqrt n`$
normalization cost. It bypasses the gradient reservation when seeking
higher precision.

For any fixed polynomial precision $`L=n^k`$, the new exponential
range eventually includes it. For example $`L=n^2`$ gives
$`D_T=O(Nn^2/b^2+n^2\log(n+2))`$, compared with the retained
additive $`n^3`$ allowance, with the same optimal T and Clifford budgets.
These are asymptotic scope statements, not practical crossover claims.

The [all-precision theorem](ALL_PRECISION_DEPTH.md) combines this scoped
result with separately proved low- and high-precision constructions and
precision splices.
