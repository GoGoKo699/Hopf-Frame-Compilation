# A gauge-fixed complex coarse circuit and state preparation

[Borrowed native compiler](BORROWED_WORKSPACE_COMPILER.md) · [State-only construction](STATE_ONLY_COMPILER.md) · [Coarse-frame decoder](COARSE_FRAME_QBP.md)

The phase-dressed Hopf family admits an actual logical native coarse circuit
with linear T-count and exact dirty return. Fixing the common phase removes
the scalar-phase synthesis problem for this state-based task. Its remaining
diagonal factors are determinant-one prefix multiplexors, so the retained
dirty reflection interpreter applies. The resulting complex residual also
fits the two-flag state construction with the existing workspace.

This is a gauge-fixed state and coarse-frame result. It does not approximate
the originally prescribed complete frame with its literal common phase,
and it does not resolve fine complete-frame compilation.

The [bounded native example](NATIVE_COMPLEX_COARSE_QBP.md) emits these
phase-prefix selections and both gradient streams for an exact two-qubit
target. It checks literal phases and arbitrary dirty inputs; the general
fine residual-table construction below remains analytic.

## 1. Target gauge and workspace

Let $`N=2^n`$, $`n\ge1`$, and let $`W_{\mathbb R}`$ be the real Hopf
frame specified by its angle tuple. Supply one real phase representative
$`\varphi_x`$ per leaf, with certified classical evaluation. Set

```math
D_\varphi=\operatorname{diag}(e^{i\varphi_x}),\qquad
\mu=\frac1N\sum_x\varphi_x,\qquad
D_0=e^{-i\mu}D_\varphi,\qquad W'=D_0W_{\mathbb R}.
```

All means below use these supplied representatives. There is no phase
unwrapping, equality test, or choice of an exact principal argument.
The gauge-fixed target state is
$`|\psi'\rangle=W'|0^n\rangle=e^{-i\mu}D_\varphi W_{\mathbb R}|0^n\rangle`$.
Use

```math
0\lt\eta\le1/64,\qquad
L=\max\{6,\lceil\log_2(1/\eta)\rceil\}\ge n,\qquad
b\ge L+n+7,\qquad
\delta_c=\min\{1/64,1/(4\sqrt N)\}.
```

Sections 2–5 construct a circuit implementing $`C\otimes I_b`$ exactly,
with no initialized work, such that

```math
\|C-W'\|\le\delta_c,\qquad T(C)=O(N),\qquad G(C)=O(Nn).
```

Here C is the actual recorded logical native unitary, not the ideal phase
diagonal or an accepted block. Its dirty return is exact on every input,
including arbitrary reference correlations. Section 6 prepares psi', or
psi' coherently with $`C|0^n\rangle`$, to error eta using two initialized
compiler flags and the same dirty pool, with $`T=O(N+L)`$, $`G=O(NL)`$.
The protocol branch and any observable work are counted separately.

## 2. Exact prefix phase factors

For each binary prefix v, let $`\mu_v`$ be the arithmetic mean of phases
in its subtree, with $`\mu_\emptyset=\mu`$ and $`\mu_x=\varphi_x`$
at a leaf. At an internal node put

```math
\gamma_v=\frac{\mu_{v1}-\mu_{v0}}2,\qquad
R_z(\gamma_v)=e^{-i\gamma_v Z}
 =\operatorname{diag}(e^{-i\gamma_v},e^{i\gamma_v}).
```

At depth $`d=0,\ldots,n-1`$, apply this two-by-two table to the next
logical bit, controlled by its d-bit prefix. Denote the layer by $`P_d`$.
It acts for **every** suffix assignment, rather than only the Hopf anchor
pair. On either child, its phase increment is
$`\mu_{v x_{d+1}}-\mu_v`$. These increments telescope along a leaf path:

```math
\prod_{d=0}^{n-1}e^{i(\mu_{x_1\cdots x_{d+1}}
                            -\mu_{x_1\cdots x_d})}
=e^{i(\varphi_x-\mu)},\qquad
D_0=P_{n-1}\cdots P_0.
```

This is the standard cascade of uniformly controlled phase rotations;
see [Möttönen et al., Section III, Eqs. (4), (5), and (7)](https://arxiv.org/pdf/quant-ph/0407010v1).
The displayed derivation fixes the signs and half-angle convention used
here. The native cost and dirty-input implementation are supplied next.

## 3. Literal determinant-one row words

Allocate each phase layer operator error

```math
\epsilon_d=\delta_c\,2^{d-n-1},\qquad
\sum_{d=0}^{n-1}\epsilon_d\lt\delta_c/2.
```

For a desired real rotation $`R_y(\theta)`$, synthesize an actual native
word Q approximating $`R_y(\theta/4)`$ up to scalar phase within
$`\epsilon_d/4`$. As in the
[borrowed construction](BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector),
form

```math
A=XQ^\dagger XQ,\qquad F=A^2.
```

Exactly $`\det A=\det F=1`$. The unknown scalar phase cancels, and
$`\|F-R_y(\theta)\|\le\epsilon_d`$. With the fixed Clifford
$`B=HS^\dagger`$, $`BYB^\dagger=Z`$, so $`BFB^\dagger`$ is a
literal determinant-one native approximation to $`R_z(\theta)`$.
Every row word has length

```math
w_d=O(1+\log_2(1/\epsilon_d))=O(n).
```

Each such actual word has an exact linear-length rewrite into
$`\{T^jHT^{-j}:0\le j\lt8\}\cup\{Z\}`$. These are involutions,
and the rewrite retains scalar phases. No row-dependent phase may be
discarded. The actual row need not be diagonal, or become diagonal after
the target Clifford; only its unitary approximation guarantee is used.

## 4. External dirty banks and the linear coarse count

At depth d there are $`S_d=2^d`$ rows. Choose the power of two

```math
\lambda_d=2^{\left\lfloor\log_2
                 \min\{\sqrt{S_d},b-d\}\right\rfloor}.
```

The reservation uses $`d-\log_2\lambda_d`$ arbitrary external selectors
and $`\lambda_d`$ arbitrary external bank bits. Conservatively
$`d+\lambda_d\le b`$; here $`b-d\ge n+8\ge1`$. Logical prefix
bits serve only as unchanged addresses. No precision core, initialized
bank, or additional predicate helper is needed for this full multiplexor.

Use the [exact reflection interpreter](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
with its optional external control omitted. For each reflection, a dirty
table load and bank routing select $`z\oplus f(x)`$; the matched second
use selects z. Their product is the prescribed reflection to power
$`f(x)`$. Complete this load/use/unload echo before the next reflection.
All banks and selectors return exactly, even though successive target
reflections need not commute. Controlled reflections and Fredkin routers
have charged exact native implementations. Thus an entire depth costs

```math
T_d=O\!\left(w_d(S_d/\lambda_d+\lambda_d+1)\right),
\qquad G_d=O(S_dw_d).
```

Let $`B_c=b-n+1`$. Since $`b-d\ge B_c`$ and rounding to a power of
two loses at most a factor two,

```math
\sum_d T_d
=O\!\left(\frac{Nn}{B_c}+n\sqrt N\right)=O(N),\qquad
\sum_d G_d=O(Nn).
```

The last equality uses $`B_c\ge n+8`$ and $`n\sqrt N=O(N)`$.
It includes the depth-zero row, which may also be emitted directly.
The constant rounding and the extra precision for the phase error split
affect word lengths, not the dirty reservation.

Compile the real frame at error $`\delta_c/2`$ by the borrowed theorem,
obtaining the exact-return logical word $`C_{\mathbb R}`$. Its precision
is $`O(n)`$, hence its T-count is $`O(N)`$ and G-count $`O(Nn)`$ in
this pool. Write $`\widehat P_d`$ for each actual native phase layer and
set

```math
C_P=\widehat P_{n-1}\cdots\widehat P_0,\qquad
C=C_PC_{\mathbb R}.
```

Unitary telescoping gives $`\|C-W'\|\le\delta_c`$. Each complete
subroutine is exactly its logical unitary tensor identity on the dirty
pool. Their composition and the actual inverse C-dagger therefore have
that same exact full-input return property. No approximately returned
phase-source work is substituted for it.

## 5. Recorded coefficients and classical application

Store the actual two-by-two row words of every $`\widehat P_d`$ and the
actual local words of $`C_{\mathbb R}`$. Their matrix entries are given
by finite native words and admit certified evaluation. A phase layer has
$`N/2`$ disjoint two-mode blocks: each prefix row is repeated across all
suffixes. Its actual row can mix the target bit. Applying the layer to a
classical amplitude vector thus costs $`O(N)`$ scalar arithmetic, and
applying C or its actual inverse costs $`O(Nn)`$.

This computation uses the recorded logical words. It neither assumes that
the approximate phase layers remain diagonal nor simulates the dirty
Hilbert space. Computing phase means, native words, coefficient enclosures,
and any required extra classical digits is separate preprocessing. No
uniform bit-complexity claim is made for arbitrary supplied computable
angles or phases.

## 6. Two-flag complex state and common-reference preparation

Define the exact residual $`\phi=C^\dagger\psi'`$. It is normalized and

```math
\|\phi-|0^n\rangle\|\le\delta_c\le1/(4\sqrt N).
```

The [state-only two-flag word](STATE_ONLY_COMPILER.md#3-two-flags-give-an-exact-half-amplitude-state)
already permits complex residual coefficients. For the common-reference
version, let unchanged protocol branch j select

```math
\phi_0=|0^n\rangle,\quad \phi_1=C^\dagger\psi',\qquad
\phi_j=a_j|0^n\rangle+v_j,\quad (v_j)_0=0.
```

On compiler target t, addressed by $`(s,j,x)`$, use

```math
U(z)=\begin{pmatrix}z&-\sqrt{1-|z|^2}\\
                    \sqrt{1-|z|^2}&\overline z\end{pmatrix},\qquad
z_{0,j,x}=a_j,\quad z_{1,j,0}=0,\quad
z_{1,j,x}=\sqrt N\,(v_j)_x\ (x\ne0).
```

The word $`Q=H_sM C_{s=1}(H^{\otimes n})H_s`$ has accepted amplitude
$`\phi_j/2`$ for both branches. The one-step state amplification
$`-Q R_{\rm init}Q^\dagger R_{\rm good}Q`$ therefore prepares either
residual with identical literal phase. The initial reflection excludes
the protocol branch; both reflections include identity on arbitrary dirty
inputs. Apply the same C afterward. No uncontrolled branch phase is
allowed: the resulting targets are specifically $`C|0^n\rangle`$ and
psi', with a relative plus sign for a plus branch input.

For $`n\ge6`$, use the certified Euler and native rotation error budget
from [the state-only construction](STATE_ONLY_COMPILER.md#5-native-table-classical-rounding-and-the-literal-dirty-budget):
$`q=L+10`$, $`\|\widehat Q-Q\|\le130\,2^{-q}`$. Fix eight of
the $`n+2`$ address bits, giving 256 invariant sectors and $`k=n-6`$
free selectors. The dirty reservation is exactly

```math
(q+1)_{\rm core}+k_{\rm selectors}
 +1_{\rm predicate\ helper}+1_{\rm signal}=L+n+7.
```

The synthesis signal is arbitrary dirty work; the two initialized flags
s and t may remain occupied. Sector errors take a maximum because inactive
sectors are exact identities. Use the actual native inverse of Q. The
complete final state-isometry error is at most
$`390\,2^{-q}\lt2^{-L}\le\eta`$. It includes flag leakage and dirty
return against arbitrary branch-and-dirty inputs and external references.
There is no added coarse error, since $`C\phi_1=\psi'`$ exactly.

For $`n\le5`$, directly compile the at most 128 rows of U(z) instead.
Certified SU(2) Euler approximation gives three determinant-one rotations
per row; phase-cancelled native real words and fixed Clifford conjugations
approximate them with length $`O(L)`$. Rewrite the actual words in the
same reflection alphabet and add the at most seven address literals by
bounded controlled-Z conjugations. The exact predicate construction uses
one arbitrary returned helper, with a bounded number of native gates per
reflection. It assumes no initialized scratch or native controlled T.
Each row and the classical completion can satisfy the same Q error budget.
The same actual-inverse amplification and common C then apply, with
$`T,G=O(L)`$ for these finitely many n. This fallback does not invoke a
two-real-state theorem for a complex target.

In both cases the total counts are $`T=O(N+L)`$ and $`G=O(NL)`$.
Only the initial zero-system subspace is promised for the fine preparation;
there is no fine complete-frame claim. Classical square-root and Euler
certification, including zero tails and $`|a_0|=1`$, follows the existing
state proof and requires no exact-zero decision.

## 7. The gauge preserves the requested energy gradients

For Hermitian O, replacing the physical target psi by
$`\psi'=e^{-i\mu}\psi`$ preserves its energy. Magnitude derivatives
satisfy $`d'_j=e^{-i\mu}d_j=a_jW'|\lambda(j)\rangle`$ because mu
depends only on the phase tuple. Their raw energy gradients are unchanged.
For a leaf-phase coordinate,

```math
\partial_{\varphi_x}\psi'=i(P_x-I/N)\psi',\qquad
P_x=|x\rangle\langle x|.
```

The extra common-phase term contributes zero to
$`2\operatorname{Re}\langle\partial_{\varphi_x}\psi'|O|\psi'\rangle`$,
since $`\langle\psi'|O|\psi'\rangle`$ is real. Equivalently, one may
freeze the chosen scalar gauge at the current tuple when forming the
derivative records. This does not differentiate the discrete native
synthesis word. Observable compilation and sampling remain separately
charged. Preparing the unrotated psi in only one interference branch would
change that branch's relative phase and is not licensed by this argument.

## 8. Additional dirty banks improve fine state preparation

The real and gauge-fixed complex state constructions, including their
coherent common-reference versions, admit a banked refinement. Retain
$`L\ge n`$, two initialized compiler flags, and the preceding complete
state-isometry error eta. If

```math
B_0=L+n+7,\qquad b\ge2B_0,
```

then

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}+n\sqrt N\right),
\qquad G=O(NL).
```

The protocol branch is still separately counted. This refines preparation
of the requested state on the initialized system; it is not a new bound
for fine compilation of its prescribed complete frame.

### Disjoint banks and the full borrowed-signal contract

For $`n\ge6`$, use the same $`q=L+10`$, core width $`m=q+1=L+11`$,
and $`k=n-6`$ free address bits as in Section 6. Eight fixed address
literals give 256 invariant sectors for the coherent two-state table.
The core, k selectors, predicate helper, and arbitrary synthesis signal
occupy exactly $`B_0`$ dirty wires. In particular $`m\le B_0`$.
Reserve the disjoint additional pool

```math
B_{\rm bank}=b-B_0\ge b/2\ge m.
```

For an S-row, m-bit programmable mask table, select a power-of-two number
of m-bit banks by

```math
\lambda=2^{\left\lfloor\log_2
 \max\{1,\min(S,\sqrt{S/m},B_{\rm bank}/m)\}\right\rfloor}.
```

Thus $`1\le\lambda\le S`$ and $`\lambda m\le B_{\rm bank}`$.
Apply the [exact whole-word dirty-bank query](OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost).
Its high-address traversal, low-address bank router, and matched XOR echo
implement precisely the original query on the arbitrary operator core,
tensored with identity on all added banks. The base selector reservation
covers the remaining address; no extra selector register is hidden in
the bank pool. The query costs

```math
T_{\rm query}=O(S/\lambda+\lambda m)
 =O\!\left(\sqrt{Sm}+m+\frac{Sm}{b}\right),
\qquad G_{\rm query}=O(Sm).
```

All banks and selectors return exactly after the query. The X- and Z-mask
queries run sequentially and reuse that returned bank pool. Their output
is the separate occupied dirty core. The core, signal, and predicate helper
are not counted again as available banks while the primitive is active.

This substitution does not assume a clean amplification signal. On the
complete input space it replaces a query by exactly the same query tensor
identity on the new banks. Hence the actual real-rotation word retains its
exact commutation with the signal's X, including its five-call amplification
and inactive-predicate action. The
[borrowed-signal extension](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
therefore retains the full-operator error $`43\,2^{-q}`$ with that signal
arbitrary and the added banks included in the dirty identity. Fixed target
Cliffords give the same conclusion for Rz. The literal scalar-phase word
is not used: the residual completion U(z) has determinant one and needs
only these three Euler rotations.

There are a constant number of masks and sources per rotation; the source
words still cost $`O(q)`$. Three Euler factors, three Q appearances, and
a fixed number of sectors consequently give fine preparation cost

```math
T_{\rm fine}=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G_{\rm fine}=O(NL).
```

Every inactive sector remains exactly identity on the entire enlarged
work space, so sector errors take a maximum. Euler certification and
actual-inverse state amplification keep the same final error
$`390\,2^{-q}\lt\eta`$. For $`n\le5`$, the direct bounded-control
Euler construction in Section 6 costs $`O(L)`$ and fits this stronger
reservation without requiring a source or a negative address length.

### The exact coarse circuit remains charged

The fine banked primitive returns its core and signal approximately. It
cannot replace the exact-return logical coarse C used to define residual
coefficients and the measurement basis. The retained real borrowed
construction and Sections 3–4 for the complex phase factors give, in the
present workspace regime,

```math
T_C=O\!\left(\frac{Nn}{n+b}+n\sqrt N\right),
\qquad G_C=O(Nn).
```

For the complex phase factors, $`b-n+1=\Theta(n+b)`$ here, so their
earlier denominator gives the same bound. All recorded coarse factors
still implement their exact logical action tensor identity on dirty work.
Since $`L\ge n`$, their first term is absorbed by $`NL/b`$.
The exact state reflections cost $`O(n^2)`$, absorbed by $`n\sqrt N`$.
Adding these charges proves the displayed banked state bound.

In particular, for $`L\ge n^2`$ the coarse term is absorbed by
$`\sqrt{NL}`$, yielding

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right).
```

For smaller L, the $`n\sqrt N`$ contribution must be retained unless
another displayed term absorbs it. The result compares constructive upper
bounds under the declared workspace; it is not a lower bound or a claim
of optimal gradient cost.
