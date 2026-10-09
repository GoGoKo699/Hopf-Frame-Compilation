# Coupled tree resolvents and constant-normalization scattering

Geometrically decreasing local defects give a height-independent norm
bound for the coupled tree inverse. Its finite power series has explicit
tail and input-support bounds, and the whole inverse has an ideal unitary
scattering realization at normalization three. The construction and its
native precision cost are distinct: the depth-ordered implementation
retains an $`O(N+nL)`$ T-count.

[Residual assembly](RESIDUAL_ASSEMBLY.md#6-what-this-assembly-resolves) ·
[Collective refinement](COLLECTIVE_PRECISION_REFINEMENT.md) ·
[Research map](../README.md)

## 1. Graded tree shifts

For a binary tree of height n with $`N=2^n`$ leaves, put

```math
\mathcal H=\bigoplus_{j=0}^n\mathcal H_j,
\qquad \dim\mathcal H_j=2^j,
\qquad P_j=\text{projection onto depth }j.
```

Let B and E increase depth by exactly one and vanish on leaf columns.
Assume $`\|B\|\le1`$ and

```math
E_j=P_{j+1}EP_j,\qquad
\|E_j\|\le d\,2^{j-n},\qquad 0\le j\lt n,
\qquad 0\le d\lt1.
```

The depth blocks may be arbitrary linear maps; no commutation or nonzero
edge assumption is imposed. Define $`A=B+E`$ and

```math
G_B=(I-B)^{-1}=\sum_{r=0}^nB^r,\qquad X=G_BE.
```

The exact finite identities are

```math
\begin{aligned}
I-A&=(I-B)(I-X),\qquad X^{n+1}=0,\\
G_A&=(I-X)^{-1}G_B,\qquad
G_A-G_B=G_AEG_B.
\end{aligned}
```

For actual tree shifts, each descendant-ancestor entry of $`G_A`$ is
the unique path product. Thus $`O(N)`$ edge data specify the
$`O(Nn)`$ possible entries. These polynomial identities remain valid
at zero edges and require no inversion oracle.

The [coarse schedule](COLLECTIVE_PRECISION_REFINEMENT.md#1-rounded-coins-and-the-actual-coarse-word)
gives one source of the defect assumption. Form A and B from the first
columns of the requested and ideal rounded comparison coins. Both shifts
have normalized outgoing vectors. At depth j, the coin-logarithm bound
$`65\,2^{-q_j}`$, with $`q_j=s+n-j+7`$, bounds the coin difference
by the skew-Hermitian exponential integral. Hence

```math
\|E_j\|\le(65/128)2^{-s}2^{j-n}\le d2^{j-n},
\qquad d=2^{-s}.
```

This comparison concerns ideal logical coins. The complete native coarse
word has dirty and rejected components and is not identified with a
strictly descending tree shift.

## 2. Power and tail bounds

**Graded inverse theorem.** For $`1\le k\le n`$,

```math
\|X^k\|\le a_k:=\frac{d^k}{\prod_{t=1}^k(2^t-1)}.
```

To prove this, expand a word of the product as
$`B^{r_k}E\cdots B^{r_1}E`$, with $`r_t\ge0`$. A nonzero input
at depth j satisfies $`j+k+\sum_t r_t\le n`$. The product of the
E-block bounds is at most

```math
\begin{aligned}
d^k2^{k(j-n)+k(k-1)/2+\sum_{t=1}^{k-1}(k-t)r_t}
&\le d^k2^{-k(k+1)/2-\sum_{t=1}^k t r_t}.
\end{aligned}
```

For this fixed word, different input depths have different output depths,
so its norm is the maximum block norm. Every B factor is contractive.
Triangle inequality and geometric summation therefore give

```math
\|X^k\|\le d^k2^{-k(k+1)/2}
 \prod_{t=1}^k\sum_{r=0}^{\infty}2^{-tr}
=\frac{d^k}{\prod_{t=1}^k(2^t-1)}.
```

The first bounds are $`d,d^2/3,d^3/21,d^4/315`$. The ratios after
the second-order term are at most $`d/7`$, giving

```math
\|(I-X)^{-1}\|
\le1+d+\frac{d^2}{3(1-d/7)}\lt\frac{43}{18}.
```

For $`K_p=\sum_{k=0}^pX^k`$ and $`0\le p\lt n`$, the same ratio
argument proves

```math
\|(I-X)^{-1}-K_p\|
\le\frac{d^{p+1}}{\prod_{t=1}^{p+1}(2^t-1)}
       \frac1{1-d/(2^{p+2}-1)}.
```

The tail is zero at $`p=n`$. For $`d=2^{-s}`$, its sufficient
error exponent is

```math
(p+1)s+\sum_{t=1}^{p+1}\log_2(2^t-1)
 +\log_2\!\left(1-\frac{2^{-s}}{2^{p+2}-1}\right).
```

Using $`2^t-1\ge2^{t-1}`$ and the denominator at least $`2/3`$
gives the simpler certificate

```math
\|(I-X)^{-1}-K_p\|
\le\frac32\,2^{-(p+1)s-p(p+1)/2}.
```

The constant bound is for the coupled inverse, not for $`G_B`$.
Multiplying the tail by $`G_B`$ costs at most another factor $`n+1`$
from its finite sum. The separate
[damping inequality](RESIDUAL_ASSEMBLY.md#finite-tree-resolvents-pack-every-insertion-order)
concerns this distinction.

## 3. Input support and descendant data

Every X increases depth, so

```math
X^k=X^kP_{\le n-k},\qquad
\mathrm{rank}\,X^k\le2^{n-k+1}-1.
```

The tail after p has input support $`P_{\le n-p-1}`$ and rank at
most $`2^{n-p}-1`$. These are computational input supports, with no
moving input basis. Tensoring every map and projector with arbitrary
dirty identity preserves the support statements and multiplies ranks
by the dirty dimension; it supplies no initialized space.

**Rank-one descendant-data witness.** Even differences of normalized
real tree shifts can have a rank-one last term with $`N/2`$ independent
descendant parameters. Set each B column to $`(1,0)^T`$ and each A
column to $`(\cos\theta_v,\sin\theta_v)^T`$ in its child pair. Then

```math
E|v\rangle=(\cos\theta_v-1)|2v\rangle
                  +\sin\theta_v|2v+1\rangle,
\qquad
\|E_j\|=\max_{\mathrm{depth}(v)=j}2|\sin(\theta_v/2)|.
```

Fix $`0\lt d\lt1`$. Choose ancestor angles in $`(0,d2^{j-n})`$, and vary the
$`N/2`$ bottom angles independently in a nonzero open subinterval
of $`(0,d/2)`$. Every defect bound holds and every ancestor difference
coefficient is nonzero. At order n, every inserted positive B power
vanishes by depth, so $`X^n=E^n`$ has only its root input column.
At a bottom parent v, its two output coordinates are

```math
a_v(\cos\theta_v-1,\sin\theta_v),\qquad a_v\ne0,
```

where $`a_v`$ depends only on ancestors. The derivative in each bottom
angle is nonzero and supported on a different child pair. The resulting
rank-one family therefore has local real dimension $`N/2`$.

Input rank alone does not count the independently specified descendant
coefficients or price their output isometry. This is not a T-count or
program-row lower bound for joint compilation: the same parameters may
be shared across stages.

## 4. Constant-normalization scattering

Now let A and B be actual edge-weighted tree shifts. Set
$`K=G_A(I-B)=(I-X)^{-1}`$, write $`y=Kx`$ and $`z=y-x`$, and
take the root incoming z to be zero. The exact local recursion is

```math
y_v=x_v+z_v,\qquad
z_{2v+b}=A_{2v+b,v}z_v+E_{2v+b,v}x_v.
```

The propagator is A, while the injected defect is E; using B as the
propagator would omit higher orders. This realizes all insertion orders
without a separate order register.

For a subtree, let $`[a_v\ M_v]`$ be its stop-output map from incoming
z and external x. The subtree obeys the same graded bound, since its
restricted depth schedule is no larger. Thus $`\|M_v\|\lt43/18`$.
With $`\sigma=3`$, define the finite quadratic maximum

```math
\rho_v=\sup_x\left\{\|a_v+M_vx\|^2-9\|x\|^2\right\}.
```

Eliminating descendant inputs and then the local scalar input gives

```math
\begin{aligned}
\beta_v&=8-\sum_b\rho_{2v+b}|E_{2v+b,v}|^2,\\
\gamma_v&=1+\sum_b\rho_{2v+b}
              \overline{A_{2v+b,v}}E_{2v+b,v},\\
\rho_v&=1+\sum_b\rho_{2v+b}|A_{2v+b,v}|^2
              +\frac{|\gamma_v|^2}{\beta_v}.
\end{aligned}
```

At leaves the sums are empty. The remaining quadratic has cross term
$`2\Re(\overline z\gamma_v x)`$ and maximum at
$`x=\overline{\gamma_v}z/\beta_v`$, proving the conjugations and
recurrence. The Schur complement of $`9I-M_v^\dagger M_v`$ gives

```math
\beta_v\ge9-(43/18)^2>3,\qquad \rho_v\ge1.
```

Indeed that positive matrix is at least $`[9-(43/18)^2]I`$;
the inverse diagonal formula gives the same lower bound for its scalar
Schur complement. Taking x zero gives $`\rho_v\ge\|a_v\|^2\ge1`$.

In semantic output order stop, left continuation, right continuation,
rejection, prescribe

```math
p_v=\frac1{\sqrt{\rho_v}}
\begin{pmatrix}
1\\ \sqrt{\rho_{2v}}A_{2v,v}\\
\sqrt{\rho_{2v+1}}A_{2v+1,v}\\
-\overline{\gamma_v}/\sqrt{\beta_v}
\end{pmatrix},\qquad
q_v=\frac13
\begin{pmatrix}
1\\ \sqrt{\rho_{2v}}E_{2v,v}\\
\sqrt{\rho_{2v+1}}E_{2v+1,v}\\
\sqrt{\beta_v}
\end{pmatrix}.
```

The recurrence gives unit norms, and
$`p_v^\dagger q_v=(\gamma_v-\gamma_v)/(3\sqrt{\rho_v})=0`$.
Complete them to a four-mode unitary. At leaves, deleting continuation
coordinates gives the fixed two-mode determinant-one matrix with columns
$`(\sqrt8,-1)^T/3`$ and $`(1,\sqrt8)^T/3`$.

The incoming continuation amplitude is $`\sqrt{\rho_v}z_v/3`$ and
the marker input is $`x_v`$. The stop output is therefore $`y_v/3`$,
and the outgoing continuations have precisely the required amplitudes.
Using the [continuation and marker slots](RESIDUAL_ASSEMBLY.md#3-exactly-one-signal-flag-supplies-the-modes)
and gathering all stop outputs produces an ideal unitary with accepted
block $`K/3`$. The root continuation is a rejected input mode, zero
only on accepted inputs. Every other input column also has its ordinary
unitary action; no outcome is discarded.

For m vertices, this network has $`2m`$ physical modes. The selected
branch interface on n logical qubits uses the internal heap of height
$`n-1`$, with $`m=N-1`$. Add a fixed two-mode block on the dummy
logical mode to reach $`2N`$ modes: one signal flag, leaving the second
flag as selector. The height-n heap instead needs $`4N-2`$ modes and
consumes both flag copies.

## 5. Native cost and feedback

The ideal scatterers are effectively approximable without real equality
tests. The separated bounds for beta and rho allow certified arithmetic
and positive square roots. At each nonterminal node, search a fixed
dense SU(4) Euler grid until its first two columns are certified within
a positive tolerance of the prescribed frame. Such a grid point exists
because every orthonormal two-frame has an SU(4) completion. Nearby
orthonormal frames can be aligned by a unitary $`O(\text{tolerance})`$
from identity; adjust its determinant on the complementary two-space.
Thus the grid word is near an exact completion, and the accepted identity
holds for every such completion.

A fixed number of two-level SU(2) factors realizes each grid point.
These use the [native bank primitive](COARSE_PREFIX_ENCODER.md#3-the-canonical-su2-native-circuit),
with literal determinant-one phases and actual inverses. The leaf
completion is already explicit. This is a finite selection rule for
generic and singular charts; generic real coefficients are not exact
Clifford+T gates.

Chronological implementation charges a fine scatterer bank at each
depth. As a scaling comparison, tolerances $`O(2^{-L}/n)`$ require
axis precision $`L+\lceil\log_2 n\rceil+O(1)`$, giving

```math
T=O(N+nL),\qquad b=O(L+n).
```

The retained fixed-mode decomposition and borrowed-signal primitives
give this sufficient workspace scaling; it does not certify the literal
additive budget $`b=N+n+7`$ for this completion. Also, K alone is not
the affine S in the selected residual interface: injection, stopping
weights and coarse transport belong in that interface's cost.

**Static feedback identity.** Give every local scatterer separate ports.
After fixed wire routing, express its one-bank map as

```math
w=Dx+Cu,\qquad v=B_0x+A_0u,\qquad u=Pv,
```

where P joins each parent continuation output to its child's incoming
port. Eliminating the internal variables gives exactly

```math
w=\left[D+CP(I-A_0P)^{-1}B_0\right]x.
```

The inverse is finite because the continuation graph is acyclic. With
the same local completions, the accepted block is $`K/3`$. Attaching
an output back to an input in this equation is not a gate in the finite
unitary circuit model. Native chronological elimination is the original
depth-ordered scatterer word; classically evaluating inverse coefficients
does not apply that inverse to an unknown quantum input.

For the selected internal heap, the separated local bank has
$`4(N/2-1)+2(N/2)=3N-4`$ ports, versus $`2N-2`$ external tree
ports. Adding the same two-mode dummy gives $`3N-2`$ versus $`2N`$;
for $`n\ge2`$ the separated realization exceeds the one-flag space.
The chronological circuit reuses its $`N-2`$ internal ports. This
count concerns that separated-port realization, not every completion.

## 6. Common conjugators and hidden prefix states

Accepted-block implementations may choose any rejected-space completion
that obeys the [selected branch promises](RESIDUAL_ASSEMBLY.md#6-what-this-assembly-resolves).
They need not preserve the strict-down filtration as complete unitaries.
Their range isometry and coefficient access nevertheless require pricing.

For one common actual native T, complete-unitary cancellation gives

```math
C_j=T^\dagger M_jT\quad\Longrightarrow\quad
C_m\cdots C_1=T^\dagger(M_m\cdots M_1)T.
```

This includes every dirty sector and uses the actual inverse. Intervening
signal phases $`D_j=e^{i\phi_jP}`$ instead leave transported projectors
$`TD_jT^\dagger=e^{i\phi_jTPT^\dagger}`$ in the middle word. For
the history conjugator $`T=SH`$, these projectors contain moving prefix
vectors. They cannot be priced as fixed flag phases. Without these
projectors, a product of unchanged-address banks is again one such bank;
the [fixed-basis witness](COLLECTIVE_PRECISION_REFINEMENT.md#8-a-fixed-basis-one-bank-obstruction)
limits the resulting ansatz.

There is an explicit prefix-state identity even when the accepted inverse
is trivial. Set E zero and fix normalized A; let $`C_0(A)`$ be the
completed, gathered network on a height-h tree. Its accepted input is
the fixed flag rotation

```math
F=\frac13\begin{pmatrix}1&-\sqrt8\\\sqrt8&1\end{pmatrix}.
```

Unitarity then implies

```math
C_0(A)=(F\otimes I)\mathrm{diag}(I,V_A)
```

for a complete unitary on the rejected sector. Here beta is 8, gamma is
1 and $`\rho_v=(h-\mathrm{depth}(v)+1)9/8`$. If $`a_A(v)`$
is the root-to-v path product, the prescribed continuation columns give

```math
V_A|\mathrm{root}\rangle
=-\frac1{\sqrt{h+1}}\sum_v a_A(v)|v\rangle.
```

Indeed the root-continuation input has stop amplitude
$`a_A(v)/\sqrt{\rho_{\rm root}}`$ and rejection amplitude
$`-a_A(v)/\sqrt{8\rho_{\rm root}}`$. Applying $`F^\dagger`$
gives the formula. Normalization follows from
$`\sum_{\mathrm{depth}(v)=j}|a_A(v)|^2=1`$ at every depth.
The column is independent of how the remaining local columns are completed.

Canceling this full physical baseline introduces a unitary carrying all
requested prefix states, although its accepted K is the identity. Its
native price does not follow from that accepted identity. Substituting
a known coarse B for A changes the formula and requires its own error
and precision bound.

## 7. Precision ledger

At $`s=\lceil N/n\rceil`$, the tail certificate in Section 2 reaches
error exponent N only with order $`\Theta(n)`$, unless the order-n
sum is used exactly. The $`O(n^2)`$ gain from the graded denominators
does not alter the leading balance with $`N=2^n`$. A separate factor
$`n+1`$ for $`G_B`$ costs $`\log_2(n+1)`$ extra error bits.

The input support suggests the following conditional doubling ledger,
with $`k_j=2^j\le O(n)`$:

```math
T_j\le c_{\rm rows}N2^{-k_j}
       +c_{\rm prec}k_js+c_{\rm route}\mathrm{poly}(n).
```

Its row charges sum to at most $`O(N)`$, and
$`s\sum_jk_j=O(N)`$. A fixed polynomial routing charge over
$`O(\log n)`$ stages is also $`O(N)`$, absorbing finitely many
small dimensions into the constant. The rank-one witness shows why
shrinking input support alone does not establish this native ledger.
Independent full-size banks or $`O(N)`$ moving encoders at every
doubling stage retain an $`O(N\log n)`$ charge. A word of p full-size
queries has the ledger $`O(pN+pq)`$, becoming $`O(nN)`$ at
$`p=\Theta(n)`$, $`q=\Theta(N)`$.

The constant inverse bound and ideal scattering construction concern
coupled expressions; they provide neither cheap separate powers nor a
constant-call native feedback operation. The selected-branch compiler
must also include the affine diagonal, occupied-selector control, actual
inverses and arbitrary dirty inputs. These component results leave the
complete-frame endpoint bounds unchanged.

## 8. Verification

The [coupled-tree resolvent checks](../../tests/test_coupled_tree_resolvent.py)
exercise the finite inverse identity,
power and tail bounds, computational input support, the normalized-shift
descendant-data witness, complete complex scattering, and the zero-defect
rejected prefix column. The all-size conclusions follow from the proofs
above; finite matrix checks diagnose their implementation.
