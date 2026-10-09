# Tree transport and the remaining endpoint cost

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Tree generators](SOURCE_REUSE_LIMITS.md#3-tree-generators-compress-the-residual-classically) · [Open endpoint](../../docs/OPEN_PROBLEM.md) · [Grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md)

The residual's linear-size classical description admits an exact sparse
transport formula, a complete logical-unitary realization of its normalized
transport, and a height-independent bound on its weighted pieces. These
establish the operator identities and normalization bounds for one joint-block candidate.
Its precision-dependent
gate cost remains unresolved: the construction below does not improve
the retained $`O(N\ell_*(n))`$ one-clean bound or close the selected
$`a=2,b=N+n+7,L=N`$ endpoint.

The candidate is to implement the tree transport coherently, multiply by
the local residual generators, and synthesize their combined operator with
one precision charge. The analysis identifies the normalized transport
exactly and tests a proposed finite-order correction around the coarse frame.
The [weighted-block continuation](WEIGHTED_TRANSPORT_BLOCK.md) now gives an
explicit one-signal-flag dilation and native synthesis using only borrowed
work, with $`T=O(N+nL)`$, $`G=O(NL)`$, and $`b\ge L+n+7`$.
Its core and borrowed amplification signal return approximately within the
full-operator error. The endpoint
$`O(N\log N)`$ T-count still misses the linear target; an alternative
$`O(L\sqrt N)`$ route retains exact dirty return at its stated width.

## 1. Sparse transport contains every path product

Use the heap-indexed tree and local words $`U_v^A`$, $`A\in\{C,W\}`$,
from the tree-generator lemma. Let $`\mathcal I=\{1,\ldots,N-1\}`$
be the internal nodes. The logical columns are indexed by the root state
$`\ast`$ and the N−1 marker nodes. This is the same common reindexing
used in that lemma, not a free physical register operation.

For an internal child $`u=2v+b`$, define

```math
(\mathcal B_A)_{u,v}=(U_v^A)_{b,0},\qquad
(\mathcal G_A)_{u,v}=(U_v^A)_{b,1},\qquad
(\mathcal G_A)_{1,\ast}=1.
```

All other entries are zero. Thus $`\mathcal B_A`$ is square on the
internal-node space and $`\mathcal G_A`$ maps the N logical columns
to that space. Both are contractions. Distinct parent columns have
disjoint child supports; the two local columns are orthonormal. In particular,

```math
\mathcal B_A^n=0,\qquad
\mathcal B_A^\dagger\mathcal G_A=0,\qquad
\mathcal B_A\mathcal B_A^\dagger+
\mathcal G_A\mathcal G_A^\dagger=I_{\mathcal I}.
```

Put

```math
\mathcal F_A=(I-\mathcal B_A)^{-1}
=\sum_{j=0}^{n-1}\mathcal B_A^j,
\qquad \mathcal P_A=\mathcal F_A\mathcal G_A.
```

Then $`(\mathcal P_A)_{u,\ast}=a_A(u)`$ and
$`(\mathcal P_A)_{u,v}=\beta_A(v,u)`$ for a strict descendant u of v;
the other marker entries vanish. Each matrix multiplication advances one
tree edge, so this is exactly the path-product formula, including singular
angles and complex coarse words.

Let $`\iota`$ embed an internal node into its logical marker column.
Write $`D_h=\mathrm{diag}(h_v)`$, $`D_k=\mathrm{diag}(k_v)`$, and
$`D=\mathrm{diag}(g_1-1,\{d_v-1\}_v)`$. The full residual is

```math
C^\dagger W-I
=D+\iota D_h\mathcal P_W+
\mathcal P_C^\dagger D_k\iota^\dagger.
```

This is an exact operator identity using $`O(N)`$ scalar generators.
It neither invokes W as a subroutine nor supplies a free quantum evaluator
for the inverse or its products.

## 2. The joint transport has a known Gram matrix

Although $`\mathcal F_A`$ can amplify by order n, its product with
$`\mathcal G_A`$ has a smaller, exactly known normalization. If d(v)
denotes the depth of node v, then

```math
\mathcal P_A^\dagger\mathcal P_A
=\mathrm{diag}\!\left(n,\{n-1-d(v)\}_v\right).
```

*Proof.* Split the rows of $`\mathcal P_A`$ by output depth
$`j=0,\ldots,n-1`$. At depth j, the nonzero columns are the root state
and the markers strictly above j. They form the complete depth-j frame,
with all other logical columns zero. Its Gram matrix is the projector onto
those columns. Summing these projectors counts n occurrences of the root
column and $`n-1-d(v)`$ occurrences of marker v. ∎

Consequently $`\|\mathcal P_A\|=\sqrt n`$. Deepest marker columns
are zero; dividing every other column by its displayed norm gives an
isometry. This treats injection and propagation jointly.

For comparison, the standalone inverse obeys

```math
\sqrt{\frac{(n+1)(2n+1)}6}
\le\|\mathcal F_A\|\le n.
```

Indeed, $`\phi_j=\mathcal B_A^j|1\rangle`$ are unit vectors on
distinct depths and hence orthonormal. On their span, the shift advances
$`\phi_j`$ to $`\phi_{j+1}`$ until the last depth. Applying its
inverse to $`n^{-1/2}\sum_j\phi_j`$ gives coefficients
$`(j+1)/\sqrt n`$, proving the lower bound; the finite geometric sum
gives the upper bound. A block encoding of this inverse alone therefore
needs normalization of order n. This does not imply that the weighted
residual above needs that normalization: the injection and h/k factors
must be analyzed together.

### The actual weighted pieces have no height penalty

Define $`\epsilon_v=1-|g_v|^2`$ and
$`\gamma=\max_{v\in\mathcal I}\epsilon_v`$; leaves have
$`\epsilon_v=0`$. The overlap recurrence gives

```math
|h_v|^2=\epsilon_v-
\sum_{b=0}^1|(U_v^W)_{b,0}|^2\epsilon_{2v+b},
```

```math
|k_v|^2=\epsilon_v-
\sum_{b=0}^1|(U_v^C)_{b,0}|^2\epsilon_{2v+b}.
```

These follow by taking the squared norm of the first column and first row,
respectively, of $`K_v=(U_v^C)^\dagger\mathrm{diag}(g_{2v},g_{2v+1})U_v^W`$.
Let $`a_A(v,u)`$ be the first-column path product from v to u, including
$`a_A(v,v)=1`$. Telescoping through the subtree gives

```math
\sum_{u\text{ internal below or equal to }v}
|a_W(v,u)|^2|h_u|^2=\epsilon_v,
```

and the same identity with C and k. Consequently

```math
\|D_h\mathcal P_W\|\le2\sqrt\gamma,
\qquad \|\mathcal P_C^\dagger D_k\|\le2\sqrt\gamma.
```

In particular, both are at most two, independent of n. If all subtree
state discrepancies are at most $`\varepsilon`$, then
$`\gamma\le\varepsilon^2`$ and both bounds are at most
$`2\varepsilon`$.

*Finite-tree proof.* Let $`\mathcal R_A`$ have rows
$`\langle q_v^A|`$ on the leaf space. Then
$`\mathcal P_A=\mathcal R_A A`$, where A denotes the complete frame;
this is an identity for proving norms, not permission to call A in a circuit.
First suppose all first-column branch amplitudes are nonzero. The root
state defines a positive leaf measure
$`\mu(x)=|q_1^A(x)|^2`$ and subtree mass $`\mu(v)=|a_A(v)|^2`$.
For a leaf vector $`\psi`$, set $`f(x)=\psi(x)/q_1^A(x)`$.
Then $`\|f\|_{L^2(\mu)}=\|\psi\|`$ and

```math
\langle q_v^A|\psi\rangle
=a_A(v)\,\langle f\rangle_{v,\mu}.
```

For $`A=W`$, the weights $`\alpha_v=\mu(v)|h_v|^2`$ obey
$`\sum_{u\subseteq v}\alpha_u\le\gamma\mu(v)`$ by the
telescoping identity. Let M be the maximum, over the finite ancestral
chain, of the subtree average of $`|f|`$. Disjoint maximal subtrees
where that average exceeds t imply

```math
\sum_v\alpha_v|\langle f\rangle_{v,\mu}|^2
\le\gamma\|M\|_2^2.
```

On the same maximal subtrees,
$`t\mu(M>t)\le\int_{\{M>t\}}|f|\,d\mu`$.
Integrating in t gives
$`\|M\|_2^2\le2\int|f|M\,d\mu\le2\|f\|_2\|M\|_2`$,
so $`\|M\|_2\le2\|f\|_2`$. This proves the forward bound.
Using C and k proves the reverse bound; adjoints and conjugating the
diagonal weights do not change it. All matrices and norms are continuous
in the local words, so approximating zero branches by nonzero branches
proves the singular cases as well. ∎

This uses the standard finite-tree Carleson embedding argument and the
L2 maximal inequality; see [Lai, arXiv:1411.5408v3, Theorems 1.6 and 1.8](https://arxiv.org/pdf/1411.5408v3).
Those inequalities are inherited. Their application here uses the exact
overlap-defect identities of the two prescribed frames. A bounded norm
does not construct a block encoding or price its gates.

### The retained borrowed compiler already supplies a cheap coarse frame

The coarse frame needed for this weighted representation does not require
a new synthesis theorem. Fix $`0\lt\varepsilon_0\le1/64`$, independently
of n, and put $`L_c=\max\{6,\lceil\log_2(1/\varepsilon_0)\rceil\}`$.
Apply the [borrowed-workspace compiler](../../docs/BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements)
at accuracy $`\varepsilon_0`$, using the endpoint's $`b=N+n+7`$ dirty
qubits and leaving both initialized flags untouched. It gives an actual
native frame C with

```math
T(C)=O\!\left(\frac{NL_c}{n+b}+L_c\sqrt N\right)=O(\sqrt N),
\qquad G(C)=O(NL_c)=O(N).
```

Every external helper is returned exactly on arbitrary inputs, including
reference correlations. Thus the circuit implements $`C\otimes I_b`$
exactly, with C close to W; the helper-return statement is not approximate.
The actual inverse has the same resource bounds.

The tree structure follows from that compiler's construction, not just its
global error theorem. Its local words are $`U_v^C=C_x^2`$, where
$`C_x=XQ_x^\dagger XQ_x`$. Scalar phases cancel literally,
$`\det C_x=1`$, and $`XC_xX=C_x^\dagger`$ makes the inactive-sector
action exactly identity. These complex native words therefore retain the
same addressed tree and prescribed local complement columns used above.

The existing per-depth allowance in its Section 3 is

```math
\epsilon_d=\varepsilon_0 2^{d-n},\qquad
\sum_{j=d}^{n-1}\epsilon_j\lt\varepsilon_0.
```

Telescoping the addressed layers within any subtree bounds that subtree's
complete-frame discrepancy by $`\varepsilon_0`$. In particular, both its
root-state and root-complement vectors differ by at most that amount.
Their overlap definitions give $`|g_v-1|,|d_v-1|\le\varepsilon_0`$,
and normalized-state overlap gives $`1-|g_v|^2\le\varepsilon_0^2`$.
Consequently

```math
\|D\|\le\varepsilon_0,\qquad
\|D_h\mathcal P_W\|,\|\mathcal P_C^\dagger D_k\|
\le2\varepsilon_0.
```

The three pieces thus have a constant total norm budget at most
$`5\varepsilon_0`$. The diagonal estimate uses literal-phase closeness;
gamma alone would not provide it. This is a reuse of the existing compiler
and the weighted norm lemma, not a new endpoint gate bound.

The unconditional C and its inverse are safe even when the two flags are
occupied, because those flags can remain untouched and the dirty helpers
are returned on their full input space. This does not price the history
unitary $`V_C`$, either weighted transport block, or an arbitrary controlled
version of C. In particular, it does not replace the target transport by
coarse transport or resolve the coefficient-filter flags.

### Frobenius-small residuals also fit the endpoint budget

At the endpoint pool $`b=N+n+7`$, fix
$`0\lt\varepsilon_0\le1/64`$ and instead use the same borrowed
compiler with

```math
\eta_c=\frac{\varepsilon_0}{\sqrt N},\qquad
L_c=\max\{6,\lceil\log_2(\sqrt N/\varepsilon_0)\rceil\}=O(n).
```

It needs no initialized work ($`a=0`$), leaving both available clean
flags untouched. Its physical circuit is exactly $`C\otimes I_b`$,
with literal-phase error $`\|C-W\|\le\eta_c`$ and exact dirty
return on arbitrary reference-entangled inputs. The retained estimates
give

```math
T(C)=O(n\sqrt N)=O(N),\qquad G(C)=O(Nn),
```

and unitary invariance and the dimension bound give

```math
\|C^\dagger W-I\|_F=\|W-C\|_F
\le\sqrt N\,\|W-C\|\le\varepsilon_0.
```

This is a specialization of the existing compiler, not a new frontier
bound. Its local-word estimate is $`O(L_c+n-d)`$, so it does not
automatically supply the original $`O(n-d+1)`$ baseline-word promise
or the grouped representation's linear table budget. Nor does a small
Frobenius norm provide a cheap coherent evaluator for the residual, its
logarithm, or its inverse Cayley transform. Those native operations and
their precision and clean-work costs remain to be supplied.

## 3. An explicit complete unitary for normalized transport

The normalized columns can be realized without an added depth register.
Use an n-qubit node basis $`0,1,\ldots,N-1`$, with zero left fixed.
At a nonterminal internal node v, put $`r=n-d(v)\ge2`$,
$`c_r=1/\sqrt r`$, and $`s_r=\sqrt{(r-1)/r}`$. On the ordered
three modes $`(v,2v,2v+1)`$, apply

```math
S_v^A=
\begin{pmatrix}
c_r&0&s_r\\
s_r(U_v^A)_{00}&(U_v^A)_{01}&-c_r(U_v^A)_{00}\\
s_r(U_v^A)_{10}&(U_v^A)_{11}&-c_r(U_v^A)_{10}
\end{pmatrix}.
```

It factors as $`\mathrm{diag}(1,U_v^A)`$ times the real orthogonal
matrix with rows $`(c_r,0,s_r)`$, $`(s_r,0,-c_r)`$, and
$`(0,1,0)`$. Thus every factor is a literal unitary on all its inputs.
Execute the factors shallow to deep; factors at the same depth act on
disjoint modes. Call their product $`V_A`$.

Let J embed the internal-node output into the n-qubit node register.
The designated columns are

```math
V_A|1\rangle=J\mathcal P_A|\ast\rangle/\sqrt n,
```

```math
V_A|2v\rangle
=J\mathcal P_A|v\rangle/\sqrt{n-1-d(v)}
\quad(d(v)\le n-2).
```

The input $`|2v\rangle`$ is untouched until its parent's factor
$`S_v^A`$; its second column launches the prescribed complementary local
vector. Every later propagation contributes $`s_r`$ and every stop
contributes $`c_r`$. The products telescope to the same inverse square
root of the remaining height at every possible stopping depth. The root
input has the same proof. The other input columns provide the complete
unitary extension; they are not projected away or assumed initialized.

For three levels, this is the actual eight-dimensional unitary obtained
chronologically from $`S_1^A`$, then the disjoint $`S_2^A,S_3^A`$.
Its columns at node inputs $`1,2,4,6`$ give the normalized root,
root-marker, left-marker, and right-marker transports. Columns at
$`0,3,5,7`$ supply the remaining completion.

### What this circuit has and has not priced

This is an explicit logical-mode construction, not an elementary
Clifford+T cost theorem. There are $`N/2-1`$ three-mode factors and
$`n-1`$ successive levels. Every $`U_v^W`$ is still a target-dependent
rotation. The fixed depth-dependent mixing angles in $`c_r,s_r`$ also
need synthesis; fixed classical data does not make their quantum gates free.
Mode routing, the marker-to-node input permutation, and actual inverses
must be compiled as well.

The existing layer compiler can price target-dependent addressed stages
separately, but then it retains a precision charge at each level. Its
retained baseline estimate is $`O(N+nL)`$, already available before this
change of representation. The displayed transport identity supplies no
joint cancellation or shared precision implementation that improves it.
There is no inference of an additive T-count lower bound from these
separately scheduled calls.

Nor is $`V_A`$ itself the required half-unitary block for W. Its accepted
columns realize normalized transport. The column norms to restore are
$`\sqrt n`$ and $`\sqrt{n-1-d(v)}`$. Simply sandwiching this unitary
between separate diagonal filters therefore retains a possible
$`\sqrt n`$ normalization factor, even though the fully weighted operator
has bounded norm. Exploiting the correlation between the transport and its
weights is a circuit-construction problem that the norm inequality does not
solve. A proposed full residual compiler
must still combine its h/k filters, diagonal term, and coarse frame with
literal phases, the two-clean lifetime ledger, and a final amplification
error proof. This bounded candidate passes the complete logical-operator
test but does not yet pass the selected linear T-count budget.
In particular, a filter flag that already contains a rejected component
cannot be reused as if initialized when calling a two-clean transport
compiler. Its coherent action on that component must be included.

## 4. Fixed-order coarse corrections miss necessary mixed terms

Could the target transport be replaced by cheap coarse transport plus a
fixed number of local corrections? Put

```math
\Delta\mathcal B=\mathcal B_W-\mathcal B_C,
\qquad \Delta\mathcal G=\mathcal G_W-\mathcal G_C,
\qquad Q=\Delta\mathcal B\mathcal F_C.
```

Each factor of Q advances at least one depth. The exact finite expansion is

```math
\mathcal P_W
=\mathcal F_C\sum_{j=0}^{n-1}Q^j\mathcal G_C
+\mathcal F_C\sum_{j=0}^{n-2}Q^j\Delta\mathcal G.
```

The perturbation-degree-r approximation $`\mathcal P_W^{[r]}`$ uses
upper limits r and $`r-1`$, respectively; its second sum is empty at
$`r=0`$. Degree $`n-1`$ is exact. For
$`\delta=\max_v\|U_v^W-U_v^C\|`$ with $`n\delta\lt1`$,

```math
\|D_h(\mathcal P_W-\mathcal P_W^{[r]})\|
\le\frac{(n+1)(n\delta)^{r+2}}{1-n\delta}.
```

To see this, use $`\|\mathcal F_C\|\le n`$,
$`\|\Delta\mathcal B\|,\|\Delta\mathcal G\|\le\delta`$,
and $`\|\mathcal G_C\|\le1`$ to sum the two geometric tails.
Their combined norm is at most
$`(n+1)(n\delta)^{r+1}/(1-n\delta)`$.
Subtree unitary telescoping gives $`\|D_h\|\le n\delta`$:
$`h_v=\langle e_v^C|q_v^W-q_v^C\rangle`$ and the remaining
subtree has at most n layers. This proves the estimate while keeping
the exact h/k and diagonal generators; it is not an error estimate for
an unspecified approximate quantum evaluator.

### A witness at every height

Take all coarse words to be the identity, an exact native circuit. Put
target rotations of angle $`\theta`$ only on a right spine, ending
at an internal node u of depth d. The exact residual entry is

```math
E_{u,\ast}=(\sin\theta)^{d+1}.
```

Every coarse right-edge first-column coefficient is zero. Reaching u
from the root therefore requires d perturbation factors. A transport
truncation of degree $`r\lt d`$ makes this entry zero, even though
the exact local generator $`h_u=\sin\theta`$ has been retained.

In the physical marker ordering, the two-level witness is
$`E_{3,0}=\sin\theta_0\sin\theta_1`$. At three levels it is
$`E_{7,0}=\sin\theta_0\sin\theta_1\sin\theta_2`$.
Thus the first-order transport correction already misses a concrete
three-level term.

Fix $`p>1`$ and take $`\theta=n^{-p}`$. For every $`r\le n-2`$,
choose $`d=r+1`$. The residual approximation then obeys

```math
\|E-\widehat E_r\|
\ge(\sin n^{-p})^{r+2}
\ge2^{-n(1+p\log_2 n)}
>2^{-2^n}
```

for sufficiently large n. The first inequality is one matrix entry;
the second uses $`\sin\theta\ge\theta/2`$ and $`r+2\le n`$.
Consequently this uniform truncation method needs full transport degree
$`n-1`$ to guarantee the endpoint accuracy from polynomially small
coarse discrepancies. Alternatively, on a fixed omitted path it must
require much finer coarse data: $`\sin\theta\le\eta^{1/(r+2)}`$.
The witness rules out the stated fixed-order substitution, not all
resummations or all two-clean compilers. It does not make source-call
costs additive.

## 5. What the next endpoint construction must change

The compact object to synthesize is the weighted, fully resummed residual
$`D+\iota D_h\mathcal P_W+\mathcal P_C^\dagger D_k\iota^\dagger`$.
Separately normalizing its inverse loses a factor of order n; treating
injection and transport together gives the exact smaller Gram matrix;
including the local overlap defects gives height-independent norms.
The logical history unitary shows that no extra depth register is needed
merely to realize its normalized columns. The unresolved requirement is
one jointly charged Clifford+T implementation, including coefficient
filters, source programming, coarse words, routing, and inverse calls.

A finite correction order cannot replace this task at the stated
high-precision endpoint. Expanding all ancestor entries retains the
old $`\Theta(nN)`$ table representation. Neither observation is a
general impossibility result. The existing grouped compiler remains
the proved resource bound while this candidate's precision cost is open.

The [revision decision](../ROUTE_HISTORY.md#revision-decision-and-next-bounded-pass)
reuses the cheap coarse frame above. The
[weighted-block construction](WEIGHTED_TRANSPORT_BLOCK.md) supplies a
forward component with constant normalization, a complete flag action, and
a native implementation. Its $`O(N\log N)`$ endpoint T-count
does not replace the grouped compiler. The [two-flag assembly](RESIDUAL_ASSEMBLY.md) now incorporates the diagonal
into an affine forward block and combines it with the reverse term, including
actual controls and the complete workspace schedule. The remaining task is
cheaper joint precision synthesis for those compatible selected branches.

### Small defects do not make the selected completion a small correction

The weighted block is stable as an operator, but the particular completion
selected in the [weighted-block construction](WEIGHTED_TRANSPORT_BLOCK.md)
need not approach its zero-defect completion. For an exact witness, fix
$`\alpha=4\varepsilon_0`$, take $`n=1`$, $`C=I`$, and
$`W(t)=R_y(t)`$, with $`-\varepsilon_0/2\lt t\lt0`$. Then

```math
h_1=\sin t,\qquad \rho_1=\sin^2t,\qquad
A(t)=\frac{\sin t}{\alpha}|1\rangle\langle0|.
```

Here A is the accepted block of the selected complete unitary $`U(t)`$.
Its terminal phase is $`\omega_1=-1`$ for negative t, whereas the
zero-defect convention chooses $`\omega_1=1`$ at zero. The root mixing
leaves the logical marker input untouched; the terminal factor therefore
sends that input to opposite signs of the same rejection output. Hence

```math
\|U(t)-U(0)\|=2,
\qquad
\|A(t)-A(0)\|=\frac{|\sin t|}{\alpha}\longrightarrow0.
```

Thus small original-frame defects do not justify a small full-unitary
correction around this zero-defect completion. This is a property of the
selected completion, not a lower bound or an obstruction to other
completions.

A native replacement may choose a different rejected-space action. For
$`F=\iota D_h\mathcal P_W`$, its component interface can require an
explicit unitary Q, its literal inverse, and

```math
\left\|J^\dagger QJ-(F/\alpha)\otimes I_b\right\|\le\delta,
```

where J initializes the allocated signal flags. It need not approximate
the selected $`U(t)`$ on rejected inputs. The retained
[amplification lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
uses an accepted-block estimate and the actual inverse, rather than
proximity to a prescribed completion. For an arbitrary replacement of this
F-only component, controlled calls and charged workspace still need proof;
it does not automatically implement the affine $`S=A+F`$ branch. The
[residual assembly](RESIDUAL_ASSEMBLY.md) already proves complete-frame
composition for its specified S/R interfaces. A replacement meeting those
interfaces can reuse that theorem with its own native resource ledger.

The [finite checks](../../tests/test_tree_transport.py) compare independent
path products and shift formulas, complete residual matrices, Gram
identities, the three-mode unitary completion, and perturbation witnesses.
They use small independent matrix checks, including singular charts;
they do not certify an elementary joint-block compiler or a new lower
bound for unrestricted frame synthesis.

## 6. Small products suggest a Cayley representation

The small complete-frame examples admit a compact skew-Hermitian
representation. It extends to an actual complex native coarse frame and
does not normalize vanishing residual vectors. This is an operator
representation, not a native source program or an improved compiler bound.
The Cayley transform and its connection with low-rank unitary/Hermitian
structure are standard; see Del Corso, Poloni, Robol, and Vandebril,
[Section 7](https://arxiv.org/abs/1811.05854). The recursion below specializes
that elementary matrix tool to the present complete tree residual, without
a priority claim for the transform or low-rank update method.

### The two-level real example

Let a, b, c be the half-angle tangents at the root, left child, and right
child. In logical order 0, 1, 2, 3, the frame applies the root rotation
on the pair (0, 2), followed by the rotations on (0, 1) and (2, 3).
Assume the three half-angle tangents are finite. Its Cayley transform is exactly

```math
K=(W-I)(W+I)^{-1}
=\begin{pmatrix}
0&-b&-a&-ac\\
b&0&-ab&-abc\\
a&ab&0&-c\\
ac&abc&c&0
\end{pmatrix}.
```

The sibling off-diagonal block is the rank-one matrix
$`-a(1,b)^{\mathsf T}(1,c)`$. The full matrix is generally dense and
full rank. Thus the useful feature is a recursive description, not a
small global rank. A three-level example can have all 28 upper-triangular entries nonzero
while retaining these recursive rank-one sibling couplings.

For a real frame with identity coarse frame, write
$`p_v=(I+K_v)e_0`$. If t is the root half-angle tangent, then

```math
K_v=\begin{pmatrix}
K_L&-t p_Lp_R^{\mathsf T}\\
t p_Rp_L^{\mathsf T}&K_R
\end{pmatrix},\qquad
p_v=\begin{pmatrix}p_L\\t p_R\end{pmatrix}.
```

Products of half-angle tangents therefore describe the coupling columns.
Using the full target W in this transform is not a uniformly small-norm
construction: W may approach an eigenvalue minus one. The residual below
uses the established coarse approximation instead.

### The actual complex-coarse residual

Let $`\mathcal R_v=\mathcal C_v^\dagger\mathcal W_v`$ be the
complete subtree residual. Denote its local coarse and target two-mode
words by C and U. The coarse C is the actual phase-calibrated native
word; it need not be real. Let E inject the two existing child-root
modes and put $`\mathsf C=I+E(C-I_2)E^\dagger`$. Then

```math
\begin{aligned}
D&=\mathcal R_L\oplus\mathcal R_R,& F_2&=UC^\dagger,\\
F&=I+E(F_2-I_2)E^\dagger,&
\mathcal R_v&=\mathsf C^\dagger DF\mathsf C.
\end{aligned}
```

Define Cayley transforms where their denominators are invertible:

```math
\begin{aligned}
K_v&=(\mathcal R_v-I)(\mathcal R_v+I)^{-1},\\
A&=K_L\oplus K_R,\qquad
\tau=(F_2-I_2)(F_2+I_2)^{-1}.
\end{aligned}
```

All these generators are skew-Hermitian. Put

```math
B=E^\dagger AE=\mathrm{diag}(\beta_L,\beta_R),\qquad
\beta_j=e_0^\dagger K_je_0,\qquad Z=(I+A)E,
```

```math
H=\tau(I_2+B\tau)^{-1}.
```

The exact recursion is

```math
K_v=\mathsf C^\dagger\bigl(A+ZHZ^\dagger\bigr)\mathsf C.
```

To prove it, set $`T=E\tau E^\dagger`$. Factoring DF plus and minus
identity gives

```math
\mathrm{Cayley}(DF)
=A+(I+A)T(I+AT)^{-1}(I-A).
```

The identity $`T(I+AT)^{-1}=EH E^\dagger`$ reduces the inverse to
two dimensions, and $`I-A=(I+A)^\dagger`$ gives the displayed result.
The push-through identity also yields
$`H^\dagger=-(I+\tau B)^{-1}\tau=-H`$, including singular tau.
No inverse of tau or of a boundary difference is used.

For real rotation words C and U and real child residuals, B is zero. The correction
then uses $`H=\tau=t\begin{pmatrix}0&-1\\1&0\end{pmatrix}`$,
where t is the half-angle tangent of the local target-minus-coarse
rotation. The complex case requires the small two-by-two inverse above;
setting B to zero there would discard genuine imaginary root entries.

### Uniform conditioning and constant local data

Use the [retained coarse approximation](#the-retained-borrowed-compiler-already-supplies-a-cheap-coarse-frame)
with a uniform subtree bound
$`\|\mathcal R_v-I\|\le\epsilon_0\lt1`$ at every subtree, as well
as $`\|F_2-I\|\le\epsilon_0`$ for the local discrepancies. Write
$`\kappa=\epsilon_0/(2-\epsilon_0)\lt1`$. Then

```math
\|K_v\|\le\kappa,\qquad \|B\|,\|\tau\|\le\kappa,
\qquad \|(I+B\tau)^{-1}\|\le\frac1{1-\kappa^2},
```

```math
\|H\|\le\frac{\kappa}{1-\kappa^2},\qquad
I_2\preceq Z^\dagger Z\preceq(1+\kappa^2)I_2.
```

The subtree hypothesis matters: separate child and local error bounds
of $`\epsilon_0`$ alone give at most twice that bound for their parent.
Use the actual uniform subtree promise, not that invalid implication.

Only constant-size local data are needed. If $`u=Ce_0`$, then

```math
\beta_v
=u^\dagger\bigl[B+(I+B)H(I-B)\bigr]u.
```

Thus the beta summaries and H matrices are computed bottom-up with a
constant number of scalar operations per node. For the coupling column,

```math
g=\bigl[I_2+H(I_2-B)\bigr]u,\qquad
p_v=\mathsf C^\dagger
\begin{pmatrix}g_0p_L\\g_1p_R\end{pmatrix}.
```

Leaves have beta zero and p equal to one. Keep these columns recursively
specified rather than expanding one dense vector at every node. The
tree therefore uses $`O(N)`$ scalar/matrix generators, with their
precision and certified classical evaluation separately charged.
After the root coarse conjugation, the sibling off-diagonal block of
each recursively defined $`K_v`$ lies in the column and row spans of
$`\{e_0,p_L\}`$ and $`\{e_0,p_R\}`$, respectively, so its rank is
at most two. These redundant bases also nest. Put

```math
\mathcal B_v=[e_0,p_v],\qquad h=(C^\dagger-I_2)(I_2+B)g.
```

The root-only coarse action gives
$`p_v|_L=h_0e_0+g_0p_L`$ and
$`p_v|_R=h_1e_0+g_1p_R`$. Consequently

```math
\mathcal B_v=
\begin{pmatrix}\mathcal B_L&0\\0&\mathcal B_R\end{pmatrix}
\begin{pmatrix}
1&h_0\\0&g_0\\0&h_1\\0&g_1
\end{pmatrix}.
```

The local difference $`K_v-A`$ is supported in these two child spans:
$`Z`$ has columns $`p_L,p_R`$, while conjugating A by the root
coarse word adds only root vectors and $`Ae_0`$. The nested identity
keeps every ancestor correction in the same spans upon restriction.
Thus every sibling off-diagonal block of the final root matrix has rank
at most two as well. This is a compact classical representation. Its
bases may lose rank at zero defects; no formula orthonormalizes or
inverts them, and no free quantum state preparation is implied.

For skew-Hermitian K and an approximation $`\widehat K`$ that is
also skew-Hermitian, inverse Cayley conversion has the stable bound

```math
\left\|(I+\widehat K)(I-\widehat K)^{-1}
 -(I+K)(I-K)^{-1}\right\|
\le2\|\widehat K-K\|.
```

Indeed the inverse transform is $`2(I-K)^{-1}-I`$, and both
resolvents have norm at most one. This is an operator perturbation
bound, not a gate-count estimate for implementing a resolvent.

### What remains to make this a compiler

The representation contains no target-frame oracle and no additional
logical modes. Its coefficients are classical data. It does not supply
the coherent loading, routing, or transport of p as free operations.
Nor does a rank-two update give a two-qubit gate: its support vectors
are generally spread over the complete child subspaces.

A useful native word must implement the jointly represented generator
and inverse Cayley action with charged controls, queries, and actual
inverses, within the declared clean/dirty allocation. A power-series
or generic block-encoding conversion must charge every use of the
generator; constant residual norm does not make fine-precision
conversion a constant-query operation. Linear classical storage alone
does not prove linear table work or one global precision charge.
The generic compiler frontier is unchanged.

The [Cayley fixtures](../../tests/test_tree_cayley.py) start with the displayed
four-mode product and extend to eight modes. They compare direct complete
matrices with the recursive data, including actual complex native coarse
words, zero defects, chart boundaries, conditioning, and inverse conversion.
They do not emit a Clifford+T implementation of the Cayley generator or
its inverse.

## 7. A native four-mode benchmark

Four real modes admit the standard magic-basis reduction to two SU(2)
factors; see [Vatan and Williams, Section III, Theorem 1](https://arxiv.org/abs/quant-ph/0308006).
The reduction is inherited. The formulas below fix its literal Clifford
convention, apply it to the Cayley data, and charge the retained native
compiler. This is a fixed-size benchmark, not a generic endpoint gain.

### An exact Clifford basis change

Use computational order 00, 01, 10, 11, with qubit 0 the first tensor
factor. Products are in matrix order. Define

```math
M=\mathrm{CNOT}_{0\to1}H_0S_0S_1\mathrm{CZ}_{01}
 =\mathrm{CNOT}_{0\to1}H_0S_0S_1H_1\mathrm{CNOT}_{0\to1}H_1.
```

Its columns are $`(\Phi^+,i\Psi^+,i\Phi^-,\Psi^-)`$, where
$`\Phi^\pm=(|00\rangle\pm|11\rangle)/\sqrt2`$ and
$`\Psi^\pm=(|01\rangle\pm|10\rangle)/\sqrt2`$. This is an
exact matrix equality, including the scalar phase. The displayed word
uses seven Clifford gates; M and its actual inverse cost fourteen in
total, with no T gates or workspace. The magic-basis identity is

```math
M\,\mathrm{SO}(4)\,M^\dagger
=\{A\otimes B:A,B\in\mathrm{SU}(2)\}.
```

The pair of factors is determined up to simultaneously changing both
signs, which preserves its literal product. The Cayley formulas below
choose the pair continuously from the identity.

### Regular factors directly from the Cayley generator

Let K be real skew-symmetric, with $`k_{ij}=K_{ij}`$ for i less
than j. Write $`\boldsymbol\sigma=(X,Y,Z)`$ and define

```math
\mathbf a=\tfrac12(-k_{01}-k_{23},-k_{03}-k_{12},-k_{02}+k_{13}),
\qquad
\mathbf b=\tfrac12(-k_{01}+k_{23},k_{03}-k_{12},-k_{02}-k_{13}).
```

Direct Clifford conjugation gives

```math
MKM^\dagger
=i\bigl[(\mathbf a\cdot\boldsymbol\sigma)\otimes I
 +I\otimes(\mathbf b\cdot\boldsymbol\sigma)\bigr].
```

Put $`r=\|\mathbf a\|`$, $`s=\|\mathbf b\|`$, and use the
positive normalization

```math
d=\sqrt{(1+(r+s)^2)(1+(r-s)^2)}\ge1.
```

Then

```math
(I+K)(I-K)^{-1}=M^\dagger(A\otimes B)M,
```

```math
A=\frac{(1-r^2+s^2)I+2i\mathbf a\cdot\boldsymbol\sigma}{d},
\qquad
B=\frac{(1+r^2-s^2)I+2i\mathbf b\cdot\boldsymbol\sigma}{d}.
```

For example, $`(1-r^2+s^2)^2+4r^2=d^2`$, so A is unitary
with determinant one; the same holds for B. On simultaneous eigenvectors
of the two commuting Pauli factors their product equals
$`(1+i(er+fs))/(1-i(er+fs))`$, for $`e,f\in\{1,-1\}`$.
This proves the inverse Cayley identity. Neither factor divides by r or
s, so zero defects and rank changes introduce no singular direction
normalization. At K zero both factors are exactly identity.

For the four-mode example in Section 6, rename the root, left, and right
half-angle tangents x, y, z. The formulas simplify to

```math
2\mathbf a=(y+z,x(y+z),x(1-yz)),\qquad
2\mathbf b=(y-z,x(y-z),x(1+yz)),
```

```math
r^2-s^2=yz,\qquad d=\sqrt{(1+x^2)(1+y^2)(1+z^2)}.
```

These recover the complete four-mode matrix, including its complement
columns. Outside this Cayley chart, an explicitly specified real plane
rotation word can instead be factored one rotation at a time and its two
SU(2) products evaluated classically. For instance,
$`|1\rangle\langle0|-|0\rangle\langle1|`$ maps to
$`i(X\otimes I+I\otimes X)/2`$. The linear conjugation formula
gives the other five plane generators. This does not assume a global
Cayley chart or discard a sign when an eigenvalue crosses minus one.

### An addressed primitive with a borrowed logical spectator

Let k address bits remain unchanged, put $`S=2^k`$ and $`n=k+2`$,
and supply factor coordinates admitting certified classical evaluation
for each real four-mode row
$`MW_xM^\dagger=A_x\otimes B_x`$. The Cayley formulas supply
such factors from finite real skew-symmetric K tables; they do not
cover a row W with an eigenvalue minus one. The complete target
is M conjugating two one-target multiplexors with that same address.

Apply the retained [one-target compiler](../../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
to each factor at $`L'=L+1`$. Its sufficient dirty reservation is
$`L'+k+9=L+k+10`$. External work of size
$`b=L+n+7=L+k+9`$, together with the other logical target,
meets this reservation: while compiling A use target 1 as the additional
dirty helper, and while compiling B use target 0. The helper returns
exactly. Arbitrary logical and reference correlations are allowed by the
primitive's dirty-input contract. Both factors reuse one initialized
flag and the external core; the endpoint's second clean flag can remain
untouched.

Let $`\mathcal A,\mathcal B`$ denote the ideal addressed factors,
and J the common initialized-flag embedding. Although the helper
assignment changes, each circuit has the same complete-isometry promise.
Therefore

```math
\begin{aligned}
\|V_BV_AJ-J\mathcal B\mathcal A\|
&\le\|V_AJ-J\mathcal A\|
  +\|V_BJ-J\mathcal B\|\\
&\le2\cdot2^{-(L+1)}=2^{-L}.
\end{aligned}
```

This norm includes every logical and dirty input, references, and actual
flag leakage; no intermediate reset is used. The exact M boundaries
preserve the bound. For $`L\ge6`$ the resulting addressed primitive has

```math
a=1,\qquad b\ge L+n+7,\qquad
T=O(S+L),\qquad G=O(SL),
```

plus the fourteen displayed Clifford boundary gates. Certified Euler
preprocessing evaluates the exact computable factor functions to its
own allocated tolerance. If K is instead rounded separately, its
inverse Cayley error of at most $`2\|\widehat K-K\|`$ must also
be allocated; the displayed budget does not include that extra error
for free. Classical preprocessing cost is not asserted to be uniform.
The source-call count is constant, not a single literal preparation and
unpreparation.

### Real residuals and the complex-coarse boundary

For real coarse $`C\in\mathrm{SO}(4)`$, the residual
$`C^\dagger W`$ remains in this family. The actual native coarse
word may instead be complex, so its residual generally does not admit
this real Cayley split. A valid fallback is the literal word

```math
C^\dagger M^\dagger V_BV_AM
```

approximating the residual $`C^\dagger W`$ of the real target W.
Charge $`T(C)`$ and $`G(C)`$ in addition
to the displayed synthesis cost. The retained unconditional coarse
frame returns its dirty helpers exactly on their full input space, so
its actual inverse preserves the error even on leaked inputs. At fixed
coarse accuracy in the selected endpoint pool $`b=N+n+7`$, its
separately proved costs are $`O(\sqrt N)`$ T
gates and $`O(N)`$ Clifford gates. This does not price a controlled
coarse word or give SO(4) closure for the complex residual itself.

The fixed-four-mode result is a corollary of the standard magic basis
and the retained compiler. Extending it to two four-mode children and
their root coupling still requires a native composition rule and a
source-call recurrence. Fixed-size $`O(L)`$ synthesis alone gives no
improvement to the generic $`O(N\ell_*(n))`$ bound.

## 8. The eight-mode root retains a controlled coupling

Use the Clifford M with columns
$`[\Phi^+,i\Psi^+,i\Phi^-,\Psi^-]`$, equivalently
$`M=\mathrm{CNOT}_{0\to1}H_0(S_0S_1\mathrm{CZ}_{01})`$.
The four-mode factorizing direction is $`MWM^\dagger`$.
Applying the same M to both children of an eight-mode frame therefore
uses $`I\otimes M`$. The root rotates the existing pair (0, 4), so
its transformed complete word is

```math
U_\theta=\exp(-i\theta Y\otimes P),\qquad
P=M|00\rangle\langle00|M^\dagger
=|\Phi^+\rangle\langle\Phi^+|
=\frac{II+XX-YY+ZZ}{4}.
```

Consequently it has the exact decomposition

```math
U_\theta
=e^{-i\theta YII/4}e^{-i\theta YXX/4}
 e^{+i\theta YYY/4}e^{-i\theta YZZ/4}.
```

The four Pauli strings commute. Each factor is a Clifford conjugate of
a single-qubit rotation, whose precision must be synthesized and charged.
This supplies a constant-size construction, not a proof that four
separate calls are necessary. At the special angle $`\theta=\pi/2`$,
three T words and one T-dagger word suffice, together with Clifford gates:
their scalar $`e^{i\pi/4}`$ is removed by
$`(S^\dagger H)^3=e^{-i\pi/4}I`$. The finite native fixture retains
that correction and every logical input column.

The root cannot instead be absorbed into independent prefix/child
factors. For $`0\le\theta\le\pi/2`$, its exact distance is

```math
\inf_{A\in U(2),\,B\in U(4)}
\|U_\theta-A\otimes B\|=2\sin(\theta/4).
```

For the lower bound, choose a unit vector chi orthogonal to Phi plus and
apply the word to the product input

```math
|0\rangle\otimes
\frac{|\Phi^+\rangle+|\chi\rangle}{\sqrt2}.
```

Use $`R_y(t)=e^{-itY}`$. The output is the equally weighted sum of
$`R_y(\theta)|0\rangle\otimes|\Phi^+\rangle`$ and
$`|0\rangle\otimes|\chi\rangle`$. Its Schmidt coefficients across
the prefix/child cut are $`\cos(\theta/2)`$ and
$`\sin(\theta/2)`$. Every product-gate output on this input is a
normalized product vector, so its distance is at least
$`\sqrt{2-2\cos(\theta/2)}=2\sin(\theta/4)`$.
For attainment, take $`A=R_y(\theta/2)`$ and $`B=I_4`$.
On both the P and its orthogonal sector the difference is a half-angle
rotation, with precisely that operator norm.

The transformed root and child families also do not commute. The
left-child root generator, originally the pair (0, 2), becomes

```math
H_L=-\frac14(I+Z)\otimes(ZI+IZ),\qquad
\|[Y\otimes P,H_L]\|=1.
```

Before Clifford conjugation their commutator is minus i times the
Hermitian generator of the pair (2, 4), which proves the norm identity.
Thus independent four-mode factorizations do not make the eight-mode
ordered word a product of globally commuting families.

These are restrictions on this particular fixed product-factor extension.
Controlled factors, a jointly programmed word, and other source boundaries
remain allowed. An eight-mode instance still has an $`O(L)`$ synthesis;
neither noncommutation nor the displayed distance makes precision costs
additive across a variable number of groups.

The [eight-mode checks](../../tests/test_eight_mode_coupling.py) verify these
complete matrices, the distance witness, and the literal four-T special
case. General angles remain charged synthesis targets.

## 9. Inverse Cayley recovers the original local scattering word

The fixed-size Cayley update has a well-conditioned inverse, but applying
Woodbury does not remove the local target word. Keep the notation of
Section 6 and write

```math
\overline K=A+ZHZ^\dagger,\qquad Q=(I-A)^{-1},\qquad
D=(I+A)(I-A)^{-1}.
```

Since $`QZ=DE`$ and $`Z^\dagger Q=E^\dagger`$, the Woodbury
identity, in a form that never inverts H, gives

```math
(I-\overline K)^{-1}
=Q+DE\,H[I_2-(I_2+B)H]^{-1}E^\dagger.
```

The remaining two-dimensional inverse cancels explicitly:

```math
\begin{aligned}
I_2-(I_2+B)H&=(I_2-\tau)(I_2+B\tau)^{-1},\\
H[I_2-(I_2+B)H]^{-1}
&=\tau(I_2-\tau)^{-1}=\frac{F_2-I_2}{2}.
\end{aligned}
```

Moreover
$`\|[I_2-(I_2+B)H]^{-1}\|\le1+\kappa^2`$, because
$`\|(I_2-\tau)^{-1}\|\le1`$. Thus inverse Cayley returns

```math
2(I-\overline K)^{-1}-I
=D[I+E(F_2-I_2)E^\dagger]=DF.
```

Including the coarse conjugation gives exactly
$`\mathcal R_v=\mathsf C^\dagger D\mathsf U`$, where
$`\mathsf U=I+E(U-I_2)E^\dagger`$. The implementation uses one
coherent child-residual call, the local target word, and the coarse inverse.
It does not introduce a new cheap implementation of that target word.

There is also a complete one-signal scattering boundary. For any unitary R
whose Cayley chart exists, put $`S=(I-K)^{-1}=(I+R)/2`$. Then

```math
\Phi(R)=
\begin{pmatrix}S&S-I\\S-I&S\end{pmatrix}
=(H_f\otimes I)\mathrm{diag}(R,I)(H_f\otimes I).
```

This specifies both occupied signal sectors and is unitary on the whole
space. It obeys $`\Phi(R_1)\Phi(R_2)=\Phi(R_1R_2)`$, hence

```math
\Phi(\mathcal R_v)
=\Phi(\mathsf C^\dagger)\Phi(D)\Phi(\mathsf U).
```

The common signal therefore composes through arbitrarily many levels with
literal inverse words. However, the displayed completion already contains
controlled R. It is not a circuit for the resolvent supplied independently
of R. Compiling its conditional local target and coarse words with the
retained layerwise primitives reproduces $`O(N+nL)`$ T count and
$`O(NL)`$ Clifford count; no joint precision program follows from this
completion alone.

Approximation must retain the full dirty-space contract. If actual unitary
child and local circuits approximate $`D\otimes I_b`$ and
$`\mathsf U\otimes I_b`$ with errors $`\delta_D,\delta_U`$,
their composed word has error at most $`\delta_D+\delta_U`$,
including dirty disturbance and reference correlations. With the child and
coarse factors exact, the error from the local word is exactly
$`\delta_U`$. Neither a small Cayley generator nor a small rank-two
correction attenuates that error.

The fixed two-by-two matrices describe the exact logical generators. A
physical approximate child circuit that disturbs dirty work can have
operator-valued dirty entries in its Cayley compression. One must not
substitute independently rounded scalar H data and claim that the exact
Woodbury identity still holds for that arbitrary approximate completion.
The actual composed unitary word and its norm bound remain valid without
such a claim.

The [Cayley checks](../../tests/test_tree_cayley.py) verify Woodbury cancellation
and all scattering ports on four/eight-mode complex-coarse fixtures,
including a native child that disturbs borrowed work. The
[four-mode native checks](../../tests/test_native_cayley.py) separately test
the magic-basis factors and the borrowed-spectator schedule. These small
fixtures check identities and interfaces; the resource bounds rely on
the proofs and retained compiler contracts.

## 10. A shared source body for changing targets

The original tree basis admits a joint native word for root and child
layers at one common source width, with one final parity correction for
arbitrary borrowed signal inputs. Its declared source-call count decreases,
but simplifying fixed masks gives the comparison words the same leading
precision cost. That cost still grows with the number of layers.

### Keep the fixed scalar word independent of the logical route

Use the [paired source and scalar words](../../docs/ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks)
at one fixed q. With f the fixed mask, define

```math
A=\mathcal S_f=\tfrac12I+X_aD_f,\qquad
D_f^\dagger=-D_f,\qquad D_f^2=-\tfrac34I.
```

A is unconditional and acts only on the signal a and precision core.
Consequently $`A^3=-I`$ and $`Z_aAZ_a=A^\dagger`$ exactly.
For an unchanged predicate h and prefix-programmed mask mu, let
$`\mathcal S_{\mu,h}`$ be the retained conditioned scalar word.
It is exactly identity on h false, including occupied signal inputs.
For a Hermitian logical Pauli P disjoint from the address, predicate,
signal, and core, set

```math
B=C_a(P)S_a\mathcal S_{\mu,h}S_a^\dagger C_a(P),\qquad
Q=A^\dagger BA,
```

where $`S_a=\mathrm{diag}(1,i)`$ and $`C_a(P)`$ is the literal
Clifford controlled Pauli. For a tree layer choose $`P=Y_t`$ on
its current target. That target is not added to the coefficient-table
address. Query selectors and the predicate helper return exactly before
the surrounding routing gates resume.

On an active address row the scalar moments satisfy
$`\{D_f,D_\mu\}=2pI`$, with $`p=s/2-\tau`$. Expansion gives

```math
B=sI+Y_aPD_\mu,\qquad K=D_\mu+2pD_f,
```

```math
Q=sI-ipZ_aP+Y_aPK,\qquad
K^\dagger=-K,\qquad K^\dagger K=(1-s^2-p^2)I.
```

For the last identity use $`D_\mu^2=-(1-s^2)I`$ and
$`D_f^2=-3I/4`$. In particular,
$`\langle0|Q|0\rangle=sI-ipP`$. On h false, B and Q are
both exactly identity; A need not be conditioned. The existing digit
rule programs $`s\approx\beta\cos\theta`$ and
$`p\approx\beta\sin\theta`$, with
$`\beta=\sin(\pi/10)`$. The same coefficient-error and
five-call amplification proof gives initialized-column error below
$`30\,2^{-q}`$ for $`e^{-i\theta P}`$.

### The full signal error has the same bound

Fix an active address row and diagonalize P and K, with eigenvalues
$`\rho\in\{1,-1\}`$ and $`i\lambda`$, respectively. Each
signal block is

```math
Q_{\rho,\lambda}=sI+i\rho(-pZ+\lambda Y)\in\mathrm{SU}(2).
```

The five-call amplified block remains in SU(2): its four literal Z
reflections contribute determinant one in total. Its comparison block
$`e^{-i\rho\theta Z}`$ also lies in SU(2). For any two SU(2)
matrices U and V,

```math
(U-V)^\dagger(U-V)=[2-\mathrm{Tr}(U^\dagger V)]I.
```

Thus the full block norm equals its signal-zero-column norm. Taking the
maximum over all sectors proves the complete-operator estimate

```math
\left\|\mathrm{Amp}(Q)
 -\exp(-i\Pi_h\theta_x Z_aP)\otimes I_{\rm dirty}\right\|
\lt30\,2^{-q}.
```

Here $`\Pi_h`$ is the predicate projector and $`\theta_x`$ the
addressed angle. The inactive sector has zero error. The diagonalization
is only a proof device, with no emitted basis preparation. The bound
includes arbitrary signal, dirty, and reference inputs. Its signal-one
target has the opposite angle; it is not the same scalar action on both
signal sectors.

### Expand and cancel the actual amplification words

Write $`Z=Z_a`$. The middle word obeys $`ZBZ=B^\dagger`$ on
the complete space. The physical amplification word is

```math
\mathrm{Amp}(Q)=QZQ^\dagger ZQZQ^\dagger ZQ.
```

Since $`AZA^\dagger=-A^\dagger Z`$, substitution gives

```math
\mathrm{Amp}(A^\dagger BA)=A^\dagger F_A(B)A,\qquad
F_A(B)=(BA^\dagger BA)^2B
=BA^\dagger BAB A^\dagger BAB.
```

The four minus signs cancel literally. The body contains five B words,
two A words, and two actual A inverses. For g ordered stages,

```math
\mathrm{Amp}(Q_g)\cdots\mathrm{Amp}(Q_1)
=A^\dagger F_A(B_g)\cdots F_A(B_1)A.
```

This is a full-unitary identity, not multiplication of accepted blocks.
Targets, prefix addresses, and suffix predicates may change because A
does not depend on them. Logical-only Clifford gates also commute through
A in this identity. The common q and fixed source are essential; a
change of width instead requires an explicitly charged bridge.

Each scalar A or B has three source calls. Before any further native
cancellation, the displayed source ledger is

```math
45g\quad\longrightarrow\quad27g+6.
```

This source-call ledger is not an optimized T count. Let U be the paired
loader, with $`2q`$ T gates. Each routed scalar is $`UR_jU^\dagger`$:
its exterior flag/logical Cliffords commute with U, while its two mask
transforms retain their internal loader pairs. Canceling adjacent exterior
pairs across the $`9g+2`$ scalar calls reduces the loader T count to
$`(72g+20)q`$, versus $`(120g+4)q`$ for the expanded baseline's
15g scalar calls. These give $`236q`$ and $`308q`$ at three and four
stages before simplifying the fixed mask.

That simplification matters for a fair comparison. Factor
$`U=G_{\rm tail}U_{\rm seed}`$, where the seed has two T gates.
The fixed $`P_f`$ commutes with every tail plane: its signs are
constant on each of the odd and even tails. Hence

```math
U^\dagger P_fU=U_{\rm seed}^\dagger P_fU_{\rm seed}.
```

Each such transform needs at most four T gates, so a hoisted fixed-A
interior needs at most eight, independently of q. The five programmable
B words per stage still retain $`8q`$ loader T gates apiece. This gives

```math
T_{\rm loader+fixed}\le
\begin{cases}
(40g+4)q+8(4g+2),&\text{shared body},\\
(40g+4)q+80g,&\text{expanded baseline}.
\end{cases}
```

Both therefore have the same leading precision coefficient. These are
declared-word bounds, not minima. Each programmable mask still appears
ten times per stage, including its actual inverses; its native queries,
source-center predicate controls, and returned helpers remain separately
charged. Fixed-mask Clifford gates are charged as well.

### Return the borrowed signal with two parity boundaries

For the original r-level real tree, let $`\Pi=\prod_{j=0}^{r-1}Z_{t_j}`$
on its r local logical bits, excluding any unchanged external prefix.
It commutes with every address and suffix predicate and reverses each
$`Y_{t_j}`$. If W is the ordered ideal tree word, its complete
signal-controlled product is

```math
|0\rangle\langle0|_a\otimes W
+|1\rangle\langle1|_a\otimes\Pi W\Pi.
```

The opposite-angle factors keep their original chronological order;
this branch is generally neither $`W^\dagger`$ nor $`W^*=W`$.
Put $`R=C_a(\Pi)`$, a product of r CZ gates. Then the complete
native word

```math
V=RA^\dagger F_A(B_r)\cdots F_A(B_1)AR
```

satisfies

```math
\|V-W\otimes I_{\rm work}\|\lt30r\,2^{-q}.
```

The signal is part of the arbitrary returned work. All logical and dirty
reference correlations are included; no intermediate reset is used.
Selectors and the predicate helper return exactly. The two parity
boundaries cost $`2r`$ Clifford CZ gates. This return proof applies
to the original tree layers. Extra interleaved logical gates would also
need to commute with Pi for this same proof; arbitrary logical Clifford
interleaving only preserves the common-A identity above. In particular,
signal-dependent lifts of complex magic-basis gates cannot be inserted
with an unchanged source ledger without their own proof and charges.

### Explicit eight- and sixteen-mode allocation

Take $`r\in\{3,4\}`$, k unchanged external prefix bits,
$`S=2^k`$, and $`n=k+r`$ total logical bits. There are g=r depth
layers. At local depth d, the table has k+d address bits and the suffix
predicate has $`r-d-1\le3`$ literals. Choose the common
$`q=L+7`$, with $`L\ge6`$. Then

```math
30r\,2^{-q}=\frac{30r}{128}\,2^{-L}\le\eta.
```

The maximum simultaneous borrowed reservation is

```math
(q+1)_{\rm core}+(k+r-1)_{\rm selectors}
+1_{\rm helper}+1_{\rm signal}=L+n+9.
```

At $`a=2,b=L+n+7`$, put the core and selectors in the b arbitrary
dirty wires and use the two available clean wires as signal and helper.
Their initial values are not used. The retained
[predicate construction](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger)
returns one separate arbitrary helper exactly, including when a source
center has the additional signal control. Its bounded predicates cost
$`O(p^2+1)`$ Toffolis. Controlled Y and the signal S conjugations
are Clifford and require no helper. This reservation borrows no logical
target.

The table sizes sum to $`S(2^r-1)=O(S)`$ at these fixed depths.
Thus $`T=O(S+L)`$ and $`G=O(SL)`$, with declared source counts
87 instead of 135 for eight modes, and 114 instead of 180 for sixteen.
These are source-appearance comparisons; the fixed-mask refinement above
aligns their leading precision cost. No generic endpoint gain follows.

For variable r, the same error sum asks for
$`q=L+\lceil\log_2(30r)\rceil`$. The unmodified reservation
then exceeds $`b+2`$ by $`\lceil\log_2(30r)\rceil-7`$ when
positive. This diagnoses this schedule, not a workspace lower bound.
Its source recurrence also remains $`27r+6`$, giving
$`O(r(L+\log r))`$ loader cost even with enough work. A joint program
with one global precision charge remains unproved.

The [joint native checks](../../tests/test_joint_source_body.py) verify
literal phases, changing-target composition, occupied signal ports, and
borrowed work in small complete matrices. Their $`q=2`$ words are
finite diagnostics; their deliberately coarse approximation errors do
not certify the stated L-bit resource allocation.

## 11. A packed Hopf scattering step and its boundary transfer

All local Hopf rotations fit into one addressed real-rotation table on
one additional mode qubit. Fixed, charged port permutations turn that
table into a unitary scattering step. The complete prescribed frame is
its boundary transfer matrix, however, rather than the step itself.
The exact identity below identifies that distinction and tests the
unchanged-coin route without assuming a target-dependent eigenbasis.

### Physical ports and the exact transfer identity

Use $`N=2^n`$ and physical labels $`(a,x)`$, with one mode bit a
and the original n-bit logical label x. The external ports $`a=0`$
keep the prescribed root/marker inputs: x zero is the root, and

```math
v=2^d+r,\qquad \lambda(v)=(2r+1)2^{n-d-1}
```

is node v's marker. Internal ports are $`(1,v)`$ for
$`2\le v\lt N`$: incoming continuation amplitudes at nonroot internal
nodes. The two remaining ports $`(1,0),(1,1)`$ are dummies.
Coin coordinates are $`(v,p)`$, with n node bits followed by one
local port bit. Define the forward input permutation P by

| Physical input | Coin input |
| --- | --- |
| $`(0,0)`$ | $`(1,0)`$ |
| $`(0,\lambda(v))`$, $`1\le v\lt N`$ | $`(v,1)`$ |
| $`(1,v)`$, $`2\le v\lt N`$ | $`(v,0)`$ |
| $`(1,0),(1,1)`$ | $`(0,0),(0,1)`$, respectively |

These are bijections of all $`2N`$ modes. Let

```math
\mathscr M=(-I_2)\oplus\bigoplus_{v=1}^{N-1}R_y(\theta_v).
```

The output permutation O sends a coin output $`(v,p)`$ to

```math
\left(1-\left\lfloor\frac{2v}{N}\right\rfloor,
                    (2v+p)\bmod N\right).
```

It is simply X on the highest bit of the coin word: children below N
are internal continuations, and children at least N are external leaves.
In the external/internal physical decomposition, write

```math
\mathscr S=O\mathscr M P
=\begin{pmatrix}D&C_{\rm out}\\ B&A\end{pmatrix}.
```

On genuine internal ports, $`A_0`$ advances one tree level, so
$`A_0^{n-1}=0`$ for $`n\ge2`$. Padding gives
$`A=(-I_2)\oplus A_0`$, with B and $`C_{\rm out}`$ zero on
the dummy coordinates. Thus A itself is not nilpotent, but $`I-A`$
is invertible, including its dummy block $`2I_2`$. Exactly,

```math
W=D+C_{\rm out}(I-A)^{-1}B
 =D+C_{\rm out}\left(\sum_{j=0}^{n-2}A_0^j\right)B,
```

where the last expression suppresses the uncoupled dummy coordinates.
For $`n=1`$, the genuine internal space is empty and $`W=D`$.

To prove the identity, feed the internal output back to its identically
labelled input. For external input x, its unique continuation vector is
$`u=Au+Bx`$, and its external output is $`y=Dx+C_{\rm out}u`$.
At node v these equations are precisely the prescribed local rotation of
the incoming amplitude and marker amplitude into its two children.
Solving downward therefore gives every column of W, with its literal
phase and marker convention. This algebraic feedback relation is not a
physical circuit instruction or an accepted-block implementation of W.
The elimination is standard block linear algebra; the local statement is
the prescribed Hopf port assignment and its charged step implementation.

### The step has a charged native implementation

The input permutation also has a short structured implementation. First
cyclically move a to the last position and flip it, giving $`(x,1-a)`$.
On port one apply the n-bit permutation

```math
H(0)=0,\qquad H(p\,1\,0^t)=0^t\,1\,p.
```

Then swap the two coin basis states $`(0,1)`$ and $`(1,0)`$.
To implement H, reverse all n bits, then reverse the suffix below the
highest one. Unroll the latter operation over its n possible highest-one
positions. Each sector has fixed prefix controls and at most n swaps;
its controls are disjoint from the swapped bits. The retained
[borrowed multi-control construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
implements each controlled swap with $`O(n^2)`$ native gates and one
returned arbitrary helper. The final two-state swap has the same
$`O(n^2)`$ bound by affine conjugation of a multi-controlled X.
Thus P costs $`O(n^4)=O(N)`$ T and Clifford gates; O is Clifford.
The helper can reuse a core wire before or after the complete coin call.

For one coin call at error $`2^{-L}`$, use the
[full-operator real-table primitive](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations).
For $`n\ge2`$, take source width $`q=L+6`$ and split two of the
n unchanged address bits into four invariant sectors. There are
$`k=n-2`$ free address bits and two predicate literals, so each sector
uses

```math
(q+1)_{\rm core}+k_{\rm selectors}
 +1_{\rm helper}+1_{\rm signal}=L+n+7
```

arbitrary dirty wires. The sector errors take a maximum, below
$`43\,2^{-q}=(43/64)2^{-L}`$. Their tables and source costs sum to
$`O(N+L)`$ T gates and $`O(NL)`$ Clifford gates. The fixed-size
$`n=1`$ case uses the retained direct borrowed-sector native words.
All port permutations are exact. Consequently one scattering step has
the full-operator guarantee and costs

```math
\|\widehat{\mathscr S}-\mathscr S\otimes I_b\|\le2^{-L},
\qquad T=O(N+L),\qquad G=O(NL),\qquad b\ge L+n+7.
```

The coin synthesis needs no initialized work. The extra mode qubit may
be initialized to select external inputs, using only one of the two
available clean wires; the bound also covers its occupied internal port.
Core and signal return are included in the approximation, and query
helpers return exactly. This prices one step, not feedback or a sequence
of finer-accuracy calls within the same reservation.

### A fixed number of unchanged-coin calls does not give the frame

Consider a coherent word making k calls to the ideal
$`\mathscr S(\theta)`$, its adjoint, or their controlled versions,
interleaved with arbitrary parameter-independent unitaries. Fixed dense
basis changes and additional initialized or dirty work are allowed.
The [original-angle query theorem](STRUCTURAL_COMPILATION_LIMITS.md#5-original-angle-parallel-coins-need-at-least-n-queries)
proves that $`k\lt n`$ has worst-case literal error at least one on
signed angle tuples, by restricting a path to Boolean parity. Here the
remaining question is the narrower canonical interval $`[0,\pi/2]`$.

Set all angles on the all-right path $`v=1,3,\ldots,N-1`$ to
one variable theta and every other angle to zero. The target entry is

```math
\langle N-1|W(\theta)|0\rangle=\sin^n\theta.
```

Each oracle entry has Laurent degree at most one in
$`z=e^{i\theta}`$. Therefore any selected complete output amplitude
$`p_k(\theta)`$ has Laurent degree at most k. For $`k\lt n`$, it
cannot equal the target on any open angle interval.

This is the elementary degree-growth argument underlying the
[quantum query polynomial method](https://homepages.cwi.nl/~rdewolf/publ/qc/polynomials.pdf),
Lemma 4.1, applied here to Laurent entries rather than Boolean variables.
The local content is the Hopf path amplitude and the specified coin model.

To quantify the error on the canonical interval, put $`d=2n`$. The polynomial

```math
P(z)=z^n\bigl(\sin^n\theta-p_k(\theta)\bigr)
```

has degree d and leading coefficient of magnitude $`2^{-n}`$.
If its modulus is at most epsilon on that arc, interpolate at
$`z_j=e^{i\pi j/(2d)}`$, $`0\le j\le d`$. Since
$`|z_j-z_l|\ge|j-l|/d`$, the leading coefficient is bounded by

```math
\varepsilon d^d\sum_{j=0}^d\frac1{j!(d-j)!}
=\varepsilon\frac{(2d)^d}{d!}.
```

The elementary estimate $`d!\ge(d/e)^d`$ follows by integrating
log x below its increasing sum. Hence

```math
\varepsilon\ge2^{-n}\frac{d!}{(2d)^d}
\ge(8e^2)^{-n}>2^{-6n}.
```

At the selected $`L=N`$ endpoint this excludes $`k\lt n`$ for
$`n\ge5`$. A complete-isometry error bound controls this amplitude
by choosing one logical input and one fixed borrowed input, so extra
work and reference correlations do not evade the conclusion. Replacing
each ideal query by a full-operator delta-approximation changes the word
by at most $`k\delta`$; the corresponding diagnostic is
$`\varepsilon_{\rm out}+k\delta\ge(8e^2)^{-n}`$.

This is a query restriction on the unchanged local-angle coin and fixed
interleaves. It does not apply to tables reprogrammed with global angle
functions, target-dependent basis changes, or other joint synthesis.
Nor does it prove additive T-counts for oracle calls. In particular,
compact tree scattering is a valid representation, but a proposed pair
of involutions must still supply the boundary transfer, rather than
silently identify the packed step with the prescribed frame.

The [Hopf scattering checks](../../tests/test_hopf_scattering.py) verify the
complete transfer matrices, port permutations, padded dummy convention,
and right-path Fourier coefficient at small dimensions. They test the
algebraic interface, not fine-precision native synthesis of feedback.
