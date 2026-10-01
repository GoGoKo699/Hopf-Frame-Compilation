# Tree transport and the remaining endpoint cost

[Tree generators](SOURCE_REUSE_LIMITS.md#3-tree-generators-compress-the-residual-classically) · [Open endpoint](OPEN_PROBLEM.md) · [Grouped compiler](CONDITIONAL_SUFFIX_COMPILER.md)

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
Apply the [borrowed-workspace compiler](BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements)
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

The [revision decision](OPEN_PROBLEM.md#revision-decision-and-next-bounded-pass)
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
[amplification lemma](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
uses an accepted-block estimate and the actual inverse, rather than
proximity to a prescribed completion. For an arbitrary replacement of this
F-only component, controlled calls and charged workspace still need proof;
it does not automatically implement the affine $`S=A+F`$ branch. The
[residual assembly](RESIDUAL_ASSEMBLY.md) already proves complete-frame
composition for its specified S/R interfaces. A replacement meeting those
interfaces can reuse that theorem with its own native resource ledger.

The [finite checks](../tests/test_tree_transport.py) compare independent
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

The [Cayley fixtures](../tests/test_tree_cayley.py) start with the displayed
four-mode product and extend to eight modes. They compare direct complete
matrices with the recursive data, including actual complex native coarse
words, zero defects, chart boundaries, conditioning, and inverse conversion.
They do not emit a Clifford+T implementation of the Cayley generator or
its inverse.
