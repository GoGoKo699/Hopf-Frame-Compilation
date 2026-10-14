# A sharp normalization limit for entry-indexed frame encodings

[Complete frame convention](../../docs/HOPF_INTERFACE.md) ·
[Sublinear-depth construction](SUBLINEAR_T_DEPTH.md#3-compress-the-complete-columns)

The balanced h-bit real Hopf frame has entrywise absolute-value norm
$`\sqrt{h+1}`$. This is exactly the best normalization attainable by
the entry-indexed row/column preparation architecture defined below.
Even a fixed operator-norm approximation error cannot change that order.
Thus changing weights within this architecture cannot remove the
worst-case square-root normalization in the current tag construction.

This is an algebraic limit on a specified encoding architecture. It is
not a lower bound on unrestricted block encodings, recovery algorithms,
or native T-depth. No circuit-cost assertion is made for the optimal
preparations used to prove sharpness.

## 1. The encoding contract

Let B be a matrix whose entries have a finite factorization

```math
B_{ij}=\sum_k\overline{r_{ijk}}c_{ijk},\qquad
\sum_{j,k}|r_{ijk}|^2\le1\quad\text{for every }i,\qquad
\sum_{i,k}|c_{ijk}|^2\le1\quad\text{for every }j.
```

This is the overlap of row and column preparations in an orthonormal
basis indexed by the shared entry $`(i,j)`$ and an auxiliary label k:
the row state for i has support on $`|i,j,k\rangle`$, and the column
state for j has support on $`|i,j,k\rangle`$. Additional orthogonal
failure components may complete the states to unit norm. Arbitrary
complex amplitudes and interference among k values are allowed.

An encoding with normalization $`\alpha\gt0`$ represents
$`\widetilde W=\alpha B`$. Throughout this note, $`|A|`$ means
entrywise absolute value, $`\|A\|`$ is the spectral/operator norm,
and $`\|A\|_*`$ is the nuclear norm.

**Normalization lemma.** Every such factorization satisfies

```math
\||B|\|\le1,\qquad \alpha\ge\||\widetilde W|\|.
```

For nonnegative unit vectors x and y, triangle inequality and
Cauchy–Schwarz give

```math
\begin{aligned}
x^{\mathsf T}|B|y
&\le\sum_{i,j,k}x_i|r_{ijk}|\,y_j|c_{ijk}|\\
&\le\left(\sum_{i,j,k}x_i^2|r_{ijk}|^2\right)^{1/2}
     \left(\sum_{i,j,k}y_j^2|c_{ijk}|^2\right)^{1/2}
\le1.
\end{aligned}
```

The spectral norm of a nonnegative matrix is obtained by maximizing
over nonnegative unit x and y: replacing arbitrary vectors by their
entrywise magnitudes cannot decrease the absolute bilinear form. This
proves the lemma. In particular, adding entry-local channels or phases
does not evade it.

## 2. The balanced frame and its complete absolute-value spectrum

Fix $`h\ge1`$ and $`Q=2^h`$. Set every angle of the h-bit Hopf
frame to $`\pi/4`$, using $`R_y(\theta)=e^{-i\theta Y}`$.
Its columns are the normalized constant vector

```math
u=Q^{-1/2}\mathbf 1
```

and one signed Haar vector $`w_{d,p}`$ for every depth
$`0\le d\le h-1`$ and d-bit prefix p. This marker column is
zero outside the subtree p, equals $`-2^{-(h-d)/2}`$ on its child
p0, and equals $`+2^{-(h-d)/2}`$ on its child p1. These are exactly
the columns indexed by the marker $`p1\,0^{h-d-1}`$ in the complete
frame convention. Denote this real orthogonal matrix by W and put
$`A=|W|`$.

Let $`g_{d,p}=|w_{d,p}|`$ be the normalized positive indicator of
subtree p, and let

```math
P_d=\sum_{p\in\{0,1\}^d}g_{d,p}g_{d,p}^{\mathsf T}.
```

Thus $`P_d`$ is the rank-$`2^d`$ orthogonal projector onto vectors
constant on each depth-d subtree, and $`P_0=uu^{\mathsf T}`$.
Summing the outer products of the absolute-value columns gives

```math
AA^{\mathsf T}=P_0+\sum_{d=0}^{h-1}P_d
              =2P_0+P_1+\cdots+P_{h-1}.
```

The projectors are nested. The constant vector has eigenvalue
$`h+1`$. A depth-d Haar vector lies in the kernel of $`P_e`$
for $`e\le d`$ and in its range for $`e\ge d+1`$. The complete
eigenvalue list is therefore

```math
h+1\quad\text{with multiplicity }1,\qquad
h-1-d\quad\text{with multiplicity }2^d,
\quad d=0,\ldots,h-1.
```

The last depth supplies the $`2^{h-1}`$ zero eigenvalues. In
particular,

```math
\||W|\|=\sqrt{h+1}.
```

Equivalently, the root-column basis vector and the h normalized uniform
vectors on the disjoint depth sectors of the column index are
orthonormal, and A sends every one of them to u. The projector calculation
also proves that no larger singular value is present.

For $`h=1`$ the projector sum is $`2P_0`$, with spectrum
$`2,0`$. Explicitly,
$`W=2^{-1/2}\begin{pmatrix}1&-1\\1&1\end{pmatrix}`$ and
$`A=2^{-1/2}\mathbf1\mathbf1^{\mathsf T}`$; no intermediate
depth is needed.

## 3. Exact sharpness within the architecture

For each leaf i there are exactly $`h+1`$ nonzero entries in its
row of W: the root column and the marker for its ancestor at every
depth. Use a single k value and set

```math
r_{ij}=\begin{cases}
1/\sqrt{h+1},&W_{ij}\ne0,\\
0,&W_{ij}=0,
\end{cases}
\qquad c_{ij}=W_{ij}.
```

Every row of r has squared norm one. Every column of c has squared
norm one because W is orthogonal. The resulting overlap is

```math
B_{ij}=\overline{r_{ij}}c_{ij}
      =\frac{W_{ij}}{\sqrt{h+1}}.
```

Together with the normalization lemma, this proves that the minimum
exact normalization in this architecture is $`\sqrt{h+1}`$.
This is a statement about allowed amplitudes and their normalization;
it does not charge or promise a native implementation of these states.

## 4. Robustness to operator-norm approximation

Suppose $`0\le\epsilon\lt1`$ and
$`\|\widetilde W-W\|\le\epsilon`$, allowing $`\widetilde W`$
and its error to be complex. Then

```math
\||\widetilde W|\|\ge(1-\epsilon)\sqrt{h+1}.
```

To prove this, keep u from Section 2 and define the nonnegative unit
vector

```math
v=\frac{A^{\mathsf T}u}{\sqrt{h+1}}.
```

Its root entry is $`1/\sqrt{h+1}`$, and its entry for any depth-d
marker is $`2^{-d/2}/\sqrt{h+1}`$. Its squared norm is
$`(1+\sum_{d=0}^{h-1}2^d2^{-d})/(h+1)=1`$, and
$`u^{\mathsf T}Av=\sqrt{h+1}`$.

Define a real dual witness by

```math
F_{ij}=u_i v_j\,\mathrm{sign}(W_{ij}),\qquad
\mathrm{sign}(0)=0.
```

It has the exact factorization $`F=WD`$, where D is positive diagonal
in the root/marker column basis, with entries

```math
D_{\rm root}=\frac1{\sqrt{h+1}},\qquad
D_{d,p}=\frac{2^{-d}}{\sqrt{h+1}}.
```

Indeed, on the support of a depth-d marker,
$`|W_{ij}|D_{d,p}=2^{-(h+d)/2}/\sqrt{h+1}=u_iv_j`$.
Both sides vanish off that support. The root column gives the same
identity with its displayed diagonal entry. Since W is orthogonal,

```math
\|F\|_*=\mathrm{Tr}\,D
=\frac{1+\sum_{d=0}^{h-1}2^d2^{-d}}{\sqrt{h+1}}
=\sqrt{h+1}.
```

Write $`E=\widetilde W-W`$. For each entry, including zeros of W,

```math
|W_{ij}+E_{ij}|
\ge |W_{ij}|+
\Re\!\left(\mathrm{sign}(W_{ij})E_{ij}\right).
```

At a zero entry the right side is zero. Multiplying by the nonnegative
weights $`u_iv_j`$ and summing proves

```math
\begin{aligned}
\||\widetilde W|\|
&\ge u^{\mathsf T}|\widetilde W|v\\
&\ge\sqrt{h+1}+\Re\,\mathrm{Tr}(F^\dagger E)\\
&\ge\sqrt{h+1}-\|F\|_*\|E\|\\
&\ge(1-\epsilon)\sqrt{h+1}.
\end{aligned}
```

Consequently an entry-indexed encoding obeying the explicit error
contract $`\|W-\alpha B\|\le\epsilon`$ must have

```math
\alpha\ge(1-\epsilon)\sqrt{h+1}.
```

This approximate bound is sharp as well: use the exact B from Section 3
and $`\alpha=(1-\epsilon)\sqrt{h+1}`$. Then
$`\alpha B=(1-\epsilon)W`$ and the operator error is exactly
$`\epsilon`$. The error here is on the represented operator
$`\alpha B`$, not on the scaled block B alone.

At $`h=1`$, $`D=I/\sqrt2`$ and $`F=W/\sqrt2`$, so the same
witness has nuclear norm $`\sqrt2`$ and proves the stated bound.

## 5. Relation to the current tag construction and the remaining freedom

The early block in the sublinear-depth theorem has
$`h=n-M`$ logical bits. Its tag source uses
$`K=2^{\lceil\log_2(h+1)\rceil}`$ values and normalization
$`\sqrt K`$. On the balanced instance it belongs to the architecture
above: a used shared label $`(i,t)`$ identifies a unique entry
$`(i,j)`$, because t specifies the root or depth and i specifies the
ancestor prefix. Padded tags contribute only unmatched row components.
Relabeling this common basis is an algebraic identification and requires
no additional circuit register.

Hence its exact normalization is within the factor

```math
1\le\frac{\sqrt K}{\sqrt{h+1}}\lt\sqrt2
```

of the architecture optimum. Arbitrary reweighting and extra entry-local
channels cannot improve its worst-case normalization to
$`o(\sqrt h)`$, even at any fixed represented-operator error below one.

The theorem does not constrain preparations whose shared basis lacks
the entry-indexed support contract, coherent changes of matrix basis,
other unitary factorizations, or direct native circuits. It also does
not supply a query lower bound for recovery from a block encoding:
such a statement requires its own specified oracle and access model.
The present conclusion identifies the limit of optimizing this
normalization, while leaving those other routes open.
