# A Gauss-sum table encoding and its missing dirty query

[Session outcome](ONE_CLEAN_FOUR_LAYER_SESSION_RESULT.md) ·
[Retained native query R91](PARALLEL_BILINEAR_WORD_QUERY.md)

The short constructive scout in the 16 October 2026 session tested a
different table mechanism: encode an arbitrary Boolean table in an
affine quadratic Gauss sum. Its native circuit has the desired count
envelope and constant T-depth, but the table sign occurs in a small
matrix element and normalized trace. It does not implement that sign
on arbitrary work inputs. The extraction obstruction below is specific
to this architecture; it is not a general T-depth lower bound.

## 1. An arbitrary table as an inverse entry

Split an a-bit address into $`(l,r)`$, with lengths
$`p=\lfloor a/2\rfloor`$, $`q=\lceil a/2\rceil`$, and write
$`H=2^p,J=2^q,Q=2^a`$. For a Boolean table $`f(l,r)`$, form a
directed acyclic graph from a binary prefix tree with source s and H
leaves u, a disjoint reversed prefix tree from J leaves v to sink t,
and cross edges $`u\to v`$ of constant weight $`f(u,v)`$.
The tree edge at address coordinate i and bit b has weight
$`[x_i=b]`$ over $`\mathbb F_2`$.

There are $`m=(2H-1)+(2J-1)=O(\sqrt Q)`$ vertices, O(m) variable
edges and at most Q constant edges. In topological order let
$`N(x)_{v,u}`$ be the edge weight and put $`A(x)=I+N(x)`$.
This matrix is affine in x and invertible for every address. Its inverse
is the finite path sum $`I+N+\cdots+N^{m-1}`$, so

```math
e_t^{\mathsf T}A(x)^{-1}e_s
=\sum_{u,v}[l=u]f(u,v)[r=v]=f(l,r).
```

No invertibility promise on the table is needed.

## 2. The native encoding and its normalization

On arbitrary m-bit work words y,z define the diagonal involution

```math
D_x|y,z\rangle
=(-1)^{y^{\mathsf T}A(x)z+e_s^{\mathsf T}y+e_t^{\mathsf T}z}|y,z\rangle.
```

Constant quadratic terms are Clifford CZ gates. Each variable tree
edge contributes one CCZ and possibly a constant CZ. All O(m) CCZ
phases use the private-dirty-helper echo of R91 §1 in eight physical
T layers. Shared controls occur only in Clifford gates; each cubic
term has its own arbitrary dirty helper, which returns exactly.
Therefore D and $`U_x=H^{\otimes2m}D_xH^{\otimes2m}`$ have

```math
T=O(\sqrt Q),\qquad G=O(Q),\qquad
b=O(\sqrt Q),\qquad D_T=O(1).
```

The work count includes y,z and the private helpers. The address is
preserved, and the identity holds with arbitrary external references.
However U is a nontrivial operator on y,z.

Summing first over y gives the exact finite-field identity

```math
\sum_{y,z}(-1)^{y^{\mathsf T}Az+b^{\mathsf T}y+c^{\mathsf T}z}
=2^m\sum_{z:Az=b}(-1)^{c^{\mathsf T}z}
=2^m(-1)^{c^{\mathsf T}A^{-1}b}.
```

For $`b=e_s,c=e_t`$, this is $`2^m(-1)^{f(x)}`$. Consequently

```math
\langle0^{2m}|U_x|0^{2m}\rangle
=2^{-2m}\mathrm{Tr}(U_x)=2^{-m}(-1)^{f(x)}.
```

Selecting this matrix element would require 2m initialized inputs
and success probability $`2^{-2m}`$. Those inputs and postselection
are outside the target model. Direct amplitude amplification of this
element would cost $`\Theta(2^m)`$ calls even if initialization were
supplied. An exact integer check of all 16 two-bit tables and all four
addresses confirms the signs and normalization; the displayed sum
proves the arbitrary-size identity.

## 3. A spectral filter cannot extract the sign on dirty work

Since U is an involution, its eigenvalue multiplicities are

```math
n_\pm(x)=\frac{2^{2m}\pm2^m(-1)^{f(x)}}2.
```

Both are positive for every address when $`m\ge1`$. The table sign
changes multiplicities rather than the set of eigenvalues.

Consider any extractor using controlled or uncontrolled calls to U,
interleaved with address-independent unitaries on an external clean
system, and no other action on y,z. On a work eigenvector of eigenvalue
$`\lambda`$, the extractor induces a fixed clean-system unitary
$`V_\lambda`$, independent of the address. For a nonconstant table,
choose two addresses with opposite signs and a +1 eigenvector at each.
The same unit vector $`V_+|0\rangle`$ must approximate both
$`e^{i\phi}|0\rangle`$ and $`-e^{i\phi}|0\rangle`$. Its two squared
distances sum to four. Thus the full initialized operator error is at
least $`\sqrt2`$, even after optimizing a common phase.

This excludes the proposed one-clean spectral filter or ordinary
controlled-U interference word. Additional address-dependent gates
on the Gaussian work lie outside the argument. Such a repair must
provide and price those gates; an uncomputed affine inverse is the
missing computation, not a free operation.

## 4. Scope check for low-degree Boolean phase branches

A separate restricted architecture selects a Boolean phase
$`(-1)^{g_s(x)}`$ with an address-independent coherent seed s and
actually unprepares the seed. Assume every $`g_s`$ has algebraic
degree at most d. Literal full initialized error $`\varepsilon\lt\sqrt2`$
implies

```math
4\Pr_s[g_s(x)\ne f(x)]\le\varepsilon^2
\quad\text{for every }x.
```

Averaging over x shows that one degree-d polynomial is within relative
Hamming distance $`\delta=\varepsilon^2/4\lt1/2`$ of f.
There are $`2^M`$ such polynomials, where
$`M=\sum_{j=0}^d\binom aj`$, and each radius-$`\lfloor\delta Q\rfloor`$
Hamming ball contains at most $`2^{H_2(\delta)Q}`$ truth tables.
Covering all Boolean tables therefore requires

```math
\sum_{j=0}^d\binom aj\ge(1-H_2(\delta))Q.
```

For fixed error this requires $`d\ge a/2-O(\sqrt a)`$. An arbitrary
function of r linear hash bits has degree at most r; XORs preserve that
bound and bilinear fusion adds the two banks' degree bounds. Thus
uniformly sublinear total hash rank cannot meet this particular
Boolean-branch contract. General coherent branches, address-dependent
seed preparation and arbitrary Clifford+T circuits are not covered.

## 5. Reopening condition

The [subsequent inverse-extraction analysis](GAUSSIAN_INVERSE_EXTRACTION.md)
gives an exact dirty-work Pauli commutator, including controlled query
extraction and a native baseline. It pays for the in-place affine map
whose inverse entry contains f. A tight affine-symplectic round bound
limits that particular repair. The later
[scalar-loop theorem](GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md) also limits
direct scalar extraction by r actual affine Gaussian phase rounds:
the returned sign has degree at most $`4r+2`$. Its path-sum and Pauli
marker invariant sees the phase that the work symplectic matrix misses.
This is still a round-architecture restriction; general native scalar
extraction remains open.

The Gauss-sum route needs a native conversion from the affine-quadratic
encoding to the scalar table phase, returning every arbitrary work
wire, with the same count envelope and sublogarithmic address T-depth.
It must handle occupied pure inputs rather than select a small matrix
element or assume maximally mixed work. No such conversion qualified
in the scout. R91, the full frame upper bound R96 and the indicator
and routing costs remain unchanged.

## 6. What the complete inverse-word query does and does not transfer

The [subsequent garbage-recursive construction](NONCOMMUTATIVE_CYCLE_CHARACTER.md#5d-a-logarithmic-depth-complete-inverse-word-query)
supplies an exact complete inverse-word XOR query on every physical
unitriangular coefficient word. It returns arbitrary work, with cubic
T-count, Clifford count and dirty width, and logarithmic T-depth in m.
Its coefficient wires change internally. This removes the matrix-element
extraction issue for a **supplied physical matrix**, but does not supply
the affine coefficient embedding or native read-only access to it.

For Section 1's $`m=\Theta(\sqrt Q)`$ table representation, the generic
substitution has $`\Theta(Q^{3/2})`$ T-count and dirty work and
$`\Theta(Q)`$ physical coefficient positions. Its stated depth guarantee
is $`O(a)`$, compared with R91's $`O(\log(a+2))`$. These are costs of
the specified generic circuit, not lower bounds on a specialized sparse
implementation. The O(m) variable edges and Q classically specified
constant edges need not be treated identically by a different compiler.

### A scoped affine-embedding capacity bound

Put $`s=\binom m2`$. Allow every strict-upper entry of U(x) to be an
arbitrary table-dependent affine Boolean function of the a address bits.
Allow the final scalar readout to be any affine linear combination of all
s strict inverse entries. There are $`2^{(a+1)s}`$ matrix specifications
and $`2^{s+1}`$ readouts. Hence this architecture represents at most

```math
2^{(a+2)s+1}
```

Boolean tables. Covering all $`2^Q`$ tables requires

```math
s\ge\frac{Q-1}{a+2}.                                  \qquad(7)
```

An additional affine function of x in the readout replaces the final
one in the exponent by a+1, without changing the asymptotic bound.
Selecting only one inverse entry is a smaller architecture. For entries
of Boolean degree at most k, the same argument replaces a+1 by
$`M_k=\sum_{j=0}^k\binom aj`$; computing those coefficients is a
separate paid task.

Thus shrinking matrix order cannot make the retained cubic generic
inversion ledger fit the square-root T/work budget within this affine
embedding. The full coefficient word alone has $`\Omega(Q/a)`$
positions. This statement does not cover arbitrary quantum Clifford
encodings, nonlinear readouts, or a sparse specialization that removes
fixed positions without initialized work.

### Completed oracle calls with affine wiring

A related boundary allows unlimited arbitrary borrowed width. Consider
k completed inverse-word queries, interleaved only with computational-basis
affine circuits, such as X, CNOT and SWAP. Call j has matrix order m_j
and exposes $`s_j=\binom{m_j}2`$ strict inverse entries. Its own
work returns. Allow all external wires to participate in the affine
circuits; this is more permissive than native read-only addressing.

Evaluate the required full-column XOR-query identity at zero output and
zero borrowed input. This is a test column, not a runtime preparation.
Initially every wire lies in the affine span of $`1,x_1,\ldots,x_a`$.
Affine wiring preserves this span. Each completed query adds at most
s_j new Boolean features; its diagonal entries add only the constant
one. Thus the physical borrowed width does not enlarge the feature
span beyond the previous query outputs.

For fixed call orders, call j's input is specified by at most
$`s_j(a+1+\sum_{i\lt j}s_i)`$ binary coefficients. With
$`M=\sum_j s_j`$, the final output uses at most a+1+M more. Therefore
all represented tables have at most

```math
M(a+2)+\sum_{i\lt j}s_i s_j+a+1
\le M(a+2)+M^2/2+a+1                                  \qquad(8)
```

parameter bits. Universal coverage requires this quantity to be at least
Q and hence $`M=\Omega(\sqrt Q)`$. At one call (7)'s stronger
$`\Omega(Q/a)`$ bound returns.

For fixed k, the **retained generic substitution ledger** charges
$`\Theta(\sum_j s_j^{3/2})`$ T gates for orders at least three.
Smaller uncontrolled queries are affine and can be absorbed into the
wiring. Convexity gives

```math
\sum_j s_j^{3/2}\ge M^{3/2}/\sqrt k,
```

so universal coverage by a fixed number of these generic calls costs
$`\Omega(Q^{3/4}/\sqrt k)`$, outside the square-root T envelope.
Allowing call orders to vary below a fixed budget adds only
$`O(k\log(M+1))`$ specification bits and leaves this conclusion
unchanged. A controlled completed query adds one affine control form
per call and replaces its new features by control-times-inverse entries.
At fixed k this adds only $`O(k(a+M))`$ parameter bits. Its charged
generic cost has the same scaling, including the nonlinear order-two
controlled call.

This is a capacity and **generic-call pricing** statement. It is not
a lower bound on arbitrary native synthesis of these functions. Arbitrary
Clifford interference, exposing internal compute garbage, nonlinear
readout, and cancellation-aware resynthesis lie outside its hypotheses.

### Why dirty loading needs an actual extraction mechanism

For fixed x, a read-only-address unitary is onto the entire arbitrary-work
space. It cannot force a separate nonempty coefficient register to U(x),
independently of every dirty input: fixing s output bits reduces that
space's dimension by $`2^s`$, even if other work may retain garbage.
The available XOR loading instead gives E+U(x).

The resulting two-query translation echo is a finite difference, not
evaluation at a clean matrix. Already at order three the corner inverse
entry is $`w+uv`$. Translating u,v by x,y gives

```math
[w+(u+x)(v+y)]+[w+uv]=xv+yu+xy.                       \qquad(9)
```

The dirty terms survive. At x=y=1, dirty inputs (u,v)=(0,0) and (1,0)
give opposite phases; the erroneous phase echo has operator distance two
from the desired $`(-1)^{xy}I`$ on one such sector. A four-corner echo
repairs this quadratic example, so (9) excludes only the naive two-call
loading claim.

## 7. An exact read-only repair with a losing sequential ledger

There is an exact way to remove a common dirty matrix offset. It does
not improve the retained read-only query baseline. For matrix order
$`n\ge2`$, let $`p=2^b\ge n`$ be the smallest such power and
$`K=\mathbb F_p`$, in a fixed binary basis containing one. For a binary
strict-upper U and an arbitrary K-valued dirty strict-upper E, define

```math
(D_tU)_{ij}=t^{j-i}U_{ij},\qquad
F_{ij}(t)=[(I+E+D_tU)^{-1}]_{ij},\qquad k=j-i.
```

Every increasing path from i to j has total index gap k. Its term has
t-degree at most k, with equality exactly when every edge chooses its U
contribution. Therefore $`[t^k]F_{ij}(t)=[(I+U)^{-1}]_{ij}`$,
independently of E. The finite-field power sum gives the exact identity

```math
\sum_{t\in K}t^{p-1-k}F_{ij}(t)=[(I+U)^{-1}]_{ij}.     \qquad(10)
```

Use $`0^0=1`$. The endpoint k=p-1 is included: the constant power sum
is p=0 in characteristic two, whereas $`\sum_t t^{p-1}=1`$.
Apply a fixed binary coordinate functional $`\lambda:K\to\mathbb F_2`$
with $`\lambda(1)=1`$ to read the binary answer.

One evaluation loads $`E\mathrel{+}=D_tU`$ by CNOTs from U, runs
the field version of the complete inverse-word query with weighted
readout $`\lambda(t^{p-1-k}F_{ij}(t))`$ into each arbitrary output,
and undoes the actual load. E and the internal work return after every
evaluation. Repeating over all t with the **same** initial E and toggling
the known diagonal ones once gives the full inverse-word XOR query.
U occurs only as a native Clifford control. Every identity holds on
arbitrary work, outputs and external references.
At n=1 the query is just the known diagonal output X, with no field,
T gate or helper; its actual inverse is charged in the fixed alphabet.

The field version has a complete native substitution. One field product
is b-squared binary products, each contributing a fixed linear form to
the output bits. The private-helper phase echo reads that linear form
without a clean parity register. This uses $`O(b^2)`$ T gates,
$`O(b^3)`$ elementary Cliffords, constant T-depth and $`O(b^2)`$
private dirty helpers per field product. The same two-shear recursion
therefore costs $`O(n^3b^2)`$ T/work, $`O(n^3b^3)`$ Cliffords and
$`O(\log n)`$ T-depth per complete evaluation. The live E word has
$`b\binom n2`$ separate bits and is included in the work bound.
The p sequential evaluations give, with zero clean work,

```math
T=O(n^4\log^2 n),\quad G=O(n^4\log^3 n),\quad
D_T=O(n\log n),\quad W=O(n^3\log^2 n).               \qquad(11)
```

The constants, field structure and unit diagonals are compiled classical
data, not initialized coefficient registers. Every actual inverse is
included. The direct phase readout and a private-helper controlled readout
obey the same asymptotic bounds. A full-space compute error delta per
evaluation, followed by its actual inverse, gives total error at most
$`2p\delta`$; the displayed native construction is exact.

For comparison, direct read-only path propagation already gives all
inverse entries with $`O(n^3)`$ T/Clifford count, $`O(n)`$ depth and
$`O(n^2)`$ arbitrary work. On an arbitrary n-by-n B, process rows from
n-2 down to zero, adding $`U_{ij}B_j`$ into row i. This realizes
$`F:B\mapsto(I+U)^{-1}B`$ using U only as native Clifford controls.
The chronological word $`F^\dagger,X_{\rm diagonal},F`$ translates
B by the full inverse matrix. Echo this translation around CNOT readout
to the desired output and use its actual inverse. B and every private
helper return. At one row the products have stable sources and can share
one enclosing Hadamard per target, so the stated count/depth/work follow
from the retained private-helper schedule. Interpolation is slower and
wider than this already available path/Weyl construction.

### The shared-mask parallelization gap

A fixed K-linear evaluation formula extracting the degree-k coefficient
from all the polynomials above needs at least k+1 distinct points. Their
span contains each lower monomial: for $`0\lt\ell\lt k`$, use one
U edge of gap ell followed by one E edge of gap k-ell; a direct E edge
gives the constant and a direct U edge gives $`t^k`$. The E coefficient
allows arbitrary K multiples of every lower monomial; the monic leading
term already suffices. If r points with
$`r\le k`$ sufficed, the monic polynomial

```math
t^{k-r}\prod_{j=1}^r(t-\alpha_j)
```

would vanish at every sampled point but have leading coefficient one.
The same contradiction holds for fixed binary-linear field readouts
selecting a coordinate with $`\lambda(1)=1`$.

Independent dirty matrices do not substitute for the common E in (10).
More generally, a product of diagonal scalar phase factors, each with
its own independent arbitrary mask, can be mask-independent only if
each factor separately is. An invertible address-independent linear
change of mask coordinates preserves their independence. This excludes
only uncoupled copies and fixed linear evaluation formulas; it leaves
coupled-mask processing, nonlinear/coherent extraction and other native
representations open.

[Bounded exact fixtures](../../tests/test_gaussian_inverse_extraction.py)
check the dirty-offset failure, finite-field extraction including the
maximal exponent, and the evaluation boundary. These support the
interfaces, not a general native lower bound. A qualifying R98 mechanism
still needs an actual table embedding and returned extraction at its
simultaneous square-root count/work envelope, native read-only address
contract and sublogarithmic address T-depth.
