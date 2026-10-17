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
limits that particular repair; direct scalar extraction remains open.

The Gauss-sum route needs a native conversion from the affine-quadratic
encoding to the scalar table phase, returning every arbitrary work
wire, with the same count envelope and sublogarithmic address T-depth.
It must handle occupied pure inputs rather than select a small matrix
element or assume maximally mixed work. No such conversion qualified
in the scout. R91, the full frame upper bound R96 and the indicator
and routing costs remain unchanged.
