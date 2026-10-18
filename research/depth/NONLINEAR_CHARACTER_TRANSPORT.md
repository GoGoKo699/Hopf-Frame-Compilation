# Nonlinear character transport: a returned native mechanism and its width cost

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Affine diagonal-source obstruction](RETURNED_DIAGONAL_ECHO_DEGREE.md) ·
[Native query baseline R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Conditional frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md)

Nonlinear dirty-work transport can expose a Boolean phase of unbounded
degree using only two diagonal source occurrences. The construction below
is a literal full-space identity and has a completely charged, slow native
implementation, including an arbitrary-output XOR version. It escapes the
fixed-affine transport hypothesis of the diagonal-source theorem.

Its apparent advantage does not survive the workspace comparison. Every
character in this unit-group representation, and more generally in every
commutative nilpotent coefficient algebra below, has Boolean degree
at most logarithmic in its physical dirty width. A degree-a target requires
at least $`2^{a-1}`$ work bits even if source synthesis and transport are
optimistically free. This excludes the representation from R98's
square-root-width regime. It is not a lower bound on other nonlinear
representations or on general native query depth. The
[noncommutative cycle character](NONCOMMUTATIVE_CYCLE_CHARACTER.md)
does escape this width obstruction, with a separately charged native
construction that does not improve the retained conjunction or R91 bounds.

## 1. A nonlinear full-return identity

For $`w\ge1`$, let

```math
R_w=\mathbb F_2[t]/(t^{w+1}),\qquad G_w=1+tR_w.
```

Every arbitrary w-bit work word z represents exactly one unit

```math
u_z=1+\sum_{i=1}^w z_i t^i.
```

Thus the group has exactly all physical work basis states. There is no
initialized word, valid-code promise, padding sector or supplied eigenstate.
Let the unchanged address be $`x\in\mathbb F_2^a`$, and choose

```math
a_x=1+\sum_{i=1}^w\ell_i(x)t^i,
```

where every $`\ell_i`$ is an affine Boolean form. Let
$`\chi:G_w\to\mathbb F_2`$ be any group homomorphism. Define

```math
D_\chi|u\rangle=(-1)^{\chi(u)}|u\rangle,\qquad
P_x|u\rangle=|a_xu\rangle.
```

Then the literal operator identity is

```math
D_\chi P_xD_\chi^\dagger P_x^\dagger
       =(-1)^{\chi(a_x)}I_W.                              (1)
```

Indeed the rightmost transport sends u to $`a_x^{-1}u`$; the two source
phases add to $`\chi(a_x^{-1}u)+\chi(u)=\chi(a_x)`$. The complete
work word returns on every input. This proves the identity on arbitrary
superpositions and with arbitrary external references, including all
literal signs and phases.

Multiplication is affine in work at a fixed address, but its complete
coordinate map is generally quadratic in the joint address and work:

```math
z_i\longmapsto z_i+\ell_i(x)
        +\sum_{j=1}^{i-1}\ell_j(x)z_{i-j}.                (2)
```

It is therefore not an allowed affine Clifford interval in the earlier
diagonal-source theorem. The inverse is multiplication by $`a_x^{-1}`$,
not a second copy of (2). Both appearances require actual native words.

## 2. All characters have logarithmic degree in the work width

Use the Euler derivation $`\delta=t\,d/dt`$, which preserves the ideal
$`(t^{w+1})`$ and hence is well defined on $`R_w`$. For odd
$`n\le w`$, define

```math
p_n(u)=[t^n]\frac{\delta u}{u}.                           (3)
```

The product rule gives $`p_n(uv)=p_n(u)+p_n(v)`$, so these are Boolean
characters. The square subgroup consists exactly of units supported on
even powers of t, and has size $`2^{\lfloor w/2\rfloor}`$.
Consequently $`G_w/G_w^2`$ has dimension $`\lceil w/2\rceil`$.
For odd i and n,

```math
p_n(1+t^i)=\begin{cases}1,&i\mid n,\\0,&i\nmid n.\end{cases}
```

Ordering the odd indices increasingly gives a triangular matrix with
ones on its diagonal. The odd $`p_n`$ therefore form a basis of **all**
$`\mathbb F_2`$-valued characters, not merely a selected source family.

Write $`u=1+b`$. In the truncated ring,

```math
u^{-1}=\prod_{j\ge0}(1+b^{2^j}),\qquad
b^{2^j}=\sum_{i=1}^w z_i t^{i2^j}.                       (4)
```

Only finitely many factors are nonidentity. The second equality uses
Frobenius and the Boolean values of the coefficient bits. Each nonconstant
factor has coefficient degree one. Choosing k distinct nonconstant factors
costs t-degree at least $`1+2+\cdots+2^{k-1}=2^k-1`$.
The numerator $`\delta u`$ adds t-degree at least one and Boolean degree
at most one. Thus

```math
\deg p_n\le1+\lfloor\log_2 n\rfloor,
\qquad
\deg\chi\le1+\lfloor\log_2 w\rfloor.                  (5)
```

Cancellation or reduction of repeated Boolean variables cannot increase
these degrees. Substituting the affine coefficient forms $`\ell_i(x)`$
also cannot increase them. Every returned phase in (1) therefore obeys

```math
\deg_x\chi(a_x)\le1+\lfloor\log_2 w\rfloor.             (6)
```

In particular, an a-bit target of degree a requires

```math
w\ge2^{a-1}=Q/2,\qquad Q=2^a.                            (7)
```

These are physical, independently arbitrary coefficient wires. For R98
with one output and extra dirty pool $`B=\Theta(\sqrt Q)`$, (7) exceeds
the available allocation for growing a. Even using the one arbitrary output
as a coefficient wire would only change the width comparison by one.
This conclusion already grants every source and transport call free cost;
optimizing their native implementations cannot remove this width requirement
within the stated representation.

### The same width barrier for every commutative nilpotent algebra

Let A be a nonzero finite-dimensional commutative associative nilpotent
$`\mathbb F_2`$ algebra, with $`w=\dim A`$ and nilpotency index
$`\nu=\min\{j:A^j=0\}`$. The unitization gives a group $`1+A`$
whose $`2^w`$ elements are exactly all coefficient words. For any
character $`\chi:1+A\to\mathbb F_2`$ and affine coefficient map

```math
u_x=1+b+\sum_{i=1}^a x_i a_i,\qquad b,a_i\in A,
```

the degree obeys

```math
\deg_x\chi(u_x)\le\lceil\log_2\nu\rceil
                 \le1+\lfloor\log_2 w\rfloor.            \qquad(7a)
```

This includes multivariate quotient algebras; it needs no classification
of their characters or square subgroups. The zero algebra has only the
trivial character and can be omitted.

For a proof, use commuting indeterminates and the universal polynomial

```math
F_d(s)=\prod_{\epsilon\in\mathbb F_2^d}
           \left(s+\sum_{j=1}^d\epsilon_jt_j\right)
       =s^{2^d}+\sum_{i=0}^{d-1}c_i(t)s^{2^i}.
```

It is additive in s, and each nonzero $`c_i`$ is homogeneous of degree
$`2^d-2^i`$. Induction follows from
$`F_{d+1}(s)=F_d(s)^2+F_d(t_{d+1})F_d(s)`$. No linear independence
of the substituted directions is assumed. Substitute $`s=1+b`$ and
arbitrary $`t_j\in A`$ in the commutative unitization. Frobenius and
the coefficient degrees give

```math
\prod_{\epsilon\in\mathbb F_2^d}
       \left(1+b+\sum_j\epsilon_jt_j\right)
    \in 1+A^{2^{d-1}}.                                  \qquad(7b)
```

When $`2^{d-1}\ge\nu`$, this product is one. Applying the character
identifies its value with the d-fold Boolean derivative of $`\chi(u_x)`$
in any address directions, so every such derivative vanishes. Set
$`d=1+\lceil\log_2\nu\rceil`$ to prove the first inequality in (7a).
The powers $`A,A^2,\ldots`$ strictly descend until zero: an equality
$`A^j=A^{j+1}`$ would persist under multiplication by A and contradict
nilpotence unless both were zero. Therefore $`\nu\le w+1`$, proving
the second inequality.

In particular, a degree-a target still requires $`w\ge2^{a-1}`$ in
every such commutative coefficient representation, before charging any
source or transport. Commutativity is essential: the cycle construction
linked above realizes degree d on $`d^2`$ coefficients using a
noncommutative path algebra. The cube argument does not apply there.

## 3. The nonlinear escape really has unbounded degree

For every $`d\ge1`$, put

```math
w_d=1+(d-1)2^{d-1},\qquad
a_x=1+\sum_{j=1}^d x_jt^{2^{j-1}},\qquad
\chi=p_{w_d}.                                           (8)
```

The returned phase $`\chi(a_x)`$ has degree exactly d. To prove this,
take its full ANF coefficient, summing over all d-bit addresses. Because
$`\delta a_x=x_1t`$, this coefficient is the coefficient of
$`t^{w_d}`$ in

```math
t\sum_{v\in V}\frac1{1+t+v},\qquad
V=\mathrm{span}_{\mathbb F_2}
          \{t^2,t^4,\ldots,t^{2^{d-1}}\}.                 (9)
```

Let $`F_V(s)=\prod_{v\in V}(s+v)`$. This is an additive polynomial:
induct on a basis of V using

```math
F_{V+\mathbb F_2b}(s)=F_V(s)^2+F_V(b)F_V(s).
```

Its derivative is therefore the constant
$`\prod_{v\in V\setminus\{0\}}v`$. The series (9) equals

```math
\frac{t\prod_{v\in V\setminus\{0\}}v}{F_V(1+t)}.
```

The denominator has constant coefficient one. For each
$`2\le j\le d`$, exactly $`2^{d-j}`$ nonzero v have least
monomial $`t^{2^{j-1}}`$. The numerator consequently has t-adic order

```math
1+\sum_{j=2}^d2^{d-j}2^{j-1}=w_d
```

and leading coefficient one. The full ANF coefficient is one, proving the
claim. This construction has only two diagonal source occurrences but
$`\Theta(d2^d)`$ physical work bits. It supplies an actual nonlinear
escape from the fixed-affine source model, not a smaller-width query.

## 4. A completely native source, transport and XOR query

Every unit has the unique factorization

```math
u=\prod_{j=1}^w(1+t^j)^{c_j},\qquad c_j\in\mathbb F_2.   (10)
```

Successively cancel the first remaining coefficient to obtain existence
and uniqueness. Multiplication adds factor exponents; each binary carry
uses $`(1+t^j)^2=1+t^{2j}`$ and reaches only an even index.
The odd factor bits $`c_j`$ are therefore characters and form another
basis of all characters.

### In-place character coordinates

A fixed reversible word F computes the factor bits on the same w dirty
wires. For $`j=1,\ldots,w`$, execute, in ascending i,

```math
z_i\longmapsto z_i+z_jz_{i-j},\qquad 2j+1\le i\le w.    (11)
```

There are exactly

```math
N_w=\sum_{j=1}^w\max\{w-2j,0\}
             =\left\lfloor\frac{(w-1)^2}{4}\right\rfloor
```

Toffolis. At the start of stage j, lower wires retain the preceding factor
bits; unprocessed wires represent a residual series whose actual lower
coefficients are zero. Its degree-j coefficient is $`c_j=z_j`$.
Dividing that residual by $`1+c_jt^j`$ leaves coefficients
$`j\lt i\le2j`$ unchanged, since contributions from the constant and
$`c_jt^j`$ cancel. For $`i\gt2j`$ the new coefficient satisfies
$`r_i=z_i+c_jr_{i-j}`$, precisely the ascending update (11).
The stored earlier factor bits are never used as residual coefficients.
This proves F on every work input. Its actual reversed native word is
its full inverse.

For a scalar character, compute F, apply Z on its selected odd factor
bits, and execute the actual inverse of F. Thus $`D_\chi`$ has an
exact native implementation with $`O(w^2)`$ T-count, elementary Clifford
count and T-depth, and no additional helper. This does not assume a
supplied nonlinear change of basis.

For q arbitrary outputs $`h_1,\ldots,h_q`$ and characters
$`\chi_1,\ldots,\chi_q`$, replace those Z gates by a fixed CZ matrix
from the outputs to the corresponding odd factor bits. The resulting
source is literally

```math
D_{\mathbf h}|u\rangle
  =(-1)^{\sum_{k=1}^q h_k\chi_k(u)}|u\rangle.
```

If the CZ matrix has R nonzero entries, $`R\le qw`$, its elementary
Clifford cost is exactly $`3R`$, using H, CNOT, H with the output as
CNOT control. The retained 7-T, 10-CNOT CCZ word, conjugated by two target
Hadamards, implements each Toffoli in (11) with 12 elementary Cliffords
and at most five T layers. Therefore one controlled source has

```math
T_D=14N_w,\qquad G_D=24N_w+3R,\qquad d_D\le10N_w.        (12)
```

Its actual inverse has the same counts for this specified H/CNOT/T word.
No output or address wire is a T target in this source.

### Paid address-dependent transport

To implement (2), visit $`i=w,w-1,\ldots,1`$. Toggle $`z_i`$ by
$`\ell_i(x)`$, then by every product
$`\ell_j(x)z_{i-j}`$ with $`j\lt i`$. The lower coefficient wires
have not yet changed, so this is exactly multiplication by $`a_x`$.
Affine loads use explicit Clifford gates. For each product, the R91
private-dirty phase echo, conjugated by target Hadamards, supplies eight
T gates in eight T layers for its linear-address part. Its phase target
is a separate arbitrary helper. A constant part adds a CNOT. Every
address bit occurs only as a Clifford control, and the private helper
returns after each product.

There are at most $`w(w-1)/2`$ such products, so this sequential transport
and its actual reverse have

```math
T_P=O(w^2),\qquad G_P,G_{P^\dagger}=O((a+1)w^2),\qquad
 d_P=O(w^2),
```

using one reused arbitrary helper. The inverse Clifford count is charged
separately: a native affine load may contain an X decomposition whose
literal inverse has a different count in the H/S/CNOT alphabet.

Substitute $`D_{\mathbf h}`$ for $`D_\chi`$ in (1). Its phase is
$`(-1)^{\sum_k h_k\chi_k(a_x)}`$ and all coefficient work returns.
Conjugating all q arbitrary outputs by Hadamards gives the complete XOR
query $`\mathbf h\mapsto\mathbf h+\boldsymbol\chi(a_x)`$.
Those $`2q`$ Hadamards are charged. The literal four-occurrence word has

```math
\begin{aligned}
T&=2T_D+2T_P=O(w^2),\\
G&=2G_D+G_P+G_{P^\dagger}+2q
                         =O((a+1)w^2+qw),\\
D_T&\le2d_D+2d_P=O(w^2).
\end{aligned}                                             (13)
```

Simultaneous extra work is $`w+1`$ arbitrary bits, excluding the original
address and q outputs; clean work is zero. The source uses only the w
coefficient wires, and the transport helper is disjoint and returned after
every use. Error and leakage are zero. Every source, control, inverse and
workspace return in (13) is explicit. These are sequential upper bounds,
with no optimality claim for this slow implementation.

### Approximate replacements retain their own source guarantees

For the scalar word, suppose supplied native sources approximate
$`D_\chi`$ and $`P_x`$ within full operator errors
$`\delta_D`$ and $`\delta_P`$, including identity return of all extra
arbitrary helpers. Their actual inverses have the same errors. The complete
word differs from (1) by at most $`2\delta_D+2\delta_P`$.
If it approximates a prescribed Boolean scalar phase with error epsilon and

```math
\varepsilon+2\delta_D+2\delta_P\lt2,
```

then that target must equal $`\chi(a_x)`$, because distinct Boolean
scalar phases have operator distance two. The degree and width conclusions
persist. The XOR version requires the corresponding full accuracy promise
for the **controlled** source $`D_{\mathbf h}`$; it is not inferred from
accuracy of an uncontrolled source. This is a forward source-error ledger,
not a constraint on arbitrary fused circuits.

## 5. Fair comparison and the next mechanism question

For every one-output table, R91 already gives
$`T=O(\sqrt Q)`$, $`G=O(Q)`$ and
$`D_T=O(\log_2(a+2))`$ with $`B=\Theta(\sqrt Q)`$ extra dirty
work. A degree-a output in the present representation needs
$`w\ge Q/2`$ before paying source or transport synthesis. Its explicit
native realization (13) gives no improvement even when that larger width
is available. The same width obstruction applies if any output of a
multi-output table has degree a.

The failed candidate nevertheless changes the constructive diagnosis.
Nonlinear transport can evade the fixed-affine diagonal-source obstruction
and can be made exactly native, with controlled phases and complete dirty
return. The obstruction here is the commutative character's logarithmic
degree in its physical coefficient width, not an assumed inability to
return nonlinear work. The noncommutative cycle construction now supplies
high-degree returned characters at polynomial width, but its paid source
and transport still give no depth improvement. The remaining question is
whether a representation can carry arbitrary tables with cheap native
controlled sources, transport and actual inverses, or whether a different
extraction identity can do so. No dimension-growing unrestricted native
lower bound or R96/R98 improvement follows from either example.

[Focused evidence](../../tests/test_nonlinear_character_transport.py) checks
the character basis, all-input factorization and actual inverse, the degree
bound and unbounded-degree witness, literal returned phases with affine
coefficient maps, and one complete native multi-output word including all
physical columns and its resource ledger. The analytic arguments establish
the arbitrary-size statements; the bounded checks do not extrapolate a
native lower bound.
