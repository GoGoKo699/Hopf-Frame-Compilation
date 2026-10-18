# Fully returned diagonal translation echoes have cubic Boolean phase

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Read-only layer results](READ_ONLY_THREE_LAYER_RIGIDITY.md)

This is a width-independent theorem for a scalar translation echo with an
arbitrary diagonal unitary source, with a fixed positive error radius for
approximation to a prescribed Boolean scalar phase. Source phases may have arbitrary
precision, including irrational multiples of pi; they need not be Boolean
signs. It does not constrain non-diagonal fixed conjugators or arbitrary
products of four Pauli representations, and it does not prove a T-depth
recurrence. The sharp native example below is a scope witness, not a resource
improvement.

## Statement

Let $`h:\mathbb F_2^w\to\mathbb R/\mathbb Z`$ be arbitrary and let
$`L,M:\mathbb F_2^a\to\mathbb F_2^w`$ be linear. Write

```math
D_h|z\rangle=e^{2\pi i h(z)}|z\rangle,\qquad
A_x=X^{Lx},\qquad B_x=X^{Mx}.
```

Define the literal commutator

```math
U_x=A_x(D_hB_xD_h^\dagger)A_x^\dagger
                     (D_hB_x^\dagger D_h^\dagger).
```

Every work input is arbitrary. Suppose for every x the complete work operator
is a scalar phase,

```math
U_x=\lambda(x)I_W,\qquad |\lambda(x)|=1.
```

Then the phase is necessarily Boolean: $`\lambda(x)=(-1)^{q(x)}`$ for a
Boolean q satisfying $`q(0)=0`$ and $`\deg_{\mathbb F_2}q\le3`$.
This bound is independent of address length, dirty width, phase precision and
the complexity of h. The conclusion is exact and includes complete work
return with arbitrary superpositions and external references.

For the additive difference $`\Delta_vh(z)=h(z+v)-h(z)`$ in
$`\mathbb R/\mathbb Z`$, direct action on every basis column gives

```math
U_x|z\rangle
=e^{-2\pi i\Delta_{Lx}\Delta_{Mx}h(z)}|z\rangle.
```

The minus sign follows from the actual right-to-left order of the displayed
commutator; the two source calls have positive phases and both actual inverse
calls have negative phases. Scalar return is therefore equivalent to
$`\Delta_{Lx}\Delta_{Mx}h(z)=c(x)`$ independent of z.
Translation by Lx negates this difference, since Lx has order at most two.
Thus $`2c(x)=0`$ in $`\mathbb R/\mathbb Z`$, including the case Lx is zero,
and uniquely $`c(x)=q(x)/2`$ with Boolean q. The sign in the exponential then
makes no difference. No initialized state, measurement, scalar convention or
intermediate work return is used.

## Proof in the additive derivative algebra

All operators below act on the abelian group of
$`\mathbb R/\mathbb Z`$-valued functions of z. For the standard address basis
$`e_i`$ define

```math
a_i=\Delta_{Le_i},\qquad b_i=\Delta_{Me_i}.
```

They commute and satisfy $`a_i^2=-2a_i`$ and $`b_i^2=-2b_i`$.
They are not assumed square-zero. The argument below only uses commutativity
and the scalar-return equations, and remains valid for overlapping images of
L and M or linearly dependent images. For sets A and B write
$`a_A=\prod_{i\in A}a_i`$ and $`b_B=\prod_{i\in B}b_i`$, with empty product I.
Expanding the translations as products of $`I+a_i`$ and $`I+b_i`$ gives the
unique multilinear expansion in the Boolean address coordinates:

```math
\Delta_{Lx}\Delta_{Mx}
  =\sum_{\varnothing\ne S\subseteq[a]}x^S E_S,
\qquad
E_S=\sum_{\substack{\varnothing\ne A,B\subseteq S\\A\cup B=S}}
                    a_A b_B.
```

The union allows A and B to overlap: $`x_i^2=x_i`$ on Boolean inputs.
Multilinear coefficients exist over any abelian group by signed Moebius
inversion. Scalar return says each $`E_Sh`$ is a constant function of z.
Consequently every nonempty derivative word R in the $`a_i,b_i`$ obeys

```math
R E_S h=0.                                                   (1)
```

First, if $`A\cap B`$ is nonempty and $`|A\cup B|\ge2`$, then
$`a_Ab_Bh=0`$. Choose i in the intersection. Since
$`E_{\{i\}}=a_ib_i`$, the monomial is $`E_{\{i\}}`$ multiplied by a nonempty
remaining derivative word, and (1) applies. This disposes of all
terms with overlapping indices below.

Fix a set S of $`n\ge3`$ distinct indices. For nonempty proper A subset S put

```math
v_A=a_A b_{S\setminus A}h.
```

For $`i\in A`$ and $`j\in S\setminus A`$ multiply $`E_{\{i,j\}}`$ by

```math
R=a_{A\setminus\{i\}}b_{S\setminus(A\cup\{j\})}.
```

R is nonempty because $`n\ge3`$. The only terms in $`E_{\{i,j\}}`$ without
an overlapping index are $`a_ib_j`$ and $`a_jb_i`$; its other terms vanish
by the preceding paragraph after multiplication by R. Equation (1) gives

```math
v_A=-v_{(A\setminus\{i\})\cup\{j\}}.                    (2)
```

The minus sign matters for non-Boolean sources. The graph of k-element
subsets joined by one exchange is connected and contains a triangle for
$`1\le k\le n-1`$ and $`n\ge3`$. For $`k\le n-2`$, fix k-1 members and
choose the remaining member in three ways; for $`k=n-1`$, use complements
of three singletons. Applying (2) around a triangle gives $`v_A=-v_A`$.
Connectivity propagates this identity, so every $`v_A`$ is two-torsion and
all values with $`|A|=k`$ agree. Denote their common value by $`v_k`$.

Now let $`n\ge4`$. Choose three indices T subset S. For each integer
$`0\le r\le n-3`$, choose R to contain one derivative on every index of
$`S\setminus T`$, with r of them a derivatives and $`n-3-r`$ of them b
derivatives. R is nonempty. In $`RE_Th=0`$ all terms with overlapping indices
vanish. There are three remaining terms having one a derivative on T and
three having two. Hence

```math
0=3v_{r+1}+3v_{r+2}=v_{r+1}+v_{r+2}.
```

The last equality uses two-torsion. Thus all $`v_1,\ldots,v_{n-1}`$ agree.
The coefficient $`E_Sh`$ is the sum of one $`v_A`$ for each nonempty proper
subset A, with overlapping terms zero:

```math
E_Sh=\sum_{\varnothing\ne A\subsetneq S}v_{|A|}
     =(2^n-2)v_1=0.
```

Every additive multilinear coefficient of c of degree at least four therefore
vanishes. The identification $`c=q/2`$ identifies this two-torsion subgroup
with $`\mathbb F_2`$ and its signed coefficients with the Boolean ANF
coefficients of q. This proves the theorem without a rank, sparsity,
phase-denominator or finite-width constant.

## Robust return at a width-independent error radius

The same degree conclusion survives a nonzero complete operator error. Let
f be a prescribed Boolean function and retain the literal echo above. If

```math
\max_x\|U_x-(-1)^{f(x)}I_W\|_{\rm op}\le\eta,
\qquad \eta\lt\frac{2}{275},
```

then $`f(0)=0`$ and $`\deg f\le3`$. The constant is a convenient analytic
bound, with no optimality claim. It is independent of work width, address
length, source precision and the size of a native source implementation.
All work inputs, including correlations with external references, remain
covered. No assumption of exact scalar return, an initialized work state,
or a rational phase denominator is used. The echo is diagonal on work and
block diagonal in the unchanged address, so the full operator error is the
maximum of these phase errors. Tensoring with an external reference does not
change that norm.

Here is a bounded four-direction proof. For an additive phase function v,
write, with $`T_u`$ denoting translation by u,

```math
\nu(v)=\sup_z|e^{2\pi i v(z)}-1|.
```

This quantity is subadditive and invariant under translation. If v is
within eta of a constant in this phase metric, applying a nonempty product
of k translation differences gives phase norm at most $`2^k\eta`$.
More generally, an integer translation polynomial P with zero coefficient
sum gives the bound $`\|P\|_1\eta`$, where $`\|P\|_1`$ is the sum of
absolute coefficients. These follow by multiplying unit complex numbers
and telescoping; they do not require choosing a real logarithm of a phase.

Restrict the address to any four directions, allowing dependence and zero
directions. Use their coordinates 1 through 4, and retain the operators
$`a_i,b_i,E_S`$ from the exact proof. Put

```math
d_x=\Delta_{Lx}\Delta_{Mx}h,\qquad d_i=a_i b_i h,
\qquad c_i=a_i+b_i+a_i b_i=T_{Le_i+Me_i}-I.
```

The target comparison is $`\nu(d_x-f(x)/2)\le\eta`$: reversing the
literal negative sign of the echo preserves this bound because the target
phase is Boolean. Also $`U_0=I`$, so $`f(0)=0`$ at the stated error.
Signed Moebius inversion gives

```math
\nu\!\left(E_Sh-\sum_{\varnothing\ne X\subseteq S}
                    (-1)^{|S|-|X|}\frac{f(X)}2\right)
\le(2^{|S|}-1)\eta.                                      (3)
```

For the four-direction coefficient, separate the disjoint terms
$`v_A=a_A b_{[4]\setminus A}h`$ from the terms with an index present in
both factors. Their overlap sum is exactly

```math
O h=\left[\prod_{i=1}^4c_i-\prod_{i=1}^4(a_i+b_i)\right]h
 =\sum_{i=1}^4
   \left[\prod_{j\lt i}(a_j+b_j)\prod_{j\gt i}c_j\right]d_i.
```

Every bracket annihilates constants. Since $`\|a_j+b_j\|_1\le4`$
and $`\|c_j\|_1\le2`$, the full-input comparison for each d_i gives

```math
\nu(Oh)\le(8+16+32+64)\eta=120\eta.                    (4)
```

Next pair two proper nonempty subsets A and B of [4] that differ by
exchanging i in A with j outside A. Their monomial sum is R times
$`a_i b_j+a_j b_i`$, where R contains the two remaining differences.
The exact translation identity is

```math
(a_i b_j+a_j b_i)h
 =d_{e_i+e_j}-T_{Le_j+Me_j}d_i
             -(T_{Le_i}+T_{Me_i}-I)d_j.                  (5)
```

The nonempty R kills each constant comparison phase and has coefficient
norm at most four. Equation (5) therefore gives
$`\nu(v_A+v_B)\le4(1+1+3)\eta=20\eta`$.
The fourteen proper nonempty subsets admit seven such pairs:

| Subset size | Pairs, with concatenated indices denoting a subset |
|---|---|
| 1 | (1,2), (3,4) |
| 2 | (12,13), (14,24), (23,34) |
| 3 | (234,134), (124,123) |

Consequently $`\nu(E_{[4]}h)\le120\eta+7(20\eta)=260\eta`$.
By (3), its distance from the Boolean fourth-difference target is at most
$`15\eta`$. A nonzero fourth difference would require
$`2\le275\eta`$, which contradicts the hypothesis.

Thus every four-direction Boolean derivative of f at zero vanishes. This
also excludes all higher algebraic degrees, not only a degree-four ANF
coefficient on one chosen coordinate set. In characteristic two the
Boolean difference operators obey

```math
\delta_u\delta_v=\delta_{u+v}+\delta_u+\delta_v.
```

A derivative of order greater than four reduces by this identity to sums
of four-direction derivatives at zero. In particular every ANF coefficient
of degree at least four is zero. This completes the robust proof without a
width-dependent continuity estimate or intermediate work-return assumption.

The relation is an analytic constraint on the actual echo. It adds no
physical calls, clean work or dirty work to a construction. If the source
is implemented natively, all source appearances, actual inverses and
translations retain their own full resource charges; the theorem does not
supply such an implementation for an arbitrary U(1) source.


### Approximate native sources remain charged

Suppose an actual fixed native source $`\widetilde D`$ approximates
$`D_h\otimes I`$ within delta in full operator norm on all its physical
work wires. Any additional source wires are arbitrary, included in that
norm, and included simultaneously in its allocation; an initialized-only
promise is not substituted here. Use the actual inverse of this same native
word at both dagger positions. Its inverse error is also at most delta.
Four-term unitary telescoping gives

```math
\|\widetilde U-U\|_{\rm op}\le4\delta.
```

Thus if the completed native echo approximates $`(-1)^{f(x)}I`$ with
full operator error epsilon and

```math
\varepsilon+4\delta\lt\frac{2}{275},
```

then f still has degree at most three. The comparison includes every source
input, every actual inverse, all returned work and all reference systems.
There is no reset, measurement or uncharged recovery.

If one source word has counts $`T_D,G_D`$, T-depth $`d_D`$, and total
physical arbitrary-work allocation W, the displayed echo has the upper
ledger

```math
T\le4T_D,\qquad D_T\le4d_D,\qquad
G\le4G_D+2\mathrm{nnz}(L)+2\mathrm{nnz}(M).
```

Here nnz counts entries equal to one in the displayed binary translation
matrices: one CNOT per entry implements a translation. Both translation
families occur twice, and all four source calls occur chronologically on
the same W wires. The simultaneous allocation is W arbitrary work and zero
clean work, excluding the read-only address. These are charged upper bounds
for this word, with no count or depth improvement over a simplified native
query baseline asserted.

## Sharp cubic example, and a literal native ledger

Use six arbitrary work bits $`p_1,p_2,p_3,q_1,q_2,q_3`$, and three read-only
address bits $`x_1,x_2,x_3`$. Take

```math
h_2(p,q)=p_1q_2q_3+p_2q_1q_3+p_3q_1q_2\quad\text{in }\mathbb F_2,
\qquad h=h_2/2,\qquad Lx=(x,0),\qquad Mx=(0,x).
```

Here $`h=h_2/2`$ uses the embedding of $`\mathbb F_2`$ into the two-torsion
subgroup of $`\mathbb R/\mathbb Z`$. With Boolean differences,

```math
\Delta_{Lx}\Delta_{Mx}h_2
 =\sum_{\mathrm{cyc}(i,j,k)}
 x_i(x_jq_k+x_kq_j+x_jx_k)
 =x_1x_2x_3.
```

The sum ranges over the three cyclic permutations of (1,2,3).
Every term linear in q occurs twice. The three cubic terms add to one.
Hence $`U_x`$ is the literal scalar $`(-1)^{x_1x_2x_3}`$ on all six work bits.
The proof checks every physical input, rather than only a prepared code.

For a completely explicit finite native implementation, let $`D_h`$ be the
product of the three CCZ gates on $`(p_i,q_j,q_k)`$. One CCZ can be built from

```math
4abc\equiv a+b+c-(a\oplus b)-(a\oplus c)-(b\oplus c)
                         +(a\oplus b\oplus c)\pmod8.
```

Apply T to a,b,c in parallel. Implement each negative pair parity by a
CNOT into one member, T-dagger there, and the inverse CNOT. Implement the
positive triple parity by two CNOTs into c, T there, and those two CNOTs
in reverse order. This exact word has 7 T/T-dagger gates, 10 elementary
CNOTs, five declared physical T layers, no helper and no residual phase.
It realizes the displayed mod-8 identity using $`T=\mathrm{diag}(1,e^{i\pi/4})`$.

Execute the three CCZ words successively. $`D_h`$ has 21 T gates, 30 elementary
Cliffords and a schedule of at most 15 physical T layers. Use the actual
reversed inverse word for every $`D_h^\dagger`$ call. The displayed commutator
has two $`D_h`$ calls and two actual $`D_h^\dagger`$ calls, plus two A and two B
words of three address-controlled CNOTs each. Its literal ledger is

| Resource | Displayed native word |
|---|---:|
| T/T-dagger count | 84 |
| Elementary Clifford count | 132 CNOTs |
| T-depth | at most 60 |
| Clean work | 0 |
| Simultaneous arbitrary work | 6 bits |
| Error / leakage | 0 / 0 |
| Work return | Exact identity on all six bits |

There are no additional outputs or helpers hidden in that allocation.
All T gates are on work; every address interaction is a CNOT control.
The interlayer address-preservation requirement is therefore literal.
The counts describe this deliberately simple explicit word, not an optimized
compiler. No saving compared with existing native scalar/query baselines is
claimed. Its role is to certify that the cubic upper degree is attained within
the full native/operator contract. The abstract four-representation form
must not be misreported as four physical T layers.

## Scope and precise unresolved step

This closes the proposed mechanism of hiding arbitrarily high degree in an
arbitrary diagonal $`D_h`$ and trying to expose it with one fully returned
double-translation commutator. Its resulting phase is necessarily Boolean
and remains cubic even when the source has arbitrary U(1) phases and its
native cost is temporarily ignored. Charging a native $`D_h`$ can only add costs.

The theorem does not address non-diagonal fixed conjugators or arbitrary
products of four representation factors. It applies to the displayed echo
with two linear translation families. The robust statement covers approximation to a prescribed Boolean scalar
phase below $`2/275`$ in full operator norm. It gives no query-degree bound
at five layers or general-frame lower bound.

The [bounded evidence checks](../../tests/test_returned_diagonal_echo_degree.py)
include the complete native word and actual inverse on every physical basis
input, the derivative expansion with overlapping and dependent maps, and
linear-algebra elimination for disjoint L,M with $`a=1,\ldots,4`$.
The elimination imposes every full-work return equation on arbitrary Boolean
truth tables h and q, rather than sampling h. In that Boolean calculation,
the maximum possible q degrees are 1,2,3,3; all degree-at-least-four coefficient
functionals belong to the equation row span. Further checks retain signed
additive derivatives, actual inverse phases, overlapping and dependent maps,
and a returned cubic source with non-Boolean phases of several denominators.
Further exact integer group-ring checks verify the bounded four-direction
identity and its error coefficient, including overlapping and dependent
translations. A five-address-bit kernel check verifies that four-direction
tests through zero exclude fifth degree as well as fourth degree. A perturbed
non-Boolean source exercises positive full-input error and the literal
inverse order. These checks support, and are not the proof of, the theorem.
