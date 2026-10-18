# Fully returned diagonal translation echoes have cubic Boolean phase

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Read-only layer results](READ_ONLY_THREE_LAYER_RIGIDITY.md)

This is an exact, width-independent theorem for a scalar translation echo
with an arbitrary diagonal unitary source. Source phases may have arbitrary
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
with two linear translation families. It also assumes exact scalar return.
No positive approximation radius, query-degree bound at five layers, or general-frame
lower bound is inferred.

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
These checks support, and are not the proof of, the theorem.
