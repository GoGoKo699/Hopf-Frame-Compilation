# Clifford output duals need a large symplectic abelianization

[Complete one-band feasibility](ONE_BAND_COPY_FEASIBILITY.md) ·
[Exact copy grading](DIRTY_COPY_GRADING.md) ·
[Commuting symplectic-orbit bound](CLIFFORD_SEED_COPY_STRUCTURE.md)

A one-band copy of an arbitrary unitary seed has Clifford output duals.
Each such dual must itself have a copyable query orbit. Its generated
symplectic group must retain almost every address label in its abelianization:
for a p-bit field query it needs at least $`2^p-p-1`$ independent binary
characters. This applies even when the dual symplectic matrices do not
commute, and places no depth or Clifford assumption on the original seed.

A separate word certificate gives quantitative exclusions. In particular,
if the symplectic product of two distinct dual orbit elements has odd
order r, then for $`p\ge3`$ the complete copy error is at least
$`1/(8r)`$. Passing these necessary conditions does not provide the full
carrier fixed point, a native seed, or a successful all-size copy. R96
and the conditional 42-band indicator target are unchanged.

## 1. One output dual gives an independently graded Clifford orbit

Let $`p\ge1`$, $`E=\mathbb F_2^p`$, $`m=2^p`$, and let Q be a literal additive
Clifford query representation on the complete carrier A:

```math
Q_xQ_y=Q_{x+y},\qquad Q_0=I.
```

All private helpers are part of A. The query may act trivially on some
of them. Let R be any unitary and $`R_x=Q_xRQ_x^\dagger`$. A proposed
physical copy C acts on A and an arbitrary m-qubit output Y, with defect

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                     -R_x\otimes X_{Y_x}\|.             \qquad(1)
```

Fix an output position t and put

```math
J=C^\dagger Z_{Y_t}C,\qquad
J_x=(Q_x^\dagger\otimes I_Y)J(Q_x\otimes I_Y).          \qquad(2)
```

For a literal one-band $`C=F D G`$, J is a Hermitian Clifford:
the actual inverse word has one physical T/T-dagger band, and sends
Paulis to Cliffords. More generally, the arguments below apply whenever
this particular output dual is Clifford, regardless of C's total depth.

The requested copied operator in (1) has sign $`(-1)^{[x=t]}`$ under
conjugation by $`Z_{Y_t}`$. Subtracting the two compared operators
before and after that conjugation, then pulling back by C and Q, gives

```math
\|(R\otimes I_Y)J_x(R^\dagger\otimes I_Y)
                      -(-1)^{[x=t]}J_x\|\le2\eta.       \qquad(3)
```

This uses only unitarity of R and J. In particular, exact copying
gives one singleton sign character of the dual orbit. Translating its
implementer by every Q gives every other singleton sign character.
For any word in the $`J_x`$ and their actual inverses, a nonzero
occurrence-parity coordinate therefore forces its trace to be zero.
By the [exact grading criterion](DIRTY_COPY_GRADING.md#1-the-exact-full-space-criterion),
the complete Clifford family $`\{J_x\}`$ admits an abstract dirty copy.

This is a necessary dual reduction, not an implementation of that copy.
It must hold separately for every physical output t. It does not
replace the other output probes or the original carrier-return equation.

The same conclusions hold for a Clifford pullback of any physical output
Pauli whose Y-Z label f has odd weight. Its signs are the translates
of f, which span every binary sign table by the
[odd-weight translation argument](NATIVE_COPY_INTERACTION.md#3-dual-support-must-carry-the-full-translation-rank).
Thus every scalar word still has zero parity. For the robust certificates
below, a nonzero parity has nonzero pairing with some translate of f,
so the same word lengths and error constants apply.

## 2. The symplectic group must retain almost all address grades

Let $`F(x)=S(Q_x\otimes I_Y)`$ have entrywise Boolean degree at most d,
where $`d\ge1`$. Define the literal Clifford group and its symplectic image

```math
K=\langle J_x:x\in E\rangle,\qquad
S=S(K),\qquad N=\ker(K\longrightarrow S).
```

The symbol S here denotes a group, not the native one-qubit phase gate.
Its abelianization is $`S_{\rm ab}=S/[S,S]`$. Since every generator
$`S(J_x)`$ is an involution, this abelianization is an elementary
abelian 2-group, and so has a well-defined binary dimension.

**Necessary abelianization theorem.** If (1) holds exactly, then

```math
\dim_{\mathbb F_2}S_{\rm ab}
 \ge 2^p-M_d(p),\qquad
M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj.                \qquad(4)
```

For the field shear, d is one and this becomes

```math
\dim S_{\rm ab}\ge m-p-1.                             \qquad(5)
```

Thus a perfect symplectic orbit group, or more generally one with
binary abelianization rank below this bound, cannot supply the dual
of an exact copy. Noncommutation is permitted, but it cannot erase
the required large abelian quotient. For p equal to one, (5) is
vacuous; no exclusion of the known one-address component follows.

### Word parity and the regular address module

Let $`\mathcal R=\mathbb F_2[E]`$, with basis $`e_x`$. Multiplication
is $`e_xe_y=e_{x+y}`$, and addition records occurrence parity. Every
scalar word in K has nonzero trace, so exact grading forces its parity
to vanish. In particular every identity word has zero parity. Two
words representing the same element consequently have the same parity,
giving a well-defined surjective homomorphism

```math
\pi:K\longrightarrow\mathcal R,\qquad
\pi(J_x)=e_x.                                         \qquad(6)
```

All scalar phases in K lie in its kernel. Conjugation by Q_y translates
the generator index, so (6) is E-equivariant for the regular translation
action on $`\mathcal R`$. Inverses contribute the same binary parity.

### Only the Pauli kernel receives a matrix-degree bound

The kernel N consists of Paulis up to scalar. Its binary label space
$`V=\ell(N)`$ is linear and Q-invariant. Scalars are killed by pi,
so its restriction factors through a linear E-equivariant map

```math
V\longrightarrow W=\pi(N)\subseteq\mathcal R.         \qquad(7)
```

Let I be the augmentation ideal of $`\mathcal R`$. The query degree
implies $`I^{d+1}V=0`$: for any d+1 directions a_i, finite differences
and the representation law give

```math
0=\Delta_{a_1}\cdots\Delta_{a_{d+1}}F(x)
 =F(x)\prod_{i=1}^{d+1}(F(a_i)-I).                    \qquad(8)
```

Therefore $`I^{d+1}W=0`$. For the field query this is also immediate
from $`F(x)=I+N_x`$ and $`N_xN_y=0`$. Pauli conjugation signs do
not alter (7)--(8), because their scalar factors have zero pi image.
No matrix-entry degree is being assigned to the nonabelian group S
or to its abelianization.

Choose address basis a_1,...,a_p and write $`t_i=e_{a_i}+e_0`$.
Then

```math
\mathcal R\cong
 \mathbb F_2[t_1,\ldots,t_p]/(t_1^2,\ldots,t_p^2),
\qquad I=(t_1,\ldots,t_p).
```

Its square-free monomials give

```math
\mathrm{Ann}_{\mathcal R}(I^{d+1})
 =\mathrm{span}\{t_A:|A|\ge p-d\},\qquad
\dim\mathrm{Ann}_{\mathcal R}(I^{d+1})=M_d(p).   \qquad(9)
```

To check necessity for a linear combination, fix a set B of d+1
indices. Multiplication by t_B kills monomials meeting B and maps
the others to distinct monomials $`t_{A\cup B}`$. There can be no
cancellation among them. A term with $`|A|\lt p-d`$ is detected
by choosing B in its complement. All larger terms are annihilated.
When $`d\ge p`$ the ideal power is zero and (9) includes the entire
regular module. For d equal to one its dimension is p+1.

Finally (6) descends to an E-equivariant surjection

```math
S=K/N\longrightarrow\mathcal R/W.                    \qquad(10)
```

The target is elementary abelian, so this map factors through
$`S_{\rm ab}`$. Equations (7)--(10) prove (4). The stronger structural
statement is that $`S_{\rm ab}`$ must have a quotient
$`\mathcal R/W`$ with $`W\subseteq\mathrm{Ann}(I^{d+1})`$.
If the Pauli-kernel label dimension is known, the bound sharpens to
$`\dim S_{\rm ab}\ge m-\min(\dim V,M_d(p))`$.

## 3. Explicit symplectic relations give robust certificates

The rank theorem is exact. A rank failure supplies some relation, but
does not bound its length uniformly in dimension, width or the chosen
symplectic group. A quantitative certificate retains that length.

Let w be a word of length ell in the dual generators and their inverses
whose symplectic product is identity. Its literal unitary value is a
Pauli up to scalar. Let $`a=\pi(w)`$ now denote the formal occurrence
parity of this specified word, without assuming exact copying or a
well-defined group parity map. Suppose d+1 address directions satisfy

```math
\sigma=\prod_{i=1}^{d+1}(1+e_{a_i})a\ne0.             \qquad(11)
```

Take the ordered product of all $`2^{d+1}`$ Q-translates of w over
their cube. The Pauli label of that product is zero by (8), so the
literal value is scalar. Its formal parity is sigma. Translate the
entire word once more so that its t-coordinate is one. The resulting
scalar word has length $`2^{d+1}\ell`$ and odd requested sign in (3).
Unitary telescoping therefore gives

```math
2\le2^{d+1}\ell\,(2\eta),\qquad
\eta\ge\frac1{2^{d+1}\ell}.                          \qquad(12)
```

The scalar need not equal +1: every scalar commutes with R, while
the requested odd sign negates it. Actual inverses and all literal
Clifford phases remain in the word. For the field query, (11) is a
two-direction test and (12) is $`\eta\ge1/(4\ell)`$.

A deficient abelianization guarantees a word satisfying (11): the
parities of symplectic identity words are the kernel of the map from
the free binary generator space to $`S_{\rm ab}`$. Indeed, a word
landing in $`[S,S]`$ can be corrected to identity by commutator lifts,
whose occurrence parity is zero. If that kernel
has dimension greater than M_d(p), it cannot lie in the annihilator
(9). This establishes existence of a finite certificate, without
claiming a uniform bound on its word length.

### Odd-order pair products are inexpensive certificates

Suppose $`z\ne0`$ and $`S(J_0)S(J_z)`$ has odd order r. Then

```math
w=(J_0J_z)^r,\qquad \ell=2r,\qquad
S(w)=I,\qquad \pi(w)=e_0+e_z.                         \qquad(13)
```

For the field query and $`p\ge3`$, choose a,b such that z,a,b are
independent. The parity in (11) is
$`(1+e_a)(1+e_b)(e_0+e_z)`$, supported on eight distinct addresses.
Consequently

```math
\eta\ge\frac1{8r}.                                   \qquad(14)
```

More generally a degree-d query admits this argument when
$`p\ge d+2`$, with floor $`1/(2^{d+2}r)`$.
For a direct copy of the Clifford family itself, rather than its use
as a dual in (3), the per-generator error is eta instead of 2 eta;
the corresponding floors in (12) and (14) double. A known shorter
literal scalar relation can of course give a stronger certificate.

This screen includes genuinely noncommuting symplectic pairs. For
example, let J be H on one u/v core pair, times $`Z_{Y_t}`$.
It is an actual one-band output dual: with

```math
V=SHTHS^\dagger,\quad VZV^\dagger=H,\qquad
B=\mathrm{CX}(u;Y_t)\mathrm{CX}(v;Y_t),\qquad
C=B(V_u^\dagger V_v^\dagger),
```

we have $`C^\dagger Z_{Y_t}C=H_uH_vZ_{Y_t}`$. C uses two
T-dagger gates in one physical band. The field address z whose
matrix is identity has Q_z equal to pairwise CNOTs. On this chosen
pair, $`J_0J_z`$ is the product of the two opposite CNOTs and has
order three. Its symplectic generators do not commute. The literal
pair word is already identity here, so this example has the stronger
dual floor $`1/6`$ wherever those two addresses are present. It is a
fixture of a failed one-band candidate, not a proposed complete copy.

### The retained commuting-orbit certificate also dualizes

If all symplectic matrices of the dual orbit commute, the
[scalar-word theorem](CLIFFORD_SEED_COPY_STRUCTURE.md#3-a-second-cube-gives-a-scalar-relation-with-nonzero-copy-parity)
supplies a word of $`2^{5d+2}`$ distinct addresses for
$`p\ge5d+2`$. Translate its support to include t and apply (3).
Thus

```math
\eta\ge2^{-5d-2},\qquad
\eta\ge1/128\quad\hbox{for the field query and }p\ge7.
```

The cited proof owns that word construction. The dual application
allows an arbitrary original unitary seed; it does not assert that
all one-band dual orbits have commuting symplectic actions.

## 4. Consequence for the constructive search

A changed one-band wall family needs more than noncommuting dual
matrices. Each output dual must avoid short forbidden-parity scalar
relations and must retain the large abelianized symplectic quotient
(4). Even then, all duals together and the complete carrier fixed
point must be satisfied by the same native seed and actual inverse.

The [exact fixtures](../../tests/test_clifford_dual_copy_grading.py)
check finite regular-module annihilators, the affine field action on
Pauli grades, a noncommuting order-three symplectic pair and its literal
one-band dual, negative scalar phases, and the eight-address parity
certificate. These finite checks support the algebraic identities;
the proofs above establish their all-size scope. No unrestricted
one-band exclusion, uniform gap from rank alone, or improved complete
frame-depth theorem follows.
