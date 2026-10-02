# Exact depth within the Majorana-linear source model

[Operator source](OPERATOR_SOURCE_COMPILER.md) · [Paired source](ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks) · [Full-frame depth](T_DEPTH_COMPILER.md)

The geometric sources have exact optimal T-depth **within a specified
Majorana-linear circuit class**. The paired source also admits a parallel
tail schedule with no additional work qubits. These results concern exact
source primitives; they neither establish unrestricted Clifford+T depth
optimality nor change the retained Hopf asymptotic bounds.

## 1. Gate class and inherited denominator bound

On every circuit wire, including any additional dirty wires, fix the
paired Jordan–Wigner Majoranas

```math
\gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j.
```

An admissible schedule alternates physical T/T-dagger layers on distinct
wires with Clifford blocks whose **complete conjugation action** is a
signed permutation of these Majoranas. Elementary Clifford gates inside
a block need not separately have that property. Arbitrary signed
permutations are allowed, with no connectivity restriction; their native
Clifford gates still contribute to G. A T layer rotates disjoint on-site
Majorana pairs by angles $`\pm\pi/4`$. Conjugating such a layer by an
allowed Clifford gives disjoint Majorana-plane rotations.

This is a discrete extension of Clifford-matchgate+T circuits by
parity-odd Clifford stages, lying in the Pin group up to global phase.
In particular, a source
$`M=\sum_i a_i\gamma_i`$ is parity odd: its conjugation matrix is
$`2aa^{\mathsf T}-I`$, with determinant minus one. Calling its synthesis
strictly even matchgate synthesis would exclude the source itself.

For entries in $`\mathbb Z[1/\sqrt2]`$, define

```math
\mathrm{sde}(x)=\min\{k\ge0:(\sqrt2)^k x\in\mathbb Z[\sqrt2]\}.
```

Take the maximum over matrix entries or vector coordinates.
Casas, Braccia, Gouzien, Cerezo, and García-Martín,
[*Matchgate synthesis via Clifford matchgates and T gates*,
Section III.2.2, equation (29)](https://arxiv.org/html/2602.05425v1#S3.SS2.SSS2),
give the denominator lower bound for even matchgate circuits: signed
permutation Clifford matrices do not increase this exponent, and one
T layer increases it by at most one. The same proof applies when
determinant-minus-one signed permutations are also allowed. Indeed,
each row of a T-layer matrix is either a coordinate row or the sum or
difference of two coordinate rows divided by $`\sqrt2`$; addition cannot
increase the maximum denominator exponent.

Consequently a depth-D admissible schedule has Majorana conjugation
matrix with $`\mathrm{sde}\le D`$. The same bound holds for the
coefficient vector obtained by conjugating one Majorana. Extra arbitrary
dirty wires cannot evade it: the desired core coefficients remain
entries of the enlarged conjugation matrix. This permits arbitrary
mixing with those wires at intermediate times and requires their actual
return for the completed source. No source state is initialized.

This is an application and Pin extension of the cited denominator
argument, not a new general depth lower-bound method. Under unrestricted
Cliffords a Pauli can spread across many physical T targets in one layer;
the Majorana denominator argument then does not apply.

## 2. Geometric source and loader

Loader exactness here means its induced Majorana transformation; a
native representative may differ from a displayed exponential by a
common scalar. Source exactness always means the literal full unitary.

The [original source](OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit)
uses $`\Gamma_j=\gamma_{2j}`$ and coefficients

```math
a_j=2^{-(j+1)/2}\quad(j\lt m-1),\qquad
a_{m-1}=2^{-(m-1)/2},\qquad m\ge2.
```

Any admissible loader taking $`\gamma_0`$ to $`M_m=\sum_j a_j\Gamma_j`$
needs depth at least $`m-1`$, from the last coefficient. The chronological
planes $`(\gamma_0,\gamma_2),(\gamma_2,\gamma_4),\ldots`$, each with
$`R_{uv}=\exp(i\pi(i\gamma_u\gamma_v)/8)`$, attain this depth.

The source itself has exact optimal depth $`2m-4`$ in the same class.
For $`m\ge3`$ its last diagonal conjugation coefficient is

```math
2a_{m-1}^2-1=2^{2-m}-1,
\qquad \mathrm{sde}=2m-4.
```

For the matching upper bound absorb the first split into the literal
Clifford $`M_2=(\Gamma_0+\Gamma_1)/\sqrt2`$, as in the
[exact source-count proof](OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers).
If $`V=R_{m-2}\cdots R_1`$, then $`M_m=V M_2V^\dagger`$.
Chronologically apply the actual inverse of V, then $`M_2`$, then V.
There are $`2(m-2)`$ T layers and $`O(m)`$ Clifford gates, with no
helper. At $`m=2`$ the source is already Clifford.

The individual even-index plane has a constant-size implementation in
the admissible class. Let $`D_j=\exp(-i\pi X_jX_{j+1}/4)`$.
Its conjugation sends $`Y_jX_{j+1}`$ to $`Z_j`$, so the chronological
word $`D_j,T_j^\dagger,D_j^\dagger`$ gives the desired plane up to a
common scalar. Each complete $`D_j`$ is a signed Majorana permutation.
All such loader scalars cancel in the actual source and controlled-source
conjugations; no relative branch phase is discarded.

## 3. Two tails in parallel

For the [paired source](ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks),
put $`m=q+1`$, $`q\ge2`$. Its head coefficient is $`1/2`$ and its
two tails have coefficients $`\sqrt{w_j/2}`$ and $`\sqrt{w_j/4}`$,
where $`w_j=2^{-j-1}`$ except for the duplicated final weight
$`w_{q-1}=2^{-(q-1)}`$. The unused final Majorana has zero coefficient.
The following list is chronological:

1. Apply $`R_{01}`$.
2. Apply $`R_{02}`$.
3. For $`j=0,\ldots,q-2`$, apply together
   $`R_{2j+1,2j+3}`$ and $`R_{2j+2,2j+4}`$.

All odd-tail planes commute with all even-tail planes, because their
Majorana indices are disjoint. Thus this interleaves the two old tails
without changing the unitary or its phase. Each pair in step 3 uses only
physical wires $`j,j+1,j+2`$, despite their shared middle wire.

Here is an explicit one-T-layer realization of that pair. On these three
wires relabel the six local Majoranas as $`\delta_0,\ldots,\delta_5`$;
the common preceding Z string cancels from all bilinears. Define the
Clifford plane swap
$`G_{uv}=\exp(i\pi(i\delta_u\delta_v)/4)`$ and apply the
chronological Clifford word

```math
G_{01},\ G_{13},\ G_{34},
\qquad D=G_{34}G_{13}G_{01}.
```

Conjugation by D sends
$`(\delta_1,\delta_3,\delta_2,\delta_4)`$ to
$`(-\delta_0,-\delta_1,\delta_2,-\delta_3)`$. Hence it sends
$`i\delta_1\delta_3=-X_jY_{j+1}`$ to $`-Z_j`$ and
$`i\delta_2\delta_4=Y_{j+1}X_{j+2}`$ to $`Z_{j+1}`$.
The desired paired block is therefore implemented chronologically by

```math
D,\quad T_jT_{j+1}^\dagger,\quad D^\dagger.
```

The two native T phases cancel each other's common scalar, so this block
is literally the desired pair. Each G is a weight-one or weight-two
Pauli $`\pi/4`$ rotation, hence a constant-size Clifford word. D needs
no helper and, as a complete block, preserves the Majorana span.

The loader now has $`D_T=q+1`$, $`T=2q`$, and $`G=O(q)`$.
Its last even-tail coefficient is $`2^{-(q+1)/2}`$, so the denominator
bound proves this depth is optimal in the specified class. The physical
seed implementations are $`T_0`$ followed by the positive
$`Y_0X_1`$ plane; their common scalars also cancel.

For the source, absorb the first seed into

```math
C_0=\frac{X_0+Y_0}{\sqrt2}=\omega^{-1}S_0X_0,
\qquad \omega=e^{i\pi/4}.
```

This is a literal Clifford: $`(HS)^3=\omega I`$ supplies its inverse
scalar using Clifford gates. It swaps $`\gamma_0,\gamma_1`$ and
negates the remaining Majoranas, so it belongs to the Pin class.
Let V be steps 2–3 above. Chronologically execute
$`V^\dagger,C_0,V`$. This gives the exact source with

```math
D_T=2q,\qquad T=4q-2,\qquad G=O(q),
```

and no additional wire. Its last even-tail diagonal conjugation entry is
$`2^{-q}-1`$, with denominator exponent $`2q`$, proving matching
restricted depth. The T-count here is an upper bound, not an asserted
new count minimum.

## 4. Scope and finite evidence

| Exact primitive | Optimal depth in the specified class |
|---|---|
| Geometric loader, $`m\ge2`$ | $`m-1`$ |
| Geometric source | $`2m-4`$ |
| Paired loader, $`q\ge2`$ | $`q+1`$ |
| Paired source | $`2q`$ |

The source equalities hold on the entire dirty Hilbert space and hence
with arbitrary references. An uncontrolled loader scalar is never
substituted for a controlled-source phase. Using an actual loader and
inverse around the central controlled Pauli still gives $`O(q)`$
controlled-source depth, but that central controlled Clifford need not
preserve the Majorana span. No controlled-source depth lower bound is
claimed here.

Arbitrary Clifford+T resynthesis, alternative approximate sources, and
joint resynthesis of multiple source calls remain outside this lower
bound. Exact denominator certificates do not establish approximate
depth lower bounds. In particular, source costs cannot be added to
deduce an unrestricted complete-frame or state-preparation lower bound.
The full-frame and state-based Hopf depth orders are unchanged.

[Finite checks](../tests/test_source_t_depth.py) compare native schedules,
their actual inverses, literal source phases, and all dirty-input columns.
They verify signed-permutation Clifford blocks and distinct physical T
targets, and use exact rational $`\mathbb Q(\sqrt2)`$ arithmetic for
the coefficient and denominator witnesses. These bounded checks support
the explicit schedules; the uniform optimality statements use the proof
above.
