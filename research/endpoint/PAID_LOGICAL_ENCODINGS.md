# Paid logical encodings and complete sector controls

[Endpoint](../../docs/OPEN_PROBLEM.md) · [Disjoint native blocks](DISJOINT_BLOCK_COMPILER.md) · [Structural limits](STRUCTURAL_COMPILATION_LIMITS.md) · [Center-phase reflections](CENTERED_PHASE_REFLECTIONS.md)

A nonlinear basis permutation can compress vacuum/one-hot data at an exact
cost linear in the number of data bits. Its tensor product handles a
branching support of exponential size, and its complete encoded Pauli
controls are Clifford. These are explicit representations of families
already covered abstractly by the disjoint-block theorem, not a stronger
unrestricted Hopf compiler. Gates and inverses below retain literal phase.
The two-clean, arbitrary-dirty/reference contract and endpoint error
$`\eta=2^{-N}`$, $`N=2^n`$, remain unchanged.

The implementation is
[paid_logical_encodings.py](../../compiler_robust_hopf/paid_logical_encodings.py).
Its gate words use chronological order and the least significant integer
bit as wire zero. It reuses the retained literal Toffoli and actual-inverse
helpers; no arbitrary-angle fine core is emitted.

## 1. A bit-optimal vacuum/one-hot permutation

For $`m\ge1`$, put $`r=\lceil\log_2(m+1)\rceil`$. Label the m data wires
by the nonempty binary subsets $`A=1,\ldots,m`$ of an r-element set.
This set of subsets is downward closed. An exact reversible permutation
on these same wires satisfies

```math
E_m|0^m\rangle=|0^r\rangle|0^{m-r}\rangle,\qquad
E_m|e_j\rangle=|j\rangle|0^{m-r}\rangle\quad(1\le j\le m).
```

The r output bits are the original singleton wires; a fixed wire permutation
can put them contiguously at zero T-count. The other inputs have a specified
reversible image, not an erasure or a promise of zero remainder.

First compute the upper subset transform using CNOTs:

```math
y_A=\bigoplus_{J\supseteq A,\ J\le m}x_J.
```

For each bit i, update $`y_A\leftarrow y_A\mathbin\oplus y_{A\cup\{i\}}`$
when i is absent from A and the larger label exists. The successive bit
passes give the displayed transform. On a one-hot input at J,
$`y_A=\mathbf1_{A\subseteq J}=\prod_{i\in A}y_{\{i\}}`$.
Now visit nonsingleton A in decreasing cardinality, choose i in A, and use
one Toffoli to perform

```math
y_A\longleftarrow y_A\mathbin\oplus y_{\{i\}}y_{A\setminus\{i\}}.
```

Its lower-degree factors have not yet been cleared. Thus the nonsingletons
become zero on the stated code, and the singleton bits encode J; vacuum
stays zero. The physical word is invertible on every input. Its exact counts are

```math
\#\mathrm{CCX}=m-r,\qquad
\#\mathrm{CNOT}=\sum_{J=1}^{m}\mathrm{popcount}(J)-r,
\qquad T\le7(m-r),\qquad \text{helpers}=0.
```

Two optimality statements have different scopes. The bit count is minimal
because the code has m+1 orthogonal states. In the exact X/CNOT/Toffoli
model, restrict the wire functions to vacuum and the m one-hot inputs and
include the constant-one function. Their initial binary span has dimension
m+1; the required final span has dimension r+1. X and CNOT preserve that
span, while replacing one wire by a Toffoli output reduces its dimension
by at most one. Hence m-r Toffolis are necessary. This is not T-optimality
for general Clifford+T circuits.

## 2. Six-bit joint synthesis and full encoded controls

In semantic order $`(a,b,c,u,v,w)`$, the full six-bit permutation is

```math
\begin{aligned}
a&=x_1\oplus x_3\oplus x_5,& b&=x_2\oplus x_3\oplus x_6,
&c&=x_4\oplus x_5\oplus x_6,\\
u&=x_3\oplus ab,&v&=x_5\oplus ac,&w&=x_6\oplus bc.
\end{aligned}
```

Its physical output order is $`(a,b,u,c,v,w)`$. The three nonlinear clears
commute. Conjugating their targets by H gives phase
$`(-1)^{abu+acv+bcw}`$. Use the exact mod-eight identity

```math
4xyz=x+y+z-(x\oplus y)-(x\oplus z)-(y\oplus z)
 +(x\oplus y\oplus z)\pmod8.
```

The three repeated linear phases on a,b,c become S gates. The remaining
odd coefficients give fifteen T/T-dagger gates. With all parity
computation and uncomputation, the complete encoder has 60 elementary
gates, fifteen T/T-dagger gates, and no helpers. Fifteen is an upper
bound, not an asserted global minimum. The module emits this word rather
than storing a copied gate-list artifact.

The inverse formulas hold on all 64 inputs:

```math
\begin{aligned}
x_1&=a\oplus u\oplus v\oplus ab\oplus ac,&
x_2&=b\oplus u\oplus w\oplus ab\oplus bc,\\
x_3&=u\oplus ab,&x_4&=c\oplus v\oplus w\oplus ac\oplus bc,\\
x_5&=v\oplus ac,&x_6&=w\oplus bc.
\end{aligned}
```

Consequently every $`Z'_j=E_6Z_jE_6^\dagger`$ is a product of Z and CZ.
For example $`Z'_1=Z_aZ_uZ_v\mathrm{CZ}_{a,b}\mathrm{CZ}_{a,c}`$ and
$`Z'_3=Z_u\mathrm{CZ}_{a,b}`$. Translations are affine: the image of
$`X_1`$ toggles a, updates $`u\leftarrow u\oplus b`$ and
$`v\leftarrow v\oplus c`$, and fixes the other bits. The implementation
supplies every translated X word and $`Y'_j=iX'_jZ'_j`$, including its
literal scalar i, entirely in H/S/CNOT. No inference about arbitrary
permutations is made merely from quadratic Boolean degree.

For an original conditional rotation let

```math
\Pi'=\prod_{l\ne j,\ {\rm specified}}
 \frac{I+(-1)^{c_l}Z'_l}{2}.
```

The commuting transported factors give the whole operator identity

```math
E_6e^{-i\theta\Pi Y_j}E_6^\dagger
=I+(\cos\theta-1)\Pi'-i\sin\theta\,\Pi'Y'_j.
```

This includes inactive inputs and crossings between remainder sectors.
It is not a native arbitrary-angle gate or permission to implement the
projector nonunitarily. For instance the legal original edge (2,3) becomes
(2,7) in physical integer labels, crossing a remainder bit. Discarding that
sector causes error at least $`|\sin\theta|`$ on its endpoint input.

## 3. A genuinely branching complete-input application

Let $`n=6k`$ and let $`\mathcal S_k`$ consist of labels whose every six-bit
block has weight at most one. It has $`7^k`$ labels and is closed under
$`p(x)=x-\mathrm{lowbit}(x)`$. Its induced computational tree has root
degree n and depth k. Keep arbitrary independent original angles on its
edges and set the others to zero. The complete target is identity outside
this set; no restriction on the actual quantum input is imposed.

The tensor encoder has $`T(E)=15k`$ and no helpers. It maps the support into
3k active bits with zero remainder. Each three-bit code is 000 through 110.
Code 111 is padding: its original preimage is the all-ones six-bit string,
not another one-hot input. The middle operation must be identity there.
The padded zero-remainder block has dimension $`8^k`$ and

```math
(8^k)^2=N.
```

Compute the complete ordered product on the support classically, then use
the [disjoint-block theorem](DISJOINT_BLOCK_COMPILER.md#1-disjoint-unitary-blocks)
once, with identity on all other labels. The emitted complete word is
$`E^\dagger\widehat M E`$, with T-count $`30k+O(N)`$ and the original
return/error contract. The remainder data bits are not newly supplied clean
workspace. Any finite overlapping program confined to this support can be
handled in the same way, independent of its length in the quantum ledger.

This support already obeyed $`(7^k)^2\le N`$ before encoding. The contribution
is the explicit inexpensive representation and its controls, not a new
asymptotic family beyond that theorem. A basis permutation followed by one
r-bit controlled block reaches at most $`2^r`$ amplitudes from a basis input.
Its error from the balanced Hopf root column is at least
$`\sqrt{1-2^r/N}`$. One such support compression therefore cannot implement
the unrestricted branching frame.

## 4. Long Pauli programs with bounded quantum algebra

There is a useful sufficient precision-sharing rule. Push Clifford gates
through a finite program of Pauli rotations and let S be the binary span of
the resulting axes. Write $`2r=\mathrm{rank}(\omega|_S)`$ for its symplectic
rank and $`c=\dim S-2r`$. Standard binary symplectic elimination supplies
an exact Clifford basis with r anticommuting generator pairs and c commuting
Z generators. The program is consequently a direct sum of $`2^r`$-dimensional
unitaries, controlled by the commuting labels, with spectators.

Multiply the entire program inside those matrices classically. Separate
scalar-on-core axes into a literal sector diagonal; the remaining sector
products are special unitary. Counting spectator copies, the squared block
budget is $`2^rN`$. For fixed r, the existing block and literal-diagonal
primitives give $`O_r(N)`$ T-count, independent of program length. All
relative sector phases are retained. This is a corollary of standard
symplectic reduction and existing native primitives, not a general compiler.
For all original Hopf levels, the Pauli span has r=n; the direct certificate
is quadratic rather than linear.

## 5. Physical common-core limitations and the price of an encoder

These statements restrict a common representation, not arbitrary native
circuits. On $`m=n+b+a`$ physical bits, allow

```math
V=C_LMC_R,\qquad
M=\bigoplus_{z\in\{0,1\}^{m-r}}M_z,\qquad M_z\in U(2^r),
```

with unrelated target-dependent physical Cliffords. The middle matrices
are unrestricted for this representation test; they are not free native
gates. It commutes with m-r independent Z Paulis. Pull them to the input,
obtaining an isotropic space S. If H is the a-dimensional initialized-code
stabilizer, its induced commuting space on the unknown n+b bits has dimension

```math
\dim(S\cap H^\perp)-\dim(S\cap H)
=\dim S-a+\dim(H\cap S^\perp)-\dim(H\cap S)
\ge n+b-r.
```

Isotropy gives the last inequality; the clean cost is not subtracted twice.
For r smaller than n, not all induced operators can live on b dirty bits.
Choose a relation with nonidentity data Pauli. Approximation error epsilon
gives an input/output Pauli intertwiner within $`2\epsilon`$. A Pauli
mapping the clean code to another syndrome gives distance $`\sqrt2`$;
differing dirty labels likewise give normalized Frobenius distance
$`\sqrt2`$, hence at least that operator distance. Define

```math
\gamma(W)=\min_{P\not\propto I,\ P,Q\ {\rm Hermitian\ Paulis}}
 \|WP-QW\|,
```

including signed and scalar Q. Therefore

```math
r\lt n\quad\Longrightarrow\quad
\|VJ_a-J_a(W\otimes I_b)\|
\ge\tfrac12\min\{\sqrt2,\gamma(W)\}.
```

No individual Clifford must preserve the initialized code. Approximate
work return suffices, and arbitrary finite helper width is allowed.
For $`E^\dagger M E`$ with a Clifford+T encoder E, each T can remove at
most one independent Pauli from the backwards-propagated commuting space:
retain its subspace commuting with the relevant Z axis. The same argument
then implies $`T(E)\ge n-r`$ whenever the target gap excludes a nontrivial
logical relation. This counts E once, not its inverse, and does not bound
an unrestricted compiler.

A legal witness is the original root star, with all its angles
$`2\arctan u`$, $`u=2^{-\ell}/(16n)`$, and
$`\ell=\lceil\sqrt N\rceil+6`$. Its distance from identity is at most
$`2^{-\ell}/8`$. For $`c=(1-u^2)/(1+u^2)`$, its root diagonal is
$`c^n`$, leaf diagonals are c, and other diagonals are one. A shifted Pauli
has a commutator entry at least $`c-c^n\ge u^2`$. A nontrivial diagonal
Pauli is detected by a marked leaf amplitude $`sc^j`$, giving at least u.
Distinct signed Paulis have separation at least $`\sqrt2`$, so proximity
to identity also gives $`\gamma(W)\ge u^2`$. For $`n\ge8`$, the resulting
$`u^2/2`$ floor exceeds $`2\eta`$. The same construction on the branching
support works because nonroot degree is at most n-5, while its root degree
is n. The common-core interface fails, although both target families are
already inexpensive complete blocks. The explicit branching encoder costs
$`5n/2`$ against the encoding-specific lower bound $`n/2`$.

## 6. Fixed logical relabeling and the regularity threshold

This extends the existing
[centered-stage theorem](STRUCTURAL_COMPILATION_LIMITS.md#arbitrary-fixed-cliffords-exactly-n-centered-diagonal-stages),
not its unrestricted scope. Let a connected spanning logical graph specify
independent two-level generators
$`H_{uv}=i(|v\rangle\langle u|-|u\rangle\langle v|)`$.
After any fixed logical basis permutation its edge displacements still span
$`\mathbb F_2^n`$, since path XORs reach every root-to-vertex difference.
Choose edges with independent displacements $`a_1,\ldots,a_n`$, binary duals
$`a_j\cdot b_i=\delta_{ij}`$, and $`z_j=b_1\oplus\cdots\oplus b_j`$.
The Hermitian Paulis $`i^{\mathrm{popcount}(a_j\&z_j)}X^{a_j}Z^{z_j}`$
anticommute pairwise and have edge-generator coefficient magnitude $`2/N`$.

A fixed physical Clifford-conjugated diagonal algebra has a commuting
Pauli basis, so it contributes at most one of these coefficients. Clean
compression preserves only I/Z clean factors; these cannot cancel logical
anticommutation. Thus an exact differentiable representation with K diagonal
slots that are Clifford-valued at the origin needs $`K\ge n`$, even with
arbitrary physical helpers and a fixed nonlinear logical relabeling. No
sufficiency claim after an arbitrary relabeling is inferred.

For a separate quantitative statement, center the factors, assume the net
fixed word acts as identity on the initialized subspace, and choose actual
Hermitian phase generators along a single-angle axis t. Put
$`q(t)=\sum_j\|H_j(t)\|`$. Scaling them together, their product has second
derivative norm at most $`q(t)^2`$. Its remainder after the linear term is
at most $`q(t)^2/2`$. For $`1\le K\lt n`$, the missing witness gives

```math
\epsilon(t)\ge\frac{2|\sin t|}{N}-\frac{q(t)^2}{2}.
```

Exact representation therefore requires $`q(t)\ge2\sqrt{|\sin t|/N}`$.
If $`\|H_j(t)\|\le B|t|^\alpha`$ on an asserted neighborhood,
$`\alpha\gt1/2`$, $`B\ge1`$, choose any positive point in that neighborhood
with

```math
t_*\le\min\{1,(NK^2B^2)^{-1/(2\alpha-1)}\}.
```

Then $`\epsilon(t_*)\ge t_*/(2N)`$. Fixed alpha, polynomial B, and an
inverse-polynomial neighborhood give an inverse-polynomial error floor.
An arbitrarily shrinking neighborhood does not. At square-root motion a
commutator of X and Z rotations generates a new order-t Y direction, so the
hypothesis cannot be silently removed. Non-Clifford cancelling baselines,
changing target-dependent skeletons, and other native words remain eligible.
The [exact center-phase audit](CENTERED_PHASE_REFLECTIONS.md) studies one
such non-Clifford baseline separately. Stage counts are not T-counts.

## 7. Evidence and attribution

[Tests](../../tests/test_paid_logical_encodings.py) check complete Boolean
permutations, promised codewords, exact Toffoli-model span loss, the full
six-bit native matrix using the retained coefficient ring, every encoded
single-Pauli control, complete encoded frames, padded branching codewords,
physical-code dimensions, and relabeled anticommuting witnesses. Bounded
matrix tests do not prove asymptotic native synthesis or all-dimension
coverage. The code emits encoders and controls, not the arbitrary-angle core.

One-hot conversion, phase-polynomial synthesis, and binary Clifford
algebra are established background: see [one-hot conversion](https://arxiv.org/abs/2206.11166),
[phase-polynomial synthesis](https://arxiv.org/abs/1601.07363), and
[binary Clifford representation](https://doi.org/10.1103/PhysRevA.68.042318).
The constructions and scope-restricted inequalities used here are given
above; these references do not supply free encoders or changed workspace
assumptions. No historical novelty certification is asserted.
