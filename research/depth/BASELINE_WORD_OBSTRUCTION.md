# Why a word cannot cancel the two native baseline families alone

[Five-layer reduction](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Current depth frontier](../../docs/OPEN_PROBLEM.md)

The Pauli-comparison proof cancels a quadratic label baseline on an
eight-vertex address cube. For a general Clifford comparison, replacing
that label sum by an operator product does not give a corresponding
universal cancellation rule. This note proves a stronger obstruction:
**no word in the comparison values and their actual inverses can both
annihilate every two-family native baseline and detect a central target
factor.** The word may have any finite length and any prescribed ordering.

This is a limitation of a proposed proof mechanism. It is not a depth
lower bound, a returned-query construction, or an obstruction to using
the additional five-layer cocycle. In particular, the two baseline
families below are not claimed to have all the shared-source and return
compatibilities of an actual five-layer query.

## 1. The exact question and theorem

Let the address group be $`V=\mathbb F_2^a`$, with $`a\ge1`$. A native
one-T Pauli family has the form

```math
A_x=C D P_xD^\dagger C^\dagger,
```

where C is a fixed Clifford, D is one physical T/T-dagger band, and
the Pauli labels of P are linear in x. The word "one-T" describes the
conjugating band; an implementation charges both D and its inverse.
Consider the baseline class

```math
K_x^0=A_x^\dagger N_x,\qquad K_0^0=I,
```

with A and N such native families. This is already a subclass of the
factorization $`A_x^\dagger P_xN_x`$ in the five-layer reduction, obtained
by setting P to identity. Here A and N may be chosen independently.

Fix a finite word

```math
\mathcal W(K)=K_{z_1}^{\epsilon_1}\cdots K_{z_s}^{\epsilon_s},
\qquad z_j\in V,\quad \epsilon_j\in\{1,-1\},
```

where exponent minus one means the actual adjoint. Set

```math
n_z=\sum_{j:z_j=z}\epsilon_j.
```

**Baseline-word obstruction.** Suppose $`\mathcal W(K^0)`$ is a scalar
multiple of identity for every finite width and every pair of native
families A,N above. Then

```math
n_z=0\quad\text{for every }z\ne0.                         \qquad(1)
```

Consequently let Q be any fixed involution commuting with every baseline
value, and let $`g:V\to\mathbb F_2`$ satisfy $`g(0)=0`$. For

```math
\widehat K_x=K_x^0Q^{g(x)},
```

the same word has exactly the same value:

```math
\mathcal W(\widehat K)=\mathcal W(K^0).                   \qquad(2)
```

A nonidentity Pauli on one disjoint arbitrary qubit is a permitted Q.
Thus such a universal word cannot detect a Boolean derivative, or any
other variation of g, through that central target factor. Equation (1)
is necessary, not sufficient, for universal baseline cancellation.

The same conclusion holds under the weaker, approximate cancellation
hypothesis

```math
\inf_{|\omega|=1}\|\mathcal W(K^0)-\omega I\|_{\rm op}
\lt\frac1{\sqrt2}
\quad\text{for every such native baseline}.              \qquad(3)
```

This is a width-independent obstruction to a proposed cancellation rule,
not a robustness theorem for five-layer returned queries.

The word contains only comparison values and their inverses. It does
not include separately known native factors, probes, the middle family
B, or a cocycle-dependent selection rule. Those are possible sources
of additional information and are outside the obstruction.

## 2. Native products of unbounded odd symplectic order

For each $`m\ge1`$, use m arbitrary work qubits and put

```math
D=T^{\otimes m},\qquad X_* = X_1\cdots X_m,
\qquad A=DX_*D^\dagger.
```

Choose the explicit Clifford

```math
C=(\mathrm{CX}_{1\to2}\mathrm{CX}_{2\to3}
 \cdots\mathrm{CX}_{m-1\to m})H^{\otimes m},
\qquad N=CDX_*D^\dagger C^\dagger.                       \qquad(4)
```

Products act right to left. Thus the CNOTs in C are executed in descending
control order, after the H gates. Its individual axis images are

```math
CZ_jC^\dagger=X_jX_{j+1}\quad(j\lt m),\qquad
CZ_mC^\dagger=X_m.
```

A and N are literal Hermitian involutions, each a native one-T image of
the Pauli $`X_*`$. Their symplectic actions are the products of
transvections along, respectively, the odd and even members of the chain

```math
c_{2j-1}=Z_j,\qquad
c_{2j}=X_jX_{j+1}\ (j\lt m),\qquad c_{2m}=X_m.              \qquad(5)
```

Here labels are used in the binary symplectic space. Each parity class
of this list is independent and isotropic. The only nonzero pairings
between different chain members are those of adjacent members.

To calculate the order, take the even-weight subspace of
$`\mathbb F_2^{2m+1}`$ with the restricted dot product. It is a
nondegenerate alternating space: its orthogonal complement in the full
space is spanned by the all-ones vector, which has odd weight. The
vectors $`e_i+e_{i+1}`$ form a basis with exactly the pairings in (5).
The transvection along $`e_i+e_{i+1}`$ swaps coordinates i and i+1.
Under this symplectic identification, the action of AN is therefore the
product of the two disjoint matchings

```math
(1\ 2)(3\ 4)\cdots(2m-1\ 2m),
\qquad
(2\ 3)(4\ 5)\cdots(2m\ 2m+1).
```

Their product is the cycle

```math
(1\ 2\ 4\ \cdots\ 2m\ 2m+1\ 2m-1\ \cdots\ 3).
```

The permutation action on the even-weight subspace is faithful for
$`2m+1\ge3`$: fixing every vector $`e_i+e_j`$ fixes every unordered pair,
and hence every coordinate. Consequently

```math
\mathrm{ord}(S_{AN})=2m+1.                         \qquad(6)
```

In particular, $`(AN)^k`$ cannot be a scalar when $`2m+1`$ does not divide
k. We need only this necessary condition; no claim about the literal
unitary order is required.

## 3. Proof of the obstruction

For any binary linear form $`\ell:V\to\mathbb F_2`$, take

```math
A_x=A^{\ell(x)},\qquad N_x=N^{\ell(x)}.
```

These are exact representations, not just projective ones, and are
native one-T Pauli families because $`P_x=X_*^{\ell(x)}`$. Since A is
Hermitian,

```math
K_x^0=(AN)^{\ell(x)}.
```

All values of this particular baseline commute. The proposed word is
therefore

```math
\mathcal W(K^0)=(AN)^{k_\ell},\qquad
k_\ell=\sum_{z\ne0}n_z\ell(z)\in\mathbb Z.
```

If this is a scalar for every m, (6) forces every positive odd integer
at least three to divide $`k_\ell`$. Hence $`k_\ell=0`$ for every linear
form. This is an equality of integers, not just a parity constraint.

The approximate hypothesis (3) gives the same conclusion. If
$`2m+1`$ does not divide $`k_\ell`$, the Clifford
$`U=(AN)^{k_\ell}`$ has nonidentity symplectic action. Some Hermitian
Pauli P is sent to a Pauli of a different binary label, so

```math
\sqrt2\le\|UPU^\dagger-P\|
\le2\|U-\omega I\|\quad\text{for every }|\omega|=1.
```

This contradicts (3). Thus the necessary divisibilities and
$`k_\ell=0`$ also follow in that case. No approximate native source or
rounding of an arbitrary unitary is used in this argument.

Index the nonzero linear forms and nonzero vectors by rows and columns
of the integer matrix $`M_{\ell,z}=\ell(z)`$. For $`a=1`$ it is the
one-by-one identity. For $`a\ge2`$, distinct nonzero linear forms are
independent, giving

```math
MM^{\mathsf T}=2^{a-2}(I+\mathbf1\mathbf1^{\mathsf T}).
```

This matrix is positive definite over the rationals, so M is invertible.
The equations $`Mn=0`$ prove (1). Finally Q commutes with every factor,
so moving all target factors to the right yields

```math
\mathcal W(\widehat K)
=\mathcal W(K^0)Q^{\sum_z n_zg(z)}
=\mathcal W(K^0),
```

including the actual inverse factors. The last equality uses (1) and
$`g(0)=0`$. This proves (2), with literal equality of the two words.

## 4. Cost, error and the return boundary

Let the nonzero form $`\ell`$ involve s address bits. Implement each
$`X_*^{\ell(x)}`$ by a CNOT from each of those s address bits to each
of the m work qubits. This uses ms elementary Cliffords, preserves each
address pointwise, and requires no parity scratch. The displayed circuit
for $`K_x^0=A_xN_x`$ has

```math
T=4m,\qquad G=4m-2+2ms,\qquad D_T=4.
```

There are m simultaneously arbitrary work qubits, zero clean qubits,
and zero further helpers. The fixed source C and its actual inverse
each cost $`2m-1`$ elementary Cliffords. Every T in D and every T-dagger
in its inverse is included. Within every T band the targets are distinct.
The actual inverse word has the same ledger; a source followed by that
inverse costs twice these displayed counts and restores all work exactly.
These are conservative costs for the displayed word, not optimality claims.

The identities are exact on every address and work column, with literal
phases and arbitrary external references. There is no precision parameter,
approximation error or leakage. **The source K itself transforms the work;
it is not a returned XOR or scalar query.** Appending Q on a disjoint qubit
in (2) is an operator comparison, not a claimed free implementation of an
arbitrary function g. No source, control or inverse is being supplied as
an uncharged compiler subroutine.

## 5. Consequence for the next invariant

The eight-factor Pauli-label cube cannot be extended merely by finding a
different universal ordering of the non-Pauli comparison values. Even
allowing an arbitrarily long fixed word does not solve that baseline-only
problem. Repeated multiplicative differences illustrate the issue: for
a direction with $`\ell(h)=1`$, successive differences of
$`(AN)^{\ell(x)}`$ give powers $`(AN)^{\pm2^{r-1}}`$, which remain
non-scalar for every r already when $`m=1`$.

The useful resumption question must retain more information: can the
actual five-layer relation between B, J and K, together with their shared
native sources, force a target-detecting relation or scalar invariant?
An invariant using probes or separate baseline factors is also not ruled
out here, but its conjugation transport and every error use must be proved.
This result neither settles the general five-layer degree question nor
changes the complete-frame upper law, selected Result D, or the scope of
the retained four-layer obstruction.

The bounded [evidence](../../tests/test_baseline_word_obstruction.py)
checks the physical source and its actual inverse on all small columns,
the axis images, exact binary symplectic orders, the incidence Gram
matrix, counterexamples to cube words, and the literal central-factor
cancellation. The arbitrary-width and all-word conclusions follow from
the analytic argument above, not from those finite checks.
