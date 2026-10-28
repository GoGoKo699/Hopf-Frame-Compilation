# A seed-degree barrier for the weak complete-probe contract

Scratch continuation after the frozen PR270 revision. No frozen repository
file is changed. This theorem strengthens the contract scope of the sector
degree barrier: the whole query endpoint need not be a wire translation.
The one complete bank probe already forces high helper dependence.

## 1. Model and robust Clifford-sector theorem

Let $`p\ge1`$ be an integer. The address is $`x\in\mathbb F_2^p`$, with
$`m=2^p`$ bank wires and q
additional arbitrary core qubits. Write $`Z_x`$ for Z on bank wire x and
identity elsewhere. The h sector helpers have a fixed common X basis r.
Assume

```math
E=\sum_r |r_X\rangle\langle r_X|\otimes E_r,
\qquad
P_x|r_X\rangle=\gamma_x(r)|(r+\Lambda x)_X\rangle,
\qquad |\gamma_x(r)|=1,
```

where $`\Lambda:\mathbb F_2^p\to\mathbb F_2^h`$ is binary-linear. The same
helper basis is fixed throughout the contract. The source core queries
$`Q_x`$ are Cliffords. Their binary symplectic label matrices $`F_x`$ have
entrywise Boolean degree at most d, with $`0\le d\le p`$.

Define the seed and its source orbit by

```math
R=E^\dagger(Z_0\otimes I_h)E,
\qquad
R_x=(Q_x\otimes P_x)R(Q_x\otimes P_x)^\dagger.
```

Assume only the complete probe comparison

```math
\max_x\|E R_xE^\dagger-Z_x\otimes I_h\|_{\rm op}
\le\varepsilon\lt\sqrt2.                         \qquad(1)
```

Every core, helper and reference input is arbitrary. Suppose every sector
$`E_r`$ is an exact Clifford up to scalar phase. Let $`M(r)`$ be its
symplectic matrix and let $`z_x`$ denote the binary label of $`Z_x`$. The seed
label is

```math
v(r)=M(r)^{-1}z_0.
```

On every active affine helper coset $`r+\mathrm{im}\Lambda`$, its Boolean
degree D satisfies

```math
D\ge p-d,\qquad \mathrm{rank}\Lambda\ge p-d.       \qquad(2)
```

If $`\mathrm{rank}\Lambda=p-d`$, the restricted degree is exactly p-d on
every active coset. There is no square-zero assumption on $`F_x-I`$, no
full-query endpoint assumption, and no assumed representation law for the
core queries. Scalar phases of the sector blocks and projective helper
queries cause no change.

## 2. Orientation, literal Pauli rounding and the top coefficient

At output helper sector r, the left-hand observable in (1) is

```math
B_{x,r}=E_r Q_x R_{r+\Lambda x}Q_x^\dagger E_r^\dagger,
\qquad R_s=E_s^\dagger Z_0E_s.                     \qquad(3)
```

The helper shift is its own inverse on binary labels. The two scalar helper
phases cancel because the seed is block diagonal. The scalar phases of E
and Q also cancel inside the displayed conjugations.

Under the Clifford-sector hypothesis, every $`B_{x,r}`$ is a signed Hermitian
Pauli. Distinct signed Hermitian Paulis have operator distance at least
$`\sqrt2`$. Equation (1) therefore forces $`B_{x,r}=Z_x`$ literally, with
positive sign. Taking labels in (3) gives

```math
F_xv(r+\Lambda x)=M(r)^{-1}z_x.                    \qquad(4)
```

This is a seed-vector identity, rather than a full symplectic intertwiner.
For fixed r, sum (4) over all addresses:

```math
\sum_xF_xv(r+\Lambda x)
 =M(r)^{-1}z_{\rm all}\ne0,
\qquad z_{\rm all}=\sum_xz_x.                     \qquad(5)
```

The sum is binary. The all-bank Z label is nonzero and $`M(r)`$ is invertible.
For any Boolean vector function of p variables, the sum of all its values is
its full-degree algebraic-normal-form coefficient. Thus the vector on the
left of (4), as a function of x, has degree p.

If v has degree D on the active helper coset, its composition with the linear
map $`x\mapsto r+\Lambda x`$ has degree at most D. Multiplying by the matrix
$`F_x`$ raises the degree by at most d. Hence $`p\le D+d`$. Since a function
on a k-dimensional coset has degree at most k, both assertions in (2) follow.

The direct observable comparison in (1) is why this rounding radius is
$`\sqrt2`$. No factor-two endpoint-to-Pauli estimate is needed. Strict
inequality is required for this rounding argument.

## 3. Affine queries and the augmented construction's leading seed coefficient

For an affine query, $`d=1`$, so the seed requires helper degree at least p-1.
Because

```math
M(r)^{-1}=J M(r)^T J,
```

with a constant binary symplectic form J, the inverse matrix has the same
maximum entry degree as the forward matrix. A product of k helper-affine
Clifford-sector stages with constant core Clifford walls therefore needs
$`k\ge p-1`$, even for the weak complete-probe contract.

For the augmented field construction, the highest address bit b spans
$`\ker\Lambda`$. Write the address as $`(a,b)`$, with p-1 lower bits, and
write

```math
F_{a,b}=F_{a,0}+bN_{\rm top},\qquad
c_C=\sum_a v(r+\Lambda a).
```

Splitting the sum in (5) over b gives the short certificate

```math
N_{\rm top}c_C=M(r)^{-1}z_{\rm all}\ne0.           \qquad(6)
```

The vector $`c_C`$ is the top coefficient on the active coset C and is
independent of its chosen base point. Consequently the pulled-back all-bank
Z label on the right of (6) is constant on that coset. In particular the
single seed pullback, not only the full tableau, has degree exactly p-1.
This is a binary-label statement; literal signs remain in (3) and (4)'s
Pauli-rounding premise.

## 4. Exact Pauli-seed variant without Clifford sector decoders

Suppose now only that every $`E_r`$ is unitary and every seed pullback
$`R_r=E_r^\dagger Z_0E_r`$ is a signed Hermitian Pauli. Keep the core queries
Clifford and require the probe equation exactly, with $`\varepsilon=0`$.

Equation (3) then forces every $`E_r^\dagger Z_xE_r`$ to be a Pauli. The map
from bank Z words to their pullback labels is linear and injective: a nonzero
bank Z word has trace zero and cannot be conjugated to a scalar Pauli. Thus
the sum of the labels of all these pullbacks is nonzero. With v now the
label of the Pauli seed, the same all-address sum proves (2).

The entire core unitary need not normalize the Pauli group. For example,
left-multiplying an exact sector decoder by a non-Clifford unitary diagonal
in the bank Z basis preserves every probe while generally destroying the
Clifford sector property. The exact seed-only theorem includes this freedom.

The uniform $`\sqrt2`$ rounding argument cannot be transferred to arbitrary
sector unitaries: their conjugated Pauli seeds need not remain Paulis. A
small non-Clifford rotation can move a Pauli continuously toward a target
without making them equal. The following distinct quantitative argument is
available instead.

## 5. A separate robust radius for arbitrary sectors with Pauli seeds

Under the Pauli-seed-only hypothesis of Section 4, fix r and define

```math
A_x=Q_xR_{r+\Lambda x}Q_x^\dagger,
\qquad B_x=E_r^\dagger Z_xE_r.
```

The $`A_x`$ are signed Hermitian Paulis, while the $`B_x`$ are commuting,
independent Hermitian involutions. The probe error says
$`\|A_x-B_x\|\le\varepsilon`$. All conclusions (2) still follow if

```math
\varepsilon\lt
\min\left\{\sqrt{2-\sqrt2},\frac2{\sqrt m}\right\}. \qquad(7)
```

First, the $`A_x`$ commute pairwise below the first threshold. If two
anticommuted, choose a common +1 eigenstate of their corresponding independent
$`B_x,B_y`$. The norm bounds imply

```math
\langle A_x\rangle,\langle A_y\rangle
 \ge1-\varepsilon^2/2.
```

But $`\|A_x+A_y\|=\sqrt2`$, so
$`2-\varepsilon^2\le\sqrt2`$, contradicting the first strict bound in (7).

If $`D+d\lt p`$, the all-address sum of the labels of the $`A_x`$ is zero.
Their commuting product is therefore $`\prod_xA_x=\sigma I`$ for
$`\sigma\in\{1,-1\}`$. Choose a joint eigenvalue pattern $`b_x`$ for the
independent $`B_x`$ such that $`\prod_xb_x=-\sigma`$, and a unit vector in
that joint eigenspace. Every joint eigenvalue pattern of the $`A_x`$ differs
from this b in at least one coordinate. Hence the exact operator inequality

```math
\sum_x(A_x-b_xI)^2\ge4I
```

holds. On the chosen joint B eigenvector its left-hand expectation is at
most $`m\varepsilon^2`$. This contradicts $`\varepsilon\lt2/\sqrt m`$.
Thus the all-address label sum is nonzero and the degree/rank argument applies.

For the standard normalization $`Q_0=I`$ up to scalar, $`A_0=B_0`$ exactly
by the seed definition. That term contributes zero to the last expectation,
so $`2/\sqrt m`$ in (7) improves to $`2/\sqrt{m-1}`$. This refinement does
not change the anticommutation threshold.

This is a width-dependent robust theorem, not the uniform Clifford-sector
radius. It uses no unjustified rounding of a non-Pauli observable.

## 6. Counting reversible nonlinear helper stages

The common-frame affine-stage bound permits fixed affine helper-coordinate
changes. Nonlinear helper motion can raise the final degree even when the
helpers are returned. It must be included in the stage ledger.

For a model with L layers of reversible helper motion of Boolean degree at
most two, and k helper-affine Clifford core stages, the final core symplectic
degree is at most $`k2^L`$. Every current helper coordinate has degree at most
$`2^L`$ in the original coordinates, and matrix multiplication adds the
degrees of the k factors. Assume the completed decoder restores the helper
sectors and satisfies the weak-probe hypotheses above. For an affine query,
k is positive and

```math
p-1\le k2^L\le2^{L+k-1},
\qquad
L+k\ge1+\left\lceil\log_2(p-1)\right\rceil
\quad(p\ge2).                                      \qquad(8)
```

This counts the specified helper-motion and core stages, including helper
uncomputation. It is not an unrestricted physical T-depth bound. Quantum
helper motion, non-Clifford intermediate sectors, or directly implemented
higher-degree stages lie outside this grammar; each still needs its own
complete native cost and work-return proof.

## 7. Bounded evidence

The scratch fixture checks a family that saturates
$`D=\mathrm{rank}\Lambda=p-d`$ for every split of the address into core-query and helper-query bits;
the actual augmented seed identity (6); complete non-Clifford Pauli-seed
examples with projective helper phases; a counterexample to inappropriate
non-Clifford Pauli rounding; the parity-distance inequality used in (7);
and a reversible helper tree that raises degree while returning its inputs.
These fixtures do not replace the all-size proof or pay a new native source.
