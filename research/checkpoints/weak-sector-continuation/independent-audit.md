# Independent audit of the weak-probe sector-degree barrier

**Verdict: valid under the stated sector and Pauli-seed hypotheses.** The
constant-radius Clifford-sector statement, the exact Pauli-seed statement,
and the stronger width-dependent robust Pauli-seed statement are distinct.
This audit concerns scratch continuation only; it changes no frozen PR270 file.

## 1. Precise sector equation and phase cancellation

The physical target bank has $`m=2^p`$ wires, disjoint from h designated
helper qubits. Decompose in the helper X basis, indexed by
$`r\in\mathbb F_2^h`$, and assume

```math
E=\bigoplus_r E_r.
```

Each sector acts on a core that contains the entire bank and may contain
further spectator qubits. The source is $`Q_x\otimes P_x`$, where the
core query is Clifford with symplectic matrix $`F_x`$ of entrywise Boolean
degree at most d. In the helper X basis, the Pauli query shifts r by
$`\Lambda x`$. It may also have an r-dependent scalar phase. Such phases
cancel in conjugating a block-diagonal probe; literal commutation of the
helper Paulis is unnecessary. Sector-dependent scalar phases of $`E_r`$
also cancel from its conjugations.

The pulled-back seed has block $`R_r=E_r^\dagger Z_0E_r`$. Thus the
physical probe comparison is equivalent to, for every r and x,

```math
\left\|E_rQ_xR_{r+\Lambda x}Q_x^\dagger E_r^\dagger-Z_x\right\|
\le\epsilon.
```

The shift has this orientation because its inverse equals itself over the
binary helper labels. The core source action is forward conjugation by
$`Q_x`$, not its adjoint. Full operator norm is the maximum of these
sector norms; restricting to one initialized helper sector would not suffice.

If every $`E_r`$ is Clifford up to scalar, the first operator in this norm
is a signed Hermitian Pauli. Distinct signed Hermitian Paulis are separated
by at least $`\sqrt2`$: the distance is $`\sqrt2`$ for anticommuting
labels and two for distinct commuting labels or opposite signs. Therefore
$`\epsilon\lt\sqrt2`$ forces exact equality with the positive $`Z_x`$.
Writing $`M(r)`$ for the core symplectic matrix of $`E_r`$ and
$`z_x`$ for the bank Z label gives exactly

```math
v(r)=M(r)^{-1}z_0,\qquad
F_xv(r+\Lambda x)=M(r)^{-1}z_x.
```

No scalar-phase assumption has been dropped in reaching this binary equation.

## 2. The address sum forces high degree on every active coset

Let $`S=\mathrm{im}\Lambda`$ and $`k=\dim S`$. On each affine coset
$`C=r+S`$, define $`D_C`$ as the maximum Boolean degree of the coordinates
of v after any affine bijection $`\mathbb F_2^k\to C`$. This definition
is independent of the chosen affine coordinates. In particular,
$`D_C\le k`$, and the address function $`x\mapsto v(r+\Lambda x)`$
has degree at most $`D_C`$.

Sum the preceding sector equation over all addresses, in the binary label
space. Its right side is

```math
M(r)^{-1}\sum_x z_x\ne0.
```

The sum of the m distinct canonical bank Z labels is the nonzero all-bank-Z
label, and the symplectic map is invertible. The left side is the full
Boolean-cube sum of a vector-valued polynomial of degree at most
$`d+D_C`$. Such a sum vanishes when that degree is below p. Hence every
active coset necessarily satisfies

```math
D_C\ge p-d,\qquad \mathrm{rank}\Lambda\ge p-d.
```

The assertion is informative when $`d\lt p`$. It does not require an
additive law for the core queries, any field-specific nilpotence identity,
or an assumption that their phases have low degree.

## 3. Exact arbitrary-sector-unitary version

For exact probes the Clifford-sector hypothesis can be weakened: require
only that every $`R_r=E_r^\dagger Z_0E_r`$ is a signed Hermitian Pauli.
Define v by these seed labels. Exact transport implies that every
$`E_r^\dagger Z_xE_r=Q_xR_{r+\Lambda x}Q_x^\dagger`$ is also Pauli.
Conjugation therefore gives a linear injective label map $`\phi_r`$
from the bank-Z subgroup, and

```math
F_xv(r+\Lambda x)=\phi_r(z_x),\qquad
\sum_x\phi_r(z_x)=\phi_r\!\left(\sum_xz_x\right)\ne0.
```

Injectivity follows because a nonzero positive Z word cannot be conjugated
to a scalar. The same polynomial argument proves the same degree and rank
bounds. The helper sectors may otherwise be arbitrary unitaries. Approximate
outputs in this class are not necessarily Pauli, so the preceding
$`\sqrt2`$ discretization cannot be applied to them.

## 4. Stronger robust statement with Pauli seeds only

The parent's additional robust argument also passes independently. Assume
the Pauli-seed hypothesis, but allow arbitrary sector unitaries. For fixed r
put

```math
A_x=Q_xR_{r+\Lambda x}Q_x^\dagger,\qquad
B_x=E_r^\dagger Z_xE_r,\qquad \|A_x-B_x\|\le\epsilon.
```

The $`A_x`$ are signed Hermitian Paulis. The $`B_x`$ are commuting
involutions with every joint sign pattern present, because they are the
simultaneous conjugates of the m independent bank-Z generators.

If $`A_i,A_j`$ anticommute, choose a common $`+1,+1`$ eigenstate of
$`B_i,B_j`$. On that state,

```math
2-2\langle A_i\rangle=\|(A_i-I)\psi\|^2\le\epsilon^2,
\qquad
2-\epsilon^2\le\langle A_i+A_j\rangle\le\sqrt2.
```

Thus anticommutation requires
$`\epsilon\ge\sqrt{2-\sqrt2}`$. Below that threshold the $`A_x`$
all commute. If $`d+D_C\lt p`$, their label sum vanishes, so
$`\prod_xA_x=\sigma I`$ with $`\sigma\in\{1,-1\}`$.
Choose a joint target sign pattern b with $`\prod_xb_x=-\sigma`$;
its joint eigenspace is nonzero. Every joint A eigenpattern differs from b
in at least one coordinate, so the operator inequality

```math
\sum_x(A_x-b_xI)^2\ge4I
```

holds. On a unit vector from the chosen B eigenspace, the left expectation
is at most $`m\epsilon^2`$. Therefore a degree defect requires
$`\epsilon\ge2/\sqrt m`$. The degree and helper-rank barrier consequently
holds under the separately scoped robust radius

```math
\epsilon\lt\min\!\left\{\sqrt{2-\sqrt2},\frac2{\sqrt m}\right\}.
```

Independence of the target bank-Z generators is essential to the sign-pattern
step. The full operator-norm comparison is essential to evaluating it on the
selected joint eigenstates. This proof does not claim approximate Pauli-domain
invariance and does not snap the approximate outputs to exact Paulis.

## 5. Affine-stage consequence and the augmented construction

If the sector Clifford is a product of s stages whose core symplectic
matrices are affine in the same helper-X coordinates, with arbitrary fixed
core Clifford walls, its full symplectic matrix has degree at most s.
Symplectic inversion preserves entrywise degree because
$`M^{-1}=JM^TJ`$. Thus $`D_C\le s`$ and necessarily

```math
s\ge p-d.
```

More generally the sum of the stage degrees must be at least $`p-d`$.
This is a restriction on a specified common-sector factorization. It is
**not** an unrestricted native T-depth lower bound. A native implementation
may mix helper X sectors internally, use stages of higher helper degree,
or have non-Pauli seed pullbacks. Completed affine conditional-Clifford
batches cannot be identified with physical T bands without a separate proof.

The retained augmented field construction satisfies the hypotheses exactly.
Its field symplectic query is affine, so $`d=1`$. Its helper shift depends
on the lower $`p-1`$ address bits and not the top bit. The largest helper's
injective field embedding already has rank $`p-1`$; smaller shared helpers
add no address rank. Hence

```math
\mathrm{rank}\Lambda=p-1,\qquad D_C=p-1
```

on every active coset. Each of the two completed quadratic-phase batches
per recursion level has an affine helper-dependent core symplectic matrix.
The present construction therefore uses $`2(p-1)`$ such stages; the bound
requires at least $`p-1`$ in this model. It rules out a bounded-stage global
batching in that model, while leaving other gauges and native factorizations
open.

A useful exact specialization writes the address as $`(a,b)`$, with the top
bit b in $`\ker\Lambda`$, and $`F_{a,b}=F_{a,0}+bN_{\rm top}`$.
If $`c_C=\sum_a v(r+\Lambda a)`$ is the top active-coset ANF coefficient,
the summed Clifford-sector equation gives

```math
N_{\rm top}c_C=M(r)^{-1}z_{\rm all}\ne0.
```

Thus the highest coefficient is nonzero, and the pulled-back all-bank-Z label
is constant throughout an active coset. This provides an additional exact
fixture identity without strengthening the scope of the depth claim.
