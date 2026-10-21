# Complete two-band sources with native read-only address

[Workspace](../../WORKSPACE.md) ·
[General two-band source screen](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) ·
[Read-only Clifford normal form](READ_ONLY_TWO_LAYER_RIGIDITY.md) ·
[Paid spectral promotion](MSS_SPECTRAL_PROMOTION.md)

The complete indicator has a stronger approximate obstruction inside
the native read-only two-band class. At amplitude $`1/\sqrt2`$ and
relative error at most $`1/8`$, such a compression source needs

```math
k\ge\left\lceil\frac{p-3}{2}\right\rceil
```

initialized clean wires. The proof permits arbitrary dirty width,
Clifford and T counts, signal mixing, clean leakage and off-code action.
It does **not** require exact flatness of the compressed block.

This closes that two-band source subclass in the sublogarithmic-clean
window of the current spectral-promotion construction. It does not close
the general moving-address approximate two-band problem, and it does
not establish an unrestricted T-depth lower bound. The distinction
between completed source return and read-only native boundaries is
essential.

## 1. Complete source and actual native scope

Let x contain p address bits, let Y contain all $`m=2^p`$ indicator
output bits, and let W be arbitrary private work. The complete target is

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle.
```

Let A contain k clean wires and let J insert $`A=0^k`$. All Y and W
inputs are arbitrary, including entanglement with an external reference.
The actual source word has two physical T/T-dagger bands,

```math
U=C_2D_2C_1D_1C_0.                                      \qquad(1)
```

Every **completed** Clifford interlayer C_j must fix the computational
address pointwise, and every physical T target must be outside x.
An elementary decomposition of C_j may temporarily move x; the
hypothesis is its completed boundary action. Neither an address that
returns only after all of U nor an algebraically reordered circuit is
enough. The identity and empty T bands are allowed.

The private-target assumption matches the current parallel-promotion
interface. It is not an additional escape from this invariant: physical
T gates directly on x, while all C_j still fix x, commute through the
word as a diagonal scalar function of x. Their phases can be absorbed
into lambda below. More generally, arbitrary diagonal address phases
do not alter the argument.

Suppose the full compressed block satisfies

```math
K=J^\dagger UJ,\qquad
\|K-\alpha I_p\|\le\varepsilon=\alpha\tau,
\qquad 0\lt\alpha\le1.                                \qquad(2)
```

The identity on W is implicit. No successful measurement or selected
input state is supplied: (2) holds on every original output and dirty
column. Arbitrary compression leakage is allowed. Literal phases are
retained below; allowing a common source phase does not weaken the
obstruction.

The address conditions make U block diagonal. Write U_x for its work
fiber, where work now includes A,Y,W, and write

```math
F_x=X_{Y_x}\otimes I_W.
```

Equation (2) implies $`\|J^\dagger U_xJ-\alpha F_x\|\le\varepsilon`$
for every x.

## 2. Three Pauli-control families, with their phases

The [address-preserving Clifford factorization](READ_ONLY_TWO_LAYER_RIGIDITY.md#1-address-preserving-clifford-interlayers)
gives

```math
C_j(x)=\lambda_j(x)P_j(x)C_j(0),\qquad
P_j(0)=I,\quad\lambda_j(0)=1.                         \qquad(3)
```

Here each P_j is a Pauli family with binary-linear labels. Its products
are projectively linear; its actual scalar multiplication cocycle is
retained. The unit-modulus lambda_j includes all address S/CZ phases
and any phases introduced by the chosen Pauli representatives.

Define the actual fixed work operators

```math
\begin{aligned}
V_*&=C_2(0)D_2,\\
S(x)&=P_1(x),\qquad A(x)=V_*S(x)V_*^\dagger,\\
W_*&=U_0,\qquad R(x)=C_0(0)^\dagger P_0(x)C_0(0),\\
P(x)&=P_2(x),\qquad\lambda(x)=\lambda_2(x)\lambda_1(x)\lambda_0(x).
\end{aligned}
```

Direct insertion into (1) gives the phase-retaining identity

```math
U_x=\lambda(x)P(x)A(x)W_*R(x).                         \qquad(4)
```

The proof below actually permits arbitrary fixed unitaries V_* and
W_* in (4). Its structural restriction is the three Pauli-control
families, not a claim that a conjugated Clifford or a compressed block
is itself a free native gate.

Let $`s_P,s_R:\mathbb F_2^p\to\mathbb F_2^k`$ be the clean X
syndromes of P and R. They are linear. Set

```math
H=\ker s_P\cap\ker s_R,\qquad
r=\dim H\ge p-2k,\qquad M=|H|=2^r.                  \qquad(5)
```

For h in H the endpoint Paulis commute with the initialized clean
projector. Their restrictions p(h),r(h) are logical Paulis with linear
labels, and

```math
J^\dagger P(h)=p(h)J^\dagger,\qquad R(h)J=Jr(h).
```

Thus (2) and (4) imply

```math
\begin{aligned}
B_h&=J^\dagger A(h)W_*J,\\
L_h&=\lambda(h)^{-1}p(h)^\dagger F_h r(h)^\dagger,\\
\|B_h-\alpha L_h\|&\le\varepsilon.                  \qquad(6)
\end{aligned}
```

Each L_h is a unit-modulus scalar times a logical Pauli. If e_h denotes
the binary label of the original output Pauli $`X_{Y_h}`$, its label is

```math
\mathrm{lab}(L_h)=e_h+\ell_0(h)                     \qquad(7)
```

for a linear map ell_0. The vectors e_h are distinct independent output
basis vectors, not merely distinct values of a one-bit predicate.
This is where the **complete bank** enters the argument.

## 3. A matrix coefficient has Fourier nuclear norm at most one

The projective cocycle of A must not be silently dropped. Introduce a
mathematical auxiliary copy of the work Pauli representation and put

```math
\Gamma(h)=A(h)\otimes\overline{S(h)}.
```

If $`S(h)S(g)=c(h,g)S(h+g)`$, the conjugate representation contributes
$`\overline{c(h,g)}`$. Consequently Gamma is an **ordinary** unitary
representation of the additive group H. Its Fourier operators

```math
\Pi_\chi=\frac1M\sum_{h\in H}\chi(h)\Gamma(h)
```

are mutually orthogonal projectors summing to identity. All characters
chi of H take values in $`\{1,-1\}`$.

Let

```math
E=J\otimes I,\qquad F=W_*J\otimes I.
```

Both are isometries with a common input dimension d. The Fourier
coefficients of $`B_h\otimes\overline{S(h)}`$ are exactly
$`E^\dagger\Pi_\chi F`$. Write
$`\|Z\|_{s,d}=(d^{-1}\mathrm{Tr}|Z|^s)^{1/s}`$ for normalized
Schatten norms. The Hilbert--Schmidt product inequality and
Cauchy--Schwarz give

```math
\begin{aligned}
\sum_\chi\|E^\dagger\Pi_\chi F\|_{1,d}
&\le\frac1d\sum_\chi
       \|\Pi_\chi E\|_{\rm HS}\|\Pi_\chi F\|_{\rm HS}\\
&\le\frac1d
 \sqrt{\mathrm{Tr}E^\dagger E\,
       \mathrm{Tr}F^\dagger F}=1.                    \qquad(8)
\end{aligned}
```

The tensor copy, Fourier projectors and Schatten calculation are proof
invariants only. No auxiliary copy, character transform, projector or
sum of operators is asserted to be an available native operation.

## 4. Complete output labels force a large Fourier nuclear norm

Put

```math
T_h=L_h\otimes\overline{S(h)},\qquad
D_\chi=\frac1M\sum_{h\in H}\chi(h)T_h.
```

Every T_h is a scalar times a Pauli unitary, and its binary label has
the form $`e_h+\ell(h)`$ for a linear ell into the enlarged Pauli
space. Parseval gives

```math
\sum_\chi\|D_\chi\|_{2,d}^2=1.                       \qquad(9)
```

In the aggregate fourth moment, character orthogonality leaves only
ordered tuples $`(a,b,c,d)`$ with $`a+b+c+d=0`$. Their linear label
terms therefore cancel. A surviving Pauli trace must also satisfy

```math
e_a+e_b+e_c+e_d=0.
```

Independence of the original output labels means each index occurs an
even number of times. There are exactly $`3M^2-2M`$ such ordered
tuples. Their actual phases can only lower the upper bound obtained
by taking absolute values. Therefore

```math
\sum_\chi\|D_\chi\|_{4,d}^4
\le\frac{3M^2-2M}{M^3}
=\frac{3M-2}{M^2}\le\frac3M.                        \qquad(10)
```

Apply Schatten interpolation to the direct sum of all D_chi, with
trace measure $`d^{-1}\sum_\chi\mathrm{Tr}`$. The normalization
is not divided by another factor M. Equations (9)--(10) yield

```math
\sum_\chi\|D_\chi\|_{1,d}
\ge\frac{M}{\sqrt{3M-2}}
\ge\sqrt{\frac M3}.                                  \qquad(11)
```

This argument retains arbitrary phases in (6). In particular, it does
not divide a matrix coefficient by an arbitrary scalar function and
then assume its Fourier nuclear norm is unchanged.

## 5. Robust complete-source bound

Let the error family be

```math
Z_h=(B_h-\alpha L_h)\otimes\overline{S(h)}.
```

Every Z_h has operator norm at most epsilon. Its Fourier coefficients
obey, by normalized Hilbert--Schmidt Parseval,

```math
\sum_\chi\|\widehat Z_\chi\|_{1,d}
\le\sqrt M\left(\sum_\chi
                  \|\widehat Z_\chi\|_{2,d}^2\right)^{1/2}
\le\sqrt M\,\varepsilon.                             \qquad(12)
```

Combining (8), (11) and (12) proves the necessary inequalities

```math
\boxed{\quad
1\ge\frac{\alpha M}{\sqrt{3M-2}}-\varepsilon\sqrt M
\ge\sqrt M\left(\frac\alpha{\sqrt3}-\varepsilon\right).
\quad}                                               \qquad(13)
```

In particular, if $`\tau\lt1/\sqrt3`$, then

```math
2^{p-2k}\le M
\le\frac1{\alpha^2(1/\sqrt3-\tau)^2}.                \qquad(14)
```

At $`\alpha=1/\sqrt2`$ and $`\tau\le1/8`$, the fact
$`1/\sqrt3>1/2`$ gives $`M\lt128/9\lt16`$. Since M is a power of two,
r is at most three, proving

```math
k\ge\left\lceil\frac{p-3}{2}\right\rceil.            \qquad(15)
```

This is an operator-norm approximation theorem, not an exact logical
action theorem. It applies whether or not $`K^\dagger K=I/2`$.
Exact flatness would improve subsequent recovery and return; it cannot
evade (13). The same proof tolerates an optimized common source phase,
which only changes the scalar coefficients of the target Paulis.

## 6. Resource frontier and the next unpaid operation

The current [spectral-promotion schedule](MSS_SPECTRAL_PROMOTION.md#1-source-contract-and-synchronized-native-stages)
requires precisely the native boundary property in Section 1. Its
amplitude recovery preserves that property when it holds for U.
Consequently a two-band source offered directly to that schedule, at
the constant accuracy of (15) or better, cannot have
$`k=o(\log(p+2))`$. This conclusion is stronger than the general
rank--norm screen **within this narrower native class**.

For the existing emitted promotion word, recovery uses k+1 clean bits
before $`h=2^{k+1}`$ private source copies are made arbitrary. Each copy
has its own q-bit old output bank, with $`q=2^p`$ for this source. Thus
(15) makes those banks alone have at least

```math
hq\ge2^{3p/2-O(1)}
```

wires. They do not fit the required $`2^p\mathrm{poly}(p)`$
source envelope; the explicit selective truth-table banks are additional.
The retained root fusion with source size $`p-O_\tau(\log(p+2))`$
does not absorb this exponential multiplicative factor. This describes
the actual supplied construction, not a lower bound on every possible
promotion or joint resynthesis. An upper bound on recovery depth is
not reinterpreted as a necessary depth cost.

The first unsupported inference in applying that wrapper to an arbitrary
moving-address source would be:

> The completed source returns its address, so all of its native
> Clifford boundaries may be treated as read-only and its private
> copies may be run in parallel on one shared address.

Neither implication follows. A remaining approximate two-band source
can move the address between its completed Clifford boundaries, but
then its compatible native sharing schedule is unpaid. Merely placing
T gates on an otherwise fixed address contributes only the scalar
phases already allowed above; that alone does not evade the invariant.
Alternatively, additional native bands or a different complete-bank
mechanism can leave the three-Pauli-control grammar.

The single-predicate zero-reflection construction is not contradicted:
it does not supply the independent output labels e_h used in (10), and
its moving native address already lies outside this class. Likewise,
none of the tensor/Fourier witnesses above supplies a native readout.
The general approximate two-band source question, the focused original
physical-pair scout, and the unrestricted complete-frame depth frontier
remain distinct open tasks.

## 7. Bounded exact checks

[`test_read_only_two_band_source_gap.py`](../../tests/test_read_only_two_band_source_gap.py)
checks the nontrivial projective cocycle cancellation with literal Pauli
phases. It also checks Parseval and the aggregate fourth-moment identity
on two- and three-bit address subgroups, including eighth-root target
phases and affine Pauli-label perturbations. The checks use exact
integer Pauli labels and cyclotomic coefficients; no dense bank matrix
or numerical source search is used. They audit the finite identities,
not the uniform proof or a native source construction.
