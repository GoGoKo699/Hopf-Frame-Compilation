# A four-band source for every balanced Clifford involution

[Clifford-interleaved lookup](CLIFFORD_INTERLEAVED_LOOKUP.md) ·
[Complete derivative contract](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) ·
[Native phase conventions](NATIVE_BAND_SOURCE_STRUCTURE.md)

Every traceless Hermitian Clifford on n qubits is the exact image of a Pauli
under a source with at most four native T bands, at most $`3n-2`$ T gates,
$`O(n^2)`$ elementary Clifford gates and no additional qubits. All n inputs
are arbitrary. Literal signs and the actual inverse are retained. This is
a useful source component, not a complete lookup construction or a stronger
frame-depth bound.

It also removes the first nonlinear block as a separate search variable in
a three-block source. After peeling two weighted cubic conjugations from the
complete target, the remaining question is simply whether the result is an
actual Clifford. The theorem constructs its first source when that test
succeeds. The two outer conjugations, their Clifford coordinates and their
complete action still have to be found and paid.

<a id="balanced-clifford-source"></a>

## 1. A traceless Clifford involution has a translation coordinate

Let K be an n-qubit Clifford with

```math
K=K^\dagger,\qquad K^2=I,\qquad \mathrm{Tr}(K)=0.
```

Its binary symplectic matrix S obeys $`S^2=I`$. Put $`N=S-I`$ over
$`\mathbb F_2`$. The image of N is isotropic, since symplecticity and
$`S^{-1}=S`$ give $`\langle Nu,Nv\rangle=0`$. Also

```math
\ker N=(\mathrm{im}\,N)^\perp.
```

For every label v in $`\ker N`$, let $`\epsilon(v)\in\{+1,-1\}`$
be the literal sign in $`K\tau(v)K^\dagger=\epsilon(v)\tau(v)`$.
Here $`\tau(v)`$ is the Hermitian Pauli with that label. These signs form
a character even if the fixed Paulis do not commute: multiplying two
conjugation identities preserves their common Pauli multiplication phase.
The trace of conjugation on the full operator space gives

```math
0=|\mathrm{Tr}(K)|^2
  =\mathrm{Tr}(\mathrm{Ad}_K)
  =\sum_{v\in\ker N}\epsilon(v).
```

Thus some fixed-label Pauli R has sign minus one. The space spanned by
$`\mathrm{im}\,N`$ and R is isotropic. Extend it to a Lagrangian L;
because it contains $`\mathrm{im}\,N`$, it lies in $`\ker N`$.
A Clifford coordinate change maps L to the computational Z algebra.
In those coordinates K sends every computational Z to itself up to sign,
so its basis permutation is a fixed translation $`z\mapsto z+a`$.
The minus sign on R forces $`a\ne0`$. An invertible linear change puts
$`a=e_0`$, the first coordinate vector.

This construction uses symplectic linear algebra on the full physical
register. It assumes no initialized coordinate. The required Clifford
coordinate changes have $`O(n^2)`$ H, S and CNOT gates: symplectic
elimination first chooses the commuting Z basis, and ordinary binary
elimination fixes its translation. Their actual native inverses have the
same asymptotic bound, with S-dagger expanded as three S gates if required.

## 2. Integrating the literal Clifford phase

Write the translated coordinates as $`(t,z)`$, with t one bit and z
containing the other $`n-1`$ bits. A monomial Clifford has one global scalar
times fourth-root phases from a quadratic polynomial modulo four. Its
Hermiticity makes that global scalar an eighth root: squaring the operator
at any basis column fixes its square as a fourth root. Accordingly,
for some integer k modulo eight and Clifford phase polynomial q modulo four,

```math
K'|t,z\rangle
 =\omega^{(1-2t)[k+2q(z)]}|t+1,z\rangle,
\qquad \omega=e^{i\pi/4}.                              \qquad(1)
```

The constant term of q is absorbed into k. The phase at t=1 is the literal
complex conjugate of the phase at t=0, because $`K'^2=I`$. Define

```math
D|t,z\rangle=\omega^{t[k+2q(z)]}|t,z\rangle.
```

Direct action on every column gives

```math
D X_t D^\dagger=K'.                                    \qquad(2)
```

The exponent is a weighted cubic phase: q has arbitrary linear coefficients
modulo four and even quadratic coefficients. Multiplication by $`2t`$
therefore gives twice-linear quadratic terms and four-times-cubic terms,
with the remaining linear term $`kt`$. This is a diagonal third-level
source. It need not be a product of cubic sign gates alone; an odd k is
already needed for an appropriate one-qubit Clifford involution.

If C denotes the accumulated coordinate Clifford, $`K=CK'C^\dagger`$.
For any specified nonidentity signed Hermitian Pauli P, choose a Clifford B
with $`BPB^\dagger=X_t`$. Then

```math
U=CDB,\qquad UPU^\dagger=K.                            \qquad(3)
```

The Clifford B absorbs the sign of P. Any global scalar in C or B cancels
in this conjugation, while k retains the literal sign of K. The emitted
word and its actual inverse are used; no projective replacement of K occurs.

## 3. A sparse quadratic basis gives four physical bands

The apparent dense quadratic phase in (1) does not require a quadratic
T count. Define the symmetric binary bilinear form b by

```math
2b(u,v)=q(u+v)-q(u)-q(v)\pmod4.
```

A nonzero diagonal entry splits off an orthogonal one-dimensional block.
Once the remaining form is alternating, a nonzero off-diagonal entry splits
off a hyperbolic pair. The remaining radical is orthogonal to every block.
This elementary elimination gives an invertible linear change of z such
that

```math
q(z)=\sum_{i=1}^{s}\varepsilon_i a_i
       +2\sum_{j=1}^{h}u_jv_j+2\sum_{\ell\in E}z_\ell\pmod4,
\quad \varepsilon_i\in\{+1,-1\},\quad s+2h\le n-1.       \qquad(4)
```

All singletons $`a_i`$ and pair coordinates $`u_j,v_j`$ are distinct.
The even-linear set E uses the remaining coordinates, including pair
coordinates when needed. Indeed the value of q on each basis vector fixes
its linear coefficient; its polar form fixes all cross coefficients. Odd
singleton values are plus or minus one modulo four. Every other basis
value is zero or two. This proves (4), including degenerate forms.
The coordinate change costs $`O(n^2)`$ CNOT gates and is included in C and
D's physical basis walls.

In these coordinates D is the product of root powers, controlled S or
S-dagger on the disjoint singletons, CCZ on the disjoint pairs, and CZ
for the even-linear terms. Use the exact integer phase identities

```math
\begin{aligned}
2ta&=t+a-(t\mathbin\oplus a),\\
4tuv&=t+u+v-(t\mathbin\oplus u)-(t\mathbin\oplus v)
       -(u\mathbin\oplus v)+(t\mathbin\oplus u\mathbin\oplus v)
       \pmod8.
\end{aligned}                                                   \qquad(5)
```

Let $`r=(k+\sum_i\varepsilon_i+h)\bmod8`$, chosen in $`\{0,\ldots,7\}`$.
Emit $`S_t^{\lfloor r/2\rfloor}`$ and, when r is odd, one $`T_t`$.
The other T gates fit in these four genuine physical bands:

| Band | Native phases on distinct physical targets |
| --- | --- |
| 1 | $`T_{a_i}^{\varepsilon_i}`$, $`T_{u_j},T_{v_j}`$, and the optional $`T_t`$ |
| 2 | $`T^{ -\varepsilon_i}`$ on $`t\oplus a_i`$; T-dagger on $`t\oplus u_j,t\oplus v_j`$ |
| 3 | T-dagger on $`u_j\oplus v_j`$ |
| 4 | T on $`t\oplus u_j\oplus v_j`$ |

Each parity is computed into its own listed non-root wire and uncomputed
immediately after that band. Bands 2, 3 and 4 require respectively
$`2s+4h`$, $`2h`$ and $`4h`$ CNOT gates including uncomputation.
Shared root controls impose no additional T band. The even-linear terms
are $`\mathrm{CZ}(t,z_\ell)`$, each paid as H, CNOT, H. No temporary
coordinate is clean or borrowed; these are invertible parity changes on
arbitrary inputs.

The exact core ledger is therefore

```math
\begin{aligned}
T_D&=2s+6h+(r\bmod2)\le3n-2,& D_T(D)&\le4,\\
G_{{\rm Cliff},D}&=2s+10h+3|E|+\lfloor r/2\rfloor,&
a_{\rm clean}&=0,\qquad b_{\rm extra}=0.
\end{aligned}                                                   \qquad(6)
```

Empty bands may be omitted. The actual inverse reverses all words and
T orientations; it has identical T count and T-depth. Replacing reversed
S gates by positive S powers can increase its Clifford count only by a
constant factor. Combining the source core with the coordinate changes
and B proves $`T_U\le3n-2`$, $`D_T(U)\le4`$ and
$`G_U=O(n^2)`$, with the same asymptotic actual-inverse ledger. All equalities
hold on the complete n-qubit space and after tensoring an arbitrary reference.

<a id="two-outer-conjugations"></a>

## 4. The exact two-outer-conjugation search

Let Q be any traceless Hermitian target involution, including the complete
lookup phase tensored with identity on all arbitrary helpers. Fix candidate
Clifford walls $`C_2,C_3`$ and diagonal weighted cubic phases $`D_2,D_3`$,
and put

```math
A=C_3D_3C_2D_2,\qquad K=A^\dagger Q A.                 \qquad(7)
```

There exists a source $`V=AU`$ with U third-level and
$`VPV^\dagger=Q`$ if and only if K is a Clifford. Necessity follows from
the definition of the third hierarchy level. For sufficiency, K is
automatically Hermitian, involutory and traceless, so (3)--(6) construct U
with at most four native bands and no additional work. In particular the
first source can always be chosen diagonal third-level in fixed Clifford
coordinates. This replaces only the observable image of U: the four-band source need
not approximate U as a full operator. Its complete derivative is identical.
The same image-set equality applies at finite error. With $`\mathcal H_n`$
the set of traceless Hermitian n-qubit Cliffords,

```math
\inf_{U\in\mathcal C_3}\|AUPU^\dagger A^\dagger-Q\|
 =\min_{L\in\mathcal H_n}\|A^\dagger Q A-L\|.         \qquad(8)
```

Every right-hand approximant has a four-band source with exactly the same
complete error. The balance and Hermiticity restrictions on L are essential;
an arbitrary nearby Clifford need not be a Pauli image. This is an image-set
identity for fixed exact outer factors, not an approximation-rounding theorem
for unrestricted moving-key walls.

If each outer nonlinear phase is a cubic-sign batch paid in eight bands,
the resulting source has at most twenty T bands and T count at most
$`3n-2+8(r_2+r_3)`$, where $`r_j`$ counts its cubic occurrences. General weighted outer
phases instead use their own paid compiler; the thirteen-band compiler
would give at most thirty bands. Clifford coordinates, phase occurrences,
private returned pools and all inverse costs remain in the count and width
ledger. The source supplied by (3) has no private helpers of its own. Here n counts
the principal space on which K acts; returned private phase-compiler helpers
can be appended afterward, with the target and source extended by identity.
They are charged to the outer batches, not silently treated as clean inputs.
For physical width $`n=2^p\mathrm{poly}(p)`$, its count and Clifford cost
fit the polynomial-overhead discovery budget. This does not certify that
the chosen outer phases or walls fit that budget or satisfy (7).

The earlier two-block and later-key-preserving three-block obstructions
continue to apply in their stated domains. Equation (7) permits arbitrary
late key mixing, but supplies neither those outer factors nor a shallow
lookup. The remaining constructive problem is to find a complete target
conjugation whose peeled image is Clifford, with qualifying native costs.
The [exact fixtures](../../tests/test_clifford_involution_source.py) check
quadratic compression, all columns of the four-band word, literal phases,
actual inverses and independently specified nonmonomial Clifford examples.
