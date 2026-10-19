# Revision: global native matching and the exact readout cost

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Quadratic-control matching](QUADRATIC_CONTROL_NATIVE_MATCHING.md) ·
[Exact dual-readout barrier](DUAL_FRAME_READOUT_BARRIER.md) ·
[Previous checkpoint](COMPUTATIONAL_PHASE_SESSION_RESULT.md)

This continuation began from verified main `4e83fe1` (PR #191) on
19 October 2026. It kept the native lists global while studying their
shared center, and tested whether the separate Walsh readout could meet
the query budget under unrestricted native resynthesis. Both questions
now have stronger answers. The complete-frame resource frontier is
unchanged.

## 1. Destination, general theorems and components

The destination remains optimal general depth for the prescribed complete
Hopf differential frame, on one circuit meeting every count, width,
Clifford and initialized-isometry error requirement.

| Role | Current conclusion |
| --- | --- |
| Strongest general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ |
| Literal domain | Integers $`n\ge1,L\ge6`$; one clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full error at most $`2^{-L}`$ |
| Selected manuscript theorem | R79 remains Result D; the selected count frontier is settled at its stated allocations |
| Unrestricted low-layer results | At most one clean and two native layers give XOR degree at most two; zero clean and three layers give sign/XOR bounds five/four |
| Retained fixed-rank result | Both fixed endpoint clean-active ranks admit a uniform bound, including arbitrary middle pairs |
| New lower components | Within the retained no-pair normal form, fixed quadratic polar rank closes another subcase; arbitrary genuine quadratic control forces a common invariant middle span and monomial Clifford endpoint images after paid normalization |
| New source obstruction | Exact standalone Walsh dual readout requires $`\Omega(m\sqrt{\log m})`$ T gates, with arbitrary Clifford interlayers and returned helpers |
| Unresolved | Growing-rank one-clean three-layer degree, a faster complete query meeting every budget, and matching general depth bounds |

At fixed accuracy and square-root dirty width, the general upper bound
is still logarithmic in n and the unrestricted lower bound is constant.
Neither a restricted degree theorem nor a count obstruction for a
separate operation closes that interval.

## 2. The actual quadratic insertion and two closed components

The retained no-middle-pair reduction, with both interior clean images
non-Pauli, gives the full logical word

```math
V=K_C^\dagger E_C^\dagger H^{b(z)}E_BE_AK_A=e^{i\theta}D_g.
```

The E lists are independent native Pauli T/T-dagger bands. The genuine
computational center R is fixed by every factor. H commutes R and the
middle axes, and b is quadratic. The extracted central cubic and real
quadratic signs remain in the literal target factorization.

If the polar rank of b is $`2r`$, its affine control is an actual
Clifford and can be absorbed into the output band and endpoint. The
remaining insertion has Pauli support-space dimension at most
$`2r+1`$. Finite insertion gives
$`\deg g\le D_{\mathrm{ins}}(0,0,0,2r+1)`$ for positive r.
For r=0, actual Clifford absorption leaves a zero-clean three-band word
and gives sign/XOR bounds five/four. These conclusions are uniform in
growing endpoint ranks when polar rank stays fixed.

A separate rigidity lemma applies to a full Clifford that is diagonal
on entire quantum fibers over an affinely spanning set of center values:
it is diagonal globally. A quadratic control of polar rank at least four
has such support. Consequently a diagonal matching-band sandwich
$`E^\dagger H^{b(z)}E`$ forces E to commute H. A literal rank-two
example shows why the support premise matters. The actual word has
different endpoint bands; this lemma alone does not identify them.

Raw polar rank is not a universal hardness parameter. An actual
three-band family with both endpoint ranks 2m and no middle pairs has
control polar rank 2m but target $`(-1)^{h\sum_jx_jy_j}`$, of degree
three. It has a shared clean flip and is already in a closed branch.
It does not witness unbounded polar rank after excluding that branch,
nor minimize rank over every equivalent representation.

## 3. Progress without restricting to a center fiber

For genuinely quadratic b, let J_u and C_u be the actual one-band
Clifford images of a computational probe at the input and output.
Diagonality gives

```math
C_u=H^{b(z)}E_BJ_uE_B^\dagger H^{b(z)}.
```

The new global argument first proves that a middle-axis image under
J_u activates no middle native axis. The non-Pauli quadratic character
cannot be hidden in that one-band image: its displacement would lie in
R, while native independence and commutation with R forbid any activated
axis there. Squaring the resulting commuting defect then gives an exact
Gram identity between the original middle list and its J_u image.
Independence on the full register forces

```math
J_uL_B=L_B.
```

Both probe families preserve
$`S=R+L_B+\langle h\rangle`$ and agree literally on its whole Pauli
algebra. This is span preservation, not permutation of the native axes;
distinct independent Gram bases can leave a nontrivial cubic phase.

The joint displacement span

```math
R_0=\mathrm{span}\{(\Phi_{J_u}-I)s:u\text{ arbitrary},\ s\in S\}
   \subseteq S\cap L_A\cap L_C
```

commutes with all retained native bands. Thus $`R+R_0`$ is a larger
actual shared center. The existing Boolean-reality argument pays its
Clifford normalization, with Hadamards only on exact identity modes.
In these coordinates the endpoint images of S are commuting monomial
Cliffords of the form

```math
|z,y\rangle\longmapsto
\gamma_s i^{\ell_s(z,y)}(-1)^{q_s(z,y)}
|z,y+a_s+L_sz\rangle.
```

The input and output families have the same computational permutation.
Their quadratic work-coordinate phases remain present. No center fiber
has been selected, no dependent projected list has been credited as an
independent native band, and no extra phase has been deleted.

## 4. The first unsupported inference and its actual scope

The monomial algebra is not automatically a Pauli algebra with a free
Clifford decoder. Gaussian elimination selected separately for each z
can have changing rank and requires a globally implemented conditional
operation. The proof supplies neither that decoder nor a bound on the
relative center phases. A direct global phase argument could avoid a
decoder; it must still use the native factorization.

There are two precise limits on shortcuts. Multiplying the target by an
arbitrary sign of the unresolved common endpoint parities preserves the
recorded central transport and covered derivatives, but need not preserve
three native bands. This proves insufficiency of those equations alone,
not a high-degree native counterexample. Sliding the insertion through
an activated endpoint Pauli replaces it by a Pauli times
$`(-1)^{b(z)}`$, generally a Clifford involution rather than a Pauli.
Its spectral T rotation is not a native Pauli-axis T gate. Calling the
conjugated list a native band leaves that operation unpaid.

The selected next lower question is therefore whether the globally
independent endpoint lists constrain the relative phases of these
matched monomial Clifford algebras, uniformly in the center size.
The new invariant span and enlarged center are the completed first
steps toward that question. A representation-dependent decoder failure
does not obstruct the unrestricted model, and a closure here would
still need the remaining growing-rank branches before giving a general
one-clean three-layer degree theorem.

## 5. A true unrestricted obstruction for the separate exact readout

For the Walsh dual readout Q on $`m=2^p`$ modes, the new bound is

```math
T(Q)\ge2L_p,\qquad
L_p=p\binom{p-1}{\lfloor(p-1)/2\rfloor}
   =\Theta(m\sqrt p).
```

An explicit downset Walsh minor has determinant
$`\pm2^{-L_p/2}`$. In Q's Pauli channel, its square occurs as a
unique coefficient with a distinct dirty-table character. Every physical
T gate increases the square-root-of-two denominator exponent by at most
one. Returned dirty work preserves the coefficient, and exact clean
return contracts channel entries by an integer sum, without a helper
normalization penalty. The proof permits arbitrary Clifford interlayers
and any finite returned clean allocation. It also survives a right unitary
gauge on the table register, because the input witness is carrier-only.

This excludes linear T-count for every exact implementation of that
standalone readout, including non-Gaussian implementations. It does not
exclude the complete indicator: R91 already implements the m-row,
m-output identity table with O(m) T gates, O(m) dirty work,
$`O(m^2)`$ Cliffords and $`O(\log(p+2))`$ depth. Joint compilation
can therefore avoid a real cost of the separated mechanism.

A new six-call wrapper also shows that a complete Q alone would suffice
for the indicator; a separate flat-source premise is unnecessary.
The wrapper charges all calls, address masks, actual adjoints and work.
It tolerates a diagonal table-phase gauge, which is narrower than the
right-unitary gauge allowed in the lower bound.

The obstruction is exact. Its very small coefficient gives no uniform
fixed-error or inverse-polynomial-error lower bound. An approximate Q
with full operator or initialized-isometry error tau gives indicator
error at most $`6\tau`$, including work-return error. Its full native
resource ledger and an approximate version of the currently exact R98
frame transfer remain unpaid. The next upper question should use a
jointly compiled complete query or a fully accounted approximate source,
not resynthesize this exact standalone Q at impossible linear T-count.

## 6. Verification and continuing contract

Both proof notes were independently audited. Seven new exact fixtures
check the paid affine-control absorption, the rank-two sandwich boundary,
the actual all-column growing-rank family, Walsh downset determinants,
a native channel coefficient, the six-call indicator with scalar gauge
and actual inverses, and clean-projector contraction. The arbitrary-width
matching and count theorems remain analytic, rather than consequences
of small numerical examples.

Every logical column, literal phase, arbitrary dirty input and reference,
actual inverse, clean leakage and work-return condition remains in scope.
The R96 frontier and selected manuscript claims are unchanged. Final
manuscript drafting remains separate from this research continuation.
