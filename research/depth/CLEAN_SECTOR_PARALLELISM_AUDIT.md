# Clean-sector parallelism: what the native schedule actually buys

[Clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Current frontier](../../docs/OPEN_PROBLEM.md)

Parallelizing the clean-sector cover looks like a way to permit more
source-clean qubits without multiplying depth by their number of sectors.
There is a sound native schedule under a strong read-only hypothesis.
That hypothesis also lets us remove the clean bank directly, at no
T-count or T-depth cost. It therefore supplies no new shallow source.

A separate rank argument shows that a universal zero-clean black-box
wrapper really needs exponentially many source occurrences in the number
of clean bits. This is an occurrence bound, not a native-depth bound;
source inspection and joint resynthesis are outside its hypothesis.
Neither assessment changes R96 or the unrestricted depth frontier.

## 1. A clean bank fixed at every native boundary is dispensable

Let A contain k source-clean qubits, and let R contain all remaining
physical wires: address, arbitrary output and arbitrary dirty work.
Suppose the actual emitted source word is

```math
V=C_d B_d C_{d-1}\cdots C_1B_1C_0,
```

where every actual Clifford boundary fixes each computational A sector,
and every physical T or T-dagger in the bands B acts on R. Equivalently,
each C commutes with every single-wire $`Z_{A_j}`$. Individual gates
inside a Clifford boundary need not preserve A. The hypothesis concerns
the actual completed boundary, including both endpoints.

Write $`C=\bigoplus_a C_a`$. Each restricted unitary C_a is itself a
Clifford on all of R, including its literal scalar. To see the Pauli
claim, conjugate any R Pauli P by C. The resulting physical Pauli
commutes with every $`Z_{A_j}`$, and hence has the form
$`Z_A^z\otimes P'`$ with its actual Pauli phase retained. Restriction
to A=a gives a Pauli on R. This map cannot send a nonidentity R Pauli
to a scalar: conjugating back a Pauli supported only on Z_A leaves that
same Z_A Pauli, because C fixes all those generators. Thus C_a normalizes
the full R Pauli group.

The literal scalar can be retained by ordinary Clifford tableau
elimination with phase tracking. The remaining scalar is an eighth root
of unity; $`(HS)^3=e^{i\pi/4}I`$ implements its generator using only
Cliffords. This is not a free global-phase correction. On n R wires,
the resulting literal restricted Clifford has an O(n squared) elementary
H, S and CNOT word, with the same order of cost for its actual inverse.

Consequently the native word

```math
V_0=C_{d,0} B_d C_{d-1,0}\cdots C_{1,0}B_1C_{0,0}
```

uses no A register at all. It has the same T-count and no greater
T-depth than V. If J inserts A=0 and the source satisfies the complete
initialized-isometry promise

```math
\|VJ-JQ_f\|\le\epsilon,
```

then exact A-sector preservation gives the full-operator promise

```math
\|V_0-Q_f\|\le\epsilon
```

on every R input, including all dirty-work and reference correlations.
The actual inverse has the same error. The simultaneous ledger is

```math
T_0=T_V,\qquad d_0\le d_V,\qquad
G_0=O((d_V+1)n^2),\qquad w_0=w_V,\qquad a_{\rm clean}=0.
```

Here w_V excludes the declared address, output and clean A registers.
No source invocation, zero-test, sector cover or error amplification is
needed. If the source address is also fixed at every boundary, its
restriction remains read-only.

For the relaxed b-bit source target, n of size
$`O_\delta(2^b(b+2)^c)`$ and qualifying
$`d_V=o(\log(b+2))`$ make this Clifford overhead
$`O_\delta(4^b(b+2)^{c_G})`$ for a fixed exponent c_G.
Thus replacing the clean bank by these restricted Cliffords stays inside
the root-fusion budgets. This is an exact compiler reduction, not a
claim that an arbitrary source with clean return satisfies its hypothesis.

## 2. Why the proposed parallel cover is dominated by restriction

Put $`M=2^k`$. Give each sector a its own old output bank, toggle,
source-private work, zero-test work and controlled-phase helpers. Share
only the address x, the arbitrary sector register A and the final output
Z. In the phase basis, assume every actual source Clifford boundary fixes
x and A, and all source T targets are private.

The translated source $`X_a V X_a`$ can be emitted without changing
shared A between boundaries: replace each actual boundary C by
$`X_a C X_a`$. The intervening X pairs cancel because every T band
acts outside A. This retains literal phases and pays the conjugations.
The zero-test has the same read-only interface. A fully read-only dirty
CCZ word pays each controlled phase on shared Z using a private helper;
all its T targets are private too.

At a fixed shared x,A,Z sector, the separate Clifford blocks then act
on disjoint private registers, times commuting scalar phases. Hence
the blocks commute and can be combined between common physical T bands.
This proves parallelism under the stated *native-boundary* hypothesis,
not just from commutativity of the completed selective operators.

The exact-sector estimate in
[clean promotion Section 4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#4-approximate-clean-return-still-gives-a-full-operator-query)
also applies with separate private banks. Each selective word is exactly
identity off its own A sector. On that sector it is
$`R_a^\dagger B_a R_a B_a`$, within $`2\epsilon`$ of the
desired phase. The full cover is a direct sum, so its error is at most
$`2\epsilon`$, with all private inputs arbitrary. The final Hadamards
on Z preserve that full-operator estimate. The sector register, toggles,
zero-test work and exact controlled-phase helpers return exactly;
possible source-private return errors are included in the norm bound.

This sound parallelization uses M private copies. More importantly,
Section 1 already deletes the source-clean A bank under the very same
native hypothesis, preserving depth and error without those copies.
The apparent new growing-clean window therefore does not enlarge the
sufficient source class within the polynomial-overhead budgets.

If only the *completed* V returns or preserves A, a Clifford boundary
may move A into private work or mix it coherently with that work. The
different sector copies then share a genuinely used quantum register.
Their local native blocks need not commute, and their T layers cannot
be interleaved on the evidence of completed selective identities alone.
Restricting such an intermediate Clifford to A=0 need not even yield
a unitary on R. The next unpaid operation is a native parallel schedule
for this genuinely used coherent clean work, with its full return
certificate. Approximate clean return does not supply that schedule.

## 3. A universal black-box wrapper needs at least all sector occurrences

Consider one fixed zero-clean unitary compiler using h occurrences of
a supplied source or its actual inverse. Its interlayers, wire locations
and query pattern are independent of the supplied source. They may be
arbitrary unitaries, move the address and use any finite arbitrary dirty
space. They must produce the literal complete promoted query on the
entire physical space. No initialized subspace is available to the
wrapper. The source itself has k promised clean bits A and at least one
old output bit t.

Compare the two source oracles

```math
V^{(0)}=I,\qquad
V^{(1)}=P_A X_t+(I-P_A)I,\qquad
P_A=|0^k\rangle\langle0^k|_A\otimes I.
```

They satisfy exact complete clean-code promises for the constant queries
zero and one, respectively. Their actions off the promised sector are
also explicit. Both sources are Hermitian involutions, so choosing an
actual-inverse occurrence does not alter this pair of examples.

Let D be the dimension of the compiler's entire physical Hilbert space,
including every arbitrary helper. Each embedded source difference has

```math
\mathrm{rank}(V^{(1)}-V^{(0)})=\frac{D}{2M}.
```

Indeed P_A selects a fraction 1/M of the space, and X_t minus identity
has rank one on the old output qubit. Tensoring arbitrary spectators,
or conjugating by any fixed wire layout, leaves that rank fraction
unchanged. Product telescoping and rank subadditivity therefore give,
for the two complete emitted unitaries,

```math
\mathrm{rank}(W^{(1)}-W^{(0)})\le\frac{hD}{2M}.
```

The literal promoted targets are identity and X_z on a new output bit,
both tensored with identity on every other physical wire. Suppose their
full-operator errors are tau_0 and tau_1, with
$`\tau_0+\tau_1\lt2`$. Then

```math
\|(W^{(1)}-W^{(0)})-(X_z-I)\|\le\tau_0+\tau_1\lt2.
```

The matrix X_z minus identity has D/2 singular values equal to two.
Any matrix of rank less than D/2 has a kernel vector in its minus-one
X_z eigenspace, and on that unit vector the displayed difference has
norm two. Thus the rank of the emitted difference is at least D/2,
which proves

```math
h\ge M=2^k.
```

This is a total source-occurrence bound. If many physical oracle calls
can actually run in one layer, they still count separately toward h.
The result does not lower-bound that layer count or unrestricted native
T-depth. It also does not apply when the compiler inspects the source
description, synthesizes restricted Clifford blocks as in Section 1,
uses source-dependent interlayers, or jointly resynthesizes a completed
word without retaining black-box source calls. The literal target phase
and full-space error are part of the theorem; they are not replaced by
a projective or initialized-sector comparison.
