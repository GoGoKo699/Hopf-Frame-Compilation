# Disjoint-axis two-band field decoders have a fixed error floor

[Field query and full contract](FIELD_SHEAR_INTERTWINER.md) ·
[One-band axis matching](ONE_BAND_AXIS_MATCHING.md) ·
[Four-address decoder](FLAT_FIELD_BANK_PILOT.md)

A full field-to-wire decoder cannot be made from two physical T bands
whose absorbed axis spans have zero intersection. The bands may be partial,
have arbitrary T/T-dagger orientations, and act jointly on arbitrary
helpers. The source field may be enlarged, and the decoded helper query
may have arbitrary Clifford action. Fixed Clifford walls may mix every register.
The full-intertwining operator error is at least

```math
c_2=\frac12\sqrt{1-\frac1{\sqrt2}}\simeq0.270598.
```

This is a restriction on the two-band arrangement, not on all two-band
words. Shared axis spans permit collisions not covered by the proof.
It does not apply to the weaker single-probe decoder contract, nor does
it give a growing depth lower bound for general frame compilation.

## 1. Scope and a stronger endpoint statement

Let p be at least two, $`m=2^p`$, and let $`Q_x`$ be a literal full-rank affine
field-shear family on 2K core qubits, with identity on h source helpers.
The field degree K is at least p; the p address scalars are independent.
Put n=2K+h and suppose n is at least m. Every nonzero address has the
same fixed Pauli-label space. Every other
Pauli label has a free address orbit of size m, contained in the binary
span of at most p+1 labels. These properties survive a fixed Clifford
change of coordinates and identity helpers.

Absorb the Clifford interlayers of a two-band unitary into

```math
E=L D_B D_A K,
```

where K,L are fixed Cliffords, the physical axes of each D are independent
commuting Hermitian Paulis, and their binary spans A,B satisfy
$`A\cap B=0`$. Either band can have inactive wires. Put
$`\widehat Q_x=KQ_xK^\dagger`$.

The stronger statement is: if every
$`D_BD_A\widehat Q_xD_A^\dagger D_B^\dagger`$ is within operator
distance strictly below $`c_2`$ of some projective Clifford, then
$`D_BD_A`$ commutes literally with every $`\widehat Q_x`$.
Consequently every endpoint is already the Clifford conjugate
$`L\widehat Q_xL^\dagger`$.

The full-intertwining bound applies to every target

```math
T_x=W_x\otimes R_x,
```

where $`W_x`$ translates m bank wires and $`R_x`$ is any Clifford on the
remaining n-m wires. The $`R_x`$ need not be identity; residual helper
queries that are later canceled by the complete probe echo are allowed.
Clifford similarity preserves the source identity $`N_xN_y=0`$ for the symplectic deviations
$`N_x=F_x-I`$. For wire translations on m bank wires, the deviations of
two independent address translations have nonzero product. The target
symplectic matrices are block diagonal between bank and residual wires,
so no residual block can cancel this bank product. Identity source helpers
also preserve the zero-product identity. Thus some source and target
endpoint symplectic actions differ. Such projective Clifford operators
have distance at least $`1/\sqrt2`$, exceeding $`c_2`$; see the separation
argument in [one-band matching](ONE_BAND_AXIS_MATCHING.md#3-exact-and-robust-proofs).

## 2. Normalize the two bands without adding a gate

Move every individual axis of $`D_A`$ that commutes with all axes of $`D_B`$
into the outer band. Denote the remaining inner band by $`D_U`$ and the
combined outer band by $`D_V`$. The moved rotations commute with both
bands, so this only reorders the original factors and preserves every
native scalar phase and T orientation.

The combined outer axes are commuting and independent because
$`A\cap B=0`$. The new spans U,V still have zero intersection.
Every remaining individual inner axis pairs nontrivially with at least
one outer axis. The reordered word remains two physical bands in fixed
Clifford coordinates, with exactly the original count and width.
If all inner axes moved, the statements below still hold with an empty
inner band.

## 3. Recover the entire outer axis list from two-term images

For a Hermitian Pauli P with label v, let r be the number of inner axes
anticommuting with P. Conjugation by $`D_U`$ produces a flat expansion on
$`v+U_v`$, where $`U_v`$ is the r-dimensional span of those active axes.
For each $`s\in U_v`$, let $`\beta(v+s)`$ be the binary activation vector
against the individual outer axes. Conjugation by the outer band then
has exactly

```math
\left|\mathrm{supp}(D_VD_U P D_U^\dagger D_V^\dagger)\right|
   =\sum_{s\in U_v}2^{\mathrm{wt}(\beta(v+s))}.          \qquad(1)
```

There are no collisions: equality between labels from different paths
would place a nonzero difference in $`U\cap V`$. A coefficient on the
branch s has magnitude
$`2^{-(r+\mathrm{wt}(\beta(v+s)))/2}`$.

If r is positive, toggling any active inner axis changes beta, by the
normalization in Section 2. For r=1 the two beta values are distinct,
so the sum in (1) is at least three. For r at least two the sum is at
least four. Thus a two-term output occurs exactly when

```math
r=0,\qquad \mathrm{wt}(\beta(v))=1.                     \qquad(2)
```

Every such output is a signed flat sum of two Hermitian Paulis with
coefficients of magnitude $`1/\sqrt2`$. Its two support labels differ
by the unique active outer axis. Conversely, each outer axis is
realized as such a difference: independence of the combined inner and
outer axes makes the commutation equations r=0 and beta equal to the
chosen unit vector simultaneously solvable.

Therefore the set of binary differences of all two-term Pauli images
is exactly the set of individual outer axes. This recovers the axis
list itself, rather than only its span. The input labels satisfying
(2) need not form a basis; no such assumption is used.

## 4. The support fingerprint is stable at fixed operator error

Use normalized Hilbert--Schmidt norm, for which Hermitian Paulis are an
orthonormal basis and a unitary has norm one. Every image that is not a
signed flat two-term sum has distance at least

```math
\delta=\sqrt{1-1/\sqrt2}=2c_2                         \qquad(3)
```

from every such sum. To see this, the overlap with a normalized flat
two-term sum is at most the sum of the two largest coefficient
magnitudes divided by $`\sqrt2`$.

If r=0 and the outer activation has weight other than one, the image
is either a single Pauli or has at least four equal-magnitude terms;
that overlap is at most $`1/\sqrt2`$. If r is at least two, every
coefficient has magnitude at most one half, giving the same bound.
For r=1, distinct beta values allow at most one coefficient of
magnitude $`1/\sqrt2`$, and every other coefficient is at most one
half. The overlap is at most $`(1+1/\sqrt2)/2`$, proving (3).

Let $`F=D_VD_U`$. Suppose a projective Clifford $`C_x`$ obeys

```math
\|F\widehat Q_xF^\dagger-C_x\|_{\rm op}\lt c_2.
```

For any two-term image $`FPF^\dagger`$, its $`C_x`$ conjugate is again a signed
flat two-term sum. Its distance from
$`F\widehat Q_xP\widehat Q_x^\dagger F^\dagger`$ is below
$`2c_2=\delta`$. Equation (3) forces the latter to be a two-term
image as well. Thus each source query permutes the intrinsic input set
of two-term images.

Moreover these two signed flat expansions must be exactly equal.
Two unequal signed flat two-term expansions have normalized
Hilbert--Schmidt distance at least one: their coefficient inner product
is at most one half. Their present distance is below delta, which is
less than one. Consequently $`C_x`$ sends the support differences of these
images onto themselves. Section 3 proves that $`C_x`$ permutes the
individual outer axes up to Hermitian signs.

## 5. Two applications of one-band matching close the argument

By the converse direction of one-band matching,
$`D_V^\dagger C_xD_V`$ is a projective Clifford. It is within $`c_2`$ of
$`D_U\widehat Q_xD_U^\dagger`$. Since
$`c_2\lt\sqrt2/4`$, robust one-band matching forces
$`\widehat Q_x`$ to permute the individual inner axes.

A nonfixed field Pauli orbit has m distinct labels in a span of
dimension at most p+1. Because $`m\gt p+1`$ for p at least two,
such an orbit cannot lie in an independent axis list. Every inner
axis is therefore fixed as a label by every query. A fixed field
Pauli label is fixed with positive literal sign, including after the
fixed Clifford coordinate change. Hence $`D_U`$ commutes with every query.

The nearby endpoint now has the form
$`D_V\widehat Q_xD_V^\dagger`$. Robust one-band matching applies
again and forces the query to permute the outer axes. The same orbit
argument fixes every outer axis literally. $`D_V`$ also commutes with every
query, proving the endpoint statement in Section 1.

The proof is independent of the source field degree and helper count,
subject to the dimensions in Section 1. The decoded residual action
can be Clifford and nontrivial. All norms concern complete operators,
so coherent addresses, dirty inputs and external references are retained. The overlap assumption on the
two absorbed spans is essential to this proof. A candidate outside it
must exploit actual shared-span support collisions; merely using
partial independent bands does not evade the result.
