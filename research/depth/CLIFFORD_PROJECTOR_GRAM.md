# Clifford projector-Gram rigidity and affine-query limits

[Accepted-source contract and recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) ·
[Invertible Clifford-fiber compression](GAUSSIAN_COMPRESSION_SOURCE_LIMITS.md) ·
[Bell-envelope acceptance](CORANK_FOURIER_BELL_SOURCE.md#bell-acceptance-charge-transfer)

This chapter concerns singular accepted Grams of complete Clifford fibers,
including fully interleaved data and auxiliary action. A strict quarter-error
promise forces the exact target Gram and at least two clean signals. A
separate degree argument limits the number of completed affine-address
queries. The first statement does not assume a Bell envelope; the second
is a query-grammar restriction, not a native T-depth lower bound.

The rigidity theorem applies to any nonidentity signed Pauli target.
For its lookup corollary, use p address bits, $`m=2^p`$, and $`n=pm`$ bits in each
carrier C, D, Y and Z. View the p-bit address as $`x\in\mathbb F_{2^p}`$.
Expand field elements in the fixed binary basis whose first element is one,
and define

```math
\begin{aligned}
v_x&=(1,x,\ldots,x^{m-1})^{\mathsf T},\\
Q_x&=X_C(v_x)X_D(v_x),& R_x&=Z_Y(v_x)Z_Z(v_x).
\end{aligned}
```

These are the complete charges from the
[corank-one Fourier construction](CORANK_FOURIER_BELL_SOURCE.md#1-the-singular-fourier-operator-and-its-unpaid-dilation).
The query count q below is distinct from the carrier Hilbert-space dimension.

<a id="interleaved-clifford-gram-rigidity"></a>

## 1. Interleaved Clifford fibers still require a growing query word

Leaving the Bell envelope does not by itself remove the charge-readout
obligation. There is a broader test for sources whose net action fixes the address
and is Clifford on each address fiber. It permits arbitrary interleaving of
data, auxiliary banks and clean signals; it does not assume the [Bell-envelope factorization](CORANK_FOURIER_BELL_SOURCE.md#bell-acceptance-charge-transfer). The earlier
[Gaussian compression result](GAUSSIAN_COMPRESSION_SOURCE_LIMITS.md) assumes
an invertible accepted block near a scaled unitary. The present statement
instead treats the singular projector-Gram promise needed by recovery.

Let $`V_x`$ be Clifford on an arbitrary bank B and k clean signals A, with
$`k\ge1`$. B includes every logical carrier and every dirty helper; all its
inputs, including reference entanglement, are retained. Write

```math
J=I_B\otimes|0^k\rangle_A,\qquad
K_x=J^\dagger V_xJ,\qquad G_x=K_x^\dagger K_x.
```

For any nonidentity signed Hermitian Pauli $`Q_x`$ on B, the strict promise

```math
\left\|G_x-\frac{I+Q_x}{4}\right\|\lt\frac14
\quad\Longrightarrow\quad
G_x=\frac{I+Q_x}{4}.                                    \qquad(36)
```

Moreover k must be at least two. With one clean signal the error is at least
one half, regardless of bank width or the length of the Clifford word. These
are statements about the net Clifford fibers; intermediate address motion
is permitted if that final structure is established. They do not cover an
arbitrary native source with one clean signal.

### Stabilizer compression proves exact Gram rigidity

Suppress x and put $`H_j=V^\dagger Z_{A_j}V`$. These are commuting signed
Hermitian Paulis. Let M be the k-by-k binary matrix whose j-th column is the
**input-clean X label** of $`H_j`$, and set $`\rho=\mathrm{rank}M`$. Thus M is
taken from the inverse tableau, with no implicit transpose of a forward
clean block. Expanding the output-clean projector gives

```math
G=2^{-k}\sum_{t\in\ker M}R_t,\qquad
R_t=J^\dagger\left(\prod_jH_j^{t_j}\right)J.           \qquad(37)
```

Terms outside the kernel vanish. Inside it, the clean Z factors act as plus
one, and the $`R_t`$ form a commuting signed Pauli group on the whole bank B.
If this group contains negative identity, its sum is zero. Otherwise its
average is a stabilizer projector $`P_S`$ and

```math
G=2^{-\rho}P_S.                                        \qquad(38)
```

If its independent data-stabilizer rank is ell, then
$`\ell\le k-\rho`$ and $`P_S`$ has rank $`2^{-\ell}\dim B`$.
Possible clean-only kernel relations change ell, not the factor in (38).

The target in (36) has eigenvalues zero and one half, each on half the bank.
A nonzero eigenvalue one in (38) is at least one half from the target spectrum;
if rho is at least two, the largest eigenvalue of G is at most one quarter.
The strict error promise therefore forces rho=1. It then gives
$`\|P_S-(I+Q)/2\|\lt1/2`$, so these projectors have equal rank and ell=1.
Consequently $`P_S=(I+Q')/2`$ for a nonidentity signed Pauli Q'. Distinct
signed Hermitian Paulis have operator distance at least square root two,
which contradicts (36) unless Q'=Q. This proves the complete equality.

For k=1, (37)--(38) allow only zero, identity, I/2, or a Pauli half-projector.
Each is at distance at least one half from the target: a nonzero projection
has expectation one on some vector whereas the target is at most one half,
and I/2 has expectation one half on the target's zero sector. Zero attains
this bound. For two signals, H on both signals instead gives G=I/4 and error
exactly one quarter. Thus the strict threshold in (36) cannot be enlarged
uniformly to include equality.

Rigidity does not prohibit the desired Gram abstractly. For a fixed Pauli
$`Q=Z_b`$, CNOT from b to a clean signal a, followed by H on a second clean
signal c, gives $`K=(I+Q)/(2\sqrt2)`$ and exactly the required Gram.
A Clifford change of bank basis handles any other fixed nonidentity Q.
These are literal Clifford words on arbitrary bank inputs. Choosing such
a basis coherently for the entire address-dependent $`Q_x`$ family still needs
its own native construction and complete cost/return proof.

### A query-count consequence for bounded clean allocation

Suppose $`V_x`$ is a word of q complete affine-address XOR shear queries, in
either direction, interleaved with arbitrary fixed address-independent
Cliffords on all bank and signal wires. The actual queries and their inverses
are retained; no principal or helper bank is initialized. Each fixed-address
query is Clifford with an affine binary symplectic matrix. Matrix
multiplication shows that every physical label entry has Boolean address
degree at most q. The inverse word has the same q queries, so the same degree
bound applies to the labels in (37). This is the physical-label argument of
the [Gaussian query grammar](GAUSSIAN_COMPRESSION_SOURCE_LIMITS.md), without
its invertible-block hypothesis.

More generally let D bound the degrees of those inverse-tableau labels.
Under (36), rho=ell=1. For one bank Pauli-label coordinate j, let $`c_j`$ be the
length-k row of that coordinate in the $`H_i`$ labels. The unique nonzero data
stabilizer label has coordinate

```math
q_j(x)=1\quad\Longleftrightarrow\quad
c_j(x)\notin\mathrm{rowspan}\,M(x).                     \qquad(39)
```

Indeed, $`c_j`$ has nonzero restriction to ker M precisely in this case, and
the image of ker M has dimension one. The membership test is a Boolean
function of the k-squared entries of M and k entries of $`c_j`$. Extending it
arbitrarily on invalid tuples gives a Boolean polynomial in these
$`k(k+1)`$ bits, hence

```math
\deg q_j\le k(k+1)D.                                  \qquad(40)
```

This bound is independent of the number of dirty wires. For k=2 it improves
to 3D. Write $`w_i=\mathrm{OR}(M_{1i},M_{2i})`$ for each of the two columns.
Since a rank-one nonzero matrix factors as $`M=u w^{\mathsf T}`$ with u
nonzero, this recovers its unique nonzero row w. Its kernel has the unique
nonzero vector $`t=(w_2,w_1)`$, so

```math
q_j=c_{j1}w_2+c_{j2}w_1,\qquad
\deg w_i\le2D,\qquad \deg q_j\le3D.                   \qquad(41)
```

The Bell target $`Q_x=X_C(v_x)X_D(v_x)`$, and its auxiliary charge $`R_x`$,
have a Pauli-label coordinate of degree p: in the fixed binary field basis
with first element one, the coordinate $`x^{2^p-1}`$ is the Boolean
nonzero-address indicator, including its value zero at x=0. Its algebraic
normal form has degree p. Thus any such q-query word satisfying the complete
strict Gram promise at every address must obey

```math
k(k+1)q\ge p;\qquad k=2\ \Longrightarrow\ 3q\ge p.
                                                               \qquad(42)
```

The fixed Vandermonde bank change does not evade this conclusion: an
invertible constant binary map preserves the maximum degree among vector
coordinates. Nor can a sign convention remove this label obstruction.

This excludes a constant number of these completed affine queries with a
bounded clean allocation, even with fully interleaved Clifford data/auxiliary
mixers. It is a query-word restriction, **not a native T-depth lower bound**.
The paid query macros retain their previously established native counts,
widths and returns; no lower bound on their jointly resynthesized depth follows
from (42). Non-Clifford fibers and intermediate address motion are outside
this query grammar. Growing clean allocation weakens (42) and is not the
bounded-clean source contract.

If an implementation is within eta of such an ideal V on the full space,
its accepted Gram is within $`2\eta`$ of G. Therefore the same exclusion
applies when its requested Gram error delta obeys
$`\delta+2\eta\lt1/4`$. For individually approximated query calls, their
full errors sum in eta. This keeps every dirty/reference return error in
the comparison and does not grant a free ideal query.

The [exact projector-Gram fixtures](../../tests/test_clifford_projector_gram.py)
retain all bank columns, literal Clifford inverses, coherent dirty input,
the strict quarter boundary, and the two-signal kernel formula. They check
these identities in finite dimensions; the all-size rigidity and query-count
statement follow from the argument above. The next accepted-source attempt
must change the fiber dynamics or the address-motion mechanism, rather than
add another fixed Clifford mixing layer around the same bounded query word.
