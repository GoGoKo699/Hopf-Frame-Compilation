# A robust obstruction to two-layer dirty source circuits

[Source primitives](OPERATOR_SOURCE_COMPILER.md) · [Restricted source depth](SOURCE_T_DEPTH.md) · [Current frontier](OPEN_PROBLEM.md)

The geometric source cannot be approximated arbitrarily well by two T
layers, even with unrestricted Clifford interlayers and arbitrarily many
returned dirty helpers. For source width at least five, the operator-norm
error is at least $`1/16`$. The same bound holds for its controlled version
at width at least four. These are constant-depth obstructions for
full-input source primitives, not growing depth lower bounds or lower
bounds for an initialized-clean complete-frame compiler.

The proof combines standard stabilizer overlap quantization with the
source's exact Pauli-transfer coefficients. The overlap fact is given in
[Aaronson–Gottesman, *Improved Simulation of Stabilizer Circuits*,
arXiv:quant-ph/0406196v5, Section III, final paragraph, p. 5](https://arxiv.org/pdf/quant-ph/0406196v5).
The observation that one T layer conjugates Paulis to Hermitian Cliffords
is also used explicitly in
[Zhang–Zhang, arXiv:2409.13809v2, Theorem III.1, equations (10)–(11)](https://arxiv.org/html/2409.13809v2#S3.SS1).
The transfer-entry reduction and source witnesses below specialize these
inherited facts; no general priority claim is made.

## 1. A discrete Pauli-transfer alphabet at T-depth at most two

Let V be a circuit on q physical qubits with T-depth at most two. Every
wire is included in its full unitary action. Clifford blocks may be
arbitrary, and a T layer may contain T or T-dagger on any set of distinct
wires. For Hermitian Pauli strings P,Q, define

```math
\widehat V_{P,Q}=2^{-q}\mathrm{Tr}(PVQV^\dagger).
```

**Lemma.** Every such coefficient belongs to

```math
\mathcal A_q=\{0\}\cup\{\pm2^{-j/2}:0\le j\le2q\}.
```

*Proof.* Write $`V=C_2D_2C_1D_1C_0`$, where each C is Clifford and
each D is one physical T layer; identity layers are allowed. Absorb the
outer Cliffords into the two Pauli strings, obtaining P',Q'. Cyclicity
of the trace gives

```math
\widehat V_{P,Q}
=2^{-q}\mathrm{Tr}\!\left[
 (D_2^\dagger P'D_2)
 C_1(D_1Q'D_1^\dagger)C_1^\dagger\right].
```

Conjugating a one-qubit Pauli by T or T-dagger gives a Hermitian Clifford:
X and Y become signed combinations $`(X\pm Y)/\sqrt2`$, while I,Z
remain Paulis. Tensor products and Clifford conjugation preserve this
property. The two displayed factors are therefore Hermitian Clifford
unitaries. Their product K is Clifford, and its trace is real because
the trace of a product of two Hermitian matrices is real.

For the normalized q-pair Bell state $`|\Omega_q\rangle`$,

```math
2^{-q}\mathrm{Tr}K
=\langle\Omega_q|(K\otimes I)|\Omega_q\rangle.
```

Both vectors are stabilizer states on 2q qubits. Their overlap magnitude
is either zero or $`2^{-j/2}`$ for an integer $`0\le j\le2q`$.
This follows directly by expanding their stabilizer projectors: the
squared overlap is zero if the common Pauli constraints conflict, and
otherwise is $`2^{s-2q}`$, where the common stabilizer subgroup has size
$`2^s`$. Reality of the trace gives the stated signs. In particular, no
assumption about the phase of a general Clifford trace is needed. ∎

This is a statement about individual full-space transfer entries. It
does not assert that a state-dependent Pauli expectation in a circuit
with two T layers belongs to the same alphabet.

## 2. Converting a transfer gap into robust operator distance

Put

```math
\mathcal A=\{0\}\cup\{\pm2^{-j/2}:j\ge0\}.
```

For any two unitaries U,V and any Hermitian Paulis P,Q on the same
physical space,

```math
|\widehat U_{P,Q}-\widehat V_{P,Q}|
\le\|UQU^\dagger-VQV^\dagger\|
\le2\|U-V\|.
```

The first step uses the normalized trace bound, and the second expands
the conjugation difference into two terms. Therefore, if one transfer
entry of U is x, every circuit V of T-depth at most two satisfies

```math
\|V-U\|\ge\frac12\mathrm{dist}(x,\mathcal A).
```

Tensoring U and both witness Paulis with identity preserves x. The
bound consequently survives any number of dirty helper wires, including
arbitrary intermediate mixing with them. The target is the full unitary
$`U\otimes I_{\rm dirty}`$; even approximate helper return is included
in this operator-norm distance. The same statement holds after optimizing
over a common global phase, since transfer coefficients ignore that phase.

## 3. Geometric source witnesses

For the existing geometric source

```math
M_m=\sum_{j=0}^{m-1}a_j\Gamma_j,\qquad
\Gamma_j=Z_0\cdots Z_{j-1}X_j,
```

write $`w_j=a_j^2`$. The weights are
$`w_j=2^{-j-1}`$ for $`j\lt m-1`$, with the duplicated last weight
$`w_{m-1}=2^{-(m-1)}`$. For any index subset S, the Pauli
$`P_S=\prod_{j\in S}Z_j`$ flips exactly those Majorana terms indexed
by S. Orthogonality of distinct Pauli strings gives

```math
2^{-m}\mathrm{Tr}(P_SM_mP_SM_m)
=1-2\sum_{j\in S}w_j.
```

At $`m\ge4`$, the single-qubit witness $`P=Z_2`$ has coefficient
$`3/4`$. Its nearest alphabet point is $`1/\sqrt2`$, proving

```math
\|V-M_m\otimes I_{\rm dirty}\|
\ge\frac{3/4-1/\sqrt2}{2}\approx0.0214466
\qquad(D_T(V)\le2).
```

At $`m\ge5`$, use $`P=Z_3`$. The coefficient is $`7/8`$, whose
nearest alphabet point is one. This yields the simpler stronger bound

```math
\|V-M_m\otimes I_{\rm dirty}\|\ge\frac1{16}.
```

Thus exact full-input implementations of these sources require at least
three T layers, as do approximations below the displayed tolerances.
No growing lower bound in m follows from this argument.
At $`m=3`$, the existing optimized source word has one T layer on each
side of its Clifford middle, so it implements the source exactly with
two T layers. The width threshold for excluding two layers is therefore
sharp, without asserting source-depth optimality at larger widths.

### Optimizing over the source's diagonal witnesses

The subset weights attain every multiple of $`2^{-(m-1)}`$ between
zero and one. Indeed, the first m−1 integer-scaled weights are the binary
digits $`2^{m-2},\ldots,1`$, and the final duplicate supplies the endpoint
$`2^{m-1}`$. Hence the displayed transfer entries cover the entire grid

```math
\mathcal G_m=\{-1,-1+h_m,\ldots,1\},\qquad h_m=2^{2-m}.
```

Let

```math
G=1-1/\sqrt2,\qquad t_*=(1+1/\sqrt2)/2,
\qquad x_m=h_m\mathrm{round}(t_*/h_m).
```

For $`m\ge5`$, the best gap among these diagonal witnesses is exactly

```math
\max_{x\in\mathcal G_m}\mathrm{dist}(x,\mathcal A)
=\frac G2-|x_m-t_*|.
```

To see optimality within this witness family, the interval
$`[1/\sqrt2,1]`$ has the largest alphabet gap. The available point
$`7/8`$ already achieves distance $`1/8`$. Every smaller positive
alphabet gap has half-width at most
$`(1/\sqrt2-1/2)/2\lt1/8`$, and negative gaps are symmetric.
Within the largest gap, the grid point nearest its midpoint is optimal.

Consequently the lower bound on operator distance is

```math
\frac G4-\frac12|x_m-t_*|
\ge\frac G4-\frac{h_m}{4}
\longrightarrow\frac{1-1/\sqrt2}{4}\approx0.0732233.
```

This optimizes a single diagonal transfer witness, not the true distance
to the set of all two-layer circuits.

## 4. Controlled source

For the controlled source

```math
C(M_m)=|0\rangle\langle0|\otimes I+
|1\rangle\langle1|\otimes M_m,
```

the witness
$`I_{\rm control}\otimes P_S`$ has diagonal transfer coefficient

```math
\frac12\left(1+1-2\sum_{j\in S}w_j\right)
=1-\sum_{j\in S}w_j.
```

These coefficients fill $`[0,1]`$ with spacing $`2^{1-m}`$.
For $`m\ge4`$, $`S=\{2\}`$ gives $`7/8`$, and therefore every
full-input circuit V of T-depth at most two satisfies

```math
\|V-C(M_m)\otimes I_{\rm dirty}\|\ge\frac1{16}.
```

The control is an arbitrary input qubit, not an initialized helper.
The optimized formula from Section 3 applies with spacing $`h_m/2`$
for $`m\ge4`$, with the same limiting constant. At $`m=3`$, the
$`Z_1`$ witness instead gives coefficient $`3/4`$ and the smaller
positive bound from Section 3.
At $`m=2`$, the existing one-layer loader and its actual inverse around
the controlled central Pauli give an exact controlled source with two
T layers, so this exclusion also starts at the smallest possible width.

## 5. Scope

These bounds allow arbitrary Clifford interlayers and arbitrarily many
dirty helpers. They therefore differ from the growing but restricted
Majorana-linear depth bounds in [the source-depth chapter](SOURCE_T_DEPTH.md).
Neither argument replaces the other.

An initialized-clean isometry contract is excluded. Restricting input
ancillas to zero inserts their stabilizer projector into the trace and
can sum many full-space transfer coefficients; that sum need not belong
to the alphabet. Closeness on those initialized columns also does not
imply full operator-norm closeness. The two clean flags in the complete
frame compiler cannot be silently treated as dirty inputs to obtain a
frame lower bound.

The result does not establish source depth optimality, a lower bound for
a source loader, an additive cost across source calls, or a growing
unrestricted complete-frame depth bound. A joint compiler may avoid
approximating this source as an isolated full-input primitive. The
complete-frame and state-based resource frontiers remain unchanged.

The [bounded checks](../tests/test_shallow_source_obstruction.py) audit
native two-layer words, exact rational subset grids, the optimized
coefficient gaps, native source witnesses with dirty extensions, and
the small-width two-layer constructions. Their finite matrices support
these identities; the uniform exclusion and robustness follow from the
proofs above.
