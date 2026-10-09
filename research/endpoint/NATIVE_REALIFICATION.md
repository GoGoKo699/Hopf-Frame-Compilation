# Native realification and exact tree-pattern restrictions

[Complete-frame contract](../../docs/HOPF_INTERFACE.md) · [One-clean grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit) · [Endpoint](../../docs/OPEN_PROBLEM.md)

Every native compiler for a real target has an exactly real native
implementation with one additional clean or dirty qubit, at most twice
the T-count, and the same full-isometry error. In contrast, exact
ring-valued orthogonal matrices that preserve the frame's mandatory
marker supports have a constant approximation gap. The distinction is
between the enlarged physical completion and an exactly preserved
logical sparsity pattern.

## 1. Realification preserves the full initialized-isometry error

For a complex, possibly rectangular matrix B, define

```math
\mathcal R(B)=
\begin{pmatrix}\mathrm{Re}B&-\mathrm{Im}B\\
\mathrm{Im}B&\mathrm{Re}B\end{pmatrix}.
```

This real-algebra homomorphism preserves products and adjoints. With
$`F=HS=2^{-1/2}\begin{pmatrix}1&i\\1&-i\end{pmatrix}`$,

```math
(F\otimes I)\mathcal R(B)(F^\dagger\otimes I)
=B\oplus\overline B,
\qquad \|\mathcal R(B)\|=\|B\|.
```

The identities on the two sides have the appropriate sizes for
rectangular B. Its singular values are therefore repeated twice by
realification, with no quotient by scalar phase.

**Theorem.** Suppose a finite native word V uses a clean qubits, b
arbitrary dirty qubits and t T or inverse-T gates, and satisfies

```math
\|VJ_a-J_a(W\otimes I_b)\|\le\varepsilon
```

for a real target W. There is an exactly real native word implementing
$`\mathcal R(V)`$, with the following alternatives:

| Additional register | Clean qubits | Dirty qubits | T-count | Full error |
|---|---:|---:|---:|---:|
| Arbitrary qubit | a | b+1 | at most 2t | at most epsilon |
| Initialized qubit | a+1 | b | at most 2t | at most epsilon |

Both variants include clean leakage, dirty return, arbitrary reference
correlations and the literal target phase.

*Proof of the error claim.* Put $`S_0=W\otimes I_b`$ and
$`E=VJ_a-J_aS_0`$. Since $`J_a`$ and $`S_0`$ are real,

```math
\mathcal R(V)(I_2\otimes J_a)
-(I_2\otimes J_a)(I_2\otimes S_0)=\mathcal R(E).
```

The norm is $`\|E\|`$ on the entire added-register space, including
complex superpositions. Consequently the additional qubit may be dirty
and entangled with any other input or inaccessible reference. Tensoring
with a reference preserves the operator norm. Restricting its input to
$`|0\rangle`$ instead gives the clean variant, with its return error
already included. ∎

### Exact native gate words

The extra realification register is r. The real gates H and CNOT act
unchanged on their original wires. Use the literal real-rotation convention

```math
R_y(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix}.
```

Then

```math
\begin{aligned}
\mathcal R(S_x)
&=C_x[R_y(\pi/2)_r]
=\mathrm{CNOT}_{x,r}\mathrm{CZ}_{x,r},\\
\mathcal R(T_x)
&=C_x[R_y(\pi/4)_r]
=C_x[H_r]\mathrm{CZ}_{x,r}.
\end{aligned}
```

These follow from $`R_y(\pi/2)=XZ`$ and
$`R_y(\pi/4)=HZ`$. Put $`A=SHTHS^\dagger`$. Direct multiplication gives
$`AZA^\dagger=H`$, hence

```math
C_x[H_r]=(I_x\otimes A_r)\mathrm{CZ}_{x,r}
                 (I_x\otimes A_r^\dagger).
```

This word contains exactly one T and one inverse T; CZ is
H–CNOT–H. Replace inverse gates by the actual inverse words. Products
then implement $`\mathcal R(V)`$ exactly, using only H, S, CNOT, T
and their native inverses. SWAPs can place r anywhere. No arbitrary
rotation primitive or supplied state enters this construction.

**Endpoint specialization.** The
[one-clean grouped compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit)
leaves the second supplied clean flag available as r. At
$`L=N`$ it therefore gives an exactly real native completion with

```math
a=2,\qquad b=N+n+7,\qquad
T=O\bigl(N[1+\log_2^*(n+2)]\bigr),\qquad
\|\mathcal R(V)J_2-J_2(W\otimes I_b)\|\le2^{-N},
```

after the fixed register reordering. Thus requiring an exactly real
enlarged native completion does not change the established asymptotic
upper bound. Realification does not by itself reduce that bound to
$`O(N)`$.

## 2. Exact ring-valued marker supports have a constant gap

Let $`R=\mathbb Z[1/\sqrt2]`$. Every exact native Clifford+T matrix has
entries in $`R[i]`$, by inspection of the generators and closure under
matrix multiplication. An exactly returned real logical unitary
therefore has entries in R, even with clean or arbitrary dirty work.

**Circle lemma.** The only $`x,y\in R`$ satisfying $`x^2+y^2=1`$ are

```math
(\pm1,0),\quad(0,\pm1),\quad
(\pm2^{-1/2},\pm2^{-1/2}),
```

with all independent signs in the last pair.

*Proof.* Write $`x=(a+b\sqrt2)/2^k`$ and
$`y=(c+d\sqrt2)/2^k`$, with integer coefficients and $`k\ge0`$.
Rational and irrational parts give

```math
a^2+c^2+2b^2+2d^2=4^k,\qquad ab+cd=0.
```

For $`k\ge1`$, a and c have the same parity. If both were odd,
the first equation modulo four would force b and d to have opposite
parity, contradicting the second equation modulo two. Thus a and c
are even. If $`k\ge2`$, set $`a=2a'`$, $`c=2c'`$. Dividing the
first equation by two shows that b and d have the same parity. If both
were odd, reduction modulo four would force $`a',c'`$ to have
opposite parity, again contradicting $`a'b+c'd=0`$ modulo two.
All four coefficients are even, so divide by two and descend in k.
At $`k=0`$ only the axis points occur; $`k=1`$ adds exactly the four
diagonal points. ∎

For each tree node $`(d,p)`$, let

```math
I_{d,p}=\{p2^{n-d},\ldots,(p+1)2^{n-d}-1\},\qquad
z_{d,p}=(2p+1)2^{n-d-1}.
```

The prescribed frame's marker column $`z_{d,p}`$ is supported in
$`I_{d,p}`$. No support restriction is imposed on column zero.

**Exact-pattern theorem.** Every real orthogonal Q with these marker
supports is a Hopf frame times a possible sign on its root column.
If Q has entries in R, every frame angle is a multiple of
$`\pi/4`$ modulo $`2\pi`$. Consequently this exact-pattern class
has at most $`2\cdot8^{N-1}`$ ring-valued matrices.

*Proof.* A bottom marker is a unit vector on coordinates
$`2p,2p+1`$. Choose the real plane rotation sending the second basis
vector to that marker. Multiplying Q by its transpose sends the marker
to $`|2p+1\rangle`$; orthogonality makes every other column zero in
that odd row. Remove all odd marker columns and rows. The remaining
matrix is orthogonal with the smaller tree's marker supports. Repeating
this operation leaves a root coordinate equal to +1 or −1 and proves
the first assertion.

If Q is ring-valued, the circle lemma makes each bottom stripping
rotation ring-valued with angle a multiple of $`\pi/4`$. The residual
matrix remains ring-valued, so induction proves the second assertion.
For an original Hopf frame the bottom marker is exactly

```math
-\sin\theta_{n-1,p}|2p\rangle+
\cos\theta_{n-1,p}|2p+1\rangle.
```

This also restricts each of its prescribed angles. No ancestor amplitude is divided out, and
singular angles and all signs are included. ∎

**Uniform approximation gap.** Set only
$`\theta_{n-1,0}=\pi/8`$ nonzero. Its marker column 1 has distance
at least

```math
2\sin(\pi/16)>0.39
```

from every ring-valued unit vector supported on coordinates 0 and 1.
Hence the operator distance from every ring-valued Q in the stated
exact-pattern class has the same lower bound, independently of n.
Additional denominator bits or exactly returned work cannot remove it.

The support hypothesis is essential: small off-pattern entries and
approximate work return lie outside this theorem. In particular, the
realified physical completion of Section 1 is not required to preserve
the logical marker pattern.

## 3. Verification

The [native-representation tests](../../tests/test_native_representation_limits.py) verify the literal gate words
in exact $`\mathbb Q(\sqrt2,i)`$ arithmetic, check the circle equations
over bounded integer ranges, and exercise the rectangular norm and
full-isometry identities. Small frame checks confirm the marker-column
gap. The all-size conclusions follow from the proofs above, independently
of these finite checks.
