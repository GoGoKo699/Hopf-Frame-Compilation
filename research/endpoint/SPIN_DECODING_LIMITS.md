# Access and decoding limits for doubled Spin lifts

[Source representations](SOURCE_REUSE_LIMITS.md#10-higher-grade-spin-identities-on-arbitrary-dirty-work) · [Complete-frame contract](../../docs/HOPF_INTERFACE.md) · [Endpoint](../../docs/OPEN_PROBLEM.md)

A Spin lift represents every complete real Hopf frame through Majorana
conjugation. Its doubled physical representation does not provide a
normalization-one coherent decoder on arbitrary dirty inputs. For
$`n\ge4`$, a single explicit frame rules out every target-dependent
inverse encoder around one doubled lift, regardless of encoding cost.
A second frame rules out unrelated target-dependent Clifford wrappers.
Both bounds allow all available physical wires to participate in the
encoding and include the prescribed exponential accuracy.

Throughout, $`N=2^n`$, $`m=N/2`$ and $`\eta=2^{-N}`$. The physical
width is $`q=n+b+2`$, with two initialized flags and b arbitrary dirty
qubits. A doubled lift requires $`q\ge N`$; write
$`s=q-N`$ for its spectator qubits. The bounds below hold for every
such b, in particular the endpoint allocation $`b=N+n+7`$.

## 1. Spin access and normalization

On m qubits, use the Hermitian Jordan–Wigner Majoranas

```math
\Gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\Gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j.
```

They square to identity and anticommute. A Givens rotation sending
$`e_a`$ to $`\cos\theta\,e_a+\sin\theta\,e_b`$ has the lift

```math
S_{ab}(\theta)=\exp[(\theta/2)\Gamma_b\Gamma_a].
```

Its conjugation sends $`\Gamma_a`$ to
$`\cos\theta\,\Gamma_a+\sin\theta\,\Gamma_b`$. Use the tree edges

```math
a=p2^{n-d},\qquad b=a+2^{n-d-1},\qquad
0\le d\lt n,\quad 0\le p\lt2^d,
```

in the frame's prescribed order, increasing in depth. Then

```math
S(W)\Gamma_jS(W)^\dagger=\sum_iW_{ij}\Gamma_i,
\qquad \rho(W)=S(W)\otimes\overline{S(W)}.
```

Normalized operator vectorization obeys

```math
\rho(W)|A\rangle\!\rangle
=|S(W)AS(W)^\dagger\rangle\!\rangle.
```

The doubled lift occupies N physical qubits. Preparing the
$`|\Gamma_j\rangle\!\rangle`$ vectors is a separate operation;
arbitrary dirty input does not supply them. Likewise, the N−1
independently specified rotations in S retain their synthesis charges.
Given a native word for S, its conjugate is native by exchanging S
with its inverse and T with its inverse, while H and CNOT remain real;
both copies are charged.

### Scalar anticommutators do not supply a normalized block

The exact scalar identity is

```math
\frac{\Gamma_iS\Gamma_jS^\dagger+
S\Gamma_jS^\dagger\Gamma_i}{2}=W_{ij}I.
```

With independent input and output address registers, an LCU flag can
select the two product orders. A uniform output-address preparation
and uniform input-address projection then give $`W/N`$, with an
additional n-qubit initialized address register. Neither this
initialization nor its amplification is part of a two-clean decoder.

The one-address alternative does not produce the same pair. For
$`Q=\sum_i|i\rangle\langle i|\otimes\Gamma_i`$ and an address
mixer A,

```math
[Q(A\otimes I)SQS^\dagger]_{ij}
=A_{ij}\Gamma_iS\Gamma_jS^\dagger.
```

Reversing the SELECT order instead gives
$`A_{ij}(S\Gamma_iS^\dagger)\Gamma_j`$, with the wrong address
index in the conjugated factor. An address partial transpose is not
a free quantum operation.

The normalization issue is also visible algebraically. Define block
matrices $`F_{ij}=\Gamma_i(S\Gamma_jS^\dagger)`$,
$`G_{ij}=(S\Gamma_jS^\dagger)\Gamma_i`$, and
$`(Q_0)_{ij}=\Gamma_i\Gamma_j`$. Then

```math
\begin{aligned}
Q_0&=NVV^\dagger,
&V|\phi\rangle&=\frac1{\sqrt N}\sum_i|i\rangle\Gamma_i|\phi\rangle,\\
F&=Q_0(W\otimes I),
&G&=(2I-Q_0)(W\otimes I).
\end{aligned}
```

Thus $`(F+G)/2=W\otimes I`$, but $`\|F\|=N`$ and
$`\|G\|=N-2`$ for $`N\ge4`$. Their scalar cancellation alone is
not a contraction or a normalized physical block.

## 2. A single target defeats target-dependent inverse encoders

**Theorem.** For $`n\ge4`$, there is an admitted frame W such that
every isometry E from the n+b logical/dirty input qubits into the q
physical qubits satisfies

```math
\| (\rho(W)\otimes I_s)E-E(W\otimes I_b)\|
\ge\Delta,\qquad
\Delta=2\sin(\pi/3^m)>2^{-N}.
```

E may depend on W and have arbitrary gate count. In particular this
excludes the architecture $`U^\dagger(\rho(W)\otimes I_s)UJ_2`$
for any target-dependent physical encoder U, by taking $`E=UJ_2`$.

*Proof.* Set all upper-layer angles to zero and the deepest angles to

```math
\theta_j=\frac{2\pi3^j}{3^m},\qquad 0\le j\lt m.
```

W is a direct sum of m real rotations with distinct eigenvalues
$`e^{\pm i\theta_j}`$. The canonical Majorana pairing makes each
Spin factor $`\exp(-i\theta_jZ_j/2)`$. Its doubled pair has phases
$`-\theta_j,0,0,+\theta_j`$. Hence a signed ternary word
$`e\in\{-1,0,1\}^m`$ gives phase and multiplicity

```math
\phi_e=\sum_je_j\theta_j,\qquad
\mathrm{mult}(e)=2^{m-|\mathrm{supp}e|}.
```

The balanced ternary sums give each residue modulo $`3^m`$ exactly
once. Therefore every desired eigenvalue $`e^{\pm i\theta_j}`$
has multiplicity $`2^{m-1}`$ in the doubled lift. With all available
spectators its eigenspace has dimension

```math
r=2^{m-1+s}=2^{b+n+1-m}\lt2^b
\qquad(n\ge4).
```

Fix its logical eigenvector u. On the $`2^b`$-dimensional input
space $`u\otimes\mathbb C^{2^b}`$, the composition of E with the
lifted eigenprojector has a nontrivial kernel. Choose a unit vector v
there. Every other lifted eigenvalue is at distance at least
$`\Delta`$ from the target eigenvalue, so its intertwining error
on v is at least $`\Delta`$. Finally,

```math
\Delta\ge\frac4{3^m}>\frac1{4^m}=2^{-N},
```

using $`\sin x\ge2x/\pi`$ on $`[0,\pi/2]`$. ∎

**Small-size and approximation scope.** At $`n=3`$ the ratio
$`r/2^b=2^{n+1-m}`$ equals one, so this theorem gives no obstruction.
The fixed-family theorem in Section 4 still applies at that size.
If a native implementation satisfies
$`\|\widetilde\rho-\rho\|\le\xi`$ on the whole physical space,
the lower bound becomes $`\Delta-\xi`$. An implementation error
larger than the witness gap cannot be ignored.

## 3. Unrelated Clifford wrappers

The inverse-encoder hypothesis matters: unrelated target-dependent
general unitaries could include W in the decoder itself. For Clifford
wrappers a separate scaled target gives an obstruction.

**Theorem.** For $`n\ge4`$, set all upper angles to zero and

```math
\theta_j=\frac{3^j}{4\cdot3^m},\qquad 0\le j\lt m.
```

For arbitrary Clifford circuits C and D on all q physical qubits,
including target-dependent choices, the full initialized-isometry
error obeys

```math
\|D(\rho(W)\otimes I_s)CJ_2-J_2(W\otimes I_b)\|
\ge g,\qquad
g=\sin\!\left(\frac1{8\cdot3^m}\right)>2^{-N}.
```

*Proof.* Suppress spectator identities and put
$`\delta=\|\rho-I\|\lt1/8`$ and
$`\omega=\|W-I\|\lt1/12`$, using the geometric sum of the
angles. Assume the displayed error $`\varepsilon\lt g`$ and
write $`E=CJ_2`$, $`F=D^\dagger J_2`$. Their ranges are equal-rank
stabilizer codes, and

```math
\|E-F(W\otimes I_b)\|\le\delta+\varepsilon\lt1/\sqrt2.
```

Distinct equal-rank stabilizer codes have a unit vector in one code
at distance at least $`1/\sqrt2`$ from the other. To see this,
expand their rank-r projectors P,Q over their signed stabilizer groups.
An opposite-sign common Pauli makes the codes orthogonal. Otherwise,
if the groups differ, their intersection is at most half of either
group and $`\mathrm{tr}(PQ)\le r/2`$. Averaging over an
orthonormal basis of the first code proves the assertion. The displayed
norm bound therefore forces the two code ranges to coincide.

Consequently $`F=EK`$, where K is the induced logical Clifford on
the n+b input qubits. Indeed DC preserves the standard initialized
code, and its induced action and inverse preserve logical Paulis.
Moreover,

```math
\|K-I\|\le\delta+\omega+\varepsilon\lt1/\sqrt2.
```

For every Pauli P,
$`\|KPK^\dagger-P\|\le2\|K-I\|\lt\sqrt2`$.
Distinct signed Paulis are at distance at least $`\sqrt2`$, so K
fixes every Pauli under conjugation and is scalar. No restriction on
its scalar phase is needed.

The scaled balanced ternary grid has minimum eigenvalue separation
$`2g`$. Its largest eigenspace, at zero phase, has multiplicity
$`2^m`$ before spectators and

```math
2^{m+s}=2^{b+n+2-m}\lt2^b\qquad(n\ge4)
```

afterwards. For any scalar-shifted desired eigenvalue, at most one
lifted eigenvalue lies strictly within distance g. Choose a dirty
direction whose encoded image has zero projection onto that one
eigenspace. Every remaining component then has error at least g,
contradicting the assumption. Finally,

```math
g\ge\frac5{48\cdot3^m}>4^{-m}=2^{-N}
\qquad(m\ge8),
```

because $`\sin x\ge(5/6)x`$ for $`0\le x\le1`$ and
$`(5/48)(4/3)^8>1`$. ∎

For an implemented lift with full-operator error at most xi, the lower
bound is $`g-\xi`$. General unrelated non-Clifford wrappers need
not preserve stabilizer codes and are outside this theorem, as are
repeated lift queries or interleaved non-Clifford access. The kernel
witnesses are already pure dirty inputs; allowing inaccessible
references cannot remove the bounds.

## 4. A fixed encoder on the even-sign family

This separate family theorem covers $`n=3`$ and permits an arbitrary
fixed encoding. Setting every Hopf angle to zero or pi gives exactly

```math
G=\{D_t=\mathrm{diag}(t_0,\ldots,t_{N-1}):
t_i\in\{-1,1\},\ \prod_it_i=1\}.
```

Each pi edge negates its two endpoints, and the tree incidence vectors
form a basis of the even-parity subspace over $`\mathbb F_2`$.
Thus every displayed element belongs to the prescribed Hopf family.

The $`2^N`$ Clifford monomials $`\Gamma_A`$ form an orthonormal
operator basis and have characters

```math
\rho(D_t)|\Gamma_A\rangle\!\rangle
=\chi_A(t)|\Gamma_A\rangle\!\rangle,
\qquad \chi_A(t)=\prod_{i\in A}t_i.
```

Two subsets give the same character precisely when they are equal or
complements. Every character therefore has multiplicity two, or
$`2^{s+1}`$ with spectators.

Suppose one fixed isometry E has intertwining error at most epsilon
for every $`D_t\in G`$. Average its errors against the coordinate
character $`\chi_{\{i\}}`$ and denote the physical character
projector by $`P_i`$. For $`n\ge3`$,

```math
\mathrm{rank}P_i=2^{s+1}=2^{b+n+3-N}\lt2^b.
```

There is a unit v in $`|i\rangle\otimes\mathbb C^{2^b}`$ with
$`P_iEv=0`$. The averaged error gives
$`\|P_iEv-Ev\|\le\varepsilon`$, whereas the left side is one.
Hence $`\varepsilon\ge1`$. If a fixed decoder is unrelated to
the encoder, the identity-target error first compares their initialized
images; the same argument then gives $`\varepsilon\ge1/2`$.
These are fixed-family statements; Section 2 supplies the stronger
single-target result for inverse encoders at $`n\ge4`$.

## 5. Verification

The [native-representation tests](../../tests/test_native_representation_limits.py) count tree-incidence rank, complement
character classes and signed-ternary eigenvalue multiplicities with
exact integer arithmetic at $`n=3,4`$. They also check both explicit
gap comparisons and the exceptional equality at $`n=3`$. The general
rank and gap proofs above do not require constructing exponentially
large physical matrices.
