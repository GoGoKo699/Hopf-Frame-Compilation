# Precision, T count and native T-depth for complete frames

[General resource frontier](../../docs/OPEN_PROBLEM.md) ·
[All-precision upper bound R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Four-band rank–precision bound](FOUR_LAYER_LINEAR_PRECISION.md)

<a id="precision-count-depth-tradeoff"></a>

## 1. The unrestricted precision tradeoff, R103

For every $`n\ge1`$, choose the complete real Hopf frame
$`W_n=I_{2^{n-1}}\otimes R_y(\theta)`$ with
$`R_y(\theta)=e^{-i\theta Y}`$ and $`\theta=\tfrac12\arccos(1/3)`$.
This sets every deepest-layer angle to the same $`\theta`$ and all
earlier angles to zero. It is also a phase-dressed complex magnitude
frame with all leaf phases zero.

Let V be any native Clifford+T word with k initialized zeros, b
arbitrary dirty qubits, T-count t and native T-depth d. On all n input
qubits and the full dirty space, assume

```math
\|VJ_k-J_k(W_n\otimes I_b)\|\le2^{-L},\qquad L\ge6.
```

An optimized common phase is allowed because the proof uses the induced
channel. The circuit may move data and work during its execution. Work
return is included in the displayed error. Then

```math
d\,t\ \ge\ \frac{(L-\log_2 25)^2}{(k+1)\ln 4}.
```

The inequality permits arbitrary Clifford interlayers, dirty width,
Clifford count and table-specific compilation. It charges physical T and
T-dagger targets in each layer. It applies to the actual native inverse
with the same error, as the ideal inverse has the same Pauli coefficient.

For one clean input and $`L\ge N=2^n`$, the simultaneous R96 count budget
is $`O(L)`$ at its stated real and complex dirty thresholds. The inequality
requires $`\Omega(L)`$ T-depth in that budget. R96 attains $`O(L)`$ depth
with $`O(L)`$ T gates and $`O(NL)`$ Cliffords, so the precision-dominated
worst-case depth is $`\Theta(L)`$ at those allocations.

At sufficient square-root dirty width and $`6\le L\le N`$, the same
one-clean count budget implies $`d=\Omega(L^{3/2}/\sqrt N)`$.
Consequently a bound of the [R102 form](COMPACT_PROGRAM_PRECISION.md),
$`O((1+L/n)\log_2(n+2))`$ at that width, can only extend to
$`L=O(N\log_2^2(n+2)/n^2)`$ in the range $`n\le L\le N`$.
Fixed-accuracy depth remains open.

## 2. Valuations of minors of a physical T band

Put $`q=n+b+k`$, $`D=4^q`$, and $`\pi=\sqrt2`$. Use the standard
Hermitian Pauli basis on q qubits, orthonormal in normalized
Hilbert–Schmidt inner product. A transfer matrix entry of V is
$`2^{-q}\mathrm{Tr}(P VQ V^\dagger)`$. Clifford transfer matrices are
signed permutations. All Clifford+T transfer entries belong to
$`\mathbb Z[\sqrt2,1/2]`$, since a single T rotation splits an
anticommuting Pauli into two terms with coefficients of modulus
$`1/\sqrt2`$ and fixes commuting Paulis.

Use the valuation $`v_\pi`$ with $`v_\pi(\pi)=1`$. The localization of
$`\mathbb Z[\sqrt2]`$ at the ideal $`(\pi)`$ is a discrete valuation ring.
Indeed, for $`u,v\in\mathbb Z`$ the two terms in $`u+v\pi`$ have distinct
valuation parity whenever both are nonzero, and factoring out the
smaller power of $`\pi`$ leaves an element nonzero modulo $`\pi`$.
Diagonal elimination over this local ring gives Smith form, with powers
of $`\pi`$ ordered by increasing exponent. It applies to transfer
matrices by first multiplying by a common power of $`\pi`$.

For a matrix B, let $`\tau_j(B)`$ denote the minimum valuation of a
$`j\times j`$ minor, with $`\tau_0(B)=0`$ and zero minors given infinite
valuation. This is the sum of the j smallest Smith exponents. The
Cauchy–Binet expansion gives, also for rectangular matrices,

```math
\tau_j(B_1 B_2)\ge\tau_j(B_1)+\tau_j(B_2).
```

An integral boundary matrix has nonnegative minor valuations.

For a one-wire T or T-dagger, multiply its transfer matrix by $`\pi`$.
Its I and Z entries are $`\pi`$, and its X,Y block has unit entries
and determinant 2. Its Smith exponents are $`0,1,1,2`$. In a physical
band with h distinct T targets, tensor the one-wire Smith forms and the
spectator identities. The exponents of the scaled band are the multiset
of a binomial variable $`X\sim\mathrm{Bin}(2h,1/2)`$, with total
multiplicity D. For the unscaled band, subtract h from each exponent.

Take the smallest j exponents, breaking ties arbitrarily, and let E be
that event in the uniform multiset. For $`0\lt j\lt D`$ its probability
is $`j/D`$. For any $`\lambda\gt0`$, conditional Jensen and the binomial
exponential moment give

```math
\frac{j}{D}\exp(\lambda(h-\mathbb E[X\mid E]))
\le \mathbb E e^{\lambda(h-X)}
=\cosh(\lambda/2)^{2h}
\le e^{h\lambda^2/4}.
```

Minimizing in $`\lambda`$ shows that an unscaled band B satisfies

```math
\tau_j(B)\ge-j\sqrt{h\ln(D/j)}.
```

The j zero and empty-band cases contribute zero. For d bands with
widths $`h_i`$, $`\sum_i h_i=t`$ and $`\sum_i\sqrt{h_i}\le\sqrt{dt}`$.
This proves the same minor bound for a full word with penalty
$`j\sqrt{dt\ln(D/j)}`$ in place of the single-band penalty. Integral
boundary restrictions and injections preserve this bound.

## 3. A large scalar block of the initialized channel

Initialize the k clean inputs with J and trace those output wires. Write
$`\mathcal E(X)=\mathrm{Tr}_{\rm clean}(VJXJ^\dagger V^\dagger)`$.
Let A be the block of its Pauli transfer matrix whose inputs and outputs
are Z on the chosen logical target and arbitrary Pauli labels on every
spectator, both logical and dirty. It is an $`r\times r`$ matrix, where

```math
s=2^{n+b-1},\qquad r=s^2=D/4^{k+1}.
```

The ideal channel block is $`cI_r`$ with $`c=1/3`$. This block uses the
complete work-return contract on every spectator operator.

The input and output maps that select A in the full physical transfer
matrix are integral. In normalized transfer coordinates the input
projector expands as a sum of all k-bit Z strings with coefficients
$`2^{-k}`$. That factor cancels the change from physical trace
normalization $`2^{-q}`$ to input-space normalization
$`2^{-(q-k)}`$. Thus one injects the sum of labels
$`P\otimes Z^u`$ at each input Pauli P, then keeps the output label
$`Q\otimes I`$ for each desired Q. In particular
$`A\in\mathbb Z[\sqrt2,1/2]^{r\times r}`$, and Section 2 gives

```math
\tau_j(A)\ge-j\sqrt{dt\ln(D/j)}\qquad (1\le j\le r).
```

### Frobenius estimates on the block

The Frobenius norm below is the Euclidean matrix norm of the indicated
transfer block, namely the square root of the sum of squared entries.
Suppose a channel on logical target and spectators has Kraus matrices
$`K_a=\sum_P P\otimes B_{a,P}`$, with P ranging over the four target
Paulis. Its Z-to-Z block is the linear map on spectator matrices

```math
X\ \,\longmapsto\ \,\sum_{a,P}\chi_Z(P)B_{a,P}X B_{a,P}^\dagger,
```

where $`\chi_Z(P)`$ is the commutation sign. This follows by tracing the
target Pauli products, which kills cross terms on distinct P. The map
$`X\mapsto BXB^\dagger`$ has transfer Frobenius norm $`\|B\|_{\mathrm F}^2`$
by vectorization. If $`\sum_a K_a^\dagger K_a\preceq\delta I`$, partial
trace over the target gives $`\sum_{a,P} B_{a,P}^\dagger B_{a,P}\preceq\delta I_s`$.
Thus the block's Frobenius norm is at most $`\delta s=\delta\sqrt r`$.

Put $`\epsilon=2^{-L}`$. The principal clean-output Kraus matrix
$`K_0=J^\dagger VJ`$ differs from the ideal unitary by at most
$`\epsilon`$ (in a common phase convention). Its channel differs from
the ideal by at most $`2\epsilon`$ in Hilbert–Schmidt induced norm, by
telescoping the two Kraus factors. Compression to the selected block
keeps this bound. All other clean-output Kraus matrices come from the
leaked part of VJ, so they satisfy
$`\sum_{a\ne0}K_a^\dagger K_a\preceq\epsilon^2 I`$.
The preceding block estimate therefore yields

```math
\|A-cI_r\|_{\mathrm F}\le(2\epsilon+\epsilon^2)\sqrt r.
```

For the real field conjugation $`\sigma(\sqrt2)=-\sqrt2`$, extend
$`\sigma`$ by $`e^{i\pi/4}\mapsto e^{3i\pi/4}`$. It sends the native
generators to unitary matrices and commutes with conjugate transpose.
The conjugate of the transfer block is therefore a signed coordinate
change of the corresponding block of a trace-preserving induced
channel. The signs only arise from Pauli Y. The same Frobenius estimate
with $`\delta=1`$ gives $`\|\sigma(A)\|_{\mathrm F}\le\sqrt r`$.
Neither estimate charges a factor growing with the clean dimension.

## 4. A nonzero determinant at the target precision

Set $`B=3A-I_r`$. The determinant is nonzero because B reduces to
$`-I_r`$ modulo 3 in $`\mathbb Z[\sqrt2,1/2]`$, where 2 is invertible.
For any nonzero $`z\in\mathbb Z[\sqrt2,1/2]`$,

```math
2^{v_\pi(z)}\le |z\,\sigma(z)|.
```

To see this, clear dyadic denominators by multiplying by a power of
$`\pi`$. The algebraic norm of the remaining nonzero integral element
after removing all factors of $`\pi`$ is a nonzero integer. Each removed
factor contributes absolute norm 2, and the cleared denominators cancel.

The product of a square matrix's singular values is at most the
corresponding power of its root-mean-square singular value. Apply this
to B and its conjugate, using the two Frobenius estimates and c=1/3.

```math
\begin{aligned}
|(\det B)\sigma(\det B)|
 &\le [12(2+\epsilon)\epsilon]^r\le(25\epsilon)^r,\\
v_\pi(\det B)&\le-r(L-\log_2 25).
\end{aligned}
```

On the other hand, expand $`\det(3A-I_r)`$ in principal minors of A,
including the j zero term. All coefficients are integral, so its
valuation is at least the minimum of $`\tau_j(A)`$ for $`0\le j\le r`$.
The function $`j\sqrt{\ln(D/j)}`$ is increasing up to r, since
$`r\le D/4`$. Section 3 yields

```math
v_\pi(\det B)\ge-r\sqrt{dt\ln(D/r)}.
```

These two bounds and $`\ln(D/r)=(k+1)\ln 4`$ prove the theorem.
The valuation argument only uses the actual circuit and its complete
initialized-isometry error. In particular it allows leakage and
approximate return at that error.
