# An exact T-count barrier for the separate Walsh dual readout

[Paid Majorana construction](MAJORANA_WALSH_QUERY.md) ·
[Complete probe interface](SYMPLECTIC_MODULE_QUERY_BARRIER.md#6-universal-pauli-probes-contain-a-complete-dirty-indicator) ·
[Native query R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Conditional frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md)

The separate dual-frame readout in the Majorana proposal cannot have
linear T-count, even after unrestricted native resynthesis. For
$`m=2^p`$ modes it requires

```math
T(Q)\ge2L_p,\qquad
L_p=p\binom{p-1}{\lfloor(p-1)/2\rfloor}
       =\Theta(m\sqrt p).                                      \qquad(1)
```

This is an **exact standalone-source** bound. It permits arbitrary
Clifford interlayers, arbitrary returned dirty helpers and any finite
initialized-and-returned clean allocation. It is not a Gaussian-circuit
restriction, a lower bound for the fully fused indicator, or an
approximation lower bound. The latter distinctions matter: R91 already
implements that indicator with O(m) T gates, and the coefficient used
in this proof is exponentially small in $`L_p`$.

The [approximate continuation](APPROXIMATE_QUERY_FRAME_TRANSFER.md)
shows that fixed frame precision needs only constant readout accuracy
as the address grows. Its robust Gaussian-conjugator obstruction has
an explicitly narrower native-decomposition hypothesis than this exact
standalone theorem.

## 1. The complete readout and its phases

Take $`p\ge2`$, $`m=2^p`$ and $`w=m/2`$ arbitrary carrier qubits.
Let $`\gamma_j`$ be the m Jordan--Wigner Majoranas and W the actual
native Walsh unitary from the linked construction, so

```math
W\gamma_jW^\dagger=A_j
 =\sum_u H_{uj}\gamma_u,\qquad
H_{uj}=m^{-1/2}(-1)^{u\cdot j},\qquad H=H^{\mathsf T}.
```

The Majorana binary Pauli labels form a basis: their symplectic Gram
matrix is $`I+\mathbf1\mathbf1^{\mathsf T}`$, which is invertible
for even m. Choose positive Hermitian physical Paulis $`\Delta_j`$
whose conjugation negates only $`\gamma_j`$. Write
$`\delta_j=W\Delta_jW^\dagger`$. These dual operators are Hermitian
involutions and anticommute pairwise. For an arbitrary m-bit register Y,
fix the matrix-product order with index zero leftmost:

```math
\begin{aligned}
C_Y&=\Delta_0^{Y_0}\Delta_1^{Y_1}\cdots\Delta_{m-1}^{Y_{m-1}},\\
Q_Y&=WC_YW^\dagger
    =\delta_0^{Y_0}\delta_1^{Y_1}\cdots\delta_{m-1}^{Y_{m-1}},\\
Q&=\sum_Y |Y\rangle\langle Y|\otimes Q_Y.
\end{aligned}                                                   \qquad(2)
```

C is a Clifford word of controlled tensor Paulis. Neither C_Y nor Q_Y
is assumed Hermitian; every inverse below is the actual adjoint in the
specified order. The target Q acts on every Y/carrier input and reference.
Its known separated implementation uses $`mp`$ T gates. The lower bound
(1) concerns every exact native implementation of this same operator,
including implementations whose intermediate gates leave the Majorana span.

## 2. A readout-only indicator interface

It is unnecessary to supply a separate shallow flat-vector source if a
complete readout is available. Let $`X_0`$ denote X on output wire
$`Y_0`$. Because index zero is leftmost in (2),

```math
Q_YQ_{Y\oplus e_0}^\dagger=\delta_0,
\qquad QX_0Q^\dagger X_0=I_Y\otimes\delta_0.          \qquad(3)
```

The source may even be supplied with a completed Y-dependent scalar:
$`\widetilde Q|Y,a\rangle=e^{i\theta(Y)}|Y\rangle Q_Y|a\rangle`$.
Its two-call extraction in (3) is then $`r(Y)\delta_0`$, where
$`r(Y)=e^{i(\theta(Y)-\theta(Y\oplus e_0))}`$. This scalar will cancel.
This transfer allows a diagonal Y-phase gauge; it does not assume that
the transfer works for an arbitrary right register unitary.

Use the coherent Clifford masks P_x of the linked construction, with
$`P_xA_0P_x^\dagger=A_x`$. They are even for $`p\ge2`$. If
$`\mathsf P=Z_0\cdots Z_{w-1}`$ is carrier parity, then
$`\Delta_j=\sigma_j i\mathsf P\gamma_j`$ for fixed signs
$`\sigma_j\in\{1,-1\}`$. W and P_x commute with parity, so
$`P_x\delta_0P_x^\dagger=\pm\delta_x`$. Define the actual word

```math
D_x=P_x\widetilde QX_0\widetilde Q^\dagger X_0P_x^\dagger.
```

It preserves Y and equals $`r(Y)`$ times a signed $`\delta_x`$.
Pairwise anticommutation therefore gives the literal identity

```math
\left(\prod_j Z_{Y_j}\right)
D_x\widetilde QD_x^\dagger\widetilde Q^\dagger
   =(-1)^{Y_x}I_{\mathrm{carrier}}.                             \qquad(4)
```

The actual $`D_x^\dagger`$ cancels both the extracted scalar and every
native phase; replacing it by a presumed Hermitian word is unnecessary.
Hadamards on all Y wires before and after (4) give
$`|x,Y,a\rangle\mapsto|x,Y\oplus e_x,a\rangle`$, with literal phase
one. Coherent addresses, arbitrary output/carrier inputs and references
are included. No carrier initialization or second control bank is used.

Suppose common forward/actual-inverse source bounds are
$`T_Q,G_Q,d_Q`$, with $`b_Q`$ returned dirty helpers. Four P_x
orientations cost at most $`20pw`$ elementary Cliffords, four X words
cost 16, and the parity correction and Hadamard endpoints cost 4m.
Thus the completely compiled wrapper has

```math
\begin{aligned}
T_{\mathrm{ind}}&\le6T_Q,&D_T(\mathrm{ind})&\le6d_Q,\\
G_{\mathrm{ind}}&\le6G_Q+20pw+4m+16,&
b_{\mathrm{ind}}&\le w+b_Q.
\end{aligned}                                                   \qquad(5)
```

The m Y wires are the indicator's outputs, not borrowed workspace.
The six source calls are sequential and reuse the returned helper pool.
The wrapper adds no clean work and inherits any clean allocation of Q;
the R98 consequence below uses a zero-clean source.
In the positive-S alphabet, common bounds must cover the actual inverse,
whose Clifford count need not equal the forward count. Consequently a
zero-clean readout with linear T/width, quadratic Clifford count and depth
d(p) would supply the complete indicator interface used in R98. The exact
linear-T premise is ruled out for this readout family by (1).

## 3. Channel denominators, including clean return

The standard channel representation and smallest-denominator exponent
are used in [Gosset, Kliuchnikov, Mosca and Russo, *An algorithm for the
T-count*, Sections 2.3 and 4, Definition 1 and equation (4.2)](https://arxiv.org/html/1308.4134v1).
The elementary bound needed here has the following self-contained proof.
For n physical qubits and Hermitian Paulis P,R, put

```math
\widehat V_{R,P}=2^{-n}\mathrm{Tr}(RVPV^\dagger),\qquad
\mathrm{sde}(v)=\min\{t\ge0:(\sqrt2)^t v\in\mathbb Z[\sqrt2]\}.
```

Define the exponent of zero to be zero. Clifford channels are signed
permutation matrices. A physical T or T-dagger leaves commuting Pauli
coordinates fixed and mixes each relevant pair with coefficients
$`\pm1/\sqrt2`$. Starting from the identity channel, a word with t
T/T-dagger gates therefore has every entry in
$`(\sqrt2)^{-t}\mathbb Z[\sqrt2]`$. Finite sums cannot increase that
denominator bound. In particular,

```math
T(V)\ge\mathrm{sde}(\widehat V_{R,P})
\quad\hbox{for every channel entry}.                            \qquad(6)
```

No restriction on Clifford mixing, circuit width or T scheduling enters
this count argument. Tensoring the target and witness Paulis with identity
on returned dirty helpers leaves the coefficient unchanged.

Initialized clean helpers do not weaken (6) either. Let V implement a
logical unitary U exactly while returning c clean qubits to zero, and put
$`\Pi=|0^c\rangle\langle0^c|`$. Then
$`V(P\otimes\Pi)V^\dagger=UPU^\dagger\otimes\Pi`$. Since
$`\Pi=2^{-c}\sum_{S\subseteq[c]}Z_S`$,

```math
\widehat U_{R,P}
 =\sum_{S\subseteq[c]}
   \widehat V_{R\otimes I_c,\,P\otimes Z_S}.                   \qquad(7)
```

This is a sum, not an average: the full-space normalization contributes
$`2^c`$, canceling the projector's $`2^{-c}`$, and
$`\mathrm{Tr}\Pi=1`$. Every summand has exponent at most the physical
T-count. Arbitrarily many returned dirty helpers can first be included
among the logical wires. The proof requires exact initialized return;
it uses no measurement or postselection.

## 4. A Walsh minor with a large exact denominator

Identify a bit word with its support subset of $`[p]`$. Let
$`K=\sqrt m\,H`$ be the integer Walsh matrix. Define the integer zeta
matrix $`Z_{u,S}=[S\subseteq u]`$. Expanding the product over
$`u\cap v`$ gives

```math
K_{uv}=\sum_{S\subseteq u\cap v}(-2)^{|S|},\qquad
K=Z\,\mathrm{diag}((-2)^{|S|})Z^{\mathsf T}.                   \qquad(8)
```

Set $`r=\lfloor(p-1)/2\rfloor`$ and
$`J=\{u:|u|\le r\}`$, with $`k=|J|`$. J is a downset, so the
restriction of (8) to its rows and columns has no terms indexed outside
J. Ordering by subset size makes $`Z_{J,J}`$ unitriangular. Hence

```math
\det H_{J,J}
 =(-1)^{\sum_{u\in J}|u|}\,
   2^{\sum_{u\in J}|u|-pk/2}
 =\pm2^{-L_p/2}.                                               \qquad(9)
```

The exponent is evaluated without an asymptotic approximation:

```math
\begin{aligned}
L_p&=\sum_{j=0}^{r}(p-2j)\binom pj\\
   &=p\sum_{j=0}^{r}
       \left[\binom{p-1}j-\binom{p-1}{j-1}\right]
     =p\binom{p-1}r.
\end{aligned}
```

The central-binomial estimate gives $`L_p=\Theta(2^p\sqrt p)`$.
Both the factorization and the attained minor are explicit; no numerical
determinant assumption is used.

## 5. The unique full-channel witness

For an ordered subset S of modes define the Hermitian Pauli

```math
\Gamma_S=i^{|S|(|S|-1)/2}\prod_{j\in S}^{\mathrm{increasing}}\gamma_j.
```

Independence of the m Majorana binary labels makes all $`\Gamma_S`$
distinct up to their fixed signs. Expanding conjugated Majorana products
and using orthogonality of H gives the exterior-power identity

```math
W\Gamma_SW^\dagger
 =\sum_{|U|=|S|}\det H_{U,S}\,\Gamma_U.                       \qquad(10)
```

Repeated-mode contributions cancel, and the common Hermitian phase for
the fixed grade appears on both sides. Thus (9) is an actual Pauli-channel
entry of W, proving $`T(W)\ge L_p`$ without a Gaussian-intermediate
assumption.

Conjugation by C_Y multiplies $`\Gamma_R`$ by
$`(-1)^{\sum_{j\in R}Y_j}`$. Using $`H=H^{\mathsf T}`$ for the
first, inverse-W action gives

```math
Q(I_Y\otimes\Gamma_J)Q^\dagger
 =\sum_{\substack{|R|=k\\|U|=k}}
   \det H_{R,J}\det H_{U,R}\,
   Z_{Y_R}\otimes\Gamma_U.                                   \qquad(11)
```

The Y-character distinguishes R and the independent Majorana labels
distinguish U. There are no collisions between the displayed Pauli terms.
The term $`R=U=J`$ has coefficient
$`(\det H_{J,J})^2=2^{-L_p}`$, whose smallest square-root-of-two
denominator exponent is $`2L_p`$. Equations (6)--(7) prove (1), including
all the helper allowances stated at the start.

A Y-dependent scalar on Q does not change (11). More generally the same
witness survives the right gauge $`Q(U_Y\otimes I_A)`$ for arbitrary
$`U_Y`$, since it commutes with the input $`I_Y\otimes\Gamma_J`$.
This also covers every separate operator with the required conjugations
$`RA_jR^\dagger=Z_{Y_j}A_j`$: $`Q^\dagger R`$ commutes with all A_j,
which generate the full carrier matrix algebra, so it is
$`U_Y\otimes I_A`$. If R preserves each computational Y sector, this
remaining freedom is only a diagonal scalar. Fixed Clifford changes of
encoding merely permute and sign the channel entries; they preserve the
same obstruction.

### The denominator method is sharp for these two targets

For any k-by-k minor of K, Cauchy--Binet applied to (8) shows divisibility
by $`2^{a_k}`$, where $`a_k`$ is the sum of the k smallest subset
weights among all p-bit words. Each term selects k diagonal factors,
so this divisibility holds term by term. A normalized H minor therefore
has exponent at most $`pk-2a_k`$. Maximizing this expression includes
all weights strictly below p/2, giving $`L_p`$. The downset minor
attains it, and (10) includes every Pauli-channel grade. Consequently

```math
\mathrm{sde}(\widehat W)=L_p,\qquad
\mathrm{sde}(\widehat Q)=2L_p.                                \qquad(12)
```

For the second upper bound, the canonical Q channel is the product of
the W channel, a signed-permutation Clifford channel and the inverse-W
channel. The witness (11) attains the resulting bound. Equation (12)
does not assert optimal T-counts. It concerns the canonical Q; arbitrary
scalar or right gauges retain the lower bound but may increase other
channel denominators. This denominator invariant alone cannot establish
the Gaussian model's stronger $`mp/2`$ source-count bound.

## 6. What joint compilation and approximation can still change

The exact scalar indicator word has the identity

```math
A_xQ_YA_xQ_Y^\dagger
 =W(\gamma_xC_Y\gamma_xC_Y^\dagger)W^\dagger
 =(-1)^{Y_x}I.
```

Its outer W cancels because the interior is scalar. A native compiler
for the whole word need not realize Q or expose the coefficient (11).
The remaining selected Pauli $`\gamma_x`$ is a coherent selector,
not a free wire translation. This algebraic cancellation alone does not
compile it faster.

There is already a concrete count separation. Apply R91 to the identity
table with m rows and m output bits, and allocate $`B=\Theta(m)`$
dirty work. Its complete indicator has
$`T=O(m)`$, $`G=O(m^2)`$, extra dirty width O(m) and
$`D_T=O(\log(p+2))`$, with zero clean work and literal phase one.
The exact separated readout instead needs $`\Omega(m\sqrt p)`$
T gates, regardless of helper supply. Thus it cannot meet R98's exact
linear-count source envelope, whereas joint compilation already avoids
that count cost. An exact linear-count improvement must bypass the
standalone Q or use another mechanism. An approximate separate Q is
not excluded.

The bound is not stable under ordinary approximation by this argument.
Its witness has magnitude $`2^{-L_p}`$, so it supplies no fixed-error or
inverse-polynomial-error exclusion for an approximate Q. A simple robust
version of the transfer makes the remaining possibility precise. Let V
approximate a fixed $`\widetilde Q`$ from Section 2 to operator error
$`\tau`$ on all Y/carrier/dirty inputs, including the identity action on
returned work. Actual adjoints have the same error. Six-call telescoping
then bounds the complete indicator's operator error by $`6\tau`$;
all Clifford wrappers are exact. For initialized clean work, the same
argument applies to the full isometry on its promised initial subspace,
including its required return. Indeed, if J embeds that subspace and
$`\|VJ-JU\|\le\tau`$, unitarity gives
$`\|V^\dagger J-JU^\dagger\|\le\tau`$. The telescoping bound includes
any return error, rather than discarding the work at the end. If V returns
work exactly, the compiled wrapper does too.

This approximate interface still needs a complete accuracy-dependent
ledger for both orientations, all work and every source call. The transfer
only permits the diagonal scalar gauge of Section 2; the arbitrary right
Y-unitary gauge in Section 5 is a lower-bound allowance. R98 as currently
stated assumes exact queries; neither an approximate-Q compiler meeting
the simultaneous budgets nor its approximate frame transfer is proved here.
No unrestricted T-depth or general query-count lower bound follows.

## 7. Exact checks

[Focused fixtures](../../tests/test_dual_frame_readout_barrier.py) check
the downset determinants and exponent identities, an actual native W
channel entry, the corresponding complete Y-character coefficient,
the six-call native readout wrapper with literal scalar gauges and actual
inverses, and the clean-projector normalization. The proofs above supply
the all-size statements; finite checks do not replace them.
