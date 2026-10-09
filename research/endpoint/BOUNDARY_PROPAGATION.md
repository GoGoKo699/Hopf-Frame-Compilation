# Boundary propagation for the regular Cayley core

The [regular fixed-tree reduction](TREE_CAYLEY_REDUCTION.md) admits an
exact preconditioner that turns dyadic cut sums into two-point boundary
incidence. Its coefficient matrix has worst-case norm $`\Theta(n)`$,
even inside the regular parameter chart. A nearest-parent propagation
operator has an explicitly known singular gap and a native two-signal
block encoding with $`O(N)`$ T count at the selected endpoint.

A complete lossless feedback word reduces the Cayley core to one prefix
encoder, an entrance diagonal, and fixed permutations. The adjoint encoder
cancels on the initialized entrance. The surviving encoder has a direct
$`O(nN)`$-T fine-precision realization. Its
[coarse realization](COARSE_PREFIX_ENCODER.md) and
[cubic collective replacement](COLLECTIVE_PRECISION_REFINEMENT.md)
cost O(N) at their stated error scales; the complete-frame endpoint bound
is unchanged.

## 1. Sparse boundary preconditioning

Use $`N=2^n`$, $`m=N-1`$, $`u=N^{-1/2}\sum_x|x\rangle`$, and the
chronological breadth-first edge order of the fixed-tree reduction. Edge j
has right-half dyadic cut $`T_j=[b_j,b_j+s_j)`$, with
$`s_j=2^{n-d_j-1}`$. Put

```math
P_\perp=I-|u\rangle\langle u|,\qquad
\chi_j=P_\perp\mathbf1_{T_j},\qquad
X=[\chi_1\ \cdots\ \chi_m],\qquad T=\mathrm{diag}(t_j).
```

Let $`P_0|x\rangle=|x-1\bmod N\rangle`$ and
$`R_0=(I-P_0)/2`$. The operator $`R_0`$ is zero on u and equals
$`(I+iA_0)^{-1}`$ on $`u^\perp`$. Direct application of decrement gives

```math
R_0\chi_j=\frac{|b_j+s_j-1\rangle-|b_j-1\rangle}{2}.
```

Define the fixed boundary-incidence matrix and its triangular coefficient by

```math
F=\sqrt2R_0X,\qquad
f_j=\frac{|b_j+s_j-1\rangle-|b_j-1\rangle}{\sqrt2},\qquad
K=X^\dagger R_0X.
```

Then

```math
\boxed{K=\frac12I+\mathrm{strictlower}(F^\dagger F).}
```

To prove triangularity, a later cut has both boundary points inside an
earlier cut or both outside it. The only possible exception would put the
left boundary of a disjoint later cut at the last point of an earlier cut.
That would make the later cut start at an even multiple of its dyadic
length, whereas every right-half cut starts at an odd multiple. Thus the
strict upper triangle is zero. Each diagonal entry is $`1/2`$. Finally,
$`R_0+R_0^\dagger=2R_0^\dagger R_0`$ implies
$`K+K^\dagger=F^\dagger F`$.

For every real parameter table, including zero entries, set

```math
H(t)=(I+2iTK)^{-1}T.
```

The inverse is defined because its matrix is lower triangular with nonzero
diagonal $`1+it_j`$. The complete logical identity is

```math
\boxed{\mathcal C_n(t)=[I-2iFH(t)F^\dagger]P_0.}
```

Indeed, $`A=A_0+2XTX^\dagger`$ on $`u^\perp`$, so the resolvent identity
gives, without an inverse of T,

```math
(I+iA)^{-1}
=R_0-2iR_0X(I+2iTK)^{-1}TX^\dagger R_0.
```

Using $`R_0=-P_0R_0^\dagger`$ yields
$`X^\dagger R_0=-(R_0X)^\dagger P_0`$ and proves the displayed
Cayley identity on $`u^\perp`$. Both sides fix u, so it holds on the
entire logical space, with the literal phase retained.

This preconditioning removes the large cut-sum vectors exactly. It does
not turn H or its physical incidence factors into free unitary operations.

## 2. An exact nearest-parent representation

Index the m edges by heap nodes $`v=1,\ldots,N-1`$. Define

```math
g_v=\frac1{1+it_v},\qquad
h_v=\frac{t_v}{1+it_v},\qquad a_v=\Re h_v,
\qquad D_h=\mathrm{diag}(h_v),\quad D_a=\mathrm{diag}(a_v).
```

On each nonterminal node define child maps

```math
\mathcal B|v\rangle=ih_v|2v\rangle+g_v|2v+1\rangle,\qquad
\mathcal G|v\rangle=-\bar g_v|2v\rangle-i\bar h_v|2v+1\rangle.
```

Both maps are zero on terminal nodes. Each pair of nonzero child columns
is orthonormal, different parents have disjoint supports, and

```math
\|\mathcal B\|,\|\mathcal G\|\le1,\qquad
\mathcal B^\dagger\mathcal G=0,\qquad \mathcal B^n=0.
```

Put $`J_{\mathcal B}=(I-\mathcal B)^{-1}`$. Then

```math
\boxed{H=D_a-iD_hJ_{\mathcal B}D_{\bar h}
                 -iD_hJ_{\mathcal B}\mathcal G D_h.}
```

For the proof, let $`y=Hx`$. Let $`p_v`$ be the inherited charge at the
right endpoint of node v's interval before inserting its own boundary
dipole. The triangular equation gives

```math
y_v=h_vx_v-ih_vp_v,\qquad p_1=0.
```

The dipole subtracts $`y_v`$ at the midpoint and adds it at the right
endpoint. Its two child charges therefore obey

```math
p_{2v}=-y_v=ih_vp_v-h_vx_v,\qquad
p_{2v+1}=p_v+y_v=g_vp_v+h_vx_v.
```

Thus $`p=\mathcal Bp+\mathcal Rx`$, where the child column of
$`\mathcal R`$ is $`(-h_v,h_v)^T`$. The local identity

```math
\begin{pmatrix}-h_v\\h_v\end{pmatrix}
=\bar h_v\begin{pmatrix}ih_v\\g_v\end{pmatrix}
 +h_v\begin{pmatrix}-\bar g_v\\-i\bar h_v\end{pmatrix}
```

gives $`\mathcal R=\mathcal BD_{\bar h}+\mathcal GD_h`$. Substitute
$`p=J_{\mathcal B}\mathcal Rx`$, use
$`J_{\mathcal B}\mathcal B=J_{\mathcal B}-I`$, and use
$`h_v+i|h_v|^2=a_v`$. This proves the formula, including all zero or
unequal parameters. The recurrence also gives the uniform bound

```math
\|H(t)\|\le1+\sqrt2n,
```

because $`\|D_h\|\le1`$, $`\|\mathcal R\|\le\sqrt2`$, and
$`J_{\mathcal B}=\sum_{j=0}^{n-1}\mathcal B^j`$.

## 3. Coefficient growth and the exact propagation gap

### The coefficient norm is of order n

For a uniform table $`t_v=\tau`$, with fixed
$`0\lt|\tau|\lt\sqrt3`$, write

```math
h=\frac\tau{1+i\tau}=a-ib,\qquad
a=\frac\tau{1+\tau^2},\qquad b=\frac{\tau^2}{1+\tau^2}.
```

The representation becomes

```math
H=aI-ibJ_{\mathcal B}-ih^2J_{\mathcal B}\mathcal G.
```

The vectors $`\phi_j=\mathcal B^j|1\rangle`$, $`0\le j\lt n`$,
have unit norm and distinct depths. They are orthonormal and span an
invariant subspace $`\mathcal L`$. Moreover,
$`P_{\mathcal L}J_{\mathcal B}\mathcal G=0`$: reducing an inner
product with powers of $`\mathcal B^\dagger`$ eventually encounters
$`\mathcal B^\dagger\mathcal G=0`$, or the root against a vector
supported below it. Consequently the compression of H to this subspace is

```math
P_{\mathcal L}H|_{\mathcal L}=aI_n-ibJ_n,
```

where $`J_n`$ has ones on and below its diagonal. Apply it to the
normalized all-ones vector. The real and imaginary output parts give

```math
\boxed{\|H(\tau\mathbf1)\|^2
 \ge a^2+b^2\frac{(n+1)(2n+1)}6.}
```

Together with the uniform upper bound, this proves

```math
\boxed{\sup_{t\in(-\sqrt3,\sqrt3)^{N-1}}\|H(t)\|=\Theta(n).}
```

The interior witness $`\tau=1`$ already gives
$`\|H\|\ge n/\sqrt{12}`$. A separate block encoding of H must therefore
have growing normalization. The boundary factor also has a degree-n row,
so $`\|F\|^2\ge n/2`$. In contrast, the complete physical expression obeys

```math
\|2FH(t)F^\dagger\|=\|I-\mathcal C_n(t)P_0^\dagger\|\le2.
```

These are restrictions on separating the factors, not on a coupled native
compiler for the complete unitary.

### Every parameter table has the same singular gap

The root vectors above form a truncated shift chain of length n for every
real table, including zeros. At depth $`d+1`$, the $`2^d`$ child-complement
columns of $`\mathcal G`$ each generate a chain of length $`n-d-1`$.
Distinct chains are orthogonal: reduce inner products with
$`\mathcal B^\dagger`$ to the first unmatched injection, where
$`\mathcal B^\dagger\mathcal G=0`$. Their total dimension is

```math
n+\sum_{d=0}^{n-2}2^d(n-d-1)=2^n-1.
```

They therefore exhaust the node space. On a chain of length r,
$`I-\mathcal B`$ is a bidiagonal difference matrix. Its squared Gram
matrix has diagonal $`(2,\ldots,2,1)`$ and off-diagonal entries -1;
the mixed endpoint sine boundary condition gives singular values

```math
2\sin\!\frac{(2k-1)\pi}{4r+2},\qquad k=1,\ldots,r.
```

Taking the longest chain gives the exact all-table statement

```math
\boxed{\sigma_{\min}(I-\mathcal B)
 =2\sin\!\frac\pi{4n+2}=\Theta(1/n).}
```

The chain spectrum is independent of the parameters. Its basis consists of
parameter-dependent propagated root and complement vectors; that basis
change is not supplied by knowing the spectrum.

## 4. A native two-signal block encoding

Pad the node space with dummy node zero, setting $`\mathcal B|0\rangle=0`$.
This section constructs a complete unitary $`U_A`$ on n node qubits and
two signal qubits a,c such that

```math
\langle0_c0_a|U_A|0_c0_a\rangle=(I-\mathcal B)/2.
```

At $`n\ge3`$ and $`\eta=2^{-N}`$, its native approximation fits
$`N+n+7`$ arbitrary dirty qubits, has $`O(N)`$ T count, and approximates
the entire unitary within $`(43/256)\eta`$. The rejected signal components
are part of the specified unitary, rather than discarded outputs.

### The child-pair bank and every dilation port

Let S map $`|v\rangle`$ to $`|2v\rangle`$ for
$`1\le v\lt N/2`$, and be zero on node zero and all terminal nodes.
Let V fix the node pair $`(0,1)`$ and have disjoint child-pair blocks

```math
Q_v=\begin{pmatrix}ih_v&-\bar g_v\\g_v&-i\bar h_v\end{pmatrix},
\qquad 1\le v\lt N/2.
```

Then $`\mathcal B=VS`$. For $`\alpha_v=\arctan t_v`$, direct
multiplication gives the literal factorization

```math
Q_v=R_y(\pi/2)R_x(\alpha_v)R_z(\alpha_v),\qquad R_y(\pi/2)=XZ.
```

Thus V needs two variable real-rotation banks with fixed Clifford target
conjugations, followed by fixed XZ at nonzero parent addresses. Set
$`\alpha_0=0`$ and omit the fixed XZ at parent zero. Certified sine and
cosine values can be obtained directly as
$`t_v/\sqrt{1+t_v^2}`$ and $`1/\sqrt{1+t_v^2}`$; no argument or equality
test is needed.

Write a mode as $`(a,v)`$. Rotate the complete $`(n+1)`$-bit word left
once, then interchange output modes $`(0,0)`$ and $`(1,1)`$. The resulting
fixed permutation $`D_S`$ has the complete mapping

```math
D_S(a,v)=
\begin{cases}
(1,1),&(a,v)=(0,0),\\
(0,0),&(a,v)=(1,N/2),\\
(\lfloor2v/N\rfloor,(2v+a)\bmod N),&\text{otherwise}.
\end{cases}
```

It is a permutation of all $`2N`$ modes and
$`\langle0_a|D_S|0_a\rangle=S`$, including the dummy and every terminal
column. Hence

```math
U_B=(I_a\otimes V)D_S,\qquad
\langle0_a|U_B|0_a\rangle=\mathcal B.
```

The exceptional output swap is CNOT from a to the node least significant
bit, followed by X on a controlled on all node bits being zero, followed
by the same CNOT. The
[borrowed multi-control construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
gives $`O(n^2)`$ exact T and Clifford cost with one returned arbitrary
helper. An outer-c-controlled wire rotation uses $`O(n)`$ Fredkins; adding
c to the exceptional swap only adds one predicate control. Thus the
controlled $`D_S`$ also has $`O(n^2)`$ exact cost. Its helper reuses a dirty
pool wire and returns before the variable banks.

### Literal sign and occupied signal sectors

Use the complete circuit word

```math
\boxed{U_A=H_cZ_cC_c(U_B)H_c.}
```

The control value is $`c=1`$, and $`Z_c`$ supplies the minus sign.
Ordinary literal Hadamards then give
$`\langle0_c0_a|U_A|0_c0_a\rangle=(I-\mathcal B)/2`$ with no omitted
scalar phase. During the controlled V bank, c is an unchanged predicate
and a is untouched. The controlled fixed XZ factor is the actual
controlled-X/controlled-Z product conditioned on nonzero parent address;
the same borrowed multi-control construction prices it in $`O(n^2)`$.

Use the
[full-operator borrowed-signal real primitive](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations).
At source width q and free address size k it has full-space error less than
$`43\,2^{-q}`$, dirty reservation $`q+k+3`$, and counts

```math
T=O(2^k+q+p^2),\qquad G=O(2^kq+q+p^2),
```

for p predicate literals. Fixed Clifford target conjugations give both
$`R_x`$ and $`R_z`$. The actual primitive is identity on each inactive
predicate sector, including occupied signal and core inputs. The
initialized-only literal phase primitive is not used.

### Repeated calls within the same width

The [self-borrowed whole-word query](../../docs/OPERATOR_SOURCE_COMPILER.md#self-borrowed-whole-word-queries)
needs no external selectors when its m output wires satisfy $`m\ge2k`$.
Use all $`k=n-1`$ parent bits as the address, and use the occupied
dilation port a as the real primitive's arbitrary signal. The full-operator
contract applies even after $`D_S`$ entangles a. The unchanged predicate
is $`c=1`$, so each source-center X has at most the two controls a,c;
no external predicate helper is needed. Fixed routing and the controlled
XZ factor borrow a core wire sequentially and return it exactly before
the next source call. Thus at source width q, with $`m=q+1`$,

```math
\|\widehat U_A-U_A\otimes I_m\|\lt 86\,2^{-q},\qquad
T=O(N+q+n^2),\qquad G=O(Nq+q+n^2).
```

This estimate includes every port, dirty input, and reference. The actual
reversed native word has the same inverse error. False predicate sectors
are exactly identity on the complete occupied workspace.

For R appearances of the block or its actual inverse, choose

```math
h=\lceil\log_2R\rceil,\qquad q=N+h+9,\qquad m=N+h+10.
```

Then $`m\ge2(n-1)`$, and a full-unitary hybrid gives

```math
\boxed{\|\widehat{\mathcal V}-\mathcal V\|
 \lt 86R\,2^{-q}\le\frac{43}{256}\eta.}
```

With $`\rho`$ additional unchanged control wires, use

```math
C_g(U_A)=H_c\,C_g(Z_c)\,C_{g,c}(U_B)\,H_c.
```

On $`g=0`$ the complete word is exactly identity. A source-center
multi-control gate may borrow another core bit as an exact returned helper;
the additional control wires remain disjoint from the core, ports, and node register.
The simultaneous external dirty allocation fits precisely under the
sufficient condition

```math
\boxed{q+1+\rho\le N+n+7
       \quad\Longleftrightarrow\quad h+\rho\le n-3.}
```

These are arbitrary occupied control wires, not supplied clean flags.
Controlled fixed routing costs
$`O(n(\rho+1)^2+(n+\rho)^2)`$ T, which is $`O(N)`$ in this range.
Consequently the direct certificate remains $`O(RN)`$ T and
$`O(RN^2)`$ Clifford. Intervening gates and errors are charged separately.
For $`\rho=0`$, every $`R\le N/8`$ fits; every fixed polynomial in n
therefore fits for sufficiently large n. Exceptional sizes need their own
construction. Saving query selectors removes a width obstruction, not the
repeated precision-source cost or the rejected amplitudes in a,c.

## 5. Complete feedback and the surviving prefix encoder

Extend the heap to all $`2N`$ modes: dummy zero, internal nodes
$`1,\ldots,N-1`$, and leaves $`N,\ldots,2N-1`$. The mode bit a and
n node bits represent this space. In this section $`\mathcal B`$ extends
the same child formula to every internal node and vanishes on the leaves
and dummy. The normalized child column and symmetric complement are

```math
b_v=(ih_v,g_v)^T,\qquad e_v=(g_v,ih_v)^T.
```

They are orthonormal. The wandering space, excluding the dummy, consists
of the root and these $`N-1`$ disjoint complement columns. A vector e
born at depth d generates the orthonormal chain
$`e,\mathcal Be,\ldots,\mathcal B^{n-d}e`$. Different chains are
orthogonal by reducing overlaps with $`\mathcal B^\dagger`$ to the
first orthogonal child columns. Their dimensions sum to $`2N-1`$.
Let R fix zero and reverse every chain. It is a complete Hermitian unitary
mapping each wandering vector to its leaf endpoint, without a separately
normalized inverse.

### Entrance phases and all logical columns

Let V be the full child bank with the SU(2) blocks $`Q_v`$ of Section 4,
now for every $`1\le v\lt N`$. Its second column becomes $`e_v`$ after
multiplication by

```math
d_v=-\frac{1-it_v}{1+it_v}.
```

Let $`\lambda(v)=(2p+1)2^{n-d-1}`$ for $`v=2^d+p`$, and let the
fixed permutation $`H_n`$ send $`\lambda(v)`$ to v and zero to zero.
Put $`v(x)=H_n(x+1\bmod N)`$, and set $`D(x)=d_{v(x)}`$ when
$`v(x)\ne0`$ and $`D(N-1)=1`$. Define the complete entrance by

```math
P_{\rm in}(a,x)=2v(x)+(1-a),\qquad
E=VP_{\rm in}(I_a\otimes D).
```

Initialized a=0 enters only odd modes. Input $`N-1`$ enters the root;
input $`\lambda(v)-1`$ enters $`e_v`$. Each local physical Cayley
factor, in incoming-anchor/fresh-marker order, is literally

```math
\begin{pmatrix}ih_v&g_v\\g_v&ih_v\end{pmatrix}.
```

Earlier edges leave each fresh marker untouched. Its completed column is
therefore its injected complement propagated to the leaves; the root
column is $`\mathcal B^n|1\rangle`$. Reversal supplies exactly these
columns, and $`X_a`$ returns the leaf modes to a=0. With $`J_a`$
inserting the initialized mode bit,

```math
X_aRE(I_a\otimes P_0)J_a=J_a\mathcal C_n(t).
```

This is an equality of complete initialized columns, including zero
rejected output rows and literal scalar phase. The second supplied clean
wire c is ideally unchanged. D is a packed literal diagonal, and V is
one two-axis bank. For routing, increment x, apply $`H_n`$, reorder bits,
and flip the appended bit. The recursion for $`H_n`$ brings the final
bit forward and repeats on the suffix only while the produced prefix is
zero. It uses $`O(n^2)`$ controlled swaps; the exact borrowed-control
construction gives a polynomial in n, hence $`O(N)`$, T cost.

### The adjoint encoder cancels on the entrance

Let $`\mathcal B_0|v\rangle=|2v\rangle`$ on internal nodes. Its
wandering roots are root 1 and the odd children. Let $`T_n`$ fix dummy
zero and map its canonical chains to the weighted chains, choosing the
second column of $`Q_v`$ for each complement. Chain-root phases cancel
in reversal, so

```math
R=T_nR_{0,n}T_n^\dagger,
\qquad \mathcal BT_n=T_n\mathcal B_0.
```

Write $`E_{\rm even}(A)`$ for A on the even-node subspace, with
$`|v\rangle\mapsto|2v\rangle`$, and identity on odd nodes. Then

```math
\boxed{T_n=V_nE_{\rm even}(T_{n-1}).}
```

On every odd node the encoder column is V's second column, so
$`V_n^\dagger T_n`$ is identity there. On an even node, intertwining
gives $`V_n^\dagger T_n|2v\rangle=\mathcal B_0T_{n-1}|v\rangle`$.
This proves the recurrence on all columns, including zero, with
$`T_0=I_2`$. Only the first $`N/2-1`$ parameters enter $`T_{n-1}`$.

Because $`P_{\rm in}J_a`$ has only odd-node support, the adjoint
embedded encoder cancels exactly at that interface. Consequently

```math
\boxed{X_aT_nR_{0,n}P_{\rm in}(I_a\otimes DP_0)J_a
       =J_a\mathcal C_n(t).}
```

Only one encoder survives. The fixed reversal sends the canonical
wandering roots to leaves, so this word uses only the final depth block
of $`T_n`$: the entire target-dependent prefix frame.

The fixed $`R_{0,n}`$ also has a polynomial-size exact word. For a nonzero
binary node write $`0^a\omega0^j`$, with $`\omega`$ beginning and
ending in 1; reversal sends it to $`0^j\omega0^a`$. Reverse all bits,
then reverse the interior between the first and last 1. Unroll the latter
over pairs of endpoint positions, using disjoint unchanged sector
predicates and controlled interior swaps. Dummy zero is fixed. Every
helper returns exactly, and the polynomial T cost is $`O(N)`$.

### Fine precision and its direct native price

At stage j of the recurrence, $`0\le j\lt n`$, a two-axis bank has
$`n-j`$ address bits and a j-bit zero-suffix predicate. The heap mode
occupies the logical register and a; c can be the real primitive's
arbitrary signal. Self-borrowed queries use only the $`q+1`$ core wires
when $`q+1\ge2n`$, and source-center predicates borrow another core bit
as their exact helper. Hence

```math
T(T_n)=\sum_{j=0}^{n-1}O(2^{n-j}+q+(n+j)^2)
       =O(N+nq+n^3),\qquad G(T_n)=O(Nq+nq+n^3).
```

Choose $`q=N+\lceil\log_2(2n)\rceil+9`$. The $`2n`$ variable
banks contribute at most $`(43/512)\eta`$ full-operator error. Their
core fits the declared dirty pool when
$`\lceil\log_2(2n)\rceil\le n-3`$, in particular for $`n\ge7`$.
D can first be compiled at error $`\eta/4`$ with the supplied two clean
wires; subsequent actual unitary words propagate that error without a
reset. Fixed routing is exact. The retained finite-dimensional compiler
handles the finitely many smaller sizes, with their constant cost absorbed
in the asymptotic certificate. Thus this complete word has a direct
$`O(nN)`$ T price and includes final signal and dirty/reference return.
It is weaker than the retained endpoint upper bound.

## 6. Scope and verification

The [coarse prefix encoder](COARSE_PREFIX_ENCODER.md) compiles this same
complete encoder at error $`(43/64)2^{-\lceil N/n\rceil}`$ with
$`O(N)`$ T gates. The
[collective refinement](COLLECTIVE_PRECISION_REFINEMENT.md) gives a cubic
replacement for its terminal frame, with a charged physical baseline and
full initialized-output return. Endpoint precision still requires a linear
total source/table ledger. An independently usable reversal and a second
encoder are unnecessary obligations; the non-small reversal does not
receive small-conjugation error suppression.

The [matrix helpers](../../compiler_robust_hopf/tree_boundary.py) and
[tests](../../tests/test_tree_boundary.py) compare all initialized feedback
columns against the original chronological edge product, including
unequal, signed, and zero parameters. They also check the recurrence,
entrance cancellation, full-unitary reversal, every dilation port and
dummy mode, literal local factors, propagation identities and spectrum,
and exact repeated-call width/error reservations. These finite matrix
diagnostics supplement the arbitrary-n proofs; they do not emit a
fine-precision linear-T compiler.
