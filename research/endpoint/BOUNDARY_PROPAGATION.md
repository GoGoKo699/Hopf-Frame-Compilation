# Boundary propagation for the regular Cayley core

The [regular fixed-tree reduction](TREE_CAYLEY_REDUCTION.md) admits an
exact preconditioner that turns dyadic cut sums into two-point boundary
incidence. Its coefficient matrix has worst-case norm $`\Theta(n)`$,
even inside the regular parameter chart. A nearest-parent propagation
operator has an explicitly known singular gap and a native two-signal
block encoding with $`O(N)`$ T count at the selected endpoint.

The block encoding is a complete unitary with rejected signal amplitudes.
It does not eliminate feedback or return both signal qubits, and therefore
does not establish a compiler for the complete Cayley core.

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

At $`n\ge3`$ and $`\eta=2^{-N}`$, its native approximation uses exactly
$`N+n+7`$ arbitrary dirty qubits, has $`O(N)`$ T count, and approximates
the entire unitary within $`(43/64)\eta`$. The rejected signal components
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

### Exact endpoint reservation and error

For each of the two variable banks choose

```math
q=N+7,\qquad k=n-3.
```

Fix two of the $`n-1`$ parent address bits as sector literals and process
their four values sequentially. The three predicate literals are those
two values and $`c=1`$. The remaining k bits are the free table address.
At every bank call the complete dirty reservation is

```math
\underbrace{N+8}_{\text{precision core}}
+\underbrace{n-3}_{\text{selectors}}
+\underbrace{1}_{\text{predicate helper}}
+\underbrace{1}_{\text{borrowed signal}}
=N+n+7.
```

The borrowed signal is an arbitrary dirty wire, distinct from a,c. The two
supplied clean wires are exactly a,c, so total physical width is
$`N+2n+9`$. The fixed routing and multi-controls run sequentially with the
banks and reuse their dirty pool.

Sector errors take a maximum because the four sectors are invariant and
the other sector words are exactly inactive. The two variable banks compose
with exact routing and literal gates, giving the complete-operator estimate

```math
\boxed{\|\widehat U_A-U_A\otimes I_b\|
 \lt86\,2^{-(N+7)}=\frac{43}{64}\,2^{-N}\lt\eta.}
```

This estimate applies to arbitrary initial a,c values, arbitrary dirty
inputs, and their references. A unitary hybrid uses the complete actual
words; it assumes no reset or intermediate exact return of approximate
work. The selectors and predicate helpers return exactly, and the core
and borrowed signal return within the displayed norm.

There are four sectors in each of two banks. Substitution into the native
cost bound gives

```math
T=8O(2^{n-3}+N+16)+O(n^2)=O(N),
```

```math
G=8O(2^{n-3}(N+7)+N+16)+O(n^2)=O(N^2).
```

All constants are independent of n and the regular parameter table. The
certified table computations terminate without a uniform classical
running-time promise.

## 5. Remaining interface and verification

The native block above has constant normalization for
$`(I-\mathcal B)/2`$. It leaves rejected amplitudes in a,c. Recovering the
physical Cayley unitary still requires a finite native feedback or matrix
function construction that includes the incidence coupling and returns
both signals on every logical/dirty input. The exact gap of order
$`1/n`$ neither supplies this operation nor proves an unrestricted
compiler lower bound.

For an elimination using Q successive approximate bank calls, a uniform
triangle-budget certificate uses per-call precision
$`N+\lceil\log_2Q\rceil+O(1)`$. The increased precision core and each
extra call must be priced within the same dirty pool; the single-block
ledger cannot be silently reused unchanged. A free chain eigenbasis or a
Schur complement expression also does not satisfy the missing interface.

The [boundary matrix helpers](../../compiler_robust_hopf/tree_boundary.py)
and [boundary tests](../../tests/test_tree_boundary.py) verify complete
logical matrices against an independently ordered edge product, unequal
and zero parameters, the triangular and nearest-parent identities, root
chain compression, the all-table gap, all dilation ports and dummy modes,
the literal local factorization, and the exact endpoint reservation.
These finite diagnostics supplement the arbitrary-n proofs; they are not
a native emitter for feedback elimination.
