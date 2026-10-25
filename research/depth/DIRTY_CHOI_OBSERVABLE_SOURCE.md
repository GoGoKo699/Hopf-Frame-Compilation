# Dirty Choi transport and a sufficient single-probe correction

[HRT character source](FULL_BANK_CHARACTER_TRANSPORT.md) ·
[Returned subset product](ALGEBRA_QUERY_EXTRACTION.md#4-fast-subset-convolution-as-a-native-returned-query) ·
[Crossed cubic kernels](CROSSED_CUBIC_KERNEL_SOURCE.md) ·
[Exact fixtures](../../tests/test_complete_observable_bases.py)

A coherent Choi-basis construction with its paid Z correction reduces to
one quadratic XOR query on the entire arbitrary-input space. It does not
supply inverse evaluation. However, the missing correction can be weakened
to one complete transported probe: that weaker source already supplies the
whole indicator through a four-corner dirty-product identity. Its remaining
phase is the Boolean Walsh dual of an explicitly cubic bent phase.

These are exact identities and a conditional native reduction. No shallow
implementation of the missing phase is supplied, and the unconditional
complete-indicator and general frame frontiers remain unchanged.

## 1. The complete dirty teleportation word has one-query normal form

Let k be any Boolean map on w bits, put $`d=2^w`$, and use three arbitrary
w-qubit registers C, A and B. Define the complete query and basis change

```math
Q_k|a,z\rangle=|a,z+k(a)\rangle,\qquad
G_{AB}=\mathrm{SWAP}_{AB}Q_k(A;B)H_A.
```

Here $`H_A`$ means Hadamards on every A wire. In rightmost-first operator
order, the basis change has the literal action

```math
G|a,z\rangle=d^{-1/2}\sum_y(-1)^{a\cdot y}|k(y)+z,y\rangle.
```

Follow G by the coherent Bell circuit $`\mathrm{CX}(C;A)`$, then $`H_C`$.
Call the resulting computational labels $`C=t,A=u,B=y`$. Apply the two
phases

```math
Z_k|t,u,y\rangle=(-1)^{t\cdot k(y)}|t,u,y\rangle,
\qquad
\mathrm{CZ}_{CA}|t,u,y\rangle=(-1)^{t\cdot u}|t,u,y\rangle.
```

The first is the paid coherent correction
$`Z_k=H_CQ_k(B;C)H_C`$; the second is a Clifford. The complete word

```math
\mathcal T=\mathrm{CZ}_{CA}Z_kH_C\mathrm{CX}(C;A)G_{AB}
```

obeys

```math
\mathcal T|c,a,z\rangle
 =d^{-1}\sum_{t,y}(-1)^{a\cdot y+t\cdot z}
             |t,c+z+k(y),y\rangle.                       \qquad(1)
```

Indeed $`u=c+z+k(y)`$. The exponent before the two corrections is
$`a\cdot y+t\cdot c`$, and adding $`t\cdot k(y)+t\cdot u`$ reduces it
to the exponent in (1). Every computational column has been retained.

Let $`P_{\rm cyc}|u,y,t\rangle=|t,u,y\rangle`$. Equation (1) is exactly

```math
\boxed{\quad
\mathcal T=P_{\rm cyc}H_BQ_k(A;C)\mathrm{CX}(B;C)H_A.
\quad}                                                   \qquad(2)
```

Thus the entire word is Clifford-equivalent to **one** out-of-place query.
No permutation property of k, clean resource, measurement or postselection
was used. For quadratic k this is one Clifford-equivalent diagonal cubic
block. The Fourier representation of its basis does not implement
$`k^{-1}`$ on an arbitrary input register.

If $`Q_k`$ has common forward/actual-inverse bounds $`t,g,d_T,h`$, with g
counting all elementary native gates, (2) has T-count t, T-depth $`d_T`$,
h returned private helpers and at most $`g+9w`$ total gates. The additional
cost is two w-Hadamard layers, w CNOTs and a register cycle implemented by
two w-wire swaps. The principal C, A and B inputs remain arbitrary; their
individual return is not asserted. Reversing (2) gives the actual inverse
with the same bounds.

## 2. The HRT correction does not disappear on the chosen probe

Use the [HRT maps](FULL_BANK_CHARACTER_TRANSPORT.md#1-the-shared-complete-source)
on the augmentation ideal of

```math
R_r=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2),
\qquad
f_r(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1,
\qquad k_r=f_r^{-1}.
```

The physical register has $`w=2^r-1`$ arbitrary coefficients. The inverse
$`k_r`$ is quadratic. The product definition gives the exact identity

```math
f_r(c+c')=f_r(c)+f_r(c')+f_r(c)f_r(c').                    \qquad(3)
```

Write $`F_r|c\rangle=|f_r(c)\rangle`$. For fixed computational resource
labels a,z, (1) gives, on an arbitrary input state,

```math
\mathcal T\bigl(|\psi\rangle|a\rangle|z\rangle\bigr)
 =d^{-1}\sum_{t,u}(-1)^{t\cdot z}
       |t,u\rangle\bigl(Z^aF_rX^{u+z}|\psi\rangle\bigr).   \qquad(4)
```

This describes the coherent output. It does not make the old dirty label z
an available computational control: z is Fourier-encoded in the t register.
Even at $`a=z=0`$, the output retains the branch $`F_rX^u`$.

The missing correction is $`R_s=F_rX^sF_r^\dagger`$, whose complete action
is $`y\mapsto y+f_r(s)+yf_r(s)`$. Let $`\lambda_r(y)=[e_{[r]}]y`$ and
let $`Z_{y,[r]}`$ be the physical Z on the top carrier coefficient. Then

```math
\begin{aligned}
R_sZ_{y,[r]}R_s^\dagger|y\rangle
 &=(-1)^{\lambda_r(y+f_r(s)+yf_r(s))}|y\rangle\\
 &=(-1)^{\lambda_r(f_r(s+k_r(y)))}|y\rangle.              \qquad(5)
\end{aligned}
```

The last equality is (3) with $`f_r(k_r(y))=y`$. The $`Z^a`$ branch phase
commutes with this computational probe; the $`R_s`$ correction does not.
For example, at y equal to zero its phase is
$`(-1)^{\lambda_r(f_r(s))}`$, whose Boolean polynomial has degree r.
This identifies the residual phase; polynomial degree alone is not a
native-depth lower bound.

## 3. One corrected probe supplies an entire nonlinear XOR query

Put $`r=p+1`$, $`m=2^p`$, with $`p\ge1`$. Define the complete phase

```math
\Xi_r|z,y\rangle
 =(-1)^{\lambda_r(f_r(z+k_r(y)))}|z,y\rangle,
\qquad
S_r|z\rangle=(-1)^{\lambda_r(f_r(z))}|z\rangle.
```

The [paid complete quadratic query](FULL_BANK_CHARACTER_TRANSPORT.md#5-the-quadratic-inverse-has-a-paid-xor-query-not-a-paid-in-place-map)
gives

```math
Q_{k_r}(y;z)\,\Xi_r\,Q_{k_r}(y;z)=S_r\otimes I_y.        \qquad(6)
```

Its two occurrences of $`k_r(y)`$ cancel in the phase. The carrier y is
arbitrary and returns; no inverse evaluation or preferred key value is
assumed. An operator-error epsilon in $`\Xi_r`$ remains epsilon in (6).

Split the whole arbitrary key as $`z=a+e_rb`$, where a has the m-1
nonconstant coefficients on the first p labels and b has all m coefficients,
including its arbitrary constant coordinate. Set
$`\mathcal F_p(a)=1+f_p(a)`$ as a full m-coordinate array. Factors containing
$`e_r`$ have pairwise zero products, so

```math
1+f_r(a+e_rb)=\mathcal F_p(a)(1+e_rb),\qquad
\lambda_r(f_r(a+e_rb))=[e_{[p]}]b\mathcal F_p(a).          \qquad(7)
```

Let $`J_B`$ reverse subset labels, $`(J_BB)_S=B_{[p]\setminus S}`$.
The last expression in (7), at $`b=J_BB`$, is the ordinary binary dot
product $`B\cdot\mathcal F_p(a)`$. Therefore

```math
\mathsf A_p=H_BJ_BS_rJ_BH_B:
|a,B\rangle\longmapsto|a,B+\mathcal F_p(a)\rangle.        \qquad(8)
```

This is a complete XOR query on arbitrary B. The constant-coordinate
toggle in $`\mathcal F_p(a)`$ is retained; it is not a prepared physical
constant. The two m-Hadamard layers and two complement permutations cost
at most $`5m`$ Clifford gates: each permutation has m/2 swaps.

## 4. Four corners remove the arbitrary key and both old banks

Let a, B and C be arbitrary extra banks of sizes m-1, m and m. For the
actual p-bit address x, put $`\iota(x)=\sum_i x_ie_i`$. Define

```math
U:B\mapsto B+\mathcal F_p(a),\qquad
V:C\mapsto C+\mathcal F_p(a+\iota(x)),\qquad
P:Y\mapsto Y+BC.
```

U is (8). V encloses the same query with the p-CNOT shift of a by
$`\iota(x)`$ and its actual inverse. P is the
[paid full unital subset-product XOR](ALGEBRA_QUERY_EXTRACTION.md#4-fast-subset-convolution-as-a-native-returned-query)
on disjoint source banks B,C and output Y. All constant coordinates are
present arbitrary wires.

Execute the following **chronological** word:

```math
P,\ U,\ P,\ V,\ P,\ U^\dagger,\ P,\ V^\dagger.           \qquad(9)
```

Writing $`u=\mathcal F_p(a)`$ and $`v=\mathcal F_p(a+\iota(x))`$, the four
product additions are

```math
BC+(B+u)C+(B+u)(C+v)+B(C+v)=uv.
```

Every a,B,C input returns. Moreover (3) gives

```math
uv=\mathcal F_p(a)\mathcal F_p(a+\iota(x))
   =\mathcal F_p(\iota(x))
   =\prod_i(1+x_ie_i),                                  \qquad(10)
```

because $`\mathcal F_p(a)^2=1`$. Thus Y gains the full monomial vector
$`v_S(x)=\prod_{i\in S}x_i`$, independently of all old bank values.

Conjugate (9) by the fixed upper subset-zeta transform on Y. It is an
involution over the binary field and the two appearances cost pm CNOTs.
The identity

```math
[x=u]=\sum_{S\supseteq u}\prod_{i\in S}x_i
```

shows that the resulting addition is exactly $`e_x`$. This proves the
complete indicator, all work return and every reference extension.
No zero key, zero coefficient bank or branch acceptance was used.

## 5. Simultaneous native ledger and the derivative-only contract

Suppose a native phase source approximates $`\Xi_r\otimes I_h`$ within
complete operator error epsilon and has common forward/actual-inverse bounds
$`t_\Xi,g_\Xi,d_\Xi,h_\Xi`$. Let the exact query $`Q_{k_r}`$ have bounds
$`t_k,g_k,d_k,h_k`$, and the exact product P have
$`t_M,g_M,d_M,h_M`$. Here **g counts every elementary native gate**.

The cited arithmetic proofs list Clifford counts G separately from T;
use $`g=G+T`$. In their notation,

```math
\rho_j=\binom{j+2}{2}2^j,\qquad s_j=j(j+1)2^{j-1},
```

so valid all-gate bounds and depth bounds are

```math
\begin{aligned}
g_k&\le(r-1)(152\rho_r+64s_r)+3(2^r-1),&d_k&\le40,\\
g_M&\le152\rho_p+64s_p+2m,&d_M&\le40.
\end{aligned}
```

The remaining count and width bounds are the cited simultaneous bounds.
The completed indicator from (6)--(10) has

```math
\begin{aligned}
T_I&\le4t_\Xi+8t_k+4t_M,\\
D_T(I)&\le4d_\Xi+8d_k+4d_M\le4d_\Xi+480,\\
G_{\rm all}(I)&\le4g_\Xi+8g_k+4g_M+pm+20m+4p,\\
w_I&\le5m-2+\max(h_\Xi,h_k,h_M),\\
\varepsilon_I&\le4\epsilon,\qquad a_{\rm clean}=0.        \qquad(11)
\end{aligned}
```

The allocated principal dirty wires are a (m-1), B and C (2m), and the
arbitrary carrier y in (6) (2m-1). The returned private pools can be reused
sequentially. The 4p CNOTs are the two key-shift pairs in V and its inverse.
Four occurrences of (8) contain four phase sources and eight quadratic
queries; (9) contains four product queries. Its four phase-source errors
telescope to $`4\epsilon`$ using actual adjoints, without assuming exact
intermediate return. The completed actual inverse has the same bounds.

A weaker source promise is enough. Let a native W on z,y and arbitrary
private work obey only

```math
\|WZ_{y,[r]}W^\dagger-\Xi_r\otimes I_{\rm work}\|_{\rm op}
 \le\epsilon.                                           \qquad(12)
```

W need not match the full correction $`R_z`$ or separately return its work.
Insert its actual $`WZ_{y,[r]}W^\dagger`$ as the phase source. For common W
bounds $`t_W,g_W,d_W,h_W`$, this yields

```math
\begin{aligned}
T_I&\le8t_W+8t_k+4t_M,\\
D_T(I)&\le8d_W+8d_k+4d_M\le8d_W+480,\\
G_{\rm all}(I)&\le8g_W+8g_k+4g_M+pm+20m+4p+8,\\
w_I&\le5m-2+\max(h_W,h_k,h_M),\\
\varepsilon_I&\le4\epsilon,\qquad a_{\rm clean}=0.        \qquad(13)
\end{aligned}
```

The eight additional gates are four literal $`Z=S^2`$ insertions. The
error in (13) is the derivative error in (12), not an assumed error of W
as a standalone correction.

The exact arithmetic overhead has constant native depth and
$`O(m\,\mathrm{poly}(p))`$ T-count, total gates and dirty width. Thus either
source with qualifying simultaneous bounds, depth $`o(\log(p+2))`$ and
every requested fixed tolerance supplies the missing relaxed indicator.
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) then
normalizes count and width. The
[fixed-precision frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md)
retains its own log-log contribution. These are sufficient reductions;
no converse from an arbitrary indicator to $`S_r`$ or $`\Xi_r`$ is claimed.

## 6. The remaining phase is a cubic bent dual

The phase in (7) has a useful explicit structure. Let $`n=m-1`$ and
write $`d_S=b_{[p]\setminus S}`$ for nonempty S. Then

```math
\lambda_r(f_r(a+e_rb))=b_{[p]}+g_p(a,d),\qquad
 g_p(a,d)=d\cdot f_p(a),                                 \qquad(14)
```

where a,d each contain n arbitrary bits. The linear term $`b_{[p]}`$ is
the top coefficient of the original key and is a single paid Z.
For binary Fourier labels u,v, summing d first gives

```math
2^{-n}\sum_{a,d}(-1)^{g_p(a,d)+u\cdot a+v\cdot d}
 =(-1)^{u\cdot k_p(v)}.                                 \qquad(15)
```

Indeed that sum vanishes unless $`f_p(a)=v`$, whose unique solution is
$`a=k_p(v)`$. Thus $`g_p`$ is bent, meaning that its normalized Boolean
Walsh transform has modulus one everywhere, and its dual is
$`g_p^*(u,v)=u\cdot k_p(v)`$. This dual has Boolean degree at most three.
Its complete diagonal phase is paid by the existing quadratic XOR query
and Hadamards on the u register.

The Walsh identity does not implement the diagonal phase $`D_{g_p}`$.
For $`\mathcal H=H^{\otimes2n}`$, double Walsh inversion instead gives

```math
\langle z'|\mathcal H D_{g_p^*}\mathcal H|z\rangle
 =2^{-n}(-1)^{g_p(z+z')}.                               \qquad(16)
```

Every matrix entry has the same nonzero modulus: this is a dense unitary
kernel, not the desired diagonal operation. The distinction is also central
to the [crossed-kernel interface](CROSSED_CUBIC_KERNEL_SOURCE.md#2-the-complete-kernel-and-collision-condition).
Equation (2) shows why the standard dirty Choi completion alone does not
remove it.

The concrete positive question is whether $`D_{g_p}`$, whose Walsh dual
$`g_p^*`$ is cubic, admits a shallow **complete diagonal implementation**,
or whether a source can meet (12) directly with the simultaneous dirty-work
and error bounds.
The four-corner reduction proves the payoff of that weaker phase target
without requiring an in-place HRT permutation or a full coherent correction.
The [exact fixtures](../../tests/test_complete_observable_bases.py) check
finite instances and ledgers; the displayed identities prove the all-size
claims.
