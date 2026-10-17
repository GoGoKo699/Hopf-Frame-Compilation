# Scalar phases of closed Gaussian circuits have bounded degree

[Original Gauss-sum encoding](GAUSS_SUM_QUERY_ENCODING.md) ·
[Paid inverse extraction](GAUSSIAN_INVERSE_EXTRACTION.md) ·
[Read-only Clifford normal form](READ_ONLY_TWO_LAYER_RIGIDITY.md#1-address-preserving-clifford-interlayers)

An affine Gaussian encoding can hide an arbitrary Boolean table in its
Gauss-sum sign. The sign alone is not visible in its symplectic matrix.
Nevertheless a circuit assembled from a bounded number of the actual
affine Gaussian phase primitives cannot return the work while extracting
an arbitrary such sign. A path-sum argument recovers the information that
ordinary symplectic tracking misses.

The restriction is to the defined Gaussian **round architecture**. It is
not an unrestricted T-depth lower bound, and it does not change R96.

## 1. Primitive, return contract, and conclusion

Let x be an a-bit read-only computational address. A Gaussian phase
round on M physical work qubits is the literal diagonal

```math
D_j(x)|w\rangle=(-1)^{q_j(x,w)}|w\rangle,
```

where $`q_j`$ is quadratic in w, and every coefficient, including
its constant term, is an affine Boolean function of x. Between r such
rounds allow arbitrary fixed native Clifford circuits on the work.
Also allow read-only Clifford circuits on address and work which
preserve every computational address, as described below. Work may
contain any number of initialized and arbitrary dirty qubits. No
address-dependent scalar with an arbitrary prescribed Boolean exponent
is supplied as a free primitive.

**Scalar-return obstruction.** Suppose this circuit has, for every x,
zero-column error

```math
\left\|V_x|0^M\rangle-(-1)^{f(x)}|0^M\rangle\right\|
\lt\kappa,
\qquad \kappa=\sqrt{2-\sqrt2}.
                                                               \qquad(1)
```

Then the algebraic normal form of f satisfies

```math
\deg f\le4r+2.                                      \qquad(2)
```

Consequently the table $`f(x)=\prod_{i=1}^a x_i`$ requires
$`r\ge\lceil(a-2)/4\rceil`$. More generally every table of degree
d requires $`r\ge\lceil(d-2)/4\rceil`$ when this lower bound is
positive. Exact return is included. A full initialized-isometry
approximation to the scalar query implies (1), since the allowed dirty
inputs include the all-zero word. Thus adding initialized helpers,
dirty helpers, or external references does not evade this scoped bound.
The proof needs only the literal zero-column phase; it assumes neither
identity on other work inputs nor a particular work dimension.

The same obstruction applies to an XOR query in this architecture.
Prepare its output in $`|{-}\rangle`$ by a fixed Clifford, execute
the query, and actually undo that preparation. This produces the
zero-column phase $`(-1)^{f(x)}`$ with the same r and the same error.
No clean-output assumption is added to the query: its full-input
contract already includes this particular output state.

Read-only Clifford intervals can be reduced to fixed work Cliffords,
affine-address X and Z controls on work, and a Clifford scalar phase on
the address. This is the separately proved read-only normal form.
Use ordered X/Z products, rather than introducing uncontrolled phase
conventions for Hermitian Paulis. These controls do not increase the
work symplectic degree. Their path-sum treatment in Section 3 introduces
only affine coefficients and, for scalar address CZ gates, a quadratic
constant term. They do not count as Gaussian rounds in (2).

The native implementation of one phase round must still be paid. If
s is the number of cubic monomials after expanding its affine
coefficients, R91's private-helper echo implements those phases in
eight T layers using $`O(s)`$ T gates, elementary Cliffords and
returned dirty helpers. Every lower-degree monomial is also an explicit
Clifford gate, including address-work terms $`x_iw_j`$, address phases,
and constant quadratic work terms. Their number can include
$`O(aM+M^2)`$ terms. Every fixed Clifford
interval is also charged. This does not assert a universal
$`O(\sqrt Q)`$ resource envelope for arbitrary dense Gaussian rounds.
Conversely, a lower bound on the number r of these completed rounds is
not a lower bound on native T-depth after an arbitrary resynthesis.

## 2. A basis-independent formula for the Arf sign

Let q be a quadratic Boolean form on $`2s`$ bits with $`q(0)=0`$.
Write $`q(v)=v^{\mathsf T}Qv`$, and suppose its alternating polar
matrix $`B=Q+Q^{\mathsf T}`$ is nonsingular. Put
$`M=B^{-1}Q`$ and define the second characteristic coefficient

```math
c_2(M)=\sum_{i\lt j}(M_{ii}M_{jj}+M_{ij}M_{ji})
\quad\text{over }\mathbb F_2.
```

Then

```math
\mathrm{Arf}(q)=c_2(B^{-1}Q)+\binom{s}{2}\pmod2.       \qquad(3)
```

To prove this without choosing an address-dependent basis in a circuit,
make a mathematical symplectic basis change with

```math
B=J=\begin{pmatrix}0&I\\I&0\end{pmatrix}.
```

Congruence $`Q\mapsto S^{\mathsf T}QS`$ conjugates $`B^{-1}Q`$
by S, so $`c_2`$ is unchanged. In this basis write

```math
Q=\begin{pmatrix}A&C\\I+C^{\mathsf T}&D\end{pmatrix},
\qquad A=A^{\mathsf T},\quad D=D^{\mathsf T}.
```

The quadratic form is
$`\sum_i u_iv_i+\sum_i A_{ii}u_i+\sum_iD_{ii}v_i`$;
its Arf invariant is $`\sum_iA_{ii}D_{ii}`$.
Meanwhile

```math
JQ=\begin{pmatrix}I+C^{\mathsf T}&D\\A&C\end{pmatrix}.
```

The block formula for its second coefficient gives

```math
c_2(JQ)=c_2(I+C)+c_2(C)
       +\mathrm{tr}(I+C)\mathrm{tr}(C)+\mathrm{tr}(DA)
       =\binom{s}{2}+\mathrm{tr}(DA).
```

Here $`c_2(I+C)=\binom{s}{2}+(s-1)\mathrm{tr}(C)+c_2(C)`$
and $`(\mathrm{tr}C)^2=\mathrm{tr}C`$ over $`\mathbb F_2`$.
The off-diagonal terms of $`\mathrm{tr}(DA)`$ cancel in pairs
because A and D are symmetric. This proves (3) for every representative
Q of q; no choice of upper-triangular representative changes it.

In particular, for a parameter family with coefficient degree at most e
and inverse-polar entry degree at most d,

```math
\deg\mathrm{Arf}(q_x)\le2(d+e).                      \qquad(4)
```

All degrees are Boolean algebraic degrees after reduction
$`x_i^2=x_i`$. The symplectic basis in the proof is only a proof
device. It is not an uncharged native transformation.

## 3. Literal zero-column return makes the path form nonsingular

First suppose all fixed work Cliffords are real, and the zero-column
return is exact. Expand them into real Clifford gates H, CNOT and
Paulis; CZ can be retained as a diagonal phase. Convert each
affine-address X control to an H, affine-address Z control, H word.
Wire values then depend affinely on the Hadamard path variables, but
not on x: all address dependence is in diagonal phases.

Let h be the total number of Hadamard gates in this expanded word.
Each Hadamard introduces one independent new path bit, and contributes
the product of its input and output bits to the phase. Starting at
$`|0^M\rangle`$ and ending with $`\langle+^M|`$ leaves every path
variable summed, with no output constraint. Thus

```math
\langle+^M|V_x|0^M\rangle
 =2^{-(h+M)/2}\sum_{u\in\mathbb F_2^h}(-1)^{q_x(u)}.  \qquad(5)
```

The path phase is quadratic in u. Its nonconstant coefficients are
affine in x. Its constant $`c(x)=q_x(0)`$ has degree at most two,
including any read-only address Clifford phase. Write
$`q_x(u)=c(x)+\bar q_x(u)`$, with $`\bar q_x(0)=0`$, and
let K be the polar matrix of $`\bar q_x`$.

Exact zero-column return makes the left side of (5) equal
$`(-1)^{f(x)}2^{-M/2}`$. The absolute Gauss sum is therefore
$`2^{h/2}`$. A quadratic form of polar rank k has Gauss sum either
zero or absolute value $`2^{h-k/2}`$. Hence $`k=h`$: h is even,
K is nonsingular for every x, and

```math
f(x)=c(x)+\mathrm{Arf}(\bar q_x).                     \qquad(6)
```

When h is zero the path form is empty and $`f=c`$, already of
degree at most two. No initialized path register or postselection is
used physically in this argument. Equation (5) is a mathematical
expansion of a complete circuit amplitude.

## 4. Pauli markers bound the inverse polar degree

For each Hadamard output path bit $`u_j`$, insert a Z at that location
when an independent formal marker $`\lambda_j`$ is one. The markers
are algebraic probes, not extra gates required of the target compiler.
The path sum acquires the linear phase $`\lambda^{\mathsf T}u`$.
Nonsingular quadratic Fourier translation gives the ratio

```math
\frac{\langle+^M|V_{x,\lambda}|0^M\rangle}
     {\langle+^M|V_x|0^M\rangle}
 =(-1)^{\bar q_x(K^{-1}\lambda)}.                    \qquad(7)
```

The polar matrix, in the marker variables, of the exponent on the
right is $`K^{-1}`$.

Move each inserted Z to the output through the **actual original
suffix** of the circuit. It becomes a real Hermitian Pauli

```math
P_j(x)=\sigma_j(x)X(p_j(x))Z(z_j(x)),
\qquad \sigma_j(x)\in\{+1,-1\}.
```

The work symplectic matrix of one Gaussian phase round is affine in x.
Fixed Clifford intervals contribute constant matrices; affine-address
Pauli controls have identity symplectic matrix. Therefore every entry
of $`p_j,z_j`$ has degree at most r, independently of M and h.
No bound on the potentially complicated individual signs
$`\sigma_j(x)`$ is assumed.

In their inherited order, the inserted word is an ordered product of
$`P_j(x)^{\lambda_j}`$ followed by the original V. Exact return
on $`|0^M\rangle`$ suffices here. Ordered real Paulis multiply by

```math
X(p)Z(z)X(p')Z(z')=(-1)^{z^{\mathsf T}p'}
                            X(p+p')Z(z+z').
```

For every p,z,
$`\langle+^M|X(p)Z(z)|0^M\rangle=2^{-M/2}`$.
Consequently the exponent of the ratio (7) is quadratic in the markers.
Its mixed coefficient for two distinct markers is a dot product
$`z_i^{\mathsf T}p_j`$, with i,j in the inherited product order.
The individual signs contribute only linear marker terms and cancel
from every mixed coefficient. Each mixed coefficient therefore has
address degree at most $`2r`$.

Comparing with (7) proves that every off-diagonal entry of
$`K^{-1}`$ has degree at most $`2r`$. Its diagonal entries are zero,
since K is nonsingular alternating. The path coefficient matrix Q has
degree at most one. Equations (3), (4) and (6) now give

```math
\deg f\le\max\{2,2(2r+1)\}=4r+2.
```

This is the step that detects a central scalar phase even though the
completed work symplectic matrix may be the identity. It neither
assumes nor computes the signs of intermediate Clifford lifts.

## 5. Complex Clifford intervals and the approximation threshold

Realify any complex work unitary by

```math
\mathcal R(U)=
\begin{pmatrix}\mathrm{Re}\,U&-\mathrm{Im}\,U\\
\mathrm{Im}\,U&\mathrm{Re}\,U\end{pmatrix}.
```

This is multiplicative. One additional rebit indexes real and imaginary
parts. H and CNOT remain real Clifford gates. An S gate on work qubit j
becomes a controlled $`J=XZ=\begin{pmatrix}0&-1\\1&0\end{pmatrix}`$
on the rebit, implemented by fixed CZ and CNOT gates with that work
control. A fixed global eighth-root phase becomes a real Clifford
rotation on the rebit, generated by HZ. Every real Gaussian phase round
acts trivially on the rebit, so the number r is unchanged.

Affine-address X/Z controls realify in the same way. A scalar
address S phase becomes an affine-address controlled J on the rebit;
an address CZ phase stays a quadratic scalar sign. Thus all nonconstant
path coefficients remain affine in x, and its constant term has degree
at most two. The rebit is a proof device; no extra physical initialized
qubit is required or granted to the compiler. Exact literal return of
the original zero column implies exact return of its realification on
the zero rebit and zero work, so Sections 3–4 apply.

Finally, for each fixed x the actual work circuit is Clifford. Its
realification sends the all-zero state to a real stabilizer state.
A real stabilizer state's computational-zero amplitude is zero or
$`\pm2^{-k/2}`$ for an integer $`k\ge0`$. Unless this state is
exactly the specified $`(-1)^{f(x)}|0\rangle`$, its Euclidean
distance from that target is at least $`\sqrt{2-\sqrt2}`$.
Realification preserves this distance. Therefore the strict bound
(1) forces exact literal zero-column return for every x, reducing the
approximate statement to the exact proof. The threshold is a statement
about this Clifford-fiber architecture, not arbitrary approximate
Clifford+T circuits.

## 6. A genuine scalar loop, and why its Fourier map must be paid

The tempting metaplectic identity is real and exact. For a nonsingular
quadratic q on $`2s`$ work bits with $`q(0)=0`$ and polar B,
define the generalized
Fourier transform

```math
(F_B)_{u,v}=2^{-s}(-1)^{u^{\mathsf T}Bv},\qquad
D_q|v\rangle=(-1)^{q(v)}|v\rangle.
```

Since B is symmetric and invertible, $`F_B^2=I`$. Completing the
quadratic Fourier sum gives

```math
(F_BD_q)^2=\gamma(q)D_qF_B,\qquad
(F_BD_q)^3=\gamma(q)I,\qquad
\gamma(q)=(-1)^{\mathrm{Arf}(q)}.                    \qquad(8)
```

Thus (8) does return every arbitrary work input with the desired Arf
sign. But $`F_B=H^{\otimes2s}C_B`$, where
$`C_B|v\rangle=|Bv\rangle`$ is an **in-place** map. Its Pauli
action includes $`B^{-\mathsf T}`$. The constant-depth doubled
shear for B does not provide this in-place map, just as in the paid
inverse-extraction note.

If a supplied literal native $`C_B`$ has resources
$`(T_B,G_B,D_B,c_B,b_B)`$, with separate clean and dirty helper
reservations and complete return, (8) makes three such
calls, three D calls, and $`6s`$ Hadamards. It has
$`T\le3T_B+3T_q`$, $`G\le3G_B+3G_q+6s`$, and
$`D_T\le3D_B+3D_q`$, with helper pools reused only at completed
boundaries. A sufficient simultaneous allocation is the $`2s`$-bit
arbitrary principal word, $`\max\{c_B,c_q\}`$ clean helpers and
$`\max\{b_B,b_q\}`$ dirty helpers. The native diagonal phase has
$`c_q=0`$; any clean requirement for the supplied in-place map remains
explicit. If each map or phase has full error at most
$`\delta`$, the six-factor replacement error is at most
$`6\delta`$, including leakage. Supplying B as a classical matrix
does not supply a native map at these resources.

In particular, trying to make both Q and $`B^{-1}`$ affine cannot
hide a high-degree table in this loop: (4) bounds its Arf degree by
four. More generally degree d in the inverse polar permits at most
$`2(d+1)`$ in its Arf sign when Q is affine. Replacing the in-place
map by a paired map on two arrays also leaves a second-array action
which must be tracked through the entire loop; it is not automatically
returned by (8).

## 7. Evidence, scope, and the remaining question

The [bounded exact fixtures](../../tests/test_gaussian_scalar_loop_obstruction.py)
check (3) for every nonsingular quadratic form on two and four bits,
larger deterministic samples, and the literal cube identity (8).
They also verify path-sum normalization and marker/inverse-polar
identification on an explicit closed affine Gaussian circuit. The
analytic proof establishes the width-independent degree and error
statements; finite examples do not certify an unrestricted depth bound.

The original inverse-entry encoding and its small trace do not escape
this restriction by closing a bounded number of the actual affine
Gaussian phase rounds into a scalar loop. Its worst-case degree a
requires linearly many such rounds even with a generous clean budget.
This stops that mechanism as a sublogarithmic query candidate. It does
not stop the project: a successful query may use address-conditioned
non-Clifford work dynamics, nonlinear coefficient functions whose
native costs are paid, or another representation outside this round
model. The general dirty-query and optimal complete-frame depth
questions remain open.
