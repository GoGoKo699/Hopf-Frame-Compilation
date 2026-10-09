# A regular fixed-tree Cayley reduction

Every complete real Hopf frame has an exact reduction to a fixed-tree
Cayley family with a bounded real parameter table. Two literal diagonal
unitaries and three explicitly charged permutations supply the surrounding
circuit. They cost $`O(N)`$ T gates at the selected endpoint; the remaining
condition is a native compiler for the Cayley family itself.

Use $`N=2^n`$, $`m=N-1`$, and the
[prescribed frame convention](../../docs/HOPF_INTERFACE.md). All statements
hold for arbitrary effectively supplied real angles, including unequal
branches and singular angles. Products act rightmost first. In particular,

```math
R_y(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix}.
```

## 1. A division-free phase gauge

List the tree edges in chronological breadth-first order. At depth d and
prefix p, the ordered endpoints are

```math
(a_{d,p},b_{d,p})=
\bigl(p2^{n-d},(2p+1)2^{n-d-1}\bigr),
\qquad 0\le d\lt n,\quad 0\le p\lt 2^d.
```

Each edge joins an existing anchor to a fresh marker. Indeed, the
2-adic valuation of the nonzero label b identifies its unique depth,
whereas a is zero or a marker of smaller depth. Pairs within a depth
are disjoint. Define

```math
v_j=(|a_j\rangle-|b_j\rangle)/\sqrt2,\qquad
z_j=e^{2i\theta_j},\qquad
F_j=I+(z_j-1)v_jv_j^\dagger,
\qquad P(\theta)=F_m\cdots F_1.
```

If $`S_j`$ is the literal basis transposition of the two endpoints, then
$`F_j=\exp(i\theta_j(I-S_j))`$. Every factor fixes the uniform vector
$`u=N^{-1/2}\sum_x|x\rangle`$.

There are effectively evaluable literal diagonals such that

```math
\boxed{W(\theta)=D_{\rm out}P(\theta)D_{\rm in}^\dagger.}
```

To construct them, set the initial and current phase at label zero to
one. At edge j assign the previously untouched marker the initial phase
$`i`$ times the current phase of its anchor. Multiply both current
endpoint phases by $`e^{-i\theta_j}`$. The initial phases form
$`D_{\rm in}`$ and the final phases form $`D_{\rm out}`$.

For the proof, let g be the current anchor phase. The two active phases
are $`g(1,i)`$, and direct multiplication gives

```math
R_y(\theta)\mathrm{diag}(g,ig)
=\mathrm{diag}(ge^{-i\theta},ige^{-i\theta})
 e^{i\theta}\begin{pmatrix}\cos\theta&-i\sin\theta\\
-i\sin\theta&\cos\theta\end{pmatrix}.
```

The last matrix, including its scalar, is F on that pair. All input
phases may be assigned in advance because earlier edges fix a future
marker. Thus $`R_jD_{j-1}=D_jF_j`$ telescopes to the boxed identity on
every column. The construction uses $`O(N)`$ phase multiplications;
certified sine/cosine enclosures suffice, with no division, argument
selection, or equality test.

At $`\theta_j=\pi/2`$, each F is exactly its transposition. The product
is modular decrement $`P_0|x\rangle=|x-1\bmod N\rangle`$. To see the
direction, write $`N=2h`$: the root swaps 0 and h, and the remaining
two disjoint subtree products decrement separately in each half. This
sends 0 to $`N-1`$, h to $`h-1`$, and every other x to $`x-1`$.

## 2. The polynomial resolvent and physical Cayley generator

Put

```math
V=[v_1\ \cdots\ v_m],\quad G=V^\dagger V=I+L+L^\dagger,
\quad L=\mathrm{strictlower}(G),\quad
D=\mathrm{diag}(z_j-1).
```

Expansion of the ordered product along increasing edge-index chains gives

```math
\boxed{P=I+VMV^\dagger,\qquad
M=(I-DL)^{-1}D=\sum_{r=0}^{m-1}(DL)^rD.}
```

The inverse is a finite polynomial because DL is strictly lower
triangular. This formula includes $`z_j=1`$ directly. Tree incidence
columns are independent: deleting a leaf forces its incident coefficient
to vanish, and induction finishes the proof. Hence unitarity of P implies

```math
M+M^\dagger+MGM^\dagger=0.
```

For the Cayley form, first assume $`z_j\ne1`$. Let R be any fixed real
orthonormal basis of $`u^\perp`$ and set $`B=R^\dagger V`$. For example,
R can consist of the nonuniform columns of the global Walsh transform.
Since $`D^{-1}=-I/2-i\mathrm{diag}(\cot\theta)/2`$,

```math
M=-2\bigl[G+L-L^\dagger+i\mathrm{diag}(\cot\theta)\bigr]^{-1}.
```

It follows that

```math
R^\dagger PR=I-2(I+iA)^{-1},\qquad
A=B^{-\dagger}\bigl[\mathrm{diag}(\cot\theta)
                         +i(L^\dagger-L)\bigr]B^{-1}.
```

A is Hermitian, so $`\|(I+iA)^{-1}\|\le1`$. Section 3 supplies
effective, nonsingular parameters for this same fixed tree for every
original angle tuple.

There is a direct physical interpretation of $`B^{-1}`$. Removing edge j
separates the dyadic subtree

```math
T_j=[b_j,b_j+2^{n-d_j-1})\cap\mathbb Z.
```

If $`Vf=x`$ and $`x\perp u`$, summing over that subtree cancels every
internal edge and leaves

```math
f_j=-\sqrt2\sum_{x'\in T_j}x_{x'}.
```

Write $`P_\perp=I-|u\rangle\langle u|`$ and
$`\chi_j=P_\perp\mathbf1_{T_j}`$. Extending the reduced generator by
zero on u gives

```math
A_{\rm full}(t)=A_0+2\sum_j t_j|\chi_j\rangle\langle\chi_j|,
\qquad t_j=\cot\theta_j,
```

where $`A_0=RB^{-\dagger}i(L^\dagger-L)B^{-1}R^\dagger`$ is fixed.
At $`t=0`$, its Cayley transform is the decrement permutation on
$`u^\perp`$. Define the full physical Cayley family by

```math
\mathcal C_n(t)=I-2P_\perp(I+iA_{\rm full}(t))^{-1}P_\perp.
```

It fixes u and is unitary on the entire logical space.

The coupled expression is essential for normalization. At the all-swap
baseline, $`M=-2(I+2L)^{-1}`$. Its principal block on the rightmost
root-to-leaf chain is $`-2J_n`$, where $`J_n`$ is lower triangular
with every entry on and below the diagonal equal to one. Consecutive
chain incidence vectors overlap by $`-1/2`$; an increasing-index overlap
path cannot leave the chain and return. Therefore

```math
\|M\|\ge 2\|J_n\|=\csc\!\left(\frac{\pi}{4n+2}\right),
\qquad \|VMV^\dagger\|=\|P_0-I\|=2.
```

The displayed singular value follows from the bidiagonal difference
matrix $`J_n^{-1}`$. Thus normalizing an independently implemented M
would introduce a growing factor even at an inexpensive permutation.

## 3. Coarse swaps give a uniformly regular fixed-tree chart

Choose bits $`e_j\in\{0,1\}`$ and put
$`\delta_j=\theta_j-e_j\pi/2`$. The local identity
$`F_j(\theta_j)=S_j^{e_j}F_j(\delta_j)`$ is literal. Moving the swaps
to the left gives

```math
P(\theta)=\widehat\Pi_m
 \prod_{j=m}^{1}\widehat\Pi_{j-1}^\dagger F_j(\delta_j)
                         \widehat\Pi_{j-1},\qquad
\widehat\Pi_j=S_j^{e_j}\cdots S_1^{e_1}.
```

The residual edge attaches the same fresh marker to
$`\widehat\Pi_{j-1}^{-1}(a_j)`$. At any depth, earlier swaps at that
depth are disjoint from the current edge.

Let $`\Pi_d`$ be the coarse permutation on d prefix bits obtained from
depths below d, with trailing zeros suppressed. In particular,
$`\Pi_n=\widehat\Pi_m`$. Define permutations recursively by

```math
\sigma_0(0)=0,\qquad
\sigma_{d+1}(2p)=2\sigma_d(p),\qquad
\sigma_{d+1}(2p+1)=2\Pi_d(\sigma_d(p))+1.
```

Set $`r_d=\Pi_d\sigma_d`$ and
$`\phi_{d,p}=\delta_{d,r_d(p)}`$. Then

```math
\boxed{P(\theta)=\Pi_n\sigma_nP(\phi)\sigma_n^\dagger.}
```

Indeed, the residual edge originally indexed p has endpoints
$`(2\Pi_d^{-1}(p),2p+1)`$. The recursion sends the canonical edge
indexed p to precisely the residual edge indexed $`r_d(p)`$. Its even
extensions preserve all earlier vertices. Within-depth factors commute,
so the within-depth parameter permutation proves the full identity.

The bits can be selected effectively so that all remaining cotangents
are bounded. Compute certified enclosures for
$`z_j=(\cos\theta_j+i\sin\theta_j)^2`$ and refine until either strict
test is certified:

```math
e_j=0\quad\text{if }\Re z_j\lt 1/2,\qquad
e_j=1\quad\text{if }\Re z_j>-1/2.
```

Use a fixed priority if both pass. These overlapping open tests cover
the unit circle, so selection terminates even at either boundary and
at every original singular input. With $`w_j=(-1)^{e_j}z_j`$,

```math
\Re w_j\lt 1/2,\qquad |w_j-1|>1,\qquad
\cot\delta_j=i\frac{w_j+1}{w_j-1}\in(-\sqrt3,\sqrt3).
```

There is no argument extraction or division by a possible zero. Permute
these real values within each depth by $`r_d`$ to obtain the table t.
Equivalently, its phases are $`e^{2i\phi_j}=(t_j+i)/(t_j-i)`$. The
resulting exact all-angle reduction is

```math
\boxed{W=D_{\rm out}\Pi_n\sigma_n
             \mathcal C_n(t)\sigma_n^\dagger D_{\rm in}^\dagger,
       \qquad |t_j|\lt \sqrt3.}
```

## 4. The permutations have exact native linear cost

A depth-k coarse layer on d prefix bits is a Boolean-table-controlled X
on bit k, conditioned on its lower suffix being zero. Let h be the suffix
predicate, z a separate arbitrary dirty helper, and c an optional extra
unchanged control. Let G toggle z by h and Q toggle the logical target
by $`zce(x)`$, omitting c when uncontrolled. The chronological word
$`G,Q,G,Q`$ toggles the target by

```math
(z\mathbin\oplus h)c e(x)\mathbin\oplus zc e(x)=hc e(x)
```

and restores z. Use the
[exact dirty-table query](../../docs/OPERATOR_SOURCE_COMPILER.md#3-dirty-programming-and-a-one-flag-scalar-block)
for Q and the
[borrowed multi-control construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
for G. During G the logical target is the returned borrowed bit; it is
restored before Q acts. Q's address and separate dirty selectors are
disjoint from its target. The controlled table has at most $`4\cdot2^k`$
rows and $`k+2`$ selectors. Suffix toggles use
$`O((d-k-1)^2)`$ Toffolis. Thus

```math
T(C_c(\Pi_d)),\ G(C_c(\Pi_d))=O(2^d+d^3),
```

with at most $`d+2`$ arbitrary external work qubits. The uncontrolled
permutation uses one fewer address bit. All words are literal exact
Clifford+T permutations, including their scalar phases and work return.

The recursion for sigma has the circuit form

```math
\sigma_{d+1}=C_{\text{new bit}=1}(\Pi_d)(\sigma_d\otimes I_2).
```

Here $`\Pi_d`$ is resynthesized on its d-bit prefix, so its suffix tests
involve only bits in that prefix. The new control is disjoint from every
target and suffix. It is incorporated into the Boolean table; no control
of an already decomposed T gate is assumed. Extra lower logical bits are
untouched. Consequently

```math
T(\sigma_n),\ G(\sigma_n)
=O\!\left(\sum_{d=1}^{n-1}(2^d+d^3)\right)
=O(N+n^4)=O(N).
```

The last absorption is uniform because $`\sup_{n\ge1}n^4/2^n\lt \infty`$.
The same linear bound holds for $`\Pi_n`$ and the actual inverse word
$`\sigma_n^\dagger`$. Controlled prefix calls have $`d\le n-1`$ and
standalone $`\Pi_n`$ uses one fewer address bit, so at most $`n+1`$
arbitrary dirty wires suffice throughout. Neither clean flag is needed.
All work returns on every input, including reference entanglement.

## 5. Two diagonal tables describe the variable generator

Let Q be the balanced complete Hopf frame, whose columns are u and the
normalized Haar wavelets on binary dyadic intervals. For a dyadic subtree
T, let $`I_T`$ be its computational support projector and
$`u_T=\mathbf1_T/\sqrt{|T|}`$. The zero-sum subspace supported in T has
projector

```math
R_T=I_T-|u_T\rangle\langle u_T|.
```

It is the sum of exactly the Haar-wavelet projectors whose intervals
are contained in T. Therefore

```math
2\sum_T t_T|\mathbf1_T\rangle\langle\mathbf1_T|
=D_a-QD_bQ^\dagger,
```

where the real diagonal entries are

```math
a_x=2\sum_{T\ni x}t_T|T|,\qquad
b_I=2\sum_{T\supseteq I}t_T|T|,\qquad b_u=0.
```

This is the identity $`|\mathbf1_T\rangle\langle\mathbf1_T|=|T|(I_T-R_T)`$
summed over the tree cuts. The physical generator is

```math
\boxed{A_{\rm full}(t)=A_0+
              P_\perp(D_a-QD_bQ^\dagger)P_\perp.}
```

Q has an exact native $`O(n^3)=O(N)`$ T-gate implementation. At each
depth all angles equal $`\pi/4`$, so its layer is
$`R_y(\pi/4)=HZ`$ controlled only on a zero suffix. Each controlled Z
is a multi-controlled Pauli. For H, the exact Clifford+T word
$`V_0=SHTHS^\dagger`$ obeys $`V_0 ZV_0^\dagger=H`$, reducing its
control to the same multi-controlled Pauli. The retained borrowed
construction costs $`O(n^2)`$ per layer and returns its helper. Thus the
fixed Haar basis is explicitly charged.

The entries of $`D_a,D_b`$ may grow with N even though t is bounded.
Their diagonal exponentials are covered by the literal diagonal compiler.
The displayed identity is a sum of Hermitian operators, not a
factorization of its Cayley transform into diagonal unitary gates.

## 6. Exact endpoint implication and verification

For $`n\ge3`$, set $`\eta=2^{-N}`$ and use two clean flags and
$`b=N+n+7`$ arbitrary dirty qubits, for physical width $`N+2n+9`$.
Suppose the family $`\mathcal C_n(t)`$, for every effectively supplied
$`t\in(-\sqrt3,\sqrt3)^{N-1}`$, admits a deterministic finite native
Clifford+T compiler in this workspace with T count $`O(N)`$ and complete
initialized-isometry error at most $`\eta/2`$. This condition includes
leakage, dirty/reference return, and all logical inputs.

Then the complete real frame admits an $`O(N)`$ T-count compiler in the
same workspace. Compile each exterior diagonal at error $`\eta/8`$
using the [packed diagonal theorem](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames),
with precision $`\ell=N+3`$. Its sufficient dirty reservation is

```math
n+1+\left\lceil\frac{N+7}{2}\right\rceil\le N+n+7.
```

Run the diagonal, permutation, and core words sequentially, reusing the
same flags and dirty pool. The permutations are exact on occupied as
well as initialized work. A unitary hybrid compares each approximate
call only on the ideal initialized input and propagates previous error
by the norm-one actual word. Hence total error is at most

```math
\eta/8+\eta/2+\eta/8=3\eta/4\lt \eta.
```

No intermediate reset or clean reinitialization is used. All supplied
angle computations and chart choices terminate by their certified
enclosures; no uniform evaluation-time bound is assumed.

The remaining hypothesis is specifically an O(N)-T native realization of
$`\mathcal C_n(t)`$. The bounded inverse norm and the fixed Haar
generator formula do not price that operation or its nonunitary
compression. Thus this reduction preserves the selected endpoint's
existing bounds while isolating a bounded, fixed-tree synthesis family.

The [boundary-propagation construction](BOUNDARY_PROPAGATION.md) gives an
exact sparse preconditioner and a native $`O(N)`$-T two-signal block for
the associated shifted propagation operator. Its proof includes the
coefficient norm, the exact all-table singular gap, and the repeated-call
workspace and error ledger. Its complete feedback identity cancels the
entrance against the inverse encoder, leaving one prefix-frame encoder
with a direct O(nN)-T fine-precision realization.
The [coarse encoder](COARSE_PREFIX_ENCODER.md) realizes that family in
linear T-count at coarse precision; linear-T fine synthesis remains the
unresolved part of the core hypothesis above.

The [exact tree phase and chart tests](../../tests/test_tree_phase_chart.py)
check all-column phase gauges, singular polynomial resolvents, decrement
order, coarse-edge conjugacy, cotangent selection, and the Haar cut
identity over rational or cyclotomic arithmetic. Finite tests supplement
the arbitrary-n proofs above; they are not a native implementation of
the conditional Cayley core.
