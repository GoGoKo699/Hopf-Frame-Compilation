# Parallel affine detuning: a paid batch and the remaining composition

[Selective-inversion source](SELECTIVE_INVERSION_DIRTY_INDICATOR.md) ·
[Exact dirty phase relocation](READ_ONLY_PHASE_RELOCATION.md) ·
[Polynomial-overhead discovery target](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

Polynomially many affine detuning presentations do have a simultaneous
native schedule. The retained phase-relocation construction implements
their entire commuting product in four T layers, with polynomial address
overhead in count and dirty width. No affine address copy is initialized
or supplied. This pays a batch, not the complete selective operation.

A single such batch, even with arbitrary fixed bank-unitary encoders,
still cannot implement the complete indicator phase: an exact cube
identity gives operator error at least one quarter. The scheduling result
and this obstruction have different roles. The former shows that logical
access count is not native T-depth; the latter locates the missing
noncommuting composition. Section 4 further proves that the read-only
Clifford-interlayer Pauli echo needs at least the ceiling of b/2 rounds
for b at least two below error square root two. These logical rounds do
not lower-bound unrestricted native depth. R96 remains unchanged.

## 1. Four native layers for a complete affine batch

Let $`b\ge1`$, $`m=2^b`$ and $`q\ge1`$. The original address
x and m-bit output bank Y are arbitrary. Give call j its own arbitrary
m-bit signal bank $`s_{j,u}`$. Choose a fixed affine presentation
$`z_j(x)=L_jx+t_j`$, a fixed permutation $`\pi_j`$ of the bank
indices, and a sign $`\sigma_j\in\{1,-1\}`$. Define

```math
\mathcal R=
\prod_{j=1}^q\prod_{u\in\mathbb F_2^b}
 \exp\!\left(\frac{i\pi\sigma_j}{4}
 (|z_j(x)\mathbin\oplus u|+1-Y_{\pi_j(u)})Z_{s_{j,u}}\right).
                                                               \qquad(1)
```

The matrix L need not be invertible: its rows will be used directly as
parity controls, not computed into a prepared address register. The
permutations are specified wire indices in those controls; no
address-dependent router is granted. All factors of (1) commute.

For a Boolean affine parity p, write T(p) for its literal phase
$`\omega^{p}`$, where $`\omega=e^{i\pi/4}`$. The retained
detuning identity gives

```math
\begin{aligned}
\mathcal R={}&\omega^{\lambda m(b+2)/2}
 \prod_{j,u}T(s_{j,u})^{-\sigma_j(b+1)}
 T(Y_{\pi_j(u)}\mathbin\oplus s_{j,u})^{-\sigma_j}\\
&\quad\cdot
 \prod_{j,u,a}T(z_{j,a}(x)\mathbin\oplus s_{j,u})
                         ^{\sigma_j(-1)^{u_a}},
\qquad \lambda=\sum_j\sigma_j .
\end{aligned}                                                   \qquad(2)
```

This includes the exact scalar. The edge signs sum to zero for each
address coordinate, regardless of its affine presentation. Equation (2)
does not assert that its displayed axes are jointly independent.

Put $`\chi=1`$ if b is even and zero otherwise. Extract one
$`T^{-\sigma_j}`$ from each signal power when chi is one, leaving
only an S power on that signal. The number of remaining odd phases is

```math
M=qm(b+1+\chi).
```

Their total linear support size H and affine-constant count C are

```math
\begin{aligned}
H&=qm(b+2+\chi)+m\sum_{j,a}\mathrm{wt}(L_{j,a}),\\
C&=m\sum_{j,a}t_{j,a}.
\end{aligned}                                                   \qquad(3)
```

The first two local types have support sizes one and two; each edge
parity has support $`1+\mathrm{wt}(L_{j,a})`$ and constant
$`t_{j,a}`$.

Apply [the retained four-layer relocation word](READ_ONLY_PHASE_RELOCATION.md#3-parallel-phase-words-and-a-sharp-scalar-example)
to these M signed affine phases. Each phase receives four private
arbitrary dirty wires. Its original parity inputs occur only as CNOT
controls during this word. Across the entire batch the T targets are
distinct, so the four physical T layers run together. Shared x, Y and
signal controls only serialize elementary Clifford gates at their
boundaries; no Clifford-depth claim is made.

Reduce each complete residual signal S power modulo four, in both the
forward word and its actual reversed word. This costs at most 3qm
positive-S gates in either direction. Emit the scalar using
$`(HS)^3=\omega I`$ on an existing private work wire, with its
exponent reduced modulo eight. The forward scalar word costs at most
42 Cliffords, and its actual inverse, expanding inverse S gates into
three positive S gates, costs at most 84. Thus common bounds for the
forward batch and its actual inverse are

```math
\begin{aligned}
T_{\mathcal R}&=15M,\qquad D_T(\mathcal R)\le4,\\
G_{\mathcal R}&\le22M+16H+64C+3qm+84,\\
b_{\rm signal}&=qm,\qquad b_{\rm extra}=4M,\qquad
a_{\rm clean}=0,\qquad\varepsilon=0.
\end{aligned}                                                   \qquad(4)
```

The allowance 3qm uses the stated reduction of entire S powers; naive
gate-by-gate positive-S expansion of an unreduced inverse could be larger.
No such reduction crosses a noncommuting gate. All relocation uncomputes
and the full inverse are actual native inverse words.

The resulting operator is exactly $`\mathcal R\otimes I_w`$
on every input, including correlations among all signals, added work
and external references. The added 4M wires return exactly. The original
signals can be entangled by the intended operator R; their return is
not supplied by (4). Every original computational label is preserved.

For polynomial q, (4) gives
$`T,b=O(qmb)`$ and $`G=O(qmb^2)`$, within the relaxed source
budgets. The original [11-layer construction](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#52-joint-resynthesis-pays-one-special-angle-detuning-call)
retains its different advantage: one call has linear T count and no
additional dirty work. The four-layer batch pays more count and width
to permit simultaneous affine presentations. It is an application of
the existing relocation theorem, not a new native primitive.

## 2. One commuting batch has an exact third-difference obstruction

Let $`\mathcal R_x`$ denote (1) on address sector x. Consider

```math
U_x=B^\dagger\mathcal R_xE,
```

where B and E are arbitrary fixed unitaries on the whole Y/signal/work
space and do not move the address. They need not be Clifford, and can
mix all banks and returned relocation helpers. For $`b\ge3`$,

```math
\begin{gathered}
\left\|\sum_x|x\rangle\langle x|\otimes U_x
            -e^{i\gamma}D_b\otimes I\right\|\ge\frac14,\\
D_b|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle .
\end{gathered}
                                                               \qquad(5)
```

To prove this, fix arbitrary computational Y and signal values only to
describe the phase of R. The sum over all u satisfies
$`\sum_u|z\mathbin\oplus u|=bm/2`$, independently of z.
Consequently the only address-dependent part of the exponent of omega
in (1) is

```math
-2\sum_{j,u,a}\sigma_j(-1)^{u_a}s_{j,u}z_{j,a}(x)
\pmod8.                                                        \qquad(6)
```

An affine Boolean parity has an integer multilinear expansion of degree
at most two modulo four: the coefficient of a product of k distinct
linear input bits is divisible by $`2^{k-1}`$, and adding the affine
constant complements the parity. Multiplication by minus two makes (6)
a polynomial of address degree at most two modulo eight. Therefore its
third mixed coordinate difference vanishes exactly. This holds for
every computational Y/signal input and is an operator identity, not
an initialized-state estimate.

Choose the three coordinate directions $`e_1,e_2,e_3`$, and order the
eight vertices by the Gray cycle

```text
000, 100, 110, 010, 011, 111, 101, 001.
```

The other address bits can be fixed arbitrarily. With these vertices
denoted $`v_0,\ldots,v_7`$, the alternating product is

```math
U_{v_0}^\dagger U_{v_1}U_{v_2}^\dagger U_{v_3}
U_{v_4}^\dagger U_{v_5}U_{v_6}^\dagger U_{v_7}=I .     \qquad(7)
```

Adjacent endpoint encoders cancel, and all R branches commute, so the
remaining expression is E-conjugate to the vanishing third difference.
The corresponding target product is the product of eight distinct
output Z Paulis. It is a nonidentity Hermitian Pauli at distance two
from I. Four positive and four negative occurrences cancel the common
comparison phase gamma. Replacing the eight unitary factors in (7)
therefore gives $`2\le8\tau`$, proving (5).

This floor is independent of q and of the ranks of the accessed affine
features. Reusing a signal bank, or taking inverse calls, does not evade
the commuting-product identity. Noncommuting interlayers between batches,
address motion and address-dependent bank routing are outside its
hypotheses. In particular the four-layer schedule in (4) does not make
a complete selective decoder free.

## 3. The complete composition still needs its own return identity

For h exact batches on a common arbitrary signal space of width n_S,
write $`M_\ell,H_\ell,C_\ell,q_\ell`$ for their parameters in
(3)–(4). Reuse a private pool of $`4\max_\ell M_\ell`$ wires,
since each batch returns that pool exactly. Let the interlayers have
common forward/actual-inverse bounds
$`T_{\rm inter},G_{\rm inter},D_{\rm inter}`$, total full-operator
error eta and additional arbitrary width w, with zero clean input.
Then the complete native word has the simultaneous bounds

```math
\begin{aligned}
T_W&\le15\sum_\ell M_\ell+T_{\rm inter},\\
G_W&\le\sum_\ell(22M_\ell+16H_\ell+64C_\ell
                         +3q_\ell m+84)+G_{\rm inter},\\
D_T(W)&\le4h+D_{\rm inter},\\
b_W&\le n_S+4\max_\ell M_\ell+w,\qquad
\|\widetilde W-W\|\le\eta.
\end{aligned}                                                   \qquad(8)
```

Here eta can be the sum of the interlayer error certificates by
telescoping; the batches themselves are exact. The interlayers must
independently justify their counts, allocation and accuracy. Equation
(8) makes no inference about return of the original signals or private
interlayer work.

For example, set $`P=\prod_{a=1}^{n_S}Z_{s_a}`$. If an ideal
word separately proves

```math
\|P W^\dagger P W-D_b\otimes I\|\le\delta_0,
```

then use the actual emitted adjoint and Hadamards on all Y bits. The
resulting complete dirty indicator satisfies

```math
\begin{aligned}
\epsilon_{\rm ind}&\le\delta_0+2\eta,\\
T_{\rm ind}&\le2T_W,\qquad D_T(\mathrm{ind})\le2D_T(W),\\
G_{\rm ind}&\le2G_W+4n_S+2m,\qquad
b_{\rm ind}\le n_S+4\max_\ell M_\ell+w .
\end{aligned}                                                   \qquad(9)
```

The two P words use 4n_S positive-S gates; the two output Hadamard words
use 2m Cliffords. This error includes every borrowed register's return,
all logical columns and arbitrary references. The first unpaid operation
is the selective identity hypothesized before (9), together with a
qualifying number of noncommuting rounds and interlayer resources.
The general two-band source obstruction and the conditional amplitude
promotion recurrence retain their separate scopes.

## 4. A Pauli echo needs linearly many read-only detuning rounds

This lemma is valid for the following completed-operator grammar. It is
an obstruction to a specified representation, not to unrestricted native
T-depth. In particular, the argument bounds a Pauli echo; it does not
bound the global scalar of an arbitrary completed Clifford-valued word.

### 4.1. Statement and precise grammar

Let x be a b-bit computational address, Y an arbitrary output bank with
one bit per address, and W an entirely arbitrary signal/private register.
Let P be a fixed Hermitian Pauli on W, acting trivially on x and Y.
Consider a word V with at most h detuning rounds, separated and surrounded
by actual Clifford interlayers preserving x and Y pointwise. Those
interlayers may mix all signal and private work, and may have the
address/bank dependence of an actual read-only Clifford. They are not
arbitrary conditionally chosen Cliffords.

Each detuning round is any finite product of full pi/4 detuning calls or
their actual inverses, with known affine presentations of x and fixed
row permutations of Y. Calls in a round are diagonal in the same signal
Z coordinates. Disjoint signal banks suffice. Sharing signal axes is
also allowed when the calls are still such a commuting diagonal product.
No x-dependent lookup or permutation of Y is hidden in the presentation.
All original signals, private work and additional implementation helpers
are arbitrary; an exact returned helper of a detuning implementation can
be included as an identity factor in its completed operator.

Define the actual-adjoint echo

```math
E=P V^\dagger P V,
\qquad
D_b|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle.
```

If for some common phase phi

```math
\|E-e^{i\phi}D_b\otimes I_W\|\lt\sqrt2,
```

then

```math
b\le\max\{2h,1\}.
```

In particular, for b at least two this requires h at least the ceiling
of b/2, irrespective of the number of parallel calls, the rank of their
affine parity features, and finite arbitrary work width. This holds even
when all calls of a round receive the simultaneous four-layer schedule
in section 1.

### 4.2. A valid one-hot test column

Fix Y to the computational table with its only one at an arbitrary
address x_*. This is a test of the full operator promise, not a runtime
preparation or an initialized helper. The selected target becomes the
scalar sign with Boolean exponent

```math
f(x)=[x=x_*],\qquad \deg f=b.
```

The word and P preserve this test subspace. For fixed x, every full
detuning call is a Clifford on W. Up to a scalar, an entire round has
the form

```math
\prod_i S_i^{t_i(x)}.
```

Here t_i is an integer sum, with either sign, of the detunings assigned
to signal i, reduced modulo four. Each detuning is a constant plus a
sum of affine Boolean functions of x. Its squarefree polynomial over
the ring of integers modulo four has degree at most two, and its
reduction modulo two is affine. These two facts survive sums and
sharing signal axes. The scalar of the round is retained in V, but
cancels exactly in every Pauli conjugation used below.

After fixing Y, each actual read-only Clifford interlayer has the form

```math
\lambda(x)X^{u(x)}Z^{v(x)}F,
```

with F fixed Clifford, u and v affine binary labels, and a retained
scalar lambda. Lambda also cancels in conjugation. The actual inverse
has the same grammar, so it suffices to analyze either orientation of
the conjugation by V.

### 4.3. Ordered Pauli degree recurrence, including the literal phase

Use ordered Paulis, not the Hermitian-section tableau convention:

```math
Q(x)=i^{g(x)}X^{a(x)}Z^{c(x)}.
```

A Boolean function of algebraic degree at most d, viewed as a
zero-or-one-valued function modulo four, has a squarefree lift of degree
at most 2d. Indeed, modulo four the XOR of its algebraic-normal-form
monomials is their integer sum minus twice the sum of their pairwise
products. Multiplying a Boolean function by two needs no pair-product
correction, so twice a Boolean degree-d function has degree at most d
over the integers modulo four.

On one signal the exact conjugation rule is

```math
S^{t(x)}X^{a(x)}Z^{c(x)}S^{-t(x)}
 =i^{t(x)a(x)}X^{a(x)}
  Z^{c(x)+(t(x)\bmod2)a(x)}.
```

The label addition is binary. There is no c-dependent extra phase in
this ordered convention. If all incoming label coordinates have degree
at most k-1, the k-th round therefore gives label degree at most k.
Its phase increment is the sum of t_i times the zero-or-one lift of a_i.
That increment has degree at most 2+2(k-1)=2k modulo four.

For completeness, a fixed Clifford F acts on an ordered Pauli with label
v by

```math
F(X^{v_X}Z^{v_Z})F^\dagger
 =i^{\ell(v)+2q(v)}X^{(Mv)_X}Z^{(Mv)_Z},
```

where M is a fixed binary symplectic matrix, ell is integer linear and
q is Boolean quadratic. This follows by multiplying the fixed images
of the ordered X and Z generators: their individual phases give ell,
and reordering their products gives the pairwise sign q. It is an exact
phase identity. When v has Boolean coordinate degree at most k, its
phase increment has degree at most 2k modulo four. Its output label
coordinates still have degree at most k.

The affine Pauli part of a read-only Clifford adds a sign whose exponent
is the symplectic pairing of an affine address label with the moving
Pauli label. Its modulo-four phase has degree at most k+1. Before any
detuning round, a fixed initial P has constant label, so the initial
interlayer gives only a degree-one phase. For k at least one, k+1 is at
most 2k. Consequently the full recurrence is

```math
\deg_{\mathbb F_2}(a),\ \deg_{\mathbb F_2}(c)\le h,
\qquad
\deg_{\mathbb Z_4}(g)\le\max\{2h,1\}.
```

Constant input Pauli phases and every inverse sign are included. No
global scalar of V is estimated or discarded from the original word;
the echo cancels that scalar exactly.

### 4.4. Scalar return and the common-phase approximation floor

For fixed x and the chosen Y, V is Clifford on all arbitrary W inputs.
Thus Q=V-dagger P V is a Hermitian Pauli and E=P Q is a phase times a
Pauli. A nonscalar phase Pauli has two opposite eigenvalues and has
operator distance at least square root two from every unit scalar.
The strict approximation premise therefore forces every E_x to be
scalar. Since P and Q are Hermitian Paulis, this scalar is plus or
minus one, and Q_x is respectively plus or minus P.

Write the fixed P in the ordered convention as

```math
P=i^{g_P}X^{a_P}Z^{c_P}.
```

Then the final labels equal a_P,c_P for every x, and g-g_P is even-valued
modulo four. Boolean-cube Möbius inversion shows that every coefficient
of its unique squarefree modulo-four polynomial is even. Dividing its
coefficients by two gives the Boolean sign exponent e(x) of E_x, with

```math
\deg e\le\max\{2h,1\}.
```

Finally, the corrected values (-1) raised to e(x)+f(x) cannot contain
both signs: the maximum distance of plus and minus one from any common
unit scalar is at least square root two. Hence e(x)+f(x) is constant
throughout the address cube. Its constant does not alter degree, and
the degree-b point indicator gives the claimed inequality.

All reductions use invariant computational test sectors and full work
operators. They include arbitrary references, complete signal/private
return and the actual adjoint. The result does not apply if P changes
the chosen Y sector, if Y moves between calls, if address motion leaves
the read-only grammar, or if non-Clifford interlayers/encoders are added.
It also does not turn a round count into an unrestricted native depth
bound: a new whole-word resynthesis can leave the displayed round
representation.
