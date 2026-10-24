# Moved zeta, order-three batches and a complete Pauli-transport interface

[Working checkpoint](../../WORKSPACE.md) ·
[Addressed-transport boundary](ENTANGLING_SELECTOR_ACCEPTANCE.md#address-transport-cost) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

Moving a fixed zeta transform with the address is not a free extension of
fixed-row readout: two uses already give the complete dirty indicator.
This note reduces moved zeta to disjoint binary order-three blocks and
identifies a weaker sufficient interface: transport one fixed Pauli,
on every input, rather than approximate the whole block permutation.
The fixed basis changes and a linear-depth native baseline are paid.
No sublogarithmic compiler for the remaining interface is supplied.

<a id="order-three-blocks"></a>

## 1. Exact binary reduction

For $`p\ge1`$ address bits let $`m=2^p`$, $`q=m/2`$, and use m arbitrary
selector bits s. All identities below include arbitrary external references.
For a binary invertible matrix Q, write $`U_Q|s\rangle=|Qs\rangle`$.
This is a literal phase-free computational-basis permutation, not a
projective Clifford representative. Define the binary matrices

```math
L=\begin{pmatrix}1&0\\1&1\end{pmatrix},\qquad
U=L^{\mathsf T},\qquad
A=UL=\begin{pmatrix}0&1\\1&1\end{pmatrix}.
```

All matrix arithmetic in this section is over $`\mathbb F_2`$.
We have $`A^3=I`$ and $`A^2+A+I=0`$. The fixed zeta is
$`M=L^{\otimes p}`$; it costs pm/2 CNOTs. Let T_x permute physical
row labels by $`u\mapsto u\oplus x`$. Since $`XLX=U`$,

```math
M_x=T_xMT_x=F(x)M,\qquad
F(x)=\prod_{i=1}^p A_i^{x_i},                         \qquad(1)
```

where A_i acts on the ith tensor coordinate. These A_i commute as
binary linear maps. This fact alone does not parallelize their native
controlled implementations on shared address and selector wires.

Let $`\omega^2+\omega+1=0`$ in $`\mathbb F_4`$. Multiplication
by omega in the binary basis (1,omega) has matrix A. Over this field,
the simultaneous eigencharacters of the A_i are indexed by
$`\sigma\in\{+1,-1\}^p`$, with eigenvalues $`\omega^{\sigma_i}`$.
Frobenius pairs sigma with minus sigma. Each pair descends to a binary
two-dimensional space, so one fixed $`S\in\mathrm{GL}(m,2)`$ gives

```math
SF(x)S^{-1}
 =\bigoplus_{[\sigma]}A^{r_\sigma(x)},\qquad
r_\sigma(x)=\sum_i\sigma_i x_i\pmod3.                   \qquad(2)
```

There are q pairs; choose, for example, representatives with first sign +1.
The residues in (2) are signed integer sums modulo three, not binary parities.
Define the coherent batch on arbitrary block inputs by

```math
\mathcal B_p|x,z\rangle
 =|x,(\bigoplus_{[\sigma]}A^{r_\sigma(x)})z\rangle.
```

Then the complete moved-zeta operator is

```math
\mathcal M_p=U_{S^{-1}}\mathcal B_p U_S U_M.             \qquad(3)
```

## 2. Constructive binary basis and its native cost

The descent supplies an explicit compiler for S. Choose a tensor-product
eigenvector v_sigma over $`\mathbb F_4`$, and pair it with its
coordinatewise square. The binary vectors

```math
\phi_\sigma(t)=t v_\sigma+t^2v_\sigma^2,
\qquad t\in\mathbb F_4,
```

are Frobenius-fixed. The two columns phi_sigma(1), phi_sigma(omega)
form a binary basis for this pair, because the two eigencharacters are
distinct. Moreover $`A_i\phi_\sigma(t)=\phi_\sigma(t\omega^{\sigma_i})`$.
Collect these columns over all pairs and invert the resulting binary
matrix to obtain S. Thus the basis is fixed independently of x and s.

Gaussian elimination emits S and its actual inverse with $`O(m^2)`$
CNOTs, zero T gates and no helpers. Binary row swaps cost three CNOTs;
reversing the emitted word pays the actual inverse with the same count.
This is a basis change on m physical selector wires. It is distinct from
an unpriced complex eigenbasis of a quantum three-cycle on two qubits.

Put $`v_0=Se_0`$. Every two-bit block component of v_0 is nonzero.
Indeed $`F(x)e_0=e_x`$ for binary addresses x, so its orbit already
spans the entire binary selector space. A zero block component would
contradict this span after applying S. Within each block an A power can
send that nonzero component to (1,0), without changing (2). Hence S may
be chosen so that

```math
v_0=\bigoplus_{[\sigma]}(1,0),\qquad
P_0=X^{v_0}=\prod_{[\sigma]}X_{\sigma,0}.                \qquad(4)
```

This normalization only changes the fixed binary Clifford walls.

## 3. A paid linear-depth native baseline and optional signal control

A block A sends $`(s,t)\mapsto(t,s\oplus t)`$ and is emitted
chronologically as CX(s,t), CX(t,s). Reversing these gates gives its
inverse. For each address bit x_i, apply the appropriate two CNOT
orientations, controlled by x_i, on every disjoint block. Each of the
two shared-root CCX batches has four physical T layers.

Let $`\chi=q\bmod2`$, $`k=\lfloor q/2\rfloor\bmod4`$ and
$`\bar k=-\lfloor q/2\rfloor\bmod4`$. A shared-root CCX batch
has T-count $`6q+\chi`$ and forward/inverse Clifford counts
$`12q+k`$, $`12q+\bar k`$, respectively. These count positive S
powers after reduction modulo four. Consequently a common ledger is

```math
T_{\mathcal B}\le p(12q+2\chi),\qquad
G_{\mathcal B}\le p(24q+6),\qquad
D_T(\mathcal B_p)\le8p.                                \qquad(5)
```

No clean or dirty helpers are used. Logical commutation of address
stages does not remove their shared physical target wires.

Full signal control is also paid, when that stronger interface is wanted.
Let $`\mathcal R=\prod_\sigma\mathrm{SWAP}(z_{\sigma,0},z_{\sigma,1})`$.
Then $`\mathcal R\mathcal B_p\mathcal R=\mathcal B_p^{-1}`$.
The chronological word

```text
B_p^-1, controlled-R, B_p, controlled-R
```

is identity at signal zero. At signal one its operator is
$`\mathcal R\mathcal B_p\mathcal R\mathcal B_p^{-1}=\mathcal B_p^{-2}=\mathcal B_p`$.
A shared-root controlled-R
Fredkin batch costs $`6q+\chi`$ T gates and four T layers.
The resulting common forward/actual-inverse ledger is

```math
\begin{aligned}
T_{\mathrm c\mathcal B}&\le(2p+1)(12q+2\chi),\\
G_{\mathrm c\mathcal B}&\le48pq+28q+8p+6,\\
D_T(\mathrm c\mathcal B)&\le16p+8.
\end{aligned}                                                \qquad(6)
```

For the Clifford bound, the two batch orientations contribute
$`48pq+2p(k+\bar k)`$ and the two Fredkins contribute at most
$`28q+6`$; use $`k+\bar k\in\{0,4\}`$. There are again no
helpers or clean assumptions. These exact words have literal phase +1.
Equations (5)--(6) are baselines, not sublogarithmic compilers.

<a id="pauli-transport-interface"></a>

## 4. The indicator needs one observable, not a whole controlled batch

Let $`F=X^{\otimes m}`$. The fixed zeta sends the all-ones binary
vector to e_0, so $`U_M F U_M^\dagger=X_{s_0}`$. Therefore

```math
\mathcal M_pF\mathcal M_p^\dagger
 =X_{s_x}=I_p,\qquad
I_p|x,s\rangle=|x,s\oplus e_x\rangle.                    \qquad(7)
```

This is the complete dirty indicator on the existing arbitrary selector
bank. No addressed flip is inserted for free, and the eight-stage
fixed-row central word is unnecessary for this reduction.

Suppose a zero-clean emitted E obeys the strong full-operator promise
$`\|E-\mathcal M_p\otimes I_W\|\le\epsilon`$ on arbitrary
private work W. Then $`EFE^\dagger`$ has indicator error at most
$`2\epsilon`$, using the actual emitted inverse. Its common bounds are

```math
T\le2T_E,\qquad G\le2G_E+4m,\qquad D_T\le2d_E.          \qquad(8)
```

The 4m Cliffords emit F using $`X=HS^2H`$; no further work is added.
Bank-Hadamard walls turn I_p into the lookup phase $`(-1)^{Y_x}`$.

Full approximation of moved zeta is sufficient but unnecessary. With
$`C=(\mathcal M_p^\dagger\otimes I_W)E`$, exactly

```math
\|EFE^\dagger-I_p\otimes I_W\|
 =\|CF-FC\|=\|CFC^\dagger F-I\|.                        \qquad(9)
```

Here CF-FC is the additive commutator, not the group commutator itself.
Thus arbitrary relative errors commuting with F disappear. A global
phase of E also cancels; a sign error in the transported observable
does not. The orientation $`EFE^\dagger`$ matters for this criterion.

In block coordinates the same exact identity becomes

```math
I_p=U_{S^{-1}}\mathcal B_pP_0\mathcal B_p^\dagger U_S.   \qquad(10)
```

The relevant mask is $`Se_0`$, not S times the all-ones vector.
It suffices to emit a zero-clean unitary E_B satisfying the weaker
complete observable promise

```math
\|E_BP_0E_B^\dagger
 -(\mathcal B_pP_0\mathcal B_p^\dagger)\otimes I_W\|
 \le\tau.                                                \qquad(11)
```

The final echo returns all work in its ideal comparison. E_B itself
need not return it or preserve the address. Every input column and
reference remains covered by (11); a prepared-state or clean-output
promise is insufficient. A full batch approximation of error epsilon
implies (11) with tau at most $`2\epsilon`$, but is stronger.
After the exact S walls the indicator has the same tau, depth at most
$`2d_{E_B}`$, T-count at most $`2T_{E_B}`$ and Clifford count at most
$`2G_{E_B}+2G_S+4q`$ using (4).

Even controlled indicator does not require controlled E_B. Replace P_0
in (10)--(11) by its signal-controlled form, q ordinary CNOTs from a
separate signal to the fixed mask. Keep E_B and both S walls unconditional.
At signal zero the actual E_B/E_B-adjoint cancel exactly. At signal one
the same observable error applies, so the complete controlled error is
still tau. This adds no T layer and uses q rather than 4q Pauli Cliffords.
If another source genuinely requires the full controlled batch or
controlled moved zeta, that is a separate obligation; (6) and the paid
controlled fixed zeta give only their linear-depth baseline.

Transporting the specified X on each two-bit block is sufficient for
(11), but the joint observable promise is weaker still. It does not
require independent blocks or an unpriced complex Fourier basis.

<a id="paid-indicator-converse"></a>

## 5. Paid converses under the quadratic Clifford budget

Conversely, a complete dirty indicator implements the full batch with
constant query overhead. The factorwise-return construction below uses
twelve calls; Section 5.1 reduces this to four by sharing cleanup, and
needs only two calls for the weaker observable source.
These constructions use arbitrary auxiliary banks, not clean
predicate bits. For any table $`f(x)\in\{0,1\}^q`$, let A_f have
column x equal to f(x), and define the Clifford map

```math
C_f:(T,b)\mapsto(T,b\oplus A_fT).
```

It costs at most qm CNOTs on an arbitrary m-bit T and q-bit b. The word
$`E_f=I_p^\dagger C_f I_p C_f`$ toggles b by f(x) and restores T.
For disjoint two-wire CNOT involutions G_j, put
$`B_G=\prod_j\mathrm c_{b_j}G_j`$. Then

```math
E_f^\dagger B_GE_fB_G=\prod_jG_j^{f_j(x)}.                \qquad(12)
```

Both b and T return. This uses four actual indicator/adjoint calls,
four C_f maps and two disjoint CCX batches. The latter cost 14q T gates,
24q Cliffords and eight T layers in total.

For each residue r_sigma let a=[r_sigma=1], c=[r_sigma=2] and
$`b=a\oplus c`$. The binary identity
$`A^{r_\sigma}=U^aL^bU^c`$ uses three such Boolean-controlled batches,
chronologically $`U^c,L^b,U^a`$. A complete indicator of error delta
therefore gives the batch with common bounds

```math
\begin{aligned}
T&\le12T_I+42q,&G&\le12G_I+12qm+72q,\\
D_T&\le12d_I+24,&w&\le m+q+w_I,\qquad\epsilon\le12\delta.
\end{aligned}                                                \qquad(13)
```

The twelve-call error follows by full-operator telescoping, including
actual adjoints; private implementation work is covered by the indicator
promise. No controlled-indicator oracle is used. Thus complete indicator,
complete batch and complete moved zeta have constant-factor reductions
within the allowed quadratic Clifford budget. This converse would not
preserve a stricter near-linear Clifford-count requirement.

Substituting the [retained exact dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md)
into (13) already gives an exact full batch, and hence moved zeta, with
$`T,w=O(m)`$, $`G=O(m^2)`$ and $`D_T=O(\log(p+2))`$.
These circuits use arbitrary returned dirty helpers and zero clean inputs.
The two-call control construction from Section 3 gives the same depth
order for the controlled batch. This improves the standalone helper-free
linear-depth baseline, but inherits its logarithmic depth from the existing
indicator; it is not an improvement to the general frame frontier.

<a id="shared-cleanup-transport"></a>

### 5.1. Two-call observable transport and four-call full batches

The observable promise permits a further saving: intermediate work need
not return after each table lookup or even after the complete source.
Keep the same arbitrary query bank T and predicate bank b. Define

```math
R_f=C_f I_p C_f,\qquad V_{f,G}=B_G R_f B_G.
```

On every computational column, including all arbitrary T,b,z,

```math
\begin{aligned}
R_f:(x,T,b)&\longmapsto(x,T\oplus e_x,b\oplus f(x)),\\
V_{f,G}:(x,T,b,z)&\longmapsto
 (x,T\oplus e_x,b\oplus f(x),\prod_jG_j^{f_j(x)}z).
\end{aligned}                                                \qquad(14)
```

The two C_f maps cancel their dependence on the initial T. The two
B_G maps apply exponents b_j and $`b_j\oplus f_j(x)`$ to the same
involution, so their product is independent of the initial b. All maps
are literal phase-free permutations; hence (14) is a full-operator
identity on coherent inputs and external references. No auxiliary
bank was initialized. One V uses one indicator call, two C_f maps and
two disjoint CCX batches.

Write $`a_\sigma=[r_\sigma=1]`$,
$`\beta_\sigma=[r_\sigma\ne0]`$ and
$`c_\sigma=[r_\sigma=2]`$, so $`\beta=a\oplus c`$.
For the first coordinate vector of one block,

```math
U^{a_\sigma}L^{\beta_\sigma}e_0=A^{r_\sigma}e_0.
```

At residues 0,1,2 the two-shear maps are I,A,L, respectively. The
residue-two map differs from $`A^2`$ by a right factor U, which fixes
e_0. Thus this is an actual change of representative in the allowed
Pauli centralizer, rather than an approximation of the full batch.

Apply $`V_{\beta,L}`$ followed by $`V_{a,U}`$, reusing T and b.
Equation (14) returns T after the two toggles, maps the principal bank
by $`\prod_\sigma U_\sigma^{a_\sigma}L_\sigma^{\beta_\sigma}`$,
and changes b by c(x). Call this full unitary E_*. The residual
address-controlled $`X_b^{c(x)}`$ commutes with the transported bank
observable, so

```math
E_*P_0E_*^\dagger=\mathcal B_pP_0\mathcal B_p^\dagger
 \otimes I_{T,b,W}.                                         \qquad(15)
```

The source itself generally changes b. Its complete Pauli echo restores
b, T and the indicator's private work W on every input. Replacing each
ideal indicator call by an emitted approximation of error delta gives
$`\|E-E_*\|\le2\delta`$ and observable error at most $`4\delta`$.
With common forward/actual-adjoint indicator bounds, the source ledger is

```math
\begin{aligned}
T_E&\le2T_I+28q,&G_E&\le2G_I+4qm+48q,\\
d_E&\le2d_I+16,&w_E&\le m+q+w_I.
\end{aligned}                                                \qquad(16)
```

The resulting indicator echo has four actual indicator/adjoint calls.
Use its actual emitted adjoint: no approximate indicator is assumed
involutive, and no intermediate approximate work return is assumed.
Full-operator telescoping already includes every return error.

There is also a saving for the stronger full batch. Apply
$`V_{c,U},V_{\beta,L},V_{a,U}`$ chronologically, then one bare
indicator on T. The three predicate shifts sum to zero. The fourth
query-bank toggle restores T. The principal action is exactly
$`U^aL^\beta U^c=A^r`$ on every block, with all private work returned.
The improved common full-batch ledger is

```math
\begin{aligned}
T_{\mathcal B}&\le4T_I+42q,&
G_{\mathcal B}&\le4G_I+6qm+72q,\\
d_{\mathcal B}&\le4d_I+24,&
w_{\mathcal B}&\le m+q+w_I,\qquad\epsilon_{\mathcal B}\le4\delta.
\end{aligned}                                                \qquad(17)
```

The fourth call may be an actual adjoint instead, since the ideal
indicator is self-adjoint. Reversing the emitted construction gives
the separately charged actual inverse with the same common bounds.
Both improvements retain zero clean work and the quadratic Clifford
budget. Substituting the existing exact indicator still gives
$`O(\log(p+2))`$ depth: these are constant-factor improvements, not
sublogarithmic constructions or a new general frame bound.

<a id="two-band-transport-screen"></a>

### 5.2. Two native T bands cannot supply the weak observable

For $`p\ge3`$, any zero-clean source E with at most two physical
T/T-dagger layers has observable error at least

```math
\tau_* = \frac{\sin(\pi/8)}{2\sqrt2}\approx0.135299          \qquad(18)
```

in (11), irrespective of its width, Clifford count or intermediate
address and work movement. This follows from the retained
[two-layer diagonal rigidity theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction),
after folding the Pauli echo at its center.

Indeed, write $`E=C_2D_2C_1D_1C_0`$ with physical T bands D_i
and arbitrary Clifford blocks C_i. A one-band encoder sends a Hermitian
Pauli to a literal Hermitian Clifford, so

```math
K=C_1D_1(C_0P_0C_0^\dagger)D_1^\dagger C_1^\dagger
```

is Clifford, including its actual phase. Consequently

```math
EP_0E^\dagger=C_2D_2K D_2^\dagger C_2^\dagger             \qquad(19)
```

has two physical T bands, not four. This folding is also proved in
[R75's two-layer Pauli-image form](FOUR_LAYER_OPERATOR_RIGIDITY.md#1-two-layer-pauli-images-and-commuting-axis-spaces).
It remains valid when P_0 is tensored with identity on every arbitrary
helper; E need not preserve any individual register.

The fixed S walls and output Hadamards turn the ideal observable into
the complete diagonal lookup sign $`(-1)^{Y_x}`$, without changing
T-depth or norm error. Its Boolean exponent has degree p+1: in

```math
Y_x=\bigoplus_u Y_u\prod_i(1\oplus x_i\oplus u_i)
```

each distinct Y_u has a nonzero full-address monomial. For completeness,
the zero-clean proof of the cited theorem applies to this multivariate
sign directly. Error below tau_* forces the same actual two-band word
to preserve every computational Z exactly by R75 operator rounding.
It is then diagonal C3. The folded observable is also Hermitian, so
its diagonal entries are literal Boolean signs, of degree at most three.
It cannot equal the target of degree p+1 greater than three. Distinct
Boolean diagonal unitaries have distance two, contradicting error below
tau_*. This argument preserves the observable's sign; a source global
phase cancels from its actual-adjoint echo.

Thus at least three source T layers are necessary at this fixed error
threshold. More generally a d-layer source echo folds to at most
$`2d-2`$ T layers for d at least one, but (18) does not exclude three
source layers or growing sublogarithmic depth. It is a native two-band
screen, distinct from the classical multiplication-degree bound below.

<a id="three-band-transport-exactness"></a>

### 5.3. Three-layer observable approximation must already be exact

The same operator rounding gives a sharper discovery condition at three
source layers. Let E be any zero-clean native circuit of T-depth at most
three, P a fixed Hermitian Pauli on its whole physical space, and D a
Boolean diagonal involution in any fixed Clifford basis. Then

```math
\|EPE^\dagger-D\|\lt\tau_*
\quad\Longrightarrow\quad EPE^\dagger=D.                    \qquad(20)
```

There is no restriction on dirty width, Clifford count, or intermediate
address and work movement. The premise is a full-operator comparison;
an initialized-sector promise is not sufficient for this argument.

To prove (20), fold the central source band as in (19). The literal
Hermitian unitary $`O=EPE^\dagger`$ has at most four T layers. Move to
the fixed Clifford basis in which D is diagonal. For every computational
Z, including every dirty-work Z, its distance from D gives

```math
\|OZO^\dagger-Z\|\le2\|O-D\|\lt2\tau_*.
```

[R75's four-layer preservation theorem](FOUR_LAYER_OPERATOR_RIGIDITY.md#5-four-layer-preservation-and-the-zero-clean-witness)
therefore makes O commute with every such Z exactly, simultaneously and
without a union bound. Thus O is diagonal. Since it is a Hermitian
involution, its entries are literal signs. Distinct Boolean diagonal
involutions have distance two, so the assumed error forces O=D.
The source's global phase cancels in its actual-adjoint echo; no target
sign or input-dependent phase is discarded.

The fixed S walls and output Hadamards put the target of (11) into this
form. Hence a source of depth at most three either meets the complete
MOD3 observable contract exactly or has error at least tau_*. Equivalently,
below that threshold its relative word
$`C=(\mathcal B_p^\dagger\otimes I_W)E`$ must commute exactly with
$`P_0\otimes I_W`$. For p at least three, Section 5.2 excludes two
source layers, while exact three-layer transport remains open. This
rounding statement neither excludes that exact construction nor bounds
growing sublogarithmic source depth. The retained constant-depth
approximate scalar zero reflection is compatible with this fixed-layer
threshold; conjunction still supplies no dimension-growing witness.

<a id="scoped-mod3-limits"></a>

## 6. Classical degree and approximation boundaries

The block matrix has the binary polynomial expansion

```math
A^{r_\sigma(x)}=\prod_i\bigl[I+(A^{\sigma_i}+I)x_i\bigr].
```

Since $`A^{\sigma_i}+I=A^{-\sigma_i}`$, its full-address coefficient
is the invertible matrix $`A^{-\sum_i\sigma_i}`$. Applied to two
arbitrary block bits, every output therefore has algebraic degree p+1.
Signal control gives degree p+2, since its matrix is
$`I+a(A^{r_\sigma(x)}-I)`$. A classical binary circuit whose only
nonlinear gates multiply two bits has multiplicative depth at least
$`\lceil\log_2(p+1)\rceil`$, or $`\lceil\log_2(p+2)\rceil`$
with signal control: each multiplication layer can at most double degree.
Products within a layer consume pre-layer values; affine interlayers are
unrestricted. Arbitrarily many clean or dirty helpers do not change this
degree argument, and it does not require their return.
This is not a native T-depth lower bound for quantum Clifford+T circuits.

One proposed dirty fanout shortcut is also excluded. Put
$`\zeta=e^{2\pi i/3}\in\mathbb C`$, distinct from the finite-field omega.
As a two-qubit quantum permutation, A has eigenvalue multiplicities
(2,1,1) for $`(1,\zeta,\zeta^2)`$. Thus $`A\otimes I_4`$ has multiplicities
(8,4,4), whereas $`A\otimes A`$ has (6,5,5). For every unitary C,
$`\|C^\dagger(A\otimes I_4)C-A\otimes A\|\ge\sqrt3`$:
the former's eight-dimensional fixed space intersects the latter's
ten-dimensional nonfixed space, where distance from identity is sqrt(3).
Extra identity helpers preserve this argument. It excludes this
conjugation-only fanout, not all dirty MOD3 arithmetic or a different
final helper action.

Independently approximating block eigenbases or phases does not pay the
complete interface automatically. A sum of block operator errors is a
sufficient global error bound; without an additional argument it is not
a necessary per-block precision bound. Joint synthesis may correlate or
cancel errors. All claims must refer to actual operators, literal phases,
arbitrary inputs and any returned implementation work.

There is an exact product-observable statement under a stated angle bound.
At a fixed address let P_j,Q_j be ideal and actual Hermitian unitary
observables, with $`\theta_j=\max_{\lambda\in\mathrm{spec}(P_jQ_j)}|\arg\lambda|`$.
Conjugation by P_j sends $`P_jQ_j`$ to its adjoint, so its extreme angles
occur in plus/minus pairs. If $`\Theta=\sum_j\theta_j\le\pi`$, then

```math
\left\|\bigotimes_jQ_j-\bigotimes_jP_j\right\|
 =2\sin(\Theta/2).
```

The product spectrum contains the sum of positive extreme angles, and
every other sum has magnitude at most Theta. This concerns independent
product observables in that angle regime, not arbitrary correlated errors,
full-gate approximation, or a native T-depth lower bound.

<a id="frontier-transfer"></a>

## 7. What a successful compiler would buy

At every fixed requested tolerance, a zero-clean E_B satisfying (11)
with $`T,w=O_\tau(m\,\mathrm{poly}(p))`$,
$`G=O_\tau(m^2\,\mathrm{poly}(p))`$ and
$`d_{E_B}=o(\log(p+2))`$ would provide the relaxed complete indicator.
The fixed S walls fit its Clifford budget. Polynomial-overhead fusion
then supplies the normalized counts, with its paid logarithmic-logarithmic
tail; the retained frame transfer carries that benefit to full frames.
At fixed precision and square-root dirty width this would improve the
current logarithmic depth contribution, but even constant source depth
retains the reductions' log-log term. No growing-precision improvement
follows from fixed-tolerance premises alone.

The first unpaid operation is sublogarithmic complete batch synthesis or
its weaker Pauli-transport interface. Neither the direct baseline nor
the indicator-based logarithmic compiler meets it, and no
balanced dirty MOD3 arithmetic is assumed available. The unconditional
R96 frontier remains unchanged.

At three source layers, (20) turns the small-fixed-error question into
an exact joint transport question. A source with more layers may instead
exploit approximation, provided its complete depth remains sublogarithmic.
The [joint-coordinate construction](WEAK_PAULI_ORBIT_COMPRESSION.md#two-coordinate-lift)
pays two coordinates at every bank size, with twelve T gates, three layers
and no extra work. It supplies the complete source at p=2 only. Retaining
its final-block commutation condition cannot extend the construction to
p at least three. The [coupled-final-block theorem](COUPLED_FINAL_ENCODER.md#bank-independent-final-boundary)
also excludes bank-independent final one-band blocks at those lengths,
despite arbitrary middle Clifford boundaries and dirty width. A separate
[full-bank X-motion word](COUPLED_FINAL_ENCODER.md#full-bank-x-motion-boundary)
also has no one-band middle completion, by its rotated-bank centralizer and
affine response code. Other bank-dependent final blocks or a sublogarithmic
source with more layers remain open. The
[four-band nested-coefficient source](WEAK_PAULI_ORBIT_COMPRESSION.md#four-band-nested-coefficient-source)
remains partial. A separate [two-band coefficient shear](WEAK_PAULI_ORBIT_COMPRESSION.md#complete-coefficient-shear)
now completes every three-address coefficient in five source layers, with
28 T gates and no helpers. Its exact all-size recurrence pays the complete
bank but retains linear emitted address depth. Coefficient coverage is
settled; a sublogarithmic joint native schedule remains unpaid. Neither
this completion nor the scoped exclusions move R96.

Exact fixtures: [binary transport and observable identities](../../tests/test_mod3_zeta_transport.py),
[native batch, actual inverses and signal controls](../../tests/test_mod3_native_batch.py),
and [shared cleanup and Pauli-echo folding](../../tests/test_mod3_observable_transport.py).
Finite fixtures corroborate these reductions and baselines; they do not
emit a sublogarithmic source or prove an unrestricted depth lower bound.
