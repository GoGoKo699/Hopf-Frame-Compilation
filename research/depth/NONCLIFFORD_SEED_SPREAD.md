# Pauli spread required of a non-Clifford seed

[Interleaved tree target](INTERLEAVED_TREE_PROBE_SOURCE.md) ·
[Clifford-orbit query rank](CLIFFORD_ORBIT_QUERY_RANK.md) ·
[Fixed-readout syndrome bound](FIXED_READOUT_SYNDROME_BOUND.md)

A genuinely non-Clifford starting reflection leaves the Clifford-orbit
rank theorem's premise, but a small Pauli-support space remains insufficient.
This note quantifies the required spread even with an arbitrary fixed
decoder, arbitrary coefficient degeneracies and arbitrary dirty width.
The bound concerns a specified queried-Clifford orbit. It is not a lower
bound on unrestricted native T depth.

Let $`p\ge1`$, $`m=2^p`$, and let $`P_x`$ be the complete pulled-back
tree family. Its m Pauli labels are independent, its commutation matrix
is nonsingular, and the $`m/2`$ addresses with final bit one form an
independent commuting family. These properties are proved in the
[target algebra](INTERLEAVED_TREE_PROBE_SOURCE.md#full-rank-probe-algebra).

## 1. The orbit fits in a fixed Pauli algebra

On w complete work wires, expand a fixed Hermitian involution R in the
Hermitian Pauli section:

```math
R=\sum_v r_v\tau(v),\qquad r_v\in\mathbb R,\qquad
\sum_v r_v^2=1,\qquad
L_0=\mathrm{span}_{\mathbb F_2}\{v:r_v\ne0\},\qquad r_0=\dim L_0.
```

Here $`r_0`$ is the **linear label-span dimension**, not the number of
support terms, their affine dimension, or their physical Pauli weight.
The identity label, if present, contributes no dimension.

Let $`C_x`$ be a complete Clifford family whose symplectic matrix
$`S(x)`$ has entrywise Boolean degree at most d. Define

```math
R_x=C_xRC_x^\dagger,\qquad
M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj.
```

For a basis $`v_1,\ldots,v_{r_0}`$ of $`L_0`$, write the vector-valued
algebraic normal forms

```math
S(x)v_j=\sum_{|A|\le d}x^A v_{j,A},\qquad
W=\mathrm{span}_{\mathbb F_2}\{v_{j,A}:j,A\},\qquad
r:=\dim W\le r_0M_d(p).                              \qquad(1)
```

Every support label of every $`R_x`$ belongs to W. Thus every orbit
member lies in the same unital Pauli algebra

```math
\mathcal A(W)=\mathrm{span}_{\mathbb C}\{\tau(v):v\in W\},
\qquad \dim_{\mathbb C}\mathcal A(W)=2^r.              \qquad(2)
```

This is an algebra because Pauli products add their labels. Pauli signs
and their address dependence do not affect (1), and no independence of
the real coefficients is assumed.

For q literal affine XOR queries interleaved with fixed Cliffords,
each query's symplectic entries are affine, so matrix multiplication
certifies $`d\le q`$. Unlike a full phase-lift argument, this label-span
bound still permits arbitrary Clifford Pauli signs. More general stages
require their own label-degree certificate.

## 2. Exact decoding requires a full quantum algebra

Write $`2s`$ for the rank of the symplectic form restricted to W and
$`t=r-2s`$ for its radical dimension. Fixed Clifford coordinates put
the algebra in the form

```math
\mathcal A(W)\cong
 M_{2^s}(\mathbb C)\otimes\mathcal D_{2^t}
 \otimes I_{2^{w-s-t}},                              \qquad(3)
```

where $`\mathcal D_{2^t}`$ is the diagonal algebra. In particular, W
has s quantum Pauli pairs and t central Z directions.

Suppose one fixed unitary E satisfies $`ER_xE^\dagger=P_x`$ for every
address. The target family generates the simple matrix algebra
$`M_{2^{m/2}}(\mathbb C)`$ with spectator multiplicity. Its conjugate
therefore embeds unitally in $`\mathcal A(W)`$. Restrict this embedding
to any central block of (3). It remains a unital representation of a
simple matrix algebra, so is injective; its defining dimension cannot
exceed $`2^s`$. Consequently

```math
s\ge m/2,\qquad r\ge m,\qquad r_0M_d(p)\ge2^p.       \qquad(4)
```

The central blocks and arbitrary spectator multiplicities cannot replace
the quantum dimension. The weaker dimension comparison
$`2^m\le\dim\mathcal A(W)=2^r`$ also proves the last two inequalities.

## 3. A robust bound from the independent commuting half

Let E be an otherwise arbitrary fixed decoder on the complete work
space, and suppose

```math
\max_x\|ER_xE^\dagger-P_x\|\le\delta,
\qquad 0\le\delta\lt\sqrt2,
\qquad \alpha(\delta)=1-H_2(\delta^2/4).               \qquad(5)
```

The original target half $`x_p=1`$ contains $`k=m/2`$ commuting,
independent Hermitian Paulis. Every joint sign string
$`z\in\{1,-1\}^k`$ therefore occurs: its joint projector has rank
$`2^{w-k}`$. Choose the normalized joint projector $`\rho_z`$ and
pull it back to $`E^\dagger\rho_zE`$.

In (3)'s coordinates, dephase the t central wires and trace out the
spectators. This channel preserves the expectation of every observable
in $`\mathcal A(W)`$, so preserves every $`R_x`$ measurement. Its
output is a state $`\sigma_z`$ on s quantum wires and t classical
binary coordinates, represented on a Hilbert space of dimension

```math
D=2^{s+t}=2^{r-s}.                                  \qquad(6)
```

For a target sign $`z_i`$, every unit vector in the corresponding joint
eigenspace obeys $`\|(ER_iE^\dagger-z_i I)\psi\|\le\delta`$.
Both compared observables are involutions, giving the sharper squared
error relation

```math
\frac{1-z_i\mathrm{Tr}(R_iE^\dagger\rho_zE)}2
 \le\frac{\delta^2}{4}.                             \qquad(7)
```

The left side is exactly the wrong-outcome probability of measuring
$`R_i`$. Compression preserves this probability. The states
$`\sigma_z`$ therefore form a quantum random-access encoding of k
independent classical bits, each recoverable with error at most
$`\delta^2/4`$.

The information bound used here is
[Nayak's Theorem 2.3](https://arxiv.org/abs/quant-ph/9904093): an encoding
on a D-dimensional Hilbert space with bit error at most e has
$`\log_2D\ge k[1-H_2(e)]`$. For completeness, its entropy argument
averages the encoding states one input bit at a time. For each binary
split, the decoding measurement and the binary Fano bound give a Holevo
information of at least $`1-H_2(e)`$. Summing these entropy increases
over k levels gives $`k[1-H_2(e)]`$; the final state's entropy is at most
$`\log_2D`$. Mixed encodings and arbitrary decoding measurements are
allowed.

Combining (1), (6) and (7) proves

```math
r_0M_d(p)\ge r\ge r-s\ge\frac{2^p}{2}\alpha(\delta).
                                                               \qquad(8)
```

The compression is a mathematical information argument. It is not a
native implementation, an initialized-work promise, or an available
source subroutine. The test states are legitimate because (5) holds on
every work input. Tensoring an arbitrary reference leaves the norm
contract unchanged.

In particular, q literal affine queries obey

```math
r_0M_q(p)\ge2^{p-1}\alpha(\delta),\qquad
q\ge\frac p2-
 \sqrt{\frac p2\ln\!\frac{2r_0}{\alpha(\delta)}}.      \qquad(9)
```

For the second inequality, if $`q\lt p/2`$, apply the binomial lower-tail
bound to $`M_q(p)/2^p`$; if $`q\ge p/2`$, it is automatic. A bounded
or polynomial-size label span therefore retains linear address-query
count at fixed error. Conversely, a fixed q requires
$`r_0=\Omega_\delta(2^p/p^q)`$. Growing physical Pauli weight by fixed
Clifford spreading does not help, because it preserves $`r_0`$.

## 4. The seed itself must supply the spread

Suppose the actual seed is $`R=V P V^\dagger`$, where P is one fixed
Pauli and V contains $`T_{\rm seed}`$ native T or T-dagger gates.
A Clifford passage preserves affine Pauli-support dimension. A single
T-axis passage keeps a label or adds that axis to it, so increases the
affine support dimension by at most one. Starting from one label gives

```math
r_0\le T_{\rm seed}+1,\qquad
(T_{\rm seed}+1)M_q(p)\ge2^{p-1}\alpha(\delta).        \qquad(10)
```

Thus a localized nonlinear gadget cannot provide the required spread
for bounded q merely by surrounding it with fixed Clifford walls and
an elaborate fixed decoder. This is a seed-count and query-grammar
constraint. An extensive seed with linear T count can still have
constant native depth and remains compatible with the discovery budget.

### An extensive seed with a paid two-band source

The following seed makes that surviving boundary concrete. On w active
wires, put

```math
A=THT^\dagger=\frac Z{\sqrt2}+\frac{X+Y}{2},\qquad
R=A^{\otimes w}\otimes I_h.
```

Its support has $`3^w`$ distinct terms, linear and affine label-span
dimension $`2w`$, and largest squared coefficient $`2^{-w}`$.
Thus w proportional to m supplies extensive spread and avoids the
dominant-coefficient premise below.

This is an actual native component, not just a prescribed reflection.
Let $`V=SHTHS^\dagger`$, so $`VZV^\dagger=H`$, and define

```math
L=(TV)^{\otimes w}\otimes I_h,\qquad
U_0=L\prod_{j=1}^{w-1}\mathrm{CX}(j;0),\qquad
U_0Z_0U_0^\dagger=R.
```

The CNOT fan-in first maps $`Z_0`$ to the product of all active Z
operators. The chronological local walls are $`S^\dagger,H,T,H,S,T`$,
with each stage applied in parallel across the w active wires. Writing
each S-dagger as three native S gates gives $`9w-1`$ total gates,
$`2w`$ T gates and two T bands. The actual inverse has the same counts
with T-dagger gates. Both words act as the identity on all h helpers;
they need no initialized work. Paying both words is required wherever
the complete source uses both directions.

### The plain local decoder has a constant error gap

The evident fixed inverse $`L^\dagger`$ has uniform decoding error at
least $`1/2`$ for the full tree family when $`p\ge2`$, even with
arbitrary fixed Clifford wrappers. More precisely, allow any Clifford
family $`C_x`$ and decoder $`E=F_{\rm out}L^\dagger F_{\rm in}`$,
where both F operators are fixed Cliffords on all active and helper
wires. Then

```math
\max_x\|EC_xRC_x^\dagger E^\dagger-P_x\|\ge\frac12.
```

There is no restriction on the number of queries in this statement.
Absorb $`F_{\rm in}`$ into $`C_x`$ and put
$`Q_x=F_{\rm out}^\dagger P_xF_{\rm out}`$. An error smaller than
$`1/2`$ would put $`S_x=LQ_xL^\dagger`$ within that distance of
$`R_x=C_xRC_x^\dagger`$ for every x.

For the single-cell unitary $`U=TV`$, direct conjugation gives

```math
UZU^\dagger=A,\qquad UXU^\dagger=-ZAZ,\qquad
UYU^\dagger=\frac{Y-X}{\sqrt2}.
```

The necessary rigidity comes from signed Pauli symmetries. Single-site
normalized traces are

```math
\mathrm{tr}(AXAX)=\mathrm{tr}(AYAY)=-\frac12,\qquad
\mathrm{tr}(AZAZ)=0,\qquad \mathrm{tr}(AIAI)=1.
```

In particular, $`AXA`$ is not a Pauli, since its Pauli coefficient
along X is $`-1/2`$. Thus R is genuinely non-Clifford for every
$`w\ge1`$. If a Pauli B acts nontrivially on any active wire,
the tensor product trace obeys $`|\mathrm{tr}(RBRB)|\le1/2`$.
The normalized Hilbert--Schmidt norm is a lower bound for operator
norm, so for either sign $`s\in\{1,-1\}`$,

```math
\|R-sBRB\|^2\ge
 \|R-sBRB\|_{2,\rm norm}^2
 =2-2s\,\mathrm{tr}(RBRB)\ge1.                       \qquad(11)
```

R's only exact signed Pauli symmetries are the helper Paulis, all with
positive sign. Their label space has dimension $`2h`$. Any other
signed Pauli relation has defect at least one; a negative relation
using only helpers has defect two. These facts survive Clifford
conjugation.

If $`BS_xB=sS_x`$, approximation would give
$`\|R_x-sBR_xB\|\lt1`$. Consequently every exact signed Pauli
symmetry of $`S_x`$ must be a positive symmetry of $`R_x`$. An active
Y factor in $`Q_x`$ is impossible because its displayed local image
has the negative symmetry Z. A nonidentity helper factor is impossible
because a helper Pauli anticommutes with it. With these factors absent,
each active identity factor gives two independent positive symmetry
directions, in addition to the $`2h`$ helper directions. Their dimension
cannot exceed $`2h`$, so active identity factors are impossible too.
Every $`Q_x`$ must therefore be, up to sign, a full X/Z string on the
active wires and identity on helpers.

Encode its choice of Z factors by $`u_x\in\mathbb F_2^w`$. The
commutation form of two such strings is

```math
\kappa(x,y)=\mathbf1\cdot u_x+\mathbf1\cdot u_y.
```

This matrix has rank at most two, whereas the tree commutation matrix
has rank $`m\ge4`$, a contradiction. An appropriately nonlocal fixed
decoder and a different extensive seed remain open.

At zero error, a simpler cardinality check already forces full active
X/Z support: $`LQL^\dagger`$ has exactly
$`3^{k_{XZ}}2^{k_Y}`$ terms, where $`k_{XZ}`$ and $`k_Y`$ count its
active X/Z and Y factors. Equality with R's $`3^w`$ support forces
$`k_{XZ}=w`$ and $`k_Y=0`$. The symmetry argument is what supplies
the constant error gap and excludes helper factors robustly.

## 5. A dominant Pauli coefficient gives a separate sharp screen

There is a stronger conclusion when one seed coefficient has squared
magnitude $`\mu\gt1/2`$. Absorb its sign into a Hermitian Pauli A and
write

```math
R=\sqrt\mu\, A+R_{\rm rest},\qquad
\sum_{v\ne a}r_v^2=1-\mu,
\qquad A_x=C_xAC_x^\dagger.
```

Suppose (5) holds with $`\delta\lt\mu-1/2`$. For a target pair, let
$`\varepsilon_{xy}\in\{1,-1\}`$ be its required commutation sign.
Four-term telescoping gives

```math
\|R_xR_y-\varepsilon_{xy}R_yR_x\|\le4\delta.          \qquad(12)
```

If the heavy Paulis $`A_x,A_y`$ had the wrong sign, their contribution
to the normal product label $`a_x+a_y`$ in (12) would have magnitude
$`2\mu`$. Every other pair contributing at that label uses a nonheavy
term on both sides: once one term is heavy, the label equation forces
the other to be heavy too. The remaining label equation pairs terms by
a partial bijection. Cauchy--Schwarz therefore bounds the magnitude of
all remaining contributions by $`2(1-\mu)`$, including every literal
Pauli phase. The coefficient, and hence the operator norm, is at least

```math
2\mu-2(1-\mu)=2(2\mu-1)\gt4\delta,
```

a contradiction. Thus the heavy-axis family has exactly the target
commutation matrix. Its nonsingularity forces m independent Pauli
labels. Their degree-d vectors span at most $`M_d(p)`$ dimensions, so

```math
M_d(p)\ge2^p,\qquad d\ge p,
\qquad q\ge p\quad\hbox{for literal affine queries}. \qquad(13)
```

If an actual reflection differs from this ideal seed by at most eta,
the same conclusion holds with $`\delta+\eta\lt\mu-1/2`$.

This coefficient screen does not exclude the intended weighted-cubic
dressing merely because its seed is non-Clifford. Indeed, if
$`R=D_0KD_0^\dagger`$, K is Clifford, $`D_0^\dagger`$ is third-level
and R is not a Pauli up to scalar, then every coefficient is a trace
of a nonscalar Clifford:

```math
r_v=\mathrm{Tr}_{\rm norm}
       \bigl((D_0^\dagger\tau(v)D_0)K\bigr),
\qquad |r_v|\le1/\sqrt2.
```

The bound follows from the retained
[flat Clifford-support fact](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension).
For weighted-cubic phases, the adjoint premise is separately certified
by negating the phase polynomial. The broad-seed question therefore
survives these necessary conditions. Its complete family relations,
actual decoder, inverses and work return remain to be constructed.

The [finite fixtures](../../tests/test_nonclifford_seed_spread.py) retain
native seed phases, local signed-symmetry gaps, actual affine-query label
growth, the central-block compression and squared error, and
heavy-coefficient product collisions.
The all-size conclusions are the analytic arguments above.
