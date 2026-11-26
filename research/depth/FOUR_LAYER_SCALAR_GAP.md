# Scalar and rank gaps at four native T layers

For a physical circuit $`V`$ of at most four native T layers, two Hermitian Paulis $`P,Q`$, and
physical dimension $`D`$, let $`e=(1-D^{-1}\mathrm{Tr}(QVPV^\dagger))/2`$. For data-bit
observables this is the mean disagreement probability on uniformly sampled computational
columns, with no clean initialization. Then

```math
e=0\quad\text{or}\quad e\ge\frac{1-\cos(\pi/12)}8.
```

There is a stronger conclusion when the exceptional columns are known exactly. If the two
observables agree on a subspace of relative dimension greater than $`1-\sin^2(\pi/24)/2`$, they
agree everywhere. In particular, agreement off one $`k`$-bit address forces agreement there as
well whenever $`k\ge7`$. The circuit may move the address and use arbitrary dirty work.

## Sparse Pauli probes

In the
[R75 midpoint frame](FOUR_LAYER_OPERATOR_RIGIDITY.md#1-two-layer-pauli-images-and-commuting-axis-spaces),
write $`A=E_AK_AE_A^\dagger`$ and $`B=E_BK_BE_B^\dagger`$ for the two transported observables.
Here each $`E_j`$ is a Clifford conjugate of one physical T band, and each $`K_j`$ is Hermitian
Clifford. Their full-space normalized Hilbert–Schmidt distance is $`2\sqrt e`$. The binary Pauli
axes of $`E_A`$ are independent because its original T gates act on distinct wires. The map
pairing a Pauli with these axes is onto. Choose a basis of its kernel and a lift of each
coordinate vector of the image. This gives a generating set of Hermitian Pauli probes $`H`$,
each anticommuting with at most one axis.

If $`H`$ meets one axis, $`E_A^\dagger H E_A`$ expands into two distinct signed Paulis with
coefficients $`1/\sqrt2`$. Conjugation by $`K_A`$ preserves Paulis, and conjugation by $`E_A`$
takes each Pauli to a Hermitian Clifford. The normalized Choi vector of $`AHA`$ is therefore
$`(U+V')/\sqrt2`$ for orthogonal unit stabilizer vectors $`U,V'`$ (absorb the signs). With no
anticommuting axis it is one stabilizer vector. On the other side $`BHB=E_B C_HE_B^\dagger`$
with $`C_H`$ Hermitian Clifford. In an $`E_B`$-diagonalizing Clifford basis, its Choi vector
$`C`$ is a stabilizer vector multiplied by the diagonal phase of $`E_B\otimes\overline E_B`$.

Each of these states is flat on its affine support. Ordinary stabilizer phases there are a fixed
eighth root times $`i`$ to a sum of affine bits times a quadratic sign. The eighth-root global
phase follows from Hermiticity by comparing a nonzero matrix entry with its transposed
conjugate, as in
[R75 §2](FOUR_LAYER_OPERATOR_RIGIDITY.md#2-the-one-phase-gap-with-arbitrary-hermitian-clifford-boundaries).
The extra diagonal phase is an eighth root raised to a signed sum of affine parity bits on the
support.

## A gap for the Choi vectors

Put $`b=2\sin(\pi/24)`$. We show that the two Choi vectors coincide or have distance at least
$`b`$. The single-stabilizer case follows directly from R75 §2. In the two-summand case, if
either overlap with $`C`$ has modulus at most $`1/2`$, orthogonality and Bessel's inequality
bound the real overlap of $`C`$ with the sum by $`(1+\sqrt3)/(2\sqrt2)=\cos(\pi/12)`$. This
gives the claimed distance.

Suppose both overlaps exceed $`1/2`$. For intersecting affine supports $`S,T`$ define
$`d(S,T)=\dim S+\dim T-2\dim(S\cap T)`$; set $`d=\infty`$ for disjoint supports. A flat-state
overlap has modulus at most $`|S\cap T|/\sqrt{|S||T|}=2^{-d/2}`$. Relative to
$`S=\mathrm{supp}C`$, each other support is therefore equal to $`S`$, one index-two subspace of
it, or one index-two extension of it.

The following phase estimate applies on every affine support or half-support below. A relative
phase $`\omega^p`$, with $`\omega=e^{i\pi/4}`$, has exponent digits of binary degrees at most
$`1,2,4`$. Indeed, for a sum of affine bits, digit $`j`$ is the elementary symmetric polynomial
of order $`2^j`$ in those bits, by the parity of the corresponding binomial coefficient. Repeat
affine bits to represent signed integer coefficients modulo eight. A quadratic sign contributes
$`4q`$ with $`\deg q\le2`$. The additional exponent $`-2q`$ occurring below contributes
$`2q+4q`$ modulo eight. Its middle digit has degree at most two and the carry from that digit to
the high digit has degree at most four. If the phase is nonconstant or a nonidentity constant,
its lowest nonzero digit takes value one on at least $`1/2,1/4,`$ or $`1/16`$ of the points,
respectively, by the minimum weight of a nonzero binary polynomial. The corresponding real
cosine losses are at least $`1-1/\sqrt2,1,2`$. Thus the mean loss is zero or at least $`1/8`$.

* **A support extends $`S`$.** If the other support lies in $`S`$, the outside contribution to
  the norm is $`1/2`$. Distinct extensions have disjoint outside cosets. If the two extensions
  coincide, both summands on $`S`$ have half the target magnitude. Since both states being
  compared have norm one, their squared distance is the sum of the two mean phase losses on
  $`S`$. This is zero or at least $`1/8`$.
* **A support is half of $`S`$.** Two noncomplementary halves leave at least a quarter of $`S`$
  uncovered. On complementary halves the squared distance is the sum of the two mean phase
  losses. If the other support equals $`S`$, orthogonality on their shared half gives masses
  $`3/4,1/4`$ for the sum on the two blocks. The target masses are $`1/2,1/2`$. The triangle
  inequality on each block gives squared distance at least
  $`(\sqrt3/2-1/\sqrt2)^2+(1/2-1/\sqrt2)^2=b^2`$.
* **All supports equal $`S`$.** The ratio $`w=V'/U`$ is a fixed eighth root times a Clifford
  phase. For odd eighth-root exponents, $`|(1+w)/\sqrt2|`$ differs from one by at least
  $`\sqrt{1+1/\sqrt2}-1>b`$ at every point. Otherwise $`w`$ takes fourth-root values. Its
  real-valued positions are either empty or at least half the support, because the low bit of
  its mod-four exponent is affine. On real positions the magnitude differs from one by at least
  $`\sqrt2-1`$. If there are no real positions, write $`w=i(-1)^q`$ with $`q`$ binary quadratic.
  Then $`(U+V')/\sqrt2=U\omega i^{-q}`$ throughout, and the phase estimate applies.

## Mean error, exceptional rank and complete operators

Every probe has $`\|AHA-BHB\|_{2,\mathrm{tr}}\le2\|A-B\|_{2,\mathrm{tr}}=4\sqrt e`$. Below $`b`$
all probe conjugations agree, so $`BA`$ commutes with a full Pauli basis and is a scalar. Since
$`A,B`$ are Hermitian involutions this gives $`A=\pm B`$, and the minus sign has $`e=1`$. This
proves the mean bound.

For an exact exceptional subspace put $`W=BA`$. Telescoping $`WHW^\dagger-H`$ bounds its rank by
$`2\mathrm{rank}(W-I)`$ and its operator norm by two, so

```math
\|AHA-BHB\|_{2,\mathrm{tr}}^2
\le 8\mathrm{rank}(W-I)/D.
```

When the right side is below $`b^2`$, the same probes force $`W`$ scalar. The fixed subspace
then gives $`W=I`$. For $`k\ge7`$ it is at most $`1/16\lt b^2`$. Positive flip probability may
be arbitrarily small on the exceptional address. Exact silence on every other computational
column means the corresponding Pauli eigenspaces agree on those columns.

For $`p\ge7`$ the complete dirty indicator $`|x,Y,w\rangle\mapsto|x,Y\oplus e_x,w\rangle`$ and
the fully returned $`p`$-controlled NOT have four-layer operator error at least
$`\tau=\sin(\pi/8)/(2\sqrt2)`$, even up to a common phase. Hadamards on the target bank
transform the exact operator into a diagonal sign. Below $`\tau`$, R75 operator rigidity forces
commutation with each computational $`Z`$, hence the implementing circuit is diagonal. Its
unit-modulus diagonal entries are eighth roots in the native phase convention by
[R75 §4](FOUR_LAYER_OPERATOR_RIGIDITY.md#4-discrete-phases-in-the-common-commutant). After
multiplying by the ideal signs, all entries lie in one disk of diameter
$`2\tau\lt 2\sin(\pi/8)`$, hence coincide. The resulting exact oracle contradicts the
single-address rank gap.

An explicit real frame sets its deepest rotation angle to $`\pi`$ at one prefix and all other
angles to zero. Its sign on that prefix acts as identity on the leaf bit. Conjugating a prefix
bit by Hadamards gives a controlled flip, up to a sign on the selected branch, when the other
$`n-2`$ prefix bits match. The diagonal freezing argument followed by the exact single-address
test still applies to that sign. Thus $`n\ge9`$ gives the same zero-clean error floor for a
complete real frame and for its complex realization. The
[canonical frame domains](../../docs/HOPF_INTERFACE.md#canonical-domains) admit these angles and
phases.

<a id="two-band-signed-gap-sources"></a>

## A two-band source cannot give a uniform sparse signed gap

The [full-space signed-gap wrapper](PRESCRIBED_SCALAR_TABLES.md#full-space-signed-gap-source)
can use a correct full-space spectral sign without first approximating a scalar
unitary. The following screen applies to that weaker source. Let V have at most
two physical T/T-dagger bands, let S be a Hermitian Pauli, and put

```math
 W=V^\dagger S V S,\qquad \mathcal H=\tfrac12(W+W^\dagger).
```

All w physical qubits are arbitrary inputs. Suppose
$`|\mathcal H|\succeq\gamma I`$, with $`\gamma\gt0`$, and let
$`0\lt\varepsilon\le1/2`$ be the smaller relative dimension of its positive
and negative eigenspaces. Then

```math
 \gamma\le\frac{\sqrt\varepsilon}{1-\sqrt\varepsilon}.\qquad\text{(1)}
```

In particular, if the sign is the prescribed scalar phase for a single marked
a-bit address, tensored with identity on every other wire, then

```math
 \gamma\le\frac{2^{-a/2}}{1-2^{-a/2}}.\qquad\text{(2)}
```

The same conclusion holds for the complementary sign. There is no width,
T-count or Clifford-count restriction, and no address-preservation assumption.

### Native normal form and the Clifford trace bound

Write the actual native word as
$`V=C_2T_2C_1T_1C_0`$, allowing empty bands, and set

```math
 \begin{aligned}
 J&=T_1C_0,\\
 A&=T_2^\dagger C_2^\dagger S C_2T_2,\qquad
 B=J S J^\dagger,\\
 K&=(C_1^\dagger A C_1)B.
 \end{aligned}
```

A physical T band sends each Pauli, under conjugation, to a Hermitian Clifford:
each active wire contributes a Pauli or $`(X\pm Y)/\sqrt2`$. Thus A and B
are literal Hermitian Cliffords, K is a literal Clifford, and

```math
 W=J^\dagger KJ,\qquad \mathcal H=J^\dagger\mathrm{Re}(K)J.\qquad\text{(3)}
```

No common phase of V survives this identity.

Let M be K's binary symplectic matrix, let
$`r=\mathrm{rank}_{\mathbb F_2}(M-I)`$, and write
$`\tau(K)=2^{-w}\mathrm{Tr}K`$. The adjoint action on the Pauli
basis is a signed permutation. Its trace is $`|\mathrm{Tr}K|^2`$;
the signs on the fixed-label space $`\ker(M-I)`$ form a character. Hence

```math
 \tau(K)=0\quad\hbox{or}\quad |\tau(K)|^2=2^{-r}.\qquad\text{(4)}
```

In the nonzero case every fixed-label Pauli commutes with K literally. The
fixed space has dimension $`2w-r`$ and symplectic radical of dimension at
most r. When $`r\lt w`$, it contains at least $`w-r`$ symplectic pairs. A Clifford
coordinate change makes K identity on those qubits, so it acts on at most r
effective qubits; for $`r\ge w`$ this bound is automatic. This is a spectral
argument, not an emitted decoder. Every
nonempty eigenspace of K consequently occupies at least a $`2^{-r}`$
fraction of the full space.

Apply this to the smaller sign sector, replacing K by $`-K`$ if necessary.
If the trace is nonzero, $`\varepsilon\ge2^{-r}`$ and therefore
$`|\tau(K)|\le\sqrt\varepsilon`$. Since every positive eigenvalue of
$`\mathrm{Re}K`$ is at least gamma and every negative one is at least
minus one,

```math
 \gamma(1-\varepsilon)-\varepsilon
 \le\mathrm{Re}\tau(K)\le\sqrt\varepsilon.
```

Rearrangement gives (1). If the trace is zero, the stronger bound
$`\gamma\le\varepsilon/(1-\varepsilon)`$ gives the same conclusion.

### A Pauli-balanced lift retains the obstruction

Let f mark one a-bit address and add an independent bit u. The balanced table
$`g(x,u)=f(x)\oplus u`$ has equal positive and negative multiplicities, but
it still cannot have a uniform two-band source gap. Assume
$`\mathrm{sgn}\mathcal H=D_g\otimes I`$, and define

```math
 D'=\mathrm{sgn}(\mathrm{Re}K)=J D_g J^\dagger,
 \qquad Q=J Z_u J^\dagger,\qquad \varepsilon=2^{-a}.
```

Identity factors are suppressed. Q is a Hermitian Clifford because J contains
only one native T band. In normalized Hilbert--Schmidt norm,

```math
 \|D'-Q\|_{2,\mathrm{tr}}=2\sqrt\varepsilon,\qquad
 \|KQ-QK\|_{2,\mathrm{tr}}\le4\sqrt\varepsilon.\qquad\text{(5)}
```

The second inequality uses $`[K,D']=0`$. Every literal nonidentity Clifford
has distance at least $`\delta_{\mathrm C}=\sqrt{2-\sqrt2}`$ from I in
this norm: (4) bounds a nonscalar's real normalized trace by
$`1/\sqrt2`$, and the scalar cases are eighth roots of unity. Therefore

```math
 4\sqrt\varepsilon\lt\delta_{\mathrm C}
 \quad\Longrightarrow\quad [K,Q]=0.
```

This holds for $`a\ge5`$. The Clifford $`L=QK`$ now satisfies
$`\mathrm{Re}L=Q\mathrm{Re}K`$ and
$`\mathrm{sgn}(\mathrm{Re}L)=QD'=J D_fJ^\dagger`$.
Its gap is unchanged and its negative fraction is epsilon. Applying (1) to L
proves (2) for the balanced table as well.

The conclusion tolerates a small error in the sign eigenspaces. If
$`\|\mathrm{sgn}\mathcal H-D_g\otimes I\|_{\mathrm{op}}\le\delta`$,
replace the right side of (5) by $`2\delta+4\sqrt\varepsilon`$. Whenever
this is below $`\delta_{\mathrm C}`$, K and Q still commute. The two
Hermitian involutions $`Q\mathrm{sgn}(\mathrm{Re}K)`$ and
$`J D_fJ^\dagger`$ are then within delta in operator norm. For
$`\delta\lt2`$ their negative ranks agree, so (2) still follows. For the
unbalanced sparse target, this rank-stability argument alone suffices.

This is a zero-clean full-space screen for a Pauli signal and at most two
native source bands. It does not cover initialized-subspace sources,
non-Pauli signals or more source bands. A general nonlinear balancing mask
need not give a Clifford Q, so the balancing argument is specific to the
Pauli lift above. No growing lower bound on general native T-depth follows.

### Why the trace-support proof does not iterate to the next band

For three source bands, (3) can leave a two-band unitary rather than a
Clifford. Even a fixed positive gap and a large trace no longer force bounded
effective spectral support. Let C fan out one bit c to m bits y, let
$`E=\prod_{i=1}^m T_{y_i}`$, and take $`m\equiv1\pmod8`$, $`m\ge9`$.
The actual two-band word has full diagonal action

```math
 K_m=E^\dagger CEC,\qquad
 K_m|c,y\rangle=\omega^{c(m-2|y|)}|c,y\rangle,
 \qquad\omega=e^{i\pi/4}.
```

Its normalized trace is $`1/2+2^{-m/2-1}`$, its real-part gap is
$`1/\sqrt2`$, and its negative fraction is
$`1/4-2^{-(m+3)/2}`$. The four nontrivial eigenvalue multiplicities are
$`2^{m-2}\pm2^{(m-3)/2}`$, each occurring twice; eigenvalue one also has
multiplicity $`2^m`$. Their greatest common divisor is
$`2^{(m-3)/2}`$. Thus a factorization into an effective operator tensored
with identity needs at least $`(m+5)/2`$ effective qubits, despite the
uniformly large trace.

The same native construction works for every positive m. At $`m=3`$ it has
gap $`1/\sqrt2`$ and negative fraction $`1/8`$, violating (1). Thus even
the displayed two-band inequality cannot be carried unchanged to three bands.

This example also occurs in the actual next source normal form. Append a
qubit q. On c, a representative target $`y_1`$, and q, the three commuting
independent axes

```math
 P_1=Z_cY_q,\qquad P_2=Z_{y_1}Y_q,\qquad
 P_3=Z_cZ_{y_1}Y_q
```

give the one-band word
$`F=\exp(i\pi(P_1+P_2-P_3)/8)`$ up to a common phase, with
$`F X_qF^\dagger=\mathrm{CZ}_{c,y_1}H_q`$. A Hadamard on
$`y_1`$ and a Clifford change of the target-bank coordinates produce a
one-band U satisfying $`U X_qU^\dagger=C H_q`$. For the three-band source
$`V=U^\dagger E U`$, its Pauli pair signal is unitarily conjugate to
$`K_m\otimes I_q`$. Its spectral sign is not asserted to be an arbitrary
prescribed diagonal table. The example blocks an iteration of the Clifford
support argument; it neither supplies the missing readout nor excludes a
different spectral-volume bound.

The [finite gap fixtures](../../tests/test_native_gap_refinements.py) check
the full normal form on the returned-helper two-band CCZ, a literal one-band
realization of the displayed Pauli image, and the next-band multiplicities.
