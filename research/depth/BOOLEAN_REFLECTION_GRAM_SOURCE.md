# Boolean spectral reshaping and a two-reflection Gram source

[Paid phase and signal word](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md) ·
[Bounded-clean query promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Complete Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Changing the actual branch spectra removes the earlier branch-moment
obstruction. This note constructs those changed branches, including an
exact eight-phase evaluator, clean removal and a literal controlled
quarter-phase on arbitrary work. The resulting source has constant native
T-depth and the polynomial intermediate budgets. Its complete Gram is an
anticommutator of two broad reflections, rather than the old phase Gram.
The original flat bank mixer still fails, for every Boolean spectral
choice. A simultaneous change of the two reflections remains open.

## 1. The changed branches and their complete accepted action

Retain p, m=2^p, n=pm, the relative phase D_x, the fixed bank mixer
$`\mathcal R`$, and $`P_0=Z_{z_0}Z_{t_0}`$ from
[the paid signal construction](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#literal-controlled-t-wall).
In particular,

```math
D_x|z,t\rangle=\omega^{h_x(z,t)}|z,t\rangle,\quad
\omega=e^{i\pi/4},\quad D_x^8=I,\quad D_0=I.
```

Choose any two fixed truth tables $`f,g:\{0,\ldots,7\}\to\{0,1\}`$,
with f(0)=g(0)=0, and define spectral reflections

```math
R_f=(-1)^{f(D_x)},\qquad R_g=(-1)^{g(D_x)},\qquad
R=\mathcal R^\dagger R_f\mathcal R,\qquad S=P_0R_g.
                                                               \qquad(1)
```

Here f(D_x) means f(h) on the eigenphase $`\omega^h`$, not a Boolean
function of a matrix entry. Both R and S are Hermitian involutions.
The high-bit table f(h)=g(h)=[h at least 4] is one concrete choice;
all the following constructions also cover every other pair of tables.
Neither reflection is supplied as an oracle.

For this high-bit choice and nonzero x, both eigenspaces of R_f have
dimension at least one quarter of the bank space. Indeed h has the
distribution of a sum of $`n-p\ge1`$ independent residues taking
0 with probability one half and each of +1,-1 with probability one
quarter, by the phase-coordinate change used in the
[branch-moment calculation](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#branch-moment-obstruction).
Every residue probability is at most one half, since convolution cannot
increase the largest probability. Symmetry gives
$`\Pr(4\le h\le7)=(1-\Pr(h=0)+\Pr(h=4))/2`$, between one
quarter and three quarters. These reflections therefore avoid the thin
accepted-support premise of the Bell call-rank obstruction.

For $`P_\pm=(I\pm R)/2`$, put

```math
E=P_+-iP_-,\qquad V=iP_++\omega P_-.
```

Use the same literal one-signal word as before, but with its two branches
changed to $`A=iE`$ and $`B=iS`$:

```math
U=T(Y_a)H_a\mathrm C_a(iS)H_a\mathrm C_a(iE)H_aS_aT(Y_a).
                                                               \qquad(2)
```

Products act rightmost first. The physical basis envelope is
$`C_E^\dagger U C_E`$, as in the original source. Only a is initialized
and accepted at zero. Substitution into the full accepted polynomial,
without assuming any commutation, gives

```math
\begin{aligned}
K&=\frac{i}{2\sqrt2}(I+S)P_+
       +\frac{1+i}{4}(I-S)P_-\\
 &=\frac1{\sqrt8}(I+SR)V,\\
K^\dagger K&=\frac14\left(I+\frac{SR+RS}{2}\right).
\end{aligned}                                                   \qquad(3)
```

Indeed the general second line with an arbitrary unitary F instead of S
has Gram $`V^\dagger[I+\mathrm{Re}(FR)]V/4`$. When S is a Hermitian
involution, $`(SR+RS)/2`$ commutes with R and hence with V, proving the
last line. At x=0, R=I and S=P_0, so (3) is the required
$`(I+P_0)/4`$. These are complete operator identities on arbitrary
coherent address, both banks, all dirty work and every reference.

<a id="paid-boolean-spectral-evaluator"></a>

## 2. Exact eight-phase evaluation with returned clean selectors

The controlled relative phase is already paid, with common forward and
actual-inverse bounds $`T_C,G_C,d_C,h_C`$ in
[equation (6)](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#literal-controlled-t-wall).
Use three clean selector bits and the exact Fourier transform

```math
F_8|k\rangle=\frac1{\sqrt8}\sum_{j=0}^7\omega^{jk}|j\rangle,
\qquad
W=F_8^\dagger\left(\prod_{j=0}^2
       \mathrm C_{a_j}(D_x)^{2^j}\right)H^{\otimes3}.
                                                               \qquad(4)
```

On an eigenvector of phase $`\omega^h`$, W maps selector zero to h
exactly. This includes every degenerate eigenvector; it is not a phase
estimation accuracy assumption. Seven serial controlled-D calls implement
the powers in (4). F_8 uses three H gates, two CS gates, one CT and one
SWAP. A literal CS is three T/T-dagger gates and two CNOTs, by
$`2ab=a+b-(a\mathbin\oplus b)`$. The retained four-dirty CT word costs
70 T and 88 Cliffords. Consequently F_8 and its actual reverse have
T at most 76, Clifford count at most 98, T-depth at most 76 and four
returned dirty helpers. Including the initial three H gates, W has
Clifford bound $`7G_C+101`$.

For an arbitrary specified output u, compute its XOR with f(h), then
apply the actual W-dagger. A three-control XOR uses four returned
private-phase Toffolis and two dirty helpers: with arbitrary d, its
chronological logical word is

```math
\mathrm{CCX}(a,b;d),\quad \mathrm{CCX}(c,d;u),\quad
\mathrm{CCX}(a,b;d),\quad \mathrm{CCX}(c,d;u).
```

The two changes to u sum to abc and d returns. Each Toffoli costs
eight T gates and at most twelve Cliffords, using one reused private
phase bit. Surround by at most six literal X words for negative
controls. Thus one minterm costs at most 32 T, 72 Cliffords and depth
32. Emit at most eight minterms for any truth table. The complete
three-clean XOR source $`V_f=W^\dagger Q_f^{\rm table}W`$ therefore has

```math
\begin{aligned}
T_V&\le14T_C+408,&G_V&\le14G_C+778,\\
d_V&\le14d_C+408,&h_V&\le\max(h_C,4).
\end{aligned}                                                   \qquad(5)
```

Every clean selector returns exactly. Its complete logical action is
$`|x,z,t,u\rangle\mapsto|x,z,t,u\mathbin\oplus f(h_x(z,t))\rangle`$.
All its implementation helpers return on every input. The off-clean
action need not have this form, and is not assumed to do so.

## 3. Clean removal and the literal controlled quarter-phase

Apply [bounded-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
with k=3 and one output. Its exact negative three-bit zero query has
bounds 32 T, 72 Cliffords, depth 32 and two dirty helpers by the same
four-Toffoli construction. The promoted Q_f is therefore an exact
zero-clean XOR on its whole space, with common forward/reverse bounds

```math
\begin{aligned}
T_Q&\le32T_V+1248,&G_Q&\le32G_V+2754,\\
d_Q&\le32d_V+1152,&h_Q&\le h_V+7.
\end{aligned}                                                   \qquad(6)
```

The former selectors and old output are included in h_Q and are
arbitrary. The same majorants hold for Q_g. This application is valid
because (5) is a complete XOR query on arbitrary output u, rather than
only a phase or a state-preparation promise.

For the controlled reflection, borrow arbitrary u and emit chronologically
Q_g(u), CZ(a,u), Q_g(u), CZ(a,u). The complete phase is
$`(-1)^{a g(h)}`$, independently of u. Thus
$`\mathrm C_a(R_g)`$ costs two Q_g calls and six Cliffords.

A quarter-phase needs a separate returned identity. Borrow arbitrary u,v
and emit the following chronological word, with f=f(h) only in this
display's scalar verification:

```math
\begin{array}{llll}
Q_f(u),&\mathrm{CS}^\dagger(a,u),&Q_f(u),&\mathrm{CS}(a,u),\\
Q_f(v),&\mathrm{CCZ}(a,u,v),&Q_f(v),&\mathrm{CCZ}(a,u,v).
\end{array}                                                     \qquad(7)
```

The first row has phase
$`i^{a(u-(u\mathbin\oplus f))}=(-i)^{af}(-1)^{auf}`$.
The second row has phase $`(-1)^{auf}`$. Their product is exactly
$`(-i)^{af}`$, with both u and v restored and no residual scalar.
Two CS words cost six T and four Cliffords; two ordinary exact CCZ
words cost fourteen T and twenty Cliffords. Equation (7) therefore
pays the raw branch $`\mathrm C_a((-i)^{f(D_x)})`$ with

```math
T\le4T_Q+20,\quad G\le4G_Q+24,\quad
d\le4d_Q+20,\quad h\le h_Q+2.                         \qquad(8)
```

The actual inverse reverses the emitted gates. In particular, a
native implementation of Q_f need not be substituted for its reverse
merely because its full logical operator is an involution. Both
orientations have the bounds above. No clean phase kickback bit and no
controlled promoted query have been assumed.

<a id="two-reflection-factorization"></a>

## 4. Complete source ledger and the exact missing factorization

Conjugate (8) by the fixed bank mixer to control E. Control S by the
Q_g reflection echo followed by the two CZ factors controlling P_0.
The two i factors are literal S_a gates. With the same elementary
Clifford convention as the original source, a common ledger for (2),
its physical envelope and its actual inverse is

```math
\begin{aligned}
T_U&\le6T_Q+22,&d_U&\le6d_Q+22,\\
G_U&\le6G_Q+2G_{\mathcal R}+4c_n+60,&
h_U&\le h_Q+2,\qquad a_{\rm clean}=1.
\end{aligned}                                                   \qquad(9)
```

The 60 includes 24 quarter-phase Cliffords, twelve controlled-S-branch
Cliffords, and 24 for the signal rotations, three H gates and the three
signal S factors in a common inverse bound. Every internal helper
returns exactly. Both principal banks remain arbitrary; their embedding
around the m-wire output adds $`(2p-1)m`$ dirty wires. Since
$`T_C=O(r+n)`$, $`G_C=O(r+s+n)`$, $`d_C\le156`$ and
$`h_C\le\max(r,4n)`$, with $`r=(m-1)p^2`$ and
$`s\le(m-1)p^3`$, (9) has constant T-depth and count/width/Clifford
budgets within polynomial-overhead fusion. All displayed primitives are
exact; no dimension-dependent synthesis accuracy is hidden.

For a general pair of Hermitian involutions R,S, (3) identifies a
precise remaining operation. Write $`Q=Q_x=Z_z(v_x)Z_t(v_x)`$.
Then

```math
K^\dagger K=(I+Q)/4
\quad\Longleftrightarrow\quad
SR=Q
\quad\Longleftrightarrow\quad
S=QR\ \hbox{and}\ [R,Q]=0.                              \qquad(10)
```

To prove the nontrivial direction, set $`A=(SR+RS)/2`$ and
$`B=(SR-RS)/(2i)`$. The unitary SR has commuting Hermitian real and
imaginary parts, with $`A^2+B^2=I`$. If A=Q then A squared is I,
so B=0 and SR=Q. Conversely SR=Q is sufficient. This is an operator
factorization obligation, not a supplied synthesis of Q.

It is robust. If the Gram error against the target is at most delta,
then

```math
\begin{aligned}
\|SR-Q\|&\le4\delta+\sqrt{8\delta},\\
\|[S,R]\|&\le\sqrt{32\delta},\\
\|[R,Q]\|,\ \|[S,Q]\|&\le8\delta.
\end{aligned}                                                   \qquad(11)
```

Indeed $`\|A-Q\|\le4\delta`$, both A and Q are contractions,
and $`\|I-A^2\|\le8\delta`$. This bounds B and hence SR-Q.
Both R and S commute with A, giving the last line. The estimates apply
also with coherent address motion and arbitrary spectator factors,
whenever the complete two-reflection Gram (3) is the actual source.

<a id="flat-spectral-reshaping-gap"></a>

## 5. Every spectral table still fails with the original flat mixer

This family changes the old branch spectra and has the correct normalized
trace. Nevertheless it cannot match the whole lookup sector. At any fixed
nonzero x, every computational column of $`\mathcal R`$ has uniform
modulus. Hence the diagonal entries of R in (1) are all
$`\mu_f=2^{-2n}\mathrm{Tr}(R_f)`$. S is diagonal. On z=0 the
relative phase is D_x=I for every t, so S=P_0 because g(0)=0.

Choose j different from zero with $`(v_x)_j=1`$; the second field
block supplies such a j. The computational vectors
$`w_+=|0,0\rangle`$ and $`w_-=|0,e_j\rangle`$ have opposite
lookup charges, but both are S eigenvectors of sign +1. Equation (3)
therefore gives

```math
\langle w_+|K^\dagger K|w_+\rangle
=\langle w_-|K^\dagger K|w_-\rangle
=\frac{1+\mu_f}{4},\qquad
\left\|K^\dagger K-\frac{I+Q_x}{4}\right\|\ge\frac14.
                                                               \qquad(12)
```

The two desired expectations are one half and zero, so the last bound
is their minimax separation. It holds for every pair of truth tables,
not just the high-bit example, and for every nonzero address and every
p at least one. Even dropping g(0)=0 merely changes both equal
expectations to $`(1+(-1)^{g(0)}\mu_f)/4`$. Flatness and the two
opposite-charge vectors, rather than the old spectral moments, are the
premises. The proof does not exclude a simultaneous change of the
spectra and the bank dynamics.

The normalized trace claim follows directly: R is X-diagonal and
therefore commutes with $`X_{z_0}`$, while S anticommutes with it because
D_x is independent of z_0 and P_0 changes sign. Thus
$`\mathrm{Tr}(SR)=0`$ and (3) has normalized trace one quarter.
Trace matching has been achieved, but it does not certify sectors.

The next concrete question in this interface is to construct two paid
broad reflections whose product approaches Q as in (10)--(11), with
genuinely changed bank dynamics. Merely changing the eight spectral
signs under the original flat mixer is paused by (12). The complete
indicator and frame-depth frontiers do not move. The proof is analytic;
the accompanying exact fixtures are finite internal checks, not external
mathematical review.
