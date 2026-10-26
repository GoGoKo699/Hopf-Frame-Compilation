# A fixed non-Gaussian readout cannot supply the missing lookup label

[Minimal tree source](ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source) ·
[Affine query boundary](GAUSSIAN_INVERSE_EXTRACTION.md#4-a-tight-affine-symplectic-round-obstruction) ·
[Complete observable contract](OBSERVABLE_TRANSPORT_SOURCE.md)

Appending one fixed carrier unitary to a low-degree Pauli source does not
give a shallow complete lookup probe. The final unitary may be arbitrary,
including non-Clifford and non-Gaussian dynamics on every carrier and dirty
helper. The obstruction follows from the number of joint Pauli syndromes
that its input family can represent. It does not assume that this fixed
unitary has a cheap native implementation.

For affine-symplectic preprocessing, exact transport still needs at least p
query stages. At any fixed operator error below the threshold below, it needs
at least p/2 minus an explicit square-root correction. These are bounds on
the stated fixed-address source grammar, not on unrestricted native T-depth.
Address motion and non-Gaussian interactions interleaved with later
address-dependent stages remain outside the theorem. R96 is unchanged.

## 1. The complete family and the fixed-readout grammar

Let p be positive and put $`m=2^p`$. The coefficient carrier has m arbitrary
bits $`c_S`$, one for each subset of $`[p]`$. Include every arbitrary helper
in the w-qubit physical carrier space B. The desired fixed-address
observables are

```math
Q_x=(-1)^{\sum_{S\subseteq[p]}c_S\prod_{i\in S}x_i}\otimes I_H
    =Z(v_x)\otimes I_H,
\qquad (v_x)_S=\prod_{i\in S}x_i.                         \qquad(1)
```

The m binary labels $`v_x`$ are independent: indexed by the support of x,
their matrix is the subset-zeta matrix, triangular with diagonal one in
an inclusion-compatible order. Thus every string of m joint eigenvalues
of the $`Q_x`$ occurs on a carrier computational state, with any helper state.
This is a test of the full operator contract, not an initialized carrier
supplied to a construction.

Fix a nonidentity Hermitian Pauli P on B. Suppose $`C_x`$ is a source whose
literal probe images are signed Hermitian Paulis,

```math
P_x=C_xPC_x^\dagger=\sigma_x\tau(a(x)),\qquad
\sigma_x\in\{+1,-1\},\qquad \deg a_j\le d.             \qquad(2)
```

Here a(x) is the full binary X/Z label, including every helper coordinate.
Only its coordinate degrees are bounded. The signs may be arbitrary
functions of x. The $`C_x`$ need not be Clifford for this general statement;
only their specified images in (2) must be Pauli.

Let U be one fixed unitary on B, independent of x. It may mix all carrier
and helper wires, but it does not act on the address. Consider the source
with fixed-address fibers $`V_x=UC_x`$ and require

```math
\max_x\|UP_xU^\dagger-Q_x\|_{\rm op}\le\delta,
\qquad 0\le\delta\lt\gamma:=\sqrt{2-\sqrt2}.          \qquad(3)
```

This is the complete operator norm, equivalent to the coherent block
condition with address retained. It includes every dirty input and any
external reference. Neither $`C_x`$ nor U must return its work separately.
There is no additional address-dependent stage after U in this grammar.
If an emitted source is within eta of this complete ideal source and has
probe error tau, unitary telescoping gives (3) with $`\delta\le\tau+2\eta`$.
Its inverse is the actual adjoint, so the same estimate includes work return.

## 2. Syndrome-covering theorem

Define

```math
M_d=\sum_{j=0}^{\min(d,p)}\binom pj,
\qquad r=\left\lfloor\frac{m\delta^2}{4}\right\rfloor.
```

**Theorem.** Under (1)--(3), the $`P_x`$ commute pairwise. Their binary-label
span has some dimension $`R\le M_d`$, and

```math
2^R\sum_{j=0}^{r}\binom mj\ \ge\ 2^m.                 \qquad(4)
```

In particular, with binary entropy $`H_2(t)=-t\log_2t-(1-t)\log_2(1-t)`$
and $`H_2(0)=0`$,

```math
M_d\ge m\bigl[1-H_2(\delta^2/4)\bigr].                \qquad(5)
```

### Small complete error first forces commutation

Choose a common positive eigenstate of all $`Q_x`$ and pull it back through U.
Equation (3), applied to that state, gives

```math
\langle P_x\rangle\ge1-\delta^2/2.
```

If two Hermitian involutions $`P_x`$ and $`P_y`$ anticommute, then
$`\langle P_x\rangle^2+\langle P_y\rangle^2\le1`$.
Indeed the square of $`sP_x+tP_y`$ is $`(s^2+t^2)I`$; choose s,t to
be the two expectations. Two expectations at least $`1-\delta^2/2`$
contradict this inequality when $`\delta\lt\sqrt{2-\sqrt2}`$.
Signed Paulis either commute or anticommute, proving pairwise commutation.
The threshold is used only for this step; no sharpness claim for the
whole transport problem follows from it.

### Arbitrary signs give an affine code, not extra label directions

Write the Boolean expansion of the label vector as

```math
a(x)=\sum_{S:\,|S|\le d}a_S\prod_{i\in S}x_i.
```

Its values lie in the span of at most $`M_d`$ coefficient vectors, so their
rank R is at most $`M_d`$, independently of physical width. A commuting
signed-Pauli family of label rank R has exactly $`2^R`$ joint syndrome
strings. They form an affine binary code $`\mathcal C\subseteq\mathbb F_2^m`$:
choose R independent generators; every other $`P_x`$ is their product with
a fixed sign. Every assignment to the independent generators occurs,
since its joint projector has rank $`2^{w-R}`$.

For any desired target syndrome $`y\in\mathbb F_2^m`$, choose a Q-family
joint eigenstate and pull it back through U. Measuring the commuting
P-family gives a random codeword $`Z\in\mathcal C`$. For each address x,

```math
4\Pr[Z_x\ne y_x]
 =\|(P_x-(-1)^{y_x}I)U^\dagger|\psi_y\rangle\|^2
 \le\delta^2.
```

Summing gives $`\mathbb E\,d_H(Z,y)\le m\delta^2/4`$. Some codeword
therefore has integer distance at most r from y. All $`2^m`$ syndromes
must be covered by the $`2^R`$ radius-r Hamming balls, proving (4).
Overlaps of these balls only make the necessary inequality stronger.
The binomial-sum bound
$`\sum_{j\le r}\binom mj\le2^{mH_2(r/m)}`$, together with monotonicity
of $`H_2`$ below one half, proves (5). This bound follows directly by expanding
$`(t+(1-t))^m`$ at $`t=r/m`$; the r=0 case is exact.

## 3. Consequences for the paid affine queries

For exact transport, r=0 in (4), so $`R\ge m`$ and $`M_d\ge2^p`$.
Consequently $`d\ge p`$. For nonzero error put

```math
\alpha_\delta=1-H_2(\delta^2/4)\gt0.
```

Since $`M_d/2^p=\Pr[\mathrm{Bin}(p,1/2)\le d]`$, the binomial lower-tail
bound gives the explicit necessary condition

```math
d\ge\frac p2-
 \sqrt{\frac p2\ln\frac1{\alpha_\delta}}.              \qquad(6)
```

For completeness, when $`d\le p/2`$, the exponential-moment estimate
$`\Pr[\mathrm{Bin}(p,1/2)\le d]\le\exp[-2(p/2-d)^2/p]`$ follows from
$`\cosh t\le\exp(t^2/2)`$ and minimization of the moment parameter.
For $`d\gt p/2`$, (6) is immediate. The exact-error conclusion d at least
p is stronger than this tail estimate at delta zero.

Now take $`C_x`$ to consist of q fixed-address Clifford query stages, each
with binary symplectic entries affine in x, interleaved with arbitrary
address-independent Clifford walls on B. The resulting symplectic matrix
is a product of q affine matrices and constant matrices. Its entries, and
therefore the labels in (2), have degree at most q. This includes both
orientations of the minimal tree's paid forward XOR query, carrier Fourier
walls, fixed bank mixing and arbitrary dirty spectators.

Equations (4)--(6) apply with d=q. Exact transport needs q at least p;
at fixed error below gamma it needs $`q\ge p/2-O_\delta(\sqrt p)`$.
Appending an arbitrary fixed non-Gaussian carrier decoder to those query
stages therefore does not supply a constant- or sublogarithmic-stage root
readout. Making its final unitary larger or more expressive does not change
the syndrome deficit.

This conclusion does not add q to the native depth of a jointly resynthesized
word. It does not cover a decoder that moves the address, further
address-dependent processing after U, non-Pauli intermediate probes before
the final fixed readout, or a different complete native representation.
Those mechanisms require their own construction and simultaneous count,
width, precision and work-return ledger. The first unpaid operation remains
such a complete joint root transport, rather than a final fixed decoder of
the available affine-query Pauli family.

The [finite fixtures](../../tests/test_fixed_readout_syndrome_bound.py) check
target-syndrome completeness, signed source-code ranks and covering radii,
the anticommuting expectation threshold and the mismatch identity on full
carrier inputs. They supplement the all-size proof; they do not synthesize
an unrestricted source or establish a general native-depth lower bound.
