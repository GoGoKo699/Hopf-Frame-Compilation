# Limits of structured frames for spectral promotion

[Parallel clean-source promotion](PARALLEL_CLEAN_PROMOTION.md#6-alternative-a-constant-spectral-gap-permits-a-paid-scalar-filter)
requires a positive whole-space gap for a sum of conjugated tensor
reflections. This note proves three separately scoped obstructions.
They use the actual marked-projector tensor structure. None excludes
general dense frames with noncommuting global Pauli bases.

Let w=k+1, M=2^k, h=M and n=hw. Partition n qubits into h blocks of
w qubits, put s=0^k1, and let p_j=|s><s| on block j. Define

```math
\begin{aligned}
D_g&=\prod_{j=1}^h(I-2p_j),\\
H_{\mathcal C}&=\sum_{i=1}^a\frac{I-C_iD_gC_i^\dagger}{2}.
\end{aligned}                                                   \qquad(1)
```

The spectral alternative needs H_C≥gamma I for a constant gamma>0.
For one block the exact Pauli coefficients are

```math
I-2|s\rangle\langle s|
=\left(1-\frac1M\right)I
 -\frac1M\sum_{0\ne v\in\mathbb F_2^w}(-1)^{v\cdot s}Z^v.
                                                               \qquad(2)
```

The coefficients of D_g tensorize. A common +1 vector for the framed
reflections, or a normalized vector with energy tending to zero in (1),
disproves the required gap.

## 1. Frames sharing a Lagrangian Pauli space

Suppose the frames are Clifford and all their marked stabilizer label
spaces lie in one common isotropic Pauli label space. Here a label is
the binary X/Z vector with its sign omitted; isotropic means that all
represented Paulis commute. In particular this hypothesis holds for
arbitrary dense monomial Clifford frames, whose marked labels are all
Z-type. It also permits a common nonstandard Lagrangian basis.

If a≤w, H_C has an exact kernel containing a stabilizer code of dimension
2^(n−ah). To prove this, let S_{i,j} be the w-dimensional label space
of the marked stabilizer group for C_i p_j C_i†. Within each fixed frame
these h spaces form a direct sum, since the original block Z labels do.
For any collection of t marked groups, one frame contributes at least
ceil(t/a) of them. Therefore

```math
\dim\sum_{(i,j)\in I}S_{i,j}
 \ge w\left\lceil\frac{|I|}{a}\right\rceil
 \ge |I|.                                                     \qquad(3)
```

The linear transversal lemma gives one globally independent nonzero
label from each S_{i,j}. For completeness, its sufficient condition is
dim(sum_{j in I} S_j)≥|I| for every I. An induction proves it. If a proper
nonempty subset is tight, choose its independent representatives first
and apply induction to the other spaces modulo its span; the condition
survives in that quotient. If there is no proper tight subset, choose
any nonzero vector from the final space and quotient by that line.
Every remaining nonempty subset had at least one unit of rank slack,
so induction again applies. The one-space case starts the induction.

Each chosen label has a definite signed marked-state stabilizer. Because
all chosen labels commute and are independent, imposing eigenvalue −1
for every such signed stabilizer is consistent and leaves dimension
2^(n−ah). Every marked projector requires eigenvalue +1 for its chosen
stabilizer, so it kills this entire code. Thus every factor in each
framed tensor reflection fixes the code, proving the claim.

The common isotropic space is essential. A linear transversal alone
does not make labels from different noncommuting Lagrangian spaces
commute. This proof does not extend to arbitrary global Clifford frames.

## 2. Bounded support cones, allowing non-Clifford gates

Now allow arbitrary unitary frames. Write p_{i,j}=C_i p_j C_i† and suppose
their physical-support overlap graph has maximum degree Delta. If

```math
e(\Delta+1)\le2^w,                                           \qquad(4)
```

then H_C has a zero eigenvalue. The transformed marks need not commute
with one another across different frames.

Each good subspace X_{i,j}=ker p_{i,j} has relative dimension 1−2^(-w),
because conjugation preserves rank. If p_{i,j} is supported on S, all
its nonneighbors act on the complementary tensor factor. Their joint
kernel therefore factors as H_S tensor K, even if those nonneighbors
overlap or do not commute. Intersecting it with X_{i,j} multiplies its
relative dimension by exactly 1−2^(-w). This is the mutual R-independence
condition for the good kernel subspaces.

Theorem 4 of Ambainis, Kempe and Sattath,
[A Quantum Lovasz Local Lemma](https://arxiv.org/pdf/0911.1696),
gives nonzero intersection for good subspaces of relative dimension
at least 1−p when their mutual R-independence dependency degree is Delta
and ep(Delta+1)≤1. Apply it with p=2^(-w). Every vector in the resulting
common kernel is fixed by all factors I−2p_{i,j}, and hence by every
framed reflection. This proves (4).

In particular, suppose each frame has ordinary depth ell in layers of
disjoint two-block gates on the fixed w-qubit partition. The gates may
be arbitrary unitaries. Each marked block reaches at most 2^ell output
blocks. The backward cone of any output block contains at most 2^ell
original blocks per frame. Thus one transformed mark overlaps at most
a 4^ell marks, including itself, and the sufficient obstruction is

```math
e a4^\ell\le2^{k+1}.                                        \qquad(5)
```

For constant a, ordinary block depth ell=o(k) therefore gives zero gap
for all sufficiently large k. This is not a native T-depth lower bound.
Dense Clifford layers can have large block cones while costing no
T-depth. A global Clifford followed by shallow non-Clifford dressing
is outside (5) if the global Clifford already makes those cones large.

## 3. Explicit global graph frames and tensor witnesses

For k≥2 split the h blocks into equal paired halves L,R, each with
r=n/2 bits. For binary u,v of length r define Bell labels by

```math
|\Phi_{u,v}\rangle
 =2^{-r/2}\sum_y(-1)^{v\cdot y}|y+u\rangle_L|y\rangle_R.
                                                               \qquad(6)
```

For a binary r by r matrix A define

```math
\begin{aligned}
Z_A|x,y\rangle&=(-1)^{x^TAy}|x,y\rangle,\\
C_A&=H^{\otimes 2r}Z_AH^{\otimes 2r}.
\end{aligned}                                                   \qquad(7)
```

These are literal Clifford frames. Dense A uses at most r^2 CZ gates
and 4r Hadamards. Decomposing each CZ into H-CNOT-H gives the explicit
native bound G≤3r^2+4r, T=d_T=0, with no helpers. The graph matrix is
B_A=[[0,A],[A^T,0]], which is full rank whenever A is invertible.

### 3.1 Alternating matrices retain an exact Bell vector

If A is symmetric with zero diagonal, x^TAx=0. The vector Phi_{0,0}
is fixed by the two-sided Hadamard and by Z_A, so C_A fixes it. Also
D_g fixes it: the two identical tensor halves contribute equal signs,
whose product is one. Consequently any number of such frames has

```math
C_AD_gC_A^\dagger|\Phi_{0,0}\rangle=|\Phi_{0,0}\rangle.
                                                               \qquad(8)
```

Invertible examples exist because r is even. For r=2d take
A=[[0,M_0],[M_0^T,0]] with any invertible d by d matrix M_0. Multiplication
matrices of nonzero elements of F_(2^d) give an explicit family. Thus
full-rank, potentially dense global graph mixing does not itself
establish the gap.

### 3.2 Symmetric trace pencils have exact orbit trial energies

Let A be symmetric, now allowing a nonzero diagonal d. Directly from
(6), two-sided Hadamard maps Phi_{u,v} to (−1)^(u·v) Phi_{v,u}, and
Z_A maps Phi_{u,v} to Phi_{u,v+d+Au}. Combining these identities gives
the complete action, with no residual phase,

```math
C_A|\Phi_{u,v}\rangle=|\Phi_{u+d+Av,v}\rangle.                  \qquad(9)
```

For symmetric frames A_1,...,A_a and a fixed v put

```math
\begin{aligned}
T_v&=\mathrm{span}_{\mathbb F_2}\{d_i+A_i v:1\le i\le a\},\\
|\psi_v\rangle
 &=|T_v|^{-1/2}\sum_{u\in T_v}|\Phi_{u,v}\rangle.
\end{aligned}                                                   \qquad(10)
```

Every frame fixes psi_v by translation of the sum. Let t(u) count the
nonzero w-bit blocks of u in the paired halves. On a zero block shift,
the two signs of the local reflection agree. On a nonzero shift they
differ at exactly the two strings s and s+u_j. Their average product
is therefore 1−4/2^w=1−2/M. Different u labels have disjoint
computational support, and D_g preserves that support. It follows that

```math
\langle\psi_v,H_{\mathcal C}\psi_v\rangle
 =\frac a2\left[1-
   \frac1{|T_v|}\sum_{u\in T_v}\left(1-\frac2M\right)^{t(u)}
  \right].                                                    \qquad(11)
```

This exact energy uses the tensor autocorrelation from (2). In
particular, average block weight o(M) in some T_v gives vanishing energy.

There are explicit full-rank trace-pencil examples with this behavior.
Choose a basis beta_1,...,beta_r of F_(2^r), and for nonzero field
parameter b define

```math
(A_b)_{ij}=\mathrm{Tr}(b\beta_i\beta_j).                        \qquad(12)
```

The trace pairing makes A_b symmetric and invertible. Its diagonal is
d_b,i=Tr(sqrt(b) beta_i). Thus b maps bijectively onto diagonal vectors,
where field elements are represented by these trace coordinates.
Choose nonzero parameters whose diagonal vectors all lie in one
physical w-bit block. At v=0 every nonzero u in T_0 then has t(u)=1,
and (11) becomes

```math
\langle\psi_0,H_{\mathcal C}\psi_0\rangle
 =\frac{a(1-2^{-\dim T_0})}{M}\longrightarrow0.                \qquad(13)
```

Each graph form remains full rank. This does not cover arbitrary trace
parameters: a globally dispersed translation span may have block weight
proportional to M. Eliminating the displayed trial vectors would also
not prove a gap, since (11) does not control states superposing different
Bell-v sectors. D_g mixes v at fixed u, while the frames shear u as in
(9); those couplings remain part of the whole-space spectral problem.

## 4. What remains unpaid

The three hypotheses are separate: one shared commuting Pauli space,
small physical support cones, and the explicit graph-frame families or
low-weight Bell orbits. None excludes dense frames with different
noncommuting Pauli bases and sufficiently dispersed label action.
No generic expander estimate, normalized trace estimate or absence of
these witnesses supplies the required least-eigenvalue bound.

The first unpaid operation in the spectral promotion alternative remains
an explicit, uniformly paid frame family with a constant whole-space
gap in (1). Its constant-clean filter and subsequent promotion are paid
conditionally in the linked proof. An actual qualifying clean source
is separately still required, so these scope results do not improve the
complete unrestricted frontier.
