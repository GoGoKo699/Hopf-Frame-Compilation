# Uniform Clifford transport from bounded Pauli deficiency

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Diagonal-source transport](NATIVE_CLIFFORD_TRANSPORT.md)

Bounded Pauli deficiency gives a source-uniform transport theorem even
when the Pauli-preserving domain is degenerate. The proof does not require
third-level membership or exact rounding of the transported endpoint. It
uses the full operator approximation first, then normalized
Hilbert--Schmidt estimates to recover the exact target Clifford's labels
and literal signs.

The deficiency bound is an additional hypothesis on the specified shared
prefix. It need not hold uniformly for general three-band prefixes, so
this result does not resolve the general six-layer degree question.

## 1. Statement and native contract

Let F be a fixed unitary on w arbitrary work qubits, and define

```math
E_F=\{p:F\tau(p)F^\dagger\text{ is a signed Hermitian Pauli}\},
\qquad k=\mathrm{codim}_{\mathbb F_2}E_F.
```

The zero label belongs to this domain. Use the tensor-factorizing
Hermitian section
$`\tau(x,z)=i^{x\cdot z}X^xZ^z`$, with the dot product in the exponent
taken over the integers. In particular
$`\tau(p_A,p_B)=\tau(p_A)\otimes\tau(p_B)`$ literally.

Let L and K be the exact native inner family and actual outer comparison
of the [six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md#2-a-clifford-valued-native-pullback).
They are Clifford up to their retained scalar phases and satisfy

```math
\|K_x-FL_xF^\dagger\|_{\rm op}\le2\eta.                 \qquad(1)
```

The native input ledger gives degree at most two for the image labels of
L on fixed probes, and degree at most six for their literal signs. Fixed
Clifford coordinate changes preserve these bounds.

**Bounded-deficiency theorem.** For $`k\ge1`$, put

```math
D(k)=(2k+1)(4k+6)=8k^2+16k+6.
```

If

```math
\eta\lt\frac{1}{4\sqrt2(2^k+1)},                         \qquad(2)
```

then

```math
\deg_x S_{K_x}\le D(k),\qquad
\deg_x\kappa_{K_x}(q)\le3D(k)
\quad\text{for every fixed Hermitian Pauli probe }q.      \qquad(3)
```

Consequently the actual zero-clean read-only six-layer query has
$`\deg f\le24k^2+48k+18`$, improving to $`8k^2+16k+6`$ when the outer
target shear is visible. The bounds are independent of dirty width for
fixed k. Neither sharpness nor attainment is asserted.

More precisely, k can be replaced in (2)--(3) by

```math
a=\frac{k+\dim\mathrm{rad}E_F}{2}\le k,             \qquad(4)
```

the active-subsystem size below. For $`k=0`$, F is Clifford. Direct
Clifford separation gives label degree two and fixed-probe sign degree
six whenever $`\eta\lt\sqrt2/4`$, with the same respective target bounds.

The complete query contract remains the one supplying (1): every address
and logical column, arbitrary output and dirty input, external reference,
literal target phase, leakage bound and full workspace return are
retained. No initialized-only source promise is substituted for it.

## 2. Signed localization of a possibly degenerate domain

E is a linear space: conjugation preserves Pauli products and their
literal scalar phases. Its induced label map preserves the symplectic
pairing. Write $`r=\dim\mathrm{rad}E_F`$. A symplectic basis of
E contains a nondegenerate subspace $`E_0`$ of dimension $`2w-k-r`$.
Its signed action extends to a fixed Clifford action: extend paired
bases in domain and image to full symplectic bases, then correct the
remaining basis signs with a fixed Pauli. Multiplication matches the
literal action on the whole subspace.

Fixed input and output Clifford coordinates put $`E_0`$ and its image
on the full Pauli algebra of a spectator subsystem B. Matching the signs
makes the remaining source commute literally with every B Pauli. The
commutant of that full algebra is therefore

```math
F=C_{\rm out}(U_A\otimes I_B)C_{\rm in},\qquad
a=w-\tfrac12\dim E_0=(k+r)/2.                            \qquad(5)
```

The radical is contained in $`E_F^\perp`$, so $`r\le k`$ and
$`a\le k`$. Its directions remain inside A; no controlled factorization
of the residual source is needed. A single T gate, for example, has
$`k=r=a=1`$ and a degenerate Pauli-preserving domain.

These fixed coordinates are analytical. They introduce no physical gate,
source call, inverse, control, helper or clean qubit. Work below in the
coordinates $`F=U_A\otimes I_B`$, with L and K conjugated accordingly.

## 3. Uniform recovery from heavy Pauli coefficients

Fix any Hermitian probe $`\tau(q_A)\otimes\tau(q_B)`$ and expand

```math
(U_A^\dagger\tau(q_A)U_A)\otimes\tau(q_B)
 =\sum_{s\in\mathbb F_2^{2a}}c_s\tau(s,q_B),
\qquad c_s\in\mathbb R,\quad\sum_sc_s^2=1.
```

There are $`m=4^a`$ possible coefficients. Set

```math
\delta=4\eta,\qquad
\rho=\frac1{\sqrt2(2^a+1)},\qquad
\mathcal H=\{s:|c_s|\gt\rho\}.                          \qquad(6)
```

The radius gives $`\delta\lt\rho`$. The retained set is nonempty because
some coefficient has magnitude at least $`2^{-a}\gt\rho`$. All these
coefficients, the retained set and choices made from it are independent
of x.

Write the literal Pauli images as

```math
\begin{aligned}
L_x\tau(s,q_B)L_x^\dagger
 &=(-1)^{\gamma_s(x)}\tau(u_s(x),v_s(x)),\\
K_x\tau(q)K_x^\dagger
 &=\sigma_x\tau(t_x,b_x),\qquad\sigma_x\in\{\pm1\}.
\end{aligned}
```

The labels u,t are on A and v,b on B. Conjugating (1) on the probe and
then by $`F^\dagger`$ yields

```math
\left\|\sum_s c_s(-1)^{\gamma_s(x)}\tau(u_s(x),v_s(x))
 -\sigma_x(U_A^\dagger\tau(t_x)U_A)\otimes\tau(b_x)\right\|_{\rm op}
 \le\delta.                                             \qquad(7)
```

Now use the normalized Hilbert--Schmidt norm

```math
\|M\|_{2,\rm n}=\sqrt{\mathrm{Tr}(M^\dagger M)/\dim M}
\le\|M\|_{\rm op}.
```

Distinct source terms have distinct full Pauli images, since L is
Clifford. If $`v_s\ne b_x`$, the target in (7) has zero coefficient on
that Pauli, so $`|c_s|\le\delta`$. Consequently

```math
v_s(x)=b_x\quad(s\in\mathcal H).                        \qquad(8)
```

Choose one fixed $`s_0\in\mathcal H`$. The entire spectator output label
is $`b_x=v_{s_0}(x)`$, of native degree at most two.

Project (7) onto the B-Pauli sector $`\tau(b_x)`$. The local coefficient
map
$`M\mapsto2^{-|B|}\mathrm{Tr}_B[(I_A\otimes\tau(b_x))M]`$
is a contraction between the normalized Hilbert--Schmidt norms. Let
$`T_x`$ be the local A contribution of the retained terms alone. The
discarded terms have squared norm at most $`m\rho^2`$ before projection,
so

```math
\|T_x-\sigma_xU_A^\dagger\tau(t_x)U_A\|_{2,\rm n}
\le\delta+2^a\rho\lt\frac1{\sqrt2}.                    \qquad(9)
```

This is a Hilbert--Schmidt tail estimate, not an operator-norm estimate.
Distinct signed operators $`\sigma U_A^\dagger\tau(t)U_A`$ have mutual
normalized Hilbert--Schmidt distance at least $`\sqrt2`$, independently
of U. Their open balls of radius $`1/\sqrt2`$ are disjoint. Thus the
local output label $`t_x`$ and its literal sign $`\sigma_x`$ are uniquely
determined by $`T_x`$. No source-specific gap to physical Paulis is used.

## 4. A small signed feature set determines the output

Choose a fixed basis $`d_1,\ldots,d_h`$ for the span of
$`\{s+s_0:s\in\mathcal H\}`$, from those differences. It has
$`h\le2a`$. By (8), multiplying the corresponding full Pauli images
shows that

```math
L_x\tau(d_j,0)L_x^\dagger
 =(-1)^{\epsilon_j(x)}\tau(z_j(x),0).                    \qquad(10)
```

Every difference image is supported entirely on A. Its sign is retained;
this is not merely an assertion about binary labels.

For each retained s, fixed binary coefficients $`e_j(s)`$ and a fixed
fourth root of unity $`\alpha_s`$ satisfy, in one fixed ordering,

```math
\tau(s,q_B)=\alpha_s
 \left[\prod_{j=1}^h\tau(d_j,0)^{e_j(s)}\right]\tau(s_0,q_B).
                                                               \qquad(11)
```

After conjugation, (11) reconstructs $`T_x`$ literally from the active
label $`u_{s_0}(x)`$ and sign $`\gamma_{s_0}(x)`$, together with the h
active labels $`z_j(x)`$ and signs $`\epsilon_j(x)`$. The one spectator
factor $`\tau(b_x)`$ factors out. No equality test over an unbounded
number of spectator bits is needed, and no multiplication phase in
$`\alpha_s`$ is discarded.

Each of the $`(h+1)2a`$ label bits has degree at most two, and each of the
$`h+1`$ sign bits has degree at most six. The unique local output label
bits and sign from (9) are Boolean functions of the realized feature
tuples. Extend each function arbitrarily to all binary feature tuples
and take its multilinear polynomial. Substitution bounds its degree by
the sum of the feature degree bounds:

```math
(h+1)(4a+6)\le(2a+1)(4a+6)=D(a).                        \qquad(12)
```

The local output labels and sign thus have degree at most $`D(a)`$;
the spectator labels still have degree two. Restoring the fixed output
Clifford preserves label degree. Its Hermitian-section sign is cubic in
the incoming label, so it adds degree at most $`3D(a)`$. This proves (3).
The argument applies directly to each fixed probe; it does not multiply
an unbounded number of output-generator signs.

For $`a=0`$, F is Clifford. Equation (1), projective Clifford separation
and $`2\eta\lt1/\sqrt2`$ make K projectively equal to $`FLF^\dagger`$.
The native degrees are then two and six, proving the separately stated
zero-deficiency case.

## 5. Target extraction, resources and limits

For visible outer target shear, the retained
[fixed linear extractor](READ_ONLY_FIVE_LAYER_REDUCTION.md#3-visible-outer-shear-a-translation-invariant-space-of-functions)
expresses $`g(x)=f(x)+f(0)`$ as a fixed linear functional of $`S_K(x)`$
plus a quadratic baseline. This gives degree at most $`D(a)`$. When the
outer shear vanishes, a fixed anticommuting probe has sign
$`p_q(x)+g(x)`$, with native baseline degree at most six. Its literal
sign bound gives degree at most $`3D(a)`$.

The proof applies to the actual comparison and complete operator error.
It neither assumes intermediate source passages are Clifford nor infers
that $`FL_xF^\dagger`$ is exactly Clifford at positive error. Actual
scalars disappear only in conjugation of Hermitian probes, after (1).
Every original T gate, elementary Clifford gate, actual inverse, control
and simultaneous work allocation remains charged. The fixed coordinate
choices, projection and interpolation are analysis, not new physical
operations, source calls, helpers or initialized states.

The active-subsystem size and error radius are uniform over sources with
bounded deficiency. However, the actual three-band prefix can have k
growing with work width; a wide parallel T band already does. The radius
in (2) decreases exponentially in k. Therefore this theorem supplies no
width-independent general six-layer result, unrestricted one-clean or
address-mixing rigidity, growing-depth recurrence, or improvement to
complete-frame resource upper bounds.

The [focused checks](../../tests/test_bounded_deficiency_transport.py)
verify degenerate signed localization, heavy-term common spectator
labels, phase-sensitive reconstruction from a small difference basis,
and genuinely approximate endpoints whose transported operators remain
non-Clifford. These are finite checks of the proof mechanisms, not
returned-query compilers or evidence that all native prefixes have
bounded deficiency. The argument above establishes the uniform theorem.
