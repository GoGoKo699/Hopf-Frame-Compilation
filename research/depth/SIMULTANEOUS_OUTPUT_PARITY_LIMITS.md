# Simultaneous output parities and native support rank

[Transported output probes](OVERLAPPING_TWO_BAND_SOURCES.md#2-transported-output-probe-for-general-encoders) ·
[Fixed Clifford-sum source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) ·
[Complete source recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#61-exact-flatness-sharpens-the-complete-error-and-return-recurrence)

These are necessary conditions and linear count constraints on actual
two-band sources. They are not unrestricted depth bounds. The separate
[fixed Clifford-sum theorem](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) excludes
the general two-band indicator source at sufficiently small fixed error
and sufficiently large address size. The results below retain useful
native structure and explain why a simpler rank argument does not give
that conclusion.

All specified output and private dirty inputs are arbitrary. Every
identity includes external references, literal phases and actual
inverses. No supplied reflection is counted as an encoder.

## 1. The scalar-commuting parities lie in a Reed–Muller space

Let $`m=2^p`$ and use the actual one-band encoders

```math
U=U_B^\dagger U_A,\quad
A=U_AZ_aU_A^\dagger,\quad B=U_BZ_aU_B^\dagger,\quad
K=J^\dagger UJ.
```

Suppose the complete source error satisfies

```math
\|K-I_p/\sqrt2\|\le\epsilon\lt\frac{\sin(\pi/8)}8.
```

This is below one tenth, so the retained rounding gives exact flatness,
$`\{A,B\}=0`$, and the complete unitary $`L=\sqrt2K`$.
For an output parity $`u\in\mathbb F_2^m`$, write

```math
C_u=U_BZ_uU_B^\dagger,\qquad
S_B=\{u:AC_uA=\pm C_u\}.
```

The products $`C_uC_v=C_{u+v}`$ show that S_B is a binary linear
subspace. Its scalar commutation signs form a character. For each u in
S_B, the transported-probe theorem supplies an actual two-band,
one-clean circuit with exact signal return and logical action
$`L^\dagger Z_uL`$. Its error is at most $`2\sqrt2\epsilon`$
against the literal Boolean diagonal sign

```math
I_p^\dagger Z_uI_p=(-1)^{u\cdot Y+u_x}.
```

Regard u as the truth table of the address function $`f_u(x)=u_x`$.
If this function has degree greater than three, the retained
[two-layer scalar-Boolean floor](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction)
is violated. Therefore

```math
S_B\subseteq\mathrm{RM}(3,p),\qquad
\dim S_B\le d_3(p):=\sum_{j=0}^3\binom pj.
\qquad\text{(S1)}
```

Here RM(3,p) denotes the space of truth tables of Boolean polynomials
of degree at most three. The restriction concerns every output parity,
not just the original individual output wires.

Since the indicator is Hermitian, the accepted block of the actual
inverse U-dagger has the same error. Interchanging the encoders proves
the same inclusion for

```math
S_A=\{u:B(U_AZ_uU_A^\dagger)B=\pm U_AZ_uU_A^\dagger\}.
\qquad\text{(S2)}
```

## 2. The exact computational-mask rank requirement

Set $`R_B=U_B^\dagger A U_B`$ and expand it in the Hermitian Pauli
basis with its actual nonzero coefficients. For a supported label v,
let $`x_{\rm out}(v)`$ be its X mask on the specified m outputs.
Fix any supported label $`v_0`$ and define

```math
H_{\rm out}(R_B)
 =\mathrm{span}\{x_{\rm out}(v)+x_{\rm out}(v_0):
                                      v\in\mathrm{supp}(R_B)\}.
```

Membership in S_B is equivalent to $`Z_uR_BZ_u=\pm R_B`$.
Pauli orthogonality makes this exactly the condition that all supported
output masks have the same dot product with u. Consequently

```math
S_B=H_{\rm out}(R_B)^\perp,\qquad
\dim H_{\rm out}(R_B)\ge m-d_3(p).
\qquad\text{(S3)}
```

This is a rank requirement on the actual computational X masks. It
does not identify coupling rank with mask span.

To keep the native count visible, write $`U_A=E_AK_A`$ and
$`U_B=E_BK_B`$, where each E is a Clifford conjugate of an actual
physical T band and each K is Clifford. Let r_A count the E_A axes
active on $`K_AZ_aK_A^\dagger`$. Their independent span is L_A,
and A has Pauli support exactly one affine coset of L_A.

Let $`L_B^A`$ be the span of the original E_B axes that fail to
commute with A as an operator, and let $`t_B(A)`$ count those axes.
Their physical independence makes this also their span dimension.
Every other E_B factor commutes with A and cancels in
$`E_B^\dagger A E_B`$. Conjugation by the retained rotations can
only add labels from $`L_B^A`$. In the common physical coordinates,

```math
\mathrm{supp}(R_B)
 \subseteq v_0+F_{K_B}^{-1}(L_A+L_B^A),\qquad
\dim H_{\rm out}(R_B)
 \le\dim(L_A+L_B^A)\le r_A+t_B(A).
```

Here $`F_{K_B}`$ is the actual symplectic label map of K_B; the
choice of affine representative does not affect the span. The mirrored
argument gives the simultaneous native constraints

```math
r_A+t_B(A)\ge m-d_3(p),\qquad
r_B+t_A(B)\ge m-d_3(p).
\qquad\text{(S4)}
```

In particular the complete source count satisfies

```math
T_A+T_B\ge m-d_3(p),\qquad
T_A+T_B+r_A+r_B\ge2\bigl(m-d_3(p)\bigr).
\qquad\text{(S5)}
```

These are linear requirements, compatible with the proposed O(m)
source budget. They do not supply a source or a growing-depth bound.

## 3. A literal one-band family with no nonzero scalar parity

For a Hermitian Pauli P, let $`T_P`$ mean a literal Clifford conjugate
of the native T gate, so

```math
T_P=e^{i\pi/8}\exp(-i\pi P/8).
```

The scalar is retained. On a signal a and m arbitrary logical wires,
define the actual Clifford and encoder

```math
C=\left(\prod_j\mathrm{CX}_{j\to a}\right)
  \left(\prod_jH_j\right),\qquad
E=C\left(\prod_jT_j\right)C^\dagger
  =\prod_jT_{X_aX_j},\qquad U_A=I,\quad U_B=EH_a.
```

C maps the independent Z_j axes to the independent commuting axes
$`X_aX_j`$. Thus this is one physical band, with common forward and
actual-inverse bounds

```math
T_U=m,\qquad G_U\le4m+1,\qquad D_T(U)=1,
\qquad b_{\rm private}=0.
\qquad\text{(S6)}
```

Only H, CNOT and T/T-dagger occur. The clean-image reflections are
$`A=Z_a`$ and $`B=X_a`$. Every logical column of the accepted block
satisfies the literal identity

```math
\sqrt2\langle0_a|U_B^\dagger|0_a\rangle
  =\prod_j(H_jT_j^\dagger H_j).
\qquad\text{(S7)}
```

Indeed the left signal bra is plus, so each $`T_{X_aX_j}^\dagger`$
acts there as the corresponding $`T_{X_j}^\dagger`$. The native
common phases agree and are not discarded. This logical unitary is a
product of one-qubit gates, not an indicator.

The complete Pauli support is also explicit:

```math
U_B^\dagger Z_aU_B
 =X_a\prod_j\frac{I-iZ_aX_j}{\sqrt2}.
\qquad\text{(S8)}
```

Every subset term is nonzero and has a distinct output X mask. Hence
$`H_{\rm out}=\mathbb F_2^m`$ and $`S_B=\{0\}`$. Actual
one-band structure, flatness and linear count therefore do not force
even one nonzero scalar-commuting output parity. This family is already
inside the stopped one-band source grammar. It diagnoses a limitation
of the rank argument rather than reopening that route.

## 4. Both parity spaces can vanish for a complete identity transition

The same limitation remains if both directions and exact Boolean
involutivity are required. Define actual one-band encoders

```math
E_A=\prod_jT_{Z_aX_j},\qquad E_B=\prod_jT_{X_aX_j},\qquad
U_A=E_A,\qquad U_B=E_BH_a.
```

For E_A the Clifford $`(\prod_jH_j)(\prod_j\mathrm{CX}_{a\to j})`$
maps Z_j to $`Z_aX_j`$; for E_B use C from Section 3. Both axis
lists are independent. The source has

```math
A=Z_a,\quad B=X_a,\quad S_A=S_B=\{0\},\qquad
L=\left(\prod_jT_{X_j}^\dagger\right)
  \left(\prod_jT_{X_j}\right)=I.
\qquad\text{(S9)}
```

The forward support calculation is (S8). For the reverse direction,
$`U_A^\dagger B U_A`$ has the same displayed subset expansion,
so its output-mask span is full as well. The common native ledger is

```math
T_U=2m,\qquad G_U\le8m+1,\qquad D_T(U)\le2,
\qquad b_{\rm private}=0.
\qquad\text{(S10)}
```

No logical wire or dirty spectator is initialized. The complete logical
identity is a Boolean involution of degree zero. It does not have the
specified indicator action; this is precisely why it satisfies these
geometric conditions without contradicting (S1) or the source gap.

## 5. Four real traces and the separate gap theorem

For logical Hermitian Paulis P,Q, put $`C=U_BPU_B^\dagger`$ and
$`D=U_AQU_A^\dagger`$. Then C commutes with B and D with A.
Let $`\mathrm{tr}_{\rm phys}`$ and $`\mathrm{tr}_{\rm logical}`$
denote normalized traces. Exact flatness gives

```math
\begin{aligned}
\alpha&=\mathrm{tr}_{\rm logical}(L^\dagger PLQ)\\
 &=\mathrm{tr}_{\rm phys}\bigl[(I+A)(I+B)CD\bigr]\\
 &=\mathrm{tr}_{\rm phys}(CD)
  +\mathrm{tr}_{\rm phys}(C\,AD)
  +\mathrm{tr}_{\rm phys}(BC\,D)
  +\mathrm{tr}_{\rm phys}(BC\,AD).
\end{aligned}
\qquad\text{(S11)}
```

All four terms are real overlaps of Hermitian Clifford involutions:
AD and BC are Hermitian because of the displayed commutations. Each
term belongs to $`\{0\}\cup\{\pm2^{-r/2}:r\ge0\}`$.
If one term has magnitude one, equality in normalized Hilbert–Schmidt
Cauchy–Schwarz makes the corresponding two involutions equal up to
sign. Substitution makes the other three terms vanish, using

```math
\mathrm{tr}_{\rm phys}(A)=\mathrm{tr}_{\rm phys}(B)
=\mathrm{tr}_{\rm phys}(AB)=0.
```

Thus that case has
$`\alpha=\pm1`$.

The four scalar alphabets alone do not exclude an apparent pattern
such as one half plus one half minus a dimension-dependent small term.
The [fixed Clifford-sum source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
resolves this issue using actual Clifford-vector relations. Its direct
defect argument also removes the need for flatness and covers bounded
clean banks. That proof and its precise accuracy and size scope have
one home there; no unrestricted depth claim is inferred from (S11).

## 6. Exact receipts

[`test_simultaneous_output_parities.py`](../../tests/test_simultaneous_output_parities.py)
checks all physical and logical columns, actual inverses, literal
normalized transitions, the subset-support identity and every output
parity for small native members. It also checks the complete identity
family in both directions. These are finite mechanism checks, not a
proof of asymptotic impossibility. The support and Reed–Muller arguments
above are analytic.
