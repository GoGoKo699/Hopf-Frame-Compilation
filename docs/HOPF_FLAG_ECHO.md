# Exact flag-echo audit for amplified Hopf sources

[Error accumulation](HOPF_ERROR_ACCUMULATION.md) · [Operator source](OPERATOR_SOURCE_COMPILER.md) · [Current frontier](OPEN_PROBLEM.md)

The proposed two-half-angle reflection echo is exactly the square of the
amplified source. It suppresses radial leakage near zero, but retains
first-order leakage at ordinary target angles. The other three diagonal
Pauli flag echoes also have exact error formulas. One has a useful exact
equal-mask identity, without supplying uniform cancellation for arbitrary
angles. These statements concern the actual source algebra on arbitrary
dirty inputs; they do not change the frame count/depth bounds.

## 1. Literal source and composite conventions

Work on one addressed row, with target operator
$`K=XZ=-iY`$, so $`K^2=-I`$. Let M, N_c, and N_s be the
Hermitian Majorana-vector involutions of the source chapter. Write

```math
\{M,N_c\}=2cI,\qquad \{M,N_s\}=2sI,\qquad
D_c=(MN_c-N_cM)/2,\quad D_s=(MN_s-N_sM)/2.
```

The real encoded coefficients satisfy $`|c|,|s|\le1`$. The exact
one-flag source operators and the two-flag query are

```math
S_c=cI+X_aD_c,\qquad S_s=sI+X_aD_s,\qquad
Q=H_b\,\mathrm{diag}_b(S_c,KS_s)\,H_b.
```

Here $`D_f^\dagger=-D_f`$ and
$`D_f^2=-(1-f^2)I`$. Target operators commute with the dirty core.
Let J initialize only flags $`(b,a)=(0,0)`$, and define

```math
P=JJ^\dagger,\qquad R=I-2P,\qquad
A=-QRQ^\dagger RQ,\qquad B=J^\dagger QJ=(cI+sK)/2,
\qquad \delta=1-c^2-s^2.
```

The inverse in A is the actual inverse of Q. The displayed scalar minus
is retained. The normalized physical source $`\widehat A=-RA`$
has the exact proposed echo

```math
R\widehat A R\widehat A=A^2.
```

Thus its second call receives the first call's leaked flags and dirty
core. No reinitialization occurs between the calls. The identities below
are operator identities tensored with arbitrary spectators and references.

## 2. Exact error of the square

Expansion of the two reflections gives

```math
A=-Q+2PQ+2QP-4QJB^\dagger J^\dagger Q,\qquad
AJ=\delta QJ+2JB.
```

Put $`D=J^\dagger Q^2J`$. A second application yields the full
initialized-isometry identity

```math
A^2J=\delta[-Q^2J+2JD+4QJ(B-B^\dagger D)]+4JB^2.
```

Since $`B^\dagger B=BB^\dagger=(1-\delta)I/4`$, its accepted
compression is

```math
J^\dagger A^2J=\delta^2D+4(1+\delta)B^2.
```

The source structure determines D, rather than leaving an arbitrary
unitary dilation. In $`Q^2=H_b\mathrm{diag}(S_c^2,-S_s^2)H_b`$,
the a-zero compression of $`S_f^2`$ is $`(2f^2-1)I`$. Hence

```math
D=tI,\qquad t=c^2-s^2,\qquad u=2cs,
\qquad J^\dagger A^2J=t(1+\delta+\delta^2)I+u(1+\delta)K.
```

Assume $`q=c^2+s^2=1-\delta\gt0`$ and define the encoded polar
rotation

```math
U=(cI+sK)/\sqrt q=R_y(\phi),\qquad \phi=\mathrm{atan2}(s,c).
```

Its doubled target is $`U^2=\cos(2\phi)I+\sin(2\phi)K`$. Thus

```math
J^\dagger A^2J=(1-\delta^3)\cos(2\phi)I
 +(1-\delta^2)\sin(2\phi)K.
```

Unitarity of A gives the error on all dirty input columns, not merely its
accepted angle. Explicitly,

```math
(A^2J-JU^2)^\dagger(A^2J-JU^2)
=2\delta^2[\sin^2(2\phi)+\delta\cos^2(2\phi)]I,
```

so

```math
\boxed{\|A^2J-JU^2\|
=\sqrt2|\delta|\sqrt{\sin^2(2\phi)+\delta\cos^2(2\phi)}.}
```

Negative radial defects are included. If $`\delta\lt0`$, the physical
bounds on c and s imply
$`|\cos(2\phi)|\le(1+\delta)/(1-\delta)`$. The radicand is
therefore nonnegative; it is at least
$`(-3\delta-\delta^2)/(1-\delta)`$. Replacing delta by its
absolute value inside the formula would change the answer.

### A fixed-angle witness and a genuine near-zero improvement

For $`c=s=r\ne0`$, the polar target is exactly $`U^2=K`$.
The accepted block is $`(1-\delta^2)K`$, while direct compression
onto rejected flags gives

```math
\langle b=1,a=0|A^2J
=\delta[(1-t^2)I+utK]=\delta I.
```

Consequently the full error is exactly $`\sqrt2|\delta|`$, and
the rejected norm is $`|\delta|\sqrt{2-\delta^2}`$.
This has no polar-angle bias to blame for the surviving first-order error.
For example, the literal dyadic source $`c=s=3/4`$ has
$`\delta=-1/8`$ and error $`\sqrt2/8`$. Nearest dyadic
approximations to $`1/\sqrt2`$ give arbitrarily fine nonzero-defect
examples at the same target angle.

There is a scoped improvement near zero. For literal rows with c=1 and
small dyadic s, $`\delta=-s^2`$, and the exact error becomes

```math
\|A^2J-JU^2\|=|s|^3\sqrt{\frac{2(3-s^2)}{1+s^2}}.
```

On an axis s=0 or c=0, the physical defect is nonnegative and the error
is $`\sqrt2\delta^{3/2}`$. These special regimes do not supply a
uniform quadratic radial-error bound at other angles.

## 3. All four diagonal Pauli flag echoes

More generally let C be a flag unitary satisfying $`CJ=J`$, and form
$`F_C=ACAC^\dagger`$. The inverse on the right is literal.
Define $`D_C=J^\dagger QCQJ`$. The same expansion gives

```math
F_CJ=\delta[-QCQJ+2JD_C+4QJ(B-B^\dagger D_C)]+4JB^2,
\qquad J^\dagger F_CJ=(1-\delta^2)U^2+\delta^2D_C.
```

No commutation assumption on D_C was used. In particular,

```math
(F_CJ-JU^2)^\dagger(F_CJ-JU^2)
=\delta^2[2I-U^{\dagger2}D_C-D_C^\dagger U^2].
```

For the common-source masks, put
$`gI=\{N_c,N_s\}/2`$. Their Majorana-vector structure makes g a
real scalar with $`|g|\le1`$, and

```math
\{D_c,D_s\}=2(cs-g)I.
```

This identity evaluates all four choices. Z_a and Z_b act on the
indicated flags, and each fixes their initialized state.

| C | Exact D_C | Exact full error $`\lVert F_CJ-JU^2\rVert`$ |
|---|---|---|
| I | $`tI`$ | $`\sqrt2\lvert\delta\rvert\sqrt{\sin^2(2\phi)+\delta\cos^2(2\phi)}`$ |
| $`Z_a`$ | 0 | $`\sqrt2\lvert\delta\rvert`$ |
| $`Z_b`$ | $`(u-g)K`$ | $`\sqrt2\lvert\delta\rvert\sqrt{1-(u-g)\sin(2\phi)}`$ |
| $`Z_aZ_b`$ | $`gK`$ | $`\sqrt2\lvert\delta\rvert\sqrt{1-g\sin(2\phi)}`$ |

For clarity, the Z_a calculation uses
$`S_fZ_aS_f=Z_a`$. Conjugating Z_b by H_b interchanges the two
query branches, so its compression involves
$`\{D_c,D_s\}`$ with a plus sign. Including Z_a reverses that
sign, producing g instead of $`u-g`$. Also $`|u-g|\le1`$:
D_C compresses the unitary QCQ. All table radicands are nonnegative.

When $`N_c=N_s`$, we have c=s and g=1. In fact the
$`Z_aZ_b`$ composite then implements K on the entire physical space.
Let $`Q_{s,c}`$ denote the query with the two full source masks
exchanged, and similarly for $`A_{s,c}`$. For $`C=Z_aZ_b`$,

```math
CQ_{c,s}C=KQ_{s,c}^\dagger,\qquad
CA_{c,s}C=KA_{s,c}^\dagger,\qquad
F_C=KA_{c,s}A_{s,c}^\dagger.
```

The first identity uses $`Z_aS_fZ_a=S_f^\dagger`$ and the
branch exchange by Z_b; the second uses commutation of K with Q and R.
Identical masks make the final product exactly K, including its phase
and all rejected flag sectors. Equality of c and s alone does not imply
that stronger full-space identity. At other angles the error obeys

```math
\|F_{Z_aZ_b}J-JU^2\|
\ge\sqrt{2(1-|\sin(2\phi)|)}\,|\delta|.
```

The same lower bound holds for Z_b. Thus their equal-angle cancellation
does not extend uniformly away from the doubled quarter-turn targets.
Z_a never suppresses the first-order defect. This classifies these four
echoes, not all Clifford flag conjugations, phase sequences, different
source pairs, or longer composite words.

## 4. What a compiler improvement would still require

The formulas compare with the encoded polar angle. A requested rotation
$`R_y(\theta)`$ additionally incurs the angle error
$`2|\sin(\phi-\theta/2)|`$ between $`U^2`$ and that target.
The equal-row witness has exactly zero such error. In general it cannot
be discarded when certifying a modified layer or complete frame.

The literal composite uses two calls to A, or six Q/Q-adjoint occurrences
before further simplification. The flag reflections and Pauli choices
are Clifford operations. Lookup banks still obey their existing exact
return contracts; source-core disturbance and flag leakage are included
in the full norm above. No supplied catalyst, reset, postselection, or
extra initialized work is assumed. Actual reversal restores all inputs.

This closes the uniform radial-cancellation question for the square and
these three Pauli variants. It does not prove a lower bound on general
source precision, fixed-accuracy frame T-depth, or the best longer word.
The separate [phased radial filter](HOPF_RADIAL_FILTER.md) treats a longer
sequence outside this four-choice audit. Its full isometry error, source
angle bias, selective-phase synthesis, and extra calls must all be charged
when composing a frame; the obstructions here do not rule out that route.

The [bounded checks](../tests/test_flag_echo.py) use small common
Clifford-algebra representations of actual dyadic source masks, preserving
both flags, the dirty core, phases, and actual inverses. They check the
square's exact formulas, its fixed-angle obstruction, and its near-zero
improvement, as well as the three Pauli variants and the exact full-space
equal-mask exception. These are algebraic matrix fixtures, not an emitted
general native composite or a numerical proof of the uniform formulas above.
