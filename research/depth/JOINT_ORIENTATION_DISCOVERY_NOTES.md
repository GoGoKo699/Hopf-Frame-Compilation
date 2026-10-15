# Algebraic side results from the joint orientation session

[R85 proof](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) · [Session result](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)

These side results are **not needed for the R85 proof**. They preserve two
exact calculations from the joint discovery session: an isolated two-level
fusion and a restriction on five fully dirty XOR shears. Neither is a
generic native lower bound, and neither supplies R85's geometric recovery
mechanism. No additional local-result label is assigned here.

## 1. Exact two-level fusion

Let $`J=-iY=\begin{pmatrix}0&-1\\1&0\end{pmatrix}`$ and
$`Q=|+\rangle\langle+|`$. The two visible branch bits span four dimensions;
deeper suffix bits are fixed to their common uniform vector. For a parent
angle $`\delta`$ relative to balance and child angles
$`\epsilon_0,\epsilon_1`$, the chronological parent-then-children operator is

```math
U=\mathrm{diag}(e^{\epsilon_0J},e^{\epsilon_1J})
  e^{\delta J\otimes Q}.
```

Define the two triples

```math
A_1=J\otimes I,\quad A_2=Z\otimes J,\quad A_3=X\otimes J,
\qquad
B_1=I\otimes J,\quad B_2=J\otimes Z,\quad B_3=J\otimes X.
```

Direct Pauli multiplication gives $`[A_j,B_k]=0`$ for every j,k.
Each triple squares to minus identity and satisfies its cyclic quaternion
relations, for example $`A_1A_2=A_3`$ and $`B_1B_2=B_3`$.
With

```math
\beta=(\epsilon_0+\epsilon_1)/2,\qquad
\gamma=(\epsilon_0-\epsilon_1)/2,\qquad d=\delta/2,
```

the two exponents are $`\beta B_1+\gamma A_2`$ and
$`d(A_1+B_3)`$. Their cross-triple commutation therefore proves the
literal identity

```math
U=\bigl(e^{\gamma A_2}e^{dA_1}\bigr)
  \bigl(e^{\beta B_1}e^{dB_3}\bigr),
```

whose two parenthesized factors commute. If the child angles have equal
magnitude and differ only in sign, one of beta and gamma vanishes. Except
at zero magnitude, exactly one vanishes.

There is a fixed Clifford basis realizing this separation. Products below
act right to left. Set

```math
M=\mathrm{CNOT}_{1\to2}(H\otimes I)\mathrm{diag}(1,i,i,1),
\qquad \mathrm{diag}(1,i,i,1)=(S\otimes S)\mathrm{CZ}.
```

Its Pauli conjugations are

```math
\begin{aligned}
MA_1M^\dagger&=iI\otimes Z,&
MA_2M^\dagger&=iI\otimes X,&
MA_3M^\dagger&=-iI\otimes Y,\\
MB_1M^\dagger&=iX\otimes I,&
MB_2M^\dagger&=iZ\otimes I,&
MB_3M^\dagger&=iY\otimes I.
\end{aligned}
```

Substitution proves, with literal phase,

```math
MUM^\dagger=
\bigl(e^{i\beta X}e^{idY}\bigr)\otimes
\bigl(e^{i\gamma X}e^{idZ}\bigr).
```

This is an exact isolated two-level identity. At two visible levels the
child orientation is a Boolean function of only one visible parent bit,
so it is already affine. External preserved prefix labels may select the
angles, but their queries and the uniform-suffix activity predicate still
require a native implementation and its resource accounting.

The identity does not close recursively on the original anchors. Indeed,

```math
M|++\rangle=e^{i\pi/4}(S\otimes I)\mathrm{CZ}|++\rangle
=\frac{e^{i\pi/4}}2(1,1,i,-i)^{\mathsf T}.
```

The two coefficient rows are orthogonal with equal norm, so either
one-qubit reduced state is $`I/2`$. Local one-qubit corrections preserve
these Schmidt coefficients and cannot restore a product anchor. An
entangling correction need not preserve the separated one-qubit factors.
Computational prefixes likewise become Bell-sector selectors, while the
new rotation axes need not be the original common real Y axis. Thus this
calculation supplies neither a closure rule for another fusion round nor
a simultaneous implementation of nested, changing-prefix gauges.

## 2. Five fully dirty XOR shears force an exponent-four target

Let V,W be finite binary vector spaces, with addition written as plus.
For arbitrary functions define

```math
A_f(x,y)=(x,y+f(x)),\qquad B_g(x,y)=(x+g(y),y).
```

Suppose five alternating shears implement
$`P_F(x,y)=(F(x),y)`$ for **every** $`x\in V,y\in W`$.
Then $`F^4=I`$, for either starting orientation. This statement requires
all of Y to be arbitrary; it does not cover an initialized coordinate.

First consider the chronological word
$`B_g,A_f,B_h,A_k,B_l`$. Put

```math
u=x+g(y),\qquad t=y+f(u),\qquad w=u+h(t).
```

As x,y vary, u,t range independently over V,W. Exact return of Y gives
$`k(u+h(t))=f(u)`$. Write $`h_0=h(0)`$ and

```math
K=\mathrm{span}\{h(t)+h_0:t\in W\}\subseteq V.
```

Comparing the equation at t and at zero shows that k is invariant under
each generator of K, hence under K. Since $`f(u)=k(u+h_0)`$, f has the
same invariance. The output equation is

```math
F(x)=x+g(y)+h\bigl(y+f(x+g(y))\bigr)+l(y).
```

For every $`v\in K`$ it implies $`F(x+v)=F(x)+v`$.
Also $`h(W)\subseteq h_0+K`$; fixing $`y=0`$ shows that
$`F(x)+x`$ lies in one fixed coset $`c+K`$. Consequently
$`G=F^2`$ is K-equivariant and $`G(x)+x\in K`$.
Writing $`G(x)=x+v`$ gives
$`G^2(x)=G(x+v)=G(x)+v=x`$, proving $`F^4=I`$.

For the chronological A-B-A-B-A case, write its operator as
$`A_lV_4=P_F`$, where $`V_4`$ is the initial four-shear word.
Each shear is an involution, and

```math
A_lP_F=P_FA_{l\circ F}.
```

Thus $`V_4A_{l\circ F}=P_F`$. The rightmost A shear executes first
and merges with the initial A shear of $`V_4`$. This produces a
four-shear chronological A-B-A-B word. Prepending the identity shear
$`B_0`$ puts it in the B-A-B-A-B case just proved. This absorption uses
the identity on all of Y and does not preserve a promised initialized-Y
subspace.

## 3. Six shears already escape that exponent

For $`X=(a,b)`$ and an arbitrary dirty bit d, execute the following six
assignments chronologically; each right-hand side uses current values:

```math
\begin{aligned}
d&\leftarrow d+a,& b&\leftarrow b+d,& d&\leftarrow d+a+b,\\
a&\leftarrow a+d,& d&\leftarrow d+b,& b&\leftarrow b+d.
\end{aligned}
```

The first three assignments give $`(a,a+b+d,a+b)`$. The fourth and
fifth give $`(b,a+b+d,d)`$, and the sixth gives

```math
(a,b,d)\longmapsto(b,a+b,d).
```

The dirty bit returns on every input. The target has square
$`(a,b)\mapsto(a+b,a)`$ and cube identity, so its order is three and
it violates the five-shear conclusion. These are linear XOR operations,
hence already Clifford operations. The example only marks the limit of
the five-shear algebra; it gives no construction for arbitrary nonlinear
permutations and no native depth obstruction.
