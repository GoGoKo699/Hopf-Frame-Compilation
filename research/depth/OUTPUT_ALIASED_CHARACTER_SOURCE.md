# Direct output readout and the weaker algebra-source contract

[Joint algebra source](FULL_BANK_CHARACTER_TRANSPORT.md) ·
[Paid quadratic query](FULL_BANK_CHARACTER_TRANSPORT.md#5-the-quadratic-inverse-has-a-paid-xor-query-not-a-paid-in-place-map) ·
[Balanced multiplier](BALANCED_ALGEBRA_PERMUTATION.md) ·
[Fast returned products](ALGEBRA_QUERY_EXTRACTION.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

The coefficient carrier need not remain disjoint from the output bank.
Its character reads exactly the nonconstant subset-zeta coordinates of
that bank. Identifying these physical wires removes the outer translation
echo: one character word, four source occurrences and only one additional
bank of carrier wires suffice. The source itself may leave arbitrary
input-dependent unitary action on its dirty helpers; the character word
cancels that action with the actual inverse.

These are complete source-interface improvements. They do not supply a
sublogarithmic native source. The paid algebra-key group operation is
constant-depth, whereas converting an ordinary coherent key into the
corresponding algebra element remains unpaid. The unrestricted depth
frontier and R96 are unchanged.

## 1. Identify the character coordinates with the output wires

Let $`b\ge1`$, $`m=2^b`$, $`r=b+1`$ and $`w=2^r-1=2m-1`$.
The b-bit address x and m-bit output Y are arbitrary. Write the output
labels as $`S\subseteq[b]`$. Use the square-zero coefficient algebra
and the retained full-input permutation

```math
R_r=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2),\qquad
F_r|c\rangle=|f(c)\rangle,\qquad
f(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1.
```

Let k be its quadratic inverse, and set

```math
\iota(x)=\sum_{i=1}^b x_i e_i,\qquad
A_x=f(\iota(x))=\prod_{i=1}^b(1+x_i e_i)-1,\qquad
\lambda(y)=[e_{[r]}]y.
```

On Y perform the binary subset-zeta transform, costing $`bm/2`$ CNOTs.
Its resulting physical values are

```math
z_S(Y)=\bigoplus_{u\subseteq S}Y_u.
```

For each nonempty $`S\subseteq[b]`$, identify the carrier wire
$`y_{[r]\setminus S}`$ with the existing output wire holding
$`z_S(Y)`$. This is a wire assignment, not a copy or a supplied value.
The $`m-1`$ assigned coordinates are distinct; they contain the dummy
generator r and exclude the full set $`[r]`$. Supply the other m carrier
coordinates as arbitrary dirty inputs. In particular, the top coordinate
$`y_{[r]}`$ and every singleton $`y_{\{i\}}`$ with $`i\le b`$ are
among those dirty wires. The output wire $`Y_\varnothing`$ stays outside
the carrier. There is no physical constant-coefficient wire.

The retained character word uses

```math
J_x=F_rX^{\iota(x)}F_r^\dagger,\qquad
C_x:\quad J_x,\ Z_{[r]},\ J_x^\dagger,\ Z_{[r]}
```

in chronological order. Its exact action is
$`C_x|y\rangle=(-1)^{\lambda(yA_x)}|y\rangle`$ on the entire carrier.
Only the assigned coordinates occur in its exponent:

```math
\begin{aligned}
\lambda(yA_x)
 &=\sum_{\varnothing\ne S\subseteq[b]}z_S(Y)\prod_{i\in S}x_i\\
 &=Y_x+Y_\varnothing.
\end{aligned}
```

The first equality follows by complementary subset multiplication; the
second is the algebraic normal form of the equality functions. Every
unassigned carrier coefficient may be arbitrary and remains present in
the source calls. The complete C word restores every carrier wire,
including the aliased output wires. Reverse the actual zeta transform,
then correct the missing constant character by $`Z_{Y_\varnothing}`$.
Conjugating the completed phase word by Hadamards on all Y wires gives

```math
I_b|x,Y,h\rangle=|x,Y\mathbin\oplus e_x,h\rangle.
```

The actual chronological word is therefore

```math
H_Y,\quad\mathrm{zeta}_Y,\quad C_x,\quad
\mathrm{zeta}_Y^\dagger,\quad Z_{Y_\varnothing},\quad H_Y.
```

Output motion inside the source is explicitly allowed. No output or
carrier coordinate is initialized, and every equality holds with
arbitrary external references.

## 2. Source helper return is unnecessary inside the character word

The sufficient ideal source need not be $`F_r\otimes I_H`$. It may be
any full unitary of the form

```math
V_0=\sum_c |f(c)\rangle\langle c|\otimes S_c,
```

where each $`S_c`$ is an arbitrary unitary on the whole dirty helper
space H. The helper action may depend on every coefficient input and
need not be computational-basis preserving. The carrier permutation is
literal; all input-dependent phases may be included in the unitaries
$`S_c`$. This is a full-space promise, not a statement about zero-filled
helpers or about the success probability of individual basis columns.

Use $`J_x=V_0X^{\iota(x)}V_0^\dagger`$. For $`y=f(c)`$, its action is

```math
|y\rangle|h\rangle\longmapsto
|f(c+\iota(x))\rangle
S_{c+\iota(x)}S_c^\dagger|h\rangle.
```

The principal permutation is still $`y\mapsto y+A_x+yA_x`$.
The top-coordinate Z reads only this principal value. In the chronological
word $`J_x,Z,J_x^\dagger,Z`$, the actual J adjoint cancels the helper
unitary exactly. The two scalar phases therefore give the same character
$`(-1)^{\lambda(yA_x)}`$, tensored with identity on H. No separate
return after either source occurrence is required.

Suppose the emitted source V obeys

```math
\|V-V_0\|_{\mathrm{op}}\le\eta
```

for such a specified ideal source, with common forward/actual-inverse
bounds $`T_V,G_V,d_V,h_V`$. Here h_V counts private dirty helper wires,
excluding the w coefficient inputs. Its actual adjoint obeys the same
norm bound relative to $`V_0^\dagger`$. The direct readout uses four
V/actual-adjoint occurrences, so unitary telescoping gives the complete
common ledger

```math
\begin{aligned}
T_I&\le4T_V,\\
G_I&\le4G_V+bm+2b+2m+6,\\
D_T(I)&\le4d_V,\\
h_I&\le m+h_V,\\
\varepsilon_I&\le4\eta,\qquad a_{\mathrm{clean}}=0.
\end{aligned}
```

The Clifford remainder includes the two zeta words, two b-CNOT address
translations, two top-coordinate Z gates, the two output-Hadamard
envelopes, and the additional constant-character Z correction. Each Z
uses $`S^2`$.
In the positive-S convention, adjacent inverse powers reduce literally
by $`S^{-2}=S^2`$, giving the same common allowance for the actual reverse.

The error includes every address, output, carrier and private-helper
return error, as well as arbitrary reference correlations. Intermediate
approximate calls need not return anything exactly. If V equals the
ideal fiber source exactly, the completed indicator returns every
extra wire exactly, despite the source's possible helper action.

## 3. The arbitrary-key correction has the same weaker return interface

The coherent correction has principal action

```math
R:\quad(z,y)\longmapsto(z,R_z(y)),\qquad
R_z(y)=y+f(z)+yf(z),\qquad R_zR_{z'}=R_{z+z'}.
```

Its key z has w arbitrary coefficient bits. An ideal source may also
apply any unitary $`S_{z,y}`$ to its separate private helper space.
Only its displayed principal permutation is prescribed. Let $`S_x`$
add $`\iota(x)`$ to z, at b CNOTs, and use the chronological echo

```math
R,\ S_x,\ R^\dagger,\ S_x^\dagger.
```

Its principal action restores z and applies $`J_x`$ to y. Its private
helper action can depend on z and y and need not return. The subsequent
character word uses this echo and its actual adjoint around the
top-coordinate Z; Section 2's cancellation then returns both the key
and every private helper. This argument includes all nonsingleton key
coefficients rather than replacing them by formal zeros.

If a native source approximates that full ideal fiber unitary within
eta, with common bounds $`T_R,G_R,d_R,h_R`$, the direct output readout
uses four R/actual-adjoint occurrences and gives

```math
\begin{aligned}
T_I&\le4T_R,\\
G_I&\le4G_R+bm+4b+2m+6,\\
D_T(I)&\le4d_R,\\
h_I&\le3m-1+h_R,\\
\varepsilon_I&\le4\eta,\qquad a_{\mathrm{clean}}=0.
\end{aligned}
```

Width counts the m unassigned carrier wires and the entire
$`2m-1`$-wire key. The four key translations explain the $`4b`$ term.
As before, the complete norm estimate does not assume exact intermediate
key or helper return. The exact fiber source gives exact final return;
an approximate source gives the full stated error bound.

## 4. Fusion and the fixed-accuracy frame consequence

Apply the retained root-fusion reduction to either complete indicator
above. It uses eight source-indicator/actual-inverse occurrences. The
combined source contribution is therefore 32 V calls, or 32 R calls,
and full error at most $`32\eta`$. With the actual large prefix
$`b_0=p-s`$ and exact suffix
$`s=\lceil a\log_2(p+2)\rceil`$ chosen by fusion, its depth is at most

```math
32d_V(b_0+1)+O_\tau(\log\log(p+4))
```

for the coefficient source, or the corresponding expression with d_R.
Request $`\eta=\tau/32`$ for final indicator tolerance tau. The source
still needs fixed polynomial address overhead in T-count and dirty
width, and at most quadratic bank size times polynomial overhead in
Clifford count. Source helpers with arbitrary fiber action are permitted
because the *completed* indicator satisfies fusion's full returned-work
promise. The exact suffix and short-address fallbacks are unchanged.

For the fixed-L frame transfer, every completed word query uses at most
four normalized indicators. Its added error is at most $`128\eta`$.
With the retained dimension-independent query bound $`K_*`$, requesting

```math
\eta_L=\frac{2^{-L}}{512K_*}
```

keeps the added frame error at most $`2^{-L}/4`$. The retained real and
complex budgets then give complete initialized-isometry error at most
$`2^{-L}`$, with the same literal one-clean and dirty allocations.
This is the same fixed-L conditional transfer, with a smaller source
constant. It is not a precision-uniform improvement when L grows.
The known logarithmic-depth source still retains the old depth order.

## 5. The paired transform and the coherent correction are equivalent interfaces

Let $`Q_k(y;z)`$ be the paid complete quadratic XOR
$`(z,y)\mapsto(z+k(y),y)`$. Its helper return is exact and its depth
is at most 40. On two arbitrary w-bit registers, execute chronologically

```math
Q_k(y;z),\quad R(z;y),\quad Q_k(y;z),\quad
\mathrm{SWAP}_{z,y}.
```

After the first query the key is $`z+k(y)`$. The correction sends the
carrier to

```math
R_{z+k(y)}(y)
 =f(k(y)+z+k(y))=f(z).
```

The second quadratic query changes the key to
$`z+k(y)+k(f(z))=k(y)`$. Thus the final swap gives

```math
P=F_r\otimes F_r^\dagger:\quad(z,y)\longmapsto(f(z),k(y)).
```

This is a full-input identity, including every old value of both banks.
Conversely, writing Q for this same quadratic query and W for SWAP,
the operator identity $`P=WQRQ`$ gives $`R=QWPQ`$. Therefore the
chronological reverse reduction is $`Q,P,W,Q`$. If R or P has arbitrary
private fiber action, the reduction has exactly the displayed principal
action with another arbitrary private fiber action. These are the weaker
source promises of Sections 2 and 3, not independently returned sources.

For either direction, the added common costs are

```math
\Delta T\le2T_{Q_k},\qquad
\Delta G\le2G_{Q_k}+3w,\qquad
\Delta D_T\le80.
```

Use a private pool of size $`\max\{h_{Q_k},h_{\mathrm{given}}\}`$,
disjoint from the two live coefficient banks. The exact Q word accepts
every arbitrary pool input, so reuse remains valid even when the given
source only returns that pool approximately or carries ideal fiber
action. With Q exact, a supplied full-operator error eta remains eta
under either reduction. Every inverse is the actual reversed word.

As a source for Section 2, the second coefficient bank occupies w
additional dirty wires; its k action is a permissible fixed helper
unitary. It returns in the completed character. It does **not** return
after the paired source alone. The identity therefore does not assert
a conversion from a quadratic XOR to a standalone in-place F, nor from
$`F\otimes F^\dagger`$ to $`F\otimes I`$ for free. It shows that a
shallow paired transform is an equivalent sufficient target for this
coherent-correction route, up to paid constant-depth overhead.

## 6. A paid algebra-key multiplication does not yet convert the key

The algebra-key group operation

```math
B:\quad(a,y)\longmapsto(a,a+y+ay)
```

already has a complete constant-depth implementation. Here a and y are
arbitrary augmentation-ideal coefficient words. Use the retained
generator-split multiplier with one bank to multiply y by $`1+a`$,
then add a by w CNOTs. For $`r\ge2`$ set

```math
\rho_r=\binom{r+1}{2}2^{r-1},\qquad
s_r=r(r-1)2^{r-2},\qquad m_0=2^{r-1}.
```

Substituting the retained fast subset-product word into that single
multiplier gives the common conservative bounds

```math
\begin{aligned}
T_B&\le336\rho_r,\\
G_B&\le978\rho_r+w,\\
D_T(B)&\le200,\\
h_B&\le9\rho_r,\qquad\varepsilon_B=0.
\end{aligned}
```

More directly, its six product instances and the original copies give
$`G_B\le576\rho_r+384s_r+18m_0-2+w`$ and
$`h_B\le m_0+6\rho_r+2s_r`$. These are safe upper counts after
omitting structural-zero coordinates. At r=1 the operation is one CNOT.

There is no supplied constant coefficient in this application. Splitting
one generator writes the carrier as $`x+e_r y'`$ and the key as
$`a_0+e_r v`$. The constant coefficient of x is structurally zero.
During the retained two shears and two copies it stays zero: multiplication
by $`a_0`$ has zero constant term, and copies from the absent x coordinate
are omitted. The subsequent correction also only uses that zero
coordinate as an omitted source. All other carrier and private mask
coordinates are physical arbitrary wires. The retained exact reverse
returns every dirty mask and correction bank.

This key is a, whereas the required coherent correction uses
$`a=f(z)`$. They must not be identified for arbitrary z. Ordinary
Hadamards on the entire carrier do not supply that conversion. Let
$`M_a`$ be multiplication by a on the augmentation-ideal coefficient
space. Since $`M_a^2=0`$ and $`M_a a=0`$, direct binary Fourier
transformation gives

```math
H_y B H_y|a,t\rangle
 =(-1)^{t\cdot a}|a,(I+M_a^{\mathsf T})t\rangle.
```

This is another quadratic principal permutation with a bilinear phase.
It contains no evaluation of f. The group characters instead have the
kernel $`(-1)^{t\cdot k(y)}`$. The corresponding single-register
group Fourier transform is $`H F_r^\dagger`$. Although its two-register
phase oracle is paid by the quadratic query, that phase oracle is not
the single-register transform. Treating it as such is precisely an
unpaid change of representation, analogous to treating an out-of-place
XOR as an in-place permutation.

## 7. The shallow Fourier kernel is not a returned scalar source

The quadratic query exposes the desired high-degree character in a
different way, but its amplitude and clean allocation must be retained.
Put $`d=2^w`$, use two arbitrary w-bit registers A and B, and define

```math
Q=Q_k(A;B),\qquad M=H_A Q H_A.
```

Its exact matrix kernel is

```math
\langle a',b'|M|a,b\rangle
 =\frac1d(-1)^{(a+a')\cdot f(b+b')}.
```

To prove this, sum over the intermediate A value y. The query requires
$`k(y)=b+b'`$, which has exactly the solution $`y=f(b+b')`$.
The two Hadamard factors supply the displayed coefficient $`1/d`$.

For arbitrary disjoint coherent keys alpha and beta, form

```math
W_{\alpha,\beta}=X_A^\alpha X_B^\beta M,\qquad
p=(-1)^{\alpha\cdot f(\beta)}.
```

Every diagonal matrix entry of W is $`p/d`$. This does not make W a
scalar on arbitrary A and B. M commutes with every X translation, as
its kernel depends only on the two input-output differences. Hence W
is a Hermitian involution. Its dimension is $`d^2`$, and its trace is
$`pd`$. For $`w\ge1`$ both eigenvalues occur, giving

```math
\|W_{\alpha,\beta}-(p/d)I\|_{\mathrm{op}}=1+1/d.
```

Thus the diagonal or trace character does not give an approximate
full-operator scalar source. In particular, a one-clean Hadamard test
compresses to $`(I+W)/2`$ and does not replace W by its trace.

### An explicit weak block and the cost of this particular scalarization

There is an exact diagonal-twirl identity on the two carrier registers:

```math
2^{-2w}\sum_{s\in\mathbb F_2^{2w}}Z^s W_{\alpha,\beta}Z^s
 =\frac{p}{d}I.
```

It has a literal unitary block encoding with $`2w`$ **clean** label
bits: Hadamards on those labels, their controlled Z word into A and B,
W, the same controlled Z word, and final label Hadamards. Compression
at label zero gives exactly the displayed block. All A and B inputs
are arbitrary. With the retained query implementation, its common
source ledger is

```math
T_{\mathrm{block}}=T_{Q_k},\qquad
d_{\mathrm{block}}\le40,\qquad
G_{\mathrm{block}}\le G_{Q_k}+20w,
```

with $`2w`$ clean labels, $`2w+h_{Q_k}`$ dirty work excluding the two
declared keys, and compressed amplitude $`2^{-w}`$. The Clifford
allowance charges both A-Hadamard banks, the coherent X translations,
all label Hadamards and both controlled-Z banks; each CZ is the literal
H-CNOT-H word. This is a compressed block, not a full returned source.

The exponential label cover is necessary for the stated diagonal-twirl
representation at any fixed useful relative error. Let $`N=d^2`$ and
$`W_s=Z^sWZ^s`$. In the normalized Hilbert-Schmidt inner product,
the N matrices $`W_s`$ are orthonormal. Indeed every entry of W has
modulus $`1/d`$, so the inner product between $`W_s`$ and $`W_t`$ is

```math
N^{-2}\sum_{i,j\in\mathbb F_2^{2w}}
             (-1)^{(s+t)\cdot(i+j)}=\mathbf1_{s=t}.
```

The inner product of each $`W_s`$ with $`(p/d)I`$ is $`1/N`$,
and the latter matrix has squared normalized norm $`1/N`$. Projection
onto any K of these orthonormal matrices therefore leaves squared error

```math
\frac1N\left(1-\frac KN\right).
```

Since operator norm dominates normalized Hilbert-Schmidt norm, an
arbitrary linear combination of K distinct mask conjugates with
operator error at most $`\epsilon/d`$, for $`0\le\epsilon\lt1`$,
must satisfy

```math
K\ge(1-\epsilon^2)N.
```

For coefficients whose sum is one, squared error is
$`\sum_s|c_s|^2-1/N\ge1/K-1/N`$, yielding the sharper
$`K\ge N/(1+\epsilon^2)`$. A single clean selector encoding K such
terms consequently needs at least $`2w-O_\epsilon(1)`$ label bits;
the exact twirl needs $`2w`$.

The exhibited block therefore does not meet the bounded-clean,
amplitude-$`1/\sqrt2`$ source route. Standard promotion of its label
bank incurs an exponential sector factor; standard amplification of
its $`2^{-w}`$ amplitude incurs an exponential occurrence cost. Neither
fits the polynomial-overhead complete-source window. These conclusions
concern this mask-twirl representation and these repairs. The mask count
is not an unrestricted native-depth lower bound, and the argument does
not exclude a different coherent scalarization or whole-word synthesis.

## 8. What remains open after this revision

The direct alias removes a needless outer echo and permits a broader
source with arbitrary unitary helper fibers. A source need not solve
its own helper-return problem before the completed character is formed.
A shallow paired transform, a coherent correction, or a direct completed
extraction can use that freedom. The paid quadratic XOR and algebra-key
multiplier are useful exact components, but their mere composition in
the retained classical architecture still retains the old cost.

The existing inverse-degree argument already permits arbitrary classical
helper garbage. An exact full-input classical map with principal output
k has inverse principal output f of degree r. A word of D completed
stages with inverse coordinate degree at most delta therefore needs
$`\delta^D\ge r`$. Likewise the coherent correction has degree r on
its full principal inputs. These are bounds on the specified completed
classical stages, not on unrestricted native T layers. Arbitrary
nonmonomial interference or genuinely quantum conditional helper
unitaries are not covered by that degree argument.

The next constructive question is consequently a sublogarithmic native
implementation of one of these weaker complete source interfaces, or
of the whole extraction, within the relaxed polynomial-overhead budgets.
No source construction or unrestricted obstruction is supplied here.
The original physical-pair lower scout remains an independent question;
none of these interface reductions bounds its unrestricted coupling rank.

## 9. Focused composition evidence

The [single exact fixture](../../tests/test_output_aliased_character.py)
uses b=2, all four arbitrary output bits, four additional arbitrary
carrier bits and one arbitrary private helper, on eleven physical wires.
It checks all 2048 phase columns of the direct construction and its
actual reverse. The zeta transform is expanded into its four actual
CNOTs. The fixture also verifies the physical alias count, the exclusion
of $`Y_\varnothing`$ from the carrier, the four source occurrences and the 26-Clifford
remainder including the analytic output-Hadamard envelopes.

The source deliberately applies an input-dependent quantum helper word:
a native controlled-H followed by a controlled Z. A checked source
column leaves that helper in a coherent superposition. Thus the test
does not merely permute classical helper values or assume source-level
helper return. The completed character restores every physical wire
with the literal selected-output phase; its actual reversed word does
the same. Output-Hadamard conjugation then gives the XOR indicator
analytically, including arbitrary reference entanglement.

F and its inverse are the retained exact coefficient permutations in
this fixture. Their native compilation and all-size cost bounds remain
in the existing proof homes. This is a composition and quantum-helper
return witness, not a new native HRT emitter or a numerical proof of
the asymptotic source reductions, Fourier-twirl bound, or unrestricted
depth frontier.
