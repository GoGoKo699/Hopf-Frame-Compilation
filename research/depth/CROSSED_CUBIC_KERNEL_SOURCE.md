# Crossed cubic kernels and complete observable transport

[Observable contract](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Output aliases](OUTPUT_ALIASED_CHARACTER_SOURCE.md) ·
[Paid quadratic query](FULL_BANK_CHARACTER_TRANSPORT.md#5-the-quadratic-inverse-has-a-paid-xor-query-not-a-paid-in-place-map) ·
[Earlier key-mixing screens](KEY_MIXING_OBSERVABLE_SCREENS.md) ·
[Involutory Gaussian loops](INVOLUTORY_GAUSS_SOURCE_LIMITS.md)

A quadratic permutation's high-degree inverse occurs in the exact
kernel of a paid non-Gaussian Fourier word. Two crossed kernels give
an explicit collision condition for complete observable transport.
The concrete cubic candidate below has zero clean work and at most
336 native T layers, with simultaneous polynomial-overhead budgets,
but fails the required full-space identity.

The reusable obstructions concern actual protected sectors and, for a
broader middle class, an unbalanced top derivative. They permit arbitrary
dirty inputs and references. They do not exclude general crossed
interference, the direct observable contract, or unrestricted native
depth. The upper and lower frontier is unchanged.

## 1. Physical inputs and the aliased target

Put $`r=p+1`$, $`m=2^p`$, $`w=2^r-1=2m-1`$ and $`T=[r]`$,
with $`p\ge1`$. Bank A has one qubit for each nonempty subset of T
and represents the augmentation ideal of

```math
R_r=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2).
```

There is no constant-coefficient qubit or initialized coefficient bank.
Use the retained quadratic inverse

```math
f(a)=\prod_{\varnothing\ne S\subseteq T}(1+a_Se_S)-1,
\qquad k(a)=a+\gamma_2(a),
\qquad k(a+u)+k(a)=k(u)+au.
```

After the fixed subset-zeta transform of Y, identify the actual output
wire holding $`z_S(Y)`$ with $`A_{T\setminus S}`$ for every
$`S\subseteq[p]`$, including the empty subset. The remaining m minus
one A coordinates, the full w-bit bank B, and all query helpers are
arbitrary private inputs. The [alias proof](OUTPUT_ALIASED_CHARACTER_SOURCE.md#1-identify-the-character-coordinates-with-the-output-wires)
gives the required observable in these coordinates:

```math
\Omega_x=Z_{A,T}
\prod_{\varnothing\ne S\subseteq[p]}
Z_{A,T\setminus S}^{\prod_{i\in S}x_i}.             \qquad(1)
```

Undoing the fixed zeta transform gives $`Z_{Y_x}\otimes I_H`$.
In particular, at singleton address $`x=e_i`$, put
$`s=T\setminus\{i\}`$; then $`\Omega_i=Z_{A,T}Z_{A,s}`$.
The formal expression $`\iota(x)=\sum_i x_ie_i`$ uses the actual p
address bits. It does not supply a zero-filled w-bit key.

## 2. The complete kernel and collision condition

Let $`\pi,\rho`$ be permutations of the w-bit labels and $`g,h`$
Boolean functions. With $`d=2^w`$ and rightmost-first operator order,
define

```math
M_{\pi,g}=H_A D_g(A)Q_\pi(A;B)H_A,
\qquad D_g|a\rangle=(-1)^{g(a)}|a\rangle,
\qquad Q_\pi|a,b\rangle=|a,b+\pi(a)\rangle.
```

The diagonal and query commute, so M is a real Hermitian involution.
Summing the intermediate A label gives the exact full kernel

```math
\langle a',b'|M_{\pi,g}|a,b\rangle
=d^{-1}(-1)^{(a+a')\cdot\pi^{-1}(b+b')
                  +g(\pi^{-1}(b+b'))}.             \qquad(2)
```

Thus the inverse genuinely occurs in an interference kernel. This does
not implement inverse evaluation on an arbitrary input bank. Two such
kernels in the same Fourier frame simply fuse:

```math
M_{\rho,h}M_{\pi,g}=H_A D_{g+h}Q_{\pi+\rho}H_A,    \qquad(3)
```

whether or not $`\pi+\rho`$ is a permutation. Instead consider the
crossed word

```math
D_{\sigma,x}|c,b\rangle=(-1)^{b\cdot\sigma_x(c)}|c,b\rangle,
\qquad D_{\sigma,x}=H_BQ_{\sigma_x}H_B,
\qquad J_x=M_{\rho,h}D_{\sigma,x}M_{\pi,g}.
```

For Fourier-B keys $`v,v'`$, set

```math
\begin{aligned}
\phi_{a,v}(c)&=d^{-1}\sum_u
(-1)^{(a+c)\cdot u+v\cdot\pi(u)+g(u)},\\
\psi_{a',v'}(c)&=d^{-1}\sum_u
(-1)^{(a'+c)\cdot u+v'\cdot\rho(u)+h(u)}.
\end{aligned}
```

Then the complete kernel in these coordinates is

```math
(J_x)_{a',v';a,v}
=\sum_{c:\,\sigma_x(c)=v+v'}
      \psi_{a',v'}(c)\phi_{a,v}(c).                \qquad(4)
```

This is an actual collision sum, without clean projection. A weak
principal source with permutation $`F_x`$ and arbitrary unitary helper
fibers $`S_{x,a}`$ is equivalent to

```math
\phi_{a,v}(c)\mathbf1_{v'=v+\sigma_x(c)}
=(S_{x,a})_{v'v}\,\psi_{F_x(a),v'}(c)
\quad\text{for every }a,v,c,v'.                  \qquad(5)
```

For example, the full arbitrary-key correction
$`F_z(a)=a+f(z)+af(z)`$, retaining every physical key bit, would meet
the sufficient [weak alias interface](OUTPUT_ALIASED_CHARACTER_SOURCE.md#3-the-arbitrary-key-correction-has-the-same-weaker-return-interface)
if (5) held. This principal-permutation promise is stronger than the
single-observable contract.

For that observable, only the forbidden parity blocks of (4) must
vanish: their output top-bit parity disagrees with the desired input
parity in (1). For Hermitian involutions $`P,\Omega`$ and unitary J,

```math
\|J^\dagger PJ-\Omega\|
=2\max\{\|\Pi_{P,+}J\Pi_{\Omega,-}\|,
         \|\Pi_{P,-}J\Pi_{\Omega,+}\|\}.          \qquad(6)
```

Every forbidden block therefore needs norm at most half the requested
observable error. Small individual kernel entries are insufficient.

## 3. A fully paid cubic candidate

Take

```math
\begin{aligned}
q(a)&=k(a)_T+a_T=\gamma_2(a)_T,\\
g(a)&=a_Tq(a),\qquad \sigma_x(a)=\iota(x)k(a),\\
M&=H_AD_gQ_kH_A,\qquad J_x=MD_{\sigma,x}M.        \qquad(7)
\end{aligned}
```

The quadratic q is independent of $`a_T`$, so g is cubic. Use the
complete exact query bounds $`T_k,G_k,d_k,h_k`$ from the
[quadratic-inverse ledger](FULL_BANK_CHARACTER_TRANSPORT.md#5-the-quadratic-inverse-has-a-paid-xor-query-not-a-paid-in-place-map),
where $`d_k\le40`$ and every native helper returns. The required actual
address multiplication query is

```math
Q_\times:|x,E,B\rangle\longmapsto|x,E,B+\iota(x)E\rangle.
```

Its monomials are $`x_iE_S`$ for nonempty S not containing i. The
[returned quadratic-query compiler](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#3-the-full-quadratic-xor-query-is-paid)
uses one arbitrary mask per monomial. With $`R=p(m-1)`$, its common
exact bounds are

```math
T_\times\le8R,\qquad d_\times\le8,\qquad
G_\times\le18R+2w,\qquad h_\times\le R.          \qquad(8)
```

There are sixteen parity CNOTs and two linear-toggle CNOTs per mask,
plus the enclosing B-Hadamard pair. Physical T targets are distinct
private masks; shared address/source/output wires are Clifford
controls. No initialized product temporary is assumed.

Reuse one additional arbitrary w-bit bank E. Chronologically,

```math
Q_k(A;E),\ CZ(A_T,E_T),\ Q_k(A;E),\ CZ(A_T,E_T),\ Z(A_T)
```

implements $`D_g`$: the first four gates give
$`a_Tk(a)_T=a_T+g(a)`$, and the last Z removes $`a_T`$. The cost is
$`(2T_k,2G_k+8,2d_k)`$, charging three elementary Cliffords per CZ
and two per $`Z=S^2`$. Likewise

```math
Q_k(A;E),\ Q_\times(x,E;B),\ Q_k(A;E),\ Q_\times(x,E;B)
```

implements $`Q_{\sigma_x}`$: the dirty E contributions cancel and E
returns. Enclose it by $`H_B`$ to obtain $`D_{\sigma,x}`$. Thus

```math
\begin{aligned}
T_J&\le8T_k+2T_\times,\\
G_J&\le8G_k+2G_\times+6w+16,\\
d_J&\le8d_k+2d_\times\le336.                    \qquad(9)
\end{aligned}
```

Beyond principal A, private allocation is
$`2w+\max(h_k,h_\times)`$: B, E, and a reusable native-helper pool.
Relative to only the m logical output wires, include the other A
coordinates, giving $`5m-3+\max(h_k,h_\times)`$ private qubits.
Clean allocation is zero. These simultaneous bounds fit
$`T,w_{\mathrm{private}}=O(2^p\mathrm{poly}(p))`$ and
$`G=O(4^p\mathrm{poly}(p))`$.

The actual reversed native word has the same bounds. Using actual
adjoints at reverse positions, complete primitive errors
$`\eta_k,\eta_\times`$ give

```math
\|\widetilde J-J\|\le8\eta_k+2\eta_\times,
\qquad
\|\widetilde J^\dagger P\widetilde J-J^\dagger PJ\|
\le16\eta_k+4\eta_\times.                       \qquad(10)
```

These are implementation errors relative to the candidate, not target
errors. The ideal E and native helpers return exactly. B and the private
A coordinates need not return after J. Their return after the probe
echo would require the complete observable identity, which fails next.

## 4. Protected sectors exclude diagonal middles

Allow arbitrary Boolean A-only chirps g,h, with
$`M_g=H_AD_gQ_kH_A`$ and $`P_u^{(g)}=M_gZ_A^uM_g`$. Direct
conjugation gives

```math
P_u^{(g)}=Z_A^uH_AD_{\Delta_ug}Q_{\Delta_uk}H_A,
\qquad \Delta_uk(a)=k(a+u)+k(a).                 \qquad(11)
```

The derivative factors are invariant under translation by u, so this
ordering is exact. With the singleton-address s from Section 1,

```math
\Delta_Tk=e_T,\qquad \Delta_sk=e_s+a_ie_T.        \qquad(12)
```

Consequently $`P_T^{(h)}`$ acts identically on $`B_s`$, whereas
$`P_s^{(g)}`$ contains one unconditional $`X_{B,s}`$ and all its other
factors act trivially on that qubit. The A-only derivative phases do
not affect these support statements. Therefore

```math
[P_T^{(h)},Z_{B,s}]=0,\qquad
\{M_g\Omega_iM_g,Z_{B,s}\}=0.                   \qquad(13)
```

Let C be any unitary satisfying $`[C,Z_{B,s}]=0`$, and set
$`J=M_hCM_g`$. Conjugating actual probe and target by $`M_g`$ gives

```math
M_g(J^\dagger Z_{A,T}J)M_g=C^\dagger P_T^{(h)}C,
\qquad M_g\Omega_iM_g=P_T^{(g)}P_s^{(g)}.        \qquad(14)
```

The first preserves computational $`B_s`$ sectors; the second flips
them. On any unit vector in one such sector their outputs are
orthogonal and both have norm one. Thus the complete observable error
is at least $`\sqrt2`$. C may move the address: use the full controlled
target and choose input address $`|e_i\rangle`$ in this proof.

Every computationally diagonal middle satisfies the premise, even if
its phase depends on $`B_s`$. Changing the support of sigma or using
distinct outer chirps cannot repair this architecture. A simple original
product witness is available because

```math
M_g|+^w\rangle_A|b\rangle_B
=(-1)^{g(0)}|+^w\rangle_A|b\rangle_B.
```

Hence $`|e_i\rangle|+^w\rangle_A|b\rangle_B`$ has squared defect
exactly two. This is a legal input proving failure, not an algorithmic
initialization promise. The statement tensors with every reference and
returned helper identity.

## 5. A plain XOR middle has an exact blind sector

Replacing $`D_\sigma`$ by $`Q_\sigma`$ is paid. For
$`\sigma_i(a)=e_ik(a)`$ it still preserves $`B_s`$: every product
coordinate contains i, whereas s does not. Section 4 already applies.

More generally, every XOR query on B is identity on the whole
Fourier-key-zero sector $`A\otimes|+^w\rangle_B`$. There

```math
M_g=(H_AD_gH_A)\otimes I_B,
\qquad M_gQ_\sigma M_g=I.                       \qquad(15)
```

The transported probe is exactly $`Z_{A,T}`$ on that sector, whereas
the singleton target is $`Z_{A,T}Z_{A,s}`$. A computational A label
with $`a_s=1`$ proves operator-norm error exactly two. This covers even
XOR queries that change $`B_s`$. It is an entire exact helper-sector
identity, not an average or a small kernel entry.

## 6. Unbalanced derivatives exclude translation-plus-phase middles

Fix the singleton address and allow any read-only-A middle

```math
C|a,b\rangle=e^{i\theta(a,b)}|a,b+\tau(a)\rangle,
\qquad CZ_{B,s}C^\dagger=D_{\tau_s}(A)Z_{B,s}.   \qquad(16)
```

The arbitrary real phase theta cancels. Neither theta nor tau needs an
efficiency bound for this screen. Write $`f=\Delta_Th`$ for the
probe-side chirp; it is independent of $`a_T`$. Then

```math
C_f=H_{A,\mathrm{proper}}D_fH_{A,\mathrm{proper}},
\qquad P_T^{(h)}=Z_{A,T}C_fX_{B,T}.
```

If f is unbalanced, $`\mathrm{Tr}(C_f)\ne0`$. One eigenspace of
$`C_f`$, of eigenvalue kappa, has dimension strictly greater than
$`N/2`$, where $`N=2^{w-1}`$. At fixed $`a_T=t`$, one eigenspace
of $`D_{\tau_s}`$, of eigenvalue epsilon, has dimension at least
$`N/2`$. Their intersection contains a unit vector u with

```math
C_fu=\kappa u,\qquad D_{\tau_s}u=\epsilon u.
```

No commutation between these operators is assumed. Put
$`\phi=|t,u\rangle_A|b\rangle_B`$. Both phi and
$`P_T^{(h)}\phi=(-1)^t\kappa|t,u\rangle|b+e_T\rangle`$
have eigenvalue $`\epsilon(-1)^{b_s}`$ under
$`D_{\tau_s}Z_{B,s}`$. Equation (16) therefore places
$`\psi=C^\dagger\phi`$ and
$`C^\dagger P_T^{(h)}C\psi`$ in the same computational
$`B_s`$ sector. The transformed target flips that sector by (13).
Its unit output is orthogonal, proving

```math
\|(M_hCM_g)^\dagger Z_{A,T}(M_hCM_g)-\Omega_i\|
\ge\sqrt2.                                      \qquad(17)
```

Equivalently the relevant anticommutator has norm two. The common-vector
proof supplies a full norm-one sector-preserving output; the conclusion
does not follow merely from a weaker bound on a diagonal block.

For the concrete q in (7), the proper coordinates form disjoint
complement pairs. Each pair contributes two to its signed sum, so
$`\mathrm{Tr}(C_q)=2^{(w-1)/2}>0`$. Thus (17) excludes every
translation-plus-phase middle for this cubic chirp, including any word
of such operations. Distinct outer chirps are allowed: only the
probe-side top derivative must be unbalanced. Balancedness is a necessary
escape from this screen, not a sufficient source construction.

## 7. A balanced repair changes the invariant, not the frontier

Adding linear terms to q cannot balance it: its proper-coordinate polar
form is nondegenerate. Removing a complementary pair and adding a linear
term on the new radical does balance it:

```math
h'(a)=a_T(q(a)+a_ua_v+a_u),\qquad v=T\setminus u.
```

This costs one extra literal CCZ and one CZ beyond the paid chirp,
preserving constant native depth and all polynomial resource envelopes.
Its Fourier conjugate contains $`X_{A,u}`$; the trace obstruction
really disappears. Sections 4 and 5 still apply for every chirp.

The natural choice $`u=\{i\}`$, $`v=s`$ has another protected sector.
Write $`q_{\mathrm{rest}}=q+a_ia_s`$, independent of both
$`a_i,a_s`$. For the symmetric kernel using this chirp,

```math
P_T'=Z_{A,T}X_{A,i}C_{\mathrm{rest}}X_{B,T},
\qquad P_s'=Z_{A,s}X_{B,s}V,
\qquad V=H_{A,i}\,\mathrm{CNOT}_{A_i\to B_T}\,H_{A,i}.
```

Whenever $`[C,Z_{A,i}]=0`$, the actual $`C^\dagger P_T'C`$ flips
$`A_i`$. On the legal $`B_T`$ Fourier-key-one sector, V restricts
to $`X_{A,i}`$ and the transformed target is

```math
P_T'P_s'=-Z_{A,T}Z_{A,s}C_{\mathrm{rest}}X_{B,s},
```

which preserves $`A_i`$. A fixed computational $`A_i`$ input with
$`B_T=|-\rangle`$ therefore gives another $`\sqrt2`$ gap. C may
have arbitrary bank dynamics and need not preserve that Fourier key.
This particular balanced repair relocates the obstruction to A; it is
not a statement about all balanced derivatives.

The candidate (7) has no unpaid primitive: it is a fully paid failed
word. Its group property $`J_xJ_y=J_{x+y}`$ is not itself an
obstruction; the known complete source $`FX^{\iota(x)}F^\dagger`$
has the same property. The first unpaid operation for a replacement is
still complete forbidden-parity cancellation in (4), or another actual
full-space identity with the qualifying native/error/return ledger.
An inverse appearing in (2) does not supply that operation.

Changed outer derivative-support frames, different balanced derivative
geometry with sufficient A motion, and bank operations beyond controlled
translations and phases remain outside the corresponding hypotheses.
These include nonlinear B permutations and genuinely nonmonomial
dynamics. No general impossibility or faster complete source follows.

<a id="normalized-forward-gram-source"></a>

## 8. A flat inverse kernel has a normalized forward-coordinate source

This trial starts with genuine Fourier interference and pays the
normalization. Its complete Gram operator nevertheless selects the
forward coordinate, not the inverse coordinate visible in a scalar
kernel. The distinction is useful when proposing a replacement.

Let $`\pi_x`$ be a permutation of n-bit strings with a supplied complete
zero-clean XOR query $`Q_\pi:(b,a)\mapsto(b,a+\pi_x(b))`$.
The two principal banks A and B are arbitrary, apart from one signal
coordinate of A specified below. The address may be coherent. For the
nilpotent pencil take $`\pi_x(b)=A_xb`$, $`A_x=N_x^{\mathsf T}`$;
its bilinear query and full native ledger are already paid in
[Section 3 of the pencil proof](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#3-the-full-quadratic-xor-query-is-paid).

Put $`d=2^n`$, let F be Hadamards on both banks, and define

```math
f_x(a,b)=a\cdot\pi_x(b),\qquad
M_x=F D_{f_x}F=H_BQ_\pi H_B.
```

With the unitary Fourier normalization $`1/d=2^{-n}`$ on 2n bits,
the Walsh sum is exactly

```math
\frac1d\sum_{a,b}(-1)^{a\cdot\pi_x(b)+u\cdot a+v\cdot b}
 =(-1)^{v\cdot\pi_x^{-1}(u)}.                         \qquad(18)
```

Consequently M is a Hermitian involution, and its full kernel is

```math
\langle a',b'|M_x|a,b\rangle
 =d^{-1}(-1)^{(b+b')\cdot\pi_x^{-1}(a+a')}.          \qquad(19)
```

This is the standard Maiorana--McFarland duality; see
[Rötteler, *Quantum algorithms for highly non-linear Boolean functions*,
Lemma 4](https://arxiv.org/abs/0811.3208).
The dual phase oracle in that paper is a separate oracle assumption.
Equation (19) supplies its phase as a matrix entry, not as a complete
diagonal gate or a supplied clean state.

Choose $`s=A_j`$ as the sole clean signal, with all other A inputs
arbitrary. Let $`R=T(Y_s)=e^{i\pi/8}R_y(\pi/4)`$, where
$`R_y(\theta)=e^{-i\theta Y_s/2}`$, and define the actual word

```math
U=R^3 M R^\dagger M^\dagger.                         \qquad(20)
```

Its rightmost occurrence is the actual adjoint. In an ideal Fourier-B
sector b, M translates A by $`\pi_x(b)`$. All translations outside s
commute through the signal rotations and cancel. With
$`q=\pi_x(b)_j`$ and $`\omega=e^{i\pi/4}`$, the complete action is

```math
R^3X_s^qR^\dagger X_s^q
 =\begin{cases}
 \omega R_y(\pi/2),&q=0,\\
 \omega R_y(\pi),&q=1.
 \end{cases}                                          \qquad(21)
```

Both ideal cases hold on either signal input. Thus every other target-bank
wire and every query helper returns exactly for the ideal query, on arbitrary inputs.
Define

```math
C_{x,j}=H_B D_{\pi_x(\cdot)_j}H_B,\qquad
\Pi_{x,j}^{\rm fwd}=\frac{I+C_{x,j}}2.
```

For the clean embedding J, (21) gives the genuinely normalized block
and its Gram operator

```math
K=J^\dagger UJ=\frac{\omega}{\sqrt2}\Pi_{x,j}^{\rm fwd},\qquad
K^\dagger K=\frac12\Pi_{x,j}^{\rm fwd}.              \qquad(22)
```

No dirty input is averaged or initialized in this calculation. It is
a complete operator identity, so it retains coherent addresses,
private work and arbitrary references. This computes a normalized projector source; its projector has not
thereby been identified with the required lookup sectors.

### 8.1 The pencil trial has the wrong column, with a full gap

For the paid pencil query, choose j=0. Then
$`C_{x,0}=X_B(N_xa_0)`$, whereas the addressed observable has label
$`\nu_x=N_x^{-1}a_0`$. To compare directly with the canonical paired
Gram contract, apply a B-Hadamard envelope and identify B with carrier C.
Conjugate also by the n CNOTs from carrier D to C. The resulting source
has exactly

```math
G_{\rm trial}(x)=\frac14\bigl(I+
 Z_C(N_xa_0)Z_D(N_xa_0)\bigr),\qquad
G_*(x)=\frac14\bigl(I+
 Z_C(N_x^{-1}a_0)Z_D(N_x^{-1}a_0)\bigr).            \qquad(23)
```

The Hadamard walls cancel in the complete source word. It is a
literal two-query signal echo, with the two carrier-pair CNOT walls.
In particular (23) is not a computational estimate of (19).

At x=0 the labels agree. For every $`p\ge2`$, at field address x=1,
their standard-coordinate preimages are respectively

```math
E N_1a_0=e_0+e_1,\qquad
E N_1^{-1}a_0=e_0+e_1+\cdots+e_{m-1},\qquad m=2^p.
```

These labels differ. Two distinct nonzero binary characters disagree
on a legitimate computational carrier assignment. Hence

```math
\|G_{\rm trial}-G_*\|_{\rm op}=\frac12
\quad(p\ge2).                                        \qquad(24)
```

The upper bound one half holds because both expressions in (23) are
one half of projectors. A private carrier coordinate may witness the
failure; those coordinates are arbitrary inputs of the contract.
Restricting them to zero would not repair the complete source. The
small p=1 pencil is involutory and does not exhibit this mismatch.

### 8.2 Native resources, precision and actual return

Let $`T_Q,G_Q,d_Q,h_Q`$ be common bounds for the supplied query and
its actual inverse, including its complete error eta. The paired source
in (23) uses one clean signal, an extra n-minus-one dirty target wires,
and the query pool. Literal $`R^3`$ and $`R^\dagger`$ need two T gates
in two layers; conjugating T by $`SH`$ implements its Y axis. Expanding
all inverse S gates gives the safe common ledger

```math
\begin{aligned}
T_U&\le2T_Q+2,&d_U&\le2d_Q+2,\\
G_U&\le2G_Q+2n+15,&b_U&\le n-1+h_Q,\\
\|\widetilde U-U\|&\le2\eta,&
\|\widetilde K^\dagger\widetilde K-G_{\rm trial}\|&\le4\eta.
\end{aligned}                                         \qquad(25)
```

Here the two n-bit carriers C,D and the address are principal inputs;
$`b_U`$ counts only additional dirty work. The source's actual reverse
has the same common bounds. Equation (21), not a statewise measurement
claim, pays all target and helper returns; the norm estimate includes
leakage of the one signal. For the exact pencil query with n=pm,
$`r=(m-1)p^2`$, $`s\le(m-1)p^3`$ and $`c_n=n^2+3n`$, this is

```math
T_U\le16r+2,\quad d_U\le18,\quad
G_U\le8c_n+8n+32r+4s+15,\quad b_U\le n-1+r.
                                                               \qquad(26)
```

All simultaneous source budgets fit the relaxed polynomial-overhead
fusion envelope. Its target Gram error is nevertheless one half before
any compilation error. Applying the paid readout and fusion to the wrong
projector cannot move the complete-indicator or frame frontier.

<a id="clifford-signal-gram-gap"></a>

## 9. A Clifford signal image cannot supply the lookup Gram

There is a broader screen for the *single* cubic Fourier trial, including
Clifford boundary changes and genuine address mixing. It concerns the
actual clean signal image and does not exclude all higher-band sources.
Let U have one clean signal and arbitrary other inputs, and suppose

```math
A=U^\dagger Z_sU
```

is a Clifford. Since A is also a Hermitian involution, put
$`H=J^\dagger AJ=2K^\dagger K-I`$. If $`Q=(I-\Omega_C\Omega_D)/2`$,
then Gram error delta from the required $`(I-Q)/2`$ is exactly
$`\|H+Q\|/2`$.

A one-clean compression of any Clifford has

```math
H^\dagger H\in\{0,I,I/2\}\quad\hbox{or}\quad
H^\dagger H=(I\pm P)/2                               \qquad(27)
```

for a nonidentity Hermitian Pauli P on the complete non-signal space.
Indeed expand $`JJ^\dagger=(I+Z_s)/2`$ between A and its adjoint.
The conjugate $`A^\dagger Z_sA`$ is a Pauli; its clean compression is
zero if its signal factor is X or Y, and otherwise is a signed Pauli.
No bound on the dirty width or on native address motion is used.

Q is a nontrivial half-rank projector. The first three alternatives in
(27) give $`\|H+Q\|\ge1/\sqrt2`$: use a vector in its range if H=0,
and a vector in its kernel if H is unitary or has all singular values
$`1/\sqrt2`$. In the remaining case put $`P_0=H^\dagger H`$.
Since H is a partial isometry with initial projector $`P_0`$,

```math
\|H+Q\|\ge\max\{\|P_0(I-Q)\|,\|Q(I-P_0)\|\}
             =\|P_0-Q\|.                              \qquad(28)
```

The two inequalities follow by taking inputs in the kernels of Q and
$`P_0`$, respectively. They do not presume that H and Q commute.

For p at least one, Q is diagonal but is not a Pauli half-projector:
its sign is the nonaffine function $`Y_x\mathbin\oplus Y'_x`$.
Every Pauli half-projector has distance at least $`1/\sqrt2`$ from Q.
For a diagonal Pauli the differing sign at a computational basis state
gives distance one. For a Pauli with nonzero X label v, decompose the
computational basis into pairs $`\{z,z+v\}`$. The Pauli projector
has rank one on each pair. If Q has two equal diagonal entries there,
the projector distance on that pair is one; if its entries differ, the
distance is $`1/\sqrt2`$. Literal Pauli phases only change the
unit-modulus off-diagonal entry and do not change either value.
Thus

```math
\boxed{\|K^\dagger K-\Pi_{\rm eq}/2\|\ge1/(2\sqrt2),
       \qquad p\ge1,\quad U^\dagger Z_sU\text{ Clifford}.}
                                                               \qquad(29)
```

A Clifford sandwich of a Boolean cubic phase satisfies this hypothesis:
its derivative in any Pauli translation is quadratic, hence its
conjugation of any Pauli is Clifford. Therefore a single cubic
Maiorana--McFarland Fourier phase, even with arbitrary Clifford input
and output changes, cannot supply the target Gram at small error.
The same proof quantifies the escape. It used only the Clifford
compression identity (27), which holds even for a non-Hermitian Clifford C.
Thus any actual source at Gram error delta satisfies

```math
\inf_{C\ {\rm Clifford}}\|U^\dagger Z_sU-C\|
 \ge1/\sqrt2-2\delta.                               \qquad(29a)
```

The clean signal image must stay a constant distance from every Clifford
at small target error. This is a statement about the actual source, not
its gate-count label or an arbitrary product of third-level gates.
The doubled source (20)
need not have a Clifford signal image and is instead rejected by its
explicit forward-sector identity (23).

## 10. What the normalized trial changes and what remains open

The trial resolves a concrete normalization question positively:
a complete, one-clean, constant-depth source of the form (22) is paid,
with every arbitrary target-bank fiber canceled. It also resolves its
lookup-sector question negatively: the inverse in the flat kernel is
not the inverse in the accepted Gram, and the pencil specialization
has the exact full-space gap (24).

A next interference word must change the accepted sector itself, for
example by a second genuinely crossed kernel or a non-Clifford signal
image with address motion. An additional Fourier wall, a relabeling of
the scalar dual, or recovery applied to (23) does not do this. The first
unpaid identity is a normalized accepted Gram equal or uniformly close
to the **inverse-column** paired lookup projector on arbitrary carriers.
No claimed source here meets that identity, and no unrestricted native
T-depth bound follows from (29).

The [exact fixtures](../../tests/test_bent_gram_source.py) check a cubic
Maiorana--McFarland phase with a higher-degree inverse dual, every
column and actual inverse of the normalized signal echo, the GF(4)
forward/inverse column mismatch, and the one-clean Clifford-compression
alternatives. They are finite guards on the analytic identities, not
an asymptotic source emitter or external mathematical review.

<a id="normalized-half-bilinear-source"></a>

## 11. A normalized half-bilinear commutator with dense bank mixing

This second discovery trial leaves the single cubic-sandwich grammar.
It uses a wide non-Clifford phase on both arbitrary banks, followed by
fixed dense Clifford mixing. Normalization is exact; finite fixtures
reject the proposed lookup sectors. No all-size impossibility is inferred
from those fixtures, and this is not a qualifying source.

### 11.1 The complete native family

Use the pencil's standard dual coordinates on two n-bit carriers z,t,
where n=pm. The fixed changes from physical carriers are
$`z=E^{-\mathsf T}C`$, $`t=E^{-\mathsf T}D`$; every carrier bit remains
arbitrary. Write

```math
 B_x=\bar N_x^{\mathsf T}-I,\qquad
 q_x(z,t)=\sum_{i,j}(B_x)_{ij}z_it_j\pmod4,\qquad
 D_x|z,t\rangle=i^{q_x(z,t)}|z,t\rangle.                 \qquad(30)
```

The coefficient $`(B_x)_{ij}`$ is its binary value, even when expressed
as an XOR of address bits. Its integer lift modulo four must not be
replaced by the integer sum of those bits. B_x has at most
$`r=(m-1)p^2`$ nonzero affine coefficient forms, each supported on at
most p address bits. The fixed binary change E is the Vandermonde
basis from the pencil proof, not a data-loading operation.

Each factor in (30) is the three-control phase $`i^{abc}`$, with
$`a=(B_x)_{ij}`$, $`b=z_i`$, $`c=t_j`$. It has an exact zero-clean
word whose physical T targets are all private. Borrow dirty bits d,e,f.
Let $`T_{ab;d}`$ be the returned eight-T Toffoli with phase helper f,
and let L(d,c) be the chronological word

```math
 T_d,\quad\mathrm{CX}(c,d),\quad T_d^\dagger,
 \quad\mathrm{CX}(c,d).
```

Its phase is $`\omega^{d-(d\oplus c)}`$, with
$`\omega=e^{i\pi/4}`$. The chronological echo
$`T_{ab;d},L,T_{ab;d},L^\dagger`$ gives
$`i^{abc}(-1)^{abcd}`$ for every old d. Cancel the residual sign by

```math
 T_{ab;e},\quad\mathrm{CCZ}(e,c,d),\quad
 T_{ab;e},\quad\mathrm{CCZ}(e,c,d).
```

Each CCZ is H on private e, the returned Toffoli with controls c,d,
and H on e. All d,e,f return literally. There are six returned
Toffolis and four private T/T-dagger gates: at most 52 T gates and
52 T layers. All original a,b,c enter only as Clifford controls.
For an affine form a, replace every control-to-private CNOT by its
literal parity CNOTs; no address wire is used as a parity target or
as a shared physical T target. A safe common Clifford bound is
$`24p+80`$ per factor, including all parity incidences. Separate
three-bit private pools let all r factors align in the same 52
native T layers, even with overlapping address and carrier controls.
Thus D and its actual inverse obey

```math
 T_D\le52r,\quad G_D\le(24p+80)r,\quad d_D\le52,
 \quad h_D\le3r.                                      \qquad(31)
```

The [exact half-bilinear fixtures](../../tests/test_half_bilinear_gram_source.py)
check the 52-T word on every input of all six wires, including dirty
values and its actual reverse. Equation (31) uses the analytic
read-only-affine-control schedule; it is not inferred from a single
numerical gate test.

For a specified dense mixer, define the invertible binary map on both
banks and its native Clifford

```math
 M(z,t)=(Ez+t,z),\qquad R=H^{\otimes2n}U_M,
 \qquad E_x=R^\dagger D_xR,
 \qquad V_x=D_x^\dagger E_x^\dagger D_xE_x.             \qquad(32)
```

R mixes both banks, and its full count is paid. No intermediate
carrier is initialized. With $`c_n=n^2+3n`$, Gaussian elimination
and literal CNOT row swaps give $`G_R\le4n^2+8n`$, also for the
actual inverse. The coherent direct sum V over x fixes the address;
it can still mix all carrier states nonmonomially. At x=0, D and V
are identity.

Let $`\Pi_0=(I+Z_{z_0}Z_{t_0})/2`$. The fixed two-T, 19-Clifford
one-signal source $`U_0`$ from the
[four-call fixture](../../tests/test_gram_four_call_readout.py) acts as
$`H_a`$ on equal z_0,t_0 and $`X_a`$ otherwise. Its accepted block is
$`\Pi_0/\sqrt2`$. Explicitly, let B be the chronological native word
$`S_a,S_a,S_a,H_a,T_a,H_a,S_a`$, whose operator is
$`B=S_aH_aT_aH_aS_a^\dagger`$. Literal conjugation gives
$`BZ_aB^\dagger=H_a`$. The controlled-H with control z_0 is the
chronological word $`B^\dagger,H_a,\mathrm{CX}(z_0,a),H_a,B`$,
where B-dagger is the actual reversed word. Emit U_0 chronologically as

```math
 \mathrm{CX}(t_0,z_0),\quad H_a,\quad
 \mathrm{CH}(z_0;a),\quad\mathrm{CX}(z_0,a),\quad
 \mathrm{CX}(t_0,z_0).
```

On even z_0,t_0 parity only the first signal H acts. On odd parity,
the two successive H gates cancel and the signal X remains. Both
carrier bits return, with no input-dependent scalar phase. The two
B orientations each cost one T gate and six Cliffords in the positive-S
alphabet, so the full word and its actual inverse have two T gates,
19 Cliffords and at most two T layers. Define the actual source, with
the physical-to-dual
change $`C_E=U_{E^{-\mathsf T}}\otimes U_{E^{-\mathsf T}}`$, by

```math
 U=C_E^\dagger V^\dagger U_0VC_E,\qquad
 K=\frac1{\sqrt2}C_E^\dagger V^\dagger\Pi_0VC_E,
 \qquad G=K^\dagger K=\frac12C_E^\dagger V^\dagger\Pi_0VC_E.
                                                               \qquad(33)
```

This is a parameterized native family, with a complete accepted
projector and its normalization derived, not just a proposed kernel.
The one source signal is the only clean input. Every other work bit
is arbitrary and the private phase pools return exactly. The common
forward/actual-inverse ledger is

```math
\begin{aligned}
 T_U&\le416r+2,&d_U&\le418,\qquad b_U\le3r,\\
 G_U&\le8(24p+80)r+8(4n^2+8n)+4c_n+19.
\end{aligned}                                          \qquad(34)
```

The address and two n-bit carriers are excluded from b_U. Their extra
physical allocation is charged by the readout ledger. Equation (34)
fits the relaxed polynomial count/width and quadratic-times-polynomial
Clifford envelopes. It is exact, so there is no hidden precision state,
per-row error sum, or supplied controlled fine-angle gate. The missing
contract is the lookup sector, not normalization or a native schedule.

### 11.2 Exact exponential sum and finite sector failures

For clarity, the whole Gram can be evaluated without interpreting any
trace as a returned scalar. Put $`Q=2^{2n}`$ and let r,s,t,w,w' here
denote binary 2n-bit carrier labels. The literal matrix entries of (32)
are

```math
\begin{aligned}
 R_{r,s}&=Q^{-1/2}(-1)^{r^{\mathsf T}Ms},\\
 (E_x)_{r,s}&=Q^{-1}\sum_w i^{q_x(w)}
                       (-1)^{w^{\mathsf T}M(r+s)},\\
 (V_x)_{r,s}&=Q^{-2}\sum_{t,w,w'}
  i^{-q_x(r)-q_x(w)+q_x(t)+q_x(w')}
  (-1)^{w^{\mathsf T}M(t+r)+{w'}^{\mathsf T}M(t+s)}.
\end{aligned}                                          \qquad(35)
```

The accepted Gram in dual coordinates is therefore

```math
G_{r,s}=\tfrac12\sum_{u:u_{z_0}=u_{t_0}}
\overline{(V_x)_{u,r}}(V_x)_{u,s}.
```

These finite sums are exact operator formulas, not an efficient classical
evaluation claim.

The required remaining identity is

```math
 V_x^\dagger Z_{z_0}Z_{t_0}V_x
   = Z_z(v_x)Z_t(v_x),\qquad
 v_x=(1,x,\ldots,x^{m-1})^{\mathsf T},                  \qquad(36)
```

with field entries expanded in the fixed binary basis. Approximate
versions need a uniform full-operator error, not only diagonal entries.
For the explicit mixer (32), the exact p=1, x=1 fixture fails (36):
on $`z=e_1,t=0`$ the transported probe is +1 while the desired probe
is -1. The complete Gram error is consequently the maximal one half.
Changing the two commutator exponents independently from +1 to -1,
or using the tested pairs (2,1), (1,2), (2,2), does not fix this input.

The p=2, x=1 fixture tests the same second-field-coordinate input
$`z=e_{p},t=0`$, using exact integer Walsh transforms on all 65,536
carrier amplitudes. The desired probe is -1. The default (+1,+1)
commutator has exact probe expectation $`19/64`$, giving Gram error
at least $`83/256`$. The sign variants (+1,-1), (-1,+1), (-1,-1)
have expectations $`39/128,39/128,19/64`$, respectively. The three
tested exponent-two pairs give $`3/32,1/4,0`$. Thus neither a literal
sign correction nor these bounded angle changes repair the tested
sectors. This is exact finite evidence, not a general lower bound or
an exclusion of a family using small-p fallbacks.

The trial is paused at (36). Reopening requires an actual changed
mixer/interference word whose complete conjugated probe meets that
identity or a dimension-uniform approximate version. An address-sensitive
nonmonomial mixer is one concrete change outside this fixed R word;
it must be specified and paid, not inferred from the presence of dense
mixing. The already paid four-call recovery cannot correct a wrong
projector. Neither this finite failure nor the Clifford-image screen
excludes general higher-band sources or moves the frame frontier.

### 11.3 A paid address-controlled shear does not repair the sectors

The next literal trial replaces the fixed mixer by

```math
 S_x=Q_{\bar N_x^{\mathsf T}}(z;t),\qquad
 S_x(z,t)=(z,t+\bar N_x^{\mathsf T}z),\qquad
 R_x=RS_x=H^{\otimes2n}U_{M S_x}.                     \qquad(37)
```

Keep the definitions of E_x, V_x and the source (33), with R_x in
place of R. This is address-controlled bank motion; it leaves the
address register unchanged. The accepted Gram remains exactly
$`C_E^\dagger V^\dagger\Pi_0VC_E/2`$ on every carrier, helper and
reference input. Its full matrix formula is (35) with M replaced by
M S_x; no dirty input is fixed to zero.

The standard-coordinate query S_x is the returned query in
[Section 3 of the pencil proof](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#3-the-full-quadratic-xor-query-is-paid)
with its four fixed basis changes omitted. Its actual reverse has the
same cost. In chronological order, reading the first row and then the
second, the emitted V_x is

```math
 \begin{gathered}
 S_x,\ R,\ D_x,\ R^\dagger,\ S_x^\dagger,\ D_x,\\
 S_x,\ R,\ D_x^\dagger,\ R^\dagger,\ S_x^\dagger,\ D_x^\dagger.
 \end{gathered}
```

Here each dagger denotes the actual reversed native word, even though
the ideal query S_x is an involution. The complete source uses eight
query calls in addition to (34). Their private pools can reuse the
existing three-r-bit pool between sequential calls. Thus the common
source/actual-inverse ledger is

```math
\begin{aligned}
 T_U&\le480r+2,&d_U&\le482,&b_U&\le3r,\\
 G_U&\le8(24p+80)r+8(4n^2+8n)+4c_n+19\\
    &\quad+8(3n+16r+2s).
\end{aligned}                                           \qquad(38)
```

The source is exact and still has one clean signal. This pays the
changed mechanism before testing its sectors; its constant depth and
polynomial overhead fit the same conditional fusion envelope.

At p=1, x=1, the two binary matrices E and bar-N-transpose coincide:
$`E=\bar N_1^{\mathsf T}=\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)`$.
Consequently $`M S_1(z,t)=(t,z)`$. The mixer is a whole-bank swap
followed by four Hadamards. D_1 is the phase $`i^{z_0t_1}`$ on
the pair z_0,t_1; $`R_1^\dagger D_1R_1`$ is its X-basis copy on
the disjoint pair t_0,z_1. They commute. Therefore V_1=I for every
integer choice of the two commutator exponents. Also V_0=I because
D_0=I. The complete p=1 V, including coherent addresses, is identity.
On x=1,z=e_1,t=0 the required probe in (36) is minus one and the
actual probe is plus one, so the complete Gram error is one half.
The exact fixture checks every column of this cancellation using the
literal binary mixer, not only that witness.

At p=2, x=1,z=e_2,t=0, the exact integer Walsh fixture gives probe
expectation $`147/512`$ for the default commutator, with norm exactly
one. Hence its Gram error is at least $`659/2048`$. The sign variants
(+1,-1), (-1,+1), (-1,-1) give $`129/512,129/512,147/512`$;
the exponent pairs (2,1), (1,2), (2,2) give $`3/128,0,0`$.
These are full-state exact receipts at a specific input, not an
all-size exclusion or an obstacle to a separately paid small-p fallback.

There is a stronger obstruction to repairing this p=1 commutator by
changing its mixer alone. For the whole p=1 circuit, including the
address control and arbitrary dirty work, D is one projector
quarter-turn, $`D=I+(i-1)P`$. Let E be any unitary and put
$`W=D^\dagger E^\dagger D E`$. If a computational input s has
$`D|s\rangle=d_s|s\rangle`$, then

```math
 \begin{aligned}
 \langle s|W|s\rangle
   &=\overline{d_s}\langle Es|D|Es\rangle
     =\overline{d_s}\bigl(1-\lambda+i\lambda\bigr),\\
 \lambda&=\langle Es|P|Es\rangle\in[0,1],\\
 |\langle s|W|s\rangle|^2
   &=(1-\lambda)^2+\lambda^2\ge1/2.
 \end{aligned}                                           \qquad(39)
```

This uses only the two eigenvalues 1 and i of D; it does not require
E to be a conjugate of D. In particular it applies to every unitary
mixer R in the present commutator, even one that mixes address,
carriers and additional dirty work.

Choose s with x=1,z=e_1,t=0 and arbitrary fixed work bits. Its initial
probe Z_z0 Z_t0 has eigenvalue plus one and its desired probe has
eigenvalue minus one. Equation (39) puts at least one half of the
final weight in the initial probe's plus sector. The final probe
expectation is therefore nonnegative, and consequently

```math
 \|G-G_*\|\ge1/4\qquad
 \text{for this p=1 quarter-turn commutator.}           \qquad(40)
```

Reversing D replaces i by minus i and leaves (39) unchanged; the
unitary E may be arbitrary, including either sign or any tested
integer power of its phase. This does not apply when the outer D
has been squared, to general higher-band sources or arbitrary longer
interference words, or to p at least two, where D_x need not have
only the two eigenvalues 1 and i. It does not exclude an asymptotic
family with a separate small-p construction.
The specific shear (37) is paused at its failed sector identity.
Reopening requires a changed interference word with correct complete
lookup sectors, rather than another normalization or readout step.

### 11.4 Genuine address Hadamards expose an explicit return failure

Two further complete words mix the address coherently. Let
$`F_X=H^{\otimes p}`$ act on the address, retain the controlled mixer
R_x from (37), and choose either
$`\mathcal R_L=F_X(\bigoplus_xR_x)`$ or
$`\mathcal R_R=(\bigoplus_xR_x)F_X`$. Put
$`D=\bigoplus_xD_x`$, $`E=\mathcal R^\dagger D\mathcal R`$ and
$`V=D^\dagger E^\dagger D E`$, then use the same actual source (33).
Normalization remains exact on the full address-and-carrier space.
The source and its actual reverse add exactly eight address-Hadamard
walls to (38), hence at most 8p Cliffords and no T gates, clean inputs,
dirty width or precision cost.

The complete kernels distinguish these from controlled-bank motion.
Write L_x=M S_x for the binary map in (37), Q=2^(2n), and use x,y
for address labels and r,s,u for carrier labels. Then

```math
 \begin{aligned}
 (\mathcal R_L)_{(y,u),(x,s)}
   &=(mQ)^{-1/2}(-1)^{x\cdot y+u^{\mathsf T}L_xs},\\
 (\mathcal R_R)_{(y,u),(x,s)}
   &=(mQ)^{-1/2}(-1)^{x\cdot y+u^{\mathsf T}L_ys},\\
 (E_L)_{(x,r),(x',s)}
   &=(mQ)^{-1}\sum_{y,u}i^{q_y(u)}
       (-1)^{y\cdot(x+x')+u^{\mathsf T}(L_xr+L_{x'}s)},\\
 (E_R)_{(x,r),(x',s)}
   &=(mQ)^{-1}\sum_{y,u}i^{q_y(u)}
       (-1)^{y\cdot(x+x')+u^{\mathsf T}L_y(r+s)}.
 \end{aligned}                                           \qquad(41)
```

Thus E_R is a convolution on the joint binary address-and-carrier
group, while E_L need not be. Neither fact identifies the required
accepted sectors or returns the address automatically.

The exact p=1 fixture checks all 32 columns. In dual coordinates and
the binary labeling $`j=z+4t+16x`$, put
$`A=V^\dagger Z_{z_0}Z_{t_0}V`$ and $`G_{\rm dual}=(I+A)/4`$.
Both choices give

```math
 A_{9,29}=-(1+i)/4,\qquad A_{18,18}=1,\qquad
 (G_{\rm dual})_{9,29}=-(1+i)/16.                     \qquad(42)
```

The off-diagonal Gram entry connects different addresses and is zero
in the required Gram. The fixed carrier basis changes do not remove
this address-sector failure. It does not imply that intermediate
address motion is forbidden.
The input 18 is x=1,z=e_1,t=0, whose desired probe is minus one;
both words consequently have the maximal Gram gap one half.

The p=2 fixture propagates all 262,144 address-and-carrier amplitudes
with integer Walsh transforms from x=1,z=e_2,t=0. Both resulting
states have norm exactly one. The left choice has probe expectation
$`20011/65536`$ and Gram gap at least $`85547/262144`$; the right
choice has expectation $`6913/16384`$ and gap at least
$`23297/65536`$. Integer amplitudes are converted to unbounded
integers before squaring, so these receipts do not use rounded norm
or phase estimates.

These two words are paused at their complete sector and address-return
failure. The exact kernels and native cost make the trials reviewable;
the finite failures do not exclude other address-moving mechanisms,
longer interference, or a family with separately paid small cases.
The first missing operation remains a complete inverse-column lookup
projector, and no complete-indicator or frame-depth frontier moves here.

<a id="normalized-xor-carry-source"></a>

## 12. An exact T-wall shear supplies a different normalized phase source

This trial replaces the expanded half-bilinear phase by a literal T
wall around the returned quadratic XOR query. It retains the XOR
carries modulo eight automatically, uses fewer private wires, and has
an exact complete accepted Gram. The tested fixed mixer and its address
Hadamard variant still select the wrong sectors. This is a paid partial
construction, not a qualifying lookup source.

### 12.1 The actual phase and native schedule

Keep the standard dual coordinates, n=pm, and
$`B_x=\bar N_x^{\mathsf T}-I`$ from Section 11. Let
$`Q_{B_x}|z,t\rangle=|z,t+B_xz\rangle`$, and put
$`P_t=\bigotimes_{j=0}^{n-1}T_{t_j}`$ and
$`\omega=e^{i\pi/4}`$. The chronological word
$`P_t^\dagger,Q_{B_x},P_t,Q_{B_x}^\dagger`$ gives

```math
 \begin{aligned}
 D_x&=Q_{B_x}^\dagger P_tQ_{B_x}P_t^\dagger,\\
 D_x|z,t\rangle&=\omega^{h_x(z,t)}|z,t\rangle,\\
 h_x(z,t)&=|t\mathbin\oplus B_xz|-|t|\pmod8.
 \end{aligned}                                             \qquad(43)
```

Here the matrix product and XOR are binary, and the bars mean integer
Hamming weight. This is the complete operator, including identity on
the returned query helpers. In particular D_0=I. If
$`a_j=(B_xz)_j`$ is its binary value, then
$`h_x=\sum_j(a_j-2t_ja_j)`$ modulo eight. Replacing a_j by an integer
sum of its contributing monomials would change the operator. For
example, two active XOR summands give a_j=0, not 2. Equation (43)
requires no synthesis of a fine controlled angle or separate carry
register.

The standard-coordinate nonlinear core of the pencil query has

```math
 T_B\le8r,\quad d_B\le8,\quad h_B\le r,\quad
 G_B\le2n+16r+2s,
 \qquad r=(m-1)p^2,\quad s\le(m-1)p^3.                  \qquad(44)
```

This is the already paid eight-layer R91 derivative schedule, with
neither the physical-coordinate E changes nor the identity shear.
Both actual query orientations return all helpers. Its physical T
targets are private. The two intervening T walls target the distinct
principal t wires and occur after query return, so they each need one
native T layer; no simultaneous reuse of a T target is inferred.
Thus D and its actual reverse have

```math
 T_D\le16r+2n,\qquad d_D\le18,\qquad
 G_D\le2G_B,\qquad h_D\le r.                            \qquad(45)
```

All these identities are exact on coherent address, both arbitrary
banks, every arbitrary query helper and references. The transposed
query B_x^T is also a returned bilinear query with the same safe bounds;
its separately tested phase uses $`|t\oplus B_x^{\mathsf T}z|-|t|`$.

This is also a paid relative phase insertion for the complete kernel
translation bank. Let P delete only the first **binary** coordinate
a_0, and put $`K_x=(I+B_x)P`$. The binary transpose convention gives
$`v_x^{\mathsf T}(I+B_x)=e_0^{\mathsf T}`$, so
$`\mathrm{im}K_x=\ker(v_x^{\mathsf T})`$, of codimension one.
The other p-1 binary directions in the first field block are retained;
they are fixed columns of K_x. Deleting that entire field block
instead would give the smaller codimension-p, field-valued kernel.
Since B_x annihilates the first field block, $`B_xP=B_x`$. Write
$`D(A)=Q_A^\dagger P_tQ_AP_t^\dagger`$ for any binary shear A.
The exact diagonal phases then give

```math
 D(K_x)D(P)^\dagger
   =Q_P D(B_x)Q_P,
 \qquad
 |t+K_xz|-|t+Pz|
   =|(t+Pz)+B_xz|-|t+Pz|.                              \qquad(45a)
```

No scalar phase is discarded. Q_P is a fixed CNOT shear using at
most n-1 gates, and it commutes with $`Z_{z_0}Z_{t_0}`$. Thus this
relative kernel echo costs only two fixed shears beyond (45); in the
commutator source they can conjugate the mixer only if the outer
fixed shears are retained. Explicitly, for the commutator V below,
$`V(Q_PD Q_P,R)=Q_PV(D,Q_PRQ_P)Q_P`$, and the corresponding dual
Gram is conjugated by Q_P as well. Q_P need not commute with the
desired lookup probe, so those outer shears cannot be dropped from
the sector check. This supplies an actual phase insertion sensitive
to the arbitrary kernel translation bank, not a central-character
readout by itself.

### 12.2 Complete normalized source, kernel and ledger

Use the fixed dense mixer from Section 11,
$`R=H^{\otimes2n}U_M`$, $`M(z,t)=(Ez+t,z)`$, with
$`G_R\le4n^2+8n`$. Set
$`F_x=R^\dagger D_xR`$ and
$`V_x=D_x^\dagger F_x^\dagger D_xF_x`$. With the same literal
H-or-X signal word U_0, the physical source is

```math
 \begin{aligned}
 U&=C_E^\dagger V^\dagger U_0VC_E,\qquad
 C_E=U_{E^{-\mathsf T}}\otimes U_{E^{-\mathsf T}},\\
 K&=\frac1{\sqrt2}C_E^\dagger V^\dagger\Pi_0VC_E,\\
 G=K^\dagger K&=\frac12C_E^\dagger V^\dagger\Pi_0VC_E,
 \qquad \Pi_0=(I+Z_{z_0}Z_{t_0})/2.
 \end{aligned}                                             \qquad(46)
```

V is controlled by the full coherent address. Equation (46) proves
normalization and rank on the entire space, with no prepared carrier
or postselection. It also holds when R is replaced by
$`\mathcal R=H_X^{\otimes p}\otimes R`$, acting on the address and
banks together; then V need not preserve the address. The required
Gram is still the fixed-address lookup projector, so intermediate
address motion adds an obligation rather than an automatic return.

For the fixed R, Q=2^(2n), and carrier labels r,s, the full finite sum is

```math
 \begin{aligned}
 (V_x)_{rs}&=Q^{-2}\sum_{t,w,w'}
 \omega^{-h_x(r)-h_x(w)+h_x(t)+h_x(w')}
 (-1)^{w^{\mathsf T}M(t+r)+w'^{\mathsf T}M(t+s)},\\
 (G_{x,\rm dual})_{rs}
 &=\frac12\sum_{u:u_{z_0}=u_{t_0}}\overline{(V_x)_{ur}}(V_x)_{us}.
 \end{aligned}
                                                               \qquad(47)
```

Every sum is over complete arbitrary carrier labels. The exact
unpaid sector identity remains
$`V_x^\dagger Z_{z_0}Z_{t_0}V_x=Z_z(v_x)Z_t(v_x)`$ in the binary
coordinates of Section 11. For an address-moving V, its right-hand
side is the corresponding controlled direct sum and its off-diagonal
address blocks must vanish.

U and its actual inverse use eight D calls, eight R calls, four fixed
carrier-basis changes and one U_0. With $`c_n=n^2+3n`$,

```math
 \begin{aligned}
 T_U&\le128r+16n+2,&d_U&\le146,&b_U&\le r,\\
 G_U&\le16G_B+8(4n^2+8n)+4c_n+19,&\varepsilon_U&=0.
 \end{aligned}                                             \qquad(48)
```

Here b_U excludes the two principal n-bit carriers and counts only
additional query work; embedding the source around the m specified
output wires also allocates the other (2p-1)m carrier wires. Each D
returns its private pool before reuse. The address-Hadamard variant
adds at most 8p Cliffords and nothing to the other bounds. These
polynomial overheads fit the fusion interface, but that interface
cannot repair a wrong projector. There is no approximation or inverse
replacement hidden in (48).

### 12.3 Exact failed sectors and the scope of the small case

At p=1, x=1, z=e_1, t=0, the desired paired character is minus one.
The fixed-mixer word instead has transported probe expectation plus
one and hence the maximal Gram gap one half. The genuine address-H
variant also has expectation plus one on that input, for either value
of the arbitrary native phase helper.

There is a useful analytic scope beyond those two mixers. On the
entire p=1 address-and-carrier space, D from (43) has spectrum
$`\{1,\omega,\omega^{-1}\}`$. For any unit vector b that is a D
eigenvector, and any unitary F on the complete non-signal space,

```math
 |\langle b|D^\dagger F^\dagger DF|b\rangle|
   =|\langle Fb|D|Fb\rangle|\ge1/\sqrt2.                \qquad(49)
```

The numerical range of D lies in the convex hull of those three
points, all with real part at least 1/sqrt(2). Choose b to be the
computational input x=1,z=e_1,t=0 above. Its initial probe is plus one
and its desired probe is minus one. Equation (49) leaves at least
one half of its final probability in the initial plus sector, so
(46) has Gram error at least one quarter. Thus changing only the
inner unitary F, including genuine address motion and arbitrary
spectator mixing, cannot solve this p=1 commutator. This is not an
all-size obstruction: D has a larger spectrum at p at least two,
and an asymptotic family may use a separate small-p construction.

The exact p=2 fixed-mixer fixture uses x=1,z=e_p,t=0. In binary
standard coordinates B_xz is z shifted right by p bits. The target
probe is again minus one, whereas the actual expectation is

```math
 \frac{57}{128}+\frac{33}{256}\sqrt2,
 \qquad
 \|G-G_*\|\ge\frac{185}{512}+\frac{33}{1024}\sqrt2.
                                                               \qquad(50)
```

Changing just the query orientation to B_x^T gives expectation
$`29/128+(21/128)\sqrt2`$ on the same input, also inconsistent with
the target. These are complete-operator witnesses, not average-input
or trace promises. They do not exclude a separate p=2 fallback.

The [exact fixtures](../../tests/test_t_wall_xor_gram_source.py) check
all coherent-address and dirty-input columns of the literal p=1 phase
and its actual inverse, the literal relative kernel-echo identity,
the full one-clean half-projector Gram, the
address-H variant on both dirty-helper inputs, and the p=2 sums.
The latter retain integer coefficients in Z[omega], check exact norm
one, and convert coefficients to unbounded integers before products.
They are finite guards and internal analytic checks, not external
mathematical review or a general native emitter.

The changed phase therefore solves its native carry and normalization
obligations, but not the complete central-character readout. Pause
these particular commutator words; the first unpaid operation is
still the conjugated-probe identity in (47), including address return
when applicable. Continuing this mechanism requires a specified
interference word that changes that identity, rather than further
extraction or precision work. No complete-indicator or frame-depth
frontier moves.

### 12.4 Charge preservation requires a different acceptance polynomial

Fix any nonzero address x. On the two arbitrary carriers, including
identity on any spectators, write
$`P_0=Z_{z_0}Z_{t_0}`$ and $`Q=Z_z(v_x)Z_t(v_x)`$.
These are distinct nonidentity diagonal Paulis. Put
$`\Pi_0=(I+P_0)/2`$, $`\Pi_{\rm eq}=(I+Q)/2`$ and
$`\Pi_{\rm neq}=(I-Q)/2`$. For any unitary V on this fixed-address
space that commutes with Q, Pauli trace orthogonality gives

```math
 \mathrm{Tr}(\Pi_{\rm neq}V^\dagger\Pi_0V)
   =\mathrm{Tr}(\Pi_{\rm neq}\Pi_0)=D_{\rm tot}/4,
 \qquad
 \left\|\tfrac12V^\dagger\Pi_0V-\tfrac12\Pi_{\rm eq}\right\|
   \ge\frac14,                                         \qquad(51)
```

where D_tot is the full dimension including the spectators. The
target vanishes on the unequal sector, while the actual positive
Gram has mean eigenvalue one quarter there. In fact the displayed
distance is exactly one half: V-dagger Pi_0 V commutes with Q, and
its unequal-sector restriction is a projector of rank D_tot/4, so
the actual Gram has eigenvalue one half where the target is zero.
Both Grams are half-projectors, so their distance is at most one half.

Kernel translations and bank-diagonal phases preserve Q. A word
using only such charge-preserving operations therefore cannot obtain
the desired Gram by conjugating the fixed Pi_0 in (46). This is a
statement about that exact fixed-address conjugation template, not
an unrestricted depth bound, an address-motion exclusion, or an
obstruction to every source that preserves Q.

The next concrete question is to build a **signal-interleaved** source
from the paid corank-one query $`Q_{K_x}`$ and its relative XOR-carry
phases: specify the signal couplings and their interleaving, derive
the actual block $`K=J^\dagger UJ`$, and test its full Gram rather than
conjugating the same fixed projector. Any signal-controlled query or
phase occurrence needs its own literal native schedule, actual
inverse and arbitrary-input return proof before a fixture treats it
as available. No qualifying interleaved word is claimed here; the
first unpaid operation remains a normalized complete lookup Gram.
