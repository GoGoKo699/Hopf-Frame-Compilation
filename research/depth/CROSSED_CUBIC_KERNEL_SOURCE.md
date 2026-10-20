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
