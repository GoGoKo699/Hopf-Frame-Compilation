# Two-layer transport of one orbit direction and its composition limit

[MOD3 transport](MOD3_ZETA_TRANSPORT.md#pauli-transport-interface) ·
[Working checkpoint](../../WORKSPACE.md) ·
[Route decisions](../ROUTE_DECISIONS.md)

The complete transport contract permits a different source representative
when its relative error commutes with the input Pauli. This freedom does
give a small native construction: one controlled swap transports the
specified Pauli through any one signed MOD3 coordinate. It does not give
the complete address orbit. The first missing identity is preservation
of the directions produced by earlier address coordinates.

<a id="one-coordinate-compression"></a>

## 1. A complete two-layer source for one coordinate

Use the normalized binary coordinates of the
[MOD3 reduction](MOD3_ZETA_TRANSPORT.md#order-three-blocks). There are
$`p\ge1`$ address bits, $`m=2^p`$ arbitrary bank bits in
$`q=m/2`$ pairs, and sign representatives $`\sigma_j\in\{+1,-1\}^p`$
whose first sign is positive. Put the all-positive representative first.
For coordinate $`i\in\{0,\ldots,p-1\}`$, define

```math
F_i=\bigoplus_{j=0}^{q-1}A^{\sigma_{j,i}},\qquad
v=\bigoplus_j(1,0),\qquad w_i=F_iv,\qquad P_0=X^v.
```

Thus each block of w_i is (0,1) at a positive sign and (1,1) at a
negative sign. The pivot block j=0 has positive sign for every i.
Let C_i be the following fixed CNOT word on the bank:

- fan the pivot first bit into every other first bit;
- fan the pivot second bit into every other second bit;
- additionally fan the pivot second bit into every negative-sign first bit.

All controls lie in the unchanged pivot pair and all targets lie outside
it. These CNOTs commute, $`C_i^2=I`$, and, as binary maps,

```math
C_i e_{0,0}=v,\qquad C_i e_{0,1}=w_i.
```

Let J swap the two pivot wires and let $`R_i=C_iJC_i`$. Emit C_i,
a swap of the pivot pair controlled by x_i, and C_i, in that chronological
order. The resulting full unitary E_i has bank action $`R_i^{x_i}`$.
In particular,

```math
E_iP_0E_i^\dagger
 =X^{F_i^{x_i}v}\otimes I_d.                         \qquad(1)
```

Here d is the arbitrary private helper used by the native swap below.
Equation (1) covers every coherent address, bank, helper and reference
input. Other address bits are arbitrary spectators. It agrees with the
complete MOD3 target on the address face where only x_i may be nonzero;
it is a partial transport statement for general p. For p=1 that face is
the entire address space, so it supplies the complete weak source.

The operator action can be written without a circuit diagram. For a
binary bank column z put $`\ell(z)=z_{0,0}\oplus z_{0,1}`$ and
$`u_i=v+w_i`$. Then

```math
R_i z=z+\ell(z)u_i,\qquad \ell(u_i)=0.                \qquad(2)
```

Thus R_i is an involution, its binary displacement rank is one, and it
exchanges v and w_i. The source generally changes bank columns; its
Pauli echo has precisely the full operator action in (1).

### Native schedule, precision and actual inverse

The retained [two-layer dirty-helper Toffoli](DIAGONAL_ENCODING_SESSION_RESULT.md#2-constructive-check-and-first-unpaid-cost)
uses eight T/T-dagger gates, sixteen CNOTs and two H gates. Its one
helper d is arbitrary and returns exactly, including literal phase.
For pivot bits a,b the chronological word

```text
CX(a,b), CCX(x_i,b;a), CX(a,b)
```

is controlled swap. Hence it costs eight T gates, twenty elementary
Cliffords and two physical T layers, with one arbitrary returned helper.
The Clifford walls contain
$`c_i=2(q-1)+n_i`$ CNOTs, where n_i is the number of negative signs.
Here $`n_0=0`$, and $`n_i=q/2`$ for i>0. The common forward and
actual-reverse source ledger is

```math
T_{E_i}=8,\qquad D_T(E_i)=2,\qquad
G_{E_i}=20+2c_i=4q+16+2n_i,\qquad w_{E_i}=1.           \qquad(3)
```

There are zero clean inputs and zero implementation error. The actual
reverse reverses every native gate and changes T to T-dagger and
conversely; no ideal involution is substituted for an emitted inverse.
Both physical T bands act on four distinct wires. Using
$`X=HS^2H`$ on the q fixed first bits, the emitted observable echo has
T-count 16, T-depth at most four, Clifford count $`2G_{E_i}+4q`$,
and the same one returned helper. Signal control of this observable
requires only q CNOTs in place of the middle Pauli; E_i remains
unconditional.

This construction does not contradict the
[order-three spectrum obstruction](MOD3_ZETA_TRANSPORT.md#scoped-mod3-limits).
That obstruction concerns conjugating a full order-three unitary into
copies of itself with identity helpers. The source here matches one
Pauli image for the two values x_i=0,1 and uses a different bank unitary.

<a id="composition-failure"></a>

## 2. Exact one-coordinate identities do not compose

Consider the candidate $`E=E_{p-1}\cdots E_0`$, with the same helper
reused after its exact return. Every source has the correct image of
P_0 individually. After the first source, however, the next source must
act on the already transported Pauli, rather than on P_0 again.
Equation (2) makes this failure explicit. All the R_i commute, since
every u_i is in the kernel of the same linear functional ell. Therefore

```math
R_{p-1}^{x_{p-1}}\cdots R_0^{x_0}v
 =v+\sum_i x_i(v+w_i).                                \qquad(4)
```

This mask is affine in the address. The required mask
$`F(x)v=\prod_iF_i^{x_i}v`$ need not be affine.
For every p at least two, set x_0=x_1=1 and all other address bits zero.
The two masks on a block j are then

| Sign $`\sigma_{j,1}`$ | Candidate mask | Required mask |
| --- | --- | --- |
| +1 | (1,0) | (1,1) |
| -1 | (0,0) | (1,0) |

Both sign classes are present. The candidate and target are distinct
commuting X Paulis, whose operator distance is two. Consequently

```math
\|EP_0E^\dagger-\mathcal B_pP_0\mathcal B_p^\dagger
                    \otimes I_d\|=2\qquad(p\ge2).      \qquad(5)
```

The equality follows from this address block and the general upper bound
two for the difference of unitaries. Private-helper return is already
exact. Additional precision cannot repair this exact sector mismatch.

<a id="prefix-orbit-rank"></a>

## 3. A prefix-orbit rank requirement for fixed coordinate factors

There is a stronger necessary condition than checking one Pauli at each
stage. Consider arbitrary fixed invertible binary maps V_0,...,V_{p-1}
on the bank and any auxiliary bits, with source action
$`V_{p-1}^{x_{p-1}}\cdots V_0^{x_0}`$. Active factors may move all
those bits. The address is the read-only selector of this stated
factorization, and a factor is identity when its selector bit is zero.
Assume the source transports v, with zero auxiliary mask, to F(x)v
with zero auxiliary mask for every address x. Then

```math
\mathrm{rank}(V_i-I)\ge2^i,\qquad
\sum_{i=0}^{p-1}\mathrm{rank}(V_i-I)\ge m-1.             \qquad(6)
```

To prove this, use $`F(x)v=Se_x`$ from the fixed binary reduction.
Setting every later address bit to zero forces the completed prefix
to send v to $`Se_y`$ for each $`0\le y\lt2^i`$. Turning on x_i
then requires

```math
(V_i-I)(Se_y)=S(e_y+e_{y+2^i})\qquad(0\le y\lt2^i),      \qquad(7)
```

with all vectors embedded by zero auxiliary coordinates. These $`2^i`$
differences are independent: their two-element supports are disjoint
before the invertible S wall. This proves (6). No commutation assumption
on the V_i is used. In the phase-free binary-permutation grammar, an
observable error strictly below two forces the same mask identity on
each address, so the same necessary rank bounds apply.

In particular, any binary conjugate of a single swap has displacement
rank one. Choosing a different C_i cannot make p such single-swap
representatives compose into the full transport when p is at least two:
the second coordinate already requires rank at least two. A stage
conjugate to k_i disjoint swaps requires $`k_i\ge2^i`$.

The label argument also extends to fixed completed controlled-Clifford
factors on the bank and arbitrary work. Use their full symplectic
matrices in (6) and embed each target X label in the symplectic space.
The same prefix differences are independent. An observable error
strictly below $`\sqrt2`$ forces equality of the Hermitian Pauli images,
since distinct signed Pauli operators have distance at least sqrt(2).
Literal phases remain an additional obligation. A physical swap has
symplectic displacement rank two, so this extension has a different
active-swap count than the binary restriction.

This is a rank requirement for the specified completed factorization.
It is not a native T-depth lower bound: the retained
[shared-root CCX batch](MOD3_ZETA_TRANSPORT.md#3-a-paid-linear-depth-native-baseline-and-optional-signal-control)
implements arbitrarily many controlled disjoint CNOTs in four T layers,
while the active bank map's displacement rank grows with their number.
The available width permits m directions. General quantum sources need
not have these fixed controlled-Clifford factors or preserve the address
at their boundaries.

## 4. Research consequence and evidence

The two-layer component pays weak transport of a single orbit direction.
The unrestricted two-band screen for the complete source, at p at least
three and the stated fixed error threshold, remains unchanged. The
first unsupported composition is replacing each coordinate action by
a representative certified only on the original Pauli. A successful
joint construction must transport the growing prefix orbit or leave this
completed-factor representation.

Pause the product of single-swap representatives. Reopen it only with
paid action on the whole required prefix orbit, or change to a source
whose non-Clifford interference supplies the joint observable directly.
Neither the component nor the rank requirement moves the complete-frame
resource frontier. The known logarithmic-depth complete indicator remains
the baseline, and the approximate scalar zero-reflection result continues
to rule out using conjunction alone as a growing fixed-accuracy witness.

The [exact fixtures](../../tests/test_weak_pauli_orbit_compression.py)
check the native source, literal Pauli echo, private helper and actual
reverse on every column in small instances. Binary fixtures check all
addresses through p=7 and the independent prefix differences. They
corroborate the analytic identities; no general quantum compiler,
unrestricted depth lower bound or external mathematical review is claimed.
