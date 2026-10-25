# Parallel cubic phases with returned dirty work

[R91 phase derivative](PARALLEL_BILINEAR_WORD_QUERY.md#1-a-phase-derivative-on-an-arbitrary-dirty-wire) ·
[Complete source contract](../../WORKSPACE.md#main-unfinished-route-and-exact-missing-operation) ·
[Coefficient shears](WEAK_PAULI_ORBIT_COMPRESSION.md#complete-coefficient-shear)

A diagonal weighted-cubic phase has a constant native T-depth compiler
even when its monomials share arbitrary logical controls. Private dirty
phase targets supply the parallelism; an exact derivative returns them.
The construction below assembles the existing R91 primitive. It is a
paid component, not a complete lookup source or a general depth improvement.

For related catalytic synthesis, [Folkertsma, Mertz, Strelchuk and Subramanian,
arXiv:2609.40144v1](https://arxiv.org/html/2609.40144v1), Lemmas 18/20 and
Theorem 23, restore arbitrary work using finite-field phase cycles at fixed
hierarchy level. Their total-depth theorem and this explicit native T-band
ledger measure different resources. The present construction assembles the
retained R91 Gray-word primitive; no priority claim for catalytic diagonal
synthesis is made. [F52](../../docs/reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries)
records the comparison and its growing-level limitation.

## 1. A complete diagonal compiler

Let z contain n arbitrary logical bits. Keep a list of Q signed quadratic
occurrences and C cubic occurrences, and define

```math
P(z)=\sum_{v=0}^{n-1}a_vz_v
       +2\sum_{j=1}^{Q}\sigma_jz_{u_j}z_{v_j}
       +4\sum_{j=1}^{C}z_{r_j}z_{s_j}z_{t_j}\pmod8,
\qquad \sigma_j\in\{+1,-1\}.                         \qquad(1)
```

Within each occurrence the indices are distinct. Repeated indices can
first be reduced using $`z_v^2=z_v`$. Repeated occurrences remain allowed,
and controls may be shared between any number of terms. In particular an
even quadratic coefficient is represented by repeated signed occurrences;
it is not silently removed from the gate ledger. Every weighted cubic
polynomial with even quadratic and four-divisible cubic coefficients has
this form. There is no constant term; any specified scalar is a separate
operation.

There is a uniform native compiler for

```math
U_P|z,w\rangle=e^{i\pi P(z)/4}|z,w\rangle
```

with Q first dirty helpers and Q+C second dirty helpers. Its emitted word
and actual native inverse obey

```math
\begin{aligned}
T&\le n+12Q+8C,&D_T&\le13,\\
G_{\rm all}&\le4n+24Q+16C,&w&=2Q+C,\\
a_{\rm clean}&=0,&\epsilon&=0.
\end{aligned}                                                   \qquad(2)
```

Here $`G_{\rm all}`$ counts **all** elementary H, S, CX, T and T-dagger
gates. It also upper-bounds the Clifford-only gate count used in the
general resource contract. The width w excludes the n logical inputs.
Every dirty input returns, including superpositions and correlations with
the logical inputs or any external reference.

### The four-band primitive and its signed use

For distinct u,v,h, the R91 identity is

```math
\begin{aligned}
p(h;u,v)&=h-(u\mathbin\oplus h)-(v\mathbin\oplus h)
                         +(u\mathbin\oplus v\mathbin\oplus h)
                         =(4h-2)uv,\\
F(u,v;h)&=e^{i\pi p(h;u,v)/4}
          =\mathrm{CS}^{\dagger}(u,v)\mathrm{CCZ}(u,v,h).
\end{aligned}                                                   \qquad(3)
```

Use the closed Gray walk h, h+u, h+u+v, h+v, h. Apply respectively
T, T-dagger, T, T-dagger at its first four vertices; one CNOT makes
each move, including the return to h. This literal word uses four T
gates and four CNOTs. For a list of private h targets, make each move
for every target, then apply the corresponding phase gates on distinct
targets. Shared logical controls therefore do not multiply the four
T bands. Changing all four phase signs emits $`F^{\dagger}`$ exactly.

For quadratic occurrence j, emit $`F(u_j,v_j;h_j)^{-\sigma_j}`$.
It contributes the desired $`\mathrm{CS}^{\sigma_j}`$ and one
$`\mathrm{CCZ}(u_j,v_j,h_j)`$ correction. The correction is the same
for both signs, since CCZ is an involution. Every h is unchanged; its
temporary input-dependent phase will be cancelled below. This stage
uses 4Q T gates, 4Q CNOTs and at most four T bands.

### Cancel the corrections together with the requested cubics

Make a list of Q+C triples: the Q correction triples $`(u_j,v_j,h_j)`$,
followed by the C requested triples. Give each triple (u,v,c) a fresh
private k. For one triple use the closed three-dimensional Gray walk.
Its visited control subsets before returning to the empty subset are

```math
\begin{gathered}
\varnothing,\ \{u\},\ \{u,v\},\ \{v\},\\
\{v,c\},\ \{u,v,c\},\ \{u,c\},\ \{c\}.
\end{gathered}                                                   \qquad(4)
```

At subset J the private target contains k plus the XOR of those control
bits. Apply T for even-sized J and T-dagger for odd-sized J. The complete
integer phase exponent satisfies

```math
\sum_{J\subseteq\{u,v,c\}}(-1)^{|J|}
       \left(k\mathbin\oplus\bigoplus_{j\in J}j\right)
       =(8k-4)uvc=4uvc\pmod8.                        \qquad(5)
```

If a control is zero, the summands cancel in pairs. If all three are one,
the four even-subset values are k and the four odd-subset values are 1-k,
giving 8k-4. Thus (5) implements CCZ(u,v,c), literally independent of k;
the final Gray move returns k.
Make corresponding moves and phase gates for all triples together.
No k occurs as a control, and all controls remain unchanged, including
any earlier dirty h. Thus the product implements all requested cubics
and cancels every earlier h-dependent correction. There are eight
physical T bands, on distinct k targets. Each triple costs eight T
gates and eight CNOTs, hence sixteen elementary native gates. Arbitrary
overlap among logical controls and earlier h controls does not change
this schedule.

Finally, each linear coefficient modulo eight is one T or T-dagger at
most, together with at most three S gates. All linear T gates fit one
band. This proves (2). Reversing the emitted word, swapping T and
T-dagger and replacing S inverses by positive S powers, gives the actual
inverse. Reduce each adjacent inverse S block modulo four; it still costs
at most three S gates per logical wire. The common bounds in (2) follow
without identifying an inverse with its forward circuit.

## 2. A paid parity class of coefficient-cube edges

Let $`p\ge2`$, $`m=2^p`$, and index the arbitrary bank bits $`A_S`$
by subsets of the p address coordinates. Put $`n_T=(I-Z_{A_T})/2`$.
For one fixed parity $`\rho\in\{0,1\}`$, consider the edge set

```math
\mathcal E_\rho=\{(S,i):i\notin S,
                                  \ |S|\equiv\rho\pmod2\},
\qquad r=|\mathcal E_\rho|=pm/4.                      \qquad(6)
```

Every tail has parity rho and every head has the opposite parity. Hence
the axes $`X_{A_S}n_{S\cup\{i\}}`$ commute, even for edges sharing a
tail or a head. The address-signed layer is

```math
L_{\rho,x}=\prod_{(S,i)\in\mathcal E_\rho}
     \exp[-i\pi Z_{x_i}X_{A_S}n_{S\cup\{i\}}/4].     \qquad(7)
```

Apply H to the parity-rho bank wires. In the resulting diagonal basis,
one edge with tail bit a and head bit b has exponent

```math
-Z_{x_i}Z_a b=-b+2x_i b+2ab-4x_iab.                  \qquad(8)
```

Summing these integer exponents gives (1), with at most m active linear
variables, Q=2r signed quadratic occurrences and C=r cubic occurrences.
Keep the occurrence lists before collecting coefficients, so coincidences
of terms do not omit any charge. The two H walls together cost m gates.
Using the slightly larger bound 2m for those walls gives convenient
common forward/actual-inverse bounds

```math
\begin{aligned}
T(L_{\rho,x})&\le m+32r,&D_T(L_{\rho,x})&\le13,\\
G_{\rm all}(L_{\rho,x})&\le6m+64r,&w&\le5r,\\
a_{\rm clean}&=0,&\epsilon&=0.
\end{aligned}                                                   \qquad(9)
```

The unsigned layer with positive angle pi/4 has per-edge polynomial
$`b-2ab`$ and fits the same upper bounds. All addresses, bank inputs and
helpers are arbitrary and returned except for the specified bank action
of the layer itself. The signed layer preserves the address exactly.

This pays a new scheduling choice: all edges in one bipartition class
can be compiled together at constant T-depth and polynomial overhead.
It does not authorize regrouping a coordinate-ordered source. Opposite
parity classes need not commute, and their complete lookup action must
be established separately. A word made from these paid layers still
needs its own complete operator and error/return proof.

The [exact fixtures](../../tests/test_dirty_cubic_phase_compiler.py) check
the primitive, full small logical/dirty matrices, shared controls,
actual native inverses, coherent helper correlations and distinct
physical T targets. These finite checks audit the compiler; (3)--(5)
give the all-size identity. No external priority or technical review is
claimed for this assembly of the retained primitive.
