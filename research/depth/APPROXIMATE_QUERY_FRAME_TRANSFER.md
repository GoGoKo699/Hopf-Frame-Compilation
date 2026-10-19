# Constant-accuracy readouts and the complete frame transfer

[Exact dual-readout barrier](DUAL_FRAME_READOUT_BARRIER.md) ·
[Complete word-query assembly](SYMPLECTIC_MODULE_QUERY_BARRIER.md#7-the-complete-r98-word-query-interface-follows-conditionally) ·
[Fixed-precision frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md) ·
[Paid Majorana source](MAJORANA_WALSH_QUERY.md)

The [joint-indicator continuation](APPROXIMATE_INDICATOR_FILTER_LIMITS.md)
tests additive polynomial phases and coherent short hashes against this
full source contract. Their complete-operator errors do not supply a
faster query; the constant-accuracy transfer proved here remains valid.

Approximate readout is a real remaining source question: the complete
R98 construction needs only constant-in-address readout accuracy at
fixed final precision. Neither the number of table rows, output batches
nor catalogue entries multiplies the source-error budget. This chapter
proves that transfer and a robust limitation on one candidate compiler.
A separately paid Majorana-preserving conjugator still needs
$`\Omega(m\log m)`$ T gates and at least $`\log_2m`$ T layers at
small constant readout error. A non-Gaussian or fused implementation is
outside that limitation. No faster native source or general-frontier
improvement is supplied here.

## 1. The complete approximate source hypothesis

Use the notation of the dual-readout chapter. For $`p\ge2`$, let
$`m=2^p`$, $`w=m/2`$, and let Y be an arbitrary m-bit register.
The carrier A has w arbitrary qubits. The canonical target is

```math
Q_p=\sum_Y|Y\rangle\langle Y|\otimes WC_YW^\dagger,
\qquad C_Y=\Delta_0^{Y_0}\cdots\Delta_{m-1}^{Y_{m-1}},
```

where the positive Hermitian physical Pauli duals and the matrix-product
order are fixed as in that chapter. Its six-call interface also permits
the target $`\widetilde Q_p=e^{i\theta_p(Y)}Q_p`$ with a fixed
completed diagonal Y-phase. It does not permit an arbitrary right
register unitary in place of this scalar gauge.

For each requested tolerance $`0\lt\tau\le1/6`$, suppose a
terminating classical compiler supplies a native unitary $`V_{p,\tau}`$
with no initialized work and b_p arbitrary dirty helpers, such that

```math
\|V_{p,\tau}-\widetilde Q_p\otimes I_{b_p}\|_{\mathrm{op}}
\le\tau.                                                       \qquad(1)
```

This is a full-operator promise on every Y/carrier/helper input, hence
also on all their superpositions and external references. Return of Y
and work is included in the norm; it is not inferred by testing basis
probabilities or discarding work. The actual inverse word is supplied.
Assume common forward/actual-inverse bounds

```math
\begin{aligned}
T(V_{p,\tau})&=O_\tau(m),&G(V_{p,\tau})&=O_\tau(m^2),\\
w+b_p&=O_\tau(m),&D_T(V_{p,\tau})&=O_\tau(d_\tau(p)),
\end{aligned}                                                   \qquad(2)
```

where $`d_\tau\ge1`$ is nondecreasing in its address argument.
Constants may depend on the requested tolerance, not on p or a table.
In particular, this is a uniform requested-accuracy interface: a larger
address does not license worse accuracy at the smaller half-addresses
later used inside a query. The inverse Clifford ledger must be paid
separately when the positive-S alphabet makes it longer.

The exact denominator obstruction does not exclude (1)--(2). Its
exponentially small coefficient is not a constant-error witness.

## 2. From six readout calls to one complete word query

Let I_p be the exact dirty indicator that adds the selected one-hot word
to its arbitrary output register. The dual-readout chapter's actual
six-call word uses $`\widetilde Q_p`$ and its adjoint, four coherent
Clifford masks, four X words, a parity correction and output Hadamards.
The completed diagonal scalar cancels. Replacing each source by
$`V_{p,\tau}`$ or its actual adjoint gives

```math
\|\widetilde I_p-I_p\otimes I\|_{\mathrm{op}}\le6\tau.
                                                               \qquad(3)
```

Unitary product telescoping proves (3); actual adjoints have the same
error as their forward words. Its common compiled resource bounds are

```math
\begin{aligned}
T(\widetilde I_p)&\le6T_V,&D_T(\widetilde I_p)&\le6d_V,\\
G(\widetilde I_p)&\le6G_V+20pw+4m+16,&
b(\widetilde I_p)&\le w+b_p.
\end{aligned}                                                   \qquad(4)
```

The m output wires are excluded from this helper count. They must be
included when the indicator is used as a borrowed bank inside a query.
For zero or one address bit, use the exact Clifford indicator directly;
it needs no approximate readout.

Now apply the R91 word-query assembly with Q rows, q output bits and
extra dirty allocation B. Its source-call audit is:

| Branch | Completed indicator appearances | Error bound |
| --- | --- | --- |
| Q equals one | None; fixed Clifford row | 0 |
| Direct full indicator | One forward and one actual inverse | $`12\tau`$ |
| Bilinear allocation with lengths H,J | Two H-indicators and two J-indicators | $`24\tau`$ |

Clifford zero/one-bit indicators can reduce these counts. In the last
row both pairs lie outside every traversal and output batch; no factor
$`\lceil q/g\rceil`$, K, Q or q multiplies the error. More generally,
indicator errors $`\epsilon_h,\epsilon_j`$ give
$`2\epsilon_h+2\epsilon_j`$. The common tolerance in (1) controls both.
Every complete call is sequential. Internal changes of a shared address
therefore require no unproved parallel interleaving rule.

The simultaneous R91 allocation is unchanged after enlarging its
indicator constant to include the carrier and helper pool in (4).
The new Clifford term is $`O_\tau(H^2+J^2)`$, bounded by
$`O_\tau(Qq)`$ in every allocation branch by the retained rounding
audit. The direct branch has $`Q^2\le Qq`$. Consequently, for
$`B\ge C_\tau(a+2)`$, where $`Q=2^a`$, the complete approximate
word query has

```math
\begin{aligned}
T&=O_\tau(\sqrt{Qq}+Qq/B),\qquad G=O_\tau(Qq),\\
D_T&=O_\tau(d_\tau(a)+Qq/B^2),\\
\|\widetilde U_f-U_f\otimes I\|_{\mathrm{op}}&\le24\tau.
\end{aligned}                                                   \qquad(5)
```

Its physical allocation is at most B, with zero clean work. No
independence of different approximation errors is assumed. All address,
output and dirty-return errors are included in (5), and its actual
inverse obeys the same estimate and common ledger.

## 3. The bounded R98 call count and final error

Replace only R98's completed long queries by (5). Keep its short
standalone tables, paired layouts, local sources and terminal fixed words
at their retained implementations. Let k_i be the number of completed
joint long-query occurrences in emitted early stage i, counting every
echo, recovery occurrence and actual inverse. The fixed stage grammar
and joint-word packing in R98 Section 2 give $`k_i\le C_E`$ for
one fixed constant C_E. This is the count after combining the former
polynomially many shared-address instances into one word query; it is
not the count of those individual instances. There are at most four
early stages.

The other long-query occurrences have explicit bounds:

| Position | Completed word-query occurrences |
| --- | ---: |
| Two initial Boolean sign phases | 4: two query/actual-inverse echoes |
| Four early stages, including all recovery | $`\sum_i k_i\le4C_E`$ |
| Whole terminal | 32 |
| Complex literal diagonal, when present | 20 |

For the terminal, the retained unreduced word has four commutator
factors, two selected swaps in each factor, two translations per swap,
and two query occurrences per translation: $`4\cdot2\cdot2\cdot2=32`$.
There is no repetition for each catalogue entry. Trivial addresses,
an omitted terminal and actual word cancellations can only reduce this
count. The complex diagonal's twenty includes both orientations.
Thus a uniform upper bound for either family is

```math
K_* =4C_E+56,
\qquad \hbox{number of approximate readout calls}\le24K_*.
                                                               \qquad(6)
```

C_E is fixed by the chronological construction, independent of n and
b. An emitted compiler may instead count its actual query occurrences
and use that smaller K. All of these are constant counts, not unspecified
factors growing with the address, precision bank or catalogue size.

Fix the final precision L, put $`\eta=2^{-L}`$, and request the same
source tolerance at every length:

```math
\tau=\frac{\eta}{96K_*}.                                    \qquad(7)
```

The added full-word error is at most $`24K_*\tau=\eta/4`$.
R98's wide real branch has retained error $`21\eta/512`$, so the
new error is at most $`149\eta/512\lt\eta`$. For the complex
branch, its real component runs at L+1 and has retained error at most
$`21\eta/1024`$; its literal diagonal contributes at most
$`\eta/2`$. Including all approximate queries gives at most
$`789\eta/1024\lt\eta`$.

These are full initialized-isometry comparisons on the frame's sole
clean flag, with arbitrary dirty inputs and references. A hybrid compares
each substituted query to its full ideal operator. Later actual
unitaries have norm one on any previously leaked components; no work
is assumed freshly initialized at a reuse boundary. Exact intermediate
return is not required of the approximants. The final bound includes
all residual address, dirty and flag error.

For completeness, a source supplied on an initialized subspace also
admits this error argument: if J embeds that subspace and
$`\|VJ-JU\|\le\epsilon`$, then
$`\|V^\dagger J-JU^\dagger\|\le\epsilon`$. The same unitary hybrid
applies when ideal intervening words preserve that subspace. Such a
source would add its initialized allocation; it is not the zero-clean
hypothesis (1)--(2) and does not automatically satisfy R98's one-clean
total budget.

At fixed L, (7) is a fixed positive constant independent of n,b.
All constants in (2), (4) and (5), including the minimum query pool,
can therefore be absorbed into R98's L-dependent stage choices and
large-dimension cutoff. Its small-dimension and narrow-width branches
retain R96 unchanged, so no stronger accuracy or width premise is
imposed at the literal thresholds.

Under (1)--(2), the resulting conditional bound is consequently

```math
\begin{aligned}
T&=O_L(\sqrt N+N/b),\qquad G=O_L(N),\\
D_T&=O_L\!\left(N/b^2+d_{\tau_L}(2n+10)
                         +\log_2\log_2(n+4)\right),
\end{aligned}                                                   \qquad(8)
```

with R98's same one-clean and family-specific dirty thresholds and full
error $`2^{-L}`$. A source with $`d_\tau(p)=o(\log(p+2))`$ for
each fixed requested tolerance would therefore suffice. No uniform
growing-L theorem is asserted. Primitive-by-primitive approximation
inside a proposed source may still incur an address-dependent precision
cost; (7) does not pay or remove that internal compilation cost.

## 4. Constant error does not repair the Gaussian conjugator

Consider a specific candidate architecture

```math
Q'_Y=V(C_Y\otimes I)V^\dagger,                          \qquad(9)
```

where V acts on the carrier and any number of arbitrary dirty helper
modes, excluding Y. Require V to have a native implementation whose
every intervening Clifford block preserves the span of all M physical
Majoranas, interspersed with physical T/T-dagger bands. A Clifford block
in this model is a signed permutation of those Majoranas. Merely having
a Gaussian endpoint is insufficient for this hypothesis. The fixed
controlled word C need not be Gaussian on the joint Y/carrier system.

The completed readout is the direct sum of (9) over Y. Suppose it has
full-operator distance at most
$`\tau\lt1/3`$ from the canonical readout, allowing a diagonal
Y-dependent scalar gauge and requiring identity on its dirty helpers.
Put $`\delta=3\tau/2\lt1/2`$ and let h_2 denote binary entropy.
Then its separately implemented conjugator satisfies

```math
T(V)\ge\frac m2\big[(1-\delta)p-h_2(\delta)\big],
\qquad D_T(V)\ge p.                                      \qquad(10)
```

For example, at $`\tau\le1/6`$ the count is at least
$`\frac m2[3p/4-h_2(1/4)]`$. This bound permits arbitrary intermediate
mixing with dirty helper modes. It is not an unrestricted approximate-Q
bound or a bound on a fused query that never implements V separately.

### Pair sectors recover the individual lines

Let O be V's real orthogonal Majorana matrix, let v_i be its ith
original-mode column, and let h_i be the corresponding Walsh column
padded with zeros on helpers. Choose any distinct original i,j.
In the Y sector supported on that pair,
$`C_Y=\pm\gamma_i\gamma_j`$. It is even, so its conjugation flips
exactly those two modes, including when helper modes are present.
The two Majorana conjugation matrices in (9) and its target are

```math
I-2P_{ij},\quad I-2P^H_{ij},\qquad
P_{ij}=v_iv_i^{\mathsf T}+v_jv_j^{\mathsf T},\quad
P^H_{ij}=h_ih_i^{\mathsf T}+h_jh_j^{\mathsf T}.
```

Unitary distance at most $`\tau`$ implies adjoint-action distance at
most $`2\tau`$; sector scalars cancel. For real v, the Majorana combination
$`\gamma(v)`$ has operator norm $`\|v\|_2`$. Hence the preceding
adjoint estimate is an operator-norm matrix estimate, giving
$`\|P_{ij}-P^H_{ij}\|\le\tau`$.

There are at least three original modes because $`p\ge2`$. For
distinct i,j,k, the identity

```math
v_iv_i^{\mathsf T}=\tfrac12(P_{ij}+P_{ik}-P_{jk})
```

and its Walsh counterpart give

```math
\|v_iv_i^{\mathsf T}-h_ih_i^{\mathsf T}\|
\le\delta,
\qquad |\langle v_i,h_i\rangle|^2\ge1-\delta^2.        \qquad(11)
```

Thus a readout approximation already forces every original column to
be close to its flat Walsh line. It is unnecessary to assume that V
approximates W itself, or to choose a common scalar for different sectors.

### Entropy and support remain large

The exact-source proof applies
[Ailon's matrix-entropy method, Theorem 2.1](https://arxiv.org/html/1305.4745v1).
Here is the robust step. Choose the sign of h_i to make its inner
product with v_i nonnegative. The total variation distance between
their squared-coordinate distributions is at most delta, since

```math
\tfrac12\sum_s|(v_i)_s^2-(h_i)_s^2|
\le\tfrac12\|v_i-h_i\|_2\|v_i+h_i\|_2
\le\delta.
```

Coarse-grain every helper coordinate into one extra outcome. The target
distribution is uniform on the m original modes. Coupling the two
distributions with mismatch probability at most delta gives an entropy
loss of at most $`h_2(\delta)+\delta\log_2m`$. Refining the helper
outcome cannot decrease entropy. Each original column therefore has
entropy at least $`(1-\delta)p-h_2(\delta)`$.

Sum these entropies over the m original columns. Signed row permutations
preserve the sum. A physical T mixes two rows; in each column their
squared mass is preserved and their entropy changes by at most that
mass. Summing over the selected columns gives at most two, since their
total mass in those rows is at most the full orthogonal rows' mass.
The identity starts at zero. This proves the count in (10), even with
any number of helper columns.

After D physical T bands, signed permutations and matching two-row
rotations give at most $`2^D`$ nonzero entries in each column. A unit
vector supported on s coordinates has squared overlap at most s/m
with h_i. Equation (11) therefore gives

```math
2^D\ge m(1-\delta^2).
```

For $`\delta\lt1/2`$ the right side exceeds m/2. Since
$`m=2^p`$ and D is an integer, $`D\ge p`$. This proves (10).

## 5. Remaining operation and verification scope

The norm transfer shows that an approximate source need not resolve
the tiny exact channel coefficient. The robust argument shows why
approximating the existing Gaussian conjugator does not supply that
source. A non-Gaussian complete compiler or a fused word remains a
different possibility, with all simultaneous native costs still unpaid.

Moving the Walsh transform to a small index register does not by itself
pay the operation either. The usual uniform-index sandwich around
$`\mathrm{SELECT}(\gamma_j)`$ has clean-index block
$`m^{-1}\sum_j\gamma_j=A/\sqrt m`$. It requires the initialized
index and the coherent selected-Pauli operation, and its normalization
is not the unitary A. This is a diagnosis of that proposed interface,
not a lower bound on every interference-based compiler.

The all-size transfer and Gaussian limitation above are analytic. The
[focused fixture](../../tests/test_approximate_query_frame_transfer.py)
uses four original modes and one arbitrary helper qubit. A native
Gaussian perturbation E genuinely mixes carrier and helper modes and
has the exact full-operator certificate

```math
(E-I)^\dagger(E-I)=(23/16-\sqrt2)I,
\qquad \tau^2\le23/4-4\sqrt2\lt1/9.
```

The compiled conjugator $`V=E(W_4\otimes I)`$ has 26 T and 534
elementary Clifford gates; its complete controlled readout has 52 T
and 1086 Cliffords. Both actual inverse words have the same counts.
The fixture checks all six pair sectors, their nonreal scalar phases
and actual adjoints, recovery of the individual line projectors, and
exact probability/support estimates with occupied helper inputs. It
checks finite native geometry, not the analytic constant-call composition
theorem, an asymptotic approximate source or an improved frontier.
