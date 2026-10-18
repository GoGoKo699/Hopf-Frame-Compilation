# Exact phase relocation onto arbitrary dirty work

[Four-layer query frontier](FOUR_LAYER_QUERY_SESSION_RESULT.md) ·
[Read-only rigidity proof](READ_ONLY_THREE_LAYER_RIGIDITY.md) ·
[Parallel bilinear query R91](PARALLEL_BILINEAR_WORD_QUERY.md)

A phase-polynomial representation can change its work-parity rank after
adding an exact null identity. The following native construction makes
this issue concrete: any finite collection of affine-parity T phases
can be implemented in four physical T layers on private arbitrary dirty
work. The original parity inputs occur only as Clifford controls. This
also attains the retained four-layer scalar degree-three boundary.

The result does not implement higher-degree Boolean scalar phases,
improve the asymptotic R91 query law, or change the complete-frame upper
law R96. No external priority claim is made for the phase identity or
its synthesis application.

## 1. A literal null identity

For a bit p and an arbitrary four-bit word w, all sums in the exponent
below are ordinary integer sums of Boolean values. Then

```math
\sum_{v\in\mathbb F_2^4}(p\oplus v\cdot w)\equiv0\pmod8.
                                                               \qquad(1)
```

If w is nonzero, the sixteen summands contain exactly eight ones.
If w is zero, the sum is 16p. With
$`\omega=e^{i\pi/4}`$ and $`\sigma\in\{1,-1\}`$, (1) gives

```math
\prod_{v\ne0}\omega^{-\sigma(p\oplus v\cdot w)}
=\omega^{\sigma p}.                                    \qquad(2)
```

This is an exact pointwise phase equality, including its global phase.
No summation over a prepared work state or averaging over unknown work
is used. Every value of w gives the same phase on p.

## 2. Four actual native layers

Let $`p(x)=\kappa\oplus\bigoplus_{j\in S}x_j`$ be a fixed
affine parity, with $`h=|S|`$ and $`\kappa\in\{0,1\}`$.
The four reserved work bits $`w_0,w_1,w_2,w_3`$ are arbitrary.
Represent a binary row by the integer whose bit i is its coefficient
on $`w_i`$; for example, row 9 is $`w_0\oplus w_3`$.
The following rows partition all fifteen nonzero work rows into four
independent sets.

| Layer | Active work rows | Actual CNOT basis word, in chronological order | Work mask for p |
|---|---|---|---:|
| 1 | 1, 2, 4 on wires 0, 1, 2 | Empty word | 7 |
| 2 | 9, 3, 5, 8 on wires 0, 1, 2, 3 | $`0\to1,\ 0\to2,\ 3\to0`$ | 14 |
| 3 | 7, 13, 6, 15 on wires 0, 1, 2, 3 | $`1\to2,\ 2\to0,\ 0\to3,\ 3\to1`$ | 4 |
| 4 | 11, 12, 14, 10 on wires 0, 1, 2, 3 | $`1\to3,\ 3\to0,\ 3\to2,\ 2\to1`$ | 8 |

Each listed basis word acts on the four work bits only. In layer 1
the inactive fourth row is 8. In every layer each active row has binary
inner product one with the displayed mask. Consequently the following
chronological word puts $`p\oplus v\cdot w`$ on its designated
physical wire, for every active row v:

1. XOR p into each work wire selected by the mask. Use one CNOT from
   each address bit in S to each selected wire, and one X per selected
   wire if $`\kappa=1`$. Each X is the literal native word
   $`HSSH`$, containing four elementary Cliffords.
2. Execute the listed work-only CNOT basis word.
3. Apply $`T^{-\sigma}`$ to each active wire, simultaneously.
4. Reverse the actual basis word, then reverse the actual word from
   step 1. The actual inverse of each X expansion is
   $`HS^{-1}S^{-1}H=HSSH`$, using the exact Clifford identity
   $`S^{-2}=S^2`$; it also costs four elementary Cliffords.

The two reversals in step 4 are actual inverses. The work returns after every complete
layer word on every computational input, and its diagonal phase is
retained. The fifteen T phases are exactly the factors of (2). Thus
the full circuit is the literal operator

```math
\sum_x\omega^{\sigma p(x)}|x\rangle\langle x|
       \otimes I_w.                                    \qquad(3)
```

The allocation is exactly four additional arbitrary dirty qubits and
zero clean qubits. The basis-word lengths sum to $`0+3+4+4=11`$;
the mask weights sum to $`3+3+1+1=8`$. Charging both computations
and both actual inverses gives the simultaneous exact ledger

```math
T=15,\qquad D_T=4,\qquad
G=22+16h+64\kappa,\qquad b=4,\qquad\varepsilon=0.
                                                               \qquad(4)
```

Here G counts the model's elementary H, S and CNOT gates; X is charged
through its four-gate native expansion. There is no Clifford-depth
claim. The explicit formula includes a constant affine term when one
is present; constants or trivial parities may of course be simplified
separately. Neither initialization nor an uncharged copy of p is used.
Every original input wire is unchanged throughout and appears only as
a CNOT control. Every physical T gate acts on work. Exact equality
(3) proves zero leakage, complete return and validity with arbitrary
external references, including correlations among all reserved work.

## 3. Parallel phase words and a sharp scalar example

For m specified affine parities $`p_j`$ and signs $`\sigma_j`$,
give each factor a separate four-bit dirty bank and run the four
corresponding T layers in parallel. The address controls may be shared.
The bank computations have disjoint targets and leave every control
unchanged, so their elementary Clifford words can be serialized at
each layer boundary. The simultaneous resources are

```math
T=15m,\quad D_T=4,\quad b=4m,\quad
G=22m+16\sum_j h_j+64\sum_j\kappa_j,\quad\varepsilon=0.  \qquad(5)
```

All T gates within a physical layer have distinct targets. The result
is the specified product of literal phases on the same circuit, with
identity on the entire dirty register. The banks need not be
independent input states. Equation (5) is an upper construction; it
does not assert optimal count or width.

Apply (5) to the standard integer parity identity

```math
a+b+c-(a\oplus b)-(a\oplus c)-(b\oplus c)
 +(a\oplus b\oplus c)=4abc.                             \qquad(6)
```

There are seven signed parities, with total support size twelve and
no affine constants. Their seven private banks therefore implement

```math
(-1)^{abc}I_w,\qquad
T=105,\quad D_T=4,\quad G=346,\quad b=28,
\quad\varepsilon=0.                                    \qquad(7)
```

All three bits a,b,c remain read-only address wires. No phase
correction, supplied gate, clean helper or evaluator is omitted.
The [zero-clean three-layer scalar bound](READ_ONLY_THREE_LAYER_RIGIDITY.md#4-literal-scalar-phases-have-degree-at-most-two)
excludes Boolean scalar degree three below error $`\sqrt2/6`$.
Consequently four layers are optimal for this scalar CCZ query within
the stated read-only interface, including that approximation radius.
The scalar degree-three bound at four layers is therefore attainable.

This construction is not a cheaper AND query than the retained
[native four-layer AND](READ_ONLY_THREE_LAYER_RIGIDITY.md#7-native-layer-fixtures-and-the-clean-qubit-distinction):
the latter uses only eight T gates, 18 elementary Cliffords and one
additional dirty helper. The difference is that all three scalar CCZ
arguments in (7) remain original address wires. Reclassifying c as an
arbitrary target and adding two Hadamards gives an AND query, but pays
105 T gates, 348 Cliffords and 28 helpers, with no layer advantage.

For a batch of $`M=qH`$ R91 bilinear phase terms
$`(-1)^{y_i z_j s_{ij}(X)}`$, apply (7) with the third argument
replaced by its specified linear row form. If its support size is
$`\ell_{ij}`$, the resulting alternative has

```math
T=105M,\quad D_T=4,\quad b=28M,\quad
G=282M+64\sum_{i,j}\ell_{ij}.                            \qquad(8)
```

For the XOR shear add the two charged Hadamards per output bit.
The retained R91 phase word instead has eight T layers, at most 8M
T gates and M private dirty wires. Equation (8) trades larger count
and simultaneous width constants for fewer layers in that component;
both already have constant depth. This does not reduce the logarithmic
address-indicator cost, change R91's asymptotic envelope, or establish
the faster query premise in R98. The production query construction and
all complete-frame statements remain unchanged.

## 4. Why the rank candidate was reformulated

The exploration first produced a 32-term phase identity with seven
dirty parity rows and a five-layer independent-row packing. Its
global density $`32/7\gt4`$ rules out four layers for that list.
It does not rule out other native representations. The construction
above replaces that list and reaches four layers, at increased count
and width. The lower-bound parameter therefore has to survive exact
null phase identities, not just Clifford changes of coordinates on a
chosen list.

For example, the sum of all fifteen nonzero linear parities on four
dirty bits is zero modulo eight: it is eight on every nonzero dirty
input and zero at the origin. These are distinct parities, so adding
this null word increases the list length by fifteen and its row rank
by four without being merely a duplicate T/T-dagger cancellation.
Disjoint copies invalidate a universal proposed inequality
$`m\gt4\,\mathrm{rank}(W)`$ based only on total term count m
and dirty row rank W. In a fixed list the stronger condition is that
every subset fit the available number of independent sets; the
present relocation identity also changes the list itself. A useful
lower invariant must be stable under both operations.

The phase mechanism does not supply a cubic XOR query at five layers.
An affine parity has the integer expansion

```math
\bigoplus_{i\in S}u_i
=\sum_{\varnothing\ne A\subseteq S}
       (-2)^{|A|-1}\prod_{i\in A}u_i.
```

Modulo eight, only terms of degree at most three survive. Therefore
every CNOT/X and diagonal-T phase word has a weighted phase polynomial
of degree at most three in all original inputs, including arbitrary
work. If its literal full-space output is a Boolean scalar phase
$`(-1)^{g(x)}`$ independent of work, restricting work to zero and
using uniqueness of the multilinear expansion modulo eight forces
$`\deg g\le3`$. Hadamards only at an input/output boundary can
convert a scalar phase to an XOR derivative, giving degree at most
two. Intermediate changes of basis beyond this diagonal mechanism
are required to reach the first open five-layer cubic-XOR test.

## 5. Evidence and limits

The [bounded native checks](../../tests/test_read_only_phase_relocation.py)
verify the sixteen-term null identity, the actual CNOT basis words and
their inverses, every four-bit dirty input for affine parity phases,
the simultaneous ledgers, and the seven-factor scalar CCZ identity.
They use exact integer phases and basis permutations rather than an
unnecessary matrix on 31 qubits. The full-space proof is (1)–(3);
finite checks do not establish a new general lower bound.

The remaining constructive question is whether a genuinely
non-diagonal native five-layer mechanism can implement XOR degree
three, or whether a further invariant excludes it. Any larger-scale
escape must also meet the count, Clifford and simultaneous-width
requirements of the complete query interface. Optimal general depth
and unrestricted one-clean four-layer rigidity remain open.
