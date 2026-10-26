# Native sharing required by a complete approximate indicator

[Current indicator target](JOINT_INDICATOR_SESSION_RESULT.md#2-the-complete-constructive-target-and-paid-transfer) ·
[Selective inversion](SELECTIVE_INVERSION_DIRTY_INDICATOR.md) ·
[Constant-accuracy frame transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md)

A complete approximate indicator must use its non-Clifford gates across
many bank observables. The statement below permits arbitrary Clifford
interlayers, temporary address and bank motion, noncommuting bases and
arbitrary dirty work. It is a necessary incidence condition, not a
dimension-growing T-depth lower bound. Its independent-row corollary
does rule out linear T-count for separately emitted row words, even if
each row is synthesized jointly and no individual row-error estimate is
assumed. Cross-row synthesis remains outside that corollary.

Throughout, let $`p\ge2`$, $`m=2^p`$, and

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle.
```

The native source U uses no initialized work and satisfies the literal
full-operator estimate

```math
\|U-I_p\otimes I_W\|_{\mathrm{op}}\le\tau\lt\tfrac12.       \qquad(1)
```

W is arbitrary and may be correlated with the address, bank and an
external reference. Approximate return is included in (1); exact return
at intermediate boundaries is not required. The actual adjoint of U is
used in every conjugation below.

## 1. An active T gate adds at most one support direction

Ignore phases only when labelling Hermitian Paulis by vectors in
$`\mathbb F_2^{2Q}`$, where Q is the total physical width. For a
nonzero operator A, let its affine Pauli support dimension mean the
dimension of the affine hull of labels with nonzero coefficients in its
unique Pauli expansion. Coefficients and their signs remain part of A.

A Clifford conjugation acts on labels by an invertible binary linear
map and preserves this dimension. Conjugation by a physical T or
T-dagger on an axis z leaves each commuting term unchanged and replaces
an anticommuting term with a linear combination on labels a and a+z.
Consequently

```math
\mathrm{supp}(TAT^\dagger)
 \subseteq\mathrm{supp}(A)\cup(\mathrm{supp}(A)+z).             \qquad(2)
```

The affine dimension increases by at most one. Cancellation can only
reduce the support. If $`[A,Z]=0`$ on that physical wire, the conjugation
does nothing at all.

Start with any one Pauli and propagate it through the actual native
word in chronological order. Call a T occurrence active for this
observable when its immediately preceding observable fails to commute
with the physical Z on that T target. If there are q active occurrences,
its final affine support dimension is at most q. In particular it is at
most the word's entire T-count. This proof imposes no restriction on
the Clifford motion between those occurrences.

## 2. Each bank observable requires p active occurrences

For each row label u, propagate $`Z_{Y_u}`$ through U and write

```math
A_u=U Z_{Y_u}U^\dagger,
\qquad B_u=I_p Z_{Y_u}I_p^\dagger
             =Z_{Y_u}(-1)^{[x=u]}.
```

Unitary telescoping gives $`\|A_u-B_u\|\le2\tau`$.
Suppose the affine support of A_u had dimension less than p. Choose a
support label a and span all its support differences. Fewer than p
binary constraints then leave a nonzero address translation
$`X_x^v`$, $`v\in\mathbb F_2^p`$, whose commutation character is zero
on every difference. It has one common commutation sign on the entire
affine support. Thus, for some $`\sigma\in\{+1,-1\}`$,

```math
X_x^v A_u X_x^v=\sigma A_u.                                  \qquad(3)
```

The desired B_u has neither sign symmetry. Its sign differs from that
of its translate at x=u and x=u+v, so the plus comparison has norm two.
At every other address the two signs agree, so the minus comparison
also has norm two. Such another address exists because p is at least
two. Hence

```math
2=\|X_x^v B_u X_x^v-\sigma B_u\|
 \le2\|A_u-B_u\|\le4\tau,
```

contradicting (1). The affine dimension is therefore at least p, and
Section 1 proves

```math
q_u\ge p\qquad\hbox{for every bank row }u.                    \qquad(4)
```

This is a robust operator-norm argument. It does not infer a support
lower bound from small individual Fourier coefficients, which could
disappear under constant-error approximation. It also explains the
restriction p at least two: the one-bit indicator is Clifford, and
the missing third address makes the minus-sign argument fail.

The same proof applies to the phase target
$`D_p|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle`$ by propagating
$`X_{Y_u}`$ instead. Hadamards on the bank relate the two native
problems without changing T-count or T-depth.

## 3. The unrestricted incidence condition

For each actual T occurrence g, let $`\iota_g`$ be the number of the
m propagated bank observables for which that occurrence is active.
Double counting (4) yields

```math
\sum_g\iota_g=\sum_u q_u\ge mp.                              \qquad(5)
```

This condition holds for the unrestricted zero-clean native model in
(1), including a word that only restores its address at the end. If its
T-count is at most $`C_\tau m`$, the average incidence of its T gates
is at least $`p/C_\tau`$. Thus a linear-count source must share native
non-Clifford action across a growing number of these observables.

An additional architectural bound $`\iota_g\le k`$ would imply
$`T\ge mp/k`$. There is no such general constant bound here: Clifford
motion can make one physical T occurrence active for many observables,
and $`\iota_g`$ may equal m. Equation (5) alone only implies
$`T\ge p`$. It gives no growing native-depth lower bound, does not
exclude the retained linear-count indicator, and does not move R96's
complete-frame frontier.

Activity is defined from actual operator commutators, not from a gate's
wire name or a classical dependency graph. Temporary address motion can
make a gate active for observables belonging to other completed blocks.
Consequently the next corollary is proved by completed-block isolation,
not by assigning a false incidence cap to individual gates.

## 4. Separately emitted bank blocks pay once per block

Partition the m bank wires into nonempty sets $`S_1,\ldots,S_b`$.
Suppose an implementation of the phase target is emitted as complete
native words

```math
U=U_b\cdots U_1,
\qquad U_j=\sum_x|x\rangle\langle x|\otimes V_{j,x},           \qquad(6)
```

where V_j acts only on bank wires in S_j and a private arbitrary work
system W_j. The W_j are disjoint. The completed U_j preserves the
computational address; it may internally mix address and bank, and it
need not preserve its bank or return its private work at its boundary.
Let t_j be the T-count of that separately emitted native word. Assume
only the global bound $`\|U-D_p\otimes I\|\le\tau\lt 1/2`$.

The completed blocks commute because they are block diagonal in the
same address and act on disjoint systems within each address sector.
For u in S_j, every other block therefore cancels from

```math
U X_{Y_u}U^\dagger=U_j X_{Y_u}U_j^\dagger.                   \qquad(7)
```

The right side is within $`2\tau`$ of
$`X_{Y_u}(-1)^{[x=u]}`$. Apply Sections 1–2 to the native word U_j
itself to obtain $`t_j\ge p`$. Thus the emitted ledger satisfies

```math
T_{\mathrm{emitted}}=\sum_{j=1}^b t_j\ge bp.
                                                               \qquad(8)
```

If every block has at most k rows, (8) is at least
$`p\lceil m/k\rceil`$. Independent row emission is the case k=1:
it needs at least mp T gates at this constant full-bank accuracy.
There is no premise that row errors are independently small, so hiding
a row-dependent phase in another row cannot evade (7).

This excludes the target $`O_\tau(m)`$ count whenever the maximum block
size is $`o(p)`$, including every fixed-size completed-block architecture.
It leaves a common joint word, cross-block resynthesis, overlapping
private work, and block allocations with $`bp=O(m)`$ available. In
particular, a uniform size cap must be $`\Omega(p)`$ to avoid this count
obstruction. A fixed nontrivial address/bank
conjugator around (6) is not automatically covered: its probe observable
must be analysed in that changed representation. The unrestricted
incidence statement (5) remains available for its actual native word.

## 5. Research consequence

The selective-inversion study improves the analytic error of a dirty
row pulse and then charges its finite native realization. Equation (8)
shows that even an ideal improvement of each separately emitted row
compiler cannot give the requested linear T-count. This is a reason to
stop that architecture, not a reason to stop the complete indicator
problem.

The next constructive question is a word sharing non-Clifford operations
across a growing set of bank observables, with both full-bank error and
the native shared operations proved. Noncommuting address/bank motion
across those boundaries is allowed. The existing exact indicator shows
that sharing can recover linear count; its remaining logarithmic
address depth is the cost to improve. The current R98 transfer still
requires all counts, widths, accuracy and depth on the same source.

The later [root-fusion reduction](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
places that sharing in an outer dirty tree: a smaller complete source
may have polynomial address overhead, which the sparse first stage
absorbs. Thus (8) continues to constrain standalone emitted banks but
does not exclude polynomial-cost independent rows inside that smaller
source. Its joint full-bank error and sublogarithmic native depth are
still unpaid. No unrestricted depth claim is added by this relaxation.

The [focused fixtures](../../tests/test_native_indicator_sharing.py)
check native support propagation and the signed-translation distinction.
The all-width incidence and block bounds are analytic statements above;
finite fixtures do not prove an asymptotic depth obstruction.

<a id="complete-probe-support-count"></a>

## 6. The complete probe requires a linear number of source T gates

A stronger count condition applies when one Pauli image is the entire lookup
phase. Let the physical width be $`n=p+m+w`$, with w arbitrary work wires,
and let P be a Hermitian Pauli. Suppose an actual native source V satisfies

```math
A=VPV^\dagger,\qquad
\|A-Q_p\|\le\epsilon,\qquad
Q_p|x,Y,W\rangle=(-1)^{Y_x}|x,Y,W\rangle.
```

There is no initialized input, and V may move every register. Return is
required by the complete target observable, not separately of V. Then

```math
\begin{aligned}
\epsilon\lt1&\quad\Longrightarrow\quad
 \text{active T occurrences for P}\ge m+p-1,\\
\epsilon\lt2\sin(\pi/8)&\quad\Longrightarrow\quad
 T(V)\ge m+p,\\
\epsilon=0&\quad\Longrightarrow\quad T(V)\ge m+p+1.       \qquad(9)
\end{aligned}
```

These are source T-count bounds. They are distinct from the
[complete indicator's unitary-nullity bound](INDICATOR_PAULI_NULLITY.md),
which applies directly to a query word. Neither statement gives growing
native depth with the permitted exponential physical width.

**Signed Pauli symmetries give a robust support inclusion.** For a nonzero
operator B, let L(B) span the differences between labels in its Pauli support.
A Hermitian Pauli R obeys $`RBR=\sigma B`$ for some
$`\sigma\in\{+1,-1\}`$ precisely when its commutation character is constant
on that support. Its label therefore belongs to $`L(B)^\perp`$.

For $`Q_p`$, write R's X mask as (a,b,d) on the address, bank and work.
Conjugation replaces its sign by $`(-1)^{Y_{x+a}+b_{x+a}}`$.
An arbitrary bank rules out every nonzero a. When a is zero, a common sign
is possible exactly for the two constant bank masks b; d and every Z mask
are free. Thus the signed-symmetry label space has dimension
$`n+w+1`$. Every forbidden signed comparison has norm exactly two, because
it compares two different diagonal sign functions. In particular,

```math
RAR=\sigma A,\quad\epsilon\lt1
\quad\Longrightarrow\quad RQ_pR=\sigma Q_p,
\qquad L(Q_p)\subseteq L(A).                              \qquad(10)
```

Otherwise telescoping would give $`2\le2\epsilon`$. The target's exact
Pauli expansion also identifies the difference space explicitly:

```math
\begin{aligned}
Q_p&=\frac1m\sum_{u,t\in\mathbb F_2^p}
 (-1)^{u\cdot t}Z_x^t Z_{Y_u},\\
L(Q_p)&=\{\text{address Z masks and even-parity bank Z masks}\},\\
\dim L(Q_p)&=p+m-1.                                    \qquad(11)
\end{aligned}
```

Its work masks vanish, and this space is isotropic. Equations (10)--(11)
and Section 1 prove the first line of (9), including arbitrary work and
references. Small individual Fourier coefficients are not used as a
continuity argument.

**Saturating that support count would make the source third-level.** Collect
the Clifford gates of the actual word into one final Clifford C, leaving
native T-angle Pauli rotations about fixed labels $`a_1,\ldots,a_t`$.
Their literal scalar factors are retained. If $`\mathcal W`$ is their binary
span, Pauli conjugation gives
$`L(A)\subseteq C\mathcal W`$, where C denotes its symplectic action on
labels. If $`t\le m+p-1`$, (10)--(11) force equality throughout:

```math
L(Q_p)=L(A)=C\mathcal W,\qquad \dim\mathcal W=m+p-1.
```

Hence every collected axis commutes with every other. A simultaneous
Clifford diagonalization makes their actual product a weighted cubic phase,
so V belongs to the third Clifford-hierarchy level and A is a Hermitian
Clifford. For p at least two, $`Y_x`$ has Boolean degree p+1 and $`Q_p`$
is not Clifford. The
[diagonal-sign threshold](CLIFFORD_INTERLEAVED_LOOKUP.md#2-two-elementary-clifford-thresholds)
therefore gives $`\|A-Q_p\|\ge2\sin(\pi/8)`$, proving the second line
of (9). The strict threshold and the distinction between total and active
T counts are retained.

For the exact last line of (9), put $`r=m+p-1`$ and suppose t=r+1.
If the collected axis span has dimension r, the same isotropic argument
applies. Otherwise all r+1 axes are independent. At each active rotation,
the previous support and its translate by the new axis are disjoint;
every old coefficient survives and the affine support dimension increases
by exactly one. An inactive rotation does nothing. Exact final dimension
r therefore leaves one inactive occurrence, which can be omitted without
changing the probe image. The resulting r-T source contradicts the bound
already proved. This argument uses exact support, not a continuity claim
for nonzero approximation error.

Support dimension itself can grow by r in a single physical band:
conjugating $`X_1\cdots X_r`$ by parallel T gates gives
$`\bigotimes_{j=1}^r(X_j+Y_j)/\sqrt2`$, with affine support dimension r.
A fixed Clifford can supply that input probe from any specified nonidentity
Pauli. Consequently (9) does not turn into a logarithmic-depth obstruction
when $`n\ge m+p`$; one must control more than support growth per band.
The [complete-probe fixtures](../../tests/test_complete_probe_support_count.py)
check the signed-symmetry gap, target support, independent-axis growth,
literal native propagation and this one-band limitation. They do not
establish an asymptotic depth bound.
