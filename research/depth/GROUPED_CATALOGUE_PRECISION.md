# Growing-precision depth from grouped dirty catalogues

[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[All-precision baseline](UNIFORM_SUBTREE_DEPTH.md) ·
[Dirty catalogue identity](DIRTY_CATALOGUE_TERMINAL.md) ·
[Executed design](PRECISION_TERMINAL_PROMPT.md)

**Local result R93.** Small complete groups can share a long-prefix
query even when their constituent logical layers do not commute. A
dirty catalogue executes the group's coarse operator. One shared binary
source corrects the residual of each entire group, with an explicit
matrix-unit dictionary and retained failure history. This improves the
inverse-polynomial-in-n precision benchmark while preserving the
one-clean contract and the original count and width bounds.

## 1. The precision range and simultaneous theorem

For integers $`n\ge1`$ and L put

```math
N=2^n,\qquad \ell=\log_2(n+2),\qquad
\eta=2^{-L},\qquad 6\le L\le2\ell,
```

```math
F_{93}(n,L)=\ell\log_2(\ell+2)
             +\frac{L\ell}{\log_2(\ell+2)}.
```

For every prescribed complete real Hopf frame W and every integer
$`b\ge17(L+n+7)`$, one external initialized qubit and at most b
arbitrary dirty qubits suffice for one coherent Clifford+T circuit V
with

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+F_{93}(n,L)\right).
\qquad(1)
```

The prescribed complex magnitude completion has the same bounds at
$`b\ge17(L+n+8)`$. Every constant in (1) is absolute, independent
of n, L, b and the prescribed angle and phase data on this precision
range. The error includes every logical column, literal phases,
initialized-work leakage and return, arbitrary dirty inputs and their
external reference correlations. Unused supplied wires remain untouched.
G counts elementary Clifford gates, including those between T layers.
There is no QRAM, measurement, reset, supplied precision state, extra
initialized bank, or free coherent evaluator. The terminating certified
evaluation interface is unchanged; classical preprocessing has no
efficiency claim.

At $`L=\Theta(\ell)`$ and $`b=\Theta(\sqrt{NL})`$, (1) gives

```math
T=O(\sqrt{NL}),\qquad G=O(NL),\qquad
D_T=O\!\left(\frac{\ell^2}{\log_2(\ell+2)}\right).
\qquad(2)
```

The accuracy here is inverse polynomial in n, the number of logical
qubits, rather than in N. This strictly improves R79's
$`O(\ell^2)`$ upper bound at that benchmark. It does not establish
the design's stronger $`O(\ell\log\ell)`$ target or a matching
lower bound. R79 remains available at every precision; R90 and the
other retained constructions may be stronger on other ranges.

We first prove the wide real construction uniformly on the slightly
larger internal range $`6\le L\le3\ell`$. This allows the
one-bit precision shift in the complex composition. Every cutoff below
is absolute after fixing the native constants. No fixed-L asymptotic
statement is used with L silently varying.

## 2. One uniform coarse reduction

Take the wide branch $`b\ge\sqrt N/n`$ and n beyond an absolute
cutoff to be fixed after the native constants. First make the retained
finite whole-family magnitude/sign reduction, at error at most
$`\eta/512`$, and pay its two exact Boolean sign phases once.
Their rows use the full n-bit logical address; the sharp one-bit query
bounds give

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\ell).
```

The finite sign data include originally singular and zero-amplitude
branches. The reduction needs no exact sign or zero decision on an
analytic input. It is the same
[full-column sign gauge](PREFIX_SELECTED_COIN_TERMINAL.md#2-exact-factorization-and-where-rounding-occurs)
used by R89 and R90. All subsequent factors refer to the resulting
positive frame until the two sign phases are restored.

Use exactly one
[R89 local causal stage](CAUSAL_COMPACT_COIN_DEPTH.md#102-one-stage-addresses-accuracy-and-resource-bounds),
starting with m equal to n and with

```math
l_0=L+\lceil\log_2(n+2)\rceil+10,\qquad
K_0=n+l_0+2,\qquad
u_0=l_0+\lceil\log_2K_0\rceil,
```

```math
m_1=\lceil A u_0^3\rceil.
\qquad(3)
```

The native constants C,d and the reserve constant A are fixed in R89
before this choice. In particular,

```math
m_1-Cu_0^3\ge(2d+80)u_0.
```

Uniformly for $`L\le3\ell`$, we have $`u_0=O(\ell)`$ and
$`m_1=O(\ell^3)`$. Hence $`m_1\le n/2`$ beyond one absolute
cutoff. The stage is therefore well defined without invoking R89's
fixed-L stopping height. Its absolute native ledger yields

```math
T_0=O(\sqrt N+N/b),\qquad G_0=O(N),\qquad
D_0=O(N/b^2+\ell).
\qquad(4)
```

For example, the large tables have at most
$`2^{n-m_1+Cu_0^3}`$ rows times fixed powers of $`K_0`$;
the reserve margin suppresses these factors by $`K_0^{-40}`$ or
more in the square-root count. Standalone work is a fixed polynomial
in n and fits $`b\ge\sqrt N/n`$ beyond the same absolute cutoff.
The sharp-query private-pool condition is polynomial in n here and
also fits. Its error is at most

```math
2^{-l_0}/4\le\frac{\eta}{4096(n+2)}\le\frac{\eta}{2048}.
\qquad(5)
```

Finish the stage's actual saved-predicate inverse before proceeding.
The remaining bottom $`m_1`$ layers retain their prescribed values.
Their upper $`n-m_1`$ logical bits form an immutable external prefix
within the next local operation. Approximate flag return is handled by
initialized-isometry telescoping, never by a reset.

## 3. A prefix-selected early decoder with local overhead

The next lemma reduces the remaining height from $`O(\ell^3)`$
to $`O(L+\log(\ell+2))`$. We give its native ledger because merely
substituting the local height for n in R79 would erase the cost of the
long prefix.

Let m layers remain, with immutable prefix P of length $`n-m`$.
Choose a protected suffix Z of M bits and let $`h=m-M`$. The exact
early factor is

```math
E=\sum_P|P\rangle\langle P|\otimes
 \left(W_{P,h}\otimes\Pi_Z+I\otimes(I-\Pi_Z)\right),
\qquad \Pi_Z=|0^M\rangle\langle0^M|.
\qquad(6)
```

Use the original
[R78 complete tagged source and subtree search](SUBTREE_PROJECTOR_DECODER.md#2-the-literal-nested-subtree-spaces),
now direct summed over P. Its sources are the original descendant
loaders with fresh local gradient preparation and its actual inverse on
every occurrence. No shared-gradient invariant is assumed.

For integer local precision l define

```math
K=m+l+2,\quad s=\lceil\log_2(h+1)\rceil,\quad
\delta=\delta_0=\frac{2^{-l}}{64(s+1)},\quad
p=l+7+\lceil\log_2(s+1)\rceil.
```

Save $`H=[Z=0^M]`$ in the sole external initialized flag. The three
s-bit tag, search and threshold records, one reused $`O(p)`$ source
bank, and fixed extra fields all occupy Z on H one. P is never copied
to initialized work, modified, routed or tested for zero. The native
conditional bank has width $`O(l+\log(m+2))`$.

For each P the nested projectors and literal comparison reflections are
exactly those of R78. A comparison reflection acts on the transported
all-zero-source-work space; it is not a logical projector tensored with
identity on arbitrary source work. Its full-norm error is at most
$`4\delta`$. Taking a direct sum over P takes a maximum, not a sum.
The complete early error therefore obeys

```math
\|V_EJ_1-J_1(E\otimes I_b)\|
\le\delta_0+8s\delta\le2^{-l}/8.
\qquad(7)
```

The physical completed word is literally identity on H zero, including
occupied source and search inputs. The final actual predicate inverse
restores H on the ideal comparison output and propagates actual leakage
unitarily. All dirty query helpers return exactly.

### 3.1 Keep the table addresses global

Put $`Q_0=2^{n-M}`$. The retained R78 source grammar gives the
following safe native ledger after adding P only to tables that depend
on its angle data.

| Operation | Safe row bound | Output or standalone local cost |
|---|---:|---:|
| Descendant Boolean phases and calibration programs | $`C(h+1)Q_0p`$ | At most $`Cp`$ outputs |
| Full-source anchor phase | $`C(h+1)Q_0p\delta^{-2}`$ | One output |
| Variable-suffix Hadamard masks | $`C(h+1)`$ | h outputs |
| Packed tag queries | $`2^h`$ and $`2^s`$ | s and h outputs |
| Local gradients and coefficient preparation | No large table | Polynomial in p, depth $`O(p)`$ |
| Thresholds and guarded tag XOR | No large table | Polynomial in s |
| Saved-suffix boundaries | No P dependence | Polynomial in M, logarithmic depth |

The anchor and packed-tag tables do not actually need P; the larger
common envelope is harmless. P occurs in the source grammar only in
Boolean phase and calibration-program queries. The variable-suffix
Hadamards, local gradients, modular adders, anchor predicates, tag
packing, threshold arithmetic and saved-Z predicates act on local
registers and do not read P.

Replace every Boolean query by the
[R88 read-only sharp interface](PREFIX_SELECTED_COIN_TERMINAL.md#4-parallel-sharp-queries-with-a-read-only-public-address).
Its shared public address enters only through Clifford controls at
non-Clifford boundaries. Parallel one-bit outputs have disjoint private
output and helper pools. Thus output bits within a table position can
run in parallel; no logical sources or spectral-search stages acting on
the same target are parallelized. This is a property of the physical
native words and their actual inverses, not an inference from eventual
preservation of P.

Since $`\delta^{-2}=4096(s+1)^2 4^l`$, fix absolute polynomial
exponents bounding the local outputs, simultaneous instances, source
appearances and standalone work. The exponential precision cost
$`4^l`$ is a table-row factor, not a polynomial instance count.
Every address has at most

```math
n-M+2l+C\log_2K+C
```

bits. If I is the chosen polynomial bound on simultaneous instances,
assign each a private fraction $`b/(C'I)`$ of the dirty pool. For a
fixed absolute degree $`d_w`$, the sufficient condition

```math
b\ge C K^{d_w}(n+l+\log_2K+2)
\qquad(8)
```

fits every sharp-query minimum pool and all separate local outputs.
This is a sufficient wide-branch condition, not a claim at the literal
minimum allocation. Protected Z is not counted as dirty query work.

A fixed absolute degree d, chosen after all those native exponents,
gives the complete envelope

```math
\begin{aligned}
T_E&\le C K^d\left(2^l\sqrt{Q_0}+4^lQ_0/b+1\right),\\
G_E&\le C K^d\left(4^lQ_0+1\right),\\
D_E&\le C\left[(s+1)
 \bigl(\log_2(n+l+\log_2K+2)+l+\log_2K\bigr)
 +K^d4^lQ_0/b^2\right].
\end{aligned}
\qquad(9)
```

There are $`O(s+1)`$ sequential source/query positions; the number
of parallel output bits is already included in the polynomial count and
pool factors. No polynomial in global n multiplies T or G in (9).
Global n remains in the table volume, logarithmic address depth and
minimum-width condition (8).

### 3.2 Reduce the residual height

Now use $`m=m_1=O(\ell^3)`$ from (3) and $`l=L+4`$. Fix $`A_E`$
large enough after (8)–(9) and the conditional clean-bank constant, and
set

```math
M=\lceil A_E[l+\log_2K]\rceil.
\qquad(10)
```

If $`M\le m/2`$, apply the decoder. Choosing $`A_E\ge2d+40`$
ensures

```math
K^d2^l2^{-M/2}\le K^{-20},\qquad
K^d4^l2^{-M}\le K^{-20}.
```

Condition (8) holds uniformly for $`b\ge\sqrt N/n`$ beyond an
absolute cutoff, since K is polynomial in ell. The standalone
$`K^d`$ costs also fit. Hence

```math
T_E=O(\sqrt N+N/b),\qquad G_E=O(N),\qquad
D_E=O\!\left(N/b^2+\ell\log_2(\ell+2)\right).
\qquad(11)
```

If $`M\gt m/2`$, skip this decoder. The failed inequality itself,
together with $`m=O(\ell^3)`$, implies
$`m=O(L+\log_2(\ell+2))`$. Thus in either case the residual
height, denoted again by M, satisfies

```math
M\le C[L+\log_2(\ell+2)].
\qquad(12)
```

When used, (7) contributes at most $`\eta/128`$ error. A skipped
decoder contributes zero. Its saved flag has completed its actual
inverse before the grouped terminal starts.

## 4. Groups and a finite coarse catalogue

Choose

```math
g=\max\left\{1,\left\lfloor\frac{\log_2(\ell+2)}8
                                      \right\rfloor\right\},
\qquad
R_0=\min\{M,\lceil A_0\log_2(\ell+2)\rceil\}.
\qquad(13)
```

Leave the deepest $`R_0`$ layers to the retained one-clean bottom
compiler. Partition the other $`M-R_0`$ consecutive layers into
R groups of at most g bits, with full groups except possibly the first,
shallowest group. If the middle is nonempty,

```math
R=\left\lceil\frac{M-R_0}{g}\right\rceil
=O\!\left(1+\frac{L}{\log_2(\ell+2)}\right).
\qquad(14)
```

For a group let X have h target bits, $`1\le h\le g`$, and let
P be the whole logical prefix above X. Let G test only logical bits
below X and above the protected bottom bank Z of $`R_0`$ bits.
It excludes X and Z. Save $`H=[Z=0^{R_0}]`$ once in the sole
external clean flag. Each local G is computed with direct H control
into separate exact zero scratch, used for that group, and uncomputed
after its complete coarse word. It is not part of retained failures.

The ideal active group $`W_P`$ is the complete h-bit real Hopf
block of the positive tuple produced in Section 2 and has determinant
one. Its angle code therefore uses the same $`[0,\pi/2]`$ domain
as R90. Its constituent layers need not commute.
P stays unchanged within this whole group, but may contain logical
outputs changed by earlier groups; every group freshly queries its own
current P. No initial value of a changed logical target is reused.

For R positive set

```math
q_\pi=\lceil\log_2(8R)\rceil,\qquad
\pi=2^{-q_\pi},\qquad \gamma=(1-\pi)^R,
```

```math
J=2^{\lceil\log_2(L+20)\rceil},\qquad
D=2^h,\qquad K_h=4D^2,
```

```math
\varepsilon_h=\frac{\pi}{(1-\pi)K_hJ},\qquad
\delta_h=\varepsilon_h/4.
\qquad(15)
```

Then $`R\pi\le1/8`$ and $`7/8\le\gamma\le1`$.
Round the group's real angles to a finite code whose whole-operator
error is below $`\delta_h/2`$. The
[R90 depth-weighted angle code](DIRTY_CATALOGUE_TERMINAL.md#8-a-finite-full-frame-program-with-linear-in-dimension-length)
uses

```math
q_h=O\!\left(2^h[h+\log_2(R+1)+\log_2(J+1)]\right)
\qquad(16)
```

bits. Its angle precisions grow with the remaining local depth, so this
bound charges the whole group's tree, including nodes reached after
other target bits have changed. Each string defines an entire rounded
$`\mathrm{SO}(D)`$ frame.

For every code v, choose the actual native R90 commutator factors
$`\widehat A_v,\widehat B_v`$ so that their commutator $`C_v`$
approximates the coded frame within $`\delta_h/2`$. The finite
certified enumeration and the ancillary-free fixed synthesis bound give

```math
S_h=O\!\left(8^h[h+\log_2(R+1)+\log_2(J+1)]\right)
\qquad(17)
```

for each factor's T-count, Clifford count and T-depth. Actual inverse
words are retained. Scalar factor phases cancel in the commutator; the
coarse group has its literal calibrated phase. The all-zero angle code
is the identity. For that entry choose the actual factors
$`\widehat A_0=\widehat B_0=I`$ exactly.

Let the queried code $`f(H,P,G)`$ be the rounded active group when
H and G are one, and zero otherwise. R90's complete selected
commutator, using two disjoint arbitrary arrays and an arbitrary dirty
program register, is exactly

```math
C_{f(H,P,G)}(X)\otimes I_{\rm arrays,program,helpers}.
\qquad(18)
```

The address excludes all of X. Equation (18) holds on every input,
including arbitrary occupied inactive work and external references.
On H zero or G zero it is literally identity because the selected
actual factors are identities. All other array factors cancel as actual
inverses. There is no additional factor from the number of catalogue
entries in its approximation error.

The exact coarse circuit $`C_i`$ for a group is the physical circuit
(18), not its rounded ideal. On the active group
$`\|C_i-W_i\|\lt\delta_h`$. Although $`W_i`$ is real,
$`C_i`$ can be complex. The correction below is defined from that
actual $`C_i`$, including its literal phases.

## 5. A residual kernel for the whole changing target

For the group target X use the fixed dictionary

```math
C_i^\dagger W_i-I
=\sum_{\lambda=(u,v,\omega)}
 a_{i,\lambda}(P,G)\,\omega|u\rangle\langle v|_X,
\qquad K_h=4D^2,
\qquad(19)
```

where $`u,v\in\{0,1\}^h`$ and
$`\omega\in\{1,-1,i,-i\}`$. Split the real and imaginary parts
of each entry into nonnegative positive and negative parts. Every
coefficient lies in $`[0,\delta_h]`$ and is literally zero on
G zero. The coefficients are diagonal in preserved P and G and commute
with the target matrix units. No endpoint address of length n is
introduced; the two endpoint labels have only h bits each.

For each finite table row use certified nonnegative rational lower
enclosures and downward grid rounding, as in the retained binary
terminal. This gives digits with

```math
\widehat a_{i,\lambda}
=\varepsilon_h\sum_{k=0}^{J-1}f_{i,\lambda,k}(P,G)2^{-k},
\qquad
0\le a_{i,\lambda}-\widehat a_{i,\lambda}
\lt2\varepsilon_h2^{1-J}.
\qquad(20)
```

Inactive rows have literal zero digits. Positive/negative-part clipping
needs no exact sign or zero decision. The finite actual coarse word and
the prescribed angle interface allow certified entry evaluation. No
coarse approximation error is added separately later: the correction
target in (19) is exactly $`C_i^\dagger W_i`$.

Use the same
[paid shared binary source](SHARED_BINARY_TERMINAL.md#3-paid-geometric-source-and-the-complete-local-kernel)
on a source register of dimension $`M_s=2J`$:

```math
g_s=\sum_{j=0}^{M_s-2}2^{-(j+1)/2}|j\rangle
       +2^{-(M_s-1)/2}|M_s-1\rangle.
```

Its exact native preparation $`P_s`$ has $`O(J)`$ count and
conservative T-depth, $`O(\log J)`$ initialized width, separate
returned preparation scratch, and a charged actual inverse. For the
downward shifts $`S_d`$,

```math
\|S_dg_s-2^{-d/2}g_s\|\le e_s\lt2\,2^{-J},
\qquad 0\le d\lt2J.
\qquad(21)
```

All used displacements are $`d=2k`$. The shift dilation retains
its underflow flag. The source is prepared once for the entire middle
stream, not once per group.

Let F contain the rare bits and decoded rare flag, term label, digit
label, loaded digit, and distinct underflow, atom-rejection and
digit-zero flags. It excludes the source, history, saved H, global
padding fields, exact preparation scratch and G scratch. A group word
$`B_i`$ prepares uniform rare, term and digit labels and XOR-loads
the digit from

```math
(H,P,G,\lambda,k).
```

It excludes X and the source. Uniform preparation and rare decoding
have direct H control, and the table has zero H-zero rows. Consequently
$`B_i`$ is literally identity on the entire H-zero space.

SELECT is identity in common mode. In rare mode, digit zero flips the
separate digit-zero failure flag. In rare mode with digit one it applies
the downward source-shift dilation, the target matrix-unit dilation and
the literal phase $`\omega`$. The target dilation first XORs
$`[X\ne v]`$ into its retained atom flag and then replaces X by
$`X\mathbin\oplus u\mathbin\oplus v`$. Its zero atom-flag
block is exactly $`|u\rangle\langle v|`$.

To implement the comparison, conjugate the zero test by the controlled
bitwise XOR of v into X, toggle the mismatch flag, and reverse the
comparison conjugation before the final XOR. Use the retained returned
conjunctions with private arbitrary helpers; no comparison output is
silently initialized. The retained atom flag is initialized only on the
accepted input column, as required by a dilation. The native action is
unitary on every occupied flag input. A safe fixed polynomial in h bounds
count, depth and returned helpers. Since $`h=O(\log\ell)`$,
this cost is smaller than the per-group $`O(\ell)`$ allowance.
All consumed comparison, XOR, source, rejection and phase operations
have direct H and bounded rare/digit controls. The retained
[literal bounded-control native implementations](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#exact-conditioning-without-a-clifford-count-penalty)
and their residual phase corrections are used.

The source and atom failure flags are independent. Their simultaneous
zero block therefore gives the product of the two contractions; failed
branches are retained. SELECT preserves P, G, all labels and the loaded
digit. Thus the actual $`B_i^\dagger`$ unloads that digit even
after X changes, on both accepted and failed branches. All query helpers
return exactly. On H zero the whole word is identity on arbitrary
occupied inputs.

The F-zero block of $`K_i=B_i^\dagger\mathrm{SELECT}_iB_i`$ is

```math
T_i=(1-\pi)I+\frac{\pi}{K_hJ}
 \sum_{\lambda,k}S_{2k}\otimes
 f_{i,\lambda,k}(P,G)\,\omega|u\rangle\langle v|_X.
\qquad(22)
```

It is a contraction. On G zero, all rare terms reject and this block is
$`(1-\pi)I`$. The digit-zero rejection must not be bypassed by G:
its common normalization is needed on every logical input. H zero has
the separate exact occupied-space identity described above.

Writing $`E_g=|g_s\rangle\otimes I`$ and
$`D_i=(1-\pi)C_i^\dagger W_i`$, (20)–(22) imply

```math
\|T_iE_g-E_gD_i\|
\le\pi e_s+\frac{2\pi}{J}2^{1-J}.
\qquad(23)
```

There are $`K_hJ`$ norm-at-most-one source terms weighted by
$`\pi/(K_hJ)`$. The sum of coefficient rounding errors is
$`(1-\pi)K_h`$ times the bound in (20), which equals the second
term of (23). The growing dictionary dimension has therefore been
charged and canceled explicitly. It is not omitted from the error sum.
After $`K_i`$ execute the actual coarse $`C_i`$ in (18), then
return G scratch. That word acts as $`C_i`$ tensored with identity
on F, source and history, with all catalogue work returned exactly.

## 6. One source, retained failures, and full return

Between consecutive groups increment a separate history iff F is
nonzero, with direct original H control. There are R minus one possible
increments. For $`R\ge2`$ use $`\lceil\log_2R\rceil`$ bits;
for R one omit history. A temporary for the nonzero test and the
arithmetic scratch are outside F and history and return before the next
group. The counter never wraps on a path beginning at zero.

Exactly as in the
[retained nonwrapping-history proof](SHARED_BINARY_TERMINAL.md#4-local-and-stream-accepted-column-estimates),
the final F-zero/history-zero block inserts F-zero projectors at every
intermediate group boundary. A path with a recorded nonzero F has a
nonzero final counter; all others have intermediate F zero. Coarse
operators commute with these projectors. The source is never projected
between groups. The complete accepted block is consequently

```math
V_{\rm acc}=(I\otimes C_R)T_R\cdots(I\otimes C_1)T_1.
\qquad(24)
```

Contractive telescoping of (23) gives

```math
\|V_{\rm acc}E_g-E_g\gamma W_{\rm mid}\|\le\xi,
```

```math
\xi\le R\pi\left(e_s+\frac2J2^{1-J}\right)
\le\frac38\,2^{-J}\le\eta/1024.
\qquad(25)
```

Thus the cap $`J=2^{\lceil\log_2(L+20)\rceil}`$ suffices
independently of the group count and dictionary dimension. Their effects
are already paid through $`\pi`$, $`\varepsilon_h`$, the
coarse precision and the number of label bits.

Apply the actual $`P_s^\dagger`$ only after the whole stream.
Use a separate paid padding source targeting
$`w_0=1/(2\gamma)\in[1/2,4/7]`$, with error
$`\epsilon_0\le\eta/2048`$ in its probability. Its preparation
and inverse cost $`O(L)`$ and use $`O(\log L)`$ initialized
bits. The rejection coordinate is a distinct bit inside Z, as in the
[R72 one-clean relocation](ONE_CLEAN_DEPTH.md#4-relocate-the-terminal-padding-flag).
The chronological half-block is

1. Prepare padding unconditionally.
2. Prepare the main source unconditionally.
3. Execute the H-guarded group stream and history.
4. Apply the actual main-source inverse.
5. Flip padding rejection iff H is one and padding mode is zero.
6. Apply the actual padding-preparation inverse.

On H zero, the stream, history and rejection are literal identities
on the full occupied-work space. The unconditional preparations cancel
in reverse nested order, with no initialization premise. On H one the
initialized block B satisfies

```math
\|2B-W_{\rm mid}\|\le2\xi+2\epsilon_0.
```

Let $`J_0`$ initialize all protected Z coordinates, including padding
rejection, unused coordinates and returned G scratch; it excludes saved
H and arbitrary dirty inputs. The physical reflection $`R_H`$ is
identity on H zero and $`I-2J_0J_0^\dagger`$ on H one. Source
inversion has already occurred at this boundary. The actual recovery

```math
\mathcal A_H=Z_HU_HR_HU_H^\dagger R_HU_H
\qquad(26)
```

has full initialized-isometry error

```math
8(\xi+\epsilon_0)\le3\eta/256\lt\eta/16.
\qquad(27)
```

The sign is the physical $`Z_H`$, not an unconditional minus sign.
Every inverse is the actual inverse of the entire finite word,
including sources, coarse catalogues, failed branches and history.
Reverse the original H predicate last. The ideal middle target preserves
Z zero and returns its work; the inverse propagates actual leakage by
unitarity. The equality and bound tensor with arbitrary external
references.

The simultaneously initialized fields have total width

```math
O\!\left(\log J+\log(R+1)+g\right)=O(\log_2(\ell+2)).
\qquad(28)
```

This includes the two sources and separate exact preparation scratch,
rare and digit fields, the two h-bit endpoint labels, retained failures,
history, padding rejection and separate exact local predicate/arithmetic
scratch. Choose $`A_0`$ after this native constant. The coarse
program and both catalogue arrays are arbitrary dirty registers outside
Z and do not consume conditional initialized work. A common maximal F
may be used for a shorter group: only its active label coordinates are
prepared uniformly; surplus coordinates remain untouched and zero on
the active initialized column. Their presence in F and the final
reflection does not alter $`K_h`$ or the accepted contraction.

## 7. Charge the tables, catalogues and simultaneous work

A group with h target bits and s lower logical bits has prefix length
$`n-s-h`$. The digit table has constant output width and at most

```math
Q_{d,i}\le C2^{n-s-h}K_hJ
         \le CNJ2^{h-s}
\qquad(29)
```

rows. Its address includes H and G; their factor is included in C.
The whole coarse-code query has bit volume

```math
Q_{c,i}q_i\le CN2^{-s}
  [g+\log_2(R+1)+\log_2(J+1)].
\qquad(30)
```

These addresses exclude the whole group target and every unrelated
live bank. They have at most $`n+O(\log\ell)`$ bits, so their
sharp-query logarithmic address depth is $`O(\ell)`$. There are
a constant number of completed queries per group, including inverses
and the constant recovery multiplicity.

Let $`\Lambda=g+\log_2(R+1)+\log_2(J+1)`$.
Consecutive full groups have lower-suffix lengths differing by g, and
the possible short group is furthest from the bottom. Therefore

```math
\begin{aligned}
\sum_iQ_{d,i}&\le CNJ2^{g-R_0},&
\sum_i\sqrt{Q_{d,i}}&\le C\sqrt{NJ2^{g-R_0}},\\
\sum_iQ_{c,i}q_i&\le CN2^{-R_0}\Lambda,&
\sum_i\sqrt{Q_{c,i}q_i}&\le C\sqrt{N2^{-R_0}\Lambda}.
\end{aligned}
\qquad(31)
```

These are geometric sums, not a sum bounded by R times its largest
term. For a nonempty middle, $`R_0=\lceil A_0\log_2(\ell+2)\rceil`$.
Increase the same fixed $`A_0`$ if necessary so that

```math
J2^{g-R_0}=O(L),\qquad 2^{-R_0}\Lambda=O(L).
\qquad(32)
```

Here $`J=O(L)`$, $`2^g=O(\ell^{1/8})`$ and
$`\Lambda=O(\log(\ell+2))`$ uniformly on the internal range.

The exact
[R90 whole-word query ledger](DIRTY_CATALOGUE_TERMINAL.md#7-the-local-catalogue-lemma-and-its-paid-query)
with v output bits gives

```math
T_Q=O(v+\sqrt{Qv}+Qv/B),\quad G_Q=O(Qv),\quad
D_Q=O(v+\log_2(a+2)+Qv/B^2),
\qquad(33)
```

provided its separate dirty pool satisfies
$`B\ge16c_I(a+2)`$. We use v equal to one for digit queries and
v equal to $`q_i`$ for coarse queries. The word-query depth's
$`q_i`$ term is retained; no improved wide-output query lemma is
needed.

By (16)–(17), the largest local code and factor costs satisfy

```math
q_g=O\!\left(\ell^{1/8}\log_2(\ell+2)\right)=o(\ell),
```

```math
S_g=O\!\left(\ell^{3/8}\log_2(\ell+2)\right)=o(\ell).
\qquad(34)
```

The largest catalogue has $`2^{q_g}`$ slots in each of two arrays,
and all its fixed factors are physically executed and charged. Its
complete width and counts obey

```math
q_g+2g2^{q_g}+2^{q_g}[S_g+gq_g]
=2^{O(\ell^{1/8}\log_2(\ell+2))}=n^{o(1)}.
\qquad(35)
```

Multiplying the counts by $`R=O(\ell)`$ keeps this
$`n^{o(1)}`$ bound. Parallel disjoint catalogue factors contribute
$`S_g`$ depth; the charged shared-control routing contributes
$`q_g`$ depth. This is a circuit on arbitrary arrays, not supplied
quantum memory. Groups may reuse the same physical dirty bank after each
complete return, with a layout for the current h; all peak sizes are
bounded by (35).

Reserve this program and both arrays simultaneously, plus all separate
local returned-helper banks. Their width is
$`o(\sqrt N/n)`$. For sufficiently large absolute n they use at
most $`b/4`$. Set $`B=\lfloor b/2\rfloor`$ for the shared
query pool. It excludes all live arrays, program, address, outputs,
protected source work and saved H. Enlarging the same cutoff ensures

```math
B\ge16c_I[n+C\log_2(\ell+2)+2],\qquad B=\Theta(b).
\qquad(36)
```

The remaining dirty fraction supplies any local returned conjunction
and arithmetic helpers. Their costs are fixed polynomials in ell,
including the local predicates on at most M logical bits, and are
negligible beside $`\sqrt N/n`$. Such local predicates have
logarithmic native T-depth. Source, history and matrix-unit arithmetic
on $`O(\log\ell)`$ fields have a safe polylogarithmic depth,
which is $`o(\ell)`$. Exact helpers are reused only after return;
no retained failure or occupied source coordinate is borrowed.

Substituting (31)–(36) into (33), then adding the catalogue and local
operations, gives

```math
T_{\rm mid}=O(\sqrt{NL}+NL/b),\qquad G_{\rm mid}=O(NL),
```

```math
D_{\rm mid}=O\!\left(NL/b^2+L+R\ell\right).
\qquad(37)
```

The $`O(L)`$ source/padding preparation and inverse cost is paid
only at the whole-stream boundaries. The constant number of recovery
appearances changes no order. Global zero-work reflections and the
saved predicate add only logarithmic depth. No g constituent layers are
executed through separate long-prefix coefficient queries: both the
coarse catalogue and (19) describe their complete group operator.

## 8. Bottom layers, errors and original-width composition

For the deepest $`R_0`$ layers use the retained
[one-clean smaller-error bottom](ONE_CLEAN_DEPTH.md#the-final-logical-layers-at-smaller-error),
with

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad Q_j=N2^{-j}.
```

Its original dirty base is $`L+n+9`$. On the wide branch the
remaining pool is $`\Theta(b)`$ and satisfies the sharp-query
minimum width beyond the same absolute cutoff. The geometric weighted
sums give

```math
T_{\rm bot}=O(\sqrt{NL}+NL/b),\qquad G_{\rm bot}=O(NL),
```

```math
D_{\rm bot}=O\!\left(NL/b^2+R_0[L+R_0+\ell]\right)
=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right).
\qquad(38)
```

Its error is below $`15\eta/64`$. Its source, guards and local
count overheads are polynomial in n and ell on this precision range and
fit the displayed count bounds. If R is zero, omit all middle machinery;
then $`M=R_0=O(\log(\ell+2))`$ and (38) applies directly.

All stages run in sequence, so their dirty peaks are maximized, not
added across time. The initialized flag is reused by complete-isometry
telescoping; actual later unitaries propagate earlier leakage. The wide
real error is bounded by

```math
\frac{\eta}{512}+\frac{\eta}{2048}
 +\frac{\eta}{128}+\frac{3\eta}{256}
 +\frac{15\eta}{64}\lt\frac\eta2.
\qquad(39)
```

The terms are the single sign reduction, the coarse stage, the local
early decoder, the recovered middle, and the bottom. Omitted stages
contribute zero. Combining (4), (11), (14), (37) and (38) proves (1)
for the wide real branch, uniformly through $`L\le3\ell`$.

For complex magnitude frames on $`L\le2\ell`$, use this real
wide construction at precision $`L+1\le3\ell`$, followed by the
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at the same precision $`L+1`$. The real error is below
$`\eta/4`$ by (39); the diagonal error is below $`\eta/2`$.
The diagonal adds

```math
T_D=O(\sqrt{NL}+L+NL/b),\qquad G_D=O(NL),\qquad
D_D=O(NL/b^2+L+\ell).
```

Both circuits fit the complex literal allocation
$`b\ge17(L+n+8)`$ on this wide branch. All supplied leaf phases,
including those at zero-amplitude leaves, are retained. The sum of
errors is below eta, and the precision shift preserves the displayed
resource orders. This is the prescribed phase-dressed magnitude family,
not arbitrary complex-unitary synthesis.

For $`b\lt\sqrt N/n`$, select the complete R79 compiler for the
respective family at the original precision L. Since

```math
NL/b^2\gt Ln^2\ge nL\ge F_{79}(n,L),
```

its entire additive depth is absorbed by the width term in (1).
It already includes the complex diagonal when appropriate; none is
appended again. For dimensions below the chosen absolute cutoff, use
that same complete original-precision compiler at any supplied width.
The displayed precision range then contains only finitely many n,L;
one absolute enlargement absorbs its bounded additive depth into
$`F_{93}(n,L)`$. These fallbacks retain the original family-specific
widths and exact stated tolerance, without substituting a stricter
internal precision into a literal reservation.

## 9. Scope of the improvement

The new step is a paid correction of an entire noncommuting group. Its
coarse code describes all h local layers, the residual matrix-unit query
excludes all h changing targets, and the later group reads its own
current preserved prefix. The returned two-array identity and the
retained source/failure/history word are used together on their complete
input spaces. A sequence of old single-target queries would retain the
unimproved $`L\ell`$ term.

At the benchmark $`L=\Theta(\ell)`$, the paid bottleneck in this
construction is $`R\ell`$ with
$`R=\Theta(\ell/\log\ell)`$. The generic fixed-unitary
synthesis bound grows as $`8^g`$ times logarithmic precision, so
$`g=\Theta(\log\ell)`$ fits within an $`O(\ell)`$
per-group depth. This explains the achieved upper bound, not a lower
bound on other mechanisms. Reaching the stronger target requires an
additional improvement to the structured group interpreter, global
prefix access, or the number of addressed groups.

The proof is uniform on its stated growing-precision range and charges
all catalogue entries, queries, native local gates, actual inverses and
simultaneous work. Finite preprocessing can be enormous, and no
practical crossover or scalable circuit emitter is asserted. Optimal
general depth remains the continuing destination.
