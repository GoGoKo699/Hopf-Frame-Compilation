# A complete bank from one nonlinear algebra source

[Current checkpoint](../../WORKSPACE.md) ·
[Paid scalar extraction and fast source](ALGEBRA_QUERY_EXTRACTION.md) ·
[Balanced HRT permutation](BALANCED_ALGEBRA_PERMUTATION.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

The [current direct output readout](OUTPUT_ALIASED_CHARACTER_SOURCE.md)
identifies the character coordinates with the output bank and uses four
source calls, with full error at most four times source error. It also
permits arbitrary quantum private-helper fibers, which cancel inside the
completed character. The eight-call disjoint-carrier construction below
is retained with its original proof and fixture.

The earlier algebra query was compared with the square-root count and
width budget for one arbitrary table. The complete-indicator relaxation
changes that comparison. One shared coefficient carrier now supplies
**all** equality outputs, with linear Clifford preprocessing of the
arbitrary output bank. It does not require one carrier per output.
The paid fast HRT source consequently fits the relaxed count and width
budgets, but retains logarithmic address depth.

This is a concrete mechanism outside the Majorana-preserving bank-motion
representation. Its next unpaid operation is a sublogarithmic native
implementation of the complete HRT source, or of the completed extraction
word. Neither the source's third-level hierarchy membership nor its
quadratic inverse pays that operation. No unrestricted depth bound or
improvement to R96 follows from the present construction.

## 1. The shared complete source

Let $`b\ge1`$, $`m=2^b`$ and $`r=b+1`$. The original address x has b bits,
and Y is an arbitrary m-bit output bank indexed by
$`u\subseteq[b]`$. All additional qubits are arbitrary dirty inputs.
Use the square-zero coefficient algebra

```math
R_r=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2).
```

Its nonconstant coefficient word has $`2^r-1=2m-1`$ physical qubits.
Every such qubit is occupied arbitrarily; the constant coefficient is a
formal algebraic constant, not a clean physical wire. Write

```math
f(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1,
\qquad F_r|c\rangle=|f(c)\rangle.
```

The [retained HRT construction](BALANCED_ALGEBRA_PERMUTATION.md#6-balanced-hrt-construction)
proves this is a permutation on the entire coefficient space and gives
its actual native inverse. Its full inverse coordinate map is quadratic:

```math
k(a)_S=a_S+
\sum_{\substack{\{I,J\}:I,J\ne\varnothing\\I\sqcup J=S}}a_Ia_J.
```

The product defining f immediately gives

```math
f(c+d)=f(c)+f(d)+f(c)f(d).
```

Each factor squares to one, and all factors commute. This identity
holds for all coefficient inputs, with no code-subspace qualification.

Define the address translation $`X^{\iota(x)}`$ by one CNOT from each
$`x_i`$ into the singleton coefficient $`c_{\{i\}}`$, where
$`\iota(x)=\sum_{i=1}^b x_i e_i`$. Set

```math
J_x=F_rX^{\iota(x)}F_r^\dagger,\qquad
A_x=\prod_{i=1}^b(1+x_i e_i)-1.
```

On every carrier input y the exact action is

```math
J_x:y\longmapsto y+A_x+yA_x.
```

This follows by writing $`y=f(c)`$ and applying the displayed product
identity. The generator $`e_r`$ does not occur in $`A_x`$.
Let $`\lambda(y)=[e_{[r]}]y`$ and let Z be the physical Z gate on
that top coefficient. The chronological word

```math
C_x:\quad J_x,\ Z,\ J_x^\dagger,\ Z
```

therefore has the exact full-space action

```math
C_x|y\rangle=(-1)^{\lambda(yA_x)}|y\rangle.
```

The carrier returns after this word. The two Z phases sum to
$`\lambda(y+A_x+yA_x)+\lambda(y)`$, and $`\lambda(A_x)=0`$.
The character is a Pauli on the carrier at each fixed address; its
coherent address dependence is not being declared Clifford.

## 2. Pay every output jointly

The equality function at u has the algebraic normal form

```math
[x=u]=\prod_{i\in u}x_i\prod_{i\notin u}(1+x_i)
     =\sum_{S\supseteq u}\prod_{i\in S}x_i.
```

Its constant coefficient is one exactly at $`u=\varnothing`$.
For nonempty S put

```math
z_S(Y)=\bigoplus_{u\subseteq S}Y_u,
\qquad
v(Y)=\sum_{\varnothing\ne S\subseteq[b]}z_S(Y)e_{[r]\setminus S}.
```

All carrier slots appearing in v contain the dummy generator and are
nonempty. Let $`L_Y:y\mapsto y+v(Y)`$, preserving x and Y.
This is an exact Clifford shear. It is paid without initialized buffers:

1. On Y, execute its in-place binary subset-zeta transform. For each
   $`i\in[b]`$ and $`S\ni i`$, add $`Y_{S\setminus\{i\}}`$ into
   $`Y_S`$, in increasing direction order. This costs $`bm/2`$ CNOTs.
2. Copy each transformed $`Y_S`$, for nonempty S, into the existing
   carrier coordinate $`y_{[r]\setminus S}`$. This costs $`m-1`$ CNOTs.
3. Reverse the actual zeta CNOT word, at another $`bm/2`$ CNOTs.

The completed L and its actual inverse each cost

```math
G_L=bm+m-1,\qquad T_L=D_T(L)=0,
```

with no extra helper. Y returns exactly before the next source word.
The temporary change of Y is allowed on every coherent bank input;
it is not a copy of Y into a clean history.

Now execute chronologically

```math
C_x,\ L_Y,\ C_x^\dagger,\ L_Y^\dagger.
```

The exact two character phases differ by

```math
\begin{aligned}
\lambda(v(Y)A_x)
 &=\sum_{\varnothing\ne S\subseteq[b]}z_S(Y)
                                     \prod_{i\in S}x_i\\
 &=\sum_{u\subseteq[b]}Y_u\bigl([x=u]+[u=\varnothing]\bigr)\\
 &=Y_x+Y_\varnothing.
\end{aligned}
```

The carrier and all source helpers return. Conjugate the completed word
by Hadamards on all Y wires, then apply X to $`Y_\varnothing`$. The
result is the literal complete indicator

```math
I_b|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle.
```

The constant correction is essential: the character construction itself
has no nonzero constant character. There is no promise that Y, any
carrier coordinate or any private helper starts at zero. The displayed
operator equality includes arbitrary references entangled with all of
them.

## 3. Complete resource, error and return recurrence

Suppose a native $`\widetilde F_r`$, acting on the carrier and $`w_F`$
private arbitrary helpers, satisfies the full-operator source promise

```math
\|\widetilde F_r-F_r\otimes I_{w_F}\|_{\mathrm{op}}\le\eta.
```

Use its **actual adjoint** wherever $`F_r^\dagger`$ occurs. Let
$`T_F,G_F,d_F`$ be simultaneous common bounds for the source and that
actual adjoint. This is a complete source promise, not an accuracy
statement on singleton-supported or zero-filled coefficient inputs.

One C word uses four source or actual-adjoint occurrences. The completed
indicator uses one C and its actual adjoint, hence eight occurrences in
total. Unitary product telescoping gives

```math
\|\widetilde I_b-I_b\otimes I_w\|_{\mathrm{op}}\le8\eta.
```

This includes every address, output, coefficient-carrier and private-work
error, and arbitrary reference entanglement. Intermediate approximate
source calls need not return the private pool exactly. Each subsequent
word is unitary on that current arbitrary state, so the hybrid argument
does not discard leaked amplitude or assume fresh work. The actual
completed inverse has the same error, because the exact indicator is an
involution.

The complete common native ledger is

```math
\begin{aligned}
T_I&\le8T_F,\\
G_I&\le8G_F+2bm+4b+4m+10,\\
D_T(I)&\le8d_F,\\
w_I&\le2m-1+w_F,\\
\varepsilon_I&\le8\eta,\qquad c_{\mathrm{clean}}=0.
\end{aligned}
```

The Clifford remainder charges four b-CNOT address translations, two
completed L words, four top-coordinate Z gates, two m-Hadamard envelopes
and one X correction. Z uses $`S^2`$ and X uses HSSH. Their reversed
adjacent powers can be reduced literally by $`S^{-2}=S^2`$, giving the
same common allowance in the positive-S alphabet. Nothing is charged as
an address-dependent free Clifford. The source-private pool is reused
across sequential source calls and stays disjoint from the live carrier,
x and Y; the carrier's $`2m-1`$ qubits are explicitly included in width.

If the source implements F exactly, the completed query is exact and
returns every carrier and private wire. An approximate source supplies
only the full error bound above. Exact private-helper return of an
approximate source alone would not prove exact carrier return of the
completed commutator; actual adjoints do not supply that missing identity.

## 4. The retained fast source meets the relaxed budgets

The [paid fast subset-product source](ALGEBRA_QUERY_EXTRACTION.md#5-a-smaller-count-complete-hrt-source)
has, for $`r\ge2`$,

```math
h_r=\lceil\log_2r\rceil,\qquad
A_r=\binom r2\,2^r,
```

```math
T_F\lt336h_rA_r,\qquad G_F\lt978h_rA_r,\qquad
 d_F\le200h_r,\qquad w_F\lt9A_r,\qquad\eta=0.
```

All these are simultaneous bounds on an actual complete native source
and its reversed inverse. The returned subset-convolution masks and
linear pre/postprocessing have already been paid in that proof; no new
arithmetic oracle is being assumed here. Substitution gives

```math
\begin{aligned}
T_I&\lt2688h_rA_r,\\
G_I&\lt7824h_rA_r+2bm+4b+4m+10,\\
D_T(I)&\le1600h_r,\\
w_I&\lt2m-1+9A_r,\qquad\varepsilon_I=0.
\end{aligned}
```

Thus the one-carrier bank has
$`T,G=O(2^b b^2\log(b+2))`$,
$`w=O(2^b b^2)`$ and $`D_T=O(\log(b+2))`$.
These counts and width fit the polynomial-overhead complete-source
window. In particular, the old square-root-width comparison for one
arbitrary table is not an obstruction to this full-bank use. The depth
remains the old logarithmic order, so inserting this exact source into
root fusion does not change the full resource frontier.

For any improved source at requested error $`\eta`$, root fusion makes
eight indicator-source calls. Its combined source contributions are
therefore 64 F/actual-adjoint occurrences, error at most $`64\eta`$,
and depth at most

```math
64d_F(b_0+1)
       +O_\tau(\log\log(p+4)).
```

Here $`b_0=p-s`$ and $`s=\lceil a\log_2(p+2)\rceil`$ are the
actual large-prefix and exact-suffix lengths chosen by root fusion,
with fixed a large enough to absorb the source overhead.

This holds provided the complete source retains fixed polynomial overhead
in count and dirty width at each requested accuracy. Choosing
$`\eta=\tau/64`$ pays the final indicator tolerance tau. The fusion's
exact small-address suffix and full work-return ledger remain as proved
in its own proof home. Its finitely many short-address fallback cases
need no absent HRT call.

## 5. The quadratic inverse has a paid XOR query, not a paid in-place map

The distinction can be made explicit. Let $`w=2^r-1`$ and write an
arbitrary coefficient input as

```math
v_i(a)=\sum_{\min S=i}a_Se_S,\qquad
q_i(a)=\sum_{\min S\gt i}a_Se_S.
```

For a disjoint unordered pair of nonempty subsets, exactly one member
has the smaller minimum label. Therefore the quadratic inverse satisfies

```math
k(a)=a+\sum_{i=1}^{r-1}v_i(a)q_i(a).
```

Each factor is a coordinate projection of the same arbitrary input a.
Apply the [returned fast subset-product phase word](ALGEBRA_QUERY_EXTRACTION.md#4-fast-subset-convolution-as-a-native-returned-query)
to each of these $`r-1`$ products, with the common arbitrary output bank z
as its phase-readout input. Give the products disjoint masks and linear
program temporaries. All a and z wires occur only as Clifford controls
inside those phase words; their physical T targets are private masks.
Consequently the corresponding forty T layers run in parallel even
though the factors share source wires and their final phase readout.
Use one enclosing Hadamard pair on each of the w output wires, then add
the linear a term by w CNOTs.

No physical constant-coordinate wire is supplied. Every constant
coefficient of the two factors is structurally zero and is omitted from
the linear forms. Their product's constant coefficient is also zero;
the absent constant output requires no Hadamards. With

```math
\rho_r=\binom{r+2}{2}\,2^r,\qquad
s_r=r(r+1)2^{r-1},
```

this gives the complete exact query

```math
Q_k:|a,z,h\rangle\longmapsto|a,z\mathbin\oplus k(a),h\rangle
```

with common forward/actual-inverse bounds

```math
\begin{aligned}
T_{Q_k}&\le56(r-1)\rho_r,\\
G_{Q_k}&\le(r-1)(96\rho_r+64s_r)+3w,\\
D_T(Q_k)&\le40,\\
w_{Q_k}&\le(r-1)(3\rho_r+s_r),\qquad\varepsilon_{Q_k}=0.
\end{aligned}
```

The w wires of z are the query's arbitrary output and are excluded from
its extra-helper count; a use of this query as a source-building
component must separately allocate them. All masks and temporaries
return exactly, including their correlations with a, z and references.
The actual word uses H, CNOT, T and T-dagger, so the reversed word has
the same simultaneous ledger. These bounds fit polynomial overhead in
$`2^r`$ and do not assume initialized arithmetic scratch.

This pays the complete out-of-place quadratic XOR. It does **not** pay
$`|a\rangle\mapsto|k(a)\rangle`$ on the original coefficient register.
The arbitrary output is $`z_{\mathrm{old}}+k(a)`$, and the original a
has not disappeared. Swapping the two registers retains both the old-z
contamination and a; the supplied query does not erase them. Calling
that swap an in-place inverse is the first unsupported shortcut in this
version of the proposal. A repair must provide a separate full-input
identity, with its native resource and return ledger. The retained
balanced source provides one logarithmic-depth repair. The argument
here is not an impossibility theorem for other conversions.

## 6. The Choi basis is paid; its coherent X correction is not

The complete quadratic query also gives a constant-depth change of basis
on **two arbitrary** w-bit registers. Put $`d=2^w`$ and define, in
operator order,

```math
G=\mathrm{SWAP}_{A,B}\,Q_k\,(H_A\otimes I_B).
```

Its literal action on every computational label pair is

```math
\begin{aligned}
G|a,z\rangle
 &=d^{-1/2}\sum_y(-1)^{a\cdot y}|k(y)+z,y\rangle\\
 &=(X^z\otimes Z^aF_r)|\Phi\rangle,\qquad
|\Phi\rangle=d^{-1/2}\sum_v|v,v\rangle.
\end{aligned}
```

The Bell vector on the right describes these basis vectors; it is not an
initialized resource assumed by the implementation. G is the explicit
all-input unitary on the left. The arbitrary private work of $`Q_k`$
returns exactly. Its common bounds are

```math
T_G=T_{Q_k},\qquad d_G\le40,\qquad
G_G\le G_{Q_k}+4w,\qquad w_G\le w_{Q_k}.
```

The additional w Hadamards and w swaps cost $`4w`$ elementary Cliffords.
All source registers a,z remain full logical inputs, rather than clean
work whose allocation is omitted. The actual G adjoint uses the actual
reverse of the query and these Clifford words.

The HRT structure identifies the two literal Pauli corrections:

```math
\begin{aligned}
F_rZ^tF_r^\dagger|y\rangle&=(-1)^{t\cdot k(y)}|y\rangle,\\
R_z|y\rangle:=F_rX^zF_r^\dagger|y\rangle
 &=|y+f(z)+yf(z)\rangle.
\end{aligned}
```

The first has a paid **coherent-key** implementation, not merely cheap
fixed-key branches. Use $`H_tQ_k(y;t)H_t`$, giving common bounds

```math
T_Z=T_{Q_k},\qquad d_Z\le40,\qquad
G_Z\le G_{Q_k}+2w,\qquad w_Z\le w_{Q_k},
\qquad\varepsilon_Z=0.
```

Every t and y is arbitrary. In contrast, the second correction requires
f(z), of degree r, in the coefficients of its affine action on y.
Each fixed-z branch is a Clifford, but the required controlled operator
is the full word

```math
R:|z,y,h\rangle\longmapsto|z,y+f(z)+yf(z),h\rangle.
```

Its joint native synthesis is not supplied by that branchwise fact.
For example, a teleportation calculation with an additional clean Bell
resource would produce branches $`Z^aF_rX^z|\psi\rangle`$ after G adjoint.
The correcting word is $`R_zZ^a`$, with the ordinary Z applied first.
Thus even that calculation still asks for this coherent R. The actual
model supplies no such clean Bell resource, and the calculation proves
no returned-work identity for arbitrary initial Bell labels. Preparing
the complete G basis alone resolves neither obligation. A fused quantum
word may use a different mechanism; the statements here do not exclude
it.

### A complete coherent X correction already supplies the indicator

The missing R is a sufficient native operation in its own right. Retain
$`r=b+1`$, $`m=2^b`$ and $`w=2m-1`$. Allocate one arbitrary w-bit key z
beside the arbitrary w-bit coefficient carrier y. Let $`S_x`$ add
$`\iota(x)`$ to the singleton coordinates of z, at b CNOTs. Execute
chronologically

```math
R,\ S_x,\ R^\dagger,\ S_x^\dagger.
```

The exact representation law
$`R_zR_{z'}=R_{z+z'}`$ gives carrier action
$`R_{z+\iota(x)}^\dagger R_z=R_{\iota(x)}=J_x`$ and restores z.
No nonsingleton key coordinate is initialized or ignored. Replace each
of the four J/actual-adjoint occurrences in Sections 1--2 by this echo
or its actual inverse. A supplied full-operator approximation to R,
with error eta and simultaneous common bounds
$`T_R,G_R,d_R,w_R`$, then gives

```math
\begin{aligned}
T_I&\le8T_R,\\
G_I&\le8G_R+2bm+8b+4m+10,\\
D_T(I)&\le8d_R,\\
w_I&\le4m-2+w_R,\\
\varepsilon_I&\le8\eta,\qquad c_{\mathrm{clean}}=0.
\end{aligned}
```

The additional carrier and key occupy $`2w=4m-2`$ arbitrary wires.
There are eight b-CNOT key translations, two completed zeta/copy shears,
four top-coordinate Z gates, two output-Hadamard envelopes and the
constant output correction. Every source occurrence is R or its actual
adjoint. Two-source telescoping bounds each J echo by $`2\eta`$;
the four completed occurrences give $`8\eta`$ without assuming exact
intermediate key or helper return. All register and reference errors
are included in that full-operator bound.

Using the retained source directly gives R with two F/actual-adjoint
calls and w controlled-X CNOTs. That implementation retains logarithmic
depth. The new reduction does not make the correction free: it identifies
exactly which coherent Clifford family would suffice if compiled at
sublogarithmic depth within the relaxed budgets. It is not a lower bound
on a whole teleportation word, on a different completion outside the
stated full-key interface, or on unrestricted native depth.

## 7. First unpaid operation and exact scope of the retained limit

The next concrete question is to compile $`F_{b+1}`$ on its **entire**
$`2^{b+1}-1`$-qubit coefficient carrier, with arbitrary returned work,
in $`o(\log(b+2))`$ native T-depth and the displayed polynomial-overhead
budgets, at every fixed requested full-operator accuracy. Equivalently
for this sufficient route, a new compiler may synthesize the whole
completed eight-call word at qualifying depth; it need not preserve
separate source occurrences. Algebraic hierarchy membership, a shallow
out-of-place quadratic XOR, or the source's quadratic inverse does not
supply either compiler.

The balanced source currently implements successive complete nonlinear
merges. Its earlier
[classical-stage inverse-degree argument](ALGEBRA_QUERY_EXTRACTION.md#deferring-arbitrary-helper-garbage-does-not-remove-inverse-degree)
remains useful but scoped. If an exact all-input computation has principal
output $`k(c)`$ and arbitrary invertible helper action, its actual inverse
has principal coordinate $`f(y)`$ of degree r. A word of D completed
classical reversible stages whose inverse coordinate degrees are at most
$`\Delta\gt1`$ must therefore obey

```math
D\ge\lceil\log_\Delta r\rceil.
```

This argument permits arbitrary affine coordinate changes and arbitrary
helper garbage within that classical promise. It does not apply to
arbitrary intervening Hadamards or to approximate quantum words, and D
is not automatically native T-depth. Whole-word interference can leave
its hypotheses. Likewise, the new algebra source is not subject merely
by name to the Majorana-preserving Clifford-boundary theorem: it uses
nonlinear physical coefficient permutations, and its actual Clifford
boundaries are not granted that grammar.

The revised assessment is therefore precise. Complete nonlinear
extraction, all-output processing and their error/return recurrence are
paid. The relaxed count and width frontier is met. Sublogarithmic native
depth remains unpaid. A faster source or completed word would improve
the conditional fixed-accuracy frame law through the retained fusion and
R98 transfer; the present paid construction does not. The focused
original physical-pair lower scout concerns a different circuit class
and remains independent of this sufficient upper-bound interface.

## 8. Focused composition evidence

The [single focused fixture](../../tests/test_full_bank_character_transport.py)
uses b=2, the four arbitrary output bits and all seven arbitrary carrier
bits. It checks all 8192 computational phase columns of the resulting
13-wire composition and its actual reverse, including the constant
correction and the complete return of Y and the carrier. The zeta/copy
shear is expanded into its eleven actual CNOTs and checked separately
inside the same all-column test. Equality of every phase column implies
the corresponding operator identity with arbitrary references; the
analytic output-Hadamard conjugation supplies the XOR form.

F and its inverse in this fixture are the retained exact coefficient
permutations. The fixture tests their new composition and output wiring;
it is not a replacement native-source receipt or a finite extrapolation
of the asymptotic theorem. The earlier balanced-source and subset-product
fixtures, together with their analytic proofs, supply the native F
implementation and resource bounds used above.
