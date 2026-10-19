# Dirty-tree compilation from complete branch differences

[Current workspace](../../WORKSPACE.md) ·
[Exact nonuniform tree](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Approximate filter limits](APPROXIMATE_INDICATOR_FILTER_LIMITS.md)

An individual tree source need not implement a complete controlled
query, or return its private dirty work after each use. Its two parent
branches need the correct **relative operator** on the entire child and
private-work space. This weaker promise has a complete error, resource,
and return recurrence. It permits exact cancellation of common work
actions and scalar sibling phases, while retaining the accumulation of
independent erroneous child actions.

The recurrence is a positive compiler component, not a new fast source.
The native branch-difference operation remains to be constructed. It
does not imply a lower bound on unrestricted native depth or supply a
reflection implementation.

## 1. A source contract on arbitrary dirty blocks

Use the tree of the [exact indicator](NONUNIFORM_DIRTY_INDICATOR.md#1-the-exact-tree-interface-permits-unequal-chunks).
The address has $`p\ge1`$ bits, the output has $`m=2^p`$ bits, and positive
chunk lengths $`b_1,\ldots,b_t`$ sum to p. Put
$`d_i=\sum_{j=1}^i b_j`$ and $`d_0=0`$. Let $`F_i`$ be the
ideal stage toggling each child by its parent times its chunk-equality
predicate. Products act right to left. Define

```math
F_{\le i}=F_i\cdots F_1,\qquad
K_i=F_{\le i}X_\emptyset F_{\le i}^\dagger,\qquad
K_0=X_\emptyset.
```

For each computational address x, $`K_i`$ flips the root and the
single selected node at each boundary through i, and is identity on
all private work. This holds on every initial tree word.

An actual stage $`V_i`$ must preserve the computational address and
every parent at boundary $`i-1`$. It may target children and private
dirty work, and must not depend on earlier ancestors. Conditional on
address x and parent bits $`r_a`$, require

```math
V_i(x,r)=
\bigotimes_{a\in\mathbb F_2^{d_{i-1}}}V_{i,a,x}(r_a).
                                                               \qquad(1)
```

The factors have disjoint child/private blocks within this stage.
Shared address and parent controls are permitted. This is an operator
contract; an actual simultaneous native schedule remains separately
required.

The selected parent is $`a_i(x)=x_{1:d_{i-1}}`$. On its block let
$`S_{i,x}`$ be X on the child selected by the current chunk, tensored
with identity on all other children and private work. Define

```math
\delta_i=\max_x
\left\|
V_{i,a_i(x),x}(1)V_{i,a_i(x),x}(0)^\dagger-S_{i,x}
\right\|_{\mathrm{op}}.                                    \qquad(2)
```

The norm includes every arbitrary child/private input. In particular,
the source may have $`V(0)=R`$ and $`V(1)=SR`$ for any unitary R;
then its differential error is zero. The order matters: common left
multiplication instead gives $`RSR^\dagger`$.

## 2. Error and complete return

Put

```math
\widetilde K_i=(V_i\cdots V_1)X_\emptyset(V_i\cdots V_1)^\dagger,
\qquad e_i=\|\widetilde K_i-K_i\|_{\mathrm{op}}.
```

Then

```math
e_0=0,\qquad e_i\le e_{i-1}+\delta_i.                         \qquad(3)
```

Insert $`V_iK_{i-1}V_i^\dagger`$ in the difference. Unitary
invariance bounds the propagated error by $`e_{i-1}`$. On each x
block the ideal $`K_{i-1}`$ flips exactly the selected parent among
this stage's controls. All other parent factors in (1) cancel literally.
Earlier ancestor flips commute through $`V_i`$. For selected-parent
input r, the remaining child action is

```math
V_{i,a_i(x),x}(r\oplus1)V_{i,a_i(x),x}(r)^\dagger.
```

Its two branches are the differential in (2) and its adjoint. Because
$`S_{i,x}`$ is Hermitian, both have distance at most $`\delta_i`$
from the ideal action. This proves (3), without summing over parents.

Let C copy every leaf into its corresponding arbitrary output bit by
CNOT. The chronological complete word is

```math
C,\quad\widetilde K_t,\quad C,\quad\widetilde K_t^\dagger.
```

Its ideal counterpart is $`I_p`$ on address/output and identity on
all tree/private wires. Replacing the two conjugates gives

```math
\epsilon_{\mathrm{ind}}
\le2e_t\le2\sum_i\delta_i.                                  \qquad(4)
```

All inverse occurrences are the actual reversed native words. Thus (4)
includes dirty-return and leakage error and holds with arbitrary
external references. If all differentials are exact, every dirty wire
returns exactly, even though the individual stages can leave arbitrary
intermediate private-work action.

A private pool may be reused between stages. The proof compares
unitaries on the full common workspace, and the ideal $`K_{i-1}`$
is identity on that pool. Its state is arbitrary at every stage; there
is no initialization or local reset premise. Disjointness of the
parent blocks within a stage is still required.

More generally, with $`E_i=F_i^\dagger V_i`$, the stage comparison
is exactly $`\|[E_i,K_{i-1}]\|`$. Equation (2) makes one useful
cancellation condition explicit rather than assuming $`V_i\approx F_i`$.

## 3. Which sibling errors cancel

Suppose a parent block further factors into row words $`W_v(r)`$,
on disjoint child/private blocks for $`v\in\mathbb F_2^{b_i}`$.
Define their differential errors by

```math
\gamma_{i,v}(x)=
\left\|W_v(1)W_v(0)^\dagger-X_v^{[x^{(i)}=v]}\right\|.
```

The parent differential is a tensor product, so

```math
\delta_i\le\max_x\sum_v\gamma_{i,v}(x).                     \qquad(5)
```

If each whole row approximates its ideal controlled query to error
$`\epsilon_i`$, then $`\gamma_{i,v}\le2\epsilon_i`$ and

```math
\epsilon_{\mathrm{ind}}\le4\sum_i2^{b_i}\epsilon_i.           \qquad(6)
```

The exponent is the chunk length $`b_i`$, not the cumulative prefix
length $`d_i`$. All other parents' errors cancel. Incorrect child
actions under the selected parent still act on arbitrary dirty children;
their ideal one-hot predicates do not give a maximum-over-children bound.

For example, suppose the row differential is its ideal child X times
$`\exp(i\theta_vX_v)`$ and $`\sum_v|\theta_v|\le\pi`$.
The exact parent-block error is

```math
2\sin\!\left(\frac12\sum_v|\theta_v|\right).                \qquad(7)
```

Independent child-X eigenvalues can choose the signs of all angles.
The extremal sum attains (7), and no other eigenvalue has a larger
distance in this angle range. Thus a square-sum or single-child estimate
is invalid. The retained [additive filter bound](APPROXIMATE_INDICATOR_FILTER_LIMITS.md#1-the-literal-additive-phase-filter-contract)
already supplies a sharper specialized obstruction for low-degree
Hamming kernels, including phase wraparound.

Scalar sibling phases can instead sum to a cheap correction. If

```math
W_v(r)=e^{i\theta r[x^{(i)}=v]}X_v^{r[x^{(i)}=v]},
```

then the complete parent block is $`e^{i\theta r}S_{i,x}^r`$,
because the equality predicates sum to one. A parent correction
$`\mathrm{diag}(1,e^{-i\theta})`$ makes that block exact.
For $`\theta=-\pi/2`$ it is one native Clifford S. This permits
a relative-phase source to be corrected jointly; it does not supply
a growing-arity source or a sublogarithmic depth bound.

## 4. Native resources, width, and the actual inverse

Stage i has $`n_i=2^{d_i}`$ rows. Suppose each row has common
forward/actual-inverse bounds $`c_{T,i},c_{G,i},w_i`$, and a supplied
simultaneous stage-depth bound $`h_i`$. Individual row depth alone
does not imply $`h_i`$ when shared controls are T targets. The retained
private-phase-target/read-only-control interface supplies a valid
schedule; another interface must prove its own schedule.

The complete echo has four source traversals. Therefore

```math
\begin{aligned}
T_{\mathrm{ind}}&\le4\sum_i n_ic_{T,i},\\
G_{\mathrm{ind}}&\le4\sum_i n_ic_{G,i}+2m+8,\\
D_{T,\mathrm{ind}}&\le4\sum_i h_i,\\
w_{\mathrm{ind}}&\le(2m-1)+\max_i n_iw_i,\\
\epsilon_{\mathrm{ind}}&\le2\sum_i\delta_i.
\end{aligned}                                                   \qquad(8)
```

The root and every live tree boundary are included in $`2m-1`$;
the original output is excluded from additional work. The final term
in G pays two leaf readouts and two literal $`X=HSSH`$ words. Use
a common maximum when forward and inverse Clifford counts differ,
or sum their separately charged orientations for a sharper bound.
Any per-parent corrections or other source work must also be included.

Local return is unnecessary under (2), but global return is exactly
the contract in (4). All resources are for actual words, with no free
controlled call or unimplemented reflection.

## 5. The unpaid operation after a sparse first stage

Choose a first chunk $`b=p-h`$ and use the retained exact nonuniform
tree for the remaining h bits. Its replicated exact suffix has
$`O(m)`$ count/width and $`O(\log(h+2))`$ depth. Suppose a row
differential primitive on b literals has costs
$`C_T(b,\eta),C_G(b,\eta),W(b,\eta)`$ and compatible batch depth
$`H(b,\eta)`$. The first stage has $`2^b`$ rows. It suffices to use

```math
\eta\le\tau/2^{b+1}
```

for differential error, or $`\epsilon\le\tau/2^{b+2}`$ for ordinary
complete-row error. The full recurrence is

```math
\begin{aligned}
T&=O(m)+O(2^{p-h}C_T),\\
G&=O(m)+O(2^{p-h}C_G),\\
w&=O(m)+O(2^{p-h}W),\\
D_T&\le4H+O(\log(h+2)),\qquad\epsilon_{\mathrm{ind}}\le\tau.
\end{aligned}                                                   \qquad(9)
```

If these row costs are polynomial in p at that precision, choose
$`h=\lceil A\log_2(p+2)\rceil`$ for a sufficiently large fixed A.
Then the first-stage count/width are absorbed into $`O(m)`$, and a
true $`H=o(\log p)`$ would give the desired depth, with an exact
$`O(\log\log p)`$ suffix. The required row precision satisfies
$`\log(1/\eta)=\Theta(p)`$. Fixed row accuracy does not suffice.
The exact retained conjunction instead has $`H=O(\log p)`$, so
substituting it recovers the current frontier.

The stronger [polynomial-overhead indicator fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
uses a whole sibling-bank source at **constant** complete error, then
absorbs its polynomial count/width overhead in this sparse first stage.
It avoids the independent-row precision demand in (9). The branch
difference recurrence explains the valid cancellation and the remaining
contract; the fusion theorem supplies the sharper sufficient source
interface.

A third source band alone supplies neither interface. The missing
operation is an actual growing-arity native word with a complete joint
differential, the required depth and compatible batch schedule. Even
the relaxed differential cannot hide a free query. For a one-parent
stage, including the sparse first stage, the two-call word
$`VX_{\mathrm{parent}}V^\dagger X_{\mathrm{parent}}`$ is within
$`\delta`$ of $`S\otimes I_{\mathrm{parent}}`$ when (2) holds.
Its two parent sectors carry the differential and its adjoint. A fixed
parent in a multiple-parent stage would require this bound at every
address, beyond its selected-address promise in (2). Thus
the relaxation removes unnecessary local return and phase obligations,
but still pays the complete controlled action.

The [two bounded fixtures](../../tests/test_dirty_tree_differential_recurrence.py)
retain literal native phases and every physical column. They check
nonreturned private work, both parent orientations, cancellation under
unselected parents, the complete root/leaf echo, and the independent
child-error accumulation in (7). They do not supply the missing source
or extrapolate a growing-depth theorem.
