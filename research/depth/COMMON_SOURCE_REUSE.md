# Common-source reuse and rejected-work diagnostics

**Status:** exact conjugation identity and scoped failures of cached-success and fresh-bank replacements; no improved complete-frame bound. Sections 1–10 refer to the [grouped-program proof](../../docs/GROUPED_PROGRAM_PREFETCH.md); original numbering is retained below.

## 11. A common-source identity and the remaining reflection

Keep the same source core and branch throughout a cached group. Put

```math
D=VH_b,\qquad
S_\ell=C_{h=u_\ell=b=1}(K_\ell)P_{h,u_\ell,W},\qquad
E_\ell=\mathrm{CZ}(h,u_\ell).
```

The controlled target word commutes with V: it acts only on flags, the
inner enable, and the logical target. The programmed mask is **not**
commuted through V. Hence, with literal operator order,

```math
Q_\ell=D^\dagger S_\ell D,\qquad
\widehat R_\ell=D R_\ell D^\dagger,\qquad
B_\ell=E_\ell S_\ell\widehat R_\ell S_\ell^\dagger
                    \widehat R_\ell S_\ell,
\qquad A_\ell=D^\dagger B_\ell D.
```

On the invariant zero-preparation-work subspace, write
$`|\omega\rangle=|+\rangle_b\otimes|\psi_m\rangle_{\rm core}`$,
where $`|\psi_m\rangle`$ is the prepared geometric source. Then

```math
\widehat R_\ell
=I-2[h=1,u_\ell=1]\otimes|\omega\rangle\langle\omega|.
```

The Section 3 selector/enable compute word $`C_\ell`$ acts on local
logical and selector work only, so it commutes with D. Define
$`\overline B_\ell=C_\ell^\dagger B_\ell C_\ell`$.
The exact two-layer identity is

```math
(C_1^\dagger A_1 C_1)(C_0^\dagger A_0 C_0)
=D^\dagger\overline B_1\overline B_0D.
```

The corresponding g-layer identity follows by telescoping. Section 10
uses incremental selector maintenance instead of paired per-layer
computations. Every one of its actual selector and enable words also
commutes with D, so the same boundary conjugation follows through that
actual interleaving; it does not require restoring the old selector
schedule. Actual inverses and every enable phase remain unchanged.
These conjugation identities hold on the full physical space, including
arbitrary inactive helpers, when $`\widehat R_\ell`$ is defined by
conjugation. The rank-one display is restricted to the invariant
zero-work subspace. No fresh source per height is needed.

This moves the preparations to the group boundaries, but leaves two
conjugated reflections per layer. Synthesizing each by its displayed
conjugation retains $`O(\log(m+2))`$ depth per reflection. The identity
is a constant-factor reorganization, not an $`O(g+\log m)`$ group-depth
theorem. It does not allow discarding the preparations while leaving an
uncharged $`\widehat R_\ell`$ gate.

### A cached success bit does not survive the programmed word

Fix an active enabled row and put
$`\Pi=|\omega\rangle\langle\omega|\otimes I_{\rm target}`$.
A hypothetical extra monitor of this projector is

```math
C_\Pi=\Pi\otimes X_a+(I-\Pi)\otimes I_a.
```

This monitor is an audit device, not an additional uncharged compiler
flag. The subspace in which a records the correct success bit is
preserved by a word acting only on source and target if and only if
that word commutes with $`\Pi`$. Computational-basis invariance of
the diagonal mask is insufficient: $`\Pi`$ is not a
computational-basis projector.

There is a legal exact row witnessing a constant failure. Choose target
angle zero, cosine mask zero, and sine mask with only its first geometric
bit set. Since that bit has squared weight one half,

```math
c=1,\qquad s=0,\qquad \zeta=0.
```

For the resulting programmed mask P, its expectation on the source is

```math
\tau=\langle\omega|P|\omega\rangle=(1+0)/2=1/2,
\qquad \|[P,\Pi]\|=\sqrt{1-\tau^2}=\sqrt3/2.
```

The commutator norm follows by resolving $`P|\omega\rangle`$ into
its parallel and orthogonal components; P is a Hermitian involution.
This holds for every $`m\ge2`$. All amplitudes and phases are exact;
it is not an inadmissible coefficient pair outside the unit circle.
The full middle word S also satisfies $`\Pi S\Pi=\Pi/2`$, because
its target compression is $`(cI+sK)/2=I/2`$.
For every normalized target state v,

```math
\|(I-\Pi)S(|\omega\rangle\otimes|v\rangle)\|=\sqrt3/2.
```

Compute the monitor from $`a=0`$, apply S without updating it, and
apply the actual monitor inverse. The rejected source component leaves
a equal to one with norm $`\sqrt3/2`$. Thus a precomputed success bit
cannot be assumed to erase, and its stale value cannot supply the next
reflection. An exact update must mix the success and failure sectors;
writing the update as a conjugation transfers the original task instead
of removing it.

### A fresh bank per invocation does not supply work return

Fresh independent banks per **complete amplified layer** are legitimate
if their conditional-zero widths are charged. Fresh banks per Q
occurrence are different. If an initialized bank is used once by Q,
and every later operation commutes with that bank's success projector,
its rejected norm is invariant. Since

```math
J^\dagger QJ=(cI+sK)/2,\qquad
(J^\dagger QJ)^\dagger(J^\dagger QJ)=(c^2+s^2)I/4,
```

the rejected norm after that one use is
$`\sqrt{1-(c^2+s^2)/4}`$. It equals $`\sqrt3/2`$ even on the exact
zero-angle row above. Later operations on other banks, together with
phases controlled by this bank's success sector, cannot return it to
zero. The repeated interactions with the same bank in amplitude
amplification remain necessary unless a different return mechanism is
supplied.

The [bounded source-reuse checks](../../tests/test_grouped_source_reuse.py)
retain all source and preparation-work input columns for the two-stage
conjugation identity, including inactive inputs. Negative controls show
that keeping the original zero reflection after moving the preparations
is incorrect. An exact zero-angle fixture with native PREP checks both the stale
monitor's rejected norm and the one-use source leakage, while proper
amplification returns that row exactly. Masks and reflections in these
fixtures are reduced operators; they do not emit a scalable reflection.

These are obstructions to two specified replacement proposals. They do
not lower-bound arbitrary source circuits, dynamically updated monitors,
alternative encodings, initialized-clean frame compilation, or
large-width T-depth. The complete-frame frontier remains unchanged. A
positive next proposal must supply a shallow literal conjugated
reflection or a different complete native identity, with every monitor
update and all work return charged.
