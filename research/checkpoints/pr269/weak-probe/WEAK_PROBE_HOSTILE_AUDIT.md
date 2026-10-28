# Independent hostile audit: weak-probe Pauli-domain continuation

Scratch only, 28 October 2026. No PR269 or repository edits.
Audited `domain-next/weak-probe/WEAK_PROBE_PAULI_DOMAIN.md` plus Root's
subsequent robust low-weight argument.

## Verdict

The exact field theorem and robust low-weight refinement are valid under the
stated complete-probe contract. No counterexample or phase/orientation defect
was found. The following scope clarification is important: the conclusions
apply when the decoded endpoint has a moving helper factor, such as
W_m tensor Q_delta. The endpoint is not required to fix those helpers; the
probe is Z_x tensor their identity. The source common-fixed-kernel assumption
is a separate condition, and q counts all target qubits beyond the m-wire bank.

## Exact theorem audit

Let n=m+q, m=2^p, R=E† Z_0 E, R_x=Q_x R Q_x†, and require
E R_x E†=Z_x tensor I_q for every x. Let H_E be the exact Pauli domain.
The earlier scalar-cube argument makes every image-domain X mask affine for
a degree-one query. Its orientation uses Q_x† P Q_x, which equals Q_x P Q_x†
for the literal involutive field family.

For K_Z={v:E†(Z^v tensor I_q)E is a Pauli}, exact covariance gives

    Q_y P_v Q_y†=P_(Pi_y v).

K_Z is a binary subspace, and its Pauli-label map phi is injective: any
nonzero bank Z word is traceless, so its unitary pullback cannot be scalar.
Scalar Pauli phases do not affect this injectivity or linearity.

The source identities N_y N_z=0 transfer through phi to second differences
of every v in K_Z. Thus v is affine. For p>=2, an affine scalar mask has a
nonzero translation stabilizer. Its source Pauli label is then fixed by one
nonzero query, hence by all queries because their fixed-label kernels coincide.
Injectivity forces v to be translation-invariant, hence constant. Therefore

    dim K_Z<=1,
    dim H_E<=p+2+2q,
    absorbed-axis rank and T(E)>=2m-p-2.

The projection of the image Pauli domain to (bank X mask, all helper labels)
has exactly K_Z as its kernel. There is no missing sign sector: signed
Hermitian Pauli pullbacks are included in the domain by definition.

For E=L(D^tensor(m/4) tensor I_q)K with arbitrary fixed joint Cliffords,
the domain contains dimension m+2q. Both outer Clifford walls preserve that
dimension, even when they mix cells and helpers. Since m>p+2 for p>=3,
this family cannot meet the exact weak-probe contract. I_q describes E's
spectator factor; Q itself need not be identity on those same physical wires.
The p=2 finite construction is consistent with the equality bound four.

## Robust low-weight argument

Assume complete probe error at most epsilon. For v in K_Z of weight w,
unitary telescoping gives, for every y,

    ||Q_y P_v Q_y†-P_(Pi_y v)|| <= 2w epsilon.             (1)

This does not assert that Pi_y v belongs to K_Z. The second operator in (1)
need not be Pauli. It is a commuting-Z pullback, which is sufficient below.

All query conjugates of the source Pauli P_v commute pairwise: for the
field shear, P commutes with every Q_x P Q_x†, and additivity reduces each
pair to this case. Their label sum on a four-address affine plane is zero,
so its ordered product is +I or -I, never ±iI.

Use a plane based at address zero. Its base term is exactly P_v; each of the
other three terms has error at most 2w epsilon in (1). The target product is
E† Z^(Delta_y Delta_z v) E. A nontrivial such Z word has distance two from
either scalar +I or -I. Hence

    epsilon < 1/(3w) => every second difference of v vanishes. (2)

So v is affine. A nonzero translation stabilizer then makes a field Pauli
and its conjugate lie within 2w epsilon<2; they commute, so the Pauli gap
forces equality. The common nonzero fixed kernel makes its label fixed under
every query. Applying (1) again, the target Z-word gap forces all translates
of v to equal v, including the literal sign. Therefore v is constant.

This proof never assumes approximate invariance of K_Z. That is the key
reason it avoids the gap in a naive robust version of exact covariance.

## A uniform weight-four improvement

For the independent-D-cell family, assume epsilon<1/4 so the image X-mask
space has dimension at most p+1. The domain dimension m+2q implies

    dim K_Z >= m-p-1,
    |F2^m/K_Z| <= 2^(p+1)=2m.

There are 1+m+binom(m,2) binary vectors of weight at most two. For m>2,

    1+m+binom(m,2)-2m = (m-1)(m-2)/2 > 0.

Two distinct such vectors therefore lie in the same quotient coset. Their
sum is a nonzero v in K_Z of weight at most four. This is a direct
pigeonhole argument; no coding-theory theorem is required.

For p>=3, four is strictly smaller than m, so this v is nonconstant. If
epsilon<1/12, then epsilon<1/(3w), contradicting (2). Thus the independent
four-wire-cell class obeys the uniform complete-probe error floor

    epsilon >= 1/12       for every p>=3.

The previously proved symplectic-rank argument remains stronger at p>=5,
where it gives epsilon>=1/4. The uniform 1/12 bound is useful already at
p=3 and p=4. It remains a restriction on independent cells with Clifford
outer walls, not on coupled non-Clifford bands or unrestricted decoders.
