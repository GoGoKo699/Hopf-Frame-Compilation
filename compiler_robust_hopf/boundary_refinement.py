"""Finite matrix diagnostics for rank-weighted blocks and canonical forests.

These routines do not emit a fine compiler, decide ranks of computable inputs,
or certify endpoint errors. Floating zero/rank decisions are diagnostics only;
proof-level algorithms use declared support/rank bounds and certified search.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import numpy as np
from .conjugate_recovery import _square


def reflection(vector):
    v = np.asarray(vector, dtype=complex)
    if v.ndim != 1 or not len(v) or not np.isfinite(v).all() or np.linalg.norm(v) == 0:
        raise ValueError('A finite nonzero vector is required.')
    v = v/np.linalg.norm(v)
    return np.eye(len(v))-2*np.outer(v, v.conj())


def real_state_extension(vector, anchor=0):
    """Full orthogonal extension taking a center to the supplied unit state."""
    v = np.asarray(vector)
    if v.ndim != 1 or np.iscomplexobj(v) or not np.isfinite(v).all():
        raise ValueError('A finite real vector is required.')
    if isinstance(anchor, bool) or not isinstance(anchor, Integral) or not 0 <= anchor < len(v):
        raise ValueError('Anchor is outside the vector.')
    if not np.isclose(v@v, 1, atol=1e-10, rtol=0):
        raise ValueError('State must be normalized.')
    e = np.eye(len(v))[:, int(anchor)]
    sign = 1 if v[anchor] >= 0 else -1
    center = np.eye(len(v))
    center[anchor, anchor] = -sign
    return np.real_if_close(reflection(e+sign*v)) @ center


def commutator(u, v):
    u, v = _square(u), _square(v)
    if u.shape != v.shape:
        raise ValueError('Commutator dimensions differ.')
    return u@v@u.conj().T@v.conj().T


def spectral_reflection_lcu(frame, phases):
    """Ideal complete LCU word with coefficient norm exactly two.

    Branch zero is the original m-dimensional input. Mode dimension is
    (rank+3)*m; this matrix construction is not a claimed native allocation.
    Returns the word, its target, and the absolute coefficient sum before pad.
    """
    frame = np.asarray(frame, dtype=complex)
    phases = np.asarray(phases, dtype=float)
    if frame.ndim != 2 or phases.ndim != 1 or frame.shape[1] != len(phases):
        raise ValueError('Frame/phase shape mismatch.')
    m, r = frame.shape
    if not m or not np.isfinite(frame).all() or not np.isfinite(phases).all():
        raise ValueError('Finite nonempty frame required.')
    if not np.allclose(frame.conj().T@frame, np.eye(r), atol=1e-10, rtol=0):
        raise ValueError('Frame columns must be orthonormal.')
    z = np.exp(1j*phases)-1
    coeff = np.r_[1+np.sum(z)/2, -z/2]
    lam = float(np.sum(np.abs(coeff)))
    if lam >= 2:
        raise ValueError('Coefficient norm must be strictly below two.')
    coeff = np.r_[coeff, (2-lam)/2, -(2-lam)/2]
    weights = np.sqrt(np.abs(coeff)/2)
    prep = np.kron(real_state_extension(weights), np.eye(m))
    select = np.eye((r+3)*m, dtype=complex)
    for a in range(r):
        start = (a+1)*m
        select[start:start+m, start:start+m] = reflection(frame[:, a])
    phase = np.ones(r+3, dtype=complex)
    nonzero = np.abs(coeff) > 0
    phase[nonzero] = coeff[nonzero]/np.abs(coeff[nonzero])
    word = prep.conj().T @ np.kron(np.diag(phase), np.eye(m)) @ select @ prep
    target = np.eye(m)+(frame*z)@frame.conj().T
    return word, target, lam


def rank_mode_bins(n_modes, sizes, ranks):
    """Exact small-support chunk/bin schedule; no quantum circuit is implied."""
    if isinstance(n_modes, bool) or not isinstance(n_modes, Integral) or n_modes < 256 or n_modes & (n_modes-1):
        raise ValueError('Use a power-of-two ambient size at least 256.')
    sizes, ranks = tuple(sizes), tuple(ranks)
    if len(sizes) != len(ranks):
        raise ValueError('Mismatched block lists.')
    chunks = []
    for j, (m, r) in enumerate(zip(sizes, ranks)):
        if any(isinstance(x, bool) or not isinstance(x, Integral) for x in (m, r)) or not 2 <= r <= m <= n_modes//128:
            raise ValueError('Only small supports with 2 <= rank <= size are accepted.')
        limit = n_modes//(8*m)-3
        q = min(2*r, m)
        chunks.append([(j, m, min(limit, q-a)) for a in range(0, q, limit)])
    bins = []
    for t in range(max(map(len, chunks), default=0)):
        current, weight = [], 0
        for rows in chunks:
            if t >= len(rows):
                continue
            row = rows[t]
            cost = row[1]*(row[2]+3)
            if current and weight+cost > n_modes//4:
                bins.append(tuple(current))
                current, weight = [], 0
            current.append(row)
            weight += cost
        if current:
            bins.append(tuple(current))
    return tuple(bins)


def forest_matching(parents):
    """Maximum matching of a topologically numbered declared forest."""
    parents = tuple(parents)
    for i, p in enumerate(parents):
        if isinstance(p, bool) or not isinstance(p, Integral) or not -1 <= p < i:
            raise ValueError('Each parent must be -1 or an earlier vertex.')
    used = set()
    edges = []
    for i in reversed(range(len(parents))):
        p = parents[i]
        if p >= 0 and i not in used and p not in used:
            edges.append((p, i))
            used.update((p, i))
    return tuple(edges)


def _real_shift(a):
    a = np.asarray(a)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or np.iscomplexobj(a) or not np.isfinite(a).all():
        raise ValueError('A finite real square shift is required.')
    if np.any(np.triu(a)):
        raise ValueError('Use a strictly outward, topologically ordered shift.')
    return a.astype(float)


def canonical_frame(shift):
    """Classical Cayley evaluation; the inverse is not a native primitive."""
    shift = _real_shift(shift)
    eye = np.eye(len(shift))
    if not len(shift):
        return eye
    r = np.linalg.solve(eye-shift, eye)
    x = r.T-r
    return np.linalg.solve((eye+x).T, (eye-x).T).T


def canonical_star_word(shift):
    """Independent local-star evaluation in topological order."""
    shift = _real_shift(shift)
    if np.any(np.count_nonzero(shift, axis=1) > 1):
        raise ValueError('A forest row has at most one parent.')
    word = np.eye(len(shift))
    for v in range(len(shift)):
        children = np.flatnonzero(shift[:, v])
        if not len(children):
            continue
        w = shift[children, v]
        size = len(w)+1
        square = float(w@w)
        block = np.eye(size)
        block[0, 0] = (1-square)/(1+square)
        block[1:, 0] = 2*w/(1+square)
        block[0, 1:] = -block[1:, 0]
        block[1:, 1:] -= 2*np.outer(w, w)/(1+square)
        labels = np.r_[v, children]
        word[labels, :] = block @ word[labels, :]
    return word


@dataclass(frozen=True)
class LeafDeflation:
    matching: np.ndarray
    basis: np.ndarray
    compact_shift: np.ndarray
    embedding: np.ndarray
    intermediate_shift: np.ndarray
    parents: tuple[int, ...]


def leaf_deflation(a, d, declared=None):
    """One exact-formula leaf-pair reduction evaluated in floating arithmetic."""
    a = _real_shift(a)
    d = np.asarray(d)
    m = len(a)
    if d.shape != (m,) or np.iscomplexobj(d) or not np.isfinite(d).all():
        raise ValueError('Matching weights must be a finite real vector.')
    mask = a != 0 if declared is None else np.asarray(declared, dtype=bool)
    if mask.shape != a.shape or np.any(np.triu(mask)) or np.any(mask.sum(axis=1) > 1) or np.any(a[~mask]):
        raise ValueError('Declared support must contain A and have unique parents.')
    parents = tuple(map(int, np.flatnonzero(np.any(mask, axis=0))))
    co = 1/np.sqrt(1+d*d)
    si = d*co
    g = np.block([[np.diag(co), np.diag(si)], [-np.diag(si), np.diag(co)]])
    b = co[:, None]*a*co[None, :]
    lt = np.block([[b, np.zeros_like(a)], [d[:, None]*b, np.zeros_like(a)]])
    basis = np.eye(2*m)
    active, gammas = [], []
    for v in parents:
        children = np.flatnonzero(mask[:, v])
        indices = [int(i) for i in children if i not in parents]+[m+int(i) for i in children]
        z = m+int(children[0])
        weights = np.r_[co*a[:, v], si*a[:, v]]
        weights[list(parents)] = 0
        gamma = float(np.linalg.norm(weights))
        if gamma:
            anchor = indices.index(z)
            local = real_state_extension(weights[indices]/gamma, anchor)
            basis[np.ix_(indices, indices)] = local
        active.extend((v, z))
        gammas.append(gamma)
    compact = np.zeros((2*len(parents), 2*len(parents)))
    for j, v in enumerate(parents):
        compact[2*j+1, 2*j] = gammas[j]*co[v]
        for k, u in enumerate(parents):
            compact[2*k, 2*j] = b[u, v]
    embedding = np.eye(2*m)[:, active]
    return LeafDeflation(g, basis, compact, embedding, lt, parents)


def recursive_leaf_frame(shift):
    """Full diagnostic reconstruction for the original binomial topology."""
    shift = _real_shift(shift)
    n = len(shift)
    if n < 1 or n & (n-1):
        raise ValueError('A positive dyadic mode dimension is required.')
    mask = np.zeros_like(shift, dtype=bool)
    for j in range(1, n):
        mask[j, j-(j & -j)] = True
    if np.any(shift[~mask]):
        raise ValueError('Shift must use the original computational tree.')
    if n == 1:
        return np.eye(1)
    order = np.r_[np.arange(0, n, 2), np.arange(1, n, 2)]
    packed = shift[np.ix_(order, order)]
    m = n//2
    out = leaf_deflation(packed[:m, :m], np.diag(packed[m:, :m]), mask[np.ix_(order[:m], order[:m])])
    core = recursive_leaf_frame(out.compact_shift) if len(out.compact_shift) else np.zeros((0, 0))
    full = np.eye(n)+out.embedding@(core-np.eye(len(core)))@out.embedding.T
    result = out.matching.T@out.basis@full@out.basis.T@out.matching.T
    restored = np.zeros_like(result)
    restored[np.ix_(order, order)] = result
    return restored
