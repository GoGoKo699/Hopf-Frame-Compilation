"""Bounded matrix diagnostics for programmed reflections, not a fine compiler.

Floating support and rank decisions are diagnostic only. The proof chapter
uses exact equations and certified finite search, not these tolerances, for
statements about computable real input. No target-dependent basis is a free gate.
"""
from __future__ import annotations

from collections import deque
from numbers import Real
import numpy as np


def _square(value):
    a = np.asarray(value, dtype=complex)
    if (a.ndim != 2 or a.shape[0] != a.shape[1] or not 2 <= len(a) <= 256
            or not np.isfinite(a).all()):
        raise ValueError('A finite square matrix of dimension 2 through 256 is required.')
    return a


def _unitary(value):
    a = _square(value)
    if not np.allclose(a.conj().T @ a, np.eye(len(a)), atol=2e-10, rtol=0):
        raise ValueError('The diagnostic requires a unitary matrix.')
    return a


def _involution(value):
    a = _unitary(value)
    if not np.allclose(a, a.conj().T, atol=2e-10, rtol=0):
        raise ValueError('A Hermitian unitary is required.')
    return a


def _phases(value, size):
    d = np.asarray(value, dtype=complex)
    if (d.shape != (size,) or not np.isfinite(d).all()
            or not np.allclose(abs(d), 1, atol=2e-10, rtol=0)):
        raise ValueError('A unit-modulus vector of the matching dimension is required.')
    return d


def _same(*matrices):
    if len({a.shape for a in matrices}) != 1:
        raise ValueError('Matrix dimensions must agree.')


def phase_reflection(mixer, center, phases):
    """Return C D (C-dagger J C) D-dagger C-dagger on all columns."""
    c, j = _unitary(mixer), _involution(center)
    _same(c, j)
    d = _phases(phases, len(c))
    a = (c*d) @ c.conj().T
    return a @ j @ a.conj().T


def programmed_pair(mixer, center, phases):
    j = _involution(center)
    return phase_reflection(mixer, j, phases) @ j


def recover_phase_gauge(base, target, *, tolerance=1e-9):
    """Recover D in target=D base D-dagger, or return None diagnostically.

    Both directed entries are checked. A small numerical edge is ignored
    only during tree construction; the final full matrix residual is still
    checked. None is not an exact nonmembership certificate for real inputs.
    """
    b, k = _unitary(base), _unitary(target)
    _same(b, k)
    if (isinstance(tolerance, bool) or not isinstance(tolerance, Real)
            or not np.isfinite(tolerance) or not 0 < tolerance < 1e-3):
        raise ValueError('Use a positive diagnostic tolerance below 1e-3.')
    if (np.max(abs(np.diag(k)-np.diag(b))) > tolerance
            or np.max(abs(abs(k)-abs(b))) > tolerance):
        return None
    size = len(b)
    adjacency = (abs(b) > tolerance) | (abs(b.T) > tolerance)
    phases = np.ones(size, dtype=complex)
    visited = np.zeros(size, dtype=bool)
    for root in range(size):
        if visited[root]:
            continue
        visited[root] = True
        queue = deque([root])
        while queue:
            i = queue.popleft()
            for j in np.flatnonzero(adjacency[i]):
                if visited[j] or i == j:
                    continue
                if abs(b[i, j]) > tolerance:
                    gain = k[i, j]/b[i, j]
                else:
                    gain = np.conj(k[j, i]/b[j, i])
                if abs(gain) == 0:
                    return None
                phases[j] = phases[i]*np.conj(gain/abs(gain))
                visited[j] = True
                queue.append(int(j))
    reconstructed = phases[:, None]*b*phases.conj()[None, :]
    if np.linalg.norm(k-reconstructed, 2) > 4*len(b)*tolerance:
        return None
    return phases


def cycle_discrepancy(base, target, cycle):
    """A rounded evaluation of the exact product-difference/p error bound."""
    b, k = _unitary(base), _unitary(target)
    _same(b, k)
    vertices = tuple(cycle)
    if (len(vertices) < 2 or len(set(vertices)) != len(vertices)
            or any(isinstance(i, bool) or not isinstance(i, (int, np.integer))
                   or not 0 <= i < len(b) for i in vertices)):
        raise ValueError('Use at least two distinct valid cycle vertices.')
    edges = tuple(zip(vertices, vertices[1:]+vertices[:1]))
    return float(abs(np.prod([k[i,j] for i,j in edges])
                     - np.prod([b[i,j] for i,j in edges]))/len(edges))


def nearest_reversible(target, center):
    """Unrestricted matrix optimum near identity; no native matrix-sign claim."""
    w, j = _unitary(target), _involution(center)
    _same(w, j)
    if np.linalg.norm(w-np.eye(len(w)), 2) >= 1:
        raise ValueError('The inertia argument assumes ||W-I|| < 1.')
    q = w @ j
    herm = (q+q.conj().T)/2
    values, vectors = np.linalg.eigh(herm)
    sign = (vectors*np.where(values > 0, 1., -1.)) @ vectors.conj().T
    omega = min(1., float(np.linalg.norm((q-q.conj().T)/(2j), 2)))
    floor = omega*np.sqrt(2/(1+np.sqrt(max(0., 1-omega*omega))))
    return sign @ j, float(floor)


def shared_relative_kernel(mixer, center, target):
    c, j, w = _unitary(mixer), _involution(center), _unitary(target)
    _same(c, j, w)
    b = c.conj().T @ j @ c
    m = c.conj().T @ w @ c
    return abs(b)**2-np.diag(np.diag(m))


def shared_pair_floor(mixer, center, target):
    """Numerical smallest singular value; certify it separately for a proof."""
    return float(np.linalg.svd(shared_relative_kernel(mixer, center, target),
                               compute_uv=False)[-1])


def recover_shared_pair(mixer, center, target, *, tolerance=1e-9):
    """Finite diagnostic for a simple relative kernel, not general coverage.

    Returns (u,v) with target=R(v)R(u), or None. Multiple numerical null
    directions are deliberately left unresolved instead of choosing one.
    """
    c, j, w = _unitary(mixer), _involution(center), _unitary(target)
    _same(c, j, w)
    if (isinstance(tolerance, bool) or not isinstance(tolerance, Real)
            or not np.isfinite(tolerance) or not 0 < tolerance < 1e-3):
        raise ValueError('Use a positive diagnostic tolerance below 1e-3.')
    matrix = shared_relative_kernel(c, j, w)
    _, singular, vh = np.linalg.svd(matrix)
    if np.count_nonzero(singular <= tolerance) != 1:
        return None
    d = vh[-1].conj()*np.sqrt(len(c))
    if np.max(abs(abs(d)-1)) > 100*tolerance:
        return None
    d /= abs(d)
    b = c.conj().T @ j @ c
    m = c.conj().T @ w @ c
    a = (b*d) @ b @ np.diag(d.conj())
    v = recover_phase_gauge(a, m, tolerance=tolerance)
    if v is None:
        return None
    u = v*d
    if np.linalg.norm(phase_reflection(c,j,v) @ phase_reflection(c,j,u)-w, 2) > 8*len(c)*tolerance:
        return None
    return u, v


def hermitian_coordinates(size):
    """Sparse basis descriptors: diagonal, then paired real/imaginary entries."""
    if isinstance(size, bool) or not isinstance(size, (int, np.integer)) or not 2 <= size <= 32:
        raise ValueError('Hermitian systems are limited to dimensions 2 through 32.')
    entries = [(i, i, 'diagonal') for i in range(size)]
    for i in range(size):
        for j in range(i+1, size):
            entries += [(i,j,'real'), (i,j,'imaginary')]
    return tuple(entries)


def reversing_moment_system(target, first_mixer, second_mixer, center):
    """Real linear necessary conditions, including all complex Hermitian H.

    Column order is hermitian_coordinates. Rows are real strict-upper
    reversing entries, imaginary upper entries, then both diagonal lists.
    Involutivity, signature and complete phase-orbit membership are additional.
    """
    w, c1, c2, j = (_unitary(target), _unitary(first_mixer),
                     _unitary(second_mixer), _involution(center))
    _same(w,c1,c2,j)
    if np.max(abs(w.imag)) > 1e-12:
        raise ValueError('This nonspectral system assumes a real target.')
    n = len(w)
    cols = []
    for a,b,kind in hermitian_coordinates(n):
        h = np.zeros((n,n),complex)
        if kind == 'diagonal':
            h[a,a] = 1
        elif kind == 'real':
            h[a,b] = h[b,a] = 1
        else:
            h[a,b], h[b,a] = 1j, -1j
        defect = w@h-h@w.T
        cols.append([defect[a,b].real for a in range(n) for b in range(a+1,n)]
                    + [defect[a,b].imag for a in range(n) for b in range(a,n)]
                    + list(np.diag(c1.conj().T@h@c1).real)
                    + list(np.diag(c2.conj().T@w@h@c2).real))
    rhs = np.r_[np.zeros(n*n), np.diag(c1.conj().T@j@c1).real,
                 np.diag(c2.conj().T@j@c2).real]
    return np.array(cols).T, rhs
