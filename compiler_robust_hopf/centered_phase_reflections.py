"""Exact-word generator and bounded diagnostics for one center-phase family.

The two diagonal occurrences are a table and its actual inverse. Matrix
functions use floating arithmetic and are not certified arbitrary-precision
compilers. Only the explicitly specified five-qubit member is emitted.
"""
from __future__ import annotations

from numbers import Integral, Real
import numpy as np

from .native_coarse_fixture import Word, inverse_word
from .paid_logical_encodings import expand_native


def _dimension(size: int) -> int:
    if isinstance(size, bool) or not isinstance(size, Integral) or not 2 <= size <= 256 or size & (size - 1):
        raise ValueError('This dense diagnostic requires a dyadic dimension from 2 to 256.')
    return int(size)


def _angle(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError('A finite real angle is required.')
    return float(value)


def _unit(vector, *, real=False):
    v = np.asarray(vector, dtype=complex)
    if v.ndim != 1:
        raise ValueError('A vector is required.')
    _dimension(len(v))
    if not np.isfinite(v).all() or not np.isclose(np.vdot(v, v).real, 1, atol=1e-10, rtol=0):
        raise ValueError('A finite unit vector is required.')
    if real and np.max(np.abs(v.imag)) > 1e-12:
        raise ValueError('A real vector is required.')
    return v.real if real else v


def walsh_matrix(size: int):
    size = _dimension(size)
    return np.array([[(-1.0) ** ((x & y).bit_count()) for y in range(size)]
                     for x in range(size)]) / np.sqrt(size)


def reflection(vector):
    v = _unit(vector)
    return np.eye(len(v)) - 2 * np.outer(v, v.conj())


def centered_word(phases):
    d = np.asarray(phases, dtype=complex)
    if d.ndim != 1:
        raise ValueError('A phase vector is required.')
    _dimension(len(d))
    if not np.isfinite(d).all() or not np.allclose(abs(d), 1, atol=1e-10, rtol=0):
        raise ValueError('All table entries must have unit modulus.')
    f = walsh_matrix(len(d))
    s = np.eye(len(d), dtype=complex)
    s[0, 0] = 1j
    r = np.eye(len(d)); r[0, 0] = -1
    q = (f * d) @ f
    return s @ q @ r @ q.conj().T @ r @ s.conj().T


def optimal_phase_table(normal):
    """A numerical attaining table and the exact-formula minimum (rounded).

    Zero transformed coordinates receive phase one. This convention is a
    finite floating diagnostic, not an exact-zero oracle for input reals.
    """
    v = _unit(normal)
    y = v.copy(); y[0] *= -1j
    y = walsh_matrix(len(v)) @ y
    d = np.ones(len(v), dtype=complex)
    nonzero = np.abs(y) > 0
    d[nonzero] = y[nonzero] / np.abs(y[nonzero])
    # For a unit y, 1-||y||_1**2/N is the squared centered-magnitude norm.
    # Evaluate that variance directly to avoid cancellation near exact coverage.
    magnitudes = abs(y)
    floor = 2 * np.linalg.norm(magnitudes - np.mean(magnitudes))
    return d, float(floor)


def flat_walsh_direction(signs):
    sigma = np.asarray(signs)
    if sigma.ndim != 1:
        raise ValueError('A sign table is required.')
    _dimension(len(sigma))
    if not np.all(np.isin(sigma, (-1, 1))) or np.sum(sigma) != 0:
        raise ValueError('The sign table must contain equally many +1 and -1 entries.')
    return walsh_matrix(len(sigma)) @ sigma.astype(float) / np.sqrt(len(sigma))


def plane_rotation(direction, angle: float):
    b, angle = _unit(direction, real=True), _angle(angle)
    if abs(b[0]) > 1e-12:
        raise ValueError('The direction must be orthogonal to the center.')
    e = np.eye(len(b))[0]
    v = np.cos(angle / 2) * e + np.sin(angle / 2) * b
    return reflection(v) @ reflection(e)


def two_coordinate_floor(angle: float) -> float:
    """Stable form of sqrt(2-2 sqrt(1-sin(angle/2)**4))."""
    angle = _angle(angle)
    t = np.sin(angle / 2) ** 2
    return float(np.sqrt(2) * t / np.sqrt(1 + np.sqrt(max(0.0, 1 - t*t))))


def ordered_hopf_witness(angle: float, size: int = 8):
    """G_(0,2) first and G_(0,1) second, identity on every other label."""
    size, angle = _dimension(size), _angle(angle)
    if size < 4:
        raise ValueError('The witness requires at least four modes.')
    c, s = np.cos(angle), np.sin(angle)
    u = np.eye(size)
    for target in (2, 1):
        u[[0, target], :] = np.array([[c, -s], [s, c]]) @ u[[0, target], :]
    return u


def ordered_hopf_floor(angle: float) -> float:
    angle = _angle(angle)
    if abs(angle) > np.pi / 2:
        raise ValueError('The stated sharp floor assumes |angle| <= pi/2.')
    return float(2 * np.sin(angle / 2)**2)


def nearest_center_reflection(angle: float, size: int = 8):
    """Attainer for the ordered witness's larger single-reflection family."""
    angle = _angle(angle)
    ordered_hopf_floor(angle)
    size = _dimension(size)
    if size < 4:
        raise ValueError('The witness requires at least four modes.')
    h = np.zeros(size); h[0] = 1
    h[1] = h[2] = np.tan(angle / 2)
    h /= np.linalg.norm(h)
    e = np.eye(size)[0]
    return reflection(h) @ reflection(e)


def boolean_direction_words() -> tuple[Word, Word]:
    """Return macro and strict native five-qubit words, at angle pi/2.

    Data are wires 0,1,2; clean work is 3,4. The target direction is
    (0,1,0,1,0,1,0,-1)/2. Work is returned on all eight logical inputs.
    No claim is made about an arbitrary-angle fine emitter.
    """
    a, b = 3, 4
    def center_phase(gate):
        return (tuple(('X', (j,)) for j in range(3))
                + (('CCX', (0, 1, a)), ('CCX', (a, 2, b)), (gate, (b,)),
                   ('CCX', (a, 2, b)), ('CCX', (0, 1, a)))
                + tuple(('X', (j,)) for j in range(3)))
    query = (('CX', (0, a)), ('CCX', (1, 2, a)))
    rotation = (('TDG', (a,)), ('X', (a,)), ('T', (a,)), ('X', (a,)))
    had = tuple(('H', (j,)) for j in range(3))
    q = had + query + rotation + inverse_word(query) + had
    macro = (center_phase('SDG') + center_phase('Z') + inverse_word(q)
             + center_phase('Z') + q + center_phase('S'))
    return macro, expand_native(macro)
