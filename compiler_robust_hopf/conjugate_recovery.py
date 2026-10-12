"""Full-sector source diagnostics and a native recovery counterexample.

Matrices below are bounded, floating diagnostics, not an arbitrary-precision
compiler. The counterexample and gate lifts are literal chronological native
words, with qubit zero the low integer bit. Fine sources are not emitted.
"""
from __future__ import annotations

from numbers import Integral
import numpy as np
from .native_ordered_banks import inverse

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def _square(a):
    a = np.asarray(a, dtype=complex)
    if a.ndim != 2 or not a.shape[0] or a.shape[0] != a.shape[1] or not np.isfinite(a).all():
        raise ValueError('A nonempty finite square matrix is required.')
    return a


def conjugate_pair(a):
    """diag(A, conjugate(A)); not A tensored with signal identity."""
    a = _square(a)
    z = np.zeros_like(a)
    return np.block([[a, z], [z, a.conj()]])


def recover(u, accepted=None):
    """Actual -U Z U-dagger Z U; accepted coordinates are first."""
    u = _square(u)
    k = len(u) // 2 if accepted is None else accepted
    if isinstance(k, bool) or not isinstance(k, Integral) or not 0 < k < len(u):
        raise ValueError('Accepted dimension must be a proper positive integer.')
    z = np.diag([1] * int(k) + [-1] * (len(u) - int(k)))
    return -u @ z @ u.conj().T @ z @ u


def recovery_blocks(u):
    """Four-block polynomial for equal-sized sectors of a unitary U."""
    u = _square(u)
    if len(u) % 2:
        raise ValueError('Equal sectors require even matrix dimension.')
    k = len(u) // 2
    a, b, c, d = u[:k, :k], u[:k, k:], u[k:, :k], u[k:, k:]
    eye = np.eye(k)
    return np.block([[3*a-4*a@a.conj().T@a, b@(eye-4*d.conj().T@d)],
                     [c@(eye-4*a.conj().T@a), 3*d-4*d@d.conj().T@d]])


def two_sector_source(e, f):
    """Exact-normalization ideal source; neither E nor F is a free gate."""
    e, f = _square(e), _square(f)
    if e.shape != f.shape:
        raise ValueError('Sector targets must have matching dimensions.')
    eye = np.eye(len(e))
    for a in (e, f):
        if not np.allclose(a.conj().T @ a, eye, atol=1e-10, rtol=0):
            raise ValueError('Sector targets must be unitary.')
    return np.block([[e/2, -np.sqrt(3)*eye/2], [np.sqrt(3)*f@e/2, f/2]])


def _strict(word):
    out = []
    for op, *bits in word:
        if op == 'x':
            out.extend([('h', bits[0]), ('s', bits[0]), ('s', bits[0]), ('h', bits[0])])
        elif op == 'z':
            out.extend([('s', bits[0]), ('s', bits[0])])
        elif op == 'cz':
            out.extend([('h', bits[1]), ('cx', *bits), ('h', bits[1])])
        else:
            out.append((op, *bits))
    return out


def lift_native_word(word, signal):
    """Emit diag(U, conjugate(U)) using a distinct arbitrary signal bit."""
    if isinstance(signal, bool) or not isinstance(signal, Integral) or signal < 0:
        raise ValueError('Signal must be a nonnegative integer wire.')
    out = []
    for gate in word:
        op, *bits = gate
        arity = 2 if op == 'cx' else 1
        if op not in ('h', 's', 'sdg', 't', 'tdg', 'x', 'cx') or len(bits) != arity:
            raise ValueError('Unsupported native gate or arity.')
        if len(set(bits)) != len(bits) or signal in bits or any(
                isinstance(b, bool) or not isinstance(b, Integral) or b < 0 for b in bits):
            raise ValueError('Data wires must be distinct from the signal.')
        if op in ('t', 'tdg'):
            pair = [('cx', bits[0], signal), ('t', signal),
                    ('cx', bits[0], signal), ('tdg', signal)]
            out.extend(pair if op == 't' else inverse(pair))
        elif op in ('s', 'sdg'):
            out.extend([(op, bits[0])] + _strict([('cz', signal, bits[0])]))
        else:
            out.append((op, *bits))
    return out


def recovery_counterexample_words():
    """Source block I/2 recovers to CZ, not identity on an occupied flag.

    Wires: 3 recovery, 2 logical, 1 and 0 dirty core. Returns strict native
    source and recovery words; no large generated JSON is stored in the repo.
    """
    f, x, b0, b1 = 3, 2, 1, 0
    d1 = [('cx', b1, b0), ('h', b0), ('h', b1), ('s', b0)]
    d2 = [('cx', b1, b0), ('h', b0), ('s', b0)]
    basis = inverse(d1)+[('tdg', b0)]+d1+inverse(d2)+[('tdg', b0)]+d2
    nword = inverse(basis)+[('x', b0)]+basis
    scalar = [('h', f), ('cx', f, b0)]+nword+[
        ('x', f), ('cx', f, b0), ('x', f), ('h', f)]
    source = scalar+[('cz', f, x)]
    recovery = source+[('z', f)]+inverse(source)+[('z', f)]+source+[
        ('x', f), ('z', f), ('x', f), ('z', f)]
    return _strict(source), _strict(recovery)


def _kron(factors):
    out = np.ones((1, 1), dtype=complex)
    for factor in factors:
        out = np.kron(out, factor)
    return out


def _mean_signs(value, q):
    weights = np.array([2.**(-j-1) for j in range(q-1)]+[2.**(1-q)])
    grid = 2.**(1-q)
    left = int(np.clip(np.floor((1-np.clip(value, -1, 1))/(2*grid)+.5), 0, 2**(q-1)))
    signs = np.ones(q)
    for j, weight in enumerate(weights):
        cost = int(round(weight/grid))
        if left >= cost:
            signs[j] = -1
            left -= cost
    return signs


def finite_phase_word(phi, q=5):
    """Full finite Majorana word, for 5 <= q <= 7 diagnostics only.

    Floating sign selection tests identities, not the certified digit rule.
    Returns Q, its five-call recovery V, Hermitian rejection B, and z.
    """
    if isinstance(q, bool) or not isinstance(q, Integral) or not 5 <= q <= 7:
        raise ValueError('This dense diagnostic supports integer q=5,6,7 only.')
    if not np.isfinite(phi) or np.iscomplexobj(phi):
        raise ValueError('Phase must be finite and real.')
    q = int(q)
    m = q+1
    gammas = [_kron([Z]*j+[p]+[I2]*(m-j-1)) for j in range(m) for p in (X, Y)]
    weights = np.array([2.**(-j-1) for j in range(q-1)]+[2.**(1-q)])
    coeff = np.zeros(2*m)
    coeff[0] = .5
    for j in range(q):
        coeff[2*j+1], coeff[2*j+2] = np.sqrt(weights[j]/2), np.sqrt(weights[j]/4)
    fixed = np.ones(2*m)
    fixed[2:2*q+1:2] = -1
    beta = (np.sqrt(5)-1)/4
    ss, tt = beta*np.cos(phi), beta*np.sin(phi)
    plus, minus = .75*ss+.5*tt, .25*ss-.5*tt
    head = 1 if plus >= 0 else -1
    signs = np.ones(2*m)
    signs[0] = head
    signs[1:2*q:2] = _mean_signs(2*plus-head/2, q)
    signs[2:2*q+1:2] = _mean_signs(4*minus, q)
    source = sum(c*g for c, g in zip(coeff, gammas))
    nf = sum(c*f*g for c, f, g in zip(coeff, fixed, gammas))
    ng = sum(c*s*g for c, s, g in zip(coeff, signs, gammas))
    eye, zero = np.eye(len(source)), np.zeros_like(source)
    had = np.kron(H, eye)
    c0 = np.block([[source, zero], [zero, eye]])
    c1 = np.block([[eye, zero], [zero, source]])
    def scalar(n):
        return had @ c0 @ np.kron(I2, n) @ c1 @ had
    a = scalar(nf)
    sg = np.kron(np.diag([1, 1j]), eye)
    block = a.conj().T @ sg @ scalar(ng) @ sg.conj().T @ a
    zflag = np.kron(Z, eye)
    v = block@zflag@block.conj().T@zflag@block@zflag@block.conj().T@zflag@block
    s = float((coeff**2) @ signs)
    tau = float((coeff**2*fixed) @ signs)
    df, dg = (source@nf-nf@source)/2, (source@ng-ng@source)/2
    b = -1j*(dg/4+df@dg@df)
    return block, v, b, s+1j*(tau-s/2)
