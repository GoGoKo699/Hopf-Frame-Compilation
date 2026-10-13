"""Signed-cycle native words and bounded spectral-chart diagnostics.

Only fixed circuits are emitted. Eigenvalue, polar, and quadrature routines
use floating arithmetic and are not certified spectral or native oracles.
Gate tuples reuse the repository's chronological little-endian convention.
"""
from __future__ import annotations

from numbers import Integral
import numpy as np

from .native_coarse_fixture import Word, inverse_word
from .paid_logical_encodings import expand_native
from .programmed_reflections import _unitary


def _integer(value, lower, upper, name):
    if isinstance(value, bool) or not isinstance(value, Integral) or not lower <= value <= upper:
        raise ValueError(f'{name} must be an integer from {lower} through {upper}.')
    return int(value)


def signed_cycle(size):
    size = _integer(size, 2, 256, 'matrix dimension')
    if size % 2:
        raise ValueError('The signed cycle requires an even dimension.')
    p = np.zeros((size,size))
    p[np.arange(1,size), np.arange(size-1)] = 1
    p[0,-1] = -1
    return p


def complement(size):
    return np.eye(len(signed_cycle(size)))[::-1]


def spectral_chart(residual):
    """Floating realization of the fixed-contour chart on its proved domain."""
    e = _unitary(residual)
    if len(e) % 2 or np.max(abs(e.imag)) > 1e-12:
        raise ValueError('A real orthogonal residual of even dimension is required.')
    e = e.real
    size = len(e)
    p, s = signed_cycle(size), complement(size)
    gap = 2*np.sin(np.pi/size)
    delta = float(np.linalg.norm(e-np.eye(size), 2))
    if delta > gap/(16*np.sqrt(size)):
        raise ValueError('Residual lies outside the proved spectral-chart radius.')
    reference = np.exp(1j*(2*np.arange(size)+1)*np.pi/size)
    basis = reference[None,:]**(-np.arange(size)[:,None])/np.sqrt(size)
    target = e@p
    values, vectors = np.linalg.eig(target)
    # In the proved domain each fixed circle contains exactly one eigenvalue.
    order = np.argmin(abs(reference[:,None]-values[None,:]),axis=1)
    if len(set(order)) != size or np.max(abs(values[order]-reference)) >= gap/3:
        raise ArithmeticError('Numerical eigenvalue labeling did not respect the fixed circles.')
    vectors = vectors[:,order]
    vectors /= np.linalg.norm(vectors,axis=0)
    q = np.array([np.outer(vectors[:,j], vectors[:,j].conj()) for j in range(size)])
    a = sum((q[j] @ np.outer(basis[:,j],basis[:,j].conj()) for j in range(size)),
            np.zeros((size,size),complex))
    if np.linalg.norm(a.imag,2) > 1e-8:
        raise ArithmeticError('Conjugate pairing lost numerical reality.')
    left, _, right = np.linalg.svd(a.real)
    o = left@right
    h1 = o@s@o.T
    h2 = target@h1
    return {'P':p, 'S':s, 'T':target, 'A':a, 'O':o, 'H1':h1, 'H2':h2,
            'projectors':q, 'reference':reference, 'gap':float(gap), 'delta':delta}


def contour_projector(target, center, radius, nodes):
    """Bounded quadrature only; caller separately proves spectral separation."""
    t = _unitary(target)
    nodes = _integer(nodes, 1, 256, 'quadrature nodes')
    if not np.isfinite(center) or not np.isfinite(radius) or not 0 < radius < 2:
        raise ValueError('A finite center and radius in (0,2) are required.')
    eye = np.eye(len(t),dtype=complex)
    result = np.zeros_like(t)
    for k in range(nodes):
        z = np.exp(2j*np.pi*k/nodes)
        result += radius*z*np.linalg.solve((center+radius*z)*eye-t,eye)/nodes
    return result


def matching_generator(size):
    """The original final-layer Hopf witness, repeated on an arbitrary bit."""
    size = _integer(size, 8, 128, 'original mode dimension')
    if size & (size-1):
        raise ValueError('Use a dyadic original dimension.')
    k = np.zeros((size,size))
    for p in range(size//2):
        a = -np.sin(2*np.pi*(2*p+1)/size)
        k[2*p+1,2*p], k[2*p,2*p+1] = a,-a
    return np.kron(np.eye(2),k)


def centered_cycle_vector(size, index):
    size = len(signed_cycle(size))
    index = _integer(index, 0, size-1, 'spectral index')
    omega = (2*index+1)*np.pi/size
    return np.exp(-1j*omega*(np.arange(size)-(size-1)/2))/np.sqrt(size)


def borrowed_mcx(controls, target, dirty) -> Word:
    """Exact Boolean echo. Helpers may be arbitrary and return exactly."""
    controls = tuple(controls)
    dirty = tuple(dirty)
    bits = controls+(target,)+dirty
    if (any(isinstance(q,bool) or not isinstance(q,Integral) or q < 0 for q in bits)
            or len(set(bits)) != len(bits)):
        raise ValueError('All controls, target and helper wires must be distinct nonnegative integers.')
    k = len(controls)
    if k == 0:
        return (('X',(target,)),)
    if k == 1:
        return (('CX',(controls[0],target)),)
    if k == 2:
        return (('CCX',controls+(target,)),)
    if len(dirty) < k-2:
        raise ValueError('Insufficient arbitrary borrowed chain.')
    ladder = tuple(('CCX',(controls[j],dirty[j-2],dirty[j-1])) for j in range(2,k-1))
    middle = (('CCX',(controls[-1],dirty[k-3],target)),)
    tail = ladder+middle+inverse_word(ladder)
    first = (('CCX',(controls[0],controls[1],dirty[0])),)
    return first+tail+first+tail


def signed_cycle_word(data_bits: int) -> Word:
    """Emit P on low data bits, identity on max(0,data_bits-3) dirty bits."""
    m = _integer(data_bits, 3, 256, 'data width')
    helpers = tuple(range(m,m+max(0,m-3)))
    word = (('H',(m-1,)),)+borrowed_mcx(tuple(range(m-1)),m-1,helpers)+(('H',(m-1,)),)
    for j in range(m-1,0,-1):
        word += borrowed_mcx(tuple(range(j)),j,helpers)
    return expand_native(word+(('X',(0,)),))


def alternative_baseline_word(data_bits: int) -> Word:
    """Exact Clifford P**(M/4) S, no helpers, on all data inputs."""
    m = _integer(data_bits, 3, 256, 'data width')
    low,high = m-2,m-1
    macro = tuple(('X',(j,)) for j in range(m))
    macro += (('CZ',(low,high)),('CX',(low,high)),('X',(low,)))
    return expand_native(macro)


def reflection_mixer_words() -> tuple[Word, Word]:
    """The two exact logical mixers rejected by the rational certificate."""
    seq = (('H',(1,)),('CCX',(2,1,0)),('CCX',(0,1,2)),('H',(1,)),
           ('H',(2,)),('CCX',(1,2,0)),('H',(0,)),('CX',(1,2)),
           ('CCX',(0,2,1)),('H',(2,)))
    other = tuple((g,tuple(2-j for j in qs)) for g,qs in reversed(seq))
    return (expand_native((('H',(0,)),('S',(0,)))+seq),
            expand_native((('H',(1,)),('S',(1,)))+other))


def hopf_corner_word() -> Word:
    """All seven n=3 angles pi/2: W|0>=|7>, W|j>=-|j-1> otherwise."""
    x = lambda q: (('X',(q,)),)
    z = lambda q: (('Z',(q,)),)
    minus = x(0)+z(0)+x(0)+z(0)
    masks = x(0)+x(1)+x(2)
    zero = masks+(('H',(2,)),('CCX',(0,1,2)),('H',(2,)))+masks
    dec = x(0)+x(1)+(('CCX',(0,1,2)),)+x(1)+x(0)+x(0)+(('CX',(0,1)),)+x(0)+x(0)
    return expand_native(minus+zero+dec)
