"""Finite diagnostics for supported dilations with disjoint rejection banks.

The affine mode maps below are exact integer permutations. The other routines
are bounded floating matrix diagnostics, not a native fine-source emitter or
an endpoint compiler. In particular, supplying a dilation matrix does not
charge or construct its precision source. Mode order is (f, e, logical),
most significant first; an arbitrary dirty factor is always last.
"""
from __future__ import annotations

from numbers import Integral, Real

import numpy as np

from .conjugate_recovery import _square


def _depth(n, depth, bank_type=None):
    if any(isinstance(x, bool) or not isinstance(x, Integral)
           for x in (n, depth)) or not 1 <= n <= 12 or not 0 <= depth < n:
        raise ValueError('Use integer 1 <= n <= 12 and 0 <= depth < n.')
    if bank_type is not None and (
            isinstance(bank_type, bool) or not isinstance(bank_type, Integral)
            or not 0 <= bank_type < 3 or depth > n-3):
        raise ValueError('Three-type banks require type 0, 1, 2 and depth <= n-3.')


def active_modes(n, depth):
    """The full rank-2^(depth+1) logical support of one Hopf layer."""
    _depth(n, depth)
    return tuple(row << (n-depth-1) for row in range(1 << (depth+1)))


def rejection_modes(n, depth, bank_type=None):
    """Dedicated mode bank, in the same (prefix, target) order as its support."""
    _depth(n, depth, bank_type)
    base = 2 << n
    kind = 0 if bank_type is None else bank_type << (n-2)
    return tuple(base + ((row & 1) << n) + kind + (1 << depth) + (row >> 1)
                 for row in range(1 << (depth+1)))


def routing_affine(n, depth, bank_type=None):
    """Destination of each non-f input bit, and the output XOR mask.

    Apply this wire permutation and XOR only when f=1. Thus it has an O(n)
    controlled-SWAP/CNOT realization, but this function emits no native word.
    The e bit has index n, and logical bit zero is least significant.
    """
    _depth(n, depth, bank_type)
    suffix = n-depth-1
    destinations = [0]*(n+1)
    for bit in range(suffix):
        destinations[bit] = depth+1+bit
    destinations[suffix] = n
    for bit in range(depth):
        destinations[n-depth+bit] = bit
    destinations[n] = depth
    mask = (1 << depth)
    if bank_type is not None:
        mask |= bank_type << (n-2)
    return tuple(destinations), mask


def routing_permutation(n, depth, bank_type=None):
    """Exact truth table: output basis label for each input basis label."""
    destinations, mask = routing_affine(n, depth, bank_type)
    flag = 1 << (n+1)
    result = np.arange(4 << n, dtype=int)
    for value in range(flag):
        moved = sum(((value >> bit) & 1) << dest
                    for bit, dest in enumerate(destinations))
        result[flag+value] = flag+(moved ^ mask)
    return result


def lift_supported_dilation(block, active, failure, mode_count, dirty_dimension=1):
    """Embed an actual unitary on (active then failure) modes; fix all others.

    Within each mode the dirty index runs fastest. Off-diagonal/rejected
    blocks may act on arbitrary dirty inputs; they are never replaced by
    scalars. Unitarity is checked numerically for this finite diagnostic.
    """
    block = _square(block)
    active, failure = tuple(active), tuple(failure)
    if (isinstance(mode_count, bool) or not isinstance(mode_count, Integral)
            or isinstance(dirty_dimension, bool)
            or not isinstance(dirty_dimension, Integral)
            or mode_count <= 0 or dirty_dimension <= 0
            or not active or len(active) != len(failure)):
        raise ValueError('Positive dimensions and equal nonempty mode lists required.')
    modes = active+failure
    if len(set(modes)) != len(modes) or any(
            isinstance(i, bool) or not isinstance(i, Integral)
            or not 0 <= i < mode_count for i in modes):
        raise ValueError('Active and failure modes must be distinct valid labels.')
    if len(block) != len(modes)*dirty_dimension or not np.allclose(
            block.conj().T @ block, np.eye(len(block)), atol=1e-10, rtol=0):
        raise ValueError('The supplied local block must be a matching full unitary.')
    indices = [mode*dirty_dimension+j for mode in modes
               for j in range(dirty_dimension)]
    out = np.eye(mode_count*dirty_dimension, dtype=complex)
    out[np.ix_(indices, indices)] = block
    return out


def nested_attenuation(logical_dimension, supports, alphas):
    """Input-column attenuation for an increasing coordinate-support chain.

    This computes only the diagonal product of supported scalar factors.
    Factoring a stage product as W D additionally requires each intervening
    unitary to preserve its declared support and fix its complement.
    """
    if (isinstance(logical_dimension, bool)
            or not isinstance(logical_dimension, Integral) or logical_dimension < 1):
        raise ValueError('Logical dimension must be a positive integer.')
    supports, alphas = tuple(supports), tuple(alphas)
    if len(supports) != len(alphas):
        raise ValueError('Supports and attenuation factors must match.')
    diagonal = np.ones(logical_dimension)
    previous = set()
    for labels, alpha in zip(supports, alphas):
        labels = tuple(labels)
        current = set(labels)
        if (not labels or len(current) != len(labels) or not previous <= current
                or any(isinstance(i, bool) or not isinstance(i, Integral)
                       or not 0 <= i < logical_dimension for i in labels)
                or isinstance(alpha, bool) or not isinstance(alpha, Real)
                or not np.isfinite(alpha) or not 0 < alpha <= 1):
            raise ValueError('Use increasing valid supports and real 0 < alpha <= 1.')
        diagonal[list(labels)] *= alpha
        previous = current
    return diagonal


def diagonal_preweight(attenuation):
    """Contractive input preweight H=(2D)^(-1); no synthesis is implied."""
    values = np.asarray(attenuation)
    if (values.ndim != 1 or not len(values) or values.dtype.kind not in 'iuf'
            or not np.isfinite(values).all() or np.any(values < .5)
            or np.any(values > 1)):
        raise ValueError('Attenuation must be a nonempty real vector in [1/2, 1].')
    return .5/values


def three_call_error_bound(accepted_error):
    """Full-isometry bound for actual -U R U-dagger R U near accepted W/2.

    For epsilon <= 1/4 this is at most (4+sqrt(14))*epsilon <= 8*epsilon;
    the last inequality is strict when epsilon is positive.
    The argument includes rejection leakage, all dirty inputs and references.
    It does not permit a separately selected inverse or omission of -I.
    """
    if (isinstance(accepted_error, bool) or not isinstance(accepted_error, Real)
            or not np.isfinite(accepted_error)
            or not 0 <= accepted_error <= .25):
        raise ValueError('Accepted-block error must be real and in [0, 1/4].')
    return (4+np.sqrt(12+8*accepted_error))*accepted_error
