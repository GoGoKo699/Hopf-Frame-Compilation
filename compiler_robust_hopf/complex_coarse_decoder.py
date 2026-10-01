"""Floating-point reconstruction for the gauged complex Hopf protocol.

Walsh butterflies retain exact Python integers. Recorded native blocks,
phase evaluation, and the final contractions use NumPy floating point.
This utility neither emits a native compiler nor provides the theorem's
certified arithmetic. Large supplied phases may require additional input
precision and certified argument reduction beyond these floating routines.
"""
from __future__ import annotations

from numbers import Real

import numpy as np

from .coarse_frame_decoder import (
    _apply_tree_blocks,
    _coarse_blocks,
    _histogram_vector,
    _real_tree_angles,
    _real_tree_pullback,
)


def _leaf_phases(values: object) -> np.ndarray:
    """Require finite real numeric phases, including in object arrays."""
    try:
        raw = np.asarray(values)
        if (raw.ndim != 1 or np.iscomplexobj(raw)
                or any(not isinstance(value, (Real, np.bool_)) for value in raw)):
            raise ValueError("leaf_phases must contain real numeric phases.")
        phases = np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("leaf_phases must be a finite one-dimensional real array.") from exc
    size = phases.size
    if size < 2 or size & (size - 1) or not np.all(np.isfinite(phases)):
        raise ValueError("leaf_phases must contain 2**n finite real phases for n >= 1.")
    return phases


def _phase_tree(phases: np.ndarray) -> tuple[float, np.ndarray]:
    size = phases.size
    means = np.zeros(2 * size)
    means[size:] = phases
    angles = np.empty(size - 1)
    for node in range(size - 1, 0, -1):
        left, right = means[2 * node], means[2 * node + 1]
        # Half first: neither a same-sign sum nor an opposite-sign
        # difference can overflow for finite floating input phases.
        means[node] = left / 2 + right / 2
        angles[node - 1] = right / 2 - left / 2
    return float(means[1]), angles


def phase_tree_angles(leaf_phases: object) -> tuple[float, np.ndarray]:
    """Return the mean phase and determinant-one prefix angles in heap order.

    Supplied phases are unwrapped real values. At each prefix p the angle
    is (mean_right-mean_left)/2. The row diag(exp(-i angle), exp(i angle))
    acts on the next bit for every suffix. In exact arithmetic their product
    is diag(exp(i*(leaf_phase-mean))). Pairwise half-sums avoid overflow in
    the mean, but do not certify accuracy for very large or ill-conditioned
    phase inputs. The supplied array is not mutated.
    """
    return _phase_tree(_leaf_phases(leaf_phases))


def apply_phase_blocks(vector: object, phase_blocks: object) -> np.ndarray:
    """Apply actual prefix-addressed rows across every suffix to a copy.

    ``vector`` has N=2**n entries, n>=1. ``phase_blocks`` has shape
    (N-1,2,2), in breadth-first heap order. Row p acts on the next bit
    throughout its prefix subtree, including every lower-bit suffix.
    Rows act chronologically in heap order; actual native approximants
    need not be diagonal or mutually commuting. This costs O(N log N)
    scalar arithmetic and O(N) storage. Unitarity is not certified.
    """
    try:
        out = np.array(vector, dtype=complex, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("vector must be a finite one-dimensional complex array.") from exc
    size = out.size
    if (out.ndim != 1 or size < 2 or size & (size - 1)
            or not np.all(np.isfinite(out))):
        raise ValueError("vector must contain 2**n finite entries for n >= 1.")
    blocks = _coarse_blocks(phase_blocks, size, "phase_blocks")
    for node, block in enumerate(blocks, start=1):
        depth = node.bit_length() - 1
        width = size >> depth
        start = (node - (1 << depth)) * width
        pairs = out[start:start + width].reshape(2, width // 2)
        out[start:start + width] = (block @ pairs).reshape(width)
    return out


def decode_complex_coarse_frame_histograms(
    theta: object,
    coarse_blocks: object,
    phase_blocks: object,
    leaf_phases: object,
    hist_x: object,
    hist_y: object,
    shots: int,
    *,
    coefficient_scale: float = 1,
) -> np.ndarray:
    """Decode the N-1 complex-state magnitude coordinates from X/Y records.

    ``theta`` and ``coarse_blocks`` describe the real Hopf tree and its
    actual native coarse tree, as in ``decode_coarse_frame_histograms``.
    ``phase_blocks`` are the actual native prefix rows from
    ``apply_phase_blocks``; together they implement C=C_P C_R.
    ``leaf_phases`` are the original unwrapped real phases. Their mean
    fixes G=diag(exp(i*(leaf_phase-mean))), so the exact target formula is
    4*Lambda*J_R.T*Re(G-dagger*C*u), u=F_N*(hist_x+i*hist_y)/shots.

    Histograms contain branch signs times observable-term signs, with
    fresh independent fair X/Y basis choices. ``coefficient_scale`` is
    the positive observable coefficient norm Lambda. Separate phase
    coordinates use their direct one-hot stream and are not returned.
    All inputs are preserved. There are O(N log N) arithmetic operations
    and O(N) stored values; integer costs depend on shots, and classical
    phase/block preprocessing and precision remain separately charged.
    This floating utility does not certify sampling or approximation error.
    """
    angles, size, height = _real_tree_angles(theta)
    blocks = _coarse_blocks(coarse_blocks, size)
    phases = _leaf_phases(leaf_phases)
    if phases.size != size:
        raise ValueError("leaf_phases must have one entry per logical leaf.")
    mean, _ = _phase_tree(phases)
    vector, scale = _histogram_vector(hist_x, hist_y, shots, coefficient_scale, size)
    vector = _apply_tree_blocks(vector, blocks, height)
    vector = apply_phase_blocks(vector, phase_blocks)
    # Separate exponentials avoid overflowing the subtraction phase-mean.
    # Their floating argument reduction is not a certified phase evaluator.
    gauge_adjoint = np.exp(-1j * phases) * np.exp(1j * mean)
    return 4 * scale * _real_tree_pullback(angles, (gauge_adjoint * vector).real)
