"""Histogram reconstruction for the coarse-frame raw-gradient protocol.

The Walsh butterflies use exact Python integers.  Coarse-block application
and the Hopf reverse traversal use NumPy floating point; this module does
not implement the theorem's certified arbitrary-precision arithmetic.
"""
from __future__ import annotations

import operator

import numpy as np

from .conventions import anchor_label, marker_label


def _integer_histogram(values: object, size: int, name: str) -> list[int]:
    try:
        array = np.asarray(values, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a one-dimensional integer histogram.") from exc
    if array.shape != (size,):
        raise ValueError(f"{name} must have exactly {size} entries.")
    result = []
    for value in array:
        try:
            count = operator.index(value)
        except TypeError:
            # Integral floating inputs denote their represented integers;
            # genuine integer inputs never pass through a float conversion.
            if (not isinstance(value, (float, np.floating))
                    or not np.isfinite(value) or value != np.floor(value)):
                raise ValueError(f"{name} must contain finite integers.") from None
            count = int(value)
        result.append(int(count))
    return result


def _integer_walsh(values: list[int]) -> list[int]:
    out = values.copy()
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            for offset in range(width):
                first, second = out[start + offset], out[start + width + offset]
                out[start + offset] = first + second
                out[start + width + offset] = first - second
        width *= 2
    return out


def decode_coarse_frame_histograms(
    theta: object,
    coarse_blocks: object,
    hist_x: object,
    hist_y: object,
    shots: int,
    coefficient_scale: float = 1,
) -> np.ndarray:
    """Decode signed X/Y histograms without forming a frame or Jacobian.

    ``theta`` contains the N-1 real Hopf angles in breadth-first heap order.
    ``coarse_blocks`` has shape (N-1, 2, 2).  Block ``node-1`` is the actual
    complex logical block on ``(anchor_label(node,n), marker_label(node,n))``;
    blocks act chronologically in node order 1,...,N-1.  Their literal phases
    must match the coarse circuit used in preparation and measurement.

    Each shot contributes its branch sign, multiplied by the selected
    observable term's sign, to one entry of ``hist_x`` or ``hist_y``.  Basis
    choices must have been independent fair X/Y draws.  ``shots`` is supplied
    separately because signed histograms lose counts through cancellation.
    ``coefficient_scale`` is the positive reflection-sum norm Lambda (one
    for a single observable).

    The return value is the empirical N-1 coordinate gradient, in heap order.
    Integer Walsh transforms copy their inputs and avoid fixed-width overflow;
    subsequent reconstruction is floating point.  There are O(N log N)
    arithmetic operations and O(N) stored values after the histograms and
    logical block coefficients exist; integer bit costs depend on shots.
    Inputs are not mutated.  This routine neither emits a native circuit nor
    certifies block unitarity, closeness to the target frame, or sampling error.
    """
    try:
        raw_angles = np.asarray(theta)
        if np.iscomplexobj(raw_angles) or any(np.iscomplexobj(value) for value in raw_angles.flat):
            raise ValueError("theta must contain real angles.")
        angles = np.asarray(theta, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("theta must be a finite one-dimensional real angle array.") from exc
    if angles.ndim != 1 or not np.all(np.isfinite(angles)):
        raise ValueError("theta must be a finite one-dimensional real angle array.")
    size = angles.size + 1
    if size < 2 or size & (size - 1):
        raise ValueError("theta must contain 2**n - 1 angles for n >= 1.")
    height = size.bit_length() - 1

    try:
        blocks = np.asarray(coarse_blocks, dtype=complex)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("coarse_blocks must contain finite complex 2-by-2 blocks.") from exc
    if blocks.shape != (size - 1, 2, 2) or not np.all(np.isfinite(blocks)):
        raise ValueError("coarse_blocks must have finite entries and shape (N-1, 2, 2).")
    try:
        count = operator.index(shots)
    except TypeError as exc:
        raise ValueError("shots must be a positive integer.") from exc
    if isinstance(shots, (bool, np.bool_)) or count <= 0:
        raise ValueError("shots must be a positive integer.")
    try:
        raw_scale = np.asarray(coefficient_scale)
        if raw_scale.ndim != 0 or np.iscomplexobj(raw_scale) or np.iscomplexobj(raw_scale.item()):
            raise ValueError("coefficient_scale must be a real scalar.")
        scale = float(coefficient_scale)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("coefficient_scale must be positive and finite.") from exc
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("coefficient_scale must be positive and finite.")

    hx = _integer_histogram(hist_x, size, "hist_x")
    hy = _integer_histogram(hist_y, size, "hist_y")
    if sum(abs(value) for value in hx) + sum(abs(value) for value in hy) > count:
        raise ValueError("The total absolute histogram count cannot exceed shots.")
    wx, wy = _integer_walsh(hx), _integer_walsh(hy)
    # Divide arbitrary-size integers before converting to complex floating
    # point.  The validated count bound makes both quotients lie in [-1,1].
    vector = np.array([complex(x / count, y / count)
                       for x, y in zip(wx, wy, strict=True)])
    for node, block in enumerate(blocks, start=1):
        indices = [anchor_label(node, height), marker_label(node, height)]
        vector[indices] = block @ vector[indices]

    cosine, sine = np.cos(angles), np.sin(angles)
    incoming = np.zeros(2 * size)
    incoming[1] = 1
    for node in range(1, size):
        incoming[2 * node] = incoming[node] * cosine[node - 1]
        incoming[2 * node + 1] = incoming[node] * sine[node - 1]
    adjoints = np.zeros(2 * size)
    adjoints[size:] = vector.real
    gradient = np.zeros(size - 1)
    for node in range(size - 1, 0, -1):
        c, s = cosine[node - 1], sine[node - 1]
        left, right = adjoints[2 * node], adjoints[2 * node + 1]
        gradient[node - 1] = incoming[node] * (-s * left + c * right)
        adjoints[node] = c * left + s * right
    return 4 * scale * gradient
