"""Exact combinatorial ledgers and finite matrix diagnostics for ordered banks.

No routine here synthesizes the fine left source or the complete endpoint.
Floating rank/positivity tests are diagnostics, not exact-real compiler oracles.
Native readout and prefix words live in ``native_ordered_banks``.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Sequence
import math

import numpy as np


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}.")
    return int(value)


@dataclass(frozen=True)
class Band:
    """One consecutive depth band [start, stop), in application order."""

    start: int
    stop: int

    def __post_init__(self) -> None:
        _integer(self.start, "start")
        _integer(self.stop, "stop", self.start + 1)

    @property
    def height(self) -> int:
        return self.stop - self.start

    @property
    def blocks(self) -> int:
        return 1 << self.start

    @property
    def block_size(self) -> int:
        return 1 << self.height

    @property
    def active_modes(self) -> int:
        return 1 << self.stop

    @property
    def squared_size(self) -> int:
        return 1 << (2 * self.stop - self.start)

    def supports(self, n: int) -> tuple[tuple[int, ...], ...]:
        n = _integer(n, "n", self.stop)
        step = 1 << (n - self.stop)
        return tuple(tuple((p * self.block_size + u) * step
                           for u in range(self.block_size))
                     for p in range(self.blocks))


def validate_bands(n: int, bands: Sequence[Band]) -> tuple[Band, ...]:
    n = _integer(n, "n", 1)
    bands = tuple(bands)
    cursor = 0
    for band in bands:
        if not isinstance(band, Band) or band.start != cursor or band.stop > n:
            raise ValueError("Bands must partition [0,n) consecutively.")
        cursor = band.stop
    if cursor != n:
        raise ValueError("Bands must cover every depth exactly once.")
    return bands


def band_schedule(n: int, singleton_tail: int = 0) -> tuple[Band, ...]:
    """Bottom-up heights 3,2,3,5,..., optionally after a singleton tail."""
    n = _integer(n, "n", 1)
    p = _integer(singleton_tail, "singleton_tail")
    if p > n:
        raise ValueError("singleton_tail cannot exceed n.")
    bottom = [1] * p
    left, j = n - p, 1
    while left:
        height = 3 if j == 1 else (1 << (j - 2)) + 1
        height = min(height, left)
        bottom.append(height)
        left -= height
        j += 1
    cursor, answer = 0, []
    for height in reversed(bottom):
        answer.append(Band(cursor, cursor + height))
        cursor += height
    return validate_bands(n, answer)


def refined_schedule(n: int) -> tuple[Band, ...]:
    """The (4+o(1))N volume schedule, using only integer logarithms."""
    n = _integer(n, "n", 3)
    ceil_log_n = (n - 1).bit_length()
    p = (ceil_log_n - 1).bit_length()
    return band_schedule(n, p)


def fixed_factor_schedule(n: int, factors: int) -> tuple[Band, ...]:
    """Rounded geometric endpoints; coincident endpoints are omitted."""
    n = _integer(n, "n", 1)
    k = _integer(factors, "factors", 1)
    denominator = (1 << k) - 1
    ends = [0]
    for j in range(1, k + 1):
        numerator = n * ((1 << k) - (1 << (k - j)))
        stop = (numerator + denominator - 1) // denominator
        if stop != ends[-1]:
            ends.append(stop)
    return validate_bands(n, [Band(a, b) for a, b in zip(ends, ends[1:])])


def ledger(n: int, bands: Sequence[Band]) -> dict[str, int]:
    bands = validate_bands(n, bands)
    N = 1 << n
    active = sum(b.active_modes for b in bands)
    return {"n": n, "logical_dimension": N, "factors": len(bands),
            "squared_size": sum(b.squared_size for b in bands),
            "history_modes": active, "shared_modes": N + active,
            "completed_modes": 2 * N + active, "available_modes": 4 * N}


def coarse_exponent(n: int) -> int:
    n = _integer(n, "n", 1)
    N = 1 << n
    root = math.isqrt(N)
    return root + (root * root < N) + 6


def modular_witness_exponent(n: int, factors: int) -> int:
    return _integer(factors, "factors", 1) * (coarse_exponent(n) + n + 2)


def ordered_product(factors: Sequence[np.ndarray]) -> np.ndarray:
    """First factor in the sequence acts first."""
    values = tuple(np.asarray(a, dtype=complex) for a in factors)
    if not values or values[0].ndim != 2 or values[0].shape[0] != values[0].shape[1]:
        raise ValueError("Supply nonempty square factors.")
    out = np.eye(len(values[0]), dtype=complex)
    for a in values:
        if a.shape != out.shape or not np.all(np.isfinite(a)):
            raise ValueError("Finite factors must have a common square shape.")
        out = a @ out
    return out


def grouped_factors(layers: Sequence[np.ndarray], bands: Sequence[Band]) -> tuple[np.ndarray, ...]:
    bands = validate_bands(len(layers), bands)
    return tuple(ordered_product(layers[b.start:b.stop]) for b in bands)


def interpolation_coefficients(target: Sequence[np.ndarray], coarse: Sequence[np.ndarray]) -> tuple[np.ndarray, ...]:
    """Coefficients of the formal ordered product, not an implemented circuit."""
    if not target or len(target) != len(coarse):
        raise ValueError("Target and coarse sequences must have equal positive length.")
    W, Q = [np.asarray(a, complex) for a in target], [np.asarray(a, complex) for a in coarse]
    ordered_product(W)
    ordered_product(Q)
    if W[0].shape != Q[0].shape:
        raise ValueError("Target and coarse factors must have the same dimension.")
    result = [np.eye(len(W[0]), dtype=complex)]
    for w, q in zip(W, Q):
        following = [np.zeros_like(q) for _ in range(len(result) + 1)]
        for k, value in enumerate(result):
            following[k] += q @ value
            following[k + 1] += (w - q) @ value
        result = following
    return tuple(result)


@dataclass(frozen=True)
class HistoryGeometry:
    n: int
    bands: tuple[Band, ...]
    heights: tuple[int, ...]
    rows: tuple[tuple[int, int], ...]
    matches: tuple[tuple[tuple[int, int], ...], ...]

    @property
    def dimension(self) -> int:
        return 1 << self.n

    @property
    def shared_modes(self) -> int:
        return len(self.rows)

    def weight(self, i: int, k: int) -> float:
        h = self.heights[i]
        if not 0 <= k <= h:
            raise ValueError("History order is outside the label's support.")
        return 2.0 ** (-(h if k == h else k + 1) / 2)

    def source(self) -> np.ndarray:
        out = np.zeros((4 * self.dimension, self.dimension), complex)
        for row, (i, k) in enumerate(self.rows):
            out[row, i] = self.weight(i, k)
        return out

    def unitary(self) -> np.ndarray:
        """Ideal full readout for bounded matrix diagnostics."""
        out = np.eye(4 * self.dimension, dtype=complex)
        for pairs in self.matches:
            for a, b in pairs:
                u, v = out[a].copy(), out[b].copy()
                out[a], out[b] = (u - v) / math.sqrt(2), (u + v) / math.sqrt(2)
        return out


def history_geometry(n: int, bands: Sequence[Band] | None = None) -> HistoryGeometry:
    n = _integer(n, "n", 1)
    bands = validate_bands(n, bands if bands is not None else [Band(d, d + 1) for d in range(n)])
    N = 1 << n
    heights = tuple(sum(i % (1 << (n - b.stop)) == 0 for b in bands) for i in range(N))
    rows = [(i, 0) for i in range(N)]
    for k in range(1, len(bands) + 1):
        rows.extend((i, k) for i in range(N) if heights[i] >= k)
    where = {pair: j for j, pair in enumerate(rows)}
    matches = tuple(tuple((where[i, k], where[i, k + 1]) for i in range(N) if heights[i] > k)
                    for k in range(len(bands)))
    if len(rows) + N > 4 * N:
        raise ValueError("History plus failure bank exceeds the available modes.")
    return HistoryGeometry(n, bands, heights, tuple(rows), matches)


def history_sources(coefficients: Sequence[np.ndarray], geometry: HistoryGeometry,
                    *, atol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Ideal completed source columns; no fine native left-source synthesis."""
    if not math.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and nonnegative.")
    C = tuple(np.asarray(c, complex) for c in coefficients)
    N, m = geometry.dimension, geometry.shared_modes
    if len(C) != len(geometry.bands) + 1 or any(c.shape != (N, N) for c in C):
        raise ValueError("Coefficient count or dimensions do not match the geometry.")
    if not all(np.all(np.isfinite(c)) for c in C):
        raise ValueError("Coefficients must be finite.")
    Q = C[0]
    if np.linalg.norm(Q.conj().T @ Q - np.eye(N), 2) > atol:
        raise ValueError("The coarse coefficient must be unitary.")
    for k in range(1, len(C)):
        outside = [i for i, h in enumerate(geometry.heights) if h < k]
        if outside and np.linalg.norm(C[k][:, outside], 2) > atol:
            raise ValueError("Coefficient violates its prescribed right support.")
    L = np.zeros((m, N), complex)
    L[:N] = Q.conj().T / math.sqrt(2)
    for row, (i, k) in enumerate(geometry.rows[N:], start=N):
        L[row] = C[k][:, i].conj() / (2 * geometry.weight(i, k))
    deficit = np.eye(N) - L.conj().T @ L
    values, vectors = np.linalg.eigh(deficit)
    if values.min() < -atol:
        raise ValueError("Left rows are not a contraction.")
    D = (vectors * np.sqrt(np.maximum(values, 0))) @ vectors.conj().T
    left = np.zeros((4 * N, N), complex)
    left[:m] = L
    left[m:m + N] = Q.conj().T @ D
    return left, geometry.source()


def prefix_matchings(n: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    n = _integer(n, "n", 1)
    offsets = [(1 << (d + 1)) - 2 for d in range(n)]
    return tuple(tuple((offsets[d] + 2 * p * (1 << (d - ell)),
                        offsets[d] + (2 * p + 1) * (1 << (d - ell)), p)
                       for d in range(ell + 1, n) for p in range(1 << ell))
                 for ell in range(n - 1))


def prefix_matrix(n: int, blocks: Sequence[Sequence[np.ndarray]]) -> np.ndarray:
    """Ideal direct sum of the actual coarse prefixes, padded to 2N modes."""
    n = _integer(n, "n", 1)
    if len(blocks) != n - 1 or any(len(row) != 1 << d for d, row in enumerate(blocks)):
        raise ValueError("Supply one row of 2**d blocks for each d < n-1.")
    out = np.eye(1 << (n + 1), dtype=complex)
    for ell, pairs in enumerate(prefix_matchings(n)):
        for a, b, p in pairs:
            block = np.asarray(blocks[ell][p], complex)
            if block.shape != (2, 2) or not np.all(np.isfinite(block)):
                raise ValueError("Each prefix block must be a finite 2x2 matrix.")
            out[[a, b]] = block @ out[[a, b]]
    return out
