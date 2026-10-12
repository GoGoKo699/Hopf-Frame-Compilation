"""Finite matrix diagnostics for the original ordered tree and its residual.

These routines do not emit native gates, supply source isometries, or certify
arbitrary-precision input enclosures. In particular, the shared-mode factors
below are not a completed block encoding. Numerical rank and singularity
tolerances are diagnostic choices, not equality oracles for a compiler.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from .conventions import anchor_label, marker_label


def _tree(dimension, ordered_edges):
    if isinstance(dimension, bool) or not isinstance(dimension, Integral) or dimension < 2:
        raise ValueError("dimension must be an integer at least two.")
    edges = tuple(tuple(edge) for edge in ordered_edges)
    if len(edges) != dimension - 1:
        raise ValueError("A tree must have dimension - 1 edges.")
    parents = list(range(dimension))
    adjacency = [[] for _ in range(dimension)]

    def root(vertex):
        while parents[vertex] != vertex:
            vertex = parents[vertex]
        return vertex

    for time, edge in enumerate(edges):
        if len(edge) != 2 or any(isinstance(v, bool) or not isinstance(v, Integral)
                                 or not 0 <= v < dimension for v in edge):
            raise ValueError("Each edge must contain two valid integer labels.")
        first, second = edge
        a, b = root(first), root(second)
        if a == b:
            raise ValueError("Edges must form a tree without loops or repeated edges.")
        parents[a] = b
        # K[first, second] = -1 in the repository's Ry convention.
        adjacency[first].append((second, time, -1))
        adjacency[second].append((first, time, 1))
    return edges, adjacency


def _unitary(matrix):
    result = np.asarray(matrix, dtype=complex)
    if result.ndim != 2 or result.shape[0] != result.shape[1] or not len(result):
        raise ValueError("Expected a nonempty square matrix.")
    if not np.all(np.isfinite(result)):
        raise ValueError("Matrix entries must be finite.")
    if not np.allclose(result.conj().T @ result, np.eye(len(result)), atol=1e-10, rtol=0):
        raise ValueError("Expected a unitary matrix to numerical tolerance.")
    return result


def hopf_tree_edges(n: int) -> tuple[tuple[int, int], ...]:
    """Original computational-label edges in breadth-first application order."""
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer.")
    return tuple((anchor_label(node, n), marker_label(node, n))
                 for node in range(1, 1 << n))


def tree_word_matrix(dimension: int, ordered_edges, blocks) -> np.ndarray:
    """Multiply two-level SU(2) blocks, first applied block listed first.

    A block is embedded in its edge's ordered pair. Complex phases are kept
    literally; no determinant correction or polar projection is performed.
    """
    edges, _ = _tree(dimension, ordered_edges)
    blocks = tuple(blocks)
    if len(blocks) != len(edges):
        raise ValueError("There must be one block per ordered edge.")
    result = np.eye(dimension, dtype=complex)
    for edge, block in zip(edges, blocks, strict=True):
        local = _unitary(block)
        if local.shape != (2, 2) or not np.isclose(np.linalg.det(local), 1, atol=1e-10, rtol=0):
            raise ValueError("Each tree block must be in SU(2).")
        indices = list(edge)
        result[indices, :] = local @ result[indices, :]
    return result


def tree_edge_cuts(dimension: int, ordered_edges):
    """Return both label sets of each edge cut, in the supplied edge order."""
    edges, adjacency = _tree(dimension, ordered_edges)
    cuts = []
    for excluded, (start, _) in enumerate(edges):
        seen, stack = {start}, [start]
        while stack:
            for neighbour, time, _ in adjacency[stack.pop()]:
                if time != excluded and neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        cuts.append((tuple(sorted(seen)), tuple(v for v in range(dimension) if v not in seen)))
    return tuple(cuts)


def tree_cayley_from_paths(dimension: int, ordered_edges, half_tangents) -> np.ndarray:
    """Build the signed path formula without multiplying or inverting W.

    The inputs are finite ``tan(theta_e / 2)`` values; singular half-angle
    charts are deliberately outside this interface. Fraction inputs retain
    exact object arithmetic, useful for finite algebraic regression checks.
    """
    edges, adjacency = _tree(dimension, ordered_edges)
    tangents = np.asarray(half_tangents)
    if (tangents.shape != (len(edges),) or np.iscomplexobj(tangents)
            or any(not isinstance(value, Real) for value in tangents)):
        raise ValueError("Supply one real half-tangent per edge.")
    if not np.all(np.isfinite(np.asarray(tangents, dtype=float))):
        raise ValueError("Half-tangents must be finite; singular charts are excluded.")
    result = np.zeros((dimension, dimension), dtype=object if tangents.dtype == object else float)
    for source in range(dimension):
        stack = [(source, -1, -1, 1)]
        while stack:
            vertex, parent, previous_time, value = stack.pop()
            for neighbour, time, orientation in adjacency[vertex]:
                if neighbour == parent:
                    continue
                comparison = 1 if previous_time < 0 or previous_time > time else -1
                next_value = value * orientation * tangents[time] * comparison
                result[source, neighbour] = next_value
                stack.append((neighbour, vertex, time, next_value))
    return result


def cayley_matrix(unitary, *, singular_atol: float = 1e-12) -> np.ndarray:
    """Numerically evaluate (W-I)(W+I)^-1, rejecting an unresolved chart."""
    if not np.isfinite(singular_atol) or singular_atol < 0:
        raise ValueError("singular_atol must be finite and nonnegative.")
    matrix = _unitary(unitary)
    identity = np.eye(len(matrix))
    plus = matrix + identity
    if np.linalg.svd(plus, compute_uv=False)[-1] <= singular_atol:
        raise ValueError("The Cayley chart is singular or unresolved at the supplied tolerance.")
    return np.linalg.solve(plus.T, (matrix - identity).T).T


@dataclass(frozen=True)
class ResidualDiagnostics:
    residual: np.ndarray
    cayley: np.ndarray
    delta: float
    cayley_norm_bound: float
    absolute_sum_bound: float


def stable_residual(target, coarse) -> ResidualDiagnostics:
    """Check the complete logical residual E=W Q†, including complex phases.

    This accepts logical matrices only. It asserts nothing about physical
    approximation errors containing operators on arbitrary dirty work.
    """
    target, coarse = _unitary(target), _unitary(coarse)
    if target.shape != coarse.shape:
        raise ValueError("Target and coarse matrices must have the same shape.")
    residual = target @ coarse.conj().T
    delta = float(np.linalg.norm(residual - np.eye(len(residual)), 2))
    if delta >= 2:
        raise ValueError("The stable residual estimate requires delta < 2.")
    return ResidualDiagnostics(residual, cayley_matrix(residual), delta,
                               delta / (2 - delta), 1 + np.sqrt(len(residual)) * delta)


@dataclass(frozen=True)
class SharedModeFactors:
    """Rectangular factors, explicitly not native source isometries."""

    a: np.ndarray
    b: np.ndarray
    gram_bound: float
    reconstruction_error: float

    @property
    def shared_modes(self) -> int:
        return self.a.shape[1]

    @property
    def separate_failure_completion_modes(self) -> int:
        """Count both additional N-mode Gram-deficit square-root banks."""
        return self.shared_modes + 2 * self.a.shape[0]


def hierarchical_shared_factors(
    residual, *, max_cut_rank: int = 2, rank_atol: float = 1e-10
) -> SharedModeFactors:
    """Numerically factor a dyadic rank-one or rank-two-cut residual.

    Rank one uses 3N-2 shared modes; rank two uses 4N-4. Each directed sibling
    block retains the requested number of SVD components (one at the bottom
    level), including zero components. Larger discarded singular values are
    rejected; the returned reconstruction error exposes numerical truncation.
    Rank one applies to W Q† when both words follow the same outward rooted
    tree order. Arbitrarily reordered tree words can require rank two.
    Orthogonal completion, inverse source words, and their native cost remain
    unimplemented. The 5N-2 or 6N-4 separate-failure counts are not dimension
    lower bounds on all possible joint completions.
    """
    if (isinstance(max_cut_rank, bool) or not isinstance(max_cut_rank, Integral)
            or max_cut_rank not in (1, 2)):
        raise ValueError("max_cut_rank must be the integer 1 or 2.")
    if not np.isfinite(rank_atol) or rank_atol < 0:
        raise ValueError("rank_atol must be finite and nonnegative.")
    matrix = _unitary(residual)
    size = len(matrix)
    if size < 2 or size & (size - 1):
        raise ValueError("The hierarchy requires a power-of-two dimension at least two.")
    n = size.bit_length() - 1
    delta = float(np.linalg.norm(matrix - np.eye(size), 2))
    shared_modes = 3 * size - 2 if max_cut_rank == 1 else 4 * size - 4
    a, b = np.zeros((size, shared_modes), complex), np.zeros((size, shared_modes), complex)
    diagonal = np.diag(matrix)
    roots = np.sqrt(np.abs(diagonal))
    a[np.arange(size), np.arange(size)] = np.divide(diagonal, roots,
                                                  out=np.zeros(size, complex), where=roots > 0)
    b[np.arange(size), np.arange(size)] = roots
    column = size
    for depth in range(n):
        width = size >> depth
        half = width // 2
        for start in range(0, size, width):
            left, right = slice(start, start + half), slice(start + half, start + width)
            for rows, columns in ((left, right), (right, left)):
                u, singular, vh = np.linalg.svd(matrix[rows, columns], full_matrices=False)
                keep = min(max_cut_rank, half)
                if len(singular) > keep and singular[keep] > rank_atol:
                    raise ValueError(f"A directed sibling block has numerical rank greater than {max_cut_rank}.")
                selected = slice(column, column + keep)
                weights = np.sqrt(singular[:keep])
                a[rows, selected] = u[:, :keep] * weights
                b[columns, selected] = vh[:keep, :].conj().T * weights
                column += keep
    assert column == shared_modes
    error = float(np.linalg.norm(a @ b.conj().T - matrix, 2))
    return SharedModeFactors(a, b, 1 + n * delta, error)
