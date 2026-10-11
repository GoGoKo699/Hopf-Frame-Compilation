"""Combinatorial certificates for the structured-forest endpoint theorems.

These helpers partition an active computational-label forest and expose the
exact chronological factorization. They do not emit native Clifford+T gates
or estimate a numerical circuit's approximation error.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Iterable


Edge = tuple[int, int]


def _height(height: int) -> int:
    if not isinstance(height, Integral) or isinstance(height, bool) or height < 1:
        raise ValueError("height must be a positive integer")
    return int(height)


def trailing_degree(vertex: int, height: int) -> int:
    """Return the number of available child bit positions; zero has degree n."""
    height = _height(height)
    if (not isinstance(vertex, Integral) or isinstance(vertex, bool)
            or not 0 <= vertex < 1 << height):
        raise ValueError("vertex must be a computational label")
    vertex = int(vertex)
    return height if vertex == 0 else (vertex & -vertex).bit_length() - 1


def computational_edges(height: int) -> tuple[Edge, ...]:
    """Return all original Hopf edges in chronological layer order."""
    height = _height(height)
    size = 1 << height
    return tuple((left, left + (1 << bit))
                 for bit in range(height - 1, -1, -1)
                 for left in range(0, size, 1 << (bit + 1)))


def _edges(height: int, edges: Iterable[Edge]) -> tuple[Edge, ...]:
    available = set(computational_edges(height))
    try:
        result = tuple(tuple(edge) for edge in edges)
    except TypeError as exc:
        raise ValueError("edges must be computational parent-child pairs") from exc
    if any(len(edge) != 2 or any(not isinstance(v, Integral) or isinstance(v, bool)
                                 for v in edge) for edge in result):
        raise ValueError("edges must be computational parent-child pairs")
    if len(set(result)) != len(result) or any(edge not in available for edge in result):
        raise ValueError("edges must be distinct computational parent-child pairs")
    return tuple((int(left), int(right)) for left, right in result)


@dataclass(frozen=True)
class ForestPartition:
    """Finite support data for U_all F_(k-1)...F_0, F_r=V_r^* K_r V_r.

    Matrices use ordinary rightmost-first multiplication. ``edges`` retain
    chronological order, and every ``cut_groups[r]`` retains its subsequence.
    Component entries are disjoint computational labels, including isolates.
    """

    height: int
    edges: tuple[Edge, ...]
    cuts: frozenset[Edge]
    components: tuple[tuple[int, ...], ...]
    component_of: tuple[int, ...]
    component_depths: tuple[int, ...]
    cut_groups: tuple[tuple[Edge, ...], ...]

    @property
    def internal_edges(self) -> tuple[Edge, ...]:
        return tuple(edge for edge in self.edges if edge not in self.cuts)

    @property
    def square_budget(self) -> int:
        return sum(len(component) ** 2 for component in self.components)

    @property
    def cut_star_square_budgets(self) -> tuple[int, ...]:
        budgets = []
        for group in self.cut_groups:
            degrees: dict[int, int] = {}
            for parent, _ in group:
                degrees[parent] = degrees.get(parent, 0) + 1
            budgets.append(sum((degree + 1) ** 2 for degree in degrees.values()))
        return tuple(budgets)

    @property
    def factor_count(self) -> int:
        """Count U_all and three factors for every nonempty quotient depth."""
        return 1 + 3 * len(self.cut_groups)


def partition_suffix_forest(
    height: int, edge_word: Iterable[Edge], cuts: Iterable[Edge]
) -> ForestPartition:
    """Partition a chronological word whose cuts form every vertex's suffix.

    The word may use a different local child order after exact same-root star
    transfer, but an incoming active edge must precede every outgoing edge.
    Invalid orders are rejected rather than silently commuting rotations.
    """
    height = _height(height)
    word = _edges(height, edge_word)
    cut_set = frozenset(_edges(height, cuts))
    if not cut_set.issubset(word):
        raise ValueError("every cut must belong to the active edge word")
    size = 1 << height
    parents: list[int | None] = [None] * size
    positions = {edge: index for index, edge in enumerate(word)}
    for left, right in word:
        parents[right] = left
    suffix_started: set[int] = set()
    for edge in word:
        left, _ = edge
        ancestor = parents[left]
        if ancestor is not None and positions[ancestor, left] >= positions[edge]:
            raise ValueError("incoming edges must precede outgoing edges")
        if edge in cut_set:
            suffix_started.add(left)
        elif left in suffix_started:
            raise ValueError("cut edges must form each outgoing suffix")

    components: list[list[int]] = []
    component_of = [-1] * size
    depths: list[int] = []
    # Every computational parent has smaller integer label than its child.
    for vertex, parent in enumerate(parents):
        if parent is None or (parent, vertex) in cut_set:
            index = len(components)
            components.append([])
            depths.append(0 if parent is None else depths[component_of[parent]] + 1)
        else:
            index = component_of[parent]
        component_of[vertex] = index
        components[index].append(vertex)
    groups: list[list[Edge]] = [[] for _ in range(max(depths, default=0))]
    for edge in word:
        if edge in cut_set:
            groups[depths[component_of[edge[0]]]].append(edge)
    return ForestPartition(
        height, word, cut_set, tuple(map(tuple, components)), tuple(component_of),
        tuple(depths), tuple(map(tuple, groups)),
    )


def partition_degree_bands(
    height: int, active_edges: Iterable[Edge], bands: int
) -> ForestPartition:
    """Cut each active component into at most ``bands`` degree windows.

    Each original component chooses its window width from its own root degree.
    The original child order is already decreasing degree, so the cuts are
    suffixes. No angle values, equality tests, or numerical rounding are used.
    The linear asymptotic square budget additionally requires the growth
    certificate stated in STRUCTURED_FOREST_COMPILERS.md.
    """
    height = _height(height)
    if not isinstance(bands, Integral) or isinstance(bands, bool) or bands < 1:
        raise ValueError("bands must be a positive integer")
    bands = int(bands)
    active = set(_edges(height, active_edges))
    word = tuple(edge for edge in computational_edges(height) if edge in active)
    size = 1 << height
    parents: list[int | None] = [None] * size
    for left, right in word:
        parents[right] = left
    roots = list(range(size))
    ordinal = [0] * size
    for vertex, parent in enumerate(parents):
        if parent is not None:
            roots[vertex] = roots[parent]
        root_degree = trailing_degree(roots[vertex], height)
        width = (root_degree + 1 + bands - 1) // bands
        ordinal[vertex] = (root_degree - trailing_degree(vertex, height)) // width
    cuts = tuple((left, right) for left, right in word
                 if ordinal[left] != ordinal[right])
    return partition_suffix_forest(height, word, cuts)
