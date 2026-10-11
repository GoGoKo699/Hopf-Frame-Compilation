"""Independent complete-matrix and graph checks for structured forest factors.

The tests verify identities and finite support budgets, not a native emitter
or the analytical asymptotic gate bound.
"""
import unittest

import numpy as np

from compiler_robust_hopf.structured_forests import (
    computational_edges,
    partition_degree_bands,
    partition_suffix_forest,
    trailing_degree,
)


ATOL = 2e-12


def _reference_edges(height):
    return tuple((prefix << (height - depth),
                  (2 * prefix + 1) << (height - depth - 1))
                 for depth in range(height) for prefix in range(1 << depth))


def _word(size, edges, angles):
    result = np.eye(size)
    for left, right in edges:
        cosine, sine = np.cos(angles[left, right]), np.sin(angles[left, right])
        first, second = result[left].copy(), result[right].copy()
        result[left] = cosine * first - sine * second
        result[right] = sine * first + cosine * second
    return result


def _factors(partition, angles):
    size = 1 << partition.height
    component_words = []
    for index, _ in enumerate(partition.components):
        edges = [edge for edge in partition.internal_edges
                 if partition.component_of[edge[0]] == index]
        component_words.append(_word(size, edges, angles))
    internal = _word(size, partition.internal_edges, angles)
    sandwiches = []
    for depth, edges in enumerate(partition.cut_groups):
        parent_word = np.eye(size)
        for index, local in enumerate(component_words):
            if partition.component_depths[index] == depth:
                parent_word = local @ parent_word
        cut_word = _word(size, edges, angles)
        sandwiches.append(parent_word.T @ cut_word @ parent_word)
    grouped = np.eye(size)
    for sandwich in sandwiches:
        grouped = sandwich @ grouped
    return internal @ grouped, sandwiches


def _bounded_support(height, bound, rng, maximal):
    children = [[] for _ in range(1 << height)]
    for left, right in _reference_edges(height):
        children[left].append(right)
    active = set()
    for parent, choices in enumerate(children):
        if maximal:
            chosen = choices[:bound]
        else:
            count = int(rng.integers(min(bound, len(choices)) + 1))
            chosen = rng.choice(choices, size=count, replace=False)
        active.update((parent, int(child)) for child in chosen)
    return tuple(edge for edge in _reference_edges(height) if edge in active)


def _stable_star_transfer(height, angles, exceptional):
    """Independent local-matrix construction of the exact input correction."""
    size = 1 << height
    children = [[] for _ in range(size)]
    for left, right in _reference_edges(height):
        children[left].append(right)
    correction = np.eye(size)
    reordered, new_angles = [], {}

    def visit(vertex):
        leaves = children[vertex]
        if not leaves:
            return
        local_edges = tuple((0, index) for index in range(1, len(leaves) + 1))
        local_angles = {(0, index): angles[vertex, child]
                        for index, child in enumerate(leaves, 1)}
        original = _word(len(leaves) + 1, local_edges, local_angles)
        root = original[:, 0]
        order = [child for child in leaves if (vertex, child) not in exceptional]
        order += [child for child in leaves if (vertex, child) in exceptional]
        remaining, replacement = 1., np.eye(len(leaves) + 1)
        for child in order:
            index = leaves.index(child) + 1
            sine = root[index] / np.sqrt(remaining)
            cosine = np.sqrt(max(0., 1. - sine * sine))
            angle = np.arctan2(sine, cosine)
            replacement = _word(len(leaves) + 1, ((0, index),),
                                {(0, index): angle}) @ replacement
            reordered.append((vertex, child))
            new_angles[vertex, child] = angle
            remaining -= root[index] ** 2
        local_correction = replacement.T @ original
        np.testing.assert_allclose(local_correction[:, 0],
                                   np.eye(len(leaves) + 1)[:, 0], atol=ATOL, rtol=0)
        correction[np.ix_(leaves, leaves)] = local_correction[1:, 1:]
        for child in leaves:
            visit(child)

    visit(0)
    return tuple(reordered), new_angles, correction


class StructuredForestTests(unittest.TestCase):
    def test_degree_band_complete_factors_include_signed_angles_and_corners(self):
        for height in range(2, 7):
            for bound, bands in ((2, 3), (3, 8)):
                for maximal in (False, True):
                    with self.subTest(height=height, bound=bound, maximal=maximal):
                        rng = np.random.default_rng(713 * height + bound)
                        edges = _bounded_support(height, bound, rng, maximal)
                        values = rng.uniform(-np.pi, np.pi, len(edges))
                        values[::4] = np.resize([0., np.pi / 2, -np.pi / 2,
                                               np.pi, -np.pi], len(values[::4]))
                        angles = dict(zip(edges, values))
                        partition = partition_degree_bands(height, edges, bands)
                        actual, _ = _factors(partition, angles)
                        expected = _word(1 << height, edges, angles)
                        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                        self.assertLessEqual(len(partition.cut_groups), bands - 1)
                        self.assertEqual(partition.square_budget,
                                         sum(len(c) ** 2 for c in partition.components))
                        self.assertEqual(sorted(v for c in partition.components for v in c),
                                         list(range(1 << height)))

    def test_stable_exceptional_transfer_keeps_the_entire_frame(self):
        for height in range(2, 7):
            rng = np.random.default_rng(987 + height)
            edges = _reference_edges(height)
            angles = dict(zip(edges, rng.uniform(-.74, .74, len(edges))))
            # Distinct intermediate magnitudes and arbitrary interleaving are
            # intentionally permitted; the algebra requires no weak cutoff.
            exceptional = {edge for edge in edges if rng.random() < .5}
            reordered, new_angles, correction = _stable_star_transfer(
                height, angles, exceptional)
            partition = partition_suffix_forest(height, reordered, exceptional)
            factored, _ = _factors(partition, new_angles)
            expected = _word(1 << height, edges, angles)
            np.testing.assert_allclose(factored @ correction, expected,
                                       atol=ATOL, rtol=0)
            retained_children = [[] for _ in range(1 << height)]
            for parent, child in reordered:
                degree = trailing_degree(parent, height)
                old_sine = abs(np.sin(angles[parent, child]))
                new_sine = abs(np.sin(new_angles[parent, child]))
                self.assertGreaterEqual(new_sine + ATOL, 2 ** (-degree / 2) * old_sine)
                self.assertLessEqual(new_sine, 2 ** (degree / 2) * old_sine + ATOL)
                self.assertGreaterEqual(np.cos(new_angles[parent, child]) + ATOL,
                                        2 ** (-degree / 2))
                if (parent, child) not in exceptional:
                    retained_children[parent].append(child)
            inherited_depth = [0] * (1 << height)
            for parent in range((1 << height) - 1, -1, -1):
                children = retained_children[parent]
                inherited_depth[parent] = max(
                    [len(children)]
                    + [index + inherited_depth[child]
                       for index, child in enumerate(children, 1)])
            self.assertLessEqual(max(inherited_depth), height)

    def test_depth_groups_must_be_reversed_in_matrix_multiplication(self):
        edges = ((0, 2), (2, 3))
        angles = {(0, 2): .37, (2, 3): -.61}
        partition = partition_suffix_forest(2, edges, edges)
        actual, factors = _factors(partition, angles)
        expected = _word(4, edges, angles)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(factors[0] @ factors[1] - expected, 2), .1)

    def test_root_degree_counts_and_fibonacci_growth_certificate(self):
        gamma, constant = 5 / 3, 25
        counts = [1, 2]
        for degree in range(2, 21):
            counts.append(1 + counts[-1] + counts[-2])
        for degree, count in enumerate(counts):
            self.assertLessEqual(count, constant * gamma ** degree)
        self.assertLess(gamma ** 4, 8)
        for height in range(1, 11):
            self.assertEqual(computational_edges(height), _reference_edges(height))
            for degree in range(height):
                self.assertEqual(sum(trailing_degree(v, height) == degree
                                     for v in range(1, 1 << height)),
                                 (1 << height) // (1 << (degree + 1)))
            rng = np.random.default_rng(height)
            active = _bounded_support(height, 2, rng, True)
            partition = partition_degree_bands(height, active, 3)
            rho = gamma ** (4 / 3) / 2
            kappa = constant ** 2 * gamma ** (1 / 3) * (1 + 1 / (2 * (1 - rho)))
            self.assertLessEqual(partition.square_budget, kappa * (1 << height))
            self.assertLessEqual(partition.factor_count, 7)

    def test_invalid_orders_and_support_are_rejected(self):
        for height, edges, cuts in (
            (2, ((2, 3), (0, 2)), ((0, 2),)),
            (2, ((0, 2), (0, 1)), ((0, 2),)),
            (2, ((0, 2),), ((2, 3),)),
            (2, ((0, 3),), ()),
            (2, ((0, 2), (0, 2)), ()),
        ):
            with self.subTest(edges=edges, cuts=cuts), self.assertRaises(ValueError):
                partition_suffix_forest(height, edges, cuts)
        for bands in (0, -1, 1.5, True):
            with self.subTest(bands=bands), self.assertRaises(ValueError):
                partition_degree_bands(2, (), bands)
        for height in (0, -1, 1.5, True):
            with self.subTest(height=height), self.assertRaises(ValueError):
                computational_edges(height)


if __name__ == "__main__":
    unittest.main()
