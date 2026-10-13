"""Independent full-matrix checks for the bounded exceptional-exit theorem.

The support envelopes are computed from edge words, including at exact corners.
These finite checks do not implement or claim a native gate emitter.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.structured_forests import partition_exceptional_exits


ATOL = 4e-12


def _degree(label, height):
    return height if label == 0 else (label & -label).bit_length() - 1


def _edges(height):
    return tuple((prefix << (height - depth),
                  (2 * prefix + 1) << (height - depth - 1))
                 for depth in range(height) for prefix in range(1 << depth))


def _word(size, edges, coefficients):
    result = np.eye(size)
    for edge in edges:
        left, right = edge
        cosine, sine = coefficients[edge]
        first, second = result[left].copy(), result[right].copy()
        result[left] = cosine * first - sine * second
        result[right] = sine * first + cosine * second
    return result


def _coefficients(edges, rng):
    return {edge: (np.cos(angle), np.sin(angle))
            for edge, angle in zip(edges, rng.uniform(-.74, .74, len(edges)))}


def _suffix_order(height, exceptional):
    children = [[] for _ in range(1 << height)]
    for left, right in _edges(height):
        children[left].append(right)
    result = []

    def visit(vertex):
        leaves = children[vertex]
        result.extend((vertex, child) for child in leaves
                      if (vertex, child) not in exceptional)
        result.extend((vertex, child) for child in leaves
                      if (vertex, child) in exceptional)
        for child in leaves:
            visit(child)

    visit(0)
    return tuple(result)


def _transfer(height, coefficients, exceptional):
    """Match each local star column and retain its complete child correction."""
    size = 1 << height
    correction, replacement_coefficients = np.eye(size), {}
    for vertex in range(size):
        leaves = [right for left, right in _edges(height) if left == vertex]
        if not leaves:
            continue
        original = _word(len(leaves) + 1,
                         tuple((0, i) for i in range(1, len(leaves) + 1)),
                         {(0, i): coefficients[vertex, child]
                          for i, child in enumerate(leaves, 1)})
        order = [child for child in leaves if (vertex, child) not in exceptional]
        order += [child for child in leaves if (vertex, child) in exceptional]
        root_column = original[:, 0]
        remaining = set(range(len(leaves) + 1))
        local_coefficients, local_word = {}, []
        for child in order:
            index = leaves.index(child) + 1
            carrier = np.linalg.norm(root_column[list(remaining)])
            remaining.remove(index)
            cosine = np.linalg.norm(root_column[list(remaining)]) / carrier
            sine = root_column[index] / carrier
            local_coefficients[0, index] = (cosine, sine)
            replacement_coefficients[vertex, child] = (cosine, sine)
            local_word.append((0, index))
        replacement = _word(len(leaves) + 1, local_word, local_coefficients)
        local_correction = replacement.T @ original
        np.testing.assert_allclose(local_correction[:, 0], np.eye(len(leaves) + 1)[:, 0],
                                   atol=ATOL, rtol=0)
        correction[np.ix_(leaves, leaves)] = local_correction[1:, 1:]
    return _suffix_order(height, exceptional), replacement_coefficients, correction


def _bounded_marks(height, bound, rng):
    component = list(range(1 << height))
    exits, cuts = {}, set()
    for left, right in _edges(height):
        parent = component[left]
        if exits.get(parent, 0) < bound and rng.random() < .6:
            cuts.add((left, right))
            exits[parent] = exits.get(parent, 0) + 1
        else:
            component[right] = parent
    return cuts


def _row_envelope(label, internal_edges):
    """Expand U^T e_label through the actual reverse word, never its values."""
    support = {label}
    for left, right in reversed(internal_edges):
        if left in support or right in support:
            support.update((left, right))
    return support


def _product(matrices, size):
    result = np.eye(size)
    for matrix in matrices:
        result = matrix @ result
    return result


def _assemble(certificate, coefficients):
    """Build rotations on actual dressed rows, independently of support metadata."""
    base = certificate.base
    size = 1 << base.height
    identity = np.eye(size)
    base_word = _word(size, base.internal_edges, coefficients)
    dressed = {}
    for edge in base.edges:
        if edge not in base.cuts:
            continue
        left, right = edge
        parent, child = base_word[left], identity[right]
        cosine, sine = coefficients[edge]
        dressed[edge] = (identity + (cosine - 1) *
                         (np.outer(parent, parent) + np.outer(child, child)) +
                         sine * (np.outer(child, parent) - np.outer(parent, child)))
    band_words = []
    for band in range(len(certificate.bands)):
        band_words.append(_product(
            (dressed[edge] for edge in certificate.internal_cuts
             if certificate.band_of[base.component_of[edge[0]]] == band), size))
    internal = _product(band_words, size)
    sandwiches, cut_words = [], []
    for depth, edges in enumerate(certificate.cut_groups):
        parent_word = _product((word for band, word in enumerate(band_words)
                                if certificate.band_depths[band] == depth), size)
        cut_word = _product((dressed[edge] for edge in edges), size)
        cut_words.append(cut_word)
        sandwiches.append(parent_word.T @ cut_word @ parent_word)
    result = base_word @ internal @ _product(sandwiches, size)
    return result, base_word, dressed, band_words, cut_words, sandwiches


class ExceptionalExitTests(unittest.TestCase):
    def _assert_supports(self, certificate, band_words, cut_words, bound):
        base = certificate.base
        size = 1 << base.height
        expected_ports = []
        for index, component in enumerate(base.components):
            root = component[0]
            support = {root}
            for left, right in base.cuts:
                if base.component_of[left] == index:
                    row = _row_envelope(left, base.internal_edges)
                    self.assertLessEqual(len(row), _degree(root, base.height) + 1)
                    self.assertTrue(row <= set(component))
                    support.update(row)
            self.assertLessEqual(len(support), 1 + bound * (_degree(root, base.height) + 1))
            expected_ports.append(tuple(sorted(support)))
        self.assertEqual(certificate.port_supports, tuple(expected_ports))
        expected_bands = [set().union(*(set(expected_ports[i]) for i in band))
                          for band in certificate.bands]
        self.assertEqual(certificate.band_supports,
                         tuple(tuple(sorted(support)) for support in expected_bands))
        self.assertEqual(certificate.square_budget,
                         sum(len(support) ** 2 for support in expected_bands))
        for word, support in zip(band_words, expected_bands):
            outside = sorted(set(range(size)) - support)
            np.testing.assert_allclose((word - np.eye(size))[outside], 0, atol=ATOL, rtol=0)
            np.testing.assert_allclose((word - np.eye(size))[:, outside], 0, atol=ATOL, rtol=0)
        for depth, (edges, word) in enumerate(zip(certificate.cut_groups, cut_words)):
            blocks = {}
            for left, right in edges:
                band = certificate.band_of[base.component_of[left]]
                blocks.setdefault(band, set(expected_bands[band])).add(right)
            supports = [blocks[band] for band in sorted(blocks)]
            self.assertEqual(certificate.cut_block_supports[depth],
                             tuple(tuple(sorted(support)) for support in supports))
            self.assertEqual(certificate.cut_square_budgets[depth],
                             sum(len(support) ** 2 for support in supports))
            for i, support in enumerate(supports):
                self.assertTrue(all(support.isdisjoint(other) for other in supports[:i]))
            self.assertLessEqual(certificate.cut_square_budgets[depth],
                                 (bound + 1) ** 2 * certificate.square_budget)
            union = set().union(*supports)
            outside = sorted(set(range(size)) - union)
            np.testing.assert_allclose((word - np.eye(size))[outside], 0, atol=ATOL, rtol=0)

    def test_complete_factorization_after_paid_star_transfer(self):
        for height in range(3, 7):
            for bands in (2, 4):
                with self.subTest(height=height, bands=bands):
                    rng = np.random.default_rng(2189 + 11 * height + bands)
                    coefficients = _coefficients(_edges(height), rng)
                    exceptional = _bounded_marks(height, 2, rng)
                    edges, transferred, correction = _transfer(height, coefficients, exceptional)
                    certificate = partition_exceptional_exits(height, edges, exceptional, bands, 2)
                    actual, base_word, dressed, band_words, cut_words, _ = _assemble(
                        certificate, transferred)
                    size = 1 << height
                    expected = _word(size, _edges(height), coefficients)
                    np.testing.assert_allclose(actual @ correction, expected, atol=ATOL, rtol=0)
                    chronological = _product((dressed[e] for e in edges if e in exceptional), size)
                    sorted_word = _product((dressed[e] for e in certificate.ordered_cuts), size)
                    np.testing.assert_allclose(sorted_word, chronological, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(base_word @ chronological,
                                               _word(size, edges, transferred), atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.T @ actual, np.eye(size), atol=ATOL, rtol=0)
                    self._assert_supports(certificate, band_words, cut_words, 2)
                    self.assertLessEqual(certificate.factor_count, 3 * bands)

    def test_exact_signed_corners_keep_conservative_coordinate_envelopes(self):
        corners = ((1., 0.), (0., 1.), (-1., 0.), (0., -1.))
        for height in range(3, 7):
            rng = np.random.default_rng(432 + height)
            exceptional = _bounded_marks(height, 2, rng)
            edges = _suffix_order(height, exceptional)
            coefficients = _coefficients(edges, rng)
            for index, edge in enumerate(edges):
                if index % 3 != 2:
                    coefficients[edge] = corners[index % len(corners)]
            certificate = partition_exceptional_exits(height, edges, exceptional, 3, 2)
            actual, _, _, band_words, cut_words, _ = _assemble(certificate, coefficients)
            np.testing.assert_allclose(actual, _word(1 << height, edges, coefficients),
                                       atol=ATOL, rtol=0)
            self._assert_supports(certificate, band_words, cut_words, 2)

    def test_distinct_ports_commute_despite_overlapping_tiny_row_coordinates(self):
        exceptional = {(4, 5), (2, 3)}
        edges = _suffix_order(3, exceptional)
        coefficients = {edge: (np.cos(.37), np.sin(.37)) for edge in edges}
        coefficients[0, 4] = (1., 1e-250)
        certificate = partition_exceptional_exits(3, edges, exceptional, 4, 2)
        _, base_word, dressed, _, _, _ = _assemble(certificate, coefficients)
        first = _row_envelope(4, certificate.base.internal_edges)
        second = _row_envelope(2, certificate.base.internal_edges)
        self.assertIn(4, second)  # An earlier sibling carrier, not an ancestor of 2.
        self.assertGreaterEqual(len(first & second), 2)
        self.assertGreater(abs(base_word[4, 0]), 0.)
        self.assertLess(abs(base_word[4, 0]), 1e-100)
        self.assertAlmostEqual(float(base_word[4] @ base_word[2]), 0., delta=ATOL)
        left, right = dressed[4, 5], dressed[2, 3]
        np.testing.assert_allclose(left @ right, right @ left, atol=ATOL, rtol=0)
        # A second exit at the SAME port is not free to commute.
        same_port = {(0, 4), (0, 2)}
        same_edges = _suffix_order(3, same_port)
        same = partition_exceptional_exits(3, same_edges, same_port, 4, 2)
        same_coefficients = {edge: (np.cos(.37), np.sin(.37)) for edge in same_edges}
        _, _, matrices, _, _, _ = _assemble(same, same_coefficients)
        self.assertGreater(np.linalg.norm(matrices[0, 4] @ matrices[0, 2] -
                                          matrices[0, 2] @ matrices[0, 4]), .01)

    def test_quotient_depth_factors_require_later_depths_on_the_left(self):
        height = 3
        exceptional = {(0, 4), (4, 6), (6, 7)}
        edges = _suffix_order(height, exceptional)
        coefficients = {edge: (np.cos(.43), np.sin(.43)) for edge in edges}
        certificate = partition_exceptional_exits(height, edges, exceptional, 4, 1)
        actual, base_word, _, band_words, _, sandwiches = _assemble(certificate, coefficients)
        self.assertEqual(len(sandwiches), 3)
        expected = _word(8, edges, coefficients)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        wrong = base_word @ _product(band_words, 8) @ _product(reversed(sandwiches), 8)
        self.assertGreater(np.linalg.norm(wrong - expected, 2), .1)

    def test_exact_dyadic_packing_repeated_degrees_and_long_chain(self):
        for height in range(3, 11):
            rng = np.random.default_rng(975 + height)
            examples = [_bounded_marks(height, 2, rng),
                        {(0, 1 << (height - 1)), (0, 1 << (height - 2))},
                        {(1 << (height - 1), (1 << (height - 1)) + 1),
                         (1 << (height - 2), (1 << (height - 2)) + 1)}]
            for exceptional in examples:
                certificate = partition_exceptional_exits(
                    height, _suffix_order(height, exceptional), exceptional, 4, 2)
                base = certificate.base
                for index, component in enumerate(base.components):
                    root, degree = component[0], _degree(component[0], height)
                    intervals = []
                    for left, right in exceptional:
                        if base.component_of[left] == index:
                            child_degree = _degree(right, height)
                            interval = set(range(right, right + (1 << child_degree)))
                            self.assertNotIn(root, interval)
                            self.assertTrue(interval <= set(range(root, root + (1 << degree))))
                            self.assertTrue(all(interval.isdisjoint(other) for other in intervals))
                            intervals.append(interval)
                    weight = sum((Fraction(len(interval), 1 << degree) for interval in intervals),
                                 Fraction())
                    self.assertLessEqual(weight, Fraction(3, 4))
            repeated = examples[-1]
            self.assertEqual(len({_degree(right, height) for _, right in repeated}), 1)
            vertex, chain = 0, set()
            for bit in reversed(range(height)):
                child = vertex + (1 << bit)
                chain.add((vertex, child))
                vertex = child
            certificate = partition_exceptional_exits(
                height, _suffix_order(height, chain), chain, 1, 1)
            self.assertEqual(max(certificate.base.component_depths), height)
            self.assertEqual(certificate.factor_count, 3)
            self.assertEqual(certificate.cut_groups, ())
            # One whole chain fits in a polynomial coordinate block, regardless
            # of the unbounded exceptional depth and branching inside its base.
            self.assertLessEqual(len(certificate.band_supports[0]),
                                 (height + 1) * (height + 4) // 2)
        # Exact strict geometric slack in the b=2, alpha=3/4, s=4 certificate.
        self.assertLess(Fraction(27, 32), 1)  # r_b**4
        self.assertLess(Fraction(3, 4) * Fraction(5, 4), 1)  # log_2(2 rho)

    def test_rejects_wrong_stable_order_and_component_exit_count(self):
        exceptional = {(4, 5), (2, 3)}
        edges = _suffix_order(3, exceptional)
        with self.assertRaisesRegex(ValueError, "exceeds.*exit bound"):
            partition_exceptional_exits(3, edges, exceptional, 4, 1)
        for value in (0, -1, True, 1.5):
            with self.subTest(value=value), self.assertRaises(ValueError):
                partition_exceptional_exits(3, edges, exceptional, 4, value)
        for cuts in (set(), {(0, 4), (0, 2)}):
            bad = list(_suffix_order(3, cuts))
            first, second = bad.index((0, 4)), bad.index((0, 2))
            bad[first], bad[second] = bad[second], bad[first]
            with self.assertRaisesRegex(ValueError, "retain decreasing degree"):
                partition_exceptional_exits(3, bad, cuts, 4, 2)


if __name__ == "__main__":
    unittest.main()
