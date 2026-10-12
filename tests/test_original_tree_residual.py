"""Independent finite checks for original ordered trees, not native synthesis.

The exact tests cover arbitrary application orders rather than the canonical
recursive chart. Complex coarse blocks are literal native-word matrices, but
the rectangular residual factors here neither emit sources nor return work.
"""
from fractions import Fraction
from itertools import combinations, permutations
import unittest

import numpy as np

from compiler_robust_hopf.frames import direct_real_frame, hopf_ry
from compiler_robust_hopf.original_tree_residual import (
    cayley_matrix, hierarchical_shared_factors, hopf_tree_edges,
    stable_residual, tree_cayley_from_paths, tree_edge_cuts, tree_word_matrix,
)
from tests.test_coupled_residual_merge import _native_small_word
from tests.test_operator_source_compiler import H, X, Z
from tests.test_tree_residual_structure import S


ATOL = 2e-10


def _dense_word(size, edges, blocks):
    """Independent full-matrix embedding and multiplication oracle."""
    result = np.eye(size, dtype=complex)
    for edge, block in zip(edges, blocks, strict=True):
        local = np.eye(size, dtype=complex)
        local[np.ix_(edge, edge)] = block
        result = local @ result
    return result


def _fraction_word(size, edges, tangents):
    result = np.eye(size, dtype=object)
    for (first, second), tangent in zip(edges, tangents, strict=True):
        cosine = (1 - tangent * tangent) / (1 + tangent * tangent)
        sine = 2 * tangent / (1 + tangent * tangent)
        local = np.eye(size, dtype=object)
        local[first, first] = local[second, second] = cosine
        local[first, second], local[second, first] = -sine, sine
        result = local @ result
    return result


def _fraction_determinant(matrix):
    rows = [[Fraction(entry) for entry in row] for row in matrix]
    answer = Fraction(1)
    for column in range(len(rows)):
        pivot = next((j for j in range(column, len(rows)) if rows[j][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        value = rows[column][column]
        answer *= value
        for j in range(column + 1, len(rows)):
            multiplier = rows[j][column] / value
            for k in range(column + 1, len(rows)):
                rows[j][k] -= multiplier * rows[column][k]
    return answer


def _complex_coarse_fixture(n):
    native = _native_small_word()
    perturbations = (native, S @ native @ S.conj().T, H @ native.conj().T @ H)
    bases = (np.eye(2), H @ Z, X @ Z, Z @ H)
    target, coarse = [], []
    for index in range((1 << n) - 1):
        base = bases[(index + 1) % len(bases)]
        target.append(base @ hopf_ry((index % 5 - 2) / (10000 * n)))
        coarse.append(base @ perturbations[index % len(perturbations)])
    return target, coarse


class OriginalTreeResidualTests(unittest.TestCase):
    def test_hopf_label_and_rotation_conventions_match_independent_complete_frame(self):
        rng = np.random.default_rng(419)
        for n in (1, 2, 3, 4):
            angles = rng.uniform(-2, 2, (1 << n) - 1)
            angles[::3] = 0
            actual = tree_word_matrix(1 << n, hopf_tree_edges(n), [hopf_ry(a) for a in angles])
            np.testing.assert_allclose(actual, direct_real_frame(n, angles), atol=ATOL, rtol=0)

    def test_exact_signed_path_formula_and_determinant_for_every_tree_order(self):
        base_edges = ((0, 1), (2, 1), (1, 3), (4, 3))
        for values in ((Fraction(2, 3), Fraction(-1, 5), Fraction(3, 7), Fraction(-2, 9)),
                       (Fraction(0), Fraction(-1, 2), Fraction(2, 5), Fraction(0))):
            for order in permutations(range(4)):
                with self.subTest(tangents=values, order=order):
                    edges = tuple(base_edges[j] for j in order)
                    tangents = tuple(values[j] for j in order)
                    word = _fraction_word(5, edges, tangents)
                    generator = tree_cayley_from_paths(5, edges, tangents)
                    identity = np.eye(5, dtype=object)
                    np.testing.assert_array_equal(generator @ (word + identity), word - identity)
                    np.testing.assert_array_equal(generator.T, -generator)
                    expected = Fraction(1 << 5)
                    for tangent in tangents:
                        expected /= 1 + tangent * tangent
                    self.assertEqual(_fraction_determinant(word + identity), expected)
                    for left, right in tree_edge_cuts(5, edges):
                        for first, second in combinations(left, 2):
                            for third, fourth in combinations(right, 2):
                                self.assertEqual(generator[first, third] * generator[second, fourth],
                                                 generator[first, fourth] * generator[second, third])
        # Changing application time really changes signs, even for one path.
        forward = tree_cayley_from_paths(3, ((0, 1), (1, 2)), (2., 3.))
        backward = tree_cayley_from_paths(3, ((1, 2), (0, 1)), (3., 2.))
        self.assertEqual(forward[0, 2], -6)
        self.assertEqual(backward[0, 2], 6)

    def test_random_ordered_trees_match_dense_cayley_and_joining_update(self):
        rng = np.random.default_rng(1842)
        for size in (2, 5, 8, 12):
            for trial in range(5):
                edges = [(int(rng.integers(vertex)), vertex) for vertex in range(1, size)]
                rng.shuffle(edges)
                edges = [edge if index % 2 else edge[::-1] for index, edge in enumerate(edges)]
                angles = rng.uniform(-1.3, 1.3, size - 1)
                angles[trial % (size - 1)] = 0
                tangents = np.tan(angles / 2)
                word = _dense_word(size, edges, [hopf_ry(a) for a in angles])
                path = tree_cayley_from_paths(size, edges, tangents)
                np.testing.assert_allclose(path, cayley_matrix(word), atol=ATOL, rtol=0)
                updated, identity = np.zeros((size, size)), np.eye(size)
                for (first, second), tangent in zip(edges, tangents, strict=True):
                    k = np.zeros((size, size))
                    k[first, second], k[second, first] = -1, 1
                    np.testing.assert_allclose(k @ updated @ k, 0, atol=ATOL, rtol=0)
                    updated += tangent * (identity - updated) @ k @ (identity + updated)
                np.testing.assert_allclose(updated, path, atol=ATOL, rtol=0)
                expected = 2. ** size * np.prod(np.cos(angles / 2) ** 2)
                self.assertAlmostEqual(float(np.linalg.det(identity + word).real) / expected, 1, places=12)

    def test_singular_corner_is_a_chart_boundary_not_a_removed_logical_column(self):
        edges = hopf_tree_edges(3)
        for position in (0, 3, 6):
            blocks = [np.eye(2) for _ in edges]
            blocks[position] = -np.eye(2)
            word = tree_word_matrix(8, edges, blocks)
            self.assertEqual(np.linalg.det(np.eye(8) + word), 0)
            np.testing.assert_allclose(word.conj().T @ word, np.eye(8), atol=ATOL, rtol=0)
            with self.assertRaisesRegex(ValueError, "singular or unresolved"):
                cayley_matrix(word)
            for offset in (-1e-3, 1e-3):
                angles = np.zeros(7)
                angles[position] = np.pi + offset
                near = _dense_word(8, edges, [hopf_ry(a) for a in angles])
                np.testing.assert_allclose(cayley_matrix(near),
                                           tree_cayley_from_paths(8, edges, np.tan(angles / 2)),
                                           atol=2e-7, rtol=0)
                expected = 256 * np.prod(np.cos(angles / 2) ** 2)
                self.assertAlmostEqual(np.linalg.det(np.eye(8) + near).real / expected, 1, places=10)

    def test_homogeneous_root_cut_and_conditioning_bounds_include_negative_tangents(self):
        for n in (1, 2, 3, 4, 5):
            size = 1 << n
            for tangent in (-.7, 0., np.sqrt(2) - 1, 1.):
                generator = tree_cayley_from_paths(size, hopf_tree_edges(n), [tangent] * (size - 1))
                cut_norm = np.linalg.norm(generator[:size // 2, size // 2:], 2)
                lower = abs(tangent) * (1 + tangent * tangent) ** (n - 1)
                self.assertAlmostEqual(cut_norm, lower, delta=ATOL)
                self.assertAlmostEqual(np.linalg.norm(generator[:, 0]) ** 2,
                                       (1 + tangent * tangent) ** n - 1, delta=ATOL)
                norm = np.linalg.norm(generator, 2)
                upper = ((1 + tangent * tangent) ** n - 1) / abs(tangent) if tangent else 0
                self.assertLessEqual(lower, norm + ATOL)
                self.assertLessEqual(norm, upper + ATOL)
                for sign in (-1, 1):
                    condition = np.linalg.cond(np.eye(size) + sign * generator)
                    self.assertGreaterEqual(condition + ATOL, np.sqrt(1 + lower ** 2) / (1 + tangent ** 2))
                    self.assertLessEqual(condition, np.sqrt(1 + upper ** 2) + ATOL)
                sign, logdet = np.linalg.slogdet(np.eye(size) - generator)
                self.assertEqual(sign, 1)
                self.assertAlmostEqual(logdet, (size - 1) * np.log1p(tangent * tangent), delta=ATOL)

    def test_outward_complex_native_coarse_residual_has_rank_one_on_every_edge(self):
        for n in (3, 4, 5):
            size, edges = 1 << n, hopf_tree_edges(n)
            target_blocks, coarse_blocks = _complex_coarse_fixture(n)
            for block in coarse_blocks:
                self.assertGreater(np.max(np.abs(block.imag)), 1e-6)
                np.testing.assert_allclose(np.linalg.det(block), 1, atol=ATOL, rtol=0)
            target = _dense_word(size, edges, target_blocks)
            coarse = tree_word_matrix(size, edges, coarse_blocks)
            report = stable_residual(target, coarse)
            self.assertLess(report.delta, 1 / (8 * np.sqrt(size)))
            self.assertGreater(np.max(np.abs(np.diag(report.residual).imag)), 1e-6)
            np.testing.assert_allclose(report.residual @ coarse, target, atol=ATOL, rtol=0)
            np.testing.assert_allclose(report.cayley + report.cayley.conj().T, 0, atol=ATOL, rtol=0)
            self.assertLessEqual(np.linalg.norm(report.cayley, 2), report.cayley_norm_bound + ATOL)
            for axis in (0, 1):
                self.assertLessEqual(np.max(np.sum(np.abs(report.residual), axis=axis)),
                                     report.absolute_sum_bound + ATOL)
                self.assertLess(report.absolute_sum_bound, 9 / 8)
            for left, right in tree_edge_cuts(size, edges):
                for rows, columns in ((left, right), (right, left)):
                    for matrix, bound in ((target, 1), (coarse, 1),
                                          (report.residual, 1), (report.cayley, 1)):
                        rank = np.linalg.matrix_rank(matrix[np.ix_(rows, columns)], tol=ATOL)
                        self.assertLessEqual(rank, bound)

    def test_reordered_complex_coarse_residual_attains_rank_two(self):
        # A joining edge with earlier nontrivial words on both sides makes
        # the rank-two upper bound sharp. Outward Hopf order need not do so.
        edges = ((0, 1), (2, 3), (1, 2))
        target_blocks, coarse_blocks = _complex_coarse_fixture(2)
        report = stable_residual(_dense_word(4, edges, target_blocks),
                                 _dense_word(4, edges, coarse_blocks))
        for matrix in (report.residual, report.cayley):
            self.assertEqual(np.linalg.matrix_rank(matrix[:2, 2:], tol=ATOL), 2)
            self.assertEqual(np.linalg.matrix_rank(matrix[2:, :2], tol=ATOL), 2)
        with self.assertRaisesRegex(ValueError, "rank greater than 1"):
            hierarchical_shared_factors(report.residual, max_cut_rank=1)
        self.assertLess(hierarchical_shared_factors(report.residual).reconstruction_error, ATOL)

    def test_shared_modes_preserve_complex_diagonal_phase_and_expose_completion_budget(self):
        for n in (3, 4):
            size, edges = 1 << n, hopf_tree_edges(n)
            target_blocks, coarse_blocks = _complex_coarse_fixture(n)
            residual = stable_residual(_dense_word(size, edges, target_blocks),
                                       _dense_word(size, edges, coarse_blocks)).residual
            factors = hierarchical_shared_factors(residual)
            self.assertEqual(factors.shared_modes, 4 * size - 4)
            self.assertEqual(factors.separate_failure_completion_modes, 6 * size - 4)
            self.assertGreater(factors.separate_failure_completion_modes, 4 * size)
            self.assertLess(factors.reconstruction_error, ATOL)
            np.testing.assert_allclose(factors.a @ factors.b.conj().T, residual, atol=ATOL, rtol=0)
            np.testing.assert_allclose(factors.a[:, :size] @ factors.b[:, :size].conj().T,
                                       np.diag(np.diag(residual)), atol=ATOL, rtol=0)
            self.assertGreater(np.max(np.abs(factors.a[:, :size].imag)), 1e-6)
            delta = np.linalg.norm(residual - np.eye(size), 2)
            self.assertAlmostEqual(factors.gram_bound, 1 + n * delta, delta=ATOL)
            self.assertLess(factors.gram_bound, 2)
            for factor in (factors.a, factors.b):
                self.assertLessEqual(np.linalg.norm(factor @ factor.conj().T, 2), factors.gram_bound + ATOL)
                column = size
                for depth in range(n):
                    count = 2 * (1 << depth) * min(2, size >> (depth + 1))
                    level = factor[:, column:column + count]
                    self.assertLessEqual(np.linalg.norm(level @ level.conj().T, 2), delta + ATOL)
                    column += count
                self.assertEqual(column, factors.shared_modes)
            rank_one = hierarchical_shared_factors(residual, max_cut_rank=1)
            self.assertEqual(rank_one.shared_modes, 3 * size - 2)
            self.assertEqual(rank_one.separate_failure_completion_modes, 5 * size - 2)
            self.assertGreater(rank_one.separate_failure_completion_modes, 4 * size)
            self.assertLess(rank_one.reconstruction_error, ATOL)
            np.testing.assert_allclose(rank_one.a @ rank_one.b.conj().T, residual, atol=ATOL, rtol=0)
            self.assertAlmostEqual(rank_one.gram_bound, factors.gram_bound, delta=ATOL)
            for factor in (rank_one.a, rank_one.b):
                self.assertLessEqual(np.linalg.norm(factor @ factor.conj().T, 2), rank_one.gram_bound + ATOL)
        # A zero residual keeps the declared slot count without guessing rank.
        identity = hierarchical_shared_factors(np.eye(8))
        np.testing.assert_array_equal(identity.a[:, 8:], 0)
        self.assertEqual(identity.gram_bound, 1)

    def test_invalid_trees_blocks_charts_and_rank_claims_are_rejected(self):
        for dimension, edges in ((3, ((0, 1),)), (3, ((0, 1), (1, 0))),
                                 (3, ((0, 0), (1, 2))), (3, ((0, 1), (1, 3))),
                                 (3, ((0, 1), (1., 2))), (True, ())):
            with self.assertRaises(ValueError):
                tree_edge_cuts(dimension, edges)
        for tangents in ((np.inf,), (np.nan,), (1j,), (), ((1,),), ('0.3',)):
            with self.assertRaises(ValueError):
                tree_cayley_from_paths(2, ((0, 1),), tangents)
        for blocks in ((), (np.diag([1, 1j]),), (np.ones((2, 2)),), (np.eye(3),)):
            with self.assertRaises(ValueError):
                tree_word_matrix(2, ((0, 1),), blocks)
        for n in (0, 1.5, True):
            with self.assertRaises(ValueError):
                hopf_tree_edges(n)
        for matrix in (np.eye(3), np.ones((4, 4)), np.array([[np.nan]])):
            with self.assertRaises(ValueError):
                hierarchical_shared_factors(matrix)
        with self.assertRaisesRegex(ValueError, "rank greater than 2"):
            hierarchical_shared_factors(np.kron(np.kron(H, H), H))
        for rank in (0, 3, 1.0, True):
            with self.assertRaises(ValueError):
                hierarchical_shared_factors(np.eye(4), max_cut_rank=rank)
        with self.assertRaises(ValueError):
            hierarchical_shared_factors(np.eye(4), rank_atol=-1)
        with self.assertRaises(ValueError):
            cayley_matrix(np.eye(2), singular_atol=np.inf)
        with self.assertRaises(ValueError):
            stable_residual(np.eye(2), np.eye(4))
        with self.assertRaises(ValueError):
            stable_residual(np.eye(2), -np.eye(2))


if __name__ == '__main__':
    unittest.main()
