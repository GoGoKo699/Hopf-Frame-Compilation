"""Small complete-residual Cayley checks, with no native synthesis claim.

The exact four-mode example uses rational arithmetic. Remaining matrices
have at most sixteen logical modes. Complex coarse words retain their actual
native phases; no extra logical modes, clean work, or target oracle is
supplied by the recursive representation.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_coupled_residual_merge import (
    _close_real_fixture, _direct_sum, _root_word, _subtree_residual,
)
from tests.test_operator_source_compiler import H, X, Z, _rotation


ATOL = 6e-11


def _right_solve(numerator, denominator):
    return np.linalg.solve(denominator.T, numerator.T).T


def _cayley(unitary):
    identity = np.eye(len(unitary))
    return _right_solve(unitary - identity, unitary + identity)


def _inverse_cayley(generator):
    identity = np.eye(len(generator))
    return _right_solve(identity + generator, identity - generator)


def _complex_fixture(height):
    coarse, target = _close_real_fixture(height)
    # Each base is itself a literal real Clifford word. Multiplication by
    # the imported native commutator keeps C complex, including at a chart
    # boundary; the requested target remains real and close to C.
    bases = (np.eye(2), H @ Z, H @ Z, X @ Z)
    angles = (0, np.pi / 4, np.pi / 4, np.pi / 2)
    for node in coarse:
        which = node % 4
        coarse[node] = bases[which] @ coarse[node]
        offset = (node + 2) / (4000 * height) if which == 1 else 0
        if which == 2:
            offset = -(node + 2) / (4000 * height)
        target[node] = _rotation(angles[which] + offset)
    return coarse, target


def _recursive_generators(height, coarse, target):
    records = {}

    def visit(remaining, node):
        if remaining == 0:
            return np.zeros((1, 1), complex), np.ones(1, complex), 0j
        left, pl, bl = visit(remaining - 1, 2 * node)
        right, pr, br = visit(remaining - 1, 2 * node + 1)
        size, half = 1 << remaining, 1 << (remaining - 1)
        identity = np.eye(2)
        a = _direct_sum(left, right)
        b = np.diag([bl, br])
        tau = _cayley(target[node] @ coarse[node].conj().T)
        denominator = identity + b @ tau
        h = _right_solve(tau, denominator)
        # Columns are built from child summaries, not from dense inversion.
        z = np.zeros((size, 2), complex)
        z[:half, 0], z[half:, 1] = pl, pr
        root = _root_word(coarse[node], size)
        k = root.conj().T @ (a + z @ h @ z.conj().T) @ root
        u = coarse[node][:, 0]
        g = (identity + h @ (identity - b)) @ u
        p = root.conj().T @ np.concatenate([g[0] * pl, g[1] * pr])
        beta = u.conj() @ (b + (identity + b) @ h @ (identity - b)) @ u
        records[node] = dict(k=k, p=p, beta=beta, h=h, z=z, g=g, b=b,
                             denominator=denominator, remaining=remaining)
        return k, p, beta

    visit(height, 1)
    return records


def _fraction_product(a, b):
    return [[sum(x * y for x, y in zip(row, column))
             for column in zip(*b)] for row in a]


def _fraction_rotation(size, low, high, tangent):
    result = [[Fraction(i == j) for j in range(size)] for i in range(size)]
    cosine = (1 - tangent * tangent) / (1 + tangent * tangent)
    sine = 2 * tangent / (1 + tangent * tangent)
    result[low][low] = result[high][high] = cosine
    result[high][low], result[low][high] = sine, -sine
    return result


class TreeCayleyTests(unittest.TestCase):
    def test_four_mode_path_products_are_exact_rational_identities(self):
        for a, b, c in ((Fraction(1, 7), Fraction(1, 11), Fraction(1, 13)),
                        (Fraction(-2, 7), Fraction(3, 11), Fraction(-4, 13)),
                        (Fraction(0), Fraction(1, 3), Fraction(0))):
            with self.subTest(half_angles=(a, b, c)):
                frame = _fraction_product(
                    _fraction_product(_fraction_rotation(4, 2, 3, c),
                                      _fraction_rotation(4, 0, 1, b)),
                    _fraction_rotation(4, 0, 2, a))
                k = [[0, -b, -a, -a * c],
                     [b, 0, -a * b, -a * b * c],
                     [a, a * b, 0, -c],
                     [a * c, a * b * c, c, 0]]
                plus = [[value + int(i == j) for j, value in enumerate(row)]
                        for i, row in enumerate(frame)]
                minus = [[value - int(i == j) for j, value in enumerate(row)]
                         for i, row in enumerate(frame)]
                self.assertEqual(_fraction_product(k, plus), minus)
                self.assertGreater(np.linalg.svd(np.array(plus, float), compute_uv=False)[-1], 1)

    def test_eight_mode_real_generator_is_dense_but_has_rank_one_sibling_blocks(self):
        height = 3
        coarse = {node: np.eye(2) for node in range(1, 8)}
        target = {node: _rotation(2 * np.arctan(1 / (node + 6))) for node in coarse}
        records = _recursive_generators(height, coarse, target)
        for node, record in records.items():
            direct = _subtree_residual(record['remaining'], coarse, target, node)
            np.testing.assert_allclose(record['k'], _cayley(direct), atol=ATOL, rtol=0)
            half = len(direct) // 2
            self.assertEqual(np.linalg.matrix_rank(record['k'][:half, half:], tol=ATOL), 1)
        k = records[1]['k']
        self.assertEqual(np.count_nonzero(np.abs(np.triu(k, 1)) > ATOL), 28)
        self.assertEqual(np.linalg.matrix_rank(k, tol=ATOL), 8)
        np.testing.assert_allclose(records[1]['p'],
                                   [1, 1 / 10, 1 / 8, 1 / 88, 1 / 7, 1 / 84, 1 / 63, 1 / 819],
                                   atol=ATOL, rtol=0)

    def test_complex_native_coarse_frames_reconstruct_every_subtree_column(self):
        for height in (2, 3, 4):
            coarse, target = _complex_fixture(height)
            records = _recursive_generators(height, coarse, target)
            for node, record in records.items():
                self.assertGreater(np.max(np.abs(coarse[node].imag)), 1e-6)
                direct = _subtree_residual(record['remaining'], coarse, target, node)
                np.testing.assert_allclose(_inverse_cayley(record['k']), direct, atol=ATOL, rtol=0)
                np.testing.assert_allclose(record['k'] + record['k'].conj().T, 0, atol=ATOL, rtol=0)
                half = len(direct) // 2
                self.assertLessEqual(np.linalg.matrix_rank(record['k'][:half, half:], tol=ATOL), 2)
                # Nested redundant bases also control ancestor contributions
                # to sibling blocks of the complete root matrix.
                child_p = [records[child]['p'] if child in records else np.ones(1)
                           for child in (2 * node, 2 * node + 1)]
                bases = [np.column_stack([np.eye(half)[:, 0], column]) for column in child_p]
                child_basis = np.zeros((2 * half, 4), complex)
                child_basis[:half, :2], child_basis[half:, 2:] = bases
                correction = (coarse[node].conj().T - np.eye(2)) @ (np.eye(2) + record['b']) @ record['g']
                transfer = np.array([[1, correction[0]], [0, record['g'][0]],
                                     [0, correction[1]], [0, record['g'][1]]])
                np.testing.assert_allclose(child_basis @ transfer,
                                           np.column_stack([np.eye(len(direct))[:, 0], record['p']]),
                                           atol=ATOL, rtol=0)
                depth = node.bit_length() - 1
                start = (node - (1 << depth)) * len(direct)
                root_block = records[1]['k'][start:start + half, start + half:start + 2 * half]
                self.assertLessEqual(np.linalg.matrix_rank(root_block, tol=ATOL), 2)
            half = 1 << (height - 1)
            self.assertEqual(np.linalg.matrix_rank(records[1]['k'][:half, half:], tol=ATOL), 2)

    def test_scalar_and_root_column_recursions_survive_zero_and_degenerate_children(self):
        for height in (2, 3):
            for case in ('all_zero', 'left_zero', 'complex'):
                coarse, target = _complex_fixture(height)
                for node in coarse:
                    # A node descends from the left child iff its second
                    # leading binary bit is zero.
                    in_left = node > 1 and (node >> (node.bit_length() - 2)) == 2
                    if case == 'all_zero' or (case == 'left_zero' and in_left):
                        coarse[node] = target[node] = np.eye(2)
                records = _recursive_generators(height, coarse, target)
                for node, record in records.items():
                    direct = _cayley(_subtree_residual(record['remaining'], coarse, target, node))
                    np.testing.assert_allclose(record['p'], (np.eye(len(direct)) + direct)[:, 0],
                                               atol=ATOL, rtol=0)
                    self.assertAlmostEqual(abs(record['beta'] - direct[0, 0]), 0, delta=ATOL)
                    self.assertAlmostEqual(record['beta'].real, 0, delta=ATOL)
                if case == 'all_zero':
                    np.testing.assert_array_equal(records[1]['k'], np.zeros((1 << height, 1 << height)))

    def test_conditioning_uses_the_actual_uniform_subtree_error(self):
        coarse, target = _complex_fixture(3)
        records = _recursive_generators(3, coarse, target)
        # Measure every subtree; a local error bound alone is insufficient.
        errors = []
        for node, record in records.items():
            residual = _subtree_residual(record['remaining'], coarse, target, node)
            errors += [np.linalg.norm(residual - np.eye(len(residual)), 2),
                       np.linalg.norm(target[node] @ coarse[node].conj().T - np.eye(2), 2)]
        delta = max(errors)
        self.assertLess(delta, .01)
        kappa = delta / (2 - delta)
        for record in records.values():
            self.assertLessEqual(np.linalg.norm(record['k'], 2), kappa + ATOL)
            self.assertLessEqual(np.linalg.norm(np.linalg.inv(record['denominator']), 2),
                                 1 / (1 - kappa * kappa) + ATOL)
            self.assertLessEqual(np.linalg.norm(record['h'], 2), kappa / (1 - kappa * kappa) + ATOL)
            np.testing.assert_allclose(record['h'] + record['h'].conj().T, 0, atol=ATOL, rtol=0)
            singular = np.linalg.svd(record['z'], compute_uv=False)
            self.assertGreaterEqual(min(singular), 1 - ATOL)
            self.assertLessEqual(max(singular), np.sqrt(1 + kappa * kappa) + ATOL)

    def test_inverse_cayley_error_bound_preserves_complete_unitarity(self):
        rng = np.random.default_rng(602)
        for size in (4, 8):
            raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
            k = (raw - raw.conj().T) / 20
            noise = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
            perturbation = (noise - noise.conj().T) / 2000
            actual, ideal = _inverse_cayley(k + perturbation), _inverse_cayley(k)
            np.testing.assert_allclose(actual.conj().T @ actual, np.eye(size), atol=ATOL, rtol=0)
            np.testing.assert_allclose(_cayley(ideal), k, atol=ATOL, rtol=0)
            self.assertLessEqual(np.linalg.norm(actual - ideal, 2),
                                 2 * np.linalg.norm(perturbation, 2) + ATOL)


if __name__ == '__main__':
    unittest.main()
