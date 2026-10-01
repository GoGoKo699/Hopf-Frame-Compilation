"""Small product-first identities and their logical-support boundary.

These compare complete matrices of dimension at most sixteen. They check
classical product coordinates, literal phases and address conventions; the
retained one-target compiler supplies the native resource theorem.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 3e-12
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULIS = (X, Y, Z)


def _quaternion_product(second, first):
    w2, v2 = second[0], second[1:]
    w1, v1 = first[0], first[1:]
    return np.r_[w2 * w1 - np.dot(v2, v1),
                 w2 * v1 + w1 * v2 + np.cross(v2, v1)]


def _quaternion_matrix(value):
    return value[0] * I2 - 1j * sum(value[j + 1] * PAULIS[j] for j in range(3))


def _axis_rotation(axis, angle):
    return np.cos(angle) * I2 - 1j * np.sin(angle) * sum(axis[j] * PAULIS[j] for j in range(3))


def _multiplexor(rows):
    result = np.zeros((2 * len(rows), 2 * len(rows)), dtype=complex)
    for index, block in enumerate(rows):
        result[2 * index:2 * index + 2, 2 * index:2 * index + 2] = block
    return result


def _pair_rotation(size, first, second, angle):
    result = np.eye(size, dtype=complex)
    result[np.ix_([first, second], [first, second])] = [
        [np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    return result


class SmallProductCompilationTests(unittest.TestCase):
    def test_two_and_three_noncommuting_rotations_keep_literal_quaternion_phase(self):
        for alpha, beta, gamma in ((.31, -.47, .23), (-.62, .19, -.41),
                                   (np.pi, 0, .27), (0, np.pi / 2, -np.pi / 2)):
            with self.subTest(alpha=alpha, beta=beta, gamma=gamma):
                ca, sa = np.cos(alpha), np.sin(alpha)
                cb, sb = np.cos(beta), np.sin(beta)
                cg, sg = np.cos(gamma), np.sin(gamma)
                rx = _axis_rotation([1, 0, 0], alpha)
                ry = _axis_rotation([0, 1, 0], beta)
                rz = _axis_rotation([0, 0, 1], gamma)
                two = np.array([ca * cb, sa * cb, ca * sb, sa * sb])
                three = np.array([cg * ca * cb - sg * sa * sb,
                                  cg * sa * cb - sg * ca * sb,
                                  cg * ca * sb + sg * sa * cb,
                                  cg * sa * sb + sg * ca * cb])
                np.testing.assert_allclose(_quaternion_matrix(two), rx @ ry, atol=ATOL, rtol=0)
                np.testing.assert_allclose(_quaternion_matrix(three), rz @ rx @ ry,
                                           atol=ATOL, rtol=0)
                recursive = _quaternion_product(np.array([cg, 0, 0, sg]), two)
                np.testing.assert_allclose(recursive, three, atol=ATOL, rtol=0)
                self.assertAlmostEqual(np.dot(three, three), 1, delta=ATOL)
        np.testing.assert_allclose(_quaternion_matrix(np.array([-1, 0, 0, 0])),
                                   -I2, atol=ATOL, rtol=0)

    def test_variable_length_words_share_one_unchanged_address_table(self):
        rng = np.random.default_rng(1307)
        lengths = (0, 1, 2, 3, 5, 8, 13, 21)
        words = []
        for row, length in enumerate(lengths):
            factors = []
            for index in range(length):
                axis = rng.normal(size=3)
                axis /= np.linalg.norm(axis)
                angle = rng.uniform(-np.pi, np.pi)
                if row == 1:
                    axis, angle = np.array([1.0, 0, 0]), np.pi
                factors.append((axis, angle))
            words.append(factors)
        # Actual stage products on all three address bits and the target.
        actual = np.eye(16, dtype=complex)
        for step in range(max(lengths)):
            rows = [_axis_rotation(*word[step]) if step < len(word) else I2 for word in words]
            actual = _multiplexor(rows) @ actual
        products = []
        for word in words:
            accumulated = np.array([1.0, 0, 0, 0])
            for axis, angle in word:
                factor = np.r_[np.cos(angle), np.sin(angle) * axis]
                accumulated = _quaternion_product(factor, accumulated)
            self.assertAlmostEqual(np.dot(accumulated, accumulated), 1, delta=ATOL)
            products.append(_quaternion_matrix(accumulated))
        compiled_table = _multiplexor(products)
        np.testing.assert_allclose(compiled_table, actual, atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual.conj().T @ actual, np.eye(16), atol=ATOL, rtol=0)
        np.testing.assert_allclose(compiled_table[:2, :2], I2, atol=ATOL, rtol=0)
        np.testing.assert_allclose(compiled_table[2:4, 2:4], -I2, atol=ATOL, rtol=0)

    def test_changing_the_address_invalidates_unchanged_row_preprocessing(self):
        first = [_axis_rotation([1, 0, 0], .37), _axis_rotation([0, 1, 0], -.41)]
        second = [_axis_rotation([0, 0, 1], .53), _axis_rotation([1, 0, 0], -.29)]
        a, b = _multiplexor(first), _multiplexor(second)
        hadamard = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        changed_address = np.kron(hadamard, I2)
        actual = b @ changed_address @ a
        unchanged_rows = _multiplexor([second[row] @ first[row] for row in range(2)])
        incorrect = changed_address @ unchanged_rows
        self.assertGreater(np.linalg.norm(actual - incorrect, 2), .3)
        # Even after removing the visible address change, an off-diagonal
        # address block remains. It cannot be a single unchanged-row table.
        residual = changed_address @ actual
        np.testing.assert_allclose(residual[:2, 2:],
                                   (second[0] - second[1]) @ first[1] / 2,
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(residual[:2, 2:], 2), .2)
        # A classical permutation is different: it permits an explicit row
        # relabeling, but still not the original unchanged-row product.
        permutation = np.kron(X, I2)
        relabeled = _multiplexor([second[row ^ 1] @ first[row] for row in range(2)])
        np.testing.assert_allclose(b @ permutation @ a, permutation @ relabeled,
                                   atol=ATOL, rtol=0)

    def test_overlapping_hopf_pairs_create_the_three_mode_mixed_path(self):
        # Existing two-qubit modes: root (0,2), right child (2,3), spectator 1.
        modes = [0, 2, 3, 1]
        for alpha, beta in ((.37, -.61), (-.23, .48), (np.pi / 2, np.pi / 2)):
            with self.subTest(alpha=alpha, beta=beta):
                first = _pair_rotation(4, 0, 2, alpha)
                second = _pair_rotation(4, 2, 3, beta)
                actual = second @ first
                ca, sa, cb, sb = np.cos(alpha), np.sin(alpha), np.cos(beta), np.sin(beta)
                expected = np.array([[ca, -sa, 0, 0],
                                     [cb * sa, cb * ca, -sb, 0],
                                     [sb * sa, sb * ca, cb, 0],
                                     [0, 0, 0, 1]], dtype=complex)
                np.testing.assert_allclose(actual[np.ix_(modes, modes)], expected,
                                           atol=ATOL, rtol=0)
                self.assertAlmostEqual(actual[3, 0], sb * sa, delta=ATOL)
                self.assertGreater(abs(actual[3, 0]), .09)
                # Reversing the two factors has no path from 0 to 3.
                self.assertAlmostEqual((first @ second)[3, 0], 0, delta=ATOL)
                np.testing.assert_allclose(actual[:, 1], np.eye(4)[:, 1], atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
