"""Four exact carrier fixtures for the involutory Gauss identities.

These checks retain literal scalar phases, every carrier column, and the
actual native inverse. They do not search parameter families or certify a
native emitter for an externally controlled matrix.
"""
from fractions import Fraction
from itertools import permutations
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, Q2, adjoint, eye, multiply, native_columns, scale, zeros,
)


def _apply(matrix, value):
    return sum((sum(bit * (value >> j & 1) for j, bit in enumerate(row)) & 1)
               << i for i, row in enumerate(matrix))


def _quadratic(matrix, value):
    n = len(matrix)
    return (sum(matrix[i][i] * (value >> i & 1) for i in range(n))
            + 2 * sum(matrix[i][j] * (value >> i & 1) * (value >> j & 1)
                      for i in range(n) for j in range(i + 1, n))) % 4


def _determinant(matrix):
    total = 0
    for perm in permutations(range(len(matrix))):
        inversions = sum(perm[i] > perm[j] for i in range(len(perm))
                         for j in range(i + 1, len(perm)))
        term = (-1) ** inversions
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        total += term
    return total


def _fourier_word(matrix):
    word = []
    for i, row in enumerate(matrix):
        if row[i]:
            word.append(('S', i))
        for j in range(i + 1, len(matrix)):
            if row[j]:
                word.extend([('H', j), ('CX', i, j), ('H', j)])
    return word + [('H', i) for i in range(len(matrix))]


def _actual_inverse(word):
    inverse = []
    for gate in reversed(word):
        inverse.extend([gate] * (3 if gate[0] == 'S' else 1))
    return inverse


class InvolutoryGaussSourceTests(unittest.TestCase):
    def test_literal_loops_inverse_and_characteristic_phase(self):
        # K4 has determinant -3 and ell=(1,1,1,1). The final example
        # is I+[[B,B],[B,B]] with B=diag(1,0), whose sixth power is -I.
        fixtures = [
            [[1, 0], [0, 1]],
            [[0, 1], [1, 0]],
            [[int(i != j) for j in range(4)] for i in range(4)],
            [[0, 0, 1, 0], [0, 1, 0, 0],
             [1, 0, 0, 0], [0, 0, 0, 1]],
        ]
        i_powers = [CO, CI, -CO, -CI]
        half_root = Q2(0, Fraction(1, 2))
        omega = CQ2(half_root, half_root)
        omega_powers = [CO]
        for _ in range(7):
            omega_powers.append(omega_powers[-1] * omega)

        for matrix in fixtures:
            n = len(matrix)
            size = 1 << n
            with self.subTest(matrix=matrix):
                self.assertTrue(all(matrix[i][j] == matrix[j][i]
                                    for i in range(n) for j in range(n)))
                self.assertEqual([_apply(matrix, _apply(matrix, 1 << i))
                                  for i in range(n)], [1 << i for i in range(n)])
                gamma = sum((i_powers[_quadratic(matrix, z)]
                             for z in range(size)), CQ2()) / (1 << (n // 2))
                diag = sum(matrix[i][i] << i for i in range(n))
                self.assertEqual(_apply(matrix, diag), diag)
                determinant = _determinant(matrix)
                characteristic_norm = sum(matrix[i][j] * (diag >> i & 1)
                                          * (diag >> j & 1)
                                          for i in range(n) for j in range(n))
                chi = (-1) ** ((determinant * determinant - 1) // 8)
                self.assertEqual(gamma, chi * omega_powers[characteristic_norm % 8])

                ell = sum((((_quadratic(matrix, _apply(matrix, 1 << i))
                             - matrix[i][i]) // 2) & 1) << i for i in range(n))
                self.assertEqual(_apply(matrix, ell), ell)
                self.assertEqual(ell.bit_count() % 2, 0)
                expected_cube = zeros(size)
                for column in range(size):
                    row = _apply(matrix, column) ^ ell
                    expected_cube[row][column] = gamma * (-1) ** ((ell & column).bit_count())

                cube_word = _fourier_word(matrix) * 3
                cube = native_columns(n, cube_word)
                inverse = native_columns(n, _actual_inverse(cube_word))
                self.assertEqual(cube, expected_cube)
                self.assertEqual(inverse, adjoint(cube))
                self.assertEqual(multiply(inverse, cube), eye(size))
                sixth = multiply(cube, cube)
                self.assertEqual(sixth, scale(eye(size), gamma * gamma))
                self.assertEqual(sixth, scale(eye(size), i_powers[_quadratic(matrix, diag)]))

        doubled_gamma = i_powers[(2 + 1 + 2) % 4]
        self.assertEqual(doubled_gamma, CI)
        self.assertEqual(gamma, doubled_gamma)
        self.assertEqual(sixth, scale(eye(16), -CO))


if __name__ == '__main__':
    unittest.main()
