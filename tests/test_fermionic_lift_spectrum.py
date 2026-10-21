"""Exact finite checks for full-occupation errors of a fermionic lift.

The phase criterion and sparse-column obstruction have analytic proofs.
These fixtures challenge phases, cross terms, and arbitrary helper return.
"""
from fractions import Fraction
from itertools import combinations, permutations
import unittest

from tests.exact_matrices import (
    CO, CZ, CQ2, Q2, adjoint, eye, multiply, native_columns, zeros,
)
from tests.test_bank_motion_walsh_indicator import _inverse


def determinant(matrix):
    size = len(matrix)
    result = CZ
    for order in permutations(range(size)):
        inversions = sum(order[i] > order[j]
                         for i in range(size) for j in range(i + 1, size))
        term = CO
        for row, column in enumerate(order):
            term *= matrix[row][column]
        result += -term if inversions & 1 else term
    return result


def exterior_lift(matrix):
    """Literal occupation basis, including vacuum and determinant sector."""
    size = len(matrix)
    result = zeros(1 << size)
    for weight in range(size + 1):
        subsets = list(combinations(range(size), weight))
        for rows in subsets:
            row = sum(1 << bit for bit in rows)
            for columns in subsets:
                column = sum(1 << bit for bit in columns)
                result[row][column] = determinant(
                    [[matrix[i][j] for j in columns] for i in rows])
    return result


def subtract(left, right):
    return [[a - b for a, b in zip(row, other)]
            for row, other in zip(left, right)]


def squared_norm(column):
    return sum((value.conj() * value for value in column), CZ)


class FermionicLiftSpectrumTests(unittest.TestCase):
    def test_determinant_one_small_phases_accumulate_on_subset(self):
        root_half = Q2(0, Fraction(1, 2))
        omega = CQ2(root_half, root_half)
        phases = [omega] * 4 + [omega.conj()] * 4
        # Every one-particle squared error is 2-sqrt(2), less than one.
        one_error = CQ2(Q2(2, -1))
        for value in phases:
            self.assertEqual((value - CO).conj() * (value - CO), one_error)
        values = []
        for mask in range(256):
            value = CO
            for bit, phase in enumerate(phases):
                if mask >> bit & 1:
                    value *= phase
            values.append(value)
        self.assertEqual(values[0], CO)       # Vacuum fixes literal phase.
        self.assertEqual(values[-1], CO)      # Determinant exactly one.
        self.assertEqual(values[15], -CO)     # All four positive phases.
        self.assertEqual((values[15] - CO).conj() * (values[15] - CO),
                         CQ2(Q2(4)))

    def test_tilted_rank_two_relative_error_has_no_dimension_factor(self):
        # Reflection about (cos(pi/8),sin(pi/8),0,0). Its entries lie in
        # Q(sqrt(2)), although the vector coordinates need a larger field.
        root_half = CQ2(Q2(0, Fraction(1, 2)))
        target, tilted = eye(4), eye(4)
        target[0][0] = -CO
        tilted[0][0], tilted[0][1] = -root_half, -root_half
        tilted[1][0], tilted[1][1] = -root_half, root_half
        relative = multiply(target, tilted)
        lifted = exterior_lift(relative)
        self.assertEqual(lifted, multiply(exterior_lift(target),
                                         exterior_lift(tilted)))
        difference = subtract(lifted, eye(16))
        square = multiply(adjoint(difference), difference)
        error_squared = CQ2(Q2(2, -1))
        for row in range(16):
            for column in range(16):
                expected = (error_squared if row == column and
                            ((column & 1) ^ ((column >> 1) & 1)) else CZ)
                self.assertEqual(square[row][column], expected)
        # Occupying neither or both active modes gives identity, including
        # arbitrary spectator occupation and the full determinant sector.
        self.assertEqual(lifted[15][15], CO)

    def test_native_two_mode_hadamard_retains_determinant_sign(self):
        root_half = CQ2(Q2(0, Fraction(1, 2)))
        single = [[root_half, root_half], [root_half, -root_half]]
        align = [('CX', 0, 1), ('H', 0), ('CX', 0, 1), ('S', 1)]
        word = ([('S', 1), ('S', 1)] + _inverse(align)
                + [('T', 0), ('TDG', 1)] + align)
        expected = exterior_lift(single)
        self.assertEqual(native_columns(2, word), expected)
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(expected))
        self.assertEqual(expected[0][0], CO)
        self.assertEqual(expected[3][3], -CO)

    def test_auxiliary_fiber_parity_decode_returns_all_dirty_columns(self):
        # Two rows, two modes per row. The second mode is arbitrary work.
        decode = [('CX', 1, 0), ('CX', 3, 2)]
        for address in range(2):
            row = 2 * address
            mark = [('S', row), ('S', row),
                    ('S', row + 1), ('S', row + 1)]
            actual = native_columns(4, decode + mark + _inverse(decode))
            for column in range(16):
                for output in range(16):
                    expected = ((-CO if column >> row & 1 else CO)
                                if output == column else CZ)
                    self.assertEqual(actual[output][column], expected)

    def test_sparse_uniform_eigenvector_cross_term_exact(self):
        # U is minus a four-cycle on data and identity on two auxiliary
        # labels. It fixes each uniform-data auxiliary vector with sign -.
        m, auxiliary = 4, 2
        size = m * auxiliary
        unitary, target = zeros(size), eye(size)
        for a in range(auxiliary):
            for u in range(m):
                unitary[((u + 1) % m) * auxiliary + a][u * auxiliary + a] = -CO
                for v in range(m):
                    target[v * auxiliary + a][u * auxiliary + a] -= CQ2(Q2(Fraction(2, m)))
        first = subtract(unitary, eye(size))
        error = subtract(unitary, target)
        frobenius_squared = CZ
        for column in range(size):
            a = column % auxiliary
            self.assertEqual(sum((first[v * auxiliary + a][column]
                                  for v in range(m)), CZ), CQ2(Q2(-2)))
            first_norm = squared_norm([row[column] for row in first])
            error_norm = squared_norm([row[column] for row in error])
            self.assertEqual(error_norm, first_norm - CQ2(Q2(Fraction(4, m))))
            self.assertEqual(error_norm, CO)
            frobenius_squared += error_norm
        # Column sparsity s=1 attains 4d(m/(s+1)-1)=8 exactly.
        self.assertEqual(frobenius_squared, CQ2(Q2(8)))


if __name__ == '__main__':
    unittest.main()
