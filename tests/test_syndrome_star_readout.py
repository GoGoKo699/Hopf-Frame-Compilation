"""Exact fixtures for a paid, but thin, syndrome-reflection acceptance.

The small native word checks the canonical two-bank sector and its literal
inverse. A separate exact matrix entry challenges the old kernel-X symmetry;
neither finite check supplies the all-size rank argument or a successful source.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, multiply, native_columns
from tests.test_two_layer_diagonal_rigidity import (
    _primitive_inverse, _reduce_adjacent_s,
)


def _inverse(word):
    return _reduce_adjacent_s(_primitive_inverse(word))


def _wide_t(bank):
    h = [('H', bit) for bit in bank]
    t = [('T', bit) for bit in bank]
    # F = T_X T_B T_X, retaining the literal T scalars.
    return h + t + h + t + h + t + h


def _x(bit):
    return [('H', bit), ('S', bit), ('S', bit), ('H', bit)]


def _p1_acceptance():
    # a=0; C=(1,2); D=(3,4); signal=5. B is D+C between the walls.
    wall = [('CX', 1, 3), ('CX', 2, 4)]
    f = _wide_t((3, 4))
    # Controlled R_S: phase -1 iff signal=1 and a+B_1=0.
    cr = ([('CX', 4, 0)] + _x(0)
          + [('H', 5), ('CX', 0, 5), ('H', 5)]
          + _x(0) + [('CX', 4, 0)])
    return (wall + [('H', 5)] + f + cr + _inverse(f) + _inverse(cr)
            + [('S', 5), ('S', 5), ('H', 5)] + list(reversed(wall)))


def _syndrome(bank, m):
    result = 0
    for row in range(m):
        if (bank >> row) & 1:
            result ^= row
    return result


def _entry(f, m, row, column):
    result = CO
    for bit in range(m):
        result *= f[(row >> bit) & 1][(column >> bit) & 1]
    return result


class SyndromeStarReadoutTests(unittest.TestCase):
    def test_native_two_bank_acceptance_all_columns_and_actual_inverse(self):
        word = _p1_acceptance()
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        f = native_columns(2, _wide_t((0, 1)))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 12)
        for a in range(2):
            signs = [(-1 if a == (b >> 1) else 1) for b in range(4)]
            rfr = [[signs[row] * value for value in entries]
                   for row, entries in enumerate(f)]
            u = multiply(adjoint(f), rfr)
            u = [[signs[row] * value for value in entries]
                 for row, entries in enumerate(u)]
            k = [[((CO if row == col else CZ) - u[row][col]) / 2
                  for col in range(4)] for row in range(4)]
            gram = multiply(adjoint(k), k)
            self.assertEqual(gram, [[
                ((2 * CO if row == col else CZ)
                 - u[row][col] - u[col][row].conj()) / 4
                for col in range(4)] for row in range(4)])
            for c in range(4):
                for d_in in range(4):
                    column = a | (c << 1) | (d_in << 3)
                    for d_out in range(4):
                        row = a | (c << 1) | (d_out << 3)
                        self.assertEqual(actual[row][column],
                                         k[d_out ^ c][d_in ^ c])
                    # Arbitrary C is a spectator after the complete two walls.
                    for row in range(32):
                        if (row & 7) != (column & 7):
                            self.assertEqual(actual[row][column], CZ)

    def test_wide_t_gram_really_breaks_a_kernel_translation_at_p3(self):
        # H(14)=1+2+3=0. The old R_S / bank-X pulse algebra must commute
        # with this translation. Two selected exact entries suffice to
        # disprove that premise for the changed accepted Gram at p=3.
        m, shift = 8, 14
        self.assertEqual(_syndrome(shift, m), 0)
        f = native_columns(1, _wide_t((0,)))
        centers = [b for b in range(1 << m) if _syndrome(b, m) == 0]

        def gram_entry(row, column):
            p_row = int(_syndrome(row, m) == 0)
            p_col = int(_syndrome(column, m) == 0)
            p_prime = sum((_entry(f, m, center, row).conj()
                           * _entry(f, m, center, column)
                           for center in centers), CZ)
            return ((CO if row == column and p_row else CZ)
                    + (1 - p_row - p_col) * p_prime)

        self.assertNotEqual(gram_entry(17, 119),
                            gram_entry(17 ^ shift, 119 ^ shift))


if __name__ == '__main__':
    unittest.main()
