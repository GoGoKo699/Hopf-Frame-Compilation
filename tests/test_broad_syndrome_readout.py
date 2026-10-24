"""Exact checks for a normalized broad-baseline syndrome source.

The native source is fully paid but does not recover the complete lookup
Gram. These fixtures check its normalization, literal inverse and the
spectator used by the analytic rank witness; they are not a depth theorem.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CQ2, CZ, Q2, adjoint, multiply, native_columns
from tests.test_syndrome_star_readout import _inverse, _syndrome, _wide_t, _x


def _ty(bit):
    # Chronological S-dagger, H, T, H, S gives literal T(Y).
    return [('S', bit)] * 3 + [('H', bit), ('T', bit),
                              ('H', bit), ('S', bit)]


def _p1_broad_acceptance(with_reflections=True):
    # a=0; C=(1,2); D=(3,4); q=5. Between walls B=C+D.
    wall = [('CX', 1, 3), ('CX', 2, 4)]
    f = _wide_t((3, 4))
    cr = ([('CX', 4, 0)] + _x(0)
          + [('H', 5), ('CX', 0, 5), ('H', 5)]
          + _x(0) + [('CX', 4, 0)])
    middle = f + cr + _inverse(f) + _inverse(cr) if with_reflections else []
    cz = [('H', 3), ('CX', 5, 3), ('H', 3)]
    return (wall + _ty(5) + middle + cz + _inverse(_ty(5))
            + [('H', 5)] + list(reversed(wall)))


def _syndrome_commutator(p, address):
    size = 1 << p
    f = native_columns(size, _wide_t(tuple(range(size))))
    signs = [-1 if _syndrome(b, size) == address else 1
             for b in range(1 << size)]
    rf = [[signs[row] * value for value in entries]
          for row, entries in enumerate(f)]
    u = multiply(adjoint(f), rf)
    return [[signs[row] * value for value in entries]
            for row, entries in enumerate(u)]


class BroadSyndromeReadoutTests(unittest.TestCase):
    def test_native_broad_block_normalization_and_actual_inverse(self):
        word = _p1_broad_acceptance()
        actual = native_columns(6, word)
        baseline = native_columns(6, _p1_broad_acceptance(False))
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 14)
        self.assertEqual(sum(g[0] in ('T', 'TDG')
                             for g in _p1_broad_acceptance(False)), 2)
        scale = CQ2(Q2(0, Fraction(1, 4)), Q2(0, 0))
        for address in range(2):
            u = _syndrome_commutator(1, address)
            k = [[scale * ((CO if row == col else CZ)
                          + (-1 if row & 1 else 1) * u[row][col])
                  for col in range(4)] for row in range(4)]
            gram = multiply(adjoint(k), k)
            self.assertEqual(gram, [[
                ((2 * CO if row == col else CZ)
                 + (-1 if row & 1 else 1) * u[row][col]
                 + (-1 if col & 1 else 1) * u[col][row].conj()) / 8
                for col in range(4)] for row in range(4)])
            for c in range(4):
                for d_in in range(4):
                    column = address | (c << 1) | (d_in << 3)
                    for row in range(32):
                        same_fiber = (row & 7) == (column & 7)
                        d_out = row >> 3
                        expected = k[d_out ^ c][d_in ^ c] if same_fiber else CZ
                        self.assertEqual(actual[row][column], expected)
                        expected_baseline = (2 * scale if row == column
                                             and not ((c ^ d_in) & 1) else CZ)
                        self.assertEqual(baseline[row][column], expected_baseline)

    def test_row_zero_spectator_and_rank_witness_dimensions(self):
        # Exact p=2 matrix checks the extra B0 spectator. Its restricted
        # rank bound, rather than the unrestricted physical-space rank,
        # is what permits the analytic failure threshold p=3.
        u = _syndrome_commutator(2, 1)
        for row in range(16):
            for column in range(16):
                if (row ^ column) & 1:
                    self.assertEqual(u[row][column], CZ)
                else:
                    self.assertEqual(u[row][column], u[row ^ 1][column ^ 1])
        size, address = 8, 1
        target_positive = [b for b in range(1 << size)
                           if b & 1 and not ((b >> address) & 1)]
        syndrome_zero_in_spectator = [
            b for b in range(1 << size)
            if b & 1 and _syndrome(b, size) == address]
        self.assertEqual(len(target_positive), 64)
        self.assertEqual(len(syndrome_zero_in_spectator), 16)
        self.assertGreater(len(target_positive),
                           2 * len(syndrome_zero_in_spectator))


if __name__ == '__main__':
    unittest.main()
