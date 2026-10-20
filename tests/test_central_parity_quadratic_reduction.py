"""Exact guards for one central-parity restriction of quadratic controls.

These small full-column fixtures retain native phases and actual inverse
words. The dependent-list guard is not a completed Boolean clean source
and does not supply an unrestricted circuit lower bound.
"""
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_phase_signal_query_limits import _ccz, _cz, _x
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _clifford_columns, _conjugate, _tau


def _block(matrix, indices):
    return [[matrix[row][column] for column in indices] for row in indices]


class CentralParityQuadraticReductionTests(unittest.TestCase):
    def test_one_relation_can_hide_nondiagonal_activated_axes(self):
        # E = product_i T_(X_t Z_i), with four genuinely independent axes.
        # Its central parity p=Z1 Z2 Z3 Z4 is their literal product.
        width = 5
        change = [('CX', 0, bit) for bit in range(1, 5)] + [('H', 0)]
        word = _inverse(change) + [('T', bit) for bit in range(1, 5)] + change
        axes = [1 | (1 << (width + bit)) for bit in range(1, 5)]
        self.assertEqual(_rank(axes), 4)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [4])
        e = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(e))
        self.assertEqual(multiply(adjoint(e), e), eye(32))

        # On p=+1, use t,z1,z2,z3 as the SAME coordinates on every column.
        # The fourth bit is their center parity z4=z1+z2+z3.
        indices = [value | ((((value >> 1) & 7).bit_count() & 1) << 4)
                   for value in range(16)]
        sector = _block(e, indices)
        sector_inverse = _block(native_columns(width, _inverse(word)), indices)
        self.assertEqual(sector_inverse, adjoint(sector))
        self.assertEqual(multiply(sector_inverse, sector), eye(16))
        for column in range(32):
            parity = ((column >> 1) & 15).bit_count() & 1
            for row in range(32):
                if (((row >> 1) & 15).bit_count() & 1) != parity:
                    self.assertEqual(e[row][column], CZ)
        for column in range(16):
            z = column >> 1
            destination = column ^ 1 if z in (0, 7) else column
            phase = CO if z == 0 else -CO if z == 7 else CI
            self.assertEqual([sector[row][column] for row in range(16)],
                             [phase if row == destination else CZ
                              for row in range(16)])

        projected_axes = [1 | (1 << (4 + bit)) for bit in range(1, 4)]
        projected_axes.append(1 | (14 << 4))
        self.assertEqual(len(set(projected_axes)), 4)
        self.assertEqual(_rank(projected_axes), 3)
        self.assertEqual(projected_axes[0] ^ projected_axes[1]
                         ^ projected_axes[2] ^ projected_axes[3], 0)
        self.assertTrue(all(axis & 15 == 1 for axis in projected_axes))

        # E_s^dagger Z_t E_s is diagonal, although every activated axis
        # has nonzero X mask. Its quadratic center sign is not a Pauli.
        dressed = _conjugate(adjoint(sector), _tau(1 << 4, 4))
        for column in range(16):
            z1, z2, z3 = ((column >> bit) & 1 for bit in range(1, 4))
            q = 1 ^ z1 ^ z2 ^ z3 ^ (z1 & z2) ^ (z1 & z3) ^ (z2 & z3)
            phase = -CO if q ^ (column & 1) else CO
            self.assertEqual([dressed[row][column] for row in range(16)],
                             [phase if row == column else CZ for row in range(16)])
        self.assertEqual(dressed[0][0] * dressed[2][2]
                         * dressed[4][4] * dressed[6][6], -CO)
        # An affine diagonal Pauli sign has product +1 on this affine square.
        # Thus E_s is non-Clifford: its actual Pauli image is non-Pauli.

    def test_rank_two_signed_insertion_becomes_literal_clifford_on_one_parity(self):
        # Retained rank-two native sandwich, now with the signed H=-Z_t.
        # b=z1(1+z2); the insertion is Clifford after fixing z1 alone.
        change = [('CX', 0, 1), ('CX', 0, 2), ('H', 0)]
        e_word = _inverse(change) + [('T', 1), ('T', 2)] + change
        f_word = _x(2) + _cz(1, 2) + _ccz((0, 1, 2)) + _x(2)
        word = e_word + f_word + _inverse(e_word)
        e, f, actual = (native_columns(3, source)
                        for source in (e_word, f_word, word))
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(8))
        self.assertEqual(actual, f)
        for column in range(8):
            t, z1, z2 = ((column >> bit) & 1 for bit in range(3))
            phase = -CO if (1 ^ t) & z1 & (1 ^ z2) else CO
            self.assertEqual([actual[row][column] for row in range(8)],
                             [phase if row == column else CZ for row in range(8)])

        omega = native_columns(1, [('T', 0)])[1][1]
        sector_change = [('CX', 0, 1), ('H', 0)]
        for s in (0, 1):
            with self.subTest(center_bit=s):
                indices = [(value & 1) | (s << 1) | ((value & 2) << 1)
                           for value in range(4)]
                es, fs, vs = (_block(matrix, indices) for matrix in (e, f, actual))
                es_word = (_inverse(sector_change)
                           + [('TDG' if s else 'T', 0), ('T', 1)]
                           + sector_change)
                es_native = native_columns(2, es_word)
                # T_(-X)=omega T_X^dagger: retain the projected scalar.
                self.assertEqual(es, scale(es_native, omega if s else CO))
                fs_word = (_x(1) + [('S', 1), ('S', 1)]
                           + _cz(1, 0) + _x(1)) if s else []
                self.assertEqual(fs, native_columns(2, fs_word))
                self.assertEqual(native_columns(2, _inverse(fs_word)), adjoint(fs))
                _clifford_columns(fs, 2)
                self.assertEqual(multiply(fs, fs), eye(4))

                # E_C'=F_s^dagger E_C F_s, K_C'=F_s^dagger. These are
                # actual Clifford conjugations; the two T targets stay native.
                ec_prime_word = _inverse(fs_word) + es_word + fs_word
                ec_prime = native_columns(2, ec_prime_word)
                self.assertEqual(ec_prime, _conjugate(adjoint(fs), es_native))
                self.assertEqual([len(band) for band in _native_t_bands(ec_prime_word)],
                                 [2])
                repacked_word = es_word + _inverse(ec_prime_word) + fs_word
                self.assertEqual(native_columns(2, repacked_word), vs)
                self.assertEqual(native_columns(2, _inverse(repacked_word)), adjoint(vs))
                self.assertEqual(multiply(adjoint(vs), vs), eye(4))


if __name__ == '__main__':
    unittest.main()
