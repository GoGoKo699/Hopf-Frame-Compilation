"""Exact source guards for the occupied Boolean common-center reduction.

The positive source and both guards are actual native one-clean words.
The quadratic-control word tests a proposed inference about finite core
support; it is not claimed to be a fixed-pair three-band Boolean source.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_original_middle_pair_insertions import _ghz, _pair_word
from tests.test_t_depth import _word
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _cz(first, second):
    return [('H', second), ('CX', first, second), ('H', second)]


def _z(wire, width):
    return _tau(1 << (width + wire), width)


def _x(wire, width):
    return _tau(1 << wire, width)


class OccupiedBooleanPairCenterTests(unittest.TestCase):
    def assert_actual_inverse_and_clean_return(self, width, word):
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(1 << width))
        self.assertEqual(_conjugate(actual, _z(0, width)), _z(0, width))
        return actual

    def test_boolean_source_with_nonpauli_images_and_affine_derivative(self):
        width = 3
        # K_A = GHZ_01 H_1, K_C = K_A CZ_12, E_A=E_C=T_0 T_1.
        # The whole physical word is CZ_12, with wire 0 returned.
        ka_word = [('H', 1), ('H', 0), ('CX', 0, 1)]
        phase_word = _cz(1, 2)
        kc_word = phase_word + ka_word
        band = [('T', 0), ('T', 1)]
        ka, kc, e = [native_columns(width, word)
                     for word in (ka_word, kc_word, band)]
        pa = _conjugate(ka, _z(0, width))
        pc = _conjugate(kc, _z(0, width))
        self.assertEqual(pa, pc)
        self.assertEqual(pa, _tau(3, width))
        j, s = _conjugate(e, pa), _conjugate(e, pc)
        self.assertEqual(j, s)
        # A nontrivial binary action rules out a Pauli clean image.
        self.assertNotIn(_conjugate(j, _x(0, width)),
                         (_x(0, width), scale(_x(0, width), -1)))
        for axis in (_z(0, width), _z(1, width)):
            self.assertEqual(_conjugate(axis, pa), scale(pa, -1))
        q = _tau(3 << width, width)  # R=<Z_0 Z_1>.
        self.assertEqual(_conjugate(q, pa), pa)
        self.assertEqual(_conjugate(q, e), e)
        word = ka_word + band + _inverse(band) + _inverse(kc_word)
        actual = self.assert_actual_inverse_and_clean_return(width, word)
        self.assertEqual(actual, native_columns(width, phase_word))
        qa, qc = _conjugate(adjoint(ka), q), _conjugate(adjoint(kc), q)
        self.assertEqual(qa, _x(1, width))
        self.assertEqual(qc, _tau(2 | (4 << width), width))
        self.assertEqual(_conjugate(actual, qa), qc)
        # The transported Boolean derivative is the actual affine sign Z_2.
        self.assertEqual(multiply(_conjugate(actual, qa), qa), _z(2, width))
        for column in range(1 << width):
            sign = -1 if column & 2 and column & 4 else 1
            self.assertEqual(actual[column][column], CO * sign)
            self.assertTrue(all(actual[row][column] == CZ
                                for row in range(1 << width) if row != column))

    def test_pauli_clean_guard_requires_nonempty_opposite_active_set(self):
        width = 2
        kc_word = [('CX', 1, 0)]
        ka_word = kc_word + [('S', 0)]
        ec_word = [('T', 0), ('T', 1)]
        ka, kc, ec = [native_columns(width, word)
                      for word in (ka_word, kc_word, ec_word)]
        pa = _conjugate(ka, _z(0, width))
        pc = _conjugate(kc, _z(0, width))
        q = _tau(3, width)
        self.assertEqual(pa, _tau(3 << width, width))
        self.assertEqual(pa, pc)
        self.assertEqual(_conjugate(ec, pc), pc)
        self.assertEqual(_conjugate(q, pa), pa)
        for axis in (_z(0, width), _z(1, width)):
            self.assertEqual(_conjugate(q, axis), scale(axis, -1))
        word = ka_word + _inverse(ec_word) + _inverse(kc_word)
        actual = self.assert_actual_inverse_and_clean_return(width, word)
        # Every occupied input is fixed literally, including the dirty bit.
        for column in (0, 2):
            self.assertEqual([actual[row][column] for row in range(4)],
                             [CO if row == column else CZ for row in range(4)])

    def test_nonboolean_guard_even_with_both_active_sets_nonempty(self):
        width = 2
        k_word = [('H', 0)]
        ea_word, ec_word = [('T', 0)], [('T', 0), ('T', 1)]
        k, ea, ec = [native_columns(width, word)
                    for word in (k_word, ea_word, ec_word)]
        p, q = _x(0, width), _x(1, width)
        self.assertEqual(_conjugate(k, _z(0, width)), p)
        self.assertEqual(_conjugate(q, p), p)
        self.assertEqual(_conjugate(q, ea), ea)
        self.assertEqual(_conjugate(ea, p), _conjugate(ec, p))
        j = _conjugate(ea, p)
        self.assertNotIn(_conjugate(j, p), (p, scale(p, -1)))
        self.assertEqual(_conjugate(q, _z(1, width)), scale(_z(1, width), -1))
        word = k_word + ea_word + _inverse(ec_word) + _inverse(k_word)
        actual = self.assert_actual_inverse_and_clean_return(width, word)
        self.assertEqual(actual, native_columns(width, [('TDG', 1)]))
        self.assertNotIn(actual[2][2] / actual[0][0], (CO, -CO))

    def test_one_quadratic_feature_has_growing_clifford_norm_distortion(self):
        width = 5
        # F=X_4^(z_0 z_1 + z_2 z_3). Borrowed controls are restored by
        # each complete Toffoli. No initialized bit or abstract gate is used.
        schedule = (_compact_toffoli(0, 1, 4, 2)
                    + _compact_toffoli(2, 3, 4, 0))
        word = _word(schedule)
        actual = native_columns(width, word)
        expected = [[CZ for _ in range(32)] for _ in range(32)]
        for column in range(32):
            q = int(column & 3 == 3) ^ int(column & 12 == 12)
            expected[column ^ (q << 4)][column] = CO
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(actual, actual), eye(32))
        dressed = _conjugate(actual, _z(4, width))
        clifford_word = [('S', 4), ('S', 4)] + _cz(0, 1) + _cz(2, 3)
        self.assertEqual(dressed, native_columns(width, clifford_word))
        # Its binary displacement has four independent central directions.
        differences = []
        for bit in range(4):
            partner = bit ^ 1
            image_label = (1 << bit) | (1 << (width + partner))
            self.assertEqual(_conjugate(dressed, _x(bit, width)),
                             _tau(image_label, width))
            differences.append(1 << (width + partner))
        self.assertEqual(_rank(differences), 4)
        # The exact CZ expansions give 16 distinct coefficients of |1/4|.
        expansion = [[CZ for _ in range(32)] for _ in range(32)]
        for mask in range(16):
            sign = (-1) ** (int(mask & 3 == 3) + int(mask & 12 == 12))
            label = (mask | 16) << width
            expansion = add(expansion, scale(_tau(label, width),
                                             Fraction(sign, 4)))
        self.assertEqual(dressed, expansion)
        self.assertEqual(16 * Fraction(1, 4), 4)
        for candidate in (word, _inverse(word)):
            count = Counter(gate[0] for gate in candidate)
            self.assertEqual(count['T'] + count['TDG'], 16)
            self.assertEqual(count['CX'], 32)
            self.assertEqual(count['H'], 4)

    def test_actual_silent_axes_cancel_the_pair_sector_odd_phase(self):
        width, active = 5, 4
        k_word = _ghz(active)
        ea_word = [('T', bit) for bit in range(width)]
        # U sends Z_t to R Z_t and fixes the four active Z axes.
        parity = [('CX', bit, 4) for bit in range(active)]
        ec_word = parity + ea_word + _inverse(parity)
        pair_word = _pair_word(active, 1, 1)
        k, ea, ec, eb = [native_columns(width, word)
                         for word in (k_word, ea_word, ec_word, pair_word)]
        pa = _conjugate(k, _z(0, width))
        j = _conjugate(ea, pa)
        b = _tau(15, width)
        d = _tau(15 | (15 << width), width)
        self.assertEqual(_conjugate(j, b), d)
        self.assertEqual(_conjugate(j, d), b)
        self.assertEqual(_conjugate(eb, j), _conjugate(ec, pa))
        self.assertEqual(_rank([1 << (width + bit) for bit in range(active)]
                               + [31 << width]), width)
        word = k_word + ea_word + pair_word + _inverse(ec_word) + _inverse(k_word)
        actual = self.assert_actual_inverse_and_clean_return(width, word)
        i = native_columns(1, [('S', 0)])[1][1]
        for column in range(0, 32, 2):
            weight = ((column >> 1) & 7).bit_count()
            if weight % 2:
                phase = i if column & 16 else CO
            else:
                phase = i if weight == 0 else CO
            self.assertEqual([actual[row][column] for row in range(32)],
                             [phase if row == column else CZ for row in range(32)])
        # The omega-sector phase is removed, but the whole occupied word
        # still fails the projective Boolean condition.
        self.assertEqual(actual[2][2], CO)
        self.assertEqual(actual[18][18], i)
        self.assertNotIn(actual[0][0] / actual[2][2], (CO, -CO))
        for candidate in (word, _inverse(word)):
            count = Counter(gate[0] for gate in candidate)
            self.assertEqual(count['T'] + count['TDG'], 12)
            self.assertEqual(count['H'] + count['CX'], 30)


if __name__ == '__main__':
    unittest.main()
