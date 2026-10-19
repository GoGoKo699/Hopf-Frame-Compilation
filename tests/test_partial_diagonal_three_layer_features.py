"""Exact fixtures for active-sector lifts and pointwise phase alphabets.

An inactive extension need not be Hermitian. Alphabet types need not
commute when they occur on different central fibers. Literal native
phases and inverse words distinguish the valid statements from both
stronger, false shortcuts.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_projected_band_matching import _cq_root
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _power(matrix, count):
    result = eye(len(matrix))
    for _ in range(count):
        result = multiply(matrix, result)
    return result


class PartialDiagonalThreeLayerFeaturesTests(unittest.TestCase):
    def test_active_sector_lift_needs_no_hermitian_inactive_extension(self):
        # Core t=0, central z=1. E=T_(Y_z) T_(Z_t Y_z)†.
        coordinates = [('CX', 1, 0), ('H', 1), ('S', 1)]
        e_word = _primitive_inverse(coordinates) + [('TDG', 0), ('T', 1)] + coordinates
        e = native_columns(2, e_word)
        z_t, z_z = _tau(4, 2), _tau(8, 2)
        x_z, y_z = _tau(2, 2), _tau(10, 2)
        plus = scale(add(eye(4), z_t), Fraction(1, 2))
        minus = scale(add(eye(4), scale(z_t, -1)), Fraction(1, 2))
        h = _conjugate(e, z_z)
        self.assertEqual(h, add(multiply(plus, z_z), multiply(minus, x_z)))
        self.assertNotEqual(multiply(h, z_z), multiply(z_z, h))
        # P=T_(Z_z) T_(Z_t Z_z) is one actual independent-axis band.
        pair_word = [('CX', 1, 0), ('T', 0), ('T', 1), ('CX', 1, 0)]
        f_word = _primitive_inverse(e_word) + pair_word + e_word
        f = native_columns(2, f_word)
        s_h = scale(add(scale(eye(4), CO + CI), scale(h, CO - CI)), Fraction(1, 2))
        expected_f = add(multiply(plus, s_h), scale(minus, _cq_root(1)))
        self.assertEqual(f, expected_f)
        self.assertEqual(multiply(f, z_z), multiply(z_z, f))
        # The actual controlled-Y correction U is identity on t=0.
        u_word = [('S', 1)] * 3 + [('CX', 0, 1), ('S', 1)]
        u = native_columns(2, u_word)
        self.assertEqual(u, add(plus, multiply(minus, y_z)))
        lifted = multiply(u, h)
        self.assertEqual(lifted, add(multiply(plus, z_z), scale(multiply(minus, z_z), -CI)))
        # Independently emitted diagonal Clifford: S_t† Z_z.
        lifted_word = [('S', 0)] * 3 + [('S', 1)] * 2
        self.assertEqual(lifted, native_columns(2, lifted_word))
        self.assertNotEqual(lifted, adjoint(lifted))
        self.assertEqual(multiply(lifted, z_z), multiply(z_z, lifted))
        self.assertEqual(multiply(plus, lifted), multiply(plus, h))
        # Use the linear S formula only after projecting to the active
        # sector. It is not asserted to be a unitary on the inactive one.
        linear_s_lift = scale(add(scale(eye(4), CO + CI), scale(lifted, CO - CI)), Fraction(1, 2))
        replacement = add(multiply(plus, linear_s_lift), scale(minus, _cq_root(1)))
        self.assertEqual(replacement, f)
        for word, actual in ((e_word, e), (f_word, f), (u_word, u), (lifted_word, lifted)):
            self.assertEqual(native_columns(2, _primitive_inverse(word)), adjoint(actual))
        self.assertEqual(adjoint(replacement), native_columns(2, _primitive_inverse(f_word)))
        # All physical columns include the inactive block with its omega.
        for column in range(4):
            exponent = 1 if column & 1 else 2 * ((column >> 1) & 1)
            for row in range(4):
                self.assertEqual(f[row][column], _cq_root(exponent) if row == column else CZ)

    def test_pointwise_commuting_pair_alphabet_retains_counts_modulo_eight(self):
        # F_X uses R=Z_z; F_Z uses R=-Z_z. Their non-scalar actions
        # occupy disjoint central fibers, although S_X and S_Z do not commute.
        coordinates_x = [('CX', 0, 1), ('H', 0)]
        coordinates_z = [('CX', 0, 1), ('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        band = [('T', 0), ('T', 1)]
        fx_word = _primitive_inverse(coordinates_x) + band + coordinates_x
        fz_word = _primitive_inverse(coordinates_z) + band + coordinates_z
        fx, fz = native_columns(2, fx_word), native_columns(2, fz_word)
        self.assertEqual(multiply(fx, fz), multiply(fz, fx))
        z_z = _tau(8, 2)
        plus = scale(add(eye(4), z_z), Fraction(1, 2))
        minus = scale(add(eye(4), scale(z_z, -1)), Fraction(1, 2))
        sx = native_columns(2, [('H', 0), ('S', 0), ('H', 0)])
        sz = native_columns(2, [('S', 0)])
        self.assertNotEqual(multiply(sx, sz), multiply(sz, sx))
        self.assertEqual(fx, add(multiply(plus, sx), scale(minus, _cq_root(1))))
        self.assertEqual(fz, add(scale(plus, _cq_root(1)), multiply(minus, sz)))
        # Modulo four would discard a real, literal global minus sign.
        self.assertEqual(_power(fx, 4), z_z)
        self.assertEqual(_power(fz, 4), scale(z_z, -1))
        four_each = (fx_word + fz_word) * 4
        self.assertEqual(native_columns(2, four_each), scale(eye(4), -1))
        self.assertEqual(_power(fx, 8), eye(4))
        self.assertEqual(_power(fz, 8), eye(4))
        # An actual interleaved native word has nine F_X and eight F_Z
        # factors. The pointwise finite-alphabet counts reduce to F_X.
        repeated_word = (fx_word + fz_word) * 8 + fx_word
        repeated = native_columns(2, repeated_word)
        self.assertEqual(repeated, fx)
        for word, actual in ((fx_word, fx), (fz_word, fz),
                             (four_each, scale(eye(4), -1)),
                             (repeated_word, repeated)):
            self.assertEqual(native_columns(2, _primitive_inverse(word)), adjoint(actual))


if __name__ == '__main__':
    unittest.main()
