"""Exact native witness for the exceptional paired occupied replacement.

The supplied three-band word and its inverse retain all physical columns.
The one-T replacement agrees on the whole initialized linear subspace,
while deliberately giving a different unused-sector extension.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_projected_band_matching import _cq_root
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


class JointNonpauliProbeAlgebraTests(unittest.TestCase):
    def test_exceptional_pair_has_literal_one_t_occupied_replacement(self):
        # v=0 is arbitrary data; t=1 is the only initialized qubit.
        # Each coordinate word sends the two independent physical Z axes
        # to {Z_v P_t, P_t}, so each following factor is one native band.
        coordinates = {
            'Z': [('CX', 1, 0)],
            'X': [('CX', 1, 0), ('H', 1)],
            'Y': [('CX', 1, 0), ('H', 1), ('S', 1)],
        }

        def band(axis, negative_first=False):
            middle = [('TDG' if negative_first else 'T', 0), ('T', 1)]
            change = coordinates[axis]
            return _primitive_inverse(change) + middle + change

        ea_word = band('X', True)   # T_Xt T_(Zv Xt)-dagger.
        b_word = band('Y', True)    # T_Yt T_(Zv Yt)-dagger.
        ec_word = band('Z')         # T_Zt T_(Zv Zt).
        ka_word, kc_word = [('H', 1), ('S', 1)], [('H', 1)]
        c0_word = [('S', 1)] * 3   # S_t-dagger, taking Y_t to X_t.
        ea, eb, ec = [native_columns(2, word)
                      for word in (ea_word, b_word, ec_word)]
        ka, kc, c0 = [native_columns(2, word)
                      for word in (ka_word, kc_word, c0_word)]

        xt = _tau(2, 2)
        zv, zt, yt = _tau(4, 2), _tau(8, 2), _tau(10, 2)
        for axis, labels in (('X', (6, 2)), ('Y', (14, 10)),
                             ('Z', (12, 8))):
            change = native_columns(2, coordinates[axis])
            for initial, expected in zip((zv, zt), labels):
                self.assertEqual(_conjugate(change, initial), _tau(expected, 2))
        self.assertEqual(_conjugate(ka, zt), yt)
        self.assertEqual(_conjugate(kc, zt), xt)

        pv0 = scale(add(eye(4), zv), Fraction(1, 2))
        pv1 = scale(add(eye(4), scale(zv, -1)), Fraction(1, 2))
        j = _conjugate(ea, yt)
        s = _conjugate(ec, xt)
        self.assertEqual(j, add(multiply(pv0, yt), multiply(pv1, zt)))
        self.assertEqual(s, add(multiply(pv0, yt), multiply(pv1, xt)))
        self.assertEqual(_conjugate(eb, j), s)
        # The middle axes really form one swapped clean orbit.
        self.assertEqual(_conjugate(j, yt), _tau(14, 2))
        self.assertEqual(_conjugate(j, _tau(14, 2)), yt)

        f0_word = ea_word + b_word + _primitive_inverse(ec_word)
        f0 = native_columns(2, f0_word)
        self.assertEqual(f0, multiply(adjoint(ec), multiply(eb, ea)))
        self.assertEqual(_conjugate(f0, yt), xt)
        self.assertEqual(_conjugate(c0, yt), xt)
        g0 = multiply(adjoint(c0), f0)
        sy_word = ([('S', 1)] * 3 + [('H', 1), ('S', 1), ('H', 1)]
                   + [('S', 1)])
        sy = native_columns(2, sy_word)  # S_Y=(S H) S (S H)-dagger.
        self.assertEqual(g0, add(pv0, scale(multiply(pv1, sy), _cq_root(-1))))
        self.assertEqual(multiply(g0, yt), multiply(yt, g0))
        self.assertEqual(multiply(g0, zv), multiply(zv, g0))

        # On P_A=Y_t's positive sector, the exact scalar is gamma=1
        # and q(Z_v)=T_v-dagger. No projective phase is discarded.
        pa_plus = scale(add(eye(4), yt), Fraction(1, 2))
        q_word = [('TDG', 0)]
        q = native_columns(2, q_word)
        self.assertEqual(multiply(g0, pa_plus), multiply(q, pa_plus))
        kc_prime = multiply(adjoint(c0), kc)
        self.assertEqual(kc_prime, ka)
        replacement = multiply(adjoint(kc_prime), multiply(q, ka))
        self.assertEqual(replacement, q)

        original_word = (ka_word + f0_word + _primitive_inverse(kc_word))
        original = native_columns(2, original_word)
        self.assertEqual(original, multiply(adjoint(kc), multiply(f0, ka)))
        # All four supplied physical columns have phase v*(2*t-1).
        for column in range(4):
            v, t = column & 1, (column >> 1) & 1
            for row in range(4):
                self.assertEqual(original[row][column],
                                 _cq_root(v * (2 * t - 1)) if row == column else CZ)
                self.assertEqual(replacement[row][column],
                                 _cq_root(-v) if row == column else CZ)

        clean_zero = scale(add(eye(4), zt), Fraction(1, 2))
        self.assertEqual(multiply(original, clean_zero),
                         multiply(replacement, clean_zero))
        self.assertNotEqual(original, replacement)
        self.assertNotEqual(original[3][3], replacement[3][3])
        # Equality of the complete occupied linear maps includes arbitrary
        # v superpositions and external reference entanglement.
        self.assertEqual(native_columns(2, original_word, columns=(0, 1)),
                         native_columns(2, q_word, columns=(0, 1)))
        for word, matrix in ((ea_word, ea), (b_word, eb), (ec_word, ec),
                             (f0_word, f0), (original_word, original),
                             (q_word, replacement)):
            self.assertEqual(native_columns(2, _primitive_inverse(word)), adjoint(matrix))
        self.assertEqual(multiply(adjoint(original), clean_zero),
                         multiply(adjoint(replacement), clean_zero))
        self.assertEqual(native_columns(2, _primitive_inverse(q_word)),
                         native_columns(2, [('T', 0)]))


if __name__ == '__main__':
    unittest.main()
