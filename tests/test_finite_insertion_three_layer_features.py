"""Exact scope fixtures for finite insertions in three-layer features.

These two identities test the need for norm distortion and an intrinsic
Pauli kernel. They are not high-degree queries or general depth bounds.
Every matrix uses literal native phases on the full physical input space.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CI, CO, CZ, Q2, adjoint, eye, native_columns
from tests.test_bounded_relation_native_matching import _rank
from tests.test_two_layer_galois_return import (
    _clifford_columns, _coefficients, _conjugate, _inverse, _product,
    _symplectic, _tau,
)


class FiniteInsertionThreeLayerFeaturesTests(unittest.TestCase):
    def test_commuting_insertion_can_break_raw_middle_rank_matching(self):
        # T_Y=(S H) T_Z (S H)-dagger, including its native scalar.
        outer_word = ([('S', 0)] * 3
                      + [('H', 0), ('T', 0), ('H', 0), ('S', 0)])
        middle_word, insertion_word = [('T', 0)], [('TDG', 0)]
        outer = native_columns(1, outer_word)
        middle = native_columns(1, middle_word)
        insertion = native_columns(1, insertion_word)
        z, x = _tau(2, 1), _tau(1, 1)
        self.assertEqual(_product(insertion, z), _product(z, insertion))
        self.assertEqual(_product(insertion, middle), eye(2))

        transported = _conjugate(outer, z)
        self.assertEqual(transported, native_columns(1, [('H', 0)]))
        self.assertEqual(_conjugate(transported, z), x)
        ranks, support_sizes = [], []
        for probe in (z, x):
            image = _conjugate(middle, probe)
            columns = _clifford_columns(image, 1)
            ranks.append(_rank([column ^ (1 << bit)
                                for bit, column in enumerate(columns)]))
            support_sizes.append(sum(value != CZ for value in _coefficients(image, 1)))
        self.assertEqual(ranks, [0, 1])
        self.assertEqual(support_sizes, [1, 2])
        # V=E_C-dagger F E_B E_A=I with E_C=E_A, despite that
        # failure of rank-profile invariance for the unmodified E_B.
        word = outer_word + middle_word + insertion_word + _inverse(outer_word)
        actual = native_columns(1, word)
        self.assertEqual(actual, eye(2))
        self.assertEqual(native_columns(1, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(1, _inverse(outer_word)), adjoint(outer))

    def test_non_c3_insertion_has_nonisotropic_support_and_zero_core_kernel(self):
        # Wire v=0 is central; wire t=1 is the entire finite core.
        # Rightmost operator factors act first: F=T_Xt T_Zt.
        width = 2
        a_word, b_word = [('TDG', 1)], [('T', 0)]
        f_word = [('T', 1), ('H', 1), ('T', 1), ('H', 1)]
        c_band = [('T', 0), ('T', 1)]
        c_word = [('H', 1)] + c_band + [('H', 1)]
        self.assertEqual({gate[1] for gate in c_band}, {0, 1})
        insertion = native_columns(width, f_word)
        x_core, z_core, z_central = 2, 8, 4
        central_axis = _tau(z_central, width)
        self.assertEqual(_product(insertion, central_axis),
                         _product(central_axis, insertion))

        omega = native_columns(1, [('T', 0)])[1][1]
        alpha, beta = (CO + omega) / 2, (CO - omega) / 2
        expected = [CZ] * 16
        expected[0] = alpha * alpha
        expected[x_core] = expected[z_core] = alpha * beta
        expected[x_core ^ z_core] = -CI * beta * beta  # XZ=-iY.
        coefficients = _coefficients(insertion, width)
        self.assertEqual(coefficients, expected)
        self.assertTrue(any(value.imag != 0 for value in coefficients))
        support = [label for label, value in enumerate(coefficients) if value != CZ]
        self.assertEqual(support, [0, x_core, z_core, x_core ^ z_core])
        self.assertEqual(_rank(support), 2)
        self.assertEqual(_symplectic(x_core, z_core, width), 1)

        # Enumerate the exact Pauli-preserving domain on all physical
        # wires: only the spectator central Pauli algebra survives.
        kernel = []
        for label in range(16):
            image = _conjugate(insertion, _tau(label, width))
            image_coefficients = _coefficients(image, width)
            if sum(value != CZ for value in image_coefficients) == 1:
                kernel.append(label)
        self.assertEqual(kernel, [0, 1, 4, 5])
        self.assertEqual(set(kernel) & {0, x_core, z_core, x_core ^ z_core}, {0})
        non_clifford_image = _conjugate(insertion, _tau(x_core, width))
        weights = {value.abs2() for value in _coefficients(non_clifford_image, width)
                   if value != CZ}
        self.assertEqual(weights, {Q2(Fraction(1, 2)), Q2(Fraction(1, 4))})
        with self.assertRaises(ValueError):
            _clifford_columns(non_clifford_image, width)

        # E_C=T_Xt T_Zv is one independent-axis band; the complete
        # returned word is identity for every v/t input, with no projection.
        word = a_word + b_word + f_word + _inverse(c_word)
        actual = native_columns(width, word)
        self.assertEqual(actual, eye(4))
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(native_columns(width, _inverse(f_word)), adjoint(insertion))


if __name__ == '__main__':
    unittest.main()
