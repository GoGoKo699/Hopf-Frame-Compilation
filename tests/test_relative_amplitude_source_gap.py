"""Exact bounded guards for the amplitude-independent source defect.

These complete native matrices check the projection, literal ordered
expansion and nonzero-compression boundary. The all-size relative gap is
analytic; no fixture is treated as an indicator source or a depth bound.
"""
from itertools import product
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_three_layer_native_mechanisms import _native_t_bands


def _conjugate(unitary, matrix):
    return multiply(multiply(unitary, matrix), adjoint(unitary))


def _product(*factors):
    result = eye(len(factors[0]))
    for factor in factors:
        result = multiply(result, factor)
    return result


def _clean_block(matrix):
    return [[matrix[row][column] for column in (0, 2)] for row in (0, 2)]


def _normalized_hs_squared(matrix):
    return sum(value.abs2() for row in matrix for value in row) / len(matrix)


class RelativeAmplitudeSourceGapTests(unittest.TestCase):
    def test_nonflat_relative_defect_retains_order_phase_and_all_columns(self):
        # Signal0 is the only clean input; logical1 is arbitrary. Each
        # encoder contains one actual T band, with overlapping logical axes.
        ua_word = [('CX', 1, 0), ('H', 1), ('T', 1)]
        ub_word = [('CX', 1, 0), ('H', 0), ('H', 1), ('T', 1), ('H', 1)]
        source = ua_word + _inverse(ub_word)
        ua, ub = native_columns(2, ua_word), native_columns(2, ub_word)
        actual = native_columns(2, source)
        self.assertEqual(actual, multiply(adjoint(ub), ua))
        self.assertEqual([len(band) for band in _native_t_bands(source)], [1, 1])
        for word, matrix in ((ua_word, ua), (ub_word, ub), (source, actual)):
            self.assertEqual(native_columns(2, _inverse(word)), adjoint(matrix))
            self.assertEqual(multiply(adjoint(matrix), matrix), eye(4))

        j = [[CO if row == 2 * column else CZ for column in range(2)]
             for row in range(4)]
        pi = multiply(j, adjoint(j))
        pi_a, pi_b = _conjugate(ua, pi), _conjugate(ub, pi)
        k = _clean_block(actual)
        gram = multiply(adjoint(k), k)
        p = native_columns(1, [('Z', 0)])
        physical_p = native_columns(2, [('Z', 1)])
        c, d = _conjugate(ub, physical_p), _conjugate(ua, physical_p)
        self.assertNotEqual(gram, scale(eye(2), (gram[0][0] + gram[1][1]) / 2))
        self.assertNotEqual(multiply(p, gram), multiply(gram, p))

        direct = add(_product(pi_a, pi_b, c, pi_a),
                     scale(_product(d, pi_a, pi_b, pi_a), -CO))
        logical = add(_product(adjoint(k), p, k), scale(multiply(p, gram), -CO))
        embedded = multiply(multiply(j, logical), adjoint(j))
        self.assertEqual(direct, _conjugate(ua, embedded))
        self.assertEqual(_product(pi_a, direct, pi_a), direct)
        self.assertNotEqual(direct, zeros(4))
        self.assertEqual(_normalized_hs_squared(direct),
                         _normalized_hs_squared(logical) / 2)

        # Keep the two three-projector expansions in their actual order.
        # Their coefficients are fixed +/-1/8, with no chosen amplitude.
        za = native_columns(2, [('Z', 0)])
        images_a = [eye(4), _conjugate(ua, za)]
        images_b = [eye(4), _conjugate(ub, za)]
        expansion, terms = zeros(4), 0
        for left, middle, right in product(range(2), repeat=3):
            positive = _product(images_a[left], images_b[middle], c, images_a[right])
            negative = _product(d, images_a[left], images_b[middle], images_a[right])
            expansion = add(expansion, scale(add(positive, scale(negative, -CO)), CO / 8))
            terms += 2
        self.assertEqual(terms, 16)
        self.assertEqual(expansion, direct)

        # A common literal source phase is emitted by a Clifford word.
        # Its actual inverse is retained, and both quadratic terms cancel it.
        phase_word = [('H', 1), ('S', 1)] * 3
        phased_word = source + phase_word
        phased = native_columns(2, phased_word)
        omega = native_columns(1, [('T', 0)])[1][1]
        self.assertEqual(phased, scale(actual, omega))
        self.assertEqual(native_columns(2, _inverse(phased_word)), adjoint(phased))
        phased_k = _clean_block(phased)
        phased_gram = multiply(adjoint(phased_k), phased_k)
        self.assertEqual(add(_product(adjoint(phased_k), p, phased_k),
                             scale(multiply(p, phased_gram), -CO)), logical)

    def test_one_signal_nonzero_norm_boundary_and_zero_block(self):
        # U_A=I and U_B=T_Y X on signal0 give A=Z0,B=-H0. The
        # arbitrary logical wire1 is a spectator. T_Y is a genuine
        # one-band native word in the D=SH Pauli coordinates.
        x_word = [('H', 0), ('S', 0), ('S', 0), ('H', 0)]
        coordinates = [('H', 0), ('S', 0)]
        ub_word = x_word + _inverse(coordinates) + [('T', 0)] + coordinates
        ub = native_columns(2, ub_word)
        self.assertEqual([len(band) for band in _native_t_bands(ub_word)], [1])
        a = native_columns(2, [('Z', 0)])
        b = _conjugate(ub, a)
        self.assertEqual(b, scale(native_columns(2, [('H', 0)]), -CO))
        actual = native_columns(2, _inverse(ub_word))
        self.assertEqual(actual, adjoint(ub))
        k = _clean_block(actual)

        half_root = native_columns(1, [('H', 0)])[0][0]
        sine_squared = (CO - half_root) / 2  # sin(pi/8)^2, exactly.
        self.assertEqual(multiply(adjoint(k), k), scale(eye(2), sine_squared))
        summed = add(a, b)
        self.assertEqual(multiply(adjoint(summed), summed),
                         scale(eye(4), 4 * sine_squared))
        # These scalar Gram identities certify ||A+B||=2||K|| and the
        # attained nonzero boundary without numerical singular values.

        orthogonal = native_columns(2, x_word)
        zero_b = _conjugate(orthogonal, a)
        self.assertEqual(zero_b, scale(a, -CO))
        self.assertEqual(add(a, zero_b), zeros(4))
        self.assertEqual(_clean_block(adjoint(orthogonal)), zeros(2))
        self.assertEqual(native_columns(2, _inverse(x_word)), adjoint(orthogonal))


if __name__ == '__main__':
    unittest.main()
