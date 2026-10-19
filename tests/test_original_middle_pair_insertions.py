"""Literal physical-pair insertion checks, not a Boolean-degree theorem.

The two-wire cases retain every signed swap and native orientation. The
four-wire source checks genuine support growth through the endpoint band;
it returns one clean wire but is not claimed as a Boolean query.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CZ, CO, CI, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _ghz(width):
    return [('H', 0)] + [('CX', 0, j) for j in range(1, width)]


def _pair_coordinates(width):
    # U Z0 U-dagger = X^width and U Z1 U-dagger = Y^width.
    sign = [('H', 1), ('S', 1), ('S', 1), ('H', 1)] if width % 4 == 2 else []
    return sign + [('CX', j, 1) for j in range(width) if j != 1] + _ghz(width)


def _pair_word(width, first, second):
    coordinates = _pair_coordinates(width)
    band = [('T' if exponent == 1 else 'TDG', wire)
            for wire, exponent in enumerate((first, second))]
    return _inverse(coordinates) + band + coordinates


def _correction_word(width, exponent):
    if exponent == 0:
        return []
    coordinates = _pair_coordinates(width)
    return (_inverse(coordinates)
            + [('S' if exponent == 2 else 'SDG', 1)] + coordinates)


class OriginalMiddlePairInsertionTests(unittest.TestCase):
    def test_signed_native_pairs_and_actual_clean_normal_form(self):
        width = 2
        size = 1 << width
        e_word = [('T', 0), ('T', 1)]
        e = native_columns(width, e_word)
        p, q = _tau(3, width), _tau(15, width)
        u = native_columns(width, _pair_coordinates(width))
        self.assertEqual(_conjugate(u, _tau(1 << width, width)), p)
        self.assertEqual(_conjugate(u, _tau(2 << width, width)), q)
        for swap_sign in (1, -1):
            # The second endpoint Pauli is XY when swap_sign is negative.
            k_word = _ghz(width) + ([('S', 1)] if swap_sign == -1 else [])
            k = native_columns(width, k_word)
            clean = _conjugate(k, _tau(1 << width, width))
            j = _conjugate(e, clean)
            self.assertEqual(_conjugate(j, p), scale(q, swap_sign))
            self.assertEqual(_conjugate(j, q), scale(p, swap_sign))
            for epsilon in (1, -1):
                for delta in (1, -1):
                    with self.subTest(s=swap_sign, epsilon=epsilon, delta=delta):
                        b_word = _pair_word(width, epsilon, delta)
                        f_word = _pair_word(width, epsilon, swap_sign * epsilon)
                        c_word = _correction_word(width, delta - swap_sign * epsilon)
                        b, f, c = [native_columns(width, word)
                                   for word in (b_word, f_word, c_word)]
                        self.assertEqual(b, multiply(c, f))
                        self.assertEqual(_conjugate(f, j), j)
                        # E_C = C E_A C-dagger, K_C = C K_A.
                        ec = _conjugate(c, e)
                        kc = multiply(c, k)
                        pc = _conjugate(kc, _tau(1 << width, width))
                        self.assertEqual(_conjugate(b, j), _conjugate(ec, pc))
                        self.assertEqual(_conjugate(adjoint(c), ec), e)
                        self.assertEqual(multiply(adjoint(c), kc), k)
                        # Chronological original word, followed by actual reverse.
                        original = (k_word + e_word + b_word + _inverse(c_word)
                                    + _inverse(e_word) + _inverse(k_word))
                        normalized = k_word + e_word + f_word + _inverse(e_word) + _inverse(k_word)
                        actual = native_columns(width, original)
                        self.assertEqual(actual, native_columns(width, normalized))
                        self.assertEqual(native_columns(width, _inverse(original)), adjoint(actual))
                        self.assertEqual(native_columns(width, _inverse(normalized)), adjoint(actual))
                        # Full physical commutation proves clean return on all data.
                        self.assertEqual(_conjugate(actual, _tau(1 << width, width)),
                                         _tau(1 << width, width))
                        self.assertEqual(multiply(adjoint(actual), actual), eye(size))
                        for word in (original, normalized):
                            self.assertEqual([len(band) for band in _native_t_bands(word)], [2, 2, 2])
                            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 6)

    def test_one_actual_pair_gains_large_endpoint_support(self):
        width = 4
        size = 1 << width
        all_bits = size - 1
        e_word = [('T', bit) for bit in range(width)]
        k_word = _ghz(width)
        f_word = _pair_word(width, 1, 1)
        e, f = [native_columns(width, word) for word in (e_word, f_word)]
        p = _tau(all_bits, width)
        q = _tau(all_bits | (all_bits << width), width)
        j = _conjugate(e, p)
        self.assertEqual(_conjugate(j, p), q)
        self.assertEqual(_conjugate(j, q), p)
        self.assertEqual(_conjugate(f, j), j)
        omega = native_columns(1, [('T', 0)])[1][1]
        alpha, beta = (CO + omega) / 2, (CO - omega) / 2
        expected = add(scale(eye(size), alpha * alpha),
                       scale(_tau(all_bits << width, width), beta * beta))
        support = [0, all_bits << width]
        # Each _tau is the literal tensor product X_(S-complement) Y_S.
        for subset in range(size):
            if subset.bit_count() % 2 == 0:
                label = all_bits | (subset << width)
                support.append(label)
                expected = add(expected, scale(_tau(label, width), alpha * beta * Fraction(1, 2)))
        dressed = _conjugate(adjoint(e), f)
        self.assertEqual(dressed, expected)
        self.assertNotEqual(alpha * beta, CZ)
        self.assertEqual(len(set(support)), 10)
        self.assertEqual(_rank(support), width)
        self.assertEqual(_rank([all_bits, all_bits | (all_bits << width)]), 2)
        word = k_word + e_word + f_word + _inverse(e_word) + _inverse(k_word)
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(_conjugate(actual, _tau(1 << width, width)),
                         _tau(1 << width, width))
        # This is a component witness, not a growing-degree Boolean family.
        self.assertEqual(actual[0][0], CI)
        self.assertEqual(actual[2][2], omega)
        self.assertEqual(actual[6][6], CO)
        self.assertNotIn(actual[2][2] / actual[6][6], (CO, -CO))
        for candidate in (word, _inverse(word)):
            count = Counter(gate[0] for gate in candidate)
            self.assertEqual(count['T'] + count['TDG'], 2 * width + 2)
            self.assertEqual(count['H'] + count['CX'], 6 * width - 2)
            self.assertEqual([len(band) for band in _native_t_bands(candidate)], [width, 2, width])


if __name__ == '__main__':
    unittest.main()
