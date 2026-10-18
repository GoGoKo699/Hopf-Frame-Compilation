"""Bounded evidence for the robust shared-prefix comparison invariant.

The native pullback, projective quadratic law and augmentation degrees
are checked separately. The abstract feature example is not a synthesis
claim for a returned query.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_one_band_axis_matching import (
    _one_qubit_cliffords, _projective_distance,
)
from tests.test_read_only_five_layer_reduction import _p
from tests.test_read_only_four_layer_rigidity import _four_layer_fixture
from tests.test_read_only_three_layer_rigidity import H, S, T, X, _cx, _degree, _wire
from tests.test_span_preserving_augmentation import (
    _dihedral_multiply, _ideal_powers, _reduce,
)


def _inverse(element, multiply, size):
    return next(candidate for candidate in range(size)
                if multiply(element, candidate) == multiply(candidate, element) == 0)


def _product(elements, multiply):
    result = 0
    for element in elements:
        result = multiply(result, element)
    return result


def _inverse_word(word):
    return [(address, -sign) for address, sign in reversed(word)]


def _defect_word(x, y):
    return [(x ^ y, 1), (y, -1), (x, -1)]


def _relation_words(x, y, z):
    defect = _defect_word(x, y)
    return [
        defect + [(z, 1)] + _inverse_word(defect) + [(z, -1)],
        defect + defect,
        (_defect_word(x ^ y, z) + _inverse_word(_defect_word(y, z))
         + _inverse_word(_defect_word(x, z))),
        (_defect_word(x, y ^ z) + _inverse_word(_defect_word(x, z))
         + _inverse_word(_defect_word(x, y))),
    ]


def _matrix_word(word, matrices):
    result = np.eye(len(matrices[0]), dtype=complex)
    for address, sign in word:
        value = matrices[address]
        result = result @ (value if sign == 1 else value.conj().T)
    return result


class SharedPrefixComparison(unittest.TestCase):
    def test_literal_native_shared_prefix_pullback(self):
        c, d = _four_layer_fixture()
        c.append(_wire(S.conj().T, 0) @ _cx(0, 1) @ _wire(H, 1))
        d.append(_wire(T, 0) @ _wire(T.conj().T, 1))
        prefix = c[4] @ d[3] @ c[3] @ d[2]
        inner = c[2] @ d[1]
        for address in range(8):
            b = prefix @ _p(2, address) @ prefix.conj().T
            j_prefix = prefix @ inner
            j = j_prefix @ _p(1, address) @ j_prefix.conj().T
            pullback = (_p(2, address) @ inner @ _p(1, address)
                        @ inner.conj().T)
            np.testing.assert_allclose(b @ j, prefix @ pullback @ prefix.conj().T,
                                       atol=3e-14)
        # This is the retained native-prefix fixture, not a returned query.
        # Its actual inverse is used; no replacement of the source by a
        # Pauli or Clifford is assumed.
        np.testing.assert_allclose(prefix.conj().T @ prefix, np.eye(4), atol=3e-14)

    def test_native_dihedral_group_defect_and_ordered_normal_form(self):
        identity = np.eye(2, dtype=complex)
        z = np.diag([1, -1]).astype(complex)
        representatives = [np.linalg.matrix_power(S, r) @ np.linalg.matrix_power(X, b)
                           for b in range(2) for r in range(4)]
        values = []
        matrices = []
        for address in range(16):
            u, v, s, t = [(address >> index) & 1 for index in range(4)]
            pauli = np.linalg.matrix_power(X, u) @ np.linalg.matrix_power(z, v)
            image = T @ np.linalg.matrix_power(X, s) @ np.linalg.matrix_power(z, t) @ T.conj().T
            matrix = pauli @ image
            coefficients = [abs(np.trace(other.conj().T @ matrix)) / 2
                            for other in representatives]
            value = int(np.argmax(coefficients))
            self.assertAlmostEqual(coefficients[value], 1)
            values.append(value)
            matrices.append(matrix)
        self.assertEqual(set(values), set(range(8)))
        multiply = _dihedral_multiply
        inverse = [_inverse(element, multiply, 8) for element in range(8)]

        def defect(x, y):
            return _product([values[x ^ y], inverse[values[y]], inverse[values[x]]],
                            multiply)

        nontrivial = False
        for x, y, zz in itertools.product(range(16), repeat=3):
            dxy = defect(x, y)
            self.assertEqual(multiply(dxy, dxy), 0)
            self.assertEqual(multiply(dxy, values[zz]), multiply(values[zz], dxy))
            self.assertEqual(defect(x ^ y, zz), multiply(defect(x, zz), defect(y, zz)))
            self.assertEqual(defect(x, y ^ zz), multiply(defect(x, y), defect(x, zz)))
            nontrivial |= multiply(values[x], values[y]) != multiply(values[y], values[x])
        self.assertTrue(nontrivial)
        for address in range(16):
            linear = [values[1 << i] for i in range(4) if address >> i & 1]
            quadratic = [defect(1 << i, 1 << j)
                         for i in range(4) for j in range(i + 1, 4)
                         if address >> i & 1 and address >> j & 1]
            self.assertEqual(_product(linear + quadratic, multiply), values[address])
        # The four Clifford relations have the claimed literal word lengths.
        for x, y, zz in itertools.product(range(0, 16, 3), repeat=3):
            words = _relation_words(x, y, zz)
            self.assertEqual([len(word) for word in words], [8, 6, 9, 9])
            for word in words:
                actual = _matrix_word(word, matrices)
                self.assertAlmostEqual(abs(np.trace(actual)) / 2, 1)
                np.testing.assert_allclose(actual.conj().T @ actual, identity, atol=2e-13)

    def test_faithful_intrinsic_features_and_quadratic_degree_bound(self):
        multiply = _dihedral_multiply
        powers = _ideal_powers(range(8), multiply, 8)
        # D3 is trivial although the complete augmentation ideal I3 is not.
        self.assertTrue(powers[3])
        self.assertEqual([element for element in range(8)
                          if _reduce((1 << element) ^ 1, powers[3]) == 0], [0])
        self.assertEqual(len({_reduce(1 << element, powers[3])
                              for element in range(8)}), 8)

        # A genuine quadratic group map into C2 x C2 can require degree-four
        # intrinsic coordinates at m=3. This illustrates the feature bound;
        # no claim is made that this abstract map is a native returned query.
        xor = lambda left, right: left ^ right
        xor_powers = _ideal_powers(range(4), xor, 8)
        quadratic = [((address & 3) == 3) + 2 * ((address & 12) == 12)
                     for address in range(16)]
        maxima = {}
        for truncation in (1, 2, 3, 8):
            features = [_reduce(1 << value, xor_powers[truncation])
                        for value in quadratic]
            degrees = [_degree([(feature >> bit) & 1 for feature in features])
                       for bit in range(4)]
            maxima[truncation] = max(degrees)
            self.assertLessEqual(max(degrees), 2 * (truncation - 1))
        self.assertEqual(maxima[3], 4)

    def test_fixed_word_error_budget_and_clifford_separation(self):
        identity = np.eye(2, dtype=complex)
        matrices = [identity, X, T @ X @ T.conj().T, S]
        angle = 0.001
        perturbation = np.cos(angle) * identity + 1j * np.sin(angle) * np.diag([1, -1])
        perturbed = [perturbation @ value for value in matrices]
        epsilon = max(np.linalg.norm(left - right, ord=2)
                      for left, right in zip(matrices, perturbed))
        for word in _relation_words(1, 2, 3):
            distance = np.linalg.norm(_matrix_word(word, perturbed)
                                      - _matrix_word(word, matrices), ord=2)
            self.assertLessEqual(distance, len(word) * epsilon + 1e-12)
        radius = np.sqrt(2) / 40
        self.assertLess(2 * 9 * radius, 1 / np.sqrt(2))
        for clifford in _one_qubit_cliffords():
            distance = _projective_distance(clifford, identity)
            if distance > 1e-10:
                self.assertGreaterEqual(distance + 1e-12, 1 / np.sqrt(2))


if __name__ == '__main__':
    unittest.main()
