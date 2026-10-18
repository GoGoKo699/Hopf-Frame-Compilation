"""Bounded checks for the outer-axis-normalizer five-layer subcase.

The proof supplies the width-independent degree bound. These checks target
its mod-four sign carry, signed native-axis correction, nonredundancy and
strict dimension thresholds; they do not enumerate general five-layer queries.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _pauli, _wire


def _tau(label, width):
    mask = (1 << width) - 1
    p, q = label & mask, label >> width
    return 1j ** (p & q).bit_count() * _pauli(p, q, width)


def _image(clifford, label, width):
    image = clifford @ _tau(label, width) @ clifford.conj().T
    coefficients = [np.trace(_tau(candidate, width) @ image) / (1 << width)
                    for candidate in range(1 << (2 * width))]
    best = int(np.argmax(np.abs(coefficients)))
    if not np.isclose(abs(coefficients[best]), 1, atol=3e-12):
        raise ValueError("The operator does not have a Clifford Pauli image")
    return best, int(coefficients[best].real < 0)


def _tableau(clifford, width):
    images = [_image(clifford, 1 << column, width)
              for column in range(2 * width)]
    matrix = np.array([[(images[column][0] >> row) & 1
                        for column in range(2 * width)]
                       for row in range(2 * width)], dtype=int)
    return matrix, [item[1] for item in images]


def _add(left, right, scale=1):
    result = dict(left)
    for monomial, coefficient in right.items():
        value = (result.get(monomial, 0) + scale * coefficient) % 4
        if value:
            result[monomial] = value
        else:
            result.pop(monomial, None)
    return result


def _multiply(left, right):
    result = {}
    for a, ca in left.items():
        for b, cb in right.items():
            result = _add(result, {a | b: ca * cb})
    return result


def _parity(variables):
    result = {variable: 1 for variable in variables}
    for left, right in itertools.combinations(variables, 2):
        result = _add(result, {left | right: 2})
    return result


def _phase_polynomial(probe, width):
    """The sign phase F modulo four, as a polynomial in variable S entries."""
    size = 2 * width
    active = [column for column in range(size) if (probe >> column) & 1]

    def variable(row, column):
        return 1 << (row * size + column)

    mask = (1 << width) - 1
    constant = ((probe & mask) & (probe >> width)).bit_count() % 4
    result = {0: constant} if constant else {}
    for wire in range(width):
        for column in active:
            result = _add(result, {variable(wire, column)
                                   | variable(wire + width, column): 1})
        for left, right in itertools.combinations(active, 2):
            result = _add(result, {variable(wire + width, left)
                                   | variable(wire, right): 2})
        p = _parity([variable(wire, column) for column in active])
        q = _parity([variable(wire + width, column) for column in active])
        result = _add(result, _multiply(p, q), -1)
    return result


def _high_bit_polynomial(phase):
    result = set()

    def toggle(monomial):
        if monomial in result:
            result.remove(monomial)
        else:
            result.add(monomial)

    low = []
    for monomial, coefficient in phase.items():
        if coefficient & 2:
            toggle(monomial)
        if coefficient & 1:
            low.append(monomial)
    for left, right in itertools.combinations(low, 2):
        toggle(left | right)
    return result


def _matrix_bits(matrix):
    return sum(int(value) << index for index, value in enumerate(matrix.flat))


def _evaluate(polynomial, values):
    return sum((values & monomial) == monomial for monomial in polynomial) % 2


class OuterAxisNormalizer(unittest.TestCase):
    def test_varying_tableau_sign_uses_mod_four_carries(self):
        fixtures = [np.eye(4, dtype=complex), _cx(0, 1),
                    _wire(H, 0) @ _cx(0, 1) @ _wire(S, 1),
                    _cx(1, 0) @ _wire(S, 0) @ _cx(0, 1) @ _wire(H, 1),
                    _wire(X, 0) @ _wire(S, 1) @ _cx(1, 0) @ _wire(H, 0)]
        carries_needed = False
        for probe in range(16):
            phase = _phase_polynomial(probe, 2)
            self.assertTrue(all(m.bit_count() <= (2 if c & 1 else 3)
                                for m, c in phase.items()))
            sign_polynomial = _high_bit_polynomial(phase)
            self.assertTrue(all(m.bit_count() <= 4 for m in sign_polynomial))
            naive = {m for m, c in phase.items() if c & 2}
            for clifford in fixtures:
                matrix, signs = _tableau(clifford, 2)
                bits = _matrix_bits(matrix)
                phase_value = sum(c for m, c in phase.items()
                                  if (bits & m) == m) % 4
                self.assertEqual(phase_value % 2, 0)
                actual_label, actual_sign = _image(clifford, probe, 2)
                vector = np.array([(probe >> j) & 1 for j in range(4)])
                expected_label = sum(int(bit) << j for j, bit in enumerate(
                    matrix @ vector % 2))
                self.assertEqual(actual_label, expected_label)
                generator_sign = sum(signs[j] for j in range(4)
                                     if (probe >> j) & 1) % 2
                self.assertEqual(actual_sign,
                                 generator_sign ^ _evaluate(sign_polynomial, bits))
                self.assertEqual(_evaluate(sign_polynomial, bits), phase_value // 2)
                carries_needed |= (_evaluate(naive, bits)
                                   != _evaluate(sign_polynomial, bits))
        self.assertTrue(carries_needed)

    def test_signed_axis_permutation_gives_the_native_band_correction(self):
        width = 3
        swap = _cx(0, 1, width) @ _cx(1, 0, width) @ _cx(0, 1, width)
        base = (_wire(H, 2, width) @ _cx(0, 2, width)
                @ _wire(S, 0, width))
        for orientation, permutation, sign_mask in itertools.product(
                itertools.product((0, 1), repeat=2), range(2), range(4)):
            clifford = base @ (swap if permutation else np.eye(8))
            clifford = _pauli(sign_mask, 0, width) @ clifford
            matrix, signs = _tableau(clifford, width)
            delta = list(orientation)
            for source in range(2):
                label = sum(int(matrix[row, width + source]) << row
                            for row in range(2 * width))
                self.assertIn(label, (1 << width, 1 << (width + 1)))
                target = label.bit_length() - 1 - width
                delta[target] ^= orientation[source] ^ signs[width + source]
            for inverse in (False, True):
                layer = np.eye(8, dtype=complex)
                correction = np.eye(8, dtype=complex)
                for wire in range(2):
                    negative = bool(orientation[wire]) ^ inverse
                    layer = _wire(T.conj().T if negative else T, wire, width) @ layer
                    if delta[wire]:
                        correction = _wire(S.conj().T if negative else S,
                                           wire, width) @ correction
                actual = layer @ clifford @ layer.conj().T
                predicted = correction @ clifford
                phase = np.trace(predicted.conj().T @ actual) / 8
                self.assertAlmostEqual(abs(phase), 1)
                np.testing.assert_allclose(actual, phase * predicted, atol=5e-14)
                _tableau(actual, width)
        # An inactive wire and its arbitrary input were retained throughout.
        # The formula is an operator identity, with no initialized helper.

    def test_returned_fixture_has_nonpauli_comparison_and_nonclifford_middle(self):
        identity = np.eye(2, dtype=complex)
        y = 1j * X @ Z
        source = T @ H @ T
        middle = source @ X @ source.conj().T
        comparison = Z @ T.conj().T @ X @ T
        _tableau(comparison, 1)
        self.assertLess(max(abs(np.trace(_tau(label, 1) @ comparison)) / 2
                            for label in range(4)), 1 - 1e-8)
        with self.assertRaises(ValueError):
            _tableau(middle, 1)
        np.testing.assert_allclose(comparison @ Z @ comparison.conj().T,
                                   -Z, atol=3e-14)
        peeled = T.conj().T @ comparison @ T
        _tableau(peeled, 1)
        inner = H @ T @ X @ T.conj().T @ H
        np.testing.assert_allclose(T @ inner @ T.conj().T, middle, atol=3e-14)
        for x, shift in itertools.product(range(2), repeat=2):
            family = [identity, peeled]
            conjugators = [identity, inner]
            np.testing.assert_allclose(
                family[x ^ shift], conjugators[shift] @ family[x]
                @ conjugators[shift].conj().T @ family[shift], atol=4e-14)
        for address in range(2):
            actual = (T @ np.linalg.matrix_power(Z, address) @ T @ H @ T
                      @ np.linalg.matrix_power(X, address)
                      @ np.linalg.matrix_power(y, address) @ T.conj().T
                      @ H @ S.conj().T)
            np.testing.assert_allclose(actual, np.linalg.matrix_power(X, address),
                                       atol=4e-14)
        self.assertEqual((4, 4, 3 + 1 + 5 + 4 + 1, 0, 0), (4, 4, 14, 0, 0))
        # Its fair simplified baseline remains one CNOT, with 0T and 1G.

    def test_span_mixing_peel_is_third_level_but_not_clifford(self):
        layer = _wire(T, 0) @ _wire(T.conj().T, 1)
        clifford = _cx(0, 1)
        matrix, _ = _tableau(clifford, 2)
        active_images = [sum(int(matrix[row, 2 + column]) << row
                             for row in range(4)) for column in range(2)]
        self.assertEqual(active_images, [4, 12])
        self.assertTrue(all(label & 3 == 0 for label in active_images))
        self.assertNotIn(active_images[1], (4, 8))
        peeled = layer.conj().T @ clifford @ layer
        with self.assertRaises(ValueError):
            _tableau(peeled, 2)
        # All Hermitian Pauli images are Cliffords, which is exactly C3.
        for label in range(16):
            image = peeled @ _tau(label, 2) @ peeled.conj().T
            _tableau(image, 2)
        # The special conjugated normalizer subgroup is closed; general
        # closure of the whole third level is never used.
        second = _wire(X, 1) @ _wire(S, 0) @ _cx(1, 0)
        second_peeled = layer.conj().T @ second @ layer
        product = peeled @ second_peeled
        np.testing.assert_allclose(product,
                                   layer.conj().T @ clifford @ second @ layer,
                                   atol=3e-14)
        for label in range(16):
            _tableau(product @ _tau(label, 2) @ product.conj().T, 2)

    def test_dimension_thresholds_are_strict_at_the_claimed_boundary(self):
        d = 122
        for n in (2 * d, 2 * d + 1, 2 * d + 2):
            dimension = (sum(math.comb(n, j) for j in range(d + 1))
                         + sum(math.comb(n, j) for j in range(d)))
            self.assertEqual(dimension < 2 ** n, n > 2 * d)
        d = 20
        for n in (2 * d + 1, 2 * d + 2, 2 * d + 3):
            dimension = 2 * sum(math.comb(n, j) for j in range(d + 1))
            self.assertEqual(dimension < 2 ** n, n > 2 * d + 1)


if __name__ == "__main__":
    unittest.main()
