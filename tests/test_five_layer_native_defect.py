"""Bounded evidence for the general read-only five-layer native defect.

The tests isolate the bilinear pullback, its one-band matching, literal
signs on moving probes, and the translation-space obstruction. Native
prefix examples below are not advertised as returned-query syntheses.
The proof, rather than these finite examples, supplies arbitrary-width
degree and robustness conclusions.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.test_one_band_axis_matching import (
    _one_qubit_cliffords, _projective_distance,
)
from tests.test_outer_axis_normalizer import _image, _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _wire


def _anf(values):
    coefficients = list(map(int, values))
    for bit in range(len(values).bit_length() - 1):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return {mask for mask, coefficient in enumerate(coefficients) if coefficient}


def _xor(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def _multiply(left, right):
    result = set()
    for a in left:
        for b in right:
            result.symmetric_difference_update({a | b})
    return result


def _substitute(polynomial, variables):
    result = set()
    for monomial in polynomial:
        term = {0}
        for index, variable in enumerate(variables):
            if monomial >> index & 1:
                term = _multiply(term, variable)
        result = _xor(result, term)
    return result


def _evaluate(polynomial, address):
    return sum(address & term == term for term in polynomial) % 2


def _degree(polynomial):
    return max((term.bit_count() for term in polynomial), default=-1)


def _projective_pauli_label(operator, width):
    coefficients = [abs(np.trace(_tau(label, width) @ operator)) / (1 << width)
                    for label in range(1 << (2 * width))]
    label = int(np.argmax(coefficients))
    if not np.isclose(coefficients[label], 1, atol=2e-12):
        raise ValueError("The operator is not projectively Pauli")
    return label


class FiveLayerNativeDefect(unittest.TestCase):
    def test_native_defect_is_a_bilinear_pauli_pullback(self):
        width = 2
        inner = _cx(0, 1) @ _wire(H, 1) @ _wire(T, 0) @ _wire(T.conj().T, 1)
        prefix = (_wire(T, 0) @ _wire(T.conj().T, 1)
                  @ _wire(H, 0) @ _cx(1, 0) @ _wire(T, 1))
        p2 = [_tau(address, width) for address in range(16)]
        # This invertible binary rearrangement makes both X and Z labels
        # vary, including anticommuting literal Pauli representatives.
        p1 = [_tau(((address << 1) & 15) | (address >> 3), width)
              for address in range(16)]
        images = [inner @ value @ inner.conj().T for value in p1]
        values = [pauli @ image for pauli, image in zip(p2, images)]
        labels = {}
        nonzero = False
        for x, y in itertools.product(range(16), repeat=2):
            defect = values[x ^ y] @ values[y].conj().T @ values[x].conj().T
            commutator = images[x] @ p2[y] @ images[x].conj().T @ p2[y].conj().T
            self.assertLess(_projective_distance(defect, commutator), 5e-12)
            labels[x, y] = _projective_pauli_label(defect, width)
            nonzero |= labels[x, y] != 0
            native = prefix @ defect @ prefix.conj().T
            dressed = [prefix @ values[address] @ prefix.conj().T
                       for address in (x ^ y, y, x)]
            np.testing.assert_allclose(native, dressed[0] @ dressed[1].conj().T
                                       @ dressed[2].conj().T, atol=8e-14)
        self.assertTrue(nonzero)
        for x, y, z in itertools.product(range(16), repeat=3):
            self.assertEqual(labels[x ^ y, z], labels[x, z] ^ labels[y, z])
            self.assertEqual(labels[x, y ^ z], labels[x, y] ^ labels[x, z])

    def test_special_defect_matches_even_when_middle_is_nonclifford(self):
        prefix = T @ H @ T
        middle = prefix @ X @ prefix.conj().T
        defect = prefix @ Z @ prefix.conj().T
        with self.assertRaises(ValueError):
            _tableau(middle, 1)
        # The defect's pulled-back Pauli is Z; the middle instead dresses X.
        # Only the former has the matching one-band input in this example.
        inner_defect = H @ T @ Z @ T.conj().T @ H
        self.assertEqual(_image(inner_defect, 2, 1)[0], 2)
        _tableau(defect, 1)
        np.testing.assert_allclose(defect, T @ X @ T.conj().T, atol=2e-14)

    def test_one_band_clifford_gap_rejects_unmatched_defects(self):
        cliffords = _one_qubit_cliffords()
        matched = unmatched = 0
        for negative, label in itertools.product((False, True), range(4)):
            layer = T.conj().T if negative else T
            inner = H @ T @ _tau(label, 1) @ T.conj().T @ H
            actual = layer @ inner @ layer.conj().T
            matching = _image(inner, 2, 1)[0] == 2
            if matching:
                matched += 1
                _tableau(actual, 1)
            else:
                unmatched += 1
                with self.assertRaises(ValueError):
                    _tableau(actual, 1)
                distance = min(_projective_distance(actual, other) for other in cliffords)
                self.assertGreaterEqual(distance + 1e-12, math.sqrt(2) / 4)
        self.assertTrue(matched and unmatched)

    def test_signed_axis_swap_and_quadratic_moving_probe(self):
        # A one-T Pauli image whose symplectic action swaps the two active
        # outer Z axes. This tests genuine mixing, orientation and signs.
        coordinate = _cx(0, 1) @ _wire(H, 0) @ _cx(0, 1)
        inner_layer = _wire(T, 0) @ _wire(T, 1)
        outer_layer = _wire(T, 0) @ _wire(T.conj().T, 1)
        native = []
        signs = []
        matrices = []
        correction_polynomials = [[], []]
        swapped = False
        for family in range(8):
            s, t, r = [(family >> bit) & 1 for bit in range(3)]
            label = s | (s << 1) | (t << 2) | (r << 3)
            inner = coordinate @ inner_layer @ _tau(label, 2) @ inner_layer.conj().T @ coordinate.conj().T
            tableau, generator_sign = _tableau(inner, 2)
            delta = [0, 1]
            for source in range(2):
                image, sign = _image(inner, 1 << (2 + source), 2)
                self.assertIn(image, (4, 8))
                target = image.bit_length() - 3
                swapped |= target != source
                delta[target] ^= source ^ sign
            correction = np.eye(4, dtype=complex)
            for target in range(2):
                correction_polynomials[target].append(delta[target])
                if delta[target]:
                    correction = _wire(S if target == 0 else S.conj().T, target) @ correction
            actual = outer_layer @ inner @ outer_layer.conj().T
            predicted = correction @ inner
            self.assertLess(_projective_distance(actual, predicted), 5e-12)
            native.append(actual)
            matrix, _ = _tableau(actual, 2)
            matrices.append(matrix)
            signs.append([_image(actual, probe, 2)[1] for probe in range(16)])
        self.assertTrue(swapped)
        self.assertTrue(any(correction_polynomials[0]))
        self.assertTrue(any(correction_polynomials[1]))
        for row, column in itertools.product(range(4), repeat=2):
            self.assertLessEqual(_degree(_anf([matrix[row, column] for matrix in matrices])), 5)

        # The three family bits and four probe bits form a seven-variable
        # truth table, obtained from literal full-operator Pauli images.
        sign_polynomial = _anf([signs[value & 7][value >> 3] for value in range(128)])
        variables = [{1 << index} for index in range(3)]
        for index in range(4):
            left, right = 3 + 2 * index, 4 + 2 * index
            variables.append({1 << left, (1 << left) | (1 << right)})
        moving_sign = _substitute(sign_polynomial, variables)
        self.assertLessEqual(_degree(moving_sign), 10)
        self.assertGreater(_degree(moving_sign), 2)
        for address in range(1 << 11):
            family = address & 7
            probe = sum(_evaluate(variables[3 + index], address) << index for index in range(4))
            self.assertEqual(_evaluate(moving_sign, address), signs[family][probe])

    def test_actual_inverse_error_budget_before_projective_rounding(self):
        exact = [H, T @ X @ T.conj().T, S]
        perturbed = []
        for index, unitary in enumerate(exact):
            angle = (index + 1) / 1000
            error = math.cos(angle) * np.eye(2) + 1j * math.sin(angle) * Z
            perturbed.append(error @ unitary)
        native_defect = exact[0] @ exact[1].conj().T @ exact[2].conj().T
        actual_defect = perturbed[0] @ perturbed[1].conj().T @ perturbed[2].conj().T
        error_sum = sum(np.linalg.norm(left - right, ord=2)
                        for left, right in zip(exact, perturbed))
        self.assertLessEqual(np.linalg.norm(actual_defect - native_defect, ord=2), error_sum + 1e-13)
        radius = math.sqrt(2) / 24
        self.assertAlmostEqual(6 * radius, math.sqrt(2) / 4)
        self.assertLess(6 * radius, 1 / math.sqrt(2))

    def test_visible_translation_space_excludes_degree_fifteen(self):
        # Full-monomial translates have full rank: each translate here is
        # the indicator of exactly one point, with a permuted support.
        variables = 6
        translates = []
        for shift in range(1 << variables):
            values = [int((address ^ shift) == (1 << variables) - 1)
                      for address in range(1 << variables)]
            translates.append(sum(value << address for address, value in enumerate(values)))
        self.assertEqual(set(translates), {1 << bit for bit in range(1 << variables)})

        def dimension_bound(degree):
            return sum(math.comb(degree, j) for j in range(8)) + sum(math.comb(degree, j) for j in range(7))

        self.assertEqual(dimension_bound(14), 1 << 14)
        for degree in range(15, 65):
            self.assertLess(dimension_bound(degree), 1 << degree)

    def test_free_outer_baseline_has_an_incompatible_high_degree_defect(self):
        # An independent outer comparison can be specified for any Boolean
        # g: S_K=(I+a N_Z)(I+(b+g)N_X). The native shared-prefix restriction
        # is essential. For a shift outside a,b its ordered defect has
        # upper-right entry Delta_y g, which has degree14 for this g.
        variables = 15
        full = (1 << variables) - 1
        shift = 1 << 2
        entries = []
        for address in range(1 << variables):
            b = (address >> 1) & 1
            target = int(address == full)
            translated = int((address ^ shift) == full)
            current_upper = b ^ target
            next_upper = b ^ translated
            # Multiply the upper row of S_K(x+y) by the right column
            # of S_K(x)^-1. S_K(y)=I for the chosen shift.
            entries.append(current_upper ^ next_upper)
            # The visible extractor recovers g even with the lower
            # shear present; it cannot supply the missing native defect.
            self.assertEqual(current_upper ^ b, target)
        self.assertEqual(_degree(_anf(entries)), 14)
        self.assertGreater(_degree(_anf(entries)), 5)


if __name__ == '__main__':
    unittest.main()
