"""Exact bounded evidence for the native span-preserving degree theorem.

The analytic proof establishes the width-independent bound. These tests
exercise its mixed-axis nilpotence, subgroup filtration distinction,
intrinsic feature action and target detection; they are not a synthesis
claim or a finite-width extrapolation.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.test_outer_axis_normalizer import _tableau
from tests.test_read_only_three_layer_rigidity import H, T, X, _cx, _wire


def _basis(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in pivots:
                vector ^= pivots[pivot]
            else:
                pivots[pivot] = vector
                break
    return pivots


def _reduce(vector, pivots):
    for pivot in sorted(pivots, reverse=True):
        if vector & (1 << pivot):
            vector ^= pivots[pivot]
    return vector


def _ring_product(left, right, multiply):
    result = 0
    while left:
        low = left & -left
        a = low.bit_length() - 1
        left ^= low
        remaining = right
        while remaining:
            low_right = remaining & -remaining
            b = low_right.bit_length() - 1
            remaining ^= low_right
            result ^= 1 << multiply(a, b)
    return result


def _ideal_powers(elements, multiply, depth):
    augmentation = _basis((1 << element) ^ 1 for element in elements)
    powers = [_basis(1 << element for element in elements), augmentation]
    while len(powers) <= depth:
        powers.append(_basis(
            _ring_product(left, right, multiply)
            for left in powers[-1].values()
            for right in augmentation.values()))
    return powers


def _dihedral_multiply(left, right):
    a, b = left % 4, left // 4
    c, d = right % 4, right // 4
    return (a + (-1 if b else 1) * c) % 4 + 4 * (b ^ d)


def _closure(generators, multiply):
    result = {0, *generators}
    while True:
        expanded = result | {multiply(a, b) for a in result for b in result}
        if expanded == result:
            return frozenset(result)
        result = expanded


def _subgroups(size, multiply):
    found = {frozenset({0})}
    todo = list(found)
    while todo:
        subgroup = todo.pop()
        for element in range(size):
            candidate = _closure(subgroup | {element}, multiply)
            if candidate not in found:
                found.add(candidate)
                todo.append(candidate)
    return found


def _symplectic_pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def _increment(axis, width):
    return np.array([
        [((axis >> row) & 1) * _symplectic_pair(axis, 1 << column, width)
         for column in range(2 * width)]
        for row in range(2 * width)], dtype=int)


class SpanPreservingAugmentation(unittest.TestCase):
    def test_native_mixed_axes_have_nonzero_cubic_but_zero_quartic_products(self):
        width = 2
        identity = np.eye(2 * width, dtype=int)
        z0, z1, x1 = 4, 8, 2
        h0, h1 = _increment(z0, width), _increment(z1, width)
        n = (_increment(x1, width) + _increment(x1 ^ z0, width)) % 2

        # The two independent outer axes commute. Their combined shear
        # preserves L=<Z0,Z1>, although each separate X1 shear does not.
        self.assertEqual(_symplectic_pair(x1, x1 ^ z0, width), 0)
        self.assertTrue(np.array_equal(n[:, 2], np.zeros(4, dtype=int)))
        self.assertTrue(np.array_equal(n[:, 3], np.array([0, 0, 1, 0])))
        self.assertFalse(np.array_equal(_increment(x1, width)[:, 3],
                                       np.zeros(4, dtype=int)))
        np.testing.assert_array_equal(n @ h1 @ n % 2, h0)
        self.assertTrue(np.any(n @ h1 @ n % 2))
        for word in itertools.product((h0, h1, n), repeat=4):
            product = identity.copy()
            for factor in word:
                product = product @ factor % 2
            self.assertFalse(np.any(product))
        symplectic = (identity + h1) @ (identity + n) % 2
        self.assertFalse(np.array_equal(symplectic @ symplectic % 2, identity))
        np.testing.assert_array_equal(
            np.linalg.matrix_power(symplectic, 4) % 2, identity)

        # Verify these are actual one-T Pauli images, retaining the fixed
        # Clifford and its actual inverse. These analytical-source fixtures
        # are not returned queries or proposed resource savings.
        layer = _wire(T, 0, width) @ _wire(T, 1, width)
        coordinates = _wire(H, 1, width) @ _cx(1, 0, width)
        source_a = layer @ _wire(X, 1, width) @ layer.conj().T
        source_n = (coordinates @ layer @ _wire(X, 0, width)
                    @ _wire(X, 1, width) @ layer.conj().T
                    @ coordinates.conj().T)
        np.testing.assert_array_equal(_tableau(source_a, width)[0], identity + h1)
        np.testing.assert_array_equal(_tableau(source_n, width)[0], identity + n)

    def test_subgroup_and_ambient_filtrations_cannot_be_identified(self):
        multiply = lambda left, right: (left + right) % 4
        ambient = _ideal_powers(range(4), multiply, 4)
        subgroup = _ideal_powers((0, 2), multiply, 3)
        target_difference = (1 << 2) ^ 1
        self.assertEqual(_reduce(target_difference, ambient[2]), 0)
        self.assertNotEqual(_reduce(target_difference, subgroup[2]), 0)
        # w=2,m=2 gives the correct ambient truncation 3, not 2.
        for vector in (0, 1, 1 << 2, target_difference):
            if _reduce(vector, ambient[3]) == 0:
                self.assertEqual(_reduce(vector, subgroup[2]), 0)

    def test_extension_lemma_for_every_subgroup_of_nonabelian_d8(self):
        multiply = _dihedral_multiply
        ambient = _ideal_powers(range(8), multiply, 8)
        self.assertNotEqual(multiply(1, 4), multiply(4, 1))
        dimensions = [
            {element for element in range(8)
             if _reduce((1 << element) ^ 1, ambient[degree]) == 0}
            for degree in (1, 2, 3)]
        self.assertEqual(dimensions, [set(range(8)), {0, 2}, {0}])
        subgroups = _subgroups(8, multiply)
        self.assertEqual(len(subgroups), 10)
        for group in subgroups:
            elements = sorted(group)
            intrinsic = _ideal_powers(elements, multiply, 4)
            for truncation in range(1, 5):
                extension_depth = 2 * (truncation - 1) + 1
                for mask in range(1 << len(elements)):
                    vector = sum(1 << element for bit, element in enumerate(elements)
                                 if mask & (1 << bit))
                    if _reduce(vector, ambient[extension_depth]) == 0:
                        self.assertEqual(_reduce(vector, intrinsic[truncation]), 0)

    def test_intrinsic_quotient_retains_automorphism_and_right_action(self):
        multiply = _dihedral_multiply
        powers = _ideal_powers(range(8), multiply, 4)

        def automorphism(element):
            a, b = element % 4, element // 4
            return (-a + b) % 4 + 4 * b

        for left, right in itertools.product(range(8), repeat=2):
            self.assertEqual(automorphism(multiply(left, right)),
                             multiply(automorphism(left), automorphism(right)))
        self.assertEqual(len({automorphism(element) for element in range(8)}), 8)
        for degree in (2, 3, 4):
            for vector in powers[degree].values():
                image = sum(1 << automorphism(element) for element in range(8)
                            if vector & (1 << element))
                self.assertEqual(_reduce(image, powers[degree]), 0)
                for right in range(8):
                    translated = _ring_product(image, 1 << right, multiply)
                    self.assertEqual(_reduce(translated, powers[degree]), 0)

    def test_native_word_degree_and_actual_inverse_target_detection(self):
        multiply = _dihedral_multiply
        powers = _ideal_powers(range(8), multiply, 3)
        factors = (4, 5, 4, 6, 2)
        self.assertTrue(all(multiply(factor, factor) == 0 for factor in factors))
        values, baselines = [], []
        for address in range(1 << len(factors)):
            word = 0
            selected = [factor for bit, factor in enumerate(factors)
                        if address & (1 << bit)]
            for factor in selected:
                word = multiply(word, factor)
            baselines.append(word)
            values.append(_reduce(1 << word, powers[3]))
            inverse = 0
            for factor in reversed(selected):
                inverse = multiply(inverse, factor)
            self.assertEqual(multiply(inverse, word), 0)
            boolean_target = int(address & 7 == 7)
            actual = multiply(word, 2 if boolean_target else 0)
            self.assertEqual(multiply(inverse, actual), 2 if boolean_target else 0)
        # Vector-valued Boolean Möbius transform checks the polynomial
        # degree after truncation, without choosing quotient coordinates.
        coefficients = values.copy()
        for bit in range(len(factors)):
            for mask in range(len(coefficients)):
                if mask & (1 << bit):
                    coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
        self.assertTrue(any(value for mask, value in enumerate(coefficients)
                            if mask.bit_count() == 2))
        self.assertTrue(all(value == 0 for mask, value in enumerate(coefficients)
                            if mask.bit_count() >= 3))
        central_target = (1 << 2) ^ 1
        self.assertEqual(_reduce(central_target, powers[2]), 0)
        self.assertNotEqual(_reduce(central_target, powers[3]), 0)

    def test_dimension_threshold_and_jennings_exponent_ledger(self):
        # gamma_1 exponent8; gamma_2,3 exponent4; gamma_4,...,7 exponent2.
        exponents = {1: 8, 2: 4, 3: 4, 4: 2, 5: 2, 6: 2, 7: 2}
        for index, exponent in exponents.items():
            power = 1
            while index * power < 8:
                power *= 2
            self.assertEqual(power % exponent, 0)
        self.assertEqual(7 * (8 - 1) + 1, 50)
        self.assertEqual((50 - 1) + (8 - 1), 56)

        def capacity(variables):
            return (sum(math.comb(variables, degree) for degree in range(57))
                    + sum(math.comb(variables, degree) for degree in range(56)))

        self.assertEqual(capacity(112), 1 << 112)
        self.assertLess(capacity(113), 1 << 113)
        for variables in range(114, 130):
            self.assertLess(capacity(variables), 1 << variables)


if __name__ == "__main__":
    unittest.main()
