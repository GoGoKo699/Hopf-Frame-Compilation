"""Small exact receipts for the non-Clifford seed spread screen.

The all-size bounds are analytic. These fixtures retain native phases,
literal affine-query label growth, and the central coordinates used only
in the information-theoretic compression argument.
"""
from fractions import Fraction
from itertools import combinations
from math import comb
import unittest

from tests.exact_matrices import (
    CO, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clifford_orbit_query_rank import (
    _anf, _normal, _one_qubit_cliffords, _rank,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _pauli(width, label):
    mask = (1 << width) - 1
    a, b = label & mask, label >> width
    phase = (a & b).bit_count() % 4
    return _normal(width, label | (phase << (2 * width)))


def _trace(matrix):
    return sum((row[i] for i, row in enumerate(matrix)), CZ)


def _conjugate(unitary, observable):
    return multiply(multiply(unitary, observable), adjoint(unitary))


def _coefficients(width, observable):
    result = {}
    for label in range(1 << (2 * width)):
        value = _trace(multiply(_pauli(width, label), observable)) / (1 << width)
        assert value.imag == 0
        if value != CZ:
            result[label] = value.real
    return result


def _affine_rank(labels):
    labels = list(labels)
    return _rank(label ^ labels[0] for label in labels)


def _spread_seed_word(width):
    # U Z_0 U dagger = tensor_j (T H T dagger)_j. No work initialization.
    return ([('CX', wire, 0) for wire in range(1, width)]
            + [('S', wire) for wire in range(width) for _ in range(3)]
            + [('H', wire) for wire in range(width)]
            + [('T', wire) for wire in range(width)]
            + [('H', wire) for wire in range(width)]
            + [('S', wire) for wire in range(width)]
            + [('T', wire) for wire in range(width)])


def _cx_label(width, label, control, target):
    a, b = label & ((1 << width) - 1), label >> width
    a ^= ((a >> control) & 1) << target
    b ^= ((b >> target) & 1) << control
    return a | (b << width)


def _symplectic(width, left, right):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def _compress_center_and_spectator(density):
    # Keep quantum wire 0, dephase central wire 1, discard spectator wire 2.
    return [[sum((density[row | (h << 2)][column | (h << 2)]
                  for h in range(2)), CZ)
             if (row >> 1) == (column >> 1) else CZ
             for column in range(4)] for row in range(4)]


class NoncliffordSeedSpreadTests(unittest.TestCase):
    def test_extensive_seed_and_every_prefix_support_growth(self):
        for width in range(1, 8):
            word = _spread_seed_word(width)
            for actual in (word, _inverse(word)):
                self.assertEqual(len(actual), 9 * width - 1)
                self.assertEqual(sum(g[0] in ('T', 'TDG') for g in actual),
                                 2 * width)
                bands = _bands(actual)
                self.assertEqual(len(bands), 2)
                for band in bands:
                    self.assertEqual({gate[1] for gate in band}, set(range(width)))
                    self.assertEqual(len(band), width)
                self.assertTrue(all(0 <= wire < width
                                    for gate in actual for wire in gate[1:]))
        for width in (1, 2):
            word = _spread_seed_word(width)
            unitary = native_columns(width, word)
            self.assertEqual(native_columns(width, _inverse(word)), adjoint(unitary))
            seed = _conjugate(unitary, _pauli(width, 1 << width))
            coeffs = _coefficients(width, seed)
            self.assertEqual(len(coeffs), 3 ** width)
            self.assertEqual(_rank(coeffs), 2 * width)
            self.assertEqual(_affine_rank(coeffs), 2 * width)
            self.assertEqual(max(value * value for value in coeffs.values()),
                             Q2(Fraction(1, 2 ** width)))
            expected = eye(1 << width)
            for wire in range(width):
                local = native_columns(width, [('TDG', wire), ('H', wire),
                                               ('T', wire)])
                expected = multiply(local, expected)
            self.assertEqual(seed, expected)
            # The complete signed-Pauli gap comes from exact local trace
            # factors, including the strict non-Clifford coefficient -1/2.
            local_dots = {0: Q2(1), 1: Q2(Fraction(-1, 2)),
                          2: Q2(0), 3: Q2(Fraction(-1, 2))}
            for label in range(1 << (2 * width)):
                pauli = _pauli(width, label)
                dot = _trace(multiply(multiply(multiply(seed, pauli), seed),
                                      pauli)) / (1 << width)
                expected_dot = Q2(1)
                for wire in range(width):
                    factor = (((label >> wire) & 1)
                              | (((label >> (width + wire)) & 1) << 1))
                    expected_dot *= local_dots[factor]
                self.assertEqual(dot, CO * expected_dot)
                if label:
                    self.assertLessEqual(abs(dot.real), Q2(Fraction(1, 2)))
                    for sign in (-1, 1):
                        self.assertGreaterEqual(2 - 2 * sign * dot.real, Q2(1))

        # The plain local inverse admits only full X/Z strings if its
        # conjugated Pauli is to retain the seed's 3**w support cardinality.
        width = 2
        local_word = [gate for gate in _spread_seed_word(width) if gate[0] != 'CX']
        local_unitary = native_columns(width, local_word)
        allowed = []
        for label in range(1 << (2 * width)):
            image = _conjugate(local_unitary, _pauli(width, label))
            coeffs = _coefficients(width, image)
            factors = [((label >> wire) & 1) | (((label >> (width + wire)) & 1) << 1)
                       for wire in range(width)]
            k_xz = sum(factor in (1, 2) for factor in factors)
            k_y = factors.count(3)
            self.assertEqual(len(coeffs), 3 ** k_xz * 2 ** k_y)
            positive, negative = [], []
            for symmetry in range(1 << (2 * width)):
                conjugated = _conjugate(_pauli(width, symmetry), image)
                if conjugated == image:
                    positive.append(symmetry)
                if conjugated == scale(image, -1):
                    negative.append(symmetry)
            self.assertEqual(bool(negative), bool(k_y))
            if not k_y:
                self.assertEqual(_rank(positive), 2 * factors.count(0))
            if len(coeffs) == 3 ** width:
                self.assertTrue(all(factor in (1, 2) for factor in factors))
                self.assertEqual(_affine_rank(coeffs), 2 * width)
                allowed.append(label)
        gram = [sum(_symplectic(width, a, b) << j for j, b in enumerate(allowed))
                for a in allowed]
        self.assertEqual(_rank(gram), 2)

        width, previous_rank, t_count = 2, 0, 0
        word = _spread_seed_word(width)
        for index, gate in enumerate(word, 1):
            prefix = native_columns(width, word[:index])
            coeffs = _coefficients(width, _conjugate(prefix, _pauli(width, 4)))
            rank = _affine_rank(coeffs)
            is_t = gate[0] in ('T', 'TDG')
            t_count += is_t
            if is_t:
                self.assertLessEqual(rank, previous_rank + 1)
            else:
                self.assertEqual(rank, previous_rank)
            self.assertLessEqual(rank, t_count)
            self.assertLessEqual(_rank(coeffs), t_count + 1)
            previous_rank = rank

    def test_two_actual_affine_queries_have_quadratic_label_span(self):
        p = 8
        pairs = list(combinations(range(p), 2))
        width = 1 + p + len(pairs)
        images = []
        z_root = 1 << width
        for address in range(1 << p):
            first = [('CX', 0, 1 + i) for i in range(p) if address >> i & 1]
            second = [('CX', 1 + i, 1 + p + k)
                      for k, (i, j) in enumerate(pairs) if address >> j & 1]
            # Each word is one literal XOR shear between disjoint wire banks.
            for query in (first, second):
                self.assertFalse({g[1] for g in query} & {g[2] for g in query})
            transformed = []
            for label in (1, z_root, 1 | z_root):
                for _, control, target in first + second:
                    label = _cx_label(width, label, control, target)
                transformed.append(label)
            expected = 1 | (address << 1)
            for k, (i, j) in enumerate(pairs):
                expected |= ((address >> i & 1) & (address >> j & 1)) << (1 + p + k)
            self.assertEqual(transformed, [expected, z_root, expected | z_root])
            images.append(expected)
        coefficients = _anf(images)
        self.assertEqual(max(mask.bit_count() for mask, value in enumerate(coefficients)
                             if value), 2)
        monomials = sum(comb(p, j) for j in range(3))
        self.assertEqual(_rank(images), monomials)
        self.assertEqual(_rank(coefficients), monomials)
        span = [1 << wire for wire in range(width)] + [z_root]
        r = _rank(images + [z_root])
        self.assertEqual(r, monomials + 1)
        gram = [sum(_symplectic(width, a, b) << j for j, b in enumerate(span))
                for a in span]
        self.assertEqual(_rank(gram), 2)
        self.assertLessEqual(r, 2 * monomials)
        self.assertLess(2 * monomials, 1 << p)

    def test_central_compression_preserves_statistics_and_squared_error(self):
        word = [('H', 0), ('H', 1), ('CX', 1, 2), ('T', 0),
                ('CX', 0, 1), ('H', 2), ('S', 1)]
        unitary = native_columns(3, word)
        for column in range(8):
            state = [[row[column]] for row in unitary]
            density = multiply(state, adjoint(state))
            compressed = _compress_center_and_spectator(density)
            self.assertEqual(_trace(compressed), CO)
            for a in (0, 1):
                for b in range(4):
                    before = _trace(multiply(_pauli(3, a | (b << 3)), density))
                    after = _trace(multiply(_pauli(2, a | (b << 2)), compressed))
                    self.assertEqual(before, after)

        target, source = _pauli(2, 4), native_columns(2, [('H', 0)])
        difference = add(source, scale(target, -1))
        squared_error = Q2(2, -1)
        self.assertEqual(multiply(difference, difference),
                         scale(eye(4), squared_error))
        for column in range(4):
            sign = 1 - 2 * (column & 1)
            error = (CO - source[column][column] * sign) / 2
            self.assertEqual(error, CO * (squared_error / 4))
            # The central Z bit is preserved, rather than silently discarded.
            self.assertEqual(_pauli(2, 8)[column][column],
                             CO * (1 - 2 * (column >> 1)))

    def test_heavy_coefficient_collisions_retain_all_literal_phases(self):
        word = [(name, 0) for name in 'HTHTHTH']
        unitary = native_columns(1, word)
        self.assertEqual(native_columns(1, _inverse(word)), adjoint(unitary))
        seed = _conjugate(unitary, _pauli(1, 2))
        self.assertEqual(_coefficients(1, seed), {
            1: Q2(Fraction(1, 2)),
            3: Q2(Fraction(-1, 2), Fraction(-1, 4)),
            2: Q2(Fraction(1, 2), Fraction(-1, 4)),
        })
        heavy = scale(_pauli(1, 3), -1)
        mu = Q2(Fraction(3, 8), Fraction(1, 4))
        self.assertGreater(mu, Q2(Fraction(1, 2)))
        orbit = [(_conjugate(c, seed), _conjugate(c, heavy))
                 for c in _one_qubit_cliffords().values()]
        self.assertEqual(len(orbit), 24)
        lower_bound = 2 * (2 * mu - 1)
        for rx, ax in orbit:
            for ry, ay in orbit:
                heavy_product = multiply(ax, ay)
                reverse = multiply(ay, ax)
                actual_sign = 1 if heavy_product == reverse else -1
                self.assertEqual(heavy_product, scale(reverse, actual_sign))
                required_sign = -actual_sign
                defect = add(multiply(rx, ry), scale(multiply(ry, rx),
                                                    -required_sign))
                coefficient = _trace(multiply(adjoint(heavy_product), defect)) / 2
                self.assertEqual(coefficient.imag, 0)
                self.assertGreaterEqual(coefficient.real, lower_bound)


if __name__ == '__main__':
    unittest.main()
