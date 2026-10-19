"""Exact bounded fixtures for two restricted approximate indicator routes.

The deleted-parity filter attains the sharp sqrt(2) complete phase error.
The hash fixture pays its native arithmetic, including actual inverses,
before comparing initialized and arbitrary seed inputs. Neither fixture
claims a lower bound for general coherent indicator circuits.
"""
from fractions import Fraction
from itertools import combinations
from math import comb
import unittest

import numpy as np

from tests.exact_matrices import CI, CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_majorana_walsh_query import _inverse as _native_inverse
from tests.test_parallel_bilinear_word_query import (
    _inverse as _phase_inverse, _monomial, _parity_phase,
)


def _interpolate(values):
    """Rational monomial coefficients, interpolating nodes 0,...,n-1."""
    result = [Fraction(0)] * len(values)
    for node, value in enumerate(values):
        numerator, denominator = [Fraction(1)], 1
        for other in range(len(values)):
            if other == node:
                continue
            expanded = [Fraction(0)] * (len(numerator) + 1)
            for power, coefficient in enumerate(numerator):
                expanded[power] -= other * coefficient
                expanded[power + 1] += coefficient
            numerator = expanded
            denominator *= node - other
        for power, coefficient in enumerate(numerator):
            result[power] += value * coefficient / denominator
    return result


def _evaluate(coefficients, argument):
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * argument + coefficient
    return value


def _native_phase(word):
    return [(('T' if gate[2] == 1 else 'TDG'), gate[1])
            if gate[0] == 'T' else gate for gate in word]


def _phase_native(word):
    return [('T', gate[1], 1 if gate[0] == 'T' else -1)
            if gate[0] in ('T', 'TDG') else gate for gate in word]


def _counts(word):
    return (sum(gate[0] in ('T', 'TDG') for gate in word),
            sum(gate[0] not in ('T', 'TDG') for gate in word))


def _hash_word(p, k):
    address = list(range(p))
    seed = [[p + row * p + column for column in range(p)] for row in range(k)]
    target = list(range(p + p * k, p + p * k + k))
    dirty = [[p + p * k + k + row * p + column for column in range(p)]
             for row in range(k)]
    layers = [[], [], [], []]
    linear = []
    for row in range(k):
        for column in range(p):
            a, t, d = seed[row][column], target[row], dirty[row][column]
            for layer, controls, power in (
                    (0, (), 1), (1, (a,), -1), (2, (t,), -1), (3, (a, t), 1)):
                layers[layer] += _parity_phase(d, controls, power)
            linear.append(('CX', address[column], d))
    relative = [gate for layer in layers for gate in layer]
    phase = _phase_inverse(relative) + linear + relative + _phase_inverse(linear)
    hadamards = [('H', wire) for wire in target]
    return (hadamards + _native_phase(phase) + hadamards,
            layers, address, seed, target, dirty)


def _zero_phase_word(k):
    """Native phase (-1)^(z [t=0]) on local target wires 0,...,k-1,z=k."""
    word = [('S', k), ('S', k)]
    for wire in range(k):
        word += [('H', k), ('CX', wire, k), ('H', k)]
    if k == 2:
        # 4abc = a+b+c-(a xor b)-(a xor c)-(b xor c)+(a xor b xor c).
        for size in (1, 2, 3):
            for subset in combinations(range(3), size):
                word += _native_phase(_parity_phase(
                    subset[-1], subset[:-1], 1 if size % 2 else -1))
    elif k != 1:
        raise ValueError('This bounded fixture only compiles one or two hash bits.')
    return word


def _hash(seed, address, p, k):
    return sum((((seed >> (row * p)) & ((1 << p) - 1) & address).bit_count() % 2)
               << row for row in range(k))


def _norm_squared(vector):
    return sum((entry.abs2() for entry in vector), 0)


class ApproximateIndicatorFilterLimitsTests(unittest.TestCase):
    def test_deleted_parity_polynomial_and_every_subset_attain_sqrt_two(self):
        for p in range(1, 5):
            m = 1 << p
            values = [Fraction(int(weight == 0)) - Fraction((-1) ** weight, m)
                      for weight in range(p + 1)]
            polynomial = _interpolate(values[:-1])
            self.assertEqual(len(polynomial), p)
            self.assertNotEqual(polynomial[-1], 0)  # degree exactly p-1
            self.assertEqual([_evaluate(polynomial, weight) for weight in range(p + 1)],
                             values)
            differences = values
            for _ in range(p):
                differences = [right - left for left, right
                               in zip(differences, differences[1:])]
            self.assertEqual(differences, [Fraction(0)])
            self.assertEqual(sum((-1) ** weight * comb(p, weight) * values[weight]
                                 for weight in range(p + 1)), 0)

            errors = [_evaluate(polynomial, leaf.bit_count()) - int(leaf == 0)
                      for leaf in range(m)]
            self.assertEqual(errors, [Fraction(-(-1) ** leaf.bit_count(), m)
                                      for leaf in range(m)])
            self.assertEqual(sum(abs(error) for error in errors), 1)
            self.assertEqual(sum((-1) ** leaf.bit_count() * error
                                 for leaf, error in enumerate(errors)), -1)
            # Enumerate every one of the 2^m arbitrary bank inputs. Integer
            # numerators avoid approximating the small phase angles.
            subset_sums = [0]
            for error in errors:
                numerator = int(m * error)
                subset_sums += [subtotal + numerator for subtotal in subset_sums]
            self.assertEqual(len(subset_sums), 1 << m)
            self.assertEqual((min(subset_sums), max(subset_sums)), (-m // 2, m // 2))
            self.assertTrue(all(abs(subtotal) <= m // 2 for subtotal in subset_sums))
            # All phase errors lie in [-pi/2,pi/2]; the two endpoint ratios
            # are exactly +/-i, certifying the attained squared error 2.
            self.assertEqual((CI - CO).abs2(), 2)
            self.assertEqual((-CI - CO).abs2(), 2)

    def test_native_hash_shear_clean_seed_dirty_seed_and_whole_bank(self):
        p = 2
        for k in (1, 2):
            word, layers, address, seed, target, dirty = _hash_word(p, k)
            width = p + 2 * p * k + k
            values = np.arange(1 << width, dtype=np.int64)
            private = {wire for row in dirty for wire in row}
            for layer in layers:
                phase_targets = [gate[1] for gate in layer if gate[0] == 'T']
                self.assertEqual(set(phase_targets), private)
                self.assertEqual(len(phase_targets), p * k)
            expected = np.zeros_like(values)
            for row in range(k):
                for column in range(p):
                    expected ^= (((values >> seed[row][column]) & 1)
                                 & ((values >> address[column]) & 1)
                                 & ((values >> target[row]) & 1))
            expected *= 4
            for orientation in (word, _native_inverse(word)):
                self.assertTrue(all(gate[0] in ('H', 'CX', 'T', 'TDG')
                                    for gate in orientation))
                self.assertEqual(_counts(orientation), (8 * p * k, 18 * p * k + 2 * k))
                # Remove only the paid H endpoints. All remaining basis
                # phases and returns are evaluated from literal gates.
                images, phases = _monomial(_phase_native(orientation[k:-k]), width)
                np.testing.assert_array_equal(images, values)
                np.testing.assert_array_equal(phases, expected)
                # Fourier-conjugate each complete target block. This checks
                # all output columns and all dirty values without a dense
                # 4096-by-4096 matrix, using the evaluated native phases.
                target_mask = sum(1 << wire for wire in target)
                for basis in range(1 << width):
                    base = basis & ~target_mask
                    incoming = sum(((basis >> wire) & 1) << row
                                   for row, wire in enumerate(target))
                    a = sum(((basis >> wire) & 1) << (row * p + column)
                            for row in range(k) for column, wire in enumerate(seed[row]))
                    x = basis & ((1 << p) - 1)
                    for outgoing in range(1 << k):
                        numerator = 0
                        for frequency in range(1 << k):
                            index = base | sum(((frequency >> row) & 1) << wire
                                               for row, wire in enumerate(target))
                            sign = (((incoming ^ outgoing) & frequency).bit_count()
                                    + int(phases[index]) // 4) % 2
                            numerator += (-1) ** sign
                        self.assertEqual(numerator, (1 << k)
                                         if outgoing == incoming ^ _hash(a, x, p, k) else 0)

            zero_word = _zero_phase_word(k)
            zero = native_columns(k + 1, zero_word)
            zero_inverse = native_columns(k + 1, _native_inverse(zero_word))
            self.assertEqual(zero_inverse, adjoint(zero))
            self.assertEqual(_counts(zero_word), (0, 5) if k == 1 else (7, 18))
            for row in range(1 << (k + 1)):
                for column in range(1 << (k + 1)):
                    self.assertEqual(zero[row][column],
                                     CO * (-1 if row == 1 << k else 1) if row == column else CZ)

            seed_wires = [wire for row in seed for wire in row]
            preparation = [('H', wire) for wire in seed_wires]
            mapped_zero = [(name, *(target[wire] if wire < k else width for wire in wires))
                           for name, *wires in zero_word]
            complete_word = (preparation + word + mapped_zero + _native_inverse(word)
                             + _native_inverse(preparation))
            for orientation in (complete_word, _native_inverse(complete_word)):
                self.assertEqual(_counts(orientation), (32, 85) if k == 1 else (71, 178))
                self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                    for gate in orientation))

            # The verified shear and its actual inverse return all hash
            # targets initialized at zero and every arbitrary private d.
            # The only remaining action is this exact small seed operator.
            seed_size = 1 << (p * k)
            h_word = [('H', wire) for wire in range(p * k)]
            h = native_columns(p * k, h_word)
            h_inverse = native_columns(p * k, _native_inverse(h_word))
            self.assertEqual(h_inverse, adjoint(h))
            self.assertEqual(multiply(h_inverse, h), eye(seed_size))
            for x in range(1 << p):
                for z in (0, 1):
                    diagonal = [zero[_hash(a, x, p, k) | (z << k)]
                                    [_hash(a, x, p, k) | (z << k)]
                                for a in range(seed_size)]
                    inverse_diagonal = [zero_inverse[_hash(a, x, p, k) | (z << k)]
                                                     [_hash(a, x, p, k) | (z << k)]
                                        for a in range(seed_size)]
                    actual = multiply(h_inverse, [[phase * entry for entry in row]
                                                  for phase, row in zip(diagonal, h)])
                    inverse = multiply(h_inverse, [[phase * entry for entry in row]
                                                   for phase, row in zip(inverse_diagonal, h)])
                    self.assertEqual(inverse, adjoint(actual))
                    ideal_phase = -1 if x == 0 and z else 1
                    clean_error = [actual[row][0] - (CO * ideal_phase if row == 0 else CZ)
                                   for row in range(seed_size)]
                    self.assertEqual(_norm_squared(clean_error),
                                     Fraction(4, 1 << k) if x and z else 0)
                    if x and z:
                        # An allowed dirty input H|A=0> becomes the all-zero
                        # bad seed after preparation. It is returned with
                        # the wrong sign, even though clean-seed error falls.
                        bad_seed = [[row[0]] for row in h_inverse]
                        output = multiply(actual, bad_seed)
                        self.assertEqual(output, [[-entry[0]] for entry in bad_seed])
                        self.assertEqual(_norm_squared([left[0] - right[0]
                                                       for left, right in zip(output, bad_seed)]), 4)

            if k < p:
                for x in range(1 << p):
                    bank_diagonal = []
                    for a in range(seed_size):
                        phase = CO
                        accepted = 0
                        for leaf in range(1 << p):
                            hashed = _hash(a, x ^ leaf, p, k)
                            accepted += hashed == 0
                            phase *= zero[hashed | (1 << k)][hashed | (1 << k)]
                        self.assertIn(accepted, (2, 4))
                        self.assertEqual(phase, CO)
                        bank_diagonal.append(phase)
                    # All Y bits are one: each even collision bucket gives
                    # +1 for every seed, whereas the exact indicator is -1.
                    bank = multiply(h_inverse, [[phase * entry for entry in row]
                                                for phase, row in zip(bank_diagonal, h)])
                    self.assertEqual(bank, eye(seed_size))
                    self.assertEqual(_norm_squared([bank[row][0] + (CO if row == 0 else CZ)
                                                    for row in range(seed_size)]), 4)


if __name__ == '__main__':
    unittest.main()
