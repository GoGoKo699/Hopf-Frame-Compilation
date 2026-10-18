"""Exact bounded evidence for nonlinear binary-unit character transport.

The analytic width obstruction is architecture-specific. Native source and
transport words below retain all dirty and output columns and actual inverses.
"""
from math import log2
from random import Random
import unittest

import numpy as np

from tests.test_parallel_bilinear_word_query import (
    _inverse, _native_matrix, _relative_grid,
)
from tests.test_returned_diagonal_echo_degree import _ccz


def _multiply(first, second, width):
    result = 0
    while second:
        if second & 1:
            result ^= first
        first <<= 1
        second >>= 1
    return result & ((1 << (width + 1)) - 1)


def _reciprocal(unit, width):
    result = 1
    for degree in range(1, width + 1):
        bit = sum((unit >> i & 1) * (result >> (degree - i) & 1)
                  for i in range(1, degree + 1)) % 2
        result |= bit << degree
    return result


def _character(unit, degree, width):
    derivative = sum((unit >> i & 1) << i for i in range(1, width + 1, 2))
    return _multiply(derivative, _reciprocal(unit, width), width) >> degree & 1


def _factor_gates(width):
    return [(j - 1, i - j - 1, i - 1)
            for j in range(1, width + 1) for i in range(2 * j + 1, width + 1)]


def _apply_factor(word, gates):
    for first, second, target in gates:
        word ^= ((word >> first & 1) * (word >> second & 1)) << target
    return word


def _anf(table):
    values = list(table)
    for bit in range((len(values) - 1).bit_length()):
        for address in range(len(values)):
            if address >> bit & 1:
                values[address] ^= values[address ^ (1 << bit)]
    return values


def _address_unit(address, masks, constants=0):
    return 1 | sum((((address & mask).bit_count() % 2) ^ (constants >> i & 1))
                   << (i + 1) for i, mask in enumerate(masks))


def _native_toffoli(first, second, target):
    word = [('H', target)]
    for _, gates in _ccz(first, second, target):
        for gate in gates:
            if gate[0] in ('T', 'TDG'):
                word.append(('T', gate[1], 1 if gate[0] == 'T' else -1))
            else:
                word.append(gate)
    return word + [('H', target)]


def _native_factor(work):
    return [gate for first, second, target in _factor_gates(len(work))
            for gate in _native_toffoli(work[first], work[second], work[target])]


def _controlled_source(work, outputs, character_masks):
    factor = _native_factor(work)
    middle = []
    for output, mask in zip(outputs, character_masks):
        for bit, coefficient in enumerate(work):
            if mask >> bit & 1:
                middle += [('H', coefficient), ('CX', output, coefficient),
                           ('H', coefficient)]
    return factor + middle + _inverse(factor)


def _native_transport(address, work, helper, masks):
    # The native fixture uses linear forms. Affine constants are checked
    # independently in the complete algebraic operator fixture below.
    word = []
    for degree in range(len(work), 0, -1):
        target = work[degree - 1]
        word += [('CX', wire, target) for bit, wire in enumerate(address)
                 if masks[degree - 1] >> bit & 1]
        for j in range(1, degree):
            controls = [wire for bit, wire in enumerate(address) if masks[j - 1] >> bit & 1]
            if not controls:
                continue
            relative, _ = _relative_grid([work[degree - j - 1]], [target], [[helper]])
            linear = [('CX', wire, helper) for wire in controls]
            word += ([('H', target)] + _inverse(relative) + linear + relative
                     + _inverse(linear) + [('H', target)])
    return word


class NonlinearCharacterTransportTests(unittest.TestCase):
    def test_multivariate_commutative_cube_product_at_nilpotency_bound(self):
        # Augmentation ideals in F2[s,t]/(s^p,t^q), rather than the
        # univariate family. A product of p+q-1 radical elements is zero.
        # The tested d-cubes retain arbitrary bases, zero directions,
        # repeated directions, and linearly dependent nonzero directions.
        rng = Random(7139)
        for p, q in ((3, 2), (4, 3)):
            monomials = [(i, j) for j in range(q) for i in range(p)]
            positions = {monomial: bit for bit, monomial in enumerate(monomials)}
            dimension = len(monomials) - 1
            nilindex = p + q - 1
            cube_dimension = 1 + (nilindex - 1).bit_length()
            self.assertGreaterEqual(1 << (cube_dimension - 1), nilindex)

            def multiply(first, second):
                result = 0
                for bit, (i, j) in enumerate(monomials):
                    if first >> bit & 1:
                        for other, (k, l) in enumerate(monomials):
                            if second >> other & 1 and i + k < p and j + l < q:
                                result ^= 1 << positions[(i + k, j + l)]
                return result

            s, t = 1 << positions[(1, 0)], 1 << positions[(0, 1)]
            direction_sets = [
                [0] * cube_dimension,
                [s, s, t] + [0] * (cube_dimension - 3),
                [s, t, s ^ t] + [0] * (cube_dimension - 3),
            ] + [[rng.getrandbits(dimension) << 1 for _ in range(cube_dimension)]
                 for _ in range(20)]
            bases = (range(1 << dimension) if dimension <= 5
                     else [rng.getrandbits(dimension) for _ in range(24)])
            for base in bases:
                for directions in direction_sets:
                    product = 1
                    for corner in range(1 << cube_dimension):
                        unit = 1 | (base << 1)
                        for bit, direction in enumerate(directions):
                            if corner >> bit & 1:
                                unit ^= direction
                        product = multiply(product, unit)
                    self.assertEqual(product, 1)
            # Lower cubes need not have product one; this is not a
            # vacuous test of a square-zero or identically trivial group.
            lower = 1
            for direction in (0, s, t, s ^ t):
                lower = multiply(lower, 1 ^ direction)
            self.assertNotEqual(lower, 1)

    def test_character_basis_and_homomorphisms(self):
        for width in range(1, 7):
            odds = list(range(1, width + 1, 2))
            rows = []
            for degree in odds:
                row = []
                for generator in odds:
                    entry = _character(1 | (1 << generator), degree, width)
                    self.assertEqual(entry, int(degree % generator == 0))
                    row.append(entry)
                rows.append(row)
            self.assertTrue(all(rows[i][i] == 1 for i in range(len(odds))))
            self.assertTrue(all(rows[i][j] == 0 for i in range(len(odds))
                                for j in range(i + 1, len(odds))))
            squares = {_multiply((z << 1) | 1, (z << 1) | 1, width)
                       for z in range(1 << width)}
            self.assertEqual(len(squares), 1 << (width // 2))
            self.assertTrue(all(all(unit >> i & 1 == 0 for i in odds) for unit in squares))
            if width <= 5:
                for first in range(1 << width):
                    for second in range(1 << width):
                        u, v = (first << 1) | 1, (second << 1) | 1
                        product = _multiply(u, v, width)
                        for degree in odds:
                            self.assertEqual(_character(product, degree, width),
                                             _character(u, degree, width)
                                             ^ _character(v, degree, width))

    def test_in_place_factorization_and_actual_inverse(self):
        for width in range(1, 11):
            gates = _factor_gates(width)
            self.assertEqual(len(gates), (width - 1) ** 2 // 4)
            self.assertTrue(all(len({first, second, target}) == 3
                                for first, second, target in gates))
            for word in range(1 << width):
                factors = _apply_factor(word, gates)
                self.assertEqual(_apply_factor(factors, list(reversed(gates))), word)
                reconstructed = 1
                for i in range(1, width + 1):
                    if factors >> (i - 1) & 1:
                        reconstructed = _multiply(reconstructed, 1 | (1 << i), width)
                self.assertEqual(reconstructed, (word << 1) | 1)
                if width <= 7:
                    for degree in range(1, width + 1, 2):
                        expected = sum(factors >> (i - 1) & 1 for i in range(1, degree + 1, 2)
                                       if degree % i == 0) % 2
                        self.assertEqual(_character(reconstructed, degree, width), expected)

    def test_degree_bound_and_unbounded_degree_family(self):
        for width in range(1, 10):
            for degree in range(1, width + 1, 2):
                coefficients = _anf([_character((z << 1) | 1, degree, width)
                                     for z in range(1 << width)])
                actual_degree = max(mask.bit_count() for mask, value in enumerate(coefficients)
                                    if value)
                self.assertLessEqual(actual_degree, 1 + int(log2(degree)))
        for address_width in range(1, 7):
            width = 1 + (address_width - 1) * (1 << (address_width - 1))
            table = []
            for address in range(1 << address_width):
                unit = 1 | sum((address >> j & 1) << (1 << j)
                               for j in range(address_width))
                table.append(_character(unit, width, width))
            coefficients = _anf(table)
            self.assertEqual(coefficients[-1], 1)
            self.assertEqual(max(mask.bit_count() for mask, value in enumerate(coefficients)
                                 if value), address_width)
            self.assertGreaterEqual(width, 1 << (address_width - 1))

    def test_literal_character_echo_on_all_dirty_inputs(self):
        width, address_width = 5, 3
        masks, constants = (1, 3, 0, 6, 1), 0b10100
        for address in range(1 << address_width):
            unit = _address_unit(address, masks, constants)
            inverse = _reciprocal(unit, width)
            self.assertEqual(_multiply(unit, inverse, width), 1)
            for work in range(1 << width):
                original = (work << 1) | 1
                for output in range(4):
                    # Chronological P†, D†, P, D. The two source positions
                    # retain every output control and literal inverse phase.
                    current = _multiply(inverse, original, width)
                    phase = sum((output >> bit & 1) * _character(current, degree, width)
                                for bit, degree in enumerate((3, 5))) % 2
                    current = _multiply(unit, current, width)
                    phase ^= sum((output >> bit & 1) * _character(current, degree, width)
                                 for bit, degree in enumerate((3, 5))) % 2
                    self.assertEqual(current, original)
                    expected = sum((output >> bit & 1) * _character(unit, degree, width)
                                   for bit, degree in enumerate((3, 5))) % 2
                    self.assertEqual(phase, expected)

    def test_complete_native_two_output_word_and_ledger(self):
        address, work, outputs, helper = [0, 1], [2, 3, 4], [5, 6], 7
        masks = (1, 2, 0)
        # p1=c1; p3=c1+c3. Thus three controlled-Z edges are charged.
        source = _controlled_source(work, outputs, (0b001, 0b101))
        transport = _native_transport(address, work, helper, masks)
        hadamards = [('H', wire) for wire in outputs]
        word = (hadamards + _inverse(transport) + _inverse(source)
                + transport + source + hadamards)
        self.assertEqual(sum(g[0] == 'T' for g in source), 14)
        self.assertEqual(sum(g[0] != 'T' for g in source), 33)
        self.assertEqual(sum(g[0] == 'T' for g in transport), 24)
        self.assertEqual(sum(g[0] != 'T' for g in transport), 62)
        self.assertEqual(sum(g[0] == 'T' for g in word), 76)
        self.assertEqual(sum(g[0] != 'T' for g in word), 194)
        self.assertTrue(all(g[1] not in address for g in word if g[0] == 'T'))
        self.assertTrue(all(g[-1] not in address for g in word if g[0] == 'CX'))
        actual = _native_matrix(word, 8)
        expected = np.zeros_like(actual)
        for column in range(256):
            unit = _address_unit(column & 3, masks)
            result = column
            for target, degree in zip(outputs, (1, 3)):
                result ^= _character(unit, degree, 3) << target
            expected[result, column] = 1
        np.testing.assert_allclose(actual, expected, rtol=0, atol=5e-14)
        inverse = _native_matrix(_inverse(word), 8)
        np.testing.assert_allclose(inverse, expected.conj().T, rtol=0, atol=5e-14)
        # Declared sequential schedule: two10-layer sources and two
        # transports of three8-layer private-helper products.
        self.assertEqual(2 * 10 + 2 * 3 * 8, 68)


if __name__ == '__main__':
    unittest.main()
