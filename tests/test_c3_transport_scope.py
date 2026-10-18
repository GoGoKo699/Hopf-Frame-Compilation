"""Exact evidence separating third-level rounding from degree transport.

The underlying subset-algebra permutation family is attributed to
He--Robitaille--Tan in the proof home. These checks retain its literal
native source/inverse and growing signed Clifford action. No three-band
source or fully returned XOR query is asserted.
"""
from __future__ import annotations

import itertools
import random
import unittest

from tests.exact_matrices import CO, CZ, native_columns
from tests.test_returned_diagonal_echo_degree import _ccz


def _support(value):
    """Algebra labels are nonempty subset masks; physical wire is label-1."""
    while value:
        bit = value & -value
        yield bit.bit_length()
        value ^= bit


def _multiply(left, right):
    result = 0
    for i in _support(left):
        for j in _support(right):
            if not i & j:
                result ^= 1 << ((i | j) - 1)
    return result


def _quadratic(value):
    result = 0
    for i, j in itertools.combinations(_support(value), 2):
        if not i & j:
            result ^= 1 << ((i | j) - 1)
    return result


def _k(value):
    return value ^ _quadratic(value)


def _f(value):
    result = 0
    for label in _support(value):
        basis = 1 << (label - 1)
        result ^= basis ^ _multiply(result, basis)
    return result


def _singleton_label(address, r):
    return sum(((address >> i) & 1) << ((1 << i) - 1) for i in range(r))


def _anf_degree(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, coefficient in enumerate(coefficients)
                if coefficient), default=-1)


def _toffolis(r):
    triples = [(i - 1, j - 1, (i | j) - 1)
               for i, j in itertools.combinations(range(1, 1 << r), 2)
               if not i & j]
    return sorted(triples, key=lambda triple: (-(triple[2] + 1).bit_count(), triple))


def _classical_word(triples, value):
    for left, right, target in triples:
        value ^= ((value >> left & 1) & (value >> right & 1)) << target
    return value


def _native_schedule(triples):
    schedule = []
    for left, right, target in triples:
        schedule += [('C', [('H', target)])]
        schedule += _ccz(left, right, target)
        schedule += [('C', [('H', target)])]
    return schedule


def _inverse(schedule):
    names = {'T': 'TDG', 'TDG': 'T', 'H': 'H', 'CX': 'CX'}
    return [(kind, [(names[name], *wires) for name, *wires in reversed(gates)])
            for kind, gates in reversed(schedule)]


class ThirdLevelTransportScope(unittest.TestCase):
    def test_nilalgebra_homomorphism_and_inverse_beyond_smallest_example(self):
        for r in (2, 3):
            limit = 1 << ((1 << r) - 1)
            for left, right in itertools.product(range(limit), repeat=2):
                product = _multiply(left, right)
                self.assertEqual(_quadratic(product), 0)
                self.assertEqual(_quadratic(left ^ right),
                                 _quadratic(left) ^ _quadratic(right) ^ product)
                self.assertEqual(_k(left ^ right ^ product), _k(left) ^ _k(right))
            for value in range(limit):
                self.assertEqual(_multiply(value, value), 0)
                self.assertEqual(_f(_k(value)), value)
                self.assertEqual(_k(_f(value)), value)
        generator = random.Random(20261018)
        # Widths 15,31,63 test the same identities beyond the seven-qubit
        # dense native fixture, without treating samples as a proof.
        for r in (4, 5, 6):
            width = (1 << r) - 1
            triples = _toffolis(r)
            self.assertEqual(len(triples), (3 ** r - 2 ** (r + 1) + 1) // 2)
            for _ in range(160):
                a, b, c = [generator.getrandbits(width) for _ in range(3)]
                self.assertEqual(_multiply(_multiply(a, b), c), _multiply(a, _multiply(b, c)))
                self.assertEqual(_quadratic(_multiply(a, b)), 0)
                self.assertEqual(_k(a ^ b ^ _multiply(a, b)), _k(a) ^ _k(b))
                self.assertEqual(_f(_k(a)), a)
                self.assertEqual(_k(_f(a)), a)
                self.assertEqual(_classical_word(triples, a), _k(a))
                self.assertEqual(_classical_word(list(reversed(triples)), a), _f(a))

    def test_images_of_pauli_generators_are_literal_affine_or_quadratic(self):
        r = 3
        width = (1 << r) - 1
        limit = 1 << width
        for basis in (1 << wire for wire in range(width)):
            translation = _f(basis)
            for value in range(limit):
                self.assertEqual(_f(_k(value) ^ basis),
                                 value ^ translation ^ _multiply(value, translation))
        for wire in range(width):
            # F Z_wire F† has the literal phase (-1)^(k(y)_wire).
            self.assertLessEqual(_anf_degree([_k(value) >> wire & 1 for value in range(limit)]), 2)
        singleton = _singleton_label(7, r)
        self.assertNotEqual(_f(singleton), _k(singleton))
        self.assertEqual(_anf_degree([_f(value) >> (width - 1) & 1
                                      for value in range(limit)]), 3)

    def test_transported_symplectic_and_literal_probe_sign_degrees_grow(self):
        for r in range(2, 7):
            width = (1 << r) - 1
            top = width - 1
            translated = [_f(_singleton_label(address, r)) for address in range(1 << r)]
            for subset in range(1, 1 << r):
                values = [value >> (subset - 1) & 1 for value in translated]
                self.assertEqual(_anf_degree(values), subset.bit_count())
            # Off-diagonal B_(top,{1}) is the coefficient of e_top
            # in e_{1}a(x), hence the product of the other address bits.
            entries = [_multiply(1, value) >> top & 1 for value in translated]
            signs = [value >> top & 1 for value in translated]
            self.assertEqual(_anf_degree(entries), r - 1)
            self.assertEqual(_anf_degree(signs), r)
            for value in translated:
                self.assertEqual(value ^ _multiply(value, value), value)
            if r == 3:
                # Verify the signed Z_top formula on every work column,
                # rather than inferring it solely from its coefficient.
                for address, value in enumerate(translated):
                    probe = sum(((((1 << wire) ^ _multiply(1 << wire, value)) >> top) & 1)
                                << wire for wire in range(width))
                    for work in range(1 << width):
                        image = _f(_k(work) ^ _singleton_label(address, r))
                        self.assertEqual(image >> top & 1,
                                         (value >> top & 1) ^ ((probe & work).bit_count() & 1))

    def test_native_source_actual_inverse_all_columns_and_resource_ledger(self):
        r = 3
        width = (1 << r) - 1
        triples = _toffolis(r)
        count = (3 ** r - 2 ** (r + 1) + 1) // 2
        schedule = _native_schedule(triples)
        inverse = _inverse(schedule)
        for actual, target in ((schedule, _k), (inverse, _f)):
            word = [gate for _, gates in actual for gate in gates]
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 7 * count)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 12 * count)
            self.assertEqual(sum(kind == 'T' for kind, _ in actual), 5 * count)
            for kind, gates in actual:
                self.assertTrue(all(all(0 <= wire < width for wire in gate[1:]) for gate in gates))
                if kind == 'T':
                    self.assertEqual(len({gate[1] for gate in gates}), len(gates))
            matrix = native_columns(width, word)
            for column in range(1 << width):
                row = target(column)
                self.assertEqual([matrix[index][column] for index in range(1 << width)],
                                 [CO if index == row else CZ for index in range(1 << width)])


if __name__ == '__main__':
    unittest.main()
