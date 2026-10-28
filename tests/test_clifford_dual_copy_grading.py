"""Exact fixtures for Clifford-dual grading and symplectic relation screens.

The all-size theorem lives in CLIFFORD_DUAL_COPY_GRADING.md.  These tests
check finite modules, literal native duals and scalar phases; they do not
infer a uniform group-relation length or an unrestricted depth bound.
"""
from collections import deque
from fractions import Fraction
from itertools import combinations
from math import comb
import unittest

from tests.exact_matrices import adjoint, eye, multiply, native_columns, scale
from tests.test_joint_lookup_pauli_source import _inverse


def _translate(table, shift):
    answer = 0
    while table:
        bit = table & -table
        answer ^= 1 << ((bit.bit_length() - 1) ^ shift)
        table ^= bit
    return answer


def _difference(table, directions):
    for direction in directions:
        table ^= _translate(table, direction)
    return table


def _anf_degree(table, p):
    coefficients = [(table >> x) & 1 for x in range(1 << p)]
    for bit in range(p):
        for mask in range(1 << p):
            if mask & (1 << bit):
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, value in enumerate(coefficients)
                if value), default=-1)


def _rank(rows):
    basis = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit in basis:
                row ^= basis[bit]
            else:
                basis[bit] = row
                break
    return len(basis)


def _identity(width):
    return tuple(1 << bit for bit in range(width))


def _compose(left, right):
    """Binary row-matrix composition, independently of native simulation."""
    output = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length() - 1]
            row ^= bit
        output.append(value)
    return tuple(output)


def _group_binary_abelianization(generators):
    """Enumerate a bounded group; edge relations determine its F2 quotient."""
    identity = _identity(len(generators[0]))
    paths = {identity: 0}
    queue = deque([identity])
    relations = []
    while queue:
        current = queue.popleft()
        for index, generator in enumerate(generators):
            product = _compose(current, generator)
            parity = paths[current] ^ (1 << index)
            if product not in paths:
                paths[product] = parity
                queue.append(product)
            else:
                relations.append(parity ^ paths[product])
    return len(paths), len(generators) - _rank(relations)


class CliffordDualCopyGradingTests(unittest.TestCase):
    def test_regular_annihilators_equal_low_degree_truth_tables(self):
        # Every vector, every ideal power, for p <= 3.  The independent
        # truth-table ANF calculation also fixes the regular-module basis.
        for p in range(1, 4):
            for d in range(p + 1):
                annihilated = 0
                directions = list(combinations([1 << bit for bit in range(p)],
                                               d + 1))
                for table in range(1 << (1 << p)):
                    actual = all(_difference(table, cube) == 0
                                 for cube in directions)
                    self.assertEqual(actual, _anf_degree(table, p) <= d)
                    annihilated += actual
                dimension = sum(comb(p, j) for j in range(d + 1))
                self.assertEqual(annihilated, 1 << dimension)
                if d == 1:
                    self.assertEqual(dimension, p + 1)

    def test_field_pauli_grade_map_is_equivariant_and_has_rank_p_plus_one(self):
        # Q_x X_u0 Q_x^dagger = X_u0 X_v^x.  These first p companion
        # columns undergo no field reduction since p <= k=2^(p-1).
        for p in range(1, 6):
            k, m = 1 << (p - 1), 1 << p

            def orbit_label(x):
                return 1 | (x << k)

            def grade(z):
                return sum(((orbit_label(x) & z).bit_count() & 1) << x
                           for x in range(m))

            generators = [grade(1 << bit) for bit in range(2*k)]
            self.assertEqual(_rank(generators), p + 1)
            for z in [0, 1, (1 << (2*k)) - 1,
                      *(1 << (k + bit) for bit in range(p))]:
                table = grade(z)
                self.assertLessEqual(_anf_degree(table, p), 1)
                for a, b in combinations(range(p), 2):
                    self.assertEqual(_difference(table, (1 << a, 1 << b)), 0)
                for y in range(m):
                    # Moving x to x+y changes only the constant Pauli
                    # pairing, by <X_v^y,Z^z>.
                    shifted_z = z ^ (((y & (z >> k)).bit_count() & 1))
                    self.assertEqual(_translate(table, y), grade(shifted_z))

    def test_literal_one_band_dual_and_actual_inverse(self):
        # Chronological word for C=B(V_u^dagger V_v^dagger),
        # V=S H T H S^dagger.  Its two T-dagger gates form one band.
        word = [('SDG', 0), ('SDG', 1), ('H', 0), ('H', 1),
                ('TDG', 0), ('TDG', 1), ('H', 0), ('H', 1),
                ('S', 0), ('S', 1), ('CX', 0, 2), ('CX', 1, 2)]
        copy = native_columns(3, word)
        reverse = native_columns(3, _inverse(word))
        self.assertEqual(reverse, adjoint(copy))
        self.assertEqual(multiply(copy, reverse), eye(8))
        self.assertEqual([gate for gate in word if gate[0] in ('T', 'TDG')],
                         [('TDG', 0), ('TDG', 1)])
        output_z = native_columns(3, [('Z', 2)])
        dual = multiply(multiply(reverse, output_z), copy)
        expected = native_columns(3, [('H', 0), ('H', 1), ('Z', 2)])
        self.assertEqual(dual, expected)
        query = native_columns(3, [('CX', 0, 1)])
        translated = multiply(multiply(query, dual), adjoint(query))
        self.assertNotEqual(multiply(dual, translated), multiply(translated, dual))
        pair = multiply(dual, translated)
        self.assertEqual(multiply(multiply(pair, pair), pair), eye(8))
        # The literal six-generator scalar relation is stronger than
        # the generic 24-generator field certificate in this example.
        self.assertEqual(Fraction(1, 6), Fraction(2, 2*6))

    def test_noncommuting_s3_pair_and_eight_address_odd_parity(self):
        # X0,X1,Z0,Z1 row coordinates; H on both wires and CNOT(0;1).
        h = (4, 8, 1, 2)
        q = (1, 3, 12, 8)
        j1 = _compose(_compose(q, h), q)
        identity = _identity(4)
        pair = _compose(h, j1)
        self.assertEqual(_compose(pair, _compose(pair, pair)), identity)
        self.assertNotEqual(pair, identity)
        self.assertNotEqual(_compose(h, j1), _compose(j1, h))
        self.assertEqual(_group_binary_abelianization((h, j1)), (6, 1))
        for p in range(3, 7):
            for z, a, b in ((1, 2, 4), (1, 3, 5)):
                parity = _difference((1 << 0) | (1 << z), (a, b))
                self.assertEqual(parity.bit_count(), 8)
                self.assertEqual(parity,
                                 sum(1 << x for x in (0, z, a, a ^ z,
                                                      b, b ^ z, a ^ b,
                                                      a ^ b ^ z)))
                for target in range(1 << p):
                    translated = _translate(parity, target)
                    self.assertEqual((translated >> target) & 1, 1)
            self.assertEqual(_difference(3, (1, 2)), 0)
        r = 3
        length = 4 * 2 * r
        self.assertEqual(Fraction(2, 2*length), Fraction(1, 8*r))

    def test_translated_pauli_cube_retains_a_negative_scalar(self):
        # An additive degree-one Clifford query; its four translated
        # Pauli values multiply to -I, not the silently assumed +I.
        pauli = native_columns(2, [('H', 0), ('Z', 0), ('H', 0),
                                   ('H', 1), ('Z', 1), ('H', 1)])
        scalar, parity = eye(4), 0
        for address in range(4):
            query = native_columns(2, [('H', bit) for bit in range(2)
                                       if (address >> bit) & 1])
            shifted = multiply(multiply(query, pauli), adjoint(query))
            scalar = multiply(scalar, shifted)
            parity ^= 1 << address
        self.assertEqual(scalar, scale(eye(4), -1))
        self.assertEqual(parity, _difference(1, (1, 2)))
        self.assertEqual(parity, 15)

    def test_odd_weight_readouts_span_all_sign_characters(self):
        for p in range(1, 4):
            m = 1 << p
            for f in range(1 << m):
                if not f.bit_count() & 1:
                    continue
                translates = [_translate(f, x) for x in range(m)]
                self.assertEqual(_rank(translates), m)
                square = 0
                for x in range(m):
                    if (f >> x) & 1:
                        square ^= translates[x]
                self.assertEqual(square, 1)


if __name__ == '__main__':
    unittest.main()
