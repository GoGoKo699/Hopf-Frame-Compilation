"""Exact fixtures for four-call Gram recovery and its two sector obligations.

All comparisons retain native scalar phases and every physical column.
These bounded identities do not construct a higher-band lookup source.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse


HALF = CQ2(Q2(Fraction(1, 2)))


def _subtract(left, right):
    return add(left, scale(right, -CO))


def _clean_z(width):
    half = 1 << (width - 1)
    return [[(CO if column < half else -CO) if row == column else CZ
             for column in range(2 * half)] for row in range(2 * half)]


def _embedding(width):
    size = 1 << (width - 1)
    return [[CO if row == column else CZ for column in range(size)]
            for row in range(2 * size)]


def _gram(unitary):
    size = len(unitary) // 2
    accepted = [row[:size] for row in unitary[:size]]
    return multiply(adjoint(accepted), accepted)


def _readout(source, signal):
    # F=(Z U† Z U)^2 in chronological order, with actual adjoints.
    z = [('S', signal), ('S', signal)]
    return (source + z + _inverse(source) + z) * 2


def _partial_gram_source(private_fiber=False):
    # Equal c,d sectors use H_a; unequal sectors use X_a. An optional
    # private H fiber demonstrates that extraction must return it.
    signal = 3 if private_fiber else 2
    rotation = ([('S', signal)] * 3
                + [('H', signal), ('T', signal),
                   ('H', signal), ('S', signal)])
    controlled_h = (_inverse(rotation)
                    + [('H', signal), ('CX', 0, signal), ('H', signal)]
                    + rotation)
    return (([('H', 2)] if private_fiber else [])
            + [('CX', 1, 0), ('H', signal)]
            + controlled_h + [('CX', 0, signal), ('CX', 1, 0)])


class GramFourCallReadoutTests(unittest.TestCase):
    def test_literal_negative_reflection_returns_all_private_fibers(self):
        source = _partial_gram_source(private_fiber=True)
        unitary = native_columns(4, source)
        gram = _gram(unitary)
        for row in range(8):
            for column in range(8):
                equal = (column & 1) == ((column >> 1) & 1)
                self.assertEqual(gram[row][column],
                                 HALF if equal and row == column else CZ)

        readout = _readout(source, 3)
        actual = native_columns(4, readout)
        self.assertEqual(native_columns(4, _inverse(readout)), adjoint(actual))
        for column in range(16):
            equal = (column & 1) == ((column >> 1) & 1)
            for row in range(16):
                self.assertEqual(actual[row][column],
                                 (-CO if equal else CO) if row == column else CZ)

    def test_noncommuting_gram_has_exact_complete_error_square(self):
        # Logical d=0 and clean a=1. The input H_d makes G fail to
        # commute with the chosen Pi=|0><0|; it cannot be scalarized.
        source = [('H', 0), ('H', 1), ('T', 1), ('H', 1), ('CX', 0, 1)]
        unitary = native_columns(2, source)
        gram = _gram(unitary)
        pi = [[CO, CZ], [CZ, CZ]]
        self.assertNotEqual(multiply(gram, pi), multiply(pi, gram))
        complement = _subtract(eye(2), pi)
        delta = _subtract(gram, scale(pi, HALF))
        f = add(_subtract(scale(multiply(gram, gram), 8), scale(gram, 8)), eye(2))
        actual = native_columns(2, _readout(source, 1))
        self.assertEqual([row[:2] for row in actual[:2]], f)

        j = _embedding(2)
        s0 = _subtract(eye(2), scale(pi, 2))
        difference = _subtract(multiply(actual, j), multiply(j, s0))
        error_square = multiply(adjoint(difference), difference)
        first = multiply(multiply(pi, multiply(delta, delta)), pi)
        second = multiply(multiply(complement,
                                   multiply(gram, _subtract(eye(2), gram))),
                          complement)
        self.assertEqual(error_square, scale(add(first, second), 16))

    def test_approximate_readout_needs_the_actual_reversed_word(self):
        source = [('H', 0), ('T', 0), ('H', 0)]
        readout = _readout(source, 0)
        actual = native_columns(1, readout)
        self.assertEqual(actual, [[CZ, -CI], [-CI, CZ]])
        self.assertNotEqual(actual, adjoint(actual))
        inverse = native_columns(1, _inverse(readout))
        self.assertEqual(inverse, adjoint(actual))
        self.assertEqual(multiply(inverse, actual), eye(2))
        self.assertEqual(multiply(actual, actual), scale(eye(2), -CO))

    def test_spectral_closure_and_sector_matching_are_distinct(self):
        identity = eye(8)
        zero = zeros(8)
        z = _clean_z(3)
        j = _embedding(3)
        pi = [[CO if row == column and column in (0, 3) else CZ
               for column in range(4)] for row in range(4)]
        for rotate_input in (False, True):
            source = (([('H', 0)] if rotate_input else [])
                      + _partial_gram_source())
            ua = native_columns(3, source)
            # UB=I is enough to test the exact diagnostic. No band or
            # asymptotic synthesis claim is inferred from this split.
            a = multiply(multiply(ua, z), adjoint(ua))
            c = multiply(a, z)
            c_squared = multiply(c, c)
            self.assertEqual(multiply(c_squared, c_squared), identity)
            self.assertEqual(multiply(add(c, identity), add(c_squared, identity)), zero)
            pa = scale(add(identity, a), HALF)
            actual_sector = scale(multiply(pa, _subtract(identity, c_squared)), HALF)
            desired_sector = multiply(multiply(multiply(multiply(ua, j), pi),
                                               adjoint(j)), adjoint(ua))
            if rotate_input:
                self.assertNotEqual(actual_sector, desired_sector)
            else:
                self.assertEqual(actual_sector, desired_sector)

        # C=I has fourth-order closure but retains a Gram-one sector.
        self.assertEqual(multiply(identity, identity), identity)
        self.assertNotEqual(multiply(add(identity, identity),
                                     add(identity, identity)), zero)


if __name__ == '__main__':
    unittest.main()
