"""Exact Pauli-orbit copy screen and the dirty address-mask gauge.

The obstruction applies to copies whose output-Z duals are Clifford.
The one-address native endpoint uses noncommuting carrier involutions.
All equalities retain scalar phases and arbitrary input/helper values.
"""
from collections import Counter
import unittest

from tests.exact_matrices import (
    CI, CO, adjoint, eye, multiply, native_columns, zeros,
)
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_field_shear_intertwiner import (
    _native_x, _permutation_matrix, _shear,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_t_depth import _word


def _scale(matrix, scalar):
    return [[entry * scalar for entry in row] for row in matrix]


def _pauli(width, label):
    """Hermitian Pauli: low bits are X, high bits are Z."""
    mask = (1 << width) - 1
    x, z = label & mask, label >> width
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    result = zeros(1 << width)
    for column in range(1 << width):
        sign = -CO if (z & column).bit_count() % 2 else CO
        result[column ^ x][column] = phase * sign
    return result


def _conjugate(unitary, operator):
    return multiply(multiply(unitary, operator), adjoint(unitary))


def _pauli_orbit(operator, offset, first, second):
    """Affine labels, with deliberately unrelated scalar query phases."""
    width = len(operator).bit_length() - 1
    orbit = []
    for address in range(4):
        label = offset
        if address & 1:
            label ^= first
        if address & 2:
            label ^= second
        phase = (CO, CI, -CO, -CI)[(address + (address >> 1)) % 4]
        query = _scale(_pauli(width, label), phase)
        orbit.append(_conjugate(adjoint(query), operator))
    return orbit


def _four_point_word(orbit):
    return multiply(multiply(multiply(orbit[1], orbit[0]), orbit[2]),
                    orbit[3])


def _one_address_copy(carrier, out0, out1):
    """Two commuting independent T axes, one native physical band."""
    canonical = [('H', out0), ('CX', out0, carrier),
                 ('CX', carrier, out1), ('H', out0)]
    diagonalizer = [('CX', carrier, out0), ('CX', out1, out0),
                   ('H', carrier), ('S', carrier),
                   ('H', out0), ('H', out1)]
    return (canonical + _inverse(diagonalizer)
            + [('T', carrier), ('TDG', out0)] + diagonalizer)


class PauliOrbitCopyTests(unittest.TestCase):
    def test_four_point_word_for_hermitian_cliffords_and_affine_paulis(self):
        wall = native_columns(
            2, [('H', 0), ('S', 1), ('CX', 0, 1), ('H', 1)])
        base_words = (
            [('H', 0)],
            [('CX', 0, 1)],
            [('H', 1), ('CX', 0, 1), ('H', 1)],
            [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)],
            [('H', 0), ('H', 1)],
            [('S', 0), ('S', 0)],
        )
        affine_maps = ((0, 1, 4), (13, 2, 4), (3, 9, 6), (5, 2, 2))
        for word in base_words:
            operator = _conjugate(wall, native_columns(2, word))
            self.assertEqual(operator, adjoint(operator))
            self.assertEqual(multiply(operator, operator), eye(4))
            for labels in affine_maps:
                with self.subTest(word=word, labels=labels):
                    orbit = _pauli_orbit(operator, *labels)
                    product = _four_point_word(orbit)
                    scalar = product[0][0]
                    self.assertIn(scalar, (CO, CI, -CO, -CI))
                    self.assertEqual(product, _scale(eye(4), scalar))
                    for member in orbit:
                        self.assertEqual(member, adjoint(member))
                        self.assertEqual(multiply(member, member), eye(4))

        # A central word can be -I: dropping its literal phase is invalid.
        orbit = _pauli_orbit(native_columns(2, [('H', 0)]), 0, 1, 4)
        central = _four_point_word(orbit)
        self.assertEqual(central, _scale(eye(4), -CO))
        # p=2 and j=3 demand one minus sign on this four-point plane.
        requested = list(orbit)
        requested[3] = _scale(requested[3], -CO)
        target = _four_point_word(requested)
        self.assertEqual(target, _scale(central, -CO))
        difference = [[left - right for left, right in zip(row, other)]
                      for row, other in zip(central, target)]
        self.assertEqual(multiply(adjoint(difference), difference),
                         _scale(eye(4), CO * 4))

    def test_native_noncommuting_one_address_copy_and_actual_inverse(self):
        word = _one_address_copy(0, 1, 2)
        actual = native_columns(3, word)
        reverse = native_columns(3, _inverse(word))
        self.assertEqual(reverse, adjoint(actual))
        self.assertNotEqual(reverse, actual)
        seed = native_columns(3, [('H', 0)])
        query = native_columns(3, [('S', 0), ('S', 0)])
        orbit = [seed, _conjugate(query, seed)]
        self.assertEqual(multiply(orbit[0], orbit[1]),
                         _scale(multiply(orbit[1], orbit[0]), -CO))
        for address, operator in enumerate(orbit):
            output = native_columns(3, _native_x(1 + address))
            self.assertEqual(_conjugate(actual, operator),
                             multiply(operator, output))
        square = native_columns(
            3, [('H', 1), ('CX', 1, 2), ('H', 1)])
        self.assertEqual(multiply(actual, actual), square)
        for native in (word, _inverse(word)):
            self.assertEqual(len(native), 20)
            counts = Counter(gate[0] for gate in native)
            self.assertEqual(counts['T'] + counts['TDG'], 2)
            bands = _bands(native)
            self.assertEqual(len(bands), 1)
            self.assertEqual(len({gate[1] for gate in bands[0]}), 2)
            self.assertEqual({wire for gate in native for wire in gate[1:]},
                             {0, 1, 2})

    def test_one_address_indicator_full_native_echo(self):
        # Address x=0, carrier u=1, outputs Y0=2 and Y1=3.
        query = [('H', 1), ('CX', 0, 1), ('H', 1)]
        source = _inverse(query) + [('H', 1)] + query
        copy = _one_address_copy(1, 2, 3)
        word = _inverse(source) + _inverse(copy) + source + copy
        expected = _permutation_matrix([
            column ^ (1 << (2 + (column & 1))) for column in range(16)])
        self.assertEqual(native_columns(4, word), expected)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(expected))
        for column in range(16):
            row = column ^ (1 << (2 + (column & 1)))
            self.assertEqual(row & 3, column & 3)
        for native in (word, _inverse(word)):
            counts = Counter(gate[0] for gate in native)
            self.assertEqual(counts['T'] + counts['TDG'], 4)
            self.assertLessEqual(len(_bands(native)), 2)
            for band in _bands(native):
                self.assertEqual(len({gate[1] for gate in band}), len(band))

    def test_elementary_dirty_pauli_mask_gauge_all_columns(self):
        # Native p=1: x,t,u,v,d, including both arbitrary dirty-helper values.
        gauge = _word(_compact_toffoli(1, 2, 3, 4))
        query = _word(_compact_toffoli(0, 2, 3, 4))
        masked_query = query + [('CX', 0, 1)]
        word = _inverse(gauge) + masked_query + gauge
        expected = native_columns(5, [('CX', 0, 1)])
        self.assertEqual(native_columns(5, word), expected)
        self.assertEqual(native_columns(5, _inverse(word)), adjoint(expected))
        for column in range(32):
            row = column ^ (2 if column & 1 else 0)
            self.assertEqual(row & 28, column & 28)
        # Independent complete p=2 field action: all x,t and all bank inputs.
        for address in range(4):
            for mask in range(4):
                for bank in range(16):
                    returned = _shear(
                        2, mask ^ address,
                        _shear(2, address, _shear(2, mask, bank)))
                    self.assertEqual(returned, bank)


if __name__ == '__main__':
    unittest.main()
