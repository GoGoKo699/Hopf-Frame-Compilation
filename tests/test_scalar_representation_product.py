"""Exact bounded evidence for the four-representation scope counterexample.

The arbitrary-degree conclusion follows from the analytic construction.
These checks retain every basis column and literal sign in small banks.
"""
from __future__ import annotations

import unittest

from tests.test_read_only_three_layer_rigidity import _degree
from tests.test_returned_diagonal_echo_degree import _ccz, _column, _inverse


def _gamma(column, label, table):
    sign = sum(table[label ^ other] * ((column >> other) & 1)
               for other in range(label)) & 1
    return column ^ (1 << label), sign


def _translation(column, shift, table):
    size = len(table)
    sign = sum(table[left ^ right] * ((column >> left) & 1)
               * ((column >> right) & 1)
               for left in range(size) for right in range(left + 1, size)
               if (left ^ shift) > (right ^ shift)) & 1
    result = sum(((column >> label) & 1) << (label ^ shift)
                 for label in range(size))
    return result, sign


def _native_source(width, table):
    """Physical chronological order; CCZ and Fredkin are exact macros."""
    gates = []
    for coordinate in range(width):
        shift = 1 << coordinate
        for left in range(len(table)):
            for right in range(left + 1, len(table)):
                if table[left ^ right] and (left ^ shift) > (right ^ shift):
                    gates.append(('CCZ', coordinate, width + left,
                                  width + right))
        for left in range(len(table)):
            if not left & shift:
                gates.append(('F', coordinate, width + left,
                              width + (left ^ shift)))
    return gates


def _apply_source(gates, column):
    sign = 0
    for kind, control, left, right in gates:
        if kind == 'CCZ':
            sign ^= ((column >> control) & (column >> left)
                     & (column >> right) & 1)
        elif column >> control & 1 and ((column >> left) ^ (column >> right)) & 1:
            column ^= (1 << left) | (1 << right)
    return column, sign


def _tables(width):
    size = 1 << width
    if width <= 2:
        return [[0] + [(mask >> index) & 1 for index in range(size - 1)]
                for mask in range(1 << (size - 1))]
    return [[0] * size,
            [int(x == size - 1) for x in range(size)],
            [x & 1 for x in range(size)],
            [int((x & 3) == 3) for x in range(size)],
            [int(x != 0) for x in range(size)],
            [0, 1, 1, 0, 0, 1, 0, 1]]


class ScalarRepresentationProductTests(unittest.TestCase):
    def test_literal_paulis_translation_and_honest_group_law(self):
        for width in range(1, 4):
            size = 1 << width
            for table in _tables(width):
                actions = [[_translation(column, x, table)
                            for column in range(1 << size)] for x in range(size)]
                for column in range(1 << size):
                    for left in range(size):
                        for right in range(size):
                            a, sa = _gamma(column, right, table)
                            a, sb = _gamma(a, left, table)
                            b, sc = _gamma(column, left, table)
                            b, sd = _gamma(b, right, table)
                            self.assertEqual((a, sa ^ sb),
                                             (b, sc ^ sd ^ table[left ^ right]))
                    for x in range(size):
                        moved, sx = actions[x][column]
                        for label in range(size):
                            # C_x gamma_label = gamma_(label+x) C_x.
                            a, sa = _gamma(column, label, table)
                            a, sb = actions[x][a]
                            b, sc = _gamma(moved, label ^ x, table)
                            self.assertEqual((a, sa ^ sb), (b, sx ^ sc))
                        for y in range(size):
                            a, sa = actions[y][column]
                            b, sb = actions[x][a]
                            self.assertEqual((b, sa ^ sb), actions[x ^ y][column])

    def test_native_source_canonical_images_and_all_column_scalar_return(self):
        for width in range(1, 4):
            size = 1 << width
            zmask = size - 1
            for table in _tables(width):
                gates = _native_source(width, table)
                actual_inverse = list(reversed(gates))
                source = [_apply_source(gates, column)
                          for column in range(1 << (width + size))]
                for column, (image, phase) in enumerate(source):
                    bank, sign = _translation(column >> width,
                                              column & zmask, table)
                    self.assertEqual((image, phase),
                                     ((bank << width) | (column & zmask), sign))
                    back, inverse_phase = _apply_source(actual_inverse, image)
                    self.assertEqual((back, phase ^ inverse_phase), (column, 0))
                    for x in range(size):
                        # U X_z^x U^dagger, using the verified actual inverse.
                        first, s1 = source[column]
                        second, s2 = source[first ^ x]
                        bank, expected_phase = _translation(column >> width,
                                                             x, table)
                        expected = (bank << width) | ((column & zmask) ^ x)
                        self.assertEqual((second, s1 ^ s2),
                                         (expected, expected_phase))
                        # Physical order is B, A, B, A. B has gamma_0 on
                        # both ends; every U/source inverse call is retained.
                        value, scalar = column, 0
                        for is_b in (True, False, True, False):
                            if is_b:
                                value ^= 1 << width
                            value, sign = source[value]
                            scalar ^= sign
                            value ^= x
                            value, sign = source[value]
                            scalar ^= sign
                            if is_b:
                                value ^= 1 << width
                        self.assertEqual((value, scalar), (column, table[x]))

    def test_source_counts_literal_ccz_and_unbounded_degree_family(self):
        schedule = _ccz(0, 1, 2)
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual((sum(map(len, layers)), len(layers)), (7, 5))
        self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'C'),
                         10)
        for column in range(8):
            self.assertEqual(_column(schedule, column),
                             (column, 4 * int(column == 7)))
            self.assertEqual(_column(schedule + _inverse(schedule), column),
                             (column, 0))
        for width in range(1, 6):
            size = 1 << width
            monomial = [int(x == size - 1) for x in range(size)]
            self.assertEqual(_degree(monomial), width)
            tables = _tables(width) if width < 4 else [monomial]
            for table in tables:
                gates = _native_source(width, table)
                swaps = sum(gate[0] == 'F' for gate in gates)
                cc_z = sum(gate[0] == 'CCZ' for gate in gates)
                self.assertEqual(swaps, width * size // 2)
                self.assertEqual(cc_z, size * sum(table) // 2)
                t_source = 7 * (swaps + cc_z)
                g_source = 14 * swaps + 10 * cc_z
                self.assertEqual(8 * t_source, 56 * (swaps + cc_z))
                self.assertEqual(8 * g_source + 4 * width + 16,
                                 112 * swaps + 80 * cc_z + 4 * width + 16)
                self.assertEqual(8 * 5 * (swaps + cc_z),
                                 40 * (swaps + cc_z))


if __name__ == '__main__':
    unittest.main()
