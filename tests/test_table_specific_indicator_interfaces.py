"""A table-specific source can mix a syndrome's arbitrary bank fiber.

Every column and the actual inverse of the paid extraction are checked.
This fixture distinguishes the smaller probe contract from a universal
indicator source; it is not an asymptotic synthesis result.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _counts, _cz, _x


class TableSpecificIndicatorInterfacesTests(unittest.TestCase):
    def test_nonlinear_row_query_with_address_and_kernel_bank_gauge(self):
        # x=(0,1), A=(2,3,4,5).  The first two C coordinates are
        # M A with M rows 1110 and 0101.  The other coordinates are
        # old A0,A3, and may be changed arbitrarily by the gauge.
        coordinate = [('CX', 2, 4), ('CX', 4, 2), ('CX', 2, 4)]
        coordinate += [('CX', 4, 2), ('CX', 3, 2), ('CX', 5, 3)]
        coordinate_matrix = native_columns(6, coordinate)
        for column in range(64):
            bits = [(column >> wire) & 1 for wire in (2, 3, 4, 5)]
            transformed = (bits[0] ^ bits[1] ^ bits[2],
                           bits[1] ^ bits[3], bits[0], bits[3])
            expected = (column & 3) | sum(bit << (2 + index)
                                          for index, bit in enumerate(transformed))
            for row in range(64):
                self.assertEqual(coordinate_matrix[row][column],
                                 CO if row == expected else CZ)

        # M e_x = (1+x0*x1, x0): a nonlinear literal row lookup.
        row_query = _x(2) + [('H', 2)] + _ccz((0, 1, 2)) + [('H', 2)]
        row_query += [('CX', 0, 3)]
        # G preserves S=(2,3), but genuinely mixes address and kernel
        # coordinates, with noncommuting bases and a native T phase.
        gauge = [('H', 0), ('CX', 0, 4), ('T', 4),
                 ('CX', 5, 1), ('H', 5)]
        source = coordinate + gauge + row_query + _inverse(coordinate)
        self.assertEqual(_counts(coordinate), (0, 6))
        self.assertEqual(_counts(row_query), (7, 17))
        for orientation in (source, _inverse(source)):
            self.assertEqual(_counts(orientation), (8, 33))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))
        source_matrix = native_columns(6, source)
        self.assertEqual(native_columns(6, _inverse(source)), adjoint(source_matrix))
        self.assertEqual(sum(source_matrix[row][0] != CZ for row in range(64)), 4)

        # V=(6,7) is arbitrary.  The Clifford probe reads only M A.
        probe = (_cz(2, 6) + _cz(3, 6) + _cz(4, 6)
                 + _cz(3, 7) + _cz(5, 7))
        output_h = [('H', 6), ('H', 7)]
        wrapper = output_h + probe + _inverse(source) + probe + source + output_h
        self.assertEqual(_counts(probe), (0, 15))
        for orientation in (wrapper, _inverse(wrapper)):
            self.assertEqual(_counts(orientation), (16, 100))
            matrix = native_columns(8, orientation)
            for column in range(256):
                x0, x1 = column & 1, (column >> 1) & 1
                expected = column ^ ((1 ^ (x0 & x1)) << 6) ^ (x0 << 7)
                for row in range(256):
                    self.assertEqual(matrix[row][column], CO if row == expected else CZ)

        # A3 is not a syndrome coordinate.  The same source does not
        # supply its individual indicator probe; the kernel freedom is
        # real, not merely a renamed bank-diagonal gauge.
        wrong_probe = _cz(5, 6)
        wrong = ([('H', 6)] + wrong_probe + _inverse(source)
                 + wrong_probe + source + [('H', 6)])
        wrong_column = native_columns(7, wrong, columns=(0,))
        self.assertTrue(any(wrong_column[row][0] != (CO if row == 0 else CZ)
                            for row in range(128)))


if __name__ == '__main__':
    unittest.main()
