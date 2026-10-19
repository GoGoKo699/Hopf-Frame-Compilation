"""An incomplete native source can still give a complete dirty indicator.

The source changes the address and private helper. Every column of the
actual two-call extraction and its native inverse is checked, retaining
literal phases. The fixture is not an asymptotic synthesis result.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, multiply, native_columns
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _counts, _cz, _x


class WeakerIndicatorInterfacesTests(unittest.TestCase):
    def test_arbitrary_address_helper_gauge_cancels_on_all_columns(self):
        # x=0, Y=(1,2), arbitrary private helper=3. G preserves Y but
        # changes x and W, including noncommuting basis changes and T.
        gauge = [('H', 0), ('H', 0)] + _ccz((1, 2, 0)) + [('H', 0)]
        gauge += [('CX', 0, 3), ('T', 3)]
        indicator = _x(1) + [('CX', 0, 1), ('CX', 0, 2)]
        source = gauge + indicator  # rightmost-first U = I G
        self.assertEqual(_counts(gauge), (8, 14))
        for orientation in (source, _inverse(source)):
            self.assertEqual(_counts(orientation), (8, 20))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))
        source_matrix = native_columns(4, source)
        gauge_matrix = native_columns(4, gauge)
        indicator_matrix = native_columns(4, indicator)
        self.assertEqual(source_matrix, multiply(indicator_matrix, gauge_matrix))
        self.assertEqual(native_columns(4, _inverse(source)), adjoint(source_matrix))

        # Input Y=01 triggers no CCX. The initial H produces a
        # superposition, then CX entangles x,W. The completed source
        # therefore does not return x or W, and is not monomial.
        column = 1 << 1
        supported = [row for row in range(16) if source_matrix[row][column] != CZ]
        self.assertEqual(len(supported), 2)
        self.assertEqual({row & 1 for row in supported}, {0, 1})
        self.assertEqual({(row >> 3) & 1 for row in supported}, {0, 1})

        # V=(4,5) is the specified arbitrary output bank. All of old Y,
        # W and x return after the full wrapper. Nothing is initialized.
        probe = _cz(1, 4) + _cz(2, 5)
        output_h = [('H', 4), ('H', 5)]
        wrapper = output_h + probe + _inverse(source) + probe + source + output_h
        for orientation in (wrapper, _inverse(wrapper)):
            self.assertEqual(_counts(orientation), (16, 56))
            actual = native_columns(6, orientation)
            for column in range(64):
                x = column & 1
                expected = column ^ (1 << (4 + x))
                for row in range(64):
                    self.assertEqual(actual[row][column], CO if row == expected else CZ)
        self.assertEqual(native_columns(6, _inverse(wrapper)),
                         adjoint(native_columns(6, wrapper)))


if __name__ == '__main__':
    unittest.main()
