"""Complete native fixtures for the dirty-tree branch-difference contract.

These bounded checks certify operator cancellation, not a growing-arity
source or an unrestricted depth theorem. Every physical column is kept.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, eye, native_columns
from tests.test_majorana_walsh_query import _inverse


def _x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


class DirtyTreeDifferentialTests(unittest.TestCase):
    def test_unreturned_private_work_and_selected_parent_cancellation(self):
        # Address 0 selects parent 1 or 2. Children are 3,4; arbitrary
        # private wires 5,6 undergo a common right action R before F.
        width = 7
        junk = [('H', 5), ('T', 3), ('CX', 3, 5),
                ('H', 6), ('TDG', 4), ('CX', 4, 6)]
        ideal_stage = [('CX', 1, 3), ('CX', 2, 4)]
        stage = junk + ideal_stage
        previous = _x(1) + [('CX', 0, 1), ('CX', 0, 2)]
        expected = previous + _x(3) + [('CX', 0, 3), ('CX', 0, 4)]
        actual = _inverse(stage) + previous + stage
        matrix = native_columns(width, stage)
        self.assertEqual(native_columns(width, _inverse(stage)), adjoint(matrix))
        self.assertEqual(native_columns(width, actual),
                         native_columns(width, expected))
        # Common left action is not interchangeable with common right action.
        wrong_stage = ideal_stage + junk
        wrong = _inverse(wrong_stage) + previous + wrong_stage
        self.assertNotEqual(native_columns(width, wrong),
                            native_columns(width, expected))
        # All-zero input has private-work superposition after the source.
        self.assertTrue(any(matrix[row][0] != CZ and row & (3 << 5)
                            for row in range(1 << width)))

        # A one-leaf complete echo: root 0, child 1, private 2, output 3.
        # It returns the root/child/private wires on every initial input.
        small_stage = [('H', 2), ('T', 1), ('CX', 1, 2), ('CX', 0, 1)]
        path = _inverse(small_stage) + _x(0) + small_stage
        copy = [('CX', 1, 3)]
        echo = copy + path + copy + _inverse(path)
        complete = native_columns(4, echo)
        self.assertEqual(complete, native_columns(4, _x(3)))
        self.assertEqual(native_columns(4, _inverse(echo)), adjoint(complete))

    def test_independent_child_angles_add_in_the_complete_operator_norm(self):
        # (HS)^3 is the literal scalar omega. Combined with S_X or its
        # actual inverse it gives exp(+-i*pi*X/4), with no scalar dropped.
        omega_word = [('H', 0), ('S', 0)] * 3
        omega = native_columns(1, [('T', 0)])[1][1]
        self.assertEqual(native_columns(1, omega_word),
                         [[omega, CZ], [CZ, omega]])
        plus = omega_word + [('H', 0), ('S', 0), ('S', 0), ('S', 0), ('H', 0)]
        minus = [(name, 1) for name, _ in _inverse(omega_word)]
        minus += [('H', 1), ('S', 1), ('H', 1)]
        word = plus + minus
        matrix = native_columns(2, word)
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(matrix))
        # The common X eigenbasis displays all four errors. Their squared
        # norms are 0,2,2,0: total |theta|=pi/2 gives error sqrt(2),
        # whereas each separate child has the smaller error sqrt(2-sqrt(2)).
        h = [('H', 0), ('H', 1)]
        diagonal = native_columns(2, h + word + h)
        self.assertTrue(all(diagonal[r][c] == CZ for r in range(4)
                            for c in range(4) if r != c))
        errors = [(diagonal[j][j] - CO).abs2() for j in range(4)]
        self.assertEqual(errors, [CO.real * value for value in (0, 2, 2, 0)])
        one = native_columns(1, [('H', 0)] + plus + [('H', 0)])
        self.assertLess((one[0][0] - CO).abs2(), CO.real * 2)
        self.assertEqual(native_columns(2, word + _inverse(word)), eye(4))


if __name__ == '__main__':
    unittest.main()
