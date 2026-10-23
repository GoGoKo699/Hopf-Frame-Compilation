"""Exact two-clean Gram recovery and complete dirty-sector promotion.

The fixtures check arbitrary rejected-code action and helper inputs.
They do not construct a normalized shallow complete lookup source.
"""
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_gram_four_call_readout import HALF, _partial_gram_source
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _counts, _cz, _x


def _subtract(left, right):
    return add(left, scale(right, -CO))


def _joint_reflection(first, second):
    return ([('S', first)] * 2 + [('S', second)] * 2
            + _cz(first, second))


def _zero_test(first, second, toggle, phase):
    return (_x(first) + _x(second)
            + _private_toffoli(first, second, toggle, phase)
            + _x(second) + _x(first))


def _readout(source, first, second):
    reflection = _joint_reflection(first, second)
    return (source + reflection + _inverse(source) + reflection) * 2


def _controlled_h(control, target):
    rotation = ([('S', target)] * 3
                + [('H', target), ('T', target), ('H', target), ('S', target)])
    return _inverse(rotation) + _cz(control, target) + rotation


class TwoCleanGramRecoveryTests(unittest.TestCase):
    def test_joint_code_reflection_and_returned_zero_test_native_counts(self):
        reflection = _joint_reflection(0, 1)
        for orientation in (reflection, _inverse(reflection)):
            self.assertEqual(_counts(orientation), (0, 7))
            actual = native_columns(2, orientation)
            for column in range(4):
                for row in range(4):
                    self.assertEqual(actual[row][column],
                                     (CO if column == 0 else -CO)
                                     if row == column else CZ)

        zero_test = _zero_test(0, 1, 2, 3)
        for orientation in (zero_test, _inverse(zero_test)):
            self.assertEqual(_counts(orientation), (8, 26))
            actual = native_columns(4, orientation)
            for column in range(16):
                target = column ^ (4 if column & 3 == 0 else 0)
                for row in range(16):
                    self.assertEqual(actual[row][column], CO if row == target else CZ)

    def test_coupled_two_signal_noncommuting_gram_and_complete_error_square(self):
        # One logical bit and two clean signals occupy unequal-size blocks.
        source = [('H', 0), ('H', 1), ('T', 1), ('CX', 0, 1),
                  ('H', 2), ('T', 2), ('CX', 1, 2), ('H', 1)]
        unitary = native_columns(3, source)
        accepted = [row[:2] for row in unitary[:2]]
        gram = multiply(adjoint(accepted), accepted)
        pi = [[CO, CZ], [CZ, CZ]]
        self.assertNotEqual(multiply(gram, pi), multiply(pi, gram))
        complement = _subtract(eye(2), pi)
        delta = _subtract(gram, scale(pi, HALF))
        expected = add(_subtract(scale(multiply(gram, gram), 8), scale(gram, 8)), eye(2))
        readout = _readout(source, 1, 2)
        actual = native_columns(3, readout)
        self.assertEqual([row[:2] for row in actual[:2]], expected)
        embedding = [[CO if row == column else CZ for column in range(2)]
                     for row in range(8)]
        target = _subtract(eye(2), scale(pi, 2))
        error = _subtract(multiply(actual, embedding), multiply(embedding, target))
        first = multiply(multiply(pi, multiply(delta, delta)), pi)
        second = multiply(multiply(complement,
                                   multiply(gram, _subtract(eye(2), gram))), complement)
        self.assertEqual(multiply(adjoint(error), error), scale(add(first, second), 16))
        self.assertEqual(native_columns(3, _inverse(readout)), adjoint(actual))

    def test_exact_joint_recovery_returns_private_fiber_but_not_off_code_query(self):
        # c,d, arbitrary private fiber, a,b.  Coupling a to b leaves the
        # accepted Gram fixed but changes rejected joint-signal action.
        source = _partial_gram_source(private_fiber=True) + _controlled_h(3, 4)
        unitary = native_columns(5, source)
        accepted = [row[:8] for row in unitary[:8]]
        gram = multiply(adjoint(accepted), accepted)
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            for row in range(8):
                self.assertEqual(gram[row][column], HALF if equal and row == column else CZ)
        readout = _readout(source, 3, 4)
        actual = native_columns(5, readout)
        self.assertEqual(native_columns(5, _inverse(readout)), adjoint(actual))
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            for row in range(32):
                self.assertEqual(actual[row][column],
                                 (-CO if equal else CO) if row == column else CZ)
        # The desired negative reflection is -1 here, but the off-code
        # action is +1.  Recovery alone does not remove the clean promise.
        self.assertEqual(actual[16][16], CO)

    def test_four_sector_promotion_cancels_fixed_phase_and_arbitrary_off_code_action(self):
        # x=0, old output=1, source-clean bank=(2,3), new output=4,
        # arbitrary toggle=5 and zero-test private phase helper=6.
        negative_identity = _x(0) + [('S', 0)] * 2 + _x(0) + [('S', 0)] * 2
        source = (negative_identity + [('CX', 0, 1)]
                  + _controlled_h(2, 0) + _controlled_h(3, 1))
        source_matrix = native_columns(4, source)
        for column in range(4):
            target = column ^ (2 if column & 1 else 0)
            for row in range(16):
                self.assertEqual(source_matrix[row][column], -CO if row == target else CZ)
        self.assertEqual(sum(source_matrix[row][4] != CZ for row in range(16)), 2)

        zero_test = _zero_test(2, 3, 5, 6)
        probe = _ccz((5, 1, 4))
        conjugated = source + probe + _inverse(source)

        def echo(word):
            return zero_test + word + _inverse(zero_test) + _inverse(word)

        selective = echo(probe) + echo(conjugated)
        # Cyclic Gray order visits 00,01,11,10 and returns to 00.
        chunks = [[('H', 4)]]
        for wire in (2, 3, 2, 3):
            chunks.extend((selective, _x(wire)))
        chunks.append([('H', 4)])
        promoted = [gate for chunk in chunks for gate in chunk]
        # Evaluate each repeated literal native block once, then compose
        # the complete word.  No ideal selective action is substituted.
        cache = {}

        def evaluate(parts):
            operator = eye(128)
            for part in parts:
                key = tuple(part)
                if key not in cache:
                    cache[key] = native_columns(7, part)
                operator = multiply(cache[key], operator)
            return operator

        actual = evaluate(chunks)
        for column in range(128):
            target = column ^ (16 if column & 1 else 0)
            for row in range(128):
                self.assertEqual(actual[row][column], CO if row == target else CZ)
        self.assertEqual(evaluate([_inverse(chunk) for chunk in reversed(chunks)]),
                         adjoint(actual))
        source_t = max(_counts(source)[0], _counts(_inverse(source))[0])
        source_g = max(_counts(source)[1], _counts(_inverse(source))[1])
        for orientation in (promoted, _inverse(promoted)):
            actual_t, actual_g = _counts(orientation)
            self.assertLessEqual(actual_t, 4 * (4 * source_t + 4 * 8 + 24 + 4))
            self.assertLessEqual(actual_g, 4 * (4 * source_g + 4 * 26 + 40 + 16) + 2)


if __name__ == '__main__':
    unittest.main()
