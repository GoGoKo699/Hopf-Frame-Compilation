"""Exact guards for joint order interference and two failed bank-moving pairs.

The controlled anchor and relative phase are literal native fixtures.
The joint matrix fixture checks the full logical control word with its
separately paid primitives. The p=1 failures are finite evidence, not an
all-size exclusion or a qualifying complete indicator.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)
from tests.test_boolean_reflection_gram_source import _cs
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_signal_interleaved_kernel_source import (
    _p1_parts, _root_powers, _two_signal_control,
)


HALF = CQ2(Fraction(1, 2))
QUARTER = CQ2(Fraction(1, 4))
NORMALIZATION = CQ2(Q2(0, Fraction(1, 4)))


def _x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _s0_word(c, l):
    return ([('S', l)] * 2 + [('S', l)] * 3
            + [('CX', c, l), ('S', l)] + [('S', c)] * 3)


def _controlled_s0(a, c, l, helper):
    return ([('H', l), ('CX', a, l), ('H', l)]
            + [('S', l)] * 3 + _private_toffoli(a, c, l, helper)
            + [('S', l)] + _inverse(_cs(a, c)))


def _compose(stages):
    result = eye(len(stages[0]))
    for stage in stages:
        result = multiply(stage, result)
    return result


def _joint_stages(r, s, minus=False):
    bank_width = (len(r) - 1).bit_length()
    a, b = bank_width, bank_width + 1
    width = bank_width + 2
    stages = [
        native_columns(width, [('H', a), ('H', b)]),
        _two_signal_control(r, 0),
        native_columns(width, _x(b) + [('CX', a, b)]),
        _two_signal_control(s, 1),
        native_columns(width, [('CX', a, b)] + _x(b) + _x(a)),
        _two_signal_control(r, 0),
        native_columns(width, _x(a)),
    ]
    if minus:
        stages.append(native_columns(width, [('S', a), ('S', a)]))
    stages.append(native_columns(width, [('CX', a, b), ('H', a)]))
    return stages


def _path_stages(r, s, minus=False):
    """Independently assemble the four paths and their Bell decoder."""
    size = len(r)
    bank_width = (size - 1).bit_length()
    a, b = bank_width, bank_width + 1
    paths = zeros(4 * size)
    # The binary label is a+2b: paths 01 and 10 both carry R.
    for label, word in enumerate((multiply(r, s), r, r, multiply(s, r))):
        for row in range(size):
            for column in range(size):
                paths[label * size + row][label * size + column] = word[row][column]
    decoder = ([('S', a), ('S', a)] if minus else [])
    decoder += [('CX', a, b), ('H', a)]
    return [native_columns(bank_width + 2, [('H', a), ('H', b)]),
            paths, native_columns(bank_width + 2, decoder)]


def _accepted(source, size):
    return [row[:size] for row in source[:size]]


def _p1_pair(address, broad=False):
    e, f = _p1_parts(address)
    p0 = native_columns(4, [('Z', 0), ('Z', 2)])
    d = multiply(adjoint(f), p0)
    if broad:
        encoder = [('CX', 0, 2), ('CX', 0, 1), ('CX', 0, 3),
                   ('CX', 1, 0), ('CX', 3, 0)]
        h = native_columns(4, [('H', wire) for wire in range(4)])
        swap = native_columns(4, [gate for first, second in ((0, 2), (1, 3))
                                 for gate in (('CX', first, second),
                                              ('CX', second, first),
                                              ('CX', first, second))])
        v, a = multiply(multiply(h, d), h), multiply(multiply(swap, d), swap)
    else:
        encoder = [('CX', 3, 0), ('CX', 3, 2), ('CX', 0, 2)]
        v, a = e, d
    r0 = native_columns(4, encoder + [('S', 0)] * 2 + _inverse(encoder))
    s0 = native_columns(4, encoder + _s0_word(2, 0) + _inverse(encoder))
    return (multiply(multiply(adjoint(v), r0), v),
            multiply(multiply(adjoint(a), s0), a))


class CoupledOrderGramSourceTests(unittest.TestCase):
    def test_literal_controlled_anchor_returns_arbitrary_helper_and_actual_reverse(self):
        word = _controlled_s0(0, 1, 2, 3)
        self.assertEqual(_counts(word), (11, 19))
        actual = native_columns(4, word)
        for column in range(16):
            a, c, l = (column >> wire & 1 for wire in range(3))
            target = column ^ (4 if a and c else 0)
            phase = -CO if a and not c and l else CO
            for row in range(16):
                self.assertEqual(actual[row][column], phase if row == target else CZ)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(actual))
        reverse_counts = _counts(_inverse(word))
        self.assertLessEqual(reverse_counts[0], 11)
        self.assertLessEqual(reverse_counts[1], 19)
        r0 = native_columns(2, [('S', 0), ('S', 0)])
        s0 = native_columns(2, _s0_word(1, 0))
        self.assertEqual(multiply(s0, s0), eye(4))
        product = multiply(r0, s0)
        charge = native_columns(2, [('Z', 1)])
        self.assertEqual(multiply(product, product), charge)
        self.assertEqual(multiply(r0, charge), multiply(charge, r0))
        self.assertEqual(multiply(s0, charge), multiply(charge, s0))
        self.assertEqual(multiply(multiply(s0, r0), s0), multiply(r0, charge))

    def test_joint_word_has_complete_order_polynomial_and_actual_reverse(self):
        r = native_columns(1, [('S', 0), ('S', 0)])
        s = native_columns(1, [('H', 0)])
        square = multiply(multiply(r, s), multiply(r, s))
        real_square = scale(add(square, adjoint(square)), HALF)
        for minus in (False, True):
            stages = _joint_stages(r, s, minus)
            actual = _compose(stages)
            self.assertEqual(actual, _compose(_path_stages(r, s, minus)))
            accepted = _accepted(actual, 2)
            sign = -CO if minus else CO
            expected = scale(add(multiply(r, s), scale(multiply(s, r), sign)),
                             NORMALIZATION)
            self.assertEqual(accepted, expected)
            self.assertEqual(multiply(adjoint(accepted), accepted),
                             scale(add(eye(2), scale(real_square, sign)), QUARTER))
            actual_reverse = _compose([adjoint(stage) for stage in reversed(stages)])
            self.assertEqual(actual_reverse, adjoint(actual))
            self.assertEqual(multiply(actual_reverse, actual), eye(8))
        signal_word = ([('H', 1), ('H', 2)] + _x(2) + [('CX', 1, 2)]
                       + [('CX', 1, 2)] + _x(2) + _x(1) + _x(1)
                       + [('CX', 1, 2), ('H', 1)])
        self.assertEqual(_counts(signal_word), (0, 22))

    def test_actual_p1_pairs_have_correct_base_and_distinct_nonzero_failures(self):
        for broad, receipt, error_squared in (
                (False, QUARTER, CQ2(2)),
                (True, CQ2(Fraction(3, 8)), CO)):
            for address in (0, 1):
                r, s = _p1_pair(address, broad)
                self.assertEqual(multiply(r, r), eye(16))
                self.assertEqual(multiply(s, s), eye(16))
                self.assertEqual(r, adjoint(r))
                self.assertEqual(s, adjoint(s))
                source = _compose(_path_stages(r, s))
                accepted = _accepted(source, 16)
                expected = scale(add(multiply(r, s), multiply(s, r)), NORMALIZATION)
                self.assertEqual(accepted, expected)
                gram = multiply(adjoint(accepted), accepted)
                self.assertEqual(multiply(gram, r), multiply(r, gram))
                self.assertEqual(multiply(gram, s), multiply(s, gram))
                square = multiply(multiply(r, s), multiply(r, s))
                self.assertEqual(gram, scale(add(eye(16),
                                                scale(add(square, adjoint(square)), HALF)),
                                            QUARTER))
                target_word = [('Z', 0), ('Z', 2)]
                if address:
                    target_word += [('Z', 1), ('Z', 3)]
                target = native_columns(4, target_word)
                target_gram = scale(add(eye(16), target), QUARTER)
                if not address:
                    self.assertEqual(square, target)
                    self.assertEqual(gram, target_gram)
                else:
                    self.assertEqual(gram[0][0], receipt)
                    self.assertEqual(target_gram[0][0], HALF)
                    difference = [square[row][0] - target[row][0] for row in range(16)]
                    self.assertEqual(sum(value.conj() * value for value in difference),
                                     error_squared)

    def test_uncontrolled_relative_phase_is_native_on_coherent_address_and_helper(self):
        # Wires: x,z0,z1,t0,t1,private helper. Neither x nor helper is fixed.
        qp = [('CX', 2, 4)]
        qb = _private_toffoli(0, 2, 3, 5)
        word = (qp + [('TDG', 3), ('TDG', 4)] + qb
                + [('T', 3), ('T', 4)] + _inverse(qb) + qp)
        self.assertEqual(_counts(word), (20, 22))
        actual = native_columns(6, word)
        powers = _root_powers()
        for column in range(64):
            x, z1 = column & 1, column >> 2 & 1
            t = column >> 3 & 3
            shifted = t ^ (z1 << 1)
            exponent = ((shifted ^ (x * z1)).bit_count()
                        - shifted.bit_count()) % 8
            phase = powers[exponent]
            for row in range(64):
                self.assertEqual(actual[row][column], phase if row == column else CZ)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))


if __name__ == '__main__':
    unittest.main()
