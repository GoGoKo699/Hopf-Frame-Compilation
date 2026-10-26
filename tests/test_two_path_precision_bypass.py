"""Exact complete-space fixtures for balanced two-path recovery bypass.

The identities and resource comparison are analytic. These bounded fixtures
retain literal phases, rejected amplitudes and every arbitrary dirty column;
they do not supply a new precision-scaling construction.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_co_moving_lookup_stage import _ccz_7t
from tests.test_operator_source_compiler import _adjoint


def _difference(left, right):
    return add(left, scale(right, -1))


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = multiply(result, matrix)
    return result


def _canonical_source(u, v):
    accepted = scale(add(u, v), Fraction(1, 2))
    rejected = scale(_difference(u, v), Fraction(1, 2))
    return ([a + b for a, b in zip(accepted, rejected)]
            + [b + a for a, b in zip(accepted, rejected)])


def _recovery(source):
    size = len(source)
    z = [[CO * (1 if i < size // 2 else -1) if i == j else CZ
          for j in range(size)] for i in range(size)]
    return scale(_product(source, z, adjoint(source), z, source), -1)


def _error_gram(output, target):
    error = _difference(output, target)
    return multiply(adjoint(error), error)


def _minus_identity(wire):
    # (Z H)^4 = -I, including its literal sign; twelve H/S Cliffords.
    return [('H', wire), ('S', wire), ('S', wire)] * 4


def _alternating_words(u, v):
    minus = _minus_identity(0)
    return u + _adjoint(v) + u + minus, v + _adjoint(u) + v + minus


class TwoPathPrecisionBypassTests(unittest.TestCase):
    def test_native_order_three_source_noncommuting_projector_and_actual_inverses(self):
        target = [('H', 0), ('T', 0), ('H', 0)]
        order_three = [('CX', 0, 1), ('CX', 1, 0)]
        minus = _minus_identity(0)
        u_word = order_three + target + minus
        v_word = order_three + order_three + target + minus
        controlled = []
        for _, control, target_wire in order_three:
            controlled += ([('H', target_wire)] + _ccz_7t(2, control, target_wire)
                           + [('H', target_wire)])
        source_word = ([('H', 2)] + order_three + controlled + target + minus
                       + [('H', 2)])
        recovered_word = (source_word + [('S', 2)] * 2 + _adjoint(source_word)
                          + [('S', 2)] * 2 + source_word + minus)
        a_word, b_word = _alternating_words(u_word, v_word)
        w, r = native_columns(2, target), native_columns(2, order_three)
        u, v = native_columns(2, u_word), native_columns(2, v_word)
        source = native_columns(3, source_word)
        self.assertEqual(source, _canonical_source(u, v))
        self.assertEqual(_product(r, r, r), eye(4))
        p = _difference(eye(4), scale(add(add(eye(4), r), multiply(r, r)),
                                     Fraction(1, 3)))
        self.assertEqual(multiply(p, p), p)
        self.assertEqual(adjoint(p), p)
        self.assertNotEqual(multiply(w, p), multiply(p, w))
        self.assertEqual(scale(add(u, v), Fraction(1, 2)),
                         multiply(w, _difference(scale(p, Fraction(3, 2)), eye(4))))
        returned = w + zeros(4, 4)
        self.assertEqual(native_columns(3, recovered_word, range(4)), returned)
        self.assertEqual(native_columns(2, a_word), w)
        self.assertEqual(native_columns(2, b_word), w)
        for width, word in ((3, source_word), (3, recovered_word),
                            (2, a_word), (2, b_word)):
            self.assertEqual(native_columns(width, word + _adjoint(word)),
                             eye(1 << width))
        count = lambda word: sum(gate[0] in ('T', 'TDG') for gate in word)
        self.assertEqual((count(source_word), count(recovered_word)), (15, 45))
        self.assertEqual((count(a_word), count(b_word)), (3, 3))
        self.assertEqual(native_columns(1, minus), scale(eye(2), -1))
        self.assertEqual(len(minus), 12)

    def test_full_error_certificate_for_different_unitary_completions(self):
        u = native_columns(2, [('H', 0), ('T', 0), ('CX', 0, 1), ('S', 1)])
        v = native_columns(2, [('T', 1), ('H', 0), ('CX', 1, 0), ('TDG', 0)])
        w = native_columns(2, [('H', 0), ('T', 0), ('H', 0)])
        a = scale(_product(u, adjoint(v), u), -1)
        b = scale(_product(v, adjoint(u), v), -1)
        certificate = scale(add(_error_gram(a, w), _error_gram(b, w)),
                            Fraction(1, 2))
        source = _canonical_source(u, v)
        reject_rotation = native_columns(2, [('H', 1), ('S', 1), ('CX', 1, 0)])
        completion_change = ([row + [CZ] * 4 for row in eye(4)]
                             + [[CZ] * 4 + row for row in reject_rotation])
        other_source = multiply(completion_change, source)
        self.assertNotEqual(source, other_source)
        for candidate in (source, other_source):
            self.assertEqual(multiply(adjoint(candidate), candidate), eye(8))
            self.assertEqual([row[:4] for row in candidate[:4]],
                             scale(add(u, v), Fraction(1, 2)))
            recovered = _recovery(candidate)
            output = [row[:4] for row in recovered]
            self.assertEqual(_error_gram(output, w + zeros(4, 4)), certificate)
            self.assertEqual([row[:4] for row in recovered[:4]],
                             scale(add(a, b), Fraction(1, 2)))
        self.assertNotEqual(_product(adjoint(u), v, adjoint(u), v, adjoint(u), v),
                            eye(4))

    def test_sharp_square_root_two_factor_retains_rejected_error(self):
        u, v = eye(2), native_columns(1, [('T', 0)])
        w = scale(adjoint(v), -1)
        a = scale(_product(u, adjoint(v), u), -1)
        b = scale(_product(v, adjoint(u), v), -1)
        self.assertEqual(a, w)
        self.assertNotEqual(b, w)
        recovered = _recovery(_canonical_source(u, v))
        output = [row[:2] for row in recovered]
        total_error = _error_gram(output, w + zeros(2, 2))
        self.assertEqual(total_error, scale(_error_gram(b, w), Fraction(1, 2)))
        accepted_error = _error_gram(output[:2], w)
        rejected_error = _product(adjoint(output[2:]), output[2:])
        self.assertEqual(accepted_error, rejected_error)
        self.assertNotEqual(rejected_error, zeros(2))
        self.assertEqual(total_error, scale(accepted_error, 2))


if __name__ == '__main__':
    unittest.main()
