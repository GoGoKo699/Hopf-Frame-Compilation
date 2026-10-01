"""Exact small-field certificate for the complex coarse fixture's error.

Each field element is a four-tuple of Fractions in the ordered basis
(1, sqrt(2), i, i*sqrt(2)). No numerical approximation or native emitter
is used to establish the K2 trace and its rational upper bound.
"""
from __future__ import annotations

from fractions import Fraction as F
import unittest


def _field(*coefficients):
    return tuple(F(value) for value in coefficients) + (F(0),) * (4 - len(coefficients))


def _add(left, right):
    return tuple(a + b for a, b in zip(left, right, strict=True))


def _negate(value):
    return tuple(-coefficient for coefficient in value)


def _product(left, right):
    result = [F(0)] * 4
    for j, a in enumerate(left):
        for k, b in enumerate(right):
            # The low basis-index bit is sqrt(2); the high bit is i.
            factor = (2 if j & k & 1 else 1) * (-1 if j & k & 2 else 1)
            result[j ^ k] += factor * a * b
    return tuple(result)


def _conjugate(value):
    return value[0], value[1], -value[2], -value[3]


def _adjoint(matrix):
    return tuple(tuple(_conjugate(matrix[k][j]) for k in range(2)) for j in range(2))


def _multiply(first, *others):
    result = first
    for other in others:
        result = tuple(tuple(_add(_product(result[j][0], other[0][k]),
                                  _product(result[j][1], other[1][k]))
                             for k in range(2)) for j in range(2))
    return result


class ComplexCoarseCertificateTests(unittest.TestCase):
    def test_exact_k2_trace_certifies_three_k3_error_contributions(self):
        zero, one = _field(), _field(1)
        identity = ((one, zero), (zero, one))
        h_entry = _field(0, F(1, 2))
        h = ((h_entry, h_entry), (h_entry, _negate(h_entry)))
        t = ((one, zero), (zero, _field(0, F(1, 2), 0, F(1, 2))))
        s = ((one, zero), (zero, _field(0, 0, 1)))

        # Independent literal K0 = T H T H T-dagger H T-dagger H.
        k = (
            (_field(F(1, 4), F(1, 2), F(1, 4), F(-1, 4)),
             _field(F(-1, 4), 0, F(-1, 4), F(1, 4))),
            (_field(F(1, 4), 0, F(-1, 4), F(1, 4)),
             _field(F(1, 4), F(1, 2), F(-1, 4), F(1, 4))),
        )
        self.assertEqual(k, _multiply(t, h, t, h, _adjoint(t), h, _adjoint(t), h))
        traces = (_field(F(1, 2), 1), _field(1, F(11, 16)),
                  _field(F(1, 2048), F(2893, 2048)))
        for level, expected_trace in enumerate(traces):
            self.assertEqual(_multiply(k, _adjoint(k)), identity)
            determinant = _add(_product(k[0][0], k[1][1]),
                               _negate(_product(k[0][1], k[1][0])))
            self.assertEqual(determinant, one)
            self.assertEqual(_add(k[0][0], k[1][1]), expected_trace)
            if level < 2:
                conjugated = _multiply(s, k, _adjoint(s))
                k = _multiply(k, conjugated, _adjoint(k), _adjoint(conjugated))

        # For SU(2), ||K2-I||^2 = 2-tr(K2).
        gap = _add(_field(2), _negate(traces[2]))
        self.assertEqual(gap, _field(F(4095, 2048), F(-2893, 2048)))
        self.assertLess(41 ** 2, 2 * 29 ** 2)  # sqrt(2) > 41/29.
        self.assertLess(gap[1], 0)
        upper_bound = gap[0] + gap[1] * F(41, 29)
        self.assertEqual(upper_bound, F(71, 29696))
        self.assertLess(upper_bound, F(1, 400))
        # ||[A,B]-I|| <= 2||A-I||||B-I|| gives delta_K3 < 1/200.
        self.assertLess(2 * upper_bound, F(1, 200))
        self.assertLess(3 * F(1, 200), F(1, 64))


if __name__ == "__main__":
    unittest.main()
