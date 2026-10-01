"""Finite rational certificates for residual-table preprocessing.

These tests certify algebraic enclosures and finite resource/error
inequalities. They do not simulate a quantum circuit or certify a
general-purpose native single-qubit synthesizer.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

from compiler_robust_hopf.residual_table_preprocessing import (
    residual_rotation_data,
    sqrt_enclosure,
)


def _p2(exponent):
    return Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << -exponent)


def _ends(interval):
    return interval.lower, interval.upper


def _add(left, right):
    return left[0] + right[0], left[1] + right[1]


def _negate(value):
    return -value[1], -value[0]


def _multiply(left, right):
    values = [a * b for a in left for b in right]
    return min(values), max(values)


def _square(value):
    lower, upper = value
    if lower <= 0 <= upper:
        return Fraction(0), max(lower * lower, upper * upper)
    squares = lower * lower, upper * upper
    return min(squares), max(squares)


class ResidualTablePreprocessingTests(unittest.TestCase):
    def assertContains(self, interval, exact):
        self.assertLessEqual(interval[0], exact)
        self.assertGreaterEqual(interval[1], exact)

    def test_integer_square_root_certificates(self):
        values = (Fraction(0), Fraction(1), Fraction(4), Fraction(1, 4),
                  Fraction(2), Fraction(5, 9), Fraction(1, 1 << 257),
                  1 - Fraction(1, 1 << 192))
        for bits in (0, 1, 7, 32, 96, 256):
            for value in values:
                with self.subTest(bits=bits, value=value):
                    interval = sqrt_enclosure(value, bits)
                    self.assertGreaterEqual(interval.lower, 0)
                    self.assertLessEqual(interval.lower ** 2, value)
                    self.assertGreaterEqual(interval.upper ** 2, value)
                    self.assertLessEqual(interval.width, _p2(-bits))
                    if interval.lower ** 2 == value:
                        self.assertEqual(interval.lower, interval.upper)

    def _check_phase_certificate(self, a0, b0, q):
        data = residual_rotation_data(a0, b0, q)
        a, b = data.interior
        norm_squared = a * a + b * b
        self.assertLess(norm_squared, 1)
        self.assertEqual(data.coefficient_error, _p2(-2 * q - 20))
        self.assertEqual(data.working_bits, data.target_bits + 4 * q + 40)
        self.assertEqual(data.completion_error_bound, _p2(-q - 4))
        # The same actual interval objects are used, not independently
        # selected half-phase branches at the two boundaries.
        self.assertIs(data.rotations[0], data.rotations[2])
        for pair in data.rotations:
            for interval in pair:
                self.assertLessEqual(interval.width, _p2(-data.target_bits))
                self.assertGreaterEqual(interval.lower, -1)
                self.assertLessEqual(interval.upper, 1)
        if data.small_radius:
            self.assertLessEqual(norm_squared, _p2(-2 * q - 12))
            self.assertEqual(tuple(tuple(_ends(x) for x in pair)
                                   for pair in data.rotations),
                             (((1, 1), (0, 0)), ((0, 0), (1, 1)),
                              ((1, 1), (0, 0))))
            return

        u = _ends(data.rotations[0][0])
        v = _negate(_ends(data.rotations[0][1]))
        radius = _ends(data.rotations[1][0])
        sine = _ends(data.rotations[1][1])
        self.assertGreaterEqual(radius[0], _p2(-q - 7))
        self.assertContains(_square(radius), norm_squared)
        self.assertContains(_square(sine), 1 - norm_squared)
        self.assertContains(_add(_square(u), _square(v)), 1)
        # These are the literal identities r*w**2=a+ib; no phase
        # alignment or floating-point trigonometric oracle is used.
        real_square = _add(_square(u), _negate(_square(v)))
        imaginary_square = _multiply((2, 2), _multiply(u, v))
        self.assertContains(_multiply(radius, real_square), a)
        self.assertContains(_multiply(radius, imaginary_square), b)

    def test_finite_disk_grid_and_exact_phase_certificates(self):
        grid = [Fraction(j, 8) for j in range(-8, 9)]
        for q in (5, 16, 40):
            for a0 in grid:
                for b0 in grid:
                    if a0 * a0 + b0 * b0 <= 1:
                        with self.subTest(q=q, real=a0, imag=b0):
                            self._check_phase_certificate(a0, b0, q)

    def test_tiny_radius_unit_boundary_and_half_phase_cut(self):
        for q in (5, 16, 64, 128):
            delta = _p2(-2 * q - 20)
            threshold = _p2(-q - 6)
            tiny_imag = _p2(-2 * q - 50)
            cases = ((0, 0), (threshold / 2, 0), (threshold, 0),
                     (2 * threshold, 0), (1, 0), (-1, 0), (0, 1),
                     (0, -1), (1 + delta, 0),
                     (Fraction(-1, 2), tiny_imag),
                     (Fraction(-1, 2), -tiny_imag),
                     (Fraction(-1, 2), 0))
            for a0, b0 in cases:
                with self.subTest(q=q, real=a0, imag=b0):
                    self._check_phase_certificate(a0, b0, q)
            self.assertTrue(residual_rotation_data(threshold, 0, q).small_radius)
            self.assertFalse(residual_rotation_data(2 * threshold, 0, q).small_radius)
            positive_cut = residual_rotation_data(Fraction(-1, 2), tiny_imag, q)
            negative_cut = residual_rotation_data(Fraction(-1, 2), -tiny_imag, q)
            exact_cut = residual_rotation_data(Fraction(-1, 2), 0, q)
            # Outer sine is -v. Its discontinuity is deliberate; both
            # copies change together, preserving the completed SU(2).
            self.assertLess(positive_cut.rotations[0][1].upper, 0)
            self.assertGreater(negative_cut.rotations[0][1].lower, 0)
            self.assertEqual(_ends(exact_cut.rotations[0][0]), (0, 0))
            self.assertEqual(_ends(exact_cut.rotations[0][1]), (-1, -1))

    def test_error_and_small_system_workspace_certificates(self):
        for q in (5, 6, 16, 64, 256, 1024):
            delta = _p2(-2 * q - 20)
            self.assertLess((1 - 2 * delta) * (1 + delta), 1)
            self.assertLessEqual(delta + 2 * delta * (1 + delta), 4 * delta)
            self.assertLessEqual(4 * delta, _p2(-q - 10))
            # sqrt(8)<3, so the following rational surrogate certifies
            # the completion's square-root perturbation term.
            self.assertLess(8, 3 ** 2)
            self.assertLess(4 * delta + 3 * _p2(-q - 10), _p2(-q - 8))
            self.assertLess(_p2(-q - 8) + 2 * _p2(-q - 6), _p2(-q - 4))
            target_bits = q + 20
            h = _p2(-(target_bits + 4 * q + 40))
            threshold = _p2(-q - 6)
            self.assertEqual(32 * h / threshold ** 4, _p2(-target_bits - 11))
        self.assertLess(Fraction(129) + Fraction(1, 16), 130)
        self.assertLess(Fraction(390, 1024), 1)
        for n in range(1, 6):
            for precision in (6, 7, 16, 64, 256):
                sectors = 1 << (n + 2)
                live_dirty = (precision + 11) + 1 + 1
                self.assertLessEqual(sectors, 128)
                self.assertEqual(live_dirty, precision + 13)
                self.assertLessEqual(live_dirty, 2 * (precision + n + 7))

    def test_rejects_uncertifiable_data_and_supports_requested_precision(self):
        with self.assertRaises(TypeError):
            residual_rotation_data(0.5, 0, 16)
        with self.assertRaises(ValueError):
            residual_rotation_data(Fraction(1, 3), 0, 16)
        with self.assertRaises(ValueError):
            residual_rotation_data(1, Fraction(1, 2), 16)
        for bits in (1, 7, 80, 256):
            data = residual_rotation_data(Fraction(-3, 4), Fraction(1, 2),
                                          16, target_bits=bits)
            self.assertEqual(data.target_bits, bits)
            self.assertTrue(all(x.width <= _p2(-bits)
                                for pair in data.rotations for x in pair))


if __name__ == "__main__":
    unittest.main()
