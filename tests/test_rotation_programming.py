"""Exact arithmetic certificates for the paired-source rotation bridge.

The oracle reconstructs both moments directly from the paired Majorana
weights. No floating point, inverse angles, or dense matrices are used.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

from compiler_robust_hopf.residual_table_preprocessing import (
    CoefficientInterval,
    residual_rotation_data,
    sqrt_enclosure,
)
from compiler_robust_hopf.rotation_programming import (
    _encode_tail_mean,
    program_rotation,
)


def _point(value):
    return CoefficientInterval(Fraction(value), Fraction(value))


def _product_bounds(left, right):
    values = [a * b for a in (left.lower, left.upper)
              for b in (right.lower, right.upper)]
    return min(values), max(values)


def _tail_weights(q):
    return tuple(Fraction(1, 1 << (j + 1)) for j in range(q - 1)) + (
        Fraction(1, 1 << (q - 1)),)


def _moments_from_signs(signs, q):
    weights = [Fraction(0)] * (2 * (q + 1))
    weights[0] = Fraction(1, 4)
    for j, weight in enumerate(_tail_weights(q)):
        weights[2 * j + 1] = weight / 2
        weights[2 * j + 2] = weight / 4
    scalar = sum(weight * (-1) ** bit for weight, bit in zip(weights, signs))
    tau = sum(weight * (-1) ** (bit + (index > 0 and index % 2 == 0))
              for index, (weight, bit) in enumerate(zip(weights, signs)))
    return scalar, scalar / 2 - tau


def _unit_vector_intervals(a, b, q):
    """Certify (a,b)/sqrt(a*a+b*b), independently of any angle solver."""
    root = sqrt_enclosure(a * a + b * b, q + 80)
    result = []
    for value in (a, b):
        endpoints = value / root.lower, value / root.upper
        result.append(CoefficientInterval(min(endpoints), max(endpoints)))
    return tuple(result)


class RotationProgrammingTests(unittest.TestCase):
    def _check_certificate(self, program):
        q = program.q
        unit = Fraction(1, 1 << q)
        self.assertEqual(len(program.signs), 2 * (q + 1))
        self.assertTrue(all(type(bit) is int and bit in (0, 1)
                            for bit in program.signs))
        self.assertEqual(program.signs[-1], 0)
        self.assertEqual((program.s, program.p), _moments_from_signs(program.signs, q))
        self.assertGreater(program.beta.lower, 0)
        self.assertLess(program.beta.upper, Fraction(1, 3))
        # beta is the positive root of 4*x*x+2*x-1, certified exactly.
        polynomial = lambda x: 4 * x * x + 2 * x - 1
        self.assertLessEqual(polynomial(program.beta.lower), 0)
        self.assertGreaterEqual(polynomial(program.beta.upper), 0)
        target_s = _product_bounds(program.beta, program.cosine)
        target_p = _product_bounds(program.beta, program.sine)
        error_s = max(abs(program.s - endpoint) for endpoint in target_s)
        error_p = max(abs(program.p - endpoint) for endpoint in target_p)
        self.assertLessEqual(error_s ** 2 + error_p ** 2,
                             program.block_error_bound ** 2)
        self.assertLess(program.block_error_bound, Fraction(5, 2) * unit)
        self.assertLessEqual(program.block_error_bound, program.beta.lower / 2)
        self.assertLessEqual(program.beta.upper + program.block_error_bound,
                             Fraction(1, 2))
        self.assertEqual(program.operator_error_bound, 43 * unit)
        self.assertLess(288 * program.block_error_bound ** 2,
                        program.operator_error_bound ** 2)

    def test_exact_circle_points_and_high_precision_moments(self):
        points = {(Fraction(1), Fraction(0)), (Fraction(-1), Fraction(0)),
                  (Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1))}
        for integer in range(-12, 13):
            t = Fraction(integer, 8)
            points.add(((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)))
        for q in (5, 8, 16, 80, 256, 1024):
            for cosine, sine in sorted(points):
                with self.subTest(q=q, cosine=cosine, sine=sine):
                    program = program_rotation(_point(cosine), _point(sine), q)
                    self._check_certificate(program)
                    self.assertEqual(program, program_rotation(_point(cosine), _point(sine), q))
        beyond_float = program_rotation(_point(Fraction(3, 5)), _point(Fraction(4, 5)), 2048)
        self._check_certificate(beyond_float)
        self.assertGreater(beyond_float.operator_error_bound, 0)
        self.assertGreater(beyond_float.operator_error_bound.denominator, 1 << 2000)

    def test_geometric_digit_boundaries_and_clipped_endpoints(self):
        for q in (5, 6, 16, 80, 512):
            weights = _tail_weights(q)
            denominator = 1 << (q - 2)
            endpoint = 1 << (q - 1)
            indices = range(endpoint) if q <= 6 else (0, 1, denominator - 1, endpoint - 1)
            for index in indices:
                # This mean is exactly halfway between grid indices k,k+1.
                middle = 1 - Fraction(2 * index + 1, 2 * denominator)
                perturbation = Fraction(1, 1 << (q + 17))
                for estimate, expected_index in (
                    (middle + perturbation, index), (middle, index + 1),
                    (middle - perturbation, index + 1),
                ):
                    bits, mean = _encode_tail_mean(estimate, q)
                    actual = sum(weight * (-1) ** bit for weight, bit in zip(weights, bits))
                    self.assertEqual(actual, mean)
                    self.assertEqual(mean, 1 - Fraction(expected_index, denominator))
                    self.assertLessEqual(abs(mean - estimate), Fraction(1, 1 << (q - 1)))
            for estimate, value in ((Fraction(2), 1), (Fraction(1), 1),
                                    (Fraction(-1), -1), (Fraction(-2), -1)):
                bits, mean = _encode_tail_mean(estimate, q)
                self.assertEqual(bits, (int(value == -1),) * q)
                self.assertEqual(mean, value)
            zero_bits, zero_mean = _encode_tail_mean(Fraction(0), q)
            self.assertEqual(zero_mean, 0)
            self.assertEqual(zero_bits, (1,) + (0,) * (q - 1))

    def test_head_zero_and_nearby_sides_need_no_exact_angle_sign(self):
        for q in (5, 80, 256):
            epsilon = Fraction(1, 1 << (q + 10))
            for displacement, expected_head in ((-epsilon, 1), (Fraction(0), 0), (epsilon, 0)):
                cosine, sine = _unit_vector_intervals(2 + displacement, Fraction(3), q)
                program = program_rotation(cosine, sine, q)
                self.assertEqual(program.signs[0], expected_head)
                self._check_certificate(program)
            # Interchanging both coefficient signs reaches the other
            # orientation without assuming a canonical angle range.
            cosine, sine = _unit_vector_intervals(Fraction(-2), Fraction(-3), q)
            self._check_certificate(program_rotation(cosine, sine, q))

    def test_residual_half_phases_endpoints_and_interval_clipping(self):
        for q in (5, 16, 80, 256):
            tiny = Fraction(1, 1 << (2 * q + 50))
            for real, imaginary in ((0, 0), (1, 0), (-1, 0), (0, 1),
                                    (Fraction(-3, 4), tiny),
                                    (Fraction(-3, 4), -tiny)):
                data = residual_rotation_data(real, imaginary, q)
                programs = tuple(program_rotation(cosine, sine, q)
                                 for cosine, sine in data.rotations)
                self.assertEqual(programs[0], programs[2])
                for program in programs:
                    self._check_certificate(program)
            width = Fraction(1, 1 << (q + 20))
            clipped = program_rotation(CoefficientInterval(Fraction(1), 1 + width),
                                       CoefficientInterval(-width, Fraction(0)), q)
            self.assertEqual(clipped.cosine, _point(1))
            self._check_certificate(clipped)
            negative = program_rotation(CoefficientInterval(-1 - width, Fraction(-1)),
                                        CoefficientInterval(Fraction(0), width), q)
            self.assertEqual(negative.cosine, _point(-1))
            self._check_certificate(negative)

    def test_rejects_noncertified_types_widths_and_inconsistent_pairs(self):
        for q in (4, -1, True, Fraction(5), 5.0):
            with self.subTest(q=q), self.assertRaises(ValueError):
                program_rotation(_point(1), _point(0), q)
        with self.assertRaises(TypeError):
            program_rotation((1, 1), _point(0), 5)
        with self.assertRaises(TypeError):
            program_rotation(CoefficientInterval(1.0, 1.0), _point(0), 5)
        with self.assertRaises(TypeError):
            program_rotation(CoefficientInterval(True, True), _point(0), 5)
        width = Fraction(1, 1 << 25)
        invalid = ((CoefficientInterval(1 - 2 * width, Fraction(1)), _point(0)),
                   (CoefficientInterval(1 + width, 1 + 2 * width), _point(0)),
                   (_point(1), _point(1)), (_point(0), _point(0)),
                   (_point(Fraction(1, 2)), _point(Fraction(1, 2))))
        for cosine, sine in invalid:
            with self.subTest(cosine=cosine, sine=sine), self.assertRaises(ValueError):
                program_rotation(cosine, sine, 5)


if __name__ == "__main__":
    unittest.main()
