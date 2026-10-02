"""Bounded actual-source checks of the phase-corrected pi/3 radial filter.

The source core is the same joint Clifford-algebra representation used by
test_flag_echo, with its literal Q, actual inverse, and arbitrary dirty
columns. These are 16-dimensional operator checks, not a new native gate
emitter. The six-rotation perturbation check verifies a continuous operator
norm budget; native phase-synthesis resources are an analytic interface.
Precision-budget inequalities below use exact rational arithmetic.
"""
from __future__ import annotations

from fractions import Fraction
import math
import unittest

import numpy as np

try:
    from .test_flag_echo import _fixture
    from .test_hopf_error_accumulation import I, Z
except ImportError:
    from test_flag_echo import _fixture
    from test_hopf_error_accumulation import I, Z


ATOL = 1e-11


def _phases(embedding):
    projector = embedding @ embedding.conj().T
    z = np.exp(1j * math.pi / 3)
    return z, np.eye(16) + (z - 1) * projector


def _flag_phase_word(errors=(0.0, 0.0, 0.0)):
    za = np.kron(I, np.kron(Z, np.eye(4)))
    zb = np.kron(Z, np.eye(8))
    result = np.eye(16, dtype=complex)
    for pauli, error in zip((za, zb, za @ zb), errors):
        angle = math.pi / 12 + error
        result = (math.cos(angle) * np.eye(16) + 1j * math.sin(angle) * pauli) @ result
    return result


class RadialFilterTests(unittest.TestCase):
    def test_actual_source_filter_controls_all_dirty_columns(self):
        for width, angle in ((4, math.pi / 4), (5, 0.23), (7, 0.51)):
            cosine, sine, _, amplified, _, embedding, polar, _ = _fixture(width, angle)
            radius_squared = cosine * cosine + sine * sine
            defect = 1 - radius_squared
            amplitude = math.sqrt(radius_squared) * (3 - radius_squared) / 2
            rejected_probability = defect * defect * (3 + defect) / 4
            self.assertAlmostEqual(1 - amplitude * amplitude, rejected_probability, delta=ATOL)
            z, phase = _phases(embedding)
            filtered = z ** -2 * amplified @ phase @ amplified.conj().T @ phase @ amplified
            output = filtered @ embedding
            np.testing.assert_allclose(embedding.conj().T @ output,
                                       amplitude * (1 + z.conjugate() * rejected_probability) * polar,
                                       atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.linalg.norm(output[4:], 2),
                                   rejected_probability ** 1.5, delta=ATOL)
            t = rejected_probability
            error_squared = t * t * (3 + t) / (2 * (1 + math.sqrt(1 - t) * (1 + t / 2)))
            actual_error = np.linalg.norm(output - embedding @ polar, 2)
            self.assertAlmostEqual(actual_error, math.sqrt(error_squared), delta=ATOL)
            self.assertLessEqual(error_squared, 2 * defect ** 4)
            np.testing.assert_allclose(filtered.conj().T @ output, embedding, atol=ATOL, rtol=0)
            # A forward call in place of the required inverse changes the
            # target even when the original source has small radial error.
            wrong = z ** -2 * amplified @ phase @ amplified @ phase @ amplified
            self.assertGreater(np.linalg.norm((wrong - filtered) @ embedding, 2), 0.1)

    def test_literal_global_phase_and_six_rotation_error_budget(self):
        _, _, _, amplified, _, embedding, _, _ = _fixture(7, 0.51)
        z, phase = _phases(embedding)
        exact = z ** -2 * amplified @ phase @ amplified.conj().T @ phase @ amplified
        tilde = _flag_phase_word()
        np.testing.assert_allclose(tilde, np.exp(-1j * math.pi / 12) * phase,
                                   atol=ATOL, rtol=0)
        uncorrected = amplified @ tilde @ amplified.conj().T @ tilde @ amplified
        np.testing.assert_allclose(-1j * uncorrected, exact, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(uncorrected - exact, 2), math.sqrt(2), delta=ATOL)
        for tolerance in (1e-3, 1 / 256):
            early = _flag_phase_word((tolerance / 6, -tolerance / 6, tolerance / 6))
            late = _flag_phase_word((-tolerance / 6, tolerance / 6, tolerance / 6))
            approximate = -1j * amplified @ late @ amplified.conj().T @ early @ amplified
            self.assertLessEqual(np.linalg.norm(approximate - exact, 2), tolerance + ATOL)

    def test_half_logarithmic_precision_budget_in_exact_arithmetic(self):
        dimensions = tuple(range(1, 34)) + (63, 64, 65, 127, 128, 129, 257, 1024)
        for n in dimensions:
            cap = 0
            while 4 ** cap < 8 * n:
                cap += 1
            old_cap = (8 * n - 1).bit_length()
            self.assertLess(4 ** (cap - 1), 8 * n)
            self.assertGreaterEqual(4 ** cap, 8 * n)
            for precision in (6, 12, 64):
                widths = [precision + 4 + min(n - depth, cap) for depth in range(n)]
                square_sum = sum((Fraction(1, 4 ** width) for width in widths), Fraction())
                self.assertLessEqual(square_sum, Fraction(11, 6144 * 4 ** precision))
                # Square only after checking the positive remaining budget:
                # 8 sqrt(2S) + 240S < 2^(-L), with no floating square root.
                remaining = Fraction(1, 2 ** precision) - 240 * square_sum
                self.assertGreater(remaining, 0)
                self.assertLess(128 * square_sum, remaining * remaining)
                for depth, width in enumerate(widths):
                    self.assertLessEqual(width, precision + n - depth + 4)
                    self.assertLessEqual(width, precision + 4 + min(n - depth, old_cap))
                stopped = min(n, cap)
                self.assertEqual(sum(widths),
                                 n * (precision + 4) + n * stopped - stopped * (stopped - 1) // 2)
                self.assertLessEqual(sum(widths), n * (precision + 4 + cap))


if __name__ == '__main__':
    unittest.main()
