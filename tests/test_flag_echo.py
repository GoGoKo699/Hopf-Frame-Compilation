"""Small literal-source algebra checks for the two-half-angle flag composite.

Every Q uses one common representation of its M, cosine-mask, and
sine-mask Majorana vectors. The squared gadget reuses that same dirty core
and both flags. This is a Clifford-algebra reduction of the actual source
identities, not a new emitted native gate word. Dense matrices are 16 by
16; all phases and the inverse are retained.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

try:
    from .test_hopf_error_accumulation import I, X, Y, Z, _amplify, _mask_vector, _query
except ImportError:
    from test_hopf_error_accumulation import I, X, Y, Z, _amplify, _mask_vector, _query


ATOL = 8e-12


def _fixture(width, angle):
    amplitudes = np.array([2.0 ** (-(j + 1) / 2) for j in range(width - 1)]
                          + [2.0 ** (-(width - 1) / 2)])
    cosine_mask, cosine = _mask_vector(math.cos(angle), amplitudes)
    sine_mask, sine = _mask_vector(math.sin(angle), amplitudes)
    vectors = np.array([amplitudes, cosine_mask, sine_mask])
    basis, _ = np.linalg.qr(vectors.T, mode='reduced')
    coordinates = vectors @ basis
    np.testing.assert_allclose(coordinates @ coordinates.T, vectors @ vectors.T,
                               atol=ATOL, rtol=0)
    source, nc, ns = [sum(coefficient * gamma for coefficient, gamma in zip(row, (X, Y, Z)))
                      for row in coordinates]
    query = _query(source, nc, ns)
    amplified, reflection = _amplify(query)
    embedding = np.eye(16, dtype=complex)[:, :4]
    radius_squared = cosine * cosine + sine * sine
    polar = np.kron((cosine * I - 1j * sine * Y) / math.sqrt(radius_squared), I)
    mask_inner_product = float(np.dot(cosine_mask, sine_mask))
    return cosine, sine, query, amplified, reflection, embedding, polar, mask_inner_product


class FlagEchoTests(unittest.TestCase):
    def test_equal_rows_leave_linear_leakage_and_zero_polar_angle_error(self):
        for width in (4, 8, 12):
            cosine, sine, _, amplified, _, embedding, polar, _ = _fixture(width, math.pi / 4)
            self.assertEqual(cosine, sine)
            defect = 1 - cosine * cosine - sine * sine
            output = amplified @ amplified @ embedding
            ideal = polar @ polar
            accepted = embedding.conj().T @ output
            np.testing.assert_allclose(accepted, (1 - defect * defect) * ideal,
                                       atol=ATOL, rtol=0)
            # The accepted block has the exact target polar angle, but its
            # omitted flag block contains a literal defect times identity.
            np.testing.assert_allclose(output[8:12], defect * np.eye(4),
                                       atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.linalg.norm(output - embedding @ ideal, 2),
                                   math.sqrt(2) * abs(defect), delta=ATOL)
            self.assertAlmostEqual(np.linalg.norm(output[4:], 2),
                                   abs(defect) * math.sqrt(2 - defect * defect), delta=ATOL)

    def test_generic_composite_accepted_rejected_blocks_and_actual_inverse(self):
        for width, angle in ((5, 0.23), (7, 0.51), (9, 0.91)):
            cosine, sine, query, amplified, reflection, embedding, polar, _ = _fixture(width, angle)
            normalized = -reflection @ amplified
            composite = reflection @ normalized @ reflection @ normalized
            np.testing.assert_allclose(composite, amplified @ amplified, atol=ATOL, rtol=0)
            output = composite @ embedding
            block = embedding.conj().T @ query @ embedding
            second = embedding.conj().T @ query @ query @ embedding
            np.testing.assert_allclose(second, (cosine * cosine - sine * sine) * np.eye(4),
                                       atol=ATOL, rtol=0)
            defect = 1 - cosine * cosine - sine * sine
            accepted = defect * defect * second + 4 * (1 + defect) * block @ block
            np.testing.assert_allclose(embedding.conj().T @ output, accepted,
                                       atol=ATOL, rtol=0)
            raw_rejected = defect * (-query @ query @ embedding
                                    + 4 * query @ embedding @ (block - block.conj().T @ second))
            np.testing.assert_allclose(output[4:], raw_rejected[4:], atol=ATOL, rtol=0)
            phi = math.atan2(sine, cosine)
            full_error = abs(defect) * math.sqrt(
                2 * (math.sin(2 * phi) ** 2 + defect * math.cos(2 * phi) ** 2))
            self.assertAlmostEqual(np.linalg.norm(output - embedding @ polar @ polar, 2),
                                   full_error, delta=ATOL)
            # Evaluate the algebraic leakage norm without subtracting two
            # nearly equal order-one numbers at fine precision.
            leakage_squared = defect * defect * (
                (2 * defect - defect ** 4) * math.cos(2 * phi) ** 2
                + (2 - defect * defect) * math.sin(2 * phi) ** 2)
            self.assertAlmostEqual(np.linalg.norm(output[4:], 2),
                                   math.sqrt(leakage_squared), delta=ATOL)
            np.testing.assert_allclose(amplified.conj().T @ amplified.conj().T @ output,
                                       embedding, atol=ATOL, rtol=0)

    def test_literal_near_zero_rows_have_only_a_scoped_improvement(self):
        for exponent in (3, 4, 5):
            sine_target = 2.0 ** -exponent
            cosine, sine, _, amplified, _, embedding, polar, _ = _fixture(
                2 * exponent + 1, math.asin(sine_target))
            self.assertEqual(cosine, 1)
            self.assertEqual(sine, sine_target)
            squared_error = np.linalg.norm(amplified @ amplified @ embedding
                                           - embedding @ polar @ polar, 2)
            expected = sine ** 3 * math.sqrt(2 * (3 - sine * sine) / (1 + sine * sine))
            self.assertAlmostEqual(squared_error, expected, delta=ATOL)
            single_error = np.linalg.norm(amplified @ embedding - embedding @ polar, 2)
            self.assertLess(squared_error, single_error)

    def test_diagonal_flag_variants_and_equal_mask_exception(self):
        za = np.kron(I, np.kron(Z, np.eye(4)))
        zb = np.kron(Z, np.eye(8))
        zt = np.kron(np.eye(4), np.kron(Z, I))
        logical_k = np.kron(-1j * Y, I)
        full_k = np.kron(np.eye(4), logical_k)
        for width, angle in ((4, math.pi / 4), (7, 0.51)):
            cosine, sine, query, amplified, _, embedding, polar, inner = _fixture(width, angle)
            defect = 1 - cosine * cosine - sine * sine
            sin_twice = math.sin(2 * math.atan2(sine, cosine))
            lower = abs(defect) * math.sqrt(2 * (1 - abs(sin_twice)))
            for flag, coefficient in ((za, 0), (zb, 2 * cosine * sine - inner),
                                      (za @ zb, inner)):
                second = embedding.conj().T @ query @ flag @ query @ embedding
                np.testing.assert_allclose(second, coefficient * logical_k, atol=ATOL, rtol=0)
                composite = amplified @ flag @ amplified @ flag.conj().T
                output = composite @ embedding
                accepted = ((1 - defect * defect) * polar @ polar
                            + defect * defect * second)
                np.testing.assert_allclose(embedding.conj().T @ output, accepted,
                                           atol=ATOL, rtol=0)
                error = np.linalg.norm(output - embedding @ polar @ polar, 2)
                self.assertAlmostEqual(error, abs(defect) * math.sqrt(
                    max(0, 2 * (1 - coefficient * sin_twice))), delta=ATOL)
                self.assertGreaterEqual(error + ATOL, lower)
            # This is an exact full-space exception, not merely clean-column
            # cancellation; a generic angle fails the same equality.
            double_flag = amplified @ za @ zb @ amplified @ za @ zb
            if cosine == sine:
                np.testing.assert_allclose(double_flag, full_k, atol=ATOL, rtol=0)
            else:
                self.assertGreater(np.linalg.norm(double_flag - full_k, 2), 0.1)
            np.testing.assert_allclose(amplified @ zt @ amplified.conj().T @ zt,
                                       amplified @ za @ amplified @ za,
                                       atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
