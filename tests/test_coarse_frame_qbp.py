"""Small ideal-matrix checks of the corrected coarse-frame QBP decoder.

Four/eight-mode fixtures use ideal complex logical coarse matrices and their
literal adjoints. They check interference and classical reconstruction, not
native gate emission, asymptotic resource bounds, or a fine frame compiler.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import direct_real_frame, real_tree_data
from tests.test_reference_state_qbp import H, _coarse_fixture, _fixtures, _observable


ATOL = 5e-12
Y_READOUT = np.array([[1., -1j], [1., 1j]]) / np.sqrt(2.)


def _walsh(height):
    result = np.array([[1.]])
    for _ in range(height):
        result = np.kron(result, H)
    return result


def _interference(coarse, target, observable, walsh):
    """Conditional X/Y distributions, each normalized before the fair setting."""
    size = len(target)
    transform = walsh @ coarse.conj().T
    prepared = np.concatenate((coarse[:, 0], target)) / np.sqrt(2.)
    controlled = np.eye(2 * size, dtype=complex)
    controlled[size:, size:] = observable
    output = np.kron(np.eye(2), transform) @ controlled @ prepared
    probabilities = np.array([
        abs(np.kron(readout, np.eye(size)) @ output).reshape(2, size) ** 2
        for readout in (H, Y_READOUT)
    ])
    return probabilities, transform


def _expected_scores(probabilities, transformed_derivatives):
    size = transformed_derivatives.shape[1]
    differences = probabilities[:, 0] - probabilities[:, 1]
    # Half the shots use each setting. The score itself has factor four.
    return 2 * np.sqrt(size) * (
        transformed_derivatives.real @ differences[0]
        + transformed_derivatives.imag @ differences[1])


def _reverse_tree(angles, leaf_weights):
    """Differentiate the fixed linear leaf functional, never C or the records."""
    size = len(angles) + 1
    incoming, adjoints = np.zeros(2 * size), np.zeros(2 * size)
    incoming[1] = 1.
    adjoints[size:] = np.real(leaf_weights)
    for node, theta in enumerate(angles, start=1):
        incoming[2 * node:2 * node + 2] = incoming[node] * np.array([
            np.cos(theta), np.sin(theta)])
    gradient = np.zeros(size - 1)
    for node in range(size - 1, 0, -1):
        c, s = np.cos(angles[node - 1]), np.sin(angles[node - 1])
        left, right = adjoints[2 * node:2 * node + 2]
        gradient[node - 1] = incoming[node] * (-s * left + c * right)
        adjoints[node] = c * left + s * right
    return gradient


class CoarseFrameQBPTests(unittest.TestCase):
    def test_randomized_quadratures_recover_all_raw_gradients(self):
        for height in (2, 3):
            size, walsh = 1 << height, _walsh(height)
            observable = _observable(size)
            np.testing.assert_allclose(observable.conj().T, observable, atol=ATOL)
            np.testing.assert_allclose(observable @ observable, np.eye(size), atol=ATOL)
            self.assertGreater(np.linalg.norm(observable.imag), .1)
            for angles in _fixtures(height):
                with self.subTest(height=height, angles=angles):
                    data = real_tree_data(angles)
                    derivatives = np.asarray(data.derivatives)
                    frame = direct_real_frame(height, angles)
                    coarse = _coarse_fixture(frame, scale=.5)
                    probabilities, transform = _interference(coarse, data.state, observable, walsh)
                    reference = transform @ coarse[:, 0]
                    response = transform @ observable @ data.state
                    np.testing.assert_allclose(reference, np.ones(size) / np.sqrt(size), atol=ATOL)
                    signs = np.array([1., -1.])[:, None]
                    np.testing.assert_allclose(probabilities[0], abs(reference + signs * response) ** 2 / 4,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(probabilities[1], abs(reference - 1j * signs * response) ** 2 / 4,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(probabilities.sum(axis=(1, 2)), 1., atol=ATOL)
                    transformed = (transform @ derivatives.T).T
                    decoded = _expected_scores(probabilities, transformed)
                    analytic = 2 * np.real(derivatives @ observable @ data.state)
                    np.testing.assert_allclose(decoded, analytic, atol=ATOL, rtol=0)

    def test_depth_records_obey_the_full_frame_distance_bound(self):
        for height in (2, 3):
            size, walsh = 1 << height, _walsh(height)
            for angles in _fixtures(height):
                data = real_tree_data(angles)
                frame = direct_real_frame(height, angles)
                coarse = _coarse_fixture(frame, scale=.5)
                error = np.linalg.norm(coarse - frame, ord=2)
                self.assertLessEqual(error, min(1 / 64, 1 / (4 * np.sqrt(size))))
                np.testing.assert_allclose(coarse.conj().T @ coarse, np.eye(size), atol=ATOL)
                derivatives = np.asarray(data.derivatives)
                for node in range(1, size):
                    np.testing.assert_allclose(derivatives[node - 1],
                                               data.incoming_amplitude[node - 1] * frame[:, marker_label(node, height)],
                                               atol=ATOL, rtol=0)
                transformed = (walsh @ coarse.conj().T @ derivatives.T).T
                for depth in range(height):
                    rows = transformed[(1 << depth) - 1:(1 << (depth + 1)) - 1]
                    bound = 4 * (1 + np.sqrt(size) * error)
                    for quadrature in (rows.real, rows.imag):
                        record_norms = 4 * np.sqrt(size) * np.linalg.norm(quadrature, axis=0)
                        self.assertLessEqual(record_norms.max(), bound + ATOL)
                    self.assertLessEqual(bound, 5 + ATOL)

    def test_histogram_actual_adjoint_and_tree_reverse_equal_record_scores(self):
        for height in (2, 3):
            size, walsh = 1 << height, _walsh(height)
            for angles in (_fixtures(height)[1], np.zeros(size - 1)):
                data = real_tree_data(angles)
                derivatives = np.asarray(data.derivatives)
                coarse = _coarse_fixture(direct_real_frame(height, angles), scale=.5)
                transform = walsh @ coarse.conj().T
                transformed = (transform @ derivatives.T).T
                records = [(setting, x, (-1 if (x + setting) % 3 == 0 else 1))
                           for x in range(size) for setting in (0, 1)]
                records += [(0, 0, 1), (1, 1, -1), (1, size - 1, 1)]
                direct, histogram = np.zeros(size - 1), np.zeros(size, dtype=complex)
                for setting, leaf, sign in records:
                    quadrature = transformed[:, leaf].real if setting == 0 else transformed[:, leaf].imag
                    direct += 4 * np.sqrt(size) * sign * quadrature
                    histogram[leaf] += sign * (1 if setting == 0 else 1j)
                direct /= len(records)
                coefficient = 4 * np.sqrt(size) / len(records)
                leaf_weights = coarse @ walsh @ histogram
                np.testing.assert_allclose(leaf_weights, transform.conj().T @ histogram, atol=ATOL)
                dense = coefficient * np.real(derivatives @ leaf_weights)
                reverse = coefficient * _reverse_tree(angles, leaf_weights)
                np.testing.assert_allclose(direct, dense, atol=ATOL, rtol=0)
                np.testing.assert_allclose(direct, reverse, atol=ATOL, rtol=0)

    def test_complex_phase_requires_y_records_and_actual_coarse_correction(self):
        # W=I, psi=e0, and the root derivative is e2. The imaginary
        # Hermitian observable has identically zero energy on real states.
        height, size, alpha = 2, 4, .01
        data = real_tree_data(np.zeros(size - 1))
        walsh = _walsh(height)
        coarse = np.eye(size, dtype=complex)
        coarse[2, 2] = np.exp(1j * alpha)
        pauli_y = np.array([[0., -1j], [1j, 0.]])
        observable = np.kron(pauli_y, np.eye(2))
        self.assertLess(np.linalg.norm(coarse - np.eye(size), ord=2),
                        min(1 / 64, 1 / (4 * np.sqrt(size))))
        probabilities, transform = _interference(coarse, data.state, observable, walsh)
        derivatives = np.asarray(data.derivatives)
        transformed = (transform @ derivatives.T).T
        np.testing.assert_allclose(_expected_scores(probabilities, transformed), 0., atol=ATOL)
        # Keep the fair setting probability but drop all Y contributions.
        only_x = 2 * np.sqrt(size) * transformed.real @ (probabilities[0, 0] - probabilities[0, 1])
        self.assertAlmostEqual(only_x[0], np.sin(2 * alpha), delta=ATOL)
        self.assertGreater(abs(only_x[0]), .01)
        # Reusing the ideal Walsh/frame scores also loses the actual phase.
        uncorrected = _expected_scores(probabilities, (walsh @ derivatives.T).T)
        self.assertAlmostEqual(uncorrected[0], 2 * np.sin(alpha), delta=ATOL)
        self.assertGreater(abs(uncorrected[0]), .01)


if __name__ == '__main__':
    unittest.main()
