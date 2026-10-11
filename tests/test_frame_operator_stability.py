"""Full-operator diagnostics for HOPF_INTERFACE.md section 4.2.

The derivative is accumulated by the independent product rule for addressed
rotation layers. It is then compared with the nested-support tangent formula.
These bounded dense-matrix checks do not model approximate source leakage or
establish any new compiler resource bound.
"""
from __future__ import annotations

import math
import unittest
from dataclasses import dataclass

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix


ATOL = 4e-12


@dataclass(frozen=True)
class Differential:
    frame: np.ndarray
    product_derivative: np.ndarray
    lower_tangent: np.ndarray
    psd: np.ndarray
    caps: np.ndarray

    @property
    def tangent(self) -> np.ndarray:
        return self.lower_tangent - self.lower_tangent.T


def _differential(n: int, angles: np.ndarray, velocities: np.ndarray) -> Differential:
    dimension = 1 << n
    product = np.eye(dimension)
    derivative = np.zeros((dimension, dimension))
    lower_tangent = np.zeros((dimension, dimension))
    psd = np.zeros((dimension, dimension))
    caps = []
    for depth in range(n):
        layer = np.eye(dimension)
        layer_derivative = np.zeros((dimension, dimension))
        edge_map = np.zeros((dimension, dimension))
        offset = (1 << depth) - 1
        caps.append(np.max(np.abs(velocities[offset:offset + (1 << depth)])))
        for prefix in range(1 << depth):
            old = prefix << (n - depth)
            new = old + (1 << (n - depth - 1))
            theta = angles[offset + prefix]
            velocity = velocities[offset + prefix]
            cosine, sine = math.cos(theta), math.sin(theta)
            pair = np.ix_([old, new], [old, new])
            layer[pair] = [[cosine, -sine], [sine, cosine]]
            layer_derivative[pair] = velocity * np.array(
                [[-sine, -cosine], [cosine, -sine]]
            )
            edge_map[new, old] = velocity
        piece = edge_map @ product
        lower_tangent += piece
        psd += piece.T @ piece
        # This computation does not use the tangent or PSD identities.
        derivative = layer_derivative @ product + layer @ derivative
        product = layer @ product
    return Differential(product, derivative, lower_tangent, psd, np.asarray(caps))


def _cap_bound(caps: np.ndarray) -> float:
    bound = 0.0
    for cap in caps:
        bound = (bound + math.sqrt(bound**2 + 4 * cap**2)) / 2
    return bound


class FrameOperatorStabilityTests(unittest.TestCase):
    def test_product_rule_derivative_and_psd_at_unrestricted_and_corner_angles(self):
        rng = np.random.default_rng(20261011)
        for n in range(1, 6):
            count = (1 << n) - 1
            for corner in (False, True):
                angles = (
                    rng.choice([0.0, math.pi / 2, -math.pi], count)
                    if corner else rng.uniform(-math.pi, math.pi, count)
                )
                result = _differential(n, angles, rng.uniform(-1, 1, count))
                with self.subTest(n=n, corner=corner):
                    # Anchor the independent layer construction to the
                    # production recursive frame and its marker convention.
                    np.testing.assert_allclose(
                        result.frame, real_frame_matrix(angles), atol=ATOL, rtol=0
                    )
                    np.testing.assert_allclose(
                        result.product_derivative,
                        result.frame @ result.tangent,
                        atol=ATOL, rtol=0,
                    )
                    np.testing.assert_allclose(
                        result.lower_tangent.T @ result.lower_tangent,
                        result.psd, atol=ATOL, rtol=0,
                    )
                    actual = np.linalg.norm(result.product_derivative, 2)
                    psd_bound = 2 * math.sqrt(np.linalg.eigvalsh(result.psd)[-1])
                    self.assertLessEqual(actual, psd_bound + ATOL)
                    self.assertLessEqual(actual, _cap_bound(result.caps) + ATOL)

    def test_uniform_depth_magnitudes_attain_full_spectral_recurrence(self):
        rng = np.random.default_rng(20261012)
        for n in range(1, 6):
            for zero_depth in (None, 0, n - 1):
                caps = rng.uniform(0.05, 1, n)
                if zero_depth is not None:
                    caps[zero_depth] = 0
                velocities = np.concatenate([
                    cap * rng.choice([-1, 1], size=1 << depth)
                    for depth, cap in enumerate(caps)
                ])
                angles = rng.uniform(-math.pi, math.pi, (1 << n) - 1)
                result = _differential(n, angles, velocities)
                spectrum = np.zeros(1)
                for cap in caps:
                    root = np.sqrt(spectrum**2 + 4 * cap**2)
                    spectrum = np.concatenate([(spectrum + root) / 2,
                                               (spectrum - root) / 2])
                with self.subTest(n=n, zero_depth=zero_depth):
                    np.testing.assert_allclose(
                        np.linalg.eigvalsh(1j * result.tangent),
                        np.sort(spectrum), atol=ATOL, rtol=0,
                    )
                    self.assertAlmostEqual(
                        np.linalg.norm(result.product_derivative, 2),
                        _cap_bound(caps), delta=ATOL,
                    )

    def test_selected_root_to_leaf_path_attains_psd_sharpness_formula(self):
        for n in range(1, 7):
            velocities = np.zeros((1 << n) - 1)
            for depth in range(n):
                velocities[2 * ((1 << depth) - 1)] = 1
            result = _differential(n, np.zeros_like(velocities), velocities)
            with self.subTest(n=n):
                np.testing.assert_allclose(
                    result.psd @ result.psd, result.psd, atol=ATOL, rtol=0
                )
                self.assertAlmostEqual(np.linalg.norm(result.psd, 2), 1, delta=ATOL)
                self.assertAlmostEqual(
                    np.linalg.norm(result.product_derivative, 2),
                    2 * math.cos(math.pi / (n + 2)), delta=ATOL,
                )

    def test_straight_path_bound_controls_complete_finite_difference(self):
        rng = np.random.default_rng(20261013)
        for n in range(1, 6):
            angles = rng.uniform(-math.pi, math.pi, (1 << n) - 1)
            original = real_frame_matrix(angles)
            for scale in (0.015, 0.25, 2.0):
                changes = rng.uniform(-scale, scale, angles.size)
                caps = np.asarray([
                    np.max(np.abs(changes[(1 << d) - 1:(1 << (d + 1)) - 1]))
                    for d in range(n)
                ])
                recurrence = _cap_bound(caps)
                rss_bound = math.sqrt(caps[0]**2 + 2 * np.sum(caps[1:]**2))
                with self.subTest(n=n, scale=scale):
                    self.assertLessEqual(recurrence, rss_bound + ATOL)
                    actual = np.linalg.norm(real_frame_matrix(angles + changes) - original, 2)
                    self.assertLessEqual(actual, min(2, recurrence) + ATOL)

    def test_isolated_angle_error_has_no_incoming_amplitude_factor(self):
        rng = np.random.default_rng(20261014)
        n = 4
        for angles in (np.zeros((1 << n) - 1),
                       rng.uniform(-math.pi, math.pi, (1 << n) - 1)):
            original = real_frame_matrix(angles)
            for coordinate in (0, 2, angles.size - 1):
                for change in (1e-5, -0.37, 2.4):
                    perturbed = angles.copy()
                    perturbed[coordinate] += change
                    difference = real_frame_matrix(perturbed) - original
                    with self.subTest(coordinate=coordinate, change=change):
                        self.assertAlmostEqual(
                            np.linalg.norm(difference, 2),
                            2 * abs(math.sin(change / 2)), delta=ATOL,
                        )
                        if np.all(angles == 0) and coordinate == 2:
                            # Node 3 has zero incoming state amplitude, but
                            # the prescribed complete frame still changes.
                            np.testing.assert_allclose(difference[:, 0], 0,
                                                       atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
