"""Bounded full-column checks for prefix-selected terminal factorization.

Native layer matrices are built from disjoint anchor/marker rotations;
the comparison frame is independently built by the library's analytic
subtree columns. These fixtures test local factorization and rounding,
not a scalable native source or the asymptotic resource allocation.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix


ATOL = 5e-12


def _layers(angles, first, stop):
    """Ordered elementary layer factors, preserving the external prefix."""
    angles = np.asarray(angles, dtype=float)
    sectors, nodes = angles.shape
    local_size = nodes + 1
    height = local_size.bit_length() - 1
    assert local_size == 1 << height
    result = np.eye(sectors * local_size)
    for depth in range(first, stop):
        layer = np.eye(len(result))
        for outer in range(sectors):
            for prefix in range(1 << depth):
                node = (1 << depth) - 1 + prefix
                left = outer * local_size + (prefix << (height - depth))
                right = left + (1 << (height - depth - 1))
                c, s = np.cos(angles[outer, node]), np.sin(angles[outer, node])
                layer[np.ix_([left, right], [left, right])] = [[c, -s], [s, c]]
        result = layer @ result
    return result


def _analytic_selected_frame(angles):
    size = len(angles[0]) + 1
    result = np.zeros((len(angles) * size, len(angles) * size))
    for prefix, values in enumerate(angles):
        selection = slice(prefix * size, (prefix + 1) * size)
        result[selection, selection] = real_frame_matrix(values).real
    return result


def _round_range(angles, first, stop, tau):
    """Finite positive fixtures; this is not an analytic sign oracle."""
    result = np.asarray(angles, dtype=float).copy()
    domain = slice((1 << first) - 1, (1 << stop) - 1)
    radix = 1 << tau
    numerators = np.clip(np.rint(radix * np.cos(result[:, domain])**2),
                         1, radix - 1)
    result[:, domain] = np.arccos(np.sqrt(numerators / radix))
    return result


class PrefixSelectedCoinTerminalTests(unittest.TestCase):
    def test_exact_selected_factorization_signed_and_singular_all_columns(self):
        height, remainder, sectors = 4, 2, 4
        angles = np.array([
            [0.0, np.pi / 2, -0.7, np.pi, -np.pi / 2, 0.31, 2.1,
             0.0, 1.2, -2.3, np.pi / 2, 0.47, -0.63, np.pi, 0.91],
            np.linspace(-2.4, 2.6, 15),
            np.linspace(0.13, 1.37, 15)[::-1],
            np.linspace(-1.1, 0.81, 15),
        ])
        early = _layers(angles, 0, height - remainder)
        terminal = _layers(angles, height - remainder, height)
        np.testing.assert_allclose(terminal @ early, _analytic_selected_frame(angles),
                                   atol=ATOL)

        # Verify the protected-sector extension independently of the layer
        # product, including every occupied protected computational input.
        expected = np.eye(sectors * (1 << height))
        for prefix in range(sectors):
            active = prefix * (1 << height) + (np.arange(4) << remainder)
            expected[np.ix_(active, active)] = real_frame_matrix(angles[prefix, :3])
        np.testing.assert_allclose(early, expected, atol=ATOL)
        occupied = [x for x in range(len(early)) if x % (1 << remainder)]
        np.testing.assert_allclose(early[:, occupied], np.eye(len(early))[:, occupied],
                                   atol=ATOL)
        # Later factors genuinely depend on nonzero external-prefix values.
        self.assertGreater(np.linalg.norm(expected[16:32, 16:32] -
                                          expected[:16, :16]), 0.1)

    def test_different_stage_denominators_preserve_full_column_target(self):
        rng = np.random.default_rng(8801)
        original = rng.uniform(0.01, np.pi / 2 - 0.01, (4, 15))
        # Singular incoming amplitudes must not suppress prescribed descendants.
        original[0, :3] = [0.0, np.pi / 2, 0.0]
        original[2, :3] = [np.pi / 2, 0.0, np.pi / 2]
        rounded = original.copy()
        stages = ((0, 2, 9), (2, 3, 11), (3, 4, 8))
        word = np.eye(64)
        errors = []
        for first, stop, tau in stages:
            changed = _round_range(original, first, stop, tau)
            start_node, stop_node = (1 << first) - 1, (1 << stop) - 1
            rounded[:, start_node:stop_node] = changed[:, start_node:stop_node]
            factor = _layers(changed, first, stop)
            exact_factor = _layers(original, first, stop)
            error = np.linalg.norm(factor - exact_factor, ord=2)
            self.assertLessEqual(error, 2 * (stop - first) * 2**(-tau / 2))
            errors.append(error)
            word = factor @ word
        np.testing.assert_allclose(word, _analytic_selected_frame(rounded), atol=ATOL)
        full_error = np.linalg.norm(word - _analytic_selected_frame(original), ord=2)
        self.assertLessEqual(full_error, sum(errors) + ATOL)
        self.assertGreater(full_error, 1e-3)

    def test_prefix_sector_errors_take_maximum_not_number_of_prefixes(self):
        original = np.array([
            [0.0, 0.71, 1.33], [0.51, np.pi / 2, 0.62],
            [1.2, 0.0, 0.37], [0.23, 1.07, 1.45],
        ])
        rounded = _round_range(original, 0, 2, 7)
        difference = _analytic_selected_frame(rounded) - _analytic_selected_frame(original)
        local_errors = [np.linalg.norm(difference[4*p:4*p+4, 4*p:4*p+4], ord=2)
                        for p in range(4)]
        self.assertAlmostEqual(np.linalg.norm(difference, ord=2), max(local_errors),
                               delta=ATOL)
        # A coherent change of basis on the preserved prefix retains the norm;
        # an arbitrary reference tensor factor retains it as well.
        hadamard = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        prefix_mix = np.kron(np.kron(hadamard, hadamard), np.eye(4))
        coherent = difference @ prefix_mix
        self.assertAlmostEqual(np.linalg.norm(np.kron(coherent, np.eye(2)), ord=2),
                               max(local_errors), delta=ATOL)


if __name__ == "__main__":
    unittest.main()
