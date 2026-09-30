"""Small independent checks of the classical tree-residual factorization.

These fixtures compare addressed circuits, recursively constructed vectors,
and scalar overlap generators. They do not implement a coherent product
evaluator or establish any asymptotic quantum gate-count improvement.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 3e-12
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
T = np.diag([1, np.exp(1j * np.pi / 4)])
S = np.diag([1, 1j])


def _rotation(angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array([[cosine, -sine], [sine, cosine]], dtype=complex)


def _su2(angle, left_phase, right_phase):
    left = np.diag([np.exp(1j * left_phase), np.exp(-1j * left_phase)])
    right = np.diag([np.exp(1j * right_phase), np.exp(-1j * right_phase)])
    return left @ _rotation(angle) @ right


def _native_word(gates):
    unitary = np.eye(2, dtype=complex)
    for gate in gates:
        unitary = gate @ unitary
    return unitary


def _marker(node, height):
    depth = node.bit_length() - 1
    prefix = node - (1 << depth)
    return (2 * prefix + 1) << (height - depth - 1)


def _addressed_frame(height, words):
    """Apply literal addressed two-row gates in shallow-to-deep order."""
    frame = np.eye(1 << height, dtype=complex)
    for depth in range(height):
        stride = 1 << (height - depth - 1)
        for prefix in range(1 << depth):
            rows = [2 * prefix * stride, (2 * prefix + 1) * stride]
            frame[rows] = words[(1 << depth) + prefix] @ frame[rows]
    return frame


def _recursive_frame(height, words):
    """Build subtree vectors independently of the addressed circuit."""
    size = 1 << height
    basis = np.eye(size, dtype=complex)
    states = {size + leaf: basis[:, leaf] for leaf in range(size)}
    complements = {}
    for node in range(size - 1, 0, -1):
        children = np.column_stack([states[2 * node], states[2 * node + 1]])
        pair = children @ words[node]
        states[node], complements[node] = pair[:, 0], pair[:, 1]
    frame = np.zeros((size, size), dtype=complex)
    frame[:, 0] = states[1]
    for node, vector in complements.items():
        frame[:, _marker(node, height)] = vector
    return frame


def _edge_product(words, ancestor, descendant, complement=False):
    """Return None off the ancestor chain; an empty path has value one."""
    route = []
    current = descendant
    while current > ancestor:
        route.append(current & 1)
        current //= 2
    if current != ancestor:
        return None
    coefficient = 1.0 + 0j
    current = ancestor
    for step, branch in enumerate(reversed(route)):
        column = int(complement and step == 0)
        coefficient *= words[current][branch, column]
        current = 2 * current + branch
    return coefficient


def _generator_residual(height, coarse, target):
    """Scalar recurrences and path products; no subtree vectors are used."""
    size = 1 << height
    overlaps = {leaf: 1.0 + 0j for leaf in range(size, 2 * size)}
    generators = {}
    for node in range(size - 1, 0, -1):
        child = np.diag([overlaps[2 * node], overlaps[2 * node + 1]])
        block = coarse[node].conj().T @ child @ target[node]
        generators[node] = block
        overlaps[node] = block[0, 0]
    result = np.zeros((size, size), dtype=complex)
    result[0, 0] = overlaps[1] - 1
    for node, block in generators.items():
        marker = _marker(node, height)
        h, k, diagonal = block[1, 0], block[0, 1], block[1, 1]
        result[marker, marker] = diagonal - 1
        result[marker, 0] = _edge_product(target, 1, node) * h
        result[0, marker] = np.conj(_edge_product(coarse, 1, node)) * k
        for ancestor in range(1, size):
            if ancestor == node:
                continue
            target_path = _edge_product(target, ancestor, node, True)
            if target_path is None:
                continue
            ancestor_marker = _marker(ancestor, height)
            coarse_path = _edge_product(coarse, ancestor, node, True)
            result[marker, ancestor_marker] = target_path * h
            result[ancestor_marker, marker] = np.conj(coarse_path) * k
    return result


class TreeResidualStructureTests(unittest.TestCase):
    def _assert_factorization(self, height, coarse, target):
        size = 1 << height
        c_frame = _addressed_frame(height, coarse)
        w_frame = _addressed_frame(height, target)
        for words, frame in ((coarse, c_frame), (target, w_frame)):
            np.testing.assert_allclose(
                _recursive_frame(height, words), frame, atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                frame.conj().T @ frame, np.eye(size), atol=ATOL, rtol=0)
        expected = c_frame.conj().T @ w_frame - np.eye(size)
        actual = _generator_residual(height, coarse, target)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        return expected

    def test_all_entries_with_actual_complex_native_coarse_words(self):
        choices = [
            _native_word([H, T, H]),
            _native_word([T, H, S]),
            _native_word([H, T.conj().T, H, T]),
        ]
        for height in (2, 3):
            size = 1 << height
            coarse = {node: choices[(node - 1) % len(choices)]
                      for node in range(1, size)}
            target = {node: _rotation(0.17 + 0.21 * node)
                      for node in range(1, size)}
            with self.subTest(height=height):
                self._assert_factorization(height, coarse, target)

    def test_all_entries_with_general_complex_su2_coarse_words(self):
        for height in (2, 3):
            size = 1 << height
            coarse = {node: _su2(0.12 * node, -0.19 * node, 0.07 + 0.04 * node)
                      for node in range(1, size)}
            target = {node: _rotation(-0.32 + 0.16 * node)
                      for node in range(1, size)}
            with self.subTest(height=height):
                self._assert_factorization(height, coarse, target)

    def test_singular_charts_need_no_division_and_keep_marker_columns(self):
        height = 3
        size = 1 << height
        coarse = {node: _native_word([H, T if node & 1 else S, H])
                  for node in range(1, size)}
        singular_angles = [0, np.pi / 2, 0, np.pi, -np.pi / 2, 0, np.pi / 2]
        target = {node: _rotation(singular_angles[node - 1])
                  for node in range(1, size)}
        self._assert_factorization(height, coarse, target)

        # Root angle zero makes the right subtree absent from column zero,
        # but its independently prescribed marker columns still vary.
        first = {1: _rotation(0), 2: _rotation(0.37), 3: _rotation(0.21)}
        second = dict(first)
        second[3] = _rotation(1.13)
        first_frame = _addressed_frame(2, first)
        second_frame = _addressed_frame(2, second)
        np.testing.assert_allclose(first_frame[:, 0], second_frame[:, 0], atol=0)
        self.assertGreater(np.linalg.norm(first_frame - second_frame, ord=2), 0.5)
        self._assert_factorization(2, first, second)

    def test_two_level_target_transport_witness_and_omitted_error(self):
        alpha, gamma, delta = 0.62, 0.55, 0.03
        coarse = {1: _rotation(gamma), 2: _rotation(0.21), 3: _rotation(-0.43)}
        target = {1: _rotation(alpha), 2: _rotation(0.21 + delta),
                  3: _rotation(0.09)}
        residual = self._assert_factorization(2, coarse, target)
        expected = np.sin(delta) * np.array([
            np.cos(alpha), -np.sin(alpha), -np.cos(gamma), np.sin(gamma)])
        actual = residual[[1, 1, 0, 2], [0, 2, 1, 1]]
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        coarse_transport = np.sin(delta) * np.array([np.cos(gamma), -np.sin(gamma)])
        omitted = np.linalg.norm(actual[:2] - coarse_transport)
        self.assertAlmostEqual(
            omitted, 2 * abs(np.sin(delta) * np.sin((alpha - gamma) / 2)),
            delta=ATOL)
        self.assertGreater(omitted, 1e-3)

    def test_three_level_target_and_coarse_path_products(self):
        alpha0, alpha1, gamma0, gamma1, delta = 0.63, -0.27, 0.44, -0.35, 0.06
        coarse = {node: _rotation(0.08 * node) for node in range(1, 8)}
        target = {node: _rotation(-0.04 + 0.11 * node) for node in range(1, 8)}
        coarse.update({1: _rotation(gamma0), 2: _rotation(gamma1),
                       4: _rotation(0.31)})
        target.update({1: _rotation(alpha0), 2: _rotation(alpha1),
                       4: _rotation(0.31 + delta)})
        residual = self._assert_factorization(3, coarse, target)
        expected = np.sin(delta) * np.array([
            np.cos(alpha0) * np.cos(alpha1), -np.sin(alpha1),
            -np.sin(alpha0) * np.cos(alpha1),
            -np.cos(gamma0) * np.cos(gamma1), np.sin(gamma1),
            np.sin(gamma0) * np.cos(gamma1)])
        actual = residual[[1, 1, 1, 0, 2, 4], [0, 2, 4, 1, 1, 1]]
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
