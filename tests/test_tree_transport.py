"""Finite checks of the tree transport identities and their limitations.

The matrices have dimension at most 16. These fixtures check path products,
the history-unitary dilation, and exact finite perturbation expansions. They
do not synthesize a Clifford+T history circuit or prove an asymptotic cost.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

from tests.test_tree_residual_structure import (
    ATOL, H, S, T, _addressed_frame, _edge_product, _marker,
    _native_word, _rotation, _su2,
)


def _path_matrix(height, words):
    """Independent path products, with heap rows and physical marker columns."""
    size = 1 << height
    result = np.zeros((size - 1, size), dtype=complex)
    for node in range(1, size):
        result[node - 1, 0] = _edge_product(words, 1, node)
        for ancestor in range(1, size):
            if ancestor == node:
                continue
            coefficient = _edge_product(words, ancestor, node, True)
            if coefficient is not None:
                result[node - 1, _marker(ancestor, height)] = coefficient
    return result


def _shift_and_injection(height, words):
    size = 1 << height
    shift = np.zeros((size - 1, size - 1), dtype=complex)
    injection = np.zeros((size - 1, size), dtype=complex)
    injection[0, 0] = 1
    for parent in range(1, size // 2):
        for branch in (0, 1):
            child = 2 * parent + branch
            shift[child - 1, parent - 1] = words[parent][branch, 0]
            injection[child - 1, _marker(parent, height)] = words[parent][branch, 1]
    return shift, injection


def _power_sum(matrix, terms):
    result = np.zeros_like(matrix)
    power = np.eye(matrix.shape[0], dtype=complex)
    for _ in range(terms):
        result += power
        power = power @ matrix
    return result


def _overlap_generators(height, coarse, target):
    size = 1 << height
    overlaps = {leaf: 1.0 + 0j for leaf in range(size, 2 * size)}
    blocks = {}
    for node in range(size - 1, 0, -1):
        child_overlap = np.diag([overlaps[2 * node], overlaps[2 * node + 1]])
        blocks[node] = coarse[node].conj().T @ child_overlap @ target[node]
        overlaps[node] = blocks[node][0, 0]
    return blocks


def _embed_rows(height, matrix):
    """Place internal-node rows at their complete-frame marker positions."""
    size = 1 << height
    result = np.zeros((size, matrix.shape[1]), dtype=complex)
    for node in range(1, size):
        result[_marker(node, height)] = matrix[node - 1]
    return result


def _history_unitary(height, words):
    size = 1 << height
    result = np.eye(size, dtype=complex)
    for depth in range(height - 1):
        remaining = height - depth
        stay = 1 / np.sqrt(remaining)
        descend = np.sqrt((remaining - 1) / remaining)
        splitter = np.array([
            [stay, 0, descend],
            [descend, 0, -stay],
            [0, 1, 0],
        ], dtype=complex)
        for node in range(1 << depth, 1 << (depth + 1)):
            local_word = np.eye(3, dtype=complex)
            local_word[1:, 1:] = words[node]
            rows = [node, 2 * node, 2 * node + 1]
            result[rows] = local_word @ splitter @ result[rows]
    return result


def _perturbation_truncation(height, coarse, target, degree):
    coarse_shift, coarse_injection = _shift_and_injection(height, coarse)
    target_shift, target_injection = _shift_and_injection(height, target)
    coarse_resolvent = _power_sum(coarse_shift, height)
    perturbation = (target_shift - coarse_shift) @ coarse_resolvent
    return coarse_resolvent @ (
        _power_sum(perturbation, degree + 1) @ coarse_injection
        + _power_sum(perturbation, degree) @ (target_injection - coarse_injection)
    )


def _word_fixtures(height):
    """Native complex words, generic complex words, and singular real charts."""
    size = 1 << height
    native_choices = (
        _native_word([H, T, H]),
        _native_word([T, H, S]),
        _native_word([H, T.conj().T, H, T]),
    )
    native = {node: native_choices[(node - 1) % len(native_choices)]
              for node in range(1, size)}
    generic = {node: _su2(0.12 * node, -0.19 * node, 0.07 + 0.04 * node)
               for node in range(1, size)}
    singular_angles = (0, np.pi / 2, np.pi, -np.pi / 2)
    singular = {node: _rotation(singular_angles[(node - 1) % len(singular_angles)])
                for node in range(1, size)}
    return native, generic, singular


class TreeTransportTests(unittest.TestCase):
    def test_path_products_resolvent_and_complete_residual(self):
        for height in (2, 3, 4):
            size = 1 << height
            native, generic, singular = _word_fixtures(height)
            for coarse, target in ((native, singular), (generic, native)):
                with self.subTest(height=height, target_singular=target is singular):
                    paths = []
                    for words in (coarse, target):
                        expected = _path_matrix(height, words)
                        shift, injection = _shift_and_injection(height, words)
                        actual = _power_sum(shift, height) @ injection
                        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                        np.testing.assert_allclose(
                            np.linalg.matrix_power(shift, height), 0, atol=0, rtol=0)
                        paths.append(actual)
                    blocks = _overlap_generators(height, coarse, target)
                    h = np.array([blocks[node][1, 0] for node in range(1, size)])
                    k = np.array([blocks[node][0, 1] for node in range(1, size)])
                    diagonal = np.zeros(size, dtype=complex)
                    diagonal[0] = blocks[1][0, 0] - 1
                    for node in range(1, size):
                        diagonal[_marker(node, height)] = blocks[node][1, 1] - 1
                    residual = (np.diag(diagonal)
                                + _embed_rows(height, h[:, None] * paths[1])
                                + _embed_rows(height, k.conj()[:, None] * paths[0]).conj().T)
                    coarse_frame = _addressed_frame(height, coarse)
                    target_frame = _addressed_frame(height, target)
                    np.testing.assert_allclose(
                        residual, coarse_frame.conj().T @ target_frame - np.eye(size),
                        atol=ATOL, rtol=0)

    def test_path_gram_diagonal_including_zero_deepest_columns(self):
        for height in (2, 3, 4):
            size = 1 << height
            expected = np.zeros(size)
            expected[0] = height
            for node in range(1, size):
                expected[_marker(node, height)] = height - node.bit_length()
            for fixture, words in enumerate(_word_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    path = _path_matrix(height, words)
                    np.testing.assert_allclose(
                        path.conj().T @ path, np.diag(expected), atol=ATOL, rtol=0)
                    self.assertAlmostEqual(np.linalg.norm(path, ord=2), np.sqrt(height),
                                           delta=ATOL)
                    for node in range(size // 2, size):
                        np.testing.assert_array_equal(path[:, _marker(node, height)], 0)

    def test_subtree_defect_telescoping_and_weighted_transport_norm(self):
        for height in (2, 3, 4):
            size = 1 << height
            native, generic, singular = _word_fixtures(height)
            close_coarse = {node: _rotation(0.13 * node) for node in range(1, size)}
            close_target = {node: _rotation(0.13 * node + float(Fraction(1, 100)))
                            for node in range(1, size)}
            for fixture, (coarse, target) in enumerate((
                    (native, generic), (native, singular), (close_coarse, close_target))):
                blocks = _overlap_generators(height, coarse, target)
                defects = {node: 1 - abs(blocks[node][0, 0]) ** 2
                           for node in range(1, size)}
                gamma = max(defects.values())
                h = np.array([blocks[node][1, 0] for node in range(1, size)])
                k = np.array([blocks[node][0, 1] for node in range(1, size)])
                with self.subTest(height=height, fixture=fixture):
                    for ancestor in range(1, size):
                        target_sum = coarse_sum = 0.0
                        for node in range(1, size):
                            target_path = _edge_product(target, ancestor, node)
                            if target_path is None:
                                continue
                            coarse_path = _edge_product(coarse, ancestor, node)
                            target_sum += abs(target_path) ** 2 * abs(h[node - 1]) ** 2
                            coarse_sum += abs(coarse_path) ** 2 * abs(k[node - 1]) ** 2
                        self.assertAlmostEqual(target_sum, defects[ancestor], delta=ATOL)
                        self.assertAlmostEqual(coarse_sum, defects[ancestor], delta=ATOL)
                    for weights, words in ((h, target), (k, coarse)):
                        norm = np.linalg.norm(weights[:, None] * _path_matrix(height, words),
                                              ord=2)
                        self.assertLessEqual(norm, 2 * np.sqrt(gamma) + ATOL)
                    if fixture == 2:
                        self.assertGreater(gamma, 0)
                        self.assertLess(gamma, 0.01)

    def test_history_scatter_unitarity_and_designated_columns(self):
        for height in (2, 3, 4):
            size = 1 << height
            for fixture, words in enumerate(_word_fixtures(height)):
                with self.subTest(height=height, fixture=fixture):
                    history = _history_unitary(height, words)
                    path = _path_matrix(height, words)
                    np.testing.assert_allclose(
                        history.conj().T @ history, np.eye(size), atol=ATOL, rtol=0)
                    np.testing.assert_array_equal(history[0], np.eye(size)[0])
                    np.testing.assert_array_equal(history[:, 0], np.eye(size)[:, 0])
                    np.testing.assert_allclose(
                        history[1:, 1], path[:, 0] / np.sqrt(height), atol=ATOL, rtol=0)
                    for node in range(1, size // 2):
                        normalization = np.sqrt(height - node.bit_length())
                        np.testing.assert_allclose(
                            history[1:, 2 * node],
                            path[:, _marker(node, height)] / normalization,
                            atol=ATOL, rtol=0)

    def test_finite_perturbation_expansion_and_small_error_estimate(self):
        for height in (2, 3, 4):
            size = 1 << height
            coarse, target, _ = _word_fixtures(height)
            coarse_shift, coarse_injection = _shift_and_injection(height, coarse)
            np.testing.assert_allclose(
                _perturbation_truncation(height, coarse, target, 0),
                _power_sum(coarse_shift, height) @ coarse_injection, atol=ATOL, rtol=0)
            np.testing.assert_allclose(
                _perturbation_truncation(height, coarse, target, height - 1),
                _path_matrix(height, target), atol=ATOL, rtol=0)

            coarse = {node: _rotation(0.13 * node) for node in range(1, size)}
            target = {node: _rotation(0.13 * node + float(Fraction(1, 100)))
                      for node in range(1, size)}
            delta = max(np.linalg.norm(target[node] - coarse[node], ord=2)
                        for node in range(1, size))
            self.assertLess(height * delta, 1)
            blocks = _overlap_generators(height, coarse, target)
            h = np.array([blocks[node][1, 0] for node in range(1, size)])
            self.assertLessEqual(max(abs(h)), height * delta + ATOL)
            path = _path_matrix(height, target)
            for degree in range(height - 1):
                actual_error = np.linalg.norm(h[:, None] * (
                    path - _perturbation_truncation(height, coarse, target, degree)), ord=2)
                bound = ((height + 1) * (height * delta) ** (degree + 2)
                         / (1 - height * delta))
                with self.subTest(height=height, degree=degree):
                    self.assertLessEqual(actual_error, bound + ATOL)

    def test_right_spine_omissions_at_two_three_and_four_levels(self):
        angle = float(Fraction(1, 16))
        for height in (2, 3, 4):
            size = 1 << height
            for degree in range(height - 1):
                coarse = {node: np.eye(2, dtype=complex) for node in range(1, size)}
                target = {node: np.eye(2, dtype=complex) for node in range(1, size)}
                node = 1
                for depth in range(degree + 2):
                    target[node] = _rotation(angle)
                    if depth < degree + 1:
                        node = 2 * node + 1
                blocks = _overlap_generators(height, coarse, target)
                truncated = _perturbation_truncation(height, coarse, target, degree)
                expected = np.sin(angle) ** (degree + 2)
                actual_residual = (_addressed_frame(height, coarse).conj().T
                                   @ _addressed_frame(height, target) - np.eye(size))
                with self.subTest(height=height, degree=degree):
                    self.assertEqual(truncated[node - 1, 0], 0)
                    self.assertAlmostEqual(actual_residual[_marker(node, height), 0],
                                           expected, delta=ATOL)
                    self.assertAlmostEqual(
                        blocks[node][1, 0] * _path_matrix(height, target)[node - 1, 0],
                        expected, delta=ATOL)
                    self.assertGreater(expected, 100 * ATOL)


if __name__ == "__main__":
    unittest.main()
