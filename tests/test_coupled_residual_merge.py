"""Bounded whole-port checks for the coupled residual boundary contract.

Logical matrices have dimension at most eight and signal/logical matrices
at most sixteen. These are algebraic checks, not native circuit emitters or
evidence of a better asymptotic compiler cost. In particular, the recursive
wrappers still contain the target local words.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_tree_residual_structure import (
    ATOL, H, S, T, _addressed_frame, _rotation, _su2,
)


SINE = np.sqrt(3) / 2


def _direct_sum(left, right):
    rows, columns = left.shape[0], right.shape[0]
    result = np.zeros((rows + columns, rows + columns), dtype=complex)
    result[:rows, :rows] = left
    result[rows:, rows:] = right
    return result


def _canonical(residual):
    """One signal flag; the displayed four blocks include rejected inputs."""
    identity = np.eye(len(residual), dtype=complex)
    return np.block([
        [residual / 2, -SINE * identity],
        [SINE * identity, residual.conj().T / 2],
    ])


def _root_word(local, size):
    """Embed on existing child-root modes 0 and size/2; add no modes."""
    result = np.eye(size, dtype=complex)
    rows = [0, size // 2]
    result[np.ix_(rows, rows)] = local
    return result


def _child_signals(left, right):
    """Direct sum child words, ordered by one shared signal then logical bits."""
    child_size = len(left) // 2
    size = 2 * child_size
    result = np.zeros((2 * size, 2 * size), dtype=complex)
    for branch, child in enumerate((left, right)):
        rows = [flag * size + branch * child_size + index
                for flag in (0, 1) for index in range(child_size)]
        result[np.ix_(rows, rows)] = child
    return result


def _merge(coarse_root, discrepancy, child_word):
    """Literal inverse wrappers; no accepted-block products are substituted."""
    size = len(coarse_root)
    identity = np.eye(size, dtype=complex)
    coarse_both = _direct_sum(coarse_root, coarse_root)
    left = _direct_sum(identity, discrepancy.conj().T)
    right = _direct_sum(discrepancy, identity)
    return left @ coarse_both.conj().T @ child_word @ coarse_both @ right


def _recursive_word(height, coarse, target, node=1):
    if height == 0:
        return _canonical(np.ones((1, 1), dtype=complex))
    left = _recursive_word(height - 1, coarse, target, 2 * node)
    right = _recursive_word(height - 1, coarse, target, 2 * node + 1)
    size = 1 << height
    coarse_root = _root_word(coarse[node], size)
    discrepancy = coarse_root.conj().T @ _root_word(target[node], size)
    return _merge(coarse_root, discrepancy, _child_signals(left, right))


def _native_small_word():
    """A finite Clifford+T commutator word, including literal scalar phases.

    Each update is K (S K S†) K† (S K† S†). The determinant-one
    commutator and its four updates are exact native words; no polar
    projection, phase normalization, or gate synthesis is performed here.
    """
    other = H @ T @ H
    word = T @ other @ T.conj().T @ other.conj().T
    for _ in range(4):
        conjugate = S @ word @ S.conj().T
        word = word @ conjugate @ word.conj().T @ conjugate.conj().T
    return word


def _close_real_fixture(height):
    """Independently changed real nodes and an actual complex native baseline."""
    native = _native_small_word()
    choices = (native, S @ native @ S.conj().T, H @ native.conj().T @ H)
    coarse = {node: choices[(node - 1) % len(choices)]
              for node in range(1, 1 << height)}
    target = {node: _rotation((2 + node) / (4000 * height))
              for node in range(1, 1 << height)}
    return coarse, target


def _subtree_residual(height, coarse, target, node):
    if height == 0:
        return np.ones((1, 1), dtype=complex)
    children = _direct_sum(
        _subtree_residual(height - 1, coarse, target, 2 * node),
        _subtree_residual(height - 1, coarse, target, 2 * node + 1),
    )
    size = 1 << height
    return (_root_word(coarse[node], size).conj().T @ children
            @ _root_word(target[node], size))


class CoupledResidualMergeTests(unittest.TestCase):
    def test_canonical_completion_has_fixed_rejected_blocks_and_exact_stability(self):
        first = _su2(0.23, -0.17, 0.08)
        second = _su2(-0.19, 0.31, -0.07)
        for residual in (first, second, np.eye(2, dtype=complex)):
            word = _canonical(residual)
            np.testing.assert_allclose(word.conj().T @ word, np.eye(4),
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(word[:2, 2:], -SINE * np.eye(2),
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(word[2:, :2], SINE * np.eye(2),
                                       atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(_canonical(first) - _canonical(second), 2),
                               np.linalg.norm(first - second, 2) / 2, places=13)

    def test_complex_native_fork_and_height_three_match_every_residual_port(self):
        for height in (2, 3):
            with self.subTest(height=height):
                coarse, target = _close_real_fixture(height)
                for node in (1, 2, 3):
                    self.assertGreater(np.max(np.abs(coarse[node].imag)), 1e-6)
                    self.assertGreater(np.linalg.norm(coarse[node] - target[node], 2), 1e-4)
                    np.testing.assert_allclose(np.linalg.det(coarse[node]), 1,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_array_equal(target[node].imag, 0)
                native_frame = _addressed_frame(height, coarse)
                target_frame = _addressed_frame(height, target)
                residual = native_frame.conj().T @ target_frame
                self.assertLess(np.linalg.norm(native_frame - target_frame, 2), 1 / 128)
                merged = _recursive_word(height, coarse, target)
                np.testing.assert_allclose(merged, _canonical(residual), atol=ATOL, rtol=0)
                np.testing.assert_allclose(merged.conj().T @ merged, np.eye(len(merged)),
                                           atol=ATOL, rtol=0)

    def test_unfolded_word_exposes_target_and_preserves_arbitrary_spectators(self):
        height = 3
        coarse, target = _close_real_fixture(height)
        native_frame = _addressed_frame(height, coarse)
        target_frame = _addressed_frame(height, target)
        size = 1 << height
        unfolded = (_direct_sum(native_frame.conj().T, target_frame.conj().T)
                    @ _canonical(np.eye(size))
                    @ _direct_sum(target_frame, native_frame))
        merged = _recursive_word(height, coarse, target)
        np.testing.assert_allclose(merged, unfolded, atol=ATOL, rtol=0)

        # Eight spectator basis states represent an idle second external
        # flag, one dirty qubit, and a reference. No spectator starts at zero.
        rng = np.random.default_rng(401)
        state = rng.normal(size=(2 * size, 8)) + 1j * rng.normal(size=(2 * size, 8))
        state /= np.linalg.norm(state)
        np.testing.assert_allclose(merged @ state, unfolded @ state, atol=ATOL, rtol=0)
        np.testing.assert_allclose(merged.conj().T @ (merged @ state), state,
                                   atol=ATOL, rtol=0)

    def test_actual_native_substitutions_keep_offdiagonal_cancellation(self):
        height = 3
        coarse, target = _close_real_fixture(height)
        base = _native_small_word()
        # A different literal native target word at every node. This is a
        # finite substitution test, not a synthesis-accuracy claim.
        substituted = {node: coarse[node] @ (base if node % 2 else base.conj().T)
                       for node in range(1, 1 << height)}
        native_frame = _addressed_frame(height, coarse)
        substituted_frame = _addressed_frame(height, substituted)
        target_frame = _addressed_frame(height, target)
        actual = _recursive_word(height, coarse, substituted)
        size = 1 << height
        np.testing.assert_allclose(actual[:size, size:], -SINE * np.eye(size),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual[size:, :size], SINE * np.eye(size),
                                   atol=ATOL, rtol=0)
        expected = _canonical(native_frame.conj().T @ substituted_frame)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        desired = _canonical(native_frame.conj().T @ target_frame)
        self.assertAlmostEqual(np.linalg.norm(actual - desired, 2),
                               np.linalg.norm(substituted_frame - target_frame, 2) / 2,
                               places=11)

    def test_merge_error_uses_max_child_error_and_matched_inverse_wrapper(self):
        coarse = _root_word(_su2(0.13, 0.19, -0.11), 4)
        discrepancy = _root_word(_rotation(0.21), 4)
        changed = _root_word(_rotation(0.2107), 4)
        children = [_canonical(_su2(0.17, 0.04, -0.09)),
                    _canonical(_su2(-0.12, 0.08, 0.15))]
        # Perturb complete child words, including their rejected blocks.
        noisy = [np.diag(np.exp(1j * np.array([0.0002, -0.0001, 0.0003, 0]))) @ children[0],
                 np.diag(np.exp(1j * np.array([0.0001, 0.0004, -0.0002, 0]))) @ children[1]]
        child_error = max(np.linalg.norm(noisy[index] - children[index], 2)
                          for index in (0, 1))
        actual = _merge(coarse, changed, _child_signals(*noisy))
        desired = _merge(coarse, discrepancy, _child_signals(*children))
        bound = child_error + np.linalg.norm(changed - discrepancy, 2) / 2
        self.assertLessEqual(np.linalg.norm(actual - desired, 2), bound + ATOL)

    def test_anchored_product_retains_parent_children_mixed_term(self):
        coarse, target = _close_real_fixture(2)
        coarse_root = _root_word(coarse[1], 4)
        descendant = _direct_sum(
            _subtree_residual(1, coarse, target, 2),
            _subtree_residual(1, coarse, target, 3),
        )
        first = coarse_root.conj().T @ descendant @ coarse_root
        second = coarse_root.conj().T @ _root_word(target[1], 4)
        identity = np.eye(4)
        anchored = _canonical(first) @ _canonical(identity).conj().T @ _canonical(second)
        defect = anchored[:4, :4] - first @ second / 2
        mixed = -3 / 8 * (first - identity) @ (second - identity)
        np.testing.assert_allclose(defect, mixed, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(mixed, 2), 1e-8)

        # The parent discrepancy acts on exactly two existing root modes.
        # Full rejected-space error also admits a small-rank factorization,
        # whose nontrivial vectors still require charged tree transport.
        first_change, second_change = first - identity, second - identity
        self.assertEqual(np.linalg.matrix_rank(second_change, tol=ATOL), 2)
        self.assertLessEqual(np.linalg.matrix_rank(mixed, tol=ATOL), 2)
        full_error = anchored - _canonical(first @ second)
        factorized = (np.vstack([SINE * first_change, first_change.conj().T / 2])
                      @ np.hstack([-SINE * second_change, second_change.conj().T / 2]) / 2)
        factorized[4:, 4:] -= second_change.conj().T @ first_change.conj().T / 2
        np.testing.assert_allclose(full_error, factorized, atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.matrix_rank(full_error, tol=ATOL), 4)
        repair = _canonical(first @ second) @ anchored.conj().T
        np.testing.assert_allclose(repair @ anchored, _canonical(first @ second),
                                   atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.matrix_rank(repair - np.eye(8), tol=ATOL), 4)
        transported = []
        for root in (0, 2):
            transported.append(np.concatenate([
                SINE * first_change[:, root], first_change.conj().T[:, root] / 2]))
            transported.append(np.concatenate([np.zeros(4), identity[:, root]]))
        moving, _ = np.linalg.qr(np.column_stack(transported))
        complement = np.eye(8) - moving @ moving.conj().T
        np.testing.assert_allclose(complement @ (repair - np.eye(8)), 0,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose((repair - np.eye(8)) @ complement, 0,
                                   atol=ATOL, rtol=0)

        # Disjoint updates are a distinct special case: their mixed term is zero.
        disjoint_first = _direct_sum(_rotation(0.14), np.eye(2))
        disjoint_second = _direct_sum(np.eye(2), _rotation(-0.27))
        disjoint = (_canonical(disjoint_first) @ _canonical(identity).conj().T
                    @ _canonical(disjoint_second))
        np.testing.assert_allclose(disjoint, _canonical(disjoint_first @ disjoint_second),
                                   atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
