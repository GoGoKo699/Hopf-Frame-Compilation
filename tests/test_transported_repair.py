"""Full-input checks of the transported rank-four commutator repair.

The candidate uses an already supplied complete child word, its actual
inverse, the local parent word and its inverse.  D(AB) is used only as an
independent reference.  No test treats a transported basis, a child call,
or a fine-precision local parent word as a free native operation.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_coupled_residual_merge import (
    SINE, _canonical, _close_real_fixture, _direct_sum, _root_word,
    _subtree_residual,
)
from tests.test_tree_residual_structure import ATOL, _rotation, _su2
from tests.test_operator_source_compiler import _word_matrix


def _repair_word(child_word, parent):
    """One signal; arbitrary whole child matrix, with its literal inverse."""
    size = len(parent)
    fixed = _canonical(np.eye(size))
    local = _direct_sum(np.eye(size), parent.conj().T)
    transport = child_word @ fixed.conj().T
    repair = local @ transport @ local.conj().T @ transport.conj().T
    anchored = transport @ _canonical(parent)
    return repair, anchored, transport, local


def _fixture(height):
    coarse, target = _close_real_fixture(height)
    children = _direct_sum(
        _subtree_residual(height - 1, coarse, target, 2),
        _subtree_residual(height - 1, coarse, target, 3),
    )
    coarse_root = _root_word(coarse[1], 1 << height)
    target_root = _root_word(target[1], 1 << height)
    return coarse_root.conj().T @ children @ coarse_root, coarse_root.conj().T @ target_root, coarse


def _moving_columns(transport, logical_size):
    """A support certificate only; this SVD is not a proposed circuit."""
    roots = (0, logical_size // 2)
    projector_columns = np.eye(2 * logical_size)[:, [logical_size + r for r in roots]]
    columns = np.concatenate((projector_columns, transport @ projector_columns), axis=1)
    left, singular, _ = np.linalg.svd(columns, full_matrices=False)
    rank = int(np.count_nonzero(singular > 1e-11))
    return projector_columns, left[:, :rank], rank


def _unitary_noise(size, angle):
    """A noncanonical signal/core perturbation, without a matrix exponential."""
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    z = np.diag([1, -1])
    generator = np.kron(x, np.kron(np.eye(size), z))
    return np.cos(angle) * np.eye(4 * size) + 1j * np.sin(angle) * generator


class TransportedRepairTests(unittest.TestCase):
    def test_complete_native_coarse_fork_and_height_three_without_product_oracle(self):
        for height in (2, 3):
            with self.subTest(height=height):
                a, b, coarse = _fixture(height)
                for node in (1, 2, 3):
                    self.assertGreater(np.max(np.abs(coarse[node].imag)), 1e-6)
                child = _canonical(a)
                repair, anchored, transport, local = _repair_word(child, b)
                identity = np.eye(len(repair))
                np.testing.assert_allclose(repair.conj().T @ repair, identity,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(repair @ anchored, _canonical(a @ b),
                                           atol=ATOL, rtol=0)
                # This is the emitted cancellation, containing only one child call.
                one_child = local @ child @ _direct_sum(b, np.eye(len(b)))
                np.testing.assert_allclose(repair @ anchored, one_child,
                                           atol=ATOL, rtol=0)
                self.assertGreater(np.linalg.norm(anchored - _canonical(a @ b), 2), 1e-8)
                np.testing.assert_allclose(
                    transport.conj().T @ transport, identity, atol=ATOL, rtol=0)

    def test_transported_modes_and_support_include_zero_and_dependent_cases(self):
        size = 4
        b = _root_word(_su2(.21, -.17, .09), size)
        generic, _, _ = _fixture(2)
        one_root = np.diag(np.exp(1j * np.array([.3, -.3, 0, 0])))
        off_roots = np.diag(np.exp(1j * np.array([0, .2, 0, -.2])))
        cases = ((np.eye(size), b, 2), (off_roots, b, 2),
                 (one_root, b, 3), (generic, b, 4), (generic, np.eye(size), 4))
        for a, parent, expected_rank in cases:
            with self.subTest(expected_rank=expected_rank, parent_identity=np.allclose(parent, np.eye(size))):
                repair, _, transport, _ = _repair_word(_canonical(a), parent)
                root_columns, basis, rank = _moving_columns(transport, size)
                self.assertEqual(rank, expected_rank)
                support = basis @ basis.conj().T
                defect = repair - np.eye(2 * size)
                np.testing.assert_allclose((np.eye(2 * size) - support) @ defect,
                                           np.zeros_like(defect), atol=ATOL, rtol=0)
                np.testing.assert_allclose(defect @ (np.eye(2 * size) - support),
                                           np.zeros_like(defect), atol=ATOL, rtol=0)
                roots = np.eye(size)[:, [0, size // 2]]
                difference_modes = np.concatenate((SINE * (a - np.eye(size)) @ roots,
                                                   .5 * (a.conj().T - np.eye(size)) @ roots), axis=0)
                np.testing.assert_allclose(transport @ root_columns - root_columns,
                                           .5 * difference_modes, atol=ATOL, rtol=0)
                if expected_rank == 2 or np.allclose(parent, np.eye(size)):
                    np.testing.assert_allclose(repair, np.eye(2 * size), atol=ATOL, rtol=0)

    def test_mixed_defect_and_weighted_stability_need_no_nonzero_mode_gap(self):
        size = 4
        a_base, _, _ = _fixture(2)
        for scale in (1, 1e-3, 1e-8, 0):
            a = _root_word(_rotation(.37 * scale), size)
            b = _root_word(_su2(.21 * scale, -.11 * scale, .07 * scale), size)
            repair, _, transport, _ = _repair_word(_canonical(a), b)
            bound = np.linalg.norm(a - np.eye(size), 2) * np.linalg.norm(b - np.eye(size), 2)
            self.assertLessEqual(np.linalg.norm(repair - np.eye(2 * size), 2), bound + ATOL)
            self.assertAlmostEqual(np.linalg.norm(transport - np.eye(2 * size), 2),
                                   .5 * np.linalg.norm(a - np.eye(size), 2), delta=ATOL)
            perturbed = _root_word(_rotation(.013), size) @ a_base
            first = _repair_word(_canonical(a_base), b)[0]
            second = _repair_word(_canonical(perturbed), b)[0]
            weighted = np.linalg.norm(b - np.eye(size), 2) * np.linalg.norm(perturbed - a_base, 2)
            self.assertLessEqual(np.linalg.norm(second - first, 2), weighted + ATOL)
            perturbed_b = _root_word(_su2(-.017, .011, .009), size) @ b
            both = _repair_word(_canonical(perturbed), perturbed_b)[0]
            joint_bound = (weighted + np.linalg.norm(perturbed - np.eye(size), 2)
                           * np.linalg.norm(perturbed_b - b, 2))
            self.assertLessEqual(np.linalg.norm(both - first, 2), joint_bound + ATOL)

    def test_noncanonical_child_and_shared_dirty_core_cancel_with_actual_inverse(self):
        a, b, _ = _fixture(2)
        size, dirty = len(a), 2
        ideal_child = np.kron(_canonical(a), np.eye(dirty))
        noise = _unitary_noise(size, .019)
        actual_child = noise @ ideal_child
        parent = np.kron(b, np.eye(dirty))
        # The matrix order remains signal, logical, dirty, with no new clean core.
        repair, anchored, _, local = _repair_word(actual_child, parent)
        emitted = local @ actual_child @ _direct_sum(parent, np.eye(size * dirty))
        np.testing.assert_allclose(repair @ anchored, emitted, atol=ATOL, rtol=0)
        reference = np.kron(_canonical(a @ b), np.eye(dirty))
        delta = np.linalg.norm(actual_child - ideal_child, 2)
        self.assertAlmostEqual(np.linalg.norm(emitted - reference, 2), delta, delta=ATOL)
        # It is a true noncanonical perturbation, so exact cancellation cannot
        # rely on canonical off-diagonal blocks of the approximate child.
        block_size = size * dirty
        self.assertGreater(np.linalg.norm(actual_child[:block_size, block_size:]
                                          + SINE * np.eye(block_size), 2), .001)
        # A maximally correlated input includes every logical/dirty column and a
        # reference; matrix equality already proves the arbitrary-reference case.
        dimension = len(emitted)
        correlated = np.eye(dimension).reshape(-1) / np.sqrt(dimension)
        difference = np.kron(repair @ anchored - emitted, np.eye(dimension)) @ correlated
        self.assertLess(np.linalg.norm(difference), ATOL)

    def test_paired_local_approximation_has_half_operator_error_on_all_ports(self):
        a, b, _ = _fixture(2)
        approximate_b = _root_word(_rotation(.023), len(b)) @ b
        repair, anchored, _, _ = _repair_word(_canonical(a), approximate_b)
        actual = repair @ anchored
        np.testing.assert_allclose(actual, _canonical(a @ approximate_b), atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(actual - _canonical(a @ b), 2),
                               .5 * np.linalg.norm(approximate_b - b, 2), delta=ATOL)

    def test_substituting_an_ideal_inverse_breaks_the_shared_core_cancellation(self):
        a, _, _ = _fixture(2)
        size = len(a)
        ideal = np.kron(_canonical(a), np.eye(2))
        actual = _unitary_noise(size, .019) @ ideal
        identity_parent = np.eye(2 * size)
        repair, anchored, transport, _ = _repair_word(actual, identity_parent)
        np.testing.assert_allclose(repair, np.eye(len(actual)), atol=ATOL, rtol=0)
        np.testing.assert_allclose(anchored, actual, atol=ATOL, rtol=0)
        ideal_transport = ideal @ _canonical(identity_parent).conj().T
        incorrect_repair = transport @ ideal_transport.conj().T
        defect = np.linalg.norm(incorrect_repair @ anchored - actual, 2)
        self.assertGreater(defect, .01)
        self.assertAlmostEqual(defect, np.linalg.norm(actual - ideal, 2), delta=ATOL)

    def test_one_actual_native_connector_and_its_inverse_preserve_cancellation(self):
        a, b, _ = _fixture(2)
        size = len(a)
        # This literal word tests the identity for an arbitrary native connector;
        # it is not presented as a fine approximation to the desired signal angle.
        signal = _word_matrix(1, [('T', 0), ('H', 0), ('TDG', 0), ('H', 0)])
        connector = np.kron(signal, np.eye(size))
        child = _direct_sum(a, np.eye(size)) @ connector @ _direct_sum(np.eye(size), a.conj().T)
        parent = _direct_sum(b, np.eye(size)) @ connector @ _direct_sum(np.eye(size), b.conj().T)
        local = _direct_sum(np.eye(size), b.conj().T)
        transport = child @ connector.conj().T
        repair = local @ transport @ local.conj().T @ transport.conj().T
        actual = repair @ transport @ parent
        reference = (_direct_sum(a @ b, np.eye(size)) @ connector
                     @ _direct_sum(np.eye(size), (a @ b).conj().T))
        np.testing.assert_allclose(actual, reference, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(actual - _canonical(a @ b), 2),
                               np.linalg.norm(connector - _canonical(np.eye(size)), 2), delta=ATOL)

    def test_matched_noisy_connector_child_and_wrappers_cancel_on_the_same_dirty_core(self):
        size = 2
        a, b = _su2(.13, -.19, .07), _su2(-.11, .17, .09)
        ideal_child = np.kron(_canonical(a), np.eye(2))
        ideal_connector = np.kron(_canonical(np.eye(size)), np.eye(2))
        parent = np.kron(b, np.eye(2))
        ideal_left = _direct_sum(np.eye(2 * size), parent.conj().T)
        ideal_right = _direct_sum(parent, np.eye(2 * size))
        child = np.exp(.017j) * _unitary_noise(size, .031) @ ideal_child
        connector = np.exp(.021j) * _unitary_noise(size, -.023) @ ideal_connector
        left = np.exp(-.034j) * _unitary_noise(size, .013) @ ideal_left
        right = np.exp(.025j) * _unitary_noise(size, -.009) @ ideal_right
        # The local completion must be this actual word. Substituting an
        # independently rounded canonical(parent) would not provide this identity.
        local_completion = left @ connector @ right
        transport = child @ connector.conj().T
        repair = left @ transport @ left.conj().T @ transport.conj().T
        anchored = transport @ local_completion
        actual = repair @ anchored
        emitted = left @ child @ right
        np.testing.assert_allclose(actual, emitted, atol=ATOL, rtol=0)
        reference = ideal_left @ ideal_child @ ideal_right
        upper = (np.linalg.norm(left - ideal_left, 2)
                 + np.linalg.norm(child - ideal_child, 2)
                 + np.linalg.norm(right - ideal_right, 2))
        self.assertLessEqual(np.linalg.norm(actual - reference, 2), upper + ATOL)
        # The connector cancels completely even though it has phase, signal,
        # and dirty-core errors. No clean return or reset is inserted between calls.
        alternate = np.exp(-.081j) * _unitary_noise(size, .047) @ ideal_connector
        alternate_transport = child @ alternate.conj().T
        alternate_repair = left @ alternate_transport @ left.conj().T @ alternate_transport.conj().T
        alternate_anchored = alternate_transport @ left @ alternate @ right
        np.testing.assert_allclose(alternate_repair @ alternate_anchored, emitted,
                                   atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
