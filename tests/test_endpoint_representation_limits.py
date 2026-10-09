"""Finite checks for the archived fixed-bank and Spin representation results.

Integer identities below are exact. Dense floating-point checks retain every
logical and dirty column; ideal Spin rotations are not native compilers.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.conventions import anchor_label, marker_label
from compiler_robust_hopf.frames import direct_real_frame
from tests.test_source_t_depth import _majoranas


ATOL = 3e-12


def _commutator(a, b):
    return a @ b - b @ a


def _bank_polynomial(a, b):
    commutator = _commutator(a, b)
    return _commutator(commutator @ commutator, a)


def _binary_rank(labels):
    pivots = {}
    for label in labels:
        while label:
            bit = label.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = label
                break
            label ^= pivots[bit]
    return len(pivots)


class EndpointRepresentationLimitsTests(unittest.TestCase):
    def test_fixed_two_by_two_bank_identity_and_exact_frame_witness(self):
        # General integer matrices check the polynomial beyond unitary pairs.
        a = np.zeros((6, 6), dtype=np.int64)
        b = np.zeros_like(a)
        for block, (left, right) in enumerate((
                ([[0, 1], [1, 0]], [[1, 0], [0, -1]]),
                ([[0, -1], [1, 0]], [[0, 1], [1, 0]]),
                ([[2, 3], [-1, 4]], [[1, -2], [5, 0]]))):
            sl = slice(2 * block, 2 * block + 2)
            a[sl, sl], b[sl, sl] = left, right
        np.testing.assert_array_equal(_bank_polynomial(a, b), 0)

        u = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        v = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        witness = _bank_polynomial(u, v)
        np.testing.assert_array_equal(witness,
                                      [[0, -2, 0], [-2, 0, -2], [0, 2, 0]])
        gram = witness.T @ witness
        # Nonzero Gram eigenvalues are exactly 8, hence norm = 2 sqrt(2).
        np.testing.assert_array_equal(gram @ gram, 8 * gram)
        self.assertEqual(np.trace(gram), 16)
        # These matrices are restrictions of two legal, overlapping tree edges.
        for node, expected in ((1, u), (2, v)):
            angles = np.zeros(7)
            angles[node - 1] = np.pi / 2
            actual = direct_real_frame(3, angles)[np.ix_([0, 4, 2], [0, 4, 2])]
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

    def test_column_support_and_sharp_total_entry_mass(self):
        for n in (2, 3, 4, 5):
            size = 1 << n
            bound = (1 + np.sqrt(2)) * size - np.sqrt(2 * size)
            balanced = np.full(size - 1, np.pi / 4)
            unequal = .19 + .31 * np.sin(np.arange(size - 1) + .4)
            singular = np.resize([0, np.pi / 2, -np.pi, .23], size - 1)
            for angles in (balanced, unequal, singular):
                frame = direct_real_frame(n, angles)
                self.assertLessEqual(np.abs(frame).sum(), bound + ATOL)
                for node in range(1, size):
                    start = anchor_label(node, n)
                    marker = marker_label(node, n)
                    end = start + 2 * (marker - start)
                    outside = np.r_[0:start, end:size]
                    np.testing.assert_allclose(frame[outside, marker], 0,
                                               atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.abs(direct_real_frame(n, balanced)).sum(),
                                   bound, delta=ATOL)
            # The global mass bound does not imply bounded row normalization.
            for depth in range(n):
                balanced[(1 << depth) - 1] = np.arccos(
                    np.sqrt((depth + 1) / (depth + 2)))
            row = direct_real_frame(n, balanced)[0]
            nonzero = np.abs(row[np.abs(row) > ATOL])
            self.assertEqual(len(nonzero), n + 1)
            np.testing.assert_allclose(nonzero ** 2, 1 / (n + 1),
                                       atol=ATOL, rtol=0)
            self.assertAlmostEqual(np.abs(row).sum(), np.sqrt(n + 1), delta=ATOL)

    def test_spin_covariance_and_overlap_on_every_dirty_column(self):
        for n in (2, 3):
            size = 1 << n
            gammas = _majoranas(size // 2)
            identity = np.eye(len(gammas[0]))
            embedding = np.vstack(gammas) / np.sqrt(size)
            angles = .37 * np.sin(np.arange(size - 1) + .2)
            lift = identity.copy()
            for node, theta in enumerate(angles, start=1):
                a, b = anchor_label(node, n), marker_label(node, n)
                plane = (np.cos(theta / 2) * identity
                         - np.sin(theta / 2) * gammas[a] @ gammas[b])
                lift = plane @ lift
            frame = direct_real_frame(n, angles)
            np.testing.assert_allclose(embedding.conj().T @ embedding, identity,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(np.kron(frame, lift) @ embedding,
                                       embedding @ lift, atol=ATOL, rtol=0)
            overlap = np.trace(frame) * identity
            for j in range(size):
                for k in range(j + 1, size):
                    overlap = overlap + (frame[j, k] - frame[k, j]) * gammas[j] @ gammas[k]
            np.testing.assert_allclose(
                embedding.conj().T @ np.kron(frame, identity) @ embedding,
                overlap / size, atol=ATOL, rtol=0)

    def test_single_plane_reflection_distance_and_overlap_singular_values(self):
        for n in (2, 3):
            size = 1 << n
            gammas = _majoranas(size // 2)
            identity = np.eye(len(gammas[0]))
            embedding = np.vstack(gammas) / np.sqrt(size)
            reflection = 2 * embedding @ embedding.conj().T - np.eye(size * len(identity))
            for theta in (0, .37, np.pi):
                angles = np.zeros(size - 1)
                angles[0] = theta
                frame = direct_real_frame(n, angles)
                a, b = anchor_label(1, n), marker_label(1, n)
                lift = (np.cos(theta / 2) * identity
                        - np.sin(theta / 2) * gammas[a] @ gammas[b])
                spin = np.kron(np.eye(size), lift)
                vector = np.kron(frame, identity)
                rotated = spin @ reflection @ spin.conj().T
                np.testing.assert_allclose(rotated,
                                           vector.conj().T @ reflection @ vector,
                                           atol=ATOL, rtol=0)
                overlap = embedding.conj().T @ vector @ embedding
                square = 1 - 4 * (size - 2) * (1 - np.cos(theta)) / size ** 2
                np.testing.assert_allclose(overlap.conj().T @ overlap,
                                           square * identity, atol=ATOL, rtol=0)
                distance = max(abs(np.linalg.eigvalsh(rotated - reflection)))
                expected = 4 * np.sqrt((size - 2) * (1 - np.cos(theta))) / size
                self.assertAlmostEqual(distance, expected, delta=ATOL)

    def test_tree_bivector_labels_are_independent_over_gf2(self):
        # Binary Pauli exponents only: no exponential spinor matrices needed.
        for n in (2, 3, 4, 5):
            size, q = 1 << n, 1 << (n - 1)
            labels = []
            for wire in range(q):
                prefix = (1 << wire) - 1
                labels += [(1 << wire) | (prefix << q),
                           (1 << wire) | ((prefix | (1 << wire)) << q)]
            self.assertEqual(_binary_rank(labels), size)
            edges = [labels[anchor_label(node, n)] ^ labels[marker_label(node, n)]
                     for node in range(1, size)]
            self.assertEqual(_binary_rank(edges), size - 1)


if __name__ == '__main__':
    unittest.main()
