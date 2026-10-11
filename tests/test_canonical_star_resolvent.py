"""Complete-operator checks for canonical stars and their sparse resolvent.

These tests compare independently composed local stars, the Cayley and
triangular formulas, and the full midpoint-reflection oracle.  They do
not implement or certify native synthesis of the feedback inverse.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.conventions import anchor_label, marker_label
from compiler_robust_hopf.frames import direct_real_frame, hopf_ry


ATOL = 8e-12


def _children(height):
    children = [[] for _ in range(1 << height)]
    for label in range(1, 1 << height):
        children[label - (label & -label)].append(label)
    return children


def _canonical_star(root_column):
    c = float(root_column[0].real)
    u = np.asarray(root_column[1:], dtype=complex)
    star = np.eye(len(root_column), dtype=complex)
    star[0, 0] = c
    star[1:, 0] = u
    star[0, 1:] = -u.conj()
    star[1:, 1:] -= np.outer(u, u.conj()) / (1 + c)
    return star


def _compose_stars(children, stars):
    size = len(children)
    complete = np.eye(size, dtype=complex)
    # Numerical label order is topological; incomparable stars commute.
    for vertex, kids in enumerate(children):
        labels = [vertex, *kids]
        embedded = np.eye(size, dtype=complex)
        embedded[np.ix_(labels, labels)] = stars[vertex]
        complete = embedded @ complete
    return complete


def _ordered_hopf_stars(height, angles):
    children = _children(height)
    stars = [np.eye(len(kids) + 1, dtype=complex) for kids in children]
    for node, theta in enumerate(angles, start=1):
        vertex = anchor_label(node, height)
        child = marker_label(node, height)
        slot = children[vertex].index(child) + 1
        rotation = np.eye(len(children[vertex]) + 1, dtype=complex)
        rotation[np.ix_([0, slot], [0, slot])] = hopf_ry(theta)
        stars[vertex] = rotation @ stars[vertex]
    return children, stars


def _shift(children, stars):
    size = len(children)
    shift = np.zeros((size, size), dtype=complex)
    carrier = np.zeros_like(shift)
    for vertex, kids in enumerate(children):
        column = stars[vertex][:, 0]
        shift[kids, vertex] = column[1:] / (1 + column[0].real)
        carrier[kids, vertex] = column[1:]
    return shift, carrier


def _complex_fixture(height):
    rng = np.random.default_rng(817 + height)
    children = _children(height)
    stars = []
    for vertex, kids in enumerate(children):
        if not kids:
            stars.append(np.eye(1, dtype=complex))
            continue
        c = (0.0, 0.02, 0.6, 1.0)[(vertex // 2) % 4]
        u = rng.normal(size=len(kids)) + 1j * rng.normal(size=len(kids))
        u *= np.sqrt(1 - c * c) / np.linalg.norm(u)
        stars.append(_canonical_star(np.r_[c, u]))
    return children, stars


class CanonicalStarResolventTests(unittest.TestCase):
    def test_same_root_transfer_preserves_the_complete_hopf_frame(self):
        for height in (1, 2, 3, 4):
            size = 1 << height
            for angles in (np.linspace(-.71, .67, size - 1),
                           np.zeros(size - 1),
                           np.where(np.arange(size - 1) % 3, .31, 0)):
                with self.subTest(height=height, angles=angles):
                    children, original = _ordered_hopf_stars(height, angles)
                    canonical = [_canonical_star(star[:, 0]) for star in original]
                    correction = np.eye(size, dtype=complex)
                    for vertex, kids in enumerate(children):
                        detail = canonical[vertex].conj().T @ original[vertex]
                        np.testing.assert_allclose(detail[:, 0],
                                                   np.eye(len(detail))[:, 0],
                                                   atol=ATOL, rtol=0)
                        correction[np.ix_(kids, kids)] = detail[1:, 1:]
                    reference = direct_real_frame(height, angles)
                    np.testing.assert_allclose(_compose_stars(children, original),
                                               reference, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(
                        _compose_stars(children, canonical) @ correction,
                        reference, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(correction.conj().T @ correction,
                                               np.eye(size), atol=ATOL, rtol=0)

    def test_cayley_and_sparse_forms_against_independent_complex_star_word(self):
        for height in (1, 2, 3, 4, 5):
            with self.subTest(height=height):
                children, stars = _complex_fixture(height)
                size = len(children)
                shift, carrier = _shift(children, stars)
                identity = np.eye(size)
                f = identity - shift
                r = np.linalg.solve(f, identity)
                x = r.conj().T - r
                d = identity + shift.conj().T @ shift
                m = d - 2 * shift
                reference = _compose_stars(children, stars)
                cayley = 2 * np.linalg.solve(identity + x, identity) - identity
                sparse = 2 * f @ np.linalg.solve(m, f.conj().T) - identity
                np.testing.assert_allclose(cayley, reference, atol=ATOL, rtol=0)
                np.testing.assert_allclose(sparse, reference, atol=ATOL, rtol=0)
                np.testing.assert_allclose(reference.conj().T @ reference,
                                           identity, atol=ATOL, rtol=0)
                for congruence in (f.conj().T @ x @ f, f @ x @ f.conj().T):
                    np.testing.assert_allclose(congruence, shift.conj().T - shift,
                                               atol=ATOL, rtol=0)
                np.testing.assert_allclose(m, (identity - carrier) @ d,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(d, np.diag(np.diag(d)), atol=ATOL, rtol=0)
                self.assertLessEqual(np.linalg.norm(x, 2), 2 * height + ATOL)
                self.assertLessEqual(np.linalg.norm(np.linalg.solve(m, identity), 2),
                                     height + 1 + ATOL)

    def test_root_included_shift_has_the_correct_nilpotency_index(self):
        height = 4
        children = _children(height)
        stars = []
        for kids in children:
            if kids:
                u = np.full(len(kids), np.sqrt(.75 / len(kids)))
                stars.append(_canonical_star(np.r_[.5, u]))
            else:
                stars.append(np.eye(1))
        shift, _ = _shift(children, stars)
        self.assertGreater(np.linalg.norm(np.linalg.matrix_power(shift, height)), 0)
        np.testing.assert_array_equal(np.linalg.matrix_power(shift, height + 1), 0)

    def test_midpoint_source_gives_the_full_local_oracle_and_literal_dummy_sign(self):
        for height in (1, 2, 3, 4):
            with self.subTest(height=height):
                children, stars = _complex_fixture(height)
                size = len(children)
                shift, _ = _shift(children, stars)
                diagonal = 1 + np.sum(np.abs(shift) ** 2, axis=0)
                v = np.vstack((np.eye(size), shift)) / np.sqrt(diagonal)
                source = np.eye(2 * size, dtype=complex)
                expected = np.eye(2 * size, dtype=complex)
                z = np.diag(np.r_[np.ones(size), -np.ones(size)])
                for vertex, kids in enumerate(children):
                    labels = [vertex, *(size + child for child in kids)]
                    expected[np.ix_(labels, labels)] = stars[vertex]
                    anchor = np.eye(2 * size)[:, vertex]
                    midpoint = anchor + v[:, vertex]
                    midpoint /= np.linalg.norm(midpoint)
                    hw = np.eye(2 * size) - 2 * np.outer(midpoint, midpoint.conj())
                    he = np.eye(2 * size) - 2 * np.outer(anchor, anchor)
                    source = hw @ he @ source
                np.testing.assert_allclose(source[:, :size], v, atol=ATOL, rtol=0)
                np.testing.assert_allclose(source.conj().T @ source,
                                           np.eye(2 * size), atol=ATOL, rtol=0)
                oracle = source @ z @ source.conj().T @ z
                np.testing.assert_allclose(oracle, expected, atol=ATOL, rtol=0)
                # The unused lower root must have +1, not a discarded -1.
                np.testing.assert_array_equal(oracle[:, size],
                                              np.eye(2 * size)[:, size])
                for vertex, kids in enumerate(children):
                    if not kids:
                        np.testing.assert_array_equal(source[:, vertex],
                                                      np.eye(2 * size)[:, vertex])

    def test_zero_angles_give_identity_in_both_the_frame_and_full_oracle(self):
        height = 3
        children, stars = _ordered_hopf_stars(height, np.zeros((1 << height) - 1))
        shift, carrier = _shift(children, stars)
        np.testing.assert_array_equal(shift, 0)
        np.testing.assert_array_equal(carrier, 0)
        np.testing.assert_array_equal(_compose_stars(children, stars), np.eye(1 << height))
        v = np.vstack((np.eye(1 << height), np.zeros((1 << height, 1 << height))))
        z = np.diag(np.r_[np.ones(1 << height), -np.ones(1 << height)])
        np.testing.assert_array_equal((2 * v @ v.T - np.eye(2 << height)) @ z,
                                      np.eye(2 << height))

    def test_leaf_omission_packs_one_bin_and_precision_core_fits(self):
        for height in range(1, 12):
            size = 1 << height
            children = _children(height)
            degrees = [len(kids) for kids in children]
            self.assertEqual(sum(degree ** 2 for degree in degrees),
                             3 * size - 2 * height - 3)
            padded = [1 << degree.bit_length() for degree in degrees if degree]
            self.assertLessEqual(sum(padded), 2 * size - 2)
            cursor = 0
            for block in sorted(padded, reverse=True):
                self.assertEqual(cursor % block, 0)
                cursor += block
            self.assertLess(cursor, 2 * size)
            if height >= 3:
                precision = size + 9
                self.assertLessEqual(precision + 1, size + height + 7)
                self.assertGreaterEqual(precision + 1, 2 * (height + 1))
                self.assertLess(30 * 8 * 2.0 ** (size - precision), 1)


if __name__ == "__main__":
    unittest.main()
