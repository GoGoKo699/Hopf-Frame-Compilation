"""Bounded algebra audits for a signed coarse-frame sparse correction.

The fixtures include the state column, singular branches, negative real
coefficients, and every marker column.  Source completions below are dense
verification matrices, not native source oracles.  The small LCU word keeps
all source and branch work through its actual inverse and literal recovery
phase.  These checks establish no resource bound or scalable emitter.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix

try:
    from .test_weighted_haar_depth import _completion, _nodes, _source_index
except ImportError:
    from test_weighted_haar_depth import _completion, _nodes, _source_index


ATOL = 5e-12


def _frames(height, *, perturbation=None):
    """Literal pair layers, with exact zero and negatively oriented branches."""
    size = 1 << height
    pairs = [(0., -1.), (-.6, .8), (-1., 0.), (.8, -.6), (1., 0.)]
    target, coarse = np.eye(size), np.eye(size)
    angles, errors = [], []
    step = perturbation if perturbation is not None else 1 / (128 * (height + 1)**2)
    for node in range(1, size):
        depth = node.bit_length() - 1
        prefix = node - (1 << depth)
        left = prefix << (height - depth)
        right = left + (1 << (height - depth - 1))
        cosine, sine = pairs[(node - 1) % len(pairs)]
        delta = step * (1 if node % 2 else -.75)
        c0 = cosine * np.cos(delta) - sine * np.sin(delta)
        s0 = sine * np.cos(delta) + cosine * np.sin(delta)
        local = np.array([[cosine, -sine], [sine, cosine]])
        rounded = np.array([[c0, -s0], [s0, c0]])
        target[[left, right]] = local @ target[[left, right]]
        coarse[[left, right]] = rounded @ coarse[[left, right]]
        angles.append(np.arctan2(sine, cosine))
        errors.append(np.linalg.norm(local - rounded, 2))
    # The independent addressed-layer fixture uses the repository's
    # prescribed marker convention, including all singular continuations.
    np.testing.assert_allclose(target, real_frame_matrix(angles), atol=ATOL, rtol=0)
    return target, coarse, max(errors)


def _signed_sources(residual):
    """One leaf, one depth tag, and two disjoint failure coordinates."""
    size = len(residual)
    height = (size - 1).bit_length()
    tags = 1 << height.bit_length()
    positive = np.abs(residual)
    doubled = np.block([[np.zeros_like(positive), positive],
                        [positive.T, np.zeros_like(positive)]])
    weights = np.linalg.solve(np.eye(2 * size) - doubled, np.ones(2 * size))
    u, v = weights[:size], weights[size:]
    row = np.zeros((size * tags * 4, size), dtype=complex)
    column = np.zeros_like(row)
    # Unlike the bounded-density correction, the root column has a live
    # overlap component.  Tag zero is reserved for this root on every leaf.
    sectors = [(0, 0, 0, size), *_nodes(size)]
    for marker, tag, start, length in sectors:
        for leaf in range(start, start + length):
            location = _source_index(size, leaf, tag, 0)
            magnitude = positive[leaf, marker]
            column[location, marker] = np.sign(residual[leaf, marker]) * np.sqrt(
                magnitude * u[leaf] / v[marker])
        column[_source_index(size, start, tag, 1), marker] = 1 / np.sqrt(v[marker])
    # Recover each marker independently from leaf and tag.  In particular,
    # this is not a second retained copy of the input marker register.
    for leaf in range(size):
        for tag in range(height + 1):
            marker = 0 if tag == 0 else (
                (2 * (leaf >> (height - tag + 1)) + 1) << (height - tag))
            row[_source_index(size, leaf, tag, 0), leaf] = np.sqrt(
                positive[leaf, marker] * v[marker] / u[leaf])
        row[_source_index(size, leaf, 0, 2), leaf] = 1 / np.sqrt(u[leaf])
    return row, column, u, v


def _lcu_word(target, coarse, *, perturbation=0.):
    """Actual unitary completions and a normalization-two Hadamard LCU."""
    row, column, _, _ = _signed_sources(target - coarse)
    left, right = _completion(row), _completion(column)
    size, dimension = len(target), len(left)
    # The source promise fixes only its initialized columns.  A generic
    # permitted action on other columns makes missing work coordinates in
    # a recovery reflection observable instead of relying on a specially
    # aligned QR completion.
    complement, _ = np.linalg.qr(np.random.default_rng(491).normal(
        size=(dimension - size, dimension - size)))
    left[:, size:] = left[:, size:] @ complement
    if perturbation:
        # A genuine source error couples main overlap to column failure;
        # retain this changed word's own adjoint in the full recovery.
        i = _source_index(size, 0, 0, 0)
        j = _source_index(size, 0, 0, 1)
        c, s = np.cos(perturbation), np.sin(perturbation)
        right[[i, j]] = np.array([[c, -s], [s, c]]) @ right[[i, j]]
    residual_word = left.conj().T @ right
    coarse_word = np.kron(np.eye(dimension // size), coarse)
    zero = np.zeros_like(coarse_word)
    selected = np.block([[coarse_word, zero], [zero, residual_word]])
    prep = np.kron(np.array([[1, 1], [1, -1]]) / np.sqrt(2), np.eye(dimension))
    return prep.conj().T @ selected @ prep, prep


def _recover(word, size, *, omit_tag_from_anchor=False, inverse=None):
    reflection = np.eye(len(word), dtype=complex)
    locations = list(range(size))
    if omit_tag_from_anchor:
        source_dimension = len(word) // 2
        locations = [index for index in range(source_dimension)
                     if (index // size) % 4 == 0]
    reflection[locations, locations] = -1
    actual_inverse = word.conj().T if inverse is None else inverse
    return -word @ reflection @ actual_inverse @ reflection @ word


class CoarseFlagSparseCorrectionTests(unittest.TestCase):
    def test_absolute_residual_bound_at_singular_negative_branches(self):
        for height in range(1, 5):
            with self.subTest(height=height):
                target, coarse, node_error = _frames(height)
                residual = target - coarse
                size = len(target)
                allowed = np.zeros((size, size), dtype=bool)
                allowed[:, 0] = True
                for marker, _, start, length in _nodes(size):
                    allowed[start:start + length, marker] = True
                self.assertTrue(np.all(allowed.sum(axis=1) == height + 1))
                np.testing.assert_array_equal(residual[~allowed], 0)
                self.assertTrue(np.any(target == 0) and np.any(target < 0))
                columns = np.linalg.norm(residual, axis=0)
                self.assertGreater(columns[0], 1e-6)
                self.assertLessEqual(columns.max(), height * node_error + ATOL)
                self.assertLessEqual(np.linalg.norm(np.abs(residual), 2),
                                     np.sqrt(height + 1) * columns.max() + ATOL)
                self.assertLess(np.linalg.norm(np.abs(residual), 2), 1 / 16)
                np.testing.assert_allclose(coarse.T @ coarse, np.eye(size),
                                           atol=ATOL, rtol=0)

    def test_signed_single_leaf_sources_include_root_tag_and_zero_columns(self):
        for height in range(1, 4):
            target, coarse, _ = _frames(height)
            for residual in (target - coarse, np.zeros_like(target)):
                with self.subTest(height=height, zero=not np.any(residual)):
                    row, column, u, v = _signed_sources(residual)
                    size = len(residual)
                    np.testing.assert_allclose(row.conj().T @ row, np.eye(size),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(column.conj().T @ column, np.eye(size),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(row.conj().T @ column, residual,
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(abs(residual) @ v, u - 1, atol=ATOL)
                    np.testing.assert_allclose(abs(residual).T @ u, v - 1, atol=ATOL)
                    self.assertTrue(np.all(u >= 1 - ATOL) and np.all(v >= 1 - ATOL))
                    if np.any(residual):
                        self.assertGreater(np.linalg.norm(column[:size, 0]), .001)
                        # Keeping the old pure-failure root silently deletes
                        # an entire nonzero state-column correction.
                        wrong = column.copy()
                        wrong[:, 0] = 0
                        wrong[_source_index(size, 0, 0, 1), 0] = 1
                        self.assertGreater(np.linalg.norm(row.conj().T @ wrong
                                                          - residual), 1e-6)
                        unsigned = np.abs(column)
                        np.testing.assert_allclose(row.conj().T @ unsigned,
                                                   abs(residual), atol=ATOL, rtol=0)
                        self.assertGreater(np.linalg.norm(abs(residual) - residual), 1e-6)

    def test_normalization_two_recovery_literal_phase_and_all_work_anchor(self):
        target, coarse, _ = _frames(2, perturbation=.012)
        size = len(target)
        word, prep = _lcu_word(target, coarse)
        embedding = np.eye(len(word), dtype=complex)[:, :size]
        np.testing.assert_allclose(embedding.conj().T @ word @ embedding,
                                   target / 2, atol=ATOL, rtol=0)
        recovered = _recover(word, size)
        wanted = embedding @ target
        np.testing.assert_allclose(recovered @ embedding, wanted, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(-recovered @ embedding - wanted, 2), 1.99)
        wrong_anchor = _recover(word, size, omit_tag_from_anchor=True)
        self.assertGreater(np.linalg.norm(wrong_anchor @ embedding - wanted, 2), .01)

        # H=0 uses actual PREP/inverse cancellation, including occupied
        # branch and source banks.  Coherent H and complex dirty/reference
        # entanglement also detect a relative minus sign on the active side.
        inactive = prep.conj().T @ np.eye(len(word)) @ prep
        np.testing.assert_allclose(inactive, np.eye(len(word)), atol=ATOL, rtol=0)
        rng = np.random.default_rng(20261016)
        state = rng.normal(size=(2, len(word), 2, 2)) + 1j * rng.normal(
            size=(2, len(word), 2, 2))
        state[1, size:] = 0
        state /= np.linalg.norm(state)
        actual = np.stack((np.einsum('ab,bdr->adr', inactive, state[0]),
                           np.einsum('ab,bdr->adr', recovered, state[1])))
        expected = state.copy()
        expected[1, :size] = np.einsum('ab,bdr->adr', target, state[1, :size])
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

    def test_perturbed_sources_use_actual_inverse_and_retain_final_leakage(self):
        target, coarse, _ = _frames(2, perturbation=.012)
        exact, _ = _lcu_word(target, coarse)
        angle = .0017
        actual, _ = _lcu_word(target, coarse, perturbation=angle)
        size = len(target)
        embedding = np.eye(len(actual), dtype=complex)[:, :size]
        result = _recover(actual, size) @ embedding
        desired = embedding @ target
        self.assertLessEqual(np.linalg.norm(result - desired, 2),
                             3 * 2 * np.sin(angle / 2) + ATOL)
        self.assertGreater(np.linalg.norm(result[size:]), 1e-6)
        np.testing.assert_allclose(result.conj().T @ result, np.eye(size),
                                   atol=ATOL, rtol=0)
        stale_inverse = _recover(actual, size, inverse=exact.conj().T) @ embedding
        self.assertGreater(np.linalg.norm(stale_inverse - result), 1e-6)


if __name__ == '__main__':
    unittest.main()
