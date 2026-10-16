"""Bounded algebraic checks for interval dilation and anchored recovery.

These fixtures retain all rejection registers and literal phases. Small
Householder completions verify interfaces; they are not native emitters.
The finite overlap-norm check does not establish the asymptotic norm bound,
and no resource assertion is inferred from these matrix calculations.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import real_frame_matrix, real_tree_data


ATOL = 4e-11


def _interval_fixture(q=2, tau=2, numerators=(1, 2, 3)):
    denominator = 1 << tau
    leaf_count, atom_count = 1 << q, 1 << (q * tau)
    intervals = {1: (0, atom_count)}
    angles = []
    for node, numerator in enumerate(numerators, 1):
        start, stop = intervals[node]
        split = start + (stop - start) * numerator // denominator
        intervals[2 * node] = (start, split)
        intervals[2 * node + 1] = (split, stop)
        angles.append(np.arccos(np.sqrt(numerator / denominator)))
    embedding = np.zeros((atom_count, leaf_count))
    membership = np.empty(atom_count, dtype=int)
    for x in range(leaf_count):
        start, stop = intervals[leaf_count + x]
        embedding[start:stop, x] = 1 / np.sqrt(stop - start)
        membership[start:stop] = x
    frame = real_frame_matrix(angles)
    balanced = real_frame_matrix(np.full(atom_count - 1, np.pi / 4))
    overlap = balanced.T @ embedding @ frame
    return intervals, embedding, membership, angles, frame, balanced, overlap


def _completion(target, phase=0.0):
    """Unitary with a real literal first column and free complex completion."""
    target = np.asarray(target, dtype=complex)
    anchor = np.zeros(len(target), dtype=complex)
    anchor[0] = 1
    difference = anchor - target
    if np.linalg.norm(difference) < 1e-13:
        unitary = np.eye(len(target), dtype=complex)
    else:
        difference /= np.linalg.norm(difference)
        unitary = np.eye(len(target)) - 2 * np.outer(difference, difference.conj())
    phases = np.exp(1j * phase * np.arange(len(target)))
    return unitary * phases[None, :]


def _apply_axes(state, matrix, axes):
    axes = tuple(axes)
    moved = np.moveaxis(state, axes, range(len(axes)))
    shape = moved.shape
    acted = matrix @ moved.reshape(matrix.shape[1], -1)
    return np.moveaxis(acted.reshape(shape), range(len(axes)), axes)


def _controlled(state, matrices, control, targets):
    result = np.empty_like(state)
    shifted = tuple(axis - (axis > control) for axis in targets)
    for value, matrix in enumerate(matrices):
        selection = [slice(None)] * state.ndim
        selection[control] = value
        result[tuple(selection)] = _apply_axes(state[tuple(selection)], matrix, shifted)
    return result


def _membership_xor(state, membership, target_axis, atom_axis):
    result = np.empty_like(state)
    for z, leaf in enumerate(membership):
        for x in range(state.shape[target_axis]):
            source = [slice(None)] * state.ndim
            destination = list(source)
            source[target_axis], source[atom_axis] = x, z
            destination[target_axis], destination[atom_axis] = x ^ int(leaf), z
            result[tuple(destination)] = state[tuple(source)]
    return result


def _rotation(angle):
    return np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle), np.cos(angle)]])


class _SquareBlock:
    """Matrix-free q=1 physical fixture with separate rejection registers.

    Axes are X,Y,Z,F,A,A_E,damping, followed by untouched optional axes.
    X and Y each have two values, Z four, F four, and the source banks two.
    The only dense primitives below act on at most 32 basis states.
    """

    shape = (2, 2, 4, 4, 2, 2, 2)

    def __init__(self, *, phase=0.0, leakage=0.0):
        _, self.embedding, self.membership, _, self.frame, self.balanced, overlap = (
            _interval_fixture(q=1, tau=2, numerators=(1,))
        )
        self.alpha = 2.0
        magnitude = np.abs(overlap)
        symmetric = np.block([[np.zeros((4, 4)), magnitude],
                              [magnitude.T, np.zeros((2, 2))]])
        weights = np.linalg.solve(np.eye(6) - symmetric / self.alpha, np.ones(6))
        f, g = weights[:4], weights[4:]
        c = np.sign(overlap) * np.sqrt(magnitude * f[:, None]
                                       / (self.alpha * g[None, :]))
        ell = np.sqrt(magnitude * g[None, :]
                      / (self.alpha * f[:, None]))
        self.columns, self.rows, self.intervals = [], [], []
        for x in range(2):
            target = np.zeros((4, 4, 2))
            target[:, 0, 0] = c[:, x]
            target[0, 1, 0] = 1 / np.sqrt(g[x])
            word = _completion(target.ravel(), phase + 0.17 * x)
            self.columns.append(np.kron(np.eye(16), _rotation(leakage)) @ word)
        for y in range(4):
            target = np.zeros((2, 4, 2))
            target[:, 0, 0] = ell[y]
            target[0, 2, 0] = 1 / np.sqrt(f[y])
            word = _completion(target.ravel(), phase - 0.11 * y)
            self.rows.append(np.kron(np.eye(8), _rotation(-leakage)) @ word)
        for x in range(2):
            target = np.zeros((4, 2))
            target[:, 0] = self.embedding[:, x]
            word = _completion(target.ravel(), phase + 0.23 * x)
            self.intervals.append(np.kron(np.eye(4), _rotation(leakage)) @ word)
        self.k = 3
        self.angle = np.pi / (4 * self.k + 2)
        self.damping = _rotation(np.arccos(self.alpha * np.sin(self.angle)))

    def anchor(self, columns=None):
        if columns is None:
            columns = np.eye(2)
        state = np.zeros(self.shape + columns.shape[1:], dtype=complex)
        state[:, 0, 0, 0, 0, 0, 0] = columns
        return state

    def interval_word(self, state, *, inverse=False):
        if inverse:
            state = _membership_xor(state, self.membership, 1, 2)
            return _controlled(state, [u.conj().T for u in self.intervals], 1, (2, 5))
        state = _controlled(state, self.intervals, 1, (2, 5))
        return _membership_xor(state, self.membership, 1, 2)

    def square(self, state, *, inverse=False):
        if inverse:
            state = np.swapaxes(state, 0, 1).copy()
            state = self.interval_word(state)
            state = _apply_axes(state, self.balanced.T, (2,))
            state = _controlled(state, self.rows, 2, (0, 3, 4))
            return _controlled(state, [u.conj().T for u in self.columns], 0, (2, 3, 4))
        state = _controlled(state, self.columns, 0, (2, 3, 4))
        state = _controlled(state, [u.conj().T for u in self.rows], 2, (0, 3, 4))
        state = _apply_axes(state, self.balanced, (2,))
        state = self.interval_word(state, inverse=True)
        return np.swapaxes(state, 0, 1).copy()

    def unitary(self, state, *, inverse=False):
        if inverse:
            state = _apply_axes(state, self.damping.T, (6,))
            return self.square(state, inverse=True)
        return _apply_axes(self.square(state), self.damping, (6,))

    @staticmethod
    def reflect(state):
        result = state.copy()
        result[:, 0, 0, 0, 0, 0, 0] *= -1
        return result

    def amplified(self, state, *, inverse=False):
        if inverse:
            for _ in range(self.k):
                state = self.unitary(state, inverse=True)
                state = self.reflect(state)
                state = self.unitary(state)
                state = -self.reflect(state)
            return self.unitary(state, inverse=True)
        state = self.unitary(state)
        for _ in range(self.k):
            state = self.reflect(state)
            state = self.unitary(state, inverse=True)
            state = self.reflect(state)
            state = -self.unitary(state)
        return state


class IntervalDilationDepthTests(unittest.TestCase):
    def test_integer_embedding_and_every_prescribed_prefix_flag(self):
        intervals, embedding, membership, angles, frame, _, _ = _interval_fixture()
        self.assertEqual([intervals[4 + x][1] - intervals[4 + x][0]
                          for x in range(4)], [2, 2, 9, 3])
        np.testing.assert_allclose(embedding.T @ embedding, np.eye(4), atol=ATOL)
        tree = real_tree_data(angles)
        for node in range(1, 8):
            start, stop = intervals[node]
            uniform = np.zeros(16)
            uniform[start:stop] = 1 / np.sqrt(stop - start)
            np.testing.assert_allclose(embedding @ tree.subtree[node], uniform, atol=ATOL)
        for node in range(1, 4):
            start, split = intervals[2 * node]
            _, stop = intervals[2 * node + 1]
            left, right = split - start, stop - split
            expected = np.zeros(16)
            expected[start:split] = -np.sqrt(right / (left * (left + right)))
            expected[split:stop] = np.sqrt(left / (right * (left + right)))
            np.testing.assert_allclose(embedding @ frame[:, marker_label(node, 2)],
                                       expected, atol=ATOL)
        self.assertTrue(np.array_equal(membership[:4], [0, 0, 1, 1]))

    def test_rectangular_overlap_root_sign_sparsity_and_norm_diagnostic(self):
        _, embedding, _, _, frame, balanced, overlap = _interval_fixture()
        np.testing.assert_allclose(overlap.T @ overlap, np.eye(4), atol=ATOL)
        np.testing.assert_allclose(balanced @ overlap, embedding @ frame, atol=ATOL)
        np.testing.assert_allclose(overlap[:, 0], np.eye(16)[:, 0], atol=ATOL)
        np.testing.assert_allclose(overlap[0, 1:], 0, atol=ATOL)
        for depth in range(4):
            rows = [marker_label((1 << depth) + p, 4) for p in range(1 << depth)]
            self.assertTrue(np.all(np.count_nonzero(np.abs(overlap[rows, 1:]) > ATOL,
                                                   axis=0) <= 3))
        alpha_bound = 1 + np.sqrt(3 * 4) + (np.sqrt(6) + 6) / (1 - 2 ** -0.5)
        self.assertLess(np.linalg.norm(np.abs(overlap), 2), alpha_bound)

    def test_membership_erasure_and_actual_inverse_ignore_source_completion(self):
        _, embedding, membership, _, _, _, _ = _interval_fixture()
        for phase in (0.0, 0.61):
            sources = [_completion(embedding[:, x], phase * (x + 1)) for x in range(4)]
            state = np.zeros((4, 16, 4), dtype=complex)
            state[:, 0] = np.eye(4)
            prepared = _controlled(state, sources, 0, (1,))
            erased = _membership_xor(prepared, membership, 0, 1)
            expected = np.zeros_like(erased)
            expected[0] = embedding
            np.testing.assert_allclose(erased, expected, atol=ATOL)
            # Test E^dagger on the entire atom space, including its orthogonal
            # complement, rather than only inverse action on prepared columns.
            atom_inputs = np.zeros((4, 16, 16), dtype=complex)
            atom_inputs[0] = np.eye(16)
            restored = _membership_xor(atom_inputs, membership, 0, 1)
            restored = _controlled(restored, [u.conj().T for u in sources], 0, (1,))
            np.testing.assert_allclose(restored[:, 0], embedding.T, atol=ATOL)

    def test_square_anchor_preserves_all_rejection_coordinates(self):
        for phase in (0.0, 0.79):
            block = _SquareBlock(phase=phase)
            output = block.square(block.anchor())
            accepted = output[:, 0, 0, 0, 0, 0, 0]
            np.testing.assert_allclose(accepted, block.frame / block.alpha, atol=ATOL)
            flat = output.reshape(-1, 2)
            np.testing.assert_allclose(flat.conj().T @ flat, np.eye(2), atol=ATOL)
            # Old-X rejection is retained in final Y; failure and source-work
            # coordinates are never silently projected out by the fixture.
            self.assertGreater(np.linalg.norm(output[:, 1:]), 0.01)
            np.testing.assert_allclose(block.square(output, inverse=True), block.anchor(),
                                       atol=ATOL)

    def test_calibrated_amplification_literal_phase_and_full_leakage_bound(self):
        exact = _SquareBlock(phase=0.47)
        output = exact.amplified(exact.anchor())
        np.testing.assert_allclose(output, exact.anchor(exact.frame), atol=ATOL)
        leakage = 2e-3
        actual = _SquareBlock(phase=0.47, leakage=leakage)
        output = actual.amplified(actual.anchor())
        error = np.linalg.norm((output - actual.anchor(actual.frame)).reshape(-1, 2), 2)
        self.assertLess(error, 8 * (2 * actual.k + 1) * leakage)
        occupied = output.copy()
        occupied[:, 0, 0, 0, 0, 0, 0] = 0
        # This symmetric source-work perturbation cancels to first order
        # in the accepted block, but its complete output still leaks.
        self.assertGreater(np.linalg.norm(occupied), leakage**2 / 10)
        np.testing.assert_allclose(actual.amplified(output, inverse=True), actual.anchor(),
                                   atol=ATOL)

    def test_fine_tag_erasure_all_bands_with_actual_inverse_and_reference(self):
        _, _, _, _, frame, _, _ = _interval_fixture()
        # A small occupied-work perturbation verifies the complete hybrid,
        # including arbitrary reference correlations and actual reversed words.
        epsilon = 3e-5
        compiler = np.kron(frame, np.eye(2))
        perturbation = np.eye(8)
        perturbation[np.ix_([0, 3], [0, 3])] = _rotation(epsilon)
        actual = perturbation @ compiler
        tags = [0, 2, 1, 2]
        state = np.zeros((8, 4, 4), dtype=complex)
        for j, tag in enumerate(tags):
            state[::2, tag, j] = frame[:, j]
        decoded = np.einsum("ab,btj->atj", actual.conj().T, state)
        erased = np.empty_like(decoded)
        for x in range(4):
            for work in range(2):
                for tag in range(4):
                    erased[2 * x + work, tag ^ tags[x]] = decoded[2 * x + work, tag]
        restored = np.einsum("ab,btj->atj", actual, erased)
        expected = np.zeros_like(restored)
        expected[::2, 0] = frame
        full_error = np.linalg.norm((restored - expected).reshape(32, 4), 2)
        self.assertGreater(full_error, epsilon / 10)
        self.assertLess(full_error, 2 * np.linalg.norm(perturbation - np.eye(8), 2) + ATOL)
        reference = np.array([[1, 2j], [-3j, 1], [2, -1j], [1j, 4]], dtype=complex)
        reference /= np.linalg.norm(reference)
        np.testing.assert_allclose(np.einsum("atj,jr->atr", restored, reference),
                                   np.einsum("atj,jr->atr", expected, reference),
                                   atol=2 * epsilon)

    def test_exact_guarded_balanced_word_variable_width_and_occupied_padding(self):
        lengths, rmax, tag_count = (1, 2, 3), 3, 4
        one = _rotation(np.pi / 4)
        product = np.kron(np.kron(one, one), one)

        def lastone(x, length):
            prefix = x >> (rmax - length)
            return 0 if prefix == 0 else length - ((prefix & -prefix).bit_length() - 1)

        def xor_word(state, active, length):
            result = np.empty_like(state)
            for x in range(8):
                for tag in range(tag_count):
                    result[x, tag ^ (lastone(x, length) if active else 0)] = state[x, tag]
            return result

        def word(state, active, length):
            state = xor_word(state, active, length)
            selected = np.empty_like(state)
            for tag in range(tag_count):
                unitary = np.ones((1, 1))
                for j in range(1, rmax + 1):
                    unitary = np.kron(unitary, one if active and tag <= j <= length
                                      else np.eye(2))
                selected[:, tag] = unitary @ state[:, tag]
            state = np.einsum("xy,ytj->xtj", product.T, selected)
            state = xor_word(state, active, length)
            return np.einsum("xy,ytj->xtj", product, state)

        for length in lengths:
            source = np.zeros((8, tag_count, 8), dtype=complex)
            source[:, 0] = np.eye(8)
            expected = np.zeros_like(source)
            balanced = real_frame_matrix(np.full((1 << length) - 1, np.pi / 4))
            expected[:, 0] = np.kron(balanced, np.eye(1 << (rmax - length)))
            np.testing.assert_allclose(word(source, True, length), expected, atol=ATOL)
            occupied = np.eye(8 * tag_count).reshape(8, tag_count, -1)
            np.testing.assert_allclose(word(occupied, False, length), occupied, atol=ATOL)


if __name__ == "__main__":
    unittest.main()
