"""Finite checks for product-coin dilation and actual-inverse flag recovery.

The source completions below are bounded dense fixtures, not native emitters.
Every rejection coordinate is retained when checking the recovery error.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import real_frame_matrix, real_tree_data


ATOL = 4e-11


def _fixture(q=2, tau=2, numerators=(1, 2, 3)):
    leaves, radix = 1 << q, 1 << tau
    atoms = radix**q
    coins = np.array([[z // radix**(q - j - 1) % radix for j in range(q)]
                      for z in range(atoms)])
    membership = np.zeros(atoms, dtype=int)
    prefixes = np.ones((atoms, q + 1), dtype=int)
    for z in range(atoms):
        node = 1
        for j in range(q):
            node = 2 * node + int(coins[z, j] >= numerators[node - 1])
            prefixes[z, j + 1] = node
        membership[z] = node - leaves
    embedding = np.zeros((atoms, leaves))
    for x in range(leaves):
        mask = membership == x
        embedding[mask, x] = 1 / np.sqrt(np.count_nonzero(mask))
    angles = np.arccos(np.sqrt(np.asarray(numerators) / radix))
    hadamard = np.ones((1, 1))
    for _ in range(q * tau):
        hadamard = np.kron(hadamard, np.array([[1, 1], [1, -1]]) / np.sqrt(2))
    last = np.max(np.where(coins != 0, np.arange(1, q + 1), 0), axis=1)
    depths = np.zeros(leaves, dtype=int)
    for node in range(1, leaves):
        depths[marker_label(node, q)] = node.bit_length()
    return {
        "q": q, "tau": tau, "radix": radix, "coins": coins,
        "membership": membership, "prefixes": prefixes, "embedding": embedding,
        "angles": angles, "frame": real_frame_matrix(angles),
        "hadamard": hadamard, "last": last, "depths": depths,
        "numerators": numerators,
    }


def _completion(target, phase):
    """Literal target column with a deliberately noncanonical completion."""
    target = np.asarray(target, dtype=complex)
    anchor = np.zeros(target.size, dtype=complex)
    anchor[0] = 1
    direction = anchor - target
    if np.linalg.norm(direction) < 1e-13:
        word = np.eye(target.size, dtype=complex)
    else:
        direction /= np.linalg.norm(direction)
        word = np.eye(target.size) - 2 * np.outer(direction, direction.conj())
    return word * np.exp(1j * phase * np.arange(target.size))[None, :]


def _apply_axes(state, word, axes):
    axes = tuple(axes)
    moved = np.moveaxis(state, axes, range(len(axes)))
    shape = moved.shape
    transformed = word @ moved.reshape(word.shape[1], -1)
    return np.moveaxis(transformed.reshape(shape), range(len(axes)), axes)


class _CoinWord:
    """Physical axes X, coins, source work, old tag, then references."""

    def __init__(self, *, phase=0.0, leakage=0.0):
        self.data = _fixture()
        self.sources = []
        self.delta = 0.0
        for x in range(4):
            target = np.zeros((16, 2))
            target[:, 0] = self.data["embedding"][:, x]
            word = _completion(target.ravel(), phase * (x + 1))
            # Mix a populated coin coordinate with a rejected source-work
            # coordinate at a different coin value. A product work rotation
            # could cancel in the decoder and would not test leakage return.
            populated = int(np.flatnonzero(target[:, 0])[0])
            pair = [2 * populated, 2 * ((populated + 9) % 16) + 1]
            rotation = np.eye(32)
            c, s = np.cos(leakage), np.sin(leakage)
            rotation[np.ix_(pair, pair)] = [[c, -s], [s, c]]
            actual = rotation @ word
            self.delta = max(self.delta, np.linalg.norm(actual[:, 0] - target.ravel()))
            self.sources.append(actual)

    def source(self, state, *, inverse=False):
        output = np.empty_like(state)
        for x, word in enumerate(self.sources):
            output[x] = _apply_axes(state[x], word.conj().T if inverse else word,
                                    (0, 1))
        return output

    def membership(self, state):
        output = np.empty_like(state)
        for x in range(4):
            for z, leaf in enumerate(self.data["membership"]):
                output[x ^ int(leaf), z] = state[x, z]
        return output

    def embedding(self, state, *, inverse=False):
        if inverse:
            return self.source(self.membership(state), inverse=True)
        return self.membership(self.source(state))

    def recover(self, state, *, active=True, offset=0):
        state = self.embedding(state)
        state = _apply_axes(state, self.data["hadamard"], (1,))
        erased = np.empty_like(state)
        for z, depth in enumerate(self.data["last"]):
            label = offset + int(depth) if active else 0
            for tag in range(state.shape[3]):
                erased[:, z, :, tag ^ label] = state[:, z, :, tag]
        state = _apply_axes(erased, self.data["hadamard"], (1,))
        return self.embedding(state, inverse=True)

    def tagged_columns(self, *, offset=0):
        tag_size = 1 << (offset + self.data["q"]).bit_length()
        state = np.zeros((4, 16, 2, tag_size, 4), dtype=complex)
        expected = np.zeros_like(state)
        for column, depth in enumerate(self.data["depths"]):
            state[:, 0, 0, offset + depth, column] = self.data["frame"][:, column]
            expected[:, 0, 0, 0, column] = self.data["frame"][:, column]
        return state, expected


class CoinDilationDepthTests(unittest.TestCase):
    def test_product_partition_isometry_and_every_scaling_flag(self):
        for data in (_fixture(), _fixture(3, 2, (1, 2, 3, 1, 3, 2, 1))):
            embedding, q = data["embedding"], data["q"]
            leaves = 1 << q
            np.testing.assert_allclose(embedding.T @ embedding, np.eye(leaves), atol=ATOL)
            tree = real_tree_data(data["angles"])
            for node in range(1, 2 * leaves):
                depth = node.bit_length() - 1
                mask = data["prefixes"][:, depth] == node
                expected = mask / np.sqrt(np.count_nonzero(mask))
                np.testing.assert_allclose(embedding @ tree.subtree[node], expected,
                                           atol=ATOL)
        self.assertEqual(np.bincount(_fixture()["membership"]).tolist(), [2, 2, 9, 3])

    def test_every_marker_has_literal_mean_zero_coin_and_uniform_suffix(self):
        for data in (_fixture(), _fixture(3, 2, (1, 2, 3, 1, 3, 2, 1))):
            q, radix = data["q"], data["radix"]
            for node in range(1, 1 << q):
                depth = node.bit_length() - 1
                numerator = data["numerators"][node - 1]
                probability = numerator / radix
                coin = np.r_[np.full(numerator, -np.sqrt((1 - probability) / numerator)),
                             np.full(radix - numerator,
                                     np.sqrt(probability / (radix - numerator)))]
                self.assertAlmostEqual(float(np.sum(coin)), 0.0, delta=ATOL)
                self.assertAlmostEqual(float(np.linalg.norm(coin)), 1.0, delta=ATOL)
                parent = data["prefixes"][:, depth] == node
                # There are K^(q-depth) suffix choices per allowed prefix.
                prefix_count = np.count_nonzero(parent) // radix**(q - depth)
                expected = parent * coin[data["coins"][:, depth]]
                expected /= np.sqrt(prefix_count * radix**(q - depth - 1))
                actual = data["embedding"] @ data["frame"][:, marker_label(node, q)]
                np.testing.assert_allclose(actual, expected, atol=ATOL)

    def test_hadamards_expose_exact_last_nonzero_coin_band_for_all_columns(self):
        for data in (_fixture(), _fixture(3, 2, (1, 2, 3, 1, 3, 2, 1))):
            transformed = data["hadamard"] @ data["embedding"] @ data["frame"]
            np.testing.assert_allclose(transformed[:, 0], np.eye(len(transformed))[:, 0],
                                       atol=ATOL)
            for column, depth in enumerate(data["depths"]):
                np.testing.assert_allclose(transformed[data["last"] != depth, column],
                                           0, atol=ATOL)
                self.assertAlmostEqual(float(np.linalg.norm(transformed[:, column])), 1.0,
                                       delta=ATOL)

    def test_physical_source_membership_and_actual_inverse(self):
        for phase in (0.0, 0.71):
            word = _CoinWord(phase=phase)
            state = np.zeros((4, 16, 2, 1, 4), dtype=complex)
            state[:, 0, 0, 0] = np.eye(4)
            expected = np.zeros_like(state)
            expected[0, :, 0, 0] = word.data["embedding"]
            output = word.embedding(state)
            np.testing.assert_allclose(output, expected, atol=ATOL)
            np.testing.assert_allclose(word.embedding(output, inverse=True), state, atol=ATOL)

    def test_full_fine_tag_erase_literal_phase_all_columns_and_offset(self):
        for phase in (0.0, 0.53):
            word = _CoinWord(phase=phase)
            for offset in (0, 3):
                state, expected = word.tagged_columns(offset=offset)
                output = word.recover(state, offset=offset)
                np.testing.assert_allclose(output, expected, atol=ATOL)
                # The trailing column axis can be an arbitrary entangled
                # external reference, rather than a classical test label.
                reference = np.array([[1, 2j], [-3j, 1], [2, -1j], [1j, 4]])
                reference /= np.linalg.norm(reference)
                np.testing.assert_allclose(output @ reference, expected @ reference,
                                           atol=ATOL)

    def test_actual_source_leakage_is_bounded_in_complete_output_norm(self):
        word = _CoinWord(phase=0.61, leakage=2e-3)
        state, expected = word.tagged_columns()
        output = word.recover(state)
        error = np.linalg.norm((output - expected).reshape(-1, 4), 2)
        self.assertGreater(error, word.delta / 10)
        self.assertLessEqual(error, 2 * word.delta + ATOL)
        rejected = output.copy()
        rejected[:, 0, 0, 0] = 0
        self.assertGreater(np.linalg.norm(rejected), word.delta / 10)
        reference = np.array([[1, 1j], [2j, -1], [3, 2], [-1j, 1]])
        reference /= np.linalg.norm(reference)
        self.assertLessEqual(np.linalg.norm((output - expected) @ reference),
                             2 * word.delta + ATOL)

    def test_inactive_word_cancels_on_arbitrary_occupied_work_and_reference(self):
        word = _CoinWord(phase=0.43, leakage=0.017)
        values = np.arange(4 * 16 * 2 * 4 * 2).reshape(4, 16, 2, 4, 2)
        state = np.exp(0.137j * values) * (1 + (values % 11) / 13)
        state /= np.linalg.norm(state)
        np.testing.assert_allclose(word.recover(state, active=False), state, atol=ATOL)


if __name__ == "__main__":
    unittest.main()
