"""Bounded algebra audits for fast capped geometric reference readout.

The matrices below are dense verification completions, not native source
oracles.  They audit the truncation error, its full-space completion, and
the original reference's complete two-block flag with occupied private
work and actual inverses.  They establish no resource bound or emitter.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 8e-12


def _capped(size, dimension=None):
    """Literal positive geometric column, with the last probability capped."""
    column = np.zeros(size if dimension is None else dimension)
    column[:size - 1] = np.exp2(-np.arange(1, size) / 2)
    column[size - 1] = np.exp2(-(size - 1) / 2)
    return column


def _completion(column, seed):
    """Fix the source column and independently mix all unspecified columns."""
    difference = column.astype(complex).copy()
    difference[0] -= 1
    word = np.eye(len(column), dtype=complex)
    if np.linalg.norm(difference):
        word -= 2 * np.outer(difference, difference.conj()) / np.vdot(
            difference, difference)
    rng = np.random.default_rng(seed)
    random = rng.normal(size=(len(column) - 1, len(column) - 1)) + 1j * rng.normal(
        size=(len(column) - 1, len(column) - 1))
    complement, _ = np.linalg.qr(random)
    word[:, 1:] = word[:, 1:] @ complement
    return word


def _plane_rotation(source, target):
    """Analysis-only nearest plane rotation; no gate for it is supplied."""
    identity = np.eye(len(source), dtype=complex)
    if np.array_equal(source, target):
        return identity
    cosine = float(np.vdot(source, target).real)
    orthogonal = target - cosine * source
    sine = np.linalg.norm(orthogonal)
    orthogonal /= sine
    return identity + (cosine - 1) * (
        np.outer(source, source.conj()) + np.outer(orthogonal, orthogonal.conj())
    ) + sine * (
        np.outer(orthogonal, source.conj()) - np.outer(source, orthogonal.conj()))


def _reference(size=8, retained=4, uniform_size=1, seed=291):
    # Order each augmented block as uniform, original coin, private bit.
    uniform = np.ones(uniform_size) / np.sqrt(uniform_size)
    source = np.kron(uniform, np.kron(_capped(retained, size), [1., 0.]))
    target = np.kron(uniform, np.kron(_capped(size), [1., 0.]))
    actual = _completion(source, seed)
    rotation = _plane_rotation(source, target)
    return actual, rotation @ actual, source, target, rotation


def _original_bands(size=8):
    """Independent original-coin basis, embedded with zero private bits."""
    original = _capped(size)
    basis = _completion(original, 715)
    basis = np.kron(basis, np.array([[1.], [0.]]))
    root, markers = basis[:, 0], basis[:, 1:]
    columns = [np.kron(root, root)]
    bands = [0]
    # Last nonzero augmented block is 1, then 2.  These size-squared columns
    # span the entire pair of original coins with initialized private work.
    columns.extend(np.kron(markers[:, j], root) for j in range(size - 1))
    bands.extend([1] * (size - 1))
    columns.extend(np.kron(basis[:, i], markers[:, j])
                   for i in range(size) for j in range(size - 1))
    bands.extend([2] * (size * (size - 1)))
    return np.column_stack(columns), bands


def _readout(word, state, *, inverse=None, enabled=True,
             omit_private=False, omit_high_coin=False):
    """Reference inverse, label XOR, reference forward, keeping all work."""
    inverse = word.conj().T if inverse is None else inverse
    transformed = np.tensordot(inverse, state, axes=(1, 0))
    if enabled:
        labels = np.arange(4)
        block_dimension = round(np.sqrt(len(word)))
        for location in range(len(word)):
            blocks = [location // block_dimension, location % block_dimension]
            if omit_private:
                blocks = [block // 2 for block in blocks]
            if omit_high_coin:
                blocks = [block % (block_dimension // 2) for block in blocks]
            band = 2 if blocks[1] else (1 if blocks[0] else 0)
            transformed[location] = transformed[location, labels ^ band].copy()
    return np.tensordot(word, transformed, axes=(1, 0))


class FastCappedGeometricReferenceTests(unittest.TestCase):
    def test_literal_cap_overlap_distance_and_precision_cutoff(self):
        for size in (2, 4, 8, 16, 32, 64, 128):
            original = _capped(size)
            for retained in (1 << j for j in range(1, size.bit_length())):
                with self.subTest(size=size, retained=retained):
                    shortened = _capped(retained, size)
                    np.testing.assert_allclose(np.linalg.norm(shortened), 1,
                                               atol=ATOL, rtol=0)
                    distance_squared = np.linalg.norm(original - shortened)**2
                    if retained == size:
                        np.testing.assert_array_equal(shortened, original)
                    else:
                        expected = (4 - 2 * np.sqrt(2)) * 2.**(-retained)
                        self.assertAlmostEqual(distance_squared, expected, places=14)
                        self.assertAlmostEqual(np.vdot(original, shortened),
                                               1 - (2 - np.sqrt(2)) * 2.**(-retained),
                                               places=14)
                        self.assertLess(distance_squared, 2.**(1 - retained))
            for precision in range(1, 17):
                cutoff = 1 << (2 * precision + 1).bit_length()
                retained = min(size, cutoff)
                self.assertLessEqual(np.linalg.norm(original - _capped(retained, size)),
                                     2.**(-precision))
        # Cutting and renormalizing the original prefix gives a different
        # column; the exact capped source has a distinct final amplitude.
        prefix = _capped(8)[:4]
        prefix /= np.linalg.norm(prefix)
        self.assertGreater(np.linalg.norm(prefix - _capped(4)), .05)

    def test_full_space_plane_completion_and_augmented_anchor(self):
        for size, retained in ((4, 2), (8, 4), (8, 8)):
            actual, ideal, source, target, rotation = _reference(
                size, retained, uniform_size=2)
            identity = np.eye(len(actual))
            with self.subTest(size=size, retained=retained):
                np.testing.assert_allclose(actual[:, 0], source, atol=ATOL, rtol=0)
                np.testing.assert_allclose(ideal[:, 0], target, atol=ATOL, rtol=0)
                np.testing.assert_allclose(ideal.conj().T @ ideal, identity,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(rotation.conj().T @ rotation, identity,
                                           atol=ATOL, rtol=0)
                error = np.linalg.norm(source - target)
                self.assertAlmostEqual(np.linalg.norm(rotation - identity, 2),
                                       error, places=12)
                self.assertAlmostEqual(np.linalg.norm(actual - ideal, 2),
                                       error, places=12)
                anchor = np.outer(ideal[:, 0], ideal[:, 0].conj())
                np.testing.assert_allclose(anchor, np.outer(target, target.conj()),
                                           atol=ATOL, rtol=0)
                # The comparison acts trivially on occupied private scratch;
                # its rotated plane lies entirely inside the zero-scratch sector.
                np.testing.assert_allclose(rotation[:, 1::2], identity[:, 1::2],
                                           atol=ATOL, rtol=0)

    def test_every_original_band_uses_all_augmented_coordinates(self):
        first, first_ideal, u, v, _ = _reference(seed=291)
        second, second_ideal, _, _, _ = _reference(seed=941)
        actual = np.kron(first, second)
        ideal = np.kron(first_ideal, second_ideal)
        embedding, bands = _original_bands()
        count, dimension = len(bands), len(actual)
        np.testing.assert_allclose(embedding.conj().T @ embedding, np.eye(count),
                                   atol=ATOL, rtol=0)
        # Test every input label as well as every root/marker column.
        initial = np.zeros((dimension, 4, 4 * count), dtype=complex)
        wanted = np.zeros_like(initial)
        for label in range(4):
            for column, band in enumerate(bands):
                initial[:, label, count * label + column] = embedding[:, column]
                wanted[:, label ^ band, count * label + column] = embedding[:, column]
        np.testing.assert_allclose(_readout(ideal, initial), wanted, atol=ATOL, rtol=0)
        result = _readout(actual, initial)
        bound = 4 * np.linalg.norm(u - v)
        self.assertLess(bound, 2)  # The sandwich estimate must be nonvacuous.
        self.assertLessEqual(np.linalg.norm((result - wanted).reshape(
            4 * dimension, 4 * count), 2), bound + ATOL)
        # Arbitrary permitted completion columns can put a marker in private
        # scratch.  Neither private scratch nor the original high coin bit
        # may be silently dropped from the transported zero test.
        for omission in ({'omit_private': True}, {'omit_high_coin': True}):
            wrong = _readout(ideal, initial, **omission)
            self.assertGreater(np.linalg.norm(wrong - wanted), .1)
        # Original tail weight classes remain positive.  The physical short
        # source deletes those amplitudes, so it cannot be substituted into
        # the mathematical matching while claiming the same exact target.
        self.assertGreater(np.linalg.norm(v[8::2]), 0)
        np.testing.assert_array_equal(u[8::2], 0)
        self.assertGreater(np.linalg.norm(_readout(actual, initial[:, :, :1])
                                         - wanted[:, :, :1]), .01)

    def test_actual_inverse_inactive_identity_and_entangled_occupied_work(self):
        first, first_ideal, u, v, _ = _reference(seed=291)
        second, second_ideal, _, _, _ = _reference(seed=941)
        actual = np.kron(first, second)
        ideal = np.kron(first_ideal, second_ideal)
        rng = np.random.default_rng(20261016)
        # H, both augmented blocks, label, dirty bit, external reference bit.
        state = rng.normal(size=(2, len(actual), 4, 2, 2)) + 1j * rng.normal(
            size=(2, len(actual), 4, 2, 2))
        state /= np.linalg.norm(state)
        result = np.stack((_readout(actual, state[0], enabled=False),
                           _readout(actual, state[1])))
        wanted = np.stack((state[0], _readout(ideal, state[1])))
        np.testing.assert_allclose(result[0], state[0], atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.norm(result - wanted),
                             4 * np.linalg.norm(u - v) * np.linalg.norm(state[1])
                             + ATOL)
        self.assertAlmostEqual(np.linalg.norm(result), 1, places=12)
        # Replacing the actual inverse by the comparison inverse destroys
        # inactive cancellation even though both have nearby source columns.
        stale = _readout(actual, state[0], enabled=False, inverse=ideal.conj().T)
        self.assertGreater(np.linalg.norm(stale - state[0]), .01)


if __name__ == '__main__':
    unittest.main()
