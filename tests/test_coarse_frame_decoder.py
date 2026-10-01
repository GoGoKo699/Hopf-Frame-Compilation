"""Small semantic checks for the integer-histogram coarse-frame decoder."""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.coarse_frame_decoder import decode_coarse_frame_histograms
from compiler_robust_hopf.frames import real_tree_data
from tests.test_operator_source_compiler import H, X, _adjoint, _word_matrix


def _recursive_coarse(blocks, node=1):
    """An independent subtree composition for small dense references."""
    if node > len(blocks):
        return np.eye(1, dtype=complex)
    left = _recursive_coarse(blocks, 2 * node)
    right = _recursive_coarse(blocks, 2 * node + 1)
    half = len(left)
    children = np.zeros((2 * half, 2 * half), complex)
    children[:half, :half], children[half:, half:] = left, right
    root = np.eye(2 * half, dtype=complex)
    root[np.ix_([0, half], [0, half])] = blocks[node - 1]
    return children @ root


class CoarseFrameDecoderTests(unittest.TestCase):
    def test_complex_literal_blocks_match_empirical_quadrature_scores(self):
        word = [('T', 0), ('H', 0), ('TDG', 0), ('H', 0)]
        words = (word, _adjoint(word), [('S', 0)] + word + [('SDG', 0)])
        for height in (2, 3):
            size, scale = 1 << height, 2.5
            blocks = np.array([_word_matrix(1, words[node % 3]) for node in range(size - 1)])
            coarse = _recursive_coarse(blocks)
            walsh = np.eye(1)
            for _ in range(height):
                walsh = np.kron(walsh, H)
            records = [(setting, leaf, -1 if (leaf + setting) % 3 == 0 else 1)
                       for leaf in range(size) for setting in (0, 1)]
            records += [(0, 0, 1), (1, 1, -1), (1, size - 1, 1)]
            hx, hy = [0] * size, [0] * size
            for setting, leaf, sign in records:
                (hx if setting == 0 else hy)[leaf] += sign
            for angles in (np.linspace(-.43, .81, size - 1), np.zeros(size - 1)):
                derivatives = np.asarray(real_tree_data(angles).derivatives)
                transformed = (walsh @ coarse.conj().T @ derivatives.T).T
                direct = np.zeros(size - 1)
                for setting, leaf, sign in records:
                    column = transformed[:, leaf]
                    direct += 4 * scale * np.sqrt(size) * sign * (
                        column.real if setting == 0 else column.imag)
                direct /= len(records)
                snapshots = angles.copy(), blocks.copy(), hx.copy(), hy.copy()
                angles.setflags(write=False)
                blocks.setflags(write=False)
                actual = decode_coarse_frame_histograms(angles, blocks, hx, hy,
                                                       len(records), scale)
                np.testing.assert_allclose(actual, direct, atol=3e-12, rtol=0)
                for value, snapshot in zip((angles, blocks, hx, hy), snapshots, strict=True):
                    np.testing.assert_array_equal(value, snapshot)
                wrong_phase = decode_coarse_frame_histograms(angles, blocks.conj(), hx, hy,
                                                            len(records), scale)
                self.assertGreater(np.linalg.norm(actual - wrong_phase), .01)

    def test_large_integer_cancellation_survives_beyond_float_and_int64_precision(self):
        large = 1 << 80
        shots = 2 * large + 1
        # Walsh component zero is exactly one. Casting before the transform
        # would erase it; X directs just that component to the theta=0 tangent.
        actual = decode_coarse_frame_histograms([0], [X], [large + 1, -large],
                                               [0, 0], shots)
        self.assertEqual(actual[0], 4 / shots)
        self.assertGreater(actual[0], 0)

    def test_explicit_shot_count_retains_cancelled_records(self):
        identity = [np.eye(2)]
        zero = decode_coarse_frame_histograms([.2], identity, [0, 0], [0, 0], 12)
        np.testing.assert_array_equal(zero, [0])
        first = decode_coarse_frame_histograms([.2], identity, [1, 0], [0, 1], 2)
        second = decode_coarse_frame_histograms([.2], identity, [1, 0], [0, 1], 8)
        np.testing.assert_allclose(second, first / 4, atol=0, rtol=0)

    def test_invalid_geometry_histograms_and_normalizations_are_rejected(self):
        valid = dict(theta=[0], coarse_blocks=[np.eye(2)], hist_x=[0, 0],
                     hist_y=[0, 0], shots=2, coefficient_scale=1)
        changes = (
            dict(theta=[]), dict(theta=[0, 0]), dict(theta=[[0]]),
            dict(theta=[np.nan]), dict(theta=[np.inf]), dict(theta=[1j]),
            dict(theta=np.array([np.complex128(1 + 2j)], dtype=object)),
            dict(coarse_blocks=np.eye(2)), dict(coarse_blocks=[np.full((2, 2), np.nan)]),
            dict(hist_x=[0]), dict(hist_y=[[0, 0]]), dict(hist_x=[.5, 0]),
            dict(hist_y=[np.nan, 0]), dict(hist_x=['1', 0]), dict(hist_y=[1j, 0]),
            dict(hist_x=[2, 0], hist_y=[0, -1]),
            dict(shots=0), dict(shots=-1), dict(shots=2.0), dict(shots=True),
            dict(coefficient_scale=0), dict(coefficient_scale=-1),
            dict(coefficient_scale=np.inf), dict(coefficient_scale=np.nan),
            dict(coefficient_scale=np.complex128(1 + 2j)),
            dict(coefficient_scale=np.array(1 + 2j, dtype=object)),
        )
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                decode_coarse_frame_histograms(**(valid | change))


if __name__ == '__main__':
    unittest.main()
